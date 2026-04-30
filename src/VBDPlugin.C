


#include <UT/UT_DSOVersion.h>
//#include <RE/RE_EGLServer.h>


#include <UT/UT_Math.h>
#include <UT/UT_Interrupt.h>
#include <GU/GU_Detail.h>
#include <GU/GU_PrimPoly.h>
#include <CH/CH_LocalVariable.h>
#include <PRM/PRM_Include.h>
#include <PRM/PRM_SpareData.h>
#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>

#include "half-edge-mesh.h"

#include <limits.h>
#include "VBDPlugin.h"

#include <cmath>
#include <string>
#include <glm.hpp>
#include <iostream>

#include <vector>
#include <unordered_set>
#include "defines.h"

//#include <CH/CH_Manager.h>
//#include <OP/OP_Director.h>
using namespace glm;
using namespace std;
using namespace HDK_Sample;

/// newSopOperator is the hook that Houdini grabs from this dll
/// and invokes to register the SOP.  In this case we add ourselves
/// to the specified operator table.
void
newSopOperator(OP_OperatorTable *table)
{
    table->addOperator(
	    new OP_Operator("VBDSolver",			// Internal name
			    "VBDSolver",			// UI name
			     SOP_VBD::myConstructor,	// How to build the SOP
			     SOP_VBD::myTemplateList,	// My parameters
			     1,				// Min # of sources
			     1,				// Max # of sources
			     SOP_VBD::myVariables,	// Local variables
			     OP_FLAG_GENERATOR)		// Flag it as generator
	    );
}

//
static PRM_Name timeScaleName("timeScale", "TimeScale");
static PRM_Name subStepsName("subSteps", "Substeps");
static PRM_Name iterationCountName("gaussSeidelIterations", "GaussSeidelIterations");

static PRM_Name physicsMaterialName("physicsMaterial", "PhysicsMaterial");

static PRM_Name springConstantName("springConstant", "SpringConstant");
static PRM_Name restLengthName("restLength", "RestLength");

static PRM_Name areaChangeResistanceName("areaChangeResistance", "AreaChangeResistance");
static PRM_Name shearResistanceName("shearResistanceName", "ShearResistanceName");

static PRM_Name constraintGroupNameName("constraintGroupName", "ConstraintGroupName");

//
static PRM_Default timeScaleDefault(1.0);
static PRM_Default subStepsDefault(1);
static PRM_Default iterationCountDefault(5);

static PRM_Default physicsMaterialDefault(0);

static PRM_Default springConstantDefault(150.0);
static PRM_Default restLengthDefault(0.3);

static PRM_Default areaChangeResistanceDefault(1.0);
static PRM_Default shearResistanceDefault(1.0);

static PRM_Default constraintGroupNameDefault(0, "ConstraintGroup");

PRM_Template
SOP_VBD::myTemplateList[] = {
   PRM_Template(PRM_FLT, PRM_Template::PRM_EXPORT_MIN, 1, &timeScaleName, &timeScaleDefault, 0),
   PRM_Template(PRM_INT, PRM_Template::PRM_EXPORT_MIN, 1, &subStepsName, &subStepsDefault, 0),
   PRM_Template(PRM_INT, PRM_Template::PRM_EXPORT_MIN, 1, &iterationCountName, &iterationCountDefault, 0),
   PRM_Template(PRM_INT, PRM_Template::PRM_EXPORT_MIN, 1, &springConstantName, &springConstantDefault, 0),
   PRM_Template(PRM_FLT, PRM_Template::PRM_EXPORT_MIN, 1, &physicsMaterialName, &physicsMaterialDefault, 0),
   PRM_Template(PRM_FLT, PRM_Template::PRM_EXPORT_MIN, 1, &restLengthName, &restLengthDefault, 0),
   PRM_Template(PRM_FLT, PRM_Template::PRM_EXPORT_MIN, 1, &areaChangeResistanceName, &areaChangeResistanceDefault, 0),
   PRM_Template(PRM_FLT, PRM_Template::PRM_EXPORT_MIN, 1, &shearResistanceName, &shearResistanceDefault, 0),
   PRM_Template(PRM_STRING, PRM_Template::PRM_EXPORT_MIN, 1, &constraintGroupNameName, &constraintGroupNameDefault, 0),
   PRM_Template()
};

enum {
	VAR_PT,		// ptnum
	VAR_NPT		// npts
};

// Local variable info, updated if looping through some pieces
CH_LocalVariable
SOP_VBD::myVariables[] = {
    { "PT",	VAR_PT, 0 },
    { "NPT",	VAR_NPT, 0 },
    { 0, 0, 0 },
};

bool
SOP_VBD::evalVariableValue(fpreal &val, int index, int thread)
{
    // myCurrPoint will be negative when we're not cooking so only try to
    // handle the local variables when we have a valid myCurrPoint index.
    if (myCurrPoint >= 0)
    {
	// Note that "gdp" may be null here, so we do the safe thing
	// and cache values we are interested in.
	switch (index)
	{
	    case VAR_PT:
		val = (fpreal) myCurrPoint;
		return true;
	    case VAR_NPT:
		val = (fpreal) myTotalPoints;
		return true;
	    default:
		/* do nothing */;
	}
    }
    // Not one of our variables, must delegate to the base class.
    return SOP_Node::evalVariableValue(val, index, thread);
}

OP_Node *
SOP_VBD::myConstructor(OP_Network *net, const char *name, OP_Operator *op)
{
    return new SOP_VBD(net, name, op);
}

SOP_VBD::SOP_VBD(OP_Network *net, const char *name, OP_Operator *op)
	: SOP_Node(net, name, op)
{
    myCurrPoint = -1;	// To prevent garbage values from being returned
}

SOP_VBD::~SOP_VBD() {}

unsigned
SOP_VBD::disableParms()
{
    return 0;
}


//// TODO: This could be put into a 'glue' class that acts as info relay between houdini sop and our sim
//// Converts input mesh into an adjacency list format
//int SOP_VBD::convertMeshToAdjacency(OP_Context &context, int inputIndex) {
//    if (lockInputs(context) >= UT_ERROR_ABORT) {
//        // Inputs won't be changed while we're looked; ensure input data doesn't change during cook
//        // Will auto unlock when context goes out of scope
//        std::cerr << "Lock Failed" << std::endl;
//        return -1;
//    }
//
//    const GU_Detail* geo = inputGeo(inputIndex, context);
//
//    uPtr<HalfEdgeMesh> mesh = mkU<HalfEdgeMesh>();
//    mesh->CreateFromGUDetail(geo);
//     
//    return 0;
//}

OP_ERROR
SOP_VBD::cookMySop(OP_Context &context)
{
    /*CH_Manager* channelManager = OPgetDirector()->getChannelManager();*/
    float fps = 24.0f;// channelManager->getSamplesPerSec();

	fpreal currTime = context.getTime();

    std::cerr << "Current Time: " << currTime << std::endl;

    // convertMeshToAdjacency(context, 0);

    int			 divisions;
    int			 xcoord =0, ycoord = 1, zcoord =2;
    UT_Interrupt	*boss;

    // Since we don't have inputs, we don't need to lock them.

    divisions  = glm::ceil(4+currTime);
    myTotalPoints = divisions;		// Set the NPT local variable value, TODO: NOT ACCURATE RN FOR THE USER!
    myCurrPoint   = 0;			// Initialize the PT local variable

    // Check to see that there hasn't been a critical error in cooking the SOP.
    if (error() < UT_ERROR_ABORT)
    {
	boss = UTgetInterrupt();
	if (divisions < 4)
	{
	    addWarning(SOP_MESSAGE, "Invalid divisions (just a test warning)");
	    divisions = 4;
	}
    if (lockInputs(context) >= UT_ERROR_ABORT) {
        // Inputs won't be changed while we're looked; ensure input data doesn't change during cook
        // Will auto unlock when context goes out of scope
        std::cerr << "Lock Failed" << std::endl;
    }
    duplicateSource(0, context);//gdp->clearAndDestroy();  // Clear all geo of this node

    if (gdp->getP()->getDataId() != inputGeoDataID) { // TODO: check for change in group and stuff
        inputGeoDataID = gdp->getP()->getDataId();

        std::cout << "New Input Mesh of Topo ID " << inputGeoDataID << ", Resetting Sim!" << std::endl;
        uPtr<HalfEdgeMesh> inputMesh = mkU<HalfEdgeMesh>();
        UT_String groupNameStr; evalString(groupNameStr, constraintGroupNameName.getToken(), 0, currTime);
        inputMesh->CreateFromGUDetail(gdp, groupNameStr);
        vbdSolver.ResetSimulation(std::move(inputMesh));
    }

	// Start the interrupt server
	if (boss->opStart("Building Sim Frame"))
	{
        vbdSolver.dt = evalFloat(timeScaleName.getToken(), 0, currTime) / fps;
        
        vbdSolver.g = vec3(0.0f, -0.98f, 0.0f); // TODO: parametrize
        vbdSolver.iterCount = evalInt(iterationCountName.getToken(), 0, currTime);
        vbdSolver.currMaterial = evalInt(physicsMaterialName.getToken(), 0, currTime);

        vbdSolver.k = evalFloat(springConstantName.getToken(), 0, currTime);
        vbdSolver.restLen = evalFloat(restLengthName.getToken(), 0, currTime);

        vbdSolver.u = evalFloat(shearResistanceName.getToken(), 0, currTime);
        vbdSolver.lambda = evalFloat(areaChangeResistanceName.getToken(), 0, currTime);

        std::cerr << "Frame " << context.getFrame() << std::endl;
        vbdSolver.SimulateUpToFrame(context.getFrame());
        vbdSolver.lastSimulatedMesh->LoadIntoExistingTopologicallySameHoudiniMesh(gdp);
	}

	// Tell the interrupt server that we've completed. Must do this
	// regardless of what opStart() returns.
	boss->opEnd();
    }

    myCurrPoint = -1;

    // This flag gets reset to false every beginning cook, so we have to set back to true
    flags().setTimeDep(true);

    return error();
}

