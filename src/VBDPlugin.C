


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
#include <OP/OP_Director.h>

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
			     2,				// Max # of sources
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

static PRM_Name collisionThresholdName("collisionThreshold", "CollisionThreshold");
static PRM_Name collisionCoefficientName("collisionCoefficient", "CollisionCoefficient");

static PRM_Name constraintGroupNameName("constraintGroupName", "ConstraintGroupName");

static PRM_Name gravityName("gravity", "Gravity");

//
static PRM_Default timeScaleDefault(1.0);
static PRM_Default subStepsDefault(1);
static PRM_Default iterationCountDefault(5);

static PRM_Default physicsMaterialDefault(0);

static PRM_Default springConstantDefault(150.0);
static PRM_Default restLengthDefault(1.0f);

static PRM_Default areaChangeResistanceDefault(1.0);
static PRM_Default shearResistanceDefault(1.0);

static PRM_Default collisionThresholdDefault(0.1f);
static PRM_Default collisionCoefficientDefault(1e6);

static PRM_Default constraintGroupNameDefault(0, "ConstraintGroup");

static PRM_Default gravityDefault[] = {
    PRM_Default(0.0f),
    PRM_Default(-0.98f),
    PRM_Default(0.0f)
};

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
   PRM_Template(PRM_FLT, PRM_Template::PRM_EXPORT_MIN, 1, &collisionThresholdName, &collisionThresholdDefault, 0),
   PRM_Template(PRM_FLT, PRM_Template::PRM_EXPORT_MIN, 1, &collisionCoefficientName, &collisionCoefficientDefault, 0),
   PRM_Template(PRM_STRING, PRM_Template::PRM_EXPORT_MIN, 1, &constraintGroupNameName, &constraintGroupNameDefault, 0),
   PRM_Template(PRM_XYZ, PRM_Template::PRM_EXPORT_MIN, 3, &gravityName, gravityDefault, 0),
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

/*static PRM_Name timeScaleName("timeScale", "TimeScale");
static PRM_Name subStepsName("subSteps", "Substeps");
static PRM_Name iterationCountName("gaussSeidelIterations", "GaussSeidelIterations");

static PRM_Name physicsMaterialName("physicsMaterial", "PhysicsMaterial");

static PRM_Name springConstantName("springConstant", "SpringConstant");
static PRM_Name restLengthName("restLength", "RestLength");

static PRM_Name areaChangeResistanceName("areaChangeResistance", "AreaChangeResistance");
static PRM_Name shearResistanceName("shearResistanceName", "ShearResistanceName");

static PRM_Name constraintGroupNameName("constraintGroupName", "ConstraintGroupName");
*/

// TODO: put in solver
SolverParams SOP_VBD::GetParams(fpreal time) {
    SolverParams params;

    CH_Manager* channelManager = OPgetDirector()->getChannelManager();
    float fps = channelManager->getSamplesPerSec();

    params.frameDt = evalFloat(timeScaleName.getToken(), 0, time) / fps;
    params.subSteps = evalInt(subStepsName.getToken(), 0, time);

    params.g = vec3(
        evalFloat(gravityName.getToken(), 0, time),
        evalFloat(gravityName.getToken(), 1, time),
        evalFloat(gravityName.getToken(), 2, time)
    );
    params.iterCount = evalInt(iterationCountName.getToken(), 0, time);
    params.currMaterial = evalInt(physicsMaterialName.getToken(), 0, time);

    params.k = evalFloat(springConstantName.getToken(), 0, time);
    params.restLen = evalFloat(restLengthName.getToken(), 0, time);

    params.u = evalFloat(shearResistanceName.getToken(), 0, time);
    params.lambda = evalFloat(areaChangeResistanceName.getToken(), 0, time);

    params.m = 1.0f; // TODO: parametrize

    params.collisionThreshold = evalFloat(collisionThresholdName.getToken(), 0, time);
    params.kc = evalFloat(collisionCoefficientName.getToken(), 0, time);

    return params;
}

OP_ERROR
SOP_VBD::cookMySop(OP_Context &context)
{
    

	fpreal currTime = context.getTime();

    UT_Interrupt* boss;
    myCurrPoint   = 0;			// Initialize the PT local variable


    // Check to see that there hasn't been a critical error in cooking the SOP.
    if (error() < UT_ERROR_ABORT)
    {
	boss = UTgetInterrupt();
    if (lockInputs(context) >= UT_ERROR_ABORT) {
        // Inputs won't be changed while we're looked; ensure input data doesn't change during cook
        // Will auto unlock when context goes out of scope
        std::cerr << "Lock Failed" << std::endl;
    }
    duplicateSource(0, context);//gdp->clearAndDestroy();  // Clear all geo of this node

    // Collision
    const GU_Detail* collisionMesh = inputGeo(1, context);
    int prevCollisionGeoDataID = collisionGeoDataID;
    collisionGeoDataID = collisionMesh ? collisionMesh->getP()->getDataId() : -1;
    // std::cerr << "Collision data: " << collisionGeoDataID << std::endl;

    // Constraint Group
    int prevConstraintGroupID = constraintGroupID;
    UT_String groupNameStr; evalString(groupNameStr, constraintGroupNameName.getToken(), 0, currTime);
    const GA_PointGroup* group = gdp->findPointGroup(groupNameStr);
    if (group == nullptr) {
        constraintGroupID = -1;
        // std::cerr << "Constraint Group \"" << groupNameStr << "\" doesn't exist!" << std::endl;
    } else {
        constraintGroupID = group->getDataId();
        // std::cerr << "Constraint Group ID: " << constraintGroupID << std::endl;
    }

    // Check for param changes
    bool paramsChanged = false;
    auto channelManager = OPgetDirector()->getChannelManager();
    int frameCount =
        channelManager->getSample(channelManager->getGlobalEnd()) -
        channelManager->getSample(channelManager->getGlobalStart()) + 1;
    if (cachedParams.size() != frameCount) {
        std::cerr << "Was size " << cachedParams.size() << " but frame count is " << frameCount << std::endl;
        paramsChanged = true;
        cachedParams.resize(frameCount);
    }
    for (int i = 0; i < frameCount; i++) {
        auto newParam = GetParams(
            channelManager->getTime(
                channelManager->getSample(channelManager->getGlobalStart()) + i
            )
        );
        if (!paramsChanged && !(cachedParams[i] == newParam)) {
            std::cerr << "New param found at " << i << std::endl;
            paramsChanged = true;
        }
        cachedParams[i] = newParam;
        // std::cerr << "putting param at " << i << " with mat " << cachedParams[i].currMaterial << std::endl;
    }

    // Invalidate Solver Cache
    if (paramsChanged) {
        // CACHE INVALIDATION
        // std::cerr << "Params changed at " << context.getFrame() << std::endl;
        vbdSolver.TruncateSimulation(context.getFrame());
        vbdSolver.SetDirty();
    }

    // Check should reset sim?
    if (
        prevConstraintGroupID != constraintGroupID || // group changed
        gdp->getP()->getDataId() != inputGeoDataID || // input geo changed
        collisionGeoDataID != prevCollisionGeoDataID // input collision changed
        ) {
        inputGeoDataID = gdp->getP()->getDataId();

        std::cerr << "Convert input mesh" << std::endl;
        uPtr<HalfEdgeMesh> inputMesh = mkU<HalfEdgeMesh>();        
        inputMesh->CreateFromGUDetail(gdp, group);

        uPtr<HalfEdgeMesh> collisionGeo = nullptr;
        if (collisionGeoDataID != -1) {
            collisionGeo = mkU<HalfEdgeMesh>();
            collisionGeo->CreateFromGUDetail(collisionMesh, nullptr);
        }

        std::cerr << "Time to reset sim" << std::endl;
        vbdSolver.ResetSimulation(std::move(inputMesh), collisionGeo ? std::move(collisionGeo) : nullptr);
    }

	// Start the interrupt server
	if (boss->opStart("Building Sim Frame"))
	{
        // std::cerr << "Frame " << context.getFrame() << std::endl;
        vbdSolver.SimulateUpToFrame(context.getFrame());
        vbdSolver.GetMesh(context.getFrame())->LoadIntoExistingTopologicallySameHoudiniMesh(gdp);
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

