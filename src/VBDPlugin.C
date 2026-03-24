


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

static PRM_Name	angleName("angle", "Angle"); // (internal, label)
static PRM_Name	stepSizeName("stepSize", "StepSize");
static PRM_Name	iterationsName("iterations", "Iterations");
static PRM_Name	fileName("grammarFilePath", "GrammarFilePath");

static PRM_Default angleDefault(30.0);
static PRM_Default stepSizeDefault(1.0);
static PRM_Default iterationsDefault(2);
static PRM_Default fileDefault(0, "UNDEFINED");

static PRM_Range angleRange(PRM_RANGE_RESTRICTED, 0, PRM_RANGE_UI, 360);
static PRM_Range stepSizeRange(PRM_RANGE_RESTRICTED, 0, PRM_RANGE_UI, 5);
static PRM_Range iterationsRange(PRM_RANGE_RESTRICTED, 0, PRM_RANGE_UI, 8);

PRM_Template
SOP_VBD::myTemplateList[] = {
   PRM_Template(PRM_FLT, PRM_Template::PRM_EXPORT_MIN, 1, &angleName, &angleDefault, 0, &angleRange),
   PRM_Template(PRM_FLT, PRM_Template::PRM_EXPORT_MIN, 1, &stepSizeName, &stepSizeDefault, 0, &stepSizeRange),
   PRM_Template(PRM_INT, PRM_Template::PRM_EXPORT_MIN, 1, &iterationsName, &iterationsDefault, 0, &iterationsRange),
   PRM_Template(PRM_FILE, 1, &fileName, &fileDefault),
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
	fpreal currTime = context.getTime();

	float angle; angle = ANGLE(currTime);
    float stepSize; stepSize = STEP_SIZE(currTime);
    uint iterationCount; iterationCount = static_cast<uint>(ITERATIONS(currTime));
    UT_String hFilePath;
    FILE_PATH(currTime, hFilePath);
    string filePath = hFilePath.toStdString();

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

    if (gdp->getP()->getDataId() != inputGeoDataID) {
        inputGeoDataID = gdp->getP()->getDataId();
        std::cout << "New Input Mesh of Topo ID " << inputGeoDataID << ", Resetting Sim!" << std::endl;
        uPtr<HalfEdgeMesh> inputMesh = mkU<HalfEdgeMesh>();
        inputMesh->CreateFromGUDetail(gdp);
        vbdSolver.ResetSimulation(std::move(inputMesh));
    }

	// Start the interrupt server
	if (boss->opStart("Building Sim Frame"))
	{
        std::cerr << "Frame " << context.getFrame() << std::endl;
        vbdSolver.SimulateUpToFrame(context.getFrame());
        vbdSolver.lastSimulatedMesh->LoadIntoExistingTopologicallySameHoudiniMesh(gdp);

     //   GA_Offset ptoff;
     //   UT_Vector3 tpos;

     //   vec3 sta = vec3(0);// b.first;
     //   vec3 end = vec3(1, 0, 0);// b.second;
     //   vec3 mid = 0.5f * (sta + end);

     //   UT_Vector3 hmid = UT_Vector3(mid[0], mid[1], mid[2]);

     //   vec3 baseDir = vec3(0.0f, 1.0f, 0.0f);
     //   vec3 dir = end - sta;
     //   float length = glm::length(dir);
     //   dir = normalize(dir);
     //   vec3 axis = normalize(glm::cross(baseDir, dir));
     //   float angle = std::acos(glm::dot(dir, baseDir));

     //   UT_Matrix4 transform;
     //   transform.identity();
     //   transform.scale(1.0f, length, 1.0f, 1.0f);
     //   transform.rotate(UT_Vector3(axis[0], axis[1], axis[2]), angle);
     //   transform.translate(hmid);

     //   float radius = 0.5f;
     //   int cols = divisions;// 4;

     //   // set pos of vertices
     //   for (int i = 0; i < cols; i++)
     //   {
     //       float angle = (float)i / cols * M_PI * 2.0f;
     //       tpos.x() = radius * cos(angle);
     //       tpos.y() = 0.5f;
     //       tpos.z() = radius * sin(angle);
     //       tpos = tpos * transform;

     //       ptoff = gdp->appendPointOffset();
     //       // std::cerr << ptoff << std::endl;
     //       gdp->setPos3(ptoff, tpos);
     //   }

     //   for (int i = 0; i < cols; i++)
     //   {
     //       float angle = (float)i / cols * M_PI * 2.0f;
     //       tpos.x() = radius * cos(angle);
     //       tpos.y() = -0.5f;
     //       tpos.z() = radius * sin(angle);
     //       tpos = tpos * transform;

     //       ptoff = gdp->appendPointOffset();
     //       // std::cerr << ptoff << std::endl;
     //       gdp->setPos3(ptoff, tpos);
     //   }

     //   for (int i = 0; i < cols; i++)
     //   {
     //       int next = (i + 1) % cols;

     //       int lsm[4] = { i, next,cols + next,cols + i };
     //       std::cerr << lsm[0] << '\t' << lsm[1] << '\t' << lsm[2] << '\t' << lsm[3] << std::endl;

     //       std::cerr << "Point Count: " << gdp->getNumPoints() << std::endl;
     //       GU_PrimPoly* poly = GU_PrimPoly::build(gdp, 4, GU_POLY_CLOSED);
     //       std::cerr << "Point Count: " << gdp->getNumPoints() << std::endl;

     //       poly->setPointOffset(0, lsm[0]);// off + i);
     //       poly->setPointOffset(1, lsm[1]);// off + next);
     //       poly->setPointOffset(2, lsm[2]);// off + cols + next);
     //       poly->setPointOffset(3, lsm[3]);// off + cols + i);

     //   }
	    //select(GU_SPrimitive);
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

