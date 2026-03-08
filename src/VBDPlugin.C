


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


#include <limits.h>
#include "VBDPlugin.h"

#include <cmath>
#include <string>
#include <glm.hpp>
#include <iostream>
using namespace glm;
using namespace std;
using namespace HDK_Sample;

//
// Help is stored in a "wiki" style text file. 
//
// See the sample_install.sh file for an example.
//
// NOTE : Follow this tutorial if you have any problems setting up your visual studio 2008 for Houdini 
//  http://www.apileofgrains.nl/setting-up-the-hdk-for-houdini-12-with-visual-studio-2008/


///
/// newSopOperator is the hook that Houdini grabs from this dll
/// and invokes to register the SOP.  In this case we add ourselves
/// to the specified operator table.
///
void
newSopOperator(OP_OperatorTable *table)
{
    table->addOperator(
	    new OP_Operator("VBDSolver",			// Internal name
			    "VBDSolver",			// UI name
			     SOP_VBD::myConstructor,	// How to build the SOP
			     SOP_VBD::myTemplateList,	// My parameters
			     0,				// Min # of sources
			     0,				// Max # of sources
			     SOP_VBD::myVariables,	// Local variables
			     OP_FLAG_GENERATOR)		// Flag it as generator
	    );
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//PUT YOUR CODE HERE
//You need to declare your parameters here
//Example to declare a variable for angle you can do like this :
static PRM_Name		angleName("angle", "Angle");
static PRM_Name		stepSizeName("stepSize", "StepSize");
static PRM_Name		iterationsName("iterations", "Iterations");
static PRM_Name		fileName("grammarFilePath", "GrammarFilePath");









//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//				     ^^^^^^^^    ^^^^^^^^^^^^^^^
//				     internal    descriptive version


// PUT YOUR CODE HERE
// You need to setup the initial/default values for your parameters here
// For example : If you are declaring the inital value for the angle parameter
static PRM_Default angleDefault(30.0);
static PRM_Default stepSizeDefault(1.0);
static PRM_Default iterationsDefault(2);
static PRM_Default fileDefault(0, "UNDEFINED");











////////////////////////////////////////////////////////////////////////////////////////

static PRM_Range angleRange(PRM_RANGE_RESTRICTED, 0, PRM_RANGE_UI, 360);
static PRM_Range stepSizeRange(PRM_RANGE_RESTRICTED, 0, PRM_RANGE_UI, 5);
static PRM_Range iterationsRange(PRM_RANGE_RESTRICTED, 0, PRM_RANGE_UI, 8);

PRM_Template
SOP_VBD::myTemplateList[] = {
// PUT YOUR CODE HERE
// You now need to fill this template with your parameter name and their default value
// EXAMPLE : For the angle parameter this is how you should add into the template
// PRM_Template(PRM_FLT,	PRM_Template::PRM_EXPORT_MIN, 1, &angleName, &angleDefault, 0),
// Similarly add all the other parameters in the template format here






/////////////////////////////////////////////////////////////////////////////////////////////

   PRM_Template(PRM_FLT, PRM_Template::PRM_EXPORT_MIN, 1, &angleName, &angleDefault, 0, &angleRange),
   PRM_Template(PRM_FLT, PRM_Template::PRM_EXPORT_MIN, 1, &stepSizeName, &stepSizeDefault, 0, &stepSizeRange),
   PRM_Template(PRM_INT, PRM_Template::PRM_EXPORT_MIN, 1, &iterationsName, &iterationsDefault, 0, &iterationsRange),
   PRM_Template(PRM_FILE, 1, &fileName, &fileDefault),
   PRM_Template()
};


// Here's how we define local variables for the SOP.
enum {
	VAR_PT,		// Point number of the star
	VAR_NPT		// Number of points in the star
};

CH_LocalVariable
SOP_VBD::myVariables[] = {
    { "PT",	VAR_PT, 0 },		// The table provides a mapping
    { "NPT",	VAR_NPT, 0 },		// from text string to integer token
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

OP_ERROR
SOP_VBD::cookMySop(OP_Context &context)
{
	fpreal currTime = context.getTime();

	// PUT YOUR CODE HERE
	// Decare the necessary variables and get always keep getting the current value in the node
	// For example to always get the current angle thats set in the node ,you need to :
	//    float angle;
	//    angle = ANGLE(now)       
    //    NOTE : ANGLE is a function that you need to use and it is declared in the header file to update your values instantly while cooking 
	float angle; angle = ANGLE(currTime);
    float stepSize; stepSize = STEP_SIZE(currTime);
    uint iterationCount; iterationCount = static_cast<uint>(ITERATIONS(currTime));
    UT_String hFilePath;
    FILE_PATH(currTime, hFilePath);
    string filePath = hFilePath.toStdString();

	///////////////////////////////////////////////////////////////////////////

	//PUT YOUR CODE HERE
	// Next you need to call your Lystem cpp functions 
	// Below is an example , you need to call the same functions based on the variables you declare
    if (filePath == "UNDEFINED") {
        std::cerr << "Current Time: " << currTime << std::endl;
        // return error();
    }

    //myplant.loadProgram(filePath);// loadProgramFromString("F\nF->F[+F]F[-F]");
    //myplant.setDefaultAngle(angle);
    //myplant.setDefaultStep(stepSize);
	





	///////////////////////////////////////////////////////////////////////////////

	// PUT YOUR CODE HERE
	// You the need call the below function for all the genrations ,so that the end points points will be
	// stored in the branches vector , you need to declare them first
	//std::vector<LSystem::Branch> branches; // array of vectors to pregenerate all iters? is that what they want? makes no sense cuz its run every cook

	//for (int i = 0; i < 2; i++) // egnerations
	//{
	//	  myplant.process(i, branches);
	//}
    //myplant.process(iterationCount, branches);





	///////////////////////////////////////////////////////////////////////////////////


	// Now that you have all the branches ,which is the start and end point of each point ,its time to render 
	// these branches into Houdini 
    

	// PUT YOUR CODE HERE
	// Declare all the necessary variables for drawing cylinders for each branch 
    float		 rad, tx, ty, tz;
    int			 divisions, plane;
    int			 xcoord =0, ycoord = 1, zcoord =2;
    float		 tmp;
    UT_Vector4		 pos;
    GU_PrimPoly		*poly;
    int			 i;
    UT_Interrupt	*boss;

    // Since we don't have inputs, we don't need to lock them.

    divisions  = glm::ceil(4+currTime);	// We need twice our divisions of points
    myTotalPoints = divisions;		// Set the NPT local variable value, TODO: NOT ACCURATE RN FOR THE USER!
    myCurrPoint   = 0;			// Initialize the PT local variable



    // Check to see that there hasn't been a critical error in cooking the SOP.
    if (error() < UT_ERROR_ABORT)
    {
	boss = UTgetInterrupt();
	if (divisions < 4)
	{
	    // With the range restriction we have on the divisions, this
	    //	is actually impossible, but it shows how to add an error
	    //	message or warning to the SOP.
	    addWarning(SOP_MESSAGE, "Invalid divisions");
	    divisions = 4;
	}
	gdp->clearAndDestroy();  // Clear all geo of this node

	// Start the interrupt server
	if (boss->opStart("Building Sim Frame"))
	{
        // PUT YOUR CODE HERE
	    // Build a polygon
	    // You need to build your cylinders inside Houdini from here
		// TIPS:
		// Use GU_PrimPoly poly = GU_PrimPoly::build(see what values it can take)
		// Also use GA_Offset ptoff = poly->getPointOffset()
		// and gdp->setPos3(ptoff,YOUR_POSITION_VECTOR) to build geometry.

        GA_Offset ptoff;
        UT_Vector3 tpos;

        vec3 sta = vec3(0);// b.first;
        vec3 end = vec3(1, 0, 0);// b.second;
        vec3 mid = 0.5f * (sta + end);

        UT_Vector3 hmid = UT_Vector3(mid[0], mid[1], mid[2]);

        vec3 baseDir = vec3(0.0f, 1.0f, 0.0f);
        vec3 dir = end - sta;
        float length = glm::length(dir);
        dir = normalize(dir);
        vec3 axis = normalize(glm::cross(baseDir, dir));
        float angle = std::acos(glm::dot(dir, baseDir));

        UT_Matrix4 transform;
        transform.identity();
        transform.scale(1.0f, length, 1.0f, 1.0f);
        transform.rotate(UT_Vector3(axis[0], axis[1], axis[2]), angle);
        transform.translate(hmid);

        float radius = 0.5f;
        int cols = divisions;// 4;

        // set pos of vertices
        for (int i = 0; i < cols; i++)
        {
            float angle = (float)i / cols * M_PI * 2.0f;
            tpos.x() = radius * cos(angle);
            tpos.y() = 0.5f;
            tpos.z() = radius * sin(angle);
            tpos = tpos * transform;

            ptoff = gdp->appendPointOffset();
            // std::cerr << ptoff << std::endl;
            gdp->setPos3(ptoff, tpos);
        }

        for (int i = 0; i < cols; i++)
        {
            float angle = (float)i / cols * M_PI * 2.0f;
            tpos.x() = radius * cos(angle);
            tpos.y() = -0.5f;
            tpos.z() = radius * sin(angle);
            tpos = tpos * transform;

            ptoff = gdp->appendPointOffset();
            // std::cerr << ptoff << std::endl;
            gdp->setPos3(ptoff, tpos);
        }

        for (int i = 0; i < cols; i++)
        {
            int off = 0 * (2 * cols + 4 * cols); // + 4 * cols because for each point offset we give for some reason appendPointOffset gets larger?!??!?!?!?!?!!?!?
            int next = (i + 1) % cols;

            int lsm[4] = { off + i, off + next,off + cols + next,off + cols + i };

            GU_PrimPoly* poly = GU_PrimPoly::build(gdp, 4, GU_POLY_CLOSED);

            poly->setPointOffset(0, lsm[0]);// off + i);
            poly->setPointOffset(1, lsm[1]);// off + next);
            poly->setPointOffset(2, lsm[2]);// off + cols + next);
            poly->setPointOffset(3, lsm[3]);// off + cols + i);

            // std::cerr << lsm[0] << ", " << lsm[1] << ", " << lsm[2] << ", " << lsm[3] << std::endl;
        }

		////////////////////////////////////////////////////////////////////////////////////////////

	    // Highlight the star which we have just generated.  This routine
	    // call clears any currently highlighted geometry, and then it
	    // highlights every primitive for this SOP. 
	    select(GU_SPrimitive);

        // This flag gets reset to false every beginning cook, so we have to set back to true
        flags().setTimeDep(true);
	}

	// Tell the interrupt server that we've completed. Must do this
	// regardless of what opStart() returns.
	boss->opEnd();
    }

    myCurrPoint = -1;
    return error();
}

