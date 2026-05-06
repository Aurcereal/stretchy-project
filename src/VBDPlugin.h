

#pragma once

//#include <GEO/GEO_Point.h>
//
#include <SOP/SOP_Node.h>
#include "VBDSolver.h"

namespace HDK_Sample {
class SOP_VBD : public SOP_Node
{
public:
    static OP_Node		*myConstructor(OP_Network*, const char *,
							    OP_Operator *);

    /// Stores the description of the interface of the SOP in Houdini.
    /// Each parm template refers to a parameter.
    static PRM_Template		 myTemplateList[];

    /// This optional data stores the list of local variables.
    static CH_LocalVariable	 myVariables[];
    // If we do evalFloat or wtver without local variable info,
    // If user types sin($PT) into a parameter, evalFloat won't know $PT and fail
    // So it needs local variable info through evalFloatInst
    // Would then need to loop over all points and separately eval 

protected:

	     SOP_VBD(OP_Network *net, const char *name, OP_Operator *op);
    virtual ~SOP_VBD();

    /// Disable parameters according to other parameters.
    virtual unsigned		 disableParms();


    /// cookMySop does the actual work of the SOP computing, in this
    /// case, a LSYSTEM
    virtual OP_ERROR		 cookMySop(OP_Context &context);

    /// This function is used to lookup local variables that you have
    /// defined specific to your SOP.
    virtual bool		 evalVariableValue(
				    fpreal &val,
				    int index,
				    int thread);
    // Add virtual overload that delegates to the super class to avoid
    // shadow warnings.
    virtual bool		 evalVariableValue(
				    UT_String &v,
				    int i,
				    int thread)
				 {
				     return evalVariableValue(v, i, thread);
				 }

private:
    SolverParams GetParams(fpreal time);
    vector<SolverParams> cachedParams;

    VBDSolver vbdSolver = VBDSolver(&cachedParams);
    GA_DataId inputGeoDataID = -1;
    GA_DataId constraintGroupID = -1;
    GA_DataId collisionGeoDataID = -1;
    GA_DataId initialVelDataID = -1;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Member variables are stored in the actual SOP, not with the geometry
    /// In this case these are just used to transfer data to the local 
    /// variable callback.
    /// Another use for local data is a cache to store expensive calculations.

	// NOTE : You can declare local variables here like this  
    int		myCurrPoint;
    int		myTotalPoints;
};
} // End HDK_Sample namespace