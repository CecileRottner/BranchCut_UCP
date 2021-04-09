#ifndef MODELEUCP
#define MODELEUCP

#include <ilcplex/ilocplex.h>

#include "InstanceUCP.h"

class ModeleUCP {

private :
    InstanceUCP* pb ;
    IloEnv env ;
    IloBoolVarArray x ;
    IloBoolVarArray u ;

    IloInt n , T ;




public:
    IloModel model ;
    ModeleUCP(InstanceUCP* inst, IloEnv envir, const IloBoolVarArray & xx, const IloBoolVarArray & uu) ;
	~ModeleUCP() ;
	
    void defineModel() ;
};

#endif /* MODELEUCP_INCLUDED */
