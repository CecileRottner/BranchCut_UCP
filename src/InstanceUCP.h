#ifndef INSTANCEUCP
#define INSTANCEUCP

#include <ilcplex/ilocplex.h>
#include <fstream>

class InstanceUCP {

public:
    IloEnv env ;
    IloInt n, T ;
    IloBoolArray Init ;
    IloIntArray L, l ;
    IloNumArray D, P, Pmax, cf, c0, cp ;

    IloBoolArray Init_ ;
    IloIntArray L_, l_ ;
    IloNumArray P_, Pmax_, cf_, c0_ , cp_;

    IloIntArray C_ ;

	IloInt SommePmax ;


public:
    InstanceUCP(IloEnv envir, const char* file) ;
    ~InstanceUCP() ;

    IloInt getn() ;
    IloInt getT() ;
    IloBool getInit(IloInt i) ;
    IloInt getL(IloInt i) ;
    IloInt getl(IloInt i) ;
    IloNum getD(IloInt i) ;
    IloNum getP(IloInt i) ;
    IloNum getPmax(IloInt i) ;
    IloNum getcf(IloInt i) ;
    IloNum getc0(IloInt i) ;
    IloNum getcp(IloInt i) ;

	IloInt getSommePmax() ;

    IloInt partition(IloNumArray const & ordre, IloIntArray & indices, IloInt p, IloInt q) ;

    void quickSort(IloNumArray const & ordre, IloIntArray & indices, IloInt p, IloInt q) ;

};

#endif /* INSTANCEUCP_INCLUDED */
