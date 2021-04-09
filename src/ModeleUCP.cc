#include "ModeleUCP.h"
#include "InstanceUCP.h"

#include <ilcplex/ilocplex.h>

using namespace std ;

ModeleUCP::ModeleUCP(InstanceUCP* inst, IloEnv envir, const IloBoolVarArray & xx, const IloBoolVarArray & uu) :
    pb (inst),
    env (envir),
    x (xx),
    u (uu) {
        n = pb->getn();
        T = pb->getT() ;
        model = IloModel(env);
}

ModeleUCP::~ModeleUCP() {
    //env.end() ;
}


void ModeleUCP::defineModel() {

    IloInt t ;
    IloInt i ;
    IloInt k ;

    IloNumVarArray pp(env, n*T, 0.0, 1000);

    // Objective Function: Minimize Cost
    IloExpr cost(env) ;
    for (t=0 ; t < T ; t++) {
        for (i=0; i<n; i++) {
            cost += x[i*T + t]*pb->getcf(i) + pb->getc0(i)*u[i*T + t] + (pp[i*T + t]+pb->getP(i)*x[i*T + t])*(pb->getcp(i)) ;

        }
    }

    model.add(IloMinimize(env, cost));

    cost.end() ;

    // Conditions initiales
     for (i=0; i<n; i++) {
         model.add(u[i*T] >= x[i*T] - pb->getInit(i) ) ;
     }

     for (i=0; i<n; i++) {
             IloExpr sum(env) ;
             for (k= 0; k < pb->getl(i) ; k++) {
                 sum += u[i*T + k] ;
             }
             model.add(sum <= 1 - pb->getInit(i) ) ;
             sum.end() ;
     }

     // Min up constraints
     for (i=0; i<n; i++) {
         for (t=pb->getL(i) -1 ; t < T ; t++) {
             IloExpr sum(env) ;
             for (k= t - pb->getL(i) + 1; k <= t ; k++) {
                 sum += u[i*T + k] ;
             }
             model.add(sum <= x[i*T + t]) ;
             sum.end() ;
         }
     }


     // Min down constraints
     for (i=0; i<n; i++) {
         for (t=pb->getl(i) ; t < T ; t++) {
             IloExpr sum(env) ;
             for (k= t - pb->getl(i) + 1; k <= t ; k++) {
                 sum += u[i*T + k] ;
             }
             model.add(sum <= 1 - x[i*T + t - pb->getl(i)]) ;
             sum.end() ;
         }
     }

     //Relation entre u et x
     for (i=0; i<n; i++) {
          for (t=1 ; t < T ; t++) {
             model.add(x[i*T + t] - x[i*T + t-1] <= u[i*T + t]);
         }
     }


     //Limite de production
     for (i=0; i<n; i++) {
          for (t=0 ; t < T ; t++) {
             model.add(pp[i*T + t] <= (pb->getPmax(i)-pb->getP(i))*x[i*T + t]);
             model.add(pp[i*T + t] >= 0);
         }
     }


     //Demande
     for (t=0; t < T ; t++) {
         IloExpr Prod(env) ;
         for (i=0; i<n; i++) {
             Prod += pp[i*T + t] + pb->getP(i)*x[i*T + t];
         }
         model.add(pb->getD(t) <= Prod);
         Prod.end() ;
     }


}
