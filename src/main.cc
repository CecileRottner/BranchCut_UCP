#include <ilcplex/ilocplex.h>
#include <ilcplex/ilocplexi.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <string>
#include <math.h>
#include <ctime>

#include "InstanceUCP.h"
#include "ModeleUCP.h"
#include "Separation.h"

using namespace std ;

ILOSTLBEGIN


// USERCUT : Cuts added with addUserCut must be real cuts, in that the solution of a MIP does not depend on whether the cuts are added or not. Instead, they are there only to strengthen the formulation. 
// LAZYCONSTRAINT :  It is an error to use addUserCut for lazy constraints, that is, constraints whose absence may potentially change the solution of the problem. 
//Use addLazyConstraints or, equivalently, addCut when you add such a constraint. 

ILOUSERCUTCALLBACK6(CtCallback, Separation &, sep, IloNumArray &, solPLx, IloNumArray &, solPLu, IloBoolVarArray &, x, IloBoolVarArray &, u, IloNumArray &, data) {


    IloInt T = sep.getT() ;

    IloInt n = sep.getn() ;

    getValues(solPLx, x) ;
    getValues(solPLu, u) ;



    clock_t start;
    start = clock();

    int nbCuts = 0 ;



    // calcul du nombre de coupes ajoutées par Cplex
    nbCuts += getNcuts(IloCplex::CutCover);
    nbCuts += getNcuts(IloCplex::CutGubCover);
    nbCuts += getNcuts(IloCplex::CutFlowCover);
    nbCuts += getNcuts(IloCplex::CutClique);
    nbCuts += getNcuts(IloCplex::CutFrac);
    nbCuts += getNcuts(IloCplex::CutMir);
    nbCuts += getNcuts(IloCplex::CutMCF);
    nbCuts += getNcuts(IloCplex::CutFlowPath);
    nbCuts += getNcuts(IloCplex::CutDisj);
    nbCuts += getNcuts(IloCplex::CutImplBd);
    nbCuts += getNcuts(IloCplex::CutZeroHalf);
    nbCuts += getNcuts(IloCplex::CutLiftProj);


    nbCuts -= data[2] ;
    data[2] += nbCuts ;

    int methode = data[4] ;
    //methode = 0 : pas de séparation
    // methode = 1 : separation des up-sets (classique)
    // methode = 2 : separation up-set classique + interval-up-set (!CoverFound et isAfterCutLoop, begin=1)

    //methode = 5 : separation up-set heuristique de Letchford
    //methode = 7 : separation up-set classique + interval-up-set (!CoverFound, begin=2)
    //methode = 8 : separation up-set classique + interval-up-set (!CoverFound, begin=2) (ajoutées pour un seul i par intervalle)



    /*  if ((data[4] >= 2) && (data[4] <= 4)) {
        suppPartielle=1 ;
    }*/

    // Test Séparation heuristique

    bool inCallback = 1;
    /*if (data[3] > 90 && data[4] == 6) {
        inCallback = 0 ;
    }*/
    /*int nbUserCuts;
        if (data[4]==5) {
                nbUserCuts=200 ;
        }
        if (data[4]==6) {
                nbUserCuts=300;
        }
        if (data[4]==7) {
                nbUserCuts=100000;
        }*/


    int CoverFound=0;

    ////////////////////  Séparation des up-sets //////////////////////////
    if ( (methode > 0) && inCallback && isAfterCutLoop() && data[3] < 300 ) {
        //cout <<"Separation up-set: " << data[3] << endl ;

        int i = 0 ;
        for (int t0=0 ; (t0 < T) ; t0++) {

            sep.computeCosts(t0, t0, solPLx, solPLu) ;
            sep.computeWeightedCosts(t0, t0, solPLx, solPLu) ;

            IloRange cons = IloRange(sep.env, -IloInfinity, + IloInfinity) ;

            if (methode != 5 && methode != 6) {

                sep.Separe(cons, i, t0, t0, 1) ;

            }

            if (data[4] == 5 || data[4] == 6) {
                sep.Separe(cons, i, t0, t0, 5) ;
            }


            if (cons.getLB() >= 0) {
                IloExpr exprCons = cons.getExpr() ;
                IloInt alpha = cons.getLB() ;
                double violation = getValue(exprCons) - alpha ;

                try {

                    add(cons, IloCplex::UseCutForce) ;
                    CoverFound=1;
                    data[3] ++ ;
                    // cout << "inégalité up-set ajoutée, itération: " << getNiterations() << ", t0: " << t0 << ", violation : " << violation <<endl ;

                    cons.end() ;

                }
                catch (...) {
                    throw;
                }

            }

        }
    }

    //////////////////////////////////////////////////////////////////////////



    //Séparation Interval-Cover si aucune cover n'est trouvée dans ce callback

    int separeIC = !CoverFound;

    double intGAP = getMIPRelativeGap() ;

    /*if ( data[6] < 2 ) {
        separeIC = 0 ;
    }*/

    int suppPartielle = 0 ;

    if (methode!=7 && methode!=8) {
        separeIC = separeIC && isAfterCutLoop() ;
    }

    if ( (methode > 1) && (methode != 5) && inCallback && getNnodes()==0 && separeIC  ) {

        IloIntArray units(sep.env, n) ;
        for (int j=0 ; j <n ; j++) {
            units[j] = j;
        }

        sep.saveIndices() ;


        double lastViolationVal = 0 ;
        double currentViolationVal = 0 ;

        int begin = 1 ;
        if (methode == 7 || methode == 8) {
            begin=2 ;
        }



        for (int t0=0 ; (t0 < T) ; t0++) {


            for (int t1 = t0+begin ; (t1 <= fmin(T-1, t0+sep.Lmax)) ; t1++) {
                //for (int t1 = t0+1 ; t1 <= fmin(T-1, t0+1) ; t1++) { // on ne cherche que les intervalles de taille 2

                //cout << "--------------------------- [t0, t1] : [" << t0 << ", " << t1 << "]-------------------------" << endl ;

                //Mise à jour des coûts
                sep.computeCosts(t0, t1, solPLx, solPLu) ;
                sep.computeWeightedCosts(t0, t1, solPLx, solPLu) ;

                if ( data[4] == 8 ) {
                    sep.permut(units) ;
                }

                int stop=0 ;

                for (int i0 = 0 ; (i0 < n) && (!stop) ; i0++) {

                    int i = units[i0] ;

                    if (sep.iOK(i, t0, t1)) {
                        IloRange cons = IloRange(sep.env, -IloInfinity, + IloInfinity) ;

                        sep.Separe(cons, i, t0, t1, 2) ;


                        if (cons.getLB() > 0 ) {
                            if ( data[4] == 8 ) {
                                stop=1 ;
                            }

                            IloExpr exprCons = cons.getExpr() ;
                            IloInt alpha = cons.getLB() ;
                            currentViolationVal = getValue(exprCons) - alpha ;
                            /*cout << "Violation : " << currentViolationVal << endl ;
                            cout << "i: " << i << endl ;*/
                            if ( (lastViolationVal - currentViolationVal == 0) && suppPartielle) {
                                cout << "non ajoutée" << endl;
                                cout << endl ;
                                cons.end() ;
                            }

                            else {
                                lastViolationVal = currentViolationVal ;


                                try {

                                    add(cons, IloCplex::UseCutForce) ;
                                    cout << "recherche u plus grand que necessaire"<< endl ;
                                    for (int i=0 ; i < n ; i++) {
                                        for (int t=1 ; t < T ; t++) {
                                            if ( (solPLu[i*T+t] > 0.000001 + solPLx[i*T+t] - solPLx[i*T+t-1]  && solPLu[i*T+t] > 0.000001  ) || (i==6 && t==21) ) {
                                                cout << "unité :" << i << ", temps " << t << endl;
                                                cout << "u_t: " <<solPLu[i*T+t] << endl ;
                                                cout << "x_t - x_t-1: " <<solPLx[i*T+t] - solPLx[i*T+t-1]<< endl ;

                                            }
                                        }
                                    }
                                    if (t0==31) {
                                        for (int t = t0 ; t <= t1 ; t++) {
                                            double sum=0 ;
                                            for (int i=0 ; i <n ; i++) {
                                                sum += solPLx[i*T+t] ;
                                                // cout << "t=21 ... unité "<< i << ": x = " << solPLx[i*T+21] << ", u = " << solPLu[i*T+21] << endl ;
                                            }
                                            cout << "sum at " << t << ": " << sum << endl ;
                                        }

                                        for (int i=0 ; i <n ; i++) {

                                            cout << "unité " << i << endl;
                                            for (int t = t0 ; t <= t1 ; t++) {

                                                cout << "t="<< t << ": x = " << solPLx[i*T+t]  << ", u = " << solPLu[i*T+t] << endl ;
                                            }
                                        }
                                    }
                                    cout << "inégalité Interval-Up-Set ajoutée, itération: " << getNiterations() << ", GAP: " << 100*getMIPRelativeGap() << ", violation : " << currentViolationVal << endl ;
                                    cout << "[" << t0 << ", " << t1 << "]" << endl ;
                                    cout << endl ;
                                    data[5] ++ ;
                                    cons.end() ;
                                }
                                catch (...) {
                                    cout << "ici" << endl ;
                                    throw;
                                }
                                // cout << endl ;
                            }
                        }

                        /*else {
                            int notOK = 0 ;
                            while (notOK < 1000) {
                                sep.Separe(cons, i, t0, t1, 4) ;

                                if (cons.getLB() > 0 ) {
                                    notOK=1000 ;
                                    IloExpr exprCons = cons.getExpr() ;
                                    IloInt alpha = cons.getLB() ;
                                    currentViolationVal = getValue(exprCons) - alpha ;
                                    cout << "Violation : " << currentViolationVal << endl ;

                                }

                                else notOK++ ;
                            }
                        }*/
                    }

                }


            }

        }


        sep.copyIndices();

    }

    data[0] += ( clock() - start ) / (double) CLOCKS_PER_SEC;



    if (getNnodes() == 0) {
        data[1] = getBestObjValue()  ;
    }

    // cout << "Cplex cuts: " << nbCuts << endl ;




}

string to_string_us(int number){
    string number_string = "";
    char ones_char = '0';
    int ones = 0;
    while(true){
        ones = number % 10;
        switch(ones){
        case 0: ones_char = '0'; break;
        case 1: ones_char = '1'; break;
        case 2: ones_char = '2'; break;
        case 3: ones_char = '3'; break;
        case 4: ones_char = '4'; break;
        case 5: ones_char = '5'; break;
        case 6: ones_char = '6'; break;
        case 7: ones_char = '7'; break;
        case 8: ones_char = '8'; break;
        case 9: ones_char = '9'; break;
        default : ; //ErrorHandling("Trouble converting number to string.");
        }
        number -= ones;
        number_string = ones_char + number_string;
        if(number == 0){
            break;
        }
        number = number/10;
    }
    return number_string;
}

string createName(int n, int T, int bloc, int demande, int symetrie, int cat, int id) {
    string s_n = to_string_us(n) + "_" ;

    string s_T = s_n + to_string_us(T) ;
    s_T = s_T + "_" ;

    string s_bloc = s_T + to_string_us(bloc);
    s_bloc = s_bloc + "_" ;

    string s_demande = s_bloc + to_string_us(demande) ;
    s_demande = s_demande + "_" ;

    string  s_symetrie = s_demande + to_string_us(symetrie) ;
    s_symetrie = s_symetrie + "_" ;

    string s_cat = s_symetrie + to_string_us(cat) ;
    s_cat = s_cat + "_" ;

    string s_id = s_cat + to_string_us(id) ;

    return s_id ;
}



int process(int n, int T, int bloc, int demande, int symetrie, int cat, int id, ofstream & fichier, double & time, string localisation, int IC) {

    string nom = createName(n, T, bloc, demande, symetrie, cat, id) ;
    string fileI = localisation + nom;
    string fileS = fileI + ".txt" ;

    cout << fileS << endl ;

    const char* file = fileS.c_str() ;

    IloEnv env ;
    InstanceUCP* inst = new InstanceUCP(env, file) ;

    IloBoolVarArray x(env,n*T);
    IloBoolVarArray u(env,n*T);

    IloNumArray solPLx(env, n*T);
    IloNumArray solPLu(env, n*T);

    ModeleUCP modeleUCP(inst,env,x,u) ;

    modeleUCP.defineModel() ;


    IloCplex cplex(modeleUCP.model) ;

    IloNum eps = cplex.getParam(IloCplex::Param::Simplex::Tolerances::Feasibility) ;

    IloNumArray data(env, 7) ;
    data[0]=0 ;
    data[1]=0 ;
    data[2]=0 ;
    data[3]=0 ;
    data[4]=IC ;
    data[5] = 0;
    data[6] = IC ;

    //data[0] = duration
    //data[1] = relax
    //data[2] = temps de séparation
    //data[3] = nbCutUPSET
    //data[4] = cplex ou nous
    //data[5] = nbCutINTERVAL
    //data[6] = gap d'intégrité quand les IUS sont séparées pour la première fois

    Separation sep(env, inst, x, u, eps) ;


    /*IloIntArray C(env,6) ;
        C[0] = 0;
        C[1] = 1 ;
        C[2] = 2 ;
        C[3] = 3 ;
        C[4] = 4 ;
        C[5] = 5 ;

        int alpha =  sep.ComputeAlpha(12, C, 569);
        cout << "Alpha ; " << alpha << endl ;*/

    cplex.use(CtCallback(env, sep, solPLx, solPLu, x, u, data)) ;

    //Paramètres
    cplex.setParam(IloCplex::Param::ClockType, 1); //1 : CPU TIME
    cplex.setParam(IloCplex::Param::Threads, 1);
    cplex.setParam(IloCplex::EpGap, 0.0000001) ;
    cplex.setParam(IloCplex::Param::TimeLimit, 3600) ;
    //cplex.setParam(IloCplex::NodeLim, 100);

    //cplex.setParam(IloCplex::EachCutLim, 10);

    //CUTS
    /*cplex.setParam(IloCplex::Cliques,-1); // ne change rien si on l'enlève
        // cplex.setParam(IloCplex::Covers,-1); //
        cplex.setParam(IloCplex::DisjCuts,-1);
        cplex.setParam(IloCplex::FlowCovers,-1);
        cplex.setParam(IloCplex::FlowPaths,-1);
        cplex.setParam(IloCplex::FracCuts,-1); //
        cplex.setParam(IloCplex::GUBCovers,-1);
        cplex.setParam(IloCplex::ImplBd,-1);
        cplex.setParam(IloCplex::MIRCuts,-1); //
        cplex.setParam(IloCplex::ZeroHalfCuts,-1); //*/


    //cplex.setParam(IloCplex::MCFCuts,-1);
    //cplex.setParam(IloCplex::MIPInterval,1);

    //désactivation de l'heuristique primale
    /* cplex.setParam(IloCplex::HeurFreq,-1);
        cplex.setParam(IloCplex::RINSHeur,-1);
        cplex.setParam(IloCplex::FPHeur,-1);*/


    //Résolution et affichage de la solution

    cplex.solve();
    /*cplex.getValues(solPLx, x) ;
        cplex.getValues(solPLu, u) ;

        cout << "x = " << solPLx << endl ;
        cout << "u = " << solPLu << endl ;*/

    double t = cplex.getCplexTime();
    //fichier << "" << endl ;

    double timeSep=data[0] ;
    data[0] = data[0]/(t-time) ;
    data[0] = 1000*100*data[0] ;
    data[0] = round(data[0]);
    data[0] = data[0] / 1000 ;

    IloNum rootGap = (cplex.getObjValue() - data[1])/cplex.getObjValue() ;

    //IloInt ite = cplex.getNiterations() ;


    int nbCuts = 0 ;
    nbCuts += cplex.getNcuts(IloCplex::CutCover);
    nbCuts += cplex.getNcuts(IloCplex::CutGubCover);
    nbCuts += cplex.getNcuts(IloCplex::CutFlowCover);
    nbCuts += cplex.getNcuts(IloCplex::CutClique);
    nbCuts += cplex.getNcuts(IloCplex::CutFrac);
    nbCuts += cplex.getNcuts(IloCplex::CutMir);
    nbCuts += cplex.getNcuts(IloCplex::CutMCF);
    nbCuts += cplex.getNcuts(IloCplex::CutFlowPath);
    nbCuts += cplex.getNcuts(IloCplex::CutDisj);
    nbCuts += cplex.getNcuts(IloCplex::CutImplBd);
    nbCuts += cplex.getNcuts(IloCplex::CutZeroHalf);
    nbCuts += cplex.getNcuts(IloCplex::CutLiftProj);

    fichier << data[4] <<  " & " << n << " & " << T  << " & " << id ;
    fichier << " & " << cplex.getObjValue()  ; //Optimal value
    fichier << " & " << data[1]; // RL
    fichier << " & " << cplex.getMIPRelativeGap() << " \\% " ; //approx gap
    fichier << " & " << data[3] ; //Up-set Cuts
    fichier << " & " << data[5] ; //Interval-up-set Cuts
    fichier <<  " & " << data[0] << " \\% " ; // Separation time
    fichier << " & " << cplex.getNnodes() ;
    fichier << " & " << nbCuts ;
    fichier << " & " << t - time ;
    fichier <<" \\\\ " << endl ;

    time = t ;

    //Destructeurs



    env.end() ;

    delete inst ;
    if (data[5] < IC*100) {
        return 1;
    }
    else return 0;

}


int
main()
{

    srand(time(NULL));
    ofstream fichier("result.txt");

    fichier << "Instance & n & T & OptVal & RootRelax & ApproxGap &  USCuts & IUSCuts & SepTime & Nodes & CplexCuts & CPU \\\\ " << endl;

    double time = 0 ;

    int T;
    int n ;
    int sym ;
    int demande ;
    int cat01 ;
    int bloc ;
    /* demande = 0 ;
    sym = 0 ;
    cat01 = 1 ;
    bloc = 2 ;



    //int process(int n, int T, int bloc, int demande, int symetrie, int cat, int id, ofstream & fichier, double & time, int IC)
     for (int id=1; id <=100 ; id++) {
        process(n, T, bloc, demande, sym, cat01, id, fichier, time, 0) ;
        process(n, T, bloc, demande, sym, cat01, id, fichier, time, 1) ;
        process(n, T, bloc, demande, sym, cat01, id, fichier, time, 2) ;
        fichier << endl ;

    }*/

    demande = 3 ;

    cat01 = 0 ;
    bloc = 1 ;
    n= 20 ;
    T= 96 ;

    string localisation ;

    localisation = "data/JOCO/Theorique/75%/";
    fichier << localisation << endl ;

    n=10 ;
    bloc=1 ;

    for (sym=0 ; sym <=1 ; sym++) {
        fichier << "Symétrie: " << sym << endl ;
        for (int id=1; id <=20 ; id++) {

            process(n, T, bloc, demande, sym, cat01, id, fichier, time, localisation, 0) ;
            process(n, T, bloc, demande, sym, cat01, id, fichier, time, localisation, 1) ;
            process(n, T, bloc, demande, sym, cat01, id, fichier, time, localisation, 2) ;
            fichier << endl ;
        }
        fichier << endl ;
    }
    fichier << endl ;





    //Littérature
    /*bloc = 1 ;

    localisation = "data/JOCO/Litt_pur/";
    n=50 ;
    fichier << "Litt_pur" << endl ;

    for (int id=1; id <=15 ; id++) {
        process(n, T, bloc, demande, sym, cat01, id, fichier, time, localisation, 0) ;
        process(n, T, bloc, demande, sym, cat01, id, fichier, time, localisation, 1) ;
        process(n, T, bloc, demande, sym, cat01, id, fichier, time, localisation, 7) ;
    }
    fichier << endl ;

    localisation = "data/JOCO/Litt_pur_Pmin/";
    n=40 ;
    fichier << "Litt_pur_Pmin" << endl ;

    for (int id=1; id <= 15 ; id++) {
        process(n, T, bloc, demande, sym, cat01, id, fichier, time, localisation, 0) ;
        process(n, T, bloc, demande, sym, cat01, id, fichier, time, localisation, 1) ;
       // process(n, T, bloc, demande, sym, cat01, id, fichier, time, localisation, 2) ;
        process(n, T, bloc, demande, sym, cat01, id, fichier, time, localisation, 7) ;
       // process(n, T, bloc, demande, sym, cat01, id, fichier, time, 8) ;
        fichier << endl ;
    }
    fichier << endl ;



    localisation = "data/JOCO/Real_pseudo01_noTAC/";
    n=10 ;
    bloc=2 ;
    fichier << "Realiste sans TAC, pseudo01" << endl ;

    for (int id=1; id <= 15 ; id++) {
        process(n, T, bloc, demande, sym, cat01, id, fichier, time, localisation, 0) ;
        process(n, T, bloc, demande, sym, cat01, id, fichier, time, localisation, 1) ;
        process(n, T, bloc, demande, sym, cat01, id, fichier, time, localisation, 7) ;
        fichier << endl ;
    }
    fichier << endl ;*/

    return 0 ;
}
