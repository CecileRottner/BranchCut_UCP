#include "InstanceUCP.h"

#include <ctime>
#include <stdlib.h> // use rand
#include <fstream> // pour lire le fichier de données

using namespace std ;


InstanceUCP::InstanceUCP(IloEnv envir, const char* file) {


    //Lecture de n et de T
    ifstream fichier(file, ios::in);

    string nom = "";
    fichier >> nom;

    while(nom!="="){
        fichier >> nom;
    }
    nom = "";
    fichier >> n;

    while(nom!="="){
        fichier >> nom;
    }
    nom = "";
    fichier >> T;

    //Lecture des données
    Init_ = IloBoolArray(env, n);
    L_ = IloIntArray(env, n);
    l_ = IloIntArray(env, n);
    P_ = IloNumArray(env, n);
    Pmax_ = IloNumArray(env, n);
    cf_ = IloNumArray(env, n);
    c0_ = IloNumArray(env, n);
    cp_ = IloNumArray(env, n);

    Init = IloBoolArray(env, n);
    L = IloIntArray(env, n);
    l = IloIntArray(env, n);
    P = IloNumArray(env, n);
    Pmax = IloNumArray(env, n);
    cf = IloNumArray(env, n);
    c0 = IloNumArray(env, n);
    cp = IloNumArray(env, n);

    D = IloNumArray(env, T);

    //Init
    while(nom!="["){
        fichier >> nom;
    }
    nom = "";
    for(IloInt j=0; j<n; j++){
        fichier >> Init_[j];
    }

    //L
    while(nom!="["){
        fichier >> nom;
    }
    nom = "";
    for(IloInt j=0; j<n; j++){
        fichier >> L_[j];
    }

    //l
    while(nom!="["){
        fichier >> nom;
    }
    nom = "";
    for(IloInt j=0; j<n; j++){
        fichier >> l_[j];
    }

    //P
    while(nom!="["){
        fichier >> nom;
    }
    nom = "";
    for(IloInt j=0; j<n; j++){
        fichier >> P_[j];
    }

    //Pmax
    while(nom!="["){
        fichier >> nom;
    }
    nom = "";
    for(IloInt j=0; j<n; j++){
        fichier >> Pmax_[j];
    }

    //cf
    while(nom!="["){
        fichier >> nom;
    }
    nom = "";
    for(IloInt j=0; j<n; j++){
        fichier >> cf_[j];
    }

    //c0
    while(nom!="["){
        fichier >> nom;
    }
    nom = "";
    for(IloInt j=0; j<n; j++){
        fichier >> c0_[j];
    }

    //cp
    while(nom!="["){
        fichier >> nom;
    }
    nom = "";
    for(IloInt j=0; j<n; j++){
        fichier >> cp_[j];
    }

    //D
    while(nom!="["){
        fichier >> nom;
    }
    nom = "";
    for(IloInt j=0; j<T; j++){
        fichier >> D[j];
    }

    C_ = IloIntArray(env, n);
    IloNumArray ordre = IloNumArray(env, n);

    for (IloInt j = 0 ; j <n ; j++ ) {
        C_[j] = j ;
        ordre[j] = 1/Pmax_[j] ;
    }

    quickSort(ordre,C_, 0, n);
    //On trie les unités par Pmax décroissantes

    for (IloInt j = 0 ; j <n ; j++ ) {
        P[j] = P_[C_[j]] ;
        Pmax[j] = Pmax_[C_[j]] ;
        L[j] = L_[C_[j]] ;
        l[j] = l_[C_[j]] ;
        c0[j] = c0_[C_[j]] ;
        cf[j] = cf_[C_[j]] ;
        cp[j] = cp_[C_[j]] ;
        Init[j] = Init_[C_[j]] ;
    }

    SommePmax = 0 ;
	for (int j = 0 ; j <n ; j++) {
		SommePmax += Pmax[j] ;
	}

    cout << "demande sur 31-36 : " ;
    for (int t=31 ; t <= 36 ; t++) {
        cout << D[t]<< " ";
    }
    cout << endl ;

    cout << "Pmax: " << Pmax << endl ;
}


InstanceUCP::~InstanceUCP() {

    Init.end() ;
    L.end() ;
    l.end() ;
    D.end() ;
    P.end() ;
    Pmax.end();
    c0.end() ;
    cf.end() ;


    Init_.end() ;
    L_.end() ;
    l_.end() ;
    P_.end() ;
    Pmax_.end();
    c0_.end() ;
    cf_.end() ;
    C_.end() ;

    env.end() ;

}

double moyenne(IloNumArray V) {
	IloInt n = V.getSize() ;
	double moyP = 0 ;
	for (IloInt j =0 ; j < n ; j++) {
		moyP += V[j];
	}
	moyP = moyP/n ;
	return moyP ;
}

double variance(IloNumArray X, IloNumArray Y) {
	IloInt n = X.getSize() ;
	double sigma = 0 ;
	double moyX = moyenne(X);
	double moyY = moyenne(Y) ;
	for (IloInt j = 0 ; j < n ; j++) {
		sigma += (X[j] - moyX)*(Y[j] - moyY);
	}
	sigma = sigma/n ;

	return sigma ;
}

IloInt InstanceUCP::getn() {
    return n ;
}

IloInt InstanceUCP::getT() {
    return T ;
}

IloBool InstanceUCP::getInit(IloInt i) {
    if ((i >= n)||(i < 0)) {
        cout << "Attention, Init[" << i << "] n'existe pas." << endl ;
    }
    return Init[i] ;
}

IloInt InstanceUCP::getL(IloInt i) {
    if ((i >= n)||(i < 0)) {
        cout << "Attention, L[" << i << "] n'existe pas." << endl ;
    }
    return L[i] ;
}

IloInt InstanceUCP::getl(IloInt i) {
    if ((i >= n)||(i < 0)) {
        cout << "Attention, l[" << i << "] n'existe pas." << endl ;
    }
    return l[i] ;
}

IloNum InstanceUCP::getD(IloInt i) {
    if ((i >= T)||(i < 0)) {
        cout << "Attention, D[" << i << "] n'existe pas." << endl ;
    }
    return D[i] ;
}

IloNum InstanceUCP::getP(IloInt i) {
    if ((i >= n)||(i < 0)) {
        cout << "Attention, P[" << i << "] n'existe pas." << endl ;
    }
    return P[i] ;
}

IloNum InstanceUCP::getPmax(IloInt i) {
    if ((i >= n)||(i < 0)) {
        cout << "Attention, Pmax[" << i << "] n'existe pas." << endl ;
    }
    return Pmax[i] ;
}

IloNum InstanceUCP::getcf(IloInt i) {
    if ((i >= n)||(i < 0)) {
        cout << "Attention, cf[" << i << "] n'existe pas." << endl ;
    }
    return cf[i] ;
}

IloNum InstanceUCP::getc0(IloInt i) {
    if ((i >= n)||(i < 0)) {
        cout << "Attention, c0[" << i << "] n'existe pas." << endl ;
    }
    return c0[i] ;
}

IloNum InstanceUCP::getcp(IloInt i) {
    if ((i >= n)||(i < 0)) {
        cout << "Attention, cp[" << i << "] n'existe pas." << endl ;
    }
    return cp[i] ;
}

IloInt InstanceUCP::getSommePmax() {
    return SommePmax ;
}

void InstanceUCP::quickSort(IloNumArray const & ordre, IloIntArray & indices, IloInt p, IloInt q) {
    IloInt r;
    if(p<q)
    {
        r=partition(ordre, indices, p,q);
        quickSort(ordre, indices,p,r);
        quickSort(ordre, indices,r+1,q);
    }
}

IloInt InstanceUCP::partition(IloNumArray const & ordre, IloIntArray & indices, IloInt p, IloInt q)
{
    IloInt x= indices[p];
    IloInt i=p;
    IloInt j;

    for(j=p+1; j<q; j++)
    {
        if(ordre[indices[j]]<= ordre[x])
        {
            i=i+1;
            swap(indices[i],indices[j]);
        }

    }

    swap(indices[i],indices[p]);
    return i;
}
