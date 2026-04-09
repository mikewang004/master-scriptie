// For linux g++ -I ./ -O2 -Wno-deprecated nematic23.cc -o nematic23
// To run ./nematic23 equil_t_088_tdot_e-3_time_24000000.txt
//***This program is very similar to nematic15.cc,   it is optimized for memory via using dynamical arrays and it does a different binning for pdf of nematic order and directors.***
// This program calculates the nematic order parameter within each grid element and does the cluster analysis and bond-bond correlations along the chains and interchains bond correlations
//additional feature: obtaining the MSID and Rg for amorphous regions. 
// // here we also calculate the g(r,z) for amorph and ordered regions
// This program is modified so that the density of amorphous and crystalline regions are calculated seperately. 
//modification Nov 2017: We modify the nematic order parameter, so that S is defined as the eigenvalue with the largest absolut value.
//This file reads data from the new version of LAMMPS 2016-2017.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include  <sys/types.h>
#include  <unistd.h>
#include <sstream>
#include <iostream>
#include<fstream>
#include <complex>
#include <algorithm>
#include <cstring>
#include <climits>
//#include <string>
// #include <nr3.h>
// #include <ran.h>
using namespace std;
/*********************************************************************************************************************/

#define  onethird   1.0/3.0
# define pi    acos(-1.)
# define pi2    2.0 * pi
# define asdf    sqrt(3.0)
#define nw  6
#define nmol 4500
#define nmo2 4500
#define nvect 1
int const Lchain=100; // number of monomers in a chain
ofstream nematicfile, clusterfile, clusterfile1;
/*********************************************************************************************************************************************************/
void orderparameter( int , double [][4], double [], double [],  double [], double [], double *,  double[]  );
void orderparameter1( int , double **, double [], double [],  double [], double [], double *, double[] ); // used for Ree and global nematic oder parameter

// ---- NEW: cluster function prototype ----
void find_clusters(
    int nx, int ny, int nz, int gridcount,
    int***  ncolor,
    double*** nn1,
    double*** nn2,
    double*** nn3
);

int &Max(int &a, int &b)
{
    return a > b ? a : b;
}

int &Min(int &a, int &b)
{
    return a <= b ? a : b;
}

/********************** 3D array ******************************************************************************/

double*** Allocate_3D_Double_Array(int x, int y, int z)
{
    double*** the_array = new double**[x];
    for(int i(0); i < x; ++i)
    {
        the_array[i] = new double*[y];

        for(int j(0); j < y; ++j)
        {
            the_array[i][j] = new double[z];

            for(int k(0); k < z; ++k)
            {
                the_array[i][j][k]= 0.;
            }
        }
    }
    return the_array;
}

void release_3D_Double_Array(double*** the_array, int x, int y, int z)
{
    for (int i = 0; i < x; ++i) 
    {
        for (int j = 0; j < y; ++j)
        {
            delete [] the_array[i][j];
        }
        delete [] the_array[i];
    }
    delete [] the_array;
}
////////////////////////////3D integer array ///////////////////////


int*** Allocate_3D_Integer_Array(int x, int y, int z)
{
    int*** the_array = new int**[x];
    for(int i(0); i < x; ++i)
    {
        the_array[i] = new int*[y];

        for(int j(0); j < y; ++j)
        {
            the_array[i][j] = new int[z];

            for(int k(0); k < z; ++k)
            {
                the_array[i][j][k]= 0.;
            }
        }
    }
    return the_array;
}

void release_3D_Integer_Array(int*** the_array, int x, int y, int z)
{
    for (int i = 0; i < x; ++i) 
    {
        for (int j = 0; j < y; ++j)
        {
            delete [] the_array[i][j];
        }
        delete [] the_array[i];
    }
    delete [] the_array;
}



/**********************************************2D array********************************************/
double** Allocate_2D_Double_Array(int x, int y)
{
    double** the_array = new double*[x];
    for(int i(0); i < x; ++i)
    {
        the_array[i] = new double[y];
    }
    return the_array;
}

void release_2D_Double_Array(double** the_array, int x, int y)
{
    for (int i = 0; i < x; ++i) 
    {
        delete [] the_array[i];
    }
    delete [] the_array;
}


//////////////////////////////////////

int** Allocate_2D_Integer_Array(int x, int y)
{
    int** the_array = new int*[x];
    for(int i(0); i < x; ++i)
    {
        the_array[i] = new int[y];
    }
    return the_array;
}

void release_2D_Integer_Array(int** the_array, int x, int y)
{
    for (int i = 0; i < x; ++i) 
    {
        delete [] the_array[i];
    }
    delete [] the_array;
}

////////////////////////////////////////////////////////////////////

int main( int argc, const char* argv[] )
{
    ifstream xyzfile, bondfile;

    double S; double director[3];
    int i, id, type, nsequenc, j,l,k;
    int Natomtype, Nbondtype, Nangletype;
    int timestep;
    int Nmon,Nbond,Nangle;
    int Nchain; //number of chains

    string str,str1,str2,str3,str4,str5,str6,str7,str8,str9,str10;
    double xlo, xhi, ylo, yhi, zlo, zhi,Lx,Ly,Lz ;
    int index, id1,id2, molid, ix,iy,iz;
    int m,n,Natom;
    int Nmax=4330000;
    int Mchain=1000;
    double xx,yy, zz;

    double xp1[5000], xp2[5000],  yp1[5000], yp2[5000],  zp1[5000], zp2[5000];

    double *xp, *yp,*zp,*xu, *yu,*zu;

    xp=new double [Nmax];
    yp=new double[Nmax];
    zp=new double [Nmax];
    xu=new double [Nmax];
    yu=new double[Nmax];
    zu=new double [Nmax];

    double cnd[4];

    cout << "New long array created successfully" << endl;

    cout << "argv[1] = "<< argv[1] << "   its length is" <<  strlen(argv[1]) << endl;

    const int fnamelength=strlen(argv[1]);

    char filename[fnamelength];  
    for(i=0; i<= fnamelength; i++)
        filename[i] = argv[1][i];

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    xyzfile.open(filename); 
    if(!xyzfile) 
    {
        cerr << "Error: file atoms could not be opened" << endl;
        exit(1);
    }
    getline (xyzfile ,str);
    cout << str<< endl;
    xyzfile >> Nmon >> str1   ;
    cout << "Nmon: " << Nmon << str1<< endl;
    xyzfile >> Natomtype >>  str1 >> str2 ;
    xyzfile >> Nbond >> str1   ;
    cout << Nbond << str1<< endl;
    xyzfile >> Nbondtype >>  str1 >> str2 ;
    xyzfile >> Nangle >> str1   ;
    cout << Nangle << str1<< endl;
    xyzfile >> Nangletype >>  str1 >> str2 ;

    xyzfile >> xlo >> xhi >> str1 >> str2 ;
    xyzfile >> ylo >> yhi  >> str1 >> str2 ;
    xyzfile >> zlo >> zhi >> str1 >> str2 ;
    cout << "xlo =" << ylo << " xhi = " << xhi << endl;
    cout << "ylo =" << ylo << " yhi = " << yhi << endl;
    cout << "zlo="<< zlo<<  " zhi=" << zhi << endl;
    Lx=xhi-xlo;
    Ly=yhi-ylo;
    Lz=zhi-zlo;  
    double   Lxhalf=Lx/2.;
    double   Lyhalf=Ly/2.;
    double   Lzhalf=Lz/2.;

    xyzfile >>  str1   ;

    // Read and discard the header line:
    // "ITEM: ATOMS id mol xu yu zu"
    std::string line;
    std::getline(xyzfile, line);

    // Now read Nmon atoms: id, molid, xu, yu, zu
    for (i = 1; i <= Nmon; i++) {
        xyzfile >> id >> molid >> xu[id] >> yu[id] >> zu[id];
        xp[id] = xu[id] - xlo;
        yp[id] = yu[id] - ylo;
        zp[id] = zu[id] - zlo;
    }

    xyzfile.close();
    cout << "reading of data file was made successfully" << endl;

    int max;
    Nchain=Nmon/Lchain;

    double rho;
    rho=Nmon/(Lx*Ly*Lz);

    cout<< "Nchain=" << Nchain << "rho=" << rho << endl;

    ///////////////////defining the bond vector//////////////////////////////////
    Nbond=Nchain*(Lchain-1);
    int nx,ny,nz,n1,n2,n3;
    double lx,ly,lz,Xl,Yl,Zl;

    double *lbond;
    lbond=new double [Nbond+1];
    int m1,m2,m3;

    double lambda[4], eevp[4],   eevm[4],  eevo[4], lambdasara[4]; 
    lambda[0]=0;
    eevp[0]=0;
    eevm[0]=0;
    eevo[0]=0;

    double *xm, *ym,*zm;

    xm=new double [Nbond+1]; // bond coordinates
    ym=new double[Nbond+1]; 
    zm=new double [Nbond+1]; 

    double **bond; //bond[Nbond+1][3]
    bond= new double*[Nbond+1 ];
    for(int i = 0; i < Nbond+1; ++i){
        bond[i] = new double [3];
        if (bond[i] == NULL) {cerr << " bond[i] Allocation problem!"; exit(1);}
    }

    cout << "Long array was created successfully" << endl;

    int count[1000], nematiccount[1000];

    std::vector<int> atom1_of_bond(Nbond + 1);
    std::vector<int> atom2_of_bond(Nbond + 1);

    for(n=0;n<1000;n++){
        count[n]=0;
        nematiccount[n]=0;
    }

    for (int n = 1; n <= Nchain; ++n) { 
        for (int i = 1; i < Lchain; ++i) { 
            int atom1 = (n - 1) * Lchain + i;
            int atom2 = atom1 + 1;
            int l     = (n - 1) * (Lchain - 1) + i;

            double bx = xu[atom2] - xu[atom1];
            double by = yu[atom2] - yu[atom1];
            double bz = zu[atom2] - zu[atom1];

            bond[l][0] = bx;
            bond[l][1] = by;
            bond[l][2] = bz;

            xm[l] = (xp[atom1] + xp[atom2]) / 2.0;
            ym[l] = (yp[atom1] + yp[atom2]) / 2.0;
            zm[l] = (zp[atom1] + zp[atom2]) / 2.0;

            double len = std::sqrt(bx*bx + by*by + bz*bz);
            bond[l][0] /= len;
            bond[l][1] /= len;
            bond[l][2] /= len;

            atom1_of_bond[l] = atom1;
            atom2_of_bond[l] = atom2;
        }
    }
    cout<< "l=" << l << endl;

    double u1[2][4];
    u1[0][0]=0;
    u1[0][1]=0;
    u1[0][2]=0;
    u1[0][3]=0;
    u1[1][0]=0;
    u1[1][1]=bond[10][0];
    u1[1][2]=bond[10][1];
    u1[1][3]=bond[10][2];

    double Dhalf2=0.25*(Lx*Lx+Ly*Ly+Lz*Lz);

    char bondvec[200];
    strcpy(bondvec,filename);
    strcat(bondvec, "_bondvecs.txt");

    // Grid assignment 

    double Lgrid;
    Lgrid=2.0;

    nx=(int) Lx/Lgrid;
    ny=(int) Ly/Lgrid;
    nz=(int) Lz/Lgrid;
    double*** SS=Allocate_3D_Double_Array(nx+1, ny+1, nz+1);
    double*** nn1=Allocate_3D_Double_Array(nx+1, ny+1, nz+1);
    double*** nn2=Allocate_3D_Double_Array(nx+1, ny+1, nz+1);
    double*** nn3=Allocate_3D_Double_Array(nx+1, ny+1, nz+1);
    int*** mm=Allocate_3D_Integer_Array(nx+1, ny+1, nz+1);
    int*** ncolor_int = Allocate_3D_Integer_Array(nx+1, ny+1, nz+1); // if you prefer int***; but below you used VLA int ncolor[nx+1][ny+1][nz+1]

    cout<< "nx=" << nx<< " ny= " << ny << " nz=  "<< nz << endl;
    lx=Lx/(nx*1.0);
    ly=Ly/(ny*1.0);
    lz=Lz/(nz*1.0);
    double vol=lx*ly*lz;
    int Ngrid=nx*ny*nz;
    cout<< "lx=" << lx<< " ly= " << ly << " lz=  "<< lz << "   Ngrids=" << nx*ny*nz << endl;
    double u[1000][4];

    ofstream bondvecfile1;
    bondvecfile1.open(bondvec);
    bondvecfile1 << "atom1 bx by bz xm ym zm nx ny nz"  << endl;



    int ncolor_arr_size_x = nx+1;
    int ncolor_arr_size_y = ny+1;
    int ncolor_arr_size_z = nz+1;

    // NOTE: your original code used: int ncolor[nx+1][ny+1][nz+1];
    // Here, I'll keep your logic but using int*** instead:
    int*** ncolor = ncolor_int;

    int idcell[10000], *amorphLabel;

    amorphLabel=new int [Nmon+1];

    cout<< " up to here fine" << endl;
    double r1,r2,r3,rr,P2Q,cosQ, P2Q1,cosQ1 , P2Q2,cosQ2;

    char nemname[200];
    strcpy(nemname,filename);
    strcat(nemname, "_nematic.dat");

    nematicfile.open(nemname);

    nematicfile << "S_glob SS  lambda_1_glob  lambda_2_glob  lambda_3_glob director_1_glob  director_2_glob  director_3_glob  "  << endl;

    char clustername[60];
    strcpy(clustername,filename);
    strcat(clustername, "_cluster.dat");
    clusterfile.open(clustername);

    char clustername1[60];
    strcpy(clustername1,filename);
    strcat(clustername1, "_clust.dat");
    clusterfile1.open(clustername1);

    max=1;
    int gridcount=0;
    int cc=0, ndomains=0;

    //**********  calculation of global nematic order parameter ***********************************//
    orderparameter1(Nbond, bond, lambda, eevp, eevm, eevo, &S, director);
    double globalNem=S;
    nematicfile << "The global nematic order parameter is: " << globalNem << endl;
    double bintheta=0.025, binS=0.02;
    int Namorph=0, tt,nS,Nss; 
    tt= (int) (pi/2.0/bintheta)+1;  cout <<"tt=" << tt << endl;
 
    double nem=0, thetaz=0, thetay=0, thetax=0;
    
    int Ns0= (int) 1.0/binS;
    Nss= Ns0*3/2+1;
    int Nbamorph=0, sumamorph=0,countS[Nss+1], counttheta1[tt+1], counttheta2[tt+1], counttheta3[tt+1];
    cout << "Nss= " << Nss << "Ns0= " << Ns0 << endl;
    for(i=0;i<=Nss;i++) countS[i]=0;
    for(i=0;i<=tt;i++) {
        counttheta1[i]=0;
        counttheta2[i]=0;
        counttheta3[i]=0;
    }
 
    int *bondamorph=new int [Nbond+1];
    int *idbond=new int [1000];

    for(i=0;i< 1000;i++) idbond[i]=0;
    for(i=0;i<=Nbond;i++)  bondamorph[i]=0;

    // Fill ncolor, SS, nn1, nn2, nn3 (as in your original code)
    for(n1=0;n1<nx;n1++)
        for(n2=0;n2<ny;n2++)
            for(n3=0;n3<nz;n3++){
                m=0;
                for(l=1; l<=Nbond;l++){
                    m1 = (int)(xm[l] / lx);
                    m2 = (int)(ym[l] / ly);
                    m3 = (int)(zm[l] / lz);

                    if( (m1==n1) && ( m2==n2) && (m3==n3) ){
                        ++m;
                        u[m][1]=bond[l][0];
                        u[m][2]=bond[l][1];
                        u[m][3]=bond[l][2];
                        idcell[m]=l/(Lchain-1)+l; 
                        idbond[m]=l;
                    }
                }
                if(m>max) max=m;
                ++count[m];
                if(m > 1) {
                    nematicfile  << "m=" << m << endl;
                    nematicfile << "n1 " << "n2 "<<" n3 "<< endl;
                    nematicfile << n1<< " " << n2 << " "<< n3 << endl;
                    orderparameter(m, u, lambda, eevp, eevm, eevo, &S,director);
                    mm[n1][n2][n3]=m;
                    SS[n1][n2][n3]=lambda[1];
                    nn1[n1][n2][n3]=eevp[1];
                    nn2[n1][n2][n3]=eevp[2];
                    nn3[n1][n2][n3]=eevp[3];
                    nem= S/binS;
                    nS= (int) nem;
                    
                    if( nS>0)   countS[nS]= countS[nS]+1;
                    else {
                        nS=abs(nS)+Ns0;
                        countS[nS]= countS[nS]+1;
                    }
                }
                if(SS[n1][n2][n3] > 0.8){ 
                    ncolor[n1][n2][n3]=1; // initial label of all the nematic grid elements
                    ++nematiccount[m];
                    ++gridcount;
                }
    }

    // ---- HERE is where you would typically call the cluster finder ----
    find_clusters(nx, ny, nz, gridcount, ncolor, nn1, nn2, nn3);

    return 0;
}//end program



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// ========================
// CLUSTER ALGORITHM START
// ========================

// Helper macro: walk union-find chain to root
#define FIND_ROOT_CLUSTER(x) while (LL[(x)] < 0) (x) = -LL[(x)]

// Helper: merge two clusters (no double-count).
inline void merge_clusters_int(int i, int j, int* LL) {
    if (i == j) return;
    if (i < j) { LL[i] += LL[j]; LL[j] = -i; }
    else        { LL[j] += LL[i]; LL[i] = -j; }
}

// Helper: compute P2Q between two director vectors
inline double compute_P2Q_cluster(double a1, double a2, double a3,
                                  double b1, double b2, double b3) {
    double cosQ = a1*b1 + a2*b2 + a3*b3;
    return 1.5*cosQ*cosQ - 0.5;
}

void find_clusters(
    int nx, int ny, int nz, int gridcount,
    int***  ncolor,
    double*** nn1,
    double*** nn2,
    double*** nn3
){
    const double SC = 0.97;

    int* LL = new int[gridcount + 2](); // enough for labels up to gridcount+1
    int label    = 1;
    int gridtest = 0;
    int n1, n2, n3;

    // ------------------------------------------------------------
    // 1D scan: first row [0][0][n3]
    // ------------------------------------------------------------
    if (ncolor[0][0][0] != 0) {
        ++gridtest;
        ++label;
        ncolor[0][0][0] = label;
        ++LL[label];
    }

    for (n3 = 1; n3 < nz; n3++) {
        if (ncolor[0][0][n3] != 0) {
            ++gridtest;

            bool has_left  = (ncolor[0][0][n3-1] != 0);
            bool conn_left = false;

            if (has_left) {
                double P2Q = compute_P2Q_cluster(
                    nn1[0][0][n3],   nn2[0][0][n3],   nn3[0][0][n3],
                    nn1[0][0][n3-1], nn2[0][0][n3-1], nn3[0][0][n3-1]);
                conn_left = (P2Q >= SC);
            }

            if (!conn_left) {
                ++label;
                ncolor[0][0][n3] = label;
                ++LL[label];
            } else {
                int i = ncolor[0][0][n3-1];
                FIND_ROOT_CLUSTER(i);
                ncolor[0][0][n3] = i;
                ++LL[i];
            }
        }
    }

    // PBC in n3 for [0][0][*]
    if (ncolor[0][0][0]    != 0 &&
        ncolor[0][0][nz-1] != 0 &&
        ncolor[0][0][0]    != ncolor[0][0][nz-1]) {

        double P2Q = compute_P2Q_cluster(
            nn1[0][0][0],    nn2[0][0][0],    nn3[0][0][0],
            nn1[0][0][nz-1], nn2[0][0][nz-1], nn3[0][0][nz-1]);

        if (P2Q >= SC) {
            int i = ncolor[0][0][0];
            int j = ncolor[0][0][nz-1];
            FIND_ROOT_CLUSTER(i); FIND_ROOT_CLUSTER(j);
            if (i != j) {
                merge_clusters_int(i, j, LL);
                int root = std::min(i, j);
                ncolor[0][0][0]    = root;
                ncolor[0][0][nz-1] = root;
            }
        }
    }

    // ------------------------------------------------------------
    // 2D scan: sheet n1 = 0
    // ------------------------------------------------------------
    for (n2 = 1; n2 < ny; n2++) {

        // leftmost [0][n2][0]
        if (ncolor[0][n2][0] != 0) {
            ++gridtest;

            bool has_n2  = (ncolor[0][n2-1][0] != 0);
            bool conn_n2 = false;

            if (has_n2) {
                double P2Q1 = compute_P2Q_cluster(
                    nn1[0][n2][0],   nn2[0][n2][0],   nn3[0][n2][0],
                    nn1[0][n2-1][0], nn2[0][n2-1][0], nn3[0][n2-1][0]);
                conn_n2 = (P2Q1 >= SC);
            }

            if (!conn_n2) {
                ++label;
                ncolor[0][n2][0] = label;
                ++LL[label];
            } else {
                int i = ncolor[0][n2-1][0];
                FIND_ROOT_CLUSTER(i);
                ncolor[0][n2][0] = i;
                ++LL[i];
            }
        }

        // interior [0][n2][n3]
        for (n3 = 1; n3 < nz; n3++) {
            if (ncolor[0][n2][n3] != 0) {
                ++gridtest;

                bool has_n3 = (ncolor[0][n2][n3-1] != 0);
                bool has_n2 = (ncolor[0][n2-1][n3]  != 0);
                double P2Q = 0, P2Q1 = 0;

                if (has_n3)
                    P2Q  = compute_P2Q_cluster(
                        nn1[0][n2][n3],   nn2[0][n2][n3],   nn3[0][n2][n3],
                        nn1[0][n2][n3-1], nn2[0][n2][n3-1], nn3[0][n2][n3-1]);
                if (has_n2)
                    P2Q1 = compute_P2Q_cluster(
                        nn1[0][n2][n3],   nn2[0][n2][n3],   nn3[0][n2][n3],
                        nn1[0][n2-1][n3], nn2[0][n2-1][n3], nn3[0][n2-1][n3]);

                bool conn_n3 = has_n3 && (P2Q  >= SC);
                bool conn_n2 = has_n2 && (P2Q1 >= SC);

                if (!conn_n3 && !conn_n2) {
                    ++label;
                    ncolor[0][n2][n3] = label;
                    ++LL[label];
                } else {
                    int i = -1, j = -1;
                    if (conn_n3) { i = ncolor[0][n2][n3-1]; FIND_ROOT_CLUSTER(i); }
                    if (conn_n2) { j = ncolor[0][n2-1][n3]; FIND_ROOT_CLUSTER(j); }

                    int root = INT_MAX;
                    if (conn_n3) root = std::min(root, i);
                    if (conn_n2) root = std::min(root, j);

                    ncolor[0][n2][n3] = root;
                    ++LL[root];

                    if (conn_n3) merge_clusters_int(root, i, LL);
                    if (conn_n2) merge_clusters_int(root, j, LL);
                }
            }
        }

        // PBC: n3 for [0][n2][*]
        if (ncolor[0][n2][0]    != 0 &&
            ncolor[0][n2][nz-1] != 0 &&
            ncolor[0][n2][0]    != ncolor[0][n2][nz-1]) {

            double P2Q = compute_P2Q_cluster(
                nn1[0][n2][0],    nn2[0][n2][0],    nn3[0][n2][0],
                nn1[0][n2][nz-1], nn2[0][n2][nz-1], nn3[0][n2][nz-1]);

            if (P2Q >= SC) {
                int i = ncolor[0][n2][0];
                int j = ncolor[0][n2][nz-1];
                FIND_ROOT_CLUSTER(i); FIND_ROOT_CLUSTER(j);
                if (i != j) {
                    merge_clusters_int(i, j, LL);
                    int root = std::min(i, j);
                    ncolor[0][n2][0]    = root;
                    ncolor[0][n2][nz-1] = root;
                }
            }
        }
    }

    // PBC in n2 for n1=0
    for (n3 = 0; n3 < nz; n3++) {
        if (ncolor[0][0][n3]    != 0 &&
            ncolor[0][ny-1][n3] != 0 &&
            ncolor[0][0][n3]    != ncolor[0][ny-1][n3]) {

            double P2Q = compute_P2Q_cluster(
                nn1[0][0][n3],    nn2[0][0][n3],    nn3[0][0][n3],
                nn1[0][ny-1][n3], nn2[0][ny-1][n3], nn3[0][ny-1][n3]);

            if (P2Q >= SC) {
                int i = ncolor[0][0][n3];
                int j = ncolor[0][ny-1][n3];
                FIND_ROOT_CLUSTER(i); FIND_ROOT_CLUSTER(j);
                if (i != j) {
                    merge_clusters_int(i, j, LL);
                    int root = std::min(i, j);
                    ncolor[0][0][n3]    = root;
                    ncolor[0][ny-1][n3] = root;
                }
            }
        }
    }

    // ------------------------------------------------------------
    // 3D scan: sheets n1 = 1..nx-1
    // ------------------------------------------------------------
    for (n1 = 1; n1 < nx; n1++) {

        // first site [n1][0][0], neighbor [n1-1][0][0]
        if (ncolor[n1][0][0] != 0) {
            ++gridtest;

            bool has_n1 = (ncolor[n1-1][0][0] != 0);
            bool conn_n1 = false;

            if (has_n1) {
                double P2Q2 = compute_P2Q_cluster(
                    nn1[n1][0][0],   nn2[n1][0][0],   nn3[n1][0][0],
                    nn1[n1-1][0][0], nn2[n1-1][0][0], nn3[n1-1][0][0]);
                conn_n1 = (P2Q2 >= SC);
            }

            if (!conn_n1) {
                ++label;
                ncolor[n1][0][0] = label;
                ++LL[label];
            } else {
                int i = ncolor[n1-1][0][0];
                FIND_ROOT_CLUSTER(i);
                ncolor[n1][0][0] = i;
                ++LL[i];
            }
        }

        // first row [n1][0][n3], neighbors: left and behind
        for (n3 = 1; n3 < nz; n3++) {
            if (ncolor[n1][0][n3] != 0) {
                ++gridtest;

                bool has_n3 = (ncolor[n1][0][n3-1] != 0);
                bool has_n1 = (ncolor[n1-1][0][n3] != 0);
                double P2Q = 0, P2Q2 = 0;

                if (has_n3)
                    P2Q = compute_P2Q_cluster(
                        nn1[n1][0][n3],   nn2[n1][0][n3],   nn3[n1][0][n3],
                        nn1[n1][0][n3-1], nn2[n1][0][n3-1], nn3[n1][0][n3-1]);
                if (has_n1)
                    P2Q2 = compute_P2Q_cluster(
                        nn1[n1][0][n3],   nn2[n1][0][n3],   nn3[n1][0][n3],
                        nn1[n1-1][0][n3], nn2[n1-1][0][n3], nn3[n1-1][0][n3]);

                bool conn_n3 = has_n3 && (P2Q  >= SC);
                bool conn_n1 = has_n1 && (P2Q2 >= SC);

                if (!conn_n3 && !conn_n1) {
                    ++label;
                    ncolor[n1][0][n3] = label;
                    ++LL[label];
                } else {
                    int i = -1, k = -1;
                    if (conn_n3) { i = ncolor[n1][0][n3-1]; FIND_ROOT_CLUSTER(i); }
                    if (conn_n1) { k = ncolor[n1-1][0][n3]; FIND_ROOT_CLUSTER(k); }

                    int root = INT_MAX;
                    if (conn_n3) root = std::min(root, i);
                    if (conn_n1) root = std::min(root, k);

                    ncolor[n1][0][n3] = root;
                    ++LL[root];

                    if (conn_n3) merge_clusters_int(root, i, LL);
                    if (conn_n1) merge_clusters_int(root, k, LL);
                }
            }
        }

        // PBC in n3 for first row of sheet n1
        if (ncolor[n1][0][0]    != 0 &&
            ncolor[n1][0][nz-1] != 0 &&
            ncolor[n1][0][0]    != ncolor[n1][0][nz-1]) {

            double P2Q = compute_P2Q_cluster(
                nn1[n1][0][0],    nn2[n1][0][0],    nn3[n1][0][0],
                nn1[n1][0][nz-1], nn2[n1][0][nz-1], nn3[n1][0][nz-1]);

            if (P2Q >= SC) {
                int i = ncolor[n1][0][0];
                int j = ncolor[n1][0][nz-1];
                FIND_ROOT_CLUSTER(i); FIND_ROOT_CLUSTER(j);
                if (i != j) {
                    merge_clusters_int(i, j, LL);
                    int root = std::min(i, j);
                    ncolor[n1][0][0]    = root;
                    ncolor[n1][0][nz-1] = root;
                }
            }
        }

        // left column [n1][n2][0], neighbors: n2-1 and n1-1
        for (n2 = 1; n2 < ny; n2++) {
            if (ncolor[n1][n2][0] != 0) {
                ++gridtest;

                bool has_n2 = (ncolor[n1][n2-1][0] != 0);
                bool has_n1 = (ncolor[n1-1][n2][0] != 0);
                double P2Q1 = 0, P2Q2 = 0;

                if (has_n2)
                    P2Q1 = compute_P2Q_cluster(
                        nn1[n1][n2][0],   nn2[n1][n2][0],   nn3[n1][n2][0],
                        nn1[n1][n2-1][0], nn2[n1][n2-1][0], nn3[n1][n2-1][0]);
                if (has_n1)
                    P2Q2 = compute_P2Q_cluster(
                        nn1[n1][n2][0],   nn2[n1][n2][0],   nn3[n1][n2][0],
                        nn1[n1-1][n2][0], nn2[n1-1][n2][0], nn3[n1-1][n2][0]);

                bool conn_n2 = has_n2 && (P2Q1 >= SC);
                bool conn_n1 = has_n1 && (P2Q2 >= SC);

                if (!conn_n2 && !conn_n1) {
                    ++label;
                    ncolor[n1][n2][0] = label;
                    ++LL[label];
                } else {
                    int j = -1, k = -1;
                    if (conn_n2) { j = ncolor[n1][n2-1][0]; FIND_ROOT_CLUSTER(j); }
                    if (conn_n1) { k = ncolor[n1-1][n2][0]; FIND_ROOT_CLUSTER(k); }

                    int root = INT_MAX;
                    if (conn_n2) root = std::min(root, j);
                    if (conn_n1) root = std::min(root, k);

                    ncolor[n1][n2][0] = root;
                    ++LL[root];

                    if (conn_n2) merge_clusters_int(root, j, LL);
                    if (conn_n1) merge_clusters_int(root, k, LL);
                }
            }

            // full 3D interior [n1][n2][n3]
            for (n3 = 1; n3 < nz; n3++) {
                if (ncolor[n1][n2][n3] != 0) {
                    ++gridtest;

                    bool has_n3 = (ncolor[n1][n2][n3-1] != 0);
                    bool has_n2 = (ncolor[n1][n2-1][n3] != 0);
                    bool has_n1 = (ncolor[n1-1][n2][n3] != 0);
                    double P2Q = 0, P2Q1 = 0, P2Q2 = 0;

                    if (has_n3)
                        P2Q  = compute_P2Q_cluster(
                            nn1[n1][n2][n3],   nn2[n1][n2][n3],   nn3[n1][n2][n3],
                            nn1[n1][n2][n3-1], nn2[n1][n2][n3-1], nn3[n1][n2][n3-1]);
                    if (has_n2)
                        P2Q1 = compute_P2Q_cluster(
                            nn1[n1][n2][n3],   nn2[n1][n2][n3],   nn3[n1][n2][n3],
                            nn1[n1][n2-1][n3], nn2[n1][n2-1][n3], nn3[n1][n2-1][n3]);
                    if (has_n1)
                        P2Q2 = compute_P2Q_cluster(
                            nn1[n1][n2][n3],   nn2[n1][n2][n3],   nn3[n1][n2][n3],
                            nn1[n1-1][n2][n3], nn2[n1-1][n2][n3], nn3[n1-1][n2][n3]);

                    bool conn_n3 = has_n3 && (P2Q  >= SC);
                    bool conn_n2 = has_n2 && (P2Q1 >= SC);
                    bool conn_n1 = has_n1 && (P2Q2 >= SC);

                    if (!conn_n3 && !conn_n2 && !conn_n1) {
                        ++label;
                        ncolor[n1][n2][n3] = label;
                        ++LL[label];
                    } else {
                        int i = -1, j = -1, k = -1;
                        if (conn_n3) { i = ncolor[n1][n2][n3-1]; FIND_ROOT_CLUSTER(i); }
                        if (conn_n2) { j = ncolor[n1][n2-1][n3]; FIND_ROOT_CLUSTER(j); }
                        if (conn_n1) { k = ncolor[n1-1][n2][n3]; FIND_ROOT_CLUSTER(k); }

                        int root = INT_MAX;
                        if (conn_n3) root = std::min(root, i);
                        if (conn_n2) root = std::min(root, j);
                        if (conn_n1) root = std::min(root, k);

                        ncolor[n1][n2][n3] = root;
                        ++LL[root];

                        if (conn_n3) merge_clusters_int(root, i, LL);
                        if (conn_n2) merge_clusters_int(root, j, LL);
                        if (conn_n1) merge_clusters_int(root, k, LL);
                    }
                }
            }

            // PBC in n3 for [n1][n2][*]
            if (ncolor[n1][n2][0]    != 0 &&
                ncolor[n1][n2][nz-1] != 0 &&
                ncolor[n1][n2][0]    != ncolor[n1][n2][nz-1]) {

                double P2Q = compute_P2Q_cluster(
                    nn1[n1][n2][0],    nn2[n1][n2][0],    nn3[n1][n2][0],
                    nn1[n1][n2][nz-1], nn2[n1][n2][nz-1], nn3[n1][n2][nz-1]);

                if (P2Q >= SC) {
                    int i = ncolor[n1][n2][0];
                    int j = ncolor[n1][n2][nz-1];
                    FIND_ROOT_CLUSTER(i); FIND_ROOT_CLUSTER(j);
                    if (i != j) {
                        merge_clusters_int(i, j, LL);
                        int root = std::min(i, j);
                        ncolor[n1][n2][0]    = root;
                        ncolor[n1][n2][nz-1] = root;
                    }
                }
            }
        }

        // PBC in n2 per sheet
        for (n3 = 0; n3 < nz; n3++) {
            if (ncolor[n1][0][n3]    != 0 &&
                ncolor[n1][ny-1][n3] != 0 &&
                ncolor[n1][0][n3]    != ncolor[n1][ny-1][n3]) {

                double P2Q = compute_P2Q_cluster(
                    nn1[n1][0][n3],    nn2[n1][0][n3],    nn3[n1][0][n3],
                    nn1[n1][ny-1][n3], nn2[n1][ny-1][n3], nn3[n1][ny-1][n3]);

                if (P2Q >= SC) {
                    int i = ncolor[n1][0][n3];
                    int j = ncolor[n1][ny-1][n3];
                    FIND_ROOT_CLUSTER(i); FIND_ROOT_CLUSTER(j);
                    if (i != j) {
                        merge_clusters_int(i, j, LL);
                        int root = std::min(i, j);
                        ncolor[n1][0][n3]    = root;
                        ncolor[n1][ny-1][n3] = root;
                    }
                }
            }
        }
    }

    // ------------------------------------------------------------
    // PBC in n1 direction between sheets 0 and nx-1
    // ------------------------------------------------------------
    for (n2 = 0; n2 < ny; n2++)
    for (n3 = 0; n3 < nz; n3++) {
        if (ncolor[0][n2][n3]    != 0 &&
            ncolor[nx-1][n2][n3] != 0 &&
            ncolor[0][n2][n3]    != ncolor[nx-1][n2][n3]) {

            double P2Q = compute_P2Q_cluster(
                nn1[0][n2][n3],    nn2[0][n2][n3],    nn3[0][n2][n3],
                nn1[nx-1][n2][n3], nn2[nx-1][n2][n3], nn3[nx-1][n2][n3]);

            if (P2Q >= SC) {
                int i = ncolor[0][n2][n3];
                int j = ncolor[nx-1][n2][n3];
                FIND_ROOT_CLUSTER(i); FIND_ROOT_CLUSTER(j);
                if (i != j) {
                    merge_clusters_int(i, j, LL);
                    int root = std::min(i, j);
                    ncolor[0][n2][n3]    = root;
                    ncolor[nx-1][n2][n3] = root;
                }
            }
        }
    }

    // ------------------------------------------------------------
    // Final relabeling: resolve union-find chains in ncolor
    // ------------------------------------------------------------
    for (n1 = 0; n1 < nx; n1++)
    for (n2 = 0; n2 < ny; n2++)
    for (n3 = 0; n3 < nz; n3++) {
        if (ncolor[n1][n2][n3] != 0) {
            int r = ncolor[n1][n2][n3];
            FIND_ROOT_CLUSTER(r);
            ncolor[n1][n2][n3] = r;
        }
    }

    // ------------------------------------------------------------
    // Readout: count independent clusters and total sites
    // ------------------------------------------------------------
    int labelcount = 0;
    int total_sites = 0;

    for (int i = 1; i <= label; i++) {
        if (LL[i] > 0) {
            ++labelcount;
            total_sites += LL[i];
        }
    }

    cout << "Total crystalline grid elements:    " << total_sites  << endl;
    cout << "Independent crystalline domains:    " << labelcount   << endl;

    delete[] LL;
}

// ========================
// CLUSTER ALGORITHM END
// ========================


 void orderparameter1( int npar, double **um, double lambda[], double eevp[4],  double eevm[4], double eevo[4], double *S, double director[3] )
{
    double a11, a12, a13, a23,a22,a33,aa11,aa22 ;
    double A[3][3], e[3][3], dummy;
    double cc1,cc0,pp,qq,phia,rr,vnorm;
    int i;     
    a11=0.;
    a12=0.;
    a13=0.;
    a22=0.;
    a23=0.;
    a33=0.;

    for(i=1; i<=npar;i++){
        a11+=um[i][0]*um[i][0];
        a12+=um[i][0]*um[i][1];
        a13+=um[i][0]*um[i][2];
        a22+=um[i][1]*um[i][1];
        a23+=um[i][1]*um[i][2];
        a33+=um[i][2]*um[i][2];
    }

    a11=a11/npar;
    a12=a12/npar;
    a13=a13/npar;
    a22=a22/npar;
    a23=a23/npar;
    a33=a33/npar;

    cc1=a12*a12+a13*a13+a23*a23-a11*a22-a11*a33-a22*a33;
    cc0=a11*a22*a33+2.0*a12*a13*a23-a12*a12*a33-a13*a13*a22-a23*a23*a11;
    pp=0.25*(1.00+3.0*cc1);
    qq=(27.0*cc0+9.0*cc1+2.0)/16.0;
   
    phia=acos(qq/sqrt(pow(pp,3)))/3.;
    rr=2.0 *sqrt(pp);

    lambda[1]=rr*cos(phia);
    lambda[2]=rr*cos(phia-2.*pi2/3.);
    lambda[3]=rr*cos(phia-pi2/3.);

    aa11=a11-onethird-2.*onethird*lambda[1];
    aa22=a22-onethird-2.*onethird*lambda[1];
    eevp[1]=a12*a23-a13*aa22;
    eevp[2]=a13*a12-a23*aa11;
    eevp[3]=aa11*aa22-a12*a12;
    vnorm=sqrt(eevp[1]*eevp[1]+eevp[2]*eevp[2]+eevp[3]*eevp[3]);
    if(vnorm > 1.0e-12) {
        eevp[1]=eevp[1]/vnorm;
        eevp[2]=eevp[2]/vnorm;
        eevp[3]=eevp[3]/vnorm;
    }
    else{
        eevp[1]=0.;
        eevp[2]=0.;
        eevp[3]=0.;
    }

    aa11=a11-onethird-2.*onethird*lambda[2];
    aa22=a22-onethird-2.*onethird*lambda[2];
    eevm[1]=a12*a23-a13*aa22;
    eevm[2]=a13*a12-a23*aa11;
    eevm[3]=aa11*aa22-a12*a12;
    vnorm=sqrt(eevm[1]*eevm[1]+eevm[2]*eevm[2]+eevm[3]*eevm[3]);
    if(vnorm > 1.0e-12) {
        eevm[1]=eevm[1]/vnorm;
        eevm[2]=eevm[2]/vnorm;
        eevm[3]=eevm[3]/vnorm;
    }
    else {
        eevm[1]=0;
        eevm[2]=0;
        eevm[3]=0;
    }

    aa11=a11-onethird-2.*onethird*lambda[3];
    aa22=a22-onethird-2.*onethird*lambda[3];
    eevo[1]=a12*a23-a13*aa22;
    eevo[2]=a13*a12-a23*aa11;
    eevo[3]=aa11*aa22-a12*a12;
    vnorm=sqrt(eevo[1]*eevo[1]+eevo[2]*eevo[2]+eevo[3]*eevo[3]);

    if(vnorm > 1.0e-12) {
        eevo[1]=eevo[1]/vnorm;
        eevo[2]=eevo[2]/vnorm;
        eevo[3]=eevo[3]/vnorm;
    }
    else{
        eevo[1]=0;
        eevo[2]=0;
        eevo[3]=0;
    }

    double S_nem, S1,S2,S3,S0=0;

    S1=abs(lambda[1]);
    S2=abs(lambda[2]);
    S3=abs(lambda[3]);
    S0=S1; director[0]=eevp[1]; director[1]=eevp[2]; director[2]=eevp[3]; S_nem=lambda[1];
    if(S2 > S0) {S0=S2; director[0]=eevm[1]; director[1]=eevm[2]; director[2]=eevm[3]; S_nem=lambda[2];}
    if(S3 > S0)  {S0=S3; director[0]=eevo[1]; director[1]=eevo[2]; director[2]=eevo[3]; S_nem=lambda[3];}

    double ss=sqrt(2.0/3.0*(lambda[1]*lambda[1]+lambda[2]*lambda[2]+lambda[3]*lambda[3])); 

    nematicfile << S_nem << " " << ss << " " << lambda[1]  << " "<< lambda[2] << " "<< lambda[3] << " " << director[0]  << " "<< director[1] << " "<< director[2] <<  endl;
    nematicfile << eevp[1]<< " " << eevp[2] << " "<< eevp[3] << endl;
    nematicfile << eevm[1]<< " " << eevm[2] << " "<< eevm[3] << endl;
    nematicfile << eevo[1]<< " " << eevo[2] << " "<< eevo[3] << endl;
    nematicfile << "ss= " << ss << endl;
    *S=S_nem;
    return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void orderparameter( int npar, double um[][4], double lambda[], double eevp[4],  double eevm[4], double eevo[4], double *S, double director[3] )
{
    double a11, a12, a13, a23,a22,a33,aa11,aa22 ;
    double A[3][3], e[3][3], dummy;
    double cc1,cc0,pp,qq,phia,rr,vnorm;
    int i;     
    a11=0.;
    a12=0.;
    a13=0.;
    a22=0.;
    a23=0.;
    a33=0.;

    for(i=1; i<=npar;i++){
        a11+=um[i][1]*um[i][1];
        a12+=um[i][1]*um[i][2];
        a13+=um[i][1]*um[i][3];
        a22+=um[i][2]*um[i][2];
        a23+=um[i][2]*um[i][3];
        a33+=um[i][3]*um[i][3];
    }

    a11=a11/npar;
    a12=a12/npar;
    a13=a13/npar;
    a22=a22/npar;
    a23=a23/npar;
    a33=a33/npar;

    cc1=a12*a12+a13*a13+a23*a23-a11*a22-a11*a33-a22*a33;
    cc0=a11*a22*a33+2.0*a12*a13*a23-a12*a12*a33-a13*a13*a22-a23*a23*a11;
    pp=0.25*(1.00+3.0*cc1);
    qq=(27.0*cc0+9.0*cc1+2.0)/16.0;
   
    phia=acos(qq/sqrt(pow(pp,3)))/3.;
    rr=2.0 *sqrt(pp);

    lambda[1]=rr*cos(phia);
    lambda[2]=rr*cos(phia-2.*pi2/3.);
    lambda[3]=rr*cos(phia-pi2/3.);

    aa11=a11-onethird-2.*onethird*lambda[1];
    aa22=a22-onethird-2.*onethird*lambda[1];
    eevp[1]=a12*a23-a13*aa22;
    eevp[2]=a13*a12-a23*aa11;
    eevp[3]=aa11*aa22-a12*a12;
    vnorm=sqrt(eevp[1]*eevp[1]+eevp[2]*eevp[2]+eevp[3]*eevp[3]);
    if(vnorm > 1.0e-12) {
        eevp[1]=eevp[1]/vnorm;
        eevp[2]=eevp[2]/vnorm;
        eevp[3]=eevp[3]/vnorm;
    }
    else{
        eevp[1]=0.;
        eevp[2]=0.;
        eevp[3]=0.;
    }

    aa11=a11-onethird-2.*onethird*lambda[2];
    aa22=a22-onethird-2.*onethird*lambda[2];
    eevm[1]=a12*a23-a13*aa22;
    eevm[2]=a13*a12-a23*aa11;
    eevm[3]=aa11*aa22-a12*a12;
    vnorm=sqrt(eevm[1]*eevm[1]+eevm[2]*eevm[2]+eevm[3]*eevm[3]);
    if(vnorm > 1.0e-12) {
        eevm[1]=eevm[1]/vnorm;
        eevm[2]=eevm[2]/vnorm;
        eevm[3]=eevm[3]/vnorm;
    }
    else {
        eevm[1]=0;
        eevm[2]=0;
        eevm[3]=0;
    }

    aa11=a11-onethird-2.*onethird*lambda[3];
    aa22=a22-onethird-2.*onethird*lambda[3];
    eevo[1]=a12*a23-a13*aa22;
    eevo[2]=a13*a12-a23*aa11;
    eevo[3]=aa11*aa22-a12*a12;
    vnorm=sqrt(eevo[1]*eevo[1]+eevo[2]*eevo[2]+eevo[3]*eevo[3]);

    if(vnorm > 1.0e-12) {
        eevo[1]=eevo[1]/vnorm;
        eevo[2]=eevo[2]/vnorm;
        eevo[3]=eevo[3]/vnorm;
    }
    else{
        eevo[1]=0;
        eevo[2]=0;
        eevo[3]=0;
    }

    double S_nem, S1,S2,S3,S0=0;

    S1=abs(lambda[1]);
    S2=abs(lambda[2]);
    S3=abs(lambda[3]);
    S0=S1; director[0]=eevp[1]; director[1]=eevp[2]; director[2]=eevp[3]; S_nem=lambda[1];
    if(S2 > S0) {S0=S2; director[0]=eevm[1]; director[1]=eevm[2]; director[2]=eevm[3]; S_nem=lambda[2];}
    if(S3 > S0)  {S0=S3; director[0]=eevo[1]; director[1]=eevo[2]; director[2]=eevo[3]; S_nem=lambda[3];}
    double ss=sqrt(2.0/3.0*(lambda[1]*lambda[1]+lambda[2]*lambda[2]+lambda[3]*lambda[3])); 

    nematicfile << S_nem << " " <<ss << " " << lambda[1]  << " "<< lambda[2] << " "<< lambda[3] << " " << director[0]  << " "<< director[1] << " "<< director[2] <<  endl;
    nematicfile << eevp[1]<< " " << eevp[2] << " "<< eevp[3] << endl;
    nematicfile << eevm[1]<< " " << eevm[2] << " "<< eevm[3] << endl;
    nematicfile << eevo[1]<< " " << eevo[2] << " "<< eevo[3] << endl;
    nematicfile << lambda[1] << " " << lambda[2] << " "<< lambda[3] << endl;

    nematicfile << "ss= " << ss << endl;
    *S=S_nem;
    return;
}