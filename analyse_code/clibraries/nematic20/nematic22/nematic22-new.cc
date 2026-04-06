// For linux g++ -I ./ -O2 -Wno-deprecated nematic22-new.cc -o nematic22-new
// To run ./nematic22-new equil_t_088_tdot_e-3_time_24000000.txt
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
{ ifstream xyzfile, bondfile;

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
//int Nmax = 720000;
//int Mchain = 100;
double xx,yy, zz;


double xp1[5000], xp2[5000],  yp1[5000], yp2[5000],  zp1[5000], zp2[5000];

double *xp, *yp,*zp,*xu, *yu,*zu;

xp=new double [Nmax]; // xp[id]



yp=new double[Nmax]; 


zp=new double [Nmax]; 


xu=new double [Nmax]; // xp[id] //unwrapped coordiantes



yu=new double[Nmax]; 


zu=new double [Nmax]; 

double cnd[4];



cout << "New long array created successfully" << endl;


cout << "argv[1] = "<< argv[1] << "   its length is" <<  strlen(argv[1]) << endl;// anouncing the name of input file for initial values

const int fnamelength=strlen(argv[1]);

char filename[fnamelength];  
for(i=0; i<= fnamelength; i++)

		filename[i] = argv[1][i];  //converting the input argument to the file name string

// // argv[1] MUST be valid; ideally ensure argc >= 2 earlier in main
// cout << "argv[1] = " << argv[1]
//      << "   its length is " << strlen(argv[1]) << endl;  // prints length

// const int fnamelength = strlen(argv[1]) + 1;  // +1 for the '\0'

// char filename[fnamelength];  // valid indices: 0 .. fnamelength-1

// // Copy all characters INCLUDING the terminating '\0'
// for (int i = 0; i < fnamelength; ++i) {
//     filename[i] = argv[1][i];
// }



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   xyzfile.open(filename); // opens the file
      
       if(!xyzfile) 
       { // file couldn't be opened
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

// for(i=1;i<=Natomtype;i++)
//  xyzfile >>  j >> l   ;

//  xyzfile >>  str1 >> str2 >> str3 >> str4  ;
// cout << str4 << endl;
// for(i=1;i<=Natomtype;i++)
// xyzfile >> xx >> yy >> zz ;

/// Bonds block
//     xyzfile >> str1;
// cout <<  str1<<  endl;

//  for(i=1;i<=Nbond; i++)
//      { xyzfile >>  j>> l >> id1 >> id2   ;
//       cout << j << " "<<  l << " "<< id1 << " " << id2 << endl;
//       //bond1[i][0]=id1;
//     // bond1[i][1]=id2; 
//     if(i==9900) cout << j << " "<<  l << " "<< id1 << " " << id2 << endl;  } 
    
/// End bond block


// Read and discard the header line:
// "ITEM: ATOMS id mol xu yu zu"
std::string line;
std::getline(xyzfile, line);   // assumes we're positioned right before this line

// Now read Nmon atoms: id, molid, xu, yu, zu
for (i = 1; i <= Nmon; i++) {
//for (i = 1; i <= 10; i ++) {
    xyzfile >> id >> molid >> xu[id] >> yu[id] >> zu[id];
    // If you still want shifted coordinates starting at (0,0,0) from xlo,ylo,zlo:
    xp[id] = xu[id] - xlo;
    yp[id] = yu[id] - ylo;
    zp[id] = zu[id] - zlo;
    //cout << xp[id] << xu[id] << xlo << endl;
}
// xyzfile >> str1;
//cout <<  str1<< endl;
// xyzfile >> str1;
//cout <<  str1<< endl;
///Velocities block 
//  xyzfile >> str1;
// cout <<  str1<<  endl;
      
//      for(i=1;i<=Nmon; i++)
//      { xyzfile >>  id >>  xx >> yy >> zz   ;
//     if(i==10000) cout << id <<  xx << " "<< yy << " " << zz << endl;  } 
/// end velocities block



 xyzfile.close();
   cout << "reading of data file was made successfully" << endl;


 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////






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
//lbond=new double[Nbond];

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
if (bond[i] == NULL) {cerr << " bond[i] Allocation problem!"; exit(1);}}


cout << "Long array was created successfully" << endl;

int count[1000], nematiccount[1000];

std::vector<int> atom1_of_bond(Nbond + 1);
std::vector<int> atom2_of_bond(Nbond + 1);

for(n=0;n<1000;n++){
count[n]=0;
nematiccount[n]=0;}


// Assumptions:
// - Atoms are ordered by chain:
//   chain 1: 1..Lchain
//   chain 2: Lchain+1..2*Lchain
//   ...
//   chain n: (n-1)*Lchain+1 .. n*Lchain
// - Nmon = Nchain * Lchain
// - Nbond = Nchain * (Lchain - 1)
// - Arrays xp,yp,zp,xu,yu,zu are valid for indices 1..Nmon
// - xm,ym,zm,bond,lbond are valid for 1..Nbond

for (int n = 1; n <= Nchain; ++n) { // n <= number of chains
    for (int i = 1; i < Lchain; ++i) { //l < length of chain 
        int atom1 = (n - 1) * Lchain + i;   // 1 .. Nmon-1
        int atom2 = atom1 + 1;              // 2 .. Nmon
        int l     = (n - 1) * (Lchain - 1) + i; // 1 .. Nbond


        double bx = xu[atom2] - xu[atom1];
        double by = yu[atom2] - yu[atom1];
        double bz = zu[atom2] - zu[atom1];

        bond[l][0] = bx;
        bond[l][1] = by;
        bond[l][2] = bz;

        xm[l] = (xp[atom1] + xp[atom2]) / 2.0;
        ym[l] = (yp[atom1] + yp[atom2]) / 2.0;
        zm[l] = (zp[atom1] + zp[atom2]) / 2.0;
        //cout << xm[l] << xp[atom1] << xp[atom2]<< endl;
        double len = std::sqrt(bx*bx + by*by + bz*bz);
        // cout << "bond length = " << len << endl;
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




// for (int n = 1; n <= Nchain; ++n) { // n <= number of chains
//     for (int i = 1; i < Lchain; ++i) { //l < length of chain 
//         int atom1 = (n - 1) * Lchain + i;   // 1 .. Nmon-1
//         int atom2 = atom1 + 1;              // 2 .. Nmon
//         int l     = (n - 1) * (Lchain - 1) + i; // 1 .. Nbond

//         bondvecfile1 << atom1 << " " << atom2 << " " << bond[l][0] << " " << bond[l][1] << " " << bond[l][2] << endl;
//     }
// }

// Grid assignment 

double*** nn1=Allocate_3D_Double_Array(nx+1, ny+1, nz+1);
double*** nn2=Allocate_3D_Double_Array(nx+1, ny+1, nz+1);
double*** nn3=Allocate_3D_Double_Array(nx+1, ny+1, nz+1);

double Lgrid;
Lgrid=2.0;

nx=(int) Lx/Lgrid;
ny=(int) Ly/Lgrid;
nz=(int) Lz/Lgrid;
cout<< "nx=" << nx<< " ny= " << ny << " nz=  "<< nz << endl;
lx=Lx/(nx*1.0);
ly=Ly/(ny*1.0);
lz=Lz/(nz*1.0);   // we have a  nx*ny*nzy grid whose mesh sides are given by lx, ly, lz
double vol=lx*ly*lz;
//lx=ly=lz=Lgrid
 int Ngrid=nx*ny*nz;
cout<< "lx=" << lx<< " ly= " << ly << " lz=  "<< lz << "   Ngrids=" << nx*ny*nz << endl;
double u[1000][4];

ofstream bondvecfile1;
bondvecfile1.open(bondvec);
bondvecfile1 << "atom1 bx by bz xm ym zm nx ny nz"  << endl;


for(n1=0;n1<nx;n1++)
  for(n2=0;n2<ny;n2++)
    for(n3=0;n3<nz;n3++){
      m=0;

//cout << n1<< " " << n2 << " "<< n3 << endl;
      for(l=1; l<=10;l++){
    //cout << "l=" << l << endl;
        m1=(int) xm[l]/lx;
        m2=(int) ym[l]/ly;
        m3=(int) zm[l]/lz;
        cout << m1 << xm[l] << lx << endl;

        if( (m1==n1) && ( m2==n2) && (m3==n3) ){

          ++m;
          u[m][1]=bond[l][0];
          u[m][2]=bond[l][1];
          u[m][3]=bond[l][2];

          bondvecfile1 << atom1_of_bond[l] << " "
                       << bond[l][0] << " "
                       << bond[l][1] << " "
                       << bond[l][2] << " "
                       << xm[l] << " "
                       << ym[l] << " "
                       << zm[l] << " "
                       << n1 << " "
                       << n2 << " "
                       << n3 << std::endl;
        }

      }
}



}//end program



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

 void orderparameter1( int npar, double **um, double lambda[], double eevp[4],  double eevm[4], double eevo[4], double *S, double director[3] )
{//cout << "up to here is fine!!" << endl;




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
 if(a11> 1.0) cout << "s.th wrong";
 if(a12> 1.0) cout << "s.th wrong";
 if(a13> 1.0) cout << "s.th wrong";
if(a22> 1.0) cout << "s.th wrong";
if(a23> 1.0) cout << "s.th wrong";
if(a33> 1.0) cout << "s.th wrong";




    cc1=a12*a12+a13*a13+a23*a23-a11*a22-a11*a33-a22*a33;
      cc0=a11*a22*a33+2.0*a12*a13*a23-a12*a12*a33-a13*a13*a22-a23*a23*a11;
      pp=0.25*(1.00+3.0*cc1);
      qq=(27.0*cc0+9.0*cc1+2.0)/16.0;
   
      phia=acos(qq/sqrt(pow(pp,3)))/3.;
      // cout << "phia=" << phia << endl;
      rr=2.0 *sqrt(pp);

   //   cout << "pp=" << pp << "  qq= " << qq<<  " rr=" << rr << " phia=" << phia << endl; 
      lambda[1]=rr*cos(phia);
      lambda[2]=rr*cos(phia-2.*pi2/3.);
      lambda[3]=rr*cos(phia-pi2/3.);
//if(npar==1)  cout << "lambda[1]=" << lambda[1] <<  "lambda[2]=" << lambda[2] <<  "lambda[3]=" << lambda[3] << endl; 
      aa11=a11-onethird-2.*onethird*lambda[1];
      aa22=a22-onethird-2.*onethird*lambda[1];
      eevp[1]=a12*a23-a13*aa22;
      eevp[2]=a13*a12-a23*aa11;
      eevp[3]=aa11*aa22-a12*a12;
      vnorm=sqrt(eevp[1]*eevp[1]+eevp[2]*eevp[2]+eevp[3]*eevp[3]);
      if(vnorm > 1.0e-12) {
      eevp[1]=eevp[1]/vnorm;
      eevp[2]=eevp[2]/vnorm;
      eevp[3]=eevp[3]/vnorm;}
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
      eevm[3]=eevm[3]/vnorm;}
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


//ordering of order parameter eigenvalues
double S_nem, S1,S2,S3,S0=0;
//double director[3];

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
//nematicfile << lambda[1] << " " << lambda[2] << " "<< lambda[3] << endl;

nematicfile << "ss= " << ss << endl;
*S=S_nem;
      return;

     }
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 void orderparameter( int npar, double um[][4], double lambda[], double eevp[4],  double eevm[4], double eevo[4], double *S, double director[3] )
{//cout << "up to here is fine!!" << endl;

   //   onethird=1./3.
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
 if(a11> 1.0) cout << "s.th wrong";
 if(a12> 1.0) cout << "s.th wrong";
 if(a13> 1.0) cout << "s.th wrong";
if(a22> 1.0) cout << "s.th wrong";
if(a23> 1.0) cout << "s.th wrong";
if(a33> 1.0) cout << "s.th wrong";




    cc1=a12*a12+a13*a13+a23*a23-a11*a22-a11*a33-a22*a33;
      cc0=a11*a22*a33+2.0*a12*a13*a23-a12*a12*a33-a13*a13*a22-a23*a23*a11;
      pp=0.25*(1.00+3.0*cc1);
      qq=(27.0*cc0+9.0*cc1+2.0)/16.0;
   
      phia=acos(qq/sqrt(pow(pp,3)))/3.;
      // cout << "phia=" << phia << endl;
      rr=2.0 *sqrt(pp);

   //   cout << "pp=" << pp << "  qq= " << qq<<  " rr=" << rr << " phia=" << phia << endl; 
      lambda[1]=rr*cos(phia);
      lambda[2]=rr*cos(phia-2.*pi2/3.);
      lambda[3]=rr*cos(phia-pi2/3.);
//if(npar==1)  cout << "lambda[1]=" << lambda[1] <<  "lambda[2]=" << lambda[2] <<  "lambda[3]=" << lambda[3] << endl; 
      aa11=a11-onethird-2.*onethird*lambda[1];
      aa22=a22-onethird-2.*onethird*lambda[1];
      eevp[1]=a12*a23-a13*aa22;
      eevp[2]=a13*a12-a23*aa11;
      eevp[3]=aa11*aa22-a12*a12;
      vnorm=sqrt(eevp[1]*eevp[1]+eevp[2]*eevp[2]+eevp[3]*eevp[3]);
      if(vnorm > 1.0e-12) {
      eevp[1]=eevp[1]/vnorm;
      eevp[2]=eevp[2]/vnorm;
      eevp[3]=eevp[3]/vnorm;}
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
      eevm[3]=eevm[3]/vnorm;}
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


//ordering of order parameter eigenvalues
double S_nem, S1,S2,S3,S0=0;
//double director[3];

S1=abs(lambda[1]);
S2=abs(lambda[2]);
S3=abs(lambda[3]);
S0=S1; director[0]=eevp[1]; director[1]=eevp[2]; director[2]=eevp[3]; S_nem=lambda[1];
if(S2 > S0) {S0=S2; director[0]=eevm[1]; director[1]=eevm[2]; director[2]=eevm[3]; S_nem=lambda[2];}
if(S3 > S0)  {S0=S3; director[0]=eevo[1]; director[1]=eevo[2]; director[2]=eevo[3]; S_nem=lambda[3];}
double ss=sqrt(2.0/3.0*(lambda[1]*lambda[1]+lambda[2]*lambda[2]+lambda[3]*lambda[3])); 

nematicfile << S_nem << " " <<ss << " " << lambda[1]  << " "<< lambda[2] << " "<< lambda[3] << " " << director[0]  << " "<< director[1] << " "<< director[2] <<  endl;
// nematicfile << S0 << " " << director[0] << " "<< director[1] << " "<< director[2] << endl;
nematicfile << eevp[1]<< " " << eevp[2] << " "<< eevp[3] << endl;
nematicfile << eevm[1]<< " " << eevm[2] << " "<< eevm[3] << endl;
nematicfile << eevo[1]<< " " << eevo[2] << " "<< eevo[3] << endl;
nematicfile << lambda[1] << " " << lambda[2] << " "<< lambda[3] << endl;

nematicfile << "ss= " << ss << endl;
*S=S_nem;
      return;

     }


 
