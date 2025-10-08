//g++ -I ./ -o nematic20-new10.exe nematic20-new.cc -O2 -Wno-deprecated
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
int const Lchain=10; // number of monomers in a chain
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



cout << "NEw long array created successfully" << endl;


cout << "argv[1] = "<< argv[1] << "   its length is" <<  strlen(argv[1]) << endl;// anouncing the name of input file for initial values

const int fnamelength=strlen(argv[1]);

char filename[fnamelength];  
for(i=0; i<= fnamelength; i++)

		filename[i] = argv[1][i];  //converting the input argument to the file name string



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
   cout << Nmon << str1<< endl;
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
   cout << "zlo="<< zlo<<  " zhi=" << zhi << endl;
Lx=xhi-xlo;
Ly=yhi-ylo;
Lz=zhi-zlo;  
double   Lxhalf=Lx/2.;
   double   Lyhalf=Ly/2.;
   double   Lzhalf=Lz/2.;


 xyzfile >>  str1   ;

for(i=1;i<=Natomtype;i++)
 xyzfile >>  j >> l   ;

 xyzfile >>  str1 >> str2 >> str3 >> str4  ;
cout << str4 << endl;
for(i=1;i<=Natomtype;i++)
xyzfile >> xx >> yy >> zz ;


xyzfile >>  str1  >> str2 >> str3 >> str4 ;

xyzfile >> xx >> yy >> zz ;

xyzfile >> str1 >> str2 >> str3   ;

cout << str3 << endl;

      for(i=1;i<=Nmon; i++)
     { xyzfile >>  id >> molid>> type >> xx >> yy >> zz >> ix >> iy >> iz  ;
        xu[id]=xx+ix*Lx;
     yu[id]=yy+iy*Ly;
     zu[id]=zz+iz*Lz;
    xp[id]=xx-xlo; //shifting coordinates so that they are all positive starting from 0
    yp[id]=yy-ylo;
    zp[id]=zz-zlo;
   
}

// xyzfile >> str1;
//cout <<  str1<< endl;
// xyzfile >> str1;
//cout <<  str1<< endl;
 xyzfile >> str1;
cout <<  str1<<  endl;
      
     for(i=1;i<=Nmon; i++)
     { xyzfile >>  id >>  xx >> yy >> zz   ;
    if(i==10000) cout << id <<  xx << " "<< yy << " " << zz << endl;  } 
    xyzfile >> str1;
cout <<  str1<<  endl;

 for(i=1;i<=Nbond; i++)
     { xyzfile >>  j>> l >> id1 >> id2   ;
      //bond1[i][0]=id1;
    // bond1[i][1]=id2; 
    if(i==9900) cout << j << " "<<  l << " "<< id1 << " " << id2 << endl;  } 
    


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

for(n=0;n<1000;n++){
count[n]=0;
nematiccount[n]=0;}


for(n=1; n <= Nchain;n++)
 for(i=1;i< Lchain;i++){
l=(n-1)*(Lchain-1)+i;
//cout << "l="<< l<< endl;

bond[l][0]=xu[l+1+(n-1)]-xu[l+(n-1)];
bond[l][1]=yu[l+1+(n-1)]-yu[l+(n-1)];
bond[l][2]=zu[l+1+(n-1)]-zu[l+(n-1)];
xm[l]=(xp[l+1+(n-1)]+xp[l+(n-1)])/2;
if(xm[l]> Lx)  cout<< "l=" << l << "xm=  "<< xm[l] << endl;
ym[l]=(yp[l+1+(n-1)]+yp[l+(n-1)])/2;
//yum[l]=(yu[l+1+(n-1)]+yu[l+(n-1)])/2;
if(ym[l]> Ly)  cout<< "l=" << l << "ym=  "<< ym[l] << endl;
zm[l]=(zp[l+1+(n-1)]+zp[l+(n-1)])/2;
//zum[l]=(zu[l+1+(n-1)]+zu[l+(n-1)])/2;
if(zm[l]> Lz)  cout<< "l=" << l << "zm=  "<< zm[l] << endl;
lbond[l]=sqrt(bond[l][0]*bond[l][0]+bond[l][1]*bond[l][1]+bond[l][2]*bond[l][2]);
if(lbond[l]> 1) cout << "l=" << l << " lbond[l]=" << lbond[l] << endl;
bond[l][0]=bond[l][0]/lbond[l];
bond[l][1]=bond[l][1]/lbond[l];
bond[l][2]=bond[l][2]/lbond[l];
//if(l==100) cout << "bond10" << bond[10][0] << " " << bond[10][1]<<  " " << bond[10][2] << endl;
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

char outname1[50];
strcpy(outname1,filename);
strcat(outname1, "_nem.tcl");

ofstream vmdfile;
vmdfile.open(outname1);
vmdfile << "mol new" << endl;
vmdfile << "color Display Background silver" << endl;

strcpy(outname1,filename);
strcat(outname1, "_max.tcl");
ofstream vmdfile1;
vmdfile1.open(outname1);
vmdfile1 << "mol new" << endl;
vmdfile1 << "color Display Background silver" << endl;

////////// calculation of R_g ////////////////////////////////////
 
int idm,idnext,idpast;
ofstream  MSIDfile;
char MSIDname[50];
strcpy(MSIDname,filename);
strcat(MSIDname, "_MSID.dat");

MSIDfile.open(MSIDname);
double xcm[Nchain], ycm[Nchain], zcm[Nchain], Rg1[Nchain], Rg2[Nchain], Rg3[Nchain], Rg=0, Rgx=0, Rgy=0, Rgz=0;

for(n=0;n< Nchain;n++){
 xcm[n]=0;
ycm[n]=0;
zcm[n]=0;
Rg1[n]=0;
Rg2[n]=0;
Rg3[n]=0;
}


for(i=0; i< Nchain;i++){
  for(j=1;j<=Lchain;j++){ //calculation of R_cm
     id1=i*Lchain+j; // id of jth atom in the chain i
  //  }} cout << "id1=" << id1 << endl;

 xcm[i]+=xu[id1]/Lchain;
     ycm[i]+=yu[id1]/Lchain;
     zcm[i]+=zu[id1]/Lchain;
}
  
  for(j=1;j<=Lchain;j++){
     id1=i*Lchain+j; // id of jth atom in the chain i
     Rg1[i]+=(xu[id1]-xcm[i])*(xu[id1]-xcm[i]);
     Rg2[i]+=(yu[id1]-ycm[i])*(yu[id1]-ycm[i]);
     Rg3[i]+=(zu[id1]-zcm[i])*(zu[id1]-zcm[i]);
}

Rgx+= Rg1[i]/Lchain;
Rgy+= Rg2[i]/Lchain;
Rgz+= Rg3[i]/Lchain;
//Rg+=sqrt((Rg1[i]+ Rg2[i]+Rg3[i])/Lchain);
Rg+=(Rg1[i]+ Rg2[i]+Rg3[i])/Lchain;
}// end for i

Rg=Rg/Nchain;
Rgx=Rgx/Nchain;
Rgy=Rgy/Nchain;
Rgz=Rgz/Nchain;


Rg=sqrt(Rg);
Rgx=sqrt(Rgx);
Rgy=sqrt(Rgy);
Rgz=sqrt(Rgz);


MSIDfile << "Rg, Rgx, Rgy, Rgz "<< Rg << " "<<  Rgx << " "<<  Rgy << " "<<  Rgz << endl;

////////// calculation of persistence length ////////////////////////////////////


double lp,lpave=0,R1;


for(n=0;n< Nchain;n++){
lp=0;
for(j=1;j<=Lchain;j++){

id1=n*Lchain+1;
id2=n*Lchain+j;
lp+=xu[id1]*xu[id2]+yu[id1]*yu[id2]+zu[id1]*zu[id2];

}
R1=xu[id1]*xu[id1]+yu[id1]*yu[id1]+zu[id1]*zu[id1];
R1=sqrt(R1);
lp=lp/R1;
lpave+=lp;
}
lp=lp/Nchain;
MSIDfile << "persistence length in the limit of N=" << Lchain << " is " << lp << endl;

/////////////////////////////////////////////////calculation of MSID ///////////////////////////////////////////////////////////



double Mr2[Lchain],rr2, MSID[Lchain],  MSIDx[Lchain], MSIDy[Lchain], MSIDz[Lchain];
for(n=1;n< Lchain;n++){
 Mr2[n]=0;
  MSIDx[n]=0;  MSIDy[n]=0; MSIDz[n]=0;
}
double dx,dy,dz,r2max;





  for(i=0; i<Nchain;i++)
  for(n=1;n<Lchain;n++)
  for(j=1; j<= (Lchain-n);j++){
   id1=i*Lchain+j; // id of jth atom in the chain i
   id2=id1+n; // id of (j+n)th atom in the chain i
   cnd[1]=xu[id1]-xu[id2];
   cnd[2]=yu[id1]-yu[id2];
   cnd[3]=zu[id1]-zu[id2];
   rr2=cnd[1]*cnd[1]+cnd[2]*cnd[2]+cnd[3]*cnd[3];
    Mr2[n]+=rr2/(Lchain-n);
    MSIDx[n]+=cnd[1]*cnd[1]/(Lchain-n);
    MSIDy[n]+=cnd[2]*cnd[2]/(Lchain-n);
    MSIDz[n]+=cnd[3]*cnd[3]/(Lchain-n);
}
 
MSIDfile << " the MSID is calculated for " << Nchain << " chains of length "<< Lchain << endl;
MSIDfile << "n" << "  "<< "MSID(n)"<< "  "<< "MSID(n)/n"<<"  "<< "MSIDx(n)"<< "  "<< "MSIDy(n)"<< "  "<< "MSIDz(n)" << endl;

for(n=1;n< Lchain;n++){
 MSID[n]=Mr2[n]/Nchain;
MSIDx[n]=MSIDx[n]/Nchain;
MSIDy[n]=MSIDy[n]/Nchain;
MSIDz[n]=MSIDz[n]/Nchain;

 MSIDfile << n << " "<< MSID[n]<< " "<< MSID[n]/n << " "<< MSIDx[n] << " "<< MSIDy[n] << " "<< MSIDz[n] <<  endl;}

/////////////////////////////////////////////////calculation of  End-to-End vector Nematic tensor ///////////////////////////////////////////////////////////
double**  Ree=Allocate_2D_Double_Array(Nchain+1,3);
double *Re;
Re=new double [Nchain+1];
 
 for(i=0; i<Nchain;i++){
   id1=i*Lchain+1; // id of firstth atom in the chain i
   id2=id1+Lchain-1; // id of Lchain_th atom in the chain i
   Ree[i+1][0]=xu[id1]-xu[id2];
   Ree[i+1][1]=yu[id1]-yu[id2];
   Ree[i+1][2]=zu[id1]-zu[id2];
   Re[i+1]=sqrt(Ree[i+1][0]*Ree[i+1][0]+Ree[i+1][1]*Ree[i+1][1]+Ree[i+1][2]*Ree[i+1][2]);
   Ree[i+1][0]=Ree[i+1][0]/Re[i+1];
   Ree[i+1][1]=Ree[i+1][1]/Re[i+1];
   Ree[i+1][2]=Ree[i+1][2]/Re[i+1];
}
char pdfseg[50];

strcpy(pdfseg,filename);
strcat(pdfseg, "_segsummary.txt");
ofstream segfile1;
segfile1.open(pdfseg);

orderparameter1(Nchain, Ree, lambda, eevp, eevm, eevo, &S, director);
double S_Ree=S;
 segfile1 << "S_Ree lambda_1_Ree  lambda_2_Ree  lambda_3_Ree director_1_Ree  director_2_Ree  director_3_Ree  "  << endl;
segfile1 << S_Ree << " " << lambda[1]  << " "<< lambda[2] << " "<< lambda[3] << " " << director[0]  << " "<< director[1] << " "<< director[2] << endl;
//////////////////////////////////////////////////////////////////// calculation of  bond-bond correlation ///////////////////////////////////////////////

strcpy(outname1,filename);
strcat(outname1, "_bondcorr.dat");

ofstream bondcorrfile;
bondcorrfile.open(outname1);
bondcorrfile << "i   P1       P2  "<< endl;

strcpy(outname1,filename);
strcat(outname1, "_intercorr.dat");

ofstream bondcorrfile1;
bondcorrfile1.open(outname1);
bondcorrfile1 << "r   P2       P2-intra  "<< endl;



int l1,l2,im, ir;

double costheta,P1Qi[Lchain], P2Qi[Lchain],P2cos;
P2Qi[0]=Nchain;  P1Qi[0]=Nchain; 


for(i=1;i< Lchain;i++) {
P2Qi[i]=0;
P1Qi[i]=0;
}

//////////////////////////////////////////////////////////////////// calculation of intra bond-bond correlation ///////////////////////////////////////////////

int ntheta, nbin=180;
double Ibintheta=nbin/pi,sumtest=0;
int histQ[nbin+1];
double theta, pdfQ[nbin+1]; // pdf of angle between two consequet bonds
for(i=0;i<=nbin;i++)  histQ[i]=0;

for(n=0; n< Nchain;n++)
  for(i=1;i < (Lchain-1);i++)
  for(j=1; j< (Lchain-i);j++){
   l1=n*(Lchain-1)+j; //id of bond connecting atom j,j+1 of chain n
   l2=n*(Lchain-1)+i+j; // bond connecting atom j+i,j+i+1 of chain n
costheta=bond[l1][0]*bond[l2][0]+bond[l1][1]*bond[l2][1]+bond[l1][2]*bond[l2][2];
P1Qi[i]+=costheta/(Lchain-i-1);
P2cos=0.5*(3.0*costheta*costheta-1.0);
P2Qi[i]+=P2cos/(Lchain-i-1);
  if(i==1){
         theta=acos(costheta);
          ir=theta*Ibintheta;
         //cout << "costheta= " << costheta <<  "  theta=" << theta << " ir=" << ir << endl;
         ++histQ[ir];
  }
}


 for(i=0; i < Lchain-1;i++){
 
 bondcorrfile <<  i << "   "  <<  P1Qi[i]/Nchain << "   "  <<  P2Qi[i]/Nchain << endl;}

strcpy(outname1,filename);
strcat(outname1, "_pdfQ.dat");

ofstream pdfQfile;
pdfQfile.open(outname1);
pdfQfile << "i   theta  nhistQ   nhist/Nbond  pdfQ  theta_degree   pdfQ_degree  "<< endl;


for(i=0; i <= nbin; i++){
 pdfQ[i]=histQ[i]*1.0/Nangle ;
sumtest+=pdfQ[i];
pdfQfile << i << " "<< i/Ibintheta <<  " "<< histQ[i]<<  " " << pdfQ[i] << " "<<  pdfQ[i]*Ibintheta <<" "<< i/Ibintheta/pi*180.0 << " "<<  pdfQ[i]*Ibintheta*pi/180.0 <<  endl;

}
cout <<  " sumtest=" << sumtest << endl;



//////////////////////////////////////////////////////////////////// calculation of interchain bond-bond correlation excluding intrachains ///////////////////////////////////////////////

double IDeltar;

r2max=min(Lx,Ly);
r2max=min(r2max,Lz)+2;
r2max=0.25*r2max*r2max;

double Deltar=1.0/0.2;

ir=r2max*Deltar;
cout << "rmax=" <<  sqrt(r2max) << " ir=" << ir << " 1/Deltar=" << 1.0/Deltar << endl;
long int Nbondcorr[ir+1],  Nbondcorr2[ir+1];
double  interCorr1[ir+1], interCorr2[ir+1],rx,ry,rz;
 interCorr1[0]=1;  interCorr2[0]=1;
/*
for(i=1;i<= ir;i++) {
interCorr1[i]=0;
interCorr2[i]=0;
Nbondcorr[i]=0;
Nbondcorr2[i]=0;
}

for(n1=1; n1< Nchain; n1++)
for(n2=n1; n2<= Nchain; n2++)
for(i=1;i < Lchain;i++)
  for(j=1;j < Lchain;j++){
 l1=(n1-1)*(Lchain-1)+i;
  l2=(n2-1)*(Lchain-1)+j;
//if(n1==1 && n2==Nchain) cout << "l1=" << l1 << " l2=" << l2 << endl;

rx=xm[l1]-xm[l2];
if (abs(rx) >   Lxhalf) {if(rx >0) rx = rx-Lx; 
                          if(rx <0) rx = rx+Lx;}
ry=ym[l1]-ym[l2];
if (abs(ry) >   Lyhalf) {if(ry >0) ry = ry-Ly; 
                          if(ry <0) ry = ry+Ly;}
rz=zm[l1]-zm[l2];
if (abs(rz) >   Lzhalf) {if(rz >0) rz = rz-Lz; 
                          if(rz <0) rz = rz+Lz;}

rr2=rx*rx+ry*ry+rz*rz; //check the way rr2 should be calculated
if(rr2 > 0.1 && rr2<r2max){
costheta=bond[l1][0]*bond[l2][0]+bond[l1][1]*bond[l2][1]+bond[l1][2]*bond[l2][2];
P2cos=0.5*(3.0*costheta*costheta-1.0);
im=(int)trunc(rr2*Deltar);
//if(n1==1 && n2==2) cout << "rr2=" << rr2 << "im=" << im << endl;
 Nbondcorr[im]++;
 interCorr1[im]+=P2cos;
 if(n2 > n1){ 
Nbondcorr2[im]++;
interCorr2[im]+=P2cos;}
}}

 IDeltar=1.0/Deltar;

    for(i=0; i<ir;i++){
       //  myfile << "ig=" << ig << "g[ig]=" << g[ig] << endl;               
      // r[i]=sqrt(i*IDeltar);
       interCorr1[i]= interCorr1[i]/Nbondcorr[i];
       interCorr2[i]= interCorr2[i]/Nbondcorr2[i];
                }

for(i=1; i < ir ;i++){
 
 bondcorrfile1 << sqrt(i*IDeltar) << "   "  <<  interCorr1[i] << "   "  <<  interCorr2[i] << endl;}

*/





//////////////////////////////////////////making the grid//////////////////////////////////////////////////////////////////////////////////////////
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
//int NN=nx*nx*ny*ny*nz*nz/2;
//each mesh of the grid is determined by its index (n1,n2,n3)

int nhist=800;
 int ip,ipx,ipy,ipz;
double Deltas= nhist/Dhalf2;
double corr[nhist+1],fQ[nhist+1], r[nhist+1], corrx[nhist+1],fQx[nhist+1], corry[nhist+1],fQy[nhist+1],corrz[nhist+1],fQz[nhist+1];
unsigned long g[nhist+1], gx[nhist+1], gy[nhist+1], gz[nhist+1];
for(ip=0;ip<=nhist;ip++){
corr[ip]=0;
corrx[ip]=0;
corry[ip]=0;
corrz[ip]=0;
g[ip]=0;
gx[ip]=0;
gy[ip]=0;
gz[ip]=0;
}


//double SS[nx+1][ny+1][nz+1],  nn1[nx+1][ny+1][nz+1],  nn2[nx+1][ny+1][nz+1],  nn3[nx+1][ny+1][nz+1]; int mm[nx+1][ny+1][nz+1]; // order parameter value & orientation  vector of grid elements

double*** SS=Allocate_3D_Double_Array(nx+1, ny+1, nz+1);
double*** nn1=Allocate_3D_Double_Array(nx+1, ny+1, nz+1);
double*** nn2=Allocate_3D_Double_Array(nx+1, ny+1, nz+1);
double*** nn3=Allocate_3D_Double_Array(nx+1, ny+1, nz+1);
int*** mm=Allocate_3D_Integer_Array(nx+1, ny+1, nz+1);

int ncolor[nx+1][ny+1][nz+1]; 
int idcell[10000], *amorphLabel;

amorphLabel=new int [Nmon+1];

cout<< " up to here fine" << endl;
double r1,r2,r3,rr,P2Q,cosQ, P2Q1,cosQ1 , P2Q2,cosQ2;



char nemname[50];
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
 int Namorph=0, tt,nS,Nss; //tt the numbr of bins for the cos angle of local nematic director.
 tt= (int) (pi/2.0/bintheta)+1;  cout <<"tt=" << tt << endl;
 
 double nem=0, thetaz=0, thetay=0, thetax=0;
    
   int Ns0= (int) 1.0/binS;
    Nss= Ns0*3/2+1; //  -0.5=< S<= 1
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
 
for(n1=0;n1<nx;n1++)
for(n2=0;n2<ny;n2++)
for(n3=0;n3<nz;n3++){
m=0;

//cout << n1<< " " << n2 << " "<< n3 << endl;
for(l=1; l<=Nbond;l++){
 //cout << "l=" << l << endl;
m1=(int) xm[l]/lx;
m2=(int) ym[l]/ly;
m3=(int) zm[l]/lz;

//cout<< "m1=" << m1<< " m2= " << m2 << " m3=  "<< m3 << endl;
 if( (m1==n1) && ( m2==n2) && (m3==n3) ){

++m;
u[m][1]=bond[l][0];
u[m][2]=bond[l][1];
u[m][3]=bond[l][2];
 idcell[m]=l/(Lchain-1)+l; //id of atom at the begining of bond vector.
 idbond[m]=l;
/*xp1[m]=xm[l]-lbond[l]/2*bond[l][0];
xp2[m]=xm[l]+lbond[l]/2*bond[l][0];

yp1[m]=ym[l]-lbond[l]/2*bond[l][1];
yp2[m]=ym[l]+lbond[l]/2*bond[l][1];

zp1[m]=zm[l]-lbond[l]/2*bond[l][2];
zp2[m]=zm[l]+lbond[l]/2*bond[l][2];*/
} //end if

}// end for l
if(m>max) max=m;
++count[m];
if(m > 1) {
//cout << "m=" << m << endl;
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
    else {//cout << "nS=" << nS << endl;
        nS=abs(nS)+Ns0;
        countS[nS]= countS[nS]+1;
   //cout << "nS_new=" << nS << "countS[nS]=" << countS[nS] << endl;
    }
}
if(SS[n1][n2][n3] > 0.8){ 
ncolor[n1][n2][n3]=1; // initial label of all the nematic grid elements
++nematiccount[m];
++gridcount;

 thetaz=acos(abs(director[2]));
	  nS= (int) (thetaz/bintheta);
	 // cout << "eevp[3]=" << eevp[3] << " thetaz=" <<  thetaz << "  nS=" << nS <<  endl;
     ++counttheta3[nS];
     thetay=acos(abs(director[1]));
	  nS= (int) (thetay/bintheta);
	      ++counttheta2[nS];
       thetax=acos(abs(director[0]));
	  nS= (int) (thetax/bintheta);
	      ++counttheta1[nS];
 } //end if S>SC

else { sumamorph+=m;
ncolor[n1][n2][n3]=0;
for(i=1;i<=m;i++){
//cout << "idcell[i]=" << idcell[i] << endl;
 bondamorph[idbond[i]]=1;
amorphLabel[idcell[i]]=1;
//Namorph++;
amorphLabel[idcell[i]+1]=1;
}
}//end else

} //end for n3
 cout <<  " sumbondamorph="<< sumamorph << endl;
cout << "nematic gridcount=" << gridcount << endl;
ndomains=gridcount;

 Namorph=0;
for(l=1;l<=Nbond;l++){
  if( bondamorph[l]==1)  ++Namorph;
 else {bondamorph[l]=-1;}
 }


 cout << "Namorph=" << Namorph << " fraction of amorph bonds=" << (Namorph*1.0)/Nbond <<  endl;

delete idbond;
 //////////////////////////////////////////////////histogram of number of bonds per cell/////////////////////////////////////////////////

    unsigned long  sumcount=0, sumcount1=0, sumcount2=0, sumcount3=0;
    double rho_nem=0, rho_amorph=0;
//double SS;

char outname[50];
strcpy(outname,filename);
strcat(outname, "_hist.txt");

ofstream histfile;
histfile.open(outname);
histfile <<  "m  " << "count[m]  count[m]*1.0/Nbond  nematic_count[m]  fraction_ nematic_cell[m] " << endl;
for(m=0; m<=max;m++){
if(count[m]>0){
histfile << m << " " << count[m] << "  "<< count[m]*1.0/Nbond << "  " << nematiccount[m] << "  " << nematiccount[m]/ (1.0*count[m])  << endl;
sumcount+=count[m]*m;
sumcount1+=count[m];
sumcount2+=nematiccount[m];
sumcount3+=nematiccount[m]*m;
rho_nem+=nematiccount[m]*m/vol;
rho_amorph+=(count[m]-nematiccount[m])*m/vol;
} //end if
}
    
    rho_amorph=rho_amorph/(Ngrid-sumcount2);
    rho_nem=rho_nem/sumcount2;
    
histfile << " total number of grid elements is:" << sumcount1 << endl;
histfile << " total number of bonds from sumcount is:" << sumcount << endl;
histfile << " total number of nematic grid elements is:" << sumcount2 << endl;
histfile << " total number of bonds in the nematic grid elements is:" << sumcount3<< endl;
histfile << " degree of crystallinity is:" << sumcount2*1.0/ sumcount1 << endl;

    
    
double crystallinity;
crystallinity=sumcount2*1.0/ sumcount1;


//ndomains=sumcount2-ndomains/3;
//histfile << " number of independent crystalline domains is:" << ndomains << endl;

 


  ///////////////////////////////////////////////////////////Histogram of local nematic order parameter and its director/////////////////////////////////////////////////////////////
 strcpy(outname,filename);
strcat(outname, "_Nemhist.txt");

ofstream Nemhistfile;
Nemhistfile.open(outname);
    Nemhistfile << "rho rho_amorph rho_nem v_mon v_amorph v_nem " << rho  << " " << rho_amorph << " " << rho_nem << " " << 1.0/rho <<  " " << 1.0/rho_amorph << " " << 1.0/rho_nem << endl;


int Ntest=0;

Nemhistfile <<  "i  " << "S "<< "count[S]" <<   " pdf[S]" << endl;
     for(i=Nss-1;i>=Ns0;i--){ Nemhistfile << i << "  " << (Ns0-i-0.5)*binS <<  " "<< countS[i] << " " <<   countS[i]*1.0/Ngrid/binS << endl; 
//Ntest+=countS[i];
}
for(i=0;i<Ns0;i++){ Nemhistfile << i << "  " << (i+0.5)*binS <<  " "<< countS[i] << " " <<   countS[i]*1.0/Ngrid/binS << endl;
//Ntest+=countS[i];
 }

//cout << "number of nematic cells is " << gridcount << "sumcountS=" << Ntest << endl;

 strcpy(outname,filename);
strcat(outname, "_directorhist.txt");

ofstream directorfile;
directorfile.open(outname);
directorfile <<  "theta "<< "count1[theta]" <<   "pdf1[theta]" << "count2[theta]" <<   "pdf2[theta]"<< "count3[theta]" <<   "pdf3[theta]"<< endl;
 for(i=0;i<= tt;i++){ directorfile << i*bintheta << "  " <<  counttheta1[i] << " " <<   counttheta1[i]*1.0/gridcount/bintheta <<  " " << counttheta2[i] << " " <<   counttheta2[i]*1.0/gridcount/bintheta <<  " " << counttheta3[i] << " " <<   counttheta3[i]*1.0/gridcount/bintheta << endl; }

 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

 
double cosphi;

for(n1=0;n1<nx;n1++)
for(n2=0;n2<ny;n2++)
for(n3=0;n3<nz;n3++){
for(m1=n1;m1<nx;m1++)
for(m2=n2;m2<ny;m2++)
for(m3=n3;m3<nz;m3++){

if(SS[n1][n2][n3] >0.8 && SS[m1][m2][m3] >0.8){
r1=(m1-n1)*lx;
r2=(m2-n2)*ly;
r3=(m3-n3)*lz;
if (r1 > Lxhalf)  r1 =Lx-r1;
if (r2 > Lyhalf)  r2 =Ly-r2;
if (r3 >Lzhalf)  r3 =Lz-r3;
rr=r1*r1+r2*r2+r3*r3;
cosphi=nn1[n1][n2][n3]*r1+nn2[n1][n2][n3]*r2+nn3[n1][n2][n3]*r3;
  cosQ=nn1[n1][n2][n3]*nn1[m1][m2][m3]+nn2[n1][n2][n3]*nn2[m1][m2][m3]+nn3[n1][n2][n3]*nn3[m1][m2][m3];
    P2Q=(3.0* cosQ*cosQ-1.0);
//**********  calculation of orientational correlations within the nematic domains ***********************************//

   ip=(int)trunc(rr*Deltas);
     g[ip]+=2; 
     corr[ip]+=P2Q;

if(0.9 < cosphi &&  cosphi < 1.0){
       gx[ip]+=2; 
     corrx[ip]+=P2Q;
}


if(0 < cosphi &&  cosphi < 0.1){
        gy[ip]+=2; 
     corry[ip]+=P2Q;
}


                          



} // if
}//end for2
} //end for 1

 



double IDelta2=1.0/Deltas;

    for(i=0; i<nhist;i++){
       //  myfile << "ig=" << ig << "g[ig]=" << g[ig] << endl;               
       r[i]=sqrt(i*IDelta2);
       fQ[i]= corr[i]/g[i];
    fQx[i]= corrx[i]/gx[i];
    fQy[i]= corry[i]/gy[i];
        
                }




ofstream P2file;
char P2name[50];
strcpy(P2name,filename);
strcat(P2name, "_P2(R).txt");
P2file.open(P2name); 





P2file  << "r   g0(r)    P2(r)    P2_par(r)     P2_per(r)   " << endl;


 for(i=0; i<nhist;i++){
 if( g[i] >0 )  
P2file <<   r[i] << "   "  <<  g[i] << "   "  <<  fQ[i]  << "   "  <<  fQx[i] << "   "  <<  fQy[i] << endl;}

release_3D_Double_Array(SS,nx+1, ny+1, nz+1);

 //////////***************** statistical analysis of tie and crystalline segments*********//////////////////////////////////////

 

 
 
 int id0,seg_sign, sum_seg, Nsegmax=500;
 //int seg[Nchain][Nsegmax];// seg[i][j] tells about the number of monomers in jth segment of chain i, its sign tells us if they are amorphous or ordered part (1 amorphous, -1 ordered)

int **seg; //seg[Nchain][Nsegmax];
seg= new int *[Nchain];
for(int i = 0; i < Nchain; ++i){
seg[i] = new int [Nsegmax];
if (seg[i] == NULL) {cerr << " seg[i] Allocation problem!"; exit(1);}}





// int Nseg[Nchain+1]; //tells the numbers of segments in chain i
int *Nseg;

Nseg= new int [Nchain+1];

for(i=0; i<Nchain;i++)
   for(n=0; n < Nsegmax;n++)
    seg[i][n]=0;
  
 
     int segmax_amorph=1, segmax_crys=1;

 for(i=0; i<Nchain;i++){
   id=i*(Lchain-1)+1;
   Nseg[i+1]=1;
   j= Nseg[i+1];
   seg[i][j]=bondamorph[id];
   seg_sign=bondamorph[id];
  for(n=1; n < Lchain-1;n++){
    //   id=i*Lchain+n; // id of nth atom in the chain i
    idnext=id+1;
   
  if(bondamorph[idnext]==bondamorph[id]){
       seg[i][j]=(abs(seg[i][j])+1)*seg_sign;  }
  else {
    Nseg[i+1]=Nseg[i+1]+1;
    // cout << amorphLabel[id] << " "<< amorphLabel[idnext] << "Nseg= "<<  Nseg[i+1]<<  endl;
    seg_sign=bondamorph[idnext];
      j= Nseg[i+1];
      seg[i][j]=seg_sign;
      
  }// end else

  id=idnext;
    
  } //end for 2
  
  }// end for 1

 Nseg[0]=Lchain/10;

/* for(i=0;i<=Nchain;i++)
   cout <<  Nseg[i] << " " ;
   cout << endl;*/




int max_Nseg=1;

for(i=1;i<=Nchain;i++){
   if(Nseg[i] > max_Nseg) max_Nseg=Nseg[i];
}

cout << " Nseg_max is "  << max_Nseg << endl;
 
cout << "The largest number of chain segments is "  << *std::max_element(Nseg,Nseg+Nchain+1) << '\n';

 cout << "The smallest number of chain segments is "  << *std::min_element(Nseg,Nseg+Nchain+1) << '\n';





 //test of seg[i][j]
 int sumseg_tot=0;
 
 for(i=0; i<Nchain;i++)
   for(j=1; j<=Nseg[i+1];j++)
     sumseg_tot+=abs(seg[i][j]);


 cout<< "  sumseg_tot=" <<  sumseg_tot << endl;
   

 
for(i=0; i<Nchain;i++){
  for(j=1;j<=Nseg[i+1];j++){
    if(seg[i][j]>0 ){
      if( seg[i][j]> segmax_amorph)
      segmax_amorph=seg[i][j];
      }
  else {if(abs(seg[i][j])> segmax_crys)       segmax_crys= abs(seg[i][j]);
    
  } //end else
  } // end for 2
 } // end for 1
 

 cout << "max length of an amorphous segment="<< segmax_amorph << " "<< "max length of crystalline segment ="<< segmax_crys << endl;
 int segmax=Max(segmax_crys, segmax_amorph);

////////////////////////////// Identification of tie chains /////////////////////////////
 
 long int Ntie=0, Nseg_tie=0; //Nseg_tie tells  the total number of segments that act as tie chains.
long int  Count_tie[segmax+1];

for(i=0;i<=segmax;i++)  Count_tie[i]=0;

 for(i=0;i<Nchain;i++)
   for(n=2;n< Nseg[i+1];n++){
     if(seg[i][n]>2 && abs(seg[i][n-1])>2 && abs(seg[i][n+1])>2 ){
       Ntie+=seg[i][n];
       j=seg[i][n]; 
   ++Count_tie[j];
   ++Nseg_tie;
   
     }
   }// end for n
	
 double ftie;  




 ftie=(Ntie*1.0)/Nbond;
 cout << "fraction of tie bonds  is" << ftie << " Nseg_tie=" << Nseg_tie << endl;
 

 double  pdf_amorph[segmax+1],  pdf_crys[segmax+1], aveseg_amorph=0, aveseg_crys=0, aveseg_tie=0, Lseg_amorph=0, Lseg_tie=0, Lseg_crys=0;
 long int Nseg_crys=0, Nseg_amorph=0;
 int Na, Nc;
long int Ntot_crys=0;


 for(i=0;i<=segmax; i++){
   pdf_amorph[i]=0;
  // pdf_tie[i]=0;
}

 for(i=0;i<=segmax;i++) pdf_crys[i]=0;


 
 for(i=0; i<Nchain;i++){
  for(j=1;j<= Nseg[i+1];j++){
    if(seg[i][j]>0){      
       Na=seg[i][j];
     ++pdf_amorph[Na];
     ++Nseg_amorph; 
      
      }
    else { Nc=abs(seg[i][j]);
	++pdf_crys[Nc]; 
       ++Ntot_crys;   
    if(Nc>2) { ++Nseg_crys; }// number of crystalline segments  
  } //end else
  } // end for 2
 } // end for 1
 
  
 for(m=1; m<=segmax;m++){
    aveseg_amorph+=m*m*pdf_amorph[m];
      Lseg_amorph+=m*pdf_amorph[m];
    
   }
    
  Nseg_tie=0;
// long int Nmoncrys=0;

 for(m=3; m<=segmax;m++){
  aveseg_crys+=m*m*pdf_crys[m];
   Lseg_crys+=m*pdf_crys[m];
   Lseg_tie+=m*Count_tie[m];
   Nseg_tie+=Count_tie[m];
    }
                              
  aveseg_crys= aveseg_crys/Lseg_crys;
 aveseg_amorph= aveseg_amorph/sumamorph;
  Lseg_crys= Lseg_crys/Nseg_crys;
 Lseg_tie= Lseg_tie/(1.0*Nseg_tie);
 Lseg_amorph= Lseg_amorph/Nseg_amorph;
   cout << "Lseg_tie" << Lseg_tie << " Nseg_tie=" << Nseg_tie << endl;  
// unsigned long  sumcount=0, sumcount1=0, sumcount2=0, sumcount3=0 ;



segfile1 << "Max_Nseg  Min_Nseg_amorph Max_Nseg_crys Max_Lseg_amorph    aveseg_crys aveseg_amorph   Lseg_crys Lseg_amorph   Lseg_tie ftie " << endl;
 segfile1 <<  *std::max_element(Nseg,Nseg+Nchain+1)  <<  " " << *std::min_element(Nseg,Nseg+Nchain+1)  << " " << segmax_crys << " "<< segmax_amorph << " "<< aveseg_crys << " "<< aveseg_amorph << " " <<  Lseg_crys << " "<< Lseg_amorph << " " << Lseg_tie << " " << ftie << endl;

strcpy(pdfseg,filename);
strcat(pdfseg, "_pdfseg.txt");
ofstream segfile;
segfile.open(pdfseg);
  
segfile <<  "m  countseg_amorph  countseg_crys  countseg_tie  pdfseg_amorph[m]  pdfseg_crys  pdfseg_tie" <<   endl;
for(m=1; m<=segmax;m++){
  //if(count[m]>0){
  segfile << m << " " << pdf_amorph[m] << " " << pdf_crys[m]<< " "<< (pdf_amorph[m]*1.0)/(Nseg_amorph*1.0) << "  "<< (pdf_crys[m]*1.0)/(1.0*Nseg_crys) << "  "<< (Count_tie[m]*1.0)/(1.0*Nseg_tie) << endl; }
//
 segfile << "the  average length of crystalline segments is="<< Lseg_crys << " "<< "of  amorphous segments is ="<< Lseg_amorph << endl;
 segfile << "the length weighted  average length of crystalline segments is="<< aveseg_crys << " "<< "of  amorphous segments is ="<< aveseg_amorph << endl;



 
 //**** test of segment recognition*****// 
 for(i=0; i<2;i++){

  cout << "i=" <<  i << " " << Nseg[i+1] << endl;
  }

 for(i=0; i<2;i++){
      sum_seg=0;
   j=1;
   while(seg[i][j]!=0){
	sum_seg+=abs(seg[i][j]);
	j++;
   }
 cout << "Nchain=" << i+1 << endl;
   cout << "i=" <<  i << "sum_seg= " << sum_seg << endl; 
for(n=1; n < Lchain;n++){
    id=i*Lchain+n;
    cout << bondamorph[id] << " ";
   
   }
 cout << "\n" <<  endl;
 }

 

delete Nseg;

for(int i = 0; i < Nchain; ++i)
delete [] seg[i];
delete [] seg;



/////////////////////////////////////////////////calculation of MSID within the amorphous domains ///////////////////////////////////////////////////////////




ofstream  MSIDfile_amorph;

strcpy(MSIDname,filename);
strcat(MSIDname, "_MSIDamorph.dat");

MSIDfile_amorph.open(MSIDname);
double Nr2[Lchain], MSIDamorph[Lchain],MSIDamorphx[Lchain],MSIDamorphy[Lchain],MSIDamorphz[Lchain];


int amorphcount[Lchain];

for(n=1;n< Lchain;n++){
 Nr2[n]=0;
amorphcount[n]=0;
MSIDamorphx[n]=0;
MSIDamorphy[n]=0;
MSIDamorphz[n]=0;

}

  for(i=0; i<Nchain;i++)
  for(n=1;n<Lchain;n++)
  for(j=1; j<= (Lchain-n);j++){
   id1=i*Lchain+j; // id of jth atom in the chain i
   id2=id1+n; // id of jth atom in the chain i

   if(amorphLabel[id1]==1 && amorphLabel[id2]==1){
   cnd[1]=xu[id1]-xu[id2];
   cnd[2]=yu[id1]-yu[id2];
   cnd[3]=zu[id1]-zu[id2];
   rr2=cnd[1]*cnd[1]+cnd[2]*cnd[2]+cnd[3]*cnd[3];
    Nr2[n]+=rr2;
    MSIDamorphx[n]+=cnd[1]*cnd[1];
    MSIDamorphy[n]+=cnd[2]*cnd[2];
    MSIDamorphz[n]+=cnd[3]*cnd[3];
   ++amorphcount[n];
}//end if
}
 
MSIDfile_amorph << " the MSID is calculated for " << Nchain << " chains of length "<< Lchain << endl;
MSIDfile_amorph << "n" << "  "<< "MSID(n)"<< "  "<< "MSID(n)/n"<<"  "<< "MSIDx(n)"<< "  "<< "MSIDy(n)"<< "  "<< "MSIDz(n)" << endl;

for(n=1;n< Lchain;n++){
 MSID[n]=Nr2[n]/amorphcount[n];
MSIDamorphx[n]=MSIDamorphx[n]/amorphcount[n];
MSIDamorphy[n]=MSIDamorphy[n]/amorphcount[n];
MSIDamorphz[n]=MSIDamorphz[n]/amorphcount[n];


 MSIDfile_amorph << n << " "<< MSID[n]<< " "<< MSID[n]/n  << " "<< MSIDamorphx[n] << " "<< MSIDamorphy[n] << " "<< MSIDamorphz[n] << endl;}


/////////////////////////////////////////////////calculation of MSID within the crystalline domains ///////////////////////////////////////////////////////////


ofstream  MSIDfile_crys;

strcpy(MSIDname,filename);
strcat(MSIDname, "_MSIDcrys.dat");

MSIDfile_crys.open(MSIDname);
double  MSIDcrys[Lchain],MSIDcrysx[Lchain],MSIDcrysy[Lchain],MSIDcrysz[Lchain];


int cryscount[Lchain];

for(n=1;n< Lchain;n++){
 Nr2[n]=0;
cryscount[n]=0;
MSIDcrysx[n]=0;
MSIDcrysy[n]=0;
MSIDcrysz[n]=0;

}

  for(i=0; i<Nchain;i++)
  for(n=1;n<Lchain;n++)
  for(j=1; j<= (Lchain-n);j++){
   id1=i*Lchain+j; // id of jth atom in the chain i
   id2=id1+n; // id of jth atom in the chain i

   if(amorphLabel[id1]!=1 && amorphLabel[id2]!=1){
   cnd[1]=xu[id1]-xu[id2];
   cnd[2]=yu[id1]-yu[id2];
   cnd[3]=zu[id1]-zu[id2];
   rr2=cnd[1]*cnd[1]+cnd[2]*cnd[2]+cnd[3]*cnd[3];
    Nr2[n]+=rr2;
    MSIDcrysx[n]+=cnd[1]*cnd[1];
    MSIDcrysy[n]+=cnd[2]*cnd[2];
    MSIDcrysz[n]+=cnd[3]*cnd[3];
   ++cryscount[n];
}//end if
}
 
MSIDfile_crys << " the MSID is calculated for " << Nchain << " chains of length "<< Lchain << endl;
MSIDfile_crys << "n" << "  "<< "MSID(n)"<< "  "<< "MSID(n)/n"<<"  "<< "MSIDx(n)"<< "  "<< "MSIDy(n)"<< "  "<< "MSIDz(n)" << endl;

for(n=1;n< Lchain;n++){
 MSID[n]=Nr2[n]/cryscount[n];
MSIDcrysx[n]=MSIDcrysx[n]/cryscount[n];
MSIDcrysy[n]=MSIDcrysy[n]/cryscount[n];
MSIDcrysz[n]=MSIDcrysz[n]/cryscount[n];


 MSIDfile_crys << n << " "<< MSID[n]<< " "<< MSID[n]/n  << " "<< MSIDcrysx[n] << " "<< MSIDcrysy[n] << " "<< MSIDcrysz[n] << endl;}


/*
//////////////////////////////////////////////////////////////////// calculation of pair correlations in amorph and ordered regions ///////////////////////////////////////////////



r2max=min(Lx,Lz)-1;

//r2max=0.25*r2max*r2max;

r2max=100;
double r2d[4];
int ngbin=600, ngbinr=600;
int ngy=600;
 Deltas= ngbinr/r2max;
int ig, igr, igy;
double ymax=6.0,rmax=6.0;
double rmax2;

rmax2=rmax*rmax;
double Deltay=ngy/ymax, Deltarho= ngbin/rmax;

cout << "rmax=" <<  sqrt(r2max) << " ngbin=" << ngbin << " 1/Deltar=" << 1.0/Deltarho << endl;
cout << "ymax=" <<  ymax << " ngy=" << ngy << " 1/Deltay=" << 1.0/Deltay << endl;


//long int Ngr[ngbin][ngy];
long int **Ngr; 
Ngr= new long int *[ngbin ];
for(int i = 0; i < ngbin; ++i){
Ngr[i] = new long int [ngy];
if (Ngr[i] == NULL) {cerr << " Ngr[i] Allocation problem!"; exit(1);}}

//  double Ag[ngbinr],Agamorph[ngbinr],Agcrys[ngbinr];



double *Ag , *Agamorph , *Agcrys ;
Ag= new double  [ngbinr];
Agamorph= new double[ngbinr];
Agcrys= new double  [ngbinr];

int *Ng, *Ngamorph, *Ngcrys;

Ng=new int [ngbinr];
Ngamorph=new int [ngbinr];
Ngcrys=new int [ngbinr];

//double  Agr[ngbin+1][ngy+1];
 double  **Agr; 
Agr= new double *[ngbin];
for(int i = 0; i < ngbin; ++i){
Agr[i] = new double [ngy];
if (Agr[i] == NULL) {cerr << " Agr[i] Allocation problem!"; exit(1);}}


for(i=0;i< ngbin;i++) 
for(j=0;j< ngy;j++)
{ Ngr[i][j]=0;
Agr[i][j]=0;
}


for(i=0;i< ngbinr;i++) 
{ Ng[i]=0;
Ag[i]=0.0;

Ngcrys[i]=0;
Agcrys[i]=0.0;

Ngamorph[i]=0;
Agamorph[i]=0.0;

}



//long int Ngrcrys[ngbin][ngy];
long int **Ngrcrys; 
Ngrcrys= new long int *[ngbin];
for(int i = 0; i < ngbin; ++i){
Ngrcrys[i] = new long int [ngy];
if (Ngrcrys[i] == NULL) {cerr << " Ngrcrys[i] Allocation problem!"; exit(1);}}


//double  Agrcrys[ngbin][ngy];
 double  **Agrcrys; 
Agrcrys= new double *[ngbin];
for(int i = 0; i < ngbin; ++i){
Agrcrys[i] = new double [ngy];
if (Agrcrys[i] == NULL) {cerr << " Agrcrys[i] Allocation problem!"; exit(1);}}


for(i=0;i< ngbin;i++) 
for(j=0;j< ngy;j++)
{ Ngrcrys[i][j]=0;
Agrcrys[i][j]=0;
}


//long int Ngramoprh[ngbin][ngy];
long int **Ngramorph; 
Ngramorph= new long int *[ngbin ];
for(int i = 0; i < ngbin; ++i){
Ngramorph[i] = new long int [ngy];
if (Ngr[i] == NULL) {cerr << " Ngramroph[i] Allocation problem!"; exit(1);}}


//double  Agramroph[ngbin][ngy];
 double  **Agramorph; 
Agramorph= new double *[ngbin];
for(int i = 0; i < ngbin; ++i){
Agramorph[i] = new double [ngy];
if (Agramorph[i] == NULL) {cerr << " Agramorph[i] Allocation problem!"; exit(1);}}


for(i=0;i< ngbin;i++) 
for(j=0;j< ngy;j++)
{ Ngramorph[i][j]=0;
Agramorph[i][j]=0.0;
}

  

cout << "Up to here fine" << endl;
     
for(n1=1; n1< Nchain; n1++)
for(n2=n1; n2<= Nchain; n2++)
for(i=1;i <= Lchain;i++)
  for(j=1;j <= Lchain;j++){
 id1=(n1-1)*Lchain+i;
  id2=(n2-1)*Lchain+j;

  
if(id1!=id2){
   rx=xp[id1]-xp[id2];
   ry=yp[id1]-yp[id2];
   rz=zp[id1]-zp[id2];



if (abs(rx) >  Lxhalf) {if(rx >0) rx = rx-Lx; 
                          else rx = rx+Lx;}

if (abs(ry) >   Lyhalf) {if(ry >0) ry = ry-Ly; 
                         else ry = ry+Ly;}

if (abs(rz) >   Lzhalf) {if(rz >0) rz = rz-Lz; 
                          else rz = rz+Lz;}


r2=rx*rx+ry*ry+rz*rz;


if(r2 < r2max  ){

 ig=(int)trunc(r2*Deltas);
  Ng[ig]+=2;


if(amorphLabel[id1]!=1 && amorphLabel[id2]!=1)   Ngcrys[ig]+=2;
  if(amorphLabel[id1]==1 && amorphLabel[id2]==1)   Ngamorph[ig]+=2;

}

rr2=rx*rx+rz*rz;

if(rr2 < rmax2 && abs(ry)< ymax ){
    
     igr=(int)trunc(sqrt(rr2)*Deltarho); 
     igy=(int)trunc(sqrt(ry*ry)*Deltay);
      Ngr[igr][igy]+=2;
    
  if(amorphLabel[id1]!=1 && amorphLabel[id2]!=1)   Ngrcrys[igr][igy]+=2;
  if(amorphLabel[id1]==1 && amorphLabel[id2]==1)   Ngramorph[igr][igy]+=2;
     }}}


double IDeltar2=1.0/Deltas,nid;
double IDeltay=1.0/Deltay;
 IDeltar=1.0/Deltarho;
double r2r[ngbinr], y[ngy];

ofstream  gfile;
char gfilename[50];
strcpy(gfilename,filename);
strcat(gfilename, "_gry.dat");

gfile.open(gfilename);

double amorphfraction=Namorph*1.0/Nmon;
for(j=0; j< ngy;j++)
for(i=0; i< ngbin;i++){
       //  myfile << "ig=" << ig << "g[ig]=" << g[ig] << endl;               
       y[i]=(i+0.50)*IDeltay;
       nid=(2*i+1)*pi*IDeltar*IDeltar*rho*2.0*IDeltay; // number of ideal gas particles in the shell
 Agr[i][j]=(Ngr[i][j]*1.0)/(Nmon * nid);
 Agramorph[i][j]=Ngramorph[i][j]/(Namorph * nid * amorphfraction);
Agrcrys[i][j]=Ngrcrys[i][j]/((Nmon-Namorph) * nid * (1.0-amorphfraction));
}
 



for(j=0; j< ngy;j++)
for(i=0; i< ngbin;i++) {
 gfile << Agr[i][j] << endl ;
}

ofstream  gfile1, gfile2,  gfile_amorph,gfile_crys;

strcpy(gfilename,filename);
strcat(gfilename, "_gr.dat");

gfile1.open(gfilename);

strcpy(gfilename,filename);
strcat(gfilename, "_gr-amorph.dat");

gfile_amorph.open(gfilename);

strcpy(gfilename,filename);
strcat(gfilename, "_gr-crys.dat");

gfile_crys.open(gfilename);

strcpy(gfilename,filename);
strcat(gfilename, "_gr-gy.dat");

gfile2.open(gfilename);



for(j=0; j< ngy;j++)
for(i=0; i< ngbin;i++) {
 gfile_crys << Agrcrys[i][j] << endl ;
}


for(j=0; j< ngy;j++)
for(i=0; i< ngbin;i++) {
 gfile_amorph << Agramorph[i][j] << endl ;
}






for(i=0; i< ngbinr;i++){
       //  myfile << "ig=" << ig << "g[ig]=" << g[ig] << endl; 
       //if(i==1) cout<< "Ng[1]=" << Ng[i] << endl;             
       r2r[i]=sqrt((i+0.5)*IDeltar2);
       nid=4.0/3.0 *pi*rho*(pow(i+1,1.5)-pow(i,1.5))*pow(IDeltar2,1.5 ); // number of ideal gas particles in the shell
       Ag[i]=Ng[i]*1.0/(Nmon * nid);
       Agamorph[i]=Ngamorph[i]*1.0/(Namorph * nid * amorphfraction);
       Agcrys[i]=Ngcrys[i]*1.0/((Nmon-Namorph)* nid * (1.0-amorphfraction));
}

 gfile2 << "i  r[i]   grho   grho-crys  grho-amorph  gy  gy-crys  gy-amorph " << endl ;

for(i=0; i<ngy;i++) gfile2 << i << " "<< y[i] << "  "<< Agr[i][0]  << " "<< Agrcrys[i][0]<< " "<< Agramorph[i][0] << " "<<  Agr[0][i] << " " << Agrcrys[0][i] << " " << Agramorph[0][i]  << endl ;

gfile1 << "i   Ng[i] r[i]   gr  gr-crys  gr-amorph " << endl ;
for(i=0; i<ngbinr;i++) gfile1 << i << " "<< Ng[i] <<" "<< r2r[i] << "  "<< Ag[i]  <<  "  "<< Agcrys[i] << "  "<< Agamorph[i]   << endl ;





delete xp,yp,zp;

delete  Ng, Ag,Agamorph, Agcrys,  Ngcrys, Ngamorph ;

for(i = 0; i < ngbin ; ++i){
delete []Agr[i];
delete []Agrcrys[i];
delete []Agramorph[i];
delete []Ngr[i];
delete []Ngrcrys[i];
delete []Ngramorph[i];
}
delete [] Agr;
delete []Agrcrys;
delete []Agramorph ;
delete []Ngr ;
delete []Ngrcrys ;
delete []Ngramorph;
*/
////////////////////////////////////////////   Cluster analysis     ///////////////////////////////////////////////////////////////
 int LL[gridcount+2], ni,nj; // label of labels
 for(i=0;i<=(gridcount+1);i++) LL[i]=0;
int gridtest=0;
 int min1; 
 int label=1;
/*cout<< "labels of first sheet" << endl;
for(n2=0;n2<ny;n2++){
for(n3=0;n3<nz;n3++)
cout << ncolor[0][n2][n3] << "  " ;
cout <<  "  "  << endl;}*/

if( ncolor[0][0][0]!=0) {
++gridtest;
++label;
++LL[label];
ncolor[0][0][0]=label; //first label will be 2.

}


/*******************checking the clusters within the first 1D-row **************************************/
double SC=0.97;
for(n3=1;n3<nz;n3++){
//cout << ncolor[0][0][n3] << "  " ;
if( ncolor[0][0][n3]!=0 ) { 
++gridtest;
if (ncolor[0][0][n3-1]==0)
{++label;
ncolor[0][0][n3]=label;
++LL[label];}
else { 
cosQ=nn1[0][0][n3]*nn1[0][0][n3-1]+nn2[0][0][n3]*nn2[0][0][n3-1]+nn3[0][0][n3]*nn3[0][0][n3-1];
P2Q=1.50* cosQ*cosQ-0.50;
if(P2Q < SC) ++label;
 ncolor[0][0][n3]=label;
 ++LL[label];
} //end else 
} //end if1
} //end for n3

// P.B.C
if( ncolor[0][0][0]!=0 && ncolor[0][0][nz-1]!=0 && ncolor[0][0][0]!=ncolor[0][0][nz-1] ) { // add the condition that and if they are not equal
cosQ=nn1[0][0][0]*nn1[0][0][nz-1]+nn2[0][0][0]*nn2[0][0][nz-1]+nn3[0][0][0]*nn3[0][0][nz-1];
if(P2Q >= SC){
i=ncolor[0][0][0];
j=ncolor[0][0][nz-1];
LL[i]+=LL[j];
LL[j]=-i;
ncolor[0][0][nz-1]=i;
}
}

//n1=0;
//cout<< "gridtest after first row=" << gridtest << endl;

/*for(n3=0;n3<nz;n3++) cout << ncolor[0][0][n3] << "  " ;

cout << endl;*/
/*******************checking the clusters within the first 2D-sheet **************************************/
for(n2=1;n2<ny;n2++){

if(ncolor[0][n2][0]!=0) { //if * checking the most left neighbor 
++gridtest;
if (ncolor[0][n2-1][0]==0){
++label;
ncolor[0][n2][0]=label;
++LL[label];}
else {  //else 1
cosQ=nn1[0][n2][0]*nn1[0][n2-1][0]+nn2[0][n2][0]*nn2[0][n2-1][0]+nn3[0][n2][0]*nn3[0][n2-1][0];
P2Q=1.50* cosQ*cosQ-0.50;

//if( n2==1 ) cout << "n1==0 && n2==1 && n3==0  P2Q=" << P2Q << endl;
if(P2Q < SC) { ++label;
 ncolor[0][n2][0]=label;
 ++LL[label];}
 else { //else 2 // add the case LL[i]<0
i=ncolor[0][n2-1][0];
while (LL[i]< 0 )   i=-LL[i];
ncolor[0][n2][0]=i;
 ++LL[i];


}//end else 2
} //end else 1
} //end if *



/////*******************checking the rows for each n2***
for(n3=1;n3<nz;n3++){
//cout << ncolor[0][n2][n3] << "  " ;
if(ncolor[0][n2][n3]!=0) {  // if1
++gridtest;
if (ncolor[0][n2][n3-1]==0 && ncolor[0][n2-1][n3]==0)
{++label;
ncolor[0][n2][n3]=label;
++LL[label];}

else{ // else *

if (ncolor[0][n2][n3-1]==0 && ncolor[0][n2-1][n3]!=0) {//lower neigbor only
cosQ1=nn1[0][n2][n3]*nn1[0][n2-1][n3]+nn2[0][n2][n3]*nn2[0][n2-1][n3]+nn3[0][n2][n3]*nn3[0][n2-1][n3];
P2Q1=1.50* cosQ1*cosQ1-0.50;



if(P2Q1 < SC){++label;
 ncolor[0][n2][n3]=label;
 ++LL[label];}

else {i=ncolor[0][n2-1][n3];
while (LL[i]<0) i=-LL[i];
ncolor[0][n2][n3]=i;
++LL[i];}


} //  end if lower neigbor only


else{ //else **
if (ncolor[0][n2][n3-1]!=0 && ncolor[0][n2-1][n3]==0) { //  if left neigbor only
cosQ=nn1[0][n2][n3]*nn1[0][n2][n3-1]+nn2[0][n2][n3]*nn2[0][n2][n3-1]+nn3[0][n2][n3]*nn3[0][n2][n3-1];
P2Q=1.50* cosQ*cosQ-0.50;

if(P2Q < SC){++label;
 ncolor[0][n2][n3]=label;
 ++LL[label];}

else {
i=ncolor[0][n2][n3-1];
while (LL[i]<0) i=-LL[i];
ncolor[0][n2][n3]=i;
++LL[i]; }


}//  end if left neigbor only 

else{ // else ***

if (ncolor[0][n2][n3-1]!=0 && ncolor[0][n2-1][n3]!=0){ //iff both left and lower neighbor exist


//check lower and right neighbors
 
cosQ=nn1[0][n2][n3]*nn1[0][n2][n3-1]+nn2[0][n2][n3]*nn2[0][n2][n3-1]+nn3[0][n2][n3]*nn3[0][n2][n3-1];
P2Q=1.50* cosQ*cosQ-0.50;

cosQ1=nn1[0][n2][n3]*nn1[0][n2-1][n3]+nn2[0][n2][n3]*nn2[0][n2-1][n3]+nn3[0][n2][n3]*nn3[0][n2-1][n3];
P2Q1=1.50* cosQ1*cosQ1-0.50;

if(P2Q < SC && P2Q1 < SC ){ ++label;
 ncolor[0][n2][n3]=label;
 ++LL[label];}

else{ //els s1
if(P2Q >= SC  && P2Q1 >= SC ) {
i=ncolor[0][n2][n3-1];
j=ncolor[0][n2-1][n3];

while (LL[i]<0) i=-LL[i];
while (LL[j]<0) j=-LL[j];

min1=min(i,j) ;
ncolor[0][n2][n3]=min1;
++LL[min1];

if (i > j) {LL[j]+=LL[i]; LL[i]=-j; } 
if (i < j) {LL[i]+=LL[j]; LL[j]=-i; } 
}
 


else{ //els s2
if(P2Q >= SC ) {
i=ncolor[0][n2][n3-1];
while (LL[i]< 0) i=-LL[i];
ncolor[0][n2][n3]=i;
++LL[i];}


else {// equiv. if(P2Q1 >= SC )    else s3 
i=ncolor[0][n2-1][n3];
while (LL[i]< 0) i=-LL[i];
ncolor[0][n2][n3]=i;
++LL[i]; 
}// end els s3

} // end els s2
} // end els s1

} //end iff

} //end els***    
} //end els**
} //end els*  


} //end if1



}//end for n3              
//P.B.C.  for n3 direction

if( ncolor[0][n2][0]!=0 && ncolor[0][n2][nz-1]!=0 && ncolor[0][n2][nz-1]!=ncolor[0][n2][0]) { //P.B.C. for n3 direction
cosQ=nn1[0][n2][0]*nn1[0][n2][nz-1]+nn2[0][n2][0]*nn2[0][n2][nz-1]+nn3[0][n2][0]*nn3[0][n2][nz-1];
if(P2Q >= SC){
i=ncolor[0][n2][0];
j=ncolor[0][n2][nz-1];

while (LL[i]<0) i=-LL[i];
while (LL[j]<0) j=-LL[j];
if (i < j) {  LL[i]+=LL[j]; LL[j]=-i; ncolor[0][n2][nz-1]=i;}
if (j < i) {  LL[j]+=LL[i]; LL[i]=-j; ncolor[0][n2][0]=j; }
}
}

} //end for n2
//P.B.C.  for n2 direction
for(n3=0;n3 < nz;n3++){
if( ncolor[0][0][n3]!=0 && ncolor[0][ny-1][n3]!=0 && ncolor[0][0][n3]!=ncolor[0][ny-1][n3]) { 
cosQ=nn1[0][0][n3]*nn1[0][ny-1][n3]+nn2[0][0][n3]*nn2[0][ny-1][n3]+nn3[0][0][n3]*nn3[0][ny-1][n3];
if(P2Q >= SC){ //cout << "P.B in n2 direction " << endl;
i=ncolor[0][0][n3];
j=ncolor[0][ny-1][n3];
while (LL[i]<0) i=-LL[i];
while (LL[j]<0) j=-LL[j];
//cout << "i=" << 2<<  " LL[2]=" << LL[2] << endl; // condition of LL[i] & LL[j] both positive should be added!!
if (i < j) {  LL[i]+=LL[j]; LL[j]=-i; ncolor[0][ny-1][n3]=i; }
if (j < i ) {  LL[j]+=LL[i]; LL[i]=-j; ncolor[0][0][n3]=j; }
}
}
}//end for n3*/

/*for(i=1;i<=label;i++){
             cout << "i=" << i<<  " LL[i]=" << LL[i] << endl; }

cout << "label value after first sheet is" << label << endl;
for(n2=0;n2<ny;n2++){
for(n3=0;n3<nz;n3++)
cout << ncolor[0][n2][n3] << "  " ;
cout <<  "  " << endl;}*/
cout<< "gridtest after first sheet=" << gridtest << endl;
////////////////////////////////////////////////////////////////////////////////////////////////////
int  dn=0,labelcount=0;
 

/*cout<< " labels of 5th sheet" << endl;
for(n2=0;n2<ny;n2++){
for(n3=0;n3<nz;n3++)
cout << ncolor[4][n2][n3] << "  " ;
cout <<  "  "  << endl;}*/
/*******************checking the clusters between 2D-sheets corresponding to different  n1 s **************************************/

for(n1=1;n1<nx;n1++){

// first point of sheet n1
if( ncolor[n1][0][0]!=0) { //if 0
++gridtest;
if (ncolor[n1-1][0][0]==0){
++label;
ncolor[n1][0][0]=label;
++LL[label];}
else {  //else 1
cosQ=nn1[n1][0][0]*nn1[n1-1][0][0]+nn2[n1][0][0]*nn2[n1-1][0][0]+nn3[n1][0][0]*nn3[n1-1][0][0];
P2Q=1.50* cosQ*cosQ-0.50;
if(P2Q < SC) { ++label;
 ncolor[n1][0][0]=label;
 ++LL[label];}
 else { //else 2
i=ncolor[n1-1][0][0];
while (LL[i]<0) i=-LL[i];
ncolor[n1][0][0]=i;
 ++LL[i];
}//end else 2
} //end else 1
} //end if 0

// *************  first row (n2=0) of sheet n1 .i.e. [n1][0][n3]******************

for(n3=1;n3<nz;n3++){

if(ncolor[n1][0][n3]!=0) {  // if1
++gridtest;
if (ncolor[n1][0][n3-1]==0 && ncolor[n1-1][0][n3]==0)
{++label;
ncolor[n1][0][n3]=label;
++LL[label];}

else{ // else *

if (ncolor[n1][0][n3-1]==0 && ncolor[n1-1][0][n3]!=0) {//lower n1 neigbor only
cosQ2=nn1[n1][0][n3]*nn1[n1-1][0][n3]+nn2[n1][0][n3]*nn2[n1-1][0][n3]+nn3[n1][0][n3]*nn3[n1-1][0][n3];
P2Q2=1.50* cosQ2*cosQ2-0.50;
if(P2Q2 < SC){ ++label;
 ncolor[n1][0][n3]=label;
 ++LL[label];}

else {
i=ncolor[n1-1][0][n3];
while (LL[i]<0) i=-LL[i];
ncolor[n1][0][n3]=i;
++LL[i];
} 

} //  end if lower neigbor only

else{ //else **
if (ncolor[n1][0][n3-1]!=0 && ncolor[n1-1][0][n3]==0) { //  if left neigbor only
cosQ=nn1[n1][0][n3]*nn1[n1][0][n3-1]+nn2[n1][0][n3]*nn2[n1][0][n3-1]+nn3[n1][0][n3]*nn3[n1][0][n3-1];
P2Q=1.50* cosQ*cosQ-0.50;

if(P2Q < SC){++label;
 ncolor[n1][0][n3]=label;
 ++LL[label];}

else {
i=ncolor[n1][0][n3-1];
while (LL[i]<0) i=-LL[i];
ncolor[n1][0][n3]=i;
++LL[i];
} //end else

}//  end if left neigbor only 

else{ // else ***

if (ncolor[n1][0][n3-1]!=0 && ncolor[n1-1][0][n3]!=0){ //iff both left and lower neighbor exist


//check lower and right neighbors
 
cosQ=nn1[n1][0][n3]*nn1[n1][0][n3-1]+nn2[n1][0][n3]*nn2[n1][0][n3-1]+nn3[n1][0][n3]*nn3[n1][0][n3-1];
P2Q=1.50* cosQ*cosQ-0.50;

cosQ2=nn1[n1][0][n3]*nn1[n1-1][0][n3]+nn2[n1][0][n3]*nn2[n1-1][0][n3]+nn3[n1][0][n3]*nn3[n1-1][0][n3];
P2Q2=1.50* cosQ2*cosQ2-0.50;

if(P2Q < SC && P2Q2 < SC ){ ++label;
 ncolor[n1][0][n3]=label;
 ++LL[label];}

else{ //els s1
if(P2Q >= SC  && P2Q2 >= SC ) {
i=ncolor[n1][0][n3-1];
j=ncolor[n1-1][0][n3];
while (LL[i]<0) i=-LL[i];
while (LL[j]<0) j=-LL[j];
min1=min(i,j) ;
ncolor[n1][0][n3]=min1;
++LL[min1];
if (j < i ) {  LL[j]+=LL[i]; LL[i]=-j;}
if (i< j) {  LL[i]+=LL[j]; LL[j]=-i;}
}


else{ //els s2
if(P2Q >= SC ) {
i=ncolor[n1][0][n3-1];
while (LL[i]<0) i=-LL[i];
ncolor[n1][0][n3]=i;
++LL[i];}

else {// equiv. if(P2Q1 >= SC )    else s3 
i=ncolor[n1-1][0][n3];
while (LL[i]<0) i=-LL[i];
ncolor[n1][0][n3]=i;
++LL[i]; 
}// end els s3

} // end els s2
} // end els s1

} //end iff

} //end els***
} //end els**
} //end els*


} //end if1

}//end for n3
//P.B.C. for n3 direction in the first row of sheet n1
if( ncolor[n1][0][0]!=0 && ncolor[n1][0][nz-1]!=0) { 
cosQ=nn1[n1][0][0]*nn1[n1][0][nz-1]+nn2[n1][0][0]*nn2[n1][0][nz-1]+nn3[n1][0][0]*nn3[n1][0][nz-1];
if(P2Q >= SC){
i=ncolor[n1][0][0];
j=ncolor[n1][0][nz-1];
while (LL[i]<0) i=-LL[i];
while (LL[j]<0) j=-LL[j];
min1=min(i,j);
if (i <j) {  LL[i]+=LL[j]; LL[j]=-i; ncolor[n1][0][nz-1]=i; }
if (j <i) {  LL[j]+=LL[i]; LL[i]=-j; ncolor[n1][0][0]=j;}
}
}



// *********************  row n2 of sheet n1 .i.e. [n1][n2][n3]******************

for(n2=1;n2<ny;n2++) {
 
if(ncolor[n1][n2][0]!=0) { //if 0, checking the most left neighbor of row n2 in sheet n1, its neigbors [n1][n2-1][0] & [n1-1][n2][0]
++gridtest;

if (ncolor[n1][n2-1][0]==0 && ncolor[n1-1][n2][0]==0){

++label;
ncolor[n1][n2][0]=label;
++LL[label];}

else{ // else *

if (ncolor[n1][n2-1][0]!=0 && ncolor[n1-1][n2][0]==0) {//lower n2 neigbor only
cosQ1=nn1[n1][n2][0]*nn1[n1][n2-1][0]+nn2[n1][n2][0]*nn2[n1][n2-1][0]+nn3[n1][n2][0]*nn3[n1][n2-1][0];
P2Q1=1.50* cosQ1*cosQ1-0.50;
if(P2Q1 < SC){++label;
ncolor[n1][n2][0]=label;
 ++LL[label];}

else {
i=ncolor[n1][n2-1][0];
while (LL[i]<0) i=-LL[i];
ncolor[n1][n2][0]=i;
++LL[i];
} 

} //  end if lower n2 neigbor only

else{ //else **
if (ncolor[n1][n2-1][0]==0 && ncolor[n1-1][n2][0]!=0) { //  if lower n1 neigbor only

cosQ2=nn1[n1][n2][0]*nn1[n1-1][n2][0]+nn2[n1][n2][0]*nn2[n1-1][n2][0]+nn3[n1][n2][0]*nn3[n1-1][n2][0];
P2Q2=1.50* cosQ2*cosQ2-0.50;

if(P2Q2 < SC){++label;
 ncolor[n1][n2][0]=label;
 ++LL[label];}

else {
i=ncolor[n1-1][n2][0];
while (LL[i]<0) i=-LL[i];
ncolor[n1][n2][0]=i;
++LL[i];
} //end else

}//  end if left neigbor only 

else{ // else *** equivalent to if (ncolor[n1-1][n2][0]!=0 && ncolor[n1][n2-1][0]!=0){ //iff both  lower neighbors exist


//check lower and right neighbors
 
cosQ1=nn1[n1][n2][0]*nn1[n1][n2-1][0]+nn2[n1][n2][0]*nn2[n1][n2-1][0]+nn3[n1][n2][0]*nn3[n1][n2-1][0];
P2Q1=1.50* cosQ1*cosQ1-0.50;

cosQ2=nn1[n1][n2][0]*nn1[n1-1][n2][0]+nn2[n1][n2][0]*nn2[n1-1][n2][0]+nn3[n1][n2][0]*nn3[n1-1][n2][0];
P2Q2=1.50* cosQ2*cosQ2-0.50;


if(P2Q1 < SC && P2Q2 < SC ){ ++label;
 ncolor[n1][n2][0]=label;
 ++LL[label];}


else{ //els s1
if(P2Q1 >= SC  && P2Q2 >= SC ) {
i=ncolor[n1-1][n2][0];
j=ncolor[n1][n2-1][0];
while (LL[i]<0) i=-LL[i];
while (LL[j]<0) j=-LL[j];
min1=min(i,j) ;
ncolor[n1][n2][0]=min1;
++LL[min1];
if (j <i) {  LL[j]+=LL[i]; LL[i]=-j;}
if (i <j) {  LL[i]+=LL[j]; LL[j]=-i;}
}


else{ //els s2
if(P2Q1 >= SC ) {
i=ncolor[n1][n2-1][0];
while (LL[i]<0) i=-LL[i];
ncolor[n1][n2][0]=i;
++LL[i];}

else {// equiv. if(P2Q2 >= SC )    else s3 
i=ncolor[n1-1][n2][0];
while (LL[i]<0) i=-LL[i];
ncolor[n1][n2][0]=i;
++LL[i]; 
}// end els s3

} // end els s2
} // end els s1


} //end els***
} //end els**
} //end els*



} //end if 0
/////*******************checking the rows for each n2***

//from here should be checked!
for(n3=1;n3<nz;n3++){

if(ncolor[n1][n2][n3]!=0) {  // if1
++gridtest;
if (ncolor[n1-1][n2][n3]==0 && ncolor[n1][n2-1][n3]==0 && ncolor[n1][n2][n3-1]==0)
{++label;
ncolor[n1][n2][n3]=label;
++LL[label];}

else{ // else 0

if (ncolor[n1-1][n2][n3]!=0 && ncolor[n1][n2-1][n3]==0 && ncolor[n1][n2][n3-1]==0) {//lower neigbor only
cosQ2=nn1[n1][n2][n3]*nn1[n1-1][n2][n3]+nn2[n1][n2][n3]*nn2[n1-1][n2][n3]+nn3[n1][n2][n3]*nn3[n1-1][n2][n3];
P2Q2=1.50* cosQ2*cosQ2-0.50;

//if(n1==2 && n2==6 && n3==5) cout << " n1==2 && n2==6 && n3==5,  P2Q=" << P2Q2 << endl;
if(P2Q2 < SC){++label;
 ncolor[n1][n2][n3]=label;
 ++LL[label];}

else {
i=ncolor[n1-1][n2][n3];
while (LL[i]<0) i=-LL[i];
ncolor[n1][n2][n3]=i;
++LL[i];
} 

} //  end if lower n1 neigbor only

else{ // else *

if (ncolor[n1-1][n2][n3]==0  && ncolor[n1][n2-1][n3]!=0 && ncolor[n1][n2][n3-1]==0) {//lower neigbor only
cosQ1=nn1[n1][n2][n3]*nn1[n1][n2-1][n3]+nn2[n1][n2][n3]*nn2[n1][n2-1][n3]+nn3[n1][n2][n3]*nn3[n1][n2-1][n3];
P2Q1=1.50* cosQ1*cosQ1-0.50;

//if(n1==1 && n2==4 && n3==2) cout << "n1==1 && n2==4 && n3==2  P2Q1=" << P2Q1 << endl;
if(P2Q1 < SC){++label;
 ncolor[n1][n2][n3]=label;
 ++LL[label];}

else {
i=ncolor[n1][n2-1][n3];
while (LL[i]<0) i=-LL[i];
ncolor[n1][n2][n3]=i;
++LL[i];
} 

} //  end if lower n2 neigbor only

else{ //else ** 

if (ncolor[n1-1][n2][n3]==0 && ncolor[n1][n2-1][n3]==0 && ncolor[n1][n2][n3-1]!=0 ) { //  if left neigbor only

cosQ=nn1[n1][n2][n3]*nn1[n1][n2][n3-1]+nn2[n1][n2][n3]*nn2[n1][n2][n3-1]+nn3[n1][n2][n3]*nn3[n1][n2][n3-1];
P2Q=1.50* cosQ*cosQ-0.50;


if(P2Q < SC){++label;
 ncolor[n1][n2][n3]=label;
 ++LL[label];}

else { 
i=ncolor[n1][n2][n3-1];
while (LL[i]<0) i=-LL[i];
ncolor[n1][n2][n3]=i;
++LL[i];
} //end else

}//  end if left  n3 neigbor only 

else{ // else ***

if (ncolor[n1-1][n2][n3]==0 && ncolor[n1][n2-1][n3]!=0 && ncolor[n1][n2][n3-1]!=0){ //iff both left and lower n2 neighbor exist


//check lower and right neighbors
 
cosQ1=nn1[n1][n2][n3]*nn1[n1][n2-1][n3]+nn2[n1][n2][n3]*nn2[n1][n2-1][n3]+nn3[n1][n2][n3]*nn3[n1][n2-1][n3];
P2Q1=1.50* cosQ1*cosQ1-0.50;

cosQ=nn1[n1][n2][n3]*nn1[n1][n2][n3-1]+nn2[n1][n2][n3]*nn2[n1][n2][n3-1]+nn3[n1][n2][n3]*nn3[n1][n2][n3-1];
P2Q=1.50* cosQ*cosQ-0.50;

if(P2Q < SC && P2Q1 < SC ){ ++label;
 ncolor[n1][n2][n3]=label;
 ++LL[label];}

else{ //els s1
if(P2Q >= SC  && P2Q1 >= SC ) {
i=ncolor[n1][n2][n3-1];
j=ncolor[n1][n2-1][n3];
while (LL[i]<0) i=-LL[i];
while (LL[j]<0) j=-LL[j];
min1=min(i,j) ;
ncolor[n1][n2][n3]=min1;
++LL[min1];
if (j <i) {  LL[j]+=LL[i]; LL[i]=-j;}
if (i< j) {  LL[i]+=LL[j]; LL[j]=-i;}
}


else{ //els s2
if(P2Q >= SC ) {
i=ncolor[n1][n2][n3-1];
while (LL[i]<0) i=-LL[i];
ncolor[n1][n2][n3]=i;
++LL[i];}

else {// equiv. if(P2Q1 >= SC )    else s3 
i=ncolor[n1][n2-1][n3];
while (LL[i]<0) i=-LL[i];
ncolor[n1][n2][n3]=i;
++LL[i]; 
}// end els s3

} // end els s2
} // end els s1

} //end iff

else{ // else ****

if (ncolor[n1-1][n2][n3]!=0 && ncolor[n1][n2-1][n3]==0 && ncolor[n1][n2][n3-1]!=0){ //iff both left and lower n2 neighbor exist


//check lower and right neighbors
 
cosQ2=nn1[n1][n2][n3]*nn1[n1-1][n2][n3]+nn2[n1][n2][n3]*nn2[n1-1][n2][n3]+nn3[n1][n2][n3]*nn3[n1-1][n2][n3];
P2Q2=1.50* cosQ2*cosQ2-0.50;

cosQ=nn1[n1][n2][n3]*nn1[n1][n2][n3-1]+nn2[n1][n2][n3]*nn2[n1][n2][n3-1]+nn3[n1][n2][n3]*nn3[n1][n2][n3-1];
P2Q=1.50* cosQ*cosQ-0.50;


if(P2Q < SC && P2Q2 < SC ){ ++label;
 ncolor[n1][n2][n3]=label;
 ++LL[label];}

else{ //els s1 
if(P2Q >= SC  && P2Q2 >= SC ) { 
//if(n1==3 && n2==8 && n3==6) cout << "n1==3 && n2==8 && n3==5  P2Q=" << P2Q  << "  P2Q2=" << P2Q2 << "  ncolor[n1-1]=" << ncolor[n1-1][n2][n3] << "  ncolor[n3-1]=" <<  ncolor[n1][n2][n3-1]<< endl;

i=ncolor[n1][n2][n3-1];
j=ncolor[n1-1][n2][n3];
while (LL[i]<0) i=-LL[i];
while (LL[j]<0) j=-LL[j];

//if(n1==3 && n2==8 && n3==6) cout << "i= " << i << " j= " << j << endl;
min1=min(i,j);
ncolor[n1][n2][n3]=min1;
++LL[min1];
if (j <i) {  LL[j]+=LL[i]; LL[i]=-j;}
if (i< j) {  LL[i]+=LL[j]; LL[j]=-i;}

//if(n1==3 && n2==8 && n3==6) cout << "i= " << i << " LL[i]= " << LL[i] << endl;
}


else{ //els s2
if(P2Q >= SC ) {
i=ncolor[n1][n2][n3-1];
while (LL[i]<0) i=-LL[i];
ncolor[n1][n2][n3]=i;
++LL[i];}

else {// equiv. if(P2Q2 >= SC )    else s3 
i=ncolor[n1-1][n2][n3];
while (LL[i]<0) i=-LL[i];
ncolor[n1][n2][n3]=i;
++LL[i]; 
}// end els s3

} // end els s2
} // end els s1

//if(n1==4 && n2==8 && n3==5) cout << "n1==4 && n2==8 && n3==5  ncolor=" << ncolor[n1][n2][n3] << endl;
} //end iff


else{ // else *****

if (ncolor[n1-1][n2][n3]!=0 && ncolor[n1][n2-1][n3]!=0 && ncolor[n1][n2][n3-1]==0){ //iff both left and lower n2 neighbor exist


//check lower and right neighbors

cosQ1=nn1[n1][n2][n3]*nn1[n1][n2-1][n3]+nn2[n1][n2][n3]*nn2[n1][n2-1][n3]+nn3[n1][n2][n3]*nn3[n1][n2-1][n3];
P2Q1=1.50* cosQ1*cosQ1-0.50; 

cosQ2=nn1[n1][n2][n3]*nn1[n1-1][n2][n3]+nn2[n1][n2][n3]*nn2[n1-1][n2][n3]+nn3[n1][n2][n3]*nn3[n1-1][n2][n3];
P2Q2=1.50* cosQ2*cosQ2-0.50;



if(P2Q1 < SC && P2Q2 < SC ){ ++label;
 ncolor[n1][n2][n3]=label;
 ++LL[label];}

else{ //els s1
if(P2Q1 >= SC  && P2Q2 >= SC ) {
i=ncolor[n1][n2-1][n3];
j=ncolor[n1-1][n2][n3];
while (LL[i]<0) i=-LL[i];
while (LL[j]<0) j=-LL[j];
min1=min(i,j) ;
ncolor[n1][n2][n3]=min1;
++LL[min1];
if (j <i) {  LL[j]+=LL[i]; LL[i]=-j;}
if (i< j) {  LL[i]+=LL[j]; LL[j]=-i;}
}


else{ //els s2
if(P2Q1 >= SC ) {
i=ncolor[n1][n2-1][n3];
while (LL[i]<0) i=-LL[i];
ncolor[n1][n2][n3]=i;
++LL[i];}

else {// equiv. if(P2Q2 >= SC )    else s3 
i=ncolor[n1-1][n2][n3];
while (LL[i]<0) i=-LL[i];
ncolor[n1][n2][n3]=i;
++LL[i]; 
}// end els s3

} // end els s2
} // end els s1

} //end iff

else{ // else 6*

if (ncolor[n1-1][n2][n3]!=0 && ncolor[n1][n2-1][n3]!=0 && ncolor[n1][n2][n3-1]!=0){ //iff all the left and lower n2 & n1 neighbors exist

cosQ=nn1[n1][n2][n3]*nn1[n1][n2][n3-1]+nn2[n1][n2][n3]*nn2[n1][n2][n3-1]+nn3[n1][n2][n3]*nn3[n1][n2][n3-1];
P2Q=1.50* cosQ*cosQ-0.50;


cosQ1=nn1[n1][n2][n3]*nn1[n1][n2-1][n3]+nn2[n1][n2][n3]*nn2[n1][n2-1][n3]+nn3[n1][n2][n3]*nn3[n1][n2-1][n3];
P2Q1=1.50* cosQ1*cosQ1-0.50; 

cosQ2=nn1[n1][n2][n3]*nn1[n1-1][n2][n3]+nn2[n1][n2][n3]*nn2[n1-1][n2][n3]+nn3[n1][n2][n3]*nn3[n1-1][n2][n3];
P2Q2=1.50* cosQ2*cosQ2-0.50;




if(P2Q1 >= SC && P2Q2>= SC  && P2Q >= SC) { //  if all 3
i=ncolor[n1][n2][n3-1];
j=ncolor[n1][n2-1][n3];
k=ncolor[n1-1][n2][n3];

while (LL[i]<0) i=-LL[i];
while (LL[j]<0) j=-LL[j];
while (LL[k]<0) k=-LL[k];

//if(n1==1 && n2==4 && n3==8) cout << "n1==1 && n2==4 && n3==8  P2Q=" << P2Q << "  P2Q1=" << P2Q1 << " P2Q2=" << P2Q2 << endl;



if (i==j && j==k) {  ncolor[n1][n2][n3]=i;   ++LL[i]; }

if (i==j && j!=k) { 

 if(i< k) { ncolor[n1][n2][n3]=i;   ++LL[i]; LL[i]+=LL[k]; LL[k]=-i; }
 if(i> k) { ncolor[n1][n2][n3]=k;   ++LL[k]; LL[k]+=LL[i]; LL[i]=-k; }
 
}


if (i==k && j!=k) { 

 if(i< j) { ncolor[n1][n2][n3]=i;   ++LL[i]; LL[i]+=LL[j]; LL[j]=-i;}
 if(j< i) { ncolor[n1][n2][n3]=j;   ++LL[j]; LL[j]+=LL[i]; LL[i]=-j;}
 
}

if (j==k && j!=i) { 

 if(i< j) { ncolor[n1][n2][n3]=i;   ++LL[i]; LL[i]+=LL[j]; LL[j]=-i; }
 if(j< i) { ncolor[n1][n2][n3]=j;   ++LL[j]; LL[j]+=LL[i]; LL[i]=-j; }
 
}


if (j!=k && j!=i && i!=k) { 

ni=min(i,j);
min1=min(ni,k);
ncolor[n1][n2][n3]=min1;
++LL[min1];

 if(i==min1) {    LL[i]+=LL[j]+LL[k]; LL[j]=-i; LL[k]=-i; }
 if(j==min1) {    LL[j]+=LL[i]+LL[k]; LL[i]=-j; LL[k]=-j; }
 if(k==min1) {    LL[k]+=LL[j]+LL[i]; LL[j]=-k; LL[i]=-k; }
}




}// end if all 3

else{ //else s1
if(P2Q1 >= SC  && P2Q2 >= SC  && P2Q < SC) {
i=ncolor[n1][n2-1][n3];
j=ncolor[n1-1][n2][n3];
while (LL[i]<0) i=-LL[i];
while (LL[j]<0) j=-LL[j];
min1=min(i,j) ;
ncolor[n1][n2][n3]=min1;
++LL[min1];
if (j <i) {  LL[j]+=LL[i]; LL[i]=-j;}
if (i< j) {  LL[i]+=LL[j]; LL[j]=-i;}
}


else{ //else s2
if(P2Q1 >= SC && P2Q2 < SC  && P2Q >= SC) {
i=ncolor[n1][n2-1][n3];
j=ncolor[n1][n2][n3-1];
while (LL[i]<0) i=-LL[i];
while (LL[j]<0) j=-LL[j];
min1=min(i,j) ;
ncolor[n1][n2][n3]=min1;
++LL[min1];
if (j <i) {  LL[j]+=LL[i]; LL[i]=-j;}
if (i< j) {  LL[i]+=LL[j]; LL[j]=-i;}
}

else{ //else s3
if(P2Q1 < SC && P2Q2 >= SC  && P2Q >= SC) {
i=ncolor[n1-1][n2][n3];
j=ncolor[n1][n2][n3-1];
while (LL[i]<0) i=-LL[i];
while (LL[j]<0) j=-LL[j];
min1=min(i,j) ;
ncolor[n1][n2][n3]=min1;
++LL[min1];
if (j <i) {  LL[j]+=LL[i]; LL[i]=-j;}
if (i< j) {  LL[i]+=LL[j]; LL[j]=-i;}
}

else {    //else s4 
if(P2Q1 < SC && P2Q2>= SC  && P2Q < SC) {
i=ncolor[n1-1][n2][n3];
while (LL[i]<0) i=-LL[i];
ncolor[n1][n2][n3]=i;
++LL[i]; }

else {    //else s5 
if(P2Q1 >= SC && P2Q2< SC  && P2Q < SC) {
i=ncolor[n1][n2-1][n3];
while (LL[i]<0) i=-LL[i];
ncolor[n1][n2][n3]=i;
++LL[i]; }

else {    //else s6 
if(P2Q1 < SC && P2Q2< SC  && P2Q >= SC) {
i=ncolor[n1][n2][n3-1];
while (LL[i]<0) i=-LL[i];
ncolor[n1][n2][n3]=i;
++LL[i]; }

else {    //else s7 

if(P2Q1 < SC && P2Q2 < SC   && P2Q < SC ){ ++label;
 ncolor[n1][n2][n3]=label;
 ++LL[label];}






}// end els s7

}// end els s6
}// end els s5
}// end els s4
}// end els s3
} // end els s2
} // end els s1

} //end iff

}//end else 6*
}//end els*****
}//end els****
} //end els***
} //end els**
} //end els*

} //end els 0
} //end if1


}//end for n3

//P.B.C. for n3 direction
if( ncolor[n1][n2][0]!=0 && ncolor[n1][n2][nz-1]!=0 && ncolor[n1][n2][nz-1]!=ncolor[n1][n2][0]) { 
cosQ=nn1[n1][n2][0]*nn1[n1][n2][nz-1]+nn2[n1][n2][0]*nn2[n1][n2][nz-1]+nn3[n1][n2][0]*nn3[n1][n2][nz-1];
if(P2Q >= SC){
i=ncolor[n1][n2][0];
j=ncolor[n1][n2][nz-1];
while (LL[i]<0) i=-LL[i];
while (LL[j]<0) j=-LL[j];
if (i < j) {  LL[i]+=LL[j]; LL[j]=-i; ncolor[n1][n2][nz-1]=i;}
if (j < i) {  LL[j]+=LL[i]; LL[i]=-j; ncolor[n1][n2][0]=j; }
}
}

} //end for n2

//P.B.C.  for n2 direction

for(n3=0;n3 < nz;n3++){
if( ncolor[n1][0][n3]!=0 && ncolor[n1][ny-1][n3]!=0 && ncolor[n1][0][n3]!=ncolor[n1][ny-1][n3]) { 
cosQ=nn1[n1][0][n3]*nn1[n1][ny-1][n3]+nn2[n1][0][n3]*nn2[n1][ny-1][n3]+nn3[n1][0][n3]*nn3[n1][ny-1][n3];
if(P2Q >= SC){ //cout << "P.B in n2 direction " << endl;
i=ncolor[n1][0][n3];
j=ncolor[n1][ny-1][n3];
while (LL[i]<0) i=-LL[i] ;
while (LL[j]<0) j=-LL[j];
//cout << "i=" << 2<<  " LL[2]=" << LL[2] << endl; // condition of LL[i] & LL[j] both positive should be added!!
if (i < j) {  LL[i]+=LL[j]; LL[j]=-i;  ncolor[n1][ny-1][n3]=i;}
if (j < i ) {  LL[j]+=LL[i]; LL[i]=-j; ncolor[n1][0][n3]=j; }
}
}
}//end for n3*/



/*if(n1==2) {cout<< "gridtest after 3rd sheet, i.e n1=4 is:  " << gridtest << endl;

for (i=0;i<=(gridcount+1);i++) acount[i]=0;
// The relevant labels start from 2

// The total number of independent crystalline domains is given by labelcount


for(i=1;i<=label;i++)
cout << "i=" << i<<  " LL[i]=" << LL[i] << endl;}
if(LL[i]>0) { ++labelcount;
              //cout << "i=" << i<<  " LL[i]=" << LL[i] << endl;
             acount[0]+=LL[i];
 }

}

cout << " The total number of independent crystalline domains is:  " << labelcount << endl;

cout << " The total number of  crystalline grid elements is:  " << acount[0] << endl;
}*/
} // end for n1


//P.B.C.  for n1 direction

for(n2=0;n2 < ny;n2++)
for(n3=0;n3 < nz;n3++){
if( ncolor[0][n2][n3]!=0 && ncolor[nx-1][n2][n3]!=0 && ncolor[0][n2][n3]!=ncolor[nx-1][n2][n3] ) { 
cosQ=nn1[0][n2][n3]*nn1[nx-1][n2][n3]+nn2[0][n2][n3]*nn2[nx-1][n2][n3]+nn3[0][n2][n3]*nn3[nx-1][n2][n3];
if(P2Q >= SC){ //cout << "P.B in n1 direction " << endl;
i=ncolor[0][n2][n3];
j=ncolor[nx-1][n2][n3];
while (LL[i]<0) i=-LL[i];
while (LL[j]<0) j=-LL[j];
//cout << "i=" << 2<<  " LL[2]=" << LL[2] << endl; // condition of LL[i] & LL[j] both positive should be added!!
if (i < j) {  LL[i]+=LL[j]; LL[j]=-i; ncolor[nx-1][n2][n3]=i; }
if (j < i ) {  LL[j]+=LL[i]; LL[i]=-j; ncolor[0][n2][n3]=j; }
}
}
}//end for n3*/


cout << " The final label value  is" << label << endl;

/*cout << " The final labels of 4th sheet:" <<  endl;
for(n2=0;n2<ny;n2++){
for(n3=0;n3<nz;n3++)
cout << ncolor[3][n2][n3] << "  " ;
cout <<  "  "  << endl;}

cout << " The final labels of 5th sheet:" <<  endl;

for(n2=0;n2<ny;n2++){
for(n3=0;n3<nz;n3++)
cout << ncolor[4][n2][n3] << "  " ;
cout <<  "  "  << endl;}*/




////////////////////////////////////////////////////Cluster distribution///////////////////////////////////



// The relevant labels start from 2

// The total number of independent crystalline domains is given by labelcount

int  nemcount=0;

int nclust=0,imax1=1;

for(i=1;i<=label;i++){
//cout << "i=" << i<<  " LL[i]=" << LL[i] << endl;
if(LL[i]>0) { ++labelcount;
              ni=LL[i]; 
             if(ni> imax1) imax1=ni;
             clusterfile1 << "i=" << i<<  " LL[i]=" << LL[i] << endl;
             nemcount+=LL[i];
 }
}

cout << " The total number of independent crystalline domains is:  " << labelcount << endl;

cout << " The total number of  crystalline grid elements is:  " << nemcount << endl;

cout<< "final gridtest =" << gridtest << endl;

cout << "imax1=" << imax1 << endl;

int jmax=imax1+1000;
double Vtot=Lx*Ly*Lz;
int ncluster[jmax]; int nmax1=0; // number of clusters with size i

for (i=0;i< jmax;i++) ncluster[i]=0;



for(i=1;i<=label;i++){
//cout << "i=" << i<<  " LL[i]=" << LL[i] << endl;
if(LL[i]>0) {  ni=LL[i];
if(ni>nmax1) nmax1=ni;
++ncluster[ni];
if(ni>1) ++nclust;
            
 }
}

cout << "nmax1=" << nmax1 << endl;



//cout << "up to hear fine" << endl;

 int bincluster[jmax];
double pdf1[jmax], pdf2[jmax];
for (i=0;i< jmax;i++) bincluster[i]=0;

for(i=1;i<=5;i++) bincluster[i]=ncluster[i];


bincluster[7]=ncluster[6]+ncluster[7];

bincluster[10]=ncluster[8]+ncluster[9]+ncluster[10];

for(i=20;i<=70;i+=10) {
for(j=i-9;j<=i;j++)
if(j<=imax1) bincluster[i]+=ncluster[j];}

for(j=71;j<=100;j++)
bincluster[100]+=ncluster[j];

for(i=200;i<=700;i+=100) {
for(j=i-99;j<=i;j++)
if(j<=imax1) bincluster[i]+=ncluster[j];}

for(j=701;j<=1000;j++)
if(j<=imax1) bincluster[1000]+=ncluster[j];

for(i=2000;i< jmax;i+=1000) {
for(j=i-999;j<=i;j++)
if(j<=imax1) bincluster[i]+=ncluster[j];}

n1=0;
bincluster[0]=1;
double avencluster=0, avevolume=0 ;

cout << "vol="<<  vol << endl;

for (i=1;i<jmax;i++) { if (ncluster[i]>0 ) {
avencluster+=ncluster[i]*i;
avevolume+=ncluster[i]*i*i;}

}

cout << " The average cluster size of crystalline domains is :  " << avencluster/labelcount << endl;
cout << " The average volume of crystalline domains is :  " << avevolume*vol/nemcount << endl;

avevolume=avevolume*vol/nemcount;

//avencluster=0; avevolume=0 ;
for (i=1;i<jmax;i++) { if (bincluster[i]>0 ) { 

n2=i;  
pdf1[i]=bincluster[i]*1.0/(n2-n1) ;
//avencluster+=bincluster[i]*(n1+n2)*0.5;
pdf2[i]=((n1+n2)*pdf1[i]*0.5)/nemcount;
//avevolume+=pdf2[i]*(n1+n2)*0.5*vol;
n1=n2;
}}





clusterfile <<  "clustersize volume  vol/Vtot" << "Ncluster[i] " << " pdf1 "<< " volume-pdf"<<  endl;
/*for (i=1;i<=imax1;i++) 
if (ncluster[i]>0) clusterfile <<  i  <<  " " <<   ncluster[i] << " "<< (ncluster[i]*1.0)/labelcount <<  " " << (i*ncluster[i]*1.0)/nemcount <<  endl;*/

//if ( bincluster[i]>0)
for (i=0;i<jmax;i++) { if (bincluster[i]>0) clusterfile <<  i  <<  " " << i*vol <<" " <<  i*vol/Vtot <<" "<<  bincluster[i] << " " <<  pdf1[i]/labelcount << " " << pdf2[i]  << " " << pdf2[i]*crystallinity << endl;}


//clusterfile <<  "the largest  crystalline cluster has " << imax1 << " nematic grid elements" <<  endl;
clusterfile <<  "the total number of clusters with at least 2 elements is " << nclust << " nematic grid elements" <<  endl;


clusterfile << " The total number of independent crystalline domains is:  " << labelcount << endl;
clusterfile << " The average cluster size of crystalline domains is :  " << avencluster/labelcount << endl;
clusterfile << " The average volume fraction of crystalline domains is :  " << avevolume/Vtot << endl;


clusterfile << " The total number of  crystalline grid elements is:  " << nemcount << endl;


clusterfile << " S_Ree globalNem Vtotal  Ngrid  crystallinity  avenclust avecclust/Ngrid avevol  avevol/Vtot " <<  endl;


clusterfile << S_Ree << " " << globalNem << " "<< Vtot << " " << sumcount1 << " " << crystallinity << " "<< avencluster/labelcount << " " << avencluster/labelcount/sumcount1  << " "<<   avevolume << " "<<  avevolume/Vtot <<  endl;


////////////////////////////////////////////////Making the vmd script for visualization of clusters + labeling particles based on the cluster they belong to. ///////////////////////////			    

int * clusterID,  * clusterIDbond;

clusterID=new int [Nmon+1]; // each monomer in the crystalline domains belongs to a a cluster. If it is in amorphous region, its clusterID=0
clusterIDbond= new int [Nbond+1]; // each bond in the crystalline domains belongs to  a cluster. If it is in amorphous region, clusterIDbond=0
 

for(i=1;i<=Nmon;i++){
clusterID[i]=0;
}

for(i=1;i<=Nbond;i++){
clusterIDbond[i]=0;
}


for(n1=0;n1<nx;n1++)
for(n2=0;n2<ny;n2++)
for(n3=0;n3<nz;n3++)
if(ncolor[n1][n2][n3]!=0) {
i=ncolor[n1][n2][n3];
while (LL[i]<0) i=-LL[i];
//if( LL[i]>1){ // we only count nematic grid elements belonging to a cluster size of at least 2.

vmdfile << "draw color "<<  i % 1024 << endl;
//vmdfile << "draw color "<<  "blue" << endl;
m=0;
for(l=1; l<=Nbond;l++){
 //cout << "l=" << l << endl;
m1=(int) xm[l]/lx;
m2=(int) ym[l]/ly;
m3=(int) zm[l]/lz;

//cout<< "m1=" << m1<< " m2= " << m2 << " m3=  "<< m3 << endl;
 if( (m1==n1) && ( m2==n2) && (m3==n3) ){

++m;

xp1[m]=xm[l]-lbond[l]/2*bond[l][0];
xp2[m]=xm[l]+lbond[l]/2*bond[l][0];

yp1[m]=ym[l]-lbond[l]/2*bond[l][1];
yp2[m]=ym[l]+lbond[l]/2*bond[l][1];

zp1[m]=zm[l]-lbond[l]/2*bond[l][2];
zp2[m]=zm[l]+lbond[l]/2*bond[l][2];

id=l/(Lchain-1)+l; //id of atom at the begining of bond vector.
clusterID[id]=i;
clusterIDbond[l]=i; // idbond=l;
} //end if

}// end for l
for(j=1;j<=mm[n1][n2][n3];j++){
//vmdfile <<  "graphics 0 line {" << xp1[j] << " "<< yp1[j]<< " " << zp1[j]<< "}   {"  << xp2[j] << " "<< yp2[j]<< " " << zp2[j] << "}  width 5 style solid" << endl;
vmdfile <<  "graphics 0 sphere {" << xp1[j] << " "<< yp1[j]<< " " << zp1[j]<< "}   radius 0.2 resolution 100" << endl;
vmdfile <<  "graphics 0 sphere {" << xp2[j] << " "<< yp2[j]<< " " << zp2[j]<< "}   radius 0.2 resolution 100" << endl;
vmdfile <<  "graphics 0 cylinder {" << xp1[j] << " "<< yp1[j]<< " " << zp1[j]<< "}   {"  << xp2[j] << " "<< yp2[j]<< " " << zp2[j] << "} radius 0.1 resolution 100" << endl;
} // end for j

if(LL[i]==imax1 ){
vmdfile1 << "draw color "<<  " red" << endl;
for(j=1;j<=mm[n1][n2][n3];j++){
//vmdfile <<  "graphics 0 line {" << xp1[j] << " "<< yp1[j]<< " " << zp1[j]<< "}   {"  << xp2[j] << " "<< yp2[j]<< " " << zp2[j] << "}  width 5 style solid" << endl;
vmdfile1 <<  "graphics 0 sphere {" << xp1[j] << " "<< yp1[j]<< " " << zp1[j]<< "}   radius 0.2 resolution 100" << endl;
vmdfile1 <<  "graphics 0 sphere {" << xp2[j] << " "<< yp2[j]<< " " << zp2[j]<< "}   radius 0.2 resolution 100" << endl;
vmdfile1 <<  "graphics 0 cylinder {" << xp1[j] << " "<< yp1[j]<< " " << zp1[j]<< "}   {"  << xp2[j] << " "<< yp2[j]<< " " << zp2[j] << "} radius 0.1 resolution 100" << endl;
} // end for j
}//end if LL[i]=max


// } //end if LL[i]>1


} //end if ncolor

/*
for(i=89000;i<=89020;i++)
cout << i <<" "<< clusterIDbond[i]<< endl;

*/





//////////////////////////////////// Identifying tie chains,loops & tails and  their length distributions  //////////////////////////////////////////////////////////////


//********************  initializaions********************************** id=l/(Lchain-1)+l; //id of atom at the begining of bond vector.***************/

long int **Id1_tie, **Ltie; //Id1_tie[Nchain][max_Nseg]; Id1_tie[i][j] tells us the ID of first atom of jth tie segment of chain i. Ltie[i][j] tells us the number of bonds of jth tie segment of chain i
Id1_tie= new long int *[Nchain];
for( i = 0; i < Nchain; ++i){
Id1_tie[i] = new long int [max_Nseg];
if (Id1_tie[i] == NULL) {cerr << " Id1_tie[i] Allocation problem!"; exit(1);}}

Ltie= new long int *[Nchain];
for( i = 0; i < Nchain; ++i){
Ltie[i] = new long int [max_Nseg];
if (Ltie[i] == NULL) {cerr << " Ltie[i] Allocation problem!"; exit(1);}}



for(i=0; i<Nchain;i++)
   for(n=0; n < max_Nseg;n++)
    Ltie[i][n]=0;
  
for(i=0; i<Nchain;i++)
   for(n=0; n < max_Nseg;n++)
    Id1_tie[i][n]=0; 


long int **Id1_loop, **Lloop; //Id1_loop[Nchain][max_Nseg]; Id1_loop[i][j] tells us the ID of first atom of jth loop segment of chain i. Lloop[i][j] tells us the number of bonds of jth loop segment of chain i
Id1_loop= new long int *[Nchain];
for(int i = 0; i < Nchain; ++i){
Id1_loop[i] = new long int [max_Nseg];
if (Id1_loop[i] == NULL) {cerr << " Id1_loop[i] Allocation problem!"; exit(1);}}

Lloop= new long int *[Nchain];
for(int i = 0; i < Nchain; ++i){
Lloop[i] = new long int [max_Nseg];
if (Lloop[i] == NULL) {cerr << " Lloop[i] Allocation problem!"; exit(1);}}



for(i=0; i<Nchain;i++)
   for(n=0; n < max_Nseg;n++)
    Lloop[i][n]=0;
  
for(i=0; i<Nchain;i++)
   for(n=0; n < max_Nseg;n++)
    Id1_loop[i][n]=0; 


  



int *loop_label, *tie_label, *tail_label ;

loop_label= new int [Nmon+1];
tie_label= new  int [Nmon+1];
tail_label= new int [Nmon+1];

long int loopcount[segmax+1], tiecount[segmax+1],tailcount[segmax+1];

for(i=0; i<Nmon;i++){
loop_label[i]=0;
tie_label[i]=0;
tail_label[i]=0;
}

for(i=0; i<= segmax;i++){
loopcount[i]=0;
tiecount[i]=0;
tailcount[i]=0;
}

  int idend, tie_max=2, loop_max=3, tail_max=1, L_tail;

long int flip[max_Nseg]; // array of ID of crystalline bonds after which bonds  flip  from -1 (crys) to amorphous (1) bonds.
long int Ntieseg, Nloopseg;

long int Nloop=0; // number of bonds belonging to loop segments

long int Ntail=0; // number of bonds belonging to tail segments

Ntie=0; 

long int Ntieseg_max=1;

cout << "initialization of tie loop variables went fine" << endl;

//********************  Characterizations of tie and chain segments **********************************/


 for(i=0; i<Nchain;i++){
  j=0;
 for(m=0;m<=max_Nseg;m++) flip[m]=0;
   /*** fidning flip points ****/
 
     for(n=1; n < Lchain;n++){
         id=i*(Lchain-1)+n; //id of nth bond of chain i
         idnext=id+1;

   
      if(bondamorph[id]==-1 && bondamorph[idnext]==1){
             ++j;
            flip[j]=id;

          }//end if


   if(bondamorph[id]==1 && bondamorph[idnext]==-1){
     ++j;
    flip[j]=idnext;

   }//end if

   
    } //end for n
 
/*** checking if chain ends are tail or not! ****/


if(j>0){ //if we have at least one flip point between ordered and amorphous segments.

id1=i*(Lchain-1)+1; // id of first bond of i_th chain.
      
 if(bondamorph[id1]==1){
    tail_label[id1]=1;
    L_tail=1;
    ++id1;
  while(bondamorph[id1]==1){
    tail_label[id1]=1;
    ++L_tail;
    ++id1; 
    }//end while
 ++tailcount[L_tail];
 if(L_tail > tail_max) tail_max=L_tail;
}//end if id1 tail

 idend=(i+1)*(Lchain-1); // id of last bond of i_th chain.
      
 if(bondamorph[idend]==1){
            tail_label[idend]=1;
             L_tail=1;
            --idend;
   while(bondamorph[idend]==1){
                              tail_label[idend]=1;
                              --idend; 
                             ++L_tail;
       }//end while

   ++tailcount[L_tail];
 if(L_tail > tail_max) tail_max=L_tail;

 }//end if idend tail

} // end if (j>0) i.e. if we have at least one flip point between ordered and amorphous segments.

/*** identifying tie and loops sequeces ****/

if(j>1){
   Ntieseg=0;
   Nloopseg=0;
   
   for(m=1; m < j;m++){
   id=flip[m];
   idnext=flip[m+1];
   l=idnext-id-1;
 
   if(clusterIDbond[id]!=clusterIDbond[idnext] &&  bondamorph[id+1]==1 &&  bondamorph[id-1]==-1 && bondamorph[idnext+1]==-1 && l>2 ) {
    {++Ntieseg;
   Id1_tie[i][Ntieseg]=id+1;
   Ltie[i][Ntieseg]=l;
   ++tiecount[l];
   for(k=id+1;k< idnext;k++)    tie_label[k]=1;
   if(l > tie_max) tie_max=l;
   } // end of if condition of being a tie chain!
}
  
 else { 
// if(clusterIDbond[id]==clusterIDbond[idnext] && clusterIDbond[id]!=0 && clusterIDbond[idnext]!=0){
if( l > 2 && bondamorph[id+1]==1 && bondamorph[id-1]==-1 && bondamorph[idnext+1]==-1 ){  
   ++Nloopseg;
   ++loopcount[l];
  for(k=id+1;k< idnext;k++)    loop_label[k]=1; 
 if(l > loop_max) loop_max=l;
   Id1_loop[i][Nloopseg]=id+1;
   Lloop[i][Nloopseg]=l;
     
   } //if condition of being a loop !

  } //end else


 }//end for m



} // end if j>1

if(Ntieseg  > Ntieseg_max) Ntieseg_max=Ntieseg;

 
  }// end for i 




cout << "loop_max=" << loop_max << " tie_max=" << tie_max << " tail_max=" << tail_max << " Ntieseg_max= "  << Ntieseg_max << endl;  

/*
for(i=0; i< Nchain; i++){

for(m=1;m <= max_Nseg;m++) 
           if( Lloop[i][m]==loop_max){              cout << "chain number=" << i << " m=" << m  << "  Lloop[i][m]="<<  Lloop[i][m] << endl;                                                                                                    } //end if
  }//end for i 
                                                
i=236;  for(m=1;m <= max_Nseg;m++)  cout << " Lloop[i][m]= "<< Lloop[i][m] << endl;
          
 for(n=1; n< Lchain; n++) {
                                                                           id= i*(Lchain-1)+n;
                                                                         cout << id <<" "<< bondamorph[id] << " "<< tail_label[id] << " " <<loop_label[id] << " " << tie_label[id] << " " << clusterIDbond[id]<< endl ;
                                                                           } */
//*******************************************checking which fraction of chain ends are in the amorphous phase*******************************************
int Ntailend=0;

double Nchainend=Nchain*2.0, fchainend;

for(i=0; i<Nchain;i++){
 id=i*(Lchain-1)+1; //id of first bond of chain i
 if(tail_label[id]==1) ++Ntailend;
 idend=id+Lchain-2;
 if( tail_label[idend]==1) ++Ntailend;
}

fchainend=Ntailend/Nchainend;


/***************************pdf of ties, loops and tails ***************************************/

double  pdf_tail[segmax+1], pdf_tie[segmax+1],   pdf_loop[segmax+1];
int  N_tiechain=0, N_loopchain=0, N_tailchain=0;
int ave_tie=0, ave_tail=0, ave_loop=0;
for(i=0;i<= segmax; i++){
    pdf_tie[i]=0;
   pdf_tail[i]=0;
   pdf_loop[i]=0;}



for(i=1;i<=Nbond;i++) {

if( tie_label[i]==1) ++Ntie;
if( loop_label[i]==1) ++Nloop;
if( tail_label[i]==1) ++Ntail;
 //if( tie_label[i]==1 && bondamorph[i]==-1) cout << "s.th wrong with labeling of ties !! for i=" << i << endl;
// if( tail_label[i]==1 && bondamorph[i]==-1 ) cout << "s.th wrong with labeling of tails !! for i=" << i << endl;
//if( loop_label[i]==1 && bondamorph[i]==-1) cout << "s.th wrong with labeling of loops!! for i=" << i << endl;
}
 
double f_tie=Ntie*1.0/Nbond;
double f_loop= Nloop*1.0/Nbond;
double f_tail= Ntail*1.0/Nbond;




strcpy(pdfseg,filename);
strcat(pdfseg, "_pdftieloop.txt");
ofstream tiefile;
tiefile.open(pdfseg);


tiefile <<  "m  tie_count pf_tie pdf_tie   loop_count pf_loop pdf_loop  tail_count pf_tail pdf_tail   " <<   endl;
int max_tiecount=1;

for(i=1;i<= segmax; i++){
ave_tie+=i*tiecount[i]; 
ave_loop+=i*loopcount[i]; 
ave_tail+=i*tailcount[i]; 
if(tiecount[i]> max_tiecount) max_tiecount=tiecount[i];
N_tiechain+=tiecount[i]; // total number of tie segments
N_loopchain+=loopcount[i]; // total number of loop segments
N_tailchain+=tailcount[i]; // total number of tail segments
 pdf_tie[i]=i*tiecount[i]*1.0/Ntie; // fraction of tie bonds belonging to a tie chain of length i.
   pdf_loop[i]=i*loopcount[i]*1.0/Nloop; // fraction of loop bonds belonging to a loop chain of length i.
   pdf_tail[i]=i*tailcount[i]*1.0/Ntail; // fraction of tail bonds belonging to a tie chain of length i.
}



cout << "N_tiechain=" << N_tiechain << " max_tiecount= " << max_tiecount  << endl;



/*** Calculating average of loops, tails and ties +fraction of chain ends being amorphous or crystalline ****/


float Ave_tie=0, Ave_tail=0, Ave_loop=0;
Ave_tie= ave_tie*1.0/N_tiechain;
 Ave_loop=ave_loop*1.0/N_loopchain;
 Ave_tail= ave_tail*1.0/N_tailchain;
for(i=1;i<= segmax; i++){
 tiefile << i << " "<< tiecount[i] <<" " <<   tiecount[i]*1.0/N_tiechain << " "<<  pdf_tie[i] << " "<<  loopcount[i] <<" " <<   loopcount[i]*1.0/N_loopchain << " "<<  pdf_loop[i] << " "<<  tailcount[i] << " "<<  tailcount[i]*1.0/N_tailchain 
<< " "<<  pdf_tail[i] << endl;
}   

cout  << " f_tail= "<< f_tail << " f_loop=" << f_loop << " f_tie= "<< f_tie << " f_tie+f_loop+f_tail= "<< f_tie+f_loop+f_tail <<  " f_amorph= " << Namorph*1.0/Nbond  <<" Ave_tie=" <<  Ave_tie <<" Ave_loop=" << Ave_loop
  <<" Ave_tail="  << Ave_tail << endl; 

tiefile  << " f_tail= "<< f_tail << " f_loop=" << f_loop << " f_tie= "<< f_tie << " f_tie+f_loop+f_tail= "<< f_tie+f_loop+f_tail <<  " f_amorph= " << Namorph*1.0/Nbond  <<" Ave_tie=" <<  Ave_tie <<" Ave_loop=" << Ave_loop
  <<" Ave_tail="  << Ave_tail << endl; 

segfile1 << "f_tail f_loop f_tie f_tie+f_loop+f_tail f_tie+f_loop+f_tail f_amorph fchainend Ave_tie Ave_loop Ave_tail N_tiechain N_loopchain N_tailchain "  <<  endl;

segfile1 << f_tail << "  " << f_loop << "  "<< f_tie << "  "<< f_tie+f_loop+f_tail <<  "  " << Namorph*1.0/Nbond  <<"  " << fchainend << "  "<<  Ave_tie <<"  " << Ave_loop
  <<"  "  << Ave_tail<<  "  " << N_tiechain << "  " << N_loopchain << "  " << N_tailchain  << endl;





/************************** calculating the end to end distance of  the tie chains  with different Lties ***************************************/
// you need to print id1 and idend for analysis of deformed sampels.


m=max_tiecount+1;
int **id1tieL=Allocate_2D_Integer_Array(segmax+1, m); // id1tieL[i][j] id of  first monomer of the jth tie chain of length i 
int **id2tieL=Allocate_2D_Integer_Array(segmax+1, m); //id2tieL[i][j] id of  last monomer of the jth tie chain of length i 

//int ntie[120];


 

double Rex2[segmax+1], Rey2[segmax+1], Rez2[segmax+1], Ree2[segmax+1];

for(k=0;k<=segmax;k++){
 Rex2[k]=0;
 Rey2[k]=0;
 Rez2[k]=0;
 Ree2[k]=0;
}



for(k=3;k<=segmax ;k++){
n=0;
if(tiecount[k]!=0){
for(i=0; i<Nchain;i++) 
for(j=1;j<=Ntieseg_max;j++){
if(Ltie[i][j]==k) {
 ++n;
id1tieL[k][n]= Id1_tie[i][j]/(Lchain-1)+Id1_tie[i][j];
 id2tieL[k][n]=id1tieL[k][n]+k; //id of the end monomer of tie chain of length 10
id1=id1tieL[k][n];
id2=id2tieL[k][n];
if(amorphLabel[id1]==0) cout << " s.th wrong with your algorithm " << endl;
if(amorphLabel[id2]==0) cout << " s.th wrong with your algorithm " << endl;

Rex2[k]+=(xu[id1]-xu[id2])*(xu[id1]-xu[id2]);
Rey2[k]+=(yu[id1]-yu[id2])*(yu[id1]-yu[id2]);
Rez2[k]+= (zu[id1]-zu[id2])*(zu[id1]-zu[id2]);
  } // end if
} //end for j

 Rex2[k]=Rex2[k]/tiecount[k];
 Rey2[k]=Rey2[k]/tiecount[k];
 Rez2[k]=Rez2[k]/tiecount[k];
 Ree2[k]=Rex2[k]+Rey2[k]+Rez2[k];

if(k==20) {cout << "tiecount[k]= " << tiecount[k] << " n=" << n << endl;

    cout << "k=" << k << " Rex2[k]= " << Rex2[k] << " Rey2[k]= " << Rey2[k] << " Rez2[k]= " << Rez2[k] << "Ree2[k]= " << Ree2[k] <<  endl;}

} // end if tiecount[k]!=0 


} //end for k

strcpy(pdfseg,filename);
strcat(pdfseg, "_Re2_tie.txt");
ofstream Refile;
Refile.open(pdfseg);

Refile <<  "k tiecount  Rex2  Rey2  Rez2  Ree2 "  <<  endl;

for(k=3;k<=segmax ;k++) 
   if(tiecount[k]!=0) Refile <<  k <<  " " << tiecount[k] << " " << Rex2[k] << " " << Rey2[k] << " " << Rez2[k] << " " << Ree2[k] <<  endl;
/////////////////// recording idfirst and idend of tie chains of different lengths ////////////////////////////////

strcpy(pdfseg,filename);
strcat(pdfseg, "_idtie_xyz.txt");
ofstream IDtiefile;
IDtiefile.open(pdfseg);

 for(k=3;k<=segmax ;k++){
   if(tiecount[k]!=0 & k<=segmax){
IDtiefile << k << " " << tiecount[k] << endl;
      for(n=1;n<=tiecount[k];n++){
IDtiefile << id1tieL[k][n]<< " " << id2tieL[k][n] << " " ;
if(n==tiecount[k]) IDtiefile << endl;
 } }// end for n


 }//end for k

/////////////////////////////////////////


 delete lbond, xm,ym, zm;



for(i = 0; i < Nbond; ++i)
delete [] bond[i];
delete [] bond;


release_3D_Double_Array(nn1,nx+1, ny+1, nz+1);
release_3D_Double_Array(nn2,nx+1, ny+1, nz+1);
release_3D_Double_Array(nn3,nx+1, ny+1, nz+1);


delete bondamorph;

m=max_tiecount+1;
release_2D_Integer_Array(id1tieL, segmax+1, m);
release_2D_Integer_Array(id2tieL, segmax+1, m);

delete xu;
delete yu;
delete zu;
delete amorphLabel;

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

 
