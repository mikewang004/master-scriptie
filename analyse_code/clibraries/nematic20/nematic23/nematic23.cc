// g++ -I ./ -O2 -Wno-deprecated nematic23.cc -o nematic23
// ./nematic23 equil_t_088_tdot_e-3_time_24000000.txt


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <complex>
#include <algorithm>
#include <cstring>
#include <climits>
#include <vector>
using namespace std;

/*********************************************************************************************************************/
#define  onethird   1.0/3.0
#define  pi         acos(-1.)
#define  pi2        2.0 * pi
#define  asdf       sqrt(3.0)
#define  nw         6
#define  nmol       4500
#define  nmo2       4500
#define  nvect      1

int const Lchain = 1000;
ofstream nematicfile, clusterfile, clusterfile1;

/*********************************************************************************************************************/
// Function prototypes
/*********************************************************************************************************************/
void orderparameter(int, double[][4], double[], double[], double[], double[], double*, double[]);
void orderparameter1(int, double**, double[], double[], double[], double[], double*, double[]);

inline double wrap_to_zero(double u, double lo, double L);
inline double wrap_to_box(double u, double lo, double L);
inline double min_image(double dx, double L);

void find_clusters(
    int nx, int ny, int nz, int gridcount,
    int***    ncolor,
    double*** nn1,
    double*** nn2,
    double*** nn3,
    int&      label_out,
    int*&     LL_out,
    int&      gridtest_out
);

void compute_cluster_distribution(
    int    label,
    int*   LL,
    int    gridtest,
    double Lx, double Ly, double Lz,
    double vol,
    double S_Ree,
    double globalNem,
    int    sumcount1,
    double crystallinity
);

/*********************************************************************************************************************/
// Utility
/*********************************************************************************************************************/
int &Max(int &a, int &b) { return a > b ? a : b; }
int &Min(int &a, int &b) { return a <= b ? a : b; }

/*********************************************************************************************************************/
// 3D array allocation helpers
/*********************************************************************************************************************/
double*** Allocate_3D_Double_Array(int x, int y, int z)
{
    double*** arr = new double**[x];
    for (int i = 0; i < x; ++i) {
        arr[i] = new double*[y];
        for (int j = 0; j < y; ++j) {
            arr[i][j] = new double[z];
            for (int k = 0; k < z; ++k)
                arr[i][j][k] = 0.0;
        }
    }
    return arr;
}

void release_3D_Double_Array(double*** arr, int x, int y, int z)
{
    for (int i = 0; i < x; ++i) {
        for (int j = 0; j < y; ++j)
            delete[] arr[i][j];
        delete[] arr[i];
    }
    delete[] arr;
}

int*** Allocate_3D_Integer_Array(int x, int y, int z)
{
    int*** arr = new int**[x];
    for (int i = 0; i < x; ++i) {
        arr[i] = new int*[y];
        for (int j = 0; j < y; ++j) {
            arr[i][j] = new int[z];
            for (int k = 0; k < z; ++k)
                arr[i][j][k] = 0;
        }
    }
    return arr;
}

void release_3D_Integer_Array(int*** arr, int x, int y, int z)
{
    for (int i = 0; i < x; ++i) {
        for (int j = 0; j < y; ++j)
            delete[] arr[i][j];
        delete[] arr[i];
    }
    delete[] arr;
}

/*********************************************************************************************************************/
// 2D array allocation helpers
/*********************************************************************************************************************/
double** Allocate_2D_Double_Array(int x, int y)
{
    double** arr = new double*[x];
    for (int i = 0; i < x; ++i)
        arr[i] = new double[y];
    return arr;
}

void release_2D_Double_Array(double** arr, int x, int y)
{
    for (int i = 0; i < x; ++i)
        delete[] arr[i];
    delete[] arr;
}

int** Allocate_2D_Integer_Array(int x, int y)
{
    int** arr = new int*[x];
    for (int i = 0; i < x; ++i)
        arr[i] = new int[y];
    return arr;
}

void release_2D_Integer_Array(int** arr, int x, int y)
{
    for (int i = 0; i < x; ++i)
        delete[] arr[i];
    delete[] arr;
}



/*********************************************************************************************************************/
// main
/*********************************************************************************************************************/
int main(int argc, const char* argv[])
{

    
    ifstream xyzfile, bondfile;

    double S;
    double director[3];
    int i, id, type, nsequenc, j, l, k;
    int Natomtype, Nbondtype, Nangletype;
    int timestep;
    int Nmon, Nbond, Nangle;
    int Nchain;

    string str, str1, str2, str3, str4, str5, str6, str7, str8, str9, str10;
    double xlo, xhi, ylo, yhi, zlo, zhi, Lx, Ly, Lz;
    int index, id1, id2, molid, ix, iy, iz;
    int m, n, Natom;
    int Nmax   = 4330000;
    int Mchain = 1000;
    double xx, yy, zz;

    double xp1[5000], xp2[5000], yp1[5000], yp2[5000], zp1[5000], zp2[5000];

    double *xp, *yp, *zp, *xu, *yu, *zu;
    xp = new double[Nmax];
    yp = new double[Nmax];
    zp = new double[Nmax];
    xu = new double[Nmax];
    yu = new double[Nmax];
    zu = new double[Nmax];

    double cnd[4];

    cout << "New long array created successfully" << endl;
    cout << "argv[1] = " << argv[1] << "   its length is " << strlen(argv[1]) << endl;

    const int fnamelength = strlen(argv[1]);
    char filename[fnamelength + 1];
    for (i = 0; i <= fnamelength; i++)
        filename[i] = argv[1][i];

    // ---- Open and parse input file ----
    xyzfile.open(filename);
    if (!xyzfile) {
        cerr << "Error: file could not be opened" << endl;
        exit(1);
    }
    getline(xyzfile, str);
    cout << str << endl;

    xyzfile >> Nmon >> str1;
    cout << "Nmon: " << Nmon << str1 << endl;
    xyzfile >> Natomtype >> str1 >> str2;
    xyzfile >> Nbond >> str1;
    cout << Nbond << str1 << endl;
    xyzfile >> Nbondtype >> str1 >> str2;
    xyzfile >> Nangle >> str1;
    cout << Nangle << str1 << endl;
    xyzfile >> Nangletype >> str1 >> str2;

    xyzfile >> xlo >> xhi >> str1 >> str2;
    xyzfile >> ylo >> yhi >> str1 >> str2;
    xyzfile >> zlo >> zhi >> str1 >> str2;
    cout << "xlo=" << xlo << " xhi=" << xhi << endl;
    cout << "ylo=" << ylo << " yhi=" << yhi << endl;
    cout << "zlo=" << zlo << " zhi=" << zhi << endl;

    Lx = xhi - xlo;
    Ly = yhi - ylo;
    Lz = zhi - zlo;
    double Lxhalf = Lx / 2.0;
    double Lyhalf = Ly / 2.0;
    double Lzhalf = Lz / 2.0;



    xyzfile >> str1;

    std::string line;
    std::getline(xyzfile, line);   // discard header line

    // for (i = 1; i <= Nmon; i++) {
    //     xyzfile >> id >> molid >> xu[id] >> yu[id] >> zu[id];
    //     xp[id] = xu[id] - xlo;
    //     yp[id] = yu[id] - ylo;
    //     zp[id] = zu[id] - zlo;
    // }

// After reading xlo,xhi,ylo,yhi,zlo,zhi and computing Lx,Ly,Lz:

    for (i = 1; i <= Nmon; i++) {
        xyzfile >> id >> molid >> xu[id] >> yu[id] >> zu[id];

        // wrap xu,yu,zu into [0,L) as C uses later
        xp[id] = wrap_to_zero(xu[id], xlo, Lx);
        yp[id] = wrap_to_zero(yu[id], ylo, Ly);
        zp[id] = wrap_to_zero(zu[id], zlo, Lz);
    }
    xyzfile.close();
    cout << "reading of data file was made successfully" << endl;

    // ---- Chain / bond setup ----
    int max_m = 1;
    Nchain = Nmon / Lchain;
    double rho = Nmon / (Lx * Ly * Lz);
    cout << "Nchain=" << Nchain << " rho=" << rho << endl;

    Nbond = Nchain * (Lchain - 1);
    int nx, ny, nz, n1, n2, n3;
    double lx, ly, lz;

    double *lbond = new double[Nbond + 1];
    int m1, m2, m3;

    double lambda[4], eevp[4], eevm[4], eevo[4], lambdasara[4];
    lambda[0] = eevp[0] = eevm[0] = eevo[0] = 0;

    double *xm = new double[Nbond + 1];
    double *ym = new double[Nbond + 1];
    double *zm = new double[Nbond + 1];

    double **bond = new double*[Nbond + 1];
    for (i = 0; i < Nbond + 1; ++i) {
        bond[i] = new double[3];
        if (!bond[i]) { cerr << "bond[i] allocation failed"; exit(1); }
    }
    cout << "Long array was created successfully" << endl;

    const int MAX_IN_CELL = 1000;
    int count[MAX_IN_CELL], nematiccount[MAX_IN_CELL];
    std::vector<int> atom1_of_bond(Nbond + 1);
    std::vector<int> atom2_of_bond(Nbond + 1);

    for (n = 0; n < MAX_IN_CELL; n++) {
        count[n] = 0;
        nematiccount[n] = 0;
    }

    // ---- Build bond vectors ----
    for (int nc = 1; nc <= Nchain; ++nc) {
        for (i = 1; i < Lchain; ++i) {
            int atom1 = (nc - 1) * Lchain + i;
            int atom2 = atom1 + 1;
            l = (nc - 1) * (Lchain - 1) + i;

            double bx = xu[atom2] - xu[atom1];
            double by = yu[atom2] - yu[atom1];
            double bz = zu[atom2] - zu[atom1];

            xm[l] = (xp[atom1] + xp[atom2]) / 2.0;
            ym[l] = (yp[atom1] + yp[atom2]) / 2.0;
            zm[l] = (zp[atom1] + zp[atom2]) / 2.0;
            bx = min_image(bx, Lx);
            by = min_image(by, Ly);
            bz = min_image(bz, Lz);

            double len = sqrt(bx*bx + by*by + bz*bz);
            bond[l][0] = bx / len;
            bond[l][1] = by / len;
            bond[l][2] = bz / len;

            atom1_of_bond[l] = atom1;
            atom2_of_bond[l] = atom2;
        }
    }
    cout << "l=" << l << endl;

    double Dhalf2 = 0.25 * (Lx*Lx + Ly*Ly + Lz*Lz);

    char bondvec[200];
    strcpy(bondvec, filename);
    strcat(bondvec, "_bondvecs.txt");

    // ---- Compute grid dimensions BEFORE allocating 3D arrays ----
    double Lgrid = 2.0;
    nx = (int)(Lx / Lgrid);
    ny = (int)(Ly / Lgrid);
    nz = (int)(Lz / Lgrid);
    cout << "nx=" << nx << " ny=" << ny << " nz=" << nz << endl;

    lx = Lx / (nx * 1.0);
    ly = Ly / (ny * 1.0);
    lz = Lz / (nz * 1.0);
    double vol   = lx * ly * lz;
    int    Ngrid = nx * ny * nz;
    cout << "lx=" << lx << " ly=" << ly << " lz=" << lz
         << "   Ngrids=" << Ngrid << endl;

    // ---- Allocate 3D arrays AFTER nx,ny,nz are known ----
    double*** SS     = Allocate_3D_Double_Array(nx, ny, nz);
    double*** nn1arr = Allocate_3D_Double_Array(nx, ny, nz);
    double*** nn2arr = Allocate_3D_Double_Array(nx, ny, nz);
    double*** nn3arr = Allocate_3D_Double_Array(nx, ny, nz);
    int***    mm     = Allocate_3D_Integer_Array(nx, ny, nz);
    int***    ncolor = Allocate_3D_Integer_Array(nx, ny, nz);

    double u[MAX_IN_CELL][4];

    ofstream bondvecfile1;
    bondvecfile1.open(bondvec);
    bondvecfile1 << "atom1 bx by bz xm ym zm nx ny nz" << endl;

    int idcell[10000];
    int *amorphLabel = new int[Nmon + 1];
    cout << " up to here fine" << endl;

    // ---- Output files ----
    char nemname[200];
    strcpy(nemname, filename);
    strcat(nemname, "_nematic.dat");
    nematicfile.open(nemname);
    nematicfile << "S_glob SS  lambda_1_glob  lambda_2_glob  lambda_3_glob "
                   "director_1_glob  director_2_glob  director_3_glob" << endl;

    char clustername[200];
    strcpy(clustername, filename);
    strcat(clustername, "_cluster.dat");
    clusterfile.open(clustername);

    char clustername1[200];
    strcpy(clustername1, filename);
    strcat(clustername1, "_clust.dat");
    clusterfile1.open(clustername1);

    int gridcount = 0;
    int cc = 0, ndomains = 0;

    // ---- Global nematic order parameter ----
    orderparameter1(Nbond, bond, lambda, eevp, eevm, eevo, &S, director);
    double globalNem = S;
    nematicfile << "The global nematic order parameter is: " << globalNem << endl;

    double bintheta = 0.025, binS = 0.02;
    int Namorph = 0, tt, nS, Nss;
    tt  = (int)(pi / 2.0 / bintheta) + 1;
    cout << "tt=" << tt << endl;

    double nem = 0;
    int Ns0 = (int)(1.0 / binS);
    Nss = Ns0 * 3 / 2 + 1;
    cout << "Nss=" << Nss << " Ns0=" << Ns0 << endl;

    int *countS     = new int[Nss + 2]();
    int *counttheta1 = new int[tt + 2]();
    int *counttheta2 = new int[tt + 2]();
    int *counttheta3 = new int[tt + 2]();

    int *bondamorph = new int[Nbond + 1]();
    int *idbond     = new int[MAX_IN_CELL]();

    // ---- Grid assignment loop ----
    for (n1 = 0; n1 < nx; n1++)
    for (n2 = 0; n2 < ny; n2++)
    for (n3 = 0; n3 < nz; n3++) {
        m = 0;
        for (l = 1; l <= Nbond; l++) {
            m1 = (int)(xm[l] / lx);
            m2 = (int)(ym[l] / ly);
            m3 = (int)(zm[l] / lz);

            if (m1 == n1 && m2 == n2 && m3 == n3) {
                if (m < MAX_IN_CELL - 1) {
                    ++m;
                    u[m][1] = bond[l][0];
                    u[m][2] = bond[l][1];
                    u[m][3] = bond[l][2];
                    idcell[m] = l / (Lchain - 1) + l;
                    idbond[m] = l;
                }
            }
        }

        if (m > max_m) max_m = m;

        int m_idx = (m < MAX_IN_CELL) ? m : (MAX_IN_CELL - 1);
        ++count[m_idx];

        if (m > 1) {
            nematicfile << "m=" << m << endl;
            nematicfile << "n1 n2 n3" << endl;
            nematicfile << n1 << " " << n2 << " " << n3 << endl;
            orderparameter(m, u, lambda, eevp, eevm, eevo, &S, director);
            mm[n1][n2][n3]     = m;
            SS[n1][n2][n3]     = lambda[1];
            // double S1 = fabs(lambda[1]);
            // double S2 = fabs(lambda[2]);
            // double S3 = fabs(lambda[3]);
            // double S0 = S1;
            // if (S2 > S0) S0 = S2;
            // if (S3 > S0) S0 = S3;

            // SS[n1][n2][n3] = S0;
            nn1arr[n1][n2][n3] = eevp[1];
            nn2arr[n1][n2][n3] = eevp[2];
            nn3arr[n1][n2][n3] = eevp[3];

            nem = S / binS;
            nS  = (int)nem;
            if (nS > 0) {
                if (nS > Nss) nS = Nss;
                countS[nS]++;
            } else {
                nS = abs(nS) + Ns0;
                if (nS > Nss) nS = Nss;
                countS[nS]++;
            }
        }

        if (SS[n1][n2][n3] > 0.8) {
            ncolor[n1][n2][n3] = 1;
            ++nematiccount[m_idx];
            ++gridcount;
        }
    }

    // ---- Cluster labeling ----
    int  label_out   = 0;
    int* LL_out      = nullptr;
    int  gridtest_out = 0;

    find_clusters(nx, ny, nz, gridcount,
                  ncolor, nn1arr, nn2arr, nn3arr,
                  label_out, LL_out, gridtest_out);

    // ---- Cluster size distribution ----
    // Set S_Ree, sumcount1, crystallinity to your actual values here
    double S_Ree       = 0.0;   // replace with your computed value
    int    sumcount1   = Ngrid; // or your actual occupied grid count
    double crystallinity = (gridcount > 0 && Ngrid > 0)
                           ? (double)gridcount / Ngrid
                           : 0.0;

    compute_cluster_distribution(
        label_out, LL_out, gridtest_out,
        Lx, Ly, Lz, vol,
        S_Ree, globalNem, sumcount1, crystallinity
    );

    delete[] LL_out;

    // ---- Cleanup ----
    release_3D_Double_Array(SS,     nx, ny, nz);
    release_3D_Double_Array(nn1arr, nx, ny, nz);
    release_3D_Double_Array(nn2arr, nx, ny, nz);
    release_3D_Double_Array(nn3arr, nx, ny, nz);
    release_3D_Integer_Array(mm,    nx, ny, nz);
    release_3D_Integer_Array(ncolor,nx, ny, nz);

    delete[] xp; delete[] yp; delete[] zp;
    delete[] xu; delete[] yu; delete[] zu;
    delete[] xm; delete[] ym; delete[] zm;
    delete[] lbond;
    delete[] amorphLabel;
    delete[] bondamorph;
    delete[] idbond;
    delete[] countS;
    delete[] counttheta1; delete[] counttheta2; delete[] counttheta3;
    for (i = 0; i < Nbond + 1; ++i) delete[] bond[i];
    delete[] bond;

    return 0;
}

/*********************************************************************************************************************/
// Cluster labeling helpers
/*********************************************************************************************************************/
#define FIND_ROOT_CLUSTER(x) while (LL[(x)] < 0) (x) = -LL[(x)]

inline void merge_clusters_int(int i, int j, int* LL)
{
    if (i == j) return;
    if (i < j) { LL[i] += LL[j]; LL[j] = -i; }
    else        { LL[j] += LL[i]; LL[i] = -j; }
}

inline double compute_P2Q_cluster(double a1, double a2, double a3,
                                   double b1, double b2, double b3)
{
    double cosQ = a1*b1 + a2*b2 + a3*b3;
    return 1.5 * cosQ * cosQ - 0.5;
}

/*********************************************************************************************************************/
// find_clusters
// Outputs: label_out  — highest label used
//          LL_out     — heap-allocated LL array (caller must delete[])
//          gridtest_out — count of occupied sites visited
/*********************************************************************************************************************/
void find_clusters(
    int nx, int ny, int nz, int gridcount,
    int***    ncolor,
    double*** nn1,
    double*** nn2,
    double*** nn3,
    int&      label_out,
    int*&     LL_out,
    int&      gridtest_out)
{
    const double SC = 0.97;

    int* LL = new int[gridcount + 2]();
    int  label    = 1;
    int  gridtest = 0;
    int  n1, n2, n3;

    // ----------------------------------------------------------
    // 1D scan: first row [0][0][n3]
    // ----------------------------------------------------------
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
                ++label; ncolor[0][0][n3] = label; ++LL[label];
            } else {
                int i = ncolor[0][0][n3-1]; FIND_ROOT_CLUSTER(i);
                ncolor[0][0][n3] = i; ++LL[i];
            }
        }
    }

    // PBC n3 for [0][0][*]
    if (ncolor[0][0][0] != 0 && ncolor[0][0][nz-1] != 0 &&
        ncolor[0][0][0] != ncolor[0][0][nz-1]) {
        double P2Q = compute_P2Q_cluster(
            nn1[0][0][0], nn2[0][0][0], nn3[0][0][0],
            nn1[0][0][nz-1], nn2[0][0][nz-1], nn3[0][0][nz-1]);
        if (P2Q >= SC) {
            int i = ncolor[0][0][0], j = ncolor[0][0][nz-1];
            FIND_ROOT_CLUSTER(i); FIND_ROOT_CLUSTER(j);
            if (i != j) {
                merge_clusters_int(i, j, LL);
                int root = min(i, j);
                ncolor[0][0][0] = root; ncolor[0][0][nz-1] = root;
            }
        }
    }

    // ----------------------------------------------------------
    // 2D scan: first sheet [0][n2][n3]
    // ----------------------------------------------------------
    for (n2 = 1; n2 < ny; n2++) {

        // leftmost [0][n2][0]
        if (ncolor[0][n2][0] != 0) {
            ++gridtest;
            bool has_n2 = (ncolor[0][n2-1][0] != 0);
            bool conn_n2 = false;
            if (has_n2) {
                double P2Q1 = compute_P2Q_cluster(
                    nn1[0][n2][0], nn2[0][n2][0], nn3[0][n2][0],
                    nn1[0][n2-1][0], nn2[0][n2-1][0], nn3[0][n2-1][0]);
                conn_n2 = (P2Q1 >= SC);
            }
            if (!conn_n2) {
                ++label; ncolor[0][n2][0] = label; ++LL[label];
            } else {
                int i = ncolor[0][n2-1][0]; FIND_ROOT_CLUSTER(i);
                ncolor[0][n2][0] = i; ++LL[i];
            }
        }

        // interior [0][n2][n3]
        for (n3 = 1; n3 < nz; n3++) {
            if (ncolor[0][n2][n3] != 0) {
                ++gridtest;
                bool has_n3 = (ncolor[0][n2][n3-1] != 0);
                bool has_n2 = (ncolor[0][n2-1][n3]  != 0);
                double P2Q = 0, P2Q1 = 0;
                if (has_n3) P2Q  = compute_P2Q_cluster(
                    nn1[0][n2][n3], nn2[0][n2][n3], nn3[0][n2][n3],
                    nn1[0][n2][n3-1], nn2[0][n2][n3-1], nn3[0][n2][n3-1]);
                if (has_n2) P2Q1 = compute_P2Q_cluster(
                    nn1[0][n2][n3], nn2[0][n2][n3], nn3[0][n2][n3],
                    nn1[0][n2-1][n3], nn2[0][n2-1][n3], nn3[0][n2-1][n3]);
                bool conn_n3 = has_n3 && (P2Q  >= SC);
                bool conn_n2 = has_n2 && (P2Q1 >= SC);
                if (!conn_n3 && !conn_n2) {
                    ++label; ncolor[0][n2][n3] = label; ++LL[label];
                } else {
                    int ii = -1, jj = -1;
                    if (conn_n3) { ii = ncolor[0][n2][n3-1]; FIND_ROOT_CLUSTER(ii); }
                    if (conn_n2) { jj = ncolor[0][n2-1][n3]; FIND_ROOT_CLUSTER(jj); }
                    int root = INT_MAX;
                    if (conn_n3) root = min(root, ii);
                    if (conn_n2) root = min(root, jj);
                    ncolor[0][n2][n3] = root; ++LL[root];
                    if (conn_n3) merge_clusters_int(root, ii, LL);
                    if (conn_n2) merge_clusters_int(root, jj, LL);
                }
            }
        }

        // PBC n3 for [0][n2][*]
        if (ncolor[0][n2][0] != 0 && ncolor[0][n2][nz-1] != 0 &&
            ncolor[0][n2][0] != ncolor[0][n2][nz-1]) {
            double P2Q = compute_P2Q_cluster(
                nn1[0][n2][0], nn2[0][n2][0], nn3[0][n2][0],
                nn1[0][n2][nz-1], nn2[0][n2][nz-1], nn3[0][n2][nz-1]);
            if (P2Q >= SC) {
                int i = ncolor[0][n2][0], j = ncolor[0][n2][nz-1];
                FIND_ROOT_CLUSTER(i); FIND_ROOT_CLUSTER(j);
                if (i != j) {
                    merge_clusters_int(i, j, LL);
                    int root = min(i, j);
                    ncolor[0][n2][0] = root; ncolor[0][n2][nz-1] = root;
                }
            }
        }
    }

    // PBC n2 for n1=0
    for (n3 = 0; n3 < nz; n3++) {
        if (ncolor[0][0][n3] != 0 && ncolor[0][ny-1][n3] != 0 &&
            ncolor[0][0][n3] != ncolor[0][ny-1][n3]) {
            double P2Q = compute_P2Q_cluster(
                nn1[0][0][n3], nn2[0][0][n3], nn3[0][0][n3],
                nn1[0][ny-1][n3], nn2[0][ny-1][n3], nn3[0][ny-1][n3]);
            if (P2Q >= SC) {
                int i = ncolor[0][0][n3], j = ncolor[0][ny-1][n3];
                FIND_ROOT_CLUSTER(i); FIND_ROOT_CLUSTER(j);
                if (i != j) {
                    merge_clusters_int(i, j, LL);
                    int root = min(i, j);
                    ncolor[0][0][n3] = root; ncolor[0][ny-1][n3] = root;
                }
            }
        }
    }

    // ----------------------------------------------------------
    // 3D scan: sheets n1 = 1..nx-1
    // ----------------------------------------------------------
    for (n1 = 1; n1 < nx; n1++) {

        // first site [n1][0][0]
        if (ncolor[n1][0][0] != 0) {
            ++gridtest;
            bool has_n1 = (ncolor[n1-1][0][0] != 0);
            bool conn_n1 = false;
            if (has_n1) {
                double P2Q2 = compute_P2Q_cluster(
                    nn1[n1][0][0], nn2[n1][0][0], nn3[n1][0][0],
                    nn1[n1-1][0][0], nn2[n1-1][0][0], nn3[n1-1][0][0]);
                conn_n1 = (P2Q2 >= SC);
            }
            if (!conn_n1) {
                ++label; ncolor[n1][0][0] = label; ++LL[label];
            } else {
                int i = ncolor[n1-1][0][0]; FIND_ROOT_CLUSTER(i);
                ncolor[n1][0][0] = i; ++LL[i];
            }
        }

        // first row [n1][0][n3]
        for (n3 = 1; n3 < nz; n3++) {
            if (ncolor[n1][0][n3] != 0) {
                ++gridtest;
                bool has_n3 = (ncolor[n1][0][n3-1] != 0);
                bool has_n1 = (ncolor[n1-1][0][n3] != 0);
                double P2Q = 0, P2Q2 = 0;
                if (has_n3) P2Q  = compute_P2Q_cluster(
                    nn1[n1][0][n3], nn2[n1][0][n3], nn3[n1][0][n3],
                    nn1[n1][0][n3-1], nn2[n1][0][n3-1], nn3[n1][0][n3-1]);
                if (has_n1) P2Q2 = compute_P2Q_cluster(
                    nn1[n1][0][n3], nn2[n1][0][n3], nn3[n1][0][n3],
                    nn1[n1-1][0][n3], nn2[n1-1][0][n3], nn3[n1-1][0][n3]);
                bool conn_n3 = has_n3 && (P2Q  >= SC);
                bool conn_n1 = has_n1 && (P2Q2 >= SC);
                if (!conn_n3 && !conn_n1) {
                    ++label; ncolor[n1][0][n3] = label; ++LL[label];
                } else {
                    int ii = -1, kk = -1;
                    if (conn_n3) { ii = ncolor[n1][0][n3-1]; FIND_ROOT_CLUSTER(ii); }
                    if (conn_n1) { kk = ncolor[n1-1][0][n3]; FIND_ROOT_CLUSTER(kk); }
                    int root = INT_MAX;
                    if (conn_n3) root = min(root, ii);
                    if (conn_n1) root = min(root, kk);
                    ncolor[n1][0][n3] = root; ++LL[root];
                    if (conn_n3) merge_clusters_int(root, ii, LL);
                    if (conn_n1) merge_clusters_int(root, kk, LL);
                }
            }
        }

        // PBC n3 for first row of sheet n1
        if (ncolor[n1][0][0] != 0 && ncolor[n1][0][nz-1] != 0 &&
            ncolor[n1][0][0] != ncolor[n1][0][nz-1]) {
            double P2Q = compute_P2Q_cluster(
                nn1[n1][0][0], nn2[n1][0][0], nn3[n1][0][0],
                nn1[n1][0][nz-1], nn2[n1][0][nz-1], nn3[n1][0][nz-1]);
            if (P2Q >= SC) {
                int i = ncolor[n1][0][0], j = ncolor[n1][0][nz-1];
                FIND_ROOT_CLUSTER(i); FIND_ROOT_CLUSTER(j);
                if (i != j) {
                    merge_clusters_int(i, j, LL);
                    int root = min(i, j);
                    ncolor[n1][0][0] = root; ncolor[n1][0][nz-1] = root;
                }
            }
        }

        // left column [n1][n2][0] + interior [n1][n2][n3]
        for (n2 = 1; n2 < ny; n2++) {

            // left column site [n1][n2][0]
            if (ncolor[n1][n2][0] != 0) {
                ++gridtest;
                bool has_n2 = (ncolor[n1][n2-1][0] != 0);
                bool has_n1 = (ncolor[n1-1][n2][0] != 0);
                double P2Q1 = 0, P2Q2 = 0;
                if (has_n2) P2Q1 = compute_P2Q_cluster(
                    nn1[n1][n2][0], nn2[n1][n2][0], nn3[n1][n2][0],
                    nn1[n1][n2-1][0], nn2[n1][n2-1][0], nn3[n1][n2-1][0]);
                if (has_n1) P2Q2 = compute_P2Q_cluster(
                    nn1[n1][n2][0], nn2[n1][n2][0], nn3[n1][n2][0],
                    nn1[n1-1][n2][0], nn2[n1-1][n2][0], nn3[n1-1][n2][0]);
                bool conn_n2 = has_n2 && (P2Q1 >= SC);
                bool conn_n1 = has_n1 && (P2Q2 >= SC);
                if (!conn_n2 && !conn_n1) {
                    ++label; ncolor[n1][n2][0] = label; ++LL[label];
                } else {
                    int jj = -1, kk = -1;
                    if (conn_n2) { jj = ncolor[n1][n2-1][0]; FIND_ROOT_CLUSTER(jj); }
                    if (conn_n1) { kk = ncolor[n1-1][n2][0]; FIND_ROOT_CLUSTER(kk); }
                    int root = INT_MAX;
                    if (conn_n2) root = min(root, jj);
                    if (conn_n1) root = min(root, kk);
                    ncolor[n1][n2][0] = root; ++LL[root];
                    if (conn_n2) merge_clusters_int(root, jj, LL);
                    if (conn_n1) merge_clusters_int(root, kk, LL);
                }
            }

            // full interior [n1][n2][n3]
            for (n3 = 1; n3 < nz; n3++) {
                if (ncolor[n1][n2][n3] != 0) {
                    ++gridtest;
                    bool has_n3 = (ncolor[n1][n2][n3-1] != 0);
                    bool has_n2 = (ncolor[n1][n2-1][n3]  != 0);
                    bool has_n1 = (ncolor[n1-1][n2][n3]  != 0);
                    double P2Q = 0, P2Q1 = 0, P2Q2 = 0;
                    if (has_n3) P2Q  = compute_P2Q_cluster(
                        nn1[n1][n2][n3], nn2[n1][n2][n3], nn3[n1][n2][n3],
                        nn1[n1][n2][n3-1], nn2[n1][n2][n3-1], nn3[n1][n2][n3-1]);
                    if (has_n2) P2Q1 = compute_P2Q_cluster(
                        nn1[n1][n2][n3], nn2[n1][n2][n3], nn3[n1][n2][n3],
                        nn1[n1][n2-1][n3], nn2[n1][n2-1][n3], nn3[n1][n2-1][n3]);
                    if (has_n1) P2Q2 = compute_P2Q_cluster(
                        nn1[n1][n2][n3], nn2[n1][n2][n3], nn3[n1][n2][n3],
                        nn1[n1-1][n2][n3], nn2[n1-1][n2][n3], nn3[n1-1][n2][n3]);
                    bool conn_n3 = has_n3 && (P2Q  >= SC);
                    bool conn_n2 = has_n2 && (P2Q1 >= SC);
                    bool conn_n1 = has_n1 && (P2Q2 >= SC);
                    if (!conn_n3 && !conn_n2 && !conn_n1) {
                        ++label; ncolor[n1][n2][n3] = label; ++LL[label];
                    } else {
                        int ii = -1, jj = -1, kk = -1;
                        if (conn_n3) { ii = ncolor[n1][n2][n3-1]; FIND_ROOT_CLUSTER(ii); }
                        if (conn_n2) { jj = ncolor[n1][n2-1][n3]; FIND_ROOT_CLUSTER(jj); }
                        if (conn_n1) { kk = ncolor[n1-1][n2][n3]; FIND_ROOT_CLUSTER(kk); }
                        int root = INT_MAX;
                        if (conn_n3) root = min(root, ii);
                        if (conn_n2) root = min(root, jj);
                        if (conn_n1) root = min(root, kk);
                        ncolor[n1][n2][n3] = root; ++LL[root];
                        if (conn_n3) merge_clusters_int(root, ii, LL);
                        if (conn_n2) merge_clusters_int(root, jj, LL);
                        if (conn_n1) merge_clusters_int(root, kk, LL);
                    }
                }
            }

            // PBC n3 for [n1][n2][*]
            if (ncolor[n1][n2][0] != 0 && ncolor[n1][n2][nz-1] != 0 &&
                ncolor[n1][n2][0] != ncolor[n1][n2][nz-1]) {
                double P2Q = compute_P2Q_cluster(
                    nn1[n1][n2][0], nn2[n1][n2][0], nn3[n1][n2][0],
                    nn1[n1][n2][nz-1], nn2[n1][n2][nz-1], nn3[n1][n2][nz-1]);
                if (P2Q >= SC) {
                    int i = ncolor[n1][n2][0], j = ncolor[n1][n2][nz-1];
                    FIND_ROOT_CLUSTER(i); FIND_ROOT_CLUSTER(j);
                    if (i != j) {
                        merge_clusters_int(i, j, LL);
                        int root = min(i, j);
                        ncolor[n1][n2][0] = root; ncolor[n1][n2][nz-1] = root;
                    }
                }
            }
        }

        // PBC n2 per sheet
        for (n3 = 0; n3 < nz; n3++) {
            if (ncolor[n1][0][n3] != 0 && ncolor[n1][ny-1][n3] != 0 &&
                ncolor[n1][0][n3] != ncolor[n1][ny-1][n3]) {
                double P2Q = compute_P2Q_cluster(
                    nn1[n1][0][n3], nn2[n1][0][n3], nn3[n1][0][n3],
                    nn1[n1][ny-1][n3], nn2[n1][ny-1][n3], nn3[n1][ny-1][n3]);
                if (P2Q >= SC) {
                    int i = ncolor[n1][0][n3], j = ncolor[n1][ny-1][n3];
                    FIND_ROOT_CLUSTER(i); FIND_ROOT_CLUSTER(j);
                    if (i != j) {
                        merge_clusters_int(i, j, LL);
                        int root = min(i, j);
                        ncolor[n1][0][n3] = root; ncolor[n1][ny-1][n3] = root;
                    }
                }
            }
        }
    }

    // PBC n1 direction
    for (n2 = 0; n2 < ny; n2++)
    for (n3 = 0; n3 < nz; n3++) {
        if (ncolor[0][n2][n3] != 0 && ncolor[nx-1][n2][n3] != 0 &&
            ncolor[0][n2][n3] != ncolor[nx-1][n2][n3]) {
            double P2Q = compute_P2Q_cluster(
                nn1[0][n2][n3], nn2[0][n2][n3], nn3[0][n2][n3],
                nn1[nx-1][n2][n3], nn2[nx-1][n2][n3], nn3[nx-1][n2][n3]);
            if (P2Q >= SC) {
                int i = ncolor[0][n2][n3], j = ncolor[nx-1][n2][n3];
                FIND_ROOT_CLUSTER(i); FIND_ROOT_CLUSTER(j);
                if (i != j) {
                    merge_clusters_int(i, j, LL);
                    int root = min(i, j);
                    ncolor[0][n2][n3] = root; ncolor[nx-1][n2][n3] = root;
                }
            }
        }
    }

    // Final relabeling pass
    for (n1 = 0; n1 < nx; n1++)
    for (n2 = 0; n2 < ny; n2++)
    for (n3 = 0; n3 < nz; n3++) {
        if (ncolor[n1][n2][n3] != 0) {
            int r = ncolor[n1][n2][n3];
            FIND_ROOT_CLUSTER(r);
            ncolor[n1][n2][n3] = r;
        }
    }

    label_out    = label;
    LL_out       = LL;      // caller owns this memory
    gridtest_out = gridtest;
}

/*********************************************************************************************************************/
// compute_cluster_distribution
/*********************************************************************************************************************/
void compute_cluster_distribution(
    int    label,
    int*   LL,
    int    gridtest,
    double Lx, double Ly, double Lz,
    double vol,
    double S_Ree,
    double globalNem,
    int    sumcount1,
    double crystallinity)
{
    int i, j;

    // ---- Pass 1: count labels, find max cluster size ----
    int labelcount = 0;
    int nemcount   = 0;
    int ni;
    int nclust = 0;
    int imax1  = 1;

    for (i = 1; i <= label; i++) {
        if (LL[i] > 0) {
            ++labelcount;
            ni = LL[i];
            if (ni > imax1) imax1 = ni;
            clusterfile1 << "i=" << i << " LL[i]=" << LL[i] << endl;
            nemcount += LL[i];
        }
    }

    cout << " The total number of independent crystalline domains is:  "
         << labelcount << endl;
    cout << " The total number of crystalline grid elements is:  "
         << nemcount << endl;
    cout << " final gridtest = " << gridtest << endl;
    cout << " imax1 = " << imax1 << endl;

    if (labelcount == 0 || nemcount == 0) {
        cerr << "compute_cluster_distribution: no clusters found." << endl;
        return;
    }

    // ---- Build histogram ----
    int jmax = imax1 + 1000;
    double Vtot = Lx * Ly * Lz;

    int*    ncluster   = new int[jmax]();
    int*    bincluster = new int[jmax]();
    double* pdf1       = new double[jmax]();
    double* pdf2       = new double[jmax]();

    int nmax1 = 0;
    for (i = 1; i <= label; i++) {
        if (LL[i] > 0) {
            ni = LL[i];
            if (ni >= jmax) continue;
            if (ni > nmax1) nmax1 = ni;
            ++ncluster[ni];
            if (ni > 1) ++nclust;
        }
    }
    cout << " nmax1 = " << nmax1 << endl;
    cout << " vol (cell volume) = " << vol << endl;

    // ---- Averages ----
    double avencluster = 0.0, avevolume = 0.0;
    for (i = 1; i < jmax; i++) {
        if (ncluster[i] > 0) {
            avencluster += ncluster[i] * i;
            avevolume   += ncluster[i] * (double)i * i;
        }
    }
    cout << " The average cluster size of crystalline domains is :  "
         << avencluster / labelcount << endl;
    avevolume = avevolume * vol / nemcount;
    cout << " The average volume of crystalline domains is :  "
         << avevolume << endl;

    // ---- Coarse binning ----
    for (i = 1; i <= 5 && i < jmax; i++)
        bincluster[i] = ncluster[i];

    if (7 < jmax) {
        if (6 < jmax) bincluster[7] += ncluster[6];
        bincluster[7] += ncluster[7];
    }
    if (10 < jmax)
        for (i = 8; i <= 10 && i < jmax; i++)
            bincluster[10] += ncluster[i];

    for (i = 20; i <= 70 && i < jmax; i += 10)
        for (j = i-9; j <= i && j <= imax1 && j < jmax; j++)
            bincluster[i] += ncluster[j];

    if (100 < jmax)
        for (j = 71; j <= 100 && j <= imax1 && j < jmax; j++)
            bincluster[100] += ncluster[j];

    for (i = 200; i <= 700 && i < jmax; i += 100)
        for (j = i-99; j <= i && j <= imax1 && j < jmax; j++)
            bincluster[i] += ncluster[j];

    if (1000 < jmax)
        for (j = 701; j <= 1000 && j <= imax1 && j < jmax; j++)
            bincluster[1000] += ncluster[j];

    for (i = 2000; i < jmax; i += 1000)
        for (j = i-999; j <= i && j <= imax1 && j < jmax; j++)
            bincluster[i] += ncluster[j];

    // ---- PDFs ----
    int n1 = 0, n2_bin = 0;
    bincluster[0] = 1;

    for (i = 1; i < jmax; i++) {
        if (bincluster[i] > 0) {
            n2_bin = i;
            double width = (n2_bin - n1) > 0 ? (double)(n2_bin - n1) : 1.0;
            pdf1[i] = bincluster[i] / width;
            pdf2[i] = ((n1 + n2_bin) * pdf1[i] * 0.5) / nemcount;
            n1 = n2_bin;
        }
    }

    // ---- Output ----
    clusterfile << "clustersize volume vol/Vtot Ncluster[i] pdf1 volume-pdf" << endl;

    for (i = 0; i < jmax; i++) {
        if (bincluster[i] > 0) {
            clusterfile << i << " "
                        << i * vol << " "
                        << i * vol / Vtot << " "
                        << bincluster[i] << " "
                        << pdf1[i] / labelcount << " "
                        << pdf2[i] << " "
                        << pdf2[i] * crystallinity << endl;
        }
    }

    clusterfile << "the total number of clusters with at least 2 elements is "
                << nclust << " nematic grid elements" << endl;
    clusterfile << " The total number of independent crystalline domains is:  "
                << labelcount << endl;
    clusterfile << " The average cluster size of crystalline domains is :  "
                << avencluster / labelcount << endl;
    clusterfile << " The average volume fraction of crystalline domains is :  "
                << avevolume / Vtot << endl;
    clusterfile << " The total number of crystalline grid elements is:  "
                << nemcount << endl;
    clusterfile << " S_Ree globalNem Vtotal Ngrid crystallinity "
                   "avenclust avenclust/Ngrid avevol avevol/Vtot" << endl;

    double avenclust = avencluster / labelcount;
    clusterfile << S_Ree << " "
                << globalNem << " "
                << Vtot << " "
                << sumcount1 << " "
                << crystallinity << " "
                << avenclust << " "
                << (sumcount1 > 0 ? avenclust / sumcount1 : 0.0) << " "
                << avevolume << " "
                << avevolume / Vtot << endl;

    delete[] ncluster;
    delete[] bincluster;
    delete[] pdf1;
    delete[] pdf2;
}

/*********************************************************************************************************************/
// orderparameter1
/*********************************************************************************************************************/
void orderparameter1(int npar, double **um, double lambda[], double eevp[4],
                     double eevm[4], double eevo[4], double *S, double director[3])
{
    double a11=0,a12=0,a13=0,a22=0,a23=0,a33=0,aa11,aa22;
    double cc1,cc0,pp,qq,phia,rr,vnorm;
    int i;

    for (i = 1; i <= npar; i++) {
        a11 += um[i][0]*um[i][0];
        a12 += um[i][0]*um[i][1];
        a13 += um[i][0]*um[i][2];
        a22 += um[i][1]*um[i][1];
        a23 += um[i][1]*um[i][2];
        a33 += um[i][2]*um[i][2];
    }
    a11/=npar; a12/=npar; a13/=npar;
    a22/=npar; a23/=npar; a33/=npar;

    cc1 = a12*a12+a13*a13+a23*a23-a11*a22-a11*a33-a22*a33;
    cc0 = a11*a22*a33+2.0*a12*a13*a23-a12*a12*a33-a13*a13*a22-a23*a23*a11;
    pp  = 0.25*(1.0+3.0*cc1);
    qq  = (27.0*cc0+9.0*cc1+2.0)/16.0;
    phia = acos(qq/sqrt(pow(pp,3)))/3.0;
    rr   = 2.0*sqrt(pp);

    lambda[1] = rr*cos(phia);
    lambda[2] = rr*cos(phia-2.0*pi2/3.0);
    lambda[3] = rr*cos(phia-pi2/3.0);

    auto computeEigvec = [&](double lam, double ev[4]) {
        aa11 = a11-onethird-2.0*onethird*lam;
        aa22 = a22-onethird-2.0*onethird*lam;
        ev[1] = a12*a23-a13*aa22;
        ev[2] = a13*a12-a23*aa11;
        ev[3] = aa11*aa22-a12*a12;
        vnorm = sqrt(ev[1]*ev[1]+ev[2]*ev[2]+ev[3]*ev[3]);
        if (vnorm > 1.0e-12) { ev[1]/=vnorm; ev[2]/=vnorm; ev[3]/=vnorm; }
        else                  { ev[1]=ev[2]=ev[3]=0.0; }
    };
    computeEigvec(lambda[1], eevp);
    computeEigvec(lambda[2], eevm);
    computeEigvec(lambda[3], eevo);

    double S0=0, S_nem=lambda[1];
    director[0]=eevp[1]; director[1]=eevp[2]; director[2]=eevp[3];
    S0 = fabs(lambda[1]);
    if (fabs(lambda[2])>S0) { S0=fabs(lambda[2]); S_nem=lambda[2]; director[0]=eevm[1]; director[1]=eevm[2]; director[2]=eevm[3]; }
    if (fabs(lambda[3])>S0) { S0=fabs(lambda[3]); S_nem=lambda[3]; director[0]=eevo[1]; director[1]=eevo[2]; director[2]=eevo[3]; }

    double ss = sqrt(2.0/3.0*(lambda[1]*lambda[1]+lambda[2]*lambda[2]+lambda[3]*lambda[3]));
    nematicfile << S_nem << " " << ss << " " << lambda[1] << " " << lambda[2] << " " << lambda[3]
                << " " << director[0] << " " << director[1] << " " << director[2] << endl;
    nematicfile << eevp[1] << " " << eevp[2] << " " << eevp[3] << endl;
    nematicfile << eevm[1] << " " << eevm[2] << " " << eevm[3] << endl;
    nematicfile << eevo[1] << " " << eevo[2] << " " << eevo[3] << endl;
    nematicfile << "ss= " << ss << endl;
    *S = S_nem;
}

/*********************************************************************************************************************/
// orderparameter
/*********************************************************************************************************************/
void orderparameter(int npar, double um[][4], double lambda[], double eevp[4],
                    double eevm[4], double eevo[4], double *S, double director[3])
{
    double a11=0,a12=0,a13=0,a22=0,a23=0,a33=0,aa11,aa22;
    double cc1,cc0,pp,qq,phia,rr,vnorm;
    int i;

    for (i = 1; i <= npar; i++) {
        a11 += um[i][1]*um[i][1];
        a12 += um[i][1]*um[i][2];
        a13 += um[i][1]*um[i][3];
        a22 += um[i][2]*um[i][2];
        a23 += um[i][2]*um[i][3];
        a33 += um[i][3]*um[i][3];
    }
    a11/=npar; a12/=npar; a13/=npar;
    a22/=npar; a23/=npar; a33/=npar;

    cc1 = a12*a12+a13*a13+a23*a23-a11*a22-a11*a33-a22*a33;
    cc0 = a11*a22*a33+2.0*a12*a13*a23-a12*a12*a33-a13*a13*a22-a23*a23*a11;
    pp  = 0.25*(1.0+3.0*cc1);
    qq  = (27.0*cc0+9.0*cc1+2.0)/16.0;
    phia = acos(qq/sqrt(pow(pp,3)))/3.0;
    rr   = 2.0*sqrt(pp);

    lambda[1] = rr*cos(phia);
    lambda[2] = rr*cos(phia-2.0*pi2/3.0);
    lambda[3] = rr*cos(phia-pi2/3.0);

    auto computeEigvec = [&](double lam, double ev[4]) {
        aa11 = a11-onethird-2.0*onethird*lam;
        aa22 = a22-onethird-2.0*onethird*lam;
        ev[1] = a12*a23-a13*aa22;
        ev[2] = a13*a12-a23*aa11;
        ev[3] = aa11*aa22-a12*a12;
        vnorm = sqrt(ev[1]*ev[1]+ev[2]*ev[2]+ev[3]*ev[3]);
        if (vnorm > 1.0e-12) { ev[1]/=vnorm; ev[2]/=vnorm; ev[3]/=vnorm; }
        else                  { ev[1]=ev[2]=ev[3]=0.0; }
    };
    computeEigvec(lambda[1], eevp);
    computeEigvec(lambda[2], eevm);
    computeEigvec(lambda[3], eevo);

    double S0=0, S_nem=lambda[1];
    director[0]=eevp[1]; director[1]=eevp[2]; director[2]=eevp[3];
    S0 = fabs(lambda[1]);
    if (fabs(lambda[2])>S0) { S0=fabs(lambda[2]); S_nem=lambda[2]; director[0]=eevm[1]; director[1]=eevm[2]; director[2]=eevm[3]; }
    if (fabs(lambda[3])>S0) { S0=fabs(lambda[3]); S_nem=lambda[3]; director[0]=eevo[1]; director[1]=eevo[2]; director[2]=eevo[3]; }

    double ss = sqrt(2.0/3.0*(lambda[1]*lambda[1]+lambda[2]*lambda[2]+lambda[3]*lambda[3]));
    nematicfile << S_nem << " " << ss << " " << lambda[1] << " " << lambda[2] << " " << lambda[3]
                << " " << director[0] << " " << director[1] << " " << director[2] << endl;
    nematicfile << eevp[1] << " " << eevp[2] << " " << eevp[3] << endl;
    nematicfile << eevm[1] << " " << eevm[2] << " " << eevm[3] << endl;
    nematicfile << eevo[1] << " " << eevo[2] << " " << eevo[3] << endl;
    nematicfile << lambda[1] << " " << lambda[2] << " " << lambda[3] << endl;
    nematicfile << "ss= " << ss << endl;
    *S = S_nem;
}


// Wrap scalar u into [0, L) (shifted box)
inline double wrap_to_zero(double u, double lo, double L)
{
    // u is e.g. xu[id], lo is xlo, L is Lx
    double x = u - lo;             // shift origin to 0
    x -= floor(x / L) * L;         // wrap into [0, L)
    return x;                      // in [0,L)
}

// If you want the wrapped coordinate back in [lo, hi):
inline double wrap_to_box(double u, double lo, double L)
{
    double x = wrap_to_zero(u, lo, L);
    return x + lo;                 // in [lo, lo+L)
}

inline double min_image(double dx, double L)
{
    if (dx >  0.5 * L) dx -= L;
    if (dx < -0.5 * L) dx += L;
    return dx;
}
