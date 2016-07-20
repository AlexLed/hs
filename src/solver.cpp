#include <QtCore>
#include "datadef.h"
#include "math.h"

//Parameters
double Theta0d = 45, //atan(100)/Pi*180,   //Angle between bisectrix and leading edge 1
       Theta0 = Theta0d/180*Pi,
       Theta1d = Theta0d,    //Angle between bisectrix and leading edge 2
       Theta1 = Theta1d/180*Pi,
       z0,
       Betad = 0, Beta,//Angle between bisectrix and U
       Gamma = 1.4,    //Cp/Cv
       Pr = 0.72,      //Prandtl number
       Hw = 0.5,       //Enthalpy on the wall
       Vw = 0,        //Suction/blow through the wall
       khi = 1000,     //parameter of interaction
       Omega = 1,      //Dependence of viscosity on temperature: mu ~ T ^ Omega
       relU=0.7, relP=0.3, relP2=0.01, a0=0.05, //relaxation coefficients
       AccuracyU=0.00001, AccuracyPressure=0.000001; //approximation accuracy  for friction and pressure
//Grid
int    imax = 81, jmax = 201, kmax = 80; //grid dimensions (x, y, z correspondingly)
double L = 1.0,//length of wing
       h = 15, //height
       dx = L/(imax-1), dy = h/(jmax-1), dz = 2.0/(kmax-1),  //grid spacing
       dt = 0.0005;  //time step
/*****************************************************************************
U,W,V,H[X][Y][Z] - profiles of longitudinal U, transverse W, normal V velocities and enthalpy H
Us,Ws,Hs[X][Y][Z] - arrays for relaxation and checking approximations
X[imax] - variable X or R [0;1]
Y[jmax] - variable Y [0;h]
Z[kmax] - variable Z or Theta [-1;1]
Delta[imax][kmax], DeltaE[imax][kmax] - displacement thickness of boundary layer
DeltaR[imax][kmax] - array for delta relaxation
P[imax][kmax], P0[imax][kmax] - pressure field
Prel[imax][kmax] - array for pressure relaxation
Pd[kmax] - modification of relaxation
AU[jmax], BU[jmax], BW[jmax], AH[jmax], BH[jmax], AP[jmax], BP[jmax] - marching coefficients
*****************************************************************************/
double *X, *Y, *Z, *Ze, **Body;
double ***U, ***W, ***H, ***V, ***T, ***Mach, ***Ro, ***Us, ***Ws, ***Hs, ***Vs;
double **Tau, **TauU, **TauW, **TauH, **DeltaDistribution, **PressureDistribution,
       **TauUDistribution, **TauWDistribution, **TauHDistribution, **DisturbancesDistribution;
double **Delta, **DeltaE, **DeltaRel, **Delta99; //arrays for displacement thickness
double **P, **P0, **Pd, **Prel; //arrays for pressure
double *AU, *BU, *BW, *AH, *BH, *AP, *BP, *Mu, *J; //progonka coefficients and others
int    i, j, k, l, n; //indexes: i -> x, r;  j -> y, eta;  k -> z, theta; n -> t
int    N = 700, LX = 0, LY = 15, LZ = 0 /*kmax*2/3*/, NIterations = 0 /*Number of iterations*/, MaxIterations = 1000;
bool   converged = false,
       thomas_method = true,
       cartesian = true,
       is_first_order_theta = false,
       is_mod_rel = true,
       draw_grid = true;
QString status;
double Simps(double, int, int), SystemC1(double P2, double alpha);
double RungeKutta2(double ***F, int i, int k, double dy, double dF0);

double ***NewArray3D(int N1, int N2, int N3) {
    double ***array;
    array = new double**[N1];
    for (i = 0; i < N1; i++) {
        array[i] = new double*[N2];
        for (j = 0; j < N2; j++) {
            array[i][j] = new double [N3];
            for (k = 0; k < N3; k++) {
                array[i][j][k] = 0;
            }
        }
    }
    return array;
}
double  **NewArray2D(int N1, int N2) {
    double **array;
    array = new double*[N1];
    for (i = 0; i < N1; i++) {
        array[i] = new double[N2];
        for (k = 0; k < N2; k++) {
            array[i][k] = 0;
        }
    }
    return array;
}
void SetZe(int imax, double dx, double Theta0) { //Form of leading edge (sea also in Edge() and Nos())
    for (i = 0; i < imax+1; i++) {
        //Ze[i] = (cbrt(i*dx*2-1)+1)*0.5+i*dx/3;
        Ze[i] = i*dx*tan(Theta0);
        //Ze[i] = i*dx*exp(i*dx)*0.5;
        //Ze[i] = log(i*dx+1);
        //Ze[i] = tan(Theta0)*atan(i*dx*20);
        //if(i*dx < 0.5) Ze[i] = i*dx*tan(Theta0);
        //else Ze[i] = (i*dx-0.5)*tan(Theta1)+tan(Theta0)*0.5;
    }
    z0 = Ze[int(imax/L)];
    //trail
    for (i = imax/L; i < imax+1; i++) {
        Ze[i] = z0;// + (i*dx-1)*tan(5/180.0*Pi);
    }
    for (i = imax/4; i < imax*3/4; i++)
        for (k = kmax/4; k < kmax*3/4; k++)
            Body[i][k] = 0;//(1+cos((X[i]+0.5)*4*Pi))*(1+cos(Z[k]*2*Pi))/10;
    //for (i = 0; i < imax; i++)
    //   for (k = 0; k < kmax; k++)
    //       Body[i][k] = pow(X[i], 0.75)*pow(1-Z[k]*Z[k], 9.25)/2;
}
void CreateArrays(bool is_interpolated, int imax, int jmax, int kmax, double dx, double dz, int N) {
    if(!is_interpolated) {
        Us       = NewArray3D(imax, jmax, kmax);
        Vs       = NewArray3D(imax, jmax, kmax);
        Ws       = NewArray3D(imax, jmax, kmax);
        Hs       = NewArray3D(imax, jmax, kmax);
        Prel     = NewArray2D(imax, kmax);
        DeltaRel = NewArray2D(imax, kmax);
    }
    U    = NewArray3D(imax+1, jmax+1, kmax+1);
    W    = NewArray3D(imax+1, jmax+1, kmax+1);
    V    = NewArray3D(imax+1, jmax+1, kmax+1);
    H    = NewArray3D(imax+1, jmax+1, kmax+1);
    T    = NewArray3D(imax+1, jmax+1, kmax+1);
    Ro   = NewArray3D(imax,   jmax,   kmax);
    Mach = NewArray3D(imax,   jmax,   kmax);

    //arrays for coordinates
    X  = new double [imax];
    Y  = new double [jmax];
    Z  = new double [kmax];
    Ze = new double [imax+1];

    Body   = NewArray2D(imax+1, kmax);
    P      = NewArray2D(imax+1, kmax+1);
    P0     = NewArray2D(imax,   kmax);
    Pd     = NewArray2D(imax,   kmax);
    Delta  = NewArray2D(imax,   kmax);
    DeltaE = NewArray2D(imax,   kmax);
    Delta99= NewArray2D(imax,   kmax);

    //temporary arrays
    Tau  = NewArray2D(jmax, kmax);
    TauU = NewArray2D(jmax, kmax);
    TauW = NewArray2D(jmax, kmax);
    TauH = NewArray2D(jmax, kmax);

    Mu = new double [jmax];
    J  = new double [kmax];
    AU = new double [jmax];
    BU = new double [jmax];
    BW = new double [jmax];
    AH = new double [jmax];
    BH = new double [jmax];
    AP = new double [kmax];
    BP = new double [kmax];

    DeltaDistribution    = NewArray2D(N, N);
    PressureDistribution = NewArray2D(N, N);
    TauUDistribution     = NewArray2D(N, N);
    TauWDistribution     = NewArray2D(N, N);
    TauHDistribution     = NewArray2D(N, N);

    DisturbancesDistribution = NewArray2D(360, kmax);

    //Coordinates
    for (i = 0; i < imax; i++) X[i] = dx*i; //X=[0;L]
    for (k = 0; k < kmax; k++) Z[k] = -1 + dz*k; //Z=[-1;1]

    //Form of leading edge
    SetZe(imax, dx, Theta0);

    //qDebug() << "Arrays mem size:" <<
    //            sizeof(double)*(imax*kmax*(jmax*11+8)+5*N*N)/1024/1024 << "MB";
}
void SetBoundaryConditions(bool is_cartesian, double Hw) {
    for (i = 0; i < imax; i++)
        for (k = 0; k < kmax; k++) {
            //Y = 0
            U[i][0][k] = 0;
            W[i][0][k] = 0;
            if( true || ((k > kmax*1/4) && (k < kmax*3/4)) || (k == 0) || (k==kmax-1))
            //if( ((k > kmax*1/4) && (k < kmax*3/4)) )
                V[i][0][k] = Vw;
            else
                V[i][0][k] = 0;
            H[i][0][k] = Hw;
            T[i][0][k] = H[i][0][k]-U[i][0][k]*U[i][0][k]-W[i][0][k]*W[i][0][k];
            //Y = inf
            if(is_cartesian) {
                U[i][jmax-1][k] = cos(Beta);
                W[i][jmax-1][k] = sin(Beta);
            }
            else {
                U[i][jmax-1][k] = cos(Theta0*Z[k]-Beta);
                W[i][jmax-1][k] = -sin(Theta0*Z[k]-Beta);
            }
            V[i][jmax-1][k] = 0;
            H[i][jmax-1][k] = 1;
            T[i][jmax-1][k] = H[i][jmax-1][k]-U[i][jmax-1][k]*U[i][jmax-1][k]-W[i][jmax-1][k]*W[i][jmax-1][k];
        }
    Mu[0] = pow(fabs(Hw), Omega-1);  Mu[jmax-1] = 0; //viscocity coefficient
    AU[0] = 0; BU[0] = 0; BW[0] = 0; AH[0] = 0; BH[0] = Hw; //marching coefficients
}
void Simm() {
    QString method = "RK"; //Progonka, RK, Ust
    double /**AW, *BW, Eu, Ew, Eh, vdy2, dW, dWs, */Integral;
    for (j = 1; j < jmax-1; j++) {
        U[0][j][0] = U[0][jmax-1][0]*pow((double)j/(jmax-1), 0.25);
        //W[0][j][0] = W[0][jmax-1][0]*pow((double)j/(jmax-1), 0.25);
        H[0][j][0] = H[0][0][0]*(1-pow((double)j/(jmax-1), 0.25))+H[0][jmax-1][0]*pow((double)j/(jmax-1), 0.25);
    }
    if (method == "RK") {
        qDebug() << "RK2 Umax =" << RungeKutta2(U,0,0,dy,(U[0][1][0]-U[0][0][0])/dy);
        //qDebug() << "RK2 Hmax =" << RungeKutta2(H,0,0,dy,(H[0][1][0]-H[0][0][0])/dy);
    }
    if (method == "Ust") {
        dt = dy*dy/4.0001;
        //System C00 without W00
        P[0][0] = 1;   Delta[0][0] = 0;
        U[0][0][0] = 0; W[0][0][0] = 0; H[0][0][0] = Hw; V[0][0][0] = 0;
        U[0][jmax-1][0] = 1; W[0][jmax-1][0] = 0; H[0][jmax-1][0] = 1;
        Us[0][2][0] = 1;
        for (NIterations = 0; (U[0][2][0] != Us[0][2][0]) && NIterations < 1000000 ; NIterations++) {
            QCoreApplication::processEvents(); //to prevent freezing of gui
            for (j = 0; j < jmax; j++) {
                Us[0][j][0] = U[0][j][0];
                Hs[0][j][0] = H[0][j][0];
            }
            Integral = H[0][0][0]-U[0][0][0]*U[0][0][0];
            for (j = 1; j < jmax-1; j+=2) Integral +=  4*(H[0][j][0]-U[0][j][0]*U[0][j][0]) + 2*(H[0][j+1][0]-U[0][j+1][0]*U[0][j+1][0]);
            Delta[0][0] = sqrt((Gamma-1)/Gamma*0.5)/P[0][0]*Integral*dy/3.0;
            P[0][0] = (1-relP)*P[0][0] + relP*(Gamma+1)*9/32.0*Delta[0][0]*Delta[0][0];
            for (j = 1; j < jmax-1; j++) {
                V[0][j][0] = V[0][j-1][0] - dy*(0.25*Us[0][j][0])/P[0][0];
                U[0][j][0] = Us[0][j][0]  - dt*(V[0][j][0]*(Us[0][j+1][0]-Us[0][j-1][0])/2/dy
                                               - (Gamma-1)/4.0/Gamma/P[0][0]*(Hs[0][j][0]-Us[0][j][0]*Us[0][j][0])
                                               - (Us[0][j+1][0]-2*Us[0][j][0]+Us[0][j-1][0])/dy/dy);
                H[0][j][0] = Hs[0][j][0]  - dt*(V[0][j][0]*(Hs[0][j+1][0]-Hs[0][j-1][0])/2/dy*Pr
                                               + (1-Pr)*(Us[0][j+1][0]*Us[0][j+1][0]-2*Us[0][j][0]*Us[0][j][0]+Us[0][j-1][0]*Us[0][j-1][0])/dy/dy
                                               - (Hs[0][j+1][0]-2*Hs[0][j][0]+Hs[0][j-1][0])/dy/dy);
            }
        }
//        qDebug() << "Delta[0][0]"<< Delta[0][0];
//        qDebug() << "P[0][0]"<< P[0][0];
//        qDebug() << "C00 is finished\n";
        //System C10 without W[1][0]
        //boundary conditions
        P[1][0] = 0;
        U[1][0][0] = 0; W[1][0][0] = 0; H[1][0][0] = 0; V[1][0][0] = 0;
        U[1][jmax-1][0] = 0; W[1][jmax-1][0] = 0; H[1][jmax-1][0] = 0;
        for (j = 1; j < jmax-1; j++) {
            U[1][j][0] = 0;
            W[1][j][0] = 0;
            H[1][j][0] = 0;
            V[1][j][0] = 0;
        }
        dt = dy*dy/4.1;
        Us[1][2][0] = 1;
        for (NIterations = 0; (U[1][2][0] != Us[1][2][0]) && NIterations < 0 ; NIterations++) {
            QCoreApplication::processEvents(); //to prevent freezing of gui
            for (j = 0; j < jmax; j++) {
                Us[1][j][0] = U[1][j][0];
                Hs[1][j][0] = H[1][j][0];
            }
            Integral = (H[1][0][0]-2*U[0][0][0]*U[1][0][0])-P[1][0]/P[0][0]*(H[0][0][0]-U[0][0][0]*U[0][0][0]);
            for (j = 1; j < jmax-1; j+=2)
                Integral += 4*((H[1][j][0]-2*U[0][j][0]*U[1][j][0])-P[1][0]/P[0][0]*(H[0][j][0]-U[0][j][0]*U[0][j][0]))
                          + 2*((H[1][j+1][0]-2*U[0][j+1][0]*U[1][j+1][0])-P[1][0]/P[0][0]*(H[0][j+1][0]-U[0][j+1][0]*U[0][j+1][0]));
            Delta[1][0] = sqrt((Gamma-1)/Gamma*0.5)/P[0][0]*Integral*dy/3.0;
            P[1][0] = (1-relP)*P[1][0]
                    + relP *(Gamma+1)/16.0*Delta[0][0]*(9*Delta[0][0]-15*Delta[1][0]);
            for (j = 1; j < jmax-1; j++) {
                V[1][j][0] = V[1][j-1][0] - dy*(0.25*Us[0][j][0]/P[0][0]*(1-P[1][0]/P[0][0])
                                          - 7.0*0.25*Us[1][j][0]/P[0][0]);
                U[1][j][0] = Us[1][j][0] - dt*(0*V[0][j][0]*(Us[1][j+1][0]-Us[1][j-1][0])/2/dy
                                         +     V[1][j][0]*(Us[0][j+1][0]-Us[0][j-1][0])/2/dy
                                         - 2.0*Us[0][j][0]*Us[1][j][0]/P[0][0]
                                         - (Gamma-1)/4.0/Gamma/P[0][0]*
                                         ((Hs[0][j][0]-Us[0][j][0]*Us[0][j][0])*(1+3.0*P[1][0]/P[0][0])
                                         +(Hs[1][j][0]-2.0*Us[0][j][0]*Us[1][j][0]))
                                         - (Us[1][j+1][0]-2*Us[1][j][0]+Us[1][j-1][0])/dy/dy);
                H[1][j][0] = Hs[1][j][0] - dt*(V[0][j][0]*(Hs[1][j+1][0]-Hs[1][j-1][0])/2/dy
                                         +     V[1][j][0]*(Hs[0][j+1][0]-Hs[0][j-1][0])/2/dy
                                         - 2.0*Us[0][j][0]*Hs[1][j][0]/P[0][0]
                                         + (1-Pr)/Pr*2.0*(Us[0][j+1][0]*Us[1][j+1][0]-2*Us[0][j][0]*Us[1][j][0]+Us[0][j-1][0]*Us[1][j-1][0])/dy/dy
                                         - (Hs[1][j+1][0]-2*Hs[1][j][0]+Hs[1][j-1][0])/dy/dy/Pr);
            }
        }
//        qDebug() << U[1][1][0];
//        qDebug() << "Delta[1][0]"<< Delta[1][0];
//        qDebug() << "P[1][0]"<< P[1][0];
//        qDebug() << "C10 is finished";
    }
    else if (method == "RK2") {
        //System C00 without W00
        P[0][0] = 1;   Delta[0][0] = 0;
        U[0][0][0] = 0; W[0][0][0] = 0; H[0][0][0] = Hw; V[0][0][0] = 0;
        U[0][jmax-1][0] = 1; W[0][jmax-1][0] = 0; H[0][jmax-1][0] = 1;
        //double Umax = 1, Hmax = 1;
        double dU0, dU0s, dU0n, Um, Ums;//, dH0, dH0s;
        for (NIterations = 0; NIterations < 1 ; NIterations++) {
            QCoreApplication::processEvents(); //to prevent freezing of gui
            Integral = H[0][0][0]-U[0][0][0]*U[0][0][0];
            for (j = 1; j < jmax-1; j+=2) Integral +=  4*(H[0][j][0]-U[0][j][0]*U[0][j][0]) + 2*(H[0][j+1][0]-U[0][j+1][0]*U[0][j+1][0]);
            Delta[0][0] = sqrt((Gamma-1)/Gamma*0.5)/P[0][0]*Integral*dy/3.0;
            P[0][0] = (1-relP)*P[0][0] + relP*(Gamma+1)*9/32.0*Delta[0][0]*Delta[0][0];

            dU0 = 0.2;
            for (int n = 0; n < 30; n++) Um = RungeKutta2(U,0,0,dy,dU0);
            qDebug() << dU0 << Um;
            dU0s = dU0;
            dU0 +=0.1;
            Ums = Um;
            for (int n = 0; n < 30; n++) Um = RungeKutta2(U,0,0,dy,dU0);
            qDebug() << dU0 << RungeKutta2(U,0,0,dy,dU0);
            for (int n = 0; n < 50; n++) {
                dU0n = dU0s+(dU0-dU0s)*(1-Ums)/(Um-Ums);
                dU0s = dU0;
                dU0 = dU0n;
                for (int k = 0; k < 30; k++) Um = RungeKutta2(U,0,0,dy,dU0);
                qDebug() << dU0 << Um;
            }
        }
        //while ((fabs(U[0][jmax-1][0]-Umax) > 0.0001) && (fabs(H[0][jmax-1][0]-Hmax) > 0.0001));
    }
}
void Edge(double Theta, int I, bool is_cartesian, int kromka) { //computing on the leading edges
    double Integral, E, Eu, Ew, Eh, vdy2; converged = false;
    double XdZe; int K;
    if(I == 0) XdZe = 1/tan(Theta0);
    else XdZe = X[I]/Ze[I];
    dt = dy*dy/4.00001;
    kromka == -1 ? K = 0 : K = kmax-1;
    //setting initial distribution and boundary conditions
    double y;
    for (j = 1; j < jmax-1; j++) {
        y = pow((double)j/(jmax-1), 0.25);
        U[I][j][K] = (double)y*U[I][jmax-1][K];
        W[I][j][K] = 0;//y*W[I][jmax-1][K];
        H[I][j][K] = H[I][0][K]*(1-y)+H[I][jmax-1][K]*y;
    }
    Prel[I][K] = P[I][K]+0.1;
    if (thomas_method || 1) //Progonka
        for (NIterations = 0; fabs(Prel[I][K]-P[I][K])>0.0000000001 && NIterations < 10000 ; NIterations++) {
        QCoreApplication::processEvents(); //to prevent freezing of gui
        Prel[I][K] = P[I][K]; //for checking approximation
//		H[I][0][K] = H[I][1][K];
        for (j = 0; j < jmax-1; j++)  T[I][j][K] = H[I][j][K]-U[I][j][K]*U[I][j][K]-W[I][j][K]*W[I][j][K];
//        Integral = T[I][0][K];
//        for (j = 1; j < jmax-1; j+=2) Integral += 4*T[I][j][K]+2*T[I][j+1][K];
        Integral = 0;
        for (j = 1; j < jmax-1; j+=2) Integral += (T[I][j-1][K]+4*T[I][j][K]+T[I][j+1][K])*dy/3.0;
        if (N%2 == 1) Integral += (T[I][j-1][K]+T[I][j][K])*dy/2.0;
        if (is_cartesian) {
            Delta[I][K] = cbrt(sqrt((Gamma-1)/Gamma*0.5)*Integral*8.0/
                              (9.0*(Gamma+1)*pow(XdZe*(Ze[I+1]-Ze[I])/dx, 2.0)));
            P[I][K] = (Gamma+1)*0.5*pow(1.5*Delta[I][K]*XdZe*(Ze[I+1]-Ze[I])/dx, 2.0);
        }
        else {
            Delta[I][K] = cbrt(sqrt((Gamma-1)/Gamma*0.5)*Integral*dy/3.0*8.0/(9.0*(Gamma+1)*pow(sin(Theta*Z[K]-Beta)/Theta, 2.0)));
            P[I][K] = 9.0/8.0*(Gamma+1)*pow(Delta[I][K]*sin(Theta*Z[K]-Beta)/Theta, 2.0);
        }
        for (j = 1; j < jmax-1; j++) Mu[j] = pow(fabs(T[I][j][K]), Omega-1); //viscocity coefficient
        for (j = 1; j < jmax-1; j++) {
            E = (Gamma-1)/Gamma*0.5/P[I][K]*dy*dy*T[I][j][K];
            if (is_cartesian) {
                Eu = z0*E*XdZe*(Ze[I+1]-Ze[I])/dx;
                Ew = -Z[K]*E*XdZe;
                V[I][j][K] = V[I][j-1][K]+Z[K]*dy/4.0/P[I][K]*
                        ((W[I][j][K]+W[I][j-1][K])*XdZe-z0*Z[K]*XdZe*(Ze[I+1]-Ze[I])/dx*(U[I][j][K]+U[I][j-1][K]));
            }
            else {
                Eu = 0;
                Ew = -Z[K]*E;
                V[I][j][K] = V[I][j-1][K]+Z[K]*dy/4.0/P[I][K]*(W[I][j][K]+W[I][j-1][K]);
            }
            Eh = (Pr-1)/Pr*(U[I][j+1][K]*U[I][j+1][K] - 2*U[I][j][K]*U[I][j][K] + U[I][j-1][K]*U[I][j-1][K] +
                            W[I][j+1][K]*W[I][j+1][K] - 2*W[I][j][K]*W[I][j][K] + W[I][j-1][K]*W[I][j-1][K]);
            if (Omega == 1) {
                vdy2 = V[I][j][K]*dy*0.5;
                AU[j] = (1-vdy2) / (2-AU[j-1]*(1+vdy2));
                BU[j] = (BU[j-1]*(1+vdy2) + Eu) / (2-AU[j-1]*(1+vdy2));
                BW[j] = (BW[j-1]*(1+vdy2) + Ew) / (2-AU[j-1]*(1+vdy2));
                AH[j] = (1/Pr-vdy2) / (2/Pr - AH[j-1] * (1/Pr+vdy2));
                BH[j] = (BH[j-1]*(1/Pr+vdy2)+Eh) / (2/Pr - AH[j-1] * (1/Pr+vdy2));
            }
            else {
                vdy2 = V[I][j][K]*dy*2;
                AU[j] = (Mu[j+1]+4*Mu[j]-Mu[j-1]-vdy2)/(8*Mu[j]+AU[j-1]*(Mu[j+1]-4*Mu[j]-Mu[j-1]-vdy2));
                BU[j] = (-BU[j-1]*(Mu[j+1]-4*Mu[j]-Mu[j-1]-vdy2)+4*Eu)/(8*Mu[j]+AU[j-1]*(Mu[j+1]-4*Mu[j]-Mu[j-1]-vdy2));
                BW[j] = (-BW[j-1]*(Mu[j+1]-4*Mu[j]-Mu[j-1]-vdy2)+4*Ew)/(8*Mu[j]+AU[j-1]*(Mu[j+1]-4*Mu[j]-Mu[j-1]-vdy2));
                AH[j] = (Mu[j+1]+4*Mu[j]-Mu[j-1]-vdy2*Pr)/(8*Mu[j]+AH[j-1]*(Mu[j+1]-4*Mu[j]-Mu[j-1]-vdy2*Pr));
                BH[j] = (-BH[j-1]*(Mu[j+1]-4*Mu[j]-Mu[j-1]-vdy2*Pr)+4*Eh*Pr)/(8*Mu[j]+AH[j-1]*(Mu[j+1]-4*Mu[j]-Mu[j-1]-vdy2*Pr));
            }
        }
        for (j = jmax-2; j >= 0; j--) {
            U[I][j][K] = U[I][j+1][K]*AU[j] + BU[j];
            W[I][j][K] = W[I][j+1][K]*AU[j] + BW[j];
            H[I][j][K] = H[I][j+1][K]*AH[j] + BH[j];
        }
    }
    else { //Ustanovlenie //no cartesian system
     for (NIterations = 0; fabs((U[I][2][K]+W[I][2][K])/(Us[I][2][K]+Ws[I][2][K]))>0.00001 && NIterations < 100000; NIterations++) {
        QCoreApplication::processEvents(); //to prevent freezing of gui
        for (j = 0; j < jmax; j++) {
            T[I][j][K]  = H[I][j][K]-U[I][j][K]*U[I][j][K]-W[I][j][K]*W[I][j][K];
            Us[I][j][K] = U[I][j][K];
            Ws[I][j][K] = W[I][j][K];
            Hs[I][j][K] = H[I][j][K];
        }
        Integral = T[I][0][K];
        for (j = 1; j < jmax-1; j+=2) Integral += 4*T[I][j][K]+2*T[I][j+1][K];
        Delta[I][K] = cbrt(sqrt((Gamma-1)/Gamma*0.5)*Integral*dy/3.0/((Gamma+1)*0.5*pow(1.5*Z[K]/Theta*sin(Theta*Z[K]-Beta), 2.0)));
        P[I][K] = (Gamma+1)/2*pow(1.5*Delta[I][K]*sin(Theta*Z[K]-Beta)/Theta, 2.0);
        for (j = 1; j < jmax-1; j++) Mu[j] = pow(fabs(T[I][j][K]), Omega-1); //viscocity coefficient
        for (j = 1; j < jmax-1; j++) {
            V[I][j][K] = V[I][j-1][K] + Z[K]*dy/4.0/P[I][K]*(Ws[I][j][K]+Ws[I][j-1][K]);
            E = (Gamma-1)/Gamma*0.5*T[I][j][K]/P[I][K];
            Eu = 0;
            Ew = -Z[K]*E;
            Eh = (Pr-1)/Pr*(Us[I][j+1][K]*Us[I][j+1][K]-2*Us[I][j][K]*Us[I][j][K]+Us[I][j-1][K]*Us[I][j-1][K]+
                            Ws[I][j+1][K]*Ws[I][j+1][K]-2*Ws[I][j][K]*Ws[I][j][K]+Ws[I][j-1][K]*Ws[I][j-1][K])/dy/dy;
            vdy2 = V[I][j][K]*dy*2;
            U[I][j][K] = Us[I][j][K] - dt*((vdy2-Mu[j+1]+Mu[j-1])*(Us[I][j+1][K]-Us[I][j-1][K])/dy/dy*0.25
                                     - Mu[j]*(Us[I][j+1][K]-2*Us[I][j][K]+Us[I][j-1][K])/dy/dy - Eu);
            W[I][j][K] = Ws[I][j][K] - dt*((vdy2-Mu[j+1]+Mu[j-1])*(Ws[I][j+1][K]-Ws[I][j-1][K])/dy/dy*0.25
                                     - Mu[j]*(Ws[I][j+1][K]-2*Ws[I][j][K]+Ws[I][j-1][K])/dy/dy - Ew);
            H[I][j][K] = Hs[I][j][K] - dt*((vdy2-Mu[j+1]/Pr+Mu[j-1]/Pr)*(Hs[I][j+1][K]-Hs[I][j-1][K])/dy/dy*0.25
                                     - Mu[j]*(Hs[I][j+1][K]-2*Hs[I][j][K]+Hs[I][j-1][K])/dy/dy/Pr - Eh);
        }
    }
 }
}
void Nos(bool is_cartesian) { //computing on the nose
    //set up initial fields
    for (k = 1; k < kmax-1; k++) P[0][k] = P[0][0] + (1+tanh(Z[k]*5))/2*((P[0][kmax-1]-P[0][0]));
    if (!is_cartesian) {
        if (Beta > 0.000001)
            for (j = 1; j < jmax-1; j++)
                for (k = 1; k < kmax-1; k++) {
                    U[0][j][k] = U[0][j][0] + (U[0][j][kmax-1] - U[0][j][0])*(U[0][jmax-1][k]-U[0][jmax-1][0])
                    / (U[0][jmax-1][kmax-1]-U[0][jmax-1][0]);
                    W[0][j][k] = W[0][j][0] + (W[0][j][kmax-1] - W[0][j][0])*(W[0][jmax-1][k]-W[0][jmax-1][0])
                    / (W[0][jmax-1][kmax-1]-W[0][jmax-1][0]);
                    H[0][j][k] = (H[0][j][0]*(kmax-k-1) + H[0][j][kmax-1]*k)/(kmax-1);
                }
        else
        for (j = 1; j < jmax-1; j++)
            for (k = 1; k < kmax-1; k++) {
                U[0][j][k] = U[0][j][kmax-1]*U[0][jmax-1][k]/U[0][jmax-1][kmax-1];
                W[0][j][k] = W[0][j][kmax-1]*W[0][jmax-1][k]/W[0][jmax-1][kmax-1];
                H[0][j][k] = (H[0][j][0]*(kmax-k-1) + H[0][j][kmax-1]*k)/(kmax-1);
            }
    }
    else {
        for (j = 1; j < jmax-1; j++)
            for (k = 1; k < kmax-1; k++) {
                U[0][j][k] = (U[0][j][0]*(kmax-k-1) + U[0][j][kmax-1]*k)/(kmax-1);
                W[0][j][k] = (W[0][j][0]*(kmax-k-1) + W[0][j][kmax-1]*k)/(kmax-1);
                H[0][j][k] = (H[0][j][0]*(kmax-k-1) + H[0][j][kmax-1]*k)/(kmax-1);
            }
    }
   double Integral, A, Aw, E, Eu, Ew, Eh, vdy2, dPdz; converged = false;
   double XdZe, dZedX;
   XdZe = 1/tan(Theta0);
   dZedX = (Ze[1]-Ze[0])/dx;
   dt = dy*dy/4.1;
   for (NIterations = 0; !converged && NIterations < 15000; NIterations++) {
       QCoreApplication::processEvents(); //to prevent freezing of gui
       //arrays for checking approximation and relaxation
       for (j = 0; j < jmax; j++) {
          memcpy(Us[0][j], U[0][j], sizeof(double)*kmax);
          memcpy(Ws[0][j], W[0][j], sizeof(double)*kmax);
          memcpy(Hs[0][j], H[0][j], sizeof(double)*kmax);
       }
      memcpy(Prel[0], P[0], sizeof(double)*kmax);
      for (k = 1; k < kmax-1; k++) {
          for (j = 0; j < jmax-1; j++)
              T[0][j][k] = H[0][j][k]-U[0][j][k]*U[0][j][k]-W[0][j][k]*W[0][j][k];
          Integral = T[0][0][k];
          for (j = 1; j < jmax-1; j+=2) Integral +=  4*T[0][j][k] + 2*T[0][j+1][k];
          Delta[0][k] = sqrt((Gamma-1)/Gamma*0.5)/P[0][k]*Integral*dy/3.0;
      }
      //Pressure
      if (is_cartesian)
          for (k = 1; k < kmax-1; k++)
              P[0][k] = (Gamma+1.0)/2.0*pow((0.75*(1-Z[k]*Z[k])*(Delta[0][k]+Body[0][k])
                        -Z[k]*(1-Z[k]*Z[k])*(Delta[0][k+1]-Delta[0][k-1]+Body[0][k+1]-Body[0][k-1])/dz*0.5
                        +1.5*Z[k]*Z[k]*(Delta[0][k]+Body[0][k]))*cos(Beta)
                        +sin(Beta)/z0*((1-Z[k]*Z[k])*(Delta[0][k+1]-Delta[0][k-1]+Body[0][k+1]-Body[0][k-1])/dz*0.5
                        -1.5*Z[k]*(Delta[0][k]+Body[0][k])),2.0);
      else
          for (k = 1; k < kmax-1; k++)
              P[0][k] = (Gamma+1.0)/2.0*pow(0.75*(1-Z[k]*Z[k])*Delta[0][k]*cos(Theta0*Z[k]-Beta)+sin(Theta0*Z[k]-Beta)/Theta0
                *(1.5*Z[k]*Delta[0][k]-(1-Z[k]*Z[k])*(Delta[0][k+1]-Delta[0][k-1])/dz/2.0),2.0);
//	  for (k = 1; k < kmax-1; k++)
//		  if(sin(Theta0*(Z[k]-Beta)) < 0) P[0][k] = (Gamma+1.0)/2.0*pow(0.75*(1-Z[k]*Z[k])*Delta[0][k]*cos(Theta0*Z[k]-Beta)+sin(Theta0*Z[k]-Beta)/Theta0
//					*(1.5*Z[k]*Delta[0][k]-(1-Z[k]*Z[k])*(Delta[0][k]-Delta[0][k-1])/dz),2);
//		  else P[0][k] = (Gamma+1.0)/2.0*pow(0.75*(1-Z[k]*Z[k])*Delta[0][k]*cos(Theta0*Z[k]-Beta)+sin(Theta0*Z[k]-Beta)/Theta0
//						*(1.5*Z[k]*Delta[0][k]-(1-Z[k]*Z[k])*(Delta[0][k+1]-Delta[0][k])/dz),2);
      if (is_mod_rel /*Modified relaxation*/) {
         AP[0] = 0; BP[0] = 0; Pd[0][0] = 0; Pd[0][kmax-1] = 0;
         for (k = 1; k < kmax-1; k++) {
            AP[k] = 1/(2+a0*dz*dz-AP[k-1]);
            BP[k] = (BP[k-1]-a0*dz*dz*(Prel[0][k]-P[0][k]))/(2+a0*dz*dz-AP[k-1]);
         }
         for (k = kmax-2; k > 0; k--) Pd[0][k] = Pd[0][k+1]*AP[k]+BP[k];
         for (k = 1; k < kmax-1; k++) P[0][k] = Prel[0][k] + relP*Pd[0][k];
      }
      else
         for (k = 1; k < kmax-1; k++) P[0][k] = relP*P[0][k] + (1-relP)*Prel[0][k];
      if (thomas_method) { //Progonka
      for (k = 1; k < kmax-1; k++) {
//		  H[0][0][k] = H[0][1][k];
          dPdz = (P[0][k+1]-P[0][k-1])/dz*0.5;
          for (j = 1; j < jmax-1; j++) Mu[j] = pow(fabs(T[0][j][k]), Omega-1); //viscocity coefficient
          A = (1.0-Z[k]*Z[k])/P[0][k];
          for (j = 1; j < jmax-1; j++) {
            //if (fabs(Aw) < 0.1) Aw = 0.1*Aw/fabs(Aw);
            E = (Gamma-1)/Gamma*0.5*T[0][j][k];
            if (is_cartesian) {
                Aw = A*(W[0][j][k]*XdZe-z0*U[0][j][k]*Z[k])*dy*dy/dz;
                V[0][j][k] = V[0][j-1][k]+dy*0.25*(Z[k]/P[0][k]*((W[0][j][k]+W[0][j-1][k])*XdZe-z0*Z[k]*XdZe*dZedX*(U[0][j][k]+U[0][j-1][k]))
                            -A*((W[0][j][k+1]-W[0][j][k-1]+W[0][j-1][k+1]-W[0][j-1][k-1])/dz*XdZe+0.5*z0*(U[0][j][k]+U[0][j-1][k])
                            -z0*Z[k]*XdZe*dZedX*(U[0][j][k+1]-U[0][j][k-1]+U[0][j-1][k+1]-U[0][j-1][k-1])/dz));
                Eu = z0*dy*dy*E/P[0][k]*(A*0.5*P[0][k]+Z[k]*Z[k]*XdZe*dZedX+Z[k]*XdZe*dZedX*A*dPdz);
                Ew = -dy*dy*E/P[0][k]*XdZe*(Z[k]+A*dPdz);
            }
            else {
                Aw = A*W[0][j][k]*dy*dy/dz;
                V[0][j][k] = V[0][j-1][k]+dy*0.25*(Z[k]/P[0][k]*(W[0][j][k]+W[0][j-1][k])
                            -A*((W[0][j][k+1]-W[0][j][k-1]+W[0][j-1][k+1]-W[0][j-1][k-1])/dz+2.5*Theta0*(U[0][j][k]+U[0][j-1][k])));
                Eu = Theta0*dy*dy*A*(W[0][j][k]*W[0][j][k]+E*0.5);
                Ew = -dy*dy*(A*Theta0*U[0][j][k]*W[0][j][k]+E/P[0][k]*(Z[k]+A*dPdz));
            }
            Eh = (Pr-1)/Pr*(U[0][j+1][k]*U[0][j+1][k]-2*U[0][j][k]*U[0][j][k]+U[0][j-1][k]*U[0][j-1][k]
                          + W[0][j+1][k]*W[0][j+1][k]-2*W[0][j][k]*W[0][j][k]+W[0][j-1][k]*W[0][j-1][k]);
            if(Omega == 1) {
                vdy2 = V[0][j][k]*dy*0.5;
                //1st order
                if (k == 1 || k == kmax-2 || is_first_order_theta || ( Beta == 0 && k == kmax/2 && kmax/2 != kmax/2.0)) {
                  AU[j] = (1-vdy2)/(fabs(Aw)+2-AU[j-1]*(1+vdy2));
                  BU[j] = (((Aw+fabs(Aw))*U[0][j][k-1]-(Aw-fabs(Aw))*U[0][j][k+1])*0.5+BU[j-1]*(1+vdy2)+Eu)
                          /(fabs(Aw)+2-AU[j-1]*(1+vdy2));
                  BW[j] = (((Aw+fabs(Aw))*W[0][j][k-1]-(Aw-fabs(Aw))*W[0][j][k+1])*0.5+BW[j-1]*(1+vdy2)+Ew)
                          /(fabs(Aw)+2-AU[j-1]*(1+vdy2));
                  AH[j] = (1/Pr-vdy2)/(fabs(Aw)+2/Pr-AH[j-1]*(1/Pr+vdy2));
                  BH[j] = (((Aw+fabs(Aw))*H[0][j][k-1]-(Aw-fabs(Aw))*H[0][j][k+1])*0.5+BH[j-1]*(1/Pr+vdy2)+Eh)
                          /(fabs(Aw)+2/Pr-AH[j-1]*(1/Pr+vdy2));
                }
                //central difference
//				if (k == 1 || k == kmax-2 || is_first_order_theta) {
//				  AU[j] = (1-vdy2)/(2-AU[j-1]*(1+vdy2));
//				  BU[j] = (-0.5*Aw*(U[0][j][k+1]-U[0][j][k-1])+BU[j-1]*(1+vdy2)+Eu)
//						  /(2-AU[j-1]*(1+vdy2));
//				  BW[j] = (-0.5*Aw*(W[0][j][k+1]-W[0][j][k-1])+BW[j-1]*(1+vdy2)+Ew)
//						  /(2-AU[j-1]*(1+vdy2));
//				  AH[j] = (1-vdy2*Pr)/(2-AH[j-1]*(1+vdy2*Pr));
//				  BH[j] = (-Pr*0.5*Aw*(H[0][j][k+1]-H[0][j][k-1])+BH[j-1]*(1+vdy2*Pr)+Eh*Pr)
//						  /(2-AH[j-1]*(1+vdy2*Pr));
//				}
                else {
                  AU[j] = (1-vdy2)/(fabs(Aw)*1.5+2-AU[j-1]*(1+vdy2));
                  BU[j] = (((Aw+fabs(Aw))*(4*U[0][j][k-1]-U[0][j][k-2])-(Aw-fabs(Aw))*(4*U[0][j][k+1]-U[0][j][k+2]))*0.25
                          +BU[j-1]*(1+vdy2)+Eu)/(fabs(Aw)*1.5+2-AU[j-1]*(1+vdy2));
                  BW[j] = (((Aw+fabs(Aw))*(4*W[0][j][k-1]-W[0][j][k-2])-(Aw-fabs(Aw))*(4*W[0][j][k+1]-W[0][j][k+2]))*0.25
                          +BW[j-1]*(1+vdy2)+Ew)/(fabs(Aw)*1.5+2-AU[j-1]*(1+vdy2));
                  AH[j] = (1/Pr-vdy2)/(fabs(Aw)*1.5-AH[j-1]*(1/Pr+vdy2)+2/Pr);
                  BH[j] = (((Aw+fabs(Aw))*(4*H[0][j][k-1]-H[0][j][k-2])-(Aw-fabs(Aw))*(4*H[0][j][k+1]-H[0][j][k+2]))*0.25
                          +BH[j-1]*(1/Pr+vdy2)+Eh)/(fabs(Aw)*1.5-AH[j-1]*(1/Pr+vdy2)+2/Pr);
                }
            }
            else {
               vdy2 = V[0][j][k]*dy*2;
               if (k == 1 || k == kmax-2 || is_first_order_theta) {
                  AU[j] = (Mu[j+1]+4*Mu[j]-Mu[j-1]-vdy2)
                          /(4*fabs(Aw)+8*Mu[j]+AU[j-1]*(Mu[j+1]-4*Mu[j]-Mu[j-1]-vdy2));
                  BU[j] = (((Aw+fabs(Aw))*U[0][j][k-1]-(Aw-fabs(Aw))*U[0][j][k+1])*2-BU[j-1]*(Mu[j+1]-4*Mu[j]-Mu[j-1]-vdy2)+4*Eu)
                          /(4*fabs(Aw)+8*Mu[j]+AU[j-1]*(Mu[j+1]-4*Mu[j]-Mu[j-1]-vdy2));
                  BW[j] = (((Aw+fabs(Aw))*W[0][j][k-1]-(Aw-fabs(Aw))*W[0][j][k+1])*2-BW[j-1]*(Mu[j+1]-4*Mu[j]-Mu[j-1]-vdy2)+4*Ew)
                          /(4*fabs(Aw)+8*Mu[j]+AU[j-1]*(Mu[j+1]-4*Mu[j]-Mu[j-1]-vdy2));
                  AH[j] = (Mu[j+1]+4*Mu[j]-Mu[j-1]-vdy2*Pr)
                          /(4*fabs(Aw)*Pr+8*Mu[j]+AH[j-1]*(Mu[j+1]-4*Mu[j]-Mu[j-1]-vdy2*Pr));
                  BH[j] = (((Aw+fabs(Aw))*H[0][j][k-1]-(Aw-fabs(Aw))*H[0][j][k+1])*2*Pr-BH[j-1]*(Mu[j+1]-4*Mu[j]-Mu[j-1]-vdy2*Pr)+4*Eh*Pr)
                          /(4*fabs(Aw)*Pr+8*Mu[j]+AH[j-1]*(Mu[j+1]-4*Mu[j]-Mu[j-1]-vdy2*Pr));
               }
               else {
                  AU[j] = (Mu[j+1]+4*Mu[j]-Mu[j-1]-vdy2)
                          /(6*fabs(Aw)+8*Mu[j]+AU[j-1]*(Mu[j+1]-4*Mu[j]-Mu[j-1]-vdy2));
                  BU[j] = (((Aw+fabs(Aw))*(4*U[0][j][k-1]-U[0][j][k-2])-(Aw-fabs(Aw))*(4*U[0][j][k+1]-U[0][j][k+2]))
                          -BU[j-1]*(Mu[j+1]-4*Mu[j]-Mu[j-1]-vdy2)+4*Eu)
                          /(6*fabs(Aw)+8*Mu[j]+AU[j-1]*(Mu[j+1]-4*Mu[j]-Mu[j-1]-vdy2));
                  BW[j] = (((Aw+fabs(Aw))*(4*W[0][j][k-1]-W[0][j][k-2])-(Aw-fabs(Aw))*(4*W[0][j][k+1]-W[0][j][k+2]))
                          -BW[j-1]*(Mu[j+1]-4*Mu[j]-Mu[j-1]-vdy2)+4*Ew)
                          /(6*fabs(Aw)+8*Mu[j]+AU[j-1]*(Mu[j+1]-4*Mu[j]-Mu[j-1]-vdy2));
                  AH[j] = (Mu[j+1]+4*Mu[j]-Mu[j-1]-vdy2*Pr)
                          /(6*fabs(Aw)*Pr+8*Mu[j]+AH[j-1]*(Mu[j+1]-4*Mu[j]-Mu[j-1]-vdy2*Pr));
                  BH[j] = (((Aw+fabs(Aw))*(4*H[0][j][k-1]-H[0][j][k-2])-(Aw-fabs(Aw))*(4*H[0][j][k+1]-H[0][j][k+2]))*Pr
                          -BH[j-1]*(Mu[j+1]-4*Mu[j]-Mu[j-1]-vdy2*Pr)+4*Eh*Pr)
                          /(6*fabs(Aw)*Pr+8*Mu[j]+AH[j-1]*(Mu[j+1]-4*Mu[j]-Mu[j-1]-vdy2*Pr));
               }
            }
         }
         for (j = jmax-2; j > 0; j--) {
            Us[0][j][k] = Us[0][j+1][k]*AU[j] + BU[j];
            Ws[0][j][k] = Ws[0][j+1][k]*AU[j] + BW[j];
            Hs[0][j][k] = Hs[0][j+1][k]*AH[j] + BH[j];
         }
      }
      //relaxation
      for (k = 1; k < kmax-1; k++)
         for (j = 1; j < jmax-1; j++) {
            U[0][j][k] = relU*Us[0][j][k] + (1-relU)*U[0][j][k];
            W[0][j][k] = relU*Ws[0][j][k] + (1-relU)*W[0][j][k];
            H[0][j][k] = relU*Hs[0][j][k] + (1-relU)*H[0][j][k];
         }
      }
      else { //Ustanovlenie
          for (k = 1; k < kmax-1; k++) {
              dPdz = (P[0][k+1]-P[0][k-1])/dz*0.5;
              for (j = 1; j < jmax-1; j++) Mu[j] = pow(fabs(T[0][j][k]), Omega-1); //viscocity coefficient
              for (j = 1; j < jmax-1; j++) {
                  A = (1.0-Z[k]*Z[k])/P[0][k];
                  Aw = A*Ws[0][j][k];
                  //if (fabs(A) < 0.0000001) A = 0;
                  V[0][j][k] = V[0][j-1][k]+dy*0.25*(Z[k]/P[0][k]*(Ws[0][j][k-1]+Ws[0][j-1][k-1]+Ws[0][j][k+1]+Ws[0][j-1][k+1]/*Ws[0][j][k]+Ws[0][j-1][k]*/)/2
                        -A*((Ws[0][j][k+1]-Ws[0][j][k-1]+Ws[0][j-1][k+1]-Ws[0][j-1][k-1])/dz
                        +2.5*Theta0*(Us[0][j][k]+Us[0][j-1][k])));
                  E = (Gamma-1)/(2*Gamma)*T[0][j][k];
                  Eu = Theta0*A*(Ws[0][j][k]*Ws[0][j][k]+E*0.5);
                  //Ew = -(Theta0*Aw*Us[0][j][k]+E/P[0][k]*(Z[k]+A*dPdz));
                  Ew = -(Theta0*Aw*Us[0][j][k]+E/P[0][k]*(Z[k]+A*dPdz));
                  Eh = (Pr-1)/Pr*(Us[0][j+1][k]*Us[0][j+1][k]-2*Us[0][j][k]*Us[0][j][k]+Us[0][j-1][k]*Us[0][j-1][k]
                                 +Ws[0][j+1][k]*Ws[0][j+1][k]-2*Ws[0][j][k]*Ws[0][j][k]+Ws[0][j-1][k]*Ws[0][j-1][k])/(dy*dy);
                  vdy2 = V[0][j][k];
                  //Crank-Nicolson
                  if (k == 1) {
                      U[0][j][k] = Us[0][j][k] - dt*( Aw*(Us[0][j][k]-Us[0][j][k-1])/dz
                          + vdy2*(Us[0][j+1][k]-Us[0][j-1][k]+Us[0][j+1][k-1]-Us[0][j-1][k-1])/dy*0.25
                          -(Mu[j+1]-Mu[j-1])*(Us[0][j+1][k]-Us[0][j-1][k]+Us[0][j+1][k-1]-Us[0][j-1][k-1])/dy/dy*0.125
                          - Mu[j]*(Us[0][j+1][k]-2*Us[0][j][k]+Us[0][j-1][k]
                                  +Us[0][j+1][k-1]-2*Us[0][j][k-1]+Us[0][j-1][k-1])/(dy*dy)*0.5 - Eu);
                      W[0][j][k] = Ws[0][j][k] - dt*( Aw*(Ws[0][j][k]-Ws[0][j][k-1])/dz
                          + vdy2*(Ws[0][j+1][k]-Ws[0][j-1][k]+Ws[0][j+1][k-1]-Ws[0][j-1][k-1])/dy*0.25
                          -(Mu[j+1]-Mu[j-1])*(Ws[0][j+1][k]-Ws[0][j-1][k]+Ws[0][j+1][k-1]-Ws[0][j-1][k-1])/dy/dy*0.125
                          - Mu[j]*(Ws[0][j+1][k]-2*Ws[0][j][k]+Ws[0][j-1][k]
                                  +Ws[0][j+1][k-1]-2*Ws[0][j][k-1]+Ws[0][j-1][k-1])/(dy*dy)*0.5
                          +(Theta0*Aw*Us[0][j][k]+E*(Z[k]+A*(P[0][k]-P[0][k-1])/dz)/P[0][k]));
                      H[0][j][k] = Hs[0][j][k] - dt*( Aw*(Hs[0][j][k]-Hs[0][j][k-1])/dz
                          + vdy2*(Hs[0][j+1][k]-Hs[0][j-1][k]+Hs[0][j+1][k-1]-Hs[0][j-1][k-1])/dy*0.25
                          -(Mu[j+1]-Mu[j-1])/Pr*(Hs[0][j+1][k]-Hs[0][j-1][k]+Hs[0][j+1][k-1]-Hs[0][j-1][k-1])/dy/dy*0.125
                          - Mu[j]*(Hs[0][j+1][k]-2*Hs[0][j][k]+Hs[0][j-1][k]
                                  +Hs[0][j+1][k-1]-2*Hs[0][j][k-1]+Hs[0][j-1][k-1])/dy/dy/Pr*0.5
                          -(Pr-1)/Pr*(Us[0][j+1][k]*Us[0][j+1][k]-2*Us[0][j][k]*Us[0][j][k]+Us[0][j-1][k]*Us[0][j-1][k]
                                     +Ws[0][j+1][k]*Ws[0][j+1][k]-2*Ws[0][j][k]*Ws[0][j][k]+Ws[0][j-1][k]*Ws[0][j-1][k]
                                     +Us[0][j+1][k-1]*Us[0][j+1][k-1]-2*Us[0][j][k-1]*Us[0][j][k-1]+Us[0][j-1][k-1]*Us[0][j-1][k-1]
                                     +Ws[0][j+1][k-1]*Ws[0][j+1][k-1]-2*Ws[0][j][k-1]*Ws[0][j][k-1]+Ws[0][j-1][k-1]*Ws[0][j-1][k-1])/(dy*dy)/2);
                  }
                  else if (k == kmax-2) {
                      U[0][j][k] = Us[0][j][k] - dt*((Aw+fabs(Aw))/2*(Us[0][j][k]-Us[0][j][k-1])/dz-(Aw-fabs(Aw))/2*(Us[0][j][k]-Us[0][j][k+1])/dz
                                      + vdy2*(Us[0][j+1][k]-Us[0][j-1][k])/dy*0.5
                                      -(Mu[j+1]-Mu[j-1])*(Us[0][j+1][k]-Us[0][j-1][k])/dy/dy*0.25
                                      - Mu[j]*(Us[0][j+1][k]-2*Us[0][j][k]+Us[0][j-1][k])/(dy*dy) - Eu);
                          W[0][j][k] = Ws[0][j][k] - dt*((Aw+fabs(Aw))/2*(Ws[0][j][k]-Ws[0][j][k-1])/dz-(Aw-fabs(Aw))/2*(Ws[0][j][k]-Ws[0][j][k+1])/dz
                                      + vdy2*(Ws[0][j+1][k]-Ws[0][j-1][k])/dy*0.5
                                      -(Mu[j+1]-Mu[j-1])*(Ws[0][j+1][k]-Ws[0][j-1][k])/dy/dy*0.25
                                      - Mu[j]*(Ws[0][j+1][k]-2*Ws[0][j][k]+Ws[0][j-1][k])/(dy*dy) - Ew);
                          H[0][j][k] = Hs[0][j][k] - dt*((Aw+fabs(Aw))/2*(Hs[0][j][k]-Hs[0][j][k-1])/dz-(Aw-fabs(Aw))/2*(Hs[0][j][k]-Hs[0][j][k+1])/dz
                                      + vdy2*(Hs[0][j+1][k]-Hs[0][j-1][k])/dy*0.5
                                      -(Mu[j+1]-Mu[j-1])/Pr*(Hs[0][j+1][k]-Hs[0][j-1][k])/dy/dy*0.25
                                      - Mu[j]*(Hs[0][j+1][k]-2*Hs[0][j][k]+Hs[0][j-1][k])/dy/dy/Pr - Eh);
          //			   U[0][j][k] = Us[0][j][k] - dt*(-Aw/2*(Us[0][j][k]-Us[0][j][k+1])/dz
          //											 + vdy2*(Us[0][j+1][k]-Us[0][j-1][k])/dy*0.5
          //											 -(Mu[j+1]-Mu[j-1])*(Us[0][j+1][k]-Us[0][j-1][k])/dy/dy*0.25
          //											 - Mu[j]*(Us[0][j+1][k]-2*Us[0][j][k]+Us[0][j-1][k])/(dy*dy) - Eu);
          //			   W[0][j][k] = Ws[0][j][k] - dt*(-Aw/2*(Ws[0][j][k]-Ws[0][j][k+1])/dz
          //											 + vdy2*(Ws[0][j+1][k]-Ws[0][j-1][k])/dy*0.5
          //											 -(Mu[j+1]-Mu[j-1])*(Ws[0][j+1][k]-Ws[0][j-1][k])/dy/dy*0.25
          //											 - Mu[j]*(Ws[0][j+1][k]-2*Ws[0][j][k]+Ws[0][j-1][k])/(dy*dy) - Ew);
          //			   H[0][j][k] = Hs[0][j][k] - dt*(-Aw/2*(Hs[0][j][k]-Hs[0][j][k+1])/dz
          //											 + vdy2*(Hs[0][j+1][k]-Hs[0][j-1][k])/dy*0.5
          //											 -(Mu[j+1]-Mu[j-1])/Pr*(Hs[0][j+1][k]-Hs[0][j-1][k])/dy/dy*0.25
          //											 - Mu[j]*(Hs[0][j+1][k]-2*Hs[0][j][k]+Hs[0][j-1][k])/dy/dy/Pr - Eh);
                      }
          //			else {
          //			   U[0][j][k] = Us[0][j][k] - dt*((Aw+fabs(Aw))/2*(3*Us[0][j][k]-4*Us[0][j][k-1]+Us[0][j][k-2])/2/dz
          //											 -(Aw-fabs(Aw))/2*(3*Us[0][j][k]-4*Us[0][j][k+1]+Us[0][j][k+2])/2/dz
          //											 + vdy2*(Us[0][j+1][k]-Us[0][j-1][k])/dy*0.5
          //											 -(Mu[j+1]-Mu[j-1])*(Us[0][j+1][k]-Us[0][j-1][k])/dy/dy*0.25
          //											 -Mu[j]*(Us[0][j+1][k]-2*Us[0][j][k]+Us[0][j-1][k])/(dy*dy) - Eu);
          //			   W[0][j][k] = Ws[0][j][k] - dt*((Aw+fabs(Aw))/2*(3*Ws[0][j][k]-4*Ws[0][j][k-1]+Ws[0][j][k-2])/2/dz
          //											 -(Aw-fabs(Aw))/2*(3*Ws[0][j][k]-4*Ws[0][j][k+1]+Ws[0][j][k+2])/2/dz
          //											 +vdy2*(Ws[0][j+1][k]-Ws[0][j-1][k])/dy*0.5
          //											 -(Mu[j+1]-Mu[j-1])*(Ws[0][j+1][k]-Ws[0][j-1][k])/dy/dy*0.25
          //											 -Mu[j]*(Ws[0][j+1][k]-2*Ws[0][j][k]+Ws[0][j-1][k])/(dy*dy) - Ew);
          //			   H[0][j][k] = Hs[0][j][k] - dt*((Aw+fabs(Aw))/2*(3*Hs[0][j][k]-4*Hs[0][j][k-1]+Hs[0][j][k-2])/2/dz
          //											 -(Aw-fabs(Aw))/2*(3*Hs[0][j][k]-4*Hs[0][j][k+1]+Hs[0][j][k+2])/2/dz
          //											 + vdy2*(Hs[0][j+1][k]-Hs[0][j-1][k])/dy*0.5
          //											 -(Mu[j+1]-Mu[j-1])/Pr*(Hs[0][j+1][k]-Hs[0][j-1][k])/dy/dy*0.25
          //											 -Mu[j]*(Hs[0][j+1][k]-2*Hs[0][j][k]+Hs[0][j-1][k])/dy/dy/Pr - Eh);
          //			}
          //			if (k == 1 || k == kmax-2 || is_first_order_theta) {
          //				U[0][j][k] = Us[0][j][k] - dt*((Aw+fabs(Aw))/2*(Us[0][j][k]-Us[0][j][k-1])/dz-(Aw-fabs(Aw))/2*(Us[0][j][k]-Us[0][j][k+1])/dz
          //											   +(vdy2+fabs(vdy2))*(Us[0][j][k]-Us[0][j-1][k])/dy*0.5
          //											   -(vdy2-fabs(vdy2))*(Us[0][j][k]-Us[0][j+1][k])/dy*0.5
          //											   -(Mu[j+1]-Mu[j-1])*(Us[0][j+1][k]-Us[0][j-1][k])/dy/dy*0.25
          //											   - Mu[j]*(Us[0][j+1][k]-2*Us[0][j][k]+Us[0][j-1][k])/(dy*dy) - Eu);
          //				W[0][j][k] = Ws[0][j][k] - dt*((Aw+fabs(Aw))/2*(Ws[0][j][k]-Ws[0][j][k-1])/dz-(Aw-fabs(Aw))/2*(Ws[0][j][k]-Ws[0][j][k+1])/dz
          //											   +(vdy2+fabs(vdy2))*(Ws[0][j][k]-Ws[0][j-1][k])/dy*0.5
          //											   -(vdy2-fabs(vdy2))*(Ws[0][j][k]-Ws[0][j+1][k])/dy*0.5
          //											   -(Mu[j+1]-Mu[j-1])*(Ws[0][j+1][k]-Ws[0][j-1][k])/dy/dy*0.25
          //											   - Mu[j]*(Ws[0][j+1][k]-2*Ws[0][j][k]+Ws[0][j-1][k])/(dy*dy) - Ew);
          //				H[0][j][k] = Hs[0][j][k] - dt*((Aw+fabs(Aw))/2*(Hs[0][j][k]-Hs[0][j][k-1])/dz-(Aw-fabs(Aw))/2*(Hs[0][j][k]-Hs[0][j][k+1])/dz
          //											   +(vdy2+fabs(vdy2))*(Hs[0][j][k]-Hs[0][j-1][k])/dy*0.5
          //											   -(vdy2-fabs(vdy2))*(Hs[0][j][k]-Hs[0][j+1][k])/dy*0.5
          //											   -(Mu[j+1]-Mu[j-1])/Pr*(Hs[0][j+1][k]-Hs[0][j-1][k])/dy/dy*0.25
          //											   - Mu[j]*(Hs[0][j+1][k]-2*Hs[0][j][k]+Hs[0][j-1][k])/dy/dy/Pr - Eh);
          //			}
                      else {
                          U[0][j][k] = Us[0][j][k] - dt*((Aw+fabs(Aw))/2*(3*Us[0][j][k]-4*Us[0][j][k-1]+Us[0][j][k-2])/2/dz
                                                        -(Aw-fabs(Aw))/2*(3*Us[0][j][k]-4*Us[0][j][k+1]+Us[0][j][k+2])/2/dz
                                                        +(vdy2+fabs(vdy2))*(Us[0][j][k]-Us[0][j-1][k])/dy*0.5
                                                        -(vdy2-fabs(vdy2))*(Us[0][j][k]-Us[0][j+1][k])/dy*0.5
                                                        -(Mu[j+1]-Mu[j-1])*(Us[0][j+1][k]-Us[0][j-1][k])/dy/dy*0.25
                                                        -Mu[j]*(Us[0][j+1][k]-2*Us[0][j][k]+Us[0][j-1][k])/(dy*dy) - Eu);
                          W[0][j][k] = Ws[0][j][k] - dt*((Aw+fabs(Aw))/2*(3*Ws[0][j][k]-4*Ws[0][j][k-1]+Ws[0][j][k-2])/2/dz
                                                        -(Aw-fabs(Aw))/2*(3*Ws[0][j][k]-4*Ws[0][j][k+1]+Ws[0][j][k+2])/2/dz
                                                        +(vdy2+fabs(vdy2))*(Ws[0][j][k]-Ws[0][j-1][k])/dy*0.5
                                                        -(vdy2-fabs(vdy2))*(Ws[0][j][k]-Ws[0][j+1][k])/dy*0.5
                                                        -(Mu[j+1]-Mu[j-1])*(Ws[0][j+1][k]-Ws[0][j-1][k])/dy/dy*0.25
                                                        -Mu[j]*(Ws[0][j+1][k]-2*Ws[0][j][k]+Ws[0][j-1][k])/(dy*dy) - Ew);
                          H[0][j][k] = Hs[0][j][k] - dt*((Aw+fabs(Aw))/2*(3*Hs[0][j][k]-4*Hs[0][j][k-1]+Hs[0][j][k-2])/2/dz
                                                        -(Aw-fabs(Aw))/2*(3*Hs[0][j][k]-4*Hs[0][j][k+1]+Hs[0][j][k+2])/2/dz
                                                        +(vdy2+fabs(vdy2))*(Hs[0][j][k]-Hs[0][j-1][k])/dy*0.5
                                                        -(vdy2-fabs(vdy2))*(Hs[0][j][k]-Hs[0][j+1][k])/dy*0.5
                                                        -(Mu[j+1]-Mu[j-1])/Pr*(Hs[0][j+1][k]-Hs[0][j-1][k])/dy/dy*0.25
                                                        -Mu[j]*(Hs[0][j+1][k]-2*Hs[0][j][k]+Hs[0][j-1][k])/dy/dy/Pr - Eh);
                      }
                   }
                }
      }
      //checking approximation
      for (k = 1; k < kmax-1; k++)
         if (fabs(1-Prel[0][k]/P[0][k]) > AccuracyPressure*relP) k = kmax;
      if(k == kmax-1) converged = true;
   }
   for (k = 0; k < kmax; k++)
       DeltaRel[0][k] = Delta[0][k];
   for (k = 0; k < kmax; k++)
       if( (U[0][jmax-1][k]-U[0][jmax-2][k])/dy > 0.005 ) {
           qDebug() << QString("Warning!!! Wrong results! (Theta0=%1,Beta=%2,Hw=%3,Vw=%4)")
                       .arg(Theta0d).arg(Betad).arg(Hw).arg(Vw);
           k = kmax;
           status = "bad";
       }
}
void FullWing(bool is_cartesian, bool is_interpolated) { //computing main field
    //set up initial field
    for (i = imax/L; i < imax; i++) {
        P[i][0] = P[int(imax/L)-1][0];
        P[i][kmax-1] = P[int(imax/L)-1][kmax-1];
        Delta[i][0] = Delta[int(imax/L)-1][0];
        Delta[i][kmax-1] = Delta[int(imax/L)-1][kmax-1];
        for (j = 0; j < jmax; j++) {
            U[i][j][0] = U[int(imax/L)-1][j][0];
            W[i][j][0] = W[int(imax/L)-1][j][0];
            H[i][j][0] = H[int(imax/L)-1][j][0];
            U[i][j][kmax-1] = U[int(imax/L)-1][j][kmax-1];
            W[i][j][kmax-1] = W[int(imax/L)-1][j][kmax-1];
            H[i][j][kmax-1] = H[int(imax/L)-1][j][kmax-1];
        }
    }
    for (i = 1; i < imax; i++)
        for (k = 1; k < kmax-1; k++) {
            if(!is_interpolated) P[i][k] = P[0][k] + (P[i][0]-P[0][0]);// - (1+cos((X[i]+0.5)*2*Pi))*(1+cos(Z[k]*Pi))/100;
            Delta[i][k] = Delta[0][k];
            for (j = 0; j < jmax; j++) {
                if(!is_interpolated) {
                    U[i][j][k] = U[0][j][k];
                    W[i][j][k] = W[0][j][k];
                    H[i][j][k] = H[0][j][k];
                }
                V[i][j][k] = V[0][j][k];
                T[i][j][k] = T[0][j][k];
                Us[i][j][k] = U[i][j][k];
                Ws[i][j][k] = W[i][j][k];
                Hs[i][j][k] = H[i][j][k];
            }
        }
    for (k = 0; k < kmax; k++) {
        if(!is_interpolated && (L == 1.0))
            P[imax-1][k] = P[0][k]+cos(Z[k]*10*Pi)/50*(1-fabs(Z[k])); //NOTE: условие на задней кромке
            //P[imax-1][k] = P[0][k]*(1-pow((1+cos(Z[k]*Pi))*0.5, 3.0));
        else
            P[imax-1][k] = P[int(imax/L)-1][k];
        P[imax][k] = P[imax-1][k];//3*P[imax-1][k]-3*P[imax-2][k]+P[imax-3][k];
        Prel[imax-1][k] = P[imax-1][k];
    }
    double Integral, Pstrong, A, Au, Aw, E, Eu, Ew, Eh, vdy, dPdx, dPdz, dDdx, dDdz, dUdx,
           AudUx, AwdUz, AudWx, AwdWz, AudHx, AwdHz, fAu, fAw;
    double k1 = sqrt((Gamma-1.0)/Gamma*0.5)*dy/3.0;
    dt = 0.000005; converged = false;
    for (NIterations = 0; !converged && NIterations < 1000/*MaxIterations*/; NIterations++) {
        QCoreApplication::processEvents(); //to prevent freezing of gui
        for (i = 1; i < imax-1; i++) memcpy(Prel[i], P[i], sizeof(double)*kmax);
        for (i = imax/L; i < imax; i++)
            for (k = 0; k < kmax; k++) {
                U[i][0][k] = U[i][1][k];
                W[i][0][k] = W[i][1][k];
                H[i][0][k] = H[i][1][k];
            }
        for (i = 1; i < imax; i++)
            for (j = 0; j < jmax-1; j++)
                for (k = 0; k < kmax; k++)
                    T[i][j][k] = H[i][j][k]-U[i][j][k]*U[i][j][k]-W[i][j][k]*W[i][j][k];
        for(i = 1; i < imax; i++)
            for(k = 0; k < kmax; k++) {
                Integral = T[i][0][k];
                for (j = 1; j < jmax-1; j+=2) Integral +=  4*T[i][j][k] + 2*T[i][j+1][k];
                Delta[i][k] = k1/P[i][k]*Integral;
            }
        for (i = 1; i < imax-1; i++)
            for (k = 1; k < kmax-1; k++) {
                if (i > 2) dDdx = (3*Delta[i][k]-4*Delta[i-1][k]+Delta[i-2][k]+3*Body[i][k]-4*Body[i-1][k]+Body[i-2][k])/dx*0.5;
                else dDdx = (Delta[i][k]-Delta[i-1][k]+Body[i][k]-Body[i-1][k])/dx;
                dDdz = (Delta[i][k+1]-Delta[i][k-1]+Body[i][k+1]-Body[i][k-1])/dz*0.5;
                if (is_cartesian) {
                    Pstrong = 0.75*(1-Z[k]*Z[k])*(Body[i][k]+Delta[i][k])+(1-Z[k]*Z[k])*X[i]*dDdx
                            -Z[k]*X[i]/Ze[i]*(Ze[i+1]-Ze[i])/dx*(dDdz*(1-Z[k]*Z[k])-1.5*Z[k]*(Delta[i][k]+Body[i][k]));
                    //P[i][k] = sqrt(X[i]*(1-Z[k]*Z[k]))/(Gamma*khi*khi)+(Gamma+1)*0.25*Pstrong*Pstrong
                        //	+Pstrong*sqrt(sqrt(X[i]*(1-Z[k]*Z[k]))/(khi*khi)+(Gamma+1)*0.25*Pstrong*(Gamma+1)*0.25*Pstrong);
                    P[i][k] = (Gamma+1.0)/2.0*pow(Pstrong*cos(Beta)+sin(Beta)/z0*((1-Z[k]*Z[k])*dDdz-1.5*Z[k]*Delta[i][k]),2.0);
                }
                else P[i][k] = (Gamma+1.0)/2.0*pow((1-Z[k]*Z[k])*(0.75*Delta[i][k]+X[i]*dDdx)
                             * cos(Theta0*Z[k]-Beta)+sin(Theta0*Z[k]-Beta)/Theta0*(1.5*Z[k]*Delta[i][k]-(1-Z[k]*Z[k])*dDdz),2);
            }
        if (is_mod_rel /*Modified relaxation*/) {
            for (i = 1; i < imax-1; i++) {
                AP[0] = 0; BP[0] = 0; Pd[i][0] = 0; Pd[i][kmax-1] = 0;
                for (k = 1; k < kmax-1; k++) {
                    AP[k] = 1/(2+a0*dz*dz-AP[k-1]);
                    BP[k] = (BP[k-1]-a0*dz*dz*(Prel[i][k]-P[i][k]))/(2+a0*dz*dz-AP[k-1]);
                }
                for (k = kmax-2; k > 0; k--) Pd[i][k] = Pd[i][k+1]*AP[k]+BP[k];
                for (k = 1; k < kmax-1; k++) P[i][k] = Prel[i][k] + relP*Pd[i][k];
            }
        }
        else
            for (i = 1; i < imax-1; i++)
                for (k = 1; k < kmax-1; k++) P[i][k] = relP2*P[i][k] + (1-relP2)*Prel[i][k];
//        for (k = 0; k < kmax; k++) {
//            P[imax-1][k] = P[imax-2][k];
//            P[imax][k] = 3*P[imax-1][k]-3*P[imax-2][k]+P[imax-3][k];
//        }
        /*
        //Modification of relaxation
        if (is_mod_rel && 0)
        {
            double A = -2/(dx*dx)-2/(dz*dz) - a0;
            double ax = 1/(dx*dx)/A;
            double az = 1/(dz*dz)/A;
            int MH = 2*(imax-2)+1; //matrix height
            int MW = (imax-2)*(kmax-2); //matrix width
            double Matrix[MH][MW];

            for (int i = 0; i < imax; i++)
                for (int k = 0; k < kmax; k++)
                    Pd[i][k] = 0;
            for (int i = 1; i < imax-1; i++)
                for (int k = 1; k < kmax-1; k++)
                    Pd[i][k] = a0/A*(Prel[i][k]-P[i][k]); // <- 1/A

            //ax = 0.5; az = 0.2;
            //Matrix filling
            for (int k = 0; k < MW; k++)
                for (int j = 0; j < MH; j++) {
                    if (j == imax-2) Matrix[j][k] = 1;
                    else if ( (j == 0 && k < (imax-2)*(kmax-3)) || (j == 2*(imax-2) && k >= imax-2)) Matrix[j][k] = az;
                    else if ((j == imax-1 && k%(imax-2) != 0) || (j == imax-3 && (k-(imax-3))%(imax-2) != 0)) Matrix[j][k] = ax;
                    else Matrix[j][k] = 0;
                }
            //Lower matrix
            double a1, a2;
            for (int k = 0; k < MW; k++)
                for (int j = 1; j <= imax-2 && (k+j) < MW; j++)
                    if (Matrix[imax-2+j][k+j] != 0) {
                        a1 = Matrix[imax-2][k];
                        a2 = Matrix[imax-2+j][k+j];
                        for (int i = 0; i <= imax-2; i++) Matrix[imax-2+j-i][k+j] -= Matrix[imax-2-i][k]/a1*a2;
                        Pd[(k+j)/(kmax-2)+1][(k+j)%(kmax-2)+1] -= Pd[k/(kmax-2)+1][k%(kmax-2)+1]/a1*a2;
                        //b[k+j] -= b[k]/a1*a2;
                    }
            //1 on the diagonal
            for (int k = MW-1; k > 0; k--) {
                Pd[k/(kmax-2)+1][k%(kmax-2)+1] /= Matrix[imax-2][k];
                //b[k] /= Matrix[imax-2][k];
                for (int j = 0; j <= imax-2 && k-j > 0; j++) Matrix[imax-2-j][k] /= Matrix[imax-2][k];
            }
            //Upper matrix
            for (int k = MW-1; k >= 0; k--)
                for (int j = 1; j <= imax-2 && k-j >= 0; j++)
                    if (Matrix[imax-2-j][k-j] != 0) {
                        Matrix[imax-2-j][k-j] -= Matrix[imax-2][k]*Matrix[imax-2-j][k-j];
                        //b[k-j] -= b[k]*Matrix[imax-2-j][k-j];
                        Pd[(k-j)/(kmax-2)+1][(k-j)%(kmax-2)+1] -= Pd[k/(kmax-2)+1][k%(kmax-2)+1]*Matrix[imax-2-j][k-j];
                    }

            for (int i = 1; i < imax-1; i++)
                for (int k = 1; k < kmax-1; k++)
                    P[i][k] = Prel[i][k] + relP*Pd[i][k];
        }
        else
            for(i = 1; i < imax-1; i++)
               for (k = 1; k < kmax-1; k++)
                  P[i][k] = relP*P[i][k] + (1-relP)*Prel[i][k];
        */
        for (i = imax/L; i < imax; i++) {
            for (j = 0; j < jmax-1; j++) {
                U[i][j][0] = U[i][j][1];
                W[i][j][0] = W[i][j][1];
                H[i][j][0] = H[i][j][1];
                U[i][j][kmax-1] = U[i][j][kmax-2];
                W[i][j][kmax-1] = W[i][j][kmax-2];
                H[i][j][kmax-1] = H[i][j][kmax-2];
            }
            P[i][0] = P[i][1];
            P[i][kmax-1] = P[i][kmax-2];
        }
    //Progonka
    for(i = 1; i < imax; i++) {
        QCoreApplication::processEvents(); //to prevent freezing of gui
        for (k = 1; k < kmax-1; k++) {
            //for (j = 1; j < jmax-1; j++) Mu[j] = pow(fabs(T[i][j][k]), Omega-1); //viscocity coefficient
            A = (1.0-Z[k]*Z[k])/P[i][k];            
            if (i == imax-1)
                dPdx = (P[i][k]-P[i-1][k])/dx;
            else
                dPdx = (P[i+1][k]-P[i-1][k])/dx*0.5;
            dPdz = (P[i][k+1]-P[i][k-1])/dz*0.5;
            if (L > 1.0) {
                if (i < imax/L) {
                    AU[0] = 0; BU[0] = 0;
                    BW[0] = 0;
                    AH[0] = 0; BH[0] = Hw;
                }
                else {
                    AU[0] = 1; BU[0] = 0;
                    BW[0] = 0;
                    AH[0] = 1; BH[0] = 0;
                }
            }
            for (j = 1; j < jmax-1; j++) {
                E = (Gamma-1)/Gamma*0.5*T[i][j][k];
                if (i != imax-1) dUdx = (U[i+1][j][k]-U[i-1][j][k]+U[i+1][j-1][k]-U[i-1][j-1][k])/dx*0.25;
                else dUdx = (U[i][j][k]-U[i-1][j][k]+U[i][j-1][k]-U[i-1][j-1][k])/dx*0.5;
                if (is_cartesian) {
                    Au = A*z0*X[i]*U[i][j][k]*dy*dy/dx;
                    Aw = A*(W[i][j][k]*X[i]/Ze[i]-z0*U[i][j][k]*Z[k]*X[i]/Ze[i]*(Ze[i+1]-Ze[i])/dx)*dy*dy/dz;
                    V[i][j][k] = V[i][j-1][k]+dy*0.25*(Z[k]/P[i][k]*((W[i][j][k]+W[i][j-1][k])*X[i]/Ze[i]-z0*X[i]/Ze[i]*(Ze[i+1]-Ze[i])/dx*Z[k]*(U[i][j][k]+U[i][j-1][k]))
                        -A*((W[i][j][k+1]-W[i][j][k-1]+W[i][j-1][k+1]-W[i][j-1][k-1])/dz*X[i]/Ze[i]+0.5*z0*(U[i][j][k]+U[i][j-1][k])
                        -z0*Z[k]*X[i]/Ze[i]*(Ze[i+1]-Ze[i])/dx*(U[i][j][k+1]-U[i][j][k-1]+U[i][j-1][k+1]-U[i][j-1][k-1])/dz+z0*X[i]*dUdx*4));
                    Eu = z0*dy*dy*E/P[i][k]*(A*(P[i][k]*0.5-X[i]*dPdx)+Z[k]*X[i]/Ze[i]*(Ze[i+1]-Ze[i])/dx*(Z[k]+A*dPdz));
                    Ew = -dy*dy*E/P[i][k]*X[i]/Ze[i]*(Z[k]+A*dPdz);
                }
                else {
                    Au = A*Theta0*X[i]*U[i][j][k]*dy*dy/dx;
                    Aw = A*W[i][j][k]*dy*dy/dz;
                    V[i][j][k] = V[i][j-1][k]+dy*0.25*(Z[k]/P[i][k]*(W[i][j][k]+W[i][j-1][k])-A*(Theta0*X[i]*dUdx*4
                        +(W[i][j][k+1]-W[i][j][k-1]+W[i][j-1][k+1]-W[i][j-1][k-1])/dz+2.5*Theta0*(U[i][j][k]+U[i][j-1][k])));
                    Eu = Theta0*dy*dy*A*(W[i][j][k]*W[i][j][k]+E*(0.5-X[i]/P[i][k]*dPdx));
                    Ew = -dy*dy*(A*Theta0*U[i][j][k]*W[i][j][k]+E*(Z[k]/P[i][k]+A/P[i][k]*dPdz));
                }
                Eh = (Pr-1)/Pr*(U[i][j+1][k]*U[i][j+1][k]-2*U[i][j][k]*U[i][j][k]+U[i][j-1][k]*U[i][j-1][k]
                              + W[i][j+1][k]*W[i][j+1][k]-2*W[i][j][k]*W[i][j][k]+W[i][j-1][k]*W[i][j-1][k]);
                vdy = V[i][j][k]*dy;
                if ( k == 1 || k == kmax-2 || is_first_order_theta ) {
                    AwdUz = ((Aw+fabs(Aw))*U[i][j][k-1]-(Aw-fabs(Aw))*U[i][j][k+1]);
                    AwdWz = ((Aw+fabs(Aw))*W[i][j][k-1]-(Aw-fabs(Aw))*W[i][j][k+1]);
                    AwdHz = ((Aw+fabs(Aw))*H[i][j][k-1]-(Aw-fabs(Aw))*H[i][j][k+1]);
                    fAw = 2*fabs(Aw);
                }
                else {
                    AwdUz = ((Aw+fabs(Aw))*(4*U[i][j][k-1]-U[i][j][k-2])-(Aw-fabs(Aw))*(4*U[i][j][k+1]-U[i][j][k+2]))*0.5;
                    AwdWz = ((Aw+fabs(Aw))*(4*W[i][j][k-1]-W[i][j][k-2])-(Aw-fabs(Aw))*(4*W[i][j][k+1]-W[i][j][k+2]))*0.5;
                    AwdHz = ((Aw+fabs(Aw))*(4*H[i][j][k-1]-H[i][j][k-2])-(Aw-fabs(Aw))*(4*H[i][j][k+1]-H[i][j][k+2]))*0.5;
                    fAw = 3*fabs(Aw);
                }
                if ( i == 1) {
                    AudUx = ((Au+fabs(Au))*U[i-1][j][k]-(Au-fabs(Au))*U[i+1][j][k]);
                    AudWx = ((Au+fabs(Au))*W[i-1][j][k]-(Au-fabs(Au))*W[i+1][j][k]);
                    AudHx = ((Au+fabs(Au))*H[i-1][j][k]-(Au-fabs(Au))*H[i+1][j][k]);
                    fAu = 2*fabs(Au);
                }
                else if (i == imax-1 || i == imax-2) {
                    AudUx = ((Au+fabs(Au))*(4*U[i-1][j][k]-U[i-2][j][k]))*0.5;
                    AudWx = ((Au+fabs(Au))*(4*W[i-1][j][k]-W[i-2][j][k]))*0.5;
                    AudHx = ((Au+fabs(Au))*(4*H[i-1][j][k]-H[i-2][j][k]))*0.5;
                    fAu = 3*fabs(Au);
                }
                else {
                    AudUx = ((Au+fabs(Au))*(4*U[i-1][j][k]-U[i-2][j][k])-(Au-fabs(Au))*(4*U[i+1][j][k]-U[i+2][j][k]))*0.5;
                    AudWx = ((Au+fabs(Au))*(4*W[i-1][j][k]-W[i-2][j][k])-(Au-fabs(Au))*(4*W[i+1][j][k]-W[i+2][j][k]))*0.5;
                    AudHx = ((Au+fabs(Au))*(4*H[i-1][j][k]-H[i-2][j][k])-(Au-fabs(Au))*(4*H[i+1][j][k]-H[i+2][j][k]))*0.5;
                    fAu = 3*fabs(Au);
                }
                AU[j] = (2-vdy)/(fAu+fAw+4-AU[j-1]*(2+vdy));
                AH[j] = (2/Pr-vdy)/(fAu+fAw+4/Pr-AH[j-1]*(2/Pr+vdy));
                BU[j] = (AudUx+AwdUz+BU[j-1]*(2+vdy)+2*Eu)/(fAu+fAw+4-AU[j-1]*(2+vdy));
                BW[j] = (AudWx+AwdWz+BW[j-1]*(2+vdy)+2*Ew)/(fAu+fAw+4-AU[j-1]*(2+vdy));
                BH[j] = (AudHx+AwdHz+BH[j-1]*(2/Pr+vdy)+2*Eh)/(fAu+fAw+4/Pr-AH[j-1]*(2/Pr+vdy));
            }
            for (j = jmax-2; j > 0; j--) {
                Us[i][j][k] = Us[i][j+1][k]*AU[j] + BU[j];
                Ws[i][j][k] = Ws[i][j+1][k]*AU[j] + BW[j];
                Hs[i][j][k] = Hs[i][j+1][k]*AH[j] + BH[j];
            }
        }
    }
    //relaxation
    for (i = 1; i < imax; i++)
        for (k = 1; k < kmax-1; k++)
            for (j = 1; j < jmax-1; j++) {
               U[i][j][k] = relU*Us[i][j][k] + (1-relU)*U[i][j][k];
               W[i][j][k] = relU*Ws[i][j][k] + (1-relU)*W[i][j][k];
               H[i][j][k] = relU*Hs[i][j][k] + (1-relU)*H[i][j][k];
            }
      //Ustanovlenie
      /*for (i = 1; i < imax-1; i++)
        {
            for (k = 1; k < kmax-1; k++)
                for (j = 1; j < jmax-1; j++) {
                    A = (1-Z[k]*Z[k])/P[i][k];
                    AU = A*Theta0*X[i]*Us[i][j][k];
                    AW = A*Ws[i][j][k];
                    //if (fabs(A) < 0.0000001) A = 0;
                    if (!is_first_order_v)
                        V[i][j][k] = V[i][j-1][k]+dy*0.25*(Z[k]/P[i][k]*(Ws[i][j][k]+Ws[i][j-1][k])
                            -A*(Theta0*X[i]*(Us[i+1][j][k]-Us[i-1][j][k]+Us[i+1][j-1][k]-Us[i-1][j-1][k])/dx
                            +(Ws[i][j][k+1]-Ws[i][j][k-1]+Ws[i][j-1][k+1]-Ws[i][j-1][k-1])/dz
                            +2.5*Theta0*(Us[i][j][k]+Us[i][j-1][k])));
                    if (is_first_order_v)
                        V[i][j][k] = V[i][j-1][k]+dy*(Z[k]/P[i][k]*Ws[i][j-1][k]/2
                            -A*(Theta0*X[i]*(Us[i+1][j-1][k]-Us[i-1][j-1][k])/(2*dx)
                            +(Ws[i][j-1][k+1]-Ws[i][j-1][k-1])/(2*dz)+1.25*Theta0*Us[i][j-1][k]));
                        E = (Gamma-1)/(2*Gamma)*T[i][j][k];
    //					if(i == imax-1)
    //						E1 = Theta0*A*(Ws[i][j][k]*Ws[i][j][k]+E*(0.5-X[i]/P[i][k]*(P[i][k]-P[i-1][k])/dx));
    //					else
                    E1 = Theta0*A*(Ws[i][j][k]*Ws[i][j][k]+E*(0.5-X[i]/P[i][k]*(P[i+1][k]-P[i-1][k])/(2*dx)));
                    E2 = -(Theta0*A*Us[i][j][k]*Ws[i][j][k]+E*(Z[k]/P[i][k]+A/P[i][k]*(P[i][k+1]-P[i][k-1])/(2*dz)));
                    E3 = (Pr-1)/Pr*( Us[i][j+1][k]*Us[i][j+1][k]-2*Us[i][j][k]*Us[i][j][k]+Us[i][j-1][k]*Us[i][j-1][k]
                    +Ws[i][j+1][k]*Ws[i][j+1][k]-2*Ws[i][j][k]*Ws[i][j][k]+Ws[i][j-1][k]*Ws[i][j-1][k])/(dy*dy);

                    U[i][j][k] = Us[i][j][k] - dt*(
                        (AU+fabs(AU))/2*(Us[i][j][k]-Us[i-1][j][k])/dx +
                        (AU-fabs(AU))/2*(Us[i+1][j][k]-Us[i][j][k])/dx +
                        (AW+fabs(AW))/2*(Us[i][j][k]-Us[i][j][k-1])/dz +
                        (AW-fabs(AW))/2*(Us[i][j][k+1]-Us[i][j][k])/dz +
                        V[i][j][k]*(Us[i][j+1][k]-Us[i][j-1][k])/(2*dy) -
                        (Us[i][j+1][k]-2*Us[i][j][k]+Us[i][j-1][k])/(dy*dy) - E1);
                    W[i][j][k] = Ws[i][j][k] - dt*(
                        (AU+fabs(AU))/2*(Ws[i][j][k]-Ws[i-1][j][k])/dx +
                        (AU-fabs(AU))/2*(Ws[i+1][j][k]-Ws[i][j][k])/dx +
                        (AW+fabs(AW))/2*(Ws[i][j][k]-Ws[i][j][k-1])/dz +
                        (AW-fabs(AW))/2*(Ws[i][j][k+1]-Ws[i][j][k])/dz +
                        V[i][j][k]*(Ws[i][j+1][k]-Ws[i][j-1][k])/(2*dy) -
                        (Ws[i][j+1][k]-2*Ws[i][j][k]+Ws[i][j-1][k])/(dy*dy) - E2);
                    H[i][j][k] = Hs[i][j][k] - dt*(
                        (AU+fabs(AU))/2*(Hs[i][j][k]-Hs[i-1][j][k])/dz +
                        (AU-fabs(AU))/2*(Hs[i+1][j][k]-Hs[i][j][k])/dz +
                        (AW+fabs(AW))/2*(Hs[i][j][k]-Hs[i][j][k-1])/dz +
                        (AW-fabs(AW))/2*(Hs[i][j][k+1]-Hs[i][j][k])/dz +
                        V[i][j][k]*(Hs[i][j+1][k]-Hs[i][j-1][k])/(2*dy) -
                        (Hs[i][j+1][k]-2*Hs[i][j][k]+Hs[i][j-1][k])/(dy*dy)/Pr - E3);
            }
        }*/
    //checking approximation
//	for (k = 1; k < kmax-1; k++)
//	   if (fabs(1-Prel[imax/2][k]/P[imax/2][k]) > AccuracyPressure*relP) k = kmax;
//	if(k == kmax-1) converged = true;
    }
}
void FullWingWithTrail(bool is_cartesian) { //computing main field and trail
    //set up initial field
    for (i = 1; i < imax; i++)
        for (k = 0; k < kmax; k++) {
            P[i][k] = P[0][k];// - (1+cos((X[i]+0.5)*2*Pi))*(1+cos(Z[k]*Pi))/100;
            Delta[i][k] = Delta[0][k];
            for (j = 0; j < jmax; j++) {
                U[i][j][k] = U[0][j][k];
                W[i][j][k] = W[0][j][k];
                H[i][j][k] = H[0][j][k];
                V[i][j][k] = V[0][j][k];
                T[i][j][k] = T[0][j][k];
                Us[i][j][k] = U[0][j][k];
                Ws[i][j][k] = W[0][j][k];
                Hs[i][j][k] = H[0][j][k];
            }
        }
    for (k = 0; k < kmax; k++) {
        P[imax-1][k] = P[0][k];//*(1-0.1*pow((1+cos(Z[k]*Pi))*0.5, 3.0));
        P[imax][k] = 3*P[imax-1][k]-3*P[imax-2][k]+P[imax-3][k];
        Prel[imax-1][k] = P[imax-1][k];
    }
    double Integral, Pstrong, A, Au, Aw, E, Eu, Ew, Eh, vdy, dPdx, dPdz, dDdx, dDdz, dUdx,
           AudUx, AwdUz, AudWx, AwdWz, AudHx, AwdHz, fAu, fAw;
    double k1 = sqrt((Gamma-1.0)/Gamma*0.5)*dy/3.0;
    dt = 0.000005; converged = false;
    for (NIterations = 0; !converged && NIterations < MaxIterations; NIterations++) {
        QCoreApplication::processEvents(); //to prevent freezing of gui
        for (i = 1; i < imax; i++) memcpy(Prel[i], P[i], sizeof(double)*kmax);
        for (i = imax/L; i < imax; i++)
            for (k = 0; k < kmax; k++) {
                U[i][0][k] = U[i][1][k];
                W[i][0][k] = W[i][1][k];
                H[i][0][k] = H[i][1][k];
            }
        for (i = 1; i < imax; i++)
            for (j = 0; j < jmax-1; j++)
                for (k = 0; k < kmax; k++)
                    T[i][j][k] = H[i][j][k]-U[i][j][k]*U[i][j][k]-W[i][j][k]*W[i][j][k];
        for(i = 1; i < imax; i++)
            for(k = 0; k < kmax; k++) {
                Integral = T[i][0][k];
                for (j = 1; j < jmax-1; j+=2) Integral +=  4*T[i][j][k] + 2*T[i][j+1][k];
                Delta[i][k] = k1/P[i][k]*Integral;
            }
        for (i = 1; i < imax-1; i++)
            for (k = 0; k < kmax; k++) {
                if (i > 2 && 0) dDdx = (3*Delta[i][k]-4*Delta[i-1][k]+Delta[i-2][k])/dx*0.5;
                else dDdx = (Delta[i][k]-Delta[i-1][k])/dx;
                dDdz = (Delta[i][k+1]-Delta[i][k-1])/dz*0.5;
                if (is_cartesian) {
                    if(i < imax) Pstrong = 0.75*(1-Z[k]*Z[k])*Delta[i][k]+(1-Z[k]*Z[k])*X[i]*dDdx
                            -Z[k]*dDdz*(1-Z[k]*Z[k])+1.5*Z[k]*Z[k]*Delta[i][k];
                    else Pstrong = 0.75*(1-Z[k]*Z[k])*Delta[i][k]+(1-Z[k]*Z[k])*X[i]*dDdx
                            -0*Z[k]*dDdz*(1-Z[k]*Z[k])+0*1.5*Z[k]*Z[k]*Delta[i][k];
                    //P[i][k] = sqrt(X[i]*(1-Z[k]*Z[k]))/(Gamma*khi*khi)+(Gamma+1)*0.25*Pstrong*Pstrong
                        //	+Pstrong*sqrt(sqrt(X[i]*(1-Z[k]*Z[k]))/(khi*khi)+(Gamma+1)*0.25*Pstrong*(Gamma+1)*0.25*Pstrong);
                    P[i][k] = (Gamma+1.0)/2.0*pow(Pstrong*cos(Beta)+sin(Beta)/z0*((1-Z[k]*Z[k])*dDdz-1.5*Z[k]*Delta[i][k]),2.0);
                }
                else P[i][k] = (Gamma+1.0)/2.0*pow((1-Z[k]*Z[k])*(0.75*Delta[i][k]+X[i]*dDdx)
                             *cos(Theta0*Z[k]-Beta)+sin(Theta0*Z[k]-Beta)/Theta0*(1.5*Z[k]*Delta[i][k]-(1-Z[k]*Z[k])*dDdz),2);
            }
        if (is_mod_rel /*Modified relaxation*/) {
            for (i = 1; i < imax; i++) {
                AP[0] = 0; BP[0] = 0; Pd[i][0] = 0; Pd[i][kmax-1] = 0;
                for (k = 1; k < kmax-1; k++) {
                    AP[k] = 1/(2+a0*dz*dz-AP[k-1]);
                    BP[k] = (BP[k-1]-a0*dz*dz*(Prel[i][k]-P[i][k]))/(2+a0*dz*dz-AP[k-1]);
                }
                for (k = kmax-2; k > 0; k--) Pd[i][k] = Pd[i][k+1]*AP[k]+BP[k];
                for (k = 1; k < kmax-1; k++) P[i][k] = Prel[i][k] + relP*Pd[i][k];
            }
        }
        else
            for (i = 1; i < imax; i++)
                for (k = 1; k < kmax-1; k++) P[i][k] = relP*P[i][k] + (1-relP)*Prel[i][k];
        for (k = 1; k < kmax-1; k++) { P[imax-1][k] = P[imax-2][k]; P[imax][k] = 3*P[imax-1][k]-3*P[imax-2][k]+P[imax-3][k]; }
    //Progonka
    for(i = 1; i < imax; i++) {
        if (i >= imax/L) {
            for (j = 0; j < jmax-1; j++) {
                U[i][j][0] = U[i][j][1];
                W[i][j][0] = W[i][j][1];
                H[i][j][0] = H[i][j][1];
                U[i][j][kmax-1] = U[i][j][kmax-2];
                W[i][j][kmax-1] = W[i][j][kmax-2];
                H[i][j][kmax-1] = H[i][j][kmax-2];
            }
            P[i][0] = P[i][1];
            P[i][kmax-1] = P[i][kmax-2];
        }
        for (k = 1; k < kmax-1; k++) {
            A = (1.0-Z[k]*Z[k])/P[i][k];
            dPdx = (P[i+1][k]-P[i-1][k])/dx*0.5;
            dPdz = (P[i][k+1]-P[i][k-1])/dz*0.5;
            if (i < imax/L) {
                AU[0] = 0; BU[0] = 0;
                BW[0] = 0;
                AH[0] = 0; BH[0] = Hw;}
            else {
                AU[0] = 1; BU[0] = 0;
                BW[0] = 0;
                AH[0] = 1; BH[0] = 0;
            }
            for (j = 1; j < jmax-1; j++) {
                E = (Gamma-1)/Gamma*0.5*T[i][j][k];
                if (i != imax-1) dUdx = (U[i+1][j][k]-U[i-1][j][k]+U[i+1][j-1][k]-U[i-1][j-1][k])/dx*0.25;
                else dUdx = (U[i][j][k]-U[i-1][j][k]+U[i][j-1][k]-U[i-1][j-1][k])/dx*0.5;
                if (is_cartesian) {
                    Au = A*z0*X[i]*U[i][j][k]*dy*dy/dx;
                    if(i < imax/L) {
                        Aw = A*(W[i][j][k]-z0*U[i][j][k]*Z[k])*dy*dy/dz;
                        V[i][j][k] = V[i][j-1][k]+dy*0.25*(Z[k]/P[i][k]*(W[i][j][k]+W[i][j-1][k]-z0*Z[k]*(U[i][j][k]+U[i][j-1][k]))
                            -A*((W[i][j][k+1]-W[i][j][k-1]+W[i][j-1][k+1]-W[i][j-1][k-1])/dz+0.5*z0*(U[i][j][k]+U[i][j-1][k])
                            -z0*Z[k]*(U[i][j][k+1]-U[i][j][k-1]+U[i][j-1][k+1]-U[i][j-1][k-1])/dz+z0*X[i]*dUdx*4));
                        Eu = z0*dy*dy*E/P[i][k]*(A*(P[i][k]*0.5-X[i]*dPdx+Z[k]*dPdz)+Z[k]*Z[k]);
                    }
                    else {
                        Aw = A*(W[i][j][k]*X[i]-z0*U[i][j][k]*Z[k]*0)*dy*dy/dz;
                        V[i][j][k] = V[i][j-1][k]+dy*0.25*
                            (Z[k]/P[i][k]*((W[i][j][k]+W[i][j-1][k])*X[i]-z0*Z[k]*0*(U[i][j][k]+U[i][j-1][k]))
                            -A*((W[i][j][k+1]-W[i][j][k-1]+W[i][j-1][k+1]-W[i][j-1][k-1])/dz+0.5*z0*(U[i][j][k]+U[i][j-1][k])
                            -0*z0*Z[k]*(U[i][j][k+1]-U[i][j][k-1]+U[i][j-1][k+1]-U[i][j-1][k-1])/dz+z0*X[i]*dUdx*4));
                        Eu = z0*dy*dy*E/P[i][k]*(A*(P[i][k]*0.5-X[i]*dPdx+0*Z[k]*dPdz)+0*Z[k]*Z[k]);
                    }
                    Ew = -dy*dy*E/P[i][k]*(Z[k]+A*dPdz);
                }
                else {
                    Au = A*Theta0*X[i]*U[i][j][k]*dy*dy/dx;
                    Aw = A*W[i][j][k]*dy*dy/dz;
                    V[i][j][k] = V[i][j-1][k]+dy*0.25*(Z[k]/P[i][k]*(W[i][j][k]+W[i][j-1][k])-A*(Theta0*X[i]*dUdx*4
                        +(W[i][j][k+1]-W[i][j][k-1]+W[i][j-1][k+1]-W[i][j-1][k-1])/dz+2.5*Theta0*(U[i][j][k]+U[i][j-1][k])));
                    Eu = Theta0*dy*dy*A*(W[i][j][k]*W[i][j][k]+E*(0.5-X[i]/P[i][k]*dPdx));
                    Ew = -dy*dy*(A*Theta0*U[i][j][k]*W[i][j][k]+E*(Z[k]/P[i][k]+A/P[i][k]*dPdz));
                }
                Eh = (Pr-1)/Pr*(U[i][j+1][k]*U[i][j+1][k]-2*U[i][j][k]*U[i][j][k]+U[i][j-1][k]*U[i][j-1][k]
                              + W[i][j+1][k]*W[i][j+1][k]-2*W[i][j][k]*W[i][j][k]+W[i][j-1][k]*W[i][j-1][k]);
                vdy = V[i][j][k]*dy;
                if ( k == 1 || k == kmax-2 || is_first_order_theta ) {
                    AwdUz = ((Aw+fabs(Aw))*U[i][j][k-1]-(Aw-fabs(Aw))*U[i][j][k+1]);
                    AwdWz = ((Aw+fabs(Aw))*W[i][j][k-1]-(Aw-fabs(Aw))*W[i][j][k+1]);
                    AwdHz = ((Aw+fabs(Aw))*H[i][j][k-1]-(Aw-fabs(Aw))*H[i][j][k+1]);
                    fAw = 2*fabs(Aw);
                }
                else {
                    AwdUz = ((Aw+fabs(Aw))*(4*U[i][j][k-1]-U[i][j][k-2])-(Aw-fabs(Aw))*(4*U[i][j][k+1]-U[i][j][k+2]))*0.5;
                    AwdWz = ((Aw+fabs(Aw))*(4*W[i][j][k-1]-W[i][j][k-2])-(Aw-fabs(Aw))*(4*W[i][j][k+1]-W[i][j][k+2]))*0.5;
                    AwdHz = ((Aw+fabs(Aw))*(4*H[i][j][k-1]-H[i][j][k-2])-(Aw-fabs(Aw))*(4*H[i][j][k+1]-H[i][j][k+2]))*0.5;
                    fAw = 3*fabs(Aw);
                }
                if ( i == 1 || i == imax-2 ) {
                    AudUx = ((Au+fabs(Au))*U[i-1][j][k]-(Au-fabs(Au))*U[i+1][j][k]);
                    AudWx = ((Au+fabs(Au))*W[i-1][j][k]-(Au-fabs(Au))*W[i+1][j][k]);
                    AudHx = ((Au+fabs(Au))*H[i-1][j][k]-(Au-fabs(Au))*H[i+1][j][k]);
                    fAu = 2*fabs(Au);
                }
                else if (i == imax-1) {
                    AudUx = ((Au+fabs(Au))*(4*U[i-1][j][k]-U[i-2][j][k]))*0.5;
                    AudWx = ((Au+fabs(Au))*(4*W[i-1][j][k]-W[i-2][j][k]))*0.5;
                    AudHx = ((Au+fabs(Au))*(4*H[i-1][j][k]-H[i-2][j][k]))*0.5;
                    fAu = 3*fabs(Au);
                }
                else {
                    AudUx = ((Au+fabs(Au))*(4*U[i-1][j][k]-U[i-2][j][k])-(Au-fabs(Au))*(4*U[i+1][j][k]-U[i+2][j][k]))*0.5;
                    AudWx = ((Au+fabs(Au))*(4*W[i-1][j][k]-W[i-2][j][k])-(Au-fabs(Au))*(4*W[i+1][j][k]-W[i+2][j][k]))*0.5;
                    AudHx = ((Au+fabs(Au))*(4*H[i-1][j][k]-H[i-2][j][k])-(Au-fabs(Au))*(4*H[i+1][j][k]-H[i+2][j][k]))*0.5;
                    fAu = 3*fabs(Au);
                }
                AU[j] = (2-vdy)/(fAu+fAw+4-AU[j-1]*(2+vdy));
                AH[j] = (2/Pr-vdy)/(fAu+fAw+4/Pr-AH[j-1]*(2/Pr+vdy));
                BU[j] = (AudUx+AwdUz+BU[j-1]*(2+vdy)+2*Eu)/(fAu+fAw+4-AU[j-1]*(2+vdy));
                BW[j] = (AudWx+AwdWz+BW[j-1]*(2+vdy)+2*Ew)/(fAu+fAw+4-AU[j-1]*(2+vdy));
                BH[j] = (AudHx+AwdHz+BH[j-1]*(2/Pr+vdy)+2*Eh)/(fAu+fAw+4/Pr-AH[j-1]*(2/Pr+vdy));
            }
            for (j = jmax-2; j > 0; j--) {
                Us[i][j][k] = Us[i][j+1][k]*AU[j] + BU[j];
                Ws[i][j][k] = Ws[i][j+1][k]*AU[j] + BW[j];
                Hs[i][j][k] = Hs[i][j+1][k]*AH[j] + BH[j];
            }
        }
    }
    //relaxation
    for (i = 1; i < imax; i++)
        for (k = 1; k < kmax-1; k++)
            for (j = 1; j < jmax-1; j++) {
               U[i][j][k] = relU*Us[i][j][k] + (1-relU)*U[i][j][k];
               W[i][j][k] = relU*Ws[i][j][k] + (1-relU)*W[i][j][k];
               H[i][j][k] = relU*Hs[i][j][k] + (1-relU)*H[i][j][k];
            }
    }
}
double SelfNumbers(double P2, double z0, double Hw) {
    QString method = "Ust"; //Progonka, RK, Ust
    double alpha = 1, alpha2 = 2;
    double Integral;
    P[0][0] = 1;
    Delta[0][0] = Delta[1][0] = Delta[1][1] = Delta[0][1] = 0;
    U[0][0][0] = 0; W[0][0][0] = 0; H[0][0][0] = Hw; V[0][0][0] = 0;
    U[0][jmax-1][0] = 1; W[0][jmax-1][0] = 0; H[0][jmax-1][0] = 1;
    U[1][0][0] = 0; W[1][0][0] = 0; H[1][0][0] = 0; V[1][0][0] = 0;
    U[1][jmax-1][0] = 0; W[1][jmax-1][0] = 0; H[1][jmax-1][0] = 0;
    for (j = 1; j < jmax-1; j++) {
        //double func = pow((double)j/(jmax-1), 0.15);
        double func = tanh((double)j*dy*2);
        U[0][j][0] = U[0][jmax-1][0]*func;
        W[0][j][0] = W[0][jmax-1][0]*func;
        H[0][j][0] = H[0][0][0]*(1-func)+H[0][jmax-1][0]*func;
        U[1][j][0] = U[1][jmax-1][0]*func;
        W[1][j][0] = W[1][jmax-1][0]*func;
        H[1][j][0] = H[1][0][0]*(1-func)+H[1][jmax-1][0]*func;
    }
    if (method == "Ust") {
        dt = dy*dy/4.0001;
        //System C0 with W1
        Delta[0][1] = 2;
        for (NIterations = 0; (fabs(Delta[0][0]-Delta[0][1]) > 0.0000001) && NIterations < 1000000 && (converged == false) ; NIterations++) {
            Delta[0][1] = Delta[0][0];
            QCoreApplication::processEvents(); //to prevent freezing of gui
            for (j = 0; j < jmax; j++) {
                Us[0][j][0] = U[0][j][0];
                Hs[0][j][0] = H[0][j][0];
                Ws[0][j][0] = W[0][j][0];
            }
            Integral = H[0][0][0]-U[0][0][0]*U[0][0][0];
            for (j = 1; j < jmax-1; j+=2) Integral +=  4*(H[0][j][0]-U[0][j][0]*U[0][j][0]) + 2*(H[0][j+1][0]-U[0][j+1][0]*U[0][j+1][0]);
            Delta[0][0] = sqrt((Gamma-1)/Gamma*0.5)/P[0][0]*Integral*dy/3.0;
            P[0][0] = (1-relP)*P[0][0] + relP*(Gamma+1)*9/32.0*Delta[0][0]*Delta[0][0];
            for (j = 1; j < jmax-1; j++) {
                V[0][j][0] = V[0][j-1][0] - dy*(W[0][j][0]+0.25*z0*Us[0][j][0])/P[0][0];
                U[0][j][0] = Us[0][j][0]  - dt*(V[0][j][0]*(Us[0][j+1][0]-Us[0][j-1][0])/2/dy
                                          - z0*(Gamma-1)/4.0/Gamma/P[0][0]*(Hs[0][j][0]-Us[0][j][0]*Us[0][j][0])
                                          - (Us[0][j+1][0]-2*Us[0][j][0]+Us[0][j-1][0])/dy/dy);
                H[0][j][0] = Hs[0][j][0]  - dt*(V[0][j][0]*(Hs[0][j+1][0]-Hs[0][j-1][0])/2/dy
                                          + (1-Pr)/Pr*(Us[0][j+1][0]*Us[0][j+1][0]-2*Us[0][j][0]*Us[0][j][0]+Us[0][j-1][0]*Us[0][j-1][0])/dy/dy
                                          - 1/Pr*(Hs[0][j+1][0]-2*Hs[0][j][0]+Hs[0][j-1][0])/dy/dy);
                W[0][j][0] = Ws[0][j][0]  - dt*(V[0][j][0]*(Ws[0][j+1][0]-Ws[0][j-1][0])/2/dy
                                          + Ws[0][j][0]/P[0][0]*(Ws[0][j][0]-z0*Us[0][j][0])
                                          + (Gamma-1)/2.0/Gamma/P[0][0]*(Hs[0][j][0]-Us[0][j][0]*Us[0][j][0])*(1+2*P2/P[0][0])
                                          - (Ws[0][j+1][0]-2*Ws[0][j][0]+Ws[0][j-1][0])/dy/dy);
            }
        }
        qDebug() << "Delta[0][0]"<< Delta[0][0];
        qDebug() << "P[0][0]"<< P[0][0];
       }
//    qDebug() << SystemC1(P2, 0.75) << Delta[0][0];
    double Da = fabs(SystemC1(P2, (alpha+alpha2)*0.5));
    bool check = false;
    while (!check) {
        if (Da < Delta[0][0]) alpha = (alpha+alpha2)*0.5;
        else alpha2 = (alpha+alpha2)*0.5;
        Da = SystemC1(P2, (alpha+alpha2)*0.5);
        if (fabs(Da-Delta[0][0]) < 0.0001) check = true;
        if (fabs(alpha-alpha2) < 0.000001) {check = true; qDebug() << "Error";}
        qDebug() << alpha << alpha2 << SystemC1(P2, alpha) << SystemC1(P2, (alpha+alpha2)*0.5) << SystemC1(P2,alpha2);
    }
//    qDebug() << "Da =" << Delta[1][0] << "  D0 =" << Delta[0][0];
//    qDebug() << "alpha =" << alpha;
    return alpha;
}
double SystemC1(double P2, double alpha) {
    double Integral;
    P[1][0] = 8/3.0*(0.75-alpha)*P[0][0];
    Delta[1][1] = 1;
    for (NIterations = 0; (fabs(Delta[1][0]-Delta[1][1]) > 0.0000001) && NIterations < 100000 && (converged == false) ; NIterations++) {
        Delta[1][1] = Delta[1][0];
        QCoreApplication::processEvents(); //to prevent freezing of gui
        for (j = 0; j < jmax; j++) {
            Us[1][j][0] = U[1][j][0];
            Hs[1][j][0] = H[1][j][0];
        }
        Integral = (H[1][0][0]-2*U[0][0][0]*U[1][0][0])-P[1][0]/P[0][0]*(H[0][0][0]-U[0][0][0]*U[0][0][0]);
        for (j = 1; j < jmax-1; j+=2)
            Integral += 4*((H[1][j][0]-2*U[0][j][0]*U[1][j][0])-P[1][0]/P[0][0]*(H[0][j][0]-U[0][j][0]*U[0][j][0]))
                      + 2*((H[1][j+1][0]-2*U[0][j+1][0]*U[1][j+1][0])-P[1][0]/P[0][0]*(H[0][j+1][0]-U[0][j+1][0]*U[0][j+1][0]));
        Delta[1][0] = sqrt((Gamma-1)/Gamma*0.5)/P[0][0]*Integral*dy/3.0;
        //P[1][0] = (1-relP)*P[1][0] + relP *(Gamma+1)*(9/16.0*Delta[0][0]*Delta[0][0]-15/16.0*Delta[0][0]*Delta[1][0]);
        for (j = 1; j < jmax-1; j++) {
            V[1][j][0] = V[1][j-1][0]- dy/P[0][0]*((alpha+1)*Ws[1][j][0]+(0.25-alpha)*Us[1][j][0]
                                     - P[1][0]/P[0][0]*(Ws[1][j][0]+0.25*z0*Us[0][j][0]));
            U[1][j][0] = Us[1][j][0] - dt*(Us[1][j][0]*(Ws[0][j][0]-z0*Us[0][j][0])*alpha/P[0][0]
                                     + V[0][j][0]*(Us[1][j+1][0]-Us[1][j-1][0])/2/dy
                                     + V[1][j][0]*(Us[0][j+1][0]-Us[0][j-1][0])/2/dy
                                     - z0*(Gamma-1)/2.0/Gamma/P[0][0]*
                                     (0.5*(Hs[1][j][0]-2.0*Us[0][j][0]*Us[1][j][0])
                                     + (Hs[0][j][0]-Us[0][j][0]*Us[0][j][0])*P[1][0]/P[0][0]*(alpha-0.5))
                                     - (Us[1][j+1][0]-2*Us[1][j][0]+Us[1][j-1][0])/dy/dy);
            H[1][j][0] = Hs[1][j][0] - dt*(Hs[1][j][0]*(Ws[0][j][0]-z0*Us[0][j][0])*alpha/P[0][0]
                                     + V[0][j][0]*(Hs[1][j+1][0]-Hs[1][j-1][0])/2/dy
                                     + V[1][j][0]*(Hs[0][j+1][0]-Hs[0][j-1][0])/2/dy
                                     + (1-Pr)/Pr*(4*(Hs[0][j+1][0]-Hs[0][j-1][0])/dy*(Us[1][j+1][0]-Us[1][j-1][0])/dy/4.0
                                                + 2*Us[0][j][0]*(Us[1][j+1][0]-2*Us[1][j][0]+Us[1][j-1][0])/dy/dy
                                                + 2*Us[1][j][0]*(Us[0][j+1][0]-2*Us[0][j][0]+Us[0][j-1][0])/dy/dy)
                                     - 1/Pr*(Hs[1][j+1][0]-2*Hs[1][j][0]+Hs[1][j-1][0])/dy/dy);
            W[1][j][0] = Ws[1][j][0] - dt*(Ws[1][j][0]*(Ws[1][j][0]-z0*Us[1][j][0])/P[0][0]
                                        + (Ws[0][j][0]-z0*Us[0][j][0])*((alpha+1)*Ws[1][j][0]-Ws[0][j][0]*P[1][0]/P[0][0])/P[0][0]
                                     + V[0][j][0]*(Ws[1][j+1][0]-Ws[1][j-1][0])/2/dy
                                     + V[1][j][0]*(Ws[0][j+1][0]-Ws[0][j-1][0])/2/dy
                                     - (Gamma-1)/2.0/Gamma/P[0][0]*(
                                        alpha*P[1][0]/P[0][0]*Ws[1][j][0]*Ws[1][j][0]-(1+2*P2/P[0][0])*(Hs[1][j][0]-2*U[0][j][0]*U[1][j][0])
                                        +(Hs[0][j][0]-Us[0][j][0]*Us[0][j][0])*(alpha*P2*P[1][0]/P[0][0]+P[1][0]*(1+2*P2/P[0][0])
                                                                               +P[1][0]*P2*(alpha+2)/P[0][0]+alpha*P[1][0])/P[0][0])
                                     - (Ws[1][j+1][0]-2*Ws[1][j][0]+Ws[1][j-1][0])/dy/dy);
        }
    }
    return Delta[1][0];
}
void Translate(int N) { //traslation to physical variables
    for (i = 0; i < imax; i++) {
        for (k = 0; k < kmax; k++) DeltaE[i][k] = Delta[i][k]*pow((1-Z[k]*Z[k]), 0.75)*pow(X[i]+0.01, 0.75);
        DeltaE[i][0] = DeltaE[i][kmax-1] = 0.00001;
        for (k = 0; k < kmax; k++) {
            TauU[i][k] = (U[i][1][k]*cos(Theta0*Z[k]-Beta)+W[i][1][k]*sin(fabs(Theta0*Z[k]-Beta))
                         -U[i][0][0]*cos(Theta0*Z[k]-Beta)-W[i][0][k]*sin(fabs(Theta0*Z[k]-Beta)))/dy;
            Tau[i][k] = sqrt(U[i][1][k]*U[i][1][k] + W[i][1][k]*W[i][1][k])/dy;
            TauW[i][k] = fabs(U[i][1][k]*sin(Theta0*Z[k]-Beta)+W[i][1][k]*cos(Theta0*Z[k]-Beta))/dy;
            TauH[i][k] = (H[i][1][k]-H[i][0][k])/dy;
        }
        for (k = 1; k < kmax-1; k++) {
            P0[i][k] = P[i][k]*pow(1-Z[k]*Z[k],-0.20)*pow(X[i]+0.1,-0.5);
            Tau[i][k]  =  Tau[i][k]*P[i][k]*pow(1-Z[k]*Z[k],-0.20)*pow(X[i]+0.1,-0.5);
            TauU[i][k] = TauU[i][k]*P[i][k]*pow(1-Z[k]*Z[k],-0.75)*pow(X[i]+0.1,-0.5);
            TauW[i][k] = TauW[i][k]*P[i][k]*pow(1-Z[k]*Z[k],-0.75)*pow(X[i]+0.1,-0.5);
            TauH[i][k] = TauH[i][k]*P[i][k]*pow(1-Z[k]*Z[k],-0.25)*pow(X[i]+0.1,-0.5);
        }
        //extrapolation
        P0[i][0] = 3*P0[i][1]-3*P0[i][2]+P0[i][3];
        P0[i][kmax-1] = 3*P0[i][kmax-2]-3*P0[i][kmax-3]+P0[i][kmax-4];
        Tau[i][0] = 3*Tau[i][1]-3*Tau[i][2]+Tau[i][3];
        Tau[i][kmax-1] = 3*Tau[i][kmax-2]-3*Tau[i][kmax-3]+Tau[i][kmax-4];
        TauU[i][0] = 3*TauU[i][1]-3*TauU[i][2]+TauU[i][3];
        TauU[i][kmax-1] = 3*TauU[i][kmax-2]-3*TauU[i][kmax-3]+TauU[i][kmax-4];
        TauH[i][0] = 3*TauH[i][1]-3*TauH[i][2]+TauH[i][3];
        TauH[i][kmax-1] = 3*TauH[i][kmax-2]-3*TauH[i][kmax-3]+TauH[i][kmax-4];
    }
    for (i = 0; i < imax; i++)
        for (j = 0; j < jmax; j++)
            for (k = 0; k < kmax; k++) {
                Ro[i][j][k] = 2*Gamma/(Gamma-1)*P[i][k]/(sqrt(X[i]*(1-Z[k]*Z[k]))+0.0001)/T[i][j][k]/100;
            }
    double x, y, h = 0.001, phi, r0, fract, fract_x, fract_y, z_e;
    if (cartesian) {
        for (int iy = 0; iy < N; iy++)
            for (int ix = -N*0.5; ix < N*0.5; ix++) {
                i = int((imax-1)*2.0*z0*iy/N);
                fract_y = (imax-1)*2.0*z0*iy/N - i;
                z_e = (Ze[i]*(1-fract_y)+Ze[i+1]*fract_y)/z0*N*0.5;
                k = int((ix/z_e+1)*(kmax-1)*0.5);
                fract_x = fabs(ix/z_e*(kmax-1)*0.5+(kmax-1)*0.5) - k;
                if(fabs(ix) < z_e && i < imax-1) {
                    PressureDistribution[ix+int(N*0.5)][iy] =
                            (P0[i][k]*(1-fract_y)+P0[i+1][k]*fract_y)*(1-fract_x)+
                            (P0[i][k+1]*(1-fract_y)+P0[i+1][k+1]*fract_y)*fract_x;
                    DeltaDistribution[ix+int(N*0.5)][iy] =
                            (DeltaE[i][k]*(1-fract_y)+DeltaE[i+1][k]*fract_y)*(1-fract_x)+
                            (DeltaE[i][k+1]*(1-fract_y)+DeltaE[i+1][k+1]*fract_y)*fract_x;
                    TauUDistribution[ix+int(N*0.5)][iy] =
                            (Tau[i][k]*(1-fract_y)+Tau[i+1][k]*fract_y)*(1-fract_x)+
                            (Tau[i][k+1]*(1-fract_y)+Tau[i+1][k+1]*fract_y)*fract_x;
                    TauWDistribution[ix+int(N*0.5)][iy] =
                            (TauW[i][k]*(1-fract_y)+TauW[i+1][k]*fract_y)*(1-fract_x)+
                            (TauW[i][k+1]*(1-fract_y)+TauW[i+1][k+1]*fract_y)*fract_x;
                    TauHDistribution[ix+int(N*0.5)][iy] =
                            (TauH[i][k]*(1-fract_y)+TauH[i+1][k]*fract_y)*(1-fract_x)+
                            (TauH[i][k+1]*(1-fract_y)+TauH[i+1][k+1]*fract_y)*fract_x;
                }
            }
    }
    else {
        for (i = 0; i < N; i++)
            for (k = 0; k < N; k++) {
                x = h*(i-N/2);
                y = h*(k-2*N/5);
                phi = atan2(x,y)/Theta0 + Beta/Theta0;
                r0 = sqrt(x*x+y*y);
                if (fabs(phi) <= 1) {
                    n = (int)floor(phi/dz) + kmax/2;
                    fract = phi/dz-floor(phi/dz);
                    DeltaDistribution[i][k] = (DeltaE[0][n]*(1-fract)+DeltaE[0][n+1]*fract)*pow(r0,0.75);
                    PressureDistribution[i][k] = (P0[0][n]*(1-fract)+P0[0][n+1]*fract)*pow(r0+0.07,-0.5);
                    TauUDistribution[i][k] = (Tau[0][n]*(1-fract)+Tau[0][n+1]*fract)*pow(r0+0.07,-0.5);
                    TauWDistribution[i][k] = TauW[0][n];
                    TauHDistribution[i][k] = (TauH[0][n]*(1-fract)+TauH[0][n+1]*fract)*pow(r0+0.07,-0.5);
                }
            }
    }
}
void Analysis() { //computing of characteristics
    //MACH NUMBERS
    for (i = 0; i < imax; i++)
        for (k = 0; k < kmax; k++)
            for (j = 0; j < jmax; j++) {
                T[i][j][k] = H[i][j][k]-U[i][j][k]*U[i][j][k]-W[i][j][k]*W[i][j][k];
                Mach[i][j][k] = sqrt((U[i][j][k]*U[i][j][k]+W[i][j][k]*W[i][j][k])/(Gamma-1)*2/T[i][j][k]);
            }
    //AERODYNAMICS COEFFICIENTS
//	double Cp;
//	Cp = 0;
//	for (i = 1; i < imax; i++)
//		for (k = 1; k < kmax-1; k++)
//			Cp += P0[i][k]*dx*dz;
//	Cp *= 2/sqrt(z0);
//	qDebug() << "Cp = " << Cp;
    //for (k = 0; k < kmax; k++) qDebug() << Z[k];
//	QFile file(QString("results/Mach=1_%1-%2,Hw=%3,Pr=%4,G=%5,%6x%7x%8.dat")
//			   .arg(Theta0d).arg(Betad).arg(Hw).arg(Pr).arg(Gamma).arg(imax).arg(jmax).arg(kmax));
//	file.open(QFile::WriteOnly | QFile::Text);
//	QTextStream out(&file);
//	out << "#Theta Mach=1\n";
    double M = 0;
    //for (float v = 0.5; v <= 10; v += 0.5) {
        //qDebug() << "\nv = " << v;
        for (k = 0; k < kmax; k++) {
            for (j = 1; j < jmax; j++) {
                M = Mach[0][j][k]-1;
                if(M > 0) {J[k] = j; /*(j-1+fabs(Ms)/(fabs(Ms)+fabs(M)))*dy;*/ j = jmax;}
            }
//			qDebug() << Z[k] << J[k] << J[k]*pow(1-Z[k]*Z[k], 0.25) << DisturbancesDistribution[180][k];
        }
    double Integral;
    for (k = 0; k < kmax; k++) {
        Integral = T[0][0][k];
        for (j = 1; j <= J[k]; j++) {
            if (j%2 == 1) Integral +=  4*T[0][j][k];
            else Integral += 2*T[0][j][k];
        }
        Delta99[0][k] = sqrt((Gamma-1)/Gamma*0.5)/P[0][k]*Integral*dy/3.0;
    }
//	qDebug() << "Hw=" << Hw;
    /*for (k = 0; k < kmax; k++) {
        for (j = 0; j < jmax; j++) {
//			Delta99[0][k] += (1-sqrt(U[0][j][k]*U[0][j][k]+W[0][j][k]*W[0][j][k]))*dy;
            if (sqrt(U[0][j][k]*U[0][j][k]+W[0][j][k]*W[0][j][k]) > 0.99) {
                Delta99[0][k] = (j-(sqrt(U[0][j][k]*U[0][j][k]+W[0][j][k]*W[0][j][k])-0.99)
                        /(sqrt(U[0][j][k]*U[0][j][k]+W[0][j][k]*W[0][j][k])
                         -sqrt(U[0][j-1][k]*U[0][j-1][k]+W[0][j-1][k]*W[0][j-1][k])))*dy;
                //qDebug() << sqrt(U[0][j][k]*U[0][j][k]+W[0][j][k]*W[0][j][k])-0.99 << sqrt(U[0][j][k]*U[0][j][k]+W[0][j][k]*W[0][j][k])-sqrt(U[0][j-1][k]*U[0][j-1][k]+W[0][j-1][k]*W[0][j-1][k]) << Delta99[0][k] << (j-1)*dy;
                j = jmax;
            }
        }
//        if (k == kmax/2) qDebug()  << "M=1 - " << J[k] << "\nDelta99 = " << Delta99[0][k];
    }*/
//	file.close();

    //CRITICAL FLOWS
//	double dk, omega1 = 0, omega2 = 0;
//	for (k = 0; k < kmax; k++) {
//		J[k] = T[0][0][k];
//		for (j = 1; j < jmax-1; j+=2)
//			J[k] +=   4*((Gamma-1)/2.0*pow(T[0][j][k]/(W[0][j][k]),2) - T[0][j][k])
//				 + 2*((Gamma-1)/2.0*pow(T[0][j+1][k]/(W[0][j+1][k]),2) - T[0][j+1][k]);
//		J[k] *= dy/3.0;
//	}
//	for (k = 0; k < kmax-1; k++) {
//		if (J[k]*J[k+1] < 0) {
//			dk = fabs(J[k])/(fabs(J[k])+fabs(J[k+1]));
//			if(J[k] < 0) omega1 = Theta0d + Betad - fabs(k+dk)*dz*Theta0*180/3.14159;
//			else omega2 = fabs(k+dk)*dz*Theta0*180/3.14159 - (Theta0d + Betad);
//			//qDebug() << "dk =" << dk << "k =" << k << "  J[k] =" << J[k];
//		}
//	}
//	}
    //}
//	for(j = 0; j < jmax; j++) qDebug() << j*dy;
//	for(j = 0; j < jmax; j++)
//		qDebug() << H[0][j][25];//sqrt(U[0][j][25]*U[0][j][25]+W[0][j][25]*W[0][j][25]);
    //qDebug() << omega2;
    //qDebug() << "k =" << k << ";   omega =" << (1-fabs(k-kmax/2)*dz)*Theta0*180/3.14159 << ";   Pi/2-Theta0 = " << 90-Theta0d;
//    QFile file("results/for_graphs.dat");
//    file.open(QFile::WriteOnly | QFile::Text);
//    QTextStream out(&file);
//    out << "V";
//    for(j = 0; j < jmax; j++)
//        out << /*QString("\n%1").arg(j*dy) <<*/ QString("\n%1").arg(sqrt(U[0][j][54]*U[0][j][54]+W[0][j][54]*W[0][j][54]));
//    out << "\n\n\nH";
//    for(j = 0; j < jmax; j++)
//        out << QString("\n%1").arg(H[0][j][54]);
//    file.close();
//	QFile file("results/pressure_contours.dat");
//	file.open(QFile::WriteOnly | QFile::Text);
//	QTextStream out(&file);
//	for(k = 0; k < N; k++)
//		for(j = 0; j < N; j++)
//			out << QString("{%1,").arg(k) << QString(" %1,").arg(j) << QString("%1},\n").arg(PressureDistribution[j][k]);
//	file.close();
//	QFile file("results/test.dat");
//	file.open(QFile::WriteOnly | QFile::Text);
//	QTextStream out(&file);
//	out << "#Pressure\n";
//	for(k = 0; k < N; k++) {
//		for(j = 0; j < N; j++)
//			out << QString("%1 ").arg(PressureDistribution[j][k]);
//		out << "\n";
//	}
//	file.close();
}
void Disturbances(double angle_min, double angle_max) {
    double af, af2, a1;
    for (k = 0; k < kmax; k++) {
        QCoreApplication::processEvents(); // to prevent gui freezing
        for (int angle = angle_min; angle <= angle_max; angle++) {
            af = 0.01; af2 = 2.0;
            a1 = Simps(af,angle,k);
            while (fabs(a1) > 0.00001 && af > 0.000001) {
                if (a1*Simps(af2,angle,k) < 0) {
                    if (a1*Simps((af+af2)/2.0,angle,k) < 0) af2 = (af+af2)/2.0;
                    else af = (af+af2)/2.0;
                  }
                else { af = af*0.5; af2 = af2 + 0.1;}
                a1 = Simps(af,angle,k);
              }
            DisturbancesDistribution[angle][k] = af;
            //angle-Theta0d*(kmax-1-k)*dz+Theta0d-180-2*Betad << "	" << af << "\n";
          }
      }
    //qDebug() << Vw << DisturbancesDistribution[180][kmax/2+int(Betad/Theta0d*kmax/2)];
    //	for (k = 0; k < kmax; k++)
    //	    qDebug() /*<< Z[k] << J[k] << J[k]*pow(1-Z[k]*Z[k], 0.5)*/ << DisturbancesDistribution[180][k];
    //	for(k = 0; k < kmax; k++) {
    //	    for (int angle = 0; angle <= 360; angle+=0.5) {
    //		af = 0.02; af2 = 1.2;
    //		    while (fabs(Simps(af,angle,k)) > 0.0001 && af > 0.00001) {
    //			if (Simps(af,angle,k)*Simps(af2,angle,k) < 0) {
    //			    if (Simps(af,angle,k)*Simps((af+af2)/2.0,angle,k) < 0) af2 = (af+af2)/2.0;
    //				else af = (af+af2)/2.0;
    //			}
    //			else { af = af/2.0; af2 = af2*1.5;}
    //		}
    //		if (fabs(angle-Theta0d*(kmax-1-k)*dz+Theta0d-90-2*Betad) < 0.25) qDebug() << Z[k] << ";" << af;
    //	}
    //}
}
double Simps(double velocity, int angle, int k) {
    double I = T[0][0][k];// alpha = Theta0*Z[k]-Beta+Pi/2;
    if (!cartesian)
        for (j = 1; j < jmax-1; j+=2)
            I += 4*((Gamma-1)/2.0*pow(T[0][j][k]/(velocity
                    -U[0][j][k]*cos(angle/180.0*Pi-(Theta0*Z[k]-Beta))-W[0][j][k]*sin(angle/180.0*Pi-(Theta0*Z[k]-Beta)))
                    /*-(U[0][j][k]*sin(alpha)+W[0][j][k]*cos(alpha))*cos(angle/180.0*Pi)
                    -(U[0][j][k]*cos(alpha)+W[0][j][k]*sin(alpha))*sin(angle/180.0*Pi))*/,2)-T[0][j][k])
               + 2*((Gamma-1)/2.0*pow(T[0][j+1][k]/(velocity
                    -U[0][j][k]*cos(angle/180.0*Pi-(Theta0*Z[k]-Beta))-W[0][j][k]*sin(angle/180.0*Pi-(Theta0*Z[k]-Beta)))
                    /*-(U[0][j][k]*sin(alpha)+W[0][j][k]*cos(alpha))*cos(angle/180.0*Pi)
                    -(U[0][j][k]*cos(alpha)+W[0][j][k]*sin(alpha))*sin(angle/180.0*Pi))*/,2)-T[0][j+1][k]);
    else
        for (j = 1; j < jmax-1; j+=2)
            I += 4*((Gamma-1)/2.0*pow(T[0][j][k]/(velocity-U[0][j][k]*cos(angle/180.0*Pi)-W[0][j][k]*sin(angle/180.0*Pi)),2)-T[0][j][k])
               + 2*((Gamma-1)/2.0*pow(T[0][j+1][k]/(velocity-U[0][j][k]*cos(angle/180.0*Pi)-W[0][j][k]*sin(angle/180.0*Pi)),2)-T[0][j+1][k]);
    return I*dy/3.0;
}
void FreeArray3D(double ***array, int N1, int N2) {
    for (i = 0; i < N1; i++) {
        for (j = 0; j < N2; j++) {
            delete []array[i][j];
        }
        delete []array[i];
    }
    delete []array;
}
void FreeArray2D(double **array, int N1) {
    for (i = 0; i < N1; i++) {
        delete []array[i];
    }
    delete []array;
}
void DeleteArrays(bool is_interpolated) {
    if (U == NULL) return;
    if (!is_interpolated) {
        FreeArray3D(Us, imax, jmax);
        FreeArray3D(Vs, imax, jmax);
        FreeArray3D(Ws, imax, jmax);
        FreeArray3D(Hs, imax, jmax);

        FreeArray2D(Prel, imax);
        FreeArray2D(DeltaRel, imax);
    }
    FreeArray3D(U, imax, jmax);
    FreeArray3D(V, imax, jmax);
    FreeArray3D(W, imax, jmax);
    FreeArray3D(H, imax, jmax);
    FreeArray3D(T, imax, jmax);
    FreeArray3D(Mach, imax, jmax);
    FreeArray3D(Ro, imax, jmax);

    FreeArray2D(P,  imax+1);
    FreeArray2D(P0, imax);
    FreeArray2D(Pd, imax);
    FreeArray2D(Delta,  imax);
    FreeArray2D(DeltaE, imax);
    FreeArray2D(Delta99, imax);
    FreeArray2D(Body, imax+1);

    FreeArray2D(Tau,  jmax);
    FreeArray2D(TauU, jmax);
    FreeArray2D(TauW, jmax);
    FreeArray2D(TauH, jmax);

    FreeArray2D(TauUDistribution, N);
    FreeArray2D(TauWDistribution, N);
    FreeArray2D(TauHDistribution, N);
    FreeArray2D(DeltaDistribution, N);
    FreeArray2D(PressureDistribution, N);
    FreeArray2D(DisturbancesDistribution, 360);

    delete []X;
    delete []Y;
    delete []Z;
    delete []Ze;
    delete []AU;
    delete []BU;
    delete []BW;
    delete []AH;
    delete []BH;
    delete []AP;
    delete []BP;
    delete []Mu;
}
void InterpolateArrays(int imaxnew, int jmaxnew, int kmaxnew, int imaxold, int jmaxold, int kmaxold) {
    if (imaxnew > imaxold) {
        for (i = 0; i < imaxnew-1; i+=2)
            for (k = 0; k < kmaxnew; k++) {
                for (j = 0; j < jmaxnew; j++) {
                    U[i][j][k] = Us[i/2][j][k];
                    U[i+1][j][k] = (Us[i/2][j][k]+Us[i/2+1][j][k])*0.5;
                    W[i][j][k] = Ws[i/2][j][k];
                    W[i+1][j][k] = (Ws[i/2][j][k]+Ws[i/2+1][j][k])*0.5;
                    H[i][j][k] = Hs[i/2][j][k];
                    H[i+1][j][k] = (Hs[i/2][j][k]+Hs[i/2+1][j][k])*0.5;
                }
                P[i][k] = Prel[i/2][k];
                P[i+1][k] = (Prel[i/2][k]+Prel[i/2+1][k])*0.5;
            }
        for (k = 0; k < kmaxnew; k++) {
            for (j = 0; j < jmaxnew; j++) {
                U[imaxnew-1][j][k] = Us[imaxold-1][j][k];
                W[imaxnew-1][j][k] = Ws[imaxold-1][j][k];
                H[imaxnew-1][j][k] = Hs[imaxold-1][j][k];
            }
            P[imaxnew-1][k] = Prel[imaxold-1][k];
            Delta[0][k] = DeltaRel[0][k];
        }
    }
    //delete old arrays
    FreeArray3D(Us, imaxold, jmaxold);
    FreeArray3D(Ws, imaxold, jmaxold);
    FreeArray3D(Hs, imaxold, jmaxold);
    FreeArray3D(Vs, imaxold, jmaxold);
    FreeArray2D(Prel, imaxold);
    FreeArray2D(DeltaRel, imaxold);
    //create new arrays
    Us = NewArray3D(imaxnew, jmaxnew, kmaxnew);
    Ws = NewArray3D(imaxnew, jmaxnew, kmaxnew);
    Hs = NewArray3D(imaxnew, jmaxnew, kmaxnew);
    Vs = NewArray3D(imaxnew, jmaxnew, kmaxnew);
    Prel = NewArray2D(imaxnew, kmaxnew);
    DeltaRel = NewArray2D(imaxnew, kmaxnew);

    for (i = 0; i < imaxnew; i++)
        for (j = 0; j < jmaxnew; j++)
            for (k = 0; k < kmaxnew; k++) {
                Us[i][j][k] = 0;
                Ws[i][j][k] = 0;
                Hs[i][j][k] = 0;
            }
}
QString SelfTest() {
    imax=2;
    jmax=141;
    kmax=100;
    dx=1;
    dy=0.04;
    dz=0.04;
    N=500;
    Hw=0.5;
    Theta0d = 125; Theta0 = Theta0d*Pi/180;
    Betad = 35; Beta = Betad*Pi/180;
    CreateArrays(false, imax, jmax, kmax, dx, dz, N);
    SetBoundaryConditions(false, Hw);
    Edge(Theta0, 0, false, -1);
    Edge(Theta0, 0, false,  1);
    Nos(false);
    return "Finished";
}
void NaSt2D() {
    double laplU, laplV, du2dx, duvdy, dv2dy, duvdx, resid = 0, Re = 1000;
    double ***F = NewArray3D(imax, jmax, 1); // 1 is for 2D case
    double ***G = NewArray3D(imax, jmax, 1);
    double **rhsP = NewArray2D(imax, jmax);
    dt = dy*dy/4.1;
    LZ = 0;
    k  = 0;
    //Boundary conditions
    for (j = 0; j < jmax; j++) {
        //left
        U[0][j][k] = 0;
        V[0][j][k] = 0;
        //right
        U[imax-1][j][k] = 0;
        V[imax-1][j][k] = 0;
    }
    for (i = 0; i < imax; i++) {
        //top
        U[i][jmax-1][k] = 2;
        V[i][jmax-1][k] = 0;
        //bottom
        U[i][0][k] = 0;
        V[i][0][k] = 0;
    }
    for (NIterations = 0; /*fabs(U[imax-1][2][0]/Us[imax-1][2][0])>0.00001*/!converged && NIterations < 1000000; NIterations++) {
        QCoreApplication::processEvents(); //to prevent freezing of gui
        //outlet boundary conditions
        //for (j = 0; j < jmax; j++) {
        //   U[imax-1][j][0] = 2*U[imax-2][j][0]-U[imax-3][j][0];
        //   V[imax-1][j][0] = 2*V[imax-2][j][0]-V[imax-3][j][0];
        //}
        dt = 0.9*qMin(Re/2/(1/dx/dx+1/dy/dy),dx/U[i][jmax-1][k]/*,dy/U[i][jmax-1][k]*/);
        for (i = 1; i < imax-1; i++)
          for (j = 1; j < jmax-1; j++) {
              laplU = (U[i+1][j][0]-2*U[i][j][0]+U[i-1][j][0])/dx/dx
                     +(U[i][j+1][0]-2*U[i][j][0]+U[i][j-1][0])/dy/dy;
              du2dx = (U[i][j][0]+U[i+1][j][0])*(U[i][j][0]+U[i+1][j][0])/4/dx
                     -(U[i-1][j][0]+U[i][j][0])*(U[i-1][j][0]+U[i][j][0])/4/dx;
              duvdy = (V[i][j][0]+V[i+1][j][0])*(U[i][j][0]+U[i][j+1][0])
                     -(V[i][j-1][0]+V[i+1][j-1][0])*(U[i][j-1][0]+U[i][j][0])/4/dy;
              F[i][j][0] = U[i][j][0] + dt*(laplU/Re-du2dx-duvdy);
          }
        for (i = 1; i < imax-1; i++)
          for (j = 1; j < jmax-1; j++) {
              laplV = (V[i+1][j][0]-2*V[i][j][0]+V[i-1][j][0])/dx/dx
                     +(V[i][j+1][0]-2*V[i][j][0]+V[i][j-1][0])/dy/dy;
              duvdx = (U[i][j][0]+U[i][j+1][0])*(V[i][j][0]+V[i+1][j][0])/4/dx
                     -(U[i-1][j][0]+U[i-1][j+1][0])*(V[i-1][j][0]+V[i][j][0])/4/dx;
              dv2dy = (V[i][j][0]+V[i][j+1][0])*(V[i][j][0]+V[i][j+1][0])
                     -(V[i][j-1][0]+V[i][j][0])*(V[i][j-1][0]+V[i][j][0])/4/dy;
              G[i][j][0] = V[i][j][0] + dt*(laplV/Re-duvdx-dv2dy);
              rhsP[i][j] = ((F[i][j][0]-F[i-1][j][0])/dx+(G[i][j][0]-G[i][j-1][0])/dy)/dt;
          }
        //SOR iterations for pressure equation
        double omeg = 1.7;
        for (int iter=1; (iter < 1000) || (resid < 0.0001); iter++) {
            for (i = 1; i < imax-1; i++)
              for (j = 1; j < jmax-1; j++) {
                  P[i][j] = (1-omeg)*P[i][j] + omeg*(-rhsP[i][j]);
                }
        }
        for (i = 1; i < imax-1; i++)
          for (j = 1; j < jmax-1; j++) {
              U[i][j][0] = F[i][j][0] - dt/dx*(P[i+1][j]-P[i][j]);
              V[i][j][0] = G[i][j][0] - dt/dy*(P[i][j+1]-P[i][j]);
            }
        //for (i = 1; i < imax-1; i++)
        //  for (j = 1; j < jmax-1; j++) {
        //      U[i][j][0] = Us[i][j][0] + dt/Re*( (Us[i+1][j][0]-2*Us[i][j][0]+Us[i-1][j][0])/dx/dx
        //                 +(Us[i][j+1][0]-2*Us[i][j][0]+Us[i][j-1][0])/dy/dy)
        //                 -dt*Us[i][j][0]*(Us[i+1][j][0]-Us[i-1][j][0])/dx
        //                 -dt*Vs[i][j][0]*(Us[i][j+1][0]-Us[i][j-1][0])/dy;
        //      V[i][j][0] = Vs[i][j][0] + dt/Re*( (Vs[i+1][j][0]-2*Vs[i][j][0]+Vs[i-1][j][0])/dx/dx
        //                 +(Vs[i][j+1][0]-2*Vs[i][j][0]+Vs[i][j-1][0])/dy/dy)
        //                 -dt*Us[i][j][0]*(Vs[i+1][j][0]-Vs[i-1][j][0])/dx
        //                 -dt*Vs[i][j][0]*(Vs[i][j+1][0]-Vs[i][j-1][0])/dy;
        //  }
    } //End of global iterations
}
double Laplacian(double ***Func, double dx, double dy) {
    return (Func[i+1][j][0]-2*Func[i][j][0]+Func[i-1][j][0])/dx/dx
          +(Func[i][j+1][0]-2*Func[i][j][0]+Func[i][j-1][0])/dy/dy;
}
double f(int j, double Fj, double dFj) {
    j++; Fj++;
    return dFj;
}
double g(int j, double Fj, double dFj, double ***F, int i, int k) {
    Fj++;
    //U00
    if ((F == U) && (i == 0) && (k == 0)) return dFj*V[0][j][0]-(Gamma-1)/4.0/Gamma/P[0][0]*(H[0][j][0]-U[0][j][0]*U[0][j][0]);
    //H00
    if ((F == H) && (i == 0) && (k == 0)) {
        if (j == 0)
            return dFj*Pr*V[0][j][0]+(1-Pr)*(U[0][j+2][0]*U[0][j+2][0]-2*U[0][j+1][0]*U[0][j+1][0]+U[0][j][0]*U[0][j][0])/dy/dy;
        else
            return dFj*Pr*V[0][j][0]+(1-Pr)*(U[0][j+1][0]*U[0][j+1][0]-2*U[0][j][0]*U[0][j][0]+U[0][j-1][0]*U[0][j-1][0])/dy/dy;
    }
    //U10
    if((F == U) && (i == 1) && (k == 0)) return dFj*V[0][j][0];
    else return 0;
}

//Solving system of ODE's:
//F'[i][j][k] = dF[j]
//dF' = C0[i][j][k] * F'[i][j][k] + C1[i][j][k]*F[i][j][k] + C2[i][j][k]
double RungeKutta2 (double ***F, int i, int k, double dy, double dF0) {
    double A1, A2, B1, B2, *dF;
    dF = new double [jmax];
    dF[0] = dF0;//(F[i][1][k]-F[i][0][k])/dy;
    for (int j = 0; j < jmax-1; j++) {
        A1 = dy*f(j, F[i][j][k], dF[j]); //f=dF[j];
        B1 = dy*g(j, F[i][j][k], dF[j], F, i, k);
        A2 = dy*f(j+1, F[i][j][k] + A1, dF[j] + B1); //f=dF[j]+B1
        B2 = dy*g(j+1, F[i][j][k] + A1, dF[j] + B1, F, i, k);
        dF[j+1] = dF[j] + 0.5*(B1+B2);
        F[i][j+1][k] = F[i][j][k] + 0.5*(A1+A2);
        //if(fabs(F[i][j+1][k]) > 3) j = jmax;
    }
    return F[i][jmax-1][k];
}
double RungeKutta4 (double ***F, int i, int k, double dy, double dF0) {
    double A1, A2, A3, A4, B1, B2, B3, B4, *dF;
    dF = new double [jmax];
    dF[0] = dF0;
    for (int j = 0; j < jmax-1; j++) {
        A1 = f(j, F[i][j][k], dF[j]); //f=dF[j];
        B1 = g(j, F[i][j][k], dF[j], F, i, k);
        A2 = f(j+0.5, F[i][j][k] + 0.5*A1, dF[j] + 0.5*B1);
        B2 = g(j+0.5, F[i][j][k] + 0.5*A1, dF[j] + 0.5*B1, F, i, k);
        A3 = f(j+0.5, F[i][j][k] + 0.5*A2, dF[j] + 0.5*B2);
        B3 = g(j+0.5, F[i][j][k] + 0.5*A2, dF[j] + 0.5*B2, F, i, k);
        A4 = f(j+1, F[i][j][k] + A3, dF[j] + B3);
        B4 = g(j+1, F[i][j][k] + A3, dF[j] + B3, F, i, k);
        dF[j+1] = dF[j] + dy*(B1+B2+B3+B4)/6.0;
        F[i][j+1][k] = F[i][j][k] + dy*(A1+A2+A3+A4)/6.0;
        if (fabs(F[i][j+1][k]) > 3) j = jmax;
    }
    return F[i][jmax-1][k];
}

/*
void Progonka (int N, double *Array) {
    double *A, *B; //coefficients
    A = new double [N];
    B = new double [N];
}
*/
