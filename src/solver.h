#ifndef SOLVER_H
#define SOLVER_H

extern double Pi, h, Theta0d, Theta0, Theta1d, Theta1, Betad, Beta, Gamma, Pr, Hw, Vw, Omega, khi, z0;
extern double relU, relP, relP2, a0, AccuracyU, AccuracyPressure;
extern int imax, jmax, kmax, N, LX, LY, LZ, NIterations, MaxIterations;
extern double dx, dy, dz, dt, L;
extern double *X, *Y, *Z, *Ze, **Body,
***U, ***W, ***H, ***V, ***T, ***Mach, ***Ro, ***Us, ***Ws, ***Hs, **Tau, **TauU, **TauW, **TauH, **TauUDistribution,
**TauWDistribution, **TauHDistribution, **DeltaDistribution, **PressureDistribution, **DisturbancesDistribution,
**Delta, **DeltaE, **Delta99, **P, **P0, **Prel, *Pd,
*AU, *BU, *BW, *AH, *BH, *AP, *BP, *Mu, *J;
extern bool converged, cartesian, is_first_order_theta, thomas_method,	is_mod_rel, draw_grid;
extern QString status;
double NewArray3D(int N1, int N2, int N3);
double NewArray2D(int N1, int N2);
void CreateArrays(bool is_interpolated, int imax, int jmax, int kmax, double dx, double dz, int N);
void SetZe(int imax, double dx, double Theta0);
void SetBoundaryConditions(bool, double);
void Edge(double, int, bool, int);
void Nos(bool);
void Simm();
double SelfNumbers(double P2, double z0, double Hw);
double SystemC1(double P2, double alpha);
double RungeKutta2(double ***F, int i, int k, double dy, double dF0);
void FullWing(bool, bool is_interpolated);
void FullWingWithTrail(bool);
void FreeArray3D(double ***array, int N1, int N2);
void FreeArray2D(double **array, int N1);
void DeleteArrays(bool is_interpolated);
void Analysis();
void Translate(int N);
void Disturbances(double angle_min, double angle_max);
void InterpolateArrays(int imax, int jmax, int kmax, int NXold, int NYold, int NZold);
QString SelfTest();

void NaSt2D();
double Laplacian(double ***Func, double dx, double dy);
#endif
