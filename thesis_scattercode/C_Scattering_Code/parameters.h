int count = 1000; /* number of iterations per data writing */
double defpote = -7; /* hydrostatic deformation potential in eV for electron, (Pikus-Bir, "a" for conduction band) */
double defpoth = 2.7; /* hydrostatic deformation potential in eV for hole, (Pikus-Bir, "a" for valence band) */
double defpoteT = 0; /* shear deformation potential in eV for electron, (Pikus-Bir, "b" and "d" for conduction band */
double defpothT = 3.8; /* shear deformation potential in eV for hole, (Pikus-Bir, ((4/5)*(b^2 + d^2/2))^1/2, GaAs --> b = 1.8, d = 5.4 */
int delc = 5; /* key for del.h
0: not sure
1: e ^ i, e is in constants.h
2: i ^ 2
3: i ^ 3
4: uniform
5: polariton i
	      */
int disp = 2; /* dispersion relationship to be used,
1 = exciton
2 = polariton
3 = renormalized
	      */
double detune = 0.00; /* polaritons --> cavity mode detuning from resonance with excitons, eV */
int dopiezo = 2; /* use piezoelectric interaction with phonons
1 = no
2 = yes
		 */
double g = 0.1; /* initial occupation number for uniform dist., not used right now */
const int Geo = 5; /* free particle dimensionality (1,2,3), harmonic potential (4), free polariton (5) */
const int GQp = 25; /* Number of GQ points */
double hw = 0.001; /*0.001 */ /* harmonic potential ground state in kb T */
#define inj 5e18 /*2E18 405.4725 1000 2000*/ /*injected gaussian density */
int initial = 9; /* initial determines the initial configuration of
particles;
0: from file
1: uniform, fermi distribution
2: bose, equilibrium distribution
3: bose phonon distribution, gaussian free boson distribution
4: gaussian distribution, free particles
gaussian distribution, harmonic potential
5: 2D bose phonon distribution, gaussian free
boson distribution
6: from near equilibrium
7: pulsed thermalized input
8: pulsed gaussian input
9: pulse with flat occupation number distribution
		 */
const double kpb = 0.25; /* base wavevector for polaritons */
double kxy = 5330000; /* polaritons --> injected in plane wavevector, m^-1 */
const double maxkp = 1.54E8; /*1.54E8; */ /* maximum k-parallel used in polaritons, m^-1 */
const double maxkT = 10; /*10.5;*/ /* maximum energy used in calculation in kT */
double maxtime = 5e-1; /* maximum time calculation will run */
double mu = -4;
double mue[1] = {-3.88}; /* chemical potential of free electrons in units of kT (eV) */
const char fname[20] = {"pol05.25.08e.d"}; /* data file name */
const char ynum[10] = {"y1value"}; /* file to check calculation progress */
const int Norf = 2; /* flag for saving density(1) or occupancy(other) */
int o = 40000; /* number of iterations */
int p = 170; /* p + 2 is number of energy points */
double pulseT = 1e-8; /* pulse length for intial = 7, 8, or 9*/
double Sa = 1; /* confinement size of system for ground state, m^2 */
int scattype = 14; /* type of scattering
1: 3D flat boson-boson
2: 2D flat boson-boson
3: 3D boson-phonon
4: 2D boson-phonon
5: gaussian quadrature testing
6: 2D polariton-polariton
7: stirctly 2D, polariton-phonon(edit before use, 7/24/04)
8: 2D polariton - 3D acoustical phonon(longitudinal and transverse)
9: 2D polariton - 2D polariton/3D acoustic phonon(longitudinal and transverse)
10: 3D flat fermi-fermi
11: 2D polariton - 2D polariton/3D acoustic and optical phonon
12: 2D boson - 2D boson/3D acoustic phonon
13: 2D polariton - free electron / 3D acoustic phonon(longitudinal and transverse)
14: 2D polariton - 2D polariton / free electron / 3D acoustic phonon(longitudinal and transverse)
15: 2D polariton - free electron
*/
int ScreenType = 2; /* flag for the way the screening is handled,
1: epsilon -> epsilon * (1 + qo / (k - k’))
2: epsilon -> epsilon * (1 + qo * aB)
		    */
int statype = 1; /* statistics used
1: Bose-Einstein
2: Boltzmann
3: Fermi-Dirac
		 */
double uprate[1] = {0.1}; /* {1E16}{0.4};*/ /* fraction of total number of particles
					       actually scattered per iteration */
double TT = 4; /* Used for setting a different initial temperature than lattice */
double TTT[1] = {10}; /* Temperature of free electrons */
