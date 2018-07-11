double D2A = 0.0001; /* 2D sample area in m^2 */
double aB = 130E-10; /* exciton Bohr radius in m */
double beta = 2; /* polariton - exciton volumic oscillator strength */
const int dg = 4; /* degeneracy */
const double dn = 1.28E20; /* particle density in particles/m^3 */
const double e = 1.2; /*2.718281828;*/ /* in energy step size. see del.h */
double EB = 0.010; /* exciton binding energy in eV */
const double eC = 1.60219E-19; /* charge on an electron in Coulombs */
#define Ex0 1.61945 /* polariton --> exciton base energy in eV */
double einf = 10.9; /* dielectric constant at infinity, set for GaAs */
double estat = 12.5 ; /* static dielectric constant */
#define em 9.1095E-31 /* electron rest mass */
#define eo 8.85E-12 /* permitivity of free space, C^2 / N m^2 */
double Exo = 1.61945;
const double gs = 0.0149; /* polariton splitting at resonance in eV, 2 times the Rabi frequency */
#define hb 6.5822E-16 /* planck’s constant in eV s */
#define Itaup 0.2e12 /* 0.2e12 */ /* inverse of photon lifetime, s^-1 */
#define Itaux 1e6 /* 1e6 */ /* inverse of exciton lifetime, s^-1 */
#define kb 0.00008617 /* boltzmann constant in eV/K */
double Lata = 5.6533E-10; /* Lattice constant, set for GaAs, meters */
double Lz = 7E-9; /* quantum well thickness, meters */
double Me = 0.067; /* exciton electron mass, in em, electron rest mass Piermarocchhi PRB 53 15834 (1996)*/
double Mh = 0.18; /* For exciton total mass, exciton hole mass, in units of em, electron rest mass */
double Mh1 = 0.08; /* Used for calculating Beta in polariton-photon interaction */
const double Pd = 6; /* 3D crystal density */
double piezo14 = 0.16; /* piezoelectic constant, set for GaAs e_(14) */
double qo[1]; /* screening parameter coefficient, e^2 / (2 * e(inf)) */
#define pi 3.14159265358979
const double Sd = 5316; /* set for GaAs */ /* 2D density kg/m^3 */
double sig = 0.5; /*0.00001; 0.000196764; */ /* injected gaussian width */
const double T = 10; /* Lattice temp in K */
#define vc 2.998E8 /* speed of light, m/s */
#define v 5.117E3 /*rough average for GaAs */ /* longitudinal speed of sound in medium, m/s */
#define vTA 3.012E3 /* transverse speed of sound in medium, m/s */
#define wLO 1.070E13 /* longitudinal optical phonon frequency, GaAs */
/* see parameters.h for count, g, o, p, uprate, and others */
