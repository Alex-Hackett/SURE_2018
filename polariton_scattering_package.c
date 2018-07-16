#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define hbar 6.582119514e-16
#define kb 8.6173303e-5
#define m0 9.10938356e-31
#define E_X_0 1.515
#define me 0.067 * m0
#define mh 0.45 * m0
#define a0 10e-9
#define n 3.43
#define hbar_omega 5e-3
#define De -8.6
#define Dh 5.7
#define u 4810
#define rho 5.3e3
#define lz 1.192970175673805e-07
#define S 1e-10
#define V 1
#define c 299792458
#define pi 3.14159

int main(int argc,char *argv[] );
double dis(double k);
double Hop_X(double k);
double F(double q);
double B(double q);
double D(double q);
double d_kk(double k, double k_prime);
double phonon_proj(double k, double k_prime);
double phonon_number(double k, double k_prime, double temp);

int main(int argc,char *argv[] )
{
  double k = atof(argv[1]);
  double k_prime = atof(argv[2]);
  double temp = atof(argv[3]);
  double u_k = Hop_X(k);
  double u_k_prime = Hop_X(k_prime);
  double delta_k_k_prime = d_kk(k, k_prime);
  double qz = phonon_proj(k, k_prime);
  double first_term = (lz * pow((u_k*u_k_prime*delta_k_k_prime),2))/(hbar * rho * V * qz * u*u);
  double B_term = B(qz);
  double B_sq = B_term*B_term;
  double D_term = D(abs(k - k_prime));
  double D_sq = D_term * D_term;
  double N_ph = phonon_number(k,k_prime,temp);
  double last_term = (delta_k_k_prime - abs(k - k_prime));
  double W = first_term * B_sq * D_sq * N_ph * last_term;
  printf("%f",W);
  return (0);
}

double dis(double k)
{
  double E_C = pow((((hbar * c)/(n)) * (pow(((pi)/(lz)),2)) + (pow(k,2))),0.5);
  double E_X = E_X_0;

  double E_L = 0.5 * (E_C + E_X - sqrt((pow(E_C-E_X,2))+pow(hbar_omega,2)));
  return (E_L);
}

double Hop_X(double k)
{
  double E_X = E_X_0;
  double E_L = dis(k);
  double u_k = 1 / (sqrt(1 + pow(((hbar_omega/2)/(E_L - E_X)),2) ));
  return (u_k); 
}

double F(double q)
{
  return pow((1 + pow(((q*a0)/(2)),2)),(-3.0/2.0));
}

double B(double q)
{
  return ((8*pi*pi)/(lz*q*((4*pi*pi)-(lz*lz * q*q)))) * (sin((lz*q)/(2)));
}

double D(double q)
{
  double Fh = F((q*mh)/(me+mh));
  double Fe = F((q*me)/(me+mh));
  return (De * Fh) - (Dh * Fe);
}

double d_kk(double k, double k_prime)
{
  double E_k_prime = dis(k_prime);
  double E_k = dis(k);
  double delta_kk = (abs(E_k_prime - E_k))/(hbar*u);
  return delta_kk;
}

double phonon_proj(double k, double k_prime)
{
  double delta_kk = d_kk(k, k_prime);
  double qz = pow(((delta_kk*delta_kk) - pow((abs(k - k_prime)),2)),0.5);
  return qz;
}

double phonon_number(double k, double k_prime, double temp)
{
  double E_k_prime = dis(k_prime);
  double E_k = dis(k);
  if(E_k_prime - E_k > 0){
    return 1/(exp((E_k_prime - E_k)/(kb*temp)) - 1);
  }
  if(E_k_prime - E_k < 0){
    return 1 + (1/(exp((E_k_prime - E_k)/(kb*temp))-1));
  }
  if(E_k_prime - E_k == 0){
    return (0);
  }
}







