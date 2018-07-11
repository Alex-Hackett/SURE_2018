#include<stdio.h>
#include<string.h>
#include<math.h>
#include "constants.h" /* file where constants are defined */
#include "variables.h" /* file where variable are defined */
#include "parameters.h" /* file where important parameters for a specific calculation are defined */
#include<stdlib.h>
#include<string.h>
/*#include "Errors.h" */
#include "polpara.h"
#include "Ecset.h"
#include "getk3.h"
#include "dEdk.h"
#include "DOS.h" /* function for calculating the density of states */
#include "rand.h"
#include "del.h" /* function for creating the mesh for the calculation */
#include "getf4.h"
#include "getb4.h"
#include "polcheck.h"
#include "polEnergy.h"
#include "polfraction.h" /* function for determining the excitonic fraction of a polariton */
#include "ktilii.h"
#include "gaussian.h" /* gaussian f vs E */
#include "boltzmann.h" /* boltzmann distribution */
/*#include "boltzmann2.h" */
#include "pulseadd.h" /* pulse input, thermalized */
#include "pulsegauss.h" /* pulse input, gaussian */
#include "pulseflat.h" /* pulse input, all states equally populated */
/*#include "pulseflat2.h" */
#include "neareq.h" /* input an instantenous nearly thermalized distribution */
#include "phononfvsE.h"
#include "updatef.h" /* function for updating occupation numbers for scattering, decay, and pumping*/
#include "polpMatx.h" /* determine the scattering matrix elements for polariton-longitudinal acoustic phonon scattering */
#include "polelMatx.h"
#include "polpTAMatx.h" /* determines the scattering matrix elements for polariton-transverse acoustic phonon scattering */
#include "polpTA_PiezoMatx.h"
#include "polpFmatx.h" /* determine the scattering matrix elements for polariton-optical phonon scattering */
#include "bospMatx.h" /* determine the scattering matric elements for exciton-exciton scattering */
#include "findf.h"
#include "3Dbosescat.h" /* 3D boson-boson scattering integral */
#include "3Dfermiscat.h" /* 3D fermi-fermi scattering integral */
#include "2Dbosescat.h" /* 2D boson-boson scattering integral */
#include "2Dpolscat.h" /* 2D polariton-polariton scattering integral */
#include "2Dpolelscat.h"
#include "2DpolFscat.h" /* 2D polariton-optical phonon scattering integral */
#include "2DxpEexchange.h" /* 2D exciton-acoutic phonon scattering integral */
#include "2DpolpEexchange.h"
#include "3DxpEexchange.h"
#include "2DpolpLAscat.h" /* 2D polariton-3D longitudinal acoustical phonon scattering integral */
#include "2DpolpTAscat.h" /* 2D polariton-3D transverse acoustical phonon scattering integral */
#include "2Dbospscat.h" /* 2D exciton- 3D acoustical phonon scattering integral */
#include "Vvsq.h"
#include "uniform.h" /* uniform f vs E */
/*#include "datafit.h" */
#include "fvsEsave.h" /* function for saving calculated values */
#include "GQtest.h"
#include "gauleg.h" /* function for generating the Gaussian Quadrature points and weights */
#include "fromfile.h" /* function for reading in values from a text file */
#include "initiate.h" /* function for directing the initialization of the calculation */
#include "scat.h" /* function for directing the type of scattering to occur */
/*#include "MinEn.h" */
#include "polReNorm.h" /* function for renormalizing the dispersion curve of polaritons */
#include "zeta.h"
int main()
{
  inj1 = inj;
  Ntau = 0;
  indexr = sqrt(einf);
  kcz = setkcz(detune, Exo, kxy); /* in polpara.h */
  Ec = setEc(kcz, kxy);
  gauleg(-1,1,yy,ww,GQp);
  mu *= kb * TT;
  /*if (initial == 1)
    {*/
  Tm = zeta(2, exp(mu / (kb * T))) / (-log(1 - exp(mu / (kb * T)))); /* zeta(dimensionality, ) */
  /*
g = log(1 - exp(mu / (kb * T))) * log(1 - exp(mu / (kb * T))) / (2 * zeta(3, exp(mu / (kb * T))));
}*/
  hw *= kb * T;
  srand( (unsigned)time( NULL ) ); /* seed the random numbers from clock*/
  m = -1; /* m = # of steps to max energy for uniformfvsE */
  /* determined in uniform.h */
  initiate(initial, del, p, DOS, m, f, Nin, Nout, N, mu, hw, expp, Tm, statype, Ec, indexr); /* initial particle distribution */
  if(initial == 0) initial = 9; /* after reading in the data a flat pulse is added in */
  Ntotal = 0;
  Etotal = 0;
  /* determine total number of particles and total energy */
  for(i = 0; i <= p; i++)
    {
      Ntotal += N[i];
      Etotal += N[i] * (f[i][0] + f[i+1][0]) / 2;
    }
  Ntotal += N[p+3];
  Nbegin = Ntotal;
  fvsEsave(del, 0, N, f, Nin, Nout, tau, taut, dataname, Etotal, Norf, fname, y, o, Ntotal, Geo, kpp, Ninpol, Ninpolel, Ninph, Avetau);
do
  {
    /* do the scattering integrals */
    scatter(scattype, f, Nin, Nout, m, zzz, del, Geo, kpp, delk, y, indexr, kcz, disp);
    /* update the occupation numbers for scattering, pumping, and recombination */
    x = updatef(Ntotal, uprate, x, dN, N, f, Nin, Nout, Geo, kpp, Stime, taut, Avetau);
    /* renormalize polariton distribution if set for in parameters.h */
    if(disp == 3)
      ReNorm(indexr);
    /* determine total number of particles and total energy */
    Etotal = 0;
    Ntotal = 0;
    dNtotal = 0;
    for(l = 0; l <= p; l++)
      {
	Ntotal += N[l];
  dNtotal += Nin[l] * del[l];
 Etotal += N[l] * (f[l][0] + f[l+1][0]) / 2;
      }
    Ntotal += N[p+3];
    /* determine time step */
    change = x / Ntotal;
    tau += change;
    taut += Stime[p+3];
    Ntau += Stime[p+3] * N[p+3];
    if(Stime[p+3] < 1e-16) uprate[0] *= 1.1;
    if(uprate[0] > 1) uprate[0] = 1;
    /* record the calculation at the interval set in parameters.h */
    if(y > (n1 * count) - 1)
      {
	fvsEsave(del, n1, N, f, Nin, Nout, tau, taut, dataname, Etotal, Norf, fname, y, o, Ntotal, Geo, kpp, Ninpol, Ninpolel, Ninph, Avetau);
	n1++;
	x = 0;
      }
    y++;
    sprintf(outfile3, ynum);
    ft = fopen(outfile3, "w");
    fprintf(ft, "y = %d\n", y);
    fclose(ft);
  }
 while((n1*count <= o) && (taut < maxtime));
 }
