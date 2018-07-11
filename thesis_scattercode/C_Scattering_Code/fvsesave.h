int fvsEsave(double del[1000], int xx, double N[1000], double f[1000][6], double Nin[1000],
	     double Nout[1000], double tau, double taut, char dataname[25], double Etotal,
	     int Norf, const char fname[12], int y, int o, double Ntotal, int Geo, double kpp[1000], double Ninpol[1000],
	     double Ninpolel[1000], double Ninph[1000], double Avetau[1])
/*, double data[200], double delk[1000]) */
{
  Tmid = Ntotal * pi * hb * hb * eC / (dg * (Me + Mh) * em * kb);
  if( y == 1 || y == 2 || y == o )
    {
      sprintf(outfile1, "fvsEdata.%d.d", xx);
      fq = fopen(outfile1, "w");
      if (scattype == 13 || scattype == 15)
	fprintf(fq, "Energy DOS Number Occ# Nin-Nout Nin Nout Ninpolel Ninph\n");
      else if (scattype == 14)
	fprintf(fq, "Energy DOS Number Occ# Nin-Nout Nin Nout Ninpol Ninph Ninpolel\n");
      else
	fprintf(fq, "Energy DOS Number Occ# Nin-Nout Nin Nout Ninpol Ninph\n");
      if (scattype == 13 || scattype == 15)
	fprintf(fq, "%e %e %e %e %e %e %e %e %e\n",
		f[0][0]/(kb*T), DOS[p+3], N[p+3], f[p+3][1], Nin[p+3] - Nout[p+3], Nin[p+3], Nout[p+3], Ninpolel[p+3], Ninph[p+3]);
      else if (scattype == 14)
	fprintf(fq, "%e %e %e %e %e %e %e %e %e %e\n",
		f[0][0]/(kb*T), DOS[p+3], N[p+3], f[p+3][1], Nin[p+3] - Nout[p+3], Nin[p+3], Nout[p+3], Ninpol[p+3], Ninph[p+3], Ninpolel[p+3]);
      else
	fprintf(fq, "%e %e %e %e %e %e %e %e %e\n",
		f[0][0]/(kb * T), DOS[p+3], N[p+3], f[p+3][1], Nin[p+3] - Nout[p+3], Nin[p+3], Nout[p+3], Ninpol[p+3], Ninph[p+3]);
      for (mmm = 0; mmm <= p; mmm++)
	{
	  if (scattype == 13 || scattype == 15)
	    fprintf(fq, "%e %e %e %e %e %e %e %e %e\n",
		    f[mmm][0] / (kb * T), DOS[mmm], N[mmm]/del[mmm], f[mmm][1], Nin[mmm] - Nout[mmm], Nin[mmm], Nout[mmm], Ninpolel[mmm], Ninph[mmm]);
	  else if (scattype == 14)
	    fprintf(fq, "%e %e %e %e %e %e %e %e %e %e\n",
		    f[mmm][0] / (kb * T), DOS[mmm], N[mmm]/del[mmm], f[mmm][1], Nin[mmm] - Nout[mmm],
		    Nin[mmm], Nout[mmm], Ninpol[mmm], Ninph[mmm], Ninpolel[mmm]);
	  else
	    fprintf(fq, "%e %e %e %e %e %e %e %e %e\n",
		    f[mmm][0] / (kb * T), DOS[mmm], N[mmm]/del[mmm], f[mmm][1], Nin[mmm] - Nout[mmm],
		    Nin[mmm], Nout[mmm], Ninpol[mmm], Ninph[mmm]);
	}
      fclose(fq);
    }
  if (xx == 0)
    {
      sprintf(outfile2, fname);
      fr = fopen(outfile2, "w");
      fprintf(fr, fname);
      fprintf(fr, "\nNtotal = %e\n", Ntotal);
      fprintf(fr, "\nave E = %f", Etotal / kb / T / Ntotal);
      fprintf(fr, "\np = %d\n", p);
      fprintf(fr, "uprate = %f\n", uprate[0]);
      fprintf(fr, "GQp = %d\n", GQp);
      fprintf(fr, "Geo = %d\n", Geo);
      fprintf(fr, "initial = %d\n", initial);
      fprintf(fr, "qo = %e Ne = %e\n", qo[0], Ne[0]);
      fprintf(fr, "count = %d\n", count);
      fprintf(fr, "delc = %d\n", delc);

  fprintf(fr, "disp = %d\n", disp);
 fprintf(fr, "iterations = %d\n", o);
 fprintf(fr, "Norf = %d\n", Norf);
 fprintf(fr, "GQp = %d\n", GQp);
 fprintf(fr, "Geo = %d\n", Geo);
 fprintf(fr, "scattype = %d\n", scattype);
 fprintf(fr, "statype = %d\n", statype);
 fprintf(fr, "T = %f, TTT = %f\n", T, TTT[0]);
 fprintf(fr, "%d ", 0);
 for ( mmm = 0; mmm <= p; mmm++)
   fprintf(fr, "%e ", del[mmm]);
 fprintf(fr, "\n");
 fprintf(fr, "0 "); /* energy points on mesh */
 for ( mmm = 0; mmm <= p; mmm++)
   fprintf(fr, "%e ", (f[mmm][0] - f[0][0])*1000); /* displays Energy in meV */
 fprintf(fr, "tau ");
 fprintf(fr, "Etotal/kb/T ");
 fprintf(fr, "Stime[0] ");
 fprintf(fr, "Stime[(p-2)/2] ");
 fprintf(fr, "Stime[p-2] ");
 fprintf(fr, "Stime[p+3] ");
 fprintf(fr, "taut ");
 fprintf(fr, "Tmid/log(1+f[0][1]) "); /* temperature for an equilibrium distribution, numerator related to n and nQ, 0.00045401*/
 fprintf(fr, "Nscat[1] Nscat[p/2] Nout[1] Nout[p/2] ");
 fprintf(fr, "Nin[0]/N[0] ");
 fprintf(fr, "Ntau ");
 fprintf(fr, "Ntotal ");
 fprintf(fr, "Uprate ");
 fprintf(fr, "Avetau ");
 if (disp == 3)
   fprintf(fr, "Every_other_line_is_Renormalized_energies");
 fprintf(fr, "\n");
 if ( Geo == 5)
   fprintf(fr,"%e ", kpp[0]);
 else
   fprintf(fr,"%e ", (f[0][0]/(kb * T)));
 for ( mmm = 0; mmm <= p; mmm++)
   {
     if ( Geo == 5)
       fprintf(fr,"%e ", kpp[mmm]);
     else
       fprintf(fr,"%e ", (f[mmm][0]/(kb * T)));
   }
 fprintf(fr,"\n");
    }
  if (Norf == 1)
    {
      fprintf(fr, "%e ", N[p+3]);
      for ( mmm = 0; mmm <= p+3; mmm++)
	/* fprintf(fr, "%e ", f[mmm][3]); */
	fprintf(fr, "%e ", N[mmm]/del[mmm]); /* N/E */
    }
  else
    {
      fprintf(fr, "%e ", (f[p+3][1]));
      for ( mmm = 0; mmm <= p; mmm++)
	fprintf(fr, "%e ", (f[mmm][1]));
    }
  fprintf(fr, "%1.6f ", tau);
  fprintf(fr, "%1.6f ", Etotal / kb / T);
  fprintf(fr, "%e ", Stime[0]);
  fprintf(fr, "%e ", Stime[(p-2)/2]);
  fprintf(fr, "%e ", Stime[p-2]);
  fprintf(fr, "%e ", Stime[p+3]);
  fprintf(fr, "%e ", taut);
  fprintf(fr, "%e ", Tmid / log (1 + f[0][1]) ); /* temperature for an equilibrium distribution, numerator related to n and nQ, 0.00045401*/
  fprintf(fr, "%e %e %e %e ", Nscat[1], Nscat[p/2], Nout[1], Nout[p/2]);
  fprintf(fr, "%e ", Nin[0] / N[0]);
  fprintf(fr, "%e ", Ntau);
  fprintf(fr, "%e ", Ntotal);
  fprintf(fr, "%e ", uprate[0]);
  fprintf(fr, "%e ", Avetau[0]);
  if (disp == 3)
    {
      fprintf(fr, "\n");
      fprintf(fr, "0 ");
      for (mmm = 0; mmm <= p; mmm++)
	fprintf(fr, "%e ", (f[mmm][0] - f[0][0]) * 1000); /* displays energy in meV */
  }
  fprintf(fr, "\n");
  /*
if ( xx == o)
{
fprintf(fr, "program time (in hours) = %f\n", (float)(clock())/3600000000); */
  if ( y == o)
    {
      fprintf(fr, "%e ", (double)(f[p+3][1]));
      for (mmm = 0; mmm <= p; mmm++)
	fprintf(fr, "%e ", (double)(f[mmm][1]));
      fprintf(fr, "final uprate = %e ", uprate[0]);
      fclose(fr);
    }
  return(1);
}
