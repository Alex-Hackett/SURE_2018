void D2polpscat(double f[1000][6], double Ninph[1000], double Noutph[1000], int m, double del[1000],
		int zzz, int Geo, double kpp[1000], double delk[1000], int y, double index, double kcz, int disp)
{
  double GQ1in, GQ1out;
  double Fq0, Fq00, GQ1in0, GQ1in00, GQ1out0, GQ1out00, Ninph0, Noutph0, r, r1;
  double f0;
  if (y == 1 || disp == 3)
    M2(W, kpp, defpote, defpoth, index, kcz, W0, Wden0, Wden00); /* initialize the scattering rate array */
  /* Dp = 4 * pi / (hb * hb * hb * v * v * v); */
  Ninph[p+3] = 0;
  Noutph[p+3] = 0;
  /* LOOP E */
  for(i = 0; i <= p; i++)
    {
      r = 0;
      k = kpp[i];
      E = f[i][0];
      Ninph[i] = 0;
      Noutph[i] = 0;
      Ninph0 = 0;
      Noutph0 = 0;
      f0 = f[i][1];
      /* dEdk0 = getdEdk(E, k); */ /* use if integrating over energy */
      /* LOOP E1 */
      for(ii = 0; ii <= p; ii++)
	{
	  r1 = 0; /* made random again 02/12/08 */
	  k1 = kpp[ii];
	  E1 = f[ii][0];
	  dEdk1 = dEdki[ii]; /* getdEdk(E1, k1, disp, index, 0, y); */
	  /* if(del[i] > (0.01 * kb * T)) 9/19/07 change */
	  /*f0 = f[i][1]; */ /* findf(i, f, delk, r);*/ /*f[i][1] + ((f[i+1][1] - f[i][1]) * r);*/
	  /*Determinef4(E, k, f, m, i, zzz, Geo, kpp, delk, del);*/
	  /* else */
	  /*f0 = f[i][1] + ((f[i+1][1] - f[i][1]) * r); */
	  /* if(del[ii] > (0.01 * kb * T)) 9/19/07 change */
	  f1 = f[ii][1]; /*findf(ii, f, delk, r1);*/ /* Determinef4(E1, k1, f, m, ii, zzz, Geo, kpp, delk, del); */
	  /*else */
	  /*f1 = f[ii][1] + ((f[ii+1][1] - f[ii][1]) * r1); */ /* 9/19/07 change */
	  /* if(f0 < 0) f0 = (f[i-1][1] + f[i][1] + f[i+1][1] + f[i+2][1]) / 4;
	     if(f1 < 0) f1 = (f[ii-1][1] + f[ii][1] + f[ii+1][1] + f[ii+2][1]) / 4; 9/19/07 change */
	  if (E != E1)
  {
    Fq = 1 / (exp(fabs(E - E1) / (kb * T)) - 1);
    if (E > E1)
      {
	GQ1in = (1 + f0) * f1 * Fq; /* replace this 6/3/07 */
	GQ1out = f0 * (1 + f1) * (1 + Fq); /* and this */
      }
    else
      {
	GQ1in = (1 + f0) * f1 * (1 + Fq); /* and this */
	GQ1out = f0 * (1 + f1) * Fq; /* and this 6/3/07 */
      }
    dfac = (Nink[i][ii] + ((Nink[i+1][ii+1] - Nink[i][ii]) * sqrt(r * r + r1 * r1) / sqrt(2))) * kpp[ii] * delk[ii] * 2;
    Ninph[i] += dfac * GQ1in;
    Noutph[i] += dfac * GQ1out;
    /* printf("dfac = %e, GQ1in = %e GQ1out = %e\n", dfac, GQ1in, GQ1out); */
  }
	  /* scattering rate in vs |q|*/ /* uncomment fclose and exit below */
	  /* if(i == 0 && ii == 0)
{
sprintf(outfile4, "Rvsq.d");
fu = fopen(outfile4, "w");
fprintf(fu, "x ");
for(j = 0 ;j <= p; j++)
fprintf(fu, "%e ", kpp[j]);
fprintf(fu, "\n");
}
if(ii == 0)
fprintf(fu, "%e ", kpp[i]);
if(E == E1 || dfac != dfac || GQ1out != GQ1out)
fprintf(fu,"0 ");
else
fprintf(fu, "%e ", dfac * GQ1out * f[i][2]);
if(ii == p)
fprintf(fu, "\n");*/
	  if(i == 0 && ii != 0)
	    Fq0 = 1 / (exp((E1 - f[0][0]) / (kb * T)) - 1);
	  else
	    Fq0 = 2 * kb * T / del[0];
	  if(ii == 0 && i != 0)
	    Fq00 = 1 / (exp((E - f[0][0]) / (kb * T)) - 1);
	  else
	    Fq00 = 2 * kb * T / del[0];
	  if(i == 0)
	    {
	      GQ1out0 = f[p+3][1] * (1 + f1) * Fq0; /* and this */
	      GQ1in0 = (1 + f[p+3][1]) * f1 * (1 + Fq0); /* and this 6/3/07 */
	      Ninph[p+3] += W0[ii] * del[ii] * GQ1in0 / (dEdk1 * Wden0[ii]);
	      Noutph[p+3] += W0[ii] * del[ii] * GQ1out0 / (dEdk1 * Wden0[ii]);
	    }
	  if(ii == 0)
	    {
	      GQ1in00 = f[p+3][1] * (1 + f0) * Fq00; /* and this 6/3/07 */
	      GQ1out00 = (1 + f[p+3][1]) * f0 * (1 + Fq00); /* and this */
	      Ninph0 = W0[i] * GQ1in00 / Wden00[i];
	      Noutph0 = W0[i] * GQ1out00 / Wden00[i];
	    }
	  if (Ninph[i] > 1E100 || Noutph[i] > 1E100 || Ninph[i] != Ninph[i])
	    {
	      printf("error in 2Dpolpscat.h\n");
	      printf("GQ1in = %e, GQ1out = %e, f0 = %e, f1 = %e\n", GQ1in, GQ1out, f0, f1);
	      printf("Fq = %e, E - E1 = %e\n", Fq, E - E1);
	      exit(1);
	    }
	  /* first E1 if statement */
	  /* if((i == 0 || i == 1) && ii == 0)
	     printf("Ninph %e Noutph %e Ninph0 %e Noutph0 %e\n", Ninph[i], Noutph[i], Ninph0, Noutph0);*/
	  if (Noutph[i] < 0)
	    {
	      printf("GQ1 out %e k %e dEdk1 %e dfac %e\n", GQ1out, k, dEdk1, dfac);
	      printf("i %d ii %d\n", i, ii);
	      printf("Ninks: %e %e %e \n", Nink[i][ii], Nink[i+1][ii+1], sqrt(r * r + r1 * r1));
	      printf("f0 = %e f1 = %e\n", f0, f1);
	      printf("f(i) = %e f(i+1) = %e f(i+2) = %e\n", f[i][1], f[i+1][1], f[i+2][1]);
	      printf("GQ1in = %e f0 = %e f1 = %e Fq = %e\n", GQ1in, f0, f1, Fq);
  printf("f(i) = %e f(i+1) = %e f(i+2) = %e r = %e\n", f[30][1], f[31][1], f[32][1], r);
 printf("E = %e f(E)= %e\n", E, 1 /(exp((E - mu - f[0][0])/(kb * T)) - 1));
 printf("E1 = %e f(E1) = %e\n", E1, 1/(exp((E1 - mu - f[0][0])/(kb * T)) - 1));
 printf("E(i) = %e E(i+1) = %e\n", f[30][0], f[31][0]);
 exit(1);
	    }
	} /* LOOP E1 */
      /* if(i == 1) exit(1);*/
      /*
printf("E = %f Eq1 = %f Eq2 = %f\n", E/kb/T, Eq1/kb/T, Eq2/kb/T);
printf("N41 = %1.12f N42 = %1.12f N5 = %1.12f N61 = %1.12f N62 = %1.12f\n", N41, N42, N5, N61, N62);
printf("qGQ1 = %f\n", qGQ1);
printf("qmax1 = %f qmin1 = %f\n", qmax1, qmin1);
printf("GQ1in = %f\n", GQ1in);
printf("qdif1 = %f qdif2 = %f phcoef = %f\n", qdif1, qdif2, phcoef);
printf("Gqtotalin1 = %f GQtotalin2 = %f\n", GQtotalin1, GQtotalin2);
printf("i = %d in = %f out = %f\n\n", i, Ninph[i], Noutph[i]);
      */
      /*printf("Nin0 = %e Nout0 = %e Nin = %e Nout = %e\n", Ninph0, Noutph0, Ninph[i], Noutph[i]);*/
      Ninph[i] += Ninph0;
      Noutph[i] += Noutph0;
      Ninph[i] *= f[i][2];
      Noutph[i] *= f[i][2];
      /*
if(initial == 2 && ( (fabs(Ninph[i] - Noutph[i])/ Ninph[i]) > 4e-15) )
{
printf("%d LA phonon difference = %e\n", i, fabs(Ninph[i] - Noutph[i])/Ninph[i]);
exit(1);
}
      */
      /*printf("%d %e %e\n", i, Ninph[i], Noutph[i]);*/
    } /* LOOP E */
  /*exit(1);*/
  Ninph[p+3] *= DOS[p+3];
  Noutph[p+3] *= DOS[p+3];
  /*fclose(fu);
    exit(1);*/
}
