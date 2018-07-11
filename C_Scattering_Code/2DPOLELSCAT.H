int D2polelscat(double f[1000][6], double Ninpolel[1000], double Noutpolel[1000], int m,
		double del[1000], int zzz, int Geo, double kpp[1000], double delk[1000], double index, int disp, int y)
{
  double Fin0, Fout0, Fin00, Fout00, dEdk3, den00, rado, dEdki2, dEdki3, kxk;
  nnn=0;
  /* include a loop for each integration variable in this list */
  /* E, E1, E2 */
  if (y == 1 || disp == 3) polelMatx(f, Ninpolel, Noutpolel, m, del, zzz, Geo, kpp, delk, index, disp, y);
  /* LOOP E */
  for (i = 0; i <= p; i++)
    {
      Ninpolel[i] = 0;
      Noutpolel[i] = 0;
      /* LOOP E1 */
      for (ii = 0; ii <= p; ii++)
	{
	  /* LOOP E2 */
	  for (iii = 0; iii <= p; iii++)
	    {
	      f5 = fe[I[i][ii][iii]] + ((fe[I[i][ii][iii]+1] - fe[I[i][ii][iii]]) * fnew[i][ii][iii]);
	      Fin = f[ii][1] * fe[iii] * (1 + f[i][1]) * (1 - f5); /* took out 1 - f5, 02/13/08*/
	      Fout = f[i][1] * f5 * (1 + f[ii][1]) * (1 - fe[iii]); /* took out 1 - fe[iii], 02/13/08 */
	      if(i == 0)
		{
		  Fin0 = f[ii][1] * fe[iii] * (1 - f5) * (1 + f[p+3][1]);
		  Fout0 = f[p+3][1] * f5 * (1 + f[ii][1]) * (1 - fe[iii]);
		}
	      if(iii == 0)
		{
		  Fin00 = f[ii][1] * fe[p+3] * (1 - f5) * (1 + f[i][1]);
		  Fout00 = f[i][1] * f5 * (1 + f[ii][1]) * fe[p+3];
		}
	      dNin = dfacel1[i][ii][iii] * Fin;
	      if (dNin != dNin || dNin > 1e100 || dNin < -1e100)
		{
		  printf("Error in 2Dpolelscat.h, #1.\n");
		  printf("df1 = %e, Fin = %e\n", dfacel1[i][ii][iii], Fin);
		  printf("f1 = %e, f2 = %e, f3 = %e, f4 = %e\n", f[i][1], f[ii][1], f[iii][1], fnew[i][ii][iii]);
		  printf("i = %d, ii = %d, iii = %d, I = %d\n", i, ii, iii, I[i][ii][iii]);
		  exit(1);
		}
	      dNout = dfacel1[i][ii][iii] * Fout;
	      Ninpolel[i] += dNin;
	      Noutpolel[i] += dNout;
	      if (Ninpolel[i] < 0 || Noutpolel[i] < 0)
		{
		  printf("error in 2Dpolelscat.h, #2\n");
		  printf("dfac %d = %e\n", i, dfacel1[i][ii][iii]);
		  printf("Fin = %e, Fout = %e\n", Fin, Fout);
		  exit(1);
		}
	      if(i == 0)
		{
		  Ninpolel[p+3] += dfacel0[ii][iii] * Fin0;
		  Noutpolel[p+3] += dfacel0[ii][iii] * Fout0;
		}
				if(iii == 0)
				  {
				    Ninpolel[i] += Fin00 * dfacel00[i][ii];
				    Noutpolel[i] += Fout00 * dfacel00[i][ii];
				  }
	    } /* loop E2 */
	} /* loop E1 */
      Ninpolel[i] *= Xp[i] * f[i][2];
      Noutpolel[i] *= Xp[i] * f[i][2];
      if(i == 0)
	{
	  Ninpolel[p+3] *= DOS[p+3];
	  Noutpolel[p+3] *= DOS[p+3];
	}
    } /* loop E */
  return(1);
}
