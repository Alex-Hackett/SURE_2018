int D2polscat(double f[1000][6], double Ninpol[1000], double Noutpol[1000], int m,
	      double del[1000], int zzz, int Geo, double kpp[1000], double delk[1000], double index, int disp, int y)
{
  double Fin0, Fout0, Fin00, Fout00, dEdk3, den00, rado, dEdki2, dEdki3, kxk;
  nnn=0;
  /* include a loop for each integration variable in this list */
  /* E, E1, E2 */
  Allconst = 1.2 * 36 * EB * EB * aB * aB * aB * aB / (hb * 32 * pi * pi * pi);
  den00 = 0;
  if (y == 1 || disp == 3)
    {
      for(i = 0; i <= p; i++)
	{
	  if(i == 0)
	    kxk = 1;
	  else
	    kxk = kpp[i] * kpp[i];
	  for(ii = 0; ii <= p; ii++)
	    {
	      Vp[i][ii] = 0;
	      if(i == 0 && ii != 0)
		kxk1 = kpp[ii];
	      else if(i != 0 && ii == 0)
		kxk1 = kpp[i];
	      else if(i == 0 && ii == 0)
		kxk1 = 1;
	      else
		kxk1 = kpp[i] * kpp[ii];
	      dEdki2 = dEdki[ii]; /* 9/19/07 change */ /* getdEdk(f[ii][0], kpp[ii], disp, index, 0, y); */
	      for(iii = 0; iii <= p; iii++)
		{
		  E3 = f[iii][0] + f[ii][0] - f[i][0];
		  dEdki3 = dEdki[iii]; /* 9/19/07 change */ /*getdEdk(f[iii][0], kpp[iii], disp, index, 0, y); */
		  if (E3 > f[0][0] && E3 < f[p+1][0])
		    {
		      k3 = Determinek3(E3, f, 0);
		      I[i][ii][iii] = polD[0];
		      k2k3 = kpp[iii] * k3;
		      dEdk3 = getdEdk(E3, k3, disp, index, 0, y);
		      fnew[i][ii][iii] = (k3 - kpp[I[i][ii][iii]]) / delk[I[i][ii][iii]];
		      /*f[I[i][ii][iii]][1] + ( (f[1 + I[i][ii][iii]][1] - f[I[i][ii][iii]][1]) * ((k3 - kpp[I[i][ii][iii]]) / delk[I[i][ii][iii]]) ); */
		      /*
if(I[i][ii][iii] > 0)
fnew[i][ii][iii] = findf(I[i][ii][iii], f, delk, (k3 - kpp[I[i][ii][iii]]) / delk[I[i][ii][iii]]);
else
fnew[i][ii][iiii] = f[0][1] + ((f[1][1] - f[0][1]) * (k3 / delk[0]));
		      */
		      if (i == 0)
			{
			  E40[ii][iii] = f[iii][0] + f[ii][0] - f[0][0];
			  k40[ii][iii] = Determinek3(E40[ii][iii], f, zzz);
			  if( (k40[ii][iii] < (kpp[ii] + kpp[iii])) && ((k40[ii][iii]*k40[ii][iii]) > (kpp[ii]-kpp[iii])*(kpp[ii]-kpp[iii])) )
			    {
  dEdk40[ii][iii] = getdEdk(E40[ii][iii], k40[ii][iii], disp, index, zzz, y);
 dfac0[ii][iii] = (0.5 * Xp[ii] * Xp[iii] * Xp[I[i][ii][iii]]) * deldel[ii][iii] * 2 * pi /
   (dEdki2 * dEdki3 * dEdk40[ii][iii] * (sqrt((2 * kpp[iii] * kpp[iii] * kpp[ii] * kpp[ii]) - (k40[ii][iii]*k40[ii][iii]*k40[ii][iii]*k40[ii][iii]) +
					      (2*k40[ii][iii]*k40[ii][iii]*(kpp[ii]*kpp[ii] + kpp[iii]*kpp[iii])) - (kpp[ii]*kpp[ii]*kpp[ii]*kpp[ii]) -
					      (kpp[iii]*kpp[iii]*kpp[iii]*kpp[iii])) ) );
			    }
			  else dfac0[ii][iii] = 0;
			}
		      if (iii == 0)
			{
			  E400[i][ii] = f[i][0] + f[ii][0] - f[0][0];
			  k400[i][ii] = Determinek3(E400[i][ii], f, zzz);
			  if( (k400[i][ii] < (kpp[i] + kpp[ii])) && ((k400[i][ii]*k400[i][ii]) > (kpp[i]-kpp[ii])*(kpp[i]-kpp[ii])) )
			    {
			      dEdk400[i][ii] = getdEdk(E400[i][ii], k400[i][ii], disp, index, zzz, y);
			      den00 = (2 * kpp[i] * kpp[i] * kpp[ii] * kpp[ii]) - (k400[i][ii] * k400[i][ii]
										   * k400[i][ii] * k400[i][ii]) + (2 * k400[i][ii] * k400[i][ii] * (kpp[ii]*kpp[ii] + kpp[i]*kpp[i])) -
				(kpp[ii]*kpp[ii]*kpp[ii]*kpp[ii]) - (kpp[i]*kpp[i]*kpp[i]*kpp[i]);
			      dfac00[i][ii] = 2 * (Xp[ii] * 0.5 * Xp[I[i][ii][iii]]) * del[ii] * 8 * pi * pi / (D2A * dEdki2 * dEdk400[i][ii] * sqrt(den00) );
			    }
			  else dfac00[i][ii] = 0;
			}
		      radp = (kxk - kpp[iii] * kpp[iii] - k2k3) / kxk1;
		      radn = (kxk - kpp[iii] * kpp[iii] + k2k3) / kxk1;
		      if(radp < -1)
			thetamax = pi;
		      else if(radp > 1)
			thetamax = 0;
		      else
			thetamax = acos(radp);
		      if (radn > 1)
			thetamin = 0;
		      else if (radn < -1)
			thetamin = pi;
		      else
			thetamin = acos(radn);
		      if (thetamax - thetamin > 1e-3) /* 1e-13 */
			{
			  dif = thetamax - thetamin;
			  ttheta = thetamax + thetamin;
			  GQ1 = 0;
			  for(j = 1; j <= GQp; j++)
			    {
			      rado = cos(((dif * yy[j]) + ttheta) / 2);
			      GQ1 += ww[j] / sqrt( (radp - rado) * (-radn + rado) );
			    }
			  GQ1 *= dif / (2 * kxk1);
			  dfac1[i][ii][iii] = (Xp[ii] * Xp[iii] * Xp[I[i][ii][iii]]) * GQ1 * deldel[ii][iii] / (dEdki2 * dEdki3 * dEdk3);
			}
		      else
			dfac1[i][ii][iii] = 0;
		    }
		  else
		    {
		      dfac1[i][ii][iii] = 0;
		      if(i == 0) dfac0[ii][iii] = 0;
		    }
		  if (dfac0[ii][iii] != dfac0[ii][iii] || dfac0[ii][iii] > 1e50 || dfac00[i][ii] != dfac00[i][ii] ||
dfac00[i][ii] > 1e50 || dfac0[ii][iii] < -1e50 || dfac00[i][ii] < -1e50 ||
		      dfac1[i][ii][iii] < - 1e50 || dfac1[i][ii][iii] > 1e50 || dfac1[i][ii][iii] != dfac1[i][ii][iii])
		    {
		      printf("error in 2Dpolscat.h\n");
		      printf("i = %d, ii = %d, iii = %d, I = %d\n", i, ii, iii, I[i][ii][iii]);
		      printf("df1 = %e df0 = %e df00 = %e \n", dfac1[i][ii][iii], dfac0[ii][iii], dfac00[i][ii]);
		      printf("GQ = %e, dEdk2 =%e, dEdk3 = %e, dEdk4 = %e, dif = %e\n", GQ1, dEdki2, dEdki3, dEdk3, dif);
		      printf("E400 = %e k400 = %e den00 = %e \n", E400[i][ii], k400[i][ii], den00);
		      printf("dEdk400 = %e, kxk1 = %e\n", dEdk400[i][ii], kxk1);
		      exit(1);
		    }
		  if(den00 < 0)
		    {
		      printf("den00 error in 2Dpolscat.h, k = %e, k1 = %e,k400 = %e, den00 = %e\n", kpp[i], kpp[ii], k400[i][ii], den00);
		      printf("i = %d, ii = %d\n", i, ii);
		      exit(1);
		    }
		}
	    }
	}
    }
  Ninpol[p+3] = 0;
  Noutpol[p+3] = 0;
  /* LOOP E */
  for (i = 0; i <= p; i++)
    {
      Ninpol[i] = 0;
      Noutpol[i] = 0;
      /* LOOP E1 */
      for (ii = 0; ii <= p; ii++)
	{
	  /* LOOP E2 */
	  for (iii = 0; iii <= p; iii++)
	    {
	      f5 = f[I[i][ii][iii]][1] + ((f[I[i][ii][iii]+1][1] - f[I[i][ii][iii]][1]) * fnew[i][ii][iii]);
	      Fin = f[ii][1] * f[iii][1] * ( 1 + f5 ) * (1 + f[i][1]);
	      Fout = f[i][1] * f5 * (1 + f[ii][1]) * (1 + f[iii][1]);
	      /* uncomment the following two lines to delete Bose effects
Fin = f[ii][1] * f[iii][1];
Fout = f5 * f[i][1];*/
	      /* printf("Fin = %e Fout = %e\n", Fin, Fout); */
	      if(i == 0)
		{
		  Fin0 = f[ii][1] * f[iii][1] * (1 + f5) * (1 + f[p+3][1]);
		  Fout0 = f[p+3][1] * f5 * (1 + f[ii][1]) * (1 + f[iii][1]);
		}
	      if(iii == 0)
		{
		  Fin00 = f[ii][1] * f[p+3][1] * (1 + f5) * (1 + f[i][1]);
		  Fout00 = f[i][1] * f5 * (1 + f[ii][1]) * (1 + f[p+3][1]);
		}
	      dNin = dfac1[i][ii][iii] * Fin;
	      if (dNin != dNin || dNin > 1e100 || dNin < -1e100)
		{
		  printf("df1 = %e, Fin = %e\n", dfac1[i][ii][iii], Fin);
		  printf("f1 = %e, f2 = %e, f3 = %e, f4 = %e\n", f[i][1], f[ii][1], f[iii][1], fnew[i][ii][iii]);
		  printf("i = %d, ii = %d, iii = %d, I = %d\n", i, ii, iii, I[i][ii][iii]);
		  exit(1);
		}
	      dNout = dfac1[i][ii][iii] * Fout;
	      Ninpol[i] += dNin;
	      Noutpol[i] += dNout;
	      if(i == 0)
		{
		  Ninpol[p+3] += dfac0[ii][iii] * Fin0;
		  Noutpol[p+3] += dfac0[ii][iii] * Fout0;
		}
	      if(iii == 0)
		{
		  Ninpol[i] += Fin00 * dfac00[i][ii];
		  Noutpol[i] += Fout00 * dfac00[i][ii];
		}
	    } /* loop E2 */
	} /* loop E1 */
      Ninpol[i] *= Allconst * Xp[i] * f[i][2];
      Noutpol[i] *= Allconst * Xp[i] * f[i][2];
  if(i == 0)
    {
      Ninpol[p+3] *= DOS[p+3];
      Noutpol[p+3] *= DOS[p+3];
    }
    } /* loop E */
  return(1);
}
