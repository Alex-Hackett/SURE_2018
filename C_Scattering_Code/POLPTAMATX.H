int polelMatx(double f[1000][6], double Ninpolel[1000], double Noutpolel[1000], int m,
	      double del[1000], int zzz, int Geo, double kpp[1000], double delk[1000], double index, int disp, int y)
{
  double Fin0, Fout0, Fin00, Fout00, dEdk3, den00, rado, dEdki2, dEdki3, kmk1, kmk12, kxk, k1xk1, qo2;
  nnn=0;
  /* include a loop for each integration variable in this list */
  /* E, E1, E2 */
  /* Allconst1 = 6 * 196 * EB * EB * aB * aB / (hb * 32 * pi * pi * pi); */ /* eC * eC / (128 * hb * eo * eo * einf * einf * pi * pi * pi); */
  Allconst1 = 1000 * EB * EB * aB * aB / (hb * pi * pi * pi * pi * pi * pi);
  qo2 = qo[0] * qo[0] * Ne[0] * Ne[0] / (kb * kb * TTT[0] * TTT[0]);
  /*
Be = Me / (Me + Mh1);
Bh = Mh1 / (Me + Mh1);
  */
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
	if(i == 0 && ii != 0)
	  kxk1 = kpp[ii];
	else if(i != 0 && ii == 0)
	  kxk1 = kpp[i];
	else if(i == 0 && ii == 0)
	  kxk1 = 1;
	else kxk1 = kpp[i] * kpp[ii];
	if(ii == 0)
	  k1xk1 = 1;
	else
	  k1xk1 = kpp[ii] * kpp[ii];
	dEdki2 = dEdki[ii];
	for(iii = 0; iii <= p; iii++)
	  {
	    E3 = f[ii][0] - f[i][0] + Eek[iii];
	    dEdki3 = dg / (4 * pi * DOSe[0]); /* DOS is a constant for the electron */
	    if (E3 > Eek[0] && E3 < Eek[p+1])
	      {
		k3 = sqrt( (2 * em * Me * (E3 - Eek[0]) ) / (hb * hb * eC) );
		jj = 0;
		while (kpp[jj] <= k3)
		  jj++;
		I[i][ii][iii] = jj--;
		k2k3 = kpp[iii] * k3;
		/* if (kpp[I[i][ii][iii]] > k3)
{
printf("error in polelMatx.h with I[i][ii][iii]\n");
printf("i = %d, ii = %d, iii = %d, I = %d\n", i, ii, iii, I[i][ii][iii]);
printf("kpp - k3 = %e\n", kpp[I[i][ii][iii]] - k3);
printf("k3 = %e\n", k3);
exit(1);
} */
		dEdk3 = dg / (4 * pi * DOSe[0]);
		fnew[i][ii][iii] = (k3 - kpp[I[i][ii][iii]]) / delk[I[i][ii][iii]];
		/* take out k = 0 ’point’ if (i == 0)
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
den00 = (2 * kpp[i] * kpp[i] * kpp[ii] * kpp[ii]) - (k400[i][ii] * k400[i][ii] * k400[i][ii] * k400[i][ii]) +
  (2 * k400[i][ii] * k400[i][ii] * (kpp[ii]*kpp[ii] + kpp[i]*kpp[i])) - (kpp[ii]*kpp[ii]*kpp[ii]*kpp[ii]) - (kpp[i]*kpp[i]*kpp[i]*kpp[i]);
 dfac00[i][ii] = 2 * (Xp[ii] * 0.5 * Xp[I[i][ii][iii]]) * del[ii] * 8 * pi * pi / (D2A * dEdki2 * dEdk400[i][ii] * sqrt(den00) );
	      }
 else dfac00[i][ii] = 0;
  } */
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
			kmk12 = kxk + k1xk1 - (2 * kxk1 * rado);
			/* if(kmk12 < 0) kmk12 = 0; */
			/*
kmk1 = sqrt(kmk12);
Bhka = Bh * kmk1 * aB / 2; */ /* variiertes sigma_h */
			Beka = kxk * aB * aB / 4;
			Bcka = kmk12 * aB * aB;
			Bhka3 = (2 + Beka) * (2 + Beka) * (2 + Beka);
			Bhka34 = (2 + Bcka + Beka) * (2 + Bcka + Beka) * (2 + Bcka + Beka);
			Beka3 = (k3 * k3) + (1 / (aB * aB)) + kmk12;
			Behka = Beka3 + qo2 + (2 * sqrt(qo2 * Beka3));
			/*
Bcka3 = 1 / ( (1 + (Bcka * Bcka)) * (1 + (Bcka * Bcka)) * (1 + (Bcka * Bcka)) );
Bhka3Beka3 = Bhka3 + Beka3 - Behka;
if (fabs(Bhka3Beka3) < 1e-15) Bhka3Beka3 = 0;
Bcka3Bhka34 = Bcka3 + Bhka34 - (2 * sqrt(Bcka3 * Bhka34));
if (fabs(Bcka3Bhka34) < 1e-15) Bcka3Bhka34 = 0;
			*/
			if(ScreenType = 1)
			  GQ1 += ww[j] / ( ((1 / (kmk12 * aB * aB)) *
					    (qo2 + kmk12 + (2 * sqrt(qo2 * kmk12)))) * sqrt( (radp - rado) * (-radn + rado) ) );
			else
			  GQ1 += ww[j] / ( Bhka3 * Bhka34 * Behka * sqrt( (radp - rado) * (-radn + rado) ) );
			/* GQ1 += (ww[j] / sqrt( (radp - rado) * (-radn + rado) )) *
			   ( (Bhka3Beka3 / (qo2 + kmk12 + (2 * sqrt(qo2 * kmk12)))) + (80 * exp(-2 * aB * aB * (qo2 + kmk12 + (2 * sqrt(qo2 * kmk12))))) ); */
			/*
GQ1 += ww[j] * ( (Bhka3Beka3 + (80 * exp(-2 * kmk12 * aB * aB)) ) /
( qo2 + kmk12 + (2 * sqrt(qo2 * kmk12))) ) / sqrt( (radp - rado) * (-radn + rado) );*/
			/* 3/21/08 removed (16 * Bcka3Bhka34 / ((1/(aB * aB)) + qo2)) ) with ( 225 * exp(-(kmk12 * aB * aB)) / qo2) */
			/* changed 225 to 80 = (14 * 2 / pi)^2 */
		      }
		    GQ1 *= dif /(2 * kxk1);
		    dfacel1[i][ii][iii] = Allconst1 * Xp[ii] * GQ1 * del[ii] * (Eek[iii+1] - Eek[iii]) / (dEdki2 * dEdki3 * dEdk3);
		  }
		else
		  dfacel1[i][ii][iii] = 0;
	      }
	    else
	      {
		dfacel1[i][ii][iii] = 0;
  if(i == 0) dfacel0[ii][iii] = 0;
	      }
	    if (dfacel0[ii][iii] != dfacel0[ii][iii] || dfacel0[ii][iii] > 1e60 || dfacel00[i][ii] != dfacel00[i][ii]
|| dfacel00[i][ii] > 1e60 || dfacel0[ii][iii] < -1e60 || dfacel00[i][ii] < -1e60 ||
		dfacel1[i][ii][iii] < - 1e60 || dfacel1[i][ii][iii] > 1e60 || dfacel1[i][ii][iii] != dfacel1[i][ii][iii])
	      {
		printf("error (2) in polelMatx.h\n");
		printf("i = %d, ii = %d, iii = %d, I = %d\n", i, ii, iii, I[i][ii][iii]);
		printf("df1 = %e df0 = %e df00 = %e \n", dfacel1[i][ii][iii], dfacel0[ii][iii], dfacel00[i][ii]);
		printf("GQ = %e, dEdk2 =%e, dEdk3 = %e, dEdk4 = %e, dif = %e\n", GQ1, dEdki2, dEdki3, dEdk3, dif);
		printf("E400 = %e k400 = %e den00 = %e \n", E400[i][ii], k400[i][ii], den00);
		printf("dEdk400 = %e, kxk1 = %e\n", dEdk400[i][ii], kxk1);
		exit(1);
	      }
	    if(den00 < 0)
	      {
		printf("den00 error in polelMatx.h, k = %e, k1 = %e,k400 = %e, den00 = %e\n", kpp[i], kpp[ii], k400[i][ii], den00);
		printf("i = %d, ii = %d\n", i, ii);
		exit(1);
	      }
	    /* printf("dfac = %e\n", dfacel1[i][ii][iii]); */
	  }
      }
  }
    }
  Ninpolel[p+3] = 0;
  Noutpolel[p+3] = 0;
  return(1);
}
