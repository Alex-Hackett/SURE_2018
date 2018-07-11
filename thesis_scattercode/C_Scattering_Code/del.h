void delSet(double del[1000], double expp, int delc, double f[1000][6], double indexr)
{
  int p1; /* p - 20 */
  switch(delc)
    {
    case 0:
      /* ----- uniform change in occupation number del ----- */
      uu = log(1 - exp(-2*g));
      fff0 = 1 / (exp(-uu) - 1);
      fffm = 1 / (exp(maxkT - uu) - 1);
      ch = 1 - exp(log(fffm/fff0)/p);
      offs = log(1 / (1 + ch));
      for (i = 0; i <= p+1; i++)
	{

  del[i] = kb * T * ( offs + log((1 + fff0 + ch * fff0)/(1 + fff0)));
 fff0 = fff0 * (1 - ch);
	}
      break;
    case 1:
      /* ----- exp(i) del ---- */
      expp = 1;
      for (j = 0; j <= p; j++)
	expp *= e;
      /*printf("%e\n", expp);*/
      dels = (expp - e) / ((e - 1) * maxkT);
      expp = 1;
      for (i = 0; i <= p + 1; i++)
	{
	  expp *= e;
	  del[i] = kb * TT * expp / dels;
	}
      break;
    case 2:
      /* ----- quadratic del ----- */
      dels = (( p * p * p ) + ( 1.5 * p * p ) + ( 0.5 * p )) / ( 3 * maxkT );
      for (i = 0; i <= p + 1; i++)
	del[i] = kb * T * (i+1) * (i+1) / dels;
      break;
    case 3:
      /* ----- cubic del ----- */
      dels = ( p * p ) * ( p * p + 2 * p + 1) / ( 4 * maxkT);
      for (i = 0; i <= p + 1; i++)
	del[i] = kb * T * (i+1) * (i+1) * (i+1) / dels;
      break;
    case 4:
      /* ----- Linear del ----- */
      for (i = 0; i <= p; i++)
	{
	  if ( i <= 50)
	    del[i] = kb * T / 100;
	  else
	    del[i] = kb * T / 100;
	}
      break;
    case 5:
      /* ----- polariton del ----- */
      /* dk quadratic in k */
      dels = (( p * p * p ) + ( 1.5 * p * p ) + ( 0.5 * p )) / (3*maxkp);
      if(disp == 1)
	{
	  kpp[0] = 0;
	  f[0][0] = 0;
	  for(i = 1; i <= p + 2; i++)
	    {
	      kpp[i] = kpp[i-1] + delk[i-1];
	      f[i][0] = hb * hb * eC * kpp[i] * kpp[i] / (2 * em * (Me + Mh));
	    }
	}
      if(disp == 2 || disp == 3)
	{
	  kpp[0] = 0;
	  Exk[0] = Ex0;
	  Eck[0] = hb * vc * kcz / indexr;
	  Elp1[0] = Exk[0] * Exk[0] + Eck[0] * Eck[0] + gs * gs;
	  f[0][0] = sqrt((Elp1[0] - sqrt(Elp1[0] * Elp1[0] - 4 * Exk[0] * Exk[0] * Eck[0] * Eck[0])) / 2);
	  f[0][4] = f[0][0];
	  delk[0] = 1 / dels;
	  for (i = 1; i <= p + 2; i++)
	    {
	      delk[i] = (1+i)*(1+i) / dels;
	      kpp[i] = kpp[i-1] + delk[i-1];
	      Exk[i] = Ex0 + (hb * hb * kpp[i] * kpp[i] * eC) / (2 * em * (Me + Mh));
	      Eck[i] = hb * (vc / indexr) * sqrt((kcz * kcz) + (kpp[i] * kpp[i]));
	      Elp1[i] = Exk[i] * Exk[i] + Eck[i] * Eck[i] + gs * gs;
	      f[i][0] = sqrt((Elp1[i] - sqrt(Elp1[i] * Elp1[i] - 4 * Exk[i] * Exk[i] * Eck[i] * Eck[i])) / 2);
	      f[i][4] = f[i][0];
	      f[i][5] = f[i][0];
	    }

  }
      for (i = 0; i <= p + 1; i++)
	{
	  del[i] = f[i+1][0] - f[i][0];
	  dEdki[i] = getdEdk(f[i][0], kpp[i], disp, indexr, i, 1);
	}
      /* for free electrons */
      DOSe[0] = 2 * em * Me / (pi * hb * hb * eC);
      Ne[0] = 0;
      qo[0] = 0; /*eC / (2 * einf * eo);*/ /* added (2 * pi) on 3/27/08 */
      for(i = 0; i <= p+1; i++)
	{
	  Eek[i] = hb * hb * kpp[i] * kpp[i] * eC / (2 * em * Me);
	  fe[i] = 1 / (exp( (Eek[i] / (kb * TTT[0])) - (mue[0]) ) + 1);
	}
      for(i = 0; i <= p; i++)
	Ne[0] += fe[i] * DOSe[0] * (Eek[i+1] - Eek[i]);
      /*printf("Ne = %e\n", Ne[0]);
	exit(1); */
      break;
    }
}
