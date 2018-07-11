double getdEdk(double EE, double KK, int disp, double indexr, int zzz, int y)
{
  int up, down;
  if(disp == 1)
    {
      EEKK = hb * hb * eC / (2 * em * (Me + Mh));
    }
  if(disp == 2 || (disp == 3 && y == 1))
    {
      Exkt = Ex0 + (hb * hb * KK * KK * eC/ (2 * em * (Me + Mh)));
      Eckt = hb * vc * sqrt(kcz * kcz + KK * KK) / indexr;
      /* E2EE = (Exkt * Exkt) + (Eck * Eck) + (gs * gs) - (2 * EE * EE);
	 E4EE = (Exkt / (Me + Mh)) + (vc * vc / (indexr * indexr));*/
      /* E3EE = 2 * Exkt * hb * hb * KK / (em * (Me + Mh));
	 dE/dk^2 ---> */
      EEKK = hb * hb * ((Exkt / (em * (Me + Mh))) * eC * (Eckt * Eckt - EE * EE) +
			(vc * vc / (indexr * indexr)) * (Exkt * Exkt - EE * EE)) / (2 * EE * ((Exkt * Exkt) + (Eckt * Eckt) + (gs * gs) - (2 * EE * EE)));
      /*
dE/dk --> EEKK = ((2 * E3EE) / (4 * EE * E2EE)) * ((Eck * Eck) - (EE * EE));
      */
      /* EEKK = hb * hb * KK * ( E4EE - (E2EE * ( (Exkt * Exkt + Eckt * Eckt + gs * gs) * E4EE - 2 * (Exkt * Exkt * Eckt / (Me + Mh) +

(vc * vc * Exkt * Exkt / (indexr * indexr))) ) / ( 2 * EE);*/
      if (EEKK < 1e-22)
	EEKK = 1e-22;
    }
    if(disp == 3 && y != 1)
      {
	Determinek3(EE, f, zzz);
	up = polD[0] + 1;
	down = polD[0] - 1;
	EEKK = (f[up][0] - f[down][0]) / ((kpp[up] * kpp[up]) - (kpp[down] * kpp[down]));
      }
  return(EEKK);
}
