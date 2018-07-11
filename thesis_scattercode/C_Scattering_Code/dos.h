/* number of states / unit volume / differential energy */
double DOSf(double DOS[1000], double f[1000][6], double indexr, int y)
{
  DOS[p+3] = dg / Sa; /* ground state density of states */
  for (jj = 0; jj <= p + 1; jj++)
    {
      if (Geo == 3)
	DOS[jj] = dg * em * Me * sqrt(em * Me * f[jj][0] / (2 * eC)) / (pi * pi * hb * hb * hb * eC);
      if (Geo == 2)
	DOS[jj] = dg * em * Me / (pi * hb * hb * eC);
      /*if (Geo == 1)
	{ DOS1 = 1 / sqrt(f[jj][0]);}*/
      if ( Geo == 4)
	DOS[jj] = (f[jj][0] / hw);
      if ( Geo == 5)
	/* DOS[jj] = kpp[jj] / (2 * pi); */
	/* DOS[jj] = (dg / (4 * pi * hb * hb)) * f[jj][0] * (Elp1[jj] - (2 * f[jj][0] * f[jj][0])) /
(((Exk[jj] / (em * (Me + Mh))) * ((Eck[jj] * Eck[jj]) - (f[jj][0] * f[jj][0]))) +
((vc * vc / (index * index)) * (Exk[jj] * Exk[jj] - f[jj][0] * f[jj][0]))); */
	DOS[jj] = (dg / (4 * pi)) / getdEdk(f[jj][0], kpp[jj], disp, indexr, 0, y);
    }
  return(1);
}
