double updatef(double Ntotal, double uprate[1], double x,
	       double dN[1000], double N[1000], double f[1000][6], double Nin[1000], double
	       Nout[1000], int Geo, double kpp[1000], double Stime[200], double taut, double Avetau[0])
{
  double Nnewtau, Nnewtau1; /* for calculating average total lifetime */
  double deltaN;
  deltaN = 0;
  Ntotalin = 0;
  Ntotalout = 0;
  for (j = 0; j <= p; j++)
    {
      Ntotalin += Nin[j];
      Ntotalout += Nout[j];
    }
  /* printf("Ntotalin %e Ntotalout %e Ntotal %e uprate %e\n", Ntotalin, Ntotalout, Ntotal, uprate[0]); */
  x = uprate[0] * Ntotal;
  for (j = 0; j <= p; j++)
    {
      f[j][3] = (Nin[j]/Ntotalin) - (Nout[j]/Ntotalout);
      /* printf ("N%d = %e changes by %e with Nin%d = %e and Nout%d = %e\n", j, N[j], x * f[j][3], j, Nin[j], j, Nout[j]); */
      N[j] += x * f[j][3];
      deltaN += f[j][3];
      if(N[j] < 0)
	{
	  N[j] -= x * f[j][3];
	  N[j] /= 2;
	}
    }
  deltaNAll[0] += deltaN;
  if(N[p+3] < 0) N[p+3] = 0;
  Stime[p+3] = 0;
  for (j = 0; j <= p; j++)
    {
      if(Nin[j] != Nout[j]*Ntotalin/Ntotalout)
	{
  Stime[j] = f[j][3] / (Nin[j] - (Nout[j]*Ntotalin/Ntotalout));
 Stime[p+3] += fabs(Stime[j]);
	}
    }
  Stime[p+3] *= x / (p+1);
  /* Stime[p+3] = 0;
for (j = 0; j <= p; j++)
{
if(Nin[j] != Nout[j])
{
Stime[p+3] += fabs(f[j][3]);
Stime[j] = f[j][3] / (Nin[j] - Nout[j]);
}
}
Stime[p+3] *= x / fabs(Ntotalin - Ntotalout);*/
  /* system losses */
  Nnewtau = 0;
  Nnewtau1 = 0;
  for(j = 0; j <= p; j++)
    {
      if(f[j][0] - f[0][0] < 0.0035) Nnewtau1 += N[j];
      N[j] -= Stime[p+3] * N[j] * ( (Itaux * Xp[j]) + (Itaup * (1 - Xp[j])) );
      /*N[j] -= Stime[p+3] * Itau * N[j] / Xp[j];*/ /* *Xp[j]*Xp[j]*Xp[j]*Xp[j]*Xp[j]*Xp[j]);*/ /*Xp? see 7/10/06 in notes.
	Should depend on photon fraction. Why Xp^7?*/
      f[j][1] = N[j] / f[j][2];
      if(f[j][0] - f[0][0] < 0.0035) Nnewtau += N[j];
    }
  Avetau[0] = (Stime[p+3] * Nnewtau1) / (Nnewtau1 - Nnewtau);
  /*printf("particle loss = %e \n", Stime[p+3] * Itau);*/
  N[p+3] -= Stime[p+3] * ( (Itaux * Xp[j]) + (Itaup * ( 1 - Xp[j])) ) * N[p+3]; /* 2^7 = 128 */
  if(N[p+3] < 0) N[p+3] = 0;
  f[p+3][1] = f[0][1]; /*N[p+3] / DOS[p+3];*/
  if(initial == 7 && (taut + Stime[p+3]) < pulseT) /* puts in a boltzmann distribution at a given temperature */
    {
      for(j = 0; j <= p; j++)
	{
	  f[j][1] += (50 * Stime[p+3] / pulseT) / (exp( ( ((kpp[j] * kpp[j] * hb * hb * eC) / (2 * em * (Me + Mh))) - mu) / (kb * TT) ) - 1);
	  /* set to put in thermalized excitons */
	  /* (exp( (f[j][0] - mu - f[0][0]) / (kb * TT) ) - 1); */
	  N[j] = f[j][1] * f[j][2];
	}
      N[p+3] += ( 50 * Stime[p+3] / pulseT) / (exp( - mu / (kb * TT) ) - 1);
      f[p+3][1] = N[p+3] / DOS[p+3];
    }
  if(initial == 8 && (taut + Stime[p+3]) < pulseT) /* puts in a gaussian distribution at a specific energy */
    {
      for(j = 0; j <= p; j++)
	{
	  f[j][1] += (Stime[p+3] / pulseT) * No[j] * 2 / (DOS[j] + DOS[j+1]);
	  N[j] = f[j][1] * f[j][2];
	}
      f[p+3][1] += (Stime[p+3] / pulseT) * No[0] * 2 / (DOS[0] + DOS[1]);
      N[p+3] = f[p+3][1] * DOS[p+3];
    }
  if(initial == 9 && (taut + Stime[p+3]) < pulseT) /* puts in a flat distribution in k space */
    {
      for(j = 0; j <= p; j++)
	{
	  f[j][1] += (Stime[p+3]/ pulseT) * g;
	  N[j] = f[j][1] * f[j][2];
	}
      f[p+3][1] += (Stime[p+3] /pulseT) * g;
      N[p+3] = f[p+3][1] * DOS[p+3];
    }
  f[p+1][1] = f[p][1] * exp( -del[p] / (kb * T) );
  f[p+2][1] = f[p+1][1] * exp( -del[p+1] / (kb * T) );
  /*printf("f(0) = %e\n", f[0][1]);
    exit(1); */
  return(x);
}
