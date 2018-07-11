double Determinef4(double E4, double KK, double f[1000][6], int m, int lll,
		   int zzz, int Geo, double fpp[1000], double delk[1000], double del[1000])
{
  if ( lll == -1)
    {
      if (E4 < f[p+1][0])
	{
	  while ( E4 > f[zzz][0])
	    {
	      zzz++;
	      if (zzz > p + 1)
		zzz = 0;
	    }
	  zzz--;
	  if( (E4 < f[zzz][0]) || (E4 > f[zzz+1][0]) )
	    {
	      printf("f4 error\n");
	      printf("E4 = %.12f, E(z) = %f, z = %d, point = %f\n", E4, f[zzz][0], zzz, point);
	      exit(1);
	    }
	  point = (E4-f[zzz][0])/ del[zzz];
	  ff4 = f[zzz][1] + (point * (f[zzz+1][1] - f[zzz][1]));
	}
      else
	ff4 = 0;
    }
  else
    {
      if ( lll == m)
	point = 0;
      else
	{
	  if ( Geo == 5)
	    point = (KK - kpp[lll]) / delk[lll];
	  else
	    point = (E4-f[lll][0])/ del[lll];
	}
      ff4 = f[lll][1] + (point * (f[lll+1][1] - f[lll][1]));
    }
  return(ff4);
}
