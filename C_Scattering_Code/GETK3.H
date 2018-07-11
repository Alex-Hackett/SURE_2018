double Determinek3(double E4, double f[1000][6], int zzz)
{
  if (E4 == f[0][0])
    zzz = 0;
  else if (E4 < f[p+1][0])
    {
      while (E4 >= f[zzz][0])
	{
	  zzz++;
	  if (zzz > p + 1)
	    zzz = 0;
	}
      zzz--;
      if (zzz < 0)
	{
	  printf("error in getk3.h, zzz < 0\n");
	  printf("E4 = %e, E(0) = %e\n", E4, f[0][0]);
	  exit(1);
	}
      if(E4 < f[zzz][0] || E4 > f[zzz+1][0])
	{
	  printf("k3 error\n");
	  printf("E4 = %.20e, E(z) = %.20e, E(z-1) = %.20e, z = %d, point = %f\n", E4, f[zzz][0], f[zzz-1][0], zzz, point);
	  exit(1);
	}
      point = (E4-f[zzz][0])/(f[zzz+1][0] - f[zzz][0]);
      ff4 = kpp[zzz] + (point * delk[zzz]);
    }
  else
    ff4 = kpp[p+1];
  polD[0] = zzz;
  return(ff4);
}
