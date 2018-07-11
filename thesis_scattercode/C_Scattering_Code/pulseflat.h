int pulseflatf(double del[1000], int p, double DOS[1000], double f[1000][6], double N[1000], double expp, double indexr, int y)
{
  delSet(del, expp, delc, f, indexr);
 for(i = 0; i <= p; i++)
   {
     for(j = 0; j <= p; j++)
       deldel[i][j] = del[i] * del[j];
   }
 if(Geo != 5)
   {
     if (Geo == 4)
       f[0][0] = hw;
     else
       f[0][0] = 0;
     for (j = 1; j <= p + 1; j++)
       f[j][0] = f[j-1][0] + del[j-1];
   }
 for (j = 0; j <= p; j++)
   f[j][1] = (4E-15/pulseT)*g;
 f[p+1][1] = f[p][1] * exp(-del[p] / (kb * T));
 f[p+2][1] = f[p+1][1] * exp(-del[p+1] / (kb * T));
 DOSf(DOS, f, indexr, y);
 for (j = 0; j <= p+1; j++)
   {
     f[j][2] = (DOS[j] + DOS[j+1]) * del[j] / 2;
     N[j] = f[j][1] * f[j][2];
   }
 N[p+3] = (4E-15 / pulseT) * DOS[p+3] * g;
 f[p+3][1] = N[p+3] / DOS[p+3];
 return(1);
}
