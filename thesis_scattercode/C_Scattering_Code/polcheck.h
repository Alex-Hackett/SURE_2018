void polprint(double kpp[1000], long double f[1000][4], double DOS[1000], double N[1000], double indexr, int y)
{
  double dddd;
  printf("i k E DOS dE/dk N\n");
  for(i=0;i<=p;i++)
    {
      printf("%d %.6e ", i, kpp[i]);
      printf("%e ", f[i][0]);
      printf("%f ", DOS[i]);
      dddd = getdEdk(f[i][0], kpp[i], disp, indexr, i, y);
      printf("%e ", dddd);
      printf("%e\n", N[i]);
    }
}
