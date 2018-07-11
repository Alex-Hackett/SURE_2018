int initiate(int initial, double del[1000], int p, double DOS[1000], int m, double f[1000][6], double Nin[1000],
	     double Nout[1000], double N[1000], double mu, double hw, double expp, double Tm, int statype, double Ec, double indexr)
{
  switch(initial)
  {
  case 0:
    fromfile(f, del, expp, indexr, 1);
    break;
  case 1:
    uniformfvsE(del, p, DOS, m, f, N, expp, Tm, indexr, 1);
    break;
  case 2:
    boltzmannfvsE(del, p, DOS, f, N, expp, statype, indexr, 1);
    break;
  case 3:
    phononfvsE(p, b);
    gaussianfvsE(del, p, f, No, expp, Ec, indexr, 1);
    break;
  case 4:
    gaussianfvsE(del, p, f, No, expp, Ec, indexr, 1);
    break;
  case 5:
    exit(1);
  case 6:
    neareq(del, p, DOS, f, N, expp, indexr, 1);
    break;
  case 7:
    pulsefvsE(del, p, DOS, f, N, expp, statype, indexr, 1);
    break;
  case 8:
    pulsegaussianfvsE(del, p, f, No, expp, Ec, indexr, 1);
    break;
  case 9:
    pulseflatf(del, p, DOS, f, N, expp, indexr, 1);
    break;
  }
  return(1);
}
