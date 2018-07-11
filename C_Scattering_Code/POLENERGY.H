double polE(double kpp1, double index, int disp)
{
  if(disp == 1)
    Elpt = hb * hb * kpp1 * kpp1 * eC / (2 * em * (Me + Mh));
  if(disp == 2)
    {
      Exkt = Ex0 + (hb * hb * kpp1 * kpp1 * eC) / (2 * em * (Me + Mh));
      Eckt = hb * vc * sqrt((kcz * kcz) + (kpp1 * kpp1)) / index;
      Elp1t = Exkt * Exkt + Eckt * Eckt + gs * gs;
      Elpt = sqrt((Elp1t - sqrt((Elp1t * Elp1t) - (4 * Exkt * Exkt * Eckt * Eckt))) / 2);
    }
  return(Elpt);
}
