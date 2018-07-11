int polfrac()
{
  for(i = 0; i <= p; i++)
    {
      if(disp == 1)
	Xp[i] = 1;
      if(disp == 2 || disp == 3)
	Xp[i] = (((Eck[i] - Exk[i]) / sqrt(Elp1[i] - 2 * Eck[i] * Exk[i])) + 1) / 2;
      ki[i] = (kpp[i] + kpp[i+1]) / 2;
      /* printf("Xp%d = %e\n", i, Xp[i]); */
    }
  return(1);
}
