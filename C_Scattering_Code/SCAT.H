int scatter(int scattype, double f[1000][6], double Nin[1000], double Nout[1000],
	    int m, int zzz, double del[1000], int Geo, double kpp[1000], double delk[1000], int y, double indexr, double kcz, int disp)
{
  switch(scattype)
    {
    case 1:
      D3bosescat(f, Nin, Nout, m, Geo, kpp, delk); /* boson-boson scattering */
      break;
    case 2:
      D2bosescat(f, Ninbos, Noutbos, m, del, zzz, Geo, kpp, delk);
      for(i = 0; i <= p; i++)
	{
	  Nin[i] = Ninbos[i];
	  Nout[i] = Noutbos[i];
	}
      Nin[p+3] = Ninbos[p+3];
      Nout[p+3] = Noutbos[p+3];
      break;
    case 3:
      D3xpEexchange(f, b, Nin, Nout, del, zzz, Geo, kpp, delk);
      break;
    case 4:
      /* D2xpEexchange(f, b, Nin, Nout, del, zzz, Geo, kpp, delk);
break;
      */
      D2bospscat(f, Ninph, Noutph, m, del, zzz, Geo, kpp, delk, y);
      for(i = 0; i <= p; i++)
	{
	  Nin[i] = Ninph[i];
	  Nout[i] = Noutph[i];
	}
      Nin[p+3] = Ninph[p+3];
      Nout[p+3] = Noutph[p+3];
      break;
    case 5:
      GQtest(f, Nin, Nout, m, del, zzz, Geo, kpp, delk);
      break;
    case 6:
      if (y == 1) polfrac();
      D2polscat(f, Ninpol, Noutpol, m, del, zzz, Geo, kpp, delk, indexr, disp, y);
      for (i = 0; i <= p; i++)
	{
	  Nin[i] = Ninpol[i];
	  Nout[i] = Noutpol[i];
	}
      Nin[p+3] = Ninpol[p+3];
      Nout[p+3] = Noutpol[p+3];
      break;
    case 7:
      if ( y == 1) polfrac();
      D2polpEexchange(f, b, Nin, Nout, del, zzz, Geo, kpp, delk, indexr, disp, y);
      break;
    case 8:
      if (y == 1) polfrac();
      D2polpscat(f, Ninph, Noutph, m, del, zzz, Geo, kpp, delk, y, indexr, kcz, disp);
      D2polpTAscat(f, NinphTA, NoutphTA, m, del, zzz, Geo, kpp, delk, y, indexr, kcz, disp);
      for (i = 0; i <= p; i++)
	{
	  Nin[i] = Ninph[i] + (2 * NinphTA[i]);
	  Nout[i] = Noutph[i] + (2 * NoutphTA[i]);
	}
      Nin[p+3] = Ninph[p+3] + NinphTA[p+3];
      Nout[p+3] = Noutph[p+3] + NoutphTA[p+3];
      break;
    case 9:
      if (y == 1) polfrac();
      D2polscat(f, Ninpol, Noutpol, m, del, zzz, Geo, kpp, delk, indexr, disp, y);
      D2polpscat(f, Ninph, Noutph, m, del, zzz, Geo, kpp, delk, y, indexr, kcz, disp);
      D2polpTAscat(f, NinphTA, NoutphTA, m, del, zzz, Geo, kpp, delk,y, indexr, kcz, disp);
      for (i = 0; i <= p; i++)
	{
	  Nin[i] = Ninph[i] + Ninpol[i] + (2 * NinphTA[i]);
	  Nout[i] = Noutph[i] + Noutpol[i] + (2 * NoutphTA[i]);
	}
      Nin[p+3] = Ninph[p+3] + Ninpol[p+3] + (2 * NinphTA[p+3]);
      Nout[p+3] = Noutph[p+3] + Noutpol[p+3] + (2 * NoutphTA[p+3]);
      break;
    case 10:
      D3fermiscat(f, Nin, Nout, m, Geo, kpp, delk); /* 3D fermi scattering */
      break;
    case 11:
      if (y == 1) polfrac();
      D2polscat(f, Ninpol, Noutpol, m, del, zzz, Geo, kpp, delk, indexr, disp, y);
      D2polpscat(f, Ninph, Noutph, m, del, zzz, Geo, kpp, delk, y, indexr, kcz, disp);
      D2polFscat(f, NinphF, NoutphF, m, del, zzz, Geo, kpp, delk, y, indexr, kcz, disp);
      for(i = 0; i <= p; i++)
	{
	  Nin[i] = Ninph[i] + Ninpol[i] + NinphF[i];
	  if(Nin[i] != Nin[i] || Nin[i] > 1e50 || Nin[i] < -1e50)
	    {
	      printf("i = %d, Ninph = %e, Ninpol = %e, NinphF = %e\n", i, Ninph[i], Ninpol[i], NinphF[i]);
	      exit(1);
	    }
	  Nout[i] = Noutph[i] + Noutpol[i] + NoutphF[i];
	}
  Nin[p+3] = Ninph[p+3] + Ninpol[p+3] + NinphF[p+3];
 Nout[p+3] = Noutph[p+3] + Noutpol[p+3] + NoutphF[p+3];
 break;
    case 12:
      D2bosescat(f, Ninbos, Noutbos, m, del, zzz, Geo, kpp, delk);
      D2bospscat(f, Ninph, Noutph, m, del, zzz, Geo, kpp, delk, y);
      for (i = 0; i <= p; i++)
	{
	  Nin[i] = Ninph[i] + Ninbos[i];
	  Nout[i] = Noutph[i] + Noutbos[i];
	}
      Nin[p+3] = Ninph[p+3] + Ninbos[p+3];
      Nout[p+3] = Noutph[p+3] + Noutbos[p+3];
      break;
    case 13:
      if (y == 1) polfrac();
      D2polelscat(f, Ninpolel, Noutpolel, m, del, zzz, Geo, kpp, delk, indexr, disp, y);
      D2polpscat(f, Ninph, Noutph, m, del, zzz, Geo, kpp, delk, y, indexr, kcz, disp);
      D2polpTAscat(f, Ninph, Noutph, m, del, zzz, Geo, kpp, delk, y, indexr, kcz, disp);
      for(i = 0; i <= p; i++)
	{
	  Nin[i] = Ninph[i] + Ninpolel[i] + (2 * NinphTA[i]);
	  Nout[i] = Noutph[i] + Noutpolel[i] + (2 * NoutphTA[i]);
	}
      Nin[p+3] = Ninph[p+3] + Ninpolel[p+3] + (2 * NinphTA[p+3]);
      Nout[p+3] = Noutph[p+3] + Noutpolel[p+3] + (2 * NoutphTA[p+3]);
      break;
    case 14:
      if (y == 1) polfrac();
      D2polelscat(f, Ninpolel, Noutpolel, m, del, zzz, Geo, kpp, delk, indexr, disp, y);
      D2polpscat(f, Ninph, Noutph, m, del, zzz, Geo, kpp, delk, y, indexr, kcz, disp);
      D2polpTAscat(f, Ninph, Noutph, m, del, zzz, Geo, kpp, delk, y, indexr, kcz, disp);
      D2polscat(f, Ninpol, Noutpol, m, del, zzz, Geo, kpp, delk, indexr, disp, y);
      for(i = 0; i <= p; i++)
	{
	  Nin[i] = Ninpol[i] + Ninpolel[i] + Ninph[i] + (2 * NinphTA[i]);
	  Nout[i] = Noutpol[i] + Noutpolel[i] + Noutph[i] + (2 * NoutphTA[i]);
	}
      Nin[p+3] = Ninpol[p+3] + Ninpolel[p+3] + Ninph[p+3] + (2 * NinphTA[p+3]);
      Nout[p+3] = Noutpol[p+3] + Noutpolel[p+3] + Noutph[p+3] + (2 * NoutphTA[p+3]);
      break;
    case 15:
      if (y == 1) polfrac();
      D2polelscat(f, Ninpolel, Noutpolel, m, del, zzz, Geo, kpp, delk,indexr, disp, y);
      for(i = 0; i <= p; i++)
	{
	  Nin[i] = Ninpolel[i];
	  Nout[i] = Noutpolel[i];
	}
      Nin[p+3] = Ninpolel[p+3];
      Nout[p+3] = Noutpolel[p+3];
      break;
    }
  return(1);
}
