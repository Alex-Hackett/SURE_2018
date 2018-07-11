int MTAPiezo2(double W[200][200][200], double kpp[1000], double defe, double defh, double index, double kcz,
	      double W0TA[1000], double Wden0TA[1000], double Wden00TA[1000])
{
  double Bp, E, E1, Ipe, Iph, Ipe1, Iph1, phi, phibad, q, q2, qz, Wden1TA;
  double Ipe10, Iph10, Ipe0, Iph0;
  long double E2;
  double defepart, defhpart, piezopart1, piezopart2e, piezopart2h, Allpartse, Allpartsh, Allpartseh;
  double phiq, qtot;
  Wden1TA = 16 * pi * pi * vTA * vTA * vTA * hb * hb * Sd / eC;
  for (j = 0; j <= p; j++)
    {
      Xj = Xp[j];
      for (jj = 0; jj <= p; jj++)
	{
	  Xjj = Xp[jj];
	  NinkTA[j][jj] = 0;
	  E = f[j][0]; /*polE(kpp[j], index);*/
	  E1 = f[jj][0]; /*polE(kpp[jj], index);*/
	  E2 = fabs(E - E1); /* E - E1 */
	  qtot = E2 / (hb * vTA);
	  phibad = ((kpp[j] * kpp[j]) + (kpp[jj] * kpp[jj]) - (E2 * E2 / (hb * hb * vTA * vTA))) / (2 * kpp[j] * kpp[jj]);
	  if (fabs(phibad) < 1)
	    {
	      phibad = acos(phibad);
	      for (i = 1; i <= p; i++)
		{
		  phi = (i - 0.5) * phibad / p; /*(i - 0.5) * pi / p;*/
		  phiq = atan( (kpp[jj] * sin(phi)) / (-kpp[j] + (kpp[jj] * cos(phi))) );
		  q2 = (kpp[j] * kpp[j]) + (kpp[jj] * kpp[jj]) - (2 * kpp[j] * kpp[jj] * cos(phi));
  /*(kpp[j] * kpp[j]) + (kpp[jj] * kpp[jj]) - (2 * kpp[j] * kpp[jj] * cos(phi));*/
  q = sqrt(q2);
 Ipe1 = 1 + (Mh1 * q * aB / (2 * (Mh1 + Me))) * (Mh1 * q * aB / (2 * (Mh1 + Me)));
 Iph1 = 1 + (Me * q * aB / (2 * (Mh1 + Me))) * (Me * q * aB / (2 * (Mh1 + Me)));
 Ipe = 1 / (sqrt(Ipe1) * Ipe1);
 Iph = 1 / (sqrt(Iph1) * Iph1);
 if ( (((E2 * E2) / (hb * hb * vTA * vTA)) - q2) > 0)
   {
     qz = sqrt(((E2 * E2) / (hb * hb * vTA * vTA)) - q2);
     Bp = 8 * pi * pi * sin(Lz * qz / 2) / (qz * Lz * ( (4 * pi * pi) - (Lz * Lz * qz * qz) ));
     defepart = defe * defe * qtot;
     defhpart = defh * defh * qtot;
     piezopart1 = piezo14 * piezo14 * q2 * ( (qz * qz) + (q2 / 8) ) /
       (16 * pi * pi * eo * eo * einf * einf * qtot * qtot * qtot * qtot * qtot);
     /*
piezopart1 = 2 * piezo14 * piezo14 * q2 * ( (qz * qz) + (q2 * cos(phiq) * cos(phiq) * sin(phiq) * sin(phiq)) ) /
(16 * pi * pi * eo * eo * einf * einf * qtot * qtot * qtot * qtot * qtot);
     */
     /* piezopart2e = defe * piezo14 * q * (qz * (cos(phiq) + sin(phiq)) + (q * cos(phiq) * sin(phiq))) /
( pi * einf * ( E2 * E2 / (hb * hb * vTA * vTA)) );
piezopart2h = defh * piezo14 * q * (qz * (cos(phiq) + sin(phiq)) + (q * cos(phiq) * sin(phiq))) /
( pi * einf * ( E2 * E2 / (hb * hb * vTA * vTA)) ); */
     piezopart2e = 0;
     piezopart2h = 0;
     Allpartse = (defepart + piezopart1 + piezopart2e) * Ipe * Ipe;
     Allpartsh = (defhpart + piezopart1 + piezopart2h) * Iph * Iph;
     Allpartseh = 2 * (-sqrt(defepart * defhpart) + piezopart1) * Ipe * Iph;
     /**** W’s assume area of cavity is 1 cm^2 *****/
     NinkTA[j][jj] += phibad * Bp * Bp * E2 * Xjj * Xj * (Allpartse + Allpartsh
							  + Allpartseh ) / (p * qz); /*Bp * Bp * hb * v * q * E2 * Xjj * Xj * (defe * Ipe - defh * Iph) * (defe * Ipe - defh * Iph) / qz; */
   }
 else
   NinkTA[j][jj] += 0;
 phi = phibad + (( i - 0.5) * (pi - phibad) / p); /*(i - 0.5) * pi / p;*/
 phiq = atan( (kpp[jj] * sin(phi)) / (-kpp[j] + (kpp[jj] * cos(phi))) );
 q2 = (kpp[j] * kpp[j]) + (kpp[jj] * kpp[jj]) - (2 * kpp[j] * kpp[jj] * cos(phi));
 /*(kpp[j] * kpp[j]) + (kpp[jj] * kpp[jj]) - (2 * kpp[j] * kpp[jj] * cos(phi));*/
 q = sqrt(q2);
 Ipe1 = 1 + (Mh1 * q * aB / (2 * (Mh1 + Me))) * (Mh1 * q * aB / (2 * (Mh1 + Me)));
 Iph1 = 1 + (Me * q * aB / (2 * (Mh1 + Me))) * (Me * q * aB / (2 * (Mh1 + Me)));
 Ipe = 1 / (sqrt(Ipe1) * Ipe1);
 Iph = 1 / (sqrt(Iph1) * Iph1);
 if ( (((E2 * E2) / (hb * hb * vTA * vTA)) - q2) > 0)
   {
     qz = sqrt(((E2 * E2) / (hb * hb * vTA * vTA)) - q2);
     Bp = 8 * pi * pi * sin(Lz * qz / 2) / (qz * Lz * ( (4 * pi * pi) - (Lz * Lz * qz * qz) ));
     defepart = defe * defe * qtot;
     defhpart = defh * defh * qtot;
     piezopart1 = piezo14 * piezo14 * q2 * ( (qz * qz) + (q2 / 8) )/
       (16 * pi * pi * eo * eo * einf * einf * qtot * qtot * qtot * qtot * qtot);
     /*
piezopart1 = 2 * piezo14 * piezo14 * q2 * ( (qz * qz) + (q2 * cos(phiq) * cos(phiq) * sin(phiq) * sin(phiq)) ) /
(16 * pi * pi * eo * eo * einf * einf * qtot * qtot * qtot * qtot * qtot);
     */
     /* piezopart2e = defe * piezo14 * q * (qz * (cos(phiq) + sin(phiq)) + (q * cos(phiq) * sin(phiq))) /
(pi * einf * ( E2 * E2 / (hb * hb * vTA * vTA)) );
piezopart2h = defh * piezo14 * q * (qz * (cos(phiq) + sin(phiq)) + (q * cos(phiq) * sin(phiq))) /
(pi * einf * ( E2 * E2 / ( hb * hb * vTA * vTA)) ); */
     piezopart2e = 0;
     piezopart2h = 0;
     Allpartse = (defepart + piezopart1 + piezopart2e) * Ipe * Ipe;
     Allpartsh = (defhpart + piezopart1 + piezopart2h) * Iph * Iph;
     Allpartseh = 2 * (-sqrt(defepart * defhpart) + piezopart1) * Ipe * Iph;
     /**** W’s assume area of cavity is 1 cm^2 *****/
     NinkTA[j][jj] += (pi - phibad) * Bp * Bp * E2 * Xjj * Xj * (Allpartse + Allpartsh + Allpartseh) / (p * qz);
   }
 else
   NinkTA[j][jj] += 0;
		}
	    }
	  else
	    {
	      for (i = 1; i <= p; i++)
		{
		  phi = ((i - 0.5) * pi / p); /*(i - 0.5) * pi / p;*/
		  phiq = atan( (kpp[jj] * sin(phi)) / ( -kpp[j] + (kpp[jj] * cos(phi))) );
		  q2 = (kpp[j] * kpp[j]) + (kpp[jj] * kpp[jj]) + (2 * kpp[j] * kpp[jj] * cos(phi));
		  /*(kpp[j] * kpp[j]) + (kpp[jj] * kpp[jj]) - (2 * kpp[j] * kpp[jj] * cos(phi));*/
		  q = sqrt(q2);
		  Ipe1 = 1 + (Mh1 * q * aB / (2 * (Mh1 + Me))) * (Mh1 * q * aB / (2 * (Mh1 + Me)));
		  Iph1 = 1 + (Me * q * aB / (2 * (Mh1 + Me))) * (Me * q * aB / (2 * (Mh1 + Me)));
		  Ipe = 1 / (sqrt(Ipe1) * Ipe1);
		  Iph = 1 / (sqrt(Iph1) * Iph1);
		  if ( (((E2 * E2) / (hb * hb * vTA * vTA)) - q2) > 0)
		    {
		      qz = sqrt(((E2 * E2) / (hb * hb * vTA * vTA)) - q2);
		      Bp = 8 * pi * pi * sin(Lz * qz / 2) / (qz * Lz * ( (4 * pi * pi) - (Lz * Lz * qz * qz) ));
		      defepart = defe * defe * qtot;
		      defhpart = defh * defh * qtot;
		      piezopart1 = piezo14 * piezo14 * q2 * ( (qz * qz) + (q2 / 8) )/
			(16 * pi * pi * eo * eo * einf * einf * qtot * qtot * qtot * qtot * qtot);
		      /*
piezopart1 = 2 * piezo14 * piezo14 * q2 * ( (qz * qz) + (q2 * cos(phiq)* cos(phiq) * sin(phiq) * sin(phiq)) ) /
(16 * pi * pi * eo * eo * einf * einf * qtot * qtot * qtot * qtot * qtot);
		      */
		      /* piezopart2e = defe * piezo14 * q * (qz * (cos(phiq) + sin(phiq)) + (q * cos(phiq) * sin(phiq))) /
(pi * einf * (E2 * E2 / ( hb * hb * vTA * vTA)) );
piezopart2h = defh * piezo14 * q * (qz * (cos(phiq) + sin(phiq)) + (q * cos(phiq) * sin(phiq))) /
(pi * einf * (E2 * E2 / ( hb * hb * vTA * vTA)) ); */
		      piezopart2e = 0;
		      piezopart2h = 0;
		      Allpartse = (defepart + piezopart1 + piezopart2e) * Ipe * Ipe;
		      Allpartsh = (defhpart + piezopart1 + piezopart2h) * Iph * Iph;
		      Allpartseh = 2 * (-sqrt(defepart * defhpart) + piezopart1) * Ipe * Iph;
		      /**** W’s assume area of cavity is 1 cm^2 *****/
		      NinkTA[j][jj] += pi * Bp * Bp * E2 * Xjj * Xj * (Allpartse + Allpartsh + Allpartseh) / (p * qz);
		    }
		  else
  NinkTA[j][jj] += 0;
		}
	    }
	  /*if(jj > 0)
	    {*/
	  if(j == 0)
	    {
	      Ipe10 = 1 + (Mh1 * kpp[jj] * aB / (2 * (Mh1 + Me))) * (Mh1 * kpp[jj] * aB / (2 * (Mh1 + Me)));
	      Iph10 = 1 + (Me * kpp[jj] * aB / (2 * (Mh1 + Me))) * (Me * kpp[jj] * aB / (2 * (Mh1 + Me)));
	      Ipe0 = 1 / (sqrt(Ipe10) * Ipe10);
	      Iph0 = 1 / (sqrt(Iph10) * Iph10);
	      if( (((f[0][0] - E1) * (f[0][0] - E1) / (hb * hb * v * v)) - (kpp[jj]*kpp[jj])) > 0 )
		{
		  qz = sqrt( ((f[0][0] - E1) * (f[0][0] - E1) / (hb * hb * vTA * vTA)) - (kpp[jj]*kpp[jj]) );
		  Bp = 8 * pi * pi * sin(Lz * qz / 2) / (qz * Lz * ( (4 * pi * pi) - (Lz * Lz * qz * qz)));
		  Wden0TA[jj] = Wden1TA * qz / (2 * pi); /* when k == 0 */
		  W0TA[jj] = Bp * Bp * (E1 - f[0][0]) * (E1 - f[0][0]) * 0.5 * Xjj *
		    (defe * Ipe0 - defh * Iph0) * (defe * Ipe0 - defh * Iph0);
		  /*Bp * Bp * hb * v * kpp[jj] * (E1 - f[0][0]) * 0.5 * Xjj * (defe * Ipe0 - defh * Iph0) * (defe * Ipe0 - defh * Iph0);*/
		  Wden00TA[jj] = Wden1TA * qz * D2A / (8 * pi * pi); /* when k != 0 */
		}
	      else
		{
		  W0TA[jj] = 0;
		  Wden0TA[jj] = 1;
		  Wden00TA[jj] = 1;
		}
	    }
	  /*}
else
{
W0[0] = 0;
Wden0[0] = 1;
Wden00[0] = 1;
}*/
	  /* printf("W = %e \n", W[j][jj][i]); */
	}
    }
for (j = 0; j <= p; j+= (p-2)/2)
  {
    for (i = 1; i <= p; i++)
      WtTA[j] += WTA[j][i][0] * delk[i] * kpp[i] * 2 * pi;
  }

for (i = 0; i <= p; i++)
  {
    for (ii = 0; ii <= p; ii++)
      NinkTA[i][ii] *= 4 / Wden1TA;
  }

/* ---- Used to check symmetry of the Nink ----*/

 for(i = 0; i <= p; i++)
   {
     for(ii = 0; ii <= p; ii++)
       {
	 if(fabs(NinkTA[i][ii] - NinkTA[ii][i]) / fabs(NinkTA[i][ii]) > 1e-14)
	   printf("Matrix element assymmetry for TA phonons, %e %e %e\n", NinkTA[i][ii], NinkTA[ii][i],fabs(NinkTA[i][ii] - NinkTA[ii][i]) /
		  fabs(NinkTA[i][ii]));
       }
   }
 NinkTA[p][p+1] = 0;
 NinkTA[p+1][p] = 0;
 /* scattering matrix element vs |q|
sprintf(outfile4, "Mvsq.d", xx);
fq = fopen(outfile4, "w");
fprintf(fq, "x ");
for(i = 0 ; i <= p; i++)
fprintf(fq, "%e ", kpp[i]);
fprintf(fq, "\n");
for ( i = 0; i <= p; i++)
{
fprintf(fq, "%e ", kpp[i]);
for(ii = 0; ii <= p; ii++)
fprintf(fq, "%e ", Nink[i][ii]);
fprintf(fq, "\n");
}
fclose(fq);
exit(1);*/
 return(1);
}


