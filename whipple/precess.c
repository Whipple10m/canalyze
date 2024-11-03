#include <stdio.h>

main()
{
  double dec, decsec, ra, rasec;
  int rahour, ramin, decdeg, decmin;
  double starting_epoch, ending_epoch;
  int i[4], j[4], k;
  char cr, cd;

  printf("Entering starting epoch: ");
  scanf("%lf",&starting_epoch);
  printf("Entering ending epoch: ");
  scanf("%lf",&ending_epoch);
  printf("Enter RA (hh mm ss.s): ");
  scanf("%d %d %lf",&rahour, &ramin, &rasec);
  printf("Enter DEC (dd mm ss.s): ");
  scanf("%d %d %lf",&decdeg, &decmin, &decsec);

  sla_dtf2r_(&rahour,&ramin,&rasec,&ra,&k);
  sla_daf2r_(&decdeg,&decmin,&decsec,&dec,&k);

  printf("RA: %lf DEC: %lf\n",ra,dec);
  sla_preces_("FK4",&starting_epoch,&ending_epoch,&ra,&dec,3);
  printf("RA: %lf DEC: %lf\n",ra,dec);

  k = 1;
  sla_dr2tf_(&k,&ra,&cr,i,1);
  sla_dr2af_(&k,&dec,&cd,j,1);

  printf("New RA, DEC = %c%d %d %d.%d, %c%d %d %d.%d\n",cr,i[0],i[1],i[2],
	 i[3],cd,j[0],j[1],j[2],j[3]);
}
