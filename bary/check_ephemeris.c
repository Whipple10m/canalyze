#include <stdio.h>
#include <math.h>

main()

{

double epoch1,epoch2,t1,t2,nu1,nu2,nudot1,nudot2;
double number;

printf("\nEnter epoch 1 in MJD: ");
scanf("%F",&epoch1);
printf("\nEnter epoch 2 in MJD: ");
scanf("%F",&epoch2);
printf("\nEnter t_jpl 1 in sec: ");
scanf("%F",&t1);
printf("\nEnter t_jpl 2 in sec: ");
scanf("%F",&t2);
printf("\nEnter nu 1 in Hz: ");
scanf("%F",&nu1);
printf("\nEnter nu 2 in Hz: ");
scanf("%F",&nu2);
printf("\nEnter nudot 1 in Hz^2: ");
scanf("%F",&nudot1);
printf("\nEnter nudot 2 in Hz^2: ");
scanf("%F",&nudot2);

epoch1 = epoch1*86400.0e+0 + t1;
epoch2 = epoch2*86400.0e+0 + t2;

/* go forward from 1 to 2 */

number = nu1*(epoch2-epoch1);
number = number + (0.5*nudot1*((epoch2-epoch1)*(epoch2-epoch1)));
printf("\nGoing forward, number = %.3f\n",number);

/* go backwards from 2 to 1 */

number = nu2*(epoch1-epoch2);
number = number + (0.5*nudot2*((epoch1-epoch2)*(epoch1-epoch2)));
printf("\nGoing backwards, number = %.3f\n",number);

}
