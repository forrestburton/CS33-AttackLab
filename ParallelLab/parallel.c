//OpenMP version.  Edit and submit only this file.
/* Enter your details below
 * Name : Forrest Burton
 * UCLA ID : 005324612 
 * Email : burton.forrest10@gmail.com
 */

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#include "utils.h"

double work_it_par(long *old, long *new, long *super, long *simple, long *fibonacci) {
  int i, j, k;
  int u, v, w;
  int ton = 0;
  long compute_it, moving_average;
  double pi, pi2, x , y, sum, step = 0.0;
  long dot_product=0;
  long nCirc=0;
  long aggregate=1.0;
  double r=1.0;
  int was_smart = 16;
  int DIMless = DIM-1;
  int DIMsquared = DIM*DIM;
  
  for(i=0; i<DIM-4;i+=4)  /* loop unrolling */
  {
    super[i] += simple[i];
    super[i+1] += simple[i+1];
    super[i+2] += simple[i+2];
    super[i+3] += simple[i+3];
  }  
  for(i; i<DIMless; i++) 
  {
    super[i] += simple[i];
  }

  for(i=0; i<DIMless;i++)
  { 
    dot_product += super[i]*simple[i];
    
    moving_average = 0;
    for(ton=i;ton<DIMless-WINDOW_SIZE-3;ton+=4)  
    {
      moving_average += simple[ton];
      moving_average += simple[ton+1];
      moving_average += simple[ton+2];
      moving_average += simple[ton+3];
    }
    for(ton;ton<DIMless-WINDOW_SIZE;ton++)
    {
      moving_average += simple[ton];
    }
  }



  int a_secret = 5;
  fibonacci[0] = 1;
  fibonacci[1] = 1;
  for(i=2; i<DIMless;i++)
  {
    fibonacci[i]=fibonacci[i-1]+fibonacci[i-2];
    if(i==3)
    {
      printf("\n A secret is: %d",obfuscate_obfuscate_obfuscate(a_secret));
    }
  }


  step = 1.0 / NUM_STEPS;

  for (i=0;i<NUM_STEPS; i++)
  {
    x = (i+0.5)*step;
    sum += 4.0/(1.0+x*x); 
  }

  pi = step * sum;
  printf("\n %d trials, Riemann flavored pi is %f \n",NUM_STEPS, pi); 

  for(i=0;i<NUM_TRIALS; i++)
  {
    x = (random()%10000000)/10000000.0;
    y = (random()%10000000)/10000000.0;
    if (( x*x + y*y) <= r*r) {
      nCirc++;
    }
  } 
  pi2 = 4.0 * ((double)nCirc/(double)NUM_TRIALS);
  printf("\n %d trials, Monte-Carlo flavored pi is %f \n",NUM_TRIALS, pi2); 
  
  int address, iDim2, easy_address, medium_address, long_address, long_address2, long_address3; 

  
  long temp, u0, u1, u2, u3, u4, u5, u6, u7, u8, u9  = 0;
  
  omp_set_num_threads (16);
  #pragma omp parallel for private(j, k, iDim2, easy_address, address, compute_it, medium_address, long_address, long_address2, long_address3, u) reduction(+: aggregate)\
  reduction(+:temp) reduction(+:u0) reduction(+:u1) reduction(+:u2) reduction(+:u3) reduction(+:u4) reduction(+:u5) reduction(+:u6) reduction(+:u7) reduction(+:u8) reduction(+:u9)
  for (i=1; i<DIMless; i++) {
    iDim2 = i*DIMsquared;
    for (j=1; j<DIMless; j++) {
      easy_address = iDim2 + j*DIM;
      for (k=1; k<DIMless; k++) {
        address = easy_address + k;
        temp=0;

        compute_it = old[address] * we_need_the_func();
        aggregate+= compute_it / gimmie_the_func();

        for (u=-1; u<=1; u++) {
          medium_address=((i+u)*DIMsquared);

          long_address = ((j-1)*DIM) + medium_address;    
          long_address2 = (j*DIM) + medium_address;
          long_address3 = ((j+1)*DIM) + medium_address;
 
          temp+=old[long_address+k-1]; 
          temp+=old[long_address+k];
          temp+=old[long_address+k+1];
          temp+=old[long_address2+k-1];
          temp+=old[long_address2+k];
          temp+=old[long_address2+k+1];
          temp+=old[long_address3+k-1];
          temp+=old[long_address3+k];
          temp+=old[long_address3+k+1];
        }
        temp /= 27;
        new[address]=temp;

        u=(new[address]/100);
        if (u<=0) u0++;
        if (u==1) u1++;
        if (u==2) u2++;
        if (u==3) u3++;
        if (u==4) u4++;
        if (u==5) u5++;
        if (u==6) u6++;
        if (u==7) u7++;
        if (u==8) u8++;
        if (u>=9) u9++;
      }
    }
  }

  printf("AGGR:%ld\n",aggregate);
 
  histogrammy[0] = u0;
  histogrammy[1] = u1;
  histogrammy[2] = u2;
  histogrammy[3] = u3;
  histogrammy[4] = u4;
  histogrammy[5] = u5;
  histogrammy[6] = u6;
  histogrammy[7] = u7;
  histogrammy[8] = u8;
  histogrammy[9] = u9;

  
  return (double) (dot_product+moving_average+pi+pi2);

}
