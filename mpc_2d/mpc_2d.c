/******* DEFINITIONS *********************************************************/
//Random number [0,1]
#define RND1   ((double)((double)rand()/(double)RAND_MAX))

//Virtual particles
#define VIRTUAL_PARTICLES

/******* LIBRARIES ***********************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/******* FILES AND MODULES ***************************************************/
#include "./modules/included_files.h"

int main(){

//READING INITIAL PARAMETERS FROM FILE (./system_parameters.h)
  read_ini_para(&mpcsteps,&mdsteps,&Lx,&Ly,
		&radius,&vis_cellsize,
		&measurement_interval,&start_flow_measurement,
		&dt,&rho,&alpha,
		&temperature,&grav,&obsMass,
		&gridshift);

  printf("%lf\n%d\n%d\n%lf\n%d\n%d\n",dt,mpcsteps,mdsteps,rho,Lx,Ly);
  printf("%lf\n%lf\n%lf\n%lf\n%d\n",alpha,temperature,grav,obsMass,radius);
  printf("%d\n%d\n%d\n",vis_cellsize,measurement_interval,start_flow_measurement);
  printf("%lf\n",gridshift);

  int step;
  initialize();
  initialPositions();
  initialVelocities();

  printf("\nStarting simulation of %i mpc-steps.\n", mpcsteps);

  for(step = 0; step < mpcsteps; step++){

   if (step % 100 == 0) thermostate();      // call the thermostate every 100 steps
   md();                                    // calculate movement of the obstacle parameters

//  MPC- routines
   stream();                                // streaming step of the fluid particles
   cells(gridshift);				//sort particles into mpc-cells
   collide();                               // collision step of the fluid and obstacle parameters

   if (step % 1000 == 0) printf("Step: %u\n", step);
   if ((step >= start_flow_measurement) && (step % measurement_interval == 0)) store_flowfield();
  }

  printf("Simulation finished!\n");
 print_flowfield();
 cleanup();

  return 0;
}
