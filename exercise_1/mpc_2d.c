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
  read_ini_para(&mpcsteps,&Lx,&Ly,
		&radius,&vis_cellsize,
		&measurement_interval,&start_flow_measurement,
		&h,&rho,&alpha,
		&temperature,&grav,&obsMass,
		&gridshift);

  initialize();

  initialPositions(N,Nobs,rx,ry);
  initialVelocities();

  printf("\nStarting simulation of %i mpc-steps.\n", mpcsteps);

//MPC-loop
  for(int step = 0; step < mpcsteps; step++){

//  Assign random velocities each timestep
    if(step % 1==0)  ramdom_vel_obst(N,Nobs,temperature,obsMass,vx,vy);

//  MPC-routines
    stream(N,grav,h,rx,ry,vx,vy);	// streaming step of the fluid particles
    cells(gridshift);			// sort particles into mpc-cells
    collide();                  	// collision step of the fluid and obstacle parameters

    if (step % 1000 == 0) printf("Step: %u\n", step);

//  Store flowfield
    if ((step >= start_flow_measurement) && (step % measurement_interval == 0)){
      cells(0.0);
      store_flowfield();
    }
  }

  printf("Simulation finished!\n");
  print_flowfield();
  cleanup();

  return 0;
}
