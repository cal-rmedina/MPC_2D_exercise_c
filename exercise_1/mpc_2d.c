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
		&h,&rho,&alpha,
		&temperature,&grav,&obsMass,
		&gridshift);

//TEST FUNCTION: Print system_parameters on screen, comment/uncomment
  print_ini_para();

  initialize();

  initialPositions(N,Nobs,rx,ry);
  initialVelocities();

  printf("\nStarting simulation of %i mpc-steps.\n", mpcsteps);

  for(int step = 0; step < mpcsteps; step++){

//  TEST FUNCTION: Kinetic energy obstacle, /output/kin_ene.py, comment/uncomment
    if((step%100==0) && (Nobs>0))  print_kin_energy(step,N,Nobs,h,obsMass,vx,vy);

/*****************************************************************************/
//TODO: Play with the numbers to activate both, if you don't want a 
//	particular routine comment/uncomment it.

//  Assign random velocities each timestep
    if(step % 1==0)  ramdom_vel_obst(N,Nobs,temperature,obsMass,vx,vy);

    if(step % mpcsteps==0)  thermostate();  //Thermostate
/*****************************************************************************/
 
//  MPC-routines
    stream(N,grav,h,rx,ry,vx,vy);	// streaming step of the fluid particles
    cells(gridshift);			// sort particles into mpc-cells
    collide();                  	// collision step of the fluid and obstacle parameters

    if (step % 1000 == 0) printf("Step: %u\n", step);
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
