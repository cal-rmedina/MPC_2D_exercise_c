/******************************************************************************
* INITIAL CONFIGURATION: SOLVENT RANDOM, VELOCITIES ACCORDING M-B.            *
* AUTHOR: CAL-RMEDINA.                                                        *
* DATE:   14.01.2021.                                                         *
******************************************************************************/

/******* LIBRARIES ***********************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/******* FILES AND MODULES ***************************************************/
#include "./modules/input.h"
#include "./modules/included_files.h"

/******* MAIN PROGRAM ********************************************************/
int main(){

//READING INITIAL PARAMETERS FROM FILE (./system_parameters.h)
  read_ini_para(&mpcsteps,&mdsteps,&Lx,&Ly,
		&radius,&vis_cellsize,
		&measurement_interval,&start_flow_measurement,
		&dt,&rho,&alpha,
		&temperature,&grav,&obsMass);

/*
  int step;
  read_parameters("parameter.ini");
  initialize();
  printf("Starting simulation of %i mpc-steps.\n", mpcsteps);
  for(step = 0; step < mpcsteps; step++)
    {
      if (step % 100 == 0) thermostate();      // call the thermostate every 100 steps
      md();                                    // calculate movement of the obstacle parameters
      stream();                                // streaming step of the fluid particles
      collide();                               // collision step of the fluid and obstacle parameters
      if (step % 1000 == 0) printf("Step: %u\n", step);
      if ((step >= start_flow_measurement) && (step % measurement_interval == 0)) store_flowfield();
    }
  printf("Simulation finished!\n");
  print_flowfield();
  cleanup();

*/

  return 0;
}
