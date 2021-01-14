0.1
20000
20
5.0
52
52
1.570796326
1.0
0.005
10.0
5
1
100
10000
//TODO: random_shift  1 0 

/******************************************************************************
 dt		mpc time-step
 mpcsteps	mpc-steps
 mdsteps	md-steps per mpc-step
 rho		density (average number of fluid particles per unit cell)
 Lx		system size in unity cells
 Ly		system size in unity cells
 alpha		SRD rotation angle of the relative velocities (radians)
 temperature	temperature of fluid
 grav		gravity; control the force pushing fluid particles in x-direction
 obsMass	mass of the fixed obstacle particles in units of fluid-molecule mass 1
 radius	radius of the obstacle in the flow; choose radius=0 for no obstacle at all
 vis_cellsize	vis_cellsize² mpc-cells are combined into one cell for the flow-field
 		visualization that means: higher value => less vector-arrows, less noise

 measurement_interval		perform flowfield measurements every measurement_interval mpc-steps

 start_flow_measurement	start measuring the flowfield after this many mpc-steps 
 				(to give the system some time for "equilibration")

-------------------------------------------------------------------------------
NOTE: Parameter file for the simulation

  Some values here are integer (e.g. "5") and some are
  floating point values (e.g. "0.05" or "10.0").
  Please don't make ints to floats or vice versa.

******************************************************************************/
