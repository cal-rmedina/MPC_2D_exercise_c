0.1
40000
20
5.0
52
52
1.570796326
1.0
0.0005
10.0
5
1
10
10000
1.0
1000.0
/******************************************************************************
REFERENCE:

0.1		h
20000		mpcsteps
20		mdsteps
5.0		rho
52		Lx
52		Ly
1.570796326	alpha
1.0		temperature
0.005		grav
10.0		obsMass
5		radius
1		vis_cellsize
100		measurement_interval
10000		start_flow_measurement
1.0		gridshift

NOTE:	Some values here are integer (e.g. "5") and some are
  	floating point values (e.g. "0.05" or "10.0").
	Please don't make ints to floats or vice versa.

Take the list above as a reference in case you don't remember the original
values or the program doesn't run with your new parameters set.

DESCRIPTION OF VARIABLES:

 h		mpc time-step
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

 gridshift			1.0 or 0.0 to activate/deactivate random shift in collision routine

******************************************************************************/
