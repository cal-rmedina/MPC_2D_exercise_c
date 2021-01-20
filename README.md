# MPC_2D_exercise_c
2D Multiparticle Collision Dynamics program in c

C code to perform 2D Multiparticle Collision Dynamics Simulations (MPC).

  -Serial C code.  
  -To compile the code, standard C compiler is needed (gcc,icc,visualstudio).
  -Make file "./Makefile" can be used to compile easier, IF YOU DO NOT KNOW 
   WHAT MAKE DOES, compile typing directly on the teminal. 
  -MIT License**

Code doesn't need to be compiled for each execution, once it is compiled,
file "./system_parameters.h" can be modified to change simulation 
parameter (rho, number of particles, steps, etc).

Once the simulation has finished, output files are created in the
directory "./output/*.dat". In the same directory there are different
python to analyze the results. VALUES OF THE SYSTEM IN PYTHON PROGRAMS 
AND SIMULATION MUST MATCH.

// MAIN PROGRAM ***************************************************************

PATH: ./mpc_2d.c

INCLUDES:

  -Libraries.
  -Files included.
  -Main program (calling all the functions executed during the routine).

// MODULES AND ROUTINES *******************************************************

/*
NOTE: This section describes the structure of the code, for the purpose of 
the exercise IT IS NOT REQUIRED TO understand the files and directories
of the program.*/

Functions		"./modules/included_files.h"

  -It includes the path of the files containing the  routines and functions
  executed in main program, if a new routine is created, its path 
  MUST BE DECLARED HERE!!!.

Variables 		"./modules/initialization/variables_declaration.h"

  -It includes the declaration of variables and pointers used in the program.
  -If new variables are used MUST BE DECLARED HERE.

Initial parameters 	"./modules/initialization/read_initial_param.h"

  -Funtion to read the initial parameters

Initialization routines "./modules/initialization/*"
  -Memory allocation, placing object
  -Initial positions
  -Initial velocities

MPC routines 		"./modules/mpc_routines/mpc_routines.h"

MD routines 		"./modules/mpc_routines/md_routines.h"

// METHOD *********************************************************************

Multiparticle Collision Dynamics in CUDA

Multiparticle collision dynamics (MPC) is a simulation method
developed to investigate the properties of mesoscale objects.
MPC is a solid alternative for the simulation of soft matter 
systems. Hydrodynamic interactions, transport of heat & mass,
as well thermal fluctuation are by construction included.

The general strategy to model complex structures with a MPC
solvent is to define an hybrid algorithm where the solvent is
simulated with the MPC technique & both, the solute description
& the solute-solvent interactions, are accounted with specific models.
