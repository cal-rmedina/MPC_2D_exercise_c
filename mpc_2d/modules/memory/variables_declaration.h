//VARIABLES

int     N;			//Number of fluid particles
int     Nobs;			//Number of obstacle particles
int     Ncell;			//Number of MPC cells
int	Lx,Ly;			//Box size x & y
int	fLx,fLy;		//Box size x & y
int     mdsteps;
int     mpcsteps;
int	radius;
int	vis_cellsize;
int	flowstorecount;
int	start_flow_measurement;
int	measurement_interval;

double	dt;
double	rho;
double	alpha;
double	grav;
double	dLx,dLy;
double	obsMass;
double	spring_force;
double	temperature;

//POINTERS

int 	*list,*head;

double	*rx,*ry;
double	*vx,*vy;
double	*flowfieldx,*flowfieldy;
double	*obsStartx,*obsStarty;

unsigned int *flowcellmass;
