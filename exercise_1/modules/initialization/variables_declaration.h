//VARIABLES

int     N;			//Number of fluid particles
int     Nobs;			//Number of obstacle particles
int     Ncell;			//Number of MPC cells
int	Lx,Ly;			//Box size x & y
int	fLx,fLy;		//Vectors in the flow-field visualization
int     mdsteps;
int     mpcsteps;
int	radius;
int	vis_cellsize;
int	flowstorecount;
int	start_flow_measurement;
int	measurement_interval;

double	h;
double	rho;
double	alpha;
double	grav;
double	dLx,dLy;
double	obsMass;
double	temperature;
double  gridshift;

//POINTERS

int 	*list,*head;

double	*rx,*ry;
double	*vx,*vy;
double	*flowfieldx,*flowfieldy;
double	*obsStartx,*obsStarty;

unsigned int *flowcellmass;
