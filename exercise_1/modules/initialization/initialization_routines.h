/*****************************************************************************/
// Initialize variables, allocate memory.
// Place the heavy obstacle particles that are fixated by springs also 
// counts the number of obstacle particles Nobs

void initialize(){

  int i;
  Ncell=Lx*(Ly+1);
  dLx=(double)Lx;
  dLy=(double)Ly;
  N=(int) ((dLx)*(dLy)*rho);    //number of fluid particles
  flowstorecount = 0;

//number of cells in the flow field visualization
  fLx = ((int)(ceil(Lx/(double)(vis_cellsize))+0.5));
  fLy = ((int)(ceil(Ly/(double)(vis_cellsize))+0.5));

// Allocate memory CPU
  obsStartx = (double *) malloc(sizeof(double)*(2*radius+1)*(2*radius+1));
  obsStarty = (double *) malloc(sizeof(double)*(2*radius+1)*(2*radius+1));

// Place obstacle
  int x,y;
  Nobs = 0;
  for (y = -(radius-1); y <= (radius-1); y++){
    for (x = -(radius-1); x <= (radius-1); x++){
      if (x*x+y*y <= (radius-1)*(radius-1)){
	obsStartx[Nobs] = (double)x;
	obsStarty[Nobs] = (double)y;
	Nobs++;
      }
    }
  }

// Shift all particles into the middle of the box:
  for (x = 0; x < Nobs; x++){
    obsStartx[x] += ((double)Lx)/2.0;
    obsStarty[x] += ((double)Ly)/2.0;
  }

// Allocate variavles
  head=(int *) malloc(sizeof(int)*Ncell);       
  list=(int *) malloc(sizeof(int)*(N+Nobs));
  rx=(double *) malloc(sizeof(double)*(N+Nobs));
  ry=(double *) malloc(sizeof(double)*(N+Nobs));
  vx=(double *) malloc(sizeof(double)*(N+Nobs));
  vy=(double *) malloc(sizeof(double)*(N+Nobs));

// fLx & fLy: number of vectors in the flow-field visualization (x & y)
  flowfieldx = (double *) malloc(sizeof(double)*fLx*fLy); 
  flowfieldy = (double *) malloc(sizeof(double)*fLx*fLy);
  flowcellmass = (unsigned int *) malloc(sizeof(unsigned int)*fLx*fLy);

  for (i = 0; i < fLx*fLy; i++){
    flowfieldx[i] = 0.0;
    flowfieldy[i] = 0.0;
    flowcellmass[i] = 0;
  }

  printf("/**********************************************************************/\n");
  printf("CPU memory allocation and placed object\n");
}
/*****************************************************************************/
// Initial random positions of the fluid particles

void initialPositions(const int np_mpc,const int np_obs,
                      double *rx_h,double *ry_h){
  int i;
  for(i=0; i<np_mpc; i++){
    rx_h[i]=RND1*dLx;        //RND1 yields uniformly distributed
    ry_h[i]=RND1*dLy;        //random numbers [0,1]
  }

  // initial positions of the obstacle particles
  for(i=np_mpc; i<np_mpc+np_obs; i++){
    rx_h[i] = obsStartx[i-np_mpc];
    ry_h[i] = obsStarty[i-np_mpc];
  }

  printf("/**********************************************************************/\n");
  printf("Initial positions\n");
}
/*****************************************************************************/
// initial velocities of the particles
// average kinetic energy of one particle in 2-dim = kT;   k = 1.0 in our units;
// therefore: <1/2 m v²> = kT;   with m = 1.0 for the fluid particles

void initialVelocities(){
  int i;
  double vxtemp,vytemp;

  double vxr;
  double vp_sol = sqrt(temperature);		//Mass_solvent= 1.0
  double vp_obs = sqrt(temperature/obsMass);

  double pi = acos(-1.0);

  vxtemp=vytemp=0.0;
  for(i=0;i<N;i++){

//  Random numbers, variance 1.0 => temperature 1.0
    vxr = sqrt(-logl(RND1));
    vx[i]=vp_sol*vxr*cos(2.0*pi*RND1);

    vxr = sqrt(-logl(RND1));
    vy[i]=vp_sol*vxr*cos(2.0*pi*RND1);

    vxtemp+=vx[i];
    vytemp+=vy[i];
    }
  vxtemp/=(double) N;
  vytemp/=(double) N;

  for(i=0;i<N;i++){
    vx[i]-=vxtemp;     // Center-of-mass velocity of the
    vy[i]-=vytemp;     // whole system should vanish
    }

//  Initial velocities of the obstacle-particles: 
  for (i = N; i < N+Nobs; i++){
    vxr = sqrt(-logl(RND1));
    vx[i] = vp_obs*vxr*cos(2.0*pi*RND1);

    vxr = sqrt(-logl(RND1));
    vy[i] = vp_obs*vxr*cos(2.0*pi*RND1);
  }

  printf("/**********************************************************************/\n");
  printf("Initial velocities, CM-velocity substracted\n");
}
/*****************************************************************************/
//Free CPU-memory (Pointers)
void cleanup(){
  free(list);
  free(head);
  free(rx);
  free(ry);
  free(vx);
  free(vy);
  free(flowfieldx);
  free(flowfieldy);
  free(flowcellmass);
  free(obsStartx);
  free(obsStarty);
}
