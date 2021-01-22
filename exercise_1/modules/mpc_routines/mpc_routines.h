/*****************************************************************************/
// Streaming step: linearly propagate the fluid particles.

void stream(const int np_mpc,
	    const double grav_h,const double h_mpc,
	    double *rx_h,double *ry_h,double *vx_h,double *vy_h){

  double time_in_wall;

//loop over MPC-particles (static obstacle)
  for(int i=0;i<np_mpc;i++){

    vx_h[i]+=grav_h*h_mpc;			// gravity x-direction

    rx_h[i]+=vx_h[i]*h_mpc;
    ry_h[i]+=vy_h[i]*h_mpc;

    rx_h[i] -= floor(rx_h[i]/dLx)*dLx;		// PBC x-direction

//  Wall (No slip boundary condition)

    if (ry_h[i] < 0){
      time_in_wall = ry_h[i] / vy_h[i];
      rx_h[i] -= 2.0 * vx_h[i] * time_in_wall;
      ry_h[i] -= 2.0 * vy_h[i] * time_in_wall;
      vx_h[i] *= -1.0;
      vy_h[i] *= -1.0;
    }
    if (ry_h[i] > dLy){
      time_in_wall = (ry_h[i]-dLy) / vy_h[i];
      rx_h[i] -= 2.0 * vx_h[i] * time_in_wall;
      ry_h[i] -= 2.0 * vy_h[i] * time_in_wall;
      vx_h[i] *= -1.0;
      vy_h[i] *= -1.0;
    }
  }
}
/*****************************************************************************/

// sort all (fluid + obstacle) particles into mpc-cells

void cells(const double gridshift){

  int icell,ix,iy;
  double xtemp,ytemp,rndx,rndy;

// gridshift: 1.0 if we want a random grid shift 0.0 IF NOT
 
  rndx=RND1*gridshift;
  rndy=RND1*gridshift;

//NOTE:
//   head[] and list[] together form Ncell "linked lists".
//     A linked list is basically a list of things (e.g. particles);
//     it is efficient to add items to such a list or remove items from it.
//     For this lesson it is not necessary to understand how they work.
//     Here we have one list for each mpc-cell, and write into
//     each list all the particles that are in the specific mpc-cell.

  for(icell=0;icell<Ncell;icell++)  head[icell]=-1;  // clear the lists

  for(int i=0;i<N+Nobs;i++){
    xtemp=rx[i]+rndx;                 // grid shift
    xtemp-=floor(xtemp/dLx)*dLx;      // periodic boundary condition
    ix=(int) xtemp;
    ytemp=ry[i]+rndy;
    iy=(int) ytemp;
    icell=ix+iy*Lx;                   // determine icell, the cell number for this particle
    list[i]=head[icell];              // add this particle to the list of cell "icell"
    head[icell]=i;
  }
}
/*****************************************************************************/
void collide(){

  int i,icell;
  double vavx,vavy,vxtemp,vytemp,vrelx,vrely,ca,sa,sin_alpha,cell_mass;
  double pi = acos(-1.0);

  ca=cos(alpha);
  sin_alpha=sin(alpha);

  for(icell=0;icell<Ncell;icell++){   // loop over all the mpc-cells
    vavx=0.0;
    vavy=0.0;
    vxtemp=0.0;
    vytemp=0.0;
    cell_mass = 0.0;

//  1. Calculate center of mass velocity of particles in this cell:

    i=head[icell];                /* loop over all      */
    while(i!=-1){                 /*  particles in cell */
      if (i < N) {             /* i is a fluid particle */
        vavx += vx[i];
	vavy += vy[i];
	cell_mass += 1.0;
      }
      else {                  /* i is a obstacle particle */
	vavx += vx[i] * obsMass;
	vavy += vy[i] * obsMass;
	cell_mass += obsMass;
      }
      i=list[i];
    }

// if you know what virtual particles are, good, if not, ask us or just ignore
// the lines in between the #ifdef and #endif statements.

#ifdef VIRTUAL_PARTICLES
    if ( ((icell < Lx) || (icell >= Ly*Lx)) && cell_mass < rho){ 

      vavx += sqrt((rho-cell_mass)*temperature) * sqrt(-logl(RND1))*cos(2.0*pi*RND1);
      vavy += sqrt((rho-cell_mass)*temperature) * sqrt(-logl(RND1))*cos(2.0*pi*RND1);
      cell_mass = rho;
    }
#endif

//  2. If there is more than one particle in the cell, rotate relative velocities:
    if(cell_mass > 1.0){
      if(RND1<0.5)sa=sin_alpha;  /* randomly choose rotation axis +/- z */
      else sa=-sin_alpha;
      vavx/=cell_mass;
      vavy/=cell_mass;
      i=head[icell];            
    
//  loop over all particles in cell
      while(i!=-1){
        vrelx=vx[i]-vavx;
        vrely=vy[i]-vavy;
        vxtemp=ca*vrelx-sa*vrely;    // 2D rotation of relative velocitiess
        vytemp=sa*vrelx+ca*vrely;
        vx[i]=vavx+vxtemp;
        vy[i]=vavy+vytemp;
        i=list[i];
      }
    }
  }
}
/*****************************************************************************/
//Set new obstacle's velocities

// np_min_h = np_mpc, np_max_h = np_obst ----> update obs-particles velocities

void ramdom_vel_obst(const int np_min_h,const int np_max_h,
                     const double temperature_h,const double obs_mass_h,
                     double *vx_h,double *vy_h){
  double vxr;
  double pi = acos(-1.0);
  double vp_obs = sqrt(temperature_h/obs_mass_h);

  for(int i=np_min_h; i<np_min_h+np_max_h; i++){
    vxr = sqrt(-logl(RND1));
    vx_h[i] = vp_obs*vxr*cos(2.0*pi*RND1);

    vxr = sqrt(-logl(RND1));
    vy_h[i] = vp_obs*vxr*cos(2.0*pi*RND1);
  }
}
