/*****************************************************************************/
#define VIRTUAL_PARTICLES

/*****************************************************************************/
// The streaming step: linearly propagate the fluid particles.
void stream()
{
  int i;
  double time_in_wall;

  for(i=0;i<N;i++)
    {
      vx[i]+=grav*dt;
      rx[i]+=vx[i]*dt;
      ry[i]+=vy[i]*dt;

      if (ry[i] < 0) {
	time_in_wall = ry[i] / vy[i];
	rx[i] -= 2.0 * vx[i] * time_in_wall;
	ry[i] -= 2.0 * vy[i] * time_in_wall;
	vx[i] *= -1.0;
	vy[i] *= -1.0;
      }
      if (ry[i] > dLy) {
	time_in_wall = (ry[i]-dLy) / vy[i];
	rx[i] -= 2.0 * vx[i] * time_in_wall;
	ry[i] -= 2.0 * vy[i] * time_in_wall;
	vx[i] *= -1.0;
	vy[i] *= -1.0;
      }
    }
}
/*****************************************************************************/
// sort all (fluid + obstacle) particles into mpc-cells
// gridshift = true if we want a random grid shift
void cells(float gridshift)
{
  int i,icell,ix,iy;
  double xtemp,ytemp,rndx,rndy;

  rndx=RND1*gridshift;
  rndy=RND1*gridshift;

  /* note:
     head[] and list[] together form Ncell "linked lists".
     A linked list is basically a list of things (e.g. particles);
     it is efficient to add items to such a list or remove items from it.
     For this lesson it is not necessary to understand how they work.
     Here we have one list for each mpc-cell, and write into
     each list all the particles that are in the specific mpc-cell.
  */
  for(icell=0;icell<Ncell;icell++)head[icell]=-1;  /* clear the lists */
  for(i=0;i<N+Nobs;i++)
    {
      xtemp=rx[i]+rndx;                 /* grid shift */
      xtemp-=floor(xtemp/dLx)*dLx;      /* periodic boundary condition */
      ix=(int) xtemp;
      ytemp=ry[i]+rndy;
      iy=(int) ytemp;
      icell=ix+iy*Lx;                   /* determine icell, the cell number for this particle */
      list[i]=head[icell];              /* add this particle to the list of cell "icell" */
      head[icell]=i;
    }
}


/*****************************************************************************/
void collide(const float gridshift)
{
  int i,icell;
  double vavx,vavy,vxtemp,vytemp,vrelx,vrely,ca,sa,sin_alpha,cell_mass;
  double pi = acos(-1.0);

  cells(gridshift);		// sort particles into mpc cells, with random grid shift
  ca=cos(alpha);
  sin_alpha=sin(alpha);

  for(icell=0;icell<Ncell;icell++)   /* loop over all the mpc-cells */
    {
      vavx=0.0;
      vavy=0.0;
      vxtemp=0.0;
      vytemp=0.0;
      cell_mass = 0.0;
      // 1. Calculate center of mass velocity of particles in this cell:
      i=head[icell];                /* loop over all      */
      while(i!=-1)                  /*  particles in cell */
	{
	  if (i < N) {              /* i is a fluid particle */
	    vavx += vx[i];
	    vavy += vy[i];
	    cell_mass += 1.0;
	  } else {                  /* i is a obstacle particle */
	    vavx += vx[i] * obsMass;
	    vavy += vy[i] * obsMass;
	    cell_mass += obsMass;
	  }
	  i=list[i];
	}
      // if you know what virtual particles are, good, if not, ask us or just ignore
      //  the lines in between the #ifdef and #endif statements.
#ifdef VIRTUAL_PARTICLES
      if ( ((icell < Lx) || (icell >= Ly*Lx)) && cell_mass < rho) { 

	vavx += sqrt((rho-cell_mass)*temperature) * sqrt(-logl(RND1))*cos(2.0*pi*RND1);
	vavy += sqrt((rho-cell_mass)*temperature) * sqrt(-logl(RND1))*cos(2.0*pi*RND1);
	cell_mass = rho;
      }
#endif
      // 2. If there is more than one particle in the cell, rotate relative velocities:
      if(cell_mass > 1.0)
	{
	  if(RND1<0.5)sa=sin_alpha;  /* randomly choose rotation axis +/- z */
	  else sa=-sin_alpha;
	  vavx/=cell_mass;
	  vavy/=cell_mass;
          i=head[icell];            /* loop over all      */
          while(i!=-1)              /*  particles in cell */
            {
	      vrelx=vx[i]-vavx;
	      vrely=vy[i]-vavy;
	      vxtemp=ca*vrelx-sa*vrely;    /* 2D rotation of      */
	      vytemp=sa*vrelx+ca*vrely;    /* relative velocities */
              vx[i]=vavx+vxtemp;
              vy[i]=vavy+vytemp;
              i=list[i];
            }
	}
    }
}
/*****************************************************************************/
/* Thermostate:
   Because of the driving force (gravity), we put more and more energy into the system over time.
   Thus, the fluid heats up. To counter this, we strongly couple a heat reservoir to the fluid,
   using a simple global thermostate. It scales down all fluid particles velocities to reach the
   desired temperature.
*/

void thermostate() {
  double av_vel2 = 0.0;
  int i;
  for (i = 0; i < N; i++) {
    av_vel2 += vx[i]*vx[i]+vy[i]*vy[i];
  }
  av_vel2 /= (double)N;
  // average velocity we want: v² = kT/m
  const double scale_factor = sqrt(temperature/(av_vel2/2.0));    /* av_vel2/2.0 is the current temperature */
  for (i = 0; i < N; i++) {
    vx[i] *= scale_factor;
    vy[i] *= scale_factor;
  }
  // print out corrected temperature deviation:
  // printf("temperature correction factor: %f\n", scale_factor);
}

/*****************************************************************************/
/* 
   molecular-dynamics part: the motion of the heavy, localized obstacle particles in their harmonic potential;
   of course, having only harmonic potentials, one could calculate the new positions
   analytically, but let's do some simple molecular dynamics (md)... 
*/
void md() {
  int step, i;
  const double mdtime = dt/((double)mdsteps);
  double fx,fy;
  for (step = 0; step < mdsteps; step++) {
    for (i = N; i < N+Nobs; i++) {
      // obsStartx[i-N] is the zero-position of this potential for obstacle particle i, 
      //  rx[i] is the actual position
      fx = spring_force*(obsStartx[i-N]-rx[i]);
      fy = spring_force*(obsStarty[i-N]-ry[i]);
      vx[i] += fx * mdtime;
      vy[i] += fy * mdtime;
      rx[i] += vx[i] * mdtime;
      ry[i] += vy[i] * mdtime;
    }
  }
}
