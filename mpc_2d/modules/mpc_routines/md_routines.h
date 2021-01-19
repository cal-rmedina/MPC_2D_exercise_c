/*****************************************************************************/
// Thermostate:
//   Because of the driving force (gravity), we put more and more energy 
//   into the system over time. 
//   Thus, the fluid heats up. To counter this, we strongly couple a heat reservoir to the fluid,
//   using a simple global thermostate. It scales down all fluid particles velocities to reach the
//   desired temperature.

void thermostate(){

  double av_vel2 = 0.0;
  int i;

  for (i = 0; i < N; i++)  av_vel2 += vx[i]*vx[i]+vy[i]*vy[i];

  printf("temperature correction factor: %f\n", av_vel2);

//average velocity we want: v² = kT/m
  av_vel2 /= (double)N;

  printf("temperature correction factor: %f\n", av_vel2);

//av_vel2/2.0 is the current temperature
  const double scale_factor = sqrt(temperature/(av_vel2/2.0));
  for (i = 0; i < N; i++) {
    vx[i] *= scale_factor;
    vy[i] *= scale_factor;
  }
// print out corrected temperature deviation:
 printf("temperature correction factor: %f\n", scale_factor);
}
/*****************************************************************************/
 
// Molecular-Dynamics: 
//   Motion of the heavy, localized obstacle particles (harmonic potential).

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
