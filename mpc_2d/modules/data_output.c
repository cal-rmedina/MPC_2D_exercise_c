
/*
  Visualization of the flowfield:
  The box is separated into visualization cells (without grid shift) that are 
   of the size vis_cellsize*vis_cellsize times the unit mpc-cell.

  store_flowfield() calculates the total momenta of the particles in the vis. cells and
   adds them up in flowfieldx[] and flowfieldy[]. It should be called after every "so many" mpc-steps.

  print_flowfield() prints out the average momenta of the fluid to a file.
   It should be called once at the end of the simulation.
  
*/

void store_flowfield() {
  cells(0.0);   /* sort particles into mpc-cells without grid shift */
  double cell_sum_vx, cell_sum_vy;
  int i, cellx, celly;
  unsigned int mass;
  for (celly = 0; celly < Ly; celly++) {
    for (cellx = 0; cellx < Lx; cellx++) {   /* loop over all mpc-cells */
      cell_sum_vx = 0.0;
      cell_sum_vy = 0.0;
      i = head[cellx+celly*Lx];
      mass = 0;
      while (i != -1) {                      /* loop over all particles in this mpc-cell,     */
	if (i < N) {                         /*  sum up the velocities of the fluid particles */
	  cell_sum_vx += vx[i];
	  cell_sum_vy += vy[i];
	  ++mass;
	}
	i = list[i];
      }
      // add the total (fluid-) momentum of this mpc-cell to the correct visualization cell
      flowfieldx[cellx/vis_cellsize+(celly/vis_cellsize)*fLx] += cell_sum_vx;  /* fLx = number of vis. cells in x-direction */
      flowfieldy[cellx/vis_cellsize+(celly/vis_cellsize)*fLx] += cell_sum_vy;
      flowcellmass[cellx/vis_cellsize+(celly/vis_cellsize)*fLx] += mass;
    }
  }
  ++flowstorecount;
}

void print_flowfield() {
  FILE* flowfield = fopen("./output/flowfield.dat", "w");
  int x,y;
  for (y = 0; y < fLy; y++) {
    for (x = 0; x < fLx; x++) {
      if (flowcellmass[x+fLx*y] > 0) 
	fprintf(flowfield, "%f %f %f %f\n", (((double)x)+0.5)*vis_cellsize, (((double)y)+0.5)*vis_cellsize, flowfieldx[x+fLx*y]/((double)flowcellmass[x+fLx*y]), flowfieldy[x+fLx*y]/((double)flowcellmass[x+fLx*y]));
      else
	fprintf(flowfield, "%f %f 0.0 0.0", (((double)x)+0.5)*vis_cellsize, (((double)y)+0.5)*vis_cellsize);
    }
  }
  fclose(flowfield);

  double flow;
  unsigned int mass;
  FILE* flowprofile = fopen("./output/flowprofile.dat", "w");
  fprintf(flowprofile, "# column 1: y position; column 2: average flow\n");
  for (y = 0; y < fLy; y++) {
    flow = 0.0;
    mass = 0;
    for (x = 0; x < fLx; x++) {
      flow += flowfieldx[x+y*fLx];
      mass += flowcellmass[x+fLx*y];
    }
    if (mass > 0) fprintf(flowprofile, "%f %f\n", (((double)y)+0.5)*vis_cellsize, flow/((double)mass));
    else fprintf(flowprofile, "%f \n", (((double)y)+0.5)*vis_cellsize);
  }
  fclose(flowprofile);  

  FILE* density = fopen("./output/density.dat", "w");
  const double flowstorecountd = (double)flowstorecount;
  for (y = 0; y < fLy; y++) {
    for (x = 0; x < fLx; x++) {
      fprintf(density, "%f %f %f\n", (((double)x)+0.5)*vis_cellsize, (((double)y)+0.5)*vis_cellsize, ((double) flowcellmass[x+fLx*y])/flowstorecountd/((double)(vis_cellsize*vis_cellsize)));
    }
    fprintf(density, "\n");
  }
  fclose(density);
}


// not used at the moment;
// call this e.g. every 100 mpc steps if you want to create "movies" from your system; watch it in VMD (probably not working on your machines)
void writeout_xyz(FILE* f, int step) {
  fprintf(f, "%u\nfluids\n", N / 50 + Nobs);
  int particle;
  for (particle = 0; particle < N; particle += 50)
    fprintf(f, "F %f %f 0 %u\n", rx[particle]-floor(rx[particle]/dLx)*dLx, ry[particle]-floor(ry[particle]/dLy)*dLy, step);
  for (particle = N; particle < N+Nobs; particle += 1)
    fprintf(f, "M %f %f 0 %u\n", rx[particle]-floor(rx[particle]/dLx)*dLx, ry[particle]-floor(ry[particle]/dLy)*dLy, step);
}
