/*****************************************************************************/
void read_ini_para(int *mpcsteps_h, int *mdsteps_h, int *Lx_h, int *Ly_h,
		   int *radius_h, int *vis_cellsize_h,
		   int *measurement_interval_h, int *start_flow_measurement_h,
		   double *dt_h, double *rho_h, double *alpha_h,
		   double *temperature_h, double *grav_h, double *obsMass_h,
		   double *gridshift_h){

//TEMPORAL VARIABLES READ FROM FILE
  int mpcsteps_t;
  int mdsteps_t;
  int Lx_t;
  int Ly_t;
  int radius_t;
  int vis_cellsize_t;
  int measurement_interval_t;
  int start_flow_measurement_t;

  double dt_t;
  double rho_t;
  double alpha_t;
  double temperature_t;
  double grav_t;
  double obsMass_t;
  double gridshift_t;
  
//CHECKING IF PARAMETERS FILE EXISTS
  FILE *fp = fopen("system_parameters.h","r"); 
    if(fp == NULL){
      printf("ERROR opening file, system_parameters.h doesn't exist!!!\n");
      exit(1);
    }

//TEMPORAL VARIABLES READ FROM FILE
  fscanf(fp,"%lf\n%d\n%d\n%lf\n%d\n%d\n",&dt_t,&mpcsteps_t,&mdsteps_t,&rho_t,&Lx_t,&Ly_t);

  fscanf(fp,"%lf\n%lf\n%lf\n%lf\n%d\n",&alpha_t,&temperature_t,&grav_t,&obsMass_t,&radius_t);

  fscanf(fp,"%d\n%d\n%d\n",&vis_cellsize_t,&measurement_interval_t,&start_flow_measurement_t);

  fscanf(fp,"%lf\n",&gridshift_t);

  fclose(fp);

//VARIABLES READ FROM FILE TO POINTERS
  *mpcsteps_h	= mpcsteps_t;
  *mdsteps_h	= mdsteps_t;
  *Lx_h 	= Lx_t;
  *Ly_h 	= Ly_t;
  *radius_h	= radius_t;
  *vis_cellsize_h = vis_cellsize_t;
  *measurement_interval_h = measurement_interval_t;
  *start_flow_measurement_h = start_flow_measurement_t;

  *dt_h	= dt_t;
  *rho_h	= rho_t;
  *alpha_h	= alpha_t;
  *temperature_h = temperature_t;
  *grav_h 	= grav_t;
  *obsMass_h	= obsMass_t;
  *gridshift_h  = gridshift_t;

  printf("/**********************************************************************/\n");
  printf("Parameters read from file\n");
}

/*****************************************************************************/
//TEST ROUTINE
void print_ini_para(){

  printf("%lf\n%d\n%d\n%lf\n%d\n%d\n",dt,mpcsteps,mdsteps,rho,Lx,Ly);
  printf("%lf\n%lf\n%lf\n%lf\n%d\n",alpha,temperature,grav,obsMass,radius);
  printf("%d\n%d\n%d\n",vis_cellsize,measurement_interval,start_flow_measurement);
  printf("%lf\n",gridshift);
}
