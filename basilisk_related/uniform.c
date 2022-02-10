/** Header file, watch out for the order. */
//#include "/home/moin/hanul/basilisk/basilisk/src/grid/multigrid3D.h"
#include "grid/multigrid3D.h"
#include "navier-stokes/centered.h"
//#include "navier-stokes/centered_connor.h" //frequency damping
#include "two-phase.h"
#include "tension.h"
//#include "tag.h"
#include "navier-stokes/conserving.h"
//#include "output_fields/output_vtu_foreach.h"
#include "output_fields/output_vtr_foreach.h"
#include "dev_input.h"
/** define number of processors in each direction. */
//case1: (12,2,2)
#define nx_proc 24
#define ny_proc 4
#define nz_proc 4
// epsilon: 0.004 -> n5low: 0.08%; n5high: 0.26%; n6: 0.31%
#define epsilon 0.001
#define Nr      200
#define M_PI    3.141592
/**
We define the radius of the jet, the initial jet length, the Reynolds
number and the surface tension coefficient. */
#define radius 0.0039
#define length 0.001
//n5low: 2.2; n5high: 1.7; n6:2.5
#define omega  2.5/0.0039
#define tzero  .09
// n5low: 0.201671; n5high: 0.655438; n6: 0.778920
#define Famp   radius * 4. * 0.778920
#define Nazm   6.
/**
Output for vtk files */
int    Nx      = 32*nx_proc; //32 for YS
double Nlen    = 0.24; // about 30 points within the radius
/** boundary conditions */
scalar f0[], ub[], ur[], ui[], vr[], vi[], wr[], wi[];
// Left face B.C, x-
#if dimension > 2
u.n[left]  = dirichlet(ub[] + epsilon*(ur[]*cos(omega*(t-tzero)) + ui[]*sin(omega*(t-tzero)))*(t >= tzero));
// u.y
u.t[left]  = dirichlet(0. +   epsilon*(vr[]*cos(omega*(t-tzero)) + vi[]*sin(omega*(t-tzero)))*(t >= tzero));
// u.z
u.r[left]  = dirichlet(0. +   epsilon*(wr[]*cos(omega*(t-tzero)) + wi[]*sin(omega*(t-tzero)))*(t >= tzero));
p[left]    = neumann(0);
f[left]    = f0[];
// Right face B.C, x+
u.n[right] = u.n[] > 0 ? neumann(0) : 0;
p[right]   = dirichlet(0);
pf[right]   = dirichlet(0);
// Front face B.C, z+
u.n[front] = u.n[] > 0 ? neumann(0) : 0;
p[front]   = dirichlet(0);
pf[front]   = dirichlet(0);
// Back face B.C, z-
u.n[back] = u.n[] < 0 ? neumann(0) : 0;
p[back]   = dirichlet(0);
pf[back]   = dirichlet(0);
// Top face B.C, y+
u.n[top] = u.n[] > 0 ? neumann(0) : 0;
p[top]   = dirichlet(0);
pf[top]   = dirichlet(0);
// Bottom face B.C, y-
u.n[bottom] = u.n[] < 0 ? neumann(0) : 0;
p[bottom]   = dirichlet(0);
pf[bottom]   = dirichlet(0);
#endif

/**
The program can take two optional command-line arguments: the maximum
level and the error threshold on velocity. */

int main (int argc, char * argv[])
{
  size(Nlen); // Domain length
  origin(0, -Nlen*ny_proc/nx_proc/2, -Nlen*nz_proc/nx_proc/2);
  dimensions(nx = nx_proc, ny = ny_proc, nz = nz_proc);
  // Domain length follows the number of processors
  N = Nx; // Number of stencils
  //n5low
//  rho1 = 998., rho2 = 1.0;
//  mu1 = 0.43e-3; mu2 = 0.47e-5;  
//  f.sigma = 1.3e-3;
  //n5high
//  rho1 = 890., rho2 = 1.0;
//  mu1 = 0.43e-3; mu2 = 0.47e-5;  
//  f.sigma = 0.34e-3;
  //n6
rho1 = 998., rho2 = 1.2;
mu1 = 0.5e-3; mu2 = 1.8e-5;  
f.sigma = 0.8e-3;

  run();
}

/**
## Initial conditions */
event init (t = 0) {
  if (!restore (file = "dump")) {
    fraction (f0, sq(radius) - sq(y) - sq(z)); 
    double iniB[Nr];
    double iniUr[Nr];
    double iniUi[Nr];
    double iniVr[Nr];
    double iniVi[Nr];
    double iniWr[Nr];
    double iniWi[Nr];
    
    // Load initial baseflow
    //FILE * fpb = fopen("../newinput/n6_base.txt","r");
    FILE * fpb = fopen("./n6_base.txt","r");
    if (fpb == NULL) {
      fprintf(stderr, "Could not find the inputfile!\n");
      return 1;
    }
    readDataNew(fpb, Nr, iniB);
    fclose(fpb);
   
    // Load initial <u> velocity
    FILE * fpu = fopen("./n6_ur.txt","r");
    if (fpu == NULL) {
      fprintf(stderr, "Could not find the inputfile!\n");
      return 1;
    }
    readDataNew(fpu, Nr, iniUr);
    fclose(fpu);
    fpu = fopen("./n6_ui.txt","r");
    if (fpu == NULL) {
      fprintf(stderr, "Could not find the inputfile!\n");
      return 1;
    }
    readDataNew(fpu, Nr, iniUi);
    fclose(fpu);
    
    // Load initial <v> velocity
    FILE * fpv = fopen("./n6_vr.txt","r");
    if (fpv == NULL) {
      fprintf(stderr, "Could not find the inputfile!\n");
      return 1;
    }
    readDataNew(fpv, Nr, iniVr);
    fclose(fpv);
    fpv = fopen("./n6_vi.txt","r");
    if (fpv == NULL) {
      fprintf(stderr, "Could not find the inputfile!\n");
      return 1;
    }
    readDataNew(fpv, Nr, iniVi);
    fclose(fpv);
    
    // Load initial <w> velocity
    FILE * fpw = fopen("./n6_wr.txt","r");
    if (fpw == NULL) {
      fprintf(stderr, "Could not find the inputfile!\n");
      return 1;
    }
    readDataNew(fpw, Nr, iniWr);
    fclose(fpw); 
    fpw = fopen("./n6_wi.txt","r");
    if (fpw == NULL) {
      fprintf(stderr, "Could not find the inputfile!\n");
      return 1;
    }
    readDataNew(fpw, Nr, iniWi);
    fclose(fpw); 
    // Allocate the initial velocity using V[].
    foreach() {
      double r = sqrt(sq(y) + sq(z));
      double the;
      if ( y > 0. ) {
        the = acos(z/r);
      } else {
        the = 2.*M_PI-acos(z/r);
      }
      // Something like, getData(r, the, iniU);
      ub[] = getDataNew(r, 3*radius, Nr, iniB); 
      ur[] = 4*cos(Nazm*the)*getDataNew(r, 3*radius, Nr, iniUr); 
      ui[] = 4*cos(Nazm*the)*getDataNew(r, 3*radius, Nr, iniUi); 
      vr[] = 4*(getDataNew(r, 3*radius, Nr, iniVr)*cos(Nazm*the)*sin(the) - getDataNew(r, 3*radius, Nr, iniWi)*sin(Nazm*the)*cos(the));
      vi[] = 4*(getDataNew(r, 3*radius, Nr, iniVi)*cos(Nazm*the)*sin(the) + getDataNew(r, 3*radius, Nr, iniWr)*sin(Nazm*the)*cos(the));
      wr[] = 4*(getDataNew(r, 3*radius, Nr, iniVr)*cos(Nazm*the)*cos(the) + getDataNew(r, 3*radius, Nr, iniWi)*sin(Nazm*the)*sin(the)); 
      wi[] = 4*(getDataNew(r, 3*radius, Nr, iniVi)*cos(Nazm*the)*cos(the) - getDataNew(r, 3*radius, Nr, iniWr)*sin(Nazm*the)*sin(the)); 
      // Check r, the, Nr, Nthe etc..      
      f[] = f0[]*(x < length); // initial liquid jet
      u.x[] = ub[]*(x < length);
    }
    fprintf(stderr, "Initialization completed! Start run()!\n");
    boundary ({f, u.x});
  }
}

/** test function for oscillating interface */
event interface (i++) {
  if ( t >= tzero ) {
    fraction (f0, sq(radius + Famp*epsilon*sin(omega*(t-tzero))*cos(Nazm*( y > 0. ? acos(z/sqrt(sq(y)+sq(z))) : 2.*M_PI-acos(z/sqrt(sq(y)+sq(z)))))) - sq(y) - sq(z)); 
  }
}

/** Outputs We log some statistics on the solver. */
event logfile (i++) {
  if (i == 0)
    fprintf (ferr,
	     "t dt mgp.i mgpf.i mgu.i grid->tn perf.t perf.speed\n");
    fprintf (ferr, "%g %g %d %d %d %ld %g %g\n", 
	   t, dt, mgp.i, mgpf.i, mgu.i, grid->tn, perf.t, perf.speed);
}

event end ( t = 10. ){
}

void backup_fields_vtr (scalar * list, vector * vlist, int i)
{
  char name[80], subname[80];
  FILE * fp;
   sprintf(name, "optimalInt_%06d_n%3.3d.vtr", i, pid());
   fp = fopen(name, "w");
   bool linear=false;
   output_vtr_ascii_foreach (list, vlist, Nlen, Nx, nx_proc, ny_proc, nz_proc, fp, linear);
   fclose (fp);

 @if _MPI
   if (pid() == 0){
     sprintf( name, "optimalInt_%06d.pvtr", i);
     sprintf( subname, "optimalInt_%06d", i);
     fp = fopen( name, "w");
     output_pvtr_ascii (list, vlist, Nlen, Nx, nx_proc, ny_proc, nz_proc, fp, subname);
     fclose(fp);
     }
   MPI_Barrier(MPI_COMM_WORLD);
 @endif
 } // End of the backup_fields function

trace
void dump3 (struct Dump p)
{
  if (npe() < 2)
    fprintf (ferr, "dump3() is aimed at MPI runs only...\n");
  FILE * fp = p.fp;
  char def[] = "dump", * file = p.file ? p.file : p.fp ? NULL : def;
  if (fp != NULL || file == NULL) {
    fprintf (ferr, "dump(): must specify a file name when using MPI\n");
    exit(1);
  }
  char name[strlen(file) + 2];
  strcpy (name, file);
  FILE * fh = fopen (name, "w");
  if (fh == NULL) {
    perror (name);
    exit (1);
  }

  scalar * dlist = dump_list (p.list ? p.list : all);
  scalar size[];
  scalar * list = list_concat ({size}, dlist); free (dlist);
  struct DumpHeader header = { t, list_len(list), iter, depth(), npe(),
                               dump_version };
#if MULTIGRID_MPI
  for (int i = 0; i < dimension; i++)
    (&header.n.x)[i] = mpi_dims[i];
  MPI_Barrier (MPI_COMM_WORLD);
#endif
  if (pid() == 0)
    dump_header (fh, &header, list);
  scalar index = {-1};
  index = new scalar;
  z_indexing (index, false);
  int cell_size = sizeof(unsigned) + header.len*sizeof(double);
  int sizeofheader = sizeof(header) + 4*sizeof(double);
  for (scalar s in list)
    sizeofheader += sizeof(unsigned) + sizeof(char)*strlen(s.name);
  subtree_size (size, false);
  long offset = 0;
  foreach_cell() {
    if (is_local(cell)) {
      if (offset == (long)0)
        offset = sizeofheader + index[]*cell_size;
      break;
    }
  }
  fseek (fh, offset, SEEK_SET); //Once
  foreach_cell() {
    if (is_local(cell)) {
      unsigned flags = is_leaf(cell) ? leaf : 0;
      fwrite (&flags, 1, sizeof(unsigned), fh);
      for (scalar s in list)
        fwrite (&s[], 1, sizeof(double), fh);
    }
    if (is_leaf(cell))
      continue;
  }
  delete ({index});
  free (list);
  fclose (fh);
} // End of dump3

event initial_vtr (i = 1100)
{
  dump3 (file = "init");//Save a snapshot file we can restart from
}
//save vtr
event vtr (i = 1100; i=i+50 ) {
  if( t > 0.10 && t < 0.11 ) {
    scalar *list = {f,p,u.x,u.y,u.z};
    vector *vlist = {NULL};
    backup_fields_vtr (list, vlist, i);
  }
  else if( t > 0.11 )
    return 1;
}

//event after_connor ( t = end ) {
//  dump3 (file = "final");//Save a snapshot file we can restart from
//  scalar *list = {f,p,u.x,u.y,u.z};
//  vector *vlist = {NULL};
//  backup_fields_vtr (list, vlist, i);
//}
