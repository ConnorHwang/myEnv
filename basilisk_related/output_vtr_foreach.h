//Add parallel output function
//Check if it works for rectangular domains
//Add binary capability


/*
This function writes one XML file which allows to read the *.vtr files generated
by output_vtr_ascii_foreach() when used in MPI. Tested in multigrid
using MPI.
*/
void output_pvtr_ascii (scalar * list, vector * vlist, double nxlen, int nx, int nx_proc, int ny_proc, int nz_proc, FILE * fp, char * subname)
{

	int nproc = nx/nx_proc;
    int ny = nx*ny_proc/nx_proc;
    int nz = nx*nz_proc/nx_proc;


    fputs ("<?xml version=\"1.0\"?>\n"
    "<VTKFile type=\"PRectilinearGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);
    fprintf (fp,"\t <PRectilinearGrid WholeExtent=\"%d %d %d %d %d %d \" GhostLevel=\"0\">\n", 0, nx, 0, ny, 0, nz);
    fputs ("\t\t\t <PCellData Scalars=\"scalars\">\n", fp);
    for (scalar s in list) {
      fprintf (fp,"\t\t\t\t <PDataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", s.name);
      fputs ("\t\t\t\t </PDataArray>\n", fp);
    }
    for (vector v in vlist) {
      fprintf (fp,"\t\t\t\t <PDataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"%s\" format=\"ascii\">\n", v.x.name);
      fputs ("\t\t\t\t </PDataArray>\n", fp);
    }
    fputs ("\t\t\t </PCellData>\n", fp);
  fputs ("\t\t\t <PCoordinates>\n", fp);
    fprintf (fp,"\t\t\t\t <PDataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", "x");
    fputs ("\t\t\t\t </PDataArray>\n", fp);
    fprintf (fp,"\t\t\t\t <PDataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", "y");
    fputs ("\t\t\t\t </PDataArray>\n", fp);
    fprintf (fp,"\t\t\t\t <PDataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", "z");
    fputs ("\t\t\t\t </PDataArray>\n", fp);
  fputs ("\t\t\t </PCoordinates>\n", fp);

    for (int i=0; i<nx_proc; i++) {
    	for (int j=0; j<ny_proc; j++) {
    		for (int k=0; k<nz_proc; k++) {
      			fprintf (fp, "<Piece Extent=\"%d %d %d %d %d %d \" Source=\"%s_n%3.3d.vtr\"/> \n", i*nproc+0,i*nproc+nproc, j*nproc+0, j*nproc+nproc, k*nproc+0, k*nz_proc+nproc, subname, i*nz_proc*ny_proc + j*nz_proc + k);
		}
	}
    }

   // for (int i = 0; i < npe(); i++)
   //   fprintf (fp, "<Piece Extent=\"%d %d %d %d %d %d \" Source=\"%s_n%06d.vtr\"/> \n", 0, nproc, 0, nproc, 0, nproc, subname, i);

    fputs ("\t </PRectilinearGrid>\n", fp);
    fputs ("</VTKFile>\n", fp);
}

/*
This function writes one XML VTK file per PID process of type rectilinear grid
(*.vtr) which can be read using Paraview. File stores scalar and vector fields
defined at the center points. Results are recorded on ASCII format. If one writes
one *.vtr file per PID process this function may be combined with
output_pvtr_ascii() above to read in parallel. Tested in multigrid
using MPI. 
*/
void output_vtr_ascii_foreach (scalar * list, vector * vlist, double nxlen, int nx, int nx_proc, int ny_proc, int nz_proc, FILE * fp, bool linear)
{
#if defined(_OPENMP)
  int num_omp = omp_get_max_threads();
  omp_set_num_threads(1);
#endif

	double dx = nxlen/nx;
	int nproc = nx/nx_proc;
	int i = pid()/(nz_proc*ny_proc);
	int j = (pid() - i*nz_proc*ny_proc)/nz_proc;
	int k = (pid() - i*nz_proc*ny_proc - j*nz_proc);


  fputs ("<?xml version=\"1.0\"?>\n"
  "<VTKFile type=\"RectilinearGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);
  fprintf (fp,"\t <RectilinearGrid WholeExtent=\"%d %d %d %d %d %d \">\n", i*nproc+0,i*nproc+nproc, j*nproc+0, j*nproc+nproc, k*nproc+0, k*nproc+nproc);
  fprintf (fp,"\t\t <Piece Extent=\"%d %d %d %d %d %d \">\n", i*nproc+0,i*nproc+nproc, j*nproc+0, j*nproc+nproc, k*nproc+0, k*nproc+nproc);
  fputs ("\t\t\t <CellData Scalars=\"scalars\">\n", fp);
  for (scalar s in list) {
    fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", s.name);
	int p=0;
	int l=0;
	int m=0;
	int n=0;
	double scalar[nproc][nproc][nproc];
    foreach(){
	l = p/(nproc*nproc);
	m = (p - l*nproc*nproc)/nproc;
	n = (p - l*nproc*nproc - m*nproc);
	scalar[l][m][n] = val(s);
	p = p+1;
    }

	for (int n=0;n<nproc;n++) {
		for (int m=0;m<nproc;m++) {
			for (int l=0;l<nproc;l++) {
   				fprintf (fp, "\t\t\t\t\t %g\n", scalar[l][m][n]);
			}
		}
	}


    fputs ("\t\t\t\t </DataArray>\n", fp);
  }
  for (vector v in vlist) {
    fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"%s\" format=\"ascii\">\n", v.x.name);
    foreach(){
#if dimension == 2
      fprintf (fp, "\t\t\t\t\t %g %g 0.\n", val(v.x), val(v.y));
#endif
#if dimension == 3
      fprintf (fp, "\t\t\t\t\t %g %g %g\n", val(v.x), val(v.y), val(v.z));
#endif
    }
    fputs ("\t\t\t\t </DataArray>\n", fp);
  }
  fputs ("\t\t\t </CellData>\n", fp);
  fputs ("\t\t\t <Coordinates>\n", fp);
    fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", "x");
	for(int l=0;l<=nproc;l++)
      		fprintf (fp, "\t\t\t\t\t %g\n",i*nproc*dx + l*dx);
    fputs ("\t\t\t\t </DataArray>\n", fp);
    fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", "y");
	for(int l=0;l<=nproc;l++)
      		fprintf (fp, "\t\t\t\t\t %g\n",j*nproc*dx + l*dx);
    fputs ("\t\t\t\t </DataArray>\n", fp);
    fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", "z");
	for(int l=0;l<=nproc;l++)
      		fprintf (fp, "\t\t\t\t\t %g\n",k*nproc*dx + l*dx);
    fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t </Coordinates>\n", fp);
  fputs ("\t\t </Piece>\n", fp);
  fputs ("\t </RectilinearGrid>\n", fp);
  fputs ("</VTKFile>\n", fp);
  fflush (fp);
#if defined(_OPENMP)
  omp_set_num_threads(num_omp);
#endif
}

///*
//This function writes one XML file which allows to read the *.vtu files generated
//by output_vtu_bin_foreach() when used in MPI. Tested in (quad- and oct-)trees
//using MPI.
//*/
//void output_pvtu_bin (scalar * list, vector * vlist, int n, FILE * fp, char * subname)
//{
//    fputs ("<?xml version=\"1.0\"?>\n"
//    "<VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);
//    fputs ("\t <PUnstructuredGrid GhostLevel=\"0\">\n", fp);
//    fputs ("\t\t\t <PCellData Scalars=\"scalars\">\n", fp);
//    for (scalar s in list) {
//      fprintf (fp,"\t\t\t\t <PDataArray type=\"Float64\" Name=\"%s\" format=\"appended\">\n", s.name);
//      fputs ("\t\t\t\t </PDataArray>\n", fp);
//    }
//    for (vector v in vlist) {
//      fprintf (fp,"\t\t\t\t <PDataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"%s\" format=\"appended\">\n", v.x.name);
//      fputs ("\t\t\t\t </PDataArray>\n", fp);
//    }
//    fputs ("\t\t\t </PCellData>\n", fp);
//    fputs ("\t\t\t <PPoints>\n", fp);
//    fputs ("\t\t\t\t <PDataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n", fp);
//    fputs ("\t\t\t\t </PDataArray>\n", fp);
//    fputs ("\t\t\t </PPoints>\n", fp);
//
//    for (int i = 0; i < npe(); i++)
//      fprintf (fp, "<Piece Source=\"%s_n%3.3d.vtu\"/> \n", subname, i);
//
//    fputs ("\t </PUnstructuredGrid>\n", fp);
//    fputs ("</VTKFile>\n", fp);
//}
//
///*
//This function writes one XML VTK file per PID process of type unstructured grid
//(*.vtu) which can be read using Paraview. File stores scalar and vector fields
//defined at the center points. Results are recorded on binary format. If one writes
//one *.vtu file per PID process this function may be combined with
//output_pvtu_bin() above to read in parallel. Tested in (quad- and oct-)trees
//using MPI. Also works with solids (when not using MPI).
//*/
//void output_vtu_bin_foreach (scalar * list, vector * vlist, int n, FILE * fp, bool linear)
//{
//#if defined(_OPENMP)
//  int num_omp = omp_get_max_threads();
//  omp_set_num_threads(1);
//#endif
//
//  vertex scalar marker[];
//  int no_points = 0, no_cells=0 ;
//  foreach_vertex(){
//    marker[] = _k;
//    no_points += 1;
//  }
//  foreach(){
//    no_cells += 1;
//  }
//
//  fputs ("<?xml version=\"1.0\"?>\n"
//  "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);
//  fputs ("\t <UnstructuredGrid>\n", fp);
//  fprintf (fp,"\t\t <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", no_points, no_cells);
//  fputs ("\t\t\t <CellData Scalars=\"scalars\">\n", fp);
//
//  int count = 0;
//  for (scalar s in list) {
//    fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" Name=\"%s\" format=\"appended\" offset=\"%d\">\n", s.name,count);
//    count += ((no_cells)+1)*8;
//    fputs ("\t\t\t\t </DataArray>\n", fp);
//  }
//  for (vector v in vlist) {
//    fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"3\"  format=\"appended\" offset=\"%d\">\n", v.x.name,count);
//    count += ((no_cells*3)+1)*8;
//    fputs ("\t\t\t\t </DataArray>\n", fp);
//  }
//  fputs ("\t\t\t </CellData>\n", fp);
//  fputs ("\t\t\t <Points>\n", fp);
//  fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\"  format=\"appended\" offset=\"%d\">\n",count);
//  count += ((no_points*3)+1)*8;
//  fputs ("\t\t\t\t </DataArray>\n", fp);
//  fputs ("\t\t\t </Points>\n", fp);
//  fputs ("\t\t\t <Cells>\n", fp);
//  fputs ("\t\t\t\t <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n", fp);
//  foreach(){
//#if dimension == 2
//    fprintf (fp, "%g %g %g %g \n", marker[], marker[1,0], marker[1,1], marker[0,1]);
//#endif
//#if dimension == 3
//    fprintf (fp, "%g %g %g %g %g %g %g %g\n", marker[], marker[1,0,0], marker[1,1,0], marker[0,1,0],marker[0,0,1], marker[1,0,1], marker[1,1,1], marker[0,1,1]);
//#endif
//  }
//  fputs ("\t\t\t\t </DataArray>\n", fp);
//  fputs ("\t\t\t\t <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n", fp);
//  for (int i = 1; i < no_cells+1; i++){
//#if dimension == 2
//    fprintf (fp, "%d \n", i*4);
//#endif
//#if dimension == 3
//    fprintf (fp, "%d \n", i*8);
//#endif
//  }
//  fputs ("\t\t\t\t </DataArray>\n", fp);
//  fputs ("\t\t\t\t <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n", fp);
//  foreach(){
//#if dimension == 2
//    fputs ("9 \n", fp);
//#endif
//#if dimension == 3
//    fputs ("12 \n", fp);
//#endif
//  }
//  fputs ("\t\t\t\t </DataArray>\n", fp);
//  fputs ("\t\t\t </Cells>\n", fp);
//  fputs ("\t\t </Piece>\n", fp);
//  fputs ("\t </UnstructuredGrid>\n", fp);
//  fputs ("\t <AppendedData encoding=\"raw\">\n", fp);
//  fputs ("_", fp);
//  unsigned long long block_len=no_cells*8;
//#if dimension == 2
//  double z=0, vz=0;
//#endif
//  for (scalar s in list) {
//    fwrite (&block_len, sizeof (unsigned long long), 1, fp);
//    foreach()
//      fwrite (&val(s), sizeof (double), 1, fp);
//  }
//  block_len=no_cells*8*3;
//  for (vector v in vlist) {
//    fwrite (&block_len, sizeof (unsigned long long), 1, fp);
//    foreach(){
//      fwrite (&val(v.x), sizeof (double), 1, fp);
//      fwrite (&val(v.y), sizeof (double), 1, fp);
//#if dimension == 2
//      fwrite (&vz, sizeof (double), 1, fp);
//#endif
//#if dimension == 3
//      fwrite (&val(v.z), sizeof (double), 1, fp);
//#endif
//    }
//  }
//  block_len=no_points*8*3;
//  fwrite (&block_len, sizeof (unsigned long long), 1, fp);
//  foreach_vertex(){
//    fwrite (&x, sizeof (double), 1, fp);
//    fwrite (&y, sizeof (double), 1, fp);
//    fwrite (&z, sizeof (double), 1, fp);
//  }
//  fputs ("\t\n", fp);
//  fputs ("\t </AppendedData>\n", fp);
//  fputs ("</VTKFile>\n", fp);
//  fflush (fp);
//#if defined(_OPENMP)
//  omp_set_num_threads(num_omp);
//#endif
//}
