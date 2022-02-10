double getData (double the, double r, const double redge, const int Nthe, const int Nr, double *inputU) {
  //int the1 = Nthe, the2 = 0, r1 = Nr, r2 = 0;
  int the1 = 0, the2 = 0, r1 = 0, r2 = 0;
  double val1, val2, val3, val4, wthe, wr;
  
  for (int i = 0; i < Nthe; i++) {
    if ( the <= i*2.0*M_PI/Nthe ) {
      the1 = i-1; the2 = i;
      for (int j = 0; j < Nr; j++) {
        if ( r <= j*redge/Nr ) {
          r1 = j-1; r2 = j;
      
          val1 = inputU[the1*Nr + r1];
          val2 = inputU[the2*Nr + r1];
          val3 = inputU[the1*Nr + r2];
          val4 = inputU[the2*Nr + r2];
          
          wthe = (the - the1*2.0*M_PI/Nthe) / ((the2 - the1)*2.0*M_PI/Nthe);
          wr   = (r - r1*redge/Nr) / (r2*redge/Nr - r1*redge/Nr);

          return (1.0-wthe)*(wr*val3 + (1.0-wr)*val1) + wthe*(wr*val4 + (1.0-wr)*val2);
        } else if ( r > (Nr-1)*redge/Nr ) {
          return 0;
        }
      }
    } else if ( the > (Nthe-1)*2.0*M_PI/Nthe ) {
      the1 = Nthe-1; the2 = 0;
      for (int j = 0; j < Nr; j++) {
        if ( r <= j*redge/Nr ) {
          r1 = j-1; r2 = j;
      
          val1 = inputU[the1*Nr + r1];
          val2 = inputU[the2*Nr + r1];
          val3 = inputU[the1*Nr + r2];
          val4 = inputU[the2*Nr + r2];
          
          wthe = (the - the1*2.0*M_PI/Nthe) / ((the2 - the1)*2.0*M_PI/Nthe);
          wr   = (r - r1*redge/Nr) / (r2*redge/Nr - r1*redge/Nr);

          return (1.0-wthe)*(wr*val3 + (1.0-wr)*val1) + wthe*(wr*val4 + (1.0-wr)*val2);
        } else if ( r > (Nr-1)*redge/Nr ) {
          return 0;
        }
      }
    }
  }
  fprintf(stderr, "Data leaked in <getData> function!\n");
  return 1;
}

void readData(FILE * fp, const int Nthe, const int Nr, double *inputU) {
  fprintf (stderr, "Reading the file...\n");
  char buf[100];
  for (int i = 0; i < Nthe; i++) {
    for (int j = 0; j < Nr; j++) {
      fscanf(fp, "%s", buf);
      inputU[i*Nr + j] = atof(buf);
    }
  }
  fprintf (stderr, "Finished reading input file\n");
}
//struct inputinfo {
//  //compulsory
//  FILE * fp;
//  //optional
//  int Nthe, Nr;
//};
