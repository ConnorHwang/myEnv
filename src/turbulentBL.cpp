#include "UgpCvCfIdealGas.hpp"
#include "oneDProfile.hpp"



// ========================================================
// Compressible turbulent boundary layer
// ========================================================


class TurbulentBL : public UgpCvCfIdealGas {

public:

  TurbulentBL(){
    if ( mpi_rank == 0 ) 
      cout << "TurbulentBL():" << endl;
  }

  ~TurbulentBL(){
    if ( mpi_rank == 0 ) 
      cout << "~TurbulentBL():" << endl;
  }

  int inflowTurbHook(double (*um_bc)[3], double (*Rd_bc)[3], double (*Rod_bc)[3],
           double* Tm_it, double* pm_it, FaZone * fazone) { 
   
    if ( mpi_rank == 0 )
      cout << "inflowTurbHook: " << fazone->getName() << endl; 


    // ==========================================
    //  U_ref   = C_infty
    //  L_ref   = delta (boundary layer thickness)
    //  T_ref   = T_infty
    //  rho_ref = rho_infty
    // ==========================================


    // Flow input parameters
    double Ma_infty   = getDoubleParam("Ma_INFINITY");
    double recovery_f = getDoubleParam("RECOVERY_FACTOR",0.82);
    double Twa        = 1.0 + recovery_f * 0.5*(eos->getGamma()-1.0) * Ma_infty*Ma_infty; 
    
    // Boundary layer input parameters (from DeGraff and Eaton, JFM 422 (2000))
    double Re_theta   = getDoubleParam("RE_THETA");
    double delta = getDoubleParam("DELTA");
    double theta = getDoubleParam("THETA");
    double Big_Delta = getDoubleParam("CLAUSER-ROTTA_DELTA");
    double H = getDoubleParam("H");
    // scale input BL such that BL thickness = 1 at downstream validation station 
    double BLgrowth = getDoubleParam("BL_GROWTH",1.5);

    // friction velocity utau/c_infty
    double utau = H*theta/Big_Delta*Ma_infty; 
    double cf =  utau *utau /(0.5*Ma_infty*Ma_infty);
    
    // convert to momentum and displacement thickness to non-dimensional units
    theta = theta/delta; // theta/delta
    double deltastar = H*theta; // delta*/delta

    // viscous wall unit
    double deltanu = sqrt(2.0/cf) * theta / Re_theta;

    // and viscosity
    double mu_infty =  Ma_infty * theta / ( Re_theta * pow(Twa,1.7) );

    // dump info on screen
    if ( mpi_rank == 0 ) {
      cout << "************************** BL Inflow Data ***************************" << endl;
      cout << " > Reynolds number based on momentum thickness  = " << Re_theta << endl;
      cout << " > free stream Mach number                      = " << Ma_infty << endl;
      cout << " > wall temperature / free stream temperature   = " << Twa << endl;
      cout << " > skin friction coefficient                    = " << cf << endl;
      cout << " > friction velocity                            = " << utau << endl;
      cout << " > displacement thickness / BL thickness        = " << deltastar << endl;
      cout << " > shape factor                                 = " << H << endl;
      cout << " > momentum thickness /  BL thickness           = " << theta << endl;
      cout << " > viscous wall unit / BL thickness             = " << deltanu << endl;
      cout << " > free stream viscosity                        = " << mu_infty << endl;
      cout << "*********************** end of BL Inflow Data ***********************" << endl;
    }
    

    //#########################################################
    //#########################################################
    //
    // Input profiles from DeGraff and Eaton, JFM 422 (2000)
    //
    //#########################################################
    // allocate memory
    int zone_nfa = fazone->nzf;

    double* Umean   = new double[zone_nfa];
    double* Uprime2 = new double[zone_nfa];
    double* Vprime2 = new double[zone_nfa];
    // Wprime2 not available in DeGraff and Eaton, JFM 422 (2000)
    double* UVprime = new double[zone_nfa];
    double* Tmean   = new double[zone_nfa];
    double* Pmean   = new double[zone_nfa];
    
    // get the ys
    double * y = new double[zone_nfa];
    FOR_FAZONE_IFA{
      y[izf] = x_fa[ifa][1]/deltanu * BLgrowth;
    }
    
    // mean velocity
    // and thermodynamics
    oneDProfile * UmeanP  = new oneDProfile;
    UmeanP->read("exp_data/DeGraff_JFM_fig2_yp_Up.dat");
    FOR_FAZONE_IFA{
      Umean[izf]   = utau * UmeanP->getData(y[izf]);
      Tmean[izf]   = Twa - recovery_f * 0.5 * (eos->getGamma()-1.0) * Umean[izf]*Umean[izf];  
      double myrho = 1.0 / Tmean[izf]; 
      eos->calcPFromRhoT(&(Pmean[izf]),&(myrho),&(Tmean[izf]),1);
    }
     
    // u'^2 velocity
    oneDProfile * Uprofile  = new oneDProfile;
    Uprofile->read("exp_data/DeGraff_JFM_fig7_yp_u2p.dat");
    FOR_FAZONE_IFA{
      Uprime2[izf] = utau * utau * Uprofile->getData(y[izf]);
    }
    delete Uprofile;
    
    // v'^2 velocity
    oneDProfile * Vprofile  = new oneDProfile;
    Vprofile->read("exp_data/DeGraff_JFM_fig14_yp_v2p.dat");
    FOR_FAZONE_IFA{
      Vprime2[izf] = utau * utau * Vprofile->getData(y[izf]);
    }
    delete Vprofile;
    
    // u'v' velocity
    oneDProfile * UVprofile  = new oneDProfile;
    UVprofile->read("exp_data/DeGraff_JFM_fig18_yp_uvp.dat");
    FOR_FAZONE_IFA{
      UVprime[izf] = -utau * utau * UVprofile->getData(y[izf]);
    }
    delete UVprofile;
    
    delete[] y;
    

    setCurrentParam(fazone->getName());
    string bc_name = getCurrentParamString();
    assert(bc_name == "INFLOW_TURB");
    
    FOR_FAZONE_IFA{
      pm_it[izf+fazone->izf_f]     = Pmean[izf];
      Tm_it[izf+fazone->izf_f]     = Tmean[izf];
      um_it[izf+fazone->izf_f][0]  = Umean[izf];
      um_it[izf+fazone->izf_f][1]  = 0.0;
      um_it[izf+fazone->izf_f][2]  = 0.0;
      Rd_it[izf+fazone->izf_f][0]  = Uprime2[izf];
      Rd_it[izf+fazone->izf_f][1]  = Vprime2[izf];
      Rd_it[izf+fazone->izf_f][2]  = Vprime2[izf]; // Wprime2 not available, use Vprime2;
      Rod_it[izf+fazone->izf_f][0] = UVprime[izf];
      Rod_it[izf+fazone->izf_f][1] = 0.0;
      Rod_it[izf+fazone->izf_f][2] = 0.0;
    }
    
    delete[] Umean;
    delete[] Uprime2;
    delete[] Vprime2;
    delete[] UVprime;
    delete[] Tmean;
    delete[] Pmean;
    //#########################################################

      
    return(0);
   

  }


  void temporalHook() {
    if (step%check_interval==0){
      if (mpi_rank==0){
        cout << " > Temporal Hook..." << endl;
      }
      double xvec[3] = {20.0, 30.0, 40.0};
      double z0 = 0.0;
      double tol_x = 0.10; //0.125 dx
      double tol_z = 0.05; //0.08333333 dz

      int mycount_z[3] = {0, 0, 0};

      double mytauw[3] = {0.0, 0.0, 0.0};
      FOR_FAZONE{
        if (fazone->getName() == "y0" ){
          double * tauwall = fazone->getR1("TAU_WALL");
          assert(tauwall != NULL);
          FOR_FAZONE_IFA{
            double xsim = x_fa[ifa][0];
            if (xvec[0]-tol_x < xsim &&  xsim < xvec[0]+tol_x){
              mycount_z[0] += 1;
              mytauw[0]    += tauwall[izf];

            }
            if (xvec[1]-tol_x < xsim &&  xsim < xvec[1]+tol_x){
              mycount_z[1] += 1;
              mytauw[1]    += tauwall[izf];

            }
            if (xvec[2]-tol_x < xsim &&  xsim < xvec[2]+tol_x){
              mycount_z[2] += 1;
              mytauw[2]    += tauwall[izf];

            }
          }
          //Compute spanwise average
          double tauw_spanavg[3];
          int count_z[3];
          MPI_Allreduce(mytauw,tauw_spanavg,3,MPI_DOUBLE,MPI_SUM,mpi_comm);
          MPI_Allreduce(mycount_z,count_z,3,MPI_INT,MPI_SUM,mpi_comm);

          FOR_I3 tauw_spanavg[i] /= count_z[i];
          if (mpi_rank == 0){
            FOR_I3 {cout << " > X = " << xvec[i] << ", TAUWALL = " << tauw_spanavg[i] << endl;}
          }
        }
      }
      

    }
  }


};


int main(int argc, char * argv[]) {
  
  try {
  
    // initialize the environment: basically 
    // initializes MPI and parameters...
    
    CTI_Init(argc,argv,"charles.in");
    
    // the solver...
    
    {
      
      // turbulent boundary layer
      TurbulentBL solver;
            
      solver.init1();
      solver.init2();
      solver.init3();
      
      solver.run();

    }

    // finalize the environment: reports parameter usage, shuts down MPI...
    
    CTI_Finalize();
    
    // rudimentary error management...

  }
  catch (int e) {
    if (e == 0)
      CTI_Finalize();
    else
      CTI_Abort();
  }
  catch(...) {
    CTI_Abort();
  }
    
  return(0);

}

