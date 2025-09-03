#include <iostream>
#include <cmath>
#include <fftw3.h>
#include <fstream>
#include <complex>
#include <string>
#include <random>
#include <mpi.h>

// these are LAMMPS include files
#include "lammps.h"
#include "input.h"
#include "atom.h"
#include "atom_vec_dipole.h"
#include "library.h"

using namespace LAMMPS_NS;

//#include <boost/timer.hpp>

using namespace std;

#define PI 3.14159265358979323846


double Lx, timestep, tfinal, Dc, alpha, c0, bf, kd, km;
long long int Nsteps;

int Ntotal, N, output_steps, read_cfile;
string op_dir, cfile_name;

int main(int narg, char** arg)
{
        
    MPI_Init(&narg,&arg);
    
    // Basic setup and checks
    string infile_name;
    if (narg == 4) {
        infile_name = std::string(arg[3]);
        
    }
    else {
        std::cout << "Format is main P in.lmp input.dat" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    
    int me,nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD,&me);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    
    int nprocs_lammps = atoi(arg[1]);
    if (nprocs_lammps > nprocs) {
        if (me == 0)
            printf("ERROR: LAMMPS cannot use more procs than available\n");
        MPI_Abort(MPI_COMM_WORLD,1);
    }
    
    int lammps;
    if (me < nprocs_lammps) lammps = 1;
    else lammps = MPI_UNDEFINED;
    MPI_Comm comm_lammps;
    MPI_Comm_split(MPI_COMM_WORLD,lammps,0,&comm_lammps);
    
    // open LAMMPS input script

    FILE *fp;
    if (me == 0) {
        fp = fopen(arg[2],"r");
        if (fp == NULL) {
            printf("ERROR: Could not open LAMMPS input script\n");
            MPI_Abort(MPI_COMM_WORLD,1);
        }
    }
    
    
  
    LAMMPS *lmp = NULL;
    if (lammps == 1) lmp = new LAMMPS(0,NULL,comm_lammps);

    int n;
    char line[1024];
    while (1) {
        if (me == 0) {
            if (fgets(line,1024,fp) == NULL) n = 0;
            else n = strlen(line) + 1;
            if (n == 0) fclose(fp);
        }
        MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
        if (n == 0) break;
        MPI_Bcast(line,n,MPI_CHAR,0,MPI_COMM_WORLD);
        if (lammps == 1) lammps_command(lmp,line);
    }
    
    
    //////////////////////////////////////////////////////
    ///// INPUT FROM DAT FILE ////////////////////////////
    
    ifstream infile(infile_name);
    
    infile >> Lx;
    infile >> N;
    infile >> read_cfile;
    infile >> cfile_name;
    infile >> Dc;
    infile >> alpha;
    infile >> bf;
    infile >> c0;
    infile >> kd;
    infile >> km;
    infile >> op_dir;
    infile >> output_steps;
    infile >> tfinal;
    infile >> timestep;

    Ntotal = lammps_get_natoms(lmp);
    Nsteps = (long long int)(tfinal / timestep);
  
    std::cout << "Ntotal = " << Ntotal << std::endl;
    std::cout << "Timestep = " << 10000*timestep << std::endl;

    string incfile = cfile_name;
    ifstream cfile(incfile);

    string opfile1 = op_dir + "/" + "particles.xyz";
    ofstream ofile1(opfile1);

    string opfile2 = op_dir + "/" + "chem_conc.xyz";
    ofstream ofile2(opfile2);

    //////////////////////////////////////////////////////
    ////// Allocate memory and variables /////////////////
    
    double *x = new double[3*Ntotal];
    double *v = new double[3*Ntotal];
    int *type = new int[Ntotal];
    double *bias = new double[Ntotal];
    double *c = new double[N];
    double *dc = new double[N];
    double *s = new double[N];    
    double *xgrid = new double[N];
    
    double dx = Lx / (N-1);
    double dx2 = dx*dx;
    double coefd = Dc/dx2;
    
    int ind1, ind2;
    int jl, jr;
        
    ///////// INITIAL CONDITION //////////////////////////
    
    for (int i = 0; i < N; i++) {
        
	if (read_cfile == 0)    
	    c[i] = c0;
        
	xgrid[i] = i*dx;
    }

    if (read_cfile == 1) {
	cfile >> N;
	for (int i = 0; i < N; i++)
	    cfile >> c[i];	
    }
        
    lammps_gather_atoms(lmp,(char *) "x",1,3,x);
    lammps_gather_atoms(lmp,(char *) "v",1,3,v);
    lammps_gather_atoms(lmp,(char *) "type",0,1,type);
    lammps_gather_atoms(lmp,(char *) "rmass",1,1,bias);
   
    MPI_Barrier(MPI_COMM_WORLD);

    //////////// BEGIN TIMESTEPPING /////////////////////////
    
    for (int i = 0; i < Nsteps; i++) {
        
        // get coordinates from lammps
        
	lammps_gather_atoms(lmp,(char *) "x",1,3,x);
	lammps_gather_atoms(lmp,(char *) "v",1,3,v);
        lammps_gather_atoms(lmp,(char *) "type",0,1,type);
	lammps_gather_atoms(lmp,(char *) "rmass",1,1,bias);
       
	MPI_Barrier(MPI_COMM_WORLD);

        ////////// Output /////////
        
        if (i % output_steps == 0) {
            
            if (me == 0) {
            
                std::cout << "Step number " << i << std::endl;
      
                ofile1 << Ntotal << "\n"; 
                ofile1 << "Time in s = " << i*timestep << "\n";
                for (int l = 0; l < Ntotal; l++)
                    ofile1 << type[l] << " " << l << " " << x[3*l] << " " << x[3*l+1] << " " << x[3*l+2] << " " << v[3*l] << " " << v[3*l+1] << " " << v[3*l+2] << " " << bias[l] << " 0.5 " << "\n";

		ofile2 << N << "\n";
		ofile2 << "Time in s = " << i*timestep << "\n";
		
		for (int l = 0; l < N; l++)
		    ofile2 << c[l] << " ";
		ofile2 << "\n";
            
		for (int l = 0; l < N; l++)
                    ofile2 << dc[l] << " ";
                ofile2 << "\n";
		
                }
        }
       

	MPI_Barrier(MPI_COMM_WORLD);
	
        if (me == 0) {
            
            ////// Add the source term //////
            for (int j = 0; j < N; j++)
                s[j] = 0.0;
            
            for (int l = 0; l < Ntotal; l++) {
                
                if (type[l] == 1) {
                    
                    ind1 = floor(x[3*l]/dx);
                    ind2 = ind1+1;
                    
                    s[ind1] += 0.5*c[ind1]/(c[ind1] + km);
                    s[ind2] += 0.5*c[ind2]/(c[ind2] + km);
                    
                }
            }
            
            // solve for the concentration field /////            
            for (int j = 0; j < N; j++) {
                
		if (j == N-1) c[j] = c0;
		else {
                    if (j == 0) jl = 1;
                    else jl = j-1;
                    jr = j+1;
                    c[j] += timestep*(coefd*(c[jr]+c[jl]-2.0*c[j]) - alpha*s[j]);
		}
            }
            
            // gradient of concentration field
            for (int j = 0; j < N; j++) {
            
                if (j == 0)
                    dc[j] = 0.0;
                else if (j == N-1)
                    dc[j] = 0.0;
                else 
                    dc[j] = 0.5*(c[j+1] - c[j-1])/dx;
            }
            
            // compute the bias
            for (int l = 0; l < Ntotal; l++) {
                
                if (type[l] == 1) {
                 
                    ind1 = floor(x[3*l]/dx);
                    ind2 = ind1+1;
		  
		    bias[l] = 0.50 + bf*(c[ind2]/(c[ind2] + kd) - c[ind1]/(c[ind1] + kd));

                    if (bias[l] > 0.9) bias[l] = 0.9;
		    if (bias[l] < 0.1) bias[l] = 0.1;
          	    	    
                }
            }    
        
        }
        
        MPI_Barrier(MPI_COMM_WORLD);
        
        MPI_Bcast(bias,Ntotal,MPI_DOUBLE,0,MPI_COMM_WORLD);
        
        MPI_Barrier(MPI_COMM_WORLD);
        
        // transmit bias to lammps
        lammps_scatter_atoms(lmp,(char *) "rmass",1,1,bias);
        
	MPI_Barrier(MPI_COMM_WORLD);
	
	// run lammps
        lmp->input->one("run 1");
        
        MPI_Barrier(MPI_COMM_WORLD);
            
    }
    
    ///////////////////////////////////////////////////////
    ///////// Deallocate //////////////////////////////////
    
    delete[] c; delete [] x; delete[] type;
    delete[] s; delete[] xgrid; delete[] dc;
    delete[] bias; delete[] v;
    
    // close down LAMMPS
    
    delete lmp;
    
    // close down MPI

    if (lammps == 1) MPI_Comm_free(&comm_lammps);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    
    // close files
   
    cfile.close(); 
    infile.close();
    ofile1.close();
    ofile2.close(); 
    
    return 0;
}
    
