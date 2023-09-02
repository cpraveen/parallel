//______________________________________________________________________________
// This program (parallelized with MPI) solves 2D advection equation on a 
// square doubly periodic domain.
//______________________________________________________________________________
// Author : Pritam Giri
// Date   : 2.10.2016
// TIFR-CAM
//______________________________________________________________________________
#include<iostream>
#include<cstdlib>
#include<fstream>
#include<ctype.h>
#include<cstring>
#include<math.h>
#include<ctime>
#include<sstream>
#include <iomanip> 
#include <mpi.h>

using namespace std;

const double PI = 4.0*atan(1.0);
const double ark[] = {0.0, 3.0/4.0, 1.0/3.0};

//______________________________________________________________________________
double initial_condition(const double x, const double y)
{
  return sin(2.0*PI*x) * sin(2.0*PI*y);
}

//______________________________________________________________________________
// Absolute value
//______________________________________________________________________________
template <class Type> 
Type Absolute(Type a)
{
	return (a > (Type)(0.0) ? a : -a);
}

//______________________________________________________________________________
// Minimum value
//______________________________________________________________________________
template <class Type> 
Type Minimum(Type a, Type b)
{
	return (a > b ? b : a);
}

//______________________________________________________________________________
// Print an error message and exit the program
//______________________________________________________________________________
void ErrorMessage(const char* message)
{
	cout << message << endl;
	exit(1);
}

//______________________________________________________________________________
// Weno reconstruction
// Return left value for face between u0, up1
// Adapted from: 
// https://github.com/cpraveen/cfdlab/blob/master/petsc/convect2d/convect.c
//______________________________________________________________________________
double weno5(double um2, double um1, double u0, double up1, double up2)
{
	double eps = 1.0e-6;
	double gamma1=1.0/10.0, gamma2=3.0/5.0, gamma3=3.0/10.0;
	double beta1, beta2, beta3;
	double u1, u2, u3;
	double w1, w2, w3;
	
	beta1 = (13.0/12.0)*pow((um2 - 2.0*um1 + u0),2) 
          + (1.0/4.0)*pow((um2 - 4.0*um1 + 3.0*u0),2);
	beta2 = (13.0/12.0)*pow((um1 - 2.0*u0 + up1),2) 
          + (1.0/4.0)*pow((um1 - up1),2);
	beta3 = (13.0/12.0)*pow((u0 - 2.0*up1 + up2),2) 
          + (1.0/4.0)*pow((3.0*u0 - 4.0*up1 + up2),2);
	
	w1 = gamma1 / pow(eps+beta1, 2);
	w2 = gamma2 / pow(eps+beta2, 2);
	w3 = gamma3 / pow(eps+beta3, 2);
	
	u1 = (1.0/3.0)*um2 - (7.0/6.0)*um1 + (11.0/6.0)*u0;
	u2 = -(1.0/6.0)*um1 + (5.0/6.0)*u0 + (1.0/3.0)*up1;
	u3 = (1.0/3.0)*u0 + (5.0/6.0)*up1 - (1.0/6.0)*up2;
	
	return (w1 * u1 + w2 * u2 + w3 * u3)/(w1 + w2 + w3);
}

//______________________________________________________________________________
// Allocate a 1D array
//______________________________________________________________________________
template <class Type> 
void Allocate(Type *&m, int d1)
{
	m = new (nothrow) Type [d1];
	if (m == 0)
		ErrorMessage("Error: Memory can not be allocated.");

	for (int i = 0; i < d1; i++)
		m[i] = (Type)(0.0);
}

//______________________________________________________________________________
// Deallocate a 1D array
//______________________________________________________________________________
template <class Type> 
void Deallocate(Type *m, int d1)
{
	delete [] m;
	m = NULL;
}

//______________________________________________________________________________
// Allocate a 2D array
//______________________________________________________________________________
template <class Type> 
void Allocate(Type **&m, int d1, int d2)
{
	m = new (nothrow) Type* [d1];
	
	if (m == 0)
		ErrorMessage("Error: Memory can not be allocated.");

	for (int i = 0; i < d1; i++)
	{
		m[i] = new (nothrow) Type [d2];
	
		if (m[i] == 0)
	   		ErrorMessage("Error: Memory can not be allocated.");
		
		for (int j = 0; j < d2; j++)
			m[i][j] = (Type)(0.0);
	}
}

//______________________________________________________________________________
// Deallocate a 2D array
//______________________________________________________________________________
template <class Type> 
void Deallocate(Type **m, int d1, int d2)
{
	for (int i = 0; i < d1; i++)
		delete [] m[i];

	delete [] m;
	m = NULL;
}

//______________________________________________________________________________
// Allocate a 2D matrix in contiguous memory locations. 
// Indexing starts with m[-xoffset][-yoffset].
//______________________________________________________________________________
template <class Type> 
void Allocate(Type **&m, Type *&mcontiguous, const int d1, const int d2, 
              const int xoffset = 0, const int yoffset = 0)
{
	Allocate(mcontiguous,d1*d2);
	
	m = new (nothrow) Type* [d1];
	
	if (m == 0)
		ErrorMessage("Error: Memory can not be allocated.");
		
	for (unsigned int i = 0; i < d1; i++)
	{
		m[i] = &(mcontiguous[d2*i]);
				
		for (unsigned int j = 0; j < d2; j++)
			m[i][j] = (Type)(0.0);
	}

	m += xoffset;

	for (int i = -xoffset; i < (d1-xoffset); i++)
		m[i] += yoffset;
}

//______________________________________________________________________________
// Deallocate a 2D matrix which has been allocated at contiguous memory 
// locations. 
// Indexing starts with m[-xoffset][-yoffset].
//______________________________________________________________________________
template <class Type> 
void Deallocate(Type **m, Type *mcontiguous, const int d1, const int d2, 
                const int xoffset = 0, const int yoffset = 0)
{
	for (int i = -xoffset; i < (d1-xoffset); i++)
		m[i] -= yoffset;
	
	m -= xoffset;
	
	delete [] m;
	m = NULL;

	Deallocate(mcontiguous,d1*d2);
}

//______________________________________________________________________________
// Periodic wrap-around neighbor search.
//______________________________________________________________________________
void SearchNeighbor(int &TopNeighbor, int &BottomNeighbor, int &LeftNeighbor, 
                    int &RightNeighbor, int &TopRightNeighbor, 
                    int &BottomLeftNeighbor, const int I, const int J, 
                    const int XDim, const int YDim)
{
	int Im = (I == 0 ? XDim-1 : I-1);
	int Ip = (I == XDim-1 ? 0 : I+1);
	int Jm = (J == 0 ? YDim-1 : J-1);
	int Jp = (J == YDim-1 ? 0 : J+1);
	
	TopNeighbor = Im*YDim+J;
	BottomNeighbor = Ip*YDim+J;

	LeftNeighbor  = I*YDim+Jm;
	RightNeighbor = I*YDim+Jp;

	TopRightNeighbor = Im*YDim+Jm;
	BottomLeftNeighbor = Ip*YDim+Jp;
}

//______________________________________________________________________________
int main(int argc, char** argv)
{
	// Initialize the MPI environment
	MPI_Init(&argc,&argv);

	// Get the number of processes
	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	// Get the rank of the process
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// Get the name of the processor
	int name_len;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	MPI_Get_processor_name(processor_name, &name_len);
	
	time_t ti, tf ;
	cout.flags( ios::dec | ios::fixed );
	cout.precision(8);
	ti = time(NULL);
	
	// XDim    : Number of processors along x direction,
	// YDim    : Number of processors along y direction,
	// Nx      : Number of cells along x direction in each processor,
	// Ny      : Number of cells along y direction in each processor,
	// XOffset : Number of layers needed from neighbors for high order reconstruction in x direction,
	// YOffset : Number of layers needed from neighbors for high order reconstruction in y direction.
	
	// Note that global size of the domain (XDim*Nx) x (YDim*Ny)
	// Number of processors required to run this code is (XDim*YDim)
	
	const int XDim = 2, YDim = 2, Nx = 50, Ny = 50, XOffset = 3, YOffset = 3, 
            XOffsetEnd = 2, YOffsetEnd = 2, SaveInterval = 100;
	const int I = rank/YDim, J = (rank - I*YDim);
	const double xmin = 0.0, xmax = 1.0, ymin = 0.0, ymax = 1.0, cfl = 0.4, 
               Tf = 10.0;
	const double dx = (xmax-xmin)/(XDim*Nx), dy = (ymax-ymin)/(YDim*Ny);

	double **Phi, *Phicontiguous, **PhiOld, **Residual;
	double t = 0.0, umax = sqrt(2.0), dt = cfl*Minimum(dx,dy)/umax, 
         lambda = dt/(dx*dy), flux;
	
	int TopNeighbor, BottomNeighbor, LeftNeighbor, RightNeighbor, 
      TopRightNeighbor, BottomLeftNeighbor;
	
	SearchNeighbor(TopNeighbor,BottomNeighbor,LeftNeighbor,RightNeighbor,
                 TopRightNeighbor,BottomLeftNeighbor,I,J,XDim,YDim);
	
	if (size != (XDim*YDim))
	{
		if (rank == 0)
			fprintf(stderr,"%d processors are needed to run this code!",XDim*YDim);
		
		MPI_Abort(MPI_COMM_WORLD,1);
	}

	Allocate(Phi,Phicontiguous,Nx+XOffset+XOffsetEnd,Ny+YOffset+YOffsetEnd,
           XOffset,YOffset);
	Allocate(PhiOld,Nx,Ny);
	Allocate(Residual,Nx,Ny);
	
	// Initialize.
	for (unsigned int i = 0; i < Nx; i++)
	{
		for (unsigned int j = 0; j < Ny; j++)
		{
			double x = xmin + dx*(I*Nx+i-0.5);
			double y = ymin + dy*(J*Ny+j-0.5);
			Phi[i][j] = initial_condition(x, y);
		}
	}
	
	// Time stepping.
	int iteration = 0;
	
	MPI_Datatype PackedColumns0, PackedColumns1;
	MPI_Status ReceiveStatus;

	MPI_Type_vector(Nx,YOffset,Ny+YOffset+YOffsetEnd,MPI_DOUBLE,&PackedColumns0);
	MPI_Type_commit(&PackedColumns0);

	MPI_Type_vector(Nx,YOffsetEnd,Ny+YOffset+YOffsetEnd,MPI_DOUBLE,
                  &PackedColumns1);
	MPI_Type_commit(&PackedColumns1);
	
	while (t < Tf)
	{
		if (t+dt > Tf)
		{
			dt = Tf - t;
			lambda = dt/(dx*dy);
		}
		
		for (unsigned int i = 0; i < Nx; i++)
			for (unsigned int j = 0; j < Ny; j++)
				PhiOld[i][j] = Phi[i][j];

		for (int rk = 0; rk < 3; rk++)
		{
			for (unsigned int i = 0; i < Nx; i++)
				for (unsigned int j = 0; j < Ny; j++)
					Residual[i][j] = 0.0;

			// Copy required portions from the neighbors.
			// message passing is needed only for multiple processors.
			if (TopNeighbor == rank)
			{
				for (int i = 0; i < XOffsetEnd; i++)
					for (int j = 0; j < Ny; j++)
						Phi[Nx+i][j] = Phi[i][j];

				for (int i = 0; i < XOffset; i++)
					for (int j = 0; j < Ny; j++)
						Phi[i-XOffset][j] = Phi[Nx-XOffset+i][j];
			}
			else
			{
				MPI_Send(&Phi[0][-YOffset] , XOffsetEnd*(Ny+YOffset+YOffsetEnd), 
                 MPI_DOUBLE, TopNeighbor   , 0, MPI_COMM_WORLD);
				MPI_Recv(&Phi[Nx][-YOffset], XOffsetEnd*(Ny+YOffset+YOffsetEnd), 
                 MPI_DOUBLE, BottomNeighbor, 0, MPI_COMM_WORLD, &ReceiveStatus);
				
				MPI_Send(&Phi[Nx-XOffset][-YOffset], XOffset*(Ny+YOffset+YOffsetEnd), 
                 MPI_DOUBLE, BottomNeighbor, 1, MPI_COMM_WORLD);
				MPI_Recv(&Phi[-XOffset][-YOffset]  , XOffset*(Ny+YOffset+YOffsetEnd), 
                 MPI_DOUBLE, TopNeighbor   , 1, MPI_COMM_WORLD, &ReceiveStatus);
			}

			if (LeftNeighbor == rank)
			{
				for (int i = 0; i < Nx; i++)
					for (int j = 0; j < YOffsetEnd; j++)
						Phi[i][Ny+j] = Phi[i][j];

				for (int i = 0; i < Nx; i++)
					for (int j = 0; j < YOffset; j++)
						Phi[i][j-YOffset] = Phi[i][Ny-YOffset+j];
			}
			else
			{
				MPI_Send(&Phi[0][0] , 1, PackedColumns1, LeftNeighbor , 2, 
                 MPI_COMM_WORLD);
				MPI_Recv(&Phi[0][Ny], 1, PackedColumns1, RightNeighbor, 2, 
                 MPI_COMM_WORLD, &ReceiveStatus);

				MPI_Send(&Phi[0][Ny-YOffset], 1, PackedColumns0, RightNeighbor, 3, 
                 MPI_COMM_WORLD);
				MPI_Recv(&Phi[0][-YOffset]  , 1, PackedColumns0, LeftNeighbor , 3, 
                 MPI_COMM_WORLD, &ReceiveStatus);
			}
			
			// x flux.
			for (int i = 0; i < Nx+1; i++)
			{
				for (int j = 0; j < Ny; j++)
				{
					flux = dy * weno5(Phi[i-3][j],Phi[i-2][j],Phi[i-1][j],Phi[i][j],
                            Phi[i+1][j]);
					
					if (i == 0)
						Residual[i][j] -= flux;
					else if ( i == Nx)
						Residual[i-1][j] += flux;
					else
					{
						Residual[i][j] -= flux;
						Residual[i-1][j] += flux;
					}
				}
			}

			// y flux.
			for (int j = 0; j < Ny+1; j++)
			{
				for (int i = 0; i < Nx; i++)
				{
					flux = dx * weno5(Phi[i][j-3],Phi[i][j-2],Phi[i][j-1],Phi[i][j],
                            Phi[i][j+1]);
					
					if (j == 0)
						Residual[i][j] -= flux;
					else if(j == Ny)
						Residual[i][j-1] += flux;
					else
					{
						Residual[i][j] -= flux;
						Residual[i][j-1] += flux;
					}
				}
			}

			// Update in each domain.
			for (unsigned int i = 0; i < Nx; i++)
				for (unsigned int j = 0; j < Ny; j++)
					Phi[i][j] = ark[rk] * PhiOld[i][j] 
                      + (1.0-ark[rk]) * (Phi[i][j] - lambda*Residual[i][j]);
		}
		
		t += dt;
		iteration++;
		
		// Calculation of error
		double LocalError = 0.0, GlobalError;

		for (unsigned int i = 0; i < Nx; i++)
				for (unsigned int j = 0; j < Ny; j++)
        {
          double x = dx*(I*Nx+i-0.5) - t;
          double y = dy*(J*Ny+j-0.5) - t;
					LocalError += Absolute(Phi[i][j] - initial_condition(x,y));
        }

		MPI_Reduce(&LocalError,&GlobalError,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		
		if (rank == 0)
			cout << "Iteration = " << iteration << ", Current time = " << t 
           << ", Error = " << GlobalError << endl;
		
		// Save solution.
		if (iteration % SaveInterval == 0)
		{
			if (iteration % SaveInterval >= 10000)
			{
				if (rank == 0)
					fprintf(stderr,"Filename counter too large!");
		
				MPI_Abort(MPI_COMM_WORLD,1);
			}
			
			// Obtain the values at top-right corner
			MPI_Send(&Phi[0][0]  , 1, MPI_DOUBLE, TopRightNeighbor  , 4, 
               MPI_COMM_WORLD);
			MPI_Recv(&Phi[Nx][Ny], 1, MPI_DOUBLE, BottomLeftNeighbor, 4, 
               MPI_COMM_WORLD, &ReceiveStatus);
			
			char *s = new char [50];	
	
			sprintf(s,"Solution-%04d-%04d.tec",iteration/SaveInterval,rank);
			ofstream FileWrite(s, ios::out);
			FileWrite.flags( ios::dec | ios::scientific);
			FileWrite.precision(8);
	
			if (!FileWrite)
				ErrorMessage("Output file couldnot be opened.");
			
			FileWrite << "TITLE = \"2D_Advection\"" << endl;
			FileWrite << "VARIABLES = \"x\", \"y\", \"phi\", \"phiExact\""<< endl;
			FileWrite << "ZONE STRANDID=1, SOLUTIONTIME=" << t << ", I=" << Ny+1 
                << ", J=" << Nx+1 << ", DATAPACKING=POINT" << endl;
	
			for (int i = 0; i <= Nx; i++)
			{
				for (int j = 0; j <= Ny; j++)
				{
					double x = xmin + dx*(I*Nx+i-0.5);
					double y = ymin + dy*(J*Ny+j-0.5);
					FileWrite << x << "\t" << y << "\t" << Phi[i][j] << "\t" 
                    << initial_condition(x-t,y-t) << endl;
				}
			}
	
			FileWrite.close();
	
			delete [] s;
		}
	}
	
	Deallocate(Phi,Phicontiguous,Nx+XOffset+XOffsetEnd,Ny+YOffset+YOffsetEnd,
             XOffset,YOffset);
	Deallocate(PhiOld,Nx,Ny);
	Deallocate(Residual,Nx,Ny);
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	if (rank == 0)
	{
		tf = time(NULL);
		cout << "Starting time :" << ctime(&ti) << endl;
		cout << "Total time elapsed : " << (tf-ti) << " seconds" << endl;
	}

	// Finalize the MPI environment.
	MPI_Finalize();

	return 0;
}
