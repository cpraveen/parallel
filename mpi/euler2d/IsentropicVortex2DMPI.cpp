//____________________________________________________________________________________
// This program (parallelized with MPI) solves 2D Euler equation on a square 
// doubly periodic domain.
//____________________________________________________________________________________
// Author : Pritam Giri
// Date   : 11.10.2016
// TIFR-CAM
//____________________________________________________________________________________
//____________________________________________________________________________________
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
const double Gamma = 1.4;
const double ark[] = {0.0, 3.0/4.0, 1.0/3.0};
const double brk[] = {1.0, 1.0/4.0, 2.0/3.0};

//____________________________________________________________________________________
// Absolute value
//____________________________________________________________________________________
template <class Type> 
Type Absolute(Type a)
{
	return (a > (Type)(0.0) ? a : -a);
}

//____________________________________________________________________________________
// Minimum value
//____________________________________________________________________________________
template <class Type> 
Type Minimum(Type a, Type b)
{
	return (a > b ? b : a);
}

//____________________________________________________________________________________
// Print an error message and exit the program
//____________________________________________________________________________________
void ErrorMessage(const char* message)
{
	cout << message << endl;
	exit(1);
}

//____________________________________________________________________________________
// Weno reconstruction
// Return left value for face between u0, up1
// Adapted from : https://github.com/cpraveen/cfdlab/blob/master/petsc/convect2d/convect.c
//____________________________________________________________________________________
double weno5(double um2, double um1, double u0, double up1, double up2)
{
	double eps = 1.0e-6;
	double Gamma1=1.0/10.0, Gamma2=3.0/5.0, Gamma3=3.0/10.0;
	double beta1, beta2, beta3;
	double u1, u2, u3;
	double w1, w2, w3;
	
	beta1 =   (13.0/12.0)*pow((um2 - 2.0*um1 + u0),2)
           + (1.0/4.0)*pow((um2 - 4.0*um1 + 3.0*u0),2);
	beta2 =   (13.0/12.0)*pow((um1 - 2.0*u0 + up1),2)
           + (1.0/4.0)*pow((um1 - up1),2);
	beta3 =   (13.0/12.0)*pow((u0 - 2.0*up1 + up2),2) 
           + (1.0/4.0)*pow((3.0*u0 - 4.0*up1 + up2),2);
	
	w1 = Gamma1 / pow(eps+beta1, 2);
	w2 = Gamma2 / pow(eps+beta2, 2);
	w3 = Gamma3 / pow(eps+beta3, 2);
	
	u1 = (1.0/3.0)*um2 - (7.0/6.0)*um1 + (11.0/6.0)*u0;
	u2 = -(1.0/6.0)*um1 + (5.0/6.0)*u0 + (1.0/3.0)*up1;
	u3 = (1.0/3.0)*u0 + (5.0/6.0)*up1 - (1.0/6.0)*up2;
	
	return (w1 * u1 + w2 * u2 + w3 * u3)/(w1 + w2 + w3);
}

//____________________________________________________________________________________
// Allocate a 1D array
//____________________________________________________________________________________
template <class Type> 
void Allocate(Type *&m, int d1)
{
	m = new (nothrow) Type [d1];
	if (m == 0)
		ErrorMessage("Error: Memory can not be allocated.");

	for (int i = 0; i < d1; i++)
		m[i] = (Type)(0.0);
}

//____________________________________________________________________________________
// Deallocate a 1D array
//____________________________________________________________________________________
template <class Type> 
void Deallocate(Type *m, int d1){
	delete [] m;
	m = NULL;
}

//____________________________________________________________________________________
// Allocate a 3D array
//____________________________________________________________________________________
template <class Type> 
void Allocate(Type ***&m, int d1, int d2, int d3)
{
	m = new (nothrow) Type** [d1];
	
	if (m == 0)
		ErrorMessage("Error: Memory can not be allocated.");
	
	for (int i = 0; i < d1; i++)
	{
      m[i] = new (nothrow) Type* [d2];
		
		if (m[i] == 0)
			ErrorMessage("Error: Memory can not be allocated.");

		for (int j = 0; j < d2; j++)
		{
			m[i][j] = new (nothrow) Type [d3];
			
			if (m[i][j] == 0)
			   ErrorMessage("Error: Memory can not be allocated.");
			
			for (int k = 0; k < d3; k++)
				m[i][j][k] = (Type)(0.0);
		}
	}
}

//____________________________________________________________________________________
// Deallocate a 3D array
//____________________________________________________________________________________
template <class Type> 
void Deallocate(Type ***m, int d1, int d2, int d3){
	for (int i = 0; i < d1; i++){
		for (int j = 0; j < d2; j++)
			delete [] m[i][j];
		delete [] m[i];
	}
	
	delete [] m;
	m = NULL;
}

//____________________________________________________________________________________
// Allocate a 3D array in contiguous memory locations. 
// Indexing starts with m[i][-xoffset][-yoffset].
//____________________________________________________________________________________
template <class Type> void Allocate(Type ***&m, Type *&mcontiguous, const int d1, const int d2, const int d3, const int xoffset = 0, const int yoffset = 0)
{
	Allocate(mcontiguous,d1*d2*d3);
	
	m = new (nothrow) Type** [d1];
	
	if (m == 0)
		ErrorMessage("Error: Memory can not be allocated.");
	
	for (int i = 0; i < d1; i++)
	{
      m[i] = new (nothrow) Type* [d2];
		
		if (m[i] == 0)
			ErrorMessage("Error: Memory can not be allocated.");

		for (int j = 0; j < d2; j++)
		{
			m[i][j] = &(mcontiguous[i*d2*d3+j*d3]);
						
			for (int k = 0; k < d3; k++)
				m[i][j][k] = (Type)(0.0);
		}
	}
	
	for (int i = 0; i < d1; i++)
		m[i] += xoffset;

	for (int i = 0; i < d1; i++)
		for (int j = -xoffset; j < (d2-xoffset); j++)
			m[i][j] += yoffset;
}

//____________________________________________________________________________________
// Deallocate a 3D array which has been allocated at contiguous memory locations. 
// Indexing starts with m[i][-xoffset][-yoffset].
//____________________________________________________________________________________
template <class Type> void Deallocate(Type ***m, Type *mcontiguous, const int d1, const int d2, const int d3, const int xoffset = 0, const int yoffset = 0)
{
	for (int i = 0; i < d1; i++)
		for (int j = -xoffset; j < (d2-xoffset); j++)
			m[i][j] -= yoffset;

	for (int i = 0; i < d1; i++)
		m[i] -= xoffset;

	for (int i = 0; i < d1; i++)
		delete [] m[i];
	
	delete [] m;
	m = NULL;

	Deallocate(mcontiguous,d1*d2*d3);
}

//____________________________________________________________________________________
// Periodic wrap-around neighbor search.
//____________________________________________________________________________________
void SearchNeighbor(int &TopNeighbor, int &BottomNeighbor, int &LeftNeighbor, int &RightNeighbor, int &TopRightNeighbor, int &BottomLeftNeighbor, const int I, const int J, const int XDim, const int YDim)
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

//____________________________________________________________________________________
// Conserved to primitive variables
//____________________________________________________________________________________
void ConservedToPrimitive(double *U, double *V)
{
  V[0] = U[0];
  V[1] = U[1]/U[0];
  V[2] = U[2]/U[0];
  V[3] = (U[3] - 0.5*V[0]*(V[1]*V[1] + V[2]*V[2]))*(Gamma-1.0);
}

//____________________________________________________________________________________
// Primitive to conserved variables
//____________________________________________________________________________________
void PrimitiveToConserved(double *U, double *V)
{
  U[0] = V[0];
  U[1] = V[0]*V[1];
  U[2] = V[0]*V[2];
  U[3] = 0.5*V[0]*(V[1]*V[1] + V[2]*V[2]) + V[3]/(Gamma-1.0);
}

//____________________________________________________________________________________
// Conserved to flux
//____________________________________________________________________________________
void ConservedToFlux(double *flux, double *U, double *V, const double nx, const double ny)
{
	ConservedToPrimitive(U,V);

	flux[0] = U[1]*nx + U[2]*ny;
	flux[1] = V[3]*nx + V[1]*flux[0];
	flux[2] = V[3]*ny + V[2]*flux[0];
	flux[3] = (U[3]+V[3])*(V[1]*nx + V[2]*ny);
}

//____________________________________________________________________________________
// Compute local time step
//____________________________________________________________________________________
double ComputeTimeStep(double ***U, double *Up, double *Vp, const double dx, const double dy, const int Nx, const int Ny)
{
	double dt = 100.0, u, a;
	
	for (unsigned int i = 0; i < Nx; i++)
	{
		for (unsigned int j = 0; j < Ny; j++)
		{
			for (unsigned int k = 0; k < 4; k++)
				Up[k] = U[k][i][j];

			ConservedToPrimitive(Up,Vp);

			u = (Vp[1]*Vp[1] + Vp[2]*Vp[2]);
			a = sqrt(Gamma*Vp[3]/Vp[0]);
			
			dt = Minimum(dt,Minimum(dx,dy)/(u+a));
		}
	}

	return dt;
}

//____________________________________________________________________________________
// Evaluate Lax Friedrich flux
//____________________________________________________________________________________
void LaxFriedrichFlux(double *flux, double *Uleft, double *Uright, double *Vleft, double *Vright, double *fluxleft, double *fluxright, const double nx, const double ny, const double lambda)
{
	ConservedToFlux(fluxleft,Uleft,Vleft,nx,ny);
	ConservedToFlux(fluxright,Uright,Vright,nx,ny);

	for (unsigned int i = 0; i < 4; i++)
		flux[i] = 0.5*(fluxleft[i]+fluxright[i]) + 0.5*(Uleft[i]-Uright[i])/lambda;
}

//____________________________________________________________________________________
// Evaluate LLF flux
//____________________________________________________________________________________
void LLFFlux(double *flux, double *Uleft, double *Uright, double *Vleft, double *Vright, double *fluxleft, double *fluxright, const double nx, const double ny)
{
	double lambda;
	
	ConservedToFlux(fluxleft,Uleft,Vleft,nx,ny);
	ConservedToFlux(fluxright,Uright,Vright,nx,ny);

	for (unsigned int i = 0; i < 4; i++)
	{
		flux[i] = 0.5*(fluxleft[i]+fluxright[i]);
		Uleft[i] = 0.5*(Uleft[i]+Uright[i]);
		Uright[i] -= Uleft[i];
	}

	ConservedToPrimitive(Uleft,Vleft);

	double u, a;

	u = Absolute(Vleft[1]*nx + Vleft[2]*ny);
	a = sqrt(Gamma*Vleft[3]/Vleft[0]);
	
	lambda = u + a;

	for (unsigned int i = 0; i < 4; i++)
		flux[i] -= lambda*Uright[i];
}

//____________________________________________________________________________________
//____________________________________________________________________________________
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
	cout.precision(4);
	ti = time(NULL);
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// Begin code
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
	//____________________________________________________________________________________
	// XDim    : Number of processors along x direction,
	// YDim    : Number of processors along y direction,
	// Nx      : Number of cells along x direction in each processor,
	// Ny      : Number of cells along y direction in each processor,
	// XOffset : Number of layers needed from neighbors for high order reconstruction in x direction,
	// YOffset : Number of layers needed from neighbors for high order reconstruction in y direction.
	
	// M       : Farfield Mach number of flow,
	// alpha   : Angle of advection (in degrees).
	
	// Note that global size of the domain (XDim*Nx) x (YDim*Ny)
	// Number of processors required to run this code is (XDim*YDim)
	
	const int XDim = 2, YDim = 2, Nx = 50, Ny = 50, XOffset = 3, YOffset = 3, XOffsetEnd = 3, YOffsetEnd = 3, variables = 4, SaveInterval = 100;
	const int I = rank/YDim, J = (rank - I*YDim);
	const double xmin = -5.0, xmax = 5.0, ymin = -5.0, ymax = 5.0, cfl = 0.4, Tf = 10.0;
	const double M = 0.5, alpha = 0.0, beta = 5.0;
	const double dx = (xmax-xmin)/(XDim*Nx), dy = (ymax-ymin)/(YDim*Ny);

	double *Ucontiguous, *Uleft, *Uright, *Vleft, *Vright, *flux, *fluxleft, *fluxright;
	double ***U, ***Uold, ***Residual;
	double t = 0.0, dt, lambda;
	
	int TopNeighbor, BottomNeighbor, LeftNeighbor, RightNeighbor, TopRightNeighbor, BottomLeftNeighbor;
	
	SearchNeighbor(TopNeighbor,BottomNeighbor,LeftNeighbor,RightNeighbor,TopRightNeighbor,BottomLeftNeighbor,I,J,XDim,YDim);
	
	if (size != (XDim*YDim))
	{
		if (rank == 0)
			fprintf(stderr,"%d processors are needed to run this code!",XDim*YDim);
		
		MPI_Abort(MPI_COMM_WORLD,1);
	}

	Allocate(U,Ucontiguous,variables,Nx+XOffset+XOffsetEnd,Ny+YOffset+YOffsetEnd,XOffset,YOffset);
	Allocate(Uold,variables,Nx,Ny);
	Allocate(Residual,variables,Nx,Ny);
	Allocate(Uleft,variables);
	Allocate(Uright,variables);
	Allocate(Vleft,variables);
	Allocate(Vright,variables);
	Allocate(flux,variables);
	Allocate(fluxleft,variables);
	Allocate(fluxright,variables);

	// Initialize.
	for (unsigned int i = 0; i < Nx; i++)
	{
		for (unsigned int j = 0; j < Ny; j++)
		{
			double x = xmin + dx*(I*Nx+i-0.5), y = ymin + dy*(J*Ny+j-0.5);

			Vleft[0] =  pow((1.0 - (Gamma-1.0)*beta*beta/(8.0*Gamma*PI*PI)*exp(1.0-x*x-y*y)),(1.0/(Gamma-1.0)));
			Vleft[1] =  M*cos(alpha*PI/180.0) - beta/(2.0*PI)*y*exp(0.5*(1.0-x*x-y*y));
			Vleft[2] =  M*sin(alpha*PI/180.0) + beta/(2.0*PI)*x*exp(0.5*(1.0-x*x-y*y));
			Vleft[3] =  pow(Vleft[0],Gamma);
			
			PrimitiveToConserved(Uleft,Vleft);

			for (unsigned int k = 0; k < variables; k++)
				U[k][i][j] = Uleft[k];
		}
	}

	// Time stepping.
	int iteration = 0;
	
	MPI_Datatype PackedColumns0, PackedColumns1;
	MPI_Status ReceiveStatus;

	MPI_Type_vector(Nx,YOffset,Ny+YOffset+YOffsetEnd,MPI_DOUBLE,&PackedColumns0);
	MPI_Type_commit(&PackedColumns0);

	MPI_Type_vector(Nx,YOffsetEnd,Ny+YOffset+YOffsetEnd,MPI_DOUBLE,&PackedColumns1);
	MPI_Type_commit(&PackedColumns1);
	
	while (t < Tf)
	{
		double Localdt = ComputeTimeStep(U,Uleft,Vleft,dx,dy,Nx,Ny);
		
		MPI_Allreduce(&Localdt,&dt,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
		dt *= cfl;

		if (t+dt > Tf)
			dt = Tf - t;
		
		lambda = dt/(dx*dy);

		for (unsigned int k = 0; k < variables; k++)
			for (unsigned int i = 0; i < Nx; i++)
				for (unsigned int j = 0; j < Ny; j++)
					Uold[k][i][j] = U[k][i][j];
		
		for (int rk = 0; rk < 3; rk++)
		{
			for (unsigned int k = 0; k < variables; k++)
				for (unsigned int i = 0; i < Nx; i++)
					for (unsigned int j = 0; j < Ny; j++)
						Residual[k][i][j] = 0.0;

			// Copy required portions from the neighbors.
			// message passing is needed only for multiple processors.
			if (TopNeighbor == rank)
			{
				for (int k = 0; k < 4; k++)
				{
					for (int i = 0; i < XOffsetEnd; i++)
						for (int j = 0; j < Ny; j++)
							U[k][Nx+i][j] = U[k][i][j];

					for (int i = 0; i < XOffset; i++)
						for (int j = 0; j < Ny; j++)
							U[k][i-XOffset][j] = U[k][Nx-XOffset+i][j];
				}
			}
			else
			{
				for (int k = 0; k < 4; k++)
				{
					MPI_Send(&U[k][0][-YOffset] , XOffsetEnd*(Ny+YOffset+YOffsetEnd), MPI_DOUBLE, TopNeighbor   , 0, MPI_COMM_WORLD);
					MPI_Recv(&U[k][Nx][-YOffset], XOffsetEnd*(Ny+YOffset+YOffsetEnd), MPI_DOUBLE, BottomNeighbor, 0, MPI_COMM_WORLD, &ReceiveStatus);
				
					MPI_Send(&U[k][Nx-XOffset][-YOffset], XOffset*(Ny+YOffset+YOffsetEnd), MPI_DOUBLE, BottomNeighbor, 1, MPI_COMM_WORLD);
					MPI_Recv(&U[k][-XOffset][-YOffset]  , XOffset*(Ny+YOffset+YOffsetEnd), MPI_DOUBLE, TopNeighbor   , 1, MPI_COMM_WORLD, &ReceiveStatus);
				}
			}

			if (LeftNeighbor == rank)
			{
				for (int k = 0; k < 4; k++)
				{
					for (int i = 0; i < Nx; i++)
						for (int j = 0; j < YOffsetEnd; j++)
							U[k][i][Ny+j] = U[k][i][j];

					for (int i = 0; i < Nx; i++)
						for (int j = 0; j < YOffset; j++)
							U[k][i][j-YOffset] = U[k][i][Ny-YOffset+j];
				}
			}
			else
			{
				for (int k = 0; k < 4; k++)
				{
					MPI_Send(&U[k][0][0] , 1, PackedColumns1, LeftNeighbor , 2, MPI_COMM_WORLD);
					MPI_Recv(&U[k][0][Ny], 1, PackedColumns1, RightNeighbor, 2, MPI_COMM_WORLD, &ReceiveStatus);

					MPI_Send(&U[k][0][Ny-YOffset], 1, PackedColumns0, RightNeighbor, 3, MPI_COMM_WORLD);
					MPI_Recv(&U[k][0][-YOffset]  , 1, PackedColumns0, LeftNeighbor , 3, MPI_COMM_WORLD, &ReceiveStatus);
				}
			}
			
			// x flux.
			for (int i = 0; i < Nx+1; i++)
			{
				for (int j = 0; j < Ny; j++)
				{
					for (int k = 0; k < 4; k++)
					{
						Uleft[k] = weno5(U[k][i-3][j],U[k][i-2][j],U[k][i-1][j],U[k][i][j],U[k][i+1][j]);
						Uright[k] = weno5(U[k][i+2][j],U[k][i+1][j],U[k][i][j],U[k][i-1][j],U[k][i-2][j]);
					}

					//LaxFriedrichFlux(flux,Uleft,Uright,Vleft,Vright,fluxleft,fluxright,1.0,0.0,dt/dx);
					LLFFlux(flux,Uleft,Uright,Vleft,Vright,fluxleft,fluxright,1.0,0.0);

					for (int k = 0; k < 4; k++)
					{
						if (i == 0)
							Residual[k][i][j] -= flux[k]*dy;
						else if ( i == Nx)
							Residual[k][i-1][j] += flux[k]*dy;
						else
						{
							Residual[k][i][j] -= flux[k]*dy;
							Residual[k][i-1][j] += flux[k]*dy;
						}
					}
				}
			}

			// y flux.
			for (int j = 0; j < Ny+1; j++)
			{
				for (int i = 0; i < Nx; i++)
				{
					for (int k = 0; k < 4; k++)
					{
						Uleft[k] = weno5(U[k][i][j-3],U[k][i][j-2],U[k][i][j-1],U[k][i][j],U[k][i][j+1]);
						Uright[k] = weno5(U[k][i][j+2],U[k][i][j+1],U[k][i][j],U[k][i][j-1],U[k][i][j-2]);
					}

					//LaxFriedrichFlux(flux,Uleft,Uright,Vleft,Vright,fluxleft,fluxright,0.0,1.0,dt/dy);
					LLFFlux(flux,Uleft,Uright,Vleft,Vright,fluxleft,fluxright,0.0,1.0);
					
					for (int k = 0; k < 4; k++)
					{
						if (j == 0)
							Residual[k][i][j] -= flux[k]*dx;
						else if(j == Ny)
							Residual[k][i][j-1] += flux[k]*dx;
						else
						{
							Residual[k][i][j] -= flux[k]*dx;
							Residual[k][i][j-1] += flux[k]*dx;
						}
					}
				}
			}

			// Update in each domain.
			for (unsigned int k = 0; k < 4; k++)
				for (unsigned int i = 0; i < Nx; i++)
					for (unsigned int j = 0; j < Ny; j++)
						U[k][i][j] = ark[rk] * Uold[k][i][j] + brk[rk] * (U[k][i][j] - lambda*Residual[k][i][j]);
		}
		
		t += dt;
		iteration++;
				
		if (rank == 0)
			cout << "Iteration = " << iteration << ", Current time = " << t << endl;
		
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
			for (unsigned int k = 0; k < variables; k++)
			{
				MPI_Send(&U[k][0][0]  , 1, MPI_DOUBLE, TopRightNeighbor  , 4, MPI_COMM_WORLD);
				MPI_Recv(&U[k][Nx][Ny], 1, MPI_DOUBLE, BottomLeftNeighbor, 4, MPI_COMM_WORLD, &ReceiveStatus);
			}


			char *s = new char [50];	
	
			sprintf(s,"Solution-%04d-%04d.tec",iteration/SaveInterval,rank);
			ofstream FileWrite(s, ios::out);
			FileWrite.flags( ios::dec | ios::scientific);
			FileWrite.precision(8);
	
			if (!FileWrite)
				ErrorMessage("Output file couldnot be opened.");
			
			FileWrite << "TITLE = \"U_t + U_x + U_y = 0\"" << endl;
			FileWrite << "VARIABLES = \"x\", \"y\", \"rho\", \"u\", \"v\", \"p\""<< endl;
			FileWrite << "ZONE STRANDID=1, SOLUTIONTIME=" << t << ", I=" << Ny+1 << ", J=" << Nx+1 << ", DATAPACKING=POINT" << endl;
			for (int i = 0; i <= Nx; i++)
			{
				for (int j = 0; j <= Ny; j++)
				{
					for (unsigned int k = 0; k < variables; k++)
						Uleft[k] = U[k][i][j];

					ConservedToPrimitive(Uleft,Vleft);
					
					FileWrite << xmin+dx*(I*Nx+i-0.5) << "\t" << ymin+dy*(J*Ny+j-0.5) << "\t" << Vleft[0] << "\t" << Vleft[1] << "\t" << Vleft[2] << "\t" << Vleft[3] << endl;
				}
			}
			
			FileWrite.close();
	
			delete [] s;
		}
	}
	
	Deallocate(U,Ucontiguous,variables,Nx+XOffset+XOffsetEnd,Ny+YOffset+YOffsetEnd,XOffset,YOffset);
	Deallocate(Uold,variables,Nx,Ny);
	Deallocate(Residual,variables,Nx,Ny);
	Deallocate(Uleft,variables);
	Deallocate(Uright,variables);
	Deallocate(Vleft,variables);
	Deallocate(Vright,variables);
	Deallocate(flux,variables);
	Deallocate(fluxleft,variables);
	Deallocate(fluxright,variables);

	//____________________________________________________________________________________

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// End code
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
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
