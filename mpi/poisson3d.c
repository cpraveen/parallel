/*
 * 3D Poisson Equation
 * -𝚫(u) = 1 in [0,1] x [0,1] x [0,1]
 *  and u = 0 on the boundary
 *
 *  Translated from Fortran version
 *  C++ by Aadi Bhure
 *  Changed to C by Praveen C
 */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#define max(a,b) ((a) > (b)) ? (a) : (b)

void CopySendBuf(int iStart, int iEnd, 
                 int jStart, int jEnd, 
                 int kStart, int kEnd, 
                 double phi[iEnd-iStart][jEnd-jStart][kEnd-kStart], 
                 int disp, int dir,
                 double *fieldSend, int MaxBufLen);

void CopyRecvBuf(int iStart, int iEnd, 
                 int jStart, int jEnd, 
                 int kStart, int kEnd, 
                 double phi[iEnd-iStart][jEnd-jStart][kEnd-kStart], 
                 int disp, int dir,
                 double *fieldRecv, int MaxBufLen);

void Jacobi_sweep(int iStart, int iEnd, int jStart, int jEnd, int kStart, int kEnd,
                  double phi[2][iEnd-iStart][jEnd-jStart][kEnd-kStart], 
                  int t0, int t1, int udim[2][3], double h,
                  double *maxdelta);

// Keep reading characters until new line is encountered.
// We use this to skip comments in the input file.
void skip(FILE *fp)
{
   char c;
   while((c=getc(fp)) != '\n')
   {
   }
}

int main(int argc, char *argv[]) {
  // Mark whether boundaries are periodic or not
  int pbc_check[3];
  // Number of cells in each dimension
  int spat_dim[3];
  // Number of processes in each dimension
  int proc_dim[3];
  // dimensions of array belonging to just this rank
  int loca_dim[3];
  // Coordinates of this rank in the MPI Grid
  int mycoord[3];
  int totmsgsize[3];
  int myid, numprocs, ierr, iter, itermax, tag;
  int myid_grid, nump_grid, tmp, t0, t1;

  MPI_Comm GRID_COMM_WORLD;
  MPI_Request req;
  MPI_Status status;

  int iStart, jStart, kStart, iEnd, jEnd, kEnd, MaxBufLen;
  int source, dest, dir, disp;

  // A 2-d array for storing the start and end indices to
  // apply the sweep on.
  // Indexed by [dir][dimension]
  int udim[2][3];

  double eps, maxdelta, h;
  double *fieldSend, *fieldRecv;

  // Initialise MPI
  MPI_Init(&argc, &argv);

  // Set numprocs to be the number of ranks, which is set at run time.
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

  // Each rank has an id. Sets myid to this particular process' id.
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  if (myid == 0) {
    printf("Reading poisson3d.in\n");
    FILE *fp;
    fp = fopen("poisson3d.in", "r");
    fscanf(fp, "%d", &tmp); skip(fp);
    fscanf(fp, "%d%d%d", &proc_dim[0], &proc_dim[1], &proc_dim[2]); skip(fp);
    fscanf(fp, "%d", &itermax); skip(fp);
    fscanf(fp, "%lf", &eps);
    fclose(fp);

    // Total number of processes = product of the number of processes along each
    // dimension.
    if (numprocs != proc_dim[0] * proc_dim[1] * proc_dim[2]) {
      printf("Total procs cannot to factorized\n");
      printf("Total procs = %d\n", numprocs);
      printf("Proc grid   = %d %d %d\n", proc_dim[0], proc_dim[1], proc_dim[2]);
      ierr = MPI_Abort(MPI_COMM_WORLD, tmp);
    }

    // Number of cells along each dimension.
    spat_dim[0] = spat_dim[1] = spat_dim[2] = tmp;

    // Don't treat boundaries as periodic along any dimension
    pbc_check[0] = pbc_check[1] = pbc_check[2] = false;
  }

  // Bcast sends the root process' value to all other processes.
  ierr = MPI_Bcast(spat_dim, 3, MPI_INTEGER, 0, MPI_COMM_WORLD);
  ierr = MPI_Bcast(proc_dim, 3, MPI_INTEGER, 0, MPI_COMM_WORLD);
  ierr = MPI_Bcast(pbc_check, 3, MPI_LOGICAL, 0, MPI_COMM_WORLD);
  ierr = MPI_Bcast(&itermax, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
  ierr = MPI_Bcast(&eps, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // Partition ranks into a 3D topology, with proc_dim ranks along each dim
  ierr = MPI_Dims_create(numprocs, 3, proc_dim);

  // Assuming dx = dy = dz, so no need for separate hx, hy, hz
  // Minus 1 because the last cell will be at the end
  h = 1.0 / (spat_dim[0] - 1);

  if (myid == 0) {
    printf("Spatial Grid: %d %d %d\n", spat_dim[0], spat_dim[1], spat_dim[2]);
    printf("MPI     Grid: %d %d %d\n", proc_dim[0], proc_dim[1], proc_dim[2]);
    printf("Spatial h   : %e\n", h);
    printf("itermax     : %d\n", itermax);
    printf("eps         : %e\n", eps);
  }

  bool reorder = true;
  // create new communicator GRID_COMM_WORLD which has the topology info
  // reorder=false would restrict it to keep old and new ids the same,
  // but reorder=true allows us to give that up in return for possibly
  // better numbering.
  ierr = MPI_Cart_create(MPI_COMM_WORLD, 3, proc_dim, pbc_check, reorder,
                         &GRID_COMM_WORLD);

  if (GRID_COMM_WORLD == MPI_COMM_NULL) {
    if (myid == 0)
      printf("Failed to create GRID_COMM_WORLD\n");
    MPI_Abort(MPI_COMM_WORLD, tmp);
  }

  // get size and rank from this new communicator
  ierr = MPI_Comm_size(GRID_COMM_WORLD, &nump_grid);
  ierr = MPI_Comm_rank(GRID_COMM_WORLD, &myid_grid);

  // get grid coordinates for this rank
  ierr = MPI_Cart_coords(GRID_COMM_WORLD, myid_grid, 3, mycoord);

  // loca_dim is the grid size assigned to each current rank
  for (int i = 0; i < 3; ++i) {
    // Number of cells along i'th dimension divided by number of ranks along
    // that dimension This is truncated division, a/b = a//b * b + a%b
    loca_dim[i] = spat_dim[i] / proc_dim[i];

    // Compensate for the truncation.
    // Increase the size of (a%b) by one.
    if (mycoord[i] < spat_dim[i] % proc_dim[i])
      loca_dim[i] += 1;
  }

  // Solution variables
  // One layer of ghost points on all sides, so 2 extra indices
  // These are the first and last indices along each dir.
  // So iEnd is the number (count) of cells along each dimension, not last
  // index.
  iStart = 0;
  iEnd = loca_dim[0] + 2;
  jStart = 0;
  jEnd = loca_dim[1] + 2;
  kStart = 0;
  kEnd = loca_dim[2] + 2;

  // TODO: iStart,jStart,kStart are always 0, we could get rid of these variables

  // Initialise the solution array: phi[i][j][k][t] with zeroes
  // where t = 0 and 1 so that we can hold one older solution as well
  // Boundary conditions are Dirichlet, zero everywhere.
  //double (*phi)[iEnd-iStart][jEnd-jStart][kEnd-kStart] = calloc(2, sizeof(*phi));
  double phi[2][iEnd-iStart][jEnd-jStart][kEnd-kStart];

  MaxBufLen = 0;

  // Size of each face
  totmsgsize[2] = loca_dim[0] * loca_dim[1];
  MaxBufLen = max(MaxBufLen, totmsgsize[2]);

  totmsgsize[1] = loca_dim[0] * loca_dim[2];
  MaxBufLen = max(MaxBufLen, totmsgsize[1]);

  totmsgsize[0] = loca_dim[1] * loca_dim[2];
  MaxBufLen = max(MaxBufLen, totmsgsize[0]);

  // Buffers to send and receive data
  fieldSend = (double*)malloc(MaxBufLen * sizeof(double));
  fieldRecv = (double*)malloc(MaxBufLen * sizeof(double));

  disp = -1;

  // udim is different for each rank, and stores which points to update
  // and which ones to ignore because they will be filled in Dirichlet BC.
  for (dir = 0; dir < 3; ++dir) {
    // We want to transfer some data to the neighbouring ranks.
    // Cart_shift with the dir and disp tells us which rank to
    // send the data to and which rank to receive from, depending on the
    // topology.
    ierr = MPI_Cart_shift(GRID_COMM_WORLD, dir, disp, &source,
                          &dest);

    // dir = 0 means x-axis;
    // disp = -1 means moving to the left i.e. smaller indices.
    // So dest is left side, and source is right side

    // If boundary isn't periodic, then shifting out of range will be
    // communicated by NULL
    if (dest != MPI_PROC_NULL) {
      // We have a neighbour on the left.
      // So we can update every cell, starting from 1, skipping the halo cell
      udim[0][dir] = 1;
    } else
      // When no neighbour on the left.
      // This means this cell is on the boundary for this dir.
      // So skip over first two cells, because halo is now meaningless,
      // and the next cell is zero by Dirichlet BC, so no need to update.
      udim[0][dir] = 2;

    if (source != MPI_PROC_NULL)
      // Neighbour on the right, so index is the last one
      // loca_dim + 2 is the number of cells, so -1 will be the index of the
      // last cell so -2 will be the last non-halo cell.
      udim[1][dir] = loca_dim[dir] + 2 - 2;
    else
      // No neighbour, Dirichlet BC, so no need to update the very last
      // non-halo.
      udim[1][dir] = loca_dim[dir] + 2 - 2 - 1;
  }

  // Make initial guess
  for(int k=kStart; k<kEnd; ++k)
   for(int j=jStart; j<jEnd; ++j)
      for(int i=iStart; i<iEnd; ++i)
      {
         phi[0][i][j][k] = 0.0;
         phi[1][i][j][k] = 0.0;
      }

  // Begin iterations
  maxdelta = 2.0 * eps;

  // just indices for phi
  t0 = 0;
  t1 = 1;

  // tags are optional, meant so that ranks can easily differentiate between
  // messages
  tag = 0;
  iter = 0;

  while(iter < itermax && maxdelta > eps) {
    for (disp=-1; disp <= 1; disp+=2) {
      for (dir = 0; dir < 3; ++dir) {
        MPI_Cart_shift(GRID_COMM_WORLD, dir, disp, &source,
                       &dest);

        if (source != MPI_PROC_NULL)
          // We have a neighbour, so we receive into fieldRecv, the total
          // msgsize amount of data in that given dir. Irecv so that
          // everybody is ready to receive before everybody sends
          // TODO Replace this with MPI_SendRecv_replace?
          MPI_Irecv(fieldRecv, totmsgsize[dir], MPI_DOUBLE_PRECISION,
                    source, tag, GRID_COMM_WORLD, &req);

        if (dest != MPI_PROC_NULL) {
          // Copy into the send buffer the data to send, which is phi_old
          CopySendBuf(iStart, iEnd, jStart, jEnd, kStart, kEnd,
                      phi[t0], disp, dir, fieldSend, MaxBufLen);

          // then send that buffer to dest.
          MPI_Send(fieldSend, totmsgsize[dir], MPI_DOUBLE_PRECISION, dest,
                   tag, GRID_COMM_WORLD);
        }

        // Consequences of the Irecv
        if (source != MPI_PROC_NULL) {
          MPI_Wait(&req, &status);
          // Copy data from the receive buffer into phi_old
          CopyRecvBuf(iStart, iEnd, jStart, jEnd, kStart, kEnd,
                      phi[t0], disp, dir, fieldRecv, MaxBufLen);
        }
      }
    }

    // All communication done, we have the latest values. Now compute new values
    Jacobi_sweep(iStart, iEnd, jStart, jEnd, kStart, kEnd,
                 phi, t0, t1, udim, h, &maxdelta);

    // Find the largest delta amongst all ranks, and set it to the max_delta inO
    // every rank
    MPI_Allreduce(MPI_IN_PLACE, &maxdelta, 1, MPI_DOUBLE_PRECISION, MPI_MAX,
                  GRID_COMM_WORLD);

    ++iter;
    if (myid == 0) {
      printf("%12d %24.16e\n", iter, maxdelta);
    }
    // New becomes old
    tmp = t0; t0 = t1; t1 = tmp;
  }

  ierr = MPI_Finalize();
  return ierr;
}

// copies a face of phi (and not the layer of halo cells)
// into a linear array.
void CopySendBuf(int iStart, int iEnd, int jStart,
                 int jEnd, int kStart, int kEnd, 
                 double phi[iEnd-iStart][jEnd-jStart][kEnd-kStart],
                 int disp, int dir,
                 double *fieldSend, int MaxBufLen) {
  int i1, i2, j1, j2, k1, k2, c;

  if (dir < 0 || dir > 2) {
    printf("CSB: dir is wrong\n");
    exit(1);
  }
  if (disp != 1 && disp != -1) {
    printf("CSB: disp is wrong\n");
    exit(1);
  }

  // i2 and iEnd are both counts, so -1 is enough
  if (dir == 2) {
    // So we are dealing with a face parallel to the z-axis
    // skip the halo
    // Additional -1 because we are converting counts to indices.
    i1 = iStart + 1;
    i2 = iEnd - 2;
    j1 = jStart + 1;
    j2 = jEnd - 2;

    if (disp == -1)
      // upper face
      k1 = k2 = 1;
    else
      // lower face
      k1 = k2 = kEnd - 2;

  } else if (dir == 1) {
    i1 = iStart + 1;
    i2 = iEnd - 2;
    k1 = kStart + 1;
    k2 = kEnd - 2;

    if (disp == -1)
      j1 = j2 = 1;
    else
      j1 = j2 = jEnd - 1;
  } else if (dir == 0) {
    j1 = jStart + 1;
    j2 = jEnd - 2;
    k1 = kStart + 1;
    k2 = kEnd - 2;

    if (disp == -1)
      i1 = i2 = 1;
    else
      i1 = i2 = iEnd - 2;
  }

  c = 0;
  for (int k = k1; k <= k2; ++k)
    for (int j = j1; j <= j2; ++j)
      for (int i = i1; i <= i2; ++i) {
        fieldSend[c] = phi[i][j][k];
        c += 1;
      }
}

// Copy into the halo cells of phi the values in the 1D receive buffer
void CopyRecvBuf(int iStart, int iEnd, int jStart,
                 int jEnd, int kStart, int kEnd, 
                 double phi[iEnd-iStart][jEnd-jStart][kEnd-kStart],
                 int disp, int dir,
                 double *fieldRecv, int MaxBufLen) {
  int i1, i2, j1, j2, k1, k2, c;

  if (dir < 0 || dir > 2) {
    printf("CRB: dir is wrong\n");
    exit(1);
  }
  if (disp != 1 && disp != -1) {
    printf("CRB: disp is wrong\n");
    exit(1);
  }

  if (dir == 2) {
    // Same logic as in Send
    i1 = iStart + 1;
    i2 = iEnd - 2;
    j1 = jStart + 1;
    j2 = jEnd - 2;

    // We are receiving into the halo cells on the correct face
    if (disp == 1)
      // receiving from above
      k1 = k2 = 0;
    else
      // receiving from below
      k1 = k2 = kEnd - 1;

  } else if (dir == 1) {
    i1 = iStart + 1;
    i2 = iEnd - 2;
    k1 = kStart + 1;
    k2 = kEnd - 2;

    if (disp == 1)
      j1 = j2 = 0;
    else
      j1 = j2 = jEnd - 1;

  } else if (dir == 0) {
    j1 = jStart + 1;
    j2 = jEnd - 2;
    k1 = kStart + 1;
    k2 = kEnd - 2;

    if (disp == 1)
      i1 = i2 = 0;
    else
      i1 = i2 = iEnd - 1;
  }

  c = 0;
  for (int k = k1; k <= k2; ++k)
    for (int j = j1; j <= j2; ++j)
      for (int i = i1; i <= i2; ++i) {
        phi[i][j][k] = fieldRecv[c];
        c += 1;
      }
}

void Jacobi_sweep(int iStart, int iEnd, int jStart, int jEnd, int kStart, int kEnd,
                  double phi[2][iEnd-iStart][jEnd-jStart][kEnd-kStart], 
                  int t0, int t1, int udim[2][3], double h,
                  double *maxdelta) {
  double rhs = 1.0;
  double one_over_six = 1.0 / 6.0;

  *maxdelta = 0.0;

  for (int k = udim[0][2]; k <= udim[1][2]; ++k)
    for (int j = udim[0][1]; j <= udim[1][1]; ++j)
      for (int i = udim[0][0]; i <= udim[1][0]; ++i) {
        phi[t1][i][j][k] =
            (phi[t0][i - 1][j][k] + phi[t0][i + 1][j][k] +
             phi[t0][i][j - 1][k] + phi[t0][i][j + 1][k] +
             phi[t0][i][j][k - 1] + phi[t0][i][j][k + 1] + h * h * rhs) *
            one_over_six;
        *maxdelta = max(*maxdelta, fabs(phi[t1][i][j][k] - phi[t0][i][j][k]));
      }
}
