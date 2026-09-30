//------------------------------------------------------------------------------
// Saves solution into files sol000.h5, sol001.h5, etc.
// You can open them in VisIt.
//------------------------------------------------------------------------------
PetscErrorCode savesol(int *c, Vec ug)
{
   char           filename[32] = "sol";
   PetscViewer    viewer;
   PetscFunctionBeginUser;
   sprintf(filename, "sol%03d.h5", *c);
   PetscCall(PetscViewerHDF5Open(PETSC_COMM_WORLD, filename, FILE_MODE_WRITE, &viewer));
   PetscCall(VecView(ug, viewer));
   PetscCall(PetscViewerDestroy(&viewer));
   ++(*c);
   PetscFunctionReturn(PETSC_SUCCESS);
}
