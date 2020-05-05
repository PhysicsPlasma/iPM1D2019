#include <petscksp.h>

int petscsolver_(double *CoeA, double *CoeB, double *CoeC, double *Source, double *Solve, int *Ns, int *ErrorHandler)
{
  Vec            x,b;              /* approx solution, RHS */
  Mat            A;                /* linear system matrix */
  KSP            ksp;              /* linear solver context */
  PC             pc;               /* preconditioner context */
  PetscErrorCode ierr;
  
  PetscInt       i,n = *Ns,its,rstart,rend,nlocal,column[3],*idx;
  PetscScalar    value[3];//, *array;
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         Compute the matrix and right-hand-side vector that define
         the linear system, Ax = b.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = VecCreate(PETSC_COMM_WORLD,&x); CHKERRQ(ierr);
  ierr = VecSetSizes(x,PETSC_DECIDE,n); CHKERRQ(ierr);
  ierr = VecSetFromOptions(x); CHKERRQ(ierr);
  ierr = VecDuplicate(x,&b); CHKERRQ(ierr);

  /* Identify the starting and ending mesh points on each
     processor for the interior part of the mesh. We let PETSc decide
     above. */

  ierr = VecGetOwnershipRange(x,&rstart,&rend); CHKERRQ(ierr);
  ierr = VecGetLocalSize(x,&nlocal); CHKERRQ(ierr);

  idx = (PetscInt*)malloc(sizeof(PetscInt) * nlocal);
  for(i=0; i<nlocal; i++)
    idx[i] = i + rstart;

  ierr = VecSetValues(b,nlocal,idx,&Source[rstart],INSERT_VALUES); CHKERRQ(ierr);

  ierr = VecAssemblyBegin(b); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(b); CHKERRQ(ierr);

  /*
     Create matrix.  When using MatCreate(), the matrix format can
     be specified at runtime.

     Performance tuning note:  For problems of substantial size,
     preallocation of matrix memory is crucial for attaining good
     performance. See the matrix chapter of the users manual for details.

     We pass in nlocal as the "local" size of the matrix to force it
     to have the same parallel layout as the vector created above.
  */
  ierr = MatCreate(PETSC_COMM_WORLD,&A); CHKERRQ(ierr);
  ierr = MatSetSizes(A,nlocal,nlocal,n,n); CHKERRQ(ierr);
  ierr = MatSetFromOptions(A); CHKERRQ(ierr);

//  MatSetUp(A);
  // preallocate matrix memory to improve performance
  // the type of matrix is AIJ
  ierr = MatSeqAIJSetPreallocation(A, 3, NULL); CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation(A, 3, NULL, 3, NULL); CHKERRQ(ierr);

  /*
     Assemble matrix.

     The linear system is distributed across the processors by
     chunks of contiguous rows, which correspond to contiguous
     sections of the mesh on which the problem is discretized.
     For matrix assembly, each processor contributes entries for
     the part that it owns locally.
  */

  
  if (!rstart) {
    rstart = 1; i = 0;
    column[0] = 0; column[1] = 1;
    value[0] = CoeB[0]; value[1] = CoeC[0];
    ierr = MatSetValues(A,1,&i,2,column,value,INSERT_VALUES); CHKERRQ(ierr);
  }

  if (rend == n) {
    rend = n-1; i = n-1;
    column[0] = n-2; column[1] = n-1;
    value[0] = CoeA[rend]; value[1] = CoeB[rend];
    ierr = MatSetValues(A,1,&i,2,column,value,INSERT_VALUES); CHKERRQ(ierr);
  }

  /* Set entries corresponding to the mesh interior */
  for (i=rstart; i<rend; ++i) {
    column[0] = i-1; column[1] = i; column[2] = i+1;
    value[0] = CoeA[i]; value[1] = CoeB[i]; value[2] = CoeC[i];
    ierr = MatSetValues(A,1,&i,3,column,value,INSERT_VALUES); CHKERRQ(ierr);
  }

  if(rstart == 1) --rstart;
  if(rend == n - 1) ++rend;
  //flawed, adjustment needed for further application
  //printf("n = %d, rstart = %d, rend = %d, nlocal = %d\n", n, rstart, rend, nlocal);  

  /* Assemble the matrix */
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  // ierr = MatView(A,PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);

  /*
     Set right-hand-side vector.
  */
  //ierr = VecSet(x,0.5); CHKERRQ(ierr);
  //for (i=0; i<n; i++)
  //  ierr = VecSetValues(b,1.0,&i,&Source[i],INSERT_VALUES); CHKERRQ(ierr);
  //ierr = VecSet(x,1.0); CHKERRQ(ierr);
  //ierr = VecView(b, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the linear solver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Create linear solver context
  */
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp); CHKERRQ(ierr);

  /*
     Set operators. Here the matrix that defines the linear system
     also serves as the preconditioning matrix.
  */
  ierr = KSPSetOperators(ksp,A,A); CHKERRQ(ierr);
  //ierr = KSPSetType(ksp,KSPBICG); CHKERRQ(ierr);
  //ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE); CHKERRQ(ierr);

  /*
     Set linear solver defaults for this problem (optional).
     - By extracting the KSP and PC contexts from the KSP context,
       we can then directly call any KSP and PC routines to set
       various options.
     - The following four statements are optional; all of these
       parameters could alternatively be specified at runtime via
       KSPSetFromOptions();
  */
/*
  ierr = KSPGetPC(ksp,&pc); CHKERRQ(ierr);
  ierr = PCSetType(pc,PCILU); CHKERRQ(ierr);
*/
  ierr = KSPSetTolerances(ksp,1.e-7,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT); CHKERRQ(ierr);
  /*
    Set runtime options, e.g.,
        -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
    These options will override those specified above as long as
    KSPSetFromOptions() is called _after_ any other customization
    routines.
  */
  ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the linear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Solve linear system
  */
  ierr = KSPSolve(ksp,b,x); CHKERRQ(ierr);

  /*
     View solver info; we could instead use the option -ksp_view to
     print this info to the screen at the conclusion of KSPSolve().
  */
/*
  printf("-----------------\n");
  printf("view the KSP info\n");
  printf("-----------------\n");
  KSPView(ksp,PETSC_VIEWER_STDOUT_SELF);
*/
  //ierr = MatView(A, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  //ierr = VecView(b, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  /*
     Pass the value of x to Solve
  */

/*  for (i=rstart; i<rend; i++)
    ierr = VecGetValues(x,1,&i,&Solve[i]); CHKERRQ(ierr);
*/
  ierr = VecGetValues(x,nlocal,idx,&Solve[rstart]); CHKERRQ(ierr);

  free(idx);

/*
  array = (PetscScalar*)malloc(sizeof(PetscScalar) * nlocal);
  ierr = VecGetArray(x,&array); CHKERRQ(ierr);
  for (i=0; i<nlocal; i++)
    Solve[rstart + i] = array[i];
  ierr = VecRestoreArray(x,&array); CHKERRQ(ierr);
  free(array);
*/
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Check solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Check the error
  */

/*
  VecAXPY(x,-1.0,u);
  VecNorm(x,NORM_2,&norm);
  if (norm > tol) {
    PetscPrintf(PETSC_COMM_WORLD,"Norm of error %g, Iterations %D\n",(double)norm,its);
  }
*/
/*
  printf("-----------------\n");
  ierr = KSPGetIterationNumber(ksp,&its); CHKERRQ(ierr);
  printf("iterations = %d\n", its);
*/
  /*
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
  */
  ierr = VecDestroy(&x); CHKERRQ(ierr);
  ierr = VecDestroy(&b); CHKERRQ(ierr);
  ierr = MatDestroy(&A); CHKERRQ(ierr);
  ierr = KSPDestroy(&ksp); CHKERRQ(ierr);

  /*
     Always call PetscFinalize() before exiting a program.  This routine
       - finalizes the PETSc libraries as well as MPI
       - provides summary and diagnostic information if certain runtime
         options are chosen (e.g., -log_view).
  */
  //ierr = PetscFinalize();
  *ErrorHandler = ierr;
  return ierr;
}

