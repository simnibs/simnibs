/* Based on ex1 from PETSc
 *
 *
*/
#include <time.h>
#include <stdio.h>
#include <petscksp.h>

static PetscErrorCode _petsc_prepare_ksp(int, char**, PetscInt, PetscInt[], PetscInt[], PetscScalar[], FILE*, KSP*);
static PetscErrorCode _print_ksp_info(KSP, FILE*);
static PetscErrorCode _petsc_solve_with_ksp(KSP, PetscInt, PetscScalar[], FILE*, PetscScalar[]);
static PetscErrorCode _dealloc(KSP);

/* Creates the ksp
 * argc, args: Arguments to be passed to PetscInitialize
 * N: size of the system
 * row_indices: csr-style row indices for the system matrix
 * column_indices: csr-style column indices for the system matrix
 * matrix_values: Value of each matrix entry for the system matrix
 * stream: Where to redirect stderr and stdout
 * ksp: KSP object to be set-up
 * n_sims: number of simulations to run
 * rhs: right-hand-side vector 
 * solution: placeholder for the solution
*/
static PetscErrorCode _petsc_prepare_ksp(int argc,char **args, PetscInt N,
			      PetscInt row_indices[], PetscInt column_indices[],
			      PetscScalar matrix_values[], FILE *stream, KSP *ksp){

  Mat            A;            /* linear system matrix */
  PetscErrorCode ierr;

  /* Initialization */
  PETSC_STDOUT=stream;
  PETSC_STDERR=stream;
  ierr = MPI_Init(&argc, &args); // We need to call MPI_init if we want to call PetscInitialize many times
  ierr = PetscInitialize(&argc,&args,NULL, NULL);CHKERRQ(ierr);
  PETSC_STDERR=stream;
  PETSC_STDOUT=stream;
  /* Set-up matrix */
  ierr = MatCreateSeqAIJWithArrays(PETSC_COMM_WORLD, (PetscInt)N, (PetscInt)N, row_indices, column_indices, matrix_values, &A);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);


  /* Set-up KSP */
  ierr = KSPCreate(PETSC_COMM_WORLD,ksp);CHKERRQ(ierr);
  ierr = KSPSetOperators(*ksp,A,A);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(*ksp);CHKERRQ(ierr);
  //ierr = KSPMonitorSet(*ksp, monitor, NULL, 0);CHKERRQ(ierr);
  ierr = KSPSetUp(*ksp);CHKERRQ(ierr);

  ierr = MatDestroy(&A);CHKERRQ(ierr);
  return ierr;
}

static PetscErrorCode _print_ksp_info(KSP ksp, FILE *stream){
  PetscErrorCode        ierr;
  KSPType               ksp_type;
  PC                    pc;
  PCType                pc_type;

  PETSC_STDERR=stream;
  PETSC_STDOUT=stream;

  ierr = KSPGetType(ksp, &ksp_type);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp, &pc);CHKERRQ(ierr);
  ierr = PCGetType(pc, &pc_type);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Running PETSc with KSP: %s PC: %s\n", ksp_type, pc_type);CHKERRQ(ierr);
  return ierr;
}

/* Uses an already existing KSP object to solve
 * ksp: KSP object
 * rhs: right-hand-side vector 
 * stream: Where to redirect stderr and stdout
 * solution: placeholder for the solution
*/
static PetscErrorCode _petsc_solve_with_ksp(KSP ksp, PetscInt N, PetscScalar rhs[], FILE *stream, PetscScalar solution[]){
  Vec                   x, b;         /* Solution, RHS*/
  PetscErrorCode        ierr;
  KSPConvergedReason    reason;
  PetscReal             rnorm;
  PetscInt              niter;


  PETSC_STDOUT=stream;
  PETSC_STDERR=stream;
  /* Set-up RHS and solution vectors */
  ierr = VecCreateSeqWithArray(PETSC_COMM_WORLD, 1, N, rhs, &b);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(b);CHKERRQ(ierr);


  ierr = VecCreateSeqWithArray(PETSC_COMM_WORLD, 1, N, solution, &x);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(x);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(x);CHKERRQ(ierr);

  /* Solve */
  ierr = KSPSolve(ksp, b, x);CHKERRQ(ierr);

  /* Information */
  ierr = KSPGetResidualNorm(ksp, &rnorm);CHKERRQ(ierr);
  ierr = KSPGetIterationNumber(ksp, &niter);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Number of iterations: %i Residual Norm: %.2e\n", niter, rnorm);CHKERRQ(ierr);
  ierr = KSPGetConvergedReason(ksp, &reason);CHKERRQ(ierr);

  if (reason < 0){
	ierr = PetscPrintf(PETSC_COMM_WORLD, "KSP Diverged with reason: %i\n", reason);CHKERRQ(ierr);
	ierr = VecDestroy(&x);CHKERRQ(ierr);
	ierr = VecDestroy(&b);CHKERRQ(ierr);
	return 1;
  } else {
	ierr = PetscPrintf(PETSC_COMM_WORLD, "KSP Converged with reason: %i\n", reason);CHKERRQ(ierr);
  }


  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = VecDestroy(&b);CHKERRQ(ierr);

  return ierr;
}

/* Destroys the KSP object and finalizes PETSc
 * ksp: KSP object to be destroyed
*/
static PetscErrorCode _dealloc(KSP ksp){
  PetscErrorCode        ierr;

  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  ierr = PetscFinalize();CHKERRQ(ierr);CHKERRQ(ierr);
  return ierr;
}
