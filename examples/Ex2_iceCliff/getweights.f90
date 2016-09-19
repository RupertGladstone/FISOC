


RECURSIVE SUBROUTINE getWeights( Model,Solver,Timestep,TransientSimulation )
  USE DefUtils
 
  IMPLICIT NONE

  !------------------------------------------------------------------------------
  !    External variables
  !------------------------------------------------------------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t), TARGET :: Solver
  LOGICAL :: TransientSimulation
  REAL(KIND=dp) :: Timestep
  !------------------------------------------------------------------------------
  !    Local variables
  !------------------------------------------------------------------------------
  TYPE(ValueList_t), Pointer :: BC
  TYPE(Variable_t), POINTER :: Var
  TYPE(Element_t),POINTER :: Element
  INTEGER, POINTER :: VarPerm(:)
  INTEGER :: VarDOFs, i, k
  REAL(KIND=dp) :: x, xmax=1000.0, xmin=500.0
  REAL(KIND=dp), POINTER :: VarValues(:)
  LOGICAL :: GotIt

  CALL INFO("getWeights","Computing mesh velocity weighting distribution", Level=1)

  Var => Solver % Variable
  IF (ASSOCIATED(Var)) THEN
     VarPerm => Var % Perm
     VarDOFs =  Var % DOFs
     VarValues => Var % Values
  ELSE
     CALL FATAL('getWeights','No Variable associated')
  END IF
  k=0
  DO i = 1,Model % NumberOfNodes
     IF (VarPerm(i) > 0) THEN
        x = Solver % Mesh % Nodes % x(i)
        IF (x < xmin) THEN
           VarValues(VarPerm(i)) = 0.0
        ELSE
           VarValues(VarPerm(i)) = (x - xmin)/(xmax - xmin)
        END IF
        PRINT *, "weight:", x, VarValues(VarPerm(i))
     END IF
  END DO
END SUBROUTINE getWeights
