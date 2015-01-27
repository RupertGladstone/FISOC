
MODULE FISOC_ISM_Wrapper

  USE ESMF
  USE ElmerSolver_mod
  USE MainUtils

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: FISOC_ISM_Wrapper_Init,  FISOC_ISM_Wrapper_Run, FISOC_ISM_Wrapper_Finalize

  ! Note that CurrentModel is shared through the Types module

CONTAINS

  ! This initialisation wrapper aims to convert the Elmer mesh and required variables 
  ! to the ESMF formats.  It also performs simple sanity/consistency checks.
  SUBROUTINE FISOC_ISM_Wrapper_Init(ESMF_ElmerFieldBundle,ESMF_ElmerMesh,ESMF_ElmerConfig)

    TYPE(ESMF_fieldBundle),INTENT(INOUT) :: ESMF_ElmerFieldBundle
    TYPE(ESMF_mesh),INTENT(INOUT)        :: ESMF_ElmerMesh
    TYPE(ESMF_config),INTENT(IN)         :: ESMF_ElmerConfig

    TYPE(Mesh_t)                         :: Elmer_Mesh

! get elmer input params, especially time stepping information

! get elmer variables list, recieve esmf required variables list (from esmf_config?), 
! check the elmer contains those variables
! (use a look up table? or lookup table embedded in esmf config object?).

!convert required variables to esmf format

!store required vars with mesh in export state 

!sanity checks: timestep size, number of time intervals to run for

    CALL ElmerSolver_init(Elmer_Mesh) ! Intended to return the mesh prior to extrusion (need to add check for this mesh)

    CALL Elmer2ESMF_mesh(Elmer_mesh,ESMF_ElmerMesh)

!note: variables tobe input (basal melt rate) should be defined (perhaps as exported vars) in the 
! sif.  These also to be checked for their presence against a list of required vars from ESMF

  END SUBROUTINE FISOC_ISM_Wrapper_Init
  
  SUBROUTINE FISOC_ISM_Wrapper_Run()

! get hold of the elmer variables for receiving inputs, and convert them here from esmf to elmer type.
    CALL ElmerSolver_run()
! get hold of list of required variables from Elmer and convert them here from elmer to esmf type.

  END SUBROUTINE FISOC_ISM_Wrapper_Run

  SUBROUTINE FISOC_ISM_Wrapper_Finalize()

    CALL ElmerSolver_finalize()

  END SUBROUTINE FISOC_ISM_Wrapper_Finalize

  SUBROUTINE Elmer2ESMF_mesh(Elmer_mesh,ESMF_ElmerMesh)
    TYPE(ESMF_mesh),INTENT(INOUT)        :: ESMF_ElmerMesh
    TYPE(Mesh_t)                         :: Elmer_Mesh

    print*,"need to convert the mesh here"

  END SUBROUTINE Elmer2ESMF_mesh

END MODULE FISOC_ISM_Wrapper
