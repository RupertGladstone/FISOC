
MODULE FISOC_AM_Wrapper

  USE ESMF
  USE FISOC_utils_MOD

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: FISOC_AM_Wrapper_Init_Phase1,  FISOC_AM_Wrapper_Init_Phase2,  &
       FISOC_AM_Wrapper_Run, FISOC_AM_Wrapper_Finalize


CONTAINS


  !--------------------------------------------------------------------------------------
  ! This dummy wrapper aims to create the dummy grid and required variables 
  ! in the ESMF formats.  
  SUBROUTINE FISOC_AM_Wrapper_Init_Phase1(FISOC_config,vm,AM_ISM_ExpFB,AM_OM_ExpFB,AM_dummyGrid,rc)

    TYPE(ESMF_fieldBundle),INTENT(INOUT)  :: AM_ISM_ExpFB, AM_OM_ExpFB
    TYPE(ESMF_config),INTENT(INOUT)       :: FISOC_config
    TYPE(ESMF_grid),INTENT(OUT)           :: AM_dummyGrid
    TYPE(ESMF_VM),INTENT(IN)              :: vm
    INTEGER,INTENT(OUT),OPTIONAL          :: rc

    INTEGER                               :: localPet ! local persistent execution thread (1:1 relationship to process)
    CHARACTER(len=ESMF_MAXSTR)            :: label
    CHARACTER(len=ESMF_MAXSTR),ALLOCATABLE:: FISOC_AM_ReqVarList(:)
    LOGICAL                               :: verbose_coupling

    CALL ESMF_ConfigGetAttribute(FISOC_config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! extract a list of required atmos variables from the configuration object
    label = 'FISOC_AM_ReqVars:'
    CALL FISOC_getListFromConfig(FISOC_config, label, FISOC_AM_ReqVarList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       PRINT*,""
       PRINT*,"******************************************************************************"
       PRINT*,"**********    AM dummy wrapper.  Init phase 1 method.    *********************"
       PRINT*,"******************************************************************************"
       PRINT*,""
       PRINT*,"Here we need to get the AM grid information into the ESMF_grid type. "
       PRINT*,"We also need to create and initialise the required variables using the "
       PRINT*,"ESMF_field type and put them into an ESMF_fieldBundle type."
       PRINT*,""
    END IF

    CALL dummyCreateGrid(AM_dummyGrid,rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

! TODO:
! we'll need to separate out a required variables list from the AM for 
! both the OM and the ISM 
    CALL FISOC_populateFieldBundle(FISOC_AM_ReqVarList,AM_ISM_ExpFB,AM_dummyGrid,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL FISOC_populateFieldBundle(FISOC_AM_ReqVarList,AM_OM_ExpFB,AM_dummyGrid,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
  END SUBROUTINE FISOC_AM_Wrapper_Init_Phase1



  SUBROUTINE FISOC_AM_Wrapper_Init_Phase2(FISOC_config,vm,AM_ISM_ImpFB,AM_OM_ImpFB, &
       AM_ISM_ExpFB,AM_OM_ExpFB,rc)

    TYPE(ESMF_fieldBundle),INTENT(INOUT)  :: AM_ISM_ImpFB, AM_OM_ImpFB, AM_ISM_ExpFB, AM_OM_ExpFB
    TYPE(ESMF_config),INTENT(INOUT)       :: FISOC_config
    TYPE(ESMF_VM),INTENT(IN)              :: vm
    INTEGER,INTENT(OUT),OPTIONAL          :: rc

    INTEGER                               :: localPet
    LOGICAL                               :: verbose_coupling

    rc = ESMF_FAILURE

    CALL ESMF_VMGet(vm, localPet=localPet, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(FISOC_config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       PRINT*,""
       PRINT*,"******************************************************************************"
       PRINT*,"**********    AM dummy wrapper.  Init phase 2 method.     ********************"
       PRINT*,"******************************************************************************"
       PRINT*,""
       PRINT*,"Here we have access to the other component's initialised fields, just in "
       PRINT*,"case the AM needs to know about these in order to complete its initialisation."
       PRINT*,""
    END IF
    
    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_AM_Wrapper_Init_Phase2


  !--------------------------------------------------------------------------------------
  SUBROUTINE FISOC_AM_Wrapper_Run(FISOC_config,vm,AM_ISM_ImpFB,AM_OM_ImpFB,AM_ISM_ExpFB, & 
       AM_OM_ExpFB,rc)

    TYPE(ESMF_config),INTENT(INOUT)               :: FISOC_config
    TYPE(ESMF_fieldbundle),INTENT(INOUT),OPTIONAL :: AM_ISM_ImpFB,AM_OM_ImpFB,AM_ISM_ExpFB,AM_OM_ExpFB
    TYPE(ESMF_VM),INTENT(IN)                      :: vm
    INTEGER,INTENT(OUT),OPTIONAL                  :: rc

    INTEGER                      :: localPet
    LOGICAL                      :: verbose_coupling

    rc = ESMF_FAILURE

    CALL ESMF_VMGet(vm, localPet=localPet, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(FISOC_config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       PRINT*,""
       PRINT*,"******************************************************************************"
       PRINT*,"*************       AM dummy wrapper.  Run method.       ********************"
       PRINT*,"******************************************************************************"
       PRINT*,""
       PRINT*,""
    END IF

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_AM_Wrapper_Run


  !--------------------------------------------------------------------------------------
  SUBROUTINE FISOC_AM_Wrapper_Finalize(FISOC_config,vm,rc)

    TYPE(ESMF_config),INTENT(INOUT)    :: FISOC_config
    TYPE(ESMF_VM),INTENT(IN)           :: vm
    INTEGER,INTENT(OUT),OPTIONAL       :: rc

    INTEGER                            :: localPet
    LOGICAL                            :: verbose_coupling

    rc = ESMF_FAILURE

    CALL ESMF_VMGet(vm, localPet=localPet, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(FISOC_config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       PRINT*,""
       PRINT*,"******************************************************************************"
       PRINT*,"************    AM dummy wrapper.  Finalise method.     *********************"
       PRINT*,"******************************************************************************"
       PRINT*,""
       PRINT*,"FISOC has taken care of clearing up ESMF types.  Here we just need to call the "
       PRINT*,"AM finalise method."
    END IF

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_AM_Wrapper_Finalize


  !--------------------------------------------------------------------------------------
  ! copy pasted from the OM...
  SUBROUTINE dummyCreateGrid(AM_dummyGrid,rc)
    
    TYPE(ESMF_grid),INTENT(INOUT)        :: AM_dummyGrid    
    INTEGER,INTENT(OUT),OPTIONAL         :: rc

    REAL(ESMF_KIND_R8),POINTER :: coordY(:),coordX(:)
    REAL(ESMF_KIND_R8)         :: Lx, Ly, dx, dy
    INTEGER                    :: ii, jj, lbnd(1), ubnd(1), nx, ny

    NULLIFY (coordY,coordX)

    ! number of grid points
    nx = 200  
    ny = 50 

    ! domain size
    Lx = 300000.
    Ly = 30000.

    ! grid cell size
    dx = Lx/FLOAT(nx-1)
    dy = Ly/FLOAT(ny-1)

    ! this next grid creation section is more or less a copy from the ref documentation example:
    ! http://www.earthsystemmodeling.org/esmf_releases/public/last/ESMF_refdoc/node5.html#SECTION05083200000000000000
    !-------------------------------------------------------------------
    ! Create the Grid:  Allocate space for the Grid object, define the
    ! topology and distribution of the Grid, and specify that it 
    ! will have global indices.  Note that here aperiodic bounds are
    ! specified by the argument name. In this call the minIndex hasn't 
    ! been set, so it defaults to (1,1,...). The default is to 
    ! divide the index range as equally as possible among the DEs
    ! specified in regDecomp. This behavior can be changed by 
    ! specifying decompFlag. 
    !-------------------------------------------------------------------
    AM_dummyGrid=ESMF_GridCreateNoPeriDim(          &
         ! Define a regular distribution
         maxIndex=(/nx,ny/), & ! define index space
         !         regDecomp=(/2,3/),  & ! define how to divide among DEs
         coordSys=ESMF_COORDSYS_CART, &
         ! Specify mapping of coords dim to Grid dim
         coordDep1=(/1/), & ! 1st coord is 1D and depends on 1st Grid dim
         coordDep2=(/2/), & ! 2nd coord is 1D and depends on 2nd Grid dim
         indexflag=ESMF_INDEX_GLOBAL, &
         rc=rc)
    
    !-------------------------------------------------------------------
    ! Allocate coordinate storage and associate it with the center
    ! stagger location.  Since no coordinate values are specified in
    ! this call no coordinate values are set yet.
    !-------------------------------------------------------------------
    CALL ESMF_GridAddCoord(AM_dummyGrid,  & 
         staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
    
    !-------------------------------------------------------------------
    ! Get the pointer to the first coordinate array and the bounds
    ! of its global indices on the local DE.   
    !-------------------------------------------------------------------
    CALL ESMF_GridGetCoord(AM_dummyGrid, coordDim=1, localDE=0, &
         staggerloc=ESMF_STAGGERLOC_CENTER, &
         computationalLBound=lbnd, computationalUBound=ubnd, &
         farrayPtr=coordX, rc=rc)
    
    !-------------------------------------------------------------------
    ! Calculate and set coordinates in the first dimension.
    !-------------------------------------------------------------------
    DO ii=lbnd(1),ubnd(1)
       coordX(ii) = (ii-1)*dx
    END DO
    
    !-------------------------------------------------------------------
    ! Get the pointer to the second coordinate array and the bounds of
    ! its global indices on the local DE.
    !-------------------------------------------------------------------
    CALL ESMF_GridGetCoord(AM_dummyGrid, coordDim=2, localDE=0, &
         staggerloc=ESMF_STAGGERLOC_CENTER, &
         computationalLBound=lbnd, computationalUBound=ubnd, &
         farrayPtr=coordY, rc=rc)
    
    !-------------------------------------------------------------------
    ! Calculate and set coordinates in the second dimension 
    !-------------------------------------------------------------------
    DO jj=lbnd(1),ubnd(1)
       coordY(jj) = (jj-1)*dy
    END DO

    NULLIFY(coordX)
    NULLIFY(coordY)

  END SUBROUTINE dummyCreateGrid
    
END MODULE FISOC_AM_Wrapper
