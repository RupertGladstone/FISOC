
MODULE FISOC_OM_Wrapper

  USE ESMF
  USE FISOC_utils_MOD

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: FISOC_OM_Wrapper_Init_Phase1,  FISOC_OM_Wrapper_Init_Phase2,  &
       FISOC_OM_Wrapper_Run, FISOC_OM_Wrapper_Finalize, OM_HandleCavity

CONTAINS

  !--------------------------------------------------------------------------------------
  ! This dummy wrapper aims to create the dummy grid and required variables 
  ! in the ESMF formats.  
  SUBROUTINE FISOC_OM_Wrapper_Init_Phase1(OM_ExpFB,OM_dummyGrid,FISOC_config,vm,rc)

    TYPE(ESMF_fieldBundle),INTENT(INOUT)  :: OM_ExpFB
    TYPE(ESMF_config),INTENT(INOUT)       :: FISOC_config
    TYPE(ESMF_grid),INTENT(OUT)           :: OM_dummyGrid
    TYPE(ESMF_VM),INTENT(IN)              :: vm
    INTEGER,INTENT(OUT),OPTIONAL          :: rc

    INTEGER                    :: mpic ! mpi comm, duplicate from the OM VM
    INTEGER                    :: localPet ! local persistent execution thread (1:1 relationship to process)
    CHARACTER(len=ESMF_MAXSTR)            :: label
    CHARACTER(len=ESMF_MAXSTR),ALLOCATABLE:: FISOC_OM_ReqVarList(:)
    LOGICAL                               :: verbose_coupling

    CALL ESMF_ConfigGetAttribute(FISOC_config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! extract a list of required ocean variables from the configuration object
    label = 'FISOC_OM_ReqVars:' ! the FISOC names for the vars
    CALL FISOC_getListFromConfig(FISOC_config, label, FISOC_OM_ReqVarList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       PRINT*,""
       PRINT*,"******************************************************************************"
       PRINT*,"**********    OM dummy wrapper.  Init phase 1 method.    *********************"
       PRINT*,"******************************************************************************"
       PRINT*,""
       PRINT*,"Here we need to get the OM grid information into the ESMF_grid type. "
       PRINT*,"We also need to create and initialise the required variables using the "
       PRINT*,"ESMF_field type and put them into an ESMF_fieldBundle type."
       PRINT*,""
    END IF

    CALL dummyCreateGrid(OM_dummyGrid,rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL FISOC_populateFieldBundle(FISOC_OM_ReqVarList,OM_ExpFB,OM_dummyGrid,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
  END SUBROUTINE FISOC_OM_Wrapper_Init_Phase1
  
  !--------------------------------------------------------------------------------------
  ! This dummy wrapper aims to create the dummy grid and required variables 
  ! in the ESMF formats.  
  SUBROUTINE FISOC_OM_Wrapper_Init_Phase2(OM_ImpFB,OM_ExpFB,FISOC_config,vm,rc)

    TYPE(ESMF_config),INTENT(INOUT)       :: FISOC_config
    TYPE(ESMF_fieldBundle),INTENT(INOUT)  :: OM_ImpFB, OM_ExpFB
    INTEGER,INTENT(OUT),OPTIONAL          :: rc
    TYPE(ESMF_VM),INTENT(IN)              :: vm

    LOGICAL                      :: verbose_coupling
    INTEGER                      :: localPet
    TYPE(ESMF_field)             :: ISM_temperature_l0
    REAL(ESMF_KIND_R8),POINTER   :: ISM_temperature_l0_ptr(:,:)
    CHARACTER(len=ESMF_MAXSTR)   :: nameList(10)

    rc = ESMF_FAILURE

    NULLIFY(ISM_temperature_l0_ptr)

    CALL ESMF_ConfigGetAttribute(FISOC_config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       PRINT*,""
       PRINT*,"******************************************************************************"
       PRINT*,"**********    OM dummy wrapper.  Init phase 2 method.    *********************"
       PRINT*,"******************************************************************************"
       PRINT*,""
       PRINT*,"Here we have access to the initialised ISM fields, just in case the OM needs "
       PRINT*,"to know about these in order to complete its initialisation."
       PRINT*,""
    END IF


    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       
       CALL ESMF_FieldBundleGet(OM_ImpFB, fieldName="ISM_temperature_l0", field=ISM_temperature_l0, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       
       CALL ESMF_FieldGet(field=ISM_temperature_l0, localDe=0, farrayPtr=ISM_temperature_l0_ptr, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       PRINT*,"Temperature field size in x direction is ",SIZE(ISM_temperature_l0_ptr(:,1))
       PRINT*,"Temperature field size in y direction is ",SIZE(ISM_temperature_l0_ptr(1,:))
       PRINT*,"Show a few rows of data... we originally set the regridding to fill with zeros where source (ISM)"
       PRINT*,"grid doesn't cover destination (OM) grid."
       PRINT*,"Row 1 data:  ",ISM_temperature_l0_ptr(1,:)
       PRINT*,"Row 2 data:  ",ISM_temperature_l0_ptr(2,:)
       PRINT*,"Row 3 data:  ",ISM_temperature_l0_ptr(3,:)
       PRINT*,"Row 11 data: ",ISM_temperature_l0_ptr(11,:)
       PRINT*,""
       
    END IF
    
  END SUBROUTINE FISOC_OM_Wrapper_Init_Phase2
  
  
  !--------------------------------------------------------------------------------------
  SUBROUTINE FISOC_OM_Wrapper_Run(FISOC_config,vm,OM_ExpFB,OM_ImpFB,rc_local)
    
    TYPE(ESMF_config),INTENT(INOUT)                :: FISOC_config
    TYPE(ESMF_fieldBundle),INTENT(INOUT),OPTIONAL  :: OM_ExpFB, OM_ImpFB 
    TYPE(ESMF_VM),INTENT(IN)                       :: vm
    INTEGER,INTENT(OUT),OPTIONAL                   :: rc_local

    INTEGER                    :: localPet, rc
    LOGICAL                    :: verbose_coupling
    TYPE(ESMF_field)           :: ISM_dTdz_l0,ISM_z_l0, OM_dBdt_l0
    REAL(ESMF_KIND_R8),POINTER :: ISM_dTdz_l0_ptr(:,:), ISM_z_l0_ptr(:,:), OM_dBdt_l0_ptr(:,:)

    rc_local = ESMF_FAILURE

    CALL ESMF_ConfigGetAttribute(FISOC_config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       PRINT*,""
       PRINT*,"******************************************************************************"
       PRINT*,"************       OM dummy wrapper.  Run method.       **********************"
       PRINT*,"******************************************************************************"
       PRINT*,""

       IF ((PRESENT(OM_ExpFB)).AND.(.NOT.(PRESENT(OM_ImpFB)))) THEN
          PRINT*,"We have no new inputs for the OM from the ISM.  We need to call the OM "
          PRINT*,"and record its output in the OM export field bundle."
       END IF

       IF ((PRESENT(OM_ExpFB)).AND.(PRESENT(OM_ImpFB))) THEN
          PRINT*,"We have new inputs for the OM from the ISM in the OM import field bundle. "
          PRINT*,"We need to send these inputs to the OM, run one timestep, and record the OM "
          PRINT*,"outputs. "
       END IF
       
       IF ((.NOT.(PRESENT(OM_ExpFB))).AND.(PRESENT(OM_ImpFB))) THEN
          PRINT*,"We have new inputs for the OM from the ISM in the OM import field bundle. "
          PRINT*,"We need to send these inputs to the OM and run one timestep. We do not "
          PRINT*,"need to collect OM outputs.  Just run the OM one timestep."
       END IF

       IF ((.NOT.(PRESENT(OM_ExpFB))).AND.(.NOT.(PRESENT(OM_ImpFB)))) THEN
          PRINT*,"We have no new inputs for the OM from the ISM, and we do not need to "
          PRINT*,"collect OM outputs.  Just run the OM one timestep."
       END IF

    END IF

    ! Lets get pointers to the depth of the ice base and the temperature gradient.  These we 
    ! get fromthe OM import state, which contains the ISM export fields on the ocean grid.
    IF (PRESENT(OM_ImpFB)) THEN

       CALL ESMF_FieldBundleGet(OM_ImpFB, fieldName="ISM_z_l0", field=ISM_z_l0, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)       
       CALL ESMF_FieldGet(field=ISM_z_l0, localDe=0, farrayPtr=ISM_z_l0_ptr, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       
       CALL ESMF_FieldBundleGet(OM_ImpFB, fieldName="ISM_dTdz_l0", field=ISM_dTdz_l0, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)       
       CALL ESMF_FieldGet(field=ISM_dTdz_l0, localDe=0, farrayPtr=ISM_dTdz_l0_ptr, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    IF (PRESENT(OM_ExpFB)) THEN
       ! Lets get a pointer to the basal melt rate.  This we get from the OM export field bundle, which 
       ! contains the OM variables to be exported to the ISM.
       CALL ESMF_FieldBundleGet(OM_ExpFB, fieldName="OM_dBdt_l0", field=OM_dBdt_l0, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)       
       CALL ESMF_FieldGet(field=OM_dBdt_l0, localDe=0, farrayPtr=OM_dBdt_l0_ptr, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF
    
    rc_local = ESMF_SUCCESS
    
  END SUBROUTINE FISOC_OM_Wrapper_Run


  !--------------------------------------------------------------------------------------
  SUBROUTINE FISOC_OM_Wrapper_Finalize(FISOC_config,localPet,rc)

    TYPE(ESMF_config),INTENT(INOUT)    :: FISOC_config
    INTEGER,INTENT(IN)                 :: localPet
    INTEGER,INTENT(OUT),OPTIONAL       :: rc

    LOGICAL                            :: verbose_coupling

    rc = ESMF_FAILURE

    CALL ESMF_ConfigGetAttribute(FISOC_config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       PRINT*,""
       PRINT*,"******************************************************************************"
       PRINT*,"************    OM dummy wrapper.  Finalise method.     **********************"
       PRINT*,"******************************************************************************"
       PRINT*,""
       PRINT*,"FISOC has taken care of clearing up ESMF types.  Here we just need to call the "
       PRINT*,"OM finalise method."
    END IF

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_OM_Wrapper_Finalize


  !--------------------------------------------------------------------------------------
  SUBROUTINE dummyCreateGrid(OM_dummyGrid,rc)
    
    TYPE(ESMF_grid),INTENT(INOUT)        :: OM_dummyGrid    
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
    OM_dummyGrid=ESMF_GridCreateNoPeriDim(          &
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
    CALL ESMF_GridAddCoord(OM_dummyGrid,  & 
         staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
    
    !-------------------------------------------------------------------
    ! Get the pointer to the first coordinate array and the bounds
    ! of its global indices on the local DE.   
    !-------------------------------------------------------------------
    CALL ESMF_GridGetCoord(OM_dummyGrid, coordDim=1, localDE=0, &
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
    CALL ESMF_GridGetCoord(OM_dummyGrid, coordDim=2, localDE=0, &
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

  SUBROUTINE OM_HandleCavity(FISOC_config, FISOC_clock, OM_ImpFB, OM_ExpFB, localPet, rc)

    TYPE(ESMF_config),INTENT(INOUT)          :: FISOC_config
    TYPE(ESMF_fieldBundle),INTENT(INOUT)     :: OM_ImpFB, OM_ExpFB
    TYPE(ESMF_Clock),INTENT(IN)              :: FISOC_clock
    INTEGER,INTENT(IN)                       :: localPet
    INTEGER,INTENT(OUT),OPTIONAL             :: rc

  END SUBROUTINE OM_HandleCavity

END MODULE FISOC_OM_Wrapper
