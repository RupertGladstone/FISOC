
!
! This is the ocean model spcific code for FISOC.  The main purpose is to transfer information 
! between ESMF structures and the OM's internal structures.
!

MODULE FISOC_OM_Wrapper

  USE ESMF
  USE FISOC_utils_MOD
  USE FISOC_types_MOD
  USE ocean_control_mod

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: FISOC_OM_Wrapper_Init_Phase1,  FISOC_OM_Wrapper_Init_Phase2,  &
       FISOC_OM_Wrapper_Run, FISOC_OM_Wrapper_Finalize

  !-----------------------------------------------------------------------
  !     Staggered grid point indices
  !     d --------- d   d --- v --- d  
  !     |           |   |           |
  !     |     c     |   u     c     u
  !     |           |   |           |
  !     d --------- d   d --- v --- d     
  !     Arakawa - B     Arakawa - C
  !     RegCM           ROMS (c = rho, d = psi)
  !-----------------------------------------------------------------------
  !
  character(len=6)   :: GRIDDES(0:4) = &
       (/'N/A   ','CROSS ','DOT   ','U     ','V     '/)
  integer, parameter :: Inan    = 0
  integer, parameter :: Icross  = 1
  integer, parameter :: Idot    = 2
  integer, parameter :: Iupoint = 3
  integer, parameter :: Ivpoint = 4
  
CONTAINS
  
  !--------------------------------------------------------------------------------------
  ! The first phase of initialisation is mainly to initialise the ocean model, and access 
  ! grid and variable initial information.
  SUBROUTINE FISOC_OM_Wrapper_Init_Phase1(OM_ExpFB,OM_grid,FISOC_config,vm,rc)

    TYPE(ESMF_config),INTENT(INOUT)       :: FISOC_config
    TYPE(ESMF_VM),INTENT(IN)              :: vm ! ESMF virtual machine (parallel context)
    TYPE(ESMF_grid),INTENT(OUT)           :: OM_grid
    TYPE(ESMF_fieldBundle),INTENT(INOUT)  :: OM_ExpFB
    INTEGER,INTENT(OUT),OPTIONAL          :: rc

    INTEGER                           :: localPet ! local persistent execution thread (1:1 relationship to process)
    INTEGER                               :: mpic ! mpi comm, duplicate from the OM VM
    CHARACTER(len=ESMF_MAXSTR)            :: label
    TYPE(ESMF_staggerLoc),ALLOCATABLE     :: OM_ReqVars_stagger(:)
    CHARACTER(len=ESMF_MAXSTR),ALLOCATABLE:: OM_ReqVarList(:),FISOC_OM_ReqVarList(:)
    CHARACTER(len=ESMF_MAXSTR)            :: OM_configFile
    LOGICAL                               :: verbose_coupling, first

    first = .TRUE.

    CALL ESMF_ConfigGetAttribute(FISOC_config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(FISOC_config, OM_configFile, label='OM_configFile:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_VMGet(vm, localPet=localPet, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL FISOC_VM_MPI_Comm_dup(vm,mpic,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF (mpic.EQ.FISOC_mpic_missing) THEN
       msg = "ERROR: not currently configured for serial ROMS simulations"
       ! TODO: check whether ROMS needs a dummy mpic in serial configuration
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    ELSE
       CALL ROMS_initialize(first,mpic,OM_configFile)
    END IF
    
    ! extract a list of required ocean variables from the configuration object
    label = 'FISOC_OM_ReqVars:' ! the FISOC names for the vars
    CALL FISOC_getStringListFromConfig(FISOC_config, label, FISOC_OM_ReqVarList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    label = 'OM_ReqVars:' ! the OM names for the vars
    CALL FISOC_getStringListFromConfig(FISOC_config, label, OM_ReqVarList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    label = 'OM_ReqVars_stagger:' ! the OM names for stagger locations corresponding to the vars
    ALLOCATE(OM_ReqVars_stagger(SIZE(OM_ReqVarList)))
    CALL FISOC_ConfigDerivedAttribute(FISOC_config, OM_ReqVars_stagger, label, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! ESMF grid creation accesses OM grid information from ROMS modules
    CALL OM_createGrid(OM_grid, localPet, verbose_coupling, rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL FISOC_populateFieldBundle(FISOC_OM_ReqVarList,OM_ExpFB,OM_grid,init_value=FISOC_missingData,fieldStagger=OM_ReqVars_stagger,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL getFieldDataFromOM(OM_ExpFB,FISOC_config,vm,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    
  END SUBROUTINE FISOC_OM_Wrapper_Init_Phase1

  
  !--------------------------------------------------------------------------------------
  SUBROUTINE FISOC_OM_Wrapper_Init_Phase2(OM_ImpFB,FISOC_config,localPet,rc)

    TYPE(ESMF_config),INTENT(INOUT)       :: FISOC_config
    TYPE(ESMF_fieldBundle),INTENT(INOUT)  :: OM_ImpFB
    INTEGER,INTENT(OUT),OPTIONAL          :: rc
    INTEGER,INTENT(IN)                    :: localPet

    LOGICAL                               :: verbose_coupling

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
       PRINT*,"**********      OM wrapper.  Init phase 2 method.        *********************"
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
  SUBROUTINE FISOC_OM_Wrapper_Run(FISOC_config,localPet,OM_ExpFB,OM_ImpFB,rc)
    
    TYPE(ESMF_config),INTENT(INOUT)                :: FISOC_config
    TYPE(ESMF_fieldBundle),INTENT(INOUT),OPTIONAL  :: OM_ExpFB, OM_ImpFB 
    INTEGER,INTENT(IN)                             :: localPet
    INTEGER,INTENT(OUT),OPTIONAL                   :: rc

    LOGICAL                    :: verbose_coupling
    TYPE(ESMF_field)           :: ISM_dTdz_l0,ISM_z_l0, OM_dBdt_l0
    REAL(ESMF_KIND_R8),POINTER :: ISM_dTdz_l0_ptr(:,:), ISM_z_l0_ptr(:,:), OM_dBdt_l0_ptr(:,:)
    INTEGER                    :: OM_dt_sec
    REAL(ESMF_KIND_R8)         :: OM_dt_sec_float


    rc = ESMF_FAILURE
    
    CALL ESMF_ConfigGetAttribute(FISOC_config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

print *,"send field data to OM"
!    IF (PRESENT(OM_ImpFB)) THEN       
!       CALL sendFieldDataToOM(OM_ImpFB,FISOC_config,vm,rc=rc)
!       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!            line=__LINE__, file=__FILE__)) &
!            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
!    END IF

    CALL ESMF_ConfigGetAttribute(FISOC_config, OM_dt_sec, label='OM_dt_sec:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    OM_dt_sec_float = REAL(OM_dt_sec,ESMF_KIND_R8)


print*,"add timestep check: check ROMS timestep is consistent with FISOC time stepping"

!alternative coupling modes re asynchronous coupling: add to issues in github

!check how to call ROMS run (just one timestep)

CALL ROMS_run(OM_dt_sec_float)
    
print*,"get field data from OM"

!    IF (PRESENT(OM_ExpFB)) THEN
!       CALL getFieldDataFromOM(OM_ExpFB,FISOC_config,vm,rc=rc)
!       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!            line=__LINE__, file=__FILE__)) &
!            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
!    END IF

print*,"put this in the ISM derived vars?"
!use ISM timestep and depth and prev depth to populate rom iceshelfvar % DDDT

    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       PRINT*,""
       PRINT*,"******************************************************************************"
       PRINT*,"************         OM wrapper.  Run method.           **********************"
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

print*,"move this stuff to subroutines"

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
       CALL  ESMF_FieldBundleWrite(OM_ExpFB, "test.nc", status=ESMF_FILESTATUS_REPLACE,&
            iofmt=ESMF_IOFMT_NETCDF, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    
       
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

print*,"check for all the places where I use the netcdf writer... remove or use a sensible netcdf file name"
    
!    OM_dBdt_l0_ptr = 12345.6
!       NULLIFY(OM_dBdt_l0_ptr)
    END IF

    rc = ESMF_SUCCESS
    
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
       PRINT*,"************       OM wrapper.  Finalise method.        **********************"
       PRINT*,"******************************************************************************"
       PRINT*,""
       PRINT*,"FISOC has taken care of clearing up ESMF types.  Here we just need to call the "
       PRINT*,"OM finalise method."
    END IF

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_OM_Wrapper_Finalize


  !--------------------------------------------------------------------------------------
  ! update the fields in the ocean export field bundle from the OM
  !--------------------------------------------------------------------------------------
  SUBROUTINE getFieldDataFromOM(OM_ExpFB,FISOC_config,vm,rc)

    USE mod_iceshelfvar, ONLY : ICESHELFVAR
    USE mod_param, ONLY       : BOUNDS, Ngrids

    IMPLICIT NONE

    TYPE(ESMF_fieldBundle),INTENT(INOUT)  :: OM_ExpFB
    TYPE(ESMF_config),INTENT(INOUT)       :: FISOC_config
    INTEGER,INTENT(OUT),OPTIONAL          :: rc
    TYPE(ESMF_VM),INTENT(IN)              :: vm

    INTEGER                               :: fieldCount, localPet, petCount
    TYPE(ESMF_FIELD)                      :: FISOC_field
    TYPE(ESMF_Field),ALLOCATABLE          :: fieldList(:)
    CHARACTER(len=ESMF_MAXSTR)            :: fieldName
    REAL(ESMF_KIND_R8),POINTER            :: ptr(:,:)
    INTEGER                               :: IstrR, IendR, JstrR, JendR ! tile start and end coords
    INTEGER                               :: ii, jj, nn
!    INTEGER                               :: LBi, UBi, LBj, UBj ! tile start and end coords including halo

    rc = ESMF_FAILURE

    IF (Ngrids .ne. 1) THEN
       msg = "ERROR: not expecting multiple ROMS grids"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    CALL ESMF_VMGet(vm, localPet=localPet, petCount=petCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
         line=__LINE__, file=__FILE__)) RETURN
  
    ! get the tile position from the OM (a tile is a rectangular domain decomposition element)
    IstrR=BOUNDS(Ngrids)%IstrR(localPet)
    IendR=BOUNDS(Ngrids)%IendR(localPet)
    JstrR=BOUNDS(Ngrids)%JstrR(localPet)
    JendR=BOUNDS(Ngrids)%JendR(localPet)
    
!    ! get the tile (including halo) position from the OM
!    LBi = BOUNDS(Ngrids)%LBi(localPet)
!    UBi = BOUNDS(Ngrids)%UBi(localPet)
!    LBj = BOUNDS(Ngrids)%LBj(localPet)
!    UBj = BOUNDS(Ngrids)%UBj(localPet)
    
    ! get a list of fields and their names form the OM export field bundle
    fieldCount = 0
    CALL ESMF_FieldBundleGet(OM_ExpFB, fieldCount=fieldCount, rc=rc)
    ALLOCATE(fieldList(fieldCount))
    CALL ESMF_FieldBundleGet(OM_ExpFB, fieldList=fieldList, rc=rc)
    
    fieldLoop: DO nn = 1,fieldCount
       
       CALL ESMF_FieldGet(fieldList(nn), name=fieldName, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
            line=__LINE__, file=__FILE__)) RETURN
       CALL ESMF_FieldGet(fieldList(nn), farrayPtr=ptr, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
            line=__LINE__, file=__FILE__)) RETURN
              
!       if (localPet == 1) then
!          print*,"STUFF", nn, size(ptr),size(ptr(:,1)),size(ptr(1,:)),localPet
!          print*,IstrR, IendR, JstrR, JendR
!          print*,size(ICESHELFVAR(Ngrids) % m),size(ICESHELFVAR(Ngrids) % m(:,1)),size(ICESHELFVAR(Ngrids) % m(1,:))
!       end if
       
       ptr = FISOC_missingData
       
       SELECT CASE (TRIM(ADJUSTL(fieldName)))
          
       CASE ('OM_dBdt_l0')
          DO jj = JstrR, JendR
             DO ii = IstrR, IendR
                ptr(ii,jj) = ICESHELFVAR(1) % m(ii,jj)
             END DO
          END DO
          
       CASE ('OM_temperature_l0')
          DO jj = JstrR, JendR
             DO ii = IstrR, IendR
                ptr(ii,jj) = ICESHELFVAR(1) % Tb(ii,jj)
             END DO
          END DO
          !ICESHELFVAR(ng) % Tb(LBi:UBi,LBj:UBj)
          
       CASE DEFAULT
          msg = "ERROR: unknown variable"
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
               line=__LINE__, file=__FILE__, rc=rc)
          CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       END SELECT
       
       IF (ASSOCIATED(ptr)) THEN
          NULLIFY(ptr)
       END IF
       
    END DO fieldLoop
    
    rc = ESMF_SUCCESS
    
  END SUBROUTINE GetFieldDataFromOM
  

  !--------------------------------------------------------------------------------------
  subroutine OM_createGrid(OM_grid, localPet, verbose_coupling, rc)
    
    use mod_grid , only : GRID
    use mod_param, only : NtileI, NtileJ, BOUNDS, Lm, Mm, Ngrids
    
    implicit none
    
    type(ESMF_grid), intent(inout) :: OM_grid
    integer, intent(in)            :: localPet 
    logical, intent(in)            :: verbose_coupling
    integer, intent(out)           :: rc
    
    integer                        :: i2, j2, ii, jj, ng, nr, tile, localDECount
    integer                        :: IstrR, IendR, JstrR, JendR
    integer                        :: IstrU, IendU, JstrU, JendU     
    integer                        :: IstrV, IendV, JstrV, JendV
    integer                        :: LBi, UBi, LBj, UBj
    integer                        :: staggerEdgeLWidth(2)
    integer                        :: staggerEdgeUWidth(2)
    integer, allocatable           :: deBlockList(:,:,:)
    real(ESMF_KIND_R8), pointer    :: ptrX(:,:), ptrY(:,:), ptrA(:,:)
    integer(ESMF_KIND_I4), pointer :: ptrM(:,:)
    character(ESMF_MAXSTR)         :: name, msgString
    
    type(ESMF_Array)               :: arrX, arrY, arrM, arrA
    type(ESMF_StaggerLoc)          :: staggerLoc
    type(ESMF_DistGrid)            :: distGrid
        
    rc = ESMF_FAILURE
    
    if (Ngrids > 1) then
       msg = 'number of nested grid is greater than 1.'//        &
            'Coupling only interacts with outermost one!'
       call ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
       ng = 1
    else
       ng = Ngrids
    end if
    
    !-----------------------------------------------------------------------
    !     Get limits of the grid arrays (based on PET and nest level)
    !-----------------------------------------------------------------------
    
    IstrR = BOUNDS(ng)%IstrR(localPet)
    IendR = BOUNDS(ng)%IendR(localPet)
    JstrR = BOUNDS(ng)%JstrR(localPet)
    JendR = BOUNDS(ng)%JendR(localPet)
    !
    IstrU = BOUNDS(ng)%Istr(localPet)
    IendU = BOUNDS(ng)%IendR(localPet)
    JstrU = BOUNDS(ng)%JstrR(localPet)
    JendU = BOUNDS(ng)%JendR(localPet)
    !
    IstrV = BOUNDS(ng)%IstrR(localPet)
    IendV = BOUNDS(ng)%IendR(localPet)
    JstrV = BOUNDS(ng)%Jstr(localPet)
    JendV = BOUNDS(ng)%JendR(localPet)
    !
    LBi = BOUNDS(ng)%LBi(localPet)
    UBi = BOUNDS(ng)%UBi(localPet)
    LBj = BOUNDS(ng)%LBj(localPet)
    UBj = BOUNDS(ng)%UBj(localPet)
    !
    if (.not.allocated(deBlockList)) then
       allocate(deBlockList(2,2,NtileI(ng)*NtileJ(ng)))
    end if
    do tile=0,NtileI(ng)*NtileJ(ng)-1
       deBlockList(1,1,tile+1)=BOUNDS(ng)%Istr(tile)
       deBlockList(1,2,tile+1)=BOUNDS(ng)%Iend(tile)
       deBlockList(2,1,tile+1)=BOUNDS(ng)%Jstr(tile)
       deBlockList(2,2,tile+1)=BOUNDS(ng)%Jend(tile)
    end do
    
    !-----------------------------------------------------------------------
    !     Create ESMF DistGrid based on ROMS model domain decomposition
    !-----------------------------------------------------------------------
    distGrid = ESMF_DistGridCreate(minIndex=(/ 1, 1 /),               &
         maxIndex=(/ Lm(ng), Mm(ng) /),     &
         deBlockList=deBlockList,           &
         rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
         line=__LINE__, file=__FILE__)) return
    
    if (allocated(deBlockList)) deallocate(deBlockList) 
    
    
    !-----------------------------------------------------------------------
    !     Set staggering type 
    !-----------------------------------------------------------------------
    staggerLoop: do ii = 1, 4 
       if (ii == Iupoint) then
          staggerLoc = ESMF_STAGGERLOC_EDGE1
          staggerEdgeLWidth = (/0,1/)
          staggerEdgeUWidth = (/1,1/)
       else if (ii == Ivpoint) then
          staggerLoc = ESMF_STAGGERLOC_EDGE2
          staggerEdgeLWidth = (/1,0/)
          staggerEdgeUWidth = (/1,1/)
       else if (ii == Icross) then
          staggerLoc = ESMF_STAGGERLOC_CENTER
          staggerEdgeLWidth = (/1,1/)
          staggerEdgeUWidth = (/1,1/)
       else if (ii == Idot) then
          staggerLoc = ESMF_STAGGERLOC_CORNER
          staggerEdgeLWidth = (/0,0/)
          staggerEdgeUWidth = (/1,1/)
       end if
       
       !-----------------------------------------------------------------------
       !     Create ESMF Grid
       !-----------------------------------------------------------------------
       if (ii == 1) then
          OM_grid = ESMF_GridCreate(distgrid=distGrid,  &
               gridEdgeLWidth=(/1,1/),                              &
               gridEdgeUWidth=(/1,1/),                              &
               indexflag=ESMF_INDEX_GLOBAL,                         &
               name="OM_grid",                                      &
               rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,  &
               line=__LINE__, file=__FILE__))                             &
               call ESMF_Finalize(endflag=ESMF_END_ABORT)
       end if
       
       !-----------------------------------------------------------------------
       !     Allocate coordinates 
       !-----------------------------------------------------------------------
       call ESMF_GridAddCoord(OM_grid,            &
            staggerLoc=staggerLoc,                &
            staggerEdgeLWidth=staggerEdgeLWidth,  &
            staggerEdgeUWidth=staggerEdgeUWidth,  &
            rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
            line=__LINE__, file=__FILE__)) return
       
       !-----------------------------------------------------------------------
       !     Allocate items for masking
       !-----------------------------------------------------------------------
       call ESMF_GridAddItem(OM_grid,       &
            staggerLoc=staggerLoc,          &
            itemflag=ESMF_GRIDITEM_MASK,    &
            rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
            line=__LINE__, file=__FILE__)) return
       
       !-----------------------------------------------------------------------
       !     Allocate items for grid area 
       !-----------------------------------------------------------------------
       call ESMF_GridAddItem(OM_grid,       &
            staggerLoc=staggerLoc,          &
            itemflag=ESMF_GRIDITEM_AREA,    &
            rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
            line=__LINE__, file=__FILE__)) return
       
       !-----------------------------------------------------------------------
       !     Get number of local DEs
       !-----------------------------------------------------------------------
       call ESMF_GridGet(OM_grid,              &
            localDECount=localDECount,         &
            rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
            line=__LINE__, file=__FILE__)) return
       
       !-----------------------------------------------------------------------
       !     Get pointers and set coordinates for the grid 
       !-----------------------------------------------------------------------
       do jj = 0, localDECount-1
          call ESMF_GridGetCoord(OM_grid,                 &
               localDE=jj,                                &
               staggerLoc=staggerLoc,                     &
               coordDim=1,                                &
               farrayPtr=ptrX,                            &
               rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
               line=__LINE__, file=__FILE__)) return
          !
          call ESMF_GridGetCoord(OM_grid,                 &
               localDE=jj,                                &
               staggerLoc=staggerLoc,                     &
               coordDim=2,                                &
               farrayPtr=ptrY,                            &
               rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
               line=__LINE__, file=__FILE__)) return
          !
          call ESMF_GridGetItem (OM_grid,                 &
               localDE=jj,                                &
               staggerLoc=staggerLoc,                     &
               itemflag=ESMF_GRIDITEM_MASK,               &
               farrayPtr=ptrM,                            &
               rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
               line=__LINE__, file=__FILE__)) return
          !
          call ESMF_GridGetItem (OM_grid,                 &
               localDE=jj,                                &
               staggerLoc=staggerLoc,                     &
               itemflag=ESMF_GRIDITEM_AREA,               &
               farrayPtr=ptrA,                            &
               rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
               line=__LINE__, file=__FILE__)) return
          !
          !-----------------------------------------------------------------------
          !     Debug: write size of pointers    
          !-----------------------------------------------------------------------
          !
          name = GRIDDES(ii)
          !
          if (verbose_coupling) then
             write(*,30) localPet, jj, adjustl("PTR/OCN/GRD/"//name), &
                  lbound(ptrX, dim=1), ubound(ptrX, dim=1),           &
                  lbound(ptrX, dim=2), ubound(ptrX, dim=2)
          end if
          !
          !-----------------------------------------------------------------------
          !     Fill the pointers    
          !-----------------------------------------------------------------------
          !
          if (ii == Idot) then
             if (verbose_coupling) then
                write(*,30) localPet, jj, adjustl("DAT/OCN/GRD/"//name),         &
                     lbound(GRID(ng)%lonp, dim=1), ubound(GRID(ng)%lonp, dim=1),     &
                     lbound(GRID(ng)%lonp, dim=2), ubound(GRID(ng)%lonp, dim=2)
             end if
             !
             do j2 = JstrV, JendR
                do i2 = IstrU, IendR
                   ptrX(i2,j2) = GRID(ng)%lonp(i2,j2)
                   ptrY(i2,j2) = GRID(ng)%latp(i2,j2)
                   ptrM(i2,j2) = int(GRID(ng)%pmask(i2,j2))
                   ptrA(i2,j2) = GRID(ng)%om_p(i2,j2)*GRID(ng)%on_p(i2,j2)
                end do
             end do
          else if (ii == Icross) then
             if (verbose_coupling) then
                write(*,30) localPet, jj, adjustl("DAT/OCN/GRD/"//name),         &
                     lbound(GRID(ng)%lonr, dim=1), ubound(GRID(ng)%lonr, dim=1),     &
                     lbound(GRID(ng)%lonr, dim=2), ubound(GRID(ng)%lonr, dim=2)
             end if
             !
             do j2 = JstrR, JendR
                do i2 = IstrR, IendR
                   ptrX(i2,j2) = GRID(ng)%lonr(i2,j2)
                   ptrY(i2,j2) = GRID(ng)%latr(i2,j2)
                   ptrM(i2,j2) = int(GRID(ng)%rmask(i2,j2))
                   ptrA(i2,j2) = GRID(ng)%om_r(i2,j2)*GRID(ng)%on_r(i2,j2)
                end do
             end do
          else if (ii == Iupoint) then
             if (verbose_coupling) then
                write(*,30) localPet, jj, adjustl("DAT/OCN/GRD/"//name),         &
                     lbound(GRID(ng)%lonu, dim=1), ubound(GRID(ng)%lonu, dim=1),     &
                     lbound(GRID(ng)%lonu, dim=2), ubound(GRID(ng)%lonu, dim=2)
             end if
             !
             do j2 = JstrU, JendU
                do i2 = IstrU, IendU
                   ptrX(i2,j2) = GRID(ng)%lonu(i2,j2)
                   ptrY(i2,j2) = GRID(ng)%latu(i2,j2)
                   ptrM(i2,j2) = int(GRID(ng)%umask(i2,j2))
                   ptrA(i2,j2) = GRID(ng)%om_u(i2,j2)*GRID(ng)%on_u(i2,j2)
                end do
             end do
          else if (ii == Ivpoint) then
             if (verbose_coupling) then
                write(*,30) localPet, jj, adjustl("DAT/OCN/GRD/"//name),         &
                     lbound(GRID(ng)%lonv, dim=1), ubound(GRID(ng)%lonv, dim=1),     &
                     lbound(GRID(ng)%lonv, dim=2), ubound(GRID(ng)%lonv, dim=2)
             end if
             !
             do j2 = JstrV, JendV
                do i2 = IstrV, IendV
                   ptrX(i2,j2) = GRID(ng)%lonv(i2,j2)
                   ptrY(i2,j2) = GRID(ng)%latv(i2,j2)
                   ptrM(i2,j2) = int(GRID(ng)%vmask(i2,j2))
                   ptrA(i2,j2) = GRID(ng)%om_v(i2,j2)*GRID(ng)%on_v(i2,j2)
                end do
             end do
          end if
          !
          !-----------------------------------------------------------------------
          !     Create temporary arrays.
          !-----------------------------------------------------------------------
          !
          if (ii == Icross) then
             arrX = ESMF_ArrayCreate(distGrid, ptrX,                         &
                  indexflag=ESMF_INDEX_DELOCAL, rc=rc) 
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,  &
                  line=__LINE__, file=__FILE__)) return
             !
             arrY = ESMF_ArrayCreate(distGrid, ptrY,                         &
                  indexflag=ESMF_INDEX_DELOCAL, rc=rc) 
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,  &
                  line=__LINE__, file=__FILE__)) return
             !
             arrM = ESMF_ArrayCreate(distGrid, ptrM,                         &
                  indexflag=ESMF_INDEX_DELOCAL, rc=rc) 
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,  &
                  line=__LINE__, file=__FILE__)) return
             !
             arrA = ESMF_ArrayCreate(distGrid, ptrA,                         &
                  indexflag=ESMF_INDEX_DELOCAL, rc=rc) 
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,  &
                  line=__LINE__, file=__FILE__)) return
          end if
          !
          !-----------------------------------------------------------------------
          !     Nullify pointers 
          !-----------------------------------------------------------------------
          !
          if (associated(ptrX)) then
             nullify(ptrX)
          end if
          if (associated(ptrY)) then
             nullify(ptrY)
          end if
          if (associated(ptrM)) then
             nullify(ptrM)
          end if
          if (associated(ptrA)) then
             nullify(ptrA)
          end if
       end do

       !-----------------------------------------------------------------------
       !     Debug: write out component grid in VTK format 
       !-----------------------------------------------------------------------
       !
       if (verbose_coupling) then
          call ESMF_GridWriteVTK(OM_grid,                     &
               filename="ocean_"//                            &
               trim(GRIDDES(ii))//                            &
               "point",                                       &
               staggerLoc=staggerLoc,                         &
               rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
               line=__LINE__, file=__FILE__)) return
       end if
    end do staggerLoop
    
30  format(" PET(",I3.3,") - DE(",I2.2,") - ", A20, " : ", 4I8)
    !
  end subroutine OM_createGrid

!       ! initialise field data to MISSING_DATA
!       do jj = 0, localDECount-1
!          !-----------------------------------------------------------------------
!          !     Get pointer to data array from field 
!          !-----------------------------------------------------------------------
!          call ESMF_FieldGet(field, localDe=jj, farrayPtr=ptr2d, rc=rc)
!          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,    &
!               line=__LINE__, file=__FILE__)) return
!          ptr2d = MISSING_R8
!          if (associated(ptr2d)) then
!             nullify(ptr2d)
!          end if

END MODULE FISOC_OM_Wrapper
