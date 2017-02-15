
MODULE FISOC_ISM_Wrapper

  USE ESMF
  USE FISOC_utils_MOD

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: FISOC_ISM_Wrapper_Init_Phase1,  FISOC_ISM_Wrapper_Init_Phase2,  &
       FISOC_ISM_Wrapper_Run, FISOC_ISM_Wrapper_Finalize

  INTERFACE FISOC_ISM_Wrapper_Init_Phase1
      MODULE PROCEDURE FISOC_ISM_Wrapper_Init_Phase1_mesh
      MODULE PROCEDURE FISOC_ISM_Wrapper_Init_Phase1_grid
   END INTERFACE

  TYPE(ESMF_config)                     :: FOOL_config 

CONTAINS

  !--------------------------------------------------------------------------------------
  ! FOOL stands for FOrcing OffLine.  This wrapper pretends to be an ISM but actually 
  ! just reads in data from file and uses it to force the ocean. 
  ! 
  SUBROUTINE FISOC_ISM_Wrapper_Init_Phase1_mesh(ISM_ReqVarList,ISM_ExpFB,ISM_Mesh,&
       FISOC_config,vm,rc)

    CHARACTER(len=ESMF_MAXSTR),INTENT(IN) :: ISM_ReqVarList(:)
    TYPE(ESMF_config),INTENT(INOUT)       :: FISOC_config
    TYPE(ESMF_VM),INTENT(INOUT)           :: vm
    TYPE(ESMF_fieldBundle),INTENT(INOUT)  :: ISM_ExpFB
    TYPE(ESMF_mesh),INTENT(OUT)           :: ISM_Mesh
    INTEGER,INTENT(OUT),OPTIONAL          :: rc

    msg = "ERROR: Dummy subroutine called probably due to ISM_gridType error"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
         line=__LINE__, file=__FILE__, rc=rc)
    CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

  END SUBROUTINE FISOC_ISM_Wrapper_Init_Phase1_mesh


  SUBROUTINE FISOC_ISM_Wrapper_Init_Phase1_grid(ISM_ReqVarList,ISM_ExpFB,ISM_Grid,&
       FISOC_config,vm,rc)

    CHARACTER(len=ESMF_MAXSTR),INTENT(IN) :: ISM_ReqVarList(:)
    TYPE(ESMF_config),INTENT(INOUT)       :: FISOC_config
    TYPE(ESMF_VM),INTENT(INOUT)           :: vm
    TYPE(ESMF_fieldBundle),INTENT(INOUT)  :: ISM_ExpFB
    TYPE(ESMF_grid),INTENT(OUT)           :: ISM_Grid
    INTEGER,INTENT(OUT),OPTIONAL          :: rc

    CHARACTER(len=ESMF_MAXSTR)            :: FOOL_configName, ISM_gridLayout
    INTEGER                               :: localPet
    LOGICAL                               :: verbose_coupling

    rc = ESMF_FAILURE

    CALL ESMF_ConfigGetAttribute(FISOC_config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! Load FOOL configuration file (get name from FISOC config)
    CALL ESMF_ConfigGetAttribute(FISOC_config, FOOL_configName, label='ISM_configFile:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
!    CALL ESMF_ConfigGetAttribute(FISOC_config, ISM_configFile_FISOC, label='ISM_configFile:', rc=rc)
!    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!         line=__LINE__, file=__FILE__)) &
!         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
!    FOOL_configName = ""
!    FOOL_configName = ISM_configFile_FISOC 
    FOOL_config = ESMF_ConfigCreate(rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_ConfigLoadFile(FOOL_config, FOOL_configName, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(FOOL_config, ISM_gridLayout, label='ISM_gridLayout:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    SELECT CASE(ISM_gridLayout)

    CASE ('isomip_plus')
       msg = "Creating ISM ISOMIP+ grid"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)

       CALL Create_ISOMIP_plus_grid(ISM_grid,rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       
    CASE DEFAULT
       msg = "ERROR: FOOL does not recognise ISM_gridLayout"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       
    END SELECT
    
    CALL ESMF_VMGet(vm, localPet=localPet, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       PRINT*,""
       PRINT*,"*****************************************************************************"
       PRINT*,"****   ISM offline forcing wrapper.  Init phase 1 method.    ****************"
       PRINT*,"*****************************************************************************"
       PRINT*,""
       PRINT*,"Here we need to get the ISM grid information into the ESMF_grid type. "
       PRINT*,"We also need to create and initialise the required variables using the "
       PRINT*,"ESMF_field type and put them into an ESMF_fieldBundle type."
       PRINT*,""
    END IF

    CALL FISOC_populateFieldBundle(ISM_ReqVarList,ISM_ExpFB,ISM_Grid,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_ISM_Wrapper_Init_Phase1_grid
  

  !--------------------------------------------------------------------------------------
  ! This dummy wrapper aims to create the dummy grid and required variables 
  ! in the ESMF formats.  
  SUBROUTINE FISOC_ISM_Wrapper_Init_Phase2(ISM_ImpFB,ISM_ExpFB,FISOC_config,vm,rc)

    TYPE(ESMF_fieldBundle),INTENT(INOUT)  :: ISM_ImpFB, ISM_ExpFB
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
       PRINT*,"**********    ISM dummy wrapper.  Init phase 2 method.    ********************"
       PRINT*,"******************************************************************************"
       PRINT*,""
       PRINT*,"Here we have access to the initialised OM fields, just in case the ISM needs "
       PRINT*,"to know about these in order to complete its initialisation."
       PRINT*,""
    END IF
    
    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_ISM_Wrapper_Init_Phase2


  !--------------------------------------------------------------------------------------
  SUBROUTINE FISOC_ISM_Wrapper_Run(FISOC_config,vm,ISM_ExpFB,ISM_ImpFB,rc)

    TYPE(ESMF_config)      :: FISOC_config
    TYPE(ESMF_fieldbundle) :: ISM_ImpFB,ISM_ExpFB
    TYPE(ESMF_VM)          :: vm
    INTEGER,INTENT(OUT),OPTIONAL :: rc

    INTEGER                      :: localPet
    TYPE(ESMF_field)             :: OM_dBdt_l0, ISM_z_l0, ISM_z_l1
    REAL(ESMF_KIND_R8),POINTER   :: OM_dBdt_l0_ptr(:),ISM_z_l0_ptr(:),ISM_z_l1_ptr(:)
    LOGICAL                      :: verbose_coupling

    rc = ESMF_FAILURE

    CALL ESMF_VMGet(vm, localPet=localPet, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! query the FISOC config
    CALL ESMF_ConfigGetAttribute(FISOC_config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       PRINT*,""
       PRINT*,"******************************************************************************"
       PRINT*,"*************       ISM dummy wrapper.  Run method.       ********************"
       PRINT*,"******************************************************************************"
       PRINT*,""
       PRINT*,"OM export fields are available.  Run the ISM and return ISM export fields "
       PRINT*,""
    END IF

    ! get import and export fields and do something with them.
    CALL ESMF_FieldBundleGet(ISM_ImpFB, fieldName="OM_dBdt_l0", field=OM_dBdt_l0, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_FieldGet(field=OM_dBdt_l0, localDe=0, farrayPtr=OM_dBdt_l0_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       PRINT*,"Lets just adjust depth coords according to basal melt rate from ocean."
       PRINT*,"Melt rates dont look quite right on the ISM grid, needs checking..."
       PRINT*,OM_dBdt_l0_ptr
       PRINT*,""
    END IF

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_z_l0", field=ISM_z_l0, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_z_l1", field=ISM_z_l1, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_FieldGet(field=ISM_z_l0, localDe=0, farrayPtr=ISM_z_l0_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_z_l1, localDe=0, farrayPtr=ISM_z_l1_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ISM_z_l0_ptr = ISM_z_l0_ptr + OM_dBdt_l0_ptr
    ISM_z_l1_ptr = ISM_z_l1_ptr + OM_dBdt_l0_ptr

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_ISM_Wrapper_Run


  !--------------------------------------------------------------------------------------
  SUBROUTINE FISOC_ISM_Wrapper_Finalize(FISOC_config,localPet,rc)

    TYPE(ESMF_config),INTENT(INOUT)    :: FISOC_config
    INTEGER,INTENT(OUT),OPTIONAL       :: rc
    INTEGER,INTENT(IN)                 :: localPet

    LOGICAL                            :: verbose_coupling

    rc = ESMF_FAILURE

    CALL ESMF_ConfigGetAttribute(FISOC_config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       PRINT*,""
       PRINT*,"******************************************************************************"
       PRINT*,"************    ISM dummy wrapper.  Finalise method.     *********************"
       PRINT*,"******************************************************************************"
       PRINT*,""
       PRINT*,"FISOC has taken care of clearing up ESMF types.  Here we just need to call the "
       PRINT*,"ISM finalise method."
    END IF

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_ISM_Wrapper_Finalize


  !--------------------------------------------------------------------------------------
  SUBROUTINE Create_ISOMIP_plus_grid(ISM_grid,rc)

    TYPE(ESMF_grid),INTENT(INOUT)  :: ISM_grid
    INTEGER,INTENT(OUT),OPTIONAL   :: rc

    CHARACTER(len=ESMF_MAXSTR)     :: fileName, fileNumber, ForcingDir, ForcingBaseName

    TYPE(ESMF_grid)                :: FOOLgrid
    TYPE(ESMF_field)               :: field, ISM_z_l0
    REAL(ESMF_KIND_R8),POINTER     :: xCoords(:), yCoords(:), field_arr(:,:), zice_arr(:,:)

    INTEGER                        :: lbnd(2), ubnd(2), year
    INTEGER                        :: nx, ny, dx, dy, ii, x1, y1

    rc = ESMF_FAILURE

    ! ISOMIP grid:
    ! 480 (x direction) by 80 (y direction) grid cells, 1000m resolution.
    ! domain goes from 320km to 800km in the x direction, 0 to 80km in y.
    ! nx = 480; ny = 80; dx = 1000; dy = 1000
    !  xCoords        = (/(ii, ii=x1, x1+nx*dx, dx)/)
    !  yCoords        = (/(ii, ii=y1, y1+ny*dy, dy)/)
    
    ! Note: we are using David G's processed netcdf files with 2km instead of 1km 
    ! resolution, so set the values accordingly here:
    nx = 240; ny = 40; dx = 2000; dy = 2000
    x1 = 321000; y1 = 1000 ! starting coords (like all) are at cell centres
    
    ! read from FOOL config file
    CALL ESMF_ConfigGetAttribute(FOOL_config, ForcingDir, label='ForcingDir:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    year = 1
!    ForcingDir      = "/short/ks3/rmg581/FISOC/examples/Ex3_offlineISM/ocn3Forcing/"
    ForcingBaseName = "isomip_plus_ocean3_"
    
!  NumForcingFiles:
!  ForcingInterval_yr:
!  ForcingStartYr:

    ! use FOOL_vars from config for varNames
    
    WRITE (fileNumber, "(I0)") year
    fileName = TRIM(ADJUSTL(ForcingDir))//TRIM(ADJUSTL(ForcingBaseName))//TRIM(ADJUSTL(fileNumber))//".nc"
    
    ! Make the grid and add coords
    FOOLgrid=ESMF_GridCreateNoPeriDim(          &
         maxIndex=(/ny,nx/), & 
         coordSys=ESMF_COORDSYS_CART, &
         coordDep1=(/1/), & ! 1st coord is 1D and depends on 1st Grid dim
         coordDep2=(/2/), & ! 2nd coord is 1D and depends on 2nd Grid dim
         indexflag=ESMF_INDEX_GLOBAL, &
         rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_GridAddCoord(FOOLgrid,  & 
         staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_GridGetCoord(FOOLgrid, 2, localDE=0, &
         staggerloc=ESMF_STAGGERLOC_CENTER, &
         farrayPtr=xCoords, &
         computationalLBound=lbnd, &
         computationalUBound=ubnd, &
         rc=rc)  
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    DO ii=lbnd(1),ubnd(1)
       xCoords(ii) = x1 - dx + (ii*dx)
    END DO
    
    CALL ESMF_GridGetCoord(FOOLgrid, 1, localDE=0, &
         staggerloc=ESMF_STAGGERLOC_CENTER, &
         farrayPtr=yCoords, &
         computationalLBound=lbnd, &
         computationalUBound=ubnd, &
         rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    DO ii=lbnd(1),ubnd(1)
       yCoords(ii) = y1 - dy + (ii*dy)
    END DO
    
    rc = ESMF_SUCCESS

  END SUBROUTINE Create_ISOMIP_plus_grid
    
END MODULE FISOC_ISM_Wrapper
