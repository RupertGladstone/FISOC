
MODULE FISOC_ISM_Wrapper

  USE ESMF
  USE Netcdf
  USE FISOC_utils_MOD

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: FISOC_ISM_Wrapper_Init_Phase1,  FISOC_ISM_Wrapper_Init_Phase2,  &
       FISOC_ISM_Wrapper_Run, FISOC_ISM_Wrapper_Finalize

  INTERFACE FISOC_ISM_Wrapper_Init_Phase1
     MODULE PROCEDURE FISOC_ISM_Wrapper_Init_Phase1_mesh
     MODULE PROCEDURE FISOC_ISM_Wrapper_Init_Phase1_grid
  END INTERFACE FISOC_ISM_Wrapper_Init_Phase1
  
  INTEGER,PARAMETER                     :: NOYEAR = -999
  
  TYPE(ESMF_config)                     :: FOOL_config 
!  REAL(ESMF_KIND_R8),PARAMETER          :: yearinsec = 365.25*24.*60.*60. 
  INTEGER                               :: year = NOYEAR

! TODO: resolve some code duplication between init 1 and run.  Setting up file names etc. for getting the var from netcdf.
! TODO: hard code timestep check instead of hard coding timestep itself.  Put expected timestep in config file.

CONTAINS

  !--------------------------------------------------------------------------------------
  ! FOOL stands for Forcing OffLine.  This wrapper pretends to be an ISM but actually 
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
    INTEGER                               :: localPet, ISM_dt_sec
    LOGICAL                               :: verbose_coupling

    CHARACTER(len=ESMF_MAXSTR)   :: fileName
    TYPE(ESMF_grid)              :: FOOLgrid
    TYPE(ESMF_field)             :: field
    REAL(ESMF_KIND_R8),POINTER   :: ptr(:,:)

    rc = ESMF_FAILURE

    CALL ESMF_ConfigGetAttribute(FISOC_config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL FISOC_ConfigDerivedAttribute(FISOC_config, ISM_dt_sec, 'ISM_dt_sec',rc=rc)
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

    ! Get the fields needed from our ISM, in this case just the netcdf file
    CALL getFieldDataFromISM(ISM_ExpFB,FISOC_config)
    
    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_ISM_Wrapper_Init_Phase1_grid
  

  !--------------------------------------------------------------------------------------
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
       PRINT*,"**********    ISM FOOL wrapper.  Init phase 2 method.    ********************"
       PRINT*,"******************************************************************************"
       PRINT*,""
       PRINT*,"Here we have access to the initialised OM fields, just in case the ISM needs "
       PRINT*,"to know about these in order to complete its initialisation."
       PRINT*,""
    END IF
    
    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_ISM_Wrapper_Init_Phase2


  !-------------------------------------------------------------------------------
  ! 
  SUBROUTINE FISOC_ISM_Wrapper_Run(FISOC_config,vm,ISM_ExpFB,ISM_ImpFB,rc)

    TYPE(ESMF_config),INTENT(INOUT)      :: FISOC_config
    TYPE(ESMF_fieldbundle),INTENT(INOUT) :: ISM_ImpFB,ISM_ExpFB
    TYPE(ESMF_VM),INTENT(IN)             :: vm
    INTEGER,INTENT(OUT),OPTIONAL         :: rc

    CHARACTER(len=ESMF_MAXSTR)   :: fileName
    INTEGER                      :: localPet
    INTEGER                      :: rank, ISM_dt_sec
    TYPE(ESMF_grid)              :: FOOLgrid
    TYPE(ESMF_field)             :: field
    LOGICAL                      :: verbose_coupling

    rc = ESMF_FAILURE

    CALL ESMF_VMGet(vm, localPet=localPet, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! We expect grids not meshes here, thus rank = dim = 2
    CALL FISOC_getFirstFieldRank(ISM_ExpFB,rank,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    IF (rank.NE.2) THEN
       msg = "FOOL wrapper expecting ISM_ExpFB rank 2"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF
    CALL FISOC_getFirstFieldRank(ISM_ImpFB,rank,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    IF (rank.NE.2) THEN
       msg = "FOOL wrapper expecting ISM_ImpFB rank 2"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    ! query the FISOC config
    CALL ESMF_ConfigGetAttribute(FISOC_config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       PRINT*,""
       PRINT*,"******************************************************************************"
       PRINT*,"*************       ISM FOOL wrapper.  Run method.       ********************"
       PRINT*,"******************************************************************************"
       PRINT*,""
       PRINT*,"OM export fields are available but for offline forcing we dont need them. "
       PRINT*,"Just read new forcing data to pass to OM."
       PRINT*,""
    END IF

    ! Get the fields needed from our ISM, in this case just the netcdf file
    CALL getFieldDataFromISM(ISM_ExpFB,FISOC_config)

  END SUBROUTINE FISOC_ISM_Wrapper_Run


  
  !--------------------------------------------------------------------------------------                                                                                                                                                     
  SUBROUTINE getFieldDataFromISM(ISM_ExpFB,FISOC_config)

    TYPE(ESMF_fieldBundle),INTENT(INOUT)     :: ISM_ExpFB 
    TYPE(ESMF_config),INTENT(INOUT)          :: FISOC_config

    TYPE(ESMF_grid)                :: FOOLgrid
    INTEGER                        :: nn, ISM_dt_sec
    REAL(ESMF_KIND_R8),POINTER     :: ptr(:,:)
    INTEGER                        :: fieldCount, rc
    TYPE(ESMF_Field),ALLOCATABLE   :: fieldList(:)
    CHARACTER(len=ESMF_MAXSTR)     :: fieldName, fileName

    CALL FISOC_ConfigDerivedAttribute(FISOC_config, ISM_dt_sec, 'ISM_dt_sec',rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL makeFileName(FOOL_config,fileName)

    ! get a list of fields and their names from the ISM export field bundle
    fieldCount = 0
    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldCount=fieldCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    ALLOCATE(fieldList(fieldCount))
    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldList=fieldList, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! We can use the first field to get hold of the grid
    CALL ESMF_FieldGet(fieldList(1), grid=FOOLgrid, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)


    fieldLoop: DO nn = 1,fieldCount
       
       ! access the FISOC version of the current field
       CALL ESMF_FieldGet(fieldList(nn), name=fieldName, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       CALL ESMF_FieldGet(fieldList(nn), farrayPtr=ptr, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       ! access the netcdf file for the current field
       SELECT CASE (TRIM(ADJUSTL(fieldName)))

       CASE ('ISM_z_l0')
          CALL readFromNC(FileName,'zice',FOOLgrid,ptr,ISM_dt_sec)
          
       CASE ('ISM_dddt')
          CALL readFromNC(FileName,'dddt',FOOLgrid,ptr,ISM_dt_sec)
          
       CASE ('ISM_dsdt')
          CALL readFromNC(FileName,'dsdt',FOOLgrid,ptr,ISM_dt_sec)

       CASE DEFAULT
          msg = "ERROR: unknown variable: "//TRIM(ADJUSTL(fieldName))
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
               line=__LINE__, file=__FILE__, rc=rc)
          CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       END SELECT

    END DO fieldLoop

    NULLIFY(ptr)

    rc = ESMF_SUCCESS

  END SUBROUTINE getFieldDataFromISM


!    CALL ESMF_ClockGet(clock, startTime, currTime
!    advanceCount, rc=rc)
!    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!         line=__LINE__, file=__FILE__)) &
!         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

!    CALL ESMF_TimeGet(time, yy, rc=rc)
!    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!         line=__LINE__, file=__FILE__)) &
!         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

  


  SUBROUTINE makeFileName(FOOL_config,fileName)

    TYPE(ESMF_config),INTENT(INOUT)        :: FOOL_config
    CHARACTER(len=ESMF_MAXSTR),INTENT(OUT) :: fileName

    INTEGER                      :: NumForcingFiles, rc
    INTEGER                      :: ForcingInterval_yr, ForcingStartYr
    CHARACTER(len=ESMF_MAXSTR)   :: ForcingBaseName 
    CHARACTER(len=ESMF_MAXSTR)   :: ForcingDir, fileNumber

    
    ! set up file name for this time interval
    CALL ESMF_ConfigGetAttribute(FOOL_config, ForcingDir, label='ForcingDir:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(FOOL_config, ForcingBaseName, label='ForcingBaseName:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(FOOL_config, NumForcingFiles, label='NumForcingFiles:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(FOOL_config, ForcingInterval_yr, label='ForcingInterval_yr:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(FOOL_config, ForcingStartYr, label='ForcingStartYr:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF (year.eq.NOYEAR) THEN
       year = ForcingStartYr
    END IF

    WRITE (fileNumber, "(I0)") year
    fileName = TRIM(ADJUSTL(ForcingDir))//TRIM(ADJUSTL(ForcingBaseName))//TRIM(ADJUSTL(fileNumber))//".nc"

    ! Increment year for next time (we assume here the ISM timestep is 1 year)
    year = year + 1

  END SUBROUTINE makeFileName



  SUBROUTINE readFromNC(FileName,VarName,FOOLgrid,ptr,ISM_dt_sec)

    REAL(ESMF_KIND_R8),POINTER,INTENT(INOUT) :: ptr(:,:)
    CHARACTER(len=ESMF_MAXSTR),INTENT(IN)    :: FileName
    CHARACTER(len=*),INTENT(IN)              :: VarName
    TYPE(ESMF_grid),INTENT(INOUT)            :: FOOLgrid
    INTEGER,INTENT(IN)                       :: ISM_dt_sec

    REAL(ESMF_KIND_R8),ALLOCATABLE :: values(:,:)
    INTEGER                        :: lbnd(2), ubnd(2), NtileI, NtileJ
    INTEGER                        :: lbx, ubx, lby, uby, nx, ny
    INTEGER                        :: status, ncid, varid, rc

    msg = "NC file: "//fileName
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

    ! Get the bounds for the local pet
    CALL ESMF_GridGetCoordBounds(FOOLgrid, 1,    &
         totalLBound=lbnd, totalUBound=ubnd,   rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    lbx = lbnd(1)
    ubx = ubnd(1)

    CALL ESMF_GridGetCoordBounds(FOOLgrid, 2,    &
         totalLBound=lbnd, totalUBound=ubnd,   rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    lby = lbnd(1)
    uby = ubnd(1)
    nx = ubx - lbx + 1
    ny = uby - lby + 1

    ALLOCATE(values(ny,nx))

    ! Extract netcdf variable.  Var names hard coded here.
    status = nf90_open(fileName, NF90_NOWRITE, ncid)
    IF(status /= nf90_NoErr) CALL handle_err(status)
    
    status = nf90_inq_varid(ncid, VarName, varid)
    IF(status /= nf90_NoErr) CALL handle_err(status)

    status = nf90_get_var(ncid, varid, values,  &
         start = (/ lby, lbx /),                      &
         count = (/ ny,  nx  /)                        )
    IF(status /= nf90_NoErr) CALL handle_err(status)

    status = nf90_close(ncid)
    IF(status /= nf90_NoErr) CALL handle_err(status)

    ptr = TRANSPOSE(values)

!    ptr = ptr / ISM_dt_sec ! convert from m/yr to m/s
ptr = ptr/31557600.0

    DEALLOCATE(values)

    rc = ESMF_SUCCESS

  CONTAINS
    SUBROUTINE handle_err(status)
      INTEGER, INTENT(IN) :: status
      
      IF(status /= nf90_noerr) THEN
         msg = "FISOC encountered NETCDF error.  Error text follows."
         CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
              line=__LINE__, file=__FILE__, rc=rc)
         msg = trim(nf90_strerror(status))           
         CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
              line=__LINE__, file=__FILE__, rc=rc)
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
      END IF
    END SUBROUTINE handle_err

  END SUBROUTINE readFromNC




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
       PRINT*,"************    ISM FOOL wrapper.  Finalise method.     *********************"
       PRINT*,"******************************************************************************"
       PRINT*,""
       PRINT*,"FISOC has taken care of clearing up ESMF types.  Here we just need to call the "
       PRINT*,"ISM finalise method."
    END IF

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_ISM_Wrapper_Finalize


  !--------------------------------------------------------------------------------------
  SUBROUTINE Create_ISOMIP_plus_grid(FOOLgrid,rc)

    TYPE(ESMF_grid),INTENT(INOUT)  :: FOOLgrid
    INTEGER,INTENT(OUT),OPTIONAL   :: rc

    CHARACTER(len=ESMF_MAXSTR)     :: fileName, fileNumber, ForcingDir, ForcingBaseName
    REAL(ESMF_KIND_R8),POINTER     :: xCoords(:), yCoords(:)
    INTEGER                        :: lbnd(2), ubnd(2), year, NtileI, NtileJ
    INTEGER                        :: nx, ny, dx, dy, ii, x1, y1, localDEcount

    rc = ESMF_FAILURE

    ! ISOMIP grid:
    ! 480 (x direction) by 80 (y direction) grid cells, 1000m resolution.
    ! domain goes from 320km to 800km in the x direction, 0 to 80km in y.
    ! nx = 480; ny = 80; dx = 1000; dy = 1000
    !  xCoords        = (/(ii, ii=x1, x1+nx*dx, dx)/)
    !  yCoords        = (/(ii, ii=y1, y1+ny*dy, dy)/)
    
    ! Note: we are using David G's processed netcdf files with 2km instead of 1km 
    ! resolution, so set the values accordingly here:
    ny = 240; nx = 40; dx = 2000; dy = 2000
    y1 = 321000; x1 = 1000 ! starting coords (like all) are at cell centres
    
    ! get the decomposition from the FOOL config
    CALL ESMF_ConfigGetAttribute(FOOL_config, NtileI, label='NtileI:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_ConfigGetAttribute(FOOL_config, NtileJ, label='NtileJ:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! read from FOOL config file
    CALL ESMF_ConfigGetAttribute(FOOL_config, ForcingDir, label='ForcingDir:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    

    ! Make the grid and add coords
    FOOLgrid=ESMF_GridCreateNoPeriDim(  &
         maxIndex=(/ny,nx/),            & 
         regDecomp=(/NtileJ,NtileI/),   &
         coordSys=ESMF_COORDSYS_CART,   &
         coordDep1=(/1/), & ! 1st coord is 1D and depends on 1st Grid dim
         coordDep2=(/2/), & ! 2nd coord is 1D and depends on 2nd Grid dim
         indexflag=ESMF_INDEX_GLOBAL,   &
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

    NULLIFY(yCoords)
    NULLIFY(xCoords)
    
    rc = ESMF_SUCCESS

  END SUBROUTINE Create_ISOMIP_plus_grid
    
END MODULE FISOC_ISM_Wrapper
