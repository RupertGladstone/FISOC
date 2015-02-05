
MODULE FISOC_utils

  USE ESMF

  IMPLICIT NONE

  PRIVATE

  PUBLIC  FISOC_getStringListFromConfig, FISOC_populateFieldBundle

  INTERFACE FISOC_populateFieldBundle
      MODULE PROCEDURE FISOC_populateFieldBundleOn2dGrid
      MODULE PROCEDURE FISOC_populateFieldBundleOnMesh
  END INTERFACE FISOC_populateFieldBundle

  CHARACTER(len=ESMF_MAXSTR) :: msg

CONTAINS


  !--------------------------------------------------------------------------------------
  SUBROUTINE FISOC_populateFieldBundleOn2dGrid(fieldNames,fieldBundle,grid,init_value,rc)

    CHARACTER(len=ESMF_MAXSTR),INTENT(IN) :: fieldNames(:)
    TYPE(ESMF_grid),INTENT(IN)            :: grid
    REAL(ESMF_KIND_R8),INTENT(IN),OPTIONAL:: init_value

    TYPE(ESMF_fieldbundle),INTENT(INOUT)  :: fieldBundle
    INTEGER,INTENT(OUT),OPTIONAL          :: rc

    INTEGER                               :: ii
    REAL(ESMF_KIND_R8)                    :: initial_value
    TYPE(ESMF_field)                      :: field
    REAL(ESMF_KIND_R8),POINTER            :: field_ptr(:,:) 

    rc = ESMF_FAILURE

    NULLIFY(field_ptr)

    IF (PRESENT(init_value)) THEN
       initial_value = init_value
    ELSE
       initial_value = 0.0
    END IF

    DO ii=1,SIZE(fieldNames)
       field = ESMF_FieldCreate(grid, typekind=ESMF_TYPEKIND_R8, name=TRIM(fieldNames(ii)), rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       CALL ESMF_FieldGet(field=field, localDe=0, farrayPtr=field_ptr, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       field_ptr(:,:) = initial_value       
       CALL ESMF_FieldBundleAdd(fieldBundle, (/field/), rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END DO

    NULLIFY(field_ptr)

    rc = ESMF_SUCCESS
    
  END SUBROUTINE FISOC_populateFieldBundleOn2dGrid
  

  !--------------------------------------------------------------------------------------
  ! Note: this subroutine is almost identical to the grid version, and should probably 
  ! be auto generated in a precompile step rather than the current hard coded duplication.
  SUBROUTINE FISOC_populateFieldBundleOnMesh(fieldNames,fieldBundle,mesh,init_value,rc)

    CHARACTER(len=ESMF_MAXSTR),INTENT(IN) :: fieldNames(:)
    TYPE(ESMF_mesh),INTENT(IN)            :: mesh
    REAL(ESMF_KIND_R8),INTENT(IN),OPTIONAL:: init_value

    TYPE(ESMF_fieldbundle),INTENT(INOUT)  :: fieldBundle
    INTEGER,INTENT(OUT),OPTIONAL          :: rc

    INTEGER                               :: ii
    REAL(ESMF_KIND_R8)                    :: initial_value
    TYPE(ESMF_field)                      :: field
    REAL(ESMF_KIND_R8),POINTER            :: field_ptr(:) 

    rc = ESMF_FAILURE

    NULLIFY(field_ptr)

    IF (PRESENT(init_value)) THEN
       initial_value = init_value
    ELSE
       initial_value = 0.0
    END IF

    DO ii=1,SIZE(fieldNames)
       field = ESMF_FieldCreate(mesh, typekind=ESMF_TYPEKIND_R8, name=TRIM(fieldNames(ii)), rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       CALL ESMF_FieldGet(field=field, localDe=0, farrayPtr=field_ptr, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       field_ptr(:) = initial_value       
       CALL ESMF_FieldBundleAdd(fieldBundle, (/field/), rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END DO

    NULLIFY(field_ptr)

    rc = ESMF_SUCCESS
    
  END SUBROUTINE FISOC_populateFieldBundleOnMesh
  

  !--------------------------------------------------------------------------------------
  SUBROUTINE FISOC_getStringListFromConfig(config,label,stringList,rc)

    CHARACTER(len=ESMF_MAXSTR),ALLOCATABLE,INTENT(INOUT) :: stringList(:)

    TYPE(ESMF_config),INTENT(INOUT)      :: config

    CHARACTER(len=ESMF_MAXSTR),INTENT(IN):: label
    INTEGER,INTENT(OUT),OPTIONAL         :: rc

    CHARACTER(len=ESMF_MAXSTR)           :: dummyString
    INTEGER                              :: listCount,ii

    rc = ESMF_FAILURE

    ! point config to start of list 
    CALL ESMF_ConfigFindLabel(config,TRIM(label),rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! how many items in list?
    listCount = 0
    DO WHILE (rc.EQ.0)
       CALL ESMF_ConfigGetAttribute(config, dummyString,rc=rc) 
       IF  (rc.EQ.ESMF_RC_NOT_FOUND) EXIT
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       listCount = listCount + 1
    END DO
    CALL ESMF_ConfigFindLabel(config, TRIM(label),rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! retrieve list items
    ALLOCATE(stringList(listCount))
    DO ii = 1,listCount
       CALL ESMF_ConfigGetAttribute(config, stringList(ii),rc=rc) 
       IF  (rc.EQ.ESMF_RC_NOT_FOUND) EXIT
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END DO

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_getStringListFromConfig

END MODULE FISOC_utils
