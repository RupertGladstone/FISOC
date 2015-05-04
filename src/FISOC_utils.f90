
MODULE FISOC_utils_MOD

  USE ESMF

  IMPLICIT NONE

  PRIVATE

  PUBLIC  FISOC_getStringListFromConfig, FISOC_populateFieldBundle, FISOC_ConfigDerivedAttribute, &
       FISOC_initCumulatorFB, FISOC_zeroBundle, FISOC_cumulateFB, FISOC_processCumulator, msg

  INTERFACE FISOC_populateFieldBundle
      MODULE PROCEDURE FISOC_populateFieldBundleOn2dGrid
      MODULE PROCEDURE FISOC_populateFieldBundleOnMesh
  END INTERFACE FISOC_populateFieldBundle

  INTERFACE FISOC_ConfigDerivedAttribute
      MODULE PROCEDURE FISOC_ConfigDerivedAttributeInteger
      MODULE PROCEDURE FISOC_ConfigDerivedAttributeStaggerLocArray
  END INTERFACE 

  CHARACTER(len=ESMF_MAXSTR) :: msg

CONTAINS

  
  !--------------------------------------------------------------------------------------
  ! set values of fields in this bundle to zero
  SUBROUTINE FISOC_zeroBundle(fieldBundle,rc)
    
    TYPE(ESMF_fieldbundle),INTENT(INOUT)  :: fieldBundle
    INTEGER,OPTIONAL,INTENT(OUT)          :: rc

    TYPE(ESMF_Field),ALLOCATABLE          :: fieldList(:)
    REAL(ESMF_KIND_R8),POINTER            :: field_ptr(:,:)
    INTEGER                               :: fieldCount, ii

    rc = ESMF_FAILURE

    ! How many fields in bundle?
    CALL ESMF_FieldBundleGet(fieldBundle, fieldCount=fieldCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! ... get list of fields from bundle.
    ALLOCATE(fieldList(fieldCount))
    CALL ESMF_FieldBundleGet(fieldBundle, fieldList=fieldList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! loop over fields, setting all values to zero
    DO ii=1,fieldCount
       CALL ESMF_FieldGet(field=fieldList(ii), localDe=0, farrayPtr=field_ptr, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       field_ptr(:,:) = 0.0

    END DO

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_zeroBundle

  !--------------------------------------------------------------------------------------
  ! add the values for the field in the bundle to the cumulator.  Later it will be divided 
  ! by the number of cumulation operations to give the average. 
  SUBROUTINE FISOC_cumulateFB(fieldBundle,FBcumulator,rc)

    TYPE(ESMF_fieldbundle),INTENT(IN)     :: fieldBundle
    TYPE(ESMF_fieldbundle),INTENT(INOUT)  :: FBcumulator
    INTEGER,OPTIONAL,INTENT(OUT)          :: rc

    INTEGER                               :: fieldCount, ii, fieldCountCum
    REAL(ESMF_KIND_R8)                    :: init_value
    TYPE(ESMF_Field),ALLOCATABLE          :: fieldList(:),fieldListCum(:)
    TYPE(ESMF_TypeKind_Flag)              :: fieldTypeKind
    CHARACTER(len=ESMF_MAXSTR)            :: fieldName, fieldNameCum
    REAL(ESMF_KIND_R8),POINTER            :: fieldCum_ptr(:,:), field_ptr(:,:) 
    TYPE(ESMF_GRID)                       :: grid

    rc = ESMF_FAILURE

    ! How many fields in bundle?
    CALL ESMF_FieldBundleGet(FBcumulator, fieldCount=fieldCountCum, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! ... get list of fields from bundle.
    ALLOCATE(fieldListCum(fieldCountCum))
    CALL ESMF_FieldBundleGet(FBcumulator, fieldList=fieldListCum, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! How many fields in bundle?
    CALL ESMF_FieldBundleGet(fieldBundle, fieldCount=fieldCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! ... get list of fields from bundle.
    ALLOCATE(fieldList(fieldCount))
    CALL ESMF_FieldBundleGet(fieldBundle, fieldList=fieldList, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF (fieldCountCum .NE. fieldCount) THEN
       msg = 'ERROR: fieldBundle and cumulator fieldBundle contain different number of fields '
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    DO ii=1,fieldCount

       CALL ESMF_FieldGet(fieldList(ii), name=fieldName, typekind=fieldTypeKind, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       CALL ESMF_FieldGet(fieldListCum(ii), name=fieldNameCum, typekind=fieldTypeKind, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       IF (TRIM(fieldNameCum) .NE. TRIM(fieldName)//"_cum") THEN
          msg = 'ERROR: cumulator field name does not match field name '
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
               line=__LINE__, file=__FILE__, rc=rc)
          CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       END IF
    
       CALL ESMF_FieldGet(field=fieldList(ii), localDe=0, farrayPtr=field_ptr, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       CALL ESMF_FieldGet(field=fieldListCum(ii), localDe=0, farrayPtr=fieldCum_ptr, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       IF (SIZE(field_ptr).NE.SIZE(fieldCum_ptr)) THEN
          msg = 'ERROR: cumulator field array pointer dims do not match '
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
               line=__LINE__, file=__FILE__, rc=rc)
          CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       END IF
    
       fieldCum_ptr(:,:) = fieldCum_ptr(:,:) + field_ptr(:,:)

    END DO

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_cumulateFB


  !--------------------------------------------------------------------------------------
  ! Calculate the average from the cumulator, and stick it back in the field bundle. 
  ! We may one day wish to do other time processing such as keeping the total intead of 
  ! averaging (cant think why though...).
  ! Note: there is a lot of code duplication here of FISOC_cumulateFB.  Is there potential 
  ! for some kind of operatore overloading for field bundle operations?  Perhaps this 
  ! already exists in ESMF?
  SUBROUTINE FISOC_processCumulator(fieldBundle,FBcumulator,FISOC_config,rc)

    TYPE(ESMF_fieldbundle),INTENT(INOUT)  :: FBcumulator, fieldBundle
    TYPE(ESMF_config),INTENT(INOUT)       :: FISOC_config
    INTEGER,OPTIONAL,INTENT(OUT)          :: rc

!    TYPE(ESMF_field)                      :: cumulatorField
    CHARACTER(len=ESMF_MAXSTR)            :: fieldName, fieldNameCum
    INTEGER                               :: OM_cum_steps, ii, fieldCount, fieldCountCum 
    TYPE(ESMF_field)                      :: fieldCum, field
    REAL(ESMF_KIND_R8),POINTER            :: fieldCum_ptr(:,:), field_ptr(:,:) 
    TYPE(ESMF_Field),ALLOCATABLE          :: fieldList(:),fieldListCum(:)
    TYPE(ESMF_TypeKind_Flag)              :: fieldTypeKind

    rc = ESMF_FAILURE

    ! get number of cumulator steps: ts_ratio / OM_outputInterval
    CALL FISOC_ConfigDerivedAttributeInteger(FISOC_config, OM_cum_steps, 'OM_cum_steps', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! How many fields in bundle?
    CALL ESMF_FieldBundleGet(FBcumulator, fieldCount=fieldCountCum, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! ... get list of fields from bundle.
    ALLOCATE(fieldListCum(fieldCountCum))
    CALL ESMF_FieldBundleGet(FBcumulator, fieldList=fieldListCum, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! How many fields in bundle?
    CALL ESMF_FieldBundleGet(fieldBundle, fieldCount=fieldCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! ... get list of fields from bundle.
    ALLOCATE(fieldList(fieldCount))
    CALL ESMF_FieldBundleGet(fieldBundle, fieldList=fieldList, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF (fieldCountCum .NE. fieldCount) THEN
       msg = 'ERROR: fieldBundle and cumulator fieldBundle contain different number of fields '
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    DO ii=1,fieldCount

       CALL ESMF_FieldGet(fieldList(ii), name=fieldName, typekind=fieldTypeKind, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       CALL ESMF_FieldGet(fieldListCum(ii), name=fieldNameCum, typekind=fieldTypeKind, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       IF (TRIM(fieldNameCum) .NE. TRIM(fieldName)//"_cum") THEN
          msg = 'ERROR: cumulator field name does not match field name '
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
               line=__LINE__, file=__FILE__, rc=rc)
          CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       END IF
    
       CALL ESMF_FieldGet(field=fieldList(ii), localDe=0, farrayPtr=field_ptr, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       CALL ESMF_FieldGet(field=fieldListCum(ii), localDe=0, farrayPtr=fieldCum_ptr, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       IF (SIZE(field_ptr).NE.SIZE(fieldCum_ptr)) THEN
          msg = 'ERROR: cumulator field array pointer dims do not match '
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
               line=__LINE__, file=__FILE__, rc=rc)
          CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       END IF
    
       field_ptr(:,:) = fieldCum_ptr(:,:) / OM_cum_steps

    END DO

    CALL FISOC_zeroBundle(FBcumulator,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_processCumulator


  !--------------------------------------------------------------------------------------
  ! populate a field bundle for cumulating outputs using the field info from an 
  ! existing field bundle.
  SUBROUTINE FISOC_initCumulatorFB(fieldBundle,FBcumulator,rc)

    TYPE(ESMF_fieldbundle),INTENT(IN)     :: fieldBundle
    TYPE(ESMF_fieldbundle),INTENT(INOUT)  :: FBcumulator
    INTEGER,OPTIONAL,INTENT(OUT)          :: rc

    INTEGER                               :: fieldCount, ii
    REAL(ESMF_KIND_R8)                    :: init_value
    TYPE(ESMF_Field),ALLOCATABLE          :: fieldList(:),fieldListCum(:)
    TYPE(ESMF_TypeKind_Flag)              :: fieldTypeKind
    CHARACTER(len=ESMF_MAXSTR)            :: fieldName, fieldNameCum
    TYPE(ESMF_GRID)                       :: grid

    rc = ESMF_FAILURE

    init_value = 0.0

    ! How many fields in bundle?
    CALL ESMF_FieldBundleGet(fieldBundle, fieldCount=fieldCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! ... get list of fields from bundle.
    ALLOCATE(fieldList(fieldCount))
    CALL ESMF_FieldBundleGet(fieldBundle, fieldList=fieldList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! loop over fields, creating a cumulator for each to add to the cumulator field bundle
    ! (cumulator is just another field in which we will sum the field over time)
    ALLOCATE(fieldListCum(fieldCount))
    DO ii=1,fieldCount

       CALL ESMF_FieldGet(fieldList(ii), name=fieldName, typekind=fieldTypeKind, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       fieldNameCum = TRIM(fieldName)//"_cum"

       CALL ESMF_FieldGet(fieldList(ii), grid=grid, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       
       fieldListCum(ii) = ESMF_FieldCreate(grid, typekind=fieldTypeKind, name=fieldNameCum, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    END DO

    CALL ESMF_FieldBundleAdd(FBcumulator, fieldListCum, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_initCumulatorFB


  !--------------------------------------------------------------------------------------
  SUBROUTINE FISOC_ConfigDerivedAttributeInteger(FISOC_config, derivedAttribute, label,rc)
    
    CHARACTER(len=*),INTENT(IN)           :: label
    TYPE(ESMF_config),INTENT(INOUT)       :: FISOC_config
    INTEGER,INTENT(OUT)                   :: derivedAttribute
    INTEGER,OPTIONAL,INTENT(OUT)          :: rc
    
    INTEGER                               :: OM_dt_sec, dt_ratio, OM_outputInterval

    rc = ESMF_FAILURE

    SELECT CASE(label)

    CASE('ISM_dt_sec')
       CALL ESMF_ConfigGetAttribute(FISOC_config, dt_ratio, label='dt_ratio:', rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__, rcToReturn=rc)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       CALL ESMF_ConfigGetAttribute(FISOC_config, OM_dt_sec, label='OM_dt_sec:', rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__, rcToReturn=rc)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       derivedAttribute = OM_dt_sec * dt_ratio
    CASE('OM_cum_steps')
       CALL ESMF_ConfigGetAttribute(FISOC_config, dt_ratio, label='dt_ratio:', rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__, rcToReturn=rc)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       CALL ESMF_ConfigGetAttribute(FISOC_config, OM_outputInterval, label='OM_outputInterval:', rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__, rcToReturn=rc)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       derivedAttribute = dt_ratio / OM_outputInterval
    CASE DEFAULT
       msg = 'ERROR: unrecognised derived config attribute label '
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END SELECT
    
    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_ConfigDerivedAttributeInteger


  !--------------------------------------------------------------------------------------
  SUBROUTINE FISOC_ConfigDerivedAttributeStaggerLocArray(FISOC_config, derivedAttribute, label, rc)
    
    CHARACTER(len=*),INTENT(IN)           :: label
    TYPE(ESMF_config),INTENT(INOUT)       :: FISOC_config
    TYPE(ESMF_staggerLoc),INTENT(INOUT)   :: derivedAttribute(:)
    INTEGER,OPTIONAL,INTENT(OUT)          :: rc
    
    CHARACTER(len=ESMF_MAXSTR),ALLOCATABLE:: attribute_stringList(:)
    INTEGER                               :: OM_dt_sec, dt_ratio, OM_outputInterval

    rc = ESMF_FAILURE

    SELECT CASE(label)

    CASE('OM_ReqVars_stagger:')
       CALL FISOC_getStringListFromConfig(FISOC_config, label, attribute_stringList,rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       CALL FISOC_OM_staggerCodes(derivedAttribute,attribute_stringList,rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CASE DEFAULT
       msg = 'ERROR: unrecognised derived config attribute label '
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END SELECT
    
    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_ConfigDerivedAttributeStaggerLocArray


  !--------------------------------------------------------------------------------------
  ! convert a list of strings describning stagger location to ESMF stagger location 
  ! integer codes (for use in ESMF operations such as creating fields on a grid).
  SUBROUTINE FISOC_OM_staggerCodes(staggerLoc,staggerChar,rc)
    
    CHARACTER(len=ESMF_MAXSTR),ALLOCATABLE,INTENT(IN) :: staggerChar(:)
    TYPE(ESMF_staggerLoc),INTENT(OUT)                 :: staggerLoc(:)
    INTEGER,INTENT(OUT),OPTIONAL                      :: rc

    INTEGER                                           :: ii

    rc = ESMF_FAILURE
    
    IF (size(staggerLoc) .ne. size(staggerChar)) THEN
       msg = 'ERROR: stagger lists must be same length'
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF
       
    DO ii = 1,size(staggerLoc)
       SELECT CASE(staggerChar(ii))
       CASE("EDGE1")
          staggerLoc(ii) = ESMF_STAGGERLOC_EDGE1
       CASE("EDGE2")
          staggerLoc(ii) = ESMF_STAGGERLOC_EDGE2
       CASE("CENTER")
          staggerLoc(ii) = ESMF_STAGGERLOC_CENTER
       CASE("CORNER")
          staggerLoc(ii) = ESMF_STAGGERLOC_CORNER
       CASE DEFAULT
          msg = 'ERROR: unrecognised staggerLoc string '
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
               line=__LINE__, file=__FILE__, rc=rc)
          CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       END SELECT
    END DO

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_OM_staggerCodes


  !--------------------------------------------------------------------------------------
  SUBROUTINE FISOC_populateFieldBundleOn2dGrid(fieldNames,fieldBundle,grid,init_value,fieldStagger,rc)

    CHARACTER(len=ESMF_MAXSTR),INTENT(IN)    :: fieldNames(:)
    TYPE(ESMF_grid),INTENT(IN)               :: grid
    REAL(ESMF_KIND_R8),INTENT(IN),OPTIONAL   :: init_value
    TYPE(ESMF_staggerLoc),INTENT(IN),OPTIONAL:: fieldStagger(:)

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

END MODULE FISOC_utils_MOD
