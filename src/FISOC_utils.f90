
MODULE FISOC_utils

  USE ESMF

  IMPLICIT NONE

  PRIVATE

  PUBLIC  FISOC_getStringListFromConfig

  CHARACTER(len=ESMF_MAXSTR) :: msg

CONTAINS

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
    CALL ESMF_ConfigFindLabel(config, 'FISOC_ISM_ReqVars:',rc=rc)
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
