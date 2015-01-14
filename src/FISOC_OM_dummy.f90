MODULE FISOC_OM
  
  USE ESMF
    
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC FISOC_OM_register
    
  CHARACTER(len=ESMF_MAXSTR) :: msg

CONTAINS
  
  SUBROUTINE FISOC_OM_register(FISOC_OM, rc)
    
    TYPE(ESMF_GridComp)  :: FISOC_OM
    INTEGER, INTENT(OUT) :: rc
    
    rc = ESMF_FAILURE

    CALL ESMF_GridCompSetEntryPoint(FISOC_OM, ESMF_METHOD_INITIALIZE, &
         userRoutine=FISOC_OM_init_phase1, phase=1, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_GridCompSetEntryPoint(FISOC_OM, ESMF_METHOD_INITIALIZE, &
         userRoutine=FISOC_OM_init_phase2, phase=2, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_GridCompSetEntryPoint(FISOC_OM, ESMF_METHOD_RUN, &
         userRoutine=FISOC_OM_run, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_GridCompSetEntryPoint(FISOC_OM, ESMF_METHOD_FINALIZE, &
         userRoutine=FISOC_OM_finalise, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    rc = ESMF_SUCCESS
 
  END SUBROUTINE FISOC_OM_register


  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_OM_init_phase1(FISOC_OM, OM_ImpSt, OM_ExpSt, FISOC_clock, rc)
    TYPE(ESMF_GridComp)    :: FISOC_OM
    TYPE(ESMF_State)       :: OM_ImpSt, OM_ExpSt 
    TYPE(ESMF_Clock)       :: FISOC_clock
    INTEGER, INTENT(OUT)   :: rc

    REAL(ESMF_KIND_R8),POINTER :: coordY(:),coordX(:)
    TYPE(ESMF_grid)        :: OM_grid
    INTEGER                :: ii, jj, lbnd(1), ubnd(1)
    TYPE(ESMF_field)       :: OM_dBdt_l0
    TYPE(ESMF_fieldBundle) :: OM_FB
    REAL(ESMF_KIND_R8),POINTER :: OM_dBdt_l0_ptr(:,:)

    rc = ESMF_FAILURE

    NULLIFY (coordY,coordX,OM_dBdt_l0_ptr)

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
    OM_grid=ESMF_GridCreateNoPeriDim(          &
         ! Define a regular distribution
         maxIndex=(/11,6/), & ! define index space
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
    CALL ESMF_GridAddCoord(OM_grid,  & 
         staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
    
    !-------------------------------------------------------------------
    ! Get the pointer to the first coordinate array and the bounds
    ! of its global indices on the local DE.   
    !-------------------------------------------------------------------
    CALL ESMF_GridGetCoord(OM_grid, coordDim=1, localDE=0, &
         staggerloc=ESMF_STAGGERLOC_CENTER, &
         computationalLBound=lbnd, computationalUBound=ubnd, &
         farrayPtr=coordX, rc=rc)
    
    !-------------------------------------------------------------------
    ! Calculate and set coordinates in the first dimension.
    !-------------------------------------------------------------------
    DO ii=lbnd(1),ubnd(1)
       coordX(ii) = (ii-1)*180000.0
    END DO
    
    !-------------------------------------------------------------------
    ! Get the pointer to the second coordinate array and the bounds of
    ! its global indices on the local DE.
    !-------------------------------------------------------------------
    CALL ESMF_GridGetCoord(OM_grid, coordDim=2, localDE=0, &
         staggerloc=ESMF_STAGGERLOC_CENTER, &
         computationalLBound=lbnd, computationalUBound=ubnd, &
         farrayPtr=coordY, rc=rc)
    
    !-------------------------------------------------------------------
    ! Calculate and set coordinates in the second dimension 
    !-------------------------------------------------------------------
    DO jj=lbnd(1),ubnd(1)
       coordY(jj) = (jj-1)*10000.0
    END DO
    
    OM_dBdt_l0 = ESMF_FieldCreate(OM_grid, typekind=ESMF_TYPEKIND_R8, name="OM_dBdt_l0", rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=OM_dBdt_l0, localDe=0, farrayPtr=OM_dBdt_l0_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    OM_dBdt_l0_ptr(:,:) = 1.0

    OM_FB = ESMF_FieldBundleCreate(name="OM fields", rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_FieldBundleAdd(OM_FB, (/OM_dBdt_l0/), rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_StateAdd(OM_ImpSt, (/OM_FB/), rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_StateAdd(OM_ExpSt, (/OM_FB/), rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    msg = "OM created grid and field and added field to import and export states"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_OM_init_phase1
  
  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_OM_init_phase2(FISOC_OM, OM_ImpSt, OM_ExpSt, FISOC_clock, rc)
    TYPE(ESMF_GridComp)    :: FISOC_OM
    TYPE(ESMF_State)       :: OM_ImpSt, OM_ExpSt 
    TYPE(ESMF_Clock)       :: FISOC_clock
    INTEGER, INTENT(OUT)   :: rc

    rc = ESMF_FAILURE

    msg = "OM initialise phase 2 allows the OM access to the ISM initial state "// &
         "(does nothing for dummy case)"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_OM_init_phase2
  
  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_OM_run(FISOC_OM, OM_ImpSt, OM_ExpSt, FISOC_clock, rc)
    TYPE(ESMF_GridComp)    :: FISOC_OM
    TYPE(ESMF_State)       :: OM_ImpSt, OM_ExpSt 
    TYPE(ESMF_Clock)       :: FISOC_clock
    INTEGER, INTENT(OUT)   :: rc

    rc = ESMF_FAILURE

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_OM_run
  
  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_OM_finalise(FISOC_OM, OM_ImpSt, OM_ExpSt, FISOC_clock, rc)
    TYPE(ESMF_GridComp)    :: FISOC_OM
    TYPE(ESMF_State)       :: OM_ImpSt, OM_ExpSt 
    TYPE(ESMF_Clock)       :: FISOC_clock
    INTEGER, INTENT(OUT)   :: rc

    rc = ESMF_FAILURE

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_OM_finalise  

END MODULE FISOC_OM
