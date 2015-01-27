MODULE FISOC_ISM
  
  USE ESMF
  USE FISOC_ISM_Wrapper
  USE FISOC_utils
    
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC FISOC_ISM_register
  
  CHARACTER(len=ESMF_MAXSTR) :: msg

CONTAINS
  
  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_ISM_register(FISOC_ISM, rc)
    
    TYPE(ESMF_GridComp)  :: FISOC_ISM
    INTEGER, INTENT(OUT) :: rc
    
    rc = ESMF_FAILURE

!    CALL ESMF_GridCompSetEntryPoint(FISOC_ISM, ESMF_METHOD_INITIALIZE, &
!         userRoutine=FISOC_ISM_init, rc=rc)
!    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!         line=__LINE__, file=__FILE__)) RETURN
    
    CALL ESMF_GridCompSetEntryPoint(FISOC_ISM, ESMF_METHOD_INITIALIZE, &
         userRoutine=FISOC_ISM_init_phase1, phase=1, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) RETURN

    CALL ESMF_GridCompSetEntryPoint(FISOC_ISM, ESMF_METHOD_INITIALIZE, &
         userRoutine=FISOC_ISM_init_phase2, phase=2, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) RETURN
    
    CALL ESMF_GridCompSetEntryPoint(FISOC_ISM, ESMF_METHOD_RUN, &
         userRoutine=FISOC_ISM_run, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) RETURN
    
    CALL ESMF_GridCompSetEntryPoint(FISOC_ISM, ESMF_METHOD_FINALIZE, &
         userRoutine=FISOC_ISM_finalise, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) RETURN

    rc = ESMF_SUCCESS
 
  END SUBROUTINE FISOC_ISM_register

  !------------------------------------------------------------------------------
  ! Initialisation is implemented in two stages.  The first stage is for independent 
  ! initialisation of the ISM and the second stage occurrs after the OM has been 
  ! initialised, to allow inter-component consistency checks or use of the OM 
  ! state to complete ISM initalisation.
  SUBROUTINE FISOC_ISM_init_phase1(FISOC_ISM, ISM_ImpSt, ISM_ExpSt, FISOC_clock, rc)
    TYPE(ESMF_GridComp)    :: FISOC_ISM
    TYPE(ESMF_State)       :: ISM_ImpSt, ISM_ExpSt
    TYPE(ESMF_Clock)       :: FISOC_clock
    INTEGER, INTENT(OUT)   :: rc

    TYPE(ESMF_config)      :: config
    TYPE(ESMF_mesh)        :: ISM_mesh
    TYPE(ESMF_fieldBundle) :: ISM_ExpFB
    CHARACTER(len=ESMF_MAXSTR) :: ISM_meshFile

    TYPE(ESMF_field)       :: ISM_temperature_l0, ISM_temperature_l1
    TYPE(ESMF_field)       :: ISM_z_l0, ISM_z_l1
    TYPE(ESMF_field)       :: ISM_dTdz_l0, ISM_z_l0_previous
    REAL(ESMF_KIND_R8),POINTER :: ISM_temperature_l0_ptr(:),ISM_temperature_l1_ptr(:) 
    REAL(ESMF_KIND_R8),POINTER :: ISM_z_l0_ptr(:),ISM_z_l1_ptr(:) 
    REAL(ESMF_KIND_R8),POINTER :: ISM_dTdz_l0_ptr(:),ISM_z_l0_previous_ptr(:) 

    INTEGER                :: numNodes, numQuadElems, numTriElems, numTotElems
    INTEGER,ALLOCATABLE    :: nodeOwners(:), nodeIds(:),elemIds(:), elemTypes(:), elemConn(:)
    REAL(ESMF_KIND_R8),ALLOCATABLE :: nodeCoords(:) 

    CHARACTER(len=ESMF_MAXSTR),ALLOCATABLE :: ISM_ReqVarList(:)
    CHARACTER(len=ESMF_MAXSTR) :: label



!    real(ESMF_KIND_R8)        :: ownedNodeCoords(:)
!    integer                   :: numOwnedElements
!    logical                   :: isMemFreed
!    type(ESMF_CoordSys_Flag)  :: coordSys
!    integer                   :: parametricDim
!    integer                   :: spatialDim
!    type(ESMF_DistGrid)       :: nodalDistgrid
!    type(ESMF_DistGrid)       :: elementDistgrid
!    integer                   :: numOwnedNodes

!    INTEGER :: fieldcount

    rc = ESMF_SUCCESS

    NULLIFY (ISM_temperature_l0_ptr,ISM_temperature_l1_ptr,ISM_z_l0_ptr, &
         ISM_z_l1_ptr,ISM_dTdz_l0_ptr,ISM_z_l0_previous_ptr)


    msg = "ISM initialise started"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

    ! Get mesh file name from the config file
    CALL ESMF_GridCompGet(FISOC_ISM, config=config, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return

    label = 'FISOC_ISM_ReqVars:'
    CALL FISOC_getStringListFromConfig(config, label, ISM_ReqVarList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(config, ISM_meshFile, label='ISM_meshFile:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__, rcToReturn=rc)) return

    CALL FISOC_ISM_Wrapper_Init(ISM_ExpFB,ISM_mesh,config)


    ! Note that currently only meshes on spherical coords can be read in from file.  So 
    ! here we use the subroutine interfaces to create simple mesh instead.
  
    ! Set number of nodes
    numNodes=9

    ! Allocate and fill the node id array.
    allocate(nodeIds(numNodes))
    nodeIds=(/1,2,3,4,5,6,7,8,9/) 

    ! Allocate and fill node coordinate array.
    ! Since this is a 2D Mesh the size is 2x the
    ! number of nodes.
    allocate(nodeCoords(2*numNodes))
    nodeCoords=(/0.0,0.0,    & ! node id 1
         1000000.0,0.0,      & ! node id 2
         2000000.0,0.0,      & ! node id 3
         0.0,20000.0,        & ! node id 4
         1000000.0,20000.0,  & ! node id 5
         2000000.0,20000.0,  & ! node id 6
         0.0,40000.0,        & ! node id 7
         1000000.0,40000.0,  & ! node id 8
         2000000.0,40000.0 /)  ! node id 9
    ! if node coords go from 0 to 2000000 and from 0 to 40000 then domain is 2000km by 40km    

    ! Allocate and fill the node owner array.
    ! Since this Mesh is all on PET 0, it's just set to all 0.
    allocate(nodeOwners(numNodes))
    nodeOwners=0 ! everything on PET 0

    ! Set the number of each type of element, plus the total number.
    numQuadElems=3
    numTriElems=2
    numTotElems=numQuadElems+numTriElems

    ! Allocate and fill the element id array.
    allocate(elemIds(numTotElems))
    elemIds=(/1,2,3,4,5/) 

    ! Allocate and fill the element topology type array.
    allocate(elemTypes(numTotElems))
    elemTypes=(/ESMF_MESHELEMTYPE_QUAD, & ! elem id 1
         ESMF_MESHELEMTYPE_TRI,  & ! elem id 2
         ESMF_MESHELEMTYPE_TRI,  & ! elem id 3
         ESMF_MESHELEMTYPE_QUAD, & ! elem id 4
         ESMF_MESHELEMTYPE_QUAD/)  ! elem id 5


    ! Allocate and fill the element connection type array.
    ! Note that entries in this array refer to the 
    ! positions in the nodeIds, etc. arrays and that
    ! the order and number of entries for each element
    ! reflects that given in the Mesh options 
    ! section for the corresponding entry
    ! in the elemTypes array. The number of 
    ! entries in this elemConn array is the
    ! number of nodes in a quad. (4) times the 
    ! number of quad. elements plus the number
    ! of nodes in a triangle (3) times the number
    ! of triangle elements. 
    allocate(elemConn(4*numQuadElems+3*numTriElems))
    elemConn=(/1,2,5,4, &  ! elem id 1
         2,3,5,   &  ! elem id 2
         3,6,5,   &  ! elem id 3
         4,5,8,7, &  ! elem id 4
         5,6,9,8/)   ! elem id 5

    ! Create Mesh structure in 1 step
    ISM_mesh = ESMF_MeshCreate(parametricDim=2,spatialDim=2, &
         nodeIds=nodeIds, nodeCoords=nodeCoords, &
         nodeOwners=nodeOwners, elementIds=elemIds,&
         elementTypes=elemTypes, elementConn=elemConn, &
         rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    
    ! After the creation we are through with the arrays, so they may be
    ! deallocated.
    deallocate(nodeIds)
    deallocate(nodeCoords)
    deallocate(nodeOwners)
    deallocate(elemIds)
    deallocate(elemTypes)
    deallocate(elemConn)
    
    ! At this point the mesh is ready to use. For example, as is 
    ! illustrated here, to have a field created on it. Note that 
    ! the Field only contains data for nodes owned by the current PET.
    ! Please see Section "Create a Field from a Mesh" under Field
    ! for more information on creating a Field on a Mesh. 
!    field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8,  rc=localrc)

! no longer creating mesh from file... currently ESMF will only do spherical coord meshes from file
!    ! create ISM mesh and use it to create zeroed fields for the ISM import and export states
!    ISM_mesh = ESMF_MeshCreate(filename=ISM_meshFile, &
!            filetypeflag=ESMF_FILEFORMAT_ESMFMESH, &
!            rc=rc)
!    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!         line=__LINE__, file=__FILE__)) &
!         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

! Note weakness: currently regridding in 2d instead of a 2d manifold in 3d space.

    ISM_temperature_l0 = ESMF_FieldCreate(ISM_mesh, typekind=ESMF_TYPEKIND_R8, name="ISM_temperature_l0", rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_temperature_l0, localDe=0, farrayPtr=ISM_temperature_l0_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    ISM_temperature_l0_ptr(:) = -1.0
    ISM_temperature_l0_ptr(1) = -10.0

    ISM_temperature_l1 = ESMF_FieldCreate(ISM_mesh, typekind=ESMF_TYPEKIND_R8, name="ISM_temperature_l1", rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_temperature_l1, localDe=0, farrayPtr=ISM_temperature_l1_ptr, rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    ISM_temperature_l1_ptr(:) = -1.1
    
    ISM_z_l0 = ESMF_FieldCreate(ISM_mesh, typekind=ESMF_TYPEKIND_R8, name="ISM_z_l0", rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_z_l0, localDe=0, farrayPtr=ISM_z_l0_ptr, rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    ISM_z_l0_ptr(:) = -100.0
    
    ISM_z_l1 = ESMF_FieldCreate(ISM_mesh, typekind=ESMF_TYPEKIND_R8, name="ISM_z_l1", rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_z_l1, localDe=0, farrayPtr=ISM_z_l1_ptr, rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    ISM_z_l1_ptr(:) = -75.0
    
    ISM_dTdz_l0 = ESMF_FieldCreate(ISM_mesh, typekind=ESMF_TYPEKIND_R8, name="ISM_dTdz_l0", rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_dTdz_l0, localDe=0, farrayPtr=ISM_dTdz_l0_ptr, rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    ISM_dTdz_l0_ptr(:) = 0.0
    
    ISM_z_l0_previous = ESMF_FieldCreate(ISM_mesh, typekind=ESMF_TYPEKIND_R8, name="ISM_z_l0_previous", rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_z_l0_previous, localDe=0, farrayPtr=ISM_z_l0_previous_ptr, rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    ISM_z_l0_previous_ptr(:) = ISM_z_l0_ptr(:)
    
    msg = "ISM created mesh and fields"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

    ISM_ExpFB = ESMF_FieldBundleCreate(name="ISM export fields", rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_FieldBundleAdd(ISM_ExpFB, (/ISM_temperature_l0/), rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_FieldBundleAdd(ISM_ExpFB, (/ISM_temperature_l1/), rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_FieldBundleAdd(ISM_ExpFB, (/ISM_z_l0/), rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_FieldBundleAdd(ISM_ExpFB, (/ISM_z_l1/), rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_FieldBundleAdd(ISM_ExpFB, (/ISM_dTdz_l0/), rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_FieldBundleAdd(ISM_ExpFB, (/ISM_z_l0_previous/), rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL FISOC_ISM_calcDerivedFields(ISM_ExpFB,config,rc)

    ! we only add the fields to the import state as a way of letting the coupler get hold of the 
    ! grid.  There must be a better way to do this.
    CALL ESMF_StateAdd(ISM_ImpSt, (/ISM_ExpFB/), rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_StateAdd(ISM_ExpSt, (/ISM_ExpFB/), rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    msg = "ISM bundled fields and added to import and export states"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

  END SUBROUTINE FISOC_ISM_init_phase1


  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_ISM_init_phase2(FISOC_ISM, ISM_ImpSt, ISM_ExpSt, FISOC_clock, rc)
    TYPE(ESMF_GridComp)    :: FISOC_ISM
    TYPE(ESMF_State)       :: ISM_ImpSt, ISM_ExpSt
    TYPE(ESMF_Clock)       :: FISOC_clock
    INTEGER, INTENT(OUT)   :: rc

    rc = ESMF_FAILURE

    msg = "ISM initialise phase 2 allows the ISM access to the OM initial state "// &
         "(does nothing for dummy case)"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_ISM_init_phase2
  

  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_ISM_run(FISOC_ISM, ISM_ImpSt, ISM_ExpSt, FISOC_clock, rc)
    TYPE(ESMF_GridComp)  :: FISOC_ISM
    TYPE(ESMF_State)     :: ISM_ImpSt, ISM_ExpSt
    TYPE(ESMF_Clock)     :: FISOC_clock
    INTEGER              :: petCount
    INTEGER, INTENT(OUT) :: rc
    
    LOGICAL                :: verbose_coupling
    TYPE(ESMF_config)      :: config
    TYPE(ESMF_fieldbundle) :: ISM_ImpFB,ISM_ExpFB
    TYPE(ESMF_field)       :: OM_dBdt_l0, ISM_z_l0, ISM_z_l1
    TYPE(ESMF_fieldBundle) :: OM_FB
    REAL(ESMF_KIND_R8),POINTER :: OM_dBdt_l0_ptr(:),ISM_z_l0_ptr(:),ISM_z_l1_ptr(:)


    rc = ESMF_FAILURE

    CALL ESMF_GridCompGet(FISOC_ISM, config=config, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)


    CALL ESMF_StateGet(ISM_ImpSt, "ISM import fields", ISM_ImpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    
    
    CALL ESMF_FieldBundleGet(ISM_ImpFB, fieldName="OM_dBdt_l0", field=OM_dBdt_l0, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_FieldGet(field=OM_dBdt_l0, localDe=0, farrayPtr=OM_dBdt_l0_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF (verbose_coupling) THEN
       PRINT*,"ISM run phase. Lets just adjust depth coords according to basal melt rate from ocean."
       PRINT*,"Melt rates dont look quite right on the ISM mesh, needs checking..."
       PRINT*,OM_dBdt_l0_ptr
    END IF

    CALL ESMF_StateGet(ISM_ExpSt, "ISM export fields", ISM_ExpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    
    
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

    msg = "ISM run complete for current timestep"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)
    
    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_ISM_run


  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_ISM_finalise(FISOC_ISM, ISM_ImpSt, ISM_ExpSt, FISOC_clock, rc)
    TYPE(ESMF_GridComp)  :: FISOC_ISM
    TYPE(ESMF_State)     :: ISM_ImpSt, ISM_ExpSt
    TYPE(ESMF_Clock)     :: FISOC_clock
    INTEGER              :: petCount
    INTEGER, INTENT(OUT) :: rc

    TYPE(ESMF_fieldbundle)       :: ISM_ImpFB, ISM_ExpFB
    INTEGER                      :: FieldCount,ii
    TYPE(ESMF_field),ALLOCATABLE :: FieldList(:)
    TYPE(ESMF_mesh)              :: ISM_mesh
    
    rc = ESMF_FAILURE

    msg = "ISM finalise: destroy fields, bundles and mesh"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)


    CALL ESMF_StateGet(ISM_ImpSt, "ISM import fields", ISM_ImpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    

    ! ...how many fields?...
    CALL ESMF_FieldBundleGet(ISM_ImpFB, fieldCount=FieldCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! ... get list of fields from bundle.
    ALLOCATE(FieldList(FieldCount))
    CALL ESMF_FieldBundleGet(ISM_ImpFB, fieldCount=FieldCount, fieldList=FieldList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_FieldGet(FieldList(1), mesh=ISM_mesh, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_MeshDestroy(ISM_mesh,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    DO ii = 1,FieldCount
       CALL ESMF_FieldDestroy(FieldList(ii), rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END DO

    DEALLOCATE(FieldList)

    CALL ESMF_FieldBundleDestroy(ISM_ImpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    



    CALL ESMF_StateGet(ISM_ExpSt, "ISM export fields", ISM_ExpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    

    ! ...how many fields?...
    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldCount=FieldCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! ... get list of fields from bundle.
    ALLOCATE(FieldList(FieldCount))
    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldCount=FieldCount, fieldList=FieldList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    DO ii = 1,FieldCount
       CALL ESMF_FieldDestroy(FieldList(ii), rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END DO

    DEALLOCATE(FieldList)

    CALL ESMF_FieldBundleDestroy(ISM_ExpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_ISM_finalise
  

  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_ISM_calcDerivedFields(ISM_ExpFB,config,rc)

    TYPE(ESMF_config),INTENT(IN) :: config
    INTEGER, INTENT(OUT)         :: rc
    TYPE(ESMF_fieldBundle)       :: ISM_ExpFB

    TYPE(ESMF_field)       :: ISM_temperature_l0, ISM_temperature_l1
    TYPE(ESMF_field)       :: ISM_z_l0, ISM_z_l1
    TYPE(ESMF_field)       :: ISM_dTdz_l0, ISM_z_l0_previous
    REAL(ESMF_KIND_R8),POINTER :: ISM_temperature_l0_ptr(:),ISM_temperature_l1_ptr(:) 
    REAL(ESMF_KIND_R8),POINTER :: ISM_z_l0_ptr(:),ISM_z_l1_ptr(:) 
    REAL(ESMF_KIND_R8),POINTER :: ISM_dTdz_l0_ptr(:),ISM_z_l0_previous_ptr(:) 

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_temperature_l0", field=ISM_temperature_l0, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_temperature_l0, localDe=0, farrayPtr=ISM_temperature_l0_ptr, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_temperature_l1", field=ISM_temperature_l1, rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_temperature_l1, localDe=0, farrayPtr=ISM_temperature_l1_ptr, rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_z_l0", field=ISM_z_l0, rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_z_l0, localDe=0, farrayPtr=ISM_z_l0_ptr, rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_z_l1", field=ISM_z_l1, rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_z_l1, localDe=0, farrayPtr=ISM_z_l1_ptr, rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_FieldBundleGet(ISM_ExpFB, fieldName="ISM_dTdz_l0", field=ISM_dTdz_l0, rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_FieldGet(field=ISM_dTdz_l0, localDe=0, farrayPtr=ISM_dTdz_l0_ptr, rc=rc)
    IF (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! the simplest approximation for the vertical temperature gradient at the ice base is to 
    ! use the lowest two levels and assume the gradient doesn't change between these two 
    ! levels.
    ISM_dTdz_l0_ptr(:) = (ISM_temperature_l1_ptr-ISM_temperature_l0_ptr) / (ISM_z_l1_ptr-ISM_z_l0_ptr)

!    PRINT*,"temperature gradient: add configurable options for more advanced calculations"

!    print*,"previous B to be calculated only at run time"

  END SUBROUTINE FISOC_ISM_calcDerivedFields

END MODULE FISOC_ISM
