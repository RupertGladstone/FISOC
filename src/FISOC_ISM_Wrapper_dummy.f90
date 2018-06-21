
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
  END INTERFACE FISOC_ISM_Wrapper_Init_Phase1

CONTAINS


  SUBROUTINE FISOC_ISM_Wrapper_Init_Phase1_grid(FISOC_config,vm,ISM_ExpFB,ISM_Grid,rc)

    TYPE(ESMF_config),INTENT(INOUT)       :: FISOC_config
    TYPE(ESMF_VM),INTENT(INOUT)           :: vm
    TYPE(ESMF_fieldBundle),INTENT(INOUT)  :: ISM_ExpFB
    TYPE(ESMF_grid),INTENT(OUT)           :: ISM_Grid
    INTEGER,INTENT(OUT),OPTIONAL          :: rc

    msg = "ERROR: Dummy subroutine called probably due to ISM_gridType error"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
         line=__LINE__, file=__FILE__, rc=rc)
    CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

  END SUBROUTINE FISOC_ISM_Wrapper_Init_Phase1_grid



  !--------------------------------------------------------------------------------------
  ! This dummy wrapper aims to create the dummy mesh and required variables 
  ! in the ESMF formats.  
  SUBROUTINE FISOC_ISM_Wrapper_Init_Phase1_mesh(FISOC_config,vm,ISM_ExpFB,ISM_dummyMesh,rc)

    TYPE(ESMF_config),INTENT(INOUT)       :: FISOC_config
    TYPE(ESMF_VM),INTENT(INOUT)           :: vm
    TYPE(ESMF_fieldBundle),INTENT(INOUT)  :: ISM_ExpFB
    TYPE(ESMF_mesh),INTENT(OUT)           :: ISM_dummyMesh
    INTEGER,INTENT(OUT),OPTIONAL          :: rc

    CHARACTER(len=ESMF_MAXSTR)            :: label
    CHARACTER(len=ESMF_MAXSTR),ALLOCATABLE:: ISM_ReqVarList(:)
    INTEGER                               :: localPet
    LOGICAL                               :: verbose_coupling

    rc = ESMF_FAILURE

    CALL ESMF_ConfigGetAttribute(FISOC_config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    label = 'FISOC_ISM_ReqVars:'
    CALL FISOC_getListFromConfig(FISOC_config, label, ISM_ReqVarList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_VMGet(vm, localPet=localPet, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       PRINT*,""
       PRINT*,"******************************************************************************"
       PRINT*,"**********    ISM dummy wrapper.  Init phase 1 method.    ********************"
       PRINT*,"******************************************************************************"
       PRINT*,""
       PRINT*,"Here we need to get the ISM mesh information into the ESMF_grid type. "
       PRINT*,"We also need to create and initialise the required variables using the "
       PRINT*,"ESMF_field type and put them into an ESMF_fieldBundle type."
       PRINT*,""
    END IF

    CALL dummyCreateMesh(ISM_dummyMesh)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL FISOC_populateFieldBundle(ISM_ReqVarList,ISM_ExpFB,ISM_dummyMesh,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_ISM_Wrapper_Init_Phase1_mesh
  

  !--------------------------------------------------------------------------------------
  ! This dummy wrapper aims to create the dummy mesh and required variables 
  ! in the ESMF formats.  
  SUBROUTINE FISOC_ISM_Wrapper_Init_Phase2(FISOC_config,vm,ISM_ImpFB,ISM_ExpFB,rc)

    TYPE(ESMF_fieldBundle),INTENT(INOUT),OPTIONAL  :: ISM_ImpFB, ISM_ExpFB
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
  SUBROUTINE FISOC_ISM_Wrapper_Run(FISOC_config,vm,ISM_ImpFB,ISM_ExpFB,rc)

    TYPE(ESMF_config),INTENT(INOUT)               :: FISOC_config
    TYPE(ESMF_fieldbundle),INTENT(INOUT),OPTIONAL :: ISM_ImpFB,ISM_ExpFB
    TYPE(ESMF_VM),INTENT(IN)                      :: vm
    INTEGER,INTENT(OUT),OPTIONAL                  :: rc

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
       PRINT*,"See the Elmer wrapper for examples of accessing ESMF vars and field bundles."
    END IF

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_ISM_Wrapper_Run


  !--------------------------------------------------------------------------------------
  SUBROUTINE FISOC_ISM_Wrapper_Finalize(FISOC_config,vm,rc)

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
       PRINT*,"************    ISM dummy wrapper.  Finalise method.     *********************"
       PRINT*,"******************************************************************************"
       PRINT*,""
       PRINT*,"FISOC has taken care of clearing up ESMF types.  Here we just need to call the "
       PRINT*,"ISM finalise method."
    END IF

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_ISM_Wrapper_Finalize


  !--------------------------------------------------------------------------------------
  SUBROUTINE dummyCreateMesh(ISM_mesh)
    TYPE(ESMF_mesh),INTENT(INOUT) :: ISM_mesh
    INTEGER,ALLOCATABLE :: nodeOwners(:), nodeIds(:),elemIds(:), elemTypes(:), elemConn(:)
    REAL(ESMF_KIND_R8),ALLOCATABLE :: nodeCoords(:)
    INTEGER :: numNodes, numQuadElems, numTriElems, numTotElems, rc
    ! Note that currently only meshes on spherical coords can be read in from file. So
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
    nodeCoords=(/0.0,0.0, & ! node id 1
         1000000.0,0.0, & ! node id 2
         2000000.0,0.0, & ! node id 3
         0.0,20000.0, & ! node id 4
         1000000.0,20000.0, & ! node id 5
         2000000.0,20000.0, & ! node id 6
         0.0,40000.0, & ! node id 7
         1000000.0,40000.0, & ! node id 8
         2000000.0,40000.0 /) ! node id 9
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
         ESMF_MESHELEMTYPE_TRI, & ! elem id 2
         ESMF_MESHELEMTYPE_TRI, & ! elem id 3
         ESMF_MESHELEMTYPE_QUAD, & ! elem id 4
         ESMF_MESHELEMTYPE_QUAD/) ! elem id 5
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
    elemConn=(/1,2,5,4, & ! elem id 1
         2,3,5, & ! elem id 2
         3,6,5, & ! elem id 3
         4,5,8,7, & ! elem id 4
         5,6,9,8/) ! elem id 5
    ! Create Mesh structure in 1 step
    ISM_mesh = ESMF_MeshCreate(parametricDim=2,spatialDim=2, &
         nodeIds=nodeIds, nodeCoords=nodeCoords, &
         nodeOwners=nodeOwners, elementIds=elemIds,&
         elementTypes=elemTypes, elementConn=elemConn, &
         coordSys=ESMF_COORDSYS_CART,                        &
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
    ! field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, rc=localrc)
  END SUBROUTINE dummyCreateMesh
    
END MODULE FISOC_ISM_Wrapper
