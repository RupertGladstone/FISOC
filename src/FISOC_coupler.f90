MODULE FISOC_coupler
  
  USE ESMF
    
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC FISOC_coupler_register
    
CONTAINS
  
  SUBROUTINE FISOC_coupler_register(FISOC_coupler, rc)
    
    TYPE(ESMF_CplComp)  :: FISOC_coupler
    INTEGER, INTENT(OUT) :: rc
    
    rc = ESMF_FAILURE

    rc = ESMF_SUCCESS
 
  END SUBROUTINE FISOC_coupler_register

END MODULE FISOC_coupler
