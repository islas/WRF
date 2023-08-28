MODULE module_scalar_tables
  USE module_driver_constants
  USE module_state_description
  USE module_domain_type, ONLY : streamrec
#include "scalar_tables.inc"
CONTAINS
  SUBROUTINE init_module_scalar_tables
     INTEGER i , j
     DO j = 1, max_domains
#include "scalar_tables_init.inc"
     END DO
  END SUBROUTINE init_module_scalar_tables
END MODULE module_scalar_tables