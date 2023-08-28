submodule( module_configure ) INIT
implicit none
  contains
    module SUBROUTINE initial_config
    !<DESCRIPTION>
    ! This routine reads in the namelist.input file and sets
    ! module_config_rec, a structure of TYPE(model_config_rec_type), which is is seen via USE association by any
    ! subprogram that uses module_configure.  The module_config_rec structure
    ! contains all namelist settings for all domains.  Variables that apply
    ! to the entire run and have only one value regardless of domain are
    ! scalars.  Variables that allow different settings for each domain are
    ! defined as arrays of dimension max_domains (defined in
    ! frame/module_driver_constants.F, from a setting passed in from
    ! configure.wrf). There is another type in WRF, TYPE(grid_config_rec_type), in which
    ! all fields pertain only to a single domain (and are all scalars). The subroutine
    ! model_to_grid_config_rec(), also in frame/module_configure.F, is used to retrieve
    ! the settings for a given domain from a TYPE(module_config_rec_type) and put them into
    ! a TYPE(grid_config_rec_type), variables of which type are often called <em>config_flags</em>
    ! in the WRF code.
    !
    ! Most of the code in this routine is generated from the Registry file
    ! rconfig entries and included from the following files (found in the inc directory):
    !
    ! <pre>
    ! namelist_defines.inc	declarations of namelist variables (local to this routine)
    ! namelist_statements.inc	NAMELIST statements for each variable
    ! namelist_defaults.inc	assignment to default values if specified in Registry
    ! config_reads.inc		read statements for each namelist record
    ! config_assigns.inc	assign each variable to field in module_config_rec
    ! </pre>
    !
    !NOTE: generated subroutines from Registry file rconfig entries are renamed nl_
    ! instead of rconfig_ due to length limits for subroutine names.
    !
    ! Note for version WRF 2.0: there is code here to force all domains to
    ! have the same mp_physics setting. This is because different mp_physics
    ! packages have different numbers of tracers but the nest forcing and
    ! feedback code relies on the parent and nest having the same number and
    ! kind of tracers. This means that the microphysics option
    ! specified on the highest numbered domain is the microphysics
    ! option for <em>all</em> domains in the run. This will be revisited.
    !
    !</DESCRIPTION>
          IMPLICIT NONE

          INTEGER              :: io_status
          INTEGER              :: i

          LOGICAL              :: nml_read_error

          CHARACTER (LEN=1024) :: nml_name

          INTEGER, PARAMETER :: nml_write_unit= 9
          INTEGER, PARAMETER :: nml_read_unit = 10

          CHARACTER (LEN=1024) :: entire_line


    ! define as temporaries
#include "namelist_defines.inc"

    ! Statements that specify the namelists
#include "namelist_statements.inc"

          OPEN ( UNIT   = nml_read_unit    ,      &
                FILE   = "namelist.input" ,      &
                FORM   = "FORMATTED"      ,      &
                STATUS = "OLD"            ,      &
                IOSTAT = io_status         )

          IF ( io_status .NE. 0 ) THEN
            CALL WRF_ERROR_FATAL ( 'ERROR OPENING namelist.input' )
          ENDIF

#ifndef NO_NAMELIST_PRINT
          OPEN ( UNIT   = nml_write_unit    ,      &
#if (DA_CORE == 1)
                FILE   = "namelist.output.da" ,      &
#else
                FILE   = "namelist.output" ,      &
#endif
                FORM   = "FORMATTED"      ,      &
                STATUS = "REPLACE"        ,      &
                IOSTAT = io_status         )

          IF ( io_status .NE. 0 ) THEN
#if (DA_CORE == 1)
            CALL WRF_ERROR_FATAL ( 'ERROR OPENING namelist.output.da' )
#else
            CALL WRF_ERROR_FATAL ( 'ERROR OPENING namelist.output' )
#endif
          ENDIF
#endif

    ! Statements that set the namelist vars to default vals
#  include "namelist_defaults.inc"

#if (DA_CORE == 1)
    ! Override the default values, because we can not assigned a arrary with different values in registry.

          as1(1:3) = (/ 0.25, 1.0, 1.5 /)
          as2(1:3) = (/ 0.25, 1.0, 1.5 /)
          as3(1:3) = (/ 0.25, 1.0, 1.5 /)
          as4(1:3) = (/ 0.25, 1.0, 1.5 /)
          as5(1:3) = (/ 0.25, 1.0, 1.5 /)
#endif

    ! Statements that read the namelist are in this file
#  include "config_reads.inc"

    ! 2004/04/28  JM (with consensus by the group of developers)
    ! This is needed to ensure that nesting will work, since
    ! different mp_physics packages have different numbers of
    ! tracers. Basically, this says that the microphysics option
    ! specified on the highest numbered domain *is* the microphysics
    ! option for the run. Not the best solution but okay for 2.0.
    !

          DO i = 1, max_dom
            mp_physics(i) = mp_physics(max_dom)
          ENDDO

    ! Statements that assign the variables to the cfg record are in this file
    ! except the namelist_derived variables where are assigned below
#undef SOURCE_RECORD
#undef DEST_RECORD
#undef SOURCE_REC_DEX
#define SOURCE_RECORD
#define DEST_RECORD model_config_rec %
#define SOURCE_REC_DEX
#include "config_assigns.inc"


          CLOSE ( UNIT = nml_read_unit , IOSTAT = io_status )

          IF ( io_status .NE. 0 ) THEN
            CALL WRF_ERROR_FATAL ( 'ERROR CLOSING namelist.input' )
          ENDIF

#ifndef NO_NAMELIST_PRINT
          CLOSE ( UNIT = nml_write_unit , IOSTAT = io_status )

          IF ( io_status .NE. 0 ) THEN
            CALL WRF_ERROR_FATAL ( 'ERROR CLOSING namelist.output' )
          ENDIF
#endif

#ifdef _WIN32
          model_config_rec%nocolons = .TRUE.   ! always no colons for Windows
#endif

          RETURN

    END SUBROUTINE initial_config
end submodule INIT