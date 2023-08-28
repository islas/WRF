  interface 
    ! Split the allocation bits into separate routines.
    !
    ! The following bit of maximal foolishness brought to you by the makers of certain compilers
    ! that can't process large files of simple branches, assignment statements and calls to very simple
    ! procedures and functions without whining to their mommies about it being too complex and insisting 
    ! you call your service representative or wait geological time to produce code that old 1-pass 
    ! compilers (or compilers for other languages) could have managed efficiently in a second or two.
    !

    module subroutine ROUTINENAME ( grid,   id, setinitval_in ,  tl_in , inter_domain_in , okay_to_alloc_in, num_bytes_allocated , &
                                    sd31, ed31, sd32, ed32, sd33, ed33, &
                                    sm31 , em31 , sm32 , em32 , sm33 , em33 , &
                                    sp31 , ep31 , sp32 , ep32 , sp33 , ep33 , &
                                    sp31x, ep31x, sp32x, ep32x, sp33x, ep33x, &
                                    sp31y, ep31y, sp32y, ep32y, sp33y, ep33y, &
                                    sm31x, em31x, sm32x, em32x, sm33x, em33x, &
                                    sm31y, em31y, sm32y, em32y, sm33y, em33y )
        USE module_domain_type
        USE module_configure, ONLY : model_config_rec, grid_config_rec_type, in_use_for_config, model_to_grid_config_rec
  !      USE module_state_description
        USE module_scalar_tables ! this includes module_state_description too

        IMPLICIT NONE

        !  Input data.

        TYPE(domain)               , POINTER          :: grid
        INTEGER , INTENT(IN)            :: id
        INTEGER , INTENT(IN)            :: setinitval_in   ! 3 = everything, 1 = arrays only, 0 = none
        INTEGER , INTENT(IN)            :: sd31, ed31, sd32, ed32, sd33, ed33
        INTEGER , INTENT(IN)            :: sm31, em31, sm32, em32, sm33, em33
        INTEGER , INTENT(IN)            :: sp31, ep31, sp32, ep32, sp33, ep33
        INTEGER , INTENT(IN)            :: sp31x, ep31x, sp32x, ep32x, sp33x, ep33x
        INTEGER , INTENT(IN)            :: sp31y, ep31y, sp32y, ep32y, sp33y, ep33y
        INTEGER , INTENT(IN)            :: sm31x, em31x, sm32x, em32x, sm33x, em33x
        INTEGER , INTENT(IN)            :: sm31y, em31y, sm32y, em32y, sm33y, em33y

        ! this argument is a bitmask. First bit is time level 1, second is time level 2, and so on.
        ! e.g. to set both 1st and second time level, use 3
        !      to set only 1st                        use 1
        !      to set only 2st                        use 2
        INTEGER , INTENT(IN)            :: tl_in
  
        ! true if the allocation is for an intermediate domain (for nesting); only certain fields allocated
        ! false otherwise (all allocated, modulo tl above)
        LOGICAL , INTENT(IN)            :: inter_domain_in, okay_to_alloc_in

        INTEGER(KIND=8) , INTENT(INOUT)         :: num_bytes_allocated
    end subroutine ROUTINENAME
  end interface