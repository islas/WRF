
   SUBROUTINE push_communicators_for_domain( id )
      USE module_dm
      INTEGER, INTENT(IN) :: id   ! if specified also does an instate for grid id
!  Only required for distrbuted memory parallel runs
#if ( defined( DM_PARALLEL ) && ( ! defined( STUBMPI ) ) )
      IF ( communicator_stack_cursor .GE. max_domains ) CALL wrf_error_fatal("push_communicators_for_domain would excede stacksize")
      communicator_stack_cursor = communicator_stack_cursor + 1

      id_stack(communicator_stack_cursor) = current_id
      local_communicator_stack( communicator_stack_cursor )    =    local_communicator
      local_communicator_periodic_stack( communicator_stack_cursor )  =    local_communicator_periodic
      local_iocommunicator_stack( communicator_stack_cursor )  =    local_iocommunicator
      local_communicator_x_stack( communicator_stack_cursor )  =    local_communicator_x
      local_communicator_y_stack( communicator_stack_cursor )  =    local_communicator_y
      ntasks_stack( communicator_stack_cursor )        =    ntasks
      ntasks_y_stack( communicator_stack_cursor )      =    ntasks_y
      ntasks_x_stack( communicator_stack_cursor )      =    ntasks_x
      mytask_stack( communicator_stack_cursor )        =    mytask
      mytask_x_stack( communicator_stack_cursor )       =    mytask_x
      mytask_y_stack( communicator_stack_cursor )       =    mytask_y

      CALL instate_communicators_for_domain( id )
#endif
   END SUBROUTINE push_communicators_for_domain
   SUBROUTINE pop_communicators_for_domain
      USE module_dm
      IMPLICIT NONE
      !  Only required for distrbuted memory parallel runs
#if ( defined( DM_PARALLEL ) && ( ! defined( STUBMPI ) ) )
      IF ( communicator_stack_cursor .LT. 1 ) CALL wrf_error_fatal("pop_communicators_for_domain on empty stack")
      current_id = id_stack(communicator_stack_cursor)
      local_communicator = local_communicator_stack( communicator_stack_cursor )
      local_communicator_periodic = local_communicator_periodic_stack( communicator_stack_cursor )
      local_iocommunicator = local_iocommunicator_stack( communicator_stack_cursor )
      local_communicator_x = local_communicator_x_stack( communicator_stack_cursor )
      local_communicator_y = local_communicator_y_stack( communicator_stack_cursor )
      ntasks = ntasks_stack( communicator_stack_cursor )
      ntasks_y = ntasks_y_stack( communicator_stack_cursor )
      ntasks_x = ntasks_x_stack( communicator_stack_cursor )
      mytask = mytask_stack( communicator_stack_cursor )
      mytask_x = mytask_x_stack( communicator_stack_cursor )
      mytask_y = mytask_y_stack( communicator_stack_cursor )
      communicator_stack_cursor = communicator_stack_cursor - 1
#endif
   END SUBROUTINE pop_communicators_for_domain
   SUBROUTINE instate_communicators_for_domain( id )
      USE module_dm
!  Only required for distrbuted memory parallel runs
#if ( defined( DM_PARALLEL ) && ( ! defined( STUBMPI ) ) )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: id
      INTEGER ierr
      current_id = id
      local_communicator          = local_communicator_store( id )
      local_communicator_periodic = local_communicator_periodic_store( id )
      local_iocommunicator        = local_iocommunicator_store( id )
      local_communicator_x        = local_communicator_x_store( id )
      local_communicator_y        = local_communicator_y_store( id )
      ntasks         = ntasks_store( id )
      mytask         = mytask_store( id )
      ntasks_x       = ntasks_x_store( id )
      ntasks_y       = ntasks_y_store( id )
      mytask_x       = mytask_x_store( id )
      mytask_y       = mytask_y_store( id )
#endif
   END SUBROUTINE instate_communicators_for_domain
   SUBROUTINE store_communicators_for_domain( id )
      USE module_dm
!  Only required for distrbuted memory parallel runs
#if ( defined( DM_PARALLEL ) && ( ! defined( STUBMPI ) ) )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: id
      local_communicator_store( id )    =    local_communicator
      local_communicator_periodic_store( id )  =    local_communicator_periodic
      local_iocommunicator_store( id )  =    local_iocommunicator
      local_communicator_x_store( id )  =    local_communicator_x
      local_communicator_y_store( id )  =    local_communicator_y
      ntasks_store( id )        =    ntasks
      ntasks_x_store( id )      =    ntasks_x
      ntasks_y_store( id )      =    ntasks_y
      mytask_store( id )        =    mytask
      mytask_x_store( id )      =    mytask_x
      mytask_y_store( id )      =    mytask_y
#endif
   END SUBROUTINE store_communicators_for_domain

!=========================================================================
! wrf_dm_patch_domain has to be outside the module because it is called
! by a routine in module_domain but depends on module domain

SUBROUTINE wrf_dm_patch_domain ( id  , domdesc , parent_id , parent_domdesc , &
                          sd1 , ed1 , sp1 , ep1 , sm1 , em1 , &
                          sd2 , ed2 , sp2 , ep2 , sm2 , em2 , &
                          sd3 , ed3 , sp3 , ep3 , sm3 , em3 , &
                                      sp1x , ep1x , sm1x , em1x , &
                                      sp2x , ep2x , sm2x , em2x , &
                                      sp3x , ep3x , sm3x , em3x , &
                                      sp1y , ep1y , sm1y , em1y , &
                                      sp2y , ep2y , sm2y , em2y , &
                                      sp3y , ep3y , sm3y , em3y , &
                          bdx , bdy )
   USE module_domain, ONLY : domain, head_grid, find_grid_by_id
   USE module_dm, ONLY : patch_domain_rsl_lite  !, push_communicators_for_domain, pop_communicators_for_domain
   IMPLICIT NONE

   INTEGER, INTENT(IN)   :: sd1 , ed1 , sd2 , ed2 , sd3 , ed3 , bdx , bdy
   INTEGER, INTENT(OUT)  :: sp1 , ep1 , sp2 , ep2 , sp3 , ep3 , &
                            sm1 , em1 , sm2 , em2 , sm3 , em3
   INTEGER               :: sp1x , ep1x , sp2x , ep2x , sp3x , ep3x , &
                            sm1x , em1x , sm2x , em2x , sm3x , em3x
   INTEGER               :: sp1y , ep1y , sp2y , ep2y , sp3y , ep3y , &
                            sm1y , em1y , sm2y , em2y , sm3y , em3y
   INTEGER, INTENT(INOUT):: id  , domdesc , parent_id , parent_domdesc

   TYPE(domain), POINTER :: parent
   TYPE(domain), POINTER :: grid_ptr

   ! this is necessary because we cannot pass parent directly into
   ! wrf_dm_patch_domain because creating the correct interface definitions
   ! would generate a circular USE reference between module_domain and module_dm
   ! see comment this date in module_domain for more information. JM 20020416

   NULLIFY( parent )
   grid_ptr => head_grid
   CALL find_grid_by_id( parent_id , grid_ptr , parent )

   CALL push_communicators_for_domain(id)

   CALL patch_domain_rsl_lite ( id  , parent, parent_id , &
                           sd1 , ed1 , sp1 , ep1 , sm1 , em1 , &
                           sd2 , ed2 , sp2 , ep2 , sm2 , em2 , &
                           sd3 , ed3 , sp3 , ep3 , sm3 , em3 , &
                                      sp1x , ep1x , sm1x , em1x , &
                                      sp2x , ep2x , sm2x , em2x , &
                                      sp3x , ep3x , sm3x , em3x , &
                                      sp1y , ep1y , sm1y , em1y , &
                                      sp2y , ep2y , sm2y , em2y , &
                                      sp3y , ep3y , sm3y , em3y , &
                           bdx , bdy )

   CALL pop_communicators_for_domain

   RETURN
END SUBROUTINE wrf_dm_patch_domain

SUBROUTINE wrf_termio_dup( comm )
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: comm
  INTEGER mytask, ntasks
#ifndef STUBMPI
  INTEGER ierr
  INCLUDE 'mpif.h'
  CALL mpi_comm_size(comm, ntasks, ierr )
  CALL mpi_comm_rank(comm, mytask, ierr )
  write(0,*)'starting wrf task ',mytask,' of ',ntasks
  CALL rsl_error_dup1( mytask, ntasks )
#else
  mytask = 0
  ntasks = 1
#endif
END SUBROUTINE wrf_termio_dup

SUBROUTINE wrf_get_myproc( myproc )
  USE module_dm , ONLY : mytask
  IMPLICIT NONE
  INTEGER myproc
  myproc = mytask
  RETURN
END SUBROUTINE wrf_get_myproc

SUBROUTINE wrf_get_nproc( nproc )
  USE module_dm , ONLY : ntasks
  IMPLICIT NONE
  INTEGER nproc
  nproc = ntasks
  RETURN
END SUBROUTINE wrf_get_nproc

SUBROUTINE wrf_get_nprocx( nprocx )
  USE module_dm , ONLY : ntasks_x
  IMPLICIT NONE
  INTEGER nprocx
  nprocx = ntasks_x
  RETURN
END SUBROUTINE wrf_get_nprocx

SUBROUTINE wrf_get_nprocy( nprocy )
  USE module_dm , ONLY : ntasks_y
  IMPLICIT NONE
  INTEGER nprocy
  nprocy = ntasks_y
  RETURN
END SUBROUTINE wrf_get_nprocy

SUBROUTINE wrf_dm_bcast_bytes ( buf , size )
   USE module_dm , ONLY : local_communicator
   IMPLICIT NONE
#ifndef STUBMPI
   INCLUDE 'mpif.h'
#endif
   INTEGER size
#ifndef NEC
   INTEGER*1 BUF(size)
#else
   CHARACTER*1 BUF(size)
#endif
#ifndef STUBMPI
   CALL BYTE_BCAST ( buf , size, local_communicator )
#endif
   RETURN
END SUBROUTINE wrf_dm_bcast_bytes

SUBROUTINE wrf_dm_bcast_string( BUF, N1 )
   IMPLICIT NONE
   INTEGER n1
! <DESCRIPTION>
! Collective operation. Given a string and a size in characters on task zero, broadcast and return that buffer on all tasks.
!
! </DESCRIPTION>
   CHARACTER*(*) buf
#ifndef STUBMPI
   INTEGER ibuf(256),i,n
   CHARACTER*256 tstr
   n = n1
   ! Root task is required to have the correct value of N1, other tasks
   ! might not have the correct value.  
   CALL wrf_dm_bcast_integer( n , 1 )
   IF (n .GT. 256) n = 256
   IF (n .GT. 0 ) then
     DO i = 1, n
       ibuf(I) = ichar(buf(I:I))
     END DO
     CALL wrf_dm_bcast_integer( ibuf, n )
     buf = ''
     DO i = 1, n
       buf(i:i) = char(ibuf(i))
     END DO
   END IF
#endif
   RETURN
END SUBROUTINE wrf_dm_bcast_string

SUBROUTINE wrf_dm_bcast_string_comm( BUF, N1, COMM )
   IMPLICIT NONE
   INTEGER n1
   INTEGER COMM
! <DESCRIPTION>
! Collective operation. Given a string and a size in characters on task zero, broadcast and return that buffer on all tasks.
!
! </DESCRIPTION>
   CHARACTER*(*) buf
#ifndef STUBMPI
   INTEGER ibuf(256),i,n
   CHARACTER*256 tstr
   n = n1
   ! Root task is required to have the correct value of N1, other tasks
   ! might not have the correct value.
   CALL BYTE_BCAST( n, IWORDSIZE, COMM )
   IF (n .GT. 256) n = 256
   IF (n .GT. 0 ) then
     DO i = 1, n
       ibuf(I) = ichar(buf(I:I))
     END DO
     CALL BYTE_BCAST( ibuf, N*IWORDSIZE, COMM )
     buf = ''
     DO i = 1, n
       buf(i:i) = char(ibuf(i))
     END DO
   END IF
#endif
   RETURN
END SUBROUTINE wrf_dm_bcast_string_comm

SUBROUTINE wrf_dm_bcast_integer( BUF, N1 )
   IMPLICIT NONE
   INTEGER n1
   INTEGER  buf(*)
   CALL wrf_dm_bcast_bytes ( BUF , N1 * IWORDSIZE )
   RETURN
END SUBROUTINE wrf_dm_bcast_integer

SUBROUTINE wrf_dm_bcast_double( BUF, N1 )
   IMPLICIT NONE
   INTEGER n1
! this next declaration is REAL, not DOUBLE PRECISION because it will be autopromoted
! to double precision by the compiler when WRF is compiled for 8 byte reals. Only reason
! for having this separate routine is so we pass the correct MPI type to mpi_scatterv
! since we were not indexing the globbuf and Field arrays it does not matter
   REAL  buf(*)
   CALL wrf_dm_bcast_bytes ( BUF , N1 * DWORDSIZE )
   RETURN
END SUBROUTINE wrf_dm_bcast_double

SUBROUTINE wrf_dm_bcast_real( BUF, N1 )
   IMPLICIT NONE
   INTEGER n1
   REAL  buf(*)
   CALL wrf_dm_bcast_bytes ( BUF , N1 * RWORDSIZE )
   RETURN
END SUBROUTINE wrf_dm_bcast_real

SUBROUTINE wrf_dm_bcast_logical( BUF, N1 )
   IMPLICIT NONE
   INTEGER n1
   LOGICAL  buf(*)
   CALL wrf_dm_bcast_bytes ( BUF , N1 * LWORDSIZE )
   RETURN
END SUBROUTINE wrf_dm_bcast_logical

SUBROUTINE write_68( grid, v , s , &
                   ids, ide, jds, jde, kds, kde, &
                   ims, ime, jms, jme, kms, kme, &
                   its, ite, jts, jte, kts, kte )
  USE module_domain, ONLY : domain
  IMPLICIT NONE
  TYPE(domain) , INTENT (INOUT) :: grid
  CHARACTER *(*) s
  INTEGER ids, ide, jds, jde, kds, kde, &
          ims, ime, jms, jme, kms, kme, &
          its, ite, jts, jte, kts, kte
  REAL, DIMENSION( ims:ime , kms:kme, jms:jme ) :: v

  INTEGER i,j,k,ierr

  logical, external :: wrf_dm_on_monitor
  real globbuf( ids:ide, kds:kde, jds:jde )
  character*3 ord, stag

  if ( kds == kde ) then
    ord = 'xy'
    stag = 'xy'
  CALL wrf_patch_to_global_real ( v, globbuf, grid%domdesc, stag, ord, &
                     ids, ide, jds, jde, kds, kde, &
                     ims, ime, jms, jme, kms, kme, &
                     its, ite, jts, jte, kts, kte )
  else

    stag = 'xyz'
    ord = 'xzy'
  CALL wrf_patch_to_global_real ( v, globbuf, grid%domdesc, stag, ord, &
                     ids, ide, kds, kde, jds, jde, &
                     ims, ime, kms, kme, jms, jme, &
                     its, ite, kts, kte, jts, jte )
  endif


  if ( wrf_dm_on_monitor() ) THEN
    WRITE(68,*) ide-ids+1, jde-jds+1 , s
    DO j = jds, jde
    DO i = ids, ide
       WRITE(68,*) globbuf(i,1,j)
    END DO
    END DO
  endif

  RETURN
END

   SUBROUTINE wrf_abort

#if ( DA_CORE != 1 )
      USE module_cpl, ONLY : coupler_on, cpl_abort
#endif

      IMPLICIT NONE
#ifndef STUBMPI
      INCLUDE 'mpif.h'
      INTEGER ierr
#if ( DA_CORE != 1 )
      IF ( coupler_on ) THEN
         CALL cpl_abort( 'wrf_abort', 'look for abort message in rsl* files' )
      ELSE
#endif
         CALL mpi_abort(MPI_COMM_WORLD,1,ierr)
#if ( DA_CORE != 1 )
      END IF
#endif
#else
      STOP
#endif
   END SUBROUTINE wrf_abort

   SUBROUTINE wrf_dm_shutdown
      IMPLICIT NONE
#ifndef STUBMPI
      INTEGER ierr
      CALL MPI_FINALIZE( ierr )
#endif
      RETURN
   END SUBROUTINE wrf_dm_shutdown

   LOGICAL FUNCTION wrf_dm_on_monitor()
      IMPLICIT NONE
#ifndef STUBMPI
      INCLUDE 'mpif.h'
      INTEGER tsk, ierr, mpi_comm_local
      CALL wrf_get_dm_communicator( mpi_comm_local )
      IF ( mpi_comm_local .NE. MPI_UNDEFINED ) THEN
        CALL mpi_comm_rank ( mpi_comm_local, tsk , ierr )
        wrf_dm_on_monitor = tsk .EQ. 0
      ELSE
        wrf_dm_on_monitor = .FALSE.
      END IF
#else
      wrf_dm_on_monitor = .TRUE.
#endif
      RETURN
   END FUNCTION wrf_dm_on_monitor

   SUBROUTINE rsl_comm_iter_init(shw,ps,pe)
      INTEGER shw, ps, pe
      INTEGER iter, plus_send_start, plus_recv_start, &
                    minus_send_start, minus_recv_start
      COMMON /rcii/ iter, plus_send_start, plus_recv_start, &
                          minus_send_start, minus_recv_start
      iter = 0
      minus_send_start = ps
      minus_recv_start = ps-1
      plus_send_start = pe
      plus_recv_start = pe+1
   END SUBROUTINE rsl_comm_iter_init

   LOGICAL FUNCTION rsl_comm_iter ( id , is_intermediate,                     &
                                    shw ,  xy , ds, de_in, ps, pe, nds,nde, &
                                    sendbeg_m, sendw_m, sendbeg_p, sendw_p,   &
                                    recvbeg_m, recvw_m, recvbeg_p, recvw_p    )
      USE module_dm, ONLY : ntasks_x, ntasks_y, mytask_x, mytask_y, minx, miny, &
                            nest_pes_x, nest_pes_y
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: id,shw,xy,ds,de_in,ps,pe,nds,nde
      LOGICAL, INTENT(IN)  :: is_intermediate  ! treated differently, coarse but with same decomp as nest
      INTEGER, INTENT(OUT) :: sendbeg_m, sendw_m, sendbeg_p, sendw_p
      INTEGER, INTENT(OUT) :: recvbeg_m, recvw_m, recvbeg_p, recvw_p
      INTEGER k, kn, ni, nj, de, Px, Py, nt, ntx, nty, me, lb, ub, ierr
      INTEGER dum
      LOGICAL went
      INTEGER iter, plus_send_start, plus_recv_start, &
                    minus_send_start, minus_recv_start
      INTEGER parent_grid_ratio, parent_start
      COMMON /rcii/ iter, plus_send_start, plus_recv_start, &
                          minus_send_start, minus_recv_start

      de = de_in
      ntx = nest_pes_x(id)
      nty = nest_pes_y(id)
      IF ( xy .EQ. 1 ) THEN  ! X/I axis
        nt = ntasks_x
        me = mytask_x
        dum = 2 * nty  ! dummy dimension length for tfp to decompose without getting div 0
        IF ( is_intermediate ) THEN
           CALL nl_get_i_parent_start(id,parent_start)
           CALL nl_get_parent_grid_ratio(id,parent_grid_ratio)
        END IF
      ELSE
        nt = ntasks_y
        me = mytask_y
        dum = 2 * ntx  ! dummy dimension length for tfp to decompose without getting div 0
        IF ( is_intermediate ) THEN
           CALL nl_get_j_parent_start(id,parent_start)
           CALL nl_get_parent_grid_ratio(id,parent_grid_ratio)
        END IF
      END IF
      iter = iter + 1

#if (DA_CORE == 0)
      went = .FALSE.
      ! send to minus
      sendw_m = 0
      sendbeg_m = 1
      IF ( me .GT. 0 ) THEN
        lb = minus_send_start
        sendbeg_m = lb-ps+1
        DO k = lb,ps+shw-1
          went = .TRUE.
          IF ( xy .eq. 1 ) THEN
            IF ( is_intermediate ) THEN
              kn =  ( k - parent_start ) * parent_grid_ratio + 1 + 1 ;
              CALL task_for_point (kn,1,nds,nde,1,dum,ntx,nty,Px,Py,minx,miny,ierr) ! modified to treat x and y separately
              IF ( ierr .NE. 0 ) CALL wrf_error_fatal('error code returned by task_for_point in module_dm.F (h)')
            ELSE
              CALL task_for_point (k,1,ds,de,1,dum,ntx,nty,Px,Py,minx,miny,ierr) ! modified to treat x and y separately
              IF ( ierr .NE. 0 ) CALL wrf_error_fatal('error code returned by task_for_point in module_dm.F (i)')
            END IF
            IF ( Px .NE. me+(iter-1) ) THEN
              exit
            END IF
          ELSE
            IF ( is_intermediate ) THEN
              kn =  ( k - parent_start ) * parent_grid_ratio + 1 + 1 ;
              CALL task_for_point (1,kn,1,dum,nds,nde,ntx,nty,Px,Py,minx,miny,ierr) ! modified to treat x and y separately
              IF ( ierr .NE. 0 ) CALL wrf_error_fatal('error code returned by task_for_point in module_dm.F (h)')
            ELSE
              CALL task_for_point (1,k,1,dum,ds,de,ntx,nty,Px,Py,minx,miny,ierr) ! modified to treat x and y separately
              IF ( ierr .NE. 0 ) CALL wrf_error_fatal('error code returned by task_for_point in module_dm.F (i)')
            END IF
            IF ( Py .NE. me+(iter-1) ) THEN
              exit
            END IF
          END IF
          minus_send_start = minus_send_start+1
          sendw_m = sendw_m + 1
        END DO
      END IF
      ! recv from minus
      recvw_m = 0
      recvbeg_m = 1
      IF ( me .GT. 0 ) THEN
        ub = minus_recv_start
        recvbeg_m = ps - ub
        DO k = minus_recv_start,ps-shw,-1
          went = .TRUE.
          IF ( xy .eq. 1 ) THEN
          IF ( is_intermediate ) THEN
            kn =  ( k - parent_start ) * parent_grid_ratio + 1 + 1 ;
              CALL task_for_point (kn,1,nds,nde,1,dum,ntx,nty,Px,Py,minx,miny,ierr) ! modified to treat x and y separately
              IF ( ierr .NE. 0 ) CALL wrf_error_fatal('error code returned by task_for_point in module_dm.F (j)')
          ELSE
              CALL task_for_point (k,1,ds,de,1,dum,ntx,nty,Px,Py,minx,miny,ierr) ! modified to treat x and y separately
              IF ( ierr .NE. 0 ) CALL wrf_error_fatal('error code returned by task_for_point in module_dm.F (k)')
          END IF
          IF ( Px .NE. me-iter ) THEN
            exit
          END IF
          ELSE
            IF ( is_intermediate ) THEN
              kn =  ( k - parent_start ) * parent_grid_ratio + 1 + 1 ;
              CALL task_for_point (1,kn,1,dum,nds,nde,ntx,nty,Px,Py,minx,miny,ierr) ! modified to treat x and y separately
              IF ( ierr .NE. 0 ) CALL wrf_error_fatal('error code returned by task_for_point in module_dm.F (j)')
            ELSE
              CALL task_for_point (1,k,1,dum,ds,de,ntx,nty,Px,Py,minx,miny,ierr) ! modified to treat x and y separately
              IF ( ierr .NE. 0 ) CALL wrf_error_fatal('error code returned by task_for_point in module_dm.F (k)')
            END IF
            IF ( Py .NE. me-iter ) THEN
              exit
            END IF
          END IF
          minus_recv_start = minus_recv_start-1
          recvw_m = recvw_m + 1
        END DO
      END IF

      ! send to plus
      sendw_p = 0
      sendbeg_p = 1
      IF ( ( xy .eq. 1 .and. me .LT. ntx-1 ) .OR. ( xy .eq. 0 .and. me .LT. nty-1 ) ) THEN
        ub = plus_send_start
        sendbeg_p = pe - ub + 1
        DO k = ub,pe-shw+1,-1
          went = .TRUE.
          IF ( xy .eq. 1 ) THEN
          IF ( is_intermediate ) THEN
            kn =  ( k - parent_start ) * parent_grid_ratio + 1 + 1 ;
              CALL task_for_point (kn,1,nds,nde,1,dum,ntx,nty,Px,Py,minx,miny,ierr) ! modified to treat x and y separately
              IF ( ierr .NE. 0 ) CALL wrf_error_fatal('error code returned by task_for_point in module_dm.F (l)')
          ELSE
              CALL task_for_point (k,1,ds,de,1,dum,ntx,nty,Px,Py,minx,miny,ierr) ! modified to treat x and y separately
              IF ( ierr .NE. 0 ) CALL wrf_error_fatal('error code returned by task_for_point in module_dm.F (m)')
          END IF
          IF ( Px .NE. me-(iter-1) ) THEN
            exit
          END IF
          ELSE
            IF ( is_intermediate ) THEN
              kn =  ( k - parent_start ) * parent_grid_ratio + 1 + 1 ;
              CALL task_for_point (1,kn,1,dum,nds,nde,ntx,nty,Px,Py,minx,miny,ierr) ! modified to treat x and y separately
              IF ( ierr .NE. 0 ) CALL wrf_error_fatal('error code returned by task_for_point in module_dm.F (l)')
            ELSE
              CALL task_for_point (1,k,1,dum,ds,de,ntx,nty,Px,Py,minx,miny,ierr) ! modified to treat x and y separately
              IF ( ierr .NE. 0 ) CALL wrf_error_fatal('error code returned by task_for_point in module_dm.F (m)')
            END IF
            IF ( Py .NE. me-(iter-1) ) THEN
              exit
            END IF
          END IF
          plus_send_start = plus_send_start - 1
          sendw_p = sendw_p + 1
        END DO
      END IF
      ! recv from plus
      recvw_p = 0
      recvbeg_p = 1
      IF ( ( xy .eq. 1 .and. me .LT. ntx-1 ) .OR. ( xy .eq. 0 .and. me .LT. nty-1 ) ) THEN
        lb = plus_recv_start
        recvbeg_p = lb - pe
        DO k = lb,pe+shw
          went = .TRUE.
          IF ( xy .eq. 1 ) THEN
          IF ( is_intermediate ) THEN
            kn =  ( k - parent_start ) * parent_grid_ratio + 1 + 1 ;
              CALL task_for_point (kn,1,nds,nde,1,dum,ntx,nty,Px,Py,minx,miny,ierr) ! modified to treat x and y separately
              IF ( ierr .NE. 0 ) CALL wrf_error_fatal('error code returned by task_for_point in module_dm.F (n)')
          ELSE

              CALL task_for_point (k,1,ds,de,1,dum,ntx,nty,Px,Py,minx,miny,ierr) ! modified to treat x and y separately

              IF ( ierr .NE. 0 ) CALL wrf_error_fatal('error code returned by task_for_point in module_dm.F (o)')
          END IF
          IF ( Px .NE. me+iter ) THEN
            exit
          END IF
          ELSE
            IF ( is_intermediate ) THEN
              kn =  ( k - parent_start ) * parent_grid_ratio + 1 + 1 ;
              CALL task_for_point (1,kn,1,dum,nds,nde,ntx,nty,Px,Py,minx,miny,ierr) ! modified to treat x and y separately
              IF ( ierr .NE. 0 ) CALL wrf_error_fatal('error code returned by task_for_point in module_dm.F (n)')
            ELSE
              CALL task_for_point (1,k,1,dum,ds,de,ntx,nty,Px,Py,minx,miny,ierr) ! modified to treat x and y separately
              IF ( ierr .NE. 0 ) CALL wrf_error_fatal('error code returned by task_for_point in module_dm.F (o)')
            END IF
            IF ( Py .NE. me+iter ) THEN
              exit
            END IF
          END IF
          plus_recv_start = plus_recv_start + 1
          recvw_p = recvw_p + 1
        END DO
      END IF
#else
      if ( iter .eq. 1 ) then
        went = .true.
      else
        went = .false.
      endif
      sendw_m = 0 ; sendw_p = 0 ; recvw_m = 0 ; recvw_p = 0
      sendbeg_m = 1 ; if ( me .GT. 0 ) sendw_m = shw ;
      sendbeg_p = 1 ; if ( me .LT. nt-1 ) sendw_p = shw
      recvbeg_m = 1 ; if ( me .GT. 0 ) recvw_m = shw ;
      recvbeg_p = 1 ; if ( me .LT. nt-1 ) recvw_p = shw ;

      ! write(0,*)'shw  ', shw , ' xy ',xy
      ! write(0,*)' ds, de, ps, pe, nds,nde ',ds, de, ps, pe, nds,nde
      ! write(0,*)'sendbeg_m, sendw_m, sendbeg_p, sendw_p, recvbeg_m, recvw_m, recvbeg_p, recvw_p '
      ! write(0,*)sendbeg_m, sendw_m, sendbeg_p, sendw_p, recvbeg_m, recvw_m, recvbeg_p, recvw_p
#endif
      !if ( went ) then
      !  write(0,*)'shw  ', shw , ' xy ',xy,' plus_send_start ',plus_send_start,' minus_send_start ', minus_send_start
      !  write(0,*)' ds, de, ps, pe, nds,nde ',ds, de, ps, pe, nds,nde
      !  write(0,*)'sendbeg_m, sendw_m, sendbeg_p, sendw_p, recvbeg_m, recvw_m, recvbeg_p, recvw_p '
      !  write(0,*)sendbeg_m, sendw_m, sendbeg_p, sendw_p, recvbeg_m, recvw_m, recvbeg_p, recvw_p
      !endif
      rsl_comm_iter = went
   END FUNCTION rsl_comm_iter

   INTEGER FUNCTION wrf_dm_monitor_rank()
      IMPLICIT NONE
      wrf_dm_monitor_rank = 0
      RETURN
   END FUNCTION wrf_dm_monitor_rank

! return the global communicator if id <= 0
   SUBROUTINE wrf_get_dm_communicator_for_id ( id, communicator )
      USE module_dm , ONLY : local_communicator_store, mpi_comm_allcompute
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: id
      INTEGER , INTENT(OUT) :: communicator
      IF ( id .le. 0 ) THEN
        communicator = mpi_comm_allcompute
      ELSE
        communicator = local_communicator_store(id)
      END IF
      RETURN
   END SUBROUTINE wrf_get_dm_communicator_for_id

   SUBROUTINE wrf_get_dm_communicator ( communicator )
      USE module_dm , ONLY : local_communicator
      IMPLICIT NONE
      INTEGER , INTENT(OUT) :: communicator
      communicator = local_communicator
      RETURN
   END SUBROUTINE wrf_get_dm_communicator

   SUBROUTINE wrf_get_dm_communicator_x ( communicator )
      USE module_dm , ONLY : local_communicator_x
      IMPLICIT NONE
      INTEGER , INTENT(OUT) :: communicator
      communicator = local_communicator_x
      RETURN
   END SUBROUTINE wrf_get_dm_communicator_x

   SUBROUTINE wrf_get_dm_communicator_y ( communicator )
      USE module_dm , ONLY : local_communicator_y
      IMPLICIT NONE
      INTEGER , INTENT(OUT) :: communicator
      communicator = local_communicator_y
      RETURN
   END SUBROUTINE wrf_get_dm_communicator_y

   SUBROUTINE wrf_get_dm_iocommunicator ( iocommunicator )
      USE module_dm , ONLY : local_iocommunicator
      IMPLICIT NONE
      INTEGER , INTENT(OUT) :: iocommunicator
      iocommunicator = local_iocommunicator
      RETURN
   END SUBROUTINE wrf_get_dm_iocommunicator

   SUBROUTINE wrf_set_dm_communicator ( communicator )
      USE module_dm , ONLY : local_communicator
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: communicator
      local_communicator = communicator
      RETURN
   END SUBROUTINE wrf_set_dm_communicator

   SUBROUTINE wrf_set_dm_iocommunicator ( iocommunicator )
      USE module_dm , ONLY : local_iocommunicator
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: iocommunicator
      local_iocommunicator = iocommunicator
      RETURN
   END SUBROUTINE wrf_set_dm_iocommunicator

   SUBROUTINE wrf_get_dm_ntasks_x ( retval )
      USE module_dm , ONLY : ntasks_x
      IMPLICIT NONE
      INTEGER , INTENT(OUT) :: retval
      retval = ntasks_x
      RETURN
   END SUBROUTINE wrf_get_dm_ntasks_x

   SUBROUTINE wrf_get_dm_ntasks_y ( retval )
      USE module_dm , ONLY : ntasks_y
      IMPLICIT NONE
      INTEGER , INTENT(OUT) :: retval
      retval = ntasks_y
      RETURN
   END SUBROUTINE wrf_get_dm_ntasks_y

! added 20151212
   SUBROUTINE wrf_set_dm_quilt_comm ( communicator )
      USE module_dm , ONLY : local_quilt_comm
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: communicator
      local_quilt_comm = communicator
      RETURN
   END SUBROUTINE wrf_set_dm_quilt_comm

   SUBROUTINE wrf_get_dm_quilt_comm ( communicator )
      USE module_dm , ONLY : local_quilt_comm
      IMPLICIT NONE
      INTEGER , INTENT(OUT) :: communicator
      communicator = local_quilt_comm
      RETURN
   END SUBROUTINE wrf_get_dm_quilt_comm


!!!!!!!!!!!!!!!!!!!!!!! PATCH TO GLOBAL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE wrf_patch_to_global_real (buf,globbuf,domdesc,stagger,ordering,&
                                       DS1,DE1,DS2,DE2,DS3,DE3,&
                                       MS1,ME1,MS2,ME2,MS3,ME3,&
                                       PS1,PE1,PS2,PE2,PS3,PE3 )
       IMPLICIT NONE
       INTEGER                         DS1,DE1,DS2,DE2,DS3,DE3,&
                                       MS1,ME1,MS2,ME2,MS3,ME3,&
                                       PS1,PE1,PS2,PE2,PS3,PE3
       CHARACTER *(*) stagger,ordering
       INTEGER fid,domdesc
       REAL globbuf(*)
       REAL buf(*)

       CALL wrf_patch_to_global_generic (buf,globbuf,domdesc,stagger,ordering,RWORDSIZE,&
                                         DS1,DE1,DS2,DE2,DS3,DE3,&
                                         MS1,ME1,MS2,ME2,MS3,ME3,&
                                         PS1,PE1,PS2,PE2,PS3,PE3 )

       RETURN
   END SUBROUTINE wrf_patch_to_global_real

   SUBROUTINE wrf_patch_to_global_double (buf,globbuf,domdesc,stagger,ordering,&
                                       DS1,DE1,DS2,DE2,DS3,DE3,&
                                       MS1,ME1,MS2,ME2,MS3,ME3,&
                                       PS1,PE1,PS2,PE2,PS3,PE3 )
       IMPLICIT NONE
       INTEGER                         DS1,DE1,DS2,DE2,DS3,DE3,&
                                       MS1,ME1,MS2,ME2,MS3,ME3,&
                                       PS1,PE1,PS2,PE2,PS3,PE3
       CHARACTER *(*) stagger,ordering
       INTEGER fid,domdesc
! this next declaration is REAL, not DOUBLE PRECISION because it will be autopromoted
! to double precision by the compiler when WRF is compiled for 8 byte reals. Only reason
! for having this separate routine is so we pass the correct MPI type to mpi_scatterv
! since we were not indexing the globbuf and Field arrays it does not matter
       REAL globbuf(*)
       REAL buf(*)

       CALL wrf_patch_to_global_generic (buf,globbuf,domdesc,stagger,ordering,DWORDSIZE,&
                                         DS1,DE1,DS2,DE2,DS3,DE3,&
                                         MS1,ME1,MS2,ME2,MS3,ME3,&
                                         PS1,PE1,PS2,PE2,PS3,PE3 )

       RETURN
   END SUBROUTINE wrf_patch_to_global_double


   SUBROUTINE wrf_patch_to_global_integer (buf,globbuf,domdesc,stagger,ordering,&
                                       DS1,DE1,DS2,DE2,DS3,DE3,&
                                       MS1,ME1,MS2,ME2,MS3,ME3,&
                                       PS1,PE1,PS2,PE2,PS3,PE3 )
       IMPLICIT NONE
       INTEGER                         DS1,DE1,DS2,DE2,DS3,DE3,&
                                       MS1,ME1,MS2,ME2,MS3,ME3,&
                                       PS1,PE1,PS2,PE2,PS3,PE3
       CHARACTER *(*) stagger,ordering
       INTEGER fid,domdesc
       INTEGER globbuf(*)
       INTEGER buf(*)

       CALL wrf_patch_to_global_generic (buf,globbuf,domdesc,stagger,ordering,IWORDSIZE,&
                                         DS1,DE1,DS2,DE2,DS3,DE3,&
                                         MS1,ME1,MS2,ME2,MS3,ME3,&
                                         PS1,PE1,PS2,PE2,PS3,PE3 )

       RETURN
   END SUBROUTINE wrf_patch_to_global_integer


   SUBROUTINE wrf_patch_to_global_logical (buf,globbuf,domdesc,stagger,ordering,&
                                       DS1,DE1,DS2,DE2,DS3,DE3,&
                                       MS1,ME1,MS2,ME2,MS3,ME3,&
                                       PS1,PE1,PS2,PE2,PS3,PE3 )
       IMPLICIT NONE
       INTEGER                         DS1,DE1,DS2,DE2,DS3,DE3,&
                                       MS1,ME1,MS2,ME2,MS3,ME3,&
                                       PS1,PE1,PS2,PE2,PS3,PE3
       CHARACTER *(*) stagger,ordering
       INTEGER fid,domdesc
       LOGICAL globbuf(*)
       LOGICAL buf(*)

       CALL wrf_patch_to_global_generic (buf,globbuf,domdesc,stagger,ordering,LWORDSIZE,&
                                         DS1,DE1,DS2,DE2,DS3,DE3,&
                                         MS1,ME1,MS2,ME2,MS3,ME3,&
                                         PS1,PE1,PS2,PE2,PS3,PE3 )

       RETURN
   END SUBROUTINE wrf_patch_to_global_logical

#ifdef DEREF_KLUDGE
# define FRSTELEM (1)
#else
# define FRSTELEM
#endif

   SUBROUTINE wrf_patch_to_global_generic (buf,globbuf,domdesc,stagger,ordering,typesize,&
                                       DS1a,DE1a,DS2a,DE2a,DS3a,DE3a,&
                                       MS1a,ME1a,MS2a,ME2a,MS3a,ME3a,&
                                       PS1a,PE1a,PS2a,PE2a,PS3a,PE3a )
       USE module_driver_constants
       USE module_timing
       USE module_wrf_error, ONLY : wrf_at_debug_level
       USE module_dm, ONLY : local_communicator, ntasks

       IMPLICIT NONE
       INTEGER                         DS1a,DE1a,DS2a,DE2a,DS3a,DE3a,&
                                       MS1a,ME1a,MS2a,ME2a,MS3a,ME3a,&
                                       PS1a,PE1a,PS2a,PE2a,PS3a,PE3A
       CHARACTER *(*) stagger,ordering
       INTEGER domdesc,typesize,ierr
       REAL globbuf(*)
       REAL buf(*)
#ifndef STUBMPI
       INTEGER                         DS1,DE1,DS2,DE2,DS3,DE3,&
                                       MS1,ME1,MS2,ME2,MS3,ME3,&
                                       PS1,PE1,PS2,PE2,PS3,PE3
       INTEGER                         ids,ide,jds,jde,kds,kde,&
                                       ims,ime,jms,jme,kms,kme,&
                                       ips,ipe,jps,jpe,kps,kpe
       LOGICAL, EXTERNAL :: wrf_dm_on_monitor, has_char

       INTEGER i, j, k,  ndim
       INTEGER  Patch(3,2), Gpatch(3,2,ntasks)
    ! allocated further down, after the D indices are potentially recalculated for staggering
       REAL, ALLOCATABLE :: tmpbuf( : )
       REAL locbuf( (PE1a-PS1a+1)*(PE2a-PS2a+1)*(PE3a-PS3a+1)/RWORDSIZE*typesize+32 )

       DS1 = DS1a ; DE1 = DE1a ; DS2=DS2a ; DE2 = DE2a ; DS3 = DS3a ; DE3 = DE3a
       MS1 = MS1a ; ME1 = ME1a ; MS2=MS2a ; ME2 = ME2a ; MS3 = MS3a ; ME3 = ME3a
       PS1 = PS1a ; PE1 = PE1a ; PS2=PS2a ; PE2 = PE2a ; PS3 = PS3a ; PE3 = PE3a

       SELECT CASE ( TRIM(ordering) )
         CASE ( 'xy', 'yx' )
           ndim = 2
         CASE DEFAULT
           ndim = 3   ! where appropriate
       END SELECT

       SELECT CASE ( TRIM(ordering) )
         CASE ( 'xyz','xy' )
            ! the non-staggered variables come in at one-less than
            ! domain dimensions, but code wants full domain spec, so
            ! adjust if not staggered
           IF ( .NOT. has_char( stagger, 'x' ) ) DE1 = DE1+1
           IF ( .NOT. has_char( stagger, 'y' ) ) DE2 = DE2+1
           IF ( ndim .EQ. 3 .AND. .NOT. has_char( stagger, 'z' ) ) DE3 = DE3+1
         CASE ( 'yxz','yx' )
           IF ( .NOT. has_char( stagger, 'x' ) ) DE2 = DE2+1
           IF ( .NOT. has_char( stagger, 'y' ) ) DE1 = DE1+1
           IF ( ndim .EQ. 3 .AND. .NOT. has_char( stagger, 'z' ) ) DE3 = DE3+1
         CASE ( 'zxy' )
           IF ( .NOT. has_char( stagger, 'x' ) ) DE2 = DE2+1
           IF ( .NOT. has_char( stagger, 'y' ) ) DE3 = DE3+1
           IF ( ndim .EQ. 3 .AND. .NOT. has_char( stagger, 'z' ) ) DE1 = DE1+1
         CASE ( 'xzy' )
           IF ( .NOT. has_char( stagger, 'x' ) ) DE1 = DE1+1
           IF ( .NOT. has_char( stagger, 'y' ) ) DE3 = DE3+1
           IF ( ndim .EQ. 3 .AND. .NOT. has_char( stagger, 'z' ) ) DE2 = DE2+1
         CASE DEFAULT
       END SELECT

     ! moved to here to be after the potential recalculations of D dims
       IF ( wrf_dm_on_monitor() ) THEN
         ALLOCATE ( tmpbuf ( (DE1-DS1+1)*(DE2-DS2+1)*(DE3-DS3+1)/RWORDSIZE*typesize+32 ), STAT=ierr )
       ELSE
         ALLOCATE ( tmpbuf ( 1 ), STAT=ierr )
       END IF
       IF ( ierr .ne. 0 ) CALL wrf_error_fatal ('allocating tmpbuf in wrf_patch_to_global_generic')

       Patch(1,1) = ps1 ; Patch(1,2) = pe1    ! use patch dims
       Patch(2,1) = ps2 ; Patch(2,2) = pe2
       Patch(3,1) = ps3 ; Patch(3,2) = pe3

       IF      ( typesize .EQ. RWORDSIZE ) THEN
         CALL just_patch_r ( buf , locbuf , size(locbuf)*RWORDSIZE/typesize, &
                                   PS1, PE1, PS2, PE2, PS3, PE3 , &
                                   MS1, ME1, MS2, ME2, MS3, ME3   )
       ELSE IF ( typesize .EQ. IWORDSIZE ) THEN
         CALL just_patch_i ( buf , locbuf , size(locbuf)*RWORDSIZE/typesize, &
                                   PS1, PE1, PS2, PE2, PS3, PE3 , &
                                   MS1, ME1, MS2, ME2, MS3, ME3   )
       ELSE IF ( typesize .EQ. DWORDSIZE ) THEN
         CALL just_patch_d ( buf , locbuf , size(locbuf)*RWORDSIZE/typesize, &
                                   PS1, PE1, PS2, PE2, PS3, PE3 , &
                                   MS1, ME1, MS2, ME2, MS3, ME3   )
       ELSE IF ( typesize .EQ. LWORDSIZE ) THEN
         CALL just_patch_l ( buf , locbuf , size(locbuf)*RWORDSIZE/typesize, &
                                   PS1, PE1, PS2, PE2, PS3, PE3 , &
                                   MS1, ME1, MS2, ME2, MS3, ME3   )
       END IF

! defined in external/io_quilt
       CALL collect_on_comm0 (  local_communicator , IWORDSIZE ,  &
                                Patch , 6 ,                       &
                                GPatch , 6*ntasks                 )

       CALL collect_on_comm0 (  local_communicator , typesize ,  &
                                locbuf , (pe1-ps1+1)*(pe2-ps2+1)*(pe3-ps3+1),   &
                                tmpbuf FRSTELEM , (de1-ds1+1)*(de2-ds2+1)*(de3-ds3+1) )

       ndim = len(TRIM(ordering))

       IF ( wrf_at_debug_level(500) ) THEN
         CALL start_timing
       END IF

       IF ( ndim .GE. 2 .AND. wrf_dm_on_monitor() ) THEN

         IF      ( typesize .EQ. RWORDSIZE ) THEN
           CALL patch_2_outbuf_r ( tmpbuf FRSTELEM , globbuf ,             &
                                   DS1, DE1, DS2, DE2, DS3, DE3 , &
                                   GPATCH                         )
         ELSE IF ( typesize .EQ. IWORDSIZE ) THEN
           CALL patch_2_outbuf_i ( tmpbuf FRSTELEM , globbuf ,             &
                                   DS1, DE1, DS2, DE2, DS3, DE3 , &
                                   GPATCH                         )
         ELSE IF ( typesize .EQ. DWORDSIZE ) THEN
           CALL patch_2_outbuf_d ( tmpbuf FRSTELEM , globbuf ,             &
                                   DS1, DE1, DS2, DE2, DS3, DE3 , &
                                   GPATCH                         )
         ELSE IF ( typesize .EQ. LWORDSIZE ) THEN
           CALL patch_2_outbuf_l ( tmpbuf FRSTELEM , globbuf ,             &
                                   DS1, DE1, DS2, DE2, DS3, DE3 , &
                                   GPATCH                         )
         END IF

       END IF

       IF ( wrf_at_debug_level(500) ) THEN
         CALL end_timing('wrf_patch_to_global_generic')
       END IF
       DEALLOCATE( tmpbuf )
#endif
       RETURN
    END SUBROUTINE wrf_patch_to_global_generic

  SUBROUTINE just_patch_i ( inbuf , outbuf, noutbuf,     &
                               PS1,PE1,PS2,PE2,PS3,PE3,  &
                               MS1,ME1,MS2,ME2,MS3,ME3   )
    IMPLICIT NONE
    INTEGER                         , INTENT(IN)  :: noutbuf
    INTEGER    , DIMENSION(noutbuf) , INTENT(OUT) :: outbuf
    INTEGER   MS1,ME1,MS2,ME2,MS3,ME3
    INTEGER   PS1,PE1,PS2,PE2,PS3,PE3
    INTEGER    , DIMENSION( MS1:ME1,MS2:ME2,MS3:ME3 ) , INTENT(IN) :: inbuf
! Local
    INTEGER               :: i,j,k,n   ,  icurs
    icurs = 1
      DO k = PS3, PE3
        DO j = PS2, PE2
          DO i = PS1, PE1
            outbuf( icurs )  = inbuf( i, j, k )
            icurs = icurs + 1
          END DO
        END DO
      END DO
    RETURN
  END SUBROUTINE just_patch_i

  SUBROUTINE just_patch_r ( inbuf , outbuf, noutbuf,     &
                               PS1,PE1,PS2,PE2,PS3,PE3,  &
                               MS1,ME1,MS2,ME2,MS3,ME3   )
    IMPLICIT NONE
    INTEGER                      , INTENT(IN)  :: noutbuf
    REAL    , DIMENSION(noutbuf) , INTENT(OUT) :: outbuf
    INTEGER   MS1,ME1,MS2,ME2,MS3,ME3
    INTEGER   PS1,PE1,PS2,PE2,PS3,PE3
    REAL    , DIMENSION( MS1:ME1,MS2:ME2,MS3:ME3 ) , INTENT(in) :: inbuf
! Local
    INTEGER               :: i,j,k   ,  icurs
    icurs = 1
      DO k = PS3, PE3
        DO j = PS2, PE2
          DO i = PS1, PE1
            outbuf( icurs )  = inbuf( i, j, k )
            icurs = icurs + 1
          END DO
        END DO
      END DO
    RETURN
  END SUBROUTINE just_patch_r

  SUBROUTINE just_patch_d ( inbuf , outbuf, noutbuf,     &
                               PS1,PE1,PS2,PE2,PS3,PE3,  &
                               MS1,ME1,MS2,ME2,MS3,ME3   )
    IMPLICIT NONE
    INTEGER                                  , INTENT(IN)  :: noutbuf
    DOUBLE PRECISION    , DIMENSION(noutbuf) , INTENT(OUT) :: outbuf
    INTEGER   MS1,ME1,MS2,ME2,MS3,ME3
    INTEGER   PS1,PE1,PS2,PE2,PS3,PE3
    DOUBLE PRECISION    , DIMENSION( MS1:ME1,MS2:ME2,MS3:ME3 ) , INTENT(in) :: inbuf
! Local
    INTEGER               :: i,j,k,n   ,  icurs
    icurs = 1
      DO k = PS3, PE3
        DO j = PS2, PE2
          DO i = PS1, PE1
            outbuf( icurs )  = inbuf( i, j, k )
            icurs = icurs + 1
          END DO
        END DO
      END DO
    RETURN
  END SUBROUTINE just_patch_d

  SUBROUTINE just_patch_l ( inbuf , outbuf, noutbuf,     &
                               PS1,PE1,PS2,PE2,PS3,PE3,  &
                               MS1,ME1,MS2,ME2,MS3,ME3   )
    IMPLICIT NONE
    INTEGER                         , INTENT(IN)  :: noutbuf
    LOGICAL    , DIMENSION(noutbuf) , INTENT(OUT) :: outbuf
    INTEGER   MS1,ME1,MS2,ME2,MS3,ME3
    INTEGER   PS1,PE1,PS2,PE2,PS3,PE3
    LOGICAL    , DIMENSION( MS1:ME1,MS2:ME2,MS3:ME3 ) , INTENT(in) :: inbuf
! Local
    INTEGER               :: i,j,k,n   ,  icurs
    icurs = 1
      DO k = PS3, PE3
        DO j = PS2, PE2
          DO i = PS1, PE1
            outbuf( icurs )  = inbuf( i, j, k )
            icurs = icurs + 1
          END DO
        END DO
      END DO
    RETURN
  END SUBROUTINE just_patch_l


  SUBROUTINE patch_2_outbuf_r( inbuf, outbuf,            &
                               DS1,DE1,DS2,DE2,DS3,DE3,  &
                               GPATCH )
    USE module_dm, ONLY : ntasks
    IMPLICIT NONE
    REAL    , DIMENSION(*) , INTENT(IN) :: inbuf
    INTEGER   DS1,DE1,DS2,DE2,DS3,DE3,GPATCH(3,2,ntasks)
    REAL    , DIMENSION( DS1:DE1,DS2:DE2,DS3:DE3 ) , INTENT(out) :: outbuf
! Local
    INTEGER               :: i,j,k,n   ,  icurs
    icurs = 1
    DO n = 1, ntasks
      DO k = GPATCH( 3,1,n ), GPATCH( 3,2,n )
        DO j = GPATCH( 2,1,n ), GPATCH( 2,2,n )
          DO i = GPATCH( 1,1,n ), GPATCH( 1,2,n )
            outbuf( i, j, k ) = inbuf( icurs )
            icurs = icurs + 1
          END DO
        END DO
      END DO
    END DO

    RETURN
  END SUBROUTINE patch_2_outbuf_r

  SUBROUTINE patch_2_outbuf_i( inbuf, outbuf,         &
                               DS1,DE1,DS2,DE2,DS3,DE3,&
                               GPATCH )
    USE module_dm, ONLY : ntasks
    IMPLICIT NONE
    INTEGER    , DIMENSION(*) , INTENT(IN) :: inbuf
    INTEGER   DS1,DE1,DS2,DE2,DS3,DE3,GPATCH(3,2,ntasks)
    INTEGER    , DIMENSION( DS1:DE1,DS2:DE2,DS3:DE3 ) , INTENT(out) :: outbuf
! Local
    INTEGER               :: i,j,k,n   ,  icurs
    icurs = 1
    DO n = 1, ntasks
      DO k = GPATCH( 3,1,n ), GPATCH( 3,2,n )
        DO j = GPATCH( 2,1,n ), GPATCH( 2,2,n )
          DO i = GPATCH( 1,1,n ), GPATCH( 1,2,n )
            outbuf( i, j, k ) = inbuf( icurs )
            icurs = icurs + 1
          END DO
        END DO
      END DO
    END DO
    RETURN
  END SUBROUTINE patch_2_outbuf_i

  SUBROUTINE patch_2_outbuf_d( inbuf, outbuf,         &
                               DS1,DE1,DS2,DE2,DS3,DE3,&
                               GPATCH )
    USE module_dm, ONLY : ntasks
    IMPLICIT NONE
    DOUBLE PRECISION    , DIMENSION(*) , INTENT(IN) :: inbuf
    INTEGER   DS1,DE1,DS2,DE2,DS3,DE3,GPATCH(3,2,ntasks)
    DOUBLE PRECISION    , DIMENSION( DS1:DE1,DS2:DE2,DS3:DE3 ) , INTENT(out) :: outbuf
! Local
    INTEGER               :: i,j,k,n   ,  icurs
    icurs = 1
    DO n = 1, ntasks
      DO k = GPATCH( 3,1,n ), GPATCH( 3,2,n )
        DO j = GPATCH( 2,1,n ), GPATCH( 2,2,n )
          DO i = GPATCH( 1,1,n ), GPATCH( 1,2,n )
            outbuf( i, j, k ) = inbuf( icurs )
            icurs = icurs + 1
          END DO
        END DO
      END DO
    END DO
    RETURN
  END SUBROUTINE patch_2_outbuf_d

  SUBROUTINE patch_2_outbuf_l( inbuf, outbuf,         &
                               DS1,DE1,DS2,DE2,DS3,DE3,&
                               GPATCH )
    USE module_dm, ONLY : ntasks
    IMPLICIT NONE
    LOGICAL    , DIMENSION(*) , INTENT(IN) :: inbuf
    INTEGER   DS1,DE1,DS2,DE2,DS3,DE3,GPATCH(3,2,ntasks)
    LOGICAL    , DIMENSION( DS1:DE1,DS2:DE2,DS3:DE3 ) , INTENT(out) :: outbuf
! Local
    INTEGER               :: i,j,k,n   ,  icurs
    icurs = 1
    DO n = 1, ntasks
      DO k = GPATCH( 3,1,n ), GPATCH( 3,2,n )
        DO j = GPATCH( 2,1,n ), GPATCH( 2,2,n )
          DO i = GPATCH( 1,1,n ), GPATCH( 1,2,n )
            outbuf( i, j, k ) = inbuf( icurs )
            icurs = icurs + 1
          END DO
        END DO
      END DO
    END DO
    RETURN
  END SUBROUTINE patch_2_outbuf_l

!!!!!!!!!!!!!!!!!!!!!!! GLOBAL TO PATCH !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE wrf_global_to_patch_real (globbuf,buf,domdesc,stagger,ordering,&
                                       DS1,DE1,DS2,DE2,DS3,DE3,&
                                       MS1,ME1,MS2,ME2,MS3,ME3,&
                                       PS1,PE1,PS2,PE2,PS3,PE3 )
       IMPLICIT NONE
       INTEGER                         DS1,DE1,DS2,DE2,DS3,DE3,&
                                       MS1,ME1,MS2,ME2,MS3,ME3,&
                                       PS1,PE1,PS2,PE2,PS3,PE3
       CHARACTER *(*) stagger,ordering
       INTEGER fid,domdesc
       REAL globbuf(*)
       REAL buf(*)

       CALL wrf_global_to_patch_generic (globbuf,buf,domdesc,stagger,ordering,RWORDSIZE,&
                                       DS1,DE1,DS2,DE2,DS3,DE3,&
                                       MS1,ME1,MS2,ME2,MS3,ME3,&
                                       PS1,PE1,PS2,PE2,PS3,PE3 )
       RETURN
    END SUBROUTINE wrf_global_to_patch_real

    SUBROUTINE wrf_global_to_patch_double (globbuf,buf,domdesc,stagger,ordering,&
                                       DS1,DE1,DS2,DE2,DS3,DE3,&
                                       MS1,ME1,MS2,ME2,MS3,ME3,&
                                       PS1,PE1,PS2,PE2,PS3,PE3 )
       IMPLICIT NONE
       INTEGER                         DS1,DE1,DS2,DE2,DS3,DE3,&
                                       MS1,ME1,MS2,ME2,MS3,ME3,&
                                       PS1,PE1,PS2,PE2,PS3,PE3
       CHARACTER *(*) stagger,ordering
       INTEGER fid,domdesc
! this next declaration is REAL, not DOUBLE PRECISION because it will be autopromoted
! to double precision by the compiler when WRF is compiled for 8 byte reals. Only reason
! for having this separate routine is so we pass the correct MPI type to mpi_scatterv
! since we were not indexing the globbuf and Field arrays it does not matter
       REAL globbuf(*)
       REAL buf(*)

       CALL wrf_global_to_patch_generic (globbuf,buf,domdesc,stagger,ordering,DWORDSIZE,&
                                       DS1,DE1,DS2,DE2,DS3,DE3,&
                                       MS1,ME1,MS2,ME2,MS3,ME3,&
                                       PS1,PE1,PS2,PE2,PS3,PE3 )
       RETURN
    END SUBROUTINE wrf_global_to_patch_double


    SUBROUTINE wrf_global_to_patch_integer (globbuf,buf,domdesc,stagger,ordering,&
                                       DS1,DE1,DS2,DE2,DS3,DE3,&
                                       MS1,ME1,MS2,ME2,MS3,ME3,&
                                       PS1,PE1,PS2,PE2,PS3,PE3 )
       IMPLICIT NONE
       INTEGER                         DS1,DE1,DS2,DE2,DS3,DE3,&
                                       MS1,ME1,MS2,ME2,MS3,ME3,&
                                       PS1,PE1,PS2,PE2,PS3,PE3
       CHARACTER *(*) stagger,ordering
       INTEGER fid,domdesc
       INTEGER globbuf(*)
       INTEGER buf(*)

       CALL wrf_global_to_patch_generic (globbuf,buf,domdesc,stagger,ordering,IWORDSIZE,&
                                       DS1,DE1,DS2,DE2,DS3,DE3,&
                                       MS1,ME1,MS2,ME2,MS3,ME3,&
                                       PS1,PE1,PS2,PE2,PS3,PE3 )
       RETURN
    END SUBROUTINE wrf_global_to_patch_integer

    SUBROUTINE wrf_global_to_patch_logical (globbuf,buf,domdesc,stagger,ordering,&
                                       DS1,DE1,DS2,DE2,DS3,DE3,&
                                       MS1,ME1,MS2,ME2,MS3,ME3,&
                                       PS1,PE1,PS2,PE2,PS3,PE3 )
       IMPLICIT NONE
       INTEGER                         DS1,DE1,DS2,DE2,DS3,DE3,&
                                       MS1,ME1,MS2,ME2,MS3,ME3,&
                                       PS1,PE1,PS2,PE2,PS3,PE3
       CHARACTER *(*) stagger,ordering
       INTEGER fid,domdesc
       LOGICAL globbuf(*)
       LOGICAL buf(*)

       CALL wrf_global_to_patch_generic (globbuf,buf,domdesc,stagger,ordering,LWORDSIZE,&
                                       DS1,DE1,DS2,DE2,DS3,DE3,&
                                       MS1,ME1,MS2,ME2,MS3,ME3,&
                                       PS1,PE1,PS2,PE2,PS3,PE3 )
       RETURN
    END SUBROUTINE wrf_global_to_patch_logical

    SUBROUTINE wrf_global_to_patch_generic (globbuf,buf,domdesc,stagger,ordering,typesize,&
                                       DS1a,DE1a,DS2a,DE2a,DS3a,DE3a,&
                                       MS1a,ME1a,MS2a,ME2a,MS3a,ME3a,&
                                       PS1a,PE1a,PS2a,PE2a,PS3a,PE3a )
       USE module_dm, ONLY : local_communicator, ntasks
       USE module_driver_constants
       IMPLICIT NONE
       INTEGER                         DS1a,DE1a,DS2a,DE2a,DS3a,DE3a,&
                                       MS1a,ME1a,MS2a,ME2a,MS3a,ME3a,&
                                       PS1a,PE1a,PS2a,PE2a,PS3a,PE3A
       CHARACTER *(*) stagger,ordering
       INTEGER domdesc,typesize,ierr
       REAL globbuf(*)
       REAL buf(*)
#ifndef STUBMPI
       INTEGER                         DS1,DE1,DS2,DE2,DS3,DE3,&
                                       MS1,ME1,MS2,ME2,MS3,ME3,&
                                       PS1,PE1,PS2,PE2,PS3,PE3
       LOGICAL, EXTERNAL :: wrf_dm_on_monitor, has_char

       INTEGER i,j,k,ord,ord2d,ndim
       INTEGER  Patch(3,2), Gpatch(3,2,ntasks)
       REAL, ALLOCATABLE :: tmpbuf( : )
       REAL locbuf( (PE1a-PS1a+1)*(PE2a-PS2a+1)*(PE3a-PS3a+1)/RWORDSIZE*typesize+32 )

       DS1 = DS1a ; DE1 = DE1a ; DS2=DS2a ; DE2 = DE2a ; DS3 = DS3a ; DE3 = DE3a
       MS1 = MS1a ; ME1 = ME1a ; MS2=MS2a ; ME2 = ME2a ; MS3 = MS3a ; ME3 = ME3a
       PS1 = PS1a ; PE1 = PE1a ; PS2=PS2a ; PE2 = PE2a ; PS3 = PS3a ; PE3 = PE3a

       SELECT CASE ( TRIM(ordering) )
         CASE ( 'xy', 'yx' )
           ndim = 2
         CASE DEFAULT
           ndim = 3   ! where appropriate
       END SELECT

       SELECT CASE ( TRIM(ordering) )
         CASE ( 'xyz','xy' )
            ! the non-staggered variables come in at one-less than
            ! domain dimensions, but code wants full domain spec, so
            ! adjust if not staggered
           IF ( .NOT. has_char( stagger, 'x' ) ) DE1 = DE1+1
           IF ( .NOT. has_char( stagger, 'y' ) ) DE2 = DE2+1
           IF ( ndim .EQ. 3 .AND. .NOT. has_char( stagger, 'z' ) ) DE3 = DE3+1
         CASE ( 'yxz','yx' )
           IF ( .NOT. has_char( stagger, 'x' ) ) DE2 = DE2+1
           IF ( .NOT. has_char( stagger, 'y' ) ) DE1 = DE1+1
           IF ( ndim .EQ. 3 .AND. .NOT. has_char( stagger, 'z' ) ) DE3 = DE3+1
         CASE ( 'zxy' )
           IF ( .NOT. has_char( stagger, 'x' ) ) DE2 = DE2+1
           IF ( .NOT. has_char( stagger, 'y' ) ) DE3 = DE3+1
           IF ( ndim .EQ. 3 .AND. .NOT. has_char( stagger, 'z' ) ) DE1 = DE1+1
         CASE ( 'xzy' )
           IF ( .NOT. has_char( stagger, 'x' ) ) DE1 = DE1+1
           IF ( .NOT. has_char( stagger, 'y' ) ) DE3 = DE3+1
           IF ( ndim .EQ. 3 .AND. .NOT. has_char( stagger, 'z' ) ) DE2 = DE2+1
         CASE DEFAULT
       END SELECT

     ! moved to here to be after the potential recalculations of D dims
       IF ( wrf_dm_on_monitor() ) THEN
         ALLOCATE ( tmpbuf ( (DE1-DS1+1)*(DE2-DS2+1)*(DE3-DS3+1)/RWORDSIZE*typesize+32 ), STAT=ierr )
       ELSE
         ALLOCATE ( tmpbuf ( 1 ), STAT=ierr )
       END IF
       IF ( ierr .ne. 0 ) CALL wrf_error_fatal ('allocating tmpbuf in wrf_global_to_patch_generic')

       Patch(1,1) = ps1 ; Patch(1,2) = pe1    ! use patch dims
       Patch(2,1) = ps2 ; Patch(2,2) = pe2
       Patch(3,1) = ps3 ; Patch(3,2) = pe3

! defined in external/io_quilt
       CALL collect_on_comm0 (  local_communicator , IWORDSIZE ,  &
                                Patch , 6 ,                       &
                                GPatch , 6*ntasks                 )
       ndim = len(TRIM(ordering))

       IF ( wrf_dm_on_monitor() .AND. ndim .GE. 2 ) THEN
         IF      ( typesize .EQ. RWORDSIZE ) THEN
           CALL outbuf_2_patch_r ( globbuf , tmpbuf FRSTELEM ,    &
                                   DS1, DE1, DS2, DE2, DS3, DE3 , &
                                   MS1, ME1, MS2, ME2, MS3, ME3 , &
                                   GPATCH                         )
         ELSE IF ( typesize .EQ. IWORDSIZE ) THEN
           CALL outbuf_2_patch_i ( globbuf , tmpbuf FRSTELEM ,    &
                                   DS1, DE1, DS2, DE2, DS3, DE3 , &
                                   GPATCH                         )
         ELSE IF ( typesize .EQ. DWORDSIZE ) THEN
           CALL outbuf_2_patch_d ( globbuf , tmpbuf FRSTELEM ,    &
                                   DS1, DE1, DS2, DE2, DS3, DE3 , &
                                   GPATCH                         )
         ELSE IF ( typesize .EQ. LWORDSIZE ) THEN
           CALL outbuf_2_patch_l ( globbuf , tmpbuf FRSTELEM ,    &
                                   DS1, DE1, DS2, DE2, DS3, DE3 , &
                                   GPATCH                         )
         END IF
       END IF

       CALL dist_on_comm0 (  local_communicator , typesize ,  &
                             tmpbuf FRSTELEM , (de1-ds1+1)*(de2-ds2+1)*(de3-ds3+1) , &
                             locbuf    , (pe1-ps1+1)*(pe2-ps2+1)*(pe3-ps3+1) )

       IF      ( typesize .EQ. RWORDSIZE ) THEN
         CALL all_sub_r ( locbuf , buf ,             &
                                   PS1, PE1, PS2, PE2, PS3, PE3 , &
                                   MS1, ME1, MS2, ME2, MS3, ME3   )

       ELSE IF ( typesize .EQ. IWORDSIZE ) THEN
         CALL all_sub_i ( locbuf , buf ,             &
                                   PS1, PE1, PS2, PE2, PS3, PE3 , &
                                   MS1, ME1, MS2, ME2, MS3, ME3   )
       ELSE IF ( typesize .EQ. DWORDSIZE ) THEN
         CALL all_sub_d ( locbuf , buf ,             &
                                   PS1, PE1, PS2, PE2, PS3, PE3 , &
                                   MS1, ME1, MS2, ME2, MS3, ME3   )
       ELSE IF ( typesize .EQ. LWORDSIZE ) THEN
         CALL all_sub_l ( locbuf , buf ,             &
                                   PS1, PE1, PS2, PE2, PS3, PE3 , &
                                   MS1, ME1, MS2, ME2, MS3, ME3   )
       END IF


       DEALLOCATE ( tmpbuf )
#endif
       RETURN
    END SUBROUTINE wrf_global_to_patch_generic

  SUBROUTINE all_sub_i ( inbuf , outbuf,              &
                               PS1,PE1,PS2,PE2,PS3,PE3,  &
                               MS1,ME1,MS2,ME2,MS3,ME3   )
    IMPLICIT NONE
    INTEGER    , DIMENSION(*) , INTENT(IN) :: inbuf
    INTEGER   MS1,ME1,MS2,ME2,MS3,ME3
    INTEGER   PS1,PE1,PS2,PE2,PS3,PE3
    INTEGER    , DIMENSION( MS1:ME1,MS2:ME2,MS3:ME3 ) , INTENT(OUT) :: outbuf
! Local
    INTEGER               :: i,j,k,n   ,  icurs
    icurs = 1
      DO k = PS3, PE3
        DO j = PS2, PE2
          DO i = PS1, PE1
            outbuf( i, j, k )  = inbuf ( icurs )
            icurs = icurs + 1
          END DO
        END DO
      END DO
    RETURN
  END SUBROUTINE all_sub_i

  SUBROUTINE all_sub_r ( inbuf , outbuf,              &
                               PS1,PE1,PS2,PE2,PS3,PE3,  &
                               MS1,ME1,MS2,ME2,MS3,ME3   )
    IMPLICIT NONE
    REAL       , DIMENSION(*) , INTENT(IN) :: inbuf
    INTEGER   MS1,ME1,MS2,ME2,MS3,ME3
    INTEGER   PS1,PE1,PS2,PE2,PS3,PE3
    REAL       , DIMENSION( MS1:ME1,MS2:ME2,MS3:ME3 ) , INTENT(OUT) :: outbuf
! Local
    INTEGER               :: i,j,k,n   ,  icurs
    icurs = 1
      DO k = PS3, PE3
        DO j = PS2, PE2
          DO i = PS1, PE1
            outbuf( i, j, k )  = inbuf ( icurs )
            icurs = icurs + 1
          END DO
        END DO
      END DO

    RETURN
  END SUBROUTINE all_sub_r

  SUBROUTINE all_sub_d ( inbuf , outbuf,              &
                               PS1,PE1,PS2,PE2,PS3,PE3,  &
                               MS1,ME1,MS2,ME2,MS3,ME3   )
    IMPLICIT NONE
    DOUBLE PRECISION    , DIMENSION(*) , INTENT(IN) :: inbuf
    INTEGER   MS1,ME1,MS2,ME2,MS3,ME3
    INTEGER   PS1,PE1,PS2,PE2,PS3,PE3
    DOUBLE PRECISION    , DIMENSION( MS1:ME1,MS2:ME2,MS3:ME3 ) , INTENT(OUT) :: outbuf
! Local
    INTEGER               :: i,j,k,n   ,  icurs
    icurs = 1
      DO k = PS3, PE3
        DO j = PS2, PE2
          DO i = PS1, PE1
            outbuf( i, j, k )  = inbuf ( icurs )
            icurs = icurs + 1
          END DO
        END DO
      END DO
    RETURN
  END SUBROUTINE all_sub_d

  SUBROUTINE all_sub_l ( inbuf , outbuf,              &
                               PS1,PE1,PS2,PE2,PS3,PE3,  &
                               MS1,ME1,MS2,ME2,MS3,ME3   )
    IMPLICIT NONE
    LOGICAL    , DIMENSION(*) , INTENT(IN) :: inbuf
    INTEGER   MS1,ME1,MS2,ME2,MS3,ME3
    INTEGER   PS1,PE1,PS2,PE2,PS3,PE3
    LOGICAL    , DIMENSION( MS1:ME1,MS2:ME2,MS3:ME3 ) , INTENT(OUT) :: outbuf
! Local
    INTEGER               :: i,j,k,n   ,  icurs
    icurs = 1
      DO k = PS3, PE3
        DO j = PS2, PE2
          DO i = PS1, PE1
            outbuf( i, j, k )  = inbuf ( icurs )
            icurs = icurs + 1
          END DO
        END DO
      END DO
    RETURN
  END SUBROUTINE all_sub_l

  SUBROUTINE outbuf_2_patch_r( inbuf, outbuf,         &
                               DS1,DE1,DS2,DE2,DS3,DE3, &
                               MS1, ME1, MS2, ME2, MS3, ME3 , &
                               GPATCH )
    USE module_dm, ONLY : ntasks
    IMPLICIT NONE
    REAL    , DIMENSION(*) , INTENT(OUT) :: outbuf
    INTEGER   DS1,DE1,DS2,DE2,DS3,DE3,GPATCH(3,2,ntasks)
    INTEGER   MS1,ME1,MS2,ME2,MS3,ME3
    REAL    , DIMENSION( DS1:DE1,DS2:DE2,DS3:DE3 ) , INTENT(IN) :: inbuf
! Local
    INTEGER               :: i,j,k,n   ,  icurs

    icurs = 1
    DO n = 1, ntasks
      DO k = GPATCH( 3,1,n ), GPATCH( 3,2,n )
        DO j = GPATCH( 2,1,n ), GPATCH( 2,2,n )
          DO i = GPATCH( 1,1,n ), GPATCH( 1,2,n )
            outbuf( icurs ) = inbuf( i,j,k )
            icurs = icurs + 1
          END DO
        END DO
      END DO
    END DO
    RETURN
  END SUBROUTINE outbuf_2_patch_r

  SUBROUTINE outbuf_2_patch_i( inbuf, outbuf,         &
                               DS1,DE1,DS2,DE2,DS3,DE3,&
                               GPATCH )
    USE module_dm, ONLY : ntasks
    IMPLICIT NONE
    INTEGER    , DIMENSION(*) , INTENT(OUT) :: outbuf
    INTEGER   DS1,DE1,DS2,DE2,DS3,DE3,GPATCH(3,2,ntasks)
    INTEGER    , DIMENSION( DS1:DE1,DS2:DE2,DS3:DE3 ) , INTENT(IN) :: inbuf
! Local
    INTEGER               :: i,j,k,n   ,  icurs
    icurs = 1
    DO n = 1, ntasks
      DO k = GPATCH( 3,1,n ), GPATCH( 3,2,n )
        DO j = GPATCH( 2,1,n ), GPATCH( 2,2,n )
          DO i = GPATCH( 1,1,n ), GPATCH( 1,2,n )
            outbuf( icurs ) = inbuf( i,j,k )
            icurs = icurs + 1
          END DO
        END DO
      END DO
    END DO
    RETURN
  END SUBROUTINE outbuf_2_patch_i

  SUBROUTINE outbuf_2_patch_d( inbuf, outbuf,         &
                               DS1,DE1,DS2,DE2,DS3,DE3,&
                               GPATCH )
    USE module_dm, ONLY : ntasks
    IMPLICIT NONE
    DOUBLE PRECISION    , DIMENSION(*) , INTENT(OUT) :: outbuf
    INTEGER   DS1,DE1,DS2,DE2,DS3,DE3,GPATCH(3,2,ntasks)
    DOUBLE PRECISION    , DIMENSION( DS1:DE1,DS2:DE2,DS3:DE3 ) , INTENT(IN) :: inbuf
! Local
    INTEGER               :: i,j,k,n   ,  icurs
    icurs = 1
    DO n = 1, ntasks
      DO k = GPATCH( 3,1,n ), GPATCH( 3,2,n )
        DO j = GPATCH( 2,1,n ), GPATCH( 2,2,n )
          DO i = GPATCH( 1,1,n ), GPATCH( 1,2,n )
            outbuf( icurs ) = inbuf( i,j,k )
            icurs = icurs + 1
          END DO
        END DO
      END DO
    END DO
    RETURN
  END SUBROUTINE outbuf_2_patch_d

  SUBROUTINE outbuf_2_patch_l( inbuf, outbuf,         &
                               DS1,DE1,DS2,DE2,DS3,DE3,&
                               GPATCH )
    USE module_dm, ONLY : ntasks
    IMPLICIT NONE
    LOGICAL    , DIMENSION(*) , INTENT(OUT) :: outbuf
    INTEGER   DS1,DE1,DS2,DE2,DS3,DE3,GPATCH(3,2,ntasks)
    LOGICAL    , DIMENSION( DS1:DE1,DS2:DE2,DS3:DE3 ) , INTENT(IN) :: inbuf
! Local
    INTEGER               :: i,j,k,n   ,  icurs
    icurs = 1
    DO n = 1, ntasks
      DO k = GPATCH( 3,1,n ), GPATCH( 3,2,n )
        DO j = GPATCH( 2,1,n ), GPATCH( 2,2,n )
          DO i = GPATCH( 1,1,n ), GPATCH( 1,2,n )
            outbuf( icurs ) = inbuf( i,j,k )
            icurs = icurs + 1
          END DO
        END DO
      END DO
    END DO
    RETURN
  END SUBROUTINE outbuf_2_patch_l


  SUBROUTINE wrf_dm_nestexchange_init
      CALL rsl_lite_nesting_reset
  END SUBROUTINE wrf_dm_nestexchange_init


!------------------------------------------------------------------

#if ( EM_CORE == 1 && DA_CORE != 1 )

!------------------------------------------------------------------

   SUBROUTINE force_domain_em_part2 ( grid, ngrid, pgrid, config_flags    &
!
#include "dummy_new_args.inc"
!
                 )
      USE module_state_description
      USE module_domain, ONLY : domain, get_ijk_from_grid
      USE module_configure, ONLY : grid_config_rec_type
      USE module_dm, ONLY : ntasks, ntasks_x, ntasks_y, local_communicator, mytask, &
                            nest_pes_x, nest_pes_y ! ,                                 &
                            !push_communicators_for_domain,pop_communicators_for_domain
      USE module_comm_nesting_dm, ONLY : halo_force_down_sub
      USE module_model_constants
      IMPLICIT NONE
!
      TYPE(domain), POINTER :: grid          ! name of the grid being dereferenced (must be "grid")
      TYPE(domain), POINTER :: ngrid
      TYPE(domain), POINTER :: pgrid         !KAL added for vertical nesting
#include "dummy_new_decl.inc"
      INTEGER nlev, msize
      INTEGER i,j,pig,pjg,cm,cn,nig,njg,retval,k,kk
      TYPE (grid_config_rec_type)            :: config_flags
      REAL xv(2000)
      INTEGER       ::          cids, cide, cjds, cjde, ckds, ckde,    &
                                cims, cime, cjms, cjme, ckms, ckme,    &
                                cips, cipe, cjps, cjpe, ckps, ckpe
      INTEGER       ::          nids, nide, njds, njde, nkds, nkde,    &
                                nims, nime, njms, njme, nkms, nkme,    &
                                nips, nipe, njps, njpe, nkps, nkpe
      INTEGER       ::          ids, ide, jds, jde, kds, kde,    &
                                ims, ime, jms, jme, kms, kme,    &
                                ips, ipe, jps, jpe, kps, kpe
      INTEGER idim1,idim2,idim3,idim4,idim5,idim6,idim7,itrace
      REAL  dummy_xs, dummy_xe, dummy_ys, dummy_ye

      !KAL variables for vertical nesting
      REAL :: p_top_m  , p_surf_m , mu_m , hsca_m , pre_c ,pre_n
      REAL, DIMENSION(pgrid%s_vert:pgrid%e_vert) :: alt_w_c
      REAL, DIMENSION(pgrid%s_vert:pgrid%e_vert+1) :: alt_u_c
      REAL, DIMENSION(ngrid%s_vert:ngrid%e_vert) :: alt_w_n
      REAL, DIMENSION(ngrid%s_vert:ngrid%e_vert+1) :: alt_u_n      
      
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: p, al
      REAL :: pfu, pfd, phm, temp, qvf, qvf1, qvf2    

      !KAL change this for vertical nesting
      ! force_domain_em_part1 packs up the interpolation onto the coarse (vertical) grid
      ! therefore the message size is based on the coarse grid number of levels
      ! here it is unpacked onto the intermediate grid
      CALL get_ijk_from_grid (  pgrid ,                   &
                                cids, cide, cjds, cjde, ckds, ckde,    &
                                cims, cime, cjms, cjme, ckms, ckme,    &
                                cips, cipe, cjps, cjpe, ckps, ckpe    )
                                
      !KAL this is the original WRF code                                
      !CALL get_ijk_from_grid (  grid ,                   &
      !                          cids, cide, cjds, cjde, ckds, ckde,    &
      !                          cims, cime, cjms, cjme, ckms, ckme,    &
      !                          cips, cipe, cjps, cjpe, ckps, ckpe    )
      CALL get_ijk_from_grid (  ngrid ,              &
                                nids, nide, njds, njde, nkds, nkde,    &
                                nims, nime, njms, njme, nkms, nkme,    &
                                nips, nipe, njps, njpe, nkps, nkpe    )

      nlev  = ckde - ckds + 1

#include "nest_interpdown_unpack.inc"

if (ngrid%vert_refine_method .NE. 0) then

      !KAL calculating the vertical coordinate for parent and nest grid (code from ndown)
      ! assume that the parent and nest have the same p_top value (as in ndown)
      
!KAL ckde is equal to e_vert of the coarse grid. There are e_vert-1 u points.  The coarse 1D grid here is e_vert+1,
!    so it is the e_vert-1 points from the coarse grid, plus a surface point plus a top point.  Extrapolation coefficients
!    are used to get the surface and top points to fill out the pro_u_c 1D array of u values from the coarse grid.     
                
      hsca_m = 6.7 !KAL scale height of the atmosphere
      p_top_m = ngrid%p_top
      p_surf_m = 1.e5
      mu_m = p_surf_m - p_top_m
!    parent
      do  k = 1,ckde
      pre_c = mu_m * pgrid%c3f(k) + p_top_m + pgrid%c4f(k)
      alt_w_c(k) =  -hsca_m * alog(pre_c/p_surf_m)
      enddo   
      do  k = 1,ckde-1
      pre_c = mu_m * pgrid%c3h(k) + p_top_m + pgrid%c4h(k)
      alt_u_c(k+1) =  -hsca_m * alog(pre_c/p_surf_m)
      enddo
      alt_u_c(1) =  alt_w_c(1)
      alt_u_c(ckde+1) =  alt_w_c(ckde)       
!    nest
      do  k = 1,nkde
      pre_n = mu_m * ngrid%c3f(k) + p_top_m + ngrid%c4f(k)
      alt_w_n(k) =  -hsca_m * alog(pre_n/p_surf_m)
      enddo
      do  k = 1,nkde-1
      pre_n = mu_m * ngrid%c3h(k) + p_top_m + ngrid%c4h(k)
      alt_u_n(k+1) =  -hsca_m * alog(pre_n/p_surf_m)
      enddo
      alt_u_n(1) =  alt_w_n(1)
      alt_u_n(nkde+1) =  alt_w_n(nkde)
        
endif   

      !KAL added this call for vertical nesting (return coarse grid dimensions to intended values)
      CALL get_ijk_from_grid (  grid ,                   &
                                cids, cide, cjds, cjde, ckds, ckde,    &
                                cims, cime, cjms, cjme, ckms, ckme,    &
                                cips, cipe, cjps, cjpe, ckps, ckpe    )
                                                      
      CALL get_ijk_from_grid (  grid ,              &
                                ids, ide, jds, jde, kds, kde,    &
                                ims, ime, jms, jme, kms, kme,    &
                                ips, ipe, jps, jpe, kps, kpe    )

      !  Vertical refinement is turned on.

      IF (ngrid%vert_refine_method .NE. 0) THEN
      
#include "nest_forcedown_interp_vert.inc"

         IF ( ngrid%this_is_an_ideal_run ) THEN
            IF ( SIZE( grid%t_init, 1 ) * SIZE( grid%t_init, 3 ) .GT. 1 ) THEN
               CALL vert_interp_vert_nesting( grid%t_init, & !CD field
                                              ids, ide, kds, kde, jds, jde, & !CD dims
                                              ims, ime, kms, kme, jms, jme, & !CD dims
                                              ips, ipe, kps, MIN( (kde-1), kpe ), jps, jpe, & !CD dims
                                              pgrid%s_vert, pgrid%e_vert, & !vertical dimension of the parent grid
                                              pgrid%cf1, pgrid%cf2, pgrid%cf3, pgrid%cfn, pgrid%cfn1, & !coarse grid extrapolation constants
                                              alt_u_c, alt_u_n ) !coordinates for parent and nest
            END IF  ! Check t_init is a fully allocated 3d array.
         END IF ! only for ideal runs


         !  Rebalance the grid on the intermediate grid.  The intermediate grid has the horizontal
         !  resolution of the parent grid, but at this point has been interpolated in the vertical
         !  to the resolution of the nest.  The base state (phb, pb, etc) from the parent grid is
         !  unpacked onto the intermediate grid every time this subroutine is called.  We need the
         !  base state of the nest, so it is recalculated here.  

         !  Additionally, we do not need to vertically interpolate the entire intermediate grid
         !  above, just the points that contribute to the boundary forcing.

         !  Base state potential temperature and inverse density (alpha = 1/rho) from
         !  the half eta levels and the base-profile surface pressure.  Compute 1/rho
         !  from equation of state.  The potential temperature is a perturbation from t0.
      
         !  Uncouple the variables moist and t_2 that are used to calculate ph_2

         DO j = MAX(jds,jps),MIN(jde-1,jpe)
            DO i = MAX(ids,ips),MIN(ide-1,ipe)
               DO k=kds,kde-1
                  grid%t_2(i,k,j) = grid%t_2(i,k,j)/((ngrid%c1h(k)*grid%mub(i,j)+ngrid%c2h(k)) + (ngrid%c1h(k)*grid%mu_2(i,j)))
                  moist(i,k,j,P_QV) = moist(i,k,j,P_QV)/((ngrid%c1h(k)*grid%mub(i,j)+ngrid%c2h(k)) + (ngrid%c1h(k)*grid%mu_2(i,j)))
               END DO
            END DO
         END DO
    
         DO j = MAX(jds,jps),MIN(jde-1,jpe)
            DO i = MAX(ids,ips),MIN(ide-1,ipe)

               DO k = 1, kpe-1
                  grid%pb(i,k,j) = ngrid%c3h(k) * grid%mub(i,j) + ngrid%c4h(k) + ngrid%p_top
             
                  !  If this is a real run, recalc t_init.
   
                  IF ( .NOT. ngrid%this_is_an_ideal_run ) THEN
                     temp = MAX ( ngrid%tiso, ngrid%t00 + ngrid%tlp*LOG(grid%pb(i,k,j)/ngrid%p00) )
                     IF ( grid%pb(i,k,j) .LT. ngrid%p_strat ) THEN
                        temp = ngrid%tiso + ngrid%tlp_strat * LOG ( grid%pb(i,k,j)/ngrid%p_strat )
                     END IF
                     grid%t_init(i,k,j) = temp*(ngrid%p00/grid%pb(i,k,j))**(r_d/cp) - t0
                  END IF
                  grid%alb(i,k,j) = (r_d/p1000mb)*(grid%t_init(i,k,j)+t0)*(grid%pb(i,k,j)/p1000mb)**cvpm
               END DO
   
               !  Integrate base geopotential, starting at terrain elevation.  This assures that
               !  the base state is in exact hydrostatic balance with respect to the model equations.
               !  This field is on full levels.
   
               grid%phb(i,1,j) = grid%ht(i,j) * g
               IF (grid%hypsometric_opt == 1) THEN
                  DO kk = 2,kpe
                     k = kk - 1
                     grid%phb(i,kk,j) = grid%phb(i,k,j) - ngrid%dnw(k)*(ngrid%c1h(k)*grid%mub(i,j)+ngrid%c2h(k))*grid%alb(i,k,j)
                  END DO
               ELSE IF (grid%hypsometric_opt == 2) THEN
                  DO k = 2,kpe
                     pfu = ngrid%c3f(k  )*grid%MUB(i,j) + ngrid%c4f(k  ) + ngrid%p_top
                     pfd = ngrid%c3f(k-1)*grid%MUB(i,j) + ngrid%c4f(k-1) + ngrid%p_top
                     phm = ngrid%c3h(k-1)*grid%MUB(i,j) + ngrid%c4h(k-1) + ngrid%p_top
                     grid%phb(i,k,j) = grid%phb(i,k-1,j) + grid%alb(i,k-1,j)*phm*LOG(pfd/pfu)
                  END DO
               ELSE
                  CALL wrf_error_fatal( 'module_dm: hypsometric_opt should be 1 or 2' )
               END IF  ! which hypsometric option
            END DO  ! i loop
         END DO  ! j loop
         !  Perturbation fields
         ALLOCATE( p (ips:ipe, kps:kpe, jps:jpe) )
         ALLOCATE( al(ips:ipe, kps:kpe, jps:jpe) )
         DO j = MAX(jds,jps),MIN(jde-1,jpe)
            DO i = MAX(ids,ips),MIN(ide-1,ipe)
               !  Integrate the hydrostatic equation (from the RHS of the bigstep vertical momentum
               !  equation) down from the top to get the pressure perturbation.  First get the pressure
               !  perturbation, moisture, and inverse density (total and perturbation) at the top-most level.
      
               kk = kpe-1
               k = kk+1
      
               qvf1 = 0.5*(moist(i,kk,j,P_QV)+moist(i,kk,j,P_QV))
               qvf2 = 1./(1.+qvf1)
               qvf1 = qvf1*qvf2
      
               p(i,kk,j) = - 0.5*((ngrid%c1f(k)*grid%Mu_2(i,j))+qvf1*(ngrid%c1f(k)*grid%Mub(i,j)+ngrid%c2f(k)))/ngrid%rdnw(kk)/qvf2
               IF ( config_flags%use_theta_m == 0) THEN
                  qvf = 1. + rvovrd*moist(i,kk,j,P_QV)
               ELSE
                  qvf = 1.
               ENDIF
               al(i,kk,j) = (r_d/p1000mb)*(grid%t_2(i,kk,j)+t0)*qvf* &
                           (((p(i,kk,j)+grid%pb(i,kk,j))/p1000mb)**cvpm) - grid%alb(i,kk,j)
      
               !  Now, integrate down the column to compute the pressure perturbation, and diagnose the two
               !  inverse density fields (total and perturbation).
      
               DO kk=kpe-2,1,-1
                  k = kk + 1
                  qvf1 = 0.5*(moist(i,kk,j,P_QV)+moist(i,kk+1,j,P_QV))
                  qvf2 = 1./(1.+qvf1)
                  qvf1 = qvf1*qvf2
                  p(i,kk,j) = p(i,kk+1,j) - ((ngrid%c1f(k)*grid%Mu_2(i,j)) + qvf1*(ngrid%c1f(k)*grid%Mub(i,j)+ngrid%c2f(k)))/qvf2/ngrid%rdn(kk+1)
                  IF ( config_flags%use_theta_m == 0) THEN
                     qvf = 1. + rvovrd*moist(i,kk,j,P_QV)
                  ELSE
                     qvf = 1.
                  ENDIF
                  al(i,kk,j) = (r_d/p1000mb)*(grid%t_2(i,kk,j)+t0)*qvf* &
                              (((p(i,kk,j)+grid%pb(i,kk,j))/p1000mb)**cvpm) - grid%alb(i,kk,j)
               END DO
      
               !  This is the hydrostatic equation used in the model after the small timesteps.  In
               !  the model, grid%al (inverse density) is computed from the geopotential.
      
               IF (grid%hypsometric_opt == 1) THEN
                  DO kk = 2,kpe
                     k = kk - 1
                     grid%ph_2(i,kk,j) = grid%ph_2(i,kk-1,j) - &
                                        ngrid%dnw(kk-1) * ( ((ngrid%c1h(k)*grid%mub(i,j)+ngrid%c2h(k))+(ngrid%c1h(k)*grid%mu_2(i,j)))*al(i,kk-1,j) &
                                        + (ngrid%c1h(k)*grid%mu_2(i,j))*grid%alb(i,kk-1,j) )
                  END DO
      
               ! Alternative hydrostatic eq.: dZ = -al*p*dLOG(p), where p is dry pressure.
               ! Note that al*p approximates Rd*T and dLOG(p) does z.
               ! Here T varies mostly linear with z, the first-order integration produces better result.
      
               ELSE IF (grid%hypsometric_opt == 2) THEN
      
                  grid%ph_2(i,1,j) = grid%phb(i,1,j)
                  DO k = 2,kpe
                     pfu = ngrid%c3f(k  )*( grid%MUB(i,j)+grid%MU_2(i,j) ) + ngrid%c4f(k  ) + ngrid%p_top
                     pfd = ngrid%c3f(k-1)*( grid%MUB(i,j)+grid%MU_2(i,j) ) + ngrid%c4f(k-1) + ngrid%p_top
                     phm = ngrid%c3h(k-1)*( grid%MUB(i,j)+grid%MU_2(i,j) ) + ngrid%c4h(k-1) + ngrid%p_top
                     grid%ph_2(i,k,j) = grid%ph_2(i,k-1,j) + (grid%alb(i,k-1,j)+al(i,k-1,j))*phm*LOG(pfd/pfu)
                  END DO
      
                  DO k = 1,kpe
                     grid%ph_2(i,k,j) = grid%ph_2(i,k,j) - grid%phb(i,k,j)
                  END DO
      
               END IF

            END DO ! i loop
         END DO ! j loop

         DEALLOCATE(p)
         DEALLOCATE(al)
      
         ! Couple the variables moist and t_2, and the newly calculated ph_2
         DO j = MAX(jds,jps),MIN(jde-1,jpe)
            DO i = MAX(ids,ips),MIN(ide-1,ipe)
               DO k=kps,kpe
               grid%ph_2(i,k,j) = grid%ph_2(i,k,j)*((ngrid%c1f(k)*grid%Mub(i,j)+ngrid%c2f(k)) + (ngrid%c1f(k)*grid%Mu_2(i,j)))
               END DO
            END DO
         END DO
         DO j = MAX(jds,jps),MIN(jde-1,jpe)
            DO i = MAX(ids,ips),MIN(ide-1,ipe)
               DO k=kps,kpe-1
               grid%t_2(i,k,j) = grid%t_2(i,k,j)*((ngrid%c1h(k)*grid%mub(i,j)+ngrid%c2h(k)) + (ngrid%c1h(k)*grid%mu_2(i,j)))
               moist(i,k,j,P_QV) = moist(i,k,j,P_QV)*((ngrid%c1h(k)*grid%mub(i,j)+ngrid%c2h(k)) + (ngrid%c1h(k)*grid%mu_2(i,j)))
               END DO
            END DO
         END DO


      END IF
                               

#include "HALO_FORCE_DOWN.inc"

      ! code here to interpolate the data into the nested domain
# include "nest_forcedown_interp.inc"

      RETURN
   END SUBROUTINE force_domain_em_part2

!------------------------------------------------------------------

   SUBROUTINE interp_domain_em_part1 ( grid, intermediate_grid, ngrid, config_flags    &
!
#include "dummy_new_args.inc"
!
                 )
      USE module_state_description
      USE module_domain, ONLY : domain, get_ijk_from_grid
      USE module_configure, ONLY : grid_config_rec_type
      USE module_dm, ONLY : ntasks, ntasks_x, ntasks_y, itrace, local_communicator, &
                            nest_task_offsets, nest_pes_x, nest_pes_y, which_kid,   &
                            intercomm_active, mpi_comm_to_kid, mpi_comm_to_mom,     &
                            mytask, get_dm_max_halo_width
      USE module_timing
      IMPLICIT NONE
!
      TYPE(domain), POINTER :: grid          ! name of the grid being dereferenced (must be "grid")
      TYPE(domain), POINTER :: intermediate_grid
      TYPE(domain), POINTER :: ngrid
#include "dummy_new_decl.inc"
      INTEGER nlev, msize
      INTEGER i,j,pig,pjg,cm,cn,nig,njg,retval,k
      INTEGER iparstrt,jparstrt,sw
      TYPE (grid_config_rec_type)            :: config_flags
      REAL xv(2000)
      INTEGER       ::          cids, cide, cjds, cjde, ckds, ckde,    &
                                cims, cime, cjms, cjme, ckms, ckme,    &
                                cips, cipe, cjps, cjpe, ckps, ckpe
      INTEGER       ::          iids, iide, ijds, ijde, ikds, ikde,    &
                                iims, iime, ijms, ijme, ikms, ikme,    &
                                iips, iipe, ijps, ijpe, ikps, ikpe
      INTEGER       ::          nids, nide, njds, njde, nkds, nkde,    &
                                nims, nime, njms, njme, nkms, nkme,    &
                                nips, nipe, njps, njpe, nkps, nkpe

      INTEGER idim1,idim2,idim3,idim4,idim5,idim6,idim7

      INTEGER icoord, jcoord, idim_cd, jdim_cd, pgr
      INTEGER thisdomain_max_halo_width
      INTEGER local_comm, myproc, nproc
      INTEGER ioffset, ierr

      CALL wrf_get_dm_communicator ( local_comm )
      CALL wrf_get_myproc( myproc )
      CALL wrf_get_nproc( nproc )

      CALL get_ijk_from_grid (  grid ,                   &
                                cids, cide, cjds, cjde, ckds, ckde,    &
                                cims, cime, cjms, cjme, ckms, ckme,    &
                                cips, cipe, cjps, cjpe, ckps, ckpe    )
      CALL get_ijk_from_grid (  intermediate_grid ,              &
                                iids, iide, ijds, ijde, ikds, ikde,    &
                                iims, iime, ijms, ijme, ikms, ikme,    &
                                iips, iipe, ijps, ijpe, ikps, ikpe    )
      CALL get_ijk_from_grid (  ngrid ,              &
                                nids, nide, njds, njde, nkds, nkde,    &
                                nims, nime, njms, njme, nkms, nkme,    &
                                nips, nipe, njps, njpe, nkps, nkpe    )

      CALL nl_get_parent_grid_ratio ( ngrid%id, pgr )
      CALL nl_get_i_parent_start ( intermediate_grid%id, iparstrt )
      CALL nl_get_j_parent_start ( intermediate_grid%id, jparstrt )
      CALL nl_get_shw            ( intermediate_grid%id, sw )
      icoord =    iparstrt - sw
      jcoord =    jparstrt - sw
      idim_cd = iide - iids + 1
      jdim_cd = ijde - ijds + 1

      nlev  = ckde - ckds + 1

      ! get max_halo_width for parent. It may be smaller if it is moad
      CALL get_dm_max_halo_width ( grid%id , thisdomain_max_halo_width )

      IF ( grid%active_this_task ) THEN
#include "nest_interpdown_pack.inc"
      END IF

      ! determine which communicator and offset to use
      IF ( intercomm_active( grid%id ) ) THEN        ! I am parent
        local_comm = mpi_comm_to_kid( which_kid(ngrid%id), grid%id )
        ioffset = nest_task_offsets(ngrid%id)
      ELSE IF ( intercomm_active( ngrid%id ) ) THEN  ! I am nest
        local_comm = mpi_comm_to_mom( ngrid%id )
        ioffset = nest_task_offsets(ngrid%id)
      END IF

      IF ( grid%active_this_task .OR. ngrid%active_this_task ) THEN
#ifndef STUBMPI
        CALL mpi_comm_rank(local_comm,myproc,ierr)
        CALL mpi_comm_size(local_comm,nproc,ierr)
#endif
        CALL rsl_lite_bcast_msgs( myproc, nest_pes_x(grid%id)*nest_pes_y(grid%id),         &
                                          nest_pes_x(ngrid%id)*nest_pes_y(ngrid%id),       &
                                          ioffset, local_comm )
      END IF

      RETURN
   END SUBROUTINE interp_domain_em_part1

!------------------------------------------------------------------

   SUBROUTINE interp_domain_em_part2 ( grid, ngrid, pgrid, config_flags    &
!
#include "dummy_new_args.inc"
!
                 )
      USE module_state_description
      USE module_domain, ONLY : domain, get_ijk_from_grid
      USE module_configure, ONLY : grid_config_rec_type
      USE module_dm, ONLY : ntasks, ntasks_x, ntasks_y, itrace, local_communicator, &
                            mytask, get_dm_max_halo_width, which_kid
                            ! push_communicators_for_domain,pop_communicators_for_domain
      USE module_comm_nesting_dm, ONLY : halo_interp_down_sub
      IMPLICIT NONE
!
      TYPE(domain), POINTER :: grid          ! name of the grid being dereferenced (must be "grid")
      TYPE(domain), POINTER :: ngrid
      TYPE(domain), POINTER :: pgrid         !KAL added for vertical nesting
#include "dummy_new_decl.inc"
      INTEGER nlev, msize
      INTEGER i,j,pig,pjg,cm,cn,nig,njg,retval,k
      TYPE (grid_config_rec_type)            :: config_flags
      REAL xv(2000)
      INTEGER       ::          cids, cide, cjds, cjde, ckds, ckde,    &
                                cims, cime, cjms, cjme, ckms, ckme,    &
                                cips, cipe, cjps, cjpe, ckps, ckpe
      INTEGER       ::          nids, nide, njds, njde, nkds, nkde,    &
                                nims, nime, njms, njme, nkms, nkme,    &
                                nips, nipe, njps, njpe, nkps, nkpe
      INTEGER       ::          ids, ide, jds, jde, kds, kde,    &
                                ims, ime, jms, jme, kms, kme,    &
                                ips, ipe, jps, jpe, kps, kpe

      INTEGER idim1,idim2,idim3,idim4,idim5,idim6,idim7

      INTEGER myproc
      INTEGER ierr
      INTEGER thisdomain_max_halo_width

      !KAL variables for vertical nesting
      REAL :: p_top_m  , p_surf_m , mu_m , hsca_m , pre_c ,pre_n
      REAL, DIMENSION(pgrid%s_vert:pgrid%e_vert) :: alt_w_c
      REAL, DIMENSION(pgrid%s_vert:pgrid%e_vert+1) :: alt_u_c
      REAL, DIMENSION(ngrid%s_vert:ngrid%e_vert) :: alt_w_n
      REAL, DIMENSION(ngrid%s_vert:ngrid%e_vert+1) :: alt_u_n


      !KAL change this for vertical nesting
      ! interp_domain_em_part1 packs up the interpolation onto the coarse (vertical) grid
      ! therefore the message size is based on the coarse grid number of levels
      ! here it is unpacked onto the intermediate grid
       CALL get_ijk_from_grid ( pgrid ,                   &
                                cids, cide, cjds, cjde, ckds, ckde,    &
                                cims, cime, cjms, cjme, ckms, ckme,    &
                                cips, cipe, cjps, cjpe, ckps, ckpe    )
      !KAL this is the original WRF code
      !CALL get_ijk_from_grid (  grid ,                   &
      !                          cids, cide, cjds, cjde, ckds, ckde,    &
      !                          cims, cime, cjms, cjme, ckms, ckme,    &
      !                          cips, cipe, cjps, cjpe, ckps, ckpe    )
      CALL get_ijk_from_grid (  ngrid ,              &
                                nids, nide, njds, njde, nkds, nkde,    &
                                nims, nime, njms, njme, nkms, nkme,    &
                                nips, nipe, njps, njpe, nkps, nkpe    )

      nlev  = ckde - ckds + 1

      CALL get_dm_max_halo_width ( ngrid%id , thisdomain_max_halo_width )

#include "nest_interpdown_unpack.inc"


if (ngrid%vert_refine_method .NE. 0) then

      !KAL calculating the vertical coordinate for parent and nest grid (code from ndown)
      ! assume that the parent and nest have the same p_top value (as in ndown)
      
!KAL ckde is equal to e_vert of the coarse grid. There are e_vert-1 u points.  The coarse 1D grid here is e_vert+1,
!    so it is the e_vert-1 points from the coarse grid, plus a surface point plus a top point.  Extrapolation coefficients
!    are used to get the surface and top points to fill out the pro_u_c 1D array of u values from the coarse grid.     
                
      hsca_m = 6.7 !KAL scale height of the atmosphere
      p_top_m = ngrid%p_top
      p_surf_m = 1.e5
      mu_m = p_surf_m - p_top_m
!    parent
      do  k = 1,ckde
      pre_c = mu_m * pgrid%c3f(k) + p_top_m + pgrid%c4f(k)
      alt_w_c(k) =  -hsca_m * alog(pre_c/p_surf_m)
      enddo   
      do  k = 1,ckde-1
      pre_c = mu_m * pgrid%c3h(k) + p_top_m + pgrid%c4h(k)
      alt_u_c(k+1) =  -hsca_m * alog(pre_c/p_surf_m)
      enddo
      alt_u_c(1) =  alt_w_c(1)
      alt_u_c(ckde+1) =  alt_w_c(ckde)       
!    nest
      do  k = 1,nkde
      pre_n = mu_m * ngrid%c3f(k) + p_top_m + ngrid%c4f(k)
      alt_w_n(k) =  -hsca_m * alog(pre_n/p_surf_m)
      enddo
      do  k = 1,nkde-1
      pre_n = mu_m * ngrid%c3h(k) + p_top_m + ngrid%c4h(k)
      alt_u_n(k+1) =  -hsca_m * alog(pre_n/p_surf_m)
      enddo
      alt_u_n(1) =  alt_w_n(1)
      alt_u_n(nkde+1) =  alt_w_n(nkde)
endif   



      !KAL added this call for vertical nesting (return coarse grid dimensions to intended values)
      CALL get_ijk_from_grid (  grid ,                   &
                                cids, cide, cjds, cjde, ckds, ckde,    &
                                cims, cime, cjms, cjme, ckms, ckme,    &
                                cips, cipe, cjps, cjpe, ckps, ckpe    )

      CALL get_ijk_from_grid (  grid ,              &
                                ids, ide, jds, jde, kds, kde,    &
                                ims, ime, jms, jme, kms, kme,    &
                                ips, ipe, jps, jpe, kps, kpe    )


if (ngrid%vert_refine_method .NE. 0) then
      
!KAL added this code (the include file) for the vertical nesting
#include "nest_interpdown_interp_vert.inc"


      !KAL finish off the 1-D variables (t_base, u_base, v_base, qv_base, and z_base) (move this out of here if alt_u_c and alt_u_n are calculated elsewhere)
      CALL vert_interp_vert_nesting_1d ( &         
                                        ngrid%t_base,                                           &    ! CD field
                                        ids, ide, kds, kde, jds, jde,                           &    ! CD dims
                                        ims, ime, kms, kme, jms, jme,                           &    ! CD dims
                                        ips, ipe, kps, MIN( (kde-1), kpe ), jps, jpe,           &    ! CD dims
                                        pgrid%s_vert, pgrid%e_vert,                             &    ! vertical dimension of the parent grid
                                        pgrid%cf1, pgrid%cf2, pgrid%cf3, pgrid%cfn, pgrid%cfn1, &    ! coarse grid extrapolation constants
                                        alt_u_c, alt_u_n)                                            ! coordinates for parent and nest
      CALL vert_interp_vert_nesting_1d ( &         
                                        ngrid%u_base,                                           &    ! CD field
                                        ids, ide, kds, kde, jds, jde,                           &    ! CD dims
                                        ims, ime, kms, kme, jms, jme,                           &    ! CD dims
                                        ips, ipe, kps, MIN( (kde-1), kpe ), jps, jpe,           &    ! CD dims
                                        pgrid%s_vert, pgrid%e_vert,                             &    ! vertical dimension of the parent grid
                                        pgrid%cf1, pgrid%cf2, pgrid%cf3, pgrid%cfn, pgrid%cfn1, &    ! coarse grid extrapolation constants
                                        alt_u_c, alt_u_n)                                            ! coordinates for parent and nest
      CALL vert_interp_vert_nesting_1d ( &         
                                        ngrid%v_base,                                           &    ! CD field
                                        ids, ide, kds, kde, jds, jde,                           &    ! CD dims
                                        ims, ime, kms, kme, jms, jme,                           &    ! CD dims
                                        ips, ipe, kps, MIN( (kde-1), kpe ), jps, jpe,           &    ! CD dims
                                        pgrid%s_vert, pgrid%e_vert,                             &    ! vertical dimension of the parent grid
                                        pgrid%cf1, pgrid%cf2, pgrid%cf3, pgrid%cfn, pgrid%cfn1, &    ! coarse grid extrapolation constants
                                        alt_u_c, alt_u_n)                                            ! coordinates for parent and nest
      CALL vert_interp_vert_nesting_1d ( &         
                                        ngrid%qv_base,                                          &    ! CD field
                                        ids, ide, kds, kde, jds, jde,                           &    ! CD dims
                                        ims, ime, kms, kme, jms, jme,                           &    ! CD dims
                                        ips, ipe, kps, MIN( (kde-1), kpe ), jps, jpe,           &    ! CD dims
                                        pgrid%s_vert, pgrid%e_vert,                             &    ! vertical dimension of the parent grid
                                        pgrid%cf1, pgrid%cf2, pgrid%cf3, pgrid%cfn, pgrid%cfn1, &    ! coarse grid extrapolation constants
                                        alt_u_c, alt_u_n)                                            ! coordinates for parent and nest
      CALL vert_interp_vert_nesting_1d ( &         
                                        ngrid%z_base,                                           &    ! CD field
                                        ids, ide, kds, kde, jds, jde,                           &    ! CD dims
                                        ims, ime, kms, kme, jms, jme,                           &    ! CD dims
                                        ips, ipe, kps, MIN( (kde-1), kpe ), jps, jpe,           &    ! CD dims
                                        pgrid%s_vert, pgrid%e_vert,                             &    ! vertical dimension of the parent grid
                                        pgrid%cf1, pgrid%cf2, pgrid%cf3, pgrid%cfn, pgrid%cfn1, &    ! coarse grid extrapolation constants
                                        alt_u_c, alt_u_n)                                            ! coordinates for parent and nest

endif
        
        CALL push_communicators_for_domain( grid%id )

#include "HALO_INTERP_DOWN.inc"

        CALL pop_communicators_for_domain

# include "nest_interpdown_interp.inc"

      RETURN
   END SUBROUTINE interp_domain_em_part2

!------------------------------------------------------------------

   SUBROUTINE interp_domain_em_small_part1 ( grid, intermediate_grid, ngrid, config_flags    &
!
#include "dummy_new_args.inc"
!
                 )
      USE module_state_description
      USE module_domain, ONLY : domain, get_ijk_from_grid
      USE module_configure, ONLY : grid_config_rec_type
      USE module_comm_dm, ONLY: halo_em_horiz_interp_sub
      USE module_dm, ONLY : ntasks, ntasks_x, ntasks_y, itrace, local_communicator, &
                            mytask, get_dm_max_halo_width,                          &
                            nest_task_offsets, mpi_comm_to_kid, mpi_comm_to_mom,    &
                            which_kid, nest_pes_x, nest_pes_y, intercomm_active
      USE module_timing
      IMPLICIT NONE
!
      TYPE(domain), POINTER :: grid          ! name of the grid being dereferenced (must be "grid")
      TYPE(domain), POINTER :: intermediate_grid
      TYPE(domain), POINTER :: ngrid
#include "dummy_new_decl.inc"
      INTEGER nlev, msize
      INTEGER i,j,pig,pjg,cm,cn,nig,njg,retval,k
      INTEGER iparstrt,jparstrt,sw
      TYPE (grid_config_rec_type)            :: config_flags
      REAL xv(2000)
      INTEGER       ::           ids,  ide,  jds,  jde,  kds,  kde,    &
                                 ims,  ime,  jms,  jme,  kms,  kme,    &
                                 ips,  ipe,  jps,  jpe,  kps,  kpe

      INTEGER       ::          cids, cide, cjds, cjde, ckds, ckde,    &
                                cims, cime, cjms, cjme, ckms, ckme,    &
                                cips, cipe, cjps, cjpe, ckps, ckpe
      INTEGER       ::          iids, iide, ijds, ijde, ikds, ikde,    &
                                iims, iime, ijms, ijme, ikms, ikme,    &
                                iips, iipe, ijps, ijpe, ikps, ikpe
      INTEGER       ::          nids, nide, njds, njde, nkds, nkde,    &
                                nims, nime, njms, njme, nkms, nkme,    &
                                nips, nipe, njps, njpe, nkps, nkpe

      INTEGER idim1,idim2,idim3,idim4,idim5,idim6,idim7

      INTEGER icoord, jcoord, idim_cd, jdim_cd, pgr
      INTEGER thisdomain_max_halo_width
      INTEGER local_comm, myproc, nproc
      INTEGER ioffset

      CALL wrf_get_dm_communicator ( local_comm )
      CALL wrf_get_myproc( myproc )
      CALL wrf_get_nproc( nproc )

      CALL get_ijk_from_grid (  grid ,                           &
                                ids, ide, jds, jde, kds, kde,    &
                                ims, ime, jms, jme, kms, kme,    &
                                ips, ipe, jps, jpe, kps, kpe     )
#ifdef DM_PARALLEL
# include "HALO_EM_HORIZ_INTERP.inc"
#endif

      CALL get_ijk_from_grid (  grid ,                   &
                                cids, cide, cjds, cjde, ckds, ckde,    &
                                cims, cime, cjms, cjme, ckms, ckme,    &
                                cips, cipe, cjps, cjpe, ckps, ckpe    )
      CALL get_ijk_from_grid (  intermediate_grid ,              &
                                iids, iide, ijds, ijde, ikds, ikde,    &
                                iims, iime, ijms, ijme, ikms, ikme,    &
                                iips, iipe, ijps, ijpe, ikps, ikpe    )
      CALL get_ijk_from_grid (  ngrid ,              &
                                nids, nide, njds, njde, nkds, nkde,    &
                                nims, nime, njms, njme, nkms, nkme,    &
                                nips, nipe, njps, njpe, nkps, nkpe    )

      CALL nl_get_parent_grid_ratio ( ngrid%id, pgr )
      CALL nl_get_i_parent_start ( intermediate_grid%id, iparstrt )
      CALL nl_get_j_parent_start ( intermediate_grid%id, jparstrt )
      CALL nl_get_shw            ( intermediate_grid%id, sw )
      icoord =    iparstrt - sw
      jcoord =    jparstrt - sw
      idim_cd = iide - iids + 1
      jdim_cd = ijde - ijds + 1

      nlev  = ckde - ckds + 1

      ! get max_halo_width for parent. It may be smaller if it is moad
      CALL get_dm_max_halo_width ( grid%id , thisdomain_max_halo_width )

      !  How many 3d arrays, so far just 3d theta-300 and geopotential perturbation,
      !  and the 2d topo elevation, three max press/temp/height fields, and three
      !  min press/temp/height fields.
   
      msize = ( 2 )* nlev + 7
   
!call wrf_debug(0,'/external/RSL_LITE/module_dm.F, calling rsl_lite_to_child')
      CALL rsl_lite_to_child_info( local_communicator, msize*RWORDSIZE     &
                              ,cips,cipe,cjps,cjpe                         &
                              ,iids,iide,ijds,ijde                         &
                              ,nids,nide,njds,njde                         &
                              ,pgr , sw                                    &
                              ,ntasks_x,ntasks_y                           &
                              ,thisdomain_max_halo_width                   &
                              ,icoord,jcoord                               &
                              ,idim_cd,jdim_cd                             &
                              ,pig,pjg,retval )
!call wrf_debug(0,'/external/RSL_LITE/module_dm.F, back from rsl_lite_to_child')
      DO while ( retval .eq. 1 )
         IF ( SIZE(grid%ph_2) .GT. 1 ) THEN
!call wrf_debug(0,'/external/RSL_LITE/module_dm.F, ph_2')
            DO k = ckds,ckde
               xv(k)= grid%ph_2(pig,k,pjg)
            END DO
            CALL rsl_lite_to_child_msg(((ckde)-(ckds)+1)*RWORDSIZE,xv)
         END IF
   
         IF ( SIZE(grid%t_2) .GT. 1 ) THEN
!call wrf_debug(0,'/external/RSL_LITE/module_dm.F, t_2')
            DO k = ckds,(ckde-1)
               xv(k)= grid%t_2(pig,k,pjg)
            END DO
            CALL rsl_lite_to_child_msg((((ckde-1))-(ckds)+1)*RWORDSIZE,xv)
         END IF
   
         IF ( SIZE(grid%ht) .GT. 1 ) THEN
!call wrf_debug(0,'/external/RSL_LITE/module_dm.F, ht')
            xv(1)= grid%ht(pig,pjg)
            CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
         END IF
   
         IF ( SIZE(grid%t_max_p) .GT. 1 ) THEN
!call wrf_debug(0,'/external/RSL_LITE/module_dm.F, t_max_p')
            xv(1)= grid%t_max_p(pig,pjg)
            CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
         END IF
   
         IF ( SIZE(grid%ght_max_p) .GT. 1 ) THEN
!call wrf_debug(0,'/external/RSL_LITE/module_dm.F, ght_max_p')
            xv(1)= grid%ght_max_p(pig,pjg)
            CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
         END IF
   
         IF ( SIZE(grid%max_p) .GT. 1 ) THEN
!call wrf_debug(0,'/external/RSL_LITE/module_dm.F, max_p')
            xv(1)= grid%max_p(pig,pjg)
            CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
         END IF
   
         IF ( SIZE(grid%t_min_p) .GT. 1 ) THEN
!call wrf_debug(0,'/external/RSL_LITE/module_dm.F, t_min_p')
            xv(1)= grid%t_min_p(pig,pjg)
            CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
         END IF
   
         IF ( SIZE(grid%ght_min_p) .GT. 1 ) THEN
!call wrf_debug(0,'/external/RSL_LITE/module_dm.F, ght_min_p')
            xv(1)= grid%ght_min_p(pig,pjg)
            CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
         END IF
   
         IF ( SIZE(grid%min_p) .GT. 1 ) THEN
!call wrf_debug(0,'/external/RSL_LITE/module_dm.F, min_p')
            xv(1)= grid%min_p(pig,pjg)
            CALL rsl_lite_to_child_msg(RWORDSIZE,xv)
         END IF
   
!call wrf_debug(0,'/external/RSL_LITE/module_dm.F, calling rsl_lite_to_child_info')
         CALL rsl_lite_to_child_info( local_communicator, msize*RWORDSIZE  &
                                     ,cips,cipe,cjps,cjpe                  &
                                     ,iids,iide,ijds,ijde                  &
                                     ,nids,nide,njds,njde                  &
                                     ,pgr , sw                             &
                                     ,ntasks_x,ntasks_y                    &
                                     ,thisdomain_max_halo_width            &
                                     ,icoord,jcoord                        &
                                     ,idim_cd,jdim_cd                      &
                                     ,pig,pjg,retval )
!call wrf_debug(0,'/external/RSL_LITE/module_dm.F, back from rsl_lite_to_child_info')
      END DO

      ! determine which communicator and offset to use
      IF ( intercomm_active( grid%id ) ) THEN        ! I am parent
        local_comm = mpi_comm_to_kid( which_kid(ngrid%id), grid%id )
        ioffset = nest_task_offsets(ngrid%id)
      ELSE IF ( intercomm_active( ngrid%id ) ) THEN  ! I am nest
        local_comm = mpi_comm_to_mom( ngrid%id )
        ioffset = nest_task_offsets(ngrid%id)
      END IF

!call wrf_debug(0,'/external/RSL_LITE/module_dm.F, calling rsl_lite_bcast')
      CALL rsl_lite_bcast_msgs( myproc, nest_pes_x(grid%id)*nest_pes_y(grid%id),         &
                                        nest_pes_x(ngrid%id)*nest_pes_y(ngrid%id),       &
                                        ioffset, local_comm )
!call wrf_debug(0,'/external/RSL_LITE/module_dm.F, back from rsl_lite_bcast')

      RETURN
   END SUBROUTINE interp_domain_em_small_part1

!------------------------------------------------------------------

   SUBROUTINE interp_domain_em_small_part2 ( grid, ngrid, config_flags    &
!
#include "dummy_new_args.inc"
!
                 )
      USE module_state_description
      USE module_domain, ONLY : domain, get_ijk_from_grid
      USE module_configure, ONLY : grid_config_rec_type
      USE module_dm, ONLY : ntasks, ntasks_x, ntasks_y, itrace, local_communicator, &
                            mytask, get_dm_max_halo_width
      USE module_comm_nesting_dm, ONLY : halo_interp_down_sub
      IMPLICIT NONE
!
      TYPE(domain), POINTER :: grid          ! name of the grid being dereferenced (must be "grid")
      TYPE(domain), POINTER :: ngrid
#include "dummy_new_decl.inc"
      INTEGER nlev, msize
      INTEGER i,j,pig,pjg,cm,cn,nig,njg,retval,k
      TYPE (grid_config_rec_type)            :: config_flags
      REAL xv(2000)
      INTEGER       ::          cids, cide, cjds, cjde, ckds, ckde,    &
                                cims, cime, cjms, cjme, ckms, ckme,    &
                                cips, cipe, cjps, cjpe, ckps, ckpe
      INTEGER       ::          nids, nide, njds, njde, nkds, nkde,    &
                                nims, nime, njms, njme, nkms, nkme,    &
                                nips, nipe, njps, njpe, nkps, nkpe
      INTEGER       ::          ids, ide, jds, jde, kds, kde,    &
                                ims, ime, jms, jme, kms, kme,    &
                                ips, ipe, jps, jpe, kps, kpe

      INTEGER idim1,idim2,idim3,idim4,idim5,idim6,idim7

      INTEGER myproc
      INTEGER ierr
      INTEGER thisdomain_max_halo_width

      CALL get_ijk_from_grid (  grid ,                   &
                                cids, cide, cjds, cjde, ckds, ckde,    &
                                cims, cime, cjms, cjme, ckms, ckme,    &
                                cips, cipe, cjps, cjpe, ckps, ckpe    )
      CALL get_ijk_from_grid (  ngrid ,              &
                                nids, nide, njds, njde, nkds, nkde,    &
                                nims, nime, njms, njme, nkms, nkme,    &
                                nips, nipe, njps, njpe, nkps, nkpe    )

      nlev  = ckde - ckds + 1

      CALL get_dm_max_halo_width ( ngrid%id , thisdomain_max_halo_width )

      CALL rsl_lite_from_parent_info(pig,pjg,retval)
      
      DO while ( retval .eq. 1 )
      
         IF ( SIZE(grid%ph_2) .GT. 1 ) THEN
            CALL rsl_lite_from_parent_msg(((ckde)-(ckds)+1)*RWORDSIZE,xv)
            DO k = ckds,ckde
               grid%ph_2(pig,k,pjg) = xv(k)
            END DO
         END IF
   
         IF ( SIZE(grid%t_2) .GT. 1 ) THEN
            CALL rsl_lite_from_parent_msg((((ckde-1))-(ckds)+1)*RWORDSIZE,xv)
            DO k = ckds,(ckde-1)
               grid%t_2(pig,k,pjg) = xv(k)
            END DO
         END IF
   
         IF ( SIZE(grid%ht) .GT. 1 ) THEN
            CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
            grid%ht(pig,pjg) = xv(1)
         END IF
   
         IF ( SIZE(grid%t_max_p) .GT. 1 ) THEN
            CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
            grid%t_max_p(pig,pjg) = xv(1)
         END IF
   
         IF ( SIZE(grid%ght_max_p) .GT. 1 ) THEN
            CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
            grid%ght_max_p(pig,pjg) = xv(1)
         END IF
   
         IF ( SIZE(grid%max_p) .GT. 1 ) THEN
            CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
            grid%max_p(pig,pjg) = xv(1)
         END IF
   
         IF ( SIZE(grid%t_min_p) .GT. 1 ) THEN
            CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
            grid%t_min_p(pig,pjg) = xv(1)
         END IF
   
         IF ( SIZE(grid%ght_min_p) .GT. 1 ) THEN
            CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
            grid%ght_min_p(pig,pjg) = xv(1)
         END IF
   
         IF ( SIZE(grid%min_p) .GT. 1 ) THEN
            CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
            grid%min_p(pig,pjg) = xv(1)
         END IF
      
         CALL rsl_lite_from_parent_info(pig,pjg,retval)
         
      END DO

      CALL get_ijk_from_grid (  grid ,              &
                                ids, ide, jds, jde, kds, kde,    &
                                ims, ime, jms, jme, kms, kme,    &
                                ips, ipe, jps, jpe, kps, kpe    )

#include "HALO_INTERP_DOWN.inc"

      CALL interp_fcn_bl ( grid%ph_2,                                           &       
                           cids, cide, ckds, ckde, cjds, cjde,                  &         
                           cims, cime, ckms, ckme, cjms, cjme,                  &         
                           cips, cipe, ckps, MIN( ckde, ckpe ), cjps, cjpe,     &         
                           ngrid%ph_2,                                          &   
                           nids, nide, nkds, nkde, njds, njde,                  &         
                           nims, nime, nkms, nkme, njms, njme,                  &         
                           nips, nipe, nkps, MIN( nkde, nkpe ), njps, njpe,     &         
                           config_flags%shw, ngrid%imask_nostag,                &         
                           .FALSE., .FALSE.,                                    &         
                           ngrid%i_parent_start, ngrid%j_parent_start,          &
                           ngrid%parent_grid_ratio, ngrid%parent_grid_ratio,    &
                           grid%ht, ngrid%ht,                                   &
                           grid%t_max_p, ngrid%t_max_p,                         &
                           grid%ght_max_p, ngrid%ght_max_p,                     &
                           grid%max_p, ngrid%max_p,                             &
                           grid%t_min_p, ngrid%t_min_p,                         &
                           grid%ght_min_p, ngrid%ght_min_p,                     &
                           grid%min_p, ngrid%min_p,                             &
                           ngrid%znw, ngrid%p_top                               )
      
      CALL interp_fcn_bl ( grid%t_2,                                            &       
                           cids, cide, ckds, ckde, cjds, cjde,                  &         
                           cims, cime, ckms, ckme, cjms, cjme,                  &         
                           cips, cipe, ckps, MIN( (ckde-1), ckpe ), cjps, cjpe, &         
                           ngrid%t_2,                                           &   
                           nids, nide, nkds, nkde, njds, njde,                  &         
                           nims, nime, nkms, nkme, njms, njme,                  &         
                           nips, nipe, nkps, MIN( (nkde-1), nkpe ), njps, njpe, &         
                           config_flags%shw, ngrid%imask_nostag,                &         
                           .FALSE., .FALSE.,                                    &         
                           ngrid%i_parent_start, ngrid%j_parent_start,          &
                           ngrid%parent_grid_ratio, ngrid%parent_grid_ratio,    &
                           grid%ht, ngrid%ht,                                   &
                           grid%t_max_p, ngrid%t_max_p,                         &
                           grid%ght_max_p, ngrid%ght_max_p,                     &
                           grid%max_p, ngrid%max_p,                             &
                           grid%t_min_p, ngrid%t_min_p,                         &
                           grid%ght_min_p, ngrid%ght_min_p,                     &
                           grid%min_p, ngrid%min_p,                             &
                           ngrid%znu, ngrid%p_top                               )

      RETURN
   END SUBROUTINE interp_domain_em_small_part2

!------------------------------------------------------------------

   SUBROUTINE feedback_nest_prep ( grid, config_flags    &
!
#include "dummy_new_args.inc"
!
)
      USE module_state_description
      USE module_domain, ONLY : domain, get_ijk_from_grid
      USE module_configure, ONLY : grid_config_rec_type
      USE module_dm, ONLY : ntasks, ntasks_x, ntasks_y, itrace, local_communicator, mytask !, &
                             !push_communicators_for_domain, pop_communicators_for_domain
      USE module_comm_nesting_dm, ONLY : halo_interp_up_sub
      IMPLICIT NONE
!
      TYPE(domain), TARGET :: grid          ! name of the grid being dereferenced (must be "grid")
      TYPE (grid_config_rec_type) :: config_flags ! configureation flags, has vertical dim of
                                                  ! soil temp, moisture, etc., has vertical dim
                                                  ! of soil categories
#include "dummy_new_decl.inc"

      INTEGER       ::          ids, ide, jds, jde, kds, kde,    &
                                ims, ime, jms, jme, kms, kme,    &
                                ips, ipe, jps, jpe, kps, kpe

      INTEGER idim1,idim2,idim3,idim4,idim5,idim6,idim7

      INTEGER       :: idum1, idum2


      CALL get_ijk_from_grid (  grid ,              &
                                ids, ide, jds, jde, kds, kde,    &
                                ims, ime, jms, jme, kms, kme,    &
                                ips, ipe, jps, jpe, kps, kpe    )

    IF ( grid%active_this_task ) THEN
      CALL push_communicators_for_domain( grid%id )

#ifdef DM_PARALLEL
#include "HALO_INTERP_UP.inc"
#endif

      CALL pop_communicators_for_domain
    END IF

   END SUBROUTINE feedback_nest_prep

!------------------------------------------------------------------

   SUBROUTINE feedback_domain_em_part1 ( grid, ngrid, config_flags    &
!
#include "dummy_new_args.inc"
!
                 )
      USE module_state_description
      USE module_domain, ONLY : domain, get_ijk_from_grid
      USE module_configure, ONLY : grid_config_rec_type, model_config_rec, model_to_grid_config_rec, set_scalar_indices_from_config
      USE module_dm, ONLY : ntasks, ntasks_x, ntasks_y, itrace, local_communicator, mytask, &
                            ipe_save, jpe_save, ips_save, jps_save,                         &
                            nest_pes_x, nest_pes_y

      IMPLICIT NONE
!
      TYPE(domain), POINTER :: grid          ! name of the grid being dereferenced (must be "grid")
      TYPE(domain), POINTER :: ngrid
#include "dummy_new_decl.inc"
      INTEGER nlev, msize
      INTEGER i,j,pig,pjg,cm,cn,nig,njg,retval,k
      TYPE(domain), POINTER :: xgrid
      TYPE (grid_config_rec_type)            :: config_flags, nconfig_flags
      REAL xv(2000)
      INTEGER       ::          cids, cide, cjds, cjde, ckds, ckde,    &
                                cims, cime, cjms, cjme, ckms, ckme,    &
                                cips, cipe, cjps, cjpe, ckps, ckpe
      INTEGER       ::          nids, nide, njds, njde, nkds, nkde,    &
                                nims, nime, njms, njme, nkms, nkme,    &
                                nips, nipe, njps, njpe, nkps, nkpe

      INTEGER idim1,idim2,idim3,idim4,idim5,idim6,idim7

      INTEGER local_comm, myproc, nproc, idum1, idum2
      INTEGER thisdomain_max_halo_width

!cyl: add variables for trajectory
      integer tjk

      INTERFACE
          SUBROUTINE feedback_nest_prep ( grid, config_flags    &
!
#include "dummy_new_args.inc"
!
)
             USE module_state_description
             USE module_domain, ONLY : domain
             USE module_configure, ONLY : grid_config_rec_type
!
             TYPE (grid_config_rec_type)            :: config_flags
             TYPE(domain), TARGET                   :: grid
#include "dummy_new_decl.inc"
          END SUBROUTINE feedback_nest_prep
      END INTERFACE
!

      CALL wrf_get_dm_communicator ( local_comm )
      CALL wrf_get_myproc( myproc )
      CALL wrf_get_nproc( nproc )

!
! intermediate grid
      CALL get_ijk_from_grid (  grid ,                                 &
                                cids, cide, cjds, cjde, ckds, ckde,    &
                                cims, cime, cjms, cjme, ckms, ckme,    &
                                cips, cipe, cjps, cjpe, ckps, ckpe    )
! nest grid
      CALL get_ijk_from_grid (  ngrid ,                                &
                                nids, nide, njds, njde, nkds, nkde,    &
                                nims, nime, njms, njme, nkms, nkme,    &
                                nips, nipe, njps, njpe, nkps, nkpe    )

      nlev  = ckde - ckds + 1

      ips_save = ngrid%i_parent_start   ! used in feedback_domain_em_part2 below
      jps_save = ngrid%j_parent_start
      ipe_save = ngrid%i_parent_start + (nide-nids+1) / ngrid%parent_grid_ratio - 1
      jpe_save = ngrid%j_parent_start + (njde-njds+1) / ngrid%parent_grid_ratio - 1

! feedback_nest_prep invokes a halo exchange on the ngrid. It is done this way
! in a separate routine because the HALOs need the data to be dereference from the
! grid data structure and, in this routine, the dereferenced fields are related to
! the intermediate domain, not the nest itself. Save the current grid pointer to intermediate
! domain, switch grid to point to ngrid, invoke feedback_nest_prep,  then restore grid
! to point to intermediate domain.

      CALL model_to_grid_config_rec ( ngrid%id , model_config_rec , nconfig_flags )
      CALL set_scalar_indices_from_config ( ngrid%id , idum1 , idum2 )
      xgrid => grid
      grid => ngrid

      CALL feedback_nest_prep ( grid, nconfig_flags    &
!
#include "actual_new_args.inc"
!
)

! put things back so grid is intermediate grid

      grid => xgrid
      CALL set_scalar_indices_from_config ( grid%id , idum1 , idum2 )

! "interp" (basically copy) ngrid onto intermediate grid

#include "nest_feedbackup_interp.inc"

      RETURN
   END SUBROUTINE feedback_domain_em_part1

!------------------------------------------------------------------

   SUBROUTINE feedback_domain_em_part2 ( grid, intermediate_grid, ngrid , config_flags    &
!
#include "dummy_new_args.inc"
!
                 )
      USE module_state_description
      USE module_domain, ONLY : domain, domain_clock_get, get_ijk_from_grid
      USE module_configure, ONLY : grid_config_rec_type, model_config_rec
      USE module_dm, ONLY : ntasks, ntasks_x, ntasks_y, itrace, local_communicator, mytask, &
                            ipe_save, jpe_save, ips_save, jps_save, get_dm_max_halo_width,  &
                            nest_pes_x, nest_pes_y,                                         &
                            intercomm_active, nest_task_offsets,                    &
                            mpi_comm_to_mom, mpi_comm_to_kid, which_kid !,            &
                             !push_communicators_for_domain, pop_communicators_for_domain

      USE module_comm_nesting_dm, ONLY : halo_interp_up_sub
      USE module_utility
      IMPLICIT NONE

!
      TYPE(domain), POINTER :: grid          ! name of the grid being dereferenced (must be "grid")
      TYPE(domain), POINTER :: intermediate_grid
      TYPE(domain), POINTER :: ngrid
      TYPE(domain), POINTER :: parent_grid

#include "dummy_new_decl.inc"
      INTEGER nlev, msize
      INTEGER i,j,pig,pjg,cm,cn,nig,njg,retval,k
      TYPE (grid_config_rec_type)            :: config_flags
      REAL xv(2000)
      INTEGER       ::          cids, cide, cjds, cjde, ckds, ckde,    &
                                cims, cime, cjms, cjme, ckms, ckme,    &
                                cips, cipe, cjps, cjpe, ckps, ckpe
      INTEGER       ::          nids, nide, njds, njde, nkds, nkde,    &
                                nims, nime, njms, njme, nkms, nkme,    &
                                nips, nipe, njps, njpe, nkps, nkpe
      INTEGER       ::          xids, xide, xjds, xjde, xkds, xkde,    &
                                xims, xime, xjms, xjme, xkms, xkme,    &
                                xips, xipe, xjps, xjpe, xkps, xkpe
      INTEGER       ::          ids, ide, jds, jde, kds, kde,    &
                                ims, ime, jms, jme, kms, kme,    &
                                ips, ipe, jps, jpe, kps, kpe

      INTEGER idim1,idim2,idim3,idim4,idim5,idim6,idim7

      INTEGER icoord, jcoord, idim_cd, jdim_cd
      INTEGER local_comm, myproc, nproc, ioffset
      INTEGER iparstrt, jparstrt, sw, thisdomain_max_halo_width
      REAL    nest_influence

      character*256 :: timestr
      integer ierr

      LOGICAL, EXTERNAL  :: cd_feedback_mask

!cyl: add variables for trajectory
      integer tjk

! On entry to this routine,
!  "grid" refers to the parent domain
!  "intermediate_grid" refers to local copy of parent domain that overlies this patch of nest
!  "ngrid" refers to the nest, which is only needed for smoothing on the parent because
!          the nest feedback data has already been transferred during em_nest_feedbackup_interp
!          in part1, above.
! The way these settings c and n dimensions are set, below, looks backwards but from the point
! of view of the RSL routine rsl_lite_to_parent_info(), call to which is included by
! em_nest_feedbackup_pack, the "n" domain represents the parent domain and the "c" domain
! represents the intermediate domain. The backwards lookingness should be fixed in the gen_comms.c
! registry routine that accompanies RSL_LITE but, just as it's sometimes easier to put up a road
! sign that says "DIP" than fix the dip,  at this point it was easier just to write this comment. JM
!
      nest_influence = 1.

      CALL domain_clock_get( grid, current_timestr=timestr )

      CALL get_ijk_from_grid (  intermediate_grid ,                   &
                                cids, cide, cjds, cjde, ckds, ckde,    &
                                cims, cime, cjms, cjme, ckms, ckme,    &
                                cips, cipe, cjps, cjpe, ckps, ckpe    )
      CALL get_ijk_from_grid (  grid ,              &
                                nids, nide, njds, njde, nkds, nkde,    &
                                nims, nime, njms, njme, nkms, nkme,    &
                                nips, nipe, njps, njpe, nkps, nkpe    )
      CALL get_ijk_from_grid (  ngrid ,              &
                                xids, xide, xjds, xjde, xkds, xkde,    &
                                xims, xime, xjms, xjme, xkms, xkme,    &
                                xips, xipe, xjps, xjpe, xkps, xkpe    )

      ips_save = ngrid%i_parent_start   ! used in feedback_domain_em_part2 below
      jps_save = ngrid%j_parent_start
      ipe_save = ngrid%i_parent_start + (xide-xids+1) / ngrid%parent_grid_ratio - 1
      jpe_save = ngrid%j_parent_start + (xjde-xjds+1) / ngrid%parent_grid_ratio - 1




IF ( ngrid%active_this_task ) THEN
!cyl add this for trajectory
    CALL push_communicators_for_domain( ngrid%id )

    do tjk = 1,config_flags%num_traj
     if (ngrid%traj_long(tjk) .eq. -9999.0) then
!       print*,'n=-9999',tjk
        ngrid%traj_long(tjk)=grid%traj_long(tjk)
        ngrid%traj_k(tjk)=grid%traj_k(tjk)
     else
!       print*,'p!=-9999',tjk
        grid%traj_long(tjk)=ngrid%traj_long(tjk)
        grid%traj_k(tjk)=ngrid%traj_k(tjk)
     endif
     if (ngrid%traj_lat(tjk) .eq. -9999.0) then
         ngrid%traj_lat(tjk)=grid%traj_lat(tjk)
         ngrid%traj_k(tjk)=grid%traj_k(tjk)
     else
         grid%traj_lat(tjk)=ngrid%traj_lat(tjk)
         grid%traj_k(tjk)=ngrid%traj_k(tjk)
     endif
    enddo
!endcyl

      CALL nl_get_i_parent_start ( intermediate_grid%id, iparstrt )
      CALL nl_get_j_parent_start ( intermediate_grid%id, jparstrt )
      CALL nl_get_shw            ( intermediate_grid%id, sw )
      icoord =    iparstrt - sw
      jcoord =    jparstrt - sw
      idim_cd = cide - cids + 1
      jdim_cd = cjde - cjds + 1

      nlev  = ckde - ckds + 1

      CALL get_dm_max_halo_width ( grid%id , thisdomain_max_halo_width )

      parent_grid => grid
      grid => ngrid
#include "nest_feedbackup_pack.inc"
      grid => parent_grid
    CALL pop_communicators_for_domain

END IF

!      CALL wrf_get_dm_communicator ( local_comm )
!      CALL wrf_get_myproc( myproc )
!      CALL wrf_get_nproc( nproc )

      ! determine which communicator and offset to use
      IF ( intercomm_active( grid%id ) ) THEN        ! I am parent
        local_comm = mpi_comm_to_kid( which_kid(ngrid%id), grid%id )
        ioffset = nest_task_offsets(ngrid%id)
      ELSE IF ( intercomm_active( ngrid%id ) ) THEN  ! I am nest
        local_comm = mpi_comm_to_mom( ngrid%id )
        ioffset = nest_task_offsets(ngrid%id)
      END IF

      IF ( grid%active_this_task .OR. ngrid%active_this_task ) THEN
#ifndef STUBMPI
        CALL mpi_comm_rank(local_comm,myproc,ierr)
        CALL mpi_comm_size(local_comm,nproc,ierr)
#endif
!call tracebackqq()
        CALL rsl_lite_merge_msgs( myproc, nest_pes_x(grid%id)*nest_pes_y(grid%id),         &
                                          nest_pes_x(ngrid%id)*nest_pes_y(ngrid%id),       &
                                          ioffset, local_comm )
      END IF

IF ( grid%active_this_task ) THEN
    CALL push_communicators_for_domain( grid%id )


#define NEST_INFLUENCE(A,B) A = B
#include "nest_feedbackup_unpack.inc"

      ! smooth coarse grid
      CALL get_ijk_from_grid (  ngrid,                           &
                                nids, nide, njds, njde, nkds, nkde,    &
                                nims, nime, njms, njme, nkms, nkme,    &
                                nips, nipe, njps, njpe, nkps, nkpe    )
      CALL get_ijk_from_grid (  grid ,              &
                                ids, ide, jds, jde, kds, kde,    &
                                ims, ime, jms, jme, kms, kme,    &
                                ips, ipe, jps, jpe, kps, kpe    )

#include "HALO_INTERP_UP.inc"

      CALL get_ijk_from_grid (  grid ,                   &
                                cids, cide, cjds, cjde, ckds, ckde,    &
                                cims, cime, cjms, cjme, ckms, ckme,    &
                                cips, cipe, cjps, cjpe, ckps, ckpe    )

#include "nest_feedbackup_smooth.inc"

    CALL pop_communicators_for_domain
END IF

      RETURN
   END SUBROUTINE feedback_domain_em_part2
#endif


!------------------------------------------------------------------

   SUBROUTINE wrf_gatherv_real (Field, field_ofst,            &
                                my_count ,                    &    ! sendcount
                                globbuf, glob_ofst ,          &    ! recvbuf
                                counts                      , &    ! recvcounts
                                displs                      , &    ! displs
                                root                        , &    ! root
                                communicator                , &    ! communicator
                                ierr )
   USE module_dm, ONLY : getrealmpitype
   IMPLICIT NONE
   INTEGER field_ofst, glob_ofst
   INTEGER my_count, communicator, root, ierr
   INTEGER , DIMENSION(*) :: counts, displs
   REAL, DIMENSION(*) :: Field, globbuf
#ifndef STUBMPI
   INCLUDE 'mpif.h'

           CALL mpi_gatherv( Field( field_ofst ),      &    ! sendbuf
                            my_count ,                       &    ! sendcount
                            getrealmpitype() ,               &    ! sendtype
                            globbuf( glob_ofst ) ,                 &    ! recvbuf
                            counts                         , &    ! recvcounts
                            displs                         , &    ! displs
                            getrealmpitype()               , &    ! recvtype
                            root                           , &    ! root
                            communicator                   , &    ! communicator
                            ierr )
#endif

   END SUBROUTINE wrf_gatherv_real

   SUBROUTINE wrf_gatherv_double (Field, field_ofst,            &
                                my_count ,                    &    ! sendcount
                                globbuf, glob_ofst ,          &    ! recvbuf
                                counts                      , &    ! recvcounts
                                displs                      , &    ! displs
                                root                        , &    ! root
                                communicator                , &    ! communicator
                                ierr )
!   USE module_dm
   IMPLICIT NONE
   INTEGER field_ofst, glob_ofst
   INTEGER my_count, communicator, root, ierr
   INTEGER , DIMENSION(*) :: counts, displs
! this next declaration is REAL, not DOUBLE PRECISION because it will be autopromoted
! to double precision by the compiler when WRF is compiled for 8 byte reals. Only reason
! for having this separate routine is so we pass the correct MPI type to mpi_scatterv
! if we were not indexing the globbuf and Field arrays it would not even matter
   REAL, DIMENSION(*) :: Field, globbuf
#ifndef STUBMPI
   INCLUDE 'mpif.h'

           CALL mpi_gatherv( Field( field_ofst ),      &    ! sendbuf
                            my_count ,                       &    ! sendcount
                            MPI_DOUBLE_PRECISION         ,               &    ! sendtype
                            globbuf( glob_ofst ) ,                 &    ! recvbuf
                            counts                         , &    ! recvcounts
                            displs                         , &    ! displs
                            MPI_DOUBLE_PRECISION                       , &    ! recvtype
                            root                           , &    ! root
                            communicator                   , &    ! communicator
                            ierr )
#endif

   END SUBROUTINE wrf_gatherv_double

   SUBROUTINE wrf_gatherv_integer (Field, field_ofst,            &
                                my_count ,                    &    ! sendcount
                                globbuf, glob_ofst ,          &    ! recvbuf
                                counts                      , &    ! recvcounts
                                displs                      , &    ! displs
                                root                        , &    ! root
                                communicator                , &    ! communicator
                                ierr )
   IMPLICIT NONE
   INTEGER field_ofst, glob_ofst
   INTEGER my_count, communicator, root, ierr
   INTEGER , DIMENSION(*) :: counts, displs
   INTEGER, DIMENSION(*) :: Field, globbuf
#ifndef STUBMPI
   INCLUDE 'mpif.h'

           CALL mpi_gatherv( Field( field_ofst ),      &    ! sendbuf
                            my_count ,                       &    ! sendcount
                            MPI_INTEGER         ,               &    ! sendtype
                            globbuf( glob_ofst ) ,                 &    ! recvbuf
                            counts                         , &    ! recvcounts
                            displs                         , &    ! displs
                            MPI_INTEGER                       , &    ! recvtype
                            root                           , &    ! root
                            communicator                   , &    ! communicator
                            ierr )
#endif

   END SUBROUTINE wrf_gatherv_integer

!new stuff 20070124
   SUBROUTINE wrf_scatterv_real (                             &
                                globbuf, glob_ofst ,          &    ! recvbuf
                                counts                      , &    ! recvcounts
                                Field, field_ofst,            &
                                my_count ,                    &    ! sendcount
                                displs                      , &    ! displs
                                root                        , &    ! root
                                communicator                , &    ! communicator
                                ierr )
   USE module_dm, ONLY : getrealmpitype
   IMPLICIT NONE
   INTEGER field_ofst, glob_ofst
   INTEGER my_count, communicator, root, ierr
   INTEGER , DIMENSION(*) :: counts, displs
   REAL, DIMENSION(*) :: Field, globbuf
#ifndef STUBMPI
   INCLUDE 'mpif.h'

           CALL mpi_scatterv(                                &
                            globbuf( glob_ofst ) ,           &    ! recvbuf
                            counts                         , &    ! recvcounts
                            displs                         , &    ! displs
                            getrealmpitype()               , &    ! recvtype
                            Field( field_ofst ),             &    ! sendbuf
                            my_count ,                       &    ! sendcount
                            getrealmpitype() ,               &    ! sendtype
                            root                           , &    ! root
                            communicator                   , &    ! communicator
                            ierr )
#endif

   END SUBROUTINE wrf_scatterv_real

   SUBROUTINE wrf_scatterv_double (                           &
                                globbuf, glob_ofst ,          &    ! recvbuf
                                counts                      , &    ! recvcounts
                                Field, field_ofst,            &
                                my_count ,                    &    ! sendcount
                                displs                      , &    ! displs
                                root                        , &    ! root
                                communicator                , &    ! communicator
                                ierr )
   IMPLICIT NONE
   INTEGER field_ofst, glob_ofst
   INTEGER my_count, communicator, root, ierr
   INTEGER , DIMENSION(*) :: counts, displs
   REAL, DIMENSION(*) :: Field, globbuf
#ifndef STUBMPI
   INCLUDE 'mpif.h'
! this next declaration is REAL, not DOUBLE PRECISION because it will be autopromoted
! to double precision by the compiler when WRF is compiled for 8 byte reals. Only reason
! for having this separate routine is so we pass the correct MPI type to mpi_scatterv
! if we were not indexing the globbuf and Field arrays it would not even matter

           CALL mpi_scatterv(                                &
                            globbuf( glob_ofst ) ,           &    ! recvbuf
                            counts                         , &    ! recvcounts
                            displs                         , &    ! displs
                            MPI_DOUBLE_PRECISION           , &    ! recvtype
                            Field( field_ofst ),             &    ! sendbuf
                            my_count ,                       &    ! sendcount
                            MPI_DOUBLE_PRECISION         ,   &    ! sendtype
                            root                           , &    ! root
                            communicator                   , &    ! communicator
                            ierr )
#endif

   END SUBROUTINE wrf_scatterv_double

   SUBROUTINE wrf_scatterv_integer (                          &
                                globbuf, glob_ofst ,          &    ! recvbuf
                                counts                      , &    ! recvcounts
                                Field, field_ofst,            &
                                my_count ,                    &    ! sendcount
                                displs                      , &    ! displs
                                root                        , &    ! root
                                communicator                , &    ! communicator
                                ierr )
   IMPLICIT NONE
   INTEGER field_ofst, glob_ofst
   INTEGER my_count, communicator, root, ierr
   INTEGER , DIMENSION(*) :: counts, displs
   INTEGER, DIMENSION(*) :: Field, globbuf
#ifndef STUBMPI
   INCLUDE 'mpif.h'

           CALL mpi_scatterv(                                &
                            globbuf( glob_ofst ) ,           &    ! recvbuf
                            counts                         , &    ! recvcounts
                            displs                         , &    ! displs
                            MPI_INTEGER                    , &    ! recvtype
                            Field( field_ofst ),             &    ! sendbuf
                            my_count ,                       &    ! sendcount
                            MPI_INTEGER         ,            &    ! sendtype
                            root                           , &    ! root
                            communicator                   , &    ! communicator
                            ierr )
#endif

   END SUBROUTINE wrf_scatterv_integer
! end new stuff 20070124

     SUBROUTINE wrf_dm_gatherv ( v, elemsize , km_s, km_e, wordsz )
      IMPLICIT NONE
      INTEGER  elemsize, km_s, km_e, wordsz
      REAL v(*)
      IF ( wordsz .EQ. DWORDSIZE ) THEN
         CALL wrf_dm_gatherv_double(v, elemsize , km_s, km_e)
      ELSE
         CALL wrf_dm_gatherv_single(v, elemsize , km_s, km_e)
      END IF
     END SUBROUTINE wrf_dm_gatherv

     SUBROUTINE wrf_dm_gatherv_double ( v, elemsize , km_s, km_e )
      IMPLICIT NONE
      INTEGER  elemsize, km_s, km_e
      REAL*8 v(0:*)
#ifndef STUBMPI
# ifndef USE_MPI_IN_PLACE
      REAL*8 v_local((km_e-km_s+1)*elemsize)
# endif
      INTEGER, DIMENSION(:), ALLOCATABLE :: recvcounts, displs
      INTEGER send_type, myproc, nproc, local_comm, ierr, i
   INCLUDE 'mpif.h'
      send_type = MPI_DOUBLE_PRECISION
      CALL wrf_get_dm_communicator ( local_comm )
      CALL wrf_get_nproc( nproc )
      CALL wrf_get_myproc( myproc )
      ALLOCATE( recvcounts(nproc), displs(nproc) )
      i = (km_e-km_s+1)*elemsize
      CALL mpi_allgather( i,1,MPI_INTEGER,recvcounts,1,MPI_INTEGER,local_comm,ierr) ;
      i = (km_s)*elemsize
      CALL mpi_allgather( i,1,MPI_INTEGER,displs,1,MPI_INTEGER,local_comm,ierr) ;
# ifdef USE_MPI_IN_PLACE
      CALL mpi_allgatherv( MPI_IN_PLACE,                                  &
# else
      DO i = 1,elemsize*(km_e-km_s+1)
        v_local(i) = v(i+elemsize*km_s-1)
      END DO
      CALL mpi_allgatherv( v_local,                                       &
# endif
                           (km_e-km_s+1)*elemsize,                        &
                           send_type,                                     &
                           v,                                             &
                           recvcounts,                                    &
                           displs,                                        &
                           send_type,                                     &
                           local_comm,                                    &
                           ierr )
      DEALLOCATE(recvcounts)
      DEALLOCATE(displs)
#endif
      return
     END SUBROUTINE wrf_dm_gatherv_double

     SUBROUTINE wrf_dm_gatherv_single ( v, elemsize , km_s, km_e )
      IMPLICIT NONE
      INTEGER  elemsize, km_s, km_e
      REAL*4 v(0:*)
#ifndef STUBMPI
# ifndef USE_MPI_IN_PLACE
      REAL*4 v_local((km_e-km_s+1)*elemsize)
# endif
      INTEGER, DIMENSION(:), ALLOCATABLE :: recvcounts, displs
      INTEGER send_type, myproc, nproc, local_comm, ierr, i
   INCLUDE 'mpif.h'
      send_type = MPI_REAL
      CALL wrf_get_dm_communicator ( local_comm )
      CALL wrf_get_nproc( nproc )
      CALL wrf_get_myproc( myproc )
      ALLOCATE( recvcounts(nproc), displs(nproc) )
      i = (km_e-km_s+1)*elemsize
      CALL mpi_allgather( i,1,MPI_INTEGER,recvcounts,1,MPI_INTEGER,local_comm,ierr) ;
      i = (km_s)*elemsize
      CALL mpi_allgather( i,1,MPI_INTEGER,displs,1,MPI_INTEGER,local_comm,ierr) ;
# ifdef USE_MPI_IN_PLACE
      CALL mpi_allgatherv( MPI_IN_PLACE,                                  &
# else
      DO i = 1,elemsize*(km_e-km_s+1)
        v_local(i) = v(i+elemsize*km_s-1)
      END DO
      CALL mpi_allgatherv( v_local,                                       &
# endif
                           (km_e-km_s+1)*elemsize,                        &
                           send_type,                                     &
                           v,                                             &
                           recvcounts,                                    &
                           displs,                                        &
                           send_type,                                     &
                           local_comm,                                    &
                           ierr )
      DEALLOCATE(recvcounts)
      DEALLOCATE(displs)
#endif
      return
     END SUBROUTINE wrf_dm_gatherv_single

      SUBROUTINE wrf_dm_decomp1d( nt, km_s, km_e )
       IMPLICIT NONE
       INTEGER, INTENT(IN)  :: nt
       INTEGER, INTENT(OUT) :: km_s, km_e
     ! local
       INTEGER nn, nnp,  na, nb
       INTEGER myproc, nproc

       CALL wrf_get_myproc(myproc)
       CALL wrf_get_nproc(nproc)
       nn = nt / nproc           ! min number done by this task
       nnp = nn
       if ( myproc .lt. mod( nt, nproc ) )   nnp = nnp + 1 ! distribute remainder

       na = min( myproc, mod(nt,nproc) ) ! Number of blocks with remainder that precede this one
       nb = max( 0, myproc - na )        ! number of blocks without a remainder that precede this one
       km_s = na * ( nn+1) + nb * nn     ! starting iteration for this task
       km_e = km_s + nnp - 1             ! ending iteration for this task
      END SUBROUTINE wrf_dm_decomp1d


SUBROUTINE wrf_dm_define_comms ( grid )
   USE module_domain, ONLY : domain
   IMPLICIT NONE
   TYPE(domain) , INTENT (INOUT) :: grid
   RETURN
END SUBROUTINE wrf_dm_define_comms

SUBROUTINE tfp_message( fname, lno )
   CHARACTER*(*) fname
   INTEGER lno
   CHARACTER*1024 mess
#ifndef STUBMPI
   WRITE(mess,*)'tfp_message: ',trim(fname),lno
   CALL wrf_message(mess)
# ifdef ALLOW_OVERDECOMP
     CALL task_for_point_message  ! defined in RSL_LITE/task_for_point.c
# else
     CALL wrf_error_fatal(mess)
# endif
#endif
END SUBROUTINE tfp_message

   SUBROUTINE set_dm_debug
    USE module_dm, ONLY : dm_debug_flag
    IMPLICIT NONE
    dm_debug_flag = .TRUE.
   END SUBROUTINE set_dm_debug
   SUBROUTINE reset_dm_debug
    USE module_dm, ONLY : dm_debug_flag
    IMPLICIT NONE
    dm_debug_flag = .FALSE.
   END SUBROUTINE reset_dm_debug
   SUBROUTINE get_dm_debug ( arg )
    USE module_dm, ONLY : dm_debug_flag
    IMPLICIT NONE
    LOGICAL arg
    arg = dm_debug_flag
   END SUBROUTINE get_dm_debug

