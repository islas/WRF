
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
