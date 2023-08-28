submodule (module_dm) CODE_0
implicit none
contains
module SUBROUTINE patch_domain_rsl_lite( id  , parent, parent_id, &
                                sd1 , ed1 , sp1 , ep1 , sm1 , em1 ,        &
                                sd2 , ed2 , sp2 , ep2 , sm2 , em2 ,        &
                                sd3 , ed3 , sp3 , ep3 , sm3 , em3 ,        &
                                      sp1x , ep1x , sm1x , em1x , &
                                      sp2x , ep2x , sm2x , em2x , &
                                      sp3x , ep3x , sm3x , em3x , &
                                      sp1y , ep1y , sm1y , em1y , &
                                      sp2y , ep2y , sm2y , em2y , &
                                      sp3y , ep3y , sm3y , em3y , &
                                bdx , bdy )

#if ( ( defined(SGIALTIX) || defined(FUJITSU_FX10) || defined(KEEP_INT_AROUND) ) && (! defined(MOVE_NESTS) ) )
      USE module_domain, ONLY : domain, head_grid, find_grid_by_id, alloc_space_field
#else
      USE module_domain, ONLY : domain, head_grid, find_grid_by_id
#endif

      IMPLICIT NONE
      INTEGER, INTENT(IN)   :: sd1 , ed1 , sd2 , ed2 , sd3 , ed3 , bdx , bdy
      INTEGER, INTENT(OUT)  :: sp1 , ep1 , sp2 , ep2 , sp3 , ep3 , &
                               sm1 , em1 , sm2 , em2 , sm3 , em3
      INTEGER, INTENT(OUT)  :: sp1x , ep1x , sp2x , ep2x , sp3x , ep3x , &
                               sm1x , em1x , sm2x , em2x , sm3x , em3x
      INTEGER, INTENT(OUT)  :: sp1y , ep1y , sp2y , ep2y , sp3y , ep3y , &
                               sm1y , em1y , sm2y , em2y , sm3y , em3y
      INTEGER, INTENT(IN)   :: id, parent_id
      TYPE(domain),POINTER  :: parent

! Local variables
      INTEGER               :: ids, ide, jds, jde, kds, kde
      INTEGER               :: ims, ime, jms, jme, kms, kme
      INTEGER               :: ips, ipe, jps, jpe, kps, kpe
      INTEGER               :: imsx, imex, jmsx, jmex, kmsx, kmex
      INTEGER               :: ipsx, ipex, jpsx, jpex, kpsx, kpex
      INTEGER               :: imsy, imey, jmsy, jmey, kmsy, kmey
      INTEGER               :: ipsy, ipey, jpsy, jpey, kpsy, kpey

      INTEGER               :: c_sd1 , c_ed1 , c_sd2 , c_ed2 , c_sd3 , c_ed3
      INTEGER               :: c_sp1 , c_ep1 , c_sp2 , c_ep2 , c_sp3 , c_ep3 , &
                               c_sm1 , c_em1 , c_sm2 , c_em2 , c_sm3 , c_em3
      INTEGER               :: c_sp1x , c_ep1x , c_sp2x , c_ep2x , c_sp3x , c_ep3x , &
                               c_sm1x , c_em1x , c_sm2x , c_em2x , c_sm3x , c_em3x
      INTEGER               :: c_sp1y , c_ep1y , c_sp2y , c_ep2y , c_sp3y , c_ep3y , &
                               c_sm1y , c_em1y , c_sm2y , c_em2y , c_sm3y , c_em3y

      INTEGER               :: c_ids, c_ide, c_jds, c_jde, c_kds, c_kde
      INTEGER               :: c_ims, c_ime, c_jms, c_jme, c_kms, c_kme
      INTEGER               :: c_ips, c_ipe, c_jps, c_jpe, c_kps, c_kpe

      INTEGER               :: idim , jdim , kdim , rem , a, b
      INTEGER               :: i, j, ni, nj, Px, Py, P

      INTEGER               :: parent_grid_ratio, i_parent_start, j_parent_start
      INTEGER               :: shw
      INTEGER               :: idim_cd, jdim_cd, ierr
      INTEGER               :: max_dom

#if (DA_CORE == 1)
      INTEGER               :: e_we, e_sn
#endif

      TYPE(domain), POINTER :: intermediate_grid
      TYPE(domain), POINTER  :: nest_grid
      CHARACTER*256   :: mess

      INTEGER parent_max_halo_width
      INTEGER thisdomain_max_halo_width
      INTEGER lats_to_mic

     lats_to_mic=0
     CALL nl_get_lats_to_mic( 1, lats_to_mic )
      IF ( lats_to_mic .GT. 0 ) THEN
        minx = -99  ! code to task_for_point to do split decomposition over MIC and host
        miny = lats_to_mic  ! number of latitudes that should be assigned to MIC
      ELSE
        minx = 1   ! normal
        miny = 1   ! normal
      END IF



      SELECT CASE ( model_data_order )
         ! need to finish other cases
         CASE ( DATA_ORDER_ZXY )
            ids = sd2 ; ide = ed2
            jds = sd3 ; jde = ed3
            kds = sd1 ; kde = ed1
         CASE ( DATA_ORDER_XYZ )
            ids = sd1 ; ide = ed1
            jds = sd2 ; jde = ed2
            kds = sd3 ; kde = ed3
         CASE ( DATA_ORDER_XZY )
            ids = sd1 ; ide = ed1
            jds = sd3 ; jde = ed3
            kds = sd2 ; kde = ed2
         CASE ( DATA_ORDER_YXZ)
            ids = sd2 ; ide = ed2
            jds = sd1 ; jde = ed1
            kds = sd3 ; kde = ed3
      END SELECT

      CALL nl_get_max_dom( 1 , max_dom )

      CALL get_dm_max_halo_width( id , thisdomain_max_halo_width )
      IF ( id .GT. 1 ) THEN
        CALL get_dm_max_halo_width( parent%id , parent_max_halo_width )
      END IF

      CALL compute_memory_dims_rsl_lite ( id, thisdomain_max_halo_width, 0 , bdx, bdy,   &
                   ids,  ide,  jds,  jde,  kds,  kde, &
                   ims,  ime,  jms,  jme,  kms,  kme, &
                   imsx, imex, jmsx, jmex, kmsx, kmex, &
                   imsy, imey, jmsy, jmey, kmsy, kmey, &
                   ips,  ipe,  jps,  jpe,  kps,  kpe, &
                   ipsx, ipex, jpsx, jpex, kpsx, kpex, &
                   ipsy, ipey, jpsy, jpey, kpsy, kpey )

     ! ensure that the every parent domain point has a full set of nested points under it
     ! even at the borders. Do this by making sure the number of nest points is a multiple of
     ! the nesting ratio. Note that this is important mostly to the intermediate domain, which
     ! is the subject of the scatter gather comms with the parent

      IF ( id .GT. 1 ) THEN
         CALL nl_get_parent_grid_ratio( id, parent_grid_ratio )
         if ( mod(ime,parent_grid_ratio) .NE. 0 ) ime = ime + parent_grid_ratio - mod(ime,parent_grid_ratio)
         if ( mod(jme,parent_grid_ratio) .NE. 0 ) jme = jme + parent_grid_ratio - mod(jme,parent_grid_ratio)
      END IF

      SELECT CASE ( model_data_order )
         CASE ( DATA_ORDER_ZXY )
            sp2 = ips ; ep2 = ipe ; sm2 = ims ; em2 = ime
            sp3 = jps ; ep3 = jpe ; sm3 = jms ; em3 = jme
            sp1 = kps ; ep1 = kpe ; sm1 = kms ; em1 = kme
            sp2x = ipsx ; ep2x = ipex ; sm2x = imsx ; em2x = imex
            sp3x = jpsx ; ep3x = jpex ; sm3x = jmsx ; em3x = jmex
            sp1x = kpsx ; ep1x = kpex ; sm1x = kmsx ; em1x = kmex
            sp2y = ipsy ; ep2y = ipey ; sm2y = imsy ; em2y = imey
            sp3y = jpsy ; ep3y = jpey ; sm3y = jmsy ; em3y = jmey
            sp1y = kpsy ; ep1y = kpey ; sm1y = kmsy ; em1y = kmey
         CASE ( DATA_ORDER_ZYX )
            sp3 = ips ; ep3 = ipe ; sm3 = ims ; em3 = ime
            sp2 = jps ; ep2 = jpe ; sm2 = jms ; em2 = jme
            sp1 = kps ; ep1 = kpe ; sm1 = kms ; em1 = kme
            sp3x = ipsx ; ep3x = ipex ; sm3x = imsx ; em3x = imex
            sp2x = jpsx ; ep2x = jpex ; sm2x = jmsx ; em2x = jmex
            sp1x = kpsx ; ep1x = kpex ; sm1x = kmsx ; em1x = kmex
            sp3y = ipsy ; ep3y = ipey ; sm3y = imsy ; em3y = imey
            sp2y = jpsy ; ep2y = jpey ; sm2y = jmsy ; em2y = jmey
            sp1y = kpsy ; ep1y = kpey ; sm1y = kmsy ; em1y = kmey
         CASE ( DATA_ORDER_XYZ )
            sp1 = ips ; ep1 = ipe ; sm1 = ims ; em1 = ime
            sp2 = jps ; ep2 = jpe ; sm2 = jms ; em2 = jme
            sp3 = kps ; ep3 = kpe ; sm3 = kms ; em3 = kme
            sp1x = ipsx ; ep1x = ipex ; sm1x = imsx ; em1x = imex
            sp2x = jpsx ; ep2x = jpex ; sm2x = jmsx ; em2x = jmex
            sp3x = kpsx ; ep3x = kpex ; sm3x = kmsx ; em3x = kmex
            sp1y = ipsy ; ep1y = ipey ; sm1y = imsy ; em1y = imey
            sp2y = jpsy ; ep2y = jpey ; sm2y = jmsy ; em2y = jmey
            sp3y = kpsy ; ep3y = kpey ; sm3y = kmsy ; em3y = kmey
         CASE ( DATA_ORDER_YXZ)
            sp2 = ips ; ep2 = ipe ; sm2 = ims ; em2 = ime
            sp1 = jps ; ep1 = jpe ; sm1 = jms ; em1 = jme
            sp3 = kps ; ep3 = kpe ; sm3 = kms ; em3 = kme
            sp2x = ipsx ; ep2x = ipex ; sm2x = imsx ; em2x = imex
            sp1x = jpsx ; ep1x = jpex ; sm1x = jmsx ; em1x = jmex
            sp3x = kpsx ; ep3x = kpex ; sm3x = kmsx ; em3x = kmex
            sp2y = ipsy ; ep2y = ipey ; sm2y = imsy ; em2y = imey
            sp1y = jpsy ; ep1y = jpey ; sm1y = jmsy ; em1y = jmey
            sp3y = kpsy ; ep3y = kpey ; sm3y = kmsy ; em3y = kmey
         CASE ( DATA_ORDER_XZY )
            sp1 = ips ; ep1 = ipe ; sm1 = ims ; em1 = ime
            sp3 = jps ; ep3 = jpe ; sm3 = jms ; em3 = jme
            sp2 = kps ; ep2 = kpe ; sm2 = kms ; em2 = kme
            sp1x = ipsx ; ep1x = ipex ; sm1x = imsx ; em1x = imex
            sp3x = jpsx ; ep3x = jpex ; sm3x = jmsx ; em3x = jmex
            sp2x = kpsx ; ep2x = kpex ; sm2x = kmsx ; em2x = kmex
            sp1y = ipsy ; ep1y = ipey ; sm1y = imsy ; em1y = imey
            sp3y = jpsy ; ep3y = jpey ; sm3y = jmsy ; em3y = jmey
            sp2y = kpsy ; ep2y = kpey ; sm2y = kmsy ; em2y = kmey
         CASE ( DATA_ORDER_YZX )
            sp3 = ips ; ep3 = ipe ; sm3 = ims ; em3 = ime
            sp1 = jps ; ep1 = jpe ; sm1 = jms ; em1 = jme
            sp2 = kps ; ep2 = kpe ; sm2 = kms ; em2 = kme
            sp3x = ipsx ; ep3x = ipex ; sm3x = imsx ; em3x = imex
            sp1x = jpsx ; ep1x = jpex ; sm1x = jmsx ; em1x = jmex
            sp2x = kpsx ; ep2x = kpex ; sm2x = kmsx ; em2x = kmex
            sp3y = ipsy ; ep3y = ipey ; sm3y = imsy ; em3y = imey
            sp1y = jpsy ; ep1y = jpey ; sm1y = jmsy ; em1y = jmey
            sp2y = kpsy ; ep2y = kpey ; sm2y = kmsy ; em2y = kmey
      END SELECT

      IF ( id.EQ.1 ) THEN
         WRITE(wrf_err_message,*)'*************************************'
         CALL wrf_message( TRIM(wrf_err_message) )
         WRITE(wrf_err_message,*)'Parent domain'
         CALL wrf_message( TRIM(wrf_err_message) )
         WRITE(wrf_err_message,*)'ids,ide,jds,jde ',ids,ide,jds,jde
         CALL wrf_message( TRIM(wrf_err_message) )
         WRITE(wrf_err_message,*)'ims,ime,jms,jme ',ims,ime,jms,jme
         CALL wrf_message( TRIM(wrf_err_message) )
         WRITE(wrf_err_message,*)'ips,ipe,jps,jpe ',ips,ipe,jps,jpe
         CALL wrf_message( TRIM(wrf_err_message) )
         WRITE(wrf_err_message,*)'*************************************'
         CALL wrf_message( TRIM(wrf_err_message) )
      END IF

      IF ( id .GT. 1 ) THEN

         CALL nl_get_shw( id, shw )
         CALL nl_get_i_parent_start( id , i_parent_start )
         CALL nl_get_j_parent_start( id , j_parent_start )
         CALL nl_get_parent_grid_ratio( id, parent_grid_ratio )

         SELECT CASE ( model_data_order )
            CASE ( DATA_ORDER_ZXY )
               idim = ed2-sd2+1
               jdim = ed3-sd3+1
               kdim = ed1-sd1+1
               c_kds = sd1                ; c_kde = ed1
            CASE ( DATA_ORDER_ZYX )
               idim = ed3-sd3+1
               jdim = ed2-sd2+1
               kdim = ed1-sd1+1
               c_kds = sd1                ; c_kde = ed1
            CASE ( DATA_ORDER_XYZ )
               idim = ed1-sd1+1
               jdim = ed2-sd2+1
               kdim = ed3-sd3+1
               c_kds = sd3                ; c_kde = ed3
            CASE ( DATA_ORDER_YXZ)
               idim = ed2-sd2+1
               jdim = ed1-sd1+1
               kdim = ed3-sd3+1
               c_kds = sd3                ; c_kde = ed3
            CASE ( DATA_ORDER_XZY )
               idim = ed1-sd1+1
               jdim = ed3-sd3+1
               kdim = ed2-sd2+1
               c_kds = sd2                ; c_kde = ed2
            CASE ( DATA_ORDER_YZX )
               idim = ed3-sd3+1
               jdim = ed1-sd1+1
               kdim = ed2-sd2+1
               c_kds = sd2                ; c_kde = ed2
         END SELECT

         idim_cd = idim / parent_grid_ratio + 1 + 2*shw + 1
         jdim_cd = jdim / parent_grid_ratio + 1 + 2*shw + 1

         c_ids = i_parent_start-shw ; c_ide = c_ids + idim_cd - 1
         c_jds = j_parent_start-shw ; c_jde = c_jds + jdim_cd - 1

#if (DA_CORE == 1)
          call nl_get_e_we( id -1, e_we )
          call nl_get_e_sn( id -1, e_sn )

         if ( c_ids .le. 0   ) c_ids = 1
         if ( c_ide .gt. e_we) c_ide = e_we
         if ( c_jds .le. 0   ) c_jds = 1
         if ( c_jde .gt. e_sn) c_jde = e_sn
#endif
         ! we want the intermediate domain to be decomposed the
         ! the same as the underlying nest. So try this:

         c_ips = -1
         nj = ( c_jds - j_parent_start ) * parent_grid_ratio + 1 + 1 ;
         ierr = 0
         DO i = c_ids, c_ide
            ni = ( i - i_parent_start ) * parent_grid_ratio + 1 + 1 ;
!jm            CALL task_for_point ( ni, nj, ids, ide, jds, jde, ntasks_x, ntasks_y, Px, Py, &
            CALL task_for_point ( ni, nj, ids, ide, jds, jde, nest_pes_x(id), nest_pes_y(id),Px,Py, &
                                  minx, miny,  ierr )
            IF ( ierr .NE. 0 ) CALL wrf_error_fatal('error code returned by task_for_point in module_dm.F (a)')
            IF ( Px .EQ. mytask_x ) THEN
               c_ipe = i
               IF ( c_ips .EQ. -1 ) c_ips = i
            END IF
         END DO
         IF ( ierr .NE. 0 ) THEN
            CALL tfp_message("<stdin>",__LINE__)
         END IF
         IF (c_ips .EQ. -1 ) THEN
            c_ipe = -1
            c_ips = 0
         END IF

         c_jps = -1
         ni = ( c_ids - i_parent_start ) * parent_grid_ratio + 1 + 1 ;
         ierr = 0
         DO j = c_jds, c_jde
            nj = ( j - j_parent_start ) * parent_grid_ratio + 1 + 1 ;
!            CALL task_for_point ( ni, nj, ids, ide, jds, jde, ntasks_x, ntasks_y, Px, Py, &
            CALL task_for_point ( ni, nj, ids, ide, jds, jde, nest_pes_x(id), nest_pes_y(id), Px, Py, &
                                  minx, miny, ierr )
            IF ( ierr .NE. 0 ) CALL wrf_error_fatal('error code returned by task_for_point in module_dm.F (b)')


            IF ( Py .EQ. mytask_y ) THEN
               c_jpe = j
               IF ( c_jps .EQ. -1 ) c_jps = j
            END IF
         END DO
         IF ( ierr .NE. 0 ) THEN
            CALL tfp_message("<stdin>",__LINE__)
         END IF
         IF (c_jps .EQ. -1 ) THEN
            c_jpe = -1
            c_jps = 0
         END IF

#if (DA_CORE == 1)
         IF (c_ipe .EQ. -1 .or. c_jpe .EQ. -1) THEN
            c_ipe = -1
            c_ips = 0
            c_jpe = -1
            c_jps = 0
         END IF


          c_kpsx = -1
          nj = ( c_jds - j_parent_start ) * parent_grid_ratio + 1 + 1 ;
          ierr = 0
          DO k = c_kds, c_kde
!             CALL task_for_point ( k, nj, kds, kde, jds, jde, ntasks_x, ntasks_y, Px, Py, &
             CALL task_for_point ( k, nj, kds, kde, jds, jde, nest_pes_x(id), nest_pes_y(id), Px, Py, &
                                   1, 1, ierr )
             IF ( Px .EQ. mytask_x ) THEN
                c_kpex = k
                IF ( c_kpsx .EQ. -1 ) c_kpsx = k
             END IF
          END DO
          IF ( ierr .NE. 0 ) THEN
             CALL tfp_message("<stdin>",__LINE__)
          END IF
          IF (c_kpsx .EQ. -1 ) THEN
             c_kpex = -1
             c_kpsx = 0
          END IF

          c_jpsx = -1
          k = c_kds ;
          ierr = 0
          DO j = c_jds, c_jde
             nj = ( j - j_parent_start ) * parent_grid_ratio + 1 + 1 ;
!             CALL task_for_point ( k, nj, kds, kde, jds, jde, ntasks_x, ntasks_y, Px, Py, &
             CALL task_for_point ( k, nj, kds, kde, jds, jde, nest_pes_x(id), nest_pes_y(id), Px, Py, &
                                   1, 1, ierr )
             IF ( Py .EQ. mytask_y ) THEN
                c_jpex = j
                IF ( c_jpsx .EQ. -1 ) c_jpsx = j
             END IF
          END DO
          IF ( ierr .NE. 0 ) THEN
             CALL tfp_message("<stdin>",__LINE__)
          END IF
          IF (c_jpsx .EQ. -1 ) THEN
             c_jpex = -1
             c_jpsx = 0
          END IF

          IF (c_ipex .EQ. -1 .or. c_jpex .EQ. -1) THEN
             c_ipex = -1
             c_ipsx = 0
             c_jpex = -1
             c_jpsx = 0
          END IF

          c_kpsy = c_kpsx   ! same as above
          c_kpey = c_kpex   ! same as above

          c_ipsy = -1
          k = c_kds ;
          ierr = 0
          DO i = c_ids, c_ide
             ni = ( i - i_parent_start ) * parent_grid_ratio + 1 + 1 ;
!             CALL task_for_point ( ni, k, ids, ide, kds, kde, ntasks_y, ntasks_x, Py, Px, &
             CALL task_for_point ( ni, k, ids, ide, kds, kde, nest_pes_y(id), nest_pes_x(id), Py, Px, &
                                   1, 1, ierr ) ! x and y for proc mesh reversed
             IF ( Py .EQ. mytask_y ) THEN
                c_ipey = i
                IF ( c_ipsy .EQ. -1 ) c_ipsy = i
             END IF
          END DO
          IF ( ierr .NE. 0 ) THEN
             CALL tfp_message("<stdin>",__LINE__)
          END IF
          IF (c_ipsy .EQ. -1 ) THEN
             c_ipey = -1
             c_ipsy = 0
          END IF
#endif


         IF ( c_ips <= c_ipe ) THEN
! extend the patch dimensions out shw along edges of domain
           IF ( mytask_x .EQ. 0 ) THEN
             c_ips = c_ips - shw
#if (DA_CORE == 1)
             c_ipsy = c_ipsy - shw  
#endif
           END IF
!           IF ( mytask_x .EQ. ntasks_x-1 ) THEN
           IF ( mytask_x .EQ. nest_pes_x(id)-1 ) THEN
             c_ipe = c_ipe + shw
#if (DA_CORE == 1)
             c_ipey = c_ipey + shw  
#endif
           END IF
           c_ims = max( c_ips - max(shw,thisdomain_max_halo_width), c_ids - bdx ) - 1
           c_ime = min( c_ipe + max(shw,thisdomain_max_halo_width), c_ide + bdx ) + 1
         ELSE
           c_ims = 0
           c_ime = 0
         END IF


! handle j dims
         IF ( c_jps <= c_jpe ) THEN
! extend the patch dimensions out shw along edges of domain
           IF ( mytask_y .EQ. 0 ) THEN
              c_jps = c_jps - shw
#if (DA_CORE == 1)
              c_jpsx = c_jpsx - shw  
#endif
           END IF
!           IF ( mytask_y .EQ. ntasks_y-1 ) THEN
           IF ( mytask_y .EQ. nest_pes_y(id)-1 ) THEN
              c_jpe = c_jpe + shw
#if (DA_CORE == 1)
              c_jpex = c_jpex + shw  
#endif
           END IF
           c_jms = max( c_jps - max(shw,thisdomain_max_halo_width), c_jds - bdx ) - 1
           c_jme = min( c_jpe + max(shw,thisdomain_max_halo_width), c_jde + bdx ) + 1
! handle k dims
         ELSE
           c_jms = 0
           c_jme = 0
         END IF
         c_kps = 1
         c_kpe = c_kde
         c_kms = 1
         c_kme = c_kde

! Default initializations
         c_sm1x = 1 ; c_em1x = 1 ; c_sm2x = 1 ; c_em2x = 1 ; c_sm3x = 1 ; c_em3x = 1
         c_sm1y = 1 ; c_em1y = 1 ; c_sm2y = 1 ; c_em2y = 1 ; c_sm3y = 1 ; c_em3y = 1

#if (DA_CORE == 1)
         c_kmsx = c_kpsx
         c_kmex = c_kpex
         c_kmsy = c_kpsy
         c_kmey = c_kpey

         IF ( c_kpsx .EQ. 0 .AND. c_kpex .EQ. -1 ) THEN  
            c_kmsx = 0
            c_kmex = 0
         END IF
         IF ( c_kpsy .EQ. 0 .AND. c_kpey .EQ. -1 ) THEN
            c_kmsy = 0
            c_kmey = 0
         END IF
         c_imsx = c_ids
         c_imex = c_ide
         c_ipsx = c_imsx
         c_ipex = c_imex

         IF ( c_ipsy .EQ. 0 .AND. c_ipey .EQ. -1 ) THEN
            c_imsy = 0
            c_imey = 0
         ELSE
            c_imsy = c_ipsy
            c_imey = c_ipey
         END IF

         c_jmsx = c_jpsx
         c_jmex = c_jpex
         c_jmsy = c_jds
         c_jmey = c_jde

         IF ( c_jpsx .EQ. 0 .AND. c_jpex .EQ. -1 ) THEN
            c_jmsx = 0
            c_jmex = 0
         ELSE
            c_jpsy = c_jmsy
            c_jpey = c_jmey
         END IF

         c_sm1x = c_imsx
         c_em1x = c_imex
         c_sm2x = c_jmsx
         c_em2x = c_jmex
         c_sm3x = c_kmsx
         c_em3x = c_kmex

         c_sm1y = c_imsy
         c_em1y = c_imey
         c_sm2y = c_jmsy
         c_em2y = c_jmey
         c_sm3y = c_kmsy
         c_em3y = c_kmey

         c_sp1x = c_ipsx
         c_ep1x = c_ipex
         c_sp2x = c_jpsx
         c_ep2x = c_jpex
         c_sp3x = c_kpsx
         c_ep3x = c_kpex

         c_sp1y = c_ipsy
         c_ep1y = c_ipey
         c_sp2y = c_jpsy
         c_ep2y = c_jpey
         c_sp3y = c_kpsy
         c_ep3y = c_kpey
#endif

         WRITE(wrf_err_message,*)'*************************************'
         CALL wrf_message( TRIM(wrf_err_message) )
         WRITE(wrf_err_message,*)'Nesting domain'
         CALL wrf_message( TRIM(wrf_err_message) )
         WRITE(wrf_err_message,*)'ids,ide,jds,jde ',ids,ide,jds,jde
         CALL wrf_message( TRIM(wrf_err_message) )
         WRITE(wrf_err_message,*)'ims,ime,jms,jme ',ims,ime,jms,jme
         CALL wrf_message( TRIM(wrf_err_message) )
         WRITE(wrf_err_message,*)'ips,ipe,jps,jpe ',ips,ipe,jps,jpe
         CALL wrf_message( TRIM(wrf_err_message) )
         WRITE(wrf_err_message,*)'INTERMEDIATE domain'
         CALL wrf_message( TRIM(wrf_err_message) )
         WRITE(wrf_err_message,*)'ids,ide,jds,jde ',c_ids,c_ide,c_jds,c_jde
         CALL wrf_message( TRIM(wrf_err_message) )
         WRITE(wrf_err_message,*)'ims,ime,jms,jme ',c_ims,c_ime,c_jms,c_jme
         CALL wrf_message( TRIM(wrf_err_message) )
         WRITE(wrf_err_message,*)'ips,ipe,jps,jpe ',c_ips,c_ipe,c_jps,c_jpe
         CALL wrf_message( TRIM(wrf_err_message) )
         WRITE(wrf_err_message,*)'*************************************'
         CALL wrf_message( TRIM(wrf_err_message) )

         SELECT CASE ( model_data_order )
            CASE ( DATA_ORDER_ZXY )
               c_sd2 = c_ids ; c_ed2 = c_ide ; c_sp2 = c_ips ; c_ep2 = c_ipe ; c_sm2 = c_ims ; c_em2 = c_ime
               c_sd3 = c_jds ; c_ed3 = c_jde ; c_sp3 = c_jps ; c_ep3 = c_jpe ; c_sm3 = c_jms ; c_em3 = c_jme
               c_sd1 = c_kds ; c_ed1 = c_kde ; c_sp1 = c_kps ; c_ep1 = c_kpe ; c_sm1 = c_kms ; c_em1 = c_kme
            CASE ( DATA_ORDER_ZYX )
               c_sd3 = c_ids ; c_ed3 = c_ide ; c_sp3 = c_ips ; c_ep3 = c_ipe ; c_sm3 = c_ims ; c_em3 = c_ime
               c_sd2 = c_jds ; c_ed2 = c_jde ; c_sp2 = c_jps ; c_ep2 = c_jpe ; c_sm2 = c_jms ; c_em2 = c_jme
               c_sd1 = c_kds ; c_ed1 = c_kde ; c_sp1 = c_kps ; c_ep1 = c_kpe ; c_sm1 = c_kms ; c_em1 = c_kme
            CASE ( DATA_ORDER_XYZ )
               c_sd1 = c_ids ; c_ed1 = c_ide ; c_sp1 = c_ips ; c_ep1 = c_ipe ; c_sm1 = c_ims ; c_em1 = c_ime
               c_sd2 = c_jds ; c_ed2 = c_jde ; c_sp2 = c_jps ; c_ep2 = c_jpe ; c_sm2 = c_jms ; c_em2 = c_jme
               c_sd3 = c_kds ; c_ed3 = c_kde ; c_sp3 = c_kps ; c_ep3 = c_kpe ; c_sm3 = c_kms ; c_em3 = c_kme
            CASE ( DATA_ORDER_YXZ)
               c_sd2 = c_ids ; c_ed2 = c_ide ; c_sp2 = c_ips ; c_ep2 = c_ipe ; c_sm2 = c_ims ; c_em2 = c_ime
               c_sd1 = c_jds ; c_ed1 = c_jde ; c_sp1 = c_jps ; c_ep1 = c_jpe ; c_sm1 = c_jms ; c_em1 = c_jme
               c_sd3 = c_kds ; c_ed3 = c_kde ; c_sp3 = c_kps ; c_ep3 = c_kpe ; c_sm3 = c_kms ; c_em3 = c_kme
            CASE ( DATA_ORDER_XZY )
               c_sd1 = c_ids ; c_ed1 = c_ide ; c_sp1 = c_ips ; c_ep1 = c_ipe ; c_sm1 = c_ims ; c_em1 = c_ime
               c_sd3 = c_jds ; c_ed3 = c_jde ; c_sp3 = c_jps ; c_ep3 = c_jpe ; c_sm3 = c_jms ; c_em3 = c_jme
               c_sd2 = c_kds ; c_ed2 = c_kde ; c_sp2 = c_kps ; c_ep2 = c_kpe ; c_sm2 = c_kms ; c_em2 = c_kme
            CASE ( DATA_ORDER_YZX )
               c_sd3 = c_ids ; c_ed3 = c_ide ; c_sp3 = c_ips ; c_ep3 = c_ipe ; c_sm3 = c_ims ; c_em3 = c_ime
               c_sd1 = c_jds ; c_ed1 = c_jde ; c_sp1 = c_jps ; c_ep1 = c_jpe ; c_sm1 = c_jms ; c_em1 = c_jme
               c_sd2 = c_kds ; c_ed2 = c_kde ; c_sp2 = c_kps ; c_ep2 = c_kpe ; c_sm2 = c_kms ; c_em2 = c_kme
         END SELECT

         ALLOCATE ( intermediate_grid )
         ALLOCATE ( intermediate_grid%parents( max_parents ) )
         ALLOCATE ( intermediate_grid%nests( max_nests ) )
         intermediate_grid%allocated=.false.
         NULLIFY( intermediate_grid%sibling )
         DO i = 1, max_nests
            NULLIFY( intermediate_grid%nests(i)%ptr )
         END DO
         NULLIFY  (intermediate_grid%next)
         NULLIFY  (intermediate_grid%same_level)
         NULLIFY  (intermediate_grid%i_start)
         NULLIFY  (intermediate_grid%j_start)
         NULLIFY  (intermediate_grid%i_end)
         NULLIFY  (intermediate_grid%j_end)
         intermediate_grid%id = id   ! these must be the same. Other parts of code depend on it (see gen_comms.c)
         intermediate_grid%num_nests = 0
         intermediate_grid%num_siblings = 0
         intermediate_grid%num_parents = 1
         intermediate_grid%max_tiles   = 0
         intermediate_grid%num_tiles_spec   = 0
#if ( EM_CORE == 1 && DA_CORE != 1 )
         intermediate_grid%active_this_task = .true.
#endif
         CALL find_grid_by_id ( id, head_grid, nest_grid )

         nest_grid%intermediate_grid => intermediate_grid  ! nest grid now has a pointer to this baby
         intermediate_grid%parents(1)%ptr => nest_grid     ! the intermediate grid considers nest its parent
         intermediate_grid%num_parents = 1

         intermediate_grid%is_intermediate = .TRUE.
         SELECT CASE ( model_data_order )
            CASE ( DATA_ORDER_ZXY )
               intermediate_grid%nids = nest_grid%sd32 ; intermediate_grid%njds = nest_grid%sd33
               intermediate_grid%nide = nest_grid%ed32 ; intermediate_grid%njde = nest_grid%sd33
            CASE ( DATA_ORDER_ZYX )
               intermediate_grid%nids = nest_grid%sd33 ; intermediate_grid%njds = nest_grid%sd32
               intermediate_grid%nide = nest_grid%ed33 ; intermediate_grid%njde = nest_grid%sd32
            CASE ( DATA_ORDER_XYZ )
               intermediate_grid%nids = nest_grid%sd31 ; intermediate_grid%njds = nest_grid%sd32
               intermediate_grid%nide = nest_grid%ed31 ; intermediate_grid%njde = nest_grid%sd32
            CASE ( DATA_ORDER_YXZ)
               intermediate_grid%nids = nest_grid%sd32 ; intermediate_grid%njds = nest_grid%sd31
               intermediate_grid%nide = nest_grid%ed32 ; intermediate_grid%njde = nest_grid%sd31
            CASE ( DATA_ORDER_XZY )
               intermediate_grid%nids = nest_grid%sd31 ; intermediate_grid%njds = nest_grid%sd33
               intermediate_grid%nide = nest_grid%ed31 ; intermediate_grid%njde = nest_grid%sd33
            CASE ( DATA_ORDER_YZX )
               intermediate_grid%nids = nest_grid%sd33 ; intermediate_grid%njds = nest_grid%sd31
               intermediate_grid%nide = nest_grid%ed33 ; intermediate_grid%njde = nest_grid%sd31
         END SELECT
         intermediate_grid%nids = ids
         intermediate_grid%nide = ide
         intermediate_grid%njds = jds
         intermediate_grid%njde = jde

         intermediate_grid%sm31x                           = c_sm1x
         intermediate_grid%em31x                           = c_em1x
         intermediate_grid%sm32x                           = c_sm2x
         intermediate_grid%em32x                           = c_em2x
         intermediate_grid%sm33x                           = c_sm3x
         intermediate_grid%em33x                           = c_em3x
         intermediate_grid%sm31y                           = c_sm1y
         intermediate_grid%em31y                           = c_em1y
         intermediate_grid%sm32y                           = c_sm2y
         intermediate_grid%em32y                           = c_em2y
         intermediate_grid%sm33y                           = c_sm3y
         intermediate_grid%em33y                           = c_em3y

#if (DA_CORE == 1)
         intermediate_grid%sp31x                           = c_sp1x
         intermediate_grid%ep31x                           = c_ep1x
         intermediate_grid%sp32x                           = c_sp2x
         intermediate_grid%ep32x                           = c_ep2x
         intermediate_grid%sp33x                           = c_sp3x
         intermediate_grid%ep33x                           = c_ep3x
         intermediate_grid%sp31y                           = c_sp1y
         intermediate_grid%ep31y                           = c_ep1y
         intermediate_grid%sp32y                           = c_sp2y
         intermediate_grid%ep32y                           = c_ep2y
         intermediate_grid%sp33y                           = c_sp3y
         intermediate_grid%ep33y                           = c_ep3y
#endif

#if ( ( defined(SGIALTIX) || defined(FUJITSU_FX10) || defined(KEEP_INT_AROUND) ) && (! defined(MOVE_NESTS) ) )
         ! allocate space for the intermediate domain
!         CALL alloc_space_field ( intermediate_grid, intermediate_grid%id , 1, 2 , .TRUE., intercomm_active( intermediate_grid%id ), &   ! use same id as nest
         CALL alloc_space_field ( intermediate_grid, intermediate_grid%id , 1, 2 , .TRUE., nest_grid%active_this_task, &   ! use same id as nest
                               c_sd1, c_ed1, c_sd2, c_ed2, c_sd3, c_ed3,       &
                               c_sm1,  c_em1,  c_sm2,  c_em2,  c_sm3,  c_em3,  &
                               c_sp1,  c_ep1,  c_sp2,  c_ep2,  c_sp3,  c_ep3,  &
                               c_sm1x, c_em1x, c_sm2x, c_em2x, c_sm3x, c_em3x, &
                               c_sm1y, c_em1y, c_sm2y, c_em2y, c_sm3y, c_em3y, &
                               c_sm1x, c_em1x, c_sm2x, c_em2x, c_sm3x, c_em3x, &   ! x-xpose
                               c_sm1y, c_em1y, c_sm2y, c_em2y, c_sm3y, c_em3y  )   ! y-xpose
#endif
         intermediate_grid%sd31                            =   c_sd1
         intermediate_grid%ed31                            =   c_ed1
         intermediate_grid%sp31                            = c_sp1
         intermediate_grid%ep31                            = c_ep1
         intermediate_grid%sm31                            = c_sm1
         intermediate_grid%em31                            = c_em1
         intermediate_grid%sd32                            =   c_sd2
         intermediate_grid%ed32                            =   c_ed2
         intermediate_grid%sp32                            = c_sp2
         intermediate_grid%ep32                            = c_ep2
         intermediate_grid%sm32                            = c_sm2
         intermediate_grid%em32                            = c_em2
         intermediate_grid%sd33                            =   c_sd3
         intermediate_grid%ed33                            =   c_ed3
         intermediate_grid%sp33                            = c_sp3
         intermediate_grid%ep33                            = c_ep3
         intermediate_grid%sm33                            = c_sm3
         intermediate_grid%em33                            = c_em3

         CALL med_add_config_info_to_grid ( intermediate_grid )

         intermediate_grid%dx = parent%dx
         intermediate_grid%dy = parent%dy
         intermediate_grid%dt = parent%dt
      END IF

      RETURN
  END SUBROUTINE patch_domain_rsl_lite






  module SUBROUTINE split_communicator
#ifndef STUBMPI
      IMPLICIT NONE
      LOGICAL mpi_inited
!      INTEGER mpi_comm_here, mpi_comm_local, comdup, comdup2, origmytask,  mytask, ntasks, ierr, io_status
      INTEGER mpi_comm_here, mpi_comm_local, comdup, comdup2, origmytask,  ierr, io_status
      INTEGER mpi_comm_me_and_mom
      INTEGER coords(3)
      INTEGER mytask_local,ntasks_local,num_compute_tasks
# if defined(_OPENMP) && defined(MPI2_THREAD_SUPPORT)
      INTEGER thread_support_provided, thread_support_requested
# endif
      INTEGER i, j, k, x, y, n_x, n_y
      INTEGER iii
      INTEGER, ALLOCATABLE :: icolor(:),icolor2(:),idomain(:)
      INTEGER comm_id
!
! Communicator definition                       Domains
!                                                       
!  6 pe Example Comm   PEs                        (1)
!        COMM_WORLD  0 1 2 3 4 5                  / !               1    0 1 2 3 4 5                (2) (3)
!               2    0 1                         |
!               3        0 1 2 3                (4)
!               4    0 1
!
!    Notes: 1. No requirement that any communicator be all tasks
!           2. A task may be a member of an arbitrary number
!              of local communicators (But you may not want to do this)
!
!
!  Namelist Split Settings (for 3 comms, 4 domains)
!  Revised namelist semantics -- no need for binding nests to separately defined communicators
!
!   (domain_id)     1    2    3    4
!   parent_id       -    1    1    2
!   comm_start      0    0    2    0
!   nest_pes_x      2    1    2    1
!   nest_pes_y      3    2    2    2
!   
!! superceded
!!  Namelist Split Settings (for 3 comms, 4 domains)
!!   (comm_id)       1    2    3    ...
!!   comm_start      0    0    2
!!   comm_pes_x      2    1    2
!!   comm_pes_y      3    2    2
!!
!!  Domain definitions
!!   (domain_id)     1    2    3    4
!!   parent_id       -    1    1    2
!!   comm_domain     1    2    3    2
!! * nest_pes_x      2    1    2    1
!! * nest_pes_y      3    2    2    2
!!
!!   [* nest_pes_x is comm_pes_x(comm_domain(domain_id))]
!

      INTEGER dims(3)
! for parallel nesting, 201408, jm
      INTEGER :: id
      INTEGER :: intercomm
      INTEGER :: domain_id,par_id,nest_id,kid_id
      INTEGER :: mytask_me_and_mom, ntasks_me_and_mom, remote_leader
      LOGICAL :: inthisone
      LOGICAL :: mytask_is_nest, mytask_is_par,isperiodic(3)
! for new quilting
      LOGICAL :: quilting_is_turned_off

!!!!! needed to sneak-peek the registry to get parent_id
! define as temporaries
#include "namelist_defines.inc"

! Statements that specify the namelists
#include "namelist_statements.inc"

      CALL MPI_INITIALIZED( mpi_inited, ierr )
      IF ( .NOT. mpi_inited ) THEN
# if defined(_OPENMP) && defined(MPI2_THREAD_SUPPORT)
        thread_support_requested = MPI_THREAD_FUNNELED
        CALL mpi_init_thread ( thread_support_requested, thread_support_provided, ierr )
        IF ( thread_support_provided .lt. thread_support_requested ) THEN
           CALL WRF_ERROR_FATAL( "failed to initialize MPI thread support")
        END IF
        mpi_comm_here = MPI_COMM_WORLD
# else
#if ( DA_CORE != 1 )
        IF ( coupler_on ) THEN
           CALL cpl_init( mpi_comm_here )
        ELSE
#endif
           CALL mpi_init ( ierr )
           mpi_comm_here = MPI_COMM_WORLD
#if ( DA_CORE != 1 )
        END IF
#endif
# endif
        CALL wrf_set_dm_communicator( mpi_comm_here )
        CALL wrf_termio_dup( mpi_comm_here )
#if (WRFPLUS == 1)
      ELSE
        CALL wrf_set_dm_communicator( local_communicator )
#endif
      END IF
! this should have been reset by init_module_wrf_quilt to be just the compute tasks
      CALL wrf_get_dm_communicator( mpi_comm_here )

      CALL MPI_Comm_rank ( mpi_comm_here, mytask_local, ierr ) ;
      CALL MPI_Comm_size ( mpi_comm_here, ntasks_local, ierr ) ;
      mpi_comm_allcompute = mpi_comm_here

      IF ( mytask_local .EQ. 0 ) THEN
        max_dom = 1
        OPEN ( unit=27, file="namelist.input", form="formatted", status="old" )
        READ ( UNIT = 27 , NML = domains , IOSTAT=io_status )
        REWIND(27)
        nio_groups = 1
        nio_tasks_per_group  = 0
        poll_servers = .false.
        READ ( 27 , NML = namelist_quilt, IOSTAT=io_status )
        CLOSE(27)
      END IF
      CALL mpi_bcast( nio_tasks_per_group  , max_domains , MPI_INTEGER , 0 , mpi_comm_here, ierr )
      CALL mpi_bcast( nio_groups , 1 , MPI_INTEGER , 0 , mpi_comm_here, ierr )
      CALL mpi_bcast( max_dom, 1 , MPI_INTEGER , 0 , mpi_comm_here, ierr )
      CALL mpi_bcast( parent_id, max_domains , MPI_INTEGER , 0 , mpi_comm_here, ierr )
      CALL quilting_disabled( quilting_is_turned_off )
      IF ( quilting_is_turned_off ) THEN
        num_io_tasks = 0
        nio_tasks_per_group  = 0
        nio_groups = 1
      ELSE
        num_io_tasks = nio_tasks_per_group(1)*nio_groups
      END IF
      CALL nl_set_max_dom(1,max_dom)  ! quilting wants to see this too

      IF ( mytask_local .EQ. 0 ) THEN
        OPEN ( unit=27, file="namelist.input", form="formatted", status="old" )
! get a sneak peek an nproc_x and nproc_y
        nproc_x = -1
        nproc_y = -1
        READ ( 27 , NML = domains, IOSTAT=io_status )
        CLOSE ( 27 )
        OPEN ( unit=27, file="namelist.input", form="formatted", status="old" )
        tasks_per_split = ntasks_local
! we need to sneak-peek the parent_id namelist setting, ,which is in the "domains" section
! of the namelist.  That namelist is registry generated, so the registry-generated information
! is #included above.
        nest_pes_x = 0    ! dimensions of communicator in X and y
        nest_pes_y = 0
        IF ( nproc_x .EQ. -1 .OR. nproc_y .EQ. -1 ) THEN
          CALL compute_mesh( ntasks_local-num_io_tasks, n_x, n_y )
        ELSE
          n_x = nproc_x
          n_y = nproc_y
        END IF
        comm_start = 0   ! make it so everyone will use same communicator if the dm_task_split namelist is not specified or is empty
        nest_pes_x(1:max_dom) = n_x
        nest_pes_y(1:max_dom) = n_y
        READ ( 27 , NML = dm_task_split, IOSTAT=io_status )
        CLOSE ( 27 )
      END IF
      CALL mpi_bcast( io_status, 1 , MPI_INTEGER , 0 , mpi_comm_here, ierr )
      IF ( io_status .NE. 0 ) THEN
! or if  dm_task_split was read but was emptly, do nothing: dm_task_split not specified, everyone uses same communicator (see above)
      END IF
      CALL mpi_bcast( tasks_per_split, 1 , MPI_INTEGER , 0 , mpi_comm_here, ierr )
      CALL mpi_bcast( nproc_x, 1 , MPI_INTEGER , 0 , mpi_comm_here, ierr )
      CALL mpi_bcast( nproc_y, 1 , MPI_INTEGER , 0 , mpi_comm_here, ierr )
      CALL mpi_bcast( comm_start, max_domains , MPI_INTEGER , 0 , mpi_comm_here, ierr )
      CALL mpi_bcast( nest_pes_x, max_domains , MPI_INTEGER , 0 , mpi_comm_here, ierr )
      CALL mpi_bcast( nest_pes_y, max_domains , MPI_INTEGER , 0 , mpi_comm_here, ierr )

      nkids = 1
      which_kid = 0
      DO i = 2, max_dom
        IF ( 1 .le. parent_id(i) .AND. parent_id(i) .LE. max_domains ) THEN
          which_kid(i) = nkids(parent_id(i))
          nkids(parent_id(i)) = nkids(parent_id(i)) + 1
        ELSE
          WRITE(wrf_err_message,*)'invalid parent id for domain ',i
          CALL wrf_error_fatal(TRIM(wrf_err_message))
        END IF
      END DO

      num_compute_tasks = -99
      DO nest_id = 1,max_dom
        IF ( nest_id .EQ. 1 ) THEN
          nest_task_offsets(nest_id) = comm_start(nest_id)
        ELSE
          IF ( comm_start(nest_id) .LT. comm_start(parent_id(nest_id)) ) THEN
            WRITE(wrf_err_message,&
        "('nest domain ',i3,'comm_start (',i3,') lt parent ',i3,' comm_start (',i3,')')") &
                   nest_id,comm_start,parent_id(nest_id),comm_start(parent_id(nest_id))
            CALL wrf_error_fatal(TRIM(wrf_err_message))
          ELSE IF ( comm_start(nest_id) .LT. &
                    comm_start(parent_id(nest_id)) &
                   +nest_pes_x(parent_id(nest_id))*nest_pes_y(parent_id(nest_id))) THEN
            nest_task_offsets(nest_id) = comm_start(nest_id)-comm_start(parent_id(nest_id))
          ELSE
            nest_task_offsets(nest_id) = nest_pes_x(parent_id(nest_id))*nest_pes_y(parent_id(nest_id))
          END IF
        END IF
        IF ((comm_start(nest_id)+nest_pes_x(nest_id)*nest_pes_y(nest_id)) .GT. num_compute_tasks ) THEN
          num_compute_tasks = (comm_start(nest_id)+nest_pes_x(nest_id)*nest_pes_y(nest_id))
        END IF
      END DO

      IF ( .TRUE. ) THEN
!jm Additional code here to set up communicator for this domain and tables
!jm mapping individual domain task IDs to the original local communicator
!jm that is unsplit over nest domains.  from now on what we are calling
!jm local_communicator will be the communicator that is used by the local
!jm nests. The communicator that spans all the nests will be renamed to
!jm intercomm_communicator.  
!jm Design note: exploring the idea of using MPI intercommunicators.  They
!jm only work in pairs so we'd have a lot of intercommunicators to set up 
!jm and keep around. We'd also have to have additional communicator arguments 
!jm to all the nesting routines in and around the RSL nesting parts.
        CALL MPI_Comm_rank ( mpi_comm_here, mytask_local, ierr ) ;
        CALL MPI_Comm_rank ( mpi_comm_here, origmytask, ierr ) ;
        CALL mpi_comm_size ( mpi_comm_here, ntasks_local, ierr ) ;
        ALLOCATE( icolor(ntasks_local) )
        ALLOCATE( icolor2(ntasks_local) )
        ALLOCATE( idomain(ntasks_local) )
        k = 0
! split off the separate local communicators

! construct list of local communicators my task is in
        comms_i_am_in = MPI_UNDEFINED
        DO i = 1, max_dom
          inthisone = .FALSE.
          icolor = 0
          DO j = comm_start(i), comm_start(i)+nest_pes_x(i)*nest_pes_y(i)-1
            IF ( j+1 .GT. ntasks_local ) THEN
              WRITE(wrf_err_message,*)"check comm_start, nest_pes_x, nest_pes_y settings in namelist for comm ",i
              CALL wrf_error_fatal(wrf_err_message)
            END IF
            icolor(j+1) = 1
          END DO
          IF ( icolor(mytask_local+1) .EQ. 1 ) inthisone = .TRUE.
          CALL MPI_Comm_dup(mpi_comm_here,comdup,ierr)
          CALL MPI_Comm_split(comdup,icolor(mytask_local+1),mytask_local,mpi_comm_local,ierr)
          IF ( inthisone ) THEN
            dims(1) = nest_pes_y(i) ! rows
            dims(2) = nest_pes_x(i)  ! columns
            isperiodic(1) = .false.
            isperiodic(2) = .false.
            CALL mpi_cart_create( mpi_comm_local, 2, dims, isperiodic, .false., comms_i_am_in(i), ierr )
          END IF
        END DO

! assign domains to communicators
        local_communicator = MPI_UNDEFINED
        CALL wrf_set_dm_quilt_comm( mpi_comm_here )   ! used by module_io_quilt_old.F
        DO i = 1, max_dom
          local_communicator_store(i) = comms_i_am_in(i)
          domain_active_this_task(i) = ( local_communicator_store(i) .NE. MPI_UNDEFINED )
          IF ( local_communicator_store(i) .NE. MPI_UNDEFINED ) THEN
             CALL MPI_Comm_size( local_communicator_store(i), ntasks_store(i), ierr )
             CALL MPI_Comm_rank( local_communicator_store(i), mytask_store(i), ierr )
             CALL mpi_cart_coords( local_communicator_store(i), mytask_store(i), 2, coords, ierr )
             IF ( ierr .NE. 0 ) CALL wrf_error_fatal('MPI_cart_coords fails ')
             mytask_y_store(i) = coords(1)   ! col task (1)
             mytask_x_store(i) = coords(2)   ! col task (x)
             CALL MPI_Comm_dup( local_communicator_store(i), comdup2, ierr )
             IF ( ierr .NE. 0 ) CALL wrf_error_fatal('MPI_Comm_dup fails ')

             CALL MPI_Comm_split(comdup2,mytask_y_store(i),mytask_store(i),local_communicator_x_store(i),ierr)
             IF ( ierr .NE. 0 ) CALL wrf_error_fatal('MPI_Comm_split fails for y ')

             CALL MPI_Comm_split(comdup2,mytask_x_store(i),mytask_store(i),local_communicator_y_store(i),ierr)
             IF ( ierr .NE. 0 ) CALL wrf_error_fatal('MPI_Comm_split fails for x ')

             CALL MPI_Comm_size( local_communicator_x_store(i), ntasks_x_store(i), ierr )
             CALL MPI_Comm_rank( local_communicator_x_store(i), mytask_x_store(i), ierr )
             CALL MPI_Comm_size( local_communicator_y_store(i), ntasks_y_store(i), ierr )
             CALL MPI_Comm_rank( local_communicator_y_store(i), mytask_y_store(i), ierr )
          END IF
        END DO

        intercomm_active  = .FALSE.
        ! iterate over parent-nest pairs
        ! split off a new communicator from the big one that includes the tasks from the parent and nest communicators
        ! starting with the parent tasks followed by the nest tasks
        ! if a task is in both (ie. the communicators overlap) set the offset at the start of the first nest task
        ! in this way, we will handle cases where the parent and nest are decomposed over the same set of tasks
        ! (in that case, the offset would be the first task of the parent-nest communicator and that communicator)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ntasks_local = num_compute_tasks
        DO nest_id = 2, max_dom
           par_id  = parent_id(nest_id)
           icolor2 = 0
           DO j = 1,ntasks_local !iterate over all the tasks in the "big" communicator
             IF ( local_communicator_store( par_id ) .NE. MPI_UNDEFINED .OR. local_communicator_store( nest_id ) .NE. MPI_UNDEFINED ) icolor2(j)=1
           END DO
        ! set mpi_comm_me_and_mom to be a communicator that has my parents tasks and mine
           icolor2 = 0
           mytask_is_nest = .FALSE.
           mytask_is_par = .FALSE.
           DO j = 1,ntasks_local

             IF ( comm_start(nest_id) .LE. j-1 .AND. j-1 .LT. comm_start(nest_id) + nest_pes_x(nest_id)*nest_pes_y(nest_id) )  THEN
               icolor2(j)=1
               if ( j-1 .EQ. mytask_local ) mytask_is_nest=.TRUE.
             END IF
             IF ( comm_start(par_id ) .LE. j-1 .AND. j-1 .LT. comm_start(par_id ) + nest_pes_x(par_id )*nest_pes_y(par_id ) )  THEN
               icolor2(j)=1
               if ( j-1 .EQ. mytask_local ) mytask_is_par=.TRUE.
             END IF
           END DO

           i = icolor2(mytask_local+1)
           CALL MPI_Comm_dup(mpi_comm_here,comdup,ierr)
           CALL MPI_Comm_split(comdup,i,origmytask,mpi_comm_me_and_mom,ierr)

           IF ( mytask_is_nest  ) THEN
              intercomm_active(nest_id)  = .TRUE.
              mpi_comm_to_mom(nest_id)   =  mpi_comm_me_and_mom
           END IF
           IF ( mytask_is_par ) THEN
              intercomm_active(par_id)              = .TRUE.
              mpi_comm_to_kid(which_kid(nest_id),par_id)  =  mpi_comm_me_and_mom
           END IF
        END DO
        DEALLOCATE( icolor )
        DEALLOCATE( icolor2 )
        DEALLOCATE( idomain )

      ELSE IF ( ( tasks_per_split .LE. ntasks_local .AND. tasks_per_split .LE. 0 ) ) THEN
        domain_active_this_task(1) = .TRUE.
        IF ( mod( ntasks_local, tasks_per_split ) .NE. 0 ) THEN
          CALL wrf_message( 'WARNING: tasks_per_split does not evenly divide ntasks. Some tasks will be wasted.' )
        END IF

        ALLOCATE( icolor(ntasks_local) )
        j = 0
        DO WHILE ( j .LT. ntasks_local / tasks_per_split )
          DO i = 1, tasks_per_split
            icolor( i + j * tasks_per_split ) = j
          END DO
          j = j + 1
        END DO

        CALL MPI_Comm_dup(mpi_comm_here,comdup,ierr)
        CALL MPI_Comm_split(comdup,icolor(mytask_local+1),mytask_local,mpi_comm_local,ierr)
        CALL wrf_set_dm_communicator( mpi_comm_local )
        CALL store_communicators_for_domain(1)
        DEALLOCATE( icolor )
      ELSE
        domain_active_this_task(1) = .TRUE.
        mpi_comm_local = mpi_comm_here
        CALL wrf_set_dm_communicator( mpi_comm_local )
        CALL store_communicators_for_domain(1)
      END IF

      CALL instate_communicators_for_domain(1)

#else
! for serial (non-MPI) builds
      IMPLICIT NONE
# if defined(_OPENMP) && defined(MPI2_THREAD_SUPPORT)
      INTEGER thread_support_provided, thread_support_requested
# endif
      INTEGER i, j, k, x, y, n_x, n_y
      INTEGER iii
      INTEGER dims(3)
! for parallel nesting, 201408, jm
      INTEGER :: id
      INTEGER :: io_status
      INTEGER :: domain_id,par_id,nest_id,kid_id

!!!!! needed to sneak-peek the registry to get parent_id
! define as temporaries
#include "namelist_defines.inc"

! Statements that specify the namelists
#include "namelist_statements.inc"

      max_dom = 1
      OPEN ( unit=27, file="namelist.input", form="formatted", status="old" )
      READ ( UNIT = 27 , NML = domains , IOSTAT=io_status )
      CLOSE(27)

      nkids = 1
      which_kid = 0
      DO i = 2, max_dom
        IF ( 1 .le. parent_id(i) .AND. parent_id(i) .LE. max_domains ) THEN
          which_kid(i) = nkids(parent_id(i))
          nkids(parent_id(i)) = nkids(parent_id(i)) + 1
        ELSE
          WRITE(wrf_err_message,*)'invalid parent id for domain ',i
          CALL wrf_error_fatal(TRIM(wrf_err_message))
        END IF
      END DO

      intercomm_active = .TRUE.
      domain_active_this_task = .TRUE.
      ntasks_stack = 1
      ntasks_y_stack = 1
      ntasks_x_stack = 1
      mytask_stack = 0
      mytask_x_stack = 0
      mytask_y_stack = 0
      ntasks_store = 1
      ntasks_y_store = 1
      ntasks_x_store = 1
      mytask_store = 0
      mytask_x_store = 0
      mytask_y_store = 0
      ntasks = 1
      ntasks_y = 1
      ntasks_x = 1
      mytask = 0
      mytask_x = 0
      mytask_y = 0
      nest_pes_x = 1
      nest_pes_y = 1
      CALL instate_communicators_for_domain(1)
#endif
   END SUBROUTINE split_communicator
end submodule CODE_0