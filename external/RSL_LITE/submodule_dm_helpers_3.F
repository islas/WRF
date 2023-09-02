#if ( EM_CORE == 1 && DA_CORE != 1 )

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

