!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief methods of the rho structure (defined in qs_rho_types)
!> \par History
!>      08.2002 created [fawzi]
!>      08.2014 kpoints [JGH]
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE qs_rho_methods
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_allocate_matrix_set, cp_dbcsr_col_block_sizes, cp_dbcsr_copy, &
       cp_dbcsr_create, cp_dbcsr_deallocate_matrix_set, &
       cp_dbcsr_distribution, cp_dbcsr_init, cp_dbcsr_p_type, &
       cp_dbcsr_row_block_sizes, cp_dbcsr_set, cp_dbcsr_type, &
       dbcsr_type_symmetric
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE input_constants,                 ONLY: use_aux_fit_basis_set,&
                                             use_orb_basis_set
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE kpoint_types,                    ONLY: get_kpoint_info,&
                                             kpoint_type
  USE lri_environment_methods,         ONLY: calculate_lri_densities
  USE lri_environment_types,           ONLY: lri_density_type,&
                                             lri_environment_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_release
  USE qs_collocate_density,            ONLY: calculate_drho_elec,&
                                             calculate_rho_elec
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_local_rho_types,              ONLY: local_rho_type
  USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
  USE qs_rho_atom_methods,             ONLY: calculate_rho_atom_coeff
  USE qs_rho_types,                    ONLY: qs_rho_clear,&
                                             qs_rho_type
  USE task_list_types,                 ONLY: task_list_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_rho_methods'

  PUBLIC :: qs_rho_update_rho, qs_rho_rebuild, duplicate_rho_type

CONTAINS

! *****************************************************************************
!> \brief rebuilds rho (if necessary allocating and initializing it)
!> \param rho the rho type to rebuild (defaults to qs_env%rho)
!> \param qs_env the environment to which rho belongs
!> \param rebuild_ao if it is necessary to rebuild rho_ao. Defaults to true.
!> \param rebuild_grids if it in necessary to rebuild rho_r and rho_g.
!>        Defaults to false.
!> \param basis_set_id ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      11.2002 created replacing qs_rho_create and qs_env_rebuild_rho[fawzi]
!> \author Fawzi Mohamed
!> \note
!>      needs updated  pw pools, s, s_mstruct and h in qs_env.
!>      The use of p to keep the structure of h (needed for the forces)
!>      is ugly and should be removed.
!>      Change so that it does not allocate a subcomponent if it is not
!>      associated and not requested?
! *****************************************************************************
  SUBROUTINE qs_rho_rebuild(rho, qs_env, rebuild_ao, rebuild_grids, basis_set_id, error)
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(in), OPTIONAL            :: rebuild_ao, rebuild_grids
    INTEGER, INTENT(IN), OPTIONAL            :: basis_set_id
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_rho_rebuild', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: headline
    INTEGER                                  :: handle, i, ic, &
                                                my_basis_set_id, ncount, &
                                                nspins, stat
    LOGICAL                                  :: dokp, failure, my_rebuild_ao, &
                                                my_rebuild_grids
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      POINTER                                :: matrixkp_s
    TYPE(cp_dbcsr_type), POINTER             :: refmatrix, tmatrix
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(kpoint_type), POINTER               :: kpoints
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool

    NULLIFY (pw_env,auxbas_pw_pool,matrix_s,dft_control)
    NULLIFY (sab_orb)
    CALL timeset(routineN,handle)
    failure=.FALSE.
    my_rebuild_ao=.TRUE.
    my_rebuild_grids=.TRUE.
    IF (PRESENT(rebuild_ao)) my_rebuild_ao=rebuild_ao
    IF (PRESENT(rebuild_grids)) my_rebuild_grids=rebuild_grids

    IF (PRESENT(basis_set_id)) THEN
      my_basis_set_id = basis_set_id
    ELSE
      my_basis_set_id = use_orb_basis_set
    END IF

    CALL get_qs_env(qs_env, kpoints=kpoints, do_kpoints=dokp, error=error)

    IF(dokp) THEN
       SELECT CASE (my_basis_set_id)
       CASE (use_orb_basis_set)
         CALL get_qs_env(qs_env, sab_orb=sab_orb, error=error)
         CALL get_kpoint_info(kpoints, matrix_s=matrixkp_s, error=error)
       CASE (use_aux_fit_basis_set)
         CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END SELECT
       refmatrix => matrixkp_s(1,1)%matrix
       CALL get_kpoint_info(kpoints, ncount=ncount, error=error)
    ELSE
       SELECT CASE (my_basis_set_id)
       CASE (use_orb_basis_set)
         CALL get_qs_env(qs_env, matrix_s=matrix_s, sab_orb=sab_orb, error=error)
       CASE (use_aux_fit_basis_set)
         CALL get_qs_env(qs_env, matrix_s_aux_fit=matrix_s, sab_aux_fit=sab_orb, error=error)
       END SELECT
       refmatrix => matrix_s(1)%matrix
    END IF

    CALL get_qs_env(qs_env, pw_env=pw_env, dft_control=dft_control, error=error)

    CALL pw_env_get(pw_env,auxbas_pw_pool=auxbas_pw_pool,error=error)
    nspins=dft_control%nspins

    IF (.NOT.ASSOCIATED(rho)) STOP routineP//" rho not associated"
    IF (.NOT.ASSOCIATED(rho%tot_rho_r)) THEN
       ALLOCATE(rho%tot_rho_r(nspins))
       rho%tot_rho_r=0.0_dp
    ENDIF

    ! rho_ao
    IF (my_rebuild_ao.OR.(.NOT.ASSOCIATED(rho%rho_ao))) THEN
       IF (.NOT.failure) THEN
          IF (ASSOCIATED(rho%rho_ao)) THEN
             CALL cp_dbcsr_deallocate_matrix_set(rho%rho_ao,error=error)
          END IF
          ! Create a new density matrix set
          CALL cp_dbcsr_allocate_matrix_set(rho%rho_ao,nspins,error=error)
          DO i=1,nspins
             ALLOCATE(rho%rho_ao(i)%matrix,STAT=stat)
             CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
             tmatrix => rho%rho_ao(i)%matrix
             IF (nspins > 1) THEN
                IF (i == 1) THEN
                   headline = "DENSITY MATRIX FOR ALPHA SPIN"
                ELSE
                   headline = "DENSITY MATRIX FOR BETA SPIN"
                END IF
             ELSE
                headline = "DENSITY MATRIX"
             END IF
             CALL cp_dbcsr_init(tmatrix, error=error)
             CALL cp_dbcsr_create(matrix=tmatrix, name=TRIM(headline),&
                  dist=cp_dbcsr_distribution(refmatrix), matrix_type=dbcsr_type_symmetric,&
                  row_blk_size=cp_dbcsr_row_block_sizes(refmatrix),&
                  col_blk_size=cp_dbcsr_col_block_sizes(refmatrix),&
                  nblks=0, nze=0, error=error)
             CALL cp_dbcsr_alloc_block_from_nbl(tmatrix,sab_orb,error=error)
             CALL cp_dbcsr_set(tmatrix,0.0_dp,error=error)
          END DO
          ! kpoints
          IF(dokp) THEN
             IF (ASSOCIATED(rho%rho_ao_kp)) THEN
                CALL cp_dbcsr_deallocate_matrix_set(rho%rho_ao_kp,error=error)
             END IF
             ! Create a new density matrix set
             CALL cp_dbcsr_allocate_matrix_set(rho%rho_ao_kp,nspins,ncount,error=error)
             DO i=1,nspins
                DO ic=1,ncount
                   ALLOCATE(rho%rho_ao_kp(i,ic)%matrix,STAT=stat)
                   CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
                   tmatrix => rho%rho_ao_kp(i,ic)%matrix
                   headline = "DENSITY MATRIX (lattice vector)"
                   CALL cp_dbcsr_init(tmatrix, error=error)
                   CALL cp_dbcsr_create(matrix=tmatrix, name=TRIM(headline),&
                        dist=cp_dbcsr_distribution(refmatrix), matrix_type=dbcsr_type_symmetric,&
                        row_blk_size=cp_dbcsr_row_block_sizes(refmatrix),&
                        col_blk_size=cp_dbcsr_col_block_sizes(refmatrix),&
                        nblks=0, nze=0, error=error)
                   CALL cp_dbcsr_alloc_block_from_nbl(tmatrix,sab_orb,error=error)
                   CALL cp_dbcsr_set(tmatrix,0.0_dp,error=error)
                END DO
             END DO
          END IF
       END IF
    END IF

    ! rho_r
    IF (my_rebuild_grids.OR..NOT.ASSOCIATED(rho%rho_r)) THEN
       IF (.NOT.failure) THEN
          IF (ASSOCIATED(rho%rho_r)) THEN
             DO i=1,SIZE(rho%rho_r)
                CALL pw_release(rho%rho_r(i)%pw,error=error)
             END DO
             DEALLOCATE(rho%rho_r,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          ALLOCATE(rho%rho_r(nspins),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF (.NOT.failure) THEN
          DO i=1,nspins
             CALL pw_pool_create_pw(auxbas_pw_pool,rho%rho_r(i)%pw,&
                  use_data=REALDATA3D,in_space=REALSPACE,error=error)
          END DO
       END IF
    END IF

    ! rho_g
    IF (my_rebuild_grids.OR..NOT.ASSOCIATED(rho%rho_g)) THEN
       IF (.NOT.failure) THEN
          IF (ASSOCIATED(rho%rho_g)) THEN
             DO i=1,SIZE(rho%rho_g)
                CALL pw_release(rho%rho_g(i)%pw,error=error)
             END DO
             DEALLOCATE(rho%rho_g,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          ALLOCATE(rho%rho_g(nspins),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF (.NOT.failure) THEN
          DO i=1,nspins
             CALL pw_pool_create_pw(auxbas_pw_pool,rho%rho_g(i)%pw,&
                  use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,error=error)
          END DO
       END IF
    END IF

    ! allocate drho_r and drho_g if xc_deriv_collocate
    IF (dft_control%drho_by_collocation) THEN
       ! drho_r
       IF (my_rebuild_grids.OR..NOT.ASSOCIATED(rho%drho_r)) THEN
          IF (.NOT.failure) THEN
             IF (ASSOCIATED(rho%drho_r)) THEN
                DO i=1,SIZE(rho%drho_r)
                   CALL pw_release(rho%drho_r(i)%pw,error=error)
                END DO
                DEALLOCATE(rho%drho_r,stat=stat)
                CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
             END IF
             ALLOCATE(rho%drho_r(3*nspins),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          IF (.NOT.failure) THEN
             DO i=1,3*nspins
                CALL pw_pool_create_pw(auxbas_pw_pool,rho%drho_r(i)%pw,&
                     use_data=REALDATA3D,in_space=REALSPACE,error=error)
             END DO
          END IF
       END IF
       ! drho_g
       IF (my_rebuild_grids.OR..NOT.ASSOCIATED(rho%drho_g)) THEN
          IF (.NOT.failure) THEN
             IF (ASSOCIATED(rho%drho_g)) THEN
                DO i=1,SIZE(rho%drho_g)
                   CALL pw_release(rho%drho_g(i)%pw,error=error)
                END DO
                DEALLOCATE(rho%drho_g,stat=stat)
                CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
             END IF
             ALLOCATE(rho%drho_g(3*nspins),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          IF (.NOT.failure) THEN
             DO i=1,3*nspins
                CALL pw_pool_create_pw(auxbas_pw_pool,rho%drho_g(i)%pw,&
                     use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,error=error)
             END DO
          END IF
       END IF
    END IF

    ! allocate tau_r and tau_g if use_kinetic_energy_density
    IF (dft_control%use_kinetic_energy_density) THEN
       ! tau_r
       IF (my_rebuild_grids.OR..NOT.ASSOCIATED(rho%tau_r)) THEN
          IF (.NOT.failure) THEN
             IF (ASSOCIATED(rho%tau_r)) THEN
                DO i=1,SIZE(rho%tau_r)
                   CALL pw_release(rho%tau_r(i)%pw,error=error)
                END DO
                DEALLOCATE(rho%tau_r,stat=stat)
                CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
             END IF
             ALLOCATE(rho%tau_r(nspins),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          IF (.NOT.failure) THEN
             DO i=1,nspins
                CALL pw_pool_create_pw(auxbas_pw_pool,rho%tau_r(i)%pw,&
                     use_data=REALDATA3D,in_space=REALSPACE,error=error)
             END DO
          END IF
       END IF

       ! tau_g
       IF (my_rebuild_grids.OR..NOT.ASSOCIATED(rho%tau_g)) THEN
          IF (.NOT.failure) THEN
             IF (ASSOCIATED(rho%tau_g)) THEN
                DO i=1,SIZE(rho%tau_g)
                   CALL pw_release(rho%tau_g(i)%pw,error=error)
                END DO
                DEALLOCATE(rho%tau_g,stat=stat)
                CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
             END IF
             ALLOCATE(rho%tau_g(nspins),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          IF (.NOT.failure) THEN
             DO i=1,nspins
                CALL pw_pool_create_pw(auxbas_pw_pool,rho%tau_g(i)%pw,&
                     use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,error=error)
             END DO
          END IF
       END IF
    END IF  ! use_kinetic_energy_density

    CALL timestop(handle)
END SUBROUTINE qs_rho_rebuild

! *****************************************************************************
!> \brief updates rho_r and rho_g to the rho%rho_ao.
!>      if use_kinetic_energy_density also computes tau_r and tau_g
!> \param rho_struct the rho structure that should be updated
!> \param qs_env the qs_env rho_struct refers to
!>        the integrated charge in r space
!> \param local_rho_set ...
!> \param task_list_external ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE qs_rho_update_rho(rho_struct, qs_env, local_rho_set, task_list_external, error)
    TYPE(qs_rho_type), POINTER               :: rho_struct
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(local_rho_type), OPTIONAL, POINTER  :: local_rho_set
    TYPE(task_list_type), OPTIONAL, POINTER  :: task_list_external
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_rho_update_rho', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, nspins
    LOGICAL                                  :: failure, gapw, gapw_xc
    REAL(KIND=dp)                            :: dum
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(lri_density_type), POINTER          :: lri_density
    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(qs_rho_type), POINTER               :: rho_xc
    TYPE(task_list_type), POINTER            :: task_list

    CALL timeset(routineN,handle)

    failure=.FALSE.

    NULLIFY(dft_control,rho_xc)
    NULLIFY(lri_env, lri_density, para_env, atomic_kind_set)

    CPPrecondition(ASSOCIATED(rho_struct),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_struct%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN

      CALL get_qs_env(qs_env,dft_control=dft_control, task_list=task_list,&
                      lri_density=lri_density, lri_env=lri_env,&
                      atomic_kind_set=atomic_kind_set,  para_env=para_env,&
                      error=error)
      IF (PRESENT(task_list_external)) task_list => task_list_external 

      nspins=dft_control%nspins
      gapw = dft_control%qs_control%gapw
      gapw_xc = dft_control%qs_control%gapw_xc

      ! get tmp grids
      IF(dft_control%qs_control%semi_empirical) THEN
        !
      ELSEIF(dft_control%qs_control%dftb) THEN
        !
      ELSEIF(dft_control%qs_control%scptb) THEN
        !
      ELSEIF(dft_control%qs_control%lrigpw) THEN
        CALL calculate_lri_densities(lri_env,lri_density,qs_env,rho_struct%rho_ao,&
                                     lri_rho_struct=rho_struct,&
                                     atomic_kind_set=atomic_kind_set,&
                                     para_env=para_env,error=error)
      ELSE
        DO ispin=1,nspins
           CALL calculate_rho_elec(matrix_p=rho_struct%rho_ao(ispin)%matrix,&
                rho=rho_struct%rho_r(ispin),&
                rho_gspace=rho_struct%rho_g(ispin),&
                total_rho=rho_struct%tot_rho_r(ispin),&
                qs_env=qs_env,soft_valid=gapw,&
                task_list_external=task_list, error=error)
        END DO
      END IF
      rho_struct%rho_r_valid=.TRUE.
      rho_struct%rho_g_valid=.TRUE.

      ! if needed compute also the gradient of the density
      IF (dft_control%drho_by_collocation) THEN
         DO ispin=1,nspins
           CPPrecondition(.NOT.PRESENT(task_list_external),cp_failure_level,routineP,error,failure)
           CALL calculate_drho_elec(matrix_p=rho_struct%rho_ao(ispin)%matrix,&
             drho=rho_struct%drho_r(3*(ispin-1)+1:3*ispin),&
             drho_gspace=rho_struct%drho_g(3*(ispin-1)+1:3*ispin),&
             qs_env=qs_env,soft_valid=gapw,error=error)
         END DO
         rho_struct%drho_r_valid=.TRUE.
         rho_struct%drho_g_valid=.TRUE.
      ENDIF

      ! if needed compute also the kinetic energy density
      IF (dft_control%use_kinetic_energy_density) THEN
         IF(dft_control%qs_control%lrigpw) THEN
            CALL stop_program(routineN,moduleN,__LINE__,&
             "LRIGPW not implemented for meta-GGAs")
         ENDIF
         DO ispin=1,nspins
            CPPrecondition(.NOT.PRESENT(task_list_external),cp_failure_level,routineP,error,failure)
            CALL calculate_rho_elec(matrix_p=rho_struct%rho_ao(ispin)%matrix,&
                 rho=rho_struct%tau_r(ispin),&
                 rho_gspace=rho_struct%tau_g(ispin),&
                 total_rho=dum, &  ! presumably not meaningful
                 qs_env=qs_env, soft_valid=gapw, &
                 compute_tau=.TRUE.,&
                 error=error)
         END DO
         rho_struct%tau_r_valid=.TRUE.
         rho_struct%tau_g_valid=.TRUE.
      ENDIF

     ! GAPW o GAPW_XC require the calculation of hard and soft local densities
      IF(gapw ) THEN
         CPPrecondition(.NOT.PRESENT(task_list_external),cp_failure_level,routineP,error,failure)
         IF(PRESENT(local_rho_set)) THEN
            CALL calculate_rho_atom_coeff(qs_env,rho_struct%rho_ao,local_rho_set%rho_atom_set,error=error)
         ELSE
            CALL calculate_rho_atom_coeff(qs_env,rho_struct%rho_ao,error=error)
         ENDIF
      ENDIF
      IF(gapw_xc ) THEN
         CPPrecondition(.NOT.PRESENT(task_list_external),cp_failure_level,routineP,error,failure)
         CALL get_qs_env(qs_env=qs_env,rho_xc=rho_xc,error=error)
         CALL calculate_rho_atom_coeff(qs_env,rho_struct%rho_ao,error=error)
         ! copy rho%rho_ao into rho_xc%rho_ao
         DO ispin=1,nspins
           CALL cp_dbcsr_copy(rho_xc%rho_ao(ispin)%matrix, rho_struct%rho_ao(ispin)%matrix, error=error)
         END DO
         DO ispin=1,nspins
           CALL calculate_rho_elec(matrix_p=rho_xc%rho_ao(ispin)%matrix,&
                rho=rho_xc%rho_r(ispin),&
                rho_gspace=rho_xc%rho_g(ispin),&
                total_rho=rho_xc%tot_rho_r(ispin),&
                qs_env=qs_env,soft_valid=gapw_xc,error=error)
         END DO
         rho_xc%rho_r_valid=.TRUE.
         rho_xc%rho_g_valid=.TRUE.
         ! if needed compute also the gradient of the density
         IF (dft_control%drho_by_collocation) THEN
            DO ispin=1,nspins
              CALL calculate_drho_elec(matrix_p=rho_xc%rho_ao(ispin)%matrix,&
                drho=rho_xc%drho_r(3*(ispin-1)+1:3*ispin),&
                drho_gspace=rho_xc%drho_g(3*(ispin-1)+1:3*ispin),&
                qs_env=qs_env,soft_valid=gapw_xc,error=error)
            END DO
            rho_xc%drho_r_valid=.TRUE.
            rho_xc%drho_g_valid=.TRUE.
         ENDIF
         ! if needed compute also the kinetic energy density
         IF (dft_control%use_kinetic_energy_density) THEN
           DO ispin=1,nspins
             CALL calculate_rho_elec(matrix_p=rho_xc%rho_ao(ispin)%matrix,&
                 rho=rho_xc%tau_r(ispin),&
                 rho_gspace=rho_xc%tau_g(ispin),&
                 total_rho=dum, &  ! presumably not meaningful
                 qs_env=qs_env, soft_valid=gapw_xc, &
                 compute_tau=.TRUE., error=error)
           END DO
           rho_xc%tau_r_valid=.TRUE.
           rho_xc%tau_g_valid=.TRUE.
         ENDIF
      ENDIF
   END IF

    CALL timestop(handle)

END SUBROUTINE qs_rho_update_rho

! *****************************************************************************
!> \brief Duplicates a pointer physically
!> \param rho_input The rho structure to be duplicated
!> \param rho_output The duplicate rho structure
!> \param qs_env The QS environment from which the auxilliary PW basis-set
!>                pool is taken
!> \param error variable to control error logging, stopping,...
!>               see module cp_error_handling
!> \par History
!>      07.2005 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!> \note
!>      Associated pointers are deallocated, nullified pointers are NOT accepted!
! *****************************************************************************
  SUBROUTINE duplicate_rho_type(rho_input, rho_output, qs_env, error)

    TYPE(qs_rho_type), POINTER               :: rho_input, rho_output
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'duplicate_rho_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, nspins
    LOGICAL                                  :: failure
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool

!   ------------------------------------------------------------------------

    CALL timeset(routineN,handle)

    failure = .FALSE.
    NULLIFY(dft_control, pw_env, auxbas_pw_pool)

    CPPrecondition(ASSOCIATED(rho_input), cp_failure_level, routineP, error, failure)
    CPPrecondition(rho_input%ref_count>0, cp_failure_level, routineP, error, failure)
    CPPrecondition(ASSOCIATED(rho_output), cp_failure_level, routineP, error, failure)
    CPPrecondition(rho_output%ref_count>0, cp_failure_level, routineP, error, failure)
    CPPrecondition(ASSOCIATED(qs_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(qs_env%ref_count>0, cp_failure_level, routineP, error, failure)

    CALL get_qs_env(qs_env=qs_env, pw_env=pw_env, dft_control=dft_control, &
                    error=error)
    CALL pw_env_get(pw_env=pw_env, auxbas_pw_pool=auxbas_pw_pool, error=error)
    nspins = dft_control%nspins

    CALL qs_rho_clear(rho_output, error=error)

    ! rho_ao
    IF (ASSOCIATED(rho_input%rho_ao)) THEN
      CALL cp_dbcsr_allocate_matrix_set(rho_output%rho_ao, nspins,error=error)
      DO i = 1, nspins
         ALLOCATE(rho_output%rho_ao(i)%matrix)
         CALL cp_dbcsr_init(rho_output%rho_ao(i)%matrix, error=error)
         CALL cp_dbcsr_copy(rho_output%rho_ao(i)%matrix,rho_input%rho_ao(i)%matrix, &
              name="myDensityMatrix_for_Spin_"//TRIM(ADJUSTL(cp_to_string(i))), &
              error=error)
         CALL cp_dbcsr_set(rho_output%rho_ao(i)%matrix,0.0_dp,error=error)
      END DO
    END IF

    ! rho_r
    IF (ASSOCIATED(rho_input%rho_r)) THEN
      ALLOCATE(rho_output%rho_r(nspins))
      DO i = 1, nspins
        CALL pw_pool_create_pw(auxbas_pw_pool, rho_output%rho_r(i)%pw, &
                                use_data=REALDATA3D, in_space=REALSPACE, &
                                error=error)
        rho_output%rho_r(i)%pw%cr3d(:,:,:) = rho_input%rho_r(i)%pw%cr3d(:,:,:)
      END DO
    END IF

    ! rho_g
    IF (ASSOCIATED(rho_input%rho_g)) THEN
      ALLOCATE(rho_output%rho_g(nspins))
      DO i = 1, nspins
        CALL pw_pool_create_pw(auxbas_pw_pool, rho_output%rho_g(i)%pw, &
                                use_data=COMPLEXDATA1D, &
                                in_space=RECIPROCALSPACE, error=error)
        rho_output%rho_g(i)%pw%cc(:)      = rho_input%rho_g(i)%pw%cc(:)
      END DO
    END IF

    ! drho_r and drho_g are only needed if calculated by collocation
    IF (dft_control%drho_by_collocation) THEN
      ! drho_r
      IF (ASSOCIATED(rho_input%drho_r)) THEN
        ALLOCATE(rho_output%drho_r(3*nspins))
        DO i = 1, 3*nspins
          CALL pw_pool_create_pw(auxbas_pw_pool, rho_output%drho_r(i)%pw, &
                                  use_data=REALDATA3D, in_space=REALSPACE, &
                                  error=error)
          rho_output%drho_r(i)%pw%cr3d(:,:,:) = rho_input%drho_r(i)%pw%cr3d(:,:,:)
        END DO
      END IF

      ! drho_g
      IF (ASSOCIATED(rho_input%drho_g)) THEN
        ALLOCATE(rho_output%drho_g(3*nspins))
        DO i = 1, 3*nspins
          CALL pw_pool_create_pw(auxbas_pw_pool, rho_output%drho_g(i)%pw, &
                                  use_data=COMPLEXDATA1D, &
                                  in_space=RECIPROCALSPACE, error=error)
          rho_output%drho_g(i)%pw%cc(:) = rho_input%drho_g(i)%pw%cc(:)
        END DO
      END IF
    END IF

    ! tau_r and tau_g are only needed in the case of Meta-GGA XC-functionals
    ! are used. Therefore they are only allocated if
    ! dft_control%use_kinetic_energy_density is true
    IF (dft_control%use_kinetic_energy_density) THEN
      ! tau_r
      IF (ASSOCIATED(rho_input%tau_r)) THEN
        ALLOCATE(rho_output%tau_r(nspins))
        DO i = 1, nspins
          CALL pw_pool_create_pw(auxbas_pw_pool, rho_output%tau_r(i)%pw, &
                                  use_data=REALDATA3D, in_space=REALSPACE, &
                                  error=error)
          rho_output%tau_r(i)%pw%cr3d(:,:,:) = rho_input%tau_r(i)%pw%cr3d(:,:,:)
        END DO
      END IF

      ! tau_g
      IF (ASSOCIATED(rho_input%tau_g)) THEN
        ALLOCATE(rho_output%tau_g(nspins))
        DO i = 1, nspins
          CALL pw_pool_create_pw(auxbas_pw_pool, rho_output%tau_g(i)%pw, &
                                  use_data=COMPLEXDATA1D, &
                                  in_space=RECIPROCALSPACE, error=error)
          rho_output%tau_g(i)%pw%cc(:) = rho_input%tau_g(i)%pw%cc(:)
        END DO
      END IF
    END IF

    rho_output%rho_g_valid = rho_input%rho_g_valid
    rho_output%rho_r_valid = rho_input%rho_r_valid
    rho_output%drho_g_valid = rho_input%drho_g_valid
    rho_output%drho_r_valid = rho_input%drho_r_valid
    rho_output%tau_r_valid = rho_input%tau_r_valid
    rho_output%tau_g_valid = rho_input%tau_g_valid
    rho_output%soft_valid  = rho_input%soft_valid
    rho_output%rebuild_each = rho_input%rebuild_each

    ! tot_rho_r
    IF (ASSOCIATED(rho_input%tot_rho_r)) THEN
      ALLOCATE(rho_output%tot_rho_r(nspins))
      DO i = 1, nspins
        rho_output%tot_rho_r(i) = rho_input%tot_rho_r(i)
      END DO
    END IF

    ! tot_rho_g
    IF (ASSOCIATED(rho_input%tot_rho_g)) THEN
      ALLOCATE(rho_output%tot_rho_g(nspins))
      DO i = 1, nspins
         rho_output%tot_rho_g(i) = rho_input%tot_rho_g(i)
      END DO
    END IF

    CALL timestop(handle)

  END SUBROUTINE duplicate_rho_type

END MODULE qs_rho_methods

