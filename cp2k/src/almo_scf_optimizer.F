!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Optimization routines for all ALMO-based SCF methods
!> \par History
!>       2011.05 created [Rustam Z Khaliullin]
!>       2014.10 as a separate file [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
MODULE almo_scf_optimizer
   USE almo_scf_diis_types,             ONLY: almo_scf_diis_extrapolate,&
                                              almo_scf_diis_init,&
                                              almo_scf_diis_push,&
                                              almo_scf_diis_release,&
                                              almo_scf_diis_type
   USE almo_scf_methods,                ONLY: &
        almo_scf_ks_blk_to_tv_blk, almo_scf_ks_to_ks_blk, almo_scf_ks_to_ks_xx, &
        almo_scf_ks_xx_to_tv_xx, almo_scf_p_blk_to_t_blk, almo_scf_t_blk_to_p, &
        almo_scf_t_blk_to_t_blk_orthonormal, almo_scf_t_to_p, apply_domain_operators, &
        apply_projector, construct_domain_preconditioner, construct_domain_r_down, &
        construct_domain_s_inv, construct_domain_s_sqrt, get_overlap, newton_grad_to_step, &
        pseudo_invert_diagonal_blk
   USE almo_scf_qs,                     ONLY: almo_scf_dm_to_ks,&
                                              almo_scf_update_ks_energy,&
                                              matrix_qs_to_almo
   USE almo_scf_types,                  ONLY: almo_scf_env_type,&
                                              optimizer_options_type
   USE cp_dbcsr_cholesky,               ONLY: cp_dbcsr_cholesky_decompose,&
                                              cp_dbcsr_cholesky_invert,&
                                              cp_dbcsr_cholesky_restore
   USE cp_dbcsr_diag,                   ONLY: cp_dbcsr_syevd
   USE cp_dbcsr_interface,              ONLY: &
        cp_dbcsr_add, cp_dbcsr_add_on_diag, cp_dbcsr_copy, cp_dbcsr_create, cp_dbcsr_desymmetrize, &
        cp_dbcsr_distribution, cp_dbcsr_filter, cp_dbcsr_finalize, cp_dbcsr_frobenius_norm, &
        cp_dbcsr_function_of_elements, cp_dbcsr_get_diag, cp_dbcsr_get_info, &
        cp_dbcsr_hadamard_product, cp_dbcsr_init, cp_dbcsr_iterator, &
        cp_dbcsr_iterator_blocks_left, cp_dbcsr_iterator_next_block, cp_dbcsr_iterator_start, &
        cp_dbcsr_iterator_stop, cp_dbcsr_multiply, cp_dbcsr_norm, cp_dbcsr_p_type, &
        cp_dbcsr_print_block_sum, cp_dbcsr_release, cp_dbcsr_reserve_block2d, cp_dbcsr_scale, &
        cp_dbcsr_set, cp_dbcsr_set_diag, cp_dbcsr_trace, cp_dbcsr_triu, cp_dbcsr_type, &
        cp_dbcsr_work_create, dbcsr_distribution_mp, dbcsr_func_artanh, dbcsr_func_dtanh, &
        dbcsr_func_inverse, dbcsr_func_tanh, dbcsr_mp_numnodes, dbcsr_norm_maxabsnorm, &
        dbcsr_type_no_symmetry
   USE cp_external_control,             ONLY: external_control
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE ct_methods,                      ONLY: analytic_line_search,&
                                              ct_step_execute,&
                                              diagonalize_diagonal_blocks
   USE ct_types,                        ONLY: ct_step_env_clean,&
                                              ct_step_env_get,&
                                              ct_step_env_init,&
                                              ct_step_env_set,&
                                              ct_step_env_type
   USE domain_submatrix_methods,        ONLY: add_submatrices,&
                                              construct_submatrices,&
                                              copy_submatrices,&
                                              init_submatrices,&
                                              maxnorm_submatrices,&
                                              release_submatrices
   USE domain_submatrix_types,          ONLY: domain_submatrix_type,&
                                              select_row
   USE input_constants,                 ONLY: &
        almo_scf_diag, almo_scf_dm_sign, cg_dai_yuan, cg_fletcher, cg_fletcher_reeves, &
        cg_hager_zhang, cg_hestenes_stiefel, cg_liu_storey, cg_polak_ribiere, cg_zero, prec_zero, &
        virt_full, xalmo_case_block_diag, xalmo_case_fully_deloc, xalmo_case_normal
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE iterate_matrix,                  ONLY: invert_Hotelling,&
                                              matrix_sqrt_Newton_Schulz
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: m_flush,&
                                              m_walltime
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'almo_scf_optimizer'

   PUBLIC :: almo_scf_block_diagonal, &
             almo_scf_xalmo_eigensolver, &
             almo_scf_xalmo_pcg

   LOGICAL, PARAMETER :: debug_mode = .FALSE.
   LOGICAL, PARAMETER :: safe_mode = .FALSE.

CONTAINS

! **************************************************************************************************
!> \brief An SCF procedure that optimizes block-diagonal ALMOs using DIIS
!> \param qs_env ...
!> \param almo_scf_env ...
!> \param optimizer ...
!> \par History
!>       2011.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE almo_scf_block_diagonal(qs_env, almo_scf_env, optimizer)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(almo_scf_env_type)                            :: almo_scf_env
      TYPE(optimizer_options_type)                       :: optimizer

      CHARACTER(len=*), PARAMETER :: routineN = 'almo_scf_block_diagonal', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iscf, ispin, nspin, nstrikes, &
                                                            ntolerate, unit_nr
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: local_nocc_of_domain
      LOGICAL                                            :: converged, prepare_to_exit, should_stop, &
                                                            use_diis, use_prev_as_guess
      REAL(KIND=dp)                                      :: energy_diff, energy_new, energy_old, &
                                                            error_norm, error_norm_ispin, &
                                                            prev_error_norm, t1, t2, &
                                                            true_mixing_fraction
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: local_mu
      TYPE(almo_scf_diis_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: almo_diis
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_ks
      TYPE(cp_dbcsr_type), ALLOCATABLE, DIMENSION(:)     :: matrix_mixing_old_blk
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%mepos == logger%para_env%source) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      ! use DIIS, it's superior to simple mixing
      use_diis = .TRUE.
      use_prev_as_guess = .FALSE.

      nspin = almo_scf_env%nspins
      ALLOCATE (local_mu(almo_scf_env%ndomains))
      ALLOCATE (local_nocc_of_domain(almo_scf_env%ndomains))

      ! init mixing matrices
      ALLOCATE (matrix_mixing_old_blk(nspin))
      ALLOCATE (almo_diis(nspin))
      DO ispin = 1, nspin
         CALL cp_dbcsr_init(matrix_mixing_old_blk(ispin))
         CALL cp_dbcsr_create(matrix_mixing_old_blk(ispin), &
                              template=almo_scf_env%matrix_ks_blk(ispin))
         CALL almo_scf_diis_init(diis_env=almo_diis(ispin), &
                                 sample_err=almo_scf_env%matrix_ks_blk(ispin), &
                                 sample_var=almo_scf_env%matrix_s_blk(1), &
                                 error_type=1, &
                                 max_length=optimizer%ndiis)
      ENDDO

      energy_old = 0.0_dp
      iscf = 0
      prepare_to_exit = .FALSE.
      true_mixing_fraction = 0.0_dp
      ! set variables that control diag/pcg switching
      nstrikes = 0
      ntolerate = 3
      error_norm = 1.0E+10_dp ! arbitrary big step

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A,A,A)') REPEAT("-", 20), &
            " Optimization of block-diagonal ALMOs ", REPEAT("-", 21)
         WRITE (unit_nr, *)
         WRITE (unit_nr, '(T2,A13,A6,A23,A14,A14,A9)') "Method", "Iter", &
            "Total Energy", "Change", "Convergence", "Time"
         WRITE (unit_nr, '(T2,A)') REPEAT("-", 79)
      ENDIF

      ! the real SCF loop
      t1 = m_walltime()
      DO

         iscf = iscf+1

         ! get a copy of the current KS matrix
         CALL get_qs_env(qs_env, matrix_ks=matrix_ks)
         DO ispin = 1, nspin
            CALL matrix_qs_to_almo(matrix_ks(ispin)%matrix, &
                                   almo_scf_env%matrix_ks(ispin), &
                                   almo_scf_env, .FALSE.)
            CALL matrix_qs_to_almo(matrix_ks(ispin)%matrix, &
                                   almo_scf_env%matrix_ks_blk(ispin), &
                                   almo_scf_env, .TRUE.)
            CALL cp_dbcsr_filter(almo_scf_env%matrix_ks(ispin), &
                                 almo_scf_env%eps_filter)
         ENDDO

         ! obtain projected KS matrix and the DIIS-error vector
         CALL almo_scf_ks_to_ks_blk(almo_scf_env)

         ! inform the DIIS handler about the new KS matrix and its error vector
         IF (use_diis) THEN
            DO ispin = 1, nspin
               CALL almo_scf_diis_push(diis_env=almo_diis(ispin), &
                                       var=almo_scf_env%matrix_ks_blk(ispin), &
                                       err=almo_scf_env%matrix_err_blk(ispin))
            ENDDO
         ENDIF

         ! get error_norm: choose the largest of the two spins
         prev_error_norm = error_norm
         DO ispin = 1, nspin
            !error_norm=cp_dbcsr_frobenius_norm(almo_scf_env%matrix_err_blk(ispin))
            CALL cp_dbcsr_norm(almo_scf_env%matrix_err_blk(ispin), &
                               dbcsr_norm_maxabsnorm, &
                               norm_scalar=error_norm_ispin)
            IF (ispin .EQ. 1) error_norm = error_norm_ispin
            IF (ispin .GT. 1 .AND. error_norm_ispin .GT. error_norm) &
               error_norm = error_norm_ispin
         ENDDO

         IF (error_norm .LT. almo_scf_env%eps_prev_guess) THEN
            use_prev_as_guess = .TRUE.
         ELSE
            use_prev_as_guess = .FALSE.
         ENDIF

         ! check convergence
         converged = .TRUE.
         IF (error_norm .GT. optimizer%eps_error) converged = .FALSE.
         ! check other exit criteria: max SCF steps and timing
         CALL external_control(should_stop, "SCF", &
                               start_time=qs_env%start_time, &
                               target_time=qs_env%target_time)
         IF (should_stop .OR. iscf >= optimizer%max_iter .OR. converged) THEN
            prepare_to_exit = .TRUE.
         ENDIF

         IF (.NOT. prepare_to_exit) THEN ! update the ALMOs and density matrix

            ! perform mixing of KS matrices
            IF (iscf .NE. 1) THEN
               IF (use_diis) THEN ! use diis instead of mixing
                  DO ispin = 1, nspin
                     CALL almo_scf_diis_extrapolate(diis_env=almo_diis(ispin), &
                                                    extr_var=almo_scf_env%matrix_ks_blk(ispin))
                  ENDDO
               ELSE ! use mixing
                  true_mixing_fraction = almo_scf_env%mixing_fraction
                  DO ispin = 1, nspin
                     CALL cp_dbcsr_add(almo_scf_env%matrix_ks_blk(ispin), &
                                       matrix_mixing_old_blk(ispin), &
                                       true_mixing_fraction, &
                                       1.0_dp-true_mixing_fraction)
                  END DO
               ENDIF
            ENDIF
            ! save the new matrix for the future mixing
            DO ispin = 1, nspin
               CALL cp_dbcsr_copy(matrix_mixing_old_blk(ispin), &
                                  almo_scf_env%matrix_ks_blk(ispin))
            ENDDO

            ! obtain ALMOs from the new KS matrix
            SELECT CASE (almo_scf_env%almo_update_algorithm)
            CASE (almo_scf_diag)

               CALL almo_scf_ks_blk_to_tv_blk(almo_scf_env)

            CASE (almo_scf_dm_sign)

               ! update the density matrix
               DO ispin = 1, nspin

                  local_nocc_of_domain(:) = almo_scf_env%nocc_of_domain(:, ispin)
                  local_mu(:) = almo_scf_env%mu_of_domain(:, ispin)
                  ! RZK UPDATE! the update algorithm is removed because
                  ! RZK UPDATE! it requires updating core LS_SCF routines
                  ! RZK UPDATE! (the code exists in the CVS version)
                  CPABORT("Density_matrix_sign has not been tested yet")
                  ! RZK UPDATE!  CALL density_matrix_sign(almo_scf_env%matrix_p_blk(ispin),&
                  ! RZK UPDATE!          local_mu,&
                  ! RZK UPDATE!          almo_scf_env%fixed_mu,&
                  ! RZK UPDATE!          almo_scf_env%matrix_ks_blk(ispin),&
                  ! RZK UPDATE!          !matrix_mixing_old_blk(ispin),&
                  ! RZK UPDATE!          almo_scf_env%matrix_s_blk(1), &
                  ! RZK UPDATE!          almo_scf_env%matrix_s_blk_inv(1), &
                  ! RZK UPDATE!          local_nocc_of_domain,&
                  ! RZK UPDATE!          almo_scf_env%eps_filter,&
                  ! RZK UPDATE!          almo_scf_env%domain_index_of_ao)
                  ! RZK UPDATE!
                  almo_scf_env%mu_of_domain(:, ispin) = local_mu(:)

               ENDDO

               ! obtain ALMOs from matrix_p_blk: T_new = P_blk S_blk T_old
               CALL almo_scf_p_blk_to_t_blk(almo_scf_env, ionic=.FALSE.)
               CALL almo_scf_t_blk_to_t_blk_orthonormal(almo_scf_env)

            END SELECT

            ! obtain density matrix from ALMOs
            CALL almo_scf_t_blk_to_p(almo_scf_env, &
                                     use_sigma_inv_guess=use_prev_as_guess)

            ! compute the new KS matrix and new energy
            CALL almo_scf_dm_to_ks(qs_env, almo_scf_env, energy_new)

         ENDIF ! prepare_to_exit

         energy_diff = energy_new-energy_old
         energy_old = energy_new
         almo_scf_env%almo_scf_energy = energy_new

         t2 = m_walltime()
         ! brief report on the current SCF loop
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T2,A13,I6,F23.10,E14.5,F14.9,F9.2)') "ALMO SCF DIIS", &
               iscf, &
               energy_new, energy_diff, error_norm, t2-t1
         ENDIF
         t1 = m_walltime()

         IF (prepare_to_exit) EXIT

      ENDDO ! end scf cycle

      IF (.NOT. converged) THEN
         IF (unit_nr > 0) THEN
            CPABORT("SCF for block-diagonal ALMOs not converged!")
         ENDIF
      ENDIF

      DO ispin = 1, nspin
         CALL cp_dbcsr_release(matrix_mixing_old_blk(ispin))
         CALL almo_scf_diis_release(diis_env=almo_diis(ispin))
      ENDDO
      DEALLOCATE (almo_diis)
      DEALLOCATE (matrix_mixing_old_blk)
      DEALLOCATE (local_mu)
      DEALLOCATE (local_nocc_of_domain)

      CALL timestop(handle)

   END SUBROUTINE almo_scf_block_diagonal

! **************************************************************************************************
!> \brief An eigensolver-based SCF to optimize extended ALMOs (i.e. ALMOs on
!>        overlapping domains)
!> \param qs_env ...
!> \param almo_scf_env ...
!> \param optimizer ...
!> \par History
!>       2013.03 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE almo_scf_xalmo_eigensolver(qs_env, almo_scf_env, optimizer)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(almo_scf_env_type)                            :: almo_scf_env
      TYPE(optimizer_options_type)                       :: optimizer

      CHARACTER(len=*), PARAMETER :: routineN = 'almo_scf_xalmo_eigensolver', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iscf, ispin, nspin, unit_nr
      LOGICAL                                            :: converged, prepare_to_exit, should_stop
      REAL(KIND=dp)                                      :: denergy_tot, energy_diff, energy_new, &
                                                            energy_old, error_norm, error_norm_0, &
                                                            spin_factor, t1, t2
      REAL(KIND=dp), DIMENSION(2)                        :: denergy_spin
      TYPE(almo_scf_diis_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: almo_diis
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_ks
      TYPE(cp_dbcsr_type)                                :: matrix_p_almo_scf_converged
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(domain_submatrix_type), ALLOCATABLE, &
         DIMENSION(:, :)                                 :: submatrix_mixing_old_blk

      CALL timeset(routineN, handle)

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%mepos == logger%para_env%source) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      nspin = almo_scf_env%nspins
      IF (nspin == 1) THEN
         spin_factor = 2.0_dp
      ELSE
         spin_factor = 1.0_dp
      ENDIF

      ! RZK-warning domain_s_sqrt and domain_s_sqrt_inv do not have spin
      ! components yet (may be used later)
      ispin = 1
      CALL construct_domain_s_sqrt( &
         matrix_s=almo_scf_env%matrix_s(1), &
         subm_s_sqrt=almo_scf_env%domain_s_sqrt(:, ispin), &
         subm_s_sqrt_inv=almo_scf_env%domain_s_sqrt_inv(:, ispin), &
         dpattern=almo_scf_env%quench_t(ispin), &
         map=almo_scf_env%domain_map(ispin), &
         node_of_domain=almo_scf_env%cpu_of_domain)
      ! TRY: construct s_inv
      !CALL construct_domain_s_inv(&
      !       matrix_s=almo_scf_env%matrix_s(1),&
      !       subm_s_inv=almo_scf_env%domain_s_inv(:,ispin),&
      !       dpattern=almo_scf_env%quench_t(ispin),&
      !       map=almo_scf_env%domain_map(ispin),&
      !       node_of_domain=almo_scf_env%cpu_of_domain)

      ! construct the domain template for the occupied orbitals
      DO ispin = 1, nspin
         ! RZK-warning we need only the matrix structure, not data
         ! replace construct_submatrices with lighter procedure with
         ! no heavy communications
         CALL construct_submatrices( &
            matrix=almo_scf_env%quench_t(ispin), &
            submatrix=almo_scf_env%domain_t(:, ispin), &
            distr_pattern=almo_scf_env%quench_t(ispin), &
            domain_map=almo_scf_env%domain_map(ispin), &
            node_of_domain=almo_scf_env%cpu_of_domain, &
            job_type=select_row)
      ENDDO

      ! init mixing matrices
      ALLOCATE (submatrix_mixing_old_blk(almo_scf_env%ndomains, nspin))
      CALL init_submatrices(submatrix_mixing_old_blk)
      ALLOCATE (almo_diis(nspin))

      ! TRY: construct block-projector
      !ALLOCATE(submatrix_tmp(almo_scf_env%ndomains))
      !DO ispin=1,nspin
      !   CALL init_submatrices(submatrix_tmp)
      !   CALL construct_domain_r_down(&
      !           matrix_t=almo_scf_env%matrix_t_blk(ispin),&
      !           matrix_sigma_inv=almo_scf_env%matrix_sigma_inv(ispin),&
      !           matrix_s=almo_scf_env%matrix_s(1),&
      !           subm_r_down=submatrix_tmp(:),&
      !           dpattern=almo_scf_env%quench_t(ispin),&
      !           map=almo_scf_env%domain_map(ispin),&
      !           node_of_domain=almo_scf_env%cpu_of_domain,&
      !           filter_eps=almo_scf_env%eps_filter)
      !   CALL multiply_submatrices('N','N',1.0_dp,&
      !           submatrix_tmp(:),&
      !           almo_scf_env%domain_s_inv(:,1),0.0_dp,&
      !           almo_scf_env%domain_r_down_up(:,ispin))
      !   CALL release_submatrices(submatrix_tmp)
      !ENDDO
      !DEALLOCATE(submatrix_tmp)

      DO ispin = 1, nspin
         ! use s_sqrt since they are already properly constructed
         ! and have the same distributions as domain_err and domain_ks_xx
         CALL almo_scf_diis_init(diis_env=almo_diis(ispin), &
                                 sample_err=almo_scf_env%domain_s_sqrt(:, ispin), &
                                 error_type=1, &
                                 max_length=optimizer%ndiis)
      ENDDO

      denergy_tot = 0.0_dp
      energy_old = 0.0_dp
      iscf = 0
      prepare_to_exit = .FALSE.

      ! the SCF loop
      t1 = m_walltime()
      DO

         iscf = iscf+1

         ! get a copy of the current KS matrix
         CALL get_qs_env(qs_env, matrix_ks=matrix_ks)
         DO ispin = 1, nspin
            CALL matrix_qs_to_almo(matrix_ks(ispin)%matrix, &
                                   almo_scf_env%matrix_ks(ispin), &
                                   almo_scf_env, .FALSE.)
            CALL matrix_qs_to_almo(matrix_ks(ispin)%matrix, &
                                   almo_scf_env%matrix_ks_blk(ispin), &
                                   almo_scf_env, .TRUE.)
            CALL cp_dbcsr_filter(almo_scf_env%matrix_ks(ispin), &
                                 almo_scf_env%eps_filter)
         ENDDO

         ! obtain projected KS matrix and the DIIS-error vector
         CALL almo_scf_ks_to_ks_xx(almo_scf_env)

         ! inform the DIIS handler about the new KS matrix and its error vector
         DO ispin = 1, nspin
            CALL almo_scf_diis_push(diis_env=almo_diis(ispin), &
                                    d_var=almo_scf_env%domain_ks_xx(:, ispin), &
                                    d_err=almo_scf_env%domain_err(:, ispin))
         ENDDO

         ! check convergence
         converged = .TRUE.
         DO ispin = 1, nspin
            !error_norm=cp_dbcsr_frobenius_norm(almo_scf_env%matrix_err_blk(ispin))
            CALL cp_dbcsr_norm(almo_scf_env%matrix_err_xx(ispin), &
                               dbcsr_norm_maxabsnorm, &
                               norm_scalar=error_norm)
            CALL maxnorm_submatrices(almo_scf_env%domain_err(:, ispin), &
                                     norm=error_norm_0)
            IF (error_norm .GT. optimizer%eps_error) THEN
               converged = .FALSE.
               EXIT ! no need to check the other spin
            ENDIF
         ENDDO
         ! check other exit criteria: max SCF steps and timing
         CALL external_control(should_stop, "SCF", &
                               start_time=qs_env%start_time, &
                               target_time=qs_env%target_time)
         IF (should_stop .OR. iscf >= optimizer%max_iter .OR. converged) THEN
            prepare_to_exit = .TRUE.
         ENDIF

         IF (.NOT. prepare_to_exit) THEN ! update the ALMOs and density matrix

            ! perform mixing of KS matrices
            IF (iscf .NE. 1) THEN
               IF (.FALSE.) THEN ! use diis instead of mixing
                  DO ispin = 1, nspin
                     CALL add_submatrices( &
                        almo_scf_env%mixing_fraction, &
                        almo_scf_env%domain_ks_xx(:, ispin), &
                        1.0_dp-almo_scf_env%mixing_fraction, &
                        submatrix_mixing_old_blk(:, ispin), &
                        'N')
                     !CALL cp_dbcsr_add(almo_scf_env%matrix_ks_blk(ispin),&
                     !        matrix_mixing_old_blk(ispin),&
                     !        almo_scf_env%mixing_fraction,&
                     !        1.0_dp-almo_scf_env%mixing_fraction)
                  END DO
               ELSE
                  DO ispin = 1, nspin
                     CALL almo_scf_diis_extrapolate(diis_env=almo_diis(ispin), &
                                                    d_extr_var=almo_scf_env%domain_ks_xx(:, ispin))
                  ENDDO
               ENDIF
            ENDIF
            ! save the new matrix for the future mixing
            DO ispin = 1, nspin
               CALL copy_submatrices( &
                  almo_scf_env%domain_ks_xx(:, ispin), &
                  submatrix_mixing_old_blk(:, ispin), &
                  copy_data=.TRUE.)
            ENDDO

            ! obtain a new set of ALMOs from the updated KS matrix
            CALL almo_scf_ks_xx_to_tv_xx(almo_scf_env)

            ! update the density matrix
            DO ispin = 1, nspin

               ! save the initial density matrix (to get the perturbative energy lowering)
               IF (iscf .EQ. 1) THEN
                  CALL cp_dbcsr_init(matrix_p_almo_scf_converged)
                  CALL cp_dbcsr_create(matrix_p_almo_scf_converged, &
                                       template=almo_scf_env%matrix_p(ispin))
                  CALL cp_dbcsr_copy(matrix_p_almo_scf_converged, &
                                     almo_scf_env%matrix_p(ispin))
               ENDIF

               ! update now
               CALL almo_scf_t_to_p( &
                  t=almo_scf_env%matrix_t(ispin), &
                  p=almo_scf_env%matrix_p(ispin), &
                  eps_filter=almo_scf_env%eps_filter, &
                  orthog_orbs=.FALSE., &
                  s=almo_scf_env%matrix_s(1), &
                  sigma=almo_scf_env%matrix_sigma(ispin), &
                  sigma_inv=almo_scf_env%matrix_sigma_inv(ispin), &
                  use_guess=.TRUE.)
               CALL cp_dbcsr_scale(almo_scf_env%matrix_p(ispin), spin_factor)

               ! obtain perturbative estimate (at no additional cost)
               ! of the energy lowering relative to the block-diagonal ALMOs
               IF (iscf .EQ. 1) THEN

                  CALL cp_dbcsr_add(matrix_p_almo_scf_converged, &
                                    almo_scf_env%matrix_p(ispin), -1.0_dp, 1.0_dp)
                  CALL cp_dbcsr_trace(almo_scf_env%matrix_ks_0deloc(ispin), &
                                      matrix_p_almo_scf_converged, &
                                      denergy_spin(ispin))

                  CALL cp_dbcsr_release(matrix_p_almo_scf_converged)

                  denergy_tot = denergy_tot+denergy_spin(ispin)

                  ! RZK-warning Energy correction can be evaluated using matrix_x
                  ! as shown in the attempt below and in the PCG procedure.
                  ! Using matrix_x allows immediate decomposition of the energy
                  ! lowering into 2-body components for EDA. However, it does not
                  ! work here because the diagonalization routine does not necessarily
                  ! produce orbitals with the same sign as the block-diagonal ALMOs
                  ! Any fixes?!

                  !CALL cp_dbcsr_init(matrix_x)
                  !CALL cp_dbcsr_create(matrix_x,&
                  !        template=almo_scf_env%matrix_t(ispin))
                  !
                  !CALL cp_dbcsr_init(matrix_tmp_no)
                  !CALL cp_dbcsr_create(matrix_tmp_no,&
                  !        template=almo_scf_env%matrix_t(ispin))
                  !
                  !CALL cp_dbcsr_copy(matrix_x,&
                  !        almo_scf_env%matrix_t_blk(ispin))
                  !CALL cp_dbcsr_add(matrix_x,almo_scf_env%matrix_t(ispin),&
                  !        -1.0_dp,1.0_dp)

                  !CALL cp_dbcsr_trace(matrix_x,&
                  !        almo_scf_env%matrix_err_xx(ispin),denergy,"T","N")

                  !denergy=denergy*spin_factor

                  !IF (unit_nr>0) THEN
                  !   WRITE(unit_nr,*) "_ENERGY-0: ", almo_scf_env%almo_scf_energy
                  !   WRITE(unit_nr,*) "_ENERGY-D: ", denergy
                  !   WRITE(unit_nr,*) "_ENERGY-F: ", almo_scf_env%almo_scf_energy+denergy
                  !ENDIF
                  !! RZK-warning update will not work since the energy is overwritten almost immediately
                  !!CALL almo_scf_update_ks_energy(qs_env,&
                  !!        almo_scf_env%almo_scf_energy+denergy)
                  !!

                  !! print out the results of the decomposition analysis
                  !CALL cp_dbcsr_hadamard_product(matrix_x,&
                  !        almo_scf_env%matrix_err_xx(ispin),&
                  !        matrix_tmp_no)
                  !CALL cp_dbcsr_scale(matrix_tmp_no,spin_factor)
                  !CALL cp_dbcsr_filter(matrix_tmp_no,almo_scf_env%eps_filter)
                  !
                  !IF (unit_nr>0) THEN
                  !   WRITE(unit_nr,*)
                  !   WRITE(unit_nr,'(T2,A)') "DECOMPOSITION OF THE DELOCALIZATION ENERGY"
                  !ENDIF

                  !mynode=dbcsr_mp_mynode(dbcsr_distribution_mp(&
                  !   cp_dbcsr_distribution(matrix_tmp_no)))
                  !WRITE(mynodestr,'(I6.6)') mynode
                  !mylogfile='EDA.'//TRIM(ADJUSTL(mynodestr))
                  !OPEN (iunit,file=mylogfile,status='REPLACE')
                  !CALL cp_dbcsr_print_block_sum(matrix_tmp_no,iunit)
                  !CLOSE(iunit)
                  !
                  !CALL cp_dbcsr_release(matrix_tmp_no)
                  !CALL cp_dbcsr_release(matrix_x)

               ENDIF ! iscf.eq.1

            ENDDO

            ! print out the energy lowering
            IF (iscf .EQ. 1) THEN
               IF (unit_nr > 0) THEN
                  WRITE (unit_nr, *)
                  WRITE (unit_nr, '(T2,A35,F25.10)') "ENERGY OF BLOCK-DIAGONAL ALMOs:", &
                     almo_scf_env%almo_scf_energy
                  WRITE (unit_nr, '(T2,A35,F25.10)') "ENERGY LOWERING:", &
                     denergy_tot
                  WRITE (unit_nr, '(T2,A35,F25.10)') "CORRECTED ENERGY:", &
                     almo_scf_env%almo_scf_energy+denergy_tot
                  WRITE (unit_nr, *)
               ENDIF
               CALL almo_scf_update_ks_energy(qs_env, &
                                              almo_scf_env%almo_scf_energy+denergy_tot)
            ENDIF

            ! compute the new KS matrix and new energy
            IF (.NOT. almo_scf_env%perturbative_delocalization) THEN
               CALL almo_scf_dm_to_ks(qs_env, almo_scf_env, energy_new)
            ENDIF

         ENDIF ! prepare_to_exit

         IF (almo_scf_env%perturbative_delocalization) THEN

            ! exit after the first step if we do not need the SCF procedure
            converged = .TRUE.
            prepare_to_exit = .TRUE.

         ELSE ! not a perturbative treatment

            energy_diff = energy_new-energy_old
            energy_old = energy_new
            almo_scf_env%almo_scf_energy = energy_new

            t2 = m_walltime()
            ! brief report on the current SCF loop
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, '(T2,A,I6,F20.9,E11.3,E11.3,E11.3,F8.2)') "ALMO SCF", &
                  iscf, &
                  energy_new, energy_diff, error_norm, error_norm_0, t2-t1
            ENDIF
            t1 = m_walltime()

         ENDIF

         IF (prepare_to_exit) EXIT

      ENDDO ! end scf cycle

      IF (.NOT. converged) THEN
         CPABORT("SCF for ALMOs on overlapping domains not converged! ")
      ENDIF

      DO ispin = 1, nspin
         CALL release_submatrices(submatrix_mixing_old_blk(:, ispin))
         CALL almo_scf_diis_release(diis_env=almo_diis(ispin))
      ENDDO
      DEALLOCATE (almo_diis)
      DEALLOCATE (submatrix_mixing_old_blk)

      CALL timestop(handle)

   END SUBROUTINE almo_scf_xalmo_eigensolver

! **************************************************************************************************
!> \brief Optimization of ALMOs using PCG-like minimizers
!> \param qs_env ...
!> \param almo_scf_env ...
!> \param optimizer   controls the optimization algorithm
!> \param quench_t ...
!> \param matrix_t_in ...
!> \param matrix_t_out ...
!> \param assume_t0_q0x ...
!> \param perturbation_only ...
!> \param special_case   to reduce the overhead special cases are implemented:
!>                       xalmo_case_normal - no special case (i.e. xALMOs)
!>                       xalmo_case_block_diag
!>                       xalmo_case_fully_deloc
!> \par History
!>       2011.11 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE almo_scf_xalmo_pcg(qs_env, almo_scf_env, optimizer, quench_t, &
                                 matrix_t_in, matrix_t_out, assume_t0_q0x, perturbation_only, &
                                 special_case)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(almo_scf_env_type)                            :: almo_scf_env
      TYPE(optimizer_options_type)                       :: optimizer
      TYPE(cp_dbcsr_type), ALLOCATABLE, DIMENSION(:)     :: quench_t, matrix_t_in, matrix_t_out
      LOGICAL, INTENT(IN)                                :: assume_t0_q0x, perturbation_only
      INTEGER, INTENT(IN), OPTIONAL                      :: special_case

      CHARACTER(len=*), PARAMETER :: routineN = 'almo_scf_xalmo_pcg', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=20)                                  :: iter_type
      INTEGER :: cg_iteration, dim0, eda_unit, fixed_line_search_niter, handle, ispin, iteration, &
         jj, line_search_iteration, max_iter, my_special_case, ncores, ndomains, occ1, &
         outer_iteration, outer_max_iter, prec_type, precond_domain_projector, unit_nr, &
         zero_neg_eiv
      LOGICAL :: converged, do_md, first_md_iteration, just_started, line_search, &
         md_in_theta_space, optimize_theta, outer_prepare_to_exit, prepare_to_exit, &
         reset_conjugator, skip_grad, use_guess, use_preconditioner
      REAL(kind=dp) :: appr_sec_der, beta, denom, e0, e1, energy_diff, energy_new, energy_old, &
         eps_skip_gradients, g0, g1, grad_norm, grad_norm_frob, kappa, kin_energy, &
         line_search_error, next_step_size_guess, numer, prec_sf_mixing_s, spin_factor, step_size, &
         t1, t2, t_norm, tau, time_step
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: evals
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_ks
      TYPE(cp_dbcsr_type) :: FTsiginv, fvo_0, grad, inv_eiv, m_theta, m_tmp_nn_1, m_tmp_no_1, &
         m_tmp_no_2, m_tmp_no_3, m_tmp_oo_1, matrix_sigma_inv_0, &
         matrix_t_0, prec_oo, prec_oo_inv, prec_vv, prev_grad, prev_minus_prec_grad, prev_step, &
         siginvTFTsiginv, ST, step, STsiginv_0, velocity
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(domain_submatrix_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: domain_r_down
      TYPE(section_vals_type), POINTER                   :: almo_print_section, input

      CALL timeset(routineN, handle)

      my_special_case = xalmo_case_normal
      IF (PRESENT(special_case)) my_special_case = special_case

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%mepos == logger%para_env%source) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, *)
         SELECT CASE (my_special_case)
         CASE (xalmo_case_block_diag)
            WRITE (unit_nr, '(T2,A,A,A)') REPEAT("-", 20), &
               " Optimization of block-diagonal ALMOs ", REPEAT("-", 21)
         CASE (xalmo_case_fully_deloc)
            WRITE (unit_nr, '(T2,A,A,A)') REPEAT("-", 20), &
               " Optimization of fully delocalized MOs ", REPEAT("-", 20)
         CASE (xalmo_case_normal)
            WRITE (unit_nr, '(T2,A,A,A)') REPEAT("-", 27), &
               " Optimization of XALMOs ", REPEAT("-", 28)
         END SELECT
         WRITE (unit_nr, *)
         WRITE (unit_nr, '(T2,A13,A6,A23,A14,A14,A9)') "Method", "Iter", &
            "Objective Function", "Change", "Convergence", "Time"
         WRITE (unit_nr, '(T2,A)') REPEAT("-", 79)
      ENDIF

      ! set local parameters using developer's keywords
      ! RZK-warning: change to normal keywords later
      do_md = almo_scf_env%logical01
      optimize_theta = almo_scf_env%logical05
      prec_sf_mixing_s = almo_scf_env%real04
      eps_skip_gradients = almo_scf_env%real01

      ! preconditioner control
      use_preconditioner = optimizer%preconditioner .NE. prec_zero
      prec_type = 4
      ! RZK-warning: prec_type here is not the same as preconditioner
      ! type in optimizer%preconditioner. change this later
      !prec_type = optimizer%preconditioner
      !if (prec_type.eq.prec_default) prec_type=prec_ks_plus_s

      ! control of the line search
      fixed_line_search_niter = 0 ! init to zero, change when eps is small enough

      ncores = dbcsr_mp_numnodes(dbcsr_distribution_mp( &
                                 cp_dbcsr_distribution(almo_scf_env%matrix_s(1))))

      ! Since it is extrememly difficult to converge the iterative
      ! procedure using T as an optimized variable, assume
      ! T = T_0 + (1-R_0)*X and optimize X for both perturbative
      ! and non-perturbative cases

      IF (almo_scf_env%nspins == 1) THEN
         spin_factor = 2.0_dp
      ELSE
         spin_factor = 1.0_dp
      ENDIF

      !!!!!! RZK-warning THIS PROCEDURE WILL WORK ONLY FOR CLOSED SHELL SYSTEMS
      !!!!!! TO ADAPT IT FOR UNRESTRICTED ORBITALS - UPDATE KS MATRIX WITH PARTIALLY
      !!!!!! OPTIMIZED ORBITALS - BOTH ALPNA AND BETA
      IF (almo_scf_env%nspins .GT. 1) THEN
         CPABORT("UNRESTRICTED ALMO SCF IS NYI(!)")
      ENDIF

      nspins = almo_scf_env%nspins

      ! init matrices for the initial guess
      DO ispin = 1, nspins
         CALL cp_dbcsr_init(m_theta)
         CALL cp_dbcsr_create(m_theta, &
                              template=matrix_t_out(ispin), &
                              matrix_type=dbcsr_type_no_symmetry)
      ENDIF

      ! create initial guess from the initial orbitals
      CALL xalmo_initial_guess(m_guess=m_theta, 
              m_t0 = almo_scf_env%matrix_t_in, &
              m_quench_t = quench_t, &
              m_overlap = almo_scf_env%matrix_s(1), &
              m_sigma_tmpl = almo_scf_env%matrix_sigma_inv, &
              nspins = nspins, &
              xalmo_history = almo_scf_env%xalmo_history, &
              assume_t0_q0x = assume_t0_q0x, &
              optimize_theta = optimize_theta, &
              envelope_amplitude = almo_scf_env%envelope_amplitude, &
              eps_filter = almo_scf_enf%eps_filter, &
              order_lanczos=almo_scf_env%order_lanczos,&
              eps_lanczos=almo_scf_env%eps_lanczos,&
              max_iter_lanczos=almo_scf_env%max_iter_lanczos)

      DO ispin = 1, nspins

         ! init temporary storage
         CALL cp_dbcsr_init(prec_vv)
         CALL cp_dbcsr_init(fvo_0)
         CALL cp_dbcsr_init(STsiginv_0)
         CALL cp_dbcsr_init(m_tmp_no_1)
         CALL cp_dbcsr_init(m_tmp_no_2)
         CALL cp_dbcsr_init(m_tmp_no_3)
         CALL cp_dbcsr_init(ST)
         CALL cp_dbcsr_init(FTsiginv)
         CALL cp_dbcsr_init(m_tmp_oo_1)
         CALL cp_dbcsr_init(m_tmp_nn_1)
         CALL cp_dbcsr_init(siginvTFTsiginv)
         CALL cp_dbcsr_init(prec_oo)
         CALL cp_dbcsr_init(prec_oo_inv)
         CALL cp_dbcsr_init(prev_grad)
         CALL cp_dbcsr_init(prev_step)
         CALL cp_dbcsr_init(grad)
         CALL cp_dbcsr_init(step)
         CALL cp_dbcsr_init(prev_minus_prec_grad)
         CALL cp_dbcsr_create(prec_vv, &
                              template=almo_scf_env%matrix_ks(ispin), &
                              matrix_type=dbcsr_type_no_symmetry)
         CALL cp_dbcsr_create(prec_oo, &
                              template=almo_scf_env%matrix_sigma(ispin), &
                              matrix_type=dbcsr_type_no_symmetry)
         CALL cp_dbcsr_create(prec_oo_inv, &
                              template=almo_scf_env%matrix_sigma(ispin), &
                              matrix_type=dbcsr_type_no_symmetry)
         CALL cp_dbcsr_create(m_tmp_oo_1, &
                              template=almo_scf_env%matrix_sigma(ispin), &
                              matrix_type=dbcsr_type_no_symmetry)
         CALL cp_dbcsr_create(siginvTFTsiginv, &
                              template=almo_scf_env%matrix_sigma(ispin), &
                              matrix_type=dbcsr_type_no_symmetry)
         CALL cp_dbcsr_create(STsiginv_0, &
                              template=matrix_t_out(ispin), &
                              matrix_type=dbcsr_type_no_symmetry)
         CALL cp_dbcsr_create(fvo_0, &
                              template=matrix_t_out(ispin), &
                              matrix_type=dbcsr_type_no_symmetry)
         CALL cp_dbcsr_create(m_tmp_no_1, &
                              template=matrix_t_out(ispin), &
                              matrix_type=dbcsr_type_no_symmetry)
         CALL cp_dbcsr_create(m_tmp_no_2, &
                              template=matrix_t_out(ispin), &
                              matrix_type=dbcsr_type_no_symmetry)
         CALL cp_dbcsr_create(m_tmp_no_3, &
                              template=matrix_t_out(ispin), &
                              matrix_type=dbcsr_type_no_symmetry)
         CALL cp_dbcsr_create(FTsiginv, &
                              template=matrix_t_out(ispin), &
                              matrix_type=dbcsr_type_no_symmetry)
         CALL cp_dbcsr_create(ST, &
                              template=matrix_t_out(ispin), &
                              matrix_type=dbcsr_type_no_symmetry)
         CALL cp_dbcsr_create(prev_grad, &
                              template=matrix_t_out(ispin), &
                              matrix_type=dbcsr_type_no_symmetry)
         CALL cp_dbcsr_create(grad, &
                              template=matrix_t_out(ispin), &
                              matrix_type=dbcsr_type_no_symmetry)
         CALL cp_dbcsr_create(prev_step, &
                              template=matrix_t_out(ispin), &
                              matrix_type=dbcsr_type_no_symmetry)
         CALL cp_dbcsr_create(step, &
                              template=matrix_t_out(ispin), &
                              matrix_type=dbcsr_type_no_symmetry)
         CALL cp_dbcsr_create(prev_minus_prec_grad, &
                              template=matrix_t_out(ispin), &
                              matrix_type=dbcsr_type_no_symmetry)

         ndomains = almo_scf_env%ndomains
         ALLOCATE (domain_r_down(ndomains))
         CALL init_submatrices(domain_r_down)

         ! create matrices to store the initial state
         CALL cp_dbcsr_init(matrix_t_0)
         CALL cp_dbcsr_init(matrix_sigma_inv_0)
         !CALL cp_dbcsr_init(matrix_sigma_0)
         !CALL cp_dbcsr_init(matrix_p_0)
         CALL cp_dbcsr_create(matrix_t_0, &
                              template=matrix_t_out(ispin), &
                              matrix_type=dbcsr_type_no_symmetry)
         CALL cp_dbcsr_create(matrix_sigma_inv_0, &
                              template=almo_scf_env%matrix_sigma_inv(ispin), &
                              matrix_type=dbcsr_type_no_symmetry)
         !CALL cp_dbcsr_create(matrix_sigma_0, &
         !                     template=almo_scf_env%matrix_sigma_inv(ispin), &
         !                     matrix_type=dbcsr_type_no_symmetry)
         !CALL cp_dbcsr_create(matrix_p_0, &
         !                     template=almo_scf_env%matrix_ks(ispin))

         CALL cp_dbcsr_set(step, 0.0_dp)

         md_in_theta_space = .FALSE. ! turn on later after several minimization steps
         IF (do_md) THEN
            CALL cp_dbcsr_init(velocity)
            CALL cp_dbcsr_create(velocity, &
                                 template=matrix_t_out(ispin))
            CALL cp_dbcsr_copy(velocity, quench_t(ispin))
            CALL cp_dbcsr_set(velocity, 0.0_dp)
            CALL cp_dbcsr_copy(prev_step, quench_t(ispin))
            CALL cp_dbcsr_set(prev_step, 0.0_dp)
            time_step = optimizer%lin_search_step_size_guess
         ENDIF

         ! invert S domains if necessary
         ! RZK-warning must be done outside the spin loop to save time
         IF (my_special_case .EQ. xalmo_case_normal) THEN
            CALL construct_domain_s_inv( &
               matrix_s=almo_scf_env%matrix_s(1), &
               subm_s_inv=almo_scf_env%domain_s_inv(:, ispin), &
               dpattern=quench_t(ispin), &
               map=almo_scf_env%domain_map(ispin), &
               node_of_domain=almo_scf_env%cpu_of_domain)
         ENDIF

         ! save the initial state and S*T_0*siginv_0
         CALL cp_dbcsr_copy(matrix_t_0, matrix_t_in(ispin))
         CALL cp_dbcsr_desymmetrize(almo_scf_env%matrix_sigma_inv(ispin), &
                                    matrix_sigma_inv_0)
         !IF (assume_t0_q0x.AND. special_case .EQ. xalmo_case_fully_deloc) THEN
            CALL cp_dbcsr_multiply("N", "N", 1.0_dp, &
                                   almo_scf_env%matrix_s(1), &
                                   matrix_t_in(ispin), &
                                   0.0_dp, ST, &
                                   filter_eps=almo_scf_env%eps_filter)
            CALL cp_dbcsr_multiply("N", "N", 1.0_dp, &
                                   ST, &
                                   matrix_sigma_inv_0, &
                                   0.0_dp, STsiginv_0, &
                                   filter_eps=almo_scf_env%eps_filter)
         !ENDIF

         ! construct domain-projector
         IF (my_special_case .EQ. xalmo_case_normal .AND. prec_type .EQ. 4) THEN
            IF (assume_t0_q0x) THEN
               CALL construct_domain_r_down( &
                  matrix_t=matrix_t_0, &
                  matrix_sigma_inv=matrix_sigma_inv_0, &
                  matrix_s=almo_scf_env%matrix_s(1), &
                  subm_r_down=domain_r_down(:), &
                  dpattern=quench_t(ispin), &
                  map=almo_scf_env%domain_map(ispin), &
                  node_of_domain=almo_scf_env%cpu_of_domain, &
                  filter_eps=almo_scf_env%eps_filter)
            ENDIF ! assume_t0_q0x
         ENDIF

         ! construct all necessary zero-order matrices
         ! note that the guess orbitals are not necessarily t_0
         IF ( perturbation_only ) THEN
            ! compute the density matrix
            CALL almo_scf_t_to_p( &
               t=matrix_t_0, &
               p=almo_scf_env%matrix_p(ispin), &
               eps_filter=almo_scf_env%eps_filter, &
               orthog_orbs=.FALSE., &
               s=almo_scf_env%matrix_s(1), &
               sigma=almo_scf_env%matrix_sigma(ispin), &
               sigma_inv=almo_scf_env%matrix_sigma_inv(ispin), &
               use_guess=.FALSE.)
            CALL cp_dbcsr_scale(almo_scf_env%matrix_p(ispin), &
                                spin_factor)
            ! RZK-warning part of the procedure must be outside the spin loop
            CALL almo_scf_dm_to_ks(qs_env, almo_scf_env, energy_new)
            CALL get_qs_env(qs_env, matrix_ks=matrix_ks)
            CALL matrix_qs_to_almo(matrix_ks(ispin)%matrix, &
                                   almo_scf_env%matrix_ks(ispin), &
                                   almo_scf_env, .FALSE.)
            CALL cp_dbcsr_filter(almo_scf_env%matrix_ks(ispin), &
                                 almo_scf_env%eps_filter)
compute initial gradient to save as fvo0                  
            CALL compute_gradient(&
                   m_grad_out=grad(ispin),&
                   m_tmp_no_2=m_tmp_no_2,&
                    m_ks=almo_scf_env%matrix_ks(ispin),&
                    m_s=almo_scf_env%matrix_s(1),&
                   m_t=matrix_t_out(ispin),&
                    m_t0=matrix_t_0,&
                    m_siginv=almo_scf_env%matrix_sigma_inv(ispin),&
                    m_quench_t=quench_t(ispin),&
                   m_FTsiginv=FTsiginv(ispin),&
                   m_siginvTFTsiginv=siginvTFTsiginv(ispin),&
                    m_ST=ST(ispin),&
                    m_STsiginv0=STsiginv_0,&
                    m_theta=m_theta,&
                    !m_sig_sqrti_ii=m_sig_sqrti_ii(ispin),&
                    domain_s_inv=almo_scf_env%domain_s_inv(:,ispin),&
                    domain_r_down=domain_r_down(:,ispin),&
                    cpu_of_domain=almo_scf_env%cpu_of_domain,&
                    domain_map=almo_scf_env%domain_map(ispin),&
                    assume_t0_q0x=assume_t0_q0x,&
                    optimize_theta=optimize_theta,&
                    perturbation_only=perturbation_only,&
                    normalize_orbitals=.FALSE.,&
                    penalty_occ_vol=.FALSE.,&
                    penalty_occ_vol_prefactor=0.0_dp,&
                    !normalize_orbitals=normalize_orbitals,&
                    !penalty_occ_vol=penalty_occ_vol,&
                    !penalty_occ_vol_prefactor=penalty_occ_vol_g_prefactor,&
                    envelope_amplitude=almo_scf_env%envelope_amplitude,&
                    eps_filter=almo_scf_env%eps_filter,&
                    spin_factor=spin_factor,&
                    special_case=my_special_case)

                  !NEW    ! save the first gradient
                  !NEW    ! it is equal to F_vo and necessary to compute
                  !NEW    ! the correction to the energy
                  NEW    CALL cp_dbcsr_copy(fvo_0, m_tmp_no_2)
                  !NEW    CALL cp_dbcsr_scale(fvo_0, 0.5_dp/spin_factor)
         ENDIF


         ! start the outer SCF loop
         outer_max_iter = optimizer%max_iter_outer_loop
         outer_prepare_to_exit = .FALSE.
         outer_iteration = 0
         grad_norm = 0.0_dp
         grad_norm_frob = 0.0_dp
         use_guess = .FALSE.

         DO

            ! start the inner SCF loop
            max_iter = optimizer%max_iter
            prepare_to_exit = .FALSE.
            line_search = .FALSE.
            converged = .FALSE.
            iteration = 0
            cg_iteration = 0
            line_search_iteration = 0
            energy_new = 0.0_dp
            energy_old = 0.0_dp
            line_search_error = 0.0_dp
            t1 = m_walltime()

            DO

               just_started = (iteration .EQ. 0) .AND. (outer_iteration .EQ. 0)

               ! switch to MD after several minimization steps
               IF (iteration .EQ. almo_scf_env%integer01 .AND. do_md) THEN
                  CALL cp_dbcsr_set(velocity, 0.0_dp)
                  CALL cp_dbcsr_set(prev_step, 0.0_dp)
                  md_in_theta_space = .TRUE.
                  first_md_iteration = .TRUE.
               ENDIF

               ! compute the MO coefficients from theta
               !NEW IF (assume_t0_q0x .AND. just_started) THEN
               !NEW    CALL cp_dbcsr_set(matrix_t_out(ispin), 0.0_dp)
               !NEW ENDIF
               IF (optimize_theta) THEN
                  CALL cp_dbcsr_norm(m_theta, &
                                     dbcsr_norm_maxabsnorm, norm_scalar=t_norm)
                  !IF (unit_nr>0) THEN
                  !   WRITE(unit_nr,*) "VAR1: |T|: ", t_norm
                  !   !WRITE(unit_nr,*) "VAR1: SIN(|T|/A), COS(|T|/A), A*SIN(|T|/A): ",&
                  !   !   SIN(t_norm/almo_scf_env%envelope_amplitude),&
                  !   !   COS(t_norm/almo_scf_env%envelope_amplitude),&
                  !   !   almo_scf_env%envelope_amplitude*SIN(t_norm/almo_scf_env%envelope_amplitude)
                  !   WRITE(unit_nr,*) "VAR1: TANH(|T|/A), D(A*TANH(|T|/A)), A*TANH(|T|/A): ",&
                  !      TANH(t_norm/almo_scf_env%envelope_amplitude),&
                  !      1.0_dp-(TANH(t_norm/almo_scf_env%envelope_amplitude))**2,&
                  !      almo_scf_env%envelope_amplitude*TANH(t_norm/almo_scf_env%envelope_amplitude)
                  !ENDIF
               ENDIF
               IF (optimize_theta) THEN
                  CALL cp_dbcsr_copy(m_tmp_no_1, m_theta)
                  CALL cp_dbcsr_function_of_elements(m_tmp_no_1, &
                                                     !func=dbcsr_func_sin,&
                                                     func=dbcsr_func_tanh, &
                                                     a0=0.0_dp, &
                                                     a1=1.0_dp/almo_scf_env%envelope_amplitude)
                  CALL cp_dbcsr_scale(m_tmp_no_1, &
                                      almo_scf_env%envelope_amplitude)
               ELSE
                  CALL cp_dbcsr_copy(m_tmp_no_1, m_theta)
               ENDIF
               CALL cp_dbcsr_hadamard_product(m_tmp_no_1, &
                                              quench_t(ispin), &
                                              matrix_t_out(ispin))
               CALL cp_dbcsr_norm(matrix_t_out(ispin), &
                                  dbcsr_norm_maxabsnorm, norm_scalar=t_norm)

               !IF (optimize_theta) THEN
               !   IF (unit_nr>0) THEN
               !      WRITE(unit_nr,*) "VAR1: |B*A*TANH(T/A)|: ", t_norm
               !   ENDIF
               !ELSE
               !   IF (unit_nr>0) THEN
               !      WRITE(unit_nr,*) "VAR1: |B*X|: ", t_norm
               !   ENDIF
               !ENDIF

               ! project out R_0
               ! NEW IF (assume_t0_q0x .AND. (.NOT. just_started)) THEN
               IF (assume_t0_q0x) THEN
                  IF (my_special_case .EQ. xalmo_case_fully_deloc) THEN
                     CALL cp_dbcsr_multiply("T", "N", 1.0_dp, &
                                            STsiginv_0, &
                                            matrix_t_out(ispin), &
                                            0.0_dp, m_tmp_oo_1, &
                                            filter_eps=almo_scf_env%eps_filter)
                     CALL cp_dbcsr_multiply("N", "N", -1.0_dp, &
                                            matrix_t_0, &
                                            m_tmp_oo_1, &
                                            1.0_dp, matrix_t_out(ispin), &
                                            filter_eps=almo_scf_env%eps_filter)
                  ELSE IF (my_special_case .EQ. xalmo_case_block_diag) THEN
                     CPABORT("cannot use projector with block-daigonal ALMOs")
                  ELSE
                     ! no special case
                     CALL apply_domain_operators( &
                        matrix_in=matrix_t_out(ispin), &
                        matrix_out=m_tmp_no_1, &
                        operator1=domain_r_down(:), &
                        operator2=almo_scf_env%domain_s_inv(:, ispin), &
                        dpattern=quench_t(ispin), &
                        map=almo_scf_env%domain_map(ispin), &
                        node_of_domain=almo_scf_env%cpu_of_domain, &
                        my_action=1, &
                        filter_eps=almo_scf_env%eps_filter, &
                        !matrix_trimmer=,&
                        use_trimmer=.FALSE.)
                     CALL cp_dbcsr_copy(matrix_t_out(ispin), &
                                        m_tmp_no_1)
                  ENDIF ! special case
                  !NEW CALL cp_dbcsr_norm(matrix_t_out(ispin), &
                  !NEW                    dbcsr_norm_maxabsnorm, norm_scalar=t_norm)
                  !IF (unit_nr>0) THEN
                  !   WRITE(unit_nr,*) "VAR1: |Q.X|: ", t_norm
                  !ENDIF
               ENDIF ! end assume_t0_q0x

               IF (assume_t0_q0x) THEN
                  CALL cp_dbcsr_add(matrix_t_out(ispin), &
                                    matrix_t_0, 1.0_dp, 1.0_dp)
                  CALL cp_dbcsr_norm(matrix_t_out(ispin), &
                                     dbcsr_norm_maxabsnorm, norm_scalar=t_norm)
                  !IF (unit_nr>0) THEN
                  !   WRITE(unit_nr,*) "VAR1: |T0+Q.X|: ", t_norm
                  !ENDIF
               ENDIF

               !! RZK-warning simple orthogonalization for block-diagonal ALMOs
               !CALL orthogonalize_mos(ket=matrix_t_out(ispin),&
               !        overlap=almo_scf_env%matrix_sigma_blk(ispin),&
               !        metric=almo_scf_env%matrix_s(1),&
               !        retain_locality=.TRUE.,&
               !        only_normalize=.TRUE.,&
               !        eps_filter=almo_scf_env%eps_filter,&
               !        order_lanczos=almo_scf_env%order_lanczos,&
               !        eps_lanczos=almo_scf_env%eps_lanczos,&
               !        max_iter_lanczos=almo_scf_env%max_iter_lanczos)

               CALL cp_dbcsr_filter(matrix_t_out(ispin), &
                                    eps=almo_scf_env%eps_filter)

               ! compute the density matrix
               CALL almo_scf_t_to_p( &
                  t=matrix_t_out(ispin), &
                  p=almo_scf_env%matrix_p(ispin), &
                  eps_filter=almo_scf_env%eps_filter, &
                  orthog_orbs=.FALSE., &
                  s=almo_scf_env%matrix_s(1), &
                  sigma=almo_scf_env%matrix_sigma(ispin), &
                  sigma_inv=almo_scf_env%matrix_sigma_inv(ispin), &
                  use_guess=use_guess)
               CALL cp_dbcsr_scale(almo_scf_env%matrix_p(ispin), &
                                   spin_factor)

               !! RZK-warning to debug lets see eigenvalues of the ALMO overlap
               !CALL cp_dbcsr_copy(almo_scf_env%matrix_sigma_blk(ispin),&
               !               almo_scf_env%matrix_sigma(ispin),&
               !               keep_sparsity=.TRUE.)
               !CALL cp_dbcsr_init(u_sigma)
               !CALL cp_dbcsr_create(u_sigma,template=almo_scf_env%matrix_sigma(ispin),&
               !        matrix_type=dbcsr_type_no_symmetry)
               !CALL cp_dbcsr_get_info(almo_scf_env%matrix_sigma_blk(ispin), nfullrows_total=occ1 )
               !ALLOCATE(evals(occ1))
               !CALL cp_dbcsr_syevd(almo_scf_env%matrix_sigma_blk(ispin),u_sigma,evals,&
               !        almo_scf_env%para_env,almo_scf_env%blacs_env)
               !WRITE(*,*) 'SIGMA_BLCK: ', evals
               !DEALLOCATE(evals)
               !CALL cp_dbcsr_create(u_sigma,template=almo_scf_env%matrix_sigma(ispin),&
               !        matrix_type=dbcsr_type_no_symmetry)
               !CALL cp_dbcsr_get_info(almo_scf_env%matrix_sigma(ispin), nfullrows_total=occ1 )
               !ALLOCATE(evals(occ1))
               !CALL cp_dbcsr_get_diag(almo_scf_env%matrix_sigma(ispin),evals)
               !WRITE(*,*) 'TRACE: ', SUM(evals)
               !CALL cp_dbcsr_syevd(almo_scf_env%matrix_sigma(ispin),u_sigma,evals,&
               !        almo_scf_env%para_env,almo_scf_env%blacs_env)
               !WRITE(*,*) 'SIGMA_FULL: ', evals
               !WRITE(*,*) 'SUMEI: ', SUM(evals)
               !DEALLOCATE(evals)
               !CALL cp_dbcsr_release(u_sigma)

               ! update the KS matrix and energy if necessary
               !NEW IF (.NOT. (perturbation_only .AND. (.NOT. just_started))) THEN
               IF ( .NOT.perturbation_only ) THEN
                  !!! RZK-warning the KS matrix must be updated outside the spin loop
                  !!! Now the code works only for restricted orbitals
                  !IF (unit_nr>0) THEN
                  !   WRITE(unit_nr,*) "....updating KS matrix...."
                  !ENDIF
                  CALL almo_scf_dm_to_ks(qs_env, almo_scf_env, energy_new)
                  CALL get_qs_env(qs_env, matrix_ks=matrix_ks)
                  CALL matrix_qs_to_almo(matrix_ks(ispin)%matrix, &
                                         almo_scf_env%matrix_ks(ispin), &
                                         almo_scf_env, .FALSE.)
                  CALL cp_dbcsr_filter(almo_scf_env%matrix_ks(ispin), &
                                       almo_scf_env%eps_filter)
               ENDIF

               !NEW ! save the initial state
               !NEW IF (just_started) THEN
               !NEW    CALL cp_dbcsr_copy(matrix_p_0, almo_scf_env%matrix_p(ispin))
               !NEW    CALL cp_dbcsr_desymmetrize(almo_scf_env%matrix_sigma(ispin), &
               !NEW                               matrix_sigma_0)
               !NEW    CALL cp_dbcsr_desymmetrize(almo_scf_env%matrix_sigma_inv(ispin), &
               !NEW                               matrix_sigma_inv_0)
               !NEW ENDIF

               !NEW IF (my_special_case .EQ. xalmo_case_normal .AND. prec_type .EQ. 4) THEN
               !NEW    ! construct domain-projector
               !NEW    IF (assume_t0_q0x .AND. just_started) THEN
               !NEW       CALL construct_domain_r_down( &
               !NEW          matrix_t=matrix_t_0, &
               !NEW          matrix_sigma_inv=matrix_sigma_inv_0, &
               !NEW          matrix_s=almo_scf_env%matrix_s(1), &
               !NEW          subm_r_down=domain_r_down(:), &
               !NEW          dpattern=quench_t(ispin), &
               !NEW          map=almo_scf_env%domain_map(ispin), &
               !NEW          node_of_domain=almo_scf_env%cpu_of_domain, &
               !NEW          filter_eps=almo_scf_env%eps_filter)
               !NEW       !CALL construct_domain_r_down(&
               !NEW       !        matrix_t=matrix_t_out(ispin),&
               !NEW       !        matrix_sigma_inv=almo_scf_env%matrix_sigma_inv(ispin),&
               !NEW       !        matrix_s=almo_scf_env%matrix_s(1),&
               !NEW       !        subm_r_down=domain_r_down(:),&
               !NEW       !        dpattern=quench_t(ispin),&
               !NEW       !        map=almo_scf_env%domain_map(ispin),&
               !NEW       !        node_of_domain=almo_scf_env%cpu_of_domain,&
               !NEW       !        filter_eps=almo_scf_env%eps_filter)
               !NEW    ENDIF ! assume_t0_q0x
               !NEW ENDIF ! debug condition

               IF (perturbation_only) THEN
                  ! calculate objective function Tr(F_0 R)
                  CALL cp_dbcsr_multiply("N", "N", 1.0_dp, &
                                         almo_scf_env%matrix_ks(ispin), &
                                         matrix_t_out(ispin), &
                                         0.0_dp, m_tmp_no_1, &
                                         filter_eps=almo_scf_env%eps_filter)
                  CALL cp_dbcsr_multiply("N", "N", 1.0_dp, &
                                         m_tmp_no_1, &
                                         almo_scf_env%matrix_sigma_inv(ispin), &
                                         0.0_dp, FTsiginv, &
                                         filter_eps=almo_scf_env%eps_filter)
                  CALL cp_dbcsr_trace(matrix_t_out(ispin), &
                                      FTsiginv, energy_new, "T", "N")
                  energy_new = energy_new*spin_factor
               ENDIF

               ! save the previous gradient to compute beta
               ! do it only if the previous grad was computed
               ! for .NOT.line_search
               IF (line_search_iteration .EQ. 0 .AND. iteration .NE. 0) &
                  CALL cp_dbcsr_copy(prev_grad, grad)

               ! compute the energy gradient if necessary
               skip_grad = (iteration .GT. 0 .AND. &
                            fixed_line_search_niter .NE. 0 .AND. &
                            line_search_iteration .NE. fixed_line_search_niter)

               IF (.NOT. skip_grad) THEN

                  !IF (unit_nr>0) THEN
                  !   WRITE(unit_nr,*) "....computing gradient...."
                  !ENDIF
                  CALL compute_gradient(&
                          m_grad_out=m_tmp_no_2,&
                          !m_grad_out=grad(ispin),&
                          m_ks=almo_scf_env%matrix_ks(ispin),&
                          m_s=almo_scf_env%matrix_s(1),&
                          m_t=matrix_t_out(ispin),&
                          m_t0=matrix_t_0,&
                          m_siginv=almo_scf_env%matrix_sigma_inv(ispin),&
                          m_quench_t=quench_t(ispin),&
                          m_FTsiginv=FTsiginv(ispin),&
                          m_siginvTFTsiginv=siginvTFTsiginv(ispin),&
                          m_ST=ST(ispin),&
                          m_STsiginv0=STsiginv_0,&
                          m_theta=m_theta,&
                          !m_sig_sqrti_ii=m_sig_sqrti_ii(ispin),&
                          domain_s_inv=almo_scf_env%domain_s_inv(:,ispin),&
                          domain_r_down=domain_r_down(:,ispin),&
                          cpu_of_domain=almo_scf_env%cpu_of_domain,&
                          domain_map=almo_scf_env%domain_map(ispin),&
                          assume_t0_q0x=assume_t0_q0x,&
                          optimize_theta=optimize_theta,&
                          perturbation_only=perturbation_only,&
                          normalize_orbitals=.FALSE.,&
                          penalty_occ_vol=.FALSE.,&
                          penalty_occ_vol_prefactor=0.0_dp,&
                          !normalize_orbitals=normalize_orbitals,&
                          !penalty_occ_vol=penalty_occ_vol,&
                          !penalty_occ_vol_prefactor=penalty_occ_vol_g_prefactor,&
                          envelope_amplitude=almo_scf_env%envelope_amplitude,&
                          eps_filter=almo_scf_env%eps_filter,&
                          spin_factor=spin_factor,&
                          special_case=my_special_case)


                  !! do d_E/d_theta first
                  !IF (.NOT. perturbation_only) THEN
                  !   CALL cp_dbcsr_multiply("N", "N", 1.0_dp, &
                  !                          almo_scf_env%matrix_ks(ispin), &
                  !                          matrix_t_out(ispin), &
                  !                          0.0_dp, m_tmp_no_1, &
                  !                          filter_eps=almo_scf_env%eps_filter)
                  !   CALL cp_dbcsr_multiply("N", "N", 1.0_dp, &
                  !                          m_tmp_no_1, &
                  !                          almo_scf_env%matrix_sigma_inv(ispin), &
                  !                          0.0_dp, FTsiginv, &
                  !                          filter_eps=almo_scf_env%eps_filter)
                  !   !CALL cp_dbcsr_multiply("N","N",1.0_dp,&
                  !   !        m_tmp_no_1,&
                  !   !        almo_scf_env%matrix_sigma_inv(ispin),&
                  !   !        0.0_dp,m_tmp_no_2,&
                  !   !        retain_sparsity=.TRUE.)
                  !ENDIF
                  !CALL cp_dbcsr_copy(m_tmp_no_2, quench_t(ispin))
                  !CALL cp_dbcsr_copy(m_tmp_no_2, &
                  !                   FTsiginv, keep_sparsity=.TRUE.)
                  !CALL cp_dbcsr_multiply("T", "N", 1.0_dp, &
                  !                       matrix_t_out(ispin), &
                  !                       FTsiginv, &
                  !                       0.0_dp, m_tmp_oo_1, &
                  !                       filter_eps=almo_scf_env%eps_filter)
                  !!CALL cp_dbcsr_multiply("N","N",1.0_dp,&
                  !!        m_tmp_oo_1,&
                  !!        almo_scf_env%matrix_sigma_inv(ispin),&
                  !!        0.0_dp,m_tmp_oo_2,&
                  !!        filter_eps=almo_scf_env%eps_filter)
                  !CALL cp_dbcsr_multiply("N", "N", 1.0_dp, &
                  !                       almo_scf_env%matrix_sigma_inv(ispin), &
                  !                       m_tmp_oo_1, &
                  !                       0.0_dp, siginvTFTsiginv, &
                  !                       filter_eps=almo_scf_env%eps_filter)
                  !CALL cp_dbcsr_multiply("N", "N", 1.0_dp, &
                  !                       almo_scf_env%matrix_s(1), &
                  !                       matrix_t_out(ispin), &
                  !                       0.0_dp, ST, &
                  !                       filter_eps=almo_scf_env%eps_filter)
                  !!NEW ! save S*T_0*siginv_0
                  !!NEW IF (assume_t0_q0x .AND. just_started &
                  !!NEW     .AND. special_case .EQ. xalmo_case_fully_deloc) THEN
                  !!NEW    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, &
                  !!NEW                           ST, &
                  !!NEW                           almo_scf_env%matrix_sigma_inv(ispin), &
                  !!NEW                           0.0_dp, STsiginv_0, &
                  !!NEW                           filter_eps=almo_scf_env%eps_filter)
                  !!NEW ENDIF
                  !CALL cp_dbcsr_multiply("N", "N", -1.0_dp, &
                  !                       ST, &
                  !                       siginvTFTsiginv, &
                  !                       1.0_dp, m_tmp_no_2, &
                  !                       retain_sparsity=.TRUE.)
                  !CALL cp_dbcsr_scale(m_tmp_no_2, &
                  !                    2.0_dp*spin_factor)
                  !CALL cp_dbcsr_filter(m_tmp_no_2, &
                  !                     eps=almo_scf_env%eps_filter)

                  !!NEW IF (perturbation_only .AND. just_started) THEN
                  !!NEW    ! save the first gradient
                  !!NEW    ! it is equal to F_vo and necessary to compute
                  !!NEW    ! the correction to the energy
                  !!NEW    CALL cp_dbcsr_copy(fvo_0, m_tmp_no_2)
                  !!NEW    CALL cp_dbcsr_scale(fvo_0, 0.5_dp/spin_factor)
                  !!NEW ENDIF

                  !! a short print-out
                  !CALL cp_dbcsr_norm(m_tmp_no_2, &
                  !                   dbcsr_norm_maxabsnorm, norm_scalar=t_norm)
                  !!IF (unit_nr>0) THEN
                  !!   WRITE(unit_nr,*) "Maximum norm of dE/dT: ", t_norm
                  !!ENDIF

                  !! project out the occupied space from the gradient
                  !IF (assume_t0_q0x) THEN
                  !   IF (my_special_case .EQ. xalmo_case_fully_deloc) THEN
                  !      CALL cp_dbcsr_copy(grad, m_tmp_no_2)
                  !      CALL cp_dbcsr_multiply("T", "N", 1.0_dp, &
                  !                             matrix_t_0, &
                  !                             grad, &
                  !                             0.0_dp, m_tmp_oo_1, &
                  !                             filter_eps=almo_scf_env%eps_filter)
                  !      CALL cp_dbcsr_multiply("N", "N", -1.0_dp, &
                  !                             STsiginv_0, &
                  !                             m_tmp_oo_1, &
                  !                             1.0_dp, grad, &
                  !                             filter_eps=almo_scf_env%eps_filter)
                  !   ELSE IF (my_special_case .EQ. xalmo_case_block_diag) THEN
                  !      ! should not be here - cannot project the zero-order space from itself
                  !      CPABORT("")
                  !   ELSE
                  !      ! no special case: normal xALMOs
                  !      CALL apply_domain_operators( &
                  !         matrix_in=m_tmp_no_2, &
                  !         matrix_out=grad, &
                  !         operator2=domain_r_down(:), &
                  !         operator1=almo_scf_env%domain_s_inv(:, ispin), &
                  !         dpattern=quench_t(ispin), &
                  !         map=almo_scf_env%domain_map(ispin), &
                  !         node_of_domain=almo_scf_env%cpu_of_domain, &
                  !         my_action=1, &
                  !         filter_eps=almo_scf_env%eps_filter, &
                  !         !matrix_trimmer=,&
                  !         use_trimmer=.FALSE.)
                  !   ENDIF ! my_special_case
                  !   CALL cp_dbcsr_copy(m_tmp_no_2, grad)
                  !ENDIF

                  !! transform d_E/d_T to d_E/d_theta
                  !IF (optimize_theta) THEN
                  !   CALL cp_dbcsr_copy(m_tmp_no_1, m_theta)
                  !   CALL cp_dbcsr_function_of_elements(m_tmp_no_1, &
                  !                                      !func=dbcsr_func_cos,&
                  !                                      func=dbcsr_func_dtanh, &
                  !                                      a0=0.0_dp, &
                  !                                      a1=1.0_dp/almo_scf_env%envelope_amplitude)
                  !   CALL cp_dbcsr_scale(m_tmp_no_1, &
                  !                       almo_scf_env%envelope_amplitude)
                  !   CALL cp_dbcsr_set(m_tmp_no_3, 0.0_dp)
                  !   CALL cp_dbcsr_filter(m_tmp_no_3, &
                  !                        eps=almo_scf_env%eps_filter)
                  !   CALL cp_dbcsr_hadamard_product(m_tmp_no_2, &
                  !                                  m_tmp_no_1, &
                  !                                  m_tmp_no_3, &
                  !                                  b_assume_value=1.0_dp)
                  !   CALL cp_dbcsr_hadamard_product(m_tmp_no_3, &
                  !                                  quench_t(ispin), &
                  !                                  grad)
                  !ELSE ! simply copy
                  !   CALL cp_dbcsr_hadamard_product(m_tmp_no_2, &
                  !                                  quench_t(ispin), &
                  !                                  grad)
                  !ENDIF
                  !CALL cp_dbcsr_filter(grad, eps=almo_scf_env%eps_filter)

               ENDIF ! skip_grad

               ! check convergence and other exit criteria
               grad_norm_frob = cp_dbcsr_frobenius_norm(grad)
               CALL cp_dbcsr_norm(grad, dbcsr_norm_maxabsnorm, &
                                  norm_scalar=grad_norm)
               converged = (grad_norm .LT. optimizer%eps_error)
               IF (converged .OR. (iteration .GE. max_iter)) THEN
                  prepare_to_exit = .TRUE.
               ENDIF
               IF (grad_norm .LT. almo_scf_env%eps_prev_guess) THEN
                  use_guess = .TRUE.
               ENDIF

               IF (md_in_theta_space) THEN

                  IF (.NOT. first_md_iteration) THEN
                     CALL cp_dbcsr_copy(prev_step, step)
                  ENDIF
                  CALL cp_dbcsr_copy(step, grad)
                  CALL cp_dbcsr_scale(step, -1.0_dp)

                  ! update velocities v(i) = v(i-1) + 0.5*dT*(a(i-1) + a(i))
                  IF (.NOT. first_md_iteration) THEN
                     CALL cp_dbcsr_add(velocity, &
                                       step, 1.0_dp, 0.5_dp*time_step)
                     CALL cp_dbcsr_add(velocity, &
                                       prev_step, 1.0_dp, 0.5_dp*time_step)
                  ENDIF
                  kin_energy = cp_dbcsr_frobenius_norm(velocity)
                  kin_energy = 0.5_dp*kin_energy*kin_energy

                  ! update positions theta(i) = theta(i-1) + dT*v(i-1) + 0.5*dT*dT*a(i-1)
                  CALL cp_dbcsr_add(m_theta, &
                                    velocity, 1.0_dp, time_step)
                  CALL cp_dbcsr_add(m_theta, &
                                    step, 1.0_dp, 0.5_dp*time_step*time_step)

                  iter_type = "MD"

                  t2 = m_walltime()
                  IF (unit_nr > 0) THEN
                     WRITE (unit_nr, '(T2,A,A2,I5,F16.7,F17.9,F17.9,F17.9,E12.3,F10.3)') &
                        "ALMO SCF ", iter_type, iteration, time_step*iteration, &
                        energy_new, kin_energy, energy_new+kin_energy, grad_norm, &
                        t2-t1
                  ENDIF
                  t1 = m_walltime()

                  IF (first_md_iteration) THEN
                     first_md_iteration = .FALSE.
                  ENDIF

               ELSE ! optimizization (not MD)

                  IF (.NOT. prepare_to_exit) THEN

                     ! check the gradient along the step direction
                     IF (iteration .NE. 0) THEN

                        IF (fixed_line_search_niter .EQ. 0) THEN

                           CALL cp_dbcsr_trace(grad, step, line_search_error, &
                                               "T", "N")
                           ! normalize the result
                           !IF (unit_nr>0) THEN
                           !   WRITE(unit_nr,*) "Angle between step/grad: ", line_search_error
                           !ENDIF
                           CALL cp_dbcsr_trace(grad, grad, denom, "T", "N")
                           !IF (unit_nr>0) THEN
                           !   WRITE(unit_nr,*) "Frobenius norm of grad:  ", SQRT(denom)
                           !ENDIF
                           line_search_error = line_search_error/SQRT(denom)
                           CALL cp_dbcsr_trace(step, step, denom, "T", "N")
                           !IF (unit_nr>0) THEN
                           !   WRITE(unit_nr,*) "Frobenius norm of step:  ", SQRT(denom)
                           !ENDIF
                           line_search_error = line_search_error/SQRT(denom)
                           IF (ABS(line_search_error) .GT. optimizer%lin_search_eps_error) THEN
                              line_search = .TRUE.
                              line_search_iteration = line_search_iteration+1
                           ELSE
                              line_search = .FALSE.
                              line_search_iteration = 0
                              IF (grad_norm .LT. eps_skip_gradients) THEN
                                 fixed_line_search_niter = ABS(almo_scf_env%integer04)
                              ENDIF
                           ENDIF

                        ELSE ! decision for fixed_line_search_niter

                           IF (.NOT. line_search) THEN
                              line_search = .TRUE.
                              line_search_iteration = line_search_iteration+1
                           ELSE
                              IF (line_search_iteration .EQ. fixed_line_search_niter) THEN
                                 line_search = .FALSE.
                                 line_search_iteration = 0
                                 line_search_iteration = line_search_iteration+1
                              ENDIF
                           ENDIF

                        ENDIF ! fixed_line_search_niter fork
                     ENDIF

                     IF (line_search) THEN
                        energy_diff = 0.0_dp
                     ELSE
                        energy_diff = energy_new-energy_old
                        energy_old = energy_new
                     ENDIF

                     ! update the step direction
                     IF (.NOT. line_search) THEN

                        !IF (unit_nr>0) THEN
                        !   WRITE(unit_nr,*) "....updating step direction...."
                        !ENDIF

                        cg_iteration = cg_iteration+1

                        IF ((just_started .AND. perturbation_only) .OR. &
                            (iteration .EQ. 0 .AND. (.NOT. perturbation_only))) THEN

                           ! compute the preconditioner
                           !IF (unit_nr>0) THEN
                           !   WRITE(unit_nr,*) "....computing preconditioner...."
                           !ENDIF

                           ! calculate (1-R)F(1-R) and S-SRS
                           ! RZK-warning take advantage: some elements will be removed by the quencher
                           ! RZK-warning S operations can be performed outside the spin loop to save time
                           ! IT IS REQUIRED THAT PRECONDITIONER DOES NOT BREAK THE LOCALITY!!!!
                           ! RZK-warning: further optimization is ABSOLUTELY NECESSARY

                           ! First S-SRS
                           !CALL cp_dbcsr_multiply("N","N",1.0_dp,&
                           !        almo_scf_env%matrix_s(1),&
                           !        matrix_t_out(ispin),&
                           !        0.0_dp,m_tmp_no_1,&
                           !        filter_eps=almo_scf_env%eps_filter)
                           CALL cp_dbcsr_multiply("N", "N", 1.0_dp, &
                                                  ST, &
                                                  almo_scf_env%matrix_sigma_inv(ispin), &
                                                  0.0_dp, m_tmp_no_3, &
                                                  filter_eps=almo_scf_env%eps_filter)
                           CALL cp_dbcsr_create(m_tmp_nn_1, &
                                                template=almo_scf_env%matrix_s(1), &
                                                matrix_type=dbcsr_type_no_symmetry)
                           CALL cp_dbcsr_desymmetrize(almo_scf_env%matrix_s(1), &
                                                      m_tmp_nn_1)
                           IF (my_special_case .EQ. xalmo_case_fully_deloc) THEN
                              ! use S instead of S-SRS
                           ELSE
                              CALL cp_dbcsr_multiply("N", "T", -1.0_dp, &
                                                     ST, &
                                                     m_tmp_no_3, &
                                                     1.0_dp, m_tmp_nn_1, &
                                                     filter_eps=almo_scf_env%eps_filter)
                           ENDIF

                           ! Second (1-R)F(1-R)
                           !CALL cp_dbcsr_multiply("N","N",1.0_dp,&
                           !        almo_scf_env%matrix_ks(ispin),&
                           !        matrix_t_out(ispin),&
                           !        0.0_dp,m_tmp_no_1,&
                           !        filter_eps=almo_scf_env%eps_filter)
                           ! re-create matrix because desymmetrize is buggy -
                           ! it will create multiple copies of blocks
                           CALL cp_dbcsr_create(prec_vv, &
                                                template=almo_scf_env%matrix_ks(ispin), &
                                                matrix_type=dbcsr_type_no_symmetry)
                           CALL cp_dbcsr_desymmetrize(almo_scf_env%matrix_ks(ispin), &
                                                      prec_vv)
                           CALL cp_dbcsr_multiply("N", "T", -1.0_dp, &
                                                  FTsiginv, &
                                                  ST, &
                                                  1.0_dp, prec_vv, &
                                                  filter_eps=almo_scf_env%eps_filter)
                           CALL cp_dbcsr_multiply("N", "T", -1.0_dp, &
                                                  ST, &
                                                  FTsiginv, &
                                                  1.0_dp, prec_vv, &
                                                  filter_eps=almo_scf_env%eps_filter)
                           !CALL cp_dbcsr_multiply("T","N",1.0_dp,&
                           !        matrix_t_out(ispin),&
                           !        m_tmp_no_1,&
                           !        0.0_dp,m_tmp_oo_1,&
                           !        filter_eps=almo_scf_env%eps_filter)
                           CALL cp_dbcsr_multiply("N", "N", 1.0_dp, &
                                                  ST, &
                                                  siginvTFTsiginv, &
                                                  0.0_dp, m_tmp_no_3, &
                                                  filter_eps=almo_scf_env%eps_filter)
                           CALL cp_dbcsr_multiply("N", "T", 1.0_dp, &
                                                  m_tmp_no_3, &
                                                  ST, &
                                                  1.0_dp, prec_vv, &
                                                  filter_eps=almo_scf_env%eps_filter)
                           CALL cp_dbcsr_add(prec_vv, m_tmp_nn_1, &
                                             1.0_dp-prec_sf_mixing_s, &
                                             prec_sf_mixing_s)
                           CALL cp_dbcsr_scale(prec_vv, 2.0_dp*spin_factor)
                           CALL cp_dbcsr_copy(m_tmp_nn_1, prec_vv)

                           ! invert using various algorithms
                           IF (my_special_case .EQ. xalmo_case_block_diag) THEN ! non-overlapping diagonal blocks

                              !precond_domain_projector=0
                              !CALL construct_domain_preconditioner(&
                              !   matrix_main=m_tmp_nn_1,&
                              !   dpattern=quench_t(ispin),&
                              !   map=almo_scf_env%domain_map(ispin),&
                              !   node_of_domain=almo_scf_env%cpu_of_domain,&
                              !   preconditioner=almo_scf_env%domain_preconditioner(:,ispin),&
                              !   use_trimmer=.FALSE.,&
                              !   my_action=precond_domain_projector)
                              CALL pseudo_invert_diagonal_blk(matrix_in=m_tmp_nn_1, &
                                                              matrix_out=prec_vv, &
                                                              nocc=almo_scf_env%nocc_of_domain(:, ispin))

                           ELSE IF (my_special_case .EQ. xalmo_case_fully_deloc) THEN ! the entire system is a block

                              ! invert using cholesky (works with S matrix, will not work with S-SRS matrix)
                              CALL cp_dbcsr_cholesky_decompose(prec_vv, &
                                                               para_env=almo_scf_env%para_env, &
                                                               blacs_env=almo_scf_env%blacs_env)
                              CALL cp_dbcsr_cholesky_invert(prec_vv, &
                                                            para_env=almo_scf_env%para_env, &
                                                            blacs_env=almo_scf_env%blacs_env, &
                                                            upper_to_full=.TRUE.)
                              CALL cp_dbcsr_filter(prec_vv, &
                                                   eps=almo_scf_env%eps_filter)
                           ELSE
                              !!! use a sophisticated domain preconditioner
                              IF (assume_t0_q0x) THEN
                                 precond_domain_projector = -1
                              ELSE
                                 precond_domain_projector = 0
                              ENDIF
                              ! for other experimental preconditioner types the inversion is
                              ! done together with applying the preconditioner
                              IF (prec_type .EQ. 4) THEN

                                 CALL construct_domain_preconditioner( &
                                    matrix_main=m_tmp_nn_1, &
                                    subm_s_inv=almo_scf_env%domain_s_inv(:, ispin), &
                                    subm_r_down=domain_r_down(:), &
                                    matrix_trimmer=quench_t(ispin), &
                                    dpattern=quench_t(ispin), &
                                    map=almo_scf_env%domain_map(ispin), &
                                    node_of_domain=almo_scf_env%cpu_of_domain, &
                                    preconditioner=almo_scf_env%domain_preconditioner(:, ispin), &
                                    use_trimmer=.FALSE., &
                                    my_action=precond_domain_projector)
                              ENDIF ! prec type
                           ENDIF

                           ! invert using cholesky (works with S matrix, will not work with S-SRS matrix)
                           !!!CALL cp_dbcsr_cholesky_decompose(prec_vv,&
                           !!!        para_env=almo_scf_env%para_env,&
                           !!!        blacs_env=almo_scf_env%blacs_env)
                           !!!CALL cp_dbcsr_cholesky_invert(prec_vv,&
                           !!!        para_env=almo_scf_env%para_env,&
                           !!!        blacs_env=almo_scf_env%blacs_env,&
                           !!!        upper_to_full=.TRUE.)
                           !!!CALL cp_dbcsr_filter(prec_vv,&
                           !!!        eps=almo_scf_env%eps_filter)
                           !!!

                           ! re-create the matrix because desymmetrize is buggy -
                           ! it will create multiple copies of blocks
                           !!!DESYM!CALL cp_dbcsr_create(prec_vv,&
                           !!!DESYM!        template=almo_scf_env%matrix_s(1),&
                           !!!DESYM!        matrix_type=dbcsr_type_no_symmetry)
                           !!!DESYM!CALL cp_dbcsr_desymmetrize(almo_scf_env%matrix_s(1),&
                           !!!DESYM!        prec_vv)
                           !CALL cp_dbcsr_multiply("N","N",1.0_dp,&
                           !        almo_scf_env%matrix_s(1),&
                           !        matrix_t_out(ispin),&
                           !        0.0_dp,m_tmp_no_1,&
                           !        filter_eps=almo_scf_env%eps_filter)
                           !CALL cp_dbcsr_multiply("N","N",1.0_dp,&
                           !        m_tmp_no_1,&
                           !        almo_scf_env%matrix_sigma_inv(ispin),&
                           !        0.0_dp,m_tmp_no_3,&
                           !        filter_eps=almo_scf_env%eps_filter)
                           !CALL cp_dbcsr_multiply("N","T",-1.0_dp,&
                           !        m_tmp_no_3,&
                           !        m_tmp_no_1,&
                           !        1.0_dp,prec_vv,&
                           !        filter_eps=almo_scf_env%eps_filter)
                           !CALL cp_dbcsr_add_on_diag(prec_vv,&
                           !        prec_sf_mixing_s)

                           !CALL cp_dbcsr_create(prec_oo,&
                           !        template=almo_scf_env%matrix_sigma(ispin),&
                           !        matrix_type=dbcsr_type_no_symmetry)
                           !CALL cp_dbcsr_desymmetrize(almo_scf_env%matrix_sigma(ispin),&
                           !        matrix_type=dbcsr_type_no_symmetry)
                           !CALL cp_dbcsr_desymmetrize(almo_scf_env%matrix_sigma(ispin),&
                           !        prec_oo)
                           !CALL cp_dbcsr_filter(prec_oo,&
                           !        eps=almo_scf_env%eps_filter)

                           !! invert using cholesky
                           !CALL cp_dbcsr_create(prec_oo_inv,&
                           !        template=prec_oo,&
                           !        matrix_type=dbcsr_type_no_symmetry)
                           !CALL cp_dbcsr_desymmetrize(prec_oo,&
                           !        prec_oo_inv)
                           !CALL cp_dbcsr_cholesky_decompose(prec_oo_inv,&
                           !        para_env=almo_scf_env%para_env,&
                           !        blacs_env=almo_scf_env%blacs_env)
                           !CALL cp_dbcsr_cholesky_invert(prec_oo_inv,&
                           !        para_env=almo_scf_env%para_env,&
                           !        blacs_env=almo_scf_env%blacs_env,&
                           !        upper_to_full=.TRUE.)

                        ENDIF

                        ! save the previous step
                        CALL cp_dbcsr_copy(prev_step, step)

                        ! compute the new step (apply preconditioner if available)
                        IF (use_preconditioner) THEN

                           !IF (unit_nr>0) THEN
                           !   WRITE(unit_nr,*) "....applying preconditioner...."
                           !ENDIF

                           SELECT CASE (prec_type)
                           CASE (1)
                              ! expensive Newton-Raphson step (the Hessian is still approximate)
                              ! RZK-warning THIS PREC HAS NOT BEEN IMPLEMENTED FOR THETA
                              IF (ncores .GT. 1) THEN
                                 CPABORT("serial code only")
                              ENDIF
                              CALL newton_grad_to_step( &
                                 matrix_grad=m_tmp_no_2, &
                                 matrix_step=m_tmp_no_1, &
                                 matrix_s=almo_scf_env%matrix_s(1), &
                                 matrix_ks=almo_scf_env%matrix_ks(ispin), &
                                 matrix_t=matrix_t_out(ispin), &
                                 matrix_sigma_inv=almo_scf_env%matrix_sigma_inv(ispin), &
                                 !matrix_ks=matrix_ks_0,&
                                 !matrix_t=matrix_t_0,&
                                 !matrix_sigma_inv=matrix_sigma_inv_0,&
                                 quench_t=quench_t(ispin), &
                                 spin_factor=spin_factor, &
                                 eps_filter=almo_scf_env%eps_filter)

                           CASE (3)

                              ! RZK-warning THIS PREC HAS NOT BEEN IMPLEMENTED FOR THETA

                              ! inversion
                              CALL cp_dbcsr_get_info(m_tmp_nn_1, nfullrows_total=dim0)
                              ALLOCATE (evals(dim0))
                              CALL cp_dbcsr_syevd(m_tmp_nn_1, prec_vv, evals, &
                                                  almo_scf_env%para_env, almo_scf_env%blacs_env)
                              ! invert eigenvalues and use eigenvectors to compute the Hessian inverse
                              ! take special care of zero eigenvalues
                              zero_neg_eiv = 0
                              CALL cp_dbcsr_get_info(almo_scf_env%matrix_sigma(ispin), nfullrows_total=occ1)
                              DO jj = 1, dim0
                                 IF (jj .LE. occ1) THEN
                                    evals(jj) = evals(jj)*0.0_dp
                                    zero_neg_eiv = zero_neg_eiv+1
                                 ELSE
                                    evals(jj) = 1.0_dp/evals(jj)
                                 ENDIF
                              ENDDO
                              IF (unit_nr > 0) THEN
                                 WRITE (*, *) 'ZERO OR NEGATIVE EIGENVALUES: ', zero_neg_eiv
                              ENDIF
                              CALL cp_dbcsr_init(inv_eiv)
                              CALL cp_dbcsr_create(inv_eiv, &
                                                   template=m_tmp_nn_1, &
                                                   matrix_type=dbcsr_type_no_symmetry)
                              CALL cp_dbcsr_add_on_diag(inv_eiv, 1.0_dp)
                              CALL cp_dbcsr_set_diag(inv_eiv, evals)
                              CALL cp_dbcsr_multiply("N", "N", 1.0_dp, &
                                                     prec_vv, &
                                                     inv_eiv, &
                                                     0.0_dp, m_tmp_nn_1, &
                                                     filter_eps=almo_scf_env%eps_filter)
                              CALL cp_dbcsr_multiply("N", "T", 1.0_dp, &
                                                     m_tmp_nn_1, &
                                                     prec_vv, &
                                                     0.0_dp, inv_eiv, &
                                                     filter_eps=almo_scf_env%eps_filter)
                              CALL cp_dbcsr_copy(prec_vv, inv_eiv)
                              CALL cp_dbcsr_release(inv_eiv)
                              DEALLOCATE (evals)

                              !!CALL cp_dbcsr_copy(step,&
                              !!        quench_t(ispin))
                              !!
                              !!CALL cp_dbcsr_multiply("N","N",1.0_dp,&
                              !!        m_tmp_no_2,&
                              !!        !grad,& - this choice is worse
                              !!        prec_oo,&
                              !!        0.0_dp,step,&
                              !!        !retain_sparsity=.TRUE.,&
                              !!        filter_eps=almo_scf_env%eps_filter)
                              !!
                              CALL cp_dbcsr_copy(m_tmp_no_1, &
                                                 quench_t(ispin))
                              !!CALL cp_dbcsr_hadamard_product(&
                              !!        quench_t(ispin),&
                              !!        step,&
                              !!        m_tmp_no_1)
                              !!
                              CALL cp_dbcsr_multiply("N", "N", -1.0_dp, &
                                                     prec_vv, &
                                                     m_tmp_no_2, &
                                                     0.0_dp, m_tmp_no_1, &
                                                     retain_sparsity=.TRUE.)

                           CASE (4)

                              IF (my_special_case .EQ. xalmo_case_block_diag .OR. &
                                  my_special_case .EQ. xalmo_case_fully_deloc) THEN

                                 CALL cp_dbcsr_multiply("N", "N", -1.0_dp, &
                                                        prec_vv, &
                                                        grad, &
                                                        0.0_dp, step, &
                                                        filter_eps=almo_scf_env%eps_filter)

                              ELSE

                                 !!! RZK-warning Currently for non-theta only
                                 IF (optimize_theta) THEN
                                    CPABORT("theta is NYI")
                                 ENDIF

                                 CALL apply_domain_operators( &
                                    matrix_in=grad, &
                                    matrix_out=step, &
                                    operator1=almo_scf_env%domain_preconditioner(:, ispin), &
                                    !operator2=,&
                                    dpattern=quench_t(ispin), &
                                    map=almo_scf_env%domain_map(ispin), &
                                    node_of_domain=almo_scf_env%cpu_of_domain, &
                                    my_action=0, &
                                    filter_eps=almo_scf_env%eps_filter)
                                 !matrix_trimmer=,&
                                 !use_trimmer=.FALSE.,&
                                 CALL cp_dbcsr_scale(step, -1.0_dp)

                                 CALL cp_dbcsr_copy(m_tmp_no_3, &
                                                    quench_t(ispin))
                                 CALL cp_dbcsr_function_of_elements(m_tmp_no_3, &
                                                                    func=dbcsr_func_inverse, &
                                                                    a0=0.0_dp, &
                                                                    a1=1.0_dp)
                                 CALL cp_dbcsr_copy(m_tmp_no_2, step)
                                 CALL cp_dbcsr_hadamard_product( &
                                    m_tmp_no_2, &
                                    m_tmp_no_3, &
                                    step)
                                 CALL cp_dbcsr_copy(m_tmp_no_3, quench_t(ispin))

                                 !CALL cp_dbcsr_create(m_tmp_oo_1,&
                                 !        template=almo_scf_env%matrix_sigma_blk(ispin),&
                                 !        matrix_type=dbcsr_type_no_symmetry)
                                 !CALL cp_dbcsr_desymmetrize(almo_scf_env%matrix_sigma_blk(ispin),m_tmp_oo_1)
                                 !CALL get_overlap(bra=matrix_t_out(ispin),&
                                 !        ket=step,&
                                 !        overlap=m_tmp_oo_1,&
                                 !        metric=almo_scf_env%matrix_s(1),&
                                 !        retain_overlap_sparsity=.TRUE.,&
                                 !        eps_filter=almo_scf_env%eps_filter)
                                 !CALL cp_dbcsr_norm(m_tmp_oo_1,&
                                 !        dbcsr_norm_maxabsnorm, norm_scalar=t_norm)
                                 !IF (unit_nr>0) THEN
                                 !   WRITE(unit_nr,*) "Step block-orthogonality error: ", t_norm
                                 !ENDIF

                              ENDIF ! special case

                           END SELECT ! preconditioner type fork

                        ELSE

                           !!!! NO PRECONDITIONER
                           CALL cp_dbcsr_copy(step, grad)
                           CALL cp_dbcsr_scale(step, -1.0_dp)

                        ENDIF

                        ! check whether we need to reset conjugate directions
                        IF (iteration .EQ. 0) THEN
                           reset_conjugator = .TRUE.
                        ENDIF

                        ! compute the conjugation coefficient - beta
                        IF (.NOT. reset_conjugator) THEN

                           SELECT CASE (optimizer%conjugator)
                           CASE (cg_hestenes_stiefel)
                              CALL cp_dbcsr_copy(m_tmp_no_1, grad)
                              CALL cp_dbcsr_add(m_tmp_no_1, prev_grad, &
                                                1.0_dp, -1.0_dp)
                              CALL cp_dbcsr_trace(m_tmp_no_1, step, numer, &
                                                  "T", "N")
                              CALL cp_dbcsr_trace(m_tmp_no_1, prev_step, denom, &
                                                  "T", "N")
                              beta = -1.0_dp*numer/denom
                           CASE (cg_fletcher_reeves)
                              CALL cp_dbcsr_trace(grad, step, numer, "T", "N")
                              CALL cp_dbcsr_trace(prev_grad, prev_minus_prec_grad, denom, "T", "N")
                              beta = numer/denom
                           CASE (cg_polak_ribiere)
                              CALL cp_dbcsr_trace(prev_grad, prev_minus_prec_grad, denom, "T", "N")
                              CALL cp_dbcsr_copy(m_tmp_no_1, grad)
                              CALL cp_dbcsr_add(m_tmp_no_1, prev_grad, 1.0_dp, -1.0_dp)
                              CALL cp_dbcsr_trace(m_tmp_no_1, step, numer, "T", "N")
                              beta = numer/denom
                           CASE (cg_fletcher)
                              CALL cp_dbcsr_trace(grad, step, numer, "T", "N")
                              CALL cp_dbcsr_trace(prev_grad, prev_step, denom, "T", "N")
                              beta = numer/denom
                           CASE (cg_liu_storey)
                              CALL cp_dbcsr_trace(prev_grad, prev_step, denom, "T", "N")
                              CALL cp_dbcsr_copy(m_tmp_no_1, grad)
                              CALL cp_dbcsr_add(m_tmp_no_1, prev_grad, 1.0_dp, -1.0_dp)
                              CALL cp_dbcsr_trace(m_tmp_no_1, step, numer, "T", "N")
                              beta = numer/denom
                           CASE (cg_dai_yuan)
                              CALL cp_dbcsr_trace(grad, step, numer, "T", "N")
                              CALL cp_dbcsr_copy(m_tmp_no_1, grad)
                              CALL cp_dbcsr_add(m_tmp_no_1, prev_grad, 1.0_dp, -1.0_dp)
                              CALL cp_dbcsr_trace(m_tmp_no_1, prev_step, denom, "T", "N")
                              beta = -1.0_dp*numer/denom
                           CASE (cg_hager_zhang)
                              CALL cp_dbcsr_copy(m_tmp_no_1, grad)
                              CALL cp_dbcsr_add(m_tmp_no_1, prev_grad, 1.0_dp, -1.0_dp)
                              CALL cp_dbcsr_trace(m_tmp_no_1, prev_step, denom, "T", "N")
                              CALL cp_dbcsr_trace(m_tmp_no_1, prev_minus_prec_grad, numer, "T", "N")
                              kappa = -2.0_dp*numer/denom
                              CALL cp_dbcsr_trace(m_tmp_no_1, step, numer, "T", "N")
                              tau = -1.0_dp*numer/denom
                              CALL cp_dbcsr_trace(prev_step, grad, numer, "T", "N")
                              beta = tau-kappa*numer/denom
                           CASE (cg_zero)
                              beta = 0.0_dp
                           CASE DEFAULT
                              CPABORT("illegal conjugator")
                           END SELECT

                           IF (beta .LT. 0.0_dp) THEN
                              IF (unit_nr > 0) THEN
                                 WRITE (unit_nr, *) "Beta is negative: ", beta
                              ENDIF
                              reset_conjugator = .TRUE.
                           ENDIF

                        ENDIF

                        IF (reset_conjugator) THEN

                           beta = 0.0_dp
                           IF (unit_nr > 0 .AND. (.NOT. just_started)) THEN
                              WRITE (unit_nr, *) "(Re)-setting conjugator to zero"
                           ENDIF
                           reset_conjugator = .FALSE.

                        ENDIF

                        ! save the preconditioned gradient (useful for beta)
                        CALL cp_dbcsr_copy(prev_minus_prec_grad, step)

                        !IF (unit_nr>0) THEN
                        !   WRITE(unit_nr,*) "....final beta....", beta
                        !ENDIF

                        ! conjugate the step direction
                        CALL cp_dbcsr_add(step, prev_step, 1.0_dp, beta)

                     ENDIF ! update the step direction

                     ! estimate the step size
                     IF (.NOT. line_search) THEN
                        e0 = energy_new
                        CALL cp_dbcsr_trace(grad, step, g0, "T", "N")
                        ! we just changed the direction and
                        ! we have only E and grad from the current step
                        ! it is not enouhg to compute step_size - just guess it
                        IF (iteration .EQ. 0) THEN
                           step_size = optimizer%lin_search_step_size_guess
                        ELSE
                           IF (next_step_size_guess .LE. 0.0_dp) THEN
                              step_size = optimizer%lin_search_step_size_guess
                           ELSE
                              ! take the last value
                              step_size = next_step_size_guess*1.05_dp
                           ENDIF
                        ENDIF
                        next_step_size_guess = step_size
                     ELSE
                        IF (fixed_line_search_niter .EQ. 0) THEN
                           e1 = energy_new
                           CALL cp_dbcsr_trace(grad, step, g1, "T", "N")
                           ! we have accumulated some points along this direction
                           ! use only the most recent g0 (quadratic approximation)
                           appr_sec_der = (g1-g0)/step_size
                           !IF (unit_nr>0) THEN
                           !   WRITE(unit_nr,'(A2,7F12.5)') &
                           !           "EG",e0,e1,g0,g1,appr_sec_der,step_size,-g1/appr_sec_der
                           !ENDIF
                           step_size = -g1/appr_sec_der
                           e0 = e1
                           g0 = g1
                        ELSE
                           ! use e0, g0 and e1 to compute g1 and make a step
                           ! if the next iteration is also line_search
                           ! use e1 and the calculated g1 as e0 and g0
                           e1 = energy_new
                           appr_sec_der = 2.0*((e1-e0)/step_size-g0)/step_size
                           g1 = appr_sec_der*step_size+g0
                           IF (unit_nr > 0) THEN
                              WRITE (unit_nr, '(A2,7F12.5)') &
                                 "EG", e0, e1, g0, g1, appr_sec_der, step_size, -g1/appr_sec_der
                           ENDIF
                           !appr_sec_der=(g1-g0)/step_size
                           step_size = -g1/appr_sec_der
                           e0 = e1
                           g0 = g1
                        ENDIF
                        next_step_size_guess = next_step_size_guess+step_size
                     ENDIF

                     ! update theta
                     CALL cp_dbcsr_add(m_theta, step, 1.0_dp, step_size)

                  ENDIF ! not.prepare_to_exit

                  IF (line_search) THEN
                     iter_type = "LS"
                  ELSE
                     iter_type = "CG"
                  ENDIF

                  t2 = m_walltime()
                  IF (unit_nr > 0) THEN
                     iter_type = TRIM("ALMO SCF "//iter_type)
                     WRITE (unit_nr, '(T2,A13,I6,F23.10,E14.5,F14.9,F9.2)') &
                        iter_type, iteration, &
                        energy_new, energy_diff, grad_norm, &
                        t2-t1
                     !WRITE(unit_nr,'(T2,A11,I6,F20.12,E12.3,E12.3,E12.3,F12.5,F10.3)') &
                     !        "ALMO SCF ",iter_type,iteration,&
                     !        energy_new,energy_diff,grad_norm,line_search_error,&
                     !        step_size,t2-t1
                  ENDIF

                  IF (my_special_case .EQ. xalmo_case_block_diag) THEN
                     almo_scf_env%almo_scf_energy = energy_new
                  ENDIF

                  t1 = m_walltime()

               ENDIF ! MD in theta space

               iteration = iteration+1
               IF (prepare_to_exit) EXIT

            ENDDO ! inner SCF loop

            IF (converged .OR. (outer_iteration .GE. outer_max_iter)) THEN
               outer_prepare_to_exit = .TRUE.
            ENDIF

            outer_iteration = outer_iteration+1
            IF (outer_prepare_to_exit) EXIT

         ENDDO ! outer SCF loop

         ! post SCF-loop calculations
         IF (converged) THEN

            ! obtain MO coefficients from final theta
            ! RZK-warning: if decide to uncomment make sure that Theta->T
            ! procedure is consistent with the procedure in the loop
            !IF (optimize_theta) THEN
            !   CALL cp_dbcsr_copy(m_tmp_no_1,m_theta)
            !   CALL cp_dbcsr_function_of_elements(m_tmp_no_1,&
            !           func=dbcsr_func_tanh,&
            !           a0=0.0_dp,&
            !           a1=1.0_dp/almo_scf_env%envelope_amplitude)
            !   CALL cp_dbcsr_hadamard_product(m_tmp_no_1,&
            !           quench_t(ispin),&
            !           matrix_t_out(ispin))
            !   CALL cp_dbcsr_scale(matrix_t_out(ispin),&
            !           almo_scf_env%envelope_amplitude)
            !ELSE
            !   CALL cp_dbcsr_copy(m_tmp_no_1,m_theta)
            !   CALL cp_dbcsr_hadamard_product(m_tmp_no_1,&
            !           quench_t(ispin),&
            !           matrix_t_out(ispin))
            !ENDIF
            !IF (perturbation_only) THEN
            !   CALL cp_dbcsr_add(matrix_t_out(ispin),&
            !           matrix_t_0,1.0_dp,1.0_dp)
            !ENDIF
            !CALL cp_dbcsr_filter(matrix_t_out(ispin),&
            !        eps=almo_scf_env%eps_filter)
            !CALL cp_dbcsr_norm(matrix_t_out(ispin),&
            !        dbcsr_norm_maxabsnorm, norm_scalar=grad_norm)
            !IF (unit_nr>0) THEN
            !   WRITE(unit_nr,*) "Maximum norm of the ALMOs: ", grad_norm
            !ENDIF

            !!!!! experiment: bump final amplitudes and get T
            !!!!IF (.NOT.optimize_theta.AND.perturbation_only) THEN
            !!!!   CALL cp_dbcsr_copy(m_tmp_no_1,m_theta)
            !!!!   CALL cp_dbcsr_hadamard_product(m_tmp_no_1,&
            !!!!           quench_t_saved,&
            !!!!           matrix_t_out(ispin))
            !!!!   IF (use_projector) THEN
            !!!!      CALL almo_scf_domain_operations(&
            !!!!              matrix_in=matrix_t_out(ispin),&
            !!!!              matrix_out=m_tmp_no_1,&
            !!!!              quench_t=quench_t(ispin),&
            !!!!              my_action=2,&
            !!!!              matrix_s=almo_scf_env%matrix_s(1),&
            !!!!              matrix_sigma_inv=matrix_sigma_inv_0,&
            !!!!              matrix_t=matrix_t_0)
            !!!!
            !!!!      CALL cp_dbcsr_copy(matrix_t_out(ispin),&
            !!!!              m_tmp_no_1)
            !!!!   ENDIF
            !!!!   IF (perturbation_only) THEN
            !!!!      CALL cp_dbcsr_add(matrix_t_out(ispin),&
            !!!!              matrix_t_0,1.0_dp,1.0_dp)
            !!!!
            !!!!   ENDIF
            !!!!ENDIF
            !!!!
            !!!!IF (.NOT.optimize_theta) THEN
            !!!!   CALL cp_dbcsr_copy(quench_t(ispin),quench_t_saved)
            !!!!ENDIF

            IF (perturbation_only) THEN

               ! use step matrix to store excitation amplitudes
               CALL cp_dbcsr_copy(step, matrix_t_0)
               CALL cp_dbcsr_add(step, matrix_t_out(ispin), &
                                 -1.0_dp, 1.0_dp)

               CALL cp_dbcsr_scale(fvo_0, spin_factor)
               CALL cp_dbcsr_trace(step, fvo_0, energy_new, "T", "N")

               !CALL cp_dbcsr_add(matrix_t_0,matrix_t_out(ispin),&
               !        -1.0_dp,1.0_dp)
               !CALL cp_dbcsr_trace(matrix_t_0,&
               !        fvo_0,energy_new,"T","N")

               ! print out the energy lowering
               IF (unit_nr > 0) THEN
                  WRITE (unit_nr, *)
                  WRITE (unit_nr, '(T2,A35,F25.10)') "ENERGY OF BLOCK-DIAGONAL ALMOs:", &
                     almo_scf_env%almo_scf_energy
                  WRITE (unit_nr, '(T2,A35,F25.10)') "ENERGY LOWERING:", &
                     energy_new
                  WRITE (unit_nr, '(T2,A35,F25.10)') "CORRECTED ENERGY:", &
                     almo_scf_env%almo_scf_energy+energy_new
                  WRITE (unit_nr, *)
               ENDIF
               CALL almo_scf_update_ks_energy(qs_env, &
                                              almo_scf_env%almo_scf_energy+energy_new)

               IF (almo_scf_env%almo_analysis%do_analysis) THEN

                  ! energy decomposition analysis (EDA)
                  IF (unit_nr > 0) THEN
                     WRITE (unit_nr, '(T2,A)') "DECOMPOSITION OF THE DELOCALIZATION ENERGY"
                  ENDIF

                  CALL cp_dbcsr_hadamard_product(step, &
                                                 fvo_0, m_tmp_no_1)
                  !CALL cp_dbcsr_hadamard_product(matrix_t_0,&
                  !        fvo_0,m_tmp_no_1)
                  CALL cp_dbcsr_filter(m_tmp_no_1, almo_scf_env%eps_filter)

                  ! open the output file, print and close
                  CALL get_qs_env(qs_env, input=input)
                  almo_print_section => section_vals_get_subs_vals(input, "DFT%ALMO_SCF%ANALYSIS%PRINT")
                  eda_unit = cp_print_key_unit_nr(logger, almo_print_section, &
                                                  "ALMO_EDA_CT", extension=".dat", local=.TRUE.)
                  CALL cp_dbcsr_print_block_sum(m_tmp_no_1, eda_unit)
                  CALL cp_print_key_finished_output(eda_unit, logger, almo_print_section, &
                                                    "ALMO_EDA_CT", local=.TRUE.)

                  ! charge transfer analysis (CTA)
                  IF (unit_nr > 0) THEN
                     WRITE (unit_nr, '(T2,A)') "DECOMPOSITION OF CHARGE TRANSFER TERMS"
                  ENDIF
                  !NEW ! first, compute [QR'R]_mu^i = [(S-SRS).X.siginv']_mu^i
                  !NEW ! a. tmp1 = S.T0
                  !NEW CALL cp_dbcsr_multiply("N", "N", 1.0_dp, &
                  !NEW                        almo_scf_env%matrix_s(1), &
                  !NEW                        matrix_t_0, &
                  !NEW                        0.0_dp, m_tmp_no_1, &
                  !NEW                        filter_eps=almo_scf_env%eps_filter)
                  !NEW ! b. tmp3 = tmp1.siginv0 = S.T0.siginv0
                  !NEW CALL cp_dbcsr_multiply("N", "N", 1.0_dp, &
                  !NEW                        m_tmp_no_1, &
                  !NEW                        matrix_sigma_inv_0, &
                  !NEW                        0.0_dp, m_tmp_no_3, &
                  !NEW                        filter_eps=almo_scf_env%eps_filter)
                  ! c. tmp1 = S.X
                  CALL cp_dbcsr_multiply("N", "N", 1.0_dp, &
                                         almo_scf_env%matrix_s(1), &
                                         step, &
                                         0.0_dp, m_tmp_no_1, &
                                         filter_eps=almo_scf_env%eps_filter)
                  ! d. tmp2 = tr(T0).tmp1 = tr(T0).S.X
                  CALL cp_dbcsr_multiply("T", "N", 1.0_dp, &
                                         matrix_t_0, &
                                         m_tmp_no_1, &
                                         0.0_dp, m_tmp_oo_1, &
                                         filter_eps=almo_scf_env%eps_filter)
                  ! e. tmp1 = tmp1 - tmp3.tmp2 = S.X - S.T0.siginv0.tr(T0).S.X
                  !         = (1-S.R0).S.X
                  CALL cp_dbcsr_multiply("N", "N", -1.0_dp, &
                  !NEW                        m_tmp_no_3, &
                                         STsiginv_0, &
                                         m_tmp_oo_1, &
                                         1.0_dp, m_tmp_no_1, &
                                         filter_eps=almo_scf_env%eps_filter)
                  ! f. tmp2 = tmp1.siginv
                  CALL cp_dbcsr_multiply("N", "N", 1.0_dp, &
                                         m_tmp_no_1, &
                                         almo_scf_env%matrix_sigma_inv(ispin), &
                                         0.0_dp, m_tmp_no_2, &
                                         filter_eps=almo_scf_env%eps_filter)
                  ! second, compute traces of blocks [RR'Q]^x_y * [X]^y_x
                  CALL cp_dbcsr_hadamard_product(step, &
                                                 m_tmp_no_2, m_tmp_no_1)
                  CALL cp_dbcsr_scale(m_tmp_no_1, spin_factor)
                  CALL cp_dbcsr_filter(m_tmp_no_1, almo_scf_env%eps_filter)

                  eda_unit = cp_print_key_unit_nr(logger, almo_print_section, &
                                                  "ALMO_CTA", extension=".dat", local=.TRUE.)
                  CALL cp_dbcsr_print_block_sum(m_tmp_no_1, eda_unit)
                  CALL cp_print_key_finished_output(eda_unit, logger, almo_print_section, &
                                                    "ALMO_CTA", local=.TRUE.)

               ENDIF ! do ALMO EDA/CTA

            ELSE

               CALL almo_scf_update_ks_energy(qs_env, &
                                              energy_new)

            ENDIF ! if perturbation only

         ENDIF ! if converged

         IF (md_in_theta_space) THEN
            CALL cp_dbcsr_release(velocity)
         ENDIF
         CALL cp_dbcsr_release(m_theta)
         CALL cp_dbcsr_release(prec_vv)
         CALL cp_dbcsr_release(prec_oo)
         CALL cp_dbcsr_release(prec_oo_inv)
         CALL cp_dbcsr_release(m_tmp_no_1)
         CALL cp_dbcsr_release(fvo_0)
         CALL cp_dbcsr_release(STsiginv_0)
         CALL cp_dbcsr_release(m_tmp_no_2)
         CALL cp_dbcsr_release(m_tmp_no_3)
         CALL cp_dbcsr_release(m_tmp_oo_1)
         CALL cp_dbcsr_release(ST)
         CALL cp_dbcsr_release(FTsiginv)
         CALL cp_dbcsr_release(siginvTFTsiginv)
         CALL cp_dbcsr_release(m_tmp_nn_1)
         CALL cp_dbcsr_release(prev_grad)
         CALL cp_dbcsr_release(prev_step)
         CALL cp_dbcsr_release(grad)
         CALL cp_dbcsr_release(step)
         CALL cp_dbcsr_release(prev_minus_prec_grad)
         !CALL cp_dbcsr_release(matrix_p_0)
         CALL cp_dbcsr_release(matrix_t_0)
         !CALL cp_dbcsr_release(matrix_sigma_0)
         CALL cp_dbcsr_release(matrix_sigma_inv_0)

         IF (.NOT. converged) THEN
            CPABORT("Optimization not converged! ")
         ENDIF

         DEALLOCATE (domain_r_down)

      ENDDO ! ispin

      CALL timestop(handle)

   END SUBROUTINE almo_scf_xalmo_pcg

! **************************************************************************************************
!> \brief Split the matrix of virtual orbitals into two:
!>        retained orbs and discarded
!> \param almo_scf_env ...
!> \par History
!>       2011.09 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE split_v_blk(almo_scf_env)

      TYPE(almo_scf_env_type)                            :: almo_scf_env

      CHARACTER(len=*), PARAMETER :: routineN = 'split_v_blk', routineP = moduleN//':'//routineN

      INTEGER                                            :: discarded_v, handle, iblock_col, &
                                                            iblock_col_size, iblock_row, &
                                                            iblock_row_size, ispin, retained_v
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: data_p, p_new_block
      TYPE(cp_dbcsr_iterator)                            :: iter

      CALL timeset(routineN, handle)

      DO ispin = 1, almo_scf_env%nspins

         CALL cp_dbcsr_work_create(almo_scf_env%matrix_v_blk(ispin), &
                                   work_mutable=.TRUE.)
         CALL cp_dbcsr_work_create(almo_scf_env%matrix_v_disc_blk(ispin), &
                                   work_mutable=.TRUE.)

         CALL cp_dbcsr_iterator_start(iter, almo_scf_env%matrix_v_full_blk(ispin))

         DO WHILE (cp_dbcsr_iterator_blocks_left(iter))

            CALL cp_dbcsr_iterator_next_block(iter, iblock_row, iblock_col, data_p, &
                                              row_size=iblock_row_size, col_size=iblock_col_size)

            IF (iblock_row .NE. iblock_col) THEN
               CPABORT("off-diagonal block found")
            ENDIF

            retained_v = almo_scf_env%nvirt_of_domain(iblock_col, ispin)
            discarded_v = almo_scf_env%nvirt_disc_of_domain(iblock_col, ispin)
            CPASSERT(retained_v .GT. 0)
            CPASSERT(discarded_v .GT. 0)

            NULLIFY (p_new_block)
            CALL cp_dbcsr_reserve_block2d(almo_scf_env%matrix_v_disc_blk(ispin), &
                                          iblock_row, iblock_col, p_new_block)
            CPASSERT(ASSOCIATED(p_new_block))
            CPASSERT(retained_v+discarded_v .EQ. iblock_col_size)
            p_new_block(:, :) = data_p(:, (retained_v+1):iblock_col_size)

            NULLIFY (p_new_block)
            CALL cp_dbcsr_reserve_block2d(almo_scf_env%matrix_v_blk(ispin), &
                                          iblock_row, iblock_col, p_new_block)
            CPASSERT(ASSOCIATED(p_new_block))
            p_new_block(:, :) = data_p(:, 1:retained_v)

         ENDDO ! iterator
         CALL cp_dbcsr_iterator_stop(iter)

         CALL cp_dbcsr_finalize(almo_scf_env%matrix_v_blk(ispin))
         CALL cp_dbcsr_finalize(almo_scf_env%matrix_v_disc_blk(ispin))

      ENDDO ! ispin

      CALL timestop(handle)

   END SUBROUTINE split_v_blk

! **************************************************************************************************
!> \brief various methods for calculating the Harris-Foulkes correction
!> \param almo_scf_env ...
!> \par History
!>       2011.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE harris_foulkes_correction(almo_scf_env)

      TYPE(almo_scf_env_type)                            :: almo_scf_env

      CHARACTER(len=*), PARAMETER :: routineN = 'harris_foulkes_correction', &
         routineP = moduleN//':'//routineN
      INTEGER, PARAMETER                                 :: cayley_transform = 1, dm_ls_step = 2

      INTEGER :: algorithm_id, handle, handle1, handle2, handle3, handle4, handle5, handle6, &
         handle7, handle8, ispin, iteration, n, nmins, nspin, opt_k_max_iter, &
         outer_opt_k_iteration, outer_opt_k_max_iter, unit_nr
      INTEGER, DIMENSION(1)                              :: fake, nelectron_spin_real
      LOGICAL :: converged, line_search, md_in_k_space, outer_opt_k_prepare_to_exit, &
         prepare_to_exit, reset_conjugator, reset_step_size, use_cubic_approximation, &
         use_quadratic_approximation
      REAL(KIND=dp) :: aa, bb, beta, conjugacy_error, conjugacy_error_threshold, &
         delta_obj_function, denom, energy_correction_final, frob_matrix, frob_matrix_base, fun0, &
         fun1, gfun0, gfun1, grad_norm, grad_norm_frob, kappa, kin_energy, line_search_error, &
         line_search_error_threshold, num_threshold, numer, obj_function, quadratic_approx_error, &
         quadratic_approx_error_threshold, safety_multiplier, spin_factor, step_size, &
         step_size_quadratic_approx, step_size_quadratic_approx2, t1, t1a, t1cholesky, t2, t2a, &
         t2cholesky, tau, time_step, x_opt_eps_adaptive, x_opt_eps_adaptive_factor
      REAL(KIND=dp), DIMENSION(1)                        :: local_mu
      REAL(KIND=dp), DIMENSION(2)                        :: energy_correction
      REAL(KIND=dp), DIMENSION(3)                        :: minima
      TYPE(cp_dbcsr_type) :: grad, k_vd_index_down, k_vr_index_down, matrix_k_central, &
         matrix_tmp1, matrix_tmp2, prec, prev_grad, prev_minus_prec_grad, prev_step, &
         sigma_oo_curr, sigma_oo_curr_inv, sigma_vv_sqrt, sigma_vv_sqrt_guess, sigma_vv_sqrt_inv, &
         sigma_vv_sqrt_inv_guess, step, t_curr, tmp1_n_vr, tmp2_n_o, tmp3_vd_vr, tmp4_o_vr, &
         tmp_k_blk, vd_fixed, vd_index_sqrt, vd_index_sqrt_inv, velocity, vr_fixed, vr_index_sqrt, &
         vr_index_sqrt_inv
      TYPE(cp_dbcsr_type), ALLOCATABLE, DIMENSION(:)     :: matrix_p_almo_scf_converged
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(ct_step_env_type)                             :: ct_step_env

!TYPE(cp_dbcsr_iterator)                  :: iter
!TYPE(cp_dbcsr_type) :: tmp11,tmp22,tmp33
!REAL(kind=dp)                            :: k_var1, k_var2
!TYPE(cp_dbcsr_type)                      :: fake_step
!
!TYPE(cp_dbcsr_type)                      :: sigma_dr, sigma_dr2, sigma_rr, sigma_rr2
!TYPE(cp_dbcsr_type)                      :: fake_a,fake_b,fake_k0
!INTEGER                                  :: retained_v,discarded_v,i_row,j_col
!TYPE(cp_dbcsr_type) :: matrix_rst0, matrix_rst1, matrix_rst2, ss_vv
!REAL(KIND=dp)       :: filter_memorize, init_filter, occ_vv
!INTEGER             :: ppp

      CALL timeset(routineN, handle)

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%mepos == logger%para_env%source) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      nspin = almo_scf_env%nspins
      energy_correction_final = 0.0_dp
      IF (nspin .EQ. 1) THEN
         spin_factor = 2.0_dp
      ELSE
         spin_factor = 1.0_dp
      ENDIF

      IF (almo_scf_env%deloc_use_occ_orbs) THEN
         algorithm_id = cayley_transform
      ELSE
         algorithm_id = dm_ls_step
      ENDIF

      t1 = m_walltime()

      SELECT CASE (algorithm_id)
      CASE (cayley_transform)

         ! rescale density matrix by spin factor
         ! so the orbitals and density are consistent with each other
         IF (almo_scf_env%nspins == 1) THEN
            CALL cp_dbcsr_scale(almo_scf_env%matrix_p(1), 1.0_dp/spin_factor)
         ENDIF

         ! transform matrix_t not matrix_t_blk (we might need ALMOs later)
         DO ispin = 1, nspin

            CALL cp_dbcsr_copy(almo_scf_env%matrix_t(ispin), &
                               almo_scf_env%matrix_t_blk(ispin))

            ! obtain orthogonalization matrices for ALMOs
            ! RZK-warning - remove this sqrt(sigma) and inv(sqrt(sigma))
            ! ideally ALMO scf should use sigma and sigma_inv in
            ! the tensor_up_down representation

            IF (unit_nr > 0) THEN
               WRITE (unit_nr, *) "sqrt and inv(sqrt) of MO overlap matrix"
            ENDIF
            CALL cp_dbcsr_init(almo_scf_env%matrix_sigma_sqrt(ispin))
            CALL cp_dbcsr_init(almo_scf_env%matrix_sigma_sqrt_inv(ispin))
            CALL cp_dbcsr_create(almo_scf_env%matrix_sigma_sqrt(ispin), &
                                 template=almo_scf_env%matrix_sigma(ispin), &
                                 matrix_type=dbcsr_type_no_symmetry)
            CALL cp_dbcsr_create(almo_scf_env%matrix_sigma_sqrt_inv(ispin), &
                                 template=almo_scf_env%matrix_sigma(ispin), &
                                 matrix_type=dbcsr_type_no_symmetry)

            CALL matrix_sqrt_Newton_Schulz(almo_scf_env%matrix_sigma_sqrt(ispin), &
                                           almo_scf_env%matrix_sigma_sqrt_inv(ispin), &
                                           almo_scf_env%matrix_sigma(ispin), &
                                           threshold=almo_scf_env%eps_filter, &
                                           order=almo_scf_env%order_lanczos, &
                                           eps_lanczos=almo_scf_env%eps_lanczos, &
                                           max_iter_lanczos=almo_scf_env%max_iter_lanczos)

            IF (safe_mode) THEN
               CALL cp_dbcsr_init(matrix_tmp1)
               CALL cp_dbcsr_create(matrix_tmp1, template=almo_scf_env%matrix_sigma(ispin), &
                                    matrix_type=dbcsr_type_no_symmetry)
               CALL cp_dbcsr_init(matrix_tmp2)
               CALL cp_dbcsr_create(matrix_tmp2, template=almo_scf_env%matrix_sigma(ispin), &
                                    matrix_type=dbcsr_type_no_symmetry)

               CALL cp_dbcsr_multiply("N", "N", 1.0_dp, almo_scf_env%matrix_sigma_sqrt_inv(ispin), &
                                      almo_scf_env%matrix_sigma(ispin), &
                                      0.0_dp, matrix_tmp1, filter_eps=almo_scf_env%eps_filter)
               CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1, &
                                      almo_scf_env%matrix_sigma_sqrt_inv(ispin), &
                                      0.0_dp, matrix_tmp2, filter_eps=almo_scf_env%eps_filter)

               frob_matrix_base = cp_dbcsr_frobenius_norm(matrix_tmp2)
               CALL cp_dbcsr_add_on_diag(matrix_tmp2, -1.0_dp)
               frob_matrix = cp_dbcsr_frobenius_norm(matrix_tmp2)
               IF (unit_nr > 0) THEN
                  WRITE (unit_nr, *) "Error for (inv(sqrt(SIG))*SIG*inv(sqrt(SIG))-I)", frob_matrix/frob_matrix_base
               ENDIF

               CALL cp_dbcsr_release(matrix_tmp1)
               CALL cp_dbcsr_release(matrix_tmp2)
            ENDIF
         ENDDO

         IF (almo_scf_env%almo_update_algorithm .EQ. almo_scf_diag) THEN

            DO ispin = 1, nspin

               t1a = m_walltime()

               line_search_error_threshold = almo_scf_env%real01
               conjugacy_error_threshold = almo_scf_env%real02
               quadratic_approx_error_threshold = almo_scf_env%real03
               x_opt_eps_adaptive_factor = almo_scf_env%real04

               !! the outer loop for k optimization
               outer_opt_k_max_iter = almo_scf_env%opt_k_outer_max_iter
               outer_opt_k_prepare_to_exit = .FALSE.
               outer_opt_k_iteration = 0
               grad_norm = 0.0_dp
               grad_norm_frob = 0.0_dp
               CALL cp_dbcsr_set(almo_scf_env%matrix_x(ispin), 0.0_dp)
               IF (almo_scf_env%deloc_truncate_virt .EQ. virt_full) outer_opt_k_max_iter = 0

               DO

                  ! obtain proper retained virtuals (1-R)|ALMO_vr>
                  CALL apply_projector(psi_in=almo_scf_env%matrix_v_blk(ispin), &
                                       psi_out=almo_scf_env%matrix_v(ispin), &
                                       psi_projector=almo_scf_env%matrix_t_blk(ispin), &
                                       metric=almo_scf_env%matrix_s(1), &
                                       project_out=.TRUE., &
                                       psi_projector_orthogonal=.FALSE., &
                                       proj_in_template=almo_scf_env%matrix_ov(ispin), &
                                       eps_filter=almo_scf_env%eps_filter, &
                                       sig_inv_projector=almo_scf_env%matrix_sigma_inv(ispin))
                  !sig_inv_template=almo_scf_env%matrix_sigma_inv(ispin),&

                  ! save initial retained virtuals
                  CALL cp_dbcsr_init(vr_fixed)
                  CALL cp_dbcsr_create(vr_fixed, &
                                       template=almo_scf_env%matrix_v(ispin))
                  CALL cp_dbcsr_copy(vr_fixed, almo_scf_env%matrix_v(ispin))

                  ! init matrices common for optimized and non-optimized virts
                  CALL cp_dbcsr_init(sigma_vv_sqrt)
                  CALL cp_dbcsr_init(sigma_vv_sqrt_inv)
                  CALL cp_dbcsr_create(sigma_vv_sqrt, &
                                       template=almo_scf_env%matrix_sigma_vv(ispin), &
                                       matrix_type=dbcsr_type_no_symmetry)
                  CALL cp_dbcsr_create(sigma_vv_sqrt_inv, &
                                       template=almo_scf_env%matrix_sigma_vv(ispin), &
                                       matrix_type=dbcsr_type_no_symmetry)
                  CALL cp_dbcsr_init(sigma_vv_sqrt_inv_guess)
                  CALL cp_dbcsr_init(sigma_vv_sqrt_guess)
                  CALL cp_dbcsr_create(sigma_vv_sqrt_inv_guess, &
                                       template=almo_scf_env%matrix_sigma_vv(ispin), &
                                       matrix_type=dbcsr_type_no_symmetry)
                  CALL cp_dbcsr_create(sigma_vv_sqrt_guess, &
                                       template=almo_scf_env%matrix_sigma_vv(ispin), &
                                       matrix_type=dbcsr_type_no_symmetry)
                  CALL cp_dbcsr_set(sigma_vv_sqrt_guess, 0.0_dp)
                  CALL cp_dbcsr_add_on_diag(sigma_vv_sqrt_guess, 1.0_dp)
                  CALL cp_dbcsr_filter(sigma_vv_sqrt_guess, almo_scf_env%eps_filter)
                  CALL cp_dbcsr_set(sigma_vv_sqrt_inv_guess, 0.0_dp)
                  CALL cp_dbcsr_add_on_diag(sigma_vv_sqrt_inv_guess, 1.0_dp)
                  CALL cp_dbcsr_filter(sigma_vv_sqrt_inv_guess, almo_scf_env%eps_filter)

                  ! do things required to optimize virtuals
                  IF (almo_scf_env%deloc_truncate_virt .NE. virt_full) THEN

                     ! project retained virtuals out of discarded block-by-block
                     ! (1-Q^VR_ALMO)|ALMO_vd>
                     ! this is probably not necessary, do it just to be safe
                     !CALL apply_projector(psi_in=almo_scf_env%matrix_v_disc_blk(ispin),&
                     !        psi_out=almo_scf_env%matrix_v_disc(ispin),&
                     !        psi_projector=almo_scf_env%matrix_v_blk(ispin),&
                     !        metric=almo_scf_env%matrix_s_blk(1),&
                     !        project_out=.TRUE.,&
                     !        psi_projector_orthogonal=.FALSE.,&
                     !        proj_in_template=almo_scf_env%matrix_k_tr(ispin),&
                     !        eps_filter=almo_scf_env%eps_filter,&
                     !        sig_inv_template=almo_scf_env%matrix_sigma_vv(ispin))
                     !CALL cp_dbcsr_copy(almo_scf_env%matrix_v_disc_blk(ispin),&
                     !        almo_scf_env%matrix_v_disc(ispin))

                     ! construct discarded virtuals (1-R)|ALMO_vd>
                     CALL apply_projector(psi_in=almo_scf_env%matrix_v_disc_blk(ispin), &
                                          psi_out=almo_scf_env%matrix_v_disc(ispin), &
                                          psi_projector=almo_scf_env%matrix_t_blk(ispin), &
                                          metric=almo_scf_env%matrix_s(1), &
                                          project_out=.TRUE., &
                                          psi_projector_orthogonal=.FALSE., &
                                          proj_in_template=almo_scf_env%matrix_ov_disc(ispin), &
                                          eps_filter=almo_scf_env%eps_filter, &
                                          sig_inv_projector=almo_scf_env%matrix_sigma_inv(ispin))
                     !sig_inv_template=almo_scf_env%matrix_sigma_inv(ispin),&

                     ! save initial discarded
                     CALL cp_dbcsr_init(vd_fixed)
                     CALL cp_dbcsr_create(vd_fixed, &
                                          template=almo_scf_env%matrix_v_disc(ispin))
                     CALL cp_dbcsr_copy(vd_fixed, almo_scf_env%matrix_v_disc(ispin))

                     !! create the down metric in the retained k-subspace
                     CALL cp_dbcsr_init(k_vr_index_down)
                     CALL cp_dbcsr_create(k_vr_index_down, &
                                          template=almo_scf_env%matrix_sigma_vv_blk(ispin), &
                                          matrix_type=dbcsr_type_no_symmetry)
                     !CALL cp_dbcsr_copy(k_vr_index_down,&
                     !        almo_scf_env%matrix_sigma_vv_blk(ispin))

                     !CALL get_overlap(bra=almo_scf_env%matrix_v_blk(ispin),&
                     !        ket=almo_scf_env%matrix_v_blk(ispin),&
                     !        overlap=k_vr_index_down,&
                     !        metric=almo_scf_env%matrix_s_blk(1),&
                     !        retain_overlap_sparsity=.FALSE.,&
                     !        eps_filter=almo_scf_env%eps_filter)

                     !! create the up metric in the discarded k-subspace
                     CALL cp_dbcsr_init(k_vd_index_down)
                     CALL cp_dbcsr_create(k_vd_index_down, &
                                          template=almo_scf_env%matrix_vv_disc_blk(ispin), &
                                          matrix_type=dbcsr_type_no_symmetry)
                     !CALL cp_dbcsr_init(k_vd_index_up)
                     !CALL cp_dbcsr_create(k_vd_index_up,&
                     !        template=almo_scf_env%matrix_vv_disc_blk(ispin),&
                     !        matrix_type=dbcsr_type_no_symmetry)
                     !CALL cp_dbcsr_copy(k_vd_index_down,&
                     !        almo_scf_env%matrix_vv_disc_blk(ispin))

                     !CALL get_overlap(bra=almo_scf_env%matrix_v_disc_blk(ispin),&
                     !        ket=almo_scf_env%matrix_v_disc_blk(ispin),&
                     !        overlap=k_vd_index_down,&
                     !        metric=almo_scf_env%matrix_s_blk(1),&
                     !        retain_overlap_sparsity=.FALSE.,&
                     !        eps_filter=almo_scf_env%eps_filter)

                     !IF (unit_nr>0) THEN
                     !   WRITE(unit_nr,*) "Inverting blocked overlap matrix of discarded virtuals"
                     !ENDIF
                     !CALL invert_Hotelling(k_vd_index_up,&
                     !        k_vd_index_down,&
                     !        almo_scf_env%eps_filter)
                     !IF (safe_mode) THEN
                     !   CALL cp_dbcsr_init(matrix_tmp1)
                     !   CALL cp_dbcsr_create(matrix_tmp1,template=k_vd_index_down,&
                     !                        matrix_type=dbcsr_type_no_symmetry)
                     !   CALL cp_dbcsr_multiply("N","N",1.0_dp,k_vd_index_up,&
                     !                          k_vd_index_down,&
                     !                          0.0_dp, matrix_tmp1,&
                     !                          filter_eps=almo_scf_env%eps_filter)
                     !   frob_matrix_base=cp_dbcsr_frobenius_norm(matrix_tmp1)
                     !   CALL cp_dbcsr_add_on_diag(matrix_tmp1,-1.0_dp)
                     !   frob_matrix=cp_dbcsr_frobenius_norm(matrix_tmp1)
                     !   IF (unit_nr>0) THEN
                     !      WRITE(unit_nr,*) "Error for (inv(SIG)*SIG-I)",&
                     !            frob_matrix/frob_matrix_base
                     !   ENDIF
                     !   CALL cp_dbcsr_release(matrix_tmp1)
                     !ENDIF

                     ! init matrices necessary for optimization of truncated virts
                     ! init blocked gradient before setting K to zero
                     ! otherwise the block structure might be lost
                     CALL cp_dbcsr_init(grad)
                     CALL cp_dbcsr_create(grad, &
                                          template=almo_scf_env%matrix_k_blk(ispin))
                     CALL cp_dbcsr_copy(grad, almo_scf_env%matrix_k_blk(ispin))

                     ! init MD in the k-space
                     md_in_k_space = almo_scf_env%logical01
                     IF (md_in_k_space) THEN
                        CALL cp_dbcsr_init(velocity)
                        CALL cp_dbcsr_create(velocity, &
                                             template=almo_scf_env%matrix_k_blk(ispin))
                        CALL cp_dbcsr_copy(velocity, almo_scf_env%matrix_k_blk(ispin))
                        CALL cp_dbcsr_set(velocity, 0.0_dp)
                        time_step = almo_scf_env%opt_k_trial_step_size
                     ENDIF

                     CALL cp_dbcsr_init(prev_step)
                     CALL cp_dbcsr_create(prev_step, &
                                          template=almo_scf_env%matrix_k_blk(ispin))

                     CALL cp_dbcsr_init(prev_minus_prec_grad)
                     CALL cp_dbcsr_create(prev_minus_prec_grad, &
                                          template=almo_scf_env%matrix_k_blk(ispin))

                     ! initialize diagonal blocks of the preconditioner to 1.0_dp
                     CALL cp_dbcsr_init(prec)
                     CALL cp_dbcsr_create(prec, &
                                          template=almo_scf_env%matrix_k_blk(ispin))
                     CALL cp_dbcsr_copy(prec, almo_scf_env%matrix_k_blk(ispin))
                     CALL cp_dbcsr_set(prec, 1.0_dp)

                     ! generate initial K (extrapolate if previous values are available)
                     CALL cp_dbcsr_set(almo_scf_env%matrix_k_blk(ispin), 0.0_dp)
                     ! matrix_k_central stores current k because matrix_k_blk is updated
                     ! during linear search
                     CALL cp_dbcsr_init(matrix_k_central)
                     CALL cp_dbcsr_create(matrix_k_central, &
                                          template=almo_scf_env%matrix_k_blk(ispin))
                     CALL cp_dbcsr_copy(matrix_k_central, &
                                        almo_scf_env%matrix_k_blk(ispin))
                     CALL cp_dbcsr_init(tmp_k_blk)
                     CALL cp_dbcsr_create(tmp_k_blk, &
                                          template=almo_scf_env%matrix_k_blk(ispin))
                     CALL cp_dbcsr_init(step)
                     CALL cp_dbcsr_create(step, &
                                          template=almo_scf_env%matrix_k_blk(ispin))
                     CALL cp_dbcsr_set(step, 0.0_dp)
                     CALL cp_dbcsr_init(t_curr)
                     CALL cp_dbcsr_create(t_curr, &
                                          template=almo_scf_env%matrix_t(ispin))
                     CALL cp_dbcsr_init(sigma_oo_curr)
                     CALL cp_dbcsr_init(sigma_oo_curr_inv)
                     CALL cp_dbcsr_create(sigma_oo_curr, &
                                          template=almo_scf_env%matrix_sigma(ispin), &
                                          matrix_type=dbcsr_type_no_symmetry)
                     CALL cp_dbcsr_create(sigma_oo_curr_inv, &
                                          template=almo_scf_env%matrix_sigma(ispin), &
                                          matrix_type=dbcsr_type_no_symmetry)
                     CALL cp_dbcsr_init(tmp1_n_vr)
                     CALL cp_dbcsr_create(tmp1_n_vr, &
                                          template=almo_scf_env%matrix_v(ispin))
                     CALL cp_dbcsr_init(tmp3_vd_vr)
                     CALL cp_dbcsr_create(tmp3_vd_vr, &
                                          template=almo_scf_env%matrix_k_blk(ispin))
                     CALL cp_dbcsr_init(tmp2_n_o)
                     CALL cp_dbcsr_create(tmp2_n_o, &
                                          template=almo_scf_env%matrix_t(ispin))
                     CALL cp_dbcsr_init(tmp4_o_vr)
                     CALL cp_dbcsr_create(tmp4_o_vr, &
                                          template=almo_scf_env%matrix_ov(ispin))
                     CALL cp_dbcsr_init(prev_grad)
                     CALL cp_dbcsr_create(prev_grad, &
                                          template=almo_scf_env%matrix_k_blk(ispin))
                     CALL cp_dbcsr_set(prev_grad, 0.0_dp)

                     !CALL cp_dbcsr_init(sigma_oo_guess)
                     !CALL cp_dbcsr_create(sigma_oo_guess,&
                     !        template=almo_scf_env%matrix_sigma(ispin),&
                     !        matrix_type=dbcsr_type_no_symmetry)
                     !CALL cp_dbcsr_set(sigma_oo_guess,0.0_dp)
                     !CALL cp_dbcsr_add_on_diag(sigma_oo_guess,1.0_dp)
                     !CALL cp_dbcsr_filter(sigma_oo_guess,almo_scf_env%eps_filter)
                     !CALL cp_dbcsr_print(sigma_oo_guess)

                  ENDIF ! done constructing discarded virtuals

                  ! init variables
                  opt_k_max_iter = almo_scf_env%opt_k_max_iter
                  iteration = 0
                  converged = .FALSE.
                  prepare_to_exit = .FALSE.
                  beta = 0.0_dp
                  line_search = .FALSE.
                  obj_function = 0.0_dp
                  conjugacy_error = 0.0_dp
                  line_search_error = 0.0_dp
                  fun0 = 0.0_dp
                  fun1 = 0.0_dp
                  gfun0 = 0.0_dp
                  gfun1 = 0.0_dp
                  step_size_quadratic_approx = 0.0_dp
                  reset_step_size = .TRUE.
                  IF (almo_scf_env%deloc_truncate_virt .EQ. virt_full) opt_k_max_iter = 0

                  ! start cg iterations to optimize matrix_k_blk
                  DO

                     CALL timeset('k_opt_vr', handle1)

                     IF (almo_scf_env%deloc_truncate_virt .NE. virt_full) THEN

                        ! construct k-excited virtuals
                        CALL cp_dbcsr_multiply("N", "N", 1.0_dp, vd_fixed, &
                                               almo_scf_env%matrix_k_blk(ispin), &
                                               0.0_dp, almo_scf_env%matrix_v(ispin), &
                                               filter_eps=almo_scf_env%eps_filter)
                        CALL cp_dbcsr_add(almo_scf_env%matrix_v(ispin), vr_fixed, &
                                          +1.0_dp, +1.0_dp)
                     ENDIF

                     ! decompose the overlap matrix of the current retained orbitals
                     !IF (unit_nr>0) THEN
                     !   WRITE(unit_nr,*) "decompose the active VV overlap matrix"
                     !ENDIF
                     CALL get_overlap(bra=almo_scf_env%matrix_v(ispin), &
                                      ket=almo_scf_env%matrix_v(ispin), &
                                      overlap=almo_scf_env%matrix_sigma_vv(ispin), &
                                      metric=almo_scf_env%matrix_s(1), &
                                      retain_overlap_sparsity=.FALSE., &
                                      eps_filter=almo_scf_env%eps_filter)
                     ! use either cholesky or sqrt
                     !! RZK-warning: strangely, cholesky does not work with k-optimization
                     IF (almo_scf_env%deloc_truncate_virt .EQ. virt_full) THEN
                        CALL timeset('cholesky', handle2)
                        t1cholesky = m_walltime()

                        ! re-create sigma_vv_sqrt because desymmetrize is buggy -
                        ! it will create multiple copies of blocks
                        CALL cp_dbcsr_create(sigma_vv_sqrt, &
                                             template=almo_scf_env%matrix_sigma_vv(ispin), &
                                             matrix_type=dbcsr_type_no_symmetry)
                        CALL cp_dbcsr_desymmetrize(almo_scf_env%matrix_sigma_vv(ispin), &
                                                   sigma_vv_sqrt)
                        CALL cp_dbcsr_cholesky_decompose(sigma_vv_sqrt, &
                                                         para_env=almo_scf_env%para_env, &
                                                         blacs_env=almo_scf_env%blacs_env)
                        CALL cp_dbcsr_triu(sigma_vv_sqrt)
                        CALL cp_dbcsr_filter(sigma_vv_sqrt, almo_scf_env%eps_filter)
                        ! apply SOLVE to compute U^(-1) : U*U^(-1)=I
                        CALL cp_dbcsr_get_info(sigma_vv_sqrt, nfullrows_total=n)
                        CALL cp_dbcsr_init(matrix_tmp1)
                        CALL cp_dbcsr_create(matrix_tmp1, template=almo_scf_env%matrix_sigma_vv(ispin), &
                                             matrix_type=dbcsr_type_no_symmetry)
                        CALL cp_dbcsr_set(matrix_tmp1, 0.0_dp)
                        CALL cp_dbcsr_add_on_diag(matrix_tmp1, 1.0_dp)
                        CALL cp_dbcsr_cholesky_restore(matrix_tmp1, n, sigma_vv_sqrt, &
                                                       sigma_vv_sqrt_inv, op="SOLVE", pos="RIGHT", &
                                                       para_env=almo_scf_env%para_env, &
                                                       blacs_env=almo_scf_env%blacs_env)
                        CALL cp_dbcsr_filter(sigma_vv_sqrt_inv, almo_scf_env%eps_filter)
                        CALL cp_dbcsr_release(matrix_tmp1)
                        IF (safe_mode) THEN
                           CALL cp_dbcsr_create(matrix_tmp1, template=almo_scf_env%matrix_sigma_vv(ispin), &
                                                matrix_type=dbcsr_type_no_symmetry)
                           CALL cp_dbcsr_desymmetrize(almo_scf_env%matrix_sigma_vv(ispin), &
                                                      matrix_tmp1)
                           CALL cp_dbcsr_multiply("T", "N", 1.0_dp, sigma_vv_sqrt, &
                                                  sigma_vv_sqrt, &
                                                  -1.0_dp, matrix_tmp1, filter_eps=almo_scf_env%eps_filter)
                           frob_matrix = cp_dbcsr_frobenius_norm(matrix_tmp1)
                           CALL cp_dbcsr_add_on_diag(matrix_tmp1, 1.0_dp)
                           frob_matrix_base = cp_dbcsr_frobenius_norm(matrix_tmp1)
                           IF (unit_nr > 0) THEN
                              WRITE (unit_nr, *) "Error for ( U^T * U - Sig )", &
                                 frob_matrix/frob_matrix_base
                           ENDIF
                           CALL cp_dbcsr_multiply("N", "N", 1.0_dp, sigma_vv_sqrt_inv, &
                                                  sigma_vv_sqrt, &
                                                  0.0_dp, matrix_tmp1, filter_eps=almo_scf_env%eps_filter)
                           frob_matrix_base = cp_dbcsr_frobenius_norm(matrix_tmp1)
                           CALL cp_dbcsr_add_on_diag(matrix_tmp1, -1.0_dp)
                           frob_matrix = cp_dbcsr_frobenius_norm(matrix_tmp1)
                           IF (unit_nr > 0) THEN
                              WRITE (unit_nr, *) "Error for ( inv(U) * U - I )", &
                                 frob_matrix/frob_matrix_base
                           ENDIF
                           CALL cp_dbcsr_release(matrix_tmp1)
                        ENDIF ! safe_mode
                        t2cholesky = m_walltime()
                        IF (unit_nr > 0) THEN
                           WRITE (unit_nr, *) "Cholesky+inverse wall-time: ", t2cholesky-t1cholesky
                        ENDIF
                        CALL timestop(handle2)
                     ELSE
                        CALL matrix_sqrt_Newton_Schulz(sigma_vv_sqrt, &
                                                       sigma_vv_sqrt_inv, &
                                                       almo_scf_env%matrix_sigma_vv(ispin), &
                                                       !matrix_sqrt_inv_guess=sigma_vv_sqrt_inv_guess,&
                                                       !matrix_sqrt_guess=sigma_vv_sqrt_guess,&
                                                       threshold=almo_scf_env%eps_filter, &
                                                       order=almo_scf_env%order_lanczos, &
                                                       eps_lanczos=almo_scf_env%eps_lanczos, &
                                                       max_iter_lanczos=almo_scf_env%max_iter_lanczos)
                        CALL cp_dbcsr_copy(sigma_vv_sqrt_inv_guess, sigma_vv_sqrt_inv)
                        CALL cp_dbcsr_copy(sigma_vv_sqrt_guess, sigma_vv_sqrt)
                        IF (safe_mode) THEN
                           CALL cp_dbcsr_init(matrix_tmp1)
                           CALL cp_dbcsr_create(matrix_tmp1, template=almo_scf_env%matrix_sigma_vv(ispin), &
                                                matrix_type=dbcsr_type_no_symmetry)
                           CALL cp_dbcsr_init(matrix_tmp2)
                           CALL cp_dbcsr_create(matrix_tmp2, template=almo_scf_env%matrix_sigma_vv(ispin), &
                                                matrix_type=dbcsr_type_no_symmetry)

                           CALL cp_dbcsr_multiply("N", "N", 1.0_dp, sigma_vv_sqrt_inv, &
                                                  almo_scf_env%matrix_sigma_vv(ispin), &
                                                  0.0_dp, matrix_tmp1, filter_eps=almo_scf_env%eps_filter)
                           CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1, &
                                                  sigma_vv_sqrt_inv, &
                                                  0.0_dp, matrix_tmp2, filter_eps=almo_scf_env%eps_filter)

                           frob_matrix_base = cp_dbcsr_frobenius_norm(matrix_tmp2)
                           CALL cp_dbcsr_add_on_diag(matrix_tmp2, -1.0_dp)
                           frob_matrix = cp_dbcsr_frobenius_norm(matrix_tmp2)
                           IF (unit_nr > 0) THEN
                              WRITE (unit_nr, *) "Error for (inv(sqrt(SIGVV))*SIGVV*inv(sqrt(SIGVV))-I)", &
                                 frob_matrix/frob_matrix_base
                           ENDIF

                           CALL cp_dbcsr_release(matrix_tmp1)
                           CALL cp_dbcsr_release(matrix_tmp2)
                        ENDIF
                     ENDIF
                     CALL timestop(handle1)

                     ! compute excitation amplitudes (to the current set of retained virtuals)
                     ! set convergence criterion for x-optimization
                     IF ((iteration .EQ. 0) .AND. (.NOT. line_search) .AND. &
                         (outer_opt_k_iteration .EQ. 0)) THEN
                        x_opt_eps_adaptive = &
                           almo_scf_env%deloc_cayley_eps_convergence
                     ELSE
                        x_opt_eps_adaptive = &
                           MAX(ABS(almo_scf_env%deloc_cayley_eps_convergence), &
                               ABS(x_opt_eps_adaptive_factor*grad_norm))
                     ENDIF
                     CALL ct_step_env_init(ct_step_env)
                     CALL ct_step_env_set(ct_step_env, &
                                          para_env=almo_scf_env%para_env, &
                                          blacs_env=almo_scf_env%blacs_env, &
                                          use_occ_orbs=.TRUE., &
                                          use_virt_orbs=.TRUE., &
                                          occ_orbs_orthogonal=.FALSE., &
                                          virt_orbs_orthogonal=.FALSE., &
                                          pp_preconditioner_full=almo_scf_env%deloc_cayley_occ_precond, &
                                          qq_preconditioner_full=almo_scf_env%deloc_cayley_vir_precond, &
                                          tensor_type=almo_scf_env%deloc_cayley_tensor_type, &
                                          neglect_quadratic_term=almo_scf_env%deloc_cayley_linear, &
                                          conjugator=almo_scf_env%deloc_cayley_conjugator, &
                                          max_iter=almo_scf_env%deloc_cayley_max_iter, &
                                          calculate_energy_corr=.TRUE., &
                                          update_p=.FALSE., &
                                          update_q=.FALSE., &
                                          eps_convergence=x_opt_eps_adaptive, &
                                          eps_filter=almo_scf_env%eps_filter, &
                                          !nspins=1,&
                                          q_index_up=sigma_vv_sqrt_inv, &
                                          q_index_down=sigma_vv_sqrt, &
                                          p_index_up=almo_scf_env%matrix_sigma_sqrt_inv(ispin), &
                                          p_index_down=almo_scf_env%matrix_sigma_sqrt(ispin), &
                                          matrix_ks=almo_scf_env%matrix_ks_0deloc(ispin), &
                                          matrix_t=almo_scf_env%matrix_t(ispin), &
                                          matrix_qp_template=almo_scf_env%matrix_vo(ispin), &
                                          matrix_pq_template=almo_scf_env%matrix_ov(ispin), &
                                          matrix_v=almo_scf_env%matrix_v(ispin), &
                                          matrix_x_guess=almo_scf_env%matrix_x(ispin))
                     ! perform calculations
                     CALL ct_step_execute(ct_step_env)
                     ! get the energy correction
                     CALL ct_step_env_get(ct_step_env, &
                                          energy_correction=energy_correction(ispin), &
                                          copy_matrix_x=almo_scf_env%matrix_x(ispin))
                     CALL ct_step_env_clean(ct_step_env)
                     ! RZK-warning matrix_x is being transformed
                     ! back and forth between orth and up_down representations
                     energy_correction(1) = energy_correction(1)*spin_factor

                     IF (opt_k_max_iter .NE. 0) THEN

                        CALL timeset('k_opt_t_curr', handle3)

                        ! construct current occupied orbitals T_blk + V_r*X
                        CALL cp_dbcsr_multiply("N", "N", 1.0_dp, &
                                               almo_scf_env%matrix_v(ispin), &
                                               almo_scf_env%matrix_x(ispin), &
                                               0.0_dp, t_curr, &
                                               filter_eps=almo_scf_env%eps_filter)
                        CALL cp_dbcsr_add(t_curr, almo_scf_env%matrix_t_blk(ispin), &
                                          +1.0_dp, +1.0_dp)

                        ! calculate current occupied overlap
                        !IF (unit_nr>0) THEN
                        !   WRITE(unit_nr,*) "Inverting current occ overlap matrix"
                        !ENDIF
                        CALL get_overlap(bra=t_curr, &
                                         ket=t_curr, &
                                         overlap=sigma_oo_curr, &
                                         metric=almo_scf_env%matrix_s(1), &
                                         retain_overlap_sparsity=.FALSE., &
                                         eps_filter=almo_scf_env%eps_filter)
                        IF (iteration .EQ. 0) THEN
                           CALL invert_Hotelling(sigma_oo_curr_inv, &
                                                 sigma_oo_curr, &
                                                 threshold=almo_scf_env%eps_filter, &
                                                 use_inv_as_guess=.FALSE.)
                        ELSE
                           CALL invert_Hotelling(sigma_oo_curr_inv, &
                                                 sigma_oo_curr, &
                                                 threshold=almo_scf_env%eps_filter, &
                                                 use_inv_as_guess=.TRUE.)
                           !CALL cp_dbcsr_copy(sigma_oo_guess,sigma_oo_curr_inv)
                        ENDIF
                        IF (safe_mode) THEN
                           CALL cp_dbcsr_init(matrix_tmp1)
                           CALL cp_dbcsr_create(matrix_tmp1, template=sigma_oo_curr, &
                                                matrix_type=dbcsr_type_no_symmetry)
                           CALL cp_dbcsr_multiply("N", "N", 1.0_dp, sigma_oo_curr, &
                                                  sigma_oo_curr_inv, &
                                                  0.0_dp, matrix_tmp1, &
                                                  filter_eps=almo_scf_env%eps_filter)
                           frob_matrix_base = cp_dbcsr_frobenius_norm(matrix_tmp1)
                           CALL cp_dbcsr_add_on_diag(matrix_tmp1, -1.0_dp)
                           frob_matrix = cp_dbcsr_frobenius_norm(matrix_tmp1)
                           !CALL cp_dbcsr_filter(matrix_tmp1,almo_scf_env%eps_filter)
                           !CALL cp_dbcsr_print(matrix_tmp1)
                           IF (unit_nr > 0) THEN
                              WRITE (unit_nr, *) "Error for (SIG*inv(SIG)-I)", &
                                 frob_matrix/frob_matrix_base, frob_matrix_base
                           ENDIF
                           CALL cp_dbcsr_release(matrix_tmp1)
                        ENDIF
                        IF (safe_mode) THEN
                           CALL cp_dbcsr_init(matrix_tmp1)
                           CALL cp_dbcsr_create(matrix_tmp1, template=sigma_oo_curr, &
                                                matrix_type=dbcsr_type_no_symmetry)
                           CALL cp_dbcsr_multiply("N", "N", 1.0_dp, sigma_oo_curr_inv, &
                                                  sigma_oo_curr, &
                                                  0.0_dp, matrix_tmp1, &
                                                  filter_eps=almo_scf_env%eps_filter)
                           frob_matrix_base = cp_dbcsr_frobenius_norm(matrix_tmp1)
                           CALL cp_dbcsr_add_on_diag(matrix_tmp1, -1.0_dp)
                           frob_matrix = cp_dbcsr_frobenius_norm(matrix_tmp1)
                           !CALL cp_dbcsr_filter(matrix_tmp1,almo_scf_env%eps_filter)
                           !CALL cp_dbcsr_print(matrix_tmp1)
                           IF (unit_nr > 0) THEN
                              WRITE (unit_nr, *) "Error for (inv(SIG)*SIG-I)", &
                                 frob_matrix/frob_matrix_base, frob_matrix_base
                           ENDIF
                           CALL cp_dbcsr_release(matrix_tmp1)
                        ENDIF

                        CALL timestop(handle3)
                        CALL timeset('k_opt_vd', handle4)

                        ! construct current discarded virtuals:
                        ! (1-R_curr)(1-Q^VR_curr)|ALMO_vd_basis> =
                        ! = (1-Q^VR_curr)|ALMO_vd_basis>
                        ! use sigma_vv_sqrt to store the inverse of the overlap
                        ! sigma_vv_inv is computed from sqrt/cholesky
                        CALL cp_dbcsr_multiply("N", "T", 1.0_dp, &
                                               sigma_vv_sqrt_inv, &
                                               sigma_vv_sqrt_inv, &
                                               0.0_dp, sigma_vv_sqrt, &
                                               filter_eps=almo_scf_env%eps_filter)
                        CALL apply_projector(psi_in=almo_scf_env%matrix_v_disc_blk(ispin), &
                                             psi_out=almo_scf_env%matrix_v_disc(ispin), &
                                             psi_projector=almo_scf_env%matrix_v(ispin), &
                                             metric=almo_scf_env%matrix_s(1), &
                                             project_out=.FALSE., &
                                             psi_projector_orthogonal=.FALSE., &
                                             proj_in_template=almo_scf_env%matrix_k_tr(ispin), &
                                             eps_filter=almo_scf_env%eps_filter, &
                                             sig_inv_projector=sigma_vv_sqrt)
                        !sig_inv_template=almo_scf_env%matrix_sigma_vv(ispin),&
                        CALL cp_dbcsr_add(almo_scf_env%matrix_v_disc(ispin), &
                                          vd_fixed, -1.0_dp, +1.0_dp)

                        CALL timestop(handle4)
                        CALL timeset('k_opt_grad', handle5)

                        ! evaluate the gradient from the assembled components
                        ! grad_xx = c0 [ (Vd_curr^tr)*F*T_curr*sigma_oo_curr_inv*(X^tr)]_xx
                        ! save previous gradient to calculate conjugation coef
                        IF (line_search) THEN
                           CALL cp_dbcsr_copy(prev_grad, grad)
                        ENDIF
                        CALL cp_dbcsr_multiply("N", "N", 1.0_dp, &
                                               almo_scf_env%matrix_ks_0deloc(ispin), &
                                               t_curr, &
                                               0.0_dp, tmp2_n_o, &
                                               filter_eps=almo_scf_env%eps_filter)
                        CALL cp_dbcsr_multiply("N", "T", 1.0_dp, &
                                               sigma_oo_curr_inv, &
                                               almo_scf_env%matrix_x(ispin), &
                                               0.0_dp, tmp4_o_vr, &
                                               filter_eps=almo_scf_env%eps_filter)
                        CALL cp_dbcsr_multiply("N", "N", 1.0_dp, &
                                               tmp2_n_o, &
                                               tmp4_o_vr, &
                                               0.0_dp, tmp1_n_vr, &
                                               filter_eps=almo_scf_env%eps_filter)
                        CALL cp_dbcsr_multiply("T", "N", 2.0_dp*spin_factor, &
                                               almo_scf_env%matrix_v_disc(ispin), &
                                               tmp1_n_vr, &
                                               0.0_dp, grad, &
                                               retain_sparsity=.TRUE.)
                        !filter_eps=almo_scf_env%eps_filter,&
                        ! keep tmp2_n_o for the next step
                        ! keep tmp4_o_vr for the preconditioner

                        ! check convergence and other exit criteria
                        grad_norm_frob = cp_dbcsr_frobenius_norm(grad)
                        CALL cp_dbcsr_norm(grad, dbcsr_norm_maxabsnorm, norm_scalar=grad_norm)
                        converged = (grad_norm .LT. almo_scf_env%opt_k_eps_convergence)
                        IF (converged .OR. (iteration .GE. opt_k_max_iter)) THEN
                           prepare_to_exit = .TRUE.
                        ENDIF
                        CALL timestop(handle5)

                        IF (.NOT. prepare_to_exit) THEN

                           CALL timeset('k_opt_energy', handle6)

                           ! compute "energy" c0*Tr[sig_inv_oo*t*F*t]
                           CALL cp_dbcsr_multiply("T", "N", spin_factor, &
                                                  t_curr, &
                                                  tmp2_n_o, &
                                                  0.0_dp, sigma_oo_curr, &
                                                  filter_eps=almo_scf_env%eps_filter)
                           delta_obj_function = fun0
                           CALL cp_dbcsr_trace(sigma_oo_curr_inv, &
                                               sigma_oo_curr, obj_function, "T", "N")
                           delta_obj_function = obj_function-delta_obj_function
                           IF (line_search) THEN
                              fun1 = obj_function
                           ELSE
                              fun0 = obj_function
                           ENDIF

                           CALL timestop(handle6)

                           ! update the step direction
                           IF (.NOT. line_search) THEN

                              CALL timeset('k_opt_step', handle7)

                              IF ((.NOT. md_in_k_space) .AND. &
                                  (iteration .GE. MAX(0, almo_scf_env%opt_k_prec_iter_start) .AND. &
                                   MOD(iteration-almo_scf_env%opt_k_prec_iter_start, &
                                       almo_scf_env%opt_k_prec_iter_freq) .EQ. 0)) THEN

                                 !IF ((iteration.eq.0).AND.(.NOT.md_in_k_space)) THEN

                                 ! compute the preconditioner
                                 IF (unit_nr > 0) THEN
                                    WRITE (unit_nr, *) "Computing preconditioner"
                                 ENDIF
                                 !CALL opt_k_create_preconditioner(prec,&
                                 !        almo_scf_env%matrix_v_disc(ispin),&
                                 !        almo_scf_env%matrix_ks_0deloc(ispin),&
                                 !        almo_scf_env%matrix_x(ispin),&
                                 !        tmp4_o_vr,&
                                 !        almo_scf_env%matrix_s(1),&
                                 !        grad,&
                                 !        !almo_scf_env%matrix_v_disc_blk(ispin),&
                                 !        vd_fixed,&
                                 !        t_curr,&
                                 !        k_vd_index_up,&
                                 !        k_vr_index_down,&
                                 !        tmp1_n_vr,&
                                 !        spin_factor,&
                                 !        almo_scf_env%eps_filter)
                                 CALL opt_k_create_preconditioner_blk(almo_scf_env, &
                                                                      almo_scf_env%matrix_v_disc(ispin), &
                                                                      tmp4_o_vr, &
                                                                      t_curr, &
                                                                      ispin, &
                                                                      spin_factor)

                              ENDIF

                              ! save the previous step
                              CALL cp_dbcsr_copy(prev_step, step)

                              ! compute the new step
                              CALL opt_k_apply_preconditioner_blk(almo_scf_env, &
                                                                  step, grad, ispin)
                              !CALL cp_dbcsr_hadamard_product(prec,grad,step)
                              CALL cp_dbcsr_scale(step, -1.0_dp)

                              ! check whether we need to reset conjugate directions
                              reset_conjugator = .FALSE.
                              ! first check if manual reset is active
                              IF (iteration .LT. MAX(almo_scf_env%opt_k_conj_iter_start, 1) .OR. &
                                  MOD(iteration-almo_scf_env%opt_k_conj_iter_start, &
                                      almo_scf_env%opt_k_conj_iter_freq) .EQ. 0) THEN

                                 reset_conjugator = .TRUE.

                              ELSE

                                 ! check for the errors in the cg algorithm
                                 !CALL cp_dbcsr_hadamard_product(prec,prev_grad,tmp_k_blk)
                                 !CALL cp_dbcsr_trace(grad,tmp_k_blk,numer,"T","N")
                                 !CALL cp_dbcsr_trace(prev_grad,tmp_k_blk,denom,"T","N")
                                 CALL cp_dbcsr_trace(grad, prev_minus_prec_grad, numer, "T", "N")
                                 CALL cp_dbcsr_trace(prev_grad, prev_minus_prec_grad, denom, "T", "N")
                                 conjugacy_error = numer/denom

                                 IF (conjugacy_error .GT. MIN(0.5_dp, conjugacy_error_threshold)) THEN
                                    reset_conjugator = .TRUE.
                                    IF (unit_nr > 0) THEN
                                       WRITE (unit_nr, *) "Lack of progress, conjugacy error is ", conjugacy_error
                                    ENDIF
                                 ENDIF

                                 ! check the gradient along the previous direction
                                 IF ((iteration .NE. 0) .AND. (.NOT. reset_conjugator)) THEN
                                    CALL cp_dbcsr_trace(grad, prev_step, numer, "T", "N")
                                    CALL cp_dbcsr_trace(prev_grad, prev_step, denom, "T", "N")
                                    line_search_error = numer/denom
                                    IF (line_search_error .GT. line_search_error_threshold) THEN
                                       reset_conjugator = .TRUE.
                                       IF (unit_nr > 0) THEN
                                          WRITE (unit_nr, *) "Bad line search, line search error is ", line_search_error
                                       ENDIF
                                    ENDIF
                                 ENDIF

                              ENDIF

                              ! compute the conjugation coefficient - beta
                              IF (.NOT. reset_conjugator) THEN

                                 SELECT CASE (almo_scf_env%opt_k_conjugator)
                                 CASE (cg_hestenes_stiefel)
                                    CALL cp_dbcsr_copy(tmp_k_blk, grad)
                                    CALL cp_dbcsr_add(tmp_k_blk, prev_grad, 1.0_dp, -1.0_dp)
                                    CALL cp_dbcsr_trace(tmp_k_blk, step, numer, "T", "N")
                                    CALL cp_dbcsr_trace(tmp_k_blk, prev_step, denom, "T", "N")
                                    beta = -1.0_dp*numer/denom
                                 CASE (cg_fletcher_reeves)
                                    !CALL cp_dbcsr_hadamard_product(prec,prev_grad,tmp_k_blk)
                                    !CALL cp_dbcsr_trace(prev_grad,tmp_k_blk,denom,"T","N")
                                    !CALL cp_dbcsr_hadamard_product(prec,grad,tmp_k_blk)
                                    !CALL cp_dbcsr_trace(grad,tmp_k_blk,numer,"T","N")
                                    !beta=numer/denom
                                    CALL cp_dbcsr_trace(grad, step, numer, "T", "N")
                                    CALL cp_dbcsr_trace(prev_grad, prev_minus_prec_grad, denom, "T", "N")
                                    beta = numer/denom
                                 CASE (cg_polak_ribiere)
                                    !CALL cp_dbcsr_hadamard_product(prec,prev_grad,tmp_k_blk)
                                    !CALL cp_dbcsr_trace(prev_grad,tmp_k_blk,denom,"T","N")
                                    !CALL cp_dbcsr_add(prev_grad,grad,-1.0_dp,1.0_dp)
                                    !CALL cp_dbcsr_hadamard_product(prec,prev_grad,tmp_k_blk)
                                    !CALL cp_dbcsr_trace(tmp_k_blk,grad,numer,"T","N")
                                    CALL cp_dbcsr_trace(prev_grad, prev_minus_prec_grad, denom, "T", "N")
                                    CALL cp_dbcsr_copy(tmp_k_blk, grad)
                                    CALL cp_dbcsr_add(tmp_k_blk, prev_grad, 1.0_dp, -1.0_dp)
                                    CALL cp_dbcsr_trace(tmp_k_blk, step, numer, "T", "N")
                                    beta = numer/denom
                                 CASE (cg_fletcher)
                                    !CALL cp_dbcsr_hadamard_product(prec,grad,tmp_k_blk)
                                    !CALL cp_dbcsr_trace(grad,tmp_k_blk,numer,"T","N")
                                    !CALL cp_dbcsr_trace(prev_grad,prev_step,denom,"T","N")
                                    !beta=-1.0_dp*numer/denom
                                    CALL cp_dbcsr_trace(grad, step, numer, "T", "N")
                                    CALL cp_dbcsr_trace(prev_grad, prev_step, denom, "T", "N")
                                    beta = numer/denom
                                 CASE (cg_liu_storey)
                                    CALL cp_dbcsr_trace(prev_grad, prev_step, denom, "T", "N")
                                    !CALL cp_dbcsr_add(prev_grad,grad,-1.0_dp,1.0_dp)
                                    !CALL cp_dbcsr_hadamard_product(prec,prev_grad,tmp_k_blk)
                                    !CALL cp_dbcsr_trace(tmp_k_blk,grad,numer,"T","N")
                                    CALL cp_dbcsr_copy(tmp_k_blk, grad)
                                    CALL cp_dbcsr_add(tmp_k_blk, prev_grad, 1.0_dp, -1.0_dp)
                                    CALL cp_dbcsr_trace(tmp_k_blk, step, numer, "T", "N")
                                    beta = numer/denom
                                 CASE (cg_dai_yuan)
                                    !CALL cp_dbcsr_hadamard_product(prec,grad,tmp_k_blk)
                                    !CALL cp_dbcsr_trace(grad,tmp_k_blk,numer,"T","N")
                                    !CALL cp_dbcsr_add(prev_grad,grad,-1.0_dp,1.0_dp)
                                    !CALL cp_dbcsr_trace(prev_grad,prev_step,denom,"T","N")
                                    !beta=numer/denom
                                    CALL cp_dbcsr_trace(grad, step, numer, "T", "N")
                                    CALL cp_dbcsr_copy(tmp_k_blk, grad)
                                    CALL cp_dbcsr_add(tmp_k_blk, prev_grad, 1.0_dp, -1.0_dp)
                                    CALL cp_dbcsr_trace(tmp_k_blk, prev_step, denom, "T", "N")
                                    beta = -1.0_dp*numer/denom
                                 CASE (cg_hager_zhang)
                                    !CALL cp_dbcsr_add(prev_grad,grad,-1.0_dp,1.0_dp)
                                    !CALL cp_dbcsr_trace(prev_grad,prev_step,denom,"T","N")
                                    !CALL cp_dbcsr_hadamard_product(prec,prev_grad,tmp_k_blk)
                                    !CALL cp_dbcsr_trace(tmp_k_blk,prev_grad,numer,"T","N")
                                    !kappa=2.0_dp*numer/denom
                                    !CALL cp_dbcsr_trace(tmp_k_blk,grad,numer,"T","N")
                                    !tau=numer/denom
                                    !CALL cp_dbcsr_trace(prev_step,grad,numer,"T","N")
                                    !beta=tau-kappa*numer/denom
                                    CALL cp_dbcsr_copy(tmp_k_blk, grad)
                                    CALL cp_dbcsr_add(tmp_k_blk, prev_grad, 1.0_dp, -1.0_dp)
                                    CALL cp_dbcsr_trace(tmp_k_blk, prev_step, denom, "T", "N")
                                    CALL cp_dbcsr_trace(tmp_k_blk, prev_minus_prec_grad, numer, "T", "N")
                                    kappa = -2.0_dp*numer/denom
                                    CALL cp_dbcsr_trace(tmp_k_blk, step, numer, "T", "N")
                                    tau = -1.0_dp*numer/denom
                                    CALL cp_dbcsr_trace(prev_step, grad, numer, "T", "N")
                                    beta = tau-kappa*numer/denom
                                 CASE (cg_zero)
                                    beta = 0.0_dp
                                 CASE DEFAULT
                                    CPABORT("illegal conjugator")
                                 END SELECT

                                 IF (beta .LT. 0.0_dp) THEN
                                    IF (unit_nr > 0) THEN
                                       WRITE (unit_nr, *) "Beta is negative, ", beta
                                    ENDIF
                                    reset_conjugator = .TRUE.
                                 ENDIF

                              ENDIF

                              IF (md_in_k_space) THEN
                                 reset_conjugator = .TRUE.
                              ENDIF

                              IF (reset_conjugator) THEN

                                 beta = 0.0_dp
                                 !reset_step_size=.TRUE.

                                 IF (unit_nr > 0) THEN
                                    WRITE (unit_nr, *) "(Re)-setting conjugator to zero"
                                 ENDIF

                              ENDIF

                              ! save the preconditioned gradient
                              CALL cp_dbcsr_copy(prev_minus_prec_grad, step)

                              ! conjugate the step direction
                              CALL cp_dbcsr_add(step, prev_step, 1.0_dp, beta)

                              CALL timestop(handle7)

                              ! update the step direction
                           ELSE ! step update
                              conjugacy_error = 0.0_dp
                           ENDIF

                           ! compute the gradient with respect to the step size in the curr direction
                           IF (line_search) THEN
                              CALL cp_dbcsr_trace(grad, step, gfun1, "T", "N")
                              line_search_error = gfun1/gfun0
                           ELSE
                              CALL cp_dbcsr_trace(grad, step, gfun0, "T", "N")
                           ENDIF

                           ! make a step - update k
                           IF (line_search) THEN

                              ! check if the trial step provides enough numerical accuracy
                              safety_multiplier = 1.0E+1_dp ! must be more than one
                              num_threshold = MAX(EPSILON(1.0_dp), &
                                                  safety_multiplier*(almo_scf_env%eps_filter**2)*almo_scf_env%ndomains)
                              IF (ABS(fun1-fun0-gfun0*step_size) .LT. num_threshold) THEN
                                 IF (unit_nr > 0) THEN
                                    WRITE (unit_nr, '(T3,A,1X,E17.7)') &
                                       "Numerical accuracy is too low to observe non-linear behavior", &
                                       ABS(fun1-fun0-gfun0*step_size)
                                    WRITE (unit_nr, '(T3,A,1X,E17.7,A,1X,E12.3)') "Error computing ", &
                                       ABS(gfun0), &
                                       " is smaller than the threshold", num_threshold
                                 ENDIF
                                 CPABORT("")
                              ENDIF
                              IF (ABS(gfun0) .LT. num_threshold) THEN
                                 IF (unit_nr > 0) THEN
                                    WRITE (unit_nr, '(T3,A,1X,E17.7,A,1X,E12.3)') "Linear gradient", &
                                       ABS(gfun0), &
                                       " is smaller than the threshold", num_threshold
                                 ENDIF
                                 CPABORT("")
                              ENDIF

                              use_quadratic_approximation = .TRUE.
                              use_cubic_approximation = .FALSE.

                              ! find the minimum assuming quadratic form
                              ! use f0, f1, g0
                              step_size_quadratic_approx = -(gfun0*step_size*step_size)/(2.0_dp*(fun1-fun0-gfun0*step_size))
                              ! use f0, f1, g1
                              step_size_quadratic_approx2 = -(fun1-fun0-step_size*gfun1/2.0_dp)/(gfun1-(fun1-fun0)/step_size)

                              IF ((step_size_quadratic_approx .LT. 0.0_dp) .AND. &
                                  (step_size_quadratic_approx2 .LT. 0.0_dp)) THEN
                                 IF (unit_nr > 0) THEN
                                    WRITE (unit_nr, '(T3,A,1X,E17.7,1X,E17.7,1X,A)') &
                                       "Quadratic approximation gives negative steps", &
                                       step_size_quadratic_approx, step_size_quadratic_approx2, &
                                       "trying cubic..."
                                 ENDIF
                                 use_cubic_approximation = .TRUE.
                                 use_quadratic_approximation = .FALSE.
                              ELSE
                                 IF (step_size_quadratic_approx .LT. 0.0_dp) THEN
                                    step_size_quadratic_approx = step_size_quadratic_approx2
                                 ENDIF
                                 IF (step_size_quadratic_approx2 .LT. 0.0_dp) THEN
                                    step_size_quadratic_approx2 = step_size_quadratic_approx
                                 ENDIF
                              ENDIF

                              ! check accuracy of the quadratic approximation
                              IF (use_quadratic_approximation) THEN
                                 quadratic_approx_error = ABS(step_size_quadratic_approx- &
                                                              step_size_quadratic_approx2)/step_size_quadratic_approx
                                 IF (quadratic_approx_error .GT. quadratic_approx_error_threshold) THEN
                                    IF (unit_nr > 0) THEN
                                       WRITE (unit_nr, '(T3,A,1X,E17.7,1X,E17.7,1X,A)') "Quadratic approximation is poor", &
                                          step_size_quadratic_approx, step_size_quadratic_approx2, &
                                          "Try cubic approximation"
                                    ENDIF
                                    use_cubic_approximation = .TRUE.
                                    use_quadratic_approximation = .FALSE.
                                 ENDIF
                              ENDIF

                              ! check if numerics is fine enough to capture the cubic form
                              IF (use_cubic_approximation) THEN

                                 ! if quadratic approximation is not accurate enough
                                 ! try to find the minimum assuming cubic form
                                 ! aa*x**3 + bb*x**2 + cc*x + dd = f(x)
                                 bb = (-step_size*gfun1+3.0_dp*(fun1-fun0)-2.0_dp*step_size*gfun0)/(step_size*step_size)
                                 aa = (gfun1-2.0_dp*step_size*bb-gfun0)/(3.0_dp*step_size*step_size)

                                 IF (ABS(gfun1-2.0_dp*step_size*bb-gfun0) .LT. num_threshold) THEN
                                    IF (unit_nr > 0) THEN
                                       WRITE (unit_nr, '(T3,A,1X,E17.7)') &
                                          "Numerical accuracy is too low to observe cubic behavior", &
                                          ABS(gfun1-2.0_dp*step_size*bb-gfun0)
                                    ENDIF
                                    use_cubic_approximation = .FALSE.
                                    use_quadratic_approximation = .TRUE.
                                 ENDIF
                                 IF (ABS(gfun1) .LT. num_threshold) THEN
                                    IF (unit_nr > 0) THEN
                                       WRITE (unit_nr, '(T3,A,1X,E17.7,A,1X,E12.3)') "Linear gradient", &
                                          ABS(gfun1), &
                                          " is smaller than the threshold", num_threshold
                                    ENDIF
                                    use_cubic_approximation = .FALSE.
                                    use_quadratic_approximation = .TRUE.
                                 ENDIF
                              ENDIF

                              ! find the step assuming cubic approximation
                              IF (use_cubic_approximation) THEN
                                 ! to obtain the minimum of the cubic function solve the quadratic equation
                                 ! 0.0*x**3 + 3.0*aa*x**2 + 2.0*bb*x + cc = 0
                                 CALL analytic_line_search(0.0_dp, 3.0_dp*aa, 2.0_dp*bb, gfun0, minima, nmins)
                                 IF (nmins .LT. 1) THEN
                                    IF (unit_nr > 0) THEN
                                       WRITE (unit_nr, '(T3,A)') &
                                          "Cubic approximation gives zero soultions! Use quadratic approximation"
                                    ENDIF
                                    use_quadratic_approximation = .TRUE.
                                    use_cubic_approximation = .TRUE.
                                 ELSE
                                    step_size = minima(1)
                                    IF (nmins .GT. 1) THEN
                                       IF (unit_nr > 0) THEN
                                          WRITE (unit_nr, '(T3,A)') &
                                             "More than one solution found! Use quadratic approximation"
                                       ENDIF
                                       use_quadratic_approximation = .TRUE.
                                       use_cubic_approximation = .TRUE.
                                    ENDIF
                                 ENDIF
                              ENDIF

                              IF (use_quadratic_approximation) THEN ! use quadratic approximation
                                 IF (unit_nr > 0) THEN
                                    WRITE (unit_nr, '(T3,A)') "Use quadratic approximation"
                                 ENDIF
                                 step_size = (step_size_quadratic_approx+step_size_quadratic_approx2)*0.5_dp
                              ENDIF

                              ! one more check on the step size
                              IF (step_size .LT. 0.0_dp) THEN
                                 CPABORT("Negative step proposed")
                              ENDIF

                              CALL cp_dbcsr_copy(almo_scf_env%matrix_k_blk(ispin), &
                                                 matrix_k_central)
                              CALL cp_dbcsr_add(almo_scf_env%matrix_k_blk(ispin), &
                                                step, 1.0_dp, step_size)
                              CALL cp_dbcsr_copy(matrix_k_central, &
                                                 almo_scf_env%matrix_k_blk(ispin))
                              line_search = .FALSE.

                           ELSE

                              IF (md_in_k_space) THEN

                                 ! update velocities v(i) = v(i-1) + 0.5*dT*(a(i-1) + a(i))
                                 IF (iteration .NE. 0) THEN
                                    CALL cp_dbcsr_add(velocity, &
                                                      step, 1.0_dp, 0.5_dp*time_step)
                                    CALL cp_dbcsr_add(velocity, &
                                                      prev_step, 1.0_dp, 0.5_dp*time_step)
                                 ENDIF
                                 kin_energy = cp_dbcsr_frobenius_norm(velocity)
                                 kin_energy = 0.5_dp*kin_energy*kin_energy

                                 ! update positions k(i) = k(i-1) + dT*v(i-1) + 0.5*dT*dT*a(i-1)
                                 CALL cp_dbcsr_add(almo_scf_env%matrix_k_blk(ispin), &
                                                   velocity, 1.0_dp, time_step)
                                 CALL cp_dbcsr_add(almo_scf_env%matrix_k_blk(ispin), &
                                                   step, 1.0_dp, 0.5_dp*time_step*time_step)

                              ELSE

                                 IF (reset_step_size) THEN
                                    step_size = almo_scf_env%opt_k_trial_step_size
                                    reset_step_size = .FALSE.
                                 ELSE
                                    step_size = step_size*almo_scf_env%opt_k_trial_step_size_multiplier
                                 ENDIF
                                 CALL cp_dbcsr_copy(almo_scf_env%matrix_k_blk(ispin), &
                                                    matrix_k_central)
                                 CALL cp_dbcsr_add(almo_scf_env%matrix_k_blk(ispin), &
                                                   step, 1.0_dp, step_size)
                                 line_search = .TRUE.
                              ENDIF

                           ENDIF

                        ENDIF ! .NOT.prepare_to_exit

                        ! print the status of the optimization
                        t2a = m_walltime()
                        IF (unit_nr > 0) THEN
                           IF (md_in_k_space) THEN
                              WRITE (unit_nr, '(T6,A,1X,I5,1X,E12.3,E16.7,F15.9,F15.9,F15.9,E12.3,F15.9,F15.9,F8.3)') &
                                 "K iter CG", iteration, time_step, time_step*iteration, &
                                 energy_correction(ispin), obj_function, delta_obj_function, grad_norm, &
                                 kin_energy, kin_energy+obj_function, beta
                           ELSE
                              IF (line_search .OR. prepare_to_exit) THEN
                                 WRITE (unit_nr, '(T6,A,1X,I3,1X,E12.3,F16.10,F16.10,E12.3,E12.3,E12.3,F8.3,F8.3,F10.3)') &
                                    "K iter CG", iteration, step_size, &
                                    energy_correction(ispin), delta_obj_function, grad_norm, &
                                    gfun0, line_search_error, beta, conjugacy_error, t2a-t1a
                                 !(flop1+flop2)/(1.0E6_dp*(t2-t1))
                              ELSE
                                 WRITE (unit_nr, '(T6,A,1X,I3,1X,E12.3,F16.10,F16.10,E12.3,E12.3,E12.3,F8.3,F8.3,F10.3)') &
                                    "K iter LS", iteration, step_size, &
                                    energy_correction(ispin), delta_obj_function, grad_norm, &
                                    gfun1, line_search_error, beta, conjugacy_error, t2a-t1a
                                 !(flop1+flop2)/(1.0E6_dp*(t2-t1))
                              ENDIF
                           ENDIF
                           CALL m_flush(unit_nr)
                        ENDIF
                        t1a = m_walltime()

                     ELSE ! opt_k_max_iter .eq. 0
                        prepare_to_exit = .TRUE.
                     ENDIF ! opt_k_max_iter .ne. 0

                     IF (.NOT. line_search) iteration = iteration+1

                     IF (prepare_to_exit) EXIT

                  ENDDO ! end iterations on K

                  IF (converged .OR. (outer_opt_k_iteration .GE. outer_opt_k_max_iter)) THEN
                     outer_opt_k_prepare_to_exit = .TRUE.
                  ENDIF

                  IF (almo_scf_env%deloc_truncate_virt .NE. virt_full) THEN

                     IF (unit_nr > 0) THEN
                        WRITE (unit_nr, *) "Updating ALMO virtuals"
                     ENDIF

                     CALL timeset('k_opt_v0_update', handle8)

                     ! update retained ALMO virtuals to restart the cg iterations
                     CALL cp_dbcsr_multiply("N", "N", 1.0_dp, &
                                            almo_scf_env%matrix_v_disc_blk(ispin), &
                                            almo_scf_env%matrix_k_blk(ispin), &
                                            0.0_dp, vr_fixed, &
                                            filter_eps=almo_scf_env%eps_filter)
                     CALL cp_dbcsr_add(vr_fixed, almo_scf_env%matrix_v_blk(ispin), &
                                       +1.0_dp, +1.0_dp)

                     ! update discarded ALMO virtuals to restart the cg iterations
                     CALL cp_dbcsr_multiply("N", "T", 1.0_dp, &
                                            almo_scf_env%matrix_v_blk(ispin), &
                                            almo_scf_env%matrix_k_blk(ispin), &
                                            0.0_dp, vd_fixed, &
                                            filter_eps=almo_scf_env%eps_filter)
                     CALL cp_dbcsr_add(vd_fixed, almo_scf_env%matrix_v_disc_blk(ispin), &
                                       -1.0_dp, +1.0_dp)

                     ! orthogonalize new orbitals on fragments
                     CALL get_overlap(bra=vr_fixed, &
                                      ket=vr_fixed, &
                                      overlap=k_vr_index_down, &
                                      metric=almo_scf_env%matrix_s_blk(1), &
                                      retain_overlap_sparsity=.FALSE., &
                                      eps_filter=almo_scf_env%eps_filter)
                     CALL cp_dbcsr_init(vr_index_sqrt_inv)
                     CALL cp_dbcsr_create(vr_index_sqrt_inv, template=k_vr_index_down, &
                                          matrix_type=dbcsr_type_no_symmetry)
                     CALL cp_dbcsr_init(vr_index_sqrt)
                     CALL cp_dbcsr_create(vr_index_sqrt, template=k_vr_index_down, &
                                          matrix_type=dbcsr_type_no_symmetry)
                     CALL matrix_sqrt_Newton_Schulz(vr_index_sqrt, &
                                                    vr_index_sqrt_inv, &
                                                    k_vr_index_down, &
                                                    threshold=almo_scf_env%eps_filter, &
                                                    order=almo_scf_env%order_lanczos, &
                                                    eps_lanczos=almo_scf_env%eps_lanczos, &
                                                    max_iter_lanczos=almo_scf_env%max_iter_lanczos)
                     IF (safe_mode) THEN
                        CALL cp_dbcsr_init(matrix_tmp1)
                        CALL cp_dbcsr_create(matrix_tmp1, template=k_vr_index_down, &
                                             matrix_type=dbcsr_type_no_symmetry)
                        CALL cp_dbcsr_init(matrix_tmp2)
                        CALL cp_dbcsr_create(matrix_tmp2, template=k_vr_index_down, &
                                             matrix_type=dbcsr_type_no_symmetry)

                        CALL cp_dbcsr_multiply("N", "N", 1.0_dp, vr_index_sqrt_inv, &
                                               k_vr_index_down, &
                                               0.0_dp, matrix_tmp1, filter_eps=almo_scf_env%eps_filter)
                        CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1, &
                                               vr_index_sqrt_inv, &
                                               0.0_dp, matrix_tmp2, filter_eps=almo_scf_env%eps_filter)

                        frob_matrix_base = cp_dbcsr_frobenius_norm(matrix_tmp2)
                        CALL cp_dbcsr_add_on_diag(matrix_tmp2, -1.0_dp)
                        frob_matrix = cp_dbcsr_frobenius_norm(matrix_tmp2)
                        IF (unit_nr > 0) THEN
                           WRITE (unit_nr, *) "Error for (inv(sqrt(SIGVV))*SIGVV*inv(sqrt(SIGVV))-I)", &
                              frob_matrix/frob_matrix_base
                        ENDIF

                        CALL cp_dbcsr_release(matrix_tmp1)
                        CALL cp_dbcsr_release(matrix_tmp2)
                     ENDIF
                     CALL cp_dbcsr_multiply("N", "N", 1.0_dp, &
                                            vr_fixed, &
                                            vr_index_sqrt_inv, &
                                            0.0_dp, almo_scf_env%matrix_v_blk(ispin), &
                                            filter_eps=almo_scf_env%eps_filter)

                     CALL get_overlap(bra=vd_fixed, &
                                      ket=vd_fixed, &
                                      overlap=k_vd_index_down, &
                                      metric=almo_scf_env%matrix_s_blk(1), &
                                      retain_overlap_sparsity=.FALSE., &
                                      eps_filter=almo_scf_env%eps_filter)
                     CALL cp_dbcsr_init(vd_index_sqrt_inv)
                     CALL cp_dbcsr_create(vd_index_sqrt_inv, template=k_vd_index_down, &
                                          matrix_type=dbcsr_type_no_symmetry)
                     CALL cp_dbcsr_init(vd_index_sqrt)
                     CALL cp_dbcsr_create(vd_index_sqrt, template=k_vd_index_down, &
                                          matrix_type=dbcsr_type_no_symmetry)
                     CALL matrix_sqrt_Newton_Schulz(vd_index_sqrt, &
                                                    vd_index_sqrt_inv, &
                                                    k_vd_index_down, &
                                                    threshold=almo_scf_env%eps_filter, &
                                                    order=almo_scf_env%order_lanczos, &
                                                    eps_lanczos=almo_scf_env%eps_lanczos, &
                                                    max_iter_lanczos=almo_scf_env%max_iter_lanczos)
                     IF (safe_mode) THEN
                        CALL cp_dbcsr_init(matrix_tmp1)
                        CALL cp_dbcsr_create(matrix_tmp1, template=k_vd_index_down, &
                                             matrix_type=dbcsr_type_no_symmetry)
                        CALL cp_dbcsr_init(matrix_tmp2)
                        CALL cp_dbcsr_create(matrix_tmp2, template=k_vd_index_down, &
                                             matrix_type=dbcsr_type_no_symmetry)

                        CALL cp_dbcsr_multiply("N", "N", 1.0_dp, vd_index_sqrt_inv, &
                                               k_vd_index_down, &
                                               0.0_dp, matrix_tmp1, filter_eps=almo_scf_env%eps_filter)
                        CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1, &
                                               vd_index_sqrt_inv, &
                                               0.0_dp, matrix_tmp2, filter_eps=almo_scf_env%eps_filter)

                        frob_matrix_base = cp_dbcsr_frobenius_norm(matrix_tmp2)
                        CALL cp_dbcsr_add_on_diag(matrix_tmp2, -1.0_dp)
                        frob_matrix = cp_dbcsr_frobenius_norm(matrix_tmp2)
                        IF (unit_nr > 0) THEN
                           WRITE (unit_nr, *) "Error for (inv(sqrt(SIGVV))*SIGVV*inv(sqrt(SIGVV))-I)", &
                              frob_matrix/frob_matrix_base
                        ENDIF

                        CALL cp_dbcsr_release(matrix_tmp1)
                        CALL cp_dbcsr_release(matrix_tmp2)
                     ENDIF
                     CALL cp_dbcsr_multiply("N", "N", 1.0_dp, &
                                            vd_fixed, &
                                            vd_index_sqrt_inv, &
                                            0.0_dp, almo_scf_env%matrix_v_disc_blk(ispin), &
                                            filter_eps=almo_scf_env%eps_filter)

                     CALL cp_dbcsr_release(vr_index_sqrt_inv)
                     CALL cp_dbcsr_release(vr_index_sqrt)
                     CALL cp_dbcsr_release(vd_index_sqrt_inv)
                     CALL cp_dbcsr_release(vd_index_sqrt)

                     CALL timestop(handle8)

                  ENDIF ! ne.virt_full

                  ! RZK-warning released outside the outer loop
                  CALL cp_dbcsr_release(sigma_vv_sqrt)
                  CALL cp_dbcsr_release(sigma_vv_sqrt_inv)
                  IF (almo_scf_env%deloc_truncate_virt .NE. virt_full) THEN
                     CALL cp_dbcsr_release(k_vr_index_down)
                     CALL cp_dbcsr_release(k_vd_index_down)
                     !CALL cp_dbcsr_release(k_vd_index_up)
                     CALL cp_dbcsr_release(matrix_k_central)
                     CALL cp_dbcsr_release(vr_fixed)
                     CALL cp_dbcsr_release(vd_fixed)
                     CALL cp_dbcsr_release(grad)
                     CALL cp_dbcsr_release(prec)
                     CALL cp_dbcsr_release(prev_grad)
                     CALL cp_dbcsr_release(tmp3_vd_vr)
                     CALL cp_dbcsr_release(tmp1_n_vr)
                     CALL cp_dbcsr_release(tmp_k_blk)
                     CALL cp_dbcsr_release(t_curr)
                     CALL cp_dbcsr_release(sigma_oo_curr)
                     CALL cp_dbcsr_release(sigma_oo_curr_inv)
                     CALL cp_dbcsr_release(step)
                     CALL cp_dbcsr_release(tmp2_n_o)
                     CALL cp_dbcsr_release(tmp4_o_vr)
                     CALL cp_dbcsr_release(prev_step)
                     CALL cp_dbcsr_release(prev_minus_prec_grad)
                     IF (md_in_k_space) THEN
                        CALL cp_dbcsr_release(velocity)
                     ENDIF

                  ENDIF

                  outer_opt_k_iteration = outer_opt_k_iteration+1
                  IF (outer_opt_k_prepare_to_exit) EXIT

               ENDDO ! outer loop for k

            ENDDO ! ispin

            ! RZK-warning update mo orbitals

         ELSE ! virtual orbitals might not be available use projected AOs

            ! compute sqrt(S) and inv(sqrt(S))
            ! RZK-warning - remove this sqrt(S) and inv(sqrt(S))
            ! ideally ALMO scf should use sigma and sigma_inv in
            ! the tensor_up_down representation
            IF (.NOT. almo_scf_env%s_sqrt_done) THEN

               IF (unit_nr > 0) THEN
                  WRITE (unit_nr, *) "sqrt and inv(sqrt) of AO overlap matrix"
               ENDIF
               CALL cp_dbcsr_init(almo_scf_env%matrix_s_sqrt(1))
               CALL cp_dbcsr_init(almo_scf_env%matrix_s_sqrt_inv(1))
               CALL cp_dbcsr_create(almo_scf_env%matrix_s_sqrt(1), &
                                    template=almo_scf_env%matrix_s(1), &
                                    matrix_type=dbcsr_type_no_symmetry)
               CALL cp_dbcsr_create(almo_scf_env%matrix_s_sqrt_inv(1), &
                                    template=almo_scf_env%matrix_s(1), &
                                    matrix_type=dbcsr_type_no_symmetry)

               CALL matrix_sqrt_Newton_Schulz(almo_scf_env%matrix_s_sqrt(1), &
                                              almo_scf_env%matrix_s_sqrt_inv(1), &
                                              almo_scf_env%matrix_s(1), &
                                              threshold=almo_scf_env%eps_filter, &
                                              order=almo_scf_env%order_lanczos, &
                                              eps_lanczos=almo_scf_env%eps_lanczos, &
                                              max_iter_lanczos=almo_scf_env%max_iter_lanczos)

               IF (safe_mode) THEN
                  CALL cp_dbcsr_init(matrix_tmp1)
                  CALL cp_dbcsr_create(matrix_tmp1, template=almo_scf_env%matrix_s(1), &
                                       matrix_type=dbcsr_type_no_symmetry)
                  CALL cp_dbcsr_init(matrix_tmp2)
                  CALL cp_dbcsr_create(matrix_tmp2, template=almo_scf_env%matrix_s(1), &
                                       matrix_type=dbcsr_type_no_symmetry)

                  CALL cp_dbcsr_multiply("N", "N", 1.0_dp, almo_scf_env%matrix_s_sqrt_inv(1), &
                                         almo_scf_env%matrix_s(1), &
                                         0.0_dp, matrix_tmp1, filter_eps=almo_scf_env%eps_filter)
                  CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1, almo_scf_env%matrix_s_sqrt_inv(1), &
                                         0.0_dp, matrix_tmp2, filter_eps=almo_scf_env%eps_filter)

                  frob_matrix_base = cp_dbcsr_frobenius_norm(matrix_tmp2)
                  CALL cp_dbcsr_add_on_diag(matrix_tmp2, -1.0_dp)
                  frob_matrix = cp_dbcsr_frobenius_norm(matrix_tmp2)
                  IF (unit_nr > 0) THEN
                     WRITE (unit_nr, *) "Error for (inv(sqrt(S))*S*inv(sqrt(S))-I)", frob_matrix/frob_matrix_base
                  ENDIF

                  CALL cp_dbcsr_release(matrix_tmp1)
                  CALL cp_dbcsr_release(matrix_tmp2)
               ENDIF

               almo_scf_env%s_sqrt_done = .TRUE.

            ENDIF

            DO ispin = 1, nspin

               CALL ct_step_env_init(ct_step_env)
               CALL ct_step_env_set(ct_step_env, &
                                    para_env=almo_scf_env%para_env, &
                                    blacs_env=almo_scf_env%blacs_env, &
                                    use_occ_orbs=.TRUE., &
                                    use_virt_orbs=almo_scf_env%deloc_cayley_use_virt_orbs, &
                                    occ_orbs_orthogonal=.FALSE., &
                                    virt_orbs_orthogonal=almo_scf_env%orthogonal_basis, &
                                    tensor_type=almo_scf_env%deloc_cayley_tensor_type, &
                                    neglect_quadratic_term=almo_scf_env%deloc_cayley_linear, &
                                    calculate_energy_corr=.TRUE., &
                                    update_p=.TRUE., &
                                    update_q=.FALSE., &
                                    pp_preconditioner_full=almo_scf_env%deloc_cayley_occ_precond, &
                                    qq_preconditioner_full=almo_scf_env%deloc_cayley_vir_precond, &
                                    eps_convergence=almo_scf_env%deloc_cayley_eps_convergence, &
                                    eps_filter=almo_scf_env%eps_filter, &
                                    !nspins=almo_scf_env%nspins,&
                                    q_index_up=almo_scf_env%matrix_s_sqrt_inv(1), &
                                    q_index_down=almo_scf_env%matrix_s_sqrt(1), &
                                    p_index_up=almo_scf_env%matrix_sigma_sqrt_inv(ispin), &
                                    p_index_down=almo_scf_env%matrix_sigma_sqrt(ispin), &
                                    matrix_ks=almo_scf_env%matrix_ks_0deloc(ispin), &
                                    matrix_p=almo_scf_env%matrix_p(ispin), &
                                    matrix_qp_template=almo_scf_env%matrix_t(ispin), &
                                    matrix_pq_template=almo_scf_env%matrix_t_tr(ispin), &
                                    matrix_t=almo_scf_env%matrix_t(ispin), &
                                    conjugator=almo_scf_env%deloc_cayley_conjugator, &
                                    max_iter=almo_scf_env%deloc_cayley_max_iter)

               ! perform calculations
               CALL ct_step_execute(ct_step_env)

               ! for now we do not need the new set of orbitals
               ! just get the energy correction
               CALL ct_step_env_get(ct_step_env, &
                                    energy_correction=energy_correction(ispin))
               !copy_da_energy_matrix=matrix_eda(ispin),&
               !copy_da_charge_matrix=matrix_cta(ispin),&

               CALL ct_step_env_clean(ct_step_env)

            ENDDO

            energy_correction(1) = energy_correction(1)*spin_factor

         ENDIF

         ! print the energy correction and exit
         DO ispin = 1, nspin

            IF (unit_nr > 0) THEN
               WRITE (unit_nr, *)
               WRITE (unit_nr, '(T2,A,I6,F20.9)') "ECORR", ispin, &
                  energy_correction(ispin)
               WRITE (unit_nr, *)
            ENDIF
            energy_correction_final = energy_correction_final+energy_correction(ispin)

            !!! print out the results of decomposition analysis
            !!IF (unit_nr>0) THEN
            !!   WRITE(unit_nr,*)
            !!   WRITE(unit_nr,'(T2,A)') "ENERGY DECOMPOSITION"
            !!ENDIF
            !!CALL cp_dbcsr_print_block_sum(eda_matrix(ispin))
            !!IF (unit_nr>0) THEN
            !!   WRITE(unit_nr,*)
            !!   WRITE(unit_nr,'(T2,A)') "CHARGE DECOMPOSITION"
            !!ENDIF
            !!CALL cp_dbcsr_print_block_sum(cta_matrix(ispin))

            ! obtain density matrix from updated MOs
            ! RZK-later sigma and sigma_inv are lost here
            CALL almo_scf_t_to_p(t=almo_scf_env%matrix_t(ispin), &
                                 p=almo_scf_env%matrix_p(ispin), &
                                 eps_filter=almo_scf_env%eps_filter, &
                                 orthog_orbs=.FALSE., &
                                 s=almo_scf_env%matrix_s(1), &
                                 sigma=almo_scf_env%matrix_sigma(ispin), &
                                 sigma_inv=almo_scf_env%matrix_sigma_inv(ispin))

            IF (almo_scf_env%nspins == 1) &
               CALL cp_dbcsr_scale(almo_scf_env%matrix_p(ispin), &
                                   spin_factor)

         ENDDO

      CASE (dm_ls_step)

         ! compute the inverse of S
         IF (.NOT. almo_scf_env%s_inv_done) THEN
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, *) "Inverting AO overlap matrix"
            ENDIF
            CALL cp_dbcsr_init(almo_scf_env%matrix_s_inv(1))
            CALL cp_dbcsr_create(almo_scf_env%matrix_s_inv(1), &
                                 template=almo_scf_env%matrix_s(1), &
                                 matrix_type=dbcsr_type_no_symmetry)
            IF (.NOT. almo_scf_env%s_sqrt_done) THEN
               CALL invert_Hotelling(almo_scf_env%matrix_s_inv(1), &
                                     almo_scf_env%matrix_s(1), &
                                     threshold=almo_scf_env%eps_filter)
            ELSE
               CALL cp_dbcsr_multiply("N", "N", 1.0_dp, almo_scf_env%matrix_s_sqrt_inv(1), &
                                      almo_scf_env%matrix_s_sqrt_inv(1), &
                                      0.0_dp, almo_scf_env%matrix_s_inv(1), &
                                      filter_eps=almo_scf_env%eps_filter)
            ENDIF

            IF (safe_mode) THEN
               CALL cp_dbcsr_init(matrix_tmp1)
               CALL cp_dbcsr_create(matrix_tmp1, template=almo_scf_env%matrix_s(1), &
                                    matrix_type=dbcsr_type_no_symmetry)
               CALL cp_dbcsr_multiply("N", "N", 1.0_dp, almo_scf_env%matrix_s_inv(1), &
                                      almo_scf_env%matrix_s(1), &
                                      0.0_dp, matrix_tmp1, &
                                      filter_eps=almo_scf_env%eps_filter)
               frob_matrix_base = cp_dbcsr_frobenius_norm(matrix_tmp1)
               CALL cp_dbcsr_add_on_diag(matrix_tmp1, -1.0_dp)
               frob_matrix = cp_dbcsr_frobenius_norm(matrix_tmp1)
               IF (unit_nr > 0) THEN
                  WRITE (unit_nr, *) "Error for (inv(S)*S-I)", &
                     frob_matrix/frob_matrix_base
               ENDIF
               CALL cp_dbcsr_release(matrix_tmp1)
            ENDIF

            almo_scf_env%s_inv_done = .TRUE.

         ENDIF

         DO ispin = 1, nspin
            ! RZK-warning the preconditioner is very important
            !       IF (.FALSE.) THEN
            !           CALL apply_matrix_preconditioner(almo_scf_env%matrix_ks(ispin),&
            !                   "forward",almo_scf_env%matrix_s_blk_sqrt(1),&
            !                   almo_scf_env%matrix_s_blk_sqrt_inv(1))
            !       ENDIF
            !CALL cp_dbcsr_filter(almo_scf_env%matrix_ks(ispin),&
            !         almo_scf_env%eps_filter)
         ENDDO

         ALLOCATE (matrix_p_almo_scf_converged(nspin))
         DO ispin = 1, nspin
            CALL cp_dbcsr_init(matrix_p_almo_scf_converged(ispin))
            CALL cp_dbcsr_create(matrix_p_almo_scf_converged(ispin), &
                                 template=almo_scf_env%matrix_p(ispin))
            CALL cp_dbcsr_copy(matrix_p_almo_scf_converged(ispin), &
                               almo_scf_env%matrix_p(ispin))
         ENDDO

         ! update the density matrix
         DO ispin = 1, nspin

            nelectron_spin_real(1) = almo_scf_env%nelectrons_spin(ispin)
            IF (almo_scf_env%nspins == 1) &
               nelectron_spin_real(1) = nelectron_spin_real(1)/2

            local_mu(1) = SUM(almo_scf_env%mu_of_domain(:, ispin))/almo_scf_env%ndomains
            fake(1) = 123523

            ! RZK UPDATE! the update algorithm is removed because
            ! RZK UPDATE! it requires updating core LS_SCF routines
            ! RZK UPDATE! (the code exists in the CVS version)
            CPABORT("CVS only: density_matrix_sign has not been updated in SVN")
            ! RZK UPDATE!CALL density_matrix_sign(almo_scf_env%matrix_p(ispin),&
            ! RZK UPDATE!                     local_mu,&
            ! RZK UPDATE!                     almo_scf_env%fixed_mu,&
            ! RZK UPDATE!                     almo_scf_env%matrix_ks_0deloc(ispin),&
            ! RZK UPDATE!                     almo_scf_env%matrix_s(1), &
            ! RZK UPDATE!                     almo_scf_env%matrix_s_inv(1), &
            ! RZK UPDATE!                     nelectron_spin_real,&
            ! RZK UPDATE!                     almo_scf_env%eps_filter,&
            ! RZK UPDATE!                     fake)
            ! RZK UPDATE!
            almo_scf_env%mu = local_mu(1)

            !IF (almo_scf_env%has_s_preconditioner) THEN
            !    CALL apply_matrix_preconditioner(&
            !             almo_scf_env%matrix_p_blk(ispin),&
            !             "forward",almo_scf_env%matrix_s_blk_sqrt(1),&
            !             almo_scf_env%matrix_s_blk_sqrt_inv(1))
            !ENDIF
            !CALL cp_dbcsr_filter(almo_scf_env%matrix_p(ispin),&
            !        almo_scf_env%eps_filter)

            IF (almo_scf_env%nspins == 1) &
               CALL cp_dbcsr_scale(almo_scf_env%matrix_p(ispin), &
                                   spin_factor)

            !CALL cp_dbcsr_trace(almo_scf_env%matrix_ks_0deloc(ispin),&
            !  almo_scf_env%matrix_p(ispin),&
            !  energy_correction(ispin))
            !IF (unit_nr>0) THEN
            !   WRITE(unit_nr,*)
            !   WRITE(unit_nr,'(T2,A,I6,F20.9)') "EFAKE",ispin,&
            !           energy_correction(ispin)
            !   WRITE(unit_nr,*)
            !ENDIF
            CALL cp_dbcsr_add(matrix_p_almo_scf_converged(ispin), &
                              almo_scf_env%matrix_p(ispin), -1.0_dp, 1.0_dp)
            CALL cp_dbcsr_trace(almo_scf_env%matrix_ks_0deloc(ispin), &
                                matrix_p_almo_scf_converged(ispin), &
                                energy_correction(ispin))

            energy_correction_final = energy_correction_final+energy_correction(ispin)

            IF (unit_nr > 0) THEN
               WRITE (unit_nr, *)
               WRITE (unit_nr, '(T2,A,I6,F20.9)') "ECORR", ispin, &
                  energy_correction(ispin)
               WRITE (unit_nr, *)
            ENDIF

         ENDDO

         DO ispin = 1, nspin
            CALL cp_dbcsr_release(matrix_p_almo_scf_converged(ispin))
         ENDDO
         DEALLOCATE (matrix_p_almo_scf_converged)

      END SELECT ! algorithm selection

      t2 = m_walltime()

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, *)
         WRITE (unit_nr, '(T2,A,F18.9,F18.9,F18.9,F12.6)') "ETOT", &
            almo_scf_env%almo_scf_energy, &
            energy_correction_final, &
            almo_scf_env%almo_scf_energy+energy_correction_final, &
            t2-t1
         WRITE (unit_nr, *)
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE harris_foulkes_correction

! **************************************************************************************************
!> \brief Computes a diagonal preconditioner for the cg optimization of k matrix
!> \param prec ...
!> \param vd_prop ...
!> \param f ...
!> \param x ...
!> \param oo_inv_x_tr ...
!> \param s ...
!> \param grad ...
!> \param vd_blk ...
!> \param t ...
!> \param template_vd_vd_blk ...
!> \param template_vr_vr_blk ...
!> \param template_n_vr ...
!> \param spin_factor ...
!> \param eps_filter ...
!> \par History
!>       2011.09 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE opt_k_create_preconditioner(prec, vd_prop, f, x, oo_inv_x_tr, s, grad, &
                                          vd_blk, t, template_vd_vd_blk, template_vr_vr_blk, template_n_vr, &
                                          spin_factor, eps_filter)

      TYPE(cp_dbcsr_type), INTENT(INOUT)                 :: prec
      TYPE(cp_dbcsr_type), INTENT(IN)                    :: vd_prop, f, x, oo_inv_x_tr, s, grad, &
                                                            vd_blk, t, template_vd_vd_blk, &
                                                            template_vr_vr_blk, template_n_vr
      REAL(KIND=dp), INTENT(IN)                          :: spin_factor, eps_filter

      CHARACTER(len=*), PARAMETER :: routineN = 'opt_k_create_preconditioner', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, p_nrows, q_nrows
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: p_diagonal, q_diagonal
      TYPE(cp_dbcsr_type)                                :: pp_diag, qq_diag, t1, t2, tmp, &
                                                            tmp1_n_vr, tmp2_n_vr, tmp_n_vd, &
                                                            tmp_vd_vd_blk, tmp_vr_vr_blk

! init diag blocks outside
! init diag blocks otside
!INTEGER                                  :: iblock_row, iblock_col,&
!                                            nblkrows_tot, nblkcols_tot
!REAL(KIND=dp), DIMENSION(:, :), POINTER  :: p_new_block
!INTEGER                                  :: mynode, hold, row, col

      CALL timeset(routineN, handle)

      ! initialize a matrix to 1.0
      CALL cp_dbcsr_init(tmp)
      CALL cp_dbcsr_create(tmp, template=prec)
      ! in order to use cp_dbcsr_set matrix blocks must exist
      CALL cp_dbcsr_copy(tmp, prec)
      CALL cp_dbcsr_set(tmp, 1.0_dp)

      ! compute qq = (Vd^tr)*F*Vd
      CALL cp_dbcsr_init(tmp_n_vd)
      CALL cp_dbcsr_create(tmp_n_vd, template=vd_prop)
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, f, vd_prop, &
                             0.0_dp, tmp_n_vd, filter_eps=eps_filter)
      CALL cp_dbcsr_init(tmp_vd_vd_blk)
      CALL cp_dbcsr_create(tmp_vd_vd_blk, &
                           template=template_vd_vd_blk)
      CALL cp_dbcsr_copy(tmp_vd_vd_blk, template_vd_vd_blk)
      CALL cp_dbcsr_multiply("T", "N", 1.0_dp, vd_prop, tmp_n_vd, &
                             0.0_dp, tmp_vd_vd_blk, &
                             retain_sparsity=.TRUE., &
                             filter_eps=eps_filter)
      ! copy diagonal elements of the result into rows of a matrix
      CALL cp_dbcsr_get_info(tmp_vd_vd_blk, nfullrows_total=q_nrows)
      ALLOCATE (q_diagonal(q_nrows))
      CALL cp_dbcsr_get_diag(tmp_vd_vd_blk, q_diagonal)
      CALL cp_dbcsr_init(qq_diag)
      CALL cp_dbcsr_create(qq_diag, &
                           template=template_vd_vd_blk)
      CALL cp_dbcsr_add_on_diag(qq_diag, 1.0_dp)
      CALL cp_dbcsr_set_diag(qq_diag, q_diagonal)
      CALL cp_dbcsr_init(t1)
      CALL cp_dbcsr_create(t1, template=prec)
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, qq_diag, tmp, &
                             0.0_dp, t1, filter_eps=eps_filter)

      ! compute pp = X*sigma_oo_inv*X^tr
      CALL cp_dbcsr_init(tmp_vr_vr_blk)
      CALL cp_dbcsr_create(tmp_vr_vr_blk, template=template_vr_vr_blk)
      CALL cp_dbcsr_copy(tmp_vr_vr_blk, template_vr_vr_blk)
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, x, oo_inv_x_tr, &
                             0.0_dp, tmp_vr_vr_blk, &
                             retain_sparsity=.TRUE., &
                             filter_eps=eps_filter)
      ! copy diagonal elements of the result into cols of a matrix
      CALL cp_dbcsr_get_info(tmp_vr_vr_blk, nfullrows_total=p_nrows)
      ALLOCATE (p_diagonal(p_nrows))
      CALL cp_dbcsr_get_diag(tmp_vr_vr_blk, p_diagonal)
      CALL cp_dbcsr_init(pp_diag)
      CALL cp_dbcsr_create(pp_diag, template=template_vr_vr_blk)
      CALL cp_dbcsr_add_on_diag(pp_diag, 1.0_dp)
      CALL cp_dbcsr_set_diag(pp_diag, p_diagonal)
      CALL cp_dbcsr_set(tmp, 1.0_dp)
      CALL cp_dbcsr_init(t2)
      CALL cp_dbcsr_create(t2, template=prec)
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, tmp, pp_diag, &
                             0.0_dp, t2, filter_eps=eps_filter)

      CALL cp_dbcsr_hadamard_product(t1, t2, prec)

      ! compute qq = (Vd^tr)*S*Vd
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, s, vd_prop, &
                             0.0_dp, tmp_n_vd, filter_eps=eps_filter)
      CALL cp_dbcsr_multiply("T", "N", 1.0_dp, vd_prop, tmp_n_vd, &
                             0.0_dp, tmp_vd_vd_blk, &
                             retain_sparsity=.TRUE., &
                             filter_eps=eps_filter)
      ! copy diagonal elements of the result into rows of a matrix
      CALL cp_dbcsr_get_diag(tmp_vd_vd_blk, q_diagonal)
      CALL cp_dbcsr_add_on_diag(qq_diag, 1.0_dp)
      CALL cp_dbcsr_set_diag(qq_diag, q_diagonal)
      CALL cp_dbcsr_set(tmp, 1.0_dp)
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, qq_diag, tmp, &
                             0.0_dp, t1, filter_eps=eps_filter)

      ! compute pp = X*sig_oo_inv*(T^tr)*F*T*sig_oo_inv*(X^tr)
      CALL cp_dbcsr_init(tmp1_n_vr)
      CALL cp_dbcsr_create(tmp1_n_vr, template=template_n_vr)
      CALL cp_dbcsr_init(tmp2_n_vr)
      CALL cp_dbcsr_create(tmp2_n_vr, template=template_n_vr)
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, t, oo_inv_x_tr, &
                             0.0_dp, tmp1_n_vr, filter_eps=eps_filter)
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, f, tmp1_n_vr, &
                             0.0_dp, tmp2_n_vr, filter_eps=eps_filter)
      CALL cp_dbcsr_multiply("T", "N", 1.0_dp, tmp1_n_vr, tmp2_n_vr, &
                             0.0_dp, tmp_vr_vr_blk, &
                             retain_sparsity=.TRUE., &
                             filter_eps=eps_filter)
      ! copy diagonal elements of the result into cols of a matrix
      CALL cp_dbcsr_get_diag(tmp_vr_vr_blk, p_diagonal)
      CALL cp_dbcsr_add_on_diag(pp_diag, 1.0_dp)
      CALL cp_dbcsr_set_diag(pp_diag, p_diagonal)
      CALL cp_dbcsr_set(tmp, 1.0_dp)
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, tmp, pp_diag, &
                             0.0_dp, t2, filter_eps=eps_filter)

      CALL cp_dbcsr_hadamard_product(t1, t2, tmp)
      CALL cp_dbcsr_add(prec, tmp, 1.0_dp, -1.0_dp)
      CALL cp_dbcsr_scale(prec, 2.0_dp*spin_factor)

      ! compute qp = X*sig_oo_inv*(T^tr)*S*Vd
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, s, vd_blk, &
                             0.0_dp, tmp_n_vd, filter_eps=eps_filter)
      CALL cp_dbcsr_multiply("T", "N", 1.0_dp, tmp_n_vd, tmp1_n_vr, &
                             0.0_dp, tmp, retain_sparsity=.TRUE., &
                             filter_eps=eps_filter)
      CALL cp_dbcsr_hadamard_product(grad, tmp, t1)
      ! gradient already contains 2.0*spin_factor
      CALL cp_dbcsr_scale(t1, -2.0_dp)

      CALL cp_dbcsr_add(prec, t1, 1.0_dp, 1.0_dp)

      CALL cp_dbcsr_function_of_elements(prec, dbcsr_func_inverse)
      CALL cp_dbcsr_filter(prec, eps_filter)

      DEALLOCATE (q_diagonal)
      DEALLOCATE (p_diagonal)
      CALL cp_dbcsr_release(tmp)
      CALL cp_dbcsr_release(qq_diag)
      CALL cp_dbcsr_release(t1)
      CALL cp_dbcsr_release(pp_diag)
      CALL cp_dbcsr_release(t2)
      CALL cp_dbcsr_release(tmp_n_vd)
      CALL cp_dbcsr_release(tmp_vd_vd_blk)
      CALL cp_dbcsr_release(tmp_vr_vr_blk)
      CALL cp_dbcsr_release(tmp1_n_vr)
      CALL cp_dbcsr_release(tmp2_n_vr)

      CALL timestop(handle)

   END SUBROUTINE opt_k_create_preconditioner

! **************************************************************************************************
!> \brief Computes a block-diagonal preconditioner for the optimization of
!>        k matrix
!> \param almo_scf_env ...
!> \param vd_prop ...
!> \param oo_inv_x_tr ...
!> \param t_curr ...
!> \param ispin ...
!> \param spin_factor ...
!> \par History
!>       2011.10 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE opt_k_create_preconditioner_blk(almo_scf_env, vd_prop, oo_inv_x_tr, &
                                              t_curr, ispin, spin_factor)

      TYPE(almo_scf_env_type), INTENT(INOUT)             :: almo_scf_env
      TYPE(cp_dbcsr_type), INTENT(IN)                    :: vd_prop, oo_inv_x_tr, t_curr
      INTEGER, INTENT(IN)                                :: ispin
      REAL(KIND=dp), INTENT(IN)                          :: spin_factor

      CHARACTER(len=*), PARAMETER :: routineN = 'opt_k_create_preconditioner_blk', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      REAL(KIND=dp)                                      :: eps_filter
      TYPE(cp_dbcsr_type)                                :: opt_k_e_dd, opt_k_e_rr, s_dd_sqrt, &
                                                            s_rr_sqrt, t1, tmp, tmp1_n_vr, &
                                                            tmp2_n_vr, tmp_n_vd, tmp_vd_vd_blk, &
                                                            tmp_vr_vr_blk

! matrices that has been computed outside the routine already

      CALL timeset(routineN, handle)

      eps_filter = almo_scf_env%eps_filter

      ! compute S_qq = (Vd^tr)*S*Vd
      CALL cp_dbcsr_init(tmp_n_vd)
      CALL cp_dbcsr_create(tmp_n_vd, template=almo_scf_env%matrix_v_disc(ispin))
      CALL cp_dbcsr_init(tmp_vd_vd_blk)
      CALL cp_dbcsr_create(tmp_vd_vd_blk, &
                           template=almo_scf_env%matrix_vv_disc_blk(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, &
                             almo_scf_env%matrix_s(1), &
                             vd_prop, &
                             0.0_dp, tmp_n_vd, filter_eps=eps_filter)
      CALL cp_dbcsr_copy(tmp_vd_vd_blk, &
                         almo_scf_env%matrix_vv_disc_blk(ispin))
      CALL cp_dbcsr_multiply("T", "N", 1.0_dp, vd_prop, tmp_n_vd, &
                             0.0_dp, tmp_vd_vd_blk, &
                             retain_sparsity=.TRUE.)

      CALL cp_dbcsr_init(s_dd_sqrt)
      CALL cp_dbcsr_create(s_dd_sqrt, &
                           template=almo_scf_env%matrix_vv_disc_blk(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
      CALL matrix_sqrt_Newton_Schulz(s_dd_sqrt, &
                                     almo_scf_env%opt_k_t_dd(ispin), &
                                     tmp_vd_vd_blk, &
                                     threshold=eps_filter, &
                                     order=almo_scf_env%order_lanczos, &
                                     eps_lanczos=almo_scf_env%eps_lanczos, &
                                     max_iter_lanczos=almo_scf_env%max_iter_lanczos)

      ! compute F_qq = (Vd^tr)*F*Vd
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, &
                             almo_scf_env%matrix_ks_0deloc(ispin), &
                             vd_prop, &
                             0.0_dp, tmp_n_vd, filter_eps=eps_filter)
      CALL cp_dbcsr_copy(tmp_vd_vd_blk, &
                         almo_scf_env%matrix_vv_disc_blk(ispin))
      CALL cp_dbcsr_multiply("T", "N", 1.0_dp, vd_prop, tmp_n_vd, &
                             0.0_dp, tmp_vd_vd_blk, &
                             retain_sparsity=.TRUE.)
      CALL cp_dbcsr_release(tmp_n_vd)

      ! bring to the blocked-orthogonalized basis
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, &
                             tmp_vd_vd_blk, &
                             almo_scf_env%opt_k_t_dd(ispin), &
                             0.0_dp, s_dd_sqrt, filter_eps=eps_filter)
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, &
                             almo_scf_env%opt_k_t_dd(ispin), &
                             s_dd_sqrt, &
                             0.0_dp, tmp_vd_vd_blk, filter_eps=eps_filter)

      ! diagonalize the matrix
      CALL cp_dbcsr_init(opt_k_e_dd)
      CALL cp_dbcsr_create(opt_k_e_dd, &
                           template=almo_scf_env%matrix_vv_disc_blk(ispin))
      CALL cp_dbcsr_release(s_dd_sqrt)
      CALL cp_dbcsr_init(s_dd_sqrt)
      CALL cp_dbcsr_create(s_dd_sqrt, &
                           template=almo_scf_env%matrix_vv_disc_blk(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
      CALL diagonalize_diagonal_blocks(tmp_vd_vd_blk, &
                                       s_dd_sqrt, &
                                       opt_k_e_dd)

      ! obtain the transformation matrix in the discarded subspace
      ! T = S^{-1/2}.U
      CALL cp_dbcsr_copy(tmp_vd_vd_blk, &
                         almo_scf_env%opt_k_t_dd(ispin))
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, &
                             tmp_vd_vd_blk, &
                             s_dd_sqrt, &
                             0.0_dp, almo_scf_env%opt_k_t_dd(ispin), &
                             filter_eps=eps_filter)
      CALL cp_dbcsr_release(s_dd_sqrt)
      CALL cp_dbcsr_release(tmp_vd_vd_blk)

      ! copy diagonal elements of the result into rows of a matrix
      CALL cp_dbcsr_init(tmp)
      CALL cp_dbcsr_create(tmp, &
                           template=almo_scf_env%matrix_k_blk_ones(ispin))
      CALL cp_dbcsr_copy(tmp, &
                         almo_scf_env%matrix_k_blk_ones(ispin))
      CALL cp_dbcsr_init(t1)
      CALL cp_dbcsr_create(t1, &
                           template=almo_scf_env%matrix_k_blk_ones(ispin))
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, &
                             opt_k_e_dd, tmp, &
                             0.0_dp, t1, filter_eps=eps_filter)
      CALL cp_dbcsr_release(opt_k_e_dd)

      ! compute S_pp = X*sigma_oo_inv*X^tr
      CALL cp_dbcsr_init(tmp_vr_vr_blk)
      CALL cp_dbcsr_create(tmp_vr_vr_blk, &
                           template=almo_scf_env%matrix_sigma_vv_blk(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
      CALL cp_dbcsr_copy(tmp_vr_vr_blk, &
                         almo_scf_env%matrix_sigma_vv_blk(ispin))
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, &
                             almo_scf_env%matrix_x(ispin), &
                             oo_inv_x_tr, &
                             0.0_dp, tmp_vr_vr_blk, &
                             retain_sparsity=.TRUE.)

      ! obtain the orthogonalization matrix
      CALL cp_dbcsr_init(s_rr_sqrt)
      CALL cp_dbcsr_create(s_rr_sqrt, &
                           template=almo_scf_env%matrix_sigma_vv_blk(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
      CALL matrix_sqrt_Newton_Schulz(s_rr_sqrt, &
                                     almo_scf_env%opt_k_t_rr(ispin), &
                                     tmp_vr_vr_blk, &
                                     threshold=eps_filter, &
                                     order=almo_scf_env%order_lanczos, &
                                     eps_lanczos=almo_scf_env%eps_lanczos, &
                                     max_iter_lanczos=almo_scf_env%max_iter_lanczos)

      ! compute F_pp = X*sig_oo_inv*(T^tr)*F*T*sig_oo_inv*(X^tr)
      CALL cp_dbcsr_init(tmp1_n_vr)
      CALL cp_dbcsr_create(tmp1_n_vr, &
                           template=almo_scf_env%matrix_v(ispin))
      CALL cp_dbcsr_init(tmp2_n_vr)
      CALL cp_dbcsr_create(tmp2_n_vr, &
                           template=almo_scf_env%matrix_v(ispin))
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, t_curr, oo_inv_x_tr, &
                             0.0_dp, tmp1_n_vr, filter_eps=eps_filter)
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, &
                             almo_scf_env%matrix_ks_0deloc(ispin), &
                             tmp1_n_vr, &
                             0.0_dp, tmp2_n_vr, filter_eps=eps_filter)
      CALL cp_dbcsr_multiply("T", "N", 1.0_dp, tmp1_n_vr, tmp2_n_vr, &
                             0.0_dp, tmp_vr_vr_blk, &
                             retain_sparsity=.TRUE.)
      CALL cp_dbcsr_release(tmp1_n_vr)
      CALL cp_dbcsr_release(tmp2_n_vr)

      ! bring to the blocked-orthogonalized basis
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, &
                             tmp_vr_vr_blk, &
                             almo_scf_env%opt_k_t_rr(ispin), &
                             0.0_dp, s_rr_sqrt, filter_eps=eps_filter)
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, &
                             almo_scf_env%opt_k_t_rr(ispin), &
                             s_rr_sqrt, &
                             0.0_dp, tmp_vr_vr_blk, filter_eps=eps_filter)

      ! diagonalize the matrix
      CALL cp_dbcsr_init(opt_k_e_rr)
      CALL cp_dbcsr_create(opt_k_e_rr, &
                           template=almo_scf_env%matrix_sigma_vv_blk(ispin))
      CALL cp_dbcsr_release(s_rr_sqrt)
      CALL cp_dbcsr_init(s_rr_sqrt)
      CALL cp_dbcsr_create(s_rr_sqrt, &
                           template=almo_scf_env%matrix_sigma_vv_blk(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
      CALL diagonalize_diagonal_blocks(tmp_vr_vr_blk, &
                                       s_rr_sqrt, &
                                       opt_k_e_rr)

      ! obtain the transformation matrix in the retained subspace
      ! T = S^{-1/2}.U
      CALL cp_dbcsr_copy(tmp_vr_vr_blk, &
                         almo_scf_env%opt_k_t_rr(ispin))
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, &
                             tmp_vr_vr_blk, &
                             s_rr_sqrt, &
                             0.0_dp, almo_scf_env%opt_k_t_rr(ispin), &
                             filter_eps=eps_filter)
      CALL cp_dbcsr_release(s_rr_sqrt)
      CALL cp_dbcsr_release(tmp_vr_vr_blk)

      ! copy diagonal elements of the result into cols of a matrix
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, &
                             tmp, opt_k_e_rr, &
                             0.0_dp, almo_scf_env%opt_k_denom(ispin), &
                             filter_eps=eps_filter)
      CALL cp_dbcsr_release(opt_k_e_rr)
      CALL cp_dbcsr_release(tmp)

      ! form the denominator matrix
      CALL cp_dbcsr_add(almo_scf_env%opt_k_denom(ispin), t1, &
                        -1.0_dp, 1.0_dp)
      CALL cp_dbcsr_release(t1)
      CALL cp_dbcsr_scale(almo_scf_env%opt_k_denom(ispin), &
                          2.0_dp*spin_factor)

      CALL cp_dbcsr_function_of_elements(almo_scf_env%opt_k_denom(ispin), &
                                         dbcsr_func_inverse)
      CALL cp_dbcsr_filter(almo_scf_env%opt_k_denom(ispin), &
                           eps_filter)

      CALL timestop(handle)

   END SUBROUTINE opt_k_create_preconditioner_blk

! **************************************************************************************************
!> \brief Applies a block-diagonal preconditioner for the optimization of
!>        k matrix (preconditioner matrices must be calculated and stored
!>        beforehand)
!> \param almo_scf_env ...
!> \param step ...
!> \param grad ...
!> \param ispin ...
!> \par History
!>       2011.10 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE opt_k_apply_preconditioner_blk(almo_scf_env, step, grad, ispin)

      TYPE(almo_scf_env_type), INTENT(INOUT)             :: almo_scf_env
      TYPE(cp_dbcsr_type), INTENT(OUT)                   :: step
      TYPE(cp_dbcsr_type), INTENT(IN)                    :: grad
      INTEGER, INTENT(IN)                                :: ispin

      CHARACTER(len=*), PARAMETER :: routineN = 'opt_k_apply_preconditioner_blk', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      REAL(KIND=dp)                                      :: eps_filter
      TYPE(cp_dbcsr_type)                                :: tmp_k

      CALL timeset(routineN, handle)

      eps_filter = almo_scf_env%eps_filter

      CALL cp_dbcsr_init(tmp_k)
      CALL cp_dbcsr_create(tmp_k, template=almo_scf_env%matrix_k_blk(ispin))

      ! transform gradient to the correct "diagonal" basis
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, &
                             grad, almo_scf_env%opt_k_t_rr(ispin), &
                             0.0_dp, tmp_k, filter_eps=eps_filter)
      CALL cp_dbcsr_multiply("T", "N", 1.0_dp, &
                             almo_scf_env%opt_k_t_dd(ispin), tmp_k, &
                             0.0_dp, step, filter_eps=eps_filter)

      ! apply diagonal preconditioner
      CALL cp_dbcsr_hadamard_product(step, &
                                     almo_scf_env%opt_k_denom(ispin), tmp_k)

      ! back-transform the result to the initial basis
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, &
                             almo_scf_env%opt_k_t_dd(ispin), tmp_k, &
                             0.0_dp, step, filter_eps=eps_filter)
      CALL cp_dbcsr_multiply("N", "T", 1.0_dp, &
                             step, almo_scf_env%opt_k_t_rr(ispin), &
                             0.0_dp, tmp_k, filter_eps=eps_filter)

      CALL cp_dbcsr_copy(step, tmp_k)

      CALL cp_dbcsr_release(tmp_k)

      CALL timestop(handle)

   END SUBROUTINE opt_k_apply_preconditioner_blk

!! **************************************************************************************************
!!> \brief Reduce the number of virtual orbitals by rotating them within
!!>        a domain. The rotation is such that minimizes the frobenius norm of
!!>        the Fov domain-blocks of the discarded virtuals
!!> \par History
!!>       2011.08 created [Rustam Z Khaliullin]
!!> \author Rustam Z Khaliullin
!! **************************************************************************************************
!  SUBROUTINE truncate_subspace_v_blk(qs_env,almo_scf_env)
!
!    TYPE(qs_environment_type), POINTER       :: qs_env
!    TYPE(almo_scf_env_type)                  :: almo_scf_env
!
!    CHARACTER(len=*), PARAMETER :: routineN = 'truncate_subspace_v_blk', &
!      routineP = moduleN//':'//routineN
!
!    INTEGER                                  :: handle, ispin, iblock_row, &
!                                                iblock_col, iblock_row_size, &
!                                                iblock_col_size, retained_v, &
!                                                iteration, line_search_step, &
!                                                unit_nr, line_search_step_last
!    REAL(KIND=dp)                            :: t1, obj_function, grad_norm,&
!                                                c0, b0, a0, obj_function_new,&
!                                                t2, alpha, ff1, ff2, step1,&
!                                                step2,&
!                                                frob_matrix_base,&
!                                                frob_matrix
!    LOGICAL                                  :: safe_mode, converged, &
!                                                prepare_to_exit, failure
!    TYPE(cp_logger_type), POINTER            :: logger
!    TYPE(cp_dbcsr_type)                      :: Fon, Fov, Fov_filtered, &
!                                                temp1_oo, temp2_oo, Fov_original, &
!                                                temp0_ov, U_blk_tot, U_blk, &
!                                                grad_blk, step_blk, matrix_filter, &
!                                                v_full_new,v_full_tmp,&
!                                                matrix_sigma_vv_full,&
!                                                matrix_sigma_vv_full_sqrt,&
!                                                matrix_sigma_vv_full_sqrt_inv,&
!                                                matrix_tmp1,&
!                                                matrix_tmp2
!
!    REAL(kind=dp), DIMENSION(:, :), POINTER  :: data_p, p_new_block
!    TYPE(cp_dbcsr_iterator)                  :: iter
!
!
!REAL(kind=dp), DIMENSION(:), ALLOCATABLE     :: eigenvalues, WORK
!REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE   :: data_copy, left_vectors, right_vectors
!INTEGER                                      :: LWORK, INFO
!TYPE(cp_dbcsr_type)                          :: temp_u_v_full_blk
!
!    CALL timeset(routineN,handle)
!
!    safe_mode=.TRUE.
!
!    ! get a useful output_unit
!    logger => cp_get_default_logger()
!    IF (logger%para_env%mepos==logger%para_env%source) THEN
!       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
!    ELSE
!       unit_nr=-1
!    ENDIF
!
!    DO ispin=1,almo_scf_env%nspins
!
!       t1 = m_walltime()
!
!       !!!!!!!!!!!!!!!!!
!       ! 0. Orthogonalize virtuals
!       !    Unfortunately, we have to do it in the FULL V subspace :(
!
!       CALL cp_dbcsr_init(v_full_new)
!       CALL cp_dbcsr_create(v_full_new,&
!               template=almo_scf_env%matrix_v_full_blk(ispin),&
!               matrix_type=dbcsr_type_no_symmetry)
!
!       ! project the occupied subspace out
!       CALL almo_scf_p_out_from_v(almo_scf_env%matrix_v_full_blk(ispin),&
!              v_full_new,almo_scf_env%matrix_ov_full(ispin),&
!              ispin,almo_scf_env)
!
!       ! init overlap and its functions
!       CALL cp_dbcsr_init(matrix_sigma_vv_full)
!       CALL cp_dbcsr_init(matrix_sigma_vv_full_sqrt)
!       CALL cp_dbcsr_init(matrix_sigma_vv_full_sqrt_inv)
!       CALL cp_dbcsr_create(matrix_sigma_vv_full,&
!               template=almo_scf_env%matrix_vv_full_blk(ispin),&
!               matrix_type=dbcsr_type_no_symmetry)
!       CALL cp_dbcsr_create(matrix_sigma_vv_full_sqrt,&
!               template=almo_scf_env%matrix_vv_full_blk(ispin),&
!               matrix_type=dbcsr_type_no_symmetry)
!       CALL cp_dbcsr_create(matrix_sigma_vv_full_sqrt_inv,&
!               template=almo_scf_env%matrix_vv_full_blk(ispin),&
!               matrix_type=dbcsr_type_no_symmetry)
!
!       ! construct VV overlap
!       CALL almo_scf_mo_to_sigma(v_full_new,&
!               matrix_sigma_vv_full,&
!               almo_scf_env%matrix_s(1),&
!               almo_scf_env%eps_filter)
!
!       IF (unit_nr>0) THEN
!          WRITE(unit_nr,*) "sqrt and inv(sqrt) of the FULL virtual MO overlap"
!       ENDIF
!
!       ! construct orthogonalization matrices
!       CALL matrix_sqrt_Newton_Schulz(matrix_sigma_vv_full_sqrt,&
!                                      matrix_sigma_vv_full_sqrt_inv,&
!                                      matrix_sigma_vv_full,&
!                                      threshold=almo_scf_env%eps_filter,&
!                                      order=almo_scf_env%order_lanczos,&
!                                      eps_lanczos=almo_scf_env%eps_lanczos,&
!                                      max_iter_lanczos=almo_scf_env%max_iter_lanczos)
!       IF (safe_mode) THEN
!          CALL cp_dbcsr_init(matrix_tmp1)
!          CALL cp_dbcsr_create(matrix_tmp1,template=matrix_sigma_vv_full,&
!                               matrix_type=dbcsr_type_no_symmetry)
!          CALL cp_dbcsr_init(matrix_tmp2)
!          CALL cp_dbcsr_create(matrix_tmp2,template=matrix_sigma_vv_full,&
!                               matrix_type=dbcsr_type_no_symmetry)
!
!          CALL cp_dbcsr_multiply("N","N",1.0_dp,matrix_sigma_vv_full_sqrt_inv,&
!                                 matrix_sigma_vv_full,&
!                                 0.0_dp,matrix_tmp1,filter_eps=almo_scf_env%eps_filter)
!          CALL cp_dbcsr_multiply("N","N",1.0_dp,matrix_tmp1,&
!                                 matrix_sigma_vv_full_sqrt_inv,&
!                                 0.0_dp,matrix_tmp2,filter_eps=almo_scf_env%eps_filter)
!
!          frob_matrix_base=cp_dbcsr_frobenius_norm(matrix_tmp2)
!          CALL cp_dbcsr_add_on_diag(matrix_tmp2,-1.0_dp)
!          frob_matrix=cp_dbcsr_frobenius_norm(matrix_tmp2)
!          IF (unit_nr>0) THEN
!             WRITE(unit_nr,*) "Error for (inv(sqrt(SIGVV))*SIGVV*inv(sqrt(SIGVV))-I)",frob_matrix/frob_matrix_base
!          ENDIF
!
!          CALL cp_dbcsr_release(matrix_tmp1)
!          CALL cp_dbcsr_release(matrix_tmp2)
!       ENDIF
!
!       ! discard unnecessary overlap functions
!       CALL cp_dbcsr_release(matrix_sigma_vv_full)
!       CALL cp_dbcsr_release(matrix_sigma_vv_full_sqrt)
!
!! this can be re-written because we have (1-P)|v>
!
!       !!!!!!!!!!!!!!!!!!!
!       ! 1. Compute F_ov
!       CALL cp_dbcsr_init(Fon)
!       CALL cp_dbcsr_create(Fon,&
!               template=almo_scf_env%matrix_v_full_blk(ispin))
!       CALL cp_dbcsr_init(Fov)
!       CALL cp_dbcsr_create(Fov,&
!               template=almo_scf_env%matrix_ov_full(ispin))
!       CALL cp_dbcsr_init(Fov_filtered)
!       CALL cp_dbcsr_create(Fov_filtered,&
!               template=almo_scf_env%matrix_ov_full(ispin))
!       CALL cp_dbcsr_init(temp1_oo)
!       CALL cp_dbcsr_create(temp1_oo,&
!               template=almo_scf_env%matrix_sigma(ispin),&
!               !matrix_type=dbcsr_type_no_symmetry)
!       CALL cp_dbcsr_init(temp2_oo)
!       CALL cp_dbcsr_create(temp2_oo,&
!               template=almo_scf_env%matrix_sigma(ispin),&
!               matrix_type=dbcsr_type_no_symmetry)
!
!       CALL cp_dbcsr_multiply("T","N",1.0_dp,almo_scf_env%matrix_t_blk(ispin),&
!               almo_scf_env%matrix_ks_0deloc(ispin),&
!               0.0_dp,Fon,filter_eps=almo_scf_env%eps_filter)
!
!       CALL cp_dbcsr_multiply("N","N",1.0_dp,Fon,&
!               almo_scf_env%matrix_v_full_blk(ispin),&
!               0.0_dp,Fov,filter_eps=almo_scf_env%eps_filter)
!
!       CALL cp_dbcsr_multiply("N","N",1.0_dp,Fon,&
!               almo_scf_env%matrix_t_blk(ispin),&
!               0.0_dp,temp1_oo,filter_eps=almo_scf_env%eps_filter)
!
!       CALL cp_dbcsr_multiply("N","N",1.0_dp,temp1_oo,&
!               almo_scf_env%matrix_sigma_inv(ispin),&
!               0.0_dp,temp2_oo,filter_eps=almo_scf_env%eps_filter)
!       CALL cp_dbcsr_release(temp1_oo)
!
!       CALL cp_dbcsr_multiply("T","N",1.0_dp,almo_scf_env%matrix_t_blk(ispin),&
!               almo_scf_env%matrix_s(1),&
!               0.0_dp,Fon,filter_eps=almo_scf_env%eps_filter)
!
!       CALL cp_dbcsr_multiply("N","N",1.0_dp,Fon,&
!               almo_scf_env%matrix_v_full_blk(ispin),&
!               0.0_dp,Fov_filtered,filter_eps=almo_scf_env%eps_filter)
!       CALL cp_dbcsr_release(Fon)
!
!       CALL cp_dbcsr_multiply("N","N",-1.0_dp,temp2_oo,&
!               Fov_filtered,&
!               1.0_dp,Fov,filter_eps=almo_scf_env%eps_filter)
!       CALL cp_dbcsr_release(temp2_oo)
!
!       CALL cp_dbcsr_multiply("N","N",1.0_dp,almo_scf_env%matrix_sigma_inv(ispin),&
!               Fov,0.0_dp,Fov_filtered,filter_eps=almo_scf_env%eps_filter)
!
!       CALL cp_dbcsr_multiply("N","N",1.0_dp,Fov_filtered,&
!               matrix_sigma_vv_full_sqrt_inv,&
!               0.0_dp,Fov,filter_eps=almo_scf_env%eps_filter)
!       !CALL cp_dbcsr_copy(Fov,Fov_filtered)
!CALL cp_dbcsr_print(Fov)
!
!       IF (safe_mode) THEN
!          CALL cp_dbcsr_init(Fov_original)
!          CALL cp_dbcsr_create(Fov_original,template=Fov)
!          CALL cp_dbcsr_copy(Fov_original,Fov)
!       ENDIF
!
!!! remove diagonal blocks
!!CALL cp_dbcsr_iterator_start(iter,Fov)
!!DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
!!
!!   CALL cp_dbcsr_iterator_next_block(iter,iblock_row,iblock_col,data_p,&
!!           row_size=iblock_row_size,col_size=iblock_col_size)
!!
!!   IF (iblock_row.eq.iblock_col) data_p(:,:)=0.0_dp
!!
!!ENDDO
!!CALL cp_dbcsr_iterator_stop(iter)
!!CALL cp_dbcsr_finalize(Fov)
!
!!! perform svd of blocks
!!!!! THIS ROUTINE WORKS ONLY ON ONE CPU AND ONLY FOR 2 MOLECULES !!!
!!CALL cp_dbcsr_init(temp_u_v_full_blk)
!!CALL cp_dbcsr_create(temp_u_v_full_blk,&
!!        template=almo_scf_env%matrix_vv_full_blk(ispin),&
!!        matrix_type=dbcsr_type_no_symmetry)
!!
!!CALL cp_dbcsr_work_create(temp_u_v_full_blk,&
!!        work_mutable=.TRUE.)
!!CALL cp_dbcsr_iterator_start(iter,Fov)
!!DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
!!
!!   CALL cp_dbcsr_iterator_next_block(iter,iblock_row,iblock_col,data_p,&
!!           row_size=iblock_row_size,col_size=iblock_col_size)
!!
!!   IF (iblock_row.ne.iblock_col) THEN
!!
!!      ! Prepare data
!!      allocate(eigenvalues(min(iblock_row_size,iblock_col_size)))
!!      allocate(data_copy(iblock_row_size,iblock_col_size))
!!      allocate(left_vectors(iblock_row_size,iblock_row_size))
!!      allocate(right_vectors(iblock_col_size,iblock_col_size))
!!      data_copy(:,:)=data_p(:,:)
!!
!!      ! Query the optimal workspace for dgesvd
!!      LWORK = -1
!!      allocate(WORK(MAX(1,LWORK)))
!!      CALL DGESVD('N','A',iblock_row_size,iblock_col_size,data_copy,&
!!              iblock_row_size,eigenvalues,left_vectors,iblock_row_size,&
!!              right_vectors,iblock_col_size,WORK,LWORK,INFO)
!!      LWORK = INT(WORK( 1 ))
!!      deallocate(WORK)
!!
!!      ! Allocate the workspace and perform svd
!!      allocate(WORK(MAX(1,LWORK)))
!!      CALL DGESVD('N','A',iblock_row_size,iblock_col_size,data_copy,&
!!              iblock_row_size,eigenvalues,left_vectors,iblock_row_size,&
!!              right_vectors,iblock_col_size,WORK,LWORK,INFO)
!!      deallocate(WORK)
!!      IF( INFO.NE.0 ) THEN
!!         CPErrorMessage(cp_failure_level,routineP,"DGESVD failed")
!!         CPPrecondition(.FALSE.,cp_failure_level,routineP,failure)
!!      END IF
!!
!!      ! copy right singular vectors into a unitary matrix
!!      NULLIFY (p_new_block)
!!      CALL cp_dbcsr_reserve_block2d(temp_u_v_full_blk,iblock_col,iblock_col,p_new_block)
!!      CPPostcondition(ASSOCIATED(p_new_block),cp_failure_level,routineP,failure)
!!      p_new_block(:,:) = right_vectors(:,:)
!!
!!      deallocate(eigenvalues)
!!      deallocate(data_copy)
!!      deallocate(left_vectors)
!!      deallocate(right_vectors)
!!
!!   ENDIF
!!ENDDO
!!CALL cp_dbcsr_iterator_stop(iter)
!!CALL cp_dbcsr_finalize(temp_u_v_full_blk)
!!!CALL cp_dbcsr_print(temp_u_v_full_blk)
!!CALL cp_dbcsr_multiply("N","T",1.0_dp,Fov,temp_u_v_full_blk,&
!!        0.0_dp,Fov_filtered,filter_eps=almo_scf_env%eps_filter)
!!
!!CALL cp_dbcsr_copy(Fov,Fov_filtered)
!!CALL cp_dbcsr_print(Fov)
!
!       !!!!!!!!!!!!!!!!!!!
!       ! 2. Initialize variables
!
!       ! temp space
!       CALL cp_dbcsr_init(temp0_ov)
!       CALL cp_dbcsr_create(temp0_ov,&
!               template=almo_scf_env%matrix_ov_full(ispin))
!
!       ! current unitary matrix
!       CALL cp_dbcsr_init(U_blk)
!       CALL cp_dbcsr_create(U_blk,&
!               template=almo_scf_env%matrix_vv_full_blk(ispin),&
!               matrix_type=dbcsr_type_no_symmetry)
!
!       ! unitary matrix accumulator
!       CALL cp_dbcsr_init(U_blk_tot)
!       CALL cp_dbcsr_create(U_blk_tot,&
!               template=almo_scf_env%matrix_vv_full_blk(ispin),&
!               matrix_type=dbcsr_type_no_symmetry)
!       CALL cp_dbcsr_add_on_diag(U_blk_tot,1.0_dp)
!
!!CALL cp_dbcsr_add_on_diag(U_blk,1.0_dp)
!!CALL cp_dbcsr_multiply("N","T",1.0_dp,U_blk,temp_u_v_full_blk,&
!!        0.0_dp,U_blk_tot,filter_eps=almo_scf_env%eps_filter)
!!
!!CALL cp_dbcsr_release(temp_u_v_full_blk)
!
!       ! init gradient
!       CALL cp_dbcsr_init(grad_blk)
!       CALL cp_dbcsr_create(grad_blk,&
!               template=almo_scf_env%matrix_vv_full_blk(ispin),&
!               matrix_type=dbcsr_type_no_symmetry)
!
!       ! init step matrix
!       CALL cp_dbcsr_init(step_blk)
!       CALL cp_dbcsr_create(step_blk,&
!               template=almo_scf_env%matrix_vv_full_blk(ispin),&
!               matrix_type=dbcsr_type_no_symmetry)
!
!       ! "retain discarded" filter (0.0 - retain, 1.0 - discard)
!       CALL cp_dbcsr_init(matrix_filter)
!       CALL cp_dbcsr_create(matrix_filter,&
!               template=almo_scf_env%matrix_ov_full(ispin))
!       ! copy Fov into the filter matrix temporarily
!       ! so we know which blocks contain significant elements
!       CALL cp_dbcsr_copy(matrix_filter,Fov)
!
!       ! fill out filter elements block-by-block
!       CALL cp_dbcsr_iterator_start(iter,matrix_filter)
!       DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
!
!          CALL cp_dbcsr_iterator_next_block(iter,iblock_row,iblock_col,data_p,&
!                  row_size=iblock_row_size,col_size=iblock_col_size)
!
!          retained_v=almo_scf_env%nvirt_of_domain(iblock_col,ispin)
!
!          data_p(:,1:retained_v)=0.0_dp
!          data_p(:,(retained_v+1):iblock_col_size)=1.0_dp
!
!       ENDDO
!       CALL cp_dbcsr_iterator_stop(iter)
!       CALL cp_dbcsr_finalize(matrix_filter)
!
!       ! apply the filter
!       CALL cp_dbcsr_hadamard_product(Fov,matrix_filter,Fov_filtered)
!
!       !!!!!!!!!!!!!!!!!!!!!
!       ! 3. start iterative minimization of the elements to be discarded
!       iteration=0
!       converged=.FALSE.
!       prepare_to_exit=.FALSE.
!       DO
!
!          iteration=iteration+1
!
!          !!!!!!!!!!!!!!!!!!!!!!!!!
!          ! 4. compute the gradient
!          CALL cp_dbcsr_set(grad_blk,0.0_dp)
!          ! create the diagonal blocks only
!          CALL cp_dbcsr_add_on_diag(grad_blk,1.0_dp)
!
!          CALL cp_dbcsr_multiply("T","N",2.0_dp,Fov_filtered,Fov,&
!                  0.0_dp,grad_blk,retain_sparsity=.TRUE.,&
!                  filter_eps=almo_scf_env%eps_filter)
!          CALL cp_dbcsr_multiply("T","N",-2.0_dp,Fov,Fov_filtered,&
!                  1.0_dp,grad_blk,retain_sparsity=.TRUE.,&
!                  filter_eps=almo_scf_env%eps_filter)
!
!          !!!!!!!!!!!!!!!!!!!!!!!
!          ! 5. check convergence
!          obj_function = 0.5_dp*(cp_dbcsr_frobenius_norm(Fov_filtered))**2
!          grad_norm = cp_dbcsr_frobenius_norm(grad_blk)
!          converged=(grad_norm.lt.almo_scf_env%truncate_v_eps_convergence)
!          IF (converged.OR.(iteration.ge.almo_scf_env%truncate_v_max_iter)) THEN
!             prepare_to_exit=.TRUE.
!          ENDIF
!
!          IF (.NOT.prepare_to_exit) THEN
!
!             !!!!!!!!!!!!!!!!!!!!!!!
!             ! 6. perform steps in the direction of the gradient
!             !    a. first, perform a trial step to "see" the parameters
!             !       of the parabola along the gradient:
!             !       a0 * x^2 + b0 * x + c0
!             !    b. then perform the step to the bottom of the parabola
!
!             ! get c0
!             c0 = obj_function
!             ! get b0 <= d_f/d_alpha along grad
!             !!!CALL cp_dbcsr_multiply("N","N",4.0_dp,Fov,grad_blk,&
!             !!!        0.0_dp,temp0_ov,&
!             !!!        filter_eps=almo_scf_env%eps_filter)
!             !!!CALL cp_dbcsr_trace(Fov_filtered,temp0_ov,b0,"T","N")
!
!             alpha=almo_scf_env%truncate_v_trial_step_size
!
!             line_search_step_last=3
!             DO line_search_step=1,line_search_step_last
!                CALL cp_dbcsr_copy(step_blk,grad_blk)
!                CALL cp_dbcsr_scale(step_blk,-1.0_dp*alpha)
!                CALL generator_to_unitary(step_blk,U_blk,&
!                        almo_scf_env%eps_filter)
!                CALL cp_dbcsr_multiply("N","N",1.0_dp,Fov,U_blk,0.0_dp,temp0_ov,&
!                        filter_eps=almo_scf_env%eps_filter)
!                CALL cp_dbcsr_hadamard_product(temp0_ov,matrix_filter,&
!                        Fov_filtered)
!
!                obj_function_new = 0.5_dp*(cp_dbcsr_frobenius_norm(Fov_filtered))**2
!                IF (line_search_step.eq.1) THEN
!                   ff1 = obj_function_new
!                   step1 = alpha
!                ELSE IF (line_search_step.eq.2) THEN
!                   ff2 = obj_function_new
!                   step2 = alpha
!                ENDIF
!
!                IF (unit_nr>0.AND.(line_search_step.ne.line_search_step_last)) THEN
!                   WRITE(unit_nr,'(T6,A,1X,I3,1X,F10.3,E12.3,E12.3,E12.3)') &
!                         "JOINT_SVD_lin",&
!                         iteration,&
!                         alpha,&
!                         obj_function,&
!                         obj_function_new,&
!                         obj_function_new-obj_function
!                ENDIF
!
!                IF (line_search_step.eq.1) THEN
!                   alpha=2.0_dp*alpha
!                ENDIF
!                IF (line_search_step.eq.2) THEN
!                   a0 = ((ff1-c0)/step1 - (ff2-c0)/step2) / (step1 - step2)
!                   b0 = (ff1-c0)/step1 - a0*step1
!                   ! step size in to the bottom of "the parabola"
!                   alpha=-b0/(2.0_dp*a0)
!                   ! update the default step size
!                   almo_scf_env%truncate_v_trial_step_size=alpha
!                ENDIF
!                !!!IF (line_search_step.eq.1) THEN
!                !!!   a0 = (obj_function_new - b0 * alpha - c0) / (alpha*alpha)
!                !!!   ! step size in to the bottom of "the parabola"
!                !!!   alpha=-b0/(2.0_dp*a0)
!                !!!   !IF (alpha.gt.10.0_dp) alpha=10.0_dp
!                !!!ENDIF
!
!             ENDDO
!
!             ! update Fov and U_blk_tot (use grad_blk as tmp storage)
!             CALL cp_dbcsr_copy(Fov,temp0_ov)
!             CALL cp_dbcsr_multiply("N","N",1.0_dp,U_blk_tot,U_blk,&
!                     0.0_dp,grad_blk,&
!                     filter_eps=almo_scf_env%eps_filter)
!             CALL cp_dbcsr_copy(U_blk_tot,grad_blk)
!
!          ENDIF
!
!          t2 = m_walltime()
!
!          IF (unit_nr>0) THEN
!             WRITE(unit_nr,'(T6,A,1X,I3,1X,F10.3,E12.3,E12.3,E12.3,E12.3,F10.3)') &
!                   "JOINT_SVD_itr",&
!                   iteration,&
!                   alpha,&
!                   obj_function,&
!                   obj_function_new,&
!                   obj_function_new-obj_function,&
!                   grad_norm,&
!                   t2-t1
!                   !(flop1+flop2)/(1.0E6_dp*(t2-t1))
!             CALL m_flush(unit_nr)
!          ENDIF
!
!          t1 = m_walltime()
!
!          IF (prepare_to_exit) EXIT
!
!       ENDDO ! stop iterations
!
!       IF (safe_mode) THEN
!          CALL cp_dbcsr_multiply("N","N",1.0_dp,Fov_original,&
!                  U_blk_tot,0.0_dp,temp0_ov,&
!                  filter_eps=almo_scf_env%eps_filter)
!CALL cp_dbcsr_print(temp0_ov)
!          CALL cp_dbcsr_hadamard_product(temp0_ov,matrix_filter,&
!                  Fov_filtered)
!          obj_function_new = 0.5_dp*(cp_dbcsr_frobenius_norm(Fov_filtered))**2
!
!          IF (unit_nr>0) THEN
!             WRITE(unit_nr,'(T6,A,1X,E12.3)') &
!                   "SANITY CHECK:",&
!                   obj_function_new
!             CALL m_flush(unit_nr)
!          ENDIF
!
!          CALL cp_dbcsr_release(Fov_original)
!       ENDIF
!
!       CALL cp_dbcsr_release(temp0_ov)
!       CALL cp_dbcsr_release(U_blk)
!       CALL cp_dbcsr_release(grad_blk)
!       CALL cp_dbcsr_release(step_blk)
!       CALL cp_dbcsr_release(matrix_filter)
!       CALL cp_dbcsr_release(Fov)
!       CALL cp_dbcsr_release(Fov_filtered)
!
!       ! compute rotated virtual orbitals
!       CALL cp_dbcsr_init(v_full_tmp)
!       CALL cp_dbcsr_create(v_full_tmp,&
!               template=almo_scf_env%matrix_v_full_blk(ispin),&
!               matrix_type=dbcsr_type_no_symmetry)
!       CALL cp_dbcsr_multiply("N","N",1.0_dp,&
!               v_full_new,&
!               matrix_sigma_vv_full_sqrt_inv,0.0_dp,v_full_tmp,&
!               filter_eps=almo_scf_env%eps_filter)
!       CALL cp_dbcsr_multiply("N","N",1.0_dp,&
!               v_full_tmp,&
!               U_blk_tot,0.0_dp,v_full_new,&
!               filter_eps=almo_scf_env%eps_filter)
!
!       CALL cp_dbcsr_release(matrix_sigma_vv_full_sqrt_inv)
!       CALL cp_dbcsr_release(v_full_tmp)
!       CALL cp_dbcsr_release(U_blk_tot)
!
!!!!! orthogonalized virtuals are not blocked
!       ! copy new virtuals into the truncated matrix
!       !CALL cp_dbcsr_work_create(almo_scf_env%matrix_v_blk(ispin),&
!       CALL cp_dbcsr_work_create(almo_scf_env%matrix_v(ispin),&
!               work_mutable=.TRUE.)
!       CALL cp_dbcsr_iterator_start(iter,v_full_new)
!       DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
!
!          CALL cp_dbcsr_iterator_next_block(iter,iblock_row,iblock_col,data_p,&
!                  row_size=iblock_row_size,col_size=iblock_col_size)
!
!          retained_v=almo_scf_env%nvirt_of_domain(iblock_col,ispin)
!
!          NULLIFY (p_new_block)
!          !CALL cp_dbcsr_reserve_block2d(almo_scf_env%matrix_v_blk(ispin),&
!          CALL cp_dbcsr_reserve_block2d(almo_scf_env%matrix_v(ispin),&
!                  iblock_row,iblock_col,p_new_block)
!          CPPostcondition(ASSOCIATED(p_new_block),cp_failure_level,routineP,failure)
!          CPPrecondition(retained_v.gt.0,cp_failure_level,routineP,failure)
!          p_new_block(:,:) = data_p(:,1:retained_v)
!
!       ENDDO ! iterator
!       CALL cp_dbcsr_iterator_stop(iter)
!       !!CALL cp_dbcsr_finalize(almo_scf_env%matrix_v_blk(ispin))
!       CALL cp_dbcsr_finalize(almo_scf_env%matrix_v(ispin))
!
!       CALL cp_dbcsr_release(v_full_new)
!
!    ENDDO ! ispin
!
!    CALL timestop(handle)
!
!  END SUBROUTINE truncate_subspace_v_blk

!UNMERGED ! *****************************************************************************
!UNMERGED !> \brief Compute MO coeffs from the main optimized variable (e.g. Theta, X) 
!UNMERGED !> \par History
!UNMERGED !>       2015.03 created [Rustam Z Khaliullin]
!UNMERGED !> \author Rustam Z Khaliullin
!UNMERGED ! *****************************************************************************
!UNMERGED   SUBROUTINE compute_xalmos_from_main_var(m_var_in, m_t_out, m_quench_t, &
!UNMERGED      m_t0, m_siginv, m_STsiginv0, m_s, m_sig_sqrti_ii_out, domain_r_down, &
!UNMERGED      domain_s_inv, domain_map, cpu_of_domain, assume_t0_q0x, just_started, &
!UNMERGED      optimize_theta, normalize_orbitals, envelope_amplitude, eps_filter, &
!UNMERGED      special_case, error)
!UNMERGED 
!UNMERGED     TYPE(cp_dbcsr_type), INTENT(IN)             :: m_var_in
!UNMERGED     TYPE(cp_dbcsr_type), INTENT(INOUT)          :: m_t_out
!UNMERGED     TYPE(cp_dbcsr_type), INTENT(INOUT)          :: m_quench_t
!UNMERGED     TYPE(cp_dbcsr_type), INTENT(IN)             :: m_t0
!UNMERGED     TYPE(cp_dbcsr_type), INTENT(INOUT)          :: m_siginv
!UNMERGED     TYPE(cp_dbcsr_type), INTENT(IN)             :: m_STsiginv0
!UNMERGED     TYPE(cp_dbcsr_type), INTENT(IN)             :: m_s
!UNMERGED     TYPE(cp_dbcsr_type), INTENT(INOUT)          :: m_sig_sqrti_ii_out
!UNMERGED     TYPE(domain_submatrix_type), DIMENSION(:),&
!UNMERGED        INTENT(IN)                               :: domain_s_inv
!UNMERGED     TYPE(domain_submatrix_type), DIMENSION(:),&
!UNMERGED        INTENT(IN)                               :: domain_r_down
!UNMERGED     INTEGER, DIMENSION(:), INTENT(IN)           :: cpu_of_domain
!UNMERGED     TYPE(domain_map_type), INTENT(IN)           :: domain_map
!UNMERGED     LOGICAL, INTENT(IN)                         :: assume_t0_q0x
!UNMERGED     LOGICAL, INTENT(IN)                         :: just_started
!UNMERGED     LOGICAL, INTENT(IN)                         :: optimize_theta
!UNMERGED     LOGICAL, INTENT(IN)                         :: normalize_orbitals
!UNMERGED     REAL(KIND=dp), INTENT(IN)                   :: envelope_amplitude
!UNMERGED     REAL(KIND=dp), INTENT(IN)                   :: eps_filter
!UNMERGED     INTEGER, INTENT(IN)                         :: special_case
!UNMERGED     TYPE(cp_error_type), INTENT(INOUT)          :: error
!UNMERGED     
!UNMERGED     CHARACTER(len=*), PARAMETER :: &
!UNMERGED       routineN = 'compute_xalmos_from_main_var', &
!UNMERGED       routineP = moduleN//':'//routineN
!UNMERGED 
!UNMERGED     INTEGER                                     :: handle
!UNMERGED     LOGICAL                                     :: failure
!UNMERGED     TYPE(cp_dbcsr_type)                         :: m_tmp_no_1,&
!UNMERGED                                                    m_tmp_oo_1
!UNMERGED     CALL timeset(routineN,handle)
!UNMERGED     
!UNMERGED     CALL cp_dbcsr_init(m_tmp_no_1,error=error)
!UNMERGED     CALL cp_dbcsr_create(m_tmp_no_1,&
!UNMERGED             template=m_quench_t,&
!UNMERGED             matrix_type=dbcsr_type_no_symmetry,&
!UNMERGED             error=error)
!UNMERGED     CALL cp_dbcsr_init(m_tmp_oo_1,error=error)
!UNMERGED     CALL cp_dbcsr_create(m_tmp_oo_1,&
!UNMERGED             template=m_siginv,&
!UNMERGED             matrix_type=dbcsr_type_no_symmetry,&
!UNMERGED             error=error)
!UNMERGED                   
!UNMERGED     ! compute the MO coefficients from theta
!UNMERGED     IF (assume_t0_q0x.AND.just_started) THEN
!UNMERGED        CALL cp_dbcsr_set(m_t_out,0.0_dp,error=error)
!UNMERGED     ENDIF
!UNMERGED     !IF (optimize_theta) THEN
!UNMERGED        !CALL cp_dbcsr_norm(m_var_in,&
!UNMERGED        !        dbcsr_norm_maxabsnorm, norm_scalar=t_norm, error=error)
!UNMERGED        !IF (unit_nr>0) THEN
!UNMERGED        !   WRITE(unit_nr,*) "VAR1: |T|: ", t_norm
!UNMERGED        !   !WRITE(unit_nr,*) "VAR1: SIN(|T|/A), COS(|T|/A), A*SIN(|T|/A): ",&
!UNMERGED        !   !   SIN(t_norm/envelope_amplitude),&
!UNMERGED        !   !   COS(t_norm/envelope_amplitude),&
!UNMERGED        !   !   envelope_amplitude*SIN(t_norm/envelope_amplitude)
!UNMERGED        !   WRITE(unit_nr,*) "VAR1: TANH(|T|/A), D(A*TANH(|T|/A)), A*TANH(|T|/A): ",&
!UNMERGED        !      TANH(t_norm/envelope_amplitude),&
!UNMERGED        !      1.0_dp-(TANH(t_norm/envelope_amplitude))**2,&
!UNMERGED        !      envelope_amplitude*TANH(t_norm/envelope_amplitude)
!UNMERGED        !ENDIF
!UNMERGED     !ENDIF
!UNMERGED     CALL cp_dbcsr_copy(m_tmp_no_1,m_var_in,error=error)
!UNMERGED     IF (optimize_theta) THEN
!UNMERGED        CALL cp_dbcsr_function_of_elements(m_tmp_no_1,&
!UNMERGED                !func=dbcsr_func_sin,&
!UNMERGED                func=dbcsr_func_tanh,&
!UNMERGED                a0=0.0_dp,&
!UNMERGED                a1=1.0_dp/envelope_amplitude,&
!UNMERGED                error=error)
!UNMERGED        CALL cp_dbcsr_scale(m_tmp_no_1,&
!UNMERGED                envelope_amplitude,error=error)
!UNMERGED     ENDIF
!UNMERGED     CALL cp_dbcsr_hadamard_product(m_tmp_no_1,m_quench_t,&
!UNMERGED             m_t_out,error=error)
!UNMERGED     
!UNMERGED     !CALL cp_dbcsr_norm(m_t_out,&
!UNMERGED     !        dbcsr_norm_maxabsnorm, norm_scalar=t_norm, error=error)
!UNMERGED     !IF (optimize_theta) THEN
!UNMERGED     !   IF (unit_nr>0) THEN
!UNMERGED     !      WRITE(unit_nr,*) "VAR1: |B*A*TANH(T/A)|: ", t_norm
!UNMERGED     !   ENDIF
!UNMERGED     !ELSE
!UNMERGED     !   IF (unit_nr>0) THEN
!UNMERGED     !      WRITE(unit_nr,*) "VAR1: |B*X|: ", t_norm
!UNMERGED     !   ENDIF
!UNMERGED     !ENDIF
!UNMERGED     
!UNMERGED     ! project out R_0
!UNMERGED     IF (assume_t0_q0x.AND.(.NOT.just_started)) THEN
!UNMERGED        IF (special_case.eq.xalmo_case_fully_deloc) THEN
!UNMERGED           CALL cp_dbcsr_multiply("T","N",1.0_dp,&
!UNMERGED                   m_STsiginv0,&
!UNMERGED                   m_t_out,&
!UNMERGED                   0.0_dp,m_tmp_oo_1,&
!UNMERGED                   filter_eps=eps_filter,&
!UNMERGED                   error=error)
!UNMERGED           CALL cp_dbcsr_multiply("N","N",-1.0_dp,&
!UNMERGED                   m_t0,&
!UNMERGED                   m_tmp_oo_1,&
!UNMERGED                   1.0_dp,m_t_out,&
!UNMERGED                   filter_eps=eps_filter,&
!UNMERGED                   error=error)
!UNMERGED        ELSE IF (special_case.eq.xalmo_case_block_diag) THEN
!UNMERGED           ! cannot use projector with block-daigonal ALMOs
!UNMERGED           CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
!UNMERGED        ELSE
!UNMERGED           ! no special case
!UNMERGED           CALL apply_domain_operators(&
!UNMERGED                   matrix_in=m_t_out,&
!UNMERGED                   matrix_out=m_tmp_no_1,&
!UNMERGED                   operator1=domain_r_down,&
!UNMERGED                   operator2=domain_s_inv,&
!UNMERGED                   dpattern=m_quench_t,&
!UNMERGED                   map=domain_map,&
!UNMERGED                   node_of_domain=cpu_of_domain,&
!UNMERGED                   my_action=1,&
!UNMERGED                   filter_eps=eps_filter,&
!UNMERGED                   use_trimmer=.FALSE.,&
!UNMERGED                   error=error)
!UNMERGED           CALL cp_dbcsr_copy(m_t_out,&
!UNMERGED                   m_tmp_no_1,error=error)
!UNMERGED        ENDIF ! special case
!UNMERGED        !CALL cp_dbcsr_norm(m_t_out,&
!UNMERGED        !        dbcsr_norm_maxabsnorm, norm_scalar=t_norm, error=error)
!UNMERGED        !IF (unit_nr>0) THEN
!UNMERGED        !   WRITE(unit_nr,*) "VAR1: |Q.X|: ", t_norm
!UNMERGED        !ENDIF
!UNMERGED     ENDIF ! end assume_t0_q0x
!UNMERGED 
!UNMERGED     IF (assume_t0_q0x) THEN
!UNMERGED        CALL cp_dbcsr_add(m_t_out,&
!UNMERGED                m_t0,1.0_dp,1.0_dp,&
!UNMERGED                error=error)
!UNMERGED        !CALL cp_dbcsr_norm(m_t_out,&
!UNMERGED        !        dbcsr_norm_maxabsnorm, norm_scalar=t_norm, error=error)
!UNMERGED        !IF (unit_nr>0) THEN
!UNMERGED        !   WRITE(unit_nr,*) "VAR1: |T0+Q.X|: ", t_norm
!UNMERGED        !ENDIF
!UNMERGED     ENDIF
!UNMERGED 
!UNMERGED     IF (normalize_orbitals) THEN
!UNMERGED        CALL orthogonalize_mos(&
!UNMERGED                ket=m_t_out,&
!UNMERGED                overlap=m_tmp_oo_1,&
!UNMERGED                metric=m_s,&
!UNMERGED                retain_locality=.TRUE.,&
!UNMERGED                only_normalize=.TRUE.,&
!UNMERGED                eps_filter=eps_filter,&
!UNMERGED                overlap_sqrti=m_sig_sqrti_ii_out,&
!UNMERGED                error=error)
!UNMERGED     ENDIF
!UNMERGED 
!UNMERGED     CALL cp_dbcsr_filter(m_t_out,eps=eps_filter,error=error)
!UNMERGED     
!UNMERGED     CALL cp_dbcsr_release(m_tmp_no_1,error=error)
!UNMERGED     CALL cp_dbcsr_release(m_tmp_oo_1,error=error)
!UNMERGED 
!UNMERGED     CALL timestop(handle)
!UNMERGED   
!UNMERGED   END SUBROUTINE compute_xalmos_from_main_var

! *****************************************************************************
!> \brief Compute the gradient wrt the main variable (e.g. Theta, X) 
!> \par History
!>       2015.03 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE compute_gradient(m_grad_out, m_ks, m_s, m_t, m_t0, &
     m_siginv, m_quench_t, m_FTsiginv, m_siginvTFTsiginv, m_ST, m_STsiginv0, &
     m_theta, domain_s_inv, domain_r_down, &
     cpu_of_domain, domain_map, assume_t0_q0x, optimize_theta, &
     perturbation_only, normalize_orbitals, penalty_occ_vol, &
     penalty_occ_vol_prefactor, envelope_amplitude, eps_filter, spin_factor, &
     special_case, m_sig_sqrti_ii)

    TYPE(cp_dbcsr_type), INTENT(INOUT)          :: m_grad_out
    !TYPE(cp_dbcsr_type), INTENT(INOUT)          :: m_tmp_no_2
    TYPE(cp_dbcsr_type), INTENT(IN)             :: m_ks
    TYPE(cp_dbcsr_type), INTENT(IN)             :: m_s
    TYPE(cp_dbcsr_type), INTENT(IN)             :: m_t
    TYPE(cp_dbcsr_type), INTENT(IN)             :: m_t0
    TYPE(cp_dbcsr_type), INTENT(IN)             :: m_siginv
    TYPE(cp_dbcsr_type), INTENT(IN)             :: m_quench_t
    TYPE(cp_dbcsr_type), INTENT(INOUT)          :: m_FTsiginv
    TYPE(cp_dbcsr_type), INTENT(INOUT)          :: m_siginvTFTsiginv
    TYPE(cp_dbcsr_type), INTENT(INOUT)          :: m_ST
    TYPE(cp_dbcsr_type), INTENT(IN)             :: m_STsiginv0
    TYPE(cp_dbcsr_type), INTENT(IN)             :: m_theta
    TYPE(domain_submatrix_type), DIMENSION(:),&
       INTENT(IN)                               :: domain_s_inv
    TYPE(domain_submatrix_type), DIMENSION(:),&
       INTENT(IN)                               :: domain_r_down
    INTEGER, DIMENSION(:), INTENT(IN)           :: cpu_of_domain
    TYPE(domain_map_type), INTENT(IN)           :: domain_map
    LOGICAL, INTENT(IN)                         :: assume_t0_q0x
    LOGICAL, INTENT(IN)                         :: optimize_theta
    LOGICAL, INTENT(IN)                         :: perturbation_only
    LOGICAL, INTENT(IN)                         :: normalize_orbitals
    LOGICAL, INTENT(IN)                         :: penalty_occ_vol
    REAL(KIND=dp), INTENT(IN)                   :: penalty_occ_vol_prefactor
    REAL(KIND=dp), INTENT(IN)                   :: envelope_amplitude
    REAL(KIND=dp), INTENT(IN)                   :: eps_filter
    REAL(KIND=dp), INTENT(IN)                   :: spin_factor
    INTEGER, INTENT(IN)                         :: special_case
    TYPE(cp_dbcsr_type), OPTIONAL, INTENT(IN)   :: m_sig_sqrti_ii
    
    CHARACTER(len=*), PARAMETER :: &
      routineN = 'compute_gradient', &
      routineP = moduleN//':'//routineN

    INTEGER                                     :: handle, dim0, ndomains
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE    :: tg_diagonal
    REAL(KIND=dp)                               :: energy_g_norm, &
                                                   penalty_occ_vol_g_norm
    LOGICAL                                     :: failure
    TYPE(cp_dbcsr_type)                         :: m_tmp_no_1,&
                                                   m_tmp_no_2,&
                                                   m_tmp_no_3,&
                                                   m_tmp_oo_1,&
                                                   m_tmp_oo_2
    TYPE(domain_submatrix_type), &
      ALLOCATABLE, DIMENSION(:)                 :: domain_r_down_current
    
    CALL timeset(routineN,handle)

    IF ( normalize_orbitals .AND. (.NOT.PRESENT(m_sig_sqrti_ii)) ) THEN
       CPABORT("Matrix must be present")
    ENDIF
    
    CALL cp_dbcsr_init(m_tmp_no_1)
    CALL cp_dbcsr_create(m_tmp_no_1,&
            template=m_quench_t,&
            matrix_type=dbcsr_type_no_symmetry)
    CALL cp_dbcsr_init(m_tmp_no_2)
    CALL cp_dbcsr_create(m_tmp_no_2,&
            template=m_quench_t,&
            matrix_type=dbcsr_type_no_symmetry)
    CALL cp_dbcsr_init(m_tmp_no_3)
    CALL cp_dbcsr_create(m_tmp_no_3,&
            template=m_quench_t,&
            matrix_type=dbcsr_type_no_symmetry)
    CALL cp_dbcsr_init(m_tmp_oo_1)
    CALL cp_dbcsr_create(m_tmp_oo_1,&
            template=m_siginv,&
            matrix_type=dbcsr_type_no_symmetry)
    CALL cp_dbcsr_init(m_tmp_oo_2)
    CALL cp_dbcsr_create(m_tmp_oo_2,&
            template=m_siginv,&
            matrix_type=dbcsr_type_no_symmetry)

    ! do d_E/d_T first
    IF (.NOT.perturbation_only) THEN
       CALL cp_dbcsr_multiply("N","N",1.0_dp,&
               m_ks,&
               m_t,&
               0.0_dp,m_tmp_no_1,&
               filter_eps=eps_filter)
       CALL cp_dbcsr_multiply("N","N",1.0_dp,&
               m_tmp_no_1,&
               m_siginv,&
               0.0_dp,m_FTsiginv,&
               filter_eps=eps_filter)
    ENDIF
    CALL cp_dbcsr_copy(m_tmp_no_2,m_quench_t)
    CALL cp_dbcsr_copy(m_tmp_no_2,m_FTsiginv,keep_sparsity=.TRUE.)
    CALL cp_dbcsr_multiply("T","N",1.0_dp,&
            m_t,&
            m_FTsiginv,&
            0.0_dp,m_tmp_oo_1,&
            filter_eps=eps_filter)
    CALL cp_dbcsr_multiply("N","N",1.0_dp,&
            m_siginv,&
            m_tmp_oo_1,&
            0.0_dp,m_siginvTFTsiginv,&
            filter_eps=eps_filter)
    CALL cp_dbcsr_multiply("N","N",1.0_dp,&
            m_s,&
            m_t,&
            0.0_dp,m_ST,&
            filter_eps=eps_filter)
    CALL cp_dbcsr_multiply("N","N",-1.0_dp,&
            m_ST,&
            m_siginvTFTsiginv,&
            1.0_dp,m_tmp_no_2,&
            retain_sparsity=.TRUE.)
    CALL cp_dbcsr_scale(m_tmp_no_2,2.0_dp*spin_factor)
    
    ! add penalty on the occupied volume: det(sigma)
    IF (penalty_occ_vol) THEN
       !RZK-warning CALL cp_dbcsr_multiply("N","N",&
       !RZK-warning         penalty_occ_vol_prefactor,&
       !RZK-warning         m_ST,&
       !RZK-warning         m_siginv,&
       !RZK-warning         1.0_dp,m_tmp_no_2,&
       !RZK-warning         retain_sparsity=.TRUE.,&
       !RZK-warning         )
CALL cp_dbcsr_copy(m_tmp_no_1,m_quench_t)
CALL cp_dbcsr_multiply("N","N",&
        penalty_occ_vol_prefactor,&
        m_ST,&
        m_siginv,&
        0.0_dp,m_tmp_no_1,&
        retain_sparsity=.TRUE.)
! this norm does not contain the normalization factors
CALL cp_dbcsr_norm(m_tmp_no_1, dbcsr_norm_maxabsnorm,&
        norm_scalar=penalty_occ_vol_g_norm)
CALL cp_dbcsr_norm(m_tmp_no_2, dbcsr_norm_maxabsnorm,&
        norm_scalar=energy_g_norm)
WRITE(*,"(A30,2F20.10)") "Energy/penalty g norms (no norm): ", energy_g_norm, penalty_occ_vol_g_norm
CALL cp_dbcsr_add(m_tmp_no_2,m_tmp_no_1,1.0_dp,1.0_dp)
    ENDIF
          
    ! take into account the factor from the normalization constraint
    IF (normalize_orbitals) THEN
       
       ! G = ( G - ST.[tr(T).G]_ii ) . [sig_sqrti]_ii
       ! this expression can be simplified to 
       ! G = ( G - c0*ST ) . [sig_sqrti]_ii
       ! where c0 = penalty_occ_vol_prefactor
       ! This is because tr(T).G_Energy = 0 and
       !                 tr(T).G_Penalty = c0*I
       CALL cp_dbcsr_copy(m_tmp_no_1,m_quench_t)
       CALL cp_dbcsr_copy(m_tmp_no_1,m_ST,keep_sparsity=.TRUE.)
       CALL cp_dbcsr_add(m_tmp_no_2,m_tmp_no_1,1.0_dp,-penalty_occ_vol_prefactor)
       CALL cp_dbcsr_copy(m_tmp_no_1,m_quench_t)
       CALL cp_dbcsr_multiply("N","N",1.0_dp,&
               m_tmp_no_2,&
               m_sig_sqrti_ii,&
               0.0_dp,m_tmp_no_1,&
               retain_sparsity=.TRUE.)

       !!! slower way of taking the norm into account
       !!CALL cp_dbcsr_copy(m_tmp_no_1,m_tmp_no_2,)
       !!CALL cp_dbcsr_multiply("N","N",1.0_dp,&
       !!        m_tmp_no_2,&
       !!        m_sig_sqrti_ii,&
       !!        0.0_dp,m_tmp_no_1,&
       !!        retain_sparsity=.TRUE.,&
       !!        )
       !!
       !!! get [tr(T).G]_ii
       !!CALL cp_dbcsr_copy(m_tmp_oo_1,m_sig_sqrti_ii,)
       !!CALL cp_dbcsr_multiply("T","N",1.0_dp,&
       !!        m_t,&
       !!        m_tmp_no_2,&
       !!        0.0_dp,m_tmp_oo_1,&
       !!        retain_sparsity=.TRUE.,&
       !!        )
       !!CALL cp_dbcsr_get_info(m_sig_sqrti_ii, nfullrows_total=dim0 )
       !!ALLOCATE(tg_diagonal(dim0))
       !!CALL cp_dbcsr_get_diag(m_tmp_oo_1,tg_diagonal,)
       !!CALL cp_dbcsr_set(m_tmp_oo_1,0.0_dp,)
       !!CALL cp_dbcsr_set_diag(m_tmp_oo_1,tg_diagonal,)
       !!DEALLOCATE(tg_diagonal)
       !!
       !!CALL cp_dbcsr_multiply("N","N",1.0_dp,&
       !!        m_sig_sqrti_ii,&
       !!        m_tmp_oo_1,&
       !!        0.0_dp,m_tmp_oo_2,&
       !!        filter_eps=eps_filter,&
       !!        )
       !!CALL cp_dbcsr_multiply("N","N",-1.0_dp,&
       !!        m_ST,&
       !!        m_tmp_oo_2,&
       !!        1.0_dp,m_tmp_no_1,&
       !!        retain_sparsity=.TRUE.,&
       !!        )
       
    ELSE
       
       CALL cp_dbcsr_copy(m_tmp_no_1,m_tmp_no_2)
    
    ENDIF ! normalize_orbitals

    ! project out the occupied space from the gradient
    IF (assume_t0_q0x) THEN
       IF (special_case.eq.xalmo_case_fully_deloc) THEN
          CALL cp_dbcsr_copy(m_grad_out,m_tmp_no_1)
          CALL cp_dbcsr_multiply("T","N",1.0_dp,&
                  m_t0,&
                  m_grad_out,&
                  0.0_dp,m_tmp_oo_1,&
                  filter_eps=eps_filter)
          CALL cp_dbcsr_multiply("N","N",-1.0_dp,&
                  m_STsiginv0,&
                  m_tmp_oo_1,&
                  1.0_dp,m_grad_out,&
                  filter_eps=eps_filter)
       ELSE IF (special_case.eq.xalmo_case_block_diag) THEN
          CPABORT("Cannot project the zero-order space from itself")
       ELSE
          ! no special case: normal xALMOs
          CALL apply_domain_operators(&
                  matrix_in=m_tmp_no_1,&
                  matrix_out=m_grad_out,&
                  operator2=domain_r_down(:),&
                  operator1=domain_s_inv(:),&
                  dpattern=m_quench_t,&
                  map=domain_map,&
                  node_of_domain=cpu_of_domain,&
                  my_action=1,&
                  filter_eps=eps_filter,&
                  !matrix_trimmer=,&
                  use_trimmer=.FALSE.)
       ENDIF ! my_special_case
       CALL cp_dbcsr_copy(m_tmp_no_1,m_grad_out)
    ENDIF

    !! check whether the gradient lies entirely in R or Q
    !CALL cp_dbcsr_multiply("T","N",1.0_dp,&
    !        m_t,&
    !        m_tmp_no_1,&
    !        0.0_dp,m_tmp_oo_1,&
    !        filter_eps=eps_filter,&
    !        )
    !CALL cp_dbcsr_multiply("N","N",1.0_dp,&
    !        m_siginv,&
    !        m_tmp_oo_1,&
    !        0.0_dp,m_tmp_oo_2,&
    !        filter_eps=eps_filter,&
    !        )
    !CALL cp_dbcsr_copy(m_tmp_no_2,m_tmp_no_1,)
    !CALL cp_dbcsr_multiply("N","N",-1.0_dp,&
    !        m_ST,&
    !        m_tmp_oo_2,&
    !        1.0_dp,m_tmp_no_2,&
    !        retain_sparsity=.TRUE.,&
    !        )
    !CALL cp_dbcsr_norm(m_tmp_no_2, dbcsr_norm_maxabsnorm,&
    !        norm_scalar=penalty_occ_vol_g_norm, )
    !WRITE(*,"(A50,2F20.10)") "Virtual-space projection of the gradient", penalty_occ_vol_g_norm
    !CALL cp_dbcsr_add(m_tmp_no_2,m_tmp_no_1,1.0_dp,-1.0_dp,)
    !CALL cp_dbcsr_norm(m_tmp_no_2, dbcsr_norm_maxabsnorm,&
    !        norm_scalar=penalty_occ_vol_g_norm, )
    !WRITE(*,"(A50,2F20.10)") "Occupied-space projection of the gradient", penalty_occ_vol_g_norm
    !CALL cp_dbcsr_norm(m_tmp_no_1, dbcsr_norm_maxabsnorm,&
    !        norm_scalar=penalty_occ_vol_g_norm, )
    !WRITE(*,"(A50,2F20.10)") "Full gradient", penalty_occ_vol_g_norm

    ! transform d_E/d_T to d_E/d_theta
    IF (optimize_theta) THEN
       CALL cp_dbcsr_copy(m_tmp_no_2,m_theta)
       CALL cp_dbcsr_function_of_elements(m_tmp_no_2,&
               !func=dbcsr_func_cos,&
               func=dbcsr_func_dtanh,&
               a0=0.0_dp,&
               a1=1.0_dp/envelope_amplitude)
       CALL cp_dbcsr_scale(m_tmp_no_2,&
               envelope_amplitude)
       CALL cp_dbcsr_set(m_tmp_no_3,0.0_dp)
       CALL cp_dbcsr_filter(m_tmp_no_3,eps=eps_filter)
       CALL cp_dbcsr_hadamard_product(m_tmp_no_1,&
               m_tmp_no_2,&
               m_tmp_no_3,&
               b_assume_value=1.0_dp)
       CALL cp_dbcsr_hadamard_product(m_tmp_no_3,&
               m_quench_t,&
               m_grad_out)
    ELSE ! simply copy
       CALL cp_dbcsr_hadamard_product(m_tmp_no_1,&
               m_quench_t,&
               m_grad_out)
    ENDIF
    CALL cp_dbcsr_filter(m_grad_out,eps=eps_filter)

    CALL cp_dbcsr_release(m_tmp_no_1)
    CALL cp_dbcsr_release(m_tmp_no_2)
    CALL cp_dbcsr_release(m_tmp_no_3)
    CALL cp_dbcsr_release(m_tmp_oo_1)
    CALL cp_dbcsr_release(m_tmp_oo_2)

    CALL timestop(handle)
  
  END SUBROUTINE compute_gradient

!FROM UNMERGED ! *****************************************************************************
!FROM UNMERGED !> \brief Compute the preconditioner matrices and invert them if necessary 
!FROM UNMERGED !> \par History
!FROM UNMERGED !>       2015.03 created [Rustam Z Khaliullin]
!FROM UNMERGED !> \author Rustam Z Khaliullin
!FROM UNMERGED ! *****************************************************************************
!FROM UNMERGED   SUBROUTINE compute_preconditioner(domain_prec_out, m_prec_out, m_ks, m_s, &
!FROM UNMERGED      m_siginv, m_quench_t,  m_FTsiginv, m_siginvTFTsiginv, m_ST, &
!FROM UNMERGED      m_STsiginv_out, m_s_vv_out, m_f_vv_out, para_env, &
!FROM UNMERGED      blacs_env, nocc_of_domain, domain_s_inv, domain_r_down, cpu_of_domain, &
!FROM UNMERGED      domain_map, assume_t0_q0x, penalty_occ_vol, penalty_occ_vol_prefactor, &
!FROM UNMERGED      eps_filter, spin_factor, special_case, error)
!FROM UNMERGED 
!FROM UNMERGED     TYPE(domain_submatrix_type), DIMENSION(:),&
!FROM UNMERGED        INTENT(INOUT)                            :: domain_prec_out
!FROM UNMERGED     TYPE(cp_dbcsr_type), INTENT(INOUT)          :: m_prec_out
!FROM UNMERGED     TYPE(cp_dbcsr_type), INTENT(INOUT)          :: m_ks
!FROM UNMERGED     TYPE(cp_dbcsr_type), INTENT(INOUT)          :: m_s
!FROM UNMERGED     TYPE(cp_dbcsr_type), INTENT(IN)             :: m_siginv
!FROM UNMERGED     TYPE(cp_dbcsr_type), INTENT(INOUT)          :: m_quench_t
!FROM UNMERGED     TYPE(cp_dbcsr_type), INTENT(IN)             :: m_FTsiginv
!FROM UNMERGED     TYPE(cp_dbcsr_type), INTENT(IN)             :: m_siginvTFTsiginv
!FROM UNMERGED     TYPE(cp_dbcsr_type), INTENT(IN)             :: m_ST
!FROM UNMERGED     TYPE(cp_dbcsr_type), INTENT(INOUT), OPTIONAL   :: m_STsiginv_out
!FROM UNMERGED     TYPE(cp_dbcsr_type), INTENT(INOUT), OPTIONAL   :: m_s_vv_out
!FROM UNMERGED     TYPE(cp_dbcsr_type), INTENT(INOUT), OPTIONAL   :: m_f_vv_out
!FROM UNMERGED     TYPE(cp_para_env_type),  POINTER            :: para_env
!FROM UNMERGED     TYPE(cp_blacs_env_type), POINTER            :: blacs_env
!FROM UNMERGED 
!FROM UNMERGED     INTEGER, DIMENSION(:), INTENT(IN)           :: nocc_of_domain 
!FROM UNMERGED     TYPE(domain_submatrix_type), DIMENSION(:),&
!FROM UNMERGED        INTENT(IN)                               :: domain_s_inv
!FROM UNMERGED     TYPE(domain_submatrix_type), DIMENSION(:),&
!FROM UNMERGED        INTENT(IN)                               :: domain_r_down
!FROM UNMERGED     INTEGER, DIMENSION(:), INTENT(IN)           :: cpu_of_domain
!FROM UNMERGED     TYPE(domain_map_type), INTENT(IN)           :: domain_map
!FROM UNMERGED     LOGICAL, INTENT(IN)                         :: assume_t0_q0x
!FROM UNMERGED     LOGICAL, INTENT(IN)                         :: penalty_occ_vol
!FROM UNMERGED     REAL(KIND=dp), INTENT(IN)                   :: penalty_occ_vol_prefactor
!FROM UNMERGED     REAL(KIND=dp), INTENT(IN)                   :: eps_filter
!FROM UNMERGED     REAL(KIND=dp), INTENT(IN)                   :: spin_factor
!FROM UNMERGED     INTEGER, INTENT(IN)                         :: special_case
!FROM UNMERGED     TYPE(cp_error_type), INTENT(INOUT)          :: error
!FROM UNMERGED     
!FROM UNMERGED     CHARACTER(len=*), PARAMETER :: &
!FROM UNMERGED       routineN = 'compute_preconditioner', &
!FROM UNMERGED       routineP = moduleN//':'//routineN
!FROM UNMERGED 
!FROM UNMERGED     INTEGER                                     :: handle,&
!FROM UNMERGED                                                    precond_domain_projector
!FROM UNMERGED     LOGICAL                                     :: failure
!FROM UNMERGED     TYPE(cp_dbcsr_type)                         :: m_tmp_nn_1,&
!FROM UNMERGED                                                    m_tmp_no_3
!FROM UNMERGED 
!FROM UNMERGED     CALL timeset(routineN,handle)
!FROM UNMERGED     
!FROM UNMERGED     CALL cp_dbcsr_init(m_tmp_nn_1,error=error)
!FROM UNMERGED     CALL cp_dbcsr_create(m_tmp_nn_1,&
!FROM UNMERGED             template=m_s,&
!FROM UNMERGED             matrix_type=dbcsr_type_no_symmetry,&
!FROM UNMERGED             error=error)
!FROM UNMERGED     CALL cp_dbcsr_init(m_tmp_no_3,error=error)
!FROM UNMERGED     CALL cp_dbcsr_create(m_tmp_no_3,&
!FROM UNMERGED             template=m_quench_t,&
!FROM UNMERGED             matrix_type=dbcsr_type_no_symmetry,&
!FROM UNMERGED             error=error)
!FROM UNMERGED 
!FROM UNMERGED     ! calculate (1-R)F(1-R) and S-SRS
!FROM UNMERGED     ! RZK-warning take advantage: some elements will be removed by the quencher
!FROM UNMERGED     ! RZK-warning S operations can be performed outside the spin loop to save time
!FROM UNMERGED     ! IT IS REQUIRED THAT PRECONDITIONER DOES NOT BREAK THE LOCALITY!!!!
!FROM UNMERGED     ! RZK-warning: further optimization is ABSOLUTELY NECESSARY
!FROM UNMERGED     
!FROM UNMERGED     ! First S-SRS
!FROM UNMERGED     CALL cp_dbcsr_multiply("N","N",1.0_dp,&
!FROM UNMERGED             m_ST,&
!FROM UNMERGED             m_siginv,&
!FROM UNMERGED             0.0_dp,m_tmp_no_3,&
!FROM UNMERGED             filter_eps=eps_filter,&
!FROM UNMERGED             error=error)
!FROM UNMERGED     CALL cp_dbcsr_desymmetrize(m_s,&
!FROM UNMERGED             m_tmp_nn_1,error=error)
!FROM UNMERGED     ! return STsiginv if necessary
!FROM UNMERGED     IF (PRESENT(m_STsiginv_out)) THEN
!FROM UNMERGED        CALL cp_dbcsr_copy(m_STsiginv_out,m_tmp_no_3,error=error)
!FROM UNMERGED     ENDIF
!FROM UNMERGED !!    IF (special_case.eq.xalmo_case_fully_deloc) THEN
!FROM UNMERGED !!       ! use S instead of S-SRS
!FROM UNMERGED !!    ELSE
!FROM UNMERGED !!       CALL cp_dbcsr_multiply("N","T",-1.0_dp,&
!FROM UNMERGED !!!       CALL cp_dbcsr_multiply("N","T",1.0_dp,&
!FROM UNMERGED !!               m_ST,&
!FROM UNMERGED !!               m_tmp_no_3,&
!FROM UNMERGED !!!               0.0_dp,m_tmp_nn_1,&
!FROM UNMERGED !!               1.0_dp,m_tmp_nn_1,&
!FROM UNMERGED !!               filter_eps=eps_filter,&
!FROM UNMERGED !!               error=error)
!FROM UNMERGED !!    ENDIF
!FROM UNMERGED     ! return S_vv = (S or S-SRS) if necessary
!FROM UNMERGED     IF (PRESENT(m_s_vv_out)) THEN
!FROM UNMERGED        CALL cp_dbcsr_copy(m_s_vv_out,m_tmp_nn_1,error=error)
!FROM UNMERGED     ENDIF
!FROM UNMERGED     
!FROM UNMERGED     ! Second (1-R)F(1-R)
!FROM UNMERGED     ! re-create matrix because desymmetrize is buggy -
!FROM UNMERGED     ! it will create multiple copies of blocks
!FROM UNMERGED     CALL cp_dbcsr_desymmetrize(m_ks,m_prec_out,error=error)
!FROM UNMERGED     CALL cp_dbcsr_multiply("N","T",-1.0_dp,&
!FROM UNMERGED             m_FTsiginv,&
!FROM UNMERGED             m_ST,&
!FROM UNMERGED             1.0_dp,m_prec_out,&
!FROM UNMERGED             filter_eps=eps_filter,&
!FROM UNMERGED             error=error)
!FROM UNMERGED     CALL cp_dbcsr_multiply("N","T",-1.0_dp,&
!FROM UNMERGED             m_ST,&
!FROM UNMERGED             m_FTsiginv,&
!FROM UNMERGED             1.0_dp,m_prec_out,&
!FROM UNMERGED             filter_eps=eps_filter,&
!FROM UNMERGED             error=error)
!FROM UNMERGED     CALL cp_dbcsr_multiply("N","N",1.0_dp,&
!FROM UNMERGED             m_ST,&
!FROM UNMERGED             m_siginvTFTsiginv,&
!FROM UNMERGED             0.0_dp,m_tmp_no_3,&
!FROM UNMERGED             filter_eps=eps_filter,&
!FROM UNMERGED             error=error)
!FROM UNMERGED     CALL cp_dbcsr_multiply("N","T",1.0_dp,&
!FROM UNMERGED             m_tmp_no_3,&
!FROM UNMERGED             m_ST,&
!FROM UNMERGED             1.0_dp,m_prec_out,&
!FROM UNMERGED             filter_eps=eps_filter,&
!FROM UNMERGED             error=error)
!FROM UNMERGED     ! return F_vv = (I-SR)F(I-RS) if necessary
!FROM UNMERGED     IF (PRESENT(m_f_vv_out)) THEN
!FROM UNMERGED        CALL cp_dbcsr_copy(m_f_vv_out,m_prec_out,error=error)
!FROM UNMERGED     ENDIF
!FROM UNMERGED 
!FROM UNMERGED #if 0
!FROM UNMERGED !penalty_only=.TRUE.
!FROM UNMERGED     WRITE(*,*) "prefactor0:", penalty_occ_vol_prefactor
!FROM UNMERGED     !IF (penalty_occ_vol) THEN
!FROM UNMERGED        CALL cp_dbcsr_desymmetrize(m_s,&
!FROM UNMERGED                m_prec_out,error=error)
!FROM UNMERGED        !CALL cp_dbcsr_scale(m_prec_out,-penalty_occ_vol_prefactor,error=error)
!FROM UNMERGED     !ENDIF
!FROM UNMERGED #else
!FROM UNMERGED     ! sum up the F_vv and S_vv terms
!FROM UNMERGED     CALL cp_dbcsr_add(m_prec_out,m_tmp_nn_1,&
!FROM UNMERGED             1.0_dp,1.0_dp,error=error)
!FROM UNMERGED     ! Scale to obtain unit step length
!FROM UNMERGED     CALL cp_dbcsr_scale(m_prec_out,2.0_dp*spin_factor,error=error)
!FROM UNMERGED     
!FROM UNMERGED     ! add the contribution from the penalty on the occupied volume
!FROM UNMERGED     IF (penalty_occ_vol) THEN
!FROM UNMERGED        CALL cp_dbcsr_add(m_prec_out,m_tmp_nn_1,&
!FROM UNMERGED                1.0_dp,penalty_occ_vol_prefactor,&
!FROM UNMERGED                error=error)
!FROM UNMERGED     ENDIF
!FROM UNMERGED #endif
!FROM UNMERGED 
!FROM UNMERGED     CALL cp_dbcsr_copy(m_tmp_nn_1,m_prec_out,error=error)
!FROM UNMERGED 
!FROM UNMERGED     ! invert using various algorithms
!FROM UNMERGED     IF (special_case.eq.xalmo_case_block_diag) THEN ! non-overlapping diagonal blocks
!FROM UNMERGED 
!FROM UNMERGED        CALL pseudo_invert_diagonal_blk(&
!FROM UNMERGED                matrix_in=m_tmp_nn_1,&
!FROM UNMERGED                matrix_out=m_prec_out,&
!FROM UNMERGED                nocc=nocc_of_domain(:),&
!FROM UNMERGED                error=error)
!FROM UNMERGED 
!FROM UNMERGED     ELSE IF (special_case.eq.xalmo_case_fully_deloc) THEN ! the entire system is a block
!FROM UNMERGED 
!FROM UNMERGED        ! invert using cholesky (works with S matrix, will not work with S-SRS matrix)
!FROM UNMERGED        CALL cp_dbcsr_cholesky_decompose(m_prec_out,&
!FROM UNMERGED                para_env=para_env,&
!FROM UNMERGED                blacs_env=blacs_env,error=error)
!FROM UNMERGED        CALL cp_dbcsr_cholesky_invert(m_prec_out,&
!FROM UNMERGED                para_env=para_env,&
!FROM UNMERGED                blacs_env=blacs_env,&
!FROM UNMERGED                upper_to_full=.TRUE.,error=error)
!FROM UNMERGED        CALL cp_dbcsr_filter(m_prec_out,&
!FROM UNMERGED                eps=eps_filter,&
!FROM UNMERGED                error=error)
!FROM UNMERGED 
!FROM UNMERGED     ELSE
!FROM UNMERGED        
!FROM UNMERGED        !!! use a true domain preconditioner with overlapping domains
!FROM UNMERGED        IF (assume_t0_q0x) THEN
!FROM UNMERGED           precond_domain_projector=-1
!FROM UNMERGED        ELSE
!FROM UNMERGED           precond_domain_projector=0
!FROM UNMERGED        ENDIF
!FROM UNMERGED 
!FROM UNMERGED !CALL copy_submatrices(domain_s_inv(:),&
!FROM UNMERGED !        domain_prec_out(:),&
!FROM UNMERGED !        copy_data=.TRUE.,&
!FROM UNMERGED !        error=error)
!FROM UNMERGED        CALL construct_domain_preconditioner(&
!FROM UNMERGED           matrix_main=m_tmp_nn_1,&
!FROM UNMERGED           subm_s_inv=domain_s_inv(:),&
!FROM UNMERGED           subm_r_down=domain_r_down(:),&
!FROM UNMERGED           matrix_trimmer=m_quench_t,&
!FROM UNMERGED           dpattern=m_quench_t,&
!FROM UNMERGED           map=domain_map,&
!FROM UNMERGED           node_of_domain=cpu_of_domain,&
!FROM UNMERGED           preconditioner=domain_prec_out(:),&
!FROM UNMERGED           use_trimmer=.FALSE.,&
!FROM UNMERGED           my_action=precond_domain_projector,&
!FROM UNMERGED           eps_zero_eigenvalues=eps_filter,&
!FROM UNMERGED           error=error)
!FROM UNMERGED 
!FROM UNMERGED     ENDIF
!FROM UNMERGED 
!FROM UNMERGED     CALL cp_dbcsr_release(m_tmp_nn_1,error=error)
!FROM UNMERGED     CALL cp_dbcsr_release(m_tmp_no_3,error=error)
!FROM UNMERGED 
!FROM UNMERGED     CALL timestop(handle)
!FROM UNMERGED   
!FROM UNMERGED   END SUBROUTINE compute_preconditioner
!FROM UNMERGED 
!FROM UNMERGED ! *****************************************************************************
!FROM UNMERGED !> \brief Compute beta for conjugate gradient algorithms 
!FROM UNMERGED !> \par History
!FROM UNMERGED !>       2015.04 created [Rustam Z Khaliullin]
!FROM UNMERGED !> \author Rustam Z Khaliullin
!FROM UNMERGED ! *****************************************************************************
!FROM UNMERGED   SUBROUTINE compute_cg_beta(beta, numer, denom, reset_conjugator, conjugator,&
!FROM UNMERGED      grad, prev_grad, step, prev_step, prev_minus_prec_grad, error) 
!FROM UNMERGED 
!FROM UNMERGED     REAL(KIND=dp), INTENT(INOUT)                :: beta
!FROM UNMERGED     REAL(KIND=dp), INTENT(INOUT), OPTIONAL      :: numer, denom
!FROM UNMERGED     LOGICAL, INTENT(INOUT)                      :: reset_conjugator
!FROM UNMERGED     INTEGER, INTENT(IN)                         :: conjugator
!FROM UNMERGED     TYPE(cp_dbcsr_type), INTENT(INOUT)          :: grad
!FROM UNMERGED     TYPE(cp_dbcsr_type), INTENT(INOUT)          :: prev_grad
!FROM UNMERGED     TYPE(cp_dbcsr_type), INTENT(INOUT)          :: step
!FROM UNMERGED     TYPE(cp_dbcsr_type), INTENT(INOUT)          :: prev_step
!FROM UNMERGED     TYPE(cp_dbcsr_type), INTENT(INOUT), OPTIONAL &
!FROM UNMERGED                                                 :: prev_minus_prec_grad
!FROM UNMERGED     TYPE(cp_error_type), INTENT(INOUT)          :: error
!FROM UNMERGED     
!FROM UNMERGED     CHARACTER(len=*), PARAMETER :: &
!FROM UNMERGED       routineN = 'compute_cg_beta', &
!FROM UNMERGED       routineP = moduleN//':'//routineN
!FROM UNMERGED 
!FROM UNMERGED     INTEGER                                     :: handle
!FROM UNMERGED     LOGICAL                                     :: failure
!FROM UNMERGED     TYPE(cp_dbcsr_type)                         :: m_tmp_no_1
!FROM UNMERGED     REAL(KIND=dp)                               :: my_numer, my_denom, &
!FROM UNMERGED                                                    kappa, tau
!FROM UNMERGED 
!FROM UNMERGED     CALL timeset(routineN,handle)
!FROM UNMERGED 
!FROM UNMERGED     IF (.NOT.PRESENT(prev_minus_prec_grad)) THEN
!FROM UNMERGED        IF (conjugator.EQ.cg_fletcher_reeves .OR. &
!FROM UNMERGED            conjugator.EQ.cg_polak_ribiere .OR. &
!FROM UNMERGED            conjugator.EQ.cg_hager_zhang ) THEN
!FROM UNMERGED           CPErrorMessage(cp_failure_level,routineP,"conjugator needs more input",error)
!FROM UNMERGED           CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
!FROM UNMERGED        ENDIF
!FROM UNMERGED     ENDIF
!FROM UNMERGED 
!FROM UNMERGED     IF (PRESENT(numer).OR.PRESENT(denom)) THEN
!FROM UNMERGED        IF (conjugator.EQ.cg_hestenes_stiefel .OR. &
!FROM UNMERGED            conjugator.EQ.cg_dai_yuan .OR. &
!FROM UNMERGED            conjugator.EQ.cg_hager_zhang ) THEN
!FROM UNMERGED           CPErrorMessage(cp_failure_level,routineP,"cannot return numer/denom",error)
!FROM UNMERGED           CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
!FROM UNMERGED        ENDIF
!FROM UNMERGED     ENDIF
!FROM UNMERGED 
!FROM UNMERGED     CALL cp_dbcsr_init(m_tmp_no_1,error=error)
!FROM UNMERGED     CALL cp_dbcsr_create(m_tmp_no_1,&
!FROM UNMERGED             template=grad,&
!FROM UNMERGED             matrix_type=dbcsr_type_no_symmetry,&
!FROM UNMERGED             error=error)
!FROM UNMERGED 
!FROM UNMERGED     SELECT CASE (conjugator)
!FROM UNMERGED     CASE (cg_hestenes_stiefel)
!FROM UNMERGED        CALL cp_dbcsr_copy(m_tmp_no_1,grad,error=error)
!FROM UNMERGED        CALL cp_dbcsr_add(m_tmp_no_1,prev_grad,&
!FROM UNMERGED                1.0_dp,-1.0_dp,error=error)
!FROM UNMERGED        CALL cp_dbcsr_trace(m_tmp_no_1,step,my_numer,&
!FROM UNMERGED                "T","N",error=error)
!FROM UNMERGED        CALL cp_dbcsr_trace(m_tmp_no_1,prev_step,my_denom,&
!FROM UNMERGED                "T","N",error=error)
!FROM UNMERGED        beta=-1.0_dp*my_numer/my_denom
!FROM UNMERGED     CASE (cg_fletcher_reeves)
!FROM UNMERGED        CALL cp_dbcsr_trace(grad,step,my_numer,"T","N",error=error)
!FROM UNMERGED        CALL cp_dbcsr_trace(prev_grad,prev_minus_prec_grad,my_denom,"T","N",error=error)
!FROM UNMERGED        beta=my_numer/my_denom
!FROM UNMERGED     CASE (cg_polak_ribiere)
!FROM UNMERGED        CALL cp_dbcsr_trace(prev_grad,prev_minus_prec_grad,my_denom,"T","N",error=error)
!FROM UNMERGED        CALL cp_dbcsr_copy(m_tmp_no_1,grad,error=error)
!FROM UNMERGED        CALL cp_dbcsr_add(m_tmp_no_1,prev_grad,1.0_dp,-1.0_dp,error=error)
!FROM UNMERGED        CALL cp_dbcsr_trace(m_tmp_no_1,step,my_numer,"T","N",error=error)
!FROM UNMERGED        beta=my_numer/my_denom
!FROM UNMERGED     CASE (cg_fletcher)
!FROM UNMERGED        CALL cp_dbcsr_trace(grad,step,my_numer,"T","N",error=error)
!FROM UNMERGED        CALL cp_dbcsr_trace(prev_grad,prev_step,my_denom,"T","N",error=error)
!FROM UNMERGED        beta=my_numer/my_denom
!FROM UNMERGED     CASE (cg_liu_storey)
!FROM UNMERGED        CALL cp_dbcsr_trace(prev_grad,prev_step,my_denom,"T","N",error=error)
!FROM UNMERGED        CALL cp_dbcsr_copy(m_tmp_no_1,grad,error=error)
!FROM UNMERGED        CALL cp_dbcsr_add(m_tmp_no_1,prev_grad,1.0_dp,-1.0_dp,error=error)
!FROM UNMERGED        CALL cp_dbcsr_trace(m_tmp_no_1,step,my_numer,"T","N",error=error)
!FROM UNMERGED        beta=my_numer/my_denom
!FROM UNMERGED     CASE (cg_dai_yuan)
!FROM UNMERGED        CALL cp_dbcsr_trace(grad,step,my_numer,"T","N",error=error)
!FROM UNMERGED        CALL cp_dbcsr_copy(m_tmp_no_1,grad,error=error)
!FROM UNMERGED        CALL cp_dbcsr_add(m_tmp_no_1,prev_grad,1.0_dp,-1.0_dp,error=error)
!FROM UNMERGED        CALL cp_dbcsr_trace(m_tmp_no_1,prev_step,my_denom,"T","N",error=error)
!FROM UNMERGED        beta=-1.0_dp*my_numer/my_denom
!FROM UNMERGED     CASE (cg_hager_zhang)
!FROM UNMERGED        CALL cp_dbcsr_copy(m_tmp_no_1,grad,error=error)
!FROM UNMERGED        CALL cp_dbcsr_add(m_tmp_no_1,prev_grad,1.0_dp,-1.0_dp,error=error)
!FROM UNMERGED        CALL cp_dbcsr_trace(m_tmp_no_1,prev_step,my_denom,"T","N",error=error)
!FROM UNMERGED        CALL cp_dbcsr_trace(m_tmp_no_1,prev_minus_prec_grad,my_numer,"T","N",error=error)
!FROM UNMERGED        kappa=-2.0_dp*my_numer/my_denom
!FROM UNMERGED        CALL cp_dbcsr_trace(m_tmp_no_1,step,my_numer,"T","N",error=error)
!FROM UNMERGED        tau=-1.0_dp*my_numer/my_denom
!FROM UNMERGED        CALL cp_dbcsr_trace(prev_step,grad,my_numer,"T","N",error=error)
!FROM UNMERGED        beta=tau-kappa*my_numer/my_denom
!FROM UNMERGED     CASE (cg_zero)
!FROM UNMERGED        beta=0.0_dp
!FROM UNMERGED     CASE DEFAULT
!FROM UNMERGED        CPErrorMessage(cp_failure_level,routineP,"illegal conjugator",error)
!FROM UNMERGED        CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
!FROM UNMERGED     END SELECT
!FROM UNMERGED     
!FROM UNMERGED     IF (beta.lt.0.0_dp) THEN
!FROM UNMERGED        !IF (unit_nr>0) THEN
!FROM UNMERGED        !   WRITE(unit_nr,*) "Beta is negative: ", beta
!FROM UNMERGED        !ENDIF
!FROM UNMERGED        reset_conjugator=.TRUE.
!FROM UNMERGED     ENDIF
!FROM UNMERGED     
!FROM UNMERGED     CALL cp_dbcsr_release(m_tmp_no_1,error=error)
!FROM UNMERGED     
!FROM UNMERGED     IF (PRESENT(numer)) THEN
!FROM UNMERGED        numer=my_numer
!FROM UNMERGED     ENDIF
!FROM UNMERGED     IF (PRESENT(denom)) THEN
!FROM UNMERGED        denom=my_denom
!FROM UNMERGED     ENDIF
!FROM UNMERGED 
!FROM UNMERGED 
!FROM UNMERGED     CALL timestop(handle)
!FROM UNMERGED   
!FROM UNMERGED   END SUBROUTINE compute_cg_beta
!FROM UNMERGED 
!FROM UNMERGED 
!FROM UNMERGED ! *****************************************************************************
!FROM UNMERGED !> \brief computes the step matrix from the gradient and Hessian using
!FROM UNMERGED !>         the Newton-Raphson method
!FROM UNMERGED !> \par History
!FROM UNMERGED !>       2015.04 created [Rustam Z. Khaliullin]
!FROM UNMERGED !> \author Rustam Z. Khaliullin
!FROM UNMERGED ! *****************************************************************************
!FROM UNMERGED   SUBROUTINE newton_grad_to_step(optimizer, m_grad, m_delta, m_s, m_ks, &
!FROM UNMERGED     m_siginv, m_quench_t, m_FTsiginv, m_siginvTFTsiginv, m_ST, m_t, &
!FROM UNMERGED     m_sig_sqrti_ii, domain_s_inv, domain_r_down, domain_map, cpu_of_domain, &
!FROM UNMERGED     nocc_of_domain, para_env, blacs_env, eps_filter, optimize_theta, &
!FROM UNMERGED     penalty_occ_vol, normalize_orbitals, penalty_occ_vol_prefactor, &
!FROM UNMERGED     penalty_occ_vol_pf2, special_case, error)
!FROM UNMERGED 
!FROM UNMERGED     TYPE(optimizer_options_type)             :: optimizer
!FROM UNMERGED     TYPE(cp_dbcsr_type), &
!FROM UNMERGED        DIMENSION(:), INTENT(IN)              :: m_grad
!FROM UNMERGED     TYPE(cp_dbcsr_type), &
!FROM UNMERGED        DIMENSION(:), INTENT(INOUT)           :: m_delta, &
!FROM UNMERGED                                                 m_s, m_ks
!FROM UNMERGED     TYPE(cp_dbcsr_type), &
!FROM UNMERGED        DIMENSION(:), INTENT(INOUT)           :: m_siginv, &
!FROM UNMERGED                                                 m_quench_t
!FROM UNMERGED     TYPE(cp_dbcsr_type), &
!FROM UNMERGED        DIMENSION(:), INTENT(IN)              :: m_FTsiginv,&
!FROM UNMERGED                                                 m_siginvTFTsiginv,&
!FROM UNMERGED                                                 m_ST, m_t, &
!FROM UNMERGED                                                 m_sig_sqrti_ii
!FROM UNMERGED     TYPE(domain_submatrix_type), &
!FROM UNMERGED        DIMENSION(:,:), INTENT(IN)            :: domain_s_inv, &
!FROM UNMERGED                                                 domain_r_down
!FROM UNMERGED     TYPE(domain_map_type), &
!FROM UNMERGED        DIMENSION(:), INTENT(IN)              :: domain_map
!FROM UNMERGED     INTEGER, DIMENSION(:), INTENT(IN)        :: cpu_of_domain
!FROM UNMERGED     INTEGER, DIMENSION(:,:), INTENT(IN)      :: nocc_of_domain 
!FROM UNMERGED     TYPE(cp_para_env_type),  POINTER         :: para_env
!FROM UNMERGED     TYPE(cp_blacs_env_type), POINTER         :: blacs_env
!FROM UNMERGED 
!FROM UNMERGED     REAL(KIND=dp), INTENT(IN)                :: eps_filter
!FROM UNMERGED     LOGICAL, INTENT(IN)                      :: optimize_theta
!FROM UNMERGED     LOGICAL, INTENT(IN)                      :: penalty_occ_vol
!FROM UNMERGED     LOGICAL, INTENT(IN)                      :: normalize_orbitals
!FROM UNMERGED     REAL(KIND=dp), INTENT(IN)                :: penalty_occ_vol_prefactor
!FROM UNMERGED     REAL(KIND=dp), INTENT(IN)                :: penalty_occ_vol_pf2
!FROM UNMERGED     INTEGER, INTENT(IN)                      :: special_case
!FROM UNMERGED     TYPE(cp_error_type), INTENT(INOUT)       :: error
!FROM UNMERGED 
!FROM UNMERGED     CHARACTER(len=*), PARAMETER :: routineN = 'newton_grad_to_step', &
!FROM UNMERGED       routineP = moduleN//':'//routineN
!FROM UNMERGED 
!FROM UNMERGED     INTEGER                                  :: handle, nspins, unit_nr,&
!FROM UNMERGED                                                 ndomains, ispin, &
!FROM UNMERGED                                                 outer_max_iter, max_iter, &
!FROM UNMERGED                                                 outer_iteration, iteration
!FROM UNMERGED     TYPE(cp_logger_type), POINTER            :: logger
!FROM UNMERGED     CHARACTER(LEN=20)                        :: iter_type
!FROM UNMERGED     REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: residue_max_norm
!FROM UNMERGED     REAL(KIND=dp)                            :: spin_factor, residue_norm, &
!FROM UNMERGED                                                 t1, t2, numer, denom, alpha, &
!FROM UNMERGED                                                 beta, numer_ispin, &
!FROM UNMERGED                                                 denom_ispin, beta_ispin, &
!FROM UNMERGED                                                 eps_error_target
!FROM UNMERGED     TYPE(cp_dbcsr_type), DIMENSION(:), &
!FROM UNMERGED        ALLOCATABLE                           :: m_residue, m_step, &
!FROM UNMERGED                                                 m_zet, m_prec,&
!FROM UNMERGED                                                 m_s_vv, m_f_vv, m_f_vo, &
!FROM UNMERGED                                                 m_Hstep, m_residue_prev, &
!FROM UNMERGED                                                 m_zet_prev, m_STsiginv
!FROM UNMERGED     LOGICAL                                  :: outer_prepare_to_exit, &
!FROM UNMERGED                                                 prepare_to_exit, converged, &
!FROM UNMERGED                                                 reset_conjugator, failure, &
!FROM UNMERGED                                                 use_preconditioner, &
!FROM UNMERGED                                                 do_exact_inversion
!FROM UNMERGED     TYPE(domain_submatrix_type), &
!FROM UNMERGED        DIMENSION(:,:), ALLOCATABLE           :: domain_prec
!FROM UNMERGED     TYPE(cp_dbcsr_type)                      :: m_tmp_oo_1, m_tmp_oo_2
!FROM UNMERGED 
!FROM UNMERGED 
!FROM UNMERGED     CALL timeset(routineN,handle)
!FROM UNMERGED 
!FROM UNMERGED     ! get a useful output_unit
!FROM UNMERGED     logger => cp_error_get_logger(error)
!FROM UNMERGED     IF (logger%para_env%mepos==logger%para_env%source) THEN
!FROM UNMERGED        unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
!FROM UNMERGED     ELSE
!FROM UNMERGED        unit_nr=-1
!FROM UNMERGED     ENDIF
!FROM UNMERGED 
!FROM UNMERGED     !!! Currently for non-theta only
!FROM UNMERGED     IF (optimize_theta) THEN
!FROM UNMERGED        CPErrorMessage(cp_failure_level,routineP,"theta is NYI",error)
!FROM UNMERGED        CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
!FROM UNMERGED     ENDIF
!FROM UNMERGED 
!FROM UNMERGED     ! set optimizer options
!FROM UNMERGED     use_preconditioner = (optimizer%preconditioner.ne.prec_zero)
!FROM UNMERGED     outer_max_iter=optimizer%max_iter_outer_loop
!FROM UNMERGED     max_iter=optimizer%max_iter
!FROM UNMERGED     eps_error_target=optimizer%eps_error
!FROM UNMERGED 
!FROM UNMERGED     ! set key dimensions
!FROM UNMERGED     nspins=SIZE(m_ks)
!FROM UNMERGED     ndomains=SIZE(domain_s_inv,1)
!FROM UNMERGED 
!FROM UNMERGED     IF (nspins == 1) THEN
!FROM UNMERGED        spin_factor = 2.0_dp
!FROM UNMERGED     ELSE
!FROM UNMERGED        spin_factor = 1.0_dp
!FROM UNMERGED     ENDIF
!FROM UNMERGED 
!FROM UNMERGED     ALLOCATE(domain_prec(ndomains,nspins))
!FROM UNMERGED     CALL init_submatrices(domain_prec,error)
!FROM UNMERGED 
!FROM UNMERGED     ! allocate matrices
!FROM UNMERGED     ALLOCATE(m_residue(nspins))
!FROM UNMERGED     ALLOCATE(m_residue_prev(nspins))
!FROM UNMERGED     ALLOCATE(m_step(nspins))
!FROM UNMERGED     ALLOCATE(m_zet(nspins))
!FROM UNMERGED     ALLOCATE(m_zet_prev(nspins))
!FROM UNMERGED     ALLOCATE(m_Hstep(nspins))
!FROM UNMERGED     ALLOCATE(m_prec(nspins))
!FROM UNMERGED     ALLOCATE(m_s_vv(nspins))
!FROM UNMERGED     ALLOCATE(m_f_vv(nspins))
!FROM UNMERGED     ALLOCATE(m_f_vo(nspins))
!FROM UNMERGED     ALLOCATE(m_STsiginv(nspins))
!FROM UNMERGED 
!FROM UNMERGED     ALLOCATE(residue_max_norm(nspins))
!FROM UNMERGED 
!FROM UNMERGED     ! initiate objects before iterations 
!FROM UNMERGED     DO ispin=1,nspins
!FROM UNMERGED 
!FROM UNMERGED        ! init matrices 
!FROM UNMERGED        CALL cp_dbcsr_init(m_residue(ispin),error=error)
!FROM UNMERGED        CALL cp_dbcsr_create(m_residue(ispin),&
!FROM UNMERGED                template=m_quench_t(ispin),&
!FROM UNMERGED                matrix_type=dbcsr_type_no_symmetry,error=error)
!FROM UNMERGED        CALL cp_dbcsr_init(m_residue_prev(ispin),error=error)
!FROM UNMERGED        CALL cp_dbcsr_create(m_residue_prev(ispin),&
!FROM UNMERGED                template=m_quench_t(ispin),&
!FROM UNMERGED                matrix_type=dbcsr_type_no_symmetry,error=error)
!FROM UNMERGED        CALL cp_dbcsr_init(m_step(ispin),error=error)
!FROM UNMERGED        CALL cp_dbcsr_create(m_step(ispin),&
!FROM UNMERGED                template=m_quench_t(ispin),&
!FROM UNMERGED                matrix_type=dbcsr_type_no_symmetry,error=error)
!FROM UNMERGED        CALL cp_dbcsr_init(m_zet_prev(ispin),error=error)
!FROM UNMERGED        CALL cp_dbcsr_create(m_zet_prev(ispin),&
!FROM UNMERGED                template=m_quench_t(ispin),&
!FROM UNMERGED                matrix_type=dbcsr_type_no_symmetry,error=error)
!FROM UNMERGED        CALL cp_dbcsr_init(m_zet(ispin),error=error)
!FROM UNMERGED        CALL cp_dbcsr_create(m_zet(ispin),&
!FROM UNMERGED                template=m_quench_t(ispin),&
!FROM UNMERGED                matrix_type=dbcsr_type_no_symmetry,error=error)
!FROM UNMERGED        CALL cp_dbcsr_init(m_Hstep(ispin),error=error)
!FROM UNMERGED        CALL cp_dbcsr_create(m_Hstep(ispin),&
!FROM UNMERGED                template=m_quench_t(ispin),&
!FROM UNMERGED                matrix_type=dbcsr_type_no_symmetry,error=error)
!FROM UNMERGED        CALL cp_dbcsr_init(m_f_vo(ispin),error=error)
!FROM UNMERGED        CALL cp_dbcsr_create(m_f_vo(ispin),&
!FROM UNMERGED                template=m_quench_t(ispin),&
!FROM UNMERGED                matrix_type=dbcsr_type_no_symmetry,error=error)
!FROM UNMERGED        CALL cp_dbcsr_init(m_STsiginv(ispin),error=error)
!FROM UNMERGED        CALL cp_dbcsr_create(m_STsiginv(ispin),&
!FROM UNMERGED                template=m_quench_t(ispin),&
!FROM UNMERGED                matrix_type=dbcsr_type_no_symmetry,error=error)
!FROM UNMERGED        CALL cp_dbcsr_init(m_f_vv(ispin),error=error)
!FROM UNMERGED        CALL cp_dbcsr_create(m_f_vv(ispin),&
!FROM UNMERGED                template=m_ks(ispin),&
!FROM UNMERGED                matrix_type=dbcsr_type_no_symmetry,error=error)
!FROM UNMERGED        CALL cp_dbcsr_init(m_s_vv(ispin),error=error)
!FROM UNMERGED        CALL cp_dbcsr_create(m_s_vv(ispin),&
!FROM UNMERGED                template=m_s(1),&
!FROM UNMERGED                matrix_type=dbcsr_type_no_symmetry,error=error)
!FROM UNMERGED        CALL cp_dbcsr_init(m_prec(ispin),error=error)
!FROM UNMERGED        CALL cp_dbcsr_create(m_prec(ispin),&
!FROM UNMERGED                template=m_ks(ispin),&
!FROM UNMERGED                matrix_type=dbcsr_type_no_symmetry,error=error)
!FROM UNMERGED 
!FROM UNMERGED        ! compute the full "gradient" - it is necessary to
!FROM UNMERGED        ! evaluate Hessian.X
!FROM UNMERGED        CALL cp_dbcsr_copy(m_f_vo(ispin),m_FTsiginv(ispin),error=error)
!FROM UNMERGED        CALL cp_dbcsr_multiply("N","N",-1.0_dp,&
!FROM UNMERGED                m_ST(ispin),&
!FROM UNMERGED                m_siginvTFTsiginv(ispin),&
!FROM UNMERGED                1.0_dp,m_f_vo(ispin),&
!FROM UNMERGED                filter_eps=eps_filter,&
!FROM UNMERGED                error=error)
!FROM UNMERGED 
!FROM UNMERGED ! RZK-warning
!FROM UNMERGED ! compute preconditioner even if we do not use it 
!FROM UNMERGED ! this is for debugging because compute_preconditioner includes 
!FROM UNMERGED ! computing F_vv and S_vv necessary for 
!FROM UNMERGED !       IF ( use_preconditioner ) THEN
!FROM UNMERGED 
!FROM UNMERGED ! domain_s_inv and domain_r_down are never used with assume_t0_q0x=FALSE          
!FROM UNMERGED           CALL compute_preconditioner(&
!FROM UNMERGED                   domain_prec_out=domain_prec(:,ispin),&
!FROM UNMERGED                   m_prec_out=m_prec(ispin),&
!FROM UNMERGED                   m_ks=m_ks(ispin),&
!FROM UNMERGED                   m_s=m_s(1),&
!FROM UNMERGED                   m_siginv=m_siginv(ispin),&
!FROM UNMERGED                   m_quench_t=m_quench_t(ispin),&
!FROM UNMERGED                   m_FTsiginv=m_FTsiginv(ispin),&
!FROM UNMERGED                   m_siginvTFTsiginv=m_siginvTFTsiginv(ispin),&
!FROM UNMERGED                   m_ST=m_ST(ispin),&
!FROM UNMERGED                   m_STsiginv_out=m_STsiginv(ispin),&
!FROM UNMERGED                   m_s_vv_out=m_s_vv(ispin),&
!FROM UNMERGED                   m_f_vv_out=m_f_vv(ispin),&
!FROM UNMERGED                   para_env=para_env,&
!FROM UNMERGED                   blacs_env=blacs_env,&
!FROM UNMERGED                   nocc_of_domain=nocc_of_domain(:,ispin),&
!FROM UNMERGED                   domain_s_inv=domain_s_inv(:,ispin),&
!FROM UNMERGED                   domain_r_down=domain_r_down(:,ispin),&
!FROM UNMERGED                   cpu_of_domain=cpu_of_domain(:),&
!FROM UNMERGED                   domain_map=domain_map(ispin),&
!FROM UNMERGED                   assume_t0_q0x=.FALSE.,&
!FROM UNMERGED                   penalty_occ_vol=penalty_occ_vol,&
!FROM UNMERGED                   penalty_occ_vol_prefactor=penalty_occ_vol_prefactor,&
!FROM UNMERGED                   eps_filter=eps_filter,&
!FROM UNMERGED                   spin_factor=spin_factor,&
!FROM UNMERGED                   special_case=special_case,&
!FROM UNMERGED                   error=error)
!FROM UNMERGED        
!FROM UNMERGED !       ENDIF ! use_preconditioner
!FROM UNMERGED        
!FROM UNMERGED        ! initial guess
!FROM UNMERGED        CALL cp_dbcsr_copy(m_delta(ispin),m_quench_t(ispin),error=error)
!FROM UNMERGED        ! in order to use cp_dbcsr_set matrix blocks must exist
!FROM UNMERGED        CALL cp_dbcsr_set(m_delta(ispin),0.0_dp,error=error)
!FROM UNMERGED        CALL cp_dbcsr_copy(m_residue(ispin),m_grad(ispin),error=error)
!FROM UNMERGED        CALL cp_dbcsr_scale(m_residue(ispin),-1.0_dp,error=error)
!FROM UNMERGED 
!FROM UNMERGED do_exact_inversion=.FALSE.
!FROM UNMERGED IF (do_exact_inversion) THEN
!FROM UNMERGED 
!FROM UNMERGED    ! copy grad to m_step temporarily
!FROM UNMERGED    ! use m_step as input to the inversion routine
!FROM UNMERGED    CALL cp_dbcsr_copy(m_step(ispin),m_grad(ispin),error=error)
!FROM UNMERGED 
!FROM UNMERGED    ! expensive "exact" inversion of the "nearly-exact" Hessian
!FROM UNMERGED    ! hopefully returns Z=-H^(-1).G
!FROM UNMERGED    CALL hessian_diag_apply(&
!FROM UNMERGED            matrix_grad=m_step(ispin),&
!FROM UNMERGED            matrix_step=m_zet(ispin),&
!FROM UNMERGED            matrix_S_ao=m_s_vv(ispin),&
!FROM UNMERGED            matrix_F_ao=m_f_vv(ispin),&
!FROM UNMERGED            !matrix_S_ao=m_s(ispin),&
!FROM UNMERGED            !matrix_F_ao=m_ks(ispin),&
!FROM UNMERGED            matrix_S_mo=m_siginv(ispin),&
!FROM UNMERGED            matrix_F_mo=m_siginvTFTsiginv(ispin),&
!FROM UNMERGED            matrix_S_vo=m_STsiginv(ispin),&
!FROM UNMERGED            matrix_F_vo=m_f_vo(ispin),&
!FROM UNMERGED            quench_t=m_quench_t(ispin),&
!FROM UNMERGED            spin_factor=spin_factor,&
!FROM UNMERGED            eps_zero=eps_filter*10.0_dp,&
!FROM UNMERGED            penalty_occ_vol=penalty_occ_vol,&
!FROM UNMERGED            penalty_occ_vol_prefactor=penalty_occ_vol_prefactor,&
!FROM UNMERGED            penalty_occ_vol_pf2=penalty_occ_vol_pf2,&
!FROM UNMERGED            m_s=m_s(1),&
!FROM UNMERGED            para_env=para_env,&
!FROM UNMERGED            blacs_env=blacs_env,&
!FROM UNMERGED            error=error)
!FROM UNMERGED     ! correct solution by the spin factor
!FROM UNMERGED     !CALL cp_dbcsr_scale(m_zet(ispin),1.0_dp/(2.0_dp*spin_factor),error=error)
!FROM UNMERGED 
!FROM UNMERGED        
!FROM UNMERGED ELSE ! use PCG to solve H.D=-G
!FROM UNMERGED 
!FROM UNMERGED        IF ( use_preconditioner ) THEN
!FROM UNMERGED 
!FROM UNMERGED           IF (special_case.eq.xalmo_case_block_diag .OR. &
!FROM UNMERGED               special_case.eq.xalmo_case_fully_deloc) THEN
!FROM UNMERGED 
!FROM UNMERGED              CALL cp_dbcsr_multiply("N","N",1.0_dp,&
!FROM UNMERGED                      m_prec(ispin),&
!FROM UNMERGED                      m_residue(ispin),&
!FROM UNMERGED                      0.0_dp,m_zet(ispin),&
!FROM UNMERGED                      filter_eps=eps_filter,&
!FROM UNMERGED                      error=error)
!FROM UNMERGED 
!FROM UNMERGED           ELSE
!FROM UNMERGED 
!FROM UNMERGED              CALL apply_domain_operators(&
!FROM UNMERGED                      matrix_in=m_residue(ispin),&
!FROM UNMERGED                      matrix_out=m_zet(ispin),&
!FROM UNMERGED                      operator1=domain_prec(:,ispin),&
!FROM UNMERGED                      dpattern=m_quench_t(ispin),&
!FROM UNMERGED                      map=domain_map(ispin),&
!FROM UNMERGED                      node_of_domain=cpu_of_domain(:),&
!FROM UNMERGED                      my_action=0,&
!FROM UNMERGED                      filter_eps=eps_filter,&
!FROM UNMERGED                      !matrix_trimmer=,&
!FROM UNMERGED                      !use_trimmer=.FALSE.,&
!FROM UNMERGED                      error=error)
!FROM UNMERGED 
!FROM UNMERGED           ENDIF ! special_case
!FROM UNMERGED    
!FROM UNMERGED        ELSE ! do not use preconditioner
!FROM UNMERGED 
!FROM UNMERGED           CALL cp_dbcsr_copy(m_zet(ispin),m_residue(ispin),error=error)
!FROM UNMERGED 
!FROM UNMERGED        ENDIF ! use_preconditioner
!FROM UNMERGED 
!FROM UNMERGED ENDIF ! do_exact_inversion
!FROM UNMERGED 
!FROM UNMERGED        CALL cp_dbcsr_copy(m_step(ispin),m_zet(ispin),error=error)
!FROM UNMERGED        
!FROM UNMERGED     ENDDO !ispin
!FROM UNMERGED 
!FROM UNMERGED 
!FROM UNMERGED     ! start the outer SCF loop
!FROM UNMERGED     outer_prepare_to_exit=.FALSE.
!FROM UNMERGED     outer_iteration=0
!FROM UNMERGED     residue_norm=0.0_dp
!FROM UNMERGED 
!FROM UNMERGED     DO
!FROM UNMERGED 
!FROM UNMERGED        ! start the inner SCF loop
!FROM UNMERGED        prepare_to_exit=.FALSE.
!FROM UNMERGED        converged=.FALSE.
!FROM UNMERGED        iteration=0
!FROM UNMERGED        t1 = m_walltime()
!FROM UNMERGED 
!FROM UNMERGED        DO
!FROM UNMERGED 
!FROM UNMERGED           ! apply hessian to the step matrix
!FROM UNMERGED           CALL apply_hessian(&
!FROM UNMERGED                   m_x_in=m_step,&
!FROM UNMERGED                   m_x_out=m_Hstep,&
!FROM UNMERGED                   m_ks=m_ks,&
!FROM UNMERGED                   m_s=m_s,&
!FROM UNMERGED                   m_siginv=m_siginv,&
!FROM UNMERGED                   m_quench_t=m_quench_t,&
!FROM UNMERGED                   m_FTsiginv=m_FTsiginv,&
!FROM UNMERGED                   m_siginvTFTsiginv=m_siginvTFTsiginv,&
!FROM UNMERGED                   m_ST=m_ST,&
!FROM UNMERGED                   m_STsiginv=m_STsiginv,&
!FROM UNMERGED                   m_s_vv=m_s_vv,&
!FROM UNMERGED                   m_ks_vv=m_f_vv,&
!FROM UNMERGED                   !m_s_vv=m_s,&
!FROM UNMERGED                   !m_ks_vv=m_ks,&
!FROM UNMERGED                   m_g_full=m_f_vo,&
!FROM UNMERGED                   m_t=m_t,&
!FROM UNMERGED                   m_sig_sqrti_ii=m_sig_sqrti_ii,&
!FROM UNMERGED                   penalty_occ_vol=penalty_occ_vol,&
!FROM UNMERGED                   normalize_orbitals=normalize_orbitals,&
!FROM UNMERGED                   penalty_occ_vol_prefactor=penalty_occ_vol_prefactor,&
!FROM UNMERGED                   eps_filter=eps_filter,&
!FROM UNMERGED                   path_num=hessian_path_reuse,&
!FROM UNMERGED                   error=error)
!FROM UNMERGED 
!FROM UNMERGED           ! alpha is computed outside the spin loop
!FROM UNMERGED           numer = 0.0_dp
!FROM UNMERGED           denom = 0.0_dp
!FROM UNMERGED           DO ispin=1,nspins
!FROM UNMERGED 
!FROM UNMERGED              CALL cp_dbcsr_trace(m_residue(ispin),m_zet(ispin),numer_ispin,&
!FROM UNMERGED                      "T","N",error=error)
!FROM UNMERGED              CALL cp_dbcsr_trace(m_step(ispin),m_Hstep(ispin),denom_ispin,&
!FROM UNMERGED                      "T","N",error=error)
!FROM UNMERGED 
!FROM UNMERGED              numer = numer + numer_ispin
!FROM UNMERGED              denom = denom + denom_ispin
!FROM UNMERGED 
!FROM UNMERGED           ENDDO !ispin
!FROM UNMERGED           
!FROM UNMERGED           alpha = numer/denom
!FROM UNMERGED 
!FROM UNMERGED           DO ispin=1,nspins
!FROM UNMERGED          
!FROM UNMERGED              ! update the variable
!FROM UNMERGED              CALL cp_dbcsr_add(m_delta(ispin),m_step(ispin),1.0_dp,alpha,&
!FROM UNMERGED                      error=error)
!FROM UNMERGED              CALL cp_dbcsr_copy(m_residue_prev(ispin),m_residue(ispin),&
!FROM UNMERGED                      error=error)
!FROM UNMERGED              CALL cp_dbcsr_add(m_residue(ispin),m_Hstep(ispin),&
!FROM UNMERGED                      1.0_dp, -1.0_dp*alpha, error=error)
!FROM UNMERGED           
!FROM UNMERGED              CALL cp_dbcsr_norm(m_residue(ispin), dbcsr_norm_maxabsnorm,&
!FROM UNMERGED                   norm_scalar=residue_max_norm(ispin), error=error)
!FROM UNMERGED 
!FROM UNMERGED           ENDDO ! ispin
!FROM UNMERGED 
!FROM UNMERGED           ! check convergence and other exit criteria
!FROM UNMERGED           residue_norm=MAXVAL(residue_max_norm)
!FROM UNMERGED           converged=(residue_norm.lt.eps_error_target)
!FROM UNMERGED           IF (converged.OR.(iteration.ge.max_iter)) THEN
!FROM UNMERGED              prepare_to_exit=.TRUE.
!FROM UNMERGED           ENDIF
!FROM UNMERGED 
!FROM UNMERGED           IF (.NOT.prepare_to_exit) THEN
!FROM UNMERGED     
!FROM UNMERGED              denom=0.0_dp
!FROM UNMERGED              numer=0.0_dp
!FROM UNMERGED              DO ispin=1,nspins
!FROM UNMERGED                 
!FROM UNMERGED                 ! save current z before the update
!FROM UNMERGED                 CALL cp_dbcsr_copy(m_zet_prev(ispin),m_zet(ispin),&
!FROM UNMERGED                         error=error)
!FROM UNMERGED     
!FROM UNMERGED                 ! compute the new step (apply preconditioner if available)
!FROM UNMERGED                 IF (use_preconditioner) THEN
!FROM UNMERGED                    
!FROM UNMERGED                    !IF (unit_nr>0) THEN
!FROM UNMERGED                    !   WRITE(unit_nr,*) "....applying preconditioner...."
!FROM UNMERGED                    !ENDIF
!FROM UNMERGED 
!FROM UNMERGED                    IF (special_case.eq.xalmo_case_block_diag .OR. &
!FROM UNMERGED                        special_case.eq.xalmo_case_fully_deloc) THEN
!FROM UNMERGED 
!FROM UNMERGED                       CALL cp_dbcsr_multiply("N","N",1.0_dp,&
!FROM UNMERGED                               m_prec(ispin),&
!FROM UNMERGED                               m_residue(ispin),&
!FROM UNMERGED                               0.0_dp,m_zet(ispin),&
!FROM UNMERGED                               filter_eps=eps_filter,&
!FROM UNMERGED                               error=error)
!FROM UNMERGED 
!FROM UNMERGED                    ELSE
!FROM UNMERGED                       
!FROM UNMERGED                       CALL apply_domain_operators(&
!FROM UNMERGED                               matrix_in=m_residue(ispin),&
!FROM UNMERGED                               matrix_out=m_zet(ispin),&
!FROM UNMERGED                               operator1=domain_prec(:,ispin),&
!FROM UNMERGED                               dpattern=m_quench_t(ispin),&
!FROM UNMERGED                               map=domain_map(ispin),&
!FROM UNMERGED                               node_of_domain=cpu_of_domain(:),&
!FROM UNMERGED                               my_action=0,&
!FROM UNMERGED                               filter_eps=eps_filter,&
!FROM UNMERGED                               !matrix_trimmer=,&
!FROM UNMERGED                               !use_trimmer=.FALSE.,&
!FROM UNMERGED                               error=error)
!FROM UNMERGED    
!FROM UNMERGED                    ENDIF ! special case
!FROM UNMERGED           
!FROM UNMERGED                 ELSE
!FROM UNMERGED 
!FROM UNMERGED                    CALL cp_dbcsr_copy(m_zet(ispin),m_residue(ispin),error=error)
!FROM UNMERGED 
!FROM UNMERGED                 ENDIF
!FROM UNMERGED              
!FROM UNMERGED                 ! compute the conjugation coefficient - beta
!FROM UNMERGED                 CALL compute_cg_beta(&
!FROM UNMERGED                         beta=beta_ispin,&
!FROM UNMERGED                         numer=numer_ispin,&
!FROM UNMERGED                         denom=denom_ispin,&
!FROM UNMERGED                         reset_conjugator=reset_conjugator,&
!FROM UNMERGED                         conjugator=cg_fletcher,&
!FROM UNMERGED                         grad=m_residue(ispin),&
!FROM UNMERGED                         prev_grad=m_residue_prev(ispin),&
!FROM UNMERGED                         step=m_zet(ispin),&
!FROM UNMERGED                         prev_step=m_zet_prev(ispin),&
!FROM UNMERGED                         error=error)
!FROM UNMERGED              
!FROM UNMERGED                 numer = numer + numer_ispin
!FROM UNMERGED                 denom = denom + denom_ispin
!FROM UNMERGED 
!FROM UNMERGED              ENDDO !ispin
!FROM UNMERGED    
!FROM UNMERGED              ! compute beta from ispin values
!FROM UNMERGED              beta = numer/denom
!FROM UNMERGED 
!FROM UNMERGED              DO ispin=1,nspins
!FROM UNMERGED              
!FROM UNMERGED                 ! conjugate the step direction
!FROM UNMERGED                 CALL cp_dbcsr_add(m_step(ispin),m_zet(ispin),beta,1.0_dp,error=error)
!FROM UNMERGED              
!FROM UNMERGED              ENDDO !ispin
!FROM UNMERGED 
!FROM UNMERGED           ENDIF ! not.prepare_to_exit
!FROM UNMERGED           
!FROM UNMERGED           t2 = m_walltime()
!FROM UNMERGED           IF (unit_nr>0) THEN
!FROM UNMERGED              !iter_type=TRIM("ALMO SCF "//iter_type)
!FROM UNMERGED              iter_type=TRIM("NR STEP")
!FROM UNMERGED              WRITE(unit_nr,'(T6,A9,I6,F14.5,F14.5,F15.10,F9.2)') &
!FROM UNMERGED                      iter_type,iteration,&
!FROM UNMERGED                      alpha,beta,residue_norm,&
!FROM UNMERGED                      t2-t1
!FROM UNMERGED           ENDIF
!FROM UNMERGED           t1 = m_walltime()
!FROM UNMERGED 
!FROM UNMERGED           iteration=iteration+1
!FROM UNMERGED           IF (prepare_to_exit) EXIT
!FROM UNMERGED 
!FROM UNMERGED        ENDDO ! inner loop
!FROM UNMERGED 
!FROM UNMERGED        IF (converged.OR.(outer_iteration.ge.outer_max_iter)) THEN
!FROM UNMERGED           outer_prepare_to_exit=.TRUE.
!FROM UNMERGED        ENDIF
!FROM UNMERGED        
!FROM UNMERGED        outer_iteration=outer_iteration+1
!FROM UNMERGED        IF (outer_prepare_to_exit) EXIT
!FROM UNMERGED 
!FROM UNMERGED     ENDDO ! outer loop
!FROM UNMERGED 
!FROM UNMERGED ! is not necessary if penalty_occ_vol_pf2=0.0
!FROM UNMERGED #if 0
!FROM UNMERGED 
!FROM UNMERGED IF (penalty_occ_vol) THEN
!FROM UNMERGED 
!FROM UNMERGED    DO ispin=1, nspins
!FROM UNMERGED    
!FROM UNMERGED       CALL cp_dbcsr_copy(m_zet(ispin),m_grad(ispin),error=error)
!FROM UNMERGED       CALL cp_dbcsr_trace(m_delta(ispin),m_zet(ispin),alpha,&
!FROM UNMERGED               "T","N",error=error)
!FROM UNMERGED       WRITE(*,*) "trace(grad.delta): ", alpha
!FROM UNMERGED       alpha=-1.0_dp/(penalty_occ_vol_pf2*alpha-1.0_dp)
!FROM UNMERGED       WRITE(*,*) "correction alpha: ", alpha
!FROM UNMERGED       CALL cp_dbcsr_scale(m_delta(ispin),alpha,error=error)
!FROM UNMERGED 
!FROM UNMERGED    ENDDO
!FROM UNMERGED 
!FROM UNMERGED ENDIF
!FROM UNMERGED 
!FROM UNMERGED #endif
!FROM UNMERGED 
!FROM UNMERGED     DO ispin=1, nspins
!FROM UNMERGED        
!FROM UNMERGED        ! check whether the step lies entirely in R or Q
!FROM UNMERGED        CALL cp_dbcsr_init(m_tmp_oo_1,error=error)
!FROM UNMERGED        CALL cp_dbcsr_create(m_tmp_oo_1,&
!FROM UNMERGED                template=m_siginv(ispin),&
!FROM UNMERGED                matrix_type=dbcsr_type_no_symmetry,error=error)
!FROM UNMERGED        CALL cp_dbcsr_init(m_tmp_oo_2,error=error)
!FROM UNMERGED        CALL cp_dbcsr_create(m_tmp_oo_2,&
!FROM UNMERGED                template=m_siginv(ispin),&
!FROM UNMERGED                matrix_type=dbcsr_type_no_symmetry,error=error)
!FROM UNMERGED        CALL cp_dbcsr_multiply("T","N",1.0_dp,&
!FROM UNMERGED                m_ST(ispin),&
!FROM UNMERGED                m_delta(ispin),&
!FROM UNMERGED                0.0_dp,m_tmp_oo_1,&
!FROM UNMERGED                filter_eps=eps_filter,&
!FROM UNMERGED                error=error)
!FROM UNMERGED        CALL cp_dbcsr_multiply("N","N",1.0_dp,&
!FROM UNMERGED                m_siginv(ispin),&
!FROM UNMERGED                m_tmp_oo_1,&
!FROM UNMERGED                0.0_dp,m_tmp_oo_2,&
!FROM UNMERGED                filter_eps=eps_filter,&
!FROM UNMERGED                error=error)
!FROM UNMERGED        CALL cp_dbcsr_copy(m_zet(ispin),m_quench_t(ispin),error=error)
!FROM UNMERGED        CALL cp_dbcsr_multiply("N","N",1.0_dp,&
!FROM UNMERGED                m_t(ispin),&
!FROM UNMERGED                m_tmp_oo_2,&
!FROM UNMERGED                0.0_dp,m_zet(ispin),&
!FROM UNMERGED                retain_sparsity=.TRUE.,&
!FROM UNMERGED                error=error)
!FROM UNMERGED        CALL cp_dbcsr_norm(m_zet(ispin), dbcsr_norm_maxabsnorm,&
!FROM UNMERGED                norm_scalar=alpha, error=error)
!FROM UNMERGED        WRITE(*,"(A50,2F20.10)") "Occupied-space projection of the step", alpha
!FROM UNMERGED        CALL cp_dbcsr_add(m_zet(ispin),m_delta(ispin),-1.0_dp,1.0_dp,error=error)
!FROM UNMERGED        CALL cp_dbcsr_norm(m_zet(ispin), dbcsr_norm_maxabsnorm,&
!FROM UNMERGED                norm_scalar=alpha, error=error)
!FROM UNMERGED        WRITE(*,"(A50,2F20.10)") "Virtual-space projection of the step", alpha
!FROM UNMERGED        CALL cp_dbcsr_norm(m_delta(ispin), dbcsr_norm_maxabsnorm,&
!FROM UNMERGED                norm_scalar=alpha, error=error)
!FROM UNMERGED        WRITE(*,"(A50,2F20.10)") "Full step", alpha 
!FROM UNMERGED        CALL cp_dbcsr_release(m_tmp_oo_1,error=error)
!FROM UNMERGED        CALL cp_dbcsr_release(m_tmp_oo_2,error=error)
!FROM UNMERGED 
!FROM UNMERGED     ENDDO
!FROM UNMERGED 
!FROM UNMERGED     ! clean up
!FROM UNMERGED     DO ispin=1,nspins
!FROM UNMERGED        CALL release_submatrices(domain_prec(:,ispin),error=error)
!FROM UNMERGED        CALL cp_dbcsr_release(m_residue(ispin),error=error)
!FROM UNMERGED        CALL cp_dbcsr_release(m_residue_prev(ispin),error=error)
!FROM UNMERGED        CALL cp_dbcsr_release(m_step(ispin),error=error)
!FROM UNMERGED        CALL cp_dbcsr_release(m_zet(ispin),error=error)
!FROM UNMERGED        CALL cp_dbcsr_release(m_zet_prev(ispin),error=error)
!FROM UNMERGED        CALL cp_dbcsr_release(m_Hstep(ispin),error=error)
!FROM UNMERGED        CALL cp_dbcsr_release(m_f_vo(ispin),error=error)
!FROM UNMERGED        CALL cp_dbcsr_release(m_f_vv(ispin),error=error)
!FROM UNMERGED        CALL cp_dbcsr_release(m_s_vv(ispin),error=error)
!FROM UNMERGED        CALL cp_dbcsr_release(m_prec(ispin),error=error)
!FROM UNMERGED        CALL cp_dbcsr_release(m_STsiginv(ispin),error=error)
!FROM UNMERGED     ENDDO !ispin
!FROM UNMERGED     DEALLOCATE(domain_prec)
!FROM UNMERGED     DEALLOCATE(m_residue)
!FROM UNMERGED     DEALLOCATE(m_residue_prev)
!FROM UNMERGED     DEALLOCATE(m_step)
!FROM UNMERGED     DEALLOCATE(m_zet)
!FROM UNMERGED     DEALLOCATE(m_zet_prev)
!FROM UNMERGED     DEALLOCATE(m_prec)
!FROM UNMERGED     DEALLOCATE(m_Hstep)
!FROM UNMERGED     DEALLOCATE(m_s_vv)
!FROM UNMERGED     DEALLOCATE(m_f_vv)
!FROM UNMERGED     DEALLOCATE(m_f_vo)
!FROM UNMERGED     DEALLOCATE(m_STsiginv)
!FROM UNMERGED     DEALLOCATE(residue_max_norm)
!FROM UNMERGED     
!FROM UNMERGED     IF (.NOT.converged)  THEN
!FROM UNMERGED        IF (unit_nr>0) WRITE(unit_nr,'(T2,A)') "Optimization not converged! "
!FROM UNMERGED        CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
!FROM UNMERGED     ENDIF
!FROM UNMERGED     
!FROM UNMERGED     ! check that the step satisfies H.step=-grad
!FROM UNMERGED 
!FROM UNMERGED     CALL timestop(handle)
!FROM UNMERGED   
!FROM UNMERGED   END SUBROUTINE newton_grad_to_step
!FROM UNMERGED 
!FROM UNMERGED ! *****************************************************************************
!FROM UNMERGED !> \brief Computes Hessian.X 
!FROM UNMERGED !> \par History
!FROM UNMERGED !>       2015.04 created [Rustam Z Khaliullin]
!FROM UNMERGED !> \author Rustam Z Khaliullin
!FROM UNMERGED ! *****************************************************************************
!FROM UNMERGED   SUBROUTINE apply_hessian(m_x_in, m_x_out, m_ks, m_s, m_siginv, &
!FROM UNMERGED      m_quench_t,  m_FTsiginv, m_siginvTFTsiginv, m_ST, m_STsiginv, m_s_vv, &
!FROM UNMERGED      m_ks_vv, m_g_full, m_t, m_sig_sqrti_ii, penalty_occ_vol, &
!FROM UNMERGED      normalize_orbitals, penalty_occ_vol_prefactor, eps_filter, path_num, error)
!FROM UNMERGED 
!FROM UNMERGED     TYPE(cp_dbcsr_type), DIMENSION(:), INTENT(INOUT) :: m_x_in
!FROM UNMERGED     TYPE(cp_dbcsr_type), DIMENSION(:), INTENT(INOUT) :: m_x_out
!FROM UNMERGED     TYPE(cp_dbcsr_type), DIMENSION(:), INTENT(INOUT) :: m_ks
!FROM UNMERGED     TYPE(cp_dbcsr_type), DIMENSION(:), INTENT(INOUT) :: m_s
!FROM UNMERGED     TYPE(cp_dbcsr_type), DIMENSION(:), INTENT(IN)    :: m_siginv
!FROM UNMERGED     TYPE(cp_dbcsr_type), DIMENSION(:), INTENT(INOUT) :: m_quench_t
!FROM UNMERGED     TYPE(cp_dbcsr_type), DIMENSION(:), INTENT(IN)    :: m_FTsiginv
!FROM UNMERGED     TYPE(cp_dbcsr_type), DIMENSION(:), INTENT(IN)    :: m_siginvTFTsiginv
!FROM UNMERGED     TYPE(cp_dbcsr_type), DIMENSION(:), INTENT(IN)    :: m_ST
!FROM UNMERGED     TYPE(cp_dbcsr_type), DIMENSION(:), INTENT(IN)    :: m_STsiginv
!FROM UNMERGED     TYPE(cp_dbcsr_type), DIMENSION(:), INTENT(INOUT) :: m_s_vv
!FROM UNMERGED     TYPE(cp_dbcsr_type), DIMENSION(:), INTENT(INOUT) :: m_ks_vv
!FROM UNMERGED     TYPE(cp_dbcsr_type), DIMENSION(:), INTENT(INOUT) :: m_g_full
!FROM UNMERGED     TYPE(cp_dbcsr_type), DIMENSION(:), INTENT(IN)    :: m_t
!FROM UNMERGED     TYPE(cp_dbcsr_type), DIMENSION(:), INTENT(IN)    :: m_sig_sqrti_ii
!FROM UNMERGED     LOGICAL, INTENT(IN)                              :: penalty_occ_vol, &
!FROM UNMERGED                                                         normalize_orbitals
!FROM UNMERGED     REAL(KIND=dp), INTENT(IN)                        :: penalty_occ_vol_prefactor
!FROM UNMERGED     REAL(KIND=dp), INTENT(IN)                        :: eps_filter
!FROM UNMERGED     INTEGER, INTENT(IN)                              :: path_num 
!FROM UNMERGED     TYPE(cp_error_type), INTENT(INOUT)               :: error
!FROM UNMERGED     
!FROM UNMERGED     CHARACTER(len=*), PARAMETER :: &
!FROM UNMERGED       routineN = 'apply_hessian', &
!FROM UNMERGED       routineP = moduleN//':'//routineN
!FROM UNMERGED 
!FROM UNMERGED     INTEGER                                     :: handle,nspins,&
!FROM UNMERGED                                                    ispin, dim0
!FROM UNMERGED     REAL(KIND=dp)                               :: spin_factor,&
!FROM UNMERGED                                                    penalty_prefactor_local
!FROM UNMERGED     REAL(KIND=dp), DIMENSION(:), ALLOCATABLE    :: tg_diagonal
!FROM UNMERGED     LOGICAL                                     :: failure
!FROM UNMERGED     TYPE(cp_dbcsr_type)                         :: m_tmp_oo_1,&
!FROM UNMERGED                                                    m_tmp_no_1,&
!FROM UNMERGED                                                    m_tmp_no_2,&
!FROM UNMERGED                                                    m_tmp_x_in
!FROM UNMERGED 
!FROM UNMERGED     CALL timeset(routineN,handle)
!FROM UNMERGED 
!FROM UNMERGED     nspins=SIZE(m_ks)
!FROM UNMERGED     
!FROM UNMERGED     IF (nspins.eq.1) THEN
!FROM UNMERGED        spin_factor = 2.0_dp
!FROM UNMERGED     ELSE
!FROM UNMERGED        spin_factor = 1.0_dp
!FROM UNMERGED     ENDIF
!FROM UNMERGED 
!FROM UNMERGED     penalty_prefactor_local=penalty_occ_vol_prefactor/(2.0_dp*spin_factor)
!FROM UNMERGED 
!FROM UNMERGED     DO ispin=1,nspins
!FROM UNMERGED 
!FROM UNMERGED        CALL cp_dbcsr_init(m_tmp_oo_1,error=error)
!FROM UNMERGED        CALL cp_dbcsr_create(m_tmp_oo_1,&
!FROM UNMERGED                template=m_siginv(ispin),&
!FROM UNMERGED                matrix_type=dbcsr_type_no_symmetry,&
!FROM UNMERGED                error=error)
!FROM UNMERGED        CALL cp_dbcsr_init(m_tmp_no_1,error=error)
!FROM UNMERGED        CALL cp_dbcsr_create(m_tmp_no_1,&
!FROM UNMERGED                template=m_quench_t(ispin),&
!FROM UNMERGED                matrix_type=dbcsr_type_no_symmetry,&
!FROM UNMERGED                error=error)
!FROM UNMERGED        CALL cp_dbcsr_init(m_tmp_no_2,error=error)
!FROM UNMERGED        CALL cp_dbcsr_create(m_tmp_no_2,&
!FROM UNMERGED                template=m_quench_t(ispin),&
!FROM UNMERGED                matrix_type=dbcsr_type_no_symmetry,&
!FROM UNMERGED                error=error)
!FROM UNMERGED        CALL cp_dbcsr_init(m_tmp_x_in,error=error)
!FROM UNMERGED        CALL cp_dbcsr_create(m_tmp_x_in,&
!FROM UNMERGED                template=m_quench_t(ispin),&
!FROM UNMERGED                matrix_type=dbcsr_type_no_symmetry,&
!FROM UNMERGED                error=error)
!FROM UNMERGED 
!FROM UNMERGED        ! transform the input X to take into account the normalization constraint
!FROM UNMERGED        IF (normalize_orbitals) THEN
!FROM UNMERGED           
!FROM UNMERGED           ! H.D = ( (H.D) - ST.[tr(T).(H.D)]_ii ) . [sig_sqrti]_ii
!FROM UNMERGED           
!FROM UNMERGED           ! get [tr(T).HD]_ii
!FROM UNMERGED           CALL cp_dbcsr_copy(m_tmp_oo_1,m_sig_sqrti_ii(ispin),error=error)
!FROM UNMERGED           CALL cp_dbcsr_multiply("T","N",1.0_dp,&
!FROM UNMERGED                   m_x_in(ispin),&
!FROM UNMERGED                   m_ST(ispin),&
!FROM UNMERGED                   0.0_dp,m_tmp_oo_1,&
!FROM UNMERGED                   retain_sparsity=.TRUE.,&
!FROM UNMERGED                   error=error)
!FROM UNMERGED           CALL cp_dbcsr_get_info(m_sig_sqrti_ii(ispin), nfullrows_total=dim0 )
!FROM UNMERGED           ALLOCATE(tg_diagonal(dim0))
!FROM UNMERGED           CALL cp_dbcsr_get_diag(m_tmp_oo_1,tg_diagonal,error=error)
!FROM UNMERGED           CALL cp_dbcsr_set(m_tmp_oo_1,0.0_dp,error=error)
!FROM UNMERGED           CALL cp_dbcsr_set_diag(m_tmp_oo_1,tg_diagonal,error=error)
!FROM UNMERGED           DEALLOCATE(tg_diagonal)
!FROM UNMERGED           
!FROM UNMERGED           CALL cp_dbcsr_copy(m_tmp_no_1,m_x_in(ispin),error=error)
!FROM UNMERGED           CALL cp_dbcsr_multiply("N","N",-1.0_dp,&
!FROM UNMERGED                   m_t(ispin),&
!FROM UNMERGED                   m_tmp_oo_1,&
!FROM UNMERGED                   1.0_dp,m_tmp_no_1,&
!FROM UNMERGED                   filter_eps=eps_filter,&
!FROM UNMERGED                   error=error)
!FROM UNMERGED           CALL cp_dbcsr_multiply("N","N",1.0_dp,&
!FROM UNMERGED                   m_tmp_no_1,&
!FROM UNMERGED                   m_sig_sqrti_ii(ispin),&
!FROM UNMERGED                   0.0_dp,m_tmp_x_in,&
!FROM UNMERGED                   filter_eps=eps_filter,&
!FROM UNMERGED                   error=error)
!FROM UNMERGED        
!FROM UNMERGED        ELSE
!FROM UNMERGED 
!FROM UNMERGED           CALL cp_dbcsr_copy(m_tmp_x_in,m_x_in(ispin),error=error)
!FROM UNMERGED 
!FROM UNMERGED        ENDIF ! normalize_orbitals
!FROM UNMERGED 
!FROM UNMERGED        IF (path_num.EQ.hessian_path_reuse) THEN
!FROM UNMERGED           
!FROM UNMERGED           ! apply pre-computed F_vv and S_vv to X
!FROM UNMERGED 
!FROM UNMERGED #if 0
!FROM UNMERGED ! RZK-warning: negative sign at penalty_prefactor_local is that 
!FROM UNMERGED ! magical fix for the negative definite problem
!FROM UNMERGED ! (since penalty_prefactor_local<0 the coeff before S_vv must
!FROM UNMERGED ! be multiplied by -1 to take the step in the right direction)
!FROM UNMERGED !CALL cp_dbcsr_multiply("N","N",-4.0_dp*penalty_prefactor_local,&
!FROM UNMERGED !        m_s_vv(ispin),&
!FROM UNMERGED !        m_tmp_x_in,&
!FROM UNMERGED !        0.0_dp,m_tmp_no_1,&
!FROM UNMERGED !        filter_eps=eps_filter,&
!FROM UNMERGED !        error=error)
!FROM UNMERGED !CALL cp_dbcsr_copy(m_x_out(ispin),m_quench_t(ispin),error=error)
!FROM UNMERGED !CALL cp_dbcsr_multiply("N","N",1.0_dp,&
!FROM UNMERGED !        m_tmp_no_1,&
!FROM UNMERGED !        m_siginv(ispin),&
!FROM UNMERGED !        0.0_dp,m_x_out(ispin),&
!FROM UNMERGED !        retain_sparsity=.TRUE.,&
!FROM UNMERGED !        error=error)
!FROM UNMERGED 
!FROM UNMERGED CALL cp_dbcsr_multiply("N","N",1.0_dp,&
!FROM UNMERGED         m_s(1),&
!FROM UNMERGED         m_tmp_x_in,&
!FROM UNMERGED         0.0_dp,m_tmp_no_1,&
!FROM UNMERGED         filter_eps=eps_filter,&
!FROM UNMERGED         error=error)
!FROM UNMERGED CALL cp_dbcsr_copy(m_x_out(ispin),m_quench_t(ispin),error=error)
!FROM UNMERGED CALL cp_dbcsr_multiply("N","N",1.0_dp,&
!FROM UNMERGED         m_tmp_no_1,&
!FROM UNMERGED         m_siginv(ispin),&
!FROM UNMERGED         0.0_dp,m_x_out(ispin),&
!FROM UNMERGED         retain_sparsity=.TRUE.,&
!FROM UNMERGED         error=error)
!FROM UNMERGED 
!FROM UNMERGED !CALL cp_dbcsr_copy(m_x_out(ispin),m_quench_t(ispin),error=error)
!FROM UNMERGED !CALL cp_dbcsr_multiply("N","N",1.0_dp,&
!FROM UNMERGED !        m_s(1),&
!FROM UNMERGED !        m_tmp_x_in,&
!FROM UNMERGED !        0.0_dp,m_x_out(ispin),&
!FROM UNMERGED !        retain_sparsity=.TRUE.,&
!FROM UNMERGED !        error=error)
!FROM UNMERGED        
!FROM UNMERGED        
!FROM UNMERGED #else
!FROM UNMERGED 
!FROM UNMERGED           ! debugging: only vv matrices, oo matrices are kronecker
!FROM UNMERGED           CALL cp_dbcsr_copy(m_x_out(ispin),m_quench_t(ispin),error=error)
!FROM UNMERGED           CALL cp_dbcsr_multiply("N","N",1.0_dp,&
!FROM UNMERGED                   m_ks_vv(ispin),&
!FROM UNMERGED                   m_tmp_x_in,&
!FROM UNMERGED                   0.0_dp,m_x_out(ispin),&
!FROM UNMERGED                   retain_sparsity=.TRUE.,&
!FROM UNMERGED                   error=error)
!FROM UNMERGED        
!FROM UNMERGED           CALL cp_dbcsr_copy(m_tmp_no_2,m_quench_t(ispin),error=error)
!FROM UNMERGED           CALL cp_dbcsr_multiply("N","N",1.0_dp,&
!FROM UNMERGED                   m_s_vv(ispin),&
!FROM UNMERGED                   m_tmp_x_in,&
!FROM UNMERGED                   0.0_dp,m_tmp_no_2,&
!FROM UNMERGED                   retain_sparsity=.TRUE.,&
!FROM UNMERGED                   error=error)
!FROM UNMERGED           CALL cp_dbcsr_add(m_x_out(ispin),m_tmp_no_2,&
!FROM UNMERGED                   1.0_dp,-4.0_dp*penalty_prefactor_local+1.0_dp,error=error)
!FROM UNMERGED #endif
!FROM UNMERGED 
!FROM UNMERGED !          ! F_vv.X.S_oo
!FROM UNMERGED !          CALL cp_dbcsr_multiply("N","N",1.0_dp,&
!FROM UNMERGED !                  m_ks_vv(ispin),&
!FROM UNMERGED !                  m_tmp_x_in,&
!FROM UNMERGED !                  0.0_dp,m_tmp_no_1,&
!FROM UNMERGED !                  filter_eps=eps_filter,&
!FROM UNMERGED !                  error=error)
!FROM UNMERGED !          CALL cp_dbcsr_copy(m_x_out(ispin),m_quench_t(ispin),error=error)
!FROM UNMERGED !          CALL cp_dbcsr_multiply("N","N",1.0_dp,&
!FROM UNMERGED !                  m_tmp_no_1,&
!FROM UNMERGED !                  m_siginv(ispin),&
!FROM UNMERGED !                  0.0_dp,m_x_out(ispin),&
!FROM UNMERGED !                  retain_sparsity=.TRUE.,&
!FROM UNMERGED !                  error=error)
!FROM UNMERGED !       
!FROM UNMERGED !          ! S_vv.X.F_oo
!FROM UNMERGED !          CALL cp_dbcsr_multiply("N","N",1.0_dp,&
!FROM UNMERGED !                  m_s_vv(ispin),&
!FROM UNMERGED !                  m_tmp_x_in,&
!FROM UNMERGED !                  0.0_dp,m_tmp_no_1,&
!FROM UNMERGED !                  filter_eps=eps_filter,&
!FROM UNMERGED !                  error=error)
!FROM UNMERGED !          CALL cp_dbcsr_copy(m_tmp_no_2,m_quench_t(ispin),error=error)
!FROM UNMERGED !          CALL cp_dbcsr_multiply("N","N",1.0_dp,&
!FROM UNMERGED !                  m_tmp_no_1,&
!FROM UNMERGED !                  m_siginvTFTsiginv(ispin),&
!FROM UNMERGED !                  0.0_dp,m_tmp_no_2,&
!FROM UNMERGED !                  retain_sparsity=.TRUE.,&
!FROM UNMERGED !                  error=error)
!FROM UNMERGED !          CALL cp_dbcsr_add(m_x_out(ispin),m_tmp_no_2,&
!FROM UNMERGED !               1.0_dp,-1.0_dp,error=error)
!FROM UNMERGED !! we have to add occ voll penalty here (the Svv termi (i.e. both Svv.D.Soo) 
!FROM UNMERGED !!  and STsiginv terms)
!FROM UNMERGED !
!FROM UNMERGED !         ! S_vo.X^t.F_vo
!FROM UNMERGED !          CALL cp_dbcsr_multiply("T","N",1.0_dp,&
!FROM UNMERGED !                  m_tmp_x_in,&
!FROM UNMERGED !                  m_g_full(ispin),&
!FROM UNMERGED !                  0.0_dp,m_tmp_oo_1,&
!FROM UNMERGED !                  filter_eps=eps_filter,&
!FROM UNMERGED !                  error=error)
!FROM UNMERGED !          CALL cp_dbcsr_copy(m_tmp_no_2,m_quench_t(ispin),error=error)
!FROM UNMERGED !          CALL cp_dbcsr_multiply("N","N",1.0_dp,&
!FROM UNMERGED !                  m_STsiginv(ispin),&
!FROM UNMERGED !                  m_tmp_oo_1,&
!FROM UNMERGED !                  0.0_dp,m_tmp_no_2,&
!FROM UNMERGED !                  retain_sparsity=.TRUE.,&
!FROM UNMERGED !                  error=error)
!FROM UNMERGED !          CALL cp_dbcsr_add(m_x_out(ispin),m_tmp_no_2,&
!FROM UNMERGED !                  1.0_dp,-1.0_dp,error=error)
!FROM UNMERGED !
!FROM UNMERGED !          ! S_vo.X^t.F_vo
!FROM UNMERGED !          CALL cp_dbcsr_multiply("T","N",1.0_dp,&
!FROM UNMERGED !                  m_tmp_x_in,&
!FROM UNMERGED !                  m_STsiginv(ispin),&
!FROM UNMERGED !                  0.0_dp,m_tmp_oo_1,&
!FROM UNMERGED !                  filter_eps=eps_filter,&
!FROM UNMERGED !                  error=error)
!FROM UNMERGED !          CALL cp_dbcsr_copy(m_tmp_no_2,m_quench_t(ispin),error=error)
!FROM UNMERGED !          CALL cp_dbcsr_multiply("N","N",1.0_dp,&
!FROM UNMERGED !                  m_g_full(ispin),&
!FROM UNMERGED !                  m_tmp_oo_1,&
!FROM UNMERGED !                  0.0_dp,m_tmp_no_2,&
!FROM UNMERGED !                  retain_sparsity=.TRUE.,&
!FROM UNMERGED !                  error=error)
!FROM UNMERGED !          CALL cp_dbcsr_add(m_x_out(ispin),m_tmp_no_2,&
!FROM UNMERGED !                  1.0_dp,-1.0_dp,error=error)
!FROM UNMERGED     
!FROM UNMERGED        ELSE IF (path_num.EQ.hessian_path_assemble) THEN
!FROM UNMERGED           
!FROM UNMERGED           ! compute F_vv.X and S_vv.X directly
!FROM UNMERGED           ! this path will be advantageous if the number
!FROM UNMERGED           ! of PCG iterations is small
!FROM UNMERGED           
!FROM UNMERGED           CPErrorMessage(cp_failure_level,routineP,"path is NYI",error)
!FROM UNMERGED           CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
!FROM UNMERGED 
!FROM UNMERGED        ELSE
!FROM UNMERGED           CPErrorMessage(cp_failure_level,routineP,"illegal path",error)
!FROM UNMERGED           CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
!FROM UNMERGED        ENDIF ! path
!FROM UNMERGED 
!FROM UNMERGED        ! transform the output to take into account the normalization constraint
!FROM UNMERGED        IF (normalize_orbitals) THEN
!FROM UNMERGED           
!FROM UNMERGED           ! H.D = ( (H.D) - ST.[tr(T).(H.D)]_ii ) . [sig_sqrti]_ii
!FROM UNMERGED           
!FROM UNMERGED           ! get [tr(T).HD]_ii
!FROM UNMERGED           CALL cp_dbcsr_copy(m_tmp_oo_1,m_sig_sqrti_ii(ispin),error=error)
!FROM UNMERGED           CALL cp_dbcsr_multiply("T","N",1.0_dp,&
!FROM UNMERGED                   m_t(ispin),&
!FROM UNMERGED                   m_x_out(ispin),&
!FROM UNMERGED                   0.0_dp,m_tmp_oo_1,&
!FROM UNMERGED                   retain_sparsity=.TRUE.,&
!FROM UNMERGED                   error=error)
!FROM UNMERGED           CALL cp_dbcsr_get_info(m_sig_sqrti_ii(ispin), nfullrows_total=dim0 )
!FROM UNMERGED           ALLOCATE(tg_diagonal(dim0))
!FROM UNMERGED           CALL cp_dbcsr_get_diag(m_tmp_oo_1,tg_diagonal,error=error)
!FROM UNMERGED           CALL cp_dbcsr_set(m_tmp_oo_1,0.0_dp,error=error)
!FROM UNMERGED           CALL cp_dbcsr_set_diag(m_tmp_oo_1,tg_diagonal,error=error)
!FROM UNMERGED           DEALLOCATE(tg_diagonal)
!FROM UNMERGED           
!FROM UNMERGED           CALL cp_dbcsr_multiply("N","N",-1.0_dp,&
!FROM UNMERGED                   m_ST(ispin),&
!FROM UNMERGED                   m_tmp_oo_1,&
!FROM UNMERGED                   1.0_dp,m_x_out(ispin),&
!FROM UNMERGED                   retain_sparsity=.TRUE.,&
!FROM UNMERGED                   error=error)
!FROM UNMERGED           CALL cp_dbcsr_copy(m_tmp_no_1,m_x_out(ispin),error=error)
!FROM UNMERGED           CALL cp_dbcsr_multiply("N","N",1.0_dp,&
!FROM UNMERGED                   m_tmp_no_1,&
!FROM UNMERGED                   m_sig_sqrti_ii(ispin),&
!FROM UNMERGED                   0.0_dp,m_x_out(ispin),&
!FROM UNMERGED                   retain_sparsity=.TRUE.,&
!FROM UNMERGED                   error=error)
!FROM UNMERGED           
!FROM UNMERGED        ENDIF ! normalize_orbitals
!FROM UNMERGED 
!FROM UNMERGED        CALL cp_dbcsr_scale(m_x_out(ispin),&
!FROM UNMERGED                2.0_dp*spin_factor,error=error)
!FROM UNMERGED        
!FROM UNMERGED        CALL cp_dbcsr_release(m_tmp_oo_1,error=error)
!FROM UNMERGED        CALL cp_dbcsr_release(m_tmp_no_1,error=error)
!FROM UNMERGED        CALL cp_dbcsr_release(m_tmp_no_2,error=error)
!FROM UNMERGED        CALL cp_dbcsr_release(m_tmp_x_in,error=error)
!FROM UNMERGED 
!FROM UNMERGED     ENDDO !ispin
!FROM UNMERGED 
!FROM UNMERGED     ! there is one more part of the hessian that comes
!FROM UNMERGED     ! from T-dependence of the KS matrix
!FROM UNMERGED     ! it is neglected here
!FROM UNMERGED 
!FROM UNMERGED     CALL timestop(handle)
!FROM UNMERGED   
!FROM UNMERGED   END SUBROUTINE apply_hessian 
!FROM UNMERGED 
!FROM UNMERGED ! *****************************************************************************
!FROM UNMERGED !> \brief Serial code that constructs an approximate Hessian
!FROM UNMERGED !> \param matrix_grad ...
!FROM UNMERGED !> \param matrix_step ...
!FROM UNMERGED !> \param matrix_S_ao ...
!FROM UNMERGED !> \param matrix_F_ao ...
!FROM UNMERGED !> \param matrix_S_mo ...
!FROM UNMERGED !> \param matrix_F_mo ...
!FROM UNMERGED !> \param matrix_S_vo ...
!FROM UNMERGED !> \param matrix_F_vo ...
!FROM UNMERGED !> \param quench_t ...
!FROM UNMERGED !> \param error ...
!FROM UNMERGED !> \par History
!FROM UNMERGED !>       2012.02 created [Rustam Z. Khaliullin]
!FROM UNMERGED !> \author Rustam Z. Khaliullin
!FROM UNMERGED ! *****************************************************************************
!FROM UNMERGED   SUBROUTINE hessian_diag_apply(matrix_grad,matrix_step,matrix_S_ao,&
!FROM UNMERGED     matrix_F_ao,matrix_S_mo,matrix_F_mo,matrix_S_vo,matrix_F_vo,quench_t,&
!FROM UNMERGED     penalty_occ_vol,penalty_occ_vol_prefactor,penalty_occ_vol_pf2,&
!FROM UNMERGED     spin_factor,eps_zero,m_s,para_env,blacs_env,error)
!FROM UNMERGED     
!FROM UNMERGED     TYPE(cp_dbcsr_type), INTENT(INOUT) :: matrix_grad, matrix_step, &
!FROM UNMERGED       matrix_S_ao, matrix_F_ao, matrix_S_mo, matrix_S_vo, &
!FROM UNMERGED       matrix_F_vo, quench_t
!FROM UNMERGED     TYPE(cp_dbcsr_type), INTENT(IN)    :: matrix_F_mo
!FROM UNMERGED     LOGICAL, INTENT(IN)                :: penalty_occ_vol
!FROM UNMERGED     REAL(KIND=dp), INTENT(IN)          :: penalty_occ_vol_prefactor
!FROM UNMERGED     REAL(KIND=dp), INTENT(IN)          :: penalty_occ_vol_pf2
!FROM UNMERGED     REAL(KIND=dp), INTENT(IN)          :: spin_factor,eps_zero
!FROM UNMERGED TYPE(cp_dbcsr_type), INTENT(IN)    :: m_s
!FROM UNMERGED TYPE(cp_para_env_type),  POINTER         :: para_env
!FROM UNMERGED TYPE(cp_blacs_env_type), POINTER         :: blacs_env
!FROM UNMERGED     TYPE(cp_error_type), INTENT(INOUT) :: error
!FROM UNMERGED 
!FROM UNMERGED     CHARACTER(len=*), PARAMETER :: routineN = 'hessian_diag_apply', &
!FROM UNMERGED       routineP = moduleN//':'//routineN
!FROM UNMERGED 
!FROM UNMERGED     INTEGER :: ao_hori_offset, ao_vert_offset, block_col, block_row, col, &
!FROM UNMERGED       copy, H_size, handle, ii, INFO, jj, lev1_hori_offset, lev1_vert_offset, &
!FROM UNMERGED       lev2_hori_offset, lev2_vert_offset, LWORK, nblkcols_tot, nblkrows_tot, &
!FROM UNMERGED       orb_i, orb_j, row, zero_neg_eiv, ncores
!FROM UNMERGED     INTEGER, ALLOCATABLE, DIMENSION(:)       :: ao_block_sizes, &
!FROM UNMERGED                                                 ao_domain_sizes, &
!FROM UNMERGED                                                 mo_block_sizes
!FROM UNMERGED     INTEGER, DIMENSION(:), POINTER           :: ao_blk_sizes, mo_blk_sizes
!FROM UNMERGED     LOGICAL                                  :: failure, found, found2, &
!FROM UNMERGED                                                 found_col, found_row
!FROM UNMERGED     REAL(KIND=dp)                            :: test_error, &
!FROM UNMERGED                                                 penalty_prefactor_local
!FROM UNMERGED     REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: eigenvalues, Grad_vec, &
!FROM UNMERGED                                                 Step_vec, tmp, tmpr, work
!FROM UNMERGED     REAL(KIND=dp), ALLOCATABLE, &
!FROM UNMERGED       DIMENSION(:, :)                        :: F_ao_block, F_mo_block, H, &
!FROM UNMERGED                                                 H1, H2, Hinv, S_ao_block, &
!FROM UNMERGED                                                 S_mo_block, test, test2
!FROM UNMERGED     REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block_p, block_p2, p_new_block
!FROM UNMERGED     TYPE(cp_dbcsr_type)                      :: matrix_F_ao_sym, &
!FROM UNMERGED                                                 matrix_F_mo_sym, &
!FROM UNMERGED                                                 matrix_S_ao_sym, &
!FROM UNMERGED                                                 matrix_S_mo_sym
!FROM UNMERGED TYPE(cp_dbcsr_type)                      :: m_prec_out, m_tmp_no_1
!FROM UNMERGED REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: test3 
!FROM UNMERGED 
!FROM UNMERGED     CALL timeset(routineN,handle)
!FROM UNMERGED     
!FROM UNMERGED     ! serial code only
!FROM UNMERGED     ncores = dbcsr_mp_numnodes(dbcsr_distribution_mp(&
!FROM UNMERGED                 cp_dbcsr_distribution(matrix_S_ao)))
!FROM UNMERGED     IF (ncores.gt.1) THEN                            
!FROM UNMERGED        CPErrorMessage(cp_failure_level,routineP,"serial code only",error)
!FROM UNMERGED        CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
!FROM UNMERGED     ENDIF
!FROM UNMERGED 
!FROM UNMERGED     nblkrows_tot = cp_dbcsr_nblkrows_total(quench_t)
!FROM UNMERGED     nblkcols_tot = cp_dbcsr_nblkcols_total(quench_t)
!FROM UNMERGED     CPPostcondition(nblkrows_tot==nblkcols_tot,cp_failure_level,routineP,error,failure)
!FROM UNMERGED     mo_blk_sizes => array_data(cp_dbcsr_col_block_sizes(quench_t))
!FROM UNMERGED     ao_blk_sizes => array_data(cp_dbcsr_row_block_sizes(quench_t))
!FROM UNMERGED     ALLOCATE(mo_block_sizes(nblkcols_tot),ao_block_sizes(nblkcols_tot))
!FROM UNMERGED     ALLOCATE(ao_domain_sizes(nblkcols_tot))
!FROM UNMERGED     mo_block_sizes(:)=mo_blk_sizes(:)
!FROM UNMERGED     ao_block_sizes(:)=ao_blk_sizes(:)
!FROM UNMERGED     ao_domain_sizes(:)=0
!FROM UNMERGED 
!FROM UNMERGED     CALL cp_dbcsr_init(matrix_S_ao_sym,error=error)
!FROM UNMERGED     CALL cp_dbcsr_create(matrix_S_ao_sym,&
!FROM UNMERGED             template=matrix_S_ao,&
!FROM UNMERGED             matrix_type=dbcsr_type_no_symmetry,error=error)
!FROM UNMERGED     CALL cp_dbcsr_desymmetrize(matrix_S_ao,matrix_S_ao_sym,error=error)
!FROM UNMERGED     CALL cp_dbcsr_scale(matrix_S_ao_sym,2.0_dp*spin_factor,error=error)
!FROM UNMERGED 
!FROM UNMERGED     CALL cp_dbcsr_init(matrix_F_ao_sym,error=error)
!FROM UNMERGED     CALL cp_dbcsr_create(matrix_F_ao_sym,&
!FROM UNMERGED             template=matrix_F_ao,&
!FROM UNMERGED             matrix_type=dbcsr_type_no_symmetry,error=error)
!FROM UNMERGED     CALL cp_dbcsr_desymmetrize(matrix_F_ao,matrix_F_ao_sym,error=error)
!FROM UNMERGED     CALL cp_dbcsr_scale(matrix_F_ao_sym,2.0_dp*spin_factor,error=error)
!FROM UNMERGED 
!FROM UNMERGED     CALL cp_dbcsr_init(matrix_S_mo_sym,error=error)
!FROM UNMERGED     CALL cp_dbcsr_create(matrix_S_mo_sym,&
!FROM UNMERGED             template=matrix_S_mo,&
!FROM UNMERGED             matrix_type=dbcsr_type_no_symmetry,error=error)
!FROM UNMERGED     CALL cp_dbcsr_desymmetrize(matrix_S_mo,matrix_S_mo_sym,error=error)
!FROM UNMERGED 
!FROM UNMERGED     CALL cp_dbcsr_init(matrix_F_mo_sym,error=error)
!FROM UNMERGED     CALL cp_dbcsr_create(matrix_F_mo_sym,&
!FROM UNMERGED             template=matrix_F_mo,&
!FROM UNMERGED             matrix_type=dbcsr_type_no_symmetry,error=error)
!FROM UNMERGED     CALL cp_dbcsr_desymmetrize(matrix_F_mo,matrix_F_mo_sym,error=error)
!FROM UNMERGED    
!FROM UNMERGED     IF (penalty_occ_vol) THEN
!FROM UNMERGED        penalty_prefactor_local = penalty_occ_vol_prefactor/(2.0_dp*spin_factor)
!FROM UNMERGED     ELSE
!FROM UNMERGED        penalty_prefactor_local = 0.0_dp 
!FROM UNMERGED     ENDIF
!FROM UNMERGED 
!FROM UNMERGED WRITE(*,*) "penalty_prefactor_local: ", penalty_prefactor_local
!FROM UNMERGED WRITE(*,*) "penalty_prefactor_2: ", penalty_occ_vol_pf2
!FROM UNMERGED     
!FROM UNMERGED     !CALL cp_dbcsr_print(matrix_grad,error=error)
!FROM UNMERGED     !CALL cp_dbcsr_print(matrix_F_ao_sym,error=error)
!FROM UNMERGED     !CALL cp_dbcsr_print(matrix_S_ao_sym,error=error)
!FROM UNMERGED     !CALL cp_dbcsr_print(matrix_F_mo_sym,error=error)
!FROM UNMERGED     !CALL cp_dbcsr_print(matrix_S_mo_sym,error=error)
!FROM UNMERGED 
!FROM UNMERGED     ! loop over domains to find the size of the Hessian
!FROM UNMERGED     H_size=0
!FROM UNMERGED     DO col = 1, nblkcols_tot
!FROM UNMERGED 
!FROM UNMERGED        ! find sizes of AO submatrices
!FROM UNMERGED        DO row = 1, nblkrows_tot
!FROM UNMERGED           
!FROM UNMERGED           CALL cp_dbcsr_get_block_p(quench_t,&
!FROM UNMERGED                   row, col, block_p, found)
!FROM UNMERGED           IF (found) THEN
!FROM UNMERGED             ao_domain_sizes(col)=ao_domain_sizes(col)+ao_blk_sizes(row)
!FROM UNMERGED           ENDIF
!FROM UNMERGED 
!FROM UNMERGED        ENDDO
!FROM UNMERGED 
!FROM UNMERGED        H_size=H_size+ao_domain_sizes(col)*mo_block_sizes(col)
!FROM UNMERGED 
!FROM UNMERGED     ENDDO
!FROM UNMERGED 
!FROM UNMERGED     ALLOCATE(H(H_size,H_size))
!FROM UNMERGED     H(:,:)=0.0_dp
!FROM UNMERGED 
!FROM UNMERGED     ! fill the Hessian matrix
!FROM UNMERGED     lev1_vert_offset=0
!FROM UNMERGED     ! loop over all pairs of fragments
!FROM UNMERGED     DO row = 1, nblkcols_tot
!FROM UNMERGED 
!FROM UNMERGED        lev1_hori_offset=0
!FROM UNMERGED        DO col = 1, nblkcols_tot
!FROM UNMERGED 
!FROM UNMERGED           ! prepare blocks for the current row-column fragment pair
!FROM UNMERGED           ALLOCATE(F_ao_block(ao_domain_sizes(row),ao_domain_sizes(col)))
!FROM UNMERGED           ALLOCATE(S_ao_block(ao_domain_sizes(row),ao_domain_sizes(col)))
!FROM UNMERGED           ALLOCATE(F_mo_block(mo_block_sizes(row),mo_block_sizes(col)))
!FROM UNMERGED           ALLOCATE(S_mo_block(mo_block_sizes(row),mo_block_sizes(col)))
!FROM UNMERGED    
!FROM UNMERGED           F_ao_block(:,:)=0.0_dp
!FROM UNMERGED           S_ao_block(:,:)=0.0_dp
!FROM UNMERGED           F_mo_block(:,:)=0.0_dp
!FROM UNMERGED           S_mo_block(:,:)=0.0_dp
!FROM UNMERGED    
!FROM UNMERGED           ! fill AO submatrices
!FROM UNMERGED           ! loop over all blocks of the AO dbcsr matrix
!FROM UNMERGED           ao_vert_offset=0
!FROM UNMERGED           DO block_row = 1, nblkcols_tot
!FROM UNMERGED              
!FROM UNMERGED              CALL cp_dbcsr_get_block_p(quench_t,&
!FROM UNMERGED                      block_row, row, block_p, found_row)
!FROM UNMERGED              IF (found_row) THEN
!FROM UNMERGED       
!FROM UNMERGED                 ao_hori_offset=0
!FROM UNMERGED                 DO block_col = 1, nblkcols_tot
!FROM UNMERGED                    
!FROM UNMERGED                    CALL cp_dbcsr_get_block_p(quench_t,&
!FROM UNMERGED                            block_col, col, block_p, found_col)
!FROM UNMERGED                    IF (found_col) THEN
!FROM UNMERGED          
!FROM UNMERGED                       CALL cp_dbcsr_get_block_p(matrix_F_ao_sym,&
!FROM UNMERGED                               block_row, block_col, block_p, found)
!FROM UNMERGED                       IF (found) THEN
!FROM UNMERGED                          ! copy the block into the submatrix
!FROM UNMERGED                          F_ao_block(ao_vert_offset+1:ao_vert_offset+ao_block_sizes(block_row),&
!FROM UNMERGED                            ao_hori_offset+1:ao_hori_offset+ao_block_sizes(block_col))&
!FROM UNMERGED                            =block_p(:,:)
!FROM UNMERGED                       ENDIF
!FROM UNMERGED          
!FROM UNMERGED                       CALL cp_dbcsr_get_block_p(matrix_S_ao_sym,&
!FROM UNMERGED                               block_row, block_col, block_p, found)
!FROM UNMERGED                       IF (found) THEN
!FROM UNMERGED                          ! copy the block into the submatrix
!FROM UNMERGED                          S_ao_block(ao_vert_offset+1:ao_vert_offset+ao_block_sizes(block_row),&
!FROM UNMERGED                            ao_hori_offset+1:ao_hori_offset+ao_block_sizes(block_col))&
!FROM UNMERGED                            =block_p(:,:)
!FROM UNMERGED                       ENDIF
!FROM UNMERGED          
!FROM UNMERGED                       ao_hori_offset=ao_hori_offset+ao_block_sizes(block_col)
!FROM UNMERGED                    
!FROM UNMERGED                    ENDIF
!FROM UNMERGED       
!FROM UNMERGED                 ENDDO
!FROM UNMERGED 
!FROM UNMERGED                 ao_vert_offset=ao_vert_offset+ao_block_sizes(block_row)
!FROM UNMERGED 
!FROM UNMERGED              ENDIF
!FROM UNMERGED    
!FROM UNMERGED           ENDDO
!FROM UNMERGED 
!FROM UNMERGED           ! fill MO submatrices
!FROM UNMERGED           CALL cp_dbcsr_get_block_p(matrix_F_mo_sym, row, col, block_p, found)
!FROM UNMERGED           IF (found) THEN
!FROM UNMERGED              ! copy the block into the submatrix
!FROM UNMERGED              F_mo_block(1:mo_block_sizes(row),1:mo_block_sizes(col))=block_p(:,:)
!FROM UNMERGED           ENDIF
!FROM UNMERGED           CALL cp_dbcsr_get_block_p(matrix_S_mo_sym, row, col, block_p, found)
!FROM UNMERGED           IF (found) THEN
!FROM UNMERGED              ! copy the block into the submatrix
!FROM UNMERGED              S_mo_block(1:mo_block_sizes(row),1:mo_block_sizes(col))=block_p(:,:)
!FROM UNMERGED           ENDIF
!FROM UNMERGED 
!FROM UNMERGED           !WRITE(*,*) "F_AO_BLOCK", row, col, ao_domain_sizes(row), ao_domain_sizes(col)
!FROM UNMERGED           !DO ii=1,ao_domain_sizes(row)
!FROM UNMERGED           !  WRITE(*,'(100F13.9)') F_ao_block(ii,:)
!FROM UNMERGED           !ENDDO
!FROM UNMERGED           !WRITE(*,*) "S_AO_BLOCK", row, col
!FROM UNMERGED           !DO ii=1,ao_domain_sizes(row)
!FROM UNMERGED           !  WRITE(*,'(100F13.9)') S_ao_block(ii,:)
!FROM UNMERGED           !ENDDO
!FROM UNMERGED           !WRITE(*,*) "F_MO_BLOCK", row, col
!FROM UNMERGED           !DO ii=1,mo_block_sizes(row)
!FROM UNMERGED           !  WRITE(*,'(100F13.9)') F_mo_block(ii,:)
!FROM UNMERGED           !ENDDO
!FROM UNMERGED           !WRITE(*,*) "S_MO_BLOCK", row, col, mo_block_sizes(row), mo_block_sizes(col)
!FROM UNMERGED           !DO ii=1,mo_block_sizes(row)
!FROM UNMERGED           !  WRITE(*,'(100F13.9)') S_mo_block(ii,:)
!FROM UNMERGED           !ENDDO
!FROM UNMERGED 
!FROM UNMERGED           ! construct tensor products for the current row-column fragment pair
!FROM UNMERGED           lev2_vert_offset=0
!FROM UNMERGED           DO orb_j=1,mo_block_sizes(row)
!FROM UNMERGED 
!FROM UNMERGED              lev2_hori_offset=0
!FROM UNMERGED              DO orb_i=1,mo_block_sizes(col)
!FROM UNMERGED IF (orb_i.EQ.orb_j.AND.row.EQ.col) THEN
!FROM UNMERGED                 H(lev1_vert_offset+lev2_vert_offset+1:lev1_vert_offset+lev2_vert_offset+ao_domain_sizes(row),&
!FROM UNMERGED                   lev1_hori_offset+lev2_hori_offset+1:lev1_hori_offset+lev2_hori_offset+ao_domain_sizes(col))&
!FROM UNMERGED != -penalty_prefactor_local*S_ao_block(:,:)
!FROM UNMERGED =F_ao_block(:,:)+S_ao_block(:,:)
!FROM UNMERGED !=S_ao_block(:,:)
!FROM UNMERGED !RZK-warning               =F_ao_block(:,:)+( 1.0_dp + penalty_prefactor_local )*S_ao_block(:,:)
!FROM UNMERGED !               =S_mo_block(orb_j,orb_i)*F_ao_block(:,:)&
!FROM UNMERGED !               -F_mo_block(orb_j,orb_i)*S_ao_block(:,:)&
!FROM UNMERGED !               +penalty_prefactor_local*S_mo_block(orb_j,orb_i)*S_ao_block(:,:)
!FROM UNMERGED ENDIF
!FROM UNMERGED                 !WRITE(*,*) row, col, orb_j, orb_i, lev1_vert_offset+lev2_vert_offset+1, ao_domain_sizes(row),&
!FROM UNMERGED                 !   lev1_hori_offset+lev2_hori_offset+1, ao_domain_sizes(col), S_mo_block(orb_j,orb_i)
!FROM UNMERGED 
!FROM UNMERGED                 lev2_hori_offset=lev2_hori_offset+ao_domain_sizes(col)
!FROM UNMERGED 
!FROM UNMERGED              ENDDO
!FROM UNMERGED 
!FROM UNMERGED              lev2_vert_offset=lev2_vert_offset+ao_domain_sizes(row)
!FROM UNMERGED 
!FROM UNMERGED           ENDDO
!FROM UNMERGED 
!FROM UNMERGED           
!FROM UNMERGED           lev1_hori_offset=lev1_hori_offset+ao_domain_sizes(col)*mo_block_sizes(col)
!FROM UNMERGED 
!FROM UNMERGED           DEALLOCATE(F_ao_block)
!FROM UNMERGED           DEALLOCATE(S_ao_block)
!FROM UNMERGED           DEALLOCATE(F_mo_block)
!FROM UNMERGED           DEALLOCATE(S_mo_block)
!FROM UNMERGED 
!FROM UNMERGED        ENDDO ! col fragment
!FROM UNMERGED 
!FROM UNMERGED        lev1_vert_offset=lev1_vert_offset+ao_domain_sizes(row)*mo_block_sizes(row)
!FROM UNMERGED 
!FROM UNMERGED     ENDDO ! row fragment
!FROM UNMERGED 
!FROM UNMERGED     CALL cp_dbcsr_release(matrix_S_ao_sym,error=error)
!FROM UNMERGED     CALL cp_dbcsr_release(matrix_F_ao_sym,error=error)
!FROM UNMERGED     CALL cp_dbcsr_release(matrix_S_mo_sym,error=error)
!FROM UNMERGED     CALL cp_dbcsr_release(matrix_F_mo_sym,error=error)
!FROM UNMERGED 
!FROM UNMERGED !!    ! Two more terms of the Hessian: S_vo.D.F_vo and F_vo.D.S_vo
!FROM UNMERGED !!    ! It seems that these terms break positive definite property of the Hessian
!FROM UNMERGED !!    ALLOCATE(H1(H_size,H_size))
!FROM UNMERGED !!    ALLOCATE(H2(H_size,H_size))
!FROM UNMERGED !!    H1=0.0_dp
!FROM UNMERGED !!    H2=0.0_dp
!FROM UNMERGED !!    DO row = 1, nblkcols_tot
!FROM UNMERGED !!
!FROM UNMERGED !!       lev1_hori_offset=0
!FROM UNMERGED !!       DO col = 1, nblkcols_tot
!FROM UNMERGED !!
!FROM UNMERGED !!          CALL cp_dbcsr_get_block_p(matrix_F_vo,&
!FROM UNMERGED !!                  row, col, block_p, found)
!FROM UNMERGED !!          CALL cp_dbcsr_get_block_p(matrix_S_vo,&
!FROM UNMERGED !!                  row, col, block_p2, found2)
!FROM UNMERGED !!
!FROM UNMERGED !!          lev1_vert_offset=0
!FROM UNMERGED !!          DO block_col = 1, nblkcols_tot
!FROM UNMERGED !!             
!FROM UNMERGED !!             CALL cp_dbcsr_get_block_p(quench_t,&
!FROM UNMERGED !!                     row, block_col, p_new_block, found_row)
!FROM UNMERGED !!
!FROM UNMERGED !!             IF (found_row) THEN
!FROM UNMERGED !!
!FROM UNMERGED !!                ! determine offset in this short loop
!FROM UNMERGED !!                lev2_vert_offset=0
!FROM UNMERGED !!                DO block_row=1,row-1
!FROM UNMERGED !!                   CALL cp_dbcsr_get_block_p(quench_t,&
!FROM UNMERGED !!                           block_row, block_col, p_new_block, found_col)
!FROM UNMERGED !!                   IF (found_col) lev2_vert_offset=lev2_vert_offset+ao_block_sizes(block_row)
!FROM UNMERGED !!                ENDDO
!FROM UNMERGED !!                !!!!!!!! short loop
!FROM UNMERGED !!
!FROM UNMERGED !!                ! over all electrons of the block
!FROM UNMERGED !!                DO orb_i=1, mo_block_sizes(col)
!FROM UNMERGED !!                
!FROM UNMERGED !!                   ! into all possible locations
!FROM UNMERGED !!                   DO orb_j=1, mo_block_sizes(block_col)
!FROM UNMERGED !!   
!FROM UNMERGED !!                      ! column is copied several times
!FROM UNMERGED !!                      DO copy=1, ao_domain_sizes(col)
!FROM UNMERGED !!   
!FROM UNMERGED !!                         IF (found) THEN
!FROM UNMERGED !!
!FROM UNMERGED !!                            !WRITE(*,*) row, col, block_col, orb_i, orb_j, copy,&
!FROM UNMERGED !!                            ! lev1_vert_offset+(orb_j-1)*ao_domain_sizes(block_col)+lev2_vert_offset+1,&
!FROM UNMERGED !!                            ! lev1_hori_offset+(orb_i-1)*ao_domain_sizes(col)+copy
!FROM UNMERGED !!
!FROM UNMERGED !!                            H1( lev1_vert_offset+(orb_j-1)*ao_domain_sizes(block_col)+lev2_vert_offset+1:&
!FROM UNMERGED !!                                lev1_vert_offset+(orb_j-1)*ao_domain_sizes(block_col)+lev2_vert_offset+ao_block_sizes(row),&
!FROM UNMERGED !!                                lev1_hori_offset+(orb_i-1)*ao_domain_sizes(col)+copy )&
!FROM UNMERGED !!                              =block_p(:,orb_i)
!FROM UNMERGED !!   
!FROM UNMERGED !!                         ENDIF ! found block in the data matrix
!FROM UNMERGED !!       
!FROM UNMERGED !!                         IF (found2) THEN
!FROM UNMERGED !!          
!FROM UNMERGED !!                            H2( lev1_vert_offset+(orb_j-1)*ao_domain_sizes(block_col)+lev2_vert_offset+1:&
!FROM UNMERGED !!                                lev1_vert_offset+(orb_j-1)*ao_domain_sizes(block_col)+lev2_vert_offset+ao_block_sizes(row),&
!FROM UNMERGED !!                                lev1_hori_offset+(orb_i-1)*ao_domain_sizes(col)+copy )&
!FROM UNMERGED !!                              =block_p2(:,orb_i)
!FROM UNMERGED !!   
!FROM UNMERGED !!                         ENDIF ! found block in the data matrix
!FROM UNMERGED !!
!FROM UNMERGED !!                      ENDDO
!FROM UNMERGED !!   
!FROM UNMERGED !!                   ENDDO
!FROM UNMERGED !!
!FROM UNMERGED !!                ENDDO
!FROM UNMERGED !!                
!FROM UNMERGED !!                !lev2_vert_offset=lev2_vert_offset+ao_block_sizes(row)
!FROM UNMERGED !!
!FROM UNMERGED !!             ENDIF ! found block in the quench matrix
!FROM UNMERGED !!
!FROM UNMERGED !!             lev1_vert_offset=lev1_vert_offset+&
!FROM UNMERGED !!                ao_domain_sizes(block_col)*mo_block_sizes(block_col)
!FROM UNMERGED !!
!FROM UNMERGED !!          ENDDO
!FROM UNMERGED !!
!FROM UNMERGED !!          lev1_hori_offset=lev1_hori_offset+&
!FROM UNMERGED !!             ao_domain_sizes(col)*mo_block_sizes(col)
!FROM UNMERGED !!
!FROM UNMERGED !!       ENDDO
!FROM UNMERGED !!
!FROM UNMERGED !!       !lev2_vert_offset=lev2_vert_offset+ao_block_sizes(row)
!FROM UNMERGED !!
!FROM UNMERGED !!    ENDDO
!FROM UNMERGED !!    H1(:,:)=H1(:,:)*2.0_dp*spin_factor
!FROM UNMERGED !!    !!!WRITE(*,*) "F_vo"
!FROM UNMERGED !!    !!!DO ii=1,H_size
!FROM UNMERGED !!    !!! WRITE(*,'(100F13.9)') H1(ii,:)
!FROM UNMERGED !!    !!!ENDDO
!FROM UNMERGED !!    !!!WRITE(*,*) "S_vo"
!FROM UNMERGED !!    !!!DO ii=1,H_size
!FROM UNMERGED !!    !!! WRITE(*,'(100F13.9)') H2(ii,:)
!FROM UNMERGED !!    !!!ENDDO
!FROM UNMERGED !!    !!!!! add terms to the hessian
!FROM UNMERGED !!    DO ii=1,H_size
!FROM UNMERGED !!       DO jj=1,H_size
!FROM UNMERGED !!! add penalty_occ_vol term
!FROM UNMERGED !!          H(ii,jj)=H(ii,jj)-H1(ii,jj)*H2(jj,ii)-H1(jj,ii)*H2(ii,jj)
!FROM UNMERGED !!       ENDDO
!FROM UNMERGED !!    ENDDO
!FROM UNMERGED !!    DEALLOCATE(H1)
!FROM UNMERGED !!    DEALLOCATE(H2)
!FROM UNMERGED              
!FROM UNMERGED !!    ! S_vo.S_vo diagonal component due to determiant constraint
!FROM UNMERGED !!    ! use grad vector temporarily
!FROM UNMERGED !!    IF (penalty_occ_vol) THEN
!FROM UNMERGED !!       ALLOCATE(Grad_vec(H_size))
!FROM UNMERGED !!       Grad_vec(:)=0.0_dp
!FROM UNMERGED !!       lev1_vert_offset=0
!FROM UNMERGED !!       ! loop over all electron blocks 
!FROM UNMERGED !!       DO col = 1, nblkcols_tot
!FROM UNMERGED !!
!FROM UNMERGED !!          ! loop over AO-rows of the dbcsr matrix
!FROM UNMERGED !!          lev2_vert_offset=0
!FROM UNMERGED !!          DO row = 1, nblkrows_tot
!FROM UNMERGED !!             
!FROM UNMERGED !!             CALL cp_dbcsr_get_block_p(quench_t,&
!FROM UNMERGED !!                     row, col, block_p, found_row)
!FROM UNMERGED !!             IF (found_row) THEN
!FROM UNMERGED !!         
!FROM UNMERGED !!                CALL cp_dbcsr_get_block_p(matrix_S_vo,&
!FROM UNMERGED !!                        row, col, block_p, found)
!FROM UNMERGED !!                IF (found) THEN
!FROM UNMERGED !!                   ! copy the data into the vector, column by column
!FROM UNMERGED !!                   DO orb_i=1, mo_block_sizes(col)
!FROM UNMERGED !!                      Grad_vec(lev1_vert_offset+ao_domain_sizes(col)*(orb_i-1)+lev2_vert_offset+1:&
!FROM UNMERGED !!                               lev1_vert_offset+ao_domain_sizes(col)*(orb_i-1)+lev2_vert_offset+ao_block_sizes(row))&
!FROM UNMERGED !!                               =block_p(:,orb_i)
!FROM UNMERGED !!                   ENDDO
!FROM UNMERGED !!
!FROM UNMERGED !!                ENDIF
!FROM UNMERGED !!          
!FROM UNMERGED !!                lev2_vert_offset=lev2_vert_offset+ao_block_sizes(row)
!FROM UNMERGED !!
!FROM UNMERGED !!             ENDIF
!FROM UNMERGED !!   
!FROM UNMERGED !!          ENDDO
!FROM UNMERGED !!
!FROM UNMERGED !!          lev1_vert_offset=lev1_vert_offset+ao_domain_sizes(col)*mo_block_sizes(col)
!FROM UNMERGED !!
!FROM UNMERGED !!       ENDDO ! loop over electron blocks 
!FROM UNMERGED !!       ! update H now
!FROM UNMERGED !!       DO ii=1,H_size
!FROM UNMERGED !!          DO jj=1,H_size
!FROM UNMERGED !!             H(ii,jj)=H(ii,jj)+penalty_occ_vol_prefactor*&
!FROM UNMERGED !!                      penalty_occ_vol_pf2*Grad_vec(ii)*Grad_vec(jj)
!FROM UNMERGED !!          ENDDO
!FROM UNMERGED !!       ENDDO
!FROM UNMERGED !!       DEALLOCATE(Grad_vec)
!FROM UNMERGED !!    ENDIF ! penalty_occ_vol
!FROM UNMERGED 
!FROM UNMERGED !S-1.G ! invert S using cholesky 
!FROM UNMERGED !S-1.G CALL cp_dbcsr_init(m_prec_out,error=error)
!FROM UNMERGED !S-1.G CALL cp_dbcsr_create(m_prec_out,&
!FROM UNMERGED !S-1.G         template=m_s,&
!FROM UNMERGED !S-1.G         matrix_type=dbcsr_type_no_symmetry,&
!FROM UNMERGED !S-1.G         error=error)
!FROM UNMERGED !S-1.G CALL cp_dbcsr_copy(m_prec_out,m_s,error=error)
!FROM UNMERGED !S-1.G CALL cp_dbcsr_cholesky_decompose(m_prec_out,&
!FROM UNMERGED !S-1.G         para_env=para_env,&
!FROM UNMERGED !S-1.G         blacs_env=blacs_env,error=error)
!FROM UNMERGED !S-1.G CALL cp_dbcsr_cholesky_invert(m_prec_out,&
!FROM UNMERGED !S-1.G         para_env=para_env,&
!FROM UNMERGED !S-1.G         blacs_env=blacs_env,&
!FROM UNMERGED !S-1.G         upper_to_full=.TRUE.,error=error)
!FROM UNMERGED !S-1.G CALL cp_dbcsr_multiply("N","N",1.0_dp,&
!FROM UNMERGED !S-1.G         m_prec_out,&
!FROM UNMERGED !S-1.G         matrix_grad,&
!FROM UNMERGED !S-1.G         0.0_dp,matrix_step,&
!FROM UNMERGED !S-1.G         filter_eps=1.0E-10_dp,&
!FROM UNMERGED !S-1.G         error=error)
!FROM UNMERGED !S-1.G !CALL cp_dbcsr_release(m_prec_out,error=error)
!FROM UNMERGED !S-1.G ALLOCATE(test3(H_size))
!FROM UNMERGED     
!FROM UNMERGED     ! convert gradient from the dbcsr matrix to the vector form
!FROM UNMERGED     ALLOCATE(Grad_vec(H_size))
!FROM UNMERGED     Grad_vec(:)=0.0_dp
!FROM UNMERGED     lev1_vert_offset=0
!FROM UNMERGED     ! loop over all electron blocks 
!FROM UNMERGED     DO col = 1, nblkcols_tot
!FROM UNMERGED 
!FROM UNMERGED        ! loop over AO-rows of the dbcsr matrix
!FROM UNMERGED        lev2_vert_offset=0
!FROM UNMERGED        DO row = 1, nblkrows_tot
!FROM UNMERGED           
!FROM UNMERGED           CALL cp_dbcsr_get_block_p(quench_t,&
!FROM UNMERGED                   row, col, block_p, found_row)
!FROM UNMERGED           IF (found_row) THEN
!FROM UNMERGED       
!FROM UNMERGED              CALL cp_dbcsr_get_block_p(matrix_grad,&
!FROM UNMERGED                      row, col, block_p, found)
!FROM UNMERGED              IF (found) THEN
!FROM UNMERGED                 ! copy the data into the vector, column by column
!FROM UNMERGED                 DO orb_i=1, mo_block_sizes(col)
!FROM UNMERGED                    Grad_vec(lev1_vert_offset+ao_domain_sizes(col)*(orb_i-1)+lev2_vert_offset+1:&
!FROM UNMERGED                             lev1_vert_offset+ao_domain_sizes(col)*(orb_i-1)+lev2_vert_offset+ao_block_sizes(row))&
!FROM UNMERGED                             =block_p(:,orb_i)
!FROM UNMERGED !WRITE(*,*) "GRAD: ", row, col, orb_i, lev1_vert_offset+ao_domain_sizes(col)*(orb_i-1)+lev2_vert_offset+1, ao_block_sizes(row)
!FROM UNMERGED                 ENDDO
!FROM UNMERGED 
!FROM UNMERGED              ENDIF
!FROM UNMERGED        
!FROM UNMERGED !S-1.G CALL cp_dbcsr_get_block_p(matrix_step,&
!FROM UNMERGED !S-1.G         row, col, block_p, found)
!FROM UNMERGED !S-1.G IF (found) THEN
!FROM UNMERGED !S-1.G    ! copy the data into the vector, column by column
!FROM UNMERGED !S-1.G    DO orb_i=1, mo_block_sizes(col)
!FROM UNMERGED !S-1.G       test3(lev1_vert_offset+ao_domain_sizes(col)*(orb_i-1)+lev2_vert_offset+1:&
!FROM UNMERGED !S-1.G                lev1_vert_offset+ao_domain_sizes(col)*(orb_i-1)+lev2_vert_offset+ao_block_sizes(row))&
!FROM UNMERGED !S-1.G                =block_p(:,orb_i)
!FROM UNMERGED !S-1.G    ENDDO
!FROM UNMERGED !S-1.G ENDIF
!FROM UNMERGED              
!FROM UNMERGED              lev2_vert_offset=lev2_vert_offset+ao_block_sizes(row)
!FROM UNMERGED 
!FROM UNMERGED           ENDIF
!FROM UNMERGED    
!FROM UNMERGED        ENDDO
!FROM UNMERGED 
!FROM UNMERGED        lev1_vert_offset=lev1_vert_offset+ao_domain_sizes(col)*mo_block_sizes(col)
!FROM UNMERGED 
!FROM UNMERGED     ENDDO ! loop over electron blocks 
!FROM UNMERGED 
!FROM UNMERGED     !WRITE(*,*) "HESSIAN"
!FROM UNMERGED     !DO ii=1,H_size
!FROM UNMERGED     ! WRITE(*,*) ii
!FROM UNMERGED     ! WRITE(*,'(20F14.10)') H(ii,:)
!FROM UNMERGED     !ENDDO
!FROM UNMERGED 
!FROM UNMERGED     ! invert the Hessian
!FROM UNMERGED     INFO=0
!FROM UNMERGED     ALLOCATE(Hinv(H_size,H_size))
!FROM UNMERGED     Hinv(:,:)=H(:,:)
!FROM UNMERGED 
!FROM UNMERGED     ! before inverting diagonalize
!FROM UNMERGED     ALLOCATE(eigenvalues(H_size))
!FROM UNMERGED     ! Query the optimal workspace for dsyev
!FROM UNMERGED     LWORK = -1
!FROM UNMERGED     ALLOCATE(WORK(MAX(1,LWORK)))
!FROM UNMERGED     CALL DSYEV('V','L',H_size,Hinv,H_size,eigenvalues,WORK,LWORK,INFO)
!FROM UNMERGED     LWORK = INT(WORK( 1 ))
!FROM UNMERGED     DEALLOCATE(WORK)
!FROM UNMERGED     ! Allocate the workspace and solve the eigenproblem
!FROM UNMERGED     ALLOCATE(WORK(MAX(1,LWORK)))
!FROM UNMERGED     CALL DSYEV('V','L',H_size,Hinv,H_size,eigenvalues,WORK,LWORK,INFO)
!FROM UNMERGED     IF( INFO.NE.0 ) THEN
!FROM UNMERGED        WRITE(*,*) 'DSYEV ERROR MESSAGE: ', INFO
!FROM UNMERGED        CPErrorMessage(cp_failure_level,routineP,"DSYEV failed",error)
!FROM UNMERGED        CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
!FROM UNMERGED     END IF
!FROM UNMERGED     DEALLOCATE(WORK)
!FROM UNMERGED 
!FROM UNMERGED     ! compute grad vector in the basis of Hessian eigenvectors
!FROM UNMERGED     ALLOCATE(Step_vec(H_size))
!FROM UNMERGED     ! Step_vec contains Grad_vec here
!FROM UNMERGED     Step_vec(:)=MATMUL(TRANSPOSE(Hinv),Grad_vec)
!FROM UNMERGED 
!FROM UNMERGED     ! compute U.tr(U)-1 = error
!FROM UNMERGED     !ALLOCATE(test(H_size,H_size))
!FROM UNMERGED     !test(:,:)=MATMUL(TRANSPOSE(Hinv),Hinv)
!FROM UNMERGED     !DO ii=1,H_size
!FROM UNMERGED     !   test(ii,ii)=test(ii,ii)-1.0_dp
!FROM UNMERGED     !ENDDO
!FROM UNMERGED     !test_error=0.0_dp
!FROM UNMERGED     !DO ii=1,H_size
!FROM UNMERGED     !   DO jj=1,H_size
!FROM UNMERGED     !      test_error=test_error+test(jj,ii)*test(jj,ii)
!FROM UNMERGED     !   ENDDO
!FROM UNMERGED     !ENDDO
!FROM UNMERGED     !WRITE(*,*) "U.tr(U)-1 error: ", SQRT(test_error)
!FROM UNMERGED     !DEALLOCATE(test)
!FROM UNMERGED     
!FROM UNMERGED     ! invert eigenvalues and use eigenvectors to compute the Hessian inverse
!FROM UNMERGED     ! project out zero-eigenvalue directions
!FROM UNMERGED     ALLOCATE(test(H_size,H_size))
!FROM UNMERGED     zero_neg_eiv=0
!FROM UNMERGED     DO jj=1, H_size
!FROM UNMERGED WRITE(*,"(I10,F20.10,F20.10)") jj, eigenvalues(jj), Step_vec(jj)
!FROM UNMERGED        IF (eigenvalues(jj).gt.eps_zero) THEN
!FROM UNMERGED           test(jj,:)=Hinv(:,jj)/eigenvalues(jj)
!FROM UNMERGED        ELSE
!FROM UNMERGED           test(jj,:)=Hinv(:,jj)*0.0_dp
!FROM UNMERGED           zero_neg_eiv=zero_neg_eiv+1
!FROM UNMERGED        ENDIF
!FROM UNMERGED     ENDDO
!FROM UNMERGED     WRITE(*,*) 'ZERO OR NEGATIVE EIGENVALUES: ', zero_neg_eiv
!FROM UNMERGED     DEALLOCATE(Step_vec)
!FROM UNMERGED     
!FROM UNMERGED     ALLOCATE(test2(H_size,H_size))
!FROM UNMERGED     test2(:,:)=MATMUL(Hinv,test)
!FROM UNMERGED     Hinv(:,:)=test2(:,:)
!FROM UNMERGED     DEALLOCATE(test,test2)
!FROM UNMERGED 
!FROM UNMERGED     !! shift to kill singularity
!FROM UNMERGED     !shift=0.0_dp
!FROM UNMERGED     !IF (eigenvalues(1).lt.0.0_dp) THEN
!FROM UNMERGED     !   CPErrorMessage(cp_failure_level,routineP,"Negative eigenvalue(s)",error)
!FROM UNMERGED     !   shift=abs(eigenvalues(1))
!FROM UNMERGED     !   WRITE(*,*) "Lowest eigenvalue: ", eigenvalues(1)
!FROM UNMERGED     !ENDIF
!FROM UNMERGED     !DO ii=1, H_size
!FROM UNMERGED     !   IF (eigenvalues(ii).gt.eps_zero) THEN
!FROM UNMERGED     !      shift=shift+min(1.0_dp,eigenvalues(ii))*1.0E-4_dp
!FROM UNMERGED     !      EXIT
!FROM UNMERGED     !   ENDIF
!FROM UNMERGED     !ENDDO
!FROM UNMERGED     !WRITE(*,*) "Hessian shift: ", shift
!FROM UNMERGED     !DO ii=1, H_size
!FROM UNMERGED     !   H(ii,ii)=H(ii,ii)+shift
!FROM UNMERGED     !ENDDO
!FROM UNMERGED     !! end shift
!FROM UNMERGED     
!FROM UNMERGED     DEALLOCATE(eigenvalues)
!FROM UNMERGED     
!FROM UNMERGED !!!!    Hinv=H
!FROM UNMERGED !!!!    INFO=0
!FROM UNMERGED !!!!    CALL DPOTRF('L', H_size, Hinv, H_size, INFO )
!FROM UNMERGED !!!!    IF( INFO.NE.0 ) THEN
!FROM UNMERGED !!!!       WRITE(*,*) 'DPOTRF ERROR MESSAGE: ', INFO
!FROM UNMERGED !!!!       CPErrorMessage(cp_failure_level,routineP,"DPOTRF failed",error)
!FROM UNMERGED !!!!       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
!FROM UNMERGED !!!!    END IF
!FROM UNMERGED !!!!    CALL DPOTRI('L', H_size, Hinv, H_size, INFO )
!FROM UNMERGED !!!!    IF( INFO.NE.0 ) THEN
!FROM UNMERGED !!!!       WRITE(*,*) 'DPOTRI ERROR MESSAGE: ', INFO
!FROM UNMERGED !!!!       CPErrorMessage(cp_failure_level,routineP,"DPOTRI failed",error)
!FROM UNMERGED !!!!       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
!FROM UNMERGED !!!!    END IF
!FROM UNMERGED !!!!    ! complete the matrix
!FROM UNMERGED !!!!    DO ii=1,H_size
!FROM UNMERGED !!!!       DO jj=ii+1,H_size
!FROM UNMERGED !!!!          Hinv(ii,jj)=Hinv(jj,ii)
!FROM UNMERGED !!!!       ENDDO
!FROM UNMERGED !!!!    ENDDO
!FROM UNMERGED 
!FROM UNMERGED     ! compute the inversion error 
!FROM UNMERGED     ALLOCATE(test(H_size,H_size))
!FROM UNMERGED     test(:,:)=MATMUL(Hinv,H)
!FROM UNMERGED     DO ii=1,H_size
!FROM UNMERGED        test(ii,ii)=test(ii,ii)-1.0_dp
!FROM UNMERGED     ENDDO
!FROM UNMERGED     test_error=0.0_dp
!FROM UNMERGED     DO ii=1,H_size
!FROM UNMERGED        DO jj=1,H_size
!FROM UNMERGED           test_error=test_error+test(jj,ii)*test(jj,ii)
!FROM UNMERGED        ENDDO
!FROM UNMERGED     ENDDO
!FROM UNMERGED     WRITE(*,*) "Hessian inversion error: ", SQRT(test_error)
!FROM UNMERGED     DEALLOCATE(test)
!FROM UNMERGED 
!FROM UNMERGED     ! prepare the output vector
!FROM UNMERGED     ALLOCATE(Step_vec(H_size))
!FROM UNMERGED     ALLOCATE(tmp(H_size))
!FROM UNMERGED     tmp(:)=MATMUL(Hinv,Grad_vec)
!FROM UNMERGED     !tmp(:)=MATMUL(Hinv,test3)
!FROM UNMERGED     Step_vec(:)=-1.0_dp*tmp(:)
!FROM UNMERGED     
!FROM UNMERGED     ALLOCATE(tmpr(H_size))
!FROM UNMERGED     tmpr(:)=MATMUL(H,Step_vec)
!FROM UNMERGED     tmp(:)=tmpr(:)+Grad_vec(:)
!FROM UNMERGED     DEALLOCATE(tmpr)
!FROM UNMERGED     WRITE(*,*) "NEWTOV step error: ", MAXVAL(ABS(tmp))
!FROM UNMERGED 
!FROM UNMERGED     DEALLOCATE(tmp)
!FROM UNMERGED 
!FROM UNMERGED     DEALLOCATE(H)
!FROM UNMERGED     DEALLOCATE(Hinv)
!FROM UNMERGED     DEALLOCATE(Grad_vec)
!FROM UNMERGED 
!FROM UNMERGED !S-1.G DEALLOCATE(test3)
!FROM UNMERGED 
!FROM UNMERGED     ! copy the step from the vector into the cp_dbcsr matrix
!FROM UNMERGED     
!FROM UNMERGED     ! re-create the step matrix to remove all blocks
!FROM UNMERGED     CALL cp_dbcsr_create(matrix_step,&
!FROM UNMERGED             template=matrix_grad,&
!FROM UNMERGED             matrix_type=dbcsr_type_no_symmetry,error=error)
!FROM UNMERGED     CALL cp_dbcsr_work_create(matrix_step,work_mutable=.TRUE.,&
!FROM UNMERGED             error=error)
!FROM UNMERGED     
!FROM UNMERGED     lev1_vert_offset=0
!FROM UNMERGED     ! loop over all electron blocks 
!FROM UNMERGED     DO col = 1, nblkcols_tot
!FROM UNMERGED 
!FROM UNMERGED        ! loop over AO-rows of the dbcsr matrix
!FROM UNMERGED        lev2_vert_offset=0
!FROM UNMERGED        DO row = 1, nblkrows_tot
!FROM UNMERGED           
!FROM UNMERGED           CALL cp_dbcsr_get_block_p(quench_t,&
!FROM UNMERGED                   row, col, block_p, found_row)
!FROM UNMERGED           IF (found_row) THEN
!FROM UNMERGED       
!FROM UNMERGED              NULLIFY (p_new_block)
!FROM UNMERGED              CALL cp_dbcsr_reserve_block2d(matrix_step,row,col,p_new_block)
!FROM UNMERGED              CPPostcondition(ASSOCIATED(p_new_block),cp_failure_level,routineP,error,failure)
!FROM UNMERGED              ! copy the data column by column
!FROM UNMERGED              DO orb_i=1, mo_block_sizes(col)
!FROM UNMERGED                 p_new_block(:,orb_i) = &
!FROM UNMERGED                    Step_vec(lev1_vert_offset+ao_domain_sizes(col)*(orb_i-1)+lev2_vert_offset+1:&
!FROM UNMERGED                             lev1_vert_offset+ao_domain_sizes(col)*(orb_i-1)+lev2_vert_offset+ao_block_sizes(row))
!FROM UNMERGED              ENDDO
!FROM UNMERGED 
!FROM UNMERGED              lev2_vert_offset=lev2_vert_offset+ao_block_sizes(row)
!FROM UNMERGED 
!FROM UNMERGED           ENDIF
!FROM UNMERGED    
!FROM UNMERGED        ENDDO
!FROM UNMERGED 
!FROM UNMERGED        lev1_vert_offset=lev1_vert_offset+ao_domain_sizes(col)*mo_block_sizes(col)
!FROM UNMERGED 
!FROM UNMERGED     ENDDO ! loop over electron blocks 
!FROM UNMERGED 
!FROM UNMERGED     DEALLOCATE(Step_vec)
!FROM UNMERGED     
!FROM UNMERGED     CALL cp_dbcsr_finalize(matrix_step,error=error)
!FROM UNMERGED     !CALL cp_dbcsr_print(matrix_step,error=error)
!FROM UNMERGED 
!FROM UNMERGED !S-1.G CALL cp_dbcsr_init(m_tmp_no_1,error=error)
!FROM UNMERGED !S-1.G CALL cp_dbcsr_create(m_tmp_no_1,&
!FROM UNMERGED !S-1.G         template=matrix_step,&
!FROM UNMERGED !S-1.G         matrix_type=dbcsr_type_no_symmetry,error=error)
!FROM UNMERGED !S-1.G CALL cp_dbcsr_multiply("N","N",1.0_dp,&
!FROM UNMERGED !S-1.G         m_prec_out,&
!FROM UNMERGED !S-1.G         matrix_step,&
!FROM UNMERGED !S-1.G         0.0_dp,m_tmp_no_1,&
!FROM UNMERGED !S-1.G         filter_eps=1.0E-10_dp,&
!FROM UNMERGED !S-1.G         error=error)
!FROM UNMERGED !S-1.G CALL cp_dbcsr_copy(matrix_step,m_tmp_no_1,error=error)
!FROM UNMERGED !S-1.G CALL cp_dbcsr_release(m_tmp_no_1,error=error)
!FROM UNMERGED !S-1.G CALL cp_dbcsr_release(m_prec_out,error=error)
!FROM UNMERGED 
!FROM UNMERGED     DEALLOCATE(mo_block_sizes,ao_block_sizes)
!FROM UNMERGED     DEALLOCATE(ao_domain_sizes)
!FROM UNMERGED     
!FROM UNMERGED     CALL cp_dbcsr_init(matrix_S_ao_sym,error=error)
!FROM UNMERGED     CALL cp_dbcsr_create(matrix_S_ao_sym,&
!FROM UNMERGED             template=quench_t,&
!FROM UNMERGED             matrix_type=dbcsr_type_no_symmetry,error=error)
!FROM UNMERGED     CALL cp_dbcsr_copy(matrix_S_ao_sym,quench_t,error=error)
!FROM UNMERGED     CALL cp_dbcsr_multiply("N","N",1.0_dp,&
!FROM UNMERGED             matrix_F_ao,&
!FROM UNMERGED             matrix_step,&
!FROM UNMERGED             0.0_dp,matrix_S_ao_sym,&
!FROM UNMERGED             retain_sparsity=.TRUE.,&
!FROM UNMERGED             error=error)
!FROM UNMERGED     CALL cp_dbcsr_init(matrix_F_ao_sym,error=error)
!FROM UNMERGED     CALL cp_dbcsr_create(matrix_F_ao_sym,&
!FROM UNMERGED             template=quench_t,&
!FROM UNMERGED             matrix_type=dbcsr_type_no_symmetry,error=error)
!FROM UNMERGED     CALL cp_dbcsr_copy(matrix_F_ao_sym,quench_t,error=error)
!FROM UNMERGED     CALL cp_dbcsr_multiply("N","N",1.0_dp,&
!FROM UNMERGED             matrix_S_ao,&
!FROM UNMERGED             matrix_step,&
!FROM UNMERGED             0.0_dp,matrix_F_ao_sym,&
!FROM UNMERGED             retain_sparsity=.TRUE.,&
!FROM UNMERGED             error=error)
!FROM UNMERGED     CALL cp_dbcsr_add(matrix_S_ao_sym,matrix_F_ao_sym,&
!FROM UNMERGED          1.0_dp,1.0_dp,error=error)
!FROM UNMERGED     CALL cp_dbcsr_scale(matrix_S_ao_sym,2.0_dp*spin_factor,error=error)
!FROM UNMERGED     CALL cp_dbcsr_add(matrix_S_ao_sym,matrix_grad,&
!FROM UNMERGED          1.0_dp,1.0_dp,error=error)
!FROM UNMERGED     CALL cp_dbcsr_norm(matrix_S_ao_sym, dbcsr_norm_maxabsnorm,&
!FROM UNMERGED          norm_scalar=test_error, error=error)
!FROM UNMERGED     WRITE(*,*) "NEWTOL step error: ", test_error
!FROM UNMERGED     CALL cp_dbcsr_release(matrix_S_ao_sym,error=error)
!FROM UNMERGED     CALL cp_dbcsr_release(matrix_F_ao_sym,error=error)
!FROM UNMERGED     
!FROM UNMERGED     CALL timestop(handle)
!FROM UNMERGED   
!FROM UNMERGED   END SUBROUTINE hessian_diag_apply 
!FROM UNMERGED 
!FROM UNMERGED ! *****************************************************************************
!FROM UNMERGED !> \brief Serial code that prints matrices readable by Mathematica 
!FROM UNMERGED !> \param matrix - matrix to print
!FROM UNMERGED !> \param error ...
!FROM UNMERGED !> \par History
!FROM UNMERGED !>       2015.05 created [Rustam Z. Khaliullin]
!FROM UNMERGED !> \author Rustam Z. Khaliullin
!FROM UNMERGED ! *****************************************************************************
!FROM UNMERGED   SUBROUTINE print_mathematica_matrix(matrix,filename,error)
!FROM UNMERGED     
!FROM UNMERGED     TYPE(cp_dbcsr_type), INTENT(INOUT) :: matrix
!FROM UNMERGED     CHARACTER(len=*)                   :: filename
!FROM UNMERGED     TYPE(cp_error_type), INTENT(INOUT) :: error
!FROM UNMERGED 
!FROM UNMERGED     CHARACTER(len=*), PARAMETER :: routineN = 'print_mathematica_matrix', &
!FROM UNMERGED       routineP = moduleN//':'//routineN
!FROM UNMERGED 
!FROM UNMERGED     INTEGER :: hori_offset, vert_offset, Ncols, Nrows, col, row, jj, &
!FROM UNMERGED       nblkcols_tot, nblkrows_tot, ncores, unit_nr, handle, fiunit
!FROM UNMERGED     CHARACTER(LEN=20)                        :: formatstr, Scols
!FROM UNMERGED     CHARACTER(LEN=200)                       :: logfile
!FROM UNMERGED     INTEGER, ALLOCATABLE, DIMENSION(:)       :: ao_block_sizes, &
!FROM UNMERGED                                                 mo_block_sizes
!FROM UNMERGED     INTEGER, DIMENSION(:), POINTER           :: ao_blk_sizes, mo_blk_sizes
!FROM UNMERGED     LOGICAL                                  :: failure, found
!FROM UNMERGED     REAL(KIND=dp), ALLOCATABLE, &
!FROM UNMERGED       DIMENSION(:, :)                        :: H
!FROM UNMERGED     REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block_p
!FROM UNMERGED     TYPE(cp_dbcsr_type)                      :: matrix_asym
!FROM UNMERGED     TYPE(cp_logger_type), POINTER            :: logger
!FROM UNMERGED 
!FROM UNMERGED     CALL timeset(routineN,handle)
!FROM UNMERGED     
!FROM UNMERGED     ! get a useful output_unit
!FROM UNMERGED     logger => cp_error_get_logger(error)
!FROM UNMERGED     IF (logger%para_env%mepos==logger%para_env%source) THEN
!FROM UNMERGED        unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
!FROM UNMERGED     ELSE
!FROM UNMERGED        unit_nr=-1
!FROM UNMERGED     ENDIF
!FROM UNMERGED 
!FROM UNMERGED     ! serial code only
!FROM UNMERGED     ncores = dbcsr_mp_numnodes(dbcsr_distribution_mp(&
!FROM UNMERGED                 cp_dbcsr_distribution(matrix)))
!FROM UNMERGED     IF (ncores.gt.1) THEN                            
!FROM UNMERGED        CPErrorMessage(cp_failure_level,routineP,"serial code only",error)
!FROM UNMERGED        CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
!FROM UNMERGED     ENDIF
!FROM UNMERGED 
!FROM UNMERGED     nblkrows_tot = cp_dbcsr_nblkrows_total(matrix)
!FROM UNMERGED     nblkcols_tot = cp_dbcsr_nblkcols_total(matrix)
!FROM UNMERGED     CPPostcondition(nblkrows_tot==nblkcols_tot,cp_failure_level,routineP,error,failure)
!FROM UNMERGED     mo_blk_sizes => array_data(cp_dbcsr_col_block_sizes(matrix))
!FROM UNMERGED     ao_blk_sizes => array_data(cp_dbcsr_row_block_sizes(matrix))
!FROM UNMERGED     ALLOCATE(mo_block_sizes(nblkcols_tot),ao_block_sizes(nblkcols_tot))
!FROM UNMERGED     mo_block_sizes(:)=mo_blk_sizes(:)
!FROM UNMERGED     ao_block_sizes(:)=ao_blk_sizes(:)
!FROM UNMERGED 
!FROM UNMERGED     CALL cp_dbcsr_init(matrix_asym,error=error)
!FROM UNMERGED     CALL cp_dbcsr_create(matrix_asym,&
!FROM UNMERGED             template=matrix,&
!FROM UNMERGED             matrix_type=dbcsr_type_no_symmetry,error=error)
!FROM UNMERGED     CALL cp_dbcsr_desymmetrize(matrix,matrix_asym,error=error)
!FROM UNMERGED 
!FROM UNMERGED     Ncols=SUM(mo_block_sizes)
!FROM UNMERGED     Nrows=SUM(ao_block_sizes)
!FROM UNMERGED     ALLOCATE(H(Nrows,Ncols))
!FROM UNMERGED     H(:,:)=0.0_dp
!FROM UNMERGED     
!FROM UNMERGED     hori_offset=0
!FROM UNMERGED     DO col = 1, nblkcols_tot
!FROM UNMERGED 
!FROM UNMERGED        vert_offset=0
!FROM UNMERGED        DO row = 1, nblkrows_tot
!FROM UNMERGED           
!FROM UNMERGED           CALL cp_dbcsr_get_block_p(matrix_asym, row, col, block_p, found)
!FROM UNMERGED           IF (found) THEN
!FROM UNMERGED       
!FROM UNMERGED              H(vert_offset+1:vert_offset+ao_block_sizes(row),&
!FROM UNMERGED                hori_offset+1:hori_offset+mo_block_sizes(col))&
!FROM UNMERGED                =block_p(:,:)
!FROM UNMERGED 
!FROM UNMERGED           ENDIF
!FROM UNMERGED           
!FROM UNMERGED           vert_offset=vert_offset+ao_block_sizes(row)
!FROM UNMERGED    
!FROM UNMERGED        ENDDO
!FROM UNMERGED 
!FROM UNMERGED        hori_offset=hori_offset+mo_block_sizes(col)
!FROM UNMERGED 
!FROM UNMERGED     ENDDO ! loop over electron blocks 
!FROM UNMERGED 
!FROM UNMERGED     CALL cp_dbcsr_release(matrix_asym,error=error)
!FROM UNMERGED    
!FROM UNMERGED     IF (unit_nr>0) THEN
!FROM UNMERGED        logfile=TRIM(ADJUSTL(filename))
!FROM UNMERGED        OPEN (fiunit,file=logfile,status='REPLACE')
!FROM UNMERGED        WRITE(Scols,"(I10)") Ncols
!FROM UNMERGED        formatstr="("//TRIM(Scols)//"E27.17)"
!FROM UNMERGED        DO jj=1, Nrows
!FROM UNMERGED           WRITE(fiunit,formatstr) H(jj,:)
!FROM UNMERGED        ENDDO
!FROM UNMERGED        CLOSE(fiunit)
!FROM UNMERGED     ENDIF
!FROM UNMERGED 
!FROM UNMERGED     DEALLOCATE(mo_block_sizes)
!FROM UNMERGED     DEALLOCATE(ao_block_sizes)
!FROM UNMERGED     DEALLOCATE(H)
!FROM UNMERGED 
!FROM UNMERGED     CALL timestop(handle)
!FROM UNMERGED   
!FROM UNMERGED   END SUBROUTINE print_mathematica_matrix 

END MODULE almo_scf_optimizer

