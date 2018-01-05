!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines for ALMO-based electron correlation methods
!> \par History
!>       2014.11 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin 
! *****************************************************************************
MODULE almo_correlation
  USE almo_scf_types,                  ONLY: almo_scf_env_type
  USE mp2,                             ONLY: mp2_main
  USE cp_dbcsr_interface,              ONLY: &
       array_data, cp_dbcsr_add, cp_dbcsr_add_on_diag, &
       cp_dbcsr_col_block_sizes, cp_dbcsr_copy, cp_dbcsr_create, &
       cp_dbcsr_desymmetrize, cp_dbcsr_distribution, cp_dbcsr_filter, &
       cp_dbcsr_finalize, cp_dbcsr_frobenius_norm, cp_dbcsr_get_block_p, &
       cp_dbcsr_get_diag, cp_dbcsr_get_info, cp_dbcsr_init, &
       cp_dbcsr_iterator, cp_dbcsr_iterator_blocks_left, &
       cp_dbcsr_iterator_next_block, cp_dbcsr_iterator_start, &
       cp_dbcsr_iterator_stop, cp_dbcsr_multiply, cp_dbcsr_nblkcols_total, &
       cp_dbcsr_nblkrows_total, cp_dbcsr_norm, cp_dbcsr_print, &
       cp_dbcsr_release, cp_dbcsr_reserve_block2d, cp_dbcsr_row_block_sizes, &
       cp_dbcsr_scale, cp_dbcsr_set, cp_dbcsr_set_diag, cp_dbcsr_transposed, &
       cp_dbcsr_type, cp_dbcsr_work_create, dbcsr_distribution_mp, &
       dbcsr_error_type, dbcsr_init_random, dbcsr_mp_group, dbcsr_mp_mynode, &
       dbcsr_mp_numnodes, dbcsr_norm_maxabsnorm, dbcsr_type_no_symmetry, &
       dbcsr_type_symmetric, cp_dbcsr_get_stored_coordinates
  USE domain_submatrix_methods,        ONLY: &
       add_submatrices, construct_dbcsr_from_submatrices, &
       construct_submatrices, copy_submatrices, copy_submatrix_data, &
       init_submatrices, multiply_submatrices, print_submatrices, &
       release_submatrices, dsyev_submatrices
  USE domain_submatrix_types,          ONLY: domain_map_type,&
                                             domain_submatrix_type,&
                                             select_row,&
                                             select_row_col
  USE bibliography,                    ONLY: Khaliullin2013,&
                                             cite_reference
  USE cp_blacs_env,                    ONLY: cp_blacs_env_release,&
                                             cp_blacs_env_retain
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE input_constants,                 ONLY: &
       almo_constraint_distance, almo_deloc_none, almo_deloc_scf, &
       almo_deloc_x, almo_deloc_x_then_scf, almo_deloc_xalmo_1diag, &
       almo_deloc_xalmo_scf, almo_deloc_xalmo_x, almo_deloc_xk, &
       almo_domain_layout_atomic, almo_domain_layout_molecular, &
       almo_mat_distr_atomic, almo_mat_distr_molecular, almo_scf_diag, &
       almo_scf_dm_sign, almo_scf_pcg, cg_hager_zhang, do_bondparm_vdw, &
       optimizer_diis, optimizer_pcg, tensor_orthogonal, virt_full, &
       virt_minimal, virt_number, virt_occ_size, xalmo_case_block_diag, &
       xalmo_case_fully_deloc, xalmo_case_normal
  USE kinds,                           ONLY: dp
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'almo_correlation'

  PUBLIC :: almo_correlation_main, compute_pseudocanonical_orbitals

CONTAINS

! *****************************************************************************
!> \brief Entry point for ALMO-based electron correlation methods
!> \param qs_env ...
!> \param almo_scf_env ...
!> \param error ...
!> \par History
!>       2014.11 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_correlation_main(qs_env,almo_scf_env,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(almo_scf_env_type)                  :: almo_scf_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'almo_correlation_main', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: unit_nr, handle
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: p_new_block
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)
    failure=.FALSE.

    ! get a useful output_unit
    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    !CALL compute_pseudocanonical_orbitals(almo_scf_env,error)
    CALL mp2_main(qs_env=qs_env,almo_scf_env=almo_scf_env,&
            calc_forces=.FALSE.,error=error)

    CALL timestop(handle)

  END SUBROUTINE almo_correlation_main

! *****************************************************************************
!> \brief Computes pseudocanonical occupied and virtual orbitals 
!> \param almo_scf_env ...
!> \param error ...
!> \par History
!>       2014.11 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE compute_pseudocanonical_orbitals(almo_scf_env,error)

    TYPE(almo_scf_env_type)                  :: almo_scf_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'compute_pseudocanonical_orbitals', &
      routineP = moduleN//':'//routineN

    TYPE(cp_logger_type), POINTER            :: logger
    LOGICAL                                  :: failure, tr
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: diagonal 
    REAL(kind=dp), DIMENSION(:, :), POINTER  :: p_new_block
    TYPE(cp_dbcsr_type)                      :: &
                                                m_tmp_no_1, &
                                                STsiginv, &
                                                matrix_no_full, &
                                                matrix_v_psc, &
                                                matrix_t_psc, &
                                                m_tmp_nn_1,&
                                                m_tmp_nn_2
    TYPE(domain_map_type)                    :: domain_map_full
    INTEGER                                  :: ndomains, idomain, jdomain, &
                                                index1_prev, ispin, dim0, &
                                                mynode, nblkrows_tot, &
                                                nblkcols_tot, row, col, &
                                                iblock_row, iblock_col, &
                                                nocc_of_block, GroupID, &
                                                hold, handle, unit_nr
    TYPE(domain_submatrix_type), &
              ALLOCATABLE, DIMENSION(:)      :: subm_v_psc, &
                                                subm_t_psc, &
                                                subm_v1
    
    CALL timeset(routineN,handle)
    failure=.FALSE.

    ! get a useful output_unit
    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    ndomains = almo_scf_env%ndomains

    ! create a neighbor map with all domains as neighbors
    ! RZK-warning: quadratic algorithm
    ALLOCATE(domain_map_full%index1(ndomains))
    ALLOCATE(domain_map_full%pairs(ndomains*ndomains,2))
    DO idomain=1, ndomains
       domain_map_full%index1(idomain)=ndomains*idomain+1
       DO jdomain=1, ndomains
          domain_map_full%pairs(jdomain+(idomain-1)*ndomains,1)=jdomain
          domain_map_full%pairs(jdomain+(idomain-1)*ndomains,2)=idomain
       ENDDO
    ENDDO
    ! check if the map is alright
    !DO idomain=1, ndomains
    !   IF (idomain.eq.1) THEN
    !      index1_prev = 1
    !   ELSE
    !      index1_prev = domain_map_full%index1(idomain-1)
    !   ENDIF
    !   IF (unit_nr>0) THEN
    !      WRITE(unit_nr,'(T2,I13,":")') idomain
    !      WRITE(unit_nr,'(T15,5I13)') &
    !         domain_map_full%pairs(index1_prev:domain_map_full%index1(idomain)-1,1) ! includes self
    !   ENDIF
    !ENDDO

    DO ispin=1,almo_scf_env%nspins

       ! RZK-warning: quadratic algorithm
       ! create a template of an NxO matrix with ALL blocks 
       CALL cp_dbcsr_init(matrix_no_full,error=error)
       CALL cp_dbcsr_create(matrix_no_full,template=almo_scf_env%matrix_t(ispin),&
               error=error)
       CALL cp_dbcsr_work_create(matrix_no_full,work_mutable=.TRUE.,&
               error=error)
       mynode = dbcsr_mp_mynode(dbcsr_distribution_mp(&
          cp_dbcsr_distribution(almo_scf_env%matrix_t(ispin))))
       nblkrows_tot = cp_dbcsr_nblkrows_total(almo_scf_env%matrix_t(ispin))
       nblkcols_tot = cp_dbcsr_nblkcols_total(almo_scf_env%matrix_t(ispin))
       DO row = 1, nblkrows_tot
          DO col = 1, nblkcols_tot
             tr = .FALSE.
             iblock_row = row
             iblock_col = col
             CALL cp_dbcsr_get_stored_coordinates(matrix_no_full,&
                     iblock_row, iblock_col, tr, hold)
             IF(hold.EQ.mynode) THEN
                NULLIFY (p_new_block)
                CALL cp_dbcsr_reserve_block2d(matrix_no_full,iblock_row,iblock_col,p_new_block)
                nocc_of_block=SIZE(p_new_block,2)
                CPPostcondition(ASSOCIATED(p_new_block),cp_failure_level,routineP,error,failure)
                CPPrecondition(nocc_of_block.gt.0,cp_failure_level,routineP,error,failure)
                p_new_block(:,:) = HUGE(1.0_dp)
             ENDIF
          ENDDO
       ENDDO
       CALL cp_dbcsr_finalize(matrix_no_full,error=error)
       CALL cp_dbcsr_filter(matrix_no_full,almo_scf_env%eps_filter,error=error)
       !CALL cp_dbcsr_print(matrix_no_full,error=error)

       ! create all necessary matrices and submatrices
       CALL cp_dbcsr_init(m_tmp_no_1,error=error)
       CALL cp_dbcsr_create(m_tmp_no_1,template=almo_scf_env%matrix_t(ispin),&
               matrix_type=dbcsr_type_no_symmetry,error=error)
       CALL cp_dbcsr_init(STsiginv,error=error)
       CALL cp_dbcsr_create(STsiginv,template=almo_scf_env%matrix_t(ispin),&
               matrix_type=dbcsr_type_no_symmetry,error=error)
       CALL cp_dbcsr_init(matrix_t_psc,error=error)
       CALL cp_dbcsr_create(matrix_t_psc,template=almo_scf_env%matrix_t(ispin),&
               matrix_type=dbcsr_type_no_symmetry,error=error)
       CALL cp_dbcsr_init(matrix_v_psc,error=error)
       CALL cp_dbcsr_create(matrix_v_psc,template=almo_scf_env%matrix_s(1),&
               matrix_type=dbcsr_type_no_symmetry,error=error)
       CALL cp_dbcsr_init(m_tmp_nn_1,error=error)
       CALL cp_dbcsr_create(m_tmp_nn_1,template=almo_scf_env%matrix_s(1),&
               matrix_type=dbcsr_type_no_symmetry,error=error)
       CALL cp_dbcsr_init(m_tmp_nn_2,error=error)
       CALL cp_dbcsr_create(m_tmp_nn_2,template=almo_scf_env%matrix_s(1),&
               matrix_type=dbcsr_type_no_symmetry,error=error)
       ALLOCATE(subm_t_psc(ndomains))
       ALLOCATE(subm_v_psc(ndomains))
       ALLOCATE(subm_v1(ndomains))
       CALL init_submatrices(subm_t_psc,error)
       CALL init_submatrices(subm_v_psc,error)
       CALL init_submatrices(subm_v1,error)
       
       ! construct occupied pseudo-canonical orbitals
       CALL cp_dbcsr_multiply("N","N",1.0_dp,almo_scf_env%matrix_t(ispin),&
               almo_scf_env%matrix_sigma_inv(ispin),0.0_dp,m_tmp_no_1,&
               filter_eps=almo_scf_env%eps_filter,error=error)
       ! convert T.siginv into submatrix using the full NO template
       CALL construct_submatrices(&
               matrix=m_tmp_no_1,&
               submatrix=subm_t_psc,&
               distr_pattern=matrix_no_full,&
               domain_map_row=domain_map_full,&
               node_of_domain=almo_scf_env%cpu_of_domain,&
               !job_type=select_row,&
               error=error)
       CALL multiply_submatrices('N','N',1.0_dp,&
               subm_t_psc,almo_scf_env%domain_ks_xx_oo_evec(:,ispin),&
               0.0_dp,almo_scf_env%domain_t_psc(:,ispin),error)

       !check whether these are truly pseudo-canonical orbitals
       !CALL construct_submatrices(&
       !        matrix=almo_scf_env%matrix_ks(ispin),&
       !        submatrix=subm_v1,&
       !        ! it is allowed to use NxO matrix as a distr_pattern
       !        ! because its size does not matter, only block distr info is used
       !        distr_pattern=matrix_no_full,&
       !        domain_map_row=domain_map_full,&
       !        domain_map_col=domain_map_full,&
       !        node_of_domain=almo_scf_env%cpu_of_domain,&
       !        !job_type=select_row_col,&
       !        error=error)
       !CALL multiply_submatrices('N','N',1.0_dp,&
       !        subm_v1,almo_scf_env%domain_t_psc(:,ispin),&
       !        0.0_dp,subm_t_psc,error)
       !CALL multiply_submatrices('T','N',1.0_dp,&
       !        almo_scf_env%domain_t_psc(:,ispin),subm_t_psc,&
       !        0.0_dp,subm_v1,error)
       !GroupID = dbcsr_mp_group(dbcsr_distribution_mp(cp_dbcsr_distribution(matrix_no_full)))
       !CALL print_submatrices(subm_v1,GroupID,error=error)
        
       ! use projected AOs as virtuals V = (I-RS)    
       CALL cp_dbcsr_multiply("N","N",1.0_dp,&
               almo_scf_env%matrix_s(1),&
               m_tmp_no_1,&
               0.0_dp,STsiginv,&
               filter_eps=almo_scf_env%eps_filter,&
               error=error)
       ! setup identity matrix
       !CALL cp_dbcsr_get_info( almo_scf_env%matrix_s(1), nfullrows_total=dim0 )
       !ALLOCATE(diagonal(dim0))
       !diagonal(:)=1.0_dp
       CALL cp_dbcsr_set(matrix_v_psc,0.0_dp,error=error)
       CALL cp_dbcsr_add_on_diag(matrix_v_psc,1.0_dp,error=error)
       !CALL cp_dbcsr_set_diag(matrix_v_psc,diagonal,error=error)
       !DEALLOCATE(diagonal)
       CALL cp_dbcsr_filter(matrix_v_psc,almo_scf_env%eps_filter,error=error)
       CALL cp_dbcsr_multiply("N","T",-1.0_dp,&
               almo_scf_env%matrix_t(ispin),&
               STsiginv,&
               1.0_dp,matrix_v_psc,&
               filter_eps=almo_scf_env%eps_filter,&
               error=error)
       
       ! convert V into submatrix using the full NO template
       CALL construct_submatrices(&
               matrix=matrix_v_psc,&
               submatrix=subm_v_psc,&
               ! it is allowed to use NxO matrix as a distr_pattern
               ! because its size does not matter, only block distr info is used
               distr_pattern=matrix_no_full,&
               domain_map_row=domain_map_full,&
               domain_map_col=almo_scf_env%domain_map(ispin),&
               node_of_domain=almo_scf_env%cpu_of_domain,&
               !job_type=select_row,&
               error=error)

!GroupID = dbcsr_mp_group(dbcsr_distribution_mp(cp_dbcsr_distribution(matrix_no_full)))
!CALL print_submatrices(subm_v_psc,GroupID,error=error)

       CALL multiply_submatrices('N','N',1.0_dp,&
               subm_v_psc,almo_scf_env%domain_s_sqrt_inv(:,ispin),&
               0.0_dp,subm_v1,error)
       CALL multiply_submatrices('N','N',1.0_dp,&
               subm_v1,almo_scf_env%domain_ks_xx_vv_evec(:,ispin),&
               0.0_dp,almo_scf_env%domain_v_psc(:,ispin),error)

       ! check whether these are truly pseudo-canonical orbitals
       !CALL construct_submatrices(&
       !        matrix=almo_scf_env%matrix_ks(ispin),&
       !        submatrix=subm_v1,&
       !        ! it is allowed to use NxO matrix as a distr_pattern
       !        ! because its size does not matter, only block distr info is used
       !        distr_pattern=matrix_no_full,&
       !        domain_map_row=domain_map_full,&
       !        domain_map_col=domain_map_full,&
       !        node_of_domain=almo_scf_env%cpu_of_domain,&
       !        !job_type=select_row_col,&
       !        error=error)
       !CALL multiply_submatrices('N','N',1.0_dp,&
       !        subm_v1,almo_scf_env%domain_v_psc(:,ispin),&
       !        0.0_dp,subm_v_psc,error)
       !CALL multiply_submatrices('T','N',1.0_dp,&
       !        almo_scf_env%domain_v_psc(:,ispin),subm_v_psc,&
       !        0.0_dp,subm_v1,error)
       !GroupID = dbcsr_mp_group(dbcsr_distribution_mp(cp_dbcsr_distribution(matrix_no_full)))
       !CALL print_submatrices(subm_v1,GroupID,error=error)
        
       ! release all necessary matrices and submatrices
       CALL cp_dbcsr_release(m_tmp_no_1,error=error)
       CALL cp_dbcsr_release(STsiginv,error=error)
       CALL cp_dbcsr_release(matrix_t_psc,error=error)
       CALL cp_dbcsr_release(matrix_v_psc,error=error)
       CALL cp_dbcsr_release(m_tmp_nn_1,error=error)
       CALL cp_dbcsr_release(m_tmp_nn_2,error=error)
       CALL cp_dbcsr_release(matrix_no_full,error=error)
       CALL release_submatrices(subm_t_psc,error)
       CALL release_submatrices(subm_v_psc,error)
       CALL release_submatrices(subm_v1,error)
       DEALLOCATE(subm_t_psc)
       DEALLOCATE(subm_v_psc)
       DEALLOCATE(subm_v1)
       
    ENDDO ! spins

    DEALLOCATE(domain_map_full%pairs)
    DEALLOCATE(domain_map_full%index1)

    CALL timestop(handle)
  
  END SUBROUTINE compute_pseudocanonical_orbitals

END MODULE almo_correlation

