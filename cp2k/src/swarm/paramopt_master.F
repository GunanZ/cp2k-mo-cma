!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Master's routines for global optimization
!> \author Ole Schuett
! *****************************************************************************
MODULE paramopt_master
  USE cmaes,                           ONLY: cmaes_input,&
                                             cmaes_type,&
                                             cmaes_run,&
                                             cmaes_init
  USE mocmaes,                         ONLY: mo_cmaes_type,&
                                             mo_cmaes_type_p,&
                                             mo_cmaes_init,&
                                             mo_cmaes_run   
  USE cp_files,                        ONLY: open_file,&
                                             close_file
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_units,                        ONLY: cp_unit_from_cp2k
  USE input_constants,                 ONLY: dump_xmol,&
                                             paramopt_do_cma,&
                                             paramopt_do_powell,&
                                             paramopt_do_grid,&
                                             paramopt_do_mocma
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_release,&
                                             section_vals_retain,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             int_8
  USE machine,                         ONLY: default_output_unit,&
                                             m_getcwd
  USE swarm_message,                   ONLY: swarm_message_get,&
                                             swarm_message_type,&
                                             swarm_message_add
  USE powell,                   ONLY: powell_optimize, &
                                             opt_state_type
  USE parallel_rng_types,              ONLY: rng_stream_type, &
                                             create_rng_stream, &
                                             delete_rng_stream,&
                                             next_random_number, &
                                             GAUSSIAN, &
                                             UNIFORM    
#include "../common/cp_common_uses.f90"

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'paramopt_master'

 PUBLIC :: paramopt_master_type
 PUBLIC :: paramopt_master_init, paramopt_master_finalize
 PUBLIC :: paramopt_master_steer


 TYPE paramopt_master_type
   PRIVATE
   INTEGER                                             :: iw = 0
   INTEGER(int_8)                                      :: count_reports = 0
   INTEGER                                             :: method
   INTEGER                                             :: num_ref_mol = 1
   INTEGER                                             :: num_input_rcvd
   INTEGER                                             :: i_iteration = 0
   TYPE(section_vals_type), POINTER                    :: paramopt_section => Null()
   TYPE(cp_error_type)                                 :: error
   !--------------------------
   REAL(KIND=dp),ALLOCATABLE,DIMENSION(:,:,:)          :: energies
   REAL(KIND=dp),ALLOCATABLE,DIMENSION(:,:,:)          :: rmsd
   REAL(KIND=dp),ALLOCATABLE,DIMENSION(:)              :: mcharge
   INTEGER,ALLOCATABLE,DIMENSION(:,:,:)                :: converged
   REAL(KIND=dp),ALLOCATABLE,DIMENSION(:,:)            :: ref_matrix
   INTEGER                                             :: database=1
   REAL(KIND=dp),DIMENSION(:,:),POINTER                :: x
   REAL(KIND=dp),DIMENSION(:,:),POINTER                :: f
   REAL(KIND=dp),DIMENSION(:),POINTER                  :: u_bound
   REAL(KIND=dp),DIMENSION(:),POINTER                  :: l_bound
   INTEGER                                             :: n,m,popsize
   REAL(KIND=dp)                                       :: step_size
   REAL(KIND=dp)                                       :: penalty=0._dp
   LOGICAL                                             :: randinput
   !--------------------------
   TYPE(cmaes_input)                                   :: cmaes_inp
   CLASS(cmaes_type),ALLOCATABLE,DIMENSION(:)          :: population
   INTEGER,ALLOCATABLE,DIMENSION(:)                    :: offspring_list
   !--------------------------
   TYPE(opt_state_type)                                :: powell_param
   !-------------------------
   LOGICAL                                             :: benchmark
   INTEGER                                             :: bench_func
   INTEGER                                             :: current_inp
 END TYPE paramopt_master_type

 CONTAINS


! *****************************************************************************
!> \brief Initializes the master of the global optimizer
!> \param this ...
!> \param para_env ...
!> \param root_section ...
!> \param input_path ...
!> \param n_walkers ...
!> \param iw ...
!> \param error ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE paramopt_master_init(this, para_env, root_section, input_path, n_walkers, iw, error)
    TYPE(paramopt_master_type), INTENT(INOUT)  :: this
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(section_vals_type), POINTER         :: powell_section
    TYPE(section_vals_type), POINTER         :: global_section
    CHARACTER(LEN=*), INTENT(IN)             :: input_path
    INTEGER, INTENT(IN)                      :: n_walkers, iw
    TYPE(cp_error_type), INTENT(inout)       :: error
    CHARACTER(len=*), PARAMETER :: routineN = 'paramopt_master_steer', &
      routineP = moduleN//':'//routineN
    TYPE(cp_logger_type), POINTER            :: logger
    INTEGER                                  :: n,m
    INTEGER                                  :: i
    INTEGER                                  :: j
    INTEGER                                  :: popsize=1
    LOGICAL                                  :: explicit_key
    REAL(KIND=dp),DIMENSION(:),POINTER       :: init_x
    !------------------
    TYPE(rng_stream_type),POINTER            :: rng_uniform_stream => Null()
    INTEGER                                  :: global_seed
    TYPE(mo_cmaes_type),DIMENSION(:),POINTER   :: bla
    TYPE(mo_cmaes_type_p),DIMENSION(:),ALLOCATABLE   :: foo
    TYPE(mo_cmaes_type_p),DIMENSION(:),ALLOCATABLE   :: moo
    REAL(kind=dp)                                    :: x,y

    NULLIFY(logger)

    this%iw = iw
    this%error = error
    !------------------------------------
    this%paramopt_section => section_vals_get_subs_vals(root_section, "SWARM%PARAM_OPT", error=error)
    global_section => section_vals_get_subs_vals(root_section, "GLOBAL", error=error)
    CALL section_vals_val_get(global_section,"SEED", i_val=global_seed, error=error)
    powell_section => section_vals_get_subs_vals(root_section, "ATOM%POWELL", error=error)
    CALL section_vals_retain(this%paramopt_section,error=error)
    CALL section_vals_val_get(this%paramopt_section,"NUM_VAR", i_val=this%n, error=error)
    CALL section_vals_val_get(this%paramopt_section,"NUM_OBJ", i_val=this%m, error=error)
    CALL section_vals_val_get(this%paramopt_section,"POPSIZE", i_val=this%popsize, error=error)
    CALL section_vals_val_get(this%paramopt_section,"REF_DATABASE", i_val=this%database, error=error)
    !------------------------------------
    n = this%n
    !popsize = 200
    !popsize = 62 
    !popsize = 2 
    popsize = this%popsize
    if (mod(popsize,2) .ne. 0) then
        STOP "popsize wrong"
    endif
    allocate(this%x(popsize,n))
    !------------------------------------
    CALL section_vals_val_get(this%paramopt_section,"UPPER_BOUNDS", r_vals=this%u_bound, error=error)
    CALL section_vals_val_get(this%paramopt_section,"LOWER_BOUNDS", r_vals=this%l_bound, error=error)
    CALL section_vals_val_get(this%paramopt_section,"METHOD", i_val=this%method, error=error)
    CALL section_vals_val_get(this%paramopt_section,"STEP_SIZE", r_val=this%step_size, error=error)
    CALL section_vals_val_get(this%paramopt_section,"BENCHMARK", l_val=this%benchmark, error=error)
    this%bench_func = 1
    !------------------------------------
    CALL section_vals_val_get(this%paramopt_section,"INITIAL_X", explicit=explicit_key, error=error)
    IF(explicit_key)THEN
        CALL section_vals_val_get(this%paramopt_section,"INITIAL_X", r_vals=init_x, error=error)
        this%x(1,:) = init_x
        this%randinput = .FALSE.
        CALL cp_assert( size(init_x) .eq. n ,cp_failure_level,cp_assertion_failed,&
                       routineP, "init_x does not match num_var",error)
    ELSE
        this%randinput = .TRUE.
        print *,'generating random x'
        CALL create_rng(rng_uniform_stream,global_seed)
        do j=1,2*popsize
            do i=1,n; this%x(j,i) = next_random_number(rng_uniform_stream,error=error); enddo 
            this%x(j,:) = var_trans(this%x(j,:),this%u_bound,this%l_bound)
         enddo
         print *,'lower bound',this%l_bound
         print *,'upper bound',this%u_bound
         print *,'done random x',this%x(1,:)
         print *,'done random x trans',var_back_trans(this%x(1,:),this%u_bound,this%l_bound)
         CALL delete_rng_stream(rng_uniform_stream,error=error)
    ENDIF
    !------------------------------------
    m = this%m
    ! this is now set from input
    !this%database = 3 ! 1 for the s22 and 2 for water_cluster and 3 for water27
    if (.not. this%benchmark) then 
        if (this%database == 1) then
            this%num_ref_mol = 3*22
            CALL read_ref_data(this%database,"s22/reference_data",this%ref_matrix) 
            allocate(this%energies(popsize,this%num_ref_mol/3,3))
            allocate(this%rmsd(popsize,this%num_ref_mol/3,3))
            allocate(this%converged(popsize,this%num_ref_mol/3,3))
            allocate(this%mcharge(popsize))
        else if (this%database == 2) then
            this%num_ref_mol = 39
            CALL read_ref_data(this%database,"water_clusters/reference_data",this%ref_matrix) 
            allocate(this%energies(popsize,this%num_ref_mol,1))
            allocate(this%rmsd(popsize,this%num_ref_mol,1))
            allocate(this%converged(popsize,this%num_ref_mol,1))
            allocate(this%mcharge(popsize))
        else if (this%database == 3) then
            this%num_ref_mol = 30 
            CALL read_ref_data(this%database,"WATER27structures/reference_data",this%ref_matrix) 
            allocate(this%energies(popsize,this%num_ref_mol,1))
            allocate(this%rmsd(popsize,this%num_ref_mol,1))
            allocate(this%converged(popsize,this%num_ref_mol,1))
            allocate(this%mcharge(popsize))
        else if (this%database == 4) then
            this%num_ref_mol = 7  
            CALL read_ref_data(this%database,"atomic_structures/reference_data",this%ref_matrix) 
            allocate(this%energies(popsize,this%num_ref_mol,1))
            allocate(this%rmsd(popsize,this%num_ref_mol,1))
            allocate(this%converged(popsize,this%num_ref_mol,1))
            allocate(this%mcharge(popsize))
        else if (this%database == 5) then
            this%num_ref_mol = 35 
            CALL read_ref_data(this%database,"atomic_structures/reference_data",this%ref_matrix) 
            allocate(this%energies(popsize,this%num_ref_mol,1))
            allocate(this%rmsd(popsize,this%num_ref_mol,1))
            allocate(this%converged(popsize,this%num_ref_mol,1))
            allocate(this%mcharge(popsize))
        endif
    else 
        allocate(this%energies(popsize,this%num_ref_mol,1))
        allocate(this%rmsd(popsize,this%num_ref_mol,1))
        allocate(this%converged(popsize,this%num_ref_mol,1))
    endif
    this%num_input_rcvd = 0
    this%current_inp = 1
    allocate(this%f(popsize,m))
    logger => cp_error_get_logger(error)
    ! ref matrix 
    ! Nr.  Dimer monomer1 monomer2 ref_energy
    !do, i=1,size(ref_matrix,dim=1)
    !     print *, ( ref_matrix(i,j), j=1,size(ref_matrix,dim=2) )
    !enddo 
    SELECT CASE (this%method)
       CASE(paramopt_do_cma)
            allocate(this%cmaes_inp%l_bound(n))
            allocate(this%cmaes_inp%u_bound(n))
            this%cmaes_inp%l_bound =0._dp
            this%cmaes_inp%u_bound =1._dp
            this%cmaes_inp%step_size = this%step_size
            this%cmaes_inp%n = n
            this%cmaes_inp%randinput = this%randinput
            this%cmaes_inp%input_seed = global_seed
            this%cmaes_inp%m = m
            this%cmaes_inp%x => this%x(1,:)
            this%cmaes_inp%f => this%f(1,:)
            !this%cmaes_inp%l_bound => this%l_bound
            !this%cmaes_inp%u_bound => this%u_bound
            ALLOCATE(cmaes_type::this%population(1))
            call this%population(1)%cmaes_init(this%cmaes_inp)
        CASE(paramopt_do_powell)
!            allocate(this%powell_param%xu(n))
!            allocate(this%powell_param%xl(n))
!            this%powell_param%xl = 0._dp
!            this%powell_param%xu = 1._dp
!            this%powell_param%xl = this%l_bound
!            this%powell_param%xu = this%u_bound
!            this%powell_param%f => this%f(1,1)
!            this%powell_param%x => this%x(1,:)
            this%powell_param%f => this%f(1,1)
            this%powell_param%xopt => this%x(1,:)
            CALL init_powell_var(this,this%powell_param,powell_section,error)
            this%x(1,:) = var_back_trans(this%x(1,:),this%u_bound,this%l_bound)
            CALL powell_optimize (this%n, this%x(1,:), &
                                  this%powell_param)
            !CALL powell_optimize (this%n, this%x(1,:) ,this%powell_param%xl, &
            !                      this%powell_param%xu,&
            !                      this%powell_param)
            this%x(1,:) = var_trans(this%x(1,:),this%u_bound,this%l_bound)
        CASE(paramopt_do_mocma)
            allocate(this%cmaes_inp%l_bound(n))
            allocate(this%cmaes_inp%u_bound(n))
            this%cmaes_inp%l_bound =0._dp
            this%cmaes_inp%u_bound =1._dp
            this%cmaes_inp%step_size = this%step_size
            this%cmaes_inp%n = n
            this%cmaes_inp%randinput = this%randinput
            this%cmaes_inp%input_seed = global_seed
            this%cmaes_inp%m = m
            this%cmaes_inp%x => this%x(1,:)
            this%cmaes_inp%f => this%f(1,:)
            !this%cmaes_inp%l_bound => this%l_bound
            !this%cmaes_inp%u_bound => this%u_bound
            ALLOCATE(mo_cmaes_type::this%population(popsize))
            allocate(this%offspring_list(popsize))
            do i=1,popsize
                this%offspring_list(i) = i
            enddo
            select type (pop => this%population)
                type is (mo_cmaes_type)
                CALL mo_cmaes_init(pop,this%cmaes_inp)
                if (.not. this%randinput) then
                    pop(1)%x = var_back_trans(pop(1)%x,this%u_bound,this%l_bound)
                    do i=1,size(pop)
                        pop(i)%f => this%f(i,:)
                        pop(i)%f = huge(1._dp)
                        pop(i)%x => this%x(i,:)
                        pop(i)%x = pop(1)%x
                        print *,'x',pop(i)%x
                    enddo
                else
                    do i=1,size(pop)
                        pop(i)%f => this%f(i,:)
                        pop(i)%f = huge(1._dp)
                        pop(i)%x => this%x(i,:)
                        pop(i)%x = var_back_trans(pop(i)%x,this%u_bound,this%l_bound)
                        print *,'x',pop(i)%x
                    enddo
                endif
            end select
       CASE(paramopt_do_grid)
            this%x(1,:) = var_back_trans(this%x(1,:),this%u_bound,this%l_bound)           
       CASE DEFAULT
        STOP "Unkown glbopt_method"
    END SELECT
    print *,'paramopt_master_init done'
 END SUBROUTINE paramopt_master_init

! *****************************************************************************
!> \brief Central steering routine of global optimizer
!> \param this ...
!> \param report ...
!> \param cmd ...
!> \param should_stop ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE paramopt_master_steer(this, report, cmd, should_stop)
    TYPE(paramopt_master_type), INTENT(INOUT)  :: this
    TYPE(swarm_message_type)                   :: report, cmd
    CHARACTER(len=default_string_length)       :: report_status
    LOGICAL, INTENT(INOUT)                     :: should_stop


    SELECT CASE (this%method)
        CASE(paramopt_do_cma,paramopt_do_powell)
            CALL do_optim(this,report,cmd,this%method,should_stop)
        CASE(paramopt_do_grid)
            CALL do_grid_search(this,report,cmd,should_stop)
        CASE(paramopt_do_mocma)
            CALL do_mocma(this,report,cmd,should_stop)
        CASE DEFAULT
            STOP "Unkown paramopt_method"
    END SELECT

   ! IF(this%E_lowest<this%E_target) THEN
   !    IF(this%iw>0)  WRITE (this%iw,"(A,I8,A)" )  &
   !       " GLBOPT| Reached E_pot < E_target after ", this%i_iteration," iterations. Quitting."
   !    should_stop = .TRUE.
   ! END IF
 END SUBROUTINE paramopt_master_steer

! *****************************************************************************
!> \brief Helper routine for glbopt_master_steer(), updates stats, etc.
!> \param this ...
!> \param report ...
!> \author Ole Schuett
! *****************************************************************************
SUBROUTINE do_grid_search(this,report,cmd,should_stop)
    TYPE(paramopt_master_type), INTENT(INOUT)  :: this
    TYPE(swarm_message_type)                   :: report, cmd
    CHARACTER(len=default_string_length)       :: report_status
    CHARACTER(len=default_string_length)       :: input_file
    LOGICAL, INTENT(INOUT)                     :: should_stop
    INTEGER                                    :: wid
    INTEGER                                    :: input
    REAL(kind=dp)                              :: myresult
    INTEGER                                    :: charge=0
    INTEGER                                    :: i,j 
    INTEGER                                    :: iounit 
    INTEGER                                    :: id=1
    INTEGER,DIMENSION(3)                       :: ind
    REAL(kind=dp),DIMENSION(:),ALLOCATABLE     :: func_array
    REAL(kind=dp),DIMENSION(:),POINTER         :: func
    allocate(func_array(size(this%energies,dim=2)))
          CALL swarm_message_get(report, "worker_id", wid)
          CALL swarm_message_get(report, "status", report_status)
          if (trim(report_status) .eq. "done") then
            this%num_input_rcvd = this%num_input_rcvd + 1
            CALL swarm_message_get(report, "pop_id",id)
            CALL swarm_message_get(report, "inp_num",input)
            CALL index2grid(input,(/1,size(this%energies,dim=2),size(this%energies,dim=3)/),ind)
            CALL swarm_message_get(report, "result", func)
            this%energies(id,ind(2),ind(3)) = func(1)
            this%rmsd(id,ind(2),ind(3))     = func(2) 
            CALL swarm_message_get(report, "converged", this%converged(id,ind(2),ind(3)))
            if (this%num_input_rcvd .eq. this%num_ref_mol) then
                CALL energy2function(this%num_ref_mol,this%energies(1,:,:), &
                                        this%ref_matrix,this%database,func_array)
                this%f(1,1) = sum(func_array)
                this%f(1,2) = sum(this%rmsd(1,:,:))
                this%f(1,3) = func_array(2)
                CALL open_file(file_name="output", file_status="UNKNOWN",&
                    file_action="WRITE", file_position="APPEND",unit_number=iounit)
                    write(iounit,*) var_trans(this%x(1,:),this%u_bound,this%l_bound), &
                                    this%f(1,:),func_array,this%rmsd(1,:,:)
                CALL close_file(iounit)
                print *,'done writing file'
                CALL grid_search(this) 
                !---reset----------
                this%current_inp = 0
                this%num_input_rcvd = 0
                this%i_iteration = this%i_iteration + 1
            endif
          endif
          if (this%current_inp .le. this%num_ref_mol .and. &
            (report_status .eq. "wait_done" .or. report_status .eq. "initial_hello")) then
            CALL swarm_message_add(cmd, "command", "energy")
            CALL swarm_message_add(cmd, "pop_id",1)
            CALL swarm_message_add(cmd, "inp_num", this%current_inp )
            CALL index2grid(this%current_inp,(/1,size(this%energies,dim=2), &
                             size(this%energies,dim=3)/),ind)
                             print *,ind
            CALL gen_input_file(ind(2),ind(3),this%database,input_file)
            CALL swarm_message_add(cmd, "inp_file", input_file)
            ! change charge if necessary
            ! w27 1-15 are 0 then 16-21 are 1 then 22-29 are -1 and 30 is 0
            if (this%database .eq. 3) then
                select case(ind(2))
                    case(1:15,30)
                        charge = 0
                    case(16:21)
                        charge = 1
                    case(22:29)
                        charge = -1
                end select
            endif
            if (this%database .eq. 5) then
                select case(ind(2))
                case(1:15,30:32,34:35)
                        charge = 0
                    case(16:21,33)
                        charge = 1
                    case(22:29)
                        charge = -1
                end select
            endif

            if (this%database .eq. 4) then
                select case(ind(2))
                    case(1:2,4)
                        charge = 0
                    case(3)
                        charge = 1
                end select
            endif
            CALL swarm_message_add(cmd, "charge", charge)
            CALL swarm_message_add(cmd, "x",&
            var_trans(this%x(1,:),this%u_bound,this%l_bound))
            this%current_inp = this%current_inp + 1
          else
            CALL swarm_message_add(cmd, "command", "wait")     
          endif
    deallocate(func_array)
END SUBROUTINE do_grid_search

! *****************************************************************************
!> \brief Helper routine for glbopt_master_steer(), updates stats, etc.
!> \param this ...
!> \param report ...
!> \author Ole Schuett
! *****************************************************************************
SUBROUTINE do_mocma(this,report,cmd,should_stop)
    TYPE(paramopt_master_type), INTENT(INOUT)  :: this
    TYPE(swarm_message_type)                   :: report, cmd
    CHARACTER(len=default_string_length)       :: report_status
    CHARACTER(len=default_string_length)       :: input_file
    LOGICAL, INTENT(INOUT)                     :: should_stop
    INTEGER                                    :: wid
    INTEGER                                    :: input
    REAL(kind=dp)                              :: myresult
    INTEGER                                    :: charge=0
    INTEGER                                    :: spinmult=1
    INTEGER                                    :: i,j 
    INTEGER                                    :: iounit 
    INTEGER                                    :: id=1 
    INTEGER                                    :: popsize,offspringsize 
    REAL(kind=dp),DIMENSION(:),POINTER         :: func
    REAL(kind=dp)                              :: penalty
    REAL(kind=dp),DIMENSION(:),ALLOCATABLE     :: func_array
    INTEGER                                    :: k,l  
    INTEGER,DIMENSION(3)                       :: ind
    popsize = size(this%population)
    offspringsize = int(0.5*popsize)
    allocate(func_array(size(this%energies,dim=2)))
    if (this%benchmark) then
          CALL swarm_message_get(report, "worker_id", wid)
          CALL swarm_message_get(report, "status", report_status)
          print *,report_status
        if (trim(report_status) .eq. "done") then
            this%num_input_rcvd = this%num_input_rcvd + 1
            CALL swarm_message_get(report, "pop_id",id)
            CALL swarm_message_get(report, "inp_num",input)
            CALL swarm_message_get(report, "result", func)
            this%f(id,:) = func
            if (this%num_input_rcvd .eq. this%num_ref_mol*offspringsize) then
            !do i=1,offspringsize
            !    print *,'this is f',this%f(i,:)
            !enddo
                !select type(mocmaes => this%population)
                ! type is (mo_cmaes_type)
                !    CALL mo_cmaes_run(mocmaes)
                !    this%offspring_list = mocmaes(1)%offspring_list
                !end select
                    do i=1,offspringsize
                        input = this%offspring_list(i)
                        CALL open_file(file_name="output", file_status="UNKNOWN",&
                        file_action="WRITE", file_position="APPEND",unit_number=iounit)
                        write(iounit,*) var_trans(this%x(input,:), &
                                    this%u_bound,this%l_bound),this%f(input,:)
                        CALL close_file(iounit)
                    enddo
                    select type(mocmaes => this%population)
                     type is (mo_cmaes_type)
                        CALL mo_cmaes_run(mocmaes)
                    this%offspring_list = mocmaes(1)%offspring_list
                    end select
                !---reset----------
                this%current_inp = 1
                this%num_input_rcvd = 0
                this%i_iteration = this%i_iteration + 1
            endif
        endif
        if (this%current_inp .le. offspringsize*this%num_ref_mol .and. &
            (report_status .eq. "wait_done" .or. report_status .eq. "initial_hello")) then
            input = this%offspring_list(this%current_inp)
            CALL swarm_message_add(cmd, "command", "function")
            CALL swarm_message_add(cmd, "inp_num", this%current_inp )
            CALL swarm_message_add(cmd, "pop_id", input )
            CALL swarm_message_add(cmd, "x",&
            var_trans(this%x(input,:),this%u_bound,this%l_bound))
            CALL swarm_message_add(cmd, "ref_func",1)
            this%current_inp = this%current_inp + 1
            print *,'send done'
        else
             CALL swarm_message_add(cmd, "command", "wait") 
        endif
    else
          CALL swarm_message_get(report, "worker_id", wid)
          CALL swarm_message_get(report, "status", report_status)
          if (trim(report_status) .eq. "done") then
            this%num_input_rcvd = this%num_input_rcvd + 1
            CALL swarm_message_get(report, "inp_num",input)
            CALL swarm_message_get(report, "pop_id",id)
            CALL index2grid(input,(/offspringsize, &
                            size(this%energies,dim=2),&
                            size(this%energies,dim=3)/),ind)
            CALL swarm_message_get(report, "result",func)
            CALL swarm_message_get(report, "converged", this%converged(id,ind(2),ind(3)))
            !if (this%converged(id,ind(2),ind(3)) .eq. 0) STOP "not converged"
            this%energies(id,ind(2),ind(3)) = func(1)
            this%rmsd(id,ind(2),ind(3))     = func(2)
            if (this%m .eq. 5) this%mcharge(id)  = this%mcharge(id) + func(5)

            if (this%num_input_rcvd .eq. offspringsize*this%num_ref_mol) then
                do i = 1,offspringsize
                    input = this%offspring_list(i)
                    CALL energy2function(this%num_ref_mol,this%energies(input,:,:), &
                                        this%ref_matrix,this%database,func_array)
                    this%f(input,1) = sum(func_array)
                    !this%f(input,2) = sum(this%rmsd(input,:,:))
                    ! this also is only first dim...so dimer of s22 yannick
                    if (this%database .eq. 1) then 
                        ! s22 is database 1
                        this%f(input,2) = sum(this%rmsd(input,:,1))
                    else if (this%database .eq. 5) then
                        this%f(input,1) = sum(func_array(1:30))
                        this%f(input,2) = sum(this%rmsd(input,:,:))
                        if (this%m .ne. 5) STOP "not enough objectives"
                        this%f(input,5) = sum(func_array(31:35))
                    else
                        this%f(input,2) = sum(this%rmsd(input,:,:))
                    endif
                    if (this%m .gt. 2) then
                      if (this%database .eq. 4) then
                        this%f(input,3) = func_array(7)
                      else
                        this%f(input,3) = func_array(2)
                      endif
                    endif
                    if (this%m .gt. 3) then
                      if (this%database .eq. 4) then
                        this%f(input,4) = this%rmsd(input,6,1)
                      else
                        this%f(input,4) = this%rmsd(input,1,1)
                      endif
                    endif
                    if (this%m .gt. 4) then
                        this%f(input,5) = this%mcharge(input)
                        print *,'end charge',this%mcharge(input)
                    endif
                    print *,'test',this%f(input,:)
                    !print *,'test',this%f(input,5)
                    !------------
                    if (any(this%f(input,:) .gt. 1e10)) this%f(input,:) = huge(1._dp)
                    CALL open_file(file_name="output", file_status="UNKNOWN",&
                        file_action="WRITE", file_position="APPEND",unit_number=iounit)
                    write(iounit,*) var_trans(this%x(input,:), &
                                    this%u_bound,this%l_bound),this%f(input,:),func_array,this%rmsd(input,:,1)
                    ! index 1 here is the first row...so the dimer of s22 ...
                    ! should be changed to be generic ?!
                    !                this%u_bound,this%l_bound),this%f(input,:),func_array,this%rmsd(input,:,:)
                    CALL close_file(iounit)
                enddo
                    select type(mocmaes => this%population)
                     type is (mo_cmaes_type)
                        CALL mo_cmaes_run(mocmaes)
                    this%offspring_list = mocmaes(1)%offspring_list
                    end select
                !---reset----------
                this%converged   = 1
                this%mcharge     = 0._dp
                this%current_inp = 1
                this%num_input_rcvd = 0
                this%i_iteration = this%i_iteration + 1
            endif
          endif
          if (this%current_inp .le.offspringsize*this%num_ref_mol .and. &
            (report_status .eq. "wait_done" .or. report_status .eq. "initial_hello")) then
            !-----------------------------------------------------------------
            CALL index2grid(this%current_inp,(/offspringsize,&
                            size(this%energies,dim=2),&
                            size(this%energies,dim=3)/),ind)
            input = this%offspring_list(ind(1))

            if (any(this%converged(input,:,:).eq.0)) then
                CALL swarm_message_add(cmd, "command", "skip")
            else
                CALL swarm_message_add(cmd, "command", "energy")
            endif
            CALL swarm_message_add(cmd, "inp_num", this%current_inp )

            print *,'iter',this%current_inp,ind
            CALL swarm_message_add(cmd, "pop_id", input )
            CALL gen_input_file(ind(2),ind(3),this%database,input_file)
            CALL swarm_message_add(cmd, "inp_file", input_file)
            ! change charge if necessary
            ! w27 1-15 are 0 then 16-21 are 1 then 22-29 are -1 and 30 is 0
            if (this%database .eq. 3) then
                select case(ind(2))
                    case(1:15,30)
                        charge = 0
                        spinmult = 1 
                    case(16:21)
                        charge = 1
                        spinmult = 1
                    case(22:29)
                        charge = -1
                        spinmult = 1 
                end select
            endif
            if (this%database .eq. 5) then
                select case(ind(2))
                    case(1:15,30)
                        charge = 0
                        spinmult = 1 
                    case(16:21)
                        charge = 1
                        spinmult = 1
                    case(22:29)
                        charge = -1
                        spinmult = 1
                    case(31)
                        charge = 0
                        spinmult = 2
                    case(32)
                        charge = 0
                        spinmult = 3
                    case(33)
                        charge = 1
                        spinmult = 4
                    case(34)
                        charge = 0
                        spinmult = 1
                    case(35)
                        charge = 0
                        spinmult = 3
                end select
            endif
            if (this%database .eq. 4) then
                select case(ind(2))
                case(1:2,4:7)
                        charge = 0
                    case(3)
                        charge = 1
                end select
                select case(ind(2))
                case(1)
                    spinmult = 2
                case(2)
                    spinmult = 3
                case(3)
                    spinmult = 4
                case(4)
                    spinmult = 1
                case(5)
                    spinmult = 3
                case(6:7)
                    spinmult = 1
                end select
            endif
            CALL swarm_message_add(cmd, "charge", charge)
            CALL swarm_message_add(cmd, "spinmult", spinmult)
            CALL swarm_message_add(cmd, "x",&
            var_trans(this%x(input,:),this%u_bound,this%l_bound))
            this%current_inp = this%current_inp + 1
          else
            CALL swarm_message_add(cmd, "command", "wait")     
          endif

      endif
      deallocate(func_array)
END SUBROUTINE do_mocma

! *****************************************************************************
!> \brief Helper routine for glbopt_master_steer(), updates stats, etc.
!> \param this ...
!> \param report ...
!> \author Ole Schuett
! *****************************************************************************
SUBROUTINE do_optim(this,report,cmd,method,should_stop)
    TYPE(paramopt_master_type), INTENT(INOUT)  :: this
    TYPE(swarm_message_type)                   :: report, cmd
    CHARACTER(len=default_string_length)       :: report_status
    CHARACTER(len=default_string_length)       :: input_file
    LOGICAL, INTENT(INOUT)                     :: should_stop
    INTEGER                                    :: wid
    INTEGER                                    :: input
    REAL(kind=dp)                              :: myresult
    INTEGER                                    :: charge=0
    INTEGER                                    :: spinmult=1
    REAL(kind=dp)                              :: penalty
    INTEGER                                    :: iounit 
    INTEGER                                    :: method 
    INTEGER                                    :: id=1 
    INTEGER,DIMENSION(3)                       :: ind
    if (this%benchmark) then
          CALL swarm_message_get(report, "worker_id", wid)
          CALL swarm_message_get(report, "status", report_status)
          if (trim(report_status) .eq. "done") then
            CALL swarm_message_get(report, "result", myresult)
            this%f = myresult
            this%x(1,:) = var_back_trans(this%x(1,:),this%u_bound,this%l_bound)
            select case(method)
                case(paramopt_do_powell)
                   CALL powell_optimize (this%n, this%x(1,:), &
                                  this%powell_param)
                   !CALL powell_optimize (this%n, this%x(1,:) ,this%powell_param%xl,&
                   !               this%powell_param%xu,&
                   !               this%powell_param)
                   if (this%powell_param%state .eq. -1) then
                    this%powell_param%state = 8
                    this%x(1,:) = var_back_trans(this%x(1,:),this%u_bound,this%l_bound)
                    CALL powell_optimize (this%n, this%x(1,:), &
                                  this%powell_param)
                    !CALL powell_optimize (this%n, this%x(1,:) ,this%l_bound,this%u_bound,&
                    !              this%powell_param)
                    this%x(1,:) = var_trans(this%x(1,:),this%u_bound,this%l_bound) 
                    IF ( this%iw > 0 ) THEN
                         WRITE(this%iw,'(" POWELL| Number of function evaluations",T71,I10)') this%powell_param%nf
                         WRITE(this%iw,'(" POWELL| Final value of function",T61,G20.10)') this%powell_param%fopt
                         WRITE(this%iw,'(" POWELL| Final x",T61,100G20.10)') this%x
                    END IF
                     STOP "powell done"
                    endif
                case(paramopt_do_cma)
                    CALL cmaes_run(this%cmaes_inp,this%population)
            end select 
            this%x(1,:) = var_trans(this%x(1,:),this%u_bound,this%l_bound) 

            CALL send_all(cmd,this%x,this%bench_func)
            !----------------------------------------------
            !CALL open_file(file_name="output", file_status="UNKNOWN",&
            !        file_action="WRITE", file_position="APPEND",unit_number=iounit)
            !        write(iounit,*) this%x(1,:),this%f(1,1)
            !CALL close_file(iounit)
          else if (report_status .eq. "wait_done" .or. report_status .eq. "initial_hello") then
            CALL send_all(cmd,this%x,this%bench_func)
          endif
    else
          CALL swarm_message_get(report, "worker_id", wid)
          CALL swarm_message_get(report, "status", report_status)
          if (trim(report_status) .eq. "done") then
            this%num_input_rcvd = this%num_input_rcvd + 1
            CALL swarm_message_get(report, "pop_id",id)
            CALL swarm_message_get(report, "inp_num",input)
            CALL index2grid(input,(/1,size(this%energies,dim=2),size(this%energies,dim=3)/),ind)
            CALL swarm_message_get(report, "result", this%energies(id,ind(2),ind(3)))
            CALL swarm_message_get(report, "converged", this%converged(id,ind(2),ind(3)))
!            CALL swarm_message_get(report, "converged", this%converged(id,ind(2),ind(3)))
!            if (this%converged(id,ind(2),ind(4)) .eq. 0) STOP "hm converged"
            !print *,'input',input,ind,this%num_input_rcvd
            if (this%num_input_rcvd .eq. this%num_ref_mol) then
!                CALL energy2function(this%energies(id,:,:),this%ref_matrix,this%database,myresult)
                !penalty = 1e3*this%penalty**2
                !this%f = myresult + penalty
                this%f = myresult 
                CALL open_file(file_name="output", file_status="UNKNOWN",&
                    file_action="WRITE", file_position="APPEND",unit_number=iounit)
                    write(iounit,*) this%x(1,:),this%f(1,1)
                CALL close_file(iounit)
                select case(method)
                    case(paramopt_do_powell)
                       this%x(1,:) = var_back_trans(this%x(1,:),this%u_bound,this%l_bound)
                       CALL powell_optimize (this%n, this%x(1,:), &
                                  this%powell_param)
                       !CALL powell_optimize (this%n, this%x(1,:) ,this%powell_param%xl,&
                       !           this%powell_param%xu,&
                       !           this%powell_param)
                       this%x(1,:) = var_trans(this%x(1,:),this%u_bound,this%l_bound) 
                       if (this%powell_param%state .eq. -1) then
                        this%powell_param%state = 8
                        this%x(1,:) = var_back_trans(this%x(1,:),this%u_bound,this%l_bound)
                        CALL powell_optimize (this%n, this%x(1,:), &
                                  this%powell_param)
                       ! CALL powell_optimize (this%n, this%x(1,:) ,this%powell_param%xl,&
                       !           this%powell_param%xu,&
                       !           this%powell_param)
                        this%x(1,:) = var_trans(this%x(1,:),this%u_bound,this%l_bound) 
                        IF ( this%iw > 0 ) THEN
                             WRITE(this%iw,'(" POWELL| Number of function evaluations",T71,I10)') this%powell_param%nf
                             WRITE(this%iw,'(" POWELL| Final value of function",T61,G20.10)') this%powell_param%fopt
                             WRITE(this%iw,'(" POWELL| Final x",T61,100G20.10)') this%x
                        END IF
                         STOP "powell done"
                        endif
                    case(paramopt_do_cma)
                        this%x(1,:) = var_back_trans(this%x(1,:),this%u_bound,this%l_bound)
                        CALL cmaes_run(this%cmaes_inp,this%population)
                        this%x(1,:) = var_trans(this%x(1,:),this%u_bound,this%l_bound) 
                end select 
                CALL penalty_fct(this%x(1,:),this%l_bound,this%u_bound,this%penalty)
                !---reset----------
                this%current_inp = 0
                this%num_input_rcvd = 0
                this%i_iteration = this%i_iteration + 1
            endif
          endif
          if (this%current_inp .le. this%num_ref_mol .and. &
            (report_status .eq. "wait_done" .or. report_status .eq. "initial_hello")) then
            CALL swarm_message_add(cmd, "command", "energy")
            CALL swarm_message_add(cmd, "pop_id", id )
            CALL swarm_message_add(cmd, "inp_num", this%current_inp )
            CALL index2grid(this%current_inp,(/1,size(this%energies,dim=2), &
                             size(this%energies,dim=3)/),ind)
            print *,'input',this%current_inp,ind
            CALL gen_input_file(ind(2),ind(3),this%database,input_file)
            CALL swarm_message_add(cmd, "inp_file", input_file)
            ! change charge if necessary
            ! w27 1-15 are 0 then 16-21 are 1 then 22-29 are -1 and 30 is 0
            if (this%database .eq. 3) then
                select case(ind(2))
                    case(1:15,30)
                        charge = 0
                    case(16:21)
                        charge = 1
                    case(22:29)
                        charge = -1
                end select
            endif
            CALL swarm_message_add(cmd, "charge", charge)
            CALL swarm_message_add(cmd, "spinmult", spinmult)
            CALL swarm_message_add(cmd, "x", this%x(1,:))
            this%current_inp = this%current_inp + 1
          else
            CALL swarm_message_add(cmd, "command", "wait")     
          endif

      endif
END SUBROUTINE do_optim
 
! *****************************************************************************
!> \brief Finalized the master of the global optimizer
!> \param this ...
!> \author Ole
! *****************************************************************************
 SUBROUTINE paramopt_master_finalize(this)
    TYPE(paramopt_master_type), INTENT(INOUT)  :: this

    TYPE(cp_logger_type), POINTER            :: logger

    NULLIFY(logger)

    SELECT CASE (this%method)
       CASE(paramopt_do_cma)
         !CALL minhop_finalize(this%minhop)
         !DEALLOCATE(this%minhop)
       CASE(2)
         !CALL mincrawl_finalize(this%mincrawl)
         !DEALLOCATE(this%mincrawl)
       CASE DEFAULT
        STOP "Unkown paramopt_method"
    END SELECT

    logger => cp_error_get_logger(this%error)

 END SUBROUTINE paramopt_master_finalize
! *****************************************************************************
!> \brief Finalizes master for Minima Hopping
!> \param this ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE energy2function(nref,energies,ref_matrix,database,func_array)
    REAL(KIND=dp),DIMENSION(:,:),INTENT(IN)         :: energies
    REAL(KIND=dp),DIMENSION(:,:),INTENT(IN)         :: ref_matrix
    INTEGER                                         :: database
    REAL(KIND=dp),DIMENSION(size(energies,dim=1))   :: func_array
    REAL(KIND=dp)                                   :: func
    REAL(KIND=dp)                                   :: kcalmol=6.27509468713739E+02_dp 
    INTEGER                                         :: i,k,n
    INTEGER                                         :: nref    
    REAL(KIND=dp),DIMENSION(:),ALLOCATABLE          :: ref_ene
    select case(database)
      case(1)
        print *,'nref',nref/3
        do i=1,nref/3
            func_array(i) = abs(energies(i,1)-(energies(i,2)+energies(i,3)))*kcalmol
            func_array(i) = abs(func_array(i) - ref_matrix(i,5))
            print *,'dimer',energies(i,1),'monomer',energies(i,2),'monomer',energies(i,3),'f',func_array(i),'ref',ref_matrix(i,5)
        enddo
        func = sum(abs(func_array))
      case(2)
        print *,'nref',nref
          do i=1,nref
            func_array(i) = abs(ref_matrix(i,2)*energies(1,1)-energies(i,1))*kcalmol
            func_array(i) = abs(func_array(i) - ref_matrix(i,5))
            print *,'cluster',energies(i,1),'monomer',energies(1,1),'x',ref_matrix(i,2), &
            ' * monomer',ref_matrix(i,2)*energies(1,1),'f',func_array(i)
        enddo
        func = sum(abs(func_array))
      case(3)
      allocate(ref_ene(4))
        ref_ene(1) = energies(1,1)  !H2O
        ref_ene(2) = energies(16,1) !H3Op
        ref_ene(3) = energies(22,1) !OHm
        ref_ene(4) = energies(11,1) !H2O8s4
        do i=1,nref
          func_array(i) = abs((ref_matrix(i,2)*energies(i,1) + dot_product(ref_matrix(i,3:6),ref_ene))*kcalmol-ref_matrix(i,7))
          print *,'cluster',energies(i,1),'x monomer',dot_product(ref_matrix(i,3:6),ref_ene),'f',func_array(i)
          if (ref_matrix(i,7) .ne. 0) func_array(i) = func_array(i)/ref_matrix(i,7)
        enddo
        func = sum(abs(func_array))
     case(4)
     allocate(ref_ene(4))
        ref_ene(1) = energies(2,1)  !O 
        ref_ene(2) = energies(4,1)  !H2
        ref_ene(3) = energies(5,1)  !O2
        ref_ene(4) = energies(6,1)  !H2O
        do i=1,nref
          print *,'energies',ref_matrix(i,2)*energies(i,1)*kcalmol,dot_product(ref_matrix(i,3:6),ref_ene)*kcalmol,ref_matrix(i,7)
          !print *,'test',abs((ref_matrix(i,2)*energies(i,1) + dot_product(ref_matrix(i,3:6),ref_ene))*kcalmol) ,ref_matrix(i,7)
          func_array(i) = abs(((ref_matrix(i,2)*energies(i,1) + dot_product(ref_matrix(i,3:6),ref_ene))*kcalmol)-ref_matrix(i,7))
          print *,'cluster ',i,energies(i,1),'x monomer',dot_product(ref_matrix(i,3:6),ref_ene),'f',func_array(i)
          if (ref_matrix(i,7) .ne. 0) func_array(i) = func_array(i)/ref_matrix(i,7)
        enddo
        func = sum(abs(func_array))
     case(5)
     allocate(ref_ene(7))
        ref_ene(1) = energies(1,1)  !H2O
        ref_ene(2) = energies(16,1) !H3Op
        ref_ene(3) = energies(22,1) !OHm
        ref_ene(4) = energies(11,1) !H2O8s4
        ref_ene(5) = energies(32,1) !O
        ref_ene(6) = energies(34,1) !H2
        ref_ene(7) = energies(35,1) !O2
        do i=1,nref
          func_array(i) = abs((ref_matrix(i,2)*energies(i,1) + dot_product(ref_matrix(i,3:6),ref_ene))*kcalmol-ref_matrix(i,7))
          print *,'cluster',energies(i,1),'x monomer',dot_product(ref_matrix(i,3:6),ref_ene),'f',func_array(i)
          if (ref_matrix(i,7) .ne. 0) func_array(i) = func_array(i)/ref_matrix(i,7)
        enddo
        func = sum(abs(func_array))
    end select
END SUBROUTINE energy2function
! *****************************************************************************
!> \brief Finalized the master of the global optimizer
!> \param this ...
!> \author Ole
! *****************************************************************************
 SUBROUTINE read_ref_data(database,ref_file,ref_matrix)
    CHARACTER(len=*)                         :: ref_file
    CHARACTER(len=default_string_length)     :: test
    INTEGER                                  :: n=0
    INTEGER                                  :: k
    INTEGER                                  :: error
    INTEGER                                  :: ref_unit
    INTEGER                                  :: database
    REAL(KIND=dp),ALLOCATABLE,DIMENSION(:,:) :: ref_matrix

       if (allocated(ref_matrix)) deallocate(ref_matrix)

       CALL open_file(file_name=ref_file,&
           file_action="READ", unit_number=ref_unit)

       DO; READ(ref_unit,*, iostat = error); IF (error == -1) EXIT; n = n + 1; END DO
       k = n - 1
     SELECT CASE(database)
        CASE(1,2)
           allocate(ref_matrix(k,5))
           REWIND(ref_unit)
           READ(ref_unit,*) test
           DO n=1,k
                READ(ref_unit,*)  ref_matrix(n,:)
           ENDDO
        CASE(3,4) 
           allocate(ref_matrix(k,7))
           REWIND(ref_unit)
           READ(ref_unit,*) test
           DO n=1,k
                READ(ref_unit,*)  ref_matrix(n,:)
           ENDDO
        CASE(5) 
           allocate(ref_matrix(k,10))
           REWIND(ref_unit)
           READ(ref_unit,*) test
           DO n=1,k
                READ(ref_unit,*)  ref_matrix(n,:)
           ENDDO
    END SELECT
           
      CALL close_file(ref_unit)

 END SUBROUTINE read_ref_data
! *****************************************************************************
!> \brief Initialize the variables for the powell optimizer
!> \param p_param ...
!> \param powell_section ...
!> \param error ...
!> \author Florian Schiffmann
! *****************************************************************************

   SUBROUTINE init_powell_var(this,p_param,powell_section,error)
      TYPE(paramopt_master_type)  :: this
      TYPE(opt_state_type)  :: p_param
      TYPE(cp_error_type), INTENT(inout)       :: error

      CHARACTER(len=*), PARAMETER :: routineN = 'init_powell_var',&
        routineP = moduleN//':'//routineN

      TYPE(section_vals_type), POINTER         :: powell_section
      INTEGER                                  :: n 
      p_param%state=0
      p_param%nf=0
      p_param%nvar=this%n
      p_param%iprint=3
      p_param%unit=this%iw
      CALL section_vals_val_get(powell_section,"ACCURACY",r_val=p_param%rhoend,error=error)
      CALL section_vals_val_get(powell_section,"STEP_SIZE",r_val=p_param%rhobeg,error=error)
      CALL cp_assert( this%step_size .eq. p_param%rhobeg ,cp_failure_level,cp_assertion_failed,&
                       routineP, "global step size different to powell step size",error)

      CALL section_vals_val_get(powell_section,"MAX_FUN",i_val=p_param%maxfun,error=error)

   END SUBROUTINE init_powell_var
! *****************************************************************************
!> \brief Finalizes master for Minima Hopping
!> \param this ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE gen_input_file(i,j,database,input_file)
    INTEGER                               :: i,j,database
    CHARACTER(len=default_string_length)  :: input_file
    CHARACTER(len=default_string_length)  :: tmp_str
 SELECT CASE(database)
  CASE(1)
    input_file="s22/config_"
    tmp_str = trim(str(i))
    if (len_trim(tmp_str) == 1 ) tmp_str = trim(str(0)) // trim(tmp_str)
    input_file = trim(input_file) // tmp_str
    SELECT CASE(j)
      CASE(1)
          input_file = trim(input_file) // ".xyz"
      CASE(2)
          input_file = trim(input_file) // "_monom1.xyz"
      CASE(3)
          input_file = trim(input_file) // "_monom2.xyz"
    END SELECT
   
  CASE(2)
    i = i - 1 
    input_file="water_clusters/"
    tmp_str = trim(str(i))
    if (len_trim(tmp_str) == 1 ) tmp_str = trim(str(0)) // trim(tmp_str)
    input_file = trim(input_file) // tmp_str
    input_file = trim(input_file) // "water.xyz"

  CASE(3)
    input_file="WATER27structures/water"
    tmp_str = trim(str(i))
    if (len_trim(tmp_str) == 1 ) tmp_str = trim(str(0)) // trim(tmp_str)
    input_file = trim(input_file) // tmp_str
    input_file = trim(input_file) // ".xyz"

  CASE(4)
    input_file="atomic_structures/water"
    tmp_str = trim(str(i))
    if (len_trim(tmp_str) == 1 ) tmp_str = trim(str(0)) // trim(tmp_str)
    input_file = trim(input_file) // tmp_str
    input_file = trim(input_file) // ".xyz"
END SELECT
END SUBROUTINE gen_input_file

! *****************************************************************************
!> \brief Finalizes master for Minima Hopping
!> \param this ...
!> \author Ole Schuett
! *****************************************************************************
SUBROUTINE grid_search(this)
    TYPE(paramopt_master_type), INTENT(INOUT)  :: this
    INTEGER                                    :: i,j,k,l,m
    INTEGER,DIMENSION(size(this%x,dim=2))            :: n
    INTEGER,DIMENSION(size(this%x,dim=2))            :: grid
    LOGICAL,DIMENSION(size(this%x,dim=2))            :: mask1,mask2
    n = int(1._dp/this%step_size)
    !n = int((this%u_bound-this%l_bound)/this%step_size)
    print *,'n',n
    if (this%i_iteration .gt. product(n)) STOP "grid sample is getting to big"
    grid(1) = modulo(this%i_iteration,n(1))  
    !this%x(1,1) = grid(1)*this%step_size+this%l_bound(1)
    this%x(1,1) = grid(1)*this%step_size
    do i=2,size(grid)
        mask1 = .FALSE. 
        mask2 = .FALSE. 
        mask1(1:i) = .TRUE.
        mask2(1:i-1) = .TRUE.
        !print *,'product',product(n,mask=mask1),product(n,mask=mask2)
        grid(i) = int(real(modulo(this%i_iteration,product(n,mask=mask1)))/real(product(n,mask=mask2)))
        !this%x(1,i) = grid(i)*this%step_size+this%l_bound(i)
        this%x(1,i) = grid(i)*this%step_size
    enddo
    print *,'grid',grid
    print *,'x',this%x(1,:)
END SUBROUTINE grid_search

! *****************************************************************************
!> \brief Finalizes master for Minima Hopping
!> \param this ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE index2grid(ind,ndim,grid_ind)
    INTEGER,INTENT(IN)                    :: ind   
    INTEGER,DIMENSION(:),INTENT(IN)       :: ndim !number of points in dimension
    INTEGER,DIMENSION(:),INTENT(OUT)      :: grid_ind !grid index
    LOGICAL,DIMENSION(size(ndim))         :: mask1,mask2
    INTEGER,DIMENSION(size(ndim))         :: mydim
    INTEGER                               :: i  
    INTEGER                               :: myind
    myind = ind - 1 ! so it starts from 0
    mydim = ndim 
    grid_ind(1) = modulo(myind,mydim(1))  
    do i=2,size(grid_ind)
        mask1 = .FALSE. 
        mask2 = .FALSE. 
        mask1(1:i) = .TRUE.
        mask2(1:i-1) = .TRUE.
        !print *,'product',product(n,mask=mask1),product(n,mask=mask2)
        grid_ind(i) = int(real(modulo(myind,product(mydim,mask=mask1)))/real(product(mydim,mask=mask2)))
    enddo
    grid_ind = grid_ind + 1
    !print *,'ndim',ndim
    !print *,'grid',grid_ind
 END SUBROUTINE index2grid

!   "Convert an integer to string."
 character(len=default_string_length) function str(k)
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
end function str


!   "penalty for x out of feasible region"
 subroutine penalty_fct(x,l_bound,u_bound,penalty)
    real(kind=dp),dimension(:), intent(inout) :: x
    real(kind=dp),dimension(:), intent(in) :: l_bound
    real(kind=dp),dimension(:), intent(in) :: u_bound
    real(kind=dp)                          :: penalty
    integer                                :: i
    penalty = 0._dp
    do i=1,size(x)
        if (x(i) .gt. u_bound(i)) then
           penalty = penalty + x(i) - u_bound(i)
            x(i) = u_bound(i) 
        endif
        if (x(i) .lt. l_bound(i)) then
           penalty = penalty - x(i) + l_bound(i) 
           x(i) = l_bound(i)
        endif
    enddo
    
end subroutine penalty_fct

! *****************************************************************************
!> \brief Central steering routine of Minima Hopping
!> \param this ...
!> \param report ...
!> \param cmd ...
!> \author Ole Schuett
! *****************************************************************************

    function var_trans(xtrans,u_bound,l_bound) result(x) 
        REAL(KIND=dp),DIMENSION(:),intent(in)              :: xtrans
        REAL(KIND=dp),DIMENSION(size(xtrans))              :: x 
        REAL(KIND=dp),DIMENSION(:)                         :: l_bound
        REAL(KIND=dp),DIMENSION(:)                         :: u_bound 
        if (any(xtrans .lt. 0._dp)) then
            print *,"var trans wrong xtrans lt",xtrans
            STOP
        endif
        if (any(xtrans .gt. 1._dp)) then
            print *, "var trans wrong xtrans gt",xtrans
            STOP
        endif
        x = xtrans*(u_bound-l_bound) + l_bound
    end function var_trans
    function var_back_trans(x,u_bound,l_bound) result(xtrans) 
        REAL(KIND=dp),DIMENSION(:),intent(in)              :: x
        REAL(KIND=dp),DIMENSION(size(x))                   :: xtrans 
        REAL(KIND=dp),DIMENSION(:)                         :: l_bound
        REAL(KIND=dp),DIMENSION(:)                         :: u_bound 
        if (any(x .lt. l_bound)) then
            print *,"var trans wrong lt",x
            STOP
        endif
        if (any(x .gt. u_bound)) then
            print *, "var trans wrong gt",x
            STOP
        endif
        xtrans = (x - l_bound)/(u_bound-l_bound)
    end function var_back_trans  

! *****************************************************************************
!> \brief Central steering routine of Minima Hopping
!> \param this ...
!> \param report ...
!> \param cmd ...
!> \author Ole Schuett
! *****************************************************************************

 SUBROUTINE create_rng(rng,input_seed)
    TYPE(rng_stream_type),POINTER            :: rng 
    REAL(kind=dp),DIMENSION(3,2)             :: seed
    INTEGER                                  :: time(3)
    INTEGER                                  :: input_seed
    TYPE(cp_error_type)                      :: error

    if (input_seed .gt. 0) then
        !deterministic seed
        seed(:,:) = &
            RESHAPE( (/input_seed*42.0_dp, &
                       input_seed*54.0_dp, &
                       input_seed*63.0_dp, &
                       input_seed*98.0_dp, &
                       input_seed*10.0_dp, &
                       input_seed*2.0_dp/),&
                     (/ 3, 2 /) )        
    else
        !random seed from itime
        call itime(time)
        seed(:,1) = REAL(time)
        seed(:,2) = (/ 1._dp, 2._dp, 3._dp /)
    endif 
    
    CALL create_rng_stream(rng,name='uniform',distribution_type=UNIFORM,seed=seed,error=error)    
 END SUBROUTINE create_rng
! *****************************************************************************
!> \brief Central steering routine of Minima Hopping
!> \param this ...
!> \param report ...
!> \param cmd ...
!> \author Ole Schuett
! *****************************************************************************

 SUBROUTINE send_all(cmd,x,bench_func,id)
    TYPE(swarm_message_type)                   :: cmd
    REAL(kind=dp),DIMENSION(:,:),POINTER       :: x
    INTEGER                                    :: bench_func
    INTEGER,OPTIONAL                           :: id

    if (present(id)) then
    CALL swarm_message_add(cmd, "command", "function")
    CALL swarm_message_add(cmd, "pop_id", id)
    CALL swarm_message_add(cmd, "x", x(id,:))
    CALL swarm_message_add(cmd, "ref_func",bench_func) 
    else
    CALL swarm_message_add(cmd, "command", "function")
    CALL swarm_message_add(cmd, "x", x(1,:))
    CALL swarm_message_add(cmd, "ref_func",bench_func) 
    endif 
 END SUBROUTINE send_all

END MODULE paramopt_master

