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
                                             pop_init
  USE cp_files,                        ONLY: open_file,&
                                             close_file
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_units,                        ONLY: cp_unit_from_cp2k
  USE input_constants,                 ONLY: dump_xmol,&
                                             paramopt_do_cma,&
                                             paramopt_do_powell,&
                                             paramopt_do_grid
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
  USE powell,                          ONLY: powell_optimize, &
                                             opt_state_type
                                             
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
   INTEGER                                             :: num_ref_mol
   INTEGER                                             :: num_input_rcvd
   INTEGER                                             :: i_iteration = 0
   TYPE(section_vals_type), POINTER                    :: paramopt_section => Null()
   TYPE(cp_error_type)                                 :: error
   !--------------------------
   REAL(KIND=dp),ALLOCATABLE,DIMENSION(:,:)            :: energies
   REAL(KIND=dp),ALLOCATABLE,DIMENSION(:,:)            :: ref_matrix
   INTEGER                                             :: database=1
   REAL(KIND=dp),DIMENSION(:),POINTER                  :: x
   REAL(KIND=dp),DIMENSION(:),POINTER                  :: f
   REAL(KIND=dp),DIMENSION(:),POINTER                  :: u_bound
   REAL(KIND=dp),DIMENSION(:),POINTER                  :: l_bound
   INTEGER                                             :: n,m
   REAL(KIND=dp)                                       :: step_size
   LOGICAL                                             :: randinput
   !--------------------------
   TYPE(cmaes_input)                                   :: cmaes_inp
   TYPE(cmaes_type),ALLOCATABLE,DIMENSION(:)           :: population
   !--------------------------
   TYPE(opt_state_type)                                :: powell_param
   !-------------------------
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
    TYPE(section_vals_type), POINTER                    :: powell_section
    CHARACTER(LEN=*), INTENT(IN)             :: input_path
    INTEGER, INTENT(IN)                      :: n_walkers, iw
    TYPE(cp_error_type), INTENT(inout)       :: error
    CHARACTER(len=*), PARAMETER :: routineN = 'paramopt_master_steer', &
      routineP = moduleN//':'//routineN
    TYPE(cp_logger_type), POINTER            :: logger
    INTEGER                                  :: n,m
    INTEGER                                  :: i
    INTEGER                                  :: j
    LOGICAL                                  :: explicit_key
    REAL(KIND=dp),DIMENSION(:),POINTER       :: init_x
    NULLIFY(logger)

    this%iw = iw
    this%error = error
    !------------------------------------
    this%paramopt_section => section_vals_get_subs_vals(root_section, "SWARM%PARAM_OPT", error=error)
    powell_section => section_vals_get_subs_vals(root_section, "ATOM%POWELL", error=error)
    CALL section_vals_retain(this%paramopt_section,error=error)
    CALL section_vals_val_get(this%paramopt_section,"NUM_VAR", i_val=this%n, error=error)
    !------------------------------------
    n = this%n
    this%bench_func = 1
    allocate(this%x(n))
    !------------------------------------
    CALL section_vals_val_get(this%paramopt_section,"INITIAL_X", explicit=explicit_key, error=error)
    IF(explicit_key)THEN
        CALL section_vals_val_get(this%paramopt_section,"INITIAL_X", r_vals=init_x, error=error)
        this%x = init_x
        this%randinput = .FALSE.
        CALL cp_assert( size(init_x) .eq. n ,cp_failure_level,cp_assertion_failed,&
                       routineP, "init_x does not match num_var",error)
    ELSE
        this%randinput = .TRUE.
    ENDIF
    CALL section_vals_val_get(this%paramopt_section,"UPPER_BOUNDS", r_vals=this%u_bound, error=error)
    CALL section_vals_val_get(this%paramopt_section,"LOWER_BOUNDS", r_vals=this%l_bound, error=error)
    CALL section_vals_val_get(this%paramopt_section,"METHOD", i_val=this%method, error=error)
    CALL section_vals_val_get(this%paramopt_section,"STEP_SIZE", r_val=this%step_size, error=error)
    !------------------------------------
    this%m = 1
    m = this%m
    this%database = 3 ! 1 for the s22 and 2 for water_cluster
    if (this%database == 1) then
        this%num_ref_mol = 3*22
        CALL read_ref_data(this%database,"s22/reference_data",this%ref_matrix) 
        allocate(this%energies(this%num_ref_mol/3,3))
    else if (this%database == 2) then
        this%num_ref_mol = 39
        CALL read_ref_data(this%database,"water_clusters/reference_data",this%ref_matrix) 
        allocate(this%energies(this%num_ref_mol,1))
    else if (this%database == 3) then
        this%num_ref_mol = 30 
        CALL read_ref_data(this%database,"WATER27structures/reference_data",this%ref_matrix) 
        allocate(this%energies(this%num_ref_mol,1))
    endif
    this%num_input_rcvd = 0
    this%current_inp = 0
    allocate(this%f(m))
    logger => cp_error_get_logger(error)
    ! ref matrix 
    ! Nr.  Dimer monomer1 monomer2 ref_energy
    !do, i=1,size(ref_matrix,dim=1)
    !     print *, ( ref_matrix(i,j), j=1,size(ref_matrix,dim=2) )
    !enddo  
    
    !this%method = paramopt_do_cma
    !this%method = paramopt_do_powell
    SELECT CASE (this%method)
       CASE(paramopt_do_cma)
            this%cmaes_inp%step_size = this%step_size
            this%cmaes_inp%n = n
            this%cmaes_inp%randinput = this%randinput
            this%cmaes_inp%input_seed = 0
            this%cmaes_inp%m = m
            this%cmaes_inp%x => this%x
            print *,'check if ok',this%cmaes_inp%x
            this%cmaes_inp%f => this%f
            this%cmaes_inp%l_bound => this%l_bound
            this%cmaes_inp%u_bound => this%u_bound
            CALL pop_init(this%population,this%cmaes_inp)
       CASE(paramopt_do_powell)
            CALL init_powell_var(this,this%powell_param,powell_section,error)
            CALL powell_optimize (this%n, this%x , this%powell_param)
       CASE(paramopt_do_grid)
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


    !CALL progress_report(this, report)
    
    SELECT CASE (this%method)
        CASE(paramopt_do_cma)
            CALL do_cma(this,report,cmd,should_stop)
        CASE(paramopt_do_powell)
            CALL do_powell(this,report,cmd,should_stop)
        CASE(paramopt_do_grid)
            CALL do_grid_search(this,report,cmd,should_stop)
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

          CALL swarm_message_get(report, "worker_id", wid)
          CALL swarm_message_get(report, "status", report_status)
          if (trim(report_status) .eq. "done") then
            this%num_input_rcvd = this%num_input_rcvd + 1
            CALL swarm_message_get(report, "inp_num",input)
            CALL index2colrow(input,size(this%energies,dim=2),i,j) 
            CALL swarm_message_get(report, "result", this%energies(i,j))
            if (this%num_input_rcvd .eq. this%num_ref_mol) then
                CALL energy2function(this%energies,this%ref_matrix,this%database,myresult)
                this%f = myresult
                CALL open_file(file_name="output", file_status="UNKNOWN",&
                    file_action="WRITE", file_position="APPEND",unit_number=iounit)
                    write(iounit,*) this%x,myresult
                CALL close_file(iounit)
                CALL grid_search(this) 
                !---reset----------
                this%current_inp = 0
                this%num_input_rcvd = 0
                this%i_iteration = this%i_iteration + 1
            endif
          endif
          if (this%current_inp .lt. this%num_ref_mol .and. &
            (report_status .eq. "wait_done" .or. report_status .eq. "initial_hello")) then
            CALL swarm_message_add(cmd, "command", "energy")
            CALL swarm_message_add(cmd, "inp_num", this%current_inp )
            CALL index2colrow(this%current_inp,size(this%energies,dim=2),i,j)
            CALL gen_input_file(i,j,this%database,input_file)
            CALL swarm_message_add(cmd, "inp_file", input_file)
            ! change charge if necessary
            ! w27 1-15 are 0 then 16-21 are 1 then 22-29 are -1 and 30 is 0
            if (this%database .eq. 3) then
                select case(i)
                    case(1:15,30)
                        charge = 0
                    case(16:21)
                        charge = 1
                    case(22:29)
                        charge = -1
                end select
            endif
            CALL swarm_message_add(cmd, "charge", charge)
            CALL swarm_message_add(cmd, "x", this%x)
            CALL swarm_message_add(cmd, "ref_func", this%bench_func)
            this%current_inp = this%current_inp + 1
          else
            CALL swarm_message_add(cmd, "command", "wait")     
          endif

END SUBROUTINE do_grid_search

! *****************************************************************************
!> \brief Helper routine for glbopt_master_steer(), updates stats, etc.
!> \param this ...
!> \param report ...
!> \author Ole Schuett
! *****************************************************************************
SUBROUTINE do_powell(this,report,cmd,should_stop)
    TYPE(paramopt_master_type), INTENT(INOUT)  :: this
    TYPE(swarm_message_type)                   :: report, cmd
    CHARACTER(len=default_string_length)       :: report_status
    CHARACTER(len=default_string_length)       :: input_file
    LOGICAL, INTENT(INOUT)                     :: should_stop
    INTEGER                                    :: wid
    INTEGER                                    :: input
    REAL(kind=dp)                              :: myresult
    INTEGER                                    :: charge=0
    REAL(kind=dp)                              :: penalty
    INTEGER                                    :: i,j 
    INTEGER                                    :: iounit 

          CALL swarm_message_get(report, "worker_id", wid)
          CALL swarm_message_get(report, "status", report_status)
          if (trim(report_status) .eq. "done") then
            this%num_input_rcvd = this%num_input_rcvd + 1
            CALL swarm_message_get(report, "inp_num",input)
            CALL index2colrow(input,size(this%energies,dim=2),i,j) 
            CALL swarm_message_get(report, "result", this%energies(i,j))
            if (this%num_input_rcvd .eq. this%num_ref_mol) then
                CALL energy2function(this%energies,this%ref_matrix,this%database,myresult)
                penalty = 1e3*penalty_fct(this%x,this%l_bound,this%u_bound)**2
                this%powell_param%f = myresult + penalty
                CALL open_file(file_name="output", file_status="UNKNOWN",&
                    file_action="WRITE", file_position="APPEND",unit_number=iounit)
                    write(iounit,*) this%x,this%powell_param%f,penalty
                CALL close_file(iounit)
                CALL powell_optimize (this%n, this%x , this%powell_param)
                if (this%powell_param%state .eq. -1) then
                     this%powell_param%state = 8
                     CALL powell_optimize (this%n, this%x , this%powell_param)
                     IF ( this%iw > 0 ) THEN
                     WRITE(this%iw,'(" POWELL| Number of function evaluations",T71,I10)') this%powell_param%nf
                     WRITE(this%iw,'(" POWELL| Final value of function",T61,G20.10)') this%powell_param%fopt
                   END IF
                    STOP "powell done"
                endif
                !---reset----------
                this%current_inp = 0
                this%num_input_rcvd = 0
                this%i_iteration = this%i_iteration + 1
            endif
          endif
          if (this%current_inp .lt. this%num_ref_mol .and. &
            (report_status .eq. "wait_done" .or. report_status .eq. "initial_hello")) then
            CALL swarm_message_add(cmd, "command", "energy")
            CALL swarm_message_add(cmd, "inp_num", this%current_inp )
            CALL index2colrow(this%current_inp,size(this%energies,dim=2),i,j)
            CALL gen_input_file(i,j,this%database,input_file)
            CALL swarm_message_add(cmd, "inp_file", input_file)
            CALL swarm_message_add(cmd, "charge", charge)
            CALL swarm_message_add(cmd, "x", this%x)
            CALL swarm_message_add(cmd, "ref_func", this%bench_func)
            this%current_inp = this%current_inp + 1
          else
            CALL swarm_message_add(cmd, "command", "wait")     
          endif
END SUBROUTINE do_powell
 
! *****************************************************************************
!> \brief Helper routine for glbopt_master_steer(), updates stats, etc.
!> \param this ...
!> \param report ...
!> \author Ole Schuett
! *****************************************************************************
SUBROUTINE do_cma(this,report,cmd,should_stop)
    TYPE(paramopt_master_type), INTENT(INOUT)  :: this
    TYPE(swarm_message_type)                   :: report, cmd
    CHARACTER(len=default_string_length)       :: report_status
    CHARACTER(len=default_string_length)       :: input_file
    LOGICAL, INTENT(INOUT)                     :: should_stop
    INTEGER                                    :: wid
    INTEGER                                    :: input
    REAL(kind=dp)                              :: myresult
    INTEGER                                    :: charge=0
    REAL(kind=dp)                              :: penalty
    INTEGER                                    :: i,j 
    INTEGER                                    :: iounit 

          CALL swarm_message_get(report, "worker_id", wid)
          CALL swarm_message_get(report, "status", report_status)
          if (trim(report_status) .eq. "done") then
            this%num_input_rcvd = this%num_input_rcvd + 1
            CALL swarm_message_get(report, "inp_num",input)
            CALL index2colrow(input,size(this%energies,dim=2),i,j) 
            CALL swarm_message_get(report, "result", this%energies(i,j))
            !print *,'input',input,'matrix i,j',i ,j
            !print *,'energy',this%energies(i,j)
            if (this%num_input_rcvd .eq. this%num_ref_mol) then
                CALL energy2function(this%energies,this%ref_matrix,this%database,myresult)
                !boundary check
                penalty = 1e3*penalty_fct(this%x,this%l_bound,this%u_bound)**2
                myresult = myresult + penalty
                this%cmaes_inp%f(1) = myresult
                CALL open_file(file_name="output", file_status="UNKNOWN",&
                    file_action="WRITE", file_position="APPEND",unit_number=iounit)
                    write(iounit,*) this%x,myresult,penalty
                CALL close_file(iounit)
                CALL cmaes_run(this%cmaes_inp,this%population)
                !---reset----------
                this%current_inp = 0
                this%num_input_rcvd = 0
                this%i_iteration = this%i_iteration + 1
            endif
          endif
          if (this%current_inp .lt. this%num_ref_mol .and. &
            (report_status .eq. "wait_done" .or. report_status .eq. "initial_hello")) then
            CALL swarm_message_add(cmd, "command", "energy")
            CALL swarm_message_add(cmd, "inp_num", this%current_inp )
            CALL index2colrow(this%current_inp,size(this%energies,dim=2),i,j)
            CALL gen_input_file(i,j,this%database,input_file)
            CALL swarm_message_add(cmd, "inp_file", input_file)
            CALL swarm_message_add(cmd, "charge", charge)
            CALL swarm_message_add(cmd, "x", this%x)
            CALL swarm_message_add(cmd, "ref_func", this%bench_func)
            this%current_inp = this%current_inp + 1
          else
            CALL swarm_message_add(cmd, "command", "wait")     
          endif
END SUBROUTINE do_cma 

! *****************************************************************************
!> \brief Helper routine for glbopt_master_steer(), updates stats, etc.
!> \param this ...
!> \param report ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE progress_report(this, report)
    TYPE(paramopt_master_type), INTENT(INOUT)  :: this
    TYPE(swarm_message_type)                 :: report

    CHARACTER(len=default_string_length)     :: status
    INTEGER                                  :: gopt_steps, md_steps, &
                                                report_worker_id
    REAL(KIND=dp)                            :: report_Epot

    this%i_iteration = this%i_iteration + 1

    CALL swarm_message_get(report, "worker_id", report_worker_id)
    CALL swarm_message_get(report, "status", status)

    IF(TRIM(status) == "ok") THEN
       CALL swarm_message_get(report, "Epot", report_Epot)
       CALL swarm_message_get(report, "md_steps", md_steps)
       CALL swarm_message_get(report, "gopt_steps", gopt_steps)
       this%count_reports    = this%count_reports + 1


       IF(this%iw>0) THEN
          WRITE (this%iw,'(A,46X,I8)')&
             " GLBOPT| Reporting worker ", report_worker_id
       END IF
    END IF
 END SUBROUTINE progress_report

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
    !CALL cp_print_key_finished_output(this%progress_traj_unit, logger, &
    !  this%glbopt_section, "PROGRESS_TRAJECTORY", error=this%error)

    !CALL section_vals_release(this%glbopt_section, error=this%error)

 END SUBROUTINE paramopt_master_finalize
! *****************************************************************************
!> \brief Finalizes master for Minima Hopping
!> \param this ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE energy2function(energies,ref_matrix,database,func)
    REAL(KIND=dp),DIMENSION(:,:),INTENT(IN)         :: energies
    REAL(KIND=dp),DIMENSION(:,:),INTENT(IN)         :: ref_matrix
    INTEGER                                         :: database
    REAL(KIND=dp),DIMENSION(size(energies,dim=1))   :: func_array
    REAL(KIND=dp)                                   :: func
    REAL(KIND=dp)                                   :: kcalmol=6.27509468713739E+02_dp 
    INTEGER                                         :: i,k,n
    REAL(KIND=dp),DIMENSION(4)                      :: ref_ene
    n=size(energies,dim=1)
    select case(database)
      case(1)
        do i=1,size(energies,dim=1)
            func_array(i) = abs(energies(i,1)-(energies(i,2)+energies(i,3)))*kcalmol
            print *,'dimer',energies(i,1),'monomer',energies(i,2),'monomer',energies(i,3),'f',func_array(i),'ref',ref_matrix(i,5)
            func_array(i) = abs(func_array(i) - ref_matrix(i,5))
        enddo
        func = sum(abs(func_array))
      case(2)
        do i=1,size(energies,dim=1)-1
            k = i + 1
            func_array(i) = ref_matrix(k,2)*energies(1,1)-energies(k,1)
            print *,'cluster',energies(k,1),'monomer',energies(1,1),'x * monomer',ref_matrix(k,2)*energies(1,1),'f',func_array(i)
        enddo
        func = sum(abs(func_array))
      case(3)
        ref_ene(1) = energies(1,1) !H2O
        ref_ene(2) = energies(16,1) !H3Op
        ref_ene(3) = energies(22,1) !OHm
        ref_ene(4) = energies(11,1) !H2O8s4
        do i=1,size(energies,dim=1)
          func_array(i) = (ref_matrix(i,2)*energies(i,1) + dot_product(ref_matrix(i,3:6),ref_ene))*kcalmol-ref_matrix(i,7)
          print *,'cluster',energies(i,1),'x monomer',dot_product(ref_matrix(i,3:6),ref_ene),'f',func_array(i)
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
        CASE(3) 
           allocate(ref_matrix(k,7))
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
    INTEGER,DIMENSION(size(this%x))            :: n
    INTEGER,DIMENSION(size(this%x))            :: grid
    LOGICAL,DIMENSION(size(this%x))            :: mask1,mask2
   
    n = int((this%u_bound-this%l_bound)/this%step_size)
    print *,'n',n
    if (this%i_iteration .gt. product(n)) STOP "grid sample is getting to big"
    grid(1) = modulo(this%i_iteration,n(1))  
    this%x(1) = grid(1)*this%step_size+this%l_bound(1)
    do i=2,size(grid)
        mask1 = .FALSE. 
        mask2 = .FALSE. 
        mask1(1:i) = .TRUE.
        mask2(1:i-1) = .TRUE.
        !print *,'product',product(n,mask=mask1),product(n,mask=mask2)
        grid(i) = int(real(modulo(this%i_iteration,product(n,mask=mask1)))/real(product(n,mask=mask2)))
        this%x(i) = grid(i)*this%step_size+this%l_bound(i)
    enddo
    print *,'grid',grid
    !    j = int(real(modulo(this%i_iteration,n(1)*n(2)))/real(n(2)))
    !    k = int(real(modulo(this%i_iteration,n(1)*n(2)*n(3)))/real(n(2)*n(3)))
    !    l = int(real(modulo(this%i_iteration,n(1)*n(2)*n(3)*n(4)))/real(n(2)*n(3)*n(4)))
    !    this%x(1) = i*this%step_size+this%l_bound(1)
    !    this%x(2) = j*this%step_size+this%l_bound(2)
    !    this%x(3) = k*this%step_size+this%l_bound(3)
    !    this%x(4) = l*this%step_size+this%l_bound(4)
END SUBROUTINE grid_search

! *****************************************************************************
!> \brief Finalizes master for Minima Hopping
!> \param this ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE index2colrow(input,indim,i,j)
    INTEGER,INTENT(IN)                    :: input,indim
    INTEGER,INTENT(OUT)                   :: j,i
    j = modulo(input,indim) + 1
    i = ceiling(real(input + 1)/real(indim)) 
 END SUBROUTINE index2colrow

!   "Convert an integer to string."
 character(len=default_string_length) function str(k)
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
end function str


!   "penalty for x out of feasible region"
 real(KIND=dp) function penalty_fct(x,l_bound,u_bound)
    real(kind=dp),dimension(:), intent(in) :: x
    real(kind=dp),dimension(:), intent(in) :: l_bound
    real(kind=dp),dimension(:), intent(in) :: u_bound
    integer                                :: i
    penalty_fct = 0._dp
    do i=1,size(x)
        if (x(i) .gt. u_bound(i)) then
           penalty_fct = penalty_fct + x(i) - u_bound(i) 
        endif
        if (x(i) .lt. l_bound(i)) then
           penalty_fct = penalty_fct - x(i) + l_bound(i) 
        endif
    enddo
    
end function penalty_fct



END MODULE paramopt_master

