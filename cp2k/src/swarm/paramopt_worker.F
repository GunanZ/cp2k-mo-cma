!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Worker routines used by global optimization schemes
!> \author Ole Schuett
! *****************************************************************************
MODULE paramopt_worker
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type,&
                                             pack_subsys_particles,&
                                             unpack_subsys_particles
  USE f77_interface,                   ONLY: calc_energy_force,&
                                             create_force_env,&
                                             destroy_force_env,&
                                             f_env_add_defaults,&
                                             f_env_type,&
                                             get_energy,&
                                             get_force,&
                                             get_natom, &
                                             get_pos, &
                                             calc_energy, &
                                             f_env_get_from_id
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE geo_opt,                         ONLY: cp_geo_opt
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: dump_xmol,&
                                             paramopt_do_cma,&
                                             paramopt_do_powell,&
                                             paramopt_do_grid
  USE input_section_types,             ONLY: section_type,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set,& 
                                             section_vals_remove_values                                             
  USE kinds,                           ONLY: default_string_length,&
                                             default_path_length,&
                                             dp
  USE md_run,                          ONLY: qs_mol_dyn
  USE mdctrl_types,                    ONLY: glbopt_mdctrl_data_type,&
                                             mdctrl_type
  USE parallel_rng_types,              ONLY: reset_to_next_rng_substream
  USE physcon,                         ONLY: angstrom,&
                                             kelvin
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE swarm_message,                   ONLY: swarm_message_add,&
                                             swarm_message_get,&
                                             swarm_message_type
#include "../common/cp_common_uses.f90"

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'paramopt_worker'

 PUBLIC :: paramopt_worker_init, paramopt_worker_finalize
 PUBLIC :: paramopt_worker_execute
 PUBLIC :: paramopt_worker_type

  TYPE paramopt_worker_type
   PRIVATE
   INTEGER                                  :: id
   INTEGER                                  :: iw
   INTEGER                                  :: group
   REAL(KIND=dp),DIMENSION(:),POINTER       :: x
   INTEGER                                  :: n
   TYPE(cp_error_type)                      :: error
   TYPE(cp_para_env_type), POINTER          :: para_env
   INTEGER                                  :: f_env_id
   INTEGER                                  :: method
   INTEGER                                  :: inp_num
   INTEGER                                  :: charge 
   TYPE(f_env_type), POINTER                :: f_env
   TYPE(force_env_type), POINTER            :: force_env
   TYPE(qs_environment_type), POINTER       :: qs_env
   TYPE(section_vals_type), POINTER         :: subsys_section
   TYPE(cp_subsys_type), POINTER            :: subsys
   TYPE(section_vals_type), POINTER         :: root_section
   TYPE(global_environment_type), POINTER   :: globenv
   TYPE(section_type), POINTER              :: input_declaration
   CHARACTER(len=default_path_length)       :: input_path
   CHARACTER(len=default_string_length)     :: inp_file
   INTEGER                                  :: ref_func = 1
   REAL(KIND=dp),DIMENSION(:),POINTER       :: D3 => Null()
 END TYPE paramopt_worker_type


 CONTAINS


! *****************************************************************************
!> \brief Initializes worker for global optimization
!> \param worker ...
!> \param input_declaration ...
!> \param para_env ...
!> \param root_section ...
!> \param input_path ...
!> \param worker_id ...
!> \param iw ...
!> \param error ...
!> \author Ole Schuett
! *****************************************************************************
   SUBROUTINE paramopt_worker_init(worker, input_declaration, para_env, root_section,&
                 input_path, worker_id, iw, error)
    TYPE(paramopt_worker_type), INTENT(INOUT)  :: worker
    TYPE(section_type), POINTER              :: input_declaration
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(section_vals_type), POINTER         :: swarm_section => Null()
    TYPE(section_vals_type), POINTER         :: paramopt_section => Null()
    CHARACTER(LEN=*), INTENT(IN)             :: input_path
    INTEGER, INTENT(in)                      :: worker_id, iw
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER                                  :: i
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: dist_in_angstrom
    TYPE(section_vals_type), POINTER         :: glbopt_section

        failure = .FALSE.

        worker%error = error
        worker%root_section => root_section
        worker%input_path = input_path
        worker%input_declaration => input_declaration
        worker%id = worker_id
        worker%iw = iw
        worker%para_env => para_env
        worker%group = worker%para_env%group
        swarm_section => section_vals_get_subs_vals(root_section, "SWARM", error=error)
        paramopt_section => section_vals_get_subs_vals(swarm_section, "PARAM_OPT", error=error)
        CALL section_vals_val_get(paramopt_section,"NUM_VAR", i_val=worker%n, error=error)
        CALL section_vals_val_get(paramopt_section,"METHOD", i_val=worker%method, error=error) 
        IF (worker%iw>0) WRITE(worker%iw,*) "paramopt_worker_init called."
        IF (worker%iw>0) WRITE(worker%iw,*) "paramopt_worker_init done."
        allocate(worker%D3(3))
  END SUBROUTINE paramopt_worker_init


! *****************************************************************************
!> \brief Central execute routine of global optimization worker
!> \param worker ...
!> \param cmd ...
!> \param report ...
!> \author Ole Schuett
! *****************************************************************************
   SUBROUTINE paramopt_worker_execute(worker, cmd, report)
    TYPE(paramopt_worker_type), INTENT(INOUT)  :: worker
    TYPE(swarm_message_type), INTENT(IN)     :: cmd
    TYPE(swarm_message_type), INTENT(INOUT)  :: report

    CHARACTER(len=default_string_length)     :: command
    REAL(kind=dp)                            :: test
   REAL(kind=dp)                            :: y


    SELECT CASE (worker%method)
       CASE(paramopt_do_cma,paramopt_do_powell,paramopt_do_grid)
            CALL swarm_message_get(cmd, "command", command)
            IF(TRIM(command) == "function") THEN
                nullify(worker%x)
                CALL swarm_message_get(cmd, "x", worker%x)
                CALL swarm_message_get(cmd, "ref_func", worker%ref_func)
                !CALL run_mdgopt(worker, cmd, report)
                !------------do calc here------------
                CALL myfun(worker%x,y,worker%ref_func)
                !------------------------------------ 
                CALL swarm_message_add(report, "status", "done") 
                CALL swarm_message_add(report, "result", y)
            ELSE IF (TRIM(command) == "energy") THEN
               nullify(worker%x)
               CALL swarm_message_get(cmd, "inp_num", worker%inp_num)
               CALL swarm_message_get(cmd, "inp_file", worker%inp_file)
               CALL swarm_message_get(cmd, "charge", worker%charge)
               CALL swarm_message_get(cmd, "x", worker%x)
               CALL swarm_message_get(cmd, "ref_func", worker%ref_func)
               !-----------------------------------------------
               CALL paramopt_energy(worker,cmd,command,report,y)
               !-----------------------------------------------
               CALL swarm_message_add(report, "status", "done")
               CALL swarm_message_add(report, "inp_num", worker%inp_num) 
               CALL swarm_message_add(report, "result", y)
            ELSE
               print *,"Worker: received unkown command:"
               STOP 
            END IF
       CASE DEFAULT
        STOP "Unkown glbopt_method"
    END SELECT
   END SUBROUTINE paramopt_worker_execute

! *****************************************************************************
!> \brief Finalizes worker for global optimization
!> \param worker ...
!> \author Ole Schuett
! *****************************************************************************
   SUBROUTINE paramopt_energy(this, cmd, command, report, energy)
    TYPE(paramopt_worker_type), &
      INTENT(INOUT)                          :: this
    CHARACTER(len=*), PARAMETER :: routineN = 'paramopt_geo_energy', &
      routineP = moduleN//':'//routineN
    CHARACTER(len=default_string_length)     :: command
    TYPE(swarm_message_type), INTENT(IN)     :: cmd
    TYPE(swarm_message_type), INTENT(INOUT)  :: report
    TYPE(section_vals_type), POINTER         :: work_section => Null()
    CHARACTER(len=default_string_length)     :: output_path, input_file
    LOGICAL                                  :: failure
    LOGICAL                                  :: converged
    INTEGER                                  :: ierr
    INTEGER                                  :: iter,i_inp,charge,i
    REAL(KIND=dp),DIMENSION(:),POINTER       :: x => Null()
    REAL(KIND=dp)                            :: energy
    REAL(KIND=dp)                            :: hb_para
    REAL(KIND=dp),DIMENSION(:),POINTER       :: dftb3_param => Null()
    REAL(KIND=dp),DIMENSION(:),POINTER       :: d3_ptr => Null()
    !CALL swarm_message_get(cmd,"input",input_file)
    !CALL swarm_message_get(cmd,"coord",coord_file)
    !CALL swarm_message_get(cmd,"num_inp",i_inp)
    !CALL swarm_message_get(cmd,"iter",iter)
    !CALL swarm_message_get(cmd,"x",this%x)

    !----for D3 optimization---------------------------
    !allocate(d3_ptr(3)) !gets deallocated automatically ??? at least it seems
    !this%D3 = (/1._dp,this%x(1:2)/)
    !!this%D3 = this%x(1:3)
    !d3_ptr = this%D3
    allocate(dftb3_param(2)) 
    dftb3_param = this%x(2:size(this%x))
    hb_para = this%x(1) 
    !d3_ptr = this%x()
    CALL section_vals_val_set(this%root_section, &
                              "FORCE_EVAL%SUBSYS%TOPOLOGY%COORD_FILE_NAME", &
                              c_val=this%inp_file, error = this%error)
    do i=1,size(dftb3_param)
        CALL section_vals_val_set(this%root_section, &
                                  i_rep_section=i,&
                                  keyword_name="FORCE_EVAL%SUBSYS%KIND%DFTB3_PARAM", &
                                  r_val=dftb3_param(i), error = this%error)
    enddo
    !CALL section_vals_val_set(this%root_section, &
    !                          "FORCE_EVAL%DFT%QS%DFTB%PARAMETER%D3_SCALING", &
    !                          r_vals_ptr=d3_ptr, error = this%error)
    
    CALL section_vals_val_set(this%root_section, &
                              "FORCE_EVAL%DFT%QS%DFTB%PARAMETER%HB_SR_PARAM", &
                              r_val=hb_para, error = this%error)
    
    CALL section_vals_val_set(this%root_section, &
                              "FORCE_EVAL%DFT%CHARGE", &
                              i_val=this%charge, error = this%error)
    !here LSD should be enabled if charge is not 0 
    if(this%charge .ne. 0) then
        CALL section_vals_val_set(this%root_section, &
                              "FORCE_EVAL%DFT%LSD", &
                              l_val=.TRUE., error = this%error)
    else
        CALL section_vals_val_set(this%root_section, &
                              "FORCE_EVAL%DFT%LSD", &
                              l_val=.FALSE., error = this%error)
    endif
    !-------------------create force env------------------------
    CALL create_force_env(this%f_env_id, &
                            input_declaration=this%input_declaration,&
                            input_path=this%input_path, &
                            input=this%root_section, &
                            output_unit=this%iw, &
                            mpi_comm=this%group, &
                            ierr=ierr)
    IF (ierr .ne. 0) THEN
        print *, "Create force env failed",trim(input_file)
        STOP
    ENDIF
    !print *,"after create force env"
    !------------------------------------------------------------
    CALL f_env_get_from_id(this%f_env_id,this%f_env,this%error)
    IF(failure) STOP "paramopt_drive_worker: f_env_add_defaults failed"
    this%force_env => this%f_env%force_env 
    this%qs_env => this%force_env%qs_env
    CALL force_env_get(this%force_env, globenv=this%globenv, subsys=this%subsys, &
                       para_env=this%para_env, error=this%error)
    this%subsys_section => section_vals_get_subs_vals(this%force_env%force_env_section, &
                                                        "SUBSYS",error=this%error)
    CALL calc_energy_force(this%f_env_id, .TRUE. ,ierr)
    CALL get_energy(this%f_env_id,energy,ierr)
    !------------------cleanup------------------------------------------------------
    CALL destroy_force_env(this%f_env_id, ierr)
    IF(ierr /= 0) STOP "paramopt_worker: destroy_force_env failed" 
    !work_section => section_vals_get_subs_vals(this%subsys_section,"COORD",error=this%error)
    !CALL section_vals_remove_values(work_section, this%error)
   END SUBROUTINE paramopt_energy
! *****************************************************************************
!> \brief Finalizes worker for global optimization
!> \param worker ...
!> \author Ole Schuett
! *****************************************************************************
   SUBROUTINE paramopt_worker_finalize(worker)
    TYPE(paramopt_worker_type), INTENT(INOUT)  :: worker

    INTEGER                                  :: ierr

  !   CALL f_env_rm_defaults(worker%f_env, worker%error)
  !   CALL destroy_force_env(worker%f_env_id, ierr)
     IF(ierr /= 0) STOP "glbopt_worker_finalize: destroy_force_env failed"
   END SUBROUTINE paramopt_worker_finalize

! *****************************************************************************
!> \brief Finalizes master for Minima Hopping
!> \param this ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE myfun(x,y,func_nr)
    REAL(KIND=dp)                        :: y
    REAL(KIND=dp)                        :: mysum,A
    REAL(KIND=dp),DIMENSION(:)           :: x
    INTEGER                              :: func_nr
    INTEGER                              :: k,i
    
       ! k = 3
       ! A = 10._dp
       ! mysum = 0._dp
       ! do i=1,k
       !     mysum = mysum + (x(i)**2 - A*cos(2_dp*pi*x(i)))
       ! enddo
       ! y = A*k + mysum
    select case(func_nr)
        case (1)
       ! spehere function
       y = sum(x**2)
        case (2)
       ! Rosenbrock function
       y = (1._dp-x(1))**2 + 100._dp*(x(2)-x(1)**2)**2
   end select
 END SUBROUTINE myfun

END MODULE paramopt_worker

