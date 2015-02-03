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
                                             paramopt_do_grid,&
                                             paramopt_do_mocma
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
   REAL(KIND=dp),DIMENSION(:),POINTER       :: func
   INTEGER                                  :: n
   INTEGER                                  :: m
   TYPE(cp_error_type)                      :: error
   TYPE(cp_para_env_type), POINTER          :: para_env
   INTEGER                                  :: f_env_id
   INTEGER                                  :: method
   INTEGER                                  :: inp_num
   INTEGER                                  :: charge 
   INTEGER                                  :: pop_id 
   INTEGER                                  :: converged=0
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

        worker%m = 3
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
        allocate(worker%func(worker%m))
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
   REAL(kind=dp)                            :: x,x1,x2,h
   REAL(kind=dp),DIMENSION(:),POINTER       :: func
   INTEGER                                  :: Id


    SELECT CASE (worker%method)
       CASE(paramopt_do_cma,paramopt_do_powell,paramopt_do_grid)
            CALL swarm_message_get(cmd, "command", command)
            IF(TRIM(command) == "function") THEN
                nullify(worker%x)
                CALL swarm_message_get(cmd, "x", worker%x)
                CALL swarm_message_get(cmd, "ref_func", worker%ref_func)
                !CALL run_mdgopt(worker, cmd, report)
                !------------do calc here------------
                CALL myfun(worker%x,worker%func(1),worker%ref_func)
                ! blablabla
                !------------------------------------ 
                CALL swarm_message_add(report, "status", "done") 
                CALL swarm_message_add(report, "result", worker%func(1))
            ELSE IF (TRIM(command) == "energy") THEN
               nullify(worker%x)
               CALL swarm_message_get(cmd, "pop_id", worker%pop_id)
               CALL swarm_message_get(cmd, "inp_num", worker%inp_num)
               CALL swarm_message_get(cmd, "inp_file", worker%inp_file)
               CALL swarm_message_get(cmd, "charge", worker%charge)
               CALL swarm_message_get(cmd, "x", worker%x)
               !-----------------------------------------------
               CALL paramopt_energy(worker,cmd,command,report)
               !-----------------------------------------------
               CALL swarm_message_add(report, "status", "done")
               CALL swarm_message_add(report, "pop_id", worker%pop_id)
               CALL swarm_message_add(report, "inp_num", worker%inp_num) 
               CALL swarm_message_add(report, "result", worker%func)
               CALL swarm_message_add(report, "converged", worker%converged)
            ELSE
               print *,"Worker: received unkown command:"
               STOP 
            END IF
       CASE (paramopt_do_mocma)
            CALL swarm_message_get(cmd, "command", command)
          IF(TRIM(command) == "function") THEN
            nullify(worker%x)
            CALL swarm_message_get(cmd, "inp_num", worker%inp_num)
            CALL swarm_message_get(cmd, "pop_id", worker%pop_id)
            CALL swarm_message_get(cmd, "x", worker%x)
            CALL swarm_message_get(cmd, "ref_func", worker%ref_func)
            write(worker%iw,*) "this is x in worker",worker%x
            !------------do calc here------------
            !!CALL myfun(worker%x,y,worker%ref_func)
            !!Zitzler¿Deb¿Thiele's function N. 1
            !worker%x = abs(worker%x)
            x1 = worker%x(1)
            x2 = worker%x(2)
            !!h = 1._dp - sqrt((x1)/(g(worker%x)))
            !h = 1._dp - (x1/g(worker%x))**2
            !worker%func(1) = (x1)
            !worker%func(2) = (g(worker%x)*h)
            worker%func(1) = 0.5_dp*(x1**2+x2**2) + sin(x1**2+x2**2)
            worker%func(2) = (3_dp*x1-2_dp*x2+4_dp)**2/8_dp + (x1-x2+1_dp)**2/27_dp + 15_dp
            worker%func(3) = 1_dp/(x1**2+x2**2+1_dp)-1.1_dp*exp(-(x1**2+x2**2))
            write(worker%iw,*) 'pop_id',worker%pop_id
            write(worker%iw,*) 'this is x',worker%x,'and f',worker%func
            !------------------------------------ 
            CALL swarm_message_add(report, "status", "done") 
            CALL swarm_message_add(report, "inp_num", worker%inp_num)
            CALL swarm_message_add(report, "pop_id", worker%pop_id)
            CALL swarm_message_add(report, "result", worker%func)
          ELSE IF (TRIM(command) == "energy") THEN
            nullify(worker%x)
            CALL swarm_message_get(cmd, "inp_num", worker%inp_num)
            CALL swarm_message_get(cmd, "pop_id", worker%pop_id)
            CALL swarm_message_get(cmd, "inp_file", worker%inp_file)
            CALL swarm_message_get(cmd, "charge", worker%charge)
            CALL swarm_message_get(cmd, "x", worker%x)
            !------------do calc here------------
            CALL paramopt_energy(worker,cmd,command,report)
            !------------------------------------ 
            CALL swarm_message_add(report, "status", "done") 
            CALL swarm_message_add(report, "inp_num", worker%inp_num)
            CALL swarm_message_add(report, "pop_id", worker%pop_id)
            CALL swarm_message_add(report, "result", worker%func)
            CALL swarm_message_add(report, "converged", worker%converged)
          ENDIF
       CASE DEFAULT
        STOP "Unkown glbopt_method"
    END SELECT
    CONTAINS 
        REAL FUNCTION g(x)
            REAL(dp),DIMENSION(:)                    :: x
            integer                                  :: n
            n = size(x)
            g = 1._dp + 9._dp/(real(n)-1._dp)*(sum(x(2:n)))
        END FUNCTION g  
   END SUBROUTINE paramopt_worker_execute

! *****************************************************************************
!> \brief Finalizes worker for global optimization
!> \param worker ...
!> \author Ole Schuett
! *****************************************************************************
   SUBROUTINE paramopt_energy(this, cmd, command, report)
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
    INTEGER                                  :: natoms
    REAL(KIND=dp),DIMENSION(:),POINTER       :: x => Null()
    REAL(KIND=dp),POINTER                    :: energy
    REAL(KIND=dp),POINTER                    :: rmsd
    REAL(KIND=dp),DIMENSION(:),POINTER       :: hb_para => Null()
    REAL(KIND=dp),DIMENSION(:),POINTER       :: dftb3_param => Null()
    REAL(KIND=dp),DIMENSION(:),POINTER       :: d3_ptr => Null()
    REAL(KIND=dp),DIMENSION(:),POINTER       :: xopt_ptr => Null()
    REAL(KIND=dp),DIMENSION(:),ALLOCATABLE   :: pos,pos2
    REAL(dp), DIMENSION(3, 3)                :: U
    REAL(dp), DIMENSION(3)                   :: center1, center2
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: g
    energy => this%func(1)
    rmsd => this%func(2)
    energy = 0._dp
    rmsd = 0._dp
    !CALL swarm_message_get(cmd,"input",input_file)
    !CALL swarm_message_get(cmd,"coord",coord_file)
    !CALL swarm_message_get(cmd,"num_inp",i_inp)
    !CALL swarm_message_get(cmd,"iter",iter)
    !CALL swarm_message_get(cmd,"x",this%x)
    write(this%iw,*) "this is x in worker",this%x
    !print *,'this is x in worker',this%x
    !print *,'this is charge', this%charge
    !print *,'this is inp',this%inp_file
    !----for D3 optimization---------------------------
!    allocate(d3_ptr(3)) !gets deallocated automatically ??? at least it seems
!    this%D3 = (/0.5_dp,this%x(4:5)/)
!!    this%D3 = (/1._dp,this%x(1:2)/)
!    !this%D3 = this%x(1:3) 
!    d3_ptr = this%D3
!    !----for dftb3 param optimization-----------------
!    ! NEEDS TO BE ADJUSTED !!!!!!!!!!!!!!!!!
!    !allocate(dftb3_param(4))
!    allocate(dftb3_param(2))
!    dftb3_param = this%x(3:4)
!    !---for hb param optimization---------------------
!    allocate(hb_para(3))
!    hb_para = 0._dp
!    hb_para(1) = this%x(1)
!    hb_para(2) = this%x(2)
!    hb_para(3) = this%x(3)
   !-------------------------------------------------
    CALL section_vals_val_set(this%root_section, &
                              "FORCE_EVAL%SUBSYS%TOPOLOGY%COORD_FILE_NAME", &
                              c_val=this%inp_file, error = this%error)
!---------------------------------------------------------------------------
!    do i=1,size(dftb3_param)
!        CALL section_vals_val_set(this%root_section, &
!                                  i_rep_section=i,&
!                                  keyword_name="FORCE_EVAL%SUBSYS%KIND%DFTB3_PARAM", &
!                                  r_val=dftb3_param(i), error = this%error)
!    enddo
!    CALL section_vals_val_set(this%root_section, &
!                              "FORCE_EVAL%DFT%QS%DFTB%PARAMETER%D3_SCALING", &
!                              r_vals_ptr=d3_ptr, error = this%error)
!    
!    CALL section_vals_val_set(this%root_section, &
!                              "FORCE_EVAL%DFT%QS%DFTB%PARAMETER%HB_SR_PARAM", &
!                              r_vals_ptr=hb_para, error = this%error)
!    
!--------------------------------------------------------------------------
    allocate(xopt_ptr(size(this%x)))
    xopt_ptr = this%x
    CALL section_vals_val_set(this%root_section, &
                              "FORCE_EVAL%DFT%QS%SCPTB%XOPT", &
                              r_vals_ptr=xopt_ptr, error= this%error)
!--------------------------------------------------------------------------                              
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
    CALL get_natom(this%f_env_id,natoms,ierr)
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
    if (this%force_env%qs_env%converged) then
        this%converged = 1
    else
        this%converged = 0
    endif
    !-----------------geo opt------------------------------------
    if (this%method .eq. paramopt_do_mocma .or. &
        this%method .eq. paramopt_do_grid) then
        if (this%converged .eq. 1) then
          ALLOCATE(pos2(3*natoms))
          ALLOCATE(pos(3*natoms))
          ALLOCATE(g(2,3*natoms))
          CALL get_pos(this%f_env_id, pos, 3*natoms, ierr)
          CALL cp_geo_opt(this%force_env, this%globenv, eval_opt_geo=.TRUE.,error=this%error)
          CALL get_energy(this%f_env_id,energy,ierr)
          CALL get_pos(this%f_env_id, pos2, 3*natoms, ierr)
          CALL get_rmsd(natoms,RESHAPE(pos,(/ 3, natoms /)), &
                                 RESHAPE(pos2,(/ 3, natoms /) ), &
                                 0, U,center2,center1, rmsd,.FALSE.,g)
          !----------------calc energy--------------------
           if (this%force_env%qs_env%converged) then
                this%converged = 1
           else
                this%converged = 0
                energy = huge(1._dp)
           endif
           deallocate(pos,pos2,g)!
        else 
            energy = huge(1._dp)
        endif
    endif
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
       if (size(x) .gt. 2) STOP "Rosenbrock selected x shold be of length 2"
       y = (1._dp-x(1))**2 + 100._dp*(x(2)-x(1)**2)**2
   end select
 END SUBROUTINE myfun
!-----------------------------------------------------------------------
SUBROUTINE get_rmsd(n, coord1, coord2, option, U, x_center, y_center, & 
     error, calc_g, g)
!-----------------------------------------------------------------------
!  This subroutine calculates the least square rmsd of two coordinate
!  sets coord1(3,n) and coord2(3,n) using a method based on quaternion.
!  If option=1, then the rotation matrix U and the centers of coord are 
!  returned.
!-----------------------------------------------------------------------
! if calc_g == .true., derivative of RMSD with respect to coord1
! is returned
!-----------------------------------------------------------------------

    INTEGER, INTENT(in)                      :: n
    REAL(dp), DIMENSION(:, :), INTENT(in)    :: coord1, coord2
    INTEGER, INTENT(in)                      :: option
    REAL(dp), DIMENSION(:, :), INTENT(out)   :: U
    REAL(dp), DIMENSION(3), INTENT(out)      :: x_center, y_center
    REAL(dp), INTENT(out)                    :: error
    LOGICAL, INTENT(in)                      :: calc_g
    REAL(dp), INTENT(out)                    :: g(:,:)

    INTEGER                                  :: i, j
    REAL(dp)                                 :: lambda, tmp(3), x_norm, y_norm
    REAL(dp), DIMENSION(3, 3)                :: Rmatrix
    REAL(dp), DIMENSION(3, n)                :: x, y
    REAL(dp), DIMENSION(4)                   :: q
    REAL(dp), DIMENSION(4, 4)                :: S
    REAL(dp), DIMENSION(n)                   :: xi, yi

! make copies of the original coordinates

  x(:,1:n) = coord1(:,1:n)
  y(:,1:n) = coord2(:,1:n)

  ! calculate the barycenters, centroidal coordinates, and the norms
  x_norm = 0.0d0
  y_norm = 0.0d0
  DO i = 1, 3
     xi(:) = x(i,:)
     yi(:) = y(i,:)
     x_center(i) = SUM(xi(1:n))/DBLE(n)
     y_center(i) = SUM(yi(1:n))/DBLE(n)
     xi(:) = xi(:) - x_center(i)
     yi(:) = yi(:) - y_center(i)
     x(i,:) = xi(:)
     y(i,:) = yi(:)
     x_norm = x_norm + DOT_PRODUCT(xi, xi)
     y_norm = y_norm + DOT_PRODUCT(yi, yi)
  END DO

  ! calculate the R matrix
  DO i = 1, 3
     DO j = 1, 3
        Rmatrix(i,j) = DOT_PRODUCT(x(i,:),y(j,:))
     END DO
  END DO

  ! S matrix
  S(1, 1) = Rmatrix(1, 1) + Rmatrix(2, 2) + Rmatrix(3, 3)
  S(2, 1) = Rmatrix(2, 3) - Rmatrix(3, 2)
  S(3, 1) = Rmatrix(3, 1) - Rmatrix(1, 3)
  S(4, 1) = Rmatrix(1, 2) - Rmatrix(2, 1)

  S(1, 2) = S(2, 1)
  S(2, 2) = Rmatrix(1, 1) - Rmatrix(2, 2) - Rmatrix(3, 3)
  S(3, 2) = Rmatrix(1, 2) + Rmatrix(2, 1)
  S(4, 2) = Rmatrix(1, 3) + Rmatrix(3, 1)

  S(1, 3) = S(3, 1)
  S(2, 3) = S(3, 2)
  S(3, 3) =-Rmatrix(1, 1) + Rmatrix(2, 2) - Rmatrix(3, 3)
  S(4, 3) = Rmatrix(2, 3) + Rmatrix(3, 2)

  S(1, 4) = S(4, 1)
  S(2, 4) = S(4, 2)
  S(3, 4) = S(4, 3)
  S(4, 4) =-Rmatrix(1, 1) - Rmatrix(2, 2) + Rmatrix(3, 3) 

  ! Calculate eigenvalues and eigenvectors, and 
  ! take the maximum eigenvalue lambda and the corresponding eigenvector q.
  CALL dstmev(S, lambda, q)

  IF (option == 1) THEN
     ! convert quaternion q to rotation matrix U
     CALL rotation_matrix(q, U)
  END IF

  ! RMS Deviation
  error = SQRT(MAX(0.0d0,((x_norm+y_norm)-2.0d0*lambda))/DBLE(n))

  IF (calc_g) THEN
     DO i = 1, n
        DO j = 1, 3
           tmp(:) = MATMUL(TRANSPOSE(U(:,:)), y(:,i))
           g(j,i) = (x(j,i) - tmp(j))/(error*DBLE(n))
        END DO
     END DO
  END IF

END SUBROUTINE get_rmsd
!-----------------------------------------------------------------------
SUBROUTINE rotation_matrix(q, U)
!-----------------------------------------------------------------------
! This subroutine constructs rotation matrix U from quaternion q.
!-----------------------------------------------------------------------

    REAL(dp), DIMENSION(:), INTENT(in)       :: q
    REAL(dp), DIMENSION(:, :), INTENT(out)   :: U

    REAL(dp)                                 :: b0, b1, b2, b3, q0, q00, q01, &
                                                q02, q03, q1, q11, q12, q13, &
                                                q2, q22, q23, q3, q33

  q0 = q(1)
  q1 = q(2)
  q2 = q(3)
  q3 = q(4)

  b0 = 2.0d0*q0
  b1 = 2.0d0*q1
  b2 = 2.0d0*q2
  b3 = 2.0d0*q3

  q00 = b0*q0-1.0d0
  q01 = b0*q1
  q02 = b0*q2
  q03 = b0*q3

  q11 = b1*q1
  q12 = b1*q2
  q13 = b1*q3  

  q22 = b2*q2
  q23 = b2*q3

  q33 = b3*q3 

  U(1,1) = q00+q11
  U(1,2) = q12-q03
  U(1,3) = q13+q02

  U(2,1) = q12+q03
  U(2,2) = q00+q22
  U(2,3) = q23-q01

  U(3,1) = q13-q02
  U(3,2) = q23+q01
  U(3,3) = q00+q33

END SUBROUTINE rotation_matrix
!-----------------------------------------------------------------------
SUBROUTINE DSTMEV(A,lambda,evec)
!-----------------------------------------------------------------------
! a simple subroutine to compute the leading eigenvalue and eigenvector
! of a symmetric, traceless 4x4 matrix A by an inverse power iteration:
! (1) the matrix is converted to tridiagonal form by 3 Givens
! rotations;  V*A*V' = T
! (2) Gershgorin's theorem is used to estimate a lower
! bound for the leading negative eigenvalue:
! lambda_1 > g=min(T11-t12,-t21+T22-t23,-t32+T33-t34,-t43+T44)
!          =
! where tij=abs(Tij)
! (3) Form the positive definite matrix 
!     B = T-gI
! (4) Use svd (algorithm svdcmp from "Numerical Recipes")
!     to compute eigenvalues and eigenvectors for SPD matrix B
! (5) Shift spectrum back and keep leading singular vector
!     and largest eigenvalue.
! (6) Convert eigenvector to original matrix A, through 
!     multiplication by V'.  
!-----------------------------------------------------------------------
    REAL(dp), DIMENSION(4, 4)                :: A
    REAL(dp)                                 :: lambda
    REAL(dp), DIMENSION(4)                   :: evec

    INTEGER                                  :: i
    INTEGER, DIMENSION(1)                    :: max_loc
    REAL(dp), DIMENSION(4)                   :: SW
    REAL(dp), DIMENSION(4, 4)                :: SV, T, V
    REAL(dp), DIMENSION(8)                   :: rv1

! must be an array
!-----------------------------------------------------------------------
!(I).   Convert to tridiagonal form, keeping similarity transform
! (a product of 3 Givens rotations)

  CALL givens4(A,T,V)

  !-----------------------------------------------------------------------
  !(II)Estimate lower bound of smallest eigenvalue by Gershgorin's theorem
  lambda=MIN(T(1,1)-ABS(T(1,2)),-ABS(T(2,1))+T(2,2)-ABS(T(2,3)),&
       -ABS(T(3,2))+T(3,3)-ABS(T(3,4)),-ABS(T(4,3))+T(4,4))
  !-----------------------------------------------------------------------
  !(III). Form positive definite matrix     T <== lambda*I - T
  DO i = 1,4
     T(i,i) = T(i,i)-lambda
  ENDDO
  !-----------------------------------------------------------------------
  !(IV). Compute singular values/vectors of SPD matrix B
  CALL svdcmp(4,T,4,4,SW,SV,rv1)
  !-----------------------------------------------------------------------
  !(V). Shift spectrum back
  max_loc = MAXLOC(SW) 
  lambda = SW(max_loc(1)) + lambda
  !lambda = SW(1) + lambda
  !-----------------------------------------------------------------------
  !(VI). Convert eigenvector to original coordinates: (V is transposed!)
  evec = MATMUL(V,SV(:,max_loc(1)))
  !write(*,*)'-----------------------------------------------------------'
  !write(*,*) 'lambda = ', lambda,'  eigenvector:  '
  !write(*,99) evec
  !write(*,*)'-----------------------------------------------------------'
  !-----------------------------------------------------------------------
  !99 format(1x,4(d19.13,1x))

END SUBROUTINE dstmev
!-----------------------------------------------------------------------
SUBROUTINE givens4(S,T,V)
!-----------------------------------------------------------------------
    REAL(dp), DIMENSION(4, 4), INTENT(in)    :: S
    REAL(dp), DIMENSION(4, 4), INTENT(out)   :: T, V

    REAL(dp)                                 :: c1, c1c2, c2, c3, r1, r2, r3, &
                                                s1, s1c2, s2, s3

!double precision :: pythag
! external        pythag
!performs givens rotations to reduce symmetric 4x4 matrix to tridiagonal

  T=S; V = 0.d0
  !-----------------------------------------------------------------------
  !Zero out entries T(4,1) and T(1,4)
  ! compute cos and sin of rotation angle in the 3-4 plane
  r1 = pythag(T(3,1),T(4,1))
  IF(r1 .NE. 0.d0) THEN
     c1 = T(3,1)/r1; s1 = T(4,1)/r1
     V(3,3) = c1  ; V(3,4) = s1
     V(4,3) =-s1  ; V(4,4) = c1
     T(3,1) = r1  ; T(4,1) = 0.d0
     T(3:4,2:4) = MATMUL(V(3:4,3:4),T(3:4,2:4))
     T(1:2,3:4) = TRANSPOSE(T(3:4,1:2))
     T(3:4,3:4) = MATMUL(T(3:4,3:4),TRANSPOSE(V(3:4,3:4)))
  ELSE
     c1 = 1.d0; s1 = 0.d0
  ENDIF
  !-----------------------------------------------------------------------
  !Zero out entries T(3,1) and T(1,3)
  ! compute cos and sin of rotation angle in the 2-3 plane
  r2 = pythag(T(3,1), T(2,1))
  IF(r2 .NE. 0.d0) THEN
     c2 = T(2,1)/r2; s2 = T(3,1)/r2
     V(2,2) = c2  ; V(2,3) = s2
     V(3,2) =-s2  ; V(3,3) = c2
     T(2,1) = r2  ; T(3,1) = 0.d0
     T(2:3,2:4) = MATMUL(V(2:3,2:3),T(2:3,2:4))
     T(1,2:3)   = T(2:3,1);  T(4,2:3) = T(2:3,4)
     T(2:3,2:3) = MATMUL(T(2:3,2:3),TRANSPOSE(V(2:3,2:3)))
  ELSE
     c2 = 1.d0; s2 = 0.d0
  ENDIF
  !-----------------------------------------------------------------------
  !Zero out entries T(4,2) and T(2,4)
  ! compute cos and sin of rotation angle in the 3-4 plane
  r3 = pythag(T(4,2), T(3,2))
  IF(r3 .NE. 0.d0) THEN
     c3 = T(3,2)/r3; s3 = T(4,2)/r3
     V(3,3) = c3  ; V(3,4) = s3
     V(4,3) =-s3  ; V(4,4) = c3
     T(3,2) = r3  ; T(4,2) = 0.d0
     T(3:4,3:4) = MATMUL(V(3:4,3:4),T(3:4,3:4))
     T(1:2,3:4) = TRANSPOSE(T(3:4,1:2))
     T(3:4,3:4) = MATMUL(T(3:4,3:4),TRANSPOSE(V(3:4,3:4)))
  ELSE
     c3 = 1.d0; s3 = 0.d0
  ENDIF
  !-----------------------------------------------------------------------
  !Compute net rotation matrix (accumulate similarity for evec. computation)
  ! To save transposing later, This is the transpose!
  V(1,1)=1.d0; V(1,2:4) = 0.d0; V(2:4,1) = 0.d0
  V(2,2) = c2;  V(3,2) = c1*s2 ; V(4,2) = s1*s2; c1c2 = c1*c2; s1c2=s1*c2
  V(2,3) = -s2*c3 ; V(3,3) = c1c2*c3-s1*s3 ; V(4,3) =  s1c2*c3+c1*s3
  V(2,4) =  s2*s3 ; V(3,4) =-c1c2*s3-s1*c3 ; V(4,4) = -s1c2*s3+c1*c3
  !-----------------------------------------------------------------------
  !write(*,*) (V(1:4,i) - W(1:4,i),i=1,4)
END SUBROUTINE givens4
!-----------------------------------------------------------------------
SUBROUTINE svdcmp(mmax,a,m,n,w,v,rv1)
!-----------------------------------------------------------------------
    INTEGER                                  :: mmax
    REAL(dp)                                 :: a(mmax,*)
    INTEGER                                  :: m, n
    REAL(dp)                                 :: w(*), v(mmax,*), rv1(*)

    INTEGER                                  :: i, its, j, jj, k, l, nm
    REAL(dp)                                 :: anorm, c, f, g, h, s, scale, &
                                                x, y, z

!,pythag

  g = 0.0d0
  scale = 0.0d0
  anorm = 0.0d0
  DO i = 1, n
     l=i+1
     rv1(i)=scale*g
     g=0.0d0
     s=0.0d0
     scale=0.0d0
     IF(i.le.m)THEN
        DO k=i,m
           scale=scale+ABS(a(k,i))
        END DO
        IF(scale.ne.0.0d0)THEN
           DO k=i,m
              a(k,i)=a(k,i)/scale
              s=s+a(k,i)*a(k,i)
           END DO
           f=a(i,i)
           g=-SIGN(SQRT(s),f)
           h=f*g-s
           a(i,i)=f-g
           DO j=l,n 
              s=0.0d0
              DO k=i,m
                 s=s+a(k,i)*a(k,j)
              END DO
              f=s/h
              DO k=i,m
                 a(k,j)=a(k,j)+f*a(k,i)
              END DO
           END DO
           DO k=i,m
              a(k,i)=scale*a(k,i)
           END DO
        ENDIF
     ENDIF
     w(i)=scale *g
     g=0.0d0
     s=0.0d0
     scale=0.0d0
     IF((i.le.m).AND.(i.ne.n))THEN
        DO k=l,n
           scale=scale+ABS(a(i,k))
        END DO
        IF(scale.ne.0.0d0)THEN
           DO k=l,n
              a(i,k)=a(i,k)/scale
              s=s+a(i,k)*a(i,k)
           END DO
           f=a(i,l)
           g=-SIGN(SQRT(s),f)
           h=f*g-s
           a(i,l)=f-g
           DO k=l,n
              rv1(k)=a(i,k)/h
           END DO
           DO j=l,m
              s=0.0d0
              DO k=l,n
                 s=s+a(j,k)*a(i,k)
              END DO
              DO k=l,n
                 a(j,k)=a(j,k)+s*rv1(k)
              END DO
           END DO
           DO k=l,n
              a(i,k)=scale*a(i,k)
           END DO
        ENDIF
     ENDIF
     anorm=MAX(anorm,(ABS(w(i))+ABS(rv1(i))))
  END DO

  DO i = n, 1, -1
     IF(i .LT. n) THEN
        IF(g.ne.0.0d0)THEN
           DO j=l,n
              v(j,i)=(a(i,j)/a(i,l))/g
           END DO
           DO j=l,n
              s=0.0d0
              DO k=l,n
                 s=s+a(i,k)*v(k,j)
              END DO
              DO k=l,n
                 v(k,j)=v(k,j)+s*v(k,i)
              END DO
           END DO
        ENDIF
        DO j=l,n
           v(i,j)=0.0d0
           v(j,i)=0.0d0
        END DO
     ENDIF
     v(i,i)=1.0d0
     g=rv1(i)
     l=i
  END DO
!-----------------------------------------------------------------------
  DO i = MIN(m,n), 1, -1
     l=i+1
     g=w(i)
     DO j=l,n
        a(i,j)=0.0d0
     END DO
     IF(g.ne.0.0d0)THEN
        g=1.0d0/g
        DO j=l,n
           s=0.0d0
           DO k=l,m
              s=s+a(k,i)*a(k,j)
           END DO
           f=(s/a(i,i))*g
           DO k=i,m
              a(k,j)=a(k,j)+f*a(k,i)
           END DO
        END DO
        DO j=i,m
           a(j,i)=a(j,i)*g
        END DO
     ELSE
        DO j= i,m
           a(j,i)=0.0d0
        END DO
     ENDIF
     a(i,i)=a(i,i)+1.0d0
  END DO
!-----------------------------------------------------------------------
  DO k=n,1,-1
     DO its=1,30
        DO l=k,1,-1
           nm=l-1
           IF((ABS(rv1(l))+anorm).eq.anorm)  GOTO 2
           IF((ABS(w(nm))+anorm).eq.anorm)  GOTO 1
        END DO
1       c=0.0d0
        s=1.0d0
        DO i=l,k
           f=s*rv1(i)
           rv1(i)=c*rv1(i)
           IF((ABS(f)+anorm).eq.anorm) GOTO 2
           g=w(i)
           h=pythag(f,g)
           w(i)=h
           h=1.0d0/h
           c= (g*h)
           s=-(f*h)
           DO j=1,m     
              y=a(j,nm)
              z=a(j,i)
              a(j,nm)=(y*c)+(z*s)
              a(j,i)=-(y*s)+(z*c)
           END DO
        END DO
2       z=w(k)
        IF(l .EQ. k)THEN
           IF(z.lt.0.0d0)THEN
              w(k)=-z
              DO j=1,n
                 v(j,k)=-v(j,k)
              END DO
           ENDIF
           GOTO 3
        ENDIF
        IF(its.eq.30) THEN
           WRITE(*,*) 'no convergence in svdcmp'
           STOP
        ENDIF
        x=w(l)
        nm=k-1
        y=w(nm)
        g=rv1(nm)
        h=rv1(k)
        f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0d0*h*y)
        g=pythag(f,1.0d0)
        f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x
        c=1.0d0
        s=1.0d0
        DO j=l,nm
           i=j+1
           g=rv1(i)
           y=w(i)
           h=s*g
           g=c*g    
           z=pythag(f,h)
           rv1(j)=z
           c=f/z
           s=h/z
           f= (x*c)+(g*s)
           g=-(x*s)+(g*c)
           h=y*s
           y=y*c
           DO jj=1,n
              x=v(jj,j)
              z=v(jj,i)
              v(jj,j)= (x*c)+(z*s)
              v(jj,i)=-(x*s)+(z*c)
           END DO
           z=pythag(f,h)
           w(j)=z
           IF(z.ne.0.0d0)THEN
              z=1.0d0/z
              c=f*z
              s=h*z
           ENDIF
           f= (c*g)+(s*y)
           x=-(s*g)+(c*y)
           DO jj=1,m
              y=a(jj,j)
              z=a(jj,i)
              a(jj,j)= (y*c)+(z*s)
              a(jj,i)=-(y*s)+(z*c)
           END DO
        END DO
        rv1(l)=0.0d0
        rv1(k)=f       
        w(k)=x
     END DO
3    CONTINUE
  END DO
!-----------------------------------------------------------------------

END SUBROUTINE svdcmp
!-----------------------------------------------------------------------
FUNCTION pythag(a,b)
!-----------------------------------------------------------------------

    REAL(dp)                                 :: a, b, pythag

    REAL(dp)                                 :: absa, absb

  absa=ABS(a)
  absb=ABS(b)
  IF(absa.gt.absb)THEN
     pythag=absa*dsqrt(1.0d0+(absb/absa)**2)
  ELSE
     IF(absb.eq.0.0d0)THEN
        pythag=0.0d0
     ELSE
        pythag=absb*dsqrt(1.0d0+(absa/absb)**2)
     ENDIF
  ENDIF

END FUNCTION pythag
END MODULE paramopt_worker

