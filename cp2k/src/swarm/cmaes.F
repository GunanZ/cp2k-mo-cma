!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines for the Minima Hopping global optimization scheme
!> \author Ole Schuett
! *****************************************************************************
MODULE cmaes
  USE bibliography,                    ONLY: Goedecker2004,&
                                             cite_reference
  USE glbopt_history,                  ONLY: history_add,&
                                             history_finalize,&
                                             history_fingerprint,&
                                             history_fingerprint_match,&
                                             history_fingerprint_type,&
                                             history_init,&
                                             history_lookup,&
                                             history_type
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mathconstants,                   ONLY: pi
  USE parallel_rng_types,              ONLY: rng_stream_type, &
                                             create_rng_stream, &
                                             delete_rng_stream,&
                                             next_random_number, &
                                             GAUSSIAN, &
                                             UNIFORM
  USE physcon,                         ONLY: kelvin
  USE swarm_message,                   ONLY: swarm_message_add,&
                                             swarm_message_get,&
                                             swarm_message_type

#include "../common/cp_common_uses.f90"

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cmaes'


 PUBLIC :: cmaes_run
 PUBLIC :: cmaes_type
 PUBLIC :: cmaes_input
 PUBLIC :: cmaes_init
 PUBLIC :: update_covariance
 PUBLIC :: update_step_size

 TYPE cmaes_type
   REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE          :: A,B,C,D,tmp_mat,BD
   REAL(KIND=dp), DIMENSION(:), POINTER              :: x
   REAL(KIND=dp), DIMENSION(:), POINTER                :: u_bound
   REAL(KIND=dp), DIMENSION(:), POINTER                :: l_bound
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE          :: xold
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE            :: z 
   REAL(KIND=dp), DIMENSION(:), POINTER              :: f
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE          :: fold
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE            :: tmp_vec
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE            :: Az
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE            :: diagD
   REAL(KIND=dp)                                       :: p_succ 
   INTEGER                                             :: p_better 
   REAL(KIND=dp)                                       :: sigma 
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE            :: pc 
   INTEGER                                             :: n
   INTEGER                                             :: iteration
   !----step size----
   REAL(KIND=dp)                                       :: damp
   REAL(KIND=dp)                                       :: p_target
   REAL(KIND=dp)                                       :: cp
   !----cov adapt----
   REAL(KIND=dp)                                       :: cc
   REAL(KIND=dp)                                       :: ccov
   REAL(KIND=dp)                                       :: p_tresh
   !------------------
   TYPE(rng_stream_type),POINTER                       :: rng_uniform_stream => Null()
   TYPE(rng_stream_type),POINTER                       :: rng_gaussian_stream => Null()
   !------------------ 
   REAL(KIND=dp)                                       :: det
   REAL(KIND=dp)                                       :: dx
 CONTAINS
   PROCEDURE                                           :: cmaes_init 
 END TYPE cmaes_type
 
 TYPE cmaes_input
   REAL(KIND=dp)                                       :: input_seed
   REAL(KIND=dp)                                       :: step_size
   REAL(KIND=dp), DIMENSION(:), POINTER                :: x
   REAL(KIND=dp), DIMENSION(:), POINTER                :: f
   INTEGER                                             :: maxiter
   INTEGER                                             :: m !num of obj
   INTEGER                                             :: n !num of dimensions
   INTEGER                                             :: iteration
   LOGICAL                                             :: randinput
   REAL(KIND=dp), DIMENSION(:), POINTER                :: u_bound
   REAL(KIND=dp), DIMENSION(:), POINTER                :: l_bound
 END TYPE cmaes_input

 CONTAINS
 
! *****************************************************************************
!> \brief Initializes master for Minima Hopping
!> \param this ...
!> \param glbopt_section ...
!> \param n_workers ...
!> \param iw ...
!> \param error ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE cmaes_init(this, input)
    CLASS(cmaes_type)                                  :: this
    TYPE(cmaes_input)                                  :: input
     
    INTEGER                                            :: n
    INTEGER                                            :: m
    INTEGER                                            :: i
    !---input init
    input%iteration = 0
    
    !---cma_inp init
    n = input%n 
    m = input%m
    this%n = n
    this%x => input%x
    this%f => input%f
    allocate(this%A(n,n))
    allocate(this%B(n,n))
    allocate(this%C(n,n))
    allocate(this%D(n,n))
    allocate(this%tmp_mat(n,n))
    allocate(this%BD(n,n))
    allocate(this%pc(n)) 
    allocate(this%xold(n)) 
    allocate(this%z(n)) 
    allocate(this%tmp_vec(n))
    allocate(this%Az(n))
    allocate(this%diagD(n))
    !------------------------
    allocate(this%fold(m))
    !-------------------------
    this%damp = 1._dp + real(n)/2._dp
    this%p_target = 1._dp/(5._dp + 1._dp/2._dp)
    this%cp = this%p_target/(2._dp+this%p_target)
    this%cc = 2._dp/(real(n)+2._dp)
    this%ccov = 2._dp/(real(n)**2+6._dp)
    this%ccov = 0.1_dp 
    this%p_tresh = 0.44_dp
    !-------------------------
    this%iteration = 0
    !------------------------
    this%tmp_vec = 0._dp
    this%p_succ = this%p_target
    this%pc = 0._dp
    if (input%step_size .eq. 0._dp) then
        this%sigma = 0.5_dp 
    else
        this%sigma = input%step_size
    endif
    this%A = 0._dp
    this%B = 0._dp
    this%C = 0._dp
    this%D = 0._dp
    this%tmp_mat = 0._dp
    this%BD = 0._dp
    do i=1,n; this%A(i,i) = 1._dp; enddo
    do i=1,n; this%B(i,i) = 1._dp; enddo
    do i=1,n; this%C(i,i) = 1._dp; enddo
    do i=1,n; this%D(i,i) = 1._dp; enddo
    do i=1,n; this%BD(i,i) = 1._dp; enddo
    CALL create_rng(this,input)
    !-------------------------
    this%l_bound => input%l_bound
    this%u_bound => input%u_bound
    !---------------------------
    !if (input%randinput) then
    !    print *,'generating random x'
    !    do i=1,n; this%x(i) = next_random_number(this%rng_uniform_stream,error=error); enddo  
    !    this%x = var_trans(this%x,this%u_bound,this%l_bound)
    !else
    this%x = input%x
    !endif
    !this%xtrans = var_back_trans(this%x,this%u_bound,this%l_bound)
    !this%xold = this%x
    !----are we doing everything in 0 ..1 ?
    this%xold = this%x
    !this%x = this%xtrans
    !-------------------------
    this%f = huge(1._dp) 
    this%fold = huge(1._dp)
  !  print *,'step_size cma',this%sigma
 END SUBROUTINE cmaes_init

! *****************************************************************************
!> \brief Initializes master for Minima Hopping
!> \param this ...
!> \param glbopt_section ...
!> \param n_workers ...
!> \param iw ...
!> \param error ...
!> \author Ole Schuett
! *****************************************************************************
! *****************************************************************************
!> \brief Central steering routine of Minima Hopping
!> \param this ...
!> \param report ...
!> \param cmd ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE cmaes_run(input,population)
    TYPE(cmaes_input)                                  :: input
    CLASS(cmaes_type),ALLOCATABLE,DIMENSION(:),TARGET             :: population
    CLASS(cmaes_type),POINTER                          :: a_k,a_k_parent => Null()
    INTEGER                                            :: iteration
    INTEGER                                            :: i
    INTEGER                                            :: j
    INTEGER                                            :: n
    TYPE(cp_error_type)                                :: error
   
    input%iteration = input%iteration + 1
    a_k => population(1)
    !a_k_parent => population(2)
    n = a_k%n
    if (a_k%f(1) .lt. epsilon(0._dp)) then
       print *, epsilon(0._dp),"target accuracy reached",input%iteration
       STOP
    endif

    !check if better
    CALL update_success(a_k%f,a_k%fold,a_k%p_better)
    CALL update_step_size(a_k)

    CALL write_out(a_k)
    if ( a_k%p_better .gt. 0) then
        !update cov
        CALL update_covariance(a_k)
        a_k%fold = a_k%f
    else
        a_k%x = a_k%xold
    endif

    !generate new population 
     do i=1,n; a_k%z(i) = next_random_number(a_k%rng_gaussian_stream,error=error); enddo
     CALL DGEMV('N',n,n,1.0d0,a_k%A,n,a_k%z,1,0.0d0,a_k%Az,1)
     a_k%xold = a_k%x
     !print *,'p-better',a_k%p_better,'f',a_k%f,'fold',a_k%fold
     !print *,'xold',a_k%x,'sigma',a_k%sigma,'Az',a_k%Az
     !print *,"C",a_k%C
     a_k%x = a_k%x+ a_k%sigma * a_k%Az
     a_k%dx = norm2(a_k%x-a_k%xold)

     !print *,'xnew',a_k%x 
     
     !bounds check
     do i=1,n
        if (a_k%x(i) .lt. 0._dp) a_k%x(i) = 0._dp 
        if (a_k%x(i) .gt. 1._dp) a_k%x(i) = 1._dp
     enddo 
     !do i=1,n
     !   if (a_k%x(i) .lt. a_k%l_bound(i)) a_k%x(i) = a_k%l_bound(i)
     !   if (a_k%x(i) .gt. a_k%u_bound(i)) a_k%x(i) = a_k%u_bound(i)
     !enddo 

 END SUBROUTINE cmaes_run

! *****************************************************************************
!> \brief Central steering routine of Minima Hopping
!> \param this ...
!> \param report ...
!> \param cmd ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE update_success(f,fold,better)
    REAL(KIND=dp),DIMENSION(:)               :: f,fold
    INTEGER,INTENT(OUT)                      :: better
    !if (all(a_k%f .lt. a_k_parent%f)) then
    if (ALL((f .lt. fold))) then
        better = 1
    else
        better = 0
    endif
 END SUBROUTINE update_success

! *****************************************************************************
!> \brief Central steering routine of Minima Hopping
!> \param this ...
!> \param report ...
!> \param cmd ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE update_step_size(this)
    CLASS(cmaes_type)                                   :: this
    this%p_succ = (1._dp - this%cp)*this%p_succ + real(this%cp*this%p_better)
    this%sigma = this%sigma*exp((this%p_succ - this%p_target)/ &
                         (this%damp*(1._dp-this%p_target)))
!    this%sigma = 0.5_dp/sqrt(minval(this%diagD))
 END SUBROUTINE update_step_size
 ! *****************************************************************************
!> \brief Central steering routine of Minima Hopping
!> \param this ...
!> \param report ...
!> \param cmd ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE update_covariance(this)
    CLASS(cmaes_type)                                   :: this
    INTEGER                                            :: n,i,j
    INTEGER                                            :: info

    n=this%n
 
    if (this%p_succ .lt. this%p_tresh) then
       this%pc = (1._dp - this%cc)*this%pc + sqrt(this%cc*(2._dp-this%cc))*(this%Az)
       CALL dgemm('N','T',n,n,1,1.0d0,this%pc,n,this%pc,n,0.0d0,this%tmp_mat,n) 
       this%C  = (1._dp - this%ccov)*this%C + this%ccov*this%tmp_mat
    else
       this%pc = (1._dp-this%cc)*this%pc
       CALL dgemm('N','T',n,n,1,1.0d0,this%pc,n,this%pc,n,0.0d0,this%tmp_mat,n)
       this%C  = (1._dp-this%ccov)*this%C + this%ccov*(this%tmp_mat + this%cc*(2._dp - this%cc)*this%C)
    endif
    
    !enforce sym
    CALL tool_symmatrix(this%C,n,this%tmp_mat)
    !-----------------------------------
    !cholesky decompose C
    this%tmp_mat = this%C
        !print *,'original' 
        !do, i=1,n
        !     write(*,*) ( this%tmp_mat(i,j), j=1,n )
        !enddo  
    CALL DPOTRF( 'L', n, this%tmp_mat, n, INFO )
    IF (info .ne. 0) STOP "eigendecomp A"
    this%A = 0._dp
    do i=1,n
        do j=1,n
            if (i .ge. j) then
                this%A(i,j) = this%tmp_mat(i,j)
            endif
        enddo
    enddo
        !this%tmp_mat = matmul(this%A,transpose(this%A))
        !print *,'cholesky'
        !do, i=1,n
        !     write(*,*) ( this%tmp_mat(i,j), j=1,n )
        !enddo  
    !-----------------------------------
    !eigenvalue decomposition
    this%tmp_mat = this%C 
    CALL tool_eigendecomp(this%tmp_mat,n,this%D,this%B,info) 
    IF (info .ne. 0) STOP "eigendecomp BD"
    do, i=1,n
         this%diagD(i) = this%D(i,i)
    enddo  
    this%BD = matmul(this%B,sqrt(this%D))

        !this%tmp_mat = matmul(this%BD,transpose(this%BD))
        !print *,'eigenvalue'
        !do, i=1,n
        !     write(*,*) ( this%tmp_mat(i,j), j=1,n )
        !enddo  

END SUBROUTINE update_covariance
! *****************************************************************************
!> \brief Central steering routine of Minima Hopping
!> \param this ...
!> \param report ...
!> \param cmd ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE create_rng(this,input)
    CLASS(cmaes_type)                                   :: this
    TYPE(cmaes_input)                                  :: input
    REAL(kind=dp),DIMENSION(3,2)                       :: seed
    INTEGER                                            :: time(3)
    TYPE(cp_error_type)                                :: error

    if (input%input_seed .gt. 0) then
        !deterministic seed
        seed(:,:) = &
            RESHAPE( (/input%input_seed*42.0_dp, &
                       input%input_seed*54.0_dp, &
                       input%input_seed*63.0_dp, &
                       input%input_seed*98.0_dp, &
                       input%input_seed*10.0_dp, &
                       input%input_seed*2.0_dp/),&
                     (/ 3, 2 /) )        
    else
        !random seed from itime
        call itime(time)
        seed(:,1) = REAL(time)
        seed(:,2) = (/ 1._dp, 2._dp, 3._dp /)
    endif 
    
    CALL create_rng_stream(this%rng_gaussian_stream,name='gaussian',distribution_type=GAUSSIAN,seed=seed,error=error)    
    CALL create_rng_stream(this%rng_uniform_stream,name='uniform',distribution_type=UNIFORM,seed=seed,error=error)    
 END SUBROUTINE create_rng

! *****************************************************************************
!> \brief Central steering routine of Minima Hopping
!> \param this ...
!> \param report ...
!> \param cmd ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE write_out(a_k)
    CLASS(cmaes_type),POINTER                           :: a_k
    INTEGER                                            :: i,j,n
    
    !print *,'x',a_k%x,a_k%xold
    !print *,"better",a_k%f,a_k%fold,a_k%p_better,a_k%p_succ,a_k%sigma,a_k%dx
    !do, i=1,n
    !     print *, ( a_k%A(i,j), j=1,n )
    !enddo
    n = size(a_k%A,dim=1) 
    open(unit=1,file="A.out") 
    do, i=1,n
         write(1,*) ( a_k%A(i,j), j=1,n )
    enddo 
    close(1)
    open(unit=1,file="C.out") 
    do, i=1,n
    enddo 
    close(1)
        open(unit=1,file="D.out") 
    do, i=1,n
         write(1,*) ( a_k%D(i,j), j=1,n )
    enddo 
    close(1)
    open(unit=1,file="B.out") 
    do, i=1,n
         write(1,*) ( a_k%B(i,j), j=1,n )
    enddo 
    close(1)
    open(unit=1,file="BD.out") 
    do, i=1,n
         write(1,*) ( a_k%BD(i,j), j=1,n )
    enddo 
    close(1)
END SUBROUTINE write_out


! *****************************************************************************
!> \brief Central steering routine of Minima Hopping
!> \param this ...
!> \param report ...
!> \param cmd ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE tool_symmatrix(C,n,triuC)
    IMPLICIT NONE
  
    !-------------------------------------------------------------------------
    !  Parameters
    !-------------------------------------------------------------------------
    INTEGER, INTENT(in)                           :: n
    REAL(dp),DIMENSION(n,n),INTENT(inout)         :: C
    REAL(dp),DIMENSION(n,n),INTENT(inout)         :: triuC

    !-------------------------------------------------------------------------
    !  Local Variables
    !-------------------------------------------------------------------------
    REAL(dp),DIMENSION(n,n)                       :: triuC1,triuC2
    LOGICAL,DIMENSION(n,n)                        :: mask
    INTEGER                                       :: i,k
    
    !-----------------------------------------------------------------------
    !  Enforce symmetry
    !-----------------------------------------------------------------------
    triuC1 = 0.
    triuC2 = 0.
    mask = .FALSE.
    DO k = 1, n       ! Create mask for upper triangular matrix
      DO i = 1, n
          IF(k .GE. i) THEN
              mask(i,k) = .TRUE.
          ELSE
              mask(i,k) = .FALSE.
          END IF
            if (isNaN(C(k,i))) then
                print *,'matrix has NaN entries'
                STOP
            endif
      END DO
    END DO
    
    WHERE(mask)
      triuC1 = C          ! Get upper triangular matrix of C
    END WHERE             ! See Matlab: triu(X)
    
    DO i = 1, n           ! set new mask, see triu(X,k), here: k=1
      mask(i,i) = .FALSE.
    END DO
    
    WHERE(mask)
      triuC2 = C          ! Get second triangular matrix
    END WHERE
    
    C = triuC1 + transpose(triuC2)    ! Enforce symmetry
      triuC = triuC1

    RETURN
  END SUBROUTINE
! *****************************************************************************
!> \brief Central steering routine of Minima Hopping
!> \param this ...
!> \param report ...
!> \param cmd ...
!> \author Ole Schuett
! *****************************************************************************

      SUBROUTINE tool_eigendecomp(A,n,eVals,eVecs,info)
      IMPLICIT NONE
      
      !-------------------------------------------------------------------------
      !  Parameters
      !-------------------------------------------------------------------------
      INTEGER, INTENT(in)							:: n
      REAL(dp),DIMENSION(n,n),INTENT(in)			:: A
      REAL(dp),DIMENSION(n,n),INTENT(out)			:: eVals
      REAL(dp),DIMENSION(n,n),OPTIONAL,INTENT(out)	:: eVecs
      INTEGER,INTENT(OUT)                           :: info
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------
      EXTERNAL 									:: dsyev 
      !-------------------------------------------------------------------------
      !  Local Variables
      !-------------------------------------------------------------------------
      CHARACTER(len=1)							:: job
      INTEGER									:: lWork
      REAL(dp),ALLOCATABLE,DIMENSION(:)			:: work
      INTEGER									:: i
      REAL(dp),DIMENSION(n)						:: eValWork
      INTEGER                                   :: allocStat
     
      !-------------------------------------------------------------------------
      !  Specify job (only compute Eigenvalues('N') or also Eigenvectors('V') )
      !-------------------------------------------------------------------------
      IF(present(eVecs)) THEN
        job = 'V'
      ELSE
        job = 'N'
      END IF
      
      !-------------------------------------------------------------------------
      !  Make a copy of A (so it isnt changed unintentionally)
      !-------------------------------------------------------------------------
      eVecs = A
      
      !-------------------------------------------------------------------------
      !  Set length of work array
      !-------------------------------------------------------------------------
      !lWork = 2*n*n+6*n+1
      lWork = -1
      ALLOCATE(work(1),stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating work(1)'
      CALL DSYEV(job,'U',n,eVecs,n,eValWork,work,-1,info)
      i = work(1)
      DEALLOCATE(work)
      
      ALLOCATE(work(i),stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating work'
      !-------------------------------------------------------------------------
      !  Call LAPACK Routine
      !-------------------------------------------------------------------------
      CALL DSYEV(job,'U',n,eVecs,n,eValWork,work,i,info)
      !-------------------------------------------------------------------------
      !  Error management
      !-------------------------------------------------------------------------
      IF(info .GT. 0) WRITE(*,*) 'EVD failed to converge'
      IF(info .LT. 0) WRITE(*,*) 'EVD: Illegal value'
      
      !-------------------------------------------------------------------------
      !  Set output variables
      !-------------------------------------------------------------------------
      DO i = 1,n
        eVals(i,i) = eValWork(i) ! Convert Vector to Matrix
      END DO


      RETURN
      END SUBROUTINE
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
        x = xtrans*(u_bound-l_bound) + l_bound
    end function var_trans
END MODULE cmaes

