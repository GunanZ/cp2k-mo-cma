!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines for the Minima Hopping global optimization scheme
!> \author Ole Schuett
! *****************************************************************************
MODULE mocmaes
  USE bibliography,                    ONLY: Goedecker2004,&
                                             cite_reference
  USE cmaes,                           ONLY: cmaes_type,&
                                             cmaes_init,&
                                             cmaes_input,&
                                             update_covariance,&
                                             update_step_size 
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
  USE cp_array_sort_r,                 ONLY: cp_1d_r_sort
#include "../common/cp_common_uses.f90"

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cmaes'

 TYPE :: mo_cmaes_type_p
    TYPE(mo_cmaes_type),POINTER                        :: p
 END TYPE

 TYPE, EXTENDS(cmaes_type)                             :: mo_cmaes_type
   REAL(kind=dp)                                       :: I_dist
   INTEGER                                             :: front
   INTEGER                                             :: d_count
   INTEGER                                             :: s_i_count
   TYPE(mo_cmaes_type_p),DIMENSION(:),ALLOCATABLE      :: s_i !should be length of pop_size
   TYPE(mo_cmaes_type),POINTER                         :: myparent=>null()
   INTEGER                                             :: popsize !offspring + parents
   LOGICAL                                             :: offspring
   INTEGER,DIMENSION(:),ALLOCATABLE                    :: offspring_list
 END TYPE mo_cmaes_type



  PUBLIC mo_cmaes_init
  PUBLIC mo_cmaes_type
  PUBLIC mo_cmaes_type_p
  PUBLIC mo_cmaes_run

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
 SUBROUTINE mo_cmaes_init(mocmaes,input)
    TYPE(cmaes_input)                                  :: input
    TYPE(mo_cmaes_type),DIMENSION(:),TARGET            :: mocmaes
    INTEGER                                            :: n
    INTEGER                                            :: m
    INTEGER                                            :: i,j 
    INTEGER                                            :: popsize
    popsize = size(mocmaes)/2 
    do i=size(mocmaes),1,-1
        mocmaes(i)%id = i
        call cmaes_init(mocmaes(i), input)
        mocmaes(i)%popsize = popsize
        allocate(mocmaes(i)%s_i(size(mocmaes)))
        do j=1,size(mocmaes)
            mocmaes(i)%s_i(j)%p => null()
        enddo
        allocate(mocmaes(i)%offspring_list(popsize/2))
        do j=1,popsize/2
            mocmaes(i)%offspring_list(j) = j
        enddo
        mocmaes(i)%s_i_count = 0 
        mocmaes(i)%d_count = 0
        mocmaes(i)%I_dist = 0._dp
        if (i .le. popsize) then
           mocmaes(i)%offspring = .true.
           mocmaes(i)%myparent => mocmaes(i+popsize)
        else
           mocmaes(i)%offspring = .false.
           mocmaes(i)%myparent => null() 
        endif
     enddo
 END SUBROUTINE mo_cmaes_init

! *****************************************************************************
!> \brief Central steering routine of Minima Hopping
!> \param this ...
!> \param report ...
!> \param cmd ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE mo_cmaes_run(mocmaes)
    implicit none 
    type(mo_cmaes_type),DIMENSION(:),TARGET             :: mocmaes
    type(mo_cmaes_type_p),DIMENSION(:),ALLOCATABLE      :: parents
    type(mo_cmaes_type_p),DIMENSION(:),ALLOCATABLE      :: offspring
    type(mo_cmaes_type),POINTER                         :: a_k,a_k_parents
    REAL(kind=dp)                                       :: lambda_succ
    REAL(kind=dp)                                       :: p_succ
    INTEGER                                             :: accepted
    INTEGER                                             :: popsize
    INTEGER                                             :: k,n,i,j
    LOGICAL                                             :: boolean
    
    !n = cmaes_opt%n
    lambda_succ = 0._dp
    p_succ = 0._dp
    boolean = .false.
    accepted = 0 
    popsize = mocmaes(1)%popsize
    n=0
    k=1
    allocate(parents(popsize),offspring(popsize))
    call update_offspring(mocmaes,offspring,parents)
    call init_selection(mocmaes)
    do i=1,popsize 
        a_k => offspring(i)%p
        a_k_parents => a_k%myparent
            call dominates(a_k,a_k_parents,boolean)
            if (boolean) then
                a_k%p_better = 1._dp
                a_k_parents%p_better = 1._dp
                accepted = accepted + 1
            else
                a_k%p_better = 0._dp
                a_k_parents%p_better = 0._dp
            endif
            !----------------------adapt sigma----------------------------------------
            CALL update_step_size(a_k_parents)  
            CALL update_step_size(a_k)
            !---------------------update covariance----------------------------------
            if (p_succ .gt. 0._dp) then
                call update_covariance(a_k)
            endif
        nullify(a_k,a_k_parents)
    enddo
    !print *,'done init selection',accepted
    call selection(mocmaes)
    call update_offspring(mocmaes,offspring,parents)
    ! not sure if update offspring needed here
    !print *,"before gen pop"
    !--------------------------generating new pop--------------------------------    
    DO j = 1, popsize
             a_k_parents => parents(j)%p 
             a_k => offspring(j)%p
             !what needs to be copied
             a_k%x = a_k_parents%x
             a_k%xold = a_k_parents%xold
             a_k%C = a_k_parents%C
             a_k%A = a_k_parents%A
             a_k%z = a_k_parents%z 
             a_k%sigma = a_k_parents%sigma
             a_k%p_succ= a_k_parents%p_succ
             a_k%pc    = a_k_parents%pc
             !---------------------
             a_k_parents%myparent => null()
             a_k%myparent => a_k_parents
             a_k%id = j
             a_k_parents%id = j+popsize
             a_k%offspring = .TRUE.
             a_k_parents%offspring = .FALSE.
             print *,'oldx',a_k%x
             CALL gen_pop(a_k)
             print *,'newx',a_k%x

     END DO
     deallocate(offspring,parents)
     !print *,"gen pop done"
 END SUBROUTINE mo_cmaes_run

!===========================================================================
SUBROUTINE init_selection(mocmaes)
    implicit none
    TYPE(mo_cmaes_type),DIMENSION(:),TARGET            :: mocmaes
    TYPE(mo_cmaes_type_p),DIMENSION(:),ALLOCATABLE     :: new_set 
    TYPE(mo_cmaes_type_p),DIMENSION(:),ALLOCATABLE     :: sel_set 
    TYPE(mo_cmaes_type_p),DIMENSION(:),ALLOCATABLE     :: tmp_set 
    INTEGER                                            :: offspring
    INTEGER                                            :: tmp_off
    INTEGER                                            :: sel_target
    REAL(kind=dp),DIMENSION(:,:),ALLOCATABLE           :: newsort
    REAL(kind=dp),DIMENSION(:),ALLOCATABLE             :: tmp_array
    INTEGER,DIMENSION(:),ALLOCATABLE                   :: idxsel
    INTEGER                                            :: i,j,k
!------------------------selection here---------------------------------------
    offspring = 0
    tmp_off = 0
    allocate(tmp_set(size(mocmaes)))
    do i = 1,size(mocmaes)
        ! reset all election varaibles
        mocmaes(i)%d_count=0
        mocmaes(i)%I_dist=0._dp
        do j =1,size(mocmaes)
            mocmaes(i)%s_i(j)%p=>null()
            mocmaes(i)%s_i_count=0
            mocmaes(i)%front=0
        enddo
    enddo
    do i=1,size(mocmaes)
        tmp_set(i)%p => mocmaes(i) !the entire population
    enddo
    !---------find pareto front and add to Q (tmp_set)----------------------------------
    CALL pareto_front(input_set=tmp_set,nondominated_set=new_set)
    CALL crowding_distance(new_set)
    !====ad front to Q
    do i = 1,size(new_set)
       new_set(i)%p%front = 1
    enddo
    tmp_off = tmp_off + size(new_set)
    offspring = offspring + size(new_set)
    !---------add more fronts to Q-----------------
    j = 1   !this is also the number of fronts! 
    do while (tmp_off .lt. size(mocmaes)) 
        ! -------copy front to selection set 
        if (allocated(sel_set)) deallocate(sel_set)
        allocate(sel_set(size(new_set)))
        do i = 1, size(new_set)
            sel_set(i)%p => new_set(i)%p
        enddo
!        print *,'size',size(tmp_set),size(sel_set)
        CALL new_front(tmp_set,sel_set,new_set)
        CALL crowding_distance(new_set)
        !print *,'size new set', size(new_set)
        tmp_off = tmp_off + size(new_set)
             do i = 1,size(new_set)
             ! -----ad to Q
             new_set(i)%p%front = j + 1
             enddo
             j = j + 1
             offspring = offspring + size(new_set)
     enddo
    deallocate(tmp_set,new_set)
    end subroutine init_selection

!===========================================================================
    SUBROUTINE pareto_front(input_set,nondominated_set)
         IMPLICIT NONE
         TYPE(mo_cmaes_type),POINTER                    :: a_k,a_k2
         TYPE(mo_cmaes_type_p),DIMENSION(:)             :: input_set
         TYPE(mo_cmaes_type_p),DIMENSION(:),ALLOCATABLE :: nondominated_set
         INTEGER                                        :: nondominated_size
         INTEGER                                        :: i,j,k,t 
         INTEGER                                        :: m
         m      = size(input_set(1)%p%f) 
         ! calculate fron with dominated count d_count and sub_set dominated s_i 
         do i =1,size(input_set)
            a_k => input_set(i)%p
            k = 1
            do j = 1,size(input_set)
                a_k2 => input_set(j)%p
                if (all(a_k%f .lt. a_k2%f)) then
                    a_k%s_i(k)%p => a_k2
                    a_k%s_i_count = a_k%s_i_count + 1
                    k = k +1
                else if (all(a_k%f .gt. a_k2%f)) then
                    a_k%d_count = a_k%d_count + 1
                endif
            enddo 
         enddo
         nondominated_size = 0
         do i =1,size(input_set)
            a_k => input_set(i)%p
            if (a_k%d_count .eq. 0) &
                nondominated_Size = nondominated_size + 1
            a_k%front = a_k%d_count + 1
            if (a_k%s_i_count .gt. size(input_set)) STOP "s_i count not possible"
            if (a_k%d_count .gt. size(input_set)) STOP "d count not possible"
        enddo
        if (allocated(nondominated_set)) deallocate(nondominated_set)
        allocate(nondominated_set(nondominated_size))
        k=1
        do i =1,size(input_set)
            a_k => input_set(i)%p
            if (a_k%d_count .eq. 0) then 
            nondominated_set(k)%p => a_k             
            k = k+1
            endif
        enddo

        END SUBROUTINE pareto_front
!===========================================================================
     SUBROUTINE crowding_distance(input_set)
         IMPLICIT NONE
         TYPE(mo_cmaes_type_p),DIMENSION(:)            :: input_Set
         REAL(kind=dp),DIMENSION(:,:),ALLOCATABLE    :: newsort
         REAL(kind=dp),DIMENSION(:),POINTER          :: fmax,fmin
         REAL(kind=dp),DIMENSION(:),ALLOCATABLE      :: maxvals,minvals
         REAL(kind=dp),DIMENSION(:),ALLOCATABLE      :: tmp
         REAL(kind=dp),DIMENSION(:),ALLOCATABLE      :: mydeviation,myIdist,maxIdist
         REAL(kind=dp)                               :: mymax,mymin,mydev
         INTEGER,DIMENSION(:),ALLOCATABLE            :: idxsel
         INTEGER                                     :: nondominated_size
         INTEGER                                     :: sel_target
         INTEGER                                     :: i,j,k,n 
         INTEGER                                     :: j_m,j_p
         LOGICAL                                     :: variance_sorting
         LOGICAL                                     :: kernel_density
         REAL(kind=dp),DIMENSION(:),ALLOCATABLE :: x,f,dens
         REAL(kind=dp)                          :: h,var,mean
         INTEGER                                :: m
         nondominated_size=size(input_set)
         m=size(input_set(1)%p%f)
         do i = 1,nondominated_size
             input_set(i)%p%I_dist = 0._dp
             !input_set(i)%p%I_dist = 1._dp
         enddo       
         !------------------------------------------------
         allocate(newsort(nondominated_size,m))
         allocate(idxsel(nondominated_size))
         allocate(x(nondominated_size))
         allocate(f(nondominated_size))
         allocate(dens(nondominated_size-2))
         ! crowding-distance-assignment
         do i = 1,nondominated_size
             newsort(i,:) = input_set(i)%p%f(:)
         enddo
        !-----------------------------------------------
         do k = 1,m     ! num objectives
           x = newsort(:,k)
           CALL cp_1d_r_sort(x,size(x),idxsel)
           f = x
           do i=1,size(x)
               x(i)=f(idxsel(i))
           enddo
           !print *,'x',x
           !print *,'f',f
           !============calculate using crowding dist NAGA-II===========
            input_set(idxsel(1))%p%I_dist = huge(1._dp) 
            input_set(idxsel(nondominated_size))%p%I_dist = huge(1._dp) 
            mymax = maxval(x)
            mymin = minval(x)
            mydev = mymax - mymin 
            do i = 2,nondominated_size-1
               j = idxsel(i)
               j_p = idxsel(i+1)
               j_m = idxsel(i-1)
               if (mydev .gt. huge(1._dp)) mydev = 1._dp
               input_set(j)%p%I_dist = input_set(j)%p%I_dist + & 
                                (input_set(j_p)%p%f(k)- &
                                 input_set(j_m)%p%f(k))/& 
                                mydev
             if (isNaN(input_set(j)%p%I_dist)) input_set(j)%p%I_dist = 0.0_dp
             if ((input_set(j)%p%I_dist) .gt. huge(1._dp)) input_set(j)%p%I_dist = 0.0_dp
             !print *,'Idist'  ,(input_set(j)%I_dist)
             enddo
           !============calculate using kde===========
           ! n = nondominated_size
           ! if (n .gt. 2) then
           !     call kde(f(2:n-1),f(2:n-1),dens,h)
           !     !mean = sum(dens)/size(dens)
           !     !var  = variance(dens)
           !     !dens = (dens-mean)/var
           !     !dens = (dens-minval(dens))/(maxval(dens)-minval(dens))
           !     !print *,'x',f(2:n-1)
           !     !print *,'dens',dens
           !     !print *,'maxval',maxval(dens),minval(dens)
           !     input_set(idxsel(1))%p%I_dist = 0.0_dp
           !     input_set(idxsel(nondominated_size))%p%I_dist = 0.0_dp
           !     do i = 1,nondominated_size-2 
           !        j = idxsel(i)
           !        input_set(j)%p%I_dist = -(abs(input_set(j)%p%I_dist*abs(dens(i))))
           !      if (isNaN(input_set(j)%p%I_dist)) input_set(j)%p%I_dist = -huge(1._dp)
           !     enddo
           ! else
           !     do i = 1,nondominated_size 
           !        input_set(i)%p%I_dist = 0._dp
           !     enddo
           ! endif
           !=================================================================
         enddo
         !do i =1,nondominated_size
         !       !print *,'kde',f(i),dens(i),h
         !       j = idxsel(i)
         !       print *,'kde',i,input_set(j)%p%I_dist
         !enddo
         deallocate(x,f,dens,newsort,idxsel) 
    END SUBROUTINE crowding_distance
!===========================================================================
    SUBROUTINE new_front(global_set,input_set,output_set)
       IMPLICIT NONE
       TYPE(mo_cmaes_type_p),DIMENSION(:),TARGET      :: input_Set
       TYPE(mo_cmaes_type_p),DIMENSION(:),INTENT(IN)  :: global_Set
       TYPE(mo_cmaes_type_p),DIMENSION(:),ALLOCATABLE :: output_Set
       TYPE(mo_cmaes_type_p),DIMENSION(:),ALLOCATABLE :: tmp_Set
       INTEGER                                      :: i,j,k,idx
       INTEGER                                      :: numrows
       TYPE(mo_cmaes_type),POINTER                        :: a_k
       TYPE(mo_cmaes_type),POINTER                        :: a_k2
       CHARACTER(len=*), PARAMETER :: routineN = 'new_front', &
                    routineP = moduleN//':'//routineN
       
       allocate(tmp_set(size(global_set)))
       numrows = size(input_set)
       if (allocated(output_set)) deallocate(output_set)
!       do i = 1,numrows
!            a_k => input_set(i)%p
!            print *,'all stuff',a_k%id,a_k%s_i_count 
!            do j=1,a_k%s_i_count; print *,'idx',a_k%s_i(j)%p%id ; enddo
!       enddo
       k = 1
       do i = 1,numrows
           a_k => input_set(i)%p
           !print *,'new_front',a_k%f,a_k%s_i,a_k%s_i_count,a_k%d_count
           do j = 1,a_k%s_i_count
               !------ go through subset ------
               a_k2 => a_k%s_i(j)%p
               a_k2%d_count = a_k2%d_count - 1
               if (a_k2%d_count .eq. 0) then
                   tmp_set(k)%p => a_k2
                   k = k + 1
               endif
           enddo
       enddo
       if (allocated(output_set)) deallocate(output_set)
       allocate(output_set(k-1))
       do i=1,size(output_set)
           output_set(i)%p => tmp_set(i)%p
       enddo 
       deallocate(tmp_set)
    END SUBROUTINE new_front
    !===========================================================================
    subroutine selection(mocmaes)
    implicit none
    TYPE(mo_cmaes_type),DIMENSION(:),TARGET                    ::mocmaes
    TYPE(mo_cmaes_type_p),DIMENSION(:),ALLOCATABLE              :: tmp_set 
    CHARACTER(len=*), PARAMETER :: routineN = 'paramopt_selection', &
      routineP = moduleN//':'//routineN
    REAL(kind=dp),DIMENSION(:),ALLOCATABLE                    :: newsort
    INTEGER,DIMENSION(:),ALLOCATABLE                          :: idxsel 
    integer                                                   :: idx
    integer                                                   :: i,j,k,n
    integer                                                   :: fsize
    integer                                                   :: psize
    integer                                                   :: front
    INTEGER                                                   :: popsize
     TYPE(cp_error_type)                                      :: error
!    !--------everything could be done more efficient with pointers!!!--------------------
    front = 1
    psize = 0
    popsize = mocmaes(1)%popsize
    k = 0
    !----reset offspring varaiable---------
    do i=1,size(mocmaes)
        mocmaes(i)%offspring = .true.
    enddo
    do while (psize .lt. popsize)
        !-------get front size-------------
        fsize = 0 
        do i=1,size(mocmaes)
          if (front .eq. mocmaes(i)%front) then
          fsize = fsize + 1
          endif
        enddo
        !-------------front size lower of what needed?----------
        if (fsize+psize .le. popsize) then
         ! add all individual from front to Q    
            do i=1,size(mocmaes)
              if (front .eq. mocmaes(i)%front) then
              !sel_set(k)%p => cmaes_opt%population(i)
              k = k + 1
              mocmaes(i)%offspring = .false.
              endif
            enddo
            psize = psize + fsize
       else
             ! check crowding distance and add n best
             allocate(tmp_set(fsize))
             !----add to temporary set
             allocate(newsort(fsize))
             allocate(idxsel(fsize))
             j=1
             do i=1,size(mocmaes)
               if (front .eq. mocmaes(i)%front) then
               tmp_set(j)%p => mocmaes(i)
               newsort(j) = mocmaes(i)%I_dist
               j = j + 1
               endif
             enddo
           !----sort according to I_dist-----------
           CALL cp_1d_r_sort(newsort,size(newsort),idxsel)
           !----add to sel set Q------------
           do i = 0,popsize-psize-1
               idx = size(idxsel)-i
               tmp_set(idxsel(idx))%p%offspring = .false.
               k = k + 1
           enddo
           psize = psize + (i+1)
           deallocate(tmp_set,newsort,idxsel)
        endif 
        front = front + 1 
     enddo
    !----------------------was selection okay?----------------------
    j=0
    do i = 1,size(mocmaes)
        if (mocmaes(i)%offspring) j = j + 1
    enddo
    !call cp_assert( k-1 .eq. popsize,cp_failure_level,cp_assertion_failed,&
    !                   routineP, cp_to_string(k)//" k ne lambda_mo "//&
    !                   cp_to_string(popsize),error)
    call cp_assert( j .eq. popsize,cp_failure_level,cp_assertion_failed,&
                       routineP, cp_to_string(j)//" j ne lambda_mo "//&
                       cp_to_string(popsize),error)
    !print *,'selection done'
    end subroutine selection
!===========================================================================
    subroutine dominates(a_k,a_k_parents,boolean)
    type(mo_cmaes_type)                             :: a_k
    type(mo_cmaes_type)                             :: a_k_parents
    LOGICAL                                         :: boolean
    IF (a_k%front .lt. a_k_parents%front) then
        boolean = .true. 
    else if (a_k%front .eq. a_k_parents%front) then
            !---depends on 
            if (a_k%I_dist .gt. a_k_parents%I_dist) then
                boolean = .true.
            else
                boolean = .false.
            endif
    else
        boolean = .false.   
    endif
    ! ugly hack 
    if (any(a_k%f .gt. 1.e100_dp) ) boolean = .false.
    end subroutine dominates
!===========================================================================
    subroutine gen_pop(a_k)
    implicit none
    type(mo_cmaes_type),POINTER                    :: a_k
    TYPE(cp_error_type),POINTER                    :: error
    REAL(kind=dp)                                  :: tmp
    INTEGER                                        :: i,n
    !generate new population 
     n=a_k%n
     do i=1,n; a_k%z(i) = next_random_number(a_k%rng_gaussian_stream,error=error); enddo
     print *,'z',a_k%z
     CALL DGEMV('N',n,n,1.0d0,a_k%A,n,a_k%z,1,0.0d0,a_k%Az,1)
     a_k%xold = a_k%x
     !print *,'p-better',a_k%p_better,'f',a_k%f,'fold',a_k%fold
     !print *,'xold',a_k%x,'sigma',a_k%sigma,'Az',a_k%Az
     !print *,"C",a_k%C
     a_k%x = a_k%x+ a_k%sigma * a_k%Az
     do i=1,size(a_k%x)
        if (a_k%x(i) .lt. a_k%l_bound(i)) a_k%x(i) = a_k%l_bound(i)
        if (a_k%x(i) .gt. a_k%u_bound(i)) a_k%x(i) = a_k%u_bound(i)
     enddo

     !!yannick specific hack DANGEROUSE 
     !if (abs(a_k%x(9)) .gt. abs(a_k%x(11))) then
     !   tmp = a_k%x(9)
     !   a_k%x(9) = a_k%x(11)
     !   a_k%x(11) = tmp
     !endif


     a_k%dx = norm2(a_k%x-a_k%xold)
    
    end subroutine gen_pop

! *****************************************************************************
!> \brief Initializes master for Minima Hopping
!> \param this ...
!> \param glbopt_section ...
!> \param n_workers ...
!> \param iw ...
!> \param error ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE update_offspring(mocmaes,offspring,parents)
    TYPE(mo_cmaes_type),DIMENSION(:),TARGET            :: mocmaes
    TYPE(mo_cmaes_type_p),DIMENSION(:)                 :: parents
    TYPE(mo_cmaes_type_p),DIMENSION(:)                 :: offspring
    INTEGER                                            :: i,j,k 
    INTEGER,DIMENSION(:),ALLOCATABLE                  :: offspring_list
    allocate(offspring_list(size(mocmaes)/2))
    do i=1,size(mocmaes)/2
       nullify(parents(i)%p,offspring(i)%p)
    enddo
    i=1
    j=1
    do k=1,size(mocmaes)
        if (mocmaes(k)%offspring) then
           offspring(i)%p => mocmaes(k) 
           offspring_list(i) = mocmaes(k)%Id
           i = i + 1
        else
           parents(j)%p => mocmaes(k)
           parents(j)%p%myparent => null()
           j = j + 1
        endif
    enddo
    do k=1,size(mocmaes)
        mocmaes(k)%offspring_list = offspring_list
    enddo
    deallocate(offspring_list)
END SUBROUTINE update_offspring

    subroutine kde(x,z,f,h)
    implicit none
    real(dp),dimension(:),intent(in)   :: x,z
    real(dp),dimension(:),allocatable  :: temp,tempf
    real(dp),dimension(:),intent(out)  :: f
    real(dp),intent(out)               :: h
    real(dp)                           :: mean,sd,xj,xi
    real(dp)                           :: from,to,step
    integer                            :: i,j,n,m,locmin(1)

    m=512
    allocate(temp(m),tempf(m))
    n = size(x)
    m = size(temp)
    from = 0.0_dp
    to   = 9.5_dp
    step = (to-from)/m
    temp(1) = from
    do i=2,m
        temp(i) = temp(i-1) + step
    enddo
    mean = sum(x)/n
    sd = sqrt(variance(x))
    !h = 1.06_dp*sd*(n**(-1/5))
    h = 0.9_dp*sd/(n**(1/5))
    if (h .eq.0 ) h = 1._dp
    !h = 0.1 
    tempf = 0._dp

    do j=1,m
        xj = temp(j)
 
        do i=1,n
            xi = x(i)
            tempf(j) = tempf(j) + K((xj-xi)/h)
        enddo

    enddo
    tempf = 1.0/(n*h)*tempf

    do i=1,m
    write(99,*) temp(i),tempf(i)
    enddo
    do i=1,size(z)
    xi = z(i)
    locmin = minloc(abs(xi-temp))
    f(i) = tempf(locmin(1))
    !print *,"minloc",locmin,xi,tempf(locmin)
    enddo
    end


    function K(u) 
        implicit none
        real(dp)    :: K,u,pi
        
        pi = 4.*atan(1.)
        K = 1.0_dp/sqrt(2.0_dp*pi)*exp(-0.5_dp*u**2)
    end
    function variance(array)
        implicit none
        real(dp),intent(in),dimension(:) :: array
        real(dp)                         :: variance,x
        x = sum(array)/size(array)
        variance = sum((array-x)**2)/size(array)
    end

END MODULE mocmaes
