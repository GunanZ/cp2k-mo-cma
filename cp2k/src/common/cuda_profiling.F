!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  routines for profiling cuda
!> \par History
!>      05.2013 created
!> \author Ole Schuett
! *****************************************************************************
MODULE cuda_profiling
  USE kinds,                           ONLY: default_string_length,&
                                             int_8

 !$ USE OMP_LIB

 IMPLICIT NONE

 PRIVATE

 PUBLIC  :: cuda_nvtx_init, cuda_nvtx_range_push, cuda_nvtx_range_pop, cuda_mem_info

#if defined( __CUDA_PROFILING )

  INTERFACE
     FUNCTION cuda_nvtx_range_push_dc(message) RESULT(level) &
       BIND(C, name="cuda_nvtx_range_push_cu") 
       USE ISO_C_BINDING
    CHARACTER(kind=C_CHAR), DIMENSION(*), &
      INTENT(IN)                             :: message
    INTEGER(KIND=C_INT)                      :: level

     END FUNCTION cuda_nvtx_range_push_dc
  END INTERFACE

  INTERFACE
     FUNCTION cuda_nvtx_range_pop_dc() RESULT(level) &
       BIND(C, name="cuda_nvtx_range_pop_cu") 
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT)                      :: level

     END FUNCTION cuda_nvtx_range_pop_dc
  END INTERFACE

  INTERFACE
    FUNCTION cudaMemGetInfo_cu (free, total) RESULT(istat) &
      BIND(C, name="cudaMemGetInfo")
      USE ISO_C_BINDING
    INTEGER(KIND=C_SIZE_T), INTENT(INOUT)    :: free, total
    INTEGER(KIND=C_INT)                      :: istat

    END FUNCTION cudaMemGetInfo_cu
  END INTERFACE

  INTERFACE
     SUBROUTINE cuda_nvtx_name_osthread_cu(name)&
          BIND(C, name="cuda_nvtx_name_osthread_cu")
        USE ISO_C_BINDING
    CHARACTER(KIND=C_CHAR), DIMENSION(*)     :: name

     END SUBROUTINE cuda_nvtx_name_osthread_cu
  END INTERFACE


#endif

 CONTAINS

#if defined( __CUDA_PROFILING )
 ! *****************************************************************************
! *****************************************************************************
!> \brief ...
! *****************************************************************************
 SUBROUTINE cuda_nvtx_init()
    USE ISO_C_BINDING
    CHARACTER(len=default_string_length)     :: threadname
    INTEGER                                  :: ithread

    !$OMP PARALLEL default (none), private (ithread,threadname)
    ithread = 0
    !$ ithread = OMP_GET_THREAD_NUM ()
    WRITE(threadname, "(I3,A,I2,A)") ithread
    CALL cuda_nvtx_name_osthread_cu(TRIM(threadname)//c_null_char)
    !$OMP END PARALLEL
  END SUBROUTINE cuda_nvtx_init

 ! *****************************************************************************
! *****************************************************************************
!> \brief ...
!> \param routineN ...
! *****************************************************************************
  SUBROUTINE cuda_nvtx_range_push(routineN)
    CHARACTER(LEN=*), INTENT(IN)             :: routineN

    INTEGER                                  :: level

    level = cuda_nvtx_range_push_dc(TRIM(routineN)//CHAR(0))
  END SUBROUTINE cuda_nvtx_range_push

 ! *****************************************************************************
! *****************************************************************************
!> \brief ...
! *****************************************************************************
  SUBROUTINE cuda_nvtx_range_pop()
    INTEGER                                  :: level

     level = cuda_nvtx_range_pop_dc()
  END SUBROUTINE cuda_nvtx_range_pop

  ! *****************************************************************************
! *****************************************************************************
!> \brief ...
!> \param free ...
!> \param total ...
! *****************************************************************************
  SUBROUTINE cuda_mem_info(free, total)
    USE ISO_C_BINDING
    INTEGER(KIND=int_8), INTENT(OUT)         :: free, total

    INTEGER(KIND=C_INT)                      :: istat
    INTEGER(KIND=C_SIZE_T)                   :: free_c, total_c

    istat = cudaMemGetInfo_cu(free_c, total_c)
    IF(istat /= 0) STOP "cuda_mem_info: failed"
    free = free_c
    total = total_c
  END SUBROUTINE cuda_mem_info

#else
 ! *****************************************************************************
! *****************************************************************************
!> \brief ...
! *****************************************************************************
  SUBROUTINE cuda_nvtx_init()
  ! nothing todo
  END SUBROUTINE cuda_nvtx_init

 ! *****************************************************************************
! *****************************************************************************
!> \brief ...
!> \param routineN ...
! *****************************************************************************
  SUBROUTINE cuda_nvtx_range_push(routineN)
    CHARACTER(LEN=*), INTENT(IN)             :: routineN

    STOP "cuda_nvtx_range_push: __CUDA_PROFILING not compiled in."
  END SUBROUTINE cuda_nvtx_range_push

 ! *****************************************************************************
! *****************************************************************************
!> \brief ...
! *****************************************************************************
  SUBROUTINE cuda_nvtx_range_pop()
     STOP "cuda_nvtx_range_push: __CUDA_PROFILING not compiled in."
  END SUBROUTINE cuda_nvtx_range_pop

  ! *****************************************************************************
! *****************************************************************************
!> \brief ...
!> \param free ...
!> \param total ...
! *****************************************************************************
  SUBROUTINE cuda_mem_info(free, total)
    INTEGER(KIND=int_8), INTENT(OUT)         :: free, total

    free  = 0
    total = 0
  END SUBROUTINE cuda_mem_info
#endif


END MODULE cuda_profiling
