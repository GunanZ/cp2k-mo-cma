!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Defines the basic variable types
!> \author vw
! *****************************************************************************
MODULE ma_kinds

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: dp
  PUBLIC :: int_4
  PUBLIC :: int_size
  PUBLIC :: default_string_length

  INTEGER, PARAMETER :: real_4 = SELECTED_REAL_KIND ( 6, 30 )
  INTEGER, PARAMETER :: real_8 = SELECTED_REAL_KIND ( 14, 200 )
  INTEGER, PARAMETER :: real_4_size = 4
  INTEGER, PARAMETER :: real_8_size = 8

  ! default integer size
  INTEGER, PARAMETER :: int_size = BIT_SIZE(0)/8

  INTEGER, PARAMETER :: int_1  = SELECTED_INT_KIND(2)
  INTEGER, PARAMETER :: int_1_size = BIT_SIZE(INT(0, int_1))/8

  INTEGER, PARAMETER :: int_2  = SELECTED_INT_KIND(4)
  INTEGER, PARAMETER :: int_2_size = BIT_SIZE(INT(0, int_2))/8

  INTEGER, PARAMETER :: int_4 = SELECTED_INT_KIND(5)
  INTEGER, PARAMETER :: int_4_size = BIT_SIZE(INT(0, int_4))/8

  INTEGER, PARAMETER :: int_8 = SELECTED_INT_KIND(10)
  INTEGER, PARAMETER :: int_8_size = BIT_SIZE(INT(0, int_8))/8

  INTEGER, PARAMETER :: sp = real_4
  INTEGER, PARAMETER :: dp = real_8

  INTEGER, PARAMETER :: default_string_length=80

END MODULE ma_kinds

