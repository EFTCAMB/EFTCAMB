!----------------------------------------------------------------------------------------
!
! This file is part of EFTCAMB.
!
! Copyright (C) 2013-2023 by the EFTCAMB authors
!
! The EFTCAMB code is free software;
! You can use it, redistribute it, and/or modify it under the terms
! of the GNU General Public License as published by the Free Software Foundation;
! either version 3 of the License, or (at your option) any later version.
! The full text of the license can be found in the file eftcamb/LICENSE at
! the top level of the EFTCAMB distribution.
!
!----------------------------------------------------------------------------------------

!> @file benchmark.F90
!! This application runs the benchmarks for a given EFTCAMB model.

program benchmarker

  use CAMB
  use omp_lib

  implicit none

  character(len=:), allocatable :: InputFile
  real(dl), allocatable :: timings(:)
  real(dl)              :: time, variance, t1, t2
  Type(TIniFile) :: Ini
  logical bad
  character(LEN=1024) :: ErrMsg, benchmark_buffer
  integer :: benchmark_count, i

  InputFile = ''
  if (GetParamCount() /= 0)  InputFile = GetParam(1)
  if (InputFile == '') error stop 'No parameter input file'

  ! EFTCAMB BENCHMARK: read from command line the number of times that the command is executed
  ! for benchmarking. More times means a more accurate measure of its timing.
  benchmark_count = 10
  benchmark_buffer = ''
  if (GetParamCount() > 1) then
      benchmark_buffer = GetParam(2)
      read(benchmark_buffer, '(I3)') benchmark_count
  end if

  call Ini%Open(InputFile, bad, .false.)
  if (bad) then
      write(*,*) 'File not found: '//trim(InputFile)
      error stop
  end if

  highL_unlensed_cl_template = Ini%Read_String_Default( &
      'highL_unlensed_cl_template', highL_unlensed_cl_template)
  call Ini%Read('number_of_threads', ThreadNum)
  call Ini%Read('DebugParam', DebugParam)
  if (Ini%HasKey('DebugMsgs')) call Ini%Read('DebugMsgs', DebugMsgs)
  Ini%Fail_on_not_found = .false.

  ! no feedback for benchmarker:
  FeedbackLevel = 0

  ! adjust options in case this is the profiler:
#ifdef PROFILE
  benchmark_count = 1
  ThreadNum = 1
#endif

  ! allocate:
  allocate( timings(benchmark_count) )
  ! call consecutively camb to get timings and variance:
  do i=1, benchmark_count
      t1 = omp_get_wtime()
      global_error_flag = 0
      if ( .not. CAMB_RunFromIni(Ini, InputFile, ErrMsg) ) global_error_flag = 1
      t2 = omp_get_wtime()
      t2 = (t2 -t1)
      timings(i) = t2
  end do

  ! compute time statistics:
  ! 1) mean and variance:
  time     = 0._dl
  if ( benchmark_count>0 ) then
      do i=1, benchmark_count
          time = time +timings(i)
      end do
      time = time/real(benchmark_count)
  end if
  ! 2) variance:
  variance = 0._dl
  if ( benchmark_count>1 ) then
      do i=1, benchmark_count
          variance = variance + (timings(i)-time)**2
      end do
      variance = sqrt( variance/real(benchmark_count-1) )
  end if

  ! print results:
  if ( benchmark_count>1 ) then
      write(*,"(F9.3,a,F9.3,a)") time, ' +- ', variance, ' (sec)'
  else if ( benchmark_count>0 .and. benchmark_count<=1 ) then
      write(*,"(F9.3,a)") time, ' (sec)'
  end if

end program benchmarker
