! Minimal analogue of InitBackground -> SolveDesignerEquations success handling (intent(out), set only on error)
module m
contains
  subroutine solve(success, fail)
    logical, intent(out) :: success
    logical, intent(in) :: fail
    if (fail) then
      success = .false.
      return
    end if
  end subroutine
  subroutine init(success, fail)
    logical, intent(inout) :: success
    logical, intent(in) :: fail
    success = .true.
    call solve(success=success, fail=fail)
  end subroutine
end module
program p
  use m
  logical :: s
  logical :: f
  integer :: i
  f = .false.
  s = transfer(-999, s)
  call init(s, f)
  print *, 'success after init (expected T): ', s, ' raw int:', transfer(s, i)
end program
