! Caller for p31_intent_out_mod.f90, mimicking fortran/results.f90 (local logical success, deferred-binding call).
module p31_caller
    use p31_mod
    implicit none
contains
    subroutine poison_stack()
        logical :: junk(256)
        junk = .false.
        call sink(junk)
    end subroutine poison_stack
    subroutine sink(j)
        logical :: j(:)
        if (j(1)) write(*, *) 'never'
    end subroutine sink
    subroutine do_background(m, flag, res)
        class(base_model), intent(inout) :: m
        integer, intent(in) :: flag
        logical, intent(out) :: res
        logical :: success
        call m%initialize_background(flag, success)
        res = success
    end subroutine do_background
end module p31_caller

program p31_test
    use p31_mod
    use p31_caller
    implicit none
    class(base_model), allocatable :: m
    logical :: res
    allocate(model3p1 :: m)
    call poison_stack()
    call do_background(m, 1, res)
    write(*, '(a,l1)') 'success returned on the success path (expected T): ', res
end program p31_test
