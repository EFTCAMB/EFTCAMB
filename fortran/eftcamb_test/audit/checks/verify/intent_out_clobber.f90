! Minimal model of EFTCAMBDesigner3p1InitBackground -> SolveDesignerEquations success handling.
module m_clobber
    implicit none
contains
    subroutine inner(success, fail)
        logical, intent(out) :: success
        logical, intent(in)  :: fail
        if (fail) then
            success = .false.
            return
        end if
    end subroutine inner

    subroutine outer(success, fail)
        logical, intent(out) :: success
        logical, intent(in)  :: fail
        success = .true.
        call inner(success, fail)
    end subroutine outer
end module m_clobber

program t
    use m_clobber
    implicit none
    logical :: s
    logical :: f
    f = .false.
    s = .false.          ! caller garbage modelled as .false. (results.f90 local success is uninitialised)
    call outer(s, f)
    write(*, '(a,l1,a)') 'success after outer (no failure) = ', s, ' (expected T)'
end program t
