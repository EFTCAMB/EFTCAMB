! Standalone replica of the call pattern in fortran/eftcamb/07f_designer_models/007p8_3p1_wDE_designer.f90:
!   InitBackground  (intent(out)) : success = .true. ; call self%solve_designer_equations(success=success)
!   SolveDesigner   (intent(out)) : only assigns success = .false. on an error path
! The caller (p31_intent_out_main.f90) mimics fortran/results.f90:584 (local logical, deferred binding call).
module p31_mod
    implicit none
    type, abstract :: base_model
    contains
        procedure(init_iface), deferred :: initialize_background
    end type base_model
    abstract interface
        subroutine init_iface(self, flag, success)
            import :: base_model
            class(base_model) :: self
            integer, intent(in) :: flag
            logical, intent(out) :: success
        end subroutine init_iface
    end interface
    type, extends(base_model) :: model3p1
    contains
        procedure :: initialize_background => init_bg
        procedure :: solve_designer_equations => solve_eq
    end type model3p1
contains
    subroutine init_bg(self, flag, success)
        class(model3p1) :: self
        integer, intent(in) :: flag
        logical, intent(out) :: success
        success = .true.
        call self%solve_designer_equations(flag, success=success)
    end subroutine init_bg
    subroutine solve_eq(self, flag, success)
        class(model3p1) :: self
        integer, intent(in) :: flag
        logical, intent(out) :: success
#ifdef FIXED
        success = .true.
#endif
        if (flag < 0) then
            success = .false.
            write(*, *) 'error path'
            return
        end if
    end subroutine solve_eq
end module p31_mod
