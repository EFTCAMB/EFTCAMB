! ----------------------------------------------------------------------------------------
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
! ----------------------------------------------------------------------------------------

!> @file 11_EFTCAMB_python.f90
!! This file contains the code to interface the EFTCAMB code to the python side
!! of CAMB. There is no physics here, just wrappers...


! ----------------------------------------------------------------------------------------
!> This module contains the code to interface the EFTCAMB code to the python side
!! of CAMB. There is no physics here, just wrappers...

!> @author Marco Raveri

    module EFTCAMB_python

    use IniObjects
    use config, only: GlobalError, global_error_flag, error_unsupported_params, error_ini
    use results
    use EFTCAMB_stability
    use MassiveNu

    implicit none

    !> Ini file used on the python path. A value that cannot be read, or a missing key without default, records an
    !! error instead of stopping the program (TIniFile ends with error stop, which kills the python interpreter):
    !! the read then returns the default (or zero) and initialize_EFTCAMB_from_py reports the first error.
    type, extends(TIniFile) :: TEFTCAMB_python_ini
        logical :: failed = .false.                       !< whether a read failed since the last reset
        character(len=:), allocatable :: error_message    !< message of the first read error
    contains
    procedure :: Error => EFTCAMB_python_ini_error
    procedure :: FailStop => EFTCAMB_python_ini_failstop
    procedure :: Read_Int => EFTCAMB_python_ini_read_int
    procedure :: Read_Double => EFTCAMB_python_ini_read_double
    procedure :: Read_Logical => EFTCAMB_python_ini_read_logical
    procedure :: reset_errors => EFTCAMB_python_ini_reset_errors
    end type TEFTCAMB_python_ini

    type(TEFTCAMB_python_ini), private :: Ini

    contains

    ! ---------------------------------------------------------------------------------------------
    !> Records an ini read error (first message kept) instead of stopping.
    subroutine EFTCAMB_python_ini_error(this, Msg, Key)

    implicit none

    class(TEFTCAMB_python_ini) :: this
    character(len=*) :: Msg
    character(len=*), optional :: Key
    character(len=:), allocatable :: message

    if (present(Key)) then
        message = 'Error for key "' // trim(Key) // '" : ' // Msg
    else
        message = 'Error : ' // Msg
    end if
    write(*, *) message
    if (.not. this%failed) this%error_message = message
    this%failed = .true.

    end subroutine EFTCAMB_python_ini_error

    ! ---------------------------------------------------------------------------------------------
    !> Records a failure instead of stopping.
    subroutine EFTCAMB_python_ini_failstop(this)

    implicit none

    class(TEFTCAMB_python_ini) :: this

    if (.not. this%failed) this%error_message = 'Error reading the EFTCAMB parameters'
    this%failed = .true.

    end subroutine EFTCAMB_python_ini_failstop

    ! ---------------------------------------------------------------------------------------------
    !> Clears the recorded read errors.
    subroutine EFTCAMB_python_ini_reset_errors(this)

    implicit none

    class(TEFTCAMB_python_ini) :: this

    this%failed = .false.
    this%error_message = ''

    end subroutine EFTCAMB_python_ini_reset_errors

    ! ---------------------------------------------------------------------------------------------
    !> As TIniFile%Read_Int; on a read error returns the default (or 0).
    function EFTCAMB_python_ini_read_int(this, Key, Default, min, max, OK) result(value)

    implicit none

    class(TEFTCAMB_python_ini) :: this
    integer :: value
    integer, intent(in), optional :: Default, min, max
    logical, intent(out), optional :: OK
    character(len=*), intent(in) :: Key
    character(len=:), pointer :: S
    integer :: status

    value = 0
    if (present(Default)) value = Default
    if (present(OK)) OK = .true.
    S => this%Read_String(Key, .not. present(Default))
    if (S == '') then
        if (present(Default)) then
            call this%ReadValues%Add(Key, Default)
        else
            call this%Error('missing key', Key)
        end if
        return
    end if
    status = 1
    if (verify(trim(S), '-+0123456789') == 0) read(S, *, iostat=status) value
    if (status /= 0) then
        value = 0
        if (present(Default)) value = Default
        if (present(OK)) then
            OK = .false.
            return
        end if
        call this%Error('error reading integer', Key)
        return
    end if
    if (present(max)) then
        if (value > max) call this%Error('value > max', Key)
    end if
    if (present(min)) then
        if (value < min) call this%Error('value < min', Key)
    end if

    end function EFTCAMB_python_ini_read_int

    ! ---------------------------------------------------------------------------------------------
    !> As TIniFile%Read_Double; on a read error returns the default (or 0).
    function EFTCAMB_python_ini_read_double(this, Key, Default, min, max) result(value)

    implicit none

    class(TEFTCAMB_python_ini) :: this
    double precision :: value
    double precision, intent(in), optional :: Default, min, max
    character(len=*), intent(in) :: Key
    character(len=:), pointer :: S
    integer :: status

    value = 0.d0
    if (present(Default)) value = Default
    S => this%Read_String(Key, .not. present(Default))
    if (S == '') then
        if (present(Default)) then
            call this%ReadValues%Add(Key, Default)
        else
            call this%Error('missing key', Key)
        end if
        return
    end if
    read(S, *, iostat=status) value
    if (status /= 0) then
        value = 0.d0
        if (present(Default)) value = Default
        call this%Error('error reading double', Key)
        return
    end if
    if (present(max)) then
        if (value > max) call this%Error('value > max', Key)
    end if
    if (present(min)) then
        if (value < min) call this%Error('value < min', Key)
    end if

    end function EFTCAMB_python_ini_read_double

    ! ---------------------------------------------------------------------------------------------
    !> As TIniFile%Read_Logical; on a read error returns the default (or false).
    function EFTCAMB_python_ini_read_logical(this, Key, Default) result(value)

    implicit none

    class(TEFTCAMB_python_ini) :: this
    logical :: value
    logical, intent(in), optional :: Default
    character(len=*), intent(in) :: Key
    character(len=:), pointer :: S
    integer :: status

    value = .false.
    if (present(Default)) value = Default
    S => this%Read_String(Key, .not. present(Default))
    if (S == '') then
        if (present(Default)) then
            call this%ReadValues%Add(Key, Default)
        else
            call this%Error('missing key', Key)
        end if
        return
    end if
    status = 1
    if (verify(trim(S), '10TF') == 0) read(S, *, iostat=status) value
    if (status /= 0) then
        value = .false.
        if (present(Default)) value = Default
        call this%Error('error reading logical', Key)
    end if

    end function EFTCAMB_python_ini_read_logical

    ! ---------------------------------------------------------------------------------------------
    !> Helper function to strip all whitespaces
    function helper_stripspaces(string)
    character(len=*) :: string
    character(len=:), allocatable :: helper_stripspaces
    integer :: stringLen
    integer :: last, actual
    stringLen = len (string)
    last = 1
    actual = 1
    do while (actual < stringLen)
        if (string(last:last) == ' ') then
            actual = actual + 1
            string(last:last) = string(actual:actual)
            string(actual:actual) = ' '
        else
            last = last + 1
            if (actual < last) &
                actual = last
        end if
    end do
    helper_stripspaces = trim(string)
    end function

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine to initialize EFTCAMB from python
    subroutine initialize_EFTCAMB_from_py_read_params(key, nkey, value, nvalue, mode) bind(c, &
        name='initialize_EFTCAMB_from_py_read_params_' )

    use iso_c_binding

    implicit none

    ! interface:
    integer(c_int), value, intent(in) :: nkey, nvalue, mode
    character(len=1), intent(in) :: key(nkey), value(nvalue)
    ! local fortran copies:
    character(len=nkey) :: f_key
    character(len=nvalue) :: f_value
    integer :: i

    ! initialize local copies:
    do i = 1, nkey
        f_key(i:i) = key(i)(1:1)
    end do
    do i = 1, nvalue
        f_value(i:i) = value(i)(1:1)
    end do

    ! if mode == 0 then create the INI (virtual) file
    if (mode == 0) then
        call Ini%Close()
        call Ini%Init(ignoreDuplicates=.true.)
        call Ini%reset_errors()
        ! if mode == 1 append settings to the INI file
    else if (mode == 1) then
        call Ini%AddString(trim(helper_stripspaces(f_key)), trim(helper_stripspaces(f_value)))
    end if

    end subroutine initialize_EFTCAMB_from_py_read_params

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine to initialize EFTCAMB from python
    subroutine initialize_EFTCAMB_from_py(P, eft_error)

    implicit none

    type(CAMBparams) :: P          !< the parameters (python passes a CAMBparams)
    integer          :: eft_error
    logical        :: reset_in_place
    type(TEFTCAMB) :: default_EFTCAMB
    type(TEFTCAMB_parameter_cache) :: default_eft_par_cache

    ! allocate, or reset to the default values in place. Keeping the objects at the same address
    ! keeps valid the python handles to them (e.g. eftcamb = pars.EFTCAMB) across initializations.
    reset_in_place = .false.
    if (allocated(P%EFTCAMB)) then
        select type (eft => P%EFTCAMB)
        type is (TEFTCAMB)
            eft = default_EFTCAMB
            reset_in_place = .true.
        end select
    end if
    if (.not. reset_in_place) then
        if (allocated(P%EFTCAMB)) deallocate(P%EFTCAMB)
        allocate(P%EFTCAMB)
    end if
    reset_in_place = .false.
    if (allocated(P%eft_par_cache)) then
        select type (cache => P%eft_par_cache)
        type is (TEFTCAMB_parameter_cache)
            cache = default_eft_par_cache
            reset_in_place = .true.
        end select
    end if
    if (.not. reset_in_place) then
        if (allocated(P%eft_par_cache)) deallocate(P%eft_par_cache)
        allocate(P%eft_par_cache)
    end if
    ! read the options in:
    eft_error = 0
    global_error_flag = 0
    call P%EFTCAMB%EFTCAMB_init_from_file(Ini)
    if (global_error_flag /= 0) eft_error = 1
    ! call the eftcamb allocation:
    if (P%EFTCAMB%EFTflag /= 0 .and. eft_error == 0 .and. .not. Ini%failed) then
        ! initialize the model from file:
        call P%EFTCAMB%EFTCAMB_init_model_from_file(Ini, eft_error)
        if (global_error_flag /= 0 .and. eft_error == 0) eft_error = 1
    end if
    ! a value that could not be read:
    if (Ini%failed) then
        eft_error = 2
        call GlobalError('EFTCAMB: ' // Ini%error_message, error_ini)
    end if
    ! on failure go back to a consistent GR state (no model), the caller raises the error:
    if (eft_error > 0) then
        P%EFTCAMB%EFTflag = 0
        if (allocated(P%EFTCAMB%model)) deallocate(P%EFTCAMB%model)
        ! keep the message of an error raised during the initialization:
        if (global_error_flag == 0) &
            call GlobalError('EFTCAMB: model initialization failed (invalid or unsupported model parameters, '// &
            'see the EFTCAMB output)', error_unsupported_params)
    end if
    ! store the parameters that have been read in the EFTCAMB object:
    call P%EFTCAMB%EFTCAMB_store_read_parameters(Ini)

    end subroutine initialize_EFTCAMB_from_py

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine to get the keys and values of the parameters read when an EFTCAMB object was initialized
    subroutine get_read_parameters(EFTCAMB, key, value, index, mode)

    implicit none

    type(TEFTCAMB) :: EFTCAMB
    character(len=200) :: key, value
    !< if mode=0 this returns the number of parameters, otherwise this is used as the parameter number to return
    integer :: index
    !< if mode=0 returns the number of parameters only, if mode=1 then returns the key, value combination
    integer :: mode

    ! if mode == 0 return the number of elements
    if (mode == 0) then
        index = 0
        if (allocated(EFTCAMB%read_param_keys)) index = size(EFTCAMB%read_param_keys)
        ! if mode == 1 return the key, value combination
    else if (mode == 1) then
        key = ''
        value = ''
        if (.not. allocated(EFTCAMB%read_param_keys)) return
        if (index < 1 .or. index > size(EFTCAMB%read_param_keys)) return
        key = EFTCAMB%read_param_keys(index)
        value = EFTCAMB%read_param_values(index)
    end if

    end subroutine get_read_parameters

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine to print EFTCAMB feedback:
    subroutine EFTCAMB_feedback(EFTCAMB, print_params)

    implicit none

    type(TEFTCAMB) :: EFTCAMB
    logical        :: print_params

    ! print the EFTCAMB header:
    call EFTCAMB%EFTCAMB_print_header()
    ! print model feedback:
    call EFTCAMB%EFTCAMB_print_model_feedback(print_params)

    end subroutine EFTCAMB_feedback

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine to get number of model parameters:
    subroutine EFTCAMB_get_num_params(EFTCAMB, num_params)

    implicit none

    type(TEFTCAMB) :: EFTCAMB
    integer        :: num_params

    num_params = 0
    if (allocated(EFTCAMB%model)) num_params = EFTCAMB%model%parameter_number

    end subroutine EFTCAMB_get_num_params

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine to get parameter names:
    subroutine EFTCAMB_get_param_names(EFTCAMB, i, name)

    implicit none

    type(TEFTCAMB)     :: EFTCAMB
    integer            :: i
    character(len=200) :: name

    name = ''
    if (allocated(EFTCAMB%model)) call EFTCAMB%model%parameter_names(i, name)

    end subroutine EFTCAMB_get_param_names

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine to get parameter labels:
    subroutine EFTCAMB_get_param_labels(EFTCAMB, i, name)

    implicit none

    type(TEFTCAMB)     :: EFTCAMB
    integer            :: i
    character(len=200) :: name

    name = ''
    if (allocated(EFTCAMB%model)) call EFTCAMB%model%parameter_names_latex(i, name)

    end subroutine EFTCAMB_get_param_labels

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine to get parameter values:
    subroutine EFTCAMB_get_param_values(EFTCAMB, i, value)

    implicit none

    type(TEFTCAMB) :: EFTCAMB
    integer        :: i
    real(dl)       :: value

    value = 0._dl
    if (allocated(EFTCAMB%model)) call EFTCAMB%model%parameter_values(i, value)

    end subroutine EFTCAMB_get_param_values

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine to get the EFTCAMB model name:
    subroutine EFTCAMB_get_model_name(EFTCAMB, model_name)

    implicit none

    type(TEFTCAMB)     :: EFTCAMB
    character(len=200) :: model_name

    if (allocated(EFTCAMB%model)) then
        model_name = EFTCAMB%model%name
    else
        model_name = ''
    end if

    end subroutine EFTCAMB_get_model_name

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine to set the EFTCAMB model name:
    subroutine EFTCAMB_set_model_name(EFTCAMB, model_name)

    implicit none

    type(TEFTCAMB)     :: EFTCAMB
    character(len=200) :: model_name

    if (allocated(EFTCAMB%model)) then
        EFTCAMB%model%name = trim(adjustl(model_name))
    end if

    end subroutine EFTCAMB_set_model_name

    ! ---------------------------------------------------------------------------------------------
    ! C-interoperable string getters/setters. The routines above take character(len=200) dummies and
    ! rely on the hidden length argument; these take every buffer length explicitly as a trailing
    ! size_t passed by value. Output buffers are filled with the (truncated) string, blank padded.
    ! The first argument is the address of the TEFTCAMB object (ctypes byref). C signatures:
    !   eftcamb_get_read_parameters_c(eft*, char* key, char* value, int* index, int* mode, size_t nkey, size_t nvalue)
    !   eftcamb_get_param_names_c(eft*, int* i, char* name, size_t n)
    !   eftcamb_get_param_labels_c(eft*, int* i, char* name, size_t n)
    !   eftcamb_get_model_name_c(eft*, char* name, size_t n)
    !   eftcamb_set_model_name_c(eft*, char* name, size_t n)

    !> Helper: copy a Fortran string into a C character buffer of length n, blank padded.
    subroutine EFTCAMB_python_to_c_buffer(string, buffer, n)

    use iso_c_binding, only: c_char, c_size_t

    implicit none

    character(len=*), intent(in)          :: string
    integer(c_size_t), intent(in)         :: n
    character(kind=c_char), intent(inout) :: buffer(n)
    integer(c_size_t) :: i

    do i = 1, n
        if (i <= len(string)) then
            buffer(i) = string(i:i)
        else
            buffer(i) = ' '
        end if
    end do

    end subroutine EFTCAMB_python_to_c_buffer

    !> Get the key and value of a parameter read at initialization (mode=1), or their number (mode=0, in index).
    subroutine EFTCAMB_get_read_parameters_c(eftcamb_ptr, key, value, index, mode, nkey, nvalue) bind(c, &
        name='eftcamb_get_read_parameters_c')

    use iso_c_binding

    implicit none

    type(c_ptr), value                    :: eftcamb_ptr
    integer(c_size_t), value              :: nkey, nvalue
    character(kind=c_char), intent(inout) :: key(nkey), value(nvalue)
    integer(c_int), intent(inout)         :: index
    integer(c_int), intent(in)            :: mode
    type(TEFTCAMB), pointer :: EFTCAMB
    character(len=200) :: f_key, f_value
    integer :: f_index, f_mode

    call c_f_pointer(eftcamb_ptr, EFTCAMB)
    f_index = index
    f_mode = mode
    call get_read_parameters(EFTCAMB, f_key, f_value, f_index, f_mode)
    if (mode == 0) then
        index = f_index
    else
        call EFTCAMB_python_to_c_buffer(f_key, key, nkey)
        call EFTCAMB_python_to_c_buffer(f_value, value, nvalue)
    end if

    end subroutine EFTCAMB_get_read_parameters_c

    !> Get the name of the i-th model parameter.
    subroutine EFTCAMB_get_param_names_c(eftcamb_ptr, i, name, n) bind(c, name='eftcamb_get_param_names_c')

    use iso_c_binding

    implicit none

    type(c_ptr), value                    :: eftcamb_ptr
    integer(c_int), intent(in)            :: i
    integer(c_size_t), value              :: n
    character(kind=c_char), intent(inout) :: name(n)
    type(TEFTCAMB), pointer :: EFTCAMB
    character(len=200) :: f_name

    call c_f_pointer(eftcamb_ptr, EFTCAMB)
    call EFTCAMB_get_param_names(EFTCAMB, int(i), f_name)
    call EFTCAMB_python_to_c_buffer(f_name, name, n)

    end subroutine EFTCAMB_get_param_names_c

    !> Get the latex label of the i-th model parameter.
    subroutine EFTCAMB_get_param_labels_c(eftcamb_ptr, i, name, n) bind(c, name='eftcamb_get_param_labels_c')

    use iso_c_binding

    implicit none

    type(c_ptr), value                    :: eftcamb_ptr
    integer(c_int), intent(in)            :: i
    integer(c_size_t), value              :: n
    character(kind=c_char), intent(inout) :: name(n)
    type(TEFTCAMB), pointer :: EFTCAMB
    character(len=200) :: f_name

    call c_f_pointer(eftcamb_ptr, EFTCAMB)
    call EFTCAMB_get_param_labels(EFTCAMB, int(i), f_name)
    call EFTCAMB_python_to_c_buffer(f_name, name, n)

    end subroutine EFTCAMB_get_param_labels_c

    !> Get the model name (empty without a model).
    subroutine EFTCAMB_get_model_name_c(eftcamb_ptr, name, n) bind(c, name='eftcamb_get_model_name_c')

    use iso_c_binding

    implicit none

    type(c_ptr), value                    :: eftcamb_ptr
    integer(c_size_t), value              :: n
    character(kind=c_char), intent(inout) :: name(n)
    type(TEFTCAMB), pointer :: EFTCAMB

    call c_f_pointer(eftcamb_ptr, EFTCAMB)
    if (allocated(EFTCAMB%model)) then
        call EFTCAMB_python_to_c_buffer(EFTCAMB%model%name, name, n)
    else
        call EFTCAMB_python_to_c_buffer('', name, n)
    end if

    end subroutine EFTCAMB_get_model_name_c

    !> Set the model name from the first n characters of name (no-op without a model).
    subroutine EFTCAMB_set_model_name_c(eftcamb_ptr, name, n) bind(c, name='eftcamb_set_model_name_c')

    use iso_c_binding

    implicit none

    type(c_ptr), value                 :: eftcamb_ptr
    integer(c_size_t), value           :: n
    character(kind=c_char), intent(in) :: name(n)
    type(TEFTCAMB), pointer :: EFTCAMB
    character(len=max(int(n), 0)) :: f_name
    integer :: i

    call c_f_pointer(eftcamb_ptr, EFTCAMB)
    do i = 1, int(n)
        f_name(i:i) = name(i)
    end do
    if (allocated(EFTCAMB%model)) EFTCAMB%model%name = trim(adjustl(f_name))

    end subroutine EFTCAMB_set_model_name_c

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine to get EFTCAMB effective w0wa values:
    subroutine EFTCAMB_get_effective_w0wa_values(EFTCAMB, w0, wa)

    implicit none

    type(TEFTCAMB) :: EFTCAMB
    real(dl)       :: w0        !< the effective w0 value (output)
    real(dl)       :: wa        !< the effective wa value (output)

    w0 = -1._dl
    wa = 0._dl
    if (EFTCAMB%EFTflag /= 0 .and. allocated(EFTCAMB%model)) call EFTCAMB%Effective_w_wa(w0, wa)

    end subroutine EFTCAMB_get_effective_w0wa_values

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine to get the EFT time evolution:
    subroutine EFTCAMB_GetEvolution(this, nq, q, neta, eta, timestep_cache)

    use GaugeInterface
    use CAMBmain

    type(CAMBdata), target :: this
    integer, intent(in)  :: nq, neta
    real(dl), intent(in) :: q(nq), eta(neta)
    type(TEFTCAMB_timestep_cache) :: timestep_cache(neta*nq)

    real(dl) :: taustart
    integer  :: q_ix
    logical  :: want_scalars
    type(EvolutionVars) :: Ev
    ! activate status:
    call SetActiveState(this)
    ! initialize:
    global_error_flag = 0
    ! check that these results carry an initialized EFTCAMB model and background:
    if (.not. EFTCAMB_python_results_ready(this, 'get_evolution', require_eft=.false.)) return
    ! get initial time and call Thermo initialization if needed:
    taustart = min(minval(eta), GetTauStart(maxval(q)))
    ! the thermal history must start before the EFT perturbations turn on (not below tau = 1e-12 Mpc, as in
    ! cmbmain InitVars; an earlier turn-on is refused by initial()):
    if (this%CP%EFTCAMB%EFTflag /= 0) then
        taustart = min(taustart, max(this%DeltaTime(0._dl, this%CP%EFTCAMB%EFTCAMB_pert_turn_on), 1.e-12_dl))
    end if
    if (.not. this%ThermoData%HasThermoData .or. taustart < this%ThermoData%tauminn) call this%ThermoData%Init(this, &
        taustart)
    ! GetNumEqns sets up the scalar hierarchy only when scalars are wanted (not after a transfer-only run):
    want_scalars = this%CP%WantScalars
    this%CP%WantScalars = .true.
    ! do the parallel loop over the requested k modes.
    !$OMP PARALLEL DO DEFAULT(SHARED), SCHEDULE(DYNAMIC), PRIVATE(EV, q_ix)
    do q_ix = 1, nq
        if (global_error_flag == 0) then
            ! reset the OpenMP-private evolution variables so that every k mode starts from the same state:
            call Ev%eft_cache%initialize()
            Ev%EFTCAMB_printing = .false.
            Ev%EFTCAMB_QS = .false.
            ! initialize evolution variables:
            Ev%q_ix = q_ix
            Ev%q = q(q_ix)
            Ev%TransferOnly = .false.
            Ev%q2 = Ev%q**2
            Ev%ThermoData => this%ThermoData
            call GetNumEqns(Ev)
            ! call utility function:
            call GetEFTOutputEvolutionFork(State, Ev, neta, eta, timestep_cache((q_ix - 1)*neta + 1:q_ix*neta))
        end if
    end do
    !$OMP END PARALLEL DO
    this%CP%WantScalars = want_scalars

    end subroutine EFTCAMB_GetEvolution

    ! ---------------------------------------------------------------------------------------------
    !> Function that checks that a CAMBdata carries an initialized EFTCAMB model and parameter cache.
    !! Otherwise it sets a global error and returns false.
    logical function EFTCAMB_python_results_ready(State, caller, require_eft)

    implicit none

    type(CAMBdata)               :: State
    character(len=*), intent(in) :: caller       !< name of the calling routine, for the error message
    logical, intent(in)          :: require_eft  !< if true EFTflag = 0 is an error

    EFTCAMB_python_results_ready = .false.
    ! nested ifs: .and. does not short-circuit.
    if (.not. allocated(State%CP%EFTCAMB)) then
        call GlobalError(caller//': EFTCAMB is not initialized for these results', error_unsupported_params)
        return
    end if
    if (State%CP%EFTCAMB%EFTflag == 0) then
        if (require_eft) then
            call GlobalError(caller//': EFTflag = 0, there is no EFTCAMB model', error_unsupported_params)
            return
        end if
    else
        if (.not. allocated(State%CP%EFTCAMB%model)) then
            call GlobalError(caller//': the EFTCAMB model is not initialized', error_unsupported_params)
            return
        end if
        if (.not. allocated(State%CP%eft_par_cache)) then
            call GlobalError(caller//': the EFTCAMB parameter cache is not initialized', error_unsupported_params)
            return
        end if
        if (.not. allocated(State%CP%eft_par_cache%grhormass)) then
            call GlobalError(caller//': nothing was computed for these results, call get_background or '// &
                'get_results with some output first', error_unsupported_params)
            return
        end if
    end if
    EFTCAMB_python_results_ready = .true.

    end function EFTCAMB_python_results_ready

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine to get the EFT time evolution k mode by kmode:
    subroutine GetEFTOutputEvolutionFork(this, EV, neta, eta, timestep_cache)

    use GaugeInterface
    use CAMBmain
    use RungeKuttaDP45Module, only: RungeKuttaDP45Settings

    type(CAMBdata) :: this
    type(EvolutionVars) :: EV
    integer, intent(in) :: neta
    real(dl), intent(in) :: eta(neta)
    type(TEFTCAMB_timestep_cache) :: timestep_cache(neta)

    type(RungeKuttaDP45Settings) :: rk_settings
    real(dl) w(EV%nvar, 9), y(EV%nvar), yprime(EV%nvar)
    real(dl) tau, tol1, tauend, taustart
    integer j, ind

    ! intialize:
    w = 0
    y = 0
    taustart = min(GetTauStart(min(500._dl, EV%q)), minval(eta))
    ! call initial conditions:
    call initial(EV, y, taustart)
    ! do the time loop:
    tau = taustart
    ind = 1
    ! tighter than the source integration: the small, fast oscillating pi sector is otherwise dominated by
    ! integration error (IntTolBoost = 8 is converged; much smaller tolerances make RK45 fail)
    tol1 = base_tol/exp(max(CP%Accuracy%AccuracyBoost*CP%Accuracy%IntTolBoost, 8._dl) - 1)
    do j = 1, neta
        tauend = eta(j)
        if (tauend < taustart) cycle
        ! evolve the preturbations:
        EV%EFTCAMB_printing = .false.
        call GaugeInterface_EvolveScal(EV, tau, y, tauend, tol1, ind, rk_settings, w)
        yprime = 0
        EV%EFTCAMB_printing = .true.
        call EV%eft_cache%initialize()
        call derivs(EV, EV%ScalEqsToPropagate, tau, y, yprime)
        ! now EV%eft_cache is filled with all quantities at the required time:
        timestep_cache(j) = EV%eft_cache
        if (global_error_flag /= 0) return
    end do

    end subroutine GetEFTOutputEvolutionFork

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine to get the EFT functions at scale factor a (no recomputation):
    subroutine EFTCAMB_GetEFTFunctions(State, na, as, etas, timestep_cache)
    ! Background quantities and EFT functions at the requested scale factors. The calculation follows
    ! the background part of derivs in equations.f90, so that the results are the same as get_evolution.
    type(CAMBdata)                :: State
    integer, intent(in)           :: na
    real(dl), intent(in)          :: as(na)
    real(dl), intent(in)          :: etas(na)
    type(TEFTCAMB_timestep_cache) :: timestep_cache(na)

    integer i, nu_i
    real(dl) :: grhonu, gpnu, grhonu_i, gpnu_i, grhormass_t, gpnudot_i, am
    real(dl) :: a, adotoa, grhorad

    global_error_flag = 0
    ! check that these results carry an initialized EFTCAMB model and background:
    if (.not. EFTCAMB_python_results_ready(State, 'get_eft_functions', require_eft=.true.)) return

    associate(EFTCAMB => State%CP%EFTCAMB, par_cache => State%CP%eft_par_cache)
        do i = 1, na
            associate(cache => timestep_cache(i))
                a = as(i)
                cache%a = a
                cache%tau = etas(i)
                ! background densities and pressure of the matter species:
                cache%grhob_t = par_cache%grhob/a
                cache%grhoc_t = par_cache%grhoc/a
                cache%grhor_t = par_cache%grhornomass/a/a
                cache%grhog_t = par_cache%grhog/a/a
                grhorad = cache%grhog_t + cache%grhor_t
                grhonu = 0._dl
                gpnu = 0._dl
                if (par_cache%Num_Nu_Massive /= 0) then
                    do nu_i = 1, par_cache%Nu_mass_eigenstates
                        grhormass_t = par_cache%grhormass(nu_i)/a/a
                        call ThermalNuBackground%rho_P(a*par_cache%nu_masses(nu_i), grhonu_i, gpnu_i)
                        grhonu = grhonu + grhormass_t*grhonu_i
                        gpnu = gpnu + grhormass_t*gpnu_i
                    end do
                end if
                cache%grhom_t = cache%grhob_t + cache%grhoc_t + grhorad + grhonu
                cache%grhoa2 = cache%grhom_t*a**2
                cache%gpresm_t = gpnu + grhorad/3._dl
                ! expansion history:
                if (.not. EFTCAMB%EFTCAMB_model_is_designer) then
                    call EFTCAMB%model%compute_background_EFT_functions(a, par_cache, cache)
                    call EFTCAMB%model%compute_adotoa(a, par_cache, cache)
                else
                    call EFTCAMB%model%compute_adotoa(a, par_cache, cache)
                end if
                adotoa = cache%adotoa
                ! massive neutrinos and first derivative of the matter pressure:
                cache%grhonu_tot = 0._dl
                cache%gpinu_tot = 0._dl
                cache%grhonudot_tot = 0._dl
                cache%gpinudot_tot = 0._dl
                cache%gpinudotdot_tot = 0._dl
                if (par_cache%Num_Nu_Massive /= 0) then
                    do nu_i = 1, par_cache%Nu_mass_eigenstates
                        grhormass_t = par_cache%grhormass(nu_i)/a/a
                        am = a*par_cache%nu_masses(nu_i)
                        call ThermalNuBackground%rho_P(am, grhonu_i, gpnu_i)
                        cache%grhonu_tot = cache%grhonu_tot + grhormass_t*grhonu_i
                        cache%gpinu_tot = cache%gpinu_tot + grhormass_t*gpnu_i
                        cache%gpinudot_tot = cache%gpinudot_tot + &
                            grhormass_t*(ThermalNuBackground%pidot(am, adotoa, gpnu_i) - 4._dl*adotoa*gpnu_i)
                        cache%grhonudot_tot = cache%grhonudot_tot + &
                            grhormass_t*(ThermalNuBackground%drho(am, adotoa) - 4._dl*adotoa*grhonu_i)
                    end do
                end if
                cache%gpresdotm_t = -4._dl*adotoa*grhorad/3._dl + cache%gpinudot_tot
                call EFTCAMB%model%compute_H_derivs(a, par_cache, cache)
                ! second derivative of the massive neutrino and matter pressure (needs Hdot):
                if (par_cache%Num_Nu_Massive /= 0) then
                    do nu_i = 1, par_cache%Nu_mass_eigenstates
                        grhormass_t = par_cache%grhormass(nu_i)/a/a
                        am = a*par_cache%nu_masses(nu_i)
                        call ThermalNuBackground%rho_P(am, grhonu_i, gpnu_i)
                        gpnudot_i = ThermalNuBackground%pidot(am, adotoa, gpnu_i)
                        cache%gpinudotdot_tot = cache%gpinudotdot_tot - &
                            4._dl*adotoa*grhormass_t*(gpnudot_i - 4._dl*adotoa*gpnu_i) + &
                            grhormass_t*(ThermalNuBackground%pidotdot(am, adotoa, cache%Hdot, gpnu_i, gpnudot_i) &
                            - 4._dl*cache%Hdot*gpnu_i - 4._dl*adotoa*gpnudot_i)
                    end do
                end if
                cache%gpresdotdotm_t = -(4._dl/3._dl)*cache%Hdot*grhorad + (16._dl/3._dl)*grhorad*adotoa**2 + &
                    cache%gpinudotdot_tot
                ! third derivative of H:
                call EFTCAMB%model%compute_H_derivs(a, par_cache, cache)
                ! remaining background quantities:
                ! (again for full models, now with the current adotoa, Hdot, Hdotdot: their background EFT
                ! functions may use them)
                call EFTCAMB%model%compute_background_EFT_functions(a, par_cache, cache)
                call EFTCAMB%model%compute_rhoQPQ(a, par_cache, cache)
                call EFTCAMB%model%compute_secondorder_EFT_functions(a, par_cache, cache)
                ! higher derivatives (EFTLambdadotdot, EFTcdotdot, EFTOmegaPPPP, Gamma PP/PPP, ...) that the model
                ! does not compute analytically, then the Hubble derivatives again: Hdotdotdot needs them (as in
                ! the stability calculation):
                call EFTCAMB%model%compute_additional_derivs(a, par_cache, cache)
                call EFTCAMB%model%compute_H_derivs(a, par_cache, cache)
            end associate
        end do
    end associate

    end subroutine EFTCAMB_GetEFTFunctions

    ! ---------------------------------------------------------------------------------------------

    end module EFTCAMB_python

! ----------------------------------------------------------------------------------------
