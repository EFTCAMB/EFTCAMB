! Track t4a-param: set_param_names / parameter_names / parameter_names_latex / init_from_file with custom names.
program set_names_check
use precision
use IniObjects
use EFTCAMB_parametrizations_1D
implicit none
class(parametrized_function_1D), allocatable :: f
type(TIniFile) :: Ini
integer :: err, i
real(dl) :: v
character(len=200) :: nm, lt
open(10, file='sn.ini', status='replace')
write(10, '(a)') 'EFTw0 = -0.85'
write(10, '(a)') 'EFTwa = 0.33'
write(10, '(a)') 'EFTwn = 2.5'
close(10)
call Ini%Open('sn.ini')
allocate(JBP_parametrization_1D :: f)
call f%set_param_names(['EFTw0', 'EFTwa', 'EFTwn'], ['w_0', 'w_a', 'n  '])
call f%set_name('EFTw', 'w')
err = 0
call f%init_from_file(Ini, err)
do i = 1, f%parameter_number
    call f%parameter_names(i, nm)
    call f%parameter_names_latex(i, lt)
    call f%parameter_value(i, v)
    write(*, '(a,i0,1x,a,1x,a,1x,es14.6)') 'P', i, trim(nm), trim(lt), v
end do
write(*, '(a,es24.16)') 'VALUE0.5 ', f%value(0.5_dl)
end program
