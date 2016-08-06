!     Code for Anisotropies in the Microwave Background
!     by Antony Lewis (http://cosmologist.info/) and Anthony Challinor
!     See readme.html for documentation. This is a sample driver routine that reads
!     in one set of parameters and produdes the corresponding output.

program driver
    use IniFile
    use CAMB
    use LambdaGeneral
    use Lensing
    use AMLUtils
    use Transfer
    use constants
    use Bispectrum
    use CAMBmain
    use omp_lib
#ifdef NAGF95
    use F90_UNIX
#endif

    ! EFTCAMB MOD START
    use EFTdef
    use compile_time_eft
    ! EFTCAMB MOD END

    implicit none

    Type(CAMBparams) P

    character(LEN=Ini_max_string_len) numstr, VectorFileName, &
        InputFile, ScalarFileName, TensorFileName, TotalFileName, LensedFileName,&
        LensedTotFileName, LensPotentialFileName,ScalarCovFileName
    integer i
    character(LEN=Ini_max_string_len) TransferFileNames(max_transfer_redshifts), &
        MatterPowerFileNames(max_transfer_redshifts), outroot, version_check
    real(dl) output_factor, nmassive
    real(dl) t1, t2

#ifdef WRITE_FITS
    character(LEN=Ini_max_string_len) FITSfilename
#endif

    logical bad

    InputFile = ''
    if (GetParamCount() /= 0)  InputFile = GetParam(1)
    if (InputFile == '') stop 'No parameter input file'

    call Ini_Open(InputFile, 1, bad, .false.)
    if (bad) stop 'Error opening parameter file'

    Ini_fail_on_not_found = .false.

    outroot = Ini_Read_String('output_root')
    if (outroot /= '') outroot = trim(outroot) // '_'

    highL_unlensed_cl_template = Ini_Read_String_Default('highL_unlensed_cl_template',highL_unlensed_cl_template)

    call CAMB_SetDefParams(P)

    P%WantScalars = Ini_Read_Logical('get_scalar_cls')
    P%WantVectors = Ini_Read_Logical('get_vector_cls',.false.)
    P%WantTensors = Ini_Read_Logical('get_tensor_cls',.false.)

    P%OutputNormalization=outNone
    output_factor = Ini_Read_Double('CMB_outputscale',1.d0)

    P%WantCls= P%WantScalars .or. P%WantTensors .or. P%WantVectors

    P%PK_WantTransfer=Ini_Read_Logical('get_transfer')

    AccuracyBoost  = Ini_Read_Double('accuracy_boost',AccuracyBoost)
    lAccuracyBoost = Ini_Read_Real('l_accuracy_boost',lAccuracyBoost)
    HighAccuracyDefault = Ini_Read_Logical('high_accuracy_default',HighAccuracyDefault)

    P%NonLinear = Ini_Read_Int('do_nonlinear',NonLinear_none)

    P%DoLensing = .false.
    if (P%WantCls) then
        if (P%WantScalars  .or. P%WantVectors) then
            P%Max_l = Ini_Read_Int('l_max_scalar')
            P%Max_eta_k = Ini_Read_Double('k_eta_max_scalar',P%Max_l*2._dl)
            if (P%WantScalars) then
                P%DoLensing = Ini_Read_Logical('do_lensing',.false.)
                if (P%DoLensing) lensing_method = Ini_Read_Int('lensing_method',1)
            end if
            if (P%WantVectors) then
                if (P%WantScalars .or. P%WantTensors) stop 'Must generate vector modes on their own'
                i = Ini_Read_Int('vector_mode')
                if (i==0) then
                    vec_sig0 = 1
                    Magnetic = 0
                else if (i==1) then
                    Magnetic = -1
                    vec_sig0 = 0
                else
                    stop 'vector_mode must be 0 (regular) or 1 (magnetic)'
                end if
            end if
        end if

        if (P%WantTensors) then
            P%Max_l_tensor = Ini_Read_Int('l_max_tensor')
            P%Max_eta_k_tensor =  Ini_Read_Double('k_eta_max_tensor',Max(500._dl,P%Max_l_tensor*2._dl))
        end if
    endif

    !  Read initial parameters.

    call DarkEnergy_ReadParams(DefIni)

    P%h0     = Ini_Read_Double('hubble')

    if (Ini_Read_Logical('use_physical',.false.)) then
        P%omegab = Ini_Read_Double('ombh2')/(P%H0/100)**2
        P%omegac = Ini_Read_Double('omch2')/(P%H0/100)**2
        P%omegan = Ini_Read_Double('omnuh2')/(P%H0/100)**2
        P%omegav = 1- Ini_Read_Double('omk') - P%omegab-P%omegac - P%omegan
    else
        P%omegab = Ini_Read_Double('omega_baryon')
        P%omegac = Ini_Read_Double('omega_cdm')
        P%omegav = Ini_Read_Double('omega_lambda')
        P%omegan = Ini_Read_Double('omega_neutrino')
    end if

    P%tcmb   = Ini_Read_Double('temp_cmb',COBE_CMBTemp)

    ! EFTCAMB MOD START

        ! 1) Initialization of EFTCAMB flags.

    if ( .not. compile_time_eftcamb ) then

        P%EFTflag = Ini_Read_Int('EFTflag',0)

        P%PureEFTmodelOmega  = Ini_Read_Int('PureEFTmodelOmega',0)
        P%PureEFTmodelGamma1 = Ini_Read_Int('PureEFTmodelGamma1',0)
        P%PureEFTmodelGamma2 = Ini_Read_Int('PureEFTmodelGamma2',0)
        P%PureEFTmodelGamma3 = Ini_Read_Int('PureEFTmodelGamma3',0)
        P%PureEFTmodelGamma4 = Ini_Read_Int('PureEFTmodelGamma4',0)
        P%PureEFTmodelGamma5 = Ini_Read_Int('PureEFTmodelGamma5',0)
        P%PureEFTmodelGamma6 = Ini_Read_Int('PureEFTmodelGamma6',0)

        P%DesignerEFTmodel = Ini_Read_Int('DesignerEFTmodel',1)
        P%AltParEFTmodel   = Ini_Read_Int('AltParEFTmodel',1)
        P%FullMappingEFTmodel = Ini_Read_Int('FullMappingEFTmodel',1)

    else if ( compile_time_eftcamb ) then

        P%EFTflag = CT_EFTflag

        P%PureEFTmodelOmega  = CT_PureEFTmodelOmega
        P%PureEFTmodelGamma1 = CT_PureEFTmodelGamma1
        P%PureEFTmodelGamma2 = CT_PureEFTmodelGamma2
        P%PureEFTmodelGamma3 = CT_PureEFTmodelGamma3
        P%PureEFTmodelGamma4 = CT_PureEFTmodelGamma4
        P%PureEFTmodelGamma5 = CT_PureEFTmodelGamma5
        P%PureEFTmodelGamma6 = CT_PureEFTmodelGamma6

        P%DesignerEFTmodel    = CT_DesignerEFTmodel
        P%AltParEFTmodel      = CT_AltParEFTmodel
        P%FullMappingEFTmodel = CT_FullMappingEFTmodel

    end if

    ! 2) Initialization of EFTCAMB model properties flags.

    if ( .not. compile_time_eftcamb ) then

        ! read the DE eos model selection flag:
        P%EFTwDE = Ini_Read_Int('EFTwDE',0)
        ! read pure EFT Horndeski model selection flag:
        P%PureEFTHorndeski = Ini_Read_Logical('PureEFTHorndeski',.false.)
        ! read RPH model selection flags:
        P%RPHmassPmodel      = Ini_Read_Int('RPHmassPmodel',0)
        P%RPHkineticitymodel = Ini_Read_Int('RPHkineticitymodel',0)
        P%RPHbraidingmodel   = Ini_Read_Int('RPHbraidingmodel',0)
        P%RPHtensormodel     = Ini_Read_Int('RPHtensormodel',0)
        ! read the Horava Solar System Free flag:
        P%HoravaSolarSystem  = Ini_Read_Logical('HoravaSolarSystem',.false.)

    else if ( compile_time_eftcamb ) then

        P%EFTwDE             = CT_EFTwDE
        P%PureEFTHorndeski   = CT_PureEFTHorndeski
        P%RPHmassPmodel      = CT_RPHmassPmodel
        P%RPHkineticitymodel = CT_RPHkineticitymodel
        P%RPHbraidingmodel   = CT_RPHbraidingmodel
        P%RPHtensormodel     = CT_RPHtensormodel
        P%HoravaSolarSystem  = CT_HoravaSolarSystem

    end if

    ! 3) Initialization of EFTCAMB model parameters.

    ! read the DE eos parameters:
    P%EFTw0  = Ini_Read_Double('EFTw0',-1._dl)
    P%EFTwa  = Ini_Read_Double('EFTwa',0._dl)
    P%EFTwn  = Ini_Read_Double('EFTwn',2._dl)
    P%EFTwat = Ini_Read_Double('EFTwat',1._dl)
    P%EFtw2  = Ini_Read_Double('EFtw2',0._dl)
    P%EFTw3  = Ini_Read_Double('EFTw3',0._dl)
    ! read pure EFT parameters:
    P%EFTOmega0    = Ini_Read_Double('EFTOmega0', 0.0_dl)
    P%EFTOmegaExp  = Ini_Read_Double('EFTOmegaExp', 0.0_dl)
    P%EFTGamma10   = Ini_Read_Double('EFTGamma10', 0.0_dl)
    P%EFTGamma1Exp = Ini_Read_Double('EFTGamma1Exp', 0.0_dl)
    P%EFTGamma20   = Ini_Read_Double('EFTGamma20', 0.0_dl)
    P%EFTGamma2Exp = Ini_Read_Double('EFTGamma2Exp', 0.0_dl)
    P%EFTGamma30   = Ini_Read_Double('EFTGamma30', 0.0_dl)
    P%EFTGamma3Exp = Ini_Read_Double('EFTGamma3Exp', 0.0_dl)
    P%EFTGamma40   = Ini_Read_Double('EFTGamma40', 0.0_dl)
    P%EFTGamma4Exp = Ini_Read_Double('EFTGamma4Exp', 0.0_dl)
    P%EFTGamma50   = Ini_Read_Double('EFTGamma50', 0.0_dl)
    P%EFTGamma5Exp = Ini_Read_Double('EFTGamma5Exp', 0.0_dl)
    P%EFTGamma60   = Ini_Read_Double('EFTGamma60', 0.0_dl)
    P%EFTGamma6Exp = Ini_Read_Double('EFTGamma6Exp', 0.0_dl)
    ! read f(R) parameters:
    P%EFTB0 = Ini_Read_Double('EFTB0', 0.0_dl)
    ! read RPH parameters:
    P%RPHmassP0        = Ini_Read_Double('RPHmassP0', 0.0_dl)
    P%RPHmassPexp      = Ini_Read_Double('RPHmassPexp', 0.0_dl)
    P%RPHkineticity0   = Ini_Read_Double('RPHkineticity0', 0.0_dl)
    P%RPHkineticityexp = Ini_Read_Double('RPHkineticityexp', 0.0_dl)
    P%RPHbraiding0     = Ini_Read_Double('RPHbraiding0', 0.0_dl)
    P%RPHbraidingexp   = Ini_Read_Double('RPHbraidingexp', 0.0_dl)
    P%RPHtensor0       = Ini_Read_Double('RPHtensor0', 0.0_dl)
    P%RPHtensorexp     = Ini_Read_Double('RPHtensorexp', 0.0_dl)
    ! read Horava parameters:
    P%Horava_xi      = Ini_Read_Double('Horava_xi', 0.0_dl)
    P%Horava_lambda  = Ini_Read_Double('Horava_lambda', 0.0_dl)
    P%Horava_eta     = Ini_Read_Double('Horava_eta', 0.0_dl)

    ! EFTCAMB MOD END

    P%yhe    = Ini_Read_Double('helium_fraction',0.24_dl)
    P%Num_Nu_massless  = Ini_Read_Double('massless_neutrinos')

    P%Nu_mass_eigenstates = Ini_Read_Int('nu_mass_eigenstates',1)
    if (P%Nu_mass_eigenstates > max_nu) stop 'too many mass eigenstates'

    numstr = Ini_Read_String('massive_neutrinos')
    read(numstr, *) nmassive
    if (abs(nmassive-nint(nmassive))>1e-6) stop 'massive_neutrinos should now be integer (or integer array)'
    read(numstr,*, end=100, err=100) P%Nu_Mass_numbers(1:P%Nu_mass_eigenstates)
    P%Num_Nu_massive = sum(P%Nu_Mass_numbers(1:P%Nu_mass_eigenstates))

    if (P%Num_Nu_massive>0) then
        P%share_delta_neff = Ini_Read_Logical('share_delta_neff', .true.)
        numstr = Ini_Read_String('nu_mass_degeneracies')
        if (P%share_delta_neff) then
            if (numstr/='') write (*,*) 'WARNING: nu_mass_degeneracies ignored when share_delta_neff'
        else
            if (numstr=='') stop 'must give degeneracies for each eigenstate if share_delta_neff=F'
            read(numstr,*) P%Nu_mass_degeneracies(1:P%Nu_mass_eigenstates)
        end if
        numstr = Ini_Read_String('nu_mass_fractions')
        if (numstr=='') then
            if (P%Nu_mass_eigenstates >1) stop 'must give nu_mass_fractions for the eigenstates'
            P%Nu_mass_fractions(1)=1
        else
            read(numstr,*) P%Nu_mass_fractions(1:P%Nu_mass_eigenstates)
        end if
    end if

    !JD 08/13 begin changes for nonlinear lensing of CMB + LSS compatibility
    !P%Transfer%redshifts -> P%Transfer%PK_redshifts and P%Transfer%num_redshifts -> P%Transfer%PK_num_redshifts
    !in the P%WantTransfer loop.
    if (((P%NonLinear==NonLinear_lens .or. P%NonLinear==NonLinear_both) .and. P%DoLensing) &
        .or. P%PK_WantTransfer) then
        P%Transfer%high_precision=  Ini_Read_Logical('transfer_high_precision',.false.)
    else
        P%transfer%high_precision = .false.
    endif

    if (P%PK_WantTransfer)  then
        P%WantTransfer  = .true.
        P%transfer%kmax          =  Ini_Read_Double('transfer_kmax')
        P%transfer%k_per_logint  =  Ini_Read_Int('transfer_k_per_logint')
        P%transfer%PK_num_redshifts =  Ini_Read_Int('transfer_num_redshifts')

        transfer_interp_matterpower = Ini_Read_Logical('transfer_interp_matterpower ', transfer_interp_matterpower)
        transfer_power_var = Ini_read_int('transfer_power_var',transfer_power_var)
        if (P%transfer%PK_num_redshifts > max_transfer_redshifts) stop 'Too many redshifts'
        do i=1, P%transfer%PK_num_redshifts
            P%transfer%PK_redshifts(i)  = Ini_Read_Double_Array('transfer_redshift',i,0._dl)
            transferFileNames(i)     = Ini_Read_String_Array('transfer_filename',i)
            MatterPowerFilenames(i)  = Ini_Read_String_Array('transfer_matterpower',i)
            if (TransferFileNames(i) == '') then
                TransferFileNames(i) =  trim(numcat('transfer_',i))//'.dat'
            end if
            if (MatterPowerFilenames(i) == '') then
                MatterPowerFilenames(i) =  trim(numcat('matterpower_',i))//'.dat'
            end if
            if (TransferFileNames(i)/= '') &
                TransferFileNames(i) = trim(outroot)//TransferFileNames(i)
            if (MatterPowerFilenames(i) /= '') &
                MatterPowerFilenames(i)=trim(outroot)//MatterPowerFilenames(i)
        end do
    else
        P%Transfer%PK_num_redshifts = 1
        P%Transfer%PK_redshifts = 0
    end if

    if ((P%NonLinear==NonLinear_lens .or. P%NonLinear==NonLinear_both) .and. P%DoLensing) then
        P%WantTransfer  = .true.
        call Transfer_SetForNonlinearLensing(P%Transfer)
    end if

    call Transfer_SortAndIndexRedshifts(P%Transfer)
    !JD 08/13 end changes

    P%transfer%kmax=P%transfer%kmax*(P%h0/100._dl)

    Ini_fail_on_not_found = .false.

    DebugParam = Ini_Read_Double('DebugParam',DebugParam)
    ALens = Ini_Read_Double('Alens',Alens)

    call Reionization_ReadParams(P%Reion, DefIni)
    call InitialPower_ReadParams(P%InitPower, DefIni, P%WantTensors)
    call Recombination_ReadParams(P%Recomb, DefIni)
    if (Ini_HasKey('recombination')) then
        i = Ini_Read_Int('recombination',1)
        if (i/=1) stop 'recombination option deprecated'
    end if

    call Bispectrum_ReadParams(BispectrumParams, DefIni, outroot)

    if (P%WantScalars .or. P%WantTransfer) then
        P%Scalar_initial_condition = Ini_Read_Int('initial_condition',initial_adiabatic)
        if (P%Scalar_initial_condition == initial_vector) then
            P%InitialConditionVector=0
            numstr = Ini_Read_String('initial_vector',.true.)
            read (numstr,*) P%InitialConditionVector(1:initial_iso_neutrino_vel)
        end if
        if (P%Scalar_initial_condition/= initial_adiabatic) use_spline_template = .false.
    end if

    if (P%WantScalars) then
        ScalarFileName = trim(outroot)//Ini_Read_String('scalar_output_file')
        LensedFileName =  trim(outroot) //Ini_Read_String('lensed_output_file')
        LensPotentialFileName =  Ini_Read_String('lens_potential_output_file')
        if (LensPotentialFileName/='') LensPotentialFileName = concat(outroot,LensPotentialFileName)
        ScalarCovFileName =  Ini_Read_String_Default('scalar_covariance_output_file','scalCovCls.dat',.false.)
        if (ScalarCovFileName/='') then
            has_cl_2D_array = .true.
            ScalarCovFileName = concat(outroot,ScalarCovFileName)
        end if
    end if
    if (P%WantTensors) then
        TensorFileName =  trim(outroot) //Ini_Read_String('tensor_output_file')
        if (P%WantScalars)  then
            TotalFileName =  trim(outroot) //Ini_Read_String('total_output_file')
            LensedTotFileName = Ini_Read_String('lensed_total_output_file')
            if (LensedTotFileName/='') LensedTotFileName= trim(outroot) //trim(LensedTotFileName)
        end if
    end if
    if (P%WantVectors) then
        VectorFileName =  trim(outroot) //Ini_Read_String('vector_output_file')
    end if

#ifdef WRITE_FITS
    if (P%WantCls) then
        FITSfilename =  trim(outroot) //Ini_Read_String('FITS_filename',.true.)
        if (FITSfilename /='') then
            inquire(file=FITSfilename, exist=bad)
            if (bad) then
                open(unit=18,file=FITSfilename,status='old')
                close(18,status='delete')
            end if
        end if
    end if
#endif


    Ini_fail_on_not_found = .false.

    !optional parameters controlling the computation

    P%AccuratePolarization = Ini_Read_Logical('accurate_polarization',.true.)
    P%AccurateReionization = Ini_Read_Logical('accurate_reionization',.false.)
    P%AccurateBB = Ini_Read_Logical('accurate_BB',.false.)
    P%DerivedParameters = Ini_Read_Logical('derived_parameters',.true.)

    version_check = Ini_Read_String('version_check')
    if (version_check == '') then
        !tag the output used parameters .ini file with the version of CAMB being used now
        call TNameValueList_Add(DefIni%ReadValues, 'version_check', version)
    else if (version_check /= version) then
        write(*,*) 'WARNING: version_check does not match this CAMB version'
    end if
    !Mess here to fix typo with backwards compatibility
    if (Ini_HasKey('do_late_rad_trunction')) then
        DoLateRadTruncation = Ini_Read_Logical('do_late_rad_trunction',.true.)
        if (Ini_HasKey('do_late_rad_truncation')) stop 'check do_late_rad_xxxx'
    else
        DoLateRadTruncation = Ini_Read_Logical('do_late_rad_truncation',.true.)
    end if

    if (HighAccuracyDefault) then
        DoTensorNeutrinos = .true.
    else
        DoTensorNeutrinos = Ini_Read_Logical('do_tensor_neutrinos',DoTensorNeutrinos )
    end if
    FeedbackLevel = Ini_Read_Int('feedback_level',FeedbackLevel)

    P%MassiveNuMethod  = Ini_Read_Int('massive_nu_approx',Nu_best)

    ThreadNum      = Ini_Read_Int('number_of_threads',ThreadNum)
    use_spline_template = Ini_Read_Logical('use_spline_template',use_spline_template)

    if (do_bispectrum) then
        lSampleBoost   = 50
    else
        lSampleBoost   = Ini_Read_Double('l_sample_boost',lSampleBoost)
    end if

    call Ini_Close

    if (.not. CAMB_ValidateParams(P)) stop 'Stopped due to parameter error'

    FeedbackLevel = 0

    write(*,"(a,a)") 'Model: ', trim(outroot)

    t1 = omp_get_wtime()

    do i=1, 10

        if (global_error_flag==0) call CAMB_GetResults(P)
        if (global_error_flag/=0) then
            write (*,*) 'Error result '//trim(global_error_message)
            stop
        endif

    end do

    t2 = omp_get_wtime()

    t2 = (t2 -t1)/10._dl

    write(*,"(a,F9.3,a)") 'Timing: ', t2, ' (sec).'

    call CAMB_cleanup
    stop

100 stop 'Must give num_massive number of integer physical neutrinos for each eigenstate'

end program driver
