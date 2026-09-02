! =================================================================================================================================
! MODULE       : control
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        "control" module contains subroutines to initialize run time control parameters. 
!!
!!\n DESCRIPTION: 
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_parameters/control.f90 $ 
!! $Date: 2026-05-15 13:34:49 +0200 (ven. 15 mai 2026) $
!! $Revision: 9546 $ 
!! \n
!_ ================================================================================================================================

MODULE control
  
  USE constantes_soil
  USE constantes_var
  USE pft_parameters
  USE dynamic_parameters
  USE vertical_soil

  IMPLICIT NONE

CONTAINS
  
!! ================================================================================================================================
!! SUBROUTINE   : control_initialize 
!!
!>\BRIEF        This subroutine reads the configuration flags which control the behaviour of the model
!!              This subroutine was previsouly named intsurf_config and located in intersurf module. 
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE control_initialize(kjpindex)

    IMPLICIT NONE

    !! 0. Variables and parameters declaration
    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                 :: kjpindex              !! Number of grid cells in the spatial domain (-)

    !! 0.4 Local variables
    LOGICAL                                    :: conflict              !! Check whether there are conflicts
    LOGICAL                                    :: temp                  !! Temp to check for conflicts  
    INTEGER(i_std)                             :: jv                    !! Local index variable
    INTEGER(i_std)                             :: ier                   !! Error handeling
    LOGICAL                                    :: hydrol_cwrr_test      !! Temporary test variable
    LOGICAL                                    :: ok_co2_test           !! Temporary test variable
    LOGICAL                                    :: ok_explicitsnow_test  !! Temporary test variable
    ! control initialisation with sechiba
    ok_sechiba = .TRUE.

    ! Start reading options from parameter file 
    !
    !Config key   = NC_RESTART_COMPRESSION 
    !Config Desc  = Restart netcdf outputs file are written in compression mode 
    !Config If    =
    !Config Def   = n
    !Config Help  = This flag allows the user to decide if the restart netcdf
    !Config         output files are compressed by default  
    !Config Units = [FLAG]
    !
    nc_restart_compression = .TRUE.
    CALL getin_p('NC_RESTART_COMPRESSION', nc_restart_compression)
    IF (printlev>=1) WRITE(numout,*) "Netcdf restart compression is : ", nc_restart_compression

    !Config Key   = SOILTYPE_CLASSIF
    !Config Desc  = Type of classification used for the map of soil types 
    !Config Def   = zobler
    !Config If    = !IMPOSE_VEG
    !Config Help  = The classification used in the file that we use here 
    !Config         There are three classification supported:  
    !Config         Zobler (7 converted to 3) and USDA (12) 
    !Config Units = [-]
    ! 
    soil_classif = 'zobler'
    CALL getin_p('SOILTYPE_CLASSIF',soil_classif)
    SELECTCASE (soil_classif)
    CASE ('zobler','none')
       nscm = nscm_usda ! owing to the fao2usda pointer
    CASE ('usda')
       nscm = nscm_usda
    CASE DEFAULT
       WRITE(numout,*) "Unsupported soil type classification: soil_classif=",soil_classif
       WRITE(numout,*) "Choose between zobler, usda and none according to the map"
       CALL ipslerr(3,'control_initialize','Bad choice of soil_classif',&
            'Choose between zobler, usda and none','')
    ENDSELECT
    IF (printlev>=1) WRITE(numout,*)'soil_classif,nscm', soil_classif, nscm

    !Config Key   = ENERGY_CONTROL
    !Config Desc  = A flag that controls severeal other flags related to  the energy budget 
    !Config Desc    scheme (enerbil or multi-layer energy budget) and water stress.
    !Config Def   = 1
    !Config If    = OK_SECHIBA
    !Config Help  = Flag that automatically controls several other flags related
    !Config Help    to multi-layering (1/2/3/4/5). 
    !Config Help    1 - DEFAULT uses the enerbil module in combination with the
    !Config Help    hydraulic architecture (ok_hydrol_arch and ok_gs_feedback
    !Config Help    true, while ok_mleb and ok_impose_canopy_structure are set
    !Config Help    to false). Uses a dynamic approach to bare soil thus 
    !Config Help    ok_bare_soil_new = false 
    !Config Help    2 - option to use enerbil module and original water stress
    !Config Help    (not hydraulic architecture). Uses a dynamic approach to bare soil thus 
    !Config Help    ok_bare_soil_new = false 
    !Config Help    3 - The energy budget is calculated using the multi-layer
    !Config Help    energy scheme with a single laye: ok_hydrol_arch,
    !Config Help    ok_gs_feedback, ok_impose_canopy_structure and ok_mleb all
    !Config Help    TRUE, but The energy budget is only calculated for a single
    !Config Help    layer (jnlvls is 1,jnlvls_under is 0,jnlvls_canopy is 1,jnlvls_over is 0).
    !Config Help    No mleb output, ok_mleb_history_file is set to FALSE. Uses a 
    !Config Help    dynamic approach to bare soil thus ok_bare_soil_new = false 
    !Config Help    4 - multi-layer energy budget: ok_hydrol_arch, ok_gs_feedback are false
    !Config Help    whereas ok_mleb is TRUE. ok_impose_canopy_structure is False, and the
    !Config Help    energy budget is calculated for multiple layers 
    !Config Help    (jnlvls is 29,jnlvls_under is 10,jnlvls_canopy is 10,jnlvls_over is 9).    
    !Config Help    No mleb output, ok_mleb_history_file is set to  FALSE. Uses an ecological 
    !Config Help    approach to bare soil thus ok_bare_soil_new = true 
    !Config Help    5 - user specific: user specific settings for these
    !Config Help    controls and layers as defined in the run.def by the user.
    !Config Units = [FLAG]
    ENERGY_CONTROL = 1 
    CALL getin_p('ENERGY_CONTROL', ENERGY_CONTROL)

    ! set flags according to the ENERGY_CONTROL flag.
    IF( ENERGY_CONTROL .EQ. 1) THEN  ! DEFAULT

       ok_hydrol_arch = .TRUE.
       ok_gs_feedback = .TRUE.
       ok_mleb = .FALSE.
       hack_can_structure = .FALSE.
       ok_mleb_history_file = .FALSE.
       ok_bare_soil_new = .FALSE.

       ! Search for possible conflicts
       conflict = .FALSE.
       temp = ok_hydrol_arch
       CALL getin_p('OK_HYDROL_ARCH', temp)
       IF (ok_hydrol_arch .NEQV. temp) conflict = .TRUE.
       temp = ok_gs_feedback
       CALL getin_p('OK_GS_FEEDBACK', temp)
       IF (ok_gs_feedback .NEQV. temp) conflict = .TRUE.
       temp = ok_mleb
       CALL getin_p('OK_MLEB', temp)
       IF (ok_mleb .NEQV. temp) conflict = .TRUE.
       temp = hack_can_structure
       CALL getin_p('OK_IMPOSE_CAN_STRUCTURE', temp)
       IF (hack_can_structure .NEQV. temp) conflict = .TRUE.
       temp = ok_mleb_history_file
       CALL getin_p('OK_MLEB_HISTORY_FILE', temp)
       IF (ok_mleb_history_file .NEQV. temp) conflict = .TRUE.
       IF (conflict) THEN
          CALL ipslerr(3,'control_initialize', &
               'Some of the parameter values set in the', &
               'run.def or set as a default are overwritten by the',&
               'values implied by the ENERGY_CONTROL flag')
       ENDIF

    ELSEIF( ENERGY_CONTROL .EQ. 2 ) THEN  ! enerbil module

       IF (printlev>=1) WRITE(numout,*) 'ENERGY_CONTROL=2 is set. '
       IF (printlev>=1) WRITE(numout,*) 'enerbil module will be used and no of the options in mleb and hydrol_arch will be used'
       ok_hydrol_arch = .FALSE.
       ok_gs_feedback = .FALSE.
       ok_mleb = .FALSE.
       hack_can_structure = .FALSE.
       ok_mleb_history_file = .FALSE.
       ok_bare_soil_new = .FALSE.
            ! Search for possible conflicts
       conflict = .FALSE.
       temp = ok_hydrol_arch
       CALL getin_p('OK_HYDROL_ARCH', temp)
       IF (ok_hydrol_arch .NEQV. temp) conflict = .TRUE.
       temp = ok_gs_feedback
       CALL getin_p('OK_GS_FEEDBACK', temp)
       IF (ok_gs_feedback .NEQV. temp) conflict = .TRUE.
       temp = ok_mleb
       CALL getin_p('OK_MLEB', temp)
       IF (ok_mleb .NEQV. temp) conflict = .TRUE.
       temp = hack_can_structure
       CALL getin_p('OK_IMPOSE_CAN_STRUCTURE', temp)
       IF (hack_can_structure .NEQV. temp) conflict = .TRUE.
       temp = ok_mleb_history_file
       CALL getin_p('OK_MLEB_HISTORY_FILE', temp)
       IF (ok_mleb_history_file .NEQV. temp) conflict = .TRUE.
       temp = ok_bare_soil_new
       CALL getin_p('OK_BARE_SOIL_NEW', temp)
       IF (ok_bare_soil_new .NEQV. temp) conflict = .TRUE.
       IF (conflict) THEN
          CALL ipslerr(3,'control_initialize', &
               'Some of the parameter values set in the', &
               'run.def or set as a default are overwritten by the',&
               'values implied by the ENERGY_CONTROL flag')
       ENDIF
    ELSEIF( ENERGY_CONTROL .EQ. 3 ) THEN  ! single layer multi-layer energy scheme

       ok_hydrol_arch = .TRUE.
       ok_gs_feedback = .TRUE.
       ok_mleb = .TRUE.
       hack_can_structure = .TRUE.
       ok_mleb_history_file = .FALSE.
       ok_bare_soil_new = .FALSE.

       ! Search for possible conflicts
       conflict = .FALSE.
       temp = ok_hydrol_arch
       CALL getin_p('OK_HYDROL_ARCH', temp)
       IF (ok_hydrol_arch .NEQV. temp) conflict = .TRUE.
       temp = ok_gs_feedback
       CALL getin_p('OK_GS_FEEDBACK', temp)
       IF (ok_gs_feedback .NEQV. temp) conflict = .TRUE.
       temp = ok_mleb
       CALL getin_p('OK_MLEB', temp)
       IF (ok_mleb .NEQV. temp) conflict = .TRUE.
       temp = hack_can_structure
       CALL getin_p('OK_IMPOSE_CAN_STRUCTURE', temp)
       IF (hack_can_structure .NEQV. temp) conflict = .TRUE.
       temp = ok_mleb_history_file
       CALL getin_p('OK_MLEB_HISTORY_FILE', temp)
       IF (ok_mleb_history_file .NEQV. temp) conflict = .TRUE.
       IF (conflict) THEN
          CALL ipslerr(3,'control_initialize', &
               'Some of the parameter values set in the', &
               'run.def or set as a default are overwritten by the',&
               'values implied by the ENERGY_CONTROL flag')
       ENDIF

    ELSEIF( ENERGY_CONTROL .EQ. 4 ) THEN  ! multi-layer

       ok_hydrol_arch = .FALSE.
       ok_gs_feedback = .FALSE.
       ok_mleb = .TRUE.
       hack_can_structure = .FALSE.
       ok_mleb_history_file = .FALSE.
       ok_bare_soil_new = .FALSE.

       ! Search for possible conflicts
       conflict = .FALSE.
       temp = ok_hydrol_arch
       CALL getin_p('OK_HYDROL_ARCH', temp)
       IF (ok_hydrol_arch .NEQV. temp) conflict = .FALSE.
       temp = ok_gs_feedback
       CALL getin_p('OK_GS_FEEDBACK', temp)
       IF (ok_gs_feedback .NEQV. temp) conflict = .FALSE.
       temp = ok_mleb
       CALL getin_p('OK_MLEB', temp)
       IF (ok_mleb .NEQV. temp) conflict = .TRUE.
       temp = hack_can_structure
       CALL getin_p('OK_IMPOSE_CAN_STRUCTURE', temp)
       IF (hack_can_structure .NEQV. temp) conflict = .TRUE.
       temp = ok_mleb_history_file
       CALL getin_p('OK_MLEB_HISTORY_FILE', temp)
       IF (ok_mleb_history_file .NEQV. temp) conflict = .TRUE.
       IF (conflict) THEN
          CALL ipslerr(3,'control_initialize', &
               'Some of the parameter values set in the', &
               'run.def or set as a default are overwritten by the',&
               'values implied by the ENERGY_CONTROL flag')
       ENDIF

    ELSEIF( ENERGY_CONTROL .EQ. 5 ) THEN 

       !Config Key   = OK_HYDROL_ARCH
       !Config Desc  = Activates the hydraulic architecture by Tuzet et al. (2017)
       !Config If    = OK_SECHIBA
       !Config Def   = y
       !Config Help  = Flag that activates the hydraulic architecture routine (true/false)
       !Config Help    The trunk version of ORCHIDEE (false) uses soil water as a 
       !Config Help    proxy for water stress and applies the stress to Vcmax.
       !Config Help    When set to true the hydraulic architecture of the vegetation
       !Config Help    is accounted for to calculate the amount of water that 
       !Config Help    can be transported through the plant given the soil and leaf
       !Config Help    potential and the conductivities of the roots, wood and 
       !Config Help    leaves. Water supply through the plant is compared against
       !Config Help    the atmospheric demand for water. If the supply is smaller
       !Config Help    then the demand, the plant experiences water stress and the
       !Config Help    stomata will be closed (water stress is now on gs rather 
       !Config Help    than Vcmax). Note that whether stomatal regulation is used or 
       !Config Help    not is controled by a separate flag: ok_gs_feedback.
       !Config Units = [FLAG]
       ok_hydrol_arch = .FALSE.
       CALL getin_p('OK_HYDROL_ARCH', ok_hydrol_arch)
       !
       !Config Key   = OK_GS_FEEDBACK
       !Config Desc  = Debug option for OK_HYDROL_ARCH  
       !Config If    = OK_SECHIBA, OK_HYDROL_ARCH
       !Config Def   = y
       !Config Help  = Flag that activates water stress on stomata (true/false)
       !Config Help    This flag is for debugging only! It allows developers
       !Config Help    to calculate GPP without any water stress. If the model is 
       !Config Help    used in production mode and ok_hydrol_arch is true this 
       !Config Help    flag should be true as well.
       !Config Units = [FLAG]
       ok_gs_feedback = .FALSE.
       CALL getin_p('OK_GS_FEEDBACK', ok_gs_feedback)
       !
       !Config Key   = OK_MLEB
       !Config Desc  = Activate multi-layer energy budget  
       !Config If    = OK_SECHIBA
       !Config Def   = y
       !Config Help  = Flag that activates the multilayer energy budget (true/false)
       !Config Help    The model uses 10 (default) canopy layers to calculate
       !Config Help    the albedo, transmittance, absorbance and GPP. These canopy
       !Config Help    layers can be combined with 10 (default) layers below and 
       !Config Help    10 layers above the canopy to calculate the energy budget
       !Config Help    (ok_mleb is y). If set to no, this flag will make the model
       !Config Help    use 10 layers for the canopy albedo, transmittance, 
       !Config Help    absorbance and GPP and just a single layer for the energy 
       !Config Help    budget. Be aware that if you wish to run with hydraulic 
       !Config Help    architechture ok_mleb needs to be se to true as well. Furthermore
       !Config Help    if you  wish to run with the original energy scheme (enerbil),
       !Config Help    set the layers for mleb to 1.
       !Config Units = [FLAG]
       ok_mleb = .TRUE.
       CALL getin_p('OK_MLEB', ok_mleb)
       !
       !+++CHECK+++
       ! This may no longer be needed. One could get the same functionality
       ! by reading an biomass map in shechiba (lai_map) and not using stomate.
       !Config Key   = OK_IMPOSE_CAN_STRUCTURE
       !Config Desc  = Debug option for OK_MLEB 
       !Config If    = OK_SECHIBA, OK_MLEB
       !Config Def   = n
       !Config Help  = This flag is for debugging only! It allows developers
       !Config Help    to use a prescribed canopy structure rather then the 
       !Config Help    structure calculate by ORCHIDEE. The flag activates the 
       !Config Help    sections of code which directly link the energy budget 
       !Config Help    scheme to the the size and LAI profile of the canopy for the 
       !Config Help    respective PFT and age class that is calculated in stomate, 
       !Config Help    for the albedo. If set to TRUE and the multi-layer budget 
       !Config Help    is activated the model takes LAI profile information and 
       !Config Help    canopy level heights from the run.def. If set to FALSE, and 
       !Config Help    and the multi-layer energy budget is used the profile 
       !Config Help    infomation and canopy levels heights comes from the
       !Config Help    PGap-based processes for calculation of stand profile 
       !Config Help    information in stomate.
       !Config Units = [FLAG]
       hack_can_structure = .FALSE.
       CALL getin_p('OK_IMPOSE_CAN_STRUCTURE', hack_can_structure)
       !++++++++++
       !
       !Config Key   = MLEB_NETCDF_FLAG
       !Config Desc  = Debug option for OK_MLEB 
       !Config If    = OK_SECHIBA, OK_MLEB
       !Config Def   = n
       !Config Help  = Flag that controls the writing of an output file with the 
       !Config Help    multi-layer energy simulations (true/false). Note that this
       !Config Help    a large file and writing it slows down the code.
       !Config Units = [FLAG]
       ok_mleb_history_file = .TRUE.
       CALL getin_p('OK_MLEB_HISTORY_FILE', ok_mleb_history_file)

       !Config Key   = OK_BARE_SOIL_NEW
       !Config Desc  = Flag that controls the view on and calculation of bare soil
       !Config If    = OK_SECHIBA or OK_STOMATE 
       !Config Def   = FALSE
       !Config Help  = Choose between the two options to calculate the bare soil.
       !Config         False uses the classic view: gaps within a canopy should be treated
       !Config         as bare soil. True uses the ecological view: gaps within a canopy are
       !Config         part of the ecosystem and should be treated as such.
       !Config         NOTE: this development needs to be finalized and checked before
       !Config         being used.
       !Config Units = [FLAG]
       ok_bare_soil_new = .FALSE.
       CALL getin_p('OK_BARE_SOIL_NEW', ok_bare_soil_new) 

    ELSE
       WRITE(numout,*) 'Current setting of ENERGY_CONTROL is not implemanted, ENERGY_CONTROL=',ENERGY_CONTROL
       CALL ipslerr_p(3,'control_initialize',&
            'ENERGY_CONTROL can only have integer values from 1-4', '','')
    ENDIF ! (ENERGY_CONTROL)

    ! Share the final settings with the user
    IF (printlev>=1) WRITE(numout,*) 'ENERGY_CONTROL is activated', ENERGY_CONTROL
    IF (printlev>=1) WRITE(numout,*) 'OK_HYDROL_ARCH is activated, ',ok_hydrol_arch
    IF (printlev>=1) WRITE(numout,*) 'OK_GS_FEEDBACK is activated, ',ok_gs_feedback    
    IF (printlev>=1) WRITE(NUmout,*) 'OK_MLEB is activated', ok_mleb
    IF (printlev>=1) WRITE(numout,*) 'OK_IMPOSE_CAN_STRUCTURE is activated', hack_can_structure
    IF (printlev>=1) WRITE(numout,*) 'OK_MLEB_HISTORY_FILE is activated',ok_mleb_history_file
    IF (printlev>=1) WRITE(numout,*) 'OK_BARE_SOIL_NEW is activated',ok_bare_soil_new
    
    ! Final consistency check of the settings
    IF (ok_hydrol_arch .AND. ok_gs_feedback)THEN

       ! This is the default set-up for running hydraulic architecture 
       ! in production mode. Using full functionality
       IF (printlev>=1) WRITE(numout,*) 'Hydraulic architecture, stomatal feedback'
       IF (printlev>=1) WRITE(numout,*) 'and the energy budget are (re)calculated'

    ELSEIF (.NOT. ok_hydrol_arch .AND. ok_gs_feedback) THEN

       ! The debug options should not be used. Not clear whether
       ! the user want to use hydraulic architecture or not
       WRITE(numout,*) 'Hydraulic architecture is not used'
       WRITE(numout,*) 'but the debug option was set to TRUE'
       WRITE(numout,*) 'Not clear whether you want to use '
       WRITE(numout,*) 'hydraulic architecture or not'
       CALL ipslerr(3,'control.f90 - run.def',&
            'setting for stomatal feedback was in conflict with',&
            'the setting of hydraulic architecture','')

    ELSEIF (ok_hydrol_arch .AND. .NOT.ok_gs_feedback) THEN

       ! Warning
       CALL ipslerr(2,'control.f90 - hydraulic architecture',&
            'is run in DEBUG mode',&
            'stomatal feedback is not used','')

    ENDIF

    IF (ok_mleb .AND. hack_can_structure) THEN

       ! Warning
       CALL ipslerr(2,'control.f90 - multi-layer energy budget',&
            'is run in debug mode',&
            'The canopy structure is fixed and imposed','')

    ELSEIF (ok_mleb .AND. .NOT.ok_mleb_history_file) THEN
    
       ! Warning
       CALL ipslerr(2,'control.f90 - multi-layer energy budget',&
            'There won be any output files for the multi-layer energy budget',&
            'If this is not intended, change the run.def','')

    ELSEIF ((.NOT.ok_mleb .AND. hack_can_structure) .OR. &
         (.NOT.ok_mleb .AND. ok_mleb_history_file)) THEN
       
       ! The debug options should not be used. Not clear whether
       ! the user want to use hydraulic architecture or not
       WRITE(numout,*) 'The multi-layer energy budget is not used'
       WRITE(numout,*) 'Debug options are set to FALSE'
       CALL ipslerr(3,'control.f90 - run.def',&
            'inconsistent settings of de debug options',&
            'because the multi-level enery budget itself is not used','')

    ENDIF

    !Config Key   = OK_HYDROL_ARCH_STORAGE
    !Config Desc  = Decides if we use the new hydraulic architecture with water storage compartments
    !Config If    = OK_HYDROL_ARCH
    !Config Def   = n
    !Config Help  = This flag allows the user to decide if the new hydraulic 
    !Config         architecture will be used with water storage inside the
    !Config         vegetation as in Tuzet et al. (2017description.
    !Config Units = [FLAG]
    ok_hydrol_arch_storage = .FALSE.
    CALL getin_p('OK_HYDROL_ARCH_STORAGE', ok_hydrol_arch_storage)
    IF (printlev>=1) WRITE(numout,*) "Tuzet et al. (2017)'s hydraulic architecture is activated with storages: ",ok_hydrol_arch_storage

    !Config Key   = OK_HYDROL_ARCH_DYN_RES
    !Config Desc  = Decides if we use the new hydraulic architecture with dynamic water resistances
    !Config If    = OK_HYDROL_ARCH 
    !Config Def   = n
    !Config Help  = Flag that activates the dynamic resistance equations used by
    !Config         Yao et al. (2021) (Under review at the time this is developed)
    !Config         instead of having parameterised hydraulic resistances as in Tuzet 
    !Config         et al. (2017), the resistances now depend on the water potential of the stage.
    !Config Units = [FLAG]
    ok_hydrol_arch_dyn_res = .FALSE.
    CALL getin_p('OK_HYDROL_ARCH_DYN_RES', ok_hydrol_arch_dyn_res)
    IF (printlev>=1) WRITE(numout,*) "Tuzet et al. (2017)'s hydraulic architecture is activated with dynamic resistances: ",ok_hydrol_arch_dyn_res

    !Config Key   = SPLIT_SOIL_PROPERTIES
    !Config Desc  = Decides if we split the soil properties between the top and inferior layers
    !Config If    = OK_HYDROL_ARCH
    !Config Def   = n
    !Config Help  = This flag allows the user to decide if the new hydraulic 
    !Config         architecture distinguishes the soil properties in the 2
    !Config         soil horizons.
    !Config         Flag that permits to distinguish properties between the superficial
    !Config         soil horizon and the inferior one in the new hydraulic architecture.
    !Config Units = [FLAG]
    split_soil_properties = .FALSE.
    CALL getin_p('SPLIT_SOIL_PROPERTIES', split_soil_properties)
    IF (printlev>=1) WRITE(numout,*) "New hydraulic architecture is activated split soil properties: ",split_soil_properties

    !Config Key   = OK_HYDROL_ARCH_MUFF
    !Config Desc  = Decides if we use the new hydraulic architecture with the muff root absorption description
    !Config If    = OK_HYDROL_ARCH
    !Config Def   = n
    !Config Help  = Flag that activates the radial resolution of Richard's equation
    !Config         around the roots. This resolution is only made when the
    !Config         hydraulic architecture is running. If set to .FALSE., the root 
    !Config         absorption is modelled by a classical resistance term.
    !Config Units = [FLAG]
    ok_hydrol_arch_muff = .FALSE.
    CALL getin_p('OK_HYDROL_ARCH_MUFF', ok_hydrol_arch_muff)
    IF (printlev>=1) WRITE(numout,*) "Tuzet et al. (2017)'s hydraulic architecture is activated with muffs: ",ok_hydrol_arch_muff

    !Config Key   = OK_TUZET_INDIV
    !Config Desc  = Puts Tuzet at tree scale
    !Config If    = OK_HYDROL_ARCH
    !Config Def   = n
    !Config Help  =
    !Config Units = [FLAG]
    ok_tuzet_indiv = .FALSE.
    CALL getin_p('OK_TUZET_INDIV', ok_tuzet_indiv)
    IF (printlev>=1) WRITE(numout,*) "Tuzet at individual tree scale : ", ok_tuzet_indiv

    !Config Key   = OK_NSL_VPD
    !Config Desc  = NSL in Tuzet
    !Config If    = OK_HYDROL_ARCH
    !Config Def   = n
    !Config Help  =
    !Config
    !Config
    !Config
    !Config Units = [FLAG]
    ok_nsl_vpd = .FALSE.
    CALL getin_p('OK_NSL_VPD', ok_nsl_vpd)
    IF (printlev>=1) WRITE(numout,*) "NSL in Tuzet : ", ok_nsl_vpd

    !Config Key   = TUZET_INF_LAYER
    !Config Desc  = Base of water potentials is inferior layer
    !Config If    = OK_HYDROL_ARCH
    !Config Def   = n
    !Config Help  =
    !Config
    !Config
    !Config
    !Config Units = [FLAG]
    tuzet_inf_layer=.FALSE.
    CALL getin_p('TUZET_INF_LAYER', tuzet_inf_layer)
    IF (printlev>=1) WRITE(numout,*) "Base of water potentials is the inferior layer : ", tuzet_inf_layer

    !Config Key   = TEST_FUNCTIONAL_ROOTS_TUZET
    !Config Desc  = Replace previously used root profile by functional
    !Config If    = OK_HYDROL_ARCH
    !Config Def   = y
    !Config Help  =
    !Config
    !Config
    !Config
    !Config Units = [FLAG]
    test_functional_roots_tuzet = .FALSE.
    CALL getin_p('TEST_FUNCTIONAL_ROOTS_TUZET', test_functional_roots_tuzet)
    IF (printlev>=1) WRITE(numout,*) "Gmin has a stress function applied : ", test_functional_roots_tuzet

    !Config Key   = NO_G0
    !Config Desc  = Puts stress function on g0
    !Config If    = OK_HYDROL_ARCH
    !Config Def   = y
    !Config Help  =
    !Config
    !Config
    !Config
    !Config Units = [FLAG]
    no_g0 = .TRUE.
    CALL getin_p('NO_G0', no_g0)
    IF (printlev>=1) WRITE(numout,*) "Gmin has a stress function applied : ", no_g0

    !Config Key   = OK_SAP_FEEDBACK
    !Config Desc  = Activate feedback of vessel mortality on sapwood conductance
    !Config If    = OK_STOMATE
    !Config Def   = FALSE
    !Config Help  =
    !Config Units = [FLAG]
    ok_sap_feedback=.FALSE.
    CALL getin_p('OK_SAP_FEEDBACK', ok_sap_feedback)
    IF (.NOT. ok_hydrol_arch_dyn_res.AND.ok_sap_feedback) THEN
       ok_sap_feedback = .FALSE.
       WRITE(numout,*) 'Sapwood mortality feedback on conductance can only happen with dynamic resistances in the hydraulic architecture'
    ENDIF

    !Config Key   = OK_LAKE_ENERGY
    !Config Desc  = Activate FLake (the energy model of lake)
    !Config If    = OK_SECHIBA
    !Config Def   = n
    !Config Help  = set to TRUE if Lake is modeled
    !Config Units = [FLAG]
    ok_lake_energy = .FALSE.
    CALL getin_p('OK_LAKE_ENERGY', ok_lake_energy)
    IF (printlev>=1) WRITE(numout,*) "LAKE is activated : ",ok_lake_energy

    !Config Key   = OK_SNOW_ICE_LAKE_PAR
    !Config Desc  = Activate FLake snow and ice albedo parametrization
    !Config If    = OK_LAKE_ENERGY
    !Config Def   = y
    !Config Help  = 
    !Config Units = [FLAG]
    ok_snow_ice_lake_par = .TRUE.
    CALL getin_p('OK_SNOW_ICE_LAKE_PAR', ok_snow_ice_lake_par)
    IF (printlev>=1) WRITE(numout,*) "Snow and ice albedo parametrization of FLake is activated : ",ok_snow_ice_lake_par

    !Config Key   = OK_BOTSED_LAKE
    !Config Desc  = Activate FLake parametrization of sediment layer
    !Config If    = OK_LAKE_ENERGY
    !Config Def   = n
    !Config Help  = 
    !Config Units = [FLAG]
    ok_botsed_lake = .FALSE.
    CALL getin_p('OK_BOTSED_LAKE', ok_botsed_lake)
    IF (printlev>=1) WRITE(numout,*) "Sediment layer parametrization of FLake is activated : ",ok_botsed_lake

    !Config Key   = RIVER_ROUTING
    !Config Desc  = Decides if we route the water or not
    !Config If    = OK_SECHIBA
    !Config Def   = y
    !Config Help  = This flag allows the user to decide if the runoff
    !Config         and drainage should be routed to the ocean
    !Config         and to downstream grid boxes.
    !Config Units = [FLAG]
    river_routing = .TRUE.
    CALL getin_p('RIVER_ROUTING', river_routing)
    IF (printlev>=1) WRITE(numout,*) "RIVER routing is activated : ",river_routing

    IF (river_routing) THEN
       !Config Key   = ROUTING_METHOD
       !Config Desc  = Choice of routing module to be used
       !Config If    = RIVER_ROUTING=T
       !Config Def   = subgrid_halfdeg
       !Config Help  = Character string used to switch between routing modules
       !Config         Possible options are
       !Config         subgrid_halfdeg, subgrid_htu,
       !Config         interp_topo_v0, interp_topo (or their old equivalents)
       !Config Units = character string
       
       routing_method='subgrid_halfdeg'
       CALL getin_p("ROUTING_METHOD",routing_method)
       
       ! Handle possibility of old flag inputs
       ! and change them to their new name
       IF      ( trim(routing_method).EQ.'standard' .OR. trim(routing_method).EQ.'subgrid_halfdeg') THEN
         routing_method='subgrid_halfdeg'
       ELSE IF ( trim(routing_method).EQ.'highres'  .OR. trim(routing_method).EQ.'subgrid_htu') THEN
         routing_method='subgrid_htu'
       ELSE IF ( trim(routing_method).EQ.'simple'   .OR. trim(routing_method).EQ.'interp_topo_v0') THEN
         routing_method='interp_topo_v0'
       ELSE IF ( trim(routing_method).EQ.'native'   .OR. trim(routing_method).EQ.'interp_topo') THEN
         routing_method='interp_topo'
       ELSE
         WRITE(numout,*) 'Following routing method is not implemented : ROUTING_METHOD=',routing_method
         CALL ipslerr_p(3,'control_initialize',&
                    'Error in ROUTING_METHOD set up', &
                    'Choose within subgrid_halfdeg, subgrid_htu, interp_topo_v0, interp_topo',&
                    'in run.def or in sechiba.card (flag ROUTING)')
       ENDIF
       IF (printlev>=1) WRITE(numout,*) "ROUTING_METHOD is set to : ",routing_method
    ELSE
       routing_method='none'
    END IF

    !Config Key   = WHICH_RT
    !Config Desc  = Choice of radiation transfer scheme
    !Config If    = 
    !Config Def   = Matrixial
    !Config Help  = Possible options are Iterative, and Matrixial
    !Config Units = character string
    which_rt = 'Matrixial'
    CALL getin_p("WHICH_RT",which_rt)
    IF (printlev>=1) WRITE(numout,*) "Radiation transfer is set to : ",which_rt
    IF (which_rt .NE. 'Iterative' .AND. which_rt .NE. 'Matrixial') THEN
            CALL ipslerr_p(3,'control_initialize',&
                    'Error in RT solver set up, choose WHICH_RT = Iterative or Matrixial', &
                    '','')
    ENDIF

    ! Control for the option HYDROL_CWRR which is not longer existing in the model. 
    ! Check here if in run.def HYDROL_CWRR=n. If that's the case then stop the model and ask the user to remove the flag from run.def
    hydrol_cwrr_test = .TRUE.
    CALL getin_p('HYDROL_CWRR', hydrol_cwrr_test)
    IF (.NOT. hydrol_cwrr_test) THEN
       CALL ipslerr_p(3,'control_initialize',&
            'HYDROL_CWRR=n is set in run.def but this option does not exist any more in ORCHIDEE', &
            'Choisnel hydrolology has been removed and CWRR is now the only hydrology module in ORCHIDEE',&
            'Remove parameter HYDROL_CWRR from run.def')
    END IF

    !Config Key   = DO_IRRIGATION
    !Config Desc  = Should we compute an irrigation flux 
    !Config If    = RIVER_ROUTING 
    !Config Def   = n
    !Config Help  = This parameters allows the user to ask the model
    !Config         to compute an irigation flux. This performed for the
    !Config         on very simple hypothesis. The idea is to have a good
    !Config         map of irrigated areas and a simple function which estimates
    !Config         the need to irrigate.
    !Config Units = [FLAG]
    do_irrigation = .FALSE.
    IF ( river_routing ) CALL getin_p('DO_IRRIGATION', do_irrigation)
    
    !Config Key   = DO_FLOODPLAINS
    !Config Desc  = Should we include floodplains 
    !Config If    = RIVER_ROUTING 
    !Config Def   = n
    !Config Help  = This parameters allows the user to ask the model
    !Config         to take into account the flood plains and return 
    !Config         the water into the soil moisture. It then can go 
    !Config         back to the atmopshere. This tried to simulate 
    !Config         internal deltas of rivers.
    !Config Units = [FLAG]  
    do_floodplains = .FALSE.
    IF ( river_routing ) CALL getin_p('DO_FLOODPLAINS', do_floodplains)

    !Config Key   = OK_SOIL_CARBON_DISCRETIZATION
    !Config Desc  = Activate soil carbon vertical discretization
    !Config If    = OK_STOMATE
    !Config Def   = FALSE
    !Config Help  = Activate soil carbon scheme with vertical discretization and vertical transport of carbon
    !Config Units = [FLAG]
    ok_soil_carbon_discretization=.FALSE.
    CALL getin_p('OK_SOIL_CARBON_DISCRETIZATION', ok_soil_carbon_discretization)
    
    !Config Key   = OK_ICE_SHEET
    !Config Desc  = Activate ice scheme on ice-sheet area
    !Config If    = OK_SECHIBA
    !Config Def   = FALSE
    !Config Help  = Activate ice scheme on ice-sheet area
    !Config Units = [FLAG]
    ok_ice_sheet = .FALSE.
    CALL getin_p('OK_ICE_SHEET', ok_ice_sheet)
    
    !Config Key  = USE_THERMO_NEWDISC
    !Config Desc = Use new subroutines for thermal capacities and conductivities (NewDisC)
    !Config If   = OK_SOIL_CARBON_DISCRETIZATION or USE_SOILC_INSULATION 
    !Config Def  = FALSE
    !Config Help = 
    !Config Units = [FLAG]
    use_thermo_newdisc = .FALSE.
    IF (ok_soil_carbon_discretization .OR. use_soilc_insulation) CALL getin_p('USE_THERMO_NEWDISC', use_thermo_newdisc)

    ! Control of option OK_EXPLICITSNOW which is not longer existing in the model. 
    ! Check here if in run.def OK_EXPLICITSNOW=n. If that's the case then stop the model and ask the user to remove the flag from run.def.
    ok_explicitsnow_test = .TRUE.
    CALL getin_p('OK_EXPLICITSNOW', ok_explicitsnow_test)
    IF (.NOT. ok_explicitsnow_test) THEN
       CALL ipslerr_p(3,'control_initialize',&
            'OK_EXPLICITSNOW=n is set in run.def but this option does not exist any more in ORCHIDEE', &
            'Explicit snow scheme is now always used in ORCHIDEE.',&
            'Remove parameter OK_EXPLICITSNOW from run.def')
    END IF

!!$    ! One of the most frequent problems is a temperature out of range
!!$    ! we provide here a way to catch that in the calling procedure. 
!!$    ! (from Jan Polcher)(true/false)
!!$    !Config Key   = DIAG_QSAT
!!$    !Config Desc  = Catch temperature out of range
!!$    !Config If    = OK_SECHIBA
!!$    !Config Def   = TRUE
!!$    !Config Help  = One of the most frequent problems is a temperature out of range
!!$    !Config Help  = we provide here a way to catch that in the calling procedure.
!!$    !Config Help  = (from Jan Polcher)(true/false)
!!$    !Config Units = [FLAG]
!!$    diag_qsat = .TRUE.
!!$    CALL getin_p('DIAG_QSAT',diag_qsat)      
    !
    !Config Key   = STOMATE_OK_STOMATE
    !Config Desc  = Activate STOMATE?
    !Config If    = OK_SECHIBA
    !Config Def   = y
    !Config Help  = set to TRUE if STOMATE is to be activated
    !Config Units = [FLAG]
    ok_stomate = .TRUE.
    CALL getin_p('STOMATE_OK_STOMATE',ok_stomate)
    IF (printlev>=1) WRITE(numout,*) 'STOMATE is activated: ',ok_stomate

    !Config Key   = OK_VESSEL_MORTALITY
    !Config Desc  = Activate death and recovery of vegetation following hydraulic failure.
    !Config If    = OK_STOMATE
    !Config Def   = FALSE
    !Config Help  = Activate death and recovery of vegetation following hydraulic failure.
    !Config Units = [FLAG]
    ok_vessel_mortality=.FALSE.
    CALL getin_p('OK_VESSEL_MORTALITY', ok_vessel_mortality)
    IF (.NOT. ok_hydrol_arch_dyn_res.AND.ok_vessel_mortality) THEN
       ok_vessel_mortality = .FALSE.
       WRITE(numout,*) 'Vessel mortality can only happen with dynamic resistances in the hydraulic architecture'
    ENDIF

    !Config Key   = READ_LAI
    !Config Desc  = Options for LAI calculation without stomate. 0: impose LAI from biomass, 1: impose LAI from phenology and LAI for site-based, 2:read biomass map 
    !Config If    = OK_STOMATE=FALSE
    !Config Def   = -1
    !Config Help  = If equal to 0: Uses a singel prescribed biomass value for grid-cell, pft and month.
    !Config         If equal to 1: Calculates a site-level LAI time-series based on day of year for phenology and 
    !Config         senescence and minimum and maximum of LAI.
    !Config         If equal to 2: read a biomass map. This option is not yet developped. 
    !Config Units = [-]
    read_lai = -1
    CALL getin_p('READ_LAI',read_lai)
    IF (ok_stomate .AND. read_lai.GE.0) THEN
       WRITE(numout,*) 'ok_stomate = ', ok_stomate
       CALL ipslerr_p(3,'slowproc_init','The read_lai option should be accompanied by ok_stomate=FALSE',&
            'Remove READ_LAI from run.def or deactivate STOMATE','')
    ELSE IF (read_lai > 2) THEN
       WRITE(numout,*) 'Following value for READ_LAI can not be used, read_lai=',read_lai
       CALL ipslerr_p(3,'slowproc_init','The read_lai option can only be -1, 0, 1, or 2. Please change in run.def','','')
    END IF
   

    ! Control for the option STOMATE_OK_CO2 which is not longer existing in the model. 
    ! Check here if in run.def STOMATE_OK_CO2=n. If that's the case then stop the model and ask the user to remove the flag from run.def
    ok_co2_test = .TRUE.
    CALL getin_p('STOMATE_OK_CO2', ok_co2_test)
    IF (.NOT. ok_co2_test) THEN
       CALL ipslerr_p(3,'control_initialize',&
            'STOMATE_OK_CO2=n is set in run.def but this option does not exist any more in ORCHIDEE', &
            'Calculation of beta coefficient using Jarvis formulation has been removed and Farquar formulation is now always used',&
            'Remove parameter STOMATE_OK_CO2 from run.def')
    END IF

    !  
    !Config Key   = DO_WOOD_HARVEST
    !Config Desc  = Activate Wood Harvest ?
    !Config If    = OK_STOMATE
    !Config Def   = n
    !Config Help  = set to TRUE if a prescribed amount of wood is to be harvested.
    !Config         Note that this approach is different from the forest management
    !Config         in which the management strategy is read from a map and the model
    !Config         left to calculate the thinnings and the harvest.
    !Config Units = [FLAG]
    do_wood_harvest = .FALSE.
    CALL getin_p('DO_WOOD_HARVEST',do_wood_harvest)

    !
    !Config Key   = STOMATE_OK_NCYCLE 
    !Config Desc  = Activate dynamic N cycle 
    !Config If    = OK_STOMATE 
    !Config Def   = y 
    !Config Help  = set to TRUE if N cycle is to be activated 
    !Config Units = [FLAG] 
    ok_ncycle = .TRUE. 
    CALL getin_p('STOMATE_OK_NCYCLE',ok_ncycle) 
    IF (printlev>=1) WRITE(numout,*) 'N cycle is activated: ',ok_ncycle  

    !
    !Config Key   = STOMATE_IMPOSE_CN
    !Config Desc  = Impose the CN ratio of leaves 
    !Config If    = OK_STOMATE 
    !Config Def   = n 
    !Config Help  = set to TRUE if IMPOSE_CN is to be activated 
    !Config Units = [FLAG] 
    impose_cn = .FALSE. 
    CALL getin_p('STOMATE_IMPOSE_CN',impose_cn) 
    IF (printlev>=1) WRITE(numout,*) 'CN ratio is imposed: ',impose_cn  
    !
    !Config Key   = RESET_IMPOSE_CN
    !Config Desc  = Reset the CN ratio of leaves 
    !Config If    = OK_STOMATE 
    !Config Def   = n 
    !Config Help  = set to TRUE if RESET_IMPOSE_CN is to be activated 
    !Config Units = [FLAG]  
    reset_impose_cn = .FALSE. 
    CALL getin_p('RESET_IMPOSE_CN',reset_impose_cn) 
    IF (printlev>=1) WRITE(numout,*) 'CN ratio is reset: ',reset_impose_cn  
    !
    !Config Key   = STOMATE_READ_CN
    !Config Desc  = Read the CN ratio of leaves 
    !Config If    = OK_STOMATE 
    !Config Def   = n 
    !Config Help  = set to TRUE if IMPOSE_CN is to be activated read
    !Config Units = [FLAG] 
    ! 
    read_cn = .FALSE. 
    CALL getin_p('STOMATE_READ_CN',read_cn) 
    IF (printlev>=1) WRITE(numout,*) 'CN ratio is read: ',read_cn  
    !
    !
    !Config Key   = STOMATE_OK_DGVM
    !Config Desc  = Activate DGVM?
    !Config If    = OK_STOMATE
    !Config Def   = n
    !Config Help  = set to TRUE if DGVM is to be activated
    !Config Units = [FLAG]
    ok_dgvm = .FALSE.
    CALL getin_p('STOMATE_OK_DGVM',ok_dgvm)
    IF ( ok_dgvm ) THEN
       CALL ipslerr(3,'control_initialize', &
            'All calls to DGVM have been commented. See ticket 504 before retoring the DGVM.',&
            'Remember that stomate now uses only veget_max but the DGVM may require rescaled veget_max (previously veget_cov_max)',&
            'This still needs to be changed before activation STOMATE_OK_DGVM')
    ENDIF

    !Config Key   = CHEMISTRY_BVOC
    !Config Desc  = Activate calculations for BVOC
    !Config If    = OK_SECHIBA
    !Config Def   = n
    !Config Help  = set to TRUE if biogenic emissions calculation is to be activated
    !Config Units = [FLAG]
    ok_bvoc = .FALSE.
    CALL getin_p('CHEMISTRY_BVOC', ok_bvoc)
    IF (printlev>=1) WRITE(numout,*) 'Biogenic emissions: ', ok_bvoc

    IF ( ok_bvoc ) THEN 
       ok_leafage         = .TRUE. 
       ok_radcanopy       = .TRUE. 
       ok_multilayer      = .TRUE.
       ok_pulse_NOx       = .TRUE.
       ok_bbgfertil_NOx   = .TRUE.
       ok_cropsfertil_NOx = .TRUE.
    ELSE
       ok_leafage         = .FALSE. 
       ok_radcanopy       = .FALSE. 
       ok_multilayer      = .FALSE.
       ok_pulse_NOx       = .FALSE.
       ok_bbgfertil_NOx   = .FALSE.
       ok_cropsfertil_NOx = .FALSE.
    ENDIF
    !
    !Config Key   = CHEMISTRY_LEAFAGE
    !Config Desc  = Activate LEAFAGE?
    !Config If    = CHEMISTRY_BVOC
    !Config Def   = n
    !Config Help  = set to TRUE if biogenic emissions calculation takes leaf age into account
    !Config Units = [FLAG]
    CALL getin_p('CHEMISTRY_LEAFAGE', ok_leafage)
    IF (printlev>=1) WRITE(numout,*) 'Leaf Age: ', ok_leafage
    !
    !Config Key   = CANOPY_EXTINCTION 
    !Config Desc  = Use canopy radiative transfer model?
    !Config If    = CHEMISTRY_BVOC 
    !Config Def   = n
    !Config Help  = set to TRUE if canopy radiative transfer model is used for biogenic emissions 
    !Config Units = [FLAG]
    CALL getin_p('CANOPY_EXTINCTION', ok_radcanopy)
    IF (printlev>=1) WRITE(numout,*) 'Canopy radiative transfer model: ', ok_radcanopy
    !
    !Config Key   = CANOPY_MULTILAYER
    !Config Desc  = Use canopy radiative transfer model with multi-layers
    !Config If    = CANOPY_EXTINCTION 
    !Config Def   = n
    !Config Help  = set to TRUE if canopy radiative transfer model is with 
    !Config         multiple layers. DO NOT CONFUSE with the settings for the
    !Config         multi-layer energy budget.  
    !Config Units = [FLAG]
    CALL getin_p('CANOPY_MULTILAYER', ok_multilayer)
    IF (printlev>=1) WRITE(numout,*) 'Multi-layer Canopy model: ', ok_multilayer
    !
    !Config Key   = NOx_RAIN_PULSE
    !Config Desc  = Calculate NOx emissions with pulse?
    !Config If    = CHEMISTRY_BVOC 
    !Config Def   = n
    !Config Help  = set to TRUE if NOx rain pulse is taken into account
    !Config Units = [FLAG]
    CALL getin_p('NOx_RAIN_PULSE', ok_pulse_NOx)
    IF (printlev>=1) WRITE(numout,*) 'Rain NOx pulsing: ', ok_pulse_NOx
    !
    !Config Key   = NOx_BBG_FERTIL
    !Config Desc  = Calculate NOx emissions with bbg fertilizing effect?
    !Config If    = CHEMISTRY_BVOC 
    !Config Def   = n
    !Config Help  = set to TRUE if NOx emissions are calculated with bbg effect 
    !Config         Fertil effect of bbg on NOx soil emissions 
    !Config Units = [FLAG]
    CALL getin_p('NOx_BBG_FERTIL', ok_bbgfertil_NOx)
    IF (printlev>=1) WRITE(numout,*) 'NOx bbg fertil effect: ', ok_bbgfertil_NOx
    !
    !Config Key   = NOx_FERTILIZERS_USE
    !Config Desc  = Calculate NOx emissions with fertilizers use?
    !Config If    = CHEMISTRY_BVOC 
    !Config Def   = n
    !Config Help  = set to TRUE if NOx emissions are calculated with fertilizers use
    !Config         Fertilizers use effect on NOx soil emissions  
    !Config Units = [FLAG] 
    !
    CALL getin_p('NOx_FERTILIZERS_USE', ok_cropsfertil_NOx)
    IF (printlev>=1) WRITE(numout,*) 'NOx Fertilizers use: ', ok_cropsfertil_NOx
    !Config Key  = Is CO2 impact on BVOC accounted for using Possell 2005 ?
    !Config Desc = In this case we use Possell 2005 parameterisation 
    !Config Desc = to take into account the impact of CO2 on biogenic emissions for 
    !Config Desc = isoprene 
    !Config Def  = n 
    !Config Help = set to TRUE if Possell parameterisation has to be considered for the CO2 impact
    !
    ok_co2bvoc_poss = .FALSE.
    CALL getin_p('CO2_FOR_BVOC_POSSELL', ok_co2bvoc_poss)
    IF (printlev>=1) WRITE(numout,*) 'CO2 impact on BVOC - Possell parameterisation: ', ok_co2bvoc_poss
    !
    !Config Key  = Is CO2 impact on BVOC accounted for using Wilkinson 2009 ? 
    !Config Desc = In this case we use Wilkinson 2009 parameterisation 
    !Config Desc = to take into account the impact of CO2 on biogenic emissions for 
    !Config Desc = isoprene 
    !Config Def  = FALSE
    !Config Help = set to TRUE if Wilkinson parameterisation has to be considered for the CO2 impact
    !
    ok_co2bvoc_wilk = .FALSE.
    CALL getin_p('CO2_FOR_BVOC_WILKINSON', ok_co2bvoc_wilk)
    IF (printlev>=1) WRITE(numout,*) 'CO2 impact on BVOC - Wilkinson parameterisation: ', ok_co2bvoc_wilk
    !
    !+++CHECK+++
    ! No longer in CN but still in CN.CAN
    !Config Key  = CONSTANT_MORTALITY
    !Config Desc = Use constant prescribed mortality or calculate as a function of last year's NPP
    !Config If   = OK_STOMATE
    !Config Def  = TRUE
    !Config Help = set to TRUE if constant mortality is to be assumed
    ok_constant_mortality = .TRUE.
    CALL getin_p('CONSTANT_MORTALITY',ok_constant_mortality)
    IF (printlev>=1) WRITE(numout,*) 'MORTALITY is assumed to be constant ',&
         '(instead of being a function of vigor): ', &
         ok_constant_mortality
    !+++++++++++
    !
    !Config Key   = OK_SNAGS
    !Config Desc  = puts a part of woody litter input in a snag litter pool
    !Config If    = OK_STOMATE
    !Config Def   = n
    !Config Help  = create snag litter pool. When set to yes, snag
    !Config         pool will be accounted for, with snag decay and snag
    !Config         fall into woody poolThe part of the woody litter
    !Config         input going to snag pool should
    !Config         depend on the nature of the disturbances 
    !Config Units = [FLAG]
    ok_snags = .FALSE.
    CALL getin_p('OK_SNAGS', ok_snags)
    !
    !Config Key   = OK_READ_FM_MAP
    !Config Desc  = Read the forest management strategy from a map  
    !Config If    = OK_STOMATE
    !Config Def   = FALSE
    !Config Help  = A logical flag determining if we read
    !Config         in the forest management strategy from a map (NetCDF file).
    !Config         This should be .TRUE. for all applications
    !Config         except for debugging and pixel-level simulations.  
    !Config         If this option is equal to TRUE, we will overwrite 
    !Config         the forest_managed_forced option, so you should be 
    !Config         careful to only use one or the other.
    !Config Units = [FLAG]
    ok_read_fm_map=.FALSE.
    CALL getin_p('OK_READ_FM_MAP',ok_read_fm_map)
    !
    !Config Key   = CONCEPT_SCALE
    !Config Desc  = Conceptual scale of the simulations for scale-dependent processes.
    !Config If    = OK_STOMATE
    !Config Def   = global 
    !Config Help  = Forest management, wind storms and bark beetle outbreak are
    !               scale dependent. Clear cuts, for example, 
    !               are common at the site-level but are unlikley at the pixel
    !               level of a 50 x 50 km2 pixel as that would imply a 2500 km2
    !               clear cut. Set the conceptual scale for forest management to:
    !               "hectare", "europe", or "global". See constantes_mtc for
    !                more details on the difference between these cases.
    !Config Units = [FLAG]
    concept_scale = 'global'
    CALL getin_p('CONCEPT_SCALE',concept_scale)
    !
    !Config Key   = FORCING_RESOLUTION
    !Config Desc  = Spatial and temporal resolution of the climate forcing
    !Config If    = OK_WINDTHROW, OK_PEST
    !Config Def   = 2degrees_6hourly 
    !Config Help  = Forest management, wind storms and bark beetle outbreaks are
    !               scale and thus forcing dependent. At the ha scale a whole forest
    !               could be destroyed by a wind storm. At a scale of 2500km2 this
    !               is highly unlikely. A handfull of storm and pest parameters
    !               depends on the temporal and spatial scale of the climate forcing.
    !               In specific applications or when testing the model, the concept_scale
    !               and forcing resolution could differ. For most applications
    !               including CMIP simulations, they should be line with each other.
    !               Note that in LMDzOR coarse spatial scales can be combined with
    !               fine temporal scales. This may require specific settings.
    !               Possible settings: "2degrees_6hourly", "1degree_daily", 0.5degree_6hourly",
    !               "8km_hourly","site", "LMDzOR_zoomed", and "LMDzOR".
    !               Use CRU_JRA parameters as the default. This setting should not
    !               affect the simulations unless OK_WINDSTORM, OK_PEST or both are set
    !               to .TRUE.
    !Config Units = [FLAG]
    forcing_resolution='2degrees_6hourly'
    CALL getin_p('FORCING_RESOLUTION',forcing_resolution)
    !
    !Config Key   = OK_READ_SP_CLEARCUT_MAP
    !Config Desc  = Read a map prescribing whether a pixel and PFT gets clearcut during spinup  
    !Config If    = OK_STOMATE
    !Config Def   = FALSE
    !Config Help  = We need this option to create the spatial heterogeneity 
    !               in forest age during spinup to mimic stochastic occurrence
    !               natural disturbance events that lead to complete forest
    !               regneration. This is also to break up the synchrony in 
    !               biomass growth during spinup for different pixels.
    !               NOTE: thechnically this feature works fine but conceptually
    !               it is still problematic. After the spinup the model restores
    !               equilibrium - this method was replaced by a post-processing
    !               approach we called "nudged spinup".
    !Config Units = [FLAG]
    ok_read_sp_clearcut_map = .FALSE.
    CALL getin_p('OK_READ_SP_CLEARCUT_MAP',ok_read_sp_clearcut_map)

    !Config Key   = OK_CHANGE_SPECIES
    !Config Desc  = Change species after a stand replacing disturbance
    !Config If    = OK_STOMATE
    !Config Def   = FALSE
    !Config Help  = Sometimes it's a good idea to change species 
    !               after a clearcut on managed forest.  If this 
    !               flag is true, we do this. If not, the same PFT 
    !               is always replanted.
    !Config Units = [FLAG]
    ok_change_species=.FALSE.
    CALL getin_p('OK_CHANGE_SPECIES',ok_change_species)
    
    !Config Key   = READ_SPECIES_CHANGE_MAP
    !Config Desc  = Read the new tree species from a species map
    !Config If    = OK_STOMATE
    !Config Def   = FALSE
    !Config Help  = A logical flag determining if we read
    !Config         in a map which changes species after a clearcut
    !Config         To be used with ok_change_species = .TRUE. or
    !Config         set species_change_force (see below).
    !Config Units = [FLAG]
    ok_read_species_change_map=.FALSE.
    CALL getin_p('OK_READ_SPECIES_CHANGE_MAP',ok_read_species_change_map)
    
    !Config Key   = OK_READ_DESIRED_FM_MAP
    !Config Desc  = Read the new FM strategu from a map
    !Config If    = OK_STOMATE, OK_CHANGE_SPECIES
    !Config Def   = FALSE
    !Config Help  = If we change the forest management after a clearcut, do 
    !               we want to read the new FM strategy from a map?
    !Config Units = [FLAG]
    ok_read_desired_fm_map=.FALSE.
    CALL getin_p('OK_READ_DESIRED_FM_MAP',ok_read_desired_fm_map)
    
    !Config Key   = OK_LITTER_RAKING
    !Config Desc  = Activite litter raking
    !Config If    = OK_STOMATE
    !Config Def   = FALSE
    !Config Help  = Check to see if we are interested in using a litter 
    !               demand map to remove litter from forest PFTs and put 
    !               it into agricultural PFTs. This simulations the 
    !               practice of litter raking.
    !               If TRUE, this flag will simulate litter raking in
    !               in grid squares.  This has the effect of moving litter
    !               once a year from forest PFTs to agricultural PFTs, if they
    !               are present on this pixel.  If TRUE, you must also provide
    !               a map with the litter demand so we know how much litter
    !               to remove for each pixel. Litter raking is a historical
    !               land use so you reconstructions to use this option.
    !Config Units = [FLAG]
    ok_litter_raking=.FALSE.
    CALL getin_p('OK_LITTER_RAKING',ok_litter_raking)
       
    !Config Key   = OK_DIMENSIONAL_PRODUCT_USE
    !Config Desc  = Product pools are based on the dimensions of the harvest
    !Config If    = OK_STOMATE
    !Config Def   = TRUE
    !Config Help  = Once the wood is harvested (through management or
    !               LCC) it ends up in wood product pools. Two options were 
    !               implemeted: (1) the product use and the longevity, of the
    !               product pools depend on the dimensions of the harvest. (2) the 
    !               dimensions are ignored and the wood is used according
    !               to fixed ratios.
    !Config Units = [FLAG]
    ok_dimensional_product_use=.TRUE.
    CALL getin_p('OK_DIMENSIONAL_PRODUCT_USE',ok_dimensional_product_use)

    !Config Key   = FORCED_CLEAR_CUT
    !Config Desc  = Use to force a clear cut at a specific year during a simulation.
    !Config If    = OK_STOMATE
    !Config Def   = FALSE
    !Config Help  = Use to force a clear cut at a specific year during a simulation.
    !               This parameter is used by the ENSEMBLE runs to ensure
    !               that the age of the simulated forest matches the age of
    !               the observations
    !Config Units = year
    forced_clear_cut= .FALSE.
    CALL getin_p('FORCED_CLEAR_CUT',forced_clear_cut)

    ! Check for consistency
    IF (ok_change_species .AND. ok_litter_raking)THEN

       ! There is no obvious conflict in combining species change
       ! and litter raking but it was never tested. Better not to
       ! combine these settings
       WRITE(numout,*) 'ERROR - conflicting settings in run.def'
       WRITE(numout,*) 'Trying to jointly use litter raking and species change'
       WRITE(numout,*) 'This should be tested first'
       CALL ipslerr(3,'ERROR: run.def',&
            'ERROR - Trying to use litter raking',&
            'at the same time of species change',&
            'The code was not developed or tested for this combination')       
    ENDIF
   
    
    IF ( ok_dgvm ) THEN

       ok_stomate = .TRUE.
       CALL ipslerr(2,'control_initialize', &
            'If you want to use the DGVM you need STOMATE',&
            'activating the flag OK_STOMATE',&
            'We set OK_STOMATE to TRUE to ensure consistency')
    ENDIF

    IF ( ok_dgvm ) ok_stomate = .TRUE.

    IF ( ok_multilayer .AND. .NOT.(ok_radcanopy) ) THEN
       ! Note that radcanopy is the radiation transfer calculated in the chemistry module
       ! it needs to be checked and confirmed that this radiation transfer approach
       ! is consistent with the radiation transfer in the albedo, transmission and
       ! absortion calculation in other parts of the model.
       ok_radcanopy  = .TRUE.
       IF (printlev>=1) WRITE(numout,*) 'You want to use the multilayer model without activating the flag CANOPY_EXTINCTION'
       IF (printlev>=1) WRITE(numout,*) 'We set CANOPY_EXTINCTION to TRUE to ensure consistency'
    ENDIF

    !
    !Config Key   = OK_C13
    !Config Desc  = Calculate C13 fractionation
    !Config If    = OK_SECHIBA 
    !Config Def   = FALSE
    !Config Help  = set to TRUE if C13 is to be calculated. C13 is
    !               calculated for leaf photosynthesis only.
    !Config Units = [FLAG]
    ok_c13 = .FALSE.
    CALL getin_p('OK_C13', ok_c13)
    IF (printlev>=1) WRITE(numout,*) 'C13 fractionation is calculated: ', ok_c13

!    IF(ok_c13) THEN
!       CALL ipslerr(3,'control.f90 - ok_c13',&
!            'There is a problem with leaf_ci recalculation after water stress',&
!            'It is not recommended to use c13 simulation until fixing that.','')   
!    ENDIF

    ! 
    !Config Key  = OK_WINDTHROW 
    !Config Desc = Activate windthrow 
    !Config If   = OK_STOMATE
    !Config Def  = FALSE
    !Config Help = Set to TRUE if storm damage needs to be accounted for. Calculates
    !              PFT-specific critical wind speeds and subsequent tree mortality from
    !              storm damage (stem breakage and uprooting)
    !Config Units = [FLAG] 
    ok_windthrow = .FALSE.
    CALL getin_p('OK_WINDTHROW',ok_windthrow) 
    IF (printlev>=1) WRITE(numout,*) 'Windthrow is activated: ',ok_windthrow
    ! 
    !Config Key  = OK_PEST 
    !Config Desc = Calculate pest outbreaks.
    !Config Def  = FALSE
    !Config If   = OK_STOMATE
    !Config Help = Set to TRUE if pest outbreaks need to be accounted for. For the moment
    !              the only pest ORCHIDEE can account for are bark beetle outbreaks. The
    !              parameters are biased toward bark beetles of Norway spruce. This is a
    !              a general flag. The pft-specific parameter beetle_pft controls whether
    !              tree mortality from beetle attacks needs to be calculated.
    !Config Units = [FLAG] 
    ok_pest = .FALSE.
    CALL getin_p('OK_PEST',ok_pest)
    IF (printlev>=1) WRITE(numout,*) 'Pest outbreak is activated: ',ok_pest

    !Config Key  = OK_PHENO 
    !Config Desc = Calculate lai and phenology.
    !Config Def  = TRUE
    !Config If   = OK_STOMATE
    !Config Help = This flag is used to switch between a calculated and prescribed lai
    !              It seems that with the new way LAI is prescribed, the flag is no 
    !              longer needed. Clean when prescribed LAI is working again.
    !Config Units = [FLAG]
    ok_pheno = .TRUE.
    CALL getin_p('OK_PHENO',ok_pheno)
    IF (printlev>=1) WRITE(numout,*) 'Phenology is activated: ',ok_pheno

    !Config Key  = USE_FLUXNET
    !Config Desc = Deactivate retrieting wind speed in stomate, gustiness in
    ! windthrow
    !Config If   = OK_STOMATE
    !Config Def  = FALSE
    !Config Help = Set to TRUE if the fluxnet forcing is used.
    !              Temperal flag for now, to separate test cases
    !Config Units = [FLAG] 
    use_fluxnet = .FALSE.
    CALL getin_p('USE_FLUXNET',use_fluxnet)
    WRITE(numout,*) 'use_fluxnet is activated: ',use_fluxnet

    !Config Key  = CAL_DEFLECTION
    !Config Desc = Calculate deflection load
    !Config If   = OK_STOMATE
    !Config Def  = FALSE
    !Config Help = This flag activates dynamic deflection load in windthrow module
    !Config Units = [FLAG] 
    cal_deflection = .FALSE.
    CALL getin_p('CAL_DEFLECTION',cal_deflection)

    !Config Key  = WIND_CAL_PFT
    !Config Desc = Calculate pft level wind speed
    !Config If   = OK_STOMATE
    !Config Def  = FALSE
    !Config Help = Enables computation of wind speed separately for each PFT, 
    !               based on its canopy height rather than pixel-level wind speed.
    !Config Units = [FLAG] 
    wind_cal_pft = .TRUE.
    CALL getin_p('WIND_CAL_PFT',wind_cal_pft)
    WRITE(numout,*) 'wind_cal_pft is activated: ',wind_cal_pft

    !Config Key  = TEST_SPACING
    !Config Desc = Manipulate spacing in windthrow
    !Config If   = OK_STOMATE
    !Config Def  = FALSE
    !Config Help = When TRUE, manipulates circ_class_n in the windthrow module 
    !               to increase tree density when trees are small.
    !Config Units = [FLAG] 
    test_spacing = .FALSE.
    CALL getin_p('TEST_SPACING',test_spacing)
    WRITE(numout,*) 'test_spacing is activated: ',test_spacing

!
    ! Configuration : number of PFTs and parameters
    !

    ! 1. Number of PFTs defined by the user

    !Config Key   = NVM
    !Config Desc  = number of PFTs  
    !Config If    = OK_SECHIBA or OK_STOMATE
    !Config Def   = 13
    !Config Help  = The number of vegetation types define by the user
    !Config Units = [-]
    nvm = 15
    CALL getin_p('NVM',nvm)
    IF (printlev>=1) WRITE(numout,*) 'The number of pfts used by the model is : ', nvm

    !
    ! TESTS
    !
    ! Tests are (final) solutions to a problem. Following testing,
    ! the code under the test_xxx flag can be accepted as final
    ! and the previous code can be removed.

    !Config Key   = TEST_DYNAMIC_ALPHA_SELF_THIN
    !Config Desc  = Use a dynamic parameter for alpha_self_thin
    !Config If    = -
    !Config Def   = n
    !Config Help  = 
    !Config Units = [FLAG]
    test_dynamic_alpha_self_thin = .TRUE.
    CALL getin_p('TEST_DYNAMIC_ALPHA_SELF_THIN',test_dynamic_alpha_self_thin)
    
    !Config Key   = TEST_Q10_ON_PLANT_UPTAKE
    !Config Desc  = Use a Q10 function to control the plant uptake
    !Config If    = -
    !Config Def   = n
    !Config Help  = 
    !Config Units = [FLAG]
    test_Q10_on_plant_uptake = .FALSE.
    CALL getin_p('TEST_Q10_ON_PLANT_UPTAKE',test_Q10_on_plant_uptake)

    !Config Key   = TEST_ADJUST_RHETERO
    !Config Desc  = Adjust heterotrophic respiration when mineralisation becomes negative
    !Config Def   = y
    !Config Help  = 
    !Config Units = [FLAG]
    test_adjust_rhetero = .FALSE.
    CALL getin_p('TEST_ADJUST_RHETERO',test_adjust_rhetero)

    !Config Key   = TEST_FDI_TYPES
    !Config Desc  = Use different Fire danger index for SPITFIRE
    !Config If    = -
    !Config Def   = n
    !Config Help  = 
    !Config Units = [FLAG]
    test_fdi_types = .TRUE.
    CALL getin_p('TEST_FDI_TYPES',test_fdi_types)
    

    ! 2. Initialize vertical discretization
    ! All initialization is done in the vertical module
    ! Calculate ngrnd and nslm
    CALL vertical_soil_init

    ! 3. Allocate and intialize the pft parameters
    CALL pft_parameters_main()
 
    ! 4. Allocate and intialize dynamic parameters
    IF (printlev>=3) THEN
       WRITE(numout,*) 'Entering dynamic parameters main'
       CALL flush(numout)
    END IF
    
    CALL dynamic_parameters_main(kjpindex)

    ! 5. Activation sub-models of ORCHIDEE
    IF (printlev>=3) THEN
       WRITE(numout,*) 'Entering activate_sub_models'
       CALL flush(numout)
    END IF

    CALL activate_sub_models()

    ! 6. Vegetation configuration
    IF (printlev>=3) THEN
       WRITE(numout,*) 'Entering veget_config'
       CALL flush(numout)
    END IF
    
    CALL veget_config

    ! 6. Read the parameters in the run.def file  according the flags
    IF (printlev>=3) THEN
       WRITE(numout,*)'In control_initialize: call config_pft_parameters'
       CALL flush(numout)
    END IF
    CALL config_pft_parameters
 
    IF ( ok_sechiba ) THEN
       IF (printlev>=3) THEN
          WRITE(numout,*)'In control_initialize: call config_sechiba_parameters'
          CALL flush(numout)
       END IF
       CALL config_sechiba_parameters

       IF (printlev>=3) THEN
          WRITE(numout,*)'In control_initialize: call config_sechiba_pft_parameters'
          CALL flush(numout)
       END IF
       CALL config_sechiba_pft_parameters()

       IF (printlev>=3) THEN
          WRITE(numout,*)'In control_initialize: call config_sechiba_dynamic_parameters'
          CALL flush(numout)
       END IF
       CALL dynamic_parameters_config_sechiba()
       
    END IF


    !! Initialize variables in constantes_soil
    IF (printlev>=3) THEN
       WRITE(numout,*)'In control_initialize: call config_soil_parameters'
       CALL flush(numout)
    END IF
    CALL config_soil_parameters()


    !! Coherence check for depth of thermosoil for long term simulation where soil thermal inertia matters
    !! ok_freeze_thermix is defined in config_soil_parameters
    IF (ok_freeze_thermix .AND. zmaxt < 11) THEN
       WRITE(numout,*) 'ERROR : Incoherence between ok_freeze_thermix activated and soil depth too small. '
       WRITE(numout,*) 'Here a soil depth of ', zmaxt, 'm is used for the soil thermodynamics'
       WRITE(numout,*) 'Set DEPTH_MAX_T=11 or higher in run.def parameter file or deactivate soil freezing'
       CALL ipslerr_p(3,'control_initialize','Too shallow soil chosen for the thermodynamic for soil freezing', &
            'Adapt run.def with at least DEPTH_MAX=11','')
    END IF
        
    IF (printlev>=3) THEN
       WRITE(numout,*)'In control_initialize: call config_co2_parameters'
       CALL flush(numout)
    END IF
    CALL config_co2_parameters
    
        
    IF ( ok_stomate ) THEN
       IF (printlev>=3) THEN
          WRITE(numout,*)'In control_initialize: call config_stomate_parameters'
          CALL flush(numout)
       END IF
       CALL config_stomate_parameters
       
       IF (printlev>=3) THEN
          WRITE(numout,*)'In control_initialize: call config_stomate_pft_parameters'
          CALL flush(numout)
       END IF
       CALL config_stomate_pft_parameters

       ! AGE_CLASS_BOUNDS_PFT (Marie 2026, design/MODULE_DESIGN_AGE_CLASS_BOUNDS_PFT.md).
       ! IMPERATIVEMENT ICI, apres config_stomate_pft_parameters : c'est elle qui lit
       ! LARGEST_TREE_DIA au namelist. Deriver les bornes plus tot (dans
       ! config_stomate_parameters, ou age_class_bound est alloue) n'aurait vu que les
       ! defauts MTC et aurait ignore sans message le 0.20 m de l'eucalyptus.
       ! No-op si OK_AGE_CLASS_BOUND_PFT=n.
       IF (printlev>=3) THEN
          WRITE(numout,*)'In control_initialize: call derive_age_class_bounds'
          CALL flush(numout)
       END IF
       CALL derive_age_class_bounds()

       IF (printlev>=3) THEN
          WRITE(numout,*)'In control_initialize: call config_stomate_dynamic_parameters'
          CALL flush(numout)
       END IF
       CALL dynamic_parameters_config_stomate()
    END IF
    
    IF ( ok_dgvm ) THEN
       IF (printlev>=3) THEN
          WRITE(numout,*)'In control_initialize: call config_dgvm_parameters'
          CALL flush(numout)
       END IF
       CALL config_dgvm_parameters
    END IF

  END SUBROUTINE control_initialize

END MODULE control
