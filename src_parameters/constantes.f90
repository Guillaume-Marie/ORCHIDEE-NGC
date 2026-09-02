! =================================================================================================================================
! MODULE       : constantes
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        "constantes" module contains subroutines to initialize most of the exernalized parameters. This module
!!              also make a use to the module constantes_var where the parameters are declared.
!!
!!\n DESCRIPTION: This module contains subroutines to initialize most of the exernalized parameters. This module
!!                also make a use to the module constantes_var where the parameters are declared.\n
!!                This module can be used to acces the subroutines and the constantes. The constantes declarations
!!                can also be used seperatly with "USE constantes_var".
!!
!! RECENT CHANGE(S): Didier Solyga : This module contains now all the externalized parameters of ORCHIDEE 
!!                   listed by modules which are not pft-dependent  
!!                   Josefine Ghattas 2013 : The declaration part has been extracted and moved to module constates_var
!!
!! REFERENCE(S)	: 
!! - Louis, Jean-Francois (1979), A parametric model of vertical eddy fluxes in the atmosphere. 
!! Boundary Layer Meteorology, 187-202.
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_parameters/constantes.f90 $
!! $Date: 2026-04-30 09:10:24 +0200 (jeu. 30 avril 2026) $
!! $Revision: 9507 $
!! \n
!_ ================================================================================================================================

MODULE constantes

  USE constantes_var
  USE pft_parameters_var, ONLY : nvm   ! AGE_CLASS_BOUNDS_PFT : age_class_bound est (nagec,nvm)
  USE defprec
  USE ioipsl_para, ONLY : getin_p, ipslerr_p
  USE mod_orchidee_para_var, ONLY : numout
  USE time, ONLY : one_day, dt_sechiba
  USE constantes_soil_var, ONLY: peat_bulk_density
  USE vertical_soil_var, ONLY: ngrnd

  IMPLICIT NONE
  
CONTAINS


!! ================================================================================================================================
!! SUBROUTINE   : activate_sub_models
!!
!>\BRIEF         This subroutine reads the flags in the configuration file to
!! activate some sub-models like routing, irrigation, fire, herbivory, ...  
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

  SUBROUTINE activate_sub_models()

    IMPLICIT NONE

    !! 0. Variables and parameters declaration

    !! 0.4 Local variables

    !_ ================================================================================================================================

     IF (ok_stomate) THEN

       !Config Key   = HERBIVORES
       !Config Desc  = herbivores allowed?
       !Config If    = OK_STOMATE 
       !Config Def   = n
       !Config Help  = With this variable, you can determine
       !Config         if herbivores are activated
       !Config Units = [FLAG]
       ok_herbivores = .FALSE.
       CALL getin_p('HERBIVORES', ok_herbivores)
       !
       !Config Key   = TREAT_EXPANSION
       !Config Desc  = treat expansion of PFTs across a grid cell?
       !Config If    = OK_STOMATE 
       !Config Def   = n
       !Config Help  = With this variable, you can determine
       !Config         whether we treat expansion of PFTs across a
       !Config         grid cell.
       !Config Units = [FLAG]
       treat_expansion = .FALSE.
       CALL getin_p('TREAT_EXPANSION', treat_expansion)

       !Config Key   = SLA_DYN
       !Config Desc  = Account for a dynamic SLA
       !Config Def   = n
       !Config if    = OK_STOMATE
       !Config Help  = If this flag is set to true (y) then the SLA
       !Config         is computed dynamically, varying with leaf biomass
       !Config Units = [FLAG]
       sla_dyn = .TRUE.
       CALL getin_p('SLA_DYN',sla_dyn)

       !Config Key   = LPJ_GAP_CONST_MORT
       !Config Desc  = Constant mortality
       !Config If    = OK_STOMATE AND NOT OK_DGVM
       !Config Def   = y/n depending on OK_DGVM
       !Config Help  = set to TRUE if constant mortality is to be activated
       !Config         
       !Config Units = [FLAG]

       ! Set Default value different if DGVM is activated.
       IF ( ok_dgvm ) THEN
          lpj_gap_const_mort=.FALSE.
       ELSE
          lpj_gap_const_mort=.TRUE.
       END IF
       CALL getin_p('LPJ_GAP_CONST_MORT', lpj_gap_const_mort)

       IF (ok_dgvm .AND. lpj_gap_const_mort) THEN
          CALL ipslerr_p(1,"activate_sub_models","Both OK_DGVM and LPJ_GAP_CONST_MORT are activated.",&
               "This combination is possible but unusual","The simulation will continue with these flags activated." )
       ELSEIF (.NOT. ok_dgvm  .AND. .NOT. lpj_gap_const_mort) THEN
           CALL ipslerr_p(3,"activate_sub_models", &
                "The combination of OK_DGVM=false and LPJ_GAP_CONST_MORT=false is not operational in this version", &
                "Some parts of the code should first be revised.","" )
       END IF

       !Config Key   = HARVEST_AGRI
       !Config Desc  = Harvest model for agricultural PFTs.
       !Config If    = OK_STOMATE 
       !Config Def   = y
       !Config Help  = Compute harvest above ground biomass for agriculture.
       !Config         Change daily turnover.
       !Config Units = [FLAG]
       harvest_agri = .TRUE. 
       CALL getin_p('HARVEST_AGRI', harvest_agri)
       !
       !Config Key   = OK_SPITFIRE
       !Config Desc  = Activate fire using SPITFIRE
       !Config If    = OK_STOMATE 
       !Config Def   = n
       !Config Help  = 
       !Config Units = [FLAG]
       ok_spitfire = .FALSE.
       CALL getin_p('OK_SPITFIRE', ok_spitfire)

       !Config Key   = OK_LIGHTNINGFIRE
       !Config Desc  = Activate lightning fire using SPITFIRE
       !Config If    = OK_SPITFIRE
       !Config Def   = y
       !Config Help  =
       !Config Units = [FLAG]
       ok_lightningfire = .TRUE.
       CALL getin_p('OK_LIGHTNINGFIRE', ok_lightningfire)


       !Config Key   = SPINUP_ANALYTIC
       !Config Desc  = Activation of the analytic resolution of the spinup.
       !Config If    = OK_STOMATE
       !Config Def   = n
       !Config Help  = Activate this option if you want to solve the spinup by the Gauss-Jordan method.
       !Config Units = BOOLEAN
       spinup_analytic = .FALSE.
       CALL getin_p('SPINUP_ANALYTIC',spinup_analytic)

       
       !Config Key   = CALCULATE_GPP_PREIND
       !Config Desc  = Calculate preindustial longterm GPP
       !Config If    = OK_STOMATE
       !Config Def   = y if SPINUP_ANALYTIC else n
       !Config Help  = Activate this flag during the spinup during the preindustrial period
       !Config Units = BOOLEAN
       IF (spinup_analytic) THEN
          calculate_gpp_preind = .TRUE.
       ELSE
          calculate_gpp_preind = .FALSE.
       END IF
       CALL getin_p('CALCULATE_GPP_PREIND',calculate_gpp_preind)
       
       
       !Config Key   = OK_DYN_GRASS_DEN
       !Config Desc  = Use continuous dynamic grassland density
       !Config If    = OK_STOMATE
       !Config Def   = n
       !Config Help  = Calculate grassland density (plants per ha) as an outcome of plant growth and
       !               resource competition. If the vegetation acquires insufficient resources to
       !               keep its resreve and labile pools up, the density will decrease. If the
       !               vegetation has plenty of resources, the density will increase.
       !Config Units = BOOLEAN
       ok_dyn_grass_den=.FALSE.
       CALL getin_p('OK_DYN_GRASS_DEN',ok_dyn_grass_den)

     ENDIF

     !
     ! HACKS TO FIX AND REMOVE
     !
     ! Hacks are temporary solutions. Often quick and dirty that
     ! need attention before the problem could be considered to
     ! be solved.
     ! 
     !Config Key   = HACK_ENERBIL_HYDROL
     !Config Desc  = Flag to skip a particular block of code in mleb.f90 
     !Config If    = -
     !Config Def   = n
     !Config Help  = For debugging only! Flag to skip a particular block of code in mleb.f90 
     !               which results in incorrect results for large scale simulations.
     !Config Units = [FLAG]
     hack_enerbil_hydrol = .FALSE. 
     CALL getin_p('HACK_ENERBIL_HYDROL',hack_enerbil_hydrol)
    
     !Config Key   = REFILLING_CAPA
     !Config Desc  = defines which fraction of embolism is considered to kill sapwood
     !Config If    = -
     !Config Def   = 1.0
     !Config Help  = The complementary fraction is thus a costless refilling of
     !Config Help  = vessels
     !Config Units = unitless
     refilling_capa = 1.0
     CALL getin_p('REFILLING_CAPA',refilling_capa)

     !Config Key   = RAIN_FAC
     !Config Desc  = remove rain or lower, or increase
     !Config If    = -
     !Config Def   = 1
     !Config Help  = -
     !Config Units = unitless
     rain_fac = 1
     CALL getin_p('RAIN_FAC',rain_fac)
 
     !Config Key   = HACK_VESSEL_LOSS
     !Config Desc  = constant vessel_loss in hydraulic_architecture
     !Config If    = OK_VESSEL_MORTALITY
     !Config Def   = -9999
     !Config Help  = When set outside the range from 0 to 1 it will not be used.
     !Config         Impose a constant vessel_loss in hydraulic_rachitecture
     !Config         when ok_vessel_mortality is TRUE. When set outside the
     !Config         range from 0 to 1 it will not be used.
     !Config Units = unitless
     hack_vessel_loss = -9999 
     CALL getin_p('HACK_VESSEL_LOSS',hack_vessel_loss)
     
     !Config Key   = HACK_KB_M1
     !Config Desc  = Use higher threshold for small lai in the calculation of kB_m1 in condveg
     !Config If    = -
     !Config Def   = n
     !Config Help  = 
     !Config Units = [FLAG]
     hack_kb_m1 = .FALSE. 
     CALL getin_p('HACK_KB_M1',hack_kb_m1)

     !
     ! HACKS TO KEEP
     !
     ! These are hacks that do not fix an problem but that
     ! enhance the functionality of the code when trusting
     ! or debugging.
     !
     !Config Key   = TRUSTING_HACK_AGE_CLASS
     !Config Desc  = Skip the calculation of redistributing biomass between age classes
     !Config If    = -
     !Config Def   = n
     !Config Help  = In age_class_distr their is a calculation that redistributes the 
     !               biomass between age classes. Although mass is conserved during this
     !               calculation it results in precision errors that propagate through
     !               the model calculations. In the trusting we have one test case that 
     !               helps to check the technical integrity of the age class distribution
     !               code: if the model is run without LCC and disturbances, A simulation
     !               without age classes should generate the exact same results as one with 
     !               age classes (but the dimensions of these simulations differ). If 
     !               biomass is recalculated this test fails. With this hack, the biomass is
     !               explicitly moved from one age class to another without any calculations
     !               and thus without any precision errors. With this hack set to .TRUE. 
     !               the model can be trusted for the age class code. This hack
     !               should not be used in simulation experiments as biomass
     !               from different age classes will NOT be merged but simply
     !               overwritten. If all goes well this should result in a mass
     !               balance error.
     !Config Units = [FLAG]
     trusting_hack_age_class = .FALSE.
     CALL getin_p('TRUSTING_HACK_AGE_CLASS',trusting_hack_age_class)
     
     !Config Key   = HACK_PGAP
     !Config Desc  = Flag to use Lambert Beer instead of Pgap 
     !Config If    = -
     !Config Def   = n
     !Config Help  = Flag to use Lambert Beer instead of Pgap to calculate veget. 
     !               Only for debugging as it may introduce inconsistencies.
     !Config Units = [FLAG]
     hack_pgap = .FALSE.
     CALL getin_p('HACK_PGAP',hack_pgap)

     !Config Key   = HACK_VEGET_MAX_NEW
     !Config Desc  = Prescribe VEGET_MAX_NEW rather than reading from a map 
     !Config If    = -
     !Config Def   = n
     !Config Help  = Read veget_max_new from run.def. This is intended to 
     !               be used in combination with a restart file. It only 
     !               works once per simulation. It can be used to debug 
     !               simplified LCC test cases. See slowproc.f90 for more
     !               details. Only for debugging as it will introduce 
     !               inconsistencies.
     !Config Units = [FLAG]
     hack_veget_max_new = .FALSE.
     CALL getin_p('HACK_VEGET_MAX_NEW',hack_veget_max_new) 
      
     
     !
     ! Check consistency (see later)
     !
!!$        IF(.NOT.(ok_routing) .AND. (doirrigation .OR. dofloodplains)) THEN
!!$           CALL ipslerr_p(2,'activate_sub_models', &
!!$               &     'Problem :you tried to activate the irrigation and floodplains without activating the routing',&
!!$               &     'Are you sure ?', &
!!$               &     '(check your parameters).')
!!$        ENDIF

!!$        IF(.NOT.(ok_stomate) .AND. (ok_herbivores .OR. treat_expansion .OR. lpj_gap_const_mort &
!!$            & .OR. harvest_agri)) THEN
!!$          CALL ipslerr_p(2,'activate_sub_models', &
!!$               &     'Problem : try to activate the following options : herbivory, treat_expansion, fire,',&
!!$               &     'harvest_agri and constant mortality without stomate activated.',&
!!$               &     '(check your parameters).')
!!$        ENDIF


  END SUBROUTINE activate_sub_models

!! ================================================================================================================================
!! SUBROUTINE   : veget_config
!!
!>\BRIEF         This subroutine reads the flags controlling the configuration for
!! the vegetation : impose_veg, veget_mpa, lai_map, etc...       
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCE(S) :
!!
!! FLOWCHART    :
!! \n
!_ ================================================================================================================================

  SUBROUTINE veget_config

    IMPLICIT NONE

    !! 0. Variables and parameters declaration

    !! 0.4 Local variables  
    CHARACTER(LEN=30)          :: veget_str         !! update frequency for landuse    
    INTEGER                    :: l
    INTEGER(i_std)             :: ier               !! Check errors in allocation
    INTEGER(i_std)             :: ilevel            !! Allocation indice

    !_ ================================================================================================================================


    !Config Key   = AGRICULTURE
    !Config Desc  = agriculture allowed?
    !Config If    = OK_SECHIBA or OK_STOMATE
    !Config Def   = y
    !Config Help  = With this variable, you can determine
    !Config         whether agriculture is allowed
    !Config Units = [FLAG]
    agriculture = .TRUE.
    CALL getin_p('AGRICULTURE', agriculture)
    !
    !Config Key   = IMPOSE_VEG
    !Config Desc  = Should the vegetation be prescribed ?
    !Config If    = OK_SECHIBA or OK_STOMATE
    !Config Def   = n
    !Config Help  = This flag allows the user to impose a vegetation distribution
    !Config         and its characteristics. It is espacially interesting for 0D
    !Config         simulations. On the globe it does not make too much sense as
    !Config         it imposes the same vegetation everywhere
    !Config Units = [FLAG]
    impveg = .FALSE.
    CALL getin_p('IMPOSE_VEG', impveg)

    !Config Key   = IMPOSE_SOILT
    !Config Desc  = Should the soil type be prescribed ?
    !Config Def   = n
    !Config If    = 
    !Config Help  = This flag allows the user to impose a soil type distribution.
    !Config         It is espacially interesting for 0D
    !Config         simulations. On the globe it does not make too much sense as
    !Config         it imposes the same soil everywhere
    !Config Units = [FLAG]
    impsoilt = .FALSE.
    CALL getin_p('IMPOSE_SOILT', impsoilt)     

    !Config Key   = IMPOSE_SLOPE
    !Config Desc  = Should reinf_slope be prescribed ?
    !Config Def   = n
    !Config If    = 
    !Config Help  = This flag allows the user to impose a uniform fraction of  
    !Config         reinfiltrated surface runoff, with value REINF_SLOPE 
    !Config Units = [FLAG]
    impslope = .FALSE. 
    CALL getin_p('IMPOSE_SLOPE', impslope)

    !Config Key   = IMPOSE_NINPUT_DEP
    !Config Desc  = Should the N inputs from atmospheric deposition be prescribed ?
    !Config Def   = n
    !Config If    = NOT IMPOSE_CN
    !Config Help  = This flag allows the user to impose N inputs from atmospheric deposition 
    !Config         It is espacially interesting for 0D
    !Config         simulations. On the globe it does not make too much sense as
    !Config         it imposes the same N inputs everywhere
    !Config Units = [FLAG]
    impose_ninput_dep = .FALSE.
    CALL getin_p('IMPOSE_NINPUT_DEP', impose_ninput_dep)     

    !Config Key   = IMPOSE_NINPUT_FERT
    !Config Desc  = Should the N inputs from fertilizer be prescribed ?
    !Config Def   = n
    !Config If    = -
    !Config Help  = This flag allows the user to impose N inputs from fertilizer application
    !Config         It is espacially interesting for 0D
    !Config         simulations. On the globe it does not make too much sense as
    !Config         it imposes the same N inputs everywhere
    !Config Units = [FLAG]
    impose_ninput_fert = .FALSE.
    CALL getin_p('IMPOSE_NINPUT_FERT', impose_ninput_fert)

    !Config Key   = IMPOSE_NINPUT_MANURE
    !Config Desc  = Should the N inputs from manure be prescribed ?
    !Config Def   = n
    !Config If    = -
    !Config Help  = This flag allows the user to impose N inputs from manure application
    !Config         It is espacially interesting for 0D
    !Config         simulations. On the globe it does not make too much sense as
    !Config         it imposes the same N inputs everywhere
    !Config Units = [FLAG]
    impose_ninput_manure = .FALSE.
    CALL getin_p('IMPOSE_NINPUT_MANURE', impose_ninput_manure)

    !Config Key   = IMPOSE_NINPUT_BNF
    !Config Desc  = Should the N inputs from biological nitrogen fixation (BNF) be prescribed ?
    !Config Def   = n
    !Config If    = -
    !Config Help  = This flag allows the user to impose N inputs from biological nitrogen fixation (BNF)
    !Config         It is espacially interesting for 0D
    !Config         simulations. On the globe it does not make too much sense as
    !Config         it imposes the same N inputs everywhere
    !Config Units = [FLAG]
    impose_ninput_bnf = .FALSE.
    CALL getin_p('IMPOSE_NINPUT_BNF', impose_ninput_bnf)

    !Config Key   = TEST_GRID
    !Config Desc  = grid cell for which extra output is written to the out_execution file
    !Config If    = OK_STOMATE
    !Config Def   = 1
    !Config Help  = Number of the grid square for which detailed output 
    !Config         If the default value is not one, this can cause crashes in debugging for small regions.
    !Config Units = [-]
    test_grid = 1
    CALL getin_p('TEST_GRID',test_grid)   

    !Config Key   = TEST_PFT
    !Config Desc  = pft for which extra output is written to the out_execution file
    !Config If    = OK_STOMATE
    !Config Def   = 6
    !Config Help  = 
    !Config Units = [-]
    test_pft = 1
    CALL getin_p('TEST_PFT',test_pft)   

    !Config Key   = MIN_CANOPY_HEIGHT_PRESCRIBED
    !Config Desc  = Starting height of the canopy
    !Config If    = READ_LAI=1
    !Config Def   = 10
    !Config Help  = 
    !Config Units = [-]
    min_canopy_height_prescribed = 10.0 
    CALL getin_p('MIN_CANOPY_HEIGHT_PRESCRIBED',min_canopy_height_prescribed)   

    !Config Key   = MAX_CANOPY_HEIGHT_PRESCRIBED
    !Config Desc  = Ending height of the canopy
    !Config If    = READ_LAI=1
    !Config Def   = 50
    !Config Help  = 
    !Config Units = [-]
    max_canopy_height_prescribed = 50.0
    CALL getin_p('MAX_CANOPY_HEIGHT_PRESCRIBED',max_canopy_height_prescribed)   

    !Config Key   = CANOPY_DENSITY_PRESCRIBED
    !Config Desc  = Prescribeddensity of the canopy in ind/m^2
    !Config If    = READ_LAI=1
    !Config Def   = 0.15
    !Config Help  = 
    !Config Units = [-]
    canopy_density_prescribed = 0.15  
    CALL getin_p('CANOPY_DENSITY_PRESCRIBED',canopy_density_prescribed)   

    !Config Key   = MAXIMUM_LAI
    !Config Desc  = Maximum LAI over the year
    !Config If    = READ_LAI=1
    !Config Def   = 5
    !Config Help  = 
    !Config Units = [-]
    maximum_lai = 5.0
    CALL getin_p('MAXIMUM_LAI',maximum_lai)   

    !Config Key   = MINIMUM_LAI
    !Config Desc  = Minimum LAI over the year
    !Config If    = READ_LAI=1
    !Config Def   = 0
    !Config Help  = 
    !Config Units = [-]
    minimum_lai = 0.0
    CALL getin_p('MINIMUM_LAI',minimum_lai)   

    !Config Key   = PHENOLOGY_DOY_START
    !Config Desc  = Starting day of the growing of the leaves
    !Config If    = READ_LAI=1
    !Config Def   = 1
    !Config Help  = 
    !Config Units = [-]
    phenology_doy_start = 1
    CALL getin_p('PHENOLOGY_DOY_START',phenology_doy_start)   

    !Config Key   = MOIGDD_ORCHIDEE_2_0
    !Config Desc  = MOIGDD from ORCHIDEE 2.0
    !Config If    = OK_STOMATE
    !Config Def   = .TRUE.
    !Config Help  = 
    !Config Units = [-]
    moigdd_orchidee_2_0 = .FALSE.
    CALL getin_p('MOIGDD_ORCHIDEE_2_0',moigdd_orchidee_2_0)

    !Config Key   = SENESCENCE_DOY_END
    !Config Desc  = Ending day of the dying of the leaves
    !Config If    = READ_LAI=1
    !Config Def   = 365
    !Config Help  = 
    !Config Units = [-]
    senescence_doy_end = 365
    CALL getin_p('SENESCENCE_DOY_END',senescence_doy_end)   

    !Config Key   = PHENOLOGY_DELTA
    !Config Desc  = Delta of doy between start and end of phenology
    !Config If    = READ_LAI=1
    !Config Def   = 30
    !Config Help  = 
    !Config Units = [-]
    phenology_delta = 30
    CALL getin_p('PHENOLOGY_DELTA',phenology_delta)

    !Config Key   = SENESCENCE_DELTA
    !Config Desc  = Delta of doy between start and end of senescence
    !Config If    = READ_LAI=1
    !Config Def   = 30
    !Config Help  = 
    !Config Units = [-]
    senescence_delta = 30 
    CALL getin_p('SENESCENCE_DELTA',senescence_delta)

    !Config Key   = NLAI
    !Config Desc  = Number of photosyntheis canopy levels 
    !Config If    = OK_SECHIBA
    !Config Def   = 10
    !Config Help  = Number of levels in the canopy used in the photosynthesis
    !Config         routine per level dictacted by nlevels. For example, if
    !Config         if nlevels = 2 and nlai = 3, the photosynthesis
    !Config         will be calculated for nlevels_to=2t*3=6 total levels.
    !Config         The number of canopy levels to be used in 
    !Config         per enegy budget canopy level photosynthesis 
    !Config Units = [-]
    nlai = 10
    CALL getin_p("NLAI",nlai)
    IF (printlev>=1) WRITE(numout,'(I4,A)') nlai,' levels are used for each canopy ' &
         // 'level to calculate photosynthesis.'
    
    ! nlevels_tot, total number of levels, is a variable that we'll pass around for convenience, 
    ! in particular when going into a subroutine. Currently, nlevels=1 is used, thus nlevels_tot=nlai.
    ! Note that when using the multi-layer budget nlai needs to be nlai=jnlvls_canopy+1
    nlevels_tot = nlevels * nlai

    IF (printlev>=2) THEN
       DO ilevel=1,nlevels
          WRITE(numout,*) 'Albedo canopy level ',ilevel,&
               ' has a lower boundary at ',z_level(ilevel),' meters.'
       ENDDO
    END IF
    
   
    ALLOCATE(lvllai_prescribed(nlevels_tot),stat=ier)
    IF (ier /= 0) THEN
       CALL ipslerr_p(3,'veget_config','Problem in allocation lvllai_prescribed','','')
    ENDIF


    !Config Key   = LVLLAI_PRESCRIBED
    !Config Desc  = LAI per level
    !Config If    = READ_LAI=1
    !Config Def   = 0.0
    !Config Help  = 
    !Config Units = [-]
    lvllai_prescribed(:) = 0.0
    CALL getin_p('LVLLAI_PRESCRIBED',lvllai_prescribed)   

    !Config Key   = VEGET_UPDATE
    !Config Desc  = Update vegetation frequency: 0Y or 1Y
    !Config If    = 
    !Config Def   = 0Y
    !Config Help  = The veget datas will be update each this time step. Must be 0Y if IMPOSE_VEG=y.
    !Config Units = [years]
    veget_update=0
    WRITE(veget_str,'(a)') '0Y'
    CALL getin_p('VEGET_UPDATE', veget_str)
    l=INDEX(TRIM(veget_str),'Y')
    READ(veget_str(1:(l-1)),"(I2.2)") veget_update

    ! Coherence test : veget_update can only be 0 or 1
    IF (veget_update /= 0 .AND. veget_update /= 1) then
       WRITE(numout,*) "Error in veget_update=", veget_update
       CALL ipslerr_p(3,'veget_config','VEGET_UPDATE can only be 0Y or 1Y.',&
            'Please correcte run.def file for VEGET_UPDATE','')
    END IF

    
    ! Coherence test for impveg and veget_update. Land use change can not be activated with impveg.
    IF (impveg .AND. veget_update > 0) THEN
       WRITE(numout,*) 'veget_update=',veget_update,' is not coeherent with impveg=',impveg
       CALL ipslerr_p(3,'slowproc_init','Incoherent values between impveg and veget_update', &
            'VEGET_UPDATE must be equal to 0Y if IMPOSE_VEG=y (impveg=true)','')
    END IF  
  

    !Config Key   = VEGETMAP_RESET
    !Config Desc  = Flag to change vegetation map without activating LAND USE change for carbon fluxes. At the same time carbon related variables are reset to zero.
    !Config If    = 
    !Config Def   = n
    !Config Help  = Use this option to change vegetation map while keeping VEGET_UPDATE=0Y
    !Config Units = [FLAG]
    vegetmap_reset = .FALSE.
    CALL getin_p('VEGETMAP_RESET', vegetmap_reset)


    !Config Key   = NINPUT_REINIT
    !Config Desc  = booleen to indicate that a new N INPUT file will be used.
    !Config If    = -
    !Config Def   = y
    !Config Help  = When set to y, the counter for the year of data grabbed in
    !Config         the Nitrogen input files will be reset to be equal to that of
    !Config         the first year present in the input file.  
    !Config         Then it is possible to change N INPUT file.
    !Config         Only seems to be
    !Config         useful set to n when the Nitrogen input file contains multiple years?
    !Config Units = [FLAG]
    ninput_reinit = .TRUE.
    CALL getin_p('NINPUT_REINIT', ninput_reinit)
    
    !Config Key   = NINPUT_YEAR
    !Config Desc  = Year of the N input map to be read
    !Config If    = -
    !Config Def   = 1
    !Config Help  = First year for N inputs vegetation
    !Config         If NINPUT_YEAR is set to 0, this means there is no time axis in the Nitrogen
    !Config         input map.  If there is a time axis, NINPUT_YEAR can be set to any four digit year.
    !Config         The code will look for an input file name corresponding to this year and
    !Config         take the data from this file.
    !Config Units = [FLAG]
    ninput_year_orig = 0
    CALL getin_p('NINPUT_YEAR', ninput_year_orig)
    
    !Config Key   = NINPUT_SUFFIX_YEAR
    !Config Desc  = Do the Ninput dataset have a 'year' suffix
    !Config If    = -
    !Config Def   = false
    !Config Help  = A flag to indicate if nitrogen input files have a year suffix (before .nc)
    !Config         If NINPUT_SUFFIX_YEAR is set to true, the code searches for a Nitrogen input file of the 
    !Config         format "filename_YEAR.nc" where YEAR is the NINPUT_YEAR
    !Config Units = [FLAG]
    ninput_suffix_year = .FALSE.
    CALL getin_p('NINPUT_SUFFIX_YEAR', ninput_suffix_year)

    !Config Key   = OK_PEAT_NODISCRETISATION
    !Config Desc  = Activate an ancien module for peatland carbon accumus 
    !Config If    = -
    !Config Def   = false
    !Config Help  = 
    !Config Units = [FLAG]
    ok_peat_NoDiscretisation = .FALSE.
    CALL getin_p('OK_PEAT_NODISCRETISATION', ok_peat_NoDiscretisation)
    
    !Config Key   = PERMA_PEAT
    !Config Desc  = Activate the carbon computation in peatland
    !Config If    = -
    !Config Def   = false
    !Config Help  = 
    !Config Units = [FLAG]
    perma_peat = .FALSE.
    CALL getin_p('PERMA_PEAT', perma_peat)
    
    IF ( perma_peat) THEN
       IF (SIZE(peat_bulk_density) .NE. ngrnd) THEN
          CALL ipslerr_p(3,"veget_config","ERROR : NGRND is 18 for peat_bulk_density par default.", &
                                                       "Need to change the peat_bulk_density.", " ")
       ENDIF
    ENDIF

    !Config Key   = FRAC1
    !Config Desc  = Fractional parameter useful for soil carbon computation on peatland
    !Config If    = PERMA_PEAT
    !Config Def   = 0.95
    !Config Help  = 
    !Config Units = []
    frac1=0.95
    CALL getin_p('FRAC1', frac1) 

    !Config Key   = FRAC2
    !Config Desc  = Fractional parameter useful for soil carbon computation on peatland
    !Config If    = PERMA_PEAT
    !Config Def   = 0.05 
    !Config Help  = 
    !Config Units = []
    frac2=0.05
    CALL getin_p('FRAC2', frac2)

    !Config Key   = AGRI_PEAT
    !Config Desc  = Activate the flag to allow having agrigculture on peatland
    !Config If    = -
    !Config Def   = false
    !Config Help  = 
    !Config Units = [FLAG]
    agri_peat = .FALSE.
    CALL getin_p('AGRI_PEAT',agri_peat)
    
    ! The following three flags are related to agri_peat
    ! They will be updated later with more details

    !Config Key   = AGRI_PEAT_PROP
    !Config Desc  = To activate different way of agri_peat
    !Config If    = AGRI_PEAT
    !Config Def   = false
    !Config Help  = 
    !Config Units = [FLAG]
    agri_peat_prop = .FALSE.
    CALL getin_p('AGRI_PEAT_PROP',agri_peat_prop)

    !Config Key   = AGRI_PEAT_MINCROP
    !Config Desc  = To activate different way of agri_peat
    !Config If    = AGRI_PEAT
    !Config Def   = false
    !Config Help  = 
    !Config Units = [FLAG]
    agri_peat_MINcrop = .FALSE.
    CALL getin_p('AGRI_PEAT_MINCROP',agri_peat_mincrop)

    !Config Key   = AGRI_PEAT_MAXCROP
    !Config Desc  = To activate different way of agri_peat
    !Config If    = AGRI_PEAT
    !Config Def   = false
    !Config Help  = 
    !Config Units = [FLAG]
    agri_peat_MAXcrop = .FALSE.
    CALL getin_p('AGRI_PEAT_MAXCROP',agri_peat_maxcrop)
        

  END SUBROUTINE veget_config


!! ================================================================================================================================
!! SUBROUTINE   : veget_config
!!
!>\BRIEF         This subroutine reads in the configuration file the imposed values of the parameters for all SECHIBA modules.  
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCE(S) :
!!
!! FLOWCHART    :
!! \n
!_ ================================================================================================================================

  SUBROUTINE config_sechiba_parameters

    IMPLICIT NONE

    !! 0. Variables and parameters declaration

    !! 0.4 Local variables
    REAL(r_std) :: nudge_tau_mc     !! Temporary variable read from run.def
    REAL(r_std) :: nudge_tau_snow   !! Temporary variable read from run.def
    INTEGER     :: ilevel           !! index

    !_ ================================================================================================================================

    ! Global : parameters used by many modules
 
    !Config Key   = NSNOW
    !Config Desc  = The number of snow levels: 3 or 12
    !Config If    = OK_SECHIBA
    !Config Def   = 3
    !Config Help  = 
    !Config Units = [-]
    nsnow=3
    CALL getin_p('NSNOW',nsnow)
    IF (nsnow /= 3 .AND. nsnow /=12 ) CALL ipslerr_p(3,"config_sechiba_parameters","ERROR : NSNOW can only be 3 or 12.", &
                                                       "Change in orchidee.def", " ")

    !Config Key   = NICE
    !Config Desc  = The number of ice levels
    !Config If    = OK_ICE_SHEET
    !Config Def   = 8
    !Config Help  = 
    !Config Units = [-]
    nice=8
    CALL getin_p('NICE',nice)
    IF ( nice /=8 ) CALL ipslerr_p(3,"config_sechiba_parameters","ERROR : NICE can only be 8 for the moment.", &
         "Change in orchidee.def", " ")

    IF (ok_ice_sheet) THEN
       !Config Key   = FRAC_MIN_ICE
       !Config Desc  = Limit of fraction of ice in the grid cell where the ice mask is considered as 1
       !Config If    = OK_ICE_SHEET AND .NOT. READ_ICE_SHEET_MASK
       !Config Def   = 0.70
       !Config Help  = This value is used as a threashold when calculation the ice_sheet_mask from frac_nobio
       !Config Units = [-]
       frac_min_ice=0.70
       CALL getin_p('FRAC_MIN_ICE',frac_min_ice)

       !Config Key   = READ_ICE_SHEET_MASK
       !Config Desc  = Read the ice sheet mask from file instead for using frac_nobio
       !Config If    = OK_ICE_SHEET
       !Config Def   = n
       !Config Help  = When this flag is set to true, FRAC_MIN_ICE is not used
       !Config Units = [-]
       read_ice_sheet_mask=.FALSE.
       CALL getin_p("READ_ICE_SHEET_MASK", read_ice_sheet_mask)

    END IF

    !Config Key   = MAXMASS_SNOW
    !Config Desc  = The maximum mass of a snow
    !Config If    = OK_SECHIBA
    !Config Def   = 3000.
    !Config Help  = 
    !Config Units = [kg/m^2]
    maxmass_snow = 3000.
    CALL getin_p('MAXMASS_SNOW',maxmass_snow)
    !
    !Config Key   = SNOWCRI
    !Config Desc  = Sets the amount above which only sublimation occures 
    !Config If    = OK_SECHIBA
    !Config Def   = 1.5
    !Config Help  = 
    !Config Units = [kg/m^2]
    snowcri = 1.5
    CALL getin_p('SNOWCRI',snowcri)
    !
    !! Initialization of sneige
    sneige = snowcri/mille
    !
    !Config Key   = MIN_WIND
    !Config Desc  = Minimum wind speed
    !Config If    = OK_SECHIBA
    !Config Def   = 0.1
    !Config Help  = 
    !Config Units = [m/s]
    min_wind = 0.1
    CALL getin_p('MIN_WIND',min_wind)
    !
    !Config Key   = MIN_QC
    !Config Desc  = Minimum qc
    !Config If    = OK_SECHIBA
    !Config Def   = 1.e-4
    !Config Help  = The minimum value for qc (qc=drag*wind) used in coupled(enerbil) and forced mode (enerbil and diffuco) 
    !Config Units = 
    min_qc = 1.e-4
    CALL getin_p('MIN_QC',min_qc)
    !
    !Config Key   = MAX_SNOW_AGE
    !Config Desc  = Maximum period of snow aging 
    !Config If    = OK_SECHIBA
    !Config Def   = 50.
    !Config Help  = 
    !Config Units = [days?]
    max_snow_age = 50._r_std
    CALL getin_p('MAX_SNOW_AGE',max_snow_age)
    !
    !Config Key   = SNOW_TRANS
    !Config Desc  = Transformation time constant for snow
    !Config If    = OK_SECHIBA
    !Config Def   = 0.2
    !Config Help  = optimized on 04/07/2016
    !Config Units = [m]
    snow_trans = 0.2_r_std
    CALL getin_p('SNOW_TRANS',snow_trans)
    !
    !Config Key   = SNOW_TRANS_NOBIO
    !Config Desc  = Transformation time constant for snow on nobio
    !Config If    = OK_SECHIBA
    !Config Def   = 1.
    !Config Help  = optimized on 18/02/2019
    !Config Units = [m]
    snow_trans_nobio = 1.0_r_std
    CALL getin_p('SNOW_TRANS_NOBIO',snow_trans_nobio)
    !
    !Config Key   = SNOWA_AGED_NOBIO_VIS
    !Config Desc  = Minimum snow albedo value for nobio type after aging (dirty old snow), visible albedo
    !Config If    = OK_SECHIBA
    !Config Def   = 0.58
    !Config Help  = optimized on 30/06/2022
    !Config Units = [-]
    snowa_aged_nobio_vis = 0.58_r_std 
    CALL getin_p('SNOWA_AGED_NOBIO_VIS',snowa_aged_nobio_vis)
    !
    !Config Key   = SNOWA_AGED_NOBIO_NIR
    !Config Desc  = Minimum snow albedo value for nobio type after aging (dirty old snow), near infrared albedo
    !Config If    = OK_SECHIBA
    !Config Def   = 0.35
    !Config Help  = optimized on 30/06/2022
    !Config Units = [-]
    snowa_aged_nobio_nir = 0.35_r_std
    CALL getin_p('SNOWA_AGED_NOBIO_NIR',snowa_aged_nobio_nir)
    !
    !Config Key   = SNOWA_DEC_NOBIO_VIS
    !Config Desc  = Decay rate of snow albedo value for nobio type, visible albedo
    !Config If    = OK_SECHIBA
    !Config Def   = 0.28
    !Config Help  = optimized on 30/06/2022
    !Config Units = [-]
    snowa_dec_nobio_vis = 0.28_r_std
    CALL getin_p('SNOWA_DEC_NOBIO_VIS',snowa_dec_nobio_vis)
    !
    !Config Key   = SNOWA_DEC_NOBIO_NIR
    !Config Desc  = Decay rate of snow albedo value for nobio type, near infrared albedo
    !Config If    = OK_SECHIBA
    !Config Def   = 0.45
    !Config Help  = optimized on 30/06/2022
    !Config Units = [-]
    snowa_dec_nobio_nir = 0.45_r_std 
    CALL getin_p('SNOWA_DEC_NOBIO_NIR',snowa_dec_nobio_nir)
    !
    !Config Key   = OMG1
    !Config Desc  = Tuning constant for snow ageing on nobio areas
    !Config If    = OK_SECHIBA
    !Config Def   = 7.
    !Config Help  = optimized on 14/10/2020
    !Config Units = [-]
    omg1 = 7.0 
    CALL getin_p('OMG1',omg1)
    !
    !Config Key   = OMG2
    !Config Desc  = Tuning constant for snow ageing on nobio areas
    !Config If    = OK_SECHIBA
    !Config Def   = 4.
    !Config Help  = optimized on 14/10/2020
    !Config Units = [-]
    omg2 = 4.0
    CALL getin_p('OMG2',omg2)
    
    !Config Key = XRHOSMIN
    !Config Desc = Fresh snow density
    !Config Def = 50.
    !Config Help = 19/02/2024
    !Config Units = kg.m-3
    xrhosmin = 50.
    CALL getin_p('XRHOSMIN',xrhosmin)
    
    !Config Key   = OK_NUDGE_MC
    !Config Desc  = Activate nudging of soil moisture
    !Config Def   = n
    !Config If    = 
    !Config Help  = 
    !Config Units = [FLAG]
    ok_nudge_mc = .FALSE.
    CALL getin_p('OK_NUDGE_MC', ok_nudge_mc)

    !Config Key   = NUDGE_TAU_MC
    !Config Desc  = Relaxation time for nudging of soil moisture expressed in fraction of the day
    !Config Def   = 1
    !Config If    = OK_NUDGE_MC
    !Config Help  = 
    !Config Units = [-]
    nudge_tau_mc = 1.0
    CALL getin_p('NUDGE_TAU_MC', nudge_tau_mc)
    IF (nudge_tau_mc < dt_sechiba/one_day) CALL ipslerr_p(3, 'hydrol_initialize', &
         'NUDGE_TAU_MC is smaller than the time step in sechiba which is not allowed.', &
         'Set NUDGE_TAU_MC higher or equal to dt_sechiba/one_day','')
    ! Calculate alpha to be used in hydrol
    alpha_nudge_mc = dt_sechiba/(one_day*nudge_tau_mc)
    IF (printlev>=2) WRITE(numout, *) 'ok_nudge_mc, nudge_tau_mc, alpha_nudge_mc =', &
         ok_nudge_mc, nudge_tau_mc, alpha_nudge_mc

    !Config Key   = OK_NUDGE_SNOW
    !Config Desc  = Activate nudging of snow variables
    !Config Def   = n
    !Config If    = 
    !Config Help  = 
    !Config Units = [FLAG]
    ok_nudge_snow = .FALSE.
    CALL getin_p('OK_NUDGE_SNOW', ok_nudge_snow)

    !Config Key   = NUDGE_TAU_SNOW
    !Config Desc  = Relaxation time for nudging of snow variables
    !Config Def   = 1
    !Config If    = OK_NUDGE_SNOW
    !Config Help  = 
    !Config Units = [-]
    nudge_tau_snow = 1.0
    CALL getin_p('NUDGE_TAU_SNOW', nudge_tau_snow)
    IF (nudge_tau_snow < dt_sechiba/one_day) CALL ipslerr_p(3, 'hydrol_initialize', &
         'NUDGE_TAU_SNOW is smaller than the time step in sechiba which is not allowed.', &
         'Set NUDGE_TAU_SNOW higher or equal to dt_sechiba/one_day','')
    ! Calculate alpha to be used in hydrol
    alpha_nudge_snow = dt_sechiba/(one_day*nudge_tau_snow)
    IF (printlev>=2) WRITE(numout, *) 'ok_nudge_snow, nudge_tau_snow, alpha_nudge_snow =', &
         ok_nudge_snow, nudge_tau_snow, alpha_nudge_snow


    !Config Key   = NUDGE_INTERPOL_WITH_XIOS
    !Config Desc  = Activate reading and interpolation with XIOS for nudging fields
    !Config Def   = n
    !Config If    = OK_NUDGE_MC or OK_NUDGE_SNOW
    !Config Help  = 
    !Config Units = [FLAG]
    nudge_interpol_with_xios = .FALSE.
    CALL getin_p('NUDGE_INTERPOL_WITH_XIOS', nudge_interpol_with_xios)

    !-
    ! condveg
    !-
    !
    !Config Key   = HEIGHT_DISPLACEMENT
    !Config Desc  = Magic number which relates the height to the displacement height.
    !Config If    = OK_SECHIBA 
    !Config Def   = 0.75
    !Config Help  = 
    !Config Units = [m]
    height_displacement = 0.66
    CALL getin_p('HEIGHT_DISPLACEMENT',height_displacement)
    !
    !Config Key   = Z0_BARE
    !Config Desc  = bare soil roughness length
    !Config If    = OK_SECHIBA 
    !Config Def   = 0.01 
    !Config Help  = 
    !Config Units = [m]
    z0_bare = 0.01
    CALL getin_p('Z0_BARE',z0_bare)
    !
    !Config Key   = Z0_ICE
    !Config Desc  = ice roughness length
    !Config If    = OK_SECHIBA 
    !Config Def   = 0.001
    !Config Help  = 
    !Config Units = [m]
    z0_ice = 0.001 
    CALL getin_p('Z0_ICE',z0_ice)

    !Config Key   = TCST_SNOWA
    !Config Desc  = Time constant of the albedo decay of snow
    !Config If    = OK_SECHIBA 
    !Config Def   = 15.0 
    !Config Help  = 
    !Config Units = [days]
    tcst_snowa = 15.0
    CALL getin_p('TCST_SNOWA',tcst_snowa)
    !
    !Config Key   = TCST_SNOWA_NOBIO
    !Config Desc  = Time constant of the albedo decay of snow for nobio
    !Config If    = OK_SECHIBA 
    !Config Def   = 2.0 
    !Config Help  = optimized on 18/02/2019
    !Config Units = [days]
    tcst_snowa_nobio = 2.0   
    CALL getin_p('TCST_SNOWA_NOBIO',tcst_snowa_nobio)
    !
    !Config Key   = SNOWCRI_ALB
    !Config Desc  = Critical value for computation of snow albedo
    !Config If    = OK_SECHIBA
    !Config Def   = 2. 
    !Config Help  = 
    !Config Units = [cm]
    snowcri_alb = 10.
    CALL getin_p('SNOWCRI_ALB',snowcri_alb)
    !
    !
    !Config Key   = VIS_DRY
    !Config Desc  = The correspondance table for the soil color numbers and their albedo 
    !Config If    = OK_SECHIBA 
    !Config Def   = 0.24, 0.22, 0.20, 0.18, 0.16, 0.14, 0.12, 0.10, 0.27
    !Config Help  = 
    !Config Units = [-]
    vis_dry = (/0.24, 0.22, 0.20, 0.18, 0.16, 0.14, 0.12, 0.10, 0.27/)
    CALL getin_p('VIS_DRY',vis_dry)
    !
    !Config Key   = NIR_DRY
    !Config Desc  = The correspondance table for the soil color numbers and their albedo 
    !Config If    = OK_SECHIBA 
    !Config Def   = 0.48, 0.44, 0.40, 0.36, 0.32, 0.28, 0.24, 0.20, 0.55
    !Config Help  = 
    !Config Units = [-]
    nir_dry = (/0.48, 0.44, 0.40, 0.36, 0.32, 0.28, 0.24, 0.20, 0.55/)
    CALL getin_p('NIR_DRY',nir_dry)
    !
    !Config Key   = VIS_WET 
    !Config Desc  = The correspondance table for the soil color numbers and their albedo
    !Config If    = OK_SECHIBA  
    !Config Def   = 0.12, 0.11, 0.10, 0.09, 0.08, 0.07, 0.06, 0.05, 0.15
    !Config Help  = 
    !Config Units = [-]
    vis_wet = (/0.12, 0.11, 0.10, 0.09, 0.08, 0.07, 0.06, 0.05, 0.15/)
    CALL getin_p('VIS_WET',vis_wet)
    !
    !Config Key   = NIR_WET
    !Config Desc  = The correspondance table for the soil color numbers and their albedo 
    !Config If    = OK_SECHIBA 
    !Config Def   = 0.24, 0.22, 0.20, 0.18, 0.16, 0.14, 0.12, 0.10, 0.31
    !Config Help  = 
    !Config Units = [-]
    nir_wet = (/0.24, 0.22, 0.20, 0.18, 0.16, 0.14, 0.12, 0.10, 0.31/)
    CALL getin_p('NIR_WET',nir_wet)
    !
    !Config Key   = ALBSOIL_VIS
    !Config Desc  = 
    !Config If    = OK_SECHIBA 
    !Config Def   = 0.18, 0.16, 0.16, 0.15, 0.12, 0.105, 0.09, 0.075, 0.25
    !Config Help  = 
    !Config Units = [-]
    albsoil_vis = (/0.18, 0.16, 0.16, 0.15, 0.12, 0.105, 0.09, 0.075, 0.25/)
    CALL getin_p('ALBSOIL_VIS',albsoil_vis)
    !
    !Config Key   = ALBSOIL_NIR 
    !Config Desc  = 
    !Config If    = OK_SECHIBA 
    !Config Def   = 0.36, 0.34, 0.34, 0.33, 0.30, 0.25, 0.20, 0.15, 0.45
    !Config Help  = 
    !Config Units = [-]
    albsoil_nir = (/0.36, 0.34, 0.34, 0.33, 0.30, 0.25, 0.20, 0.15, 0.45/)
    CALL getin_p('ALBSOIL_NIR',albsoil_nir)
    !-
    !Config Key   = ALB_DEADLEAF 
    !Config Desc  = albedo of dead leaves, VIS+NIR 
    !Config If    = OK_SECHIBA 
    !Config Def   = 0.12, 0.35
    !Config Help  = 
    !Config Units = [-]
    alb_deadleaf = (/ .12, .35/)
    CALL getin_p('ALB_DEADLEAF',alb_deadleaf)
    !
    !Config Key   = ALB_ICE
    !Config Desc  = albedo of ice, VIS+NIR
    !Config If    = OK_SECHIBA
    !Config Def   = 0.60, 0.20
    !Config Help  = 
    !Config Units = [-]
    alb_ice = (/ .60, .20/)
    CALL getin_p('ALB_ICE',alb_ice)
    !
    !Config Key   = ZSICOEF1
    !Config Desc  = Thickness of ice levels
    !Config If    = OK_ICE_SHEET
    !Config Def   = 0.01, 0.05, 0.15, 0.5, 1., 5., 10., 50.
    !Config Help  = 
    !Config Units = [m]
    ZSICOEF1  = (/0.01, 0.05, 0.15, 0.5, 1., 5., 10., 50./)
    CALL getin_p('ZSICOEF1',zsicoef1)
    
    !
    ! Get the fixed snow albedo if needed
    !
    !Config Key   = CONDVEG_SNOWA
    !Config Desc  = The snow albedo used by SECHIBA
    !Config Def   = 1.E+20
    !Config if    = OK_SECHIBA
    !Config Help  = This option allows the user to impose a snow albedo.
    !Config         Default behaviour is to use the model of snow albedo
    !Config         developed by Chalita (1993).
    !Config Units = [-]
    fixed_snow_albedo = undef_sechiba
    CALL getin_p('CONDVEG_SNOWA',fixed_snow_albedo)
    !
    !Config Key   = ALB_BARE_MODEL
    !Config Desc  = Switch bare soil albedo dependent (if TRUE) on soil wetness
    !Config Def   = n
    !Config if    = OK_SECHIBA
    !Config Help  = If TRUE, the model for bare soil albedo is the old formulation.
    !Config         Then it depend on the soil dry or wetness. If FALSE, it is the 
    !Config         new computation that is taken, it is the mean of soil albedo.
    !Config Units = [FLAG]
    alb_bare_model = .FALSE.
    CALL getin_p('ALB_BARE_MODEL',alb_bare_model)
    !
    !Config Key   = ALB_BG_MODIS
    !Config Desc  = Read bare soil albedo from file with background MODIS data only if not found in restart file
    !Config Def   = n
    !Config if    = OK_SECHIBA
    !Config Help  = If TRUE, the bare soil albedo is read from file
    !Config         based on background MODIS data.  
    !Config         If FALSE, computaion depends on ALB_BARE_MODEL
    !Config Units = [FLAG]
    alb_bg_modis = .FALSE.
    CALL getin_p('ALB_BG_MODIS',alb_bg_modis)

    IF (alb_bg_modis) THEN
       !Config Key   = ALB_BG_MODIS_REINIT
       !Config Desc  = Impose reading of bare soil albedo from file even if it is in the restart file
       !Config Def   = n
       !Config if    = ALB_BG_MODIS
       !Config Help  = If TRUE, the bare soil albedo is read from file
       !Config         based on background MODIS data only if ALB_BG_MOIDS=y
       !Config         If FALSE, file will still be read if not found in the restart file
       !Config Units = [FLAG]
       alb_bg_modis_reinit = .FALSE.
       CALL getin_p('ALB_BG_MODIS_REINIT',alb_bg_modis_reinit)

       IF (alb_bg_modis_reinit) THEN
          IF (printlev>=1) WRITE(numout,*)'Read bare soil albedo from file when starting the model. '
          IF (printlev>=1) WRITE(numout,*)'The albedo file will always be read even if the variables is found in the restart file'
       ELSE
          IF (printlev>=1) WRITE(numout,*)'Read bare soil albedo from file when starting the model only if no restart file is found or the variable is missing. '
       END IF

    ELSE
       IF (printlev>=1) WRITE(numout,*)'Reading of background albedo is not activated.'
       IF (printlev>=1) WRITE(numout,*)'ALB_BG_MODIS=',alb_bg_modis, ' ALB_BARE_MODEL=', alb_bare_model
    END IF
    !
    !Config Key   = IMPOSE_AZE
    !Config Desc  = Should the surface parameters be prescribed
    !Config Def   = n
    !Config if    = OK_SECHIBA
    !Config Help  = This flag allows the user to impose the surface parameters
    !Config         (Albedo Roughness and Emissivity). It is espacially interesting for 0D
    !Config         simulations. On the globe it does not make too much sense as
    !Config         it imposes the same vegetation everywhere
    !Config Units = [FLAG]
    impaze = .FALSE.
    CALL getin_p('IMPOSE_AZE',impaze)
    !
    IF(impaze) THEN
       !
       !Config Key   = CONDVEG_Z0
       !Config Desc  = Surface roughness
       !Config Def   = 0.15
       !Config If    = IMPOSE_AZE
       !Config Help  = Surface rougness to be used on the point if a 0-dim version
       !Config         of SECHIBA is used. Look at the description of the forcing  
       !Config         data for the correct value.
       !Config Units = [m]
       z0_scal = 0.15 
       CALL getin_p('CONDVEG_Z0', z0_scal)  
       !
       !Config Key   = ROUGHHEIGHT
       !Config Desc  = Height to be added to the height of the first level
       !Config Def   = 0.0
       !Config If    = IMPOSE_AZE
       !Config Help  = ORCHIDEE assumes that the atmospheric level height is counted
       !Config         from the zero wind level. Thus to take into account the roughness
       !Config         of tall vegetation we need to correct this by a certain fraction
       !Config         of the vegetation height. This is called the roughness height in
       !Config         ORCHIDEE talk.
       !Config Units = [m]
       roughheight_scal = zero 
       CALL getin_p('ROUGHHEIGHT', roughheight_scal)
       ! 
       !Config Key   = CONDVEG_ALBVIS
       !Config Desc  = SW visible albedo for the surface
       !Config Def   = 0.25
       !Config If    = IMPOSE_AZE
       !Config Help  = Surface albedo in visible wavelengths to be used 
       !Config         on the point if a 0-dim version of SECHIBA is used. 
       !Config         Look at the description of the forcing data for 
       !Config         the correct value.
       !Config Units = [-]
       albedo_scal = (/ 0.25, 0.25 /) 
       CALL getin_p('CONDVEG_ALBVIS', albedo_scal(ivis))
       !
       !Config Key   = CONDVEG_ALBNIR
       !Config Desc  = SW near infrared albedo for the surface
       !Config Def   = 0.25
       !Config If    = IMPOSE_AZE
       !Config Help  = Surface albedo in near infrared wavelengths to be used 
       !Config         on the point if a 0-dim version of SECHIBA is used. 
       !Config         Look at the description of the forcing data for 
       !Config         the correct value.
       !Config Units = [-]
       CALL getin_p('CONDVEG_ALBNIR', albedo_scal(inir))
       !
       !Config Key   = CONDVEG_EMIS
       !Config Desc  = Emissivity of the surface for LW radiation
       !Config Def   = 1.0
       !Config If    = IMPOSE_AZE
       !Config Help  = The surface emissivity used for compution the LE emission
       !Config         of the surface in a 0-dim version. Values range between 
       !Config         0.97 and 1.. The GCM uses 0.98.
       !Config Units = [-]
       emis_scal = 1.0
       CALL getin_p('CONDVEG_EMIS', emis_scal)
    ENDIF

    new_watstress = .FALSE.
    CALL getin_p('NEW_WATSTRESS',new_watstress)
    IF(new_watstress) THEN
       alpha_watstress = 1.
       CALL getin_p('ALPHA_WATSTRESS',alpha_watstress)
    ENDIF

    !
    !Config Key   = ROUGH_DYN
    !Config Desc  = Account for a dynamic roughness height
    !Config Def   = y
    !Config if    = OK_SECHIBA
    !Config Help  = If this flag is set to true (y) then the roughness
    !Config         height is computed dynamically, varying with LAI
    !Config Units = [FLAG]
    rough_dyn = .TRUE. 
    CALL getin_p('ROUGH_DYN',rough_dyn)

    IF ( rough_dyn ) THEN
       
       !Config Key   = USE_RATIO_Z0M_Z0H
       !Config Desc  = To impose a constant ratio as in ROUGH_DYN=F
       !Config Def   = TRUE
       !Config If    = ROUGH_DYN
       !Config Help  = 
       !Config Units = [-]
       use_ratio_z0m_z0h = .TRUE.
       CALL getin_p('USE_RATIO_Z0M_Z0H', use_ratio_z0m_z0h)
       !
       !Config Key   = C1
       !Config Desc  = Constant used in the formulation of the ratio of 
       !Config         the ratio of friction velocity to the wind speed
       !Config         at the canopy top 
       !Config         See Ershadi et al. (2015) for more info
       !Config Def   = 0.32
       !Config If    = ROUGH_DYN
       !Config Help  = 
       !Config Units = [-]
       c1 = 0.32
       CALL getin_p('C1', c1)
       !
       !Config Key   = C2
       !Config Desc  = Constant used in the formulation of the ratio of 
       !Config         the ratio of friction velocity to the wind speed
       !Config         at the canopy top 
       !Config         See Ershadi et al. (2015) for more info
       !Config Def   = 0.264
       !Config If    = ROUGH_DYN
       !Config Help  = 
       !Config Units = [-]
       c2 = 0.264
       CALL getin_p('C2', c2)
       !
       !Config Key   = C3
       !Config Desc  = Constant used in the formulation of the ratio of 
       !Config         the ratio of friction velocity to the wind speed
       !Config         at the canopy top 
       !Config         See Ershadi et al. (2015) for more info
       !Config Def   = 15.1
       !Config If    = ROUGH_DYN
       !Config Help  = 
       !Config Units = [-]
       c3 = 15.1 
       CALL getin_p('C3', c3)
       !
       !Config Key   = Cdrag_foliage
       !Config Desc  = Drag coefficient of the foliage
       !Config         See Ershadi et al. (2015) and Su et al. (2001) 
       !Config         for more info
       !Config Def   = 0.2
       !Config If    = ROUGH_DYN
       !Config Help  = 
       !Config Units = [-]
       Cdrag_foliage = 0.2  
       CALL getin_p('CDRAG_FOLIAGE', Cdrag_foliage)
       !
       !Config Key   = Ct
       !Config Desc  = Heat transfer coefficient of the leaf
       !Config         See Ershadi et al. (2015) and Su et al. (2001) 
       !Config         for more info
       !Config Def   = 0.01
       !Config If    = ROUGH_DYN
       !Config Help  = 
       !Config Units = [-]
       Ct = 0.01
       CALL getin_p('CT', Ct)
       !
       !Config Key   = Prandtl
       !Config Desc  = Prandtl number used in the calculation of Ct*
       !Config         See Su et al. (2001) for more info
       !Config Def   = 0.71
       !Config If    = ROUGH_DYN
       !Config Help  = 
       !Config Units = [-]
       Prandtl = 0.71 
       CALL getin_p('PRANDTL', Prandtl)
    ENDIF
    !- 
    ! Variables related to the explicitsnow module
    !-
    !Config Key = xansmax 
    !Config Desc = maximum snow albedo
    !Config If = OK_SECHIBA
    !Config Def = 0.85
    !Config Help = 
    !Config Units = [-]
    xansmax = 0.85
    CALL getin_p('XANSMAX',xansmax)
    !
    !Config Key = xansmin 
    !Config Desc = minimum snow albedo
    !Config If = OK_SECHIBA
    !Config Def = 0.50
    !Config Help = 
    !Config Units = [-]
    xansmin = 0.50
    CALL getin_p('XANSMIN',xansmin)
    !
    !Config Key = xans_todry 
    !Config Desc = albedo decay rate for the dry snow 
    !Config If = OK_SECHIBA
    !Config Def = 0.008
    !Config Help = 
    !Config Units = [S-1]
    xans_todry = 0.008
    CALL getin_p('XANSDRY',xans_todry)
    !
    !Config Key = xans_t 
    !Config Desc = albedo decay rate for the wet snow 
    !Config If = OK_SECHIBA
    !Config Def = 0.24
    !Config Help = 
    !Config Units = [S-1]
    xans_t = 0.240
    CALL getin_p('XANS_T',xans_t)
    !
    !Config Key = xrhosmax 
    !Config Desc = maximum snow density 
    !Config If = OK_SECHIBA
    !Config Def = 750
    !Config Help = 
    !Config Units = [kg m-3]
    xrhosmax = 750.
    CALL getin_p('XRHOSMAX',xrhosmax)
    !
    !Config Key = xrhosmin
    !Config Desc = minimum snow density 
    !Config If = OK_SECHIBA
    !Config Def = 50
    !Config Help = 
    !Config Units = [kg m-3]
    xrhosmin = 50.
    CALL getin_p('XRHOSMIN',xrhosmin)
    !
    !Config Key = xwsnowholdmax1
    !Config Desc = snow holding capacity 1
    !Config If = OK_SECHIBA
    !Config Def = 0.03
    !Config Help = 
    !Config Units = [-]
    xwsnowholdmax1 = 0.03
    CALL getin_p('XWSNOWHOLDMAX1',xwsnowholdmax1)
    !
    !Config Key = xwsnowholdmax2
    !Config Desc = snow holding capacity 2
    !Config If = OK_SECHIBA
    !Config Def = 0.10
    !Config Help = 
    !Config Units = [-]
    xwsnowholdmax2 = 0.10
    CALL getin_p('XWSNOWHOLDMAX2',xwsnowholdmax2)
    !
    !Config Key = xsnowrhohold 
    !Config Desc = snow density 
    !Config If = OK_SECHIBA
    !Config Def = 200.0
    !Config Help = 
    !Config Units = [kg/m3]
    xsnowrhohold = 200.0
    CALL getin_p('XSNOWRHOHOLD',xsnowrhohold)
    ! 
    !Config Key = ZSNOWTHRMCOND1
    !Config Desc = Thermal conductivity Coef 1
    !Config If = OK_SECHIBA
    !Config Def = 0.02 
    !Config Help = 
    !Config Units = [W/m/K]
    ZSNOWTHRMCOND1 = 0.02
    CALL getin_p('ZSNOWTHRMCOND1',ZSNOWTHRMCOND1)
    !
    !Config Key = ZSNOWTHRMCOND2
    !Config Desc = Thermal conductivity Coef 2
    !Config If = OK_SECHIBA
    !Config Def = 2.5E-6
    !Config Help = 
    !Config Units = [W m5/(kg2 K)]
    ZSNOWTHRMCOND2 = 2.5e-6
    CALL getin_p('ZSNOWTHRMCOND2',ZSNOWTHRMCOND2)
    !
    !Config Key = ZSNOWTHRMCOND_AVAP
    !Config Desc = Thermal conductivity Coef 1 water vapor
    !Config If = OK_SECHIBA
    !Config Def = -0.06023
    !Config Help = 
    !Config Units = [W/m/K]
    ZSNOWTHRMCOND_AVAP = -0.06023
    CALL getin_p('ZSNOWTHRMCOND_AVAP',ZSNOWTHRMCOND_AVAP)
    !
    !Config Key = ZSNOWTHRMCOND_BVAP
    !Config Desc = Thermal conductivity Coef 2 water vapor
    !Config If = OK_SECHIBA
    !Config Def = -2.5425
    !Config Help = 
    !Config Units = [W/m]
    ZSNOWTHRMCOND_BVAP = -2.5425
    CALL getin_p('ZSNOWTHRMCOND_BVAP',ZSNOWTHRMCOND_BVAP)
    !
    !Config Key = ZSNOWTHRMCOND_CVAP
    !Config Desc = Thermal conductivity Coef 3 water vapor
    !Config If = OK_SECHIBA
    !Config Def = -289.99
    !Config Help = 
    !Config Units = [K]
    ZSNOWTHRMCOND_CVAP = -289.99
    CALL getin_p('ZSNOWTHRMCOND_CVAP',ZSNOWTHRMCOND_CVAP)

    !Snow compaction factors
    !Config Key = ZSNOWCMPCT_RHOD
    !Config Desc = Snow compaction coefficent
    !Config If = OK_SECHIBA
    !Config Def = 150.0
    !Config Help = 
    !Config Units = [kg/m3]
    ZSNOWCMPCT_RHOD = 150.0 
    CALL getin_p('ZSNOWCMPCT_RHOD',ZSNOWCMPCT_RHOD)

    !Config Key = ZSNOWCMPCT_ACM
    !Config Desc = Coefficent for the thermal conductivity
    !Config If = OK_SECHIBA
    !Config Def = 2.8e-6
    !Config Help = 
    !Config Units = [1/s]
    ZSNOWCMPCT_ACM = 2.8e-6
    CALL getin_p('ZSNOWCMPCT_ACM',ZSNOWCMPCT_ACM)

    !Config Key = ZSNOWCMPCT_BCM
    !Config Desc = Coefficent for the thermal conductivity
    !Config If = OK_SECHIBA
    !Config Def = 0.04
    !Config Help = 
    !Config Units = [1/K]
    ZSNOWCMPCT_BCM = 0.04
    CALL getin_p('ZSNOWCMPCT_BCM',ZSNOWCMPCT_BCM)

    !Config Key = ZSNOWCMPCT_CCM
    !Config Desc = Coefficent for the thermal conductivity
    !Config If = OK_SECHIBA
    !Config Def = 460.
    !Config Help = 
    !Config Units = [m3/kg]
    ZSNOWCMPCT_CCM = 460.
    CALL getin_p('ZSNOWCMPCT_CCM',ZSNOWCMPCT_CCM)

    !Config Key = eta_0
    !Config Desc = eta_0
    !Config If = OK_SECHIBA
    !Config Def = 7.62237e6
    !Config Help =
    !Config Units = [Pa/s]
    eta_0 = 7.62237e6 
    CALL getin_p('eta_0',eta_0)

    !Config Key = c_eta
    !Config Desc = c_eta vionnet viscosity
    !Config If = OK_SECHIBA
    !Config Def = 358.
    !Config Help =
    !Config Units = [kg.m-3]
    c_eta = 358.0
    CALL getin_p('c_eta',c_eta)

    !Config Key = a_eta
    !Config Desc = a_eta Anderson/Vionnet viscosity
    !Config If = OK_SECHIBA
    !Config Def = 0.1
    !Config Help =
    !Config Units =
    a_eta = 0.1
    CALL getin_p('a_eta',a_eta)

    !Config Key = b_eta
    !Config Desc = b_eta Anderson/Vionnet viscosity
    !Config If = OK_SECHIBA
    !Config Def = 0.023
    !Config Help =
    !Config Units =
    b_eta = 0.023
    CALL getin_p('b_eta',b_eta)
    
    !Config Key = gs_vionnet
    !Config Desc = gs_vionnet
    !Config If = OK_SECHIBA
    !Config Def = 0.1
    !Config Help =
    !Config Units =
    gs_vionnet = 0.1
    CALL getin_p('gs_vionnet',gs_vionnet)

    !Config Key = RHOS_WS
    !Config Desc = FLAG rhos wind formula
    !Config If = OK_SECHIBA
    !Config Def = .FALSE.
    !Config Help =
    !Config Units = [-]
    RHOS_WS = .FALSE.
    CALL getin_p('RHOS_WS',RHOS_WS)

    !Config Key = INIT_SNOWPACK
    !Config Desc = FLAG init snowpack
    !Config If = OK_SECHIBA
    !Config Def = .FALSE.
    !Config Help =
    !Config Units = [-]
    INIT_SNOWPACK = .FALSE.
    CALL getin_p('INIT_SNOWPACK',INIT_SNOWPACK)

    !Config Key = SURFDENSITY_XRHOSMIN
    !Config Desc = FLAG init surfdensity
    !Config If = OK_SECHIBA
    !Config Def = .FALSE.
    !Config Help =
    !Config Units = [-]
    SURFDENSITY_XRHOSMIN = .FALSE.
    CALL getin_p('SURFDENSITY_XRHOSMIN',SURFDENSITY_XRHOSMIN)

    !Config Key = ZSNOWCMPCT_V0
    !Config Desc = Vapor coefficent for the thermal conductivity
    !Config If = OK_SECHIBA
    !Config Def = 3.7e7
    !Config Help = 
    !Config Units = [Pa/s]
    ZSNOWCMPCT_V0 = 3.7e7 
    CALL getin_p('ZSNOWCMPCT_V0',ZSNOWCMPCT_V0)

    !Config Key = ZSNOWCMPCT_VT
    !Config Desc = Vapor coefficent for the thermal conductivity
    !Config If = OK_SECHIBA
    !Config Def = 0.081
    !Config Help = 
    !Config Units = [1/K]
    ZSNOWCMPCT_VT = 0.081 
    CALL getin_p('ZSNOWCMPCT_VT',ZSNOWCMPCT_VT)

    !Config Key = ZSNOWCMPCT_VR
    !Config Desc = Vapor coefficent for the thermal conductivity
    !Config If = OK_SECHIBA
    !Config Def = 0.018
    !Config Help = 
    !Config Units = [m3/kg]
    ZSNOWCMPCT_VR = 0.018
    CALL getin_p('ZSNOWCMPCT_VR',ZSNOWCMPCT_VR)

    !Config Key = ZMAX
    !Config Desc = Snow grid parameter : max thickness of each layer
    !Config If = OK_SECHIBA
    !Config Def = (/0.01, 0.05, 0.15, 0.5, 1., 0., 0., 0., 1., 0.5, 0.1, 0.02/)
    !Config Help =
    !Config Units = [m]
    Zmax = (/0.01, 0.05, 0.15, 0.5, 1., 0., 0., 0., 1., 0.5, 0.1, 0.02/)
    CALL getin_p('ZMAX',Zmax)
    
    !Config Key = xrhosref
    !Config Desc = Reference density used in the viscosity calculation
    !Config If = OK_SECHIBA
    !Config Def = 250
    !Config Help = 
    !Config Units = [kg/m3]
    xrhosref = 250.
    CALL getin_p('XRHOSREF',xrhosref)

    !Config Key = COMPACT_SNOW_METHOD
    !Config Desc = Flag to choose the way to calculate snow compaction
    !Config If = OK_SECHIBA
    !Config Def = 1
    !Config Help = 1:original 2:follows compaction Decharme et al 2016
    !Config Units = [-]
    COMPACT_SNOW_METHOD = 1
    CALL getin_p('COMPACT_SNOW_METHOD',COMPACT_SNOW_METHOD)

    !Config Key = VISCOSITY_METHOD
    !Config Desc = Flag to choose the way to calculate snow viscosity
    !Config If = OK_SECHIBA
    !Config Def = 1
    !Config Help = 1:Anderson 2:Vionnet 2012 (CROCUS)
    !Config Units = [-]
    VISCOSITY_METHOD = 1
    CALL getin_p('VISCOSITY_METHOD',VISCOSITY_METHOD)

    !Config Key = xrhosmob
    !Config Desc = Reference density for the calculation of the mobility index
    !Config If = OK_SECHIBA
    !Config Def = 295
    !Config Help = 
    !Config Units = [kg/m3]
    xrhosmob = 295.
    CALL getin_p('XRHOSMOB',xrhosmob)

    !Config Key = xrhoswind
    !Config Desc = Density threshold for wind-driven compaction 
    !Config If = OK_SECHIBA
    !Config Def = 350
    !Config Help = 
    !Config Units = [kg/m3]
    xrhoswind = 350.
    CALL getin_p('XRHOSWIND',xrhoswind)

    !Config Key = ZSNOWCMPCT_X0
    !Config Desc = Viscosity coefficient 
    !Config If = OK_SECHIBA
    !Config Def = 7.82237e6
    !Config Help = 
    !Config Units = [Pa/s]
    ZSNOWCMPCT_X0 = 7.62237e6
    CALL getin_p('ZSNOWCMPCT_X0',ZSNOWCMPCT_X0)

    !Config Key = ZSNOWCMPCT_XT
    !Config Desc = Viscosity coefficient 
    !Config If = OK_SECHIBA
    !Config Def = 0.1
    !Config Help = 
    !Config Units = [1/K]
    ZSNOWCMPCT_XT = 0.1
    CALL getin_p('ZSNOWCMPCT_XT',ZSNOWCMPCT_XT)

    !Config Key = ZSNOWCMPCT_XR
    !Config Desc = Viscosity coefficient 
    !Config If = OK_SECHIBA
    !Config Def = 0.023
    !Config Help = 
    !Config Units = [m3/kg]
    ZSNOWCMPCT_XR = 0.023
    CALL getin_p('ZSNOWCMPCT_XR',ZSNOWCMPCT_XR)

    !Config Key = DTref
    !Config Desc = temperature reference used for the viscosity calculation 
    !Config If = OK_SECHIBA
    !Config Def = 5.0
    !Config Help = 
    !Config Units = [K]
    DTref = 5.
    CALL getin_p('DTREF',DTref)

    !Config Key = ACMPCT
    !Config Desc = Coefficient used for the calculatio of the wind-driven compaction index 
    !Config If = OK_SECHIBA
    !Config Def = 2.868
    !Config Help = 
    !Config Units = [-]
    ACMPCT = 2.868
    CALL getin_p('ACMPCT',ACMPCT)

    !Config Key = BCMPCT
    !Config Desc = Coefficient used for the calculatio of the wind-driven compaction index 
    !Config If = OK_SECHIBA
    !Config Def = 0.085
    !Config Help = 
    !Config Units = [s/m]
    BCMPCT = 0.085
    CALL getin_p('BCMPCT',BCMPCT)

    !Config Key = KCMPCT
    !Config Desc = Coefficient used for the calculatio of the wind-driven compaction index 
    !Config If = OK_SECHIBA
    !Config Def = 1.25
    !Config Help = 
    !Config Units = [-]
    KCMPCT = 1.25 
    CALL getin_p('KCMPCT',KCMPCT)

    !Config Key = a_tau
    !Config Desc = Coefficient used for the calculatio of the compaction rate
    !Config If = OK_SECHIBA
    !Config Def = 10
    !Config Help = 
    !Config Units = [-]
    a_tau = 10 
    CALL getin_p('A_TAU',a_tau)

    !Config Key = b_tau
    !Config Desc = Coefficient used for the calculatio of the compaction rate
    !Config If = OK_SECHIBA
    !Config Def = 3.25
    !Config Help = 
    !Config Units = [-]
    b_tau = 3.25
    CALL getin_p('B_TAU',b_tau)

    !Config Key = Amob
    !Config Desc = Coefficient used for the calculatio of the mobility index
    !Config If = OK_SECHIBA
    !Config Def = 1.25
    !Config Help = 
    !Config Units = [-]
    Amob = 1.25 
    CALL getin_p('AMOB',Amob)
    
    !Config Key = ZICETHRMCOND1
    !Config Desc = Ice Thermal conductivity Coef 1
    !Config If = OK_SECHIBA
    !Config Def = 6.727
    !Config Help = 
    !Config Units = [W/m/K]
    ZICETHRMCOND1 = 6.727 
    CALL getin_p('ZICETHRMCOND1',ZICETHRMCOND1)
    
    !Config Key = ZICETHRMCOND2
    !Config Desc = Ice Thermal conductivity Coef 2
    !Config If = OK_SECHIBA
    !Config Def = -0.0041
    !Config Help = 
    !Config Units = [W/m/K]
    ZICETHRMCOND2 = -0.0041
    CALL getin_p('ZICETHRMCOND2',ZICETHRMCOND2)
    
    !Config Key = ZICETHRMHEAT1
    !Config Desc = Ice Thermal heat capacity Coef 1
    !Config If = OK_SECHIBA
    !Config Def = 2115.3
    !Config Help = 
    !Config Units = [W/m/K]
    ZICETHRMHEAT1 = 2115.3
    CALL getin_p('ZICETHRMHEAT1',ZICETHRMHEAT1)
    
    !Config Key = ZICETHRMHEAT2
    !Config Desc = Ice Thermal heat capacity Coef 2
    !Config If = OK_SECHIBA
    !Config Def = 7.79293
    !Config Help = 
    !Config Units = [W/m/K]
    ZICETHRMHEAT2 = 7.79293
    CALL getin_p('ZICETHRMHEAT2',ZICETHRMHEAT2)


    !
    ! Variables related to crop irrigation
    !  

    !Config Key   = IRRIG_DOSMAX
    !Config Desc  = The maximum irrigation water injected per hour (kg.m^{-2}/hour)
    !Config If    = DO_IRRIGATION AND old_irrig_scheme = FALSE
    !Config Def   = 3.
    !Config Help  =
    !Config Units = [kg.m^{-2}/hour]
    irrig_dosmax = 3. 
    CALL getin_p('IRRIG_DOSMAX',irrig_dosmax)

    !Config Key   = CUM_DH_THR
    !Config Desc  = Threshold depth to define root zone, and calculate water deficit for irrigation
    !Config If    = DO_IRRIGATION AND OLD_IRRIG_SCHEME = FALSE
    !Config Def   = 0.64
    !Config Help  =
    !Config Units = [m]
    cum_dh_thr = 0.64
    CALL getin_p('CUM_DH_THR',cum_dh_thr)

    !Crop irrigation
    !Config Key   = IRRIGATED_SOILTILE
    !Config Desc  = Do we introduce a new soil tile for irrigated croplands?
    !Config If    = DO_IRRIGATION AND old_irrig_scheme = FALSE
    !Config Def   = n
    !Config Help  =
    !Config Units = [FLAG]
    irrigated_soiltile = .FALSE.
    CALL getin_p('IRRIGATED_SOILTILE',irrigated_soiltile)

    !
    !Config Key   = OLD_IRRIG_SCHEME
    !Config Desc  = Do we run with the old irrigation scheme?
    !Config If    = DO_IRRIGATION
    !Config Def   = n
    !Config Help  =
    !Config Units = [FLAG]
    old_irrig_scheme = .FALSE.  
    CALL getin_p('OLD_IRRIG_SCHEME',old_irrig_scheme)

    !   To run with new irrigation scheme
    !Config Key   = IRRIG_ST
    !Config Desc  = Which is the soil tile with irrigation flux
    !Config If    = DO_IRRIGATION and old_irrig_scheme = FALSE and IRRIGATED_SOILTILE = FALSE
    !Config Def   = 3
    !Config Help  =
    !Config Units = []
    irrig_st = 3 
    CALL getin_p('IRRIG_ST',irrig_st)

    !  To run with NEW irrigation scheme
    !Config Key   = AVAIL_RESERVE
    !Config Desc  = Max. fraction of routing reservoir volume that can be used for irrigation
    !Config If    = DO_IRRIGATION and old_irrig_scheme = FALSE
    !Config Def   = 0.9,0.9,0.9
    !Config Help  = Available water from routing reservoirs, to withdraw for irrigation
    !Config Help  = IMPORTANT: As the routing model uses 3 reservoirs, dimension is set to 3
    !Config Help  = IMPORTANT: Order of available water must be in this order: streamflow, fast, and slow reservoir
    !Config Units = []
    avail_reserve = (/0.9,0.9,0.9/)
    CALL getin_p('AVAIL_RESERVE',avail_reserve)

    !  To run with NEW irrigation scheme
    !Config Key   = BETA_IRRIG
    !Config Desc  = Threshold multiplier of Target SM to calculate root deficit
    !Config If    = DO_IRRIGATION and OLD_IRRIG_SCHEME = FALSE
    !Config Def   = 0.9
    !Config Help  =
    !Config Units = []
    beta_irrig = 0.9
    CALL getin_p('BETA_IRRIG',beta_irrig)
    !
    !   To run with new irrigation scheme
    !Config Key   = LAI_IRRIG_MIN
    !Config Desc  = Min. LAI (Leaf Area Index) to trigger irrigation
    !Config If    = DO_IRRIGATION and old_irrig_scheme = FALSE and IRRIGATED_SOILTILE = FALSE
    !Config Def   = 0.1
    !Config Help  =
    !Config Units = [m2/m2]
    lai_irrig_min = 0.1
    CALL getin_p('LAI_IRRIG_MIN',lai_irrig_min)

    !   To run with new irrigation scheme
    !Config Key   = IRRIG_MAP_DYNAMIC_FLAG
    !Config Desc  = Do we use a dynamic irrig map?
    !Config If    = DO_IRRIGATION
    !Config Def   = n
    !Config Help  =
    !Config Units = [FLAG]
    irrig_map_dynamic_flag = .FALSE.
    CALL getin_p('IRRIG_MAP_DYNAMIC_FLAG',irrig_map_dynamic_flag)

    !   To run with new irrigation scheme
    !Config Key   = SELECT_SOURCE_IRRIG
    !Config Desc  = Do we use the new priorization scheme, based on maps of equipped area with surface water?
    !Config If    = DO_IRRIGATION
    !Config Def   = n
    !Config Help  =
    !Config Units = [FLAG]
    select_source_irrig = .FALSE.
    CALL getin_p('SELECT_SOURCE_IRRIG',select_source_irrig)

    !Config Key   = REINFILTR_IRRIGFIELD
    !Config Desc  = Reinfiltrate all/part of runoff in crop soil tile
    !Config If    = DO_IRRIGATION
    !Config Def   = n
    !Config Help  =
    !Config Units = [FLAG]
    Reinfiltr_IrrigField = .FALSE.
    CALL getin_p('REINFILTR_IRRIGFIELD',Reinfiltr_IrrigField)

    !   Externalized reinf_slope for reinfiltration in irrigated cropland
    !Config Key   = reinf_slope_cropParam
    !Config Desc  = DO_IRRIGATION
    !Config If    = DO_IRRIGATION = T, REINFILTR_IRRIGFIELD = T
    !Config Def   = 0.8
    !Config Help  =
    !Config Units = [FRACTION 0-1]
    reinf_slope_cropParam = 0.8 
    CALL getin_p('REINF_SLOPE_CROPPARAM',reinf_slope_cropParam)

    !   Externalized for available volume to adduction
    !Config Key   = a_stream_adduction  
    !Config Desc  = DO_IRRIGATION
    !Config If    = DO_IRRIGATION = T
    !Config Def   = 0.05
    !Config Help  =
    !Config Units = [FRACTION 0-1]
    a_stream_adduction = 0.05
    CALL getin_p('A_STREAM_ADDUCTION',a_stream_adduction)


    !Surface resistance
    !
    !Config Key = CB
    !Config Desc = Constant of the Louis scheme 
    !Config If = OK_SECHIBA
    !Config Def = 5.0
    !Config Help = 
    !Config Units = [-]
    cb = 5._r_std  
    CALL getin_p('CB',cb)
    !
    !Config Key = CC
    !Config Desc = Constant of the Louis scheme 
    !Config If = OK_SECHIBA
    !Config Def = 5.0
    !Config Help = 
    !Config Units = [-]
    cc = 5._r_std
    CALL getin_p('CC',cc)
    !
    !Config Key = CD
    !Config Desc = Constant of the Louis scheme 
    !Config If = OK_SECHIBA
    !Config Def = 5.0
    !Config Help = 
    !Config Units = [-]
    cd = 5._r_std 
    CALL getin_p('CD',cd)
    !
    !Config Key = RAYT_CSTE
    !Config Desc = Constant in the computation of surface resistance  
    !Config If = OK_SECHIBA
    !Config Def = 125
    !Config Help = 
    !Config Units = [W.m^{-2}]
    rayt_cste = 125.
    CALL getin_p('RAYT_CSTE',rayt_cste)
    !
    !Config Key = DEFC_PLUS
    !Config Desc = Constant in the computation of surface resistance  
    !Config If = OK_SECHIBA
    !Config Def = 23.E-3
    !Config Help = 
    !Config Units = [K.W^{-1}]
    defc_plus = 23.e-3 
    CALL getin_p('DEFC_PLUS',defc_plus)
    !
    !Config Key = DEFC_MULT
    !Config Desc = Constant in the computation of surface resistance  
    !Config If = OK_SECHIBA
    !Config Def = 1.5
    !Config Help = 
    !Config Units = [K.W^{-1}]
    defc_mult = 1.5
    CALL getin_p('DEFC_MULT',defc_mult)


    ! Sets the levels for the energy calculations based on the
    ! ENERGY_CONTROL flag. 
    IF( ENERGY_CONTROL .LE. 3 ) THEN        ! single layer
       
       jnlvls = 1
       jnlvls_under = 0
       jnlvls_canopy = 1
       jnlvls_over = 0    
       
    ELSEIF( ENERGY_CONTROL .EQ. 4 ) THEN   ! multi-layer
       
       jnlvls = 29
       jnlvls_under = 10
       jnlvls_canopy = 10
        jnlvls_over = 9
        
     ELSEIF( ENERGY_CONTROL .EQ. 5 ) THEN  ! run.def specified
        !Config Key   = JNLVLS
        !Config Desc  = number of photosyntheis canopy levels 
        !Config If    = OK_SECHIBA
        !Config Def   = 29
        !Config Help  = The total number of layers in the 
        !Config       = multi-layer energy budget calculations
        !Config Units = [-]
        jnlvls = 29
        CALL getin_p("JNLVLS",jnlvls)
        
        !Config Key   = JNLVLS_UNDER
        !Config Desc  = number of energy layers under the canopy 
        !Config If    = OK_SECHIBA
        !Config Def   = 10
        !Config Help  = Number of energy layers under the canopy
        !Config       = used for multi-layer energy budget calculations
        !Config Units = [-]
        jnlvls_under = 10
        CALL getin_p("JNLVLS_UNDER",jnlvls_under)

        !Config Key   = JNLVLS_CANOPY
        !Config Desc  = number of energy layers in the canopy 
        !Config If    = OK_SECHIBA
        !Config Def   = 10
        !Config Help  = Number of energy, albedo and photosynthesis 
        !Config       = layers in the canopy used for multi-layer 
        !Confog       = energy budget calculations
        !Config Units = [-]
        jnlvls_canopy = 9
        CALL getin_p("JNLVLS_CANOPY",jnlvls_canopy)

        !Config Key   = JNLVLS_OVER
        !Config Desc  = number of energy layers over the canopy 
        !Config If    = OK_SECHIBA
        !Config Def   = 10
        !Config Help  = Number of energy layers over the canopy
        !Config       = used for multi-layer energy budget calculations
        !Config Units = [-]
        jnlvls_over = 10
        CALL getin_p("JNLVLS_OVER",jnlvls_over)
     ENDIF ! (ENERGY_CONTROL)

     IF (printlev>=1) WRITE(numout,*)'The set JNLVLS levels', jnlvls, jnlvls_under,jnlvls_canopy,jnlvls_over 

     !Config Key   = NLEV_TOP
     !Config Desc  = Maximum number of canopy levels at the top.
     
     !Config If    = OK_SECHIBA
     !Config Def   = 10
     !Config Help  = Maximum number of canopy levels that are 
     !               used to construct the "top" layer of the 
     !               canopy. The top layer is used in the 
     !               calculation transpiration. 
     !               Should not exceed nlai
     !Config Units = [-]
     nlev_top = 10
     CALL getin_p('NLEV_TOP',nlev_top)        
     IF(nlev_top .GT. nlai) THEN
        nlev_top = nlai
        IF (printlev>=1) WRITE(numout,*) 'The numbers of levels in the "top" was '
        IF (printlev>=1) WRITE(numout,*) '  larger than the total number of levels'
        IF (printlev>=1) WRITE(numout,*) '  AUTO CORRECT: nlev_top now = ',nlev_top
     ENDIF
     !+++++++++++

     !
     ! diffuco 
     !
     !Config Key   = LAIMAX
     !Config Desc  = Maximum LAI
     !Config If    = OK_SECHIBA
     !Config Def   = 
     !Config Help  = 
     !Config Units = [m^2/m^2]
     laimax = 12.
     CALL getin_p('LAIMAX',laimax)
     !
     !Config Key   = DEW_VEG_POLY_COEFF
     !Config Desc  = coefficients of the polynome of degree 5 for the dew
     !Config If    = OK_SECHIBA
     !Config Def   = 0.887773, 0.205673, 0.110112, 0.014843, 0.000824, 0.000017 
     !Config Help  = 
     !Config Units = [-]
     dew_veg_poly_coeff = (/ 0.887773, 0.205673, 0.110112, 0.014843,  0.000824,  0.000017 /)
     CALL getin_p('DEW_VEG_POLY_COEFF',dew_veg_poly_coeff)
     !
     ! DOWNREGULATION_CO2 is deactivated in ORCHIDEE. 
     ! The code is left in because related to vcmax which is kept as a diagnostic variable.
!     !Config Key   = DOWNREGULATION_CO2
!     !Config Desc  = Activation of CO2 downregulation.
!     !Config If    = OK_SECHIBA
!     !Config Def   = n
!     !Config Help  = 
!     !Config Units = [FLAG]
     downregulation_co2 = .FALSE.
!     CALL getin_p('DOWNREGULATION_CO2',downregulation_co2)
     
!     !Config Key   = DOWNREGULATION_CO2_BASELEVEL
!     !Config Desc  = CO2 base level
!     !Config If    = OK_SECHIBA 
!     !Config Def   = 380.
!     !Config Help  = 
!     !Config Units = [ppm]
     downregulation_co2_baselevel = 380.
!     CALL getin_p('DOWNREGULATION_CO2_BASELEVEL',downregulation_co2_baselevel)
     
     !Config Key   = GB_REF
     !Config Desc  = Leaf bulk boundary layer resistance
     !Config If    = 
     !Config Def   = 1./25.
     !Config Help  = 
     !Config Units = [s m-1]
     gb_ref = 1./25.
     CALL getin_p('GB_REF',gb_ref)
     !
     !Config Key   = BULK_DEFAULT
     !Config Desc  = default bulk density
     !Config If    = OK_SECHIBA 
     !Config Def   = 1000.0
     !Config Help  = The bulk density is the weight of soil in a
     !Config         given volume.  This default is used if no other value
     !Config         is found in the restart file.
     !Config Units = [kg/m3]
     bulk_default = 1000.0  
     CALL getin_p('BULK_DEFAULT',bulk_default)
     !
     !Config Key   = PH_DEFAULT
     !Config Desc  = default soil pH
     !Config If    = OK_SECHIBA 
     !Config Def   = 5.5
     !Config Help  = Gives the value of the soil pH if a value is not
     !Config         found in the restart file.
     !Config Units = [-]
     ph_default = 5.5 
     CALL getin_p('PH_DEFAULT',ph_default)
     !
     !Config Key   = MIN_VEGFRAC 
     !Config Desc  = Minimal fraction of mesh a vegetation type can occupy 
     !Config If    = OK_SECHIBA 
     !Config Def   = 0.001 
     !Config Help  = 
     !Config Units = [-]
     min_vegfrac = 0.001 
     CALL getin_p('MIN_VEGFRAC',min_vegfrac)

    IF (min_vegfrac.LT.min_stomate) THEN
       
       ! In slowproc_adjust_delta_veget_max several series of IF-statements implicitly
       ! assume that min_vegfarc > min_stomate. If this is not the case strange things
       ! may happen. These strange things woulD probably surface as a mass balance problem
       ! or a problem with conserving the surface areas of a pixel. It is a very easy
       ! inconsistency to catch at this point. It may be much harder to understand 
       ! what is going wrong if the model is used with an inconstent min_vegfrac value.
       WRITE(numout,*) 'min_vegfrac and min_stomate ', min_vegfrac, min_stomate
       CALL ipslerr_p(3,'Check your run.def','min_vegfrac should be larger than min_stomate',&
            'If not several of the IF-statements in slowproc_adjust_delta_veget_max',&
            'may not work as intended')

    END IF
    !
    !Config Key   = FRAC_NOBIO_FIXED_TEST1 
    !Config Desc  = Value for frac_nobio for tests in 0-dim simulations (0-1, unitless)
    !Config If    = OK_SECHIBA 
    !Config Def   = 0.0
    !Config Help  = 
    !Config Units = [-]
    frac_nobio_fixed_test_1 = 0.0
    CALL getin_p('FRAC_NOBIO_FIXED_TEST1',frac_nobio_fixed_test_1)
    !
    !Config Key   = STEMPDIAG_BID 
    !Config Desc  = only needed for an initial LAI if there is no restart file
    !Config If    = OK_SECHIBA 
    !Config Def   = 280.
    !Config Help  = 
    !Config Units = [K]
    stempdiag_bid = 280. 
    CALL getin_p('STEMPDIAG_BID',stempdiag_bid)
 
    !
    !Config Key   = MIN_N
    !Config Desc  = Minimum allowable n_mineralisation in som_dynamics
    !Config If    = OK_STOMATE
    !Config Def   = 0.0001
    !Config Help  =
    !Config Units = gNH4-N/m^2/day
    min_n = 0.0001
    CALL getin_p('MIN_N',min_n)    
    !
    !Config Key   = MAX_CN
    !Config Desc  = Maximum allowable ratio of som_input_total(:,icarbon)
    !               to som_input_total(:,initrogen).
    !Config If    = OK_STOMATE
    !Config Def   = 250
    !Config Help  =
    !Config Units = [-]
    max_cn = 250 
    CALL getin_p('MAX_CN',max_cn)
    !
    !Config Key   = SNC
    !Config Desc  = Structural nitrogen concentration
    !Config If    = OK_STOMATE
    !Config Def   = 0.004
    !Config Help  = Structural nitrogen [gN gC-1] based on C:N of dead wood 
    !               (White et al., 2000) assuming carbon dry matter ratio of 0.5
    !Config Units = [gN gC-1]
    snc = 0.004 
    CALL getin_p('SNC',snc)

    !Config Key   = SUGAR_LOAD_MIN
    !Config Desc  = Lower bound for sugar loading when used to regulate NUE
    !Config If    = OK_STOMATE
    !Config Def   = 0.0
    !Config Help  = Sugar loading is a ratio that will be contained between
    !               sugar_load_min and sugar_load_max. If it is 1 all is perfect
    !               Lower bound for sugar loading when used to regulate NUE. 
    !               Sugar loafding results in a strong reduction of GPP. By
    !               taking a rather high lower limit the impact of sugar loading
    !               is limited. This will result in high C-reserves unless 
    !               leaching is implemented. 
    !Config Units = [-]
    sugar_load_min = 0.05 
    CALL getin_p('SUGAR_LOAD_MIN',sugar_load_min)
    !
    !Config Key   = SUGAR_LOAD_MAX
    !Config Desc  = Upper bound for sugar loading when used to regulate NUE
    !Config If    = OK_STOMATE
    !Config Def   = 1.0
    !Config Help  = Sugar loading is a ratio that will be contained between
    !               sugar_load_min and sugar_load_max. If it is 1 all is perfect
    !Config Units = [-]
    sugar_load_max = 1.0
    CALL getin_p('SUGAR_LOAD_MAX',sugar_load_max)
    !-
    !Config Key   = NCIRC
    !Config Desc  = Number of basal area classes in allocation scheme
    !               circ classes could be considered as cohorts within a stand 
    !Config If    = OK_STOMATE, OK_SECHIBA 
    !Config Def   = 2 
    !Config Help  = 
    !Config Units = [-]
    ncirc = 1 
    CALL getin_p('NCIRC',ncirc)
    !
    !Config Key   = SLOPE_RA
    !Config Desc  = Reduction factor to make resp_maint less temperature sensitive
    !Config If    = OK_STOMATE
    !Config Def   = 1.
    !Config Help  =
    !Config Units = [-]
    slope_ra = 40
    CALL getin_p('SLOPE_RA',slope_ra)
    !
    !Config Key   = LAIEFF_SOLAR_ANGLE
    !Config Desc  = The solar zenith angle for effective LAI
    !Config If    = OK_SECHIBA
    !Config Def   = 60
    !Config Help  = The solar zenith angle used in the calculation of the 
    !               effective LAI. Pinty recommends a value of 60 degrees.
    !Config Units = [degrees]
    laieff_solar_angle=60.0_r_std
    CALL getin_p('LAIEFF_SOLAR_ANGLE',laieff_solar_angle) 
    IF (printlev>=1) WRITE(numout,*) 'Maximum number of optimization steps'//&
         'tried for albedo n-layer optimization [degrees]: ',laieff_solar_angle
    ! convert to radians
    laieff_solar_angle=laieff_solar_angle/180.0_r_std*pi 
    !
    !Config Key   = LAIEFF_ZERO_CUTOFF
    !Config Desc  = Cutoff for effective lai values
    !Config If    = OK_SECHIBA
    !Config Def   = 0.0000001
    !Config Help  = This is an arbitrary cutoff to make sure we don't pass 
    !               zero values of crown diameter and trunk diameter
    !               to a subroutine that will choke on them
    !Config Units = [-]
    laieff_zero_cutoff=0.0000001_r_std
    CALL getin_p('LAIEFF_ZERO_CUTOFF',laieff_zero_cutoff) 
    IF (printlev>=1) WRITE(numout,*) 'A minimum threshold for trunk and crown '//&
         'diameters for numerical stability: ',laieff_zero_cutoff
    !
    !Config Key   = DIRECT_LIGHT_WEIGHT
    !Config Desc  = The weighting factor to weight different sources of light
    !Config If    = OK_SECHIBA
    !Config Def   = 0.5
    !Config Help  = The weighting factor used in the calculation of the
    !               albedo and canopy absorbed light from different sources.
    !               Currently LDMZ and forcing files don't have information
    !               on the amount of solar radiation hitting the canopy 
    !               directly from the sun (collimated) or first reflected off
    !               clouds and aerosols.  But these values are calculated in
    !               the albedo routines.  We combine them into a single
    !               value with a simple weighting controlled by this variable.
    !               We select a value of 0.5 because we don't know what will
    !               be generally applicable.
    !Config Units = [degrees]
    direct_light_weight=0.5_r_std
    CALL getin_p('direct_light_weight',direct_light_weight) 
    IF (printlev>=1) WRITE(numout,*) 'Weighting fraction for calculating the total'//&
         'albedo and absorbed light from direct and diffuse sources [-]: ',direct_light_weight

    !Config Key   = MAINT_RESP_CONTROL
    !Config Desc  = Sets the approach to calculate Rm
    !Config If    = OK_SECHIBA
    !Config Def   = 'trunk'
    !Config Help  =  Choose between three options to calculate the maint respiration. If
    !                set to 'nitrogen' the plant nitrogen pools and temperature will 
    !                drive maint_resp. If set to 'cn', maint_resp will be adjusted for 
    !                the expected c/n ratio of the plant biomass. Also the temperature 
    !                control itself is calculated according to Krinner et al 2005 instead
    !                of Sitch et al 2003. This is the exact approach and parameters used
    !                in orchidee 3. If set to 'trunk' the approach is very similar to 'cn'
    !                except for the Rm when T is below zero and how labile and reserve 
    !                carbon contribute to the Rm. The largest difference between 'cn'
    !                and 'trunk' is in the parameter values.
    !Config Units = [-]
    maint_resp_control='trunk'
    CALL getin_p('MAINT_RESP_CONTROL',maint_resp_control) 
    IF (printlev>=1) WRITE(numout,*) 'Maintenance respiration is based on ',maint_resp_control
    !
    !Config Key   = CROWN_PACKING
    !Config Desc  = Packing efficiency of the crowns within the canopy space
    !Config If    = OK_SECHIBA
    !Config Def   = 1.
    !Config Help  = The crowns are assumed to be ellipsoids. The sum of all the
    !               individual crown volumes may exceed the entire canopys space
    !               given by the surface*height. If this happens the crown dimensions
    !               (the axes of the ellipsoids) will be recalculated. This 
    !               calculation assume a packing efficiency. Note that close-packing
    !               of spheres has a maximum efficiency of 0.74 (listed as some mental 
    !               guidance). Lower packing efficiencies will result in more gaps and
    !               has consequences for all variables that rely on veget.
    !Config Units = [-]
    crown_packing = 1.0 
    CALL getin_p('CROWN_PACKING',crown_packing)

    !
    !Config Key   = ALPHA_HEIGHT_PRECIP
    !Config Desc  = alpha parameter in height-precipitation relationship
    !Config If    = 
    !Config Def   = 3.22
    !Config Help  = alpha parameter in height-precipitation relationship
    !               Fitted on the GEDI and CRUJRA data by using
    !               svn/TOOLS/FIT_TREEHEIGHT_PRECIPITATION/dynamic_tree_height.py
    !Config Units = [m]
    alpha_height_precip = 3.22 
    CALL getin_p('ALPHA_HEIGHT_PRECIP',alpha_height_precip)
    !
    !Config Key   = BETA_HEIGHT_PRECIP
    !Config Desc  = beta parameter in height-precipitation relationship
    !Config If    = 
    !Config Def   = 0.33877
    !Config Help  = alpha parameter in height-precipitation relationship
    !               Fitted on the GEDI and CRUJRA data by using
    !               svn/TOOLS/FIT_TREEHEIGHT_PRECIPITATION/dynamic_tree_height.py
    !Config Units = [-]
    beta_height_precip = 0.33877
    CALL getin_p('BETA_HEIGHT_PRECIP',beta_height_precip)

    !Config Key   = N_ITER_TUZET
    !Config Desc  = number of iterations of the simili-implicit scheme in Tuzet
    !Config If    = OK_HYDROL_ARCH
    !Config Def   = 25
    !Config Help  =
    !Config Units = unitless
    n_iter_tuzet = 25 
    CALL getin_p('N_ITER_TUZET',n_iter_tuzet)
    

    

    !-
    ! lake
    !-
    IF (ok_lake_energy) THEN
       !Config Key   = NLAKE
       !Config Desc  = Maximum number of lake type in a grid cell
       !Config If    = OK_LAKE_ENERGY
       !Config Def   = 3
       !Config Help  = 
       !Config Units = [-]
       CALL getin_p('NLAKE',nlake)
              
       !Config Key   = ALBEDO_WATER_LAKE
       !Config Desc  = free water lake albedo
       !Config If    = OK_LAKE_ENERGY 
       !Config Def   = 0.07
       !Config Help  = 
       !Config Units = [-]
       CALL getin_p('ALBEDO_WATER_LAKE',albedo_water_lake)
       
       !Config Key   = LIGHT_EXT_SNOW_LAKE
       !Config Desc  = snow lake light extinction
       !Config If    = OK_LAKE_ENERGY
       !Config Def   = 15
       !Config Help  = 
       !Config Units = [m-1]
       CALL getin_p('LIGHT_EXT_SNOW_LAKE', light_ext_snow_lake)
       
       !Config Key   = LIGHT_EXT_ICE_LAKE
       !Config Desc  = ice lake light extinction
       !Config If    = OK_LAKE_ENERGY
       !Config Def   = 8
       !Config Help  = 
       !Config Units = [m-1]
       CALL getin_p('LIGHT_EXT_ICE_LAKE', light_ext_ice_lake)
       
       !Config Key   = HEIGHT_TQ_LAKE
       !Config Desc  = Height where temperature and humidity are measured
       !Config If    = OK_LAKE_ENERGY
       !Config Def   = 2
       !Config Help  = 
       !Config Units = [m]
       CALL getin_p('HEIGHT_TQ_LAKE',height_tq_lake)
       
       !Config Key   = HEIGHT_U_LAKE
       !Config Desc  = snow lake light extinction
       !Config If    = OK_LAKE_ENERGY 
       !Config Def   = 2
       !Config Help  = 
       !Config Units = [m]
       CALL getin_p('HEIGHT_U_LAKE', height_u_lake)
       
       !Config Key   = T_CLIM_LAKE 
       !Config Desc  = Climate temperature for initialization if constant method is used
       !Config If    = OK_LAKE_ENERGY
       !Config Def   = 277.15
       !Config Help  = 
       !Config Units = [K]
       CALL getin_p('T_CLIM_LAKE', T_clim_lake)

       !Config Key   = LAKE_DEPTH_CONS
       !Config Desc  = Uniform constant value for lake depth instead of spatial values from file 
       !Config If    = OK_LAKE_ENERGY
       !Config Def   = 0
       !Config Help  = This value is used only if different from zero
       !Config Units = [m]
       CALL getin_p('LAKE_DEPTH_CONS', lake_depth_cons)
       
       !Config Key   = OK_DEPTH_LAKE_MAX
       !Config Desc  = Truncation of the depth lake for FLake use (Should be used)
       !Config If    = OK_LAKE_ENERGY
       !Config Def   = y 
       !Config Help  = 
       !Config Units = [-]
       CALL getin_p('OK_DEPTH_LAKE_MAX', ok_depth_lake_max)
       
       IF (ok_depth_lake_max) THEN
          ! Coherence test: ok_depth_lake_max and lake_depth_cons can not be at the same time
          IF (lake_depth_cons > 0)  CALL ipslerr_p (3,'constantes.f90', &
               'OK_DEPTH_LAKE_MAX=y and LAKE_DEPTH_CONS>0 are incompatible', &
               'Change run.def for one of these options','')       

          !Config Key   = DEPTH_LAKE_MAX
          !Config Desc  = Maximul lake depth
          !Config If    = OK_LAKE_ENERGY and OK_DEPTH_LAKE_MAX
          !Config Def   = 60 
          !Config Help  = 
          !Config Units = [m]
          CALL getin_p('DEPTH_LAKE_MAX', depth_lake_max)
       END IF


       IF (.NOT. ok_snow_ice_lake_par) THEN
          !Config Key   = ALBEDO_ICE_LAKE
          !Config Desc  = ice lake albedo
          !Config If    = OK_LAKE_ENERGY and OK_SNOW_ICE_LAKE_PAR=n
          !Config Def   = 0.4
          !Config Help  = 
          !Config Units = [-]
          CALL getin_p('ALBEDO_ICE_LAKE', albedo_ice_lake)
          
          !Config Key   = ALBEDO_SNOW_LAKE
          !Config Desc  = snow lake albedo
          !Config If    = ok_lake_energy and OK_SNOW_ICE_LAKE_PAR=n
          !Config Def   = 0.6
          !Config Help  = 
          !Config Units = [-]
          CALL getin_p('ALBEDO_SNOW_LAKE', albedo_snow_lake)
       END IF
       
       IF (ok_botsed_lake) THEN
          !Config Key   = DEPTH_SED_LAKE
          !Config Desc  = depth layer of sediment layer
          !Config If    = OK_BOTSED_LAKE
          !Config Def   = 10
          !Config Help  = 
          !Config Units = [m]
          CALL getin_p('DEPTH_SED_LAKE', depth_sed_lake)
                    
          !Config Key   = T_SED_LAKE
          !Config Desc  = Sediment bottom temperature
          !Config If    = OK_BOTSED_LAKE
          !Config Def   = 277.15
          !Config Help  = 
          !Config Units = [K]
          CALL getin_p('T_SED_LAKE', t_sed_lake)
       END IF
    END IF

  END SUBROUTINE config_sechiba_parameters


!! ================================================================================================================================
!! SUBROUTINE   : config_co2_parameters 
!!
!>\BRIEF        This subroutine reads in the configuration file all the parameters
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): None 
!!
!! REFERENCE(S) :
!!
!! FLOWCHART    :
!! \n
!_ ================================================================================================================================

  SUBROUTINE config_co2_parameters

    IMPLICIT NONE

    !! 0. Variables and parameters declaration

    !! 0.4 Local variables

    !_ ================================================================================================================================

    !
    !Config Key   = LAI_LEVEL_DEPTH
    !Config Desc  = 
    !Config If    = 
    !Config Def   = 0.15
    !Config Help  = 
    !Config Units = [-]
    lai_level_depth = 0.15 
    CALL getin_p('LAI_LEVEL_DEPTH',lai_level_depth)
    !
    !Config Key   = Oi
    !Config Desc  = Intercellular oxygen partial pressure
    !Config If    = 
    !Config Def   = 210000.
    !Config Help  = See Legend of Figure 6 of Yin et al. (2009)
    !Config Units = [ubar]
    Oi=210000.
    CALL getin_p('Oi',Oi)

    !Config Key   = THRESHOLD_C13_ASSIM
    !Config Desc  = If assimilation falls below this threshold the delta_c13 is set to zero 
    !Config If    = OK_C13
    !Config Def   = 0.01 
    !Config Help  = 
    !Config Units = [-]
    threshold_c13_assim = 0.01
    CALL getin_p('THRESHOLD_C13_ASSIM',threshold_c13_assim)
    !
    !Config Key   = C13_A
    !Config Desc  = Coefficient for fractionation occurring due to diffusion in air 
    !Config If    = OK_C13
    !Config Def   = 0.01 
    !Config Help  = 
    !Config Units = [-]
    c13_a = 4.4 
    CALL getin_p('C13_A',c13_a)
    !
    !Config Key   = C13_B
    !Config Desc  = Coefficient for fractionation caused by carboxylation 
    !Config If    = OK_C13
    !Config Def   = 0.01 
    !Config Help  = 
    !Config Units = [-]
    c13_b = 27.
    CALL getin_p('C13_B',c13_b)

  END SUBROUTINE config_co2_parameters


!! ================================================================================================================================
!! SUBROUTINE   : config_stomate_parameters 
!!
!>\BRIEF        This subroutine reads in the configuration file all the parameters 
!! needed when stomate is activated (ie : when OK_STOMATE is set to true).
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCE(S) :
!!
!! FLOWCHART    :
!! \n
!_ ================================================================================================================================

  SUBROUTINE config_stomate_parameters

    IMPLICIT NONE

    !! 0. Variables and parameters declaration

    !! 0.4 Local variables
    LOGICAL                          :: l_error     !! Check errors in allocation
    INTEGER(i_std)                   :: ier         !! Check errors in allocation
    INTEGER(i_std)                   :: iout        !! Index
    CHARACTER(LEN=2)                 :: striout     !! string suffix indicating output diameter class
    INTEGER(i_std)                   :: jv_loc      !! PFT index (AGE_CLASS_BOUNDS_PFT)
    INTEGER(i_std)                   :: ic_loc      !! Age class index (AGE_CLASS_BOUNDS_PFT)
    REAL(r_std)                      :: acifsum     !! Sum of AGE_CLASS_INIT_FRAC, used to renormalise it
    REAL(r_std), ALLOCATABLE, DIMENSION(:) :: age_class_bound_glo !! Global bounds read from the namelist, broadcast to every PFT

    !_ ================================================================================================================================

    l_error = .FALSE.

    !Config Key   = EXP_KF
    !Config Desc  = Exponential of the sensitivity of k_latosa to tree height
    !Config If    = OK_STOMATE
    !Config Def   = 3.0
    !Config Help  = Exponential used to tune the sensitivity of k_latosa to tree 
    !               height. 
    !Config Units = [-]
    exp_kf = 3.0
    CALL getin_p('EXP_KF',exp_kf)

    !-
    ! constraints_parameters
    !-
    !
    !Config Key   = TOO_LONG 
    !Config Desc  = longest sustainable time without regeneration (vernalization)
    !Config If    = OK_STOMATE
    !Config Def   = 5.
    !Config Help  = 
    !Config Units = [years]
    too_long = 5.
    CALL getin_p('TOO_LONG',too_long)

    !-
    ! fire parameters
    !-
    !
    !Config Key   = TAU_FIRE 
    !Config Desc  = Time scale for memory of the fire index (days). Validated for one year in the DGVM. 
    !Config If    = OK_STOMATE 
    !Config Def   = 30.
    !Config Help  = 
    !Config Units = [days]
    tau_fire = 30.
    CALL getin_p('TAU_FIRE',tau_fire)
    !
    !Config Key   = CR_FIRE
    !Config Desc  = VPD fire-danger precip-suppression coefficient: FDI = VPD*exp(-CR_FIRE*precip_month). Was hardcoded 2.0; exposed as a calibration knob.
    !Config If    = OK_SPITFIRE
    !Config Def   = 2.0
    !Config Help  = Lower -> less rain suppression -> higher FDI -> more fire. Compensates the scale-dependence of the 0.5deg/weekly VPD.
    !Config Units = [-]
    cr_fire = 2.0
    CALL getin_p('CR_FIRE',cr_fire)
    !
    !Config Key   = A_ND_SCALE
    !Config Desc  = Multiplicative scaling of the human-ignition parameter a_nd (human_ign = exp(-0.5*sqrt(popd))*a_nd*A_ND_SCALE*popd). Calibration knob.
    !Config If    = OK_SPITFIRE
    !Config Def   = 1.0
    !Config Help  = Higher -> more human ignitions -> more fire. Absorbs the scale-dependence of anthropogenic ignition density (0.5deg/daily).
    !Config Units = [-]
    a_nd_scale = 1.0
    CALL getin_p('A_ND_SCALE',a_nd_scale)
    !
    !Config Key   = LITTER_CRIT
    !Config Desc  = Critical litter quantity for fire
    !Config If    = OK_STOMATE 
    !Config Def   = 200.
    !Config Help  = 
    !Config Units = [gC/m^2]
    litter_crit = 200. 
    CALL getin_p('LITTER_CRIT',litter_crit)

    !Config Key   = FIRE_RESIST_LIGNIN
    !Config Desc  = 
    !Config If    = OK_STOMATE 
    !Config Def   = 0.5
    !Config Help  = 
    !Config Units = [-]
    fire_resist_lignin = 0.5
    CALL getin_p('FIRE_RESIST_LIGNIN',fire_resist_lignin)

    !Config Key   = CO2FRAC
    !Config Desc  = What fraction of a burned plant compartment goes into the atmosphere
    !Config If    = OK_STOMATE 
    !Config Def   = 0.95, 0.95, 0., 0.3, 0., 0., 0.95, 0.95
    !Config Help  = 
    !Config Units = [-]
    co2frac = (/ .95, .95, 0., 0.3, 0., 0., .95, .95, .95 /) 
    CALL getin_p('CO2FRAC',co2frac)
    !
    !Config Key   = BCFRAC_COEFF
    !Config Desc  = 
    !Config If    = OK_STOMATE 
    !Config Def   = 0.3, 1.3, 88.2 
    !Config Help  = 
    !Config Units = [-]
    bcfrac_coeff = (/ .3,  1.3,  88.2 /) 
    CALL getin_p('BCFRAC_COEFF',bcfrac_coeff)
    
    !
    !Config Key   = ALLOC_FIREFUEL
    !Config Desc  = 
    !Config If    = OK_STOMATE 
    !Config Def   = 0.045,0.075,0.21,0.67
    !Config Help  = 
    !Config Units = [-]  
    CALL getin_p('ALLOC_FIREFUEL',alloc_firefuel)
    
    !Config Key   = FIREFRAC_COEFF 
    !Config Desc  = 
    !Config If    = OK_STOMATE 
    !Config Def   = 0.45, 0.8, 0.6, 0.13
    !Config Help  = 
    !Config Units = [-]
    firefrac_coeff = (/ 0.45, 0.8, 0.6, 0.13 /)
    CALL getin_p('FIREFRAC_COEFF',firefrac_coeff)

    !Config Key   = REF_GREFF
    !Config Desc  = Asymptotic maximum mortality rate
    !Config If    = OK_STOMATE 
    !Config Def   = 0.035
    !Config Help  = Set asymptotic maximum mortality rate from Sitch 2003
    !Config         (they use 0.01) (year^{-1})
    !Config Units = [1/year]
    ref_greff = 0.035
    CALL getin_p('REF_GREFF',ref_greff)

    !Config Key   = OK_AED_WIND
    !Config Desc  = Flag for opening fragmentation effects on wind
    !Config If    = .NOT. FIRE_DISABLE
    !Config Def   = n
    !Config Help  = fragmentation wind effect
    !Config Units = [FLAG]
    ok_aed_wind = .FALSE.
    CALL getin_p('OK_AED_WIND', ok_aed_wind)
    !
    !Config Key   = OK_AED_SIZE
    !Config Desc  = Flag for opening fragmentation effects on fire size
    !Config If    = .NOT. FIRE_DISABLE
    !Config Def   = n
    !Config Help  = flag that allows fragmentation fire size effect
    !Config Units = [FLAG]
    ok_aed_size = .FALSE.
    CALL getin_p('OK_AED_SIZE', ok_aed_size)
    !
    !Config Key   = OK_AED_HUMIGN
    !Config Desc  = Flag for opening fragmentation effects on human ignition
    !Config If    = .NOT. FIRE_DISABLE
    !Config Def   = n
    !Config Help  = flag that allows fragmentation human ignition effect
    !Config Units = [FLAG]
    ok_aed_humign = .FALSE.
    CALL getin_p('OK_AED_HUMIGN', ok_aed_humign)
    !
    !Config Key   = OK_AED_FUEL
    !Config Desc  = Flag for opening fragmentation effects on feul wetness
    !Config If    = .NOT. FIRE_DISABLE
    !Config Def   = n
    !Config Help  = flag that allows fragmentation fuel wetness effect
    !Config Units = [FLAG]
    ok_aed_fuel = .FALSE.
    CALL getin_p('OK_AED_FUEL', ok_aed_fuel)
    !
    !Config Key   = CF_TRESHOLD
    !Config Desc  = Crown fire happens when ignited_fraction is higher than cf_threshold
    !Config If    = 
    !Config Def   = 0.5 (Could be optimized)
    !Config Help  = 
    !Config Units = [-] 
    cf_threshold = 0.5  
    CALL getin_p('CF_TRESHOLD',cf_threshold)
    !
    !Config Key   = PATCH_MULTIPLIER_GRASS
    !Config Desc  = Multiplier to get max patch size for grass
    !Config If    = 
    !Config Def   = 1.2 (Bowring et al.. 2024, Table 1)
    !Config Help  = 
    !Config Units = [-] 
    patch_multiplier_grass = 1.2
    CALL getin_p('PATCH_MULTIPLIER_GRASS',patch_multiplier_grass)
    !
    !Config Key   = PATCH_MULTIPLIER_FOREST
    !Config Desc  = Multiplier to get max patch size for forest
    !Config If    = 
    !Config Def   = 1 (Bowring et al.. 2024, Table 1)
    !Config Help  = 
    !Config Units = [-]  
    patch_multiplier_forest = 1
    CALL getin_p('PATCH_MULTIPLIER_FOREST',patch_multiplier_forest)
    !
    !Config Key   = WETNESS_DIFF
    !Config Desc  = Relative fuel wetness difference between edge and interior
    !Config If    = 
    !Config Def   = 0.25 (Bowring et al.. 2024)
    !Config Help  = 
    !Config Units = [-]  
    wetness_diff = 0.25 
    CALL getin_p('WETNESS_DIFF',wetness_diff)
    !
    !Config Key   = EDGE_IGNITION_FACTOR
    !Config Desc  = Relative fuel wetness difference between edge and interior
    !Config If    = 
    !Config Def   = 1.0 (Bowring et al.. 2024)
    !Config Help  = 
    !Config Units = [-]  
    edge_ignition_factor = 1.0
    CALL getin_p('EDGE_IGNITION_FACTOR',edge_ignition_factor)
    !
    !Config Key   = EDGE_WIND_FACTOR
    !Config Desc  = Relative fuel wetness difference between edge and interior
    !Config If    = 
    !Config Def   = 1.0 (Bowring et al.. 2024)
    !Config Help  = 
    !Config Units = [-]  
    edge_wind_factor = 1.0
    CALL getin_p('EDGE_WIND_FACTOR',edge_wind_factor)
    !
    !Config Key   = EDGE_DISTANCE_IGNITION
    !Config Desc  = The edge distance for increasing human ignition
    !Config If    = 
    !Config Def   = 1.0 (Bowring et al.. 2024, Table 1)
    !Config Help  = 
    !Config Units = [m]  
    edge_distance_ignition = 1.0
    CALL getin_p('EDGE_DISTANCE_IGNITION',edge_distance_ignition)
    !
    !Config Key   = EDGE_DISTANCE_WIND
    !Config Desc  = The edge distance for increasing wind speed
    !Config If    = 
    !Config Def   = 16.0 (Bowring et al.. 2024, Table 1)
    !Config Help  = 
    !Config Units = [m]  
    edge_distance_wind = 16.0
    CALL getin_p('EDGE_DISTANCE_WIND',edge_distance_wind)
    !
    !Config Key   = EDGE_DISTANCE_WETNESS
    !Config Desc  = The edge distance for decreasing feul moisture
    !Config If    =
    !Config Def   = 20.0 (Bowring et al.. 2024, Table 1)
    !Config Help  =
    !Config Units = [m]
    edge_distance_wetness = 20.0
    CALL getin_p('EDGE_DISTANCE_WETNESS',edge_distance_wetness)
    !
    !Config Key   = OK_AED_LIGHT
    !Config Desc  = Flag for opening fragmentation effects on canopy light regime (Pgap, LAIeff)
    !Config If    = OK_STOMATE
    !Config Def   = n
    !Config Help  = Activates the edge effect on canopy porosity (Pgap +alpha) and effective LAI
    !Config         (LAIeff *(1+beta), beta<0). Non-directional formulation. Requires Edge_length.nc
    !Config         input (same file as the SPITFIRE fragmentation flags).
    !Config         [CATEGORIE : OPTION SCIENTIFIQUE -- choix de modelisation, permanent]
    !Config Units = [FLAG]
    ok_aed_light = .FALSE.
    CALL getin_p('OK_AED_LIGHT', ok_aed_light)
    !
    !Config Key   = EDGE_PGAP_FACTOR
    !Config Desc  = Relative change of Pgap at the edge vs interior
    !Config If    = OK_AED_LIGHT
    !Config Def   = 0.225 (Marie 2026, edge_effect.pdf slide 16)
    !Config Help  = Multiplies (1 + pgap_edge_factor * fedge_pgap(AED_light)) onto Pgap_cumul
    !Config         and Pgap_perlevel computed by effective_lai. Positive: edge canopies are more porous.
    !Config Units = [-]
    pgap_edge_factor = 0.225
    CALL getin_p('EDGE_PGAP_FACTOR', pgap_edge_factor)
    !
    !Config Key   = EDGE_LAIEFF_FACTOR
    !Config Desc  = Relative change of LAIeff at the edge vs interior
    !Config If    = OK_AED_LIGHT
    !Config Def   = -0.36 (Marie 2026, edge_effect.pdf slide 16)
    !Config Help  = Multiplies (1 + laieff_edge_factor * fedge_pgap(AED_light)) onto laieff
    !Config         computed by effective_lai. Negative: edge canopies have lower effective LAI.
    !Config Units = [-]
    laieff_edge_factor = -0.36
    CALL getin_p('EDGE_LAIEFF_FACTOR', laieff_edge_factor)
    !
    !Config Key   = EDGE_DISTANCE_LIGHT
    !Config Desc  = The edge distance for the canopy light effect
    !Config If    = OK_AED_LIGHT
    !Config Def   = 20.0 (Marie 2026 ; Laurance 2002 ; Schmidt 2017)
    !Config Help  = Penetration depth of the edge perturbation on canopy light regime.
    !Config         Used in the shared fragmentation form fedge_X = 1 - (AED-edge_distance_X/2)^2/AED^2.
    !Config Units = [m]
    edge_distance_light = 20.0
    CALL getin_p('EDGE_DISTANCE_LIGHT', edge_distance_light)
    !
    !Config Key   = OK_AED_FEEDBACK
    !Config Desc  = Flag for perturbation -> edge_length feedback (Marie 2026)
    !Config If    = OK_STOMATE
    !Config Def   = n
    !Config Help  = Activates a dynamic edge_length updated each timestep by
    !Config         fire/storm/pest-induced area perturbations, with exponential
    !Config         recovery toward the static Edge_length.nc reference. The dynamic
    !Config         edge_length then feeds AED_fire / AED_light via calculate_aed.
    !Config         [CATEGORIE : OPTION SCIENTIFIQUE -- choix de modelisation, permanent]
    !Config Units = [FLAG]
    ok_aed_feedback = .FALSE.
    CALL getin_p('OK_AED_FEEDBACK', ok_aed_feedback)
    !
    !Config Key   = EDGE_FEEDBACK_ALPHA_FIRE
    !Config Desc  = Edge length generated per m² burnt
    !Config If    = OK_AED_FEEDBACK
    !Config Def   = 0.0237 (moyenne pan-EU de la carte EFDA, 2026-07-30)
    !Config Help  = Linear sensitivity of edge_length to burnt area. Calibrated
    !Config         against the European Forest Disturbance Atlas (Senf et al. 2025)
    !Config         on AL/IT/FR/DE/SE (1985-2023): per-agent perimeter-to-area ratio
    !Config         dEL/dA for fire patches. Country-by-country range 0.018-0.029.
    !Config Units = [m/m²]
    !Config         Pan-EU forest-area mean of the EFDA calibration.
    alpha_fire_edge = 0.0237
    CALL getin_p('EDGE_FEEDBACK_ALPHA_FIRE', alpha_fire_edge)
    !
    !Config Key   = EDGE_FEEDBACK_ALPHA_STORM
    !Config Desc  = Edge length generated per m² stormed
    !Config If    = OK_AED_FEEDBACK
    !Config Def   = 0.0241
    !Config Help  = Storm patches are larger and more compact than dispersed
    !Config         bark-beetle mortality → lower perimeter/area ratio.
    !Config         Decoupling: EFDA groups wind+beetle (median α = 0.040 m/m²);
    !Config         Patacca 2023 DFDE 2000-2019 gives f_storm = 63%, f_beetle = 37%
    !Config         in timber volume; geometric hypothesis α_pest = 2·α_storm
    !Config         (beetle patches ~4× smaller area, 2× higher periphery/surface)
    !Config         then yields α_storm ≈ 0.029, α_pest ≈ 0.058, recombining to 0.040.
    !Config         Guillaume M. -- The values in force are LOWER than that derivation
    !Config         (0.0241 and 0.0481, same 2:1 ratio). The derivation above is kept as
    !Config         the method; which of the two is intended has not been re-established.
    !Config Units = [m/m²]
    !Config         Pan-EU forest-area mean of the EFDA calibration.
    alpha_storm_edge = 0.0241
    CALL getin_p('EDGE_FEEDBACK_ALPHA_STORM', alpha_storm_edge)
    !
    !Config Key   = EDGE_FEEDBACK_ALPHA_PEST
    !Config Desc  = Edge length generated per m² pest-killed
    !Config If    = OK_AED_FEEDBACK
    !Config Def   = 0.0481
    !Config Help  = Bark-beetle mortality is dispersed; many small individual
    !Config         patches yield a HIGHER perimeter/area ratio than compact
    !Config         storm or fire scars. See α_storm help for the EFDA+Patacca
    !Config         decoupling procedure.
    !Config Units = [m/m²]
    !Config         Pan-EU forest-area mean of the EFDA calibration.
    alpha_pest_edge = 0.0481
    CALL getin_p('EDGE_FEEDBACK_ALPHA_PEST', alpha_pest_edge)
    !
    !Config Key   = EDGE_FEEDBACK_ALPHA_HARVEST
    !Config Desc  = Edge length generated per m² harvested or LCC-cut
    !Config If    = OK_AED_FEEDBACK
    !Config Def   = 0.0219
    !Config Help  = Commercial cuts (clearcut, thinning, coppice) and LCC-driven
    !Config         clearings are typically rectangular elongated parcels with a
    !Config         high perimeter-to-area ratio. Calibrated as the median dEL/dA
    !Config         observed for the "harvest" agent across AL/IT/FR/DE/SE in EFDA
    !Config         1985-2023. In Europe, harvest dominates total disturbed area
    !Config         (~73 % of pixel-years over the five-country sample) and thus
    !Config         contributes the bulk of the perturbation-driven edge production.
    !Config         Sourced from harvest_area (= veget_max·area·contfrac, m²) summed
    !Config         across PFTs and across iharvest+ilcc components, accumulated in
    !Config         dA_harvest by stomate_lpj.
    !Config Units = [m/m²]
    !Config         Pan-EU forest-area mean of the EFDA calibration.
    alpha_harvest_edge = 0.0219
    CALL getin_p('EDGE_FEEDBACK_ALPHA_HARVEST', alpha_harvest_edge)
    !
    !Config Key   = EDGE_FEEDBACK_TAU_RECOVERY
    !Config Desc  = Edge-length recovery time toward Edge_length.nc reference
    !Config If    = OK_AED_FEEDBACK
    !Config Def   = 3.4715e8 (11 yr in seconds, EFDA 5-country spinup §10bis)
    !Config Help  = Exponential relaxation timescale of edge_length_dyn toward
    !Config         edge_length_ref (Bowring 2024 climatology). Calibrated by
    !Config         AED_SPINUP test §10bis on AL/IT/FR/DE/SE: τ_rec=11 yr gives
    !Config         steady-state endpoint bias <5% vs EFDA EL_obs(2023). The
    !Config         historical default of 20 yr (6.3072e8 s) is deprecated.
    !Config Units = [s]
    tau_rec_edge = 3.4715e8
    CALL getin_p('EDGE_FEEDBACK_TAU_RECOVERY', tau_rec_edge)
    !
    !Config Key   = EDGE_FEEDBACK_REF_FRAC
    !Config Desc  = Fraction of the EFDA edge reference kept as the feedback baseline
    !Config If    = OK_AED_FEEDBACK
    !Config Def   = 0.0 (pure perturbation-driven edge, no EFDA double-count)
    !Config Help  = edge_length_ref (Edge_length.nc) is the OBSERVED EFDA edge, which
    !Config         already integrates the historical perturbations. Using it as the
    !Config         additive baseline / relaxation target of the feedback makes
    !Config         EL_dyn = EFDA + source > EFDA by construction (double-count, see
    !Config         REPORT_AED_SPINUP_G_EFDA.md §4bis). The feedback baseline is now
    !Config         EDGE_FEEDBACK_REF_FRAC·edge_length_ref : 0 = continuous-forest
    !Config         pristine baseline (EL_dyn built purely from simulated disturbances,
    !Config         validated against EFDA) ; 1 = legacy EL_dyn = EFDA + source.
    !Config         edge_length_ref keeps its role as the static edge when the feedback
    !Config         is OFF and as the validation reference.
    !Config Units = [-]
    aed_edge_ref_frac = 0.0
    CALL getin_p('EDGE_FEEDBACK_REF_FRAC', aed_edge_ref_frac)
    !
    !Config Key   = EDGE_MIN_PATCH_AREA
    !Config Desc  = Minimum forest/disturbance patch area used for the edge-length cap
    !Config If    = OK_AED_FEEDBACK
    !Config Def   = 4000. (0.4 ha, EFDA-inferred patch size AL/IT/FR/DE/SE)
    !Config Help  = Upper physical bound on edge_length_dyn is 4·(area·contfrac)/
    !Config         sqrt(EDGE_MIN_PATCH_AREA) = 4·sqrt(area)·sqrt(N_patch), the
    !Config         perimeter of a cell fully tiled by square patches of this size
    !Config         (N_patch = area/EDGE_MIN_PATCH_AREA). Replaces the single-patch
    !Config         cap 4·sqrt(area) which under-bounds fragmented landscapes by
    !Config         2-3 orders of magnitude. Smaller value → higher cap.
    !Config Units = [m²]
    edge_min_patch_area = 4000.
    CALL getin_p('EDGE_MIN_PATCH_AREA', edge_min_patch_area)
    !
    ! ------------------------------------------------------------------------
    ! Guillaume M. -- Management intensity: one concept replacing the six earlier levers
    ! OK_HARVEST_FRAGMENTATION, EDGE_SHAPE_COEF, CLEARCUT_AREA_BY_FMCLASS,
    ! LANDSCAPE_SCENARIO, HARVEST_CLEARCUT_AREA and CLEARCUT_AREA_FILE.
    ! See design/MODULE_DESIGN_MANAGEMENT_INTENSITY.md.
    !Config Key   = OK_MANAGEMENT_INTENSITY
    !Config Desc  = Single management-intensity concept driving rotation, cut diameter and clearcut size
    !Config If    = OK_STOMATE
    !Config Cat   = SCIENTIFIC OPTION
    !Config Def   = n
    !Config Help  = Replaces six earlier and partly overlapping levers by one intensity class per
    !Config         grid cell, read from MANAGEMENT_INTENSITY_FILE. The intensity sets the rotation
    !Config         length, the cut diameter and the clearcut size, hence the edge geometry.
    !Config         See design/MODULE_DESIGN_MANAGEMENT_INTENSITY.md for the parameter table.
    !Config Units = [FLAG]
    ok_management_intensity = .FALSE.
    CALL getin_p('OK_MANAGEMENT_INTENSITY', ok_management_intensity)
    !
    !Config Key   = MANAGEMENT_INTENSITY_FILE
    !Config Desc  = NetCDF map of management-intensity class fractions (fields "f_class1".."f_class5")
    !Config If    = OK_MANAGEMENT_INTENSITY
    !Config Def   = NONE
    !Config Help  = Read annually and interpolated to the model grid (like alpha_map.nc).
    !Config         Fractions of forest area per Scherpenhuijzen 2025 class:
    !Config         1 Unmanaged, 2 Close-to-nature, 3 Combined-objective, 4 Intensive,
    !Config         5 Very-intensive. Built by scripts/build_management_intensity.py
    !Config         (historical 1990-2014 from EFDA-driven decision tree; SSP 2015-2100).
    !Config Units = [FILE]
    management_intensity_file = 'NONE'
    CALL getin_p('MANAGEMENT_INTENSITY_FILE', management_intensity_file)
    !
    !Config Key   = MI_CLEARCUT_SIZE
    !Config Desc  = Clearcut patch size (m2) per management-intensity class 1..5
    !Config If    = OK_MANAGEMENT_INTENSITY
    !Config Def   = 1e9,10000,30000,90000,250000
    !Config Help  = Intensive management leads to LARGE clearcuts, close-to-nature to SMALL
    !Config         ones (landscape-aesthetics constraint). Range 1 ha .. 25 ha; class 1
    !Config         (unmanaged) maps to a huge area -> alpha ~ 0. This is the realisation
    !Config         geometry feeding edge = alpha*cleared_area with alpha = k/sqrt(S), NOT
    !Config         a cut area.
    !Config Units = [m2]
    mi_clearcut_size = (/ 1.e9, 10000., 30000., 90000., 250000. /)
    CALL getin_p('MI_CLEARCUT_SIZE', mi_clearcut_size)
    !
    !Config Key   = MI_DIA_FACTOR
    !Config Desc  = Multiplier on the clearcut diameter, per management-intensity class
    !Config If    = OK_MANAGEMENT_INTENSITY
    !Config Def   = 1., 1., 1., 1., 1.
    !Config Help  = The management intensity modulates the TARGET DIAMETER, not the
    !Config         rotation length (decision 2026-07-27). A MULTIPLIER on
    !Config         largest_tree_dia(pft) rather than an absolute diameter, so that the
    !Config         difference between species is preserved (eucalyptus 0.20 m vs oak
    !Config         0.38 m) while management still shifts the target. Because the age
    !Config         class bounds are FRACTIONS of the clearcut diameter, this rescales the
    !Config         whole maturity ladder of the pixel: an intensively managed cell has
    !Config         its four classes compressed towards small diameters, reaches its cut
    !Config         diameter sooner, produces more age-class-1 area and therefore more
    !Config         edge -- the management/edge coupling becomes structural instead of a
    !Config         calibrated coefficient. All values at 1 = bit-identical.
    !Config         [CATEGORIE : OPTION SCIENTIFIQUE -- choix de modelisation, permanent]
    !Config Units = [-]
    mi_dia_factor(:) = un
    CALL getin_p('MI_DIA_FACTOR', mi_dia_factor)

    !Config Key   = MI_ROTATION_FACTOR
    !Config Desc  = Multiplicateur de la rotation de reference, par classe d'intensite
    !Config If    = OK_MANAGEMENT_INTENSITY
    !Config Def   = 20.0, 1.60, 1.00, 0.65, 0.35
    !Config Help  = rotation(pft,classe) = ROTATION_REF(pft) * MI_ROTATION_FACTOR(classe),
    !Config         le taux de recolte en etant l'inverse. STRICTEMENT symetrique de
    !Config         MI_DIA_FACTOR, qui module le diametre de la meme maniere.
    !Config         /!\ ORDRE DES CLASSES : 1 = non gere, 2 = close-to-nature,
    !Config         3 = multifonctionnel, 4 = intensif, 5 = tres intensif. La classe 3
    !Config         vaut 1.0 par construction : c'est ELLE la reference de ROTATION_REF.
    !Config         Le non gere recoit 20.0 -- soit plusieurs siecles, donc pas de
    !Config         recolte -- plutot que zero, qui rendrait le taux infini.
    !Config         Verification croisee sur les itineraires publies : chene x
    !Config         multifonctionnel = 110 ans (CNPF 100-150), hetre x intensif = 71
    !Config         (CNPF 70-90), eucalyptus x tres intensif = 10.5 (PROF 10-14),
    !Config         eucalyptus x multifonctionnel = 30 (PROF futaie 25-35).
    !Config         Sans effet si ROTATION_REF <= 0 : la rotation reste indefinie, et
    !Config         seuls les PFT arbres lisent cette carte.
    !Config         [CATEGORIE : OPTION SCIENTIFIQUE -- choix de modelisation, permanent]
    !Config Units = [-]
    !Config Key   = N_CLASS0
    !Config Desc  = Duree de la phase de regeneration apres coupe rase (an)
    !Config If    = OK_AGE_CLASS_BOUND_PFT
    !Config Def   = 0
    !Config Help  = Entre la coupe rase et un peuplement installe il s'ecoule 3 a 5 ans :
    !Config         preparation du terrain, plantation, degagements. Le modele TELEPORTE --
    !Config         la coupe produit immediatement un peuplement a qmd_init = 2.5 cm, ce qui
    !Config         correspond deja a un plant de pepiniere de 3 ans (fiche Douglas CNPF :
    !Config         "racines nues de 3 ans maximum, 1+1, 1+2 ou 2+1"). Le convoyeur de classe 0
    !Config         represente ces annees au lieu de les compenser par AGE_STAND_ESTAB.
    !Config         /!\ VALEUR NEUTRE = LE DEFAUT : 0 => convoyeur inactif, bit-neutre.
    !Config         Guillaume M. -- /!\ RESIDENCE TIME IS N_CLASS0 - 1 YEARS: the block runs
    !Config         restitute -> feed -> age, so area fed in year Y lands in slot 2 and is
    !Config         only read back from slot N_CLASS0 in year Y+N_CLASS0-1. Set 5 for a
    !Config         4-year regeneration phase.
    !Config Units = [an]
    n_class0 = 0
    CALL getin_p('N_CLASS0', n_class0)
    IF (n_class0 < 0) THEN
       CALL ipslerr_p(3,'constantes','N_CLASS0 < 0','une duree ne peut pas etre negative','')
    ENDIF

    !Config Key   = AC12_H_RATIO_START
    !Config Desc  = Rapport de hauteur h1/h_ref ou commence la montee de la classe 1 vers la 2
    !Config If    = OK_MATURITY_TRANSFER
    !Config Def   = 0.85
    !Config Help  = Ce qui cree la lisiere est le DIFFERENTIEL de hauteur avec le peuplement
    !Config         voisin, pas une hauteur absolue : un peuplement de 6 m a cote d'un de 7 m
    !Config         n'est pas une lisiere, a cote d'un de 25 m si. h_ref est la hauteur de la
    !Config         PREMIERE classe non vide parmi 2, 3, 4 ; si les trois sont vides on retombe
    !Config         sur la borne absolue AGE_CLASS_FRAC(1)*LARGEST_TREE_DIA.
    !Config         /!\ Calcule en DIAMETRE et jamais en hauteur : H = pipe_tune2*D^pipe_tune3
    !Config         donc r = (D1/D_ref)^pipe_tune3 et pipe_tune2 SE SIMPLIFIE. C'est decisif,
    !Config         pipe_tune2 etant le parametre defaillant (8.6-9.5 m pour un semis de 2.5 cm).
    !Config         0.85 correspond a D1/D_ref = 0.70 a 0.73 selon le MTC : quasi invariant.
    !Config         /!\ NE PAS ABAISSER sans reverifier la retroaction : l'aire qui monte
    !Config         abaisse h_ref, donc augmente r. A 0.85 l'effet est negligeable (h_ref chute
    !Config         de 0.13 % la 1re annee, la fusion etant ponderee par la BIOMASSE et la
    !Config         classe 1 ne pesant que 3.3 % de la classe 2 a surface egale). Plus bas,
    !Config         la boucle se resserre.
    !Config Units = [-]
    ac12_h_ratio_start = 0.85_r_std
    CALL getin_p('AC12_H_RATIO_START', ac12_h_ratio_start)

    !Config Key   = OK_DISTURB_VIA_CLASS0
    !Config Desc  = Feu, tempete et scolyte transitent par le registre de classe 0
    !Config If    = N_CLASS0 > 0
    !Config Def   = n
    !Config Help  = Constat 2026-08-07 : seule la recolte progressive alimentait le registre,
    !Config         donc les perturbations naturelles sautaient la phase de regeneration.
    !Config         Les trois agents sont credites ENSEMBLE, sur le taux somme, comme la
    !Config         condition de bascule les traite deja. Defaut n => bit-neutre.
    !Config         [CATEGORIE : OPTION SCIENTIFIQUE -- choix de modelisation, permanent]
    !Config Units = [-]
    ok_disturb_via_class0 = .FALSE.
    CALL getin_p('OK_DISTURB_VIA_CLASS0', ok_disturb_via_class0)

    !Config Key   = OK_CLASS0_REPLANT
    !Config Desc  = La restitution du convoyeur de classe 0 replante la fraction rendue
    !Config If    = N_CLASS0 > 0
    !Config Def   = n
    !Config Help  = Sans cela l'aire rase rendue au creneau de classe 1 n'est JAMAIS
    !Config         regarnie : `establish` exige un creneau entierement vide de biomasse et
    !Config         `recruite` est interdit sur 96,8 % de l'aire (ifm_thin_clearcut). La
    !Config         densite de classe 1 s'effondre (2700 -> 50 tiges/ha) pendant que son QMD
    !Config         monte a 88 cm : quelques arbres restes seuls dans le creneau.
    !Config         Defaut n => bit-neutre.
    !Config         [CATEGORIE : OPTION SCIENTIFIQUE -- choix de modelisation, permanent]
    !Config Units = [-]
    ok_class0_replant = .FALSE.
    CALL getin_p('OK_CLASS0_REPLANT', ok_class0_replant)

    !Config Key   = AC12_DIA_RATIO_START
    !Config Desc  = Rapport max_dia/borne a partir duquel les classes 2 et 3 montent
    !Config If    = OK_AGE_CLASS_BOUND_PFT
    !Config Def   = 0.60
    !Config Help  = Au-dela de la classe 1 la bascule n'est plus une fermeture de couvert
    !Config         mais une dimension d'exploitabilite : elle se juge en DIAMETRE, sans
    !Config         l'exposant allometrique pipe_tune3. Consequence a connaitre : ce dernier
    !Config         vaut 0.42 a 0.52 selon le MTC, donc le critere en hauteur etait
    !Config         implicitement plus severe pour certaines essences que pour d'autres
    !Config         (0.85 en hauteur = 0.66 a 0.72 en diametre). Le critere en diametre est
    !Config         uniforme. 0.60 est DELIBEREMENT sous l'equivalent iso-comportement
    !Config         (~0.70) : les classes 2 et 3 se plaquent contre leur borne, mediane
    !Config         mesuree 0.67 en diametre, et un seuil a 0.70 les y laisserait.
    !Config Units = [-]
    ac12_dia_ratio_start = 0.60_r_std
    CALL getin_p('AC12_DIA_RATIO_START', ac12_dia_ratio_start)

    !Config Key   = AC12_DIA_RATIO_LAST
    !Config Desc  = Seuil propre au dernier passage, celui qui alimente la classe terminale
    !Config If    = OK_AGE_CLASS_BOUND_PFT
    !Config Def   = -1 (repli sur AC12_DIA_RATIO_START, donc bit-neutre)
    !Config Help  = La classe terminale est la seule dont le receveur est aussi l'ensemble
    !Config         coupable : elle est remplie au 31 decembre puis rasee dans la meme
    !Config         sequence. Son bilan est negatif sur 10 pixels sur 12, la sortie valant
    !Config         1.2 a 4.9 fois l'alimentation, et les deux seuls pixels excedentaires
    !Config         sont ceux qui alimentent 13 et 15 annees sur 21 contre 3 a 11 ailleurs.
    !Config         Le discriminant est donc la CADENCE d'alimentation. Un seuil distinct
    !Config         permet de l'accelerer sans toucher au passage 2->3, dont la cadence est
    !Config         deja de 19 a 21 annees sur 21 sur plusieurs pixels.
    !Config Units = [-]
    ! Guillaume M. -- Read AFTER ac12_dia_ratio_start: the fallback copies it, so reading
    ! it first would capture the default instead of the namelist value.
    ac12_dia_ratio_last = -un
    CALL getin_p('AC12_DIA_RATIO_LAST', ac12_dia_ratio_last)
    IF (ac12_dia_ratio_last <= zero) ac12_dia_ratio_last = ac12_dia_ratio_start

    !Config Key   = OK_RDI_RAMP
    !Config Desc  = La consigne de RDI descend pendant la classe d'amelioration
    !Config If    = OK_CLEARCUT_LAST_CLASS
    !Config Def   = n
    !Config Help  = Sans rampe, la bande RDI est CONSTANTE sur toute la classe : le
    !Config         peuplement subirait un echelon a l'entree au lieu d'une eclaircie
    !Config         progressive. La rampe interpole entre RDI_BAND_LOW (entree) et
    !Config         RDI_RAMP_OUT (sortie), en fonction de l'avancement du peuplement mesure
    !Config         par son DIAMETRE entre les deux bornes de sa classe.
    !Config         Interpolation GEOMETRIQUE : on interpole un rapport, et le nombre de
    !Config         tiges decroit multiplicativement (taux de retrait constant), pas par
    !Config         soustraction constante.
    !Config         Le diametre plutot que l'age : c'est deja lui qui definit l'appartenance
    !Config         a la classe, et il porte la fertilite de la maille -- une station riche
    !Config         progresse plus vite dans la rampe, ce qui est l'effet recherche.
    !Config         [CATEGORIE : OPTION SCIENTIFIQUE -- choix de modelisation, permanent]
    !Config Units = [FLAG]
    ok_rdi_ramp = .FALSE.
    CALL getin_p('OK_RDI_RAMP', ok_rdi_ramp)

    !Config Key   = OK_FORCED_THINNING
    !Config Desc  = Eclaircie ANNUELLE forcee sur les classes de maturite jeunes
    !Config If    = OK_STOMATE
    !Config Def   = n
    !Config Help  = Le declencheur d'eclaircie ordinaire est un franchissement de bande RDI.
    !Config         Il ne mord pas sur un peuplement dont la droite d'auto-eclaircie tolere
    !Config         beaucoup de tiges : le RDI reste sous la bande et l'eclaircie n'a jamais
    !Config         lieu. Ce drapeau rend l'eclaircie ANNUELLE et INCONDITIONNELLE sur les
    !Config         classes jusqu'a FORCED_THIN_LAST_CLASS, a intensite FORCED_THIN_INTENSITY.
    !Config         Sur ces classes elle REMPLACE le declencheur RDI, elle ne s'y ajoute pas.
    !Config         Voir design/MODULE_DESIGN_FORCED_THINNING.md
    !Config         [CATEGORIE : OPTION SCIENTIFIQUE -- choix de modelisation, permanent]
    !Config Units = [FLAG]
    ok_forced_thinning = .FALSE.
    CALL getin_p('OK_FORCED_THINNING', ok_forced_thinning)

    !Config Key   = FORCED_THIN_LAST_CLASS
    !Config Desc  = Derniere classe de maturite soumise a l'eclaircie forcee
    !Config If    = OK_FORCED_THINNING
    !Config Def   = 2
    !Config Help  = Les classes 1 a FORCED_THIN_LAST_CLASS sont eclaircies chaque annee.
    !Config         Au-dela, le declencheur RDI ordinaire reprend la main. La classe
    !Config         terminale reste hors eclaircie par construction (thin_this_slot).
    !Config         [CATEGORIE : OPTION SCIENTIFIQUE -- choix de modelisation, permanent]
    !Config Units = [-]
    forced_thin_last_class = 2
    CALL getin_p('FORCED_THIN_LAST_CLASS', forced_thin_last_class)

    !Config Key   = OK_MIN_DENSITY_RDI
    !Config Desc  = CONDITION 1 juge l'effondrement sur le RDI, pas sur un nombre de tiges
    !Config If    = OK_MIN_DENSITY_RESET
    !Config Def   = n
    !Config Help  = CONDITION 1 rase un peuplement dont la densite tombe sous DENS_TARGET
    !Config         (200 tiges/ha). Son commentaire reconnait la faiblesse : "density alone
    !Config         cannot tell a stand that is open BY DESIGN from one that has collapsed",
    !Config         et le discriminant retenu est la TAILLE. Or un peuplement d'arbres
    !Config         d'avenir, c'est peu de tiges DEJA EPAISSES -- exactement la signature
    !Config         attribuee a l'effondrement. Le discriminant en taille ne separe donc
    !Config         plus rien des qu'on ouvre volontairement un peuplement murissant.
    !Config         A y, l'effondrement se juge sur le RDI : un peuplement conduit a sa
    !Config         consigne n'est pas effondre, quel que soit son nombre de tiges.
    !Config         /!\ Le RDI lu est celui de l'ANNEE PRECEDENTE (rdi_stand) : rdi n'est
    !Config         calcule qu'apres CONDITION 1 dans sapiens_forestry. Un decalage d'un an
    !Config         sur un detecteur d'effondrement est sans consequence.
    !Config         [CATEGORIE : OPTION SCIENTIFIQUE -- choix de modelisation, permanent]
    !Config Units = [FLAG]
    ok_min_density_rdi = .FALSE.
    CALL getin_p('OK_MIN_DENSITY_RDI', ok_min_density_rdi)

    !Config Key   = MIN_DENSITY_RDI_FLOOR
    !Config Desc  = RDI sous lequel un peuplement est tenu pour effondre
    !Config If    = OK_MIN_DENSITY_RDI
    !Config Def   = 0.02
    !Config Help  = Doit rester STRICTEMENT SOUS la consigne de densite la plus basse que
    !Config         l'on compte imposer, sinon le garde-fou rase le peuplement qu'on
    !Config         cherche a construire. Repere : une conduite en arbres d'avenir vise
    !Config         ~59 tiges/ha a 16 cm, soit RDI ~ 0.04 ; le plancher doit donc etre
    !Config         nettement en dessous. 0.02 laisse un facteur deux de marge.
    !Config Units = [-]
    min_density_rdi_floor = 0.02_r_std
    CALL getin_p('MIN_DENSITY_RDI_FLOOR', min_density_rdi_floor)

    ! Guillaume M. -- COPPICE_2CLASS: a coppice has no maturity ladder in diameter.
    ! See design/MODULE_DESIGN_COPPICE_2CLASS.md.
    !
    !Config Key   = OK_COPPICE_2CLASS
    !Config Desc  = Deux classes de maturite pour le taillis, graduees en HAUTEUR
    !Config If    = OK_AGE_CLASS_BOUND_PFT
    !Config Def   = n
    !Config Help  = L'echelle de maturite en diametre est un objet de futaie reguliere :
    !Config         elle encode "bois de qualite croissante jusqu'a la dimension
    !Config         marchande". Un taillis produit beaucoup de brins fins et ne franchit
    !Config         jamais le premier barreau -- mesure sur la sclerophylle : 4 a 5 cm de
    !Config         diametre dominant contre une borne de classe 1 a 7.5 cm, soit 100 %
    !Config         de l'aire sous la borne toutes les annees. La branche de declassement
    !Config         renvoie donc en permanence chaque creneau mur vers la classe 1, aire
    !Config         et biomasse conservees, et la structure d'age s'effondre en un an.
    !Config         A y, les creneaux en taillis n'ont plus que la classe 1 et la classe
    !Config         terminale, et le passage entre les deux se juge en hauteur.
    !Config         [CATEGORIE : OPTION SCIENTIFIQUE -- choix de modelisation, permanent]
    !Config Units = [FLAG]
    ok_coppice_2class = .FALSE.
    CALL getin_p('OK_COPPICE_2CLASS', ok_coppice_2class)

    !Config Key   = COP_R_UP
    !Config Desc  = Rapport de hauteurs a partir duquel la classe 1 alimente la terminale
    !Config If    = OK_COPPICE_2CLASS
    !Config Def   = 1.0
    !Config Help  = r = (max_dia/borne de classe 1)^pipe_tune3 est un rapport de HAUTEURS :
    !Config         H = pipe_tune2 * D^pipe_tune3, donc pipe_tune2 se simplifie et le
    !Config         critere ne depend pas de ce parametre defaillant. r = 1 signifie que
    !Config         les plus gros brins ont atteint la taille qui definit la fermeture du
    !Config         couvert -- le seul barreau qui garde un sens en taillis.
    !Config Units = [-]
    cop_r_up = un
    CALL getin_p('COP_R_UP', cop_r_up)

    !Config Key   = COP_R_DOWN
    !Config Desc  = Rapport de hauteurs sous lequel un creneau en taillis retombe en classe 1
    !Config If    = OK_COPPICE_2CLASS
    !Config Def   = 0.4
    !Config Help  = HYSTERESIS OBLIGATOIRE, et c'est la raison d'etre de ce parametre. Un
    !Config         seuil de descente egal a celui de montee fait osciller le peuplement
    !Config         contre sa borne au lieu de le faire passer. Le seuil doit aussi rester
    !Config         STRICTEMENT POSITIF : la descente est le seul canal par lequel une
    !Config         tempete ou un scolyte renvoie un peuplement rase en classe 1, la
    !Config         supprimer casserait le rajeunissement par perturbation.
    !Config Units = [-]
    cop_r_down = 0.4_r_std
    CALL getin_p('COP_R_DOWN', cop_r_down)

    ! Guillaume M. -- Explicit guard rather than a silent clamp: inverted thresholds do not
    ! degrade the output, they make a stand climb and fall back through the same slot every
    ! year, which reads as a model defect instead of a configuration error.
    IF (ok_coppice_2class .AND. cop_r_down >= cop_r_up) THEN
       WRITE(numout,*) 'COP_R_DOWN, COP_R_UP: ', cop_r_down, cop_r_up
       CALL ipslerr_p(3, 'constantes', 'COP_R_DOWN must be strictly below COP_R_UP', &
            'equal or inverted thresholds make the stand oscillate against its bound', &
            'lower COP_R_DOWN or raise COP_R_UP')
    ENDIF

    !Config Key   = OK_MATURITY_TRANSFER
    !Config Desc  = Promotion de maturite par COPIE INTACTE du peuplement, cadence MAT_F_REF
    !Config If    = OK_AGE_CLASS_BOUND_PFT
    !Config Def   = n
    !Config Help  = When y, three things change at once: (a) the moving stand is an intact
    !Config         copy of the donor -- every circumference class travels, none is selected,
    !Config         both sides keep their per-m2 and per-tree properties unchanged; (b) the
    !Config         area fraction moved per year is MAT_F_REF(iagec), a prescribed rate, not
    !Config         a function of the size distribution; (c) merge_biomass_pfts merges the
    !Config         receiver and the entrant by diameter proximity instead of re-drawing the
    !Config         Weibull, so the receiver structure is history, not a law.
    !Config         The eligibility thresholds (AC12_*_RATIO_*) still decide IF a stand may
    !Config         be promoted; this flag decides HOW MUCH and WHAT.
    !Config         /!\ NEUTRAL VALUE = THE DEFAULT (n) => historical paths, bit-neutral.
    !Config Units = [-]
    ok_maturity_transfer = .FALSE.
    CALL getin_p('OK_MATURITY_TRANSFER', ok_maturity_transfer)
    ! Guillaume M. -- Hard dependency of COPPICE_2CLASS, checked here because
    ! OK_MATURITY_TRANSFER is read after it. Coppice promotion rides the partial-transfer
    ! branch: without it no slot has any path to its terminal class.
    IF (ok_coppice_2class .AND. .NOT. ok_maturity_transfer) THEN
       CALL ipslerr_p(3, 'constantes', 'OK_COPPICE_2CLASS=y requires OK_MATURITY_TRANSFER=y', &
            'the coppice promotion rides the partial-transfer branch; without it no slot &
            &can ever reach its terminal class', '')
    ENDIF
    !
    !Config Key   = MAT_BAND_LOW
    !Config Desc  = Lower edge of the regulator dead band on the receiver filling ratio
    !Config If    = OK_MATURITY_TRANSFER
    !Config Def   = 0.7
    !Config Help  = The promotion rate stays at MAT_F_REF while the receiver filling ratio
    !Config         A_rec / (MAT_TARGET_FRAC * A_group) sits inside [MAT_BAND_LOW,
    !Config         MAT_BAND_HIGH]. The band exists so the regulator corrects structure,
    !Config         not noise: a target is a range, never a number to hit exactly.
    !Config Units = [-]
    mat_band_low = 0.7_r_std
    CALL getin_p('MAT_BAND_LOW', mat_band_low)
    !
    !Config Key   = MAT_BAND_HIGH
    !Config Desc  = Upper edge of the regulator dead band on the receiver filling ratio
    !Config If    = OK_MATURITY_TRANSFER
    !Config Def   = 1.3
    !Config Units = [-]
    mat_band_high = 1.3_r_std
    CALL getin_p('MAT_BAND_HIGH', mat_band_high)
    IF (mat_band_low < zero .OR. mat_band_high < mat_band_low) THEN
       CALL ipslerr_p(3,'constantes','MAT_BAND_LOW/HIGH inconsistent', &
            'need 0 <= MAT_BAND_LOW <= MAT_BAND_HIGH','')
    ENDIF
    !
    !Config Key   = MAT_GAIN
    !Config Desc  = Exponent of the regulator cross correction outside the dead band
    !Config If    = OK_MATURITY_TRANSFER
    !Config Def   = 1.0
    !Config Help  = Outside the dead band f = MAT_F_REF * (r_don/r_rec)**MAT_GAIN: it
    !Config         pushes when the donor is full AND the receiver empty, holds back in
    !Config         the opposite case. 0 disables the regulator exactly (f = MAT_F_REF
    !Config         everywhere), which is the A/B lever against the unregulated transfer.
    !Config         The growth time constant is decadal: a high gain will oscillate.
    !Config Units = [-]
    mat_gain = un
    CALL getin_p('MAT_GAIN', mat_gain)
    IF (mat_gain < zero) THEN
       CALL ipslerr_p(3,'constantes','MAT_GAIN < 0', &
            'a negative gain inverts the correction; use 0 to disable the regulator','')
    ENDIF
    !
    !Config Key   = MAT_F_MIN
    !Config Desc  = Floor of the regulated promotion rate
    !Config If    = OK_MATURITY_TRANSFER
    !Config Def   = 0.0
    !Config Units = [1/year]
    mat_f_min = zero
    CALL getin_p('MAT_F_MIN', mat_f_min)
    !
    !Config Key   = MAT_F_MAX
    !Config Desc  = Ceiling of the regulated promotion rate
    !Config If    = OK_MATURITY_TRANSFER
    !Config Def   = 0.15
    !Config Help  = Hard bound on what the regulator may demand in one year. An empty
    !Config         receiver drives the raw correction to infinity; this ceiling is what
    !Config         actually limits it, so it must stay well under one to keep the donor
    !Config         from emptying in a few years.
    !Config Units = [1/year]
    mat_f_max = 0.15_r_std
    CALL getin_p('MAT_F_MAX', mat_f_max)
    IF (mat_f_min < zero .OR. mat_f_max > un .OR. mat_f_max < mat_f_min) THEN
       CALL ipslerr_p(3,'constantes','MAT_F_MIN/MAX inconsistent', &
            'need 0 <= MAT_F_MIN <= MAT_F_MAX <= 1','')
    ENDIF
    !
    !Config Key   = MAT_SHIFT
    !Config Desc  = Maturity step of the moving stand (relative diameter offset)
    !Config If    = OK_MATURITY_TRANSFER
    !Config Def   = 0.05
    !Config Help  = The moving stand is displaced FORWARD along the self-thinning
    !Config         trajectory (diameter * (1+s), stems * (1+s)**(1/beta)) and the donor
    !Config         backward, so donor and receiver never touch: this is what creates the
    !Config         gap between adjacent maturity classes. Measured over 100 years, the
    !Config         area regulator alone lets classes 2-4 converge within 6 cm of QMD and
    !Config         the mature class shrink to 12 percent: the gap is a hard requirement
    !Config         of the regulator, not an option. 0 restores the intact copy (A/B).
    !Config Units = [-]
    mat_shift = 0.05_r_std
    CALL getin_p('MAT_SHIFT', mat_shift)
    IF (mat_shift < zero .OR. mat_shift > 0.5_r_std) THEN
       CALL ipslerr_p(3,'constantes','MAT_SHIFT out of [0, 0.5]', &
            'a maturity step is a small relative diameter offset','')
    ENDIF
    ! Guillaume M. -- The regulator without the gap is the measured degeneracy mode
    ! (classes converge, then the maturity recall vanishes for good). Legitimate only
    ! for an A/B isolation run, never for production.
    IF (ok_maturity_transfer .AND. mat_gain > zero .AND. .NOT. (mat_shift > zero)) THEN
       CALL ipslerr_p(2,'constantes','area regulator active with MAT_SHIFT = 0', &
            'adjacent maturity classes will converge and the maturity scale degenerate', &
            'acceptable for an A/B test only')
    ENDIF
    !
    !Config Key   = OK_HARVEST_FEED_COUPLING
    !Config Desc  = Per-species harvest quotas, demand modulated by the mature-class filling ratio
    !Config If    = OK_PROGRESSIVE_HARVEST
    !Config Def   = n
    !Config Help  = The progressive-harvest quota is an ABSOLUTE area flux A_forest/R taken
    !Config         from the mature class whatever its state: the smaller the class, the
    !Config         larger the fraction cut (an anti-regulation; the feed side is regulated,
    !Config         the cut side is open loop -- design doc 4.4). When y, two things change:
    !Config         (a) the quota is computed PER age-class species group, demand and supply
    !Config         both filtered on thinning-managed slots (one allowable cut per working
    !Config         circle -- the pixel-level pot pays a deficit of mature beech with mature
    !Config         spruce); (b) the demand is multiplied by g = r4**HARV_COUPLING_GAIN
    !Config         outside the MAT_BAND dead band, r4 being the mature-class filling ratio
    !Config         against MAT_TARGET_FRAC(nagec): deficit defers the cut (saving cuttings),
    !Config         surplus anticipates it (recovery), inside the band the rotation rules
    !Config         alone. g is clamped to the named constants harv_g_min/max = 0.5/1.5.
    !Config         /!\ NEUTRAL VALUE = THE DEFAULT (n) => historical pixel quota, bit-neutral.
    !Config         [CATEGORIE : OPTION SCIENTIFIQUE -- choix de modelisation, permanent]
    !Config Units = [-]
    ok_harvest_feed_coupling = .FALSE.
    CALL getin_p('OK_HARVEST_FEED_COUPLING', ok_harvest_feed_coupling)
    !
    !Config Key   = HARV_COUPLING_GAIN
    !Config Desc  = Exponent of the harvest demand modulation outside the dead band
    !Config If    = OK_HARVEST_FEED_COUPLING
    !Config Def   = 1.0
    !Config Help  = g = r4**gain, clamped to [harv_g_min, harv_g_max]. 0 gives g = 1
    !Config         everywhere: per-group quotas WITHOUT modulation, the isolation lever
    !Config         of the granularity change (which is not bit-neutral by construction).
    !Config Units = [-]
    harv_coupling_gain = un
    CALL getin_p('HARV_COUPLING_GAIN', harv_coupling_gain)
    IF (harv_coupling_gain < zero) THEN
       CALL ipslerr_p(3,'constantes','HARV_COUPLING_GAIN < 0', &
            'a negative gain inverts the correction; use 0 for per-group quotas only','')
    ENDIF

    mi_rotation_factor(:) = (/ 20.0_r_std, 1.60_r_std, 1.00_r_std, 0.65_r_std, 0.35_r_std /)
    CALL getin_p('MI_ROTATION_FACTOR', mi_rotation_factor)
    IF (ANY(mi_rotation_factor(:) <= zero)) THEN
       WRITE(numout,*) 'MI_ROTATION_FACTOR doit etre strictement positif : ', mi_rotation_factor(:)
       CALL ipslerr_p(3,'constantes.f90','MI_ROTATION_FACTOR must be > 0', &
            'une rotation nulle ou negative rendrait le taux infini ou negatif','')
    ENDIF
    IF (ANY(mi_dia_factor(:) <= zero)) THEN
       WRITE(numout,*) 'MI_DIA_FACTOR must be strictly positive : ', mi_dia_factor(:)
       CALL ipslerr_p(3,'constantes.f90','MI_DIA_FACTOR must be > 0',&
            'it multiplies the clearcut diameter and the age class bounds','')
    ENDIF
    !
    !Config Key   = MI_ROTATION_AGE
    !Config Desc  = Target rotation age (yr) per management-intensity class 1..5
    !Config If    = OK_MANAGEMENT_INTENSITY
    !Config Def   = 999,180,120,90,20
    !Config Help  = Class 1 (unmanaged) = 999 yr, i.e. never reaches the rotation criterion
    !Config         (de facto protection). This is the ONLY channel through which management
    !Config         intensity acts on the model dynamics (via AGE_ROTATION).
    !Config         CALIBRATION (2026-07-26, design §2.6) — literature, NOT the EFDA
    !Config         renewal-rate estimator, which is biased HIGH (Landsat misses partial
    !Config         fellings; the forest mask includes unproductive/set-aside area):
    !Config           class 5 = 20 yr  : Duncker et al. 2012 define short-rotation forestry
    !Config                              as <=20 yr; Scherpenhuijzen calibrated the
    !Config                              disturbance-frequency threshold on that definition.
    !Config                              Eucalyptus S. Europe 8-12 yr.
    !Config           class 4 = 90 yr  : S. Sweden pine/spruce ~80 yr, Picea N. Europe ~100 yr,
    !Config                              Central Europe spruce 70-110 yr.
    !Config           classes 2-3 (180/120) : expert judgement — the literature reports
    !Config                              rotations by species and country, not by FMA, and
    !Config                              continuous-cover forestry has no true rotation.
    !Config         NB the edge/rotation first-order compensation is ROBUST to this table
    !Config         (product within +-3 % across very different tables) : what it changes is
    !Config         the HARVEST REGIME (pan-EU mean rotation 94 -> 124 yr), not the edge.
    !Config Units = [yr]
    !Config /!\ MI_ROTATION_AGE supprimee : la rotation vient de ROTATION_REF par PFT,
    !Config     mise a l'echelle par MI_ROTATION_FACTOR selon la classe d'intensite.
    !
    ! Guillaume M. -- MI_HARVEST_RATE removed: the rotation of every tree PFT comes
    ! from ROTATION_REF scaled by MI_ROTATION_FACTOR, which overrode the rate-based
    ! value everywhere it was read. The rate only ever filled entries no tree reads.
    !
    !Config Key   = MI_EDGE_SHAPE_COEF
    !Config Desc  = Geometric shape factor k in alpha_c = k/sqrt(MI_CLEARCUT_SIZE(c))
    !Config If    = OK_MANAGEMENT_INTENSITY
    !Config Def   = 3.800
    !Config Help  = Calibrated under the aggregate constraint: the pan-EU forest-area
    !Config         weighted mean of alpha must equal 0.0219, the value reproducing the
    !Config         EFDA edge observation (12.73 Mm). Re-derive on the REALISED alpha when
    !Config         MI_AGE_MOD_STRENGTH /= 0 (S then depends on model state).
    !Config Units = [-]
    mi_edge_shape_coef = 3.800
    CALL getin_p('MI_EDGE_SHAPE_COEF', mi_edge_shape_coef)
    !
    !Config Key   = MI_AGE_MOD_STRENGTH
    !Config Desc  = Strength of the clearcut-size modulation by the stand maturity index J
    !Config If    = OK_MANAGEMENT_INTENSITY
    !Config Def   = 1.0
    !Config Help  = Where mature forest forms a large block, real clearcuts are large;
    !Config         where it is scattered or young they are small. Carrier J = mean
    !Config         normalised maturity over the 4 age classes (0 = all young, 1 = all
    !Config         mature). NOT a concentration index: a Herfindahl is blind to direction
    !Config         (all-mature and all-just-cut both give HHI=1) and measured flat over 52
    !Config         model years, whereas J shows the ageing/harvest phase alternation.
    !Config         Makes S a PHASE DIAGNOSTIC: S grows while stands age, shrinks during
    !Config         intensive harvest. 0 = modulation off (per-class table only).
    !Config Units = [-]
    mi_age_mod_strength = 1.0
    CALL getin_p('MI_AGE_MOD_STRENGTH', mi_age_mod_strength)
    !
    !Config Key   = MI_CLEARCUT_SIZE_MIN
    !Config Desc  = Lower bound of the clearcut patch size after age modulation
    !Config If    = OK_MANAGEMENT_INTENSITY
    !Config Def   = 10000. (1 ha)
    !Config Units = [m2]
    mi_clearcut_size_min = 10000.
    CALL getin_p('MI_CLEARCUT_SIZE_MIN', mi_clearcut_size_min)
    !
    !Config Key   = MI_CLEARCUT_SIZE_MAX
    !Config Desc  = Upper bound of the clearcut patch size after age modulation
    !Config If    = OK_MANAGEMENT_INTENSITY
    !Config Def   = 250000. (25 ha)
    !Config Units = [m2]
    mi_clearcut_size_max = 250000.
    CALL getin_p('MI_CLEARCUT_SIZE_MAX', mi_clearcut_size_max)
    !
    !Config Key   = PEST_BIOMASS_REF
    !Config Desc  = Standing biomass density used to convert beetle kill (gC/m²) to area
    !Config If    = OK_AED_FEEDBACK
    !Config Def   = 10000.0 (temperate mature stand)
    !Config Help  = Used in stomate_pest to map B_beetles_kill (gC m-2 yr-1) to
    !Config         a m² area affected per m² of forest, before applying alpha_pest_edge.
    !Config Units = [gC/m²]
    pest_biomass_ref = 10000.0
    CALL getin_p('PEST_BIOMASS_REF', pest_biomass_ref)
    !Config Key   = AED_BEETLE_EPIDEMIC_THRESHOLD
    !Config Desc  = Mass-attack index (i_beetles_massattack, 0-1) above which beetle
    !Config         kills feed the AED edge budget (epidemic phase). Below this the
    !Config         kills are endemic = background mortality of isolated individuals,
    !Config         which does not open a forest/non-forest edge -> excluded from dA_pest.
    !Config Units = [-]
    aed_beetle_epidemic_threshold = 0.5
    CALL getin_p('AED_BEETLE_EPIDEMIC_THRESHOLD', aed_beetle_epidemic_threshold)
    !
    !Config Key   = OK_DIA_STAGGER
    !Config Desc  = Cold start: establish each age class at its own diameter bin
    !Config If    = OK_STOMATE
    !Config Cat   = SCIENTIFIC OPTION
    !Config Def   = n
    !Config Help  = A cold start otherwise gives every age class the same structure, so all
    !Config         classes reach the clearcut threshold together and produce a synchronised
    !Config         harvest and edge wave. This flag anchors each class on its age_class_bound
    !Config         diameter and splits heartwood/sapwood from that diameter instead of the
    !Config         legacy 0.01 floor. It fixes the TIMING of the age-class machinery, not the
    !Config         stocking level, and NOT the area split (see AGE_CLASS_INIT_FRAC).
    !Config Units = [FLAG]
    ok_dia_stagger = .FALSE.
    CALL getin_p('OK_DIA_STAGGER', ok_dia_stagger)
    !
    ! Guillaume M. -- OK_EDGE_3TERMS removed: the meso term is now UNCONDITIONAL. It was a
    ! validation scaffold whose exit criterion is met. Keeping it would have left the
    ! configurations that switched it off with no mature-forest edge at all, once
    ! AED_OPEN_W1..W4 were removed.
    !
    !Config Key   = ROAD_LENGTH_FILE
    !Config Desc  = Road-length map for the MACRO edge term; 'NONE' switches the term off
    !Config If    = OK_AED_FEEDBACK
    !Config Def   = NONE
    !Config Help  = Restores the ORIGINAL AED channel: Bowring, Li, Mouillot, Rosan & Ciais
    !Config         (2024, Nat. Commun. 15, 9176) drive fragmentation in
    !Config         ORCHIDEE-MICT-SPITFIRE from GRIP road density; this project had replaced
    !Config         it by the EFDA forest/non-forest edge. The two are ADDITIVE: on 5249
    !Config         internal gaps of the EFDA masks only 7 have a road-like shape (<0.1 % of
    !Config         the area) -- at 30 m a forest road hides under a canopy that closes over
    !Config         it. Built by scripts/build_road_edge_panEU.py. The file holds a ROAD
    !Config         length (m per cell), NOT an edge length: the model derives
    !Config         EL_macro = 2 * RL * ROAD_GRIP_CORRECTION. 'NONE' = bit-neutral.
    !Config         See design/MODULE_DESIGN_FRAGMENTATION_3TERMES.md section 3.3.
    !Config Units = [FILE]
    road_length_file = 'NONE'
    CALL getin_p('ROAD_LENGTH_FILE', road_length_file)
    !
    !Config Key   = ROAD_GRIP_CORRECTION
    !Config Desc  = GRIP under-reporting correction applied to the road length
    !Config If    = ROAD_LENGTH_FILE
    !Config Def   = 2.
    !Config Help  = GRIP under-represents actual road length, chiefly unpaved roads. Bowring
    !Config         et al. (2024) assume a factor 2, which puts total global road length
    !Config         between the World Road Statistics and CIA estimates. Set to 1 for raw
    !Config         GRIP.
    !Config         /!\ This is NOT the factor 2 printed in their AED equation
    !Config         AED = 2*Area*f/sum(RL). THAT one is GEOMETRIC -- it falls out of giving
    !Config         each circular patch its FULL perimeter without halving shared boundaries
    !Config         (n = A/(pi r^2) patches, 2*pi*r each => RL = 2A/r). There are two
    !Config         factors 2 in that paper and they are not in the same place. Ours is
    !Config         hard-coded in update_edge_length (a road has two sides) and is not a knob.
    !Config Units = [-]
    road_grip_correction = 2._r_std
    CALL getin_p('ROAD_GRIP_CORRECTION', road_grip_correction)
    !
    !Config Key   = OK_ROAD_FOREST_MASK
    !Config Desc  = The road map is already restricted to forest (variable road_length_forest)
    !Config If    = ROAD_LENGTH_FILE
    !Config Def   = y
    !Config Help  = TRUE reads 'road_length_forest': road length already weighted by the
    !Config         forest fraction computed at 5 arcmin (~8 km) from the EFDA 30 m masks.
    !Config         Much finer than the model's 1-degree fraction, hence far less biased --
    !Config         but STATIC, so it does not follow land-use change.
    !Config         FALSE reads 'road_length_total' and applies the model's own forest
    !Config         fraction instead: it then follows LUC for free, but OVERESTIMATES forest
    !Config         road access, since roads are denser outside forest. The gap between the
    !Config         two is printed by scripts/build_road_edge_panEU.py -- measure it, do not
    !Config         assume it.
    !Config         [CATEGORIE : OPTION SCIENTIFIQUE -- masque statique fin vs suivi du
    !Config         changement d'usage : arbitrage de modelisation, permanent]
    !Config         Guillaume M. -- The y default is deliberate: every configuration runs
    !Config         on it without setting it. Under a moving land cover, n becomes the
    !Config         defensible choice, at the price of a high bias.
    !Config Units = [FLAG]
    ok_road_forest_mask = .TRUE.
    CALL getin_p('OK_ROAD_FOREST_MASK', ok_road_forest_mask)
    !
    !Config Key   = EDGE_GRAIN_DEFAULT
    !Config Desc  = Fallback sub-grid fragmentation grain for the meso edge term (m)
    !Config If    = OK_AED_FEEDBACK
    !Config Def   = 384.
    !Config Help  = Pan-EU median of the per-cell inversion against Edge_length_1deg.nc.
    !Config         Used for MTCs with no calibrated grain -- tropical MTC 2 and 3, absent
    !Config         from the EFDA map. Per-MTC values live in edge_grain_mtc_of
    !Config         (stomate_data). Must stay in 1e2-1e3 m: outside that range it is a
    !Config         correction factor, not a length.
    !Config Units = [m]
    edge_grain_default = 384.
    CALL getin_p('EDGE_GRAIN_DEFAULT', edge_grain_default)
    !
    !
    ! Guillaume M. -- AGECLASS_FRAC_FILE removed: the cold-start class split is set by
    ! AGE_CLASS_INIT_FRAC only. Initialising from an observed map and then validating
    ! against that same observation is a tautology; an even split is a neutral reference.
    !
    ! Guillaume M. -- /!\ AGE_CLASS_AGE_FRAC MOVED from here to just after the NAGEC read.
    ! At this point nagec still holds its default 1, so the array was sized to 1 and
    ! classes 2..nagec were read OUT OF BOUNDS. Allocate NOTHING here that depends on nagec.
    !
    ! Guillaume M. -- The whole maturity scale is derived from LARGEST_TREE_DIA per PFT;
    ! no external diameter map is involved.
    !
    !
    !
    ! Guillaume M. -- AGE_RECRUIT_OFFSET stays REMOVED. Its role -- offsetting the
    ! seedling-to-sapling lag -- is taken over below by AGE_STAND_ESTAB, on a verifiable
    ! ground: `prescribe` installs seedlings at a non-zero DIAMETER, rather than an
    ! adjustment of an age prescribed at initialisation.
    !
    !Config Key   = AGE_STAND_ESTAB
    !Config Desc  = Age attribue a un peuplement fraichement etabli (plancher de age_stand_bm)
    !Config If    = OK_STOMATE
    !Config Def   = 5.
    !Config Help  = `prescribe` n'installe pas des graines mais des semis d'un DIAMETRE
    !Config         prescrit non nul : le peuplement nait donc deja age de quelques annees.
    !Config         Sans ce plancher, une classe 1 regeneree apres coupe rase afficherait 0
    !Config         puis 1, 2, 3 ans, soit un decalage permanent de AGE_STAND_ESTAB annees.
    !Config         La refonte du 2026-07-25 avait retire AGE_RECRUIT_OFFSET en pariant que
    !Config         « la dynamique absorbe le decalage » -- elle ne l'absorbe pas, le
    !Config         diametre initial reste (constat utilisateur 2026-07-29).
    !Config         Purement DIAGNOSTIQUE : age_stand_bm n'est lu par aucune decision.
    !Config Units = [an]
    age_stand_estab = 5.
    CALL getin_p('AGE_STAND_ESTAB', age_stand_estab)
    IF (age_stand_estab < zero) THEN
       WRITE(numout,*) 'AGE_STAND_ESTAB negatif : ', age_stand_estab
       CALL ipslerr_p(3,'constantes','AGE_STAND_ESTAB doit etre >= 0','','')
    ENDIF
    !
    !
    !
    !Config Key   = DIA_ROTATION_TOL
    !Config Desc  = Lower tolerance (fraction) on the clearcut diameter
    !Config If    = OK_STOMATE
    !Config Def   = 0.175
    !Config Units = [-]
    dia_rotation_tol = 0.175
    CALL getin_p('DIA_ROTATION_TOL', dia_rotation_tol)
    !
    !
    ! Guillaume M. -- TARGET_ROTATION_AGE and TARGET_ROTATION_AGE_FILE are REMOVED.
    ! target_rotation_age is now filled by set_management_intensity from the class
    ! fractions: one source, one file, no competing scalar branch.
    IF (ok_management_intensity .AND. TRIM(management_intensity_file) == 'NONE') THEN
       CALL ipslerr_p(3, 'constantes', 'OK_MANAGEMENT_INTENSITY=y requires MANAGEMENT_INTENSITY_FILE', &
            'no class-fraction map to derive the management parameters from', '')
    ENDIF
    !
    ! Guillaume M. -- PROGRESSIVE_HARVEST: the clearcut withdraws a rotation-paced FRACTION
    ! of the mature slot instead of emptying it at once.
    ! See design/MODULE_DESIGN_PROGRESSIVE_HARVEST.md.
    !
    !Config Key   = OK_PROGRESSIVE_HARVEST
    !Config Desc  = Clearcut a rotation-paced fraction of the mature slot instead of emptying it
    !Config If    = OK_MANAGEMENT_INTENSITY
    !Config Def   = n
    !Config Help  = In a forest managed on a rotation R, a fraction 1/R of the area is
    !Config         renewed each year. When TRUE, the CONDITION-2 clearcut withdraws
    !Config         f_cut = MIN(1, (A_forest/R) / A_mature) of the eligible slot area
    !Config         instead of the whole slot; the remainder keeps its age, diameter and
    !Config         biomass and is eligible again next year. The trigger criterion itself
    !Config         is unchanged -- only the action becomes fractional. Corrects together
    !Config         the age-class patchwork, the structural over-harvest (HARVEST_AREA
    !Config         16-26 %/yr vs 1/R ~ 1.2 %/yr) and the uncontrolled rotation.
    !Config         FALSE = all-or-nothing (bit-neutral).
    !Config         REQUIRES OK_MANAGEMENT_INTENSITY=y, the source of R via
    !Config         target_rotation_age. This is enforced below.
    !Config         [CATEGORIE : OPTION SCIENTIFIQUE -- choix de modelisation, permanent]
    !Config Units = [FLAG]
    ok_progressive_harvest = .FALSE.
    CALL getin_p('OK_PROGRESSIVE_HARVEST', ok_progressive_harvest)
    !
    ! Guillaume M. -- f_cut needs no adjustable floor or cap: it is bounded by MIN(un, ...)
    ! where it is computed, and it is a ratio of positive areas. Saturation towards 1 is
    ! watched by test 5.5 of the validation suite, which measures it rather than clipping it.
    !
    ! Guillaume M. -- The dependency is hard: without OK_MANAGEMENT_INTENSITY there is no
    ! target rotation R to pace the harvest from. A silent combination of management flags
    ! is what once inflated the edge length by a factor 2.5, so it is refused here.
    IF (ok_progressive_harvest .AND. .NOT. ok_management_intensity) THEN
       CALL ipslerr_p(3, 'constantes', 'OK_PROGRESSIVE_HARVEST=y requires OK_MANAGEMENT_INTENSITY=y', &
            'the harvest pace comes from target_rotation_age (R)', '')
    ENDIF
    ! Guillaume M. -- Placed HERE and not with its own block: ok_progressive_harvest is
    ! read only at this point (order trap). Without the pre-pass the coupling has no
    ! quota to modulate and would be silently inert -- refused, same policy as above.
    IF (ok_harvest_feed_coupling .AND. .NOT. ok_progressive_harvest) THEN
       CALL ipslerr_p(3, 'constantes', 'OK_HARVEST_FEED_COUPLING=y requires OK_PROGRESSIVE_HARVEST=y', &
            'the coupling modulates the progressive-harvest quota (design doc 4.4)', '')
    ENDIF
    !
    !
    !
    !
    !
    !
    !
    !
    !
    !
    !
    !
    !
    !Config Key   = CC_N_FLOOR
    !Config Desc  = Effectif sous lequel une classe de circonference n'a plus sa propre taille
    !Config If    = OK_STOMATE
    !Config Def   = 0.0
    !Config Help  = La normalisation de la biomasse par individu, dans mortality_clean, divise
    !Config         par l'effectif de la classe. Les poids de weibull_class_dist laissent
    !Config         presque rien dans les classes extremes, donc le diviseur marche vers zero
    !Config         et la biomasse par arbre explose : `[MESURE]` une classe portant
    !Config         0,02 tige/ha atteignait 94 cm de diametre, en croissant de 1,9 cm/an quand
    !Config         le peuplement croit de 0,5. Le tri qui suit rangeait ces valeurs proprement
    !Config         en haut de grille, ce qui les faisait lire comme une queue de peuplement.
    !Config         Au-dessous du seuil la classe ADOPTE la taille de sa voisine.
    !Config         /!\ Ses tiges restent en place. Vider la classe a ete essaye et casse
    !Config         l'invariant que mortality_clean existe pour tenir : growth_fun_all divise
    !Config         b_inc_tot par circ_class_n CLASSE PAR CLASSE (expression tableau), garde
    !Config         seulement sur la SOMME, donc une classe exactement vide rend Inf x 0 = NaN
    !Config         sur tout le PFT des le lendemain.
    !Config         Ordre de grandeur utile : 1e-4 ind m-2, soit 1 tige/ha.
    !Config         /!\ VALEUR NEUTRE = LE DEFAUT (0.0) => aucune classe n'est sous le seuil.
    !Config         [CATEGORIE : OPTION SCIENTIFIQUE -- choix de modelisation, permanent]
    !Config Units = [ind m-2]
    cc_n_floor = zero
    CALL getin_p('CC_N_FLOOR', cc_n_floor)
    IF (cc_n_floor < zero) THEN
       CALL ipslerr_p(3,'constantes','CC_N_FLOOR ne peut pas etre negatif', &
            'un effectif plancher negatif ne veut rien dire','')
    ENDIF
    !
    !Config Key   = OK_EDGE_FROM_AGE_CLASS
    !Config Desc  = Derive the forest edge length from the young age-class open area
    !Config If    = OK_AED_FEEDBACK
    !Config Cat   = SCIENTIFIC OPTION
    !Config Def   = n
    !Config Help  = Chooses HOW the open area that generates forest edge is accounted for.
    !Config         TRUE  : edge comes from A_open, the young forest age-class area
    !Config                 (sum of o_c * veget_max). The age-class machinery conserves
    !Config                 area, so A_open <= forest area by construction: no double
    !Config                 counting, and the recovery time is the age-class residence
    !Config                 time rather than a prescribed constant.
    !Config         FALSE : disturbance adds edge which then decays exponentially with
    !Config                 EDGE_FEEDBACK_TAU_RECOVERY (phenomenological relaxation).
    !Config         Does not modify veget_max. Requires OK_AED_FEEDBACK=y.
    !Config         See design/MODULE_DESIGN_AED_REGROWTH_V2.md.
    !Config Units = [FLAG]
    ok_edge_from_age_class = .FALSE.
    CALL getin_p('OK_EDGE_FROM_AGE_CLASS', ok_edge_from_age_class)
    !
    !Config Key   = AED_RECRUIT_TIME
    !Config Desc  = Years to canopy reclosure (class-0 conveyor length)
    !Config If    = OK_EDGE_FROM_AGE_CLASS
    !Config Def   = 11
    !Config Help  = Residence time of a disturbed (stand-replaced) patch in the
    !Config         class-0 open pool before it graduates back to forest. Replaces
    !Config         tau_rec. Default 11 yr from the EFDA box-model pre-check and
    !Config         AED_SPINUP §10bis (both converge on ~11 yr). Boreal stands may
    !Config         need a longer value.
    !Config Units = [yr]
    n_recruit = 11
    CALL getin_p('AED_RECRUIT_TIME', n_recruit)
    !
    !Config Key   = AED_RECRUIT_TIME_MAX
    !Config Desc  = Max class-0 conveyor length (= max per-PFT recovery time)
    !Config If    = OK_EDGE_FROM_AGE_CLASS
    !Config Def   = 30
    !Config Help  = Dimensions class0_area (3rd dim). The per-PFT recovery times
    !Config         recruit_time_pft (built at runtime from MTC: sclerophyll 20,
    !Config         boreal 30, deciduous 10, ...) must not exceed this. Default 30
    !Config         covers the slowest (boreal). Cells graduate at their veget_max-
    !Config         weighted effective time n_recruit_eff(ji), not this maximum.
    !Config Units = [yr]
    n_recruit_max = 30
    CALL getin_p('AED_RECRUIT_TIME_MAX', n_recruit_max)
    n_recruit_max = MAX(n_recruit_max, n_recruit)
    !
    !
    !Config Key   = AED_MIX_TAU
    !Config Desc  = v2 leaky agent-flux memory (yr) for the per-agent edge mix
    !Config If    = OK_EDGE_FROM_AGE_CLASS
    !Config Def   = 11
    !Config Units = [yr]
    aed_mix_tau = 11._r_std
    CALL getin_p('AED_MIX_TAU', aed_mix_tau)
    aed_mix_tau = MAX(aed_mix_tau, 1._r_std)
    !
    ! Guillaume M. -- AED_OPEN_W1..W4 removed. W2/W3/W4 in fact described the
    ! forest/non-forest boundary of a mature stand, which the meso term now produces
    ! geometrically (A*4*f*(1-f)/r) with no tuned parameter, and W1 was 1 by construction.
    ! A_open therefore keeps age class 1 only, the sole genuine post-disturbance gap.
    ! The ipslerr guard forbidding W3/W4 > 0 under OK_EDGE_3TERMS goes with them.
    !
    !Config Key   = AED_BG_FRAC
    !Config Desc  = v2 background-edge fraction: young/open forest (A_open) keeps a
    !Config         standing-fragmentation edge even with no current disturbance flux
    !Config If    = OK_EDGE_FROM_AGE_CLASS
    !Config Def   = 1.0
    !Config Help  = Closes the temperate EL_dyn=0 hole (piste B): cells with A_open>0 but
    !Config         flux_tot~0 get edge = aed_bg_frac * A_open * alpha_bg (climatological
    !Config         pan-EU agent mix). 0 reproduces the legacy flux-gated behaviour.
    !Config Units = [-]
    aed_bg_frac = 1.0_r_std
    CALL getin_p('AED_BG_FRAC', aed_bg_frac)
    aed_bg_frac = MAX(aed_bg_frac, zero)
    !
    ! Guillaume M. -- OK_AED_EDGE_TWO_TERM and AED_ALPHA_FRAG removed: the landscape term
    ! and the MESO term of edge_length_base describe the SAME forest/non-forest interface,
    ! both proportional to f*(1-f)*area and both feeding edge_length_dyn, so the interface
    ! was counted twice. The meso term, anchored on a physical grain per MTC, is the one
    ! kept. See design/MODULE_DESIGN_AED_EDGE_TWO_TERM.md.
    !
    !Config Key   = AED_ALPHA_REG
    !Config Desc  = Weight of the regeneration-gap term (age class 1) in A_open
    !Config If    = OK_EDGE_FROM_AGE_CLASS
    !Config Def   = 1.0
    !Config Help  = A_open = AED_ALPHA_REG * A_class1, optionally weighted by the stand
    !Config         openness (OK_AED_EDGE_RDI_WEIGHT). Age class 1 is the only true
    !Config         internal gap, so it keeps weight 1 (the value the removed AED_OPEN_W1
    !Config         carried). Only the AREA is set here -- the area-to-length conversion by
    !Config         the agent-weighted alpha is unchanged. TO BE CALIBRATED against EFDA
    !Config         together with MI_EDGE_SHAPE_COEF, which shares the same anchor.
    !Config Units = [-]
    aed_alpha_reg = 1.0_r_std
    CALL getin_p('AED_ALPHA_REG', aed_alpha_reg)
    !
    ! Guillaume M. -- Coherence guards for AED_REGROWTH: it consumes the dA_* accumulators
    ! and edge_length_dyn owned by AED_FEEDBACK, and needs at least a 1-year conveyor.
    IF (ok_edge_from_age_class .AND. .NOT. ok_aed_feedback) THEN
       CALL ipslerr_p(3, 'config_stomate_parameters', &
            'OK_EDGE_FROM_AGE_CLASS=y requires OK_AED_FEEDBACK=y', &
            'edge_length_dyn and the dA_* accumulators are owned by AED_FEEDBACK.', '')
    ENDIF
    IF (ok_edge_from_age_class .AND. n_recruit < 1) THEN
       CALL ipslerr_p(3, 'config_stomate_parameters', &
            'AED_RECRUIT_TIME must be >= 1 year when OK_EDGE_FROM_AGE_CLASS=y', '', '')
    ENDIF
    !
    ! Guillaume M. -- calculate_aed always consumes edge_length_dyn, which only
    ! update_edge_length fills and which belongs to AED_FEEDBACK. Without it the array stays
    ! zero and the edge modules would silently run as if there were no edge. Refuse it.
    IF ((ok_aed_light .OR. ok_aed_humign .OR. ok_aed_fuel .OR. ok_aed_size .OR. ok_aed_wind) &
         .AND. .NOT. ok_aed_feedback) THEN
       CALL ipslerr_p(3, 'config_stomate_parameters', &
            'les modules OK_AED_* requierent OK_AED_FEEDBACK=y', &
            'calculate_aed consomme edge_length_dyn, alimentee par AED_FEEDBACK seul', '')
    ENDIF
    !-
    ! allocation parameters
    !-
    !Config Key   = RESERVE_TIME_TREE 
    !Config Desc  = maximum time during which reserve is used (trees) 
    !Config If    = OK_STOMATE 
    !Config Def   = 30.
    !Config Help  = 
    !Config Units = [days]
    reserve_time_tree = 30.
    CALL getin_p('RESERVE_TIME_TREE',reserve_time_tree)
    !
    !Config Key   = RESERVE_TIME_GRASS 
    !Config Desc  = maximum time during which reserve is used (grasses) 
    !Config If    = OK_STOMATE 
    !Config Def   = 20. 
    !Config Help  = 
    !Config Units = [days]
    reserve_time_grass = 20.
    CALL getin_p('RESERVE_TIME_GRASS',reserve_time_grass)

    !-
    ! data parameters
    !
    !Config Key   = PRECIP_CRIT 
    !Config Desc  = minimum precip
    !Config If    = OK_STOMATE 
    !Config Def   = 100.
    !Config Help  = 
    !Config Units = [mm/year]
    precip_crit = 100. 
    CALL getin_p('PRECIP_CRIT',precip_crit)
    !
    !Config Key   = GDD_CRIT_ESTAB
    !Config Desc  = minimum gdd for establishment of saplings
    !Config If    = OK_STOMATE 
    !Config Def   = 150. 
    !Config Help  = 
    !Config Units = [-]
    gdd_crit_estab = 150.
    CALL getin_p('GDD_CRIT_ESTAB',gdd_crit_estab)
    !
    !Config Key   = SPRING_DAYS_MAX
    !Config Desc  = Maximum number of days during which we watch for possible spring frost damage
    !Config If    = OK_STOMATE 
    !Config Def   = 40
    !Config Help  = 
    !Config Units = [days]
    spring_days_max = 40
    CALL getin_p('SPRING_DAYS_MAX',spring_days_max)
    !
    !Config Key   = FPC_CRIT
    !Config Desc  = critical fpc, needed for light competition and establishment
    !Config If    = OK_STOMATE 
    !Config Def   = 0.95
    !Config Help  = 
    !Config Units = [-]
    fpc_crit = 0.95
    CALL getin_p('FPC_CRIT',fpc_crit)
    !
    !Config Key   = ALPHA_GRASS
    !Config Desc  = sapling characteristics : alpha's
    !Config If    = OK_STOMATE 
    !Config Def   = 0.5
    !Config Help  = 
    !Config Units = [-]
    alpha_grass = 0.5
    CALL getin_p('ALPHA_GRASS',alpha_grass)
    !
    !Config Key   = ALPHA_TREE
    !Config Desc  = sapling characteristics : alpha's 
    !Config If    = OK_STOMATE 
    !Config Def   = 1.
    !Config Help  = 
    !Config Units = [-]
    alpha_tree = 1.
    CALL getin_p('ALPHA_TREE',alpha_tree)
    !
    !Config Key   = STRUCT_TO_LEAVES
    !Config Desc  = Fraction of structural carbon in grass and crops as a share of the leaf
    ! carbon pool. Only used for grasses and crops (thus NOT for trees)
    !Config If    = OK_STOMATE 
    !Config Def   = 0.05
    !Config Help  = NOTE: the line using this variable is
    !Config         commented out in r5976, and thus this variable is
    !Config         not used.
    !Config Units = [-]
    struct_to_leaves = 0.05
    CALL getin_p(' STRUCT_TO_LEAVES',struct_to_leaves)
    !
    !Config Key   = LABILE_TO_TOTAL
    !Config Desc  = Fraction of the labile pool in trees, grasses and crops
    !Config If    = OK_STOMATE 
    !Config Def   = 0.01
    !Config Help  = Fraction of the labile pool in trees, grasses and crops as a share of
    !Config Help  = the total carbon pool (accounting for the N-content of the different tissues).
    !Config Units = [-]
    labile_to_total = 0.01 
    CALL getin_p('LABILE_TO_TOTAL',labile_to_total)
    ! 
    !Config Key   = TEST_BIOMASS_INIT
    !Config Desc  = Use initial biomass estimate to determine size of plant at establishment
    !Config If    = OK_STOMATE 
    !Config Def   = n
    !Config Help  = If set to true init_biomass used to determine size of
    !               plants at establishment. If set to false, init_lai is used.
    !Config Units = [FLAG]
    test_biomass_init = .TRUE.
    CALL getin_p('TEST_BIOMASS_INIT',test_biomass_init)

    !
    !Config Key   = TAU_WEEK
    !Config Desc  = time scales for phenology and other processes
    !Config If    = OK_STOMATE 
    !Config Def   = 7.
    !Config Help  = 
    !Config Units = [days]
    tau_week = 7.
    CALL getin_p('TAU_WEEK',tau_week)
    !
    !Config Key   = TAU_HUM_MONTH
    !Config Desc  = time scales for phenology and other processes
    !Config If    = OK_STOMATE 
    !Config Def   = 20. 
    !Config Help  = 
    !Config Units = [days]
    tau_hum_month = 20. 
    CALL getin_p('TAU_HUM_MONTH',tau_hum_month)
    !
    !Config Key   = TAU_HUM_WEEK
    !Config Desc  = time scales for phenology and other processes
    !Config If    = OK_STOMATE 
    !Config Def   = 7.
    !Config Help  = 
    !Config Units = [days]
    tau_hum_week = 7.
    CALL getin_p('TAU_HUM_WEEK',tau_hum_week)
    !
    !Config Key   = TAU_PRECIP_MONTH
    !Config Desc  = time scales for phenology and other processes
    !Config If    = OK_STOMATE 
    !Config Def   = 30.
    !Config Help  = 
    !Config Units = [days]
    tau_precip_month = 20.
    CALL getin_p('TAU_PRECIP_MONTH',tau_precip_month)
    !
    !Config Key   = TAU_T2M_MONTH
    !Config Desc  = time scales for phenology and other processes
    !Config If    = OK_STOMATE 
    !Config Def   = 20.
    !Config Help  = 
    !Config Units = [days]
    tau_t2m_month = 20.
    CALL getin_p('TAU_T2M_MONTH',tau_t2m_month)
    !
    !Config Key   = TAU_T2M_WEEK
    !Config Desc  = time scales for phenology and other processes
    !Config If    = OK_STOMATE 
    !Config Def   = 7.
    !Config Help  = 
    !Config Units = [days]
    tau_t2m_week = 7.
    CALL getin_p('TAU_T2M_WEEK',tau_t2m_week)
    !
    !Config Key   = TAU_TSOIL_MONTH 
    !Config Desc  = time scales for phenology and other processes
    !Config If    = OK_STOMATE 
    !Config Def   = 20. 
    !Config Help  = 
    !Config Units = [days]
    tau_tsoil_month = 20.
    CALL getin_p('TAU_TSOIL_MONTH',tau_tsoil_month)
    !
    !Config Key   = TAU_GPP_WEEK 
    !Config Desc  = time scales for phenology and other processes
    !Config If    = OK_STOMATE 
    !Config Def   = 7. 
    !Config Help  = 
    !Config Units = [days]
    tau_gpp_week = 7. 
    CALL getin_p('TAU_GPP_WEEK',tau_gpp_week)
    !
    !Config Key   = TAU_GPP_YEAR 
    !Config Desc  = time scales to average gpp over a year
    !Config If    = OK_STOMATE 
    !Config Def   = 365
    !Config Help  = gpp_year is averaged over tau_gpp_year
    !Config Units = [days]
    tau_gpp_year = 365.
    CALL getin_p('TAU_GPP_YEAR',tau_gpp_year)
    !
    !Config Key   = TAU_GPP_SELF_THIN 
    !Config Desc  = time scales to average gpp for self thinning
    !Config If    = OK_STOMATE 
    !Config Def   = 30
    !Config Help  = pre_indust_ref_gpp is averaged over
    !               tau_gpp_self_thin and used to adjust
    !               alpha_self_thinning (which is a proxy for
    !               resource availability
    !Config Units = [years]
    tau_gpp_self_thin = 30.
    CALL getin_p('TAU_GPP_SELF_THIN',tau_gpp_self_thin)
    !
    !Config Key   = TAU_SUGARLOAD_WEEK 
    !Config Desc  = time scales for calculating sugar loading
    !Config If    = OK_STOMATE 
    !Config Def   = 7. 
    !Config Help  = 
    !Config Units = [days]
    tau_sugarload_week = 7. 
    CALL getin_p('TAU_SUGARLOAD_WEEK',tau_sugarload_week)
    !
    !Config Key   = TAU_GDD
    !Config Desc  = time scales for phenology and other processes
    !Config If    = OK_STOMATE 
    !Config Def   = 40. 
    !Config Help  = 
    !Config Units = [days]
    tau_gdd = 40.
    CALL getin_p('TAU_GDD',tau_gdd)
    !
    !Config Key   = TAU_NGD
    !Config Desc  = time scales for phenology and other processes
    !Config If    = OK_STOMATE 
    !Config Def   = 43.
    !Config Help  = 
    !Config Units = [days]
    tau_ngd = 43.
    CALL getin_p('TAU_NGD',tau_ngd)
    !
    !Config Key   = COEFF_TAU_LONGTERM
    !Config Desc  = time scales for phenology and other processes
    !Config If    = OK_STOMATE 
    !Config Def   = 3. 
    !Config Help  = 
    !Config Units = [days]
    coeff_tau_longterm = 3.
    CALL getin_p('COEFF_TAU_LONGTERM',coeff_tau_longterm)
    !-
    !
    !Config Key   = BM_SAPL_CARBRES 
    !Config Desc  = 
    !Config If    = OK_STOMATE 
    !Config Def   = 5. 
    !Config Help  = 
    !Config Units = [-]
    bm_sapl_carbres = 5.
    CALL getin_p('BM_SAPL_CARBRES',bm_sapl_carbres)
    !
    !Config Key   = BM_SAPL_SAPABOVE
    !Config Desc  = 
    !Config If    = OK_STOMATE
    !Config Def   = 0.5 
    !Config Help  = 
    !Config Units = [-]
    bm_sapl_sapabove = 0.5 
    CALL getin_p('BM_SAPL_SAPABOVE',bm_sapl_sapabove)
    !
    !Config Key   = BM_SAPL_HEARTABOVE 
    !Config Desc  = 
    !Config If    = OK_STOMATE 
    !Config Def   = 2.
    !Config Help  = 
    !Config Units = [-]
    bm_sapl_heartabove = 2.
    CALL getin_p('BM_SAPL_HEARTABOVE',bm_sapl_heartabove)
    !
    !Config Key   = BM_SAPL_HEARTBELOW 
    !Config Desc  = 
    !Config If    = OK_STOMATE 
    !Config Def   = 2. 
    !Config Help  = 
    !Config Units = [-]
    bm_sapl_heartbelow = 2.
    CALL getin_p('BM_SAPL_HEARTBELOW',bm_sapl_heartbelow)

    !Config Key   = BM_SAPL_LABILE 
    !Config Desc  = 
    !Config If    = OK_STOMATE 
    !Config Def   = 5. 
    !Config Help  = 
    !Config Units = [-]
    bm_sapl_labile = 5.
    CALL getin_p('BM_SAPL_LABILE',bm_sapl_labile)

    !Config Key   = INIT_SAPL_MASS_LABILE
    !Config Desc  = 
    !Config If    = OK_STOMATE 
    !Config Def   = 5. 
    !Config Help  = 
    !Config Units = [-]
    init_sapl_mass_labile = 5.
    CALL getin_p('INIT_SAPL_MASS_LABILE',init_sapl_mass_labile)

    !Config Key   = INIT_SAPL_MASS_LEAF_NAT
    !Config Desc  = 
    !Config If    = OK_STOMATE 
    !Config Def   = 0.1 
    !Config Help  = 
    !Config Units = [-]
    init_sapl_mass_leaf_nat = 0.1
    CALL getin_p('INIT_SAPL_MASS_LEAF_NAT',init_sapl_mass_leaf_nat)
    !
    !Config Key   = INIT_SAPL_MASS_LEAF_AGRI
    !Config Desc  = 
    !Config If    = OK_STOMATE 
    !Config Def   = 1. 
    !Config Help  = 
    !Config Units = [-]
    init_sapl_mass_leaf_agri = 1. 
    CALL getin_p('INIT_SAPL_MASS_LEAF_AGRI',init_sapl_mass_leaf_agri)
    !
    !Config Key   = INIT_SAPL_MASS_CARBRES
    !Config Desc  = 
    !Config If    = OK_STOMATE 
    !Config Def   = 5. 
    !Config Help  = 
    !Config Units = [-]
    init_sapl_mass_carbres = 5. 
    CALL getin_p('INIT_SAPL_MASS_CARBRES',init_sapl_mass_carbres)
    !
    !Config Key   = INIT_SAPL_MASS_ROOT
    !Config Desc  = 
    !Config If    = OK_STOMATE 
    !Config Def   = 0.1 
    !Config Help  = 
    !Config Units = [-]
    init_sapl_mass_root = 0.1
    CALL getin_p('INIT_SAPL_MASS_ROOT',init_sapl_mass_root)
    !
    !Config Key   = INIT_SAPL_MASS_FRUIT
    !Config Desc  = 
    !Config If    = OK_STOMATE 
    !Config Def   = 0.3 
    !Config Help  = 
    !Config Units = [-]
    init_sapl_mass_fruit = 0.3 
    CALL getin_p('INIT_SAPL_MASS_FRUIT',init_sapl_mass_fruit)
    !
    !Config Key   = CN_SAPL_INIT 
    !Config Desc  = 
    !Config If    = OK_STOMATE 
    !Config Def   = 0.5 
    !Config Help  = 
    !Config Units = [-]
    cn_sapl_init = 0.5 
    CALL getin_p('CN_SAPL_INIT',cn_sapl_init)
    !
    !Config Key   = MIGRATE_TREE 
    !Config Desc  = 
    !Config If    = OK_STOMATE 
    !Config Def   = 10000.
    !Config Help  = 
    !Config Units = [m/year]
    migrate_tree = 10.*1.E3 
    CALL getin_p('MIGRATE_TREE',migrate_tree)
    !
    !Config Key   = MIGRATE_GRASS
    !Config Desc  = 
    !Config If    = OK_STOMATE 
    !Config Def   = 10000.
    !Config Help  = 
    !Config Units = [m/year]
    migrate_grass = 10.*1.E3 
    CALL getin_p('MIGRATE_GRASS',migrate_grass)
    !
    !Config Key   = LAI_INITMIN_TREE
    !Config Desc  = 
    !Config If    = OK_STOMATE 
    !Config Def   = 0.3
    !Config Help  = 
    !Config Units = [m^2/m^2]
    lai_initmin_tree = 0.3 
    CALL getin_p('LAI_INITMIN_TREE',lai_initmin_tree)
    !
    !Config Key   = LAI_INITMIN_GRASS 
    !Config Desc  = 
    !Config If    = OK_STOMATE 
    !Config Def   = 0.1
    !Config Help  = 
    !Config Units = [m^2/m^2]
    lai_initmin_grass = 0.1 
    CALL getin_p('LAI_INITMIN_GRASS',lai_initmin_grass)
    !
    !Config Key   = DIA_COEFF
    !Config Desc  = 
    !Config If    = OK_STOMATE 
    !Config Def   = 4., 0.5
    !Config Help  = 
    !Config Units = [-]
    dia_coeff = (/ 4., 0.5 /)    
    CALL getin_p('DIA_COEFF',dia_coeff)
    !
    !Config Key   = BM_SAPL_LEAF
    !Config Desc  = 
    !Config If    = OK_STOMATE 
    !Config Def   = 4., 4., 0.8, 5. 
    !Config Help  = 
    !Config Units = [-]
    bm_sapl_leaf = (/ 4., 4., 0.8, 5./)
    CALL getin_p('BM_SAPL_LEAF',bm_sapl_leaf)

    !-
    ! litter parameters
    !-

    !Config Key   = CN
    !Config Desc  = C/N ratio
    !Config If    = OK_STOMATE 
    !Config Def   = 40., 40., 40., 40., 40., 40., 40., 40.
    !Config Help  = 
    !Config Units = [-]
    CN_fix = (/ 40., 40., 40., 40., 40., 40., 40., 40., 40. /) 
    CALL getin_p('CN',CN_fix)

    !Config Key   = FRAC_SOIL_STRUCT_SUA
    !Config Desc  = frac_soil(istructural,isurface,iabove)
    !Config If    = OK_STOMATE 
    !Config Def   = 0.55
    !Config Help  = 
    !Config Units = [-]
    frac_soil_struct_sua = 0.4
    CALL getin_p('FRAC_SOIL_STRUCT_SUA',frac_soil_struct_sua)

    !Config Key   = FRAC_SOIL_STRUCT_AA
    !Config Desc  = frac_soil(istructural,iactive,iabove)
    !Config If    = OK_STOMATE 
    !Config Def   = 0.4
    !Config Help  = 
    !Config Units = [-]
    frac_soil_struct_aa = 0.4
    CALL getin_p('FRAC_SOIL_STRUCT_AA',frac_soil_struct_aa)

    !Config Key   = FRAC_SOIL_METAB_SUA 
    !Config Desc  = frac_soil(imetabolic,isurface,iabove) 
    !Config If    = OK_STOMATE 
    !Config Def   = 0.4 
    !Config Help  = 
    !Config Units = [-]
    frac_soil_metab_sua = 0.4
    CALL getin_p('FRAC_SOIL_METAB_SUA',frac_soil_metab_sua)

    !Config Key   = FRAC_SOIL_METAB_AA 
    !Config Desc  = frac_soil(imetabolic,iactive,iabove) 
    !Config If    = OK_STOMATE 
    !Config Def   = 0.4 
    !Config Help  = 
    !Config Units = [-]
    frac_soil_metab_aa = 0.4
    CALL getin_p('FRAC_SOIL_METAB_AA',frac_soil_metab_aa)

    !Config Key   = TURN_METABOLIC
    !Config Desc  = 
    !Config If    = OK_STOMATE 
    !Config Def   = 0.066
    !Config Help  = 
    !Config Units = [years]
    turn_metabolic = 15
    CALL getin_p('TURN_METABOLIC',turn_metabolic)

    !Config Key   = TURN_STRUCT 
    !Config Desc  = 
    !Config If    = OK_STOMATE 
    !Config Def   = 0.245 
    !Config Help  = 
    !Config Units = [years]
    turn_struct = 4
    CALL getin_p('TURN_STRUCT',turn_struct)

    !Config Key   = TURN_WOODY
    !Config Desc  = 
    !Config If    = OK_STOMATE 
    !Config Def   = 1.33
    !Config Help  = unknown where this one comes from (not Parton1993)
    !Config Units = [years]
    turn_woody = 1.33  
    CALL getin_p('TURN_WOODY',turn_woody)

    !Config Key   = SNAG_FRAC
    !Config Desc  = Fraction of woody litter input going to snags
    !Config If    = OK_STOMATE 
    !Config Def   = 0.3
    !Config Help  = This should be refined with _cc objects for snags, counting
    !Config Help  = snags per death scenario, each cause of tree death producing
    !Config Help  = different proportions of snags
    !Config Help  = Proportion of the woody litter input going into snag 
    !Config Help  = the rest going in the iwoody pool, drawn from a synthesis :
    !Config Help  = Malysheva et al. 2019 over Russia,
    !Config Help  = Motta et al. 2006 over Italian Alps,
    !Config Help  = Rahman et al. 2008 over an austrian forest,
    !Config Help  = Caron et al. 2009 over Fennoscandia
    !Config Units = [-]
    snag_frac = 0.3225
    CALL getin_p('SNAG_FRAC',snag_frac)
    
    !Config Key   = TURN_SNAG
    !Config Desc  = Decomposition / Turnover of snag pool to soil organic matter
    !Config If    = OK_STOMATE 
    !Config Def   = 0.00133
    !Config Help  = Turnover of the snag litter pool. This value should be lower
    !Config Help  = than the TURN_SNAG_FALL. If not, the interest of introducing
    !Config       = a snag pool disappears
    !Config Units = [years]
    turn_snag = 0.0133 
    CALL getin_p('TURN_SNAG',turn_snag)
    
    !Config Key   = TURN_SNAG_FALL
    !Config Desc  = Snag fall / Turnover of snag pool to woody litter pool
    !Config If    = OK_STOMATE 
    !Config Def   = 0.133
    !Config Help  = This is the rate for a flux from litter to litter 
    !Config       = (snag to woody)
    !Config Units = [years]
    turn_snag_fall = 0.133
    CALL getin_p('TURN_SNAG_FALL',turn_snag_fall)

    !Config Key   = METABOLIC_REF_FRAC
    !Config Desc  =
    !Config If    = OK_STOMATE 
    !Config Def   = 0.85  
    !Config Help  = 
    !Config Units = [-]
    metabolic_ref_frac = 0.85
    CALL getin_p('METABOLIC_REF_FRAC',metabolic_ref_frac)

    !Config Key   = Z_DECOMP
    !Config Desc  = scaling depth for soil activity
    !Config If    = OK_STOMATE 
    !Config Def   = 0.2
    !Config Help  = 
    !Config Units = [m]
    z_decomp = 0.2 
    CALL getin_p('Z_DECOMP',z_decomp)

    !Config Key   = FRAC_SOIL_STRUCT_A 
    !Config Desc  = frac_soil(istructural,iactive,ibelow)
    !Config If    = OK_STOMATE 
    !Config Def   = 0.45
    !Config Help  = 
    !Config Units = [-]
    frac_soil_struct_ab = 0.45 
    CALL getin_p('FRAC_SOIL_STRUCT_AB',frac_soil_struct_ab)
    !
    !Config Key   = FRAC_SOIL_STRUCT_SA
    !Config Desc  = frac_soil(istructural,islow,iabove)
    !Config If    = OK_STOMATE
    !Config Def   = 0.7  
    !Config Help  = 
    !Config Units = [-]
    frac_soil_struct_sa = 0.7
    CALL getin_p('FRAC_SOIL_STRUCT_SA',frac_soil_struct_sa)
    !
    !Config Key   = FRAC_SOIL_STRUCT_SB
    !Config Desc  = frac_soil(istructural,islow,ibelow) 
    !Config If    = OK_STOMATE 
    !Config Def   = 0.7  
    !Config Help  = 
    !Config Units = [-]
    frac_soil_struct_sb = 0.7 
    CALL getin_p('FRAC_SOIL_STRUCT_SB',frac_soil_struct_sb)
    !
    !Config Key   = FRAC_SOIL_METAB_AB 
    !Config Desc  = frac_soil(imetabolic,iactive,ibelow)
    !Config If    = OK_STOMATE 
    !Config Def   = 0.45  
    !Config Help  = 
    !Config Units = [-]
    frac_soil_metab_ab = 0.45
    CALL getin_p('FRAC_SOIL_METAB_AB',frac_soil_metab_ab)
    !
    !Config Key   = METABOLIC_LN_RATIO
    !Config Desc  = 
    !Config If    = OK_STOMATE 
    !Config Def   = 0.018  
    !Config Help  = 
    !Config Units = [-]
    metabolic_LN_ratio = 0.018 
    CALL getin_p('METABOLIC_LN_RATIO',metabolic_LN_ratio) 
    !
    !Config Key   = SOIL_Q10
    !Config Desc  = 
    !Config If    = OK_STOMATE 
    !Config Def   = 0.69 (=ln2)
    !Config Help  = 
    !Config Units = [-]
    soil_Q10 = 0.69 
    CALL getin_p('SOIL_Q10',soil_Q10)
    !
    !Config Key   = SOIL_Q10_UPTAKE
    !Config Desc  = 
    !Config If    = OK_STOMATE 
    !Config Def   = 0.69 (=ln2)
    !Config Help  = 
    !Config Units = [-]
    soil_Q10_uptake = 0.45 
    CALL getin_p('SOIL_Q10_UPTAKE',soil_Q10_uptake)
    !
    !Config Key   = TSOIL_REF
    !Config Desc  = 
    !Config If    = OK_STOMATE 
    !Config Def   = 30. 
    !Config Help  = 
    !Config Units = [C]
    tsoil_ref = 30. 
    CALL getin_p('TSOIL_REF',tsoil_ref)
    !
    !Config Key   = LITTER_STRUCT_COEF 
    !Config Desc  = 
    !Config If    = OK_STOMATE 
    !Config Def   = 3. 
    !Config Help  = 
    !Config Units = [-]
    litter_struct_coef = 3.
    CALL getin_p('LITTER_STRUCT_COEF',litter_struct_coef)
    !
    !Config Key   = MOIST_COEFF
    !Config Desc  = 
    !Config If    = OK_STOMATE 
    !Config Def   = 1.1, 2.4, 0.29
    !Config Help  = 
    !Config Units = [-]
    moist_coeff = (/ 1.1,  2.4,  0.29 /)
    CALL getin_p('MOIST_COEFF',moist_coeff)
    !
    !Config Key   = MOISTCONT_MIN
    !Config Desc  = minimum soil wetness to limit the heterotrophic respiration
    !Config If    = OK_STOMATE 
    !Config Def   = 0.25
    !Config Help  = 
    !Config Units = [-]
    moistcont_min = 0.25 
    CALL getin_p('MOISTCONT_MIN',moistcont_min)
    !
    !Config Key   = FINEROOTDEPTH_RATIO
    !Config Desc  = the ratio of fine root to overall root e-folding depth (for C inputs)
    !Config If    = OK_STOMATE and OK_SOIL_CARBON_DISCRETIZATION 
    !Config Def   = 0.5
    !Config Help  = 
    !Config Units = [-]
    finerootdepthratio = 0.5
    CALL getin_p('FINEROOTDEPTH_RATIO',finerootdepthratio)
    !
    !Config Key   = ALTROOT_RATIO
    !Config Desc  = the maximum ratio of fine root depth to active layer thickness (for C inputs)
    !Config If    = OK_STOMATE and OK_SOIL_CARBON_DISCRETIZATION 
    !Config Def   = 0.5
    !Config Help  = 
    !Config Units = [-]
    altrootratio = 0.5
    CALL getin_p('ALTROOT_RATIO',altrootratio)
    !
    !
    !Config Key   = FRAC_WOODY
    !Config Desc  = Coefficient for determining the lignin fraction of woody litter
    !Config If    = OK_STOMATE
    !Config Def   = 0.65
    !Config Help  =
    !Config Units = [-]
    frac_woody = 0.65
    CALL getin_p('FRAC_WOODY',frac_woody)
    !
    !Config Key   = FRAC_SNAG
    !Config Desc  = Coefficient for determining the lignin fraction of snag litter
    !Config If    = OK_STOMATE
    !Config Def   = 0.65
    !Config Help  =
    !Config Units = [-]
    frac_snag = 0.65 
    CALL getin_p('FRAC_SNAG',frac_snag)
      

    !-
    ! lpj parameters
    !-
    !
    !Config Key   = FRAC_TURNOVER_DAILY 
    !Config Desc  = 
    !Config If    = OK_STOMATE 
    !Config Def   = 0.55
    !Config Help  = 
    !Config Units = [-]
    frac_turnover_daily = 0.55 
    CALL getin_p('FRAC_TURNOVER_DAILY',frac_turnover_daily)   

    !-
    ! npp parameters
    !-
    !
    !Config Key   = TAX_MAX
    !Config Desc  = maximum fraction of allocatable biomass used for maintenance respiration
    !Config If    = OK_STOMATE 
    !Config Def   = 0.8
    !Config Help  = 
    !Config Units = [-]
    tax_max = 0.8
    CALL getin_p('TAX_MAX',tax_max) 

    !-
    ! phenology parameters
    !-
    !Config Key   = MIN_GROWTHINIT_TIME 
    !Config Desc  = minimum time since last beginning of a growing season
    !Config If    = OK_STOMATE 
    !Config Def   = 300. 
    !Config Help  = 
    !Config Units = [days]
    min_growthinit_time = 300.
    CALL getin_p('MIN_GROWTHINIT_TIME',min_growthinit_time)
    !
    !Config Key   = GDDNCD_REF 
    !Config Desc  = 
    !Config If    = OK_STOMATE 
    !Config Def   = 964. 
    !Config Help  = 
    !Config Units = [-]
    gddncd_ref = 964.
    CALL getin_p('GDDNCD_REF',gddncd_ref)
    !
    !Config Key   = GDDNCD_CURVE
    !Config Desc  = 
    !Config If    = OK_STOMATE 
    !Config Def   = 0.0058 
    !Config Help  = 
    !Config Units = [-]
    gddncd_curve = 0.0058
    CALL getin_p('GDDNCD_CURVE',gddncd_curve)
    !
    !Config Key   = GDDNCD_OFFSET
    !Config Desc  = 
    !Config If    = OK_STOMATE 
    !Config Def   = 12.8. 
    !Config Help  = 
    !Config Units = [-]
    gddncd_offset = 12.8   
    CALL getin_p('GDDNCD_OFFSET',gddncd_offset)
    !-
    ! respiration parameters
    !-
    !
    !Config Key   = MAINT_RESP_MIN_VMAX
    !Config Desc  = 
    !Config If    = OK_STOMATE 
    !Config Def   = 0.3
    !Config Help  = 
    !Config Units = [-]
    maint_resp_min_vmax = 0.3
    CALL getin_p('MAINT_RESP_MIN_VMAX',maint_resp_min_vmax)  
    !
    !Config Key   = MAINT_RESP_COEFF 
    !Config Desc  = 
    !Config If    = OK_STOMATE 
    !Config Def   = 1.4 
    !Config Help  = 
    !Config Units = [-]
    maint_resp_coeff = 1.4 
    CALL getin_p('MAINT_RESP_COEFF',maint_resp_coeff)

    !-
    ! soilcarbon parameters 
    !-
    !Config Key   = ACTIVE_TO_PASS_CLAY_FRAC
    !Config Desc  = 
    !Config if    = OK_STOMATE 
    !Config Def   = 0.032 
    !Config Help  =
    !Config Units = [-]
    active_to_pass_clay_frac = 0.032
    CALL getin_p('ACTIVE_TO_PASS_CLAY_FRAC',active_to_pass_clay_frac)
    !
    !Config Key   = ACTIVE_TO_PASS_REF_FRAC
    !Config Desc  = Fixed fraction from Active to Passive pool
    !Config if    = OK_STOMATE 
    !Config Def   = 0.003
    !Config Help  =
    !Config Units = [-]
    active_to_pass_ref_frac = 0.003
    CALL getin_p('ACTIVE_TO_PASS_REF_FRAC',active_to_pass_ref_frac)  
    !
    !Config Key   = SURF_TO_SLOW_REF_FRAC
    !Config Desc  = Fixed fraction from Surface to Slow pool
    !Config if    = OK_STOMATE 
    !Config Def   = 0.4
    !Config Help  =
    !Config Units = [-]
    surf_to_slow_ref_frac = 0.4
    CALL getin_p('SURF_TO_SLOW_REF_FRAC',surf_to_slow_ref_frac)  
    !
    !Config Key   = ACTIVE_TO_CO2_REF_FRAC
    !Config Desc  = Fixed fraction from Active pool to CO2 emission
    !Config if    = OK_STOMATE 
    !Config Def   = 0.85
    !Config Help  =
    !Config Units = [-]
    active_to_CO2_ref_frac = 0.85 
    CALL getin_p('ACTIVE_TO_CO2_REF_FRAC',active_to_co2_ref_frac)  
    !
    !Config Key   = SLOW_TO_PASS_REF_FRAC
    !Config Desc  = Fixed fraction from Slow to Passive pool
    !Config if    = OK_STOMATE 
    !Config Def   = 0.003
    !Config Help  =
    !Config Units = [-]
    slow_to_pass_ref_frac = 0.03
    CALL getin_p('SLOW_TO_PASS_REF_FRAC',slow_to_pass_ref_frac)  
    !
    !Config Key   = SLOW_TO_CO2_REF_FRAC
    !Config Desc  = Fixed fraction from Slow pool to CO2 emission
    !Config if    = OK_STOMATE 
    !Config Def   = 0.55
    !Config Help  =
    !Config Units = [-]
    slow_to_CO2_ref_frac = 0.55
    CALL getin_p('SLOW_TO_CO2_REF_FRAC',slow_to_co2_ref_frac)  
    !
    !Config Key   = PASS_TO_ACTIVE_REF_FRAC
    !Config Desc  = Fixed fraction from Passive to Active pool
    !Config if    = OK_STOMATE 
    !Config Def   = 0.45
    !Config Help  =
    !Config Units = [-]
    pass_to_active_ref_frac = 0.45
    CALL getin_p('PASS_TO_ACTIVE_REF_FRAC',pass_to_active_ref_frac)  
    !
    !Config Key   = PASS_TO_SLOW_REF_FRAC
    !Config Desc  = Fixed fraction from Passive to Slow pool
    !Config if    = OK_STOMATE 
    !Config Def   = 0.
    !Config Help  =
    !Config Units = [-]
    pass_to_slow_ref_frac = 0.0
    CALL getin_p('PASS_TO_SLOW_REF_FRAC',pass_to_slow_ref_frac)  
    !
    !Config Key   = ACTIVE_TO_CO2_CLAY_SILT_FRAC
    !Config Desc  = Clay-Silt-dependant fraction from Active pool to CO2 emission
    !Config if    = OK_STOMATE 
    !Config Def   = 0.68
    !Config Help  =
    !Config Units = [-]
    active_to_CO2_clay_silt_frac = 0.68
    CALL getin_p('ACTIVE_TO_CO2_CLAY_SILT_FRAC',active_to_co2_clay_silt_frac)  
    !
    !Config Key   = SLOW_TO_PASS_CLAY_FRAC
    !Config Desc  = Clay-dependant fraction from Slow to Passive pool
    !Config if    = OK_STOMATE 
    !Config Def   = -0.009
    !Config Help  =
    !Config Units = [-]
    slow_to_pass_clay_frac = 0
    CALL getin_p('SLOW_TO_PASS_CLAY_FRAC',slow_to_pass_clay_frac)  
    !
    !Config Key   = SOM_TURN_ISURFACE
    !Config Desc  = turnover in surface pool
    !Config if    = OK_STOMATE 
    !Config Def   = 6.0
    !Config Help  =
    !Config Units =  [year-1]
    som_turn_isurface = 6.0 
    CALL getin_p('SOM_TURN_ISURFACE',som_turn_isurface)
    !
    !Config Key   = SOM_TURN_IACTIVE
    !Config Desc  = turnover in active pool
    !Config if    = OK_STOMATE 
    !Config Def   = 7.3
    !Config Help  =
    !Config Units =  [year-1]
    som_turn_iactive = 7.3
    CALL getin_p('SOM_TURN_IACTIVE',som_turn_iactive)
    !
    !Config Key   = SOM_TURN_ISLOW
    !Config Desc  = turnover in slow pool
    !Config if    = OK_STOMATE 
    !Config Def   = 0.2
    !Config Help  =
    !Config Units = [year-1]
    som_turn_islow = 0.1 
    CALL getin_p('SOM_TURN_ISLOW',som_turn_islow)
    !
    !Config Key   = FSLOW
    !Config Desc  = converting factor from active to slow pool turnover
    !Config if    = OK_STOMATE and OK_SOIL_CARBON_DISCRETIZATION
    !Config Def   = 37
    !Config Help  = 
    !Config Units = [-]
    fslow = 37.0 
    CALL getin_p('FSLOW',fslow)
    !
    !Config Key   = FPASSIVE
    !Config Desc  = converting factor from active to slow pool turnover
    !Config if    = OK_STOMATE and OK_SOIL_CARBON_DISCRETIZATION
    !Config Def   = 1617.45 
    !Config Help  = 
    !Config Units = [-]
    IF (ok_soil_carbon_discretization) THEN
       fpassive = 6469.8
    ELSE
       fpassive = 1617.45
    ENDIF
    CALL getin_p('FPASSIVE',fpassive)
    !
    !Config Key   = STOMATE_TAU
    !Config Desc  = turnover of the active pool 
    !Config if    = OK_STOMATE and OK_SOIL_CARBON_DISCRETIZATION
    !Config Def   = 4.699E6
    !Config Help  = 
    !Config Units = [seconds]
    stomate_tau = 4.699E6
    CALL getin_p('STOMATE_TAU',stomate_tau)
    !
    !Config Key   = DEPTH_MODIFIER
    !Config Desc  = turnover rate modifier depending on depth
    !Config if    = OK_STOMATE and OK_SOIL_CARBON_DISCRETIZATION
    !Config Def   = 1.E6
    !Config Help  = e-folding depth of turnover rates,following Koven et al.,2013,Biogeosciences. A very large value means no depth modification
    !Config Units = [-]
    depth_modifier = 1.E6
    CALL getin_p('DEPTH_MODIFIER',depth_modifier)
    !
    !Config Key   = SOM_TURN_IACTIVE_CLAY_FRAC
    !Config Desc  = clay-dependant parameter impacting on turnover rate of active pool - Tm parameter of Parton et al. 1993 (-)
    !Config if    = OK_STOMATE 
    !Config Def   = 0.75
    !Config Help  = 
    !Config Units = [-]
    som_turn_iactive_clay_frac = 0.75
    CALL getin_p('SOM_TURN_IACTIVE_CLAY_FRAC',som_turn_iactive_clay_frac)
    !
    !Config Key   = SOM_INIT_ACTIVE
    !Config Desc  = Initial active SOM carbon 
    !Config if    = OK_STOMATE
    !Config Def   = 1000
    !Config Help  = Putting some carbon in the soil speeds-up the spinup. As 
    !               this carbon comes with some soil nitrogen it also helps
    !               to grow vegetation when starting from scratch
    !Config Units = [g m-2]
    som_init_active = 1000
    CALL getin_p('SOM_INIT_ACTIVE',som_init_active)
    !
    !Config Key   = SOM_INIT_SLOW
    !Config Desc  = Initial slow SOM carbon
    !Config if    = OK_STOMATE
    !Config Def   = 3000
    !Config Help  = Putting some carbon in the soil speeds-up the spinup. As 
    !               this carbon comes with some soil nitrogen it also helps
    !               to grow vegetation when starting from scratch
    !Config Units = [g m-2]
    som_init_slow = 3000
    CALL getin_p('SOM_INIT_SLOW',som_init_slow)
    !
    !Config Key   = SOM_INIT_PASSIVE
    !Config Desc  = Initial passive SOM carbon
    !Config if    = OK_STOMATE
    !Config Def   = 3000
    !Config Help  = Putting some carbon in the soil speeds-up the spinup. As 
    !               this carbon comes with some soil nitrogen it also helps
    !               to grow vegetation when starting from scratch
    !Config Units = [g m-2]
    som_init_passive = 5000
    CALL getin_p('SOM_INIT_PASSIVE',som_init_passive)
    !
    !Config Key   = SOM_INIT_SURFACE
    !Config Desc  = Initial surface SOM carbon
    !Config if    = OK_STOMATE
    !Config Def   = 1000
    !Config Help  =
    !Config Units = [g m-2]
    som_init_surface = 1000 
    CALL getin_p('SOM_INIT_SURFACE',som_init_surface)
    !
    !Config Key   = CN_TARGET_IACTIVE_REF
    !Config Desc  = CN target ratio of active pool for soil min N = 0
    !Config if    = OK_STOMATE 
    !Config Def   = 15.
    !Config Help  = Putting some carbon in the soil speeds-up the spinup. As 
    !               this carbon comes with some soil nitrogen it also helps
    !               to grow vegetation when starting from scratch
    !Config Units = [-]
    CN_target_iactive_ref = 15.
    CALL getin_p('CN_TARGET_IACTIVE_REF',CN_target_iactive_ref)
    !
    !Config Key   = CN_TARGET_ISLOW_REF
    !Config Desc  = CN target ratio of slow pool for soil min N = 0
    !Config if    = OK_STOMATE 
    !Config Def   = 20.
    !Config Help  = 
    !Config Units = [-]
    CN_target_islow_ref = 20.
    CALL getin_p('CN_TARGET_ISLOW_REF',CN_target_islow_ref)
    !
    !Config Key   = CN_TARGET_IPASSIVE_REF
    !Config Desc  = CN target ratio of passive pool for soil min N = 0
    !Config if    = OK_STOMATE 
    !Config Def   = 10.
    !Config Help  = 
    !Config Units = [-]
    IF (ok_soil_carbon_discretization) THEN
       CN_target_ipassive_ref = 250.
    ELSE
       CN_target_ipassive_ref = 10.
    ENDIF
    CALL getin_p('CN_TARGET_IPASSIVE_REF',CN_target_ipassive_ref)
    !
    !Config Key   = CN_TARGET_ISURFACE_REF
    !Config Desc  = CN target ratio of surface pool for soil min N = 0
    !Config if    = OK_STOMATE 
    !Config Def   = 20.
    !Config Help  = 
    !Config Units = [-]
    CN_target_isurface_ref = 20.
    CALL getin_p('CN_TARGET_ISURFACE_REF',CN_target_isurface_ref)
    !
    !Config Key   = CN_TARGET_IACTIVE_NMIN
    !Config Desc  = CN target ratio change per mineral N unit (g m-2) for active pool 
    !Config if    = OK_STOMATE 
    !Config Def   = -6.
    !Config Help  = 
    !Config Units = [(g m-2)-1]
    CN_target_iactive_Nmin = -6. 
    CALL getin_p('CN_TARGET_IACTIVE_NMIN',CN_target_iactive_Nmin)
    !
    !Config Key   = CN_TARGET_ISLOW_NMIN
    !Config Desc  = CN target ratio change per mineral N unit (g m-2) for slow pool 
    !Config if    = OK_STOMATE 
    !Config Def   = -4.
    !Config Help  = 
    !Config Units = [(g m-2)-1]
    CN_target_islow_Nmin = -4.
    CALL getin_p('CN_TARGET_ISLOW_NMIN',CN_target_islow_Nmin)
    !
    !Config Key   = CN_TARGET_IPASSIVE_NMIN
    !Config Desc  = CN target ratio change per mineral N unit (g m-2) for passive pool 
    !Config if    = OK_STOMATE 
    !Config Def   = -1.5
    !Config Help  = 
    !Config Units = [(g m-2)-1]
    IF (ok_soil_carbon_discretization) THEN
        CN_target_ipassive_Nmin = -0.5
    ELSE
        CN_target_ipassive_Nmin = -1.5
    ENDIF
    CALL getin_p('CN_TARGET_IPASSIVE_NMIN',CN_target_ipassive_Nmin)
    !
    !Config Key   = CN_TARGET_ISURFACE_PNC
    !Config Desc  = CN target ratio change per plant nitrogen content unit (%) for surface pool
    !Config if    = OK_STOMATE 
    !Config Def   = -5.
    !Config Help  = 
    !Config Units = [-]
    CN_target_isurface_pnc = -5.
    CALL getin_p('CN_TARGET_ISURFACE_PNC',CN_target_isurface_pnc)
    
    ! soil nitrogen dynamic parameters
    !-
    !Config Key   = H_SAXTON
    !Config Desc  = Coefficient h for computing soil moisture content at saturation
    !Config If    = OK_STOMATE 
    !Config Def   = 0.332
    !Config Help  = Used for soil porosity when calculating the maximum
    !Config         pore volume of the soil in to calculate the 
    !Config         volumetric fraction of aneorobic microsites (ANVF) in 
    !Config         the soil, which gives an idea of how much aneorobic 
    !Config         bacteria can be transforming soil nitrogen.
    !Config Units = [m^3/m^3]
    h_saxton = 0.332
    CALL getin_p('H_SAXTON',h_saxton)
    !-
    !Config Key   = J_SAXTON
    !Config Desc  = Coefficient j for computing soil moisture content at saturation
    !Config If    = OK_STOMATE 
    !Config Def   = -7.251*1e-4 
    !Config Help  = Used for soil porosity when calculating the maximum
    !Config         pore volume of the soil in to calculate the 
    !Config         volumetric fraction of aneorobic microsites (ANVF) in 
    !Config         the soil, which gives an idea of how much aneorobic 
    !Config         bacteria can be transforming soil nitrogen
    !Config Units = [m^3/m^3]
    j_saxton = -7.251*1e-4 
    CALL getin_p('J_SAXTON',j_saxton)
    !-
    !Config Key   = K_SAXTON
    !Config Desc  = Coefficient k for computing soil moisture content at saturation
    !Config If    = OK_STOMATE 
    !Config Def   = O.1276
    !Config Help  = Used for soil porosity when calculating the maximum
    !Config         pore volume of the soil in to calculate the 
    !Config         volumetric fraction of aneorobic microsites (ANVF) in 
    !Config         the soil, which gives an idea of how much aneorobic 
    !Config         bacteria can be transforming soil nitrogen
    !Config Units = [m^3/m^3]
    k_saxton = 0.1276 
    CALL getin_p('K_SAXTON',k_saxton)
    !-
    !Config Key   = DIFFUSIONO2_POWER_1
    !Config Desc  = Power used in the equation defining the diffusion of oxygen in soil
    !Config If    = OK_STOMATE 
    !Config Def   = 3.33
    !Config Help  = Diffusion of oxygen determines how well anerobic bacteria
    !Config         bacteria can live in the soil, which impacts the nitrogen
    !Config         dynamics.  This is taken from the literature.
    !Config Units = [-]
    diffusionO2_power_1 = 3.33
    CALL getin_p('DIFFUSIONO2_POWER_1',diffusionO2_power_1)
    !-
    !Config Key   = DIFFUSIONO2_POWER_2
    !Config Desc  = Power used in the equation defining the diffusion of oxygen in soil
    !Config If    = OK_STOMATE 
    !Config Def   = 2.0
    !Config Help  = Diffusion of oxygen determines how well anerobic bacteria
    !Config         bacteria can live in the soil, which impacts the nitrogen
    !Config         dynamics.  This is taken from the literature.
    !Config Units = [-]
    diffusionO2_power_2 = 2.0
    CALL getin_p('DIFFUSIONO2_POWER_2',diffusionO2_power_2)
    !-
    !Config Key   = F_NOFROST
    !Config Desc  = Temperature-related Factor impacting on Oxygen diffusion rate
    !Config If    = OK_STOMATE 
    !Config Def   = 1.2
    !Config Help  = Diffusion of oxygen determines how well anerobic bacteria
    !Config         bacteria can live in the soil, which impacts the nitrogen
    !Config         dynamics.  This is taken from Table 2 of Li et al, 2000
    !Config Units = [-]
    F_nofrost = 1.2 
    CALL getin_p('F_NOFROST',F_nofrost)
    !-
    !Config Key   = F_FROST
    !Config Desc  = Temperature-related Factor impacting on Oxygen diffusion rate
    !Config If    = OK_STOMATE 
    !Config Def   = 0.8
    !Config Help  = Diffusion of oxygen determines how well anerobic bacteria
    !Config         bacteria can live in the soil, which impacts the nitrogen
    !Config         dynamics.  This is taken from Table 2 of Li et al, 2000
    !Config Units = [-]
    F_frost = 0.8 
    CALL getin_p('F_FROST',F_frost)
    !-
    !Config Key   = A_ANVF
    !Config Desc  = Coefficient used in the calculation of Volumetric fraction of anaerobic microsites
    !Config If    = OK_STOMATE 
    !Config Def   = 0.85
    !Config Help  = Anerobic bacteria grow in soil microsites, which impact 
    !Config         nitrogen dynamics.  The equation using these parameters
    !Config         is from the literature, but no values are given in the
    !Config         paper.  This value is taken from a previous version of
    !Config         the code.
    !Config Units = [-]
    a_anvf = 0.85
    CALL getin_p('A_ANVF',a_anvf)
    !-
    !Config Key   = B_ANVF
    !Config Desc  = Coefficient used in the calculation of Volumetric fraction of anaerobic microsites
    !Config If    = OK_STOMATE 
    !Config Def   = 1.
    !Config Help  = Anerobic bacteria grow in soil microsites, which impact 
    !Config         nitrogen dynamics.  The equation using these parameters
    !Config         is from the literature, but no values are given in the
    !Config         paper.  This value is taken from a previous version of
    !Config         the code.
    !Config Units = [-]
    b_anvf = 1.
    CALL getin_p('B_ANVF',b_anvf)
    !-
    !Config Key   = A_FIXNH4
    !Config Desc  = Coefficient used in the calculation of the Fraction of adsorbed NH4+
    !Config If    = OK_STOMATE 
    !Config Def   = 0.41
    !Config Help  = In particular, this seems to be for the calculation of
    !Config         adsorption onto soil clays.  Taken from the literature.
    !Config Units = [-]
    a_FixNH4 = 0.41
    CALL getin_p('A_FIXNH4',a_FixNH4)
    !-
    !Config Key   = B_FIXNH4
    !Config Desc  = Coefficient used in the calculation of the Fraction of adsorbed NH4+
    !Config If    = OK_STOMATE 
    !Config Def   = -0.47
    !Config Help  = In particular, this seems to be for the calculation of
    !Config         adsorption onto soil clays.  Taken from the literature.
    !Config Units = [-]
    b_FixNH4 = -0.47
    CALL getin_p('B_FIXNH4',b_FixNH4)
    !-
    !Config Key   = CLAY_MAX
    !Config Desc  = Coefficient used in the calculation of the Fraction of adsorbed NH4+
    !Config If    = OK_STOMATE 
    !Config Def   = 0.63
    !Config Help  = In particular, this seems to be for the calculation of
    !Config         adsorption onto soil clays.  Taken from the literature.
    !Config Units = [-]
    clay_max = 0.63
    CALL getin_p('CLAY_MAX',clay_max)
    !-
    !Config Key   = FW_NIT_0
    !Config Desc  = Coefficient used in the calculation of the Response of Nitrification to soil moisture
    !Config If    = OK_STOMATE 
    !Config Def   = -0.0243
    !Config Help  = Taken from the literature.
    !Config Units = [-]
    fw_nit_0 = -0.0243
    CALL getin_p('FW_NIT_0',fw_nit_0)
    !-
    !Config Key   = FW_NIT_1
    !Config Desc  = Coefficient used in the calculation of the Response of Nitrification to soil moisture
    !Config If    = OK_STOMATE 
    !Config Def   = 0.9975
    !Config Help  = Taken from the literature.
    !Config Units = [-]
    fw_nit_1 = 0.9975
    CALL getin_p('FW_NIT_1',fw_nit_1)
    !-
    !Config Key   = FW__NIT_2
    !Config Desc  = Coefficient used in the calculation of the Response of Nitrification to soil moisture
    !Config If    = OK_STOMATE 
    !Config Def   = -5.5368
    !Config Help  = Taken from the literature.
    !Config Units = [-]
    fw_nit_2 = -5.5368
    CALL getin_p('FW_NIT_2',fw_nit_2)
    !-
    !Config Key   = FW_NIT_3
    !Config Desc  = Coefficient used in the calculation of the Response of Nitrification to soil moisture
    !Config If    = OK_STOMATE 
    !Config Def   = 17.651
    !Config Help  = Taken from the literature.
    !Config Units = [-]
    fw_nit_3 = 17.651
    CALL getin_p('FW_NIT_3',fw_nit_3)
    !-
    !Config Key   = FW_NIT_4
    !Config Desc  = Coefficient used in the calculation of the Response of Nitrification to soil moisture
    !Config If    = OK_STOMATE 
    !Config Def   = -12.904
    !Config Help  = Taken from the literature.
    !Config Units = [-]
    fw_nit_4 = -12.904 
    CALL getin_p('FW_NIT_4',fw_nit_4)
    !-
    !Config Key   = FT_NIT_0
    !Config Desc  = Coefficient used in the calculation of the Response of Nitrification to Temperature
    !Config If    = OK_STOMATE 
    !Config Def   = -0.0233
    !Config Help  =
    !Config Units = [-]
    ft_nit_0 = -0.0233
    CALL getin_p('FT_NIT_0',ft_nit_0)
    !-
    !Config Key   = FT_NIT_1
    !Config Desc  = Coefficient used in the calculation of the Response of Nitrification to Temperature
    !Config If    = OK_STOMATE 
    !Config Def   = 0.3094
    !Config Help  = Taken from the literature.  NOTE: Zhang et al 2002 fold 
    !Config         in the factor 0.1 with the parameter for the term linear 
    !Config         in temperature (ft_nit_1), while we group it with the soil
    !Config         temperature as per the rest of the terms.  Checking 
    !Config         the parameter values, this leads to a first impression 
    !Config         that they differ by a factor of 10, when in reality the 
    !Config         same overall result is calculated.
    !Config Units = [-]
    ft_nit_1 = 0.3094
    CALL getin_p('FT_NIT_1',ft_nit_1)
    !-
    !Config Key   = FT_NIT_2
    !Config Desc  = Coefficient used in the calculation of the Response of Nitrification to Temperature
    !Config If    = OK_STOMATE 
    !Config Def   = -0.2234
    !Config Help  =
    !Config Units = [-]
    ft_nit_2 = -0.2234
    CALL getin_p('FT_NIT_2',ft_nit_2)
    !-
    !Config Key   = FT_NIT_3
    !Config Desc  = Coefficient used in the calculation of the Response of Nitrification to Temperature
    !Config If    = OK_STOMATE 
    !Config Def   = 0.1566
    !Config Help  =
    !Config Units = [-]
    ft_nit_3 = 0.1566
    CALL getin_p('FT_NIT_3',ft_nit_3)
    !-
    !Config Key   = FT_NIT_4
    !Config Desc  = Coefficient used in the calculation of the Response of Nitrification to Temperature
    !Config If    = OK_STOMATE 
    !Config Def   = -0.0272
    !Config Help  =
    !Config Units = [-]
    ft_nit_4 = -0.0272
    CALL getin_p('FT_NIT_4',ft_nit_4)
    !-
    !Config Key   = FPH_0
    !Config Desc  = Coefficient used in the calculation of the Response of Nitrification to pH
    !Config If    = OK_STOMATE 
    !Config Def   = -1.2314
    !Config Help  = Taken from the literature.
    !Config Units = [-]
    fph_0 = -1.2314
    CALL getin_p('FPH_0',fph_0)
    !-
    !Config Key   = FPH_1
    !Config Desc  = Coefficient used in the calculation of the Response of Nitrification to pH
    !Config If    = OK_STOMATE 
    !Config Def   = 0.7347
    !Config Help  = Taken from the literature.
    !Config Units = [-]
    fph_1 = 0.7347
    CALL getin_p('FPH_1',fph_1)
    !-
    !Config Key   = FPH_2
    !Config Desc  = Coefficient used in the calculation of the Response of Nitrification to pH
    !Config If    = OK_STOMATE 
    !Config Def   = -0.0604
    !Config Help  = Taken from the literature.
    !Config Units = [-]
    fph_2 = -0.0604
    CALL getin_p('FPH_2',fph_2)
    !-
    !Config Key   = FTV_0
    !Config Desc  = Coefficient used in the calculation of the response of NO2 or NO production during nitrificationof to Temperature
    !Config If    = OK_STOMATE 
    !Config Def   = 2.72
    !Config Help  =
    !Config Units = [-]
    ftv_0 = 2.72 
    CALL getin_p('FTV_0',ftv_0)
    !-
    !Config Key   = FTV_1
    !Config Desc  = Coefficient used in the calculation of the response of NO2 or NO production during nitrificationof to Temperature
    !Config If    = OK_STOMATE 
    !Config Def   = 34.6
    !Config Help  =
    !Config Units = [-]
    ftv_1 = 34.6
    CALL getin_p('FTV_1',ftv_1)
    !-
    !Config Key   = FTV_2
    !Config Desc  = Coefficient used in the calculation of the response of NO2 or NO production during nitrificationof to Temperature
    !Config If    = OK_STOMATE 
    !Config Def   = 9615.
    !Config Help  =
    !Config Units = [-]
    ftv_2 = 9615.
    CALL getin_p('FTV_2',ftv_2)
    !-
    !Config Key   = K_NITRIF
    !Config Desc  = Nitrification rate at 20 ◦C and field capacity
    !Config If    = OK_STOMATE 
    !Config Def   = 2.0
    !Config Help  = The literature value is 0.2, Schmid et al., 2001
    !Config         (https://doi.org/10.1023/A:1012694218748)
    !Config         However, the value in OCN appears to be 2.0.  We keep
    !Config         the OCN value (see stomate_soilcarbon.f90 for more info).
    !Config Units = [day**-1]
    IF (ok_soil_carbon_discretization) THEN
        k_nitrif = 0.5 
    ELSE
        k_nitrif = 2.0
    ENDIF
    CALL getin_p('K_NITRIF',k_nitrif)
    !-
    !Config Key   = N2O_NITRIF_P
    !Config Desc  = Reference n2o production per N-NO3 produced g N-N2O
    !Config If    = OK_STOMATE 
    !Config Def   = 0.0006
    !Config Help  = Taken from Zhang et al., 2002 - Appendix A p. 102
    !Config Units = [gN-N2O (gN-NO3)-1]
    IF (ok_soil_carbon_discretization) THEN
       n2o_nitrif_p = 0.0001
    ELSE
       n2o_nitrif_p = 0.0006
    ENDIF
    CALL getin_p('N2O_NITRIF_P',n2o_nitrif_p)
    !-
    !Config Key   = NO_NITRIF_P
    !Config Desc  = Reference NO production per N-NO3 produced g N-N2O
    !Config If    = OK_STOMATE 
    !Config Def   = 0.0025
    !Config Help  = Taken from Zhang et al., 2002 - Appendix A p. 102
    !Config Units = [gN-NO (gN-NO3)-1]
    IF (ok_soil_carbon_discretization) THEN
        no_nitrif_p = 0.0005
    ELSE
        no_nitrif_p = 0.0025
    ENDIF    
    CALL getin_p('NO_NITRIF_P',no_nitrif_p)
    !-
    !Config Key   = CHEMO_T0
    !Config Desc  = Coefficient used in the calculation of the Response of NO production from chemodenitrification to Temperature
    !Config If    = OK_STOMATE 
    !Config Def   = -31494
    !Config Help  =
    !Config Units = [-]
    chemo_t0 = -31494.
    CALL getin_p('CHEMO_T0',chemo_t0)
    !-
    !Config Key   = CHEMO_PH0
    !Config Desc  = Coefficient used in the calculation of the Response of NO production from chemodenitrification to pH
    !Config If    = OK_STOMATE 
    !Config Def   = -1.62
    !Config Help  =
    !Config Units = [-]
    chemo_ph0 = -1.62
    CALL getin_p('CHEMO_PH0',chemo_ph0)
    !-
    !Config Key   = CHEMO_0
    !Config Desc  = Coefficient used in the calculation of NO production from chemodenitrification
    !Config If    = OK_STOMATE 
    !Config Def   = 30.
    !Config Help  =
    !Config Units = [-]
    IF (ok_soil_carbon_discretization) THEN
         chemo_0 = 5.
    ELSE
         chemo_0 = 30.
    ENDIF
    CALL getin_p('CHEMO_0',chemo_0)
    !-
    !Config Key   = CHEMO_1
    !Config Desc  = Coefficient used in the calculation of NO production from chemodenitrification
    !Config If    = OK_STOMATE 
    !Config Def   = 16565
    !Config Help  =
    !Config Units = [-]
    chemo_1 = 16565.
    CALL getin_p('CHEMO_1',chemo_1)
    !-
    !Config Key   = FT_DENIT_0
    !Config Desc  = Coefficient used in the response of relative growth rate of total denitrifiers to Temperature
    !Config If    = OK_STOMATE 
    !Config Def   = 2.
    !Config Help  =
    !Config Units = [-]
    ft_denit_0 = 2.
    CALL getin_p('FT_DENIT_0',ft_denit_0)
    !-
    !Config Key   = FT_DENIT_1
    !Config Desc  = Coefficient used in the response of relative growth rate of total denitrifiers to Temperature
    !Config If    = OK_STOMATE 
    !Config Def   = 22.5
    !Config Help  =
    !Config Units = [-]
    ft_denit_1 = 22.5
    CALL getin_p('FT_DENIT_1',ft_denit_1)
    !-
    !Config Key   = FT_DENIT_2
    !Config Desc  = Coefficient used in the response of relative growth rate of total denitrifiers to Temperature
    !Config If    = OK_STOMATE 
    !Config Def   = 10
    !Config Help  =
    !Config Units = [-]
    ft_denit_2 = 10.
    CALL getin_p('FT_DENIT_2',ft_denit_2)
    !-
    !Config Key   = FPH_NO3_0
    !Config Desc  = Coefficient used in the response of relative growth rate of NO3 denitrifiers to pH
    !Config If    = OK_STOMATE 
    !Config Def   = 4.25
    !Config Help  =
    !Config Units = [-]
    fph_no3_0 = 4.25
    CALL getin_p('FPH_NO3_0',fph_no3_0)
    !-
    !Config Key   = FPH_NO3_1
    !Config Desc  = Coefficient used in the response of relative growth rate of NO3 denitrifiers to pH
    !Config If    = OK_STOMATE 
    !Config Def   = 0.5
    !Config Help  =
    !Config Units = [-]
    fph_no3_1 = 0.5  
    CALL getin_p('FPH_NO3_1',fph_no3_1)
    !-
    !Config Key   = FPH_NO_0
    !Config Desc  = Coefficient used in the response of relative growth rate of NO denitrifiers to pH
    !Config If    = OK_STOMATE 
    !Config Def   = 5.25
    !Config Help  =
    !Config Units = [-]
    fph_no_0 = 5.25
    CALL getin_p('FPH_NO_0',fph_no_0)
    !-
    !Config Key   = FPH_NO_1
    !Config Desc  = Coefficient used in the response of relative growth rate of NO denitrifiers to pH
    !Config If    = OK_STOMATE 
    !Config Def   = 1.
    !Config Help  =
    !Config Units = [-]
    fph_no_1  = 1. 
    CALL getin_p('FPH_NO_1',fph_no_1)
    !-
    !Config Key   = FPH_N2O_0
    !Config Desc  = Coefficient used in the response of relative growth rate of N2O denitrifiers to pH
    !Config If    = OK_STOMATE 
    !Config Def   = 6.25
    !Config Help  =
    !Config Units = [-]
    fph_n2o_0 = 6.25
    CALL getin_p('FPH_N2O_0',fph_n2o_0)
    !-
    !Config Key   = FPH_N2O_1
    !Config Desc  = Coefficient used in the response of relative growth rate of N2O denitrifiers to pH
    !Config If    = OK_STOMATE 
    !Config Def   = 1.5
    !Config Help  =
    !Config Units = [-]
    fph_n2o_1 = 1.5
    CALL getin_p('FPH_N2O_1',fph_n2o_1)
    !-
    !Config Key   = KN
    !Config Desc  = Half Saturation of N oxydes 
    !Config If    = OK_STOMATE 
    !Config Def   = 0.083
    !Config Help  =
    !Config Units = [kgN/m**3]
    Kn = 0.083 
    CALL getin_p('KN',Kn)
    !-
    !Config Key   = CTE_BACT
    !Config Desc  = Denitrification activiy of bacteria     
    !Config If    = OK_STOMATE
    !Config Def   = 0.00005
    !Config Help  = 
    !Config Units = [-]
    IF (ok_soil_carbon_discretization) THEN
       cte_bact = 0.00001
    ELSE
       cte_bact = 0.00015
    ENDIF
    CALL getin_p('CTE_BACT',cte_bact)
    !-
    !Config Key   = MU_NO3_MAX
    !Config Desc  = Maximum Relative growth rate of NO3 denitrifiers
    !Config If    = OK_STOMATE 
    !Config Def   = 0.67
    !Config Help  =
    !Config Units = [hour**-1]
    mu_no3_max = 0.67
    CALL getin_p('MU_NO3_MAX',mu_no3_max)
    !-
    !Config Key   = MU_NO_MAX
    !Config Desc  = Maximum Relative growth rate of NO denitrifiers
    !Config If    = OK_STOMATE 
    !Config Def   = 0.34
    !Config Help  =
    !Config Units = [hour**-1]
    mu_no_max = 0.67
    CALL getin_p('MU_NO_MAX',mu_no_max)
    !-
    !Config Key   = MU_N2O_MAX
    !Config Desc  = Maximum Relative growth rate of N2O denitrifiers
    !Config If    = OK_STOMATE 
    !Config Def   = 0.34
    !Config Help  =
    !Config Units = [hour**-1]
    mu_n2o_max = 0.67
    CALL getin_p('MU_N2O_MAX',mu_n2o_max)
    !-
    !Config Key   = Y_NO3
    !Config Desc  = Maximum growth yield of NO3 denitrifiers on N oxydes
    !Config If    = OK_STOMATE 
    !Config Def   = 0.401
    !Config Help  =
    !Config Units = [kgC / kgN]
    Y_no3 = 0.401
    CALL getin_p('Y_NO3',Y_no3)
    !-
    !Config Key   = Y_NO
    !Config Desc  = Maximum growth yield of NO denitrifiers on N oxydes
    !Config If    = OK_STOMATE 
    !Config Def   = 0.428
    !Config Help  =
    !Config Units = [kgC / kgN]
    Y_no = 0.428 
    CALL getin_p('Y_NO',Y_no)
    !-
    !Config Key   = Y_N2O
    !Config Desc  = Maximum growth yield of N2O denitrifiers on N oxydes
    !Config If    = OK_STOMATE 
    !Config Def   = 0.151
    !Config Help  =
    !Config Units = [kgC / kgN]
    Y_n2o = 0.151
    CALL getin_p('Y_N2O',Y_n2O)
    !-
    !Config Key   = M_NO3
    !Config Desc  = Maintenance coefficient on NO3
    !Config If    = OK_STOMATE 
    !Config Def   = 0.09
    !Config Help  =
    !Config Units = [kgN / kgC / hour]
    M_no3 = 0.09
    CALL getin_p('M_NO3',M_no3)
    !-
    !Config Key   = M_NO
    !Config Desc  = Maintenance coefficient on NO
    !Config If    = OK_STOMATE 
    !Config Def   = 0.035
    !Config Help  =
    !Config Units = [kgN / kgC / hour]
    M_no = 0.035
    CALL getin_p('M_NO',M_no)
    !-
    !Config Key   = M_N2O
    !Config Desc  = Maintenance coefficient on N2O
    !Config If    = OK_STOMATE 
    !Config Def   = 0.079
    !Config Help  =
    !Config Units = [kgN / kgC / hour]
    M_n2o = 0.079
    CALL getin_p('M_N2O',M_n2o)
    !-
    !Config Key   = MAINT_C
    !Config Desc  = Maintenance coefficient of carbon
    !Config If    = OK_STOMATE 
    !Config Def   = 0.0076
    !Config Help  =
    !Config Units = [kgC / kgC / hour]
    Maint_c = 0.0076
    CALL getin_p('MAINT_C',Maint_c)
    !-
    !Config Key   = YC
    !Config Desc  = Maximum growth yield on soluble carbon
    !Config If    = OK_STOMATE 
    !Config Def   = 0.503
    !Config Help  =
    !Config Units = [kgC / kgC ]
    Yc = 0.503 
    CALL getin_p('YC',Yc)
    !-
    !Config Key   = F_CLAY_0
    !Config Desc  = Coefficient used in the eq. defining the response of N-emission to clay fraction
    !Config If    = OK_STOMATE 
    !Config Def   = 0.13
    !Config Help  =
    !Config Units = [-]
    F_clay_0 = 0.13
    CALL getin_p('F_CLAY_0',F_clay_0)
    !-
    !Config Key   = F_CLAY_1
    !Config Desc  = Coefficient used in the eq. defining the response of N-emission to clay fraction
    !Config If    = OK_STOMATE 
    !Config Def   = -0.079
    !Config Help  =
    !Config Units = [-]
    F_clay_1 = -0.079
    CALL getin_p('F_CLAY_1',F_clay_1)
    !-
    !Config Key   = RATIO_NH4_FERT
    !Config Desc  = Proportion of ammonium in the fertilizers (ammo-nitrate) 
    !Config If    = OK_STOMATE 
    !Config Def   = 0.875
    !Config Help  =
    !Config Units = [-]
    ratio_nh4_fert = 0.5
    CALL getin_p('RATIO_NH4_FERT',ratio_nh4_fert)
    !-
    !Config Key   = BOOST_FERT
    !Config Desc  = use for fertilisation scenario while using inputs maps 
    !Config If    = OK_STOMATE 
    !Config Def   = 0.0
    !Config Help  =
    !Config Units = [-]
    boost_fert = 0.0
    CALL getin_p('BOOST_FERT',boost_fert)

    !-
    !Config Key   = MANURE_TYPE
    !Config Desc  = C:N ratio of organic fertilizers coming from Fuchs,et al,
    !Effets agronomiques attendus de l’épandage des Mafor sur les écosystèmes
    !agricoles et forestiers, Valoris. des matières Fertil. d’origine résiduaire
    !sur les sols à usage Agric. ou For., 364–567 [online] manure type table 3-1-1 Fuchs et al. 2014
    !Config If    = OK_STOMATE 
    !Config Def   = CASE (e.g dairy_cow, feedlots)
    !Config Help  =
    !Config Units = [-]
    manure_type = 'average'
    CALL getin_p('MANURE_TYPE',manure_type)

    SELECT CASE (manure_type)
    CASE ('average','none','forced')
        cn_ratio_manure = 13.7           !! C:N ratio of organic fertilizer (average over table 3-1-1 from Fuchs et al. 2014) 
        CALL getin_p('CN_RATIO_MANURE',cn_ratio_manure)
    CASE ('cow')
        cn_ratio_manure = 10.28          !! C:N ratio of organic fertilizer of cow manure (table 3-1-1 from Fuchs et al. 2014)
    CASE ('dairy')
        cn_ratio_manure = 18.93          !! C:N ratio of organic fertilizer of dairy cow manure (table 3-1-1 from Fuchs et al. 2014)
    CASE ('dairyres')
        cn_ratio_manure = 24.63          !! C:N ratio of organic fertilizer of dairy cow + residues (table 3-1-1 from Fuchs et al. 2014)
    CASE ('feedlot')
        cn_ratio_manure = 9.57           !! C:N ratio of organic fertilizer of cow feedlots (table 3-1-1 from Fuchs et al. 2014)
    CASE ('comp')
        cn_ratio_manure = 11.99          !! C:N ratio of organic fertilizer of composted manure cow (table 3-1-1 from Fuchs et al. 2014)
    CASE ('pig')
        cn_ratio_manure = 4.08           !! C:N ratio of organic fertilizer of pig (table 3-1-1 from Fuchs et al. 2014)
    CASE ('pigsolid')
        cn_ratio_manure = 11.65          !! C:N ratio of  pigsolid organic fertilizer (table 3-1-1 from Fuchs et al. 2014)
    CASE ('goat')
        cn_ratio_manure = 11.36          !! C:N ratio of organic fertilizer of goat manure (table 3-1-1 from Fuchs et al. 2014)
    CASE ('poultry')
        cn_ratio_manure = 6.06           !! C:N ratio of organic fertilizer of poultry (table 3-1-1 from Fuchs et al. 2014)
    CASE DEFAULT
        WRITE(numout,*) "Choose a manure  type"
        STOP 'constantes.f90'
    ENDSELECT
    IF (printlev>=2) WRITE(numout,*)'manure type,cn_ratio_manure', manure_type, cn_ratio_manure

    !-
    !-
    ! Arrays
    !-
    !Config Key   = K_N_MIN
    !Config Desc  = [NH4+] and [NO3-] for which the Nuptake equals vmax/2. 
    !Config If    = OK_STOMATE 
    !Config Def   = 30. 30.
    !Config Help  =
    !Config Units = [umol per litter]
    K_N_min = (/ 30., 30. /)
    CALL getin_p('K_N_min',K_N_min)
    !-
    !Config Key   = LOW_K_N_MIN
    !Config Desc  = Rate of N uptake not associated with Michaelis- Menten Kinetics for Ammonium
    !Config If    = OK_STOMATE 
    !Config Def   = 0.0002 0.0002
    !Config Help  =
    !Config Units = [umol**-1]
    low_K_N_min = (/ 0.0002, 0.0002 /)
    CALL getin_p('LOW_K_N_min',low_K_N_min)

    !Config Key   = EMM_FAC
    !Config Desc  = Factor for reducing NH3 emission  
    !Config If    = OK_NCYCLE
    !Config Def   = 0.2
    !Config Help  =
    !Config Units = [-]
    IF (ok_soil_carbon_discretization) THEN
        emm_fac = 0.05
    ELSE
        emm_fac = 0.2
    ENDIF        
    CALL getin_p('EMM_FAC',emm_fac)

    !Config Key   = FACT_KN_NO
    !Config Desc  = Factor for adusting kn constant for NOx production
    !Config If    = OK_NCYCLE
    !Config Def   = 0.012
    !Config Help  =
    !Config Units = [-]
    fact_kn_no = 0.012 
    CALL getin_p('FACT_KN_NO',fact_kn_no)

    !Config Key   = FACT_KN_N2O
    !Config Desc  = Factor for adusting kn constant for N2O production
    !Config If    = OK_NCYCLE
    !Config Def   = 0.04
    !Config Help  =
    !Config Units = [-]
    fact_kn_n2o = 0.04 
    CALL getin_p('FACT_KN_N2O',fact_kn_n2o)

    !Config Key   = KFWDENIT
    !Config Desc  = Factor for adjusting sensitivity of denitrification to water content
    !Config If    = OK_NCYCLE
    !Config Def   = -5.
    !Config Help  =
    !Config Units = [-]
    kfwdenit = -5. 
    CALL getin_p('KFWDENIT',kfwdenit)

    !Config Key   = FWDENITFC
    !Config Desc  = Value at field capacity of the sensitivity function of denitrification to water content
    !Config If    = OK_NCYCLE
    !Config Def   = 0.05
    !Config Help  =
    !Config Units = [-]
    fwdenitfc = 0.05 
    CALL getin_p('FWDENITFC',fwdenitfc)

    !Config Key   = FRACN_DRAINAGE
    !Config Desc  = Fraction of NH3/NO3 loss by drainage 
    !Config If    = OK_NCYCLE
    !Config Def   = 1.0
    !Config Help  =
    !Config Units = [-]
    fracn_drainage = 1.0 
    CALL getin_p('FRACN_DRAINAGE',fracn_drainage)

    !Config Key   = FRACN_RUNOFF
    !Config Desc  = Fraction of NH3/NO3 loss by runoff
    !Config If    = OK_NCYCLE
    !Config Def   = 0.3
    !Config Help  =
    !Config Units = [-]
    fracn_runoff = 1.0
    CALL getin_p('FRACN_RUNOFF',fracn_runoff)
    !-
    !Config Key   = LEAF_N_DMAX
    !Config Desc  = Maximal elasticity of foliage N concentrations
    !Config If    = OK_STOMATE 
    !Config Def   = 0.25
    !Config Help  = Maximal elasticity of foliage N concentrations (Zaehle et al 2010, SI2, Table 1)
    !Config Units = ???
    Dmax = 0.25 
    CALL getin_p('LEAF_N_DMAX',DMAX)
    !-
    !Config Key   = P_N_UPTAKE
    !Config Desc  = Minimum value of the correction factor for plant N uptake
    !               case of enough N in the reserve
    !Config If    = OK_STOMATE 
    !Config Def   = 0.6
    !Config Help  =
    !Config Units = [-]
    p_n_uptake = 0.6  
    CALL getin_p('P_N_UPTAKE',p_n_uptake)
    
    
    !-
    ! growth_fun_all
    !
    !Config Key   = SYNC_THRESHOLD
    !Config Desc  = The threshold value for a warning when we sync biomass
    !Config If    = OK_STOMATE 
    !Config Def   = 0.1 
    !Config Help  = 
    !Config Units = [-]
    sync_threshold = 0.0001
    CALL getin_p('SYNC_THRESHOLD',sync_threshold)
    !
    !Config Key   = MAX_DELTA_KF
    !Config Desc  = Maximum change in KF from one time step to another
    !Config If    = OK_STOMATE 
    !Config Def   = 0.1 
    !Config Help  = Maximum change in KF from one time step to another (m)
    !Config Help  = This is a bit arbitrary.
    !Config Units = [m]
    max_delta_KF = 0.1
    CALL getin_p('MAX_DELTA_KF',max_delta_KF)
    !
    !Config Key   = MAINT_FROM_GPP
    !Config Desc  = Some carbon needs to remain to support the growth, hence, 
    !               respiration will be limited. In this case resp_maint 
    !               (gC m-2 dt-1) should not be more than 80% (::maint_from_gpp) 
    !               of the GPP (gC m-2 s-1)
    !Config If    = OK_STOMATE 
    !Config Def   = 0.8
    !Config Help  = Some carbon needs to remain to support the growth, hence, 
    !Config Help  = respiration will be limited. In this case resp_maint 
    !Config Help  = (gC m-2 dt-1) should not be more than 80% (::maint_from_gpp) 
    !Config Help  = of the GPP (gC m-2 s-1)
    !Config Units = [-]
    maint_from_gpp = 0.8 
    CALL getin_p('MAINT_FROM_GPP',maint_from_gpp)


    !-
    ! turnover parameters
    !-
    !
    !Config Key   = NEW_TURNOVER_TIME_REF
    !Config Desc  = 
    !Config If    = OK_STOMATE 
    !Config Def   = 20. 
    !Config Help  = 
    !Config Units = [days]
    new_turnover_time_ref = 20.
    CALL getin_p('NEW_TURNOVER_TIME_REF',new_turnover_time_ref)

    !-
    ! vmax parameters
    !-
    !
    !Config Key   = VMAX_OFFSET 
    !Config Desc  = offset (minimum relative vcmax)
    !Config If    = OK_STOMATE 
    !Config Def   = 0.3
    !Config Help  = offset (minimum vcmax/vmax_opt)
    !Config Units = [-]
    vmax_offset = 0.3
    CALL getin_p('VMAX_OFFSET',vmax_offset)
    !
    !Config Key   = LEAFAGE_FIRSTMAX
    !Config Desc  = leaf age at which vmax attains vcmax_opt (in fraction of critical leaf age)
    !Config If    = OK_STOMATE 
    !Config Def   = 0.03 
    !Config Help  = relative leaf age at which vmax attains vcmax_opt
    !Config Units = [-]
    leafage_firstmax = 0.03
    CALL getin_p('LEAFAGE_FIRSTMAX',leafage_firstmax)
    !
    !Config Key   = LEAFAGE_LASTMAX 
    !Config Desc  = leaf age at which vmax falls below vcmax_opt (in fraction of critical leaf age) 
    !Config If    = OK_STOMATE 
    !Config Def   = 0.5 
    !Config Help  = relative leaf age at which vmax falls below vcmax_opt
    !Config Units = [-]
    leafage_lastmax = 0.5
    CALL getin_p('LEAFAGE_LASTMAX',leafage_lastmax)
    !
    !Config Key   = LEAFAGE_OLD 
    !Config Desc  = leaf age at which vmax attains its minimum (in fraction of critical leaf age)
    !Config If    = OK_STOMATE 
    !Config Def   = 1.
    !Config Help  = relative leaf age at which vmax attains its minimum
    !Config Units = [-]
    leafage_old = 1.
    CALL getin_p('LEAFAGE_OLD',leafage_old)
    !
    !-
    ! season parameters
    !-
    !
    !Config Key   = GPPFRAC_DORMANCE 
    !Config Desc  = rapport maximal GPP/GGP_max pour dormance
    !Config If    = OK_STOMATE 
    !Config Def   = 0.2 
    !Config Help  = 
    !Config Units = [-]
    gppfrac_dormance = 0.2
    CALL getin_p('GPPFRAC_DORMANCE',gppfrac_dormance)
    !
    !Config Key   = TAU_CLIMATOLOGY
    !Config Desc  = tau for "climatologic variables 
    !Config If    = OK_STOMATE 
    !Config Def   = 20
    !Config Help  =
    !Config Units = [days]
    tau_climatology = 20.
    CALL getin_p('TAU_CLIMATOLOGY',tau_climatology)
    !
    !Config Key   = HVC1 
    !Config Desc  = parameters for herbivore activity
    !Config If    = OK_STOMATE 
    !Config Def   = 0.019
    !Config Help  = 
    !Config Units = [-]
    hvc1 = 0.019
    CALL getin_p('HVC1',hvc1)
    !
    !Config Key   = HVC2 
    !Config Desc  = parameters for herbivore activity 
    !Config If    = OK_STOMATE 
    !Config Def   = 1.38
    !Config Help  = 
    !Config Units = [-]
    hvc2 = 1.38
    CALL getin_p('HVC2',hvc2)
    !
    !Config Key   = LEAF_FRAC_HVC
    !Config Desc  = parameters for herbivore activity 
    !Config If    = OK_STOMATE 
    !Config Def   = 0.33
    !Config Help  = 
    !Config Units = [-]
    leaf_frac_hvc = 0.33 
    CALL getin_p('LEAF_FRAC_HVC',leaf_frac_hvc)
    !
    !Config Key   = TLONG_REF_MAX
    !Config Desc  = maximum reference long term temperature 
    !Config If    = OK_STOMATE 
    !Config Def   = 303.1
    !Config Help  = 
    !Config Units = [K]
    tlong_ref_max = 303.1
    CALL getin_p('TLONG_REF_MAX',tlong_ref_max)
    !
    !Config Key   = TLONG_REF_MIN 
    !Config Desc  = minimum reference long term temperature 
    !Config If    = OK_STOMATE 
    !Config Def   = 253.1
    !Config Help  = 
    !Config Units = [K]
    tlong_ref_min = 253.1
    CALL getin_p('TLONG_REF_MIN',tlong_ref_min)
    !
    !Config Key   = NCD_MAX_YEAR
    !Config Desc  = 
    !Config If    = OK_STOMATE 
    !Config Def   = 3. 
    !Config Help  = NCD : Number of Chilling Days
    !Config Units = [days]
    ncd_max_year = 3.
    CALL getin_p('NCD_MAX_YEAR',ncd_max_year)
    !
    !Config Key   = NGD_THRESHOLD 
    !Config Desc  = 
    !Config If    = OK_STOMATE 
    !Config Def   = -5. 
    !Config Help  = NGD : Growing Degree-Day since midwinter
    !Config Units = [days]
    ngd_threshold = -5.
    CALL getin_p('NGD_THRESHOLD',ngd_threshold)
    !
    !Config Key   = GREEN_AGE_EVER 
    !Config Desc  = 
    !Config If    = OK_STOMATE 
    !Config Def   = 2. 
    !Config Help  = 
    !Config Units = [-]
    green_age_ever = 2.
    CALL getin_p('GREEN_AGE_EVER',green_age_ever)
    !
    !Config Key   = GREEN_AGE_DEC
    !Config Desc  = 
    !Config If    = OK_STOMATE 
    !Config Def   = 0.5 
    !Config Help  = 
    !Config Units = [-]
    green_age_dec = 0.5
    CALL getin_p('GREEN_AGE_DEC',green_age_dec)

    !Config Key   = NGD_MIN_DORMANCE
    !Config Desc  = Minimum length (days) of the dormance period for species with the ngd phenology type
    !Config If    = OK_STOMATE 
    !Config Def   = 55.
    !Config Help  = 
    !Config Units = [days]
    ngd_min_dormance = 55.
    CALL getin_p('NGD_MIN_DORMANCE',ngd_min_dormance)
 
    
    !Config Key   = NAGEC
    !Config Desc  = Number of age classes
    !Config If    = OK_STOMATE 
    !Config Def   = 1 
    !Config Help  = Number of age classes in forestry and lcchange
    !               age classes could be considered age classes across stands 
    !               in the same pixel. They help to describe landscape heterogeneity
    !               they are most useful when land cover change is used.
    !Config Units = [-]
    nagec = 1
    CALL getin_p('NAGEC',nagec)
    !
    ! Guillaume M. -- AGE_CLASS_AGE_FRAC: age of each class as a fraction of the target
    ! rotation. /!\ MUST stay AFTER the NAGEC read. Earlier in the file nagec still holds 1,
    ! so the array is sized to 1, the ]0,1] guard only tests element 1 and classes 2..nagec
    ! are read OUT OF BOUNDS. Same ordering trap as LARGEST_TREE_DIA.
    ALLOCATE(age_class_age_frac(MAX(nagec,1)),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for age_class_age_frac. We stop.',nagec
       STOP 'constantes.f90'
    ENDIF
    IF (nagec == 4) THEN
       ! Guillaume M. -- SATURATING ladder, not the mid-point of the diameter bin: diameter
       ! saturates with age (Chapman-Richards), so at median diameter a stand is markedly
       ! OLDER than half its rotation. Mid-bin values (0.125/0.375/0.625/0.875) would
       ! underestimate the age of the intermediate classes.
       age_class_age_frac(1) = 0.10
       age_class_age_frac(2) = 0.45
       age_class_age_frac(3) = 0.70
       age_class_age_frac(4) = 0.92
    ELSE
       DO ic_loc = 1, MAX(nagec,1)
          age_class_age_frac(ic_loc) = (REAL(ic_loc,r_std) - 0.5_r_std) / REAL(MAX(nagec,1),r_std)
       ENDDO
    ENDIF
    !Config Key   = AGE_CLASS_AGE_FRAC
    !Config Desc  = Stand age of each age class as a fraction of the target rotation age
    !Config If    = OK_STOMATE .AND. NAGEC > 1
    !Config Def   = 0.10, 0.45, 0.70, 0.92
    !Config Help  = age_stand_bm(class) = AGE_CLASS_AGE_FRAC(class) * target_rotation_age.
    !Config         Structured cold start: the age is derived from the rotation the
    !Config         pixel is managed for, so it stays
    !Config         consistent with the management-intensity map. Saturating ladder rather
    !Config         than mid-diameter-bin, because diameter saturates with age.
    !Config Units = [-]
    CALL getin_p('AGE_CLASS_AGE_FRAC', age_class_age_frac)
    IF (ANY(age_class_age_frac(:) <= zero) .OR. ANY(age_class_age_frac(:) > un)) THEN
       WRITE(numout,*) 'AGE_CLASS_AGE_FRAC out of ]0,1] : ',age_class_age_frac(:)
       CALL ipslerr_p(3,'constantes.f90','AGE_CLASS_AGE_FRAC must lie in ]0,1]',&
            'it is a fraction of the target rotation age','')
    ENDIF
    !Config Key   = STAND_AGE_ROTATION_MAX
    !Config Desc  = Plafond de la rotation utilisee pour DERIVER l'age au demarrage a froid
    !Config If    = OK_STOMATE .AND. NAGEC > 1
    !Config Def   = 200.
    !Config Help  = OPTION SCIENTIFIQUE (permanente). age_stand_bm(c) =
    !Config         AGE_CLASS_AGE_FRAC(c) * MIN(target_rotation_age, ce plafond).
    !Config         La carte d'intensite de gestion porte une sentinelle "non gere" a 999
    !Config         ans ; sans plafond, 16,3 % de la surface forestiere recevait des ages de
    !Config         100/450/699/919 ans (mesure spin51v2chk a froid, 2026-07-27), ce qui
    !Config         rendait l'age NON MONOTONE entre classes (classe 3 a 197,7 ans contre
    !Config         classe 4 a 194,0) alors que le convoyeur suppose l'age croissant avec
    !Config         l'indice de classe. Un peuplement non gere n'a pas de rotation : 0,92*999
    !Config         n'est pas un age. Le plafond ne touche QUE la derivation d'age ;
    !Config         target_rotation_age garde 999 pour la cadence de recolte.
    !Config Units = [years]
    stand_age_rotation_max = 200.
    CALL getin_p('STAND_AGE_ROTATION_MAX', stand_age_rotation_max)
    IF (stand_age_rotation_max <= zero) THEN
       WRITE(numout,*) 'STAND_AGE_ROTATION_MAX = ',stand_age_rotation_max
       CALL ipslerr_p(3,'constantes.f90','STAND_AGE_ROTATION_MAX doit etre > 0', &
            'c est un plafond en annees sur la rotation servant a deriver l age','')
    ENDIF
    !
    ! Guillaume M. -- Size guard: makes the ordering trap impossible to reintroduce silently.
    IF (SIZE(age_class_age_frac) /= MAX(nagec,1)) THEN
       WRITE(numout,*) 'SIZE(age_class_age_frac), nagec = ',SIZE(age_class_age_frac),nagec
       CALL ipslerr_p(3,'constantes.f90', &
            'age_class_age_frac alloue avant la lecture de NAGEC', &
            'deplacer son ALLOCATE apres CALL getin_p(NAGEC)','')
    ENDIF
    !

    ! Guillaume M. -- AGE_CLASS_BOUNDS_PFT: age_class_bound is now (nagec, nvm). The global
    ! namelist vector is broadcast unchanged to every column, so the OFF path is strictly
    ! identical. The per-PFT derivation happens later, in
    ! pft_parameters:derive_age_class_bounds, because LARGEST_TREE_DIA is not read yet here.
    ALLOCATE(age_class_bound(nagec,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for age_class_bound. We stop. We need nagec x nvm words',nagec,nvm
       STOP 'constantes.f90'
    ENDIF
    ALLOCATE(age_class_bound_glo(nagec),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for age_class_bound_glo. We stop. We need nagec words',nagec
       STOP 'constantes.f90'
    ENDIF
    !
    !Config Key   = AGE_CLASS_BOUND
    !Config Desc  = Boundaries of the age classes
    !Config If    = OK_STOMATE
    !Config Def   = 5.0
    !Config Help  = The number of age class bounds should be identical
    !               to NAGEC. Two sets of default values are provided.
    !Config Units = [m]
    IF(nagec == 4)THEN
       age_class_bound_glo(1) = 0.07
       age_class_bound_glo(2) = 0.20
       age_class_bound_glo(3) = 0.40
       age_class_bound_glo(4) = 5.00
    ELSEIF(nagec == 3)THEN
       age_class_bound_glo(1) = 0.07
       age_class_bound_glo(2) = 0.20
       age_class_bound_glo(3) = 5.00
    ELSEIF(nagec == 2)THEN
       age_class_bound_glo(1) = 0.09
       age_class_bound_glo(2) = 5.00
    ELSEIF(nagec == 1)THEN
       age_class_bound_glo(1) = 5.00
    ELSE
       age_class_bound_glo(:) = -9999.
    ENDIF
    CALL getin_p('AGE_CLASS_BOUND',age_class_bound_glo)

    IF (age_class_bound_glo(1) == -9999.) THEN
       WRITE(numout,*) 'The code does not contain default values for age_class_bound'
       WRITE(numout,*) 'for this number of age classes (nagec) ',nagec
       CALL ipslerr_p (3,'constantes.f90', 'Define age class bounds',&
            'add default values in constantes.f90',&
            'or add values in the orchidee_pft.def')
    END IF

    ! Guillaume M. -- Broadcast the global vector to every PFT (default state = historical behaviour).
    DO jv_loc = 1, nvm
       age_class_bound(:,jv_loc) = age_class_bound_glo(:)
    ENDDO
    DEALLOCATE(age_class_bound_glo)

    !
    !Config Key   = MAT_F_REF
    !Config Desc  = Reference promotion rate per age class under OK_MATURITY_TRANSFER
    !Config If    = OK_MATURITY_TRANSFER
    !Config Def   = 0.065, 0.057, 0.035, 0.0 (nagec=4) ; 0.05 otherwise
    !Config Help  = Fraction of the donor area promoted to the next maturity class per
    !Config         year when the eligibility threshold is met. The nagec=4 defaults are
    !Config         the inverse residence times MEASURED on a converged reference run
    !Config         (Little's law over the last 50 years), not tuned values. The terminal
    !Config         entry is never read: that class leaves through the harvest.
    !Config Units = [1/year]
    ! Guillaume M. -- Allocated HERE and not with the module constants: mat_f_ref is
    ! dimensioned by nagec, which getin_p only fixes a few lines above. An earlier
    ! allocation would silently size it to the compile-time default.
    ALLOCATE(mat_f_ref(nagec),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for mat_f_ref. We stop. We need nagec words',nagec
       STOP 'constantes.f90'
    ENDIF
    IF (nagec == 4) THEN
       mat_f_ref(1) = 0.065_r_std
       mat_f_ref(2) = 0.057_r_std
       mat_f_ref(3) = 0.035_r_std
       mat_f_ref(4) = zero
    ELSE
       mat_f_ref(:) = 0.05_r_std
       mat_f_ref(nagec) = zero
    ENDIF
    CALL getin_p('MAT_F_REF', mat_f_ref)
    IF (ok_maturity_transfer .AND. &
         (MINVAL(mat_f_ref) < zero .OR. MAXVAL(mat_f_ref) > un)) THEN
       CALL ipslerr_p(3,'constantes','MAT_F_REF out of [0,1]', &
            'a promotion rate is a fraction of the donor area per year','')
    ENDIF

    !
    !Config Key   = MAT_TARGET_FRAC
    !Config Desc  = Target area share of each age class within its group (regulator setpoint)
    !Config If    = OK_MATURITY_TRANSFER
    !Config Def   = 1/nagec for every class
    !Config Help  = The regulator measures each class filling ratio against
    !Config         MAT_TARGET_FRAC * A_group. A uniform target assumes equal residence
    !Config         times, which diameter growth does not guarantee; a target derived from
    !Config         the class diameter bounds would be sounder -- documented as open in the
    !Config         design. Must sum to one.
    !Config Units = [-]
    ALLOCATE(mat_target_frac(nagec),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for mat_target_frac. We stop. We need nagec words',nagec
       STOP 'constantes.f90'
    ENDIF
    mat_target_frac(:) = un / REAL(nagec,r_std)
    CALL getin_p('MAT_TARGET_FRAC', mat_target_frac)
    ! Guillaume M. -- The harvest-feed coupling reads the same setpoint (its r4 uses the
    ! terminal entry), so the sum check must hold whenever either consumer is active.
    IF ((ok_maturity_transfer .OR. ok_harvest_feed_coupling) .AND. &
         ABS(SUM(mat_target_frac) - un) > 1.E-6_r_std) THEN
       CALL ipslerr_p(3,'constantes','MAT_TARGET_FRAC must sum to one', &
            'the entries are area shares of the age-class group','')
    ENDIF

    !
    ! Guillaume M. -- AGE_CLASS_BOUNDS_PFT: see design/MODULE_DESIGN_AGE_CLASS_BOUNDS_PFT.md
    !
    !Config Key   = OK_AGE_CLASS_BOUND_PFT
    !Config Desc  = Derive age class bounds from each PFT largest_tree_dia
    !Config If    = OK_STOMATE
    !Config Def   = n
    !Config Help  = When y, age_class_bound(:,jv) = AGE_CLASS_FRAC(:)*largest_tree_dia(jv)
    !Config         for tree PFTs, so that every species reaches its last age class just
    !Config         before its clearcut diameter. n = single global vector for all PFTs
    !Config         (bit-identical to the historical behaviour).
    !Config         [CATEGORIE : ECHAFAUDAGE DE VALIDATION -- retirer quand le test §6.2
    !Config          du design AGE_CLASS_BOUNDS_PFT (les 4 classes atteignables) sera passe]
    !Config Units = [-]
    ok_age_class_bound_pft = .FALSE.
    CALL getin_p('OK_AGE_CLASS_BOUND_PFT', ok_age_class_bound_pft)

    ALLOCATE(age_class_frac(MAX(nagec-1,1)),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for age_class_frac. We stop.',nagec
       STOP 'constantes.f90'
    ENDIF
    IF (nagec == 4) THEN
       age_class_frac(1) = 0.25
       age_class_frac(2) = 0.50
       age_class_frac(3) = 0.75
    ELSE
       ! Guillaume M. -- Evenly spaced bounds for any other nagec.
       DO ic_loc = 1, MAX(nagec-1,1)
          age_class_frac(ic_loc) = REAL(ic_loc,r_std) / REAL(MAX(nagec,1),r_std)
       ENDDO
    ENDIF
    !Config Key   = AGE_CLASS_FRAC
    !Config Desc  = Age class bounds as a fraction of largest_tree_dia
    !Config If    = OK_AGE_CLASS_BOUND_PFT
    !Config Def   = 0.25, 0.50, 0.75
    !Config Help  = nagec-1 strictly increasing values in ]0,1[. /!\ LA DERNIERE VALEUR
    !Config         N'EST PLUS LIBRE depuis le 2026-07-28 : elle est DERIVEE de
    !Config         DIA_ROTATION_TOL (voir ci-dessous). La renseigner au namelist est sans
    !Config         effet sur le dernier element.
    !Config Units = [-]
    CALL getin_p('AGE_CLASS_FRAC', age_class_frac)

    ! Guillaume M. -- Cold-start area split between age classes. Allocated HERE, after NAGEC
    ! has been read: a vector sized before that read silently ends up with one element and
    ! only its first value is ever used. Default is an even split, so a cold start with no
    ! namelist entry builds a normal forest instead of putting everything in class 1.
    ALLOCATE(age_class_init_frac(MAX(nagec,1)),stat=ier)
    IF (ier /= 0) THEN
       WRITE(numout,*) 'Memory allocation error for age_class_init_frac. We stop.',nagec
       STOP 'constantes.f90'
    ENDIF
    age_class_init_frac(:) = un / REAL(MAX(nagec,1),r_std)
    !Config Key   = AGE_CLASS_INIT_FRAC
    !Config Desc  = Cold-start area fraction of each age class
    !Config If    = OK_STOMATE
    !Config Cat   = SCIENTIFIC OPTION
    !Config Def   = 0.25, 0.25, 0.25, 0.25 (an even split, i.e. 1/nagec)
    !Config Help  = nagec values, renormalised to 1. Sets how a cold start spreads each
    !Config         species area over its age classes. AGE_CLASS_INIT_FRAC = 1,0,0,0 puts
    !Config         everything in the youngest class, which is the historical behaviour.
    !Config Units = [-]
    CALL getin_p('AGE_CLASS_INIT_FRAC', age_class_init_frac)
    ! Guillaume M. -- A vector that does not sum to 1 is renormalised; an all-zero one is
    ! refused, because it would silently wipe the vegetation area at the cold start.
    acifsum = SUM(age_class_init_frac(1:MAX(nagec,1)))
    IF (acifsum <= zero .OR. ANY(age_class_init_frac(1:MAX(nagec,1)) < zero)) THEN
       WRITE(numout,*) 'AGE_CLASS_INIT_FRAC : ',age_class_init_frac(:)
       CALL ipslerr_p(3,'constantes','AGE_CLASS_INIT_FRAC must be >= 0 with a positive sum',&
            'a null or negative split would wipe the vegetation area at cold start','')
    ENDIF
    age_class_init_frac(:) = age_class_init_frac(:) / acifsum
    WRITE(numout,*) 'AGE_CLASS_INIT_FRAC (renormalise) : ',age_class_init_frac(:)

    ! Guillaume M. -- The lower bound of the LAST age class is DERIVED from DIA_ROTATION_TOL
    ! so that the class coincides exactly with the cuttable set. Kept independent, the two
    ! settings left a dead zone "mature but not cuttable" and needed a coherence guard.
    ! /!\ Narrowing the class shrinks A_mature, so f_cut = (A_forest/R)/A_mature rises; if
    ! it saturates at 1 the harvest reverts to all-or-nothing. Watched by validation test 5.5.
    IF (nagec > 1) THEN
       age_class_frac(nagec-1) = un - dia_rotation_tol
    ENDIF

    ! Guillaume M. -- Opening of age class 1 at establishment.
    ALLOCATE(prescribe_rdi_frac(MAX(nagec,1)),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for prescribe_rdi_frac. We stop.',nagec
       STOP 'constantes.f90'
    ENDIF
    ! Guillaume M. -- Negative sentinel = disabled: rdi_use keeps rdi_start(ivm), bit-neutral.
    prescribe_rdi_frac(:) = -un
    !Config Key   = PRESCRIBE_RDI_FRAC
    !Config Desc  = Establishment RDI as a fraction of rdi_max, per age class
    !Config If    = OK_STOMATE
    !Config Def   = -1.
    !Config Help  = nagec values in ]0,1]. Sets the stand density installed by prescribe
    !Config         to prescribe_rdi_frac(class)*rdi_max(pft,management) instead of
    !Config         rdi_start(pft). Typical: 0.20, 0.85, 0.95, 1.00 -- age class 1 far
    !Config         from the self-thinning ceiling (a genuine canopy gap), classes 2-4
    !Config         close to it. A fraction, not an absolute RDI, because rdi_max depends
    !Config         on the management type. Any value <= 0 disables the whole feature
    !Config         (bit-identical to the historical behaviour).
    !Config Units = [-]
    CALL getin_p('PRESCRIBE_RDI_FRAC', prescribe_rdi_frac)

    !Config Key   = OK_MIN_DENSITY_RESET
    !Config Desc  = Keep the minimum-density resets driven by dens_target
    !Config If    = OK_STOMATE
    !Config Def   = y
    !Config Help  = y = historical behaviour: a stand below dens_target is clearcut
    !Config         (sapiens_forestry CONDITION 1). n = that clearcut is disabled.
    !Config         Scope: CONDITION 1 only. The slow-death block of stomate_kill also
    !Config         keys on dens_target but has its own flag (OK_SLOW_DEATH_DENSITY),
    !Config         because the two are NOT of the same nature: CONDITION 1 resets the
    !Config         stand at once, the other is a progressive mortality sink.
    !Config         Measured on a 1-year panEU run: CONDITION 1 fires 407 times against
    !Config         21 for the rotation clearcut, 95 % of them in age class 1. This reset
    !Config         is also what forced the desync workaround in stomate_prescribe
    !Config         (a low establishment RDI falls below dens_target and is removed at
    !Config         once), so disabling it is a prerequisite for PRESCRIBE_RDI_FRAC to
    !Config         have any effect on age class 1.
    !Config         [CATEGORIE : ECHAFAUDAGE DE VALIDATION -- retirer avec
    !Config          OK_SLOW_DEATH_DENSITY, apres leur controle de signe pluriannuel]
    !Config Units = [-]
    ok_min_density_reset = .TRUE.
    CALL getin_p('OK_MIN_DENSITY_RESET', ok_min_density_reset)

    !Config Key   = OK_SLOW_DEATH_DENSITY
    !Config Desc  = Keep the slow death of stands below dens_target
    !Config If    = OK_STOMATE
    !Config Def   = n
    !Config Help  = y = historical behaviour: a stand below dens_target loses ndying stems
    !Config         per day (stomate_kill, circ_class_n -= ndying) and dies over
    !Config         ndying_year, which spreads the wood litter input over several years.
    !Config         n = that mortality sink is removed. MEASURED CONSEQUENCE: removing it
    !Config         raises the age-class-1 RDI from 0.325 to 0.400 over 4280 slots on a
    !Config         1-year panEU run, i.e. it makes age class 1 DENSER -- the opposite of
    !Config         the opening sought by PRESCRIBE_RDI_FRAC. Kept as a separate knob from
    !Config         OK_MIN_DENSITY_RESET precisely so the two effects stay attributable.
    !Config         [CATEGORIE : ECHAFAUDAGE DE VALIDATION -- garder SEPARE de
    !Config          OK_MIN_DENSITY_RESET : c'est cette separation qui a permis d'imputer les
    !Config          +0.074 de RDI classe 1 a ce bloc et non au RDI d'etablissement]
    !Config Units = [-]
    ok_slow_death_density = .FALSE.
    CALL getin_p('OK_SLOW_DEATH_DENSITY', ok_slow_death_density)

    !Config Key   = OK_AED_EDGE_RDI_WEIGHT
    !Config Desc  = Weight the AED age-class-1 edge source by the stand openness
    !Config If    = OK_EDGE_FROM_AGE_CLASS
    !Config Def   = n
    !Config Help  = When y, the regeneration term of the edge becomes
    !Config         alpha_reg * SUM(veget_max(class 1)) * MAX(0, 1-MIN(1, RDI)) instead of
    !Config         the raw area, so a regenerating stand stops producing edge as it
    !Config         stocks up. The recovery dynamics become emergent and subsume
    !Config         tau_rec / Trecruit. RDI is used rather than veget/veget_max because
    !Config         the latter does NOT respond to the clearcut (measured: median gap
    !Config         change +0.037, 66 % of 1109 cut events show no opening at all) and
    !Config         collapses every winter for deciduous PFTs. Requires
    !Config         OK_EDGE_FROM_AGE_CLASS: the weight applies to A_open, which only that
    !Config         path computes. Re-anchor AED_ALPHA_REG and MI_EDGE_SHAPE_COEF
    !Config         afterwards (median class-1 weight is about 0.675).
    !Config         [CATEGORIE : ECHAFAUDAGE DE VALIDATION -- retirer quand le controle de
    !Config          signe pluriannuel du poids (1-RDI) sera fait sur aedstandage]
    !Config Units = [-]
    ok_aed_edge_rdi_weight = .FALSE.
    CALL getin_p('OK_AED_EDGE_RDI_WEIGHT', ok_aed_edge_rdi_weight)

    IF (ok_aed_edge_rdi_weight .AND. .NOT. ok_edge_from_age_class) THEN
       CALL ipslerr_p(3,'constantes.f90',&
            'OK_AED_EDGE_RDI_WEIGHT requires OK_EDGE_FROM_AGE_CLASS',&
            'the weight applies to A_open, which only that path computes,',&
            'without it the flag would be silently inert')
    ENDIF

    IF (ANY(prescribe_rdi_frac(:) > zero)) THEN
       DO ic_loc = 1, nagec
          IF (prescribe_rdi_frac(ic_loc) <= zero .OR. prescribe_rdi_frac(ic_loc) > un) THEN
             WRITE(numout,*) 'PRESCRIBE_RDI_FRAC out of ]0,1] : ',prescribe_rdi_frac(:)
             CALL ipslerr_p(3,'constantes.f90',&
                  'PRESCRIBE_RDI_FRAC must lie in ]0,1] for every age class',&
                  'a fraction above 1 would exceed rdi_max and be thinned back at once',&
                  'set every class, or leave the whole array <= 0 to disable')
          ENDIF
       ENDDO
       IF (ok_min_density_reset) THEN
          CALL ipslerr_p(2,'constantes.f90',&
               'PRESCRIBE_RDI_FRAC with OK_MIN_DENSITY_RESET still on',&
               'an open age class 1 falls below dens_target and is clearcut at once,',&
               'so the prescribed opening will not survive its first year')
       ENDIF
    ENDIF

    !Config Key   = OK_CLEARCUT_LAST_CLASS
    !Config Desc  = Restrict management clearcut to the last age class
    !Config If    = OK_STOMATE
    !Config Def   = n
    !Config Help  = When y, the management clearcut triggers of sapiens_forestry_main
    !Config         (CONDITION 2/2bis diameter-age, CONDITION 3/3bis PAI<MAI) are refused
    !Config         outside the last age class of the group. Stand-failure resets
    !Config         (CONDITION 1), forced/spinup clearcuts, species change and all natural
    !Config         disturbances are NOT affected: they must be able to reset any class.
    !Config         [CATEGORIE : ECHAFAUDAGE DE VALIDATION -- retirer avec
    !Config          OK_AGE_CLASS_BOUND_PFT, meme test]
    !Config Units = [-]
    ok_clearcut_last_class = .FALSE.
    CALL getin_p('OK_CLEARCUT_LAST_CLASS', ok_clearcut_last_class)
    !
    !Config Key   = OK_RDI_BAND_CLASS1
    !Config Desc  = Age class 1 joins the thinned classes and gets the RDI band
    !Config If    = OK_CLEARCUT_LAST_CLASS
    !Config Def   = n
    !Config Help  = Guillaume M. -- OK_CLEARCUT_LAST_CLASS thins classes 2..nagec-1 only:
    !Config         class 1 is left RDI-FREE as a post-disturbance opening. That closes a
    !Config         trap on fertile-poor cells -- a dense class 1 can neither be thinned nor
    !Config         grow past age_class_bound(1), so it never reaches class 2 and the stand
    !Config         stalls (observed boreal, 40 yr at about 6000 stems/ha). When y, class 1
    !Config         is thinned like the intermediate classes and reads RDI_BAND_LOW/HIGH,
    !Config         which are per-slot hence per-age-class. The last class stays unthinned.
    !Config         /!\ Costs the edge signal its openness: class 1 is what AED_ALPHA_REG
    !Config         converts into edge, and thinning it lowers RDI, hence raises the
    !Config         (1-RDI) weight of OK_AED_EDGE_RDI_WEIGHT. Re-anchor before any EFDA
    !Config         comparison.
    !Config         [CATEGORIE : OPTION SCIENTIFIQUE -- choix de modelisation, permanent]
    !Config Units = [-]
    ok_rdi_band_class1 = .FALSE.
    CALL getin_p('OK_RDI_BAND_CLASS1', ok_rdi_band_class1)
    !
    !Config Key   = AC12_DIA_RATIO_DOWN
    !Config Desc  = Demotion threshold, as a fraction of the class-1 diameter bound
    !Config If    = OK_STOMATE
    !Config Def   = 1.0
    !Config Help  = Guillaume M. -- HYSTERESIS between promotion and demotion. Promotion
    !Config         fires on a HEIGHT ratio (AC12_H_RATIO_START = 0.85), which is 0.66 to
    !Config         0.72 of the bound in DIAMETER since pipe_tune3 spans 0.42-0.52, while
    !Config         demotion fired at max_dia < bound, a ratio of 1.0. The two windows
    !Config         overlapped, so a stand promoted below the bound was demoted the same
    !Config         year, whole slot -- measured 26 years of frozen area on a boreal pine.
    !Config         Set it BELOW the promotion ratio in diameter (so under about 0.66) to
    !Config         open a gap. 1.0 reproduces the previous behaviour exactly.
    !Config         /!\ Judged in DIAMETER, never in height: converting 0.85 by hand has
    !Config         already inverted one criterion in this repository.
    !Config         /!\ Never disable the demotion branch: it is the ONLY path by which
    !Config         storm and beetle send a wiped-out stand back to age class 1.
    !Config Units = [-]
    ac12_dia_ratio_down = un
    CALL getin_p('AC12_DIA_RATIO_DOWN', ac12_dia_ratio_down)
    !
    !Config Key   = ROTATION_GROWTH_WEIGHT
    !Config Desc  = Weight of site conditions against the published rotation itinerary
    !Config If    = OK_MANAGEMENT_INTENSITY
    !Config Def   = 0.
    !Config Help  = Guillaume M. -- The applied rotation becomes
    !Config            (1-w)*R_recommendation + w*R_growth
    !Config         w = 0 : the published itinerary alone, whatever the site can deliver.
    !Config                 This is the previous behaviour, hence bit-neutral.
    !Config         w = 1 : perfectly site-adapted management, the recommendation carries
    !Config                 no weight at all.
    !Config         Motive: target_rotation_age derives from rotation_ref, which is ONE
    !Config         value per PFT for the whole of Europe -- site fertility enters nowhere.
    !Config         Measured on a hemiboreal pixel: 49.6 yr prescribed against 97 yr for
    !Config         the model to reach the cut diameter and 84-87 yr observed, so the
    !Config         cadence claims stands that do not exist and empties the mature class.
    !Config         Convex, so the result never leaves the interval of the two terms.
    !Config         [CATEGORIE : OPTION SCIENTIFIQUE -- choix de modelisation, permanent]
    !Config Units = [-]
    rotation_growth_weight = zero
    CALL getin_p('ROTATION_GROWTH_WEIGHT', rotation_growth_weight)
    IF (rotation_growth_weight < zero .OR. rotation_growth_weight > un) THEN
       CALL ipslerr_p(3, 'constantes', 'ROTATION_GROWTH_WEIGHT must lie in [0,1]', &
            'it is the weight of a CONVEX combination', '')
    ENDIF
    !
    !Config Key   = DIA_GROWTH_TAU
    !Config Desc  = Memory of the diameter-increment integrator feeding R_growth
    !Config If    = ROTATION_GROWTH_WEIGHT > 0
    !Config Def   = 30.
    !Config Help  = Guillaume M. -- The increment is averaged over decades so the cadence
    !Config         does not chase interannual noise: a rotation is a slow quantity, and
    !Config         feeding it a noisy growth signal would make the harvest oscillate on
    !Config         top of the age structure it is meant to stabilise.
    !Config Units = [yr]
    dia_growth_tau = 30._r_std
    CALL getin_p('DIA_GROWTH_TAU', dia_growth_tau)
    dia_growth_tau = MAX(dia_growth_tau, un)
    !
    !Config Key   = DIA_GROWTH_AREA_TOL
    !Config Desc  = Relative area change above which a slot's increment is rejected
    !Config If    = ROTATION_GROWTH_WEIGHT > 0
    !Config Def   = 0.01
    !Config Help  = Guillaume M. -- A slot that exchanged area this year did not GROW, it
    !Config         was mixed: promotions bring in trees and merge_biomass_pfts dilutes the
    !Config         per-tree biomass, so max_dia barely moves. Measured, class 3 came out at
    !Config         0.097 cm/yr against 0.577 in class 1, giving a single 126-yr leg and a
    !Config         rotation of 165 yr for 85-87 observed. Only slots whose area held still
    !Config         carry a cohort, so only they feed the integrator. Raising this admits
    !Config         more samples at the cost of reintroducing the dilution bias.
    !Config Units = [-]
    dia_growth_area_tol = 0.01_r_std
    CALL getin_p('DIA_GROWTH_AREA_TOL', dia_growth_area_tol)
    dia_growth_area_tol = MAX(dia_growth_area_tol, zero)
    !
    !Config Key   = ROTATION_GROWTH_MIN
    !Config Desc  = Lower bound of the growth-derived rotation
    !Config If    = ROTATION_GROWTH_WEIGHT > 0
    !Config Def   = 10.
    !Config Help  = Guillaume M. -- Guards against an absurd rotation built on a barely
    !Config         populated slot. An increment near zero would otherwise send R_growth
    !Config         to infinity, and a spuriously large one would collapse it.
    !Config Units = [yr]
    rotation_growth_min = 10._r_std
    CALL getin_p('ROTATION_GROWTH_MIN', rotation_growth_min)
    !
    !Config Key   = ROTATION_GROWTH_MAX
    !Config Desc  = Upper bound of the growth-derived rotation
    !Config If    = ROTATION_GROWTH_WEIGHT > 0
    !Config Def   = 400.
    !Config Units = [yr]
    rotation_growth_max = 400._r_std
    CALL getin_p('ROTATION_GROWTH_MAX', rotation_growth_max)
    IF (ac12_dia_ratio_down <= zero .OR. ac12_dia_ratio_down > un) THEN
       CALL ipslerr_p(3, 'constantes', 'AC12_DIA_RATIO_DOWN must lie in ]0,1]', &
            'above 1 the demotion would fire on a stand that has reached its bound', '')
    ENDIF
    !
    ! Guillaume M. -- Without the last-class restriction the thinning interval is not
    ! defined, so the flag would be silently inert -- the trap this repository has paid for.
    IF (ok_rdi_band_class1 .AND. .NOT. ok_clearcut_last_class) THEN
       CALL ipslerr_p(3, 'constantes', 'OK_RDI_BAND_CLASS1=y requires OK_CLEARCUT_LAST_CLASS=y', &
            'it widens the thinned age-class interval, which only that flag defines', '')
    ENDIF
    !

    ! Guillaume M. -- Guards active only under the flag: with the flag off the fractions are
    ! unused and an existing configuration must not be invalidated.
    IF (ok_age_class_bound_pft) THEN
       IF (nagec < 2) THEN
          CALL ipslerr_p(3,'constantes.f90','OK_AGE_CLASS_BOUND_PFT needs NAGEC >= 2',&
               'with a single age class the bounds carry no function','')
       ENDIF
       DO ic_loc = 1, nagec-1
          IF (age_class_frac(ic_loc) <= zero .OR. age_class_frac(ic_loc) >= un) THEN
             WRITE(numout,*) 'AGE_CLASS_FRAC out of ]0,1[ : ',age_class_frac(:)
             CALL ipslerr_p(3,'constantes.f90','AGE_CLASS_FRAC must lie in ]0,1[','','')
          ENDIF
          IF (ic_loc > 1) THEN
             IF (age_class_frac(ic_loc) <= age_class_frac(ic_loc-1)) THEN
                WRITE(numout,*) 'AGE_CLASS_FRAC not strictly increasing : ',age_class_frac(:)
                CALL ipslerr_p(3,'constantes.f90','AGE_CLASS_FRAC must be strictly increasing','','')
             ENDIF
          ENDIF
       ENDDO
       ! Guillaume M. -- The clearcut trigger is largest_tree_dia*(1-dia_rotation_tol) and
       ! must fall INSIDE the last age class. Both sides are now equal by construction, so
       ! this guard no longer checks a tunable constraint but the derivation itself: it can
       ! only fail if someone removes the derivation and forgets this test.
       IF (ABS(age_class_frac(nagec-1) - (un - dia_rotation_tol)) > 1.e-12) THEN
          WRITE(numout,*) 'AGE_CLASS_FRAC(nagec-1), 1-DIA_ROTATION_TOL : ',&
               age_class_frac(nagec-1), un - dia_rotation_tol
          CALL ipslerr_p(3,'constantes.f90',&
               'AGE_CLASS_FRAC(nagec-1) doit EGALER 1-DIA_ROTATION_TOL',&
               'la derniere classe d age doit coincider avec l ensemble coupable',&
               'la derivation a-t-elle ete retiree ?')
       ENDIF
       IF (.NOT. ok_clearcut_last_class) THEN
          CALL ipslerr_p(1,'constantes.f90',&
               'OK_AGE_CLASS_BOUND_PFT without OK_CLEARCUT_LAST_CLASS',&
               'bounds are consistent but nothing forbids a clearcut in another class','')
       ENDIF
    ELSE
       IF (ok_clearcut_last_class) THEN
          CALL ipslerr_p(2,'constantes.f90',&
               'OK_CLEARCUT_LAST_CLASS without OK_AGE_CLASS_BOUND_PFT',&
               'with global bounds no PFT ever reaches the last age class,',&
               'so every management clearcut would be refused and the forest would freeze')
       ENDIF
    ENDIF

    ! Check for inconsistent configuration
    IF (spinup_analytic .AND. veget_update .GT.0 .AND. nagec .GT.1) THEN
       WRITE(numout,*) 'Use the analytical spinup, ',spinup_analytic
       WRITE(numout,*) 'Use land cover changes, ',veget_update
       WRITE(numout,*) 'Number of age classes, ',nagec
       CALL ipslerr_p(3,'The spinup cannot handle land cover changes',&
            'and age classes at the same time',&
            'If you really need a spinup with lcc you will have to',&
            'further develop the code in age_class_distr (sapiens_lcchange.f90)')
    END IF
    !
    !Config Key   = MIN_WATER_STRESS
    !Config Desc  = Minimal value for wstress_fac
    !Config If    = OK_STOMATE 
    !Config Def   = 0.1 
    !Config Help  = 
    !Config Units = [-]
    min_water_stress = 0.1
    CALL getin_p('MIN_WATER_STRESS',min_water_stress)
    !
    !Config Key   = NDIA_HARVEST
    !Config Desc  = Number of basal area classes in which the harvest is stored
    !Config If    = OK_STOMATE, OK_DIMENSIONAL_PRODUCT_USE 
    !Config Def   = 5 
    !Config Help  = Number of basal area classes in which the harvest
    !               is stored. This is useful when harvest is further
    !               used in a wood-use module
    !Config Units = [-]
    IF (ok_dimensional_product_use) THEN
       ndia_harvest = 7
    ELSE
       ndia_harvest = 1
    END IF
    CALL getin_p('NDIA_HARVEST',ndia_harvest)
    !
    !Config Key   = MAX_HARVEST_DIA
    !Config Desc  = The maximum diamter of tree which can be harvested
    !Config If    = OK_STOMATE 
    !Config Def   = 1.0
    !Config Help  = The maximum diamter of tree which can be 
    !               harvested.  Notice that we will create a class 
    !               that is one size larger than this to make sure
    !               we keep track of all the wood.
    !Config Units = [m]
    max_harvest_dia = 0.7
    CALL getin_p('MAX_HARVEST_DIA',max_harvest_dia) 

    !Config Key   = N_PAI
    !Config Desc  = Number of years used for the calculation of the periodic annual increment
    !Config If    = OK_STOMATE 
    !Config Def   = 5 
    !Config Help  = The number of years used for the cumulative
    !Config Help  = averages of the periodic annual increment. By
    !Config Help  = setting this value to 1000, this functionality is not used.
    !Config Help  = This way of management is extremly rational. Owners are
    !Config Help  = to follow yield tables than to monitor stand growth every 5 years.  
    !Config Units = [-]
    n_pai = 1000 
    CALL getin_p('N_PAI',n_pai)
    !
    !Config Key   = NTREES_PROFIT
    !Config Desc  = Number of trees below which the forest will be cut and replanted
    !Config If    = FOREST_MANAGEMENT  
    !Config Def   = 100 
    !Config Help  = The number of trees over which the average
    !               height is calculated to determine if the
    !               stand will be profitable to thin.
    !               Fewer trees than ntrees_profit is considered no
    !               longer profitable. Hence the stand will be cut 
    !               and replaced. It mimics disturbance in unmanaged forests
    !               and it prevents forest from becoming unrealistically old.
    !Config Units = [number of trees]
    ntrees_profit=1
    CALL getin_p('NTREES_PROFIT',ntrees_profit)
    !
    ! If we don't read the new species from a map, we will use veget_max
    ! instead. So we will replant with the current species. This will
    ! only happen if species_change_force equals -9999. If 
    ! species_change_force has a different value in the run.def that
    ! value will be used. Note that species_change_force is intended for
    ! testing and debugging.
    !Config Key   = SPECIES_CHANGE_FORCE
    !Config Desc  = New species after a final cut for testing and debugging only
    !Config If    = OK_STOMATE
    !Config Def   = -9999
    !Config Help  = If we don't read the new species from a map, we will use veget_max
    !               instead. So we will replant with the current species. This will
    !               only happen if species_change_force equals -9999. If 
    !               species_change_force has a different value in the run.def that
    !               value will be used. Note that species_change_force is intended for
    !               testing and debugging.
    !               This is the PFT number which is replanted after a
    !               clearcut, if such a thing is being done.
    !               To be used with lchange_species = .TRUE. and
    !               lread_species_change_map = .FALSE. The
    !               forced value is mainly useful for debugging
    !Config Units = [PFT number]
    species_change_force=-9999
    CALL getin_p('SPECIES_CHANGE_FORCE',species_change_force)

    !Config Key   = FM_CHANGE_FORCE
    !Config Desc  = New management after a final cut for testing and debugging only
    !Config If    = OK_STOMATE, LCHANGE_SPECIES
    !Config Def   = ifm_none
    !Config Help  = If we don't read the new speciesFM strategies from a map, we 
    !               will force it with this variable. Following a harvest all
    !               management will be set to unmanaged.
    !               This is the FM strategy which is used for the replant
    !               after a clearcut, if such a thing is being done.
    !               To be used with lchange_species = .TRUE. and
    !               lread_desired_fm_map = .FALSE. The forced value is
    !               mainly useful for debugging
    !Config Units = [1, 2, 3 or 4; unitless]
    fm_change_force = ifm_none
    CALL getin_p('FM_CHANGE_FORCE',fm_change_force)      

    ! Bark beetle attack module
    !Config Key   = nb_years_bgi
    !Config Desc  = numbers of years over which bark beetle generation index is calculated
    !Config If    = OK_PEST, OK_STOMATE
    !Config Def   = 3
    !Config Help  = numbers we have to average to calculate the bark beetle
    !generation index. use in beetle damage module
    !Config Units = [years]
    nb_years_bgi = 3
    CALL getin_p('NB_YEARS_BGI',nb_years_bgi)

    !-
    ! windthrow
    !-
    !Config Key   = DAILY_MAX_TUNE
    !Config Desc  = Non linear tuning factor for daily maximum wind speed used in windthrow module
    !Config If    = OK_WINDTHROW, stomate main program
    !Config Def   = 1.000
    !Config Help  =
    !Config Units = [-]
    daily_max_tune = 1.0
    CALL getin_p('DAILY_MAX_TUNE',daily_max_tune)

    !Config Key   = ELEVATE_WIND
    !Config Desc  = Non linear tuning factor for daily maximum wind speed used
    !in windthrow module
    !Config If    = OK_WINDTHROW, stomate main program
    !Config Def   = 10.0
    !Config Help  =
    !Config Units = [-]
    elevate_wind = 10.0
    CALL getin_p('ELEVATE_WIND',elevate_wind)

    !Config Key   = H_THRESHOLD_WIND
    !Config Desc  = Ratio of height threshold to pipe_tune2 for truncating wind damage
    !Config If    = OK_WINDTHROW, stomate main program
    !Config Def   = 0.3
    !Config Help  =
    !Config Units = [-]
    h_threshold_wind = 0.3
    CALL getin_p('H_THRESHOLD_WIND',h_threshold_wind)

    !Config Key   = D_THRESHOLD_WIND
    !Config Desc  = Fraction of the largest tree diameter used as reference threshold
    !Config If    = OK_WINDTHROW, stomate main program
    !Config Def   = 0.7
    !Config Help  =
    !Config Units = [-]
    d_threshold_wind = 0.7
    CALL getin_p('D_THRESHOLD_WIND',d_threshold_wind)

    !Config Key   = CLEAR_CUT_MAX
    !Config Desc  = The maximum contiguous area clearfelled used to calculate edge length
    !Config If    = OK_WINDTHROW, OK_STOMATE
    !Config Def   = 20000.
    !Config Help  =
    !Config Units = m2
    clear_cut_max = 20000.0
    CALL getin_p('CLEAR_CUT_MAX',clear_cut_max)
    
    !Config Key   = LEGACY_YEARS_WIND
    !Config Desc  = The years used to calculate the maximum/mean/sum related to wind speed
    !Config If    = OK_WINDTHROW, OK_STOMATE
    !Config Def   = 5
    !Config Help  =
    !Config Units = years
    legacy_years_wind = 5
    CALL getin_p('LEGACY_YEARS_WIND',legacy_years_wind)

    !Config Key   = LEGACY_YEARS
    !Config Desc  = The years used to calculate the total damage from disturbances
    !Config If    = OK_WINDTHROW, OK_STOMATE
    !Config Def   = 3
    !Config Help  =
    !Config Units = years
    legacy_years = 3
    CALL getin_p('LEGACY_YEARS',legacy_years)

    !Config Key   = LEGACY_YEARS_WOOD
    !Config Desc  = The years used to calculate the wood leftover legacy (bark beetles)
    !Config If    = OK_PEST, OK_STOMATE
    !Config Def   = 
    !Config Help  =
    !Config Units = years
    legacy_years_wood = 1 
    CALL getin_p('LEGACY_YEARS_WOOD',legacy_years_wood)

    !Config Key   = BEETLE_LEGACY
    !Config Desc  = The years used to calculate wood leftover
    !Config If    = OK_PEST, OK_STOMATE
    !Config Def   = 6
    !Config Help  =
    !Config Units = years
    beetle_legacy = 6 
    CALL getin_p('BEETLE_LEGACY',beetle_legacy)
    
    !Config Key   = WIND_DAYS
    !Config Desc  = the number of days to save sum of wind speed overtime  
    !Config If    = OK_WINDTHROW, stomate main program
    !Config Def   = 3
    !Config Help  =
    !Config Units = day
    wind_days = 3
    CALL getin_p('WIND_DAYS',wind_days)

    !Config Key   = GAP_THRESHOLD
    !Config Desc  = The threshold for the basal area loss to be considered to create the gap with edges
    !Config If    = OK_WINDTHROW, stomate main program
    !Config Def   = 0.3
    !Config Help  =
    !Config Units = -
    gap_threshold = 0.3
    CALL getin_p('GAP_THRESHOLD',gap_threshold)

    !Config Key   = NB_DAYS_STORM
    !Config Desc  = the number of days at which the max wind speed is less than wind_speed_storm_thr
    !Config If    = OK_WINDTHROW, stomate main program
    !Config Def   = 5
    !Config Help  =
    !Config Units = days
    nb_days_storm = 5
    CALL getin_p('NB_DAYS_STORM',nb_days_storm)

    !! Scale dependent parameters to detect wind storms
    !Config Key   = WIND_SUM_THRESHOLD
    !Config Desc  = Threshold for the sum of wind_days of wind speed for all
    ! time stpe to calculate wind damage 
    !Config If    = OK_WINDTHROW, stomate main program
    !Config Def   = 200
    !Config Help  =
    !Config Units = m/s
    
    !Config Key   = WIND_RATIO_THRESHOLD
    !Config Desc  = Threshold for the ratio of wind to longterm wind for each ime stpe to sum ratios
    !Config If    = OK_WINDTHROW, stomate main program
    !Config Def   = 3.0
    !Config Help  =
    !Config Units = -

    !Config Key   = WIND_SPEED_STORM_THR
    !Config Desc  = the wind speed threshold above which is_storm flag is set to TRUE 
    !Config If    = OK_WINDTHROW, stomate main program
    !Config Def   = 20.000
    !Config Help  =
    !Config Units = meter per second     
    SELECT CASE (forcing_resolution)
    CASE('2degrees_6hourly')
       ! Parameters for CRU_JRA (Jina Jeong)
       wind_ratio_threshold = 3
       wind_sum_threshold = 200
       wind_speed_storm_thr = 20.
    CASE('1degree_daily')
       ! Parameters for ISIMIP (Pengyi Zhang)
       wind_ratio_threshold = 1.8
       wind_sum_threshold = 190
       wind_speed_storm_thr = 7.8
    CASE('0.5degree_6hourly')
       ! Parameters for CRU_JRA (Jina Jeong)
       wind_ratio_threshold = 3
       wind_sum_threshold = 200
       wind_speed_storm_thr = 20.
    CASE('8km_hourly','site')
       ! Parameters for SAFRAN and FLUXNET (Guillaume Marie & Tom Bastiaan)
       wind_ratio_threshold = 3
       wind_sum_threshold = 200
       wind_speed_storm_thr = 16
       !+++DELETE+++
       ! Delete once the parameters have been set
       WRITE(numout,*) 'Parameters based on a multiplication of mas woind speed and duration'
       CALL ipslerr_p(3,'Remember to specify the resolution of the forcing', &
            'Cannot select the parameters to detect wind storms','','')
       !++++++++++++
    CASE('LMDzOR_zoomed')
       ! Parameters for zoom over Europe
       ! grid: 142*144 - 15 minute time steps (Marine Lanet & Pengyi Zhang)
       wind_ratio_threshold = -9999
       wind_sum_threshold = -9999
       wind_speed_storm_thr = -9999
       !+++DELETE+++
       ! Delete once the parameters have been set
       CALL ipslerr_p(3,'Remember to specify the resolution of the forcing', &
            'Cannot select the parameters to detect wind storms','','')
       !++++++++++++
    CASE('LMDzOR')
       ! Parameters for global LMDzOR 142*144 - 15 minute time steps
       wind_ratio_threshold = -9999
       wind_sum_threshold = -9999
       wind_speed_storm_thr = -9999
       !+++DELETE+++
       ! Delete once the parameters have been set
       CALL ipslerr_p(3,'Remember to specify the resolution of the forcing', &
            'Cannot select the parameters to detect wind storms','','')
       !++++++++++++
     CASE DEFAULT
       CALL ipslerr_p(3,'Remember to specify the resolution of the forcing', &
            'Cannot select the parameters to detect wind storms','','')
    ENDSELECT
    CALL getin_p('WIND_RATIO_THRESHOLD',wind_ratio_threshold)
    ! JJ 2026: Gompertz weight params for the storm wind-ratio damage weighting (replaces
    ! OK_WIND_ACCLIM). Defaults a=1/b=100/c=0 -> weight ~= 1 (near no-op); calibrate b,c.
    !Config Key   = WIND_RATIO_A
    !Config Desc  = Asymptote of the Gompertz weight applied to storm damage
    !Config If    = OK_STOMATE
    !Config Cat   = SCIENTIFIC OPTION
    !Config Def   = 1.0
    !Config Help  = The three WIND_RATIO_* parameters shape a Gompertz weight
    !Config         w = a * exp(-exp(-b * (ratio - c))) applied to the storm damage
    !Config         (stomate_windthrow, weight_ratio). Defaults a=1, b=100, c=0 give a
    !Config         steep step at ratio ~= 0, hence w ~= 1 for any positive ratio (near
    !Config         no-op). The weight can only REDUCE damage (0..a); it cannot create
    !Config         damage where the modified wind stays below the CWS.
    !Config Units = [-]
    wind_ratio_a = 1.0
    !Config Key   = WIND_RATIO_B
    !Config Desc  = Displacement of the Gompertz weight applied to storm damage
    !Config If    = OK_STOMATE
    !Config Cat   = SCIENTIFIC OPTION
    !Config Def   = 100.0
    !Config Help  = Larger b shifts the weight towards low values at small wind ratios.
    !Config         See WIND_RATIO_A for the full expression.
    !Config Units = [-]
    wind_ratio_b = 100.0
    !Config Key   = WIND_RATIO_C
    !Config Desc  = Growth rate of the Gompertz weight applied to storm damage
    !Config If    = OK_STOMATE
    !Config Cat   = SCIENTIFIC OPTION
    !Config Def   = 0.0
    !Config Help  = c = 0 makes the weight independent of the wind ratio, which is the
    !Config         neutral default. See WIND_RATIO_A for the full expression.
    !Config Units = [-]
    wind_ratio_c = 0.0
    CALL getin_p('WIND_RATIO_A',wind_ratio_a)
    CALL getin_p('WIND_RATIO_B',wind_ratio_b)
    CALL getin_p('WIND_RATIO_C',wind_ratio_c)
    CALL getin_p('WIND_SUM_THRESHOLD',wind_sum_threshold)
    CALL getin_p('WIND_SPEED_STORM_THR',wind_speed_storm_thr)

    !Config Key   = REPLACE_THR_WIND
    !Config Desc  = the threshold for the ratio of biomss revmoved from storm event to replace stand 
    !               to the youngest age class
    !Config If    = OK_WINDTHROW, stomate main program
    !Config Def   = 0.25
    !Config Help  = Indiviudal mortality vs stand replacing disturbance
    !Config Units = unitless
    replace_threshold_wind = 0.25
    CALL getin_p('REPLACE_THR_WIND',replace_threshold_wind)

    !Config Key   = REPLACE_THR_BEETLE
    !Config Desc  = the threshold for the ratio of biomss revmoved from bark
    !               beetle outbreak to replace stand to the youngest age class
    !Config If    = OK_PEST, stomate main program
    !Config Def   = 0.05
    !Config Help  = Indiviudal mortality vs stand replacing disturbance
    !Config Units = unitless
    replace_threshold_beetle = 0.05
    CALL getin_p('REPLACE_THR_BEETLE',replace_threshold_beetle)
   
    !Config Key   = REPLACE_THR_FIRE
    !Config Desc  = the threshold for the ratio of biomss revmoved from fire
    !               to replace stand to the youngest age class
    !Config If    = OK_SPITFIRE, stomate main program
    !Config Def   = 0.1
    !Config Help  = Indiviudal mortality vs stand replacing disturbance
    !Config Units = unitless
    replace_threshold_fire = 0.1
    CALL getin_p('REPLACE_THR_FIRE',replace_threshold_fire)    
     
    !Config Key   = FORCED_CLEAR_CUT
    !Config Desc  = Use to force a clear cut at a specific year during a simulation.
    !Config If    = OK_STOMATE
    !Config Def   = .FALSE.
    !Config Help  = 
    !Config Units = year
    forced_clear_cut= .FALSE.
    CALL getin_p('FORCED_CLEAR_CUT',forced_clear_cut)

    !Config Key   = NOUTDIACLASS
    !Config Desc  = Number of diameter classes in the output files
    !Config If    = OK_STOMATE
    !Config Def   = 16
    !Config Help  = ncirc is the number of circumference classes
    !               ORCHIDEE. For tree demography it could be more 
    !               convenient to get the results in diameter
    !               classes with fixed boundaries.
    !Config Units = [-]
    noutdiaclass = 16 
    CALL getin_p('NOUTDIACLASS',noutdiaclass)

    ALLOCATE(out_dia_class(noutdiaclass+1),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for out_dia_class. We stop. We need nagec words',noutdiaclass
       STOP 'constantes.f90'
    ENDIF

    !Config Key   = OUT_DIA_CLASS
    !Config Desc  = Boundaries of diameter classes in the output files
    !Config If    = OK_STOMATE
    !Config Def   = zero, .01, .05, .1, .15, .2, .3, .4, .5, .6, .7, .8, .9, 1., 1.5, 2., 9.9
    !Config Help  = Boundaries of the diameter classes used in an
    !               additional output file. Given the way this part is
    !               coded, this list should start with zero and end with
    !               a very high and thus unlikely diameter. The length of
    !               the list is thus noutdiaclass+1.
    !Config Units = [m]  
    IF (noutdiaclass.EQ.16) THEN
       out_dia_class =   (/  0.0, 0.01, 0.05, 0.1, 0.15, 0.2, 0.3, &      
                 0.4,  0.5,  0.6, 0.7,  0.8, 0.9, 1.0, &  
                 1.5,  2.0,  9.9 /)   
    ELSEIF(noutdiaclass.GT.16) THEN
       CALL ipslerr_p(3,'not enough output variables have been defined', &
            'only 16 diameter classes have been defined for xios', &
            'have a look at src_xml and define more classes in the', &
            'field and file definition files')
    ELSE
       out_dia_class(:) = -9999
    END IF
    CALL getin_p('OUT_DIA_CLASS',out_dia_class)
   
    IF (out_dia_class(1).EQ.-9999) THEN
       WRITE(numout,*) 'The code does not contain default values for out_dia_class'
       WRITE(numout,*) 'for this number of output diameter classes (noutdiaclass) ', noutdiaclass
       CALL ipslerr_p (3,'constantes.f90', 'Define bounds of output diameter classes',&
            'add default values in constantes.f90','or add values in the orchidee_pft.def')
    END IF
    DO iout = 2,noutdiaclass+1
       WRITE(striout,'(I2.2)') iout-1
       WRITE(numout,*) 'OUT_CLASS '//TRIM(striout), &
            ' in XIOS contains trees with diameters' 
       WRITE(numout,*) 'between ',out_dia_class(iout-1),' and ', &
            out_dia_class(iout)
    END DO

    !Config Key   = USE_HEIGHT_DOM
    !Config Desc  = Use the dominant vegetation height instead of the average height when calculating roughness length
    !Config If    = OK_STOMATE
    !Config Def   = .FALSE.
    !Config Help  = This is a somewhat odd switch, since it impacts sechiba
    !               but can only be used when stomate is activated.
    !               A logical flag determining if we use  the dominant height of
    !               the vegetation (.TRUE.) or the quadratic mean height (.FALSE.)
    !Config Units = [-]
    use_height_dom= .FALSE.
    CALL getin_p('USE_HEIGHT_DOM',use_height_dom)

    !Config Key   = ERR_ACT
    !Config Desc  = Action following an error
    !Config If    = OK_STOMATE
    !Config Def   = 1
    !Config Help  = The code distinguishes between two options to check for mass balance
    !               problems:
    !               ERR_ACT = 1 is recommended when running global
    !               long-term simulations. Under this option, mass balance closure is
    !               checked for all biogeochemical processes but only at the highest level
    !               thus stomate.f90 and stomate_lpj.f90. Although the mass balance checks
    !               are not very expensive in terms of computer time, skipping the numerous
    !               lower level checks is expected to save some time. Under this option the
    !               mass balance error is only written to the history file. No information
    !               is provided in which subroutine the problem occurred.
    !               ERR_ACT = 2
    !               is recommended when developing and testing the model. Now the mass
    !               balance is explicitly checked in stomate.f90, stomate_lpj.f90 and all
    !               its subroutines. Under this option the mass balance error is written to
    !               the history file and if the mass balance is not closed, the warning
    !               message will indicate in which subroutine the problem likely
    !               originated. 
    !               ERR_ACT = 3 is recommended when having a problem with
    !               mass balance closure. The mass balance is explicitly checked in
    !               stomate.f90, stomate_lpj.f90 and all its subroutines. If a mass balance
    !               occurs, the model is stopped. 
    !               ERR_ACT = 4 is for mass balance enthousiasts.
    !               If the model crashed on a mass balance issue in stomate_growth_fun_all
    !               this setting will activate intermediate checks. Intermediate checks may 
    !               help to narrow down in which part of the code the problem occurs. It also
    !               activates nbp consistency checks for those interested in better defining
    !               nep, nbp and tracing all the C and N in stomate_lpj.f90
    !Config Units = [1: write to history file, 2: warn and write to history file, and 3&4: stop the model]
    err_act = 1
    CALL getin_p('ERR_ACT',err_act)

    SELECT CASE (err_act)
    CASE(1)
        IF (printlev>=2) WRITE(numout,*) 'ERR_ACT=1: Many ipserr messages related to conservation problems will be skipped'
        plev = 0
    CASE(2)
        IF (printlev>=2) WRITE(numout,*) 'ERR_ACT=2: In case of a conservation problem the user will be warned'
        plev = 2
    CASE(3, 4)
        IF (printlev>=2) WRITE(numout,*) 'ERR_ACT=3 or 4: In case of a conservation problem the model will be stopped'
        plev = 3
    CASE DEFAULT
        CALL ipslerr_p(3,'Error in constantes.f90',&
            'Unknown value for ERR_ACT',&
            'Not clear how to treat errors','')   
    ENDSELECT

    !Config Key   = OK_FORCE_PHENO
    !Config Desc  = Use to force phenology when the conditions are not suitable
    !Config If    = OK_STOMATE
    !Config Def   = .TRUE.
    !Config Help  = Temperature phenology is very predictable and happens every year but
    !               moisture-driven phenology is less stable because the conditions
    !               may not be satisfied during 12 months resulting in a year without
    !               a canopy. If the model has passed the average budbreak day by a
    !               prescribed PFT-specific offset (::force_pheno), budbreak will be forced
    !               (given the condition that buds are avaialble)
    !Config Units = [-]
    ok_force_pheno = .TRUE.
    CALL getin_p('OK_FORCE_PHENO',ok_force_pheno)

    !Config Key   = OK_CROPS_IPRESENESCENCE
    !Config Desc  = Flag to allow crops to enter ipresenesence. Allows
    !allocation to 'fruits'.
    !Config If    = OK_STOMATE
    !Config Def   = .FALSE.
    !Config Help  =  This forces crops to spend a period of time, before 
    !                senescence, allocating carbon to 'fruits' instead of
    !                leaves and therefore decrease foliage and LAI.
    !Config Units = [-]
    ok_crops_ipresenescence = .FALSE.
    CALL getin_p('OK_CROPS_IPRESENESCENCE',ok_crops_ipresenescence)
    !
    !Config Key   = NPP_LONGTERM_INIT
    !Config Desc  = 
    !Config If    = OK_STOMATE
    !Config Def   = 10.
    !Config Help  = 
    !Config Units = [gC/m^2/year]
    npp_longterm_init = 10.
    CALL getin_p('NPP_LONGTERM_INIT',npp_longterm_init)
    !
    !Config Key   = SPAT_MOD_SELF_THIN
    !Config Desc  = Spatial modifier of alpha_self_thinning
    !Config If    = OK_STOMATE
    !Config Def   = 0.15
    !Config Help  = Spatial modifier of alpha_self_thinning during the spinup phase. This parameter truncates
    !               the min/max change of alpha_self_thinning 
    !Config Units = -
    spat_mod_self_thin = 0.15
    CALL getin_p('SPAT_MOD_SELF_THIN',spat_mod_self_thin)
    !
    !Config Key   = SPAT_EXP_SELF_THIN
    !Config Desc  = Spatial exponent of alpha_self_thinning
    !Config If    = OK_STOMATE
    !Config Def   = 0.33
    !Config Help  = Spatial exponent of alpha_self_thinning during the spinup phase. This parameter account
    !               for geomteric scaling between the diameter (in the self thinning relationship) and
    !               biomass (as a proxy for carrying capacity)
    !Config Units = -
    spat_exp_self_thin = 0.33
    CALL getin_p('SPAT_EXP_SELF_THIN',spat_exp_self_thin)


    !Config Key   = FLUX_TOT_COEFF
    !Config Desc  = coeff related to decomposition on peatland, esp agri_peat
    !Config if    = OK_STOMATE and agri_peat
    !Config Def   = 1.2, 1.4,.75
    !Config Help  =
    !Config Units = [days] 
    flux_tot_coeff = (/ 1.2, 1.4, .75/)  
    CALL getin_p('FLUX_TOT_COEFF',flux_tot_coeff)

    !Config Key   = P_A
    !Config Desc  = acrotelm density (g/m3)
    !Config if    = OK_STOMATE and ok_peat_NoDiscretisation
    !Config Def   = 3.5E+4
    !Config Help  =
    !Config Units = [(g/m3)] 
    p_A=3.5E+4
    CALL getin_p('P_A',p_A)

    !Config Key   = P_C
    !Config Desc  = catotelm density (g/m3)
    !Config if    = OK_STOMATE and ok_peat_NoDiscretisation
    !Config Def   = 9.1E+4 
    !Config Help  =
    !Config Units = [(g/m3)] 
    p_C=9.1E+4 
    CALL getin_p('P_C',p_C)

    !Config Key   = CF_A
    !Config Desc  = carbon fraction in catotelm peat
    !Config if    = OK_STOMATE and ok_peat_NoDiscretisation
    !Config Def   = 0.50 
    !Config Help  =
    !Config Units = [] 
    cf_A=0.50 
    CALL getin_p('CF_A',cf_A)

    !Config Key   = CF_C
    !Config Desc  = carbon fraction in acrotelm peat
    !Config if    = OK_STOMATE and ok_peat_NoDiscretisation
    !Config Def   = 0.52 
    !Config Help  =
    !Config Units = [] 
    cf_C=0.52 
    CALL getin_p('CF_C',cf_C)

    !Config Key   = V_RATIO
    !Config Desc  = ratio anaerobic-aerobic CO2
    !Config if    = OK_STOMATE and ok_peat_NoDiscretisation
    !Config Def   = 0.35 
    !Config Help  =
    !Config Units = [] 
    v_ratio=0.35 
    CALL getin_p('V_RATIO',v_ratio)

    !Config Key   = KA_INI
    !Config Desc  = acrotelm decomposition rate (an-1)
    !Config if    = OK_STOMATE and ok_peat_NoDiscretisation
    !Config Def   = 0.067
    !Config Help  =
    !Config Units = [year-1] 
    KA_ini=0.067
    CALL getin_p('KA_INI',KA_ini)

    !Config Key   = KP_INI
    !Config Desc  = carbon fraction in acrotelm peat
    !Config if    = OK_STOMATE and ok_peat_NoDiscretisation
    !Config Def   = 1.91E-2
    !Config Help  =
    !Config Units = [] 
    KP_ini=1.91E-2
    CALL getin_p('KP_INI',KP_ini)

    !Config Key   = KC_INI
    !Config Desc  = catotelm decomposition rate (an-1)
    !Config if    = OK_STOMATE and ok_peat_NoDiscretisation
    !Config Def   = 3.35E-5
    !Config Help  =
    !Config Units = [year-1] 
    KC_ini=3.35E-5
    CALL getin_p('KC_INI',KC_ini)
    

  END SUBROUTINE config_stomate_parameters

!! ================================================================================================================================
!! SUBROUTINE   : config_dgvm_parameters 
!!
!>\BRIEF        This subroutine reads in the configuration file all the parameters 
!! needed when the DGVM model is activated (ie : when ok_dgvm is set to true).
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCE(S) :
!!
!! FLOWCHART    :
!! \n
!_ ================================================================================================================================

  SUBROUTINE config_dgvm_parameters   

    IMPLICIT NONE

    !! 0. Variables and parameters declaration

    !! 0.4 Local variables

    !_ ================================================================================================================================    

    !-
    ! establish parameters
    !-
    !
    !Config Key   = ESTAB_MAX_TREE
    !Config Desc  = Maximum tree establishment rate 
    !Config If    = OK_DGVM
    !Config Def   = 0.12 
    !Config Help  = 
    !Config Units = [-]
    estab_max_tree = 0.12
    CALL getin_p('ESTAB_MAX_TREE',estab_max_tree)
    !
    !Config Key   = ESTAB_MAX_GRASS
    !Config Desc  = Maximum grass establishment rate
    !Config If    = OK_DGVM
    !Config Def   = 0.12 
    !Config Help  = 
    !Config Units = [-]
    estab_max_grass = 0.12
    CALL getin_p('ESTAB_MAX_GRASS',estab_max_grass)
    !
    !Config Key   = ESTABLISH_SCAL_FACT
    !Config Desc  = 
    !Config If    = OK_DGVM 
    !Config Def   = 5.
    !Config Help  = 
    !Config Units = [-]
    establish_scal_fact = 5.
    CALL getin_p('ESTABLISH_SCAL_FACT',establish_scal_fact)
    !
    !Config Key   = MAX_TREE_COVERAGE 
    !Config Desc  = 
    !Config If    = OK_DGVM 
    !Config Def   = 0.98
    !Config Help  = 
    !Config Units = [-]
    max_tree_coverage = 0.98
    CALL getin_p('MAX_TREE_COVERAGE',max_tree_coverage)
    !
    !Config Key   = IND_0_ESTAB
    !Config Desc  = 
    !Config If    = OK_DGVM 
    !Config Def   = 0.2
    !Config Help  = 
    !Config Units = [-]
    ind_0_estab = 0.2
    CALL getin_p('IND_0_ESTAB',ind_0_estab)

    !-
    ! light parameters
    !-
    !
    !Config Key   = ANNUAL_INCREASE
    !Config Desc  = for diagnosis of fpc increase, compare today's fpc to last year's maximum (T) or to fpc of last time step (F)?
    !Config If    = OK_DGVM
    !Config Def   = y
    !Config Help  = 
    !Config Units = [FLAG]
    annual_increase = .TRUE.
    CALL getin_p('ANNUAL_INCREASE',annual_increase)
    !
    !Config Key   = MIN_COVER 
    !Config Desc  = For trees, minimum fraction of crown area occupied 
    !Config If    = OK_DGVM
    !Config Def   = 0.05 
    !Config Help  = 
    !Config Units = [-]
    min_cover = 0.05
    CALL getin_p('MIN_COVER',min_cover)

    !-
    ! pftinout parameters
    !
    !Config Key   = IND_0 
    !Config Desc  = initial density of individuals
    !Config If    = OK_DGVM
    !Config Def   = 0.02 
    !Config Help  = 
    !Config Units = [-]
    ind_0 = 0.02
    CALL getin_p('IND_0',ind_0)
    !
    !Config Key   = MIN_AVAIL
    !Config Desc  = minimum availability
    !Config If    = OK_DGVM
    !Config Def   = 0.01
    !Config Help  = 
    !Config Units = [-]
    min_avail = 0.01   
    CALL getin_p('MIN_AVAIL',min_avail)
    !
    !Config Key   = RIP_TIME_MIN
    !Config Desc  = 
    !Config If    = OK_DGVM
    !Config Def   = 1.25 
    !Config Help  = 
    !Config Units = [year]
    RIP_time_min = 1.25
    CALL getin_p('RIP_TIME_MIN',RIP_time_min) 
    !
    !Config Key   = EVERYWHERE_INIT
    !Config Desc  = 
    !Config If    = OK_DGVM
    !Config Def   = 0.05 
    !Config Help  = 
    !Config Units = [-]
    everywhere_init = 0.05
    CALL getin_p('EVERYWHERE_INIT',everywhere_init)

  END SUBROUTINE config_dgvm_parameters

END MODULE constantes
