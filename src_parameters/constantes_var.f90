!!=================================================================================================================================
! MODULE       : constantes_var
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        constantes_var module contains most constantes like pi, Earth radius, etc...
!!              and all externalized parameters except pft-dependent constants.
!!
!!\n DESCRIPTION: This module contains most constantes and the externalized parameters of ORCHIDEE which 
!!                are not pft-dependent.\n
!!                In this module, you can set the flag diag_qsat in order to detect the pixel where the
!!                temperature is out of range (look qsatcalc and dev_qsatcalc in qsat_moisture.f90).\n
!!                The Earth radius is approximated by the Equatorial radius.The Earth's equatorial radius a,
!!                or semi-major axis, is the distance from its center to the equator and equals 6,378.1370 km.
!!                The equatorial radius is often used to compare Earth with other planets.\n
!!                The meridional mean is well approximated by the semicubic mean of the two axe yielding 
!!                6367.4491 km or less accurately by the quadratic mean of the two axes about 6,367.454 km
!!                or even just the mean of the two axes about 6,367.445 km.\n
!!                This module is already USE in module constantes. Therefor no need to USE it seperatly except
!!                if the subroutines in module constantes are not needed.\n
!!                
!! RECENT CHANGE(S):
!!
!! REFERENCE(S)	: 
!! - Louis, Jean-Francois (1979), A parametric model of vertical eddy fluxes in the atmosphere. 
!! Boundary Layer Meteorology, 187-202.\n
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_parameters/constantes_var.f90 $
!! $Date: 2026-05-13 15:11:27 +0200 (mer. 13 mai 2026) $
!! $Revision: 9540 $
!! \n
!_ ================================================================================================================================

MODULE constantes_var

  USE defprec

  IMPLICIT NONE
!-

                         !-----------------------!
                         !  ORCHIDEE CONSTANTS   !
                         !-----------------------!

  !
  ! FLAGS 
  !
  ! Flags represent different option in the code. Both options of a flag
  ! (true and false) should be well developped. A flag is an option we
  ! want to keep for a long time. Options that are no longer used may
  ! be removed and with them their flag will be removed.
  !
 
  INTEGER, SAVE       :: energy_control            !! Flag that automatically controls several other flags related to 
                                                   !! multi-layering (1/2/3).
!$OMP THREADPRIVATE(energy_control)
  
  LOGICAL, SAVE       :: ok_hydrol_arch            !! Flag that activates the hydraulic architecture routine (true/false)
!$OMP THREADPRIVATE(ok_hydrol_arch)

  LOGICAL       :: ok_hydrol_arch_muff       !! Flag that activates the radial resolution of Richard's equation around the roots. 
!$OMP THREADPRIVATE(ok_hydrol_arch_muff)

  LOGICAL       :: ok_hydrol_arch_dyn_res    !! Flag that activates the dynamic resistance equations used by Yao et al. (2021) (Under review at the time this is developed)
!$OMP THREADPRIVATE(ok_hydrol_arch_dyn_res)

  LOGICAL       :: split_soil_properties           !! Flag that permits to distinguish properties between the superficial and inferios soil horizon
!$OMP THREADPRIVATE(split_soil_properties)

  LOGICAL       :: ok_hydrol_arch_storage    !! Flag that activates the resolution of water storage inside the vegetation as in Tuzet et al. (2017).
!$OMP THREADPRIVATE(ok_hydrol_arch_storage)


  LOGICAL, SAVE       :: ok_gs_feedback            !! Flag that deactivates the stress applied to transpiration                                                              
!$OMP THREADPRIVATE(ok_gs_feedback)

  LOGICAL, SAVE       :: ok_mleb                   !! Flag that activates the multilayer energy budget (true/false)                                                              
!$OMP THREADPRIVATE(ok_mleb)
  
  LOGICAL, SAVE       :: ok_mleb_history_file      !! Flag that controls the writing of an aditional output file (true/false).
!$OMP THREADPRIVATE(ok_mleb_history_file)

  LOGICAL, SAVE :: ok_read_fm_map                  !! Thal to read in the forest management strategy from a map.                                                             
!$OMP THREADPRIVATE(ok_read_fm_map)

  CHARACTER(LEN=30), SAVE :: concept_scale         !! Set the scale for scale dependent parameters and processes.
!$OMP THREADPRIVATE(concept_scale)

  CHARACTER(LEN=30), SAVE :: forcing_resolution    !! Spatial and temporal resolution of the climate forcing
!$OMP THREADPRIVATE(forcing_resolution)

  LOGICAL, SAVE :: ok_read_sp_clearcut_map         !! Read a map prescribing whether a pxiel and PFT gets clearcut during spinup
!$OMP THREADPRIVATE(ok_read_sp_clearcut_map)

  LOGICAL, SAVE :: ok_change_species               !! A logical flag determining if we change species after a clearcut.
!$OMP THREADPRIVATE(ok_change_species)

  LOGICAL, SAVE :: ok_read_species_change_map      !! Read the new tree species from a species map
!$OMP THREADPRIVATE(ok_read_species_change_map)

  LOGICAL, SAVE :: ok_read_desired_fm_map          !! Read the new tree species from a species map
!$OMP THREADPRIVATE(ok_read_desired_fm_map)

  LOGICAL, SAVE :: ok_litter_raking                !! Activite litter raking
!$OMP THREADPRIVATE(ok_litter_raking)
 
  LOGICAL       :: ok_dimensional_product_use      !! Product pools based on the dimensions of the harvest 
!$OMP THREADPRIVATE(ok_dimensional_product_use)

  LOGICAL       :: ok_constant_mortality           !! Use constant mortality or calculate mortality as a function of last years´s NPP
!$OMP THREADPRIVATE(ok_constant_mortality)

  LOGICAL, SAVE :: ok_c13                          !! Activate carbon isotope concetration of biomass
!$OMP THREADPRIVATE(ok_c13)

  LOGICAL, SAVE :: ok_windthrow                    !! Activate the wind throw module.
!$OMP THREADPRIVATE(ok_windthrow)

  LOGICAL, SAVE :: ok_pest                         !! Activate pest  module. Trees will be killed by bark beetle
 !$OMP THREADPRIVATE(ok_pest)

  LOGICAL, SAVE :: ok_bare_soil_new                !! Flag that controls the view on and calculation of bare soil
!$OMP THREADPRIVATE(ok_bare_soil_new)

  CHARACTER(LEN=30), SAVE :: maint_resp_control    !! Choose an approach to calculate the maint respiration
!$OMP THREADPRIVATE(maint_resp_control) 

  LOGICAL, SAVE :: ok_force_pheno                  !! Use to force phenology when the conditions are not suitable
!$OMP THREADPRIVATE(ok_force_pheno)

 LOGICAL, SAVE :: moigdd_orchidee_2_0              !! Flag to use the ORCHIDEE 2.0 moigdd
                                                   !! phenology function
 !$OMP THREADPRIVATE(moigdd_orchidee_2_0) 

 LOGICAL, SAVE :: ok_crops_ipresenescence          !! Flag to allow crops to enter ipresenesence. Allows allocation to 'fruits'.
!$OMP THREADPRIVATE(ok_crops_ipresenescence)          

 LOGICAL, SAVE :: test_biomass_init                !! Use initial biomass estimate to determine size of plant at establishment.
!$OMP THREADPRIVATE(test_biomass_init)

 LOGICAL, SAVE :: ok_dyn_grass_den                 !! Use continuous dynamic grass density
!$OMP THREADPRIVATE(ok_dyn_grass_den)

  LOGICAL, SAVE :: nc_restart_compression          !! activate netcdf restart compression
!$OMP THREADPRIVATE(nc_restart_compression)

  CHARACTER(LEN=80), SAVE :: which_rt              !! Flag for radiative transfer: 'Iterative', 'Matrixial'
!$OMP THREADPRIVATE(which_rt)

  LOGICAL :: river_routing                         !! activate river routing
!$OMP THREADPRIVATE(river_routing)

  CHARACTER(LEN=80), SAVE :: routing_method        !! Character string used to switch between routing modules
!$OMP THREADPRIVATE(routing_method) 

  LOGICAL, SAVE :: ok_nudge_mc                     !! Activate nudging of soil moisture 
!$OMP THREADPRIVATE(ok_nudge_mc)

  LOGICAL, SAVE :: ok_nudge_snow                   !! Activate nudging of snow variables
!$OMP THREADPRIVATE(ok_nudge_snow)

  LOGICAL, SAVE :: nudge_interpol_with_xios        !! Activate reading and interpolation with XIOS for nudging fields
!$OMP THREADPRIVATE(nudge_interpol_with_xios)

  LOGICAL :: do_floodplains                        !! activate flood plains
!$OMP THREADPRIVATE(do_floodplains)

  LOGICAL :: do_irrigation                         !! activate computation of irrigation flux
!$OMP THREADPRIVATE(do_irrigation)

  LOGICAL :: ok_sechiba                            !! activate physic of the model
!$OMP THREADPRIVATE(ok_sechiba)

  LOGICAL :: ok_stomate                            !! activate carbon cycle
!$OMP THREADPRIVATE(ok_stomate)

  LOGICAL :: ok_ncycle                             !! activate nitrogen cycle
!$OMP THREADPRIVATE(ok_ncycle)

  LOGICAL :: impose_cn                             !! impose the CN ratio of leaves
!$OMP THREADPRIVATE(impose_cn)

  LOGICAL :: reset_impose_cn                       !! reset the CN ratio of leaves
!$OMP THREADPRIVATE(reset_impose_cn)

  LOGICAL :: read_cn                               !! read the CN ratio of leaves
!$OMP THREADPRIVATE(read_cn)

  LOGICAL :: ok_dgvm                               !! activate dynamic vegetation
!$OMP THREADPRIVATE(ok_dgvm)

  LOGICAL :: ok_lake_energy                        !! activate FLake
!$OMP THREADPRIVATE(ok_lake_energy)

  LOGICAL :: ok_snow_ice_lake_par                  !! activate FLake parametrization of snow and ice albedo
!$OMP THREADPRIVATE(ok_snow_ice_lake_par)

  LOGICAL :: ok_botsed_lake                        !! activate FLake parametrization of sediment layer
!$OMP THREADPRIVATE(ok_botsed_lake)

  LOGICAL :: do_wood_harvest                       !! activate wood harvest
!$OMP THREADPRIVATE(do_wood_harvest)

  LOGICAL :: ok_pheno                              !! activate the calculation of lai using stomate rather than a prescription
!$OMP THREADPRIVATE(ok_pheno)

  LOGICAL :: ok_bvoc                               !! activate biogenic volatile organic coumpounds
!$OMP THREADPRIVATE(ok_bvoc)

  LOGICAL :: ok_leafage                            !! activate leafage
!$OMP THREADPRIVATE(ok_leafage)

  LOGICAL :: ok_snags                              !! puts a part of woody litter input in a snag litter pool           
!$OMP THREADPRIVATE(ok_snags)

  LOGICAL :: ok_radcanopy                          !! use canopy radiative transfer model (BVOC module)
!$OMP THREADPRIVATE(ok_radcanopy)

  LOGICAL :: ok_multilayer                         !! use canopy radiative transfer model with multi-layers
!$OMP THREADPRIVATE(ok_multilayer)

  LOGICAL :: ok_pulse_NOx                          !! calculate NOx emissions with pulse
!$OMP THREADPRIVATE(ok_pulse_NOx)

  LOGICAL :: ok_bbgfertil_NOx                      !! calculate NOx emissions with bbg fertilizing effect
!$OMP THREADPRIVATE(ok_bbgfertil_NOx)

  LOGICAL :: ok_cropsfertil_NOx                    !! calculate NOx emissions with fertilizers use
!$OMP THREADPRIVATE(ok_cropsfertil_NOx)

  LOGICAL :: ok_co2bvoc_poss                       !! CO2 inhibition on isoprene activated following Possell et al. (2005) model
!$OMP THREADPRIVATE(ok_co2bvoc_poss)

  LOGICAL :: ok_co2bvoc_wilk                       !! CO2 inhibition on isoprene activated following Wilkinson et al. (2006) model
!$OMP THREADPRIVATE(ok_co2bvoc_wilk)

  LOGICAL :: use_fluxnet                           !! If true, disables retrieval of 10 m wind speed and gustiness multiplication
!$OMP THREADPRIVATE(use_fluxnet)
  LOGICAL :: cal_deflection                        !! If true, includes deflection load when calculating the critical wind speed in the wind module
!$OMP THREADPRIVATE(cal_deflection)
  LOGICAL :: wind_cal_pft                          !! If true, computes wind speed separately for each PFT based on canopy height
!$OMP THREADPRIVATE(wind_cal_pft)
  LOGICAL :: test_spacing                          !! If true, applies modified tree spacing in windthrow calculations
!$OMP THREADPRIVATE(test_spacing)


  CHARACTER(LEN=80), SAVE     :: restname_in       !! Input Restart files name for Sechiba component  
!$OMP THREADPRIVATE(restname_in)

  CHARACTER(LEN=80), SAVE :: restname_out          !! Output Restart files name for Sechiba component
!$OMP THREADPRIVATE(restname_out)

  CHARACTER(LEN=80), SAVE :: stom_restname_in      !! Input Restart files name for Stomate component
!$OMP THREADPRIVATE(stom_restname_in)

  CHARACTER(LEN=80), SAVE :: stom_restname_out     !! Output Restart files name for Stomate component
!$OMP THREADPRIVATE(stom_restname_out)

  !
  ! HACKS TO FIX AND REMOVE
  !
  ! Hacks are temporary solutions. Often quick and dirty that
  ! need attention before the problem could be considered to
  ! be solved.
  ! 
  LOGICAL, SAVE :: hack_enerbil_hydrol                    !! For debugging only!
!$OMP THREADPRIVATE(hack_enerbil_hydrol)

  REAL(r_std), SAVE :: hack_vessel_loss                   !! Constant vessel_loss in hydraulic_architecture
!$OMP THREADPRIVATE(hack_vessel_loss)
  
  REAL(r_std), SAVE :: rain_fac                           ! rain control factor
!$OMP THREADPRIVATE(rain_fac)

  REAL(r_std), SAVE :: refilling_capa                     ! refilling capacity, costless in terms of C
!$OMP THREADPRIVATE(refilling_capa)

  LOGICAL, SAVE :: hack_kb_m1                             !! Use higher threshold for small lai in the calculation of kB_m1 in condveg
!$OMP THREADPRIVATE(hack_kb_m1)


  !
  ! HACKS TO KEEP
  !
  ! These are hacks that do not fix an problem but that
  ! enhance the functionality of the code when trusting
  ! or debugging.
  !
  LOGICAL, SAVE :: trusting_hack_age_class                !! Skip the calculation of the biomass redistribution of the age classes
!$OMP THREADPRIVATE(trusting_hack_age_class)

  LOGICAL, SAVE :: hack_pgap                              !! For debugging only! If =.TRUE. the model uses Lambert Beer to calculate veget instead of Pgap.
!$OMP THREADPRIVATE(hack_pgap)

  LOGICAL, SAVE :: hack_veget_max_new                     !! Prescribe veget_max from run.def rather than a PFT map. Useful feature to debug simplified LCCs. 
                                                          !! See slowproc.f90 for details on this feature.
!$OMP THREADPRIVATE(hack_veget_max_new) 
  
  LOGICAL, SAVE       :: hack_can_structure               !! Overwrite the calculate canopy structure with a preset structure. Ony to be used for debugging !
!$OMP THREADPRIVATE(hack_can_structure)


  !
  ! TESTS
  !
  ! Tests are (final) solutions to a problem. Following testing,
  ! the code under the test_xxx flag can be accepted as final
  ! and the previous code can be removed.
  ! 
  LOGICAL, SAVE :: test_dynamic_alpha_self_thin          !! Use a dynamic parameter for alpha_selfthinning
  !$OMP THREADPRIVATE(test_dynamic_alpha_self_thin)
  LOGICAL, SAVE :: test_Q10_on_plant_uptake              !! Use a Q10 function to control the plant uptake.
!$OMP THREADPRIVATE(test_Q10_on_plant_uptake)
  LOGICAL, SAVE :: test_adjust_rhetero                   !! Adjust heterotrophic respiration when mineralisation becomes negative
!$OMP THREADPRIVATE(test_adjust_rhetero)

  LOGICAL, SAVE :: test_fdi_types                        !! Use vpd-based fire danger index.   
  !$OMP THREADPRIVATE(test_fdi_types)
  
  !
  ! PARAMETERS
  !
  REAL(r_std), SAVE :: min_n                             !! Minimum allowable n_mineralisation when truncating som_input_total(:,initrogen) in stomate_litter. 
!$OMP THREADPRIVATE(min_n)
  REAL(r_std), SAVE :: max_cn                            !! Maximum allowable ratio of som_input_total(:,icarbon) to som_input_total(:,initrogen).
!$OMP THREADPRIVATE(max_cn)
  REAL(r_std), SAVE :: snc                               !! Structural nitrogen  [gN g-1C] based on C:N of dead wood (White et al., 2000)
                                                         !! assuming carbon dry matter ratio of 0.5  
!$OMP THREADPRIVATE(snc) 
  REAL(r_std), SAVE :: slope_ra                          !! Reduction factor to make resp_maint less temperature sensitive
!$OMP THREADPRIVATE(slope_ra)
  INTEGER(i_std), SAVE :: spring_days_max                !! Maximum number of days during which we watch for possible spring frost damage  
!$OMP THREADPRIVATE(spring_days_max)
  REAL(r_std), SAVE :: alpha_nudge_mc                    !! Nudging constant for soil moisture 
!$OMP THREADPRIVATE(alpha_nudge_mc)
  REAL(r_std), SAVE :: alpha_nudge_snow                  !! Nudging constant for snow variables
!$OMP THREADPRIVATE(alpha_nudge_snow)
  REAL(r_std), SAVE :: spat_mod_self_thin                !! Spatial modifier of alpha_self_thinning during the spinup phase
!$OMP THREADPRIVATE(spat_mod_self_thin)
  REAL(r_std), SAVE :: spat_exp_self_thin                !! Spatial exponent of alpha_self_thinning during the spinup phase
!$OMP THREADPRIVATE(spat_exp_self_thin)
  
  !
  ! SPECIAL VALUES - These variables come without a getin_p. That is the reason
  ! why their default value is set here
  !
  INTEGER(i_std), PARAMETER :: undef_int = 999999999     !! undef integer for integer arrays (unitless)

  LOGICAL, SAVE :: OFF_LINE_MODE = .FALSE.               !! ORCHIDEE detects if it is coupled with a GCM or 
                                                         !! just use with one driver in OFF-LINE. (true/false)
!$OMP THREADPRIVATE(OFF_LINE_MODE)
  INTEGER, SAVE :: printlev=2                            !! Standard level for text output [0, 1, 2, 3]
!$OMP THREADPRIVATE(printlev)
  REAL(r_std), SAVE :: val_exp = 999999.                 !! Specific value if no restart value  (unitless)
!$OMP THREADPRIVATE(val_exp)
  REAL(r_std), PARAMETER :: undef = -9999.               !! Special value for stomate (unitless)
  
  REAL(r_std), PARAMETER :: min_sechiba = 1.E-8_r_std    !! Epsilon to detect a near zero floating point (unitless)
  REAL(r_std), PARAMETER :: undef_sechiba = 1.E+20_r_std !! The undef value used in SECHIBA (unitless)
  
  REAL(r_std), PARAMETER :: min_stomate = 1.E-8_r_std    !! Epsilon to detect a near zero floating point (unitless)
  REAL(r_std), PARAMETER :: large_value = 1.E33_r_std    !! some large value (for stomate) (unitless)

  
  !
  !  DIMENSIONING AND INDICES PARAMETERS  
  ! 
  INTEGER(i_std), PARAMETER :: istart = 1                !! Index to store values at the start
  INTEGER(i_std), PARAMETER :: iend = 2                  !! Index to store values at the end
  INTEGER(i_std), PARAMETER :: ibare_sechiba = 1         !! Index for bare soil in Sechiba (unitless)
  INTEGER(i_std), PARAMETER :: ivis = 1                  !! index for albedo in visible range (unitless)
  INTEGER(i_std), PARAMETER :: inir = 2                  !! index for albeod i near-infrared range (unitless) 
  INTEGER(i_std), PARAMETER :: n_spectralbands=2         !! number of spectral bands
  INTEGER(i_std), PARAMETER :: nnobio = 1                !! Number of other surface types: land ice (lakes,cities, ...) (unitless)
  INTEGER(i_std), PARAMETER :: iice = 1                  !! Index for land ice (see nnobio) (unitless)
  !-
  !! Soil
  INTEGER(i_std), PARAMETER :: classnb = 9               !! Levels of soil colour classification (unitless)
  !-
  INTEGER(i_std), PARAMETER :: nleafages = 4             !! leaf age discretisation ( 1 = no discretisation )(unitless)
  !-
  !! litter fractions: indices (unitless)
  INTEGER(i_std), PARAMETER :: ileaf = 1                 !! Index for leaf compartment (unitless)
  INTEGER(i_std), PARAMETER :: isapabove = 2             !! Index for sapwood above compartment (unitless)
  INTEGER(i_std), PARAMETER :: isapbelow = 3             !! Index for sapwood below compartment (unitless)
  INTEGER(i_std), PARAMETER :: iheartabove = 4           !! Index for heartwood above compartment (unitless)
  INTEGER(i_std), PARAMETER :: iheartbelow = 5           !! Index for heartwood below compartment (unitless)
  INTEGER(i_std), PARAMETER :: iroot = 6                 !! Index for roots compartment (unitless)
  INTEGER(i_std), PARAMETER :: ifruit = 7                !! Index for fruits compartment (unitless)
  INTEGER(i_std), PARAMETER :: icarbres = 8              !! Index for reserve compartment (unitless)
  INTEGER(i_std), PARAMETER :: ilabile = 9               !! Index for reserve compartment (unitless) 
  INTEGER(i_std), PARAMETER :: nparts = 9                !! Number of biomass compartments (unitless)
  !-
  !! indices for assimilation parameters 
  INTEGER(i_std), PARAMETER :: ivcmax = 1                !! Index for vcmax (assimilation parameters) (unitless)
  INTEGER(i_std), PARAMETER :: inue = 2                  !! Index for nue (assimilation parameters) (unitless)
  INTEGER(i_std), PARAMETER :: ileafN = 3                !! Index for leaf N (assimilation parameters) (unitless)
  INTEGER(i_std), PARAMETER :: npco2 = 3                 !! Number of assimilation parameters (unitless)

  !-
  !! trees and litter: indices for the parts of heart-
  !! and sapwood above and below the ground 
  INTEGER(i_std), PARAMETER :: iabove = 1               !! Index for above part (unitless)
  INTEGER(i_std), PARAMETER :: ibelow = 2               !! Index for below part (unitless)
  INTEGER(i_std), PARAMETER :: nlevs = 2                !! Number of levels for trees and litter (unitless)
  !-
  !! litter: indices for metabolic and structural part
  INTEGER(i_std), PARAMETER :: imetabolic = 1           !! Index for metabolic litter (unitless)
  INTEGER(i_std), PARAMETER :: istructural = 2          !! Index for structural litter (unitless)
  INTEGER(i_std), PARAMETER :: iwoody = 3               !! Index for woody litter (unitless)
  INTEGER(i_std), PARAMETER :: isnag = 4                !! Index for snag litter (unitless)
  INTEGER(i_std), PARAMETER :: nlitt = 4                !! Number of levels for litter compartments (unitless)
  !-
  !! litterfuel: indices for different types of fuel according to the scale of hours
  !! needed to reach equilibrium with the ambient moisture.
  INTEGER(i_std), PARAMETER :: ihour1 = 1       !! Index for metabolic litter (unitless)
  INTEGER(i_std), PARAMETER :: ihour10 = 2      !! Index for structural litter (unitless)
  INTEGER(i_std), PARAMETER :: ihour100 = 3     !! Index for woody litter (unitless)
  INTEGER(i_std), PARAMETER :: ihour1000 = 4    !! Index for woody litter (unitless)
  INTEGER(i_std), PARAMETER :: nhour = 4        !! Number of levels for litter compartments (unitless)
  !-
  !! ntypefire: indices indicating surface or crown fuel burning.
  INTEGER(i_std), PARAMETER :: ifiresurface = 1 !! Index for surface burning (unitless)
  INTEGER(i_std), PARAMETER :: ifirecrown = 2   !! Index for crown burning (unitless)
  INTEGER(i_std), PARAMETER :: nfiretype = 2    !! Number of vegetation fire types (unitless)
  !-

  !! carbon pools: indices
  INTEGER(i_std), PARAMETER :: iactive = 1              !! Index for active carbon pool (unitless)
  INTEGER(i_std), PARAMETER :: islow = 2                !! Index for slow carbon pool (unitless)
  INTEGER(i_std), PARAMETER :: ipassive = 3             !! Index for passive carbon pool (unitless)
  INTEGER(i_std), PARAMETER :: isurface = 4             !! Index for passive carbon pool (unitless)
  INTEGER(i_std), PARAMETER :: ncarb = 4                !! Number of soil carbon pools (unitless)
  !-
  !! For isotopes and nitrogen
  INTEGER(i_std), PARAMETER :: nelements = 2            !! Number of elements considered
  INTEGER(i_std), PARAMETER :: icarbon = 1              !! Index for carbon 
  INTEGER(i_std), PARAMETER :: initrogen = 2            !! Index for nitrogen 
  !! N-cycle : indices
  INTEGER(i_std), PARAMETER :: iammonium = 1            !! Index for Ammonium 
  INTEGER(i_std), PARAMETER :: initrate  = 2            !! Index for Nitrate
  INTEGER(i_std), PARAMETER :: inox      = 3            !! Index for NOX
  INTEGER(i_std), PARAMETER :: initrous  = 4            !! Index for N2O
  INTEGER(i_std), PARAMETER :: idinitro  = 5            !! Index for N2
  INTEGER(i_std), PARAMETER :: nnspec    = 5            !! Number of N-species considered
  INTEGER(i_std), PARAMETER :: nionspec  = 2            !! Number of ions form considered (ammonium, nitrate)
  
  !! For N-deposition
  INTEGER(i_std), PARAMETER :: iatm_ammo = 1            !! Index for N input from Ammonium N atmospheric deposition
  INTEGER(i_std), PARAMETER :: iatm_nitr = 2            !! Index for N input from Nitrate N atmospheric deposition
  INTEGER(i_std), PARAMETER :: ibnf      = 3            !! Index for N input from BNF
  INTEGER(i_std), PARAMETER :: imanure   = 4            !! Index for N input from Manure
  INTEGER(i_std), PARAMETER :: ifert_ammo = 5           !! Index for N input from Fertilisation
  INTEGER(i_std), PARAMETER :: ifert_nitr = 6           !! Index for N input from Fertilisation
  INTEGER(i_std), PARAMETER :: ninput    = 6            !! Number of N-input considered  
 

  INTEGER(i_std), PARAMETER :: i_nh4_to_no3 = 1         !! Index for NO3 production
  INTEGER(i_std), PARAMETER :: i_nh4_to_no  = 2         !! Index for NO production
  INTEGER(i_std), PARAMETER :: i_nh4_to_n2o = 3         !! Index for N2O production
  INTEGER(i_std), PARAMETER :: n_nh4_to_x = 3           !! Number of NH4 pathways

  INTEGER(i_std), PARAMETER :: i_no3_to_nox = 1         !! Index for NO3 consumption
  INTEGER(i_std), PARAMETER :: i_nox_to_n2o  = 2        !! Index for NO/Nox consumption
  INTEGER(i_std), PARAMETER :: i_n2o_to_n2 = 3          !! Index for N2O consumption
  INTEGER(i_std), PARAMETER :: n_n_to_x = 3             !! Number of N pathways

  INTEGER(i_std), PARAMETER :: nmonth = 12              !! Months in a year; used for input .nc files with monthly arrays
  
  !! Updates for the mass balance closure in stomate_lpj
  INTEGER(i_std), PARAMETER :: ibeg      = 1            !! At the begining of the routine
  INTEGER(i_std), PARAMETER :: ipre      = 2            !! After precribe
  INTEGER(i_std), PARAMETER :: iphe      = 3            !! After phenology
  INTEGER(i_std), PARAMETER :: igro      = 4            !! After growth functional allocation
  INTEGER(i_std), PARAMETER :: iage      = 5            !! After Age class distribution
  INTEGER(i_std), PARAMETER :: iluc      = 6            !! After land cover change
  INTEGER(i_std), PARAMETER :: idis      = 7            !! After land cover change after disturbances
  INTEGER(i_std), PARAMETER :: icle      = 8            !! After mortality_clean
  INTEGER(i_std), PARAMETER :: irec      = 9            !! After recruitment
  INTEGER(i_std), PARAMETER :: ispc      = 10           !! After Species Change
  INTEGER(i_std), PARAMETER :: nupdate1  = 10           !! Number of step in stomate_lpj where veget_max and atm_to_bm is update 
  INTEGER(i_std), PARAMETER :: ihar      = 11           !! After wood harvest
  INTEGER(i_std), PARAMETER :: imor      = 12           !! After mortality
  INTEGER(i_std), PARAMETER :: itur      = 13           !! After wood harvest
  INTEGER(i_std), PARAMETER :: ilpj      = 14           !! After all subroutines in lpj
  INTEGER(i_std), PARAMETER :: nupdate2  = 14           !! Number of step in stomate_lpj where wood volume is update

  ! These next sets of parameters are now used for both circ_class_kill and
  ! for the harvest_pool.  One source of confusion is what to do with trees that
  ! die from self-thinning or forest dieoffs.  These happen in all forests, regardless
  ! of management strategy.  I decided to put death of this kind into ifm_none, since
  ! it is the only type of mortality found in an unmanaged forest.  If the mortality
  ! does not kill the whole forest (e.g. self thinning), it goes into icut_thin.  If it
  ! does (forest dieoff), it goes into icut_clear.  The biomass is killed in lpj_gap.

  !! Indices used for forest management strategies
  INTEGER(i_std), PARAMETER :: nfm_types = 7             !! The total number of forest management 
                                                         !! strategies we can use
  INTEGER(i_std), PARAMETER :: ifm_none = 1              !! No human intervention in the forests.
  INTEGER(i_std), PARAMETER :: ifm_thin = 2              !! Regular thinning and harvesting of 
                                                         !! wood based on RDI.
  INTEGER(i_std), PARAMETER :: ifm_cop = 3               !! Coppicing for fuelwood.
  INTEGER(i_std), PARAMETER :: ifm_src = 4               !! Short rotation coppices for biomass 
                                                         !! production.
  INTEGER(i_std), PARAMETER :: ifm_uneven = 5            !! Uneven-aged management
  INTEGER(i_std), PARAMETER :: ifm_grass = 6             !! Grazing or cutting
  INTEGER(i_std), PARAMETER :: ifm_crop = 7              !! Crop harvest

  !! Parameter for thinning
  REAL(r_std), PARAMETER    :: r_max = 0.66              !! At r_max*diameter_max the thinning regime switches
                                                         !! from thinning from below to thinning from above.
  
  !! flag to trigger clearcut during spinup 
  INTEGER(i_std), PARAMETER :: flag_spinup_clearcut = 1  !! 

  !! Indices used for harvest pools
  INTEGER(i_std), PARAMETER :: ncut_times = 13           !! The total number of times when trees 
                                                         !! are cut and wood harvested.
  INTEGER(i_std), PARAMETER :: icut_clear = 1            !! A clearcut where all biomass is removed.
  INTEGER(i_std), PARAMETER :: icut_thin = 2             !! Thinning of biomass to reduce the 
                                                         !! number of trees.
  INTEGER(i_std), PARAMETER :: icut_lcc_wood = 3         !! Wood harvest following land cover 
                                                         !! change (LCC)
  INTEGER(i_std), PARAMETER :: icut_lcc_res = 4          !! Site clearing, removal of the stumps 
                                                         !! and branches following LCC
  INTEGER(i_std), PARAMETER :: icut_crop = 5             !! Crop harvest
  INTEGER(i_std), PARAMETER :: icut_grass = 6            !! Grazing or cutting
  INTEGER(i_std), PARAMETER :: icut_cop1 = 7             !! The first coppice cut
  INTEGER(i_std), PARAMETER :: icut_cop2 = 8             !! The second (and subsequent) coppice cut
  INTEGER(i_std), PARAMETER :: icut_cop3 = 9             !! The last coppice cut (only for SRC)
  INTEGER(i_std), PARAMETER :: icut_storm_break = 10     !! Stem breakage due to storm
  INTEGER(i_std), PARAMETER :: icut_storm_uproot = 11    !! Tee uprooting due to storm
  INTEGER(i_std), PARAMETER :: icut_beetle = 12          !! Tree killed due to bark beetle attacks
  INTEGER(i_std), PARAMETER :: icut_fire = 13            !! Tree killed due to fire disturbance

  !! Indices used to define the product pools
  !! Numbers based on Eggers 2008 - EFI report
  INTEGER(i_std), PARAMETER :: nshort = 1                !! Length in years of the short-lived product pool (GE 1)
  INTEGER(i_std), PARAMETER :: nmedium =17               !! Length in years of the medium-lived product pool (GT 4)
  INTEGER(i_std), PARAMETER :: nlong = 50                !! Length in years of the long-lived product pool (GT 4)

  !! Indices used for land use output variables
  INTEGER(i_std), PARAMETER :: nlanduse = 2              !! Total number of different land uses
  INTEGER(i_std), PARAMETER :: iharvest = 1              !! Biomass removals due to all sorts of harvest
  INTEGER(i_std), PARAMETER :: ilcc = 2                  !! Biomass removals due to all sorts of land cover changes

  !! Indices used to check the mass balance closure
  INTEGER(i_std), PARAMETER :: nmbcomp = 5               !! The total number of components in 
                                                         !! our mass balance check
  INTEGER(i_std), PARAMETER :: iatm2land = 1             !! atmosphere to land fluxes such as GPP 
                                                         !! and co2_2_bm
  INTEGER(i_std), PARAMETER :: iland2atm = 2             !! land to atmosphere fluxes such as Rh, 
                                                         !! Ra and product decomposition
  INTEGER(i_std), PARAMETER :: ilat2out = 3              !! outgoing lateral flux i.e. DOC leaching 
                                                         !! for the litter routine
  INTEGER(i_std), PARAMETER :: ilat2in = 4               !! incoming lateral flux i.e. N deposition 
                                                         !! for the land
  INTEGER(i_std), PARAMETER :: ipoolchange = 5           !! change in pool size i.e. change in biomass

  !! Indices used for warning tracking
  INTEGER(i_std), PARAMETER :: nwarns = 1                !! The total number of warnings we track
  INTEGER(i_std), PARAMETER :: iwphoto = 1               !! A warning about division by zero in photosynthesis

  !! Indices used for wind damage  
  INTEGER(i_std), PARAMETER :: ibreakage = 1             !! The index for stem breakage dur to wind damage
  INTEGER(i_std), PARAMETER :: ioverturning = 2          !! The index for the tree overtuning due to wind damage

  !! Indices for orphan fluxes
  INTEGER(i_std), PARAMETER :: norphans = 8              !! Total number of orphan fluxes (unitless)
  INTEGER(i_std), PARAMETER :: ivegold = 1               !! Index for veget_max before LCC
  INTEGER(i_std), PARAMETER :: ivegnew = 2               !! Index for veget_max before LCC (includes veget_max of orphan fluxes)
  INTEGER(i_std), PARAMETER :: igpp = 3                  !! Index for gpp_daily
  INTEGER(i_std), PARAMETER :: ico2bm = 4                !! Index for co2_to_bm
  INTEGER(i_std), PARAMETER :: irmain = 5                !! Index for maintenance respiration
  INTEGER(i_std), PARAMETER :: irgrow = 6                !! Index for growth respiration 
  INTEGER(i_std), PARAMETER :: inpp = 7                  !! Index for npp_daily
  INTEGER(i_std), PARAMETER :: irhet = 8                 !! Index for total heterotrophic respiration
  !

  !! Indices for phenology
  ! The variable ::plant_status replaces several variables (senescence, 
  ! begin_leaves and allow_phenoinit) that describe the phenological status of the plant
  ! by storing these different aspects of phenology in a single variable, inconsistencies
  ! become impossible or at least easier to check.
  ! When the model starts from scratch the status is set to iprescribe this allows us
  ! to grow leaves from the first year onwards. The plant should then go through the 
  ! different growth phases: ibudsavail, ibudbreak, icanopy, isenescent, idormant and 
  ! finally idead. Following idormant the status should return ibudsavail to initiate 
  ! another cycle in the subsequent growing season. Following idead new vegetation 
  ! should be prescribed.
  INTEGER(i_std), PARAMETER :: inone = 0                 !! No plants thus no status
  INTEGER(i_std), PARAMETER :: iprescribe = 1            !! Prescribe a PFT
  INTEGER(i_std), PARAMETER :: ibudsavail = 2            !! Buds are present
  INTEGER(i_std), PARAMETER :: ibudbreak = 3             !! Day that the buds break and leaf 
                                                         !! on-set begins
  INTEGER(i_std), PARAMETER :: icanopy = 4               !! Canopy is present
  INTEGER(i_std), PARAMETER :: ipresenescence = 5        !! Wood growth cessation
  INTEGER(i_std), PARAMETER :: isenescent = 6            !! The plant is senescent
  INTEGER(i_std), PARAMETER :: idormant = 7              !! The plant is dormant
  INTEGER(i_std), PARAMETER :: idead = 8                 !! The plant was killed

  !
  !! Indices used for analytical spin-up
  INTEGER(i_std), PARAMETER :: nbpools = 12             !! Total number of carbon pools (unitless)
  INTEGER(i_std), PARAMETER :: istructural_above = 1    !! Index for structural litter above (unitless)
  INTEGER(i_std), PARAMETER :: istructural_below = 2    !! Index for structural litter below (unitless)
  INTEGER(i_std), PARAMETER :: imetabolic_above = 3     !! Index for metabolic litter above (unitless)
  INTEGER(i_std), PARAMETER :: imetabolic_below = 4     !! Index for metabolic litter below (unitless)
  INTEGER(i_std), PARAMETER :: iwoody_above = 5         !! Index for woody litter above (unitless)
  INTEGER(i_std), PARAMETER :: iwoody_below = 6         !! Index for woody litter below (unitless)
  INTEGER(i_std), PARAMETER :: isnag_above = 7          !! Index for snag litter above (unitless)
  INTEGER(i_std), PARAMETER :: isnag_below = 8          !! Index for snag litter below (unitless)
  INTEGER(i_std), PARAMETER :: iactive_pool = 9         !! Index for active carbon pool (unitless)
  INTEGER(i_std), PARAMETER :: islow_pool   = 10        !! Index for slow carbon pool (unitless)
  INTEGER(i_std), PARAMETER :: ipassive_pool = 11       !! Index for passive carbon pool (unitless)
  INTEGER(i_std), PARAMETER :: isurface_pool = 12       !! Index for surface carbon pool (unitless)

  !
  !! Indices for the burried carbon and nitrogen following LUC
  INTEGER(i_std), PARAMETER :: nbury=6                  !! Total number of buried pools
  INTEGER(i_std), PARAMETER :: ilitter=1                !! Index for litter
  INTEGER(i_std), PARAMETER :: ifreshlitter=2           !! Index for fresh litter
  INTEGER(i_std), PARAMETER :: ifreshsom=3              !! Index for fresh som
  INTEGER(i_std), PARAMETER :: ibact=4                  !! Index for bacteria
  INTEGER(i_std), PARAMETER :: isom=5                   !! Index for som
  INTEGER(i_std), PARAMETER :: iminnitrogen=6           !! Index for mineral soil nitrogen

  !
  !! Indicies used for output variables on Landuse tiles defined according to LUMIP project
  !! Note that ORCHIDEE do not represent pasture and urban land. Therefor the variables will have 
  !! val_exp as missing value for these tiles. 
  INTEGER(i_std), PARAMETER :: nlut=4                   !! Total number of landuse tiles according to LUMIP
  INTEGER(i_std), PARAMETER :: id_psl=1                 !! Index for primary and secondary land
  INTEGER(i_std), PARAMETER :: id_pst=2                 !! Index for pasture land
  INTEGER(i_std), PARAMETER :: id_crp=3                 !! Index for crop land
  INTEGER(i_std), PARAMETER :: id_urb=4                 !! Index for urban land

  !! Indices used for canopy structure (Pgap & eff lai) 
  INTEGER(i_std),PARAMETER   :: ndist_types=6           !! the number of distributions we need in the LAI effective routines
  INTEGER(i_std),PARAMETER   :: iheight=1               !! the tree height distribution
  INTEGER(i_std),PARAMETER   :: idiameter=2             !! the trunk diameter distribution
  INTEGER(i_std),PARAMETER   :: icnvol=3                !! the crown volume distribution
  INTEGER(i_std),PARAMETER   :: icnarea=4               !! the crown area distribution
  INTEGER(i_std),PARAMETER   :: icndiaver=5             !! the verticle crown diameter distribution
  INTEGER(i_std),PARAMETER   :: icndiahor=6             !! the horizontal crown diameter distribution

  !! indices used in AR5 output
  INTEGER(i_std),PARAMETER   :: nlctypes=3              !! Number of land cover types
  INTEGER(i_std),PARAMETER   :: iforest=1               !! forest
  INTEGER(i_std),PARAMETER   :: igrass=2                !! grass
  INTEGER(i_std),PARAMETER   :: icrop=3                 !! crops

  !! Indices used in the root profile
  INTEGER(i_std),PARAMETER   :: nroot_prof=2            !! Number of different root profiles
  INTEGER(i_std),PARAMETER   :: istruc=1                !! Structural root profile
  INTEGER(i_std),PARAMETER   :: ifunc=2                 !! Functional root profile for water uptake

  !! Indices used to defile the vertical root rofile
  INTEGER(i_std),PARAMETER   :: ndepths=2               !! Number of components in the definition of a vertical soil profile
  INTEGER(i_std),PARAMETER   :: inode=1                 !! Node
  INTEGER(i_std),PARAMETER   :: iinterface=2            !! Interface


  !
  ! NUMERICAL AND PHYSICS CONSTANTS
  !
  !-
  ! 1. Mathematical and numerical constants
  !-
  REAL(r_std), PARAMETER :: pi = 3.141592653589793238   !! pi souce : http://mathworld.wolfram.com/Pi.html (unitless)
  REAL(r_std), PARAMETER :: euler = 2.71828182845904523 !! e source : http://mathworld.wolfram.com/e.html (unitless)
  REAL(r_std), PARAMETER :: zero = 0._r_std             !! Numerical constant set to 0 (unitless)
  REAL(r_std), PARAMETER :: undemi = 0.5_r_std          !! Numerical constant set to 1/2 (unitless)
  REAL(r_std), PARAMETER :: un = 1._r_std               !! Numerical constant set to 1 (unitless)
  REAL(r_std), PARAMETER :: moins_un = -1._r_std        !! Numerical constant set to -1 (unitless)
  REAL(r_std), PARAMETER :: deux = 2._r_std             !! Numerical constant set to 2 (unitless)
  REAL(r_std), PARAMETER :: trois = 3._r_std            !! Numerical constant set to 3 (unitless)
  REAL(r_std), PARAMETER :: quatre = 4._r_std           !! Numerical constant set to 4 (unitless)
  REAL(r_std), PARAMETER :: cinq = 5._r_std             !![DISPENSABLE] Numerical constant set to 5 (unitless)
  REAL(r_std), PARAMETER :: six = 6._r_std              !![DISPENSABLE] Numerical constant set to 6 (unitless)
  REAL(r_std), PARAMETER :: huit = 8._r_std             !! Numerical constant set to 8 (unitless)
  REAL(r_std), PARAMETER :: mille = 1000._r_std         !! Numerical constant set to 1000 (unitless)

  !-
  ! 2 . Physics - these are constantes that cannot be overwritten by a getin_p. Therefore
  !     their default value is set here.
  !-
  REAL(r_std), PARAMETER :: R_Earth = 6378000.              !! radius of the Earth : Earth radius ~= Equatorial radius (m)
  REAL(r_std), PARAMETER :: mincos  = 0.0001                !! Minimum cosine value used for interpolation (unitless) 
  REAL(r_std), PARAMETER :: pb_std = 1013.                  !! standard pressure (hPa)
  REAL(r_std), PARAMETER :: ZeroCelsius = 273.15            !! 0 degre Celsius in degre Kelvin (K)
  REAL(r_std), PARAMETER :: tp_00 = 273.15                  !! 0 degre Celsius in degre Kelvin (K)
  REAL(r_std), PARAMETER :: chalsu0 = 2.8345E06             !! Latent heat of sublimation (J.kg^{-1})
  REAL(r_std), PARAMETER :: chalev0 = 2.5008E06             !! Latent heat of evaporation (J.kg^{-1}) 
  REAL(r_std), PARAMETER :: chalfu0 = chalsu0-chalev0       !! Latent heat of fusion (J.kg^{-1}) 
  REAL(r_std), PARAMETER :: c_stefan = 5.6697E-8            !! Stefan-Boltzman constant (W.m^{-2}.K^{-4})
  REAL(r_std), PARAMETER :: cp_air = 1004.675               !! Specific heat of dry air (J.kg^{-1}.K^{-1}) 
  REAL(r_std), PARAMETER :: cte_molr = 287.05               !! Specific constant of dry air (kg.mol^{-1}) 
  REAL(r_std), PARAMETER :: kappa = cte_molr/cp_air         !! Kappa : ratio between specific constant and specific heat 
                                                            !! of dry air (unitless)
  REAL(r_std), PARAMETER :: msmlr_air = 28.964E-03          !! Molecular weight of dry air (kg.mol^{-1})
  REAL(r_std), PARAMETER :: msmlr_h2o = 18.02E-03           !! Molecular weight of water vapor (kg.mol^{-1}) 
  REAL(r_std), PARAMETER :: cp_h2o = &                      !! Specific heat of water vapor (J.kg^{-1}.K^{-1}) 
       & cp_air*(quatre*msmlr_air)/( 3.5_r_std*msmlr_h2o) 
  REAL(r_std), PARAMETER :: cte_molr_h2o = cte_molr/quatre  !! Specific constant of water vapor (J.kg^{-1}.K^{-1}) 
  REAL(r_std), PARAMETER :: retv = msmlr_air/msmlr_h2o-un   !! Ratio between molecular weight of dry air and water 
                                                            !! vapor minus 1(unitless)  
  REAL(r_std), PARAMETER :: rvtmp2 = cp_h2o/cp_air-un       !! Ratio between specific heat of water vapor and dry air
  
  REAL(r_std), PARAMETER :: rho_h2o= 0.9991_r_std           !! Density of water at 15°C (g cm-3)     
  REAL(r_std), PARAMETER :: cepdu2 = (0.1_r_std)**2         !! Squared wind shear (m^2.s^{-2}) 
  REAL(r_std), PARAMETER :: ct_karman = 0.40_r_std          !! Van Karmann Constant (unitless)
  REAL(r_std), PARAMETER :: cte_grav = 9.80665_r_std        !! Acceleration of the gravity (m.s^{-2})
  REAL(r_std), PARAMETER :: pa_par_hpa = 100._r_std         !! Transform pascal into hectopascal (unitless)
  REAL(r_std), PARAMETER :: RR = 8.314                      !! Ideal gas constant (J.mol^{-1}.K^{-1})
  REAL(r_std), PARAMETER :: Sct = 1370.                     !! Solar constant (W.m^{-2}) 

  REAL(r_std), PARAMETER :: omega_earth = 7.2921E-05        !! Rotation rate of the Earth (rad.s^{-1})
  REAL(r_std), PARAMETER :: mm_m = 1000._r_std              !! conversion from milimeters to meters

  !-
  ! 3. Climatic constants
  !-
  !! Constantes of the Louis scheme 
  REAL(r_std), SAVE :: cb                         !! Constant of the Louis scheme (unitless);
                                                  !! reference to Louis (1979)
!$OMP THREADPRIVATE(cb)
  REAL(r_std), SAVE :: cc                         !! Constant of the Louis scheme (unitless);
                                                  !! reference to Louis (1979)
!$OMP THREADPRIVATE(cc)
  REAL(r_std), SAVE :: cd                         !! Constant of the Louis scheme (unitless);
                                                  !! reference to Louis (1979)
!$OMP THREADPRIVATE(cd)
  REAL(r_std), SAVE :: rayt_cste                  !! Constant in the computation of surface resistance (W.m^{-2})
!$OMP THREADPRIVATE(rayt_cste)
  REAL(r_std), SAVE :: defc_plus                  !! Constant in the computation of surface resistance (K.W^{-1})
!$OMP THREADPRIVATE(defc_plus)
  REAL(r_std), SAVE :: defc_mult                  !! Constant in the computation of surface resistance (K.W^{-1})
!$OMP THREADPRIVATE(defc_mult)

  !-
  ! 4. Soil thermodynamics constants
  !-
  ! Look at constantes_soil.f90

  !
  ! OPTIONAL PARTS OF THE MODEL
  !
  LOGICAL, PARAMETER     :: diag_qsat=.TRUE.      !! CHECK - should this be externalized? One of the most frequent problems is a temperature out of range
                                                  !! we provide here a way to catch that in the calling procedure. 
                                                  !! (from Jan Polcher)(true/false) 
  LOGICAL, SAVE     :: almaoutput                 !! Selects the type of output for the model.(true/false)
                                                  !! Value is read from run.def in intersurf_history
!$OMP THREADPRIVATE(almaoutput)

  !
  ! DIVERSE
  !
  CHARACTER(LEN=100), SAVE :: stomate_forcing_name='NONE'  !! Name of STOMATE forcing file (unitless)
                                                           ! Compatibility with Nicolas Viovy driver.
!$OMP THREADPRIVATE(stomate_forcing_name)
  CHARACTER(LEN=100), SAVE :: stomate_Cforcing_name='NONE' !! Name of soil forcing file (unitless)
                                                           ! Compatibility with Nicolas Viovy driver.
!$OMP THREADPRIVATE(stomate_Cforcing_name)
  CHARACTER(LEN=100), SAVE :: stomate_Cforcing_discretization_name='NONE' !! Name of soil carbon discretization forcing file (unitless)
!$OMP THREADPRIVATE(stomate_Cforcing_discretization_name)
  INTEGER(i_std), SAVE :: forcing_id                 !! Index of the forcing file (unitless)
!$OMP THREADPRIVATE(forcing_id)
  LOGICAL, SAVE :: allow_forcing_write=.TRUE.        !! Allow writing of stomate_forcing file. 
                                                     !! This variable will be set to false for teststomate.
!$OMP THREADPRIVATE(allow_forcing_write)
  LOGICAL, SAVE :: ok_sort_lon_lat                   !! Activate sort in increasing order of latitudes and longitudes read from the forcing file
!$OMP THREADPRIVATE(ok_sort_lon_lat)



                         !------------------------!
                         !  SECHIBA PARAMETERS    !
                         !------------------------!
 

  !
  ! GLOBAL PARAMETERS   
  !
  REAL(r_std), SAVE :: min_wind            !! The minimum wind (m.s^{-1})
!$OMP THREADPRIVATE(min_wind)
  REAL(r_std), SAVE :: min_qc              !! The minimum value for qc (qc=drag*wind) used in coupled(enerbil) and forced mode (enerbil and diffuco) 
  REAL(r_std), SAVE :: snowcri             !! Sets the amount above which only sublimation occurs (kg.m^{-2})
!$OMP THREADPRIVATE(snowcri)


  !
  ! FLAGS ACTIVATING SUB-MODELS
  !
  LOGICAL, SAVE :: treat_expansion                     !! Do we treat PFT expansion across a grid point after introduction? (true/false)
!$OMP THREADPRIVATE(treat_expansion)
  LOGICAL, SAVE :: ok_herbivores                       !! flag to activate herbivores (true/false)
!$OMP THREADPRIVATE(ok_herbivores)
  LOGICAL, SAVE :: harvest_agri                        !! flag to harvest aboveground biomass from agricultural PFTs)(true/false)
!$OMP THREADPRIVATE(harvest_agri)
  LOGICAL, SAVE :: lpj_gap_const_mort                  !! constant moratlity (true/false). Default value depend on OK_DGVM.
!$OMP THREADPRIVATE(lpj_gap_const_mort)
  LOGICAL, SAVE :: ok_spitfire                         !! Activate fires with spitfire
!$OMP THREADPRIVATE(ok_spitfire)
  LOGICAL, SAVE :: ok_lightningfire                    !! Activate lightning fires with spitfire
!$OMP THREADPRIVATE(ok_lightningfire)
  CHARACTER(LEN=100), SAVE :: road_length_file           !! FRAG_3TERMS macro: road-length NetCDF; 'NONE' = term off
!$OMP THREADPRIVATE(road_length_file)
  REAL(r_std), SAVE :: road_grip_correction              !! FRAG_3TERMS macro: GRIP under-reporting factor (Bowring et al. 2024: 2.0)
!$OMP THREADPRIVATE(road_grip_correction)
  LOGICAL, SAVE :: ok_road_forest_mask                   !! FRAG_3TERMS macro: T = map already restricted to forest (5 arcmin EFDA);
                                                         !! F = map is total road length, apply the model's own forest fraction
!$OMP THREADPRIVATE(ok_road_forest_mask)
  REAL(r_std), SAVE :: edge_grain_default                !! FRAG_3TERMS: fallback sub-grid fragmentation grain (m)
!$OMP THREADPRIVATE(edge_grain_default)
  LOGICAL, SAVE :: ok_aed_wind                         !! flag that allows fragmentation wind effect (true/false)
!$OMP THREADPRIVATE(ok_aed_wind)
  LOGICAL, SAVE :: ok_aed_size                         !! flag that allows fragmentation fire size effect (true/false)
!$OMP THREADPRIVATE(ok_aed_size)
  LOGICAL, SAVE :: ok_aed_humign                       !! flag that allows fragmentation human ignition effect (true/false)
!$OMP THREADPRIVATE(ok_aed_humign)
  LOGICAL, SAVE :: ok_aed_fuel                         !! flag that allows fragmentation fuel wetness effect (true/false)
!$OMP THREADPRIVATE(ok_aed_fuel)
  LOGICAL, SAVE :: ok_aed_light                        !! flag that allows fragmentation light/canopy effect (true/false)
!$OMP THREADPRIVATE(ok_aed_light)
  LOGICAL, SAVE :: ok_aed_feedback                     !! flag for perturbation -> edge_length feedback (Marie 2026, true/false)
!$OMP THREADPRIVATE(ok_aed_feedback)
  LOGICAL, SAVE :: ok_edge_from_age_class                     !! edge length derived from the young age-class open area A_open, instead of a tau_rec relaxation
!$OMP THREADPRIVATE(ok_edge_from_age_class)
  LOGICAL, SAVE :: stand_age_bm_in_restart = .FALSE.   !! STAND_AGE: TRUE when age_stand_bm was present in the restart -> cold-start fallback to age_stand must NOT override it
!$OMP THREADPRIVATE(stand_age_bm_in_restart)
  REAL(r_std), SAVE :: dia_rotation_tol                !! AGE_ROTATION: lower tolerance on the clearcut diameter (fraction, ~0.15-0.20)
!$OMP THREADPRIVATE(dia_rotation_tol)
  LOGICAL, SAVE :: ok_progressive_harvest              !! PROGRESSIVE_HARVEST: clearcut a rotation-paced FRACTION of the mature slot
                                                       !! instead of emptying it. Default FALSE = all-or-nothing (bit-neutral)
!$OMP THREADPRIVATE(ok_progressive_harvest)
  REAL(r_std), SAVE :: cc_n_floor           !! Effectif sous lequel une classe de circonference
                                            !! n'existe pas (ind m-2). Sa biomasse TOTALE et ses tiges
                                            !! descendent d'une classe. Defaut 0 => historique (bit-neutre).
!$OMP THREADPRIVATE(cc_n_floor)
  LOGICAL, SAVE :: ok_rdi_ramp                         !! RDI consigne DESCEND pendant la classe d'amelioration (rampe intra-classe)
!$OMP THREADPRIVATE(ok_rdi_ramp)
  LOGICAL, SAVE :: ok_forced_thinning                  !! Eclaircie ANNUELLE forcee sur les classes de maturite jeunes
!$OMP THREADPRIVATE(ok_forced_thinning)
  INTEGER(i_std), SAVE :: forced_thin_last_class       !! Derniere classe de maturite soumise a l'eclaircie forcee
!$OMP THREADPRIVATE(forced_thin_last_class)
  LOGICAL, SAVE :: ok_min_density_rdi                  !! CONDITION 1 judges collapse on RDI, not on a fixed stem count:
                                                       !! a stand thinned to its target is sparse BY DESIGN, not collapsed
!$OMP THREADPRIVATE(ok_min_density_rdi)
  REAL(r_std), SAVE :: min_density_rdi_floor           !! RDI below which a stand is held to have collapsed (-)
!$OMP THREADPRIVATE(min_density_rdi_floor)
  LOGICAL, SAVE :: ok_coppice_2class                   !! COPPICE_2CLASS: a coppice has two maturity classes, not four, and
                                                       !! graduates on HEIGHT. Default FALSE = current behaviour (bit-neutral)
!$OMP THREADPRIVATE(ok_coppice_2class)
  REAL(r_std), SAVE :: cop_r_up                        !! COPPICE_2CLASS: height ratio above which class 1 feeds the terminal class (-)
!$OMP THREADPRIVATE(cop_r_up)
  REAL(r_std), SAVE :: cop_r_down                      !! COPPICE_2CLASS: height ratio below which a slot falls back to class 1 (-)
!$OMP THREADPRIVATE(cop_r_down)
  REAL(r_std), SAVE :: age_stand_estab                 !! Age d'un peuplement fraichement etabli (an), plancher de age_stand_bm
!$OMP THREADPRIVATE(age_stand_estab)
  REAL(r_std), SAVE :: aed_mix_tau                     !! AED_REGROWTH v2: memory (yr) of the leaky agent-flux integrator (agent mix)
!$OMP THREADPRIVATE(aed_mix_tau)
  ! Guillaume M. -- A_open carries the regeneration gap ONLY. The four opening weights
  ! AED_OPEN_W1..W4 and the landscape term AED_ALPHA_FRAG both described the mature-forest
  ! edge, which the meso term produces from the landscape geometry (EL_meso = A*4*f*(1-f)/r).
  ! Keeping either alongside the meso term counts the same interface twice.
  REAL(r_std), SAVE :: aed_alpha_reg                   !! Weight of the regeneration gap (age class 1), the only true internal gap
!$OMP THREADPRIVATE(aed_alpha_reg)
  REAL(r_std), SAVE :: aed_bg_frac                     !! AED_REGROWTH v2: background-edge fraction for young/open forest with no current disturbance flux (piste B)
!$OMP THREADPRIVATE(aed_bg_frac)
  LOGICAL, SAVE :: spinup_analytic                     !! Flag to activate analytical resolution for spinup (true/false)
!$OMP THREADPRIVATE(spinup_analytic)
  LOGICAL, SAVE :: calculate_gpp_preind                !! Flag to activate calculation of gpp longterm preindustrial reference values
!$OMP THREADPRIVATE(calculate_gpp_preind)
  LOGICAL, SAVE :: ok_soil_carbon_discretization       !! Flag to activate soil carbon discretization (vertical carbon and soil carbon thermal insulation)
!$OMP THREADPRIVATE(ok_soil_carbon_discretization)
  LOGICAL, SAVE :: new_carbinput_intdepzlit            !! Flag to activate a different method to redistribute the litter input into the different soil layers.
!$OMP THREADPRIVATE(new_carbinput_intdepzlit)
  LOGICAL, SAVE :: ok_forcesoil_write                  !! Flag to activate soil carbon discretization output file to be used for forcesoil program
!$OMP THREADPRIVATE(ok_forcesoil_write)
  LOGICAL, SAVE :: ok_vessel_mortality                 !! Flag to activate death and recovery of vegetation following hydraulic failure.
!$OMP THREADPRIVATE(ok_vessel_mortality)
  LOGICAL, SAVE :: ok_sap_feedback                     !! Flag to activate sapwood mortality feedback on sapwood resistance in Tuzet
!$OMP THREADPRIVATE(ok_sap_feedback)
  LOGICAL, SAVE :: no_g0                               !! Flag that deactivates leaf residual conductance
!$OMP THREADPRIVATE(no_g0)  
  LOGICAL, SAVE :: test_functional_roots_tuzet         !! Flag that uses functional roots profile in tuzet
!$OMP THREADPRIVATE(test_functional_roots_tuzet)  
  LOGICAL, SAVE :: ok_nsl_vpd                          !! Flag that activates VPD control over gs (in addition to water stress control)
!$OMP THREADPRIVATE(ok_nsl_vpd)   
  LOGICAL, SAVE :: ok_tuzet_indiv                      !! Flag that activates Tuzet at indiv scale
!$OMP THREADPRIVATE(ok_tuzet_indiv)  
  LOGICAL, SAVE :: tuzet_inf_layer                     !! Flag that triggers using the inferior layer as base of plant water potentials
!$OMP THREADPRIVATE(tuzet_inf_layer)
  LOGICAL, SAVE :: ok_ice_sheet                        !! Flag to activate ice scheme on ice-sheet area
!$OMP THREADPRIVATE(ok_ice_sheet)
  LOGICAL, SAVE :: use_thermo_newdisc                  !! Activate new equations (NewDisC) for capacity and conductivities
!$OMP THREADPRIVATE(use_thermo_newdisc)
  LOGICAL, SAVE :: impose_soc_bd                       !! Impose soil carbon(SOC), bulk density(BDOD) and coarse matter fraction(CFVO) in thermosoil for site simulation
!$OMP THREADPRIVATE(impose_soc_bd)


  !
  ! CONFIGURATION VEGETATION
  !
  LOGICAL, SAVE :: agriculture             !! allow agricultural PFTs (true/false)
!$OMP THREADPRIVATE(agriculture)
  LOGICAL, SAVE :: impveg                  !! Impose vegetation ? (true/false)
!$OMP THREADPRIVATE(impveg)
  LOGICAL, SAVE :: impsoilt                !! Impose soil ? (true/false)
!$OMP THREADPRIVATE(impsoilt)
  LOGICAL, SAVE :: impslope                !! Impose reinf_slope ? (true/false)
!$OMP THREADPRIVATE(impslope)
  LOGICAL, SAVE :: impose_ninput_dep       !! Impose N input values from the atmosphere ? (true/false)
!$OMP THREADPRIVATE(impose_ninput_dep)
  LOGICAL, SAVE :: impose_ninput_fert      !! Impose N input values from fertilizer ? (true/false)        
!$OMP THREADPRIVATE(impose_ninput_fert)    
  LOGICAL, SAVE :: impose_ninput_manure    !! Impose N input values from manure? (true/false)        
!$OMP THREADPRIVATE(impose_ninput_manure)                     
  LOGICAL, SAVE :: impose_ninput_bnf       !! Impose N input values from biological nitrogen fixation (BNF) ? (true/false)        
!$OMP THREADPRIVATE(impose_ninput_bnf) 
  INTEGER, SAVE :: read_lai                !! Option for LAI if STOMATE is not activated (0: no LAI reading, 1: site-based LAI reading, 2: read LAI map)
!$OMP THREADPRIVATE(read_lai)
  LOGICAL, SAVE :: vegetmap_reset          !! Reset the vegetation map and reset carbon related variables
!$OMP THREADPRIVATE(vegetmap_reset)
  INTEGER(i_std) , SAVE :: veget_update    !! Update frequency in years for landuse (nb of years)
!$OMP THREADPRIVATE(veget_update)
  LOGICAL, SAVE :: ninput_reinit           !! To change N INPUT file in a run. (true/false)
!$OMP THREADPRIVATE(ninput_reinit)
  
  !
  ! PARAMETERS CONTROLED BY THE TIME STEP IN THE CODE
  ! NO GETIN NEEDED - INITIALIZED HERE
  !
  LOGICAL, SAVE :: do_now_stomate_lcchange = .FALSE.    !! Time to call lcchange in stomate_lpj
!$OMP THREADPRIVATE(do_now_stomate_lcchange)
  LOGICAL, SAVE :: do_now_stomate_woodharvest = .FALSE. !! Time to call woodharvest in stomate_lpj
!$OMP THREADPRIVATE(do_now_stomate_woodharvest)
  LOGICAL, SAVE :: done_stomate_lcchange = .FALSE.      !! If true, call lcchange in stomate_lpj has just been done. 
!$OMP THREADPRIVATE(done_stomate_lcchange)
  LOGICAL, SAVE :: do_now_recruit = .FALSE.             !! Time to call prescribe to calculate recruits in stomate_lpj
!$OMP THREADPRIVATE(do_now_recruit)
 
  !
  ! PARAMETERS USED BY BOTH HYDROLOGY MODELS
  !
  REAL(r_std), SAVE :: max_snow_age                     !! Maximum period of snow aging (days)
!$OMP THREADPRIVATE(max_snow_age)
  REAL(r_std), SAVE :: snow_trans                       !! Transformation time constant for snow (m), reduced from the value 0.3 (04/07/2016)
!$OMP THREADPRIVATE(snow_trans)
  REAL(r_std), SAVE :: snow_trans_nobio                 !! Transformation time constant for snow on nobio (m), (18/02/2019)
!$OMP THREADPRIVATE(snow_trans_nobio)
  REAL(r_std), SAVE :: snowa_aged_nobio_vis             !! Minimum snow albedo value for nobio type after aging (dirty old snow), visible albedo (09/12/2022)
!$OMP THREADPRIVATE(snowa_aged_nobio_vis)
  REAL(r_std), SAVE :: snowa_aged_nobio_nir             !! Minimum snow albedo value for nobio type after aging (dirty old snow), near infrared albedo (09/12/2022)
!$OMP THREADPRIVATE(snowa_aged_nobio_nir)
  REAL(r_std), SAVE :: snowa_dec_nobio_vis              !! Decay rate of snow albedo value for nobio type, visible albedo (09/12/2022)
!$OMP THREADPRIVATE(snowa_dec_nobio_vis) 
  REAL(r_std), SAVE :: snowa_dec_nobio_nir              !! Decay rate of snow albedo value for nobio type, near infrared albedo (09/12/2022)
!$OMP THREADPRIVATE(snowa_dec_nobio_nir)
  REAL(r_std), SAVE :: sneige                           !! Lower limit of snow amount (kg.m^{-2})
!$OMP THREADPRIVATE(sneige)
  REAL(r_std), SAVE :: maxmass_snow                     !! The maximum mass of snow (kg.m^{-2})
!$OMP THREADPRIVATE(maxmass_snow)

  !! Heat capacity
  REAL(r_std), PARAMETER :: rho_water = 1000.           !! Density of water (kg/m3)
  REAL(r_std), PARAMETER :: rho_ice = 920.              !! Density of ice (kg/m3)
  REAL(r_std), PARAMETER :: rho_soil = 2700.            !! Density of soil particles (kg/m3), value from Peters-Lidard et al. 1998
  REAL(r_std), PARAMETER :: rho_orga = 1300.            !! Particle density of organic matter (kg/m3), value from Farouki 1981
  REAL(r_std), PARAMETER :: alpha_VB = 0.58             !! Van Bemmelen parameter (conversion factor from SOM to SOC), Heaton 2016

  !! Thermal conductivities
  REAL(r_std), PARAMETER :: cond_water = 0.6            !! Thermal conductivity of liquid water (W/m/K)
  REAL(r_std), PARAMETER :: cond_ice = 2.2              !! Thermal conductivity of ice (W/m/K)
  REAL(r_std), PARAMETER :: cond_solid = 2.32           !! Thermal conductivity of mineral soil particles (W/m/K)

  !! Time constant of long-term soil humidity (s) 
  REAL(r_std), PARAMETER :: lhf = 0.3336*1.E6           !! Latent heat of fusion (J/kg)

  INTEGER(i_std), SAVE      :: nsnow                    !! Number of levels in the snow for explicit snow scheme
!$OMP THREADPRIVATE(nsnow)
  INTEGER(i_std), SAVE      :: nice                     !! Number of levels in ice (for ice-sheet area)
!$OMP THREADPRIVATE(nice)
  REAL(r_std), SAVE         :: frac_min_ice             !! Limit of fraction of ice in the grid cell where the ice mask is considered as 1
!$OMP THREADPRIVATE(frac_min_ice)
  LOGICAL, SAVE             :: read_ice_sheet_mask      !! If true, read the ice sheet mask from file
!$OMP THREADPRIVATE(read_ice_sheet_mask)

  REAL(r_std), PARAMETER    :: XMD    = 28.9644E-3 
  REAL(r_std), PARAMETER    :: XBOLTZ      = 1.380658E-23 
  REAL(r_std), PARAMETER    :: XAVOGADRO   = 6.0221367E+23 
  REAL(r_std), PARAMETER    :: XRD    = XAVOGADRO * XBOLTZ / XMD 
  REAL(r_std), PARAMETER    :: XCPD   = 7.* XRD /2. 
  REAL(r_std), PARAMETER    :: phigeoth = 0.057 ! 0. DKtest 
  REAL(r_std), PARAMETER    :: thick_min_snow = .01 

  !! The maximum snow density and water holding characterisicts 
  REAL(r_std), SAVE         :: xrhosmax                            !! (kg m-3) 
!$OMP THREADPRIVATE(xrhosmax)
  REAL(r_std), SAVE         :: xwsnowholdmax1                      !! (-) 
!$OMP THREADPRIVATE(xwsnowholdmax1)
  REAL(r_std), SAVE         :: xwsnowholdmax2                      !! (-) 
!$OMP THREADPRIVATE(xwsnowholdmax2)
  REAL(r_std), SAVE         :: xsnowrhohold                        !! (kg/m3) 
!$OMP THREADPRIVATE(xsnowrhohold)
  REAL(r_std), SAVE         :: xrhosmin                            !! (kg m-3) 
!$OMP THREADPRIVATE(xrhosmin)
  REAL(r_std), PARAMETER    :: xci = 2.106e+3 
  REAL(r_std), PARAMETER    :: xrv = 6.0221367e+23 * 1.380658e-23 /18.0153e-3 !! ISBA-ES Critical snow depth at which snow grid thicknesses constant 
  REAL(r_std), PARAMETER    :: xsnowcritd = 0.03                   !! (m) 

  !! The threshold of snow depth used for preventing numerical problem in thermal calculations
  REAL(r_std), PARAMETER    :: snowcritd_thermal = 0.01            !! (m)  
  
  !! ISBA-ES CROCUS (Pahaut 1976): snowfall density coefficients: 
  REAL(r_std), PARAMETER       :: snowfall_a_sn = 109.0            !! (kg/m3) 
  REAL(r_std), PARAMETER       :: snowfall_b_sn =   6.0            !! (kg/m3/K) 
  REAL(r_std), PARAMETER       :: snowfall_c_sn =  26.0            !! [kg/(m7/2 s1/2)] 

  REAL(r_std), PARAMETER       :: dgrain_new_max=  2.0e-4          !! (m) : Maximum grain size of new snowfall 
  
  !! Used in explicitsnow to prevent numerical problems as snow becomes vanishingly thin. 
  REAL(r_std), PARAMETER                :: psnowdzmin = .0001      !! m 
  REAL(r_std), PARAMETER                :: xsnowdmin = .000001     !! m 

  REAL(r_std), PARAMETER                :: ph2o = 1000.            !! Water density [kg/m3] 
  
  ! ISBA-ES Thermal conductivity coefficients from Anderson (1976): 
  ! see Boone, Meteo-France/CNRM Note de Centre No. 70 (2002) 
  REAL(r_std), SAVE                     :: ZSNOWTHRMCOND1          !! [W/m/K] 
!$OMP THREADPRIVATE(ZSNOWTHRMCOND1)
  REAL(r_std), SAVE                     :: ZSNOWTHRMCOND2          !! [W m5/(kg2 K)] 
!$OMP THREADPRIVATE(ZSNOWTHRMCOND2)
  
  ! Thermal ice conductivity coefficients Yin CHao Yen 1981 (Temp > 195K): 
  REAL(r_std), SAVE                     :: ZICETHRMCOND1           !! [W/m/K] 
  REAL(r_std), SAVE                     :: ZICETHRMCOND2           !! [W/m/K]
  ! Thermal ice heat capacity coefficients : 
  REAL(r_std), SAVE                     :: ZICETHRMHEAT1           !! [W/m/K]
  REAL(r_std), SAVE                     :: ZICETHRMHEAT2           !! [W/m/K]
  
  ! ISBA-ES Thermal conductivity: Implicit vapor diffn effects 
  ! (sig only for new snow OR high altitudes) 
  ! from Sun et al. (1999): based on data from Jordan (1991) 
  ! see Boone, Meteo-France/CNRM Note de Centre No. 70 (2002) 
  ! 
  REAL(r_std), SAVE                       :: ZSNOWTHRMCOND_AVAP    !! (W/m/K) 
!$OMP THREADPRIVATE(ZSNOWTHRMCOND_AVAP)
  REAL(r_std), SAVE                       :: ZSNOWTHRMCOND_BVAP    !! (W/m) 
!$OMP THREADPRIVATE(ZSNOWTHRMCOND_BVAP)
  REAL(r_std), SAVE                       :: ZSNOWTHRMCOND_CVAP    !! (K) 
!$OMP THREADPRIVATE(ZSNOWTHRMCOND_CVAP)
  REAL(r_std),SAVE :: xansmax                                      !! Maxmimum snow albedo
!$OMP THREADPRIVATE(xansmax)
  REAL(r_std),SAVE :: xansmin                                      !! Miniumum snow albedo
!$OMP THREADPRIVATE(xansmin)
  REAL(r_std),SAVE :: xans_todry                                   !! Albedo decay rate for dry snow
!$OMP THREADPRIVATE(xans_todry)
  REAL(r_std),SAVE :: xans_t                                       !! Albedo decay rate for wet snow
!$OMP THREADPRIVATE(xans_t)

  ! ISBA-ES Thermal conductivity coefficients from Anderson (1976):
  ! see Boone, Meteo-France/CNRM Note de Centre No. 70 (2002)
  REAL(r_std), PARAMETER                  :: XP00 = 1.E5

  ! ISBA-ES Thermal conductivity: Implicit vapor diffn effects
  ! (sig only for new snow OR high altitudes)
  ! from Sun et al. (1999): based on data from Jordan (1991)
  ! see Boone, Meteo-France/CNRM Note de Centre No. 70 (2002)
  !
  REAL(r_std), SAVE          :: ZSNOWCMPCT_RHOD                   !! (kg/m3)
!$OMP THREADPRIVATE(ZSNOWCMPCT_RHOD)
  REAL(r_std), SAVE          :: ZSNOWCMPCT_ACM                    !! (1/s)
!$OMP THREADPRIVATE(ZSNOWCMPCT_ACM)
  REAL(r_std), SAVE          :: ZSNOWCMPCT_BCM                    !! (1/K) 
!$OMP THREADPRIVATE(ZSNOWCMPCT_BCM)
  REAL(r_std), SAVE          :: ZSNOWCMPCT_CCM                    !! (m3/kg)
!$OMP THREADPRIVATE(ZSNOWCMPCT_CCM)
  REAL(r_std), SAVE          :: ZSNOWCMPCT_V0                     !! (Pa/s) 
!$OMP THREADPRIVATE(ZSNOWCMPCT_V0)
  REAL(r_std), SAVE          :: ZSNOWCMPCT_VT                     !! (1/K)
!$OMP THREADPRIVATE(ZSNOWCMPCT_VT)
  REAL(r_std), SAVE          :: ZSNOWCMPCT_VR                     !! (m3/kg)
!$OMP THREADPRIVATE(ZSNOWCMPCT_VR) 
  REAL(r_std), SAVE          :: eta_0                             !! (Pa/s)
!$OMP THREADPRIVATE(eta_0)
  REAL(r_std), SAVE          :: c_eta
!$OMP THREADPRIVATE(c_eta)
  REAL(r_std), SAVE          :: a_eta 
!$OMP THREADPRIVATE(a_eta)
  REAL(r_std), SAVE          :: b_eta
!$OMP THREADPRIVATE(b_eta)
  REAL(r_std), SAVE          :: gs_vionnet
!$OMP THREADPRIVATE(gs_vionnet)
  LOGICAL, SAVE              :: RHOS_WS
!$OMP THREADPRIVATE(RHOS_WS)
  LOGICAL, SAVE              :: INIT_SNOWPACK
!$OMP THREADPRIVATE(INIT_SNOWPACK)
  LOGICAL, SAVE              :: SURFDENSITY_XRHOSMIN
!$OMP THREADPRIVATE(SURFDENSITY_XRHOSMIN)
  REAL(r_std), DIMENSION(12), SAVE :: Zmax                        !! (m)
!$OMP THREADPRIVATE(Zmax)
  INTEGER(i_std), SAVE       :: COMPACT_SNOW_METHOD               !! Method to calculate snow compaction
!$OMP THREADPRIVATE(COMPACT_SNOW_METHOD)
  INTEGER(i_std), SAVE       :: VISCOSITY_METHOD                  !! Method to calculate snow viscosity
!$OMP THREADPRIVATE(VISCOSITY_METHOD)
  
  ! Constants for use of the compaction scheme from Decharme et al.(2016)
  REAL(r_std), SAVE          :: xrhosref                          !!  (kg/m3)
  REAL(r_std), SAVE          :: xrhosmob                          !!  (kg/m3)
  REAL(r_std), SAVE          :: xrhoswind                         !!  (kg/m3)

  REAL(r_std), SAVE          :: ZSNOWCMPCT_X0                     !! (Pa/s)

  REAL(r_std), SAVE          :: ZSNOWCMPCT_XT                     !! (1/K)
  REAL(r_std), SAVE          :: ZSNOWCMPCT_XR                     !! (m3/kg)
  REAL(r_std), SAVE          :: DTref                             !! (K)

  REAL(r_std), SAVE          :: ACMPCT           
  REAL(r_std), SAVE          :: BCMPCT                            !! (s/m)
  REAL(r_std), SAVE          :: KCMPCT                            !! (unitless)

  REAL(r_std), SAVE          :: a_tau                             !! 
  REAL(r_std), SAVE          :: b_tau                             !!

  REAL(r_std), SAVE          :: Amob                              !! (unitless) 

  !
  ! PARAMETERS USED IN TUZET
  !
  INTEGER(i_std), SAVE       :: n_iter_tuzet                      !! (unitless)


  !
  ! PARAMETERS USED FOR CANOPY LAYERS (Albedo, photosynthesis, energy budget)
  !

  INTEGER(i_std), PARAMETER   :: nlevels = 1                  !! Originally the number of levels in the canopy used in  
                                                              !! calculation of the energy budget. After the
                                                              !! mleb calculations have been implemented in enerbil, 
                                                              !! the jnlvls levels are determining the levels used 
                                                              !! within the multi-layer energy budget calculations. However, 
                                                              !! nlevels are still used in calculate_z_level_photo. Cannot 
                                                              !! be deleted before any decision regarding the vertical 
                                                              !! layering has been made.
  INTEGER(i_std), SAVE        :: nlai                         !! Number of levels in the canopy used in the photosynthesis
!$OMP THREADPRIVATE(nlai)
  INTEGER(i_std), SAVE        :: nlevels_tot                  !! Total number of levels, nlevels*nlai
!$OMP THREADPRIVATE(nlevels_tot)
  INTEGER(i_std), SAVE        :: jnlvls                       !! Number of levels in the multilayer energy budget scheme
!$OMP THREADPRIVATE(jnlvls)
  INTEGER(i_std), SAVE        :: jnlvls_under                 !! Number of levels in the understorey of the multilayer energy budget scheme
!$OMP THREADPRIVATE(jnlvls_under)
  INTEGER(i_std), SAVE        :: jnlvls_canopy                !! Number of levels in the canopy of the multilayer energy budget scheme
!$OMP THREADPRIVATE(jnlvls_canopy)
  INTEGER(i_std), SAVE        :: jnlvls_over                  !! Number of levels in the overstorey of the multilayer energy budget scheme
!$OMP THREADPRIVATE(jnlvls_over)


  INTEGER(i_std), SAVE        :: nlev_top                     !! Maximum number of canopy levels at the "top"
!$OMP THREADPRIVATE(nlev_top) 
  REAL(r_std), PARAMETER, &
         DIMENSION (nlevels) :: z_level = (/ 0.0 /)           !! The height of the bottom of each canopy layer                                              
                                                              !! @tex $(m)$ @endtex


  !
  ! Parameters for determining the effective LAI for use in Pinty's albedo scheme
  !
  REAL(r_std), SAVE          ::  laieff_solar_angle           !! the zenith angle for effective LAI
!$OMP THREADPRIVATE(laieff_solar_angle)
  REAL(r_std), SAVE          ::  laieff_zero_cutoff           !! Cutoff for effective lai values
!$OMP THREADPRIVATE(laieff_zero_cutoff)
  REAL(r_std), SAVE          ::  direct_light_weight          !! The weighting factor to weight different sources of light
!$OMP THREADPRIVATE(direct_light_weight)

  !
  ! PARAMETERS FOR HYDRAULIC ARCHITECTURE
  !
  REAL(r_std), PARAMETER, DIMENSION(2)    :: a_viscosity = (/0.556,0.022/) !! Empirical parameters to adjust the resistance of fine
                                                                           !! root and sapwood to the temperature dependency of the
                                                                           !! viscosity of water Cochard et al 2000
  !
  ! BVOC : Biogenic activity  for each age class
  !
  REAL(r_std), SAVE, DIMENSION(nleafages) :: iso_activity                  !! Biogenic activity for each 
                                                                           !! age class : isoprene (unitless)
!$OMP THREADPRIVATE(iso_activity)
  REAL(r_std), SAVE, DIMENSION(nleafages) :: methanol_activity             !! Biogenic activity for each
                                                                           !! age class : methanol (unnitless)
!$OMP THREADPRIVATE(methanol_activity)


 !
 ! Parameters for irrigation scheme
 !
  REAL(r_std), SAVE :: irrig_dosmax                   !! The maximum irrigation water injected per hour (kg.m^{-2}/hour)
!$OMP THREADPRIVATE(irrig_dosmax)
  REAL(r_std), SAVE :: cum_dh_thr                     !! Cumulated nroot threshoold to define root zone, and calculate water deficit for irrigation (-)
!$OMP THREADPRIVATE(cum_dh_thr)
  LOGICAL, SAVE :: irrigated_soiltile                 !! Do we introduce a new soil tile for irrigated croplands? (true/false)
!$OMP THREADPRIVATE(irrigated_soiltile)
  LOGICAL, SAVE :: old_irrig_scheme                   !! Do we run with the old irrigation scheme? (true/false)  , add to compatiblity
!$OMP THREADPRIVATE(old_irrig_scheme)
  INTEGER, SAVE :: irrig_st                           !! Which is the soil tile with irrigation flux
!$OMP THREADPRIVATE(irrig_st)
  REAL(r_std), SAVE, DIMENSION(3) :: avail_reserve    !! Available water from routing reservoirs, to withdraw for irrigation
                                                      !! IMPORTANT: As the routing model uses 3 reservoirs, dimension is set to 3
                                                      !! IMPORTANT: Order of available water must be in this order: streamflow, fast, and slow reservoir
!$OMP THREADPRIVATE(avail_reserve)
  REAL(r_std), SAVE :: beta_irrig                     !! Threshold multiplier of Target SM to calculate root deficit(unitless)
!$OMP THREADPRIVATE(beta_irrig)
  REAL(r_std), SAVE :: lai_irrig_min                  !! Minimum LAI to trigger irrigation (kg.m^{-2}/hour)
!$OMP THREADPRIVATE(lai_irrig_min)
  LOGICAL, SAVE :: irrig_map_dynamic_flag             !! Do we use a dynamic irrig map?
!$OMP THREADPRIVATE(irrig_map_dynamic_flag)
  LOGICAL, SAVE :: select_source_irrig                !! Do we use the new priorization scheme, based on maps of equipped area with surface water?
!$OMP THREADPRIVATE(select_source_irrig)
  LOGICAL, SAVE :: Reinfiltr_IrrigField               !! Do we reinfiltrate all runoff from crop soil tile?O
!$OMP THREADPRIVATE(Reinfiltr_IrrigField)
  REAL, SAVE :: reinf_slope_cropParam                 !! Externalized for irrigated cropland, when Reinfiltr_IrrigField=.TRUE.
                                                      !! Max value of reinf_slope in irrig_st  
!$OMP THREADPRIVATE(reinf_slope_cropParam)
  REAL, SAVE :: a_stream_adduction                    !! Externalized for available volume to adduction
!$OMP THREADPRIVATE(a_stream_adduction)



  !
  ! condveg.f90
  !

  ! 1. Scalar

  ! 1.1 Flags used inside the module

  LOGICAL, SAVE :: alb_bare_model           !! Switch for choosing values of bare soil 
                                            !! albedo (see header of subroutine)
                                            !! (true/false)
!$OMP THREADPRIVATE(alb_bare_model)
  LOGICAL, SAVE :: alb_bg_modis             !! Switch for choosing values of bare soil background albedo read from file (true/false)
!$OMP THREADPRIVATE(alb_bg_modis)
  LOGICAL, SAVE :: alb_bg_modis_reinit      !! Impose reading from file even if variable is found in the restart file(true/false)
!$OMP THREADPRIVATE(alb_bg_modis_reinit)
  LOGICAL, SAVE :: impaze                   !! Switch for choosing surface parameters
                                            !! (see header of subroutine).  
                                            !! (true/false)
!$OMP THREADPRIVATE(impaze)
  LOGICAL, SAVE :: rough_dyn                !! Chooses between two methods to calculate the 
                                            !! the roughness height : static or dynamic (varying with LAI)
                                            !! (true/false)
!$OMP THREADPRIVATE(rough_dyn)
  LOGICAL, SAVE :: use_ratio_z0m_z0h        !! To impose a constant ratio as done in ROUGH_DYN=F
!$OMP THREADPRIVATE(use_ratio_z0m_z0h)

  LOGICAL, SAVE :: sla_dyn                  !! Chooses between two methods to calculate the 
                                            !! specific leaf area: static or dynamic (varying with LAI or biomass)
                                            !! (true/false)
!$OMP THREADPRIVATE(sla_dyn)

  LOGICAL, SAVE :: new_watstress
!$OMP THREADPRIVATE(new_watstress)

  REAL(r_std), SAVE :: alpha_watstress
!$OMP THREADPRIVATE(alpha_watstress)
  LOGICAL, SAVE :: use_height_dom           !! Use the dominant vegetation height instead of the average height when calculating roughness length
!$OMP THREADPRIVATE(use_height_dom)

  ! 1.2 Others 

  REAL(r_std), SAVE :: height_displacement               !! Factor to calculate the zero-plane displacement
                                                         !! height from vegetation height (m)
!$OMP THREADPRIVATE(height_displacement)
  REAL(r_std), SAVE :: crown_packing                     !! Parameter to describe how much of the canopy space could 
                                                         !! be filled with crowns (think: close packing of equal spheres)
!$OMP THREADPRIVATE(crown_packing)
  REAL(r_std), SAVE :: alpha_height_precip               !! alpha parameter in height-precipitation relationship
!$OMP THREADPRIVATE(alpha_height_precip)
  REAL(r_std), SAVE :: beta_height_precip                !! beta parameter in height-precipitation relationship
!$OMP THREADPRIVATE(beta_height_precip)
  REAL(r_std), SAVE :: z0_bare                           !! bare soil roughness length (m)
!$OMP THREADPRIVATE(z0_bare)
  REAL(r_std), SAVE :: z0_ice                            !! ice roughness length (m) --> Standard value
!$OMP THREADPRIVATE(z0_ice)
  REAL(r_std), SAVE :: tcst_snowa                        !! Time constant of the albedo decay of snow (days)
!$OMP THREADPRIVATE(tcst_snowa)
  REAL(r_std), SAVE :: tcst_snowa_nobio                  !! Time constant of the albedo decay of snow on nobio (days) (18/02/2018)
!$OMP THREADPRIVATE(tcst_snowa_nobio)
REAL(r_std), SAVE :: omg1                                !! First tuning constant for snow ageing on nobio areas (14/10/2020)
!$OMP THREADPRIVATE(omg1)
REAL(r_std), SAVE :: omg2                                !! Second tuning constant for snow ageing on nobio areas (14/10/2020)
!$OMP THREADPRIVATE(omg2)
  REAL(r_std), SAVE :: snowcri_alb                       !! Critical value for computation of snow albedo (cm) 
!$OMP THREADPRIVATE(snowcri_alb)
  REAL(r_std), SAVE :: fixed_snow_albedo                 !! To choose a fixed snow albedo value (unitless)
!$OMP THREADPRIVATE(fixed_snow_albedo)
  REAL(r_std), SAVE :: z0_scal                           !! Surface roughness height imposed (m)
!$OMP THREADPRIVATE(z0_scal)
  REAL(r_std), SAVE :: roughheight_scal                  !! Effective roughness Height depending on zero-plane 
                                                         !! displacement height (m) (imposed)
!$OMP THREADPRIVATE(roughheight_scal)
  REAL(r_std), SAVE :: emis_scal                         !! Surface emissivity imposed (unitless)
!$OMP THREADPRIVATE(emis_scal)

  REAL(r_std), SAVE :: c1                                !! Constant used in the formulation of the ratio of 
!$OMP THREADPRIVATE(c1)                                  !! friction velocity to the wind speed at the canopy top
                                                         !! see Ershadi et al. (2015) for more info
  REAL(r_std), SAVE :: c2                                !! Constant used in the formulation of the ratio of 
!$OMP THREADPRIVATE(c2)                                  !! friction velocity to the wind speed at the canopy top
                                                         !! see Ershadi et al. (2015) for more info
  REAL(r_std), SAVE :: c3                                !! Constant used in the formulation of the ratio of 
!$OMP THREADPRIVATE(c3)                                  !! friction velocity to the wind speed at the canopy top
                                                         !! see Ershadi et al. (2015) for more info
  REAL(r_std), SAVE :: Cdrag_foliage                     !! Drag coefficient of the foliage
!$OMP THREADPRIVATE(Cdrag_foliage)                       !! See Ershadi et al. (2015) and Su et. al (2001) for more info
  REAL(r_std), SAVE :: Ct                                !! Heat transfer coefficient of the leaf
!$OMP THREADPRIVATE(Ct)                                  !! See Ershadi et al. (2015) and Su et. al (2001) for more info
  REAL(r_std), SAVE :: Prandtl                           !! Prandtl number used in the calculation of Ct_star
!$OMP THREADPRIVATE(Prandtl)                             !! See Su et. al (2001) for more info



  ! 2. Arrays

  
! 2. Arrays
 
  REAL(r_std), SAVE, DIMENSION(n_spectralbands) :: alb_deadleaf  !! albedo of dead leaves, VIS+NIR (unitless)
!$OMP THREADPRIVATE(alb_deadleaf)
  REAL(r_std), SAVE, DIMENSION(n_spectralbands) :: alb_ice       !! albedo of ice, VIS+NIR (unitless)
!$OMP THREADPRIVATE(alb_ice)
  REAL(r_std), SAVE, DIMENSION(n_spectralbands) :: albedo_scal   !! Albedo values for visible and near-infrared 
                                                                 !! used imposed (unitless) 
!$OMP THREADPRIVATE(albedo_scal)
  REAL(r_std) , SAVE, DIMENSION(classnb) :: vis_dry              !! Soil albedo values to soil colour classification:
                                                                 !! dry soil albedo values in visible range
!$OMP THREADPRIVATE(vis_dry)
  REAL(r_std), SAVE, DIMENSION(classnb) :: nir_dry               !! Soil albedo values to soil colour classification:
                                                                 !! dry soil albedo values in near-infrared range 
!$OMP THREADPRIVATE(nir_dry)
  REAL(r_std), SAVE, DIMENSION(classnb) :: vis_wet               !! Soil albedo values to soil colour classification:
                                                                 !! wet soil albedo values in visible range 
!$OMP THREADPRIVATE(vis_wet)
  REAL(r_std), SAVE, DIMENSION(classnb) :: nir_wet               !! Soil albedo values to soil colour classification:
                                                                 !! wet soil albedo values in near-infrared range
!$OMP THREADPRIVATE(nir_wet)
  REAL(r_std), SAVE, DIMENSION(classnb) :: albsoil_vis           !! Soil albedo values to soil colour classification:
                                                                 !! Averaged of wet and dry soil albedo values
                                                                 !! in visible and near-infrared range
!$OMP THREADPRIVATE(albsoil_vis) 
  REAL(r_std), SAVE, DIMENSION(classnb) :: albsoil_nir           !! Soil albedo values to soil colour classification:
                                                                 !! Averaged of wet and dry soil albedo values
                                                                 !! in visible and near-infrared range
!$OMP THREADPRIVATE(albsoil_nir)
  REAL(r_std), PARAMETER                :: alb_threshold = 0.0000000001_r_std !! A threshold for the iteration of the
                                                                 !! multilevel albedo.  Could be externalised.
                                                                 !! Fairly arbitrary, although if a level has
                                                                 !! no LAI the absorption often ends up being
                                                                 !! equal to this value, so it should not
                                                                 !! be high.
  REAL(r_std), SAVE, DIMENSION(8) :: ZSICOEF1                    !! Ice grid parameters 8 layers : thickness of ice levels
!$OMP THREADPRIVATE(zsicoef1)

  !
  ! diffuco.f90
  !
  
  ! 0. Constants

  REAL(r_std), PARAMETER :: Tetens_1 = 0.622         !! Ratio between molecular weight of water vapor and molecular weight  
                                                     !! of dry air (unitless)
  REAL(r_std), PARAMETER :: Tetens_2 = 0.378         !!
  REAL(r_std), PARAMETER :: ratio_H2O_to_CO2 = 1.6   !! Ratio of water vapor diffusivity to the CO2 diffusivity (unitless)
  REAL(r_std), PARAMETER :: mol_to_m_1 = 0.0244      !!
  REAL(r_std), PARAMETER :: RG_to_PAR = 0.5          !!
  REAL(r_std), PARAMETER :: W_to_mol = 4.6           !! W_to_mmol * RG_to_PAR = 2.3

  ! 1. Scalar

  LOGICAL, SAVE :: ldq_cdrag_from_gcm                !! Set to .TRUE. if you want q_cdrag coming from GCM
!$OMP THREADPRIVATE(ldq_cdrag_from_gcm)
  REAL(r_std), SAVE :: laimax                        !! Maximal LAI used for splitting LAI into N layers (m^2.m^{-2})
!$OMP THREADPRIVATE(laimax)
  LOGICAL, SAVE :: downregulation_co2                !! Set to .TRUE. if you want CO2 downregulation.
!$OMP THREADPRIVATE(downregulation_co2)
  REAL(r_std), SAVE :: downregulation_co2_baselevel  !! CO2 base level (ppm)
!$OMP THREADPRIVATE(downregulation_co2_baselevel)
  REAL(r_std), SAVE :: gb_ref                        !! Leaf bulk boundary layer resistance (s m-1)
!$OMP THREADPRIVATE(gb_ref)

  ! 3. Coefficients of equations

  REAL(r_std), SAVE :: lai_level_depth               !!
!$OMP THREADPRIVATE(lai_level_depth)
  REAL(r_std), SAVE               :: c13_a           !! fractionation against during diffusion
!$OMP THREADPRIVATE(c13_a)
  REAL(r_std), SAVE               :: c13_b           !! fractionation against during carboxylation
!$OMP THREADPRIVATE(c13_b)
  REAL(r_std), SAVE               :: threshold_c13_assim !! If assimilation falls below this threshold
                                                     !! the delta_c13 is set to zero
!$OMP THREADPRIVATE(threshold_c13_assim)
  !
  REAL(r_std), SAVE, DIMENSION(6) :: dew_veg_poly_coeff !! coefficients of the 5 degree polynomomial used
                                                     !! in the equation of coeff_dew_veg
!$OMP THREADPRIVATE(dew_veg_poly_coeff)
!
  REAL(r_std), SAVE               :: Oi              !! Intercellular oxygen partial pressure (ubar)
!$OMP THREADPRIVATE(Oi)
  !
  ! slowproc.f90 
  !

  ! 1. Scalar

  INTEGER(i_std), SAVE :: ninput_year_orig          !!  first year for N inputs (number)
!$OMP THREADPRIVATE(ninput_year_orig)
  LOGICAL, SAVE :: ninput_suffix_year               !! Do the Ninput datasets have a 'year' suffix ? (y/n)  
!$OMP THREADPRIVATE(ninput_suffix_year)
  REAL(r_std), SAVE :: bulk_default                 !! Default value for bulk density of soil (kg/m3)
!$OMP THREADPRIVATE(bulk_default)
  REAL(r_std), SAVE :: ph_default                   !! Default value for pH of soil (-)
!$OMP THREADPRIVATE(ph_default)

  REAL(r_std), SAVE :: min_vegfrac                  !! Minimal fraction of mesh a vegetation type can occupy (0-1, unitless)
!$OMP THREADPRIVATE(min_vegfrac)
  REAL(r_std), SAVE :: frac_nobio_fixed_test_1      !! Value for frac_nobio for tests in 0-dim simulations (0-1, unitless)
!$OMP THREADPRIVATE(frac_nobio_fixed_test_1)
  
  REAL(r_std), SAVE :: stempdiag_bid                !! only needed for an initial LAI if there is no restart file
!$OMP THREADPRIVATE(stempdiag_bid)




  !
  ! lake.f90 
  !
  ! Flag used in the lake module
  LOGICAL, SAVE                        :: ok_depth_lake_max = .TRUE.       !! Activation of Mean lake depth Truncation
  !$OMP THREADPRIVATE(ok_depth_lake_max) 

  ! Scalar to use in lake module
  INTEGER(i_std), SAVE                :: nlake = 3                !! Number of lake types  
!$OMP THREADPRIVATE(nlake)
  REAL(r_std), SAVE                   :: T_clim_lake = 277.15     !! Temperature to initialize the Temperature column 
!$OMP THREADPRIVATE(T_clim_lake) 
  REAL(r_std), SAVE                   :: depth_lake_max = 60.     !! Truncation of lake depth 
!$OMP THREADPRIVATE(depth_lake_max)
  REAL(r_std), SAVE                   :: albedo_water_lake = 0.07 !! Water surface albedo with respect to the solar radiation
!$OMP THREADPRIVATE(albedo_water_lake)
  REAL(r_std), SAVE                   :: albedo_ice_lake = 0.4    !! Ice surface albedo on the lake with respect to the solar radiation
!$OMP THREADPRIVATE(albedo_ice_lake)
  REAL(r_std), SAVE                   :: albedo_snow_lake = 0.6   !! Snow surface albedo on the lake with respect to the solar radiation
!$OMP THREADPRIVATE(albedo_snow_lake)

!  Optical characteristics for water, ice and snow.
!  The simplest one-band approximation is used as a reference. 

  REAL(r_std), SAVE                   :: light_ext_ice_lake = 8.      
!$OMP THREADPRIVATE(light_ext_ice_lake)
  REAL(r_std), SAVE                   :: light_ext_snow_lake = 15.
!$OMP THREADPRIVATE(light_ext_snow_lake)

  REAL(r_std), SAVE                   :: t_sed_lake = 277.15      !! Temperature at the outer edge of 
                                                                  !! the thermally active layer of the bottom sediments [K]
!$OMP THREADPRIVATE(t_sed_lake)
  REAL(r_std), SAVE                   :: depth_sed_lake = 10      !! Depth of the thermally active layer of the bottom sediments [m]
!$OMP THREADPRIVATE(depth_sed_lake)

  REAL(r_std), SAVE                   :: lake_depth_cons = 0      !! Depth of the lake if there is not lake map to initialize [m]
!$OMP THREADPRIVATE(lake_depth_cons)

  REAL(r_std), PARAMETER :: h_Snow_min_flk = 1.0E-5        !! Minimum snow thickness [m]
  REAL(r_std), PARAMETER :: h_Ice_min_flk  = 1.0E-9        !! Minimum ice thickness [m]
  REAL(r_std), PARAMETER :: h_ML_min_flk   = 1.0E-2        !! Minimum mixed-layer depth [m]
  REAL(r_std), PARAMETER :: h_ML_max_flk   = 1.0E+3        !! Maximum mixed-layer depth [m]   
  REAL(r_std), PARAMETER :: H_B1_min_flk   = 1.0E-3        !! Minimum thickness of the upper layer of bottom sediments [m]
  REAL(r_std), PARAMETER :: u_star_min_flk = 1.0E-6        !! Minimum value of the surface friction velocity [m s^{-1}]

  REAL(r_std), PARAMETER :: b1_vap   = 610.78              !! Coefficient [N m^{-2} = kg m^{-1} s^{-2}]
  REAL(r_std), PARAMETER :: b3_vap   = 273.16              !! Triple point [K]
  REAL(r_std), PARAMETER :: b2w_vap  = 17.2693882          !! Coefficient (water)
  REAL(r_std), PARAMETER :: b2i_vap  = 21.8745584          !! Coefficient (ice) 
  REAL(r_std), PARAMETER :: b4w_vap  = 35.86               !! Coefficient (temperature) [K]
  REAL(r_std), PARAMETER :: b4i_vap  = 7.66                !! Coefficient (temperature) [K]		

  !  Derived thermodynamic parameters 
  REAL(r_std), PARAMETER :: num_1o3_sf = 1./3.             !! 1/3
  REAL(r_std), PARAMETER :: tpl_grav          = 9.81       !! Acceleration due to gravity [m s^{-2}]
  REAL(r_std), PARAMETER :: tpl_T_r           = 277.13     !! Temperature of maximum density of fresh water [K]
  REAL(r_std), PARAMETER :: tpl_T_f           = 273.15     !! Fresh water freezing point [K]
  REAL(r_std), PARAMETER :: tpl_a_T           = 1.6509E-05 !! Constant in the fresh-water equation of state [K^{-2}]
  REAL(r_std), PARAMETER :: tpl_rho_w_r       = 1.0E+03    !! Maximum density of fresh water [kg m^{-3}]
  REAL(r_std), PARAMETER :: tpl_rho_I         = 9.1E+02    !! Density of ice [kg m^{-3}]
  REAL(r_std), PARAMETER :: tpl_rho_S_min     = 1.0E+02    !! Minimum snow density [kg m^{-3}]
  REAL(r_std), PARAMETER :: tpl_rho_S_max     = 4.0E+02    !! Maximum snow density [kg m^{-3}]
  REAL(r_std), PARAMETER :: tpl_Gamma_rho_S   = 2.0E+02    !! Empirical parameter [kg m^{-4}] (value from Mironov 2012 before 2.0E+02) 
                                                           !! in the expression for the snow density 
  REAL(r_std), PARAMETER :: tpl_L_f           = 3.3E+05    !! Latent heat of fusion [J kg^{-1}]
  REAL(r_std), PARAMETER :: tpl_c_w           = 4.2E+03    !! Specific heat of water [J kg^{-1} K^{-1}]
  REAL(r_std), PARAMETER :: tpl_c_I           = 2.1E+03    !! Specific heat of ice [J kg^{-1} K^{-1}]
  REAL(r_std), PARAMETER :: tpl_c_S           = 2.1E+03    !! Specific heat of snow [J kg^{-1} K^{-1}]
  REAL(r_std), PARAMETER :: tpl_kappa_w       = 5.46E-01   !! Molecular heat conductivity of water [J m^{-1} s^{-1} K^{-1}]
  REAL(r_std), PARAMETER :: tpl_kappa_I       = 2.29       !! Molecular heat conductivity of ice [J m^{-1} s^{-1} K^{-1}]
  REAL(r_std), PARAMETER :: tpl_kappa_S_min   = 0.2        !! Minimum molecular heat conductivity of snow [J m^{-1} s^{-1} K^{-1}] (value from Mironov 2012 before 0.2)
  REAL(r_std), PARAMETER :: tpl_kappa_S_max   = 1.5        !! Maximum molecular heat conductivity of snow [J m^{-1} s^{-1} K^{-1}]
  REAL(r_std), PARAMETER :: tpl_Gamma_kappa_S = 1.3        !! Empirical parameter [J m^{-2} s^{-1} K^{-1}] 
                                                           !! in the expression for the snow heat conductivity  
  REAL(r_std), PARAMETER :: C_T_min       = 0.65           !! Minimum value of the shape factor C_T (thermocline)
  REAL(r_std), PARAMETER :: C_T_max       = 0.8            !! Maximum value of the shape factor C_T (thermocline)
  REAL(r_std), PARAMETER :: Phi_T_pr0_1   = 40./3.         !! Constant in the expression for the T shape-function derivative 
  REAL(r_std), PARAMETER :: Phi_T_pr0_2   = 20./3.         !! Constant in the expression for the T shape-function derivative 
  REAL(r_std), PARAMETER :: C_TT_1        = 11./18.        !! Constant in the expression for C_TT (thermocline)
  REAL(r_std), PARAMETER :: C_TT_2        = 7./45.         !! Constant in the expression for C_TT (thermocline)
  REAL(r_std), PARAMETER :: C_B1          = 2./3.          !! Shape factor (upper layer of bottom sediments)
  REAL(r_std), PARAMETER :: C_B2          = 3./5.          !! Shape factor (lower layer of bottom sediments)
  REAL(r_std), PARAMETER :: Phi_B1_pr0    = 2.             !! B1 shape-function derivative 
  REAL(r_std), PARAMETER :: C_S_lin       = 0.5            !! Shape factor (linear temperature profile in the snow layer)
  REAL(r_std), PARAMETER :: Phi_S_pr0_lin = 1.             !! S shape-function derivative (linear profile) 
  REAL(r_std), PARAMETER :: C_I_lin       = 0.5            !! Shape factor (linear temperature profile in the ice layer)
  REAL(r_std), PARAMETER :: Phi_I_pr0_lin = 1.             !! I shape-function derivative (linear profile) 
  REAL(r_std), PARAMETER :: Phi_I_pr1_lin = 1.             !! I shape-function derivative (linear profile) 
  REAL(r_std), PARAMETER :: Phi_I_ast_MR  = 2.             !! Constant in the MR2004 expression for I shape factor
  REAL(r_std), PARAMETER :: C_I_MR        = 1./12.         !! Constant in the MR2004 expression for I shape factor
  REAL(r_std), PARAMETER :: H_Ice_max     = 3.             !! Maximum ice tickness in 
                                                           !! the Mironov and Ritter (2004, MR2004) ice model [m] 
  REAL(r_std), PARAMETER :: u_wind_min_sf = 1.0E-02        !! Minimum wind speed [m s^{-1}]
  REAL(r_std), PARAMETER :: u_star_min_sf = 1.0E-04        !! Minimum value of friction velocity [m s^{-1}]
  REAL(r_std), PARAMETER :: c_accur_sf    = 1.0E-07        !! A small number (accuracy)
  REAL(r_std), PARAMETER :: c_small_sf    = 1.0E-04        !! A small number (used to compute fluxes)
  REAL(r_std), PARAMETER :: c_small_flk   = 1.0E-10        !! A small number  
  REAL(r_std), PARAMETER :: c_cbl_1       = 0.17           !! Constant in the CBL entrainment equation
  REAL(r_std), PARAMETER :: c_cbl_2       = 1.             !! Constant in the CBL entrainment equation
  REAL(r_std), PARAMETER :: c_sbl_ZM_n    = 0.5            !! Constant in the ZM1996 equation for the equilibrium SBL depth
  REAL(r_std), PARAMETER :: c_sbl_ZM_s    = 10.            !! Constant in the ZM1996 equation for the equilibrium SBL depth
  REAL(r_std), PARAMETER :: c_sbl_ZM_i    = 20.            !! Constant in the ZM1996 equation for the equilibrium SBL depth
  REAL(r_std), PARAMETER :: c_relax_h     = 0.030          !! Constant in the relaxation equation for the SBL depth
  REAL(r_std), PARAMETER :: c_relax_C     = 0.0030         !! Constant in the relaxation equation for the shape factor


  REAL(r_std), PARAMETER :: c_Karman      = 0.40      !! The von Karman constant 
  REAL(r_std), PARAMETER :: Pr_neutral    = 1.0       !! Turbulent Prandtl number at neutral static stability
  REAL(r_std), PARAMETER :: Sc_neutral    = 1.0       !! Turbulent Schmidt number at neutral static stability
  REAL(r_std), PARAMETER :: c_MO_u_stab   = 5.0       !! Constant of the MO theory (wind, stable stratification)
  REAL(r_std), PARAMETER :: c_MO_t_stab   = 5.0       !! Constant of the MO theory (temperature, stable stratification)
  REAL(r_std), PARAMETER :: c_MO_q_stab   = 5.0       !! Constant of the MO theory (humidity, stable stratification)
  REAL(r_std), PARAMETER :: c_MO_u_conv   = 15.0      !! Constant of the MO theory (wind, convection)
  REAL(r_std), PARAMETER :: c_MO_t_conv   = 15.0      !! Constant of the MO theory (temperature, convection)
  REAL(r_std), PARAMETER :: c_MO_q_conv   = 15.0      !! Constant of the MO theory (humidity, convection)
  REAL(r_std), PARAMETER :: c_MO_u_exp    = 0.25      !! Constant of the MO theory (wind, exponent)
  REAL(r_std), PARAMETER :: c_MO_t_exp    = 0.5       !! Constant of the MO theory (temperature, exponent)
  REAL(r_std), PARAMETER :: c_MO_q_exp    = 0.5       !! Constant of the MO theory (humidity, exponent)
  REAL(r_std), PARAMETER :: z0u_ice_rough = 1.0E-03   !! Aerodynamic roughness of the ice surface [m] (rough flow)
  REAL(r_std), PARAMETER :: c_z0u_smooth  = 0.1       !! Constant in the expression for z0u (smooth flow) 
  REAL(r_std), PARAMETER :: c_z0u_rough   = 1.23E-02  !! The Charnock constant in the expression for z0u (rough flow)
  REAL(r_std), PARAMETER :: c_z0u_rough_L = 1.00E-01  !! An increased Charnock constant (used as the upper limit)
  REAL(r_std), PARAMETER :: c_z0u_ftch_f  = 0.70      !! Factor in the expression for fetch-dependent Charnock parameter
  REAL(r_std), PARAMETER :: c_z0u_ftch_ex = 0.3333333 !! Exponent in the expression for fetch-dependent Charnock parameter
  REAL(r_std), PARAMETER :: c_z0t_rough_1 = 4.0       !! Constant in the expression for z0t (factor) 
  REAL(r_std), PARAMETER :: c_z0t_rough_2 = 3.2       !! Constant in the expression for z0t (factor)
  REAL(r_std), PARAMETER :: c_z0t_rough_3 = 0.5       !! Constant in the expression for z0t (exponent) 
  REAL(r_std), PARAMETER :: c_z0q_rough_1 = 4.0       !! Constant in the expression for z0q (factor)
  REAL(r_std), PARAMETER :: c_z0q_rough_2 = 4.2       !! Constant in the expression for z0q (factor)
  REAL(r_std), PARAMETER :: c_z0q_rough_3 = 0.5       !! Constant in the expression for z0q (exponent)
  REAL(r_std), PARAMETER :: c_z0t_ice_b0s = 1.250     !! Constant in the expression for z0t over ice
  REAL(r_std), PARAMETER :: c_z0t_ice_b0t = 0.149     !! Constant in the expression for z0t over ice
  REAL(r_std), PARAMETER :: c_z0t_ice_b1t = -0.550    !! Constant in the expression for z0t over ice
  REAL(r_std), PARAMETER :: c_z0t_ice_b0r = 0.317     !! Constant in the expression for z0t over ice
  REAL(r_std), PARAMETER :: c_z0t_ice_b1r = -0.565    !! Constant in the expression for z0t over ice
  REAL(r_std), PARAMETER :: c_z0t_ice_b2r = -0.183    !! Constant in the expression for z0t over ice
  REAL(r_std), PARAMETER :: c_z0q_ice_b0s = 1.610     !! Constant in the expression for z0q over ice
  REAL(r_std), PARAMETER :: c_z0q_ice_b0t = 0.351     !! Constant in the expression for z0q over ice
  REAL(r_std), PARAMETER :: c_z0q_ice_b1t = -0.628    !! Constant in the expression for z0q over ice
  REAL(r_std), PARAMETER :: c_z0q_ice_b0r = 0.396     !! Constant in the expression for z0q over ice
  REAL(r_std), PARAMETER :: c_z0q_ice_b1r = -0.512    !! Constant in the expression for z0q over ice
  REAL(r_std), PARAMETER :: c_z0q_ice_b2r = -0.180    !! Constant in the expression for z0q over ice
  REAL(r_std), PARAMETER :: Re_z0s_ice_t  = 2.5       !! Threshold value of the surface Reynolds number used to compute z0t and z0q over ice (Andreas 2002)
  REAL(r_std), PARAMETER :: Re_z0u_thresh = 0.1       !! Threshold value of the roughness Reynolds number [value from Zilitinkevich, Grachev, and Fairall (200)

  !  Dimensionless constants 
  REAL(r_std), PARAMETER :: c_free_conv   = 0.14      !! Constant in the expressions for fluxes in free convection
  !  Dimensionless constants 
  REAL(r_std), PARAMETER :: c_lwrad_emis  = 1.        !! Surface emissivity with respect to the long-wave radiation

    !  Thermodynamic parameters

  REAL(r_std), PARAMETER :: tpsf_C_StefBoltz = 5.67E-08    !! The Stefan-Boltzmann constant [W m^{-2} K^{-4}]
  REAL(r_std), PARAMETER :: tpsf_R_dryair    = 2.8705E+02  !! Gas constant for dry air [J kg^{-1} K^{-1}]
  REAL(r_std), PARAMETER :: tpsf_R_watvap    = 4.6151E+02  !! Gas constant for water vapour [J kg^{-1} K^{-1}]
  REAL(r_std), PARAMETER :: tpsf_c_a_p       = 1.005E+03   !! Specific heat of air at constant pressure [J kg^{-1} K^{-1}]
  REAL(r_std), PARAMETER :: tpsf_L_evap      = 2.501E+06   !! Specific heat of evaporation [J kg^{-1}]
  REAL(r_std), PARAMETER :: tpsf_nu_u_a      = 1.50E-05    !! Kinematic molecular viscosity of air [m^{2} s^{-1}]
  REAL(r_std), PARAMETER :: tpsf_kappa_t_a   = 2.20E-05    !! Molecular temperature conductivity of air [m^{2} s^{-1}]
  REAL(r_std), PARAMETER :: tpsf_kappa_q_a   = 2.40E-05    !! Molecular diffusivity of air for water vapour [m^{2} s^{-1}]
  
  !  Derived thermodynamic parameters
  REAL(r_std), PARAMETER :: tpsf_Rd_o_Rv  = tpsf_R_dryair/tpsf_R_watvap        !! Ratio of gas constants (Rd/Rv)
  REAL(r_std), PARAMETER :: tpsf_alpha_q  = (1.-tpsf_Rd_o_Rv)/tpsf_Rd_o_Rv     !! Diemsnionless ratio 
  
  !  Thermodynamic parameters
  REAL(r_std), PARAMETER :: P_a_ref = 1.0E+05                   !! Reference pressure [N m^{-2} = kg m^{-1} s^{-2}]

  INTEGER(i_std), PARAMETER :: nband_optic_max = 10
  REAL(r_std), PARAMETER  :: opticpar_water_ref = 3.            !! water ref
  REAL(r_std), PARAMETER  :: opticpar_whiteice_ref = 17.1       !! White ice
  REAL(r_std), PARAMETER  :: opticpar_blueice_ref = 8.4         !! Blue ice
  REAL(r_std), PARAMETER  :: opticpar_drysnow_ref = 25.0        !! Dry snow 
  REAL(r_std), PARAMETER  :: opticpar_meltingsnow_ref = 15.0    !! Melting snow
  REAL(r_std), PARAMETER  :: opticpar_snow_opaque_ref = 1.0E+07 !! Opaque ice
  REAL(r_std), PARAMETER  :: opticpar_ice_opaque_ref = 1.0E+07  !! Opaque snow
  REAL(r_std), PARAMETER  :: albedo_water_lake_ref  = 0.07      !! Water
  REAL(r_std), PARAMETER  :: albedo_whiteice_ref    = 0.50      !! White ice
  REAL(r_std), PARAMETER  :: albedo_blueice_ref     = 0.15      !! Blue ice
  REAL(r_std), PARAMETER  :: albedo_drysnow_ref     = 0.60      !! Dry snow 
  REAL(r_std), PARAMETER  :: albedo_meltingsnow_ref = 0.10      !! Melting snow
  REAL(r_std), PARAMETER  :: albedo_snow_max        = 0.87      !! Dry snow 
  REAL(r_std), PARAMETER  :: albedo_snow_min        = 0.50      !! Melting snow
  REAL(r_std), PARAMETER  :: c_albice_MR = 95.6                 !! Constant in the interpolation formula for 
                                                                !! the ice albedo (Mironov and Ritter 2004)		
  REAL(r_std), SAVE :: height_u_lake  = 2                       !! Height above the lake surface where the wind speed is measured [m]
  REAL(r_std), SAVE :: height_tq_lake = 2                       !! Height where temperature and humidity are measured [m]

   
                           !-----------------------------!
                           !  STOMATE AND LPJ PARAMETERS !
                           !-----------------------------!

  ! Allometric relationships
  REAL(r_std), SAVE :: exp_kf                       !! Tuning parameter for the relationship between tree height and k_latosas
!$OMP THREADPRIVATE(exp_kf)
  !
  ! lpj_constraints.f90
  !
  
  ! 1. Scalar

  REAL(r_std), SAVE  :: too_long                    !! longest sustainable time without 
                                                    !! regeneration (vernalization) (years)
!$OMP THREADPRIVATE(too_long)


  !
  ! lpj_establish.f90
  !

  ! 1. Scalar

  REAL(r_std), SAVE :: estab_max_tree               !! Maximum tree establishment rate (ind/m2/dt_stomate)
!$OMP THREADPRIVATE(estab_max_tree)
  REAL(r_std), SAVE :: estab_max_grass              !! Maximum grass establishment rate (ind/m2/dt_stomate)
!$OMP THREADPRIVATE(estab_max_grass)
  
  ! 3. Coefficients of equations

  REAL(r_std), SAVE :: establish_scal_fact          !!
!$OMP THREADPRIVATE(establish_scal_fact)
  REAL(r_std), SAVE :: max_tree_coverage            !! (0-1, unitless)
!$OMP THREADPRIVATE(max_tree_coverage)
  REAL(r_std), SAVE :: ind_0_estab                  !! = ind_0 * 10.
!$OMP THREADPRIVATE(ind_0_estab)


  !
  ! lpj_fire.f90
  !

  ! 1. Scalar

  REAL(r_std), SAVE :: tau_fire                     !! Time scale for memory of the fire index (days).
!$OMP THREADPRIVATE(tau_fire)
  REAL(r_std), SAVE :: cr_fire                      !! VPD fire-danger precip-suppression coeff (FDI = VPD*exp(-cr_fire*precip_month)); calibration knob, was hardcoded 2.0.
!$OMP THREADPRIVATE(cr_fire)
  REAL(r_std), SAVE :: a_nd_scale                   !! Multiplicative scaling of the human-ignition parameter a_nd (default 1.0); fire knob: higher -> more human ignitions -> more fire. Absorbs scale-dependence of anthropogenic ignition density.
!$OMP THREADPRIVATE(a_nd_scale)
  REAL(r_std), SAVE :: litter_crit                  !! Critical litter quantity for fire
                                                    !! below which iginitions extinguish 
                                                    !! @tex $(gC m^{-2})$ @endtex
!$OMP THREADPRIVATE(litter_crit)
  REAL(r_std), SAVE :: fire_resist_lignin           !!
!$OMP THREADPRIVATE(fire_resist_lignin)
  ! 2. Arrays

  REAL(r_std), SAVE, DIMENSION(nparts) :: co2frac   !! The fraction of the different biomass 
                                                    !! compartments emitted to the atmosphere 
!$OMP THREADPRIVATE(co2frac)                                                         !! when burned (unitless, 0-1)  

  ! 3. Coefficients of equations

  REAL(r_std), SAVE, DIMENSION(3) :: bcfrac_coeff   !! (unitless)
!$OMP THREADPRIVATE(bcfrac_coeff)
  REAL(r_std), SAVE, DIMENSION(4) :: firefrac_coeff !! (unitless)

!$OMP THREADPRIVATE(firefrac_coeff)
  REAL(r_std), SAVE, DIMENSION(nhour) :: alloc_firefuel = (/0.045, 0.075, 0.21, 0.67/)  !! Allocating fraction
                                                          !! of forest aboveground biomass pools, except for
                                                          !leaf and fruit, to different dead fuel
                                                          !categories for fire modelling.
!$OMP THREADPRIVATE(alloc_firefuel)

  !
  ! stomate_spitfire.f90
  !

  REAL(r_std), SAVE          :: cf_threshold               !! Crown fire happens when ignited_fraction is higher than cf_threshold
!$OMP THREADPRIVATE(cf_threshold)
  REAL(r_std), SAVE          :: patch_multiplier_grass     !! Multiplier to get max patch size for grass
!$OMP THREADPRIVATE(patch_multiplier_grass)
  REAL(r_std), SAVE          :: patch_multiplier_forest    !! Multiplier to get max patch size for trees
!$OMP THREADPRIVATE(patch_multiplier_forest)
  REAL(r_std), SAVE          :: edge_ignition_factor       !! The multiplier for increasing ignition due to edge (unitless) (0-1)
!$OMP THREADPRIVATE(edge_ignition_factor)
  REAL(r_std), SAVE          :: edge_wind_factor           !! The multiplier for increasing wind due to edge (unitless) (0-1)
!$OMP THREADPRIVATE(edge_wind_factor)
  REAL(r_std), SAVE          :: wetness_diff               !! Relative fuel wetness difference between edge and interior (unitless) (0-1)
!$OMP THREADPRIVATE(wetness_diff)
  REAL(r_std), SAVE          :: edge_distance_ignition     !! The edge distance for increasing human ignition (unit: m)
!$OMP THREADPRIVATE(edge_distance_ignition)
  REAL(r_std), SAVE          :: edge_distance_wind         !! The edge distance for increasing wind speed (unit: m)
!$OMP THREADPRIVATE(edge_distance_wind)
  REAL(r_std), SAVE          :: edge_distance_wetness      !! The edge distance for decreasing feul moisture (unit: m)
!$OMP THREADPRIVATE(edge_distance_wetness)
  REAL(r_std), SAVE          :: edge_distance_light        !! The edge distance for the canopy light effect (unit: m)
!$OMP THREADPRIVATE(edge_distance_light)
  REAL(r_std), SAVE          :: pgap_edge_factor           !! Relative change of Pgap at the edge vs interior (unitless)
!$OMP THREADPRIVATE(pgap_edge_factor)
  REAL(r_std), SAVE          :: laieff_edge_factor         !! Relative change of LAIeff at the edge vs interior (unitless, negative)
!$OMP THREADPRIVATE(laieff_edge_factor)
  REAL(r_std), SAVE          :: alpha_fire_edge            !! Edge generated per m² burnt (m/m², Marie 2026 AED_FEEDBACK)
!$OMP THREADPRIVATE(alpha_fire_edge)
  REAL(r_std), SAVE          :: alpha_storm_edge           !! Edge generated per m² stormed (m/m²)
!$OMP THREADPRIVATE(alpha_storm_edge)
  REAL(r_std), SAVE          :: alpha_pest_edge            !! Edge generated per m² pest-killed (m/m²)
!$OMP THREADPRIVATE(alpha_pest_edge)
  REAL(r_std), SAVE          :: alpha_harvest_edge         !! Edge generated per m² harvested or LCC-cut (m/m²). Repli scalaire du chemin
                                                           !! legacy (OK_MANAGEMENT_INTENSITY=n et pas d'alpha_map) : CONSERVE pour la
                                                           !! bit-neutralite OFF, au meme titre que alpha_{fire,storm,pest}_edge
!$OMP THREADPRIVATE(alpha_harvest_edge)
  REAL(r_std), SAVE          :: tau_rec_edge               !! Edge-length recovery time toward the feedback baseline (s)
!$OMP THREADPRIVATE(tau_rec_edge)
  REAL(r_std), SAVE          :: aed_edge_ref_frac          !! Fraction of edge_length_ref (EFDA) kept as the continuous-forest
                                                           !! baseline of the edge feedback. 0 = pure perturbation-driven
                                                           !! edge (no EFDA double-count, default) ; 1 = legacy EL_dyn=EFDA+source
!$OMP THREADPRIVATE(aed_edge_ref_frac)
  REAL(r_std), SAVE          :: edge_min_patch_area        !! Min disturbance/forest patch area for the edge cap (m²)
!$OMP THREADPRIVATE(edge_min_patch_area)
  !! --- Intensite de gestion : concept unique dont derivent tous les parametres de
  !!     gestion par maille (design/MODULE_DESIGN_MANAGEMENT_INTENSITY.md, 2026-07-25).
  !!     Remplace le levier fragmentation (ok_harvest_fragmentation, edge_shape_coef,
  !!     clearcut_area_by_fmclass, landscape_scenario, harvest_clearcut_area,
  !!     clearcut_area_file) ET la source spatiale d'AGE_ROTATION (target_rotation_age*).
  INTEGER(i_std), PARAMETER  :: nmiclass = 5               !! Nb de classes d'intensite (Scherpenhuijzen 2025)
  LOGICAL, SAVE              :: ok_management_intensity    !! Derive rotation / alpha_harvest / clearcut_area des fractions de classe
!$OMP THREADPRIVATE(ok_management_intensity)
  CHARACTER(LEN=100), SAVE   :: management_intensity_file  !! NetCDF annuel, fractions f_class1..5 ; "NONE" = desactive
!$OMP THREADPRIVATE(management_intensity_file)
  REAL(r_std), SAVE, DIMENSION(nmiclass) :: mi_clearcut_size !! Taille de parcelle de realisation par classe (m2), 1 ha .. 25 ha
!$OMP THREADPRIVATE(mi_clearcut_size)
  REAL(r_std), SAVE, DIMENSION(nmiclass) :: mi_dia_factor    !! Facteur multiplicatif sur le diametre de coupe par classe d'intensite (-), 1 = inchange
!$OMP THREADPRIVATE(mi_dia_factor)
  INTEGER(i_std), SAVE :: n_class0          !! Guillaume M. -- MATURITY_CONVEYOR: number of class-0
                                            !! slots. Residence time is n_class0 - 1 years.
                                            !! 0 => convoyeur inactif, comportement historique (bit-neutre).
!$OMP THREADPRIVATE(n_class0)
  REAL(r_std), SAVE :: ac12_h_ratio_start   !! MATURITY_CONVEYOR : rapport de hauteur h1/h_ref a partir duquel
                                            !! la classe 1 commence a monter en classe 2 (-)
!$OMP THREADPRIVATE(ac12_h_ratio_start)
  LOGICAL, SAVE :: ok_disturb_via_class0   !! MATURITY_CONVEYOR : feu, tempete et scolyte creditent
                                            !! le registre de classe 0 au lieu de redescendre
                                            !! directement en classe 1. Defaut FALSE => bit-neutre.
!$OMP THREADPRIVATE(ok_disturb_via_class0)
  LOGICAL, SAVE :: ok_class0_replant       !! MATURITY_CONVEYOR : la restitution du convoyeur de classe 0
                                            !! REPLANTE la fraction rendue au creneau de classe 1.
                                            !! Defaut FALSE => bit-neutre.
!$OMP THREADPRIVATE(ok_class0_replant)
  REAL(r_std), SAVE :: ac12_dia_ratio_start !! MATURITY_CONVEYOR : rapport de DIAMETRE max_dia/borne a partir
                                            !! duquel les classes 2 et 3 montent (-). Distinct du critere de
                                            !! classe 1 : la bascule 1->2 est une fermeture de couvert, donc
                                            !! une affaire de hauteur ; au-dela c'est une dimension
                                            !! d'exploitabilite, qui se juge en diametre.
!$OMP THREADPRIVATE(ac12_dia_ratio_start)
  REAL(r_std), SAVE :: ac12_dia_ratio_last  !! MATURITY_CONVEYOR : meme rapport, mais pour le DERNIER passage,
                                            !! celui qui alimente la classe terminale (-). Elle seule a pour
                                            !! receveur l'ensemble coupable : remplie au 31 decembre puis rasee
                                            !! dans la meme sequence, son bilan est negatif si la cadence
                                            !! d'alimentation ne depasse pas celle de la coupe.
                                            !! <= 0 => repli sur ac12_dia_ratio_start (bit-neutre).
!$OMP THREADPRIVATE(ac12_dia_ratio_last)
  LOGICAL, SAVE :: ok_maturity_transfer     !! MATURITY_TRANSFER : the moving stand is an intact COPY of the
                                            !! donor (no circumference-class selection), its rate comes from
                                            !! MAT_F_REF, and the receiver merges by diameter proximity
                                            !! instead of re-drawing the Weibull. .FALSE. => historical paths.
!$OMP THREADPRIVATE(ok_maturity_transfer)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: mat_f_ref !! MATURITY_TRANSFER : reference promotion rate per age
                                                            !! class (fraction of donor area per year). Allocated to
                                                            !! nagec AFTER NAGEC is read; the terminal entry is unused
                                                            !! (that class leaves through the harvest, not upward).
!$OMP THREADPRIVATE(mat_f_ref)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: mat_target_frac !! MATURITY_TRANSFER : target area share of each age
                                                                  !! class within its group (sums to one). The regulator
                                                                  !! measures filling ratios against it. Allocated to
                                                                  !! nagec AFTER NAGEC is read.
!$OMP THREADPRIVATE(mat_target_frac)
  REAL(r_std), SAVE :: mat_band_low          !! MATURITY_TRANSFER : lower edge of the regulator dead band on the
                                             !! receiver filling ratio; inside the band f = MAT_F_REF (-)
!$OMP THREADPRIVATE(mat_band_low)
  REAL(r_std), SAVE :: mat_band_high         !! MATURITY_TRANSFER : upper edge of the dead band (-)
!$OMP THREADPRIVATE(mat_band_high)
  REAL(r_std), SAVE :: mat_gain              !! MATURITY_TRANSFER : exponent of the cross correction
                                             !! (r_don/r_rec)**gain outside the dead band; 0 = regulator off (-)
!$OMP THREADPRIVATE(mat_gain)
  REAL(r_std), SAVE :: mat_f_min             !! MATURITY_TRANSFER : floor of the regulated rate (1/year)
!$OMP THREADPRIVATE(mat_f_min)
  REAL(r_std), SAVE :: mat_f_max             !! MATURITY_TRANSFER : ceiling of the regulated rate (1/year)
!$OMP THREADPRIVATE(mat_f_max)
  REAL(r_std), SAVE :: mat_shift             !! MATURITY_TRANSFER : maturity step of the moving stand, as a
                                             !! relative DIAMETER offset along the self-thinning trajectory.
                                             !! Creates the gap between adjacent classes that the area
                                             !! regulator erodes; 0 = intact copy (no gap) (-)
!$OMP THREADPRIVATE(mat_shift)
  LOGICAL, SAVE :: ok_harvest_feed_coupling !! HARVEST_FEED_COUPLING : the progressive-harvest quota is computed
                                            !! PER SPECIES GROUP and its demand is modulated by the mature-class
                                            !! filling ratio through the MAT_BAND dead band (design doc 4.4).
                                            !! .FALSE. => historical pixel-level quota (bit-neutral).
!$OMP THREADPRIVATE(ok_harvest_feed_coupling)
  REAL(r_std), SAVE :: harv_coupling_gain    !! HARVEST_FEED_COUPLING : exponent of the demand modulation
                                             !! g = r4**gain outside the dead band; 0 = per-group quota
                                             !! without modulation (-)
!$OMP THREADPRIVATE(harv_coupling_gain)
  ! Guillaume M. -- Well-posedness bounds of the demand modulation, NOT calibration
  ! targets: contracts and sanitary fellings prevent a full stop, mobilization
  ! capacity prevents a doubling.
  REAL(r_std), PARAMETER :: harv_g_min = 0.5_r_std  !! floor of the demand modulation g (-)
  REAL(r_std), PARAMETER :: harv_g_max = 1.5_r_std  !! ceiling of the demand modulation g (-)
  ! Guillaume M. -- PHANTOM COHORT PURGE (always-on guard, design doc 4.5): a floor-area
  ! slot whose QMD exceeds this factor times largest_tree_dia is not a stand but a
  ! numerical runaway (per-tree growth in 1/n); its biomass goes to litter and the slot
  ! is reset. Named constant, not a knob: consistent with validation criterion 5.
  REAL(r_std), PARAMETER :: phantom_dia_factor = 2.0_r_std  !! QMD/largest_tree_dia ratio above which a floor-area slot is purged (-)
  ! Guillaume M. -- /!\ Residual slots FLOAT slightly above min_vegfrac, they do not sit
  ! on it (measured 1.009e-5 for a 1e-5 floor: a <= min_vegfrac test misses them by
  ! 0.9 percent). Real stands are >= 100x the floor; 10x separates the two populations.
  REAL(r_std), PARAMETER :: phantom_area_factor = 10.0_r_std !! veget_max/min_vegfrac ratio below which a slot counts as floor-area (-)
  REAL(r_std), SAVE, DIMENSION(nmiclass) :: mi_rotation_factor !! Multiplicateur de la rotation de reference, par classe d'intensite (-).
                                                               !! STRICTEMENT symetrique de mi_dia_factor :
                                                               !!   diametre de coupe = largest_tree_dia(pft) * mi_dia_factor(c)
                                                               !!   rotation          = rotation_ref(pft)     * mi_rotation_factor(c)
                                                               !! Actif seulement si rotation_ref(pft) > 0 ; sinon la rotation
                                                               !! reste indefinie (zero), valeur qu'aucun arbre ne lit.
!$OMP THREADPRIVATE(mi_rotation_factor)
  REAL(r_std), SAVE          :: mi_edge_shape_coef         !! k dans alpha_c = k/sqrt(S_c) ; cale sur la lisiere EFDA (alpha pan-EU moyen 0.0219)
!$OMP THREADPRIVATE(mi_edge_shape_coef)
  REAL(r_std), SAVE          :: mi_age_mod_strength        !! Force de la modulation de S par la concentration d'age (HHI) ; 0 = desactivee
!$OMP THREADPRIVATE(mi_age_mod_strength)
  REAL(r_std), SAVE          :: mi_clearcut_size_min       !! Borne basse de S apres modulation (m2)
!$OMP THREADPRIVATE(mi_clearcut_size_min)
  REAL(r_std), SAVE          :: mi_clearcut_size_max       !! Borne haute de S apres modulation (m2)
!$OMP THREADPRIVATE(mi_clearcut_size_max)
  REAL(r_std), SAVE          :: pest_biomass_ref           !! Standing biomass density used to convert beetle kill to area (gC/m²)
!$OMP THREADPRIVATE(pest_biomass_ref)
  REAL(r_std), SAVE          :: aed_beetle_epidemic_threshold !! Mass-attack index (0-1) above which beetle kills feed the AED edge budget (epidemic phase); below = endemic / background mortality, excluded
!$OMP THREADPRIVATE(aed_beetle_epidemic_threshold)
  INTEGER(i_std), SAVE       :: n_recruit                  !! AED_REGROWTH: years to canopy reclosure (scalar fallback / non-forest cells)
!$OMP THREADPRIVATE(n_recruit)
  INTEGER(i_std), SAVE       :: n_recruit_max              !! AED_REGROWTH: max conveyor length = max recruit time over PFTs (dim of class0_area)
!$OMP THREADPRIVATE(n_recruit_max)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: recruit_time_pft   !! AED_REGROWTH: per-PFT canopy-reclosure time (yr), mapped from MTC
!$OMP THREADPRIVATE(recruit_time_pft)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:) :: n_recruit_eff   !! AED_REGROWTH: per-point effective graduation bin (veget_max-weighted)
!$OMP THREADPRIVATE(n_recruit_eff)
  LOGICAL, SAVE              :: ok_dia_stagger             !! Cold-start: establish each age class at its own diameter bin (staggered qmd)
!$OMP THREADPRIVATE(ok_dia_stagger)
  ! /!\ dia_inv_file, dia_inv_dmax_scale/min/max et dia_inv_dbar_default SUPPRIMES le
  ! 2026-07-30 : l'echelle de maturite est entierement derivee de LARGEST_TREE_DIA par PFT,
  ! la carte externe de diametre moyen D_bar n'a plus de role. Voir CALL_GRAPH.md.

  !
  ! lpj_gap.f90
  !

  ! 1. Scalar

  REAL(r_std), SAVE :: ref_greff                    !! Asymptotic maximum mortality rate
                                                    !! @tex $(year^{-1})$ @endtex
!$OMP THREADPRIVATE(ref_greff)

  !               
  ! lpj_light.f90 
  !              

  ! 1. Scalar
  
  LOGICAL, SAVE :: annual_increase                  !! for diagnosis of fpc increase, compare today's fpc to last year's maximum (T) or
                                                    !! to fpc of last time step (F)? (true/false)
!$OMP THREADPRIVATE(annual_increase)
  REAL(r_std), SAVE :: min_cover                    !! For trees, minimum fraction of crown area occupied
                                                    !! (due to its branches etc.) (0-1, unitless)
                                                    !! This means that only a small fraction of its crown area
                                                    !! can be invaded by other trees.
!$OMP THREADPRIVATE(min_cover)
  !
  ! lpj_pftinout.f90 
  !

  ! 1. Scalar

  REAL(r_std), SAVE :: min_avail                    !! minimum availability
!$OMP THREADPRIVATE(min_avail)
  REAL(r_std), SAVE :: ind_0                        !! initial density of individuals
!$OMP THREADPRIVATE(ind_0)
  ! 3. Coefficients of equations
  
  REAL(r_std), SAVE :: RIP_time_min                 !! test whether the PFT has been eliminated lately (years)
!$OMP THREADPRIVATE(RIP_time_min)
  REAL(r_std), SAVE :: npp_longterm_init            !! Initialisation value for npp_longterm (gC.m^{-2}.year^{-1})
!$OMP THREADPRIVATE(npp_longterm_init)
  REAL(r_std), SAVE :: everywhere_init              !!
!$OMP THREADPRIVATE(everywhere_init)




  !
  ! stomate_data.f90 
  !

  ! 1. Scalar 

  ! 1.2 climatic parameters 

  REAL(r_std), SAVE :: precip_crit               !! minimum precip, in (mm/year)
!$OMP THREADPRIVATE(precip_crit)
  REAL(r_std), SAVE :: gdd_crit_estab            !! minimum gdd for establishment of saplings
!$OMP THREADPRIVATE(gdd_crit_estab)
  REAL(r_std), SAVE :: fpc_crit                  !! critical fpc, needed for light competition and establishment (0-1, unitless)
!$OMP THREADPRIVATE(fpc_crit)

  ! 1.3 sapling characteristics

  REAL(r_std), SAVE :: alpha_grass               !! alpha coefficient for grasses (unitless)
!$OMP THREADPRIVATE(alpha_grass)
  REAL(r_std), SAVE :: alpha_tree                !! alpha coefficient for trees (unitless)
!$OMP THREADPRIVATE(alpha_tree)
  REAL(r_std), SAVE :: struct_to_leaves          !! Fraction of structural carbon in grass and crops as a share of the leaf
                                                 !! carbon pool.
!$OMP THREADPRIVATE(struct_to_leaves)

  REAL(r_std), SAVE :: labile_to_total           !! Fraction of the labile pool in trees, grasses and crops 
!$OMP THREADPRIVATE(labile_to_total)



  ! 1.4  time scales for phenology and other processes (in days)
  REAL(r_std), SAVE :: tau_week                  !! (days)  
!$OMP THREADPRIVATE(tau_week)
  REAL(r_std), SAVE :: tau_hum_month             !! (days)       
!$OMP THREADPRIVATE(tau_hum_month)
  REAL(r_std), SAVE :: tau_precip_month             !! (days)       
!$OMP THREADPRIVATE(tau_precip_month)
  REAL(r_std), SAVE :: tau_hum_week              !! (days)  
!$OMP THREADPRIVATE(tau_hum_week)
  REAL(r_std), SAVE :: tau_t2m_month             !! (days)      
!$OMP THREADPRIVATE(tau_t2m_month)
  REAL(r_std), SAVE :: tau_t2m_week              !! (days)  
!$OMP THREADPRIVATE(tau_t2m_week)
  REAL(r_std), SAVE :: tau_tsoil_month           !! (days)     
!$OMP THREADPRIVATE(tau_tsoil_month)
  REAL(r_std), SAVE :: tau_gpp_week              !! (days)  
!$OMP THREADPRIVATE(tau_gpp_week)
  REAL(r_std), SAVE :: tau_gpp_year              !! (days)  
!$OMP THREADPRIVATE(tau_gpp_year)
  REAL(r_std), SAVE :: tau_gpp_self_thin         !! (years)  
!$OMP THREADPRIVATE(tau_gpp_self_thin)
  REAL(r_std), SAVE :: tau_sugarload_week        !! (days)
!$OMP THREADPRIVATE(tau_sugarload_week)
  REAL(r_std), SAVE :: tau_gdd                   !! (days)  
!$OMP THREADPRIVATE(tau_gdd)
  REAL(r_std), SAVE :: tau_ngd                   !! (days)  
!$OMP THREADPRIVATE(tau_ngd)
  REAL(r_std), SAVE :: coeff_tau_longterm        !! (unitless)
!$OMP THREADPRIVATE(coeff_tau_longterm)
  !! Moyenne temporelle des precipitations totales (pluie + neige) du forcage, par point
  !! de terre local, en mm/jour. Remplie par forcing_meanprecip (mode force uniquement),
  !! consommee par slowproc_init au DEMARRAGE A FROID pour initialiser precip_longterm
  !! maille par maille. NON allouee en mode couple : slowproc_init retombe alors sur le
  !! scalaire PRECIP_LONGTERM_INIT. Le test ALLOCATED() EST le drapeau.
  !! Raison d'etre : sans cela precip_longterm vaut zero au premier pas de temps, pipe_tune2
  !! est rabattu sur min_pipe_tune2, et prescribe dimensionne les peuplements dans une
  !! allometrie qui n'a pas cours -> RDI d'etablissement divise par ~2,25 (2026-07-29).
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:) :: precip_forcing_mean   !! (mm/jour)
!$OMP THREADPRIVATE(precip_forcing_mean)
  REAL(r_std), SAVE :: tau_longterm_max          !! (days)  
!$OMP THREADPRIVATE(tau_longterm_max)

  ! 3. Coefficients of equations

  REAL(r_std), SAVE :: bm_sapl_carbres           !!
!$OMP THREADPRIVATE(bm_sapl_carbres)
  REAL(r_std), SAVE :: bm_sapl_labile            !!
!$OMP THREADPRIVATE(bm_sapl_labile)
  REAL(r_std), SAVE :: bm_sapl_sapabove          !!
!$OMP THREADPRIVATE(bm_sapl_sapabove)
  REAL(r_std), SAVE :: bm_sapl_heartabove        !!
!$OMP THREADPRIVATE(bm_sapl_heartabove)
  REAL(r_std), SAVE :: bm_sapl_heartbelow        !!
!$OMP THREADPRIVATE(bm_sapl_heartbelow)
  REAL(r_std), SAVE :: init_sapl_mass_leaf_nat   !!
!$OMP THREADPRIVATE(init_sapl_mass_leaf_nat)
  REAL(r_std), SAVE :: init_sapl_mass_leaf_agri  !!
!$OMP THREADPRIVATE(init_sapl_mass_leaf_agri)
  REAL(r_std), SAVE :: init_sapl_mass_carbres    !!
!$OMP THREADPRIVATE(init_sapl_mass_carbres)
  REAL(r_std), SAVE :: init_sapl_mass_labile     !!
!$OMP THREADPRIVATE(init_sapl_mass_labile)
  REAL(r_std), SAVE :: init_sapl_mass_root       !!
!$OMP THREADPRIVATE(init_sapl_mass_root)
  REAL(r_std), SAVE :: init_sapl_mass_fruit      !!  
!$OMP THREADPRIVATE(init_sapl_mass_fruit)
  REAL(r_std), SAVE :: cn_sapl_init              !!
!$OMP THREADPRIVATE(cn_sapl_init)
  REAL(r_std), SAVE :: migrate_tree              !!
!$OMP THREADPRIVATE(migrate_tree)
  REAL(r_std), SAVE :: migrate_grass             !!
!$OMP THREADPRIVATE(migrate_grass)
  REAL(r_std), SAVE :: lai_initmin_tree          !!
!$OMP THREADPRIVATE(lai_initmin_tree)
  REAL(r_std), SAVE :: lai_initmin_grass         !!
!$OMP THREADPRIVATE(lai_initmin_grass)
  REAL(r_std), SAVE, DIMENSION(2) :: dia_coeff   !!
!$OMP THREADPRIVATE(dia_coeff)
  REAL(r_std), SAVE, DIMENSION(4) :: bm_sapl_leaf!!
!$OMP THREADPRIVATE(bm_sapl_leaf)


  !
  ! stomate_litter.f90 
  !

  ! 0. Constants

  REAL(r_std), PARAMETER :: Q10 = 10.            !!

  ! 1. Scalar

  REAL(r_std), SAVE :: z_decomp                  !!  Maximum depth for soil decomposer's activity (m)
!$OMP THREADPRIVATE(z_decomp)

  ! 2. Arrays

  REAL(r_std), SAVE :: frac_soil_struct_sua      !! corresponding to frac_soil(istructural,isurface,iabove) 
!$OMP THREADPRIVATE(frac_soil_struct_sua)
  REAL(r_std), SAVE :: frac_soil_struct_aa       !! corresponding to frac_soil(istructural,iactive,iabove) 
!$OMP THREADPRIVATE(frac_soil_struct_aa)
 REAL(r_std), SAVE :: frac_soil_struct_ab        !! corresponding to frac_soil(istructural,iactive,ibelow)
!$OMP THREADPRIVATE(frac_soil_struct_ab)
  REAL(r_std), SAVE :: frac_soil_struct_sa       !! corresponding to frac_soil(istructural,islow,iabove)
!$OMP THREADPRIVATE(frac_soil_struct_sa)
  REAL(r_std), SAVE :: frac_soil_struct_sb       !! corresponding to frac_soil(istructural,islow,ibelow)
!$OMP THREADPRIVATE(frac_soil_struct_sb)
  REAL(r_std), SAVE :: frac_soil_metab_sua       !! corresponding to frac_soil(imetabolic,isurface,iabove)
!$OMP THREADPRIVATE(frac_soil_metab_sua)
  REAL(r_std), SAVE :: frac_soil_metab_aa        !! corresponding to frac_soil(imetabolic,iactive,iabove)
!$OMP THREADPRIVATE(frac_soil_metab_aa)
  REAL(r_std), SAVE :: frac_soil_metab_ab        !! corresponding to frac_soil(imetabolic,iactive,ibelow)
!$OMP THREADPRIVATE(frac_soil_metab_ab)
  REAL(r_std), SAVE :: frac_woody                !! Coefficient for determining the lignin fraction of woody litter
!$OMP THREADPRIVATE(frac_woody)
  REAL(r_std), SAVE :: frac_snag                 !! Coefficient for determining the lignin fraction of snag litter (unitless)
!$OMP THREADPRIVATE(frac_snag)
  REAL(r_std), SAVE, DIMENSION(nparts) :: CN_fix !! C/N ratio of each plant pool (0-100, unitless)

!$OMP THREADPRIVATE(CN_fix)

  ! 3. Coefficients of equations

  REAL(r_std), SAVE :: metabolic_ref_frac        !! used by litter and soilcarbon (0-1, unitless)
!$OMP THREADPRIVATE(metabolic_ref_frac)
  REAL(r_std), SAVE :: metabolic_LN_ratio        !! (0-1, unitless)   
!$OMP THREADPRIVATE(metabolic_LN_ratio)
  ! Exponential decay rates (yr-1) - From Parton et al., 1993
  REAL(r_std), SAVE :: turn_metabolic            !! Exponential decay rate of metabolic litter pool [y-1]
!$OMP THREADPRIVATE(turn_metabolic)
  REAL(r_std), SAVE :: turn_struct               !! Exponential decay rate of structural litter pool [y-1]
!$OMP THREADPRIVATE(turn_struct)
  REAL(r_std), SAVE :: turn_woody                !! unknown where this one comes from (not Parton1993)
!$OMP THREADPRIVATE(turn_woody)
  REAL(r_std), SAVE :: turn_snag                 !! Exponential decay rate of snag pool to soil organic matter
                                                 !! [years-1] This one derives from Boulanger et al. 2011
!$OMP THREADPRIVATE(turn_snag)
  REAL(r_std), SAVE :: turn_snag_fall            !! Exponential fall rate of snag pool to woody litter pool
                                                 !! [years-1]
!$OMP THREADPRIVATE(turn_snag_fall)
  REAL(r_std), SAVE :: snag_frac                 !! Proportion of the woody litter input going into snag 
                                                 !! the rest going in the iwoody pool, drawn from a synthesis : Malysheva et al. 2019 over Russia,
                                                 !! Motta et al. 2006 over Italian Alps, Rahman et al. 2008 over an austrian forest, Caron et al. 2009 over Fennoscandia
                                                  
!$OMP THREADPRIVATE(snag_frac)
  REAL(r_std), SAVE :: soil_Q10                  !!= ln 2
!$OMP THREADPRIVATE(soil_Q10)
  REAL(r_std), SAVE :: soil_Q10_uptake           !!
!$OMP THREADPRIVATE(soil_Q10_uptake)
  REAL(r_std), SAVE :: tsoil_ref                 !!
!$OMP THREADPRIVATE(tsoil_ref)
  REAL(r_std), SAVE :: litter_struct_coef        !! 
!$OMP THREADPRIVATE(litter_struct_coef)
  REAL(r_std), SAVE, DIMENSION(3) :: moist_coeff !!
!$OMP THREADPRIVATE(moist_coeff)
  REAL(r_std), SAVE :: moistcont_min             !! minimum soil wetness to limit the heterotrophic respiration
!$OMP THREADPRIVATE(moistcont_min)
  REAL(r_std), SAVE :: finerootdepthratio        !! the ratio of fine root to overall root e-folding depth (for C inputs)
!$OMP THREADPRIVATE(finerootdepthratio)
  REAL(r_std), SAVE :: altrootratio              !! the maximum ratio of fine root depth to active layer thickness (for C inputs)
!$OMP THREADPRIVATE(altrootratio)

  !
  ! stomate_lpj.f90
  !

  ! 1. Scalar

  REAL(r_std), SAVE :: frac_turnover_daily       !! (0-1, unitless)
!$OMP THREADPRIVATE(frac_turnover_daily)


  !
  ! stomate_npp.f90 
  !

  ! 1. Scalar

  REAL(r_std), SAVE :: tax_max                   !! Maximum fraction of allocatable biomass used 
                                                 !! for maintenance respiration (0-1, unitless)
!$OMP THREADPRIVATE(tax_max)


  !
  ! stomate_phenology.f90
  !

  ! 1. Scalar

  REAL(r_std), SAVE :: min_growthinit_time       !! minimum time since last beginning of a growing season (days)
!$OMP THREADPRIVATE(min_growthinit_time)

  ! 3. Coefficients of equations
  
  REAL(r_std), SAVE :: gddncd_ref                !! Empirical parameter to determine budbreak for ncdgdd phenology type. Reprensents gdd.
!$OMP THREADPRIVATE(gddncd_ref)
  REAL(r_std), SAVE :: gddncd_curve              !! Empirical parameter to determine budbreak for ncdgdd phenology type
!$OMP THREADPRIVATE(gddncd_curve)
  REAL(r_std), SAVE :: gddncd_offset             !! Empirical parameter to determine budbreak for ncdgdd phenology type
!$OMP THREADPRIVATE(gddncd_offset)


  !
  ! stomate_resp.f90
  !

  ! 3. Coefficients of equations

  REAL(r_std), SAVE :: maint_resp_min_vmax      !!
!$OMP THREADPRIVATE(maint_resp_min_vmax)
  REAL(r_std), SAVE :: maint_resp_coeff         !!
!$OMP THREADPRIVATE(maint_resp_coeff)


  !
  ! stomate_som_dynamics.f90 (in stomate_soilcarbon.f90 or
  ! stomate_soil_carbon_discretization.f90)   
  !

  ! 2. Arrays 

  ! 2.1 Fixed fraction from one pool to another (or to CO2 emission)

  REAL(r_std), SAVE :: active_to_pass_ref_frac  !! from active pool: depends on clay content  (0-1, unitless)
                                                !! corresponding to frac_carb(:,iactive,ipassive)
!$OMP THREADPRIVATE(active_to_pass_ref_frac)
  REAL(r_std), SAVE :: surf_to_slow_ref_frac    !! from surface pool
                                                !! corresponding to frac_carb(:,isurf,islow)
!$OMP THREADPRIVATE(surf_to_slow_ref_frac)

  REAL(r_std), SAVE :: active_to_CO2_ref_frac   !! from active pool: depends on clay content  (0-1, unitless)
                                                !! corresponding to frac_resp(:,iactive)
!$OMP THREADPRIVATE(active_to_CO2_ref_frac)
  REAL(r_std), SAVE :: slow_to_pass_ref_frac    !! from slow pool: depends on clay content  (0-1, unitless) 
                                                !! corresponding to frac_carb(:,islow,ipassive)
!$OMP THREADPRIVATE(slow_to_pass_ref_frac)
  REAL(r_std), SAVE :: slow_to_CO2_ref_frac     !! from slow pool (0-1, unitless) 
                                                !! corresponding to frac_resp(:,islow)
!$OMP THREADPRIVATE(slow_to_CO2_ref_frac)
  REAL(r_std), SAVE :: pass_to_active_ref_frac  !! from passive pool (0-1, unitless)
                                                !! corresponding to frac_carb(:,ipassive,iactive)
!$OMP THREADPRIVATE(pass_to_active_ref_frac)
  REAL(r_std), SAVE :: pass_to_slow_ref_frac    !! from passive pool (0-1, unitless)
                                                !! corresponding to frac_carb(:,ipassive,islow)
!$OMP THREADPRIVATE(pass_to_slow_ref_frac)

  ! 2.2 som carbon pools

  REAL(r_std), SAVE :: som_init_active          !! Initial active SOM carbon (g m-2) 
!$OMP THREADPRIVATE(som_init_active)
  REAL(r_std), SAVE :: som_init_slow            !! Initial slow SOM carbon (g m-2)
!$OMP THREADPRIVATE(som_init_slow)
  REAL(r_std), SAVE :: som_init_passive         !! Initial passive SOM carbon (g m-2)
!$OMP THREADPRIVATE(som_init_passive)
  REAL(r_std), SAVE :: som_init_surface         !! Initial surface SOM carbon (g m-2)
!$OMP THREADPRIVATE(som_init_surface)

  ! 3. Define Variable fraction from one pool to another (function of silt and clay fraction)
  REAL(r_std), SAVE :: active_to_pass_clay_frac
!$OMP THREADPRIVATE(active_to_pass_clay_frac)
  REAL(r_std), SAVE :: active_to_CO2_clay_silt_frac
!$OMP THREADPRIVATE(active_to_CO2_clay_silt_frac)
  REAL(r_std), SAVE :: slow_to_pass_clay_frac 
!$OMP THREADPRIVATE(slow_to_pass_clay_frac)

  ! C to N target ratios of differnt pools
  REAL(r_std), SAVE ::  CN_target_iactive_ref   !! CN target ratio of active pool for soil min N = 0
!$OMP THREADPRIVATE(CN_target_iactive_ref)
  REAL(r_std), SAVE ::  CN_target_islow_ref     !! CN target ratio of slow pool for soil min N = 0
!$OMP THREADPRIVATE(CN_target_islow_ref)
  REAL(r_std), SAVE ::  CN_target_ipassive_ref  !! CN target ratio of passive pool for soil min N = 0
!$OMP THREADPRIVATE(CN_target_ipassive_ref)
  REAL(r_std), SAVE ::  CN_target_isurface_ref  !! CN target ratio of surface pool for litter nitrogen content = 0
!$OMP THREADPRIVATE(CN_target_isurface_ref)
  REAL(r_std), SAVE ::  CN_target_iactive_Nmin  !! CN target ratio change per mineral N unit (g m-2) for active pool
!$OMP THREADPRIVATE(CN_target_iactive_Nmin) 
  REAL(r_std), SAVE ::  CN_target_islow_Nmin    !! CN target ratio change per mineral N unit (g m-2) for slow pool 
!$OMP THREADPRIVATE(CN_target_islow_Nmin)
  REAL(r_std), SAVE ::  CN_target_ipassive_Nmin !! CN target ratio change per mineral N unit (g m-2) for passive pool 
!$OMP THREADPRIVATE(CN_target_ipassive_Nmin)
  REAL(r_std), SAVE ::  CN_target_isurface_pnc  !! CN target ratio change per plant nitrogen content unit (%) for surface pool
  !! Turnover in SOM pools (year-1)
!$OMP THREADPRIVATE(CN_target_isurface_pnc)
  REAL(r_std), SAVE :: som_turn_isurface        !! turnover of surface pool (year-1)
!$OMP THREADPRIVATE(som_turn_isurface)
  REAL(r_std), SAVE :: som_turn_iactive         !! turnover of active pool (year-1)
!$OMP THREADPRIVATE(som_turn_iactive)
  REAL(r_std), SAVE :: som_turn_islow           !! turnover of slow pool (year-1)
!$OMP THREADPRIVATE(som_turn_islow)
  REAL(r_std), SAVE :: fslow                    !! convertiing factor to go from active pool turnover to slow pool turnover from Guimberteau et al 2018 GMD
!$OMP THREADPRIVATE(fslow)
  REAL(r_std), SAVE :: fpassive                 !! convertiing factor to go from active pool turnover to passive pool turnover from Guimberteau et al 2018 GMD
!$OMP THREADPRIVATE(fpassive)
  REAL(r_std), SAVE :: stomate_tau 
!$OMP THREADPRIVATE(stomate_tau)
  REAL(r_std), SAVE :: depth_modifier           !! e-folding depth of turnover rates,following Koven et al.,2013,Biogeosciences. A very large value means no depth modification
!$OMP THREADPRIVATE(depth_modifier)
  REAL(r_std), SAVE :: som_turn_iactive_clay_frac !! clay-dependant parameter impacting on turnover rate of active pool 
                                                !! Tm parameter of Parton et al. 1993 (-)
!$OMP THREADPRIVATE(som_turn_iactive_clay_frac)

  !
  ! stomate_turnover.f90
  !

  ! 3. Coefficients of equations

  REAL(r_std), SAVE :: new_turnover_time_ref    !!(days)
!$OMP THREADPRIVATE(new_turnover_time_ref)


  !
  ! stomate_vmax.f90
  !
 
  ! 1. Scalar

  REAL(r_std), SAVE :: vmax_offset              !! minimum leaf efficiency (unitless)
!$OMP THREADPRIVATE(vmax_offset)
  REAL(r_std), SAVE :: leafage_firstmax         !! relative leaf age at which efficiency
                                                !! reaches 1 (unitless)
!$OMP THREADPRIVATE(leafage_firstmax)
  REAL(r_std), SAVE :: leafage_lastmax          !! relative leaf age at which efficiency
                                                !! falls below 1 (unitless)
!$OMP THREADPRIVATE(leafage_lastmax)
  REAL(r_std), SAVE :: leafage_old              !! relative leaf age at which efficiency
                                                !! reaches its minimum (vmax_offset) 
                                                !! (unitless)
!$OMP THREADPRIVATE(leafage_old)
  REAL(r_std), SAVE :: sugar_load_min           !! Lower bound for sugar loading when used to regulate NUE. 
                                                !! Sugar loafding results in a strong reduction of GPP. By
                                                !! taking a rather high lower limit the impact of sugar loading
                                                !! is limited. This will result in high C-reserves unless 
                                                !! leaching is implemented. 
!$OMP THREADPRIVATE(sugar_load_min)
  REAL(r_std), SAVE :: sugar_load_max           !! Upper bound for sugar loading when used to regulate NUE
!$OMP THREADPRIVATE(sugar_load_max)

  !
  ! nitrogen_dynamics (in stomate_soilcarbon.f90) 
  !

  ! 0. Constants
  REAL(r_std), PARAMETER :: D_air = 1.73664     !! Oxygen diffusion rate in the air = 0.07236 m2/h
                                               !! from Table 2 of Li et al, 2000
                                               !! (m**2/day)

  REAL(r_std), PARAMETER :: C_molar_mass = 12  !! Carbon Molar mass (gC mol-1) 

  REAL(r_std), PARAMETER :: Pa_to_hPa    = 0.01      !! Conversion factor from Pa to hPa (-)
  REAL(r_std), PARAMETER :: V_O2         = 0.209476  !! Volumetric fraction of O2 in air (-)

  REAL(r_std), PARAMETER :: pk_NH4 = 9.25      !! The negative logarithm of the acid dissociation constant K_NH4     
                                               !! See Table 4 of Li et al. 1992 and Appendix A of Zhang et al. 2002    

                                     
  ! 1. Scalar

  ! Coefficients for defining maximum porosity
  ! From Saxton, K.E., Rawls, W.J., Romberger, J.S., Papendick, R.I., 1986
  ! Estimationg generalized soil-water characteristics from texture. 
  ! Soil Sci. Soc. Am. J. 50, 1031-1036
  ! Cited in Table 5 (page 444) of
  ! Y. Pachepsky, W.J. Rawls
  ! Development of Pedotransfer Functions in Soil Hydrology
  ! Elsevier, 23 nov. 2004 - 542 pages
  ! http://books.google.fr/books?id=ar_lPXaJ8QkC&printsec=frontcover&hl=fr#v=onepage&q&f=false
  REAL(r_std), SAVE :: h_saxton               !! h coefficient 
!$OMP THREADPRIVATE(h_saxton)
  REAL(r_std), SAVE :: j_saxton               !! j coefficient 
!$OMP THREADPRIVATE(j_saxton)
  REAL(r_std), SAVE :: k_saxton               !! k coefficient
!$OMP THREADPRIVATE(k_saxton)

  ! Values of the power used in the equation defining the diffusion of oxygen in soil
  ! from Table 2 of Li et al, 2000
  REAL(r_std), SAVE :: diffusionO2_power_1    !! (unitless)
!$OMP THREADPRIVATE(diffusionO2_power_1)
  REAL(r_std), SAVE :: diffusionO2_power_2    !! (unitless) 
!$OMP THREADPRIVATE(diffusionO2_power_2)

  ! Temperature-related Factors impacting on Oxygen diffusion rate
  ! From eq. 2 of Table 2 (Li et al, 2000)
  REAL(r_std), SAVE ::   F_nofrost            !! (unitless)
!$OMP THREADPRIVATE(F_nofrost)
  REAL(r_std), SAVE ::   F_frost              !! (unitless)
!$OMP THREADPRIVATE(F_frost)

  ! Coefficients used in the calculation of Volumetric fraction of anaerobic microsites 
  ! a and b constants are not specified in Li et al., 2000
  ! S. Zaehle used a=0.85 and b=1 without mention to any publication
  REAL(r_std), SAVE ::   a_anvf               !! (-)
!$OMP THREADPRIVATE(a_anvf)
  REAL(r_std), SAVE ::   b_anvf               !! (-)
!$OMP THREADPRIVATE(b_anvf)

  ! Coefficients used in the calculation of the Fraction of adsorbed NH4+
  ! Li et al. 1992, JGR, Table 4
  REAL(r_std), SAVE ::   a_FixNH4             !! (-)
!$OMP THREADPRIVATE(a_FixNH4)
  REAL(r_std), SAVE ::   b_FixNH4             !! (-)
!$OMP THREADPRIVATE(b_FixNH4)
  REAL(r_std), SAVE ::   clay_max             !! (-)
!$OMP THREADPRIVATE(clay_max)

  ! Coefficients used in the calculation of the Response of Nitrification
  ! to soil moisture
  ! Zhang et al. 2002, Ecological Modelling, appendix A, page 101
  REAL(r_std), SAVE ::   fw_nit_0             !! (-)
!$OMP THREADPRIVATE(fw_nit_0)
  REAL(r_std), SAVE ::   fw_nit_1             !! (-)
!$OMP THREADPRIVATE(fw_nit_1)
  REAL(r_std), SAVE ::   fw_nit_2             !! (-)
!$OMP THREADPRIVATE(fw_nit_2)
  REAL(r_std), SAVE ::   fw_nit_3             !! (-)
!$OMP THREADPRIVATE(fw_nit_3)
  REAL(r_std), SAVE ::   fw_nit_4             !! (-)
!$OMP THREADPRIVATE(fw_nit_4)

  ! Coefficients used in the calculation of the Response of Nitrification
  ! to Temperature
  ! Zhang et al. 2002, Ecological Modelling, appendix A, page 101
  REAL(r_std), SAVE ::   ft_nit_0             !! (-)
!$OMP THREADPRIVATE(ft_nit_0)
  REAL(r_std), SAVE ::   ft_nit_1             !! (-)
!$OMP THREADPRIVATE(ft_nit_1)
  REAL(r_std), SAVE ::   ft_nit_2             !! (-)
!$OMP THREADPRIVATE(ft_nit_2)
  REAL(r_std), SAVE ::   ft_nit_3             !! (-)
!$OMP THREADPRIVATE(ft_nit_3)
  REAL(r_std), SAVE ::   ft_nit_4             !! (-)
!$OMP THREADPRIVATE(ft_nit_4)

  ! Coefficients used in the calculation of the Response of Nitrification
  ! to pH
  ! Zhang et al. 2002, Ecological Modelling, appendix A, page 101
  REAL(r_std), SAVE ::   fph_0                !! (-)
!$OMP THREADPRIVATE(fph_0)
  REAL(r_std), SAVE ::   fph_1                !! (-)
!$OMP THREADPRIVATE(fph_1)
  REAL(r_std), SAVE ::   fph_2                !! (-)
!$OMP THREADPRIVATE(fph_2)

  ! Coefficients used in the calculation of the response of NO2 or NO 
  ! production during nitrificationof to Temperature
  ! Zhang et al. 2002, Ecological Modelling, appendix A, page 102
  REAL(r_std), SAVE ::   ftv_0                !! (-)
!$OMP THREADPRIVATE(ftv_0)
  REAL(r_std), SAVE ::   ftv_1                !! (-)
!$OMP THREADPRIVATE(ftv_1)
  REAL(r_std), SAVE ::   ftv_2                !! (-) 
!$OMP THREADPRIVATE(ftv_2)

  REAL(r_std), SAVE ::   k_nitrif             !! Nitrification rate at 20 ◦C and field capacity (day-1)
                                              !! Schmid et al., 2001 give a value of 0.2 per day.
                                              !! (https://doi.org/10.1023/A:1012694218748)
                                              !! OCN used a value of 2.0.  We keep the OCN value.
!$OMP THREADPRIVATE(k_nitrif)

  REAL(r_std), SAVE ::   n2o_nitrif_p         !! Reference n2o production per N-NO3 produced g N-N2O  (g N-NO3)-1
                                              !! From Zhang et al., 2002 - Appendix A p. 102
!$OMP THREADPRIVATE(n2o_nitrif_p)
  REAL(r_std), SAVE ::   no_nitrif_p          !! Reference NO production per N-NO3 produced g N-NO  (g N-NO3)-1
                                              !! From Zhang et al., 2002 - Appendix A p. 102
!$OMP THREADPRIVATE(no_nitrif_p)

  ! NO production from chemodenitrification
  ! based on Kesik et al., 2005, Biogeosciences
  ! Coefficients used in the calculation of the Response to Temperature
  REAL(r_std), SAVE ::   chemo_t0             !! (-)
!$OMP THREADPRIVATE(chemo_t0)
  ! Coefficients use in the calculation of the Response to pH
  REAL(r_std), SAVE ::   chemo_ph0            !! (-)
!$OMP THREADPRIVATE(chemo_ph0)
  ! Coefficients used in the calculation of NO production from chemodenitrification
  REAL(r_std), SAVE ::   chemo_0              !! (-)
!$OMP THREADPRIVATE(chemo_0)
  REAL(r_std), SAVE ::   chemo_1              !! (-)
!$OMP THREADPRIVATE(chemo_1)

  ! Denitrification processes
  ! Li et al, 2000, JGR Table 4 eq 1, 2 and 4
  !
  ! Coefficients used in the Temperature response of 
  ! relative growth rate of total denitrifiers - Eq. 2 Table 4 of Li et al., 2000
  REAL(r_std), SAVE ::   ft_denit_0           !! (-)
!$OMP THREADPRIVATE(ft_denit_0)
  REAL(r_std), SAVE ::   ft_denit_1           !! (-)
!$OMP THREADPRIVATE(ft_denit_1)
  REAL(r_std), SAVE ::   ft_denit_2           !! (-)
!$OMP THREADPRIVATE(ft_denit_2)
  !
  ! Coefficients used in the pH response of 
  ! relative growth rate of total denitrifiers - Eq. 2 Table 4 of Li et al., 2000
  REAL(r_std), SAVE ::   fph_no3_0            !! (-) 
!$OMP THREADPRIVATE(fph_no3_0)
  REAL(r_std), SAVE ::   fph_no3_1            !! (-)
!$OMP THREADPRIVATE(fph_no3_1)
  REAL(r_std), SAVE ::   fph_no_0             !! (-)
!$OMP THREADPRIVATE(fph_no_0)
  REAL(r_std), SAVE ::   fph_no_1             !! (-)
!$OMP THREADPRIVATE(fph_no_1)
  REAL(r_std), SAVE ::   fph_n2o_0            !! (-)
!$OMP THREADPRIVATE(fph_n2o_0)
  REAL(r_std), SAVE ::   fph_n2o_1            !! (-)
!$OMP THREADPRIVATE(fph_n2o_1)

  REAL(r_std), SAVE ::   Kn                   !! Half Saturation of N oxydes (kgN/m3)
                                              !! Table 4 of Li et al., 2000
!$OMP THREADPRIVATE(Kn)

  REAL(r_std), SAVE ::   cte_bact
!$OMP THREADPRIVATE(cte_bact)

  ! Maximum Relative growth rate of Nox denitrifiers
  ! Eq.1 Table 4 Li et al., 2000 
  REAL(r_std), SAVE ::   mu_no3_max          !! (hour-1)
!$OMP THREADPRIVATE(mu_no3_max)
  REAL(r_std), SAVE ::   mu_no_max           !! (hour-1)
!$OMP THREADPRIVATE(mu_no_max)
  REAL(r_std), SAVE ::   mu_n2o_max          !! (hour-1)
!$OMP THREADPRIVATE(mu_n2o_max)

  ! Maximum growth yield of NOx denitrifiers on N oxydes
  ! Table 4 Li et al., 2000
  REAL(r_std), SAVE ::   Y_no3               !! (kgC / kgN)
!$OMP THREADPRIVATE(Y_no3)
  REAL(r_std), SAVE ::   Y_no                !! (kgC / kgN)
!$OMP THREADPRIVATE(Y_no)
  REAL(r_std), SAVE ::   Y_n2o               !! (kgC / kgN)
!$OMP THREADPRIVATE(Y_n2o)

  ! Maintenance coefficient on N oxyde
  ! Table 4 Li et al., 2000
  REAL(r_std), SAVE ::   M_no3               !! (kgN / kgC / hour)
!$OMP THREADPRIVATE(M_no3)
  REAL(r_std), SAVE ::   M_no                !! (kgN / kgC / hour)
!$OMP THREADPRIVATE(M_no)
  REAL(r_std), SAVE ::   M_n2o               !! (kgN / kgC / hour)
!$OMP THREADPRIVATE(M_n2o)

        
  REAL(r_std), SAVE ::   Maint_c             !! Maintenance coefficient of carbon (kgC/kgC/h)
                                             !! Table 4 Li et al., 2000
!$OMP THREADPRIVATE(Maint_c)
  REAL(r_std), SAVE ::   Yc                  !! Maximum growth yield on soluble carbon (kgC/kgC)
                                             !! Table 4 Li et al., 2000
!$OMP THREADPRIVATE(Yc)

  !! Coefficients used in the eq. defining the response of N-emission to clay fraction (-)
  !! from  Table 4, Li et al. 2000
  REAL(r_std), SAVE ::   F_clay_0   
!$OMP THREADPRIVATE(F_clay_0) 
  REAL(r_std), SAVE ::   F_clay_1 
!$OMP THREADPRIVATE(F_clay_1)


  REAL(r_std), SAVE ::   ratio_nh4_fert      !! Proportion of ammonium in the fertilizers (ammo-nitrate) 
                                           
!$OMP THREADPRIVATE(ratio_nh4_fert)

  REAL(r_std), SAVE ::   boost_fert      !!  
!$OMP THREADPRIVATE(boost_fert)


  CHARACTER(LEN=30),SAVE :: manure_type      !! type of manure based on Fuchs et al. 2014 table 3-1-1

!$OMP THREADPRIVATE(manure_type)

  REAL(r_std), SAVE ::   cn_ratio_manure     !! C:N ratio of organic fertilizer (average over table 3-1-1 from Fuchs et al. 2014) 
 
!$OMP THREADPRIVATE(cn_ratio_manure)

  ! 2. Arrays
  REAL(r_std), SAVE, DIMENSION(2) :: K_N_min !! [NH4+] (resp. [NO3-]) for which the Nuptake 
                                             !! equals vmax/2.   (umol per litter)
                                             !! from Kronzucker, 1995
!$OMP THREADPRIVATE(K_N_min)

  REAL(r_std), SAVE, DIMENSION(2) :: low_K_N_min !! Rate of N uptake not associated with 
                                             !! Michaelis- Menten Kinetics for Ammonium 
                                             !! (ind.1) and Nitrate (ind.2)
                                             !! from Kronzucker, 1995 ((umol)-1)
!$OMP THREADPRIVATE(low_K_N_min)

  REAL(r_std), SAVE      :: emm_fac          !! Factor for reducing NH3 emission  
!$OMP THREADPRIVATE(emm_fac)
  REAL(r_std), SAVE      :: fact_kn_no       !! Factor for adusting kn constant for NOx production
!$OMP THREADPRIVATE(fact_kn_no)
  REAL(r_std), SAVE      :: fact_kn_n2o      !! Factor for adusting kn constant for N2O production
!$OMP THREADPRIVATE(fact_kn_n2o)
  REAL(r_std), SAVE      :: kfwdenit         !! Factor for adjusting sensitivity of denitrification to water content                     
!$OMP THREADPRIVATE(kfwdenit)  
 REAL(r_std), SAVE      :: fwdenitfc         !! Value at field capacity of the sensitivity function of denitrification to water content                     
!$OMP THREADPRIVATE(fwdenitfc)          
  REAL(r_std), SAVE      :: fracn_drainage   !! Fraction of NH3/NO3 loss by drainage 
!$OMP THREADPRIVATE(fracn_drainage)
  REAL(r_std), SAVE      :: fracn_runoff     !! Fraction of NH3/NO3 loss by runoff
!$OMP THREADPRIVATE(fracn_runoff)


  !! Other N-related parameters
  REAL(r_std), SAVE      :: Dmax             !! Maximal elasticity of foliage N concentrations (Zaehle et al 2010, SI2, Table 1)
!$OMP THREADPRIVATE(Dmax)

  REAL(r_std), SAVE :: reserve_time_tree     !! Maximum number of days during which
                                             !! carbohydrate reserve may be used for 
                                             !! trees (days)
!$OMP THREADPRIVATE(reserve_time_tree)
  
  REAL(r_std), SAVE :: reserve_time_grass    !! Maximum number of days during which
                                             !! carbohydrate reserve may be used for 
                                             !! grasses (days)
!$OMP THREADPRIVATE(reserve_time_grass)
  REAL(r_std), SAVE :: p_n_uptake            !! The minimum correction factor for nitrogen uptake in case of 
                                             !! enough nitrogen in the reserve 
!$OMP THREADPRIVATE(p_n_uptake)


  !! Scalars for slowproc_impose_site_level_lai (read_lai == 1)
  REAL(r_std), SAVE      :: max_canopy_height_prescribed        !! Maximum height of the canopy
!$OMP THREADPRIVATE(max_canopy_height_prescribed)
  REAL(r_std), SAVE      :: min_canopy_height_prescribed        !! Starting height of the canopy
!$OMP THREADPRIVATE(min_canopy_height_prescribed)
  REAL(r_std), SAVE      :: canopy_density_prescribed           !! Prescribed minimum LAI 
!$OMP THREADPRIVATE(canopy_density_prescribed)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: lvllai_prescribed !! Prescribed LAI per level
!$OMP THREADPRIVATE(lvllai_prescribed)
  REAL(r_std), SAVE      :: maximum_lai                         !! Prescribed maximum LAI 
!$OMP THREADPRIVATE(maximum_lai)
  REAL(r_std), SAVE      :: minimum_lai                         !! Prescribed minimum LAI 
!$OMP THREADPRIVATE(minimum_lai)
  INTEGER(i_std), SAVE   :: phenology_doy_start                 !! Day of Year at which phenology starts (leaves appear)
!$OMP THREADPRIVATE(phenology_doy_start)
  INTEGER(i_std), SAVE   :: senescence_doy_end                  !! Day of Year at which senescence ends (leaves fall)
                                                                !!(between 1 and 365, higher than phenology_doy_start)
!$OMP THREADPRIVATE(senescence_doy_end)
  INTEGER(i_std), SAVE   :: phenology_delta                     !! Delta of Doy between start and end of phenology
!$OMP THREADPRIVATE(phenology_delta)
  INTEGER(i_std), SAVE   :: senescence_delta                    !! Delta of Doy between start and end of senescence
!$OMP THREADPRIVATE(senescence_delta)


  !
  ! stomate_windthrow.f90 / stomtate_pest.f90
  !

  ! 0. Constants
  REAL(r_std), SAVE :: clear_cut_max                  !! The maximum contiguous area allowed to be clearfelled
                                                      !! @tex $(m^{2})$ @endtex
!$OMP THREADPRIVATE(clear_cut_max)
  REAL(r_std), PARAMETER :: c_surface = 0.003         !! Surface Drag Coefficient (Raupach 1994) (unitless)
  REAL(r_std), PARAMETER :: c_drag = 0.3              !! Element Drag Coefficient (Raupach 1994) (unitless)
  REAL(r_std), PARAMETER :: c_displacement = 7.5      !! Used by Raupach to calculate the zero-plane displacement (Raupach 1994) (unitless)
  REAL(r_std), PARAMETER :: c_roughness = 2.0         !! Used by Raupach to calculate the surface roughness length (Raupach 1994) (unitless)
  REAL(r_std), PARAMETER :: air_density = 1.2226      !! The value of air density (kg*m-3). If needed, this can be derived dynamically from
                                                      !! other modules of ORCHIDEE, but considering the range of values it can hold, it is probably
                                                      !! not worth additional calculations for being used in WINDTHROW.

  INTEGER(i_std), SAVE :: legacy_years_wind           !! The years used to calculate the maximum/mean/sum related to wind speed
!$OMP THREADPRIVATE(legacy_years_wind)
  INTEGER(i_std), SAVE :: wind_days                   !! The number of timesteps to caluclate mean wind speed to determine calcualte wind damage or not
!$OMP THREADPRIVATE(wind_days)
  
  REAL(r_std), SAVE :: wind_sum_threshold             !! The threshold for the sum of wind_days wind speed for all time step to calculate wind damage
!$OMP THREADPRIVATE(wind_sum_threshold)
  REAL(r_std), SAVE :: wind_ratio_threshold           !! The threshold for the ratio of wind to longterm wind for each time step to sum ratios
!$OMP THREADPRIVATE(wind_ratio_threshold)
  REAL(r_std), SAVE :: wind_ratio_a                    !! JJ 2026: Gompertz asymptote for the storm wind-ratio damage weight
!$OMP THREADPRIVATE(wind_ratio_a)
  REAL(r_std), SAVE :: wind_ratio_b                    !! JJ 2026: Gompertz steepness for the storm wind-ratio damage weight
!$OMP THREADPRIVATE(wind_ratio_b)
  REAL(r_std), SAVE :: wind_ratio_c                    !! JJ 2026: Gompertz inflection for the storm wind-ratio damage weight
!$OMP THREADPRIVATE(wind_ratio_c)
  REAL(r_std), SAVE :: gap_threshold                  !! The threshold for the basal area loss to be considered to create the gap with edges
!$OMP THREADPRIVATE(gap_threshold)
  INTEGER(i_std), SAVE :: legacy_years                !! The years used to calculate the total damage from disturbances with in legacy_years and the default is 3 years.
!$OMP THREADPRIVATE(legacy_years)
  INTEGER(i_std), SAVE :: legacy_years_wood           !! The years used to calculate the wood leftover legacy.
!$OMP THREADPRIVATE(legacy_years_wood)
  INTEGER(i_std), SAVE :: beetle_legacy               !! The years used to calculate wood leftover.
!$OMP THREADPRIVATE(beetle_legacy)

! 1. Scalar
  REAL(r_std), SAVE :: wind_speed_storm_thr           !! the wind speed threshold above which is_storm flag is set to TRUE 
!$OMP THREADPRIVATE(wind_speed_storm_thr)
  REAL(r_std), SAVE :: daily_max_tune                 !! This is a linear tunning factor to adjust the calculated daily maximum wind speed from forcing dataset.
!$OMP THREADPRIVATE(daily_max_tune)
  REAL(r_std), SAVE :: elevate_wind                   !! Height conversion of critical wind speeds above zero-plane displacement
!$OMP THREADPRIVATE(elevate_wind)
  REAL(r_std), SAVE :: h_threshold_wind               !! Ratio of height threshold to pipe_tune2 for truncating wind damage
!$OMP THREADPRIVATE(h_threshold_wind)
  REAL(r_std), SAVE :: d_threshold_wind               !! Fraction of the largest tree diameter used as reference threshold
!$OMP THREADPRIVATE(d_threshold_wind)
  INTEGER(i_std), SAVE :: nb_days_storm               !! the number of days at which the max wind speed is less than wind_speed_storm_thr
!$OMP THREADPRIVATE(nb_days_storm)
  INTEGER(i_std), SAVE :: nb_years_bgi                !! the number of years needed to average the bark beetle generation index 
!$OMP THREADPRIVATE(nb_years_bgi)
  REAL(r_std), SAVE :: replace_threshold_wind         !! the threshold for the ratio of biomss revmoved from storm event
                                                      !! to replace stand to the youngest age class
!$OMP THREADPRIVATE(replace_threshold_wind)
  REAL(r_std), SAVE :: replace_threshold_beetle       !! the threshold for the ratio of biomss revmoved from bark beetle outbreak
                                                      !! to replace stand to the youngest age class
!$OMP THREADPRIVATE(replace_threshold_beetle)

  REAL(r_std), SAVE :: replace_threshold_fire         !! the threshold for the ratio of biomss revmoved from fire
                                                      !! to replace stand to the youngest age class
!$OMP THREADPRIVATE(replace_threshold_fire)

  ! stomate_season.f90 
  !

  ! 1. Scalar

  REAL(r_std), SAVE :: gppfrac_dormance        !! report maximal GPP/GGP_max for dormance (0-1, unitless)
!$OMP THREADPRIVATE(gppfrac_dormance)
  REAL(r_std), SAVE :: tau_climatology         !! tau for "climatologic variables (years)
!$OMP THREADPRIVATE(tau_climatology)
  REAL(r_std), SAVE :: hvc1                    !! parameters for herbivore activity (unitless)
!$OMP THREADPRIVATE(hvc1)
  REAL(r_std), SAVE :: hvc2                    !! parameters for herbivore activity (unitless)
!$OMP THREADPRIVATE(hvc2)
  REAL(r_std), SAVE :: leaf_frac_hvc           !! leaf fraction (0-1, unitless)
!$OMP THREADPRIVATE(leaf_frac_hvc)
  REAL(r_std), SAVE :: tlong_ref_max           !! maximum reference long term temperature (K)
!$OMP THREADPRIVATE(tlong_ref_max)
  REAL(r_std), SAVE :: tlong_ref_min           !! minimum reference long term temperature (K)
!$OMP THREADPRIVATE(tlong_ref_min)

  ! 3. Coefficients of equations

  REAL(r_std), SAVE :: ncd_max_year
!$OMP THREADPRIVATE(ncd_max_year)
  REAL(r_std), SAVE :: ngd_threshold
!$OMP THREADPRIVATE(ngd_threshold)
  REAL(r_std), SAVE :: green_age_ever
!$OMP THREADPRIVATE(green_age_ever)
  REAL(r_std), SAVE :: green_age_dec
!$OMP THREADPRIVATE(green_age_dec)
  
  REAL(r_std), SAVE :: ngd_min_dormance
!$OMP THREADPRIVATE(ngd_min_dormance)

  !
  ! sapiens_forestry.f90
  !

  INTEGER(i_std), SAVE      :: ncirc                      !! Number of circumference classes used to calculate C allocation. This creates within stand heterogeneity (i.e., structured canopies).
!$OMP THREADPRIVATE(ncirc)
  INTEGER(i_std), SAVE      :: nagec                      !! Number of age classes. This creates landscape-level heterogeneity (i.e. the effects of deforestation and harvesting are longer seen in the energy budget)
!$OMP THREADPRIVATE(nagec)
  ! AGE_CLASS_BOUNDS_PFT (Marie 2026, design/MODULE_DESIGN_AGE_CLASS_BOUNDS_PFT.md) :
  ! les bornes deviennent PAR PFT (nagec, nvm). Motif : largest_tree_dia est par PFT
  ! (0.20 m pour l'eucalyptus, 0.38 m ailleurs) alors que les bornes etaient globales,
  ! si bien que la coupe rase tombait en classe 2 ou 3 selon l'essence et que la classe 4
  ! n'etait jamais atteinte. Dans l'approche AED les classes sont FONCTIONNELLES
  ! (classe 1 = ouvert post-perturbation, seule assiette d'A_open ; classe nagec = fin de
  ! rotation), donc la fonction doit etre invariante par PFT.
  ! Gate OFF : toutes les colonnes recoivent le vecteur global => bit-neutre.
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:) :: age_class_bound !! Boundaries of the age classes defined as threshold diameters (nagec, nvm) (m)
!$OMP THREADPRIVATE(age_class_bound)
  LOGICAL, SAVE             :: ok_age_class_bound_pft     !! Derive age class bounds from largest_tree_dia (per PFT) instead of the global vector
!$OMP THREADPRIVATE(ok_age_class_bound_pft)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: age_class_frac !! Age class bounds as a fraction of largest_tree_dia (nagec-1) (-)
!$OMP THREADPRIVATE(age_class_frac)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: age_class_init_frac !! Cold-start area fraction of each age class (nagec) (-)
!$OMP THREADPRIVATE(age_class_init_frac)
  LOGICAL, SAVE             :: ok_clearcut_last_class     !! Management clearcut (CONDITION 2/3) restricted to the last age class
!$OMP THREADPRIVATE(ok_clearcut_last_class)
  LOGICAL, SAVE             :: ok_rdi_band_class1         !! Age class 1 joins the thinned classes and gets the RDI band
!$OMP THREADPRIVATE(ok_rdi_band_class1)
  REAL(r_std), SAVE         :: ac12_dia_ratio_down        !! Demotion threshold as a fraction of age_class_bound(1); 1.0 = no hysteresis (-)
!$OMP THREADPRIVATE(ac12_dia_ratio_down)
  ! Guillaume M. -- ROTATION_GROWTH: weight of the SITE against the published silvicultural
  ! itinerary in the applied rotation. See design/MODULE_DESIGN_ROTATION_GROWTH.md.
  REAL(r_std), SAVE         :: rotation_growth_weight     !! Weight of site conditions in the rotation; 0 = recommendation alone (-)
!$OMP THREADPRIVATE(rotation_growth_weight)
  REAL(r_std), SAVE         :: dia_growth_tau             !! Memory of the diameter-increment integrator (yr)
!$OMP THREADPRIVATE(dia_growth_tau)
  REAL(r_std), SAVE         :: dia_growth_area_tol        !! Relative area change above which a slot's increment is rejected (-)
!$OMP THREADPRIVATE(dia_growth_area_tol)
  REAL(r_std), SAVE         :: rotation_growth_min        !! Lower bound of the growth-derived rotation (yr)
!$OMP THREADPRIVATE(rotation_growth_min)
  REAL(r_std), SAVE         :: rotation_growth_max        !! Upper bound of the growth-derived rotation (yr)
!$OMP THREADPRIVATE(rotation_growth_max)
  ! Bande de RDI cible FIXE pour les classes d'age intermediaires (2026-07-28).
  ! Court-circuite les polynomes calculate_rdi_boundaries, cales par PFT et par TYPE de
  ! gestion mais AVEUGLES a l'intensite : a un diametre donne ils plafonnent trop bas et le
  ! peuplement ne peut jamais se densifier. Calee sur Pucher (RDI observe median 0,686).
  ! Ouverture de la classe 1 a l'etablissement (2026-07-27, demande utilisateur).
  ! Mesure : sur 1109 arrivees d'aire en classe 1, la coupe ne degrade PAS la fermeture
  ! radiative (ecart median +0.037, 66 % des evenements sans ouverture) alors qu'elle
  ! produit bien la discontinuite structurelle (hauteur 13.8 -> 5.1 m). La classe 1 est
  ! donc "jeune" mais pas "ouverte", et le terme AED alpha_reg*A_class1 sur-compte.
  ! Exprime en FRACTION de rdi_max(ivm, forest_managed) et non en RDI absolu : rdi_max
  ! depend de la gestion (0.9 non gere / 0.7 rotationnel), un RDI absolu de 0.8 en
  ! classe 4 depasserait le plafond d'une foret eclaircie qui l'eclaircirait aussitot.
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: prescribe_rdi_frac !! Establishment RDI as a fraction of rdi_max, per age class (nagec) (-)
!$OMP THREADPRIVATE(prescribe_rdi_frac)
  LOGICAL, SAVE             :: ok_min_density_reset       !! Keep the dens_target-driven clearcut (sapiens_forestry CONDITION 1)
!$OMP THREADPRIVATE(ok_min_density_reset)
  LOGICAL, SAVE             :: ok_slow_death_density      !! Keep the dens_target-driven slow death of sparse stands (stomate_kill)
!$OMP THREADPRIVATE(ok_slow_death_density)
  LOGICAL, SAVE             :: ok_aed_edge_rdi_weight     !! Weight the AED age-class-1 edge source by the stand openness (1 - RDI)
!$OMP THREADPRIVATE(ok_aed_edge_rdi_weight)
  ! Verrou 1 de la desynchronisation (2026-07-27). Le split d'aire sur les classes d'age
  ! de slowproc_readvegetmax etait gate sur restname_in=='NONE', c'est-a-dire l'absence de
  ! TOUT restart. Or le spinup frais repart d'un restart PARTIEL (sol conserve, veget_max
  ! largue, carte PFT relue) : la condition etait fausse, la branche all-to-youngest
  ! s'executait, et 100 % de l'aire atterrissait en classe 1 (mesure spin51check).
  ! Cold-start structure (2026-07-27, cadrage utilisateur) : UN SEUL champ force de
  ! l'exterieur (les fractions d'aire par classe, issues de Pucher), le reste derive en
  ! interne -- diametre par les bornes de classe (deja par PFT), RDI par
  ! prescribe_rdi_frac, age par la rotation cible.
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: age_class_age_frac !! Stand age of each age class as a fraction of the target rotation age (nagec) (-)
!$OMP THREADPRIVATE(age_class_age_frac)
  REAL(r_std), SAVE :: stand_age_rotation_max                        !! Cap on the rotation used to DERIVE the cold-start stand age only (years)
!$OMP THREADPRIVATE(stand_age_rotation_max)
  LOGICAL, SAVE             :: forced_clear_cut           !! flag at wich the stand is force to clear-cut
!$OMP THREADPRIVATE(forced_clear_cut)
  INTEGER(i_std), SAVE      :: noutdiaclass               !! Number of diameter classes in the output files. ncirc is the number of circumference classes
                                                          !! ORCHIDEE. For tree demography it could be more convenient to get the results in diameter
                                                          !! classes with fixed boundaries.
!$OMP THREADPRIVATE(noutdiaclass)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: out_dia_class !! Boundaries of the diameter classes (no longer circumference as in
                                                          !! ORCHIEE!) (m). Given how the code is setup, this list should start with zero and end with
                                                          !! a very high and unlikely diameter. The length is thus noutdiaclass+1.
!$OMP THREADPRIVATE(out_dia_class)
  
  INTEGER(i_std), SAVE      :: test_pft                   !! Number of PFT for which detailed output 
!$OMP THREADPRIVATE(test_pft)
  INTEGER(i_std), SAVE      :: test_grid                  !! Number of the grid square for which detailed output 
!$OMP THREADPRIVATE(test_grid)
  
  INTEGER(i_std), SAVE      :: ndia_harvest             !! The number of diameter classes used for the wood harvest pools.
!$OMP THREADPRIVATE(ndia_harvest)

  REAL(r_std), SAVE         :: max_harvest_dia          !! The largest diameter for the harvest pools to
                                                        !! keep track of harvested wood from forests.
!$OMP THREADPRIVATE(max_harvest_dia)

  INTEGER(i_std), SAVE      :: n_pai                    !! Number of years used for the calculation of the periodic annual increment
!$OMP THREADPRIVATE(n_pai)

  INTEGER(i_std), SAVE      :: ntrees_profit            !! Number of trees below which the forest will be cut and replanted
!$OMP THREADPRIVATE(ntrees_profit)

  INTEGER, SAVE             :: species_change_force     !! New species after a final cut for testing and debugging only
  !$OMP THREADPRIVATE(species_change_force)

  INTEGER, SAVE             :: fm_change_force          !! New management after a final cut for testing and debugging only
!$OMP THREADPRIVATE(fm_change_force)

 REAL(r_std), SAVE        :: min_water_stress          !! Minimal value for wstress_fac (unitless, 0-1)
!$OMP THREADPRIVATE(min_water_stress)

REAL(r_std), SAVE         :: max_delta_KF              !! Maximum change in KF from one time step to another (m)
                                                       !! This is a bit arbitrary. 
!$OMP THREADPRIVATE(max_delta_KF)

REAL(r_std), SAVE         :: maint_from_gpp            !! Some carbon needs to remain to support the growth, hence, 
                                                       !! respiration will be limited. In this case resp_maint 
                                                       !! (gC m-2 dt-1) should not be more than 80% (::maint_from_gpp) 
                                                       !! of the GPP (gC m-2 s-1)
!$OMP THREADPRIVATE(maint_from_gpp)
  REAL(r_std), PARAMETER :: m2_to_km2 = 1000000.       !! Conversion from m2 to km2 
  REAL(r_std), PARAMETER :: m2_to_ha = 10000.          !! Conversion from m2 to hectares
  REAL(r_std), PARAMETER :: ha_to_m2 = 0.0001          !! Conversion from hectares (forestry) to m2 (rest of the code)
  REAL(r_std), PARAMETER :: m_to_cm = 100.             !! Conversion from m to cm
  REAL(r_std), PARAMETER :: cm_to_m = 0.01             !! Conversion from cm to m
  REAL(r_std), PARAMETER :: peta_to_unit = 1.0E15      !! Convert Peta to unit
  REAL(r_std), PARAMETER :: tera_to_unit = 1.0E12      !! Convert Tera to unit 
  REAL(r_std), PARAMETER :: giga_to_unit = 1.0E09      !! Convert Giga to unit
  REAL(r_std), PARAMETER :: mega_to_unit = 1.0E06      !! Convert Mega to unit
  REAL(r_std), PARAMETER :: kilo_to_unit = 1.0E03      !! Convert Kilo to unit
  REAL(r_std), PARAMETER :: centi_to_unit = 1.0E02     !! Convert centi to unit
  REAL(r_std), PARAMETER :: milli_to_unit = 1.0E-03    !! Convert milli to unit
  REAL(r_std), PARAMETER :: carbon_to_kilo = 2.0E-03   !! Convert g carbon to kilo biomass

  !
  ! Debugging
  !
  INTEGER(i_std), SAVE :: err_act                      !! There are three levels of error checking 
                                                       !! see constantes.f90 for more details
!$OMP THREADPRIVATE(err_act)
  INTEGER(i_std), SAVE :: plev                         !! print level of the subroutine ipslerr_p 
                                                       !! (1:note, 2: warn and 3:stop) 
!$OMP THREADPRIVATE(plev) 
  REAL(r_std), SAVE    :: sync_threshold               !! The threshold above which a warning is generated when the biomass is synchronized
!$OMP THREADPRIVATE(sync_threshold) 

  ! Soil carbon related variables
  LOGICAL, SAVE                   :: soilc_isspinup = .FALSE.            !! Activate soil carbon spinup. This value might be overwritten in force soil - no getin needed.
!$OMP THREADPRIVATE(soilc_isspinup)
  REAL(r_std), PARAMETER          :: O2_init_conc = 298.3                !! gO2/m**3 mean for Cherskii
  REAL(r_std), PARAMETER          :: CH4_init_conc = 0.001267            !! gCH4/m**3 mean for Cherskii
  REAL(r_std), PARAMETER          :: z_root_max = 0.5                    !! Depth at which litter carbon input decays e-fold (root depth); 0.5 for compar w/WH
  REAL(r_std), PARAMETER          :: diffO2_air = 1.596E-5               !! oxygen diffusivity in air (m**2/s)
  REAL(r_std), PARAMETER          :: diffO2_w = 1.596E-9                 !! oxygen diffusivity in water (m**2/s)
  REAL(r_std), PARAMETER          :: O2_surf = 0.209                     !! oxygen concentration in surface air (molar fraction)
  REAL(r_std), PARAMETER          :: diffCH4_air = 1.702E-5              !! methane diffusivity in air (m**2/s)
  REAL(r_std), PARAMETER          :: diffCH4_w = 2.0E-9                  !! methane diffusivity in water (m**2/s)
  REAL(r_std), PARAMETER          :: CH4_surf = 1700.E-9                 !! methane concentration in surface air (molar fraction)
  REAL(r_std), PARAMETER          :: tetasat =  .5                       !! volumetric water content at saturation (porosity)
  REAL(r_std), PARAMETER          :: BunsenO2 = 0.038                    !! Bunsen coefficient for O2 (10C, 1bar)
  REAL(r_std), PARAMETER          :: BunsenCH4 = 0.043                   !! Bunsen coefficient for CH4 (10C, Wiesenburg et Guinasso, Jr., 1979)
  REAL(r_std), PARAMETER          :: ebuthr = 0.9                        !! Soil humidity threshold for ebullition
  REAL(r_std), PARAMETER          :: wCH4 = 16.                          !! molar weight of CH4 (g/mol)
  REAL(r_std), PARAMETER          :: wO2 = 32.                           !! molar weight of O2 (g/mol)
  REAL(r_std), PARAMETER          :: wC = 12.                            !! molar weight of C (g/mol)
  REAL(r_std), PARAMETER          :: avm = .01                           !! minimum air volume (m**3 air/m**3 soil)
  REAL(r_std), PARAMETER          :: hmin_tcalc = .001                   !! minimum total snow layer thickness below which we ignore diffusion across the snow layer
  

  ! Peatland and tides
  !
  LOGICAL :: ok_peat_hydro                   !! flag to activate peat in sechiba
!$OMP THREADPRIVATE(ok_peat_hydro)
  LOGICAL :: tides                        !! flag to activate mangroves, if ok_peat_hydro 
!$OMP THREADPRIVATE(tides)
   LOGICAL, SAVE :: perma_peat            !! flag to activate soil carbon computation over peatland, available if ok_peat_hydro 
!$OMP THREADPRIVATE(perma_peat)
  LOGICAL, SAVE :: agri_peat              !! flag to activate agriculture development on peatland, avaialble if ok_peat_hydro 
!$OMP THREADPRIVATE(agri_peat)
  LOGICAL, SAVE :: ok_peat_NoDiscretisation                !! Activate an ancien module for peatland carbon accumus, available if ok_peat_hydro 
!$OMP THREADPRIVATE(ok_peat_NoDiscretisation)
  LOGICAL, SAVE :: agri_peat_prop         !! To activate different way of agri_peat
!$OMP THREADPRIVATE(agri_peat_prop)
  LOGICAL, SAVE :: agri_peat_MAXcrop      !! To activate different way of agri_peat
!$OMP THREADPRIVATE(agri_peat_MAXcrop)
  LOGICAL, SAVE :: agri_peat_MINcrop      !! To activate different way of agri_peat
!$OMP THREADPRIVATE(agri_peat_MINcrop)
  REAL(r_std), SAVE :: frac1     !!fraction related to soil carbon and decomposition in peatland
!$OMP THREADPRIVATE(frac1)
  REAL(r_std), SAVE :: frac2     !!fraction related to soil carbon and decomposition in peatland
!$OMP THREADPRIVATE(frac2)
  REAL(r_std), SAVE, DIMENSION(3) :: flux_tot_coeff !! coeff related to decomposition rate on peatland, esp agri_peat
!$OMP THREADPRIVATE(flux_tot_coeff)
  
  ! the following variables are related to ok_peat_NoDiscretisation, to be removed later xw
  REAL(r_std), SAVE :: p_A                !!acrotelm density (g/m3)
!$OMP THREADPRIVATE(p_A)
  REAL(r_std), SAVE :: p_C                !!catotelm density (g/m3)
!$OMP THREADPRIVATE(p_C)
  REAL(r_std), SAVE :: cf_A               !!carbon fraction in acrotelm peat
!$OMP THREADPRIVATE(cf_A)
  REAL(r_std), SAVE :: cf_C               !!carbon fraction in catotelm peat
!$OMP THREADPRIVATE(cf_C)
  REAL(r_std), SAVE :: v_ratio            !!ratio anaerobic-aerobic CO2
!$OMP THREADPRIVATE(v_ratio)
  REAL(r_std), SAVE :: KA_ini             !!acrotelm decomposition rate (an-1)
!$OMP THREADPRIVATE(KA_ini )
  REAL(r_std), SAVE :: KP_ini             !!catotelm formation rate (an-1)
!$OMP THREADPRIVATE(KP_ini )
  REAL(r_std), SAVE :: KC_ini             !!catotelm decomposition rate (an-1)
!$OMP THREADPRIVATE(KC_ini )
  


END MODULE constantes_var
