! =================================================================================================================================
! MODULE       : stomate_soil_carbon_discretization
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see
! ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       Calculate permafrost soil carbon dynamics following POPCRAN by Dmitry Khvorstyanov
!!      
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_stomate/stomate_soil_carbon_discretization.f90 $ 
!! $Date: 2026-05-11 11:11:46 +0200 (lun. 11 mai 2026) $
!! $Revision: 9526 $
!! \n
!_
!================================================================================================================================

MODULE stomate_soil_carbon_discretization
  
  ! modules used:
  USE ioipsl_para  
  USE constantes_soil_var
  USE constantes_soil
  USE constantes_var
  USE pft_parameters
  USE vertical_soil
  USE stomate_data
  USE grid
  USE mod_orchidee_para
  USE xios_orchidee
  USE function_library, ONLY : check_mass_balance
  USE time, ONLY: FirstTsYear


  IMPLICIT NONE
  PRIVATE
  PUBLIC stomate_soil_carbon_discretization_deep_somcycle, stomate_soil_carbon_discretization_clear, &
       stomate_soil_carbon_discretization_microactem, calc_vert_int_som

  INTEGER(i_std), SAVE                                 :: printlev_loc   !! Local level of text output for current module
!$OMP THREADPRIVATE(printlev_loc)  
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)         :: zf_soil        !! depths of full levels (m)
!$OMP THREADPRIVATE(zf_soil)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)         :: zi_soil        !! depths of intermediate levels (m)
!$OMP THREADPRIVATE(zi_soil)
  REAL(r_std), SAVE                                    :: mu_soil
!$OMP THREADPRIVATE(mu_soil)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:)     :: alphaO2_soil
!$OMP THREADPRIVATE(alphaO2_soil)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:)     :: betaO2_soil
!$OMP THREADPRIVATE(betaO2_soil)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:)     :: alphaCH4_soil
!$OMP THREADPRIVATE(alphaCH4_soil)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:)     :: betaCH4_soil
!$OMP THREADPRIVATE(betaCH4_soil)
  
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:) 	:: heights_snow      !! total thickness of snow levels (m)
!$OMP THREADPRIVATE(heights_snow)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:,:)     	:: zf_snow           !! depths of full levels (m)
!$OMP THREADPRIVATE(zf_snow)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:,:)     	:: zi_snow           !! depths of intermediate levels (m)
!$OMP THREADPRIVATE(zi_snow)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)     	:: zf_snow_nopftdim  !! depths of full levels (m)
!$OMP THREADPRIVATE(zf_snow_nopftdim)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)     	:: zi_snow_nopftdim  !! depths of intermediate levels (m)
!$OMP THREADPRIVATE(zi_snow_nopftdim)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  :: zf_coeff_snow
!$OMP THREADPRIVATE(zf_coeff_snow)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  :: zi_coeff_snow
!$OMP THREADPRIVATE(zi_coeff_snow)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)    :: mu_snow
!$OMP THREADPRIVATE(mu_snow)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  :: alphaO2_snow
!$OMP THREADPRIVATE(alphaO2_snow)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  :: betaO2_snow
!$OMP THREADPRIVATE(betaO2_snow)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  :: alphaCH4_snow
!$OMP THREADPRIVATE(alphaCH4_snow)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  :: betaCH4_snow
!$OMP THREADPRIVATE(betaCH4_snow)

  REAL(r_std), DIMENSION(:,:,:,:), ALLOCATABLE, SAVE  :: alpha_a
!$OMP THREADPRIVATE(alpha_a)
  REAL(r_std), DIMENSION(:,:,:,:), ALLOCATABLE, SAVE  :: alpha_s
!$OMP THREADPRIVATE(alpha_s)
  REAL(r_std), DIMENSION(:,:,:,:), ALLOCATABLE, SAVE  :: alpha_p
!$OMP THREADPRIVATE(alpha_p)
  REAL(r_std), DIMENSION(:,:),   ALLOCATABLE, SAVE  :: mu_soil_rev
!$OMP THREADPRIVATE(mu_soil_rev)

  REAL(r_std), DIMENSION(:,:,:,:), ALLOCATABLE, SAVE  :: beta_a
!$OMP THREADPRIVATE(beta_a)
  REAL(r_std), DIMENSION(:,:,:,:), ALLOCATABLE, SAVE  :: beta_s
!$OMP THREADPRIVATE(beta_s)
  REAL(r_std), DIMENSION(:,:,:,:), ALLOCATABLE, SAVE  :: beta_p
!$OMP THREADPRIVATE(beta_p)
  LOGICAL, DIMENSION(:,:), ALLOCATABLE, SAVE        :: cryoturb_location
!$OMP THREADPRIVATE(cryoturb_location)
  LOGICAL, DIMENSION(:,:), ALLOCATABLE, SAVE        :: bioturb_location
!$OMP THREADPRIVATE(bioturb_location)
  REAL(r_std), DIMENSION(:,:,:), ALLOCATABLE, SAVE      :: airvol_soil
!$OMP THREADPRIVATE(airvol_soil)
  REAL(r_std), DIMENSION(:,:,:), ALLOCATABLE, SAVE      :: totporO2_soil              !! total oxygen porosity in the soil
!$OMP THREADPRIVATE(totporO2_soil)
  REAL(r_std), DIMENSION(:,:,:), ALLOCATABLE, SAVE      :: totporCH4_soil             !! total methane porosity in the soil
!$OMP THREADPRIVATE(totporCH4_soil)
  REAL(r_std), DIMENSION(:,:,:), ALLOCATABLE, SAVE      :: diffO2_soil                !! oxygen diffusivity in the soil (m**2/s)
!$OMP THREADPRIVATE(diffO2_soil)
  REAL(r_std), DIMENSION(:,:,:), ALLOCATABLE, SAVE      :: diffCH4_soil               !! methane diffusivity in the soil (m**2/s)
!$OMP THREADPRIVATE(diffCH4_soil)
  REAL(r_std), DIMENSION(:,:,:),ALLOCATABLE, SAVE       :: airvol_snow
!$OMP THREADPRIVATE(airvol_snow)
  REAL(r_std), DIMENSION(:,:,:),ALLOCATABLE, SAVE       :: totporO2_snow              !! total oxygen porosity in the snow 
!$OMP THREADPRIVATE(totporO2_snow)
  REAL(r_std), DIMENSION(:,:,:),ALLOCATABLE, SAVE       :: totporCH4_snow             !! total methane porosity in the snow
!$OMP THREADPRIVATE(totporCH4_snow)
  REAL(r_std), DIMENSION(:,:,:),ALLOCATABLE, SAVE       :: conduct_snow
!$OMP THREADPRIVATE(conduct_snow)
  REAL(r_std), DIMENSION(:,:,:), ALLOCATABLE, SAVE      :: diffCH4_snow               !! methane diffusivity in the snow (m**2/s)
!$OMP THREADPRIVATE(diffCH4_snow)
  REAL(r_std), DIMENSION(:,:,:), ALLOCATABLE, SAVE      :: diffO2_snow                !! oxygen diffusivity in the snow (m**2/s)
!$OMP THREADPRIVATE(diffO2_snow)
  REAL(r_std), DIMENSION(:,:), ALLOCATABLE, SAVE        :: alt
!$OMP THREADPRIVATE(alt)
  INTEGER(i_std), DIMENSION(:,:), ALLOCATABLE, SAVE     :: alt_ind                    !! active layer thickness  
!$OMP THREADPRIVATE(alt_ind)
  REAL(r_std), DIMENSION(:,:),ALLOCATABLE, SAVE         :: z_root                     !! Rooting depth
!$OMP THREADPRIVATE(z_root)
  INTEGER(i_std), DIMENSION(:,:),ALLOCATABLE, SAVE      :: rootlev                    !! The deepest model level within the rooting depth
!$OMP THREADPRIVATE(rootlev)
  LOGICAL,DIMENSION(:,:),ALLOCATABLE,  SAVE             :: veget_mask_2d              !! whether there is vegetation 
!$OMP THREADPRIVATE(veget_mask_2d)
  REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:,:), SAVE    :: fc                         !! flux fractions within carbon pools
!$OMP THREADPRIVATE(fc)
  REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:), SAVE      :: fr                         !! fraction of decomposed carbon that goes into the atmosphere
!$OMP THREADPRIVATE(fr)

  REAL(r_std), ALLOCATABLE, DIMENSION(:,:), SAVE        :: deep_rhC_a_pftmean !!rh from active C pool per layer (gC/m**3/s)
  !$OMP THREADPRIVATE(deep_rhC_a_pftmean)
  REAL(r_std), ALLOCATABLE, DIMENSION(:,:), SAVE        :: deep_rhC_s_pftmean !!rh from slow C pool per layer (gC/m**3/s)
  !$OMP THREADPRIVATE(deep_rhC_s_pftmean)
  REAL(r_std), ALLOCATABLE, DIMENSION(:,:), SAVE        :: deep_rhC_p_pftmean !!rh from passive C pool per layer (gC/m**3/s)
  !$OMP THREADPRIVATE(deep_rhC_p_pftmean)
  REAL (r_std), ALLOCATABLE, DIMENSION(:), SAVE         :: rhC_a_pftmean !!rh from active C pool integrated over depth (gC/m**2/s)
  !$OMP THREADPRIVATE(rhC_a_pftmean)
  REAL (r_std), ALLOCATABLE, DIMENSION(:), SAVE         :: rhC_s_pftmean !!rh from slow C pool integrated over depth (gC/m**2/s)
  !$OMP THREADPRIVATE(rhC_s_pftmean)
  REAL (r_std), ALLOCATABLE, DIMENSION(:), SAVE         :: rhC_p_pftmean !!rh from passive C pool integrated over depth (gC/m**2/s)
  !$OMP THREADPRIVATE(rhC_p_pftmean)
  REAL(r_std),ALLOCATABLE,SAVE, DIMENSION(:,:,:)        :: frac_carb          !! Flux fractions between carbon pools 
                                                                              !! (second index=origin, third index=destination) 
                                                                              !! (unitless, 0-1)
  !$OMP THREADPRIVATE(frac_carb)
CONTAINS

!!
!================================================================================================================================
!! SUBROUTINE     : stomate_soil_carbon_discretization_deep_somcycle
!!
!>\BRIEF          Recalculate vegetation cover and LAI
!!
!!\n DESCRIPTION : 
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): None 
!!
!! REFERENCE(S)   : None
!! 
!! FLOWCHART :
!_
!================================================================================================================================

  SUBROUTINE stomate_soil_carbon_discretization_deep_somcycle(kjpindex, index, itau, time_step, lalo, clay, &
       silt, tsurf, tprof, hslong_in, snow, heat_Zimov, pb, &  
       sfluxCH4_deep, sfluxCO2_deep, &
       deepSOM_a, deepSOM_s, deepSOM_p, soil_n_min, O2_soil, CH4_soil, O2_snow, CH4_snow, &
       depth_organic_soil, som_input, veget_max, &
       altmax, altmax_ind, altmax_lastyear, altmax_ind_lastyear, &
       som, som_surf, resp_hetero_soil, fbact, CN_target, fixed_cryoturbation_depth, &
       snowdz,snowrho, n_mineralisation, root_depth, MatrixA, &
       CN_som_litter_longterm, tau_CN_longterm, &
       deepSOM_peat, soiltile)


!! 0. Variable and parameter declaration    

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                                  :: kjpindex
    REAL(r_std), INTENT(in)                                     :: time_step                  !! time step in seconds
    INTEGER(i_std), intent(in)                                  :: itau                       !! time step number
    REAL(r_std),DIMENSION(:,:),INTENT(in)                       :: lalo                       !! Geogr. coordinates (latitude,longitude) (degrees)
    REAL(r_std), DIMENSION(:), INTENT(in)                       :: pb                         !! surface pressure [pa]
    REAL(r_std), DIMENSION(:), INTENT(in)                       :: clay                       !! clay content
    REAL(r_std), DIMENSION(:), INTENT(in)                       :: silt                       !! clay content
    INTEGER(i_std),DIMENSION(:),INTENT(in)                      :: index                      !! Indeces of the points on the map
    REAL(r_std), DIMENSION(:), INTENT (in)                      :: snow                       !! Snow mass [Kg/m^2]
    REAL(r_std), DIMENSION(:,:), INTENT(in)                     :: snowdz                     !! Snow depth [m]
    REAL(r_std), DIMENSION(:,:), INTENT(in)                     :: snowrho                    !! snow density  (Kg/m^3) 
    REAL(r_std), DIMENSION(:,:,:),INTENT (in)                   :: tprof                      !! deep temperature profile
    REAL(r_std), DIMENSION(:,:,:),INTENT (in)                   :: hslong_in                  !! deep long term soil humidity profile
    REAL(r_std), DIMENSION(:,:,:,:),INTENT(in)                  :: som_input                  !! carbon going into carbon pools  [gC/(m**2 of ground)/day]
    REAL(r_std), DIMENSION(:,:),INTENT(in)                      :: veget_max                  !! Maximum vegetation fraction
    REAL(r_std), DIMENSION(:), INTENT(in)                       :: tsurf                      !! skin temperature  [K]
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)                   :: fbact                      !! turnover constant (day)
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)                   :: CN_target                  !! C to N ratio of SOM flux from one pool to another (gN m-2 dt-1) 
    REAL(r_std), DIMENSION (:,:,:), INTENT(in)                  :: root_depth                 !! Node and interface numbers at which the deepest roots
                                                                                              !! occur (1 to nslm, unitless)
    REAL(r_std), DIMENSION(:,:,:),INTENT(in)                    :: soil_n_min                 !! mineral nitrogen in the soil (gN/m**2)
 
   !! 0.2 Output variables
    REAL(r_std), DIMENSION(:), INTENT(out)                      :: sfluxCH4_deep              !! total CH4 flux [g CH4 / m**2 / s]
    REAL(r_std), DIMENSION(:), INTENT(out)                      :: sfluxCO2_deep              !! total CO2 flux [g C / m**2 / s]
    REAL(r_std), DIMENSION(:,:), INTENT(out)                    :: resp_hetero_soil           !! soil heterotrophic respiration (first in gC/day/m**2 of ground )
    REAL(r_std), DIMENSION(:,:,:), INTENT (out)                 :: heat_Zimov                 !! Heating associated with decomposition  [W/m**3 soil]
    REAL(r_std), DIMENSION(:,:,:,:), INTENT (out)               :: som                        !! vertically-integrated (diagnostic) soil carbon pool: active, slow, 
                                                                                              !! or passive, (gC/(m**2 of ground))
    REAL(r_std), DIMENSION(:,:,:,:), INTENT (out)               :: som_surf                   !! vertically-integrated (diagnostic) soil carbon pool: active, slow, 
                                                                                              !! or passive, (gC/(m**2 of ground))

    !! 0.3 Modified variables
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)              :: deepSOM_a                  !! Active soil carbon (g/m**3)
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)              :: deepSOM_s                  !! Slow soil carbon (g/m**3) 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)              :: deepSOM_p                  !! Passive soil carbon (g/m**3)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)                :: O2_snow                    !! oxygen in the snow (g O2/m**3 air)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)                :: O2_soil                    !! oxygen in the soil (g O2/m**3 air)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)                :: CH4_snow                   !! methane in the snow (g CH4/m**3 air)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)                :: CH4_soil                   !! methane in the soil (g CH4/m**3 air)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)                  :: altmax                     !! active layer thickness (m)
    INTEGER(i_std), DIMENSION(kjpindex,nvm), INTENT(inout)      :: altmax_ind
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(inout)         :: altmax_lastyear
    INTEGER(i_std), DIMENSION(kjpindex,nvm), INTENT(inout)      :: altmax_ind_lastyear
    REAL(r_std), DIMENSION(:,:),INTENT(inout)                   :: fixed_cryoturbation_depth  !! depth to hold cryoturbation to for fixed runs
    REAL(r_std), DIMENSION(:,:), INTENT(inout)                  :: n_mineralisation           !! net nitrogen mineralisation of decomposing SOM
    REAL(r_std), DIMENSION(:),   INTENT (inout)                 :: depth_organic_soil         !! depth to organic soil
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)              :: MatrixA                    !! Matrix containing the fluxes between the carbon pools
                                                                                              !! per sechiba time step 
                                                                                              !! @tex $(gC.m^2.day^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)                :: CN_som_litter_longterm     !! Longterm CN ratio of litter and som pools (gC/gN)
    REAL(r_std), INTENT(inout)                                  :: tau_CN_longterm            !! Counter used for calculating the longterm CN ratio 
                                                                                              !! of SOM and litter pools (seconds)

    ! peatland
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm, nelements)                  :: deepSOM_pt      !! peat carbon concentration (g/m3)
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm, nelements), INTENT(inout)   :: deepSOM_peat    !! soil carbon for peat (g/m2)
    REAL(r_std), DIMENSION(kjpindex,nvm, nelements)                        :: peat_OLT        !! peat organic layer thickness (m)
    REAL(r_std),DIMENSION (:,:), INTENT(in)                                :: soiltile        !! Fraction of soil tile within vegtot (0-1, unitless)


    !! 0.4 Local variables
    REAL(r_std)                                                    :: dt                      !! Time step \f$(dt_sechiba one_day^{-1})$\f
    REAL(r_std), DIMENSION(kjpindex)                               :: overburden
    REAL(r_std), DIMENSION(kjpindex,nvm)                           :: fluxCH4,febul
    REAL(r_std), DIMENSION(kjpindex,nvm)                           :: sfluxCH4    
    REAL(r_std), DIMENSION(kjpindex,nvm)                           :: flupmt
    REAL(r_std), DIMENSION(kjpindex,nvm)                           :: MT                              !! depth-integrated methane consumed in methanotrophy
    REAL(r_std), DIMENSION(kjpindex,nvm)                           :: MG                              !! depth-integrated methane released in methanogenesis
    REAL(r_std), DIMENSION(kjpindex,nvm)                           :: CH4i                            !! depth-integrated methane
    REAL(r_std), DIMENSION(kjpindex,nvm)                           :: CH4ii                           !! depth-integrated initial methane
    REAL(r_std), DIMENSION(kjpindex,nvm)                           :: dC1i                            !! depth-integrated oxic decomposition carbon

    REAL(r_std), DIMENSION(kjpindex,nvm)                           :: Tref                            !! Ref. temperature for growing season caluculation (C)	
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm)                     :: deltaCH4g, deltaCH4
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm,nelements)           :: deltaSOM1_a,deltaSOM1_s,deltaSOM1_p,deltaSOM2,deltaSOM3
    REAL(r_std), DIMENSION(kjpindex,ngrnd)                         :: deltaSOM1_a_pftmean,deltaSOM1_s_pftmean,deltaSOM1_p_pftmean
    REAL(r_std), DIMENSION(kjpindex)                               :: deltaSOM1_a_pftmean_stock,deltaSOM1_s_pftmean_stock,deltaSOM1_p_pftmean_stock
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm)                     :: CH4ini_soil
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm)                     :: hslong                          !! deep long term soil humidity profile
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm,nelements)           :: fluxout_a                       !! fluxes leaving the active pool (gC_N m-2 dt-1) 
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm,nelements)           :: fluxout_s                       !! fluxes leaving the slow pool (gC_N m-2 dt-1)   
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm,nelements)           :: fluxout_p                       !! fluxes leaving the passive pool (gC_N m-2 dt-1) 
    REAL(r_std), DIMENSION(kjpindex,nvm)                           :: fluxout_ave_a                   !! weighted average of thefluxes leaving the active pool (gC_N m-2 dt-1) 
    REAL(r_std), DIMENSION(kjpindex,nvm)                           :: fluxout_ave_s                   !! weighted average of thefluxes leaving the slow pool (gC_N m-2 dt-1)   
    REAL(r_std), DIMENSION(kjpindex,nvm)                           :: fluxout_ave_p                   !! weighted average of thefluxes leaving the passive pool (gC_N m-2 dt-1)  
    REAL(r_std), DIMENSION(ncarb,ncarb,ngrnd,nelements)            :: somflux                         !! fluxes between soil orgnanic matter reservoirs
    INTEGER(i_std)                   :: ip, il, itz, iz
    REAL(r_std), SAVE, DIMENSION(3)  :: lhc                       !! specific heat of soil organic matter oxidation (J/kg carbon)
!$OMP THREADPRIVATE(lhc)
    REAL(r_std), SAVE                :: O2m                       !! oxygen concentration [g/m3] below which there is anoxy 
!$OMP THREADPRIVATE(O2m)
    LOGICAL, SAVE                    :: ok_methane                !! Is Methanogenesis and -trophy taken into account?
!$OMP THREADPRIVATE(ok_methane)
    LOGICAL, SAVE                    :: ok_cryoturb               !! cryoturbate the carbon?
!$OMP THREADPRIVATE(ok_cryoturb)
    REAL(r_std), SAVE                :: cryoturbation_diff_k_in   !! input time constant of cryoturbation (m^2/y)
!$OMP THREADPRIVATE(cryoturbation_diff_k_in)
    REAL(r_std), SAVE                :: bioturbation_diff_k_in    !! input time constant of bioturbation (m^2/y)
!$OMP THREADPRIVATE(bioturbation_diff_k_in)
    REAL(r_std), SAVE                :: tau_CH4troph              !! time constant of methanetrophy (s)
!$OMP THREADPRIVATE(tau_CH4troph)
    REAL(r_std), SAVE                :: fbactratio                !! time constant of methanogenesis (ratio to that of oxic)
!$OMP THREADPRIVATE(fbactratio)
    LOGICAL, SAVE                    :: firstcall = .TRUE.        !! first call?
!$OMP THREADPRIVATE(firstcall)
    REAL(r_std), SAVE, DIMENSION(2)  :: lhCH4                     !! specific heat of methane transformation  (J/kg) (/ 3.1E6, 9.4E6 /)
!$OMP THREADPRIVATE(lhCH4)
    LOGICAL, SAVE                    :: oxlim                     !! O2 limitation taken into account
!$OMP THREADPRIVATE(oxlim)
    REAL(r_std), PARAMETER           :: refdep = 0.20_r_std       !! Depth to compute reference temperature for the growing season (m). WH2000 use 0.50
    REAL(r_std), PARAMETER           :: Tgr  = 5.                 !! Temperature when plant growing starts and this becomes constant
    REAL(r_std)                      :: scnd
    REAL(r_std)                      :: organic_layer_thickness
    INTEGER(i_std)                   :: ier, iv, m, jv, iele
    CHARACTER(80)                    :: yedoma_map_filename
    REAL(r_std)                      :: yedoma_depth, yedoma_cinit_act, yedoma_cinit_slo, yedoma_cinit_pas
    LOGICAL                          :: reset_yedoma_carbon
    LOGICAL, SAVE                    :: MG_useallCpools = .true.  !! Do we allow all three C pools to feed methanogenesis?
!$OMP THREADPRIVATE(MG_useallCpools)
    CHARACTER(LEN=10)                :: part_str                  !! string suffix indicating an index
    REAL(r_std), SAVE                :: max_shum_value = 1.0      !! maximum saturation degree on the thermal axes
!$OMP THREADPRIVATE(max_shum_value)
    REAL(r_std)                      :: delta_n_min               !! Amount of nitrogen that was overspend 
                                                                  !! @tex $(N gN m^{-2})$ @endtex
    REAL(r_std)                      :: reduction_factor          !! Factor to reduce the som pools in case of 
                                                                  !! overspending (unitless)
    REAL(r_std), DIMENSION(kjpindex,nvm)                  :: net_som_n_flux   !! Total flux between all soil pools. Used 
                                                                              !! to adjust som if overspending 
                                                                              !! @tex $(N gN m^{-2})$ @endtex

    REAL(r_std), DIMENSION (kjpindex,nvm)                         :: veget_max_bg
    INTEGER(i_std)                                        :: imbc, igrnd,ivm  !! Indices (unitless)
    INTEGER(i_std)                                        :: inbpools, icarb  !! Indices (unitless)
    REAL(r_std), DIMENSION(kjpindex,nvm,nmbcomp,nelements):: check_intern     !! Contains the components of the internal
                                                                              !! mass balance check for this routine
                                                                              !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(kjpindex,nvm,nelements)        :: closure_intern   !! Check closure of internal mass balance
                                                                              !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(kjpindex,nvm,nelements)        :: pool_start       !! Start and end pool of this routine 
                                                                              !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(kjpindex,nvm,nelements)        :: pool_end         !! Start and end pool of this routine 
                                                                              !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex

!_ ================================================================================================================================

    IF (printlev>=3) WRITE(*,*) 'Entering  stomate_soil_carbon_discretization_deep_somcycle'
    
    !! 0. first call 
    IF ( firstcall ) THEN                

       overburden(:)=1.
       !
       !Config Key   = organic_layer_thickness
       !Config Desc  = The thickness of organic layer
       !Config Def   = 0.0 
       !Config If    = OK_SOIL_CARBON_DISCRETIZATION
       !Config Help  = This parameters allows the user to prescibe the organic
       !Config         layer thickness 
       !Config Units = [-]
       !
       organic_layer_thickness = 0.
       CALL getin_p('organic_layer_thickness', organic_layer_thickness)
       depth_organic_soil(:) = overburden(:)*organic_layer_thickness
       
       !Config Key   = OK_METHANE
       !Config Desc  = Is Methanogenesis and methanotrophy taken into account?
       !Config Def   = n
       !Config If    = OK_SOIL_CARBON_DISCRETIZATION
       !Config Help  = 
       !Config Units = [FLAG]
       !
       ok_methane = .FALSE.
       CALL getin_p('OK_METHANE',ok_methane)
       
       IF (ok_methane) THEN
          CALL ipslerr_p(3,'stomate_soil_carbon_discretization',&
               'there is some code to calculate CH4 but mass conservation', &
               'has not yet been checked.',&
               'confirm mass conservation before using this code!')
       END IF
       !
       !Config Key   = HEAT_CO2_ACT
       !Config Desc  = specific heat of soil organic matter oxidation for active carbon
       !Config Def   = 40.0E6 
       !Config If    = OK_SOIL_CARBON_DISCRETIZATION
       !Config Help  = 
       !Config Units = [J/Kg C]
       !
       lhc(iactive) = 40.0e6
       CALL getin_p('HEAT_CO2_ACT',lhc(iactive))
       !
       !Config Key   = HEAT_CO2_SLO
       !Config Desc  = specific heat of soil organic matter oxidation for slow carbon pool 
       !Config Def   = 30.0E6 
       !Config If    = OK_SOIL_CARBON_DISCRETIZATION
       !Config Help  = 
        !Config Units = [J/Kg C]
       !
       lhc(islow) = 30.0E6
       CALL getin_p('HEAT_CO2_SLO',lhc(islow))
       !
       !Config Key   = HEAT_CO2_PAS
       !Config Desc  = specific heat of soil organic matter oxidation for passive carbon pool 
       !Config Def   = 10.0E6 
       !Config If    = OK_SOIL_CARBON_DISCRETIZATION
       !Config Help  = 
       !Config Units = [J/Kg C]
       !
       lhc(ipassive) = 10.0e6
       CALL getin_p('HEAT_CO2_PAS',lhc(ipassive))
       !
       !Config Key   = TAU_CH4_TROPH
       !Config Desc  = time constant of methanetrophy
       !Config Def   = 432000 
       !Config If    = OK_SOIL_CARBON_DISCRETIZATION
       !Config Help  = 
       !Config         
       !Config Units = [s]
       ! 
       tau_CH4troph = 432000
       CALL getin_p('TAU_CH4_TROPH',tau_CH4troph)
       !
       !Config Key   = TAU_CH4_GEN_RATIO
       !Config Desc  = time constant of methanogenesis (ratio to that of oxic)
       !Config Def   = 9.0 
       !Config If    = OK_SOIL_CARBON_DISCRETIZATION
       !Config Help  = 
       !Config Units = [-]
       !  
       fbactratio = 9.0
       CALL getin_p('TAU_CH4_GEN_RATIO',fbactratio)
       !
       !Config Key   = O2_SEUIL_MGEN
       !Config Desc  = oxygen concentration below which there is anoxy 
       !Config Def   = 3.0 
       !Config If    = OK_SOIL_CARBON_DISCRETIZATION
       !Config Help  = 
       !Config Units = [g/m3]
       !  
       O2m = 3.0
       CALL getin_p('O2_SEUIL_MGEN',O2m)
       !
       !Config Key   = HEAT_CH4_GEN
       !Config Desc  = specific heat of methanogenesis 
       !Config Def   = 0 
       !Config If    = OK_SOIL_CARBON_DISCRETIZATION
       !Config Help  = 
       !Config Units = [J/kgC]
       !  
       lhCH4(1) = 0
       CALL getin_p('HEAT_CH4_GEN',lhCH4(1))
       !
       !Config Key   = HEAT_CH4_TROPH
       !Config Desc  = specific heat of methanotrophy 
       !Config Def   = 0 
       !Config If    = OK_SOIL_CARBON_DISCRETIZATION
       !Config Help  = 
       !Config         
       !Config Units = [J/kgC]
       !  
       lhCH4(2) = 0
       CALL getin_p('HEAT_CH4_TROPH',lhCH4(2))
       !
       !Config Key   = O2_LIMIT
       !Config Desc  = O2 limitation taken into account
       !Config Def   = n
       !Config If    = OK_SOIL_CARBON_DISCRETIZATION
       !Config Help  = 
       !Config Units = [flag]
       ! 
       oxlim=.FALSE. 
       CALL getin_p('O2_LIMIT',oxlim)
       !
       !Config Key   = cryoturbate
       !Config Desc  = Do we allow for cyoturbation?
       !Config Def   = y 
       !Config If    = OK_SOIL_CARBON_DISCRETIZATION
       !Config Help  = 
       !Config Units = [flag]
       ! 
       ok_cryoturb=.TRUE.
       CALL getin_p('cryoturbate',ok_cryoturb)
       !
       !Config Key   = cryoturbation_diff_k_in
       !Config Desc  = diffusion constant for cryoturbation 
       !Config Def   = 0.001 
       !Config If    = OK_SOIL_CARBON_DISCRETIZATION
       !Config Help  = 
       !Config Units = [m2/year]
       !  
       cryoturbation_diff_k_in = .001
       CALL getin_p('cryoturbation_diff_k',cryoturbation_diff_k_in)
       !
       !Config Key   = bioturbation_diff_k_in
       !Config Desc  = diffusion constant for bioturbation 
       !Config Def   = 0.0
       !Config If    = OK_SOIL_CARBON_DISCRETIZATION
       !Config Help  = 
       !Config Units = [m2/year]
       !  
       bioturbation_diff_k_in = 0.0001
       CALL getin_p('bioturbation_diff_k',bioturbation_diff_k_in)
       !
       !Config Key   = MG_useallCpools
       !Config Desc  = Do we allow all three C pools to feed methanogenesis?
       !Config Def   = y 
       !Config If    = OK_SOIL_CARBON_DISCRETIZATION
       !Config Help  = 
       !Config Units = [flag]
       !    
       MG_useallCpools = .TRUE.
       CALL getin_p('MG_useallCpools', MG_useallCpools)
       !
       !Config Key   = max_shum_value
       !Config Desc  = maximum saturation degree on the thermal axes
       !Config Def   = 1 
       !Config If    = OK_SOIL_CARBON_DISCRETIZATION
       !Config Help  = 
       !Config Units = [-]
       !   
       max_shum_value=1.0 
       CALL getin_p('max_shum_value',max_shum_value)
       hslong(:,:,:) = MAX(MIN(hslong_in(:,:,:),max_shum_value),zero)
       !

       !!  Arrays allocations

       ALLOCATE (veget_mask_2d(kjpindex,nvm),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'stomate_soil_carbon_discretization_deep_somcycle', &
            'Pb in alloc for veget_mask_2d','','')
  
       ALLOCATE (alt(kjpindex,nvm),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'stomate_soil_carbon_discretization_deep_somcycle', &
            'Pb in alloc for alt','','')
       
       ALLOCATE (alt_ind(kjpindex,nvm),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'stomate_soil_carbon_discretization_deep_somcycle', &
            'Pb in alloc for alt_ind','','')

       ALLOCATE (z_root(kjpindex,nvm),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'stomate_soil_carbon_discretization_deep_somcycle', &
            'Pb in alloc for z_root','','')

       ALLOCATE (rootlev(kjpindex,nvm),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'stomate_soil_carbon_discretization_deep_somcycle', &
            'Pb in alloc for rootlev','','')

       ALLOCATE (heights_snow(kjpindex,nvm),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'stomate_soil_carbon_discretization_deep_somcycle', &
            'Pb in alloc for heights_snow','','')

       ALLOCATE (zf_soil(0:ngrnd),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'stomate_soil_carbon_discretization_deep_somcycle', &
            'Pb in alloc for zf_soil','','')
       
       ALLOCATE (zi_soil(ngrnd),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'stomate_soil_carbon_discretization_deep_somcycle', &
            'Pb in alloc for zi_soil','','')

       ALLOCATE (zf_snow(kjpindex,0:nsnow,nvm),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'stomate_soil_carbon_discretization_deep_somcycle', &
            'Pb in alloc for zf_snow','','')

       ALLOCATE (zi_snow(kjpindex,nsnow,nvm),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'stomate_soil_carbon_discretization_deep_somcycle', &
            'Pb in alloc for zi_snow','','')

       ALLOCATE (zf_snow_nopftdim(kjpindex,0:nsnow),stat=ier)   
       IF (ier /= 0) CALL ipslerr_p(3,'stomate_soil_carbon_discretization_deep_somcycle', &
            'Pb in alloc for zf_snow_nopftdim','','')
       
       ALLOCATE (zi_snow_nopftdim(kjpindex,nsnow),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'stomate_soil_carbon_discretization_deep_somcycle', &
            'Pb in alloc for zi_snow_nopftdim','','')

       ALLOCATE (airvol_soil(kjpindex,ngrnd,nvm),stat=ier)       
       IF (ier /= 0) CALL ipslerr_p(3,'stomate_soil_carbon_discretization_deep_somcycle', &
            'Pb in alloc for airvol_soil','','')
       
       ALLOCATE (totporO2_soil(kjpindex,ngrnd,nvm),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'stomate_soil_carbon_discretization_deep_somcycle', &
            'Pb in alloc for totporO2_soil','','')

       ALLOCATE (totporCH4_soil(kjpindex,ngrnd,nvm),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'stomate_soil_carbon_discretization_deep_somcycle', &
            'Pb in alloc for totporCH4_soil','','')

       ALLOCATE (diffO2_soil(kjpindex,ngrnd,nvm),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'stomate_soil_carbon_discretization_deep_somcycle', &
            'Pb in alloc for diffO2_soil','','')

       ALLOCATE (diffCH4_soil(kjpindex,ngrnd,nvm),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'stomate_soil_carbon_discretization_deep_somcycle', &
            'Pb in alloc for diffCH4_soil','','')

       ALLOCATE (airvol_snow(kjpindex,nsnow,nvm),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'stomate_soil_carbon_discretization_deep_somcycle', &
            'Pb in alloc for airvol_snow','','')

       ALLOCATE (totporO2_snow(kjpindex,nsnow,nvm),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'stomate_soil_carbon_discretization_deep_somcycle', &
            'Pb in alloc for totporO2_snow','','')

       ALLOCATE (totporCH4_snow(kjpindex,nsnow,nvm),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'stomate_soil_carbon_discretization_deep_somcycle', &
            'Pb in alloc for totporCH4_snow','','')

       ALLOCATE (conduct_snow(kjpindex,nsnow,nvm),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'stomate_soil_carbon_discretization_deep_somcycle', &
            'Pb in alloc for conduct_snow','','')

       ALLOCATE (diffO2_snow(kjpindex,nsnow,nvm),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'stomate_soil_carbon_discretization_deep_somcycle', &
            'Pb in alloc for diffO2_snow','','')

       ALLOCATE (diffCH4_snow(kjpindex,nsnow,nvm),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'stomate_soil_carbon_discretization_deep_somcycle', &
            'Pb in alloc for diffCH4_snow','','')

       ALLOCATE (deep_rhC_a_pftmean(kjpindex,ngrnd),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'stomate_soil_carbon_discretization_deep_somcycle', &
            'Pb in alloc for deep_rhC_a_pftmean','','')

       ALLOCATE (deep_rhC_s_pftmean(kjpindex,ngrnd),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'stomate_soil_carbon_discretization_deep_somcycle', &
            'Pb in alloc for deep_rhC_s_pftmean','','')

       ALLOCATE (deep_rhC_p_pftmean(kjpindex,ngrnd),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'stomate_soil_carbon_discretization_deep_somcycle', &
            'Pb in alloc for deep_rhC_p_pftmean','','')

       ALLOCATE (rhC_a_pftmean(kjpindex),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'stomate_soil_carbon_discretization_deep_somcycle', &
            'Pb in alloc for rhC_a_pftmean','','')

       ALLOCATE (rhC_s_pftmean(kjpindex),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'stomate_soil_carbon_discretization_deep_somcycle', &
            'Pb in alloc for rhC_s_pftmean','','')

       ALLOCATE (rhC_p_pftmean(kjpindex),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'stomate_soil_carbon_discretization_deep_somcycle', &
            'Pb in alloc for rhC_p_pftmean','','')

       veget_max_bg(:,2:nvm) = veget_max(:,2:nvm)
       veget_max_bg(:,1) = MAX((un - SUM(veget_max(:,2:nvm), 2)), zero)
!!       veget_mask_2d(:,:) = veget_max_bg .GT. EPSILON(zero)
!!       WHERE( ALL((.NOT. veget_mask_2d(:,:)), dim=2) )
!!          veget_mask_2d(:,1) = .TRUE.
!!       END WHERE
       veget_mask_2d(:,:) = .TRUE.

       alt(:,:) = 0
       alt_ind(:,:) = 0
       z_root(:,:) = 0
       rootlev(:,:) = 0
       febul(:,:) = 0
       flupmt(:,:) = 0
       ! make sure gas concentrations where not defined by veget_mask are equal
       !to initial conditions
       IF (ok_methane) THEN
          DO iv = 1, ngrnd
             WHERE ( .NOT. veget_mask_2d(:,:) )
                O2_soil(:,iv,:) = O2_init_conc
                CH4_soil(:,iv,:) = CH4_init_conc
             END WHERE
          END DO
          DO iv = 1, nsnow
             WHERE ( .NOT. veget_mask_2d(:,:) )
                O2_snow(:,iv,:) = O2_surf
                CH4_snow(:,iv,:) = CH4_surf
             END WHERE
          END DO
       ELSE
          DO iv = 1, ngrnd
             WHERE ( .NOT. veget_mask_2d(:,:) )
                O2_soil(:,iv,:) = zero
                CH4_soil(:,iv,:) = zero
             END WHERE
          END DO
          DO iv = 1, nsnow
             WHERE ( .NOT. veget_mask_2d(:,:) )
                O2_snow(:,iv,:) = zero
                CH4_snow(:,iv,:) = zero
             END WHERE
          END DO
       ENDIF
       

       heights_snow(:,:) = zero
       zf_soil(:) = zero
       zi_soil(:) = zero
       zf_snow(:,:,:) = zero
       zi_snow(:,:,:) = zero
       zf_snow_nopftdim(:,:) = zero
       zi_snow_nopftdim(:,:) = zero
       airvol_soil(:,:,:) = zero
       totporO2_soil(:,:,:) = zero
       totporCH4_soil(:,:,:) = zero
       diffO2_soil(:,:,:) = zero
       diffCH4_soil(:,:,:) = zero
       airvol_snow(:,:,:) = zero
       totporO2_snow(:,:,:) = zero
       totporCH4_snow(:,:,:) = zero
       conduct_snow(:,:,:) = zero
       diffO2_snow(:,:,:) = zero
       diffCH4_snow(:,:,:) = zero

       ! get snow and soil levels
       DO iv = 1, nvm
          heights_snow(:,iv) = SUM(snowdz(:,1:nsnow), 2)
       ENDDO
       ! Calculating intermediate and full depths for snow
       call snowlevels (kjpindex, snowdz, zi_snow, zf_snow, veget_max_bg)

       ! here we need to put the shallow and deep soil levels together to make the complete soil levels. 
       ! This requires pulling in the indices from thermosoil and deepsoil_freeze.
       zi_soil(:) = znt(:)
       zf_soil(1:ngrnd) = zlt(:)
       zf_soil(0) = 0.

       !    allocate arrays for gas diffusion        !
       !    get diffusion coefficients: heat capacity,
       !    conductivity, and oxygen diffusivity
       
       CALL get_gasdiff (kjpindex,hslong,tprof,snow,airvol_snow, &
            totporO2_snow,totporCH4_snow,diffO2_snow,diffCH4_snow, &
            airvol_soil,totporO2_soil,totporCH4_soil,diffO2_soil,diffCH4_soil, snowrho)

       !
       !    initialize soil temperature calculation
       !
       CALL soil_gasdiff_main (kjpindex,time_step,index,'initialize', &	
            pb,tsurf,tprof,diffO2_snow,diffCH4_snow, &
            totporO2_snow,totporCH4_snow,O2_snow,CH4_snow,diffO2_soil,diffCH4_soil, &
            totporO2_soil,totporCH4_soil,O2_soil,CH4_soil, zi_snow, zf_snow)
       
       IF (.NOT.ok_methane) THEN
          totporCH4_snow(:,:,:) = zero
          CH4_snow(:,:,:) = zero
          totporCH4_soil(:,:,:) = zero
          CH4_soil(:,:,:) = zero
       END IF

       !
       !    calculate the coefficients
       !
       CALL soil_gasdiff_main (kjpindex,time_step,index,'coefficients', &
            pb,tsurf,tprof,diffO2_snow,diffCH4_snow, &
            totporO2_snow,totporCH4_snow,O2_snow,CH4_snow,diffO2_soil,diffCH4_soil, &
            totporO2_soil,totporCH4_soil,O2_soil,CH4_soil, zi_snow, zf_snow)
       

       IF (printlev>=3 ) THEN
          WRITE(*,*) 'stomate_soil_carbon_discretization_deep_somcycle: finished firstcall calcs'
       ENDIF

       ! reset
       !
       !Config Key   = reset_yedoma_carbon
       !Config Desc  = Do we reset carbon concentrations for yedoma region?
       !Config Def   = n 
       !Config If    = OK_SOIL_CARBON_DISCRETIZATION
       !Config Help  = 
       !Config Units = [flag]
       ! 
       reset_yedoma_carbon = .false.
       CALL getin_p('reset_yedoma_carbon',reset_yedoma_carbon)

       IF (reset_yedoma_carbon) THEN
          !
          !Config Key   = yedoma_map_filename
          !Config Desc  = The filename for yedoma map
          !Config Def   = yedoma_map.nc 
          !Config If    = OK_SOIL_CARBON_DISCRETIZATION
          !Config Help  = 
          !Config Units = []
          yedoma_map_filename = 'NONE'
          CALL getin_p('yedoma_map_filename', yedoma_map_filename)
          !
          !Config Key   = yedoma_depth
          !Config Desc  = The depth for soil carbon in yedoma
          !Config Def   = 20 
          !Config If    = OK_SOIL_CARBON_DISCRETIZATION
          !Config Help  = 
          !Config Units = [m]
          yedoma_depth = zero 
          CALL getin_p('yedoma_depth', yedoma_depth)
          !
          !Config Key   = deepC_a_init
          !Config Desc  = Carbon concentration for active soil C pool in yedoma
          !Config Def   = 1790.1  
          !Config If    = OK_SOIL_CARBON_DISCRETIZATION
          !Config Help  = 
          !Config Units = [gC/?]
          yedoma_cinit_act = zero 
          CALL getin_p('deepC_a_init', yedoma_cinit_act)
          !
          !Config Key   = deepC_s_init
          !Config Desc  = Carbon concentration for slow soil C pool in yedoma
          !Config Def   = 14360.8
          !Config If    = OK_SOIL_CARBON_DISCRETIZATION
          !Config Help  = 
          !Config Units = [gC/?]
          yedoma_cinit_slo = zero
          CALL getin_p('deepC_s_init', yedoma_cinit_slo)
          !
          !Config Key   = deepC_p_init
          !Config Desc  = Carbon concentration for passive soil C pool in yedoma
          !Config Def   = 1436
          !Config If    = OK_SOIL_CARBON_DISCRETIZATION
          !Config Help  = 
          !Config Units = [gC/>]
          yedoma_cinit_pas = zero
          CALL getin_p('deepC_p_init', yedoma_cinit_pas)
          ! intialize the yedoma carbon stocks
          CALL initialize_yedoma_carbonstocks(kjpindex, lalo, deepSOM_a, deepSOM_s, deepSOM_p, &
               yedoma_map_filename, yedoma_depth, yedoma_cinit_act,yedoma_cinit_slo, yedoma_cinit_pas, altmax_ind)
       ENDIF


    ENDIF ! firstcall

    !! Initialize check for mass balance closure
    !  The mass balance is calculated at the end of this routine.
    IF (err_act.GT.1) THEN
       dt = dt_sechiba/one_day
       pool_start(:,:,:) = zero

       DO iele = 1,nelements
          DO igrnd = 1,ngrnd
             pool_start(:,:,iele) = pool_start(:,:,iele) + &
                  (deepSOM_a(:,igrnd,:,iele) + deepSOM_s(:,igrnd,:,iele) + deepSOM_p(:,igrnd,:,iele)) * &
                  (zf_soil(igrnd)-zf_soil(igrnd-1)) * veget_max(:,:)
          END DO
       ENDDO

       DO iele = 1,nelements
          DO icarb = 1,ncarb
             pool_start(:,:,iele) = pool_start(:,:,iele) + &
                  som_input(:,icarb,:,iele) * dt * veget_max(:,:)
          ENDDO
       ENDDO
       
       ! Initial nitrogen pool not included in the loops above
       pool_start(:,:,initrogen) = pool_start(:,:,initrogen) + &
            ( n_mineralisation(:,:) ) * veget_max(:,:)

    END IF ! err_act.GT.1


    !! Prepare values for arrays
    veget_max_bg(:,2:nvm) = veget_max(:,2:nvm)
    veget_max_bg(:,1) = MAX((un - SUM(veget_max(:,2:nvm), 2)), zero)

    IF ( ANY(rootlev(:,:) .GT. ngrnd) ) THEN
       WRITE(*,*) 'problems with rootlev:', rootlev
       CALL ipslerr_p(3,'stomate_soil_carbon_discretization_deep_somcycle',&
            'Problems with rootlev','','')
    ENDIF

    DO iv = 1, nvm
       heights_snow(:,iv) = SUM(snowdz(:,1:nsnow), 2)
    ENDDO
    !
    ! define initial CH4 value (before the time step)
    CH4ini_soil(:,:,:) = CH4_soil(:,:,:)

    ! apply maximum soil wetness criteria to prevent soils from turning to wetlands where they aren't supposed to
    hslong(:,:,:) = MAX(MIN(hslong_in(:,:,:),max_shum_value),zero)


    ! update the gas profiles
    !
    CALL soil_gasdiff_main (kjpindex, time_step, index, 'diffuse', &
         pb,tsurf,tprof,diffO2_snow,diffCH4_snow, &
         totporO2_snow,totporCH4_snow,O2_snow,CH4_snow,diffO2_soil,diffCH4_soil, &
         totporO2_soil,totporCH4_soil,O2_soil,CH4_soil, zi_snow, zf_snow)

    IF (.NOT.ok_methane) THEN
       totporCH4_snow(:,:,:) = zero
       CH4_snow(:,:,:) = zero
       totporCH4_soil(:,:,:) = zero
       CH4_soil(:,:,:) = zero
    END IF

    ! get new snow levels and interpolate gases on these levels
    !
    CALL snow_interpol (kjpindex,O2_snow, CH4_snow, zi_snow, zf_snow, veget_max_bg, snowdz)

    CALL altcalc (kjpindex, time_step, scnd, tprof, zi_soil, alt, alt_ind, altmax, altmax_ind, &
         altmax_lastyear, altmax_ind_lastyear)     

    ! Rooting depth is calculated in hydrol_root_profile. In hydrol.f90
    ! altmax has not been updated yet and thus has the same value as
    ! altmax_lastyear in stomate_soil_carbon_discretization.
    ! In hydrol_root_profile it is checked that the rooting depth is within 
    ! the active layer 
    WHERE ( altmax_lastyear(:,:) .LT. z_root_max .and. veget_mask_2d(:,:) )
       z_root(:,:) = altmax_lastyear(:,:)
       rootlev(:,:) = altmax_ind_lastyear(:,:) 
    ELSEWHERE ( veget_mask_2d(:,:) )
       z_root(:,:) = z_root_max
       rootlev(:,:) = altmax_ind_lastyear(:,:)  
    ENDWHERE

    !
    ! Carbon input into the soil
    !
    !+++CHECK+++
    ! The name z_root is confusing. z_root describes the depth of the
    ! som_input rather than the actual depth of the roots (which is done
    ! by root_depth).
    CALL sominput(kjpindex, time_step, itau*time_step, tprof, tsurf, hslong, &
         z_root, altmax_lastyear, deepSOM_a, deepSOM_s, deepSOM_p, &
         som_input, veget_max_bg, root_depth)
    !+++++++++++

    CALL permafrost_decomp (kjpindex, time_step, tprof, airvol_soil, &
         oxlim, tau_CH4troph, ok_methane, fbactratio, O2m, &
         totporO2_soil, totporCH4_soil, hslong, clay, silt, &
         deepSOM_a, deepSOM_s, deepSOM_p, &
         deltaCH4g, deltaCH4, deltaSOM1_a, deltaSOM1_s, deltaSOM1_p, deltaSOM2, &
         deltaSOM3, deltaSOM1_a_pftmean, deltaSOM1_s_pftmean, deltaSOM1_p_pftmean, &
         deltaSOM1_a_pftmean_stock, deltaSOM1_s_pftmean_stock, deltaSOM1_p_pftmean_stock, &
         O2_soil, CH4_soil, fbact, CN_target, MG_useallCpools, veget_max, veget_max_bg, &
         fluxout_a,fluxout_s,fluxout_p,somflux, &
         deepSOM_pt,  deepSOM_peat,  peat_OLT,  soiltile)
    
    
    DO ip = 1, kjpindex
       DO iv = 1, nvm
          IF ( veget_mask_2d(ip,iv) ) THEN
             ! oxic decomposition
             heat_Zimov(ip,:,iv) = lhc(iactive)*1.E-3*deltaSOM1_a(ip,:,iv,icarbon) + &
                  lhc(islow)*1.E-3*deltaSOM1_s(ip,:,iv,icarbon) + &
                  lhc(ipassive)*1.E-3*deltaSOM1_p(ip,:,iv,icarbon)
             !
             ! methanogenesis
             heat_Zimov(ip,:,iv) = heat_Zimov(ip,:,iv) + lhCH4(1)*1.E-3*deltaSOM2(ip,:,iv,icarbon)
             !
             ! methanotrophy
             heat_Zimov(ip,:,iv) = heat_Zimov(ip,:,iv) + lhCH4(2)*1.E-3*deltaCH4(ip,:,iv) *  &
                  totporCH4_soil(ip,:,iv)
             !
             heat_Zimov(ip,:,iv) = heat_Zimov(ip,:,iv)/time_step

             !
             fluxCH4(ip,iv) = zero
          ELSE
             heat_Zimov(ip,:,iv) = zero
             fluxCH4(ip,iv) = zero
          ENDIF
       ENDDO
    ENDDO

    IF  ( .NOT. firstcall) THEN 
       !
       ! Plant-mediated CH4 transport     
       !
       CALL traMplan(CH4_soil,O2_soil,kjpindex,time_step,totporCH4_soil,totporO2_soil,z_root, &
            rootlev,Tgr,Tref,hslong,flupmt, &
            refdep, zi_soil, tprof)
       !	flupmt=zero
       !
       ! CH4 ebullition
       !

       CALL ebullition (kjpindex,time_step,tprof,totporCH4_soil,hslong,CH4_soil,febul)
       !
    ENDIF

    !
    MT(:,:)=zero   
    MG(:,:)=zero    
    CH4i(:,:)=zero  
    CH4ii(:,:)=zero 
    dC1i(:,:)=zero 
    net_som_n_flux(:,:) = zero
    !
    ! Processes contributing to the CH4 fluxes and pools
    IF (ok_methane) THEN
       DO ip = 1, kjpindex
          DO iv = 1, nvm
             IF (  veget_mask_2d(ip,iv) ) THEN
                DO il=1,ngrnd
                   MG(ip,iv) = MG(ip,iv) + deltaCH4g(ip,il,iv)*totporCH4_soil(ip,il,iv) * &
                        ( zf_soil(il) - zf_soil(il-1) )
                   CH4i(ip,iv) = CH4i(ip,iv) + CH4_soil(ip,il,iv)*totporCH4_soil(ip,il,iv) * &
                        (zf_soil(il)-zf_soil(il-1))
                   CH4ii(ip,iv) = CH4ii(ip,iv) +  &
                        CH4ini_soil(ip,il,iv)*totporCH4_soil(ip,il,iv) * &
                        (zf_soil(il)-zf_soil(il-1))
                END DO
             END IF
          END DO
       END DO
    END IF

    ! Processes contributing to the CO2 fluxes and pools
    IF ( test_adjust_rhetero) THEN
       DO ip = 1, kjpindex
          DO iv = 1, nvm
             IF (  veget_mask_2d(ip,iv) ) THEN
                ! There is enough nitrogen in the mineralised pool
                ! to allow for the immobiolization (if any) in 
                ! soil_carbon [gN m-2]. This is the same as what is
                ! done in the trunk.
                DO il=1,ngrnd 
                   net_som_n_flux(ip,iv) = net_som_n_flux(ip,iv) + &
                        (deltaSOM1_a(ip,il,iv,initrogen) + &
                        deltaSOM1_s(ip,il,iv,initrogen) + &
                        deltaSOM1_p(ip,il,iv,initrogen)) * &
                        ( zf_soil(il) - zf_soil(il-1) )
                ENDDO
                
                IF ((soil_n_min(ip,iv,initrate)+soil_n_min(ip,iv,iammonium)) + n_mineralisation(ip,iv) + &
                     net_som_n_flux(ip,iv) .GE. zero) THEN
                   DO il=1,ngrnd        
                      dC1i(ip,iv) = dC1i(ip,iv) + &
                           (deltaSOM1_a(ip,il,iv,icarbon)+deltaSOM1_s(ip,il,iv,icarbon)+deltaSOM1_p(ip,il,iv,icarbon)) * &
                           ( zf_soil(il) - zf_soil(il-1) )
                      n_mineralisation(ip,iv)=n_mineralisation(ip,iv)+ &
                           (deltaSOM1_a(ip,il,iv,initrogen)+deltaSOM1_s(ip,il,iv,initrogen)+deltaSOM1_p(ip,il,iv,initrogen))* &
                           ( zf_soil(il) - zf_soil(il-1) )
                   END DO
                ELSE
                   ! There is not enough soil nitrogen to immobilize later on. 
                   ! Need to decrease the som decomposition. Calculate the 
                   ! amount of nitrogen [gN m-2] that was overspend [gN m-2]
                   delta_n_min =  net_som_n_flux(ip,iv) + soil_n_min(ip,iv,initrate) + soil_n_min(ip,iv,iammonium) + n_mineralisation(ip,iv)
                   
                   ! We calculate a reduction factor for decomposition accordingly 
                   reduction_factor = (net_som_n_flux(ip,iv)-delta_n_min)/net_som_n_flux(ip,iv)
                   IF ((reduction_factor .GT. un) .OR. (reduction_factor .LT. zero)) THEN
                      WRITE (numout,*) ' reduction_factor is lower than zero or greater than one. We stop.'
                      CALL ipslerr_p(3,'stomate_soil_carbon_discretization_deep_somcycle',&
                           'Reduction_factor is lower than zero or greater than one.','','')
                   ELSE
                      DO il=1,ngrnd
                         dC1i(ip,iv) = dC1i(ip,iv) + &
                              reduction_factor * (deltaSOM1_a(ip,il,iv,icarbon)+deltaSOM1_s(ip,il,iv,icarbon)+deltaSOM1_p(ip,il,iv,icarbon)) * &
                              ( zf_soil(il) - zf_soil(il-1) )
                         n_mineralisation(ip,iv)=n_mineralisation(ip,iv)+ &
                              reduction_factor * (deltaSOM1_a(ip,il,iv,initrogen)+deltaSOM1_s(ip,il,iv,initrogen)+deltaSOM1_p(ip,il,iv,initrogen))* &
                              ( zf_soil(il) - zf_soil(il-1) )
                         DO iele=1,nelements
                            IF (iele .EQ. icarbon) THEN
                               deepSOM_a(ip,il,iv,iele) = deepSOM_a(ip,il,iv,iele) + (un-reduction_factor) * (deltaSOM1_a(ip,il,iv,iele)/fr(ip,iactive,iv))
                               deepSOM_s(ip,il,iv,iele) = deepSOM_s(ip,il,iv,iele) + (un-reduction_factor) * (deltaSOM1_s(ip,il,iv,iele)/fr(ip,islow,iv)) 
                               deepSOM_p(ip,il,iv,iele) = deepSOM_p(ip,il,iv,iele) + (un-reduction_factor) * (deltaSOM1_p(ip,il,iv,iele)/fr(ip,ipassive,iv)) 
                            ELSEIF (iele .EQ. initrogen) THEN
                               deepSOM_a(ip,il,iv,iele) = deepSOM_a(ip,il,iv,iele) + (un-reduction_factor) * (deltaSOM1_a(ip,il,iv,iele) + &
                                    somflux(iactive,ipassive,il,iele) + somflux(iactive,islow,il,iele))
                               deepSOM_s(ip,il,iv,iele) = deepSOM_s(ip,il,iv,iele) + (un-reduction_factor) * (deltaSOM1_s(ip,il,iv,iele) + &
                                    somflux(islow,ipassive,il,iele) + somflux(islow,iactive,il,iele))
                               deepSOM_p(ip,il,iv,iele) = deepSOM_p(ip,il,iv,iele) + (un-reduction_factor) * (deltaSOM1_p(ip,il,iv,iele) + &
                                    somflux(ipassive,iactive,il,iele) + somflux(islow,ipassive,il,iele))
                            ELSE                             
                               WRITE (numout,*) ' deltaSOM. We we have an elements which is neither carbon nor nitrogen. We stop.'
                               CALL ipslerr_p(3,'stomate_soil_carbon_discretization_deep_somcycle', &
                                    'We we have an elements which is neither carbon nor nitrogen','','')
                            ENDIF
                         ENDDO
                      END DO
                   ENDIF
                   ! Error checking
                   IF ((reduction_factor.LT. -min_stomate) .OR. (reduction_factor.GT. un)) THEN
                      WRITE(numout,*) 'ERROR: inconsistency with nitrogen reduction',&
                           reduction_factor
                      CALL ipslerr_p(3,'stomate_soil_carbon_discretization_deep_somcycle',&
                           'inconsistency in calculating nitrogen overspending','','')
                   ENDIF
                ENDIF
             ENDIF
          ENDDO
       ENDDO

    ELSE
       
       ! test_adjust_rhetero=.FALSE.
       DO ip = 1, kjpindex
          DO iv = 1, nvm
             IF (  veget_mask_2d(ip,iv) ) THEN
                DO il=1,ngrnd
                   dC1i(ip,iv) = dC1i(ip,iv) + &
                        (deltaSOM1_a(ip,il,iv,icarbon)+deltaSOM1_s(ip,il,iv,icarbon)+deltaSOM1_p(ip,il,iv,icarbon)) * &
                        ( zf_soil(il) - zf_soil(il-1) )
                   n_mineralisation(ip,iv)=n_mineralisation(ip,iv)+ &
                        (deltaSOM1_a(ip,il,iv,initrogen)+deltaSOM1_s(ip,il,iv,initrogen)+deltaSOM1_p(ip,il,iv,initrogen))* &
                        ( zf_soil(il) - zf_soil(il-1) )
                END DO
             END IF
          END DO
       END DO
       
    END IF ! if test_adjust_rhetero

    DO ip = 1, kjpindex
       ! Total CH4 flux
       sfluxCH4_deep(ip) = SUM(veget_max_bg(ip,:)*( CH4ii(ip,:)-CH4i(ip,:)+MG(ip,:)-MT(ip,:) ))/time_step
       ! TotalCO2 flux
       sfluxCO2_deep(ip) = SUM(veget_max_bg(ip,:)*( dC1i(ip,:) + MT(ip,:)*(12./16.) ) )/time_step
    END DO

    ! Not all calculations above account for veget_max. Fluxes and pools
    ! might be calculated for PFTs that do not exist. This results in 
    ! problems with the consistency check of NBP.
    WHERE (veget_max(:,:).GT.min_stomate)
       resp_hetero_soil(:,:) = ( dC1i(:,:) + MT(:,:)*(12./16.) ) *one_day/time_step
       sfluxCH4(:,:) = ( CH4ii(:,:)-CH4i(:,:)+MG(:,:)-MT(:,:) ) *one_day/time_step
    ELSEWHERE
       resp_hetero_soil(:,:) = zero
       sfluxCH4(:,:) = zero
    END WHERE

    ! Not all calculations above account for veget_max. Fluxes and pools
    ! might be calculated for PFTs that do not exist. This results in 
    ! problems with the consistency check of NBP.
    DO il=1,ngrnd
       DO iele=1,nelements
          WHERE (veget_max(:,:).LT.min_stomate)
             deepSOM_a(:,il,:,iele) = zero
             deepSOM_s(:,il,:,iele) = zero
             deepSOM_p(:,il,:,iele) = zero
          END WHERE
       END DO
    END DO

    ! define pft-mean soil C sent to the atmosphere as CO2
    deltaSOM1_a_pftmean(:,:) = 0._r_std
    deltaSOM1_s_pftmean(:,:) = 0._r_std
    deltaSOM1_p_pftmean(:,:) = 0._r_std
    DO iv = 1, nvm
       DO il=1,ngrnd
         deltaSOM1_a_pftmean(:,il)  = deltaSOM1_a_pftmean(:,il)  + &
              deltaSOM1_a(:,il,iv,icarbon) * veget_max_bg(:,iv)
         deltaSOM1_s_pftmean(:,il)    = deltaSOM1_s_pftmean(:,il)    + &
              deltaSOM1_s(:,il,iv,icarbon) * veget_max_bg(:,iv)
         deltaSOM1_p_pftmean(:,il) = deltaSOM1_p_pftmean(:,il) + &
              deltaSOM1_p(:,il,iv,icarbon) * veget_max_bg(:,iv)
       END DO
    END DO

    deltaSOM1_a_pftmean_stock(:) = 0._r_std
    deltaSOM1_s_pftmean_stock(:) = 0._r_std
    deltaSOM1_p_pftmean_stock(:) = 0._r_std
    DO il = 1, ngrnd
      deltaSOM1_a_pftmean_stock(:) = deltaSOM1_a_pftmean_stock(:) + deltaSOM1_a_pftmean(:,il)*(zf_soil(il)-zf_soil(il-1))
      deltaSOM1_s_pftmean_stock(:) = deltaSOM1_s_pftmean_stock(:) + deltaSOM1_s_pftmean(:,il)*(zf_soil(il)-zf_soil(il-1))
      deltaSOM1_p_pftmean_stock(:) = deltaSOM1_p_pftmean_stock(:) + deltaSOM1_p_pftmean(:,il)*(zf_soil(il)-zf_soil(il-1))
    ENDDO
    
    deep_rhC_a_pftmean(:,:) = deltaSOM1_a_pftmean(:,:) / time_step
    deep_rhC_s_pftmean(:,:) = deltaSOM1_s_pftmean(:,:) / time_step
    deep_rhC_p_pftmean(:,:) = deltaSOM1_p_pftmean(:,:) / time_step
    rhC_a_pftmean(:) = deltaSOM1_a_pftmean_stock(:) / time_step
    rhC_s_pftmean(:) = deltaSOM1_s_pftmean_stock(:) / time_step
    rhC_p_pftmean(:) = deltaSOM1_p_pftmean_stock(:) / time_step    

    ! calculate coefficients for cryoturbation calculation
    IF (ok_cryoturb) CALL cryoturbate(kjpindex, time_step, altmax_ind_lastyear, deepSOM_a, deepSOM_s, deepSOM_p, &
         'coefficients', cryoturbation_diff_k_in/(one_day*one_year),bioturbation_diff_k_in/(one_day*one_year), &
         altmax_lastyear, fixed_cryoturbation_depth, veget_max)

    IF (ok_cryoturb) CALL cryoturbate(kjpindex, time_step, altmax_ind_lastyear, deepSOM_a, deepSOM_s, deepSOM_p, &
         'diffuse', cryoturbation_diff_k_in/(one_day*one_year), bioturbation_diff_k_in/(one_day*one_year), &
         altmax_lastyear, fixed_cryoturbation_depth, veget_max)


    ! calculate the coefficients for the next timestep:
    !
    ! get diffusion coefficients: heat capacity,
    !    conductivity, and oxygen diffusivity
    !
    CALL get_gasdiff (kjpindex,hslong,tprof,snow,airvol_snow, &
         totporO2_snow,totporCH4_snow,diffO2_snow,diffCH4_snow, &
         airvol_soil,totporO2_soil,totporCH4_soil,diffO2_soil,diffCH4_soil, snowrho)

    !
    ! calculate the coefficients for the next time step
    !
    CALL soil_gasdiff_main (kjpindex,time_step,index,'coefficients', &
         pb,tsurf,tprof,diffO2_snow,diffCH4_snow, &
         totporO2_snow,totporCH4_snow,O2_snow,CH4_snow,diffO2_soil,diffCH4_soil, &
         totporO2_soil,totporCH4_soil,O2_soil,CH4_soil, zi_snow, zf_snow)

    call calc_vert_int_som(kjpindex, deepSOM_a, deepSOM_s, deepSOM_p, som, som_surf, zf_soil)
    IF (printlev>=3) WRITE(*,*) 'after calc_vert_int_som'

    
    ! XIOS history output
    IF ( .NOT. soilc_isspinup ) THEN
       
       !CALL xios_orchidee_send_field ('fluxCH4',  sfluxCH4)
       !CALL xios_orchidee_send_field ('febul', (febul*one_day))
       !CALL xios_orchidee_send_field ('flupmt',  (flupmt*one_day))
       CALL xios_orchidee_send_field ( 'alt', alt )
       CALL xios_orchidee_send_field ( 'altmax', altmax)
       !CALL xios_orchidee_send_field ( 'sfluxCH4_deep', sfluxCH4_deep)
       !CALL xios_orchidee_send_field ( 'sfluxCO2_deep', sfluxCO2_deep)
       CALL xios_orchidee_send_field ( 'pb', pb)
      
       CALL xios_orchidee_send_field ( 'deep_rhC_a_pftmean', deep_rhC_a_pftmean ) 
       CALL xios_orchidee_send_field ( 'deep_rhC_s_pftmean', deep_rhC_s_pftmean)
       CALL xios_orchidee_send_field ( 'deep_rhC_p_pftmean', deep_rhC_p_pftmean)
       CALL xios_orchidee_send_field ( 'rhC_a_pftmean', rhC_a_pftmean)
       CALL xios_orchidee_send_field ( 'rhC_s_pftmean', rhC_s_pftmean)
       CALL xios_orchidee_send_field ( 'rhC_p_pftmean', rhC_p_pftmean)
       !CALL xios_orchidee_send_field ( 'O2_soil', O2_soil)
       !CALL xios_orchidee_send_field ( 'CH4_soil', CH4_soil)
       !CALL xios_orchidee_send_field ('O2_snow', O2_snow)
       !CALL xios_orchidee_send_field ( 'CH4_snow', CH4_snow)

       !CALL xios_orchidee_send_field ( 'deltaCH4g',  deltaCH4g)
       !CALL xios_orchidee_send_field ( 'deltaCH4',  deltaCH4)

       !CALL xios_orchidee_send_field ( 'heat_Zimov',  heat_Zimov)

       !CALL xios_orchidee_send_field ( 'totporO2_soil', totporO2_soil)
       !CALL xios_orchidee_send_field ( 'diffO2_soil', diffO2_soil)
       !CALL xios_orchidee_send_field ( 'alphaO2_soil',  alphaO2_soil)
       !CALL xios_orchidee_send_field ( 'betaO2_soil',  betaO2_soil)
       !CALL xios_orchidee_send_field ( 'totporCH4_soil', totporCH4_soil)
       !CALL xios_orchidee_send_field ( 'diffCH4_soil', diffCH4_soil)
       !CALL xios_orchidee_send_field ('alphaCH4_soil', alphaCH4_soil)
       !CALL xios_orchidee_send_field ( 'betaCH4_soil', betaCH4_soil)

       !
       IF (perma_peat) THEN
          CALL xios_orchidee_send_field ( 'deepSOM_peat', deepSOM_peat)
          CALL xios_orchidee_send_field ( 'peat_OLT', peat_OLT)
          CALL xios_orchidee_send_field ( 'deepSOM_pt', deepSOM_pt)
       ENDIF

    ENDIF
 !! 5. (Quasi-)Analytical Spin-up

    !! 5.2 Finish to fill MatrixA with fluxes between soil pools

    IF (spinup_analytic) THEN

       ! During the spinup PFT1 (bare soil) should still be empty
       ! for the moment the only way to get carbon/nitrogen in PFT1
       ! is after a LCC. LCC should not be used during spinup.

       fluxout_ave_a(:,:) = zero
       fluxout_ave_s(:,:) = zero
       fluxout_ave_p(:,:) = zero

       DO ivm = 2,nvm

          ! We calculated weigheted average of the fluxes over the soil profile
          ! only for carbon
          DO igrnd= 1, ngrnd
             fluxout_ave_a(:,ivm) = fluxout_ave_a(:,ivm) + fluxout_a(:,igrnd,ivm,icarbon) * &
                                    (zf_soil(igrnd)-zf_soil(igrnd-1))/zf_soil(ngrnd)
             fluxout_ave_s(:,ivm) = fluxout_ave_s(:,ivm) + (fluxout_s(:,igrnd,ivm,icarbon) * &
                                    (zf_soil(igrnd)-zf_soil(igrnd-1))/zf_soil(ngrnd)) 
             fluxout_ave_p(:,ivm) = fluxout_ave_p(:,ivm) + (fluxout_p(:,igrnd,ivm,icarbon) * &
                                    (zf_soil(igrnd)-zf_soil(igrnd-1))/zf_soil(ngrnd))

          ENDDO

          ! flux leaving the active pool
          MatrixA(:,ivm,iactive_pool,iactive_pool) = moins_un * &
                fluxout_ave_a(:,ivm) 
          ! flux received by the active pool from the slow pool
          MatrixA(:,ivm,iactive_pool,islow_pool) = fc(:,islow,iactive,ivm)*fluxout_ave_s(:,ivm) 

          ! flux received by the active pool from the passive pool
          MatrixA(:,ivm,iactive_pool,ipassive_pool) =  fc(:,ipassive,iactive,ivm)*fluxout_ave_p(:,ivm)  

          ! flux leaving the slow pool
          MatrixA(:,ivm,islow_pool,islow_pool) = moins_un * &
                 fluxout_ave_s(:,ivm)
          ! flux received by the slow pool from the active pool
          MatrixA(:,ivm,islow_pool,iactive_pool) = fc(:,iactive,islow,ivm)*fluxout_ave_a(:,ivm)    

          ! flux received by the slow pool from the passive pool
          MatrixA(:,ivm,islow_pool,ipassive_pool) = fc(:,ipassive,islow,ivm)*fluxout_ave_p(:,ivm)  

          ! flux leaving the passive pool
          MatrixA(:,ivm,ipassive_pool,ipassive_pool) =  moins_un * &
                 fluxout_ave_p(:,ivm)
          ! flux received by the passive pool from the active pool
          MatrixA(:,ivm,ipassive_pool,iactive_pool) = fc(:,iactive,ipassive,ivm)*fluxout_ave_a(:,ivm)

          ! flux received by the passive pool from the slow pool
          MatrixA(:,ivm,ipassive_pool,islow_pool) = fc(:,islow,ipassive,ivm)*fluxout_ave_s(:,ivm) 

          IF (printlev_loc>=4) WRITE(numout,*)'Finish to fill MatrixA'

!!$          ! Error checking
!!$          ! If there is no veget_max, the litter pools should be empty
!!$          error_count(:) = zero
!!$          WHERE ( SUM(ABS(som(:,:,ivm,initrogen)),2) .GT. min_stomate .AND. &
!!$               veget_max(:,ivm) .LT. min_stomate)
!!$             error_count(:) = un
!!$          END WHERE
!!$          IF (SUM(error_count(:)).GT.zero) THEN
!!$             DO ipts = 1,npts
!!$                IF (error_count(ipts).GT.zero) THEN
!!$                   WRITE(numout,*) 'There is litter but no veget_max in pixel and PFT, ', ipts, ivm
!!$                   WRITE(numout,*) 'veget_max, ', veget_max(ipts,ivm)
!!$                   DO icarb = 1,ncarb
!!$                      WRITE(numout,*) 'som C, ', ipts, icarb, ivm, som(ipts,icarb,ivm,icarbon)
!!$                      WRITE(numout,*) 'som N, ', ipts, icarb, ivm, som(ipts,icarb,ivm,initrogen)
!!$                   END DO
!!$                   CALL ipslerr_p(3, 'stomate_soilcarbon','Problem when preparing for the analytical spinup',&
!!$                        'veget_max indicates that the PFT is not available','however there is soilcarbon')
!!$                END IF
!!$             END DO
!!$          END IF
!!$          !-

          WHERE (SUM(deepSOM_a(:,:,ivm,initrogen),DIM=2) .GT. min_stomate)
             CN_som_litter_longterm(:,ivm,iactive_pool)  = ( CN_som_litter_longterm(:,ivm,iactive_pool) * (tau_CN_longterm-dt) &
                  + SUM(deepSOM_a(:,:,ivm,icarbon),DIM=2)/ SUM(deepSOM_a(:,:,ivm,initrogen),DIM=2) * dt)/ (tau_CN_longterm)
          ENDWHERE

          WHERE(SUM(deepSOM_s(:,:,ivm,initrogen),DIM=2) .GT. min_stomate)
             CN_som_litter_longterm(:,ivm,islow_pool)    = ( CN_som_litter_longterm(:,ivm,islow_pool) * (tau_CN_longterm-dt) &
                  + SUM(deepSOM_s(:,:,ivm,icarbon),DIM=2)/SUM(deepSOM_s(:,:,ivm,initrogen),DIM=2) * dt)/ (tau_CN_longterm)
          ENDWHERE

          WHERE(SUM(deepSOM_p(:,:,ivm,initrogen),DIM=2) .GT. min_stomate)
             CN_som_litter_longterm(:,ivm,ipassive_pool) =  ( CN_som_litter_longterm(:,ivm,ipassive_pool) * (tau_CN_longterm-dt) &
                  + SUM(deepSOM_p(:,:,ivm,icarbon),DIM=2)/SUM(deepSOM_p(:,:,ivm,initrogen),DIM=2) * dt)/ (tau_CN_longterm)
          ENDWHERE

       ENDDO ! Loop over # PFTS

       ! Add Identity for each submatrix(7,7) 
       DO inbpools = 1,nbpools
          MatrixA(:,:,inbpools,inbpools) = MatrixA(:,:,inbpools,inbpools) + un
       END DO
    ENDIF ! (spinup_analytic)

    ! Error checking
    ! When running the model without CH4 processes these pools
    ! and fluxes should be zero. If not the mass balance will
    ! not be closed.
    IF (.NOT. ok_methane) THEN
       IF (SUM(sfluxCH4_deep(:)).GT.zero) THEN
          CALL ipslerr_p(3,'stomate_soil_carbon_discretization',&
               'sfluxCH4_deep is not zero','this should not be the case',&
               'because ok_methane is .FALSE.')
       END IF
       IF (SUM(SUM(sfluxCH4(:,:),2)).GT.zero) THEN
          CALL ipslerr_p(3,'stomate_soil_carbon_discretization',&
               'sfluxCH4 is not zero','this should not be the case',&
               'because ok_methane is .FALSE.')
       END IF
       IF (SUM(SUM(SUM(CH4_snow(:,:,:),2),2)).GT.zero) THEN
          CALL ipslerr_p(3,'stomate_soil_carbon_discretization',&
               'CH4_snow is not zero','this should not be the case',&
               'because ok_methane is .FALSE.')
       END IF
       IF (SUM(SUM(SUM(CH4_soil(:,:,:),2),2)).GT.zero) THEN
          CALL ipslerr_p(3,'stomate_soil_carbon_discretization',&
               'CH4_soil is not zero','this should not be the case',&
               'because ok_methane is .FALSE.')
       END IF
    END IF

    !! Check mass balance closure
    IF (err_act.GT.1) THEN

       !! 4.1.1 Calculate components of the mass balance
       pool_end(:,:,:) = zero 
       pool_end(:,:,initrogen) = pool_end(:,:,initrogen) + &
            (n_mineralisation(:,:) * veget_max(:,:) ) 
       DO iele = 1,nelements
          DO igrnd = 1,ngrnd
             pool_end(:,:,iele) = pool_end(:,:,iele) + &
                  (deepSOM_a(:,igrnd,:,iele) + deepSOM_s(:,igrnd,:,iele) + deepSOM_p(:,igrnd,:,iele)) * &
                  (zf_soil(igrnd)-zf_soil(igrnd-1)) * veget_max(:,:)
          END DO
       END DO

       !! 4.1.2 Calculate mass balance
       check_intern(:,:,:,:) = zero
       check_intern(:,:,iatm2land,:) = zero
       check_intern(:,:,iland2atm,icarbon) = -un * (resp_hetero_soil(:,:) + sfluxCH4(:,:)) * veget_max(:,:) * dt 
       check_intern(:,:,iland2atm,initrogen) = -un * zero
       check_intern(:,:,ilat2out,:) = -un *zero
       check_intern(:,:,ilat2in,:) =  zero
       check_intern(:,:,ipoolchange,:) = -un * (pool_end(:,:,:) - &
            pool_start(:,:,:))
       closure_intern = zero
       DO imbc = 1,nmbcomp
          DO iele = 1,nelements
             closure_intern(:,:,iele) = closure_intern(:,:,iele) + &
                  check_intern(:,:,imbc,iele)
          END DO
       END DO

       CALL check_mass_balance("stomate_soil_carbon_discretization_deep_somcycle", closure_intern, &
            kjpindex, pool_end, pool_start, veget_max, 'pft')

    END IF ! err_act.gt.1

    IF (printlev>=3) WRITE(*,*) 'cdk: leaving  stomate_soil_carbon_discretization_deep_somcycle'
    
    IF ( firstcall )  firstcall = .FALSE.


  END SUBROUTINE stomate_soil_carbon_discretization_deep_somcycle
  
!!
!================================================================================================================================
!! SUBROUTINE   : altcalc
!!
!>\BRIEF        This routine calculate active layer thickness
!!
!! DESCRIPTION : 
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : alt
!! 
!! REFERENCE(S) : None
!!
!! FLOWCHART11    : None
!! \n
!_
!================================================================================================================================  
  SUBROUTINE altcalc (kjpindex,time_step,scnd, temp, zprof, alt, alt_ind, altmax, altmax_ind, &
        altmax_lastyear, altmax_ind_lastyear)

  !! 0. Variable and parameter declaration

    !! 0.1  Input variables
    INTEGER(i_std), INTENT(in) 			         :: kjpindex
    REAL(r_std), INTENT(in)                              :: time_step           !! time step in seconds
    REAL(r_std), INTENT(in)			         :: scnd                !! model time & time step
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)            :: temp                !! soil temperature
    REAL(r_std), DIMENSION(:), INTENT(in)   	         :: zprof               !! soil depths (m)

    !! 0.2 Output variables

    REAL(r_std), DIMENSION(:,:), INTENT(out)	         :: alt	                !! active layer thickness  
    INTEGER, DIMENSION(:,:), INTENT(out)	         :: alt_ind	        !! active layer index  
    
    !! 0.3 Modified variables

    REAL(r_std), DIMENSION(:,:),INTENT(inout)	         :: altmax_lastyear     !! Maximum active-layer thickness
    REAL(r_std), DIMENSION(:,:),INTENT(inout)            :: altmax              !! Maximum active-layer thickness
    INTEGER(i_std), DIMENSION(:,:),INTENT(inout)	 :: altmax_ind          !! Maximum over the year active-layer index
    INTEGER(i_std), DIMENSION(:,:),INTENT(inout)         :: altmax_ind_lastyear !! Maximum over the year active-layer index

    !! 0.4 Local variables

    INTEGER				                 :: ix,iz,il,iv         !! grid indices
    LOGICAL, SAVE 				         :: firstcall = .TRUE.
!$OMP THREADPRIVATE(firstcall)
    INTEGER(i_std), SAVE  			         :: id, id2
!$OMP THREADPRIVATE(id)
!$OMP THREADPRIVATE(id2)
    LOGICAL, SAVE			                 :: check = .FALSE.
!$OMP THREADPRIVATE(check)
    LOGICAL, SAVE			                 :: newaltcalc
!$OMP THREADPRIVATE(newaltcalc)
    LOGICAL, DIMENSION(kjpindex,nvm)                     :: inalt, bottomlevelthawed
    INTEGER                                              :: lev

    
    IF ( firstcall )  THEN
       firstcall = .FALSE.

       !Config Key   = newaltcalc
       !Config Desc  = calculate alt ?
       !Config Def   = n
       !Config If    = OK_SOIL_CARBON_DISCRETIZATION
       !Config Help  = 
       !Config Unit  = [flag]
       newaltcalc = .FALSE.
       CALL getin_p('newaltcalc', newaltcalc)
    ENDIF !firstcall

       ! all other timesteps
       IF ( .NOT. newaltcalc ) THEN
          DO ix = 1, kjpindex
             DO iv = 1, nvm
                IF ( veget_mask_2d(ix,iv) ) THEN
                   iz = 1
                   DO WHILE( temp(ix,iz,iv) > ZeroCelsius .AND. iz < ngrnd )
                      iz = iz + 1	  
                   END DO
                   IF( iz == 1 ) THEN 
                      ! it means that all is frozen
                      alt(ix,iv) = zero
                   ELSE
                      alt(ix,iv) = zprof(iz-1)
                   END IF
                   alt_ind(ix,iv) = iz-1
                END IF
             END DO
          END DO
       ELSE          
          ! initialize for pfts that don't exist
          alt(:,:) = zprof(ngrnd)  
          bottomlevelthawed(:,:) = .FALSE.
          ! start from bottom and work up instead
          WHERE (temp(:,ngrnd,:) > ZeroCelsius ) 
             bottomlevelthawed(:,:) = .TRUE.
             alt(:,:) = zprof(ngrnd)
             alt_ind(:,:) = ngrnd
          END WHERE
          inalt(:,:) = .FALSE.
          DO iz = 1, ngrnd - 1
             lev = ngrnd - iz
             WHERE ( temp(:,lev,:) > ZeroCelsius .AND. .NOT. inalt(:,:) .AND. .NOT. bottomlevelthawed(:,:) )
                inalt(:,:) = .TRUE.
                alt(:,:) = zprof(lev)
                alt_ind(:,:) = lev
             ELSEWHERE ( temp(:,lev,:) <= ZeroCelsius .AND. inalt(:,:) .AND. .NOT. bottomlevelthawed(:,:) )
                inalt(:,:) = .FALSE.
             END WHERE
          END DO
          WHERE ( .NOT. inalt .AND. .NOT. bottomlevelthawed(:,:) ) 
             alt(:,:) = zero
             alt_ind(:,:) = 0
          END WHERE
       ENDIF

       ! debug
       IF ( check ) THEN
          IF (ANY(alt(:,:) .GT. zprof(ngrnd))) THEN
             WRITE(*,*) 'error: alt greater than soil depth.'
          ENDIF
       ENDIF

       ! Maximum over the year active layer thickness
       WHERE ( ( alt(:,:) .GT. altmax(:,:) ) .AND. veget_mask_2d(:,:)  ) 
          altmax(:,:) = alt(:,:)
          altmax_ind(:,:) = alt_ind(:,:)
       ENDWHERE
       
       IF ( .NOT. soilc_isspinup ) THEN
             ! do it on the second timestep, that way when we are writing restart files it is not done before that!
             ! now we are doing daily permafrost calcs, so just run it on the second day.
             ! JG: I don't see any reason to use the second day. This should probably be changed into FirstTsYear (variable in module time)
          IF (FirstTsYear) THEN
            ! Reinitialize ALT_max
             altmax_lastyear(:,:) = altmax(:,:)
             altmax_ind_lastyear(:,:) = altmax_ind(:,:)
             altmax(:,:) = alt(:,:)
             altmax_ind(:,:) = alt_ind(:,:)
          END IF
       ELSE

             ! for spinup, best to set altmax_lastyear to altmax, and not boter to reset since every year is the same, 
             ! and if you try to do so, it doesn't work properly --  06 may 2010
          altmax_lastyear(:,:) = altmax(:,:)
          altmax_ind_lastyear(:,:) = altmax_ind(:,:)
       END IF

    IF (printlev>=3) WRITE(*,*) 'leaving  altcalc'
  END SUBROUTINE altcalc
  
!!
!================================================================================================================================
!! SUBROUTINE   : soil_gasdiff_main
!!
!>\BRIEF        This routine calculate oxygen and methane in the snow/soil medium 
!!
!! DESCRIPTION : 
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!! 
!! REFERENCE(S) : None
!!
!! FLOWCHART11    : None
!! \n
!_
!================================================================================================================================   
  SUBROUTINE soil_gasdiff_main( kjpindex,time_step,index,cur_action, &
       psol,tsurf,tprof,diffO2_snow,diffCH4_snow, &
       totporO2_snow,totporCH4_snow,O2_snow,CH4_snow,diffO2_soil,diffCH4_soil, &
       totporO2_soil,totporCH4_soil,O2_soil,CH4_soil, zi_snow, zf_snow)

  !! 0. Variable and parameter declaration

    !! 0.1  Input variables

    INTEGER(i_std), INTENT(in)                    :: kjpindex           !! number of grid points 
    REAL(r_std), INTENT(in)                       :: time_step          !! time step in seconds
    CHARACTER(LEN=*), INTENT(in)                  :: cur_action         !! what to do 
    REAL(r_std), DIMENSION(:), INTENT(in)         :: psol               !! surface pressure (Pa)
    REAL(r_std), DIMENSION(:), INTENT(in)         :: tsurf              !! Surface temperature (K)
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)     :: tprof              !! Soil temperature (K)
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)     :: diffO2_snow        !! oxygen diffusivity (m**2/s)
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)     :: diffCH4_snow       !! methane diffusivity (m**2/s)
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)     :: totporO2_snow      !! total O2 porosity (Tans, 1998)
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)     :: totporCH4_snow     !! total CH4 porosity (Tans, 1998)
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)     :: diffO2_soil        !! oxygen diffusivity (m**2/s)
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)     :: diffCH4_soil       !! methane diffusivity (m**2/s)
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)     :: totporO2_soil      !! total O2 porosity (Tans, 1998)
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)     :: totporCH4_soil     !! total CH4 porosity (Tans, 1998)
    INTEGER(i_std),DIMENSION(:),INTENT(in)        :: index              !! Indeces of permafrost points on the map
    REAL(r_std), DIMENSION(kjpindex,0:nsnow,nvm), INTENT(in)   :: zf_snow            !! depths of full levels (m)
    REAL(r_std), DIMENSION(kjpindex,nsnow,nvm), INTENT(in)     :: zi_snow            !! depths of intermediate levels (m)

    !! 0.2  Output variables

    !! 0.3  Modified variables

    REAL(r_std), DIMENSION(kjpindex,nsnow,nvm), INTENT(inout)  :: O2_snow            !! oxygen (g O2/m**3 air)
    REAL(r_std), DIMENSION(kjpindex,nsnow,nvm), INTENT(inout)  :: CH4_snow           !! methane (g CH4/m**3 air)
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm), INTENT(inout)  :: O2_soil            !! oxygen (g O2/m**3 air)
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm), INTENT(inout)  :: CH4_soil           !! methane (g CH4/m**3 air)
  
    !! 0.4 local variables

    CHARACTER(LEN=50), SAVE        :: last_action = 'not called'
!$OMP THREADPRIVATE(last_action)
    
    
    ! 1. ensure that we do not repeat actions
    !
    IF ( TRIM(cur_action) .EQ. TRIM(last_action) ) THEN
       !
       CALL ipslerr_p(3, 'soil_gasdiff_main','CANNOT TAKE THE SAME ACTION TWICE: ',cur_action, last_action)
       !
    ENDIF
    !
    ! 2. decide what to do
    !
    IF ( TRIM(cur_action) .EQ. 'initialize' ) THEN
       !
       ! 2.1 initialize
       !
       IF ( TRIM(last_action) .NE. 'not called' ) THEN
          !
          CALL ipslerr_p(3, 'soil_gasdiff_main','SOIL MODEL CANNOT BE INITIALIZED TWICE.', '', '')
          !
       ENDIF
       !
       CALL soil_gasdiff_alloc( kjpindex )
       !
    ELSEIF ( TRIM(cur_action) .EQ. 'diffuse' ) THEN
       !
       ! 2.2 calculate soil temperatures
       !
       CALL soil_gasdiff_diff( kjpindex,time_step,index,psol,tsurf, O2_snow, CH4_snow, O2_soil, CH4_soil)
       !
    ELSEIF ( TRIM(cur_action) .EQ. 'coefficients' ) THEN
       !
       ! 2.3 calculate coefficients (heat flux and apparent surface heat capacity)
       !
       CALL soil_gasdiff_coeff( kjpindex,time_step,tprof,O2_snow,CH4_snow, &
            diffO2_snow,diffCH4_snow,totporO2_snow,totporCH4_snow,O2_soil,CH4_soil, &
            diffO2_soil,diffCH4_soil,totporO2_soil,totporCH4_soil, zi_snow, zf_snow)
       !
    ELSE
       !
       ! 2.4 do not know this action
       !
       CALL ipslerr_p(3,'soil_gasdiff_main', 'This action does not exists and must be implemented', &
                        'or its wrong:',cur_action)
       !
    ENDIF
    !
    ! 2.5 keep last action in mind
    !
    last_action = TRIM(cur_action)
    
    IF (printlev>=3) WRITE(*,*) 'leaving  soil_gasdiff_main'
  END SUBROUTINE soil_gasdiff_main
 
!!
!================================================================================================================================
!! SUBROUTINE   : soil_gasdiff_alloc
!!
!>\BRIEF        This routine allocate arrays related to oxygen and methane in the snow/soil medium 
!!
!! DESCRIPTION : 
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!! 
!! REFERENCE(S) : None
!!
!! FLOWCHART11    : None
!! \n
!_
!================================================================================================================================   
  SUBROUTINE soil_gasdiff_alloc( kjpindex )
   
  !! 0. Variable and parameter declaration

    !! 0.1  Input variables
 
    INTEGER(i_std), INTENT(in)                             :: kjpindex

    !! 0.2 Output variables

    !! 0.3 Modified variables
    
    !! 0.4 local variables

    INTEGER(i_std)                                         :: ier
    
    ! Allocate the variables that need to be saved after soil_gasdiff_coeff

      ALLOCATE (alphaO2_soil(kjpindex,ngrnd,nvm),stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3,'soil_gasdiff_alloc', 'Pb in alloc for alphaO2_soil','','')

      ALLOCATE (betaO2_soil(kjpindex,ngrnd,nvm),stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3,'soil_gasdiff_alloc', 'Pb in alloc for betaO2_soil','','')

      ALLOCATE (alphaCH4_soil(kjpindex,ngrnd,nvm),stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3,'soil_gasdiff_alloc', 'Pb in alloc for alphaCH4_soil','','')

      ALLOCATE (betaCH4_soil(kjpindex,ngrnd,nvm),stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3,'soil_gasdiff_alloc', 'Pb in alloc for betaCH4_soil','','')

      ALLOCATE (alphaO2_snow(kjpindex,nsnow,nvm),stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3,'soil_gasdiff_alloc', 'Pb in alloc for alphaO2_snow','','')

      ALLOCATE (betaO2_snow(kjpindex,nsnow,nvm),stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3,'soil_gasdiff_alloc', 'Pb in alloc for betaO2_snow','','')

      ALLOCATE (alphaCH4_snow(kjpindex,nsnow,nvm),stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3,'soil_gasdiff_alloc', 'Pb in alloc for alphaCH4_snow','','')

      ALLOCATE (betaCH4_snow(kjpindex,nsnow,nvm),stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3,'soil_gasdiff_alloc', 'Pb in alloc for betaCH4_snow','','')

      ALLOCATE (zf_coeff_snow(kjpindex,0:nsnow,nvm),stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3,'soil_gasdiff_alloc', 'Pb in alloc for zf_coeff_snow','','')

      ALLOCATE (zi_coeff_snow(kjpindex,nsnow,nvm),stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3,'soil_gasdiff_alloc', 'Pb in alloc for zi_coeff_snow','','')

      ALLOCATE (mu_snow(kjpindex,nvm),stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3,'soil_gasdiff_alloc', 'Pb in alloc for mu_snow','','')


      alphaO2_soil(:,:,:) = zero
      betaO2_soil(:,:,:) = zero
      alphaCH4_soil(:,:,:) = zero
      betaCH4_soil(:,:,:) = zero
      alphaO2_snow(:,:,:) = zero
      betaO2_snow(:,:,:) = zero
      alphaCH4_snow(:,:,:) = zero
      betaCH4_snow(:,:,:) = zero
      zf_coeff_snow(:,:,:) = zero
      zi_coeff_snow(:,:,:) = zero
      mu_snow(:,:) = zero
    
  END SUBROUTINE soil_gasdiff_alloc
 
!!
!================================================================================================================================
!! SUBROUTINE   : soil_gasdiff_coeff
!!
!>\BRIEF        This routine calculate coeff related to gas diffuvisity
!!
!! DESCRIPTION : 
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!! 
!! REFERENCE(S) : None
!!
!! FLOWCHART11    : None
!! \n
!_
!================================================================================================================================   
  
  SUBROUTINE soil_gasdiff_coeff( kjpindex,time_step,tprof,O2_snow,CH4_snow, &
       diffO2_snow,diffCH4_snow,totporO2_snow,totporCH4_snow,O2_soil,CH4_soil, &
       diffO2_soil,diffCH4_soil,totporO2_soil,totporCH4_soil, zi_snow, zf_snow)


  !! 0. Variable and parameter declaration

    !! 0.1  Input variables

    INTEGER(i_std), INTENT(in)                    :: kjpindex            !! number of grid points 
    REAL(r_std), INTENT(in)                       :: time_step           !! time step in seconds
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)     :: tprof               !! Soil temperature (K)
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)     :: diffO2_snow         !! oxygen diffusivity (m**2/s)
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)     :: diffCH4_snow        !! methane diffusivity (m**2/s)
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)     :: totporO2_snow       !! total O2 porosity (Tans, 1998)
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)     :: totporCH4_snow      !! total CH4 porosity (Tans, 1998)
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)     :: diffO2_soil         !! oxygen diffusivity (m**2/s)
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)     :: diffCH4_soil        !! methane diffusivity (m**2/s)
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)     :: totporO2_soil       !! total O2 porosity (Tans, 1998)
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)     :: totporCH4_soil      !! total CH4 porosity (Tans, 1998)
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)     :: O2_snow             !! oxygen (g O2/m**3 air)
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)     :: CH4_snow            !! methane (g CH4/m**3 air)
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)     :: O2_soil             !! oxygen (g O2/m**3 air)
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)     :: CH4_soil            !! methane (g CH4/m**3 air)
    REAL(r_std), DIMENSION(:,0:,:), INTENT(in)    :: zf_snow             !! depths of full levels (m)
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)     :: zi_snow             !! depths of intermediate levels (m)

    !! 0.2  Output variables

    !! 0.3  Modified variables

    !! 0.4 local variables

    REAL(r_std), DIMENSION(kjpindex,nsnow,nvm)         :: xcO2_snow,xdO2_snow
    REAL(r_std), DIMENSION(kjpindex,nsnow,nvm)         :: xcCH4_snow,xdCH4_snow
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm)         :: xcO2_soil,xdO2_soil
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm)         :: xcCH4_soil,xdCH4_soil
    INTEGER(i_std)                                     :: il
    REAL(r_std), DIMENSION(kjpindex,nvm)               :: xeO2,xeCH4
    LOGICAL, DIMENSION(kjpindex,nvm)                   :: snow_height_mask_2d


    ! loop over materials (soil, snow), beginning at the bottom
    !
    ! 1. define useful variables linked to geometry and physical properties
    !
    ! 1.1 normal levels
    !
    ! default value if inexistent
    xcO2_snow(:,:,:) = 0
    xdO2_snow(:,:,:) = 0
    xcCH4_snow(:,:,:) = 0
    xdCH4_snow(:,:,:) = 0
    xcO2_soil(:,:,:) = 0
    xdO2_soil(:,:,:) = 0
    xcCH4_soil(:,:,:) = 0
    xdCH4_soil(:,:,:) = 0
    xeO2 = 0
    xeCH4 = 0
    !
    snow_height_mask_2d(:,:) = ( heights_snow(:,:) .GT. hmin_tcalc )
    !
    DO il = 1,nsnow-1
       !
       WHERE ( snow_height_mask_2d(:,:) .AND. veget_mask_2d(:,:) )
          !
          xcO2_snow(:,il,:) = ( zf_snow(:,il,:) - zf_snow(:,il-1,:) ) * &
               totporO2_snow(:,il,:) / time_step
          xcCH4_snow(:,il,:) = ( zf_snow(:,il,:) - zf_snow(:,il-1,:) ) * &
               totporCH4_snow(:,il,:) / time_step
          !
          xdO2_snow(:,il,:) = diffO2_snow(:,il,:) /  &
               (zi_snow(:,il+1,:)-zi_snow(:,il,:))
          xdCH4_snow(:,il,:) = diffCH4_snow(:,il,:) /  &
               (zi_snow(:,il+1,:)-zi_snow(:,il,:))
          !
       ENDWHERE
    END DO
    !
    DO il = 1,ngrnd-1
       !
       WHERE ( veget_mask_2d(:,:) )
          !
          xcO2_soil(:,il,:) = ( zf_soil(il) - zf_soil(il-1) ) * &
               totporO2_soil(:,il,:) / time_step
          xcCH4_soil(:,il,:) = ( zf_soil(il) - zf_soil(il-1) ) * &
               totporCH4_soil(:,il,:) / time_step
          !
          xdO2_soil(:,il,:) = diffO2_soil(:,il,:) /  &
               (zi_soil(il+1)-zi_soil(il))
          xdCH4_soil(:,il,:) = diffCH4_soil(:,il,:) /  &
               (zi_soil(il+1)-zi_soil(il))
          !
       ENDWHERE
       !
    ENDDO
    !
    ! 1.2 for the lower boundary, define a similar geometric variable.
    !
    !snow
    !
    WHERE ( snow_height_mask_2d(:,:) .AND. veget_mask_2d(:,:) ) 
       xcO2_snow(:,nsnow,:) = ( zf_snow(:,nsnow,:) -  &
            zf_snow(:,nsnow-1,:) ) *  &
            totporO2_snow(:,nsnow,:) / time_step
       xdO2_snow(:,nsnow,:) = diffO2_snow(:,nsnow,:) /  &
            ( zi_soil(1) +  &
            zf_snow(:,nsnow,:) - zi_snow(:,nsnow,:) )
       xcCH4_snow(:,nsnow,:) = ( zf_snow(:,nsnow,:) -  &
            zf_snow(:,nsnow-1,:) ) * &
            totporCH4_snow(:,nsnow,:) / time_step
       xdCH4_snow(:,nsnow,:) = diffCH4_snow(:,nsnow,:) /  &
            ( zi_soil(1) +  &
            zf_snow(:,nsnow,:) - zi_snow(:,nsnow,:) )
    ENDWHERE
    !
    ! soil
    ! 
    WHERE (  veget_mask_2d(:,:) ) ! removed heights_soil logic
       xcO2_soil(:,ngrnd,:) =  &
            ( zf_soil(ngrnd) - zf_soil(ngrnd-1) ) * &
            totporO2_soil(:,ngrnd,:) / time_step
       xdO2_soil(:,ngrnd,:) = diffO2_soil(:,ngrnd,:) /  &
            ( zf_soil(ngrnd) - zi_soil(ngrnd) )
       xcCH4_soil(:,ngrnd,:) =  &
            ( zf_soil(ngrnd) - zf_soil(ngrnd-1) ) * &
            totporCH4_soil(:,ngrnd,:) / time_step
       xdCH4_soil(:,ngrnd,:) = diffCH4_soil(:,ngrnd,:) /  &
            ( zf_soil(ngrnd) - zi_soil(ngrnd) )
    ENDWHERE    
    !
    ! 1.3 extrapolation factor from first levels to surface
    !
    WHERE ( snow_height_mask_2d(:,:)  .AND. veget_mask_2d(:,:) )
       mu_snow(:,:) = zi_snow(:,1,:) / ( zi_snow(:,2,:) - zi_snow(:,1,:) )
    ELSEWHERE ( veget_mask_2d(:,:) )
       mu_snow(:,:) = .5 ! any value
    ENDWHERE
    !
    mu_soil = zi_soil(1) / ( zi_soil(2) - zi_soil(1) )
    !
    ! 2. bottom level: treatment depends on lower boundary condition
    !
    ! soil
    !
    WHERE ( veget_mask_2d(:,:) ) ! removed heights_soil logic
       !
       xeO2(:,:) = xcO2_soil(:,ngrnd,:) + xdO2_soil(:,ngrnd-1,:)
       xeCH4(:,:) = xcCH4_soil(:,ngrnd,:) + xdCH4_soil(:,ngrnd-1,:)
       !
       alphaO2_soil(:,ngrnd-1,:) = xdO2_soil(:,ngrnd-1,:) / xeO2(:,:)
       alphaCH4_soil(:,ngrnd-1,:) = xdCH4_soil(:,ngrnd-1,:)  &
            / xeCH4(:,:)
       !
       betaO2_soil(:,ngrnd-1,:) =  &
            (xcO2_soil(:,ngrnd,:)*O2_soil(:,ngrnd,:))/xeO2(:,:)
       betaCH4_soil(:,ngrnd-1,:) =  &
            (xcCH4_soil(:,ngrnd,:)*CH4_soil(:,ngrnd,:))/xeCH4(:,:)
       !
    ENDWHERE
    !
    !snow
    !
    WHERE ( snow_height_mask_2d(:,:) .AND. veget_mask_2d(:,:) )
       !
       ! dernier niveau
       !
       xeO2(:,:) = xcO2_soil(:,1,:) + &
            (1.-alphaO2_soil(:,1,:))*xdO2_soil(:,1,:) +  &
            xdO2_snow(:,nsnow,:)
       xeCH4(:,:) = xcCH4_soil(:,1,:) + &
            (1.-alphaCH4_soil(:,1,:))*xdCH4_soil(:,1,:) + &
            xdCH4_snow(:,nsnow,:)
       !
       alphaO2_snow(:,nsnow,:) = xdO2_snow(:,nsnow,:)/xeO2(:,:)
       alphaCH4_snow(:,nsnow,:) = xdCH4_snow(:,nsnow,:) &
            /xeCH4(:,:)
       !
       betaO2_snow(:,nsnow,:) =  &
            ( xcO2_soil(:,1,:)*O2_soil(:,1,:) + &
            xdO2_soil(:,1,:)*betaO2_soil(:,1,:) ) &
            / xeO2(:,:)
       betaCH4_snow(:,nsnow,:) =  &
            ( xcCH4_soil(:,1,:)*CH4_soil(:,1,:) + &
            xdCH4_soil(:,1,:)*betaCH4_soil(:,1,:) ) &
            / xeCH4(:,:)
       !
       ! avant-dernier niveau
       !
       xeO2(:,:) = xcO2_snow(:,nsnow,:) + &
            (1.-alphaO2_snow(:,nsnow,:))*xdO2_snow(:,nsnow,:) + &
            xdO2_snow(:,nsnow-1,:)
       xeCH4(:,:) = xcCH4_snow(:,nsnow,:) + &
            (1.-alphaCH4_snow(:,nsnow,:))*xdCH4_snow(:,nsnow,:) &
            + xdCH4_snow(:,nsnow-1,:)
       !
       alphaO2_snow(:,nsnow-1,:) =  &
            xdO2_snow(:,nsnow-1,:) / xeO2(:,:)
       alphaCH4_snow(:,nsnow-1,:) =  &
            xdCH4_snow(:,nsnow-1,:) / xeCH4(:,:)
       !
       betaO2_snow(:,nsnow-1,:) = &
            ( xcO2_snow(:,nsnow,:)*O2_snow(:,nsnow,:) + &
            xdO2_snow(:,nsnow,:)*betaO2_snow(:,nsnow,:) ) &
            / xeO2(:,:)
       betaCH4_snow(:,nsnow-1,:) = &
            ( xcCH4_snow(:,nsnow,:)*CH4_snow(:,nsnow,:) + &
            xdCH4_snow(:,nsnow,:)*betaCH4_snow(:,nsnow,:) ) &
            / xeCH4(:,:)
       !
    ELSEWHERE ( veget_mask_2d(:,:) )
       !
       alphaO2_snow(:,nsnow,:) = 1.
       alphaCH4_snow(:,nsnow,:) = 1.
       betaO2_snow(:,nsnow,:) = zero
       betaCH4_snow(:,nsnow,:) = zero
       !
       alphaO2_snow(:,nsnow-1,:) = 1.
       alphaCH4_snow(:,nsnow-1,:) = 1.
       betaO2_snow(:,nsnow-1,:) = zero
       betaCH4_snow(:,nsnow-1,:) = zero
       !
    ENDWHERE
    !    
            
    !
    ! 3. the other levels
    !
    DO il = nsnow-2,1,-1 !snow
       !
       WHERE ( snow_height_mask_2d(:,:) .AND. veget_mask_2d(:,:) )
          !
          xeO2(:,:) = xcO2_snow(:,il+1,:) +  &
               (1.-alphaO2_snow(:,il+1,:))*xdO2_snow(:,il+1,:) + xdO2_snow(:,il,:)
          xeCH4(:,:) = xcCH4_snow(:,il+1,:) +  &
               (1.-alphaCH4_snow(:,il+1,:))*xdCH4_snow(:,il+1,:) +  &
               xdCH4_snow(:,il,:)
          !
          alphaO2_snow(:,il,:) = xdO2_snow(:,il,:) / xeO2(:,:)
          alphaCH4_snow(:,il,:) = xdCH4_snow(:,il,:) / xeCH4(:,:)
          !
          betaO2_snow(:,il,:) =  &
               ( xcO2_snow(:,il+1,:)*O2_snow(:,il+1,:) +  &
               xdO2_snow(:,il+1,:)*betaO2_snow(:,il+1,:) ) / xeO2(:,:)
          betaCH4_snow(:,il,:) =  &
               ( xcCH4_snow(:,il+1,:)*CH4_snow(:,il+1,:) +  &
               xdCH4_snow(:,il+1,:)*betaCH4_snow(:,il+1,:) ) / xeCH4(:,:)
          !
       ELSEWHERE ( veget_mask_2d(:,:) )
          !
          alphaO2_snow(:,il,:) = 1.
          alphaCH4_snow(:,il,:) = 1.
          !
          betaO2_snow(:,il,:) = zero
          betaCH4_snow(:,il,:) = zero
          !
       ENDWHERE
       !
    ENDDO
    !
    DO il = ngrnd-2,1,-1 !soil
       !
       WHERE ( veget_mask_2d(:,:) ) !removed heights_soil logic
          !
          xeO2(:,:) = xcO2_soil(:,il+1,:) +  &
               (1.-alphaO2_soil(:,il+1,:))*xdO2_soil(:,il+1,:) + xdO2_soil(:,il,:)
          xeCH4(:,:) = xcCH4_soil(:,il+1,:) +  &
               (1.-alphaCH4_soil(:,il+1,:))*xdCH4_soil(:,il+1,:) +  &
               xdCH4_soil(:,il,:)
          !
          alphaO2_soil(:,il,:) = xdO2_soil(:,il,:) / xeO2(:,:)
          alphaCH4_soil(:,il,:) = xdCH4_soil(:,il,:) / xeCH4(:,:)
          !
          betaO2_soil(:,il,:) =  &
               ( xcO2_soil(:,il+1,:)*O2_soil(:,il+1,:) +  &
               xdO2_soil(:,il+1,:)*betaO2_soil(:,il+1,:) ) / xeO2(:,:)
          betaCH4_soil(:,il,:) =  &
               ( xcCH4_soil(:,il+1,:)*CH4_soil(:,il+1,:) +  &
               xdCH4_soil(:,il+1,:)*betaCH4_soil(:,il+1,:) ) / xeCH4(:,:)
          !
       ENDWHERE
       !
    ENDDO
    !
    ! 4. store thickness of the different levels for all soil types (for security) 
    !
    zf_coeff_snow(:,:,:) = zf_snow(:,:,:)
    zi_coeff_snow(:,:,:) = zi_snow(:,:,:)


  END SUBROUTINE soil_gasdiff_coeff
 
!!
!================================================================================================================================
!! SUBROUTINE   : soil_gasdiff_diff
!!
!>\BRIEF        This routine update oxygen and methane in the snow and soil 
!!
!! DESCRIPTION : 
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!! 
!! REFERENCE(S) : None
!!
!! FLOWCHART11    : None
!! \n
!_
!================================================================================================================================    
  
  SUBROUTINE soil_gasdiff_diff( kjpindex,time_step,index,pb,tsurf, O2_snow, CH4_snow, O2_soil, CH4_soil)
   
  !! 0. Variable and parameter declaration

    !! 0.1  Input variables

    INTEGER(i_std), INTENT(in)                          :: kjpindex             !! number of grid points 
    REAL(r_std), INTENT(in)                             :: time_step            !! time step in seconds
    REAL(r_std), DIMENSION(:), INTENT(in)               :: pb                   !! Surface pressure
    REAL(r_std), DIMENSION(:), INTENT(in)               :: tsurf                !! Surface temperature
    INTEGER(i_std),DIMENSION(:),INTENT(in)              :: index                !! Indeces of the points on the map 
    !! 0.2  Output variables

    !! 0.3  Modified variables

    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)        :: O2_snow              !! oxygen (g O2/m**3 air)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)        :: CH4_snow             !! methane (g CH4/m**3 air)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)        :: O2_soil              !! oxygen (g O2/m**3 air)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)        :: CH4_soil             !! methane (g CH4/m**3 air)

    !! 0.4 local variables
 
    INTEGER(i_std)                                             :: it, ip, il, iv
    LOGICAL, DIMENSION(kjpindex,nvm)                           :: snowtop
    REAL(r_std), DIMENSION(kjpindex,nvm)                       :: O2sa, CH4sa
    
    !
    ! 1.1 Determine which is the first existing soil type.
    !
    snowtop(:,:) = .FALSE.  
    !
    !ignore snow for now...
    WHERE ( heights_snow(:,:) .GT. hmin_tcalc )
       snowtop(:,:) = .TRUE.
    ENDWHERE
    !
    ! 2.gas diffusion
    !
    ! 2.1 top level
    !
    ! 2.1.1 non-existing
    !
    DO iv = 1, nvm
       O2sa(:,iv) = pb(:)/(RR*tsurf(:)) * O2_surf * wO2
       CH4sa(:,iv) = pb(:)/(RR*tsurf(:)) * CH4_surf * wCH4
    ENDDO
    !
    WHERE ( (.NOT. snowtop(:,:)) .AND. veget_mask_2d(:,:) ) ! it equals 1 (snow) but there is no snow...
       !
       O2_snow(:,1,:) = O2sa(:,:)
       CH4_snow(:,1,:) = CH4sa(:,:)
       !
       O2_soil(:,1,:) = ( O2sa(:,:) + mu_soil*betaO2_soil(:,1,:) ) / &
            ( 1. + mu_soil*(1.-alphaO2_soil(:,1,:)) )
       CH4_soil(:,1,:) = ( CH4sa(:,:) + mu_soil*betaCH4_soil(:,1,:) ) / &
            ( 1. + mu_soil*(1.-alphaCH4_soil(:,1,:)) )
       !
    ENDWHERE
    !
    ! 2.1.2 first existing soil type
    !
    WHERE ( snowtop(:,:) .AND. veget_mask_2d(:,:) )
       !
       O2_snow(:,1,:) = ( O2sa(:,:) + mu_snow(:,:)*betaO2_snow(:,1,:) ) / &
            ( 1. + mu_snow(:,:)*(1.-alphaO2_snow(:,1,:)) )
       CH4_snow(:,1,:) = ( CH4sa(:,:) + mu_snow(:,:)*betaCH4_snow(:,1,:) ) / &
            ( 1. + mu_snow(:,:)*(1.-alphaCH4_snow(:,1,:)) )
       !
       O2_soil(:,1,:) =  &
            alphaO2_snow(:,nsnow,:) * O2_snow(:,nsnow,:) + &
            betaO2_snow(:,nsnow,:)
       CH4_soil(:,1,:) =  &
            alphaCH4_snow(:,nsnow,:) * CH4_snow(:,nsnow,:) + &
            betaCH4_snow(:,nsnow,:)
       ! debug: need to check for weird numbers here!
    ENDWHERE
    !
    ! 2.2 other levels
    !
    DO il = 2, nsnow
       
       WHERE ( veget_mask_2d(:,:) )
          !
          O2_snow(:,il,:) =  &
               alphaO2_snow(:,il-1,:) * O2_snow(:,il-1,:) + &
               betaO2_snow(:,il-1,:)
          CH4_snow(:,il,:) =  &
               alphaCH4_snow(:,il-1,:) * CH4_snow(:,il-1,:) + &
               betaCH4_snow(:,il-1,:)
       END WHERE
    ENDDO
    DO il = 2, ngrnd
       
       WHERE ( veget_mask_2d(:,:)  )
          !
          O2_soil(:,il,:) =  &
               alphaO2_soil(:,il-1,:) * O2_soil(:,il-1,:) + &
               betaO2_soil(:,il-1,:)
          CH4_soil(:,il,:) =  &
               alphaCH4_soil(:,il-1,:) * CH4_soil(:,il-1,:) + &
               betaCH4_soil(:,il-1,:)
       END WHERE
    ENDDO

  END SUBROUTINE soil_gasdiff_diff
 
!!
!================================================================================================================================
!! SUBROUTINE   : get_gasdiff
!!
!>\BRIEF        This routine update oxygen and methane in the snow and soil 
!!
!! DESCRIPTION : 
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!! 
!! REFERENCE(S) : None
!!
!! FLOWCHART11    : None
!! \n
!_
!================================================================================================================================    
  SUBROUTINE get_gasdiff (kjpindex,hslong,tprof,snow,airvol_snow, &
       totporO2_snow,totporCH4_snow,diffO2_snow,diffCH4_snow, &
       airvol_soil,totporO2_soil,totporCH4_soil,diffO2_soil,diffCH4_soil, snowrho)
   
  !! 0. Variable and parameter declaration

    !! 0.1  Input variables

    INTEGER(i_std), INTENT(in)                    :: kjpindex          !! number of grid points 
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)     :: hslong            !! deep long term soil humidity profile 
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)     :: tprof             !! Soil temperature (K)      
    REAL(r_std), DIMENSION(:,:), INTENT(in)       :: snowrho           !! snow density (Kg/m^3)
    REAL(r_std), DIMENSION(:),  INTENT (in)       :: snow              !! Snow mass [Kg/m^2]
    !! 0.2  Output variables

    REAL(r_std), DIMENSION(:,:,:), INTENT(out)    :: airvol_soil
    REAL(r_std), DIMENSION(:,:,:), INTENT(out)    :: totporO2_soil     !! total O2 porosity (Tans, 1998)
    REAL(r_std), DIMENSION(:,:,:), INTENT(out)    :: totporCH4_soil    !! total CH4 porosity 
    REAL(r_std), DIMENSION(:,:,:), INTENT(out)    :: diffO2_soil       !! oxygen diffusivity (m**2/s)
    REAL(r_std), DIMENSION(:,:,:), INTENT(out)    :: diffCH4_soil      !! methane diffusivity (m**2/s)
    REAL(r_std), DIMENSION(:,:,:), INTENT(out)    :: airvol_snow
    REAL(r_std), DIMENSION(:,:,:), INTENT(out)    :: totporO2_snow     !! total O2 porosity (Tans, 1998)
    REAL(r_std), DIMENSION(:,:,:), INTENT(out)    :: totporCH4_snow    !! total CH4 porosity (Tans, 1998)
    REAL(r_std), DIMENSION(:,:,:), INTENT(out)    :: diffO2_snow       !! oxygen diffusivity (m**2/s)
    REAL(r_std), DIMENSION(:,:,:), INTENT(out)    :: diffCH4_snow      !! methane diffusivity (m**2/s)
   
    !! 0.3  Modified variables

    !! 0.4 local variables
 
    REAL(r_std), DIMENSION(kjpindex,nsnow,nvm)                 :: density_snow
    REAL(r_std), DIMENSION(kjpindex,nsnow,nvm)                 :: porosity_snow
    REAL(r_std), DIMENSION(kjpindex,nsnow,nvm)                 :: tortuosity_snow
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm)                 :: density_soil
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm)                 :: porosity_soil
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm)                 :: tortuosity_soil
    INTEGER(i_std)                                             :: it,ip, il, iv
    REAL(r_std)                                                :: x, rho_iw
    REAL(r_std)                                                :: csat, fng
    REAL(r_std),  SAVE                                         :: cond_fact 
!$OMP THREADPRIVATE(cond_fact)
    LOGICAL, SAVE                                              :: firstcall_get_gasdiff=.TRUE.
!$OMP THREADPRIVATE(firstcall_get_gasdiff)
    

    IF (firstcall_get_gasdiff) THEN
       cond_fact=1.
       CALL getin_p('COND_FACT',cond_fact) 
       WRITE(*,*) 'COND_FACT=',cond_fact
       firstcall_get_gasdiff=.FALSE.
    ENDIF
    
    !
    ! 1. Three-layers snow model with snow density resolved at each snow layer
    !
    DO iv = 1, nvm 
       density_snow(:,:,iv) = snowrho(:,:)
    ENDDO
    porosity_snow(:,:,:) = (1. - density_snow(:,:,:)/rho_ice )
    tortuosity_snow(:,:,:) = porosity_snow(:,:,:)**(1./3.)     ! based on Sommerfeld et al., GBC, 1996
    diffO2_snow(:,:,:) = diffO2_air * porosity_snow(:,:,:) * tortuosity_snow(:,:,:)
    diffCH4_snow(:,:,:) = diffCH4_air * porosity_snow(:,:,:) * tortuosity_snow(:,:,:)
    airvol_snow(:,:,:) = MAX(porosity_snow(:,:,:),avm)
    totporO2_snow(:,:,:) = airvol_snow(:,:,:)
    totporCH4_snow(:,:,:) = airvol_snow(:,:,:)
    !
    ! 2. soil: depends on temperature and soil humidity
    !
    DO ip = 1, kjpindex
       !
       DO iv = 1, nvm
          !
          IF ( veget_mask_2d(ip,iv) ) THEN
             !
             DO il = 1, ngrnd
                !
                ! 2.1 soil dry density, porosity, and dry heat capacity
                !
                porosity_soil(ip,il,iv) = tetasat
                !
                !
                ! 2.2 heat capacity and density as a function of
                !     ice and water content
                !  removed these as we are calculating thermal evolution in the sechiba subroutines
               
                !
                ! 2.3 oxygen diffusivity: soil can get waterlogged,
                !     therefore take soil humidity into account
                !
                tortuosity_soil(ip,il,iv) = 2./3. ! Hillel, 1980
                airvol_soil(ip,il,iv) = porosity_soil(ip,il,iv)*(1.-hslong(ip,il,iv))  
                totporO2_soil(ip,il,iv) = airvol_soil(ip,il,iv) + porosity_soil(ip,il,iv)*BunsenO2*hslong(ip,il,iv)  
                totporCH4_soil(ip,il,iv) = airvol_soil(ip,il,iv) + porosity_soil(ip,il,iv)*BunsenCH4*hslong(ip,il,iv)  
                diffO2_soil(ip,il,iv) = (diffO2_air*airvol_soil(ip,il,iv) + & 
                     diffO2_w*BunsenO2*hslong(ip,il,iv)*porosity_soil(ip,il,iv))*tortuosity_soil(ip,il,iv)  
                diffCH4_soil(ip,il,iv) = (diffCH4_air*airvol_soil(ip,il,iv) + & 
                     diffCH4_w*BunsenCH4*hslong(ip,il,iv)*porosity_soil(ip,il,iv))*tortuosity_soil(ip,il,iv)  
                !
          END DO
       ELSE
          tortuosity_soil(ip,:,iv) = EPSILON(0.)
          airvol_soil(ip,:,iv) =  EPSILON(0.)
          totporO2_soil(ip,:,iv) =  EPSILON(0.)
          totporCH4_soil(ip,:,iv) = EPSILON(0.)
          diffO2_soil(ip,:,iv) = EPSILON(0.)
          diffCH4_soil(ip,:,iv) =  EPSILON(0.)
       END IF
    ENDDO
 ENDDO

END SUBROUTINE get_gasdiff
  
!!
!================================================================================================================================
!! SUBROUTINE   : traMplan
!!
!>\BRIEF        This routine calculates plant-mediated transport of methane
!!
!! DESCRIPTION : 
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!! 
!! REFERENCE(S) : None
!!
!! FLOWCHART11    : None
!! \n
!_
!================================================================================================================================    
  SUBROUTINE traMplan(CH4,O2,kjpindex,time_step,totporCH4,totporO2,z_root,rootlev,Tgr,Tref,hslong,flupmt, &
       refdep, zi_soil, tprof)

  !! 0. Variable and parameter declaration

    !! 0.1  Input variables
    
    INTEGER(i_std), INTENT(in)			     :: kjpindex    
    REAL(r_std), INTENT(in)                          :: time_step      !! time step in seconds
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)        :: totporO2       !! total oxygen porosity
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)        :: totporCH4      !! total methane porosity
    REAL(r_std), DIMENSION(:,:,:),INTENT(in)         :: tprof          !! soil temperature (K)
    INTEGER(i_std),DIMENSION(:,:),INTENT(in)	     :: rootlev        !! the deepest model level within the rooting depth 
    REAL(r_std), DIMENSION(:,:),INTENT(in)	     :: z_root         !! the rooting depth
    REAL(r_std), INTENT(in) 			     :: Tgr            !! Temperature at which plants begin to grow (C)
    REAL(r_std), DIMENSION(ngrnd), INTENT(in)        :: zi_soil        !!  depths at intermediate levels 
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)        :: hslong         !! deep soil humidity

    !! 0.2 Output variables

    REAL(r_std), DIMENSION(:,:), INTENT(out)         :: flupmt         !! plant-mediated methane flux (g m-2 s-1)

    !! 0.3 Modified variables

    REAL(r_std), DIMENSION(:,:),INTENT(inout)        :: Tref           !! Ref. temperature for growing season caluculation (C)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)     :: O2
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)     :: CH4

    !! 0.4 local variables
    REAL(r_std), DIMENSION(kjpindex,nvm)   			  :: CH4atm         !! CH4 atm concentration
    REAL(r_std), DIMENSION(kjpindex,nvm)                          :: dCH4           !! delta CH4 per m3 air
    REAL(r_std), DIMENSION(kjpindex,nvm)   			  :: dO2            !! O2 change
    REAL(r_std), DIMENSION(kjpindex,nvm)	       		  :: fgrow          !! Plant growing state (maturity index)
    REAL(r_std)                            			  :: froot          !! vertical distribution of roots
    REAL(r_std)					                  :: Tmat           !! Temperature at which plants reach maturity (C)
    REAL(r_std), PARAMETER 				          :: La_min = zero
    REAL(r_std), PARAMETER                                        :: La = 4.
    REAL(r_std), PARAMETER 				          :: La_max = La_min + La
    REAL(r_std), PARAMETER 				          :: Tveg = 10      !! Vegetation type control on the plant-mediated transport, Adjustable parameter,
                                                                                    !! but we start from 10 following Walter et al (2001) tundra value 
    REAL(r_std), PARAMETER 				          :: Pox = 0.5      !! fraction of methane oxydized near the roots
    LOGICAL, SAVE 				                  :: firstcall=.TRUE.
!$OMP THREADPRIVATE(firstcall)
    INTEGER(i_std)       				          :: il,ip, iv
    REAL(r_std), INTENT(in)			                  :: refdep         !! Depth to compute reference temperature for the growing season (m)
    INTEGER(i_std), SAVE				          :: reflev = 0     !! Level closest to reference depth refdep
!$OMP THREADPRIVATE(reflev)
   

    IF (firstcall) THEN
       firstcall = .FALSE.

       ! Find the level closest to refdep
       DO il=1,ngrnd
          IF (zi_soil(il) .GT. refdep .AND. reflev.EQ.0) reflev = il-1
       ENDDO
       IF (reflev.EQ.0) reflev = ngrnd

    ENDIF

     ! Update seasonal reference temperature trace record  
     WHERE ( veget_mask_2d(:,:) )          
        Tref(:,:) = tprof(:,reflev,:) - ZeroCelsius
     END WHERE

    Tmat = Tgr + 10._r_std
    flupmt(:,:) = zero
    CH4atm(:,:) = zero  


    ! Plant growing state (maturity index)
    WHERE (Tref(:,:).LE.Tgr .AND. veget_mask_2d(:,:) )
       fgrow(:,:) = La_min
    ELSEWHERE (Tref(:,:).GE.Tmat .AND. veget_mask_2d(:,:) )
       fgrow(:,:) = La_max
    ELSEWHERE   ( veget_mask_2d(:,:))
       fgrow(:,:) = La_min + La * (1 - ((Tmat - Tref(:,:))/(Tmat - Tgr))**2)
    ENDWHERE

    DO ip=1,kjpindex
       DO iv = 1, nvm
          IF ( (z_root(ip,iv) .GT. 0.) .AND. veget_mask_2d(ip,iv) ) THEN ! added this to prevent pmt calcs when soil frozen
             DO il=1,rootlev(ip,iv)
                ! vertical distribution of roots
                froot = MAX( 2 * (z_root(ip,iv) - REAL( zi_soil(il) )) / z_root(ip,iv), zero) 
                ! Methane removal from a given depth. We assume that the methane
                ! in air pores is always in equilibrium with that dissolved 
                ! in water-filled pores. If soil humidity is low,
                ! with root water as well
                ! We assume that PMT is proportional to soil humidity
                dCH4(ip,iv) = 0.01_r_std * Tveg * froot * fgrow(ip,iv) * hslong(ip,il,iv) * (CH4(ip,il,iv) - CH4atm(ip,iv))  
                ! No transport if soil concentration is less than atmospheric
                IF (dCH4(ip,iv).LT.CH4atm(ip,iv)) dCH4(ip,iv) = zero 
                ! Strange thing in WH 2001: 0.01*Tveg*froot*fgrow > 1
                ! at Tveg=15, froot&fgrow=max, i.e. more CH4 is taken than available
                ! So need to impose a limitation:
                IF (dCH4(ip,iv).GT.CH4(ip,il,iv)) dCH4(ip,iv) = CH4(ip,il,iv)
                ! Methane concentration is decreased within the root layer:
                
                CH4(ip,il,iv) = CH4(ip,il,iv) - dCH4(ip,iv)
                ! O2 concentration is decreased in reaction with
                ! dCH4*Pox*time_step
                dO2(ip,iv) = dCH4(ip,iv)*Pox * wO2/wCH4 * totporCH4(ip,il,iv)/totporO2(ip,il,iv)
                IF ( dO2(ip,iv).LT.O2(ip,il,iv) ) O2(ip,il,iv) = O2(ip,il,iv) - dO2(ip,iv)
                
                ! CO2 concentration is increased by dCH4(:)*Pox
                
                ! Integration	 
                flupmt(ip,iv) = flupmt(ip,iv) + dCH4(ip,iv)*totporCH4(ip,il,iv)/time_step * (1 - Pox) * &
                     ( zf_soil(il) - zf_soil(il-1) )
             ENDDO
          END IF
       ENDDO
    ENDDO
    
    
  END SUBROUTINE traMplan
  
!!
!================================================================================================================================
!! SUBROUTINE   : ebullition
!!
!>\BRIEF        This routine calculates CH4 ebullition
!!
!! DESCRIPTION : 
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!! 
!! REFERENCE(S) : None
!!
!! FLOWCHART11    : None
!! \n
!_
!================================================================================================================================     
  SUBROUTINE ebullition (kjpindex,time_step,tprof,totporCH4_soil,hslong,Ch4_soil,febul)
   
  !! 0. Variable and parameter declaration

    !! 0.1  Input variables 

    INTEGER(i_std), INTENT(in)                     :: kjpindex
    REAL(r_std), INTENT(in)                        :: time_step      !! time step in seconds
    REAL(r_std), DIMENSION(:,:,:),INTENT(in)       :: tprof          !! soil temperature (K)
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)      :: totporCH4_soil !! total methane porosity
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)      :: hslong         !! deep soil humidity

    !! 0.2 Output variables
    
    REAL(r_std), DIMENSION(:,:), INTENT(out)       :: febul          !! CH4 ebullition

    !! 0.3 Modified variables

    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)   :: Ch4_soil       !! methane 

    !! 0.4 Local variables
    REAL(r_std)                                                 :: dCH4, CH4d
    INTEGER(i_std)                                              :: ip, il, iv
    REAL(r_std)                                                 :: dz
    REAL(r_std), PARAMETER                                      :: tortuosity=2./3. 
    REAL(r_std), PARAMETER                                      :: wsize=0.01
    REAL(r_std), PARAMETER                                      :: CH4wm = 12.E-3 !! CH4 concentration threshold for ebullition (8-16 mg/m3 in Walter&Heimann 2000)
    REAL(r_std)                                                 :: hum

    febul(:,:)=zero

    DO ip=1,kjpindex
       DO iv = 1, nvm
          IF ( veget_mask_2d(ip,iv) ) THEN
             IF (hslong(ip,1,iv).GT.ebuthr) THEN
                DO il = ngrnd, 1, -1
                   CH4d = Ch4_soil(ip,il,iv) - CH4wm/BunsenCH4
                   IF (CH4d .GT. EPSILON(0.)) THEN  
                      IF (il.GT.1) THEN
                         dz = zi_soil(il) - zi_soil(il-1)
                         hum = ( hslong(ip,il,iv) + hslong(ip,il-1,iv) ) / 2 
                      ELSE
                         dz = zi_soil(1)
                         hum = hslong(ip,1,iv)
                      ENDIF
                      
                      dCH4 = hum**( dz/wsize/tortuosity ) * CH4d
                      dCH4 = CH4d 
                      
                      Ch4_soil(ip,il,iv) = Ch4_soil(ip,il,iv) - dCH4
                      
                   
                      febul(ip,iv) = febul(ip,iv) + dCH4 * totporCH4_soil(ip,il,iv) *  &
                           ( zf_soil(il) - zf_soil(il-1) ) / time_step
                      
                   ENDIF 
                ENDDO
             ENDIF 
          END IF
       ENDDO
    ENDDO
  END SUBROUTINE ebullition
  
!!
!================================================================================================================================
!! FUNCTION   : stomate_soil_carbon_discretization_microactem
!!
!>\BRIEF        This function calculates parameters describing bacterial activity (time constant tau[s]) as a function of temperature
!!
!! DESCRIPTION : mc_peat, soiltile is added by xw following the function microactem in mict peat
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!! 
!! REFERENCE(S) : None
!!
!! FLOWCHART11    : None
!! \n
!_
!================================================================================================================================     
  FUNCTION stomate_soil_carbon_discretization_microactem &
       ( temp, frozen_respiration_func, moist_in, &
       i_ind, j_ind, k_ind, zi_soil, mc_peat, soiltile  ) RESULT ( fbact )

    ! tdeep_celsius, frozen_respiration_func, hsdeep, kjpindex, ngrnd, nvm, znt,&
    !       mc_peat, soiltile
    !! 0. Variable and parameter declaration

    !! 0.1  Input variables 
    
    INTEGER(i_std), INTENT(in)                        :: i_ind       !! kjpindex
    INTEGER(i_std), INTENT(in)                        :: j_ind       !! ngrnd
    INTEGER(i_std), INTENT(in)                        :: k_ind       !! nvm
    INTEGER(i_std), INTENT(in)                        :: frozen_respiration_func !! Method for soil decomposition function
    REAL, DIMENSION(i_ind, j_ind, k_ind), INTENT(in)  :: moist_in    !! deep long term soil humidity profile (unitless)
    REAL, DIMENSION(i_ind, j_ind, k_ind), INTENT(in)  :: temp        !! deep temperature profile celsius (C)
    REAL, DIMENSION(j_ind), INTENT(in)                :: zi_soil     !! soil layers

    !! 0.2 Output variables 

    !! 0.3 Modified variables

    !! 0.4 Local variables

    REAL, DIMENSION(i_ind, j_ind, k_ind)              :: fbact            !!   turnover constant (day)
    REAL, DIMENSION(i_ind, j_ind, k_ind)              :: tempfunc_result  !! temperature profile in C
    REAL, DIMENSION(i_ind, j_ind, k_ind)              :: temp_kelvin      !! tmperature profile in K
    INTEGER                                           :: ii, ij, ik
    REAL, DIMENSION(i_ind, j_ind, k_ind)              :: moistfunc_result !! moisture  profile 
    logical, parameter                                :: limit_decomp_moisture = .true.  !! to be removed later.

    !peat
    REAL(r_std), DIMENSION (i_ind,j_ind), INTENT(in)  :: mc_peat    !! soil moisture for peat 
    REAL(r_std),DIMENSION (:,:), INTENT(in)           :: soiltile   !! soil tile fraction
    REAL, DIMENSION(j_ind,k_ind)                      :: peat_tau   !! peat turn over time
    REAL, DIMENSION(i_ind, j_ind, k_ind)              :: moistfunc_result_peat !! moisture control function for soil heteo respiration over peat
    REAL(r_std),DIMENSION(dim_moyanno_moist)          :: mc         !! relative humidity used for Moyano et al., 2012, with 0.02 interval used    
    REAL(r_std),DIMENSION(dim_moyanno_moist)          :: pcsr       !! response of soil heteo respiration to soil moisture 
    REAL(r_std),DIMENSION(dim_moyanno_moist)          :: sr         !! normalized response  of soil heteo respiration to soil moisture by using the max
    REAL(r_std),DIMENSION(dim_moyanno_moist)          :: corgmat    !! normalizedresponse of soil heteo respiration to soil moisture by using the min-max
    INTEGER(i_std)                                    :: ind        !! the index for the max of sr 
    INTEGER(i_std)                                    :: mc_ind     !! index of orc modeled soil moisture corresponding to mc array used for moyano
                                                                    !! when peat activated
    INTEGER(i_std)                                    :: agri_mc_ind !! index of orc modeled soil moisture corresponding to mc array used for moyano
                                                                    !! when agri_peat activated
    INTEGER(i_std)                                    :: istm
    
    temp_kelvin(:,:,:) = temp(:,:,:) + ZeroCelsius
    SELECT CASE(frozen_respiration_func)

    CASE(0) ! this is the standard ORCHIDEE state

       tempfunc_result(:,:,:) = EXP( soil_Q10 * ( temp_kelvin(:,:,:) - (ZeroCelsius+30.) ) / 10. )
       tempfunc_result(:,:,:) = MIN( 1._r_std, tempfunc_result(:,:,:) )

    CASE(1)  ! cutoff respiration when T < -1C
       WHERE (temp_kelvin(:,:,:) .GT. ZeroCelsius ) ! normal as above
          tempfunc_result(:,:,:) = EXP( soil_Q10 * ( temp_kelvin(:,:,:) - (ZeroCelsius+30.) ) / 10. )
       ELSEWHERE (temp_kelvin(:,:,:) .GT. ZeroCelsius - 1. )  ! linear dropoff to zero
          tempfunc_result(:,:,:) = (temp_kelvin(:,:,:) - (ZeroCelsius - 1.)) * &
               EXP( soil_Q10 * ( ZeroCelsius - (ZeroCelsius+30.) ) / 10. )
       ELSEWHERE  ! zero
          tempfunc_result(:,:,:) = EPSILON(0.)
       endwhere

       tempfunc_result(:,:,:) = MAX(MIN( 1._r_std, tempfunc_result(:,:,:) ), EPSILON(0.))

    CASE(2)  ! cutoff respiration when T < -3C
       WHERE (temp_kelvin(:,:,:) .GT. ZeroCelsius ) ! normal as above
          tempfunc_result(:,:,:) = EXP( soil_Q10 * ( temp_kelvin(:,:,:) - (ZeroCelsius+30.) ) / 10. )
       ELSEWHERE (temp_kelvin(:,:,:) .GT. ZeroCelsius - 3. )  ! linear dropoff to zero
          tempfunc_result(:,:,:) = ((temp_kelvin(:,:,:) - (ZeroCelsius - 3.))/3.) * &
               EXP( soil_Q10 * ( ZeroCelsius - (ZeroCelsius+30.) ) / 10. )
       ELSEWHERE  ! zero
          tempfunc_result(:,:,:) = EPSILON(0.)
       endwhere

    CASE(3)  ! q10 = 100 when below zero
       WHERE (temp_kelvin(:,:,:) .GT. ZeroCelsius ) ! normal as above
          tempfunc_result(:,:,:) = EXP( soil_Q10 * ( temp_kelvin(:,:,:) - (ZeroCelsius+30.) ) / 10. )
       ELSEWHERE 
          tempfunc_result(:,:,:) = EXP( log(100.) * ( temp_kelvin(:,:,:) - (ZeroCelsius) ) / 10. ) * &
               EXP( soil_Q10 * ( -30. ) / 10. )
       endwhere

    CASE(4)  ! q10 = 1000 when below zero
       WHERE (temp_kelvin(:,:,:) .GT. ZeroCelsius ) ! normal as above
          tempfunc_result(:,:,:) = EXP( soil_Q10 * ( temp_kelvin(:,:,:) - (ZeroCelsius+30.) ) / 10. )
       ELSEWHERE 
          tempfunc_result(:,:,:) = EXP( log(1000.) * ( temp_kelvin(:,:,:) - (ZeroCelsius) ) / 10. ) * &
               EXP( soil_Q10 * ( -30. ) / 10. )
       endwhere

    CASE DEFAULT
       CALL ipslerr_p(3,'stomate_soil_carbon_discretization_microactem', &
            'frozen_respiration_func can only be 0, 1, 2, 3 or 4','','')
    END SELECT
    tempfunc_result(:,:,:) = MAX(MIN( 1._r_std, tempfunc_result(:,:,:) ), EPSILON(0.))

    !---- stomate residence times: -----!
    ! residence times in carbon pools (days)
    !carbon_tau(iactive) = .149 * one_year        !!!!???? 1.5 years
    !carbon_tau(islow) = 5.48 * one_year          !!!!???? 25 years
    !carbon_tau(ipassive) = 241. * one_year       !!!!???? 1000 years
    !-----------------------------------!
    IF ( limit_decomp_moisture ) THEN
       ! stomate moisture control function
       moistfunc_result(:,:,:) = -1.1 * moist_in(:,:,:) * moist_in(:,:,:) + 2.4 * moist_in(:,:,:) - 0.29
       moistfunc_result(:,:,:) = max( 0.25_r_std, min( 1._r_std, moistfunc_result(:,:,:) ) )
    ELSE
       moistfunc_result(:,:,:) = 1._r_std
    ENDIF

    DO ij = 1, ngrnd
      fbact(:,ij,:) = stomate_tau/(moistfunc_result(:,ij,:) * tempfunc_result(:,ij,:)) / EXP(-zi_soil(ij)/depth_modifier)
    ENDDO

    ! chaoyue: We tentatively increase the turnover of soil C in croplands,
    ! as shown here to decrease its tau -- the residence time.
    DO ik = 1,nvm
        fbact(:,:,ik) = fbact(:,:,ik)/ decomp_factor(ik)
    ENDDO

    !!! Moyano et al., 2012,for organic soils
    !! volumetric moisture, 0.02 interval
    !! mc below is between 0.01 et 0.89 if dim_moyano_moist is 45
    IF (perma_peat) THEN
       DO ii=1,dim_moyanno_moist
          mc(ii)=0.01+0.02*(ii-1)
       ENDDO
    ENDIF
    
    !!calculate pcsr according to equation2 in  Moyano et al., 2012
    !!for orgainc soil, bd=1.2 g/cm3, clay=0.3 fraction, organic carbon 0.05 g/g
    ! proportional response of soil respiration to soil moisture
    IF (perma_peat) THEN
       pcsr(:)=0.97509-0.48212*mc(:)+1.83997*(mc(:)**2)-1.56379*(mc(:)**3)+ &
               0.09867*1.2+1.39944*0.05+0.17938*0.3-0.30307*mc(:)*1.2-0.30885*mc(:)*0.3
       ! Moyano equation 3: 20250902
       !pcsr(:)=1.134-0.67*mc(:) + 1.08*(mc(:)**2)-0.5*(mc(:)**3)
    ENDIF  
    
    !!relative respiration
    IF (perma_peat) THEN
       DO ii=1,dim_moyanno_moist
          IF (ii==1) THEN
             sr(ii) = pcsr(ii)
          ELSE
             sr(ii)=sr(ii-1)* pcsr(ii) 
          ENDIF
       ENDDO
       sr(:)=sr(:)/ MAXVAL(sr)
    ENDIF

     !!!rescaling respiration from 0 to 1 in the range of 0 to optimum
    IF (perma_peat ) THEN
        corgmat(:)=sr(:)
        ind= MAXLOC(corgmat,1)  
           corgmat(1:ind)=corgmat(1:ind)-MINVAL(corgmat(1:ind))
           corgmat(1:ind)=corgmat(1:ind)/MAXVAL(corgmat(1:ind))
    ENDIF

    !!! find corgmat value corresponding to current volumetric moisture
    IF (perma_peat ) THEN
       DO ii=1,i_ind
          DO ij=1,j_ind
             DO ik=1,k_ind
                istm = pref_soil_veg(ik)
                IF (is_peat(ik) .AND. soiltile(ii,istm) .GT. zero ) THEN
                   ! To be verified 202603 xw
                   mc_ind = MIN(dim_moyanno_moist, MAX(1, INT(mc_peat(ii,ij)/0.02)+1))
                   moistfunc_result_peat(ii,ij,ik)= corgmat(mc_ind) 
                   moistfunc_result_peat(ii,ij,ik)= MIN(un,MAX(EPSILON(0.),moistfunc_result_peat(ii,ij,ik)))
                ELSE
                   moistfunc_result_peat(ii,ij,ik) = moistfunc_result(ii,ij,ik)                      
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    
    IF (agri_peat) THEN 
       CALL ipslerr_p(3,'stomate_soil_carbon_discretization_microactem', &
            'agri_peat is activated ','the below code should be modified firstly by using is_croppeat True', &
            'instead of using hard-coded pft index.')
       DO ii=1,i_ind
          DO ij=1,j_ind
             DO ik=1,k_ind
                istm = pref_soil_veg(ik)
                ! from mict-peat: pft16 wheat on peatland; pft1 natural peat
                ! in the trunk, we have not decided how many different agricol vegetations will be used. 
                ! To verify if we should use the 45 or dim_moyano_moist, xw 202603
                IF (ik==15 .OR. ik==16 .AND. soiltile(ii,istm) .GT. zero) THEN  !!crops on peatland
                   agri_mc_ind = MIN(dim_moyanno_moist, MAX(1, INT(mc_peat(ii,ij)/0.02)+1))
                   moistfunc_result(ii,ij,ik)= corgmat(agri_mc_ind)
                   moistfunc_result(ii,ij,ik)= MIN(un,MAX(EPSILON(0.),moistfunc_result(ii,ij,ik)))
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    
    !!! peat turnover time increase with depth
    ! 
    peat_tau(:,:) = 0.0
    IF (perma_peat) THEN
       DO ik=1,k_ind
          DO ij=1, j_ind
             IF (is_peat(ik)) THEN
                IF (ij .LE. 12) THEN
                   peat_tau(ij,ik)= tau_peat*EXP(zi_soil(ij)/z_tau)
                ELSE
                   peat_tau(ij,ik)= tau_peat*EXP(zi_soil(12)/z_tau)
                ENDIF
             ENDIF
          ENDDO
       ENDDO
       
       DO ii=1,i_ind
          DO ij=1,j_ind
             DO ik=1,k_ind
                istm = pref_soil_veg(ik)
                IF (is_peat(ik) .AND. soiltile(ii,istm) .GT. zero) THEN
                   fbact(ii,ij,ik)= peat_tau(ij,ik)/(moistfunc_result_peat(ii,ij,ik)* tempfunc_result(ii,ij,ik))    
                   fbact(ii,ij,ik) = fbact(ii,ij,ik)/ decomp_factor(ik)
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    
  END FUNCTION stomate_soil_carbon_discretization_microactem
  
  
!!
!================================================================================================================================
!! SUBROUTINE   : snowlevels
!!
!>\BRIEF        This routine calculates depths of full levels and intermediate
!!              levels related to snow pack
!!
!! DESCRIPTION : 
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!! 
!! REFERENCE(S) : None
!!
!! FLOWCHART11    : None
!! \n
!_
!================================================================================================================================     
 
  SUBROUTINE snowlevels( kjpindex, snowdz, zi_snow, zf_snow, veget_max )

  !! 0. Variable and parameter declaration

    !! 0.1  Input variables     

    INTEGER(i_std), INTENT(in)                               :: kjpindex
    REAL(r_std), DIMENSION(:,:),INTENT(in)                   :: veget_max     !! maximum vegetation fraction
    REAL(r_std), DIMENSION(:,:),INTENT(in)                   :: snowdz        !! snow depth [m]

    !! 0.2 Output variables

    !! 0.3 Modified variables

    REAL(r_std), DIMENSION(:,0:,:), INTENT(inout)            :: zf_snow       !! depths of full levels (m)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)             :: zi_snow       !! depths of intermediate levels (m)

    !! 0.4 Local variables

    REAL(r_std), DIMENSION(kjpindex,nvm)                                :: z_alpha       !! parameter of the geometric series
    INTEGER(i_std)                                                      :: il,it, ix, iv
    INTEGER(i_std)                                                      :: it_beg,it_end
    INTEGER(i_std), PARAMETER                                           :: niter = 10
    REAL(r_std), DIMENSION(kjpindex)                                    :: dxmin
    INTEGER(i_std), DIMENSION(kjpindex)                                 :: imin
    INTEGER(i_std)                                                      :: i,j
    REAL(r_std), DIMENSION(kjpindex,nvm)                                :: xi, xf
    REAL(r_std), DIMENSION(kjpindex,nsnow,nvm)                          :: snowdz_pft 

    snowdz_pft(:,:,:) = 0.0 
    DO il = 1,nsnow 
       DO iv = 1, nvm
          WHERE ( veget_mask_2d(:,iv) ) 
                snowdz_pft(:,il,iv) = snowdz(:,il)
          ENDWHERE
       ENDDO 
    ENDDO
    !
    ! calculate snow discretisation
    !
    WHERE ( veget_mask_2d(:,:) )
       zf_snow(:,0,:) = 0.
    END WHERE
    !
    DO il = 1, nsnow
       IF ( il .EQ. 1 ) THEN
          WHERE ( veget_mask_2d(:,:) )
             
             zi_snow(:,il,:) = snowdz_pft(:,1,:) / 2. 
             
             zf_snow(:,il,:) = snowdz_pft(:,1,:)
        
          END WHERE
       ENDIF

       IF ( il .GT. 1 ) THEN
          WHERE ( veget_mask_2d(:,:) )
             
             zi_snow(:,il,:) = zf_snow(:,il-1,:) + snowdz_pft(:,il,:) / 2 
             
             zf_snow(:,il,:) = SUM(snowdz_pft(:,1:il,:),2)

          END WHERE
       ENDIF

    ENDDO
    
    DO ix = 1, kjpindex
       DO il = 1, nsnow
          zi_snow_nopftdim(ix,il) = SUM(zi_snow(ix,il,:)*veget_max(ix,:))
          zf_snow_nopftdim(ix,il) = SUM(zf_snow(ix,il,:)*veget_max(ix,:))
       END DO
    END DO
    
  END SUBROUTINE snowlevels
 
!!
!================================================================================================================================
!! SUBROUTINE   : snow_interpol
!!
!>\BRIEF        This routine interpolates oxygen and methane into snow layers
!!
!! DESCRIPTION : 
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!! 
!! REFERENCE(S) : None
!!
!! FLOWCHART11    : None
!! \n
!_
!================================================================================================================================     
 
  SUBROUTINE snow_interpol (kjpindex,snowO2, snowCH4, zi_snow, zf_snow, veget_max, snowdz)
    
  !! 0. Variable and parameter declaration

    !! 0.1  Input variables     

    INTEGER(i_std), INTENT(in)                       :: kjpindex
    REAL(r_std), DIMENSION(:,:), INTENT(in)          :: snowdz       !! snow depth at each layer [m]
    REAL(r_std), DIMENSION(:,:),INTENT(in)           :: veget_max    !! maximum vegetation fraction                                                           

    !! 0.2 Output variables                                     
                                                                
    !! 0.3 Modified variables                                   

    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)     :: snowO2       !! snow oxygen (g O2/m**3 air)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)     :: snowCH4      !! snow methane (g CH4/m**3 air), needed just for num. scheme
    REAL(r_std), DIMENSION(:,0:,:), INTENT(inout)    :: zf_snow      !! depths at full levels
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)     :: zi_snow      !! depths at intermediate levels

    !! 0.4 Local variables
    REAL(r_std), DIMENSION(kjpindex,nsnow,nvm)                  :: isnow        !! index of first old layer that is deeper
    INTEGER(i_std), DIMENSION(kjpindex,nsnow,nvm)               :: i1,i2        !! indices of the layers used for the inter- or extrapolation
    REAL(r_std), DIMENSION(kjpindex,nsnow,nvm)                  :: snowO2o      !! initial snow oxygen (g O2/m**3 air)
    REAL(r_std), DIMENSION(kjpindex,nsnow,nvm)                  :: snowCH4o     !! initial snow methane (g CH4/m**3 air)
    REAL(r_std), DIMENSION(kjpindex,nvm)                        :: dzio         !! initial distance between two levels
    INTEGER(i_std)                                      :: il, it, ip, ill, iv  !! indices
    REAL(r_std), DIMENSION(kjpindex,0:nsnow,nvm)              :: zfo            !! initial depths at full levels
    REAL(r_std), DIMENSION(kjpindex,nsnow,nvm)                :: zio            !! initial depths at intermediate levels    
    
    
    
    ! 1. save old discretisation and temperatures
    
    zio(:,:,:) = zi_snow(:,:,:)
    
    zfo(:,:,:) = zf_snow(:,:,:)
    
    snowO2o(:,:,:) = snowO2(:,:,:)
    snowCH4o(:,:,:) = snowCH4(:,:,:)
    
    ! 2. new discretisation
    
    CALL snowlevels( kjpindex, snowdz, zi_snow, zf_snow, veget_max)
    
    ! 3. for each new intermediate layer, look for the first old intermediate 
    !   layer that is deeper
    
    DO il = 1, nsnow
       
       isnow(:,il,:) = -1
       
       DO ill = nsnow,1,-1
          
          WHERE ( zio(:,ill,:) .GT. zi_snow(:,il,:) .AND. veget_mask_2d(:,:) )
             
             isnow(:,il,:) = ill
             
          ENDWHERE
          
       ENDDO
       
    ENDDO
   
    ! 4. determine which levels to take for the inter- or extrapolation
    

    DO ip = 1, kjpindex
       DO iv = 1, nvm
          IF ( veget_mask_2d(ip,iv) ) THEN
             DO il = 1, nsnow
                !
                IF ( isnow(ip,il,iv) .EQ. 1  ) THEN
                   !
                   ! 4.1 first old layer is below new layer:
                   !       extrapolation from layers 1 and 2
                   !
                   i1(ip,il,iv) = 1
                   i2(ip,il,iv) = 2
                   !
                ELSEIF ( isnow(ip,il,iv) .EQ. -1 ) THEN
                   !
                   ! 4.2 new layer is below last old layer:
                   !       extrapolation from layers nsnow-1 and nsnow
                   !
                   i1(ip,il,iv) = nsnow-1
                   i2(ip,il,iv) = nsnow
                   !
                ELSE
                   !
                   ! 4.3 new layer is between two old layers: interpolation
                   !
                   i1(ip,il,iv) = isnow(ip,il,iv)-1
                   i2(ip,il,iv) = isnow(ip,il,iv)
                   !
                ENDIF
                
             ENDDO
          ENDIF
       ENDDO
    ENDDO
    
    ! 5. inter- or extrapolate
    
    DO ip = 1, kjpindex
       DO iv = 1, nvm
          IF ( veget_mask_2d(ip,iv) ) THEN
             DO il = 1, nsnow
                dzio(ip,iv) = zio(ip,i2(ip,il,iv),iv) - zio(ip,i1(ip,il,iv),iv)
                
                IF ( dzio(ip,iv) .GT. min_stomate ) THEN
                   
                   snowO2(ip,il,iv) =  snowO2o(ip,i1(ip,il,iv),iv) + &
                        ( zi_snow(ip,il,iv) - zio(ip,i1(ip,il,iv),iv) ) / dzio(ip,iv) * &
                        ( snowO2o(ip,i2(ip,il,iv),iv) - snowO2o(ip,i1(ip,il,iv),iv)  )
                   snowCH4(ip,il,iv) =  snowCH4o(ip,i1(ip,il,iv),iv) + &
                        ( zi_snow(ip,il,iv) - zio(ip,i1(ip,il,iv),iv) ) / dzio(ip,iv) * &
                        ( snowCH4o(ip,i2(ip,il,iv),iv) - snowCH4o(ip,i1(ip,il,iv),iv)  )
                   
                ELSE
                   
                   snowO2(ip,il,iv) = snowO2o(ip,i1(ip,il,iv),iv) 
                   snowCH4(ip,il,iv) = snowCH4o(ip,i1(ip,il,iv),iv) 
                   
                ENDIF
                
             ENDDO
          ENDIF
       ENDDO
       
    ENDDO
  END SUBROUTINE snow_interpol
 
!!
!================================================================================================================================
!! SUBROUTINE   : stomate_soil_carbon_discretization_clear
!!
!>\BRIEF        
!!
!! DESCRIPTION : 
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!! 
!! REFERENCE(S) : None
!!
!! FLOWCHART11    : None
!! \n
!_
!================================================================================================================================     
  SUBROUTINE stomate_soil_carbon_discretization_clear()
    IF (ALLOCATED(veget_mask_2d)) DEALLOCATE(veget_mask_2d)
    IF (ALLOCATED(heights_snow)) DEALLOCATE(heights_snow)
    IF (ALLOCATED(zf_soil)) DEALLOCATE(zf_soil)
    IF (ALLOCATED(zi_soil)) DEALLOCATE(zi_soil)
    IF (ALLOCATED(zf_snow)) DEALLOCATE(zf_snow)
    IF (ALLOCATED(zi_snow)) DEALLOCATE(zi_snow)
    IF (ALLOCATED(alphaO2_soil )) DEALLOCATE(alphaO2_soil )
    IF (ALLOCATED(betaO2_soil )) DEALLOCATE(betaO2_soil )
    IF (ALLOCATED(alphaCH4_soil )) DEALLOCATE(alphaCH4_soil )
    IF (ALLOCATED(betaCH4_soil )) DEALLOCATE(betaCH4_soil )
    IF (ALLOCATED(alphaO2_snow )) DEALLOCATE(alphaO2_snow )
    IF (ALLOCATED(betaO2_snow )) DEALLOCATE(betaO2_snow )
    IF (ALLOCATED(alphaCH4_snow )) DEALLOCATE(alphaCH4_snow )
    IF (ALLOCATED(betaCH4_snow )) DEALLOCATE(betaCH4_snow )
    IF (ALLOCATED(zf_coeff_snow )) DEALLOCATE(zf_coeff_snow )
    IF (ALLOCATED(zi_coeff_snow )) DEALLOCATE(zi_coeff_snow )
    IF (ALLOCATED(mu_snow )) DEALLOCATE(mu_snow )
    
  END SUBROUTINE stomate_soil_carbon_discretization_clear

!!
!================================================================================================================================
!! SUBROUTINE   : initialize_yedoma_carbonstocks
!!
!>\BRIEF        This routine intialize soil carbon in yedoma region
!!
!! DESCRIPTION : 
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!! 
!! REFERENCE(S) : None
!!
!! FLOWCHART11    : None
!! \n
!_
!================================================================================================================================     

  SUBROUTINE initialize_yedoma_carbonstocks(kjpindex, lalo, soilsom_a, soilsom_s, soilsom_p, &
       yedoma_map_filename, yedoma_depth, yedoma_cinit_act, yedoma_cinit_slo, yedoma_cinit_pas, altmax_ind)
   
  !! 0. Variable and parameter declaration

    !! 0.1  Input variables     
 
    INTEGER(i_std), INTENT(in)                                :: kjpindex            !! domain size
    REAL(r_std), DIMENSION(:,:), INTENT(in)                   :: lalo                !! geographic lat/lon
    CHARACTER(LEN=80), INTENT (in)                            :: yedoma_map_filename !! yedoma map
    REAL(r_std), INTENT(in)                                   :: yedoma_depth        !! depth of yedoma carbon stock
    REAL(r_std), INTENT(in)                                   :: yedoma_cinit_act    !! initial active soil C concentration 
    REAL(r_std), INTENT(in)                                   :: yedoma_cinit_slo    !! initial slow soil C concentration 
    REAL(r_std), INTENT(in)                                   :: yedoma_cinit_pas    !! initial passive soil C concentration 
    INTEGER(i_std), DIMENSION(:,:),INTENT(in)	              :: altmax_ind          !! Maximum over the year active-layer index

    !! 0.2 Output variables

    !! 0.3 Modified variables

    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)            :: soilsom_a             !! active soil C concentration
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)            :: soilsom_s             !! slow soil C concentration
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)            :: soilsom_p             !! passive soil C concentration

    !! 0.4 Local variables
    REAL(r_std), DIMENSION(kjpindex)                                 :: yedoma
    INTEGER(i_std)                                                   :: il, ils, ip, ix, iy, imin, jmin, ier, iv
    REAL(r_std)                                                      :: dlon, dlonmin, dlat, dlatmin
    INTEGER(i_std)                                                   :: iml, jml, lml, tml, fid
    REAL(r_std),ALLOCATABLE,DIMENSION(:,:)                           :: xx,yy, yedoma_file
    REAL(r_std),ALLOCATABLE,DIMENSION(:)                             :: x,y
    REAL(r_std)                                                      :: lev(1), date, dt
    INTEGER(i_std)                                                   :: itau(1)
    INTEGER(i_std)                                                   :: yedoma_depth_index, iz
    REAL(r_std), DIMENSION(kjpindex)                                 :: NC_yedoma_act                   !! NC ratio of the active pool for yedoma  
    REAL(r_std), DIMENSION(kjpindex)                                 :: NC_yedoma_slo                   !! NC ratio of the slow pool for yedoma
    REAL(r_std), DIMENSION(kjpindex)                                 :: NC_yedoma_pas                   !! NC ratio of the passive pool for yedoma 
    ! plus bas, on prend la temperature lue dans un fichier climato si celui-ci existe

    IF ( yedoma_map_filename .EQ. "NONE" ) THEN
       yedoma(:) = zero
    ELSE IF ( yedoma_map_filename .EQ. "EVERYWHERE" ) THEN
       yedoma(:) = 1.
    ELSE
       CALL flininfo(yedoma_map_filename,iml, jml, lml, tml, fid)

       ALLOCATE (yy(iml,jml),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'initialize_yedoma_carbonstocks', 'Pb in alloc for yy','','')

       ALLOCATE (xx(iml,jml),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'initialize_yedoma_carbonstocks', 'Pb in alloc for xx','','')

       ALLOCATE (x(iml),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'initialize_yedoma_carbonstocks', 'Pb in alloc for x','','')

       ALLOCATE (y(jml),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'initialize_yedoma_carbonstocks', 'Pb in alloc for y','','')

       ALLOCATE (yedoma_file(iml,jml),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'initialize_yedoma_carbonstocks', 'Pb in alloc for yedoma_file','','')

       CALL flinopen (yedoma_map_filename, .FALSE., iml, jml, lml, &
            xx, yy, lev, tml, itau, date, dt, fid)
       CALL flinget (fid, 'yedoma', iml, jml, lml, tml, &
            1, 1, yedoma_file)
       CALL flinclo (fid)
       ! On suppose que le fichier est regulier.
       ! Si ce n'est pas le cas, tant pis. Les temperatures seront mal
       ! initialisees et puis voila. De toute maniere, il faut avoir
       ! l'esprit mal tourne pour avoir l'idee de faire un fichier de
       ! climatologie avec une grille non reguliere.
       x(:) = xx(:,1)
       y(:) = yy(1,:)
       ! prendre la valeur la plus proche
       DO ip = 1, kjpindex
          dlonmin = HUGE(1.)
          DO ix = 1,iml
             dlon = MIN( ABS(lalo(ip,2)-x(ix)), ABS(lalo(ip,2)+360.-x(ix)), ABS(lalo(ip,2)-360.-x(ix)) )
             IF ( dlon .LT. dlonmin ) THEN
                imin = ix
                dlonmin = dlon
             ENDIF
          ENDDO
          dlatmin = HUGE(1.)
          DO iy = 1,jml
             dlat = ABS(lalo(ip,1)-y(iy))
             IF ( dlat .LT. dlatmin ) THEN
                jmin = iy
                dlatmin = dlat
             ENDIF
          ENDDO
          yedoma(ip) = yedoma_file(imin,jmin)
       ENDDO
       DEALLOCATE (yy)
       DEALLOCATE (xx)
       DEALLOCATE (x)
       DEALLOCATE (y)
       DEALLOCATE (yedoma_file)
    ENDIF
    
    yedoma_depth_index = 0
    DO iz = 1, ngrnd
       IF (znt(iz) .LE. yedoma_depth ) yedoma_depth_index = yedoma_depth_index + 1
    END DO
    WRITE(*,*) 'yedoma_depth_index ', yedoma_depth_index, ' at depth ', yedoma_depth

    IF ( yedoma_depth_index .GT. 0) THEN
       DO ix = 1, kjpindex
          DO iv = 2, nvm  !!! no yedoma carbon for PFT zero.
             IF ( veget_mask_2d(ix,iv) ) THEN
                DO iz = 1, yedoma_depth_index
                   IF (yedoma(ix) .GT. 0.)  THEN
                      IF ( iz .GE. altmax_ind(ix,iv) ) THEN  !!! only put yedoma carbon at base of and below the active layer

!MERGE: provisory fix for NC ratio should be done in a cleaner way later
                            NC_yedoma_act(:)=0.1
                            NC_yedoma_slo(:)=0.1
                            NC_yedoma_pas(:)=0.1
                            CALL getin_p('NC_yedoma_act',NC_yedoma_act)
                            CALL getin_p('NC_yedoma_slo',NC_yedoma_slo)
                            CALL getin_p('NC_yedoma_pas',NC_yedoma_pas)

                         soilsom_a(ix, iz,iv,icarbon) = yedoma_cinit_act
                         soilsom_s(ix, iz,iv,icarbon) = yedoma_cinit_slo
                         soilsom_p(ix, iz,iv,icarbon) = yedoma_cinit_pas
                         soilsom_a(ix, iz,iv,initrogen) = yedoma_cinit_act * NC_yedoma_act(ix)
                         soilsom_s(ix, iz,iv,initrogen) = yedoma_cinit_slo * NC_yedoma_slo(ix)
                         soilsom_p(ix, iz,iv,initrogen) = yedoma_cinit_pas * NC_yedoma_pas(ix)
                      ELSE
                         soilsom_a(ix, iz,iv,:) = zero
                         soilsom_s(ix, iz,iv,:) = zero
                         soilsom_p(ix, iz,iv,:) = zero
                      ENDIF
                   ELSE
                      soilsom_a(ix, iz,iv,:) = zero
                      soilsom_s(ix, iz,iv,:) = zero
                      soilsom_p(ix, iz,iv,:) = zero
                   END IF
                END DO
             ENDIF
          ENDDO
       ENDDO
    ENDIF

  END SUBROUTINE initialize_yedoma_carbonstocks
!!
!================================================================================================================================
!! SUBROUTINE   : sominput
!!
!>\BRIEF        This routine calculate carbon input to the soil
!!
!! DESCRIPTION : 
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!! 
!! REFERENCE(S) : None
!!
!! FLOWCHART11    : None
!! \n
!_
!================================================================================================================================     
  SUBROUTINE sominput(kjpindex,time_step,time,tprof,tsurf,hslong, z_root,altmax, &
       soilsom_a, soilsom_s, soilsom_p,som_input, veget_max, root_depth)
   
  !! 0. Variable and parameter declaration

    !! 0.1  Input variables     
 
    INTEGER(i_std), INTENT(in)				:: kjpindex         !! domain size
    REAL(r_std), INTENT(in)                             :: time_step        !! time step in seconds
    REAL(r_std), INTENT(in)				:: time  
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)           :: tprof            !! Soil temperature (K)
    REAL(r_std), DIMENSION(:), INTENT(in)		:: tsurf	         !! Surface temperature (K)
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)	        :: hslong           !! deep soil humidity
    REAL(r_std), DIMENSION(:,:),INTENT(in)		:: z_root           !! the rooting depth 
    REAL(r_std), DIMENSION(:,:),INTENT(in)		:: altmax           !! Maximum over the year active-layer thickness
    REAL(r_std), DIMENSION(:,:),INTENT(in)              :: veget_max        !! Maximum fraction of vegetation type
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)         :: som_input        !! quantity of organic matter  going into the SOM  pools from litter decomposition (gC/(m**2 of ground)/day)
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)           :: root_depth       !! Node and interface numbers at which the deepest roots
                                                                            !! occur (1 to nslm, unitless)

    !! 0.2 Output variables


    !! 0.3 Modified variables

    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)      :: soilsom_a          !! active soil organic matter
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)      :: soilsom_s          !! slow soil organic matter
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)      :: soilsom_p          !! passive soil organic matter

    !! 0.4 Local variables

    REAL(r_std)                                                 :: dt               !! Time step \f$(dt_sechiba one_day^{-1})$\f
    REAL(r_std), DIMENSION(kjpindex,ncarb,nvm,nelements)                       :: dsom_litter        !! depth-integrated carbon input due to litter decomposition
    REAL(r_std), DIMENSION(kjpindex,ncarb,nvm,nelements)                       :: som_input_finite
    REAL(r_std), DIMENSION(kjpindex,nvm)                                       :: intdep           !! integral depth of carbon deposition   
    REAL(r_std), DIMENSION(kjpindex,ncarb,nvm,nelements)                       :: sominp_correction
    REAL(r_std), DIMENSION(kjpindex,ncarb,nvm,nelements)                       :: som_input_TS
    LOGICAL, SAVE  			                             :: firstcall = .TRUE.
!$OMP THREADPRIVATE(firstcall)
    REAL(r_std), DIMENSION(kjpindex,nvm)                             :: z_lit            !! litter input e-folding depth
    INTEGER                                                          :: il,ic,iv,ix,iele
    INTEGER                                                          :: ipts, ivm
    INTEGER(i_std), SAVE  			                     :: id, id2, id3, id4
!$OMP THREADPRIVATE(id)
!$OMP THREADPRIVATE(id2)
!$OMP THREADPRIVATE(id3)
!$OMP THREADPRIVATE(id4)
    INTEGER                                                          :: recn
    LOGICAL, SAVE                                                    :: correct_carboninput_vertprof
!$OMP THREADPRIVATE(correct_carboninput_vertprof)
    LOGICAL, SAVE                                                    :: new_carbinput_intdepzlit
!$OMP THREADPRIVATE(new_carbinput_intdepzlit)
    REAL(r_std), DIMENSION(ngrnd)                                    :: z_thickness
    REAL(r_std), DIMENSION(ngrnd)                                    :: root_prof
    REAL(r_std), SAVE                                                :: finerootdepthratio = 0.5  !! the ratio of fine root to overall root e-folding depth (for C inputs)
!$OMP THREADPRIVATE(finerootdepthratio)
    REAL(r_std), SAVE                                                :: altrootratio = 0.5        !! the maximum ratio of fine root depth to active layer thickness (for C inputs)
!$OMP THREADPRIVATE(altrootratio)
    INTEGER(i_std)                                                   :: imbc, igrnd      !! Indices (unitless)
    INTEGER(i_std)                                                   :: icarb            !! Indices (unitless)
    REAL(r_std), DIMENSION(kjpindex,nvm,nmbcomp,nelements)           :: check_intern     !! Contains the components of the internal
                                                                                         !! mass balance check for this routine
                                                                                         !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(kjpindex,nvm,nelements)                   :: closure_intern   !! Check closure of internal mass balance
                                                                                         !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(kjpindex,nvm,nelements)                   :: pool_start       !! Start and end pool of this routine 
                                                                                         !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(kjpindex,nvm,nelements)                   :: pool_end         !! Start and end pool of this routine 
                                                                                         !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex

    REAL(r_std), DIMENSION(ngrnd)                                    :: som_profile      !! vertical profile (0-1, unitless)
    INTEGER(r_std)                                                   :: best_layer       !! number of the layer of which the depth
                                                                                         !! best matches a specific target depth 
                                                                                         !! (0-ngrnd, unitless)
    REAL(r_std), DIMENSION(kjpindex,ncarb,ngrnd,nvm,nelements)       :: dsom_litter_z    !! depth_dependent organic matter input due to litter
!_ ================================================================================================================================

    
       IF (firstcall) THEN

          !Config Key   = new_carbinput_intdepzlit
          !Config Desc  = ???
          !Config Def   = n
          !Config If    = OK_SOIL_CARBON_DISCRETIZATION
          !Config Help  = 
          !Config Units = [flag]
          new_carbinput_intdepzlit = .FALSE.
          CALL getin_p('new_carbinput_intdepzlit', new_carbinput_intdepzlit)

          !
          !Config Key   = correct_carboninput_vertprof
          !Config Desc  = ???
          !Config Def   = n
          !Config If    = OK_SOIL_CARBON_DISCRETIZATION
          !Config Help  = 
          !Config Units = [flag]
          correct_carboninput_vertprof = .TRUE.
          CALL getin_p('correct_carboninput_vertprof', correct_carboninput_vertprof)

          firstcall = .FALSE.
          !
       ENDIF ! firstcall
       
       ! Calculate thickness of the soil layers for later use. Could be put into
       ! firstcall loop but it will then needs to be defined as SAVE
       DO il = 1, ngrnd
          z_thickness(il) = zf_soil(il) - zf_soil(il-1) 
       END DO

       !! Initialize check for mass balance closure
       !  Carbon available at the start of this subroutine.
       IF (err_act.GT.1) THEN

          dt = dt_sechiba/one_day
          pool_start(:,:,:) = zero
          DO iele = 1,nelements
             DO igrnd = 1,ngrnd
                pool_start(:,:,iele) = pool_start(:,:,iele) + &
                     (soilsom_a(:,igrnd,:,iele) + soilsom_s(:,igrnd,:,iele) + soilsom_p(:,igrnd,:,iele)) * &
                     z_thickness(igrnd) * veget_max(:,:)
             END DO
          ENDDO
          
          DO iele = 1,nelements
             DO icarb = 1,ncarb
                pool_start(:,:,iele) = pool_start(:,:,iele) + &
                     som_input(:,icarb,:,iele) * dt * veget_max(:,:)
             ENDDO
          ENDDO

       END IF ! err_act.GT.1

       !
       ! 1. Litter input and decomposition
       !
       ! add up the soil carbon from all veg pools, and change units from (gC/(m**2 of ground)/day) to gC/m^2 per timestep      
       som_input_TS(:,:,:,:) = som_input(:,:,:,:)*time_step/one_day
       
       ! The objective is to distribute the litter inputs over the vertical
       ! soil layers. There are two sources that matter: (1) aboveground litter
       ! and (2) belowground litter. The vertical distribution of the som_input
       ! in the soil should thus follow a vertical profile described by the 
       ! weighted average of the vertical profile of the belowground som input
       ! and the belowground som_input. The belowground litter comes from the roots
       ! and could simply follow the structural root distribution. The problem is
       ! that for the moment the variable som_input does not distinguish between
       ! aboveground and belowground sources. The weighting factor is therefore 
       ! not available yet. Either this separation is added to the code or some 
       ! assumptions have to be made here. Because the above ground litter 
       ! decomposes at the top of the soil, most of its som_input is expeceted in 
       ! the top soil layers. Because Earth worms, ants and other mesofauna that 
       ! could transport litter to deeper soil layers are not yet accounted for, 
       ! assumptions will be required.

       ! The som_inputs of the belowground litter follow
       ! and exponentially decreasing profile whereas the aboveground som_inputs
       ! will also follow an exponentially decreasing profile. The weighted mean
       ! of these profiles is expecetd to be an exponentially decreasing profile
       ! with its exponential to be less than humcste as we expect that the som_input
       ! in the top leayers will exceed the input from just the roots in the top
       ! layers because the aboveground litter pool is expected to exceed the below
       ! ground litter pool. This reasoning is reflected in the calculation of z_lit
       ! below. z_lit is the parameter that will determine the shape of the 
       ! exponentially decreasing som_inputs with increasing soil depth.

       ! Note that the calculation of z_lit is based on soil depths. Given z_lit
       ! is a bit of an arbitrary parameter, this is not wrong but it is confusing.
       ! The depth of the integration is determined by intdep and does not depend
       ! on z_lit. If z_lit is set to 0.2 m but intdep is set at 2 m, there will still
       ! be litter input at the deeper layers. Z_lit determines the shape (and through
       ! the shape where the majority of the inputs end up) but it has no strict
       ! control over the depth. Which is suggested by defining z_lit as a depth.
       ! If z_lit1 > z_lit2, th eresulting profile of z_lit1 will have more som inputs
       ! at deeper layers than z_lit2. From that point of view z_lit could be defined
       ! as a proxy for depth.
       
       !
       ! 2. Carbon input e-folding depth. We distribute with e-depth = min(z_root,intdep) 
       !    and integral depth = min(altmax,z_org)
       !    e-folding depth cannot be greater than integral depth
       IF ( .NOT. new_carbinput_intdepzlit ) THEN
          ! change to make intdep equal to z_root alone.
          z_lit(:,:) = z_root(:,:)
          intdep(:,:) = z_root(:,:)
       ELSE
          !change to separate e-folding depths for roots from total depth over which to integrate
          ! z_lit is the e-folding depth
          DO ivm=1,nvm
             z_lit(:,ivm) = MIN(finerootdepthratio/humcste(ivm), altmax(:,ivm)*altrootratio)
          END DO
          DO ipts = 1,kjpindex
             DO ivm = 1,nvm
                ! intdep is the maximum depth of integration; Use zdr to be consistent with
                ! root_depth. Both variables assume that the indices go from 0:nslm
                intdep(ipts,ivm) = zdr(NINT(root_depth(ipts,ivm,iinterface)))
             END DO
          END DO
       ENDIF
       
       !
       ! 3. Carbon input. 
       !
       dsom_litter_z(:,:,:,:,:) = zero
       dsom_litter(:,:,:,:)=zero
       DO ipts = 1,kjpindex
          DO ivm = 1,nvm

             ! The original code was using zi_soil. If we use zi_soil here we introduce 
             ! inconsistencies in the mass balance because we will distribute the 
             ! carbon over zi_soil but then integrate it over zf_soil. Given the top 
             ! layers are rather thin and the number of layers seems arbitrary anyway, 
             ! it seems acceptable to move from zi_soil to zf_soil.
             ! We always need to integrate to the exact depth of one of the layers, 
             ! otherwise we will add a lot of complexity to close the mass balance. If 
             ! intdep(:,:) .LT. zi_soil(2), we set the depth to an exact layer. A solution
             ! was added in case intdep(:,:) .GT. zi_soil(2). In that case the nearest 
             ! exact layer was used to distribute the carbon and nitrogen.
             IF ( intdep(ipts,ivm) .LT. zf_soil(2) ) THEN
                ! Litter is decomposed somehow (?) even when alt == 0. To avoid carbon loss,
                ! we distribute this carbon within the first 2 soil layers when alt == 0
                ! Adding EPSILON(zero) avoids problems with the next IF-statement where
                ! zf_soil is compared against intdep.
                intdep(ipts,ivm) = zf_soil(2) + EPSILON(0.)
             ELSE
                ! Find the layer that best represents the depth of interest. Note
                ! that before intdep was a depth in m, here it is being transfered 
                ! into a layer number.
                best_layer = MINLOC(ABS(zf_soil(:) - intdep(ipts,ivm)),dim=1)
                intdep(ipts,ivm) = zf_soil(best_layer) + EPSILON(0.)
             ENDIF
             IF ( z_lit(ipts,ivm) .LT. zf_soil(2) ) THEN
                ! Litter is decomposed somehow (?) even when alt == 0. To avoid carbon loss,
                ! we distribute this carbon within the first 2 soil layers when alt == 0
                z_lit(ipts,ivm) = zf_soil(2)
             ELSE
                ! Find the layer that best represents the depth of interest
                best_layer = MINLOC(ABS(zf_soil(:) - z_lit(ipts,ivm)),dim=1)
                z_lit(ipts,ivm) = zf_soil(best_layer)
             ENDIF

             ! Initialize
             som_profile(:) = zero
  
             ! Rewritten the calculation of the redistribution factor. Irrespective
             ! of whether zi_soil or zf_soil is used, the original code appeared 
             ! bugged to me. It did not account for the first layer. It divided by
             ! som_input_TS by z_lit which implies it redistributed the input
             ! twice (the division is a uniform redistribution, the EXP an exponential
             ! redistribution). Division by z_lit results in a change in the units
             ! this conversion is now taken care of later in the code.
             DO il = 1, ngrnd
                ! Calculate the som profile
                IF ( zf_soil(il) .LT. intdep(ipts,ivm) .AND. veget_mask_2d(ipts,ivm) ) THEN
                   ! The first part of the equation, i.e., 1/(1-EXP(-intdep/z_lit) 
                   ! is a scaling factor to ensure that the profile down to 
                   ! intdep will add up to 1 (note that z_lit determines the shape of
                   ! exponential function that will be used to redistribute the som
                   ! over the vertical layers).
                   som_profile(il) = un / ( un - EXP( -intdep(ipts,ivm) / z_lit(ipts,ivm) ) ) * &
                        ( EXP(-zf_soil(il-1)/z_lit(ipts,ivm))  - &
                        EXP( -zf_soil(il)/z_lit(ipts,ivm) ) )
                ELSE
                   ! This layer is too far down and is no longer accounted for
                   som_profile(il) = zero
                END IF
             END DO

             ! Check for errors. If none, overcome possible precision issues
             IF (veget_mask_2d(ipts,ivm)) THEN
                IF (ABS(SUM(som_profile(:))-un).GT.min_stomate) THEN
                   ! The above calculation is less precise than expeceted
                   IF (err_act.GT.1) THEN
                      WRITE(numout,*) 'ipts, ivm, ', ipts,ivm
                      WRITE(numout,*) 'zf_soil, ', zf_soil(:)
                      WRITE(numout,*) 'intdep, z_lit,', intdep(ipts,ivm),z_lit(ipts,ivm)
                      WRITE(numout,*) 'som_profile, ', som_profile(:)
                      WRITE(numout,*) 'Sum of som_profile, ', SUM(som_profile(:))
                      CALL ipslerr_p(err_act,'som_input','The sum of the som profile differs from 1',&
                           'mass balance will not be closed','')
                   ENDIF
                ELSE
                   ! The above calculation should result in a som_profile of 0.9999999 which is
                   ! pretty good but not good enough to keep the mass balance closed. Put the 
                   ! residual in the top layer (just because that layer should always be present)
                   som_profile(1) = un - SUM(som_profile(2:ngrnd))
                END IF
             END IF

             ! Use the profile to redistribute the som_input
             DO ic = 1, ncarb
                DO iele = 1, nelements
                   DO il = 1, ngrnd
                      ! Divide by the tickness of the layer, not the tickness of the whole 
                      ! profile (which is the case when using z_lit). The unit of som_input_Ts
                      ! is gC/m^2 per timestep. We want dsom_litter_z to be in gC/m^3 per
                      ! timestep so divide by the layer depth over which the som input is to 
                      ! be distributed.
                      dsom_litter_z(ipts,ic,il,ivm,iele) = som_input_TS(ipts,ic,ivm,iele) * &
                           som_profile(il)/ z_thickness(il)
                      ! The original code already used zf_soil. This may have resulted in an
                      ! inconsitency. Given that zi_soil has now been changed to zf_soil. The
                      ! original code is consistent with the new approach that uses zf_soil
                      dsom_litter(ipts,ic,ivm,iele) = dsom_litter(ipts,ic,ivm,iele) + &
                           dsom_litter_z(ipts,ic,il,ivm,iele) * z_thickness(il)
                   END DO
                END DO
             END DO

          END DO
       END DO
        
       ! Update the active, slow and passive soilsom pools
       DO il = 1, ngrnd
          DO iele = 1, nelements
             WHERE ( veget_mask_2d(:,:) ) 
                soilsom_a(:,il,:,iele) = soilsom_a(:,il,:,iele) + dsom_litter_z(:,iactive,il,:,iele)
                soilsom_s(:,il,:,iele) = soilsom_s(:,il,:,iele) + dsom_litter_z(:,islow,il,:,iele)
                soilsom_p(:,il,:,iele) = soilsom_p(:,il,:,iele) + dsom_litter_z(:,ipassive,il,:,iele)
             END WHERE
          ENDDO
       END DO
       
       !! Check mass balance closure
       IF (err_act.GT.1) THEN
        
          !! 4.1.1 Calculate components of the mass balance
          pool_end(:,:,:) = zero 
          DO iele = 1,nelements
             DO igrnd = 1,ngrnd
                pool_end(:,:,iele) = pool_end(:,:,iele) + &
                     (soilsom_a(:,igrnd,:,iele) + soilsom_s(:,igrnd,:,iele) + soilsom_p(:,igrnd,:,iele)) * &
                     z_thickness(igrnd) * veget_max(:,:)
             END DO
          END DO

          !! 4.1.2 Calculate mass balance
          check_intern(:,:,:,:) = zero
          check_intern(:,:,iatm2land,:) = zero
          check_intern(:,:,iland2atm,icarbon) = -un * zero 
          check_intern(:,:,iland2atm,initrogen) = -un * zero
          check_intern(:,:,ilat2out,:) = -un *zero
          check_intern(:,:,ilat2in,:) =  zero
          check_intern(:,:,ipoolchange,:) = -un * (pool_end(:,:,:) - &
               pool_start(:,:,:))
          closure_intern = zero
          DO imbc = 1,nmbcomp
             DO iele = 1,nelements
                closure_intern(:,:,iele) = closure_intern(:,:,iele) + &
                     check_intern(:,:,imbc,iele)
             END DO
          ENDDO

          CALL check_mass_balance("sominput", closure_intern, kjpindex, pool_end, &
               pool_start, veget_max, 'pft')

       END IF ! err_act.GT.1
      
  END SUBROUTINE sominput

!!
!================================================================================================================================
!! SUBROUTINE   : cryoturbate
!!
!>\BRIEF        This routine calculates cryoturbation process
!!
!! DESCRIPTION : 
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!! 
!! REFERENCE(S) : None
!!
!! FLOWCHART11    : None
!! \n
!_
!================================================================================================================================     
  
  SUBROUTINE cryoturbate(kjpindex, time_step, altmax_ind, deepSOM_a, deepSOM_s, deepSOM_p, &
       action, diff_k_const, bio_diff_k_const, altmax_lastyear, fixed_cryoturbation_depth, &

       veget_max)

  !! 0. Variable and parameter declaration

    !! 0.1  Input variables     

    INTEGER(i_std), INTENT(in)			      :: kjpindex         !! domain size
    REAL(r_std), INTENT(in)                           :: time_step        !! time step in seconds
    INTEGER(i_std), DIMENSION(:,:),INTENT(in)	      :: altmax_ind       !! Maximum over the year active-layer index
    REAL(r_std), DIMENSION(:,:),INTENT(in)            :: altmax_lastyear  !! Maximum over the year active-layer thickness
    CHARACTER(LEN=*), INTENT(in)                      :: action           !! what to do
    REAL(r_std), INTENT(in)                           :: diff_k_const
    REAL(r_std), INTENT(in)                           :: bio_diff_k_const
    REAL(r_std), DIMENSION(:,:),INTENT(in)            :: veget_max         !! Maximum vegetation fraction

    !! 0.2 Output variables 

    !! 0.3 Modified variables
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)   :: deepSOM_a          !! soil soil organic matter (g/m**3) active
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)   :: deepSOM_s          !! soil soil organic matter (g/m**3) slow
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)   :: deepSOM_p          !! soil soil organic matter (g/m**3) passive
    REAL(r_std), DIMENSION(:,:),INTENT(inout)        :: fixed_cryoturbation_depth  !! depth to hold cryoturbation to for fixed runs

    !! 0.4 Local variables
    LOGICAL, SAVE			                       :: firstcall = .TRUE.
!$OMP THREADPRIVATE(firstcall)
    LOGICAL, SAVE			                       :: use_new_cryoturbation
!$OMP THREADPRIVATE(use_new_cryoturbation)
    INTEGER, SAVE			                       :: cryoturbation_method
!$OMP THREADPRIVATE(cryoturbation_method)
    REAL(r_std), DIMENSION(kjpindex,nvm,nelements)                       :: altSOM_a_old       !! soil organic matter (g/m**2) active integrated over active layer before cryoturbation
    REAL(r_std), DIMENSION(kjpindex,nvm,nelements)                       :: altSOM_s_old       !! soil organic matter (g/m**2) slow integrated over active layer before cryoturbation
    REAL(r_std), DIMENSION(kjpindex,nvm,nelements)                       :: altSOM_p_old       !! soil organic matter (g/m**2) passive integrated over active layer before cryoturbation
    REAL(r_std), DIMENSION(kjpindex,nvm,nelements)                       :: altSOM_a
    REAL(r_std), DIMENSION(kjpindex,nvm,nelements)                       :: altSOM_s
    REAL(r_std), DIMENSION(kjpindex,nvm,nelements)                       :: altSOM_p
    INTEGER(i_std), PARAMETER                                  :: n_totakefrom = 3 !! how many surface layers to subtract from in mass balance
    REAL(r_std), DIMENSION(kjpindex,nvm,nelements)                       :: surfSOM_totake_a   !! active soil organic matter to subtract from surface layers to maintain mass balance (g/m**3) 
    REAL(r_std), DIMENSION(kjpindex,nvm,nelements)                       :: surfSOM_totake_s   !! slow soil organic matter to subtract from surface layers to maintain mass balance (g/m**3) 
    REAL(r_std), DIMENSION(kjpindex,nvm,nelements)                       :: surfSOM_totake_p   !! passive soil organic matter to subtract from surface layers to maintain mass balance (g/m**3) 
    REAL(r_std), DIMENSION(kjpindex,nvm)                       :: error_a
    REAL(r_std), DIMENSION(kjpindex,nvm)                       :: error_s
    REAL(r_std), DIMENSION(kjpindex,nvm)                       :: error_p
    INTEGER(i_std)                                             :: ip, il, ier, iv, iele
    CHARACTER(LEN=20), SAVE                                    :: last_action = 'not called'
!$OMP THREADPRIVATE(last_action)
    INTEGER(i_std)                                             :: cryoturb_date
    REAL(r_std), SAVE                                          :: max_cryoturb_alt
!$OMP THREADPRIVATE(max_cryoturb_alt)
    REAL(r_std), SAVE                                          :: min_cryoturb_alt
!$OMP THREADPRIVATE(min_cryoturb_alt)
    REAL(r_std), SAVE                                          :: bioturbation_depth
!$OMP THREADPRIVATE(bioturbation_depth)
    LOGICAL, SAVE			                       :: reset_fixed_cryoturbation_depth
!$OMP THREADPRIVATE(reset_fixed_cryoturbation_depth)
    LOGICAL, SAVE			                       :: use_fixed_cryoturbation_depth
!$OMP THREADPRIVATE(use_fixed_cryoturbation_depth)
    REAL(r_std), DIMENSION(kjpindex,nvm)	               :: cryoturbation_depth
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm)	               :: diff_k               !! Diffusion constant (m^2/s)
    REAL(r_std), DIMENSION(kjpindex,nvm)	               :: xe_a
    REAL(r_std), DIMENSION(kjpindex,nvm)	               :: xe_s
    REAL(r_std), DIMENSION(kjpindex,nvm)	               :: xe_p
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm)	               :: xc_cryoturb
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm)	               :: xd_cryoturb
    INTEGER(i_std)                                        :: imbc, igrnd      !! Indices (unitless)
    INTEGER(i_std)                                        :: inbpools         !! Indices (unitless)
    REAL(r_std), DIMENSION(kjpindex,nvm,nmbcomp,nelements):: check_intern     !! Contains the components of the internal
                                                                              !! mass balance check for this routine
                                                                              !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(kjpindex,nvm,nelements)        :: closure_intern   !! Check closure of internal mass balance
                                                                              !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(kjpindex,nvm,nelements)        :: pool_start       !! Start and end pool of this routine 
                                                                              !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(kjpindex,nvm,nelements)        :: pool_end         !! Start and end pool of this routine 
                                                                              !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex


!_ ================================================================================================================================

    
    ! 1. ensure that we do not repeat actions
    !
    IF ( action .EQ. last_action ) THEN
       !
       WRITE(*,*) 'CANNOT TAKE THE SAME ACTION TWICE: ',TRIM(action)
       CALL ipslerr_p(3,'cryoturbate', 'Cannot take same action twice','action equal last_action','')
       !
    ENDIF
    
    IF (firstcall) THEN
          
          ! 2. faire les trucs du debut
          
          ! 2.1 allocation des variables
          ALLOCATE (alpha_a(kjpindex,ngrnd,nvm,nelements),stat=ier)
          IF (ier /= 0) CALL ipslerr_p(3,'cryoturbate', 'Pb in alloc for alpha_a','','')
	  alpha_a(:,:,:,:)=0.

          ALLOCATE (alpha_s(kjpindex,ngrnd,nvm,nelements),stat=ier)
          IF (ier /= 0) CALL ipslerr_p(3,'cryoturbate', 'Pb in alloc for alpha_s','','')
	  alpha_s(:,:,:,:)=0.

          ALLOCATE (alpha_p(kjpindex,ngrnd,nvm,nelements),stat=ier)
          IF (ier /= 0) CALL ipslerr_p(3,'cryoturbate', 'Pb in alloc for alpha_p','','')
	  alpha_p(:,:,:,:)=0.

          ALLOCATE (mu_soil_rev(kjpindex,nvm),stat=ier)
          IF (ier /= 0) CALL ipslerr_p(3,'cryoturbate', 'Pb in alloc for mu_soil_rev','','')
	  mu_soil_rev(:,:)=0.

          ALLOCATE (beta_a(kjpindex,ngrnd,nvm,nelements),stat=ier)
          IF (ier /= 0) CALL ipslerr_p(3,'cryoturbate', 'Pb in alloc for beta_a','','')
	  beta_a(:,:,:,:)=0.

          ALLOCATE (beta_s(kjpindex,ngrnd,nvm,nelements),stat=ier)
          IF (ier /= 0) CALL ipslerr_p(3,'cryoturbate', 'Pb in alloc for beta_s','','')
	  beta_s(:,:,:,:)=0.
          
          ALLOCATE (beta_p(kjpindex,ngrnd,nvm,nelements),stat=ier)
          IF (ier /= 0) CALL ipslerr_p(3,'cryoturbate', 'Pb in alloc for beta_p','','')
 	  beta_p(:,:,:,:)=0.
        
          ALLOCATE (cryoturb_location(kjpindex,nvm),stat=ier)
          IF (ier /= 0) CALL ipslerr_p(3,'cryoturbate', 'Pb in alloc for cryoturb_location','','')

          ALLOCATE (bioturb_location(kjpindex,nvm),stat=ier)
          IF (ier /= 0) CALL ipslerr_p(3,'cryoturbate', 'Pb in alloc for bioturb_location','','')
           
          cryoturb_location(:,:) = .false.
          
          !
          !Config Key   = use_new_cryoturbation
          !Config Desc  = use new scheme to calculate cryoturbation
          !Config Def   = n
          !Config If    = OK_SOIL_CARBON_DISCRETIZATION
          !Config Help  = 
          !Config Units = [flag]
          use_new_cryoturbation = .false.
          CALL getin_p('use_new_cryoturbation', use_new_cryoturbation)
          !
          !Config Key   = cryoturbation_method
          !Config Desc  = Which method should be used to calculate cryoturbation
          !Config Def   = 1
          !Config If    =  OK_SOIL_CARBON_DISCRETIZATION 
          !Config Help  = 1: linear dropoff to zero between alt and 2*alt;
          !               2: exponential dropoff with e-folding up to a distance equal 
          !               to alt, below the active layer; 3: exponential dropoff with 
          !               e-folding up to a distance equal to alt; 4: starting at surface
          !               linear dropoff to zero between alt and 3*alt
          !Config Units = []
          cryoturbation_method = 4
          CALL getin_p('cryoturbation_method', cryoturbation_method)
          !
          !Config Key   = max_cryoturb_alt
          !Config Desc  = ???
          !Config Def   = 1
          !Config If    = OK_SOIL_CARBON_DISCRETIZATION
          !Config Help  = 
          !Config Units = [???]
          max_cryoturb_alt = 3.
          CALL getin_p('max_cryoturb_alt',max_cryoturb_alt)
          !
          !Config Key   = min_cryoturb_alt
          !Config Desc  = ???
          !Config Def   = 1
          !Config If    = OK_SOIL_CARBON_DISCRETIZATION
          !Config Help  = 
          !Config Units = [???]
          min_cryoturb_alt = 0.01
          CALL getin_p('min_cryoturb_alt',min_cryoturb_alt)
          !
          !Config Key   = reset_fixed_cryoturbation_depth
          !Config Desc  = reset fixed cryoturbation depth
          !Config Def   = n
          !Config If    = OK_SOIL_CARBON_DISCRETIZATION
          !Config Help  = 
          !Config Units = [flag]
          reset_fixed_cryoturbation_depth = .FALSE.
          CALL getin_p('reset_fixed_cryoturbation_depth',reset_fixed_cryoturbation_depth)
          IF (reset_fixed_cryoturbation_depth) THEN
             fixed_cryoturbation_depth = altmax_lastyear
          ENDIF
          !
          !Config Key   = use_fixed_cryoturbation_depth
          !Config Desc  = use fixed cryoturbation depth
          !Config Def   = n
          !Config If    = OK_SOIL_CARBON_DISCRETIZATION
          !Config Help  = 
          !Config Units = [flag]
          use_fixed_cryoturbation_depth = .FALSE.
          CALL getin_p('use_fixed_cryoturbation_depth',use_fixed_cryoturbation_depth)
          bioturb_location(:,:) = .false.
          !
          !Config Key   = bioturbation_depth
          !Config Desc  = maximum bioturbation depth 
          !Config Def   = 2
          !Config If    = OK_SOIL_CARBON_DISCRETIZATION
          !Config Help  = 
          !Config Units =[m]
          bioturbation_depth = 2.
          CALL getin_p('bioturbation_depth',bioturbation_depth)
          
          firstcall = .FALSE.
    ENDIF

    !! Initialize check for mass balance closure
    IF (err_act.GT.1) THEN

       pool_start(:,:,:) = zero
       
       DO iele = 1,nelements
          DO igrnd = 1,ngrnd
             pool_start(:,:,iele) = pool_start(:,:,iele) + &
                  (deepSOM_a(:,igrnd,:,iele) + deepSOM_s(:,igrnd,:,iele) + deepSOM_p(:,igrnd,:,iele)) * &
                  (zf_soil(igrnd)-zf_soil(igrnd-1)) * veget_max(:,:)
          END DO
       END DO

    END IF ! err_act.gt.1
       
    IF ( action .EQ. 'diffuse' ) THEN
          ! 1. calculate the total soil organic matter in the active layer
          altSOM_a_old(:,:,:) = zero
          altSOM_s_old(:,:,:) = zero
          altSOM_p_old(:,:,:) = zero
          altSOM_a(:,:,:) = zero
          altSOM_s(:,:,:) = zero
          altSOM_p(:,:,:) = zero

          DO ip = 1, kjpindex
             DO iv = 1, nvm
               DO iele = 1, nelements
                IF ( cryoturb_location(ip,iv) .OR. bioturb_location(ip,iv) )THEN 
                   ! 1. calculate the total soil organic matter 
                   DO il = 1, ngrnd
                      altSOM_a_old(ip,iv,iele) = altSOM_a_old(ip,iv,iele) + deepSOM_a(ip,il,iv,iele)*(zf_soil(il)-zf_soil(il-1))
                      altSOM_s_old(ip,iv,iele) = altSOM_s_old(ip,iv,iele) + deepSOM_s(ip,il,iv,iele)*(zf_soil(il)-zf_soil(il-1))
                      altSOM_p_old(ip,iv,iele) = altSOM_p_old(ip,iv,iele) + deepSOM_p(ip,il,iv,iele)*(zf_soil(il)-zf_soil(il-1))
                   ENDDO
                   
                   ! 2. diffuse the soil organic matter                 
                   deepSOM_a(ip,1,iv,iele) = (deepSOM_a(ip,1,iv,iele)+mu_soil_rev(ip,iv)*beta_a(ip,1,iv,iele)) / &
                        (1.+mu_soil_rev(ip,iv)*(1.-alpha_a(ip,1,iv,iele)))
                   deepSOM_s(ip,1,iv,iele) = (deepSOM_s(ip,1,iv,iele)+mu_soil_rev(ip,iv)*beta_s(ip,1,iv,iele)) / &
                        (1.+mu_soil_rev(ip,iv)*(1.-alpha_s(ip,1,iv,iele)))
                   deepSOM_p(ip,1,iv,iele) = (deepSOM_p(ip,1,iv,iele)+mu_soil_rev(ip,iv)*beta_p(ip,1,iv,iele)) / &
                        (1.+mu_soil_rev(ip,iv)*(1.-alpha_p(ip,1,iv,iele)))

                   DO il = 2, ngrnd
                      deepSOM_a(ip,il,iv,iele) = alpha_a(ip,il-1,iv,iele)*deepSOM_a(ip,il-1,iv,iele) + beta_a(ip,il-1,iv,iele)
                      deepSOM_s(ip,il,iv,iele) = alpha_s(ip,il-1,iv,iele)*deepSOM_s(ip,il-1,iv,iele) + beta_s(ip,il-1,iv,iele)
                      deepSOM_p(ip,il,iv,iele) = alpha_p(ip,il-1,iv,iele)*deepSOM_p(ip,il-1,iv,iele) + beta_p(ip,il-1,iv,iele)
                   ENDDO

                   ! 3. recalculate the total soil organic matter
                   DO il = 1, ngrnd
                      altSOM_a(ip,iv,iele) = altSOM_a(ip,iv,iele) + deepSOM_a(ip,il,iv,iele)*(zf_soil(il)-zf_soil(il-1))
                      altSOM_s(ip,iv,iele) = altSOM_s(ip,iv,iele) + deepSOM_s(ip,il,iv,iele)*(zf_soil(il)-zf_soil(il-1))
                      altSOM_p(ip,iv,iele) = altSOM_p(ip,iv,iele) + deepSOM_p(ip,il,iv,iele)*(zf_soil(il)-zf_soil(il-1))
                   ENDDO

                   
                   IF ( altSOM_a_old(ip,iv,iele) > min_stomate .AND. &
                        (ABS(altSOM_a(ip,iv,iele)-altSOM_a_old(ip,iv,iele))/altSOM_a_old(ip,iv,iele).GT. min_stomate) ) THEN
                      WRITE (numout,*) 'DZ warn: cryoturbate: total C or N not conserved','iele=',iele, 'ip=',ip,'iv=',iv, &
                           'A,diff=',altSOM_a(ip,iv,iele),altSOM_a_old(ip,iv,iele),altSOM_a(ip,iv,iele)-altSOM_a_old(ip,iv,iele), &
                           (altSOM_a(ip,iv,iele)-altSOM_a_old(ip,iv,iele))/altSOM_a_old(ip,iv,iele)
                      CALL ipslerr_p (3,'cryoturbate','','','')
                   ENDIF

                   IF ( altSOM_s_old(ip,iv,iele) > min_stomate .AND. &
                        (ABS(altSOM_s(ip,iv,iele)-altSOM_s_old(ip,iv,iele))/altSOM_s_old(ip,iv,iele).GT. min_stomate) ) THEN
                      WRITE (numout,*) 'DZ warn: cryoturbate: total C or N not conserved','iele=',iele,'ip=',ip,'iv=',iv, &
                           'S,diff=',altSOM_s(ip,iv,iele),altSOM_s_old(ip,iv,iele),altSOM_s(ip,iv,iele)-altSOM_s_old(ip,iv,iele), &
                           (altSOM_s(ip,iv,iele)-altSOM_s_old(ip,iv,iele))/altSOM_s_old(ip,iv,iele)
                      CALL ipslerr_p (3,'cryoturbate','','','')
                   ENDIF

                   IF ( altSOM_p_old(ip,iv,iele) > min_stomate .AND. &
                        (ABS(altSOM_p(ip,iv,iele)-altSOM_p_old(ip,iv,iele))/altSOM_p_old(ip,iv,iele).GT. min_stomate) ) THEN
                      WRITE (numout,*) 'DZ warn: cryoturbate: total C or N not conserved','iele=',iele, 'ip=',ip,'iv=',iv, &
                           'P,diff=',altSOM_p(ip,iv,iele),altSOM_p_old(ip,iv,iele),altSOM_p(ip,iv,iele)-altSOM_p_old(ip,iv,iele), &
                           (altSOM_p(ip,iv,iele)-altSOM_p_old(ip,iv,iele))/altSOM_p_old(ip,iv,iele)
                      CALL ipslerr_p (3,'cryoturbate','','','')
                   ENDIF

                   ! 4. subtract the organic matter in the top layer(s) so that the total organic matter content of the active layer is conserved.             
                   ! for now remove this correction term...
!                   surfC_totake_a(ip,iv) = (altC_a(ip,iv)-altC_a_old(ip,iv))/(zf_soil(altmax_ind(ip,iv))-zf_soil(0))
!                   surfC_totake_s(ip,iv) = (altC_s(ip,iv)-altC_s_old(ip,iv))/(zf_soil(altmax_ind(ip,iv))-zf_soil(0))
!                   surfC_totake_p(ip,iv) = (altC_p(ip,iv)-altC_p_old(ip,iv))/(zf_soil(altmax_ind(ip,iv))-zf_soil(0))
!                   deepC_a(ip,1:altmax_ind(ip,iv),iv) = deepC_a(ip,1:altmax_ind(ip,iv),iv) - surfC_totake_a(ip,iv)
!                   deepC_s(ip,1:altmax_ind(ip,iv),iv) = deepC_s(ip,1:altmax_ind(ip,iv),iv) - surfC_totake_s(ip,iv)
!                   deepC_p(ip,1:altmax_ind(ip,iv),iv) = deepC_p(ip,1:altmax_ind(ip,iv),iv) - surfC_totake_p(ip,iv)
!
!                   ! if negative values appear, we don't subtract the delta-C from top layers
!                   IF (ANY(deepC_a(ip,1:altmax_ind(ip,iv),iv) .LT. zero) ) THEN
!                      deepC_a(ip,1:altmax_ind(ip,iv),iv)=deepC_a(ip,1:altmax_ind(ip,iv),iv)+surfC_totake_a(ip,iv)
!                      IF (altC_a(ip,iv) .GT. zero) THEN
!                         deepC_a(ip,:,iv)=deepC_a(ip,:,iv)*altC_a_old(ip,iv)/altC_a(ip,iv)
!                      ENDIF
!                   ENDIF
!                   IF (ANY(deepC_s(ip,1:altmax_ind(ip,iv),iv) .LT. zero) ) THEN
!                      deepC_s(ip,1:altmax_ind(ip,iv),iv)=deepC_s(ip,1:altmax_ind(ip,iv),iv)+surfC_totake_s(ip,iv)
!                      IF (altC_s(ip,iv) .GT. zero) THEN
!                         deepC_s(ip,:,iv)=deepC_s(ip,:,iv)*altC_s_old(ip,iv)/altC_s(ip,iv)
!                      ENDIF
!                   ENDIF
!                   IF (ANY(deepC_p(ip,1:altmax_ind(ip,iv),iv) .LT. zero) ) THEN
!                      deepC_p(ip,1:altmax_ind(ip,iv),iv)=deepC_p(ip,1:altmax_ind(ip,iv),iv)+surfC_totake_p(ip,iv)
!                      IF (altC_p(ip,iv) .GT. zero) THEN
!                         deepC_p(ip,:,iv)=deepC_p(ip,:,iv)*altC_p_old(ip,iv)/altC_p(ip,iv)
!                      ENDIF
!                   ENDIF

          	   ! Consistency check. Potentially add to STRICT_CHECK flag
                   IF ( ANY(deepSOM_a(ip,:,iv,iele) .LT. zero) ) THEN
                      WRITE (numout,*) 'cryoturbate: deepSOM_a<0','iele=',iele, &
                           'ip=',ip,'iv=',iv,'deepSOM_a=',deepSOM_a(ip,:,iv,iele)
                      CALL ipslerr_p (3,'cryoturbate','','','')                            
                   ENDIF
                   IF ( ANY(deepSOM_s(ip,:,iv,iele) .LT. zero) ) THEN
                      WRITE (numout,*) 'cryoturbate: deepSOM_s<0','iele=',iele, &
                           'ip=',ip,'iv=',iv,'deepSOM_s=',deepSOM_s(ip,:,iv,iele)         
                      CALL ipslerr_p (3,'cryoturbate','','','')                            
                   ENDIF
                   IF ( ANY(deepSOM_p(ip,:,iv,iele) .LT. zero) ) THEN
                      WRITE (numout,*) 'cryoturbate: deepSOM_p<0','iele=',iele, &
                           'ip=',ip,'iv=',iv,'deepSOM_p=',deepSOM_p(ip,:,iv,iele)         
                     CALL ipslerr_p (3,'cryoturbate','','','')                            
                   ENDIF

                ENDIF
              ENDDO !End loop over nelements
             ENDDO
          ENDDO


    ELSEIF ( action .EQ. 'coefficients' ) THEN
       IF (firstcall) THEN
          WRITE(*,*) 'error: initilaizations have to happen before coefficients calculated. we stop.'
          CALL ipslerr_p (3,'cryoturbate','Initilaizations have to happen before coefficients calculated','','')
       ENDIF

       cryoturb_location(:,:) =  ( altmax_lastyear(:,:) .LT. max_cryoturb_alt ) &
!In the former vertical discretization scheme the first level was at 0.016 cm; now it's only 0.00048 so we set an equivalent threshold directly as a fixed depth of 1 cm,
            .AND. ( altmax_lastyear(:,:) .GE. min_cryoturb_alt ) .AND. veget_mask_2d(:,:)
       IF (use_fixed_cryoturbation_depth) THEN
          cryoturbation_depth(:,:) = fixed_cryoturbation_depth(:,:)
       ELSE
          cryoturbation_depth(:,:) = altmax_lastyear(:,:)
       ENDIF

       bioturb_location(:,:) = ( ( altmax_lastyear(:,:) .GE. max_cryoturb_alt ) .AND. veget_mask_2d(:,:) )

       DO ip = 1, kjpindex
          DO iv = 1,nvm
             IF ( cryoturb_location(ip,iv) ) THEN
                !
                IF (use_new_cryoturbation) THEN
                   SELECT CASE(cryoturbation_method)
                   CASE(1)
                      !
                      DO il = 1, ngrnd ! linear dropoff to zero between alt and 2*alt
                         IF ( zi_soil(il) .LE. cryoturbation_depth(ip,iv) ) THEN
                            diff_k(ip,il,iv) = diff_k_const
                         ELSE
                            diff_k(ip,il,iv) = diff_k_const*(un-MAX(MIN((zi_soil(il)/cryoturbation_depth(ip,iv))-un,un),zero))
                         ENDIF
                      END DO
                      !
                   CASE(2)
                      !
                      DO il = 1, ngrnd ! exponential dropoff with e-folding distace = alt, below the active layer
                         IF ( zi_soil(il) .LE. cryoturbation_depth(ip,iv) ) THEN
                            diff_k(ip,il,iv) = diff_k_const
                         ELSE
                            diff_k(ip,il,iv) = diff_k_const*(EXP(-MAX((zi_soil(il)/cryoturbation_depth(ip,iv)-un),zero)))
                         ENDIF
                      END DO
                      !
                   CASE(3)
                      !
                      ! exponential dropoff with e-folding distace = alt, starting at surface
                      diff_k(ip,:,iv) = diff_k_const*(EXP(-(zi_soil(:)/cryoturbation_depth(ip,iv))))
                      !
                   CASE(4)
                      !
                      DO il = 1, ngrnd ! linear dropoff to zero between alt and 3*alt
                         IF ( zi_soil(il) .LE. cryoturbation_depth(ip,iv) ) THEN
                            diff_k(ip,il,iv) = diff_k_const
                         ELSE
                            diff_k(ip,il,iv) = diff_k_const*(un-MAX(MIN((zi_soil(il)-cryoturbation_depth(ip,iv))/ &
                                 (2.*cryoturbation_depth(ip,iv)),un),zero))
                         ENDIF
                         IF ( zf_soil(il) .GT. max_cryoturb_alt ) THEN
                            diff_k(ip,il,iv) = zero
                         ENDIF
                      END DO
                      ! 
                      IF (printlev>=3) WRITE(*,*) 'cryoturb method 4: ip, iv, diff_k(ip,:,iv): ', ip, iv, diff_k(ip,:,iv)
                   CASE(5)
                      !
                      DO il = 1, ngrnd ! linear dropoff to zero between alt and 3m
                         IF ( zi_soil(il) .LE. cryoturbation_depth(ip,iv) ) THEN
                            diff_k(ip,il,iv) = diff_k_const
                         ELSE
                            diff_k(ip,il,iv) = diff_k_const*(un-MAX(MIN((zi_soil(il)-cryoturbation_depth(ip,iv))/ &
                                 (3.-cryoturbation_depth(ip,iv)),un),zero))
                         ENDIF
                      END DO
                      !
                      IF (printlev>=3) WRITE(*,*) 'cryoturb method 5: ip, iv, diff_k(ip,:,iv): ', ip, iv, diff_k(ip,:,iv)
                   END SELECT
                   
                   ELSE ! old cryoturbation scheme 
                   !
                   diff_k(ip,1:altmax_ind(ip,iv),iv) = diff_k_const
                   diff_k(ip, altmax_ind(ip,iv)+1,iv) = diff_k_const/10.
                   diff_k(ip, altmax_ind(ip,iv)+2,iv) = diff_k_const/100.
                   diff_k(ip,(altmax_ind(ip,iv)+3):ngrnd,iv) = zero
                ENDIF
             ELSE IF ( bioturb_location(ip,iv) ) THEN
                DO il = 1, ngrnd
                   IF ( zi_soil(il) .LE. bioturbation_depth ) THEN
                      diff_k(ip,il,iv) = bio_diff_k_const
                   ELSE
                      diff_k(ip,il,iv) = zero
                   ENDIF
                END DO
	     ELSE
                diff_k(ip,:,iv) = zero
             END IF
          END DO
       END DO

       mu_soil_rev=diff_k(:,1,:)*time_step/(zf_soil(1)-zf_soil(0))/(zi_soil(2)-zi_soil(1))
       
       DO il = 1,ngrnd-1
          WHERE ( cryoturb_location(:,:) .OR. bioturb_location(:,:) )
             xc_cryoturb(:,il,:) = (zf_soil(il)-zf_soil(il-1))  / time_step
             xd_cryoturb(:,il,:) = diff_k(:,il,:) / (zi_soil(il+1)-zi_soil(il))
          endwhere
       ENDDO


       DO iele = 1,nelements        
          WHERE ( cryoturb_location(:,:) .OR. bioturb_location(:,:)  )
             xc_cryoturb(:,ngrnd,:) = (zf_soil(ngrnd)-zf_soil(ngrnd-1))  / time_step
          
             !bottom
             xe_a(:,:) = xc_cryoturb(:,ngrnd,:)+xd_cryoturb(:,ngrnd-1,:)
             xe_s(:,:) = xc_cryoturb(:,ngrnd,:)+xd_cryoturb(:,ngrnd-1,:)
             xe_p(:,:) = xc_cryoturb(:,ngrnd,:)+xd_cryoturb(:,ngrnd-1,:)
             alpha_a(:,ngrnd-1,:,iele) = xd_cryoturb(:,ngrnd-1,:) / xe_a(:,:)
             alpha_s(:,ngrnd-1,:,iele) = xd_cryoturb(:,ngrnd-1,:) / xe_s(:,:)
             alpha_p(:,ngrnd-1,:,iele) = xd_cryoturb(:,ngrnd-1,:) / xe_p(:,:)
             beta_a(:,ngrnd-1,:,iele) = xc_cryoturb(:,ngrnd,:)*deepSOM_a(:,ngrnd,:,iele) / xe_a(:,:)
             beta_s(:,ngrnd-1,:,iele) = xc_cryoturb(:,ngrnd,:)*deepSOM_s(:,ngrnd,:,iele) / xe_s(:,:)
             beta_p(:,ngrnd-1,:,iele) = xc_cryoturb(:,ngrnd,:)*deepSOM_p(:,ngrnd,:,iele) / xe_p(:,:)
          END WHERE

          !other levels
          DO il = ngrnd-2,1,-1
             WHERE ( cryoturb_location(:,:) .OR. bioturb_location(:,:) )
                xe_a(:,:) = xc_cryoturb(:,il+1,:) + (1.-alpha_a(:,il+1,:,iele))*xd_cryoturb(:,il+1,:) + xd_cryoturb(:,il,:)
                xe_s(:,:) = xc_cryoturb(:,il+1,:) + (1.-alpha_s(:,il+1,:,iele))*xd_cryoturb(:,il+1,:) + xd_cryoturb(:,il,:)
                xe_p(:,:) = xc_cryoturb(:,il+1,:) + (1.-alpha_p(:,il+1,:,iele))*xd_cryoturb(:,il+1,:) + xd_cryoturb(:,il,:)
                alpha_a(:,il,:,iele) = xd_cryoturb(:,il,:) / xe_a(:,:)
                alpha_s(:,il,:,iele) = xd_cryoturb(:,il,:) / xe_s(:,:)
                alpha_p(:,il,:,iele) = xd_cryoturb(:,il,:) / xe_p(:,:)
                beta_a(:,il,:,iele) = (xc_cryoturb(:,il+1,:)*deepSOM_a(:,il+1,:,iele) + &
                     xd_cryoturb(:,il+1,:)*beta_a(:,il+1,:,iele)) / xe_a(:,:)
                beta_s(:,il,:,iele) = (xc_cryoturb(:,il+1,:)*deepSOM_s(:,il+1,:,iele) + &
                     xd_cryoturb(:,il+1,:)*beta_s(:,il+1,:,iele)) / xe_s(:,:)
                beta_p(:,il,:,iele) = (xc_cryoturb(:,il+1,:)*deepSOM_p(:,il+1,:,iele) + &
                     xd_cryoturb(:,il+1,:)*beta_p(:,il+1,:,iele)) / xe_p(:,:)
             END WHERE
          ENDDO
       ENDDO !End lop over nelements
    ELSE
       !
       ! do not know this action
       !
       CALL ipslerr_p(3, 'cryoturbate', 'DO NOT KNOW WHAT TO DO:', TRIM(action), '')
       !
    ENDIF
    
    ! keep last action in mind
    !
    last_action = action

    !! Check mass balance closure
    IF (err_act.GT.1) THEN

       !! 4.1.1 Calculate components of the mass balance
       pool_end(:,:,:) = zero 
       DO iele = 1,nelements
          DO igrnd = 1,ngrnd
             pool_end(:,:,iele) = pool_end(:,:,iele) + &
                  (deepSOM_a(:,igrnd,:,iele) + deepSOM_s(:,igrnd,:,iele) + deepSOM_p(:,igrnd,:,iele)) * &
                  (zf_soil(igrnd)-zf_soil(igrnd-1)) * veget_max(:,:)
          END DO
       ENDDO
       
       !! 4.1.2 Calculate mass balance
       check_intern(:,:,:,:) = zero
       check_intern(:,:,iatm2land,:) = zero
       check_intern(:,:,iland2atm,icarbon) = -un * zero
       check_intern(:,:,iland2atm,initrogen) = -un * zero
       check_intern(:,:,ilat2out,:) = -un *zero
       check_intern(:,:,ilat2in,:) =  zero
       check_intern(:,:,ipoolchange,:) = -un * (pool_end(:,:,:) - &
            pool_start(:,:,:))
       closure_intern = zero
       DO imbc = 1,nmbcomp
          DO iele = 1,nelements
             closure_intern(:,:,iele) = closure_intern(:,:,iele) + &
                  check_intern(:,:,imbc,iele)
          END DO
       ENDDO

       CALL check_mass_balance("cryoturbate", closure_intern, &
            kjpindex, pool_end, pool_start, veget_max, 'pft')

    END IF ! err_act.GT.1
    
  END  SUBROUTINE cryoturbate

!!
!================================================================================================================================
!! SUBROUTINE   : permafrost_decomp
!!
!>\BRIEF        This routine calculates soil organic matter decomposition
!! DESCRIPTION : 
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!! 
!! REFERENCE(S) : None
!!
!! FLOWCHART11    : None
!! \n
!_
!================================================================================================================================     

  SUBROUTINE permafrost_decomp (kjpindex, time_step, tprof, airvol_soil, &
       oxlim, tau_CH4troph, ok_methane, fbactratio, O2m, &
       totporO2_soil, totporCH4_soil, hslong, clay, silt, &
       deepSOM_a, deepSOM_s, deepSOM_p, &
       deltaCH4g, deltaCH4, deltaSOM1_a, deltaSOM1_s, deltaSOM1_p, deltaSOM2, &
       deltaSOM3, deltaSOM1_a_pftmean, deltaSOM1_s_pftmean, deltaSOM1_p_pftmean, & 
       deltaSOM1_a_pftmean_stock, deltaSOM1_s_pftmean_stock, deltaSOM1_p_pftmean_stock, &
       O2_soil, CH4_soil, fbact_out, CN_target, MG_useallCpools, veget_max, veget_max_bg,&
       fluxout_a,fluxout_s,fluxout_p,somflux, &
       deepSOM_pt,deepSOM_peat,peat_OLT, soiltile)

  !! 0. Variable and parameter declaration

    !! 0.1  Input variables     

    INTEGER(i_std), INTENT(in)			    :: kjpindex        !! domain size
    REAL(r_std), DIMENSION(:,:),INTENT(in)          :: veget_max       !! maximum vegetation fraction
    REAL(r_std),DIMENSION(kjpindex,nvm),INTENT(in)  :: veget_max_bg    
    REAL(r_std), INTENT(in)                         :: time_step       !! time step in seconds
    REAL(r_std), DIMENSION(:,:,:),   INTENT(in)     :: tprof           !! deep temperature profile
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)       :: airvol_soil
    LOGICAL, INTENT(in)                             :: oxlim           !! O2 limitation taken into account
    REAL(r_std), INTENT(in)                         :: tau_CH4troph    !! time constant of methanetrophy (s)
    LOGICAL, INTENT(in)                             :: ok_methane      !! Is Methanogenesis and -trophy taken into account? 
    REAL(r_std), INTENT(in)                         :: fbactratio      !! time constant of methanogenesis (ratio to that of oxic)
    REAL(r_std), INTENT(in)                         :: O2m             !! oxygen concentration [g/m3] below which there is anoxy 
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)       :: totporO2_soil   !! total O2 porosity (Tans, 1998)
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)       :: totporCH4_soil  !! total CH4 porosity (Tans, 1998)
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)       :: hslong          !! deep soil humidity
    REAL(r_std), DIMENSION(:), INTENT(in)           :: clay            !! clay content
    REAL(r_std), DIMENSION(:), INTENT(in)           :: silt            !! silt content
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)       :: fbact_out
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)       :: CN_target       !! C to N ratio of SOM flux from one pool to another (gN m-2 dt-1)  
    LOGICAL, INTENT(in)                             :: MG_useallCpools !! Do we allow all three C pools to feed methanogenesis?
    REAL(r_std),DIMENSION (:,:), INTENT(in)         :: soiltile        !! Fraction of soil tile within vegtot (0-1, unitless)

    !! 0.2 Output variables
                                                                                                       !! Analytic spinup  (gC or N m-2 dt-1)
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm,nelements), INTENT(out)             :: fluxout_a         !! fluxes leaving the active pool (gC or N m-2 dt-1)
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm,nelements), INTENT(out)             :: fluxout_s         !! fluxes leaving the slow pool  (gC or N m-2 dt-1)
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm,nelements), INTENT(out)             :: fluxout_p         !! fluxes leaving the passive pool (gC or N m-2 dt-1)
    REAL(r_std), DIMENSION(ncarb,ncarb,ngrnd,nelements), INTENT(out)              :: somflux           !! fluxes between soil orgnanic matter reservoirs

    !! 0.3 Modified variables

    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm,nelements), INTENT(inout)  :: deepSOM_a         !! soil organic matter (g/m**3) active
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm,nelements), INTENT(inout)  :: deepSOM_s         !! soil organic matter (g/m**3) slow
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm,nelements), INTENT(inout)  :: deepSOM_p         !! soil organic matter (g/m**3) passive
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm), INTENT(inout)            :: deltaCH4
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm), INTENT(inout)            :: deltaCH4g
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm,nelements), INTENT(inout)  :: deltaSOM1_a
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm,nelements), INTENT(inout)  :: deltaSOM1_s
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm,nelements), INTENT(inout)  :: deltaSOM1_p

    REAL(r_std), DIMENSION(kjpindex,ngrnd), INTENT(inout)                :: deltaSOM1_a_pftmean
    REAL(r_std), DIMENSION(kjpindex,ngrnd), INTENT(inout)                :: deltaSOM1_s_pftmean
    REAL(r_std), DIMENSION(kjpindex,ngrnd), INTENT(inout)                :: deltaSOM1_p_pftmean
    REAL(r_std), DIMENSION(kjpindex), INTENT(inout)                      :: deltaSOM1_a_pftmean_stock
    REAL(r_std), DIMENSION(kjpindex), INTENT(inout)                      :: deltaSOM1_s_pftmean_stock
    REAL(r_std), DIMENSION(kjpindex), INTENT(inout)                      :: deltaSOM1_p_pftmean_stock
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm,nelements), INTENT(inout)  :: deltaSOM2
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm,nelements), INTENT(inout)  :: deltaSOM3
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm), INTENT(inout)            :: O2_soil         !! oxygen (g O2/m**3 air)
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm), INTENT(inout)            :: CH4_soil        !! methane (g CH4/m**3 air)

    !peatland
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm,nelements), INTENT(out)    :: deepSOM_pt      !! peat carbon concentration (g/m3)
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm,nelements), INTENT(inout)  :: deepSOM_peat    !! soil carbon for peat
    REAL(r_std), DIMENSION(kjpindex,nvm,nelements),INTENT(out)           :: peat_OLT        !! peat organic layer thickness (in m)
    REAL(r_std), DIMENSION(ngrnd)                                        :: peat_BD         !! bulk density of the soil
    REAL(r_std), DIMENSION(ngrnd)                                        :: peat_SOC        !! soil carbon concentration of the soil
    REAL(r_std), ALLOCATABLE, DIMENSION(:),SAVE                          :: SOMmax          !! maximum allowed carbon content at each soil layer
!$OMP THREADPRIVATE(SOMmax)
    REAL(r_std)                                                          :: excessSOM       !! the excess of soil carbon to be transferred to lower layer
    REAL(r_std)                                                          :: max_vsreal      !! fraction of soil carbon after remove the excesss
    REAL(r_std)                                                          :: trans_flux      !! soil carbon transferred from the excess
    REAL(r_std)                                                          :: SOMthick        !! the thickness of soil carbon
    
    !! 0.4 Local variables

    INTEGER(i_std)                                             :: ier, istm
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm)                 :: nadd_soil       !! number of moles created / m**3 of air
    REAL(r_std)                                                :: fbact_a,fbact_s, fbact_p,temp
    REAL(r_std)                                                :: fbactCH4_a, fbactCH4_s, fbactCH4_p
    REAL(r_std), DIMENSION(nelements)                          :: dC
    REAL(r_std)                                                :: dCm
    REAL(r_std)                                                :: dCH4,dCH4m,dO2
    INTEGER(i_std)                                             :: il, ip, iv, iele
    LOGICAL, SAVE                                              :: firstcall = .TRUE.        !! first call?
!$OMP THREADPRIVATE(firstcall)
    INTEGER(i_std)                                             :: imbc, igrnd      !! Indices (unitless)
    INTEGER(i_std)                                             :: inbpools, icarb  !! Indices (unitless)
    REAL(r_std), DIMENSION(kjpindex,nvm,nmbcomp,nelements)     :: check_intern     !! Contains the components of the internal
                                                                                   !! mass balance check for this routine   
                                                                                   !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(kjpindex,nvm,nelements)             :: closure_intern   !! Check closure of internal mass balance
                                                                                   !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(kjpindex,nvm,nelements)             :: pool_start       !! Start and end pool of this routine 
                                                                                   !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(kjpindex,nvm,nelements)             :: pool_end         !! Start and end pool of this routine 
                                                                                   !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex

!_ ================================================================================================================================

    IF (firstcall) THEN

       ALLOCATE (fc(kjpindex,3,3,nvm),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'permafrost_decomp', 'Pb in alloc for fc','','')

       ALLOCATE (fr(kjpindex,3,nvm),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'permafrost_decomp', 'Pb in alloc for fr','','')

       ALLOCATE (frac_carb(kjpindex,ncarb,ncarb), stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'permafrost_decomp', 'Pb in alloc for frac_carb','','')

       !! peatland
       ALLOCATE (SOMmax(ngrnd),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'permafrost_decomp', 'Pb in alloc for SOMmax','','')
 
       !
       ! calculate soil organic matter flux fractions
       !
       DO iv =1,nvm
          fc(:,iactive,iactive,iv) = 0.0_r_std
          fc(:,iactive,ipassive,iv) = 0.004_r_std
          fc(:,iactive,islow,iv) = 1._r_std - (.85-.68*clay(:)) - fc(:,iactive,ipassive,iv)
          !
          fc(:,islow,islow,iv) = .0_r_std
          fc(:,islow,iactive,iv) = .42_r_std
          fc(:,islow,ipassive,iv) = .03_r_std
          !
          fc(:,ipassive,ipassive,iv) = .0_r_std
          fc(:,ipassive,iactive,iv) = .45_r_std
          fc(:,ipassive,islow,iv) = .0_r_std
          !
          fr(:,:,iv) = 1._r_std-fc(:,:,iactive,iv)-fc(:,:,islow,iv)-fc(:,:,ipassive,iv)
          firstcall = .FALSE.
       END DO
!BG to be tested the same parameters than the soil without discretization
!       !
!       ! calculate soil organic matter flux fractions
!       !
!       DO iv =1,nvm
!          fc(:,iactive,iactive,iv) = 0.0_r_std
!          fc(:,iactive,ipassive,iv) = active_to_pass_ref_frac + active_to_pass_clay_frac*clay(:)
!          fc(:,iactive,islow,iv) = un - frac_carb(:,iactive,ipassive) - (active_to_co2_ref_frac - &
!             active_to_co2_clay_silt_frac*(clay(:)+silt(:)))
!          !
!          fc(:,islow,islow,iv) = .0_r_std
!          fc(:,islow,ipassive,iv) = slow_to_pass_ref_frac + slow_to_pass_clay_frac*clay(:)
!          fc(:,islow,iactive,iv) = un - frac_carb(:,islow,ipassive) - slow_to_co2_ref_frac
!          !
!          fc(:,ipassive,ipassive,iv) = .0_r_std
!          fc(:,ipassive,iactive,iv) = pass_to_active_ref_frac
!          fc(:,ipassive,islow,iv) = pass_to_slow_ref_frac
!          !
!          fr(:,:,iv) = un-fc(:,:,iactive,iv)-fc(:,:,islow,iv)-fc(:,:,ipassive,iv)
!          firstcall = .FALSE.
!       END DO

       IF (printlev>=3) THEN
          DO ip = 1,kjpindex
             WRITE(*,*) 'cdk: permafrost_decomp: i, fraction respired gridcell(i) :', ip, fr(ip,:,1)
          END DO
       ENDIF
    ENDIF

    ! bulk density of the soil
    peat_BD(:) = peat_bulk_density(:)   !in g/cm3
    ! soil organic carbon concentration
    peat_SOC(:)=1._r_std/((0.4*peat_BD(:)+0.13)**2.19)
    peat_SOC(:) = peat_SOC(:)*0.01    
    ! the maximum allowed carbon content per soil layer 
    ! zf_soil,zi_soil: in m
    ! peat_BD in g/cm3  
    DO il=1,ngrnd
       SOMmax(il) =  peat_BD(il)*1.E6*peat_SOC(il)*(zf_soil(il)-zf_soil(il-1))   !SOMmax in g/m**2
    ENDDO

    !! Initialize check for mass balance closure
    IF (err_act.GT.1) THEN

       pool_start(:,:,:) = zero
       
       DO iele = 1,nelements
          DO igrnd = 1,ngrnd
             pool_start(:,:,iele) = pool_start(:,:,iele) + &
                  (deepSOM_a(:,igrnd,:,iele) + deepSOM_s(:,igrnd,:,iele) + deepSOM_p(:,igrnd,:,iele)) * &
                  (zf_soil(igrnd)-zf_soil(igrnd-1)) * veget_max(:,:)
          END DO
       END DO

    END IF ! err_act.gt.1

    !
    ! calculate som consumption
    !
    nadd_soil(:,:,:) = zero
    somflux(:,:,:,:) = zero

    deltaSOM1_a(:,:,:,:) = zero
    deltaSOM1_s(:,:,:,:) = zero
    deltaSOM1_p(:,:,:,:) = zero
    deltaCH4(:,:,:) = zero
    deltaCH4g(:,:,:) = zero
    deltaSOM2(:,:,:,:) = zero
    deltaSOM3(:,:,:,:) = zero   
    fluxout_a(:,:,:,:) = zero
    fluxout_s(:,:,:,:) = zero
    fluxout_p(:,:,:,:) = zero

    DO ip = 1, kjpindex
       !
       DO iv = 1, nvm
          !
          IF (  veget_mask_2d(ip,iv) ) THEN
             !
             DO il = 1, ngrnd
                !
                ! 1 function that gives soil organic matter residence time as a function of
                !     soil temperature (in seconds)
                !
                temp = tprof(ip,il,iv) - ZeroCelsius
                fbact_a = fbact_out(ip,il,iv)
                fbact_a = MAX(fbact_a,time_step)
                !
                IF ( fbact_a/HUGE(1.) .GT. .1 ) THEN
                   fbact_s = fbact_a
                   fbact_p = fbact_a
                ELSE
                   fbact_s = fbact_a * fslow
                   fbact_p = fbact_a * fpassive
                ENDIF
                !
                ! methanogenesis: first guess, 10 times (fbactratio) slower than oxic
                ! decomposition
                IF ( fbact_a/HUGE(1.) .GT. .1 ) THEN 
                   fbactCH4_a = fbact_a
                   fbactCH4_s = fbact_s
                   fbactCH4_p = fbact_p
                ELSE
                   fbactCH4_a = fbact_a * fbactratio
                   IF ( MG_useallCpools ) THEN
                      fbactCH4_s = fbact_s * fbactratio
                      fbactCH4_p = fbact_p * fbactratio
                   ELSE
                      fbactCH4_s = HUGE(1.0)
                      fbactCH4_p = HUGE(1.0)
                   ENDIF
                ENDIF
                !
                ! 2 oxic decomposition: carbon and oxygen consumption
                !
                ! 2.1 active
                !
             DO iele = 1, nelements

                IF (oxlim) THEN
                   dCm = O2_soil(ip,il,iv)*airvol_soil(ip,il,iv)*wC/wO2
                   dC(iele) = MIN(deepSOM_a(ip,il,iv,iele) * time_step/fbact_a,dCm)
                ELSE
                   dC(iele) = deepSOM_a(ip,il,iv,iele) * time_step/fbact_a
                ENDIF

                ! pour actif
                dC(iele) = dC(iele) * ( un - som_turn_iactive_clay_frac * clay(ip) )

                ! keep track of the decomposition controlling the flux  out for Analytcial spinup
                IF (spinup_analytic) THEN
                    fluxout_a(ip,il,iv,iele) = time_step/fbact_a
                ENDIF

                ! flux vers les autres reservoirs
                IF (iele .EQ. icarbon) THEN
                somflux(iactive,ipassive,il,iele) = fc(ip,iactive,ipassive,iv) * dC(iele) 
                somflux(iactive,islow,il,iele) = fc(ip,iactive,islow,iv) * dC(iele) 
                ELSEIF (iele .EQ. initrogen) THEN
                somflux(iactive,ipassive,il,iele) = fc(ip,iactive,ipassive,iv) * dC(icarbon) / CN_target(ip,iv, ipassive) 
                somflux(iactive,islow,il,iele) = fc(ip,iactive,islow,iv) * dC(icarbon) / CN_target(ip,iv, islow)
                ELSE
                   IF (ier /= 0) CALL ipslerr_p(3,'permafrost_decomp', 'We we have an elements which is neither carbon nor nitrogen','','')
                ENDIF
                !
                deepSOM_a(ip,il,iv,iele) = deepSOM_a(ip,il,iv,iele) - dC(iele)
                dO2 = wO2/wC *dC(icarbon)*fr(ip,iactive,iv) / totporO2_soil(ip,il,iv)
                O2_soil(ip,il,iv) = MAX( O2_soil(ip,il,iv) - dO2, zero)
                ! keep delta C * fr in memory (generates energy)
                IF (iele .EQ. icarbon) THEN
                deltaSOM1_a(ip,il,iv,iele) = dC(iele)*fr(ip,iactive,iv) !!this line!!!
                ELSEIF (iele .EQ. initrogen) THEN
                deltaSOM1_a(ip,il,iv,iele) = dC(iele) - (somflux(iactive,ipassive,il,iele)+somflux(iactive,islow,il,iele))
                ELSE
                   IF (ier /= 0) CALL ipslerr_p(3,'permafrost_decomp', 'We we have an elements which is neither carbon nor nitrogen','','')
                ENDIF
                !
                ! 2.2 slow       
                !
                IF (oxlim) THEN
                   dCm = O2_soil(ip,il,iv)*airvol_soil(ip,il,iv)*wC/wO2
                   dC(iele) = MIN(deepSOM_s(ip,il,iv,iele) * time_step/fbact_s,dCm)
                ELSE
                   dC(iele) = deepSOM_s(ip,il,iv,iele) * time_step/fbact_s
                ENDIF
                ! keep track of the decomposition controlling the flux out for Analytcial spinup
                IF (spinup_analytic) THEN
                    fluxout_s(ip,il,iv,iele) = time_step/fbact_s
                ENDIF

                ! flux vers les autres reservoirs
                IF (iele .EQ. icarbon) THEN
                somflux(islow,iactive,il,iele) = fc(ip,islow,iactive,iv) * dC(iele)
                somflux(islow,ipassive,il,iele) = fc(ip,islow,ipassive,iv) * dC(iele)
                ELSEIF (iele .EQ. initrogen) THEN
                somflux(islow,iactive,il,iele) = fc(ip,islow,iactive,iv) * dC(icarbon) / CN_target(ip,iv, iactive)
                somflux(islow,ipassive,il,iele) = fc(ip,islow,ipassive,iv) * dC(icarbon) / CN_target(ip,iv, ipassive)
                ELSE
                   IF (ier /= 0) CALL ipslerr_p(3,'permafrost_decomp', 'We we have an elements which is neither carbon nor nitrogen','','')
                ENDIF
                !
                deepSOM_s(ip,il,iv,iele) = deepSOM_s(ip,il,iv,iele) - dC(iele)
                dO2 = wO2/wC * dC(iele)*fr(ip,islow,iv) / totporO2_soil(ip,il,iv)
                O2_soil(ip,il,iv) = MAX( O2_soil(ip,il,iv) - dO2, zero)
                ! keep delta C * fr in memory (generates energy)
                IF (iele .EQ. icarbon) THEN
                deltaSOM1_s(ip,il,iv,iele) = dC(iele)*fr(ip,islow,iv) !!this line!!!
                ELSEIF (iele .EQ. initrogen) THEN
                deltaSOM1_s(ip,il,iv,iele) = dC(iele) - (somflux(islow,iactive,il,iele)+somflux(islow,ipassive,il,iele))
                ELSE
                   IF (ier /= 0) CALL ipslerr_p(3,'permafrost_decomp', 'We we have an elements which is neither carbon nor nitrogen','','')
                ENDIF
                !
                ! 2.3 passive
                !
                IF (oxlim) THEN
                   dCm = O2_soil(ip,il,iv)*airvol_soil(ip,il,iv)*wC/wO2
                   dC(iele) = MIN(deepSOM_p(ip,il,iv,iele) * time_step/fbact_p,dCm)
                ELSE
                   dC(iele) = deepSOM_p(ip,il,iv,iele) * time_step/fbact_p
                ENDIF
                ! keep track of the decomposition controlling the flux out for Analytcial spinup
                IF (spinup_analytic) THEN 
                    fluxout_p(ip,il,iv,iele) = time_step/fbact_p
                ENDIF

                ! flux vers les autres reservoirs
                IF (iele .EQ. icarbon) THEN
                somflux(ipassive,iactive,il,iele) = fc(ip,ipassive,iactive,iv) * dC(iele)
                somflux(ipassive,islow,il,iele) = fc(ip,ipassive,islow,iv) * dC(iele)
                ELSEIF (iele .EQ. initrogen) THEN
                somflux(ipassive,iactive,il,iele) = fc(ip,ipassive,iactive,iv) * dC(icarbon) / CN_target(ip,iv, iactive)
                somflux(ipassive,islow,il,iele) = fc(ip,ipassive,islow,iv) * dC(icarbon) / CN_target(ip,iv, islow)
                ELSE
                   IF (ier /= 0) CALL ipslerr_p(3,'permafrost_decomp', 'We we have an elements which is neither carbon nor nitrogen','','')
                ENDIF
                !
                deepSOM_p(ip,il,iv,iele) = deepSOM_p(ip,il,iv,iele) - dC(iele)
                dO2 = wO2/wC * dC(iele)*fr(ip,ipassive,iv) / totporO2_soil(ip,il,iv)
                O2_soil(ip,il,iv) = MAX( O2_soil(ip,il,iv) - dO2, zero)
                ! keep delta C * fr in memory (generates energy)
                IF (iele .EQ. icarbon) THEN
                deltaSOM1_p(ip,il,iv,iele) = dC(iele)*fr(ip,ipassive,iv) !!this line!!!
                ELSEIF (iele .EQ. initrogen) THEN
                deltaSOM1_p(ip,il,iv,iele) = dC(iele)- (somflux(ipassive,iactive,il,iele)+somflux(ipassive,islow,il,iele))
                ELSE
                   IF (ier /= 0) CALL ipslerr_p(3,'permafrost_decomp', 'We we have an elements which is neither carbon nor nitrogen','','')
                ENDIF
                !
                !
                ! 3 methanogenesis or methanotrophy
                !   
                !
                IF ((ok_methane) .AND. (iele .EQ. icarbon)) THEN
                   !
                   !
                   ! 3.1 active pool methanogenesis
                   dC(iele) = deepSOM_a(ip,il,iv,iele) * time_step / fbactCH4_a * EXP(-O2_soil(ip,il,iv)*(1+hslong(ip,il,iv) * &
                        (BunsenO2-1.)) / O2m ) !DKtest: when commented, no ox lim for MG
                   ! pour actif
                   dC(iele) = dC(iele) * ( 1. - .75 * clay(ip) )
                   dCH4 = dC(iele)*fr(ip,iactive,iv) * wCH4/wC / totporCH4_soil(ip,il,iv)
                   !
                   !
                   ! flux vers les autres reservoirs
                   somflux(iactive,ipassive,il,iele)=somflux(iactive,ipassive,il,iele)+fc(ip,iactive,ipassive,iv)*dC(iele)
                   somflux(iactive,islow,il,iele)=somflux(iactive,islow,il,iele)+fc(ip,iactive,islow,iv)*dC(iele)
                   !
                   deepSOM_a(ip,il,iv,iele) = deepSOM_a(ip,il,iv,iele) - dC(iele)
                   !
                   deltaCH4g(ip,il,iv) = dCH4
                   !
                   CH4_soil(ip,il,iv) = CH4_soil(ip,il,iv) + dCH4
                   ! keep delta C*fr in memory (generates energy)
                   deltaSOM2(ip,il,iv,iele) = dC(iele)*fr(ip,iactive,iv)
                   !
                   ! how many moles of gas / m**3 of air did we generate?
                   ! (methanogenesis generates 1 molecule net if we take
                   !  B -> B' + CH4 ) 
                   nadd_soil(ip,il,iv) = nadd_soil(ip,il,iv) + dCH4/wCH4
                   !
                   !
                   IF ( MG_useallCpools ) THEN
                      !
                      ! 3.2 slow pool methanogenesis  cdk: adding this to allow other carbon pools to participate in MG
                      dC(iele) = deepSOM_s(ip,il,iv,iele) * time_step / fbactCH4_s * EXP(-O2_soil(ip,il,iv)*(1+hslong(ip,il,iv) * &
                           (BunsenO2-1.)) / O2m ) !DKtest: when commented, no ox lim for MG
                      dCH4 = dC(iele)*fr(ip,islow,iv) * wCH4/wC / totporCH4_soil(ip,il,iv)
                      !
                      ! flux vers les autres reservoirs
                      somflux(islow,ipassive,il,iele)=somflux(islow,ipassive,il,iele)+fc(ip,islow,ipassive,iv)*dC(iele)
                      somflux(islow,iactive,il,iele)=somflux(islow,iactive,il,iele)+fc(ip,islow,iactive,iv)*dC(iele)
                      !
                      deepSOM_s(ip,il,iv,iele) = deepSOM_s(ip,il,iv,iele) - dC(iele)
                      !
                      deltaCH4g(ip,il,iv) = deltaCH4g(ip,il,iv) + dCH4
                      CH4_soil(ip,il,iv) = CH4_soil(ip,il,iv) + dCH4
                      ! keep delta C*fr in memory (generates energy)
                      deltaSOM2(ip,il,iv,iele) = deltaSOM2(ip,il,iv,iele) + dC(iele)*fr(ip,islow,iv)
                      !
                      ! how many moles of gas / m**3 of air did we generate?
                      ! (methanogenesis generates 1 molecule net if we take
                      !  B -> B' + CH4 ) 
                      nadd_soil(ip,il,iv) = nadd_soil(ip,il,iv) + dCH4/wCH4
                      !	      
                      !
                      !
                      ! 3.3 passive pool methanogenesis  cdk: adding this to allow other carbon pools to participate in MG
                      dC(iele) = deepSOM_p(ip,il,iv,iele) * time_step / fbactCH4_p * EXP(-O2_soil(ip,il,iv)*(1+hslong(ip,il,iv) * &
                          (BunsenO2-1.)) / O2m ) !DKtest: when commented, no ox lim for MG
                      dCH4 = dC(iele)*fr(ip,ipassive,iv) * wCH4/wC / totporCH4_soil(ip,il,iv)
                      !
                      ! flux vers les autres reservoirs
                      somflux(ipassive,islow,il,iele)=somflux(ipassive,islow,il,iele)+fc(ip,ipassive,islow,iv)*dC(iele)
                      somflux(ipassive,iactive,il,iele)=somflux(ipassive,iactive,il,iele)+fc(ip,ipassive,iactive,iv)*dC(iele)
                      !
                      deepSOM_p(ip,il,iv,iele) = deepSOM_p(ip,il,iv,iele) - dC(iele)
                      !
                      deltaCH4g(ip,il,iv) = deltaCH4g(ip,il,iv) + dCH4
                      CH4_soil(ip,il,iv) = CH4_soil(ip,il,iv) + dCH4
                      ! keep delta C*fr in memory (generates energy)
                      deltaSOM2(ip,il,iv,iele) = deltaSOM2(ip,il,iv,iele) + dC(iele)*fr(ip,ipassive,iv)
                      !
                      ! how many moles of gas / m**3 of air did we generate?
                      ! (methanogenesis generates 1 molecule net if we take
                      !  B -> B' + CH4 ) 
                      nadd_soil(ip,il,iv) = nadd_soil(ip,il,iv) + dCH4/wCH4
                      !	      
                      !
                   ENDIF
                   !
                   ! trophy: 
                   ! no temperature dependence except that T>0��C (Price et
                   ! al, GCB 2003; Koschorrek and Conrad, GBC 1993).
                   ! tau_CH4troph is such that we fall between values of
                   ! soil methane oxidation flux given by these authors.
                   ! 
                   IF ( temp .GE. zero ) THEN
                      !
                      dCH4m = O2_soil(ip,il,iv)/2. * wCH4/wO2 * totporO2_soil(ip,il,iv)/totporCH4_soil(ip,il,iv)
                      !      		dCH4m = CH4_soil(ip,il,iv)  !DKtest - no ox lim to trophy
                      dCH4 = MIN( CH4_soil(ip,il,iv) * time_step/MAX(tau_CH4troph,time_step), dCH4m )
                      CH4_soil(ip,il,iv) = CH4_soil(ip,il,iv) - dCH4
                      dO2 = 2.*dCH4 * wO2/wCH4 * totporCH4_soil(ip,il,iv)/totporO2_soil(ip,il,iv)
                      O2_soil(ip,il,iv) = MAX( O2_soil(ip,il,iv) - dO2, zero)
                      ! keep delta CH4 in memory (generates energy)
                      deltaCH4(ip,il,iv) = dCH4
                      ! carbon (g/m3 soil) transformed to CO2
                      deltaSOM3(ip,il,iv,icarbon)=dCH4/wCH4*wC*totporCH4_soil(ip,il,iv)
                      ! how many moles of gas / m**3 of air did we generate?
                      ! (methanotrophy consumes 2 molecules net if we take
                      !  CH4 + 2 O2 -> CO2 + 2 H2O )
                      nadd_soil(ip,il,iv) = nadd_soil(ip,il,iv)-2.*dCH4/wCH4
                      !
                   ENDIF
                   
                ENDIF
               
                ! 4 add fluxes between reservoirs
                
                deepSOM_a(ip,il,iv,iele)=deepSOM_a(ip,il,iv,iele)+somflux(islow,iactive,il,iele)+somflux(ipassive,iactive,il,iele)
                deepSOM_s(ip,il,iv,iele)=deepSOM_s(ip,il,iv,iele)+somflux(iactive,islow,il,iele)+somflux(ipassive,islow,il,iele)
                deepSOM_p(ip,il,iv,iele)=deepSOM_p(ip,il,iv,iele)+somflux(iactive,ipassive,il,iele)+somflux(islow,ipassive,il,iele)
          
             ENDDO
             
            ENDDO ! End loop over nelements
            
          ENDIF

       ENDDO
       
    ENDDO

    IF (perma_peat) THEN
       deepSOM_pt(:,:,:,:)=zero
       deepSOM_peat(:,:,:,:)=zero
    ENDIF
    
    IF (perma_peat) THEN
       DO ip = 1, kjpindex
          DO il = 1, ngrnd
             DO iv = 1, nvm
                istm = pref_soil_veg(iv)
                DO iele = 1, nelements
                   IF (is_peat(iv) .AND. veget_mask_2d(ip,iv) .AND. soiltile(ip,istm) .GT. zero) THEN
                      !!total carbon in each layer (sum of active,slow,passive)
                      deepSOM_peat(ip,il,iv,iele)=deepSOM_a(ip,il,iv, iele)+deepSOM_s(ip,il,iv,iele)+deepSOM_p(ip,il,iv,iele) !g/m^3
                      deepSOM_peat(ip,il,iv,iele)=deepSOM_peat(ip,il,iv,iele)*(zf_soil(il)-zf_soil(il-1))!g/m^2
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDIF

    ! compare deepSOM_peat with SOMmax, the excess will be transferred to lower layer
    IF (perma_peat) THEN
       DO ip = 1, kjpindex
          DO il = 1, ngrnd-1 
             DO iv = 1, nvm
                istm = pref_soil_veg(iv)
                DO iele = 1, nelements
                   IF (is_peat(iv) .AND. veget_mask_2d(ip,iv) .AND. soiltile(ip,istm) .GT. zero) THEN
                      IF (deepSOM_peat(ip,il,iv,iele) .GT. frac1*SOMmax(il)) THEN
                         excessSOM=deepSOM_peat(ip,il,iv,iele)*frac2
                         max_vsreal= (deepSOM_peat(ip,il,iv,iele)-excessSOM)/deepSOM_peat(ip,il,iv,iele)
                         deepSOM_peat(ip,il,iv,iele)=deepSOM_peat(ip,il,iv,iele)-excessSOM
                         deepSOM_a(ip,il,iv,iele)= deepSOM_a(ip,il,iv,iele)* max_vsreal
                         
                         deepSOM_s(ip,il,iv,iele)= deepSOM_s(ip,il,iv,iele)* max_vsreal
                         deepSOM_p(ip,il,iv,iele)= deepSOM_p(ip,il,iv,iele)* max_vsreal
                         
                         trans_flux=excessSOM*(zf_soil(il)-zf_soil(il-1))/(zf_soil(il+1)-zf_soil(il))
                         deepSOM_a(ip,il+1,iv,iele)= deepSOM_a(ip,il+1,iv,iele)+(deepSOM_a(ip,il,iv,iele)/deepSOM_peat(ip,il,iv,iele))*trans_flux
                         deepSOM_s(ip,il+1,iv,iele)= deepSOM_s(ip,il+1,iv,iele)+(deepSOM_s(ip,il,iv,iele)/deepSOM_peat(ip,il,iv,iele))*trans_flux
                         deepSOM_p(ip,il+1,iv,iele)= deepSOM_p(ip,il+1,iv,iele)+(deepSOM_p(ip,il,iv,iele)/deepSOM_peat(ip,il,iv,iele))*trans_flux
                         deepSOM_peat(ip,il+1,iv,iele)= deepSOM_peat(ip,il+1,iv,iele)+ trans_flux                           
                      ENDIF
                      deepSOM_pt(ip,il,iv,iele)=deepSOM_peat(ip,il,iv,iele)/(zf_soil(il)-zf_soil(il-1))
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDDO 
    ENDIF

    IF (perma_peat) THEN
       DO ip = 1, kjpindex
          DO iv=1,nvm
             istm = pref_soil_veg(iv)
             DO iele = 1, nelements
                IF (veget_mask_2d(ip,iv)) THEN
                   peat_OLT(ip,iv,iele) = zero
                   IF (is_peat(iv) .AND. soiltile(ip,istm) .GT. zero) THEN
                      il=1
                      DO WHILE ((deepSOM_peat(ip,il,iv,iele) .GT. min_stomate) .AND. (il<ngrnd))   
                         SOMthick = (zf_soil(il)-zf_soil(il-1))*deepSOM_peat(ip,il,iv,iele)/SOMmax(il)
                         peat_OLT(ip,iv,iele)=zf_soil(il-1)+SOMthick
                         il=il+1
                      ENDDO
                      IF ((il==ngrnd) .AND. (deepSOM_peat(ip,il,iv,iele) .GT. min_stomate))THEN
                         SOMthick =(zf_soil(il)-zf_soil(il-1))*deepSOM_peat(ip,il,iv,iele)/SOMmax(il)
                         peat_OLT(ip,iv,iele)=zf_soil(il-1)+SOMthick
                         peat_OLT(ip,iv,iele)= MIN(zf_soil(il),peat_OLT(ip,iv, iele))
                      ENDIF
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDDO

    ENDIF
    !
    !! Check mass balance closure
    IF (err_act.GT.1) THEN

       !! 4.1.1 Calculate components of the mass balance
       pool_end(:,:,:) = zero 
       DO iele = 1,nelements
          DO igrnd = 1,ngrnd
             pool_end(:,:,iele) = pool_end(:,:,iele) + &
                  (deepSOM_a(:,igrnd,:,iele) + deepSOM_s(:,igrnd,:,iele) + deepSOM_p(:,igrnd,:,iele) + &
                  deltaSOM1_a(:,igrnd,:,iele) + deltaSOM1_s(:,igrnd,:,iele) + deltaSOM1_p(:,igrnd,:,iele)) * &
                  (zf_soil(igrnd)-zf_soil(igrnd-1)) * veget_max(:,:)
          END DO
       ENDDO
       
       !! 4.1.2 Calculate mass balance
       check_intern(:,:,:,:) = zero
       check_intern(:,:,iatm2land,:) = zero
       check_intern(:,:,iland2atm,:) = -un * zero
       check_intern(:,:,ilat2out,:) = -un *zero
       check_intern(:,:,ilat2in,:) =  zero
       check_intern(:,:,ipoolchange,:) = -un * (pool_end(:,:,:) - &
            pool_start(:,:,:))
       closure_intern = zero
       DO imbc = 1,nmbcomp
          DO iele = 1,nelements
             closure_intern(:,:,iele) = closure_intern(:,:,iele) + &
                  check_intern(:,:,imbc,iele)
          END DO
       ENDDO

       CALL check_mass_balance("permafrost_decomp", closure_intern, &
            kjpindex, pool_end, pool_start, veget_max, 'pft')

    END IF ! err_act.GT.1

  END SUBROUTINE permafrost_decomp


!!
!================================================================================================================================
!! SUBROUTINE   : calc_vert_int_som
!!
!>\BRIEF        This routine calculates carbon decomposition
!!
!! DESCRIPTION : 
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!! 
!! REFERENCE(S) : None
!!
!! FLOWCHART11    : None
!! \n
!_
!================================================================================================================================     

  SUBROUTINE calc_vert_int_som(kjpindex, deepSOM_a, deepSOM_s, deepSOM_p, som, som_surf, zf_soil)

  !! 0. Variable and parameter declaration

    !! 0.1  Input variables     

    INTEGER(i_std), INTENT(in)                     :: kjpindex     !! domain size
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)    :: deepSOM_a    !! active pool deepsom
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)    :: deepSOM_s    !! slow pool deepsom
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)    :: deepSOM_p    !! passive pool deepsom
    REAL(r_std), DIMENSION(0:), INTENT(in)         :: zf_soil      !! depths at full levels
   
    !! 0.2 Output variables

    REAL(r_std), DIMENSION(:,:,:,:),INTENT (out)   :: som          !! vertically-integrated som pool: active, slow, or passive, (gC/(m**2 of ground))
    REAL(r_std), DIMENSION(:,:,:,:),INTENT (out)   :: som_surf     !! vertically-integrated som pool to 1 meter: active, slow, or passive,(gC/(m**2 of ground))

    !! 0.3 Modified variables

    !! 0.4 Local variables
    INTEGER(i_std)                                            :: il,iele, ivm
    real(r_std), parameter                                    :: maxdepth=2.!! depth to which we intergrate the carbon for som_surf calculation                              

    som(:,:,:,:) = zero
    DO il = 1, ngrnd
     DO iele = 1, nelements
       WHERE ( veget_mask_2d(:,:) ) 
          som(:,iactive,:,iele) = som(:,iactive,:,iele) + deepSOM_a(:,il,:,iele)*(zf_soil(il)-zf_soil(il-1))
          som(:,islow,:,iele) = som(:,islow,:,iele) + deepSOM_s(:,il,:,iele)*(zf_soil(il)-zf_soil(il-1))
          som(:,ipassive,:,iele) = som(:,ipassive,:,iele) + deepSOM_p(:,il,:,iele)*(zf_soil(il)-zf_soil(il-1))
       END WHERE
     ENDDO
    ENDDO

    som_surf(:,:,:,:) = zero
    DO il = 1, ngrnd
     DO iele = 1, nelements
       IF (zf_soil(il-1) .lt. maxdepth ) THEN
          WHERE ( veget_mask_2d(:,:) ) 
             som_surf(:,iactive,:,iele) = som_surf(:,iactive,:,iele) + &
                  deepSOM_a(:,il,:,iele)*(min(maxdepth,zf_soil(il))-zf_soil(il-1))
             som_surf(:,islow,:,iele) = som_surf(:,islow,:,iele) + &
                  deepSOM_s(:,il,:,iele)*(min(maxdepth,zf_soil(il))-zf_soil(il-1))
             som_surf(:,ipassive,:,iele) = som_surf(:,ipassive,:,iele) + &
                  deepSOM_p(:,il,:,iele)*(min(maxdepth,zf_soil(il))-zf_soil(il-1))
          END WHERE
       ENDIF
     ENDDO
    ENDDO

  END SUBROUTINE calc_vert_int_som

END MODULE stomate_soil_carbon_discretization
