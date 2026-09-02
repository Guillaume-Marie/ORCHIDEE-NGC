! =================================================================================================================================/
! MODULE       : stomate_io
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       Module for read and write of restart files for all stomate modules.
!!
!!\n DESCRIPTION : This module contains the subroutines readstart and writerestart. All variables that will be read or written
!!                 are passed as argument to the subroutines. The subroutine readstart is called from stomate_initialize and 
!!                 writerestart is called from stomate_finalize.
!!                 Note: Not all variables saved in the start files are absolutely necessary. However, Sechiba's and Stomate's 
!!                 PFTs are not necessarily identical, and for that case this information needs to be saved.
!!
!!
!! RECENT CHANGE(S) : None
!!
!! REFERENCE(S)	: None
!!
!! SVN :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_stomate/stomate_io.f90 $
!! $Date: 2026-04-08 12:19:11 +0200 (mer. 08 avril 2026) $
!! $Revision: 9466 $
!! \n
!_ ================================================================================================================================
MODULE stomate_io
  USE stomate_data
  USE constantes
  USE constantes_mtc
  USE pft_parameters_var
  USE pft_parameters
  USE constantes_soil
  USE grid
  USE mod_orchidee_para
  USE ioipsl_para
  USE structures
  USE interpol_help
  USE time, ONLY : one_year
  USE xios_orchidee
  USE sapiens_forestry,   ONLY : sapiens_forestry_set_fm, sapiens_forestry_set_species_change, sapiens_forestry_set_desired_fm

  IMPLICIT NONE
  
  PRIVATE
  PUBLIC readrestart, writerestart

  

  REAL(r_std),ALLOCATABLE,DIMENSION(:),SAVE       :: trefe         !! reference temperature (K)
!$OMP THREADPRIVATE(trefe)
CONTAINS


!! ================================================================================================================================
!! SUBROUTINE   : readrestart
!!
!>\BRIEF        Read all variables for stomate from restart file. 
!!
!! DESCRIPTION  : Read all variables for stomate from restart file. 
!!                Initialize the variables if they were not found in the restart file or if there was no restart file.
!!                
!! \n
!_ ================================================================================================================================

  SUBROUTINE readrestart &
       & (npts, index, lalo, temp_air, dt_days, date_loc, &
       &  adapted, regenerate, vegstress_day, gdd_init_date, litterhum_daily, &
       &  t2m_daily, t2m_min_daily, tsurf_daily, t2m_max_daily, tsoil_daily, &
       &  precip_daily, vpd_daily_mean, vpd_daily_max,vpd_mean_week, vpd_max_week, &
       &  gpp_daily, npp_daily, turnover_daily, turnover_resid, &
       &  vegstress_month, vegstress_week, vegstress_season, &
       &  t2m_longterm, tau_longterm, t2m_month, t2m_week, &
       &  tsoil_month,precip_month, fireindex, firelitter, &
       &  maxvegstress_lastyear, maxvegstress_thisyear, &
       &  minvegstress_lastyear, minvegstress_thisyear, &
       &  maxgppweek_lastyear, maxgppweek_thisyear, &
       &  gdd0_lastyear, gdd0_thisyear, &
       &  precip_lastyear, precip_thisyear, &
       &  gdd_m5_dormance,  gdd_from_growthinit, gdd_midwinter, &
       &  ncd_dormance, ngd_minus5, &
       &  PFTpresent, npp_longterm, croot_longterm, n_reserve_longterm, lm_lastyearmax, &
       &  lm_thisyearmax, maxfpc_lastyear, maxfpc_thisyear, &
       &  turnover_longterm, gpp_week, gpp_year, gpp_decade, &
       &  resp_maint_part, resp_maint_week, &
       &  leaf_age, leaf_frac, leaf_age_crit, plant_status, when_growthinit, age, &
       &  resp_hetero, resp_maint, resp_growth, co2_fire, &
       &  veget_lastlight, everywhere, need_adjacent, RIP_time, &
       &  time_hum_min, hum_min_dormance, &
       &  litter, dead_leaves, &
       &  som, lignin_struc, lignin_wood, lignin_snag, turnover_time, &
       &  fco2_lu, fco2_wh, fco2_ha, &
       &  prod_s, prod_m, prod_l, flux_s, flux_m, flux_l, &
       &  fDeforestToProduct, fLulccResidue, fHarvestToProduct, &
       &  bm_to_litter, bm_to_litter_resid, tree_bm_to_litter, &
       &  tree_bm_to_litter_resid, carb_mass_total, &
       &  Tseason, Tseason_length, Tseason_tmp, & 
       &  Tmin_spring_time, &
       &  global_years, ok_equilibrium, nbp_accu_flux, &
       &  nbp_pool_start, &
       &  MatrixV, VectorU, previous_stock, current_stock, &
       &  assim_param, CN_som_litter_longterm, &
       &  tau_CN_longterm, KF, k_latosa_adapt, &
       &  rue_longterm, cn_leaf_min_season, nstress_season, &
       &  soil_n_min, p_O2, bact, &
       &  forest_managed, &
       &  species_change_map, fm_change_map, lpft_replant, lai_per_level, &
       &  laieff_fit, wstress_season, wstress_month, &
       &  age_stand, age_stand_bm, rotation_n, last_cut, mai, pai, &
       &  previous_wood_volume, mai_count, coppice_dens, &
       &  light_tran_to_floor_season, daylight_count, veget_max, gap_area_save, &
       &  deepSOM_a, deepSOM_s, deepSOM_p, O2_soil, CH4_soil, O2_snow, CH4_snow, &
       &  heat_Zimov, altmax, altmax_ind, altmax_lastyear, altmax_ind_lastyear, &
       &  depth_organic_soil, fixed_cryoturbation_depth, &
       &  cn_leaf_init_2D, sugar_load, harvest_cut, &
       &  harvest_pool_acc, harvest_area_acc, burried_litter, burried_fresh_ltr, &
       &  burried_fresh_som, burried_bact, &
       &  burried_min_nitro,burried_som, &
       &  burried_deepSOM_a, burried_deepSOM_s, burried_deepSOM_p, &
       &  wood_leftover_legacy, i_beetles_activity_legacy,season_drought_legacy, &
       &  P_beetles_attacked_legacy, beetle_diapause, sumTeff, woody_litter_to_use, &
       &  i_beetles_generation, B_beetles_kill_legacy, max_wind_speed_storm, &
       &  max_wind_ratio_storm, wind_ratio_max_save, wind_ratio_sum_save, & ! JJ 2026: added companions
       &  wind_speed_max_save, wind_mean_save, wind_sum, is_storm, count_storm, &
       &  biomass_init_drought, kill_vessels, &
       &  ni_acc, litterfuel, &
       &  vessel_loss_previous, grow_season_len, doy_start_gs, doy_end_gs, &
       &  mean_start_gs, total_ba_init, & 
       &  carbon_acro, carbon_cato, height_acro, deepSOM_peat)


    IMPLICIT NONE

    ! 0 declarations
    !-
    ! 0.1 input
    !-
    INTEGER(i_std),INTENT(in)                              :: npts                     !! Domain size
    INTEGER(i_std),DIMENSION(:),INTENT(in)                 :: index                    !! Indices of the points on the map
    REAL(r_std),DIMENSION(:,:),INTENT(in)                  :: lalo                     !! Geogr. coordinates (latitude,longitude) (degrees)
    REAL(r_std),DIMENSION(:),INTENT(in)                    :: temp_air                 !! Air temperature from forcing file or coupled model (K)
    REAL(r_std),DIMENSION(:,:), INTENT(in)                 :: cn_leaf_init_2D          !! initial leaf C/N ratio 
    REAL(r_std),DIMENSION(:,:),INTENT(in)                  :: veget_max                !! Maximum fraction of vegetation type including 
                                                                                       !! non-biological fraction (unitless)
 
    !-
    ! 0.2 output
    !-
    REAL(r_std),INTENT(out)                                :: dt_days                  !! time step of STOMATE in days
    INTEGER(i_std),INTENT(out)                             :: date_loc                 !! date_loc (d)
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: adapted                  !! Winter too cold? between 0 and 1
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: regenerate               !! Winter sufficiently cold? between 0 and 1
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: vegstress_day         !! daily moisture availability
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: gdd_init_date            !! date for beginning of gdd count
    REAL(r_std),DIMENSION(:),INTENT(out)                   :: litterhum_daily          !! daily litter humidity
    REAL(r_std),DIMENSION(:),INTENT(out)                   :: t2m_daily                !! daily 2 meter temperatures (K)
    REAL(r_std),DIMENSION(:),INTENT(out)                   :: t2m_min_daily            !! daily minimum 2 meter temperatures (K)
    REAL(r_std),DIMENSION(:),INTENT(out)                   :: t2m_max_daily            !! daily maximum 2 meter temperatures (K)
    REAL(r_std),DIMENSION(:),INTENT(out)                   :: tsurf_daily              !! daily surface temperatures (K)
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: tsoil_daily              !! daily soil temperatures (K)
    REAL(r_std),DIMENSION(:),INTENT(out)                   :: precip_daily             !! daily precipitations (mm/day) (for phenology)
    REAL(r_std),DIMENSION(:),INTENT(out)                   :: vpd_daily_mean           !! daily mean vpd (kPa/day) (for spitfire)
    REAL(r_std),DIMENSION(:),INTENT(out)                   :: vpd_daily_max            !! daily max vpd (kPa/day) (for spitfire)
    REAL(r_std),DIMENSION(:),INTENT(out)                   :: vpd_mean_week            !! weekly mean vpd (kPa/day) (for spitfire)
    REAL(r_std),DIMENSION(:),INTENT(out)                   :: vpd_max_week             !! weekly max vpd (kPa/day) (for spitfire)
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: gpp_daily                !! daily gross primary productivity (gC/m**2/day)
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: npp_daily                !! daily net primary productivity (gC/m**2/day)
    REAL(r_std),DIMENSION(:,:,:,:),INTENT(out)             :: turnover_daily           !! daily turnover rates (gC/m**2/day)
    REAL(r_std),DIMENSION(:,:,:,:),INTENT(out)             :: turnover_resid           !! The turnover left from turnover_daily at any given time step  
                                                                                       !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex 
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: vegstress_month       !! "monthly" moisture availability
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: vegstress_week        !! "weekly" moisture availability
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: vegstress_season      !! mean growing season moisture availability (used for allocation response)
    REAL(r_std),DIMENSION(:),INTENT(out)                   :: t2m_longterm             !! "long term" 2 meter temperatures (K)
    REAL(r_std), INTENT(out)                               :: tau_longterm             !! "tau_longterm"
    REAL(r_std),DIMENSION(:),INTENT(out)                   :: t2m_month                !! "monthly" 2 meter temperatures (K)
    REAL(r_std),DIMENSION(:),INTENT(out)                   :: Tseason                  !! "seasonal" 2 meter temperatures (K) 
    REAL(r_std),DIMENSION(:),INTENT(out)                   :: Tseason_length           !! temporary variable to calculate Tseason
    REAL(r_std),DIMENSION(:),INTENT(out)                   :: Tseason_tmp              !! temporary variable to calculate Tseason
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: Tmin_spring_time         !!
    REAL(r_std),DIMENSION(:),INTENT(out)                   :: t2m_week                 !! "weekly" 2 meter temperatures (K)
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: tsoil_month              !! "monthly" soil temperatures (K)
    REAL(r_std),DIMENSION(:),INTENT(out)                   :: precip_month             !! "monthly" precipitation 
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: fireindex                !! Probability of fire
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: firelitter               !! Longer term total litter above the ground, gC/m**2 of ground
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: maxvegstress_lastyear !! last year's maximum moisture availability
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: maxvegstress_thisyear !! this year's maximum moisture availability
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: minvegstress_lastyear !! last year's minimum moisture availability
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: minvegstress_thisyear !! this year's minimum moisture availability
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: maxgppweek_lastyear      !! last year's maximum weekly GPP
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: maxgppweek_thisyear      !! this year's maximum weekly GPP
    REAL(r_std),DIMENSION(:),INTENT(out)                   :: gdd0_lastyear            !! last year's annual GDD0
    REAL(r_std),DIMENSION(:),INTENT(out)                   :: gdd0_thisyear            !! this year's annual GDD0
    REAL(r_std),DIMENSION(:),INTENT(out)                   :: precip_lastyear          !! last year's annual precipitation (mm/year)
    REAL(r_std),DIMENSION(:),INTENT(out)                   :: precip_thisyear          !! this year's annual precipitation (mm/year)
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: gdd_m5_dormance          !! growing degree days, threshold -5 deg C (for phenology)
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: gdd_from_growthinit      !! growing degree days, from begin of season
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: gdd_midwinter            !! growing degree days since midwinter (for phenology)
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: ncd_dormance             !! number of chilling days since leaves were lost (for phenology)
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: ngd_minus5               !! number of growing days, threshold -5 deg C (for phenology)
    LOGICAL,DIMENSION(:,:),INTENT(out)                     :: PFTpresent               !! PFT exists (equivalent to fpc_max > 0 for natural PFTs)
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: npp_longterm             !! "long term" net primary productivity (gC/m**2/year)
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: croot_longterm           !! "long term" root carbon mass (gC/m**2/year)
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: n_reserve_longterm       !! "long term" actual to potential N reserve pool (unitless)
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: lm_lastyearmax           !! last year's maximum leaf mass, for each PFT (gC/m**2)
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: lm_thisyearmax           !! this year's maximum leaf mass, for each PFT (gC/m**2)
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: maxfpc_lastyear          !! last year's maximum fpc for each natural PFT, on ground
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: maxfpc_thisyear          !! this year's maximum fpc for each PFT, on *total* ground (see stomate_season)   
    REAL(r_std),DIMENSION(:,:,:,:),INTENT(out)             :: turnover_longterm        !! "long term" turnover rate (gC/m**2/year)
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: gpp_week                 !! "weekly" GPP (gC/day/(m**2 covered)
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: gpp_year                 !! "annual" GPP (gC/day/(m**2 covered)
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: gpp_decade               !! "decadal" GPP (gC/day/(m**2 covered)
    REAL(r_std),DIMENSION(:,:,:),INTENT(out)               :: resp_maint_part          !! maintenance resp (gC/m**2)
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: resp_maint_week          !! "weekly" maintenance respiration (gC/day/(m**2 covered)
    REAL(r_std),DIMENSION(:,:,:),INTENT(out)               :: leaf_age                 !! leaf age (days)
    REAL(r_std),DIMENSION(:,:,:),INTENT(out)               :: leaf_frac                !! fraction of leaves in leaf age class
    REAL(r_std), DIMENSION(:,:),INTENT(out)                :: leaf_age_crit            !! critical leaf age (days)
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: plant_status             !! Growth and phenological status of the plant
                                                                                       !! The different stati are defined in constantes
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: when_growthinit          !! how many days ago was the beginning of the growing season
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: age                      !! mean age (years)
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: resp_hetero              !! heterotrophic respiration (gC/day/m**2)
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: resp_maint               !! maintenance respiration (gC/day/m**2)
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: resp_growth              !! growth respiration (gC/day/m**2)
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: co2_fire                 !! carbon emitted into the atmosphere by fire (living and dead biomass)
                                                                                       !! (in gC/m**2/time step)
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: veget_lastlight          !! vegetation fractions (on ground) after last light competition
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: everywhere               !! is the PFT everywhere in the grid box or very localized (after its introduction)
    LOGICAL,DIMENSION(:,:),INTENT(out)                     :: need_adjacent            !! in order for this PFT to be introduced, does it have to be present in an adjacent grid box?
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: RIP_time                 !! How much time ago was the PFT eliminated for the last time (y)
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: time_hum_min             !! time elapsed since strongest moisture availability (d)
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: hum_min_dormance         !! minimum moisture during dormance
    REAL(r_std),DIMENSION(:,:,:,:,:),INTENT(out)           :: litter                   !! fraction of litter above the ground belonging to different PFTs
                                                                                       !! separated for natural and agricultural PFTs.
    REAL(r_std),DIMENSION(:,:,:),INTENT(out)               :: dead_leaves              !! dead leaves on ground, per PFT, metabolic and structural in gC/(m**2 of ground)
    REAL(r_std),DIMENSION(:,:,:,:),INTENT(out)             :: som                      !! Soil Organic Matter pool: active, slow, or passive, (gC (or N)/m**2)
    REAL(r_std),DIMENSION(:,:,:),INTENT(out)               :: lignin_struc             !! ratio Lignine/Carbon in structural litter, above and below ground,(gC/m**2)
    REAL(r_std),DIMENSION(:,:,:),INTENT(out)               :: lignin_wood              !! ratio Lignine/Carbon in woody litter, above and below ground,(gC/m**2)
    REAL(r_std),DIMENSION(:,:,:),INTENT(out)               :: lignin_snag              !! ratio Lignine/Carbon in snag litter, above and below ground,(gC/m**2)
    REAL(r_std),DIMENSION(:,:,:),INTENT(out)               :: turnover_time            !!
    INTEGER(i_std), INTENT(out)                            :: global_years             !! for spinup matrix  
    LOGICAL, DIMENSION(:), INTENT(out)                     :: ok_equilibrium           !! for spinup matrix  
    REAL(r_std), DIMENSION(:,:), INTENT(out)               :: nbp_accu_flux            !! accumulated Net Biospheric Production over the whole simulationm (gC/N m-2)
    REAL(r_std), DIMENSION(:,:), INTENT(out)               :: nbp_pool_start           !! C an dN stocks at previous time step (gC/N m-2)
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(out)           :: MatrixV                  !!
    REAL(r_std), DIMENSION(:,:,:), INTENT(out)             :: VectorU                  !!
    REAL(r_std), DIMENSION(:,:,:), INTENT(out)             :: previous_stock           !!
    REAL(r_std), DIMENSION(:,:,:), INTENT(out)             :: current_stock            !!
    REAL(r_std), DIMENSION(:,:,:), INTENT(out)             :: CN_som_litter_longterm   !! Longterm CN ratio of litter and som pools (gC/gN)
    REAL(r_std), INTENT(out)                               :: tau_CN_longterm          !! Counter used for calculating the longterm CN ratio of SOM and litter pools (seconds)
    REAL(r_std), DIMENSION(:,:,:), INTENT(out)             :: assim_param              !!
    REAL(r_std), DIMENSION(:,:), INTENT(out)               :: KF                       !! Scaling factor to convert sapwood mass into leaf mass (m)
    REAL(r_std), DIMENSION(:,:), INTENT(out)               :: k_latosa_adapt           !! Leaf to sapwood area adapted for water stress. Adaptation takes place at the end of the year (m)
    REAL(r_std), DIMENSION(:,:), INTENT(out)               :: mai                      !! The mean annual increment @tex $(m**3 / m**2 / year)$ @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(out)               :: pai                      !! The period annual increment @tex $(m**3 / m**2 / year)$ @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(out)               :: previous_wood_volume     !! The volume of the tree trunks in a stand for the previous year. @tex $(m**3 / m**2 )$ @endtex
    INTEGER(i_std), DIMENSION(:,:),INTENT(out)             :: mai_count                !! The number of times we've calculated the volume increment for a stand
    REAL(r_std), DIMENSION(:,:),INTENT(out)                :: coppice_dens             !! The density of a coppice at the first cutting. @tex $( 1 / m**2 )$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(out)               :: rue_longterm             !! longterm radiation use efficiency
    INTEGER(i_std), DIMENSION(:,:), INTENT(out)            :: age_stand                !! Age of stand (years)
    REAL(r_std), DIMENSION(:,:), INTENT(out)               :: age_stand_bm             !! Biomass-weighted conserved mean stand age (years) - STAND_AGE
    INTEGER(i_std), DIMENSION(:,:), INTENT(out)            :: rotation_n               !! Rotation number (number of rotation since pft is managed)
    INTEGER(i_std), DIMENSION(:,:), INTENT(out)            :: last_cut                 !! Years since last thinning (years)
    REAL(r_std), DIMENSION(:,:), INTENT(out)               :: cn_leaf_min_season       !! Seasonal min CN ratio of leaves 
    REAL(r_std), DIMENSION(:,:), INTENT(out)               :: nstress_season           !! N-related seasonal stress (used for allocation) 
    REAL(r_std), DIMENSION(:,:,:), INTENT(out)             :: soil_n_min               !! mineral nitrogen in the soil (gN/m**2) (first index=npts, second index=nvm, third index=nnspec) 
    REAL(r_std), DIMENSION(:,:), INTENT(out)               :: p_O2                     !! partial pressure of oxigen in the soil (hPa)(first index=npts, second index=nvm)
    REAL(r_std), DIMENSION(:,:), INTENT(out)               :: bact                     !! denitrifier biomass (gC/m**2) (first index=npts, second index=nvm)
    INTEGER(i_std), DIMENSION(:,:), INTENT(out)            :: forest_managed           !! forest management flag
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(out)         :: prod_s                   !!
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(out)         :: prod_m                   !!
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(out)         :: prod_l                   !!
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(out)         :: flux_s                   !!
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(out)         :: flux_m                   !!
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(out)         :: flux_l                   !!
    INTEGER(i_std), DIMENSION(:,:), INTENT(out)            :: species_change_map       !! A map which gives the PFT number that each PFT will be replanted as in case of a clearcut.
                                                                                       !! (1-nvm,unitless)
    INTEGER(i_std), DIMENSION(:,:), INTENT(out)            :: fm_change_map            !! A map which gives the desired FM strategy when the PFT will be replanted after a clearcut.
                                                                                       !! (1-nvm,unitless)
    LOGICAL, DIMENSION(:,:), INTENT(out)                   :: lpft_replant             !! Indicates if this PFT has either died this year or been clearcut/coppiced.  If it has, it is not
                                                                                       !! replanted until the end of the year.
    REAL(r_std), DIMENSION(:,:,:),INTENT(out)              :: lai_per_level            !! The amount of LAI in each physical canopy level. @tex $( m**2 / m**2 )$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:),INTENT(out)            :: deepSOM_a                !!
    REAL(r_std), DIMENSION(:,:,:,:),INTENT(out)            :: deepSOM_s                !!
    REAL(r_std), DIMENSION(:,:,:,:),INTENT(out)            :: deepSOM_p                !!
    REAL(r_std), DIMENSION(:,:,:),INTENT(out)              :: O2_soil                  !!
    REAL(r_std), DIMENSION(:,:,:),INTENT(out)              :: CH4_soil                 !!
    REAL(r_std), DIMENSION(:,:,:),INTENT(out)              :: O2_snow                  !!
    REAL(r_std), DIMENSION(:,:,:),INTENT(out)              :: CH4_snow                 !!
    REAL(r_std), DIMENSION(:,:,:), INTENT(out)             :: heat_Zimov               !! heating associated with decomposition [W/m**3 soil]
    REAL(r_std), DIMENSION(:,:), INTENT(out)               :: altmax                   !! Active layer thickness (m) 
    INTEGER(i_std), DIMENSION(npts,nvm), INTENT(out)       :: altmax_ind
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)          :: altmax_lastyear 
    INTEGER(i_std), DIMENSION(npts,nvm), INTENT(out)       :: altmax_ind_lastyear
    REAL(r_std), DIMENSION(:),INTENT(out)                  :: depth_organic_soil       !! Depth at which there is still organic matter (m)
    REAL(r_std), DIMENSION(:,:), INTENT(out)               :: fixed_cryoturbation_depth!! Depth to hold cryoturbation to for fixed runs  
    REAL(r_std), DIMENSION(:,:), INTENT(out)               :: sugar_load               !! Relative sugar loading of the labile pool (unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(out)               :: harvest_cut              !! Type of cutting that was used for the harvest (unitless) 
    TYPE(laieff_type),DIMENSION (:,:,:),INTENT(out)        :: laieff_fit               !! Fitted parameters for the effective LAI
    REAL(r_std), DIMENSION(:,:), INTENT(out)               :: wstress_season           !! Water stress factor, based on hum_rel_daily (unitless, 0-1)
    REAL(r_std), DIMENSION(:,:), INTENT(out)               :: wstress_month            !! Water stress factor, based on hum_rel_daily (unitless, 0-1)
    REAL(r_std), DIMENSION(:,:), INTENT(out)               :: fDeforestToProduct       !!
    REAL(r_std), DIMENSION(:,:), INTENT(out)               :: fLulccResidue            !!
    REAL(r_std), DIMENSION(:,:), INTENT(out)               :: fHarvestToProduct        !!
    REAL(r_std),DIMENSION(:,:,:,:),INTENT(out)             :: bm_to_litter             !! Background (not senescence-driven) mortality of biomass
                                                                                       !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex
    REAL(r_std),DIMENSION(:,:,:,:),INTENT(out)             :: bm_to_litter_resid       !! Left over bm_to_litter at any specific time step
                                                                                       !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex
    REAL(r_std),DIMENSION(:,:,:,:),INTENT(out)             :: tree_bm_to_litter        !! Conversion of biomass to litter 
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std),DIMENSION(:,:,:,:),INTENT(out)             :: tree_bm_to_litter_resid  !! Left over bm_to_litter_resid. Written here, used in stomate.f90 
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std),DIMENSION(:),INTENT(out)                   :: carb_mass_total          !!
    REAL(r_std), DIMENSION(:,:), INTENT(out)               :: light_tran_to_floor_season !! Mean seasonal fraction of light transmitted to the forest floor (unitless, 0-1)
    REAL(r_std), DIMENSION(:,:), INTENT(out)               :: daylight_count           !! Time steps dt_radia during daylight and when there is growth (gpp>0)
    REAL(r_std), DIMENSION(:,:,:), INTENT(out)             :: gap_area_save            !! Total gap area created by more than 30% basal area loss
                                                                                       !! in the last 5 years (m^{2})      
    REAL(r_std), DIMENSION(:,:), INTENT(out)               :: total_ba_init            !! Total basal area saved at the first day of the year (m^{2}/m^{2})
    REAL(r_std), DIMENSION(:,:,:), INTENT(out)             :: season_drought_legacy    !! mean growing season moisture availability
    REAL(r_std), DIMENSION(:,:,:), INTENT(out)             :: wood_leftover_legacy     !!
    REAL(r_std), DIMENSION(:,:), INTENT(out)               :: P_beetles_attacked_legacy!!
    INTEGER(i_std), DIMENSION(:,:), INTENT(out)            :: beetle_diapause          !!
    REAL(r_std), DIMENSION(:,:), INTENT(out)               :: sumTeff                  !! sum of temp of beetle phenology
    REAL(r_std), DIMENSION(:,:), INTENT(out)               :: woody_litter_to_use      !! Yearly sum of woody litter pool to use for wood_leftover in the pest module
                                                                                       !! (gC m^{-2})
    REAL(r_std), DIMENSION(:,:,:),INTENT(out)              :: i_beetles_activity_legacy        !! biomass of tree from the same species that was infected during the previous timestep 
    REAL(r_std), DIMENSION(:,:,:), INTENT(out)             :: B_beetles_kill_legacy    !!
    REAL(r_std), DIMENSION(:,:,:), INTENT(out)             :: i_beetles_generation  !! number of generation that BB can achieved in one year
    LOGICAL, DIMENSION(:),INTENT(out)                      :: is_storm
    INTEGER(i_std), DIMENSION(:),INTENT(out)               :: count_storm
    REAL(r_std), DIMENSION(:), INTENT(out)                 :: max_wind_speed_storm     !! Daily maximum wind speed at 2 meter (ms-1) during a storm event
    REAL(r_std), DIMENSION(:), INTENT(out)                 :: max_wind_ratio_storm     !! JJ 2026: Daily maximum wind ratio during a storm event
    REAL(r_std), DIMENSION(:,:), INTENT(out)               :: wind_ratio_sum_save      !! Sum of wind speed ratio over longterm wind speed when this ratio exceeded certain threshold,
                                                                                       !! stored for wind_days
    REAL(r_std), DIMENSION(:,:), INTENT(out)               :: wind_ratio_max_save      !! JJ 2026: Daily maximum wind ratio, stored for wind_days (default 3)
    REAL(r_std), DIMENSION(:,:), INTENT(out)               :: wind_speed_max_save      !! Daily maximum wind speed, stored for wind_days (default 3) (m s-1)
    REAL(r_std), DIMENSION(:,:), INTENT(out)               :: wind_mean_save           !! Annual average wind speed, stored for legacy years (m s-1)
    REAL(r_std), DIMENSION(:), INTENT(out)                 :: wind_sum                 !! Sum of wind speed for all time step (m s-1)
    REAL(r_std),DIMENSION(:,:,:,:,:),INTENT(out)           :: harvest_pool_acc         !! Records the quantity of wood harvested and thinned due to forest management and LCC.
    REAL(r_std), DIMENSION(:,:,:), INTENT(out)             :: harvest_area_acc         !! Harvested area (m^{2})
    REAL(r_std),DIMENSION(:),INTENT(out)                   :: fco2_lu                  !!
    REAL(r_std),DIMENSION(:),INTENT(out)                   :: fco2_wh                  !!
    REAL(r_std),DIMENSION(:),INTENT(out)                   :: fco2_ha                  !!
    REAL(r_std),DIMENSION(:,:,:,:),INTENT(out)             :: burried_litter           !! Litter burried under non-biological land uses (gC orNm-2)
    REAL(r_std),DIMENSION(:,:,:),INTENT(out)               :: burried_fresh_ltr        !! Fresh litter burried under non-biological land uses (gC orN m-2)
    REAL(r_std),DIMENSION(:,:,:),INTENT(out)               :: burried_fresh_som        !! Fresh som burried under non-biological land uses (gC or Nm-2)
    REAL(r_std),DIMENSION(:),INTENT(out)                   :: burried_bact             !! Bacteria burried under non-biological land uses (gC m-2)
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: burried_min_nitro        !! Mineral nitrogen burried under non-biological land uses(gC or N m-2)
    REAL(r_std),DIMENSION(:,:,:),INTENT(out)               :: burried_som              !! Som burried under non-biological land uses (gC or N m-2)
    REAL(r_std),DIMENSION(:,:,:),INTENT(out)               :: burried_deepSOM_a        !! Som burried under non-biological land uses (gC or N m-2)
    REAL(r_std),DIMENSION(:,:,:),INTENT(out)               :: burried_deepSOM_s        !! Som burried under non-biological land uses (gC or N m-2)
    REAL(r_std),DIMENSION(:,:,:),INTENT(out)               :: burried_deepSOM_p        !! Som burried under non-biological land uses (gC or N m-2)
    REAL(r_std),DIMENSION(:,:,:,:,:),INTENT(out)           :: biomass_init_drought     !! Biomass of heartwood or sapwood before onset of drought. 
                                                                                       !! Used to compute turnover on same reference biomass in 
                                                                                       !! stomate_turnover.f90. Should remain the same along one 
                                                                                       !! entire drought episode and be updated inbetween 
                                                                                       !! droughts (gCor N tree-1).
    LOGICAL,DIMENSION(:,:,:),INTENT(out)                    :: kill_vessels             !! Flag to kill vessels at the end of the day when there is embolism.
    REAL(r_std),DIMENSION(:),INTENT(out)                    :: ni_acc                   !! Nesterov index (square of degree Celcius)
    REAL(r_std), DIMENSION(:,:,:,:,:),INTENT(out)           :: litterfuel               !! Dead litter fuel above ground. (gC m^{-2})
    REAL(r_std),DIMENSION(:,:,:),INTENT(out)               :: vessel_loss_previous     !! vessel loss at the previous time step
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: grow_season_len          !! growing season length in days for deciduous PFTs. 
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: doy_start_gs             !! growing season starting day of year (DOY) for deciduous PFTs.
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: doy_end_gs               !! growing season end day of year (DOY) for deciduous PFTs.
    REAL(r_std),DIMENSION(:,:),INTENT(out)                 :: mean_start_gs            !! mean growing season starting day for deciduous PFTs.

    ! peatland 
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out)            :: carbon_acro              !! carbon of acrotelm (gC/m**2), if ok_peat_NoDiscretisation activated. 
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out)            :: carbon_cato              !! carbon of catotelm (gC/m**2)if ok_peat_NoDiscretisation activated. 
    REAL(r_std),DIMENSION(npts),INTENT(out)                :: height_acro              !! height of acrotelm if ok_peat_NoDiscretisation activated. 
    REAL(r_std), DIMENSION(npts,ngrnd,nvm, nelements),INTENT(out) ::deepSOM_peat       !! soil carbon in peatland
                                                      

  !! 0.4 Local variables
    
    REAL(r_std)                                                         :: date_real
    REAL(r_std),DIMENSION(npts,nvm)                                     :: PFTpresent_real         !! PFT exists (equivalent to fpc_max > 0 for natural PFTs), real
    REAL(r_std),DIMENSION(npts)                                         :: is_storm_real
    REAL(r_std),DIMENSION(npts)                                         :: count_storm_real

    REAL(r_std),DIMENSION(npts,nvm)                                     :: need_adjacent_real      !! in order for this PFT to be introduced,
                                                                                                   !! does it have to be present in an adjacent grid box? - real
    CHARACTER(LEN=80)                                                   :: var_name                !! To store variables names for I/O
    CHARACTER(LEN=10)                                                   :: part_st                 !! string suffix indicating an index
    CHARACTER(LEN=10)                                                   :: circ_str                !! string suffix indicating an index            
    REAL(r_std),DIMENSION(1)                                            :: xtmp                    !! temporary storage
    INTEGER(i_std)                                                      :: j,k,l,m,ivma            !! index
    CHARACTER(LEN=2),DIMENSION(nelements)                               :: element_str             !! string suffix indicating element
    CHARACTER(LEN=6), DIMENSION(nbpools)                                :: pools_str
    REAL(r_std), DIMENSION(npts)                                        :: ok_equilibrium_real    
    INTEGER                                                             :: n,ilev,ipts, igrn       !! Indices
    INTEGER                                                             :: ivm,icarb,iele,ilit     !! Indices
    INTEGER                                                             :: alloc_err               !! Error indice from allocate
    REAL(r_std), DIMENSION(npts,nvm,ncirc)                              :: temp_real2              !! temporary real to allow restget 
                                                                                                   !! to work on multi-dimensional integers
    REAL(r_std), DIMENSION(npts,nvm)                                    :: temp_real               !! temporary real to allow restget 
                                                                                                   !! to work on multi-dimensional integers
    REAL(r_std), DIMENSION(npts,nvm)                                    :: r_replant               !! Getting logical values from the restart
                                                                                                   !! is not possible, so this is a temporary
                                                                                                   !! array where 1.0 is TRUE.
    REAL(r_std),DIMENSION(npts,nvm,nlevels_tot,nparams_laieff)          :: temp_array              !! To store structure values for I/O
    CHARACTER(LEN=10)                                                   :: part_str2               !! string suffix indicating an index
    CHARACTER(LEN=2), DIMENSION(legacy_years_wind+1)                    :: wyear_str               !! string suffix indicating wind year index
    CHARACTER(LEN=2)                                                    :: pyear_str               !! string suffix indicating year index for pest module 
    CHARACTER(LEN=10), DIMENSION(nlctypes)                              :: lctype_str              !! string suffix for the land cover type
    CHARACTER(LEN=10), DIMENSION(nlanduse)                              :: luse_str                !! string suffix for the land use type
    REAL(r_std), DIMENSION(0:ngrnd)                                     :: zf_soil
    CHARACTER(LEN=255)                                                  :: filename, name          !! Name of file for SOC initialization
    REAL(r_std), DIMENSION(nvm)                                         :: temp                    !! temporary variable enabling to use getin for variables
                                                                              !! with other dimensions than nvm

    ! Permafrost carbon processes
    LOGICAL :: read_input_deepC_a
    LOGICAL :: read_input_deepC_s
    LOGICAL :: read_input_deepC_p
    LOGICAL :: read_input_thawed_humidity
    LOGICAL :: read_input_depth_organic_soil

    ! AED_FEEDBACK / AED_SPINUP scratch (Marie 2026)
    REAL(r_std) :: class0_year_clock_loc
    REAL(r_std) :: n_class0_file                !! N_CLASS0 stored in the restart file (-)

!_ ================================================================================================================================

    IF (printlev >= 3) WRITE(numout,*) 'Entering readrestart'
    !-
    ! 1 string definitions
    !-
    DO l=1,nlctypes
       IF (l == iforest) THEN
          lctype_str(l) = '_forest'
       ELSEIF (l == igrass) THEN
          lctype_str(l) = '_grass'
       ELSEIF (l == icrop) THEN
          lctype_str(l) = '_crop'
       ELSE
          CALL ipslerr_p(3,'stomate_io readrestart','Define lctype_str(l)','','')
       ENDIF   
    END DO 
    !-
    DO l=1,nelements
       IF (l == icarbon) THEN
          element_str(l) = '_c'
       ELSEIF (l == initrogen) THEN
          element_str(l) = '_n'
       ELSE
          CALL ipslerr_p(3,'stomate_io readrestart','Define element_str','','')
       ENDIF
    ENDDO
    !-
    DO l=1,nlanduse
       IF (l == iharvest) THEN
          luse_str(l) = '_harvest'
       ELSEIF (l == ilcc) THEN
          luse_str(l) = '_lcc'
       ELSE
          CALL ipslerr_p(3,'stomate_io readrestart','Define luse_str','','')
       ENDIF
    ENDDO
    !-
    pools_str(1:nbpools) =(/'str_ab ','str_be ','met_ab ','met_be ','wood_ab','wood_be',& 
         & 'snag_ab', 'snag_be', 'actif  ','slow   ','passif ','surface'/)

    ! Vertical soil layers
    zf_soil(1:ngrnd) = zlt(:)
    zf_soil(0) = 0.

    !-
    ! 2 run control
    !-
    ! 2.2 time step of STOMATE in days
    !-    If the variable is not in the restart file, then un will be used as default value
    CALL restget_p(rest_id_stomate, 'dt_days', itime, .TRUE., un, dt_days)
    !-
    ! 2.3 date
    !-    If the variable is not in the restart file, then zero will be used as default value
    CALL restget_p (rest_id_stomate, 'date', itime, .TRUE., zero, date_real)
    date_loc = NINT(date_real)
    !-
    ! 3 daily meteorological variables
    !-
    vegstress_day(:,:) = val_exp
    var_name = 'vegstress_day'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., vegstress_day, 'gather', nbp_glo, index_g)
    IF (ALL(vegstress_day(:,:) == val_exp)) vegstress_day(:,:) = zero
    !-
    gdd_init_date(:,:) = val_exp
    var_name = 'gdd_init_date'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 2, 1, itime, &
         &              .TRUE., gdd_init_date, 'gather', nbp_glo, index_g)
    ! Keep val_exp as initial value for gdd_init_date(:,2)
    IF (ALL(gdd_init_date(:,1) == val_exp)) gdd_init_date(:,1) = 365.

    !-
    litterhum_daily(:) = val_exp
    var_name = 'litterhum_daily'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &              .TRUE., litterhum_daily, 'gather', nbp_glo, index_g)
    IF (ALL(litterhum_daily(:) == val_exp)) litterhum_daily(:) = zero
    !-
    t2m_daily(:) = val_exp
    var_name = 't2m_daily'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &                .TRUE., t2m_daily, 'gather', nbp_glo, index_g)
    IF (ALL(t2m_daily(:) == val_exp)) t2m_daily(:) = zero
    !-
    t2m_min_daily(:) = val_exp
    var_name = 't2m_min_daily'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &                .TRUE., t2m_min_daily, 'gather', nbp_glo, index_g)
    IF (ALL(t2m_min_daily(:) == val_exp)) t2m_min_daily(:) = large_value
    !-
    t2m_max_daily(:) = val_exp
    var_name = 't2m_max_daily'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &                .TRUE., t2m_max_daily, 'gather', nbp_glo, index_g)
    IF (ALL(t2m_max_daily(:) == val_exp)) THEN
        t2m_max_daily(:) = zero
    ENDIF

    !-
    tsurf_daily(:) = val_exp
    var_name = 'tsurf_daily'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &                .TRUE., tsurf_daily, 'gather', nbp_glo, index_g)
    ! The initial value is set to the current temperature at 2m
    IF (ALL(tsurf_daily(:) == val_exp)) tsurf_daily(:) = temp_air(:)
    !-
    tsoil_daily(:,:) = val_exp
    var_name = 'tsoil_daily'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo,   nslm, 1, itime, &
         &                .TRUE., tsoil_daily, 'gather', nbp_glo, index_g)
    IF (ALL(tsoil_daily(:,:) == val_exp)) tsoil_daily(:,:) = zero
    !-
    precip_daily(:) = val_exp
    var_name = 'precip_daily'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &                .TRUE., precip_daily, 'gather', nbp_glo, index_g)
    IF (ALL(precip_daily(:) == val_exp)) precip_daily(:) = zero

    vpd_daily_mean(:) = val_exp
    var_name = 'vpd_daily_mean'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &                .TRUE., vpd_daily_mean, 'gather', nbp_glo, index_g)
    IF (ALL(vpd_daily_mean(:) == val_exp)) vpd_daily_mean(:) = zero

    vpd_daily_max(:) = val_exp
    var_name = 'vpd_daily_max'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &                .TRUE., vpd_daily_max, 'gather', nbp_glo, index_g)
    IF (ALL(vpd_daily_max(:) == val_exp)) vpd_daily_max(:) = zero

    vpd_mean_week(:) = val_exp
    var_name = 'vpd_mean_week'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &                .TRUE., vpd_mean_week, 'gather', nbp_glo, index_g)
    IF (ALL(vpd_mean_week(:) == val_exp)) vpd_mean_week(:) = zero

    vpd_max_week(:) = val_exp
    var_name = 'vpd_max_week'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &                .TRUE., vpd_max_week, 'gather', nbp_glo, index_g)
    IF (ALL(vpd_max_week(:) == val_exp)) vpd_max_week(:) = zero

    !-
    ! 4 productivities
    !-
    gpp_daily(:,:) = val_exp
    var_name = 'gpp_daily'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., gpp_daily, 'gather', nbp_glo, index_g)
    IF (ALL(gpp_daily(:,:) == val_exp)) gpp_daily(:,:) = zero
    !-
    npp_daily(:,:) = val_exp
    var_name = 'npp_daily'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., npp_daily, 'gather', nbp_glo, index_g)
    IF (ALL(npp_daily(:,:) == val_exp)) npp_daily(:,:) = zero
    !-
    turnover_daily(:,:,:,:) = val_exp
    CALL restget_p (rest_id_stomate, 'turnover_daily', nbp_glo, nvm, nparts, nelements, itime, & 
        &                .TRUE., turnover_daily, 'gather', nbp_glo, index_g) 
    IF (ALL(turnover_daily  == val_exp)) turnover_daily(:,:,:,:) = zero
    !-
    CALL restget_p (rest_id_stomate, 'turnover_resid', nbp_glo, nvm, nparts, nelements, itime, & 
        &                .TRUE., turnover_resid, 'gather', nbp_glo, index_g) 
    IF (ALL(turnover_resid  == val_exp)) turnover_resid(:,:,:,:) = zero

    !-
    ! 5 monthly meteorological variables
    !-
    ! The following variables vegstress_month and vegstress_week
    ! are not initialized here if they were not found in the restart file. Initalization will be done in 
    ! stomate_season with daily variables calculated in stomate.
    vegstress_month(:,:) = val_exp
    var_name = 'vegstress_month'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., vegstress_month, 'gather', nbp_glo, index_g)
    
    vegstress_week(:,:) = val_exp
    var_name = 'vegstress_week'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., vegstress_week, 'gather', nbp_glo, index_g)

    ! vegstress_season is intialized to 1 if not found in restart file in oposite of the variable above
    vegstress_season(:,:) = val_exp
    var_name = 'vegstress_season'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., vegstress_season, 'gather', nbp_glo, index_g)
    IF (ALL(vegstress_season(:,:) == val_exp)) THEN
       vegstress_season(:,:) = un
       WHERE (veget_max(:,:).LT.min_stomate)
          vegstress_season(:,:) = zero
       END WHERE
    END IF



    !
    ! Longterm temperature at 2m
    !
    var_name = 't2m_longterm'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &              .TRUE., t2m_longterm, 'gather', nbp_glo, index_g)

    IF (ALL(t2m_longterm(:) == val_exp)) THEN
       ! t2m_longterm is not in restart file
       ! The initial value for the reference temperature is set to the current temperature
       t2m_longterm(:)=temp_air(:)
       ! Set the counter to 2 time steps
       tau_longterm=2
    ELSE
       ! t2m_longterm was in the restart file
       ! Now read tau_longterm
       ! tau_longterm is a scalar, therefor only master process read this value
       CALL restget_p (rest_id_stomate, 'tau_longterm', itime, .TRUE., val_exp, tau_longterm)
       IF (tau_longterm == val_exp) THEN
             ! tau_longterm is not found in restart file. 
             ! This is not normal as t2m_longterm was in restart file. Write a warning and initialize it to tau_longterm_max
          CALL ipslerr_p(2, 'stomate_io readrestart','tau_longterm was not in restart file',&
               'But t2m_longterm was in restart file','')
          tau_longterm = tau_longterm_max
       END IF

    END IF
    !-
    t2m_month(:) = val_exp
    var_name = 't2m_month'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &              .TRUE., t2m_month, 'gather', nbp_glo, index_g)
    IF (ALL(t2m_month(:) == val_exp)) t2m_month(:) = temp_air(:)
    
    CALL restget_p (rest_id_stomate, 'Tseason', nbp_glo, 1     , 1, itime, &
         .TRUE., Tseason, 'gather', nbp_glo, index_g)
    IF (ALL(Tseason(:) == val_exp)) Tseason(:) = temp_air(:)
    
    CALL restget_p (rest_id_stomate,'Tseason_length', nbp_glo, 1     , 1, itime, &
         .TRUE., Tseason_length, 'gather', nbp_glo, index_g)
    IF (ALL(Tseason_length(:) == val_exp)) Tseason_length(:) = zero
    
    CALL restget_p (rest_id_stomate, 'Tseason_tmp', nbp_glo, 1     , 1, itime, &
         .TRUE., Tseason_tmp, 'gather', nbp_glo, index_g)
    IF (ALL(Tseason_tmp(:) == val_exp)) Tseason_tmp(:) = zero

    CALL restget_p (rest_id_stomate, 'Tmin_spring_time', nbp_glo, nvm, 1, itime, &
         .TRUE., Tmin_spring_time, 'gather', nbp_glo, index_g)
    IF (ALL(Tmin_spring_time(:,:) == val_exp)) Tmin_spring_time(:,:) = zero

    t2m_week(:) = val_exp
    var_name = 't2m_week'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &              .TRUE., t2m_week, 'gather', nbp_glo, index_g)
    ! The initial value is set to the current temperature
    IF (ALL(t2m_week(:) == val_exp)) t2m_week(:) = temp_air(:)
    
    tsoil_month(:,:) = val_exp
    var_name = 'tsoil_month'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo,   nslm, 1, itime, &
         &              .TRUE., tsoil_month, 'gather', nbp_glo, index_g)

    ! The initial value is set to the current temperature
    IF (ALL(tsoil_month(:,:) == val_exp)) THEN
       DO l=1,nslm
          tsoil_month(:,l) = temp_air(:)
       ENDDO
    ENDIF
    !-
    precip_month(:) = val_exp
    var_name = 'precip_month'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo,   1, 1, itime, &
         &              .TRUE., precip_month, 'gather', nbp_glo, index_g)
    IF (ALL(precip_month(:) == val_exp))  precip_month(:) = zero
    !-
    wstress_season(:,:) = val_exp
    var_name = 'wstress_season'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                .TRUE., wstress_season, 'gather', nbp_glo, index_g)
    IF (ALL(wstress_season(:,:) == val_exp)) THEN
       wstress_season(:,:) = un
       WHERE(veget_max(:,:).LT.min_stomate)
          wstress_season(:,:) = zero
       END WHERE
    ENDIF
    !-
    wstress_month(:,:) = val_exp
    var_name = 'wstress_month'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                .TRUE., wstress_month, 'gather', nbp_glo, index_g)

    IF (ALL(wstress_month(:,:) == val_exp)) THEN
       wstress_month(:,:) = un
       WHERE(veget_max(:,:).LT.min_stomate)
          wstress_month(:,:) = zero
       END WHERE
    ENDIF
    !-
    ! 6 fire probability
    !-
    fireindex(:,:) = val_exp
    var_name = 'fireindex'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &              .TRUE., fireindex, 'gather', nbp_glo, index_g)
    IF (ALL(fireindex(:,:) == val_exp)) fireindex(:,:) = zero
    !-
    firelitter(:,:) = val_exp
    var_name = 'firelitter'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &              .TRUE., firelitter, 'gather', nbp_glo, index_g)
    IF (ALL(firelitter(:,:) == val_exp)) firelitter(:,:) = zero
    !-
    ! 7 maximum and minimum moisture availabilities for tropic phenology
    !-
    maxvegstress_lastyear(:,:) = val_exp
    var_name = 'maxmoistr_last'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., maxvegstress_lastyear, 'gather', nbp_glo, index_g)
    IF (ALL(maxvegstress_lastyear(:,:) == val_exp)) &
         &     maxvegstress_lastyear(:,:) = zero
    !-
    maxvegstress_thisyear(:,:) = val_exp

    var_name = 'maxmoistr_this'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., maxvegstress_thisyear, 'gather', nbp_glo, index_g)
    IF (ALL(maxvegstress_thisyear(:,:) == val_exp)) &
         &     maxvegstress_thisyear(:,:) = zero
    !-
    minvegstress_lastyear(:,:) = val_exp
    var_name = 'minmoistr_last'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., minvegstress_lastyear, 'gather', nbp_glo, index_g)
    IF (ALL(minvegstress_lastyear(:,:) == val_exp)) THEN
       minvegstress_lastyear(:,:) = un
       WHERE(veget_max(:,:).LT.min_stomate)
          minvegstress_lastyear(:,:) = zero
       END WHERE
    END IF
    !-
    minvegstress_thisyear(:,:) = val_exp
    var_name = 'minmoistr_this'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., minvegstress_thisyear, 'gather', nbp_glo, index_g)
    IF (ALL( minvegstress_thisyear(:,:) == val_exp)) &
         &     minvegstress_thisyear(:,:) = large_value
    !-
    ! 8 maximum "weekly" GPP
    !-
    maxgppweek_lastyear(:,:) = val_exp
    var_name = 'maxgppweek_lastyear'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., maxgppweek_lastyear, 'gather', nbp_glo, index_g)
    IF (ALL(maxgppweek_lastyear(:,:) == val_exp)) &
         &     maxgppweek_lastyear(:,:) = zero
    !-
    maxgppweek_thisyear(:,:) = val_exp
    var_name = 'maxgppweek_thisyear'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., maxgppweek_thisyear, 'gather', nbp_glo, index_g)
    IF (ALL(maxgppweek_thisyear(:,:) == val_exp)) &
         &     maxgppweek_thisyear(:,:) = zero
    !-
    ! 9 annual GDD0
    !-
    gdd0_thisyear(:) = val_exp
    var_name = 'gdd0_thisyear'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &              .TRUE., gdd0_thisyear, 'gather', nbp_glo, index_g)
    IF (ALL(gdd0_thisyear(:) == val_exp)) gdd0_thisyear(:) = zero
    !-
    gdd0_lastyear(:) = val_exp
    var_name = 'gdd0_lastyear'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &              .TRUE., gdd0_lastyear, 'gather', nbp_glo, index_g)
    IF (ALL(gdd0_lastyear(:) == val_exp)) gdd0_lastyear(:) = gdd_crit_estab
    !-
    ! 10 annual precipitation
    !-
    precip_thisyear(:) = val_exp
    var_name = 'precip_thisyear'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &              .TRUE., precip_thisyear, 'gather', nbp_glo, index_g)
    IF (ALL(precip_thisyear(:) == val_exp)) precip_thisyear(:) = zero
    !-
    precip_lastyear(:) = val_exp
    var_name = 'precip_lastyear'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &              .TRUE., precip_lastyear, 'gather', nbp_glo, index_g)
    IF (ALL(precip_lastyear(:) == val_exp)) &
         &     precip_lastyear(:) = precip_crit
    !-
    ! 11 derived "biometeorological" variables
    !-
    gdd_m5_dormance(:,:) = val_exp
    var_name = 'gdd_m5_dormance'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., gdd_m5_dormance, 'gather', nbp_glo, index_g)
    IF (ALL(gdd_m5_dormance(:,:) == val_exp)) THEN
       gdd_m5_dormance(:,:) = undef
       WHERE (veget_max(:,:).LT.min_stomate)
          gdd_m5_dormance(:,:) = zero
       END WHERE
    END IF
    !-
    gdd_from_growthinit(:,:) = val_exp
    var_name = 'gdd_from_growthinit'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., gdd_from_growthinit, 'gather', nbp_glo, index_g)
    IF (ALL(gdd_from_growthinit(:,:) == val_exp)) &
         &     gdd_from_growthinit(:,:) = zero
    !-
    gdd_midwinter(:,:) = val_exp
    var_name = 'gdd_midwinter'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., gdd_midwinter, 'gather', nbp_glo, index_g)
    IF (ALL(gdd_midwinter(:,:) == val_exp)) THEN
       gdd_midwinter(:,:) = undef
       WHERE(veget_max(:,:).LT.min_stomate)
          gdd_midwinter(:,:) = zero
       END WHERE
    END IF
    !-
    ncd_dormance(:,:) = val_exp
    var_name = 'ncd_dormance'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., ncd_dormance, 'gather', nbp_glo, index_g)
    IF (ALL(ncd_dormance(:,:) == val_exp)) THEN
       ncd_dormance(:,:) = undef
       WHERE(veget_max(:,:).LT.min_stomate)
          ncd_dormance(:,:) = zero
       END WHERE
    END IF
    !-
    ngd_minus5(:,:) = val_exp
    var_name = 'ngd_minus5'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., ngd_minus5, 'gather', nbp_glo, index_g)
    IF (ALL(ngd_minus5(:,:) == val_exp)) ngd_minus5(:,:) = zero
    !-
    time_hum_min(:,:) = val_exp
    var_name = 'time_hum_min'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., time_hum_min, 'gather', nbp_glo, index_g)
    IF (ALL(time_hum_min(:,:) == val_exp)) THEN
       time_hum_min(:,:) = undef
       WHERE(veget_max(:,:).LT.min_stomate)
          time_hum_min(:,:) = zero
       END WHERE
    END IF
    !-
    hum_min_dormance(:,:) = val_exp
    var_name = 'hum_min_dormance'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., hum_min_dormance, 'gather', nbp_glo, index_g)
    IF (ALL(hum_min_dormance(:,:) == val_exp)) THEN
       hum_min_dormance(:,:) = undef
       WHERE(veget_max(:,:).LT.min_stomate)
          hum_min_dormance(:,:) = zero
       END WHERE
    END IF
    !-
    ! 12 Plant status
    !-
    CALL restget_p (rest_id_stomate, 'PFTpresent', nbp_glo, nvm  , 1, itime, &
         &              .TRUE., PFTpresent_real, 'gather', nbp_glo, index_g)
    IF (ALL(PFTpresent_real(:,:) == val_exp)) PFTpresent_real(:,:) = zero
    WHERE (PFTpresent_real(:,:) >= .5)
       PFTpresent = .TRUE.
    ELSEWHERE
       PFTpresent = .FALSE.
    ENDWHERE

    !CALL restget_p (rest_id_stomate, 'ind', nbp_glo, nvm  , 1, itime, &
    !     &              .TRUE., ind, 'gather', nbp_glo, index_g)
    !IF (ALL(ind(:,:) == val_exp)) ind(:,:) = zero
    !-
    CALL restget_p (rest_id_stomate, 'adapted', nbp_glo, nvm  , 1, itime, &
         &              .TRUE., adapted, 'gather', nbp_glo, index_g)
    IF (ALL(adapted(:,:) == val_exp)) adapted(:,:) = zero
    !-
    CALL restget_p (rest_id_stomate, 'regenerate', nbp_glo, nvm  , 1, itime, &
         &              .TRUE., regenerate, 'gather', nbp_glo, index_g)
    IF (ALL(regenerate(:,:) == val_exp)) regenerate(:,:) = zero
    !-
    CALL restget_p (rest_id_stomate, 'npp_longterm', nbp_glo, nvm  , 1, itime, &
         &              .TRUE., npp_longterm, 'gather', nbp_glo, index_g)
    IF (ALL(npp_longterm(:,:) == val_exp)) npp_longterm(:,:) = zero
    !-
    CALL restget_p (rest_id_stomate, 'croot_longterm', nbp_glo, nvm  , 1, itime, &
         &              .TRUE., croot_longterm, 'gather', nbp_glo, index_g)
    IF (ALL(croot_longterm(:,:) == val_exp)) croot_longterm(:,:) = zero
    !-
    CALL restget_p (rest_id_stomate, 'n_reserve_longterm', nbp_glo, nvm  , 1, itime, &
         &              .TRUE., n_reserve_longterm, 'gather', nbp_glo, index_g)
    IF (ALL(n_reserve_longterm(:,:) == val_exp)) THEN
       n_reserve_longterm(:,:) = un
       WHERE(veget_max(:,:).LT.min_stomate)
          n_reserve_longterm(:,:) = zero
       END WHERE
    END IF
    !-
    CALL restget_p (rest_id_stomate, 'lm_lastyearmax', nbp_glo, nvm  , 1, itime, &
         &              .TRUE., lm_lastyearmax, 'gather', nbp_glo, index_g)
    IF (ALL(lm_lastyearmax(:,:) == val_exp)) lm_lastyearmax(:,:) = zero
    !-
    CALL restget_p (rest_id_stomate, 'lm_thisyearmax', nbp_glo, nvm  , 1, itime, &
         &              .TRUE., lm_thisyearmax, 'gather', nbp_glo, index_g)
    IF (ALL(lm_thisyearmax(:,:) == val_exp)) lm_thisyearmax(:,:) = zero
    !-
    CALL restget_p (rest_id_stomate, 'maxfpc_lastyear', nbp_glo, nvm  , 1, itime, &
         &              .TRUE., maxfpc_lastyear, 'gather', nbp_glo, index_g)
    IF (ALL(maxfpc_lastyear(:,:) == val_exp)) maxfpc_lastyear(:,:) = zero
    !-
    CALL restget_p (rest_id_stomate, 'maxfpc_thisyear', nbp_glo, nvm  , 1, itime, &
         &              .TRUE., maxfpc_thisyear, 'gather', nbp_glo, index_g)
    IF (ALL(maxfpc_thisyear(:,:) == val_exp)) maxfpc_thisyear(:,:) = zero
    !-
    CALL restget_p (rest_id_stomate, 'turnover_time', nbp_glo, nvm  , nparts, itime, &
         &              .TRUE., turnover_time, 'gather', nbp_glo, index_g)
    IF ( ALL( turnover_time(:,:,:) == val_exp)) turnover_time(:,:,:) = 100.
    !-
    CALL restget_p (rest_id_stomate, 'turnover_longterm', nbp_glo, nvm, nparts, nelements, itime, & 
         &              .TRUE., turnover_longterm, 'gather', nbp_glo, index_g) 
    IF (ALL(turnover_longterm == val_exp)) turnover_longterm(:,:,:,:) = zero
    !-
    CALL restget_p (rest_id_stomate, 'gpp_week', nbp_glo, nvm  , 1, itime, &
         &              .TRUE., gpp_week, 'gather', nbp_glo, index_g)
    IF (ALL(gpp_week(:,:) == val_exp)) gpp_week(:,:) = zero
    !-
    CALL restget_p (rest_id_stomate, 'resp_maint_week', nbp_glo, nvm  , 1, itime, &
         &              .TRUE., resp_maint_week, 'gather', nbp_glo, index_g)
    IF (ALL(resp_maint_week(:,:) == val_exp)) resp_maint_week(:,:) = zero
    !-
    CALL restget_p (rest_id_stomate, 'lai_per_level', nbp_glo, nvm, nlevels_tot, itime, &
         &               .TRUE., lai_per_level, 'gather', nbp_glo, index_g)
    IF (ALL(lai_per_level == val_exp)) lai_per_level(:,:,:) = zero
    !
    ! This is ugly, since the restart routines were not developed to deal with structures. 
    ! But all I have to do is copy some things to arrays.
    CALL restget_p (rest_id_stomate, 'laieff_fit', nbp_glo, nvm , nlevels_tot, nparams_laieff, itime, &
          .TRUE., temp_array, 'gather', nbp_glo, index_g)
    IF (ALL(temp_array == val_exp)) temp_array = zero

    DO ipts=1,npts
       DO ivm=1,nvm
          DO ilev=1,nlevels_tot
             laieff_fit(ipts,ivm,ilev)%a=temp_array(ipts,ivm,ilev,1)
             laieff_fit(ipts,ivm,ilev)%b=temp_array(ipts,ivm,ilev,2)
             laieff_fit(ipts,ivm,ilev)%c=temp_array(ipts,ivm,ilev,3)
             laieff_fit(ipts,ivm,ilev)%d=temp_array(ipts,ivm,ilev,4)
             laieff_fit(ipts,ivm,ilev)%e=temp_array(ipts,ivm,ilev,5)
          ENDDO
       ENDDO
    ENDDO
    !
    CALL restget_p (rest_id_stomate, 'maint_resp', nbp_glo, nvm, nparts, itime, &
         &                   .TRUE., resp_maint_part, 'gather', nbp_glo, index_g)
    IF (ALL(resp_maint_part == val_exp)) resp_maint_part(:,:,:) = zero
    !-
    CALL restget_p (rest_id_stomate, 'leaf_age', nbp_glo, nvm, nleafages, itime, &
         &                   .TRUE., leaf_age(:,:,:), 'gather', nbp_glo, index_g)
    IF (ALL(leaf_age == val_exp)) leaf_age(:,:,:) = zero
    !-
    CALL restget_p (rest_id_stomate, 'leaf_frac', nbp_glo, nvm, nleafages, itime, &
         &                  .TRUE., leaf_frac(:,:,:), 'gather', nbp_glo, index_g)
    IF (ALL(leaf_frac(:,:,:) == val_exp)) leaf_frac(:,:,:) = zero
    !-
    ! leaf_age_crit depends on t2m_longterm and longevity_leaf. At 
    ! the initialization phase t2m_longterm is not yet known so, leag_age_crit 
    ! cannot yet be calculated (see stomate_season.f90). A work around is used 
    ! to guess leaf_age_crit for the first time step.
    CALL restget_p (rest_id_stomate, 'leaf_age_crit', nbp_glo, nvm, 1, itime, &
         &                  .TRUE., leaf_age_crit(:,:), 'gather', nbp_glo, index_g)
    IF (ALL(leaf_age_crit(:,:) == val_exp)) THEN
       leaf_age_crit(:,1) = zero
       DO ivm = 2,nvm
          leaf_age_crit(:,:) = longevity_leaf(ivm)
       END DO
    END IF
    !-
    CALL restget_p (rest_id_stomate, 'plant_status', nbp_glo, nvm  , 1, itime, &
         &                .TRUE., plant_status, 'gather', nbp_glo, index_g)
    IF (ALL(plant_status(:,:) == val_exp)) THEN
       plant_status(:,:) = iprescribe
       WHERE(veget_max(:,:).LT.min_stomate)
          plant_status(:,:) = inone
       END WHERE
       plant_status(:,ibare_sechiba) = inone
    ENDIF

    !-
    CALL restget_p (rest_id_stomate, 'when_growthinit', nbp_glo, nvm  , 1, itime, &
         &                .TRUE., when_growthinit, 'gather', nbp_glo, index_g)
    IF (ALL(when_growthinit(:,:) == val_exp)) when_growthinit(:,:) = 240.
    !-
    CALL restget_p (rest_id_stomate, 'age', nbp_glo, nvm  , 1, itime, &
         &                .TRUE., age, 'gather', nbp_glo, index_g)
    IF (ALL(age(:,:) == val_exp)) age(:,:) = zero
    !-
    ! 13 CO2
    !-
    CALL restget_p (rest_id_stomate, 'resp_hetero', nbp_glo, nvm, 1, itime, &
         &                .TRUE., resp_hetero, 'gather', nbp_glo, index_g)
    IF (ALL(resp_hetero(:,:) == val_exp)) resp_hetero(:,:) = zero
    !-
    CALL restget_p (rest_id_stomate, 'resp_maint', nbp_glo, nvm  , 1, itime, &
         &                .TRUE., resp_maint, 'gather', nbp_glo, index_g)
    IF (ALL(resp_maint(:,:) == val_exp)) resp_maint(:,:) = zero
    !-
    CALL restget_p (rest_id_stomate, 'resp_growth', nbp_glo, nvm  , 1, itime, &
         &                .TRUE., resp_growth, 'gather', nbp_glo, index_g)
    IF (ALL(resp_growth(:,:) == val_exp)) resp_growth(:,:) = zero
    !-
    co2_fire(:,:) = val_exp
    var_name = 'co2_fire'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &                .TRUE., co2_fire, 'gather', nbp_glo, index_g)
    IF (ALL(co2_fire(:,:) == val_exp)) co2_fire(:,:) = zero
    !-
    ! 14 vegetation distribution after last light competition
    !-
    var_name = 'veget_lastlight'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &                .TRUE., veget_lastlight, 'gather', nbp_glo, index_g)
    IF (ALL(veget_lastlight(:,:) == val_exp)) veget_lastlight(:,:) = zero
    !-
    ! 15 establishment criteria
    !-
    var_name = 'everywhere'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &                .TRUE., everywhere, 'gather', nbp_glo, index_g)
    IF (ALL(everywhere(:,:) == val_exp)) everywhere(:,:) = zero
    !-
    var_name = 'need_adjacent'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &                .TRUE., need_adjacent_real, 'gather', nbp_glo, index_g)
    IF (ALL(need_adjacent_real(:,:) == val_exp)) &
         &     need_adjacent_real(:,:) = zero
    WHERE ( need_adjacent_real(:,:) >= .5 )
       need_adjacent = .TRUE.
    ELSEWHERE
       need_adjacent = .FALSE.
    ENDWHERE
    !-
    var_name = 'RIP_time'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &                .TRUE., RIP_time, 'gather', nbp_glo, index_g)
    IF (ALL(RIP_time(:,:) == val_exp)) RIP_time(:,:) = large_value
    !-
    ! 17 litter
    !-
    !litter(:,:,:,:,:) = val_exp
    CALL restget_p (rest_id_stomate, 'litter_c', nbp_glo, nlitt, nvm, nlevs, itime, &
         &                     .TRUE., litter(:,:,:,:,icarbon), 'gather', nbp_glo, index_g)
    IF (ALL(litter(:,:,:,:,icarbon) == val_exp)) THEN
       litter(:,:,:,:,icarbon) = zero
    ENDIF

    CALL restget_p (rest_id_stomate, 'litter_n', nbp_glo, nlitt, nvm, nlevs,itime, &
         &                     .TRUE., litter(:,:,:,:,initrogen), 'gather', nbp_glo,index_g)
    IF (ALL(litter(:,:,:,:,initrogen) == val_exp)) THEN
       litter(:,:,:,:,initrogen) = zero
    ENDIF
    
    !spitfire
    litterfuel(:,:,:,:,:) = val_exp
    CALL restget_p (rest_id_stomate, 'litterfuel', nbp_glo, nlitt, nvm, nhour, nelements, itime, &
         &                     .TRUE., litterfuel(:,:,:,:,:), 'gather', nbp_glo, index_g)
    IF (ALL(litterfuel(:,:,:,:,:) == val_exp)) THEN
       litterfuel(:,:,:,:,:)= zero
    ENDIF
    
    CALL restget_p (rest_id_stomate, 'dead_leaves', nbp_glo, nvm, nlitt, itime, &
          &                   .TRUE., dead_leaves, 'gather', nbp_glo, index_g)
    IF (ALL(dead_leaves == val_exp)) dead_leaves = zero

    CALL restget_p (rest_id_stomate, 'soil_carbon', nbp_glo, ncarb, nvm, itime, &
         &                   .TRUE., som(:,:,:,icarbon), 'gather', nbp_glo, index_g) 
    ! Initialise the soil pools. 
    IF (ALL(som(:,:,:,icarbon) == val_exp)) THEN
       ! We set the pools to zero for all PFTs
       som(:,:,:,icarbon) = zero
       ! Initialize all PFTs that are present. Note that PFT1 is
       ! left at zero.
       DO ivm = 2,nvm
          WHERE (veget_max(:,ivm) .GT. min_stomate) 
             ! Set the initial values for PFTs that are present.  
             som(:,iactive,ivm,icarbon) = som_init_active 
             som(:,isurface,ivm,icarbon) = som_init_surface 
             som(:,islow,ivm,icarbon) = som_init_slow
             som(:,ipassive,ivm,icarbon) = som_init_passive
          END WHERE
       ENDDO  
    ENDIF

    CALL restget_p (rest_id_stomate, 'soil_nitrogen', nbp_glo, ncarb, nvm, itime, & 
         &                   .TRUE., som(:,:,:,initrogen), 'gather', nbp_glo, index_g) 
    IF (ALL(som(:,:,:,initrogen) == val_exp)) THEN 
       som(:,iactive,:,initrogen) = som(:,iactive,:,icarbon) / CN_target_iactive_ref 
       som(:,isurface,:,initrogen) = som(:,isurface,:,icarbon) / CN_target_isurface_ref 
       som(:,islow,:,initrogen) = som(:,islow,:,icarbon) / CN_target_islow_ref 
       som(:,ipassive,:,initrogen) =  som(:,ipassive,:,icarbon) / CN_target_ipassive_ref 
    ENDIF
    CALL restget_p (rest_id_stomate, 'lignin_struc', nbp_glo, nvm, nlevs, itime, &
         &     .TRUE., lignin_struc, 'gather', nbp_glo, index_g)
    IF (ALL(lignin_struc == val_exp)) lignin_struc(:,:,:) = zero

    CALL restget_p (rest_id_stomate, 'lignin_wood', nbp_glo, nvm, nlevs, itime, &
         &     .TRUE., lignin_wood, 'gather', nbp_glo, index_g)
    IF (ALL(lignin_wood == val_exp)) lignin_wood(:,:,:) = zero
    
    CALL restget_p (rest_id_stomate, 'lignin_snag', nbp_glo, nvm, nlevs, itime, &
         &     .TRUE., lignin_snag, 'gather', nbp_glo, index_g)
    IF (ALL(lignin_snag == val_exp)) lignin_snag(:,:,:) = zero


    ! 18 Product use and LCC
    !-
    var_name = 'fco2_lu'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
         .TRUE., fco2_lu, 'gather', nbp_glo, index_g)
    IF (ALL(fco2_lu(:) == val_exp)) fco2_lu(:) = zero

    var_name = 'fco2_wh'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
         .TRUE., fco2_wh, 'gather', nbp_glo, index_g)
    IF (ALL(fco2_wh(:) == val_exp)) fco2_wh(:) = zero

    var_name = 'fco2_ha'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
         .TRUE., fco2_ha, 'gather', nbp_glo, index_g)
    IF (ALL(fco2_ha(:) == val_exp)) fco2_ha(:) = zero

       
    IF (vegetmap_reset) THEN

       ! Reset vegetation map related variables instead of reading from restart file
       ! vegetmap_reset is an option to change vegetation map without activating LAND 
       ! USE change for carbon fluxes. At the same time carbon related variables are 
       ! reset to zero. Use this option to change vegetation map while keeping 
       ! VEGET_UPDATE=0Y.

       ! ChaoYue: Be very careful for these following lines though. The First 6
       ! variables all carry legacy information on harvest or land-use-change 
       ! related pools. Simply reseting them to zero might lead to mass balance
       ! problem.
       prod_s(:,:,:,:,:) = zero
       prod_m(:,:,:,:,:) = zero
       prod_l(:,:,:,:,:) = zero
       flux_s(:,:,:,:,:) = zero
       flux_m(:,:,:,:,:) = zero
       flux_l(:,:,:,:,:) = zero

    ELSE

       DO l = 1,nlctypes

          var_name = 'prod_s'//TRIM(lctype_str(l))
          CALL restget_p (rest_id_stomate, var_name, nbp_glo, nshort+1, nelements, nlanduse, itime, &
               &                .TRUE., prod_s(:,:,:,:,l), 'gather', nbp_glo, index_g)
          IF (ALL(prod_s(:,:,:,:,l) == val_exp)) prod_s(:,:,:,:,l) = zero

          var_name = 'prod_m'//TRIM(lctype_str(l))
          CALL restget_p (rest_id_stomate, var_name, nbp_glo, nmedium+1, nelements, nlanduse, itime, &
               &                .TRUE., prod_m(:,:,:,:,l), 'gather', nbp_glo, index_g)
          IF (ALL(prod_m(:,:,:,:,l) == val_exp)) prod_m(:,:,:,:,l) = zero
          
          var_name = 'prod_l'//TRIM(lctype_str(l))
          CALL restget_p (rest_id_stomate, var_name, nbp_glo, nlong+1, nelements, nlanduse, itime, &
               &                .TRUE., prod_l(:,:,:,:,l), 'gather', nbp_glo, index_g)
          IF (ALL(prod_l(:,:,:,:,l) == val_exp)) prod_l(:,:,:,:,l) = zero
          
          var_name = 'flux_s'//TRIM(lctype_str(l))
          CALL restget_p (rest_id_stomate, var_name, nbp_glo, nshort+1 , nelements, nlanduse, itime, &
               &                .TRUE., flux_s(:,:,:,:,l), 'gather', nbp_glo, index_g)
          IF (ALL(flux_s(:,:,:,:,l) == val_exp)) flux_s(:,:,:,:,l) = zero
          
          var_name = 'flux_m'//TRIM(lctype_str(l))
          CALL restget_p (rest_id_stomate, var_name, nbp_glo, nmedium+1, nelements, nlanduse, itime, &
            &                .TRUE., flux_m(:,:,:,:,l), 'gather', nbp_glo, index_g)
          IF (ALL(flux_m(:,:,:,:,l) == val_exp)) flux_m(:,:,:,:,l) = zero
          
          var_name = 'flux_l'//TRIM(lctype_str(l))
          CALL restget_p (rest_id_stomate, var_name, nbp_glo, nlong+1, nelements, nlanduse, itime, &
               &                .TRUE., flux_l(:,:,:,:,l), 'gather', nbp_glo, index_g)
          IF (ALL(flux_l(:,:,:,:,l) == val_exp)) flux_l(:,:,:,:,l) = zero

       END DO

    END IF  ! vegetmap_reset


    fDeforestToProduct(:,:) = val_exp
    var_name = 'fDeforestToProduct'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &   .TRUE., fDeforestToProduct, 'gather', nbp_glo, index_g)
    IF (ALL(fDeforestToProduct(:,:) ==val_exp)) fDeforestToProduct(:,:) = zero

    fLulccResidue(:,:) = val_exp
    var_name = 'fLulccResidue'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &   .TRUE., fLulccResidue, 'gather', nbp_glo, index_g)
    IF (ALL(fLulccResidue(:,:) ==val_exp)) fLulccResidue(:,:) = zero

    fHarvestToProduct(:,:) = val_exp
    var_name = 'fHarvestToProduct'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &   .TRUE., fHarvestToProduct, 'gather', nbp_glo, index_g)
    IF (ALL(fHarvestToProduct(:,:) ==val_exp)) fHarvestToProduct(:,:) = zero

    !-
    bm_to_litter(:,:,:,:) = val_exp
    CALL restget_p (rest_id_stomate, 'bm_to_litter', nbp_glo, nvm, nparts, nelements, itime, &
         &                .TRUE., bm_to_litter, 'gather', nbp_glo, index_g)
    IF (ALL(bm_to_litter == val_exp)) bm_to_litter(:,:,:,:) = zero
    !-
    bm_to_litter_resid(:,:,:,:) = val_exp
    CALL restget_p (rest_id_stomate, 'bm_to_litter_resid', nbp_glo, nvm, nparts, nelements, itime, &
         &                .TRUE., bm_to_litter_resid, 'gather', nbp_glo, index_g)
    IF (ALL(bm_to_litter_resid == val_exp)) bm_to_litter_resid(:,:,:,:) = zero
    !-
    tree_bm_to_litter(:,:,:,:) = val_exp
    CALL restget_p (rest_id_stomate, 'tree_to_litter', nbp_glo, nvm, nparts, nelements, itime, &
          &                .TRUE., tree_bm_to_litter, 'gather', nbp_glo, index_g)
    IF (ALL(tree_bm_to_litter == val_exp)) tree_bm_to_litter = zero
    !-
    tree_bm_to_litter_resid(:,:,:,:) = val_exp
    CALL restget_p (rest_id_stomate, 'tree_to_litter_resid', nbp_glo, nvm, nparts, nelements, itime, &
          &                .TRUE., tree_bm_to_litter_resid, 'gather', nbp_glo, index_g)
    IF (ALL(tree_bm_to_litter_resid == val_exp)) tree_bm_to_litter_resid = zero
    !- 
    carb_mass_total(:) = val_exp
    var_name = 'carb_mass_total'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &              .TRUE., carb_mass_total, 'gather', nbp_glo, index_g)
    IF (ALL(carb_mass_total(:) == val_exp)) carb_mass_total(:) = zero
    !-
    ! Harvest_pool_acc is a 5D variable. Write separatly for c and N
    harvest_pool_acc(:,:,:,:,:) = val_exp
    harvest_area_acc(:,:,:) = val_exp
    DO l = 1,nlanduse
       ! 
       var_name ='harvest_pool_acc'//TRIM(luse_str(l))
       CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm, ndia_harvest+1, nelements, itime, &
         &    .TRUE., harvest_pool_acc(:,:,:,:,l), 'gather', nbp_glo, index_g)
       IF ( ALL(harvest_pool_acc(:,:,:,:,l) == val_exp) ) harvest_pool_acc(:,:,:,:,l) = 0.
       !-
       var_name ='harvest_area_acc'//TRIM(luse_str(l))
       CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
            &                .TRUE., harvest_area_acc(:,:,l), 'gather', nbp_glo, index_g)
       IF (ALL(harvest_area_acc(:,:,l) == val_exp)) harvest_area_acc(:,:,l) = zero
    ENDDO
    !-
    IF ( ok_soil_carbon_discretization ) THEN
       deepSom_a(:,:,:,:) = val_exp
       CALL restget_p (rest_id_stomate, 'deepSOM_a', nbp_glo, ngrnd, nvm, nelements,  itime, &
            .TRUE., deepSOM_a, 'gather', nbp_glo, index_g)

       deepSom_s(:,:,:,:) = val_exp
       CALL restget_p (rest_id_stomate, 'deepSOM_s', nbp_glo, ngrnd, nvm, nelements, itime, &
            .TRUE., deepSOM_s, 'gather', nbp_glo, index_g)
       deepSom_p(:,:,:,:) = val_exp
       CALL restget_p (rest_id_stomate, 'deepSOM_p', nbp_glo, ngrnd, nvm, nelements, itime, &
            .TRUE., deepSOM_p, 'gather', nbp_glo, index_g)

       ! Initialize if not found in the restart file
       IF (ALL(deepSOM_a == val_exp) .OR. ALL(deepSOM_s == val_exp) .OR. ALL(deepSOM_p == val_exp)) THEN
          ! One or several of the variables were not found in the restart file
          
          IF (use_initsom) THEN
             IF (printlev >= 2) WRITE(numout,*) 'Reading initSOM file for initializing soil organic matter'
             CALL read_initsoc_file(npts,lalo,neighbours, resolution, contfrac, veget_max, deepSOM_a(:,:,:,icarbon), deepSOM_s(:,:,:,icarbon), deepSOM_p(:,:,:,icarbon)) 
             CALL read_initson_file(npts,lalo,neighbours, resolution, contfrac, veget_max, deepSOM_a(:,:,:,icarbon), deepSOM_s(:,:,:,icarbon), deepSOM_p(:,:,:,icarbon), &
                                    deepSOM_a(:,:,:,initrogen), deepSOM_s(:,:,:,initrogen), deepSOM_p(:,:,:,initrogen))
          ELSE
             IF (printlev >= 2) WRITE(numout,*) 'Initializing soil organic matter using scalar values'
             deepSOM_a(:,:,:,icarbon) = 10.
             deepSOM_s(:,:,:,icarbon) = 30.
             deepSOM_p(:,:,:,icarbon) = 50.
             deepSOM_a(:,:,:,initrogen) = deepSOM_a(:,:,:,icarbon) / CN_target_iactive_ref
             deepSOM_s(:,:,:,initrogen) = deepSOM_s(:,:,:,icarbon) / CN_target_islow_ref
             deepSOM_p(:,:,:,initrogen) = deepSOM_p(:,:,:,icarbon) / CN_target_ipassive_ref
          ENDIF
       END IF

       CALL restget_p (rest_id_stomate, 'O2_soil', nbp_glo, ngrnd, nvm, itime, &
            .TRUE., O2_soil, 'gather', nbp_glo, index_g)
       IF (ALL(O2_soil == val_exp)) O2_soil = O2_init_conc

       CALL restget_p (rest_id_stomate,'CH4_soil', nbp_glo, ngrnd, nvm, itime, &
            .TRUE., CH4_soil, 'gather', nbp_glo, index_g)
       IF (ALL(CH4_soil == val_exp)) CH4_soil = CH4_init_conc

       CALL restget_p (rest_id_stomate, 'O2_snow', nbp_glo, nsnow, nvm, itime, &
            .TRUE., O2_snow, 'gather', nbp_glo, index_g)
       IF (ALL(O2_snow == val_exp)) O2_snow = O2_init_conc

       CH4_snow(:,:,:) = val_exp
       CALL restget_p (rest_id_stomate,'CH4_snow', nbp_glo, nsnow, nvm, itime, &
            .TRUE., CH4_snow, 'gather', nbp_glo, index_g)
       IF (ALL(CH4_snow == val_exp)) CH4_snow = CH4_init_conc
   
       CALL restget_p (rest_id_stomate,'heat_Zimov', nbp_glo, ngrnd, nvm, itime, &
            .TRUE., heat_Zimov, 'gather', nbp_glo, index_g)
       IF (ALL(heat_Zimov == val_exp)) THEN
          heat_Zimov(:,:,:) = 0.0
       ENDIF

       CALL restget_p (rest_id_stomate,'altmax_ind', nbp_glo, nvm, 1, itime, &
            .TRUE., altmax_ind, 'gather', nbp_glo, index_g)
       IF (ALL(altmax_ind(:,:) == val_exp)) THEN
          altmax_ind(:,:) = 0.0
       ENDIF

       CALL restget_p (rest_id_stomate,'altmax_lastyear', nbp_glo, nvm, 1, itime, &
            .TRUE., altmax_lastyear, 'gather', nbp_glo, index_g)
       IF (ALL(altmax_lastyear(:,:) == val_exp)) THEN
          altmax_lastyear(:,:) = 0.0
       ENDIF

       CALL restget_p (rest_id_stomate,'altmax_ind_lastyear', nbp_glo, nvm, 1, itime, &
            .TRUE., altmax_ind_lastyear, 'gather', nbp_glo, index_g)
       IF (ALL(altmax_ind_lastyear(:,:) == val_exp)) THEN
          altmax_ind_lastyear(:,:) = 0.0
       ENDIF

       CALL restget_p (rest_id_stomate,'depth_organic_soil', nbp_glo, 1, 1, itime, &
           .TRUE., depth_organic_soil(:), 'gather', nbp_glo, index_g)
       IF (ALL(depth_organic_soil(:) == val_exp)) THEN
           depth_organic_soil(:) = 0.0
           read_input_depth_organic_soil = .TRUE.
       ENDIF

       fixed_cryoturbation_depth(:,:) = val_exp
       CALL restget_p (rest_id_stomate,'fixed_cryoturb_depth', nbp_glo, nvm, 1, itime, &
            .TRUE., fixed_cryoturbation_depth, 'gather', nbp_glo, index_g)
       IF (ALL(fixed_cryoturbation_depth(:,:) == val_exp)) THEN
          fixed_cryoturbation_depth(:,:) = 0.0
       ENDIF
       
    ENDIF

    CALL restget_p (rest_id_stomate,'altmax', nbp_glo, nvm, 1, itime, &
         .TRUE., altmax, 'gather', nbp_glo, index_g)
    IF (ALL(altmax(:,:) == val_exp)) THEN
       IF ( ok_soil_carbon_discretization ) THEN
          ! altmax will be calculated in stomate_soil_carbon_discretization
          altmax(:,:) = 0.0
       ELSE
          ! altmax will not be calculated in this configuration but it is used
          ! in hydrol_root_profile. Hence, altmax needs to get a value that does
          ! not affect the calculation of the root profile.
          altmax(:,:) = zdr(nslm)
       ENDIF
    ENDIF
    
    CALL restget_p (rest_id_stomate, 'nbp_accu_flux', nbp_glo, nelements, 1, itime, &
         &     .TRUE., nbp_accu_flux, 'gather', nbp_glo, index_g)
    IF (ALL(nbp_accu_flux == val_exp)) THEN
       ! There is no restart value which means that the model is at its very first
       ! time step. Most likley the user has seeded C and N in the soil to speed
       ! up the spinup. This C and N basically bypasses GPP and biological nitrogen
       ! fixation. As it once entered the ecosystems (although in the model-world
       ! we don really how and when) it should be accounted for in the consistency
       ! and mass balance checks.
       IF (ok_soil_carbon_discretization) THEN
          DO iele = 1,nelements
             nbp_accu_flux(:,iele) = zero
             DO igrn = 1,ngrnd
                nbp_accu_flux(:,iele) = nbp_accu_flux(:,iele) + &
                     SUM((deepSOM_a(:,igrn,:,iele) + deepSOM_s(:,igrn,:,iele) + &
                     deepSOM_p(:,igrn,:,iele)) * &
                     (zf_soil(igrn)-zf_soil(igrn-1)) * veget_max(:,:),2)
             END DO
          END DO
       ELSE
          DO iele = 1,nelements
             nbp_accu_flux(:,iele) = SUM(SUM(som(:,:,:,iele),2)*veget_max(:,:),2)
          END DO
       ENDIF
    ENDIF

    ! nbp_pool_start and nbp_accu_flux should be identical. They are
    ! kept as independent variables such that this identity can be checked for
    ! in stomate.f90. If we would use a single variable the flux based nbp would
    ! contain some pool-based values or the other way around. We want to have
    ! two completely independent ways of calculating nbp.
    nbp_pool_start(:,:) = val_exp
    var_name = 'nbp_pool_start'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nelements, 1, itime, &
         &              .TRUE., nbp_pool_start(:,:), 'gather', nbp_glo, index_g)
    IF (ALL(nbp_pool_start(:,:) == val_exp))THEN
       ! There is no restart value which means that the model is at its very first
       ! time step. Most likley the user has seeded C and N in the soil to speed
       ! up the spinup. This C and N basically bypasses GPP and biological nitrogen
       ! fixation. As it once entered the ecosystems (although in the model-world
       ! we don really how and when) it should be accounted for in the consistency
       ! and mass balance checks.
       IF (ok_soil_carbon_discretization) THEN
          DO iele = 1,nelements
             nbp_pool_start(:,iele)=zero
             DO igrn = 1,ngrnd
                nbp_pool_start(:,iele) = nbp_pool_start(:,iele) + &
                     SUM((deepSOM_a(:,igrn,:,iele) + deepSOM_s(:,igrn,:,iele) + &
                     deepSOM_p(:,igrn,:,iele)) * &
                     (zf_soil(igrn)-zf_soil(igrn-1)) * veget_max(:,:),2)
             END DO
          END DO
       ELSE
          DO iele = 1,nelements
             nbp_pool_start(:,iele) = SUM(SUM(som(:,:,:,iele),2)*veget_max(:,:),2)
          END DO
       ENDIF
    ENDIF

    sugar_load(:,:) = val_exp
    var_name = 'sugar_load'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &   .TRUE., sugar_load, 'gather', nbp_glo, index_g)
    IF (ALL(sugar_load(:,:) ==val_exp)) THEN
       sugar_load(:,:) = un
       WHERE(veget_max(:,:).LT.min_stomate)
          sugar_load(:,:) = zero
       END WHERE
    END IF

    harvest_cut(:,:) = val_exp
    var_name = 'harvest_cut'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &   .TRUE., harvest_cut, 'gather', nbp_glo, index_g)
    IF (ALL(harvest_cut(:,:) ==val_exp)) harvest_cut(:,:) = zero
    
    burried_litter(:,:,:,:) = val_exp
    var_name = 'burried_litter'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nlitt, nlevs, nelements, itime, &
         &   .TRUE., burried_litter, 'gather', nbp_glo, index_g)
    IF (ALL(burried_litter(:,:,:,:) ==val_exp)) burried_litter(:,:,:,:) = zero

    burried_fresh_ltr(:,:,:) = val_exp
    var_name = 'burried_fresh_ltr'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nparts, nelements, itime, &
         &   .TRUE., burried_fresh_ltr, 'gather', nbp_glo, index_g)
    IF (ALL(burried_fresh_ltr(:,:,:) ==val_exp)) burried_fresh_ltr(:,:,:) = zero

    burried_fresh_som(:,:,:) = val_exp
    var_name = 'burried_fresh_som'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nparts, nelements, itime, &
         &   .TRUE., burried_fresh_som, 'gather', nbp_glo, index_g)
    IF (ALL(burried_fresh_som(:,:,:) ==val_exp)) burried_fresh_som(:,:,:) = zero

    burried_bact(:) = val_exp
    var_name = 'burried_bact'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
         &   .TRUE., burried_bact, 'gather', nbp_glo, index_g)
    IF (ALL(burried_bact(:) ==val_exp)) burried_bact(:) = zero
    
    !spitfire
    ni_acc(:) = val_exp
    var_name = 'ni_acc'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
         &   .TRUE., ni_acc, 'gather', nbp_glo, index_g)
    IF (ALL(ni_acc(:) == val_exp)) ni_acc(:) = zero
    
    burried_min_nitro(:,:) = val_exp
    var_name = 'burried_min_nitro'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nnspec, 1, itime, &
         &   .TRUE., burried_min_nitro, 'gather', nbp_glo, index_g)
    IF (ALL(burried_min_nitro(:,:) ==val_exp)) burried_min_nitro(:,:) = zero

    burried_som(:,:,:) = val_exp
    var_name = 'burried_som'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, ncarb, nelements, itime, &
         &   .TRUE., burried_som, 'gather', nbp_glo, index_g)
    IF (ALL(burried_som(:,:,:) ==val_exp)) burried_som(:,:,:) = zero

    burried_deepSOM_a(:,:,:) = val_exp
    var_name = 'burried_deepSOM_a'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, ngrnd, nelements, itime, &
         &   .TRUE., burried_deepSOM_a, 'gather', nbp_glo, index_g)
    IF (ALL(burried_deepSOM_a(:,:,:) ==val_exp)) burried_deepSOM_a(:,:,:) = zero

    burried_deepSOM_s(:,:,:) = val_exp
    var_name = 'burried_deepSOM_s'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, ngrnd, nelements, itime, &
         &   .TRUE., burried_deepSOM_s, 'gather', nbp_glo, index_g)
    IF (ALL(burried_deepSOM_s(:,:,:) ==val_exp)) burried_deepSOM_s(:,:,:) = zero

    burried_deepSOM_p(:,:,:) = val_exp
    var_name = 'burried_deepSOM_p'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, ngrnd, nelements, itime, &
         &   .TRUE., burried_deepSOM_p, 'gather', nbp_glo, index_g)
    IF (ALL(burried_deepSOM_p(:,:,:) ==val_exp)) burried_deepSOM_p(:,:,:) = zero

    ! If the variable is not in the restart file, then zero will be used as default value
    CALL restget_p (rest_id_stomate, 'Global_years', itime, .TRUE., zero, global_years)

    !-
    ! 19. Spinup
    !-
    !MatrixV(:,:,:,:) = 1.0
    !VectorU(:,:,:) = 0.0
    !previous_stock(:,:,:) = 0.0
    !current_stock(:,:,:) = 0.0
    !tau_CN_longterm = 0.0
    !total_ba_init(:,:) = 0.0
    !CN_som_litter_longterm(:,:,:)= 0.0

    IF (spinup_analytic) THEN

       !-
       ok_equilibrium_real(:) = val_exp
       var_name = 'ok_equilibrium'
       CALL restget_p (rest_id_stomate, var_name, nbp_glo , 1  , 1, itime, &
            &                .TRUE., ok_equilibrium_real,'gather', nbp_glo, index_g)
       IF (ALL(ok_equilibrium_real(:) == val_exp)) ok_equilibrium_real(:) = zero
       WHERE(ok_equilibrium_real(:) >= 0.5) 
          ok_equilibrium = .TRUE.
       ELSEWHERE
          ok_equilibrium = .FALSE.
       ENDWHERE
       !-
       MatrixV(:,:,:,:) = val_exp
       CALL restget_p (rest_id_stomate, 'MatrixV', nbp_glo, nvm, nbpools, nbpools, itime, &
                  &                     .TRUE., MatrixV, 'gather', nbp_glo, index_g)
       ! If nothing is found in the restart file, we initialize each submatrix by identity
       IF (ALL(MatrixV(:,:,:,:) == val_exp))  THEN 
          MatrixV(:,:,:,:) = zero
          DO l = 1,nbpools
             MatrixV(:,:,l,l) = un
          END DO
       END IF

       VectorU(:,:,:)  = val_exp
       CALL restget_p &
            &    (rest_id_stomate, 'Vector_U', nbp_glo, nvm, nbpools, itime, &
            &     .TRUE., VectorU, 'gather', nbp_glo, index_g)
       IF (ALL(VectorU == val_exp))  VectorU = zero

       previous_stock(:,:,:)  = val_exp
       CALL restget_p &
            &    (rest_id_stomate, 'previous_stock', nbp_glo, nvm, nbpools, itime, &
            &     .TRUE., previous_stock, 'gather', nbp_glo, index_g)
       IF (ALL(previous_stock == val_exp))  previous_stock = undef_sechiba

       current_stock(:,:,:)  = val_exp
       CALL restget_p &
            &    (rest_id_stomate, 'current_stock', nbp_glo, nvm, nbpools, itime, &
            &     .TRUE., current_stock, 'gather', nbp_glo, index_g)
       IF (ALL(current_stock == val_exp))  current_stock = zero

       CN_som_litter_longterm(:,:,:)  = val_exp
       CALL restget_p &
            &    (rest_id_stomate, 'CN_longterm', nbp_glo, nvm, nbpools, itime, &
            &     .TRUE., CN_som_litter_longterm, 'gather', nbp_glo, index_g)
       IF (ALL(CN_som_litter_longterm == val_exp))  CN_som_litter_longterm = zero

       ! If the variable is not in the restart file, then dt_sechiba/one_day will 
       ! be used as default value
       CALL restget_p(rest_id_stomate, 'tau_CN_longterm', itime, &
            .TRUE., dt_sechiba/one_day, tau_CN_longterm)
 
    ENDIF ! spinup_matrix_method

    KF(:,:) = val_exp
    CALL restget_p (rest_id_stomate, 'KF', nbp_glo, nvm  , 1, itime, &
         &                .TRUE., KF, 'gather', nbp_glo, index_g)
    IF (ALL(KF(:,:) == val_exp)) KF(:,:) = zero

    k_latosa_adapt(:,:) = val_exp
    CALL restget_p (rest_id_stomate, 'k_latosa_adapt', nbp_glo, nvm  , 1, itime, &
         &                .TRUE., k_latosa_adapt, 'gather', nbp_glo, index_g)
    DO m = 1,nvm
       IF (ALL(k_latosa_adapt(:,m) == val_exp)) THEN
          k_latosa_adapt(:,m) = k_latosa_min(m)
          WHERE (veget_max(:,m).LT.min_stomate)
             k_latosa_adapt(:,m) = zero
          END WHERE
       END IF
    ENDDO

    rue_longterm(:,:) = val_exp
    var_name = 'rue_longterm'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm , 1, itime, &
         &        .TRUE., rue_longterm(:,:), 'gather', nbp_glo, index_g)
    IF (ALL(rue_longterm(:,:) == val_exp)) THEN
       rue_longterm(:,:) = un
       WHERE(veget_max(:,:).LT.min_stomate)
          rue_longterm(:,:) = zero
       END WHERE
    END IF

    cn_leaf_min_season(:,:) = val_exp 
    var_name = 'cn_leaf_min_season' 
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, & 
         &              .TRUE., cn_leaf_min_season, 'gather', nbp_glo, index_g) 
    IF ( ALL(cn_leaf_min_season(:,:) == val_exp) ) THEN 
       cn_leaf_min_season(:,:) = cn_leaf_init_2D(:,:) 
       WHERE(veget_max(:,:).LT.min_stomate)
          cn_leaf_min_season(:,:) = zero
       END WHERE
    ENDIF
    
    nstress_season(:,:) = val_exp 
    var_name = 'nstress_season' 
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, & 
         &              .TRUE., nstress_season, 'gather', nbp_glo, index_g) 
    IF ( ALL(nstress_season(:,:) == val_exp) ) THEN
       nstress_season(:,:) = un
       WHERE(veget_max(:,:).LT.min_stomate)
          nstress_season(:,:) = zero
       END WHERE
    END IF
    
    CALL restget_p (rest_id_stomate, 'soil_n_min', nbp_glo, nvm, nnspec, itime, & 
         &              .TRUE., soil_n_min, 'gather', nbp_glo, index_g) 
!    IF ( ALL(soil_n_min == val_exp) ) soil_n_min(:,:,:)=100000000. 
    IF ( ALL(soil_n_min == val_exp) ) soil_n_min(:,:,:)=0. 

    p_O2(:,:) = val_exp 
    var_name = 'p_O2'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, & 
         &              .TRUE., p_O2(:,:), 'gather', nbp_glo, index_g) 
    IF ( ALL(p_O2(:,:) == val_exp) ) THEN
       p_O2(:,:)=200.
       WHERE(veget_max(:,:).LT.min_stomate)
          p_O2(:,:) = zero
       END WHERE
    END IF

    bact(:,:) = val_exp 
    var_name = 'bact'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, & 
         &              .TRUE., bact(:,:), 'gather', nbp_glo, index_g) 
    IF ( ALL(bact(:,:) == val_exp) ) THEN
       bact(:,:)=10
       WHERE(veget_max(:,:).LT.min_stomate)
          bact(:,:) = zero
       END WHERE
    END IF

    var_name = 'age_stand'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., temp_real, 'gather', nbp_glo, index_g)
    IF ( ALL(temp_real(:,:) == val_exp) ) THEN
       age_stand(:,:) = 0
    ELSE
       age_stand = NINT(temp_real)
    ENDIF
    !-
    ! STAND_AGE (Marie 2026): conserved biomass-weighted mean stand age (real).
    ! Absent in restart -> fallback to age_stand, then derived from the target rotation
    ! in stomate_initialize.
    var_name = 'age_stand_bm'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., age_stand_bm, 'gather', nbp_glo, index_g)
    IF ( ALL(age_stand_bm(:,:) == val_exp) ) THEN
       age_stand_bm(:,:) = REAL(age_stand(:,:), r_std)
       stand_age_bm_in_restart = .FALSE.   ! absent -> file-init (mode file) may set it
    ELSE
       stand_age_bm_in_restart = .TRUE.    ! present -> do NOT override (keep the aged value)
    ENDIF
    !-
    ! PROGRESSIVE_HARVEST (Marie 2026): the fractional clearcut is marked by
    ! sapiens_forestry in year Y and consumed by age_class_distr in year Y+1 -- the same
    ! one-year lag the whole-slot clearcut already has, age_class_distr running before
    ! sapiens_forestry inside the annual step. libIGCM runs one-year periods, so that
    ! hand-off crosses a restart boundary: without persisting ph_cut_frac the trees would
    ! be harvested and accounted while the matching area never moved down to age class 1.
    ! Allocated here (rather than only in sapiens_forestry_main) because readrestart runs
    ! first. Absent from an older restart -> zero, i.e. no pending split.
    IF (ok_progressive_harvest) THEN
       IF (.NOT. ALLOCATED(ph_cut_frac)) ALLOCATE(ph_cut_frac(npts,nvm))
       ph_cut_frac(:,:) = val_exp
       CALL restget_p (rest_id_stomate, 'ph_cut_frac', nbp_glo, nvm, 1, itime, &
            &              .TRUE., ph_cut_frac, 'gather', nbp_glo, index_g)
       IF ( ALL(ph_cut_frac(:,:) == val_exp) ) ph_cut_frac(:,:) = zero
    ENDIF
    !-
    temp_real(:,:) = val_exp
    var_name = 'rotation_n'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., temp_real, 'gather', nbp_glo, index_g)
    IF ( ALL(temp_real(:,:) == val_exp) ) THEN
       rotation_n(:,:) = 1
    ELSE
       rotation_n(:,:) = NINT(temp_real(:,:))
    ENDIF
    !-
    temp_real(:,:) = val_exp
    var_name = 'last_cut'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., temp_real, 'gather', nbp_glo, index_g)
    IF ( ALL(temp_real(:,:) == val_exp) ) THEN
       last_cut(:,:) = 0
    ELSE
       last_cut(:,:) = NINT(temp_real(:,:))
    ENDIF
    !-
    mai(:,:) = val_exp
    var_name = 'mai'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., mai, 'gather', nbp_glo, index_g)
    IF ( ALL(mai(:,:) == val_exp) ) mai(:,:) = zero
    !-
    pai(:,:) = val_exp
    var_name = 'pai'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., pai, 'gather', nbp_glo, index_g)
    IF ( ALL(pai(:,:) == val_exp) ) pai(:,:) = zero
    !-
    previous_wood_volume(:,:) = val_exp
    var_name = 'previous_wood_volume'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., previous_wood_volume, 'gather', nbp_glo, index_g)
    IF ( ALL(previous_wood_volume(:,:) == val_exp) ) previous_wood_volume(:,:) = zero
    !-
    temp_real(:,:) = val_exp
    var_name = 'mai_count'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., temp_real, 'gather', nbp_glo, index_g)
    IF ( ALL(temp_real(:,:) == val_exp) ) THEN
       mai_count(:,:) = 0
    ELSE
       mai_count = NINT(temp_real)
    ENDIF
    !-
    coppice_dens(:,:) = val_exp
    var_name = 'coppice_dens'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., coppice_dens, 'gather', nbp_glo, index_g)
    IF ( ALL(coppice_dens(:,:) == val_exp) ) coppice_dens(:,:) = zero

    
    ! Read and initialize forest_managed
    var_name = 'forest_managed'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., temp_real, 'gather', nbp_glo, index_g)

    IF ( ALL(temp_real(:,:) == val_exp) .OR. vegetmap_reset) THEN
       ! forest_managed is not found in the restart file or option to re-initialize is set.

       ! Initalize forest_managed by reading map or run.def
       CALL sapiens_forestry_set_fm(npts, lalo, neighbours, resolution, contfrac, forest_managed)
    ELSE
       ! The variable was found in restart file, transform to integer
       forest_managed = NINT(temp_real)
    ENDIF

    ! 2026-06-14 (AED): global = continuous-cover forestry. sapiens_forestry_set_fm
    ! reinterprets ifm_thin->ifm_uneven (no rotation clearfell) ONLY when (re)reading the
    ! FM map; on a restart, forest_managed is read as-is (e.g. ifm_thin from a europe
    ! spinup), bypassing it. Re-apply the reinterpretation here so the global regime is
    ! robust to restart (mirrors sapiens_forestry_set_fm ~l.2702).
    IF (concept_scale=='global') THEN
       WHERE(forest_managed==ifm_thin)
          forest_managed = ifm_uneven
       END WHERE
    END IF


    ! Read and initialize species_change_map
    IF (ok_change_species) THEN
       var_name = 'species_change_map'
       CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm , 1, itime, &
            .TRUE., r_replant(:,:), 'gather', nbp_glo, index_g)

       IF ( ALL(r_replant(:,:) == val_exp) ) THEN
          ! r_replant was not found in the restart file

          ! Initialize species_change_map by reading map or run.def
          CALL sapiens_forestry_set_species_change(npts, lalo, neighbours, &
               resolution, contfrac, species_change_map)
       ELSE
          ! r_replant was found in the restart file. Transform to INTEGER values.
          DO ipts=1,npts
             DO ivm=1,nvm
                IF (r_replant(ipts,ivm) .LT. nvm*2) THEN
                   species_change_map(ipts,ivm)=NINT(r_replant(ipts,ivm))
                ELSE
                   species_change_map(ipts,ivm)=0
                ENDIF
             ENDDO
          ENDDO
       END IF
       
    END IF ! end ok_change_species

    
    ! Read and initialize fm_change_map
    var_name = 'fm_change_map'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm , 1, itime, &
         &        .TRUE., r_replant(:,:), 'gather', nbp_glo, index_g)

    IF ( ALL(r_replant(:,:) == val_exp) ) THEN
       ! The variable was not found in the restart file. Now do initialization.

       ! Initialize fm_change_map by reading file or run.def
       CALL sapiens_forestry_set_desired_fm(npts, lalo, neighbours, &
            resolution, contfrac, forest_managed, fm_change_map)

    ELSE
       ! r_replant was found in the restart file. Now transform to integer.
       DO ipts=1,npts
          DO ivm=1,nvm
             IF (r_replant(ipts,ivm) .LT. nvm*2) THEN
                fm_change_map(ipts,ivm)=NINT(r_replant(ipts,ivm))
             ELSE
                fm_change_map(ipts,ivm)=0
             ENDIF
          ENDDO
       ENDDO
    END IF ! end ALL(r_replant==val_exp)
    
    
    IF(ok_windthrow) THEN
       CALL restget_p (rest_id_stomate, 'is_storm', nbp_glo, 1, 1, itime, &
            &              .TRUE., is_storm_real, 'gather', nbp_glo, index_g)
       IF (ALL(is_storm_real(:) == val_exp)) is_storm_real(:) = zero
       WHERE (is_storm_real(:) >= .5)
          is_storm = .TRUE.
       ELSEWHERE
          is_storm = .FALSE.
       ENDWHERE

       CALL restget_p (rest_id_stomate, 'count_storm', nbp_glo, 1, 1, itime, &
                        .TRUE., count_storm_real,'gather', nbp_glo, index_g)
       count_storm = NINT(count_storm_real)
       !-
       CALL restget_p (rest_id_stomate, 'max_wind_speed_storm', nbp_glo,  1, 1, itime, &
            &           .TRUE., max_wind_speed_storm, 'gather', nbp_glo, index_g)
       IF (ALL(max_wind_speed_storm(:) == val_exp)) max_wind_speed_storm(:) = zero
       !-
       ! JJ 2026: companion of max_wind_speed_storm
       CALL restget_p (rest_id_stomate, 'max_wind_ratio_storm', nbp_glo,  1, 1, itime, &
            &           .TRUE., max_wind_ratio_storm, 'gather', nbp_glo, index_g)
       IF (ALL(max_wind_ratio_storm(:) == val_exp)) max_wind_ratio_storm(:) = un
       !-
       CALL restget_p (rest_id_stomate, 'wind_speed_max_save', nbp_glo,  wind_days, 1, itime, &
            &           .TRUE., wind_speed_max_save, 'gather', nbp_glo, index_g)
       IF (ALL(wind_speed_max_save(:,:) == val_exp)) wind_speed_max_save(:,:) = zero
       !-
       gap_area_save(:,:,:) = val_exp
       var_name = 'gap_area_save'
       CALL restget_p &
            &    (rest_id_stomate, var_name, nbp_glo, nvm, legacy_years_wind, itime, &
            &     .TRUE., gap_area_save(:,:,:), 'gather', nbp_glo,index_g)
       IF (ALL(gap_area_save(:,:,:) == val_exp)) gap_area_save(:,:,:)= zero

       total_ba_init(:,:) = val_exp
       var_name = 'total_ba_init'
       CALL restget_p &
            &    (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
            &     .TRUE., total_ba_init(:,:), 'gather', nbp_glo,index_g)
       IF (ALL(total_ba_init(:,:) == val_exp)) total_ba_init(:,:)= zero
    ENDIF

    CALL restget_p (rest_id_stomate, 'wind_sum', nbp_glo,  1, 1, itime, &
         &           .TRUE., wind_sum, 'gather', nbp_glo, index_g)
    IF (ALL(wind_sum(:) == val_exp)) wind_sum(:) = zero
    !-
    CALL restget_p (rest_id_stomate, 'wind_mean_save', nbp_glo, legacy_years_wind, 1, itime, &
         &           .TRUE., wind_mean_save, 'gather', nbp_glo, index_g)
    IF (ALL(wind_mean_save(:,:) == val_exp)) wind_mean_save(:,:) = un
    !-
    CALL restget_p (rest_id_stomate, 'wind_ratio_sum_save', nbp_glo, wind_days, 1, itime, &
         &           .TRUE., wind_ratio_sum_save, 'gather', nbp_glo, index_g)
    IF (ALL(wind_ratio_sum_save(:,:) == val_exp)) wind_ratio_sum_save(:,:) = zero
    !-
    ! JJ 2026: companion of wind_ratio_sum_save
    CALL restget_p (rest_id_stomate, 'wind_ratio_max_save', nbp_glo, wind_days, 1, itime, &
         &           .TRUE., wind_ratio_max_save, 'gather', nbp_glo, index_g)
    IF (ALL(wind_ratio_max_save(:,:) == val_exp)) wind_ratio_max_save(:,:) = un

 

      ! Read legacy variables for pest_damage module
      wood_leftover_legacy(:,:,:) = val_exp
      DO l=1,legacy_years_wood
        WRITE(pyear_str,'(I2.2)') l
        var_name = 'wood_leftover_'//pyear_str
        CALL restget_p &
            &    (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
            &     .TRUE., wood_leftover_legacy(:,:,l), 'gather', nbp_glo,index_g)
        IF (ALL(wood_leftover_legacy(:,:,l) == val_exp)) wood_leftover_legacy(:,:,l)= zero
      ENDDO

      season_drought_legacy(:,:,:) = val_exp
      DO l=1,legacy_years
        WRITE(pyear_str,'(I2.2)') l
        var_name = 'season_drought_'//pyear_str
        CALL restget_p &
            &    (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
            &     .TRUE., season_drought_legacy(:,:,l), 'gather',nbp_glo,index_g)
        IF (ALL(season_drought_legacy(:,:,l) == val_exp)) THEN
           DO ivma = 1,nvmap
              j = start_index(ivma)
              season_drought_legacy(:,j,l)= un
              season_drought_legacy(:,j+1:,l)= zero
           END DO
        END IF
      ENDDO

      i_beetles_generation(:,:,:) = val_exp
      DO l=1,legacy_years
        WRITE(pyear_str,'(I2.2)') l
        var_name = 'i_beetles_generation_'//pyear_str
        CALL restget_p &
            &    (rest_id_stomate, var_name, nbp_glo, nvm,1, itime, &
            &     .TRUE., i_beetles_generation(:,:,l),'gather',nbp_glo,index_g)
        IF (ALL(i_beetles_generation(:,:,l) == val_exp)) i_beetles_generation(:,:,l)= un
      ENDDO

      sumTeff(:,:) = val_exp
      CALL restget_p &
            &    (rest_id_stomate, 'sumTeff', nbp_glo, nvm,1, itime, &
            &     .TRUE.,sumTeff(:,:),'gather',nbp_glo,index_g)
      IF (ALL(sumTeff(:,:) == val_exp)) sumTeff(:,:)= zero

      woody_litter_to_use(:,:) = val_exp
      CALL restget_p &
            &    (rest_id_stomate, 'woody_litter_to_use', nbp_glo, nvm,1, itime, &
            &     .TRUE.,woody_litter_to_use(:,:),'gather',nbp_glo,index_g)
      IF (ALL(woody_litter_to_use(:,:) == val_exp)) woody_litter_to_use(:,:)= zero

      beetle_diapause(:,:) = val_exp
      CALL restget_p &
            &    (rest_id_stomate, 'beetle_diapause', nbp_glo, nvm,1, itime, &
            &     .TRUE., beetle_diapause(:,:),'gather',nbp_glo,index_g)
      IF (ALL(beetle_diapause(:,:) == val_exp)) beetle_diapause(:,:)= zero

      i_beetles_activity_legacy(:,:,:) = val_exp
      DO l=1,legacy_years
        WRITE(pyear_str,'(I2.2)') l
        var_name = 'i_beetles_activity_'//pyear_str
        CALL restget_p &
           &    (rest_id_stomate, var_name, nbp_glo, nvm,1, itime, &
           &     .TRUE., i_beetles_activity_legacy(:,:,l),'gather',nbp_glo,index_g)
        IF (ALL(i_beetles_activity_legacy(:,:,l) == val_exp)) i_beetles_activity_legacy(:,:,l)= zero
      ENDDO

      B_beetles_kill_legacy(:,:,:) = val_exp
      DO l=1,beetle_legacy
        WRITE(pyear_str,'(I2.2)') l
        var_name = 'B_beetles_kill_legacy_'//pyear_str
        CALL restget_p &
            &    (rest_id_stomate, var_name, nbp_glo, nvm,1, itime, &
            &     .TRUE.,B_beetles_kill_legacy(:,:,l),'gather',nbp_glo,index_g)
        IF (ALL(B_beetles_kill_legacy(:,:,l) == val_exp)) B_beetles_kill_legacy(:,:,l)=zero
      ENDDO

      P_beetles_attacked_legacy(:,:) = val_exp
      var_name = 'P_beetles_attacked_legacy'
      CALL restget_p &
            &    (rest_id_stomate, var_name, nbp_glo, nvm,1, itime, &
            &     .TRUE.,P_beetles_attacked_legacy(:,:),'gather',nbp_glo,index_g)
        IF (ALL(P_beetles_attacked_legacy(:,:) == val_exp)) P_beetles_attacked_legacy(:,:)= zero

    ! The routines don't like to store logical variables.  So instead
    ! we create a real array.  If the value of the real array is one,
    ! we assign a value of TRUE to the logical array.  Anything
    ! else is .FALSE.
    r_replant(:,:) = val_exp
    var_name = 'lpft_replant'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm , 1, itime, &
         &        .TRUE., r_replant(:,:), 'gather', nbp_glo, index_g)
    DO ipts=1,npts
       DO ivm=1,nvm
          IF(ABS(r_replant(ipts,ivm) - un) .LT. 0.001)THEN
             lpft_replant(ipts,ivm)=.TRUE.
          ELSE
             lpft_replant(ipts,ivm)=.FALSE.
          ENDIF
       ENDDO
    ENDDO
     
    ! Read assim_param from restart file. The initialization of assim_param will 
    ! be done in stomate_var_init if the variable is not in the restart file.
    assim_param(:,:,:)  = val_exp
    CALL restget_p &
         &    (rest_id_stomate, 'assim_param', nbp_glo, nvm, npco2, itime, &
         &     .TRUE., assim_param, 'gather', nbp_glo, index_g)

    light_tran_to_floor_season(:,:) = val_exp 
    var_name = 'light_tran_season' 
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, & 
         &     .TRUE., light_tran_to_floor_season, 'gather', nbp_glo, index_g) 
    IF (ALL(light_tran_to_floor_season(:,:) == val_exp)) light_tran_to_floor_season = zero

    CALL restget_p (rest_id_stomate, 'daylight_count', nbp_glo, nvm, 1,itime, &
         .TRUE., daylight_count(:,:), 'gather', nbp_glo, index_g)
    IF ( ALL(daylight_count(:,:) == val_exp) ) daylight_count(:,:) = 0.

    ! Drought mortality
    IF (ok_hydrol_arch) THEN
       DO l = 1,nelements
          var_name = 'bio_ini_drought'//TRIM(element_str(l))
          CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm, ncirc, nparts, itime, &
               &    .TRUE., biomass_init_drought(:,:,:,:,l), 'gather', nbp_glo, index_g)
          IF ( ALL(biomass_init_drought(:,:,:,:,l) == val_exp) ) biomass_init_drought(:,:,:,:,l) = 0.
       ENDDO
       
       CALL restget_p (rest_id_stomate, 'kill_vessels', nbp_glo, nvm, ncirc , itime, &
            &   .TRUE., temp_real2, 'gather', nbp_glo, index_g)
       IF (ALL(temp_real2(:,:,:) == val_exp)) temp_real2(:,:,:) = zero
       WHERE (temp_real2(:,:,:) >= .5)
          kill_vessels = .TRUE.
       ELSEWHERE
          kill_vessels = .FALSE.
       ENDWHERE
       
       CALL restget_p (rest_id_stomate, 'vessel_loss_previous', nbp_glo, nvm, ncirc , itime, &
            &  .TRUE., vessel_loss_previous, 'gather', nbp_glo, index_g)
       IF (ALL(vessel_loss_previous(:,:,:) == val_exp)) vessel_loss_previous(:,:,:) = zero
       
    END IF
       
    ! Phenological variables
    CALL restget_p (rest_id_stomate, 'grow_season_len', nbp_glo, nvm  , 1, itime, &
         &  .TRUE., grow_season_len, 'gather', nbp_glo, index_g)
    IF (ALL(grow_season_len(:,:) == val_exp)) grow_season_len(:,:) = zero

    CALL restget_p (rest_id_stomate, 'doy_start_gs', nbp_glo, nvm  , 1, itime, &
         &  .TRUE., doy_start_gs, 'gather', nbp_glo, index_g)
    IF (ALL(doy_start_gs(:,:) == val_exp)) doy_start_gs(:,:) = zero

    CALL restget_p (rest_id_stomate, 'doy_end_gs', nbp_glo, nvm  , 1, itime, &
         &  .TRUE., doy_end_gs, 'gather', nbp_glo, index_g)
    IF (ALL(doy_end_gs(:,:) == val_exp)) doy_end_gs(:,:) = zero

    CALL restget_p (rest_id_stomate, 'mean_start_gs', nbp_glo, nvm  , 1, itime, &
         &  .TRUE., mean_start_gs, 'gather', nbp_glo, index_g)
    IF (ALL(mean_start_gs(:,:) == val_exp)) mean_start_gs(:,:) = zero

    !! peatland
    IF ( ok_peat_NoDiscretisation) THEN
       height_acro(:) = val_exp
       var_name = 'height_acro'
       CALL restget_p (rest_id_stomate,var_name, nbp_glo, 1, 1, itime, &
            &               .TRUE., height_acro, 'gather', nbp_glo,index_g)
       IF (ALL(height_acro(:) == val_exp)) THEN
          height_acro(:) = zero
       ENDIF
       
       carbon_acro(:,:) = val_exp
       var_name = 'carbon_acro'
       CALL restget_p (rest_id_stomate,var_name, nbp_glo, nvm, 1, itime, &
            &               .TRUE., carbon_acro, 'gather', nbp_glo,index_g)
       IF (ALL(carbon_acro(:,:) == val_exp))THEN
          carbon_acro(:,:) = zero
       ENDIF
       
       carbon_cato(:,:) = val_exp
       var_name = 'carbon_cato'
       CALL restget_p (rest_id_stomate,var_name, nbp_glo, nvm, 1, itime, &
            &               .TRUE., carbon_cato, 'gather', nbp_glo,index_g)
       IF (ALL(carbon_cato(:,:) == val_exp)) THEN
          carbon_cato(:,:) = zero
       ENDIF
       
       deepSOM_peat(:,:,:,:) = val_exp
       var_name = 'deepSOM_peat'
       CALL restget_p (rest_id_stomate, var_name, nbp_glo, ngrnd, nvm,nelements, itime, &
            &               .TRUE., deepSOM_peat, 'gather', nbp_glo, index_g)
       IF (ALL(deepSOM_peat == val_exp)) deepSOM_peat(:,:,:,:) = zero 
    ENDIF
       
    ! Dynamic parameters. They are global so they are not passed through
    ! the argument list of the SUBROUTINE.
    CALL restget_p (rest_id_stomate, 'pre_indust_ref_gpp', nbp_glo, nvm  , 1, itime, &
         &  .TRUE., pre_indust_ref_gpp, 'gather', nbp_glo, index_g)
       
    IF (ALL(pre_indust_ref_gpp(:,:) == val_exp)) THEN

       ! First assign the default values. The dimensions don't match. Put it
       ! into a loop.
       DO ivm = 1,nvm
          pre_indust_ref_gpp(:,ivm) = init_pre_indust_ref_gpp(ivm)
       END DO
  
    END IF
    !-
    CALL restget_p (rest_id_stomate, 'gpp_year', nbp_glo, nvm  , 1, itime, &
         &              .TRUE., gpp_year, 'gather', nbp_glo, index_g)

    IF (ALL(gpp_year(:,:) == val_exp)) THEN
       
       ! If gpp_year would be set to zero it would take long to compensate
       ! for the initial 364 zeros. Therefore gpp_year should be initialized
       ! as good as possible. In case the model starts from scratch the
       ! gpp_year will start with a fixed prescribed value of GPP.
       gpp_year(:,:) = pre_indust_ref_gpp(:,:)
       
    END IF

    CALL restget_p (rest_id_stomate, 'gpp_decade', nbp_glo, nvm  , 1, itime, &
         &              .TRUE., gpp_decade, 'gather', nbp_glo, index_g)

    IF (ALL(gpp_decade(:,:) == val_exp)) THEN
       
       ! If gpp_decade would be set to zero it would take very long to
       ! compensate for the initial 3649 zeros. Therefore gpp_decade should
       ! be initialized as good as possible.
       ! 1. In case the model starts from scratch the gpp decade will start
       ! with a fixed prescribed value of GPP.
       ! 2. During the spinup pre_indust_ref_gpp and gpp_decade are equal
       ! 3. After the spinup pre_indust_ref_gpp is no longer updated but
       ! gpp_decade is continuously updated.
       gpp_decade(:,:) = pre_indust_ref_gpp(:,:)

    END IF

    !-
    ! AED_FEEDBACK (Marie 2026) — restore the dynamic edge_length state from the
    ! restart. If no restart entry is found (val_exp), keep the freshly-allocated
    ! zero state : stomate.f90 rafraichira alors edge_length_dyn depuis
    ! Edge_length.nc a la premiere lecture.
    ! (Le test de convergence AED_SPINUP et ses scalaires ont ete retires le
    !  2026-07-26 ; les vieux restarts qui les contiennent restent lisibles.)
    !-
    IF (ok_aed_feedback) THEN
       edge_length_dyn(:) = val_exp
       var_name = 'edge_length_dyn'
       CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
            &              .TRUE., edge_length_dyn, 'gather', nbp_glo, index_g)
       IF (ALL(edge_length_dyn(:) == val_exp)) edge_length_dyn(:) = zero

    ENDIF

    !-
    ! ROTATION_GROWTH (2026-08-17). dia_growth porte une memoire de plusieurs DECENNIES et
    ! circ_dia_last le diametre de chaque classe de circonference au passage precedent :
    ! libIGCM relance le binaire a
    ! chaque periode, donc sans ecriture en reprise l'integrateur serait remis a zero
    ! CHAQUE ANNEE et n'accumulerait jamais rien. Repli sur zero = "pas d'historique",
    ! ce que update_rotation_growth traite comme non estimable.
    !-
    IF (ALLOCATED(dia_growth)) THEN
       dia_growth(:,:) = val_exp
       var_name = 'dia_growth'
       CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
            &              .TRUE., dia_growth, 'gather', nbp_glo, index_g)
       IF (ALL(dia_growth(:,:) == val_exp)) dia_growth(:,:) = zero

       circ_dia_last(:,:,:) = val_exp
       CALL restget_p (rest_id_stomate, 'circ_dia_last', nbp_glo, nvm, ncirc, itime, &
            &              .TRUE., circ_dia_last, 'gather', nbp_glo, index_g)
       IF (ALL(circ_dia_last == val_exp)) circ_dia_last(:,:,:) = zero

       circ_n_last(:,:,:) = val_exp
       CALL restget_p (rest_id_stomate, 'circ_n_last', nbp_glo, nvm, ncirc, itime, &
            &              .TRUE., circ_n_last, 'gather', nbp_glo, index_g)
       IF (ALL(circ_n_last == val_exp)) circ_n_last(:,:,:) = zero

       ! /!\ rotation_growth est DERIVEE, mais elle doit quand meme etre persistee : elle
       ! n'est recalculee que le 31 decembre par age_class_distr, alors qu'une periode vaut
       ! un an. Sans reprise elle vaut ZERO du 1er janvier au 30 decembre, donc le melange
       ! n'agit que sur les dernieres heures de l'annee -- inerte en pratique, et la sortie
       ! annuelle n'en montre qu'une moyenne diluee (mesure : 50,2 affiche pour 165 calcule).

       rotation_growth(:,:) = val_exp
       var_name = 'rotation_growth'
       CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
            &              .TRUE., rotation_growth, 'gather', nbp_glo, index_g)
       IF (ALL(rotation_growth(:,:) == val_exp)) rotation_growth(:,:) = zero
    ENDIF

    !-
    ! AED_REGROWTH (Marie 2026) — restore the class-0 conveyor (disturbed-area
    ! cohorts + annual clock). Trecruit(PFT) (2026-06-18): the conveyor 3rd dim is
    ! n_recruit_max, which may change between runs. The restart variable was
    ! renamed 'class0_area' -> 'class0_area_v2' so a pre-extension (11-bin) restart
    ! is simply absent under the new name: restget returns val_exp and the conveyor
    ! re-inits cold (a re-spinup re-equilibrates it). See MODULE_DESIGN_AED_REGROWTH.md §9bis.
    !-
    IF (ok_edge_from_age_class .AND. ALLOCATED(class0_area)) THEN
       class0_area(:,:,:) = val_exp
       CALL restget_p (rest_id_stomate, 'class0_area_v2', nbp_glo, nagent, n_recruit_max, itime, &
            &              .TRUE., class0_area, 'gather', nbp_glo, index_g)
       IF (ALL(class0_area == val_exp)) class0_area(:,:,:) = zero

       ! MATURITY_CONVEYOR : convoyeur de regeneration. INDISPENSABLE en fichier de reprise :
       ! libIGCM redemarre a chaque annee, sans persistance l'aire en attente disparaitrait a
       ! chaque periode et la regeneration ne reviendrait jamais.
       ! Guillaume M. -- N_CLASS0 must not change between two runs: class0_veg carries it
       ! as a dimension, so restget would abort on a generic "Incompatibility for I2
       ! dimension" that names no cause. Read the stored slot count first and say it plainly.
       n_class0_file = val_exp
       CALL restget_p (rest_id_stomate, 'n_class0', itime, .TRUE., val_exp, n_class0_file)
       IF (n_class0_file /= val_exp .AND. NINT(n_class0_file) /= n_class0) THEN
          WRITE(numout,*) 'N_CLASS0 in restart file =', NINT(n_class0_file), &
               ' but the namelist asks for', n_class0
          CALL flush(numout)
          CALL ipslerr_p(3,'readrestart', &
               'Attention: class 0 dimension changed, it is not allowed!', &
               'class0_veg is stored with n_class0 slots in the restart file', &
               'keep N_CLASS0 unchanged, or start from a restart without class0_veg')
       ENDIF

       IF (n_class0 > 0 .AND. ALLOCATED(class0_veg)) THEN
          class0_veg(:,:,:) = val_exp
          CALL restget_p (rest_id_stomate, 'class0_veg', nbp_glo, nvmap, n_class0, itime, &
               &              .TRUE., class0_veg, 'gather', nbp_glo, index_g)
          ! INSTRUMENTATION 2026-08-03 : le registre repart a zero chaque annee alors que
          ! le fichier de reprise le contient bien (case 2 non nulle verifiee). Le repli
          ! ci-dessous transforme un echec de relecture en zeros SANS AUCUN MESSAGE --
          ! motif deja rencontre dans ce projet. On mesure donc ce que restget a rendu.
          WRITE(numout,*) '[CLASSE0-RESTART] val_exp restants =', &
               COUNT(class0_veg == val_exp), ' / ', SIZE(class0_veg), &
               ' | somme relue =', SUM(class0_veg, MASK=(class0_veg /= val_exp))
          CALL flush(numout)
          IF (ALL(class0_veg == val_exp)) class0_veg(:,:,:) = zero
       ENDIF

       class0_year_clock_loc = val_exp
       CALL restget_p (rest_id_stomate, 'class0_year_clock_s', itime,            &
            &              .TRUE., val_exp, class0_year_clock_loc)
       IF (class0_year_clock_loc /= val_exp) class0_year_clock_s = class0_year_clock_loc

       ! AED_REGROWTH v2 leaky agent-flux integrator (npts,nagent)
       IF (ALLOCATED(agent_flux)) THEN
          agent_flux(:,:) = val_exp
          CALL restget_p (rest_id_stomate, 'aed_agent_flux', nbp_glo, nagent, 1, itime, &
               &              .TRUE., agent_flux, 'gather', nbp_glo, index_g)
          IF (ALL(agent_flux == val_exp)) agent_flux(:,:) = zero
       ENDIF
    ENDIF

    IF (printlev >= 4) WRITE(numout,*) 'Leaving readrestart ', rest_id_stomate
    !-----------------------
  END SUBROUTINE readrestart


!! ================================================================================================================================
!! SUBROUTINE   : writerestart
!!
!>\BRIEF        Write all variables for stomate from restart file. 
!!
!! DESCRIPTION  : Write all variables for stomate from restart file. 
!!                
!! \n
!_ ================================================================================================================================

  SUBROUTINE writerestart &
       & (npts, index, dt_days, date_loc, veget_max, &
       &  adapted, regenerate, vegstress_day, gdd_init_date, litterhum_daily, &
       &  t2m_daily, t2m_min_daily,t2m_max_daily, tsurf_daily, tsoil_daily, &
       &  precip_daily,vpd_daily_mean,vpd_daily_max,vpd_mean_week, vpd_max_week, gpp_daily, npp_daily, &
       &  turnover_daily, turnover_resid, vegstress_month, &
       &  vegstress_week, vegstress_season, &
       &  t2m_longterm, tau_longterm, t2m_month, t2m_week, &
       &  tsoil_month,precip_month, fireindex, firelitter, &
       &  maxvegstress_lastyear, maxvegstress_thisyear, &
       &  minvegstress_lastyear, minvegstress_thisyear, &
       &  maxgppweek_lastyear, maxgppweek_thisyear, &
       &  gdd0_lastyear, gdd0_thisyear, &
       &  precip_lastyear, precip_thisyear, &
       &  gdd_m5_dormance, gdd_from_growthinit, gdd_midwinter, ncd_dormance, ngd_minus5, &
       &  PFTpresent, npp_longterm, croot_longterm, n_reserve_longterm, lm_lastyearmax, &
       &  lm_thisyearmax, maxfpc_lastyear, maxfpc_thisyear, &
       &  turnover_longterm, gpp_week, gpp_year, gpp_decade, &
       &  resp_maint_part, resp_maint_week, &
       &  leaf_age, leaf_frac, leaf_age_crit, plant_status, when_growthinit, age, &
       &  resp_hetero, resp_maint, resp_growth, &
       &  co2_fire, &
       &  veget_lastlight, everywhere, need_adjacent, RIP_time, &
       &  time_hum_min, hum_min_dormance, &
       &  litter, dead_leaves, &
       &  som, lignin_struc, lignin_wood, lignin_snag, turnover_time, &
       &  fco2_lu, fco2_wh, fco2_ha,  &
       &  prod_s, prod_m, prod_l, &
       &  flux_s, flux_m, flux_l, &
       &  fDeforestToProduct, fLulccResidue, fHarvestToProduct, &
       &  bm_to_litter, bm_to_litter_resid, tree_bm_to_litter, &
       &  tree_bm_to_litter_resid, carb_mass_total, &
       &  Tseason, Tseason_length, Tseason_tmp, & 
       &  Tmin_spring_time, &
       &  global_years, ok_equilibrium, nbp_accu_flux, &
       &  nbp_pool_start, &
       &  MatrixV, VectorU, previous_stock, current_stock, &
       &  assim_param, CN_som_litter_longterm, &
       &  tau_CN_longterm, KF, k_latosa_adapt, &
       &  rue_longterm, cn_leaf_min_season, nstress_season, soil_n_min, p_O2, bact, &
       &  forest_managed, &
       &  species_change_map, fm_change_map, lpft_replant, lai_per_level, &
       &  laieff_fit, wstress_season, wstress_month,&
       &  age_stand, age_stand_bm, rotation_n, last_cut, mai, pai, &
       &  previous_wood_volume, mai_count, coppice_dens, &
       &  light_tran_to_floor_season,daylight_count, gap_area_save, &
       &  deepSOM_a, deepSOM_s, deepSOM_p, O2_soil, CH4_soil, O2_snow, CH4_snow, &
       &  heat_Zimov, altmax, altmax_ind, altmax_lastyear, altmax_ind_lastyear, &
       &  depth_organic_soil, fixed_cryoturbation_depth, &
       &  sugar_load, harvest_cut, & 
       &  harvest_pool_acc, harvest_area_acc, burried_litter, burried_fresh_ltr, &
       &  burried_fresh_som, burried_bact, &
       &  burried_min_nitro,burried_som, &
       &  burried_deepSOM_a, burried_deepSOM_s, burried_deepSOM_p,&
       &  wood_leftover_legacy,i_beetles_activity_legacy,season_drought_legacy,&
       &  P_beetles_attacked_legacy, beetle_diapause, sumTeff, woody_litter_to_use, &
       &  i_beetles_generation, B_beetles_kill_legacy, max_wind_speed_storm, &
       &  max_wind_ratio_storm, wind_ratio_max_save, wind_ratio_sum_save, & ! JJ 2026: added companions
       &  wind_speed_max_save, wind_mean_save, wind_sum, is_storm, count_storm, &
       &  biomass_init_drought, kill_vessels, &
       &  ni_acc, litterfuel, &
       &  vessel_loss_previous, grow_season_len, doy_start_gs, doy_end_gs, &
       &  mean_start_gs, total_ba_init, & 
       &  carbon_acro, carbon_cato, height_acro,deepSOM_peat)

    ! 0 declarations
    !-
    ! 0.1 input
    !-
    INTEGER(i_std),INTENT(in)                             :: npts                     !! Domain size
    INTEGER(i_std),DIMENSION(:),INTENT(in)                :: index                    !! Indices of the points on the map
    REAL(r_std),INTENT(in)                                :: dt_days                  !! time step of STOMATE in days
    INTEGER(i_std),INTENT(in)                             :: date_loc                 !! date_loc (d)
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: veget_max                !! Maximum fraction of vegetation type including non-biological fraction (unitless)
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: adapted                  !! Winter too cold? between 0 and 1
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: regenerate               !! Winter sufficiently cold? between 0 and 1
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: vegstress_day         !! daily moisture availability
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: gdd_init_date            !! date for beginning of gdd count
    REAL(r_std),DIMENSION(:),INTENT(in)                   :: litterhum_daily          !! daily litter humidity
    REAL(r_std),DIMENSION(:),INTENT(in)                   :: t2m_daily                !! daily 2 meter temperatures (K)
    REAL(r_std),DIMENSION(:),INTENT(in)                   :: t2m_min_daily            !! daily minimum 2 meter temperatures (K)
    REAL(r_std),DIMENSION(:),INTENT(in)                   :: t2m_max_daily            !! daily maximum 2 meter temperatures (K)
    REAL(r_std),DIMENSION(:),INTENT(in)                   :: tsurf_daily              !! daily surface temperatures (K)
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: tsoil_daily              !! daily soil temperatures (K)
    REAL(r_std),DIMENSION(:),INTENT(in)                   :: precip_daily             !! daily precipitations (mm/day) (for phenology)
    REAL(r_std),DIMENSION(:),INTENT(in)                   :: vpd_daily_mean           !! daily mean vpd (kPa/day) (for spitfire)
    REAL(r_std),DIMENSION(:),INTENT(in)                   :: vpd_daily_max            !! daily max vpd (kPa/day) (for spitfire)
    REAL(r_std),DIMENSION(:),INTENT(in)                   :: vpd_mean_week            !! Weekly mean vpd (kPa/day) (for spitfire)
    REAL(r_std),DIMENSION(:),INTENT(in)                   :: vpd_max_week             !! Weekly max vpd (kPa/day) (for spitfire)
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: gpp_daily                !! daily gross primary productivity (gC/m**2/day)
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: npp_daily                !! daily net primary productivity (gC/m**2/day)
    REAL(r_std),DIMENSION(:,:,:,:),INTENT(in)             :: turnover_daily           !! daily turnover rates (gC/m**2/day)
    REAL(r_std),DIMENSION(:,:,:,:),INTENT(in)             :: turnover_resid           !! The turnover left from turnover_daily at any given time step  
                                                                                      !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex 
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: vegstress_month       !! "monthly" moisture availability
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: vegstress_week        !! "weekly" moisture availability
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: vegstress_season      !! mean growing season moisture availability (used for allocation response)
    REAL(r_std),DIMENSION(:),INTENT(in)                   :: t2m_longterm             !! "long term" 2 meter temperatures (K)
    REAL(r_std), INTENT(in)                               :: tau_longterm             !! "tau_longterm"
    REAL(r_std),DIMENSION(:),INTENT(in)                   :: t2m_month                !! "monthly" 2 meter temperatures (K)
    REAL(r_std),DIMENSION(:),INTENT(in)                   :: Tseason                  !! "seasonal" 2 meter temperatures (K) 
    REAL(r_std),DIMENSION(:),INTENT(in)                   :: Tseason_length           !! temporary variable to calculate Tseason
    REAL(r_std),DIMENSION(:),INTENT(in)                   :: Tseason_tmp              !! temporary variable to calculate Tseason
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: Tmin_spring_time         !!
    REAL(r_std),DIMENSION(:),INTENT(in)                   :: t2m_week                 !! "weekly" 2 meter temperatures (K)
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: tsoil_month              !! "monthly" soil temperatures (K)
    REAL(r_std),DIMENSION(:),INTENT(in)                   :: precip_month             !! "monthly" precipitation
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: fireindex                !! Probability of fire
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: firelitter               !! Longer term total litter above the ground, gC/m**2 of ground
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: maxvegstress_lastyear !! last year's maximum moisture availability
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: maxvegstress_thisyear !! this year's maximum moisture availability
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: minvegstress_lastyear !! last year's minimum moisture availability
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: minvegstress_thisyear !! this year's minimum moisture availability
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: maxgppweek_lastyear      !! last year's maximum weekly GPP
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: maxgppweek_thisyear      !! this year's maximum weekly GPP
    REAL(r_std),DIMENSION(:),INTENT(in)                   :: gdd0_lastyear            !! last year's annual GDD0
    REAL(r_std),DIMENSION(:),INTENT(in)                   :: gdd0_thisyear            !! this year's annual GDD0
    REAL(r_std),DIMENSION(:),INTENT(in)                   :: precip_lastyear          !! last year's annual precipitation (mm/year)
    REAL(r_std),DIMENSION(:),INTENT(in)                   :: precip_thisyear          !! this year's annual precipitation (mm/year)
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: gdd_m5_dormance          !! growing degree days, threshold -5 deg C (for phenology)
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: gdd_from_growthinit      !! growing degree days, from begin of season
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: gdd_midwinter            !! growing degree days since midwinter (for phenology)
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: ncd_dormance             !! number of chilling days since leaves were lost (for phenology)
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: ngd_minus5               !! number of growing days, threshold -5 deg C (for phenology)
    LOGICAL,DIMENSION(:,:),INTENT(in)                     :: PFTpresent               !! PFT exists (equivalent to fpc_max > 0 for natural PFTs)
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: npp_longterm             !! "long term" net primary productivity (gC/m**2/year)
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: croot_longterm           !! "long term" root carbon mass (gC/m**2/year)
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: n_reserve_longterm       !! "long term" actual to potential N reserve pool (unitless)
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: lm_lastyearmax           !! last year's maximum leaf mass, for each PFT (gC/m**2)
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: lm_thisyearmax           !! this year's maximum leaf mass, for each PFT (gC/m**2)
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: maxfpc_lastyear          !! last year's maximum fpc for each natural PFT, on ground
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: maxfpc_thisyear          !! this year's maximum fpc for each PFT, on *total* ground (see stomate_season)   
    REAL(r_std),DIMENSION(:,:,:,:),INTENT(in)             :: turnover_longterm        !! "long term" turnover rate (gC/m**2/year)
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: gpp_week                 !! "weekly" GPP (gC/day/(m**2 covered)
    REAL(r_std),DIMENSION(:,:),INTENT(inout)              :: gpp_year                 !! "annual" GPP (gC/day/(m**2 covered)
    REAL(r_std),DIMENSION(:,:),INTENT(inout)              :: gpp_decade               !! "decadal" GPP (gC/day/(m**2 covered)
    REAL(r_std),DIMENSION(:,:,:),INTENT(in)               :: resp_maint_part          !! maintenance resp (gC/m**2)
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: resp_maint_week          !! "weekly" maintenance respiration (gC/day/(m**2 covered)
    REAL(r_std),DIMENSION(:,:,:),INTENT(in)               :: leaf_age                 !! leaf age (days)
    REAL(r_std),DIMENSION(:,:,:),INTENT(in)               :: leaf_frac                !! fraction of leaves in leaf age class
    REAL(r_std), DIMENSION(:,:),INTENT(in)                :: leaf_age_crit            !! critical leaf age (days)
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: plant_status             !! Growth and phenological status of the plant
                                                                                       !! The different stati are defined in constantes
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: when_growthinit          !! how many days ago was the beginning of the growing season
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: age                      !! mean age (years)
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: resp_hetero              !! heterotrophic respiration (gC/day/m**2)
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: resp_maint               !! maintenance respiration (gC/day/m**2)
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: resp_growth              !! growth respiration (gC/day/m**2)
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: co2_fire                 !! carbon emitted into the atmosphere by fire (living and dead biomass)
                                                                                       !! (in gC/m**2/time step)
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: veget_lastlight          !! vegetation fractions (on ground) after last light competition
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: everywhere               !! is the PFT everywhere in the grid box or very localized (after its introduction)
    LOGICAL,DIMENSION(:,:),INTENT(in)                     :: need_adjacent            !! in order for this PFT to be introduced, does it have to be present in an adjacent grid box?
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: RIP_time                 !! How much time ago was the PFT eliminated for the last time (y)
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: time_hum_min             !! time elapsed since strongest moisture availability (d)
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: hum_min_dormance         !! minimum moisture during dormance
    REAL(r_std),DIMENSION(:,:,:,:,:),INTENT(in)           :: litter                   !! fraction of litter above the ground belonging to different PFTs
                                                                                       !! separated for natural and agricultural PFTs.
    REAL(r_std),DIMENSION(:,:,:),INTENT(in)               :: dead_leaves              !! dead leaves on ground, per PFT, metabolic and structural in gC/(m**2 of ground)
    REAL(r_std),DIMENSION(:,:,:,:),INTENT(in)             :: som                      !! Soil Organic Matter pool: active, slow, or passive, (gC (or N)/m**2)
    REAL(r_std),DIMENSION(:,:,:),INTENT(in)               :: lignin_struc             !! ratio Lignine/Carbon in structural litter, above and below ground,(gC/m**2)
    REAL(r_std),DIMENSION(:,:,:),INTENT(in)               :: lignin_wood              !! ratio Lignine/Carbon in woody litter, above and below ground,(gC/m**2)
    REAL(r_std),DIMENSION(:,:,:),INTENT(in)               :: lignin_snag              !! ratio Lignine/Carbon in snag litter, above and below ground,(gC/m**2)
    REAL(r_std),DIMENSION(:,:,:),INTENT(in)               :: turnover_time            !!
    INTEGER(i_std), INTENT(in)                            :: global_years             !! for spinup matrix  
    LOGICAL, DIMENSION(:), INTENT(in)                     :: ok_equilibrium           !! for spinup matrix  
    REAL(r_std), DIMENSION(:,:), INTENT(in)               :: nbp_accu_flux            !! accumulated Net Biospheric Production over the whole simulationm (gC/N m-2)
    REAL(r_std), DIMENSION(:,:), INTENT(in)               :: nbp_pool_start           !! C an dN stocks at previous time step (gC/N m-2)
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)           :: MatrixV                  !!
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)             :: VectorU                  !!
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)             :: previous_stock           !!
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)             :: current_stock            !!
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)             :: CN_som_litter_longterm   !! Longterm CN ratio of litter and som pools (gC/gN)
    REAL(r_std), INTENT(in)                               :: tau_CN_longterm          !! Counter used for calculating the longterm CN ratio of SOM and litter pools (seconds)
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)             :: assim_param              !!
    REAL(r_std), DIMENSION(:,:), INTENT(in)               :: KF                       !! Scaling factor to convert sapwood mass into leaf mass (m)
    REAL(r_std), DIMENSION(:,:), INTENT(in)               :: k_latosa_adapt           !! Leaf to sapwood area adapted for water stress. Adaptation takes place at the end of the year (m)
    REAL(r_std), DIMENSION(:,:), INTENT(in)               :: mai                      !! The mean annual increment @tex $(m**3 / m**2 / year)$ @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(in)               :: pai                      !! The period annual increment @tex $(m**3 / m**2 / year)$ @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(in)               :: previous_wood_volume     !! The volume of the tree trunks in a stand for the previous year. @tex $(m**3 / m**2 )$ @endtex
    INTEGER(i_std), DIMENSION(:,:),INTENT(in)             :: mai_count                !! The number of times we've calculated the volume increment for a stand
    REAL(r_std), DIMENSION(:,:),INTENT(in)                :: coppice_dens             !! The density of a coppice at the first cutting. @tex $( 1 / m**2 )$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(in)               :: rue_longterm             !! longterm radiation use efficiency
    INTEGER(i_std), DIMENSION(:,:), INTENT(in)            :: age_stand                !! Age of stand (years)
    REAL(r_std), DIMENSION(:,:), INTENT(in)               :: age_stand_bm             !! Biomass-weighted conserved mean stand age (years) - STAND_AGE
    INTEGER(i_std), DIMENSION(:,:), INTENT(in)            :: rotation_n               !! Rotation number (number of rotation since pft is managed)
    INTEGER(i_std), DIMENSION(:,:), INTENT(in)            :: last_cut                 !! Years since last thinning (years)
    REAL(r_std), DIMENSION(:,:), INTENT(in)               :: cn_leaf_min_season       !! Seasonal min CN ratio of leaves 
    REAL(r_std), DIMENSION(:,:), INTENT(in)               :: nstress_season           !! N-related seasonal stress (used for allocation) 
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)             :: soil_n_min               !! mineral nitrogen in the soil (gN/m**2) (first index=npts, second index=nvm, third index=nnspec) 
    REAL(r_std), DIMENSION(:,:), INTENT(in)               :: p_O2                     !! partial pressure of oxigen in the soil (hPa)(first index=npts, second index=nvm)
    REAL(r_std), DIMENSION(:,:), INTENT(in)               :: bact                     !! denitrifier biomass (gC/m**2) (first index=npts, second index=nvm)
    INTEGER(i_std), DIMENSION(:,:), INTENT(in)            :: forest_managed           !! forest management flag
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(in)         :: prod_s                   !!
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(in)         :: prod_m                   !!
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(in)         :: prod_l                   !!
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(in)         :: flux_s                   !!
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(in)         :: flux_m                   !!
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(in)         :: flux_l                   !!
    INTEGER(i_std), DIMENSION(:,:), INTENT(in)            :: species_change_map       !! A map which gives the PFT number that each PFT will be replanted as in case of a clearcut.
                                                                                       !! (1-nvm,unitless)
    INTEGER(i_std), DIMENSION(:,:), INTENT(in)            :: fm_change_map            !! A map which gives the desired FM strategy when the PFT will be replanted after a clearcut.
                                                                                       !! (1-nvm,unitless)
    LOGICAL, DIMENSION(:,:), INTENT(in)                   :: lpft_replant             !! Indicates if this PFT has either died this year or been clearcut/coppiced.  If it has, it is not
                                                                                       !! replanted until the end of the year.
    REAL(r_std), DIMENSION(:,:,:),INTENT(in)              :: lai_per_level            !! The amount of LAI in each physical canopy level. @tex $( m**2 / m**2 )$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:),INTENT(in)            :: deepSOM_a                !!
    REAL(r_std), DIMENSION(:,:,:,:),INTENT(in)            :: deepSOM_s                !!
    REAL(r_std), DIMENSION(:,:,:,:),INTENT(in)            :: deepSOM_p                !!
    REAL(r_std), DIMENSION(:,:,:),INTENT(in)              :: O2_soil                  !!
    REAL(r_std), DIMENSION(:,:,:),INTENT(in)              :: CH4_soil                 !!
    REAL(r_std), DIMENSION(:,:,:),INTENT(in)              :: O2_snow                  !!
    REAL(r_std), DIMENSION(:,:,:),INTENT(in)              :: CH4_snow                 !!
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)             :: heat_Zimov               !! heating associated with decomposition [W/m**3 soil]
    REAL(r_std), DIMENSION(:,:), INTENT(in)               :: altmax                   !! Active layer thickness (m) 
    INTEGER(i_std), DIMENSION(npts,nvm), INTENT(in)       :: altmax_ind
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)          :: altmax_lastyear
    INTEGER(i_std), DIMENSION(npts,nvm), INTENT(in)       :: altmax_ind_lastyear
    REAL(r_std), DIMENSION(:),INTENT(in)                  :: depth_organic_soil       !! Depth at which there is still organic matter (m)
    REAL(r_std), DIMENSION(:,:), INTENT(in)               :: fixed_cryoturbation_depth!! Depth to hold cryoturbation to for fixed runs  
    REAL(r_std), DIMENSION(:,:), INTENT(in)               :: sugar_load               !! Relative sugar loading of the labile pool (unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(in)               :: harvest_cut              !! Type of cutting that was used for the harvest (unitless) 
    TYPE(laieff_type),DIMENSION (:,:,:),INTENT(in)        :: laieff_fit               !! Fitted parameters for the effective LAI
    REAL(r_std), DIMENSION(:,:), INTENT(in)               :: wstress_season           !! Water stress factor, based on hum_rel_daily (unitless, 0-1)
    REAL(r_std), DIMENSION(:,:), INTENT(in)               :: wstress_month            !! Water stress factor, based on hum_rel_daily (unitless, 0-1)
    REAL(r_std), DIMENSION(:,:), INTENT(in)               :: fDeforestToProduct       !!
    REAL(r_std), DIMENSION(:,:), INTENT(in)               :: fLulccResidue            !!
    REAL(r_std), DIMENSION(:,:), INTENT(in)               :: fHarvestToProduct        !!
    REAL(r_std),DIMENSION(:,:,:,:),INTENT(in)             :: bm_to_litter             !! Background (not senescence-driven) mortality of biomass
                                                                                      !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex
    REAL(r_std),DIMENSION(:,:,:,:),INTENT(in)             :: bm_to_litter_resid       !! Left over bm_to_litter at any specific time step
                                                                                      !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex
    REAL(r_std),DIMENSION(:,:,:,:),INTENT(in)             :: tree_bm_to_litter        !! Conversion of biomass to litter 
                                                                                      !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std),DIMENSION(:,:,:,:),INTENT(in)             :: tree_bm_to_litter_resid  !! Left over bm_to_litter_resid. Written here, used in stomate.f90 
                                                                                      !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std),DIMENSION(:),INTENT(in)                   :: carb_mass_total          !!
    REAL(r_std), DIMENSION(:,:), INTENT(in)               :: light_tran_to_floor_season !! Mean seasonal fraction of light transmitted to the forest floor (unitless, 0-1)
    REAL(r_std), DIMENSION(:,:), INTENT(in)               :: daylight_count           !! Time steps dt_radia during daylight and when there is growth (gpp>0)
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)             :: gap_area_save            !! Total gap area created by more than 30% basal area loss in the last 5 years (m^{2})
    REAL(r_std), DIMENSION(:,:), INTENT(in)               :: total_ba_init            !! Total basal area saved at the first day of the year (m^{2}/m^{2})
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)             :: season_drought_legacy    !! mean growing season moisture availability 
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)             :: wood_leftover_legacy     !!
    REAL(r_std), DIMENSION(:,:), INTENT(in)               :: P_beetles_attacked_legacy        !!
    INTEGER(i_std), DIMENSION(:,:), INTENT(in)            :: beetle_diapause          !!
    REAL(r_std), DIMENSION(:,:), INTENT(in)               :: sumTeff                  !! sum of temperazture for beetle phenology
    REAL(r_std), DIMENSION(:,:), INTENT(in)               :: woody_litter_to_use      !! Yearly sum of woody litter pool to use for wood_leftover in the pest module
                                                                                      !! (gC m^{-2})
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)             :: i_beetles_activity_legacy        !! biomass of tree from the same species that was infected during the previous timestep 
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)             :: B_beetles_kill_legacy
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)             :: i_beetles_generation  !! number of generation that BB can achieved in one year
    LOGICAL,DIMENSION(:),INTENT(in)                       :: is_storm
    INTEGER(i_std),DIMENSION(:),INTENT(in)                :: count_storm
    REAL(r_std), DIMENSION(:), INTENT(in)                 :: max_wind_speed_storm     !! Daily maximum wind speed at 2 meter (ms-1) during a storm event
    REAL(r_std), DIMENSION(:), INTENT(in)                 :: max_wind_ratio_storm     !! JJ 2026: Daily maximum wind ratio during a storm event
    REAL(r_std), DIMENSION(:,:), INTENT(in)               :: wind_ratio_sum_save      !! Sum of wind speed ratio over longterm wind speed when this ratio exceeded certain threshold,
                                                                                      !! stored for wind_days
    REAL(r_std), DIMENSION(:,:), INTENT(in)               :: wind_ratio_max_save      !! JJ 2026: Daily maximum wind ratio, stored for wind_days (default 3)
    REAL(r_std), DIMENSION(:,:), INTENT(in)               :: wind_speed_max_save      !! Daily maximum wind speed, stored for wind_days (default 3) (m s-1)
    REAL(r_std), DIMENSION(:,:), INTENT(in)               :: wind_mean_save           !! Annual average wind speed, stored for legacy years (m s-1)
    REAL(r_std), DIMENSION(:), INTENT(in)                 :: wind_sum                 !! Sum of wind speed for all time step (m s-1)
    REAL(r_std),DIMENSION(:,:,:,:,:),INTENT(in)           :: harvest_pool_acc         !! Records the quantity of wood harvested and thinned due to forest management and LCC.
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)             :: harvest_area_acc         !! Harvested area (m^{2})
    REAL(r_std),DIMENSION(:),INTENT(in)                   :: fco2_lu                  !!
    REAL(r_std),DIMENSION(:),INTENT(in)                   :: fco2_wh                  !!
    REAL(r_std),DIMENSION(:),INTENT(in)                   :: fco2_ha                  !!
    REAL(r_std),DIMENSION(:,:,:,:),INTENT(in)             :: burried_litter           !! Litter burried under non-biological land uses (gC orNm-2)
    REAL(r_std),DIMENSION(:,:,:),INTENT(in)               :: burried_fresh_ltr        !! Fresh litter burried under non-biological land uses (gC orN m-2)
    REAL(r_std),DIMENSION(:,:,:),INTENT(in)               :: burried_fresh_som        !! Fresh som burried under non-biological land uses (gC or Nm-2)
    REAL(r_std),DIMENSION(:),INTENT(in)                   :: burried_bact             !! Bacteria burried under non-biological land uses (gC m-2)
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: burried_min_nitro        !! Mineral nitrogen burried under non-biological land uses(gC or N m-2)
    REAL(r_std),DIMENSION(:,:,:),INTENT(in)               :: burried_som              !! Som burried under non-biological land uses (gC or N m-2)
    REAL(r_std),DIMENSION(:,:,:),INTENT(in)               :: burried_deepSOM_a        !! Som burried under non-biological land uses (gC or N m-2)
    REAL(r_std),DIMENSION(:,:,:),INTENT(in)               :: burried_deepSOM_s        !! Som burried under non-biological land uses (gC or N m-2)
    REAL(r_std),DIMENSION(:,:,:),INTENT(in)               :: burried_deepSOM_p        !! Som burried under non-biological land uses (gC or N m-2)
    REAL(r_std),DIMENSION(:,:,:,:,:),INTENT(in)           :: biomass_init_drought     !! Biomass of heartwood or sapwood before onset of drought. 
                                                                                       !! Used to compute turnover on same reference biomass in 
                                                                                       !! stomate_turnover.f90. Should remain the same along one 
                                                                                       !! entire drought episode and be updated inbetween 
                                                                                       !! droughts (gCor N tree-1).
    REAL(r_std),DIMENSION(:),INTENT(in)                  :: ni_acc                   !! Nesterov index (square of degree Celcius)
    REAL(r_std), DIMENSION(:,:,:,:,:),INTENT(in)         :: litterfuel               !! Dead litter fuel above ground. (gC m^{-2})
    LOGICAL,DIMENSION(:,:,:),INTENT(in)                   :: kill_vessels             !! Flag to kill vessels at the end of the day when there is embolism.
    REAL(r_std),DIMENSION(:,:,:),INTENT(in)               :: vessel_loss_previous     !! vessel loss at the previous time step
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: grow_season_len          !! growing season length in days for deciduous PFTs. 
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: doy_start_gs             !! growing season starting day of year (DOY) for deciduous PFTs.
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: doy_end_gs               !! growing season end day of year (DOY) for deciduous PFTs.
    REAL(r_std),DIMENSION(:,:),INTENT(in)                 :: mean_start_gs            !! mean growing season starting day for deciduous PFTs.
    !! peatland
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in)            :: carbon_acro              !! carbon of acrotelm
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in)            :: carbon_cato              !! carbon of catotelm
    REAL(r_std),DIMENSION(npts),INTENT(in)                :: height_acro              !! height of acrotelm
    REAL(r_std), DIMENSION(npts,ngrnd,nvm,nelements),INTENT(in) :: deepSOM_peat       !! soil carbon in peat
    
!! 0.4 Local variables
    
    REAL(r_std)                                                         :: date_real               !! date, real
    REAL(r_std),DIMENSION(npts,nvm)                                     :: PFTpresent_real         !! PFT exists (equivalent to fpc_max > 0 for natural PFTs), real
    REAL(r_std),DIMENSION(npts)                                         :: is_storm_real
    REAL(r_std),DIMENSION(npts)                                         :: count_storm_real
    REAL(r_std),DIMENSION(npts,nvm)                                     :: need_adjacent_real      !! in order for this PFT to be introduced,
                                                                                                   !! does it have to be present in an adjacent grid box? - real 
    CHARACTER(LEN=80)                                                   :: var_name                !! To store variables names for I/O
    CHARACTER(LEN=10)                                                   :: circ_str                !! string suffix indicating an index
    CHARACTER(LEN=10)                                                   :: part_str                !! string suffix indicating an index
    REAL(r_std), DIMENSION(npts)                                        :: ok_equilibrium_real     
    REAL(r_std),DIMENSION(1)                                            :: xtmp                    !! temporary storage
    INTEGER(i_std)                                                      :: j,k,l,m                 !! index
    CHARACTER(LEN=2),DIMENSION(nelements)                               :: element_str             !! string suffix indicating element
    CHARACTER(LEN=6),DIMENSION(nbpools)                                 :: pools_str
    INTEGER                                                             :: n,ilev,ipts,ivm         !! Indices
    REAL(r_std), DIMENSION(npts,nvm)                                    :: temp_real               !! temporary real to allow restget 
                                                                                                   !! to work on multi-dimensional integers
    REAL(r_std), DIMENSION(npts,nvm,ncirc)                              :: temp_real2              !! temporary real to allow restget 
                                                                                                   !! to work on multi-dimensional integers
    REAL(r_std), DIMENSION(npts,nvm)                                    :: r_replant               !! Getting logical values from the restart
                                                                                                   !! is not possible, so this is a temporary
                                                                                                   !! array where 1.0 is TRUE.
    REAL(r_std),DIMENSION(npts,nvm,nlevels_tot,nparams_laieff)          :: temp_array              !! To store structure values for I/O
    CHARACTER(LEN=10)                                                   :: part_str2               !! string suffix indicating an index
    
    CHARACTER(LEN=2), DIMENSION(legacy_years_wind+1)                    :: wyear_str               !! string suffix indicating wind year index
    CHARACTER(LEN=2)                                                    :: pyear_str               !! string suffix indicating year index for pest module
    CHARACTER(LEN=10), DIMENSION(nlctypes)                              :: lctype_str              !! string suffix for the land cover type
    CHARACTER(LEN=10), DIMENSION(nlanduse)                              :: luse_str                !! string suffix for the land use type
!!_================================================================================================================================   

    IF (printlev >= 3) WRITE(numout,*) 'Entering writerestart'

    !MatrixV(:,:,:,:) = 1.0
    !VectorU(:,:,:) = 0.0
    !previous_stock(:,:,:) = 0.0
    !current_stock(:,:,:) = 0.0
    !tau_CN_longterm = 0.0
    !total_ba_init(:,:) = 0.0
    !CN_som_litter_longterm(:,:,:)= 0.0

    !-
    ! 1 string definitions
    !-
    DO l=1,nlctypes
       IF (l == iforest) THEN
          lctype_str(l) = '_forest'
       ELSEIF (l == igrass) THEN
          lctype_str(l) = '_grass'
       ELSEIF (l == icrop) THEN
          lctype_str(l) = '_crop'
       ELSE
          CALL ipslerr_p(3,'stomate_io writerestart','Define lctype_str(l)','','')
       ENDIF   
    END DO
    !-
    DO l=1,nelements
       IF     (l == icarbon) THEN
          element_str(l) = '_c'
       ELSEIF (l == initrogen) THEN
          element_str(l) = '_n'
       ELSE
          CALL ipslerr_p(3,'stomate_io writerestart','Define element_str','','')
       ENDIF
    ENDDO
    !-
    DO l=1,nlanduse
       IF (l == iharvest) THEN
          luse_str(l) = '_harvest'
       ELSEIF (l == ilcc) THEN
          luse_str(l) = '_lcc'
       ELSE
          CALL ipslerr_p(3,'stomate_io readrestart','Define luse_str','','')
       ENDIF
    ENDDO
    !-
    pools_str(1:nbpools) =(/'str_ab ','str_be ','met_ab ','met_be ','wood_ab','wood_be', & 
         & 'snag_ab', 'snag_be', 'actif  ','slow   ','passif ','surface'/) 
    !-
    IF (is_root_prc) THEN
       CALL ioconf_setatt_p ('UNITS','-')
       CALL ioconf_setatt_p ('LONG_NAME',' ')
    ENDIF
    
    !-
    ! 2.2 time step of STOMATE in days
    !-
    CALL restput_p (rest_id_stomate, 'date', itime, date_loc)
    !-
    CALL restput_p (rest_id_stomate, 'dt_days', itime, dt_days)
    !-
    ! 2.3 date
    !-
    CALL restput_p (rest_id_stomate, 'date', itime, date_loc)
    !-
    ! 3 daily meteorological variables
    !-
    var_name = 'vegstress_day'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                vegstress_day, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'gdd_init_date'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,    2, 1, itime, &
         &              gdd_init_date, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'litterhum_daily'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,    1, 1, itime, &
         &                litterhum_daily, 'scatter', nbp_glo, index_g)
    !-
    var_name = 't2m_daily'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,    1, 1, itime, &
         &                t2m_daily, 'scatter', nbp_glo, index_g)
    !-
    var_name = 't2m_min_daily'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,    1, 1, itime, &
         &                t2m_min_daily, 'scatter', nbp_glo, index_g)
    !-
    var_name = 't2m_max_daily'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,    1, 1, itime, &
         &                    t2m_max_daily, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'tsurf_daily'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,    1, 1, itime, &
         &                tsurf_daily, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'tsoil_daily'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nslm, 1, itime, &
         &                tsoil_daily, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'precip_daily'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,    1, 1, itime, &
         &                precip_daily, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'vpd_daily_mean'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,    1, 1, itime, &
         &                vpd_daily_mean, 'scatter', nbp_glo, index_g)
    var_name = 'vpd_daily_max'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,    1, 1, itime, &
         &                vpd_daily_max, 'scatter', nbp_glo, index_g)
    var_name = 'vpd_mean_week'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,    1, 1, itime, &
         &                vpd_mean_week, 'scatter', nbp_glo, index_g)
    var_name = 'vpd_max_week'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,    1, 1, itime, &
         &                vpd_max_week, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'wstress_season'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                wstress_season, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'wstress_month'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                wstress_month, 'scatter', nbp_glo, index_g)
    !-
    ! 4 productivities
    !-
    var_name = 'gpp_daily'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                gpp_daily, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'npp_daily'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                npp_daily, 'scatter', nbp_glo, index_g)
    !-
    CALL restput_p (rest_id_stomate, 'turnover_daily', nbp_glo, nvm, nparts, nelements, itime, &
         &                   turnover_daily, 'scatter', nbp_glo, index_g)
    !-
    CALL restput_p (rest_id_stomate, 'turnover_resid', nbp_glo, nvm, nparts, nelements, itime, &
         &                   turnover_resid, 'scatter', nbp_glo, index_g)

    !-
    ! 5 monthly meteorological variables
    !-
    var_name = 'vegstress_month'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                vegstress_month, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'vegstress_week'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                vegstress_week, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'vegstress_season'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                vegstress_season, 'scatter', nbp_glo, index_g)
    !-
    var_name = 't2m_longterm'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,    1, 1, itime, &
         &                t2m_longterm, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'tau_longterm'
    CALL restput_p (rest_id_stomate, var_name, itime, tau_longterm)
    !-   
    var_name = 't2m_month'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,    1, 1, itime, &
                         t2m_month, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'Tseason'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,    1, 1, itime, &
         Tseason, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'Tseason_length'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,    1, 1, itime, &
         Tseason_length, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'Tseason_tmp'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,    1, 1, itime, &
         Tseason_tmp, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'Tmin_spring_time'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         Tmin_spring_time, 'scatter', nbp_glo, index_g)
    !-
    var_name = 't2m_week'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,    1, 1, itime, &
         &                t2m_week, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'tsoil_month'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nslm, 1, itime, &
         &                tsoil_month, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'precip_month'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,    1, 1, itime, &
         &                precip_month, 'scatter', nbp_glo, index_g)  
    !-
    ! 6 fire probability
    !-
    var_name = 'fireindex'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                fireindex, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'firelitter'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                firelitter, 'scatter', nbp_glo, index_g)
    !-
    ! 7 maximum and minimum moisture availabilities for tropic phenology
    !-
    var_name = 'maxmoistr_last'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                maxvegstress_lastyear, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'maxmoistr_this'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                maxvegstress_thisyear, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'minmoistr_last'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                minvegstress_lastyear, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'minmoistr_this'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                minvegstress_thisyear, 'scatter', nbp_glo, index_g)
    !-
    ! 8 maximum "weekly" GPP
    !-
    var_name = 'maxgppweek_lastyear'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                maxgppweek_lastyear, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'maxgppweek_thisyear'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                maxgppweek_thisyear, 'scatter', nbp_glo, index_g)
    !-
    ! 9 annual GDD0
    !-
    var_name = 'gdd0_thisyear'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
         &                gdd0_thisyear, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'gdd0_lastyear'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
         &                gdd0_lastyear, 'scatter', nbp_glo, index_g)
    !-
    ! 10 annual precipitation
    !-
    var_name = 'precip_thisyear'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
         &                precip_thisyear, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'precip_lastyear'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
         &                precip_lastyear, 'scatter', nbp_glo, index_g)
    !-
    ! 11 derived "biometeorological" variables
    !-
    var_name = 'gdd_m5_dormance'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                gdd_m5_dormance, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'gdd_from_growthinit'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &              gdd_from_growthinit, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'gdd_midwinter'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                gdd_midwinter, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'ncd_dormance'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                ncd_dormance, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'ngd_minus5'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                ngd_minus5, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'time_hum_min'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                time_hum_min, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'hum_min_dormance'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                hum_min_dormance, 'scatter', nbp_glo, index_g)

    !-
    ! 12 Plant status
    !-
    var_name = 'PFTpresent'
    WHERE ( PFTpresent(:,:) )
       PFTpresent_real = un
    ELSEWHERE
       PFTpresent_real = zero
    ENDWHERE
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                PFTpresent_real, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'turnover_time'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, nparts, itime, &
         &                turnover_time, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'adapted'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                adapted, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'regenerate'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                regenerate, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'npp_longterm'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                npp_longterm, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'croot_longterm'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                croot_longterm, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'n_reserve_longterm'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                n_reserve_longterm, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'lm_lastyearmax'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                lm_lastyearmax, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'lm_thisyearmax'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                lm_thisyearmax, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'maxfpc_lastyear'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                maxfpc_lastyear, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'maxfpc_thisyear'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                maxfpc_thisyear, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'turnover_longterm'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, nparts, nelements, itime, &
         &                turnover_longterm, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'gpp_week'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                gpp_week, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'maint_resp'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, nparts, itime, &
         &                   resp_maint_part, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'resp_maint_week'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                resp_maint_week, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'leaf_age'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, nleafages, itime, &
         &                  leaf_age, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'leaf_frac'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, nleafages, itime, &
         &                   leaf_frac, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'leaf_age_crit'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                   leaf_age_crit, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'lai_per_level'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, nlevels_tot, itime, &
         &                   lai_per_level, 'scatter', nbp_glo, index_g)
    !
    DO ipts=1,npts
       DO ivm=1,nvm
          DO ilev=1,nlevels_tot
             temp_array(ipts,ivm,ilev,1)=laieff_fit(ipts,ivm,ilev)%a
             temp_array(ipts,ivm,ilev,2)=laieff_fit(ipts,ivm,ilev)%b
             temp_array(ipts,ivm,ilev,3)=laieff_fit(ipts,ivm,ilev)%c
             temp_array(ipts,ivm,ilev,4)=laieff_fit(ipts,ivm,ilev)%d
             temp_array(ipts,ivm,ilev,5)=laieff_fit(ipts,ivm,ilev)%e
          ENDDO
       ENDDO
    ENDDO
    CALL restput_p (rest_id_stomate, 'laieff_fit', nbp_glo, nvm, nlevels_tot, nparams_laieff, itime, &
               &                   temp_array, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'plant_status'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                plant_status, 'scatter', nbp_glo, index_g)
    !- 
    var_name = 'when_growthinit'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                when_growthinit, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'age'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &                age, 'scatter', nbp_glo, index_g)
    !-
    ! 13 CO2
    !-
    var_name = 'resp_hetero'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                resp_hetero, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'resp_maint'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                resp_maint, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'resp_growth'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                resp_growth, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'co2_fire'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,  nvm, 1, itime, &
         &                co2_fire, 'scatter', nbp_glo, index_g)
    !-
    ! 14 vegetation distribution after last light competition
    !-
    var_name = 'veget_lastlight'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                veget_lastlight, 'scatter', nbp_glo, index_g)
    !-
    ! 15 establishment criteria
    !-
    var_name = 'everywhere'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                everywhere, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'need_adjacent'
    WHERE (need_adjacent(:,:))
       need_adjacent_real = un
    ELSEWHERE
       need_adjacent_real = zero
    ENDWHERE
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                need_adjacent_real, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'RIP_time'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                RIP_time, 'scatter', nbp_glo, index_g)
    !-
    ! 17 litter
    !-
    var_name = 'litter_c'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nlitt, nvm, nlevs, itime, &
        &                 litter(:,:,:,:,icarbon), 'scatter', nbp_glo, index_g)
    !-
    var_name = 'litter_n'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nlitt, nvm, nlevs,itime, &
        &                 litter(:,:,:,:,initrogen), 'scatter', nbp_glo, index_g)
    !-
    CALL restput_p (rest_id_stomate, 'litterfuel', nbp_glo, nlitt, nvm, nhour,nelements, itime, &
         &                 litterfuel(:,:,:,:,:), 'scatter', nbp_glo, index_g)
    !-
    var_name = 'dead_leaves'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,  nvm, nlitt, itime, &
        &                   dead_leaves, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'soil_carbon'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, ncarb, nvm, itime, &
         &                   som(:,:,:,icarbon), 'scatter', nbp_glo, index_g)
    !-
    var_name = 'soil_nitrogen'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, ncarb, nvm, itime, & 
         &                   som(:,:,:,initrogen), 'scatter', nbp_glo, index_g)
    !- 
    var_name = 'lignin_struc'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, nlevs, itime, &
         &                   lignin_struc, 'scatter', nbp_glo, index_g)
    !- 
    var_name = 'lignin_wood'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, nlevs, itime, &
         &                   lignin_wood, 'scatter', nbp_glo, index_g)
    !- 
    var_name = 'lignin_snag'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, nlevs, itime, &
         &                   lignin_snag, 'scatter', nbp_glo, index_g)
    !-
    ! 18 land cover change
    !-
    var_name = 'fco2_lu'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
         &                fco2_lu, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'fco2_wh'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
         &                fco2_wh, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'fco2_ha'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
         &                fco2_ha, 'scatter', nbp_glo, index_g)
    !-
    DO l = 1,nlctypes

       var_name = 'prod_s'//TRIM(lctype_str(l))
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, nshort+1, nelements, nlanduse, itime, &
            &                prod_s(:,:,:,:,l), 'scatter', nbp_glo, index_g)
       !-
       var_name = 'prod_m'//TRIM(lctype_str(l))
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, nmedium+1, nelements, nlanduse, itime, &
            &                prod_m(:,:,:,:,l), 'scatter', nbp_glo, index_g)
       !-
       var_name = 'prod_l'//TRIM(lctype_str(l))   
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, nlong+1, nelements, nlanduse, itime, &
            &                prod_l(:,:,:,:,l), 'scatter', nbp_glo, index_g)
       !-
       var_name = 'flux_s'//TRIM(lctype_str(l))
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, nshort+1, nelements, nlanduse, itime, &
            &                flux_s(:,:,:,:,l), 'scatter', nbp_glo, index_g)
       !-
       var_name = 'flux_m'//TRIM(lctype_str(l))
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, nmedium+1, nelements, nlanduse, itime, &
            &                flux_m(:,:,:,:,l), 'scatter', nbp_glo, index_g)
       !-
       var_name = 'flux_l'//TRIM(lctype_str(l))
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, nlong+1, nelements, nlanduse, itime, &
            &                flux_l(:,:,:,:,l), 'scatter', nbp_glo, index_g)
    END DO
    !-
    var_name = 'bm_to_litter'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, nparts, nelements, itime, &
         &                bm_to_litter, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'bm_to_litter_resid'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, nparts, nelements, itime, &
         &                bm_to_litter_resid, 'scatter', nbp_glo, index_g)
    !-
    CALL restput_p (rest_id_stomate, 'tree_to_litter', nbp_glo, nvm, nparts, nelements, itime, &
          &                tree_bm_to_litter, 'scatter', nbp_glo, index_g)
    !-
    CALL restput_p (rest_id_stomate, 'tree_to_litter_resid', nbp_glo, nvm, nparts, nelements, itime, &
          &                tree_bm_to_litter_resid, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'fDeforestToProduct'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              fDeforestToProduct, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'fLulccResidue'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              fLulccResidue, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'fHarvestToProduct'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &              fHarvestToProduct, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'carb_mass_total'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
         &              carb_mass_total, 'scatter', nbp_glo, index_g)
    !-
    ! Harvest pool
    DO l = 1,nlanduse
       var_name = 'harvest_pool_acc'//TRIM(luse_str(l))
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, ndia_harvest+1, nelements, itime, &
            harvest_pool_acc(:,:,:,:,l), 'scatter', nbp_glo, index_g)
       !-
       var_name = 'harvest_area_acc'//TRIM(luse_str(l))
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
            &                harvest_area_acc(:,:,l), 'scatter', nbp_glo, index_g)
    ENDDO
    
    !! C and N burried during land cover change 
    var_name = 'burried_litter'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nlitt, nlevs, nelements, itime, &
         &                burried_litter, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'burried_fresh_ltr'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nparts, nelements, itime, &
         &                burried_fresh_ltr, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'burried_fresh_som'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nparts, nelements, itime, &
         &                burried_fresh_som, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'burried_bact'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
         &                burried_bact, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'ni_acc'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
         &                ni_acc, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'burried_min_nitro'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nnspec, 1, itime, &
         &                burried_min_nitro, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'burried_som'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, ncarb, nelements, itime, &
         &                burried_som, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'burried_deepSOM_a'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, ngrnd, nelements, itime, &
         &                burried_deepSOM_a, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'burried_deepSOM_s'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, ngrnd, nelements, itime, &
         &                burried_deepSOM_s, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'burried_deepSOM_p'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, ngrnd, nelements, itime, &
         &                burried_deepSOM_p, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'area'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
         &                area, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'contfrac'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
         &                contfrac, 'scatter', nbp_glo, index_g)
    !-

    IF (ok_windthrow) THEN

       var_name = 'is_storm'
       WHERE ( is_storm(:) )
           is_storm_real = un
       ELSEWHERE
           is_storm_real = zero
       ENDWHERE
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
             &                is_storm_real, 'scatter', nbp_glo, index_g)
     
       var_name = 'count_storm'
       count_storm_real = REAL(count_storm,r_std)
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
            &                count_storm_real, 'scatter', nbp_glo, index_g)
       var_name = 'max_wind_speed_storm'
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
            &                max_wind_speed_storm, 'scatter', nbp_glo, index_g)
       ! JJ 2026: companion of max_wind_speed_storm
       var_name = 'max_wind_ratio_storm'
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
            &                max_wind_ratio_storm, 'scatter', nbp_glo, index_g)
       var_name = 'wind_speed_max_save'
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, wind_days, 1, itime, &
            &                wind_speed_max_save, 'scatter', nbp_glo, index_g)
       var_name = 'gap_area_save'
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, legacy_years_wind, itime, &
            &             gap_area_save, 'scatter', nbp_glo, index_g)
       !-
       var_name = 'total_ba_init'
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
            &             total_ba_init, 'scatter', nbp_glo, index_g)
    ENDIF
    ! Below variables are outside of ok_windthorw flag due to the restartability.
    var_name = 'wind_sum'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
         &                wind_sum, 'scatter', nbp_glo, index_g)
    var_name = 'wind_mean_save'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, legacy_years_wind, 1, itime, &
         &                wind_mean_save, 'scatter', nbp_glo, index_g)
    var_name = 'wind_ratio_sum_save'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, wind_days, 1, itime, &
         &                wind_ratio_sum_save, 'scatter', nbp_glo, index_g)
    ! JJ 2026: companion of wind_ratio_sum_save
    var_name = 'wind_ratio_max_save'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, wind_days, 1, itime, &
         &                wind_ratio_max_save, 'scatter', nbp_glo, index_g)

    !-
       var_name = 'sumTeff'
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm,1,itime, &
            &             sumTeff, 'scatter', nbp_glo, index_g)
       !-
       var_name = 'woody_litter_to_use'
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm,1,itime, &
            &             woody_litter_to_use, 'scatter', nbp_glo, index_g)
       !-
       var_name = 'beetle_diapause'
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm,1,itime, &
            &             beetle_diapause, 'scatter', nbp_glo, index_g)
       !-
       var_name = 'P_beetles_attacked_legacy'
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm,1, itime, &
               &     P_beetles_attacked_legacy, 'scatter', nbp_glo,index_g)

       DO l=1,legacy_years

          WRITE(pyear_str,'(I2.2)') l
          var_name = 'season_drought_'//pyear_str
          CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
               &     season_drought_legacy(:,:,l), 'scatter', nbp_glo, index_g)
          !-
          var_name = 'i_beetles_generation_'//pyear_str
          CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm,1, itime, &
               &     i_beetles_generation(:,:,l), 'scatter', nbp_glo, index_g)
          !-
          var_name = 'i_beetles_activity_'//pyear_str
          CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm,1, itime, &
               &     i_beetles_activity_legacy(:,:,l), 'scatter', nbp_glo, index_g)
          !-
       ENDDO

       DO l=1,legacy_years_wood
          WRITE(pyear_str,'(I2.2)') l
          var_name = 'wood_leftover_'//pyear_str
          CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
               &     wood_leftover_legacy(:,:,l), 'scatter', nbp_glo, index_g)
          !-
       ENDDO

       !-
       DO l=1,beetle_legacy
          WRITE(pyear_str,'(I2.2)') l
          var_name = 'B_beetles_kill_legacy_'//pyear_str
          CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm,1, itime, &
               &     B_beetles_kill_legacy(:,:,l), 'scatter', nbp_glo, index_g)
       ENDDO

    var_name = 'sugar_load'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm , 1, itime, &
         &              sugar_load, 'scatter', nbp_glo, index_g)

    var_name = 'harvest_cut'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm , 1, itime, &
         &              harvest_cut, 'scatter', nbp_glo, index_g)

    IF ( ok_soil_carbon_discretization ) THEN
       var_name= 'deepSOM_a'
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, ngrnd, nvm, nelements, itime, &
            deepSOM_a, 'scatter', nbp_glo, index_g)
       !-
       var_name= 'deepSOM_s'
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, ngrnd, nvm, nelements, itime, &
            deepSOM_s, 'scatter', nbp_glo, index_g)
       !-
       var_name= 'deepSOM_p' 
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, ngrnd, nvm, nelements, itime, &
            deepSOM_p, 'scatter', nbp_glo, index_g)
       !-
       var_name= 'O2_soil' 
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, ngrnd, nvm, itime, &
            O2_soil, 'scatter', nbp_glo, index_g)
       !-
       var_name= 'CH4_soil'
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, ngrnd, nvm, itime, &
            CH4_soil, 'scatter', nbp_glo, index_g)
       !-
       var_name= 'O2_snow'
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, nsnow, nvm, itime, &
            O2_snow, 'scatter', nbp_glo, index_g)
       !-
       var_name= 'CH4_snow'
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, nsnow, nvm, itime, &
            CH4_snow, 'scatter', nbp_glo, index_g)
       !-
       var_name= 'heat_Zimov'
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, ngrnd, nvm, itime, &
            heat_Zimov, 'scatter', nbp_glo, index_g)
       !-
       var_name= 'depth_organic_soil'
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
            depth_organic_soil, 'scatter', nbp_glo, index_g)
       !-
       CALL restput_p (rest_id_stomate, 'altmax_ind', nbp_glo, nvm, 1, itime, &
            altmax_ind, 'scatter', nbp_glo, index_g)
       CALL restput_p (rest_id_stomate, 'altmax_lastyear', nbp_glo, nvm, 1, itime, &
            altmax_lastyear, 'scatter', nbp_glo, index_g)
       CALL restput_p (rest_id_stomate, 'altmax_ind_lastyear', nbp_glo, nvm, 1, itime, &
            altmax_ind_lastyear, 'scatter', nbp_glo, index_g)
       !-
       var_name = 'fixed_cryoturb_depth'
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
            fixed_cryoturbation_depth, 'scatter', nbp_glo, index_g)
    ENDIF
    !-
    ! Altmax is used in hydrol also when ok_soil_carbon_discretization=n. Hence,
    ! it should be outside the IF (ok_soil_carbon_discretization) THEN
    var_name= 'altmax' 
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         altmax, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'nbp_accu_flux'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nelements, 1, itime, &
         &              nbp_accu_flux, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'nbp_pool_start'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nelements, 1, itime, &
         &              nbp_pool_start, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'Global_years'
    CALL restput_p (rest_id_stomate, var_name, itime, global_years)
    
    !-
    ! 19. Spinup
    !-
    IF (spinup_analytic) THEN

       var_name = 'ok_equilibrium'
       WHERE(ok_equilibrium(:))
          ok_equilibrium_real = un
       ELSEWHERE
          ok_equilibrium_real = zero
       ENDWHERE
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
            &               ok_equilibrium_real, 'scatter', nbp_glo, index_g)
       !-
       var_name = 'MatrixV'
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, nbpools, nbpools, itime, &
            &                MatrixV, 'scatter', nbp_glo, index_g)
       !-
       var_name = 'Vector_U'
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, nbpools, itime, &
            &                VectorU, 'scatter', nbp_glo, index_g)
       !-
       var_name = 'previous_stock'
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, nbpools, itime, &
            &                previous_stock, 'scatter', nbp_glo, index_g)
       !-
       var_name = 'current_stock'
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, nbpools, itime, &
            &                current_stock, 'scatter', nbp_glo, index_g)
       !-
       var_name = 'CN_longterm'
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, nbpools, itime, &
            &                CN_som_litter_longterm, 'scatter', nbp_glo, index_g)
       !-
       var_name = 'tau_CN_longterm'
       CALL restput_p (rest_id_stomate, var_name, itime, tau_CN_longterm)
       
    ENDIF !(spinup_analytic)

    !-
    var_name = 'KF'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
          &              KF(:,:), 'scatter', nbp_glo, index_g)
    !-
    var_name = 'k_latosa_adapt'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
          &              k_latosa_adapt(:,:), 'scatter', nbp_glo, index_g)
    !-
    var_name = 'rue_longterm'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
          &              rue_longterm(:,:), 'scatter', nbp_glo, index_g)
    !-
    var_name = 'cn_leaf_min_season' 
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, & 
         &              cn_leaf_min_season(:,:), 'scatter', nbp_glo, index_g) 
    !- 
    var_name = 'nstress_season' 
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, & 
         &              nstress_season(:,:), 'scatter', nbp_glo, index_g)
    !- 
    var_name = 'soil_n_min' 
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, nnspec, itime, & 
         &              soil_n_min, 'scatter', nbp_glo, index_g) 
    !- 
    var_name = 'p_O2' 
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, & 
         &              p_O2(:,:), 'scatter', nbp_glo, index_g) 
    !- 
    var_name = 'bact' 
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, & 
         &              bact(:,:), 'scatter', nbp_glo, index_g)
    !-
    var_name = 'age_stand'
    temp_real = REAL(age_stand,r_std)
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &              temp_real, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'age_stand_bm'                            ! STAND_AGE (Marie 2026), real
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &              age_stand_bm, 'scatter', nbp_glo, index_g)
    !-
    ! PROGRESSIVE_HARVEST (Marie 2026): pending fractional clearcut, consumed by
    ! age_class_distr at the start of next year. See the readrestart counterpart.
    IF (ok_progressive_harvest .AND. ALLOCATED(ph_cut_frac)) THEN
       CALL restput_p (rest_id_stomate, 'ph_cut_frac', nbp_glo, nvm, 1, itime, &
            &              ph_cut_frac, 'scatter', nbp_glo, index_g)
    ENDIF
    !-
    var_name = 'rotation_n'
    temp_real = REAL(rotation_n,r_std)  
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &              temp_real, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'last_cut'
    temp_real = REAL(last_cut,r_std)
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &              temp_real, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'mai'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &              mai, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'pai'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &              pai, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'coppice_dens'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &              coppice_dens, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'previous_wood_volume'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &              previous_wood_volume, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'mai_count'
    temp_real = REAL(mai_count,r_std)
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &              temp_real, 'scatter', nbp_glo, index_g)

    var_name = 'forest_managed'
    temp_real = REAL(forest_managed,r_std)
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &              temp_real, 'scatter', nbp_glo, index_g)

    IF (ok_change_species) THEN
       r_replant(:,:)=REAL(species_change_map(:,:),r_std)
       var_name = 'species_change_map'
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
            &              r_replant, 'scatter', nbp_glo, index_g)
    END IF
    
    r_replant(:,:)=REAL(fm_change_map(:,:),r_std)
    var_name = 'fm_change_map'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &              r_replant, 'scatter', nbp_glo, index_g)
    !-
    DO ipts=1,npts
       DO ivm=1,nvm
          IF(lpft_replant(ipts,ivm))THEN
             r_replant(ipts,ivm)=un
          ELSE
             r_replant(ipts,ivm)=zero
          ENDIF
       END DO
    END DO
    var_name = 'lpft_replant'
    CALL restput_p(rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &              r_replant, 'scatter', nbp_glo, index_g)
    !-
    CALL restput_p (rest_id_stomate, 'assim_param', nbp_glo, nvm, npco2, itime, &
        &                assim_param, 'scatter', nbp_glo, index_g)
    
    ! 21 Seasonal mean transmitted light for recruitment in DOFOCO  
    var_name = 'light_tran_season' 
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, & 
        &                light_tran_to_floor_season, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'daylight_count'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime,&
    &                   daylight_count, 'scatter', nbp_glo, index_g)
    !-
    ! Drought mortality
    IF (ok_hydrol_arch) THEN
       DO l = 1,nelements
          var_name = 'bio_ini_drought'//TRIM(element_str(l))
          CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, ncirc, nparts, itime, &
               biomass_init_drought(:,:,:,:,l), 'scatter', nbp_glo, index_g)
       ENDDO
       !-
       temp_real2(:,:,:) = zero
       WHERE ( kill_vessels(:,:,:) )
          temp_real2 = un
       ENDWHERE
       var_name = 'kill_vessels'
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, ncirc, itime, &
            &   temp_real2, 'scatter', nbp_glo, index_g)
       !-
       var_name = 'vessel_loss_previous'
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, ncirc, itime, &
            &  vessel_loss_previous(:,:,:), 'scatter', nbp_glo, index_g)
       !-
    END IF
    !-    
    var_name = 'grow_season_len'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &  grow_season_len(:,:), 'scatter', nbp_glo, index_g)
    !-
    var_name = 'doy_start_gs'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &  doy_start_gs(:,:), 'scatter', nbp_glo, index_g)
    !-
    var_name = 'doy_end_gs'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &  doy_end_gs(:,:), 'scatter', nbp_glo, index_g)
    !-
    var_name = 'mean_start_gs'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &  mean_start_gs(:,:), 'scatter', nbp_glo, index_g)

    !-
    IF ( ok_peat_NoDiscretisation ) THEN
    
       var_name = 'carbon_acro'
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
            &              carbon_acro, 'scatter', nbp_glo, index_g)
       
       var_name = 'carbon_cato'
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
            &              carbon_cato, 'scatter', nbp_glo, index_g)
       
       var_name = 'height_acro'
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
            &              height_acro, 'scatter', nbp_glo, index_g)
       
       var_name = 'deepSOM_peat'
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, ngrnd, nvm, nelements, itime, &
            deepSOM_peat, 'scatter', nbp_glo, index_g)
    ENDIF

    !-
    var_name = 'gpp_daily'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                gpp_daily, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'gpp_year'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                gpp_year, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'gpp_decade'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                gpp_decade, 'scatter', nbp_glo, index_g)
    !-
    ! Dynamic parameters. They are global so they are not passed through
    ! the argument list of the SUBROUTINE.
    var_name = 'pre_indust_ref_gpp'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &  pre_indust_ref_gpp(:,:), 'scatter', nbp_glo, index_g)

    !-
    ! AED_FEEDBACK (Marie 2026) — write the dynamic edge_length state as a 1D field.
    IF (ok_aed_feedback) THEN
       var_name = 'edge_length_dyn'
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
            &              edge_length_dyn, 'scatter', nbp_glo, index_g)
    ENDIF

    !-
    ! ROTATION_GROWTH (2026-08-17) — ecrire l'integrateur d'accroissement. Sans cette
    ! ecriture la lecture ne trouverait jamais rien : une periode vaut un an, donc la
    ! memoire de plusieurs decennies serait perdue a chaque pas. rotation_growth n'est PAS
    ! ecrite, elle est derivee de ces deux-la a chaque passage annuel.
    IF (ALLOCATED(dia_growth)) THEN
       var_name = 'dia_growth'
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
            &              dia_growth, 'scatter', nbp_glo, index_g)
       CALL restput_p (rest_id_stomate, 'circ_dia_last', nbp_glo, nvm, ncirc, itime, &
            &              circ_dia_last, 'scatter', nbp_glo, index_g)
       CALL restput_p (rest_id_stomate, 'circ_n_last', nbp_glo, nvm, ncirc, itime, &
            &              circ_n_last, 'scatter', nbp_glo, index_g)
       var_name = 'rotation_growth'
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
            &              rotation_growth, 'scatter', nbp_glo, index_g)
    ENDIF

    !-
    ! AED_REGROWTH (Marie 2026) — write the class-0 conveyor state.
    !-
    IF (ok_edge_from_age_class .AND. ALLOCATED(class0_area)) THEN
       CALL restput_p (rest_id_stomate, 'class0_area_v2', nbp_glo, nagent, n_recruit_max, itime, &
            &              class0_area, 'scatter', nbp_glo, index_g)
       IF (n_class0 > 0 .AND. ALLOCATED(class0_veg)) THEN
          CALL restput_p (rest_id_stomate, 'class0_veg', nbp_glo, nvmap, n_class0, itime, &
               &              class0_veg, 'scatter', nbp_glo, index_g)
          ! Guillaume M. -- slot count stored alongside, so readrestart can refuse a
          ! changed N_CLASS0 with a clear message instead of a generic restget abort.
          CALL restput_p (rest_id_stomate, 'n_class0', itime, REAL(n_class0,r_std))
       ENDIF
       CALL restput_p (rest_id_stomate, 'class0_year_clock_s', itime, class0_year_clock_s)
       IF (ALLOCATED(agent_flux)) &
            CALL restput_p (rest_id_stomate, 'aed_agent_flux', nbp_glo, nagent, 1, itime, &
            &              agent_flux, 'scatter', nbp_glo, index_g)
    ENDIF

    !--------------------------
  END SUBROUTINE writerestart
  !-
  !===
  !-

!! ================================================================================================================================
!! SUBROUTINE   : read_initsoc_file
!!
!>\BRIEF          
!!
!! DESCRIPTION : Read of SOM map to initialize C pools. This subroutine is called only when use_inisom=true and no restart file 
!!               was found.
!!                
!!
!! RECENT CHANGE(S) : None
!! 
!! MAIN OUTPUT VARIABLE(S): initSOC_a, initSOC_s, initSOC_p : SOC from observations
!!                          
!! REFERENCE(S) :
!!
!! FLOWCHART    : None 
!! \n 
!_ ================================================================================================================================

  SUBROUTINE read_initsoc_file(nbpt, lalo, neighbours, resolution, contfrac, veget_max, initSOC_a, initSOC_s, initSOC_p)

    !! 0. Variable and parameter declaration
    
    !! 0.1 Input variables
    
    INTEGER(i_std), INTENT(in)                    :: nbpt                  !! Number of points for which the data needs to be interpolated (unitless)             
    REAL(r_std), INTENT(in)                       :: lalo(nbpt,2)          !! Vector of latitude and longitudes (degree)        
    INTEGER(i_std), INTENT(in)                    :: neighbours(nbpt,NbNeighb)!! Vector of neighbours for each grid point (1=N,2=E,3=S,4=W)  
    REAL(r_std), INTENT(in)                       :: resolution(nbpt,2)    !! The size of each grid cell in X and Y (km)
    REAL(r_std), INTENT(in)                       :: contfrac(nbpt)        !! Fraction of land in each grid cell (unitless)   
    REAL(r_std),DIMENSION(nbpt,nvm),INTENT(in)    :: veget_max             !! Maximum fraction of vegetation type including non-biological fraction (unitless)

    !! 0.3 Modified variables
    REAL(r_std), DIMENSION(:,:,:), INTENT(out)    :: initSOC_a
    REAL(r_std), DIMENSION(:,:,:), INTENT(out)    :: initSOC_s
    REAL(r_std), DIMENSION(:,:,:), INTENT(out)    :: initSOC_p

    !! 0.4 Local variables

    INTEGER(i_std)                                :: nbvmax                !! nbvmax for interpolation (unitless) 
    CHARACTER(LEN=80)                             :: filename
    INTEGER(i_std)                                :: iml, jml, lml, tml    !! Indices
    INTEGER(i_std)                                :: fid, ib, ip, jp, fopt !! Indices
    INTEGER(i_std)                                :: ilf, ks               !! Indices
    INTEGER(i_std)                                :: ji, jg, jv                 
    REAL(r_std)                                   :: totarea               !! Help variable to compute average SOC
    REAL(r_std), ALLOCATABLE, DIMENSION(:)        :: lat_lu, lon_lu        !! Latitudes and longitudes read from input file
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)      :: lat_rel, lon_rel      !! Help variable to read file data and allocate memory
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)      :: mask_lu               !! Help variable to read file data and allocate memory
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:)   :: mask
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)      :: initSOC
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:)      :: initSOC_file        !! Help variable to read file data and allocate memory
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)      :: sub_area              !! Help variable to read file data and allocate memory
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:,:) :: sub_index             !! Help variable to read file data and allocate memory
    CHARACTER(LEN=30)                             :: callsign              !! Help variable to read file data and allocate memory
    CHARACTER(LEN=100)                            :: str                   !! Temporary string var 
    LOGICAL                                       :: ok_interpol           !! Optional return of aggregate_2d
    INTEGER                                       :: ALLOC_ERR             !! Help varialbe to count allocation error
!_
!================================================================================================================================

    !! 1. Open file and allocate memory

    ! Open file with SOM map

    !Config Key   = INITSOC_FILE
    !Config Desc  = File with soil carbon stocks for initialization
    !Config If    = OK_SOIL_CARBON_DISCRETIZATION, USE_INITSOM
    !Config Def   = initSOC.nc
    !Config Help  =
    !Config Units = [FILE]
    filename = 'initSOC.nc'
    CALL getin_p('INITSOC_FILE',filename)
    
    ! Read data from file
    IF (is_root_prc) CALL flininfo(filename, iml, jml, lml, tml, fid)
    CALL bcast(iml)
    CALL bcast(jml)
    CALL bcast(lml)
    CALL bcast(tml)
    
    IF (lml .NE. ngrnd) THEN
       WRITE(numout, *) 'In read_initsoc_file, ngrnd=', ngrnd, ', depth found in file=', lml
       CALL ipslerr_p(3, 'read_initSOC_file',    &
            'depth from the file must be the same as ngrnd', &
            filename,'' )
    ENDIF
    
    ALLOCATE(lon_lu(iml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'read_initSOC_file','Problem in allocation of variable lon_lu','','')
    
    ALLOCATE(lat_lu(jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'read_initSOC_file','Problem in allocation of variable lat_lu','','')
    
    ALLOCATE(mask_lu(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'read_initSOC_file','Pb in allocation for mask_lu','','')
    
    ALLOCATE(initSOC_file(iml,jml,lml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'read_initSOC_file','Pb in allocation for initSOC_file','','')
    
    ALLOCATE(initSOC(nbpt,lml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'read_initSOC_file','Pb in allocation for initSOC_file','','')
    
    IF (is_root_prc) THEN
       CALL flinget(fid, 'longitude', iml, 0, 0, 0, 1, 1, lon_lu)
       CALL flinget(fid, 'latitude', jml, 0, 0, 0, 1, 1, lat_lu)
       CALL flinget(fid, 'mask', iml, jml, 0, 0, 1, 1, mask_lu)
       CALL flinget(fid, 'soil_organic_carbon', iml, jml, lml, tml, 1, 1, initSOC_file)
       
       CALL flinclo(fid)
    ENDIF
    
    CALL bcast(lon_lu)
    CALL bcast(lat_lu)
    CALL bcast(mask_lu)
    CALL bcast(initSOC_file)
    
    ALLOCATE(lon_rel(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'read_initSOC_file','Pb in allocation for lon_rel','','')
    
    ALLOCATE(lat_rel(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'read_initSOC_file','Pb in allocation for lat_rel','','')
    
    ALLOCATE(mask(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'read_initSOC_file','Problem in allocation of variable mask','','')
    
    DO jp=1,jml
       lon_rel(:,jp) = lon_lu(:)
    ENDDO
    DO ip=1,iml
       lat_rel(ip,:) = lat_lu(:)
    ENDDO
    
    mask(:,:) = zero
    WHERE (mask_lu(:,:) > zero )
       mask(:,:) = un
    ENDWHERE
    
    ! Set nbvmax to 200 for interpolation
    ! This number is the dimension of the variables in which we store 
    ! the list of points of the source grid which fit into one grid box of the
    ! target. 
    nbvmax = 200
    callsign = 'soil organic carbon'
    
    ! Start interpolation
    ok_interpol=.FALSE.
    DO WHILE ( .NOT. ok_interpol )
       WRITE(numout,*) "Projection arrays for ",callsign," : "
       WRITE(numout,*) "nbvmax = ",nbvmax
       
       ALLOCATE(sub_area(nbpt,nbvmax), STAT=ALLOC_ERR)
       IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'read_initSOC_file','Pb in allocation for sub_area','','')
       sub_area(:,:)=zero
       
       ALLOCATE(sub_index(nbpt,nbvmax,2), STAT=ALLOC_ERR)
       IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'read_initSOC_file','Pb in allocation for sub_index','','')
       sub_index(:,:,:)=0
       
       CALL aggregate_p(nbpt, lalo, neighbours, resolution, contfrac, &
            iml, jml, lon_rel, lat_rel, mask, callsign, &
            nbvmax, sub_index, sub_area, ok_interpol)
       
       IF ( .NOT. ok_interpol ) THEN
          DEALLOCATE(sub_area)
          DEALLOCATE(sub_index)
          nbvmax = nbvmax * 2
       ENDIF
    ENDDO
    
    ! Compute the average
    initSOC(:,:) = zero
    DO ib = 1, nbpt
       fopt = COUNT(sub_area(ib,:) > zero)
       IF ( fopt > 0 ) THEN
          totarea = zero
          DO ilf = 1, fopt
             ip = sub_index(ib,ilf,1)
             jp = sub_index(ib,ilf,2)
             initSOC(ib,:) = initSOC(ib,:) + initSOC_file(ip,jp,:) * sub_area(ib,ilf)
             totarea = totarea + sub_area(ib,ilf)
          ENDDO
          ! Normalize
          initSOC(ib,:) = initSOC(ib,:)/totarea
       ELSE
          ! Set defalut value for points where the interpolation fail
          WRITE(numout,*) 'On point ', ib, ' no points were found for interpolation data. The grid cell is initialized to zero.'
          WRITE(numout,*) 'Location : ', lalo(ib,2), lalo(ib,1)
          initSOC(ib,:) = 0.
       ENDIF
    ENDDO

    !! Split initSOC into active, slow and passive pools
    DO jv = 1,nvm
       DO ji =1,nbpt
          DO jg =1,ngrnd
             initSOC_a(ji,jg,jv)=initSOC(ji,jg) * initSOM_frac_a 
             initSOC_s(ji,jg,jv)=initSOC(ji,jg) * initSOM_frac_s 
             initSOC_p(ji,jg,jv)=initSOC(ji,jg) * initSOM_frac_p 
          ENDDO
       ENDDO
    ENDDO

    DEALLOCATE (lat_lu)
    DEALLOCATE (lat_rel)
    DEALLOCATE (lon_lu)
    DEALLOCATE (lon_rel)
    DEALLOCATE (mask_lu)
    DEALLOCATE (mask)
    DEALLOCATE (initSOC_file)
    DEALLOCATE (initSOC)
    DEALLOCATE (sub_area)
    DEALLOCATE (sub_index)

  END SUBROUTINE read_initsoc_file


!! ================================================================================================================================
!! SUBROUTINE   : read_initson_file
!!
!>\BRIEF          
!!
!! DESCRIPTION : Read of SOM map to initialize C pools. This subroutine is called only when use_inisom=true and no restart file 
!!               was found.
!!                
!!
!! RECENT CHANGE(S) : None
!! 
!! MAIN OUTPUT VARIABLE(S): initSON_a, initSON_s, initSON_p : SON from observations
!!                          
!! REFERENCE(S) :
!!
!! FLOWCHART    : None 
!! \n 
!_ ================================================================================================================================

  SUBROUTINE read_initson_file(nbpt, lalo, neighbours, resolution, contfrac, veget_max, initSOC_a, initSOC_s, initSOC_p, initSON_a, initSON_s, initSON_p)

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                    :: nbpt                  !! Number of points for which the data needs to be interpolated (unitless)             
    REAL(r_std), INTENT(in)                       :: lalo(nbpt,2)          !! Vector of latitude and longitudes (degree)        
    INTEGER(i_std), INTENT(in)                    :: neighbours(nbpt,NbNeighb)!! Vector of neighbours for each grid point (1=N,2=E,3=S,4=W)  
    REAL(r_std), INTENT(in)                       :: resolution(nbpt,2)    !! The size of each grid cell in X and Y (km)
    REAL(r_std), INTENT(in)                       :: contfrac(nbpt)        !! Fraction of land in each grid cell (unitless)   
    REAL(r_std),DIMENSION(nbpt,nvm),INTENT(in)    :: veget_max             !! Maximum fraction of vegetation type including non-biological fraction (unitless)
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)     :: initSOC_a
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)     :: initSOC_s
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)     :: initSOC_p

    !! 0.3 Modified variables
    REAL(r_std), DIMENSION(:,:,:), INTENT(out)    :: initSON_a
    REAL(r_std), DIMENSION(:,:,:), INTENT(out)    :: initSON_s
    REAL(r_std), DIMENSION(:,:,:), INTENT(out)    :: initSON_p

    !! 0.4 Local variables

    INTEGER(i_std)                                :: nbvmax                !! nbvmax for interpolation (unitless) 
    CHARACTER(LEN=80)                             :: filename
    INTEGER(i_std)                                :: iml, jml, lml, tml    !! Indices
    INTEGER(i_std)                                :: fid, ib, ip, jp, fopt !! Indices
    INTEGER(i_std)                                :: ilf, ks               !! Indices
    INTEGER(i_std)                                :: ji, jg, jv
    REAL(r_std)                                   :: totarea               !! Help variable to compute average SOC
    REAL(r_std), ALLOCATABLE, DIMENSION(:)        :: lat_lu, lon_lu        !! Latitudes and longitudes read from input file
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)      :: lat_rel, lon_rel      !! Help variable to read file data and allocate memory
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)      :: mask_lu               !! Help variable to read file data and allocate memory
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:)   :: mask
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)      :: initSON
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:)    :: initSON_file          !! Help variable to read file data and allocate memory
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)      :: sub_area              !! Help variable to read file data and allocate memory
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:,:) :: sub_index             !! Help variable to read file data and allocate memory
    CHARACTER(LEN=30)                             :: callsign              !! Help variable to read file data and allocate memory
    CHARACTER(LEN=100)                            :: str                   !! Temporary string var 
    LOGICAL                                       :: ok_interpol           !! Optional return of aggregate_2d
    INTEGER                                       :: ALLOC_ERR             !! Help varialbe to count allocation error
!_
!================================================================================================================================

    !! 1. Open file and allocate memory

    ! Open file with SON map

    !Config Key   = INITSON_FILE
    !Config Desc  = File with soil nitrogen stocks
    !Config If    = OK_SOIL_CARBON_DISCRETIZATION, USE_INITSOM
    !Config Def   = initSON.nc
    !Config Help  =
    !Config Units = [FILE]
    filename = 'initSON.nc'

    CALL getin_p('INITSON_FILE',filename)

    ! Read data from file
    IF (is_root_prc) CALL flininfo(filename, iml, jml, lml, tml, fid)
    CALL bcast(iml)
    CALL bcast(jml)
    CALL bcast(lml)
    CALL bcast(tml)

    IF (lml .NE. ngrnd) THEN
       WRITE(numout, *) 'In read_initson_file, ngrnd=', ngrnd, ', depth found in file=', lml
       CALL ipslerr_p(3, 'read_initSON_file',    &
            'depth from the file must be the same as ngrnd', &
            filename,'' )
    ENDIF

    ALLOCATE(lon_lu(iml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'read_initSON_file','Problem in allocation of variable lon_lu','','')

    ALLOCATE(lat_lu(jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'read_initSON_file','Problem in allocation of variable lat_lu','','')

    ALLOCATE(mask_lu(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'read_initSON_file','Pb in allocation for mask_lu','','')

    ALLOCATE(initSON_file(iml,jml,lml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'read_initSON_file','Pb in allocation for initSON_file','','')

    ALLOCATE(initSON(nbpt,lml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'read_initSON_file','Pb in allocation for initSON','','')

    IF (is_root_prc) THEN
       CALL flinget(fid, 'longitude', iml, 0, 0, 0, 1, 1, lon_lu)
       CALL flinget(fid, 'latitude', jml, 0, 0, 0, 1, 1, lat_lu)
       CALL flinget(fid, 'mask', iml, jml, 0, 0, 1, 1, mask_lu)
       CALL flinget(fid, 'soil_organic_nitrogen', iml, jml, lml, tml, 1, 1, initSON_file)

       CALL flinclo(fid)
    ENDIF

    CALL bcast(lon_lu)
    CALL bcast(lat_lu)
    CALL bcast(mask_lu)
    CALL bcast(initSON_file)

    ALLOCATE(lon_rel(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'read_initSON_file','Pb in allocation for lon_rel','','')

    ALLOCATE(lat_rel(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'read_initSON_file','Pb in allocation for lat_rel','','')

    ALLOCATE(mask(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'read_initSON_file','Problem in allocation of variable mask','','')

    DO jp=1,jml
       lon_rel(:,jp) = lon_lu(:)
    ENDDO
    DO ip=1,iml
       lat_rel(ip,:) = lat_lu(:)
    ENDDO

    mask(:,:) = zero
    WHERE (mask_lu(:,:) > zero )
       mask(:,:) = un
    ENDWHERE

   ! Set nbvmax to 200 for interpolation
    ! This number is the dimension of the variables in which we store 
    ! the list of points of the source grid which fit into one grid box of the
    ! target. 
    nbvmax = 200 
    callsign = 'soil organic nitrogen'

    ! Start interpolation
    ok_interpol=.FALSE.
    DO WHILE ( .NOT. ok_interpol )
       IF (printlev >= 2) WRITE(numout,*) "Projection arrays for ",callsign," : "
       IF (printlev >= 2) WRITE(numout,*) "nbvmax = ",nbvmax

       ALLOCATE(sub_area(nbpt,nbvmax), STAT=ALLOC_ERR)
       IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'read_initSON_file','Pb in allocation for sub_area','','')
       sub_area(:,:)=zero

       ALLOCATE(sub_index(nbpt,nbvmax,2), STAT=ALLOC_ERR)
       IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'read_initSON_file','Pb in allocation for sub_index','','')
       sub_index(:,:,:)=0

       CALL aggregate_p(nbpt, lalo, neighbours, resolution, contfrac, &
            iml, jml, lon_rel, lat_rel, mask, callsign, &
            nbvmax, sub_index, sub_area, ok_interpol)

       IF ( .NOT. ok_interpol ) THEN
          DEALLOCATE(sub_area)
          DEALLOCATE(sub_index)
          nbvmax = nbvmax * 2
       ENDIF
    ENDDO

    ! Compute the average
    initSON(:,:) = zero
    DO ib = 1, nbpt
       fopt = COUNT(sub_area(ib,:) > zero)
       IF ( fopt > 0 ) THEN
          totarea = zero
          DO ilf = 1, fopt
             ip = sub_index(ib,ilf,1)
             jp = sub_index(ib,ilf,2)
             initSON(ib,:) = initSON(ib,:) + initSON_file(ip,jp,:) * sub_area(ib,ilf)
             totarea = totarea + sub_area(ib,ilf)
          ENDDO
          ! Normalize
          initSON(ib,:) = initSON(ib,:)/totarea
       ELSE
          ! Set defalut value for points where the interpolation fail
          IF (printlev >= 2) WRITE(numout,*) 'On point ', ib, ' no points were found for interpolation data. The grid cell is initialized to the carbon values divided by the CN target of each pool.'
          IF (printlev >= 2) WRITE(numout,*) 'Location : ', lalo(ib,2), lalo(ib,1)
          initSON(ib,:) = 0.
       ENDIF
    ENDDO


    !! Split initSON into active, slow and passive pools

    DO jv = 1,nvm
       DO ji =1,nbpt
          DO jg =1,ngrnd
             IF (initSON(ji,jg) .EQ. zero) THEN
                initSON_a(ji,jg,jv)=initSOC_a(ji,jg,jv)/ CN_target_iactive_ref
                initSON_s(ji,jg,jv)=initSOC_s(ji,jg,jv)/ CN_target_islow_ref
                initSON_p(ji,jg,jv)=initSOC_p(ji,jg,jv)/ CN_target_ipassive_ref
             ELSE
                initSON_a(ji,jg,jv)=initSON(ji,jg) * initSOM_frac_a  !! Same fraction as for SOC
                initSON_s(ji,jg,jv)=initSON(ji,jg) * initSOM_frac_s
                initSON_p(ji,jg,jv)=initSON(ji,jg) * initSOM_frac_p
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    DEALLOCATE (lat_lu)
    DEALLOCATE (lat_rel)
    DEALLOCATE (lon_lu)
    DEALLOCATE (lon_rel)
    DEALLOCATE (mask_lu)
    DEALLOCATE (mask)
    DEALLOCATE (initSON_file)
    DEALLOCATE (initSON)
    DEALLOCATE (sub_area)
    DEALLOCATE (sub_index)

  END SUBROUTINE read_initson_file  

END MODULE stomate_io
