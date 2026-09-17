! =================================================================================================================================
! MODULE       : stomate
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       Groups the subroutines that: (1) initialize all variables in 
!! stomate, (2) read and write forcing files of stomate and the soil component,
!! (3) aggregates and convert variables to handle the different time steps 
!! between sechiba and stomate, (4) call subroutines that govern major stomate
!! processes (litter,\ soil, and vegetation dynamics) and (5) structures these tasks 
!! in stomate_main
!!
!!\n DESCRIPTION : None
!!
!! RECENT CHANGE(S) : None
!!
!! REFERENCE(S)	: None
!!
!! SVN :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_stomate/stomate.f90 $
!! $Date: 2026-04-24 17:38:54 +0200 (ven. 24 avril 2026) $
!! $Revision: 9501 $
!! \n
!_ ================================================================================================================================

MODULE stomate

  ! Modules used:
  USE netcdf
  USE defprec
  USE grid
  USE time, ONLY : one_day, one_year, dt_sechiba, &
                   dt_stomate, ts_annual_proc, LastTsMonth, &
                   FirstTsDay, FirstTsYear, year_end, month_end, &
                   day_end, sec_end
  USE constantes
  USE constantes_soil
  USE vertical_soil_var
  USE pft_parameters
  USE dynamic_parameters
  USE structures
  USE sapiens_agriculture,ONLY : sapiens_agriculture_initialize
  USE sapiens_forestry,   ONLY : sapiens_forestry_set_fm,  &
                                 sapiens_forestry_set_species_change, &
                                 sapiens_forestry_set_desired_fm
  USE stomate_io
  USE stomate_data
  USE stomate_season
  USE stomate_lpj
  USE stomate_litter
  USE stomate_vmax
  USE stomate_som_dynamics
  USE stomate_resp
  USE mod_orchidee_para
  USE ioipsl_para 
  USE xios_orchidee
  USE function_library,   ONLY : cc_to_lai, wood_to_qmheight, wood_to_height, &
                                 check_vegetation_area, check_mass_balance, &
                                 check_pixel_area, get_printlev, weibull_class_dist, &
                                 fm_polynomial_parameter_fitting
  USE matrix_resolution
  USE utils
  USE stomate_soil_carbon_discretization  
  USE stomate_io_soil_carbon_discretization   
  USE stomate_laieff
  USE stomate_spitfire

  IMPLICIT NONE

  ! Private & public routines

  PRIVATE
  PUBLIC stomate_main,stomate_clear, stomate_initialize, stomate_finalize

  INTERFACE stomate_accu
     MODULE PROCEDURE stomate_accu_r1d, stomate_accu_r2d, stomate_accu_r3d, stomate_accu_r4d
  END INTERFACE

  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:,:):: som_surf             !! Carbon pool integrated to over surface soils: active, slow, or passive
!$OMP THREADPRIVATE(som_surf)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: age                  !! Age of PFT it normalized by biomass - can increase and
                                                                         !! decrease - (years)
!$OMP THREADPRIVATE(age)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: adapted              !! Winter too cold for PFT to survive (0-1, unitless)
!$OMP THREADPRIVATE(adapted)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: regenerate           !! Winter sufficiently cold to produce viable seeds 
                                                                         !! (0-1, unitless)
!$OMP THREADPRIVATE(regenerate)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: everywhere           !! Is the PFT everywhere in the grid box or very localized 
                                                                         !! (after its intoduction)
!$OMP THREADPRIVATE(everywhere)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: fireindex            !! Probability of fire (unitless)
!$OMP THREADPRIVATE(fireindex)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: veget_lastlight      !! Vegetation fractions (on ground) after last light 
                                                                         !! competition (unitless) 
!$OMP THREADPRIVATE(veget_lastlight)
  REAL(r_std), ALLOCATABLE,SAVE,DIMENSION(:,:)   :: fpc_max              !! "maximal" coverage fraction of a grid box (LAI -> 
                                                                         !! infinity) on ground. [??CHECK??] It's set to zero here, 
                                                                         !! and then is used once in lpj_light.f90 to test if 
                                                                         !! fpc_nat is greater than it. Something seems missing
!$OMP THREADPRIVATE(fpc_max)
  LOGICAL,ALLOCATABLE,SAVE,DIMENSION(:,:)        :: PFTpresent           !! PFT exists (equivalent to veget > 0 for natural PFTs)
!$OMP THREADPRIVATE(PFTpresent)
  REAL,ALLOCATABLE,SAVE,DIMENSION(:,:)           :: plant_status         !! Growth and phenological status of the plant
                                                                         !! Different stati defined in constantes
!$OMP THREADPRIVATE(plant_status)

  LOGICAL,ALLOCATABLE,SAVE,DIMENSION(:,:)        :: need_adjacent        !! This PFT needs to be in present in an adjacent gridbox 
                                                                         !! if it is to be introduced in a new gridbox
!$OMP THREADPRIVATE(need_adjacent)
!--
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: vegstress_day     !! Daily plant available water -root profile weighted 
                                                                         !! (0-1, unitless)
!$OMP THREADPRIVATE(vegstress_day)

  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: stressed_daily       !! Accumulated proxy for stressed ecosystem functioning
                                                                         !! see variable stressed defined in sechiba
!$OMP THREADPRIVATE(stressed_daily)

  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: unstressed_daily     !! Accumulated proxy for unstressed ecosystem functioning
                                                                         !! see variable stressed defined in sechiba 
!$OMP THREADPRIVATE(unstressed_daily)

  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:,:,:) :: biomass_init_drought !! Biomass of heartwood or sapwood before onset of drought. 
                                                                          !! Used to compute turnover on same reference biomass in 
                                                                          !! stomate_turnover.f90. Should remain the same along one 
                                                                          !! entire drought episode and be updated inbetween droughts.
!$OMP THREADPRIVATE(biomass_init_drought)

  LOGICAL,ALLOCATABLE,SAVE,DIMENSION(:,:,:)      :: kill_vessels         !! Flag to kill vessels at the end of the day when there is embolism.
 !$OMP THREADPRIVATE(kill_vessels)

  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:)  :: vessel_loss_previous !! Proportion of conductivity lost due to cavitation, accumulated
                                                                         !!  on the previous day (no unit). Used to compute vl_diff_daily.
!$OMP THREADPRIVATE(vessel_loss_previous)
  
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:)  :: vessel_loss_daily    !! Proportion of conductivity lost due to cavitation in the xylem, 
                                                                         !! accumulated per day (no unit). See variable vessel_loss defined 
                                                                         !! in sechiba.f90.
!$OMP THREADPRIVATE(vessel_loss_daily)

  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: daylight             !! Time steps dt_radia during daylight
!$OMP THREADPRIVATE(daylight)

 REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)     :: daylight_count    !! Time steps dt_radia during daylight and when there is growth (gpp>0)
!$OMP THREADPRIVATE(daylight_count)

  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: transpir_daily       !! Daily demand of water for transpiration @tex $(mm dt^{-1})$ @endtex
!$OMP THREADPRIVATE(transpir_daily) 


  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: vegstress_week    !! "Weekly" plant available water -root profile weighted
                                                                         !! (0-1, unitless)
!$OMP THREADPRIVATE(vegstress_week)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: vegstress_month   !! "Monthly" plant available water -root profile weighted
                                                                         !! (0-1, unitless)
!$OMP THREADPRIVATE(vegstress_month)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: vegstress_season  !! Mean growing season moisture availability (used for 
                                                                         !! allocation response)
!$OMP THREADPRIVATE(vegstress_season)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: maxvegstress_lastyear   !! Last year's max plant available water -root profile 
                                                                         !! weighted (0-1, unitless)
!$OMP THREADPRIVATE(maxvegstress_lastyear)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: maxvegstress_thisyear   !! This year's max plant available water -root profile 
                                                                         !! weighted (0-1, unitless) 
!$OMP THREADPRIVATE(maxvegstress_thisyear)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: minvegstress_lastyear   !! Last year's min plant available water -root profile 
                                                                         !! weighted (0-1, unitless)  
!$OMP THREADPRIVATE(minvegstress_lastyear)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: minvegstress_thisyear   !! This year's minimum plant available water -root profile
                                                                         !! weighted (0-1, unitless)
!$OMP THREADPRIVATE(minvegstress_thisyear)
!---  
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: t2m_daily            !! Daily air temperature at 2 meter (K)
!$OMP THREADPRIVATE(t2m_daily)

  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: Tseason              !! "seasonal" 2 meter temperatures (K)
!$OMP THREADPRIVATE(Tseason)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: Tseason_length       !! temporary variable to calculate Tseason
!$OMP THREADPRIVATE(Tseason_length)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: Tseason_tmp          !! temporary variable to calculate Tseason
!$OMP THREADPRIVATE(Tseason_tmp)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: Tmin_spring_time     !! Number of days after begin_leaves (leaf onset) 
!$OMP THREADPRIVATE(Tmin_spring_time)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: t2m_week             !! Mean "weekly" (default 7 days) air temperature at 2 
                                                                         !! meter (K)  
!$OMP THREADPRIVATE(t2m_week)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: t2m_month            !! Mean "monthly" (default 20 days) air temperature at 2 
                                                                         !! meter (K)
!$OMP THREADPRIVATE(t2m_month)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: t2m_longterm         !! Mean "Long term" (default 3 years) air temperature at 
                                                                         !! 2 meter (K) 
!$OMP THREADPRIVATE(t2m_longterm)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: t2m_min_daily        !! Daily minimum air temperature at 2 meter (K)
!$OMP THREADPRIVATE(t2m_min_daily)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: t2m_max_daily        !! Daily maximum air temperature at 2 meter (K)
!$OMP THREADPRIVATE(t2m_max_daily)

  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: tsurf_daily          !! Daily surface temperatures (K)
!$OMP THREADPRIVATE(tsurf_daily)

!---
! variables added for windthrow module  ---
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)       :: wind_max_daily      !! Daily maximum wind speed at 2 meter (ms-1)
!$OMP THREADPRIVATE(wind_max_daily)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)       :: wind_ratio_sum      !! Sum of all the ratio of current wind to the longterm wind
!$OMP THREADPRIVATE(wind_ratio_sum)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)       :: wind_sum            !! Sum of all the wind speed of the day
!$OMP THREADPRIVATE(wind_sum)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)       :: max_wind_speed_storm!! Daily maximum wind speed at 2 meter during a storm event (ms-1)
!$OMP THREADPRIVATE(max_wind_speed_storm)
! Guillaume M. -- Companion of max_wind_speed_storm: storm state, kept in the restart file.
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)       :: max_wind_ratio_storm!! Daily maximum wind ratio during a storm event
!$OMP THREADPRIVATE(max_wind_ratio_storm)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)     :: wind_ratio_sum_save !! Saved values of wind_ratio_sum for wind_days
!$OMP THREADPRIVATE(wind_ratio_sum_save)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)     :: wind_speed_max_save !! Saved values of wind_max_daily for wind_days
!$OMP THREADPRIVATE(wind_speed_max_save)
! Guillaume M. -- Companion of wind_ratio_sum_save: rolling buffer, kept in the restart file.
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)     :: wind_ratio_max_save !! Saved values of wind_ratio_max for wind_days
!$OMP THREADPRIVATE(wind_ratio_max_save)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)     :: wind_mean_save !! Saved values of annual wind average for legacy years
!$OMP THREADPRIVATE(wind_mean_save)
  INTEGER(i_std),ALLOCATABLE,SAVE,DIMENSION(:)    :: count_storm         !! Number of days after a storm 
!$OMP THREADPRIVATE(count_storm)
  LOGICAL,ALLOCATABLE,SAVE,DIMENSION(:)           :: is_storm            !! Are we in a storm event ? 
!$OMP THREADPRIVATE(is_storm)

  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)       :: wind_max            !! Temporary daily maximum speed used to calculate wind_max_daily (ms-1)
!$OMP THREADPRIVATE(wind_max)
! Guillaume M. -- Transient daily accumulator, companion of wind_max, buffered into wind_ratio_max_save.
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)       :: wind_ratio_max      !! Temporary daily maximum wind ratio used to calculate wind_ratio_max_save
!$OMP THREADPRIVATE(wind_ratio_max)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)       :: soil_temp_daily     !! Daily maximum soil temperature at 0.8 meter below ground(K)
!$OMP THREADPRIVATE(soil_temp_daily)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)       :: soil_max_daily      !! Temporary daily maximum soil temperature used to calculate soil_temp_speed_daily (ms-1)
!$OMP THREADPRIVATE(soil_max_daily)
!---
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: precip_daily         !! Daily precipitations sum @tex $(mm day^{-1})$ @endtex
!$OMP THREADPRIVATE(precip_daily)

  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: vpd_daily_mean       !! Daily mean VPD @tex $(mm day^{-1})$ @endtex
!$OMP THREADPRIVATE(vpd_daily_mean)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: vpd_daily_max       !! Daily max VPD @tex $(mm day^{-1})$ @endtex
!$OMP THREADPRIVATE(vpd_daily_max)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: vpd_mean_week       !! Weekly mean Daily mean VPD
!$OMP THREADPRIVATE(vpd_mean_week)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: vpd_max_week        !! Weekly mean Daily maximum VPD 
!$OMP THREADPRIVATE(vpd_max_week)

  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: precip_lastyear      !! Last year's annual precipitation sum                                                                                                                                                 !! @tex $??(mm year^{-1})$ @endtex
!$OMP THREADPRIVATE(precip_lastyear)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: precip_thisyear      !! This year's annual precipitation sum 
                                                                         !! @tex $??(mm year^{-1})$ @endtex 
!$OMP THREADPRIVATE(precip_thisyear)
!---
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: tsoil_daily          !! Daily soil temperatures (K)
!$OMP THREADPRIVATE(tsoil_daily)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: tsoil_month          !! Soil temperatures at each soil layer integrated over a
                                                                         !! month (K) 
!$OMP THREADPRIVATE(tsoil_month)
!--- 
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: precip_month         !! Last month's precipitation sum (mm month^{-1})
!$OMP THREADPRIVATE(precip_month)

  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: litterhum_daily      !! Daily litter humidity (0-1, unitless)
!$OMP THREADPRIVATE(litterhum_daily)
!---
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:)  :: control_moist        !! Moisture control of heterotrophic respiration 
                                                                         !! (0-1, unitless)
!$OMP THREADPRIVATE(control_moist)
 REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:)   :: drainage             !! Fraction of water lost from the soil column by leaching (-)  
!$OMP THREADPRIVATE(drainage) 
 REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)     :: drainage_daily       !! Daily Fraction of water lost from the soil column by leaching (-)  
!$OMP THREADPRIVATE(drainage_daily) 
 REAL(r_std),ALLOCATABLE,SAVE, DIMENSION(:,:)    :: n_mineralisation_d   !! net nitrogen mineralisation of decomposing SOM (gN/m**2/day), assumed to be NH4  
!$OMP THREADPRIVATE(n_mineralisation_d) 
 REAL(r_std),ALLOCATABLE,SAVE, DIMENSION(:,:,:)  :: plant_n_uptake_daily !! Uptake of soil N by plants (gN/m**2/day)  
!$OMP THREADPRIVATE(plant_n_uptake_daily) 
 REAL(r_std), ALLOCATABLE,SAVE,DIMENSION(:,:)    :: atm_to_immob_daily   !! Nitrogen taken from the atmosphere to support immobilisation
!$OMP THREADPRIVATE(atm_to_immob_daily)
 REAL(r_std), ALLOCATABLE,SAVE,DIMENSION(:,:,:)  :: leaching_daily       !! Mineral nitrogen leached from the soil(g/m**2/day)
!$OMP THREADPRIVATE(leaching_daily) 
 REAL(r_std), ALLOCATABLE,SAVE,DIMENSION(:,:,:)  :: emission_daily       !! Volatile losses of nitrogen (gN/m**2/day)
!$OMP THREADPRIVATE(emission_daily)  
 REAL(r_std), ALLOCATABLE,SAVE,DIMENSION(:,:,:)  :: n_input_daily        !! Fertilizer, deposition and biological fixation of nitrogen (gN/m**2/day)
!$OMP THREADPRIVATE(n_input_daily)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:)  :: control_temp         !! Temperature control of heterotrophic respiration at the
                                                                         !! different soil levels (0-1, unitless)
!$OMP THREADPRIVATE(control_temp)
!---
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: gdd_init_date        !! inital date for gdd count 
!$OMP THREADPRIVATE(gdd_init_date)

  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: gdd_from_growthinit  !! gdd from beginning of season (C)
!$OMP THREADPRIVATE(gdd_from_growthinit)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: gdd0_lastyear        !! Last year's annual Growing Degree Days,
                                                                         !! threshold 0 deg C (K) 
!$OMP THREADPRIVATE(gdd0_lastyear)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: gdd0_thisyear        !! This year's annual Growing Degree Days,
                                                                         !! threshold 0 deg C (K)
!$OMP THREADPRIVATE(gdd0_thisyear)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: gdd_m5_dormance      !! Growing degree days for onset of growing season, 
                                                                         !! threshold -5 deg C (K)
!$OMP THREADPRIVATE(gdd_m5_dormance)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: gdd_midwinter        !! Growing degree days for onset of growing season, 
                                                                         !! since midwinter (K)
!$OMP THREADPRIVATE(gdd_midwinter)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: ncd_dormance         !! Number of chilling days since leaves were lost (days) 
!$OMP THREADPRIVATE(ncd_dormance)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: ngd_minus5           !! Number of growing days, threshold -5 deg C (days)
!$OMP THREADPRIVATE(ngd_minus5)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: hum_min_dormance     !! Minimum moisture during dormance (0-1, unitless) 
!$OMP THREADPRIVATE(hum_min_dormance)
!---
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: gpp_daily            !! Daily gross primary productivity per ground area 
                                                                         !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex
!$OMP THREADPRIVATE(gpp_daily) 
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: resp_maint_week      !! Mean "weekly" (default 7 days) maintenance respiration
                                                                         !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex
!$OMP THREADPRIVATE(resp_maint_week) 
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: gpp_week             !! Mean "weekly" (default 7 days) GPP  
                                                                         !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex
!$OMP THREADPRIVATE(gpp_week)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: gpp_year             !! Mean "annual" (default 365 days) GPP  
                                                                         !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex
!$OMP THREADPRIVATE(gpp_year)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: gpp_decade           !! Mean "decadal" (default 10 years) GPP  
                                                                         !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex
!$OMP THREADPRIVATE(gpp_decade)
  
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: maxgppweek_lastyear  !! Last year's maximum "weekly" GPP  
                                                                         !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex 
!$OMP THREADPRIVATE(maxgppweek_lastyear)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: maxgppweek_thisyear  !! This year's maximum "weekly" GPP  
                                                                         !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex  
!$OMP THREADPRIVATE(maxgppweek_thisyear)
!---
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: npp_daily            !! Daily net primary productivity per ground area 
                                                                         !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex 
!$OMP THREADPRIVATE(npp_daily)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: npp_longterm         !! "Long term" (default 3 years) net primary productivity 
                                                                         !! per ground area  
                                                                         !! @tex $(gC m^{-2} year^{-1})$ @endtex   
!$OMP THREADPRIVATE(npp_longterm)
!---
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: croot_longterm       !! "Long term" (default 3 years) root carbon mass
                                                                         !! per ground area
                                                                         !! @tex $(gC m^{-2} year^{-1})$ @endtex
!$OMP THREADPRIVATE(croot_longterm)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: n_reserve_longterm   !! "Long term" (default 3 years) actual to potential N
                                                                         !! reserve pool (0-1, unitless)
!$OMP THREADPRIVATE(n_reserve_longterm)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:)  :: resp_maint_part_radia!! Maintenance respiration of different plant parts per 
                                                                         !! total ground area at Sechiba time step  
                                                                         !! @tex $(gC m^{-2} dt_sechiba^{-1})$ @endtex
!$OMP THREADPRIVATE(resp_maint_part_radia)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:)  :: resp_maint_part      !! Maintenance respiration of different plant parts per
                                                                         !! total ground area at Stomate time step 
                                                                         !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex
!$OMP THREADPRIVATE(resp_maint_part)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: resp_maint_radia     !! Maintenance respiration per ground area at Sechiba time
                                                                         !! step
                                                                         !! @tex $(gC m^{-2} dt_sechiba^{-1})$ @endtex
!$OMP THREADPRIVATE(resp_maint_radia)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: resp_maint_d         !! Maintenance respiration per ground area at Stomate time 
                                                                         !! step  
                                                                         !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex
!$OMP THREADPRIVATE(resp_maint_d)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: resp_growth_d        !! Growth respiration per ground area 
                                                                         !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex
!$OMP THREADPRIVATE(resp_growth_d)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: resp_hetero_d        !! Heterotrophic respiration per ground area 
                                                                         !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex
!$OMP THREADPRIVATE(resp_hetero_d)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: resp_hetero_litter_d !! Heterotrophic respiration from litter per ground area 
                                                                         !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex
!$OMP THREADPRIVATE(resp_hetero_litter_d)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: resp_hetero_soil_d   !! Heterotrophic respiration from soil per ground area 
                                                                         !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex
!$OMP THREADPRIVATE(resp_hetero_soil_d)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: resp_hetero_radia    !! Heterothrophic respiration pe\r ground area at Sechiba
                                                                         !! time step 
                                                                         !! @tex $(gC m^{-2} dt_sechiba^{-1})$ @endtex 
!$OMP THREADPRIVATE(resp_hetero_radia)
!---
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:)   :: turnover_time       !! Turnover time of grasses 
                                                                         !! @tex $(dt_stomate^{-1})$ @endtex 
!$OMP THREADPRIVATE(turnover_time)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:,:) :: turnover_daily      !! Senescence-driven turnover (better: mortality) of 
                                                                         !! leaves and roots  
                                                                         !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex
!$OMP THREADPRIVATE(turnover_daily)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:,:) :: turnover_resid      !! The turnover left from turnover_daily at any given time step  
                                                                         !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex
!$OMP THREADPRIVATE(turnover_resid)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:,:) :: turnover_littercalc !! Senescence-driven turnover (better: mortality) of 
                                                                         !! leaves and roots at Sechiba time step 
                                                                         !! @tex $(gC m^{-2} dt_sechiba^{-1})$ @endtex 
!$OMP THREADPRIVATE(turnover_littercalc)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:,:) :: turnover_longterm   !! "Long term" (default 3 years) senescence-driven 
                                                                         !! turnover (better: mortality) of leaves and roots 
                                                                         !! @tex $(gC m^{-2} year^{-1})$ @endtex
!$OMP THREADPRIVATE(turnover_longterm)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:,:) :: bm_to_litter        !! Background (not senescence-driven) mortality of biomass
                                                                         !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex
!$OMP THREADPRIVATE(bm_to_litter)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:,:) :: bm_to_litter_resid  !! Left over bm_to_litter at any specific time step
                                                                         !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex
!$OMP THREADPRIVATE(bm_to_litter_resid)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:,:) :: tree_bm_to_litter   !! Background (not senescence-driven) mortality of biomass
                                                                         !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex
!$OMP THREADPRIVATE(tree_bm_to_litter)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:,:) :: tree_bm_to_litter_resid !! Left over tree_bm_to_litter at any specific time step
                                                                         !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex
!$OMP THREADPRIVATE(tree_bm_to_litter_resid)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:,:) :: bm_to_littercalc    !! conversion of biomass to litter per ground area at 
                                                                         !! Sechiba time step 
                                                                         !! @tex $(gC m^{-2} dt_sechiba^{-1})$ @endtex 
!$OMP THREADPRIVATE(bm_to_littercalc)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:,:) :: tree_bm_to_littercalc    !! conversion of biomass to litter per ground area at 
                                                                         !! Sechiba time step 
                                                                         !! @tex $(gC m^{-2} dt_sechiba^{-1})$ @endtex 
!$OMP THREADPRIVATE(tree_bm_to_littercalc)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:)  :: dead_leaves          !! Metabolic and structural pools of dead leaves on ground
                                                                         !! per PFT @tex $(gC m^{-2})$ @endtex 
!$OMP THREADPRIVATE(dead_leaves)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:,:,:):: litter             !! Above and below ground metabolic and structural litter 
                                                                         !! per ground area 
                                                                         !! @tex $(gC m^{-2})$ @endtex 
!$OMP THREADPRIVATE(litter)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:,:,:):: litterfuel         !! Aboveground deadfuel consisting of metabolic and structural 
                                                                         !! litter, classified into 4 categories depending on the scale
                                                                         !! of number hours to reach equilibrium with the ambient moisture
                                                                         !! per ground area
                                                                         !! @tex$(gC m^{-2})$ @endtex
!$OMP THREADPRIVATE(litterfuel)

  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: firelitter           !! Total litter above the ground that could potentially 
                                                                         !! burn @tex $(gC m^{-2})$ @endtex 
!$OMP THREADPRIVATE(firelitter)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:,:):: carbon_input         !! Quantity of carbon going into carbon pools from litter
                                                                         !! decomposition per ground area  at Sechiba time step 
                                                                         !! @tex $(gC m^{-2} dt_sechiba^{-1})$ @endtex 
!$OMP THREADPRIVATE(carbon_input)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:,:):: nitrogen_input       !! Quantity of nitrogen going into nitrogen pools from litter 
                                                                         !! decomposition per ground area  at Sechiba time step  
                                                                         !! @tex $(gC m^{-2} dtradia^{-1})$ @endtex  
!$OMP THREADPRIVATE(nitrogen_input) 
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:,:)  :: som_input_daily    !! Daily quantity of carbon going into carbon pools from
                                                                         !! litter decomposition per ground area 
                                                                         !! @tex $(gC m^{-2} day^{-1})$ @endtex 
!$OMP THREADPRIVATE(som_input_daily)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:,:)  :: som                !! Soil organic matter  pools per ground area: active, slow, or 
                                                                         !! passive, @tex $(gC or N m^{-2})$ @endtex 
!$OMP THREADPRIVATE(som)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:,:)  :: burried_litter     !! Litter burried under non-biological land uses (gC or N m-2)
!$OMP THREADPRIVATE(burried_litter)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:)    :: burried_fresh_ltr  !! Fresh litter burried under non-biological land uses (gC or N m-2)
!$OMP THREADPRIVATE(burried_fresh_ltr)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:)    :: burried_fresh_som  !! Fresh som burried under non-biological land uses (gC or N m-2)
!$OMP THREADPRIVATE(burried_fresh_som)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)        :: burried_bact       !! Bacteria burried under non-biological land uses (gC m-2)
!$OMP THREADPRIVATE(burried_bact)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)      :: burried_min_nitro  !! Mineral nitrogen burried under non-biological land uses (gC or N m-2)
!$OMP THREADPRIVATE(burried_min_nitro)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:)    :: burried_som        !! Som burried under non-biological land uses (gC or N m-2)
!$OMP THREADPRIVATE(burried_som)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:)    :: burried_deepSOM_a  !! Som burried under non-biological land uses (gC or N m-2)
!$OMP THREADPRIVATE(burried_deepSOM_a)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:)    :: burried_deepSOM_s  !! Som burried under non-biological land uses (gC or N m-2)
!$OMP THREADPRIVATE(burried_deepSOM_s)  
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:)    :: burried_deepSOM_p  !! Som burried under non-biological land uses (gC or N m-2)
!$OMP THREADPRIVATE(burried_deepSOM_p)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:)  :: lignin_struc         !! Ratio Lignine/Carbon in structural litter for above and
                                                                         !! below ground compartments (unitless)
!$OMP THREADPRIVATE(lignin_struc)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:)  :: lignin_wood          !! Ratio Lignine/Carbon in woody litter for above and
                                                                         !! below ground compartments (unitless)	
!$OMP THREADPRIVATE(lignin_wood)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:)  :: lignin_snag          !! Ratio Lignine/Carbon in snag litter for above and
                                                                         !! below ground compartments (unitless)	
!$OMP THREADPRIVATE(lignin_snag)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: lm_lastyearmax       !! Last year's maximum leaf mass per ground area for each
                                                                         !! PFT @tex $(gC m^{-2})$ @endtex  
!$OMP THREADPRIVATE(lm_lastyearmax)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: lm_thisyearmax       !! This year's maximum leaf mass per ground area for each
                                                                         !! PFT @tex $(gC m^{-2})$ @endtex  
!$OMP THREADPRIVATE(lm_thisyearmax)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: maxfpc_lastyear      !! Last year's maximum fpc for each natural PFT, on ground
                                                                         !! [??CHECK] fpc but this ones look ok (computed in 
                                                                         !! season, used in light)?? 
!$OMP THREADPRIVATE(maxfpc_lastyear)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: maxfpc_thisyear      !! This year's maximum fpc for each PFT, on ground (see 
                                                                         !! stomate_season), [??CHECK] fpc but this ones look ok 
                                                                         !! (computed in season, used in light)??
!$OMP THREADPRIVATE(maxfpc_thisyear)
!---
  REAL(r_std), ALLOCATABLE,SAVE,DIMENSION(:,:,:)  :: leaf_age            !! Age of different leaf classes (days)
!$OMP THREADPRIVATE(leaf_age)
  REAL(r_std), ALLOCATABLE,SAVE,DIMENSION(:,:,:)  :: leaf_frac           !! PFT fraction of leaf mass in leaf age class (0-1, 
                                                                         !! unitless) 
!$OMP THREADPRIVATE(leaf_frac)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: when_growthinit      !! Days since beginning of growing season (days)
!$OMP THREADPRIVATE(when_growthinit)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: herbivores           !! Time constant of probability of a leaf to be eaten by a
                                                                         !! herbivore (days)
!$OMP THREADPRIVATE(herbivores)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: RIP_time             !! How much time ago was the PFT eliminated for the last 
                                                                         !! time (year)
!$OMP THREADPRIVATE(RIP_time)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: time_hum_min         !! Time elapsed since strongest moisture limitation (days) 
!$OMP THREADPRIVATE(time_hum_min)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: drain_daily          !! daily fraction of water lost from the soil column by leaching (-)  
!$OMP THREADPRIVATE(drain_daily)

 REAL(r_std), ALLOCATABLE,SAVE,DIMENSION(:,:)    :: cn_leaf_min_season   !! Seasonal min CN ratio of leaves 
!$OMP THREADPRIVATE(cn_leaf_min_season)
 REAL(r_std), ALLOCATABLE,SAVE,DIMENSION(:,:)    :: nstress_season       !! N-related seasonal stress (used for allocation) 
!$OMP THREADPRIVATE(nstress_season)
 REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:)   :: soil_n_min           !! mineral nitrogen in the soil (gN/m**2)   
                                                                         !! (first index=kjpindex, second index=nvm, third index=nnspec)   
!$OMP THREADPRIVATE(soil_n_min)
 REAL(r_std), ALLOCATABLE,SAVE,DIMENSION(:,:)    :: p_O2                 !! partial pressure of oxigen in the soil (hPa)(first index=kjpindex, second index=nvm)
!$OMP THREADPRIVATE(p_O2)                      
 REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)  :: bact                 !! denitrifier biomass (gC/m**2)
                                                                         !! (first index=npts, second index=nvm)
!$OMP THREADPRIVATE(bact)   
!---
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: co2_fire             !! Carbon emitted to the atmosphere by burning living 
                                                                         !! and dead biomass 
                                                                         !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex
!$OMP THREADPRIVATE(co2_fire)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:,:)  :: emissions_fire     !! Carbon and other elements emitted to the atmosphere by burning living 
                                                                         !! and dead biomass
!$OMP THREADPRIVATE(emissions_fire)
 
!!$  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: co2_to_bm_dgvm       !! Psuedo-photosynthesis,C used to provide seedlings with
!!$                                                                         !! an initial biomass, arbitrarily removed from the 
!!$                                                                         !! atmosphere  
!!$                                                                         !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex 
!!$!$OMP THREADPRIVATE(co2_to_bm_dgvm)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)  :: atm_to_immob           !! N taken from the atmosphere to support immobilisation
                                                                         !! @tex $(gN m^{-2} dt_stomate^{-1})$ @endtex 
!$OMP THREADPRIVATE(atm_to_immob)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:,:,:)  :: prod_s           !! Wood products remaining in the 1 year-turnover pool 
!$OMP THREADPRIVATE(prod_s)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:,:,:)  :: prod_m           !! Wood products remaining in the 10 year-turnover pool 
                                                                         !! after the annual release for each compartment 
                                                                         !! @tex $(gC m^{-2})$ @endtex    
                                                                         !! (1:11 input from year of land cover change),
                                                                         !! dimension(#pixels,1:11 years
!$OMP THREADPRIVATE(prod_m)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:,:,:)  :: prod_l           !! Wood products remaining in the 100 year-turnover pool
                                                                         !! after the annual release for each compartment
                                                                         !! @tex $(gC m^{-2})$ @endtex  
                                                                         !! (1:101 input from year of land cover change), 
                                                                         !! dimension(#pixels,1:101 years)
!$OMP THREADPRIVATE(prod_l)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:,:,:)  :: flux_s           !! Wood decomposition from the 1 year-turnover pool 
                                                                         !! compartments 
                                                                         !! @tex $(gC m^{-2} year^{-1})$ @endtex 
                                                                         !! dimension(#pixels,1:2)  
!$OMP THREADPRIVATE(flux_s)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:,:,:)  :: flux_m           !! Wood decomposition from the 10 year-turnover pool 
                                                                         !! compartments 
                                                                         !! @tex $(gC m^{-2} year^{-1})$ @endtex 
                                                                         !! dimension(#pixels,1:11)  
!$OMP THREADPRIVATE(flux_m)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:,:,:)  :: flux_l           !! Wood decomposition from the 100 year-turnover pool 
                                                                         !! compartments 
                                                                         !! @tex $(gC m^{-2} year^{-1})$ @endtex
                                                                         !! dimension(#pixels,1:101)
!$OMP THREADPRIVATE(flux_l)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:,:)    :: flux_prod_s      !! Release during first year following land cover change 
                                                                         !! (paper, burned, etc...) 
                                                                         !! @tex $(gC m^{-2} year^{-1})$ @endtex  
!$OMP THREADPRIVATE(flux_prod_s)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:,:)    :: flux_prod_m      !! Total annual release from the 10 year-turnover pool
                                                                         !! sum of flux_m  
                                                                         !! @tex $(gC m^{-2} year^{-1})$ @endtex
!$OMP THREADPRIVATE(flux_prod_m)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:,:,:)    :: flux_prod_l      !! Total annual release from the 100 year-turnover pool 
                                                                         !! sum of flux_l 
                                                                         !! @tex $(gC m^{-2} year^{-1})$ @endtex
!$OMP THREADPRIVATE(flux_prod_l)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: fco2_lu              !! CO2 flux between atmosphere and biosphere from land-use 
                                                                         !! (without forest management)
                                                                         !! @tex $(gC m^{-2} one_day^{-1})$ @endtex
!$OMP THREADPRIVATE(fco2_lu)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: fco2_wh              !! CO2 Flux to Atmosphere from Wood Harvesting (positive from atm to land)
                                                                         !! @tex $(gC m^{-2} one_day^{-1})$ @endtex
!$OMP THREADPRIVATE(fco2_wh)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: fco2_ha              !! CO2 Flux to Atmosphere from Crop Harvesting (positive from atm to land)
                                                                         !! @tex $(gC m^{-2} one_day^{-1})$ @endtex
!$OMP THREADPRIVATE(fco2_ha)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: fDeforestToProduct   !! Deforested biomass into product pool due to anthropogenic
                                                                         !! land use change

!$OMP THREADPRIVATE(fDeforestToProduct)   
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: fLulccResidue        !! Carbon mass flux into soil and litter due to anthropogenic land use or land cover change                                                                          
!$OMP THREADPRIVATE(fLulccResidue)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)    :: fHarvestToProduct    !! Deforested biomass into product pool due to anthropogenic
                                                                         !! land use
!$OMP THREADPRIVATE(fHarvestToProduct)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:):: woodharvestpft       !! New year wood harvest per  PFT
!$OMP THREADPRIVATE(woodharvestpft)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)      :: carb_mass_total      !! Total on-site and off-site C pool 
                                                                         !! @tex $(gC m^{-2})$ @endtex                        
!$OMP THREADPRIVATE(carb_mass_total)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:)   :: deepSOM_a      !! deep active SOM profile (g/m**3)
!$OMP THREADPRIVATE(deepSOM_a)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:)   :: deepSOM_s      !! deep slow SOM profile (g/m**3)
!$OMP THREADPRIVATE(deepSOM_s)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:)   :: deepSOM_p      !! deep passive SOM profile (g/m**3)
!$OMP THREADPRIVATE(deepSOM_p)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:,:)    :: altmax_ind    
!$OMP THREADPRIVATE(altmax_ind)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)       :: altmax_lastyear
!$OMP THREADPRIVATE(altmax_lastyear)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:,:)    :: altmax_ind_lastyear
!$OMP THREADPRIVATE(altmax_ind_lastyear)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: O2_soil          !! deep oxygen
!$OMP THREADPRIVATE(O2_soil)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: CH4_soil         !! deep methane
!$OMP THREADPRIVATE(CH4_soil)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: O2_snow          !! snow oxygen
!$OMP THREADPRIVATE(O2_snow)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: CH4_snow         !! snow methane
!$OMP THREADPRIVATE(CH4_snow)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: tdeep_daily      !! daily t profile (K)
!$OMP THREADPRIVATE(tdeep_daily)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: hsdeep_daily     !! daily humidity profile (unitless)
!$OMP THREADPRIVATE(hsdeep_daily)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)       :: temp_sol_daily   !! daily soil surface temp (K)
!$OMP THREADPRIVATE(temp_sol_daily)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)       :: pb_pa_daily      !! daily surface pressure [Pa]
!$OMP THREADPRIVATE(pb_pa_daily)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)       :: snow_daily       !! daily snow mass
!$OMP THREADPRIVATE(snow_daily)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: fbact            !! turnover constant for soil carbon discretization (day)
!$OMP THREADPRIVATE(fbact)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: decomp_rate      !! decomposition constant for soil carbon discretization (day-1)
!$OMP THREADPRIVATE(decomp_rate)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: decomp_rate_daily!! decomposition constant for soil carbon discretization (day)
!$OMP THREADPRIVATE(decomp_rate_daily)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)     :: fixed_cryoturbation_depth  !! depth to hold cryoturbation to for fixed runs
!$OMP THREADPRIVATE(fixed_cryoturbation_depth)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)     :: snowdz_daily       !! daily snow depth profile [m]
!$OMP THREADPRIVATE(snowdz_daily)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)     :: snowrho_daily      !! daily snow density profile (Kg/m^3)
!$OMP THREADPRIVATE(snowrho_daily)

  ! Below are the variables needed to be written to the soil carbon discretization spinup file
  REAL(r_std),DIMENSION(:,:,:,:,:),ALLOCATABLE  :: som_input_2pfcforcing   !! quantity of carbon going into carbon pools from
                                                                           !! litter decomposition per ground area 
                                                                           !! @tex $(gC m^{-2} day^{-1})$ @endtex for forcesoil
!$OMP THREADPRIVATE(som_input_2pfcforcing)
  REAL(r_std),DIMENSION(:,:),ALLOCATABLE      :: pb_2pfcforcing            !! surface pressure [Pa] for forcesoil
!$OMP THREADPRIVATE(pb_2pfcforcing)
  REAL(r_std),DIMENSION(:,:),ALLOCATABLE      :: snow_2pfcforcing          !! snow mass for forcesoil
!$OMP THREADPRIVATE(snow_2pfcforcing)
  REAL(r_std),DIMENSION(:,:,:,:),ALLOCATABLE  :: tprof_2pfcforcing         !! Soil temperature (K) for forcesoil
!$OMP THREADPRIVATE(tprof_2pfcforcing)
  REAL(r_std),DIMENSION(:,:,:,:),ALLOCATABLE  :: fbact_2pfcforcing         !! turnover constant for forcesoil (day)
!$OMP THREADPRIVATE(fbact_2pfcforcing)
  REAL(r_std),DIMENSION(:,:,:,:),ALLOCATABLE  :: hslong_2pfcforcing        !! Soil humiditity (-) for forcesoil
!$OMP THREADPRIVATE(hslong_2pfcforcing)
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE    :: veget_max_2pfcforcing     !! Vegetation coverage taking into account non-biological
                                                                           !! coverage (unitless) for forcesoil
!$OMP THREADPRIVATE(veget_max_2pfcforcing)
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE    :: rprof_2pfcforcing         !! Coefficient of the exponential functions that 
                                                                           !! relates root density to soil depth (unitless), for forcesoil  
!$OMP THREADPRIVATE(rprof_2pfcforcing)
  REAL(r_std),DIMENSION(:,:),ALLOCATABLE      :: tsurf_2pfcforcing         !! Surface temperatures (K), for forcesoil 
!$OMP THREADPRIVATE(tsurf_2pfcforcing)
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE    :: snowdz_2pfcforcing        !! Snow depth profile [m], for forcesoil
!$OMP THREADPRIVATE(snowdz_2pfcforcing)
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE    :: snowrho_2pfcforcing       !! Snow density profile (Kg/m^3), for forcesoil
!$OMP THREADPRIVATE(snowrho_2pfcforcing)
  REAL(r_std),DIMENSION(:,:,:,:),ALLOCATABLE  :: CN_target_2pfcforcing     !! 
!$OMP THREADPRIVATE(CN_target_2pfcforcing)
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE  :: n_mineralisation_2pfcforcing     !! 
!$OMP THREADPRIVATE(n_mineralisation_2pfcforcing)

!---
  REAL(r_std), SAVE                              :: tau_longterm
!$OMP THREADPRIVATE(tau_longterm)
  REAL(r_std),SAVE                               :: dt_days=zero         !! Time step of STOMATE (days) 
!$OMP THREADPRIVATE(dt_days)
  INTEGER(i_std),SAVE                            :: days_since_beg=0     !! Number of full days done since the start of the simulation
!$OMP THREADPRIVATE(days_since_beg)
  INTEGER(i_std),ALLOCATABLE,SAVE,DIMENSION(:)   :: nforce               !! Number of states calculated for the soil forcing 
                                                                         !! variables (unitless), dimension(::nparan*::nbyear) both 
                                                                         !! given in the run definition file    
!$OMP THREADPRIVATE(nforce)
  INTEGER(i_std), SAVE                           :: spinup_period        !! Period of years used to calculate the resolution of the system for spinup analytic. 
                                                                         !! This period correspond in most cases to the period of years of forcing data used
!$OMP THREADPRIVATE(spinup_period)
  INTEGER,PARAMETER                              :: r_typ = nf90_real4   !! Specify data format (server dependent)
!---
  LOGICAL, SAVE                                  :: do_slow=.FALSE.      !! Flag that determines whether stomate_accu calculates
                                                                         !! the sum(do_slow=.FALSE.) or the mean 
                                                                         !! (do_slow=.TRUE.)
!$OMP THREADPRIVATE(do_slow)
  LOGICAL, SAVE                                  :: l_first_stomate = .TRUE.!! Is this the first call of stomate?
!$OMP THREADPRIVATE(l_first_stomate)
!--- 
  INTEGER(i_std), SAVE                               :: global_years        !! Global counter of years (year)
!$OMP THREADPRIVATE(global_years)
  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:)           :: ok_equilibrium      !! Logical array marking the points where the resolution is ok 
                                                                            !! (true/false)
!$OMP THREADPRIVATE(ok_equilibrium)
  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:)           :: carbon_eq           !! Logical array to mark the carbon pools at equilibrium ? 
                                                                            !! If true, the job stops. (true/false)
!$OMP THREADPRIVATE(carbon_eq)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)     :: nbp_accu_flux       !! Accumulated Net Biospheric Production over the year (gC.m^2 )
!$OMP THREADPRIVATE(nbp_accu_flux)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)     :: nbp_pool_start      !! Biomass pool as calculated from the 
                                                                            !! previous time step (gC/N m-2)
!$OMP THREADPRIVATE(nbp_pool_start)
  REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:,:)       :: matrixA             !! matrix containing the fluxes between the carbon pools
                                                                            !! per sechiba time step 
                                                                            !! @tex $(gC.m^2.day^{-1})$ @endtex
!$OMP THREADPRIVATE(matrixA)
  REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:)         :: vectorB             !! vector containing the litter increase per sechiba time step
                                                                            !! @tex $(gC m^{-2})$ @endtex
!$OMP THREADPRIVATE(vectorB)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) :: matrixV             !! matrix containing the accumulated values of matrixA 
!$OMP THREADPRIVATE(matrixV)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: vectorU             !! matrix containing the accumulated values of vectorB
!$OMP THREADPRIVATE(vectorU)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) :: matrixW             !! matrix containing the opposite of matrixA
!$OMP THREADPRIVATE(matrixW)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: previous_stock      !! Array containing the carbon stock calculated by the analytical
                                                                            !! method in the previous resolution
!$OMP THREADPRIVATE(previous_stock)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: current_stock       !! Array containing the carbon stock calculated by the analytical
                                                                            !! method in the current resolution 
!$OMP THREADPRIVATE(current_stock)
  REAL(r_std), SAVE                                  :: eps_carbon          !! Stopping criterion for carbon pools (unitless,0-1)
!$OMP THREADPRIVATE(eps_carbon)
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)           :: sigma            !! Threshold for indivudal tree growth (m, trees whose
                                                                         !! circumference is smaller than sigma don't grow much)
!$OMP THREADPRIVATE(sigma)

  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:,:)  :: age_stand        !! Age of stand (years)
!$OMP THREADPRIVATE(age_stand)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)     :: age_stand_bm     !! Biomass-weighted mean stand age, conserved through the diameter age-class conveyor (years) - STAND_AGE, Marie 2026
!$OMP THREADPRIVATE(age_stand_bm)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)     :: age_stand_area   !! AREA-weighted mean stand age, strict mirror of age_stand_bm to expose the weighting-basis gap (years) - STAND_AGE
!$OMP THREADPRIVATE(age_stand_area)
  
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:,:)  :: rotation_n       !! Rotation number (number of rotation since pft is managed)
!$OMP THREADPRIVATE(rotation_n)
  
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:,:)  :: last_cut         !! Years since last thinning (years)
!$OMP THREADPRIVATE(last_cut)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  :: CN_som_litter_longterm !! Longterm CN ratio of litter and som pools (gC/gN)
!$OMP THREADPRIVATE(CN_som_litter_longterm)
  REAL(r_std), SAVE                                 :: tau_CN_longterm      !! Counter used for calculating the longterm CN ratio of SOM and litter pools (seconds)    
!$OMP THREADPRIVATE(tau_CN_longterm)

  ! Functional Allocation

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)     :: KF               !! Scaling factor to convert sapwood mass
                                                                         !! into leaf mass (m). The initial value is calculated
                                                                         !! in prescribe and updated during allocation
!$OMP THREADPRIVATE(KF)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)     :: k_latosa_adapt   !! Leaf to sapwood area adapted for waterstress. 
                                                                         !! Adaptation takes place at the end of the year
                                                                         !! (m)
!$OMP THREADPRIVATE(k_latosa_adapt)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:,:) :: harvest_pool_acc !! the accumulative value of harvest_pool throughout everyday.
!$OMP THREADPRIVATE(harvest_pool_acc)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)     :: harvest_type     !! Type of management that resulted
                                                                         !! in the harvest (unitless)
!$OMP THREADPRIVATE(harvest_type)  

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)     :: harvest_cut      !! Type of cutting that was used for the harvest
                                                                         !! (unitless)
!$OMP THREADPRIVATE(harvest_cut)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: harvest_area_acc !! Harvested area (m^{2})
!$OMP THREADPRIVATE(harvest_area_acc)  

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: gap_area_save    !! Total gap area created by more than 30% basal area loss 
                                                                         !!  in the last 5 years (m^{2})
!$OMP THREADPRIVATE(gap_area_save)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)     :: total_ba_init    !! Total basal area saved at the first day of the year per PFT 
                                                                         !! (m^{2}/m^{2})
!$OMP THREADPRIVATE(total_ba_init)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)       :: harvest_pool_bound       !! The boundaries of the diameter classes
                                                                                 !! in the wood harvest pools
                                                                                 !! @tex $(m)$ @endtex 
!$OMP THREADPRIVATE(harvest_pool_bound)

!! START : stomate_pest module (bark beetle outbreak)

INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: beetle_diapause           !! A beetle phenology 
                                                                                 !! stage which trigger 
                                                                                 !! the reproduction (binary) 
!$OMP THREADPRIVATE(beetle_diapause)
! Variables related to bark beetle module
REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  :: wood_leftover_legacy        !! woody litter used in the 
                                                                                 !! calculation of the windthrow                                                                                  !! susceptibility (gC.m-2)
!$OMP THREADPRIVATE(wood_leftover_legacy)
REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  :: season_drought_legacy       !! poxy of the tree healthiness
                                                                                 !! use in the calculation of 
                                                                                 !! the beetle outbreak 
                                                                                 !! susceptibility (unitless)  
!$OMP THREADPRIVATE(season_drought_legacy)
REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:)    :: i_beetles_generation      !! number of beetle generations
                                                                                 !! per year 
!$OMP THREADPRIVATE(i_beetles_generation)
REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)      :: P_beetles_attacked_legacy !! use to trigger the 
                                                                                 !! epidemic flag (unitless)
!$OMP THREADPRIVATE(P_beetles_attacked_legacy)
REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:)    :: B_beetles_kill_legacy      !! Stored beetle damage use 
                                                                                 !!in the estimation of the 
                                                                                 !!woody leftover (gC.m-2)
!$OMP THREADPRIVATE(B_beetles_kill_legacy)
REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  :: i_beetles_activity_legacy   !! Proxy of the size of the 
                                                                                 !! beetle population (unitless)
!$OMP THREADPRIVATE(i_beetles_activity_legacy)
REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: sumTeff                      !! Sum of air temperature used                                                                                   
                                                                                 !! in the calculation of the 
                                                                                 !! beetle diapause (°C)
!$OMP THREADPRIVATE(sumTeff)
REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: woody_litter_to_use          !! Yearly sum of woody litter pool to use for 
                                                                                 !! wood_left_over in the pest module (gC m^{-2})
!$OMP THREADPRIVATE(woody_litter_to_use)
 
!! END : stomate_pest module

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)     :: mai              !! The mean annual increment used in
                                                                         !! forestry.  It is the average change
                                                                         !! in the wood volume of of the trunk
                                                                         !! over the lifetime of the forest.
                                                                         !! @tex $(m**3 / m**2 / year)$ @endtex 
!$OMP THREADPRIVATE(mai)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)     :: pai              !! The period annual increment used in
                                                                         !! forestry.  It is the average change
                                                                         !! in the wood volume of of the trunk
                                                                         !! over the past n_pai years of the forest,
                                                                         !! where n_pai is defined in constants.f90.
                                                                         !! @tex $(m**3 / m**2 / year)$ @endtex 
!$OMP THREADPRIVATE(pai)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)     :: previous_wood_volume !! The volume of the tree trunks
                                                                         !! in a stand for the previous year.
                                                                         !! @tex $(m**3 / m**2 )$ @endtex 
!$OMP THREADPRIVATE(previous_wood_volume)

  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:,:)  :: mai_count        !! The number of times we've
                                                                         !! calculated the volume increment
                                                                         !! for a stand
!$OMP THREADPRIVATE(mai_count)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)     :: coppice_dens     !! The density of a coppice at the first
                                                                         !! cutting.
                                                                         !! @tex $( 1 / m**2 )$ @endtex 
!$OMP THREADPRIVATE(coppice_dens)
  REAL(r_std), ALLOCATABLE,SAVE,DIMENSION(:,:)       :: rue_longterm     !! Longterm radiation use efficiency (??units??)
!$OMP THREADPRIVATE(rue_longterm)
  REAL(r_std), ALLOCATABLE,SAVE,DIMENSION(:,:)       :: leaf_age_crit    !! critical leaf age (days)
!$OMP THREADPRIVATE(leaf_age_crit)
  REAL(r_std), ALLOCATABLE,SAVE,DIMENSION(:,:)       :: leaf_classes     !! Width of each leaf age class (days)
!$OMP THREADPRIVATE(leaf_classes)
  REAL(r_std), ALLOCATABLE,SAVE,DIMENSION(:,:,:,:,:) :: bm_sapl_2D 
!$OMP THREADPRIVATE(bm_sapl_2D)
  REAL(r_std),SAVE                                   :: dt_forcesoil        !! Time step of soil forcing file (days)
!$OMP THREADPRIVATE(dt_forcesoil)
  INTEGER(i_std),PARAMETER                           :: nparanmax=366       !! Maximum number of time steps per year for forcesoil
  INTEGER(i_std),SAVE                                :: nparan              !! Number of time steps per year for forcesoil read from run definition (unitless) 
!$OMP THREADPRIVATE(nparan)
  INTEGER(i_std),SAVE                                :: nbyear=1            !! Number of years saved for forcesoil (unitless) 
!$OMP THREADPRIVATE(nbyear)
  INTEGER(i_std),SAVE                                :: iatt                !! Time step of forcing of soil processes (iatt = 1 to ::nparan*::nbyear) 
!$OMP THREADPRIVATE(iatt)
  INTEGER(i_std),SAVE                                :: iatt_old=1          !! Previous ::iatt
!$OMP THREADPRIVATE(iatt_old)
  CHARACTER(LEN=100), SAVE                           :: Cforcing_discretization_name !! Name of forcing file 2
!$OMP THREADPRIVATE(Cforcing_discretization_name)
  INTEGER(i_std), SAVE                               :: frozen_respiration_func  !! Method for soil decomposition function
!$OMP THREADPRIVATE(frozen_respiration_func)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:,:)  :: forest_managed     !! forest management flag (is the forest managed?) (0-4,unitless)
!$OMP THREADPRIVATE(forest_managed)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:,:)  :: species_change_map   !! A map which gives the PFT number that each
                                                                         !! PFT will be replanted as in case of a clearcut.
                                                                         !! (1-nvm,unitless)
!$OMP THREADPRIVATE(species_change_map)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:,:)  :: fm_change_map    !! A map which gives the desired FM strategy when
                                                                         !! the PFT will be replanted after a clearcut.
                                                                         !! (1-nvm,unitless)
!$OMP THREADPRIVATE(fm_change_map)
  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:,:)         :: lpft_replant     !! Indicates if this PFT has either died this year
                                                                         !! or been clearcut/coppiced.  If it has, it is not
                                                                         !! replanted until the end of the year.
!$OMP THREADPRIVATE(lpft_replant)
REAL(r_std), ALLOCATABLE, SAVE,DIMENSION(:,:)      :: wstress_season   !! Water stress factor, based on hum_rel_daily
                                                                         !! (unitless, 0-1)
!$OMP THREADPRIVATE(wstress_season)

  REAL(r_std), ALLOCATABLE, SAVE,DIMENSION(:,:)     :: wstress_month    !! Water stress factor, based on hum_rel_daily
                                                                         !! (unitless, 0-1)
!$OMP THREADPRIVATE(wstress_month)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)    :: light_tran_to_floor_season  !! Mean seasonal fraction of light transmitted 
                                                                                   !! to canopy levels
!$OMP THREADPRIVATE(light_tran_to_floor_season)
  INTEGER(i_std), SAVE                              :: printlev_loc                !! Local level of text output for current module
!$OMP THREADPRIVATE(printlev_loc)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)    :: sugar_load                  !! Relative sugar loading of the labile pool (unitless) 
!$OMP THREADPRIVATE(sugar_load)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)    :: grow_season_len             !! growing season length in days for deciduous PFTs.
!$OMP THREADPRIVATE(grow_season_len) 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)    :: doy_start_gs                !! growing season starting day of year (DOY) for 
                                                                                   !! deciduous PFTs.
!$OMP THREADPRIVATE(doy_start_gs)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)    :: doy_end_gs                  !! growing season end day of year (DOY) for 
                                                                                   !! deciduous PFTs.
!$OMP THREADPRIVATE(doy_end_gs)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)    :: mean_start_gs               !! mean growing season starting day for 
                                                                                   !! deciduous PFTs.
!$OMP THREADPRIVATE(mean_start_gs)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)      :: est_co2flux_d               !! Accumulated estimate of co2flux within a single time step. The value 
                                                                                   !! is not exact because several fluxes are missing (gC m-2 s-1)
!$OMP THREADPRIVATE(est_co2flux_d)

!spitfire
  REAL(r_std),ALLOCATABLE,SAVE, DIMENSION(:)        :: ni_acc                      !! Nesterov index (square of degree Celcius)
!$OMP THREADPRIVATE(ni_acc)

  REAL(r_std),ALLOCATABLE,SAVE, DIMENSION(:,:)      :: lightning                   !! Number of daily lightning (flashes km^{-2} day^{-1})
!$OMP THREADPRIVATE(lightning)

  REAL(r_std),ALLOCATABLE,SAVE, DIMENSION(:)        :: pop_dens                    !! Population density
!$OMP THREADPRIVATE(pop_dens)

  REAL(r_std),ALLOCATABLE,SAVE, DIMENSION(:)        :: a_nd                        !! Nesterov index (square of degree Celcius)
!$OMP THREADPRIVATE(a_nd)
  REAL(r_std),ALLOCATABLE,SAVE, DIMENSION(:)        :: edge_length                 !! Edge length from road observations (m) (urban road should be excluded)
!$OMP THREADPRIVATE(edge_length)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)         :: daily_wspeed_fire           !! Daily wind speed at 2 meter (ms-1)
!$OMP THREADPRIVATE(daily_wspeed_fire)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)       :: mbc_daily                   !! Cumulated daily mass balance closure g/m2/day
!$OMP THREADPRIVATE(mbc_daily)

  ! Parameters (from Chunjing) useful for the soil carbon without discretisation. To be removed later.
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)         :: height_acro
!$OMP THREADPRIVATE(height_acro)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)         :: height_cato
!$OMP THREADPRIVATE(height_cato)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)       :: carbon_acro
!$OMP THREADPRIVATE(carbon_acro)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)       :: carbon_cato
!$OMP THREADPRIVATE(carbon_cato)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:,:)       :: resp_acro_oxic !!respiration of acrotelm( oxic )
!$OMP THREADPRIVATE(resp_acro_oxic)
  REAL(r_std),ALLOCATABLE,SAVE, DIMENSION(:,:)      :: resp_acro_anoxic !!respiration of acrotelm( anoxic )
!$OMP THREADPRIVATE(resp_acro_anoxic)
  REAL(r_std),ALLOCATABLE,SAVE, DIMENSION(:,:)      :: resp_cato !!respiration of catotelm
!$OMP THREADPRIVATE(resp_cato)
  REAL(r_std),ALLOCATABLE,SAVE, DIMENSION(:,:)      :: acro_to_cato
!$OMP THREADPRIVATE(acro_to_cato)
  REAL(r_std),ALLOCATABLE,SAVE, DIMENSION(:,:)      :: litter_to_acro
!$OMP THREADPRIVATE(litter_to_acro)
  REAL(r_std),ALLOCATABLE,SAVE, DIMENSION(:)        :: tcarbon_acro
!$OMP THREADPRIVATE(tcarbon_acro)
  REAL(r_std),ALLOCATABLE,SAVE, DIMENSION(:)        :: tcarbon_cato
!$OMP THREADPRIVATE(tcarbon_cato) 
  REAL(r_std),ALLOCATABLE,DIMENSION(:,:)            :: wtp_pt    ! above surface water table for peatland
!$OMP THREADPRIVATE(wtp_pt)
  REAL(r_std),ALLOCATABLE,DIMENSION(:)              :: wtp_daily ! above surface water table for peatland, daily basis
!$OMP THREADPRIVATE(wtp_daily)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:):: deepSOM_peat ! soil carbon when peatland is activated
!$OMP THREADPRIVATE(deepSOM_peat)
  
PUBLIC  dt_days, days_since_beg, do_slow

CONTAINS
  

!! ================================================================================================================================
!! SUBROUTINE 	: stomate_initialize
!!
!>\BRIEF        Initialization routine for stomate module. 
!!
!! DESCRIPTION  : Initialization routine for stomate module. Read options from parameter file, allocate variables, read variables 
!!                from restart file and initialize variables if necessary. 
!!                
!! \n
!_ ================================================================================================================================

SUBROUTINE stomate_initialize &
        (kjit,           kjpij,             kjpindex,                        &
         rest_id_stom,   hist_id_stom,      hist_id_stom_IPCC,               &
         index,          lalo,              neighbours,   resolution,        &
         contfrac,       clay,              silt,                            &
         bulk,           temp_air,                                           &
         veget,          veget_max,                                          &
         deadleaf_cover, assim_param,      circ_class_biomass, circ_class_n, &
         lai_per_level,  laieff_fit,       temp_growth,                      &
         som_total,      heat_Zimov,       altmax, depth_organic_soil,       & 
         cn_leaf_init_2D)

    IMPLICIT NONE
    !! 0. Variable and parameter declaration
    !! 0.1 Input variables
    INTEGER(i_std),INTENT(in)                       :: kjit              !! Time step number (unitless)
    INTEGER(i_std),INTENT(in)                       :: kjpij             !! Total size of the un-compressed grid (unitless)
    INTEGER(i_std),INTENT(in)                       :: kjpindex          !! Domain size - terrestrial pixels only (unitless)
    INTEGER(i_std),INTENT(in)                       :: rest_id_stom      !! STOMATE's _Restart_ file identifier (unitless)
    INTEGER(i_std),INTENT(in)                       :: hist_id_stom      !! STOMATE's _history_ file identifier (unitless)
    INTEGER(i_std),INTENT(in)                       :: hist_id_stom_IPCC !! STOMATE's IPCC _history_ file identifier(unitless) 
    INTEGER(i_std),DIMENSION(:),INTENT(in)          :: index             !! The indices of the terrestrial pixels only (unitless) 
    REAL(r_std),DIMENSION(:,:),INTENT(in)           :: lalo              !! Geographical coordinates (latitude,longitude) for pixels (degrees) 
    INTEGER(i_std),DIMENSION(:,:),INTENT(in)        :: neighbours        !! Neighoring grid points if land for the DGVM (unitless) 
    REAL(r_std),DIMENSION(:,:),INTENT(in)           :: resolution        !! Size in x an y of the grid (m) - surface area of the gridbox 
    REAL(r_std),DIMENSION (:), INTENT (in)          :: contfrac          !! Fraction of continent in the grid cell (unitless)
    REAL(r_std),DIMENSION(:),INTENT(in)             :: clay              !! Clay fraction of soil (0-1, unitless)
    REAL(r_std),DIMENSION(:),INTENT(in)             :: silt              !! Silt fraction of soil (0-1, unitless) 
    REAL(r_std),DIMENSION(:),INTENT(in)             :: bulk              !! Bulk density (kg/m**3) 
    REAL(r_std),DIMENSION(:),INTENT(in)             :: temp_air          !! Air temperature at first atmospheric model layer (K)
    REAL(r_std),DIMENSION(:,:),INTENT(in)           :: veget             !! Fraction of vegetation type including 
                                                                         !! non-biological fraction (unitless) 
    REAL(r_std),DIMENSION(:,:),INTENT(in)           :: veget_max         !! Maximum fraction of vegetation type including 
                                                                         !! non-biological fraction (unitless) 
    REAL(r_std),DIMENSION(:,:), INTENT(in)          :: cn_leaf_init_2D   !! initial leaf C/N ratio 

    !! 0.2 Output variables

    REAL(r_std),DIMENSION(:),INTENT(out)            :: deadleaf_cover    !! Fraction of soil covered by dead leaves (unitless)
    REAL(r_std),DIMENSION(:,:,:),INTENT(out)        :: assim_param       !! min+max+opt temperatures (K) & vmax for photosynthesis  
                                                                         !! @tex $(\mu mol m^{-2}s^{-1})$ @endtex  
    REAL(r_std),DIMENSION(:),INTENT(out)            :: temp_growth       !! Growth temperature (Â°C)  
                                                                         !! Is equal to t2m_month 
    REAL(r_std), DIMENSION(:,:,:), INTENT (out)     :: heat_Zimov        !! heating associated with decomposition [W/m**3 soil]
    REAL(r_std),DIMENSION(:,:), INTENT(out)         :: altmax            !! Maximul active layer thickness (m). Be careful, here active means non frozen.
                                                                         !! Not related with the active soil carbon pool.
    REAL(r_std), DIMENSION(:), INTENT (out)         :: depth_organic_soil!! Depth at which there is still organic matter (m)

    !! 0.3 Modified variables
    REAL(r_std),DIMENSION(:,:,:,:,:),INTENT(inout)  :: circ_class_biomass!! Biomass per circumference class @tex $(gC tree^{-1})$ @endtex
    REAL(r_std),DIMENSION(:,:,:),INTENT(inout)      :: circ_class_n      !! Number of trees within each circumference
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)    :: lai_per_level     !! This is the LAI per vertical level
                                                                         !! @tex $(m^{2} m^{-2})$ @endtex
    TYPE(laieff_type),DIMENSION(:,:,:),INTENT(inout):: laieff_fit        !! Fitted parameters for the effective LAI

    REAL(r_std), DIMENSION(:,:,:,:), INTENT (inout) :: som_total         !! total soil carbon for use in thermal calcs (g/m**3)
 
    !! 0.4 Local variables
    REAL(r_std)                                   :: dt_days_read         !! STOMATE time step read in restart file (days)
    INTEGER(i_std)                                :: l,k,ji,jv,i,j,ipts   !! indices
    INTEGER(i_std)                                :: ivm,icir             !! indices
    REAL(r_std),PARAMETER                         :: max_dt_days = 5.     !! Maximum STOMATE time step (days)
    REAL(r_std),DIMENSION(kjpindex,nvm)           :: gpp_daily_x          !! "Daily" gpp for teststomate  
                                                                          !! @tex $(??gC m^{-2} dt_stomate^{-1})$ @endtex 
    INTEGER(i_std)                                :: ier                  !! Check errors in netcdf call (unitless)
    INTEGER(i_std)                                :: max_totsize          !! Memory management - maximum memory size (Mb)
    INTEGER(i_std)                                :: totsize_1step        !! Memory management - memory required to store one 
                                                                          !! time step on one processor (Mb) 
    INTEGER(i_std)                                :: totsize_tmp          !! Memory management - memory required to store one 
                                                                          !! time step on all processors(Mb) 
    INTEGER(i_std)                                :: vid                  !! Variable identifer of netCDF (unitless)
    INTEGER(i_std)                                :: nneigh               !! Number of neighbouring pixels
    INTEGER(i_std)                                :: direct               !! 
    LOGICAL                                       :: l_error              !! error flag
    REAL(r_std)                                   :: temp_total           !! Used for renormalizing
    CHARACTER(LEN=200)                            :: temp_str 
    CHARACTER(LEN=255)                            :: filename
    CHARACTER(LEN=255)                            :: field_name
    INTEGER(i_std)                                :: imi                  !! indice de classe d'intensite de gestion (1..nmiclass)
     INTEGER(i_std)                                :: iage_st              !! Indice de classe d'age d'un creneau PFT (init STAND_AGE derive)
 !================================================================================================================================
    
    !! 1. Initialize variable

    !! Initialize local printlev
      ! Guillaume M. -- get_printlev calls getin_p, a COLLECTIVE bcast. Called here,
      ! inside the compute loop under a local firstcall, it hangs the run as soon as the
      ! ranks stop reaching it together. Use printlev, read once at startup by every rank.
    printlev_loc=printlev
    !! Initialize module stomate_laieff
    CALL stomate_laieff_initialize()

    !! Initialize module sapiens_agriculture
    CALL sapiens_agriculture_initialize()

    !! Update flag
    l_first_stomate = .FALSE.
    
    !! 1.1 Store current time step in a common variable
    itime = kjit
    
!!$    !! 1.3.1 Set lai 
!!$    lai(:,ibare_sechiba) = zero
!!$    DO i = 1, kjpindex
!!$       DO j = 2,nvm
!!$          lai(i,j) = cc_to_lai(circ_class_biomass(i,j,:,ileaf,icarbon),circ_class_n(i,j,:),j)
!!$       ENDDO
!!$    ENDDO

    !! 1.4.0 Parameters for spinup
    !
    eps_carbon = 0.01
    !Config Key   = EPS_CARBON
    !Config Desc  = Allowed error on carbon stock
    !Config If    = SPINUP_ANALYTIC
    !Config Def   = 0.01
    !Config Help  = 
    !Config Units = [%]   
    CALL getin_p('EPS_CARBON',eps_carbon)       
    
    
    !Config Key   = SPINUP_PERIOD
    !Config Desc  = Period to calulcate equilibrium during spinup analytic
    !Config If    = SPINUP_ANALYTIC
    !Config Def   = -1
    !Config Help  = Period corresponds in most cases to the number of years of forcing data used in the spinup.
    !Config Units = [years]   
    spinup_period = -1
    CALL getin_p('SPINUP_PERIOD',spinup_period)       
    
    ! Check spinup_period values. 
    ! For periods uptil 6 years, to obtain equilibrium, a bigger period have to be used 
    ! and therefore spinup_period is adjusted to 10 years. 
    IF (spinup_analytic) THEN
       IF (spinup_period <= 0) THEN
          WRITE(numout,*) 'Error in parameter spinup_period. This parameter must be > 0 : spinup_period=',spinup_period
          CALL ipslerr_p (3,'stomate_initialize', &
               'Parameter spinup_period must be set to a positive integer.', &
               'Set this parameter to the number of years of forcing data used for the spinup.', &
               '')
       END IF
       IF (printlev >=1) WRITE(numout,*) 'Spinup analytic is activated using eps_carbon=',&
            eps_carbon, ' and spinup_period=',spinup_period
    END IF
    

    !! 1.4.1 Allocate memory for all variables in stomate
    ! Allocate memory for all variables in stomate, build new index
    ! tables accounting for the PFTs, read and check flags and set file
    ! identifier for restart and history files.
    CALL stomate_init (kjpij, kjpindex, index, lalo, &
         rest_id_stom, hist_id_stom, hist_id_stom_IPCC)
    
    !! 1.4.2 Initialization some parameters
    ! Note that lots of this is now taken care off in 
    ! constantes_mtc.f90 and pft_parameters. f90. 
    CALL data
    
    !! 1.4.3 Initial conditions
    
    !! 1.4.3.1 Read initial values for STOMATE's variables from the _restart_ file

    !! Estimate polynomial parameter for the diameter/rdi curve based on forest
    !management type
    CALL fm_polynomial_parameter_fitting()

    IF (ok_spitfire) THEN
       filename='lightning.nc'
       !Config Key   = LIGHTNING_FILE
       !Config Desc  = The file name from which lightning numbers will be read
       !Config If    = OK_SPITFIRE
       !Config Help  = The file name from which lightning numbers will be read
       CALL getin_p('LIGHTNING_FILE',filename)
       field_name='lightning'
       !Config Key   = LIGHTNING_FIELD_NAME
       !Config Desc  = The field name representing the number of lightnings
       !Config If    = OK_SPITFIRE
       !Config Help  = The field name representing the number of lightnings
       CALL getin_p('LIGHTNING_FIELD_NAME',field_name)
       CALL spitfire_month_input(kjpindex, lalo, neighbours, resolution, &
            contfrac, lightning, filename, field_name)
       CALL xios_orchidee_send_field("interp_diag_lightning",lightning)

       !Config Key   = POPDENS_FILE
       !Config Desc  = The file name from which population density will be read
       !Config If    = OK_SPITFIRE
       !Config Help  = The file name from which population density will be read
       filename='popdens.nc'
       CALL getin_p('POPDENS_FILE',filename)
       field_name='popdens'
       !Config Key   = POPDENS_FIELD_NAME
       !Config Desc  = The field name representing population density
       !Config If    = OK_SPITFIRE
       !Config Help  = The field name representing population density
       CALL getin_p('POPDENS_FIELD_NAME',field_name)
       CALL spitfire_annual_input(kjpindex, lalo, neighbours, resolution, &
            contfrac, pop_dens, filename, field_name)
       CALL xios_orchidee_send_field("interp_diag_popdens",pop_dens)

       !Config Key   = HUMANIGN_FILE
       !Config Desc  = The file name from which the parameter values of a_nd will be read
       !Config If    = OK_SPITFIRE
       !Config Help  = The file name from which the parameter values of a_nd will be read
       filename='humign.nc'
       CALL getin_p('HUMANIGN_FILE',filename)
       field_name='humign'
       !Config Key   = HUMANIGN_FIELD_NAME
       !Config Desc  = The field name representing the a_nd parameter for calculating human ignitions
       !Config If    = OK_SPITFIRE
       !Config Help  = The field name representing the a_nd parameter for calculating human ignitions
       CALL getin_p('HUMANIGN_FIELD_NAME',field_name)
       CALL spitfire_annual_input(kjpindex, lalo, neighbours, resolution, &
            contfrac, a_nd, filename, field_name)
       CALL xios_orchidee_send_field("interp_diag_humign",a_nd)
    ENDIF

    ! Guillaume M. -- Edge_length.nc feeds the SPITFIRE fragmentation flags, the canopy
    ! light edge effect (OK_AED_LIGHT, independent of OK_SPITFIRE) and the recovery target
    ! of OK_AED_FEEDBACK. This block therefore sits outside the IF (ok_spitfire) block.
    IF (ok_aed_humign .OR. ok_aed_fuel .OR. ok_aed_size .OR. ok_aed_wind &
         .OR. ok_aed_light .OR. ok_aed_feedback) THEN
      !Config Key   = EDGE_LENGTH_FILE
      !Config Desc  = The file name from which the parameter values of edge_length will be read
      !Config If    = OK_AED_* (any fragmentation flag, incl. OK_AED_LIGHT)
      !Config Help  = The file name from which the parameter values of edge_length will be read
      filename='Edge_length.nc'
      CALL getin_p('EDGE_LENGTH_FILE',filename)
      field_name='edge_length'
      !Config Key   = EDGE_LENGTH_FIELD_NAME
      !Config Desc  = The field name representing the edge_length parameter
      !Config If    = OK_AED_* (any fragmentation flag, incl. OK_AED_LIGHT)
      !Config Help  = The field name representing the edge_length parameter
      CALL getin_p('EDGE_LENGTH_FIELD_NAME',field_name)
      CALL spitfire_annual_input(kjpindex, lalo, neighbours, resolution, &
            contfrac, edge_length, filename, field_name)
      CALL xios_orchidee_send_field("interp_diag_edge_length",edge_length)

      ! Guillaume M. -- AED_FEEDBACK: refresh the recovery target every time
      ! Edge_length.nc is read, and initialise the dynamic state on first read.
      IF (ok_aed_feedback) THEN
         IF (ALLOCATED(edge_length_ref)) edge_length_ref(:) = edge_length(:)
         ! Guillaume M. -- Seed from the observation ONLY off the v2 path. Under REGROWTH_V2
         ! the edge is fully recomputed each year from A_open and is seeded further down from
         ! the simulated structure: the EFDA observation stays the verification target only.
         IF (ALLOCATED(edge_length_dyn) .AND. .NOT. ok_edge_from_age_class) THEN
            IF (ALL(edge_length_dyn(:) == zero)) edge_length_dyn(:) = edge_length(:)
         ENDIF

         ! Guillaume M. -- No alpha map is read: three of its four fields are near uniform
         ! and are now fixed parameters, and the fourth (harvest) comes from
         ! set_management_intensity.
      ENDIF
    ENDIF

    ! Guillaume M. -- No D_bar map (DIA_INV_FILE) is read: the maturity scale comes entirely
    ! from LARGEST_TREE_DIA per PFT. OK_DIA_STAGGER still exists and requires no file.

    ! Guillaume M. -- Management intensity (design/MODULE_DESIGN_MANAGEMENT_INTENSITY.md):
    ! one annual read only, the class fractions f_class1..5 (Scherpenhuijzen 2025),
    ! interpolated onto the model grid. One concept, one file: set_management_intensity
    ! then derives rotation, alpha_harvest and clearcut_area from it.
    IF (ok_management_intensity .AND. TRIM(management_intensity_file) /= 'NONE') THEN
       IF (ALLOCATED(mi_frac)) THEN
          filename = management_intensity_file
          DO imi = 1, nmiclass
             WRITE(field_name,'(A,I1)') 'f_class', imi
             CALL spitfire_annual_input(kjpindex, lalo, neighbours, resolution, &
                  contfrac, mi_frac(:,imi), filename, field_name)
          ENDDO
       ENDIF
    ENDIF

    ! Guillaume M. -- Regional reference rotation (design/MODULE_DESIGN_ROTATION_REGION.md):
    ! static map, read ONCE here, before the first set_management_intensity below, which
    ! substitutes it to ROTATION_REF(pft) where it carries a value. NONE = never allocated.
    IF (ok_management_intensity .AND. TRIM(rotation_ref_file) /= 'NONE') THEN
       CALL rotation_ref_from_file(kjpindex, lalo, neighbours, resolution, contfrac)
    ENDIF

    ! Guillaume M. -- MACRO edge term (design/MODULE_DESIGN_FRAGMENTATION_3TERMES.md 3.3):
    ! road length per grid cell, read ONCE since the map is static. Two fields are available:
    ! road_length_forest, already restricted to forest and finer than the model fraction but
    ! frozen; road_length_total, raw, to which the model applies its own forest fraction,
    ! which follows land-use change but overestimates (roads are denser outside forest).
    IF (TRIM(road_length_file) /= 'NONE') THEN
       IF (ALLOCATED(road_length_map)) THEN
          filename = road_length_file
          IF (ok_road_forest_mask) THEN
             field_name = 'road_length_forest'
          ELSE
             field_name = 'road_length_total'
          ENDIF
          CALL spitfire_annual_input(kjpindex, lalo, neighbours, resolution, &
               contfrac, road_length_map(:), filename, field_name)
          ! Guillaume M. -- A length cannot be negative; the interpolation can be.
          road_length_map(:) = MAX(zero, road_length_map(:))
          IF (printlev >= 2) THEN
             WRITE(numout,*) 'AED macro : ', TRIM(field_name), ' lue depuis ', TRIM(filename)
             WRITE(numout,*) 'AED macro : longueur de route (m) min/moy/max = ', &
                  MINVAL(road_length_map), &
                  SUM(road_length_map)/REAL(MAX(kjpindex,1),r_std), &
                  MAXVAL(road_length_map)
          ENDIF
       ENDIF
    ENDIF

    ! Guillaume M. -- FRAG_3TERMS permanent: lakes, coastline and rivers. The map already
    ! carries its geometric factors (a river's two banks) and is already restricted to forest,
    ! so it is added AS IS in calculate_aed -- no factor 2, no forest weighting.
    IF (TRIM(edge_hydro_file) /= 'NONE') THEN
       IF (ALLOCATED(edge_hydro_map)) THEN
          filename   = edge_hydro_file
          field_name = 'edge_hydro'
          CALL spitfire_annual_input(kjpindex, lalo, neighbours, resolution, &
               contfrac, edge_hydro_map(:), filename, field_name)
          ! Guillaume M. -- A length cannot be negative; the interpolation can be.
          edge_hydro_map(:) = MAX(zero, edge_hydro_map(:))
          IF (printlev >= 2) THEN
             WRITE(numout,*) 'AED permanent : ', TRIM(field_name), ' lue depuis ', TRIM(filename)
             WRITE(numout,*) 'AED permanent : lisiere eau+rivieres (m) min/moy/max = ', &
                  MINVAL(edge_hydro_map), &
                  SUM(edge_hydro_map)/REAL(MAX(kjpindex,1),r_std), &
                  MAXVAL(edge_hydro_map)
          ENDIF
       ENDIF
    ENDIF

    ! Get values from _restart_ file. Note that only ::kjpindex, ::index, ::lalo
    ! and ::resolution are input variables, all others are output variables.
    CALL readrestart &
         (kjpindex, index, lalo, temp_air, &
         dt_days_read, days_since_beg, &
         adapted, regenerate, &
         vegstress_day, gdd_init_date, litterhum_daily, &
         t2m_daily, t2m_min_daily,t2m_max_daily, tsurf_daily, tsoil_daily, &
         precip_daily,vpd_daily_mean,vpd_daily_max, vpd_mean_week, vpd_max_week, &
         gpp_daily, npp_daily, turnover_daily, turnover_resid, &
         vegstress_month, vegstress_week, vegstress_season,&
         t2m_longterm, tau_longterm, t2m_month, t2m_week, &
         tsoil_month,precip_month, fireindex, firelitter, &
         maxvegstress_lastyear, maxvegstress_thisyear, &
         minvegstress_lastyear, minvegstress_thisyear, &
         maxgppweek_lastyear, maxgppweek_thisyear, &
         gdd0_lastyear, gdd0_thisyear, &
         precip_lastyear, precip_thisyear, &
         gdd_m5_dormance,  gdd_from_growthinit, gdd_midwinter, &
         ncd_dormance, ngd_minus5, &
         PFTpresent, npp_longterm, croot_longterm, n_reserve_longterm, &
         lm_lastyearmax, lm_thisyearmax, &
         maxfpc_lastyear, maxfpc_thisyear, &
         turnover_longterm, gpp_week, gpp_year, gpp_decade, resp_maint_part, resp_maint_week, &
         leaf_age, leaf_frac, leaf_age_crit, plant_status, when_growthinit, age, &
         resp_hetero_d, resp_maint_d, resp_growth_d, co2_fire, &
         veget_lastlight, everywhere, need_adjacent, RIP_time, &
         time_hum_min, hum_min_dormance, litter, dead_leaves, &
         som, lignin_struc, lignin_wood, lignin_snag, turnover_time,&
         fco2_lu, fco2_wh, fco2_ha, &
         prod_s, prod_m, prod_l, flux_s, flux_m, flux_l, &
         fDeforestToProduct, fLulccResidue,fHarvestToProduct, &
         bm_to_litter, bm_to_litter_resid, tree_bm_to_litter, &
         tree_bm_to_litter_resid, carb_mass_total, &
         Tseason, Tseason_length, Tseason_tmp, Tmin_spring_time, &
         global_years, ok_equilibrium, nbp_accu_flux, nbp_pool_start, &
         matrixV, vectorU, previous_stock, current_stock, &
         assim_param, CN_som_litter_longterm, &
         tau_CN_longterm, KF, k_latosa_adapt, &
         rue_longterm, cn_leaf_min_season, nstress_season, &
         soil_n_min, p_O2, bact, forest_managed, &
         species_change_map, fm_change_map, lpft_replant, lai_per_level, &
         laieff_fit, wstress_season, wstress_month, &
         age_stand, age_stand_bm, age_stand_area, rotation_n, last_cut, mai, pai, &
         previous_wood_volume, mai_count, coppice_dens, &
         light_tran_to_floor_season,daylight_count, veget_max, gap_area_save, &
         deepSOM_a, deepSOM_s, deepSOM_p, O2_soil, CH4_soil, O2_snow, CH4_snow, & 
         heat_Zimov, altmax, altmax_ind, altmax_lastyear, altmax_ind_lastyear, &
         depth_organic_soil,fixed_cryoturbation_depth, &
         cn_leaf_init_2D, sugar_load, harvest_cut, &
         harvest_pool_acc, harvest_area_acc, burried_litter, burried_fresh_ltr, &
         burried_fresh_som, burried_bact, &
         burried_min_nitro, burried_som, &
         burried_deepSOM_a, burried_deepSOM_s, burried_deepSOM_p,&
         wood_leftover_legacy, i_beetles_activity_legacy, season_drought_legacy,&
         P_beetles_attacked_legacy, beetle_diapause, sumTeff, woody_litter_to_use, &
         i_beetles_generation, B_beetles_kill_legacy, &
         max_wind_speed_storm, max_wind_ratio_storm, wind_ratio_max_save, & ! JJ 2026: added companions
         wind_ratio_sum_save, wind_speed_max_save, wind_mean_save, wind_sum, &
         is_storm, count_storm, biomass_init_drought, kill_vessels, &
         ni_acc, litterfuel, &
         vessel_loss_previous,grow_season_len, doy_start_gs, doy_end_gs, &
         mean_start_gs, total_ba_init, &
         carbon_acro, carbon_cato, height_acro, deepSOM_peat)

    ! Guillaume M. -- Cold-start stand age is derived, not prescribed: age_stand_bm(c) =
    ! age_class_age_frac(c) * target rotation, so one maturity scale governs area, diameter,
    ! RDI and age together. set_management_intensity is called here because
    ! target_rotation_age is filled only later in stomate_main; the call is idempotent.
    ! /!\ stomate_initialize runs every period: without the restart guard the age is frozen.
    IF (.NOT. stand_age_bm_in_restart .AND. nagec > 1) THEN
       CALL set_management_intensity(kjpindex, veget_max)
       IF (ALLOCATED(target_rotation_age)) THEN
          DO j = 2, nvm
             IF (.NOT. is_tree(j)) CYCLE
             iage_st = j - start_index(agec_group(j)) + 1
             IF (iage_st < 1 .OR. iage_st > nagec) CYCLE
             DO ji = 1, kjpindex
                IF (veget_max(ji,j) > min_stomate .AND. &
                    target_rotation_age(ji,j) > min_stomate) THEN
                   ! Guillaume M. -- The rotation is capped HERE, not in target_rotation_age,
                   ! which must keep its "unmanaged" sentinel for the harvest pace. Uncapped,
                   ! that sentinel gives ages of several hundred years, non-monotonic across
                   ! classes, while the conveyor assumes age grows with class index.
                   age_stand_bm(ji,j) = age_class_age_frac(iage_st) &
                                        * MIN(target_rotation_age(ji,j), &
                                              stand_age_rotation_max)
                   age_stand_area(ji,j) = age_stand_bm(ji,j)
                ENDIF
             ENDDO
          ENDDO
          IF (printlev >= 1) WRITE(numout,*) &
               '[STAND_AGE] age_stand_bm derive de la rotation cible, fractions = ', &
               age_class_age_frac(:)
       ENDIF
    ENDIF

    ! Guillaume M. -- AED: seed edge_length_dyn from the SIMULATED structure. Under
    ! REGROWTH_V2 the edge is fully recomputed each year from A_open, so calling
    ! update_edge_length here duplicates no formula and keeps the seed consistent.
    ! dt = zero: the v2 branch ignores dt, the historical branch stays harmless.
    ! /!\ Must come AFTER set_management_intensity: alpha_harvest and clearcut_area need it.
    IF (ok_aed_feedback .AND. ok_edge_from_age_class .AND. ALLOCATED(edge_length_dyn)) THEN
       IF (ALL(edge_length_dyn(:) == zero)) THEN
          CALL update_edge_length(kjpindex, zero, area, contfrac, veget_max)
          IF (printlev >= 1) WRITE(numout,*) &
               '[AED] edge_length_dyn amorce depuis la structure simulee, moyenne = ', &
               SUM(edge_length_dyn(:)) / MAX(1, kjpindex)
       ENDIF
    ENDIF

    IF (reset_impose_cn) THEN
       DO ipts = 1,kjpindex
         circ_class_biomass(ipts,1,:,ileaf,initrogen) = circ_class_biomass(ipts,1,:,ileaf,icarbon) / &
             cn_leaf_init_2D(ipts,1)
         DO j=2,nvm
            circ_class_biomass(ipts,j,:,ileaf,initrogen) = circ_class_biomass(ipts,j,:,ileaf,icarbon) / &
               cn_leaf_init_2D(ipts,j)
            circ_class_biomass(ipts,j,:,iroot,initrogen) = circ_class_biomass(ipts,j,:,iroot,icarbon) / &
               cn_leaf_init_2D(ipts,j)*fcn_root(j)
            circ_class_biomass(ipts,j,:,ifruit,initrogen) = circ_class_biomass(ipts,j,:,ifruit,icarbon) / &
               cn_leaf_init_2D(ipts,j)*fcn_root(j)
            circ_class_biomass(ipts,j,:,isapabove,initrogen) = circ_class_biomass(ipts,j,:,isapabove,icarbon) / &
               cn_leaf_init_2D(ipts,j)*fcn_wood(j)
            circ_class_biomass(ipts,j,:,isapbelow,initrogen) = circ_class_biomass(ipts,j,:,isapbelow,icarbon) / &
               cn_leaf_init_2D(ipts,j)*fcn_wood(j)
            circ_class_biomass(ipts,j,:,iheartabove,initrogen) = circ_class_biomass(ipts,j,:,iheartabove,icarbon) / &
               cn_leaf_init_2D(ipts,j)*fcn_wood(j)
            circ_class_biomass(ipts,j,:,iheartbelow,initrogen) = circ_class_biomass(ipts,j,:,iheartbelow,icarbon) / &
               cn_leaf_init_2D(ipts,j)*fcn_wood(j)
         END DO
       END DO
    ENDIF

    !! Calculate som_total to be used when USE_SOILC_INSULATION=y and USE_REFSOC=n
    IF (ok_soil_carbon_discretization) THEN
       som_total(:,:,:,:) = deepSOM_a(:,:,:,:) + deepSOM_s(:,:,:,:) + deepSOM_p(:,:,:,:)
    ELSE
       som_total(:,:,:,:) = zero
    END IF


    !! 1.4.5 Check time step
       
    !! 1.4.5.1 Allow STOMATE's time step to change although this is dangerous
    IF (dt_days /= dt_days_read) THEN
       WRITE(numout,*) 'slow_processes: STOMATE time step changes:', &
            & dt_days_read,' -> ',dt_days
    ENDIF
    
    !! 1.4.5.2 Time step has to be a multiple of a full day
    IF ( ( dt_days-REAL(NINT(dt_days),r_std) ) > min_stomate ) THEN
       WRITE(numout,*) 'slow_processes: STOMATE time step is not a mutiple of a full day:', &
            & dt_days,' days.'
       STOP
    ENDIF
    
    !! 1.4.5.3 upper limit to STOMATE's time step
    IF ( dt_days > max_dt_days ) THEN
       WRITE(numout,*) 'slow_processes: STOMATE time step exceeds the maximum value:', &
            & dt_days,' days > ', max_dt_days, ' days.'  
       STOP
    ENDIF
    
    !! 1.4.5.4 STOMATE time step must not be less than the forcing time step
    IF ( dt_sechiba > dt_days*one_day ) THEN
       WRITE(numout,*) &
            & 'slow_processes: STOMATE time step ::dt_days smaller than forcing time step ::dt_sechiba'
       STOP
    ENDIF
    
    !! 1.4.5.6 Final message on time step
    IF (printlev >=2) WRITE(numout,*) 'Slow_processes, STOMATE time step (days): ', dt_days
    

    ! 1.4.7b Write forcing file for the soil carbon discretization module 
    IF ( ok_soil_carbon_discretization ) THEN
       !Config  Key  = OK_FORCESOIL_WRITE
       !Config  Desc = .TRUE. or .FALSE. If .TRUE.,
       !the file Cforcing_duscretization_name will be written.
       !Config  If   = OK_SOIL_CARBON_DISCRETIZATION
       !Config  Def  = .FALSE.
       !Config  Help = Flag to write Forcesoil forcing file
       ok_forcesoil_write = .FALSE.
       CALL getin_p('OK_FORCESOIL_WRITE', ok_forcesoil_write)

       IF ( ok_forcesoil_write) THEN

          !Config  Key  = STOMATE_CFORCING_NAME
          !Config  Desc = Name of STOMATE's carbon forcing file or NONE. If NONE the file will not be written.
          !Config  If   = OK_SOIL_CARBON_DISCRETIZATION
          !Config  Def  = stomate_cforcing.nc
          !Config  Help = Name that will be given to STOMATE's carbon soil discretization 
          !Config         offline forcing file
          Cforcing_discretization_name = 'stomate_cforcing.nc'
          CALL getin ('STOMATE_CFORCING_NAME', Cforcing_discretization_name)
          IF (printlev >=2) WRITE(numout,*) 'Writing of forcing file for forcesoil will be done in file:', Cforcing_discretization_name


         ! Time step of forcesoil
         !Config Key   = FORCESOIL_STEP_PER_YEAR
         !Config Desc  = Number of time steps per year for carbon spinup.
         !Config If    = STOMATE_CFORCING_NAME and OK_STOMATE and OK_SOIL_CARBON_DISCRETIZATION
         !Config Def   = 365 (366, ...)
         !Config Help  = Number of time steps per year for carbon spinup.
         !Config Units = [days, months, year]
         nparan = 365 !year_length_in_days
         CALL getin_p('FORCESOIL_STEP_PER_YEAR', nparan)
         
         ! Correct if setting is out of bounds 
         IF ( nparan < 1 ) THEN
            WRITE(temp_str, *) "Value found:", nparan
            CALL ipslerr_p(3, 'stomate_initialize', &
                  'Invalid value for FORCESOIL_STEP_PER_YEAR ', &
                  'Expected value is > 0', temp_str)
         ENDIF

         !Config Key   = FORCESOIL_NB_YEAR
         !Config Desc  = Number of years saved for carbon spinup.
         !Config If    = STOMATE_CFORCING_NAME and OK_STOMATE
         !Config Def   = 1
         !Config Help  = Number of years saved for carbon spinup. If internal parameter cumul_Cforcing is TRUE in stomate.f90
         !Config         Then this parameter is forced to one.
         !Config Units = [years]
         nbyear = 1
         CALL getin_p('FORCESOIL_NB_YEAR', nbyear)

         ! Make use of ::nparan to calculate ::dt_forcesoil
         dt_forcesoil = zero
         nparan = nparan+1
         DO WHILE ( dt_forcesoil < dt_stomate/one_day )
            nparan = nparan-1
            IF ( nparan < 1 ) THEN
               CALL ipslerr_p(3,'stomate_initialize','Problem with number of soil forcing time steps','nparan < 1','')
            ENDIF
            dt_forcesoil = one_year/REAL(nparan,r_std)
         ENDDO
         IF ( nparan > nparanmax ) THEN
           CALL ipslerr_p(3,'stomate_initialize','Problem with number of soil forcing time steps','nparan > nparanmax','')
         ENDIF
         WRITE(numout,*) 'Time step of soil forcing (d): ',dt_forcesoil

         IF (is_root_prc) CALL SYSTEM ('rm -f '//TRIM(Cforcing_discretization_name))
       
         ALLOCATE( nforce(nparan*nbyear), stat=ier)
         IF (ier /= 0) CALL ipslerr_p(3, 'stomate_initialize', 'Problem allocating nforce', 'Error code=', ier)
         ALLOCATE(som_input_2pfcforcing(kjpindex,ncarb,nvm,nelements,nparan*nbyear))
         ALLOCATE(pb_2pfcforcing(kjpindex,nparan*nbyear))
         ALLOCATE(snow_2pfcforcing(kjpindex,nparan*nbyear))
         ALLOCATE(tprof_2pfcforcing(kjpindex,ngrnd,nvm,nparan*nbyear))
         ALLOCATE(fbact_2pfcforcing(kjpindex,ngrnd,nvm,nparan*nbyear))
         ALLOCATE(hslong_2pfcforcing(kjpindex,ngrnd,nvm,nparan*nbyear))
         ALLOCATE(veget_max_2pfcforcing(kjpindex,nvm,nparan*nbyear))
         ALLOCATE(rprof_2pfcforcing(kjpindex,nvm,nparan*nbyear))
         ALLOCATE(tsurf_2pfcforcing(kjpindex,nparan*nbyear))
         ALLOCATE(snowdz_2pfcforcing(kjpindex,nsnow,nparan*nbyear))
         ALLOCATE(snowrho_2pfcforcing(kjpindex,nsnow,nparan*nbyear))
         ALLOCATE(CN_target_2pfcforcing(kjpindex,nvm,ncarb,nparan*nbyear))
         ALLOCATE(n_mineralisation_2pfcforcing(kjpindex,nvm,nparan*nbyear))
         nforce(:) = zero
         som_input_2pfcforcing(:,:,:,:,:) = zero
         pb_2pfcforcing(:,:) = zero
         snow_2pfcforcing(:,:) = zero
         tprof_2pfcforcing(:,:,:,:) = zero
         fbact_2pfcforcing(:,:,:,:) = zero
         hslong_2pfcforcing(:,:,:,:) = zero
         veget_max_2pfcforcing(:,:,:) = zero
         rprof_2pfcforcing(:,:,:) = zero
         tsurf_2pfcforcing(:,:) = zero
         snowdz_2pfcforcing(:,:,:) = zero
         snowrho_2pfcforcing(:,:,:) = zero
         CN_target_2pfcforcing(:,:,:,:) = zero
         n_mineralisation_2pfcforcing(:,:,:) = zero

       ENDIF ! TRIM(Cforcing_discretization_name) /= 'NONE'
    ENDIF ! ok_soil_carbon_discretization
   
    IF (ok_peat_NoDiscretisation) THEN
       ALLOCATE(wtp_pt(kjpindex,nparan*nbyear), stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'stomate_initialize','Pb in alloc for wtp_pt','','')
       wtp_pt(:,:)=zero
    ENDIF

    !! 1.4.9 Initialize non-zero variables
    CALL stomate_var_init &
         (kjpindex, veget_max, leaf_age, leaf_frac, &
         leaf_age_crit, dead_leaves, &
         veget, deadleaf_cover, assim_param, &
         circ_class_biomass, circ_class_n, sugar_load)
          
    ! Initialize temp_growth
    temp_growth(:)=t2m_month(:)-tp_00 

   !Config Key   = FROZEN_RESPIRATION_FUNC 
   !Config Desc  = Method for soil decomposition function 
   !Config If    = OK_SOIL_CARBON_DISCRETIZATION 
   !Config Def   = 1
   !Config Help  = 
   !Config Units = [1]
   frozen_respiration_func = 0
   CALL getin_p('FROZEN_RESPIRATION_FUNC',frozen_respiration_func)
   IF (printlev >=2) WRITE(numout, *)' frozen soil respiration function:  ', frozen_respiration_func
      
  END SUBROUTINE stomate_initialize
  

!! ================================================================================================================================
!! SUBROUTINE 	: stomate_main
!!
!>\BRIEF        Manages variable initialisation, reading and writing forcing 
!! files, aggregating data at stomate's time step (dt_stomate), aggregating data
!! at longer time scale (i.e. for phenology) and uses these forcing to calculate
!! CO2 fluxes (NPP and respirations) and C-pools (litter, soil, biomass, ...)
!!
!! DESCRIPTION  : The subroutine manages 
!! divers tasks:
!! (1) Initializing all variables of stomate (first call)
!! (2) Reading and writing forcing data (last call)
!! (3) Adding CO2 fluxes to the IPCC history files
!! (4) Converting the time steps of variables to maintain consistency between
!! sechiba and stomate
!! (5) Use these variables to call stomate_lpj, maint_respiration, littercalc,
!! som. The called subroutines handle: climate constraints 
!! for PFTs, PFT dynamics, Phenology, Allocation, NPP (based on GPP and
!! authothropic respiration), fire, mortality, vmax, assimilation temperatures,
!! all turnover processes, light competition, sapling establishment, lai,  
!! land cover change and litter and soil dynamics.
!! (6) Use the spin-up method developed by Lardy (2011)(only if SPINUP_ANALYTIC 
!! is set to TRUE).
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): deadleaf_cover, assim_param, veget, 
!! veget_max, resp_maint, resp_hetero, resp_growth, 
!! fco2_flux_out, fco2_lu_out, fco2_wh_out, fco2_ha_out.
!!
!! REFERENCES	: 
!! - Lardy, R, et al., A new method to determine soil organic carbon equilibrium,
!! Environmental Modelling & Software (2011), doi:10.1016|j.envsoft.2011.05.016
!!
!! FLOWCHART    : 
!! \latexonly 
!! \includegraphics[scale=0.5]{stomatemainflow.png}
!! \endlatexonly
!! \n
!_ ================================================================================================================================
  
  SUBROUTINE stomate_main &
       & (kjit, kjpij, kjpindex, njsc, &
       &  index, lalo, neighbours, resolution, contfrac, frac_nobio, &
       &  totfrac_nobio, clay, &
       &  silt, bulk, temp_air, temp_sol, stempdiag, precip_longterm, &
       &  vegstress, humrel, & 
       &  shumdiag, litterhumdiag, precip_rain, precip_snow, &
       &  tmc_pft, drainage_pft, runoff_pft, swc_pft, gpp, VPD, deadleaf_cover, &
       &  assim_param, qsintveg, &
       &  frac_age, veget, veget_max, &
       &  veget_max_new, loss_gain, frac_nobio_new, fraclut, &
       &  rest_id_stom, hist_id_stom, hist_id_stom_IPCC, &
       &  fco2_flux_out, fco2_lu_out, fco2_wh_out, fco2_ha_out, &
       &  resp_maint,resp_hetero,resp_growth,temp_growth, &
       &  soil_pH, pb, n_input_raw, month, &
       &  tdeep, hsdeep, snow, heat_Zimov, sfluxCH4_deep, sfluxCO2_deep, &  
       &  som_total, snowdz, snowrho, altmax, depth_organic_soil, cn_leaf_min_2D, cn_leaf_max_2D, cn_leaf_init_2D, &
       &  circ_class_biomass, circ_class_n, lai_per_level, &
       &  laieff_fit, Light_Abs_Tot, Light_Tran_Tot, &
       &  laieff_isotrop, z_array_out, max_height_store, &
       &  transpir, transpir_mod, &
       &  coszang,stressed, unstressed, &
       &  u, v, mcs_hydrol, &
       &  mcfc_hydrol, vessel_loss, root_profile, root_depth, us, &
       &  Pgap_cumul, z0m,&
       &  wtp,  shumdiag_peat,  mc_peat_above,&
       &  liqwt_ratio,  shumdiag_croppeat, &
       &  mc_croppeat_above, &
       &  shumdiag_man,  mc_man_above, soiltile)


    IMPLICIT NONE


    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std),INTENT(in)                       :: kjit              !! Time step number (unitless)
    INTEGER(i_std),INTENT(in)                       :: kjpindex          !! Domain size - terrestrial pixels only (unitless)
    INTEGER(i_std),INTENT(in)                       :: kjpij             !! Total size of the un-compressed grid (unitless)
    INTEGER(i_std),DIMENSION(:), INTENT (in)        :: njsc              !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless)
    INTEGER(i_std),INTENT(in)                       :: rest_id_stom      !! STOMATE's _Restart_ file identifier (unitless)
    INTEGER(i_std),INTENT(in)                       :: hist_id_stom      !! STOMATE's _history_ file identifier (unitless)
    INTEGER(i_std),INTENT(in)                       :: hist_id_stom_IPCC !! STOMATE's IPCC _history_ file identifier 
                                                                         !! (unitless) 
    INTEGER(i_std),DIMENSION(:),INTENT(in)          :: index             !! Indices of the pixels on the map. Stomate uses a 
                                                                         !! reduced grid excluding oceans. ::index contains 
                                                                         !! the indices of the terrestrial pixels only 
                                                                         !! (unitless) 
    INTEGER(i_std),DIMENSION(:,:),INTENT(in)        :: neighbours        !! Neighoring grid points if land for the DGVM 
                                                                         !! (unitless) 
    REAL(r_std),DIMENSION(:,:),INTENT(in)           :: lalo              !! Geographical coordinates (latitude,longitude) 
                                                                         !! for pixels (degrees) 
    REAL(r_std),DIMENSION(:,:),INTENT(in)           :: resolution        !! Size in x an y of the grid (m) - surface area of 
                                                                         !! the gridbox 
    REAL(r_std),DIMENSION (:), INTENT (in)          :: contfrac          !! Fraction of continent in the grid cell (unitless)
    REAL(r_std),DIMENSION (:), INTENT (in)          :: totfrac_nobio     !! Total fraction of ice+lakes+cities etc. in the mesh
    REAL(r_std),DIMENSION(:),INTENT(in)             :: clay              !! Clay fraction of soil (0-1, unitless)
    REAL(r_std),DIMENSION(:),INTENT(in)             :: silt              !! Silt fraction of soil (0-1, unitless) 
    REAL(r_std),DIMENSION(:),INTENT(in)             :: bulk              !! Bulk density (kg/m**3) 
    REAL(r_std),DIMENSION(:,:),INTENT(inout)        :: vegstress         !! Relative soil moisture (0-1, unitless)
    REAL(r_std),DIMENSION(:,:),INTENT (inout)       :: humrel            !! Relative humidity - not used in stomate (needed in age_class_distr)
    REAL(r_std),DIMENSION(:),INTENT(in)             :: temp_air          !! Air temperature at first atmosperic model layer (K)
    REAL(r_std),DIMENSION(:),INTENT(in)             :: temp_sol          !! Surface temperature (K)
    REAL(r_std),DIMENSION(:,:),INTENT(in)           :: stempdiag         !! Soil temperature (K)
    REAL(r_std),DIMENSION(:,:),INTENT(in)           :: shumdiag          !! Relative soil moisture (0-1, unitless)
    REAL(r_std),DIMENSION(:),INTENT(in)             :: litterhumdiag     !! Litter humidity (0-1, unitless)
    REAL(r_std),DIMENSION(:,:),INTENT(in)           :: transpir          !! transpiration @tex $(kg m^{-2} timestep^{-1})$
    REAL(r_std),DIMENSION(:,:),INTENT(in)           :: transpir_mod      !! transpir divided by veget_max
    REAL(r_std),DIMENSION(:),INTENT(in)             :: precip_rain       !! Rain precipitation  
                                                                         !! @tex $(mm dt_stomate^{-1})$ @endtex 
    REAL(r_std),DIMENSION(:),INTENT(in)             :: precip_snow       !! Snow precipitation  
                                                                         !! @tex $(mm dt_stomate^{-1})$ @endtex 
    REAL(r_std),DIMENSION(:),INTENT (in)            :: u                 !! Lowest level wind speed in direction u (m/s)
    REAL(r_std),DIMENSION(:),INTENT (in)            :: v                 !! Lowest level wind speed in direction v (m/s)
    REAL(r_std), DIMENSION (:,:), INTENT(in)        :: tmc_pft           !! Total soil water per PFT (mm/m2) 
    REAL(r_std), DIMENSION (:,:), INTENT(in)        :: drainage_pft      !! Drainage per PFT (mm/m2)   
    REAL(r_std), DIMENSION (:,:), INTENT(in)        :: runoff_pft        !! Runoff per PFT (mm/m2)   
    REAL(r_std), DIMENSION (:,:), INTENT(in)        :: swc_pft           !! Relative Soil water content [tmcr:tmcs] per pft (-)     
    REAL(r_std),DIMENSION(:,:),INTENT(in)           :: gpp               !! GPP of total ground area  
                                                                         !! @tex $(gC m^{-2} time step^{-1})$ @endtex 
                                                                         !! Calculated in sechiba, account for vegetation 
                                                                         !! cover and effective time step to obtain ::gpp_d
    REAL(r_std), DIMENSION (:), INTENT(in)          :: VPD               !! Vapor Pressure Deficit (kPa) 
    REAL(r_std),DIMENSION(:,:),INTENT(inout)        :: frac_nobio_new    !! New fraction of nobio per gridcell
    REAL(r_std),DIMENSION(:),INTENT(in)             :: soil_pH           !! soil pH 	 
    REAL(r_std),DIMENSION(:), INTENT(in)            :: pb                !! Air pressure (hPa) 
    REAL(r_std),DIMENSION(:,:,:,:),INTENT(in)       :: n_input_raw       !! Nitrogen inputs into the soil  (gN/m**2/timestep) 
    REAL(r_std),DIMENSION(:,:), INTENT(in)          :: cn_leaf_min_2D    !! minimal leaf C/N ratio 
    REAL(r_std),DIMENSION(:,:), INTENT(in)          :: cn_leaf_max_2D    !! maximal leaf C/N ratio 
    REAL(r_std),DIMENSION(:,:), INTENT(in)          :: cn_leaf_init_2D   !! initial leaf C/N ratio 
    REAL(r_std),DIMENSION(:), INTENT(in)            :: mcs_hydrol        !! Saturated volumetric water content output to be used in stomate_soilcarbon
    REAL(r_std),DIMENSION(:), INTENT(in)            :: mcfc_hydrol       !! Volumetric water content at field capacity output to be used in stomate_soilcarbon
    REAL(r_std),DIMENSION(:,:),INTENT(in)           :: fraclut           !! Fraction of landuse tiles
    REAL(r_std),DIMENSION(:), INTENT(in)            :: coszang           !! the cosine of the zenith angle
    REAL(r_std), DIMENSION(:,:), INTENT(inout)      :: loss_gain         !! losses and gains due to LCC distributed over all
                                                                         !! age classes and thus taking the age-classes into 
                                                                         !! account (unitless, 0-1) 
    REAL(r_std),DIMENSION(:,:),INTENT(inout)        :: veget_max_new     !! New "maximal" coverage fraction of a PFT: only if 
                                                                         !! vegetation is updated in slowproc
    ! Variables for soil carbon discretization
    REAL(r_std), DIMENSION(:,:), INTENT(in)         :: snowdz            !! snow depth [m]
    REAL(r_std), DIMENSION(:,:), INTENT(in)         :: snowrho           !! snow density (Kg/m^3)
    REAL(r_std), DIMENSION(:,:,:), INTENT (in)      :: tdeep             !! deep temperature profile (K)
    REAL(r_std), DIMENSION(:,:,:), INTENT (in)      :: hsdeep            !! deep long term soil humidity profile (unitless)
    REAL(r_std), DIMENSION(:,:,:), INTENT (out)     :: heat_Zimov        !! heating associated with decomposition [W/m**3 soil]
    REAL(r_std), DIMENSION(:), INTENT (out)         :: sfluxCH4_deep     !! surface flux of CH4 to atmosphere from soil
    REAL(r_std), DIMENSION(:), INTENT (out)         :: sfluxCO2_deep     !! surface flux of CO2 to atmosphere from soil
    REAL(r_std), DIMENSION(:), INTENT (in)          :: snow              !! Snow mass [Kg/m^2]
    REAL(r_std), DIMENSION(:,:,:,:), INTENT (inout) :: som_total         !! total soil carbon for use in thermal calcs (g/m**3)
    REAL(r_std), DIMENSION(:,:),INTENT(inout)       :: altmax            !! Maximul active layer thickness (m). Be careful, here active means non frozen.
                                                                         !! Not related with the active soil carbon pool.
    REAL(r_std), DIMENSION(:),   INTENT (inout)     :: depth_organic_soil!! Depth at which there is still organic matter (m) 
    REAL(r_std), DIMENSION(:,:,:),   INTENT (in)    :: vessel_loss       !! Proportion of conductivity lost due to cavitation in the xylem (no unit).
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)     :: root_profile      !! Normalized root mass/length fraction in each soil layer 
                                                                         !! (0-1, unitless)
    REAL(r_std), DIMENSION (:,:,:), INTENT(in)      :: root_depth        !! Node and interface numbers at which the deepest roots
                                                                         !! occur (1 to nslm, unitless)
    INTEGER(i_std), INTENT(in)                      :: month             !! month number required for n_input (1-12)
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)       :: Pgap_cumul        !! The probability of finding a gap in the
                                                                         !! in canopy from the top of the canopy
                                                                         !! to a given level.
                                                                         !! (unitless, between 0-1)
    REAL(r_std), DIMENSION(:), INTENT(in)           :: z0m               !! Surface roughness for momentum (m)

    !! peatland 
    REAL(r_std),DIMENSION (kjpindex,nslm),INTENT(in) :: shumdiag_peat    !! soil moisture content useful for the decompistion when peat is activated
    REAL(r_std),DIMENSION (kjpindex,nslm),INTENT(in) :: shumdiag_croppeat!! soil moisture content useful for the decompistion when agri_peat is activated
    REAL(r_std),DIMENSION(kjpindex),INTENT(in)       :: liqwt_ratio      !! liquid water ratio in peat, compared to the satured one
    REAL(r_std),DIMENSION(kjpindex),INTENT(in)       :: wtp              !! water table position when peat is activated
    REAL(r_std),DIMENSION(kjpindex),INTENT(in)       :: mc_peat_above    !! liquid water content within the first 4 layers when peat is activated
    REAL(r_std),DIMENSION(kjpindex),INTENT(in)       :: mc_croppeat_above!! liquid water content within the first 4 layers when agri_peat activated
    REAL(r_std),DIMENSION(kjpindex,ngrnd)            :: mc_peat          !! layer 1-11: equivalent to shumdiag_peat at the corresponing level
                                                                         !! other layers: the value at last layer of shumdiag_peat
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm)       :: hsdeep_new       !! deep long term soil humidity profile
    ! tides
    REAL(r_std),DIMENSION (kjpindex,nslm),INTENT(in) :: shumdiag_man     !! soil moisture content useful for the decompistion when tides activated
    REAL(r_std),DIMENSION(kjpindex,ngrnd)            :: mc_man           !! equivalent to shumdiag_man
    REAL(r_std),DIMENSION(kjpindex),INTENT(in)       :: mc_man_above      !! liquid water content within the first 4 layers when tides  activated
    REAL(r_std),DIMENSION (:,:), INTENT(in)             :: soiltile       !! Fraction of soil tile within vegtot (0-1, unitless)
    
    !! 0.2 Output variables

    REAL(r_std),DIMENSION(:),INTENT(out)            :: fco2_flux_out     !! CO2 flux between atmosphere and biosphere per 
                                                                         !! average ground area 
                                                                         !! @tex $(gC m^{-2} dt_sechiba^{-1})$ @endtex  
    REAL(r_std),DIMENSION(:),INTENT(out)            :: fco2_lu_out       !! CO2 flux between atmosphere and biosphere from 
                                                                         !! land-use (without forest management) (gC/m2/dt_stomate)
    REAL(r_std),DIMENSION(:),INTENT(out)            :: fco2_wh_out       !! CO2 Flux to Atmosphere from Wood Harvesting (gC/m2/dt_stomate)
    REAL(r_std),DIMENSION(:),INTENT(out)            :: fco2_ha_out       !! CO2 Flux to Atmosphere from Crop Harvesting (gC/m2/dt_stomate)
    REAL(r_std),DIMENSION(:,:),INTENT(out)          :: resp_maint        !! Maitenance component of autotrophic respiration in 
                                                                         !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex 
    REAL(r_std),DIMENSION(:,:),INTENT(out)          :: resp_growth       !! Growth component of autotrophic respiration in 
                                                                         !! @tex ($gC m^{-2} dt_stomate^{-1}$) @endtex
    REAL(r_std),DIMENSION(:,:),INTENT(out)          :: resp_hetero       !! Heterotrophic respiration in  
                                                                         !! @tex $(gC m^{-2} dt_stomate^{-1})$ @endtex  
    REAL(r_std),DIMENSION(:),INTENT(out)            :: temp_growth       !! Growth temperature (Â°C)  
                                                                         !! Is equal to t2m_month 
    REAL(r_std),DIMENSION(:,:),INTENT(out)          :: max_height_store  !! ???


    !! 0.3 Modified

    REAL(r_std),DIMENSION(:,:),INTENT(in)           :: veget             !! Fraction of vegetation type including 
                                                                         !! non-biological fraction (unitless) 
    REAL(r_std),DIMENSION(:,:),INTENT(inout)        :: veget_max         !! Maximum fraction of vegetation type including 
                                                                         !! non-biological fraction (unitless) 
    REAL(r_std),DIMENSION(:,:,:),INTENT(inout)      :: assim_param       !! vmax, nue and leaf N for photosynthesis
                                                                         !! @tex $(\mu mol m^{-2}s^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)      :: qsintveg          !! Water on vegetation due to interception @tex $(kg m^{-2})$ @endtex
    REAL(r_std),DIMENSION(:),INTENT(inout)          :: deadleaf_cover    !! Fraction of soil covered by dead leaves 
                                                                         !! (unitless) 
    REAL(r_std),DIMENSION(:,:,:),INTENT(inout)      :: frac_age          !! Age efficacity from STOMATE     
    REAL(r_std), DIMENSION(:,:,:,:,:),INTENT(inout) :: circ_class_biomass!! Biomass per circumference class @tex $(gC tree^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:,:),INTENT(inout)     :: circ_class_n      !! Number of trees within each circumference - includes grasses 
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)    :: lai_per_level     !! This is the LAI per vertical level
                                                                         !! @tex $(m^{2} m^{-2})$ @endtex
    TYPE(laieff_type),DIMENSION (:,:,:),INTENT(inout) :: laieff_fit      !! Fitted parameters for the effective LAI
    REAL(r_std),DIMENSION (:,:,:), INTENT (inout)   :: Light_Abs_Tot     !! Absorbed radiation per level for photosynthesis
    REAL(r_std),DIMENSION (:,:,:), INTENT (inout)   :: Light_Tran_Tot    !! Transmitted radiation per level for photosynthesis
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)    :: laieff_isotrop    !! Effective LAI
    REAL(r_std),DIMENSION(:,:),INTENT(inout)        :: stressed          !! adjusted ecosystem functioning. Takes the unit of the variable
                                                                         !! used as a proxy for waterstress (assigned in sechiba).
    REAL(r_std),DIMENSION(:,:),INTENT(inout)        :: unstressed        !! initial ecosystem functioning after the first calculation and 
                                                                         !! before any recalculations. Takes the unit of the variable used
                                                                         !! as a proxy for unstressed (assigned in sechiba).
    REAL(r_std),DIMENSION(:,:,:,:), INTENT(inout)   :: z_array_out       !! An output of h_array, to use in sechiba
    REAL(r_std),DIMENSION(:,:),INTENT(inout)        :: frac_nobio        !! Fraction of grid cell covered by lakes, land 
                                                                         !! ice, cities, ... (unitless)
    REAL(r_std),DIMENSION(:,:,:,:), INTENT(inout)   :: us                !! Water stress index for transpiration
                                                                         !! (by soil layer and PFT) (0-1, unitless)
    REAL(r_std),DIMENSION(:),INTENT(inout)          :: precip_longterm   !! Longterm annual precipitation sum mm year^{-1}


    !! 0.4 local variables

    CHARACTER(LEN=10), DIMENSION(nelements)       :: element_str              !! string suffix indicating element    
    REAL(r_std)                                   :: dt_days_read             !! STOMATE time step read in restart file (days)
    INTEGER(i_std)                                :: l,k,ji, jv, i, j, m, ip      !! indices
    INTEGER(i_std)                                :: igrn                     !! indices
    REAL(r_std),PARAMETER                         :: max_dt_days = 5.         !! Maximum STOMATE time step (days)
    REAL(r_std)                                   :: hist_days                !! Writing frequency for history file (days)
    REAL(r_std),DIMENSION(0:nslm)                 :: z_soil                   !! Variable to store depth of the different soil 
    !! layers (m)
    REAL(r_std),DIMENSION(kjpindex)               :: cvegtot                  !! Total "vegetation" cover (unitless)
    REAL(r_std),DIMENSION(kjpindex)               :: precip                   !! Total liquid and solid precipitation  
    !! @tex $(??mm dt_stomate^{-1})$ @endtex 
    REAL(r_std),DIMENSION(kjpindex,nvm)           :: gpp_d                    !! Gross primary productivity per ground area 
    !! @tex $(??gC m^{-2} dt_stomate^{-1})$ @endtex  
    REAL(r_std),DIMENSION(kjpindex,nvm)           :: gpp_daily_x              !! "Daily" gpp for teststomate  
    !! @tex $(??gC m^{-2} dt_stomate^{-1})$ @endtex 
    REAL(r_std),DIMENSION(kjpindex,nvm)           :: resp_hetero_litter       !! Litter heterotrophic respiration per ground area 
    !! @tex $(gC m^{-2} day^{-1})$ @endtex  
    !! ??Same variable is also used to 
    !! store heterotrophic respiration per ground area 
    !! over ::dt_sechiba?? 
    REAL(r_std),DIMENSION(nvm)                    :: ld_redistribute          !! logical set to redistribute som and litter
    REAL(r_std),DIMENSION(kjpindex,nvm)           :: resp_hetero_soil         !! soil heterotrophic respiration  
    !! @tex $(gC m^{-2} day^{-1})$ @endtex
    REAL(r_std),DIMENSION(kjpindex,nlevs)         :: control_moist_inst       !! Moisture control of heterotrophic respiration 
    !! (0-1, unitless) 
    REAL(r_std),DIMENSION(kjpindex,nlevs)         :: control_temp_inst        !! Temperature control of heterotrophic 
    !! respiration, above and below (0-1, unitless) 
    REAL(r_std),DIMENSION(kjpindex,ncarb,nvm,nelements) :: som_input_inst     !! Quantity of carbon going into carbon pools from 
    !! litter decomposition 
    !! @tex $(gC m^{-2} day^{-1})$ @endtex
    INTEGER(i_std)                                :: ier                      !! Check errors in netcdf call (unitless)
    REAL(r_std)                                   :: sf_time                  !! Intermediate variable to calculate current time step 
    REAL(r_std), DIMENSION(kjpindex)              :: vartmp                   !! Temporary variable
    INTEGER(i_std)                                :: nneigh                   !! Number of neighbouring pixels
    REAL(r_std), DIMENSION(kjpindex,nvm)           :: longevity_eff_root                 !! Effective root turnover time that accounts
    !! waterstress (days)
    REAL(r_std), DIMENSION(kjpindex,nvm)           :: longevity_eff_sap                  !! Effective sapwood turnover time that accounts
    !! waterstress (days)
    REAL(r_std), DIMENSION(kjpindex,nvm)           :: longevity_eff_leaf                 !! Effective leaf turnover time that accounts
    !! waterstress (days)
    REAL(r_std), DIMENSION(kjpindex,nvm)           :: wstress_adapt                !! Factor to account for a long acclimation of
    !! of the PFT to the long-term waterstress in
    !! the pixel
    REAL(r_std), DIMENSION(kjpindex,nvm,nionspec)  :: leaching                     !! mineral nitrogen leached from the soil
    REAL(r_std), DIMENSION(kjpindex,nvm,nnspec)    :: emission                     !! volatile losses of nitrogen (gN/m**2/timestep)
    REAL(r_std), DIMENSION(kjpindex,nvm,nmbcomp,nelements) &
         :: check_intern                 !! Contains the components of the internal
    !! mass balance chech for this routine
    !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(kjpindex,nvm,nelements) :: closure_intern               !! Check closure of internal mass balance
                                                                                   !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(kjpindex,nvm,nelements) :: pool_start                   !! Start and end pool of this routine 
                                                                                   !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(kjpindex,nvm,nelements) :: pool_end                     !! Start and end pool of this routine 
                                                                                   !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(kjpindex,nvm)           :: veget_max_begin              !! veget_max at the start of the routine.
                                                                                   !! Used for consistency checks
    INTEGER(i_std)                                 :: inspec, ininput, inionspec   !! Indices
    INTEGER(i_std)                                 :: ipts, ivm, ilitt, ilev, icir !! Indices
    INTEGER(i_std)                                 :: icarb, ipar, iele, imbc      !! Indices
    INTEGER(i_std)                                 :: iday, iyr, ifm               !! Indices
    REAL(r_std),DIMENSION(kjpindex,nvm)           :: count_daylight                !! Time steps dt_radia during daylight
    REAL(r_std), DIMENSION(kjpindex,nvm,nbpools)  :: carbon_stock                  !! Array containing the carbon stock for each pool
                                                                                   !! used by ORCHIDEE
    REAL(r_std), DIMENSION(kjpindex,nvm)          :: n_mineralisation              !! net nitrogen mineralisation of decomposing SOM 
                                                                                   !!   (gN/m**2/day), supposed to be NH4 
    REAL(r_std), DIMENSION(kjpindex,nvm,nionspec) :: plant_n_uptake                !! Uptake of soil N by plants  
                                                                                   !! (gN/m**2/timestep)  
    REAL(r_std),DIMENSION(kjpindex,nvm)           :: resp_total_soil               !! soil heterotrophic respiration (gC/day/m**2 by PFT) 
    REAL(r_std), DIMENSION(kjpindex,nvm,ncarb)    :: CN_target                     !! C to N ratio of SOM flux from one pool to another (gN m-2 dt-1) 
    REAL(r_std)                                   :: weight_spinup                 !! How do we account for spinup computation (0-1)
    LOGICAL                                       :: partial_spinup                !! in order to spinup only slow and passive pools
    LOGICAL                                       :: nitrogen_spinup               !! in order to spinup only carbon pools
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm)    :: tdeep_celsius                 !! deep temperature profile celsius (C)
    REAL(r_std), DIMENSION(kjpindex)              :: tsoil_decomp                  !! Temperature used for decompostition in soil (K)
    REAL(r_std), DIMENSION(kjpindex,nvm,nlevs,nelements) :: snag_to_wood_flux      !! Snag fall to woody litter pool (gC m-2)
    REAL(r_std), DIMENSION(kjpindex)              :: wind_speed_actual             !! Actualwind speed calculated from the actual half hourly
                                                                                   !! values u and v as given in the driver (ms-1)
    REAL(r_std), DIMENSION(kjpindex,wind_days)    :: tmp_wind                      !! Temporal variable to move window for wind_ratio_sum/max_save
    REAL(r_std), DIMENSION(kjpindex,legacy_years_wind) :: tmp_wind_yr              !! Temporal variable to move window for wind_mean
    REAL(r_std), DIMENSION(kjpindex)              :: u_raw, v_raw                  !! Retrieved wind speed before transformation (ms-1)
    REAL(r_std), DIMENSION(ngrnd)                 :: zi_soil                       !! depths of intermediate soil levels (m)
    REAL(r_std), DIMENSION(0:ngrnd)               :: zf_soil                       !! depths of full soil levels (m)
    REAL(r_std), DIMENSION(kjpindex,nvm)          :: tmp                           !! dummy variable to calculate xios output
    REAL(r_std), DIMENSION(kjpindex,nvm,ncirc)    :: tmp2                          !! dummy variable to calculate xios output
    REAL(r_std), DIMENSION(kjpindex,nvm,ncirc)    :: vl_diff_daily                 !! Difference in conductivity lost since the day before.
    REAL(r_std), DIMENSION(kjpindex,nvm,ncirc)    :: vessel_mortality_daily        !! Proportion of daily vessel mortality due to cavitation in the xylem.
    REAL(r_std), DIMENSION(kjpindex)              :: est_co2flux                   !! estimate of co2flux within a single time step. The value is not exact
                                                                                   !! because several fluxes are missing (gC m-2 s-1)
    REAL(r_std), DIMENSION(kjpindex,nelements)    :: nbp_daily_flux                !! NBP at the end of stomate_lpj.f90 (gC m-2 d-1)
    REAL(r_std), DIMENSION(kjpindex,nvm)          :: fco2_flux                     !! NEE for ESM (gC m^{-2} one_day^{-1})
    REAL(r_std), DIMENSION(kjpindex)              :: nbp_intersurf                 !! NBP as it will be calculated in intersurf (gC m^{-2} one_day^{-1})
    REAL(r_std), DIMENSION(kjpindex)              :: error_count
    REAL(r_std), DIMENSION(kjpindex,nvm,12,ninput):: n_input                       !! This n_input_raw but corrected for empty PFTs and avoiding double 
                                                                                   !! counting (gN/m**2/timestep) 
    REAL(r_std), DIMENSION(kjpindex,nvm,0:nfm_types-1)  :: forest_area             !! Forest area: total and splt by managmenet strategy (km^2) 
    REAL(r_std), DIMENSION(kjpindex, ngrnd,nvm,nelements) :: rel_deepSOM_a         !! Distribution of the active pool within the profile (unitless)
    REAL(r_std), DIMENSION(kjpindex, ngrnd,nvm,nelements) :: rel_deepSOM_s         !! Distribution of the slow pool within the profile (unitless)
    REAL(r_std), DIMENSION(kjpindex, ngrnd,nvm,nelements) :: rel_deepSOM_p         !! Distribution of the passive pool within the profile (unitless)
    REAL(r_std), DIMENSION(kjpindex, ngrnd,nvm,nelements) :: deepSOM_a_stock       !! deep active SOM expressed in stock (g/m**2) 
    REAL(r_std), DIMENSION(kjpindex, ngrnd,nvm,nelements) :: deepSOM_s_stock       !! deep slow SOM expressed in stock (g/m**2)
    REAL(r_std), DIMENSION(kjpindex, ngrnd,nvm,nelements) :: deepSOM_p_stock       !! deep passive SOM expressed in stock (g/m**2)
    REAL(r_std), DIMENSION(kjpindex,nvm)          :: qc_negative_gpp               !! Count number of days that the gpp_daily is negative
    REAL(r_std), DIMENSION(kjpindex, nelements)           :: mbc_stomate_lpj       !! Mass balance closure from stomate_lpj
    INTEGER(i_std)                                        :: istm

    !_ ================================================================================================================================

    !! 1. Initialize variables

    !! 1.1 Store current time step in a common variable
    itime = kjit

    !! 1.4 Initialize first call
    resp_growth(:,:) = zero
    resp_maint(:,:) = zero
    resp_hetero(:,:) = zero
    atm_to_immob(:,:) = zero
    plant_n_uptake(:,:,:) = zero
    n_mineralisation(:,:) = zero
    fco2_lu(:) = zero
    fco2_wh(:) = zero
    fco2_ha(:) = zero
    qc_negative_gpp(:,:) = zero

    ! Daily initialization
    IF (FirstTsDay) THEN
       
       ! Reset daily variables
       est_co2flux_d(:) = zero

       ! Dynamic parameters that are calculated once per day
       CALL dynamic_parameters_calculate(kjpindex, precip_longterm)
    
    END IF

    ! Check that initialization is done
    IF (l_first_stomate) CALL ipslerr_p(3,'stomate_main','Initialization not yet done.','','')
    
    IF (printlev_loc >= 4) THEN 
       WRITE(numout,*) 'stomate_main: date=',days_since_beg, &
            ' ymds=', year_end, month_end, day_end, sec_end, &
            ' itime=', itime, ' do_slow=',do_slow
    ENDIF


    ! Yearly update of forest management variables
    IF ( FirstTsYear ) THEN

       IF (veget_update > 0) THEN
          !! In this case the model makes use of LCC, as veget_max is changing,
          !! the matching fm maps should be read to make sure that the management
          !! of newly planted forests is prescribed. Update forest_managed by
          !! reading the map for the matching year
          CALL sapiens_forestry_set_fm(kjpindex, lalo, neighbours, resolution, contfrac, forest_managed)
       END IF

       ! Species change and management change are typically used in simulations
       ! experiments that branch off from a single spinup, transient and historical
       ! simulation. The spinup, transient and historical simulation (S/T/H) are ran without
       ! these changes. Following S/T/H, the restart contains values (-999) for the desired
       ! species and the desired management. In the experiment run we want to overwrite
       ! these values by the desired values reflecting the experiment. This is done here.
       ! It would be enough to only do this at the start of the experiment. For the moment
       ! it is done every year. Note that the name of the flag "ok_change_species" is not
       ! well chosen as it controls both species and management changes. A better name could
       ! be "ok_forestry_experiments"
       IF (ok_change_species) THEN
          !! Update species change to replant by reading map or run.def
          CALL sapiens_forestry_set_species_change(kjpindex, lalo, neighbours, &
               resolution, contfrac, species_change_map)
          
          !! Update fm_change_map by reading file or run.def
          CALL sapiens_forestry_set_desired_fm(kjpindex, lalo, neighbours, &
               resolution, contfrac, forest_managed, fm_change_map)
       END IF

    END IF

    ! Guillaume M. -- These are diagnostics only, sent on the end-of-day cadence (do_slow),
    ! the same one VEGET_MAX and HARVEST_LCC_AREA_ACC use. A field sent once a year is never
    ! captured by XIOS in a monthly file and comes out entirely masked. veget_max also drifts
    ! within the year through progressive harvest, so the daily value carries information.
    IF (do_slow) THEN

       !! Calculate the global extent (km2) of the different management strategies
       ! The last two ifm strategies are related to agriculture hence
       ! nfm_types-2
       forest_area(:,:,:) = zero
       DO ivm = 2,nvm

          DO ifm = 0,nfm_types-2
             ! Calculate the area (km2) of each FM strategy (including undefined)
             WHERE (forest_managed(:,ivm) == ifm)
                forest_area(:,ivm,ifm) = veget_max(:,ivm) * area(:) * &
                     contfrac(:) / m2_to_km2
             END WHERE
          END DO

          !! Calculate the total forest area (km2) per grid cell
          ! Indices 0 to 5 are used for ifm_strategies. nfm_types-1 is
          ! used for the total forest area.
          IF (is_tree(ivm)) THEN
             forest_area(:,ivm,nfm_types-1) = veget_max(:,ivm) * area(:) * &
                  contfrac(:) / m2_to_km2
          END IF

       END DO

       ! Send forest area and area under management to XIOS
       CALL xios_orchidee_send_field("AREA_UNDEFINED", forest_area(:,:,0))
       CALL xios_orchidee_send_field("AREA_UNMANAGED", forest_area(:,:,ifm_none))
       CALL xios_orchidee_send_field("AREA_ROTATIONAL", forest_area(:,:,ifm_thin))
       CALL xios_orchidee_send_field("AREA_COPPICE", forest_area(:,:,ifm_cop))
       CALL xios_orchidee_send_field("AREA_SRC", forest_area(:,:,ifm_src))
       CALL xios_orchidee_send_field("AREA_UNEVEN", forest_area(:,:,ifm_uneven))
       CALL xios_orchidee_send_field("FOREST_AREA", forest_area(:,:,nfm_types-1))
       ! Guillaume M. -- AC12_R and AC12_F_REAL are NOT sent here: age_class_distr fills
       ! them on the last time step of the year, after this block, so a send placed here
       ! only ever ships zeros. They are sent from stomate_lpj right after that call.

    END IF


    !! 3. Special treatment for some input arrays.

    !! 3.1 Sum of liquid and solid precipitation
    precip(:) = ( precip_rain(:) + precip_snow(:) )*one_day/dt_sechiba

    !! 3.3 Adjust time step of GPP
    ! Units: gpp [gC m-2 veget_max dt_sechiba-1], gpp_d [gC m-2 day-1],
    ! and gpp_daily [gC m-2 days-1]. Note that in diffuco_trans_co2 gpp
    ! was multiplied by veget_max. That is why is has to divided by it
    ! here.
    ! No GPP for bare soil
    gpp_d(:,1) = zero
    ! GPP per PFT
    DO j = 2,nvm   
       WHERE (veget_max(:,j) > min_stomate)
          ! The PFT is available on the pixel convert to gC m-2 day-1
          gpp_d(:,j) =  gpp(:,j)/ veget_max(:,j) * one_day/dt_sechiba
       ELSEWHERE
          ! The PFT is absent on the pixel
          gpp_d(:,j) = zero
       ENDWHERE
    ENDDO

    !!! peatland
    !! Use peatland mc to calculate moisture inhibition function (Moyano et al., 2012)
    !! layer 1-11: mc_peat=mc(soiltile4)
    !! layer 12-32: mc_peat=median value of the top 11 layers. To avoid extreme dry condition in layers below 2m
    ! 
    IF (perma_peat) THEN
       DO ji=1,kjpindex
          DO j = 1, ngrnd
             IF (j .LE. nslm) THEN 
                mc_peat(ji,j) = shumdiag_peat(ji,j)
             ELSE
                mc_peat(ji,j) = shumdiag_peat(ji,nslm) 
             ENDIF
          ENDDO
       ENDDO
    ENDIF
    !!! tides
    IF (perma_peat .AND. tides ) THEN
       DO ji=1,kjpindex
          DO j = 1, ngrnd
             IF (j .LE. nslm) THEN 
                mc_man(ji,j) = shumdiag_man(ji,j)
             ELSE
                mc_man(ji,j) = shumdiag_man(ji,nslm) 
             ENDIF
          ENDDO
       ENDDO
    ENDIF
    IF (perma_peat) THEN
      DO ji=1,kjpindex
         DO jv=1,nvm
            istm = pref_soil_veg(jv)
            IF (is_peat(jv) .AND. soiltile(ji,istm) .GT. zero) THEN
               hsdeep_new(ji,:,jv)=mc_peat(ji,:)
            ELSE
               hsdeep_new(ji,:,jv)=hsdeep(ji,:,jv)
            ENDIF
         ENDDO
      ENDDO 
      
      ! agri_peat 
      IF (agri_peat) THEN
        DO j=1,ngrnd
           IF (j .LE. nslm) THEN
             hsdeep_new(:,j,15)=shumdiag_croppeat(:,j)
             hsdeep_new(:,j,16)=shumdiag_croppeat(:,j)
           ELSE
             hsdeep_new(:,j,15)=shumdiag_croppeat(:,nslm)
             hsdeep_new(:,j,16)=shumdiag_croppeat(:,nslm)
           ENDIF
        ENDDO
      ENDIF
    ENDIF
    !

    !! 4. Calculate variables for dt_stomate (i.e. "daily")

    ! Note: If dt_days /= 1, then variables 'xx_daily' (eg. half-daily or bi-daily) are by definition
    ! not expressed on a daily basis. This is not a problem but could be
    ! confusing

    !! 4.1. Calculate water stress accounting for hydraulic architecture
    !  If hydraulic architecture is used, vegstress_day (water stress) for use in stomate is 
    !  not determined from vegstress as calculated in hydrol.f90, but it is calculated as the 
    !  ratio between a proxy for stressed and unstressed ecosystem functioning. 
    !  Nightvalues are exclude. On the first day of the simulation (no biomass)
    !  the stressed and unstressed proxies are probably both equal to the initilazed value (=zero, 
    !  initialized in sechiba.f90), to avoid numerical issues in this case we set vegstress_day 
    !  to 1. If hydraulic architecture is not used vegstress_day is calculated from vegstress as 
    !  determined in hydrol.
    IF (ok_hydrol_arch) THEN
       
       IF (ok_vessel_mortality) THEN
          ! Initialize values.
          vl_diff_daily(:,:,:) = zero
          vessel_mortality_daily(:,:,:) = zero

          ! Accumulates half hourly values of vessel_loss into daily values.
          DO icir = 1,ncirc
             DO ipts = 1,kjpindex
                DO ivm = 1,nvm
                   IF (veget_max(ipts,ivm) .LE. min_stomate) THEN
                      ! Where there is no vegetation on the tile, PLC is zero.
                      vessel_loss_daily(ipts,ivm,icir) = zero
                   ELSE
                      ! Where there is vegetation, vessel_loss_daily is updated to take
                      ! the maximum daily value of vessel_loss.  
                      vessel_loss_daily(ipts,ivm,icir) = MAX(vessel_loss(ipts,ivm,icir), &
                           vessel_loss_daily(ipts,ivm,icir))
                   ENDIF
                END DO
             END DO
          END DO

       ELSE

          ! When ok_vessel_loss is not used it still need
          ! values to keep xios happy
          biomass_init_drought(:,:,:,:,:) = zero
          kill_vessels(:,:,:) = .FALSE.
          vessel_loss_previous(:,:,:) = zero
          vessel_loss_daily(:,:,:) = zero
          vl_diff_daily(:,:,:) = zero
          vessel_mortality_daily(:,:,:) = zero

       END IF   

       ! Accumulate the half hourly values into a daily value. Accumulate the
       ! stressed and unstressed proxy first and then take the average. We 
       ! first accumulate and then take the ratio because that way we better 
       ! account for the night time values. If we do it in the other order 
       ! we need to assign a value to the ratio during the night. Whether we
       ! take zero or one, this will bias our waterstress number because 
       ! the number of half hours during the night is different for the 
       ! different pixels. So although water stress could be higher in the 
       ! south than in the north. During the growing season, day are shorter
       ! in the south so if we set the ratio to 1 during the night, our daily 
       ! water stress in the south may be less than in the north because we
       ! more 1's in the daily time series.

       ! +++ CHECK +++
       ! The description below shows the following is based on an old hydr arch scheme

       ! The order of the calculation may depend on the proxy used. For example, 
       ! Accumulate transpir_supply and transpir first and then calculate ratio
       ! By doing this water stress is buffered (it is assumed that if 
       ! ::transpir_supply is larger than ::transpir at one timestep it can buffer 
       ! potential water stress in the next timestep. In reality this is not the 
       ! case: transpir_supply is a potential value not a realized one.
       ! +++++++++++++
       DO ipts = 1,kjpindex

          ! The cosine of the zenith angle, is used to 
          ! identify night values.
          IF (coszang(ipts) .LT. min_stomate) THEN

             ! Redundant because the value will never be used
             ! For all pixels
             stressed(ipts,:) = zero
             unstressed(ipts,:) = zero

             ! For PFT1 under LCC
             ! vir_stressed = zero
             ! vir_unstressed = zero

          ELSE

             ! No vegetation present so ecosystem functioning
             ! was not defined
             WHERE (veget_max(ipts,:) .LT. min_stomate)

                ! To avoid uninitialized values
                stressed(ipts,:) = zero
                unstressed(ipts,:) = zero

                ! Update the values to avoid uninitialized fields
                ! in ::stressed_daily and ::unstressed_daily
                stressed_daily(ipts,:) = stressed_daily(ipts,:) + &
                     stressed(ipts,:)
                unstressed_daily(ipts,:) = unstressed_daily(ipts,:) + &
                     unstressed(ipts,:)

                ! For PFT1 under LCC
                ! vir_stressed = zero
                ! vir_unstressed = zero

             ELSEWHERE

                ! The pixel and pft contain vegetation to calculate the stress
                stressed_daily(ipts,:) = stressed_daily(ipts,:) + &
                     stressed(ipts,:)
                unstressed_daily(ipts,:) = unstressed_daily(ipts,:) + &
                     unstressed(ipts,:)

                ! For PFT1 under LCC
                ! vir_stressed = zero
                ! vir_unstressed = zero

             ENDWHERE

          ENDIF

       ENDDO

       ! Calculate values at the end of the day
       IF (do_slow) THEN

          ! Calculates difference between current embolism and that of the
          ! day before in the variable vl_diff_daily. The current embolism is
          ! stored in vessel_loss_daily. The embolism of the day before is
          ! stored in vessel_loss_previous.
          vl_diff_daily(:,:,:) = vessel_loss_daily(:,:,:) - vessel_loss_previous(:,:,:)

          ! There are three scenarios regarding whether effect of embolism on
          ! turnover should be accounted for, and what the value of vessel
          ! mortality should be. Turnover from sapwood to heartwood is
          ! calculated in the MODULE stomate_turnover.f90.
          WHERE (vl_diff_daily(:,:,:) .GT. min_stomate)

             ! First scenario: embolism is increasing, so vl_diff_daily
             ! is positive. This means drought is ongoing. Effect of embolism on
             ! turnover should be accounted for, therefore kill_vessels is set
             ! to .TRUE. and vl_diff_daily are the vessels that started 
             ! malfunctioning in the current time step.
             kill_vessels = .TRUE.

             ! Vessel mortality is calculated as a fraction of embolism
             ! because embolized vessels do not die instantly, and
             ! might eventually recover as they get refilled. Default is 80%
             vessel_mortality_daily(:,:,:) = refilling_capa*vl_diff_daily(:,:,:)
             
          ELSEWHERE (ABS(vl_diff_daily(:,:,:)) .LT. min_stomate .AND. &
               vessel_loss_daily(:,:,:) .GT. zero) 
              
             ! Second scenario: embolism is stagnating or decreasing
             ! vl_diff_daily is null, but vessel_loss_daily is
             ! positive. This means drought is ongoing. Some of the already
             ! embolized vessels might die. Effect of embolism on turnover
             ! should be accounted for, therefore kill_vessels is set to .TRUE.
             kill_vessels(:,:,:) = .TRUE.
                      
             ! Vessel mortality is calculated as a fraction of total
             ! embolism to get a mortality value different from zero.
             vessel_mortality_daily(:,:,:) = &
                  0.01 * vessel_loss_daily(:,:,:)
       
          ELSEWHERE
             
             ! Third scenario: embolism is decreasing, which means drought
             ! is ending, or embolism is stagnating, and vessel_loss_daily
             ! is null, which means there is no drought. There is no effect
             ! of embolism on turnover, therefore kill_vessels is set to
             ! .FALSE.
             kill_vessels(:,:,:) = .FALSE.
             vessel_mortality_daily(:,:,:) = zero
             
          END WHERE
          
          ! Calculate the biomass as the start of a drought
          DO ipts=1,kjpindex
             DO ivm=1,nvm
                DO icir=1,ncirc

                   IF (.NOT. kill_vessels(ipts,ivm,icir)) THEN

                      ! Calculate the biomass as the start of a drought. This is the
                      ! Reference biomass for heartwood and sapwood. This variable
                      ! is used to avoid overestimating the sapwood mortality
                      ! since ::vessel_mortality_daily(:,ivm) is calculated as a
                      ! proportion of the model tree sapwood. Reference biomass is
                      ! recalculated whenever ::kill_vessels(:,ivm) is FALSE, that
                      ! is to say inbetween droughts.
                      biomass_init_drought(ipts,ivm,icir,:,:) = &
                           circ_class_biomass(ipts,ivm,icir,:,:)

                   ELSEIF ( SUM(SUM(biomass_init_drought(ipts,ivm,icir,:,:),1)).EQ. zero) THEN

                      ! Calculate the reference biomass also at the start of a simulation
                      ! We also calculate it when total sapwood biomass is zero so
                      ! that it already has a value on the first day of simulation.
                      biomass_init_drought(ipts,ivm,icir,:,:) = &
                           circ_class_biomass(ipts,ivm,icir,:,:)

                   END IF

                END DO
             END DO
          END DO
          
          ! Calculate the mean waterstress value at the end of each day
          WHERE (unstressed_daily(:,:) .LT. min_stomate &
               .OR. stressed_daily(:,:) .LT. min_stomate)

             ! No ecosystem function thus no data, we assume
             ! there is no water stress
             vegstress_day(:,:) = un

          ELSEWHERE

             ! Calculate water stress. If we first calculate
             ! the daily sum and then the ratio there is no need
             ! to divide by daylight. If we accumulate ratios
             ! we will have to divide by daylight. This is
             ! arbitrary. To do this properly we should calculate
             ! the turgor in the cells and calculate growth
             ! based on that see i.e. Fatichi et al 2013, New
             ! Phytologist. We use a simple numerical construct
             ! (the sqrt of the ratio) to overcome that complexity.
             vegstress_day(:,:) = stressed_daily(:,:) / &
                  unstressed_daily(:,:)

          ENDWHERE

          !---TEMP---
          IF(printlev_loc>=4)THEN 
             DO ipts=1,kjpindex
                DO ivm=1,nvm
                   IF( vegstress_day(ipts,ivm) .lt. 1.) THEN
                      WRITE(numout,*)'ivm,stressed_daily,unstressed_daily,diff,vegstress_day'

                      WRITE(numout,*) ivm, stressed_daily(ipts,ivm),unstressed_daily(ipts,ivm),&
                           & stressed_daily(ipts,ivm)-unstressed_daily(ipts,ivm), vegstress_day(ipts,ivm) 
                   ENDIF
                ENDDO
             ENDDO
          ENDIF
          !----------

          ! Set to zero to start accumulating for the next day
          stressed_daily(:,:) = zero
          unstressed_daily(:,:) = zero

       ELSE

          ! Set to a large value so that it is easy to detect problems
          ! This value should never be used. It should only be passed
          ! to other routine at the end of the day when do_slow is false
          vegstress_day(:,:) = large_value

       ENDIF !(do_slow)

    ELSE

       ! No hydrological architecture
       CALL stomate_accu (do_slow, vegstress, vegstress_day)
       
       ! When hydraulic architecture is not used, ok_vessel_loss
       ! cannot be used either so give its key variables values
       ! to keep xios and other parts of the code happy
       biomass_init_drought(:,:,:,:,:) = zero
       kill_vessels(:,:,:) = .FALSE.
       vessel_loss_previous(:,:,:) = zero
       vessel_loss_daily(:,:,:) = zero
       vl_diff_daily(:,:,:) = zero
       vessel_mortality_daily(:,:,:) = zero

    ENDIF ! (ok_hydrol_arch)

    IF (do_slow) THEN
      ! Set to zero to start accumulating for the next day
      daylight(:,:) = zero
    ELSE
      DO ipts = 1,kjpindex
        IF (coszang(ipts) .GE. min_stomate) THEN
          daylight(ipts,:) = daylight(ipts,:) + 1
        ENDIF
      ENDDO
    ENDIF

    IF (do_slow) THEN

       ! Biomass drought is a diagnostic variable to check whether all goes
       ! well. Its values were aggregated to keep the output simple. 
       CALL xios_orchidee_send_field("SAP_INIT_DROUGHT_c", &
            biomass_init_drought(:,:,:,isapabove,icarbon) + &
            biomass_init_drought(:,:,:,isapbelow,icarbon)) 
       tmp2(:,:,:) = zero
       WHERE(kill_vessels(:,:,:))
          tmp2 = un
       ENDWHERE
       CALL xios_orchidee_send_field("KILL_VESSELS",tmp2(:,:,:))
       CALL xios_orchidee_send_field("VESSEL_LOSS_PREVIOUS",vessel_loss_previous(:,:,:))
       CALL xios_orchidee_send_field("VESSEL_LOSS_DAILY",vessel_loss_daily(:,:,:))
       CALL xios_orchidee_send_field("VESSEL_LOSS",vessel_loss(:,:,:))
       CALL xios_orchidee_send_field("VL_DIFF_DAILY",vl_diff_daily(:,:,:))
       CALL xios_orchidee_send_field("VESSEL_MORTALITY_DAILY",vessel_mortality_daily(:,:,:))

       ! At the end of the day, vessel_loss_previous takes the value of
       ! vessel_loss_daily, before vessel_loss_daily is recalculated, the
       ! next day.
       vessel_loss_previous(:,:,:) = vessel_loss_daily(:,:,:)

       ! Set to zero to start accumulating again at the beginning of the
       ! next day
       vessel_loss_daily(:,:,:) = zero

    ELSE

       ! In stomate.f90 XIOS will be called 48 times per day and it will 
       ! store the average. That is not at all what we want for, e.g. 
       ! vessel_loss_daily. Therefore we will write NaN 47 times and only the 
       ! daily value we are interested in at the end of the day.
       ! It is not the cae for vessel_loss which is calculated half-hourly
 
       ! It is not the cae for vessel_loss which is calculated half-hourly
       tmp(:,:) = xios_default_val
       tmp2(:,:,:) = xios_default_val
       CALL xios_orchidee_send_field("SAP_INIT_DROUGHT_c",tmp2(:,:,:)) 
       CALL xios_orchidee_send_field("KILL_VESSELS",tmp2(:,:,:))
       CALL xios_orchidee_send_field("VESSEL_LOSS_PREVIOUS",tmp2(:,:,:))
       CALL xios_orchidee_send_field("VESSEL_LOSS_DAILY",tmp2(:,:,:))
       CALL xios_orchidee_send_field("VL_DIFF_DAILY",tmp2(:,:,:))
       CALL xios_orchidee_send_field("VESSEL_MORTALITY_DAILY",tmp2(:,:,:))
       CALL xios_orchidee_send_field("VESSEL_LOSS",vessel_loss(:,:,:))

    END IF

    ! Debug
    IF (printlev_loc>=4) THEN
       IF(do_slow) THEN
          WRITE(numout,*) 'CHECK: vegstress_day after stomate_accu',&
               vegstress_day(:,:)
       ENDIF
    ENDIF
    !-

    ! Here we add some calculations for daily max/min for windthrow module 
    ! Below processes were inside of the flag ok_windthrow before. But in case
    ! of restarting the wind simulation from simulations without wind
    ! disturbance, first five years over estimate wind damage because mean wind
    ! speed is not added to the restart. So it is outisde of flag.
    !IF (ok_windthrow) THEN

    ! Retrieve the wind speed to the raw value before transformation. 
    ! In dim2driver u,v is transformed to traslocate 10m to 2m above zhd. But because
    ! 1) this transformation decreases wind speed a lot especially when the wind speed 
    ! is high and  2) we need to calculate gustiness in the windthrow module, the raw wind 
    ! speed, (exactly speaking, 6-hourly wind speed read from forcing and interpolated
    ! to half-hourly time scale) will be retrieved.
    ! To be comparable with critical wind speed calculate in the windthrow
    ! module, U and V must be recalculated at the top of the canopy. This way
    ! the wind speed remain independant from the canopy height. Do to that we
    ! fist need to retreive the pixel average height (including all PFT exept
    ! PFT 1 ??)   

    u_raw(:) = u(:)
    v_raw(:) = v(:)
 
   ! Pixel height is no longer used here because it complicates the parameterization 
   ! of storm detection variables. The input wind speed would otherwise vary with stand 
   ! structure, deviating from the raw forcing values. 
   ! Instead, PFT-level wind speed calculations have been incorporated directly 
   ! into the windthrow module.
    IF(.NOT. use_fluxnet) THEN
           u_raw(:) = u(:) * LOG(10.0/z0m(:)) / LOG(2.0/z0m(:))
           v_raw(:) = v(:) * LOG(10.0/z0m(:)) / LOG(2.0/z0m(:))
    ENDIF

    ! Calculate the actual wind speed
    wind_speed_actual(:) = SQRT(u_raw(:)**2 + v_raw(:)**2)
    ! Detecting wind storm for the large scale is not easy. Especially since
    ! there are some areas that normaly has a high wind speed (ex. coastal area), 
    ! using the same wind speed threshold for the global can result in the continuous
    ! storm in those ares. The best solution would be developing the adjustment of
    ! trees to strong winds, (such as deeper roots ,change in crown shapes, and
    ! change in d/h relationship), bur unfortunately, there are technical limits to 
    ! achieve that. As another approach, here we calculate ratio of current wind
    ! speed over longterm average wind speed of the pixel (5yrs). If the ratio is
    ! higher than threshold for some time steps, then the storm will be detected in 
    ! the windthrow module.
    ! +++PLEASE note that this ratio is sensitive to the RESOLUTION of the forcing. 
    WHERE(wind_speed_actual(:)/ (SUM(wind_mean_save(:,:),2)/legacy_years_wind) &
              .GT. wind_ratio_threshold) 
       wind_ratio_sum(:) = wind_ratio_sum(:) + &
             wind_speed_actual(:)/(SUM(wind_mean_save(:,:),2)/legacy_years_wind)
    ENDWHERE
    IF(printlev_loc .GE. 4 .AND. ok_windthrow) THEN
       WRITE(numout,*) 'u,v entered',u(test_grid),v(test_grid)
       WRITE(numout,*) 'u,v raw',u_raw(test_grid),v_raw(test_grid)
       WRITE(numout,*) 'wind_speed_current',wind_speed_actual(test_grid)
       WRITE(numout,*) 'wind_mean', SUM(wind_mean_save(test_grid,:))/legacy_years_wind
       WRITE(numout,*) 'wind_ratio_sum',wind_ratio_sum(test_grid)
    ENDIF

    ! Accumulate the half hourly values into a daily value. Accumulate the
    ! daily wind speed first and then take the average. grnd_80 is the index 
    ! of ngrnd closest to 80 cm depth. stempdiag is discretized along diaglev
    ! (see control.f90) meaning that stempdiag gives the temperature for the
    ! node as specified in znh (see vertical_soil.f90).
    soil_max_daily(:) = MAX(stempdiag(:,grnd_80),soil_max_daily(:))

    ! Calculate daily maximum wind speed
    wind_max(:) = MAX(wind_speed_actual(:),wind_max(:))
    ! Guillaume M. -- Daily maximum of the wind ratio, companion of wind_max.
    wind_ratio_max(:) = MAX(wind_ratio_max(:), &
         wind_speed_actual(:)/(SUM(wind_mean_save(:,:),2)/legacy_years_wind))
    ! wind_sum to calculate average wind. This will be used for longterm
    ! average wind speed.
    wind_sum(:) = wind_sum(:) + wind_speed_actual(:)
    daily_wspeed_fire = zero
    IF (do_slow) THEN
       ! Moving the variables for saving wind speed
       tmp_wind(:,:) = zero
       DO iday = 1, wind_days - 1
          tmp_wind(:,iday+1) = wind_ratio_sum_save(:,iday)
       ENDDO
       tmp_wind(:,1) = wind_ratio_sum(:)
       wind_ratio_sum_save = tmp_wind

       ! Daily max wind speed is saved for several days. Since the storm can
       ! lasts several days, to calculate damage once, the module calculate
       ! start and the end of the storm. Maximum wind speed over that period
       ! will be used for damage calculation. 
       tmp_wind(:,:) = zero
       DO iday = 1, wind_days - 1
         tmp_wind(:,iday+1) = wind_speed_max_save(:,iday)
       ENDDO
       tmp_wind(:,1) = wind_max(:)
       wind_speed_max_save = tmp_wind

       ! Guillaume M. -- Same rolling-buffer shift, applied to the wind ratio maximum.
       tmp_wind(:,:) = zero
       DO iday = 1, wind_days - 1
         tmp_wind(:,iday+1) = wind_ratio_max_save(:,iday)
       ENDDO
       tmp_wind(:,1) = wind_ratio_max(:)
       wind_ratio_max_save = tmp_wind


       IF (printlev_loc .GE. 4 .AND. ok_windthrow) THEN
          WRITE(numout,*) 'ratio_sum',wind_ratio_sum(test_grid)
          WRITE(numout,*) 'save_ratio_sum',wind_ratio_sum_save(test_grid,:)
          WRITE(numout,*) 'sum',wind_sum(test_grid)
          WRITE(numout,*) 'save_max', wind_speed_max_save(test_grid,:)
       ENDIF

       IF(ts_annual_proc) THEN

          IF(printlev_loc .GE. 4) THEN
             WRITE(numout,*) 'year_mean',wind_sum(test_grid)/(one_year * one_day/dt_sechiba)
             WRITE(numout,*) 'wind_mean_all_bf',wind_mean_save(test_grid,:)
          ENDIF
          ! 
          tmp_wind_yr(:,:) = zero
          DO iyr = 1, legacy_years_wind - 1
             tmp_wind_yr(:,iyr+1) = wind_mean_save(:,iyr)
          ENDDO             
          tmp_wind_yr(:,1) = wind_sum(:)/(one_year * one_day /dt_sechiba)
          wind_mean_save = tmp_wind_yr
          wind_sum(:) = zero
       ENDIF ! ts_annual_proc

       ! wind_max_daily/soil_temp_daily will be passed to the wind_damage
       ! module. wind_max_daily/soil_max_daily reset to zero so they can be 
       ! used again for the next day
       wind_max_daily(:) = wind_max(:)
       wind_max(:) = zero
       wind_ratio_sum(:) = zero
       wind_ratio_max(:) = zero  ! JJ 2026: reset daily accumulator
       soil_temp_daily(:) = soil_max_daily(:)
       soil_max_daily(:) = zero

    ENDIF !do_slow
    !ENDIF

      DO ipts = 1,kjpindex
        IF (coszang(ipts) .GE. min_stomate) THEN
          ! 3.4 Calculate photoperiod and beetle diapause status
          DO ivm=2,nvm
            IF (daylight(ipts,ivm)*0.5 > diapause_thres_daylength(ivm)) THEN
              beetle_diapause(ipts,ivm) = 1
            ENDIF
          ENDDO
        ENDIF
      ENDDO

    !! Calculate the light that reaches each canopy layer
    !  Compute seasonal daytime transmitted light to canopy levels 
    !  This quantity is used to calculate how much recruitment can occur
    !  underneath the canopy. Recruitment is simulated in stomate_prescribe.f90
    DO ipts=1,kjpindex

       DO ivm=1,nvm

          ! If we are at a daytime time step and there is growth (gpp>0) then
          ! accumulate instantaneous transmitted light to each canopy level
          IF ( coszang(ipts) .GT. min_stomate .AND. &
               gpp(ipts,ivm) .GT. min_stomate .AND. &
               Pgap_cumul(ipts,ivm,1).NE.zero) THEN

             ! Use the light that was transmitted through all the layers 
             ! to reach the forest floor (= level 1)
             daylight_count(ipts,ivm) = daylight_count(ipts,ivm) + 1
             light_tran_to_floor_season(ipts,ivm) = &
                  light_tran_to_floor_season(ipts,ivm) + &
                  Pgap_cumul(ipts,ivm,1)

             ! Debug
             IF (printlev_loc.GT.4 .AND. ivm.EQ. test_pft)THEN
                WRITE(numout,*) 'It is daytime and growth occurs, '&
                     &'daylight_count(ipts,ivm)= ', daylight_count(ipts,ivm)
                WRITE(numout,*) 'gpp(ipts,ivm)= ', gpp(ipts,ivm)
                WRITE(numout,*) 'Absolute transmitted light '&
                     &'to the forest floor ',ivm,ipts,&
                     light_tran_to_floor_season(ipts,ivm)
             ENDIF
             !-

          ENDIF   ! daytime and growing season

       ENDDO

       ! Calculate average transmitted light at the end of the year
       ! Note that ts_annual_proc is true if current time-step corresponds to the moment when the annual
       ! processes should be done. By default, this corresponds to the last time-step of the year.
       IF (ts_annual_proc) THEN

          ! Calculate average transmitted light at the end of the year
          ! Only account for days during the growing season. Note that
          ! light_tran_to_floor_season is only a correctly calculated
          ! the last day of the year. That is OK because recruitment is
          ! calculated only the last day of the year as well.
          DO ivm=2,nvm

             IF (daylight_count(ipts,ivm) .GT. zero) THEN

                light_tran_to_floor_season(ipts,ivm) = &
                     light_tran_to_floor_season(ipts,ivm) / &
                     daylight_count(ipts,ivm)

             ELSE

                IF (veget_max(ipts,ivm).GT.min_stomate) THEN
                   ! There was no GPP during this year in this PFT. So,
                   ! all the light is transmitted to the ground level.
                   light_tran_to_floor_season(ipts,ivm) = un
                ELSE
                   ! There is no vegetation so the variable is not defined
                   light_tran_to_floor_season(ipts,ivm) = zero
                END IF

             ENDIF

             ! Debug
             IF(printlev_loc>=4)THEN
                WRITE(numout,*) 'DEBUG in stomate.f90 it is end of the ', &
                     'year in stomate.f90'
                WRITE(numout,*) 'daylight_count to divide by here is, ', &
                     daylight_count(ipts,ivm)
                WRITE(numout,*) 'transmitted light, ',&
                     light_tran_to_floor_season(ipts,ivm)
             ENDIF
             !- 
          ENDDO

          ! Reset the counter for the next year. Note that
          ! light_tran_to_floor_season will be reset after
          ! it was send to XIOS and after it was used in 
          ! stomate_prsecribe.f90 for calculating recruitment
          daylight_count(ipts,:) = zero

       ENDIF ! ts_annual_proc

    ENDDO ! ipts=1,kjpindex
    

    IF (.NOT.ts_annual_proc) THEN
       ! The correct value can only be calculated at the end
       ! of the year. Send an NaN to XIOS.
       tmp(:,:) = xios_default_val
       CALL xios_orchidee_send_field("LIGHT_TRAN_SEASON",tmp(:,:))
    ELSE
       ! At the last day of the year we are sending the correct
       ! value to XIOS.
       CALL xios_orchidee_send_field("LIGHT_TRAN_SEASON",light_tran_to_floor_season(:,:))
    ENDIF
    
    !! 4.1 Accumulate instantaneous variables (do_slow=.FALSE.) 
    ! Accumulate instantaneous variables (do_slow=.FALSE.) and eventually 
    ! calculate daily mean value (do_slow=.TRUE.)
    CALL stomate_accu(do_slow, litterhumdiag, litterhum_daily)
    CALL stomate_accu(do_slow, temp_air,      t2m_daily)
    CALL stomate_accu(do_slow, VPD,           vpd_daily_mean)
    CALL stomate_accu(do_slow, temp_sol,      tsurf_daily)
    CALL stomate_accu(do_slow, stempdiag,     tsoil_daily)
    CALL stomate_accu(do_slow, precip,        precip_daily)
    CALL stomate_accu(do_slow, gpp_d,         gpp_daily)
    CALL stomate_accu(do_slow, drainage_pft, drainage_daily) 
    CALL stomate_accu(do_slow, tdeep, tdeep_daily)
    CALL stomate_accu(do_slow, hsdeep, hsdeep_daily)
    CALL stomate_accu(do_slow, decomp_rate,decomp_rate_daily)
    CALL stomate_accu(do_slow, snow, snow_daily)
    CALL stomate_accu(do_slow, pb * 100., pb_pa_daily)
    CALL stomate_accu(do_slow, temp_sol, temp_sol_daily)
    CALL stomate_accu(do_slow, snowdz, snowdz_daily)
    CALL stomate_accu(do_slow, snowrho, snowrho_daily)
    CALL stomate_accu(do_slow, wind_speed_actual,daily_wspeed_fire)

    ! Monitor gpp_daily
    IF (do_slow) THEN
       WHERE (gpp_daily(:,:) .LT. zero)
          ! Quality check. Record for which pixels and PFTs
          ! the daily gpp is negative. XIOS will accumulate over a year.
          qc_negative_gpp(:,:) = un
       END WHERE
    END IF
    CALL xios_orchidee_send_field("QC_NEGATIVE_GPP",qc_negative_gpp)
    
    !! 4.2 Daily minimum temperature
    t2m_min_daily(:) = MIN( temp_air(:), t2m_min_daily(:) )    
    t2m_max_daily(:) = MAX( temp_air(:), t2m_max_daily(:) )

    !! Daily maximum VPD
    vpd_daily_max(:) = MAX( VPD, vpd_daily_max(:) )
    
    !! 4.3 Calculate maintenance respiration
    ! Note: lai is passed as output argument to overcome previous problems with 
    ! natural and agricultural vegetation types.
    CALL maint_respiration &
         & (kjpindex, temp_air, t2m_longterm, stempdiag, root_profile, &
         & circ_class_n, circ_class_biomass,resp_maint_part_radia, cn_leaf_init_2D)

    ! Maintenance respiration separated by plant parts
    resp_maint_part(:,:,:) = resp_maint_part(:,:,:) &
         & + resp_maint_part_radia(:,:,:)
     ! Aggregate maintenance respiration across the different plant parts
    resp_maint_radia(:,:) = zero
    DO ivm=2,nvm
       DO k= 1, nparts
          resp_maint_radia(:,ivm) = resp_maint_radia(:,ivm) &
               & + resp_maint_part_radia(:,ivm,k)
       ENDDO
    ENDDO

    !! Calculate how much bm will be added to the litter during this
    !  time step. Needs to be done before the mass balance check because
    !  the values calculated here have to be used in the mass balance 
    ! check
    !  Including: litter update, lignin content, PFT parts, litter decay,
    !  litter heterotrophic respiration, dead leaf soil cover.
    !  Note: there is no vertical discretisation in the soil for litter decay.
    n_mineralisation(:,:) = zero
    IF (do_slow) THEN
       ! Use the residual to achieve a higher precision of the calculations
       turnover_littercalc(:,:,:,:) = turnover_resid(:,:,:,:)
       bm_to_littercalc(:,:,:,:) = bm_to_litter_resid(:,:,:,:)
       tree_bm_to_littercalc(:,:,:,:) = tree_bm_to_litter_resid(:,:,:,:)
    ELSE
       ! Use 1/48th of the daily turnover and bm_to_litter.
       turnover_littercalc(:,:,:,:) = turnover_daily(:,:,:,:) * dt_sechiba/one_day
       bm_to_littercalc(:,:,:,:) = bm_to_litter(:,:,:,:) * dt_sechiba/one_day
       tree_bm_to_littercalc(:,:,:,:) = tree_bm_to_litter(:,:,:,:) * dt_sechiba/one_day   
    ENDIF
    
    !! 4.4 Initialize check for mass balance closure
    !  Mass balance closure for the half-hourly (dt_sechiba)
    !  processes in stomate.f90. This test is always performed.
    !  If err_act.EQ.1 then the value of the mass balance error
    !  -if any- is written to the history file. 
    check_intern(:,:,:,:) = zero
    pool_start(:,:,:) = zero
    DO iele = 1,nelements

       ! Biomass pool (gC m-2)*(m2 m-2). 
       ! Note that we only check where the bm_to_litter and turnover_daily
       ! are going to be processed during this time step. With every time step the
       ! litter pool will increase but the values of turnover_daily and 
       ! bm_to_litter remain constant in stomate.lpj. The values of 
       ! bm_to_litter_resid and turnover)resid are changing with every time
       ! step.  
       DO ipar = 1,nparts
          pool_start(:,:,iele) = pool_start(:,:,iele) + &
               (turnover_littercalc(:,:,ipar,iele) + &
               bm_to_littercalc(:,:,ipar,iele)) * veget_max(:,:)
       ENDDO

       ! Litter pool (gC m-2)*(m2 m-2) 
       DO ilitt = 1,nlitt
          DO ilev = 1,nlevs
             pool_start(:,:,iele) = pool_start(:,:,iele) + &
                  litter(:,ilitt,:,ilev,iele) * veget_max(:,:)
          ENDDO
       ENDDO

       IF (ok_soil_carbon_discretization) THEN
          ! Define the soil layers
          zf_soil(:) = zero
          zi_soil(:) = zero
          zi_soil(:) = znt(:)
          zf_soil(1:ngrnd) = zlt(:)
          zf_soil(0) = 0.
          ! Soil carbon (gC m-3) * (m2 m-2)
          DO igrn = 1,ngrnd
             pool_start(:,:,iele) = pool_start(:,:,iele) + &
                  (deepSOM_a(:,igrn,:,iele) + deepSOM_s(:,igrn,:,iele) + &
                  deepSOM_p(:,igrn,:,iele)) * &
                  (zf_soil(igrn)-zf_soil(igrn-1)) * veget_max(:,:)
          END DO
       ELSE
          ! Soil carbon (gC m-2) *  (m2 m-2)
          DO icarb = 1,ncarb
             pool_start(:,:,iele) = pool_start(:,:,iele) + &
                  som(:,icarb,:,iele) * veget_max(:,:)
          ENDDO
       ENDIF

       DO ivm = 1, nvm
          pool_start(:,ivm,iele) = pool_start(:,ivm,iele) + &
               SUM(SUM(harvest_pool_acc(:,ivm,:,iele,:),3),2)/(area(:)*contfrac(:))
       END DO

    ENDDO ! # nelements

    ! Account for the N-pool in the soil 
    DO inspec = 1,nnspec
       pool_start(:,:,initrogen) = pool_start(:,:,initrogen) + &
            soil_n_min(:,:,inspec) * veget_max(:,:)
    ENDDO

    ! Atm_to_immob is in gC m-2 dt_sechiba-1
    check_intern(:,:,iatm2land,initrogen) = - un * &
         atm_to_immob(:,:) * veget_max(:,:) * dt_sechiba/one_day
    
    ! Initialize check for area conservation
    veget_max_begin(:,:) = veget_max(:,:)

    ! The n_input we get from slowproc is the raw data in the sense that there is
    ! input for PFTs without veget_max. If the model is ran without updating the
    ! n_input (for example during the spinup), a copy of the raw data is needed
    ! in the case the model makes use of age classes. The same applied to land cover 
    ! change but land cover change is typically not used during the spinup. The
    ! key problem is that with age classes new PFTs are created so we need n_inputs
    ! for those PFTs. Checking the mass balance for n in stomate and stomate_lpj
    ! is easier if a clean n_input variable could be used. In this case clean means
    ! that there are only values for PFTs that have a veget_max. Both needs can be
    ! combined by using different variables in slowproc and stomate and its 
    ! subroutines.
    n_input(:,:,:,:) = n_input_raw(:,:,:,:)

    CALL littercalc (kjpindex, &
         turnover_littercalc, bm_to_littercalc, tree_bm_to_littercalc, &
         veget_max, temp_sol, stempdiag, shumdiag, litterhumdiag, som, &
         clay, silt, soil_n_min, n_input, month, harvest_pool_acc, litter, dead_leaves, &
         lignin_struc, lignin_wood, lignin_snag, n_mineralisation, deadleaf_cover, resp_hetero_litter, &
         litterfuel, som_input_inst, control_temp_inst, control_moist_inst, &
         matrixA, vectorB, CN_target, CN_som_litter_longterm, tau_CN_longterm, &
         circ_class_biomass, circ_class_n, tsoil_decomp, snag_to_wood_flux, &
         shumdiag_peat,mc_peat_above,shumdiag_croppeat, mc_croppeat_above, &
         shumdiag_man,mc_man_above, soiltile)

    ! Calculate the residuals which are then used at the last time step. This
    ! approach helps to get a higher precision in the consistency cross-checks
    ! for nbp.
    IF (do_slow) THEN
       ! At the last time step the residual should be exactly zero as calculated
       ! above
       turnover_resid(:,:,:,:) = zero
       bm_to_litter_resid(:,:,:,:) = zero
       tree_bm_to_litter_resid(:,:,:,:) = zero
    ELSE
       ! Update the remaining turnover and bm_to_litter
       turnover_resid(:,:,:,:) = turnover_resid(:,:,:,:) - &
            turnover_littercalc(:,:,:,:)
       bm_to_litter_resid(:,:,:,:) = bm_to_litter_resid(:,:,:,:) - &
            bm_to_littercalc(:,:,:,:)
       tree_bm_to_litter_resid(:,:,:,:) = tree_bm_to_litter_resid(:,:,:,:) - &
            tree_bm_to_littercalc(:,:,:,:) 
    ENDIF

    ! Heterothropic litter respiration during time step ::dt_sechiba 
    ! @tex $(gC m^{-2})$ @endtex
    resp_hetero_litter(:,:) = resp_hetero_litter(:,:) * dt_sechiba/one_day

    IF ( ok_soil_carbon_discretization ) THEN
       ! BG 201902: I commented the following line since it seems that in MICT pb is
       ! converted by error in hpa whereas it is already in hpa
       !          pb_pa = pb * 100.

       !permafrost:  get the residence time for soil carbon
       IF ( printlev>=3 ) WRITE(*,*) 'cdk debug stomate: prep to calc fbact'
       tdeep_celsius(:,:,:) = 0
       tdeep_celsius = tdeep - ZeroCelsius
       fbact = stomate_soil_carbon_discretization_microactem ( &
            tdeep_celsius, frozen_respiration_func, hsdeep, kjpindex, ngrnd, nvm, znt,&
            mc_peat, soiltile) 
       decomp_rate = 1./fbact
       heat_Zimov = zero
       ! should input daily-averaged values here
       !temp_sol -> tsurf daily, tdeep, hsdeep, stempdiag, shumdiag,
       !profil_froz_diag, snow, pb_pa...
       
       CALL stomate_soil_carbon_discretization_deep_somcycle(kjpindex, index, itime, &
            dt_sechiba, lalo, clay, silt, temp_sol, tdeep, hsdeep, snow, heat_Zimov, pb, &
            sfluxCH4_deep, sfluxCO2_deep, deepSOM_a, deepSOM_s, deepSOM_p, soil_n_min, O2_soil, &
            CH4_soil, O2_snow, CH4_snow, depth_organic_soil, som_input_inst, &
            veget_max, altmax, altmax_ind, altmax_lastyear, altmax_ind_lastyear, &
            som, som_surf, resp_hetero_soil, &
            fbact, CN_target, fixed_cryoturbation_depth, snowdz, snowrho, &
            n_mineralisation, root_depth, matrixA, &
            CN_som_litter_longterm, tau_CN_longterm, &
            deepSOM_peat, soiltile)

       resp_hetero_soil(:,:) = resp_hetero_soil(:,:) * dt_sechiba/one_day

       ! Total heterothrophic respiration during time step ::dt_sechiba 
       ! @tex $(gC m^{-2})$ @endtex
       resp_hetero_radia(:,:) = resp_hetero_litter(:,:) + resp_hetero_soil(:,:)
       resp_hetero_d(:,:) = resp_hetero_d(:,:) + resp_hetero_radia(:,:)
       resp_hetero_litter_d(:,:) = resp_hetero_litter_d(:,:) + resp_hetero_litter(:,:)
       resp_hetero_soil_d(:,:) = resp_hetero_soil_d(:,:) + resp_hetero_soil(:,:)

       ! Sum heterotrophic and autotrophic respiration in soil
       resp_total_soil(:,:) = resp_hetero_radia(:,:) + & 
            resp_maint_part_radia(:,:,isapbelow) + resp_maint_part_radia(:,:,iroot)

       som_total(:,:,:,:) = deepSOM_a(:,:,:,:) + deepSOM_s(:,:,:,:) + deepSOM_p(:,:,:,:)
       som(:,isurface,:,:) = zero
       som(:,iactive,:,:) = SUM(deepSOM_a(:,:,:,:),2)
       som(:,islow,:,:) = SUM(deepSOM_s(:,:,:,:),2)
       som(:,ipassive,:,:) = SUM(deepSOM_p(:,:,:,:),2)

       rel_deepSOM_a(:,:,:,:) =zero
       rel_deepSOM_s(:,:,:,:) =zero
       rel_deepSOM_p(:,:,:,:) =zero

       DO igrn = 1,ngrnd
          deepSOM_a_stock(:,igrn,:,:) = deepSOM_a(:,igrn,:,:) * &
               (zf_soil(igrn)-zf_soil(igrn-1))
          deepSOM_s_stock(:,igrn,:,:) = deepSOM_s(:,igrn,:,:) * &
               (zf_soil(igrn)-zf_soil(igrn-1)) 
          deepSOM_p_stock(:,igrn,:,:) = deepSOM_p(:,igrn,:,:) * &
               (zf_soil(igrn)-zf_soil(igrn-1)) 
       END DO

       DO iele = 1,nelements
          DO ivm = 1, nvm
             DO igrn = 1,ngrnd
                DO ipts = 1, kjpindex 
                   IF ((veget_max(ipts,ivm).GT.min_stomate) .AND. &
                       (SUM(deepSOM_a_stock(ipts,:,ivm,iele),1).GT.min_stomate) .AND. &
                       (SUM(deepSOM_s_stock(ipts,:,ivm,iele),1).GT.min_stomate) .AND. &
                       (SUM(deepSOM_p_stock(ipts,:,ivm,iele),1).GT.min_stomate)) THEN
                      rel_deepSOM_a(ipts,igrn,ivm,iele)=deepSOM_a_stock(ipts,igrn,ivm,iele)/SUM(deepSOM_a_stock(ipts,:,ivm,iele),1)
                      rel_deepSOM_s(ipts,igrn,ivm,iele)=deepSOM_s_stock(ipts,igrn,ivm,iele)/SUM(deepSOM_s_stock(ipts,:,ivm,iele),1)
                      rel_deepSOM_p(ipts,igrn,ivm,iele)=deepSOM_p_stock(ipts,igrn,ivm,iele)/SUM(deepSOM_p_stock(ipts,:,ivm,iele),1)
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDDO

    ELSE

       !! 4.5 Soil carbon dynamics and soil heterotrophic respiration
       ! Note: there is no vertical discretisation in the soil for litter decay.
       CALL som_dynamics (kjpindex, lalo, clay, silt, veget_max,&
            som_input_inst, control_temp_inst, control_moist_inst, drainage_pft,&
            CN_target, som, soil_n_min, resp_hetero_soil, matrixA, &
            n_mineralisation, CN_som_litter_longterm, tau_CN_longterm, &
            height_acro, height_cato, carbon_acro, carbon_cato, &
            tcarbon_acro, tcarbon_cato, resp_acro_oxic,&
            resp_acro_anoxic, resp_cato, acro_to_cato, litter_to_acro, wtp)
 
       ! The computation of resp_acro_oxic_d etc by using the sum is removed from the trunk, when ok_peat_NoDiscretisation T.
       ! such computation can be realised by using xios later if need.

       ! Initialize variables for soil carbon discretization
       som_surf(:,:,:,:) = som(:,:,:,:)
       som_total(:,:,:,:) = zero
       heat_Zimov = zero    

       ! Heterothropic soil respiration during time step ::dt_sechiba 
       ! @tex $(gC m^{-2})$ @endtex 
       resp_hetero_soil(:,:) = resp_hetero_soil(:,:) * dt_sechiba/one_day
        
       ! Total heterothrophic respiration during time step ::dt_sechiba 
       ! @tex $(gC m^{-2})$ @endtex
       resp_hetero_radia(:,:) = resp_hetero_litter(:,:) + resp_hetero_soil(:,:)
       resp_hetero_d(:,:) = resp_hetero_d(:,:) + resp_hetero_radia(:,:)
       resp_hetero_litter_d(:,:) = resp_hetero_litter_d(:,:) + resp_hetero_litter(:,:)
       resp_hetero_soil_d(:,:) = resp_hetero_soil_d(:,:) + resp_hetero_soil(:,:)       
       
       ! Sum heterotrophic and autotrophic respiration in soil. Note that this is an
       ! estimate because in stomate resp_maint is recalculated while accounting
       ! for the available C. Not all respiration estimated in resp_maint_part_radia
       ! will always happen.
       resp_total_soil(:,:) = resp_hetero_radia(:,:) + & 
            resp_maint_part_radia(:,:,isapbelow) + resp_maint_part_radia(:,:,iroot)
              
       IF (printlev>=3) WRITE (numout,*) '4.5'
       IF (printlev>=3) WRITE (numout,*) 'resp_hetero_litter(test_grid,test_pft):', &
            resp_hetero_litter(test_grid,test_pft)
       IF (printlev>=3) WRITE (numout,*) 'resp_hetero_soil(test_grid,test_pft):', &
            resp_hetero_soil(test_grid,test_pft)
       IF (printlev>=3) WRITE (numout,*) 'resp_maint_part_radia(test_grid,test_pft,isapbelow):', &
            resp_maint_part_radia(test_grid,test_pft,isapbelow)
       IF (printlev>=3) WRITE (numout,*) 'resp_maint_part_radia(test_grid,test_pft,iroot):', &
            resp_maint_part_radia(test_grid,test_pft,iroot)
       
    ENDIF ! End of if (ok_soil_carbon_discretization)

    IF (ok_ncycle) THEN
       
       CALL nitrogen_dynamics(kjpindex, njsc, clay, MAX(zero, un - silt - clay), & 
            tsoil_decomp, tmc_pft, drainage_pft, runoff_pft, swc_pft, veget_max, &
            resp_total_soil, som, &  
            n_input, month, soil_ph, n_mineralisation, pb, &  
            plant_n_uptake, bulk, soil_n_min, p_O2, bact, atm_to_immob, &
            leaching, emission, ld_redistribute, circ_class_biomass, &
            circ_class_n, cn_leaf_min_2D, cn_leaf_max_2D, cn_leaf_init_2D, &
            mcs_hydrol, mcfc_hydrol, croot_longterm, n_reserve_longterm, &
            sugar_load)  
      
    ENDIF

    ! Accumulate over the day
    plant_n_uptake_daily(:,:,:) = plant_n_uptake_daily(:,:,:) + plant_n_uptake(:,:,:)
    atm_to_immob_daily(:,:) = atm_to_immob_daily(:,:) + atm_to_immob(:,:)
    emission_daily(:,:,:) = emission_daily(:,:,:) + emission(:,:,:)
    leaching_daily(:,:,:) = leaching_daily(:,:,:) + leaching(:,:,:)
    n_input_daily(:,:,:) = n_input_daily(:,:,:) + n_input(:,:,month,:)

    IF (ok_peat_NoDiscretisation) THEN
       CALL stomate_accu (do_slow, wtp, wtp_daily)
    ENDIF
    !
    !! 4.7 Accumulate instantaneous variables (do_slow=.FALSE.) 
    ! Accumulate instantaneous variables (do_slow=.FALSE.) and eventually 
    ! calculate daily mean value (do_slow=.TRUE.) 
    CALL stomate_accu (do_slow, som_input_inst, som_input_daily)
    CALL stomate_accu (do_slow, n_mineralisation, n_mineralisation_d)

    !! 4.8 Check numerical consistency of this routine
    !  These checks only check the processes that happen
    !  every half-hour (dt_radia). This test is always 
    !  performed. If err_act.EQ.1 then the value of the 
    !  mass balance error -if any- is written to the 
    !  history file.

    !  Check surface area
    CALL check_vegetation_area("stomate dt_sechiba", kjpindex, veget_max_begin, &
         veget_max,'pixel')

    ! 4.8.2 Mass balance closure (dt_radia)
    ! Calculate final carbon and nitrogen pools
    pool_end(:,:,:) = zero
    DO iele = 1,nelements

       ! Litter pool
       DO ilitt = 1,nlitt
          DO ilev = 1,nlevs
             pool_end(:,:,iele) = pool_end(:,:,iele) + &
                  litter(:,ilitt,:,ilev,iele) * veget_max(:,:)
          ENDDO
       ENDDO

       IF (ok_soil_carbon_discretization) THEN
          ! Soil carbon
          DO igrn = 1,ngrnd
             pool_end(:,:,iele) = pool_end(:,:,iele) + &
                  (deepSOM_a(:,igrn,:,iele) + deepSOM_s(:,igrn,:,iele) + &
                  deepSOM_p(:,igrn,:,iele)) * &
                  (zf_soil(igrn)-zf_soil(igrn-1)) * veget_max(:,:)
          END DO
       ELSE
          ! Soil carbon (gC m-2) *  (m2 m-2)
          DO icarb = 1,ncarb
             pool_end(:,:,iele) = pool_end(:,:,iele) + &
                  som(:,icarb,:,iele) * veget_max(:,:)
          ENDDO
       ENDIF
       
       DO ivm = 1, nvm
          pool_end(:,ivm,iele) = pool_end(:,ivm,iele) + &
               SUM(SUM(harvest_pool_acc(:,ivm,:,iele,:),3),2)/&
               (area(:)*contfrac(:))
       END DO
    ENDDO ! # nelements

    ! The nitrogen pool in the soil may have changed
    DO inspec = 1,nnspec
       pool_end(:,:,initrogen) = pool_end(:,:,initrogen) + &
            soil_n_min(:,:,inspec) * veget_max(:,:)
    ENDDO

    ! Calculate mass balance  
    DO iele = 1,nelements
       check_intern(:,:,ipoolchange,iele) = &
            -un * (pool_end(:,:,iele) - pool_start(:,:,iele))
    ENDDO

    check_intern(:,:,iatm2land,initrogen) = check_intern(:,:,iatm2land,initrogen) + &
         atm_to_immob(:,:) * veget_max(:,:) * dt_sechiba/one_day
    
    check_intern(:,:,iland2atm,icarbon) = -un * (resp_hetero_litter(:,:) + &
         resp_hetero_soil(:,:)) * veget_max(:,:)

    DO ininput = 1,ninput
       check_intern(:,:,iatm2land,initrogen) = &
            check_intern(:,:,iatm2land,initrogen) + &
            n_input(:,:,month,ininput)*dt_sechiba/one_day * veget_max(:,:)
    ENDDO

    DO inspec= 1, nnspec
       check_intern(:,:,iland2atm,initrogen) = &
            check_intern(:,:,iland2atm,initrogen) &
            -un * (emission(:,:,inspec) * veget_max(:,:))
    ENDDO

    DO inionspec = 1, nionspec
       check_intern(:,:,ilat2out,initrogen) = &
            check_intern(:,:,ilat2out,initrogen) &
            -un * ( plant_n_uptake(:,:,inionspec) + &
            leaching(:,:,inionspec) ) * veget_max(:,:)
    ENDDO

    closure_intern(:,:,:) = zero
    DO imbc = 1,nmbcomp
       DO iele = 1,nelements
          ! Debug
          IF (printlev_loc>=4) WRITE(numout,*) &
               'check_intern, imbc, iele, ', imbc, &
               iele, SUM(check_intern(:,:,imbc,iele),2)
          !-
          closure_intern(:,:,iele) = closure_intern(:,:,iele) + &
               check_intern(:,:,imbc,iele)
       ENDDO
    ENDDO

    CALL check_mass_balance("stomate dt_sechiba", closure_intern, kjpindex, &
         pool_end, pool_start, veget_max, 'pft')

    ! Cumulate mass balance closure over the day at each time step
    mbc_daily(:,:) = mbc_daily(:,:) + SUM(closure_intern(:,:,:),2)
    
    
    !! 5. Daily processes - performed at the end of the day
    IF (do_slow) THEN
       !+++CHECK+++
       ! No longer needed. lai is no longer passed
       ! circ_class_biomass and circ_class_n are now the
       ! prognostic variables.
!!$       !! 5.1 Update lai
!!$       ! Use lai from stomate
!!$       ! ?? check if this is the only time ok_pheno is used??
!!$       ! ?? Looks like it is the only time. But this variables probably is defined 
!!$       ! in stomate_constants or something, in which case, it is difficult to track.
!!$       IF (ok_pheno) THEN
!!$          !! 5.1.1 Update LAI 
!!$          ! Set lai of bare soil to zero
!!$          lai(:,ibare_sechiba) = zero
!!$          ! lai for all PFTs
!!$          DO ipts = 1, kjpindex
!!$             DO j = 2, nvm
!!$                lai(ipts,j) = cc_to_lai(circ_class_biomass(ipts,j,:,ileaf,icarbon),&
!!$                     circ_class_n(ipts,j,:),j)
!!$             ENDDO
!!$          ENDDO
!!$          frac_age(:,:,:) = leaf_frac(:,:,:)
!!$       ELSE 
!!$          ! 5.1.2 Use a prescribed lai
!!$          ! WARNING: code in setlai is effectively the same as the lines above
!!$          ! Update subroutine if LAI should be prescribed. This is a bit an optimistic
!!$          ! function. It is rather difficult to force an lai with a dynamic allocation 
!!$          ! and a dynamic nitrogen cycle. This will quickly result in inconsistencies. 
!!$          CALL  setlai(kjpindex, lai, circ_class_biomass,circ_class_n) 
!!$          frac_age(:,:,:) = zero
!!$       ENDIF
       !++++++++++++
       !! 5.2 Calculate long-term "meteorological" and biological parameters
       ! mainly in support of calculating phenology. If ts_annual_proc=.TRUE.
       ! annual values are update (i.e. last time-step of the year by default).
       CALL season_pre_disturbance &
            &          (kjpindex, dt_days, &
            &           veget, veget_max, &
            &           vegstress_day, t2m_daily, tsoil_daily, lalo, &
            &           precip_daily, npp_daily, circ_class_biomass, circ_class_n, &
            &           turnover_daily, gpp_daily, when_growthinit, &
            &           SUM(resp_maint_part,3), resp_maint_week, &
            &           maxvegstress_lastyear, maxvegstress_thisyear, &
            &           minvegstress_lastyear, minvegstress_thisyear, &
            &           maxgppweek_lastyear, maxgppweek_thisyear, &
            &           gdd0_lastyear, gdd0_thisyear, &
            &           precip_lastyear, precip_thisyear, precip_longterm, &
            &           lm_lastyearmax, lm_thisyearmax, &
            &           maxfpc_lastyear, maxfpc_thisyear, &
            &           vegstress_month, vegstress_week, t2m_longterm, tau_longterm, &
            &           vpd_daily_mean, vpd_daily_max, vpd_mean_week, vpd_max_week, &
            &           t2m_month, t2m_week, tsoil_month, precip_month, &
            &           npp_longterm, croot_longterm, turnover_longterm, gpp_week, &
            &           gpp_year, gpp_decade, plant_status, &
            &           gdd_m5_dormance, gdd_midwinter, ncd_dormance, ngd_minus5, &
            &           time_hum_min, hum_min_dormance, gdd_init_date, &
            &           gdd_from_growthinit, herbivores, &
            &           Tseason, Tseason_length, Tseason_tmp, &
            &           Tmin_spring_time, t2m_min_daily,t2m_max_daily, &
            &           cn_leaf_min_season,nstress_season, &
            &           vegstress_season,rue_longterm, cn_leaf_init_2D, &
            &           litter, leaf_age_crit, leaf_classes,&
            &           wood_leftover_legacy, season_drought_legacy,i_beetles_generation, &
            &           woody_litter_to_use, sumTeff, beetle_diapause,ni_acc)

       !! 5.4 Waterstress

       !  The waterstress factor varies between 0.1 and 1 and is calculated
       !  from ::vegstress_season. The latter is only used in the allometric 
       !  allocation and its time integral is determined by longevity_sap for trees
       !  (see constantes_mtc.f90 for longevity_sap and see pft_constantes.f90 for 
       !  the definition of tau_hum_growingseason). The time integral for 
       !  grasses and crops is a prescribed constant (see constantes.f90). By
       !  having ::wstress_fac working on the turnover, water stress is progated
       !  into LF.
       !  Because the calculated values for ::wstress_fac are too low for this purpose
       !  Sonke Zhaele multiply it by two in the N-branch. This approach maintains
       !  the physiological basis of KF while combining it with a simple 
       !  multiplicative factor for water stress. Clearly after multiplication with 
       !  2, wstress is closer to 1 and will thus result in a KF values closer to 
       !  the physiologically expected KF. 
       !  In this implementation we take the sqrt (this is done in stomate where 
       !  ::vegstress_day is calculated from ::stressed and ::unstressed. The
       !  transformation from the ratio between stressed an unstressed gpp into a 
       !  numerical value that is used in the allocation and turnover is arbitrairy.
       !  A more physiological approach accounting for turgor would be needed to de
       !  fundamentally better.
       !  Note that the current implementation allows for the plants to adapt to drought 
       !  by adjusting its allocation. This is a long-term effect and it is long-term
       !  because ::vegstress_season integrates over ::longevity_sap. For the moment no
       !  short term effects to drought are implemented. Short-term effects should be 
       !  implmented on mortality (through loss of turgor, heat stress, carbon starvation).
       wstress_season(:,1) = zero
       wstress_month(:,1) = zero

       DO jv = 2,nvm

          ! Calculate waterstress
          ! Water stress used in stomate
          ! Set wstress to 1-vegstress so that the value is consistent with its
          ! meaning. wstress=0 indicates no stress, wstress=1 indicates stress
          WHERE (veget_max(:,jv).GT.min_stomate) 
             wstress_season(:,jv) = un - MAX(vegstress_season(:,jv), min_water_stress)
             wstress_month(:,jv) = un - MAX(vegstress_month(:,jv), min_water_stress)
          END WHERE

          !+++CHECK+++
          ! The reduction of the leaf longevity should probably depend on
          ! the leaf skin temperature that will become available through
          ! the multi-layer energy budget. For the moment we don't have
          ! water stress on the leaves. A long term adaption could be 
          ! through ::sla
          longevity_eff_leaf(:,jv) = longevity_leaf(jv) * un

          ! The reduction of the root longevity depends on the soil moisture
          ! stress which we believe is reasonably well captured by our 
          ! proxy for wstress. We need to produce more roots to take the water
          ! from those layers that have water.
          ! feedback to c-allocation has been switched off
          longevity_eff_root(:,jv) = longevity_root(jv) *  un

          ! The reduction of sapwood longevity should depend on the cavitation
          ! which is calculated in hydraulic_arch module. Should be linked
          ! once the memory of cavitation is calculated
          longevity_eff_sap(:,jv) = longevity_sap(jv) *  un
          !+++++++++

       ENDDO

       ! Add to history files
       ! If the soil-based wstress is used, the variables vegstress_xxx 
       ! reflect the moisture in the soil. If the hydraulic architecture
       ! is used vegstress_xxx reflect the ratio between the potential
       ! and actual gpp. wstress_xxx variables are basically 1-vegstress_xxx 
       ! with a minimal value. The xxx_month and xxx_season series of both
       ! variables are almost identical. Hence they were given a very
       ! different output level. 
       CALL histwrite_p (hist_id_stomate, 'VEGSTRESS_WEEK', itime, &
            vegstress_week, kjpindex*nvm, horipft_index)

       CALL xios_orchidee_send_field("VEGSTRESS_DAY",vegstress_day(:,:))
       CALL xios_orchidee_send_field("VEGSTRESS_WEEK",vegstress_week(:,:))
       CALL xios_orchidee_send_field("VEGSTRESS_MONTH",vegstress_month(:,:))

       !! 5.3 Use all processes included in stomate
       !! 5.3.1  Activate stomate processes
       ! Activate stomate processes (the complete list of processes depends
       ! on whether the DGVM is used or not). Processes include: climate constraints
       ! for PFTs, PFT dynamics, Phenology, Allocation, NPP (based on GPP and
       ! autotrophic respiration), fire, mortality, vmax, assimilation temperatures,
       ! all turnover processes, light competition, sapling establishment, lai and
       ! land cover change.

       ! Guillaume M. -- Refresh the fragmentation Average Edge Distance (AED_fire,
       ! AED_light) before any consumer; no-op when all OK_AED_* flags are FALSE.
       ! One edge length only, edge_length_dyn, evolved by update_edge_length at the end of
       ! the previous step: feeding AED_light from a static map while the disturbance
       ! feedback uses the dynamic state would make the two describe different landscapes.
       CALL calculate_aed(kjpindex, veget_max, area, contfrac, edge_length_dyn)

       CALL stomate_lpj_vegetation (kjpindex, dt_days, &
            &             neighbours, resolution, herbivores, &
            &             tsurf_daily, tsoil_daily, t2m_daily, t2m_min_daily, vpd_daily_mean, vpd_daily_max, &
            &             vpd_mean_week, vpd_max_week, &
            &             litterhum_daily, vegstress, humrel, &
            &             maxvegstress_lastyear, maxvegstress_thisyear, &
            &             minvegstress_lastyear, minvegstress_thisyear, &
            &             gdd0_lastyear, &
            &             precip_lastyear, precip_thisyear, &
            &             vegstress_month, vegstress_week, &
            &             t2m_longterm, t2m_month, t2m_week, tau_longterm, &
            &             tsoil_month,precip_month, &
            &             gdd_m5_dormance, gdd_from_growthinit, gdd_midwinter, ncd_dormance, ngd_minus5, &
            &             turnover_longterm, gpp_daily, gpp_week, gpp_year, gpp_decade, resp_maint_week, &
            &             time_hum_min, hum_min_dormance, &
            &             maxfpc_lastyear, maxfpc_thisyear, resp_maint_part,&
            &             PFTpresent, age, fireindex, firelitter, &
            &             leaf_age, leaf_frac, adapted, regenerate, &
            &             plant_status,when_growthinit, litter, &
            &             dead_leaves, som, som_surf, lignin_struc, lignin_wood, lignin_snag, &
            &             veget_max, veget_max_new, veget, fraclut, npp_longterm, croot_longterm, &
            &             lm_lastyearmax, lm_thisyearmax, &
            &             veget_lastlight, everywhere, need_adjacent, RIP_time, &
            &             npp_daily, turnover_daily, turnover_resid, turnover_time,&
            &             control_moist_inst, control_temp_inst, som_input_daily, &
            &             atm_to_immob_daily, co2_fire, emissions_fire, &
            &             resp_hetero_d, resp_hetero_litter_d, resp_hetero_soil_d, resp_maint_d, resp_growth_d, &
            &             deadleaf_cover, assim_param, qsintveg, &
            &             bm_to_litter, bm_to_litter_resid, tree_bm_to_litter, tree_bm_to_litter_resid, &
            &             prod_s, prod_m, prod_l, flux_s, flux_m, flux_l, &
            &             flux_prod_s, flux_prod_m, flux_prod_l, carb_mass_total, &
            &             fpc_max, MatrixA, MatrixV, VectorB, VectorU, &
            &             deepSOM_a, deepSOM_s, deepSOM_p, &
            &             Tseason, Tmin_spring_time, KF, k_latosa_adapt,&
            &             cn_leaf_min_season, nstress_season, vegstress_season, soil_n_min, &
            &             rue_longterm, plant_n_uptake_daily, &
            &             circ_class_n, circ_class_biomass, forest_managed, &
            &             longevity_eff_leaf, longevity_eff_sap, longevity_eff_root, &
            &             species_change_map, fm_change_map, lpft_replant, &
            &             age_stand, age_stand_bm, age_stand_area, rotation_n, last_cut, mai, pai, &
            &             previous_wood_volume, mai_count, coppice_dens, &
            &             harvest_pool_bound, harvest_pool_acc, & 
            &             harvest_type, harvest_cut, harvest_area_acc, &
            &             lai_per_level, laieff_fit, Light_Abs_Tot, Light_Tran_Tot, &
            &             laieff_isotrop, z_array_out, max_height_store, &
            &             wstress_month, wstress_season, &
            &             light_tran_to_floor_season, p_O2, bact, &
            &             CN_som_litter_longterm, &
            &             max_wind_speed_storm, max_wind_ratio_storm, count_storm, is_storm, & ! JJ 2026: threaded companions, removed wind_max_daily & wind_longterm
            &             wind_ratio_max_save, wind_ratio_sum_save, wind_speed_max_save, &
            &             soil_temp_daily, gap_area_save, &
            &             woodharvestpft, &
            &             fDeforestToProduct, fLulccResidue,fHarvestToProduct,  &
            &             cn_leaf_min_2D, cn_leaf_max_2D, cn_leaf_init_2D,bm_sapl_2D, & 
            &             sugar_load, daylight_count, &
            &             n_reserve_longterm, loss_gain, frac_nobio, frac_nobio_new, &
            &             burried_litter, burried_fresh_ltr, &
            &             burried_fresh_som, burried_bact, &
            &             burried_min_nitro, &
            &             burried_som, burried_deepSOM_a, burried_deepSOM_s, &
            &             burried_deepSOM_p,&
            &             i_beetles_generation, &
            &             season_drought_legacy, wood_leftover_legacy, &
            &             i_beetles_activity_legacy, P_beetles_attacked_legacy, &
            &             B_beetles_kill_legacy, &
            &             beetle_diapause, sumTeff, woody_litter_to_use, &
            &             kill_vessels, vessel_mortality_daily, vessel_loss_previous, &
            &             pop_dens, a_nd, edge_length,  &
            &             daily_wspeed_fire, ni_acc, litterfuel, lightning(:,month_end), &
            &             biomass_init_drought, leaf_age_crit, leaf_classes, &
            &             grow_season_len, doy_start_gs, doy_end_gs, mean_start_gs, &
            &             emission_daily, leaching_daily, n_input, &
            &             n_input_daily, fco2_flux, &
            &             nbp_accu_flux, nbp_pool_start, nbp_daily_flux, &
            &             root_profile, total_ba_init, &
            &             maxgppweek_lastyear, maxgppweek_thisyear, us, &
            &             snag_to_wood_flux,   mbc_stomate_lpj, &
            &             tcarbon_acro, tcarbon_cato, carbon_acro, carbon_cato, height_acro, height_cato)

       ! Guillaume M. -- Management intensity: derive rotation, alpha_harvest and
       ! clearcut_area from the class fractions before integrating the edge. veget_max
       ! carries the age structure that modulates plot size. No-op when
       ! ok_management_intensity=FALSE. See design/MODULE_DESIGN_MANAGEMENT_INTENSITY.md.
       CALL set_management_intensity(kjpindex, veget_max)

       ! Guillaume M. -- AED_FEEDBACK: integrate edge_length_dyn from the perturbation
       ! accumulators (dA_fire / dA_storm / dA_pest) filled during stomate_lpj_vegetation.
       ! No-op when ok_aed_feedback=FALSE. See MODULE_DESIGN_AED_FEEDBACK.md.
       CALL update_edge_length(kjpindex, dt_days * one_day, area, contfrac, veget_max)

       !! 5.5.2 Long term adaptation of allocation to water stress
       ! ::wstress season is calculated as the seasonal mean of the
       ! ratio between the stressed and unstressed GPP. If the plant
       ! experiences a short spell of drought, leaves will be killed
       ! (see stomate_turnover). However, when the drought stress is 
       ! maintained during the season, it is assumed that the plant 
       ! will economise its canopy and therefore adjust its allocation
       ! factors to grow less leaves. To avoid that the canopy can only
       ! shrink the plants will try to grow more leaves when they do
       ! not experience any water stress. Extending the LAI, however,
       ! increases the chances that the plant will experience drough 
       ! stress in the future. These feedbacks should stabilize the LAI.
       wstress_adapt(:,:) = zero
       WHERE (wstress_season(:,:) .LT. 0.01)

          ! Increase the leaf allocation by 5% over the whole
          ! year. If there was no water stress during the whole 
          ! year, the following year more C will be allocated
          ! to the leaves.
          wstress_adapt(:,:) = 1.05

       ELSEWHERE

          wstress_adapt(:,:) = wstress_season(:,:)

       ENDWHERE

       DO j = 2,nvm

          WHERE ( k_latosa_adapt(:,j) .GE. k_latosa_max(j) )

             k_latosa_adapt(:,j) = k_latosa_max(j)

          ENDWHERE

       ENDDO

       CALL histwrite_p (hist_id_stomate, 'K_LATOSA_ADAPT', itime, &
            k_latosa_adapt(:,:), kjpindex*nvm, horipft_index)
       CALL xios_orchidee_send_field("K_LATOSA_ADAPT",k_latosa_adapt)


       !! Outputs from Stomate
       
       ! CO2 from wood harvest
       fco2_wh(:) = flux_prod_s(:,icarbon,iharvest,iforest) + &
                    flux_prod_m(:,icarbon,iharvest,iforest) + &
                    flux_prod_l(:,icarbon,iharvest,iforest) + &
                    flux_prod_s(:,icarbon,ilcc,iforest) + &
                    flux_prod_m(:,icarbon,ilcc,iforest) + &
                    flux_prod_l(:,icarbon,ilcc,iforest)
       
       ! CO2 from crop and grassland
       fco2_ha(:) = flux_prod_s(:,icarbon,iharvest,icrop) + &
                    flux_prod_m(:,icarbon,iharvest,icrop) + &
                    flux_prod_l(:,icarbon,iharvest,icrop) + &
                    flux_prod_s(:,icarbon,ilcc,icrop) + &
                    flux_prod_m(:,icarbon,ilcc,icrop) + &
                    flux_prod_l(:,icarbon,ilcc,icrop) + &
                    flux_prod_s(:,icarbon,iharvest,igrass) + &
                    flux_prod_m(:,icarbon,iharvest,igrass) + &
                    flux_prod_l(:,icarbon,iharvest,igrass) + &
                    flux_prod_s(:,icarbon,ilcc,igrass) + &
                    flux_prod_m(:,icarbon,ilcc,igrass) + &
                    flux_prod_l(:,icarbon,ilcc,igrass)
       
       ! Calculate the total CO2 flux from land use change
       fco2_lu(:) = fco2_wh(:) + fco2_ha(:)
       
       !! 5.6b update forcing variables for soil carbon in soil
       IF ( ok_soil_carbon_discretization .AND. ok_forcesoil_write ) THEN

          ! NOTE: This is currently working only for calendrier with 365days 
          ! and not for gregorian calendrier, see ticket 550
          sf_time = MODULO(REAL(days_since_beg,r_std)-1,one_year*REAL(nbyear,r_std))
          iatt=FLOOR(sf_time/dt_forcesoil)+1
          IF ((iatt < 1) .OR. (iatt > nparan*nbyear)) THEN
                WRITE(numout,*) 'Error with days_since_beg=',days_since_beg
                WRITE(numout,*) 'Error with nbyear=',nbyear
                WRITE(numout,*) 'Error with nparan=',nparan
                WRITE(numout,*) 'Error with sf_time=',sf_time
                WRITE(numout,*) 'Error with dt_forcesoil=',dt_forcesoil
                WRITE(numout,*) 'Error with iatt=',iatt
                CALL ipslerr_p (3,'stomate', &
                     &          'Error with iatt.', '', &
                     &          '(Problem with dt_forcesoil ?)')
          ENDIF

          iatt_old=iatt
             
          nforce(iatt) = nforce(iatt) + 1
          som_input_2pfcforcing(:,:,:,:,iatt) = som_input_2pfcforcing(:,:,:,:,iatt) + &
               som_input_daily(:,:,:,:)
          pb_2pfcforcing(:,iatt) = pb_2pfcforcing(:,iatt) + pb_pa_daily(:)
          snow_2pfcforcing(:,iatt) = snow_2pfcforcing(:,iatt) + snow_daily(:)
          tprof_2pfcforcing(:,:,:,iatt) = tprof_2pfcforcing(:,:,:,iatt) + tdeep_daily(:,:,:)
          !cdk treat fbact differently so that we take the mean rate, not the mean
          !residence time
          fbact_2pfcforcing(:,:,:,iatt) = fbact_2pfcforcing(:,:,:,iatt) + decomp_rate_daily(:,:,:)
          hslong_2pfcforcing(:,:,:,iatt) = hslong_2pfcforcing(:,:,:,iatt) + hsdeep_daily(:,:,:)
          veget_max_2pfcforcing(:,:,iatt) = veget_max_2pfcforcing(:,:,iatt) + veget_max(:,:)          ! no need to accum, it is fixed
          DO j=1,nvm
             rprof_2pfcforcing(:,j,iatt) = rprof_2pfcforcing(:,j,iatt) + 1./humcste(j)
          END DO
          tsurf_2pfcforcing(:,iatt) = tsurf_2pfcforcing(:,iatt) + temp_sol_daily(:)
          !adding two snow forcings
          snowdz_2pfcforcing(:,:,iatt) = snowdz_2pfcforcing(:,:,iatt) + snowdz_daily(:,:)
          CN_target_2pfcforcing(:,:,:,iatt) = CN_target_2pfcforcing(:,:,:,iatt) + CN_target(:,:,:)
          n_mineralisation_2pfcforcing(:,:,iatt) = n_mineralisation_2pfcforcing(:,:,iatt) +&
               n_mineralisation_d(:,:)

          !wtp_pt: level of water above surface
          IF (ok_peat_NoDiscretisation) THEN
             wtp_pt(:,iatt)=wtp_pt(:,iatt)+wtp_daily(:)
          ENDIF
          

       ENDIF ! ok_soil_carbon_discretization .AND. ok_forcesoil_write


       ! Cumulate the mass balance closure from stomate and stomate_lpj to be output for analyse
       mbc_daily = mbc_daily + mbc_stomate_lpj
       CALL xios_orchidee_send_field("MBC_daily_c",mbc_daily(:,icarbon))
       CALL xios_orchidee_send_field("MBC_daily_n",mbc_daily(:,initrogen))
       
       !! Reset daily variables
       vegstress_day(:,:) = zero
       litterhum_daily(:) = zero
       t2m_daily(:) = zero
       t2m_min_daily(:) = large_value
       t2m_max_daily(:) = zero
       tsurf_daily(:) = zero
       tsoil_daily(:,:) = zero
       precip_daily(:) = zero
       vpd_daily_mean(:) = zero
       vpd_daily_max(:) = zero
       gpp_daily(:,:) = zero
       resp_maint_part(:,:,:)=zero
       resp_hetero_d=zero
       resp_hetero_litter_d=zero
       resp_hetero_soil_d=zero
       drainage_daily(:,:) = zero 
       plant_n_uptake_daily(:,:,:)=zero 
       atm_to_immob_daily(:,:)=zero
       leaching_daily(:,:,:)=zero
       emission_daily(:,:,:)=zero
       n_input_daily(:,:,:)=zero
       n_mineralisation_d(:,:)=zero 
       tdeep_daily=zero
       hsdeep_daily=zero
       decomp_rate_daily=zero
       snow_daily=zero
       pb_pa_daily=zero
       temp_sol_daily=zero
       snowdz_daily=zero
       snowrho_daily=zero
       daily_wspeed_fire(:) = zero
       mbc_daily(:,:) = zero
       
       IF (printlev_loc >= 3) THEN
          WRITE(numout,*) 'stomate_main: daily processes done'
       ENDIF       

    END IF ! do_slow
    !! Estimate the half hourly (dt_sechiba) CO2 flux. 
    ! Comment written assuming that dt_sechiba is 1800s. For the first
    ! 47 time steps a first estimate for est_co2flux will be made such
    ! that LMDz gets an CO2 flux at each time step. The last time step
    ! of the day (t=48) we have calculated the exact NEE in
    ! stomate_lpj.f90 and will use this variable to calculate the
    ! residual between est_co2flux_d which is the SUM for t=1 to 47 of
    ! est_co2flux and the exact NEE as calculated in stomate_lpj.f90.
    ! This residual will be assigned to est_co2flux (t=48).
    IF (do_slow) THEN

       ! Last time step of the day (t=48). Estimate of the carbon
       ! land_atmosphere fluxes as the residual of the exact nbp and the
       ! estimates of est_co2flux at the previous time steps. The unit of
       ! this residual are [gC m-2 s-1]. The initial units of est_co2flux_d
       ! are [gC m-2 s-1] but is contains the cummulative value over 47
       ! time steps. NEE is in [gC m-2 d-1]. At the start of the day
       ! est_co2flux_d is set to zero at the beginning of this module.
       ! difference = 48 half hourly blocks - 47 half hourly blocks.
       ! The remainder represents the value of the last half hourly
       ! block (=48th) and therefore needs to be divided by dt_sechiba
       ! to get a value in [gC m-2 pixel s-1]
       est_co2flux(:) = SUM(fco2_flux(:,:),2)/one_day - est_co2flux_d(:)
       ! By plotting the ESTIMATED_NEE against the exact NEE calculated and
       ! outputed in stomate_lpj against each other, one can see how good/bad
       ! this estimate is. The smaller the difference between both fluxes
       ! the better. A difference is unavoidable.
       CALL xios_orchidee_send_field("ESTIMATED_NEE",est_co2flux_d*dt_sechiba)

       ! We are about to pass several variables calculated in stomate_lpj 
       ! all the way up to intersurf.f90. Before passing these variables 
       ! and converting their units we need to make sure the carbon 
       ! balance is closed. fco2_flux, fco2_lu, fco2_wh and fco2_ha are in 
       ! gC m-2 dt_stomate-1 and fco2_lu = fco2_wh + fco2_ha
       ! NBP = fco2_nep + fco2_lu
       nbp_intersurf(:) = SUM(fco2_flux(:,:),2) + fco2_lu(:)

       ! Send the outcome of the mass balance check to the history files
       CALL xios_orchidee_send_field("MBC_NBP3_c", &
            nbp_intersurf(:)+nbp_daily_flux(:,icarbon))

       ! Error checking
       IF (err_act.GE.2) THEN
          error_count(:) = zero
          ! This test can only be done at the of the day when the exact nbp has
          ! calculated in stomate_lpj. Hence, its position within a DO_SLOW loop.
          WHERE (ABS(nbp_intersurf(:)+nbp_daily_flux(:,icarbon)).GT.min_stomate)
             error_count(:) = un
          END WHERE

          ! Write error message if needed
          IF (SUM(error_count).GT.zero) THEN
             DO ipts = 1,kjpindex
                IF (error_count(ipts).GT.zero) THEN
                   
                   ! All variables are written in gC m-1 dt_stomate-1
                   WRITE(numout,*) 'fco2_wh, ', fco2_wh(:)
                   WRITE(numout,*) 'fco2_ha, ', fco2_ha(:)
                   WRITE(numout,*) 'fco2_flux, ', SUM(fco2_flux(:,:),2)
                   WRITE(numout,*) 'NBP_intersurf, ', nbp_intersurf(:)
                   WRITE(numout,*) 'NBP_stomate_lpj, ', nbp_daily_flux(:,icarbon)

                   ! The calculation of NBP in stomate_lpj is checked
                   ! make sure to have err_act=3 when you run into this
                   ! error. If not, rerun with err_act=3. Err_act=3 checks 
                   ! for mass balance closure in each subroutine as well 
                   ! as a consistency check between pool-based and flux-based 
                   ! NBP calculations. If all those checks are passed, and this
                   ! error persists, the best place to start debugging is the 
                   ! calculation of nbp_intersurf. Check whether the same fluxes 
                   ! are included, i.e., fire, leaching, bvoc's, ... 
                   ! If a missing flux is indeed the cause of the crash, add 
                   ! these fluxes and remember to transfer them all the way 
                   ! up to intersurf because nbp is recalculted there before it 
                   ! is being send to LMDZ.
                   CALL ipslerr_p(plev,'stomate.f90', &
                        'NBP calculated as it will be done in intersurf', &
                        'differs from the NBP calculated in stomate_lpj', &
                        'Look in the code for more details on this problem')
                END IF
             END DO
          END IF
       END IF
       !-
       
    ELSE
       ! All time steps (t=1 to 47) excep the last time step of the day.
       ! Estimate of the land_atmosphere fluxes [gC m-2 s-1]. This is an
       ! estimate because resp_maint_part is estimated in stomate but
       ! can be reduced in stomate_lpj. Also fluxes from harvest, land
       ! cover change, product use, etc ... are only calculated at the
       ! end of the year. Units of several of the variables have to be
       ! converted. Initial units: resp_hetero_radia [gC m-2 dt_sechiba-1],
       ! gpp_d [gC m-2 day-1], est_co2flux [gC m-2 pixel s-1], 
       ! and resp_maint_part_radia [gC m-2 dt_sechiba-1]
       ! NOTE that atm_to_bm is not used in stomate.f90. If added be
       ! be careful with its units. The units are [gC m-2 dt_sechiba-1]
       ! for the first 47 time steps but they become [gC m-2 dt_stomate-1]
       ! after the returned from stomate_lpj.f90.
       est_co2flux(:) = SUM( ( ( gpp_d(:,:)/one_day - & 
            SUM(resp_maint_part_radia(:,:,:),3)/dt_sechiba - &
            resp_hetero_radia(:,:)/dt_sechiba ) * &
            veget_max(:,:) ),2)
       est_co2flux_d(:) = est_co2flux_d(:) + est_co2flux(:)
    END IF

    !! Prepare module variables for slowproc
    ! Update some more tricky variables. fco2_flux is calculates
    ! only once per day but it should be send every hal-hour to
    ! sechiba and the orchideedriver. Use the value from the
    ! restart file for the first 47 time steps. When stomate_lpj
    ! is called, the fco2_flux will be recalculated and stored in
    ! the restart. fco2_flux_out [gC m-2 s-1]
    fco2_flux_out(:)=est_co2flux(:)*(1-totfrac_nobio)
    fco2_lu_out(:)=fco2_lu(:)*(1-totfrac_nobio)
    fco2_wh_out(:)=fco2_wh(:)*(1-totfrac_nobio)
    fco2_ha_out(:)=fco2_ha(:)*(1-totfrac_nobio)

       !+++CHECK+++
       ! CHECK whether veget_max is correct when lcc is used. most likely this
       ! would only result in a small error but it is better to avoid this
       ! small error in the first place. Multiply every day with veget_max ?
       !! Respiration and fluxes
       ! In stomate_lpj only part of estimated resp_maint for the different plant
       ! parts may get used. resp_maint can therefore be less than resp_maint_part.
       ! Use the value for resp_maint calculated in stomate_lpj (growth_fun_all.f90)
       resp_maint(:,:) = resp_maint_radia(:,:) * veget_max(:,:)
       resp_maint(:,ibare_sechiba) = zero
       resp_growth(:,:) = resp_growth_d(:,:) * veget_max(:,:) * &
            dt_sechiba / one_day
       resp_growth(:,ibare_sechiba) = zero
       resp_hetero(:,:) = resp_hetero_radia(:,:) * veget_max(:,:)
       temp_growth(:)=t2m_month(:)-tp_00


    ! Write to sechiba history file for debugging purpose
    ! Convert the units from gC/m2/s to kgC/m2/s as for the other
    ! C-fluxes in sechiba_history
    CALL xios_orchidee_send_field("netco2flux",est_co2flux(:)/1.e3)
    ! Convert the units from gC/m2/day to kgC/m2/s as for the other
    ! C-fluxes in sechiba_history.
    CALL xios_orchidee_send_field("fco2_lu",fco2_lu(:)/1.e3/one_day)
    CALL xios_orchidee_send_field("fco2_wh",fco2_wh(:)/1.e3/one_day)
    CALL xios_orchidee_send_field("fco2_ha",fco2_ha(:)/1.e3/one_day)

    ! Count how many years have been passsed since the start of the
    ! simulation. Note that global_years is written to the restart
    ! files so it is cumulative since the start of the spinup.
    IF (ts_annual_proc) THEN

       ! Increase the years counter every ts_annual_proc which is by default 
       ! the last sechiba time step of each year
       global_years = global_years + 1 

    END IF

    !! 7. Analytical spinup
    IF (spinup_analytic) THEN

       tau_CN_longterm = tau_CN_longterm + dt_sechiba/one_day 
       !! 7.1. Update V and U at sechiba time step
       DO m = 2,nvm
          DO j = 1,kjpindex
             ! V <- A * V
             matrixV(j,m,:,:) = MATMUL(matrixA(j,m,:,:),matrixV(j,m,:,:))
             ! U <- A*U + B
             vectorU(j,m,:) = MATMUL(matrixA(j,m,:,:),vectorU(j,m,:)) + vectorB(j,m,:)
          ENDDO ! loop pixels
       ENDDO ! loop PFTS
       
       IF (ts_annual_proc) THEN

          ! 7.2.3 Is global_years is a multiple of the period time ?
          ! 3.2.1 When global_years is a multiple of the spinup_period, we calculate :
          !       1) the mean nbp flux over the period. This value is restarted
          !       2) we solve the matrix system by Gauss Jordan method
          !       3) We test if a point is at equilibrium : if yes, we mark the 
          !          point (ok_equilibrium array)
          !       4) Then we reset the matrix 
          !       5) We erase the carbon_stock calculated by ORCHIDEE by the one 
          !          found by the method
          IF( MOD(global_years, spinup_period) == 0 ) THEN

             WRITE(numout,*) 'Spinup analytic : Calculate if system is in equlibrium. global_years=', &
                  global_years

             ! Tag 2.1 and ORCHIDEE 3.0 calculate an nbp (called nbp_accu) but it seems
             ! that this nbp is not used in any calculations neither is it written to a
             ! history file. Given that this version of ORCHIDEE already has two nbps and
             ! several related variables, it was decided not to add yet another nbp 
             ! variable that appears to be purly diagnostic in the first place.

             carbon_stock(:,ibare_sechiba,:) = zero
             ! Prepare the matrix for the resolution
             ! Add a temporary matrix W which contains I-matrixV
             ! we should take the opposite of matrixV and add the 
             ! identitiy : we solve (I-matrixV)*C = vectorU
             matrixW(:,:,:,:) = moins_un * matrixV(:,:,:,:)
             DO jv = 1,nbpools
                matrixW(:,:,jv,jv) =  matrixW(:,:,jv,jv) + un
             ENDDO
             ! Since the surface pool is not considered when the soil carbon
             ! discretization is activated we force the matrix to be one.
             IF (ok_soil_carbon_discretization) THEN
                matrixW(:,:,isurface_pool,isurface_pool) = un
             ENDIF 
             carbon_stock(:,:,:) = vectorU(:,:,:)

             !  Solve the linear system
             DO m = 2,nvm
                DO j = 1,kjpindex
                   ! the solution will be stored in vectorU : so it should be 
                   ! restarted before loop over kjpindex and nvm, so we solved 
                   ! kjpindex*(nvm-1) (7,7) linear systems
                   CALL gauss_jordan_method(nbpools,matrixW(j,m,:,:),carbon_stock(j,m,:))
                ENDDO ! loop pixels
             ENDDO ! loop PFTS

             ! Reset temporary matrixW
             matrixW(:,:,:,:) = zero 

             previous_stock(:,:,:) = current_stock(:,:,:)
             current_stock(:,:,:) = carbon_stock(:,:,:)
  
             ! The relative error is calculated over the passive carbon pool 
             ! (sum over the pfts) over the pixel.
             CALL error_L1_passive(kjpindex,nvm, nbpools, current_stock, &
                  previous_stock, veget_max, eps_carbon, carbon_eq)   

             !! ok_equilibrium is saved,
             WHERE( carbon_eq(:) .AND. .NOT.(ok_equilibrium(:)) )
                ok_equilibrium(:) = .TRUE.  
             ENDWHERE

             IF (printlev_loc .GT. 4) THEN
                WRITE(numout,*) 'current_stock actif:', &
                     current_stock(test_grid,test_pft,iactive)
                WRITE(numout,*) 'current_stock slow:',&
                     current_stock(test_grid,test_pft,islow)
                WRITE(numout,*) 'current_stock passif:', &
                     current_stock(test_grid,test_pft,ipassive)
                WRITE(numout,*) 'current_stock surface:', &
                     current_stock(test_grid,test_pft,isurface)
             END IF

             ! Reset matrixV for the pixel to the identity matrix and vectorU to zero
             matrixV(:,:,:,:) = zero
             vectorU(:,:,:) = zero
             DO jv = 1,nbpools
                matrixV(:,:,jv,jv) = un
             END DO

             IF (printlev >= 2) WRITE(numout,*) 'Reset for matrixV and VectorU done'    

             !! Write the values found in the standard outputs of ORCHIDEE 
             litter(:,istructural,:,iabove,icarbon) = carbon_stock(:,:,istructural_above)
             litter(:,istructural,:,ibelow,icarbon) = carbon_stock(:,:,istructural_below)
             litter(:,imetabolic,:,iabove,icarbon)  = carbon_stock(:,:,imetabolic_above)
             litter(:,imetabolic,:,ibelow,icarbon)  = carbon_stock(:,:,imetabolic_below)
             litter(:,iwoody,:,iabove,icarbon)      = carbon_stock(:,:,iwoody_above)
             litter(:,iwoody,:,ibelow,icarbon)      = carbon_stock(:,:,iwoody_below)
             litter(:,isnag,:,iabove,icarbon)       = carbon_stock(:,:,isnag_above)
             litter(:,isnag,:,ibelow,icarbon)       = carbon_stock(:,:,isnag_below)
             IF (ok_soil_carbon_discretization) THEN
                DO igrn = 1, ngrnd 
                       deepSOM_a(:,igrn,:,icarbon) = carbon_stock(:,:,iactive_pool)* &
                                                    rel_deepSOM_a(:,igrn,:,icarbon)/ &
                                                    (zf_soil(igrn)-zf_soil(igrn-1))
                       deepSOM_s(:,igrn,:,icarbon) = carbon_stock(:,:,islow_pool)* &
                                                    rel_deepSOM_s(:,igrn,:,icarbon)/ &
                                                    (zf_soil(igrn)-zf_soil(igrn-1))
                       deepSOM_p(:,igrn,:,icarbon) = carbon_stock(:,:,ipassive_pool)* &
                                                    rel_deepSOM_p(:,igrn,:,icarbon)/ &
                                                    (zf_soil(igrn)-zf_soil(igrn-1))
                ENDDO
             ELSE
                som(:,iactive,:,icarbon)               = carbon_stock(:,:,iactive_pool)
                som(:,isurface,:,icarbon)              = carbon_stock(:,:,isurface_pool)
                som(:,islow,:,icarbon)                 = carbon_stock(:,:,islow_pool)
                som(:,ipassive,:,icarbon)              = carbon_stock(:,:,ipassive_pool) 
             ENDIF

             WHERE( CN_som_litter_longterm(:,:,istructural_above) .GT. min_stomate)
                litter(:,istructural,:,iabove,initrogen) = &
                     litter(:,istructural,:,iabove,icarbon) &
                     / CN_som_litter_longterm(:,:,istructural_above)
             ENDWHERE
   
             WHERE( CN_som_litter_longterm(:,:,istructural_below) .GT. min_stomate)
                litter(:,istructural,:,ibelow,initrogen) = &
                     litter(:,istructural,:,ibelow,icarbon) &
                     / CN_som_litter_longterm(:,:,istructural_below)
             ENDWHERE
   
             WHERE( CN_som_litter_longterm(:,:,imetabolic_above) .GT. min_stomate)
                litter(:,imetabolic,:,iabove,initrogen) = &
                     litter(:,imetabolic,:,iabove,icarbon) &
                     / CN_som_litter_longterm(:,:,imetabolic_above)
             ENDWHERE
   
             WHERE( CN_som_litter_longterm(:,:,imetabolic_below) .GT. min_stomate)
                litter(:,imetabolic,:,ibelow,initrogen) = &
                     litter(:,imetabolic,:,ibelow,icarbon)  &
                     / CN_som_litter_longterm(:,:,imetabolic_below)
             ENDWHERE
             
             WHERE( CN_som_litter_longterm(:,:,iwoody_above) .GT. min_stomate)
                litter(:,iwoody,:,iabove,initrogen) = &
                     litter(:,iwoody,:,iabove,icarbon)    &
                     / CN_som_litter_longterm(:,:,iwoody_above)    
             ENDWHERE
   
             WHERE( CN_som_litter_longterm(:,:,iwoody_below) .GT. min_stomate)
                litter(:,iwoody,:,ibelow,initrogen) =  &
                     litter(:,iwoody,:,ibelow,icarbon)     &
                     / CN_som_litter_longterm(:,:,iwoody_below)    
             ENDWHERE
             
             WHERE( CN_som_litter_longterm(:,:,isnag_above) .GT. min_stomate)
                litter(:,isnag,:,iabove,initrogen) = &
                     litter(:,isnag,:,iabove,icarbon)    &
                     / CN_som_litter_longterm(:,:,isnag_above)    
             ENDWHERE
   
             WHERE( CN_som_litter_longterm(:,:,isnag_below) .GT. min_stomate)
                litter(:,isnag,:,ibelow,initrogen) =  &
                     litter(:,isnag,:,ibelow,icarbon)     &
                     / CN_som_litter_longterm(:,:,isnag_below)    
             ENDWHERE

            IF (ok_soil_carbon_discretization) THEN
                DO igrn= 1, ngrnd
                  WHERE(CN_som_litter_longterm(:,:,iactive_pool) .GT. min_stomate)
                         deepSOM_a(:,igrn,:,initrogen) = &
                             deepSOM_a(:,igrn,:,icarbon) &
                             / CN_som_litter_longterm(:,:,iactive_pool)
                  ENDWHERE
                ENDDO

                DO igrn= 1, ngrnd
                  WHERE(CN_som_litter_longterm(:,:,islow_pool) .GT. min_stomate)
                         deepSOM_s(:,igrn,:,initrogen) = &
                             deepSOM_s(:,igrn,:,icarbon) &
                             / CN_som_litter_longterm(:,:,iactive_pool)
                  ENDWHERE
                ENDDO

                DO igrn= 1, ngrnd
                  WHERE(CN_som_litter_longterm(:,:,ipassive_pool) .GT. min_stomate)
                         deepSOM_p(:,igrn,:,initrogen) = &
                             deepSOM_p(:,igrn,:,icarbon) &
                             / CN_som_litter_longterm(:,:,iactive_pool)
                  ENDWHERE
                ENDDO
            ELSE
                WHERE(CN_som_litter_longterm(:,:,iactive_pool) .GT. min_stomate)
                   som(:,iactive,:,initrogen) = &
                        som(:,iactive,:,icarbon)    &
                        / CN_som_litter_longterm(:,:,iactive_pool)    
                ENDWHERE
                
                WHERE(CN_som_litter_longterm(:,:,isurface_pool) .GT. min_stomate)
                    som(:,isurface,:,initrogen) = &
                        som(:,isurface,:,icarbon)   &
                        / CN_som_litter_longterm(:,:,isurface_pool)   
                ENDWHERE
   
                WHERE(CN_som_litter_longterm(:,:,islow_pool) .GT. min_stomate)
                    som(:,islow,:,initrogen) = &
                         som(:,islow,:,icarbon)       &
                         / CN_som_litter_longterm(:,:,islow_pool)      
                ENDWHERE
   
                WHERE(CN_som_litter_longterm(:,:,ipassive_pool) .GT. min_stomate)
                    som(:,ipassive,:,initrogen) = &
                        som(:,ipassive,:,icarbon)     &
                         / CN_som_litter_longterm(:,:,ipassive_pool)     
                ENDWHERE
             ENDIF

             CN_som_litter_longterm(:,:,:) = zero
             tau_CN_longterm = dt_sechiba/one_day
             ! Final step, test if all points at the local domain are at equilibrium
             ! The simulation can be stopped when all local domains have 
             ! reached the equilibrium
             IF (printlev >=1) THEN
                IF (ALL(ok_equilibrium)) THEN
                   WRITE(numout,*) 'Spinup analytic : Equilibrium for carbon &
                        &pools is reached for current local domain'
                ELSE
                   WRITE(numout,*) 'Spinup analytic : Equilibrium for carbon &
                        &pools is not yet reached for current local domain'
                END IF
             END IF

          ENDIF ! ( MOD(global_years,spinup_period) == 0)

       ENDIF ! (ts_annual_proc)
       
    ENDIF !(spinup_analytic)
   
    !! Consistency cross-checking (in stomate_lpj.f90)
    IF (do_slow .AND. spinup_analytic .AND. ts_annual_proc .AND. &
         MOD(global_years, spinup_period) .EQ. 0) THEN
       
       ! During this time step soil carbon was recalculated by 
       ! making use of the analytical spinup. This recalculation
       ! violates mass conservation. Cross-checks will thus fail.
       ! Recalculate nbp_accu_flux and nbp_pool_start. 
       CALL calculate_nbp_pool(kjpindex, veget_max, litter, deepSOM_a, &
            deepSOM_s, deepSOM_p, zf_soil, som, bm_to_litter, &
            turnover_daily, circ_class_biomass, circ_class_n, &
            harvest_pool_acc, prod_s, prod_m, prod_l, soil_n_min, &
            nbp_pool_start)

       ! Make sure that at the next time step the cross-check starts
       ! with the pools as updated in the spinup
       nbp_accu_flux(:,:) = nbp_pool_start(:,:)

    END IF ! do_slow

    ! Error checking
    IF(err_act.GT.1)THEN
   
       ! All initial checks should be done in slowproc right after the map
       ! is being read. If vegetation fractions or frac_nobio is adjusted
       ! afterwards, mass balance problems are unavoidable. Check whether
       ! veget_max and frac_nobio are still consistent.
       
       ! Quality check. It is still expected that the different vegetation 
       ! fractions in each pixel sums up to exactly one.
       CALL check_pixel_area("End of stomate", kjpindex, veget_max, frac_nobio)
       
       ! Note that the other check can only be performed the day of the change
       
    END IF ! err_act.GT.1

    IF (printlev >= 4) WRITE(numout,*) 'Leaving stomate_main'

  END SUBROUTINE stomate_main

!! ================================================================================================================================
!! SUBROUTINE 	: stomate_finalize
!!
!>\BRIEF        Write variables to restart file
!!
!! DESCRIPTION  : Write variables to restart file
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCES	: 
!!
!! \n
!_ ================================================================================================================================

  SUBROUTINE stomate_finalize (kjit, kjpindex, index, clay, silt, bulk, assim_param , &
       heat_Zimov, altmax, depth_organic_soil, circ_class_biomass, circ_class_n, &
       lai_per_level, laieff_fit, veget_max)
     
    IMPLICIT NONE
    
    !! 0. Variable and parameter declaration
    !! 0.1 Input variables
    INTEGER(i_std),INTENT(in)                       :: kjit              !! Time step number (unitless)
    INTEGER(i_std),INTENT(in)                       :: kjpindex          !! Domain size - terrestrial pixels only (unitless)
    INTEGER(i_std),DIMENSION(:),INTENT(in)          :: index             !! Indices of the terrestrial pixels only (unitless)
    REAL(r_std),DIMENSION(:),INTENT(in)             :: clay              !! Clay fraction of soil (0-1, unitless)
    REAL(r_std),DIMENSION(:),INTENT(in)             :: silt              !! Silt fraction of soil (0-1, unitless) 
    REAL(r_std),DIMENSION(:),INTENT(in)             :: bulk              !! Bulk density (kg/m**3)
    REAL(r_std),DIMENSION(:,:,:),INTENT(in)         :: assim_param       !! min+max+opt temperatures (K) & vmax for 
                                                                         !! photosynthesis   
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(in)   :: circ_class_biomass!!  Biomass components of the model tree  
                                                                         !! within a circumference class
                                                                         !! class @tex $(g C ind^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)       :: circ_class_n      !! Number of trees within each circumference

                                                                         !! class @tex $(m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)       :: lai_per_level     !! This is the LAI per vertical level
                                                                         !! @tex $(m^{2} m^{-2})$
    TYPE(laieff_type),DIMENSION (:,:,:),INTENT(in) :: laieff_fit         !! Fitted parameters for the effective LAI
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)      :: heat_Zimov         !! heating associated with decomposition [W/m**3 soil]
    REAL(r_std),DIMENSION(:,:),INTENT(in)          :: altmax             !! Maximul active layer thickness (m). Be careful, here active means non frozen.
                                                                         !! Not related with the active soil carbon pool.
    REAL(r_std), DIMENSION(:), INTENT(in)          :: depth_organic_soil !! Depth at which there is still organic matter (m)
    REAL(r_std),DIMENSION(:,:),INTENT(in)          :: veget_max          !! Maximum fraction of vegetation type including non-biological fraction (unitless)
    
    !! 0.2 Output variables
    
    !! 0.3 Modified variables
    
    !! 0.4 Local variables
    REAL(r_std)                                   :: dt_days_read             !! STOMATE time step read in restart file (days)
    INTEGER(i_std)                                :: l,k,ji, jv, i, j, m      !! indices    
    REAL(r_std),PARAMETER                         :: max_dt_days = 5.         !! Maximum STOMATE time step (days)
    REAL(r_std)                                   :: hist_days                !! Writing frequency for history file (days)
    REAL(r_std),DIMENSION(0:nslm)                 :: z_soil                   !! Variable to store depth of the different soil layers (m)
    REAL(r_std),DIMENSION(kjpindex)               :: cvegtot                  !! Total "vegetation" cover (unitless)
    REAL(r_std),DIMENSION(kjpindex)               :: precip                   !! Total liquid and solid precipitation  
                                                                              !! @tex $(??mm dt_stomate^{-1})$ @endtex 
    REAL(r_std),DIMENSION(kjpindex,nvm)           :: gpp_d                    !! Gross primary productivity per ground area 
                                                                              !! @tex $(??gC m^{-2} dt_stomate^{-1})$ @endtex  
    REAL(r_std),DIMENSION(kjpindex,nvm)           :: gpp_daily_x              !! "Daily" gpp for teststomate  
                                                                              !! @tex $(??gC m^{-2} dt_stomate^{-1})$ @endtex 
    REAL(r_std),DIMENSION(kjpindex,nvm)           :: vcmax                    !! Maximum rate of carboxylation
                                                                              !! @tex $(\mumol m^{-2} s^{-1})$ @endtex
    REAL(r_std),DIMENSION(kjpindex,nlevs)         :: control_moist_inst       !! Moisture control of heterotrophic respiration 
                                                                              !! (0-1, unitless) 
    REAL(r_std),DIMENSION(kjpindex,nlevs)         :: control_temp_inst        !! Temperature control of heterotrophic 
                                                                              !! respiration, above and below (0-1, unitless) 
    INTEGER(i_std)                                :: ier                      !! Check errors in netcdf call (unitless)
    REAL(r_std)                                   :: sf_time                  !! Intermediate variable to calculate current time 
                                                                              !! step 
    REAL(r_std), DIMENSION(kjpindex)              :: vartmp                   !! Temporary variable
    INTEGER(i_std)                                :: direct                   !! ??
    REAL(r_std), DIMENSION(kjpindex,nvm,nbpools)  :: carbon_stock                  !! Array containing the carbon stock for each pool
                                                                                   !! used by ORCHIDEE
    

!_ ================================================================================================================================
    
    !! 1. Write restart file for stomate
    IF (printlev>=3) WRITE (numout,*) 'Write restart file for STOMATE' 
    CALL writerestart &
         (kjpindex, index, dt_days, days_since_beg, veget_max, &
         adapted, regenerate, &
         vegstress_day, gdd_init_date, litterhum_daily, &
         t2m_daily, t2m_min_daily, t2m_max_daily, tsurf_daily, tsoil_daily, &
         precip_daily,vpd_daily_mean,vpd_daily_max, vpd_mean_week, vpd_max_week, &
         gpp_daily, npp_daily, turnover_daily, turnover_resid, &
         vegstress_month, vegstress_week, vegstress_season, &
         t2m_longterm, tau_longterm, t2m_month, t2m_week, &
         tsoil_month,precip_month, fireindex, firelitter, &
         maxvegstress_lastyear, maxvegstress_thisyear, &
         minvegstress_lastyear, minvegstress_thisyear, &
         maxgppweek_lastyear, maxgppweek_thisyear, &
         gdd0_lastyear, gdd0_thisyear, &
         precip_lastyear, precip_thisyear, &
         gdd_m5_dormance, gdd_from_growthinit, gdd_midwinter, ncd_dormance, ngd_minus5, &
         PFTpresent, npp_longterm, croot_longterm, n_reserve_longterm, &
         lm_lastyearmax, lm_thisyearmax, &
         maxfpc_lastyear, maxfpc_thisyear, &
         turnover_longterm, gpp_week, gpp_year, gpp_decade, resp_maint_part, resp_maint_week, &
         leaf_age, leaf_frac, leaf_age_crit, &
         plant_status, when_growthinit, age, &
         resp_hetero_d, resp_maint_d, resp_growth_d, &
         co2_fire, &
         veget_lastlight, everywhere, need_adjacent, RIP_time, &
         time_hum_min, hum_min_dormance, &
         litter, dead_leaves, &
         som, lignin_struc, lignin_wood, lignin_snag, turnover_time,&
         fco2_lu, fco2_wh, fco2_ha, &
         prod_s, prod_m, prod_l, &
         flux_s, flux_m, flux_l, &
         fDeforestToProduct, fLulccResidue,fHarvestToProduct, &
         bm_to_litter, bm_to_litter_resid, tree_bm_to_litter, &
         tree_bm_to_litter_resid, carb_mass_total, &
         Tseason, Tseason_length, Tseason_tmp, &
         Tmin_spring_time, &
         global_years, ok_equilibrium, nbp_accu_flux, &
         nbp_pool_start, &
         matrixV, vectorU, previous_stock, current_stock, &
         assim_param, CN_som_litter_longterm, &
         tau_CN_longterm, KF, k_latosa_adapt, rue_longterm, &
         cn_leaf_min_season, &
         nstress_season, soil_n_min, p_O2, bact, &
         forest_managed, &
         species_change_map, fm_change_map, lpft_replant, lai_per_level, &
         laieff_fit, wstress_season, wstress_month, &
         age_stand, age_stand_bm, age_stand_area, rotation_n, last_cut, mai, pai, &
         previous_wood_volume, mai_count, coppice_dens, &
         light_tran_to_floor_season, daylight_count, gap_area_save, &
         deepSOM_a, deepSOM_s, deepSOM_p, O2_soil, CH4_soil, O2_snow, CH4_snow, &
         heat_Zimov, altmax, altmax_ind, altmax_lastyear, altmax_ind_lastyear, &
         depth_organic_soil, fixed_cryoturbation_depth, &
         sugar_load, harvest_cut, &
         harvest_pool_acc, harvest_area_acc, burried_litter, burried_fresh_ltr, &
         burried_fresh_som, burried_bact, &
         burried_min_nitro, burried_som, &
         burried_deepSOM_a, burried_deepSOM_s, burried_deepSOM_p,&
         wood_leftover_legacy, i_beetles_activity_legacy, season_drought_legacy,&
         P_beetles_attacked_legacy, beetle_diapause, sumTeff, woody_litter_to_use, &
         i_beetles_generation, B_beetles_kill_legacy, &
         max_wind_speed_storm, max_wind_ratio_storm, wind_ratio_max_save, & ! JJ 2026: added companions
         wind_ratio_sum_save, wind_speed_max_save, wind_mean_save, wind_sum, &
         is_storm, count_storm, biomass_init_drought, kill_vessels, &
         ni_acc, litterfuel, &
         vessel_loss_previous,grow_season_len, doy_start_gs, doy_end_gs, &
         mean_start_gs, total_ba_init, & 
         carbon_acro, carbon_cato, height_acro, deepSOM_peat ) 

    !! 3. Collect variables that force the soil processes in stomate

    !! Write the soil carbon forcing file
    IF ( ok_soil_carbon_discretization .AND. ok_forcesoil_write ) THEN
      WRITE(numout,*) &
           'stomate: writing the forcing file for permafrost carbon spinup'
      !
      DO iatt = 1, nparan*nbyear
         IF ( nforce(iatt) > 0 ) THEN
            som_input_2pfcforcing(:,:,:,:,iatt) = &
                 som_input_2pfcforcing(:,:,:,:,iatt)/REAL(nforce(iatt),r_std)
            pb_2pfcforcing(:,iatt) = &
                 pb_2pfcforcing(:,iatt)/REAL(nforce(iatt),r_std)
            snow_2pfcforcing(:,iatt) = &
                 snow_2pfcforcing(:,iatt)/REAL(nforce(iatt),r_std)
            tprof_2pfcforcing(:,:,:,iatt) = &
                 tprof_2pfcforcing(:,:,:,iatt)/REAL(nforce(iatt),r_std)
            fbact_2pfcforcing(:,:,:,iatt) = &
                 1./(fbact_2pfcforcing(:,:,:,iatt)/REAL(nforce(iatt),r_std))
            !!!cdk invert this so we take the mean decomposition rate rather than the mean
            !residence time
            hslong_2pfcforcing(:,:,:,iatt) = &
                 hslong_2pfcforcing(:,:,:,iatt)/REAL(nforce(iatt),r_std)
            veget_max_2pfcforcing(:,:,iatt) = &
                 veget_max_2pfcforcing(:,:,iatt)/REAL(nforce(iatt),r_std)
            rprof_2pfcforcing(:,:,iatt) = &
                 rprof_2pfcforcing(:,:,iatt)/REAL(nforce(iatt),r_std)
            tsurf_2pfcforcing(:,iatt) = &
                 tsurf_2pfcforcing(:,iatt)/REAL(nforce(iatt),r_std)
            ! Adding another two snow forcing
            snowdz_2pfcforcing(:,:,iatt) = &
                 snowdz_2pfcforcing(:,:,iatt)/REAL(nforce(iatt),r_std)
            snowrho_2pfcforcing(:,:,iatt) = &
                 snowrho_2pfcforcing(:,:,iatt)/REAL(nforce(iatt),r_std)
            CN_target_2pfcforcing(:,:,:,iatt) = &
                 CN_target_2pfcforcing(:,:,:,iatt)/REAL(nforce(iatt),r_std)
            n_mineralisation_2pfcforcing(:,:,iatt) = &
                 n_mineralisation_2pfcforcing(:,:,iatt)/REAL(nforce(iatt),r_std)
   
            IF (ok_peat_NoDiscretisation) THEN
                wtp_pt(:,iatt)= wtp_pt(:,iatt)/REAL(nforce(iatt),r_std)
            ENDIF

         ELSE
            WRITE(numout,*) &
                 &         'We have no soil carbon forcing data for this time step:', &
                 &         iatt
            WRITE(numout,*) ' -> we set them to zero'
            !soilcarbon_input(:,:,:,iatt) = zero
            !control_moist(:,:,iatt) = zero
            !control_temp(:,:,iatt) = zero
            som_input_2pfcforcing(:,:,:,:,iatt) = zero
            pb_2pfcforcing(:,iatt) = zero
            snow_2pfcforcing(:,iatt) = zero
            tprof_2pfcforcing(:,:,:,iatt) = zero
            fbact_2pfcforcing(:,:,:,iatt) = zero
            hslong_2pfcforcing(:,:,:,iatt) = zero
            veget_max_2pfcforcing(:,:,iatt) = zero
            rprof_2pfcforcing(:,:,iatt) = zero
            tsurf_2pfcforcing(:,iatt) = zero
            snowdz_2pfcforcing(:,:,iatt) = zero
            snowrho_2pfcforcing(:,:,iatt) = zero
            CN_target_2pfcforcing(:,:,:,iatt) = zero
            n_mineralisation_2pfcforcing(:,:,iatt) = zero
            
            IF (ok_peat_NoDiscretisation) THEN
                wtp_pt(:,iatt)=zero
            ENDIF
         ENDIF
      ENDDO
      
      IF (printlev >=3) WRITE (numout,*) 'Create Cforcing file : ',TRIM(Cforcing_discretization_name)
      CALL stomate_io_soil_carbon_discretization_write( Cforcing_discretization_name,                 &
                nbp_glo,            nbp_mpi_para_begin(mpi_rank),   nbp_mpi_para(mpi_rank),     nparan,         &
                nbyear,             index_g,                                                                    &
                clay,               depth_organic_soil,             lalo,                                       &
                snowdz_2pfcforcing, snowrho_2pfcforcing,            som_input_2pfcforcing,                      &
                tsurf_2pfcforcing,  pb_2pfcforcing,                 snow_2pfcforcing,                           &
                tprof_2pfcforcing,  fbact_2pfcforcing,              veget_max_2pfcforcing,                      &
                rprof_2pfcforcing,  hslong_2pfcforcing,             CN_target_2pfcforcing,                      &
                n_mineralisation_2pfcforcing,  wtp_pt)

   ENDIF ! ok_soil_carbon_discretization .AND. ok_forcesoil_write
 
  END SUBROUTINE stomate_finalize


!! ================================================================================================================================
!! SUBROUTINE 	: stomate_init
!!
!>\BRIEF        The routine is called only at the first simulation. At that 
!! time settings and flags are read and checked for internal consistency and 
!! memory is allocated for the variables in stomate.
!!
!! DESCRIPTION  : The routine reads the 
!! following flags from the run definition file:
!! -ipd (index of grid point for online diagnostics)\n
!! -ok_herbivores (flag to activate herbivores)\n
!! -treat_expansion (flag to activate PFT expansion across a pixel\n
!! -harvest_agri (flag to harvest aboveground biomass from agricultural PFTs)\n
!! \n
!! Check for inconsistent setting between the following flags:
!! -ok_stomate\n
!! -ok_dgvm\n
!! \n
!! Memory is allocated for all the variables of stomate and new indexing tables 
!! are build. New indexing tables are needed because a single pixel can conatin 
!! several PFTs. The new indexing tables have separate indices for the different 
!! PFTs. Similar index tables are build for land use cover change.\n
!! \n
!! Several global variables and land cover change variables are initialized to 
!! zero.\n
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): Strictly speaking the subroutine has no output 
!! variables. However, the routine allocates memory and builds new indexing 
!! variables for later use.\n 
!!
!! REFERENCE(S)	: None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE stomate_init &
       &  (kjpij, kjpindex, index, lalo, &
       &   rest_id_stom, hist_id_stom, hist_id_stom_IPCC)

  !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std),INTENT(in)                    :: kjpij             !! Total size of the un-compressed grid, including 
                                                                      !! oceans (unitless) 
    INTEGER(i_std),INTENT(in)                    :: kjpindex          !! Domain size - number of terrestrial pixels 
                                                                      !! (unitless) 
    INTEGER(i_std),INTENT(in)                    :: rest_id_stom      !! STOMATE's _Restart_ file identifier
    INTEGER(i_std),INTENT(in)                    :: hist_id_stom      !! STOMATE's _history_ file identifier
    INTEGER(i_std),INTENT(in)                    :: hist_id_stom_IPCC !! STOMATE's IPCC _history_ file identifier 
    INTEGER(i_std),DIMENSION(:),INTENT(in)       :: index             !! Indices of the terrestrial pixels on the global 
                                                                      !! map 
    REAL(r_std),DIMENSION(:,:),INTENT(in)        :: lalo              !! Geogr. coordinates (latitude,longitude) (degrees)
   
    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables

    LOGICAL                                      :: l_error           !! Check errors in netcdf call
    INTEGER(i_std)                               :: ier               !! Check errors in netcdf call
    INTEGER(i_std)                               :: ji,j,ipd,l        !! Indices
    INTEGER(i_std)                               :: idia              !! indices
!_ ================================================================================================================================
    
  !! 1. Online diagnostics

    IF ( kjpindex > 0 ) THEN
       !Config Key  = STOMATE_DIAGPT
       !Config Desc = Index of grid point for online diagnostics
       !Config If    = OK_STOMATE
       !Config Def  = 1
       !Config Help = This is the index of the grid point which
       !               will be used for online diagnostics.
       !Config Units = [-]
       ! By default ::ipd is set to 1
       ipd = 1
       ! Get ::ipd from run definition file
       CALL getin_p('STOMATE_DIAGPT',ipd)
       ipd = MIN( ipd, kjpindex )
       IF ( printlev >=3 ) THEN
          WRITE(numout,*) 'Stomate: '
          WRITE(numout,*) '  Index of grid point for online diagnostics: ',ipd
          WRITE(numout,*) '  Lon, lat:',lalo(ipd,2),lalo(ipd,1)
          WRITE(numout,*) '  Index of this point on GCM grid: ',index(ipd)
       END IF
    ENDIF

  !! 2. Check consistency of flags

    IF ( ( .NOT. ok_stomate ) .AND. ok_dgvm ) THEN
       WRITE(numout,*) 'Cannot do dynamical vegetation without STOMATE.'
       WRITE(numout,*) 'Inconsistency between ::ok_stomate and ::ok_dgvm'
       WRITE(numout,*) 'Stop: fatal error'
       STOP
    ENDIF

    IF (printlev >=2) THEN
       WRITE(numout,*) 'stomate first call - overview of the activated flags:'
       WRITE(numout,*) '  STOMATE: ', ok_stomate
       WRITE(numout,*) '  LPJ: ', ok_dgvm
    END IF

  !! 4. Allocate memory for STOMATE's variables

    l_error = .FALSE.

    ALLOCATE(adapted(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for adapted. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(regenerate(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for regenerate. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(vegstress_day(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for vegstress_day. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(stressed_daily(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for stressed_daily. We stop. We need kjpindex*nvm words',kjpindex,nvm
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF
    stressed_daily = zero

    ALLOCATE(unstressed_daily(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for unstressed_daily. We stop. We need kjpindex*nvm words',kjpindex,nvm
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF
    unstressed_daily(:,:) = zero

    ALLOCATE(biomass_init_drought(kjpindex,nvm,ncirc,nparts,nelements),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for biomass_init_drought. We stop. We need kjpindex*nvm*ncirc*nparts*nelements words',kjpindex,nvm,ncirc,nparts,nelements
       CALL ipslerr_p (3, 'stomate_init', 'Memory allocation issue','','')
    ENDIF

    ALLOCATE(kill_vessels(kjpindex,nvm,ncirc),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for kill_vessels. We stop. We need kjpindex*nvm words.',kjpindex,nvm,ncirc
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(vessel_loss_previous(kjpindex,nvm,ncirc),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for vessel_loss_previous. We stop. We need kjpindex*nvm words',kjpindex,nvm,ncirc
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF

    ALLOCATE(vessel_loss_daily(kjpindex,nvm,ncirc),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for vessel_loss_daily. We stop. We need kjpindex*nvm words',kjpindex,nvm,ncirc
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF
    vessel_loss_daily(:,:,:) = zero

    ALLOCATE(daylight(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for daylight. We stop. We need kjpindex*nvm words',kjpindex,nvm
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF
    daylight(:,:) = zero

    ALLOCATE(light_tran_to_floor_season(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for light_tran_to_floor_season. We stop. We need kjpindex*nvm*nlevels_tot words', &
            kjpindex,nvm,nlevels_tot
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF
    light_tran_to_floor_season(:,:) = zero

    ALLOCATE(daylight_count(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for daylight_count. We stop. We need kjpindex*nvm words',kjpindex,nvm
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF
    daylight_count(:,:) = zero

    ALLOCATE(transpir_daily(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for transpir_daily. We stop. We need kjpindex*nvm words',kjpindex,nvm
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF
    transpir_daily=zero


    ALLOCATE(litterhum_daily(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for litterhum_daily. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(t2m_daily(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for t2m_daily. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(t2m_min_daily(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for t2m_min_daily. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(t2m_max_daily(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       CALL ipslerr_p(3,'sechiba_init','Pb in alloc for t2m_max_daily','','')
    ENDIF

    ALLOCATE(wind_max_daily(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for wind_max_daily','','')
    wind_max_daily(:) = zero

    ALLOCATE(wind_ratio_sum(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for wind_ratio_sum','','')
    wind_ratio_sum(:)=zero

    ALLOCATE(wind_sum(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for wind_sum','','')
    wind_sum(:) = zero

    ALLOCATE(max_wind_speed_storm(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for max_wind_speed_storm','','')
    max_wind_speed_storm(:) = zero

    ! Guillaume M. -- Companion of max_wind_speed_storm.
    ALLOCATE(max_wind_ratio_storm(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for max_wind_ratio_storm','','')
    max_wind_ratio_storm(:) = un

    ALLOCATE(wind_ratio_sum_save(kjpindex,wind_days),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for wind_ratio_sum_save','','')
    wind_ratio_sum_save(:,:) = zero

    ! Guillaume M. -- Companion of wind_ratio_sum_save.
    ALLOCATE(wind_ratio_max_save(kjpindex,wind_days),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for wind_ratio_max_save','','')
    wind_ratio_max_save(:,:) = un

    ALLOCATE(wind_speed_max_save(kjpindex,wind_days),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for wind_speed_max_save','','')
    wind_speed_max_save(:,:) = zero

    ALLOCATE(wind_mean_save(kjpindex,legacy_years_wind),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for wind_mean_save','','')
    wind_mean_save(:,:) = un

    ALLOCATE(is_storm(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for is_storm','','')
    is_storm(:) = .FALSE.

    ALLOCATE(count_storm(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for count_storm','','')
    count_storm(:) = NINT(zero)

    ALLOCATE(wind_max(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for wind_max','','')
    wind_max(:)=zero

    ! Guillaume M. -- Transient daily accumulator, companion of wind_max.
    ALLOCATE(wind_ratio_max(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for wind_ratio_max','','')
    wind_ratio_max(:)=zero

    ALLOCATE(soil_temp_daily(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for soil_temp_daily. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(soil_max_daily(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for soil_max_daily. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF
    soil_max_daily(:) = zero

    ALLOCATE(tsurf_daily(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for tsurf_daily. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(tsoil_daily(kjpindex,nslm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for tsoil_daily. We stop. We need kjpindex*nslm words',kjpindex,nslm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(precip_daily(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for precip_daily. We stop. We need kjpindex words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(vpd_daily_mean(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for vpd_daily_mean. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(vpd_daily_max(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for vpd_daily_max. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(vpd_mean_week(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for vpd_mean_week. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(vpd_max_week(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for vpd_max_week. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF


    ALLOCATE(gpp_daily(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for gpp_daily. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(npp_daily(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for npp_daily. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(turnover_daily(kjpindex,nvm,nparts,nelements),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for turnover_daily. We stop. We need kjpindex*nvm*nparts*nelements words', &
       &   kjpindex,nvm,nparts,nelements
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(turnover_resid(kjpindex,nvm,nparts,nelements),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for turnover_resid. We stop. We need kjpindex*nvm*nparts*nelements words', &
       &   kjpindex,nvm,nparts,nelements
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(turnover_littercalc(kjpindex,nvm,nparts,nelements),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for turnover_littercalc. We stop. We need kjpindex*nvm*nparts*nelements words', & 
        &  kjpindex,nvm,nparts,nelements
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(vegstress_month(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for vegstress_month. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(vegstress_week(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for vegstress_week. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(vegstress_season(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for vegstress_season. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(t2m_longterm(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for t2m_longterm. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(t2m_month(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for t2m_month. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(Tseason(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for Tseason. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(Tseason_length(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for Tseason_length. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(Tseason_tmp(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for Tseason_tmp. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(Tmin_spring_time(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for Tmin_spring_time. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(t2m_week(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for t2m_week. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(tsoil_month(kjpindex,nslm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for tsoil_month. We stop. We need kjpindex*nslm words',kjpindex,nslm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(precip_month(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for precip_month. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(fireindex(kjpindex,nvm),stat=ier) 
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for fireindex. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(firelitter(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for firelitter. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(maxvegstress_lastyear(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for maxvegstress_lastyear. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(maxvegstress_thisyear(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for maxvegstress_thisyear. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(minvegstress_lastyear(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for minvegstress_lastyear. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(minvegstress_thisyear(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for minvegstress_thisyear. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(maxgppweek_lastyear(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for maxgppweek_lastyear. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(maxgppweek_thisyear(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for maxgppweek_thisyear. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(gdd0_lastyear(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for gdd0_lastyear. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(gdd0_thisyear(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for gdd0_thisyear. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(gdd_init_date(kjpindex,2),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for gdd_init_date. We stop. We need kjpindex*2 words',kjpindex,2
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(gdd_from_growthinit(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for gdd_from_growthinit. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(precip_lastyear(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for precip_lastyear. We stop. We need kjpindex*nvm words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(precip_thisyear(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for precip_thisyear. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(gdd_m5_dormance(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for gdd_m5_dormance. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(gdd_midwinter(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for gdd_midwinter. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(ncd_dormance(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for ncd_dormance. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(ngd_minus5(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for ngd_minus5. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(PFTpresent(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for PFTpresent. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(npp_longterm(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for npp_longterm. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(croot_longterm(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for croot_longterm. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

   ALLOCATE(n_reserve_longterm(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for n_reserve_longterm. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(lm_lastyearmax(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for lm_lastyearmax. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(lm_thisyearmax(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for lm_thisyearmax. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(maxfpc_lastyear(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for maxfpc_lastyear. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(maxfpc_thisyear(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for maxfpc_thisyear. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(turnover_longterm(kjpindex,nvm,nparts,nelements),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for turnover_longterm. We stop. We need kjpindex*nvm*nparts*nelements words', & 
       &    kjpindex,nvm,nparts,nelements
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(gpp_week(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for gpp_week. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(gpp_year(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for gpp_year. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(gpp_decade(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for gpp_decade. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(resp_maint_week(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for resp_maint_week. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(plant_status(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for plant_status. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(when_growthinit(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for when_growthinit. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(age(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for age. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(resp_hetero_d(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for resp_hetero_d. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(resp_hetero_litter_d(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for resp_hetero_litter_d. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(resp_hetero_soil_d(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for resp_hetero_soil_d. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(resp_hetero_radia(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for resp_hetero_radia. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(resp_maint_d(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for resp_maint_d. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(resp_growth_d(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for resp_growth_d. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(co2_fire(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for co2_fire. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(emissions_fire(kjpindex,nvm,nelements,nfiretype),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for emissions_fire. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(atm_to_immob(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for atm_to_immob. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(veget_lastlight(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for veget_lastlight. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(everywhere(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for everywhere. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(need_adjacent(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for need_adjacent. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(leaf_age(kjpindex,nvm,nleafages),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for leaf_age. We stop. We need kjpindex*nvm*nleafages words', & 
       &      kjpindex,nvm,nleafages
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(leaf_frac(kjpindex,nvm,nleafages),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for leaf_frac. We stop. We need kjpindex*nvm*nleafages words', & 
       &      kjpindex,nvm,nleafages
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(RIP_time(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for RIP_time. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(time_hum_min(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for time_hum_min. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(hum_min_dormance(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for hum_min_dormance. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF


    ALLOCATE(litter(kjpindex,nlitt,nvm,nlevs,nelements),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for litter. We stop. We need kjpindex*nlitt*nvm*nlevs*nelements words', & 
       &    kjpindex,nlitt,nvm,nlevs,nelements
       STOP 'stomate_init'
    ENDIF
    
    ALLOCATE(litterfuel(kjpindex,nlitt,nvm,nhour,nelements),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for litterfuel. We stop. We need kjpindex*nlitt*nvm*nhour*nelements words', &
       &    kjpindex,nlitt,nvm,nhour,nelements
       STOP 'stomate_init'
    ENDIF
    
    ALLOCATE(dead_leaves(kjpindex,nvm,nlitt),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for dead_leaves. We stop. We need kjpindex*nvm*nlitt words', & 
       &   kjpindex,nvm,nlitt
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(som(kjpindex,ncarb,nvm,nelements),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for som. We stop. We need kjpindex*ncarb*nvm*nelements words',&
            kjpindex,ncarb,nvm,nelements
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(burried_litter(kjpindex,nlitt,nlevs,nelements),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for burried_litter. We stop. We need kjpindex*nlitt*nlevs*nelements words',&
            kjpindex,nlitt,nlevs,nelements
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(burried_fresh_ltr(kjpindex,nparts,nelements),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for burried_fresh_ltr. We stop. We need kjpindex*nparts*nelements words',&
            kjpindex,nparts,nelements
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(burried_fresh_som(kjpindex,nparts,nelements),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for burried_fresh_som. We stop. We need kjpindex*nparts*nelements words',&
            kjpindex,nparts,nelements
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(burried_bact(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for burried_bact. We stop. We need kjpindex words',&
            kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(burried_min_nitro(kjpindex,nnspec),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for burried_min_nitro. We stop. We need kjpindex*nnspec words',&
            kjpindex,nnspec
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(burried_deepSOM_a(kjpindex,ngrnd,nelements),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for burried_deepSOM_a. We stop. We need kjpindex*ngrnd*ncarb*nelements words',&
            kjpindex,ngrnd,ncarb,nelements
       STOP 'stomate_init'
    ENDIF
    
    ALLOCATE(burried_deepSOM_s(kjpindex,ngrnd,nelements),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for burried_deepSOM_s. We stop. We need kjpindex*ngrnd*ncarb*nelements words',&
            kjpindex,ngrnd,ncarb,nelements
       STOP 'stomate_init'
    ENDIF
    
    ALLOCATE(burried_deepSOM_p(kjpindex,ngrnd,nelements),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for burried_deepSOM_p. We stop. We need kjpindex*ngrnd*ncarb*nelements words',&
            kjpindex,ngrnd,ncarb,nelements
       STOP 'stomate_init'
    ENDIF
    
    ALLOCATE(burried_som(kjpindex,ncarb,nelements),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for burried_som. We stop. We need kjpindex*ncarb*nelements words',&
            kjpindex,ncarb,nelements
       STOP 'stomate_init'
    ENDIF
       
    ALLOCATE(som_surf(kjpindex,ncarb,nvm,nelements),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'stomate_init', 'Pb in alloc for som_surf','','')

    ALLOCATE(lignin_struc(kjpindex,nvm,nlevs),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for lignin_struc. We stop. We need kjpindex*nvm*nlevs words',&
            kjpindex,nvm,nlevs
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(lignin_wood(kjpindex,nvm,nlevs),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for lignin_wood. We stop. We need kjpindex*nvm*nlevs words',&
            kjpindex,nvm,nlevs
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(lignin_snag(kjpindex,nvm,nlevs),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for lignin_snag. We stop. We need kjpindex*nvm*nlevs words',&
            kjpindex,nvm,nlevs
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(turnover_time(kjpindex,nvm,nparts),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for turnover_time. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(bm_to_litter(kjpindex,nvm,nparts,nelements),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for bm_to_litter. We stop. We need kjpindex*nvm*nparts*nelements words', & 
       &    kjpindex,nvm,nparts,nelements
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(bm_to_litter_resid(kjpindex,nvm,nparts,nelements),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for bm_to_litter_resid. We stop. We need kjpindex*nvm*nparts*nelements words', & 
       &    kjpindex,nvm,nparts,nelements
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(tree_bm_to_litter(kjpindex,nvm,nparts,nelements),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for tree_bm_to_litter. We stop. We need kjpindex*nvm*nparts*nelements words', & 
       &    kjpindex,nvm,nparts,nelements
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(tree_bm_to_litter_resid(kjpindex,nvm,nparts,nelements),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for tree_bm_to_litter_resid. We stop. We need kjpindex*nvm*nparts*nelements words', & 
       &    kjpindex,nvm,nparts,nelements
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(bm_to_littercalc(kjpindex,nvm,nparts,nelements),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for bm_to_littercalc. We stop. We need kjpindex*nvm*nparts*nelements words', &
       &   kjpindex,nvm,nparts,nelements
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(tree_bm_to_littercalc(kjpindex,nvm,nparts,nelements),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for tree_bm_to_littercalc. We stop. We need kjpindex*nvm*nparts*nelements words', &
       &   kjpindex,nvm,nparts,nelements
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(herbivores(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for herbivores. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(resp_maint_part_radia(kjpindex,nvm,nparts),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for resp_maint_part_radia. We stop. We need kjpindex*nvm*nparts words', &
       &  kjpindex,nvm,nparts
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(resp_maint_radia(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for resp_maint_radia. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(resp_maint_part(kjpindex,nvm,nparts),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for resp_maint_part. We stop. We need kjpindex*nvm*nparts words', &
       &    kjpindex,nvm,nparts
       STOP 'stomate_init'
    ENDIF
    resp_maint_part(:,:,:) = zero

    ALLOCATE(hori_index(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for hori_index. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(horipft_index(kjpindex*nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for horipft_index. We stop. We need kjpindex*nvm words',kjpindex*nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(horican_index(kjpindex*nlevels_tot),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for horican_index. We stop. We need kjpindex*nlevels_tot words',&
            kjpindex*nlevels_tot
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF

    ALLOCATE(horicut_index(kjpindex*ncut_times),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for horicut_index. We stop. We need kjpindex*ncut_times words',&
            kjpindex*ncut_times
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF

    ALLOCATE (horip_s_index(kjpindex*(nshort+1)), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for horip_s_index. We stop. We need kjpindex*10 words',kjpindex,nshort
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF

    ALLOCATE (horip_m_index(kjpindex*(nmedium+1)), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for horip_m_index. We stop. We need kjpindex*10 words',kjpindex,nmedium
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF

    ALLOCATE (horip_l_index(kjpindex*(nlong+1)), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for horip_l_index. We stop. We need kjpindex*100 words',kjpindex,nlong
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF
 
    ALLOCATE (horip_ss_index(kjpindex*(nshort+1)), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for horip_ss_index. We stop. We need kjpindex*11 words',kjpindex,nshort+1
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF 

    ALLOCATE (horip_mm_index(kjpindex*(nmedium+1)), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for horip_mm_index. We stop. We need kjpindex*11 words',kjpindex,nmedium+1
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF

    ALLOCATE (horip_ll_index(kjpindex*(nlong+1)), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for horip_ll_index. We stop. We need kjpindex*101 words',kjpindex,nlong+1
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF

    ALLOCATE (prod_s(kjpindex,1:nshort+1,nelements,nlanduse,nlctypes), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for prod_s. We stop. We need kjpindex*(nshort+1)*nelements*nlanduse  words', &
            kjpindex,nshort+1,nelements,nlanduse,nlctypes
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF  
    prod_s(:,:,:,:,:) = zero
    

    ALLOCATE (prod_m(kjpindex,1:nmedium+1,nelements,nlanduse,nlctypes), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for prod_m. We stop. We need kjpindex*(nmedium+1)*nelements*nlanduse words', &
            kjpindex,nmedium+1,nelements,nlanduse,nlctypes
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF 
    prod_m(:,:,:,:,:) = zero

    ALLOCATE (prod_l(kjpindex,1:nlong+1,nelements,nlanduse,nlctypes), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for prod_l. We stop. We need kjpindex*(nlong+1)*nelements*nlanduse* words', &
            kjpindex,nlong+1,nelements,nlanduse,nlctypes
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF
    prod_l(:,:,:,:,:) = zero

    ALLOCATE (flux_s(kjpindex,1:nshort+1,nelements,nlanduse,nlctypes), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for flux_s. We stop. We need kjpindex*nshort*nelements*nlanduse words', &
            kjpindex,nshort,nelements,nlanduse,nlctypes
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF
    flux_s(:,:,:,:,:) = zero

    ALLOCATE (flux_m(kjpindex,1:nmedium+1,nelements,nlanduse,nlctypes), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for flux_m. We stop. We need kjpindex*nmedium*nlanduse words', &
            kjpindex,nmedium,nelements,nlanduse,nlctypes
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF
    flux_m(:,:,:,:,:) = zero

    ALLOCATE (flux_l(kjpindex,1:nlong+1,nelements,nlanduse,nlctypes), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for flux_l. We stop. We need kjpindex*nlong*nelements*nlanduse words', &
            kjpindex,nlong,nelements,nlanduse,nlctypes
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF
    flux_l(:,:,:,:,:) = zero

    ALLOCATE (flux_prod_s(kjpindex,nelements,nlanduse,nlctypes), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for flux_prod_s. We stop. We need kjpindex*nelements*nlanduse words', &
            kjpindex,nelements,nlanduse,nlctypes
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF
    flux_prod_s(:,:,:,:) = zero

    ALLOCATE (flux_prod_m(kjpindex,nelements,nlanduse,nlctypes), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for flux_prod_m. We stop. We need kjpindex*nelements*nlanduse words', &
            kjpindex,nelements,nlanduse,nlctypes
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF
    flux_prod_m(:,:,:,:) = zero

    ALLOCATE (flux_prod_l(kjpindex,nelements,nlanduse,nlctypes), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for flux_prod_l. We stop. We need kjpindex*nelements*nlanduse words', &
            kjpindex,nelements,nlanduse,nlctypes
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF
    flux_prod_l(:,:,:,:) = zero

    ALLOCATE (fco2_lu(kjpindex), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for fco2_lu. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (fco2_wh(kjpindex), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for fco2_wh. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (fco2_ha(kjpindex), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for fco2_ha. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (woodharvestpft(kjpindex,nvm), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for woodharvestpft. We stop. We need kjpindex*nvm words',kjpindex*nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (fDeforestToProduct(kjpindex,nvm), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for fDeforestToProduct. We stop. We need kjpindex*nvm words',kjpindex*nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (fLulccResidue(kjpindex,nvm), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for fLulccResidue. We stop. We need kjpindex*nvm words',kjpindex*nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (fHarvestToProduct(kjpindex,nvm), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for fHarvestToProduct. We stop. We need kjpindex*nvm words',kjpindex*nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (carb_mass_total(kjpindex), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for carb_mass_total. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (som_input_daily(kjpindex,ncarb,nvm,nelements), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for som_input_daily. We stop. We need kjpindex*ncarb*nvm*nelements words', & 
       &    kjpindex,ncarb,nvm,nelements
       STOP 'stomate_init'
    ENDIF

    ALLOCATE (fpc_max(kjpindex,nvm), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for fpc_max. We stop. We need kjpindex*nvm words',kjpindex,nvm
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(cn_leaf_min_season(kjpindex,nvm),stat=ier) 
    l_error = l_error .OR. (ier /= 0) 
    IF (l_error) THEN 
       WRITE(numout,*) 'Memory allocation error for cn_leaf_min_season. We stop. We need kjpindex*nvm words',kjpindex,nvm 
       STOP 'stomate_init' 
    ENDIF 
 
    ALLOCATE(nstress_season(kjpindex,nvm),stat=ier) 
    l_error = l_error .OR. (ier /= 0) 
    IF (l_error) THEN 
       WRITE(numout,*) 'Memory allocation error for nstress_season. We stop. We need kjpindex*nvm words',kjpindex,nvm 
       STOP 'stomate_init' 
    ENDIF 
 
    ALLOCATE(soil_n_min(kjpindex,nvm,nnspec),stat=ier) 
    l_error = l_error .OR. (ier /= 0) 
    IF (l_error) THEN 
       WRITE(numout,*) 'Memory allocation error for soil_n_min. We stop. We need kjpindex*nvm words',kjpindex,nvm,nnspec 
       STOP 'stomate_init' 
    ENDIF 

    ALLOCATE(p_O2(kjpindex,nvm),stat=ier) 
    l_error = l_error .OR. (ier /= 0) 
    IF (l_error) THEN 
       WRITE(numout,*) 'Memory allocation error for p_O2. We stop. We need kjpindex*nvm words',kjpindex,nvm 
       STOP 'stomate_init' 
    ENDIF 

    ALLOCATE(bact(kjpindex,nvm),stat=ier) 
    l_error = l_error .OR. (ier /= 0) 
    IF (l_error) THEN 
       WRITE(numout,*) 'Memory allocation error for bact. We stop. We need kjpindex*nvm words',kjpindex,nvm 
       STOP 'stomate_init' 
    ENDIF 

    ALLOCATE(ok_equilibrium(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0) 
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for ok_equilibrium. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(drainage_daily(kjpindex,nvm),stat=ier) 
    l_error = l_error .OR. (ier /= 0) 
    IF (l_error) THEN 
       WRITE(numout,*) ' Memory allocation error for drainage_daily. We stop. We need kjpindex*nvm words = ',kjpindex, nvm
       STOP 'drainage_daily' 
    ENDIF 
 
    ALLOCATE (plant_n_uptake_daily(kjpindex,nvm,nionspec), stat=ier) 
    l_error = l_error .OR. (ier.NE.0) 
    IF (l_error) THEN 
       WRITE(numout,*) ' Memory allocation error for plant_n_uptake_daily. We stop. We need kjpindex words = ',kjpindex*nvm*nionspec 
       STOP 'plant_n_uptake_daily' 
    ENDIF 
 
    ALLOCATE (n_mineralisation_d(kjpindex,nvm), stat=ier) 
    l_error = l_error .OR. (ier.NE.0) 
    IF (l_error) THEN 
       WRITE(numout,*) ' Memory allocation error for n_mineralisation_d. We stop. We need kjpindex words = ',kjpindex*nvm 
       STOP 'n_mineralisation_d' 
    ENDIF 

    ALLOCATE (atm_to_immob_daily(kjpindex,nvm), stat=ier)
    l_error = l_error .OR. (ier.NE.0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for atm_to_immob_daily. We stop. We need kjpindex words = ',kjpindex*nvm
       STOP 'atm_to_immob_daily'
    ENDIF

    ALLOCATE (leaching_daily(kjpindex,nvm,nionspec), stat=ier)
    l_error = l_error .OR. (ier.NE.0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for leaching_daily. We stop. We need kjpindex words = ',kjpindex*nvm*nionspec
       STOP 'leaching_daily'
    ENDIF

    ALLOCATE (emission_daily(kjpindex,nvm,nnspec), stat=ier)
    l_error = l_error .OR. (ier.NE.0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for emission_daily. We stop. We need kjpindex words = ',kjpindex*nvm*nnspec
       STOP 'emission_daily'
    ENDIF

    ALLOCATE (n_input_daily(kjpindex,nvm,ninput), stat=ier)
    l_error = l_error .OR. (ier.NE.0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for n_input_daily. We stop. We need kjpindex words = ',kjpindex*nvm*ninput
       STOP 'n_input_daily'
    ENDIF

    ALLOCATE(carbon_eq(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for carbon_eq. We stop. We need kjpindex words',kjpindex
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(nbp_accu_flux(kjpindex,nelements),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for nbp_accu_flux. We stop. We need kjpindex*nelements words',kjpindex*nelements
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(nbp_pool_start(kjpindex,nelements),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for nbp_pool_start. We stop. We need kjpindex*nelements words',kjpindex*nelements
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(matrixA(kjpindex,nvm,nbpools,nbpools),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for matrixA. We stop. We need kjpindex*nvm*nbpools*nbpools words',  & 
       &     kjpindex, nvm, nbpools, nbpools
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(vectorB(kjpindex,nvm,nbpools),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for vectorB. We stop. We need kjpindex*nvm*nbpools words',  & 
       &     kjpindex, nvm, nbpools
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(vectorU(kjpindex,nvm,nbpools),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for vectorU. We stop. We need kjpindex*nvm*nbpools words',  & 
       &     kjpindex, nvm, nbpools
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(matrixV(kjpindex,nvm,nbpools,nbpools),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for matrixV. We stop. We need kjpindex*nvm*nbpools*nbpools words',  & 
       &     kjpindex, nvm, nbpools, nbpools
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(matrixW(kjpindex,nvm,nbpools,nbpools),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for matrixW. We stop. We need kjpindex*nvm*nbpools*nbpools words',  & 
       &     kjpindex, nvm, nbpools, nbpools
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(previous_stock(kjpindex,nvm,nbpools),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for previous_stock. We stop. We need kjpindex*nvm*nbpools words',  & 
       &     kjpindex, nvm, nbpools
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(current_stock(kjpindex,nvm,nbpools),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for current_stock. We stop. We need kjpindex*nvm*nbpools words',  & 
       &     kjpindex, nvm, nbpools
       STOP 'stomate_init'
    ENDIF

    ALLOCATE(CN_som_litter_longterm(kjpindex,nvm,nbpools),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for CN_som_litter_longterm. We stop. We need kjpindex*nvm*nbpools words',  & 
       &     kjpindex, nvm, nbpools
       STOP 'stomate_init'
    ENDIF
    
    ALLOCATE(KF(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for KF. We stop. We need nvm words = ',kjpindex*nvm
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF
    KF(:,:) = zero ! Is there a better place in the code for this?

    ALLOCATE(k_latosa_adapt(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for k_latosa_adapt. We stop. We need nvm words = ',kjpindex*nvm
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF

    ALLOCATE(harvest_pool_acc(kjpindex,nvm,ndia_harvest+1,nelements,nlanduse),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for harvest_pool_acc. We stop. We need many words = ',&
            kjpindex*nvm*(ndia_harvest+1)*nelements*nlanduse
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF
    harvest_pool_acc(:,:,:,:,:) = zero

    ALLOCATE(harvest_type(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for harvest_type. We stop. We need many words = ',&
            kjpindex*nvm
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF
    harvest_type(:,:) = zero

    ALLOCATE(harvest_cut(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for harvest_cut. We stop. We need many words = ',&
            kjpindex*nvm
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF
    harvest_cut(:,:) = zero

    ALLOCATE(harvest_area_acc(kjpindex,nvm,nlanduse),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for harvest_area_acc. We stop. We need many words = ',&
            kjpindex*nvm*nlanduse
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF
    harvest_area_acc(:,:,:) = zero

    ALLOCATE(gap_area_save(kjpindex,nvm,legacy_years_wind),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for gap_area_save. We stop. We need many words = ',&
            kjpindex*nvm*legacy_years_wind
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF

    ALLOCATE(total_ba_init(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for total_ba_init. We stop. We need many words = ',&
            kjpindex*nvm
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF

    ALLOCATE(harvest_pool_bound(0:ndia_harvest+1),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for harvest_pool_bound. ' // &
            'We stop. We need ndia_harvest+2 words = ',&
            ndia_harvest+2
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF

    ! Here we can initialize the values of this array, too. They 
    ! should never change over the course of the simulation.
    harvest_pool_bound(ndia_harvest+1) = val_exp
    DO idia = 0,ndia_harvest
       harvest_pool_bound(idia) = max_harvest_dia * &
            REAL(idia,r_std) / REAL(ndia_harvest,r_std)
    ENDDO

    ALLOCATE(sumTeff(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF
    sumTeff(:,:)=zero


    ALLOCATE(beetle_diapause(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF
    beetle_diapause(:,:)=zero

    ALLOCATE(mai(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for mai. ' // &
            'We stop. We need kjpindex*nvm words = ',&
            kjpindex*nvm
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF

    ALLOCATE(pai(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for pai. ' // &
            'We stop. We need kjpindex*nvm words = ',&
            kjpindex*nvm
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF

    ALLOCATE(previous_wood_volume(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for previous_wood_volume. ' // &
            'We stop. We need kjpindex*nvm words = ',&
            kjpindex*nvm
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF
    
    ALLOCATE(mai_count(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for mai_count. We stop. We need kjpindex*nvm words',  & 
       &     kjpindex, nvm
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF

    ALLOCATE(coppice_dens(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for coppice_dens. ' // &
            'We stop. We need kjpindex*nvm words = ',&
            kjpindex*nvm
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF

    ALLOCATE (rue_longterm(kjpindex,nvm), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for rue_longterm. We stop. We need kjpindex*nlevs words',kjpindex,nvm
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF
    rue_longterm(:,:) = un

    ALLOCATE (leaf_age_crit(kjpindex,nvm), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for leaf_age_crit. We stop. We need kjpindex*nlevs words',kjpindex,nvm
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF

    ALLOCATE (leaf_classes(kjpindex,nvm), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for leaf_classes. We stop. We need kjpindex*nlevs words',kjpindex,nvm
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF

    ALLOCATE(forest_managed(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for forest_managed. We stop. We need kjpindex*nvm words',  & 
       &     kjpindex, nvm
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF

    ALLOCATE(species_change_map(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for species_change_map. We stop. We need kjpindex*nvm words',  & 
       &     kjpindex, nvm
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF
    species_change_map(:,:)=0

    ALLOCATE(fm_change_map(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for fm_change_map. We stop. We need kjpindex*nvm words',  & 
       &     kjpindex, nvm
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF
    fm_change_map(:,:)=0

    ALLOCATE(lpft_replant(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for lpft_replant. We stop. We need kjpindex*nvm words',  & 
       &     kjpindex, nvm
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF
    lpft_replant(:,:)=.FALSE.

    ALLOCATE(age_stand(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for age_stand. We stop. We need kjpindex*nvm words',  &
       &     kjpindex, nvm
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF

    ! Guillaume M. -- Conserved biomass-weighted stand age. Always allocated (cheap) so that
    ! restart and deallocation stay unconditional; the diagnostic is always computed.
    ALLOCATE(age_stand_bm(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for age_stand_bm. We stop. We need kjpindex*nvm words',  &
       &     kjpindex, nvm
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF

    ALLOCATE(age_stand_area(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for age_stand_area. We stop. We need kjpindex*nvm words',  &
       &     kjpindex, nvm
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF

    ALLOCATE(rotation_n(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for rotation_n. We stop. We need kjpindex*nvm words',  & 
       &     kjpindex, nvm
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF

    ALLOCATE(last_cut(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for last_cut. We stop. We need kjpindex*nvm words',  & 
       &     kjpindex, nvm
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF

    ALLOCATE(sigma(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for sigma. We stop. We need kjpindex*nvm words',  & 
       &     kjpindex, nvm
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF

    ALLOCATE(wstress_season(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for wstress_season. ' // &
            'We stop. We need kjpindex*nvm words = ',&
            kjpindex, nvm
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF

    ALLOCATE(wstress_month(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for wstress_month. ' // &
            'We stop. We need kjpindex*nvm words = ',&
            kjpindex, nvm
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF

    ALLOCATE (deepSOM_a(kjpindex, ngrnd,nvm,nelements), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'stomate_init', 'Pb in alloc for deepSOM_a','','')
    
    ALLOCATE (deepSOM_s(kjpindex, ngrnd,nvm,nelements), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'stomate_init', 'Pb in alloc for deepSOM_s','','')
    
    ALLOCATE (deepSOM_p(kjpindex, ngrnd,nvm,nelements), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'stomate_init', 'Pb in alloc for deepSOM_p','','')
    
    ALLOCATE (altmax_ind(kjpindex,nvm), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'stomate_init', 'Pb in alloc for altmax_ind','','')

    ALLOCATE (altmax_lastyear(kjpindex,nvm), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'stomate_init', 'Pb in alloc for altmax_lastyear','','')

    ALLOCATE (altmax_ind_lastyear(kjpindex,nvm), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'stomate_init', 'Pb in alloc for altmax_ind_lastyear','','')
 
    ALLOCATE (O2_soil(kjpindex, ngrnd,nvm), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'stomate_init', 'Pb in alloc for O2_soil','','')
    
    ALLOCATE (CH4_soil(kjpindex, ngrnd,nvm), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'stomate_init', 'Pb in alloc for CH4_soil','','')
    
    ALLOCATE (O2_snow(kjpindex, nsnow,nvm), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'stomate_init', 'Pb in alloc for O2_snow','','')
    
    ALLOCATE (CH4_snow(kjpindex, nsnow,nvm), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'stomate_init', 'Pb in alloc for CH4_snow','','')
    
    ALLOCATE (tdeep_daily(kjpindex, ngrnd,nvm), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'stomate_init', 'Pb in alloc for tdeep_daily','','')
    
    ALLOCATE (fbact(kjpindex, ngrnd,nvm), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'stomate_init', 'Pb in alloc for fbact','','')

    ALLOCATE (decomp_rate(kjpindex, ngrnd,nvm), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'stomate_init', 'Pb in alloc for decomp_rate','','')
    decomp_rate=0.0
    
    ALLOCATE (decomp_rate_daily(kjpindex, ngrnd,nvm), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'stomate_init', 'Pb in alloc for decomp_rate_daily','','')
    
    ALLOCATE (hsdeep_daily(kjpindex, ngrnd,nvm), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'stomate_init', 'Pb in alloc for hsdeep_daily','','')
    
    ALLOCATE (temp_sol_daily(kjpindex), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'stomate_init', 'Pb in alloc for temp_sol_daily','','')
    
    ALLOCATE (snow_daily(kjpindex), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'stomate_init', 'Pb in alloc for snow_daily','','')

    ALLOCATE (pb_pa_daily(kjpindex), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'stomate_init', 'Pb in alloc for pb_pa_daily','','')
    
    ALLOCATE(fixed_cryoturbation_depth(kjpindex,nvm),stat=ier )
    IF (ier /= 0) CALL ipslerr_p(3,'stomate_init', 'Pb in alloc for fixed_cryoturbation_depth','','')
    
    ALLOCATE (snowdz_daily(kjpindex,nsnow), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'stomate_init', 'Pb in alloc for snowdz_daily','','')
    
    ALLOCATE (snowrho_daily(kjpindex,nsnow), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'stomate_init', 'Pb in alloc for snowrho_daily','','')    

    tdeep_daily=zero
    hsdeep_daily=zero
    decomp_rate_daily=zero
    snow_daily=zero
    pb_pa_daily=zero
    temp_sol_daily=zero
    snowdz_daily=zero
    snowrho_daily=zero

    ALLOCATE (bm_sapl_2D(kjpindex,nvm,ncirc,nparts,nelements), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN 
       WRITE(numout,*) 'Memory allocation error for bm_sapl_2D. We stop. ',kjpindex,nvm,nparts,nelements
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF
    bm_sapl_2D(:,:,:,:,:) = zero

    ALLOCATE(sugar_load(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for sugar_load. ' // &
            'We stop. We need kjpindex*nvm words = ',&
            kjpindex, nvm
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF
    sugar_load(:,:) = un

    ALLOCATE(grow_season_len(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for grow_season_len. ' // &
            'We stop. We need kjpindex*nvm words = ',&
            kjpindex, nvm
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF
    
    ALLOCATE(doy_start_gs(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for doy_start_gs. ' // &
            'We stop. We need kjpindex*nvm words = ',&
            kjpindex, nvm
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF

    ALLOCATE(doy_end_gs(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for doy_end_gs. ' // &
            'We stop. We need kjpindex*nvm words = ',&
            kjpindex, nvm
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF

    ALLOCATE(mean_start_gs(kjpindex,nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for mean_start_gs. ' // &
            'We stop. We need kjpindex*nvm words = ',&
            kjpindex, nvm
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF

    ALLOCATE(est_co2flux_d(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for est_co2flux_d. ' // &
            'We stop. We need kjpindex words = ',&
            kjpindex
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF
    
    ALLOCATE(ni_acc(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for ni_acc. ' // &
            'We stop. We need kjpindex words = ',&
            kjpindex, nvm
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF
    ni_acc(:) = zero

    ALLOCATE(lightning(kjpindex,12),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for lightning. ' // &
            'We stop. We need kjpindex words = ',&
            kjpindex, nvm
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF
    lightning(:,:) = zero

    ALLOCATE(pop_dens(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for pop_dens. ' // &
            'We stop. We need kjpindex words = ',&
            kjpindex, nvm
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF

    ALLOCATE(a_nd(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for a_nd. ' // &
            'We stop. We need kjpindex words = ',&
            kjpindex, nvm
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF

    ALLOCATE(edge_length(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for edge_length. ' // &
            'We stop. We need kjpindex words = ',&
            kjpindex, nvm
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF
    ! Guillaume M. -- Default to zero so a plain fire run (OK_SPITFIRE=y, all OK_AED_* off),
    ! which does not read Edge_length.nc, gets AED_fire = 50000 m in calculate_aed (no edge
    ! effect, historical SPITFIRE default) instead of an undefined value. The Edge_length.nc
    ! read overwrites it when any OK_AED_* flag is active, so AED runs are unaffected.
    edge_length(:) = zero

    ! Guillaume M. -- AED_FEEDBACK: allocate the dynamic edge_length state and the per-step
    ! perturbation accumulators owned by stomate_data. Init to zero; edge_length_dyn and
    ! edge_length_ref are populated from edge_length(:) on the first Edge_length.nc read.
    ALLOCATE(edge_length_dyn(kjpindex), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    ALLOCATE(edge_length_ref(kjpindex), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    ALLOCATE(dA_fire(kjpindex),  stat=ier)
    l_error = l_error .OR. (ier /= 0)
    ALLOCATE(dA_storm(kjpindex), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    ALLOCATE(dA_pest(kjpindex),  stat=ier)
    l_error = l_error .OR. (ier /= 0)
    ALLOCATE(dA_harvest(kjpindex), stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for AED_FEEDBACK arrays. ' // &
            'kjpindex = ', kjpindex
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF
    edge_length_dyn(:) = zero
    edge_length_ref(:) = zero
    dA_fire(:)    = zero
    dA_storm(:)   = zero
    dA_pest(:)    = zero
    dA_harvest(:) = zero

    ! Guillaume M. -- Management intensity: one allocation for the whole concept, the class
    ! fractions read annually plus the three fields set_management_intensity derives from
    ! them. mi_frac starts at zero => a cell is unmanaged until the file is read (no harvest
    ! edge, rotation = class 1 = never).

    ! Guillaume M. -- AED_EDGE_RDI_WEIGHT is allocated as soon as the weight is requested,
    ! outside the ok_management_intensity block: the two concepts are independent.
    IF (ok_aed_edge_rdi_weight) THEN
       ALLOCATE(rdi_stand(kjpindex, nvm), stat=ier)
       IF (ier /= 0) THEN
          WRITE(numout,*) 'Memory allocation error for rdi_stand. We stop. kjpindex x nvm = ',&
               kjpindex, nvm
          CALL ipslerr_p(3,'stomate_init','Memory allocation error for rdi_stand','','')
       ENDIF
       rdi_stand(:,:) = zero
    ENDIF

    ! Guillaume M. -- MACRO edge term: road length per grid cell. Allocated as soon as the
    ! file is provided, independently of ok_management_intensity -- the two maps share the
    ! reader and nothing else.
    IF (TRIM(road_length_file) /= 'NONE') THEN
       ALLOCATE(road_length_map(kjpindex), stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'stomate_init','Pb alloc road_length_map','','')
       road_length_map(:) = zero
    ENDIF

    ! Guillaume M. -- FRAG_3TERMS permanent: lake, coast and river edge. Same pattern as the
    ! road map -- static, read once, added as is.
    IF (TRIM(edge_hydro_file) /= 'NONE') THEN
       ALLOCATE(edge_hydro_map(kjpindex), stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'stomate_init','Pb alloc edge_hydro_map','','')
       edge_hydro_map(:) = zero
    ENDIF

    IF (ok_management_intensity) THEN
       ! Guillaume M. -- ROTATION_GROWTH output, allocated with the management it blends
       ! into (it used to hide under OK_AED_EDGE_RDI_WEIGHT, an unrelated flag). The blend
       ! itself is gated by rotation_growth_weight. mat_age_entry is allocated only when a
       ! consumer asks for it, so that the default restart keeps its variable list.
       ALLOCATE(rotation_growth(kjpindex, nvm), stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (ier == 0) rotation_growth(:,:) = zero
       IF (mat_a3_online .OR. rotation_growth_weight > zero) THEN
          ALLOCATE(mat_age_entry(kjpindex, nvm), stat=ier)
          l_error = l_error .OR. (ier /= 0)
          IF (ier == 0) mat_age_entry(:,:) = -un
       ENDIF
       ALLOCATE(mi_frac(kjpindex, nmiclass), stat=ier)
       l_error = l_error .OR. (ier /= 0)
       ALLOCATE(clearcut_area(kjpindex), stat=ier)
       l_error = l_error .OR. (ier /= 0)
       ALLOCATE(alpha_harvest(kjpindex), stat=ier)
       l_error = l_error .OR. (ier /= 0)
       ALLOCATE(dia_factor(kjpindex), stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'stomate_init','Pb alloc dia_factor','','')
       dia_factor(:) = un
       ALLOCATE(target_rotation_age(kjpindex,nvm), stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for management-intensity arrays. ' // &
               'kjpindex = ', kjpindex, ' nmiclass = ', nmiclass
          CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
       ENDIF
       mi_frac(:,:)                = zero
       clearcut_area(:)            = mi_clearcut_size_max
       alpha_harvest(:)        = zero
       ! Guillaume M. -- No management map here, so every stand takes the unmanaged
       ! sentinel: rotation_ref scaled by the unmanaged intensity factor. Zero elsewhere
       ! means "no rotation defined"; only tree PFTs ever read this map.
       target_rotation_age(:,:)  = zero
       DO j = 1, nvm
          IF (rotation_ref(j) > zero) &
               target_rotation_age(:,j) = rotation_ref(j) * mi_rotation_factor(1)
       ENDDO
    ENDIF

    ! Guillaume M. -- AED_REGROWTH class-0 conveyor: disturbed-area cohorts by agent and
    ! age, allocated under OK_EDGE_FROM_AGE_CLASS. It never touches veget_max; it is a pure
    ! area register feeding the edge length and the FRAC_CLASS0 diagnostic.

    ! Guillaume M. -- The maturity conveyor is allocated on N_CLASS0 ALONE. It is a
    ! different conveyor from the edge one: indexed by GROUP, not by agent, and it moves
    ! real area. Nesting it inside the edge flag made N_CLASS0 silently inert whenever
    ! that flag was off, with no message.
    IF (n_class0 > 0) THEN
       ALLOCATE(class0_veg(kjpindex, nvmap, n_class0), stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'stomate_init','Pb alloc class0_veg','','')
       class0_veg(:,:,:) = zero
       ! Guillaume M. -- Allocated HERE and not under ok_edge_from_age_class: same reason as
       ! class0_veg above, this conveyor is gated on n_class0 and nothing else.
       ALLOCATE(class0_restitute(kjpindex, nvm), stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'stomate_init','Pb alloc class0_restitute','','')
       class0_restitute(:,:) = zero
    ENDIF

    IF (ok_edge_from_age_class) THEN
       ALLOCATE(class0_area(kjpindex, nagent, n_recruit_max), stat=ier)
       l_error = l_error .OR. (ier /= 0)
       ALLOCATE(class0_recovered(kjpindex), stat=ier)
       l_error = l_error .OR. (ier /= 0)
       ALLOCATE(agent_flux(kjpindex, nagent), stat=ier)   ! AED_REGROWTH v2 (A_open agent mix)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for AED_REGROWTH class0 arrays. ' // &
               'kjpindex = ', kjpindex, ' nagent = ', nagent, ' n_recruit_max = ', n_recruit_max
          CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
       ENDIF
       class0_area(:,:,:)  = zero
       class0_recovered(:) = zero
       class0_year_clock_s = zero
       agent_flux(:,:)     = zero
    ENDIF

    ! Guillaume M. -- There is no automatic AED spin-up convergence test: convergence is
    ! judged on the diagnostic figures.

    ALLOCATE(daily_wspeed_fire(kjpindex),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) 'Memory allocation error for daily_wspeed_fire. We stop. We need kjpindex words',kjpindex
       CALL ipslerr_p (3,'stomate_init', 'Memory allocation issue','','')
    ENDIF
    
    ! bark beetle module. Initialize the variable here so they all have a value 
    ! of zero when the bark beetle module is not used. 
    ALLOCATE (wood_leftover_legacy(kjpindex,nvm,legacy_years_wood),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'stomate_init','Pb in alloc for wood_leftover_legacy','','')
    wood_leftover_legacy(:,:,:) = zero

    ALLOCATE (season_drought_legacy(kjpindex,nvm,legacy_years),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'stomate_init','Pb in alloc for season_drought_legacy','','')
    season_drought_legacy(:,:,:) = zero

    ALLOCATE (i_beetles_activity_legacy(kjpindex,nvm,legacy_years),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'stomate_init','Pb in alloc for i_beetles_activity_legacy','','')
    i_beetles_activity_legacy(:,:,:) = zero

    ALLOCATE (B_beetles_kill_legacy(kjpindex,nvm,beetle_legacy),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'stomate_init','Pb in alloc for B_beetles_kill_legacy','','')
    B_beetles_kill_legacy(:,:,:) = zero

    ALLOCATE (i_beetles_generation(kjpindex,nvm,legacy_years),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'stomate_init','Pb in alloc for beetle_generation_index','','')
    i_beetles_generation(:,:,:) = zero

    ALLOCATE (P_beetles_attacked_legacy(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'stomate_init','Pb in alloc for risk_index','','')
    P_beetles_attacked_legacy(:,:) = zero

    ALLOCATE (woody_litter_to_use(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'stomate_init','Pb in alloc for woody_litter_to_use','','')
    woody_litter_to_use(:,:) = zero

    ALLOCATE (mbc_daily(kjpindex,nelements),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'stomate_init','Pb in alloc for mbc_daily','','')
    mbc_daily(:,:) = zero


  !! 5. File definitions

    ! Store history and restart files in common variables
    hist_id_stomate = hist_id_stom
    hist_id_stomate_IPCC = hist_id_stom_IPCC
    rest_id_stomate = rest_id_stom
    
    ! In STOMATE reduced grids are used containing only terrestrial pixels.
    ! Build a new indexing table for the vegetation fields separating 
    ! between the different PFTs. Note that ::index has dimension (kjpindex) 
    ! wheras ::indexpft has dimension (kjpindex*nvm). 

    hori_index(:) = index(:)

    DO j = 1, nvm
       DO ji = 1, kjpindex
          horipft_index((j-1)*kjpindex+ji) = index(ji)+(j-1)*kjpij + offset_omp - offset_mpi
       ENDDO
    ENDDO

     DO j = 1, nlevels_tot
       DO ji = 1, kjpindex
          horican_index((j-1)*kjpindex+ji) = index(ji)+(j-1)*kjpij + offset_omp - offset_mpi
       ENDDO
    ENDDO

    DO j = 1, ncut_times
       DO ji = 1, kjpindex
          horicut_index((j-1)*kjpindex+ji) = index(ji)+(j-1)*kjpij + offset_omp - offset_mpi
       ENDDO
    ENDDO

    ! Similar index tables are build for the wood use
    DO j = 1, nshort
       DO ji = 1, kjpindex
          horip_s_index((j-1)*kjpindex+ji) = &
               index(ji)+(j-1)*kjpij + offset_omp - offset_mpi
       ENDDO
    ENDDO

    DO j = 1, nmedium
       DO ji = 1, kjpindex
          horip_m_index((j-1)*kjpindex+ji) = &
               index(ji)+(j-1)*kjpij + offset_omp - offset_mpi
       ENDDO
    ENDDO
    
    DO j = 1, nlong
       DO ji = 1, kjpindex
          horip_l_index((j-1)*kjpindex+ji) = &
               index(ji)+(j-1)*kjpij + offset_omp - offset_mpi
       ENDDO
    ENDDO

    DO j = 1, nshort+1
       DO ji = 1, kjpindex
          horip_ss_index((j-1)*kjpindex+ji) = &
               index(ji)+(j-1)*kjpij + offset_omp - offset_mpi
       ENDDO
    ENDDO

    DO j = 1, nmedium+1
       DO ji = 1, kjpindex
          horip_mm_index((j-1)*kjpindex+ji) = &
               index(ji)+(j-1)*kjpij + offset_omp - offset_mpi
       ENDDO
    ENDDO

    DO j = 1, nlong+1
       DO ji = 1, kjpindex
          horip_ll_index((j-1)*kjpindex+ji) = &
               index(ji)+(j-1)*kjpij + offset_omp - offset_mpi
       ENDDO
    ENDDO
  

  !! 6. Initialization of global and land cover change variables. 

    ! All variables are cumulative variables. bm_to_litter is not and is therefore
    ! excluded
    turnover_daily(:,:,:,:) = zero
    resp_hetero_d(:,:) = zero
    resp_hetero_litter_d(:,:) = zero
    resp_hetero_soil_d(:,:) = zero
    som_input_daily(:,:,:,:) = zero
    drainage_daily(:,:) = zero 
    atm_to_immob_daily(:,:) = zero
    leaching_daily(:,:,:) = zero
    emission_daily(:,:,:) = zero
    n_input_daily(:,:,:) = zero
    woodharvestpft(:,:) = zero
    fpc_max(:,:)=zero
   
    ! n variables
    nstress_season(:,:) = zero 
    soil_n_min(:,:,:) = zero
    plant_n_uptake_daily(:,:,:) = zero 
    n_mineralisation_d(:,:) = zero 

    fDeforestToProduct(:,:)=zero
    fLulccResidue(:,:)=zero
    fHarvestToProduct(:,:)=zero

    ! peatland
    ALLOCATE (wtp_daily(kjpindex), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'stomate_init','Pb in alloc for wtp_daily','','')
    ALLOCATE(carbon_acro(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'stomate_init','Pb in alloc for carbon_acro','','')
    ALLOCATE(carbon_cato(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'stomate_init','Pb in alloc for carbon_cato','','')
    ALLOCATE(height_acro(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'stomate_init','Pb in alloc for height_acro','','')
    ALLOCATE(height_cato(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'stomate_init','Pb in alloc for height_cato','','')
    ALLOCATE(tcarbon_acro(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'stomate_init','Pb in alloc for tcarbon_acro','','')
    ALLOCATE(tcarbon_cato(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'stomate_init','Pb in alloc for tcarbon_cato','','')
    ALLOCATE(resp_acro_oxic(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'stomate_init','Pb in alloc for resp_acro_oxic','','')
    ALLOCATE(resp_acro_anoxic(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'stomate_init','Pb in alloc for resp_acro_anoxic','','')
    ALLOCATE(resp_cato(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'stomate_init','Pb in alloc for resp_cato','','')
    ALLOCATE(litter_to_acro(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'stomate_init','Pb in alloc for litter_to_acro','','')
    ALLOCATE(acro_to_cato(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'stomate_init','Pb in alloc for acro_to_cato','','')
    ALLOCATE (deepSOM_peat(kjpindex, ngrnd, nvm, nelements), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'stomate_init','Pb in alloc for deepSOM_peat','','')

  END SUBROUTINE stomate_init


!! ================================================================================================================================
!! SUBROUTINE 	: stomate_clear
!!
!>\BRIEF        Deallocate memory of the stomate variables.
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCES	: None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  
  SUBROUTINE stomate_clear

  !! Deallocate all dynamics variables
    IF (ALLOCATED(adapted)) DEALLOCATE(adapted)
    IF (ALLOCATED(regenerate)) DEALLOCATE(regenerate)
    IF (ALLOCATED(vegstress_day)) DEALLOCATE(vegstress_day)
    IF (ALLOCATED(transpir_daily)) DEALLOCATE(transpir_daily)
    IF (ALLOCATED(gdd_init_date)) DEALLOCATE(gdd_init_date)
    IF (ALLOCATED(litterhum_daily)) DEALLOCATE(litterhum_daily)
    IF (ALLOCATED(t2m_daily))  DEALLOCATE(t2m_daily)
    IF (ALLOCATED(t2m_min_daily))  DEALLOCATE(t2m_min_daily)
    IF (ALLOCATED(t2m_max_daily))  DEALLOCATE(t2m_max_daily)
    IF (ALLOCATED(tsurf_daily))  DEALLOCATE(tsurf_daily)
    IF (ALLOCATED(tsoil_daily)) DEALLOCATE(tsoil_daily)
    IF (ALLOCATED(precip_daily)) DEALLOCATE(precip_daily)
    IF (ALLOCATED(vpd_daily_mean)) DEALLOCATE(vpd_daily_mean)
    IF (ALLOCATED(vpd_daily_max)) DEALLOCATE(vpd_daily_max)
    IF (ALLOCATED(vpd_mean_week)) DEALLOCATE(vpd_mean_week)
    IF (ALLOCATED(vpd_max_week)) DEALLOCATE(vpd_max_week)
    IF (ALLOCATED(gpp_daily)) DEALLOCATE(gpp_daily)
    IF (ALLOCATED(npp_daily)) DEALLOCATE(npp_daily)
    IF (ALLOCATED(turnover_daily)) DEALLOCATE(turnover_daily)
    IF (ALLOCATED(turnover_resid)) DEALLOCATE(turnover_resid)
    IF (ALLOCATED(turnover_littercalc)) DEALLOCATE(turnover_littercalc)
    IF (ALLOCATED(vegstress_month)) DEALLOCATE(vegstress_month)
    IF (ALLOCATED(vegstress_week)) DEALLOCATE(vegstress_week)
    IF (ALLOCATED(vegstress_season)) DEALLOCATE(vegstress_season)
    IF (ALLOCATED(t2m_longterm)) DEALLOCATE(t2m_longterm)
    IF (ALLOCATED(t2m_month)) DEALLOCATE(t2m_month)
    IF (ALLOCATED(Tseason)) DEALLOCATE(Tseason)
    IF (ALLOCATED(Tseason_length)) DEALLOCATE(Tseason_length)
    IF (ALLOCATED(Tseason_tmp)) DEALLOCATE(Tseason_tmp)
    IF (ALLOCATED(Tmin_spring_time)) DEALLOCATE(Tmin_spring_time)
    IF (ALLOCATED(t2m_week)) DEALLOCATE(t2m_week)
    IF (ALLOCATED(tsoil_month)) DEALLOCATE(tsoil_month)
    IF (ALLOCATED(precip_month)) DEALLOCATE(precip_month)
    IF (ALLOCATED(fireindex)) DEALLOCATE(fireindex)
    IF (ALLOCATED(firelitter)) DEALLOCATE(firelitter)
    IF (ALLOCATED(maxvegstress_lastyear)) DEALLOCATE(maxvegstress_lastyear)
    IF (ALLOCATED(maxvegstress_thisyear)) DEALLOCATE(maxvegstress_thisyear)
    IF (ALLOCATED(minvegstress_lastyear)) DEALLOCATE(minvegstress_lastyear)
    IF (ALLOCATED(minvegstress_thisyear)) DEALLOCATE(minvegstress_thisyear)
    IF (ALLOCATED(maxgppweek_lastyear)) DEALLOCATE(maxgppweek_lastyear)
    IF (ALLOCATED(maxgppweek_thisyear)) DEALLOCATE(maxgppweek_thisyear)
    IF (ALLOCATED(gdd0_lastyear)) DEALLOCATE(gdd0_lastyear)
    IF (ALLOCATED(gdd0_thisyear)) DEALLOCATE(gdd0_thisyear)
    IF (ALLOCATED(precip_lastyear)) DEALLOCATE(precip_lastyear)
    IF (ALLOCATED(precip_thisyear)) DEALLOCATE(precip_thisyear)
    IF (ALLOCATED(gdd_m5_dormance)) DEALLOCATE(gdd_m5_dormance)
    IF (ALLOCATED(gdd_from_growthinit)) DEALLOCATE(gdd_from_growthinit)
    IF (ALLOCATED(gdd_midwinter)) DEALLOCATE(gdd_midwinter)
    IF (ALLOCATED(ncd_dormance)) DEALLOCATE(ncd_dormance)
    IF (ALLOCATED(ngd_minus5))  DEALLOCATE(ngd_minus5)
    IF (ALLOCATED(PFTpresent)) DEALLOCATE(PFTpresent)
    IF (ALLOCATED(wind_ratio_sum)) DEALLOCATE(wind_ratio_sum)
    IF (ALLOCATED(wind_sum)) DEALLOCATE(wind_sum)
    IF (ALLOCATED(max_wind_speed_storm)) DEALLOCATE(max_wind_speed_storm)
    IF (ALLOCATED(max_wind_ratio_storm)) DEALLOCATE(max_wind_ratio_storm)  ! JJ 2026
    IF (ALLOCATED(wind_ratio_sum_save)) DEALLOCATE(wind_ratio_sum_save)
    IF (ALLOCATED(wind_ratio_max_save)) DEALLOCATE(wind_ratio_max_save)  ! JJ 2026
    IF (ALLOCATED(wind_speed_max_save)) DEALLOCATE(wind_speed_max_save)
    IF (ALLOCATED(wind_mean_save)) DEALLOCATE(wind_mean_save)
    IF (ALLOCATED(wind_max)) DEALLOCATE(wind_max)
    IF (ALLOCATED(wind_ratio_max)) DEALLOCATE(wind_ratio_max)  ! JJ 2026
    IF (ALLOCATED(is_storm)) DEALLOCATE(is_storm)
    IF (ALLOCATED(count_storm)) DEALLOCATE(count_storm)
    IF (ALLOCATED(npp_longterm)) DEALLOCATE(npp_longterm)
    IF (ALLOCATED(croot_longterm)) DEALLOCATE(croot_longterm)
    IF (ALLOCATED(n_reserve_longterm)) DEALLOCATE(n_reserve_longterm)
    IF (ALLOCATED(lm_lastyearmax)) DEALLOCATE(lm_lastyearmax)
    IF (ALLOCATED(lm_thisyearmax)) DEALLOCATE(lm_thisyearmax)
    IF (ALLOCATED(maxfpc_lastyear)) DEALLOCATE(maxfpc_lastyear)
    IF (ALLOCATED(maxfpc_thisyear)) DEALLOCATE(maxfpc_thisyear)
    IF (ALLOCATED(turnover_longterm)) DEALLOCATE(turnover_longterm)
    IF (ALLOCATED(gpp_week)) DEALLOCATE(gpp_week)
    IF (ALLOCATED(gpp_year)) DEALLOCATE(gpp_year)
    IF (ALLOCATED(gpp_decade)) DEALLOCATE(gpp_decade)
    IF (ALLOCATED(resp_maint_week)) DEALLOCATE(resp_maint_week)
    IF (ALLOCATED(plant_status)) DEALLOCATE(plant_status)
    IF (ALLOCATED(when_growthinit)) DEALLOCATE(when_growthinit)
    IF (ALLOCATED(age))  DEALLOCATE(age)
    IF (ALLOCATED(resp_hetero_d)) DEALLOCATE(resp_hetero_d)
    IF (ALLOCATED(resp_hetero_litter_d)) DEALLOCATE(resp_hetero_litter_d)
    IF (ALLOCATED(resp_hetero_soil_d)) DEALLOCATE(resp_hetero_soil_d)
    IF (ALLOCATED(resp_hetero_radia)) DEALLOCATE(resp_hetero_radia)
    IF (ALLOCATED(resp_maint_d)) DEALLOCATE(resp_maint_d)
    IF (ALLOCATED(resp_growth_d)) DEALLOCATE(resp_growth_d)
    IF (ALLOCATED(co2_fire)) DEALLOCATE(co2_fire)
    IF (ALLOCATED(emissions_fire)) DEALLOCATE(emissions_fire)
    IF (ALLOCATED(atm_to_immob)) DEALLOCATE(atm_to_immob)
    IF (ALLOCATED(veget_lastlight)) DEALLOCATE(veget_lastlight)
    IF (ALLOCATED(everywhere)) DEALLOCATE(everywhere)
    IF (ALLOCATED(need_adjacent)) DEALLOCATE(need_adjacent)
    IF (ALLOCATED(leaf_age)) DEALLOCATE(leaf_age)
    IF (ALLOCATED(leaf_frac)) DEALLOCATE(leaf_frac)
    IF (ALLOCATED(RIP_time)) DEALLOCATE(RIP_time)
    IF (ALLOCATED(time_hum_min)) DEALLOCATE(time_hum_min)
    IF (ALLOCATED(hum_min_dormance)) DEALLOCATE(hum_min_dormance)
    IF (ALLOCATED(litter)) DEALLOCATE(litter)
    IF (ALLOCATED(dead_leaves)) DEALLOCATE(dead_leaves)
    IF (ALLOCATED(som)) DEALLOCATE(som)
    IF (ALLOCATED(som_surf)) DEALLOCATE(som_surf)
    IF (ALLOCATED(lignin_struc)) DEALLOCATE(lignin_struc)
    IF (ALLOCATED(burried_litter)) DEALLOCATE(burried_litter)
    IF (ALLOCATED(burried_fresh_ltr)) DEALLOCATE(burried_fresh_ltr)
    IF (ALLOCATED(burried_fresh_som)) DEALLOCATE(burried_fresh_som)
    IF (ALLOCATED(burried_bact)) DEALLOCATE(burried_bact)
    IF (ALLOCATED(burried_min_nitro)) DEALLOCATE(burried_min_nitro)
    IF (ALLOCATED(burried_som)) DEALLOCATE(burried_som)
    IF (ALLOCATED(burried_deepSOM_a)) DEALLOCATE(burried_deepSOM_a)
    IF (ALLOCATED(burried_deepSOM_s)) DEALLOCATE(burried_deepSOM_s)
    IF (ALLOCATED(burried_deepSOM_p)) DEALLOCATE(burried_deepSOM_p)
    IF (ALLOCATED(lignin_wood)) DEALLOCATE(lignin_wood)
    IF (ALLOCATED(lignin_snag)) DEALLOCATE(lignin_snag)
    IF (ALLOCATED(turnover_time)) DEALLOCATE(turnover_time)
    IF (ALLOCATED(bm_to_litter)) DEALLOCATE(bm_to_litter)
    IF (ALLOCATED(bm_to_litter_resid)) DEALLOCATE(bm_to_litter_resid)
    IF (ALLOCATED(tree_bm_to_litter)) DEALLOCATE(tree_bm_to_litter)
    IF (ALLOCATED(tree_bm_to_litter_resid)) DEALLOCATE(tree_bm_to_litter_resid)
    IF (ALLOCATED(bm_to_littercalc)) DEALLOCATE(bm_to_littercalc)
    IF (ALLOCATED(tree_bm_to_littercalc)) DEALLOCATE(tree_bm_to_littercalc)
    IF (ALLOCATED(herbivores)) DEALLOCATE(herbivores)
    IF (ALLOCATED(resp_maint_part_radia)) DEALLOCATE(resp_maint_part_radia)
    IF (ALLOCATED(resp_maint_radia)) DEALLOCATE(resp_maint_radia)
    IF (ALLOCATED(resp_maint_part)) DEALLOCATE(resp_maint_part)
    IF (ALLOCATED(hori_index)) DEALLOCATE(hori_index)
    IF (ALLOCATED(horipft_index)) DEALLOCATE(horipft_index)
    IF (ALLOCATED(horican_index)) DEALLOCATE(horican_index)
    IF (ALLOCATED(horicut_index)) DEALLOCATE(horicut_index)
    IF (ALLOCATED(horip_s_index)) DEALLOCATE (horip_s_index)
    IF (ALLOCATED(horip_m_index)) DEALLOCATE (horip_m_index)
    IF (ALLOCATED(horip_l_index)) DEALLOCATE (horip_l_index)
    IF (ALLOCATED(horip_ss_index)) DEALLOCATE (horip_ss_index)
    IF (ALLOCATED(horip_mm_index)) DEALLOCATE (horip_mm_index)
    IF (ALLOCATED(horip_ll_index)) DEALLOCATE (horip_ll_index)
    !
    IF (ALLOCATED(ok_equilibrium)) DEALLOCATE(ok_equilibrium)
    IF (ALLOCATED(carbon_eq)) DEALLOCATE(carbon_eq)
    IF (ALLOCATED(matrixA)) DEALLOCATE(matrixA)
    IF (ALLOCATED(vectorB)) DEALLOCATE(vectorB)
    IF (ALLOCATED(matrixV)) DEALLOCATE(matrixV)
    IF (ALLOCATED(vectorU)) DEALLOCATE(vectorU)
    IF (ALLOCATED(matrixW)) DEALLOCATE(matrixW)
    IF (ALLOCATED(previous_stock)) DEALLOCATE(previous_stock)
    IF (ALLOCATED(current_stock)) DEALLOCATE(current_stock)
    IF (ALLOCATED(sigma)) DEALLOCATE (sigma) 
    IF (ALLOCATED(age_stand)) DEALLOCATE (age_stand)
    IF (ALLOCATED(age_stand_bm)) DEALLOCATE (age_stand_bm)
    IF (ALLOCATED(age_stand_area)) DEALLOCATE (age_stand_area)
    IF (ALLOCATED(rotation_n)) DEALLOCATE (rotation_n)
    IF (ALLOCATED(last_cut)) DEALLOCATE (last_cut)
    IF (ALLOCATED(CN_som_litter_longterm)) DEALLOCATE(CN_som_litter_longterm) 
    IF (ALLOCATED(KF)) DEALLOCATE (KF)
    IF (ALLOCATED(k_latosa_adapt)) DEALLOCATE (k_latosa_adapt)
    IF (ALLOCATED(harvest_pool_acc)) DEALLOCATE (harvest_pool_acc)
    IF (ALLOCATED(harvest_type)) DEALLOCATE (harvest_type)
    IF (ALLOCATED(harvest_cut)) DEALLOCATE (harvest_cut)
    IF (ALLOCATED(harvest_area_acc)) DEALLOCATE (harvest_area_acc)
    IF (ALLOCATED(gap_area_save)) DEALLOCATE (gap_area_save)
    IF (ALLOCATED(total_ba_init)) DEALLOCATE (total_ba_init)
    IF (ALLOCATED(harvest_pool_bound)) DEALLOCATE (harvest_pool_bound)
    IF (ALLOCATED(sumTeff)) DEALLOCATE (sumTeff)
    IF (ALLOCATED(beetle_diapause)) DEALLOCATE (beetle_diapause)
    IF (ALLOCATED(prod_s)) DEALLOCATE (prod_s)
    IF (ALLOCATED(prod_m)) DEALLOCATE (prod_m)
    IF (ALLOCATED(prod_l)) DEALLOCATE (prod_l)
    IF (ALLOCATED(flux_s)) DEALLOCATE (flux_s)
    IF (ALLOCATED(flux_m)) DEALLOCATE (flux_m)
    IF (ALLOCATED(flux_l)) DEALLOCATE (flux_l)
    IF (ALLOCATED(flux_prod_s)) DEALLOCATE (flux_prod_s)
    IF (ALLOCATED(flux_prod_m)) DEALLOCATE (flux_prod_m)
    IF (ALLOCATED(flux_prod_l)) DEALLOCATE (flux_prod_l)

    IF (ALLOCATED(mai)) DEALLOCATE (mai)
    IF (ALLOCATED(pai)) DEALLOCATE (pai)
    IF (ALLOCATED(previous_wood_volume)) DEALLOCATE (previous_wood_volume)
    IF (ALLOCATED(mai_count)) DEALLOCATE (mai_count)
    IF (ALLOCATED(coppice_dens)) DEALLOCATE (coppice_dens)
    IF (ALLOCATED(wstress_season)) DEALLOCATE (wstress_season)
    IF (ALLOCATED(wstress_month)) DEALLOCATE (wstress_month)
    IF (ALLOCATED(rue_longterm)) DEALLOCATE (rue_longterm)
    IF (ALLOCATED(bm_sapl_2D)) DEALLOCATE (bm_sapl_2D)
    IF (ALLOCATED(sugar_load)) DEALLOCATE (sugar_load)
    IF (ALLOCATED(nbp_accu_flux)) DEALLOCATE(nbp_accu_flux)
    IF (ALLOCATED(nbp_pool_start)) DEALLOCATE(nbp_pool_start)
    IF (ALLOCATED(nforce)) DEALLOCATE(nforce)
    IF (ALLOCATED(control_moist)) DEALLOCATE(control_moist)
    IF (ALLOCATED(control_temp)) DEALLOCATE(control_temp)
    IF (ALLOCATED(carbon_input)) DEALLOCATE(carbon_input)
    IF (ALLOCATED(nitrogen_input)) DEALLOCATE(nitrogen_input)
    IF ( ALLOCATED (fco2_lu)) DEALLOCATE (fco2_lu)
    IF ( ALLOCATED (fco2_wh)) DEALLOCATE (fco2_wh)
    IF ( ALLOCATED (fco2_ha)) DEALLOCATE (fco2_ha)
    IF ( ALLOCATED (woodharvestpft)) DEALLOCATE (woodharvestpft)
    IF ( ALLOCATED (fDeforestToProduct)) DEALLOCATE (fDeforestToProduct)
    IF ( ALLOCATED (fLulccResidue)) DEALLOCATE (fLulccResidue)
    IF ( ALLOCATED (fHarvestToProduct)) DEALLOCATE (fHarvestToProduct)
    IF ( ALLOCATED (som_input_daily)) DEALLOCATE (som_input_daily)

    IF ( ALLOCATED (drainage_daily)) DEALLOCATE(drainage_daily) 
    IF ( ALLOCATED (plant_n_uptake_daily)) DEALLOCATE(plant_n_uptake_daily) 
    IF ( ALLOCATED (n_mineralisation_d)) DEALLOCATE(n_mineralisation_d)
    IF ( ALLOCATED (atm_to_immob_daily)) DEALLOCATE(atm_to_immob_daily)
    IF ( ALLOCATED (emission_daily)) DEALLOCATE(emission_daily)
    IF ( ALLOCATED (leaching_daily)) DEALLOCATE(leaching_daily)
    IF ( ALLOCATED (n_input_daily)) DEALLOCATE(n_input_daily)
    IF ( ALLOCATED (cn_leaf_min_season)) DEALLOCATE (cn_leaf_min_season) 
    IF ( ALLOCATED (nstress_season)) DEALLOCATE (nstress_season) 
    IF ( ALLOCATED (soil_n_min)) DEALLOCATE (soil_n_min) 
    IF ( ALLOCATED (p_O2)) DEALLOCATE (p_O2) 
    IF ( ALLOCATED (bact)) DEALLOCATE (bact) 
    IF ( ALLOCATED (fpc_max)) DEALLOCATE (fpc_max)

    IF (ALLOCATED(forest_managed)) DEALLOCATE (forest_managed)
    IF (ALLOCATED(species_change_map)) DEALLOCATE (species_change_map)
    IF (ALLOCATED(fm_change_map)) DEALLOCATE (fm_change_map) 
    IF (ALLOCATED(lpft_replant)) DEALLOCATE (lpft_replant)
    IF (ALLOCATED(grow_season_len)) DEALLOCATE (grow_season_len)
    IF (ALLOCATED(doy_start_gs)) DEALLOCATE (doy_start_gs)
    IF (ALLOCATED(doy_end_gs)) DEALLOCATE (doy_end_gs)
    IF (ALLOCATED(mean_start_gs)) DEALLOCATE (mean_start_gs)
    IF (ALLOCATED(est_co2flux_d)) DEALLOCATE (est_co2flux_d)
    IF ( ALLOCATED (wood_leftover_legacy)) DEALLOCATE (wood_leftover_legacy)
    IF ( ALLOCATED (season_drought_legacy)) DEALLOCATE (season_drought_legacy)
    IF ( ALLOCATED (i_beetles_generation)) DEALLOCATE(i_beetles_generation)
    IF ( ALLOCATED (P_beetles_attacked_legacy))DEALLOCATE(P_beetles_attacked_legacy)
    IF ( ALLOCATED (B_beetles_kill_legacy))DEALLOCATE(B_beetles_kill_legacy)
    IF ( ALLOCATED (i_beetles_activity_legacy))DEALLOCATE(i_beetles_activity_legacy)
    IF ( ALLOCATED (woody_litter_to_use))DEALLOCATE(woody_litter_to_use)
    IF (ALLOCATED(litterfuel)) DEALLOCATE (litterfuel)
    IF (ALLOCATED(ni_acc))  DEALLOCATE (ni_acc)
    IF (ALLOCATED(lightning)) DEALLOCATE (lightning)
    IF (ALLOCATED(pop_dens)) DEALLOCATE (pop_dens)
    IF (ALLOCATED(a_nd)) DEALLOCATE (a_nd)
    IF (ALLOCATED(edge_length)) DEALLOCATE (edge_length)
    IF (ALLOCATED(daily_wspeed_fire)) DEALLOCATE (daily_wspeed_fire)
    ! Guillaume M. -- Fragmentation AED variables owned by stomate_data.
    IF (ALLOCATED(AED_fire))    DEALLOCATE (AED_fire)
    IF (ALLOCATED(AED_light))   DEALLOCATE (AED_light)
    IF (ALLOCATED(area_forest)) DEALLOCATE (area_forest)
    ! Guillaume M. -- AED_FEEDBACK dynamic state and accumulators.
    IF (ALLOCATED(edge_length_dyn)) DEALLOCATE (edge_length_dyn)
    IF (ALLOCATED(edge_length_ref)) DEALLOCATE (edge_length_ref)
    IF (ALLOCATED(dA_fire))    DEALLOCATE (dA_fire)
    IF (ALLOCATED(dA_storm))   DEALLOCATE (dA_storm)
    IF (ALLOCATED(dA_pest))    DEALLOCATE (dA_pest)
    IF (ALLOCATED(dA_harvest)) DEALLOCATE (dA_harvest)
    ! Guillaume M. -- AED_REGROWTH class-0 conveyor.
    IF (ALLOCATED(class0_area))      DEALLOCATE (class0_area)
    IF (ALLOCATED(class0_recovered)) DEALLOCATE (class0_recovered)
    IF (ALLOCATED(class0_restitute)) DEALLOCATE (class0_restitute)
    IF (ALLOCATED(agent_flux))       DEALLOCATE (agent_flux)         ! AED_REGROWTH v2
    IF (ALLOCATED(recruit_time_pft)) DEALLOCATE (recruit_time_pft)   ! Trecruit(PFT) 2026-06-18
    IF (ALLOCATED(n_recruit_eff))    DEALLOCATE (n_recruit_eff)
 !! 2. reset l_first

    l_first_stomate=.TRUE.

 !! 3. call to clear functions

    CALL season_pre_disturbance_clear
    CALL season_post_disturbance_clear
    CALL stomate_lpj_clear
    CALL littercalc_clear
    CALL vmax_clear
    CALL stomate_soil_carbon_discretization_clear

    IF ( ALLOCATED (deepSOM_a)) DEALLOCATE(deepSOM_a)
    IF ( ALLOCATED (deepSOM_s)) DEALLOCATE(deepSOM_s)
    IF ( ALLOCATED (deepSOM_p)) DEALLOCATE(deepSOM_p)
    IF ( ALLOCATED (altmax_ind)) DEALLOCATE(altmax_ind)
    IF ( ALLOCATED (altmax_lastyear)) DEALLOCATE(altmax_lastyear)
    IF ( ALLOCATED (altmax_ind_lastyear)) DEALLOCATE(altmax_ind_lastyear)
    IF ( ALLOCATED (O2_soil)) DEALLOCATE(O2_soil)
    IF ( ALLOCATED (CH4_soil)) DEALLOCATE(CH4_soil)
    IF ( ALLOCATED (O2_snow)) DEALLOCATE(O2_snow)
    IF ( ALLOCATED (CH4_snow)) DEALLOCATE(CH4_snow)
    IF ( ALLOCATED (tdeep_daily)) DEALLOCATE(tdeep_daily)
    IF ( ALLOCATED (fbact)) DEALLOCATE(fbact)
    IF ( ALLOCATED (decomp_rate)) DEALLOCATE(decomp_rate)
    IF ( ALLOCATED (decomp_rate_daily)) DEALLOCATE(decomp_rate_daily)
    IF ( ALLOCATED (hsdeep_daily)) DEALLOCATE(hsdeep_daily)
    IF ( ALLOCATED (temp_sol_daily)) DEALLOCATE(temp_sol_daily)
    IF ( ALLOCATED (som_input_daily)) DEALLOCATE(som_input_daily)
    IF ( ALLOCATED (pb_pa_daily)) DEALLOCATE(pb_pa_daily)
    IF ( ALLOCATED (snow_daily)) DEALLOCATE(snow_daily)
    IF ( ALLOCATED (fixed_cryoturbation_depth)) DEALLOCATE(fixed_cryoturbation_depth)
    IF ( ALLOCATED (snowdz_daily)) DEALLOCATE(snowdz_daily)
    IF ( ALLOCATED (snowrho_daily)) DEALLOCATE(snowrho_daily) 

    IF (ALLOCATED(carbon_acro)) DEALLOCATE(carbon_acro)
    IF (ALLOCATED(carbon_cato)) DEALLOCATE(carbon_cato)
    IF (ALLOCATED(height_acro)) DEALLOCATE(height_acro)
    IF (ALLOCATED(height_cato)) DEALLOCATE(height_cato)
    IF (ALLOCATED(tcarbon_acro)) DEALLOCATE(tcarbon_acro)
    IF (ALLOCATED(tcarbon_cato)) DEALLOCATE(tcarbon_cato)
    IF (ALLOCATED(resp_acro_oxic)) DEALLOCATE(resp_acro_oxic)
    IF (ALLOCATED(resp_acro_anoxic)) DEALLOCATE(resp_acro_anoxic)
    IF (ALLOCATED(resp_cato)) DEALLOCATE(resp_cato)
    IF (ALLOCATED(litter_to_acro)) DEALLOCATE(litter_to_acro)
    IF (ALLOCATED(acro_to_cato)) DEALLOCATE(acro_to_cato)
    IF (ALLOCATED(wtp_daily)) DEALLOCATE(wtp_daily)
    IF (ALLOCATED(wtp_pt)) DEALLOCATE(wtp_pt)
    IF (ALLOCATED(deepSOM_peat)) DEALLOCATE(deepSOM_peat)
    
  END SUBROUTINE stomate_clear


!! ================================================================================================================================
!! SUBROUTINE 	: stomate_var_init
!!
!>\BRIEF        Initialize variables of stomate with a none-zero initial value.
!! Subroutine is called only if ::ok_stomate = .TRUE. STOMATE diagnoses some 
!! variables for SECHIBA : assim_param, deadleaf_cover, etc. These variables can 
!! be recalculated from STOMATE's prognostic variables.
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): leaf age (::leaf_age) and fraction of leaves in leaf 
!! age class (::leaf_frac). The maximum water on vegetation available for 
!! interception, fraction of soil covered by dead leaves
!! (::deadleaf_cover) and assimilation parameters (:: assim_param).
!!
!! REFERENCE(S)	: None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  
  SUBROUTINE stomate_var_init &
       &  (kjpindex, veget_max, leaf_age, leaf_frac, &
       &   leaf_age_crit, dead_leaves, &
       &   veget, deadleaf_cover, assim_param, &
       &   circ_class_biomass, circ_class_n, sugar_load)


  !! 0. Variable and parameter declaration

    !! 0.1 Input variables
 
    INTEGER(i_std),INTENT(in)                                :: kjpindex           !! Domain size - terrestrial pixels only
    REAL(r_std),DIMENSION(:,:),INTENT(in)                    :: veget              !! Fraction of pixel covered by PFT. Fraction 
                                                                                   !! accounts for none-biological land covers 
                                                                                   !! (unitless) 
    REAL(r_std),DIMENSION(:,:),INTENT(in)                    :: veget_max          !! Fractional coverage: maximum share of the pixel 
                                                                                   !! covered by a PFT (unitless) 
    REAL(r_std),DIMENSION(:,:,:),INTENT(in)                  :: dead_leaves        !! Metabolic and structural fraction of dead leaves 
                                                                                   !! per ground area 
                                                                                   !! @tex $(gC m^{-2})$ @endtex  
    REAL(r_std),DIMENSION(:,:,:,:,:), INTENT(in)             :: circ_class_biomass !! @tex $(gC m^{-2})$ @endtex

    REAL(r_std),DIMENSION(:,:,:),INTENT(in)                  :: circ_class_n       !! Number of tree in each circumference class
    REAL(r_std),DIMENSION(:,:,:),INTENT(in)                  :: leaf_age           !! Age of different leaf classes per PFT (days)
    REAL(r_std),DIMENSION(:,:,:),INTENT(in)                  :: leaf_frac          !! Fraction of leaves in leaf age class per PFT 
                                                                                   !! (unitless; 1)
    REAL(r_std),DIMENSION(:,:),INTENT(in)                    :: sugar_load         !! Relative sugar loading of the labile pool (unitless)

    !! 0.2 Modified variables
    
    REAL(r_std),DIMENSION(:,:,:), INTENT(inout)              :: assim_param        !! min+max+opt temperatures (K) & vmax for 
                                                                                   !! photosynthesis  
                                                                                   !! @tex $(\mumol m^{-2} s^{-1})$ @endtex
    REAL(r_std),DIMENSION(:,:), INTENT(inout)                :: leaf_age_crit      !! critical leaf age (days)
   

    !! 0.3 Output variables

    REAL(r_std),DIMENSION(:), INTENT (out)                   :: deadleaf_cover     !! Fraction of soil covered by dead leaves 
                                                                                   !! (unitless) 

    ! 0.4 Local variables
   
    REAL(r_std),PARAMETER                                 :: dt_0 = zero     !! Dummy time step, must be zero
    REAL(r_std),DIMENSION(kjpindex,nvm,nleafages)         :: leaf_age_tmp    !! Temporary variable
    REAL(r_std),DIMENSION(kjpindex,nvm,nleafages)         :: leaf_frac_tmp   !! Temporary variable
                                                                             !! (unitless; 1)     
    INTEGER(i_std)                                        :: j               !! Index (untiless)
    
!_ ================================================================================================================================   

    ! Only if stomate is activated
    IF (printlev>=4) WRITE(numout,*) 'Entering stomate_var_init'
 
    !! 1. photosynthesis parameters

    !! 1. Calculate assim_param if it was not found in the restart file
    IF (ALL(assim_param(:,:,:)==val_exp)) THEN
       ! Use temporary leaf_age_tmp and leaf_frac_tmp to preserve the input variables from being modified by the subroutine vmax.
       leaf_age_tmp(:,:,:)=leaf_age(:,:,:)
       leaf_frac_tmp(:,:,:)=leaf_frac(:,:,:)

       !! 1.1 Calculate a temporary vcmax (stomate_vmax.f90)
       CALL vmax (kjpindex, dt_0, leaf_age_tmp, leaf_frac_tmp, assim_param, &
            circ_class_biomass, circ_class_n, sugar_load, leaf_age_crit, &
            leaf_classes, veget_max)
    END IF

    !! 2. Dead leaf cover (stomate_litter.f90)
    CALL deadleaf (kjpindex, veget_max, dead_leaves, deadleaf_cover)     
    
  END SUBROUTINE stomate_var_init


!! ================================================================================================================================
!! INTERFACE 	: stomate_accu
!!
!>\BRIEF        Accumulate a variable for the time period specified by 
!! dt_sechiba or calculate the mean value over the period of dt_stomate
!! 
!! DESCRIPTION : Accumulate a variable for the time period specified by 
!! dt_sechiba or calculate the mean value over the period of dt_stomate.
!! stomate_accu interface can be used for variables having 1, 2 or 3 dimensions.
!! The corresponding subruoutine stomate_accu_r1d, stomate_accu_r2d or
!! stomate_accu_r3d will be selected through the interface depending on the number of dimensions.
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): accumulated or mean variable ::field_out:: 
!!
!! REFERENCE(S)	: None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
    SUBROUTINE stomate_accu_r1d (ldmean, field_in, field_out)
    
  !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    LOGICAL,INTENT(in)                     :: ldmean    !! Flag to calculate the mean over
    REAL(r_std),DIMENSION(:),INTENT(in)    :: field_in  !! Field that needs to be accumulated
    
    !! 0.2 Modified variables
    REAL(r_std),DIMENSION(:),INTENT(inout) :: field_out !! Accumulated or mean field

!_ ================================================================================================================================

  !! 1. Accumulate field

    field_out(:) = field_out(:)+field_in(:)*dt_sechiba
   
  !! 2. Mean fields

    IF (ldmean) THEN
       field_out(:) = field_out(:)/dt_stomate
    ENDIF

  END SUBROUTINE stomate_accu_r1d

  SUBROUTINE stomate_accu_r2d (ldmean, field_in, field_out)
    
  !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    LOGICAL,INTENT(in)                       :: ldmean    !! Flag to calculate the mean over
    REAL(r_std),DIMENSION(:,:),INTENT(in)    :: field_in  !! Field that needs to be accumulated
    
    !! 0.2 Modified variables
    REAL(r_std),DIMENSION(:,:),INTENT(inout) :: field_out !! Accumulated or mean field

!_ ================================================================================================================================

  !! 1. Accumulate field

    field_out(:,:) = field_out(:,:)+field_in(:,:)*dt_sechiba
   
  !! 2. Mean fields

    IF (ldmean) THEN
       field_out(:,:) = field_out(:,:)/dt_stomate
    ENDIF

  END SUBROUTINE stomate_accu_r2d

  SUBROUTINE stomate_accu_r3d (ldmean, field_in, field_out)
    
  !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    LOGICAL,INTENT(in)                         :: ldmean    !! Flag to calculate the mean over
    REAL(r_std),DIMENSION(:,:,:),INTENT(in)    :: field_in  !! Field that needs to be accumulated
    
    !! 0.2 Modified variables
    REAL(r_std),DIMENSION(:,:,:),INTENT(inout) :: field_out !! Accumulated or mean field

!_ ================================================================================================================================

  !! 1. Accumulate field

    field_out(:,:,:) = field_out(:,:,:)+field_in(:,:,:)*dt_sechiba
   
  !! 2. Mean fields

    IF (ldmean) THEN
       field_out(:,:,:) = field_out(:,:,:)/dt_stomate
    ENDIF

  END SUBROUTINE stomate_accu_r3d


  SUBROUTINE stomate_accu_r4d (ldmean, field_in, field_out)
    
  !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    LOGICAL,INTENT(in)                         :: ldmean    !! Flag to calculate the mean over
    REAL(r_std),DIMENSION(:,:,:,:),INTENT(in)    :: field_in  !! Field that needs to be accumulated
    
    !! 0.2 Modified variables
    REAL(r_std),DIMENSION(:,:,:,:),INTENT(inout) :: field_out !! Accumulated or mean field

!_ ================================================================================================================================

  !! 1. Accumulate field

    field_out(:,:,:,:) = field_out(:,:,:,:)+field_in(:,:,:,:)*dt_sechiba
   
  !! 2. Mean fields

    IF (ldmean) THEN
       field_out(:,:,:,:) = field_out(:,:,:,:)/dt_stomate
    ENDIF

  END SUBROUTINE stomate_accu_r4d

END MODULE stomate
