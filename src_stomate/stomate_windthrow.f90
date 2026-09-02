! =================================================================================================================================
! MODULE        : stomate_wind
!
! CONTACT	: yiyingchen@gate.sinica.edu.tw
!
! LICENCE      	: IPSL (2006). This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        : Calculates critical wind speed for tree damage for each half-hourly period of the simulation.
!                 This wind speed can be compared to the actual wind speed to decide whether a tree damage occurs.
!!
!!\streamlining_n DESCRIPTION: None
!!
!! RECENT CHANGE(S): None
!!
!! SVN           :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_stomate/stomate_windthrow.f90 $ 
!! $Date: 2026-02-05 13:04:02 +0100 (jeu. 05 févr. 2026) $
!! $Revision: 9369 $
!! \streamlining_n
!_ ================================================================================================================================

MODULE stomate_windthrow

  ! modules used:
  USE ioipsl_para
  USE grid 
  USE xios_orchidee
  USE time, ONLY : dt_sechiba
  USE stomate_data
  USE constantes
  USE constantes_soil
  USE pft_parameters
  USE dynamic_parameters
  USE pft_parameters_var
  USE function_library, ONLY: wood_to_dia, wood_to_height, &
       wood_to_cv, cc_to_lai, cc_to_biomass, biomass_to_lai, &
       wood_to_cn_dia, wood_to_qmheight, get_printlev,&
       wood_to_qmdia, Nmax, wood_to_ba, cc_kill_to_area


  IMPLICIT NONE

  ! private & public routines

  PRIVATE 
  PUBLIC wind_damage
!!$  wind_clear, streamlining, calculate_force, calculate_bending_moment, &
!!$  calculate_wind_process, elevate

  LOGICAL, SAVE                   :: firstcall_windthrow = .TRUE.     !! First call
!$OMP THREADPRIVATE(firstcall_windthrow)
  LOGICAL, SAVE                   :: ok_wind_product = .FALSE.        !! Single-threshold storm trigger
                                                                      !! (wind_sum*daily_wind >= K) vs legacy AND of
                                                                      !! 2 thresholds. .FALSE. -> bit-neutral
!$OMP THREADPRIVATE(ok_wind_product)
  REAL(r_std), SAVE               :: wind_storm_product_thr = 13000.  !! Combined storm-trigger threshold K (ok_wind_product)
!$OMP THREADPRIVATE(wind_storm_product_thr)
  ! Guillaume M. -- Wind acclimation is carried by the Gompertz damage weight
  ! (wind_ratio_a/b/c, from constantes) applied to the damage rate. There is no
  ! separate acclimation knob acting on the critical wind speed.
  INTEGER(i_std), SAVE            :: printlev_loc                     !! Local level of text output for current module
!$OMP THREADPRIVATE(printlev_loc)

CONTAINS


!! ================================================================================================================================
!!  SUBROUTINE   : wind_clear 
!!
!>\BRIEF        Set the firstcall flag to .TRUE. and activate initialization
!! 
!_ ================================================================================================================================

  SUBROUTINE wind_clear
    firstcall_windthrow = .TRUE.
  END SUBROUTINE wind_clear


!! ================================================================================================================================
!! SUBROUTINE     : wind_damage
!!
!>\BRIEF          This subroutine calculates the critical wind speeds for both uprooting and stem breakage
!!                for each half-hourly period of the simulation. The algorithm follows Barry Gardiner's wind 
!!                damage risk model called GALES (Hale et al. 2015).
!! 
!! DESCRIPTION	: The purpose of this module is to describe the actual sensitivity of a given PFT to the wind that is present
!! at the same pixel. For each circumference class in each PFT, two critical wind speed values (CWS) are calculated in every
!! half-hourly time step by this subroutine, following the approach of differentiating two types of wind damage of trees:
!! - uprooting (the whole tree is removed from the soil together with its root plate)
!! - stem breakage (the root plate and the bottom part of the stem remains attached to the ground).
!! Whichever threshold (i.e. CWS) is reached first for a tree, the according damage type occurs in the model.
!!
!! The module completes the following tasks:
!! 1. Calculates...
!! 
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): :: Critical wind speed for stem breakage 10 m above the zero-plane displacement (m/s): U_Break_10
!!                             Critical wind speed for tree overturn 10 m above the zero-plane displacement (m/s): U_Overturn_10
!!
!! REFERENCES    :
!! - Hale, S., Gardiner, B., Nicoll, B., Taylor, P., and Pizzirani, S. (2015). Comparison and validation of three versions
!! of a forest wind risk model. Environmental Modelling and Software. In Press.
!! ...
!!
!! FLOWCHART     : None
!! \streamlining_n 
!_ ================================================================================================================================
  
  SUBROUTINE wind_damage (npts, nlevels_tot, circ_class_biomass, &
       veget_max, circ_class_n, plant_status, &
       circ_class_kill, gap_area_save, forest_managed, &
       soil_temp_daily, root_profile, max_wind_speed_storm, max_wind_ratio_storm, &
       count_storm, is_storm, wind_ratio_max_save, wind_ratio_sum_save, &
       wind_speed_max_save)


  IMPLICIT NONE

 !! 0. Variable declarations

 !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                       :: npts                     !! Domain size @tex $(unitless)$ @endtex
    INTEGER(i_std), INTENT(in)                       :: nlevels_tot              !! Total level numbers in photosynthesis
    INTEGER(i_std), DIMENSION(:,:), INTENT(in)       :: forest_managed           !! forest management flag (is the forest
                                                                                 !! being managed?)
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(in)    :: circ_class_biomass       !! Biomass in of an individual tree in each circ class
                                                                                 !! @tex $(m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)        :: circ_class_n             !! Number of trees within each circumference
    REAL(r_std), DIMENSION(:,:), INTENT(in)          :: veget_max                !! "maximal" coverage fraction of a PFT on the ground
                                                                                 !! May sum to less than unity if the pixel has nobio 
                                                                                 !! area. (unitless, 0-1)
    REAL(r_std), DIMENSION(:,:), INTENT(in)          :: plant_status             !! Growth and phenological status of the plant
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)        :: gap_area_save            !! Gap area created by more than 30% of the basal area loss
                                                                                 !! in the last 5 years
                                                                                 !! Dimension(npts,nvm,legacy_years_wind)
    REAL(r_std), DIMENSION(:), INTENT(in)            :: soil_temp_daily          !! Daily maximum soil temperature at 0.8 meter below ground (K)
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)      :: root_profile             !! Normalized root mass/length fraction in each soil layer (0-1, unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(in)          :: wind_ratio_sum_save      !! Sum of wind speed ratio to longterm wind speed when this ratio exceeded
                                                                                 !! certain threshold, stored for wind_days
    REAL(r_std), DIMENSION(:,:), INTENT(in)          :: wind_ratio_max_save      !! Daily maximum wind ratio (vs 5yr mean), stored for wind_days - JJ 2026
    REAL(r_std), DIMENSION(:,:), INTENT(in)          :: wind_speed_max_save      !! Daily maximum wind speed, stored for wind_days (default 3) (m s-1)

 !! 0.2 Output variables
    
 !! 0.3 Modified variables

    REAL(r_std), DIMENSION(:,:,:,:,:), &
                                   INTENT(inout) :: circ_class_kill              !! Number of trees within a circ that needs
                                                                                 !! to be killed @tex $(ind m^{-2})$ @endtex  
    REAL(r_std), DIMENSION(:), INTENT(inout)     :: max_wind_speed_storm         !! Daily maximum wind speed at 2 meter (ms-1) during a storm event
    REAL(r_std), DIMENSION(:), INTENT(inout)     :: max_wind_ratio_storm         !! Daily maximum wind ratio (vs 5yr mean) during a storm event - JJ 2026

    INTEGER(i_std), DIMENSION(:), INTENT(inout)  :: count_storm                  !! Number of days after a storm

    LOGICAL, DIMENSION(:), INTENT(inout)         :: is_storm                     !! Are we in a storm event ? 

 !! 0.4 Local variables
    INTEGER(i_std)                               :: islm,ipts,ivm                !! Index for the loops
    INTEGER(i_std)                               :: icir,i,iele, iyear           !! Index for the loops
    REAL(r_std), DIMENSION(npts,nvm)             :: area_damaged_storm           !! Forest area killed by storm, cc-kill basal-area route (m2)
    CHARACTER(30)                                :: var_name                     !! To store variable names for I/O
    INTEGER(i_std), DIMENSION(npts)              :: soil_type                    !! Soil types: 
                                                                                 !!  1. Free-draining mineral soils;
                                                                                 !!  2. Gleyed mineral soils; 
                                                                                 !!  3. Peaty mineral soils; and 
                                                                                 !!  4. Deep Peats
    REAL(r_std), DIMENSION(npts)                 :: wind_speed_daily             !! Daily wind speed to be used to calculated wind damage (m s-1)
    REAL(r_std), DIMENSION(npts,nvm,ncirc)       :: virtual_circ_class_n         !! Virtual number of trees within each circumference
    REAL(r_std), DIMENSION(ncirc)                :: mean_dbh                     !! Diameter at breast height (dbh) of the mean tree
                                                                                 !! of the circumference-class @tex $(m)$ @endtex
    REAL(r_std), DIMENSION(ncirc)                :: ccn_inv                      !! JJ 2026: circ_class_n restricted to inventory-visible
                                                                                 !! trees (dbh >= dia_thresh_inv), used for spacing/RDI (ind m-2)
    REAL(r_std), DIMENSION(ncirc)                :: mean_height                  !! Height of the mean tree of the given dbh-class 
                                                                                 !! @tex $(m)$ @endtex
    REAL(r_std)                                  :: lai                          !! Leaf area index @tex $(unitless)$ @endtex
    REAL(r_std)                                  :: canopy_breadth_dom           !! Canopy breadth (width of canopy) @tex $(m)$ @endtex
    REAL(r_std)                                  :: canopy_depth_dom             !! Canopy depth (height of canopy) @tex $(m)$ @endtex
    REAL(r_std)                                  :: canopy_height_dom            !! The height of the starting point of the canopy from 
                                                                                 !! below @tex $(m)$ @endtex
    REAL(r_std), DIMENSION(ncirc)                :: crown_width                  !! Crown width of the mean tree (m)  
    REAL(r_std), DIMENSION(ncirc)                :: crown_depth                  !! Crown width of the mean tree (m)  
    REAL(r_std), DIMENSION(ncirc)                :: crown_volume                 !! Crown volume of the mean tree @tex $(m^{3})$ @endtex
    REAL(r_std), DIMENSION(ncirc)                :: crown_mass                   !! Mass of the mean tree crown @tex $(kg)$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)             :: stem_mass                    !! Mass of the mean tree stem @tex $(kg)$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)             :: current_spacing              !! Spacing between trees in the given dbh-class 
                                                                                 !! @tex $(m)$ @endtexi
    INTEGER(i_std)                               :: rooting_depth                !! Rooting depth of the mean tree (three categories;
                                                                                 !! 1: shallow, 2: deep, 3:average 
                                                                                 !! See constantes_soil_var.f90 for the definitions
    REAL(r_std)                                  :: streamlining_c               !! Streamlining parameter. @tex $(unitless)$ @endtex. 
                                                                                 !! Streamlining is the change of shape of the crowns 
                                                                                 !! due to wind.
    REAL(r_std)                                  :: streamlining_n               !! Streamlining parameter. @tex $(unitless)$ @endtex. 
                                                                                 !! Streamlining is the change of shape of the crowns 
                                                                                 !! due to wind.
    REAL(r_std)                                  :: area_timber_removals_5_years !! Area of the total timber removal in the
                                                                                 !! previous 5 year @tex $(m^{2})$ @endtex
    REAL(r_std)                                  :: tree_heights_from_edge       !! Width of the forest edge area, i.e. area_closer (in number of tree heights)
                                                                                 !! @tex $(unitless)$ @endtex
    REAL(r_std)                                  :: n_gap                        !! The number of gaps with an area of clear_cut_max @tex $(unitless)$ @endtex
    REAL(r_std)                                  :: length_side_gap              !! The length of one side of the square-shaped gap @tex $(m)$ @endtex
    REAL(r_std)                                  :: area_around_gap              !! The framing edge area around the gap with different gustiness
                                                                                 !! @tex $(m^{2})$ @endtex
    REAL(r_std)                                  :: area_total_closer            !! The total area of forest edge areas affected by wind @tex $(m^{2})$ @endtex
    REAL(r_std)                                  :: area_total_further           !! The total area of forest without edge areas @tex $(m^{2})$ @endtex
    REAL(r_std)                                  :: ratio_spacing_height         !! Ratio of tree height and tree spacing @tex $(unitless)$ @endtex
    REAL(r_std)                                  :: gap_size                     !! The length (width of the cross-section) of the gap @tex $(m)$ @endtex
    REAL(r_std)                                  :: mean_gap_factor              !! Variable used for calculating gustiness 
                                                                                 !! @tex $(unitless)$ @endtex
    REAL(r_std)                                  :: max_gap_factor               !! Variable used for calculating gustiness 
                                                                                 !! @tex $(unitless)$ @endtex
    REAL(r_std)                                  :: term_a                       !! Variable used for calculating gustiness 
                                                                                 !! @tex $(unitless)$ @endtex
    REAL(r_std)                                  :: term_b                       !! Variable used for calculating gustiness 
                                                                                 !! @tex $(unitless)$ @endtex
    REAL(r_std)                                  :: new_gust_closer              !! Variable used for calculating gustiness 
                                                                                 !! @tex $(unitless)$ @endtex
    REAL(r_std)                                  :: new_gust_further             !! Variable used for calculating gustiness 
                                                                                 !! @tex $(unitless)$ @endtex
    REAL(r_std)                                  :: new_gust_edge                !! Variable used for calculating gustiness 
                                                                                 !! @tex $(unitless)$ @endtex
    REAL(r_std)                                  :: mean_bm_gust_factor_closer   !! Mean bending moment of the forest edge area 
                                                                                 !! (for calculating the edge factor)
                                                                                 !! @tex $(Nm/kg)$ @endtex
    REAL(r_std)                                  :: mean_bm_gust_factor_further  !! Mean bending moment of the inner forest area (for calculating
                                                                                 !! the edge factor) @tex $(Nm/kg)$ @endtex
    REAL(r_std)                                  :: edge_factor                  !! The ratio of the mean bending moment in the forest edge area 
                                                                                 !! and mean bending moment in the inner forest area
                                                                                 !! @tex $(unitless)$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)             :: gust_factor_closer           !! Gust factor for calculating the critical wind speed of the forest
                                                                                 !! edge area @tex $(unitless)$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)             :: gust_factor_further          !! Gust factor for calculating the critical wind speed of the inner
                                                                                 !! forest area @tex $(unitless)$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)             :: g_closer                     !! Gustiness of the forest edge area @tex $(unitless)$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)             :: g_further                    !! Gustiness of the inner forest area @tex $(unitless)$ @endtex
    REAL(r_std)                                  :: overturning_moment_multiplier!! Overturning moment multiplier representing the effects of species,
                                                                                 !! soil type, soil depths and senescence @tex $(Nm/kg)$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)             :: max_overturning_moment       !! The maximum overturning moment the mean tree can withstand
                                                                                 !! @tex $(Nm/kg)$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)             :: max_breaking_moment          !! The maximum breaking moment the mean tree can withstand
                                                                                 !! @tex $(Nm/kg)$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)             :: overturning_moment_closer    !! The maximum overturning moment the mean tree can withstand in the
                                                                                 !! forest edge area with gustiness taken into account @tex $(Nm/kg)$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)             :: overturning_moment_further   !! The maximum overturning moment the mean tree can withstand in the
                                                                                 !! inner forest area with gustiness taken into account @tex $(Nm/kg)$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)             :: breaking_moment_closer       !! The maximum breaking moment the mean tree can withstand in the
                                                                                 !! forest edge area with gustiness taken into account @tex $(Nm/kg)$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)             :: breaking_moment_further      !! The maximum breaking moment the mean tree can withstand in the
                                                                                 !! inner forest area with gustiness taken into account @tex $(Nm/kg)$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)             :: cws_break_closer             !! Critical wind speed for stem breakage in the forest edge area
                                                                                 !! @tex $(m/s; mean hourly wind speed)$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)             :: cws_break_closer_10          !! Critical wind speed for stem breakage in the forest edge area
                                                                                 !! @tex $(m/s; mean hourly wind speed)$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)             :: cws_overturn_closer          !! Critical wind speed for overturning in the forest edge area
                                                                                 !! @tex $(m/s; mean hourly wind speed)$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)             :: cws_overturn_closer_10       !! Critical wind speed for overturning in the forest edge area
                                                                                 !! @tex $(m/s; mean hourly wind speed)$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)             :: cws_final_closer             !! The critical wind speed of the final form of damage (breakage / uprooting)
                                                                                 !! in the forest edge area @tex $(m/s; mean hourly wind speed)$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)             :: d_h_ratio                    !! Array for storing the stand spacing and tree height ratio for output (-)
    REAL(r_std), DIMENSION(npts)                 :: wind_speed_actual            !! The actual wind speed over the pixel
                                                                                 !! @tex $(m/s; mean hourly wind speed)$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)             :: wind_damage_rate_closer      !! Damage rate (%) in the forest edge area in case of damage
                                                                                 !! @tex $(0-1)$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)             :: cws_break_further            !! Critical wind speed for stem breakage in the inner forest area
                                                                                 !! @tex $(m/s; mean hourly wind speed)$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)             :: cws_break_further_10         !! Critical wind speed for stem breakage in the inner forest area
                                                                                 !! @tex $(m/s; mean hourly wind speed)$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)             :: cws_overturn_further         !! Critical wind speed for overturning in the inner forest area
                                                                                 !! @tex $(m/s; mean hourly wind speed)$ @endtex
!Jina: Seems it if for 10 m above the z0 but has the same description.
    REAL(r_std), DIMENSION(npts,nvm)             :: cws_overturn_further_10      !! Critical wind speed for overturning in the inner forest area
                                                                                 !! @tex $(m/s; mean hourly wind speed)$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)             :: cws_final_further            !! The critical wind speed of the final form of damage (breakage / uprooting)
                                                                                 !! in the inner forest area @tex $(m/s; mean hourly wind speed)$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)             :: wind_damage_rate_further     !! Damage rate (%) in the forest edge area in case
                                                                                 !! of damage in the inner forest area @tex $(0-1)$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)             :: z0_wind                      !! surface roughness calculated from the subroutine streamlining()
    REAL(r_std), DIMENSION(npts,nvm)             :: zhd_wind                     !! zero plane displacement height calculated from the subroutine streamlining(m) 
    REAL(r_std), DIMENSION(npts,nvm)             :: gamma_solved_bk_closer
    REAL(r_std), DIMENSION(npts,nvm)             :: gamma_solved_ov_closer
    REAL(r_std), DIMENSION(npts,nvm)             :: gamma_solved_bk_further
    REAL(r_std), DIMENSION(npts,nvm)             :: gamma_solved_ov_further      !! the ratio between the critical wind speed at canopy top and the friction
                                                                                 !! velocity (uh / u*)

    INTEGER(i_std), DIMENSION(npts,nvm)          :: wind_damage_type_closer      !! Integer store the wind damage type for closer area  ibreakage=1, ioverturning=2
    INTEGER(i_std), DIMENSION(npts,nvm)          :: wind_damage_type_further     !! Integer store the wind damage type for further area
    REAL(r_std), DIMENSION(npts,nvm,ncirc)       :: kill_break                   !! Biomass loss through stem breakage from storm 
                                                                                 !! @tex $(ind m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,ncirc)       :: kill_uproot                  !! Number of trees killed through uproot from storm
                                                                                 !! @tex $(ind m^{-2})$ @endtex
    REAL(r_std)                                  :: crown_top                    !! Temporal variable for storing the tallest tree crown height in different PFTs  
    REAL(r_std)                                  :: crown_bottom                 !! Temporal variable for storing the smallest tree crown height in different PFTs 
    REAL(r_std)                                  :: distance_from_force          !! The distance from the force was set to 1.3 for breakage and 0.0 for overturning
    REAL(r_std), DIMENSION(npts)                 :: closer_area_fraction         !! Closer area fraction 
    REAL(r_std), DIMENSION(npts)                 :: further_area_fraction        !! Further area fraction 
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements) :: biomass  
    REAL(r_std), DIMENSION(npts)                 :: wind_in                      !! To save the daily maximum wind speed at 2 meter (ms-1)
    REAL(r_std), DIMENSION(npts,nvm)             :: DLF_breaking_closer
    REAL(r_std), DIMENSION(npts,nvm)             :: DLF_overturning_closer
    REAL(r_std), DIMENSION(npts,nvm)             :: DLF_breaking_further
    REAL(r_std), DIMENSION(npts,nvm)             :: DLF_overturning_further
    REAL(r_std)                                  :: eta                          !! Ratio of friction velocity to the wind speed at the canopy top. Also see condveg.f90
    REAL(r_std), DIMENSION(npts)                 :: zhd_pixel                    !! pixel level zpd
    REAL(r_std)                                  :: z0m_pft                      !! Roughness height for momentum for a specific PFT
    REAL(r_std)                                  :: pft_height_dom
    REAL(r_std), DIMENSION(npts,nvm)             :: max_wind_speed_pft_further
    REAL(r_std), DIMENSION(npts,nvm)             :: max_wind_speed_pft_closer
    REAL(r_std), DIMENSION(npts,nvm)             :: wind_damage_all              !! Damage from wind throw (gC m-2) 
    REAL(r_std), DIMENSION(npts,nvm)             :: wind_damage_rate             !! Damage rate of biomass from wind throw (gC m-2) 
    REAL(r_std), DIMENSION(npts,nvm)             :: wood_damage_all              !! Damage from wind throw (gC m-2) 

    REAL(r_std), DIMENSION(npts,nvm)             :: height
    REAL(r_std)                                  :: ave_height
    REAL(r_std)                                  :: z0_wind_bk
    REAL(r_std)                                  :: z0_wind_ov
    REAL(r_std)                                  :: zhd_wind_bk
    REAL(r_std)                                  :: zhd_wind_ov
    REAL(r_std), DIMENSION(npts,nvm)             :: qm_dia
    REAL(r_std), DIMENSION(npts,nvm)             :: rdi
    REAL(r_std)                                  :: actual_maxheight
    REAL(r_std), DIMENSION(npts)                 :: weight_ratio                 !! JJ 2026: Gompertz damage weight from storm wind-ratio (0..wind_ratio_a)
    REAL(r_std), DIMENSION(npts)                 :: wind_ratio_max               !! JJ 2026: yesterday's daily max wind ratio (diagnostic)

    ! Temp output Variables for GALES model intercomparision
    REAL(r_std), DIMENSION(npts,nvm)             :: canopy_depth_out
    REAL(r_std), DIMENSION(npts,nvm)             :: tree_h_edge_out
    REAL(r_std), DIMENSION(npts,nvm)             :: mean_gap_f_out
    REAL(r_std), DIMENSION(npts,nvm)             :: gap_size_out
    REAL(r_std), DIMENSION(npts,nvm)             :: edge_factor_out
    REAL(r_std)                                  :: roots_below                  !! temporary variable containing the share of roots below 80 cm depth (0-1, unitless)
    REAL(r_std), DIMENSION(npts)                 :: is_kill
    REAL(r_std), DIMENSION(npts,nvm)             :: wind_damage_wood_vol         !! Wood volume loss from wind throw @tex $(m^3 m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)             :: wind_profile_closer          !! Effective wind multiplier (log-profile, or 1) for the closer logistic
    REAL(r_std), DIMENSION(npts,nvm)             :: wind_profile_further         !! Effective wind multiplier (log-profile, or 1) for the further logistic

!_ ============================================================================================================================================================================

    IF (firstcall_windthrow) THEN
       !! Initialize local printlev
      ! Guillaume M. -- get_printlev calls getin_p, a COLLECTIVE bcast. Called here,
      ! inside the compute loop under a local firstcall, it hangs the run as soon as the
      ! ranks stop reaching it together. Use printlev, read once at startup by every rank.
       printlev_loc=printlev
       !Config Key   = OK_WIND_PRODUCT
       !Config Desc  = Single-threshold storm trigger combining wind magnitude and anomaly
       !Config If    = OK_WINDTHROW
       !Config Def   = n
       !Config Help  = single-threshold storm trigger wind_sum*daily_wind >= WIND_STORM_PRODUCT_THR
       !Config         (combines absolute magnitude + relative/acclimation anomaly).
       !Config         Bit-neutral when .FALSE.
       !Config         [CATEGORIE : OPTION SCIENTIFIQUE -- choix de modelisation, permanent]
       !Config Units = [FLAG]
       ok_wind_product = .FALSE.
       CALL getin_p('OK_WIND_PRODUCT', ok_wind_product)
       !Config Key   = WIND_STORM_PRODUCT_THR
       !Config Desc  = Threshold of the single-product storm trigger
       !Config If    = OK_WIND_PRODUCT
       !Config Cat   = SCIENTIFIC OPTION
       !Config Def   = 13000.
       !Config Help  = The trigger fires when wind_sum * daily_wind reaches this value, so one
       !Config         number carries both the absolute magnitude and the relative anomaly.
       !Config         Only read when OK_WIND_PRODUCT is true; the storm module is bit-neutral
       !Config         otherwise.
       !Config Units = [m2/s2]
       wind_storm_product_thr = 13000.
       CALL getin_p('WIND_STORM_PRODUCT_THR', wind_storm_product_thr)

       ! The code is specific for CRU-NCEP. If a different forcing
       ! is used, the Max Wind Ratio parameters should be adjusted
    END IF
    
    IF (printlev_loc>=2) WRITE(numout,*) 'Entering stomate_windthrow.f90'
    
    !! 0. Initialized the values of variables
    wind_damage_rate_closer = zero
    wind_damage_rate_further = zero
    wind_damage_rate = zero
    cws_final_further = zero
    cws_final_closer = zero
    d_h_ratio = zero
    kill_uproot = zero
    kill_break = zero
    cws_break_further_10 = zero
    cws_break_closer_10 = zero
    cws_overturn_further_10 = zero
    cws_overturn_closer_10 = zero
    cws_break_further = zero
    cws_break_closer = zero
    cws_overturn_further = zero
    cws_overturn_closer = zero
    virtual_circ_class_n = zero
    gamma_solved_bk_closer = zero
    gamma_solved_ov_closer = zero
    gamma_solved_bk_further = zero
    gamma_solved_ov_further = zero
    current_spacing = zero
    zhd_wind = zero
    breaking_moment_closer = zero
    breaking_moment_further = zero
    overturning_moment_closer = zero
    overturning_moment_further = zero
    is_kill = un
    z0_wind = zero
    wind_damage_all = zero
    wood_damage_all = zero
    wind_damage_wood_vol = zero
    ! Guillaume M. -- Storm emulator: stand state independent of the flags, reset at each call.
    wind_profile_closer = un
    wind_profile_further = un
    DLF_breaking_closer = zero
    DLF_overturning_closer = zero
    DLF_breaking_further = zero
    DLF_overturning_further = zero
    max_wind_speed_pft_closer = zero
    max_wind_speed_pft_further = zero
    height = zero
    rdi = zero
    qm_dia = zero
    closer_area_fraction = zero
    further_area_fraction = zero

    ! Guillaume M. -- wind_in is set per pixel from wind_speed_max_save in the storm loop below.
    !-Debug-
    IF (printlev_loc .GT. 4) THEN
        DO ipts = 1, npts
            WRITE(numout,*) 'Daily maximum windspeed passed from &
                stomate to stomate_winthrow:', wind_speed_daily(ipts)
            WRITE(numout,*) 'Daily maximum soil temperature at 80cm below ground:',&
                soil_temp_daily(ipts)
        ENDDO
    ENDIF
    !-
    
    !! 1. Getting the input variables in the required form for the actual subroutine    
    ! calculate the total biomass     
    biomass(:,:,:,:) = cc_to_biomass(npts,nvm,circ_class_biomass,circ_class_n)

    !! 1.1 Convert 6-hr CRU-NCEP into half-hourly wind speed     
    ! The max_wind_ratio parameters were determined by relating fluxnet
    ! half-hourly data to 6h CRU-NCEP data. For this reason the
    ! parameters should be considered specific to the CRU-NCEP forcing.
    ! The max_wind_ratio parameters were also found to depend on the 
    ! spatial aggregation, hence, they are specific to the 0.5 degree
    ! CRU-NCEP data (see Chen et al for details). 
   
    ! The following relationship to convert CRU-NCEP into half-hourly
    ! wind speeds was found:
    !wind_speed_daily(:) = -5.922*wind_speed_daily(:)  &
    !        + 2.051*wind_speed_daily(:)**2.  &
    !        - 0.191*wind_speed_daily(:)**3. &
    !        + 0.006*wind_speed_daily(:)**4.
    !wind_speed_daily(:) = 1.9751 + 0.4783 * wind_speed_daily(:)  !raw
    ! Based on the comparison between CRUJRA and 1) the multiple FLUXNET sites 
    ! (See Chen 2018) and 2) maximum wind speed from during the strom recorded 
    ! in XMS datasets, it is decided to use CRUJRA wind without manipulation. 

    ! In case this convertion did not result in the observed wind
    ! damage, we kept a linear tuning parameter. Ideally, no tuning
    ! should be applied and therefore daily_max_tune should be 1.0
    ! wind_speed_daily(:) = wind_max_daily(:) * daily_max_tune 
    ! Now maximum wind speed over wind_days is used to calculated wind damage
    IF(use_fluxnet) THEN
        wind_speed_daily(:) = MAXVAL(wind_speed_max_save(:,:),DIM=2)
    ELSE
        wind_speed_daily(:) = MAXVAL(wind_speed_max_save(:,:),DIM=2)*daily_max_tune
    ENDIF

    ! Storms do not last for exactly 24 hours. 
    ! Storms that pass overnight would be accounted for twice 
    ! unless we add criteria to tell ORCHIDEE this is just a
    ! single storm. The value wind_speed_damage_thr
    ! corresponds to the wind speed threshold above which is_storm
    ! flag is set the TRUE and nb_days_storm (set to 5.0 days) correspond to the 
    ! number of days at which the max wind speed is less than the
    ! threshold which means that the storm is over and hence 
    ! is_storm flag to FALSE. As a consequence the windstorm 
    ! damage is calculated nb_days_storm (set to 5 days) 
    ! after the end of a storm. One main consequences of this
    ! new approach is that damages due to wind is only accounted
    ! when the flag is_storm is set to TRUE.
    IF (ok_wind_product) THEN
        ! Guillaume M. -- Single-threshold onset: the accumulated relative-wind anomaly
        ! (wind_ratio_sum, acclimation to the local long-term wind) times the absolute daily
        ! wind, tested against one combined threshold. Magnitude and local anomaly are thus
        ! carried by a single metric. See DESIGN_CALIB_PERTURBATION.md.
        WHERE (SUM(wind_ratio_sum_save(:,:),2) * wind_speed_daily(:) .GE. wind_storm_product_thr &
               .AND. .NOT. is_storm(:))
            is_storm(:) = .TRUE.
            count_storm(:) = 0
        ENDWHERE
    ELSE
        WHERE (SUM(wind_ratio_sum_save(:,:),2) .GE. wind_sum_threshold &
               .AND. wind_speed_daily(:) .GE. wind_speed_storm_thr &
               .AND. .NOT. is_storm(:))
            is_storm(:) = .TRUE.
            count_storm(:) = 0
        ENDWHERE
    ENDIF

    IF(printlev_loc.ge.4) THEN
        WRITE(numout,*) 'test grid',test_grid
        WRITE(numout,*) 'wind sum',SUM(wind_ratio_sum_save(test_grid,:))
        WRITE(numout,*) 'threshold', wind_sum_threshold
        WRITE(numout,*) 'is_storm bf set',is_storm(test_grid)
        WRITE(numout,*) 'count_storm bf set',count_storm(test_grid)
        WRITE(numout,*) 'max_wind',wind_speed_daily(test_grid)
    ENDIF

    DO ipts = 1, npts

       ! Guillaume M. -- Per-pixel daily maximum wind and wind ratio, used by the
       ! diagnostics and by the Gompertz damage weight below.
       wind_in(ipts) = wind_speed_max_save(ipts,1)
       wind_ratio_max(ipts) = wind_ratio_max_save(ipts,1)

       max_wind_speed_storm(ipts)= &
         max(max_wind_speed_storm(ipts), wind_speed_daily(ipts))
       max_wind_ratio_storm(ipts) = &
         max(max_wind_ratio_storm(ipts), MAXVAL(wind_ratio_max_save(ipts,:)))

       ! Guillaume M. -- Gompertz weight on the storm wind ratio (against the 5-year mean):
       ! scales the damage magnitude with how anomalous the storm wind is locally.
       weight_ratio(ipts) = &
           wind_ratio_a * &
           exp(-exp(-wind_ratio_b * (max_wind_ratio_storm(ipts) - wind_ratio_c)))

       IF(is_storm(ipts)) THEN
!          max_wind_speed_storm(ipts)= &
!             max(max_wind_speed_storm(ipts), &
!             wind_speed_daily(ipts))

          ! Guillaume M. -- The storm decays (count++) once the trigger metric falls back
          ! below its threshold, using the same metric as the onset test above.
          IF ( ( ok_wind_product .AND. &
                 SUM(wind_ratio_sum_save(ipts,:))*wind_speed_daily(ipts) .LT. wind_storm_product_thr ) &
               .OR. ( (.NOT.ok_wind_product) .AND. &
                 SUM(wind_ratio_sum_save(ipts,:)) .LT. wind_sum_threshold ) ) THEN
             count_storm(ipts) =  count_storm(ipts) + 1
          ENDIF
          IF(count_storm(ipts) .LT. nb_days_storm) THEN
!             CYCLE
              is_kill(ipts) = zero
          ELSE
             count_storm(ipts) = 0
             is_storm(ipts) = .FALSE.
          ENDIF
       ELSE
!         CYCLE
          is_kill(ipts) = zero
       ENDIF

      IF(printlev_loc.ge.4 .AND. ipts .EQ. test_grid) THEN
          WRITE(numout,*) 'wind throw, ipts',ipts
          WRITE(numout,*) 'is_storm af set',is_storm(ipts)
          WRITE(numout,*) 'is_kill',is_kill(ipts)
          WRITE(numout,*) 'count_storm af set',count_storm(ipts)
          WRITE(numout,*) 'max_wind',wind_speed_daily(ipts)
          WRITE(numout,*) 'wind_mod',wind_speed_daily(ipts)
          WRITE(numout,*) 'wind_ratio_sum', SUM(wind_ratio_sum_save(ipts,:))
      ENDIF

      !! 1.3 Defining the soil type
      ! This sub-routine uses 4 types of soil for calculating the 
      ! critical wind speed, so 12 different soil types present in 
      ! other modules of ORCHIDEE need to be grouped into these
      ! 4 generic types.

      ! +++CHECK+++
      ! For the moment only 3 soil types are distinguished.
      ! Needs to be checked with all 12 USDA soil types. 
      ! If more soil type are used it will become possible
      ! to make use of the different soils distinguished in
      ! the windthrow module. The code could then look as
      ! follow:
!!$       SELECT CASE()
!!$       CASE (x:y)
!!$          soil_type(ipts) = ifreedrainage  
!!$       CASE (x:y)
!!$          soil_type(ipts) = igleyed
!!$       CASE (x:y)
!!$          soil_type(ipts) = ipeaty
!!$       CASE (x:y)
!!$          soil_type(ipts) = ipeat 
!!$       CASE DEFAULT
!!$          CALL ipslerr_p(3, 'stomate_windthrow.f90',&
!!$               'soil type is not defined','','')
!!$       END SELECT

      ! For the moment all soils are assumed to be freely draining
      soil_type(ipts) = ifree_draining
      ! +++++++++++++

      ave_height = zero
      DO ivm = 2,nvm
         height(ipts,ivm) = wood_to_qmheight(circ_class_biomass(ipts,ivm,:,:,icarbon), &
                     circ_class_n(ipts,ivm,:), ivm, pipe_tune2(ipts,ivm))

         lai = cc_to_lai(circ_class_biomass(ipts,ivm,:,ileaf,icarbon),&
              circ_class_n(ipts,ivm,:),ivm)
         IF (lai.GT.min_sechiba .AND. height(ipts,ivm).LT.min_sechiba) THEN
            height(ipts,ivm) = 0.1
         ENDIF
         ave_height =ave_height + veget_max(ipts,ivm)*height(ipts,ivm)
      ENDDO

      zhd_pixel(ipts) = ave_height *  height_displacement

      DO ivm = 2, nvm   
        ! Initialize
        area_timber_removals_5_years = zero
        area_timber_removals_5_years = zero
        n_gap = zero

        ! Don't do the calculations if there is no vegetation
        ! or the vegetation is bare soil, grass or crop
        IF ((veget_max(ipts,ivm).LT.min_stomate) .OR. .NOT.is_tree(ivm)) THEN
            CYCLE
        ENDIF
       
        !! 1.4 Calculating the rooting depth
        ! If root density (proportion of root mass) is higher than 40% below 80 cm below
        ! the soil surface, then we treat rooting depth as "deep". Otherwise, it's
        ! "shallow". Note that when using the static root profile deep vs shallow will be
        ! the same for all pixels belonhing to the same PFT because it solely depends on 
        ! humcste. The threshold of 40% was chosen to ensure that all trees are shallow
        ! rooting with the static root profile.
        IF (grnd_80.EQ.nslm) THEN
            ! Handle the case where the deepest soil layer is less
            ! than 80 cm. All roots must then be above 80 cm.
            roots_below = un
        ELSE
            ! Calculate the share of roots below 80 cm
            roots_below = zero
            DO islm = grnd_80,nslm
                roots_below = roots_below + root_profile(ipts,ivm,islm,istruc)
            ENDDO
        ENDIF
        IF (roots_below > 0.40) THEN
            rooting_depth = ideep
        ELSE
            rooting_depth = ishallow
        END IF
         
        !! 1.5 Set crown parameters
        !  To set crown parameters leaf mass needs to be available
        IF (SUM(circ_class_biomass(ipts,ivm,:,ileaf,icarbon)) .GT. min_stomate) THEN
            ! There is a canopy
            streamlining_c = streamlining_c_leaf(ivm)
            streamlining_n = streamlining_n_leaf(ivm)
        ELSE
            ! If the tree is leafless, the following parameters become different. 
            ! The weight of leaves will be absent.
            ! Streamlining changes as the surface of leaves is absent.
            streamlining_c = streamlining_c_leafless(ivm)
            streamlining_n = streamlining_n_leafless(ivm)
        END IF

        !-Debug-
        IF(printlev_loc .GT. 4 .AND. ivm .EQ. test_pft .AND. ipts .EQ. test_grid)THEN
            WRITE(numout,*) 'Initialisation of critical windspeed'
            WRITE(numout,*) 'ipts, ivm, ',ipts, ivm
            WRITE(numout,*) 'soil_type(ipts), ',soil_type(ipts)
            WRITE(numout,*) 'rooting_depth_type,(1:deep, 2:shallow,3:average), ',&
                rooting_depth
            WRITE(numout,*) 'streamlining_c, ', streamlining_c
            WRITE(numout,*) 'streamlining_n, ', streamlining_n
        ENDIF

        !! 1.6 Calculate stand characteristics
        !  Calculating diameter at breast height (dbh) in m of the
        !  mean tree of all circ classes.
        mean_dbh(:) = wood_to_dia(circ_class_biomass(ipts,ivm,:,:,icarbon),&
             ivm,pipe_tune2(ipts,ivm))

        ! Calculating the height in m of the mean tree of the given dbh-classes
        mean_height(:) = wood_to_height(circ_class_biomass(ipts,ivm,:,:,icarbon),&
             ivm,pipe_tune2(ipts,ivm))
          
        ! Crown_volume is calculated as the woody biomass of an individual tree 
        crown_volume(:) = wood_to_cv( circ_class_biomass(ipts,ivm,:,:,icarbon), &
              circ_class_n(ipts,ivm,:),ivm,pipe_tune2(ipts,ivm))

        ! Crown_width is calculated as(crown_vloume * 6.0/pi)^(1/3)
        !crown_width(:) =  ((crown_volume(:)*6.0)/pi)**(1.0/3.0)      
        CALL wood_to_cn_dia(circ_class_biomass(ipts,ivm,:,:,icarbon), &
               circ_class_n(ipts,ivm,:),ivm,crown_width, &
               crown_depth,pipe_tune2(ipts,ivm))

        ! ORCHIDEE assumes that the crowns are spherical. Hence, width breath and depth
        ! are all the same. Set canopy_breadth as the crown_width for all circumference 
        ! classes. Not anymore! 
        ! After rXXXX, only the biggest class is used for calculating critical wind speed         
        canopy_breadth_dom = crown_width(ncirc) 
        canopy_depth_dom = crown_depth(ncirc)  
        canopy_height_dom = mean_height(ncirc)

        ! Calculate gaps and the surrounding areas 
        ! IMPORTANT: AS WE DO NOT HAVE INFORMATION ON STAND EDGES IN THE ORCHIDEE 
        ! SIMULATIONS, A POSSIBLE PROXY TO USE HERE IS A VARIABLE THAT REPRESENTS EARLIER 
        ! GAP AREAS IN THE PREVIOUS 5 YEARS. AFTER UCH A PERIOD, THE FOREST EDGE IS
        ! ASSUMED TO GET ACCLIMATIZED TO THE NEW CIRCUMSTANCES, AND IT'S VULNERABILITY
        ! DIMINISHES (see Persson, 1975; Valinger and Fridman, 2011) 
        ! Note that gap area for the current year is not calculated yet (which is
        ! calculated at the end of the year in stomate_season), so the gap  created t-1 
        ! to t-5 is applied.
        ! Calculating the total area of timber removals in the previous five years 
        ! (set by the variable ::legacy_years_wind)
        ! +++NOTE After adding land_cover_change in relation to the disturbances, there
        ! becomes a very small chance to create a gap. In the future, Need to think about
        ! adding notion gap/edge to the land cover change 

        ! +++NOTE: Improvements needed for the calculation of closer/further area.
        ! Cases have been identified where the timber removal area exceeds
        ! the total area due to veget_max movement. The area for timber
        ! removal is calculated based on changes in basal area, incorporating
        ! veget_max. As a result, if a large portion of veget_max is moved,
        ! the calculated removal area can exceed the remaining total area.
        !
        ! The closer/further area assumes that timber removal is conducted within
        ! the current cohort and calculates the gap accordingly. However, this
        ! does not apply when the removed area is transferred to another PFT. If
        ! the timber removal area exceeds the total area, it could be
        ! interpreted as the entire area of the PFT being "closer to gap." While
        ! this might seem plausible, it does not feel like an accurate
        ! calculation and requires further consideration.
 
        area_timber_removals_5_years = zero
        DO iyear = 1, legacy_years_wind
            area_timber_removals_5_years = area_timber_removals_5_years + &
                gap_area_save(ipts,ivm,iyear)
        ENDDO

        ! Number of gaps with the area of the maximum allowed clearfelling:
        n_gap  = area_timber_removals_5_years / clear_cut_max 
          
        ! The length of the sides of the squared-shaped gaps:
        length_side_gap = SQRT(clear_cut_max)

        !! 1.6 Calculating the current tree spacing in the PFT
        ! The current tree spacing between tree stems are calculated from the number of
        ! stems per hectare. circ_class_n is the number of stems per m2; it's converted 
        ! to stems per hectare with m2_to_ha which has a value of 10000.
        ! The virtual tree spacing assumes that the LAI remains constant and the
        ! calculates how many trees of a given circumference class would be needed to
        ! result in the set LAI. Because lai changes throughout the season because of 
        ! leaf fall, wood mass was used instead. Because LAI is related to wood mass
        ! (see allocation). The principle remains the same 
        DO icir = 1,ncirc
            IF (  (circ_class_biomass(ipts,ivm,icir,isapabove,icarbon) + &
               circ_class_biomass(ipts,ivm,icir,iheartabove,icarbon)) .GT. zero ) THEN
                virtual_circ_class_n(ipts,ivm,icir) =  &
                (biomass(ipts,ivm,isapabove,icarbon) + &
                biomass(ipts,ivm,iheartabove,icarbon)) / &
                (circ_class_biomass(ipts,ivm,icir,isapabove,icarbon) + &
                circ_class_biomass(ipts,ivm,icir,iheartabove,icarbon))
            ELSE                 
                virtual_circ_class_n(ipts,ivm,icir) = zero 
            ENDIF

            ! Guillaume M. -- Inventory-visible density keeps only the circumference
            ! classes whose dbh reaches dia_thresh_inv.
            ccn_inv(icir) = circ_class_n(ipts,ivm,icir)
            IF(mean_dbh(icir) .LT. dia_thresh_inv(ivm) ) THEN
                ccn_inv(icir) = zero
            ENDIF
        ENDDO

        ! Guillaume M. -- When every circumference class sits below the inventory
        ! threshold, SUM(ccn_inv) is zero and the spacing below divides by zero:
        ! fall back to the largest class.
        IF(SUM(ccn_inv(:)) == zero) THEN
            ccn_inv(ncirc) = circ_class_n(ipts,ivm,ncirc)
        ENDIF

        ! Now the spacing is calculated based on the relative density index based on the
        ! virtual circ class to better deal with unmanaged forest and heterogenous canopy.      
        IF (virtual_circ_class_n(ipts,ivm,ncirc) .GT. min_stomate) THEN
            ! The original equation is current_spacing = 100.0 / &
            ! SQRT(m2_to_ha * virtual_circ_class_n(ipts,ivm,icir))
            ! given that m2_to_ha = 10000 it reduces to:
            !     current_spacing(ipts,ivm) = 1.0 / &
            !            SQRT(virtual_circ_class_n(ipts,ivm,ncirc))
            ! The above spacing calculation worked with previous version of
            ! ORCHIDEE, but with the developments, e.g., dynamic height and
            ! alpha self-thinning this spacing calculation gives too low CWS
            ! with high damage. Thus calculate current_spacing based on all
            ! number of trees not only biggest class.
            ! The previous spacing calculation (above) was compatible with earlier versions of
            ! ORCHIDEE. However, due to recent developments (e.g., dynamic tree height and alpha
            ! self-thinning), it now underestimates the critical wind speed (CWS), leading to
            ! high damage. To address this, we now compute `current_spacing` using the total number of
            ! trees, not just those in the largest diameter class.
            ! Guillaume M. -- /!\ The gate above tests virtual_circ_class_n but the divisor
            ! is SUM(ccn_inv): two DIFFERENT arrays. A slot whose ccn_inv sums to zero
            ! (e.g. freshly purged phantom slot) slipped through and put an Inf in the
            ! CURRENT_SPACING field, which XIOS refuses to write (numeric conversion not
            ! representable). The MAX floor is bit-neutral whenever the sum is >= min_stomate.
            current_spacing(ipts,ivm) = 1.0 / &
                        SQRT(MAX(SUM(ccn_inv(:)), min_stomate))
        ELSE
            ! In principle, we should never enter this condition. The spacing between 
            ! 2 trees equals the size of the pixel. It is thus empty. We defined this
            ! case such that if we do end here, a very large cws will be calculated and 
            ! the code should not suffer from divide by zero.
            current_spacing(ipts,ivm) = SQRT(area(ipts))
            WRITE(numout,*) 'WARNING from WINDFALL module:'
            WRITE(numout,*) 'the current spacing set to pixel spacing:', current_spacing(ipts,ivm)
            CALL ipslerr_p(2,'Wind throw module',&
                'current spacing was set to pixel dimensions',&
                'it is therefore assumed to be empty','') 
        ENDIF

        qm_dia(ipts,ivm) = wood_to_qmdia(circ_class_biomass(ipts,ivm,:,:,icarbon),&
                 circ_class_n(ipts,ivm,:), ivm, pipe_tune2(ipts,ivm))
        rdi(ipts,ivm) = (SUM(circ_class_n(ipts,ivm,:))* &
             m2_to_ha)/Nmax(qm_dia(ipts,ivm)*m_to_cm,alpha_self_thinning(ipts,ivm),&
             ivm, forest_managed(ipts,ivm))


        IF(test_spacing) THEN
           ! Temporary adjustment of circ_class_n to represent young forest
           ! dynamics.
           ! - In principle, CWS is higher (i.e. stands are less vulnerable to
           ! wind damage) 
           !   when they are young, due to:
           !   (1) closer spacing between trees,
           !   (2) lower wind speeds at lower canopy height,
           !   (3) greater stem flexibility (not represented in ForestGALES).
           ! - Current wind module (after rXXXX) accounts for (2), but initializes young
           !   stands too sparsely. This accelerates early growth, but makes
           !   young trees too vulnerable to wind.
           ! - A better long-term solution would be to:
           !   * initialize stands with more individuals at establishment, and
           !   * increase the RDI target for young stands.
           !   However, implementing this would require deeper structural 
           !   changes in forest, which are difficult to tune and would strongly affect
           !   overall model performance.
           !
           ! Workaround:
           ! - circ_class_n is modified here to artificially increase the
           !   number of trees when stands are small. This enforces higher
           !   early stand density, leading to higher CWS and more
           !   wind resistance for young forests.
           ! Guillaume M. -- Same MAX floor as above: this branch OVERWRITES the guarded
           ! value and would reintroduce the Inf on an empty slot.
           current_spacing(ipts,ivm) = 1.0/ SQRT(MAX(SUM(ccn_inv(:)), min_stomate)** &
                        MIN(mean_dbh(ncirc)/(largest_tree_dia(ivm)*d_threshold_wind),un))
        ENDIF
        
        IF(printlev_loc.GE.4 .AND. ivm .EQ. test_pft .AND. ipts .EQ. test_grid) THEN
            WRITE(numout,*) 'ivm',ivm
            WRITE(numout,*) 'ccn', circ_class_n(ipts,ivm,:)
            WRITE(numout,*) 'v_ccn', virtual_circ_class_n(ipts,ivm,:)
            WRITE(numout,*) 'spacing_virtual:',1.0 / SQRT(virtual_circ_class_n(ipts,ivm,:))
            WRITE(numout,*) 'rooting_depth_type,(1:deep,2:shallow,3:average)',rooting_depth
            WRITE(numout,*) 'biomass',circ_class_biomass(ipts,ivm,ncirc,:,icarbon)
            WRITE(numout,*) 'wind_in',wind_in(ipts)
            WRITE(numout,*) 'wind_mod', wind_speed_daily(ipts)
            WRITE(numout,*) 'wind_max',max_wind_speed_storm(ipts)
            WRITE(numout,*) 'current_spacing',current_spacing(ipts,ivm)
            WRITE(numout,*) 'rdi',rdi
        ENDIF 
        
        !! 1.7 Calculate green stem mass
        ! Stem mass is calculated as the aboveground woody biomass minus the branches. 
        ! circ_class_biomass is the biomass of an individual tree (gC/tree)
        ! Here, we convert the unit to (kgC/tree) by dividing kilo_to_unit (1000)
        stem_mass(ipts,ivm) = (circ_class_biomass(ipts,ivm,ncirc,isapabove,icarbon) + &
            circ_class_biomass(ipts,ivm,ncirc,iheartabove,icarbon)) * &
            (un - branch_ratio(ivm)) / kilo_to_unit

        !-Debug-
        IF (printlev_loc .GT. 4 .AND. ipts .EQ. test_grid .AND. ivm .EQ. test_pft) THEN
            WRITE(numout,*) 'stem_mass,',stem_mass(ipts,ivm)
            WRITE(numout,*) 'green_density,', green_density(ivm)
            WRITE(numout,*) 'pipe_sensity,', pipe_density(ivm)
            WRITE(numout,*) 'kilo to unit,', kilo_to_unit
            WRITE(numout,*) 'virtual_circ_class', virtual_circ_class_n(ipts,ivm,ncirc)
            WRITE(numout,*) 'sla,',sla(ivm)
        ENDIF

        ! Wind throw occurs on living trees so the stem mass should be converted from dry
        ! biomass in gC to wet biomass in kg to be used in the GALES equations. 
        ! Stem_mass thus converted as follows: (green_density / (pipe_density/1000)) 
        ! from dry wood mass (kgC/tree) to wet wood mass (Kg/trees)
        stem_mass(ipts,ivm) = stem_mass(ipts,ivm) * &
            green_density(ivm)/(pipe_density(ivm)/kilo_to_unit)
   
        !! 1.8 Calculating leaf area index 
        ! We need to use the virtual biomass in each circumference class to calculate lai, 
        ! and tree_heights_from_edge was calculated as function of lai. 
        ! As CWS is calculated only for the biggest circ_class, lai for the biggest circ
        ! class is calculated using the function biomass_to_lai to account for a dynamic 
        ! lai calculation.
        lai = biomass_to_lai(circ_class_biomass(ipts,ivm,ncirc,ileaf,icarbon)* &
                circ_class_n(ipts,ivm,ncirc),ivm)

        !-Debug-
        IF (printlev_loc .GT. 4 .AND. ipts .EQ. test_grid .AND. ivm .EQ. test_pft) THEN           
            WRITE(numout,*) 'check spacing and circ_number in a biggest class'
            WRITE(numout,*) 'virtual_tree_numbers:', &
                virtual_circ_class_n(ipts,ivm,ncirc)
            WRITE(numout,*) 'current spacing: ',current_spacing(ipts,ivm)
            WRITE(numout,*) 'area:', area(ipts)
            WRITE(numout,*) 'ipts, ',ipts,'ivm, ',ivm
            WRITE(numout,*) 'mean dbh, ', mean_dbh(:)
            WRITE(numout,*) 'mean height, ', mean_height(:)
            WRITE(numout,*) 'stem_mass, ', stem_mass(ipts,ivm)
            WRITE(numout,*) 'crown_volume, ',crown_volume
            WRITE(numout,*) 'crown_mass, ',crown_mass
            WRITE(numout,*) 'lai, ',lai
            WRITE(numout,*) 'virtual_circ_n, ', virtual_circ_class_n(ipts,ivm,:) 
            WRITE(numout,*) 'real_circ_n, ', circ_class_n(ipts,ivm,:) 
            WRITE(numout,*) 'End of Section 1:  Calculated canopy structure &
                     & in difffernt circumferences for each PFT' 
        ENDIF

        !! 2. Calculating the gustiness of wind in the PFTs

        !! 2.1 Calculate gaps and the surrounding areas 
        ! IMPORTANT: AS WE DO NOT HAVE INFORMATION ON STAND EDGES IN THE ORCHIDEE 
        ! SIMULATIONS, A POSSIBLE PROXY TO USE HERE IS A VARIABLE THAT REPRESENTS 
        ! EARLIER HARVESTS IN THE PREVIOUS 5 YEARS. AFTER SUCH A PERIOD, THE FOREST
        ! EDGE IS ASSUMED TO GET ACCLIMATIZED TO THE NEW CIRCUMSTANCES, AND IT'S  
        ! VULNERABILITY DIMINISHES (see Persson, 1975; Valinger and Fridman, 2011)
        ! There is a big increase in gustiness after a certain distance from the edge of
        ! the forest, hence we divide the area of the PFT into two areas:
        !  - the area which is further than X tree heights from any edge that was created 
        !    in the previous five years (area_total_further)
        !  - the area which is closer than X tree heights from any edge that was created 
        !    in the previous five years (area_total_closer).
        ! Of course, these smaller areas depend on the fragmentation of the gaps (create 
        ! in the last five years). As we have no information on this, we assume a fixed 
        ! gap size representing the most likely size of clear cuts (clear_cut_max) and 
        ! that the gaps are square-shaped. Then, we divide the area of harvest in the 
        ! previous five years (area_timber_removals_5_years) by clear_cut_max. 
        ! Hence, from the number of gaps we can calculate the area_total_closer. 
        ! As each circumference class will have its forest area divided, we calculate two 
        ! gustiness values and hence, two critical wind speeds.
        ! The value of X depends on the leaf area index of the PFT and is calculated 
        ! as follows: tree_heights_from_edge is defined as "X/h", X sets to 28/LAI and h
        ! is set to meant_tree_height.
        ! This comes from the equation of calculating canopy penetration depth in
        ! Pan et al. 2017?. Ying Pan work on the large eddy simulation and is now at 
        ! Pennsylvania State University   

        ! Here, lai for each diameter class is calculated by virtual_class_n 
        ! a threshold "0.1" is used to avoid the issue of divided by small lai < 0.1 
        ! for the cases of deciduous tree species 
        IF (lai .GT. 0.1) THEN
            tree_heights_from_edge = (28.0/lai) / mean_height(ncirc) 
        ELSE
            tree_heights_from_edge = (28.0/0.1) / mean_height(ncirc) 
        ENDIF
         
        ! As suggested by Barry (the father of GALES), the tree height from edge 
        ! should be limited to 9.0     
        IF (tree_heights_from_edge .GT. 9.0)  tree_heights_from_edge = 9.0
         
        ! The frame-shaped area around the clear cut gap with X as the width of the frame:
        area_around_gap = ((length_side_gap + &
            tree_heights_from_edge * mean_height(ncirc)*2.0)**2.0) - &
            clear_cut_max

        ! Total area of the "frames" in the particular dbh-class of the particular PFT: 
        ! the multiplication by 0.25 is because of the assumption that wind direction
        ! does not change during the half-hourly time step, therefore only one of the 
        ! four sides of the gaps are affected.
        ! Use MIN to make sure area_total_closer and area_total_further are positive.
        !area_total_closer = MAX(MIN(n_gap * area_around_gap * 0.25, &
        !                      area(ipts)*veget_max(ipts,ivm) - (n_gap * clear_cut_max)),&
        !                      zero)

        ! Previously, the maximum possible area_closer was assumed to be 
        ! the total area minus the harvested area. However, as mentioned earlier, 
        ! if the timber removal area exceeds the total area due to veget_max 
        ! movement, area_closer can become negative. 
        !
        ! In such cases, the total area could instead be considered as area_closer. 
        ! (See the note above the calculation of n_gap) 
        ! For now, the maximum area is set to the total area, but this approach 
        ! requires further refinement and improvement.
        area_total_closer = MAX(MIN(n_gap * area_around_gap * 0.25, &
                              area(ipts)*veget_max(ipts,ivm) ),&
                              zero)

        ! The total area which is nor gap neither "frame":
        ! In principle, area_total_further should be greater than zero. However, in
        ! runs with multiple age classes, area can change due to veget_max movement,
        ! while area_timber_removal remains unchanged. Consider whether area_timber_removal
        ! should also be adjusted along with veget_max.
        area_total_further = MAX(area(ipts) * veget_max(ipts,ivm) - &
                            area_total_closer - (n_gap * clear_cut_max),zero)

        IF (area(ipts).GT.zero) THEN
            ! Get the area fraction
            ! Area fraction is used to divide biomass for removal in the module `cal_storm_circ_class_kill`. 
            ! As 1) circ_class_kill is based on the proportion of killing biomass to standing biomass;  
            ! 2) It is assumed there is no biomass in gaps, therefore, the sum of the fractions should equal one.  
            ! Otherwise, it can undereste biomass to kill at the circ level, which  can lead to all individuals 
            ! being removed (see the module `cal_storm_circ_class_kill`).  
            ! To address this, a subtraction of the timber removal area has been added.  
            closer_area_fraction(ipts) = MIN(un,area_total_closer / (area(ipts)*veget_max(ipts,ivm) - &
                 (n_gap * clear_cut_max)))
            further_area_fraction(ipts) = MAX(zero, un - closer_area_fraction(ipts))
        ELSE

            WRITE(numout,*) 'ERROR: the total area of the pixel is zero'
            CALL ipslerr_p(3,'stomate_windthrow',&
                    'area of the pixel is zero','','')
        ENDIF

        !-Debug-
        IF(printlev_loc .GE. 4 .AND. ipts .EQ. test_grid .AND. ivm .EQ. test_pft)THEN 
            WRITE(numout,*) 'Section 2.0: Calculating the gustiness of wind &
                      & in the PFT and area fractions in the grid.'
            WRITE(numout,*) 'n_gap',n_gap
            WRITE(numout,*) 'crown_width, ', crown_width(:)
            WRITE(numout,*) 'canopy_height, ', mean_height(:)
            WRITE(numout,*) 'total area',(area(ipts)*veget_max(ipts,ivm))
            WRITE(numout,*) 'gap_area',gap_area_save(ipts,ivm,:)
            WRITE(numout,*) 'timber removal', n_gap * clear_cut_max
            WRITE(numout,*) 'area_closer_fraction, ', closer_area_fraction(ipts)
            WRITE(numout,*) 'area_further_fraction, ', further_area_fraction(ipts)
            WRITE(numout,*) 'area_total_closer',area_total_closer
            WRITE(numout,*) 'area_total_further',area_total_further
        ENDIF

        !! 2.2 Calculate ratio to spacing height 
        ! The ratio of spacing to height describes the functional (for wind) stand
        ! density. If there are many tall trees, the spacing is low and thus ratio_
        ! spacing_height will be low as well. Values of the ratio of mean tree spacing
        ! and height are limited between 0.075 and 0.45. Bigger values don't make any 
        ! difference, smaller values are not really possible in real life.
        ! The D/h  value are presented in Hale et al., 2015
        ratio_spacing_height = MAX(MIN(current_spacing(ipts,ivm) / &
                                       mean_height(ncirc),0.45),0.075)

        ! Copy the ratio_spacing_height to the D_H_ratio for output 
        d_h_ratio(ipts,ivm) = ratio_spacing_height
         
        !! 2.3 Calculate functional gap size
        ! Calculate whether the gap is large enough to increase the sensitivity of the
        ! remaining trees to wind throw. GALES uses the length of the edge of the 
        ! square-shaped gaps calculated in the previous section.
        gap_size = length_side_gap
 
        ! If the gap is wider than 10 times the tree height, then wind loading on the 
        ! stand edge does not increase any more. 
        ! The threshold of 10 is used in the equations below
        IF (gap_size > 10.0 * mean_height(ncirc)) THEN
            mean_gap_factor = 1.
            max_gap_factor = 1. 
        ELSE
            mean_gap_factor = (0.001 * (gap_size / &
                mean_height(ncirc))**0.562) / (0.001 * 10.**0.562)
            max_gap_factor = (0.0064 * (gap_size / &
                mean_height(ncirc))**0.3467) / (0.0064 * 10.**0.3467)
        END IF
         
        ! In the following section, there are three types of terms for forest edges:
        ! the edge (edge), the area around the edge (closer)  and the area further away
        ! from the edge (further). Only closer and further are used in the calculation
        ! of cws. Gustiness is different in both areas. 

        ! Calculating term_a and term_b for the new gust method. The definition of the
        ! term_a and term_b are in the Hale  et al. (2015) eq(7) 
        ! D/h is the ratio_spacing_height
        ! x/h is the tree_heights_from_edge
        term_a = MAX(-2.1 * ratio_spacing_height + 0.91, zero) 
        term_b = 1.0611 * LOG(ratio_spacing_height) + 4.2
        
        ! New gust factors:
        ! The gustiness away from the edge is always higher than at the edge.
        ! We applied term_a*tree_heights_from_edge  + term_b for gustiness calculation.
        ! The gustiness away from the edge is always higher than at the edge.
        ! For the further area, tree_height_from_edge is always greater than 9.0 thus,
        ! the above calculation will produce an even higher value for gustiness
        ! calculation. So, We set tree_heights_from_edge as 9.0.
        ! The gustiness at the further area will than be higher then the gustiness at 
        ! the closer area.        
        new_gust_closer  = term_a * tree_heights_from_edge + term_b
        new_gust_further = term_a * 9.0  + term_b
        new_gust_edge    = term_a * zero + term_b
       
        ! Calculating gustiness:
        g_closer(ipts,ivm) = new_gust_closer
        g_further(ipts,ivm) = new_gust_further

        ! The following section is used to calculate the edge_factor. This section 
        ! calculates how the mean loading on a tree changes as a function of tree height, 
        ! spacing, distance from edge and size of any upwind gap. An upwind gap greater 
        ! than 10 tree heights is assumed to be an infinite gap. The output (edge_factor)
        ! is used to modify the calculated wind loading on the tree, which assumes a 
        ! steady wind and a position well inside the forest.
        mean_bm_gust_factor_closer = (0.68 * ratio_spacing_height - &
            0.0385) + (-0.68 * ratio_spacing_height + 0.4785) * ((1.7239 * &
            ratio_spacing_height + 0.0316)**tree_heights_from_edge) * &
            mean_gap_factor
        mean_bm_gust_factor_further = (0.68 * ratio_spacing_height - &
            0.0385) + (-0.68 * ratio_spacing_height + 0.4785) * ((1.7239 * &
            ratio_spacing_height + 0.0316)**9.0) * mean_gap_factor
        edge_factor = mean_bm_gust_factor_closer / mean_bm_gust_factor_further
        ! Copy array for output 
        edge_factor_out(ipts,ivm)  = edge_factor
        ! Calculating gust_factor accounting for gustiness and edge_factor:
        gust_factor_closer(ipts,ivm) = g_closer(ipts,ivm) * edge_factor
           
        ! Suggested by Barry
        ! For further we force edge_factor = 1. This is because everything is
        ! compared to wind load on tree in the middle of the forest
        edge_factor = 1.0 
        gust_factor_further(ipts,ivm) = g_further(ipts,ivm) *  edge_factor
         
        !-Debug-
        IF(printlev_loc .GT. 4 .AND. ipts .EQ. test_grid .AND. ivm .EQ. test_pft)THEN 
            WRITE(numout,*) 'End of Section 2.0: Calculating the &
                      & gustiness of wind in the PFT'
            WRITE(numout,*) 'term_a, ', term_a
            WRITE(numout,*) 'term_b, ', term_b
            WRITE(numout,*) 'new_gust_closer, ', new_gust_closer
            WRITE(numout,*) 'ratio_spacing_height, ', ratio_spacing_height
            WRITE(numout,*) 'new_gust_further, ', new_gust_further
            WRITE(numout,*) 'new_gust_edge, ', new_gust_edge 
            WRITE(numout,*) 'tree_heights_from_edge, ', tree_heights_from_edge
            WRITE(numout,*) 'Correction factor for the forest &
                      & edge dimension in number of circumference'
            WRITE(numout,*) 'edge_factor, ', edge_factor
            WRITE(numout,*) 'gustness_closer, ', g_closer(ipts,ivm)
            WRITE(numout,*) 'gustness_futher, ', g_further(ipts,ivm) 
            WRITE(numout,*) 'gust_factor_closer, ' , gust_factor_closer(ipts,ivm)
            WRITE(numout,*) 'gust_factor_further, ', gust_factor_further(ipts,ivm)
            WRITE(numout,*) 'Mean bm_gust_factor_closer:', mean_bm_gust_factor_closer
            WRITE(numout,*) 'Mean bm_gust_factor_further:', mean_bm_gust_factor_further
            WRITE(numout,*) 'Mean gap_factor:', mean_gap_factor
            WRITE(numout,*) 'plant_status',plant_status(ipts,ivm)
        ENDIF
                  
        !! 3. Calculating critical moments (tree mechanics section)
        ! overturning_moment_multiplier: this comes form the species parameters file and 
        ! is related to the depth and type of the soil. Rooting depth has three categories: 
        ! Shallow (less then 0.8 m deep); Deep (deeper than 0.8 m); Average (average values 
        ! to be used when rooting depth is unknown). Soil type has four categories: 
        ! Free-draining mineral soils; Gleyed mineral soils; Peaty mineral soils; Deep peats. 
        ! As these are generic soil types, a modified soil map of Europe would be needed 
        ! for the optimal use of WINDTHROW.
        IF (soil_type(ipts) == ifree_draining .AND. rooting_depth == ishallow &
                  .AND. (plant_status(ipts,ivm).EQ.icanopy .OR. &
                         plant_status(ipts,ivm).EQ.ipresenescence)) THEN
            overturning_moment_multiplier = overturning_free_draining_shallow(ivm)
        ELSE IF (soil_type(ipts) == ifree_draining .AND. rooting_depth == ishallow &
                  .AND. (plant_status(ipts,ivm).NE.icanopy .OR. &
                         plant_status(ipts,ivm).NE.ipresenescence)) THEN
            overturning_moment_multiplier = overturning_free_draining_shallow_leafless(ivm)
        ELSE IF (soil_type(ipts) == ifree_draining .AND. rooting_depth == ideep &
                  .AND. (plant_status(ipts,ivm).EQ.icanopy .OR.&
                         plant_status(ipts,ivm).EQ.ipresenescence)) THEN
            overturning_moment_multiplier = overturning_free_draining_deep(ivm)
        ELSE IF (soil_type(ipts) == ifree_draining .AND. rooting_depth == ideep &
                  .AND. (plant_status(ipts,ivm).NE.icanopy .OR. &
                         plant_status(ipts,ivm).NE.ipresenescence)) THEN
            overturning_moment_multiplier = overturning_free_draining_deep_leafless(ivm)
        ELSE IF (soil_type(ipts) == ifree_draining .AND. rooting_depth == iaverage &
                  .AND. (plant_status(ipts,ivm).EQ.icanopy .OR. &
                         plant_status(ipts,ivm).EQ.ipresenescence)) THEN                
            overturning_moment_multiplier = overturning_free_draining_average(ivm)
        ELSE IF (soil_type(ipts) == ifree_draining .AND. rooting_depth == iaverage &
                  .AND. (plant_status(ipts,ivm).NE.icanopy .OR. &
                         plant_status(ipts,ivm).NE.ipresenescence)) THEN
            overturning_moment_multiplier = overturning_free_draining_average_leafless(ivm)
        ELSE IF (soil_type(ipts) == igleyed .AND. rooting_depth == ishallow &
                  .AND. (plant_status(ipts,ivm).EQ.icanopy .OR. &
                         plant_status(ipts,ivm).EQ.ipresenescence)) THEN                
            overturning_moment_multiplier = overturning_gleyed_shallow(ivm)
        ELSE IF (soil_type(ipts) == igleyed .AND. rooting_depth == ishallow &
                  .AND. (plant_status(ipts,ivm).NE.icanopy .OR. &
                         plant_status(ipts,ivm).NE.ipresenescence)) THEN
            overturning_moment_multiplier = overturning_gleyed_shallow_leafless(ivm)
        ELSE IF (soil_type(ipts) == igleyed .AND. rooting_depth == ideep &
                  .AND. (plant_status(ipts,ivm).EQ.icanopy .OR. &
                         plant_status(ipts,ivm).EQ.ipresenescence)) THEN                
            overturning_moment_multiplier = overturning_gleyed_deep(ivm)
        ELSE IF (soil_type(ipts) == igleyed .AND. rooting_depth == ideep &
                  .AND. (plant_status(ipts,ivm).NE.icanopy .OR. &
                         plant_status(ipts,ivm).NE.ipresenescence)) THEN
            overturning_moment_multiplier = overturning_gleyed_deep_leafless(ivm)
        ELSE IF (soil_type(ipts) == igleyed .AND. rooting_depth == iaverage &
                  .AND. (plant_status(ipts,ivm).EQ.icanopy .OR. &
                         plant_status(ipts,ivm).EQ.ipresenescence)) THEN
            overturning_moment_multiplier = overturning_gleyed_average(ivm)
        ELSE IF (soil_type(ipts) == igleyed .AND. rooting_depth == iaverage &
                  .AND. (plant_status(ipts,ivm).NE.icanopy .OR. &
                         plant_status(ipts,ivm).NE.ipresenescence)) THEN
            overturning_moment_multiplier = overturning_gleyed_average_leafless(ivm)
        ELSE IF (soil_type(ipts) == ipeaty .AND. rooting_depth == ishallow &
                  .AND. (plant_status(ipts,ivm).EQ.icanopy .OR. &
                         plant_status(ipts,ivm).EQ.ipresenescence)) THEN
            overturning_moment_multiplier = overturning_peaty_shallow(ivm)
        ELSE IF (soil_type(ipts) == ipeaty .AND. rooting_depth == ishallow &
                  .AND. (plant_status(ipts,ivm).NE.icanopy .OR. &
                         plant_status(ipts,ivm).NE.ipresenescence)) THEN
            overturning_moment_multiplier = overturning_peaty_shallow_leafless(ivm)
        ELSE IF (soil_type(ipts) == ipeaty .AND. rooting_depth == ideep &
                  .AND. (plant_status(ipts,ivm).EQ.icanopy .OR. &
                         plant_status(ipts,ivm).EQ.ipresenescence)) THEN
            overturning_moment_multiplier = overturning_peaty_deep(ivm)
        ELSE IF (soil_type(ipts) == ipeaty .AND. rooting_depth == ideep &
                  .AND. (plant_status(ipts,ivm).NE.icanopy .OR. &
                         plant_status(ipts,ivm).NE.ipresenescence)) THEN
            overturning_moment_multiplier = overturning_peaty_deep_leafless(ivm)
        ELSE IF (soil_type(ipts) == ipeaty .AND. rooting_depth == iaverage &
                  .AND. (plant_status(ipts,ivm).EQ.icanopy .OR. &
                         plant_status(ipts,ivm).EQ.ipresenescence)) THEN
            overturning_moment_multiplier = overturning_peaty_average(ivm)
        ELSE IF (soil_type(ipts) == ipeaty .AND. rooting_depth == iaverage &
                  .AND. (plant_status(ipts,ivm).NE.icanopy .OR. &
                         plant_status(ipts,ivm).NE.ipresenescence)) THEN
            overturning_moment_multiplier = overturning_peaty_average_leafless(ivm)
        ELSE IF (soil_type(ipts) == ipeat .AND. rooting_depth == ishallow &
                  .AND. (plant_status(ipts,ivm).EQ.icanopy .OR. &
                         plant_status(ipts,ivm).EQ.ipresenescence)) THEN
            overturning_moment_multiplier = overturning_peat_shallow(ivm)
        ELSE IF (soil_type(ipts) == ipeat .AND. rooting_depth == ishallow &
                  .AND. (plant_status(ipts,ivm).NE.icanopy .OR. &
                         plant_status(ipts,ivm).NE.ipresenescence)) THEN
            overturning_moment_multiplier = overturning_peat_shallow_leafless(ivm)
        ELSE IF (soil_type(ipts) == ipeat .AND. rooting_depth == ideep &
                  .AND. (plant_status(ipts,ivm).EQ.icanopy .OR. &
                         plant_status(ipts,ivm).EQ.ipresenescence)) THEN
            overturning_moment_multiplier = overturning_peat_deep(ivm)
        ELSE IF (soil_type(ipts) == ipeat .AND. rooting_depth == ideep &
                  .AND. (plant_status(ipts,ivm).NE.icanopy .OR. &
                         plant_status(ipts,ivm).NE.ipresenescence)) THEN
            overturning_moment_multiplier = overturning_peat_deep_leafless(ivm)
        ELSE IF (soil_type(ipts) == ipeat .AND. rooting_depth == iaverage &
                  .AND. (plant_status(ipts,ivm).EQ.icanopy .OR. &
                         plant_status(ipts,ivm).EQ.ipresenescence)) THEN
            overturning_moment_multiplier = overturning_peat_average(ivm)
        ELSE IF (soil_type(ipts) == ipeat .AND. rooting_depth == iaverage &
                  .AND. (plant_status(ipts,ivm).NE.icanopy .OR. &
                         plant_status(ipts,ivm).NE.ipresenescence)) THEN
            overturning_moment_multiplier = overturning_peat_average_leafless(ivm)
        ELSE
            !average value if parameters are missing
            overturning_moment_multiplier = 125.0
        END IF
             
        ! Maximum moment calculated for tree overturning.
        max_overturning_moment(ipts,ivm) = overturning_moment_multiplier * &
            stem_mass(ipts,ivm)
           
        ! Maximum moment calculated for stem breakage.
        max_breaking_moment(ipts,ivm) = modulus_rupture(ivm) * f_knot(ivm) * pi * &
            (mean_dbh(ncirc)**3.0) / 32.0
           
        ! This is where the gustiness of wind is accounted for so that half hourly mean 
        ! wind speeds can be used in the simulations.
        ! In the end, moments and hence, critical wind speeds are calculated for both
        ! areas (closer & further) with different gustiness.
        overturning_moment_closer(ipts,ivm) = max_overturning_moment(ipts,ivm) / &
            gust_factor_closer(ipts,ivm)
        overturning_moment_further(ipts,ivm) = max_overturning_moment(ipts,ivm) / &
            gust_factor_further(ipts,ivm)
        breaking_moment_closer(ipts,ivm) = max_breaking_moment(ipts,ivm) / &
            gust_factor_closer(ipts,ivm)
        breaking_moment_further(ipts,ivm) = max_breaking_moment(ipts,ivm) / &
            gust_factor_further(ipts,ivm)
      
        !-Debug-
        !IF(printlev_loc .GE. 4 .AND. ipts .EQ. test_grid .AND. ivm .EQ. test_pft)THEN
        !    WRITE(numout,*) 'SENSCENCE type,', plant_status(ipts,ivm)
        !    WRITE(numout,*) 'Overturning multiplier for free drainage average root leafless all PFTs,', &  
        !                        overturning_free_draining_average_leafless(ivm)
        !    WRITE(numout,*) 'Overturning multiplier for free drainage average root leaf all PFTs,', &  
        !                        overturning_free_draining_average(ivm)
        !    WRITE(numout,*) 'Overturning multiplier for free drainage shallow root leafless all PFTs,', &  
        !                        overturning_free_draining_shallow_leafless(ivm)
        !    WRITE(numout,*) 'Overturning multiplier for free drainage shallow root leaf all PFTs,', &  
        !                        overturning_free_draining_shallow(ivm)
        !    WRITE(numout,*) 'Overturning multiplier for free drainage deep root leafless all PFTs,', &  
        !                        overturning_free_draining_deep_leafless(ivm)
        !    WRITE(numout,*) 'Overturning multiplier for free drainage deep root leaf all PFTs,', &  
        !                        overturning_free_draining_deep(ivm)
        !    WRITE(numout,*) 'End of Section 3 : Calculating critical moments (tree mechanics section)'
        !    WRITE(numout,*) 'Overturning moment multiplier, ', overturning_moment_multiplier
        !    WRITE(numout,*) 'maximum moment for overturning, ', max_overturning_moment(ipts,ivm)
        !    WRITE(numout,*) 'maximum moment for stem breakag, ', max_breaking_moment(ipts,ivm)
        !    WRITE(numout,*) 'overturning closer,', overturning_moment_closer(ipts,ivm)
        !    WRITE(numout,*) 'breaking closer, ', breaking_moment_closer(ipts,ivm)
        !    WRITE(numout,*) 'overturning further, ', overturning_moment_further(ipts,ivm)
        !    WRITE(numout,*) 'breaking further, ', breaking_moment_further(ipts,ivm)
        !ENDIF
      
        !! 4. Calculating critical wind speed and damage rate

        !! 4.1 Calculating critical wind speed in the forest and 10 m above the zero-plane displacement   
        ! +++CHECK+++
        ! Ideally ORCHIDEE should use only one surface roughness (z0) and displacement 
        ! height. This is not straightforward because here surface roughness and 
        ! displacement height are calculated as a function of wind speed to account for 
        ! streamlining of the canopy under strong winds. This is done in the subroutine 
        ! streamlining(...) in the calculate_wind_process(). Nevertheless, the subroutine 
        ! streamlining uses the canopy structure from ORCHIDEE and when storm damage
        ! occurs, the canopy structure is adjusted and this new structure will be used 
        ! in condveg.f90 . The inconsistency is thus between the roughness length 
        ! calculations in condveg.f90 and stomate_windthrow.f90. 
        ! +++++++++++
        ! Critical wind speed for stem breakage in the edge area (m/s)
        distance_from_force = 0.0 
        call calculate_wind_process(breaking_moment_closer(ipts,ivm), stem_mass(ipts,ivm), &
                mean_dbh(ncirc), current_spacing(ipts,ivm), &
                canopy_depth_dom, canopy_breadth_dom, canopy_height_dom, &
                streamlining_n, streamlining_c, cws_break_closer(ipts,ivm), &                      
                z0_wind_bk, zhd_wind_bk, distance_from_force, &
                gamma_solved_bk_closer(ipts,ivm), ivm, ipts, DLF_breaking_closer(ipts,ivm)) 
              
        ! Critical wind speed for stem breakage 10 m above 
        ! the zero-plane displacement in the edge area(m/s)
        call elevate(cws_break_closer(ipts,ivm), z0_wind_bk, zhd_wind_bk, & 
                canopy_height_dom, cws_break_closer_10(ipts,ivm))
       
        ! Critical wind speed for tree overturning in the edge area (m/s)
        distance_from_force = 0.0
        call calculate_wind_process(overturning_moment_closer(ipts,ivm), stem_mass(ipts,ivm), &
                mean_dbh(ncirc), current_spacing(ipts,ivm), &
                canopy_depth_dom, canopy_breadth_dom, canopy_height_dom, &
                streamlining_n, streamlining_c, cws_overturn_closer(ipts,ivm), &                      
                z0_wind_ov, zhd_wind_ov, distance_from_force, &
                gamma_solved_ov_closer(ipts,ivm), ivm, ipts, DLF_overturning_closer(ipts,ivm)) 
                             
        ! Critical wind speed for tree overturning 10 m above
        ! the zero-plane displacement in the edge area (m/s)
        call elevate(cws_overturn_closer(ipts,ivm), z0_wind_ov, zhd_wind_ov, & 
                canopy_height_dom, cws_overturn_closer_10(ipts,ivm))         

        ! IF the soil at 80 below ground is frozen, the module only allows the stem breakage
        IF ( soil_temp_daily(ipts) .GT. 273.15 ) THEN 
        ! Calculating damage severity in case any of the critical wind speeds is reached 
        ! in the area around gaps (forest edge area). Whichever of the two damage types 
        ! (overturning, stem breakage) has the lowest critical wind speed, it will be the 
        ! occurring one.
            IF (cws_overturn_closer_10(ipts,ivm) < cws_break_closer_10(ipts,ivm)) THEN
                cws_final_closer(ipts,ivm) = cws_overturn_closer_10(ipts,ivm)
                wind_damage_type_closer(ipts,ivm) = ioverturning
                z0_wind(ipts,ivm) = z0_wind_ov
                zhd_wind(ipts,ivm) = zhd_wind_ov
            ELSE
                cws_final_closer(ipts,ivm) = cws_break_closer_10(ipts,ivm)
                wind_damage_type_closer(ipts,ivm) = ibreakage
                z0_wind(ipts,ivm) = z0_wind_bk
                zhd_wind(ipts,ivm) = zhd_wind_bk
            END IF
        ELSE
            ! when the soil at 80 cm  is frozen, the only damage type is limited to the
            ! stem_breakage
                cws_final_closer(ipts,ivm) = cws_break_closer_10(ipts,ivm)
                wind_damage_type_closer(ipts,ivm) = ibreakage
                z0_wind(ipts,ivm) = z0_wind_bk
                zhd_wind(ipts,ivm) = zhd_wind_bk
        ENDIF

        ! Guillaume M. -- Wind acclimation does not scale the critical wind speed here: it
        ! enters as the Gompertz weight (weight_ratio) on the closer damage rate below.

        !max_wind_speed_pft_closer(ipts,ivm) = max_wind_speed_storm(ipts) * &
        !    LOG((elevate_wind + zhd_wind(ipts,ivm) - zhd_pixel(ipts))/z0_wind(ipts,ivm))/ &
        !    LOG(elevate_wind/z0_wind(ipts,ivm))
        actual_maxheight = pipe_tune2(ipts,ivm) * largest_tree_dia(ivm) ** &
            (pipe_tune3(ivm))
        IF(actual_maxheight .LT. canopy_height_dom) THEN
           actual_maxheight = canopy_height_dom
        ENDIF

!        max_wind_speed_pft_closer(ipts,ivm) = max_wind_speed_storm(ipts) * &
!            LOG(canopy_height_dom*height_displacement/z0_wind(ipts,ivm))/ &
!            LOG(actual_maxheight*height_displacement/z0_wind(ipts,ivm))
        wind_profile_closer(ipts,ivm) = &
            LOG((canopy_height_dom - canopy_height_dom*height_displacement)/z0_wind(ipts,ivm))/ &
            LOG((actual_maxheight - canopy_height_dom*height_displacement)/z0_wind(ipts,ivm))
        max_wind_speed_pft_closer(ipts,ivm) = max_wind_speed_storm(ipts) * wind_profile_closer(ipts,ivm)
        ! Guillaume M. -- The closer logistic uses the profile-adjusted wind only when
        ! wind_cal_pft is true, raw max_wind_speed_storm otherwise. The stored multiplier
        ! must mirror that choice.
        IF (.NOT. wind_cal_pft) wind_profile_closer(ipts,ivm) = un

        ! Based on observations of wind damage in Scotland and how much we had to adjust
        ! the critical wind speed to get agreement between the wind speed above the 
        ! canopy, the calculated critical wind speed and the presence of damage (Hale
        ! et al, 2015). I am suggesting the relationship below to calculate the percentage 
        ! of damage at a grid point as a function of the meteorological wind speed at 10m
        ! above the zero-plane displacement (prezhd_wind), and the critical wind speed 
        ! for damage (crit_wind) calculated by ForestGALES. I would set maximum to 80% 
        ! for the moment (sets the maximum total amount of damage = 80%) and the scaling 
        ! factor (s) to 0.8 (this affects the rate at which damage percentage increases
        ! with increasing wind speed). We can adjust these if the seem to be wrong or 
        ! I find other data that I can use to make a better assessment.
        ! We apply a sigmodial function to tune the "damage level" by using a function 
        ! which is dependent on the difference between the actual wind speed and the 
        ! final_cws (the lower one of the cws from overturning or stem_breakage. 
        IF(wind_cal_pft) THEN
           wind_damage_rate_closer(ipts,ivm) = max_damage_closer(ivm) * ( &
            (1./(1.+exp(-((max_wind_speed_pft_closer(ipts,ivm) -  &
            cws_final_closer(ipts,ivm)) / sfactor_closer(ivm))))) - &
            (1./(1.+exp(cws_final_closer(ipts,ivm)/sfactor_closer(ivm)))) )
        ELSE
           wind_damage_rate_closer(ipts,ivm) = max_damage_closer(ivm) * ( &
            (1./(1.+exp(-((max_wind_speed_storm(ipts) - cws_final_closer(ipts,ivm)) / &
            sfactor_closer(ivm))))) - &
            (1./(1.+exp(cws_final_closer(ipts,ivm)/sfactor_closer(ivm)))) ) 
        ENDIF
        ! Guillaume M. -- Apply the Gompertz wind-ratio weight to the closer damage
        ! magnitude, capped at max_damage_closer.
        wind_damage_rate_closer(ipts,ivm) = MIN( &
            max_damage_closer(ivm), wind_damage_rate_closer(ipts,ivm) * weight_ratio(ipts))
        !wind_damage_rate_closer(ipts,ivm) = max_damage_closer(ivm)
        !-Debug-
        IF (printlev_loc .GE. 4 .AND. ivm .EQ. test_pft .AND. ipts .EQ.test_grid) THEN 
            WRITE(numout,*) 'max_damage_closer:', max_damage_closer(ivm)
            WRITE(numout,*) 'sfactor_closer:', sfactor_closer(ivm)
            WRITE(numout,*) 'max_wind_speed_storm (m/s):', max_wind_speed_storm(ipts)
            WRITE(numout,*) 'max_wind_speed_pft_closer (m/s):', max_wind_speed_pft_closer(ipts,ivm)
            WRITE(numout,*) 'cws_break_closer',cws_break_closer_10(ipts,ivm)
            WRITE(numout,*) 'cws_overturn_closer',cws_overturn_closer_10(ipts,ivm)
            WRITE(numout,*) 'cws_closer (m/s)', cws_final_closer(ipts,ivm)
            WRITE(numout,*) 'wind_damage_rate:', wind_damage_rate_closer(ipts,ivm)
        ENDIF
  
        !! 4.2 Calculating critical wind speed in the forest and 10 m above the 
        !! zero-plane displacement for the further forest area
      
        ! Whichever of the two damage types (overturning, stem breakage) has the lowest 
        ! critical wind speed, it will be the occurring one.
        ! Critical wind speed for stem breakage in the inner forest area (m/s)
        distance_from_force = 0.0
        call calculate_wind_process(breaking_moment_further(ipts,ivm), stem_mass(ipts,ivm), &
                mean_dbh(ncirc), current_spacing(ipts,ivm), &
                canopy_depth_dom, canopy_breadth_dom, canopy_height_dom, &
                streamlining_n, streamlining_c, cws_break_further(ipts,ivm), &
                z0_wind_bk, zhd_wind_bk, distance_from_force, &
                gamma_solved_bk_further(ipts,ivm), ivm, ipts, DLF_breaking_further(ipts,ivm))
           
        ! Critical wind speed for stem breakage 10 m above the zero-plane 
        ! displacement in the inner forest area (m/s)
        call elevate(cws_break_further(ipts,ivm), z0_wind_bk, zhd_wind_bk, &
                canopy_height_dom, cws_break_further_10(ipts,ivm))
             
        ! Critical wind speed for tree overturning in the inner forest area (m/s)
        distance_from_force = 0.0
        call calculate_wind_process(overturning_moment_further(ipts,ivm), stem_mass(ipts,ivm), &
                mean_dbh(ncirc), current_spacing(ipts,ivm), &
                canopy_depth_dom, canopy_breadth_dom, canopy_height_dom, &
                streamlining_n, streamlining_c, cws_overturn_further(ipts,ivm), &
                z0_wind_ov, zhd_wind_ov, distance_from_force, &
                gamma_solved_ov_further(ipts,ivm),ivm, ipts, DLF_overturning_further(ipts,ivm))
   
        ! Critical wind speed for tree overturning 10 m above the zero-plane 
        ! displacement in the inner forest area (m/s)
        call elevate(cws_overturn_further(ipts,ivm), z0_wind_ov, &
                zhd_wind_ov, canopy_height_dom, cws_overturn_further_10(ipts,ivm))
       
        ! Whichever of the two damage types (overturning, stem breakage) has the 
        ! lowest critical wind speed, it will be the occurring one.
        IF (cws_overturn_further_10(ipts,ivm) < cws_break_further_10(ipts,ivm)) THEN
                cws_final_further(ipts,ivm) = cws_overturn_further_10(ipts,ivm)
            wind_damage_type_further(ipts,ivm) = ioverturning
            z0_wind(ipts,ivm) = z0_wind_ov
            zhd_wind(ipts,ivm) = zhd_wind_ov

        ELSE
            cws_final_further(ipts,ivm) = cws_break_further_10(ipts,ivm)
            wind_damage_type_further(ipts,ivm) = ibreakage
            z0_wind(ipts,ivm) = z0_wind_bk
            zhd_wind(ipts,ivm) = zhd_wind_bk
        END IF

        ! Guillaume M. -- Wind acclimation does not scale the critical wind speed here: it
        ! enters as the Gompertz weight (weight_ratio) on the further damage rate below.

        !max_wind_speed_pft_further(ipts,ivm) = max_wind_speed_storm(ipts) * &
        !    LOG((elevate_wind + zhd_wind(ipts,ivm) - zhd_pixel(ipts))/z0_wind(ipts,ivm))/ &
        !    LOG(elevate_wind/z0_wind(ipts,ivm))

        wind_profile_further(ipts,ivm) = &
            LOG((canopy_height_dom - canopy_height_dom*height_displacement)/z0_wind(ipts,ivm))/ &
            LOG((actual_maxheight - canopy_height_dom*height_displacement)/z0_wind(ipts,ivm))
        max_wind_speed_pft_further(ipts,ivm) = max_wind_speed_storm(ipts) * wind_profile_further(ipts,ivm)
        ! Guillaume M. -- The stored multiplier mirrors the wind_cal_pft choice made by
        ! the further logistic.
        IF (.NOT. wind_cal_pft) wind_profile_further(ipts,ivm) = un

        ! Wind damage level is considered as a logistic function with sigmoidal
        ! shape, which is calculated as function of a critical wind speed, an actual
        ! wind speed, a scaling factor(s)=0.8 and a maximum damage rate(80%)=0.8.
        IF(wind_cal_pft) THEN
           wind_damage_rate_further(ipts,ivm) = max_damage_further(ivm) * ( &
            (1./(1.+exp(-((max_wind_speed_pft_further(ipts,ivm) - &
            cws_final_further(ipts,ivm)) / sfactor_further(ivm))))) - &
            (1./(1.+exp(cws_final_further(ipts,ivm)/sfactor_further(ivm)))) )
        ELSE
           wind_damage_rate_further(ipts,ivm) = max_damage_further(ivm) * ( &
            (1./(1.+exp(-((max_wind_speed_storm(ipts)- cws_final_further(ipts,ivm)) / &
            sfactor_further(ivm))))) - &
            (1./(1.+exp(cws_final_further(ipts,ivm)/sfactor_further(ivm)))) ) 
        ENDIF
        ! Guillaume M. -- Apply the Gompertz wind-ratio weight to the further damage
        ! magnitude. The cap is max_damage_further, not max_damage_closer.
        wind_damage_rate_further(ipts,ivm) = MIN( &
            max_damage_further(ivm), wind_damage_rate_further(ipts,ivm) * weight_ratio(ipts))
        !wind_damage_rate_further(ipts,ivm) = max_damage_further(ivm)
        !-Debug-
        IF (printlev_loc .GE. 4 .AND. ivm .EQ. test_pft .AND. ipts .EQ.test_grid) THEN 
            WRITE(numout,*) 'max_damage_further:', max_damage_further(ivm)
            WRITE(numout,*) 'sfactor_further:', sfactor_further(ivm)
            WRITE(numout,*) 'max_wind_speed_storm (m/s):', max_wind_speed_storm(ipts)
            WRITE(numout,*) 'cws_break_further',cws_break_further_10(ipts,ivm)
            WRITE(numout,*) 'cws_overturn_further',cws_overturn_further_10(ipts,ivm)
            WRITE(numout,*) 'cws_further (m/s)', cws_final_further(ipts,ivm)
            WRITE(numout,*) 'max_wind_speed_pft_further:', max_wind_speed_pft_further(ipts,ivm)
            WRITE(numout,*) 'wind_damage_rate:', wind_damage_rate_further(ipts,ivm)  
        ENDIF
     
        !! 5. Determine which trees are killed by wind damage
        ! We only do the damage for the virtual_circ_class_n GT min_stomate 
        !+++JINA virtual or ccn?       
        IF (virtual_circ_class_n(ipts,ivm,ncirc) .GT. min_stomate) THEN
  
            !! 5.1. Kill trees in closer area  
            IF (wind_damage_type_closer(ipts,ivm) == ibreakage) THEN
                IF(printlev_loc .GE. 4 .AND. ivm .EQ. test_pft .AND. ipts .EQ. test_grid) THEN
                    WRITE(numout,*) 'Closer: kill_breakage'
                ENDIF
                CALL cal_storm_circ_class_kill( &
                    SUM(circ_class_biomass(ipts,ivm,:,:,icarbon),2), &
                    circ_class_n(ipts,ivm,:),wind_damage_rate_closer(ipts,ivm), &
                    closer_area_fraction(ipts), kill_break(ipts,ivm,:), ivm, ipts ) 

                ! We assume that there is Only one damage type will occur due to
                ! its lower limiting wind speed
                kill_uproot(ipts,ivm,:) = zero
            ELSE
                IF(printlev_loc .GE. 4 .AND. ivm .EQ. test_pft .AND. ipts .EQ. test_grid) THEN
                    WRITE(numout,*) 'Closer: kill_uproot'
                ENDIF
                CALL cal_storm_circ_class_kill( &
                    SUM(circ_class_biomass(ipts,ivm,:,:,icarbon),2), &
                    circ_class_n(ipts,ivm,:),wind_damage_rate_closer(ipts,ivm), &
                    closer_area_fraction(ipts), kill_uproot(ipts,ivm,:), ivm, ipts ) 
                ! We assume that there is Only one damage type will occur due to
                ! its lower limited wind speed
                kill_break(ipts,ivm,:) = zero
            ENDIF
                 
            !! 5.2 Kill trees in further area  
            IF (wind_damage_type_further(ipts,ivm) == ibreakage) THEN
                ! Total damage = closer + further
                IF(printlev_loc .GE. 4 .AND. ivm .EQ. test_pft .AND. ipts .EQ. test_grid) THEN
                    WRITE(numout,*) 'Further: kill_breakage'
                ENDIF
                CALL cal_storm_circ_class_kill( &
                    SUM(circ_class_biomass(ipts,ivm,:,:,icarbon),2), &
                    circ_class_n(ipts,ivm,:),wind_damage_rate_further(ipts,ivm), &
                    further_area_fraction(ipts), kill_break(ipts,ivm,:), ivm, ipts ) 
                ! We assume that there is only one damage type will occur due to
                ! its lower limited wind speed
                kill_uproot(ipts,ivm,:) = zero
            ELSE
                IF(printlev_loc .GE. 4 .AND. ivm .EQ. test_pft .AND. ipts .EQ. test_grid) THEN
                    WRITE(numout,*) 'Further: kill_uproot'
                ENDIF

                CALL cal_storm_circ_class_kill( &
                    SUM(circ_class_biomass(ipts,ivm,:,:,icarbon),2), &
                    circ_class_n(ipts,ivm,:),wind_damage_rate_further(ipts,ivm), &
                    further_area_fraction(ipts), kill_uproot(ipts,ivm,:), ivm, ipts ) 
                ! We assume that there is only one damage type will occur due to
                ! its lower limited wind speed
                kill_break(ipts,ivm,:) = zero
            ENDIF

        ELSE

            ! set no damage for the no vegetation area 
            kill_break(ipts,ivm,:) = zero
            kill_uproot(ipts,ivm,:) = zero
        END IF ! if ccn > min_stomate
              
        ! For the tree height below five meters (7.0m), we assumed that the damage will
        ! not occur. At the moment, we don't have good solution for the calculation of
        ! CWS for the small trees, so we just simply keep the small trees. 
        ! The 7-meter threshold was too rigid, and we need to consider PFT variabilities.
        ! Height threshold to truncate wind damage was determined by multiplying
        ! h_threshold_wind (0.3) to pipe_tune2 (maximum height when the diameter is 1m).
        ! 0.3 is quit arbitral decision based on simulations of MTC7 to avoid young
        ! trees keep dying from wind.
        IF ( mean_height(ncirc) .LT. pipe_tune2(ipts,ivm)* h_threshold_wind ) THEN
            kill_break(ipts,ivm,:) = zero
            kill_uproot(ipts,ivm,:) = zero
        ENDIF

        IF ( is_kill(ipts) == zero ) THEN
            kill_break(ipts,ivm,:) = zero
            kill_uproot(ipts,ivm,:) = zero
        ENDIF  

        ! Debug
        IF(printlev_loc .GE. 4 .AND. ivm .EQ. test_pft .AND. ipts .EQ. test_grid) THEN
            WRITE(numout,*) 'damage_further',wind_damage_rate_further(ipts,ivm)
            WRITE(numout,*) 'damage_closer',wind_damage_rate_closer(ipts,ivm)
            WRITE(numout,*) 'damage_type',wind_damage_type_further(ipts,ivm), &
                                               wind_damage_type_closer(ipts,ivm)
            WRITE(numout,*) 'ccn',circ_class_n(ipts,ivm,:)
            WRITE(numout,*) 'h_threshold',pipe_tune2(ipts,ivm)* h_threshold_wind
            WRITE(numout,*) 'height',mean_height(ncirc)
            WRITE(numout,*) 'kill_break',SUM(kill_break(ipts,ivm,:))
            WRITE(numout,*) 'kill_uproot',SUM(kill_uproot(ipts,ivm,:))
            WRITE(numout,*) 'is_kill',is_kill(ipts)
        ENDIF


        ! Save mortality to circ_class_cill
        circ_class_kill(ipts,ivm,:,forest_managed(ipts,ivm),icut_storm_break) = &
            kill_break(ipts,ivm,:)
        circ_class_kill(ipts,ivm,:,forest_managed(ipts,ivm),icut_storm_uproot) = &
            kill_uproot(ipts,ivm,:)
        DO icir = 1,ncirc
           wind_damage_all(ipts,ivm) = wind_damage_all(ipts,ivm) + &
                SUM(circ_class_biomass(ipts,ivm,icir,:,icarbon)) * &
                (kill_uproot(ipts,ivm,icir) + kill_break(ipts,ivm,icir))
           wood_damage_all(ipts,ivm) = wood_damage_all(ipts,ivm) + &
                (circ_class_biomass(ipts,ivm,icir,isapabove,icarbon) + &
                  circ_class_biomass(ipts,ivm,icir,iheartabove,icarbon) )*(un-branch_ratio(ivm))*  &
                (kill_uproot(ipts,ivm,icir) + kill_break(ipts,ivm,icir))

           IF (printlev_loc .GE. 4 .AND. ivm .EQ. test_pft .AND. ipts .EQ. test_grid) THEN
              WRITE(numout,*) 'icir',icir
              WRITE(numout,*) 'wind_damage_all',wind_damage_all(ipts,ivm)
              WRITE(numout,*) 'kill',kill_uproot(ipts,ivm,icir) + kill_break(ipts,ivm,icir)
              WRITE(numout,*) 'biomass',SUM(circ_class_biomass(ipts,ivm,icir,:,icarbon))
           ENDIF
        ENDDO

        wind_damage_wood_vol(ipts,ivm) = wood_damage_all(ipts,ivm)/ (pipe_density(ivm))
        wind_damage_rate(ipts,ivm) = wind_damage_all(ipts,ivm)/ &
                SUM(SUM(circ_class_biomass(ipts,ivm,:,:,icarbon),2))

        !-Debug-
        IF (printlev_loc .GE. 4 .AND. ipts==test_grid .AND. ivm==test_pft) THEN
            WRITE(numout,*) 'End of Section 4: Calculating critical wind speed and damage rate'
            WRITE(numout,*) 'critical wind speed further, ', cws_final_further(ipts,ivm)
            WRITE(numout,*) 'critical wind speed closer, ',  cws_final_closer(ipts,ivm)
            WRITE(numout,*) 'wind damage type further, ', wind_damage_type_further(ipts,ivm)
            WRITE(numout,*) 'wind damage type closer, ',  wind_damage_type_closer(ipts,ivm)
            WRITE(numout,*) 'forest managed type, ', forest_managed(ipts,ivm) 
            WRITE(numout,*) 'circ_class_n(ipts,ivm,:),',  circ_class_n(ipts,ivm,:) 
            WRITE(numout,*) 'circ_class_kill(ipts,ivm,:,:,icut_storm_break,', &
                circ_class_kill(ipts,ivm,:,forest_managed(ipts,ivm),icut_storm_break) 
            WRITE(numout,*) 'circ_class_kill(ipts,ivm,:,:,icut_storm_uproot,', &
                circ_class_kill(ipts,ivm,:,forest_managed(ipts,ivm),icut_storm_uproot) 
            WRITE(numout,*) 'circ_class_kill, for icut_1:9_:'
            WRITE(numout,*) circ_class_kill(ipts,ivm,:,forest_managed(ipts,ivm),1:9)
            WRITE(numout,*) 'wind_damage_rate',wind_damage_rate(ipts,ivm)
        ENDIF
      
        ! Copy some values for the output files
        canopy_depth_out(ipts,ivm) = canopy_depth_dom
        tree_h_edge_out(ipts,ivm) = tree_heights_from_edge
        mean_gap_f_out(ipts,ivm) = mean_gap_factor
        gap_size_out(ipts,ivm) = gap_size  
      

      END DO ! loop nvm
      
      IF (is_kill(ipts) == un) THEN
           max_wind_speed_storm(ipts) = zero
           max_wind_ratio_storm(ipts) = zero
      ENDIF

    END DO ! loop npts
   
    !! 5. Writing out history variables

    ! Write the CWS to stomate_history_file
    !CWS at further area  
    CALL histwrite_p(hist_id_stomate, 'CWS_FURTHER', itime, &
             cws_final_further(:,:), npts*nvm, horipft_index)
    !CWS at closer area 
    CALL histwrite_p(hist_id_stomate, 'CWS_CLOSER', itime, &
             cws_final_closer(:,:), npts*nvm, horipft_index)
    !CWS at closer area for stem breakage
    CALL histwrite_p(hist_id_stomate, 'CWS_CLOSER_BK', itime, &
             cws_break_closer_10(:,:), npts*nvm, horipft_index)
    !CWS at closer area for overturning
    CALL histwrite_p(hist_id_stomate, 'CWS_CLOSER_OV', itime, &
             cws_overturn_closer_10(:,:), npts*nvm, horipft_index)
    !CWS at further  area for stem break
    CALL histwrite_p(hist_id_stomate, 'CWS_FURTHER_BK', itime, &
             cws_break_further_10(:,:), npts*nvm, horipft_index)
    !CWS at further area for overturning
    CALL histwrite_p(hist_id_stomate, 'CWS_FURTHER_OV', itime, &
             cws_overturn_further_10(:,:), npts*nvm, horipft_index)
    !CWS at closer area for stem breakage
    CALL histwrite_p(hist_id_stomate, 'CWS_TOP_CLO_BK', itime, &
             cws_break_closer(:,:), npts*nvm, horipft_index)
    !CWS at closer area in each diameter class for overturning
    CALL histwrite_p(hist_id_stomate, 'CWS_TOP_CLO_OV', itime, &
             cws_overturn_closer(:,:), npts*nvm, horipft_index)
    !CWS at further  area for stem break
    CALL histwrite_p(hist_id_stomate, 'CWS_TOP_FUR_BK', itime, &
             cws_break_further(:,:), npts*nvm, horipft_index)
    !CWS at further area for overturning
    CALL histwrite_p(hist_id_stomate, 'CWS_TOP_FUR_OV', itime, &
             cws_overturn_further(:,:), npts*nvm, horipft_index)
    !D_H_ratio  for each diameter class
    CALL histwrite_p(hist_id_stomate, 'D_H_RATIO', itime, &
             d_h_ratio(:,:), npts*nvm, horipft_index)
    ! V_IND virtual individual tree numbers
    DO icir =1,ncirc        
        ! V_IND virtual individual tree numbers for each diameter class
        WRITE(var_name,'(A,I3.3)') 'V_IND_',icir
        CALL histwrite_p(hist_id_stomate, var_name, itime, &
             virtual_circ_class_n(:,:,icir), npts*nvm, horipft_index)
    ENDDO
    ! G_CLOSER gustiness for closer area 
    CALL histwrite_p(hist_id_stomate, 'G_CLOSER', itime, &
             gust_factor_closer(:,:), npts*nvm, horipft_index)
    ! G_FURTHER gustiness for further area 
    CALL histwrite_p(hist_id_stomate, 'G_FURTHER', itime, &
             gust_factor_further(:,:), npts*nvm, horipft_index)
    ! G_FURTHER gustiness for further area
    CALL histwrite_p(hist_id_stomate, 'STEM_MASS', itime, &
             stem_mass(:,:), npts*nvm, horipft_index)
    ! roughness length
    CALL histwrite_p(hist_id_stomate, 'WIND_Z0', itime, &
             z0_wind(:,:), npts*nvm, horipft_index)
    ! displacement height
    CALL histwrite_p(hist_id_stomate, 'WIND_D', itime, &
             zhd_wind(:,:), npts*nvm, horipft_index)
    ! canopy depth
    CALL histwrite_p(hist_id_stomate, 'CANOPY_DEPTH', itime, &
             canopy_depth_out(:,:), npts*nvm, horipft_index)
    ! tree heights from edge 
    CALL histwrite_p(hist_id_stomate, 'TREE_H_EDGE', itime, &
             tree_h_edge_out(:,:), npts*nvm, horipft_index)
    ! mean gap factor
    CALL histwrite_p(hist_id_stomate, 'MEAN_GAP_F', itime, &
             mean_gap_f_out(:,:), npts*nvm, horipft_index)
     ! gap size
     CALL histwrite_p(hist_id_stomate, 'GAP_SIZE', itime, &
          gap_size_out(:,:), npts*nvm, horipft_index)
    ! edge factor
        CALL histwrite_p(hist_id_stomate, 'EDGE_FACTOR', itime, &
             edge_factor_out(:,:), npts*nvm, horipft_index)
    ! maximum moment for overturning
    CALL histwrite_p(hist_id_stomate, 'MAX_OV_MOMENT', itime, &
             max_overturning_moment(:,:), npts*nvm, horipft_index)
    ! maximum moment breakage 
    CALL histwrite_p(hist_id_stomate, 'MAX_BK_MOMENT', itime, &
             max_breaking_moment(:,:), npts*nvm, horipft_index)
    ! Daily Maximum wind speed
    CALL histwrite_p(hist_id_stomate, 'DAILY_MAX_WIND', itime, &
          wind_speed_daily, npts, hori_index)
   
    CALL xios_orchidee_send_field("CWS_FURTHER", cws_final_further)
    CALL xios_orchidee_send_field("CWS_CLOSER", cws_final_closer)
    CALL xios_orchidee_send_field("CWS_CLOSER_BK", cws_break_closer_10)
    CALL xios_orchidee_send_field("CWS_CLOSER_OV", cws_overturn_closer_10)
    CALL xios_orchidee_send_field("CWS_FURTHER_BK", cws_break_further_10)
    CALL xios_orchidee_send_field("CWS_FURTHER_OV", cws_overturn_further_10)
    CALL xios_orchidee_send_field("GAP_AREA",gap_area_save(:,:,1))
    CALL xios_orchidee_send_field("DAILY_MAX_WIND_MOD", wind_speed_daily)
    CALL xios_orchidee_send_field("DAILY_MAX_WIND", wind_in)
    CALL xios_orchidee_send_field("WIND_SUM",SUM(wind_ratio_sum_save(:,:),2))
    CALL xios_orchidee_send_field("WIND_DAMAGE_RATE",wind_damage_rate)
    CALL xios_orchidee_send_field("IS_KILL_STORM", is_kill)
    CALL xios_orchidee_send_field("AREA_CLOSER_F", closer_area_fraction)
    CALL xios_orchidee_send_field("AREA_FURTHER_F", further_area_fraction)
    CALL xios_orchidee_send_field("GAMMA_FURTHER_BK",gamma_solved_bk_further)
    CALL xios_orchidee_send_field("GAMMA_FURTHER_OV",gamma_solved_ov_further)
    CALL xios_orchidee_send_field("GAMMA_CLOSER_BK",gamma_solved_bk_closer)
    CALL xios_orchidee_send_field("GAMMA_CLOSER_OV",gamma_solved_ov_closer)
    CALL xios_orchidee_send_field("CURRENT_SPACING",current_spacing)
    CALL xios_orchidee_send_field("ZHD_WIND",zhd_wind)
    CALL xios_orchidee_send_field("Z0_WIND",z0_wind)
    CALL xios_orchidee_send_field("MOMENT_CLOSER_BK",breaking_moment_closer)
    CALL xios_orchidee_send_field("MOMENT_CLOSER_OV",overturning_moment_closer)
    CALL xios_orchidee_send_field("MOMENT_FURTHER_BK",breaking_moment_further)
    CALL xios_orchidee_send_field("MOMENT_FURTHER_OV",overturning_moment_further)
    CALL xios_orchidee_send_field("V_IND", virtual_circ_class_n)
    CALL xios_orchidee_send_field("KILL_STORM_BK", kill_break)
    CALL xios_orchidee_send_field("KILL_STORM_OV", kill_uproot)
    CALL xios_orchidee_send_field("WIND_DAMAGE_BIOMASS", wind_damage_all)
    CALL xios_orchidee_send_field("WIND_DAMAGE_WOODVOL", wind_damage_wood_vol)
    CALL xios_orchidee_send_field("DLF_CLOSER_BK",DLF_breaking_closer)
    CALL xios_orchidee_send_field("DLF_CLOSER_OV",DLF_overturning_closer)
    CALL xios_orchidee_send_field("DLF_FURTHER_BK",DLF_breaking_further)
    CALL xios_orchidee_send_field("DLF_FURTHER_OV",DLF_overturning_further)
    CALL xios_orchidee_send_field("DAILY_MAX_WIND_PFT",max_wind_speed_pft_further)

    ! Guillaume M. -- Storm-perturbed forest area through the shared cc_kill_to_area
    ! basal-area route, same unit (m2) as AREA_DAMAGED_{FIRE,PEST,HARVEST}; GAP_AREA
    ! (gap geometry) is kept alongside. n_kill = stem breakage + uprooting, the
    ! forest-management dimension collapsed by SUM(.,DIM=3). Daily increment sent with
    ! XIOS accumulate: this is a gross area, never sum it across years.
    area_damaged_storm(:,:) = zero
    DO ivm = 1,nvm
       IF ( is_tree(ivm) ) THEN
          CALL cc_kill_to_area(npts, ivm,                                                   &
               SUM(circ_class_kill(:,ivm,:,:,icut_storm_break)                              &
                 + circ_class_kill(:,ivm,:,:,icut_storm_uproot), DIM=3),                    &
               circ_class_biomass(:,ivm,:,:,:), circ_class_n(:,ivm,:), veget_max(:,ivm),    &
               area(:)*contfrac(:), pipe_tune2(:,ivm), area_damaged_storm(:,ivm))
       ENDIF
    ENDDO
    CALL xios_orchidee_send_field("AREA_DAMAGED_STORM", area_damaged_storm)

    ! Guillaume M. -- Feed the storm-only damaged area to the edge-length accumulator,
    ! as dA_fire and dA_pest do. ALLOCATED(dA_storm) is equivalent to ok_aed_feedback:
    ! the accumulators exist only then. Filled here, before update_edge_length
    ! (stomate.f90) consumes and resets it. Do not feed gap_area_save instead: it counts
    ! all causes, harvest and background mortality included.
    IF (ALLOCATED(dA_storm)) THEN
       dA_storm(:) = dA_storm(:) + SUM(area_damaged_storm(:,:), DIM=2)
    ENDIF

   ! CALL xios_orchidee_send_field("KILL_STORM_OV", SUM(kill_uproot,DIM=3))

  END SUBROUTINE wind_damage


  !! =============================================================================================
  !! SUBROUTINE      : calculate_wind_process
  !!
  !>\BRIEF           : This subroutine calculates the breaking and overturning critical wind
  !!                   speeds
  !! 
  !! DESCRIPTION	   : None
  !!
  !! RECENT CHANGE(S): None
  !!
  !! RETURN VALUE    : 
  !!
  !! REFERENCES      : Gardiner et al. 2010, hale et al. 2012, Raupach 1994.
  !!
  !! FLOWCHART       : None
  !! 
  !_ =============================================================================================
     
  SUBROUTINE calculate_wind_process(critical_moment, stem_mass, mean_dbh, &
               current_spacing, canopy_depth_dom, canopy_breadth_dom, canopy_height_dom, &
               streamlining_n, streamlining_c, &
               guess_wind_speed, z0_wind, zhd_wind, distance_from_force, gamma_solved, ivm, & 
               ipts, DLF_calc)
  
  IMPLICIT NONE
  
  !! 0. Variable and parameter declaration
  
    !! 0.1 Input variables
    REAL, INTENT(in)  :: current_spacing
    REAL, INTENT(in)  :: canopy_depth_dom
    REAL, INTENT(in)  :: canopy_breadth_dom
    REAL, INTENT(in)  :: canopy_height_dom
    REAL, INTENT(in)  :: critical_moment
    REAL, INTENT(in)  :: streamlining_n
    REAL, INTENT(in)  :: streamlining_c
    REAL, INTENT(in)  :: distance_from_force
    REAL, INTENT(in)  :: stem_mass
    REAL, INTENT(in)  :: mean_dbh    
    INTEGER(i_std),INTENT(in)            :: ivm, ipts

    !! 0.2 Output variables
  
    REAL, INTENT(out)   :: guess_wind_speed         !! A guess for the critical wind speed
    REAL, INTENT(out)   :: z0_wind
    REAL, INTENT(out)   :: zhd_wind                       
    REAL, INTENT(inout) :: gamma_solved             !! the ratio between the critical wind speed at canopy
                                                    !!  top and the friction velocity (uh / u*)
    REAL, INTENT(inout) :: DLF_calc

    !! 0.3 Modified variables
  
    !! 0.4 Local variables
    INTEGER           :: loop_count                 !! Counting the iteration
    INTEGER           :: max_count=20               !! Max iteration times set to 20 iterations              
    REAL              :: delta                      !! Difference between the guessed
                                                    !! and the calculated wind speed
    REAL              :: wind_precision             !! How close we should
                                                    !! get with the iterations
    REAL              :: lambda                     !!
    REAL              :: lambda_capital 
    REAL              :: psih
    REAL              :: critical_moment_cws

    REAL, PARAMETER   :: air_density = 1.2250
    REAL, PARAMETER   :: f_crown_mass = 1.136
    REAL, PARAMETER   :: f_knots = 1.000
      
    !! 0.5 Variables used in  calculate_force subroutine and  calculate_bending_moment subroutine                  
    !! (Both subroutines are not used in the current code)
    REAL              :: gammaSolved
    REAL              :: d
    REAL              :: force_on_tree
    REAL              :: height_of_force
    REAL              :: bending_moment  
    REAL              :: porosity
    REAL              :: crit_wind_speed
      
!_ =============================================================================================
  
  !! 1. Subroutine
  
    ! This bit belongs to the original form of doing the iterations. For more info, see the notes below.
    !  guess_wind_speed = 64.0
    !  delta = guess_wind_speed / 2.0
      
    guess_wind_speed = 25.0
    delta = guess_wind_speed / 2.0 
    wind_precision = 0.01
    loop_count = 0
    DLF_calc = 1.136

    IF(cal_deflection) THEN
       CALL dlf_fun_orc(DLF_calc, critical_moment, canopy_height_dom, canopy_depth_dom, &
                         canopy_breadth_dom, stem_mass, mean_dbh)
       DLF_calc = MIN(MAX(0.5,DLF_calc),2.5)
    ENDIF

    DO WHILE (delta > wind_precision .OR. &
         zhd_wind - distance_from_force .LT. min_stomate)
        CALL streamlining(guess_wind_speed, current_spacing, &
              canopy_breadth_dom, canopy_depth_dom, canopy_height_dom, & 
              streamlining_n, streamlining_c, gamma_solved, zhd_wind, z0_wind)


        ! Calculate critical wind speed to compare with wind speeds used for calculating
        ! streamlining in the iterative calculation of the final critical wind speed)
        IF (zhd_wind - distance_from_force .LT. min_stomate) THEN
           ! Impossible to calculate the critical wind speed
           crit_wind_speed = 100
        ELSE
           !crit_wind_speed = (gamma_solved / current_spacing) * (critical_moment / &
           !     (f_crown_mass * air_density * DLF_calc * (zhd_wind - distance_from_force) ))**0.5 
           crit_wind_speed = (gamma_solved / current_spacing) * (critical_moment / &
                (air_density * DLF_calc * (zhd_wind - distance_from_force) ))**0.5

        ENDIF

        ! We use the absolute error as the criteria for stopping the iteration
        ! loop    
        delta = abs( guess_wind_speed - crit_wind_speed ) 
         
        guess_wind_speed = crit_wind_speed
         
        loop_count = loop_count  +1      
        ! Add the maximum iteration times condition for CWS calculation
        IF (loop_count .GE. max_count) THEN
            WRITE(numout,*) 'WARANNING from stomate_windthrow.f90: The maximum iteration times for CWS claculation!!'  
            EXIT
        ENDIF
    END DO
     
  END SUBROUTINE calculate_wind_process
  

  !! ===========================================================================================
  !! SUBROUTINE      : streamlining
  !!
  !>\BRIEF           : This subroutine calculates the streamlining of tree crowns. Streamlining is
  !!                   the change of shape of the crowns due to wind.
  !! 
  !! DESCRIPTION	   : None
  !!
  !! RECENT CHANGE(S): None
  !!
  !! RETURN VALUE    : Surface roughness and zero-plane displacement among other variables that
  !!                   are used in the other subroutines
  !! REFERENCES      :
  !!                   Raupach (1994) Simplified expressions for vegetation roughness length and
  !!                   zero-plane displacement as functions of canopy height and area index, Boundary-Layer Meteorology 
  !!                   October 1994, Volume 71, Issue 1, pp 211-216, doi: 10.1007/BF00709229
  !!
  !! FLOWCHART       : None
  !! 
  !_ ============================================================================================
     
  SUBROUTINE streamlining(guess_wind_speed, current_spacing, &
               canopy_breadth_dom, canopy_depth_dom, canopy_height_dom, &
               streamlining_n, streamlining_c, gamma_solved, zhd_wind, z0_wind)
  
  IMPLICIT NONE

  !! 0. Variable and parameter declaration
  
    !! 0.1 Input variables
  
    REAL,  INTENT(in)  :: guess_wind_speed  !! the guess wind speed (m s{-1})
    REAL,  INTENT(in)  :: current_spacing   !! the average spacing between trees (m)
    REAL,  INTENT(in)  :: canopy_breadth_dom!! the maximum width of the canopy   (m)
    REAL,  INTENT(in)  :: canopy_depth_dom  !! the length of the live tree crown (m)
    REAL,  INTENT(in)  :: canopy_height_dom !! the dominant canopy height (m) 
    REAL,  INTENT(in)  :: streamlining_n    !! streamlining parameter n (-) see Eq. A12 in Hale et al 2015
    REAL,  INTENT(in)  :: streamlining_c    !! streamlining parameter c (-) see Eq. A12 in Hale et al 2015  
 
    !! 0.2 Output variables
  
    REAL, INTENT(out)  :: gamma_solved      !! a function for solving surface roughness (-)   
    REAL, INTENT(out)  :: zhd_wind          !! zero plane of displacement height (m)
    REAL, INTENT(out)  :: z0_wind           !! surface roughness (m)
  
    !! 0.3 Modified variables
  
    !! 0.4 Local variables

    REAL               :: lambda            !! roughness density, the frontal
                                            !! area of roughness elements
    REAL               :: lambda_capital    !! canopy area index 
    REAL               :: psih              !! stability correction function for heat transfer (-)  
    REAL               :: porosity          !! the effective drag ?? (-)
                                            !! see equation A12 in the supplementary material of
                                            !! Hale et al (2015)

    !! The aerodynamic parameter values used here are based on the reference of Raupach (1994)
    !! --- TEMP ----
    !! Need to be discussed for the consistency in the ORCHIDEE
    !! ------------         
    REAL, PARAMETER :: c_displacement = 7.5   
    REAL, PARAMETER :: c_surface = 0.003             
    REAL, PARAMETER :: c_drag = 0.3               
    REAL, PARAMETER :: c_roughness = 2.0      
    REAL, PARAMETER :: lambda_capital_max = 0.6
    REAL, PARAMETER :: wind_streamlining_max = 25.0 
    REAL, PARAMETER :: wind_streamlining_min = 10.0
   !_ =============================================================================================
  
    !! 1. Subroutine
  
    IF (guess_wind_speed > wind_streamlining_max) THEN
        ! Above 25 m/s, tree crowns cannot streamline any further.
        ! The porosity set to its minimum value 2.35*25^-0.51 ~ 0.46
        porosity = streamlining_c * wind_streamlining_max**(-streamlining_n)            
    ELSE IF ( guess_wind_speed < wind_streamlining_min ) THEN
        ! Below 10 m/s, tree crowns do not streamline and porosity set to its
        ! maximum value 2.35*10^-0.51 ~ 0.73
        porosity = streamlining_c * wind_streamlining_min**(-streamlining_n)            
    ELSE
        ! Wind speed between 10 m/s and 25 m/s,  tree crown do the streamlining 
        porosity = streamlining_c * guess_wind_speed**(-streamlining_n)
    END IF
 
    ! lambda = (canopy_breadth/2.0) * canopy_depth *  porosity / current_spacing **2.0
    ! Here, we directly calculated the mean width of canopy from the value of "canopy_breadth"  
    ! The shape of tree crown is not rhomboid anymore. In the ORCHIDEE-CAN
    ! the shape of tree crown is circle, so the frontal area should be
    ! calculated as (pi/4 * canopy_breadth * canopy_depth). 
    lambda = canopy_breadth_dom * canopy_depth_dom *  porosity / current_spacing **2.0
    lambda = 2*(pi/4 * canopy_breadth_dom * canopy_depth_dom *  porosity / current_spacing**2.0)
      
    lambda_capital = lambda
      
    IF (lambda_capital > lambda_capital_max) THEN
        ! Gamma is needed for the calculation of the roughness length.
        gamma_solved = 1.0 / ((c_surface + c_drag * lambda_capital_max / 2.0)**0.5)
    ELSE
        gamma_solved = 1.0 / ((c_surface + c_drag * lambda_capital / 2.0)**0.5)
    END IF
      
    ! Calculation of the zero-plane displacement: 
    ! see equations from A7 to A11 in the supplementary material of Hale et al. (2015)  
    zhd_wind = canopy_height_dom * (1.0 - ((1.0 - EXP(-(c_displacement * &
        lambda_capital)**0.5)) / (c_displacement * lambda_capital)**0.5))
      
    ! Calculation of the aerodynamic roughness length:
    ! psih = "Profile influence function" at h (height of the roughness elements).
    psih = LOG(c_roughness) -1.0 + c_roughness**(-1.0)
      
    z0_wind = (canopy_height_dom - zhd_wind) * EXP(psih - ct_karman * gamma_solved)
      
  END SUBROUTINE streamlining


  !!
  !============================================================================================
  !! SUBROUTINE      : calculate_dlf
  !!
  !>\BRIEF           : This subroutine calculates the deflection loading factor
  !! 
  !! RETURN VALUE    : DLF_calc
  !!
  !! REFERENCES      : Gardiner 1989 'Mechanical characteristics of Sitka
  !                         Spruce' Eq. 5: deflection at distance x from top of tree. 
  !                    https://www.forestresearch.gov.uk/tools-and-resources/fthr/forestgales/
  !                    fgr R package
  !_
  !============================================================================================

  SUBROUTINE dlf_fun_orc(DLF_calc, critical_moment, canopy_height_dom, canopy_depth_dom, canopy_breadth_dom, &
                         stem_mass, mean_dbh)

  IMPLICIT NONE

  !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    REAL,  INTENT(in)  :: critical_moment    
    REAL,  INTENT(in)  :: canopy_height_dom    
    REAL,  INTENT(in)  :: canopy_depth_dom    
    REAL,  INTENT(in)  :: canopy_breadth_dom   
    REAL,  INTENT(in)  :: stem_mass    
    REAL,  INTENT(in)  :: mean_dbh    


    !! 0.2 Output variables
    REAL, INTENT(out)  :: DLF_calc               !! calculated deflection loading factor

    !! 0.3 Modified variables

    !! 0.4 Local variables
    REAL               :: force_of_wind          
    REAL               :: deflection1            !! Displacement at mid-crwon height
    REAL               :: deflection2            !! Deflection ¾ of the way down the stem.  
                                                 !! Assumption that centre of stem mass is at this height
    
    !! --- TEMP ----
    !! Need to be discussed for the consistency in the ORCHIDEE
    !! ------------         
    REAL, PARAMETER :: crown_density = 2.5
    REAL, PARAMETER :: snow_density = 150
    REAL, PARAMETER :: snow_depth = 0
   !_
   !=============================================================================================

    force_of_wind = critical_moment / (canopy_height_dom - canopy_depth_dom/2)
    CALL deflection_fun(deflection1, canopy_depth_dom/2, canopy_height_dom, force_of_wind, &
         mean_dbh, canopy_height_dom, canopy_depth_dom, (canopy_height_dom - canopy_depth_dom/2))
    CALL deflection_fun(deflection2, (0.75 * canopy_height_dom), canopy_height_dom, force_of_wind, &
         mean_dbh, canopy_height_dom, canopy_depth_dom, (canopy_height_dom - canopy_depth_dom/2))

    DLF_calc = 1 / (1 - ((deflection1 * &
      ! weight of crown and snow
      ((((pi * canopy_depth_dom * (canopy_breadth_dom/2)**2)/3) * crown_density ) + &
      (pi * (canopy_breadth_dom/2)**2 * snow_depth *snow_density)) * cte_grav) + &
      ! weight of stem  
      deflection2 * stem_mass * cte_grav) / critical_moment)
    !Note that moment is in both numerator and denominator, so the magnitude of
    !the moment does not have any impact on DLF. It was decided to be kept in
    !the equation for completeness.

  END SUBROUTINE

  !!
  !============================================================================================
  !! SUBROUTINE      : deflection_fun
  !!
  !>\BRIEF           : This subroutine calculates the deflection loading factor
  !! 
  !! RETURN VALUE    : deflection
  !!
  !! REFERENCES      : Gardiner 1989 'Mechanical characteristics of Sitka
  !                         Spruce' Eq. 5: deflection at distance x from top of
  !                         tree. 
  !                    https://www.forestresearch.gov.uk/tools-and-resources/fthr/forestgales/
  !                    fgr R package
  !_
  !============================================================================================


  SUBROUTINE deflection_fun(deflection, x, lever_arm, fow, dbh, ht, cr_depth, pull_height) 

  IMPLICIT NONE

  !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    REAL,  INTENT(in)  :: x                  !! Height along the tree stem where deflection is to be calculated (m)
    REAL,  INTENT(in)  :: lever_arm          !! Length of the lever arm (m).
    REAL,  INTENT(in)  :: fow                !! Applied wind loading in the function.
    REAL,  INTENT(in)  :: dbh                !! Diameter of the stem (m)
    REAL,  INTENT(in)  :: ht                 !! Tree height (m)
    REAL,  INTENT(in)  :: cr_depth           !! Length of the tree crown (m)
    REAL,  INTENT(in)  :: pull_height        !! Height along the tree stem where the wind loading is applied (m)

    !! 0.2 Output variables
    REAL, INTENT(out)  :: deflection

    !! 0.3 Modified variables

    !! 0.4 Local variables
    REAL               :: r_value            !! Ratio of distance from top of tree where force is applied to tree height
    REAL               :: i_value            !! Second area moment of Inertia calculated at tree base
    REAL               :: hh                 !! Constant

    !! Need to be discussed for the consistency in the ORCHIDEE
    !! ------------         
    REAL, PARAMETER :: snow_depth = 0
    REAL, PARAMETER :: an = -1.4
    REAL, PARAMETER :: bn = -0.4
    REAL, PARAMETER :: cn = 0.6
    REAL, PARAMETER :: moe = 7.6e+09
   !_
   !=============================================================================================

    hh = an * bn * cn
    i_value = pi * (dbh**4 / 64)
    r_value = (lever_arm - pull_height) / lever_arm

    deflection = ((fow*lever_arm ** 3) / (moe * i_value * hh)) * (an * (x/lever_arm) ** cn - &
           r_value * cn*(x/lever_arm) ** bn + &
           r_value * bn*cn*(x/lever_arm) - an * cn * (x/lever_arm) + an * bn - &
           r_value * an * cn)

  END SUBROUTINE

  !!
  !=============================================================================================
  !! SUBROUTINE      : elevate
  !!
  !>\BRIEF           : This subroutine converts the critical wind speeds
  !calculated above to the
  !!                   wind speed 10 meters above the zero-plane displacement
  !! 
  !! DESCRIPTION           : None
  !!
  !! RECENT CHANGE(S): None
  !!
  !! RETURN VALUE    : 
  !!
  !! REFERENCES      :
  !!
  !! FLOWCHART       : None
  !!  
  !_
  !=============================================================================================

  SUBROUTINE elevate(uh_speed, z0_wind, zhd_wind, mean_height,elevate_result)

    IMPLICIT NONE

  !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    REAL,  INTENT(in)  :: uh_speed
    REAL,  INTENT(in)  :: z0_wind
    REAL,  INTENT(in)  :: zhd_wind
    REAL,  INTENT(in)  :: mean_height

    !! 0.2 Output variables

    REAL, INTENT(out)  :: elevate_result

    !! 0.3 Modified variables

    !! 0.4 Local variables

  !_
  !==============================================================================================

    !! 1. Subroutine
    elevate_result = uh_speed * log(elevate_wind / z0_wind) / &
           log((mean_height - zhd_wind) / z0_wind)

  END SUBROUTINE elevate


  !!
  !=============================================================================================
  !! SUBROUTINE      : cal_storm_circ_class_kill
  !!
  !>\BRIEF           : This subroutine calculates number of trees per circ class based on damage rate
  !! 
  !! DESCRIPTION           : None
  !!
  !! RECENT CHANGE(S): None
  !!
  !! RETURN VALUE    : 
  !!
  !! REFERENCES      :
  !!
  !! FLOWCHART       : None
  !!  
  !_
  !=============================================================================================

  SUBROUTINE cal_storm_circ_class_kill(biomass_temp, circ_class_n_temp, damage_rate, &
                area_frac, kill_temp, ivm, ipts)

        IMPLICIT NONE
  !! 0. Variable and parameter declaration
        !! 0.1 Input variables
	REAL(r_std),DIMENSION(:),INTENT(in)  :: biomass_temp        !! total carbon mass per diameter class
	REAL(r_std),DIMENSION(:),INTENT(in)  :: circ_class_n_temp   !! circ_class_n for pixel, pft
	REAL(r_std),INTENT(in)               :: damage_rate         !! calculated storm damage rate 
	REAL(r_std),INTENT(in)               :: area_frac           !! fraction of area for closer or further to the edge
        INTEGER(i_std),INTENT(in)            :: ivm, ipts
        !! 0.2 Output variables
        
        !! 0.3 Modified variables
	REAL(r_std),DIMENSION(:),INTENT(out) :: kill_temp           !! storm mortality for each diameter class
        
        !! 0.4 Local variables
	REAL(r_std)                          :: biomass_to_kill     !! target biomass to remove
	REAL(r_std)                          :: living_biomass      !! living carbon mass per closer/further area
	INTEGER(i_std)                       :: icir
 
  !==============================================================================================
  
    ! Now the damage rate is based on the total biomass and number of trees
    ! Biomass is removed from the largest circ class to smallest one.
    biomass_to_kill = SUM(biomass_temp * circ_class_n_temp) * damage_rate * area_frac
    
    IF(printlev_loc .GE. 4 .AND. ivm .EQ. test_pft .AND. ipts .EQ. test_grid) THEN
        WRITE(numout,*) 'biomass_temp',biomass_temp
        WRITE(numout,*) 'ccn_temp',circ_class_n_temp
        WRITE(numout,*) 'damage_rate',damage_rate
        WRITE(numout,*) 'biomass_to_kill',biomass_to_kill
        WRITE(numout,*) 'kill_temp',kill_temp
    ENDIF
    
    DO icir=ncirc,1,-1

        living_biomass= &
            biomass_temp(icir) * circ_class_n_temp(icir) * area_frac

        IF(printlev_loc .GE. 4 .AND. ivm .EQ. test_pft .AND. ipts .EQ. test_grid) THEN
            WRITE(numout,*) 'icir',icir
            WRITE(numout,*) 'living_biomass',living_biomass
            WRITE(numout,*) 'area_frac',area_frac
        ENDIF

        IF(living_biomass .EQ. zero) THEN
           CYCLE
        ENDIF

        IF(living_biomass .LE. biomass_to_kill)THEN
            ! biomass to kill exceeds the living mass. Remove all by the fraction of area.
            biomass_to_kill = biomass_to_kill - living_biomass
            kill_temp(icir) = kill_temp(icir) + circ_class_n_temp(icir)

        ELSE
            ! biomass to kill is less than the living mass. Remove number of individuals
            ! proportional to the biomass
            kill_temp(icir)= kill_temp(icir) + biomass_to_kill/living_biomass * &
                circ_class_n_temp(icir)
            biomass_to_kill = zero

        ENDIF ! living_biomass .LE. biomass_to_kill
        IF(printlev_loc .GE. 4 .AND. ivm .EQ. test_pft .AND. ipts .EQ. test_grid) THEN
            WRITE(numout,*) 'after calculation'
            WRITE(numout,*) 'kill_temp',kill_temp(icir)
            WRITE(numout,*) 'biomass_to_kill',biomass_to_kill
        ENDIF
    ENDDO ! loop over circ classes

    WHERE(kill_temp(:) < min_stomate)
        kill_temp(:) = zero
    ENDWHERE

  END SUBROUTINE 

  !==============================================================================================
  ! Below subroutines are not used for now
  !==============================================================================================
   
  !! ============================================================================================
  !! SUBROUTINE      : calculate_force
  !!
  !>\BRIEF           : This subroutine calculates the force affecting the tree.
  !! 
  !! DESCRIPTION	   : None
  !!
  !! RECENT CHANGE(S): None
  !!
  !! RETURN VALUE    : 
  !!
  !! REFERENCES      :
  !!
  !! FLOWCHART       : None
  !! 
  !_ ============================================================================================
     
    SUBROUTINE calculate_force(guess_wind_speed, a_gamma_solved, &
                               a_zhd_wind, current_spacing, a_mean_height, & 
                               force_on_tree, height_of_force)
  
    IMPLICIT NONE
  
  !! 0. Variable and parameter declaration
  
      !! 0.1 Input variables
  
      REAL, INTENT(in)  :: guess_wind_speed
      REAL, INTENT(in)  :: a_gamma_solved
      REAL, INTENT(in)  :: a_zhd_wind
      REAL, INTENT(in)  :: current_spacing
      REAL, INTENT(in)  :: a_mean_height
      !! 0.2 Output variables
  
      REAL, INTENT(out) :: force_on_tree
      REAL, INTENT(out) :: height_of_force
  
      !! 0.3 Modified variables
  
      !! 0.4 Local variables
      REAL,PARAMETER    :: air_density = 1.2226 !!This should call from constant.f90 (kg m^{-3}) 
  !_ =============================================================================================
  
      !! 1. Subroutine
      
      force_on_tree = air_density * (guess_wind_speed * &
           current_spacing / a_gamma_solved)**2.0

      !! --- TEMP ---
      !!   I checked with the reference the force of wind can be calculated as 
      !!   air_density*ustar^{2}*D^{2} as the gamma is defined as  U/ustar , D is the current spacing
      !!   the force of wind should be as 
      !!   force_on_tree = air_density * (U/gamma)^{2} * D^{2}
      !! ------------      
      height_of_force = a_zhd_wind
      
      ! The maximum of the zero-plane displacement height is locked at 0.8 tree height.
      ! IF (height_of_force > 0.8 * a_mean_height) THEN
      !    height_of_force = 0.8 * a_mean_height
      ! END IF
      
    END subroutine calculate_force
    
  
  !! =============================================================================================
  !! SUBROUTINE      : calculate_bending_moment
  !!
  !>\BRIEF           : Subroutine to calculate the bending moment
  !!                   (to compare with critical moments in the iterative calculation of the critical
  !!                   wind speeds)
  !! 
  !! DESCRIPTION	   : None
  !!
  !! RECENT CHANGE(S): None
  !!
  !! RETURN VALUE    : 
  !!
  !! REFERENCES      :
  !!
  !! FLOWCHART       : None
  !! 
  !_ =============================================================================================
     
    SUBROUTINE calculate_bending_moment (a_force_on_tree, a_height_of_force, bending_moment)
  
    IMPLICIT NONE
  
  !! 0. Variable and parameter declaration
  
      !! 0.1 Input variables
  
      REAL, INTENT(in)  :: a_force_on_tree 
      REAL, INTENT(in)  :: a_height_of_force
  
      !! 0.2 Output variables
  
      REAL, INTENT(out) :: bending_moment
  
      !! 0.3 Modified variables
   
      !! 0.4 Local variables
      REAL              :: f_crown_mass 

  !_ =============================================================================================
      !! 1. Subroutine
  
      bending_moment = a_force_on_tree * a_height_of_force * f_crown_mass
      
    END SUBROUTINE calculate_bending_moment
    

END MODULE stomate_windthrow
