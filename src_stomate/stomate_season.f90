! =================================================================================================================================
! MODULE        : stomate_season
!
! CONTACT       : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE       : IPSL (2006). This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       This module calculates long-term meteorological parameters from daily temperatures
!! and precipitations (essentially for phenology). 
!!      
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_stomate/stomate_season.f90 $ 
!! $Date: 2026-02-20 13:02:42 +0100 (ven. 20 févr. 2026) $
!! $Revision: 9386 $
!! \n
!_ =================================================================================================================================

MODULE stomate_season

  ! modules used:
  USE xios_orchidee
  USE ioipsl_para
  USE stomate_data
  USE constantes
  USE constantes_soil
  USE constantes_mtc
  USE pft_parameters
  USE dynamic_parameters
  USE solar
  USE grid
  USE function_library,  ONLY: lai_to_biomass,biomass_to_lai, &
                               wood_to_ba, get_printlev
  USE time,              ONLY : one_year, julian_diff, ts_annual_proc

  IMPLICIT NONE

  ! private & public routines

  PRIVATE
  PUBLIC season_pre_disturbance, season_post_disturbance, &
       season_pre_disturbance_clear, season_post_disturbance_clear

  LOGICAL, SAVE              :: firstcall_season_pre_disturbance = .TRUE.  !! first call (true/false)
!$OMP THREADPRIVATE(firstcall_season_pre_disturbance)
  LOGICAL, SAVE              :: firstcall_season_post_disturbance = .TRUE.  !! first call (true/false)
!$OMP THREADPRIVATE(firstcall_season_post_disturbance)

  INTEGER(i_std), SAVE       :: printlev_loc                                !! Local level of text output for current module
!$OMP THREADPRIVATE(printlev_loc)

CONTAINS

!! ================================================================================================================================
!! SUBROUTINE   : season_clear_post_disturbance
!!
!>\BRIEF          Flag setting 
!!
!! DESCRIPTION  : This subroutine sets flags ::firstcall_season, to .TRUE., and therefore activates   
!!                section 1.1 of the ::season subroutine which writes messages to the output. \n
!!                This subroutine is called at the end of the subroutine ::stomate_clear, in the 
!!                module ::stomate.
!!
!! RECENT CHANGE(S):None
!!
!! MAIN OUTPUT VARIABLE(S): ::firstcall_season_pre_disturbance
!!
!! REFERENCE(S)  : None 
!!
!! FLOWCHART     : None
!! \n             
!_ =================================================================================================================================

  SUBROUTINE season_pre_disturbance_clear
    firstcall_season_pre_disturbance =.TRUE.
  END SUBROUTINE season_pre_disturbance_clear


!! ================================================================================================================================
!! SUBROUTINE   : season_clear_post_disturbance
!!
!>\BRIEF          Flag setting 
!!
!! DESCRIPTION  : This subroutine sets flags ::firstcall_season, to .TRUE., and therefore activates   
!!                section 1.1 of the ::season subroutine which writes messages to the output. \n
!!                This subroutine is called at the end of the subroutine ::stomate_clear, in the 
!!                module ::stomate.
!!
!! RECENT CHANGE(S):None
!!
!! MAIN OUTPUT VARIABLE(S): ::firstcall_season_post_disturbance
!!
!! REFERENCE(S)  : None 
!!
!! FLOWCHART     : None
!! \n             
!_ =================================================================================================================================

  SUBROUTINE season_post_disturbance_clear
    firstcall_season_post_disturbance =.TRUE.
  END SUBROUTINE season_post_disturbance_clear


!! ================================================================================================================================
!! SUBROUTINE   : season_pre_disturbance
!!
!>\BRIEF          This subroutine calculates many of the long-term biometeorological variables
!!                needed in the phenology subroutines and in the calculation of long-term vegetation
!!                dynamics in the LPJ DGVM. 
!!
!! DESCRIPTION    This subroutine is called by the module ::stomate before LPJ, and mainly deals  
!!                with the calculation of long-term meteorological variables, carbon fluxes and 
!!                vegetation-related variables that are used to calculate vegetation dynamics 
!!                in the stomate modules relating to the phenology and to the longer-term 
!!                changes in vegetation type and fractional cover in the LPJ DGVM modules. \n
!!                In sections 2 to 5, longer-term meteorological variables are calculated. The
!!                long-term moisture availabilities are used in the leaf onset and senescence
!!                phenological models ::stomate_phenology and ::stomate_turnover that require
!!                a moisture condition. The long term temperatures are also required for phenology
!!                but in addition they are used in calculations of C flux and the presence and
!!                establishment of vegetation patterns on a longer timescale in the LPJ DGVM modules.
!!                Finally the monthly soil humidity/relative soil moisture is used in the C
!!                allocation module ::stomate_alloc. \n
!!                Sections 12 to 14 also calculate long-term variables of C fluxes, including NPP,
!!                turnover and GPP. The first two are used to calculate long-term vegetation 
!!                dynamics and land cover change in the LPJ DVGM modules. The weekly GPP is used
!!                to determine the dormancy onset and time-length, as described below. \n
!!                The long-term variables described above are are used in certain vegetation
!!                dynamics processes in order to maintain consistency with the basic hypotheses of the 
!!                parameterisations of LPJ, which operates on a one year time step (Krinner et al., 2005).
!!                In order to reduce the computer memory requirements, short-term variables (e.g. daily
!!                temperatures) are not stored for averaging over a longer period. Instead
!!                the long-term variables (e.g. monthly temperature) are calculated at each time step
!!                using a linear relaxation method, following the equation:
!!                \latexonly
!!                \input{season_lin_relax_eqn1.tex}
!!                \endlatexonly
!!                \n 
!!                The long-term variables are therefore updated on a daily basis. This method allows for 
!!                smooth temporal variation of the long-term variables which are used to calculate 
!!                vegetation dynamics (Krinner et al., 2005). \n
!!                Sections 6 to 11 calculate the variables required for to determine leaf onset in
!!                the module ::stomate_phenology. 
!!                These include :
!!                - the dormance onset and time-length, when GPP is below a certain threshold \n
!!                - the growing degree days (GDD), which is the sum of daily temperatures greater than
!!                -5 degrees C, either since the onset of the dormancy period, or since midwinter \n
!!                - the number of chilling days, which is the number of days with a daily temperature
!!                lower than a PFT-dependent threshold since the beginning of the dormancy period \n
!!                - the number of growing days, which is the the number of days with a temperature
!!                greater than -5 degrees C since the onset of the dormancy period \n
!!                - the time since the minimum moisture availability in the dormancy period. \n
!!                These variables are used to determine the conditions needed for the start of the
!!                leaf growing season. The specific models to which they correspond are given below. \n
!!                Sections 15 to 20 are used to update the maximum/minimum or the sum of various 
!!                meteorological and biological variables during the year that are required for
!!                calculating either leaf onset or longer-term vegetation dynamics the following year. \n
!!                At the end of the year, these variables are updated from "thisyear" to "lastyear", 
!!                in Section 21 of this subroutine, for use the following year. \n
!!                Finally the probably amount of herbivore consumption is calculated in Section 22,
!!                following McNaughton et al. (1989).
!!
!! RECENT CHANGE(S): None
!!                
!! MAIN OUTPUT VARIABLE(S): :: herbivores,
!!                        :: maxvegstress_lastyear, :: maxvegstress_thisyear, 
!!                        :: minvegstress_lastyear, :: minvegstress_thisyear, 
!!                        :: maxgppweek_lastyear, :: maxgppweek_thisyear, 
!!                        :: gdd0_lastyear, :: gdd0_thisyear, 
!!                        :: precip_lastyear, :: precip_thisyear, 
!!                        :: lm_lastyearmax, :: lm_thisyearmax,  
!!                        :: maxfpc_lastyear, :: maxfpc_thisyear, 
!!                        :: vegstress_month, :: vegstress_week, 
!!                        :: t2m_longterm, :: t2m_month, :: t2m_week, 
!!                        :: tsoil_month, ::  
!!                        :: npp_longterm, :: turnover_longterm, :: gpp_week, 
!!                        :: gdd_m5_dormance, :: gdd_midwinter, 
!!                        :: ncd_dormance, :: ngd_minus5, :: time_lowgpp, 
!!                        :: time_hum_min, :: hum_min_dormance 
!!      
!! REFERENCES   : 
!! - Krinner, G., N. Viovy, N. de Noblet-Ducoudre, J. Ogee, J. Polcher, P. 
!! Friedlingstein, P. Ciais, S. Sitch and I.C. Prentice (2005), A dynamic global
!! vegetation model for studies of the coupled atmosphere-biosphere system, Global
!! Biogeochemical Cycles, 19, doi:10.1029/2003GB002199.
!! - McNaughton, S.J., M. Oesterheld, D.A. Frank and K.J. Williams (1989), 
!! Ecosystem-level patterns of primary productivity and herbivory in terrestrial
!! habitats, Nature, 341, 142-144.
!!                       
!! FLOWCHART    : 
!! \latexonly
!! \includegraphics[scale = 1]{season_flowchart_part1.png}
!! \includegraphics[scale = 1]{season_flowchart_part2.png}
!! \includegraphics[scale = 1]{season_flowchart_part3.png}
!! \endlatexonly
!! \n   
!_ =================================================================================================================================

  SUBROUTINE season_pre_disturbance (npts, dt, veget, veget_max, &
       vegstress_day, t2m_daily, tsoil_daily, lalo,  &
       precip_daily, npp_daily, circ_class_biomass, circ_class_n, &
       turnover_daily, gpp_daily, when_growthinit, &
       resp_maint_daily, resp_maint_week, &
       maxvegstress_lastyear, maxvegstress_thisyear, &
       minvegstress_lastyear, minvegstress_thisyear, &
       maxgppweek_lastyear, maxgppweek_thisyear, &
       gdd0_lastyear, gdd0_thisyear, &
       precip_lastyear, precip_thisyear, precip_longterm, &
       lm_lastyearmax, lm_thisyearmax, &
       maxfpc_lastyear, maxfpc_thisyear, &
       vegstress_month, vegstress_week, t2m_longterm, &
       tau_longterm, vpd_daily_mean, vpd_daily_max, vpd_mean_week, vpd_max_week, &
       t2m_month, t2m_week, &
       tsoil_month, precip_month,&
       npp_longterm, croot_longterm, turnover_longterm, gpp_week, &
       gpp_year, gpp_decade, plant_status, &
       gdd_m5_dormance, gdd_midwinter, ncd_dormance, ngd_minus5, &
       time_hum_min, hum_min_dormance, gdd_init_date , & 
       gdd_from_growthinit, herbivores, &
       Tseason, Tseason_length, Tseason_tmp, &
       Tmin_spring_time, t2m_min_daily, t2m_max_daily,&
       cn_leaf_min_season,nstress_season, vegstress_season, &
       rue_longterm, cn_leaf_init_2D, litter, leaf_age_crit, leaf_classes, &
       wood_leftover_legacy, season_drought_legacy, beetle_generation_index, &
       woody_litter_to_use, sumTeff, beetle_diapause, ni_acc)


    !! 0. Variable and parameter declaration

    !
    !! 0.1 Input variables
    !
    INTEGER(i_std), INTENT(in)                             :: npts              !! Domain size - number of grid cells (unitless)
    REAL(r_std), INTENT(in)                                :: dt                !! time step in days (dt_days)
    REAL(r_std), DIMENSION(:,:), INTENT(in)                :: veget             !! coverage fraction of a PFT. Here: fraction of 
                                                                                !! total ground. (0-1, unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(in)                :: veget_max         !! "maximal" coverage fraction of a PFT (for LAI -> 
                                                                                !! infinity) (0-1, unitless)Here: fraction of 
                                                                                !! total ground. 
    REAL(r_std), DIMENSION(:,:), INTENT(in)                :: vegstress_day  !! Daily moisture availability (0-1, unitless)
    REAL(r_std), DIMENSION(:), INTENT(in)                  :: t2m_daily         !! Daily 2 meter temperature (K)
    REAL(r_std), DIMENSION(:,:), INTENT(in)                :: tsoil_daily       !! Daily soil temperature (K)
    REAL(r_std), DIMENSION(:,:), INTENT(in)                :: lalo              !!  array of lat/lon
    REAL(r_std), DIMENSION(:), INTENT(in)                  :: precip_daily      !! Daily mean precipitation @tex ($mm day^{-1}$) 
                                                                                !! @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(in)                :: npp_daily         !! daily net primary productivity @tex ($gC m^{-2} 
                                                                                !! day^{-1}$) @endtex
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(in)          :: circ_class_biomass!! biomass @tex ($gC m^{-2} of ground$) @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)              :: circ_class_n      !! biomass @tex ($gC m^{-2} of ground$) @endtex

    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)            :: turnover_daily    !! Turnover rates @tex ($gC m^{-2} day^{-1}$) 
                                                                                !! @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(in)                :: gpp_daily         !! daily gross primary productivity  
                                                                                !! (Here: @tex $gC m^{-2} of total ground 
                                                                                !! day^{-1}$) @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(in)                :: resp_maint_daily  !! daily maintenance respiration
                                                                                !! (Here: @tex $gC m^{-2} of total ground 
                                                                                !! day^{-1}$) @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(in)                :: when_growthinit   !! how many days ago was the beginning of the 
                                                                                !! growing season (days) 
    REAL(r_std), DIMENSION(:), INTENT(in)                  :: t2m_min_daily     !! Daily minimum 2-meter temperature (K)
    REAL(r_std), DIMENSION(:), INTENT(in)                  :: t2m_max_daily     !! Daily maximum 2-meter temperature (K)  
    REAL(r_std), DIMENSION(:), INTENT(in)                  :: vpd_daily_mean    !! Daily mean VPD
    REAL(r_std), DIMENSION(:), INTENT(in)                  :: vpd_daily_max     !! Daily max VPD
    REAL(r_std), DIMENSION(:,:), INTENT(inout)             :: plant_status      !! Growth and phenological status of the plant       
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)           :: wood_leftover_legacy
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)           :: season_drought_legacy
    REAL(r_std), DIMENSION(:,:), INTENT(inout)             :: sumTeff !! Yearly thermal sum that is effective for bark beetles development
    REAL(r_std), DIMENSION(:,:), INTENT(inout)             :: woody_litter_to_use     !! woody litter pool to use for wood_left_over in the pest module
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)           :: beetle_generation_index !! Number of generation that BB can achieved in one year 
    INTEGER(i_std), DIMENSION(:,:), INTENT(inout)          :: beetle_diapause         !! 
    REAL(r_std), DIMENSION(:,:), INTENT(in)                :: cn_leaf_init_2D   !! initial leaf C/N ratio 
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(in)          :: litter

    !
    !! 0.2 Output variables 
    ! (diagnostic)
    !
    REAL(r_std), DIMENSION(:,:), INTENT(out)               :: herbivores        !! time constant of probability of a leaf to be 
                                                                                !! eaten by a herbivore (days) 

    !
    !! 0.3 Modified variables
    !
    REAL(r_std), DIMENSION(:,:), INTENT(inout)             :: maxvegstress_lastyear  !! last year's maximum moisture 
                                                                                        !! availability (0-1, unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)             :: maxvegstress_thisyear  !! this year's maximum moisture 
                                                                                        !! availability (0-1, unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)             :: minvegstress_lastyear  !! last year's minimum moisture 
                                                                                        !! availability (0-1, unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)             :: minvegstress_thisyear  !! this year's minimum moisture 
                                                                                        !! availability (0-1, unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)             :: maxgppweek_lastyear       !! last year's maximum weekly GPP 
                                                                                        !! @tex ($gC m^{-2} week^{-1}$) @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)             :: maxgppweek_thisyear       !! this year's maximum weekly GPP 
                                                                                        !! @tex ($gC m^{-2} week^{-1}$) @endtex 
    REAL(r_std), DIMENSION(:), INTENT(inout)               :: gdd0_lastyear             !! last year's annual GDD0 (C)
    REAL(r_std), DIMENSION(:), INTENT(inout)               :: gdd0_thisyear             !! this year's annual GDD0 (C)
    REAL(r_std), DIMENSION(:), INTENT(inout)               :: precip_lastyear           !! last year's annual precipitation 
                                                                                        !! @tex ($mm year^{-1}$) @endtex 
    REAL(r_std), DIMENSION(:), INTENT(inout)               :: precip_thisyear           !! this year's annual precipitation 
                                                                                        !! @tex ($mm year^{-1}$) @endtex
    REAL(r_std), DIMENSION(:), INTENT(inout)               :: precip_longterm           !! longterm annual precipitation sum 
                                                                                        !! @tex ($mm year^{-1}$) @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)             :: lm_lastyearmax            !! last year's maximum leaf mass, for each 
                                                                                        !! PFT @tex ($gC m^{-2}$) @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)             :: lm_thisyearmax            !! this year's maximum leaf mass, for each 
                                                                                        !! PFT @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)             :: maxfpc_lastyear           !! last year's maximum fpc for each PFT, on 
                                                                                        !! ground (0-1, unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)             :: maxfpc_thisyear           !! this year's maximum fpc for each PFT, on 
                                                                                        !! ground (0-1, unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)             :: vegstress_month        !! "monthly" moisture availability 
                                                                                        !! (0-1, unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)             :: vegstress_week         !! "weekly" moisture availability 
                                                                                        !! (0-1, unitless)
    REAL(r_std), DIMENSION(:), INTENT(inout)               :: t2m_longterm              !! "long term" 2-meter temperatures (K)
    REAL(r_std), INTENT(inout)                             :: tau_longterm     
    REAL(r_std), DIMENSION(:), INTENT(inout)               :: t2m_month                 !! "monthly" 2-meter temperatures (K)
    REAL(r_std), DIMENSION(:), INTENT(inout)               :: t2m_week                  !! "weekly" 2-meter temperatures (K)
    REAL(r_std), DIMENSION(:), INTENT(inout)               :: Tseason                   !! "seasonal" 2-meter temperatures (K), 
                                                                                        !! used to constrain boreal treeline
    REAL(r_std), DIMENSION(:), INTENT(inout)               :: Tseason_tmp               !! temporary variable to calculate Tseason 
    REAL(r_std), DIMENSION(:), INTENT(inout)               :: Tseason_length            !! temporary variable to calculate Tseason: 
                                                                                        !! number of days when t2m_week is higher than 0 degree 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)             :: Tmin_spring_time          !! Number of days after begin_leaves (leaf onset)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)             :: tsoil_month               !! "monthly" soil temperatures (K)
    REAL(r_std), DIMENSION(:), INTENT(inout)               :: precip_month              !! "monthly" preciptation ($mm month^{-1}$)
    REAL(r_std), DIMENSION(:), INTENT(inout)               :: vpd_mean_week             !! Weekly mean VPD
    REAL(r_std), DIMENSION(:), INTENT(inout)               :: vpd_max_week              !! Weekly mean maximum VPD
    REAL(r_std), DIMENSION(:,:), INTENT(inout)             :: npp_longterm              !! "long term" net primary productivity 
                                                                                        !! @tex ($gC m^{-2} year^{-1}$) @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)             :: croot_longterm            !! "long term" root carbon mass  
 	                                                                                !! @tex ($gC m^{-2}) @endtex
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)         :: turnover_longterm         !! "long term" turnover rate 
                                                                                        !! @tex ($gC m^{-2} year^{-1}$) @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)             :: gpp_week                  !! "weekly" GPP @tex ($gC m^{-2} day^{-1}$)@endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)             :: gpp_year                  !! "annual" GPP @tex ($gC m^{-2} day^{-1}$)@endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)             :: gpp_decade                !! "decadal" GPP @tex ($gC m^{-2} day^{-1}$)@endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)             :: resp_maint_week           !! "weekly" maintenance respiration @tex ($gC m^{-2} day^{-1}$) 
                                                                                        !! @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)             :: gdd_m5_dormance           !! growing degree days above threshold -5 
                                                                                        !! deg. C (C) 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)             :: gdd_midwinter             !! growing degree days since midwinter (C)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)             :: ncd_dormance              !! number of chilling days since leaves 
                                                                                        !! were lost (days) 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)             :: ngd_minus5                !! number of growing days (days)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)             :: time_hum_min              !! time elapsed since strongest moisture 
                                                                                        !! availability (days) 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)             :: hum_min_dormance          !! minimum moisture during dormance 
                                                                                        !! (0-1, unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)             :: gdd_init_date             !! inital date for gdd count
    REAL(r_std), DIMENSION(:,:), INTENT(inout)             :: gdd_from_growthinit       !! growing degree days, threshold 0 deg. C 
                                                                                        !! since beginning of season
    REAL(r_std), DIMENSION(:,:), INTENT(inout)             :: cn_leaf_min_season        !! seasonal minimal leaf nitrogen concentratation CN ratio  
	                                                                                !! ((gC) (gN)-1) 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)             :: nstress_season            !! Seasonal nitrogen stress 
 	   
    REAL(r_std), DIMENSION(:,:), INTENT(inout)             :: vegstress_season       !! soil moisture during growing season
                                                                                        !! (used for allometry&allocation)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)             :: rue_longterm              !! Longterm radiation use efficiency 
                                                                                        !! (units?)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)             :: leaf_age_crit             !! critical leaf age (days)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)             :: leaf_classes              !! width of each leaf age class (days) 
    REAL(r_std), DIMENSION(:), INTENT(inout)               :: ni_acc                    !! Nesterov index (square of degree Celcius)
    !
    !! 0.4 Local variables
    !
    INTEGER(i_std)                                          :: j, ipts, iyear                 !! indices (unitless)
    REAL(r_std)                                             :: ncd_max                  !! maximum ncd (to avoid floating point 
                                                                                        !! underflows) (days) 
    REAL(r_std), DIMENSION(npts,nvm)                        :: is_growing_season         !! temporal mask for XIOS
    REAL(r_std), DIMENSION(npts,nvm)                        :: is_pot_growing_season     !! temporal mask for XIOS (potential)
    REAL(r_std), DIMENSION(npts)                            :: sumfpc_nat               !! sum of natural fpcs (0-1, unitless)
                                                                                        !! [DISPENSABLE] 
    REAL(r_std), DIMENSION(npts,nvm)                        :: weighttot                !! weight of biomass @tex ($gC m^{-2}$) 
                                                                                        !! @endtex 
    REAL(r_std), DIMENSION(npts,nvm)                        :: nlflong_nat              !! natural long-term leaf NPP 
                                                                                        !! @tex ($gC m^{-2} year^{-1}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm)                        :: green_age                !! residence time of green tissue (years)
    REAL(r_std), DIMENSION(npts)                            :: consumption              !! herbivore consumption 
                                                                                        !! @tex ($gC m^{-2} day^{-1}$) @endtex
    REAL(r_std), DIMENSION(npts)                            :: fracnat                  !! fraction of each gridcell occupied by 
                                                                                        !! natural vegetation (0-1, unitless)
    REAL(r_std), DIMENSION(npts)                            :: solad                    !! maximal radation during current day 
                                                                                        !! (clear sky condition)
    REAL(r_std), DIMENSION(npts)                            :: solai                    !! maximal radation during current day
                                                                                        !! (clear sky condition)
    REAL(r_std), DIMENSION(npts)                            :: cloud                    !! cloud fraction
    REAL(r_std), DIMENSION(npts)                            :: nn                       !! Temporary variable for storing biomass(initrogen) for the sum of labile+root+leaf pools (gN m-2) 
    REAL(r_std), DIMENSION(npts)                            :: cn                       !! CN ration for  labile+root+leaf pools (gC (gN)-1) 
    REAL(r_std), DIMENSION(npts,nvm,legacy_years)           :: tmp_legacy
    REAL(r_std), DIMENSION(npts,nvm,legacy_years_wood)      :: tmp_legacy2
    REAL(r_std), DIMENSION(npts)                            :: vartmp                   !! Temporary variable
    REAL(r_std), DIMENSION(npts,nvm)                        :: gtemp 
!_ =================================================================================================================================

    IF (firstcall_season_pre_disturbance) THEN
      !! Initialize local printlev
      ! Guillaume M. -- get_printlev calls getin_p, a COLLECTIVE bcast. Called here,
      ! inside the compute loop under a local firstcall, it hangs the run as soon as the
      ! ranks stop reaching it together. Use printlev, read once at startup by every rank.
      printlev_loc=printlev
    END IF

    IF (printlev_loc>=2) WRITE(numout,*) 'Entering season'

    !! 1. Initializations
    !! 1.1 Calculate ::ncd_max - the maximum possible NCD (number of chilling days) as:
    !!     \latexonly
    !!     \input{season_ncdmax_eqn2.tex}
    !!     \endlatexonly
    !!     \n
    !!     where one_year is 1 year in seconds (defined in ::constantes).
    !
    ncd_max = ncd_max_year * one_year

    IF ( firstcall_season_pre_disturbance ) THEN

       !
       !! 1.2 first call - output message giving the setting of the ::gppfrac_dormance,
       !!     ::hvc1, ::hvc2 and ::leaf_frac (as a percentage) parameters which are 
       !!     defined at the beginning of this subroutine. Also outputs the value of 
       !!     ::ncd_max.
       !

       IF ( printlev_loc>=3) THEN

          WRITE(numout,*) 'season: '

          WRITE(numout,*) '   > maximum GPP/GGP_max ratio for dormance (::gppfrac_dormance) :',gppfrac_dormance !*

          WRITE(numout,*) '   > maximum possible ncd(days) (::ncd_max) : ',ncd_max !*

          WRITE(numout,*) '   > herbivore consumption C (gC/m2/d) as a function of NPP (gC/m2/d): (::hvc1) (::hvc2)' !*
          WRITE(numout,*) '     C=',hvc1,' * NPP^',hvc2
          WRITE(numout,*) '   > for herbivores, suppose that (::leaf_frac_hvc)',leaf_frac_hvc*100., &
               '% of NPP is allocated to leaves'

       ENDIF

       !
       !! 1.3 Initialize monthly and weekly variables with daily values if they were not in the restart file
       !      This initalization is done here and not in readreastart(stomate_io.f90) because the daily variables 
       !      are calculated during the first day in stomate, after readrestart.

       IF ( ALL(vegstress_month(:,:) == val_exp) )  THEN
          ! The variable is not found in restart file
          IF(printlev_loc>=2) WRITE(numout,*) 'Warning! We have to initialize the ''monthly'' moisture availabilities. '
          vegstress_month(:,:) = vegstress_day(:,:)
       ENDIF

       IF ( ALL(vegstress_week(:,:) == val_exp) ) THEN
          ! The variable is not found in restart file
          IF(printlev_loc>=2) WRITE(numout,*) 'Warning! We have to initialize the ''weekly'' moisture availabilities. '
          vegstress_week(:,:) = vegstress_day(:,:)
       ENDIF

       !! 1.4 reset firstcall_season flag
       firstcall_season_pre_disturbance = .FALSE.

    ENDIF

    ! determine min yearly clear sky radiation (solstice) as beginning for gdd count
    cloud(:) = zero
    CALL downward_solar_flux (npts, lalo(:,1),julian_diff,12.,cloud,1,solad,solai)

    ! Shortest day of the year
    WHERE (solad(:) .LT. gdd_init_date(:,2))
       gdd_init_date(:,1)= julian_diff
       gdd_init_date(:,2)= solad(:)
    ENDWHERE

    !
    !! NOTE: Sections 2. to 5. compute slowly-varying, "long-term" (i.e. weekly/monthly)
    !! input variables using a linear relaxation method, following the equation:
    !! \latexonly
    !! \input{season_lin_relax_eqn1.tex}
    !! \endlatexonly
    !! \n 
    !! as described in the introduction to this subroutine (see Krinner et al., 2005).  
    !! The time constant in the above equation is given for each of the variables
    !! described.
    !

    !
    !! 2. Moisture availability (relative soil moisture, including the root profile) :
    !!    The time constants (as in the above equation) for these calculations are 
    !!    given by the parameters ::tau_hum_month (for the monthly
    !!    calculation) or ::tau_hum_week (for the weekly),
    !!    which are set in ::stomate_data to be 20 and 7 days respectively.
    !!    If the moisture availability is less than zero, it is set to zero. 
    !!    These variables are mostly used in the phenological leaf onset and senescence
    !!    models in the modules ::stomate_phenology and ::stomate_turnover, 
    !!    respectively. They are used in models which require a moisture limitation
    !!    condition. These are the 'hum', 'moi', 'humgdd' and 'moigdd' onset models,
    !!    and the 'mixed' and 'dry' (weekly only) senescence models. In addition,
    !!    the weekly moisture availability is used to calculate the limitation
    !!    on the fraction of NPP allocated to different compartments in the module
    !!    ::stomate_alloc. 
    !

    !
    !! 2.1 "monthly" (::vegstress_month)
    !
    vegstress_month = ( vegstress_month * ( tau_hum_month - dt ) + &
         vegstress_day * dt ) / tau_hum_month

    DO j = 2,nvm  ! Loop over # PFTs
       WHERE ( ABS(vegstress_month(:,j)) .LT. EPSILON(zero) )
          vegstress_month(:,j) = zero
       ENDWHERE
    ENDDO

    !
    !! 2.2 "weekly" (::vegstress_week)
    !

    vegstress_week = ( vegstress_week * ( tau_hum_week - dt ) + &
         vegstress_day * dt ) / tau_hum_week

    DO j = 2,nvm ! Loop over # PFTs
       WHERE ( ABS(vegstress_week(:,j)) .LT. EPSILON(zero) ) 
          vegstress_week(:,j) = zero
       ENDWHERE
    ENDDO

    !
    !! 3. 2-meter temperatures
    !!    The time constants for the "long-term", "monthly" and "weekly" 2-meter
    !!    temperatures are given by the parameters ::tau_longterm_max, 
    !!    ::tau_t2m_month, and ::tau_t2m_week,
    !!    which are set in ::stomate_data to be 3 * one year (in seconds, as 
    !!    described above) and 20 and 7 days respectively.
    !!    If the temperature is less than zero, it is set to zero.
    !!    In addition the "long term" temperature is written to the history file. \n
    !!    These variables are used in many different modules of STOMATE. 
    !!    The longterm t2m is limied to the intervall [tlong_ref_min, tlong_ref_max]
    !!    using the equation:
    !!      \latexonly
    !!      \input{season_t2m_long_ref_eqn3.tex}
    !!      \endlatexonly
    !!      \n
    !!    The "monthly" and "weekly" temperature is used in the modules 
    !!    ::stomate_phenology and ::stomate_turnover for the onset and senescence 
    !!    phenological models which require a temperature condition. In addition
    !!    the "monthly" temperature is used in ::lpj_constraints to determine 
    !!    the presence and regeneration of the vegetation and in the module
    !!    ::stomate_assimtemp to calculate photosynthesis temperatures.
        
    ! Update tau_longterm
    tau_longterm = MIN(tau_longterm+dt,tau_longterm_max)
    
    ! Recalculate the reference temperature using the old reference temperature and the current temperature
    t2m_longterm(:) = ( t2m_longterm(:) * ( tau_longterm - dt ) + &
         t2m_daily(:) * dt ) / tau_longterm

    ! The longterm reference is not allowed to go outside the interval [tlong_ref_min, tlong_ref_max]
    t2m_longterm(:) = MAX( tlong_ref_min, MIN( tlong_ref_max, t2m_longterm(:) ) )

    CALL xios_orchidee_send_field("t2m_longterm",t2m_longterm)

    CALL histwrite_p (hist_id_stomate, 'T2M_LONGTERM', itime, &
         t2m_longterm, npts, hori_index)

    !
    !! 3.3 "monthly" (::t2m_month)
    !

    t2m_month = ( t2m_month * ( tau_t2m_month - dt ) + &
         t2m_daily * dt ) / tau_t2m_month

    WHERE ( ABS(t2m_month(:)) .LT. EPSILON(zero) )
       t2m_month(:) = zero
    ENDWHERE
    
    ! Calculate "seasonal" temperature
    WHERE ( t2m_week (:) .GT. ZeroCelsius )
       Tseason_tmp(:) = Tseason_tmp(:) + t2m_week(:)
       Tseason_length(:)=Tseason_length(:) + dt
    ENDWHERE

    ! calculate Tmin in spring (after onset)
    DO j=1, nvm
       IF (leaf_tab(j)==1 .AND. pheno_type(j)==2) THEN
          ! leaf_tab=broadleaf and pheno_typ=summergreen
          ! Only treat broadleag and summergreen pfts
          
          WHERE ( Tmin_spring_time(:,j)>0 .AND. (Tmin_spring_time(:,j)<spring_days_max) )
             Tmin_spring_time(:,j)=Tmin_spring_time(:,j)+1
          ELSEWHERE
             Tmin_spring_time(:,j)=0
          ENDWHERE
          
          WHERE ( plant_status(:,j) .EQ. ibudbreak )
             Tmin_spring_time(:,j)=1
          ENDWHERE
       END IF
    END DO

    !
    !! 3.4 "weekly" (::t2m_week)
    !
    t2m_week = ( t2m_week * ( tau_t2m_week - dt ) + &
         t2m_daily * dt ) / tau_t2m_week

    WHERE ( ABS(t2m_week(:)) .LT. EPSILON(zero) )
       t2m_week(:) = zero
    ENDWHERE
    
    !
    !! 3.5 "weekly" (::vpd_mean_week) (::vpd_max_week)
    !
    vpd_mean_week = (vpd_mean_week * (tau_week - dt) + &
             vpd_daily_mean * dt ) / tau_week
    WHERE ( ABS(vpd_mean_week(:)) .LT. EPSILON(zero) )
       vpd_mean_week(:)=zero
    ENDWHERE
    
    vpd_max_week = (vpd_max_week * (tau_week - dt) + &
             vpd_daily_max * dt ) / tau_week
    WHERE ( ABS(vpd_max_week(:)) .LT. EPSILON(zero) )
       vpd_max_week(:)=zero
    ENDWHERE    
    
    !
    !! 4. "monthly" soil temperatures (::tsoil_month)
    !!    The time constant is given by the parameter ::tau_tsoil_month, 
    !!    which is set in ::stomate_data to be 20 days.
    !!    If the monthly soil temperature is less than zero, it is set to zero. \n
    !!    This variable is used in the ::stomate_allocation module.
    !

    tsoil_month = ( tsoil_month * ( tau_tsoil_month - dt ) + &
         tsoil_daily(:,:) * dt ) / tau_tsoil_month

    WHERE ( ABS(tsoil_month(:,:)) .LT. EPSILON(zero) )
       tsoil_month(:,:) = zero
    ENDWHERE


    !! "monthly" precipitation
    precip_month = ( precip_month * ( tau_precip_month - dt ) + &
         precip_daily(:) * dt ) / tau_precip_month
    WHERE ( ABS(precip_month(:)) .LT. EPSILON(zero) )
        precip_month(:) = zero
    ENDWHERE

 
    !! 24. Calculate the critical leaf age
    ! longevity_leaf is a prescribed parameter for the longevity of a 
    ! typical leaf/needle at the average temperature for that PFT. For
    ! PFTs where a large range in longevity has been observed, this is
    ! accounted for in the calculation of leaf_age_crit. Leaf_age_crit is 
    ! thus the location-specific longevity. It the age of a leaf/needle
    ! exceeds crit_leaf_age, leaf turnover will start (see stomate_turnover).
    leaf_age_crit(:,1) = zero
    leaf_classes(:,1) = zero

    DO j = 2,nvm 

       !! 4.1 Calculate critical leaf age       
       !  Critical leaf age depends on long-term temperature
       !  generally, lower turnover in cooler climates. For several
       !  PFTs the parameters are set such that leaf_age_crit = longevity_leaf 
       leaf_age_crit(:,j) = &
            MIN( longevity_leaf(j) * leaf_age_crit_coeff1(j) , &
            MAX( longevity_leaf(j) * leaf_age_crit_coeff2(j) , &
            longevity_leaf(j) - leaf_age_crit_coeff3(j) * &
            ( t2m_longterm(:) - ZeroCelsius - leaf_age_crit_tref(j) ) ) )

       ! Calculation how many days each leaf age class contains
       leaf_classes(:,j) = leaf_age_crit(:,j) / REAL(nleafages,r_std)

    END DO

    !
    !! 6. Calculate the dormance time-length (::time_lowgpp).
    !!    The dormancy time-length is increased by the stomate time step when   
    !!    the weekly GPP is below one of two thresholds, 
    !!    OR if the growing season started more than 2 years ago AND the amount of biomass
    !!    in the carbohydrate reserve is 4 times as much as the leaf biomass.
    !!    i.e. the plant is declared dormant if it is accumulating carbohydrates but not 
    !!    using them, so that the beginning of the growing season can be detected .
    !     This last condition was added by Nicolas Viovy.
    !!    - NV: special case (3rd condition)
    !!    Otherwise, set ::time_lowgpp to zero. \n
    !!    The weekly GPP must be below either the parameter ::min_gpp_allowed, which is set at 
    !!    the beginning of this subroutine to be 0.3gC/m^2/year, 
    !!    OR it must be less than last year's maximum weekly GPP multiplied by the
    !!    maximal ratio of GPP to maximum GPP (::gppfrac_dormance), which is set to 0.2 at the 
    !!    beginning of this subroutine. \n
    !!    This variable is used in the ::stomate_phenology module to allow the detection of the
    !!    beginning of the vegetation growing season for different onset models.
    !!    Each PFT is looped over.
    !
    !NVMODIF
!!$    DO j = 2,nvm ! Loop over # PFTs
!!$       WHERE ( ( gpp_week(:,j) .LT. min_gpp_allowed ) .OR. & 
!!$            ( gpp_week(:,j) .LT. gppfrac_dormance * maxgppweek_lastyear(:,j) ) .OR. &
!!$            ( ( when_growthinit(:,j) .GT. 2.*one_year ) .AND. &
!!$            ( biomass(:,j,icarbres,icarbon) .GT. biomass(:,j,ileaf,icarbon)*4. ) ) )
!!$          !       WHERE ( ( gpp_week(:,j) .EQ. zero ) .OR. & 
!!$          !            ( gpp_week(:,j) .LT. gppfrac_dormance * maxgppweek_lastyear(:,j) ) .OR. &
!!$          !            ( ( when_growthinit(:,j) .GT. 2.*one_year ) .AND. &
!!$          !            ( biomass(:,j,icarbres,icarbon) .GT. biomass(:,j,ileaf,icarbon)*4. ) ) )
!!$       
!!$          time_lowgpp(:,j) = time_lowgpp(:,j) + dt
!!$          
!!$       ELSEWHERE
!!$          
!!$          time_lowgpp(:,j) = zero
!!$
!!$       ENDWHERE
!!$    ENDDO

    !
    !! 7. Calculate the growing degree days (GDD - ::gdd_m5_dormance).
    !!    This variable is the GDD sum of the temperatures higher than -5 degrees C.
    !!    It is inititalised to 0 at the beginning of the dormancy period (i.e.
    !!    when ::time_lowgpp>0), and is set to "undef" when there is GPP 
    !!    (::time_lowgpp=0), as it is not used during the growing season. \n 
    !!    ::gdd_m5_dormance is further scaled by (::tau_gdd -dt / ::tau_gdd), 
    !!    - ::tau_gdd is set to 40. in the module ::stomate_constants. \n
    !     Nicolas Viovy - not sure why this is...
    !!    This variable is used in the ::stomate_phenology module for the 
    !!    'humgdd' and 'moigdd' leaf onset models. \n
    !!    Each PFT is looped over but ::gdd_m5_dormance is only calculated for 
    !!    those PFTs for which a critical GDD is defined, i.e. those which are 
    !!    assigned to the 'humgdd' or 'moigdd' models. \n
    !!    Finally if GDD sum is less than zero, then it is set to zero.
    DO j = 2,nvm ! Loop over # PFTs

       IF (.NOT. natural(j)) THEN
          ! reset counter: start of the growing season
          WHERE ((when_growthinit(:,j) .EQ. zero))
	     gdd_from_growthinit(:,j) = zero
	  ENDWHERE
          ! increase gdd counter
          WHERE ( t2m_daily(:) .GT. (ZeroCelsius) .AND. &
               veget_max(:,j) .GT. min_stomate)
             gdd_from_growthinit(:,j) = gdd_from_growthinit(:,j) + &
                                 dt * ( t2m_daily(:) - (ZeroCelsius) )
          ENDWHERE
       ELSE
          gdd_from_growthinit(:,j) = zero
       ENDIF

       ! only for PFTs for which critical gdd is defined
       ! gdd_m5_dormance is set to 0 at the end of the growing season. It is set to undef
       ! at the beginning of the growing season.
       IF ( ALL(pheno_gdd_crit(j,:) .NE. undef) ) THEN

          !
          !! 7.1 set to zero if undefined and there is no GPP
          !
          ! When the shorest day is used to start counting global
          ! maps of the start (and in turn end) of the growing season
          ! show a sudden shift in growing season by 180 days. This
          ! could be seen by a straight line along the equator.
          WHERE (gdd_init_date(:,1) .EQ. julian_diff)

             gdd_m5_dormance(:,j) = zero

          ENDWHERE

          !
          !! 7.2 set to undef if there is GPP
          !

          WHERE ( when_growthinit(:,j) .EQ. zero )

             gdd_m5_dormance(:,j) = undef

          ENDWHERE

          !
          !! 7.3 normal update as described above where ::gdd_m5_dormance 
          !  is defined (not set to "undef").
          !
          WHERE ( ( t2m_daily(:) .GT. ( ZeroCelsius + gdd_threshold(j) ) ) .AND. &
               ( gdd_m5_dormance(:,j) .NE. undef )           )
             gdd_m5_dormance(:,j) = gdd_m5_dormance(:,j) + &
                  dt * ( t2m_daily(:) - ( ZeroCelsius + gdd_threshold(j) ) )
          ENDWHERE

       ENDIF

    ENDDO

    !
    !! 7.4 Set to zero if GDD is less than zero.

    DO j = 2,nvm ! Loop over # PFTs
       WHERE ( ABS(gdd_m5_dormance(:,j)) .LT. EPSILON(zero) )
          gdd_m5_dormance(:,j) = zero
       ENDWHERE
    ENDDO

    !
    !! 8. Calculate the growing degree days (GDD) since midwinter (::gdd_midwinter)
    !!    This variable represents the GDD sum of temperatures higher than a PFT-dependent
    !!    threshold (::ncdgdd_temp), since midwinter.
    !!    Midwinter is detected if the monthly temperature (::t2m_month) is lower than the weekly 
    !!    temperature (::t2m_week) AND if the monthly temperature is lower than the long-term 
    !!    temperature (t2m_longterm). These variables were calculated earlier in this subroutine. \n
    !!    ::gdd_midwinter is initialised to 0.0 when midwinter is detected, and from then on 
    !!    increased with each temperature greater than ::ncdgdd_temp, which is defined
    !!    in the module ::stomate_constants. \n
    !!    ::gdd_midwinter is set to "undef" when midsummer is detected, following the opposite
    !!    conditions to those used to define midwinter. \n
    !!    The variable is used in the ::stomate_phenology module for the leaf onset model 'ncdgdd'.
    !!    Each PFT is looped over but the ::gdd_midwinter is only calculated for those
    !!    PFTs for which a critical 'ncdgdd' temperature is defined, i.e. those which are 
    !!    assigned to the 'ncdgdd' model. \n
    !

    DO j = 2,nvm ! Loop over # PFTs

       ! only for PFTs for which ncdgdd_crittemp is defined

       IF ( ncdgdd_temp(j) .NE. undef ) THEN

          !
          !! 8.1 set to 0 if undef and if we detect "midwinter"
          !

!!$          WHERE ( ( gdd_midwinter(:,j) .EQ. undef ) .AND. &
!!$               ( t2m_month(:) .LT. t2m_week(:) ) .AND. &
!!$               ( t2m_month(:) .LT. t2m_longterm(:) )    )
          WHERE (gdd_init_date(:,1) .EQ. julian_diff)

             gdd_midwinter(:,j) = zero
          ENDWHERE

          !
          !! 8.2 set to undef if we detect "midsummer"
          !

!!$          WHERE ( ( t2m_month(:) .GT. t2m_week(:) ) .AND. &
!!$               ( t2m_month(:) .GT. t2m_longterm(:) )    )
!!$
!!$             gdd_midwinter(:,j) = undef
!!$
!!$          ENDWHERE
          ! When_growthinit can be zero because the PFT
          ! does not exist and it can be zero because the
          ! PFT is dormant. The WHERE statements should
          ! distinguish between those cases.
          WHERE ( when_growthinit(:,j) .EQ. zero .AND. &
               veget_max(:,j) .GT. min_stomate)

             gdd_midwinter(:,j) = undef

          ENDWHERE

          !
          !! 8.3 normal update as described above
          !

          WHERE ( gdd_midwinter(:,j) .NE. undef .AND. &
               t2m_daily(:) .GT. ncdgdd_temp(j)+ZeroCelsius .AND. &
               veget_max(:,j) .GT. min_stomate )

             gdd_midwinter(:,j) = &
                  gdd_midwinter(:,j) + &
                  dt * ( t2m_daily(:) - ( ncdgdd_temp(j)+ZeroCelsius ) )

          ENDWHERE

       ENDIF

    ENDDO

    !
    !! 9. Calculate the number of chilling days (NCD) since leaves were lost (::ncd_dormance).
    !!    This variable is initialised to 0 at the beginning of the dormancy period (::time_lowgpp>0)
    !!    and increased by the stomate time step when the daily temperature is lower than the 
    !!    PFT-dependent threshold ::ncdgdd_temp, which is defined for each PFT 
    !!    in a table (::ncdgdd_temp_tab) in the module ::stomate_data. \n
    !!    It is set to "undef" when there is GPP (::time_lowgpp=0) as it is not needed during
    !!    the growing season. \n
    !!    The variable is used in the ::stomate_phenology module for the leaf onset model 'ncdgdd'.
    !!    Each PFT is looped over but the ::ncd_dormance is only calculated for those
    !!    PFTs for which a critical 'ncdgdd' temperature is defined, i.e. those which are 
    !!    assigned to the 'ncdgdd' model.
    !

    DO j = 2,nvm ! Loop over # PFTs

       IF ( ncdgdd_temp(j) .NE. undef ) THEN

          !
          !! 9.1 set to zero if undefined and there is no GPP
          !

          WHERE (gdd_init_date(:,1) .EQ. julian_diff .AND. &
               veget_max(:,j) .GT. min_stomate)

             ncd_dormance(:,j) = zero

          ENDWHERE

          !
          !! 9.2 set to undef if there is GPP
          !

          WHERE ( when_growthinit(:,j) .LE. un .AND. &
               veget_max(:,j) .GT. min_stomate)

             ncd_dormance(:,j) = undef

          ENDWHERE

          !
          !! 9.3 normal update, as described above, where ::ncd_dormance is defined
          !

          WHERE ( ncd_dormance(:,j) .NE. undef .AND. &
               t2m_daily(:) .LE. ncdgdd_temp(j)+ZeroCelsius .AND. &
               veget_max(:,j) .GT. min_stomate )

             ncd_dormance(:,j) = MIN( ncd_dormance(:,j) + dt, ncd_max )

          ENDWHERE

       ENDIF

    ENDDO

    !
    !! 10. Calculate the number of growing days (NGD) since leaves were lost (::ngd_minus5).
    !!     This variable is initialised to 0 at the beginning of the dormancy period (::time_lowgpp>0)
    !!     and increased by the stomate time step when the daily temperature is higher than the threshold
    !!     -5 degrees C. \n    
    !!     ::ngd_minus5 is further scaled by (::tau_ngd -dt / ::tau_ngd), 
    !!     - ::tau_ngd is set to 50. in the module ::stomate_constants. \n
    !      Nicolas Viovy - not sure why this is...
    !!     The variable is used in the ::stomate_phenology module for the leaf onset model 'ngd'.
    !!     Each PFT is looped over.
    !!     If the NGD is less than zero, it is set to zero.
    !

    DO j = 2,nvm ! Loop over # PFTs

       !
       !! 10.1 Where there is GPP (i.e. ::time_lowgpp=0), set NGD to 0. 
       !!      This means that we only take into account NGDs when the leaves are off
       !

       WHERE (gdd_init_date(:,1) .EQ. julian_diff)
          ngd_minus5(:,j) = zero
       ENDWHERE

       !
       !! 10.2 normal update, as described above.
       !

       WHERE ( t2m_daily(:) .GT. (ZeroCelsius + ngd_threshold) .AND. &
            veget_max(:,j) .GT. min_stomate)
          ngd_minus5(:,j) = ngd_minus5(:,j) + dt
       ENDWHERE
       
       !+++CHECK+++
       ! ngd_minus5 is further scaled by (::tau_ngd -dt / ::tau_ngd), 
       ! tau_ngd is set to 50. in the module ::stomate_constants. \n
       ! Nicolas Viovy - not sure why this is...
       ngd_minus5(:,j) = ngd_minus5(:,j) * ( tau_ngd - dt ) / tau_ngd
       !+++++++++++

    ENDDO

    WHERE ( ngd_minus5(:,:) .LT. zero )
       ngd_minus5(:,:) = zero
    ENDWHERE
    

    !
    !! 11. Calculate the minimum humidity/relative soil moisture since dormance began (::hum_min_dormance) 
    !!     and the time elapsed since this minimum (::time_hum_min). \n
    !!     The minimum moisture availability occurs when the monthly moisture availability, which is updated 
    !!     daily earlier in this subroutine as previously described, is at a minimum during the dormancy period.
    !!     Therefore the ::hum_min_dormance is initialised to the monthly moisture availability 
    !!     at the beginning of the dormancy period (i.e. if it was previously set to undefined and the
    !!     ::time_lowgpp>0) AND then whenever the monthly moisture availability is less than it has previously
    !!     been. \n
    !!     Consequently, the time counter (::time_hum_min) is initialised to 0 at the beginning of the dormancy  
    !!     period (i.e. if it was previously set to undefined and the ::time_lowgpp>0) AND when the minimum 
    !!     moisture availability is reached, and is increased by the stomate time step every day throughout
    !!     the dormancy period. \n
    !!     ::time_hum_min is used in the ::stomate_phenology module for the leaf onset models 'moi' and 'moigdd'.
    !!     Each PFT is looped over but the two variables are only calculated for those
    !!     PFTs for which the critical parameter ::hum_min_time is defined, i.e.   
    !!     those which are assigned to the 'moi' or 'moigdd' models.
    !

    DO j = 2,nvm ! Loop over # PFTs

       IF ( hum_min_time(j) .NE. undef ) THEN

          !
          !! 11.1 initialize if undefined and there is no GPP
          !

          WHERE (when_growthinit(:,j) .EQ. zero)

             time_hum_min(:,j) = zero
             hum_min_dormance(:,j) = vegstress_month(:,j)

          ENDWHERE

          !
          !! 11.3 normal update, as described above, where ::time_hum_min and ::hum_min_dormance are defined
          !

          !! 11.3.1 increase time counter by the stomate time step

          WHERE ( hum_min_dormance(:,j) .NE. undef )

             time_hum_min(:,j) = time_hum_min(:,j) + dt

          ENDWHERE

          !! 11.3.2 set time counter to zero if minimum is reached

          WHERE (( hum_min_dormance(:,j) .NE. undef ) .AND. &
               ( vegstress_month(:,j) .LE. hum_min_dormance(:,j) ) )

             hum_min_dormance(:,j) = vegstress_month(:,j)
             time_hum_min(:,j) = zero

          ENDWHERE

       ENDIF

    ENDDO

    !! NOTE: Sections 12. to 14. compute slowly-varying, "long-term" (i.e. weekly/monthly)
    !! C fluxes (NPP, turnover, GPP) using the same linear relaxation method as in 
    !! Sections 2. and 5. and described in the introduction to this section, as given 
    !! by the following equation:
    !! \latexonly
    !! \input{season_lin_relax_eqn1.tex}
    !! \endlatexonly
    !! \n 
    !! The following variables are calculated using the above equation, and the time constant 
    !! is given for each. 
    !

    !
    !! 12. Update the "long term" NPP. (::npp_daily in gC/m^2/day, ::npp_longterm in gC/m^2/year.)
    !!     The time constant is given by the parameter ::tau_longterm, 
    !!     which is set in ::stomate_data to be 3 * one year (in seconds, as 
    !!     described above). If the ::npp_longterm is less than zero then set it to zero. \n
    !!     ::npp_longterm is used in ::stomate_lpj in the calculation of long-term
    !!     vegetation dynamics and in ::stomate_lcchange in the calculation of 
    !!     land cover change. It is also used to calculate diagnose the hebivory activity in 
    !!     Section 22 of this subroutine.
    !

    npp_longterm = ( npp_longterm * ( tau_longterm - dt ) + &
         (npp_daily*one_year) * dt                          ) / &
         tau_longterm

    WHERE ( npp_longterm(:,:) .LT. EPSILON(zero) )
       npp_longterm(:,:) = zero
    ENDWHERE
    

    WHERE(SUM(circ_class_biomass(:,:,:,iroot,icarbon),3) .GT. min_stomate) 
       croot_longterm(:,:) = ( croot_longterm(:,:) * ( tau_longterm - dt ) + & 
            SUM(circ_class_biomass(:,:,:,iroot,icarbon) * circ_class_n(:,:,:),3) * dt) / & 
            tau_longterm 
    ENDWHERE

    !
    !! 13. Update the "long term" turnover rates (in gC/m^2/year).
    !!     The time constant is given by the parameter ::tau_longterm, 
    !!     which is set in ::stomate_data to be 3 * one year (in seconds, as 
    !!     described above). If the ::turnover_longterm is less than zero then set it to zero.\n
    !!     ::turnover_longterm is used in ::stomate_lpj and :: lpg_gap in the calculation 
    !!     of long-term vegetation dynamics.
    !

    turnover_longterm(:,:,:,:) = ( turnover_longterm(:,:,:,:) * ( tau_longterm - dt ) + &
         (turnover_daily(:,:,:,:)*one_year) * dt                          ) / &
         tau_longterm

    WHERE ( turnover_longterm(:,:,:,:) .LT. zero )
       turnover_longterm(:,:,:,:) = zero
    ENDWHERE
    

    !
    !! 14. Update the "weekly" GPP (where there is vegetation), otherwise set to zero.
    !!     The time constant is given by the parameter ::tau_gpp_week, 
    !!     which is set in ::stomate_data to be 7 days. If the ::gpp_week is 
    !!     less than zero then set it to zero. \n
    !!     ::gpp_week is used to update the annual maximum weekly GPP (::maxgppweek_thisyear)
    !!     in Section 16 of this subroutine, which is then used to update the variable
    !!     ::maxgppweek_lastyear in Section 21 of this subroutine. Both ::gpp_week and 
    !!     ::maxgppweek_lastyear are used in Section 6 of this subroutine to calculate
    !!     the onset and time-length of the dormancy period. 
    !      Note: Used to be weekly GPP divided by veget_max, i.e. per ground covered, but not anymore.
    !

    WHERE ( veget_max .GT. zero )

       gpp_week = ( gpp_week * ( tau_gpp_week - dt ) + &
            gpp_daily * dt ) / tau_gpp_week

    ELSEWHERE

       gpp_week = zero

    ENDWHERE

    ! Use min_stomate as a threshold because we divide by gpp_week later in the code
    WHERE ( gpp_week(:,:) .LT. min_stomate )
       gpp_week(:,:) = zero
    ENDWHERE

    ! Always calculate long-term (annual and decadal)
    ! Use a nested approach to surpress the seasonal variation. Calculate the
    ! annual gpp and use this annual value once per year to calculate a decadal
    ! gpp.
    DO j = 2,nvm

       IF (is_tree(j)) THEN
           
          WHERE ( veget_max(:,j) .GT. zero )
          
             gpp_year(:,j) = (gpp_year(:,j) * ( tau_gpp_year - dt ) + &
                  gpp_daily(:,j) * dt ) / tau_gpp_year
       
          ELSEWHERE

             ! Set to a number to avoid undefined values for gpp_decade
             gpp_year(:,j) = zero

          ENDWHERE
 
          IF ( ts_annual_proc ) THEN

             ! This is a gapfilled variable don't use ELSEWHERE
             ! to set gpp_decade to zero where veget_max = 0 
             WHERE ( veget_max(:,j) .GT. zero )
          
                gpp_decade(:,j) = (gpp_decade(:,j) * ( tau_gpp_self_thin - 1 ) + &
                     gpp_year(:,j)) / tau_gpp_self_thin
       
             ENDWHERE
             
          END IF   
          
       END IF
       
    END DO
    
    IF (calculate_gpp_preind .AND. &
         ts_annual_proc) THEN

       ! Calculate longterm mean GPP during the spinup. This GPP
       ! functions as the reference GPP during the spinup, transient,
       ! and historical simulations.
       DO j = 2,nvm

          IF (is_tree(j)) THEN

             ! Do not set pre_indust_ref_gpp to zero where veget_max = 0.
             ! We need a reference gpp even if the PFT is not yet present
             ! in the pixel.
             WHERE ( veget_max(:,j) .GT. zero )

                WHERE ( gpp_year(:,j) .GE. zero )

                   ! Calculate the reference gpp if gpp >= zero. This
                   ! should be the dominant case. If gpp is frequently
                   ! less then zero the PFT should die.
                   pre_indust_ref_gpp(:,j) = ( pre_indust_ref_gpp(:,j) * ( tau_gpp_self_thin - dt ) + &
                        gpp_year(:,j) * dt ) / tau_gpp_self_thin

                ELSEWHERE

                   ! Calculate the reference gpp if gpp < zero. Don't
                   ! use the negative gpp because that could result in the
                   ! reference gpp getting negative. Truncate gpp to zero
                   ! instead. Given that gpp_daily is now set to zero, this
                   ! term can be removed from the calculation.
                   pre_indust_ref_gpp(:,j) = ( pre_indust_ref_gpp(:,j) * ( tau_gpp_self_thin - dt ) ) / &
                        tau_gpp_self_thin
                   
                ENDWHERE
          
             ENDWHERE
             
          END IF

       END DO

    END IF ! calculate_gpp_preind .AND. ts_annual_proc


    !! 14.1 longterm radiation use efficiency
    ! We add an ELSEWHERE condition for the biomass check of longterm radiation use efficiency
    ! If we don't have the biomass, use the previous value. The additional statements are needed
    ! to make sure that rue_longterm always has a value (even during severe droughtstress).
    DO j=2,nvm

       WHERE(SUM(circ_class_biomass(:,j,:,ileaf,icarbon)*circ_class_n(:,j,:),2).GT. min_stomate)

          rue_longterm(:,j) =  ( rue_longterm(:,j) * ( one_year - dt ) + &
               gpp_daily(:,j) / &
               (1. - exp(-0.5 * biomass_to_lai(SUM(circ_class_biomass(:,j,:,ileaf,icarbon)*&
                circ_class_n(:,j,:),2),npts,j))) * dt ) / (one_year)
         
       ELSEWHERE

          rue_longterm(:,j) = rue_longterm(:,j)

       ENDWHERE

    ENDDO

    !! 14.2 compute weekly maintenance respiration. We simply use the tau_gpp_week.
    WHERE ( veget_max .GT. zero )

       resp_maint_week = ( resp_maint_week * ( tau_gpp_week - dt ) + &
            resp_maint_daily * dt ) / tau_gpp_week

    ELSEWHERE

       resp_maint_week = zero

    ENDWHERE

    DO j = 2,nvm ! Loop over # PFTs
       WHERE ( ABS(resp_maint_week(:,j)) .LT. EPSILON(zero) )
          resp_maint_week(:,j) = zero
       ENDWHERE
    ENDDO

    !CYmark: There seems no good reason to calculate resp_maint before 
    ! any truncation is made in stomate_growth_fun_all.f90. Need to think more
    ! carefully on how to properly calculate resp_maint when using its value
    ! to define Ipresenescence. For the moment we just remove any negative value.
    DO j = 2,nvm ! Loop over # PFTs
       WHERE ( resp_maint_week(:,j) .LT. zero )
          resp_maint_week(:,j) = zero
       ENDWHERE
    ENDDO
    !
    !! 15. Update the maximum and minimum moisture availabilities (::maxvegstress_thisyear 
    !!     and ::minvegstress_thisyear). If the daily moisture availability, ::vegstress_day,
    !!     which was calculated earlier in this subroutine, is greater or less than the current
    !!     value of the maximum or minimum moisture availability, then set the value for the maximum
    !!     or minimum to the current daily value respectively. \n
    !!     ::maxvegstress_thisyear and ::minvegstress_thisyear are used to update the variables 
    !!     ::maxvegstress_lastyear and ::minvegstress_lastyear in Section 21 of this subroutine, 
    !!     which are used in the module ::stomate_phenology for the leaf onset models 'hum' and 
    !!     'humgdd', and in the module ::stomate_turnover for the leaf senescence models 'dry' 
    !!     and 'mixed'.
    !

    WHERE ( vegstress_day .GT. maxvegstress_thisyear )
       maxvegstress_thisyear = vegstress_day
    ENDWHERE

    WHERE ( vegstress_day .LT. minvegstress_thisyear )
       minvegstress_thisyear = vegstress_day
    ENDWHERE

    !
    !! 16. Update the annual maximum weekly GPP (::maxgppweek_thisyear), if the weekly GPP
    !!     is greater than the current value of the annual maximum weekly GPP.
    !!     The use of this variable is described in Section 14.
    !

    WHERE ( gpp_week .GT. maxgppweek_thisyear )
       maxgppweek_thisyear = gpp_week
    ENDWHERE

    !
    !! 17. Update the annual GDD0 by adding the current daily 2-meter temperature 
    !!     (multiplied by the stomate time step, which is one day), if it is greater
    !!     than zero degrees C. \n
    !!     This variable is mostly used in the module ::lpj_establish.
    !

    WHERE ( t2m_daily .GT. ZeroCelsius )
       gdd0_thisyear = gdd0_thisyear + dt * ( t2m_daily - ZeroCelsius )
    ENDWHERE

    !
    !! 18. Update the annual precipitation by adding the current daily precipitation
    !!     amount (multiplied by the stomate time step, which is usually one day). \n
    !
    precip_thisyear = precip_thisyear + dt * precip_daily

    !
    !! 18. Update the longterm annual precipitation by adding the current daily precipitation
    !!     amount (multiplied by the stomate time step, which is usually one day). \n
    !
    precip_longterm = (precip_longterm * (tau_longterm - dt) + precip_daily * dt) / tau_longterm
    
    !
    !! 19. Update the annual maximum leaf mass for each PFT (::lm_thisyearmax) and the maximum fractional
    !!     plant cover if the LPJ DGVM is activated.       
    !
    
    !
    !! 19.1 If the LPJ DGVM is activated first the fraction of natural vegetation (::fracnat), i.e. non-
    !!      agricultural vegetation (PFTs 2-11) is calculated for each PFT. Each PFT is looped over.
    !
    IF(ok_dgvm ) THEN

       fracnat(:) = un
       DO j = 2,nvm ! Loop over # PFTs
          IF ( .NOT. natural(j) ) THEN
             fracnat(:) = fracnat(:) - veget_max(:,j)
          ENDIF
       ENDDO

    ENDIF

    !
    !! 19.2 If LPJ and STOMATE are activated, first the maximum fractional plant cover needs to be updated, 
    !!      and then this year's leaf biomass. Each PFT is looped over.
    !!      Both are updated according to the linear relaxation method described above. 
    !!      \latexonly
    !!      \input{season_lin_relax_eqn1.tex}
    !!      \endlatexonly
    !!      \n
    !!      The time constant for this process is set as one year (in seconds) divided by the variable
    !!      ::leaf_age_crit, which gives the leaf lifetime in days and is set for each PFT in the module 
    !!      ::stomate_constants. Each PFT is looped over. \n
    !

    IF ( ok_stomate ) THEN
       IF(ok_dgvm ) THEN
          DO j=2,nvm ! Loop over # PFTs

             !
             !! 19.2.1 Calculate maximum fractional plant cover (::maxfpc_lastyear).
             !!        If natural vegetation is present in the grid cell, and the leaf
             !!        biomass is greater than three-quarters of last year's leaf biomass, the maximum fractional plant
             !!        cover for last year is updated. \n 
             !!        The short-term variable (Xs in the above equation) that is being used to update the long-term 
             !!        maximum fractional plant cover is the fractional cover of natural vegetation, specified as 
             !!        ::veget/::fracnat. Last year's value is then set to be this year's value.
             !

             IF ( natural(j) .AND. ok_dgvm ) THEN

                WHERE ( fracnat(:) .GT. min_stomate .AND. SUM(circ_class_biomass(:,j,:,ileaf,icarbon)*&
                    circ_class_n(:,j,:),2).GT. lm_lastyearmax(:,j)*0.75 )
                   maxfpc_lastyear(:,j) = ( maxfpc_lastyear(:,j) * (one_year/(leaf_age_crit(:,j)/365) - dt ) + &
                        veget(:,j) / fracnat(:) * dt ) / (one_year/(leaf_age_crit(:,j)/365))
                ENDWHERE
                maxfpc_thisyear(:,j) = maxfpc_lastyear(:,j) ! just to initialise value

             ENDIF

!NV : correct initialization
!!$             WHERE(biomass(:,j,ileaf,icarbon).GT. lm_lastyearmax(:,j)*0.75)
!!$                lm_lastyearmax(:,j) = ( lm_lastyearmax(:,j) * ( one_year/(leaf_age_crit(:,j)/365)- dt ) + &
!!$                     biomass(:,j,ileaf,icarbon) * dt ) / (one_year/(leaf_age_crit(:,j)/365))
!!$             ENDWHERE
!!$             lm_thisyearmax(:,j)=lm_lastyearmax(:,j) ! just to initialise value

             
             !
             !! 19.2.2 Update this year's leaf biomass (::lm_thisyearmax).
             !!        The short-term variable (Xs in the above equation) that is being used to update the long-term 
             !!        this year's leaf biomass is the leaf biomass pool (::biomass(i,j,ileaf).
             !
             WHERE (lm_thisyearmax(:,j) .GT. min_stomate)
                WHERE( SUM(circ_class_biomass(:,j,:,ileaf,icarbon)*&
                    circ_class_n(:,j,:),2).GT. lm_thisyearmax(:,j)*0.75)
                   lm_thisyearmax(:,j) = ( lm_thisyearmax(:,j) * ( one_year/(leaf_age_crit(:,j)/365) - dt ) + &
                         SUM(circ_class_biomass(:,j,:,ileaf,icarbon)*&
                         circ_class_n(:,j,:),2) * dt ) / (one_year/(leaf_age_crit(:,j)/365))
                ENDWHERE
             ELSEWHERE
                lm_thisyearmax(:,j) = SUM(circ_class_biomass(:,j,:,ileaf,icarbon)*&
                    circ_class_n(:,j,:),2)
             ENDWHERE

          ENDDO

       ELSE
          
          !
          !! 19.3 If LPJ DGVM is not activated but STOMATE is, the maximum leaf mass is set to be the same  
          !!      as the leaf biomass (::biomass(i,j,ileaf), without a change in the maximum fractional plant cover. 
          !
          DO j = 2,nvm ! Loop over # PFTs
             WHERE (  SUM(circ_class_biomass(:,j,:,ileaf,icarbon)*&
                    circ_class_n(:,j,:),2) .GT. lm_thisyearmax(:,j) )
                lm_thisyearmax(:,j) = SUM(circ_class_biomass(:,j,:,ileaf,icarbon)*&
                    circ_class_n(:,j,:),2)
             ENDWHERE
          ENDDO

       ENDIF !(ok_dgvm)
    ELSE

          !
          !! 19.4 If STOMATE is not activated, the maximum leaf biomass is set to be the maximum possible 
          !!      LAI of the PFT.
          !

       DO j = 2,nvm ! Loop over # PFTs
          lm_thisyearmax(:,j) = lai_to_biomass(lai_max(j),j)
       ENDDO

    ENDIF  !(ok_stomate)

    !
    !! 20. Update the annual maximum fractional plant cover for each PFT if the current fractional cover 
    !!     (::veget), is larger than the current maximum value.
    !!     ::veget is defined as fraction of total ground. Therefore, maxfpc_thisyear has the same unit.
    !

    WHERE ( veget(:,:) .GT. maxfpc_thisyear(:,:) )
       maxfpc_thisyear(:,:) = veget(:,:)
    ENDWHERE
 
    ! 20.X nitrogen stress for carbon allocation 
    IF ( Dmax.GT.0.0 ) THEN 
       DO j=2,nvm 
          nn(:) = SUM(circ_class_biomass(:,j,:,ileaf,initrogen)*&             
                    circ_class_n(:,j,:),2) + SUM(circ_class_biomass(:,j,:,iroot,initrogen)*&
                    circ_class_n(:,j,:),2) + SUM(circ_class_biomass(:,j,:,ilabile,initrogen)*&
                    circ_class_n(:,j,:),2)   
          WHERE (nn(:) .GT. min_stomate)  
             cn(:) = (  SUM(circ_class_biomass(:,j,:,ileaf,icarbon)*&
                    circ_class_n(:,j,:),2) + SUM(circ_class_biomass(:,j,:,iroot,icarbon)*&
                    circ_class_n(:,j,:),2) + SUM(circ_class_biomass(:,j,:,ilabile,icarbon)*&
                    circ_class_n(:,j,:),2) ) / nn(:)   
          ELSEWHERE 
             cn(:) = zero 
          ENDWHERE 
 
          WHERE (cn(:) .GT. min_stomate)  
             nstress_season(:,j)= ( nstress_season(:,j) * ( one_year - dt ) + & 
                                    MAX(MIN( (cn_leaf_init_2D(:,j)/cn(:)),1.),0.3 ) & 
                  * dt ) / one_year  
          ENDWHERE 
       ENDDO 
 
    ELSE 
       nstress_season(:,:)=1.0 
    ENDIF 
       
    CALL histwrite_p (hist_id_stomate, 'NSTRESS_SEASON', itime, &
         nstress_season, npts, horipft_index)
    CALL xios_orchidee_send_field("NSTRESS_SEASON",nstress_season)
    CALL xios_orchidee_send_field("GPP_WEEK",gpp_week)
    CALL xios_orchidee_send_field("MAINT_RESP_WEEK",resp_maint_week)

    !
    ! 2.3 "growing season"
    ! is_growing_season and is_pot_growing_season are are calculated according to 
    ! defferent definition of the growing season. Both are used in xios as temporal 
    ! mask for the calculation of Water stress but maybe also used in NPP, LAI etc ...

    ! xios_default_val is used here because this value is not used in the averaging 
    ! operation by xios so by mutiplying is_growing_season and a variable such as 
    ! vegstress_day in field_def_orchidee.xml we let xios do the calculation and 
    ! release ORCHIDEE from this duty.  
    is_growing_season(:,:) = xios_default_val
    is_pot_growing_season(:,:) = xios_default_val
    gtemp(:,:) = zero

    DO j=2,nvm
       ! same equation as used in stomate_growth_fun_all.f90 for the activity of 
       ! labile carbon pool 
       WHERE ((t2m_week(:)-ZeroCelsius) .GT. tmin_labile(j))
         gtemp(:,j) = EXP((e0_labile(j))*(1.0/(tref_labile(j)-tmin_labile(j)) - &
                      1.0/(t2m_week(:)-ZeroCelsius-tmin_labile(j))))
       ENDWHERE
       IF(is_tree(j))THEN
         WHERE( SUM(circ_class_biomass(:,j,:,ileaf,icarbon)*&
                    circ_class_n(:,j,:),2).GT.min_stomate)
            is_growing_season(:,j) = un
            vegstress_season(:,j) = ( vegstress_season(:,j) * ( longevity_sap(j) - dt ) + &
                        MIN(1.,2.*vegstress_day(:,j)) * dt ) / longevity_sap(j)
         ENDWHERE
       ELSE
         WHERE( SUM(circ_class_biomass(:,j,:,ileaf,icarbon)*&
                    circ_class_n(:,j,:),2).GT.min_stomate)
            is_growing_season(:,j) = un
            vegstress_season(:,j) = ( vegstress_season(:,j) * ( 30.  - dt ) + &
                        MIN(1.,2.*vegstress_day(:,j)) * dt ) / 30.
         ENDWHERE
       ENDIF
       WHERE(gtemp(:,j) .GE. un)
          is_pot_growing_season(:,j) = un
       ENDWHERE
    ENDDO

    CALL xios_orchidee_send_field("GROWING_SEASON",is_growing_season)
    CALL xios_orchidee_send_field("POT_GROWING_SEASON",is_pot_growing_season)

    ! 20.X nitrogen stress for photosynthesis  
    !    DO j=2,nvm 
    !       WHERE (biomass(:,j,ileaf,initrogen) .GT. min_stomate)  
    !          nstress2_season(:,j)= ( nstress2_season(:,j) * ( one_year/(leaf_age_crit(:,j)/365) - dt ) + & 
    !               MAX(1.-(nstress2_daily(:,j)-1.),0.3) *dt )/ (one_year/(leaf_age_crit(:,j)/365) )  
    !       ENDWHERE 
    !    ENDDO        

    !
    !! 21. At the end of the every year, last year's maximum and minimum moisture availability,
    !!     annual GDD0, annual precipitation, annual max weekly GPP, and maximum leaf mass are 
    !!     updated with the value calculated for the current year, and then their value reset.
    !!     Either to zero for the maximum variables, or to ::large_value (set to be 1.E33 in ::
    !!     the module ::stomate_constants) for the minimum variable ::minvegstress_thisyear.
    !!
    !!     Note that ts_annual_proc is true if current time-step corresponds to the moment when the annual
    !!     processes should be done. By default, this corresponds to the last time-step of the year.

    IF ( ts_annual_proc ) THEN

       !
       !! 21.1 Update last year's values. \n
       !!      The variables ::maxvegstress_lastyear, ::minvegstress_lastyear and ::maxgppweek_lastyear  
       !!      are updated using the relaxation method:
       !!      \latexonly
       !!      \input{season_lin_relax_eqn1.tex}
       !!      \endlatexonly
       !!      \n
       !!      where Xs is this year's value, dt (Delta-t) is set to 1 (year), and the time constant (tau) !??
       !!      is set by the parameter ::tau_climatology, which is set to 20 years at the beginning 
       !!      of this subroutine. \n
       !!      The other variables (::gdd0_lastyear, ::precip_lastyear, ::lm_lastyearmax and 
       !!      ::maxfpc_lastyear) are just replaced with this year's value.
       !
       !NVMODIF
       maxvegstress_lastyear(:,:) = (maxvegstress_lastyear(:,:)*(tau_climatology-1)+ maxvegstress_thisyear(:,:))/tau_climatology
       minvegstress_lastyear(:,:) = (minvegstress_lastyear(:,:)*(tau_climatology-1)+ minvegstress_thisyear(:,:))/tau_climatology
       maxgppweek_lastyear(:,:) =( maxgppweek_lastyear(:,:)*(tau_climatology-1)+ maxgppweek_thisyear(:,:))/tau_climatology
       !       maxvegstress_lastyear(:,:) = maxvegstress_thisyear(:,:)
       !       minvegstress_lastyear(:,:) = minvegstress_thisyear(:,:)
       !       maxgppweek_lastyear(:,:) = maxgppweek_thisyear(:,:)
       
       gdd0_lastyear(:) = gdd0_thisyear(:)

       precip_lastyear(:) = precip_thisyear(:)

       lm_lastyearmax(:,:) = lm_thisyearmax(:,:)

       maxfpc_lastyear(:,:) = maxfpc_thisyear(:,:)

       ! Calculate Tseason
       Tseason(:) = zero
       WHERE ( Tseason_length(:) .GT. min_sechiba )
          Tseason(:) = Tseason_tmp(:) / Tseason_length(:)
       ENDWHERE

       !
       !! 21.2 Reset new values for the "this year" variables. \n
       !!      The maximum variables are set to zero and the minimum variable ::minvegstress_thisyear
       !!      is set to ::large_value (set to be 1.E33 in the module ::stomate_constants).
       !

       maxvegstress_thisyear(:,:) = zero
       minvegstress_thisyear(:,:) = large_value

       maxgppweek_thisyear(:,:) = zero

       gdd0_thisyear(:) = zero

       precip_thisyear(:) = zero

       lm_thisyearmax(:,:) = zero

       maxfpc_thisyear(:,:) = zero

       Tseason_tmp(:) = zero
       Tseason_length(:) =zero
       Tmin_spring_time(:,:)=zero

       !
       ! 21.3 Special treatment for maxfpc. !?? 
       !! 21.3 Set the maximum fractional plant cover for non-natural vegetation 
       !!      (i.e. agricultural C3 and C4 - PFT 12 and 13) vegetation for last year to be zero.
       !

       !
       ! 21.3.1 Only take into account natural PFTs
       !
       DO j = 2,nvm ! Loop over # PFTs
          IF ( .NOT. natural(j) ) THEN
             maxfpc_lastyear(:,j) = zero
          ENDIF
       ENDDO

       ! 21.3.2 In Stomate, veget is defined as a fraction of ground, not as a fraction 
       !        of total ground. maxfpc_lastyear will be compared to veget in lpj_light.
       !        Therefore, we have to transform maxfpc_lastyear.


       ! 21.3.3 The sum of the maxfpc_lastyear for natural PFT must not exceed fpc_crit (=.95).
       !        However, it can slightly exceed this value as not all PFTs reach their maximum 
       !        fpc at the same time. Therefore, if sum(maxfpc_lastyear) for the natural PFTs
       !        exceeds fpc_crit, we scale the values of maxfpc_lastyear so that the sum is
       !        fpc_crit.

!!$       ! calculate the sum of maxfpc_lastyear
!!$       sumfpc_nat(:) = zero
!!$       DO j = 2,nvm ! Loop over # PFTs
!!$          sumfpc_nat(:) = sumfpc_nat(:) + maxfpc_lastyear(:,j)
!!$       ENDDO
!!$
!!$       ! scale so that the new sum is fpc_crit
!!$       DO j = 2,nvm ! Loop over # PFTs 
!!$          WHERE ( sumfpc_nat(:) .GT. fpc_crit )
!!$             maxfpc_lastyear(:,j) = maxfpc_lastyear(:,j) * (fpc_crit/sumfpc_nat(:))
!!$          ENDWHERE
!!$       ENDDO

       ! 3. Calculate legacies for bark beetle suceptibilities
       !! Legacies that determine the vulnerability of bark beetle outbreaks
       ! Drop the last year, keep the other years but move them
       ! One postion in the array

         tmp_legacy2(:,:,:) = zero
         IF(legacy_years_wood .GT. un) THEN
            DO iyear = 1, legacy_years_wood - 1
               ! year time shift   
               tmp_legacy2(:,:,iyear+1) = wood_leftover_legacy(:,:,iyear)
            ENDDO
         ENDIF
         ! Add this year on the first position
         tmp_legacy2(:,:,1) = woody_litter_to_use(:,:)
         ! Update 
         wood_leftover_legacy(:,:,:) = tmp_legacy2(:,:,:)
         ! Reset the yearly variable
         woody_litter_to_use(:,:) = zero
         ! for windtrow damage

         tmp_legacy(:,:,:) = 1
        ! Drop the last year, keep the other years but move them
         ! One postion in the array
         DO iyear = 1, legacy_years-1
           ! year time shift   
           tmp_legacy(:,:,iyear+1) = season_drought_legacy(:,:,iyear)
         ENDDO ! end of iyear loop
         ! Add this year on the first position
         DO j = 2,nvm ! Loop over # PFTs
           ! this equation and it parameter values are coming from
           ! Temperli et al. 2013
           IF(beetle_pft(j)) THEN
             WHERE(SUM(SUM(circ_class_biomass(:, j, :, :, icarbon), 3)*circ_class_n(:, j, :), 2) .GT. min_stomate)
               tmp_legacy(:,j,1)= vegstress_season(:,j)
             ENDWHERE
           ELSE
             tmp_legacy(:,j,1) = 1
           ENDIF
         ENDDO
         tmp_legacy(:,:,1) = vegstress_season(:,:)       ! Update
         season_drought_legacy(:,:,:) = tmp_legacy(:,:,:)

       ! for beetle generation index
         tmp_legacy(:,:,:) = 1
         DO iyear = 1, legacy_years-1
           tmp_legacy(:,:,iyear+1) = beetle_generation_index(:,:,iyear)
         ENDDO ! end of iyear loop

         DO j = 2,nvm ! Loop over # PFTs
            ! this equation and it parameter values are coming from
            ! Temperli et al. 2013
            IF(beetle_pft(j)) THEN

              WHERE(SUM(SUM(circ_class_biomass(:, j, :, :, icarbon), 3) * circ_class_n(:, j, :), 2) .GT. min_stomate)
                tmp_legacy(:,j,1)= 1/(1+exp(-beetle_generation_a(j)*(sumTeff(:,j)/ &
                    beetle_generation_b(j)-beetle_generation_c(j))))
              ENDWHERE
            ELSE
              tmp_legacy(:,j,1) = zero
            ENDIF
         ENDDO

         IF(printlev_loc .GE. 4) THEN
            WRITE(numout,*) 'sumTeff',sumTeff(test_grid,:)
            WRITE(numout,*) 'beetle_generation_legacy',tmp_legacy(test_grid,:,1)
         ENDIF
         ! Update and reset for next year.
         beetle_generation_index(:,:,:) = tmp_legacy(:,:,:)
         sumTeff(:,:) = zero
         beetle_diapause(:,:) = 0


    ENDIF  ! ts_annual_proc

    !
    !! 22. Diagnose herbivore activity (::herbivores), determined as probability for a leaf to be
    !!     eaten in a day (follows McNaughton et al., 1989). \n
    !!     The amount of herbivore activity is used in the modules ::lpj_establish  
    !!     and ::stomate_turnover.
    !

    !
    !! 22.1 First calculate the mean long-term leaf NPP in grid box, mean residence
    !!      time (years) of green tissue to give the biomass available for herbivore consumption.
    !!      Each PFT is looped over, though the calculation is only made for
    !!      natural vegetation (PFTs 2-11).
    !

    nlflong_nat(:,:) = zero
    weighttot(:,:) = zero
    green_age(:,:) = zero
    !
    DO j = 2,nvm ! Loop over # PFTs
       !
       IF ( natural(j) ) THEN
          !
          !! 22.1.1 Calculate the total weight of the leaves (::weighttot) as last year's leaf biomass
          !
          weighttot(:,j) = lm_lastyearmax(:,j)
          
          !
          !! 22.1.2 Calculate the mean long-term leaf NPP as the long-term NPP calculated in Section 12 of 
          !!        this subroutine weighted by the leaf fraction (::leaf_frac), which is defined to be 0.33
          !!        at the beginning of this subroutine.
          !
          nlflong_nat(:,j) = npp_longterm(:,j) * leaf_frac_hvc
          
          !
          !! 22.1.3 Calculate the mean residence time of the green tissue (::green_age) 
          !!        This is calculated as the sum of 6 months 
          !!        for natural seasonal vegetation (i.e. PFTs 3, 6, 8-11), and 2 years for evergreen 
          !!        (PFTs 2, 4, 5, 8), multiplied by last year's leaf biomass for each PFT, divided by the 
          !!        total weight of the leaves for all PFTs.
          !         This is a crude approximation.  !!?? By whom?
          !!        The difference between seasonal and evergreen vegetation is determined by the parameter
          !!        ::pheno_model, which specifies the onset model of the PFT.
 
          !
          IF ( pheno_model(j) .EQ. 'none' ) THEN
             green_age(:,j) = green_age_ever * lm_lastyearmax(:,j)
          ELSE
             green_age(:,j) = green_age_dec * lm_lastyearmax(:,j)
          ENDIF
          !
       ENDIF
       !
    ENDDO
    !
    WHERE ( weighttot(:,:) .GT. min_sechiba )
       green_age(:,:) = green_age(:,:) / weighttot(:,:)
    ELSEWHERE
       green_age(:,:) = un
    ENDWHERE

    !
    !! 22.2 McNaughton et al. (1989) give herbivore consumption as a function of mean, long-term leaf NPP.
    !!      as it gives an estimate of the edible biomass. The consumption of biomass by herbivores and the 
    !!      resultant herbivore activity are calculated following the equations:
    !!      \latexonly
    !!      \input{season_consumption_eqn4.tex}
    !!      \endlatexonly
    !!      \n 
    !!      and 
    !!      \latexonly
    !!      \input{season_herbivore_eqn5.tex}
    !!      \endlatexonly
    !!      \n 
    !       

    DO j = 2,nvm ! Loop over # PFTs
       !
       IF ( natural(j) ) THEN
          !
          WHERE ( nlflong_nat(:,j) .GT. zero )
             consumption(:) = hvc1 * nlflong_nat(:,j) ** hvc2
             herbivores(:,j) = one_year * green_age(:,j) * nlflong_nat(:,j) / consumption(:)
          ELSEWHERE
             herbivores(:,j) = 100000.
          ENDWHERE
          !
       ELSE
          !
          herbivores(:,j) = 100000.
          !
       ENDIF
       !
    ENDDO
    herbivores(:,ibare_sechiba) = zero

    !! 23. Update plant_status
    ! ibudbreak was kept for one day. It can now be
    ! upgraded to the next phenological stage, i.e.,
    ! icanopy. If this is not done the model will
    ! try to grow new leaves where there are already 
    ! leaves
    WHERE ( plant_status(:,:) .EQ. ibudbreak )

       plant_status(:,:) = icanopy

    ENDWHERE

    ! Calculate cn_leaf_min_season only the day after the buds break (hence the 2 in the IF-statement).
    WHERE((when_growthinit(:,:) .EQ. 2) .AND. (SUM(circ_class_biomass(:,:,:,ileaf,initrogen), DIM=3) .GT. min_stomate))
       cn_leaf_min_season(:,:) = SUM((circ_class_biomass(:,:,:,ileaf,icarbon)*circ_class_n(:,:,:)), DIM=3) / &
            SUM((circ_class_biomass(:,:,:,ileaf,initrogen)*circ_class_n(:,:,:)), DIM=3)
    ENDWHERE
    
    !! 24. Calculate Nesterov Index used for open vegetation fire simulation.

    ! Calculate Nesterov Index following Eq.(5) in Thonicke et al. (2010); 
    ! Thonicke et al. (2010) required the cumulaitve summing condition as
    ! precipitation>3mm,
    ! Not sure where the second limiation below comes from.
    vartmp(:) = zero
    WHERE ( (precip_daily(:) .LE. 3.) .AND. ( (t2m_min_daily(:)-4.) .GE. ZeroCelsius) )
      vartmp(:) = (t2m_max_daily(:)-ZeroCelsius)*( t2m_max_daily(:) - (t2m_min_daily(:)-4.) )
    ENDWHERE

    WHERE (vartmp(:) .GT. zero)
      ni_acc(:) = ni_acc(:) + vartmp(:)
    ELSEWHERE
      ni_acc(:) = zero
    ENDWHERE

    IF (printlev>=4) WRITE(numout,*) 'Leaving season_pre_disturbance'
        
  END SUBROUTINE season_pre_disturbance


!! ================================================================================================================================
!! SUBROUTINE   : season_post_disturbance
!!
!>\BRIEF          This subroutine calculates long-term variables related to harvest, 
!!                wind, pest and fire disturbances.
!!
!! DESCRIPTION    This module is called in stomate_lpj after all disturbances have been 
!!                calculated to ensure that the disturbance variables are updated
!!                before the end of the year and that the updated variables are stored
!!                in the restart files
!!
!! RECENT CHANGE(S): None
!!                
!! MAIN OUTPUT VARIABLE(S):  
!!      
!! REFERENCES   : 
!! - Chen et al 2018. Geoscientific model development.
!! - Temperli et al 2013. 
!!                       
!! FLOWCHART    : None
!!   
!_ =================================================================================================================================

  SUBROUTINE season_post_disturbance (npts, dt, t2m_daily, tau_longterm, &
       circ_class_biomass, litter, circ_class_n, & 
       gap_area_save, sumTeff, beetle_diapause, n_reserve_balance, &
       n_reserve_longterm, doy_start_gs, doy_end_gs, mean_start_gs, &
       valid_start_gs, total_ba_init, veget_max, woody_litter_by_cut, &
       woody_litter_to_use)



    !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                                   :: npts                    !! Domain size - number of grid cells (unitless)
    REAL(r_std), INTENT(in)                                      :: dt                      !! time step in days (dt_days)
    REAL(r_std), INTENT(in)                                      :: tau_longterm
    REAL(r_std), DIMENSION(:), INTENT(in)                        :: t2m_daily               !! Daily 2 meter temperatures (K)
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(in)                :: circ_class_biomass      !! Biomass per circumference class 
                                                                                            !! @tex ($gC m^{-2} of ground$) @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)                    :: circ_class_n            !! Number of indiviudals per circumference class
    REAL(r_std),DIMENSION(:,:,:,:,:), INTENT(in)                 :: litter                  !! Above and below ground metabolic and structural
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)                    :: woody_litter_by_cut     !! Saved woody litter by icut not by iparts
    REAL(r_std), DIMENSION(:,:), INTENT(inout)                   :: woody_litter_to_use     !! woody litter pool to use for wood_left_over in the pest module
    REAL(r_std), DIMENSION(:,:), INTENT(inout)                   :: sumTeff                 !! Yearly thermal sum that is effective for bark beetles development
    INTEGER(i_std), DIMENSION(:,:), INTENT(inout)                :: beetle_diapause         !!
 
    REAL(r_std), DIMENSION(:,:), INTENT(in)                      :: n_reserve_balance       !! Actual to potential N reserve pool (unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(in)                      :: total_ba_init           !! basal area calculated at the start of the year
                                                                                            !! (m^{2}/m^{2})
    REAL(r_std), DIMENSION(:,:), INTENT(in)                      :: veget_max               !! "Maximal" coverage fraction of a PFT (LAI -> infinity)
                                                                                            !! on ground
    LOGICAL, DIMENSION(:,:), INTENT(in)                          :: valid_start_gs          !! Doy was calculated by ORCHIDEE. It is therefore valid 
                                                                                            !! and should be used to update mean_start_gs
    
    !! 0.3 Modified variables
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)                 :: gap_area_save           !! 5 (legacy_years_wind) years of the gap area created by 
                                                                                            !! more than 30% basal area loss
                                                                                            !! Used in windthrow calculations
    REAL(r_std), DIMENSION(:,:), INTENT(inout)                   :: n_reserve_longterm      !! "longer term" actual to potential  N reserve pool
                                                                                            !! (0-1, unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)                   :: doy_start_gs            !! growing season starting day of year (DOY) for 
                                                                                            !! deciduous PFTs.
    REAL(r_std), DIMENSION(:,:), INTENT(inout)                   :: doy_end_gs              !! growing season end day of year (DOY) for 
                                                                                            !! deciduous PFTs.
    REAL(r_std), DIMENSION(:,:), INTENT(inout)                   :: mean_start_gs           !! mean growing season starting day for 
                                                                                            !! deciduous PFTs.
    
    !! 0.4 Local variables
    INTEGER(i_std)                                               :: ipts,j,iyear            !! Indices (unitless)
    CHARACTER(30)                                                :: var_name
    REAL(r_std)                                                  :: Teff                    !! Intermediate calculation in order to estimate the 
                                                                                            !! number of beetle generation
    REAL(r_std), DIMENSION(npts)                                 :: tbark_daily             !! Temperature inside the bark
    REAL(r_std)                                                  :: total_ba_loss           !! Basal area loss compared to the start of the year 
                                                                                            !! (m^{2}/m^{2})

!_ =================================================================================================================================

    IF (firstcall_season_post_disturbance) THEN

       ! Initialize local printlev
      ! Guillaume M. -- get_printlev calls getin_p, a COLLECTIVE bcast. Called here,
      ! inside the compute loop under a local firstcall, it hangs the run as soon as the
      ! ranks stop reaching it together. Use printlev, read once at startup by every rank.
       printlev_loc=printlev
       ! Reset firstcall
       firstcall_season_post_disturbance = .FALSE.

    END IF

    IF (printlev_loc>=3) WRITE(numout,*) 'Entering season_post_disturbance'

    ! 1. Bark beetle phenology in response to temperature 
    ! Sum of temperature use in the calculation of the number of beetle
    ! generation all the parameter value are coming from Temperli et al. 2013
    ! To account that inside the bark, the temperature is 2°C highe
    
    tbark_daily(:)= t2m_daily(:) - ZeroCelsius + 2.
    DO ipts = 1, npts
         DO j = 2,nvm ! Loop over # PFTs
            IF (beetle_pft(j))THEN
               IF ( tbark_daily(ipts) .GT. min_temp_beetle(j) .AND. tbark_daily(ipts) &
                      .LT. max_temp_beetle(j) .AND. beetle_diapause(ipts,j) == 0) THEN
    
                  Teff= (opt_temp_beetle(j) - min_temp_beetle(j))*exp(eff_temp_beetle_a(j)* &
                        tbark_daily(ipts)) - exp(eff_temp_beetle_a(j)* &
                        eff_temp_beetle_b(j)-(eff_temp_beetle_b(j) - tbark_daily(ipts))/ &
                        eff_temp_beetle_c(j)) - eff_temp_beetle_d(j)
                  sumTeff(ipts,j)= sumTeff(ipts,j) + Teff
               ENDIF
            ELSE
              sumTeff(ipts,j)= zero
            ENDIF
            IF(printlev_loc .GE. 4 .AND. ipts == test_grid .AND. j == test_pft) THEN
               WRITE(numout,*) 'Teff',Teff
               WRITE(numout,*) 'SumTeff',SumTeff(ipts,j)
               WRITE(numout,*) 'tbark_daily',tbark_daily(ipts)
            ENDIF
         ENDDO
    ENDDO
    ! Daily sum of litter to use for the pest wood_leftover
    ! Woody litter was collected during mortality processes (stomate_kill,
    ! sapience_kill,sapience_lcc) and here the woody litter pool to use is selected.
    woody_litter_to_use(:,:) = woody_litter_to_use(:,:) + &
          woody_litter_by_cut(:,:,icut_storm_break) +  woody_litter_by_cut(:,:,icut_storm_uproot) + &
          woody_litter_by_cut(:,:,icut_lcc_wood)
    IF(printlev_loc .GE. 4) THEN
         WRITE(numout,*) 'post_season'
         WRITE(numout,*) 'woody_litter_to_use',woody_litter_to_use(:,4)
         DO ipts = 1,ncut_times
            WRITE(numout,*) 'icut, woody_litter_by_cut',ipts, woody_litter_by_cut(:,4,ipts)
         ENDDO
    ENDIF

    IF (ts_annual_proc) THEN

       ! 2. Calculate the gap created by basal area loss
       ! If trees are exposed to the newly formed edge, the trees get more
       ! susceptible to windthrow. To apply this, the windthrow module creates
       ! virtual gap area to calculate areas for trees closer to the edge. Here, we
       ! assume that if the basal area loss is more than 30% (gap_threshold), this loss
       ! is intensive enough to create trees exposed to new edges. After legacy_years_wind
       ! (5 years now) we consider that the trees are acclimated to the edge.
 
       ! Move values for gap_area_save with a one-year time shift. 
       ! Drop the oldest value (t + legacy_years_wind) and replace it with the value for 
       ! t + (legacy_years_wind - 1). Move t0 to t+1. Empty the t0 and will be filled.
       IF (ok_windthrow) THEN
          DO iyear = legacy_years_wind,2,-1
             gap_area_save(:,:,iyear) = gap_area_save(:,:,iyear-1)
          ENDDO
          gap_area_save(:,:,1) = zero

          ! Calculate basal area difference at the end of the year. 
          DO j = 2,nvm
             IF ( is_tree(j) ) THEN
                DO ipts = 1,npts
                   total_ba_loss = total_ba_init(ipts,j) -  &
                        SUM(wood_to_ba(circ_class_biomass(ipts,j,:,:,icarbon),j, &
                        pipe_tune2(ipts,j)) * circ_class_n(ipts,j,:))

                   ! If the difference is larger than tha gap_threshold (rate) of the
                   ! initial basal area then save the basal area loss to gap_area_save
                   ! Currently gap_threshold is set to 0.3 since thinning is
                   ! considered as strong if it is over 30% of the stand
                   IF( total_ba_loss .GE. gap_threshold * total_ba_init(ipts,j) .AND. &
                        total_ba_init(ipts,j) .GT. min_stomate) THEN
                      gap_area_save(ipts,j,1) = total_ba_loss/ &
                         total_ba_init(ipts,j) * veget_max(ipts,j) * area(ipts)
                   ENDIF

                   IF (SUM(gap_area_save(ipts,j,:)) .GT. area(ipts)* &
                                                    veget_max(ipts,j) )THEN
                        ! Think should not happen. need to stop the model?
                   ENDIF
                ENDDO
             ENDIF
          ENDDO

          ! AED_FEEDBACK (Marie 2026): the storm feed to dA_storm MOVED to
          ! stomate_windthrow (wind_damage), where it now uses the STORM-specific
          ! area_damaged_storm (cc_kill_to_area on icut_storm_break+uproot) instead of
          ! gap_area_save. gap_area_save was total >30% basal-area loss (all causes) and
          ! double-counted harvest + background mortality. gap_area_save is kept for its
          ! own role (windthrow edge-vulnerability legacy, legacy_years_wind buffer).
       ENDIF ! ok_windthrow
    END IF

    ! 4. Nitrogen reserves for root_uptake
    ! Update the "long term" n_reserve_balance. 
    ! (::n_reserve_balance in /day, ::n_reserve_longterm in /year.)
    ! The time constant is given by the parameter ::tau_longterm, 
    ! which is set in ::stomate_data to be 3 * one year (in seconds, as 
    ! described above). This calculation should be done after stomate_lpj
    ! so it is placed in the post_disturbance. 
    ! n_reserve_longterm is used in ::nitrogen_dynamics to calculate
    ! feedback of N uptake based on the longterm ratio of actual N reserve 
    ! to N reserve target. This calculation is placed post_disturbance
    ! since n_reserve_balance is calculated in stomate_lpj.
    WHERE (veget_max(:,:).GT.min_stomate) 
       n_reserve_longterm(:,:) = (n_reserve_longterm(:,:) * (tau_longterm - dt ) + &
            n_reserve_balance(:,:) * dt ) / tau_longterm
    END WHERE

    ! 5. Phenology
    ! Calculating the 7 year average budbreak days
    ! Start the average with real values instead of zeros. If not it will take 
    ! much longer before the average values represent real start doy
    WHERE (mean_start_gs(:,:) .EQ. zero)
       mean_start_gs(:,:) = doy_start_gs(:,:)
    END WHERE

    ! doy_start_gs gets one value per year. The other 364 days their
    ! value is zero. On those days we don't want to update the average. Note
    ! that tau_week is used. This means that we will average over seven years instead
    ! of seven days because this calculation is done only once per year. Only use 
    ! calculated values to update the mean. If the mean was used as a fail safe option
    ! for phenology or senescence this DOY should NOT be used to update the mean.
    WHERE (mean_start_gs(:,:).GT.zero .AND. &
         doy_start_gs(:,:).GT.zero .AND. &
         valid_start_gs(:,:))

       mean_start_gs(:,:) = (mean_start_gs(:,:) * (tau_week - dt) + &
            doy_start_gs(:,:) * dt ) / tau_week

    END WHERE

    IF (printlev>=4) WRITE(numout,*) 'Leaving season_post_disturbance'

  END SUBROUTINE season_post_disturbance

END MODULE stomate_season
