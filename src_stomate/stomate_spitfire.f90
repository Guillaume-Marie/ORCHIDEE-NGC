!  ==============================================================================================================================
!  MODULE		 			: stomate_spitfire
!
!  CONTACT		 			: orchidee-help _at_ ipsl.jussieu.fr
!
!  LICENCE	 			        : IPSL (2006)
!  This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF				       This module simulates fires in
!! burnable vegetations in a prognostic way. 
!! 
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/DOC/ORCHIDEE/src_sechiba/enerbil.f90 $
!! $Date: 2012-00-00 12:22:59 +0100  $
!! $Revision: ??? $
!! \n
!_ ================================================================================================================================

!!---------------------------------------------------------------------------------
!! Below are some remarks that further document model behavior or existing minor
!! issues in the model, they're minor in nature and do not affect large-scale
!! model performance.
!! ----------------------- Remarks ------------------------------------------------
!! Date 2023/05/09
! 1. There are few places where clarifications are needed, mainmy minor, they're
! marked as ??, or [UNCLEAR]
! 2. Sources of some parameters are not clear and are marked as "SOURCE"
! 3. Remaining questions are marked as QUESTION
!! --------------------------------------------------------------------------------


MODULE stomate_spitfire

  ! modules used:
  USE netcdf
  USE ioipsl_para
  USE constantes
  USE grid
  USE xios_orchidee
  USE stomate_data
  USE pft_parameters
  USE pft_parameters_var
  USE interpol_help,        ONLY: aggregate_p

  USE function_library, ONLY: wood_to_dia, wood_to_height, wood_to_cn_dia, &
                              check_mass_balance, cc_kill_to_area

  IMPLICIT NONE

  ! private & public routines
  PRIVATE
  PUBLIC spitfire, spitfire_clear, spitfire_month_input, spitfire_annual_input

  ! first call
  LOGICAL, SAVE             :: firstcall_spitfire = .TRUE.
!$OMP THREADPRIVATE(firstcall_spitfire)  
  INTEGER(i_std), SAVE          :: printlev_loc                  !! Local level of text output for current module
!$OMP THREADPRIVATE(printlev_loc)

CONTAINS

!! ================================================================================================================================
!! SUBROUTINE    : spitfire_clear
!!
!>\BRIEF        Set the flag ::firstcall_sapiens_forestry to .TRUE. and as such activate section 2 of the subroutine forestry (see below).
!! 
!_ ================================================================================================================================

  SUBROUTINE spitfire_clear
    firstcall_spitfire = .TRUE.
  END SUBROUTINE spitfire_clear

!! ================================================================================================================================
!! SUBROUTINE     : spitfire
!!
!>\BRIEF          A prognostic fire module. It uses daily temperature and
!! precipitation to calculate a fire danger index, which is combined with ignition sources to
!! derive number of fires. Simulation of fire size is based on a fire spread model under
!! the assumption of equilibrium combustion state. Then burned area is derived by combining
!! fire size and the number of fires. Fuel consumption and tree crown
!! scorching, tree mortality were last calculated. 
!!    Note that the module is fully empirical in a sense that it describes fire spread and intensity
!! corresponding to certain conditions of fuel moisture, fuel compaction, and the amount of fuel being
!! burned, rather than explicitly simulating the energy processes during combustion such as 
!! condutance and convective heating of unburned fuels, and radiative and convective dissipation of
!! energy released during combustion.
!!
!! DESCRIPTION	: 
!!
!! 1. Calculate fire danger index from daily temperature and precipitation using
!! Nesterov Index (NI). The fuel mositure and moisture of extinction are compared
!! to derive a daily fire danger index (FDI). Fuel moiture is calculated empirically using
!! Nesterov Index by weighting among differnet types of fuels. The fuel type with
!! a bigger surface-area-to-volume ratio receives less weight in determining the
!! fuel moisture state.
!! 2. Use daily FDI and ignition source data to derive daily fire numbers.
!! 3. Apply Rothermel's equation to calculate fire spread rate for both forward and backward
!! directions. Rate of spread is determined in a way responding to fire
!! reaction intensity(+), propagating flux ratio(+), wind speed(+), fuel bulk
!! density(-),effective heating number(-), and heat of pre-ignition(-).
!! 4. Using daily FDI to derive fire duration time, combined with rate of spread to
!! derive the mean fire size. Combine fire size and fire numbers to derive burned area.
!! 5. Fuel mositure calculated during the calculation of daily FDI is used to
!! determin consumption fraction of different types of fuels.
!! 6. Calculate fire frontline intensity using forward fire spread rate, fuel consumption, 
!! and burned area, thus further to determin flame height and fire damage to tree crown. 
!! 7. Fire caused tree mortality is determined based on two considerations, crown damage due to 
!! flame combustion, and damage of fire residence time to tree barks. Crowning scorching results in
!! live biomass combustion as emissions, while unburned dead tree biomass is transfered to litter pool.
!! 
!! In summary, fire weather (temperature, preciptation, wind speed) plays a central 
!! role in this module, it determines fuel moisture, which futher determines
!! daily fire danger index (FDI) --> nubmer of fires --> fire duration --> fire
!! size --> burned area --> fire intensity, residence time --> tree demage
!! 
!! Important Note:
!! 1. With a single day and and within a single pixel, fire size is considered as uniform, hence
!! burned ara = fire size X fire number. 
!! 2. No multi-day fire duration is considered. Fire duration consisting of multiple days is 
!! hence simulated as independent fires within consecutive days as long as fuel moisture conditions
!! allows meaningful burning (i.e., with a reasonalbe frontline intensity).
!! 3. Most of the equations are two from documents: Thonicke et al. (2010) and a supplemenary 
!! SPITFIRE technical document. When not otherwise mentioned, all cited Equations are from
!! Thonicke et al. (2010).
!! 
!! RECENT CHANGE(S): None
!!
!! REFERENCES    :
!! - K.Thonicke, A.Spessa, I.C. Prentice, S.P.Harrison, L.Dong, and
!! C.Carmona-Moreno, 2010: The influence of vegetation, fire spread and fire
!! behaviour on biomass burning and trace gas emissions: results from a
!! process-based model. Biogeosciences, 7, 1991-2011.
!! - Rothermel, R. C.: A Mathematical Model for Predicting Fire
!! Spread in Wildland Fuels, Intermountain Forest and Range Experiment 
!! Station, Forest Service, US Dept. of Agriculture, Ogden, UtahUSDA 
!! Forest Service General Technical Report INT-115, 1972.
!! - Pyne, S. J., Andrews, P. L., and Laven, R. D.: Introduction to wildland
!! fire, 2 edition, Wiley, New York, 769 pp., 1996
!!
!! FLOWCHART     : None
!! \n 
!_ ================================================================================================================================


   SUBROUTINE spitfire(npts, dt_days, lalo, veget, veget_max, resolution, contfrac,   &
                       wind_speed_actual, soilhum_daily, ni_acc, forest_managed,         &
                       vpd_daily_mean, vpd_daily_max, vpd_mean_week, vpd_max_week, precip_month,    &
                       popd, a_nd, edge_length, lightning, circ_class_biomass, circ_class_n,       &
                       circ_class_kill, bm_to_litter, litter, litterfuel,        &
                       emissions_fire)


    
    ! Local parameters
    INTEGER, PARAMETER     :: ntrace=6              !! Number of different trace gas types included
    REAL(r_std), PARAMETER :: MINER_TOT=0.055       !! Mineral content of fuel mass. Source: Thonicke et al., 2010. p1010 Table A1 Row 5
    REAL(r_std), PARAMETER :: H=18000.0             !! Heat content of fuel mass (kJ kg^{-1}). Source: Thonicke et al. (2010) Appendix A
    REAL(r_std), PARAMETER :: sigma_1hr=66.0        !! Surface-area-to-volume ratio (cm^{2} cm^{-3}). Source:  (USDA 1978, Harvey et al. 1997);
                                                    !! These two sources are given in a working doc by K. Thonicke, the exact sources unknown.
    REAL(r_std), PARAMETER :: sigma_10hr=3.58       !!
    REAL(r_std), PARAMETER :: sigma_100hr=0.98      !!
    REAL(r_std), PARAMETER :: me_livegrass=0.2      !! Moisture of extinction for live grass fuel (unitless) 
    REAL(r_std), PARAMETER :: carbon_to_drymass= 2  !! Assuming a fraction of 0.5 for carbon to dry mass

    ! [SOURCE] this two parameters are said to follow Brown et al. 1981 but I don't
    ! have the paper.
    REAL(r_std), PARAMETER :: fbd_a = 1.2           !! Ratio of fuel bulk density of 10hr fuels to 1hr fuel. (unitless) 
    REAL(r_std), PARAMETER :: fbd_b = 1.4           !! Ratio of fuel bulk density of 100hr fuels to 1hr fuel. (unitless)

    ! Input variables
    INTEGER, INTENT(in)                             :: npts                        !! Domain size
    REAL(r_std), INTENT(in)                         :: dt_days                     !! Time step of stomate (days)
    REAL(r_std),DIMENSION(:,:), INTENT(in)          :: veget                       !! Fractional coverage: actually share of the pixel 
                                                                                   !! covered by a PFT taking into account LAI
    REAL(r_std), DIMENSION(:,:)                     :: veget_max                   !! "Maximal" coverage fraction of a PFT (LAI -> infinity) on ground.
    REAL(r_std), DIMENSION(npts,2), INTENT(in)      :: lalo                        !! Vector of latitude and longitudes (beware of the order !)
    REAL(r_std), DIMENSION(:,:,:,:,:),INTENT(in)    :: circ_class_biomass          !! Biomass per circumference class @tex $(gC tree^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)       :: circ_class_n                !! Number of individuals in each circ class
                                                                                   !! @tex $(m^{-2})$ @endtex
    INTEGER(i_std), DIMENSION(:,:), INTENT(in)      :: forest_managed              !! forest management flag (is the forest
                                                                                   !! being managed?)
    REAL(r_std), DIMENSION(npts,2), INTENT(in)      :: resolution                  !! Resolution at each grid point in m (1=E-W, 2=N-S)
    REAL(r_std),DIMENSION (npts), INTENT (in)       :: contfrac                    !! Fraction of continent in the grid
    REAL(r_std), DIMENSION(npts), INTENT(in)        :: wind_speed_actual            !! Daily wind speed (m s^{-1})
    REAL(r_std), DIMENSION(npts), INTENT(in)        :: lightning                   !! Number of daily lightning (flashes km^{-2} day^{-1}) 
    
    REAL(r_std),DIMENSION (npts)       :: baadjust_ratio              !! Scaling factors used to adjust the simulated burned area. (unitless)
    !REAL(r_std),DIMENSION (npts), INTENT (in)       :: baadjust_ratio              !! Scaling factors used to adjust the simulated burned area. (unitless)
    !                                                                               !! This factor is used to scale simulated burned area to fit satellite observations. 
    REAL(r_std), DIMENSION(npts), INTENT(in)        :: popd                        !! Human population density (person km^{-2}) 
    REAL(r_std), DIMENSION(npts), INTENT(in)        :: a_nd                        !! Parameter for potential human-caused ignitions; Used in Eq. (3)
    REAL(r_std), DIMENSION(npts), INTENT(in)        :: edge_length                 !! Edge length from road observations (m) (urban road should be excluded)
    REAL(r_std), DIMENSION(npts), INTENT(in)        :: soilhum_daily               !! Daily top layer soil humidity, here we need only the top soil 
                                                                                   !! layer humdity, and this is addressed by 
                                                                                   !! passing soilhum_daily(:,1) as soilhum_daily when we 
                                                                                   !! call SPITFIRE in stomate_lpj.f90
    REAL(r_std), DIMENSION(npts),INTENT(in)         :: ni_acc                      !! Nesterov index (square of degree Celcius)
    REAL(r_std), DIMENSION(:),INTENT(in)            :: vpd_daily_mean              !! Daily mean VPD
    REAL(r_std), DIMENSION(:),INTENT(in)            :: vpd_daily_max               !! Daily max VPD
    REAL(r_std), DIMENSION(:),INTENT(in)            :: vpd_mean_week
    REAL(r_std), DIMENSION(:),INTENT(in)            :: vpd_max_week
    REAL(r_std), DIMENSION(npts),INTENT(in)         :: precip_month                !! Monthly precipitation

    ! Modified variables
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)           :: litter           !! Metabolic and structural litter, above and belowground. (gC m^{-2})
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)           :: litterfuel       !! Dead litter fuel above ground. (gC m^{-2})
    REAL(r_std), DIMENSION(:,:,:,:),INTENT(inout)              :: emissions_fire   !! Emissions from fire. (gC m^{-2})
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)             :: bm_to_litter     !! Conversion of biomass to litter (gC m^{-2} dtslow^{-1})
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)           :: circ_class_kill  !! Number of trees within a circ class that needs
                                                                                   !! to be killed @tex $(ind m^{-2})$ @endtex


    ! Local variables
    REAL(r_std), DIMENSION(npts)                    :: dfm_livegrass             !! Daily fuel moisture for live grass. (unitless) (0,1)
    REAL(r_std), DIMENSION(npts)                    :: dfm                       !! Daily fuel moisture of dead fuels on ground excluding 1000hr fuel. (0,1)
    REAL(r_std), DIMENSION(npts)                    :: ignition_efficiency       !! Ignition efficiency as a function of fuel load, (0,1), unitless.
    REAL(r_std), DIMENSION(npts)                    :: dfm_1hr                   !! Daily fuel moisture for 1hr-fuel. (unitless) (0,1)
    REAL(r_std), DIMENSION(npts)                    :: vd                        !! Vegetation density for calculating vpd-based fire danger index
    REAL(r_std), DIMENSION(npts,nhour,nvm)          :: litterfuel_class_pft      !! PFT-level aboveground litter fuel of different fuel classes
                                                                                 !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)                :: litterfuel_pft            !! PFT-level aboveground litter fuel by integrating all fuel classes.
                                                                                 !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nhour)              :: litterfuel_class_grid     !! Gridcell value of aboveground litter fuel of different fuel classes
                                                                                 !! @tex $(gC m^{-2})$ @endtex per ground area
    REAL(r_std), DIMENSION(npts)                    :: litterfuel_grid           !! Total aboveground litter fuel over the grid cell, based on per area of
                                                                                 !! of burnable PFTs.
                                                                                 !! @tex $(gC m^{-2})$ @endtex per ground area
    REAL(r_std), DIMENSION(npts)                    :: litterfuel_fine           !! Sum of 1hr/10hr/100hr litter fuel for all burnable PFTs. (gC m^{-2})
    REAL(r_std), DIMENSION(npts,nvm)                :: livegrass_pft             !! PFT-level live grass biomass.
                                                                                 !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts)                    :: livegrass_grid            !! Area-weighted sum of live aboveground biomass of burnable herbaceous PFTs. (gC m^{-2})
    REAL(r_std), DIMENSION(npts)                    :: netfuel_grid              !! Total net fuel (mineral content corrected) by summing the dead fine fuel 
                                                                                 !! and live grass. Note that 1000h dead fuel is excluded.
                                                                                 !! (kg dry mass m^{-2}).
    REAL(r_std), DIMENSION(npts)                    :: me_litterfuel_grid        !! Grid level moisture of extinction weighted by the 
                                                                                 !! dead groud litter fuel of different burnable PFTs
                                                                                 !! @tex $(unitless)$ @endtex
    REAL(r_std), DIMENSION(npts)                    :: me_grid                   !! The moisture of extinction weighted by ground dead fuel and livegrass,
                                                                                 !! m_{e} in Eq.(8). (unitless; should be understood as 
                                                                                 !! mass_water/(mass_water + mass_dry_fuel))
    REAL(r_std), DIMENSION(npts,nvm)                :: bulk_density_pft          !! PFT-level fuel bulk density. (kg dry mass m^{-3})
    REAL(r_std), DIMENSION(npts)                    :: bulkdensity_deadfinefuel  !! Fuel bulk density for dead fine fuel weighted by mass of different PFTs.
                                                                                 !! (kg dry mass m^{-3}) 
    REAL(r_std), DIMENSION(npts)                    :: bulkdensity_livegrass     !! Fuel bulk density for live grass weighted by mass of different burnable grass PFTs.
                                                                                 !! (kg dry mass m^{-3}) 
    REAL(r_std), DIMENSION(npts)                    :: bulkdensity_grid          !! Fuel bulk density weighted by mass of dead fine fuel and live grass.
                                                                                 !! (kg dry mass m^{-3}) 
                                                                                 !! including only burnable PFTs, in unit of kg dry mass per m3. 
    REAL(r_std), DIMENSION(npts)                    :: fdi                       !! Fire danger index. (unitless) (0-1)
    REAL(r_std), DIMENSION(npts)                    :: fdi_niacc
    REAL(r_std), DIMENSION(npts)                    :: fdi_vpd_daily_mean
    REAL(r_std), DIMENSION(npts)                    :: fdi_vpd_daily_max
    REAL(r_std), DIMENSION(npts)                    :: fdi_vpd_mean_week
    REAL(r_std), DIMENSION(npts)                    :: fdi_vpd_max_week
    REAL(r_std), DIMENSION(npts)                    :: wind_speed                !! U_{forward} in Eq. (A5) (i.e., wind speed adjusted by fractions of 
                                                                                 !! herbaceous and tree covers). (m min^{-1})
    REAL(r_std), DIMENSION(npts)                    :: wind_forward              !! It represents 3.281*U_{forward} in Eq. (A5). (m min^{-1})
    REAL(r_std), DIMENSION(npts)                    :: ros_f                     !! Forward fire spred rate (m min^{-1})
    REAL(r_std), DIMENSION(npts)                    :: ros_b                     !! Backward fire spread rate (m min^{-1})
    REAL(r_std), DIMENSION(npts)                    :: lb                        !! Length-to-breadth ratio of the fire ellipse, unitless, range[1,8]
    REAL(r_std), DIMENSION(npts)                    :: df                        !! Diameter at the forward direction for the fire ellipse (m)
    REAL(r_std), DIMENSION(npts)                    :: db                        !! Diameter at the backward direction for the fire ellipse (m)
    REAL(r_std), DIMENSION(npts)                    :: fire_durat                !! Fire duration (minutes) 
    REAL(r_std), DIMENSION(npts)                    :: sigma                     !! Surface-area-to-volume ratio weighted by mass for 1h, 10h and 100 fuel types.
                                                                                 !! (cm^{-1}, given that area in unit of cm^{2} and volume in cm^{3})
    REAL(r_std), DIMENSION(npts)                    :: fpc_tree_total            !! total tree fraction
    REAL(r_std), DIMENSION(npts)                    :: fpc_grass_total           !! total grass fraction (including burnable herbaceous PFTs)
    REAL(r_std), DIMENSION(npts)                    :: wetness                   !! Ratio of dead fuel moisture to the moisture of extinction. (unitless)
    REAL(r_std), DIMENSION(npts)                    :: frontline_intensity       !! Surface fire frontline intensity (kW m^{-1})
                                                                                 !! this variable is PFT vcmax weighted sum of above one
    REAL(r_std), DIMENSION(npts)                    :: area_burnt_orig           !! Daily burned area before adjustment by using empirical scaling coefficient (ha).
    REAL(r_std), DIMENSION(npts)                    :: area_burnt                !! Daily burned area (ha)
    REAL(r_std), DIMENSION(npts)                    :: numfire                   !! Number of fires (1)
    REAL(r_std), DIMENSION(npts)                    :: numfire_human             !! Number of fires by human ignitions (1)
    REAL(r_std), DIMENSION(npts)                    :: numfire_lightning         !! Number of fires by lightning ignitions (1)
    REAL(r_std), DIMENSION(npts)                    :: fire_frac                 !! Fraction of area burned against burnable ground area only. (unitless)
    REAL(r_std), DIMENSION(npts)                    :: fire_frac_grid            !! Fraction of area burned against grid cell land area. (unitless)
    REAL(r_std), DIMENSION(npts)                    :: area_land                 !! Land area within gridcel excluding water body (ha, i.e., 10000 m2)
    REAL(r_std), DIMENSION(npts)                    :: area_burnable_veg         !! Land area covered with burnable vegetation within gridcell (ha)
    REAL(r_std), DIMENSION(npts)                    :: human_ign                 !! Human ignitions (1 day^{-1} km^{-2})
    REAL(r_std), DIMENSION(npts)                    :: lightn_ign                !! Lightning ignitions (1 day^{-1} km^{-2})
    REAL(r_std), DIMENSION(npts)                    :: lightning_cg_ratio        !! Fractions of cloud-to-ground lightnings of total lightnings (unitless)
    REAL(r_std), DIMENSION(npts)                    :: human_suppression         !! Huamn suppression effects on lightning ignitions (unitless) (0-1)
    REAL(r_std), DIMENSION(npts,nhour)              :: cf_litterfuel             !! Combustion fractions for different surface fuel categories. (unitless) (0,1)
    REAL(r_std), DIMENSION(npts)                    :: cf_finefuel               !! Combustion fraction for surface fine litter fuel, including 1hr/10hr/100hr fuel.
                                                                                 !! (unitless) (0-1)
    REAL(r_std), DIMENSION(npts)                    :: cf_lg                     !! Combustion fraction for livegrass biomass (unitless) (0-1)
    REAL(r_std), DIMENSION(npts)                    :: var_gamma                 !! A variable used in calculating fire residence time (tau_l) 
                                                                                 !! for postfire mortality.
    REAL(r_std), DIMENSION(npts,nvm)                :: tau_l                     !! Residence time of fire (minutes)
    REAL(r_std), DIMENSION(npts,nvm,ncirc)          :: tau_c                     !! Critical time for cambial damage (minutes)
    REAL(r_std), DIMENSION(npts,nvm,ncirc)          :: diameter                  !! Diameter at breast height for representative tree (m)
    REAL(r_std), DIMENSION(npts,nvm,ncirc)          :: bark_thickess             !! bark thickness (cm)
    REAL(r_std), DIMENSION(npts,nvm,ncirc)          :: canopy_height             !! Tree height (m)
    REAL(r_std), DIMENSION(npts,nvm,ncirc)          :: crown_length              !! Crown length, i.e. crown vertical diameter (m)
    REAL(r_std), DIMENSION(npts,nvm,ncirc)          :: crown_width               !! Crown width, i.e. crown horizontal diameter (m)
    REAL(r_std), DIMENSION(npts,nvm,ncirc)          :: crown_base_height         !! Height of the base of tree crown (m)
    REAL(r_std), DIMENSION(npts,nvm)                :: scorching_height          !! Fire flame scorching height (m)
    REAL(r_std), DIMENSION(npts,nvm,ncirc)          :: ck                        !! Proportion of crown scorched by fire. (unitless) (0-1)
    REAL(r_std), DIMENSION(npts,nvm,ncirc)          :: pm_tau                    !! Tree mortality due to cambial damage (unitless) (0-1)
    REAL(r_std), DIMENSION(npts,nvm,ncirc)          :: pm_ck                     !! Tree mortality due to crown damage (unitless) (0-1)
    REAL(r_std), DIMENSION(npts,nvm,ncirc)          :: postf_mort                !! Fire-caused forest mortality. (unitless) (0-1)

    REAL(r_std), DIMENSION(npts,nvm)                :: dcflux_fire_pft           !! fire carbon flux to atmosphere; including both crown emssions 
                                                                                 !! and litter consumption emissions
    REAL(r_std), DIMENSION(npts,nvm,ntrace)         :: emissions_trace_gas       !! Trace gas emissions from fire (g m^{-2})
    REAL(r_std), DIMENSION(npts)                    :: mean_fire_size_original   !! mean fire size before the correction of surface fire intensity
    REAL(r_std), DIMENSION(npts)                    :: mean_fire_size            !! mean fire size
    REAL(r_std), DIMENSION(npts)                    :: burnable_vegfrac          !! Fraction of the grid cell that are occupied by burnable vegetation.
    REAL(r_std), DIMENSION(npts,nvm)                :: relative_cover            !! Relative ground coverage for a given PFT against the total ground area of burnable
                                                                                 !! PFTs only. (unitless) (0-1)
    REAL(r_std), DIMENSION(npts,nlitt,nvm,nhour,nelements) :: litterfuel_consumed !! Litter fuel consumed during fire on per PFT area basis. (gC m^{-2})
    REAL(r_std), DIMENSION(npts,nlitt,nvm,nhour,nelements) :: litterfuel_fc_ba    !! Litter fuel consumed during fire on the area basis of burned arae only (gC m^{-2}).
    REAL(r_std), DIMENSION(npts)                    :: finefuel_fc_ba            !! Litter fuel consumed during fire on the area basis of burned arae only (gC m^{-2}).
    REAL(r_std), DIMENSION(nvm,ntrace)              :: ef_trace                  !! Emission factors for trace gases (g (kg dry mass)^{-1})
    REAL(r_std), DIMENSION(npts)                    :: vartmp                    !!
    REAL(r_std), DIMENSION(npts,nvm)                :: vartmp_pft                !!
    REAL(r_std), DIMENSION(npts,nvm)                :: area_damaged_fire         !! Forest area killed by fire, cc-kill basal-area route (m2)
    INTEGER :: ivm,x,m,k,n,i,icir,ipts,ihour,iele,ilit,ipar,ilev,imbc
    INTEGER :: itestpft,itestipts

                                                                                 !! of each month.
    REAL(r_std)                                     :: fuel_low_bound=200        !! Lower bound of total fuel load below which ignition efficiency is 0 
                                                                                 !! @tex $gC m^{-2}$ @endtex 
    REAL(r_std)                                     :: fuel_high_bound=1000      !! Upper bound of total fule load above which ignition efficiency is 1
                                                                                 !! @tex $gC m^{-2}$ @endtex 
                                                                                 !! Sources: Arora & Boer (2005)

    REAL(r_std)                                     :: surface_threshold         !! A threshold in frontline intensity beyond which tree crown scorching
                                                                                 !! is considered.
    REAL(r_std), DIMENSION(npts,nvm)                :: disturb_crown             !! a temporary variable, used to store the C emissions 
                                                                                 !! from crown scorch of live tree crown.

    REAL(r_std), DIMENSION(npts)                    :: fc_crown                  !! PFT weighted carbon emissions caused by crown scorching (gCm^{-2}).
                                                                                 !! only for history file writing.
    REAL(r_std), DIMENSION(npts,nvm,nmbcomp,nelements)&
                                                    :: check_intern              !! Contains the components of the internal
                                                                                 !! mass balance chech for this routine
                                                                                 !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)      :: closure_intern            !! Check closure of internal mass balance
                                                                                 !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)      :: pool_start                !! Start pool of this routine 
                                                                                 !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)      :: pool_end                  !! End pool of this routine 
                                                                                 !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    CHARACTER(LEN=8), DIMENSION(nhour)              :: fueltype_str               !! string suffix indicating fuel category

    ! AED (Average Edge Distance, m) is now provided by stomate_data as AED_fire,
    ! computed by stomate_data.calculate_aed (see stomate.f90, called per-day).
    REAL(r_std), DIMENSION(npts)                    :: fedge_wind                 !! Fraction of edge area affected by wind (unitless) (0-1)
    REAL(r_std), DIMENSION(npts)                    :: max_patch_fire_size        !! Maximum fire size
    REAL(r_std), DIMENSION(npts,nvm)                :: max_patch_fire_size_PFT    !! Maximum fire size for each PFT
    REAL(r_std), DIMENSION(npts)                    :: Patch_Area                 !! Patch area for each circle
    REAL(r_std), DIMENSION(npts,nvm)                :: grass_biomass_wet          !! total wet mass of grass biomass (fuel) in (tons/Ha)
    LOGICAL, DIMENSION(npts,nvm)                    :: crown_fire                 !! Crown fire flag
    REAL(r_std), DIMENSION(npts,nvm)                :: ignited_weight             !! The biomass of trees that suffer from crown fire (gC/m2)
    REAL(r_std), DIMENSION(npts,nvm)                :: total_weight               !! The biomass of trees (gC/m2)
    REAL(r_std), DIMENSION(npts,nvm)                :: ignited_fraction           !! The Fraction of trees that suffer from crown fire (unitless) (0-1)
    REAL(r_std)                                     :: grass_wind_red_factor=0.6  !! Grass wind reduction factor (unitless) (0-1)
    REAL(r_std)                                     :: tree_wind_red_factor=0.4   !! Forest wind reduction factor (unitless) (0-1)
    REAL(r_std), DIMENSION(npts)                    :: mean_fire_size_temp        !! mean fire size
!_ ===============================================================================================================================

    IF (firstcall_spitfire) THEN   
       ! Initialize local printlev
      ! Guillaume M. -- get_printlev calls getin_p, a COLLECTIVE bcast. Called here,
      ! inside the compute loop under a local firstcall, it hangs the run as soon as the
      ! ranks stop reaching it together. Use printlev, read once at startup by every rank.
       printlev_loc=printlev
       firstcall_spitfire=.FALSE.
    END IF

    IF (printlev .GE. 2) WRITE(numout,*) 'Entering spitfire module'
    
    ! Initialise variables
    litterfuel_consumed(:,:,:,:,:) = zero
    bulkdensity_grid(:) = zero
    disturb_crown(:,:)=zero
    fc_crown(:)=zero
    surface_threshold = zero
    litterfuel_class_grid(:,:) =zero
    area_burnable_veg(:) = zero


    !!++TEMP++
    baadjust_ratio(:) = 1.
    !! Set values for some variables.
    ! We assume 0.03 of lightning flashes input has the power to ignite a fire.
    ! This was introduced in Yue et al. (2014) GMD Sect. 2.1.1 on Page 2749 following
    ! Prentice et al. (2011).
    lightning_cg_ratio(:)=0.03
    CALL calculate_ground_lightn_ratio(npts,lalo,lightning_cg_ratio)

    ! Assign emission factors for trace gas emissions which convert dry mass emissions
    ! into trace gas emissions. This was maintained from Thonicke et al. (2010) but 
    ! might be subjected to revisions by people interested fire air chemistry.
    DO ivm=1,nvm
       ef_trace(ivm,1)=ef_CO2(ivm)
       ef_trace(ivm,2)=ef_CO(ivm)
       ef_trace(ivm,3)=ef_CH4(ivm)
       ef_trace(ivm,4)=ef_VOC(ivm)
       ef_trace(ivm,5)=ef_TPM(ivm)
       ef_trace(ivm,6)=ef_NOx(ivm)
    ENDDO
    ! Initialise the test PFT
    itestipts=1
    itestpft=2
    IF (printlev_loc .GE. 4) THEN
      WRITE(numout,*) 'Debuging for PFT: ',itestpft, 'pixel: ',itestipts
      WRITE(numout,*) 'veget_max: ', veget_max(itestipts,itestpft)
      DO iele=1,nelements
        WRITE(numout,*) 'Element: ', iele
        ! output only the compartment biomass by summing over the axis of circ class
        WRITE(numout,*) 'Biomass: ', SUM(circ_class_biomass(itestipts,itestpft,:,:,iele),DIM=1)
        WRITE(numout,*) 'litter: ', litter(itestipts,:,itestpft,iabove,iele)
        ! output by summing ove the axis of fuel category
        WRITE(numout,*) 'litterfuel: ', SUM(litterfuel(itestipts,:,itestpft,:,iele),DIM=2)
      ENDDO
    ENDIF


 !! 1. Initialize check for mass balance closure
  
    !! 1.2 Initialize check for mass balance closure 
    IF (err_act .GT. 1) THEN

       pool_start(:,:,:) = zero
       DO iele = 1,nelements

          DO ipar = 1,nparts
             DO icir = 1,ncirc
                !  Initial biomass
                pool_start(:,:,iele) = pool_start(:,:,iele) + &
                     (circ_class_biomass(:,:,icir,ipar,iele) * &
                     circ_class_n(:,:,icir) * veget_max(:,:)) 
             ENDDO
          ENDDO

          ! The litter carbon pools need to be included because fires can burn aboveground litter.
          ! Litter pool (gC m-2) *  (m2 m-2) 
          DO ilit = 1,nlitt
             DO ilev = 1,nlevs
                pool_start(:,:,iele) = pool_start(:,:,iele) + &
                     litter(:,ilit,:,ilev,iele) * veget_max(:,:)
             ENDDO
          ENDDO

       ENDDO ! iele

    ENDIF ! err_act.GT.1

    !! 0. Prepare some variables for the simulation of fires.

    ! We separate dead fuel into different classes according to their scale of hours needed
    ! to reach moisture equilibrium with the ambient moisture, 1h, 10h, 100h and 1000h. But
    ! only the first three classes (1h, 10h, and 100h, collectively called as fine fuel) 
    ! contribute to determining fire behavior. This is
    ! reflected in the calculation of a few key variables that determine fire behavior.
    ! These variables include: surface-area-to-volume ratio; fuel bulk density; fuel load.

    ! Note:
    ! (1) For all variables with the name ending with '_grid', if they measure the mass
    ! on per area basis, then the area is based on per burnable ground area.

    ! 0.1 calculate the fraction of area occupied by burnable PFTs.
    burnable_vegfrac(:) = zero
    DO ivm = 1,nvm
      IF (burnable(ivm)) THEN 
           burnable_vegfrac(:)=burnable_vegfrac(:)+veget_max(:,ivm)
      ENDIF
    ENDDO

    ! Calculate pixel area.
    area_land(:) = area(:)*contfrac(:)/10000.0  ! unit:hectar (10000 m2)
    ! [QUESTION] Not sure how to deal with nonbio
    ! Chao: 2021-12-29: this relies on a critical assumption, nonbiofrac=0. 
    ! This variable provides the maximal area that could be burned.
    area_burnable_veg(:) = area_land(:)*burnable_vegfrac(:)  !unit:hectar (10000 m2)

    ! AED_fire is now computed in stomate_data.calculate_aed (called from stomate
    ! once per day, ahead of spitfire). The local recomputation, range check, and
    ! ipslerr_p fatal guard have been moved to that centralised routine to allow
    ! re-use by the canopy LIGHT edge effect (stomate_laieff).
    IF (printlev_loc .GE. 5 .AND. &
        (ok_aed_humign .OR. ok_aed_fuel .OR. ok_aed_size .OR. ok_aed_wind)) THEN
      WRITE(numout,*) 'AED_fire ', AED_fire(:)
    ENDIF

    relative_cover(:,:) = zero
    DO ipts = 1,npts
      IF (burnable_vegfrac(ipts) .GT. min_stomate)  THEN
        DO ivm = 1,nvm
          IF (burnable(ivm)) THEN 
               relative_cover(ipts,ivm)=veget_max(ipts,ivm)/burnable_vegfrac(ipts)
          ENDIF
        ENDDO
      ENDIF
    ENDDO
 
    ! Calculate different variables of aboveground litter fuel 
    litterfuel_class_pft(:,:,:) = zero
    DO ivm=1,nvm
      IF (burnable(ivm)) THEN 
        DO ihour=1,nhour
          litterfuel_class_pft(:,ihour,ivm) = SUM(litterfuel(:,:,ivm,ihour,icarbon),DIM=2)
        ENDDO
      ENDIF
    ENDDO
    litterfuel_pft(:,:) = SUM(litterfuel_class_pft(:,:,:),DIM=2)

    ! Note that litte fuel density should be measured by ground surface that's burnable.
    DO ivm=1,nvm
      IF (burnable(ivm)) THEN 
        DO ihour=1,nhour
          litterfuel_class_grid(:,ihour) = litterfuel_class_grid(:,ihour) + &
             litterfuel_class_pft(:,ihour,ivm)*relative_cover(:,ivm)
        ENDDO
      ENDIF
    ENDDO

    litterfuel_grid(:) = SUM(litterfuel_class_grid,DIM=2)
    litterfuel_fine(:) = SUM(litterfuel_class_grid(:,ihour1:ihour100),DIM=2)

    ! 0.2 Calculate the moisture of extinction weighted by the fuel load of different burnable PFTs.
    ! me: moisture of extinction; me is defined in constantes_mtc.f90 
    ! as 'flammability threshold' and is PFT-specific
    me_litterfuel_grid(:)=zero  !PFT weighted moisture of extinction.
    DO ipts=1,npts
      IF (litterfuel_grid(ipts) .GT. min_stomate) THEN
        DO ivm=1,nvm
          IF (burnable(ivm)) THEN
            me_litterfuel_grid(ipts) = me_litterfuel_grid(ipts) + me(ivm)* &
              litterfuel_pft(ipts,ivm)*relative_cover(ipts,ivm)/litterfuel_grid(ipts)
          ENDIF
        ENDDO
      ENDIF
    ENDDO

    ! Calculate tree and herbaceous foliar projective covers respectively. The fraction of pasture
    ! could be included, depending on how the 'burnable' variable is parameterized for pasture PFT.
    fpc_tree_total(:)=zero
    fpc_grass_total(:)=zero
    DO ivm=1,nvm
      IF (burnable(ivm)) THEN
        IF (is_tree(ivm)) THEN
          fpc_tree_total(:)=fpc_tree_total(:)+relative_cover(:,ivm)
        ELSE
          fpc_grass_total(:)=fpc_grass_total(:)+relative_cover(:,ivm) 
        ENDIF
      ENDIF
    ENDDO

    ! 0.3 Introduce reduction factor for wind speed with respect to forest or grass ground coverage.
    ! wind_speed is $U_{forward}$ in Thonicke et al. (2010).
    ! Note that the unit of input wind speed (wind_speed_actual) should be $(m s^{-1})$, 
    ! so that ROS is in unit of $(m min^{-1})$.
    ! Here we mulitply by 60.0 to convert the unit of wind speed from meter-per-second to 
    ! meter-per-minute. 0.4 and 0.6 are the wind speed reduction factor before the average 
    ! as weighted by foliar projective cover is calculated. (see texts right
    ! above Eq.10 in Thonicke et al. 2010). But I am not sure about its theoretical basis.

    ! For the edge depth of wind infltration, we set the value at 16 m. Effective EDW IND is actually 4 m, 
    ! since at any time in any patch, we assume the wind can only come from one of four idealised wind
    ! directions so that at the grid-scale average edge depth for wind is 4 m. 
    ! The edge depth for wind (edge_distance_wind=16.0 m) is changable.
    fedge_wind(:) = 1.0-((AED_fire(:)-edge_distance_wind/4)**2)/(AED_fire(:)**2)
    wind_speed(:)=(fpc_tree_total(:)*wind_speed_actual(:)*60.0*tree_wind_red_factor)+ &
    (fpc_grass_total(:)*wind_speed_actual(:)*60.0*grass_wind_red_factor)

    IF (ok_aed_wind) THEN
      IF (printlev_loc .GE. 5) THEN
        WRITE(numout,*) 'wind_speed before:', wind_speed(:)
        WRITE(numout,*) 'wind_speed multiplier', edge_wind_factor*(1+fedge_wind(:))
      ENDIF
      wind_speed(:)=(fpc_tree_total(:)*wind_speed_actual(:)*60.0*tree_wind_red_factor*edge_wind_factor*(1+fedge_wind(:)))+ &
      (fpc_grass_total(:)*wind_speed_actual(:)*60.0*grass_wind_red_factor)
      IF (printlev_loc .GE. 5) THEN
        WRITE(numout,*) 'wind_speed after:', wind_speed(:)
      ENDIF
    ENDIF


    wind_forward(:)=3.281*wind_speed(:) !wind_forward is used in Eq(A5) as 3.281*Uforward     

    ! 0.4 Calculate the length-to-breadth ratio of the elliptical fire scar.
    ! the minimum value of lb is 1 and the maximum value is 8
    ! [UNCLEAR] Notice that for Eq.(12) and Eq.(13) in Thonicke 2010, 
    ! 0.06*wind_speed is used here as U_{forward}. Compared with original Eq.(12),
    ! an extra 0.06 was multiplied. Not yet know why this 0.06 is used, it's
    ! in the original code given to me by Patricia.
    ! [UNCLEAR] Not sure where this 16.67 threshold comes from
    WHERE (wind_speed(:) .LT. 16.67)
      lb(:) = 1    
    ELSEWHERE 
      lb(:) = MIN( 8.0, fpc_tree_total(:) * &
                       (1.0+(8.729*((1.0-(EXP(-0.03*0.06*wind_speed(:))))**2.155))) + &
                       (fpc_grass_total(:)*(1.1+((0.06*wind_speed(:))**0.0464))) )
    ENDWHERE

    ! -- We explain the method to calculate an integrated fuel bulk density --
    ! To calculate fuel bulk density for the whole grid cell, we need to first calculate bulk density
    ! for different burnable PFTs, and then calculate for the whole grid cell. 
    ! Derivation of the equation used to calculate bulk density for a given PFT is below:
    ! We have the fuel mass for 1h, 10h and 100h fuel as M1,M10,M100, bulk density for 1h fuel as rho_1,
    ! The ratios of bulk density for 10h and 100h fuel to 1h fuel as `a`,`b`, then
    ! MeanDensity = (M1+M10+M100)/( M1/rho_1 + M10/(rho_1*a) + M100/(rho_1*b) )
    !             = rho_1*(M1+M10+M100) / (M1+M10/a+M100/b)
    ! 
    ! Grid-level fuel bulk density that includes all burnable PFTs could be derived similarly,
    ! but `M1`,`M10`, etc. in this case should be `M_{ivm}`, based on per ground area rather than
    ! per PFT area. The total amount of litter fuel including all PFTs should also be based on per ground area.  

    ! bulkdensity_deadfuel is PFT-specific bulk density for 1h fuel and defined in constantes_mtc.f90, with the
    ! unit of $(kg dry mass m^{-3})$. Note here all variables of litter fuel is defined in unit of carbon,
    ! but it has no influence on the calculation process if we assume a constant coefficient to 
    ! convert carbon to dry mass. 
    bulk_density_pft(:,:) = 1.0
    vartmp(:) = zero
    DO ivm=1,nvm
      IF (burnable(ivm)) THEN
        WHERE (litterfuel_pft(:,ivm) .GT. min_stomate)
          bulk_density_pft(:,ivm) = litterfuel_pft(:,ivm)/ &
             ( litterfuel_class_pft(:,ihour1,ivm) + litterfuel_class_pft(:,ihour10,ivm)/fbd_a + &
               litterfuel_class_pft(:,ihour100,ivm)/fbd_b ) * bulkdensity_deadfuel(ivm)
         
          vartmp(:) = vartmp(:) + (litterfuel_pft(:,ivm)*relative_cover(:,ivm))/bulk_density_pft(:,ivm)
        ENDWHERE
      ENDIF
    ENDDO

    bulkdensity_deadfinefuel(:)=zero
    WHERE(litterfuel_fine(:) .GT. min_stomate)
      bulkdensity_deadfinefuel(:) = litterfuel_fine(:)/vartmp(:)
    ENDWHERE


    ! Calculate livegrass biomass. I assume that live grass biomass contains leaf,
    ! aboveground sapwood and heartwood and fruit.
    livegrass_pft(:,:) = zero
    DO ivm=1,nvm
      IF ( burnable(ivm) .AND. (.NOT. is_tree(ivm)) ) THEN
        livegrass_pft(:,ivm) = &
          SUM( (circ_class_biomass(:,ivm,:,ileaf,icarbon) + &
                circ_class_biomass(:,ivm,:,isapabove,icarbon) + &
                circ_class_biomass(:,ivm,:,iheartabove,icarbon) + &
                circ_class_biomass(:,ivm,:,ifruit,icarbon)) *circ_class_n(:,ivm,:),DIM=2 )
      ENDIF
    ENDDO
    livegrass_grid(:) = SUM( livegrass_pft(:,:)*relative_cover(:,:),DIM=2 )

    
    ! Calculate bulk density for live grass fuel.
    ! bulkdensity_livefuel is PFT-dependent live grass fuel bulk density and could be parameterized.
    ! For the equation used please refer to the calculation of litter fuel bulk density above.
    vartmp(:) = zero
    DO ivm=1,nvm
      IF ( burnable(ivm) .AND. (.NOT. is_tree(ivm)) ) THEN
        WHERE (livegrass_pft(:,ivm) .GT. min_stomate)
          vartmp(:) = vartmp(:) + (livegrass_pft(:,ivm)*relative_cover(:,ivm))/bulkdensity_livefuel(ivm)
        ENDWHERE
      ENDIF
    ENDDO

    bulkdensity_livegrass(:)=zero
    WHERE(livegrass_grid(:) .GT. min_stomate)
      bulkdensity_livegrass(:) = livegrass_grid(:)/vartmp(:)
    ENDWHERE

    ! Calculate grid cell level bulk density.
    bulkdensity_grid(:) = zero
    vartmp(:) = litterfuel_fine(:) + livegrass_grid(:)
    DO ipts=1,npts
      IF (bulkdensity_deadfinefuel(ipts) .LT. min_stomate) THEN
        IF (bulkdensity_livegrass(ipts) .LT. min_stomate) THEN
          bulkdensity_grid(ipts) = zero     ! Both bulk densities are too small. 
                      ! We expect very little fuel on this grid.
                      ! Thus we set bulkdensity_grid as zero.
        ELSE
          bulkdensity_grid(ipts) = bulkdensity_livegrass(ipts)
        ENDIF

      ELSE
        IF (bulkdensity_livegrass(ipts) .LT. min_stomate) THEN
          bulkdensity_grid(ipts) = bulkdensity_deadfinefuel(ipts)
        ELSE
          bulkdensity_grid(ipts) = vartmp(ipts)/( litterfuel_fine(ipts)/bulkdensity_deadfinefuel(ipts) + &
                                                  livegrass_grid(ipts)/bulkdensity_livegrass(ipts) ) 
          
        ENDIF

      ENDIF 
    ENDDO
    
    ! Calculate deal fuel and livegrass weighted average moisture of extinction.
    ! me_grid is m_{e} in Eq(8); moisture of extinction. 
    WHERE ( (litterfuel_fine(:) .GT. min_stomate) .OR. (livegrass_grid(:) .GT. min_stomate) )
      vartmp(:) = livegrass_grid(:) / ( litterfuel_fine(:)+livegrass_grid(:) )
      me_grid(:) = me_livegrass * vartmp(:) + me_litterfuel_grid(:) * (1-vartmp(:))
    ELSEWHERE
      me_grid(:) = zero
    ENDWHERE

    !! 1. Calculate daily fire burned area


    !! Burnable vegetation density calculation
    vd(:) = zero
    DO ivm=1,nvm
      IF (burnable(ivm) ) THEN
        vd(:) = vd(:) + veget(:,ivm)
      ENDIF
    ENDDO
    
    !! 1.1 Calculate fire danger index
    CALL fire_danger_index(npts, nvm, veget, me_grid, litterfuel_class_grid, ni_acc, &
                           vpd_daily_mean,vpd_daily_max,vpd_mean_week,vpd_max_week, precip_month, vd,   &
                           dfm, dfm_1hr, fdi_niacc, fdi_vpd_daily_mean,fdi_vpd_daily_max, fdi_vpd_mean_week, &
                           fdi_vpd_max_week)
    IF (test_fdi_types) THEN
        ! Use the vpd-based fire danger index    
        fdi(:) = fdi_vpd_max_week(:)
    ELSE 
        ! Use the original temperature-based fire danger index    
        fdi(:) = fdi_niacc(:)
    ENDIF
    
    
    ! Calculate human ignitions
    CALL calculate_human_ignition(npts,popd,a_nd,AED_fire,human_ign)

    ! Apply the ignition efficiency as adjusted by fuel load.
    ! This was introduced by Yue et al. (2014) GMD Page 2750.
    IF (printlev_loc .GE. 5) THEN
      WRITE(numout,*) 'Before ignition_efficiency'
      WRITE(numout,*) 'litterfuel_grid',litterfuel_grid(:)
      WRITE(numout,*) 'livegrass_grid',livegrass_grid(:)
    ENDIF
    vartmp(:) = litterfuel_grid(:) + livegrass_grid(:)
    WHERE (vartmp(:) .LE. fuel_low_bound) 
      ignition_efficiency(:) = zero
    ELSEWHERE (vartmp(:) .GT. fuel_low_bound .AND. vartmp(:) .LT. fuel_high_bound)
      ignition_efficiency(:) = (vartmp(:) - fuel_low_bound)/(fuel_high_bound - fuel_low_bound)
    ELSEWHERE
      ignition_efficiency(:) = 1.0
    ENDWHERE

    !! 1.2 Calculate number of fires from lightning and human ignitions

    ! Dervie potential fire number by lightnings. We also account for human suppression in 
    ! lightning ignitions. Note that human suppression in human-ignited fires is implicit.
    CALL human_suppression_func(npts,popd,human_suppression)
    lightn_ign(:)=lightning(:)*lightning_cg_ratio(:)*(1-human_suppression(:))

    ! Here we apply Eq.(2), fdi: unitless; lightning ignition: 1/day/km**2 (as in the original nc file); 
    ! human_ign: unit as 1/day/km**2
    ! 0.01*area_land: grid cell total land area in km**2; 
    ! final unit: 1/day in the concered grid cell 

    numfire_lightning(:)=fdi(:)*ignition_efficiency(:)* lightn_ign(:) * (0.01 * area_burnable_veg(:))  
    numfire_human(:)=fdi(:)*ignition_efficiency(:)* human_ign(:) * (0.01 * area_burnable_veg(:))
    IF (ok_lightningfire) THEN
       ! Account for both the human and natural fires
       numfire(:)=numfire_lightning(:) + numfire_human(:)
    ELSE
       ! Ignore the fires from lightening
       numfire(:)=  numfire_human(:)
    ENDIF

    !! 1.3 Calculate daily area burnt

    ! Calculate net fuel load on the groud. Note that only fine fule and livegrass contribute to
    ! fire propogation.
    netfuel_grid(:) = (litterfuel_fine(:)+livegrass_grid(:)) * carbon_to_drymass &
                      *(1.0-MINER_TOT)/1000.0 !  kg dry mass per squre meter

    ! Calculate mass weighted surface-area-to-volume ratio for fine ground litter fuel.
    WHERE (litterfuel_fine(:).GT.min_stomate)
      sigma(:) = ( litterfuel_class_grid(:,ihour1) * sigma_1hr + &
                   litterfuel_class_grid(:,ihour10) * sigma_10hr + &
                   litterfuel_class_grid(:,ihour100) * sigma_100hr ) / litterfuel_fine(:)
    ELSEWHERE
      sigma(:)=0.00001
    ENDWHERE

    ! Calculate fire spread rate.
    ! Output variables: ros_f(ROS_{f,surface} in Eq.9); var_gamma, wetness (dfm/me).
    CALL rate_of_spread(npts, netfuel_grid, wind_forward, sigma,   &
                        dfm, me_grid, H, bulkdensity_grid,         &
                        ros_f, var_gamma, wetness)

    ! Calculate backward surface fire spread rate ROS_{b,surface} in Eq(10); 
    ! wind_speed is U_{forward} in Eq.(10) 
    ros_b(:) = ros_f(:) * EXP(-0.012 * wind_speed(:)) 
    WHERE (ros_f(:) .LT. 0.05) 
      ros_b(:) = zero 
    ENDWHERE 

    ! Calculate fire duration as a function of fdi; Eq.(14); unit: min
    WHERE (fdi(:) .GT. min_stomate)
      fire_durat(:)=241.0/(1.0+(240.0*EXP(-11.06*fdi(:)))) 
    ELSEWHERE
      fire_durat(:)=zero
    ENDWHERE

    db(:)=ros_b(:)*fire_durat(:) !unit: m/min * min = m
    df(:)=ros_f(:)*fire_durat(:)
    
    !! Calculate area burnt. 
    mean_fire_size(:)=(pi/(4.0*lb(:)))*((df(:)+db(:))**2.0)/10000.0 !mean_fire_size with unit as ha
    
    area_burnt_orig(:)=numfire(:)* mean_fire_size(:) 

    ! This baadjust_ratio was used to calibrate burned area with GFED4 burned area data
    area_burnt(:)=area_burnt_orig(:)*baadjust_ratio(:) 
    
    ! Ensure burned area within the limit of burnable vegetated area
    WHERE (area_burnt(:) .GT. area_burnable_veg(:))
      area_burnt(:)=area_burnable_veg(:)
    ENDWHERE

    ! AED_FEEDBACK (Marie 2026): accumulate burnt area for edge_length feedback.
    ! area_burnt is in ha; convert to m² before adding to the per-step accumulator.
    IF (ok_aed_feedback .AND. ALLOCATED(dA_fire)) THEN
      dA_fire(:) = dA_fire(:) + area_burnt(:) * 10000.0_r_std
    ENDIF

    !WRITE(numout,*) 'CheckoutB', area_burnt(itestipts),baadjust_ratio(itestipts),&
   !mean_fire_size(itestipts),db(itestipts),df(itestipts),fire_durat(itestipts),&
   !wind_speed(itestipts),numfire(itestipts),ignition_efficiency(itestipts)
    ! Adjust numfire to conserve the equation of 'numfire*mean_fire_size = burnt_area',
    ! after adjusting the burned area with adjusting ratio.
    WHERE (area_burnt_orig(:) .GT. min_stomate)
      numfire(:) = numfire(:) * area_burnt(:) / area_burnt_orig(:)
      numfire_lightning(:) = numfire_lightning(:) * area_burnt(:) / area_burnt_orig(:)
      numfire_human(:) = numfire_human(:) * area_burnt(:) / area_burnt_orig(:)
    ENDWHERE
    
    
    WHERE(area_burnable_veg(:) .GT. min_stomate)
      fire_frac(:) = area_burnt(:)/area_burnable_veg(:)
    ELSEWHERE
      fire_frac(:) = zero
    ENDWHERE

    !! 2. Calculate daily fire emissions

    ! Calculate livegrass fuel moisture
    WHERE (livegrass_grid(:) .GT. min_stomate)
      ! Eq(B2); soilhum_daily is $\omiga_{s,l}$ in Eq.(B2)
      dfm_livegrass(:) = max(zero,((10.0/9.0)*soilhum_daily(:)-(1.0/9.0)))  
    ELSEWHERE
      dfm_livegrass(:)=1.0
    ENDWHERE
    CALL combustion_fraction(npts, livegrass_grid, litterfuel_class_grid,    &
                 dfm_1hr, dfm_livegrass, wetness, me_grid, cf_litterfuel,    &
                 cf_lg, AED_fire)

    ! Calculate litter fuel consumed based on burned area only and per PFT area basis,
    ! respectively.
    litterfuel_fc_ba(:,:,:,:,:) = zero
    DO ipts = 1,npts
      DO ivm = 1,nvm 
        IF (burnable(ivm)) THEN
          DO ihour = 1,nhour
            litterfuel_fc_ba(ipts,:,ivm,ihour,:) = litterfuel(ipts,:,ivm,ihour,:) &
                     *cf_litterfuel(ipts,ihour)
            litterfuel_consumed(ipts,:,ivm,ihour,:) = litterfuel(ipts,:,ivm,ihour,:) &
                     *cf_litterfuel(ipts,ihour)*fire_frac(ipts)
          ENDDO
        ENDIF
      ENDDO
    ENDDO

    ! Calculate surface fire frontline intensity.
    ! We calculate fine fuel consumed during fire, which include fuels of 1h, 10h and 100h.
    ! Note the fine fuel consumption should integrate all burnable PFTs by using their 
    ! relative fractions of burned area as the weighting factor.
    ! Only fine fuel consumed contributes to the calculation of fire frontline intensity.
    ! The first SUM over DIM=2 sums over the dimension of litter type, and the second SUM
    ! over DIM=3 sums over the dimension of fuel type, thus yielding the two remaining 
    ! dimensions as (npts,nvm). The third SUM sums over the dimension of nvm, thus leaving 
    ! the last dimension as npts. Note we have to convert first to dry mass and then subtract
    ! the mineral content. Note `finefuel_fc_ba` means fine fuel consumption on burned area basis.
    ! We have to integrate all burnable PFTs, whose fractions in total burned area are assumed
    ! in proportion to their `relative_cover`.
    finefuel_fc_ba(:) = zero
    finefuel_fc_ba(:) = SUM( SUM(SUM(litterfuel_fc_ba(:,:,:,ihour1:ihour100,icarbon),DIM=2),DIM=3) * &
                                relative_cover(:,:), DIM=2) 

    ! We have to convert fine fuel consumption (based on burned area) from C to dry mass, excluding
    ! the mineral content. This new variable (g dry mass m^{-2}) is then used to calculate 
    ! fire frontline intensity, following Eq. (15).
    frontline_intensity(:) = zero
    vartmp(:) = finefuel_fc_ba(:) * carbon_to_drymass * (1.0-MINER_TOT)
    WHERE ( fire_frac(:) .GT. min_stomate .AND. vartmp(:) .GT. min_stomate )
      frontline_intensity(:)=H*(vartmp(:)/1000.0/fire_frac(:))*(ros_f(:)/60.0)  
    ENDWHERE

    cf_finefuel(:) = zero
    WHERE(litterfuel_fine(:) .GT. min_stomate)
      cf_finefuel(:) = finefuel_fc_ba(:)/litterfuel_fine(:)
    ENDWHERE

    DO ipts = 1,npts
      IF (cf_finefuel(ipts) .GT. 1) THEN
        WRITE(numout,*) 'Grid cell number: ', ipts
        CALL ipslerr_p (3,'spitfire',&
             'combustion fraction for fine fuel exceeds one','','')
      ENDIF
    ENDDO

    ! Apply consequences of litter fuel combustion on existing litter and fuel variables.
    !emissions_fire(:,:,:,ifiresurface) = SUM(SUM(litterfuel_consumed,DIM=2),DIM=3)/dt_days 
    !litterfuel(:,:,:,:,:) = litterfuel(:,:,:,:,:) - litterfuel_consumed(:,:,:,:,:)
    !litter(:,:,:,iabove,:) = litter(:,:,:,iabove,:) - SUM(litterfuel_consumed(:,:,:,:,:),DIM=4)
    
    !litterfuel(:,:,:,:,icarbon) = litterfuel(:,:,:,:,icarbon) - litterfuel_consumed(:,:,:,:,icarbon)
    !litter(:,:,:,iabove,icarbon) = litter(:,:,:,iabove,icarbon) - SUM(litterfuel_consumed(:,:,:,:,icarbon),DIM=4)
    !! 3. Calculate fire-driven forest mortality

    !! 3.1 Prepare stand structure information to derive tree mortality
    canopy_height(:,:,:) = zero
    crown_length(:,:,:) = zero
    DO ivm = 1,nvm
      IF(is_tree(ivm)) THEN
        DO ipts = 1,npts
          ! Calculating diameter at breast height (in meter) of the
          ! representative tree for all circ classes.
          diameter(ipts,ivm,:) = wood_to_dia(circ_class_biomass(ipts,ivm,:,:,icarbon),&
             ivm,pipe_tune2(ipts,ivm))

          ! Calculating the height in m for the mean tree of the given circumference class 
          canopy_height(ipts,ivm,:) = wood_to_height(circ_class_biomass(ipts,ivm,:,:,icarbon),&
             ivm,pipe_tune2(ipts,ivm))
          
          ! Vertical and horizontal crown diameter
          CALL wood_to_cn_dia(circ_class_biomass(ipts,ivm,:,:,icarbon), &
               circ_class_n(ipts,ivm,:),ivm,crown_width(ipts,ivm,:), &
               crown_length(ipts,ivm,:),pipe_tune2(ipts,ivm))
        ENDDO
      ENDIF
    ENDDO
    
    crown_base_height(:,:,:) = zero
    crown_base_height(:,:,:) = canopy_height(:,:,:) - crown_length(:,:,:)

    !! 3.2 Calculate tree mortality due to fire
    ck(:,:,:)=zero
    pm_ck(:,:,:)=zero
     
    DO ipts = 1,npts

      ! Fire can only cause damage when surface fire intensity is above a threshold
      ! Calculate scorching height: the fact that the relationship between scorching height
      ! and fire intensity depends on PFT is counter-intuitive, but we will keep using this
      ! empirical relationship in Thonicke et al. (2010)
      IF (frontline_intensity(ipts) .GE. surface_threshold) THEN
        scorching_height(ipts,:) = f_sh(:)*(frontline_intensity(ipts)**0.667)  ! applying Eq.(16); 
                                                 ! f_sh(ivm) is defined in constantes_mtc.f90
      ELSE
        scorching_height(ipts,:) = zero       
        tau_l(ipts,:) = zero
      ENDIF
    ENDDO

    ! Here we calculate the maximum fire size constrained by fragmentation
    ! We need tree information for accounting the fragmentation effects on fire size
    ! In the MICT version, we don't have circumfirence classes and thus only one crown_base_height. 
    ! If scorching_height is higher than crown_base_height, crown fire happens and no limitation on fire size.
    ! In the TRUNK, if we have many circumfirence classes (e.g. 20), This introduces small trees 
    ! with low crown base heights, which can be easily exceeded by the scorching height.
    ! To address this, we use the fraction of ignited canopy biomass as a proxy for crown fire initiation.
    ! When this fraction exceeds a prescribed threshold (cf_threshold), crown fire is assumed to occur.
    
    IF (ok_aed_size) THEN
      IF (printlev_loc .GE. 5) THEN
        WRITE(numout,*) 'Before ok_aed_size'
        WRITE(numout,*) 'numfire_human',numfire_human(:)
        WRITE(numout,*) 'frontline_intensity',frontline_intensity(:)
        WRITE(numout,*) 'area_burnt', area_burnt(:)
        WRITE(numout,*) 'area_burnt_orig',area_burnt_orig(:)
        WRITE(numout,*) 'litterfuel_consumed',litterfuel_consumed(:,:,:,:,:)
        WRITE(numout,*) 'litterfuel',litterfuel(:,:,:,:,:)
        WRITE(numout,*) 'SUM(litterfuel(1,:,2,:,1),DIM=2)',SUM(litterfuel(1,:,2,:,1),DIM=2)
        WRITE(numout,*) 'cf_litterfuel',cf_litterfuel(:,:)
        mean_fire_size_temp(:) = mean_fire_size(:)
      ENDIF
      !3.2.1. Recalculate mean fire size
      ignited_weight(:,:) = 0.0
      total_weight(:,:)   = 0.0
      
      DO ipts = 1,npts
        DO ivm = 1,nvm
          IF(is_tree(ivm)) THEN
            DO icir = 1,ncirc
                total_weight(ipts, ivm) = total_weight(ipts, ivm) + &
                    SUM(circ_class_biomass(ipts,ivm,:,icir,icarbon),1)
                IF (scorching_height(ipts,ivm) .GE. crown_base_height(ipts, ivm, icir)) THEN
                  ignited_weight(ipts, ivm) = ignited_weight(ipts, ivm) + &
                    SUM(circ_class_biomass(ipts, ivm,:,icir,icarbon),1)
                ENDIF
            ENDDO
          ENDIF !IF(is_tree(ivm)) THEN
        ENDDO
      ENDDO !DO icir = 1,ncirc

      DO ipts = 1,npts
        DO ivm = 1,nvm
          IF (total_weight(ipts,ivm) .GT. min_stomate)  THEN
            ignited_fraction(ipts,ivm) = ignited_weight(ipts,ivm)/total_weight(ipts,ivm)
          ELSE
            ignited_fraction(ipts,ivm) = 0.0
          ENDIF
        ENDDO
      ENDDO

      WHERE (ignited_fraction(:,:) .GE. cf_threshold)
        crown_fire(:,:) = .TRUE.
      ELSEWHERE
        crown_fire(:,:) = .FALSE.
      ENDWHERE
      
      max_patch_fire_size_PFT(:,:)=0.0
      max_patch_fire_size(:)=0.0
      grass_biomass_wet(:,:)=0.0
    
      DO ivm = 1,nvm
        IF( burnable(ivm) .AND. (.NOT. is_tree(ivm)) ) THEN
          DO ihour=1,nhour
            grass_biomass_wet(:,ivm)=grass_biomass_wet(:,ivm)+litterfuel_class_pft(:,ihour,ivm)*(1/0.45)*0.01
          ENDDO
          grass_biomass_wet(:,ivm)=grass_biomass_wet(:,ivm)+livegrass_pft(:,ivm)*(1/0.45)*0.01
          !(1/0.45 converts back to wet ; 0.01 =g/m2 --> ton/ha)
        ENDIF
      ENDDO

      Patch_Area(:)=pi*(AED_fire(:)**2)/10000.0 ! m2 --> ha
      DO ivm=1,nvm
        IF (burnable(ivm)) THEN
          IF (is_tree(ivm)) THEN
            WHERE(crown_fire(:,ivm))
              max_patch_fire_size_PFT(:,ivm) = area_land(:)*veget_max(:,ivm)
            ELSEWHERE
              max_patch_fire_size_PFT(:,ivm) = Patch_Area(:)*patch_multiplier_forest
            ENDWHERE
          ELSE !is grass
            WHERE (grass_biomass_wet(:,ivm) .GE. 2.4) ! 
              max_patch_fire_size_PFT(:,ivm) = area_land(:)*veget_max(:,ivm)
            ELSEWHERE
              max_patch_fire_size_PFT(:,ivm) = Patch_Area(:)*patch_multiplier_grass
            ENDWHERE
          ENDIF !IF (is_tree(ivm)) THEN
        ENDIF !IF (burnable(:)) THEN
      ENDDO !DO ivm=1,nvm

      DO ivm = 1,nvm
        IF( burnable(ivm) ) THEN
          max_patch_fire_size(:) = max_patch_fire_size(:) + max_patch_fire_size_PFT(:,ivm)*veget_max(:,ivm)
        ENDIF
      ENDDO

      DO ipts = 1,npts
        IF (burnable_vegfrac(ipts) .GT. min_stomate)  THEN
          max_patch_fire_size(ipts) = max_patch_fire_size(ipts)/burnable_vegfrac(ipts)
        ELSE
          max_patch_fire_size(ipts) = area_land(ipts)
          CALL ipslerr_p(2,'spitfire AED', '', '',&
          &                 'burnable_vegfrac too small') ! Warning error
        ENDIF
      ENDDO

      IF (printlev_loc .GE. 5) THEN
        WRITE(numout,*) 'grass_biomass_wet',grass_biomass_wet(:,:)
        WRITE(numout,*) 'crown_fire',crown_fire(:,:)
        WRITE(numout,*) 'area_land',area_land(:)
        WRITE(numout,*) 'mean_fire_size without fragmentation constrain',mean_fire_size(:)
        WRITE(numout,*) 'max_patch_fire_size',max_patch_fire_size(:)
        IF (ANY(mean_fire_size(:) .GT. max_patch_fire_size(:))) THEN
          WRITE(numout,*) 'mean_fire_size_change'
        ENDIF
      ENDIF
      mean_fire_size(:) = MIN(mean_fire_size(:),max_patch_fire_size(:))

      
      ! We recalculate everythin using the updated mean_fire_size
      area_burnt_orig(:)=numfire(:)* mean_fire_size(:) 

      ! This baadjust_ratio was used to calibrate burned area with GFED4 burned area data
      area_burnt(:)=area_burnt_orig(:)*baadjust_ratio(:) 
      
      ! Ensure burned area within the limit of burnable vegetated area
      WHERE (area_burnt(:) .GT. area_burnable_veg(:)) 
        area_burnt(:)=area_burnable_veg(:)
      ENDWHERE
  
      !WRITE(numout,*) 'CheckoutB', area_burnt(itestipts),baadjust_ratio(itestipts),&
     !mean_fire_size(itestipts),db(itestipts),df(itestipts),fire_durat(itestipts),&
     !wind_speed(itestipts),numfire(itestipts),ignition_efficiency(itestipts)
      ! Adjust numfire to conserve the equation of 'numfire*mean_fire_size = burnt_area',
      ! after adjusting the burned area with adjusting ratio.
      WHERE (area_burnt_orig(:) .GT. min_stomate)
        numfire(:) = numfire(:) * area_burnt(:) / area_burnt_orig(:)
        numfire_lightning(:) = numfire_lightning(:) * area_burnt(:) / area_burnt_orig(:)
        numfire_human(:) = numfire_human(:) * area_burnt(:) / area_burnt_orig(:)
      ENDWHERE
      
      
      WHERE(area_burnable_veg(:) .GT. min_stomate)
        fire_frac(:) = area_burnt(:)/area_burnable_veg(:)
      ELSEWHERE
        fire_frac(:) = zero
      ENDWHERE
  
      !! 3.2.2. Recalculate daily fire emissions
  
      ! Calculate livegrass fuel moisture
      !WHERE (livegrass_grid(:) .GT. min_stomate)
        ! Eq(B2); soilhum_daily is $\omiga_{s,l}$ in Eq.(B2)
      !  dfm_livegrass(:) = max(zero,((10.0/9.0)*soilhum_daily(:)-(1.0/9.0)))  
      !ELSEWHERE
      !  dfm_livegrass(:)=1.0
      !ENDWHERE
      CALL combustion_fraction(npts, livegrass_grid, litterfuel_class_grid,    &
                   dfm_1hr, dfm_livegrass, wetness, me_grid, cf_litterfuel,    &
                   cf_lg, AED_fire)

      ! Calculate litter fuel consumed based on burned area only and per PFT area basis,
      ! respectively.
      litterfuel_fc_ba(:,:,:,:,:) = zero
      litterfuel_consumed(:,:,:,:,:) = zero
      DO ipts = 1,npts
        DO ivm = 1,nvm 
          IF (burnable(ivm)) THEN
            DO ihour = 1,nhour
              litterfuel_fc_ba(ipts,:,ivm,ihour,:) = litterfuel(ipts,:,ivm,ihour,:) &
                       *cf_litterfuel(ipts,ihour)
              litterfuel_consumed(ipts,:,ivm,ihour,:) = litterfuel(ipts,:,ivm,ihour,:) &
                       *cf_litterfuel(ipts,ihour)*fire_frac(ipts)
            ENDDO
          ENDIF
        ENDDO
      ENDDO
  
      ! Calculate surface fire frontline intensity.
      ! We calculate fine fuel consumed during fire, which include fuels of 1h, 10h and 100h.
      ! Note the fine fuel consumption should integrate all burnable PFTs by using their 
      ! relative fractions of burned area as the weighting factor.
      ! Only fine fuel consumed contributes to the calculation of fire frontline intensity.
      ! The first SUM over DIM=2 sums over the dimension of litter type, and the second SUM
      ! over DIM=3 sums over the dimension of fuel type, thus yielding the two remaining 
      ! dimensions as (npts,nvm). The third SUM sums over the dimension of nvm, thus leaving 
      ! the last dimension as npts. Note we have to convert first to dry mass and then subtract
      ! the mineral content. Note `finefuel_fc_ba` means fine fuel consumption on burned area basis.
      ! We have to integrate all burnable PFTs, whose fractions in total burned area are assumed
      ! in proportion to their `relative_cover`.
      finefuel_fc_ba(:) = zero
      finefuel_fc_ba(:) = SUM( SUM(SUM(litterfuel_fc_ba(:,:,:,ihour1:ihour100,icarbon),DIM=2),DIM=3) * &
                                  relative_cover(:,:), DIM=2) 
  
      ! We have to convert fine fuel consumption (based on burned area) from C to dry mass, excluding
      ! the mineral content. This new variable (g dry mass m^{-2}) is then used to calculate 
      ! fire frontline intensity, following Eq. (15).
      frontline_intensity(:) = zero
      vartmp(:) = finefuel_fc_ba(:) * carbon_to_drymass * (1.0-MINER_TOT)
      WHERE ( fire_frac(:) .GT. min_stomate .AND. vartmp(:) .GT. min_stomate )
        frontline_intensity(:)=H*(vartmp(:)/1000.0/fire_frac(:))*(ros_f(:)/60.0)  
      ENDWHERE
  
      cf_finefuel(:) = zero
      WHERE(litterfuel_fine(:) .GT. min_stomate)
        cf_finefuel(:) = finefuel_fc_ba(:)/litterfuel_fine(:)
      ENDWHERE
  
      DO ipts = 1,npts
        IF (cf_finefuel(ipts) .GT. 1) THEN
          WRITE(numout,*) 'Grid cell number: ', ipts
          CALL ipslerr_p (3,'spitfire',&
               'combustion fraction for fine fuel exceeds one','','')
        ENDIF
      ENDDO
       
      DO ipts = 1,npts
  
        ! Fire can only cause damage when surface fire intensity is above a threshold
        ! Calculate scorching height: the fact that the relationship between scorching height
        ! and fire intensity depends on PFT is counter-intuitive, but we will keep using this
        ! empirical relationship in Thonicke et al. (2010)
        IF (frontline_intensity(ipts) .GE. surface_threshold) THEN
          scorching_height(ipts,:) = f_sh(:)*(frontline_intensity(ipts)**0.667)  ! applying Eq.(16); 
                                                   ! f_sh(ivm) is defined in constantes_mtc.f90
        ELSE
          scorching_height(ipts,:) = zero       
          tau_l(ipts,:) = zero
        ENDIF
      ENDDO
      IF (printlev_loc .GE. 5) THEN
        WRITE(numout,*) 'After ok_aed_size'
        WRITE(numout,*) 'numfire_human',numfire_human(:)
        WRITE(numout,*) 'frontline_intensity',frontline_intensity(:)
        WRITE(numout,*) 'litterfuel_consumed',litterfuel_consumed(:,:,:,:,:)
        WRITE(numout,*) 'litterfuel',litterfuel(:,:,:,:,:)
        WRITE(numout,*) 'SUM(litterfuel(1,:,2,:,1),DIM=2)',SUM(litterfuel(1,:,2,:,1),DIM=2)
        WRITE(numout,*) 'cf_litterfuel',cf_litterfuel(:,:)
        IF (ANY(ABS(mean_fire_size_temp(:) - mean_fire_size(:)) .GT. 1e-20)) THEN
          WRITE(numout,*) 'Difference detected'
          WRITE(numout,*) 'mean_fire_size', mean_fire_size(:)
          WRITE(numout,*) 'mean_fire_size_temp', mean_fire_size_temp(:)
          WRITE(numout,*) 'area_burnt', area_burnt(:)
          WRITE(numout,*) 'area_burnt_orig', area_burnt_orig(:)
        ENDIF
      ENDIF
    ENDIF !IF (ok_aed_size) THEN
    ! Apply consequences of litter fuel combustion on existing litter and fuel variables.
    ! This should normally be done in Part 2. However, when considering fragmentation effects,
    ! the behavior may change, so this section is placed after ok_aed_size.
    emissions_fire(:,:,:,ifiresurface) = SUM(SUM(litterfuel_consumed,DIM=2),DIM=3)/dt_days 
    litterfuel(:,:,:,:,:) = litterfuel(:,:,:,:,:) - litterfuel_consumed(:,:,:,:,:)
    litter(:,:,:,iabove,:) = litter(:,:,:,iabove,:) - SUM(litterfuel_consumed(:,:,:,:,:),DIM=4)

    DO ipts = 1,npts
      ! calculate fire residence time
      IF (var_gamma(ipts) .GT. min_stomate) THEN
        tau_l(ipts,:) = 2.0*(cf_finefuel(ipts)/var_gamma(ipts)) ! Equation(46) in SPITFIRE technical document.
      ELSE
        tau_l(ipts,:) = zero
      ENDIF

      DO ivm = 1,nvm
        IF (burnable(ivm) .AND. is_tree(ivm)) THEN
          DO icir = 1,ncirc

            ! calculate fraction of tree crown scorched by fire
            IF (scorching_height(ipts,ivm) .GE. canopy_height(ipts,ivm,icir)) THEN
              ck(ipts,ivm,icir) = 1.0 
            ELSE IF (scorching_height(ipts,ivm) .LE. crown_base_height(ipts,ivm,icir)) THEN
              ck(ipts,ivm,icir) = zero
            ELSE
              ck(ipts,ivm,icir) = (scorching_height(ipts,ivm) - crown_base_height(ipts,ivm,icir)) &
                           /crown_length(ipts,ivm,icir)
            ENDIF

            ! post-fire mortality from crown scorching
            IF (ck(ipts,ivm,icir) .GT. min_stomate) THEN
              pm_ck(ipts,ivm,icir) = r_ck(ivm)*(ck(ipts,ivm,icir)**p_ck(ivm)) ! applying Eq.(22); 
                                               ! r_ck & p_ck are defined in constantes_mtc.f90
            ENDIF

            ! post-fire mortality from cambial damage
            ! calculate bark thickness and critical residence time
            bark_thickess(ipts,ivm,icir) = BTpar1(ivm) * diameter(ipts,ivm,icir) *100 + BTpar2(ivm)   ! applying Eq. (21)
                                ! Note that Eq. (21) requires the unit of DBH as centimeter.
            tau_c(ipts,ivm,icir) = 2.9 * (bark_thickess(ipts,ivm,icir)**2.0)   ! applying Eq. (20)

            ! Calculate mortality due to cambial damage
            ! tau_l/tau_c is the ratio of fire residence time against the critical residence time 
            ! for cambial damage.
            IF (tau_c(ipts,ivm,icir) .GE. min_stomate) THEN
              IF ( (tau_l(ipts,ivm)/tau_c(ipts,ivm,icir)) .GE. 2.0) THEN
                pm_tau(ipts,ivm,icir) = 1.0
              ELSE IF ((tau_l(ipts,ivm)/tau_c(ipts,ivm,icir)) .GT. 0.22) THEN
                ! original value for the intercept in Thonicke is 0.125, which will yield negative value for pm_tau.
                ! This causes mass balance error in the model. I hence change to 0.12386 which is 0.563*0.22
                pm_tau(ipts,ivm,icir) = (0.563*tau_l(ipts,ivm)/tau_c(ipts,ivm,icir))-0.12386
              ELSE IF ((tau_l(ipts,ivm)/tau_c(ipts,ivm,icir)) .LE. 0.22) THEN
                pm_tau(ipts,ivm,icir) = zero
              ENDIF
            ELSE
              pm_tau(ipts,ivm,icir) = zero    
            ENDIF
             
            ! Calculate post-fire mortality from crown scorching and cambial kill; Eq. (18)
            postf_mort(ipts,ivm,icir) = pm_tau(ipts,ivm,icir)+pm_ck(ipts,ivm,icir)  &
                                         -(pm_tau(ipts,ivm,icir)*pm_ck(ipts,ivm,icir))

            circ_class_kill(ipts,ivm,icir,forest_managed(ipts,ivm),icut_fire) = &
               circ_class_n(ipts,ivm,icir)*postf_mort(ipts,ivm,icir)*fire_frac(ipts)

          ENDDO ! icir = 1,ncirc
        ENDIF 
      ENDDO ! ivm = 1,nvm
    ENDDO ! ipts = 1,npts

    ! Prepare variables for final output



    !!++TODO: part of dead trees are crown-scorched and generate crown emissions. 
    !! this needs to be subtracted out from mortality.

    !!++TODO: combustion and mortality of live grass needs to be treated.

    !       ! Here explains how consequence of fire-killed trees are handled.
    !       ! the fraction of trees that are kill in fire is `postf_mort*fire_frac`,
    !       ! Of the killed trees, `ck` of them have been crown-scroached by fire
    !       ! and the 1-,10- and 100-h live fuel and 5% of the 1000-h fuel are combusted. 
    !       ! the biomass of trees that are crown-scoarched but not consumed in fire 
    !       ! would go to bm_to_litter. Then, the killed but not crown-scroached tree
    !       ! biomass (1-ck) go directly to bm_to_litter.
    !       
    !       ! In summary, all the biomass of killed tress should be removed from 
    !       ! live biomass pool as they are either combusted or transferred to bm_to_litter. 
    !       ! Of the fraction (postf_mort * fire_frac) that are killed in fire, of which
    !       !    --- ck (fraction) have been crown scorached, of which
    !       !        -- 1-,10- and 100-h live fule and 5% of the 1000-h fuel are combusted;
    !       !        -- 95% 1000-h fuel go to bm_to_litter;  
    !       !    --- (1-ck) (fraction) go to bm_to_litter.

    !       tree_kill_frac(:)=postf_mort(:,ivm)*fire_frac(:)
    !       !disturb_crown is carbon emssions from crown scorching.
    !       disturb_crown(:,ivm) = tree_kill_frac(:) * ck(:,ivm) *                    & 
    !        (biomass(:,ivm,ileaf) + biomass(:,ivm,ifruit) +  &
    !         0.0450 * biomass(:,ivm,isapabove) + 0.0450 * biomass(:,ivm,iheartabove) +     &
    !         0.0450 * biomass(:,ivm,icarbres) +                                         &
    !         0.0750 * biomass(:,ivm,isapabove) + 0.0750 * biomass(:,ivm,iheartabove) +     &
    !         0.0750 * biomass(:,ivm,icarbres) +                                         &
    !         0.2100 * biomass(:,ivm,isapabove) + 0.2100 * biomass(:,ivm,icarbres) +      &
    !         0.2100 * biomass(:,ivm,iheartabove) +                                       & 
    !         0.6700 * 0.050 * biomass(:,ivm,icarbres) + 0.6700 * 0.050 * biomass(:,ivm,iheartabove) + &
    !         0.6700 * 0.050 * biomass(:,ivm,isapabove)                               & 
    !        ) 

    !       !fire flux to atmosphere, now biomass emissions from trees are included.
    !       dcflux_fire_pft(:,ivm)=dcflux_fire_pft(:,ivm)+disturb_crown(:,ivm)
    !       fc_crown(:)=fc_crown(:)+disturb_crown(:,ivm)*veget_max(:,ivm)  

    !       ! biomass to litter transfer; (the killed but not crown scorching) AND (crown scorching but not combusted)
    !       ! 1-ck(:,ivm)---killed but not crown scorching;
    !       ! ck(:,ivm)*(1-0.045-0.075-0.21-0.67*0.05)---crown scorching but not combusted. 
    !       bm_to_litter(:,ivm,ileaf) = bm_to_litter(:,ivm,ileaf) +                           & 
    !                                 (tree_kill_frac(:) * (1-ck(:,ivm)) * biomass(:,ivm,ileaf)) 
    !       bm_to_litter(:,ivm,ifruit) = bm_to_litter(:,ivm,ifruit) +                         & 
    !                                 (tree_kill_frac(:) * (1-ck(:,ivm)) * biomass(:,ivm,ifruit)) 
    !       bm_to_litter(:,ivm,isapabove) = bm_to_litter(:,ivm,isapabove) +                   & 
    !                        (tree_kill_frac(:) * (1-ck(:,ivm)+ck(:,ivm)*  &
    !                        (1-0.045-0.075-0.21-0.67*0.05))* biomass(:,ivm,isapabove)) 
    !       bm_to_litter(:,ivm,icarbres) = bm_to_litter(:,ivm,icarbres) +                     &
    !                        (tree_kill_frac(:) * (1-ck(:,ivm)+ck(:,ivm)*  &
    !                        (1-0.045-0.075-0.21-0.67*0.05))* biomass(:,ivm,icarbres)) 
    !       bm_to_litter(:,ivm,iheartabove) = bm_to_litter(:,ivm,iheartabove) +               &
    !                        (tree_kill_frac(:) * (1-ck(:,ivm)+ck(:,ivm)*  &
    !                        (1-0.045-0.075-0.21-0.67*0.05))* biomass(:,ivm,iheartabove)) 
    !       ! move root/sapbelow/heartbelow biomass to belowground litter pool.
    !       bm_to_litter(:,ivm,iroot) = bm_to_litter(:,ivm,iroot) +                           & 
    !                                 (tree_kill_frac(:) * biomass(:,ivm,iroot)) 
    !       bm_to_litter(:,ivm,isapbelow) = bm_to_litter(:,ivm,isapbelow) +           & 
    !                                 (tree_kill_frac(:) * biomass(:,ivm,isapbelow)) 
    !       bm_to_litter(:,ivm,iheartbelow) = bm_to_litter(:,ivm,iheartbelow) +           & 
    !                                 (tree_kill_frac(:) * biomass(:,ivm,iheartbelow)) 

    !       ! biomass update after combustion and litter transfer.
    !       biomass(:,ivm,1) = biomass(:,ivm,1) - (tree_kill_frac(:) * biomass(:,ivm,1)) 
    !       biomass(:,ivm,2) = biomass(:,ivm,2) - (tree_kill_frac(:) * biomass(:,ivm,2)) 
    !       biomass(:,ivm,3) = biomass(:,ivm,3) - (tree_kill_frac(:) * biomass(:,ivm,3)) 
    !       biomass(:,ivm,4) = biomass(:,ivm,4) - (tree_kill_frac(:) * biomass(:,ivm,4)) 
    !       biomass(:,ivm,5) = biomass(:,ivm,5) - (tree_kill_frac(:) * biomass(:,ivm,5)) 
    !       biomass(:,ivm,6) = biomass(:,ivm,6) - (tree_kill_frac(:) * biomass(:,ivm,6)) 
    !       biomass(:,ivm,7) = biomass(:,ivm,7) - (tree_kill_frac(:) * biomass(:,ivm,7)) 
    !       biomass(:,ivm,8) = biomass(:,ivm,8) - (tree_kill_frac(:) * biomass(:,ivm,8)) 


    !       ! Update number of individuals surviving fire, ready for the next day. 
    !       ind(:,ivm)= ind(:,ivm) - nind_kill(:,ivm)

    !     ELSE !grass

    !       WHERE (frontline_intensity(:).ge.surface_threshold .AND. PFTpresent(:,ivm)) 
    !         !grass leaf & fruit consumption are also included in litter consumption
    !         litter_consump_pft(:,ivm)=litter_consump_pft(:,ivm)+fc_lg(:,ivm)+fc_lf(:,ivm) !leaf & fruit consumption in fire for GRASS also included in litter consumption
    !         biomass(:,ivm,ileaf)=biomass(:,ivm,ileaf)-fc_lg(:,ivm)
    !         biomass(:,ivm,ifruit)=biomass(:,ivm,ifruit)-fc_lf(:,ivm)
    !       ENDWHERE
    !     ENDIF   !IF (tree(ivm))

    !     ! add ground fuel combustion to the output variable of fire carbon flux per PFT
    !     ! dcflux_fire_pft:total fire emissions, including crown fuels and ground litter consumption
    !     dcflux_fire_pft(:,ivm)=dcflux_fire_pft(:,ivm) + litter_consump_pft(:,ivm)
    !     litter_consump(:)=litter_consump(:)+(litter_consump_pft(:,ivm) *veget_max(:,ivm)/dt_days)
    !   ENDIF   !IF burnable
    ! ENDDO  !fire effects per PFTj
  
   
    IF (printlev_loc .GE. 4) THEN
      WRITE(numout,*) 'Debuging for PFT: ',itestpft, 'pixel: ',itestipts
      DO iele=1,nelements
        WRITE(numout,*) 'Element: ', iele
        ! output only the compartment biomass by summing over the axis of circ class
        WRITE(numout,*) 'Biomass: ', SUM(circ_class_biomass(itestipts,itestpft,:,:,iele),DIM=1)
        write(numout,*) 'litter: ', litter(itestipts,:,itestpft,iabove,iele)
        ! output by summing ove the axis of fuel category
        write(numout,*) 'litterfuel: ', SUM(litterfuel(itestipts,:,itestpft,:,iele),DIM=2)
        write(numout,*) 'veget_max: ', veget_max(itestipts,itestpft)
      ENDDO
    ENDIF

    IF (printlev_loc .GE. 5) THEN
      WRITE(numout,*) 'Check vars: '
      WRITE(numout,*) 'ros_f: ',ros_f(:)
      WRITE(numout,*) 'wetness: ',wetness(:)
      WRITE(numout,*) 'me_grid: ',me_grid(:)
      WRITE(numout,*) 'dfm_livegrass: ',dfm_livegrass(:)
      WRITE(numout,*) 'dfm: ',dfm(:)
      WRITE(numout,*) 'dfm_1hr: ',dfm_1hr(:)
      WRITE(numout,*) 'ni_acc: ',ni_acc(:)
      DO ihour=1,nhour
        WRITE(numout,*) 'ihour,cf_litterfuel : ',ihour,cf_litterfuel(:,ihour)
      ENDDO
    ENDIF

    !! 3. Check numerical consistency of this routine
    IF (err_act .GT. 1) THEN

       ! 3.1. Check surface area
       ! No need to check surface area because nothing is done to change veget_max
       ! CALL check_vegetation_area("spitfire", npts, veget_max_begin, &
       !      veget_max,'pft')

       ! 3.2. Calculate final biomass
       pool_end(:,:,:) = zero
       DO iele = 1,nelements

          ! Biomass pool
          DO ipar = 1,nparts
             DO icir =1, ncirc
                pool_end(:,:,iele) = pool_end(:,:,iele) + &
                     (circ_class_biomass(:,:,icir,ipar,iele) * &
                     circ_class_n(:,:,icir)* veget_max(:,:))
             ENDDO

          ENDDO

          ! Litter pool (gC m-2) *  (m2 m-2) 
          DO ilit = 1,nlitt
             DO ilev = 1,nlevs
                pool_end(:,:,iele) = pool_end(:,:,iele) + &
                     litter(:,ilit,:,ilev,iele) * veget_max(:,:)
             ENDDO
          ENDDO

       ENDDO


       ! 3.3. Calculate components of the mass balance
       check_intern(:,:,iatm2land,:) = zero
       DO iele=1,nelements
          check_intern(:,:,iland2atm,iele) = un * emissions_fire(:,:,iele,ifiresurface) * veget_max(:,:)
       ENDDO
       check_intern(:,:,ilat2out,:) = zero
       check_intern(:,:,ilat2in,:) = -un * zero
       check_intern(:,:,ipoolchange,:) = &
            un * (pool_end(:,:,:) - pool_start(:,:,:))

       closure_intern(:,:,:) = zero
       DO imbc = 1,nmbcomp
          !DO iele=1,1
          DO iele=1,nelements
             ! Debug
             IF (printlev_loc>=4) THEN
                WRITE(numout,*) 'check_intern, ivm, imbc, iele, ', imbc, &
                     iele, check_intern(:,test_pft,imbc,iele)
             ENDIF
             !-
             closure_intern(:,:,iele) = closure_intern(:,:,iele) + &
                  check_intern(:,:,imbc,iele)
          ENDDO
       ENDDO

       !! 3.4. Check mass balance closure
       CALL check_mass_balance("spitfire", closure_intern, npts, &
            pool_end, pool_start, veget_max,'pft')

    ENDIF ! IF(err_act .GT. 1)


    !! 4. Prepare variables for output and make adjustments for time steps

    fire_frac_grid(:) = fire_frac(:) * burnable_vegfrac(:) / dt_days 
    area_burnt(:) = area_burnt(:)/dt_days
    numfire(:) = numfire(:)/dt_days
    numfire_lightning(:) = numfire_lightning(:)/dt_days
    numfire_human(:) = numfire_human(:)/dt_days

    ! Calculate trace gas emissions
    DO x=1,ntrace !trace_species
      DO ivm=1,nvm
        emissions_trace_gas(:,ivm,x)=emissions_fire(:,ivm,icarbon,ifiresurface) &
                                     * carbon_to_drymass/1000. * ef_trace(ivm,x)
      ENDDO
    ENDDO

    IF (ok_spitfire) THEN
        DO ihour=1,nhour
          IF (ihour == ihour1) THEN
                  fueltype_str(ihour) = '1hr'
          ELSEIF (ihour == ihour10) THEN
                  fueltype_str(ihour) = '10hr'
          ELSEIF (ihour == ihour100) THEN
                  fueltype_str(ihour) = '100hr'
          ELSE
                  fueltype_str(ihour) = '1000hr'
          ENDIF
          CALL xios_orchidee_send_field("LITTERFUEL_c_"//TRIM(fueltype_str(ihour)),SUM(litterfuel(:,:,:,ihour,icarbon),DIM=2))
        ENDDO

        !CALL xios_orchidee_send_field("CO2_FIRE",emissions_fire(:,:,icarbon,ifiresurface))
        CALL xios_orchidee_send_field("AREA_BURNABLE",area_burnable_veg)
        CALL xios_orchidee_send_field("AREA_BURNT",area_burnt)
        CALL xios_orchidee_send_field("FIRE_FRACTION",fire_frac_grid)

        DO ivm = 1,nvm
          vartmp_pft(:,ivm) = area_burnt(:) * relative_cover(:,ivm)
        ENDDO
        
        CALL xios_orchidee_send_field("AREA_BURNT_PFT",vartmp_pft)
        CALL xios_orchidee_send_field("FDI",fdi)
        CALL xios_orchidee_send_field("NUMFIRE",numfire)
        CALL xios_orchidee_send_field("NUMFIRE_LIGHTNING",numfire_lightning)
        CALL xios_orchidee_send_field("NUMFIRE_HUMAN",numfire_human)
        CALL xios_orchidee_send_field("POPDENS",popd)
        CALL xios_orchidee_send_field("A_ND",a_nd)
        CALL xios_orchidee_send_field("FIRESIZE",mean_fire_size)
        CALL xios_orchidee_send_field("IGNITION_EFFICIENCY",ignition_efficiency)
        CALL xios_orchidee_send_field("FDI_NIACC", fdi_niacc)
        CALL xios_orchidee_send_field("FDI_VPD_DAILY_MEAN", fdi_vpd_daily_mean)
        CALL xios_orchidee_send_field("FDI_VPD_DAILY_MAX",  fdi_vpd_daily_max)
        CALL xios_orchidee_send_field("FDI_VPD_MEAN_WEEK",  fdi_vpd_mean_week)
        CALL xios_orchidee_send_field("FDI_VPD_MAX_WEEK",   fdi_vpd_max_week)
        CALL xios_orchidee_send_field("CCFIRE_MOR_c",   SUM(circ_class_kill(:,:,:,:,icut_fire),DIM=4))

        !! Harmonised perturbed forest surface from fire kills (circ-class basal-area route,
        !! shared cc_kill_to_area). Same metric/unit (m2) as AREA_DAMAGED_{STORM,PEST,HARVEST},
        !! and distinct from AREA_BURNT (fire-spread footprint, ha/day). Daily per-event
        !! increment ; XIOS operation="accumulate" yields a GROSS area over the output interval
        !! (uniform-pixel, double-counts repeated disturbance, may exceed stand area). NOT a
        !! bounded cumulative fraction -- do not sum across years (see
        !! DEV_HARMONIZE_DISTURBANCE_AREA.md, Jina's note). The forest-management dimension is
        !! collapsed by SUM(.,DIM=3): only the forest_managed slot is non-zero (cf. CCFIRE_MOR_c).
        area_damaged_fire(:,:) = zero
        DO ivm = 1,nvm
          IF ( is_tree(ivm) ) THEN
            CALL cc_kill_to_area(npts, ivm, SUM(circ_class_kill(:,ivm,:,:,icut_fire), DIM=3), &
                 circ_class_biomass(:,ivm,:,:,:), circ_class_n(:,ivm,:), veget_max(:,ivm),    &
                 area(:)*contfrac(:), pipe_tune2(:,ivm), area_damaged_fire(:,ivm))
          ENDIF
        ENDDO
        CALL xios_orchidee_send_field("AREA_DAMAGED_FIRE", area_damaged_fire)

        IF (ok_aed_humign .OR. ok_aed_fuel .OR. ok_aed_size .OR. ok_aed_wind) THEN
          CALL xios_orchidee_send_field("AED",AED_fire)
        ENDIF
        !CALL xios_orchidee_send_field("AREA_BURNT_HUMAN",numfire_human(:)*mean_fire_size(:))
        !CALL xios_orchidee_send_field("AREA_BURNT_LIGHTNING",numfire_lightning(:)*mean_fire_size(:))
        !CALL xios_orchidee_send_field("LIGHTN_IGN",lightn_ign)
        !CALL xios_orchidee_send_field("EMI_FIRE_CO",emissions_trace_gas(:,:,2))
        !CALL xios_orchidee_send_field("EMI_FIRE_CH4",emissions_trace_gas(:,:,3))
        !CALL xios_orchidee_send_field("EMI_FIRE_VOC",emissions_trace_gas(:,:,4))
        !CALL xios_orchidee_send_field("EMI_FIRE_TPM",emissions_trace_gas(:,:,5))
        !CALL xios_orchidee_send_field("EMI_FIRE_NOx",emissions_trace_gas(:,:,6))
    ENDIF
    END SUBROUTINE spitfire


!! ================================================================================================================================
!! SUBROUTINE     : fire_danger_index
!!                                                                              
!>\BRIEF        Calculation of fire danger index using the Nesterov Index      
!!                                                                              
!! DESCRIPTION   : We use an empirical relationship linking fuels of different surface-to-volume ratios
!!                 and weather conditions to fuel moisture, which is further compared with the moisture
!!                 of extinction to get the probability of successful ignitions. 
!!                                                                              
!! REFERENCE(S)   : 
!! - K.Thonicke, A.Spessa, I.C. Prentice, S.P.Harrison, L.Dong, and
!! C.Carmona-Moreno, 2010: The influence of vegetation, fire spread and fire
!! behaviour on biomass burning and trace gas emissions: results from a
!! process-based model. Biogeosciences, 7, 1991-2011.
!!                                                                              
!! \n                                                                           
!_ ================================================================================================================================

    SUBROUTINE  fire_danger_index(npts, nvm, veget, me_grid, litterfuel_class_grid, ni_acc, &
                                  vpd_daily_mean, vpd_daily_max, vpd_mean_week, vpd_max_week, precip_month, vd, & 
                                  dfm, dfm_1hr, fdi_niacc, fdi_vpd_daily_mean, fdi_vpd_daily_max, fdi_vpd_mean_week, &
                                  fdi_vpd_max_week)


    IMPLICIT NONE

    ! Input variables
    INTEGER, INTENT(in)                                  :: npts                    !! Domain size
    INTEGER, INTENT(in)                                  :: nvm
    REAL(r_std),DIMENSION(npts,nvm), INTENT(in)          :: veget                   !! Fractional coverage: actually share of the pixel covered 
                                                                                    !! by a PFT taking into account LAI
    REAL(r_std), DIMENSION(npts), INTENT(in)             :: me_grid                 !! Moisture of extinction for the whole grid (unitless)
    REAL(r_std), DIMENSION(npts,nhour), INTENT(in)       :: litterfuel_class_grid   !! Gridcell aboveground litter fuel for different fuel classes
                                                                                    !! @tex $(gC m^{-2})$ @endtex per ground area
    REAL(r_std), DIMENSION(npts),INTENT(in)              :: ni_acc                  !! Nesterov index (square of degree Celcius)  
    REAL(r_std), DIMENSION(npts),INTENT(in)              :: vpd_daily_mean          !! Daily VPD
    REAL(r_std), DIMENSION(npts),INTENT(in)              :: vpd_daily_max
    REAL(r_std), DIMENSION(npts),INTENT(in)              :: vpd_mean_week
    REAL(r_std), DIMENSION(npts),INTENT(in)              :: vpd_max_week
    REAL(r_std), DIMENSION(npts),INTENT(in)              :: precip_month            !! Monthly precipitation
    REAL(r_std), DIMENSION(npts),INTENT(in)              :: vd                      !! Vegetation density

    ! Output variables
    REAL(r_std), DIMENSION(npts), INTENT(out)            :: dfm                     !! Daily fuel moisture (0,1), unitless. This fuel moisture should be understood
                                                                                    !! as (mass_water)/(mass_water + mass_dry_fuel)
    REAL(r_std), DIMENSION(npts), INTENT(out)            :: dfm_1hr                 !! Fuel moisture for the 1hr fuel, used only in calculating the dfm_lg_d1hr
    REAL(r_std), DIMENSION(npts), INTENT(out)            :: fdi_niacc               !! Fire danger index, unitless
    REAL(r_std), DIMENSION(npts), INTENT(out)            :: fdi_vpd_daily_mean
    REAL(r_std), DIMENSION(npts), INTENT(out)            :: fdi_vpd_daily_max
    REAL(r_std), DIMENSION(npts), INTENT(out)            :: fdi_vpd_mean_week
    REAL(r_std), DIMENSION(npts), INTENT(out)            :: fdi_vpd_max_week


    ! Local variables
    ! Alpha values are in proportion to the surface-to-volume ratios of different fuel classes.
    ! They're used as weighting factors to derive an overall drying speed of fuel bed.
    ! The alpha values for 1h/10h/100h fuel were used and were provided in Thonicke et al. 2010
    ! (Eq.6).  
    REAL(r_std), PARAMETER                           :: alpha_1hr=0.001            !! @tex $(Celcius^{-2})$ @endtex
    REAL(r_std), PARAMETER                           :: alpha_10hr=5.42E-5
    REAL(r_std), PARAMETER                           :: alpha_100hr=1.49E-5
    ! c_r (VPD precip-suppression coeff) is now the calibratable module getin cr_fire (CR_FIRE, default 2.0) -- fire knob
    REAL(r_std), DIMENSION(npts)                     :: alpha_fuel
    REAL(r_std), DIMENSION(npts)                     :: vartmp
    REAL(r_std), DIMENSION(npts)                     :: vartmp2
    REAL(r_std), DIMENSION(npts)                     :: Ft
    INTEGER                                          :: ivm 
    ! Initialise
    alpha_fuel(:)=zero
    fdi_niacc(:)=zero
    fdi_vpd_daily_mean(:)=zero
    fdi_vpd_daily_max(:)=zero
    fdi_vpd_mean_week(:)=zero
    fdi_vpd_max_week(:)=zero

    vartmp2(:)=zero
    ! Calcualte daily litter fuel moisture, using an empirical relationship. (Thonicke et al. 2010)
    ! (Eq. 6). The influence of 1000h fuel was ignored because we assume its moisture does not influence
    ! the likelihood of fire ignition.
    vartmp(:) = litterfuel_class_grid(:,ihour1) + litterfuel_class_grid(:,ihour10) + & 
                litterfuel_class_grid(:,ihour100)
    WHERE (vartmp(:) .GT. min_stomate) 
      alpha_fuel(:) = (alpha_1hr*litterfuel_class_grid(:,ihour1)+ &
                       alpha_10hr*litterfuel_class_grid(:,ihour10)+ &
                       alpha_100hr*litterfuel_class_grid(:,ihour100))/vartmp(:)
    ELSEWHERE
      ! There is almost no (fine) fuel to start a fire, we thus give a alpha_fuel as zero,
      ! so that fuel moisture would be 1 and fire danger index would be zero. I.e., there will
      ! be no fire.
      alpha_fuel(:) = zero
    ENDWHERE
    
    dfm(:)=EXP(-alpha_fuel(:)*ni_acc(:))  
                                          
    ! Moisture for the 1h fuel. Given that ni_acc is always equal to or larger than zero,
    ! dfm_1hr has a maximum value of 1 and will approaches to zero when ni_acc is very big.
    dfm_1hr(:) = EXP(-alpha_1hr * ni_acc(:)) 

        
    ! VPD-based FDI
    ! The vegetation density was excluded from this calculation.
    ! Because we want a pure climate-driven fire danger index.
    fdi_vpd_daily_mean(:) = max(zero,(vpd_daily_mean(:) * exp(-cr_fire * precip_month(:))) )
    fdi_vpd_daily_max(:) = max(zero,(vpd_daily_max(:) * exp(-cr_fire * precip_month(:))) )
    fdi_vpd_mean_week(:) = max(zero,(vpd_mean_week(:) * exp(-cr_fire * precip_month(:))) )
    fdi_vpd_max_week(:) = max(zero,(vpd_max_week(:) * exp(-cr_fire * precip_month(:))) )


    !vartmp2(:)=zero        
    ! Accumulate weighted fire danger based different vegetation types
    !DO ivm = 1,nvm
    !   IF (burnable(ivm)) THEN
    !      vartmp2(:)= vartmp2(:) + veget(:,ivm)*fdi_sf_pft(ivm)*ft(:)
    !   ENDIF
    !ENDDO
        
    ! Normalize by total vegetation cover if burnable vegetation
    ! exists 
    ! WHERE (vd(:) .LE. min_stomate)
    !    fdi(:)=zero
    ! ELSEWHERE
    !    fdi(:)=vartmp2(:)/vd(:)
    ! ENDWHERE

    ! Calculate fire danger index, following Eq.(8) in Thonicke et al. 2010.
    ! me_grid acts as the 'me' in Eq.(8)
    WHERE (me_grid(:) .LE. min_stomate)
      ! Moisture of extinction usually ranges 0.2 to 0.35 for different vegetative PFTs,
      ! so if we get a very small moisture of extinction for the whole grid cell, likely
      ! for this grid there is very little to almost no fuel. We can assume its fire danger
      ! index as zero. There is no risk for fire.
      fdi_niacc(:)=zero
    ELSEWHERE
      fdi_niacc(:)=max(zero,(1.0 - dfm(:)/me_grid(:) ))
    ENDWHERE

    END SUBROUTINE fire_danger_index




!! ================================================================================================================================
!! SUBROUTINE     : calculate_human_ignition
!!                                                                              
!>\BRIEF        Calculation human ignitions
!!                                                                              
!! DESCRIPTION   : Calculation human ignitions based on an emipical equation by Thonicke et al. 2010
!!                                                                              
!! REFERENCE(S)   : 
!! - K.Thonicke, A.Spessa, I.C. Prentice, S.P.Harrison, L.Dong, and
!! C.Carmona-Moreno, 2010: The influence of vegetation, fire spread and fire
!! behaviour on biomass burning and trace gas emissions: results from a
!! process-based model. Biogeosciences, 7, 1991-2011.
!!                                                                              
!! \n                                                                           
!_ ================================================================================================================================

    SUBROUTINE calculate_human_ignition(npts,popd,a_nd,AED,human_ign) 

    IMPLICIT NONE

    ! Input Variables
    INTEGER, INTENT(in)                              ::  npts  !! Domain size
    REAL(r_std), DIMENSION(npts), INTENT(in)         ::  popd  !! Population density (individuals km^{-2})
    REAL(r_std), DIMENSION(npts), INTENT(in)         ::  a_nd  !! Ignitions individual^{-1} day^{-1}
    REAL(r_std), DIMENSION(npts), INTENT(in)         ::  AED   !! Average edge distance (m)
    ! Output Variables
    REAL(r_std), DIMENSION(npts), INTENT(out)        ::  human_ign   !! Human ignitions (ignitions day^{-1} km^{-2})
    ! Local Variables
    REAL(r_std), DIMENSION(npts)                     :: fedge_ign    !! Fraction of edge area influenced by human ignition (unitless) (0-1)
    WHERE(popd(:).ge.min_stomate)
      ! compared with Eq(3) in Thonicke et al. (2010), we have to multiply further
      ! by 100 to change the unit from ignitions day^{-1}ha^{-1} to ignitions day^{-1}km^{-2}
      ! I don't know whether the '/100' at the end of Eq.(3) is for unit conversion or not.
      human_ign(:)=(EXP(-0.5*(popd(:)**0.5)))*a_nd(:)*a_nd_scale*popd(:)/100.0 * 100.
    ELSEWHERE
      human_ign(:)=zero
    ENDWHERE
    IF (ok_aed_humign) THEN
      ! AED Input
      ! Here, we assume that the probability of human ignitions increases as a direct function of the fragmentation perimeter
      ! as a fraction of total area, assuming that the fragmentation edge depth is 1m.  Thus an AED of 20m yields a grid area 
      ! coverage of 10% of a 0.5 degree resolution grid cell. This probability is scaled to the ignitions/person/km2/day as a 
      ! constant (/10000).  It results in significantly increased ignitions at low and high population density when fragmentation 
      ! is at a maximum (AED 20m or less), which decreases exponentially as fragmentation decreases (AED increases). 
      ! This is the clearest at high population densities, where the suppression effect of high population is counteracted by fragmentation.
      ! However, this does not mean that fire spread occurs (it just describes the probability that fire ignites). 
      ! The edge depth for human ignition (edge_distance_ignition=1.0 m) is changable.
      fedge_ign(:) = 1.0-((AED(:)-edge_distance_ignition)**2)/(AED(:)**2)
      IF (printlev_loc .GE. 5) THEN
        WRITE(numout,*) 'human_ign before:', human_ign(:)
        WRITE(numout,*) 'human_ign multiplier', edge_ignition_factor*fedge_ign(:)
      ENDIF
      WHERE(human_ign(:) .gt. min_stomate)
        human_ign(:)= human_ign(:)+edge_ignition_factor*fedge_ign(:)
      ELSEWHERE
        human_ign(:)=human_ign(:)
      ENDWHERE
      IF (printlev_loc .GE. 5) THEN
        WRITE(numout,*) 'human_ign after:', human_ign(:)
      ENDIF
    ENDIF

    END SUBROUTINE calculate_human_ignition

!! ================================================================================================================================
!! SUBROUTINE     : calculate_ground_lightn_ratio
!!                                                                              
!>\BRIEF        Calculation the ratio of cloud-to-ground lightning against the total lightning numbers,
!!              i.e., including both cloud-to-ground and cloud-to-cloud.
!!                                                                              
!! REFERENCE(S)   : 
!! - Prentice, S. A. and Mackerras, D.: The ratio of cloud to cloud- ground lightning flashes in 
!!   thunderstorms, J. Appl. Meteorol., 16, 545–550, 1977.
!! - F. Li, X. D. Zeng, and S. Levis.: A process-based fire parameterization of intermediate 
!!   complexity in a Dynamic Global Vegetation Model Biogeosciences, 9, 2761–2780, 2012.
!!                                                                              
!! \n                                                                           
!_ ================================================================================================================================

    SUBROUTINE calculate_ground_lightn_ratio(npts,lalo,lightning_cg_ratio) 

    ! Input Variables
    INTEGER, INTENT(in)                              :: npts   !! Domain size
    REAL(r_std), DIMENSION(npts,2), INTENT(in)       :: lalo   !! Vector of latitude and longitudes (beware of the order !)
    ! Output Variables
    REAL(r_std), DIMENSION(npts), INTENT(out)        :: lightning_cg_ratio   !! Ratio of cloud-to-ground to
                                                                             !! total lightnings.
    REAL(r_std), DIMENSION(npts)                     :: llat
    INTEGER :: ipts

    ! latitude
    llat = lalo(:,1)

    ! The original Prentice & Mackerras paper only fits the empirical equation for up to 60 degree.
    ! Hence we set 60 degree as a saturating degree.
    DO ipts = 1,npts
      IF (llat(ipts) .GT. 60) THEN
        llat(ipts) = 60.
      ELSE
        IF (llat(ipts) .LT. -60) THEN
          llat(ipts) = -60.
        ENDIF
      ENDIF
    ENDDO

    llat  = lalo(:,1)*pi/180.
    ! hour angle is zero at noon, hence its cos = 1 and the term can be
    ! dropped from the equation
    lightning_cg_ratio(:) = 1/(5.16+2.16*COS(3*llat))
    !WRITE(numout,*) 'lightning_cg_ratio: ',lightning_cg_ratio(:)
    END SUBROUTINE calculate_ground_lightn_ratio

!! ================================================================================================================================
!! SUBROUTINE     : human_suppression_func
!!                                                                              
!>\BRIEF        Account for human suppression for lightning ignited fires.
!!                                                                              
!! REFERENCE(S)   : 
!! Li et al. (2011) A process-based fire parameterization of intermediate 
!! complexity in a Dynamic Global Vegetation Model.
!!                                                                              
!! FLOWCHART   : None                                                           
!! \n                                                                           
!_ ================================================================================================================================

    SUBROUTINE human_suppression_func(npts,popd,human_suppression)

     IMPLICIT NONE

     ! Input Variables
     INTEGER, INTENT(in)                              :: npts
     REAL(r_std), DIMENSION(npts), INTENT(in)         ::  popd  ! Population  density (individuals km^{-2})
     ! Output variables
     REAL(r_std), DIMENSION(npts), INTENT(out)        ::  human_suppression

     WHERE(popd(:) .GT. min_stomate) 
       ! Source: Li et al. (2011) A process-based fire parameterization of intermediate complexity in a Dynamic Global
       ! Vegetation Model. Page 2766
       human_suppression(:) = 0.99-0.98*EXP(-0.025*popd(:)) 
     ELSEWHERE
       human_suppression(:)=zero
     ENDWHERE

    END SUBROUTINE human_suppression_func       


!! ================================================================================================================================
!! SUBROUTINE     : rate_of_spread
!!                                                                              
!>\BRIEF        Calculation fire spread rate, This is the application of Eq.(9) in
!!     Thonicke et al. (2010)
!!                                                                              
!! REFERENCE(S)   : 
!! - K.Thonicke, A.Spessa, I.C. Prentice, S.P.Harrison, L.Dong, and
!! C.Carmona-Moreno, 2010: The influence of vegetation, fire spread and fire
!! behaviour on biomass burning and trace gas emissions: results from a
!! process-based model. Biogeosciences, 7, 1991-2011.
!!                                                                              
!! FLOWCHART   : None                                                           
!! \n                                                                           
!_ ================================================================================================================================

    SUBROUTINE rate_of_spread(npts, netfuel_grid, wind_forward, sigma,   &
                              dfm, me_grid, H, bulkdensity_grid,         &
                              ros_f, var_gamma, wetness)

    IMPLICIT NONE    

    ! Input variables
    INTEGER, INTENT(in)                              :: npts
    REAL(r_std), DIMENSION(npts), INTENT(in)         :: wind_forward 
    REAL(r_std), DIMENSION(npts), INTENT(in)         :: bulkdensity_grid  !! Fuel bulk density weighted by mass of dead fine fuel and live grass.
                                                                          !! (kg dry mass m^{-3}) 
    REAL(r_std), DIMENSION(npts), INTENT(in)         :: sigma             !! Mass weighted mean fuel surface-area-to-volume ratio. (cm^{-1})
    REAL(r_std), DIMENSION(npts), INTENT(in)         :: netfuel_grid      !! Total net fuel (mineral content corrected) by summing the 
                                                                          !! dead fuel and live grass. (kg dry mass m^{-2})
    REAL(r_std), DIMENSION(npts), INTENT(in)         :: dfm               !! Daily fuel moisture (0,1), unitless.
    REAL(r_std), DIMENSION(npts), INTENT(in)         :: me_grid           !! The moisture of extinction for fuel over the whole grid. (0,1) (unitless)
    REAL(r_std), INTENT(in)                          :: H                 !! heat content of the fuel (18000kJ/kg)

    ! Output variables       
    REAL(r_std), DIMENSION(npts), INTENT(out)        :: ros_f             !! Forward fire spread rate (m min^{-1}), ROS_{f,surface} in Thonicke et al. 2010 
    REAL(r_std), DIMENSION(npts), INTENT(out)        :: var_gamma
    REAL(r_std), DIMENSION(npts), INTENT(out)        :: wetness           !! Ratio of dead fuel moisture to the moisture of extinction. (unitless)

    ! Local variables
    REAL(r_std), DIMENSION(npts)                     :: beta              !! Packing ratio (unitless). The ratio of fuel_bulk_density against 
                                                                          !! oven-dry particle density. It measures how lossely the fuel is packed.
    REAL(r_std), DIMENSION(npts)                     :: beta_op           !! Optimum packing ratio (unitless)
    REAL(r_std), DIMENSION(npts)                     :: q_ig              !! Heat of pre-ignition, i.e., the amount of heat required to ignite a
                                                                          !! given fuel mass (kJ kg^{-1})
    REAL(r_std), DIMENSION(npts)                     :: eps               !! Effective heating number, i.e. the proportion of a fuel particle that 
                                                                          !! is heated to ignition temperature at the time flaming combustion starts.
                                                                          !! (unitless)
    REAL(r_std), DIMENSION(npts)                     :: phi_wind          !! A multiplier that accounts for the effect of wind in increasing the 
                                                                          !! effectiveness of propogating flux heating. cf. Eq. (9) in Thoncke et al. 2010
    REAL(r_std), DIMENSION(npts)                     :: xi                !! Propagating flux ratio (unitless). It measures the proportion of energy
                                                                          !! released during fuel combustion used to heat adjacent fuel. 
    REAL(r_std), DIMENSION(npts)                     :: var_B             !! Temporary variable. Eq. (A7)
    REAL(r_std), DIMENSION(npts)                     :: var_C             !! Temporary variable. Eq. (A8)
    REAL(r_std), DIMENSION(npts)                     :: var_E             !! Temporary variable. Eq. (A9)
    REAL(r_std), DIMENSION(npts)                     :: var_A             !! Temporary variable. cf. Table A1 Row 4
    REAL(r_std), DIMENSION(npts)                     :: gamma_max         !! Maximum reaction velocity (min^{-1})
    REAL(r_std), DIMENSION(npts)                     :: gamma_aptr        !! Optimum reaction velocity (min^{-1})
    REAL(r_std), DIMENSION(npts)                     :: ir                !! Reaction intensity, i.e., energy release rate per unit area of fire front
                                                                          !! (kJ m^{-2} min^{-1})
    REAL(r_std), DIMENSION(npts)                     :: ratio_beta        !! Ratio between fuel packing ratio and optimum packing ratio
    REAL(r_std), DIMENSION(npts)                     :: moist_damp        !! Moisture dampening coefficient to adjust reaction intensity. (unitless)
    REAL(r_std), DIMENSION(npts)                     :: vartmp


    REAL(r_std), PARAMETER :: MINER_DAMP=0.41739 ! Mineral dampening coefficient; Thonicke 2010 Table A1 Row 7
    REAL(r_std), PARAMETER :: part_dens=513.0    ! Oven-dry fuel partical density (kg m-3); 
                                                 ! Thonicke 2010 Table A1 Row 3
    INTEGER :: m
    INTEGER :: itestpft,itestipts

    itestipts=5

    ! Initialise
    xi(:)=zero

    beta(:) = bulkdensity_grid(:) / part_dens ! packing ratio
    beta_op(:) = 0.200395*(sigma(:)**(-0.8189)) ! Eq.(A6)
    ratio_beta(:) = beta(:) / beta_op(:)  ! prepare for Eq.(A5)

    q_ig(:)=581.0+2594.0*dfm(:)  ! Heat of pre-ignition; Eq. (A4).
    eps(:)=EXP(-4.528/sigma(:)) ! Effective heating number; Eq. (A3)

    !! Calculating the correction factor due to wind effect on the propogating heat flux
    var_B(:)=0.15988*(sigma(:)**0.54)  ! Eq.(A7)
    var_C(:)=7.47*(EXP(-0.8711*(sigma(:)**0.55)))  ! Eq.(A8)
    var_E(:)=0.7515*(EXP(-0.01094*sigma(:)))   ! Eq.(A9)

    WHERE( ratio_beta(:) .GT. 0.00001 .AND. wind_forward(:) .GT. min_stomate )
      phi_wind(:)=var_C(:)*(wind_forward(:)**var_B(:))*(ratio_beta(:)**(-var_E(:)))  ! Eq.(A5)
    ENDWHERE
    

    ! chaoyue: I guess this code is given to me by Patricia. But I don't know the source.
    WHERE(ratio_beta(:) .LE. 0.00001)
      phi_wind(:)=var_C(:)*(wind_forward(:)**var_B(:))*0.00001  !??UNIMPORTANT not sure this will pose a problem 
                                                        !but when ratio_beta --> 0+, ratio_beta**(-e) can be very big, which 
                                                        !is on the contrary with the simplied form.
    ENDWHERE
    
    WHERE(wind_forward(:) .LE. min_stomate)
      phi_wind(:)=zero
    ENDWHERE

    !! Calculate propogating flux ratio
    WHERE (sigma(:) .LE. 0.00001) 
      xi(:)=zero
    ELSEWHERE
      xi(:) = (EXP((0.792 + 3.7597 * (sigma(:)**0.5)) * &  ! Eq.(A2)
         (beta(:) + 0.1))) / (192 + 7.9095 * sigma(:))
    ENDWHERE

    !! Calculate optimum reaction intensity
    var_A(:)=8.9033*(sigma(:)**(-0.7913))  ! Source: Table A1 Row 4
    WHERE ( var_A(:) .LE. 0.00001 .OR. ratio_beta(:) .LE. 0.00001 )
       vartmp(:)=zero
    ELSEWHERE
       vartmp(:)=EXP(var_A(:)*(1.0-ratio_beta(:)))
    ENDWHERE
      
    gamma_max(:)=1.0/(0.0591+2.926*(sigma(:)**(-1.5)))  !Thonicke 2010 Table A1 Row 2

    WHERE(ratio_beta(:) .GT. 0.0001)
      gamma_aptr(:)=gamma_max(:)*(ratio_beta(:)**var_A(:))*vartmp(:)
    ELSEWHERE
      gamma_aptr(:)=gamma_max(:)*(0.0001)*vartmp(:)
    ENDWHERE
    
    !! Calculate reaction intensity
    ! wn/me in Eq. row 6 in Table A1 in Thonicke 2010.
    WHERE (me_grid(:) .GT. min_stomate)
      wetness(:)=dfm(:)/me_grid(:)
    ELSEWHERE
      wetness(:)=1.0
    ENDWHERE

    WHERE (wetness(:) .GT. min_stomate)
      moist_damp(:)=max(zero,(1.0-(2.59*wetness(:))+ &
        (5.11*(wetness(:)**2.0))-(3.52*(wetness(:)**3.0))))
    ELSEWHERE
      moist_damp(:)=zero
    ENDWHERE

    ir(:)=gamma_aptr(:) * netfuel_grid(:) * H * moist_damp(:) * MINER_DAMP !Eq.(A1)
    !WRITE(numout,*) 'CheckoutC', gamma_aptr(itestipts),netfuel_grid(itestipts),moist_damp(itestipts),&
    !me_grid(itestipts),wetness(itestipts),dfm(itestipts)

    !! For use in calculating fire residence time (tau_l), further used in calculating tree mortality.
    ! Eq. (43) in SPITFIRE technical document.
    var_gamma(:)=gamma_aptr(:)*moist_damp(:)*MINER_DAMP
      
    !! Calculate forward fire spread rate
    WHERE ( bulkdensity_grid(:) .LE. min_stomate .OR. &
            eps(:) .LE. min_stomate .OR. q_ig(:) .LE. min_stomate )
      ros_f(:)=zero
    ELSEWHERE
      ros_f(:)=(ir(:) * xi(:) * (1.0 + phi_wind(:))) / &
          (bulkdensity_grid(:) * eps(:) * q_ig(:))
    ENDWHERE

!    WRITE(numout,*) 'CheckoutD', ros_f(itestipts),ir(itestipts),phi_wind(itestipts),&
!   bulkdensity_grid(itestipts),eps(itestipts),q_ig(itestipts)

    END SUBROUTINE rate_of_spread

!! ================================================================================================================================
!! SUBROUTINE     : combustion_fraction
!!                                                                              
!>\BRIEF        Calculation fuel combustion fractions during fire.
!!                                                                              
!! DESCRIPTION   : None                            
!!                                                                              
!! REFERENCE(S)   : 
!! - K.Thonicke, A.Spessa, I.C. Prentice, S.P.Harrison, L.Dong, and
!! C.Carmona-Moreno, 2010: The influence of vegetation, fire spread and fire
!! behaviour on biomass burning and trace gas emissions: results from a
!! process-based model. Biogeosciences, 7, 1991-2011.
!!                                                                              
!_ ================================================================================================================================

    SUBROUTINE combustion_fraction(npts, livegrass_grid, litterfuel_class_grid,    &
                 dfm_1hr, dfm_livegrass, wetness, me_grid, cf_litterfuel,          &
                 cf_lg, AED)


    IMPLICIT NONE

    integer :: d,ivm,k

    !! Input variables
    INTEGER, INTENT(in)                                         :: npts
    REAL(r_std), DIMENSION(npts), INTENT(in)                    :: livegrass_grid   !! Area-weighted sum of live aboveground biomass of burnable herbaceous PFTs. (gC m^{-2})

    REAL(r_std), DIMENSION(npts,nhour), INTENT(in)              :: litterfuel_class_grid  !! Gridcell value of aboveground litter fuel of different fuel classes
                                                                                          !! (gC m^{-2}) over burnable ground area
    REAL(r_std), DIMENSION(npts), INTENT(in)                    :: dfm_1hr          !! Daily fuel moisture for 1hr fuel, unitless, (0,1)
    REAL(r_std), DIMENSION(npts), INTENT(in)                    :: dfm_livegrass    !! Daily live grass fuel moisture, unitless, (0,1)
    REAL(r_std), DIMENSION(npts), INTENT(in)                    :: wetness          !! Ratio of dead fuel moisture to the moisture of extinction. (unitless)
    REAL(r_std), DIMENSION(npts), INTENT(in)                    :: me_grid          !! The moisture of extinction for fuel over the whole grid. (0,1) (unitless)
    REAL(r_std), DIMENSION(npts), INTENT(in)                    :: AED              !! Average edge distance (m)


    !! Output variables
    REAL(r_std), DIMENSION(npts,nhour), INTENT(out)             :: cf_litterfuel    !! Combustion fractions for different surface fuel categories. (unitless) (0,1)
    REAL(r_std), DIMENSION(npts), INTENT(out)                   :: cf_lg            !! Combustion fraction for live grass fuel. (unitless) (0,1)

    !! Local variables
    REAL(r_std), DIMENSION(npts)                                :: dfm_lg_d1hr      !! Fuel moisture for 1hr dead fuel and livegrass moisture combined. (unitless) (0,1)
    REAL(r_std), DIMENSION(npts)                                :: cf_1hr           !! Combustion fraction for 1h fuel. (unitless) (0,1)
    REAL(r_std), DIMENSION(npts)                                :: cf_10hr          !! Combustion fraction for 10h fuel. (unitless) (0,1)
    REAL(r_std), DIMENSION(npts)                                :: cf_100hr         !! Combustion fraction for 100h fuel. (unitless) (0,1)
    REAL(r_std), DIMENSION(npts)                                :: cf_1000hr        !! Combustion fraction for 1000h fuel. (unitless) (0,1)
    REAL(r_std), DIMENSION(npts)                                :: wetness_lg       !! Wetness of live grass. (unitless) (0-1)
    REAL(r_std), DIMENSION(npts)                                :: wetness_1hr      !! Wetness of 1h fuel. (unitless) (0-1) 
    REAL(r_std), PARAMETER :: me_livegrass=0.2                                      !! Moisture of extinction for live grass (unitless) 
    REAL(r_std), DIMENSION(npts)                                :: wetness_edge     !! Wetness modified by edge effects. (unitless)
    REAL(r_std), DIMENSION(npts)                                :: fedge_fuel       !! Fraction of edge area where fuel moisture is modified (unitless) (0-1)


    ! Initialise
    wetness_lg(:)=zero
    wetness_1hr(:)=zero
    wetness_edge(:)=wetness(:)
    cf_1hr(:)=zero
    cf_10hr(:)=zero
    cf_100hr(:)=zero
    cf_1000hr(:)=zero
    !! Combustion fraction for live grass and 1h dead fuel.

    ! Caculate wetness for live grass and 1h dead fuel. The wetness of 1h dead fuel is influenced by
    ! the moisture of live grass.
    ! dfm_livegrass is the $\omega_{lg}$ in Eq.(B2) and Eq.(B3)
    WHERE (litterfuel_class_grid(:,ihour1) .GT. min_stomate)
      dfm_lg_d1hr(:)=dfm_1hr(:)+(dfm_livegrass(:)*(livegrass_grid(:)/litterfuel_class_grid(:,ihour1))) !Eq.(B3)
    ELSEWHERE
      dfm_lg_d1hr(:)=1.0
    ENDWHERE
 
    WHERE (me_grid(:) .GT. min_stomate)
      wetness_lg(:)=dfm_livegrass(:)/me_livegrass
      wetness_1hr(:)=dfm_lg_d1hr(:)/me_grid(:)
    ELSEWHERE
      wetness_lg(:)=1.0
      wetness_1hr(:)=1.0
    ENDWHERE

    ! Implementation based on Eq. (9) in Bowring et al. (2024).
    ! Equation (10) is not used here, as modifying only the threshold of the piecewise
    ! function would introduce a discontinuity.
    ! We decreaded the fuel moisture but only in the combustion_fraction. 
    ! In future, should also be consistently adjusted in the computation of FDI and ROS.
    ! We assume a edge depth of 20 m , with an effective mean value of 10 m.
    ! The edge depth for changing wetness (edge_distance_wetness=20.0 m) is changable.
    IF (ok_aed_fuel) THEN
      fedge_fuel(:) = 1.0-((AED(:)-edge_distance_wetness/2)**2)/(AED(:)**2)
      IF (printlev_loc .GE. 5) THEN
        WRITE(numout,*) 'wetness_lg before:', wetness_lg(:)
        WRITE(numout,*) 'wetness_1hr before:', wetness_1hr(:)
        WRITE(numout,*) 'wetness_edge before:', wetness_edge(:)
        WRITE(numout,*) 'wetness multiplier', (1.0-(0.25/2)*fedge_fuel(:))
      ENDIF
      wetness_lg(:) = wetness_lg(:)*(1.0-(0.25/2)*fedge_fuel(:))
      wetness_1hr(:) = wetness_1hr(:)*(1.0-(0.25/2)*fedge_fuel(:))
      wetness_edge(:) = wetness(:)*(1.0-(0.25/2)*fedge_fuel(:))
      IF (printlev_loc .GE. 5) THEN
        WRITE(numout,*) 'wetness_lg after:', wetness_lg(:)
        WRITE(numout,*) 'wetness_1hr after:', wetness_1hr(:)
        WRITE(numout,*) 'wetness_edge after:', wetness_edge(:)
      ENDIF
    ENDIF

    ! Combustion fraction for live grass biomass, by applying Eq. (B1)
    WHERE (wetness_lg(:) .LE. 0.18)
       cf_lg(:)=1.0
    ELSEWHERE ( wetness_lg(:) .GT. 0.18 .AND. wetness_lg(:) .LE. 0.73 )
       cf_lg(:)=1.11-0.62*wetness_lg(:)
    ELSEWHERE ( wetness_lg(:) .GT. 0.73 .AND. wetness_lg(:) .LE. 1.0)
       cf_lg(:)=2.45-2.45*wetness_lg(:)
    ELSEWHERE (wetness_lg(:) .GT. 1.0)
       cf_lg(:)=zero
    ENDWHERE
    
    ! Combustion fraction for 1h fuel, applying Eq. (B1)
    WHERE (wetness_1hr(:) .LE. 0.18)
       cf_1hr(:)=1.0  
    ELSEWHERE ( wetness_1hr(:) .GT. 0.18 .AND. wetness_1hr(:) .LE. 0.73 )
       cf_1hr(:)=1.11-0.62*wetness_1hr(:)
    ELSEWHERE ( wetness_1hr(:) .GT. 0.73 .AND. wetness_1hr(:) .LE. 1.0 )
       cf_1hr(:)=2.45-2.45*wetness_1hr(:)
    ELSEWHERE (wetness_1hr(:) .GT. 1.0)
       cf_1hr(:)=zero
    ENDWHERE

    WHERE (cf_1hr(:).GT.0.9)
      cf_1hr(:) = 0.9
    ENDWHERE
    
    !! Combustion fraction for 10h fuel, applying Eq. (B4)
    IF (.not. ok_aed_fuel) THEN
      WHERE (wetness(:) .LE. 0.12)
        cf_10hr(:)=1.0
      ELSEWHERE ( wetness(:) .GT. 0.12 .AND. wetness(:) .LE. 0.51 )
        cf_10hr(:)=1.09-0.72*wetness(:)
      ELSEWHERE ( wetness(:) .GT. 0.51 .AND. wetness(:) .LE. 1.0 )
        cf_10hr(:)=1.47-1.47*wetness(:)
      ELSEWHERE (wetness(:) .GT. 1.0)
        cf_10hr(:)=zero
      ENDWHERE
    ELSE ! ok_aed_fuel
      WHERE (wetness_edge(:) .LE. 0.12)
        cf_10hr(:)=1.0
      ELSEWHERE ( wetness_edge(:) .GT. 0.12 .AND. wetness_edge(:) .LE. 0.51 )
        cf_10hr(:)=1.09-0.72*wetness(:)
      ELSEWHERE ( wetness_edge(:) .GT. 0.51 .AND. wetness_edge(:) .LE. 1.0 )
        cf_10hr(:)=1.47-1.47*wetness_edge(:)
      ELSEWHERE (wetness_edge(:) .GT. 1.0)
        cf_10hr(:)=zero
      ENDWHERE
    ENDIF ! IF (.not. ok_aed_fuel) THEN

    
    ! I forget why we need this. This is a little arbitrary though.
    WHERE (cf_10hr(:) .GT. 0.9)
      cf_10hr(:) = 0.9
    ENDWHERE
    
    !! Combustion fraction for 100h fuel
    WHERE (wetness(:) .LE. 0.38)
      cf_100hr(:)=0.98-0.85*wetness(:)
    ELSEWHERE ( wetness(:) .GT. 0.38 .AND. wetness(:) .LE. 1.0 )
      cf_100hr(:)=1.06-1.06*wetness(:)
    ELSEWHERE (wetness(:) .GT. 1.0)
      cf_100hr(:)=zero
    ENDWHERE

    !! Combustion fraction for 1000h fuel
    cf_1000hr(:) = -0.8*wetness(:)+0.8
    
    ! This correction is added in Yue et al. (2014), GMD Sect. 2.2 on Page
    ! 2750.It seems have some problems. Because the fire_max_cf_100hr is at the
    ! PFT levels. But the cf_100hr only at the pixel scale. Hence, PZ here take
    ! the average value for fire_max_cf_100hr and fire_max_cf_1000hr.

    !WHERE (cf_100hr(:) .GT. fire_max_cf_100hr(ivm))
    !  cf_100hr(:) = fire_max_cf_100hr(ivm)
    !ENDWHERE

    !WHERE (cf_1000hr(:) .GT. fire_max_cf_1000hr(ivm))
    !  cf_1000hr(:) = fire_max_cf_1000hr(ivm)
    !ENDWHERE

    WHERE (cf_100hr(:) .GT. 0.7)
      cf_100hr(:) = 0.7
    ENDWHERE
    
    WHERE (cf_1000hr(:) .GT. 0.4)
      cf_1000hr(:) = 0.4
    ENDWHERE    

    WHERE (cf_1000hr(:) .LT. 0.0)
      cf_1000hr(:) = 0.0
    ENDWHERE    
    
    ! Write to variable
    cf_litterfuel(:,ihour1) = cf_1hr(:)   
    cf_litterfuel(:,ihour10) = cf_10hr(:)   
    cf_litterfuel(:,ihour100) = cf_100hr(:)   
    cf_litterfuel(:,ihour1000) = cf_1000hr(:)   

    END SUBROUTINE combustion_fraction


!! ================================================================================================================================
!! SUBROUTINE     : deforestation_fire
!!                                                                              
!>\BRIEF       Simulate deforestation processes from land cover change.
!!                                                                              
!! DESCRIPTION   : None                            
!!                                                                              
!! RECENT CHANGE(S): None                                                       
!!                                                                              
!! RETURN VALUE :  Nonw
!! 
!! REFERENCE(S)   : 
!!                                                                              
!! FLOWCHART   : None                                                           
!! \n                                                                           
!_ ================================================================================================================================

    SUBROUTINE deforestation_fire(npts,fdi,lcc,biomass,fuel_1hr,&
        fuel_10hr,fuel_100hr,fuel_1000hr,bafrac_deforest,&
        emideforest_fuel_1hr,emideforest_fuel_10hr,emideforest_fuel_100hr,&
        emideforest_fuel_1000hr,emideforest_biomass)

    IMPLICIT NONE

    !! Input variables
    ! Domain size
    INTEGER, INTENT(in)                                     :: npts
    REAL(r_std), DIMENSION(npts)                            :: fdi             !! daily fire danger index,climatic fire risk, unitless (0,1)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)            :: lcc               !! gross forest cover loss
    REAL(r_std), DIMENSION(npts,nvm,nparts), INTENT(in)     :: biomass
    REAL(r_std), DIMENSION(npts,nvm,nlitt), INTENT(in)      :: fuel_1hr          !! [gC^{m-2}]
    REAL(r_std), DIMENSION(npts,nvm,nlitt), INTENT(in)      :: fuel_10hr         !! [gC^{m-2}]
    REAL(r_std), DIMENSION(npts,nvm,nlitt), INTENT(in)      :: fuel_100hr        !! [gC^{m-2}]
    REAL(r_std), DIMENSION(npts,nvm,nlitt), INTENT(in)      :: fuel_1000hr       !! [gC^{m-2}]
     
    !! Output variables
    REAL(r_std), DIMENSION(npts,nvm),INTENT(inout)          :: bafrac_deforest         !! Deforestation fire burned fraction,unitless
    REAL(r_std), DIMENSION(npts,nvm,nlitt),INTENT(inout)    :: emideforest_fuel_1hr    !! fuel consumption for 1hr fuel (gCm^{-2})
    REAL(r_std), DIMENSION(npts,nvm,nlitt),INTENT(inout)    :: emideforest_fuel_10hr   !! 
    REAL(r_std), DIMENSION(npts,nvm,nlitt),INTENT(inout)    :: emideforest_fuel_100hr  !! 
    REAL(r_std), DIMENSION(npts,nvm,nlitt),INTENT(inout)    :: emideforest_fuel_1000hr !! 
    REAL(r_std), DIMENSION(npts,nvm,nparts),INTENT(inout)   :: emideforest_biomass     !! live biomass emissions from deforestation fires [gCm^{-2}]
    
    !! Local variables
    REAL(r_std), DIMENSION(npts)                            :: deforest_factor         !!Factor to scale the deforestation fire fraction, unitless
    REAL(r_std), DIMENSION(npts)                            :: vartmp
    REAL(r_std), DIMENSION(npts)                            :: cc_litter               !! Fire consumption fraction for litter, leaf,fruit and reserve biomass
    REAL(r_std), DIMENSION(npts)                            :: cc_wood                 !! Fire consumption fraction for heart and sapwood
    INTEGER(i_std)                                          :: ilit,ivm                  !! Indices (unitless)


    deforest_factor(:) = 0.017
    cc_litter(:) = 0.7 + 0.3*fdi(:)
    cc_wood(:) = 0.3 + fdi(:)*0.3  

    DO ivm=2,nvm
      IF (is_tree(ivm)) THEN
        WHERE (lcc(:,ivm) > min_stomate)
          !The equation of y=0.58*x-0.006 is derived by regressing GFED3 reported
          !deforestation fire burned fraction to the gross forest loss by Hansen 2010
          !for the closed canopy Amazonian forest
          vartmp(:) = MAX(0.0001,lcc(:,ivm)*0.58-0.006)
          bafrac_deforest(:,ivm) = deforest_factor(:) * vartmp(:) * fdi(:) 

          !We calculate fire emissions from litter and biomass
          emideforest_fuel_1hr(:,ivm,imetabolic) = cc_litter(:) * bafrac_deforest(:,ivm) * fuel_1hr(:,ivm,imetabolic)
          emideforest_fuel_1hr(:,ivm,istructural) = cc_litter(:) * bafrac_deforest(:,ivm) * fuel_1hr(:,ivm,istructural)
          emideforest_fuel_10hr(:,ivm,imetabolic) = cc_litter(:) * bafrac_deforest(:,ivm) * fuel_10hr(:,ivm,imetabolic)
          emideforest_fuel_10hr(:,ivm,istructural) = cc_litter(:) * bafrac_deforest(:,ivm) * fuel_10hr(:,ivm,istructural)
          emideforest_fuel_100hr(:,ivm,imetabolic) = cc_litter(:) * bafrac_deforest(:,ivm) * fuel_100hr(:,ivm,imetabolic)
          emideforest_fuel_100hr(:,ivm,istructural) = cc_litter(:) * bafrac_deforest(:,ivm) * fuel_100hr(:,ivm,istructural)
          emideforest_fuel_1000hr(:,ivm,imetabolic) = cc_litter(:) * bafrac_deforest(:,ivm) * fuel_1000hr(:,ivm,imetabolic)
          emideforest_fuel_1000hr(:,ivm,istructural) = cc_litter(:) * bafrac_deforest(:,ivm) * fuel_1000hr(:,ivm,istructural)

          emideforest_biomass(:,ivm,icarbres) = cc_litter(:) * bafrac_deforest(:,ivm) * biomass(:,ivm,icarbres)
          emideforest_biomass(:,ivm,ifruit) = cc_litter(:) * bafrac_deforest(:,ivm) * biomass(:,ivm,ifruit)
          emideforest_biomass(:,ivm,ileaf) = cc_litter(:) * bafrac_deforest(:,ivm) * biomass(:,ivm,ileaf)
          emideforest_biomass(:,ivm,iheartabove) = cc_wood(:) * bafrac_deforest(:,ivm) * biomass(:,ivm,iheartabove)
          emideforest_biomass(:,ivm,isapabove) = cc_wood(:) * bafrac_deforest(:,ivm) * biomass(:,ivm,isapabove)
        ENDWHERE 
      END IF
    END DO
   
   END SUBROUTINE deforestation_fire
   
!! ================================================================================================================================
!! SUBROUTINE     : deforestation_fire_proxy
!!                                                                              
!>\BRIEF       Simulate deforestation processes from land cover change.
!!                                                                              
!! DESCRIPTION   : None                            
!!                                                                              
!! RECENT CHANGE(S): None                                                       
!!                                                                              
!! RETURN VALUE :  Nonw
!! 
!! REFERENCE(S)   : 
!!                                                                              
!! FLOWCHART   : None                                                           
!! \n                                                                           
!_ ================================================================================================================================

    SUBROUTINE deforestation_fire_proxy(npts,fdi,lcc,                     &
        deforest_biomass_remain,                     &
        def_fuel_1hr_remain,def_fuel_10hr_remain,                           &
        def_fuel_100hr_remain,def_fuel_1000hr_remain,                       &
        bafrac_deforest, emideforest_fuel_1hr,emideforest_fuel_10hr,        &
        emideforest_fuel_100hr, emideforest_fuel_1000hr,emideforest_biomass)


    IMPLICIT NONE

    !! Input variables
    ! Domain size
    INTEGER, INTENT(in)                                     :: npts
    REAL(r_std), DIMENSION(npts)                            :: fdi             !! daily fire danger index,climatic fire risk, unitless (0,1)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)            :: lcc               !! gross forest cover loss
    REAL(r_std), DIMENSION(npts,nvm,nlitt), INTENT(in)        :: def_fuel_1hr_remain     !! [gC^{m-2}]
    REAL(r_std), DIMENSION(npts,nvm,nlitt), INTENT(in)        :: def_fuel_10hr_remain    !! [gC^{m-2}]
    REAL(r_std), DIMENSION(npts,nvm,nlitt), INTENT(in)        :: def_fuel_100hr_remain   !! [gC^{m-2}]
    REAL(r_std), DIMENSION(npts,nvm,nlitt), INTENT(in)        :: def_fuel_1000hr_remain  !! [gC^{m-2}]
    REAL(r_std), DIMENSION(npts,nvm,nparts), INTENT(in)      :: deforest_biomass_remain 
     
    !! Output variables
    REAL(r_std), DIMENSION(npts,nvm),INTENT(inout)          :: bafrac_deforest         !! Deforestation fire burned fraction,unitless
    REAL(r_std), DIMENSION(npts,nvm,nlitt),INTENT(inout)    :: emideforest_fuel_1hr    !! fuel consumption for 1hr fuel (gCm^{-2})
    REAL(r_std), DIMENSION(npts,nvm,nlitt),INTENT(inout)    :: emideforest_fuel_10hr   !! 
    REAL(r_std), DIMENSION(npts,nvm,nlitt),INTENT(inout)    :: emideforest_fuel_100hr  !! 
    REAL(r_std), DIMENSION(npts,nvm,nlitt),INTENT(inout)    :: emideforest_fuel_1000hr !! 
    REAL(r_std), DIMENSION(npts,nvm,nparts),INTENT(inout)   :: emideforest_biomass     !! live biomass emissions from deforestation fires [gCm^{-2}]
    
    !! Local variables
    REAL(r_std), DIMENSION(npts,nvm)                        :: bafrac_net              !! Deforestation fire burned fraction,unitless
    REAL(r_std), DIMENSION(npts)                            :: deforest_factor         !!Factor to scale the deforestation fire fraction, unitless
    REAL(r_std), DIMENSION(npts)                            :: vartmp
    REAL(r_std), DIMENSION(npts)                            :: cc_litter               !! Fire consumption fraction for litter, leaf,fruit and reserve biomass
    REAL(r_std), DIMENSION(npts)                            :: cc_wood                 !! Fire consumption fraction for heart and sapwood
    INTEGER(i_std)                                          :: ilit,ivm                  !! Indices (unitless)


    ! the 7.2 is tested. The test result could be found at 
    ! /home/cyue/Documents/ORCHIDEE/TEST_MICTv6bis/Deforestation_test/Deforestation_fire_annual_2001_2010_55e0e3e_T2.nc
    ! with the notebook be found at: /home/cyue/Documents/ORCHIDEE/TEST_MICTv6bis/Test_tag_MICTv6bisT2
    ! the readme file is at: /home/orchidee01/ychao/MICTv6/bin/readme.txt
    deforest_factor(:) = 0.017*7.2
    cc_litter(:) = 0.7 + 0.3*fdi(:)
    cc_wood(:) = 0.3 + fdi(:)*0.3  

    DO ivm=2,nvm
      IF (is_tree(ivm)) THEN
        WHERE (lcc(:,ivm) > min_stomate)
          !The equation of y=0.58*x-0.006 is derived by regressing GFED3 reported
          !deforestation fire burned fraction to the gross forest loss by Hansen 2010
          !for the closed canopy Amazonian forest
          vartmp(:) = MAX(0.0001,lcc(:,ivm)*0.58-0.006)
          bafrac_deforest(:,ivm) = deforest_factor(:) * vartmp(:) * fdi(:) 
          bafrac_net(:,ivm)=bafrac_deforest(:,ivm)/lcc(:,ivm)

          !We calculate fire emissions from litter and biomass
          !note here because cc_litter*bafrac_net < 1, the fuel will never be 
          !exhausted.
          emideforest_fuel_1hr(:,ivm,imetabolic) = cc_litter(:) * bafrac_net(:,ivm) * def_fuel_1hr_remain(:,ivm,imetabolic)
          emideforest_fuel_1hr(:,ivm,istructural) = cc_litter(:) * bafrac_net(:,ivm) * def_fuel_1hr_remain(:,ivm,istructural)
          emideforest_fuel_10hr(:,ivm,imetabolic) = cc_litter(:) * bafrac_net(:,ivm) * def_fuel_10hr_remain(:,ivm,imetabolic)
          emideforest_fuel_10hr(:,ivm,istructural) = cc_litter(:) * bafrac_net(:,ivm) * def_fuel_10hr_remain(:,ivm,istructural)
          emideforest_fuel_100hr(:,ivm,imetabolic) = cc_litter(:) * bafrac_net(:,ivm) * def_fuel_100hr_remain(:,ivm,imetabolic)
          emideforest_fuel_100hr(:,ivm,istructural) = cc_litter(:) * bafrac_net(:,ivm) * def_fuel_100hr_remain(:,ivm,istructural)
          emideforest_fuel_1000hr(:,ivm,imetabolic) = cc_litter(:) * bafrac_net(:,ivm) * def_fuel_1000hr_remain(:,ivm,imetabolic)
          emideforest_fuel_1000hr(:,ivm,istructural) = cc_litter(:) * bafrac_net(:,ivm) * def_fuel_1000hr_remain(:,ivm,istructural)

          emideforest_biomass(:,ivm,icarbres) = cc_litter(:) * bafrac_net(:,ivm) * deforest_biomass_remain(:,ivm,icarbres)
          emideforest_biomass(:,ivm,ifruit) = cc_litter(:) * bafrac_net(:,ivm) * deforest_biomass_remain(:,ivm,ifruit)
          emideforest_biomass(:,ivm,ileaf) = cc_litter(:) * bafrac_net(:,ivm) * deforest_biomass_remain(:,ivm,ileaf)
          emideforest_biomass(:,ivm,iheartabove) = cc_wood(:) * bafrac_net(:,ivm) * deforest_biomass_remain(:,ivm,iheartabove)
          emideforest_biomass(:,ivm,isapabove) = cc_wood(:) * bafrac_net(:,ivm) * deforest_biomass_remain(:,ivm,isapabove)
        ENDWHERE 
      END IF
    END DO
   
   END SUBROUTINE deforestation_fire_proxy

!! ================================================================================================================================
!! SUBROUTINE	: spitfire_month_input
!!
!>\BRIEF        Reads in the maps with monthly time step but being used in spitfire on daily time step.
!!
!!
!! DESCRIPTION  : The information is read in for a single year for all pixels present in the simulation, and
!!                interpolated to the resolution being used for the current run. This subroutine was adapted
!!                from slowproc_Ninput subroutine. For easy understanding, some variable names were kept as
!!                original.
!!
!! RECENT CHANGE(S): 
!!
!! MAIN OUTPUT VARIABLE(S): output_data
!!
!! REFERENCE(S)	: None.
!! 
!! FLOWCHART : None.
!! \n
!_ ================================================================================================================================

  SUBROUTINE spitfire_month_input(nbpt, lalo, neighbours, resolution, contfrac, output_data, &
                                  filename, field_name)

    !
    !! 0. Variable and parameter declaration
    !

    !
    !! 0.1 Input variables
    !
    INTEGER(i_std), INTENT(in)                           :: nbpt           !! Number of points for which the data needs to be interpolated
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in)           :: lalo           !! Vector of latitude and longitudes (beware of the order !)
    INTEGER(i_std), DIMENSION(nbpt,8), INTENT(in)        :: neighbours     !! Vector of neighbours for each grid point 
    ! (1=N, 2=NE, 3=E, 4=SE, 5=S, 6=SW, 7=W, 8=NW)
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in)           :: resolution     !! The size in km of each grid-box in X and Y
    REAL(r_std), DIMENSION(nbpt), INTENT(in)             :: contfrac       !! Fraction of continent in the grid
    CHARACTER(LEN=*),INTENT(in)                          :: filename       !! Input file name
    CHARACTER(LEN=*),INTENT(in)                          :: field_name     !! Name of the field 

    !
    !! 0.2 Modified variables
    !

    !
    !! 0.3 Output variables
    !
    REAL(r_std)                                          :: output_data (nbpt,12)  !! Final output data.

    !! 0.4 Local variables
    !
    CHARACTER(LEN=30)                                    :: callsign
    INTEGER(i_std)                                       :: iml, jml, lml, tml, fid, ib, ip, jp, vid, l, im
    INTEGER(i_std)                                       :: idi, idi_last, nbvmax
    REAL(r_std)                                          :: coslat
    REAL(r_std), DIMENSION(12)                           :: Ninput_val 
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:)          :: mask
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:,:)        :: sub_index
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)             :: lat_rel, lon_rel 
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:)           :: Ninput_map    !! the variable immediately filled by reading the file
    REAL(r_std), ALLOCATABLE, DIMENSION(:)               :: lat_lu, lon_lu
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)             :: sub_area
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:)           :: resol_lu

    INTEGER(i_std)                                       :: nix, njx, iv, i
    !
    LOGICAL                                              :: ok_interpol = .FALSE.      !! optionnal return of aggregate_2d
    !
    INTEGER                                              :: ALLOC_ERR
    LOGICAL                                              :: latitude_exists, longitude_exists !! Test existence of variables in the input files
!_ ================================================================================================================================




    !Config Key   = FM_FILE
    !Config Desc  = Name of file from which the forest management map is to be read
    !Config If    = OK_STOMATE
    !Config Def   = FMmap.nc
    !Config Help  = The name of the file to be opened to read a forest management
    !Config         map (including a layer for every PFT) is given here. 
    !Config Units = [FILE]
    !
    IF((TRIM(filename) .NE. 'NONE') .AND. (TRIM(filename) .NE. 'none')) THEN

       IF (xios_interpolation) THEN
          ! Read and interpolate with XIOS
          ! There should be 12 time step for the concerning variable in the file
          CALL xios_orchidee_recv_field(TRIM(field_name), output_data)

       ELSE
          ! Read with IOIPSL and interpolate with aggregate

          IF (is_root_prc) CALL flininfo(filename, iml, jml, lml, tml, fid)
          CALL bcast(iml)
          CALL bcast(jml)
          CALL bcast(lml)
          CALL bcast(tml)
          
          ALLOCATE(lat_lu(jml), STAT=ALLOC_ERR)
          IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'spitfire_month_input','Problem in allocation of variable lat_lu','','')
          
          ALLOCATE(lon_lu(iml), STAT=ALLOC_ERR)
          IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'spitfire_month_input','Problem in allocation of variable lon_lu','','')
          
          ALLOCATE(Ninput_map(iml,jml,tml), STAT=ALLOC_ERR)
          IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'spitfire_month_input','Problem in allocation of variable Ninput_map','','')
          
          ALLOCATE(resol_lu(iml,jml,2), STAT=ALLOC_ERR)
          IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'spitfire_month_input','Problem in allocation of variable resol_lu','','')
          
          WRITE(numout,*) 'Reading the file'
          
          IF (is_root_prc) THEN
             CALL flinquery_var(fid, 'longitude', longitude_exists)
             IF(longitude_exists)THEN
                CALL flinget(fid, 'longitude', iml, 0, 0, 0, 1, 1, lon_lu)
             ELSE
                CALL flinget(fid, 'lon', iml, 0, 0, 0, 1, 1, lon_lu)
             ENDIF
             CALL flinquery_var(fid, 'latitude', latitude_exists)
             IF(latitude_exists)THEN
                CALL flinget(fid, 'latitude', jml, 0, 0, 0, 1, 1, lat_lu)
             ELSE
                CALL flinget(fid, 'lat', jml, 0, 0, 0, 1, 1, lat_lu)
             ENDIF
             CALL flinget(fid, field_name, iml, jml, 0, tml, 1, tml, Ninput_map)
             !
             CALL flinclo(fid)
          ENDIF
          CALL bcast(lon_lu)
          CALL bcast(lat_lu)
          CALL bcast(Ninput_map)
          
          
          ALLOCATE(lon_rel(iml,jml), STAT=ALLOC_ERR)
          IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'slowproc_slope','Problem in allocation of variable lon_rel','','')
          
          ALLOCATE(lat_rel(iml,jml), STAT=ALLOC_ERR)
          IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'slowproc_slope','Problem in allocation of variable lat_rel','','')
          
          DO ip=1,iml
             lat_rel(ip,:) = lat_lu(:)
          ENDDO
          DO jp=1,jml
             lon_rel(:,jp) = lon_lu(:)
          ENDDO
          !
          !
          ! Mask of permitted variables.
          !
          ALLOCATE(mask(iml,jml), STAT=ALLOC_ERR)
          IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'slowproc_slope','Problem in allocation of variable mask','','')
          
          mask(:,:) = zero
          DO ip=1,iml
             DO jp=1,jml
                IF (ANY(Ninput_map(ip,jp,:) .GE. 0.)) THEN
                   mask(ip,jp) = un
                ENDIF
                !
                ! Resolution in longitude
                !
                coslat = MAX( COS( lat_rel(ip,jp) * pi/180. ), mincos )     
                IF ( ip .EQ. 1 ) THEN
                   resol_lu(ip,jp,1) = ABS( lon_rel(ip+1,jp) - lon_rel(ip,jp) ) * pi/180. * R_Earth * coslat
                ELSEIF ( ip .EQ. iml ) THEN
                   resol_lu(ip,jp,1) = ABS( lon_rel(ip,jp) - lon_rel(ip-1,jp) ) * pi/180. * R_Earth * coslat
                ELSE
                   resol_lu(ip,jp,1) = ABS( lon_rel(ip+1,jp) - lon_rel(ip-1,jp) )/2. * pi/180. * R_Earth * coslat
                ENDIF
                !
                ! Resolution in latitude
                !
                IF ( jp .EQ. 1 ) THEN
                   resol_lu(ip,jp,2) = ABS( lat_rel(ip,jp) - lat_rel(ip,jp+1) ) * pi/180. * R_Earth
                ELSEIF ( jp .EQ. jml ) THEN
                   resol_lu(ip,jp,2) = ABS( lat_rel(ip,jp-1) - lat_rel(ip,jp) ) * pi/180. * R_Earth
                ELSE
                   resol_lu(ip,jp,2) =  ABS( lat_rel(ip,jp-1) - lat_rel(ip,jp+1) )/2. * pi/180. * R_Earth
                ENDIF
                !
             ENDDO
          ENDDO
          !
          !
          ! The number of maximum vegetation map points in the GCM grid is estimated.
          ! Some lmargin is taken.
          !
          IF (is_root_prc) THEN
             nix=INT(MAXVAL(resolution_g(:,1))/MAXVAL(resol_lu(:,:,1)))+2
             njx=INT(MAXVAL(resolution_g(:,2))/MAXVAL(resol_lu(:,:,2)))+2
             nbvmax = nix*njx
          ENDIF
          CALL bcast(nbvmax)
          !
          callsign="Reading monthly maps for spitfire"
          ok_interpol = .FALSE.
          DO WHILE ( .NOT. ok_interpol )
             !
             WRITE(numout,*) "Projection arrays for ",callsign," : "
             WRITE(numout,*) "nbvmax = ",nbvmax
             
             ALLOCATE(sub_index(nbpt,nbvmax,2), STAT=ALLOC_ERR)
             IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'spitfire_month_input','Problem in allocation of variable sub_index','','')
             sub_index(:,:,:)=0
             
             ALLOCATE(sub_area(nbpt,nbvmax), STAT=ALLOC_ERR)
             IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'spitfire_month_input','Problem in allocation of variable sub_area','','')
             sub_area(:,:)=zero
             
             CALL aggregate_p(nbpt, lalo, neighbours, resolution, contfrac, &
                  &                iml, jml, lon_rel, lat_rel, mask, callsign, &
                  &                nbvmax, sub_index, sub_area, ok_interpol)
             
             IF (.NOT. ok_interpol ) THEN
                IF (printlev_loc>=3) WRITE(numout,*) 'nbvmax will be increased from ',nbvmax,' to ', nbvmax*2
                DEALLOCATE(sub_area)
                DEALLOCATE(sub_index)
                nbvmax = nbvmax * 2
             END IF
          END DO
          !
          !
          DO ib = 1, nbpt
             !-
             !- Reinfiltration coefficient due to the slope: Calculation with parameteres maxlope_ro 
             !-
             Ninput_val(:) = zero
             
             ! Initialize last index to the highest possible 
             idi_last=nbvmax
             DO idi=1, nbvmax
                ! Leave the do loop if all sub areas are treated, sub_area <= 0
                IF ( sub_area(ib,idi) <= zero ) THEN
                   ! Set last index to the last one used
                   idi_last=idi-1
                   ! Exit do loop
                   EXIT
                END IF
                
                ip = sub_index(ib,idi,1)
                jp = sub_index(ib,idi,2)
                
                IF(tml == 12) THEN
                   Ninput_val(:) = Ninput_val(:) + Ninput_map(ip,jp,:) * sub_area(ib,idi)
                ELSE
                   Ninput_val(:) = Ninput_val(:) + Ninput_map(ip,jp,1) * sub_area(ib,idi)
                ENDIF
             ENDDO
             
             IF ( idi_last >= 1 ) THEN
                output_data(ib,:) = Ninput_val(:) / SUM(sub_area(ib,1:idi_last)) 
             ELSE
                CALL ipslerr_p(2,'spitfire_month_input', '', '',&
                     &                 'No information for a point') ! Warning error
                output_data(ib,:) = 0.
             ENDIF
          ENDDO
          !
          
          DEALLOCATE(Ninput_map)
          DEALLOCATE(sub_index)
          DEALLOCATE(sub_area)
          DEALLOCATE(mask)
          DEALLOCATE(lon_lu)
          DEALLOCATE(lat_lu)
          DEALLOCATE(lon_rel)
          DEALLOCATE(lat_rel)
          
       END IF ! xios_interpolation

       !
       WRITE(numout,*) 'Interpolation Done in spitfire_month_input for ',TRIM(field_name)
       !
       !
    ELSE
       output_data(:,:)=zero
    ENDIF

  END SUBROUTINE spitfire_month_input

!! ================================================================================================================================
!! SUBROUTINE	: spitfire_annual_input
!!
!>\BRIEF        Reads in the maps with monthly time step but being used in spitfire on daily time step.
!!
!!
!! DESCRIPTION  : The information is read in for a single year for all pixels present in the simulation, and
!!                interpolated to the resolution being used for the current run. This subroutine was adapted
!!                from slowproc_Ninput subroutine. For easy understanding, some variable names were kept as
!!                original.
!!
!! RECENT CHANGE(S): 
!!
!! MAIN OUTPUT VARIABLE(S): output_data
!!
!! REFERENCE(S)	: None.
!! 
!! FLOWCHART : None.
!! \n
!_ ================================================================================================================================

  SUBROUTINE spitfire_annual_input(nbpt, lalo, neighbours, resolution, contfrac, output_data, &
                                   filename, field_name)

    !
    !! 0. Variable and parameter declaration
    !

    !
    !! 0.1 Input variables
    !
    INTEGER(i_std), INTENT(in)                           :: nbpt           !! Number of points for which the data needs to be interpolated
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in)           :: lalo           !! Vector of latitude and longitudes (beware of the order !)
    INTEGER(i_std), DIMENSION(nbpt,8), INTENT(in)        :: neighbours     !! Vector of neighbours for each grid point 
    ! (1=N, 2=NE, 3=E, 4=SE, 5=S, 6=SW, 7=W, 8=NW)
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in)           :: resolution     !! The size in km of each grid-box in X and Y
    REAL(r_std), DIMENSION(nbpt), INTENT(in)             :: contfrac       !! Fraction of continent in the grid
    CHARACTER(LEN=*),INTENT(in)                          :: filename       !! Input file name
    CHARACTER(LEN=*),INTENT(in)                          :: field_name     !! Name of the field 

    !
    !! 0.2 Modified variables
    !

    !
    !! 0.3 Output variables
    !
    REAL(r_std)                                          :: output_data (nbpt)  !! Final output data.

    !! 0.4 Local variables
    !
    CHARACTER(LEN=30)                                    :: callsign
    INTEGER(i_std)                                       :: iml, jml, lml, tml, fid, ib, ip, jp, vid, l, im
    INTEGER(i_std)                                       :: idi, idi_last, nbvmax
    REAL(r_std)                                          :: coslat
    REAL(r_std)                                          :: Ninput_val 
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:)          :: mask
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:,:)        :: sub_index
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)             :: lat_rel, lon_rel 
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:)           :: Ninput_map    !! the variable immediately filled by reading the file
    REAL(r_std), ALLOCATABLE, DIMENSION(:)               :: lat_lu, lon_lu
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)             :: sub_area
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:)           :: resol_lu

    INTEGER(i_std)                                       :: nix, njx, iv, i
    !
    LOGICAL                                              :: ok_interpol = .FALSE.      !! optionnal return of aggregate_2d
    !
    INTEGER                                              :: ALLOC_ERR
    LOGICAL                                              :: latitude_exists, longitude_exists !! Test existence of variables in the input files
!_ ================================================================================================================================




    !Config Key   = FM_FILE
    !Config Desc  = Name of file from which the forest management map is to be read
    !Config If    = OK_STOMATE
    !Config Def   = FMmap.nc
    !Config Help  = The name of the file to be opened to read a forest management
    !Config         map (including a layer for every PFT) is given here. 
    !Config Units = [FILE]
    !
    IF((TRIM(filename) .NE. 'NONE') .AND. (TRIM(filename) .NE. 'none')) THEN

       IF (xios_interpolation) THEN
          ! Read and interpolate with XIOS
          ! There should be 12 time step for the concerning variable in the file
          CALL xios_orchidee_recv_field(TRIM(field_name), output_data)

       ELSE
          ! Read with IOIPSL and interpolate with aggregate

          IF (is_root_prc) CALL flininfo(filename, iml, jml, lml, tml, fid)
          CALL bcast(iml)
          CALL bcast(jml)
          CALL bcast(lml)
          CALL bcast(tml)

          IF (tml .GT. 1) THEN
             CALL ipslerr_p(3,'spitfire_annual_input','The input time length is bigger than 1','','')
          ENDIF
          
          ALLOCATE(lat_lu(jml), STAT=ALLOC_ERR)
          IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'spitfire_annual_input','Problem in allocation of variable lat_lu','','')
          
          ALLOCATE(lon_lu(iml), STAT=ALLOC_ERR)
          IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'spitfire_annual_input','Problem in allocation of variable lon_lu','','')
          
          ALLOCATE(Ninput_map(iml,jml,tml), STAT=ALLOC_ERR)
          IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'spitfire_annual_input','Problem in allocation of variable Ninput_map','','')
          
          ALLOCATE(resol_lu(iml,jml,2), STAT=ALLOC_ERR)
          IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'spitfire_annual_input','Problem in allocation of variable resol_lu','','')
          
          WRITE(numout,*) 'Reading the file'
          
          IF (is_root_prc) THEN
             CALL flinquery_var(fid, 'longitude', longitude_exists)
             IF(longitude_exists)THEN
                CALL flinget(fid, 'longitude', iml, 0, 0, 0, 1, 1, lon_lu)
             ELSE
                CALL flinget(fid, 'lon', iml, 0, 0, 0, 1, 1, lon_lu)
             ENDIF
             CALL flinquery_var(fid, 'latitude', latitude_exists)
             IF(latitude_exists)THEN
                CALL flinget(fid, 'latitude', jml, 0, 0, 0, 1, 1, lat_lu)
             ELSE
                CALL flinget(fid, 'lat', jml, 0, 0, 0, 1, 1, lat_lu)
             ENDIF
             CALL flinget(fid, field_name, iml, jml, 0, tml, 1, tml, Ninput_map)
             !
             CALL flinclo(fid)
          ENDIF
          CALL bcast(lon_lu)
          CALL bcast(lat_lu)
          CALL bcast(Ninput_map)
          
          
          ALLOCATE(lon_rel(iml,jml), STAT=ALLOC_ERR)
          IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'slowproc_slope','Problem in allocation of variable lon_rel','','')
          
          ALLOCATE(lat_rel(iml,jml), STAT=ALLOC_ERR)
          IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'slowproc_slope','Problem in allocation of variable lat_rel','','')
          
          DO ip=1,iml
             lat_rel(ip,:) = lat_lu(:)
          ENDDO
          DO jp=1,jml
             lon_rel(:,jp) = lon_lu(:)
          ENDDO
          !
          !
          ! Mask of permitted variables.
          !
          ALLOCATE(mask(iml,jml), STAT=ALLOC_ERR)
          IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'slowproc_slope','Problem in allocation of variable mask','','')
          
          mask(:,:) = zero
          DO ip=1,iml
             DO jp=1,jml
                IF (ANY(Ninput_map(ip,jp,:) .GE. 0.)) THEN
                   mask(ip,jp) = un
                ENDIF
                !
                ! Resolution in longitude
                !
                coslat = MAX( COS( lat_rel(ip,jp) * pi/180. ), mincos )     
                IF ( ip .EQ. 1 ) THEN
                   resol_lu(ip,jp,1) = ABS( lon_rel(ip+1,jp) - lon_rel(ip,jp) ) * pi/180. * R_Earth * coslat
                ELSEIF ( ip .EQ. iml ) THEN
                   resol_lu(ip,jp,1) = ABS( lon_rel(ip,jp) - lon_rel(ip-1,jp) ) * pi/180. * R_Earth * coslat
                ELSE
                   resol_lu(ip,jp,1) = ABS( lon_rel(ip+1,jp) - lon_rel(ip-1,jp) )/2. * pi/180. * R_Earth * coslat
                ENDIF
                !
                ! Resolution in latitude
                !
                IF ( jp .EQ. 1 ) THEN
                   resol_lu(ip,jp,2) = ABS( lat_rel(ip,jp) - lat_rel(ip,jp+1) ) * pi/180. * R_Earth
                ELSEIF ( jp .EQ. jml ) THEN
                   resol_lu(ip,jp,2) = ABS( lat_rel(ip,jp-1) - lat_rel(ip,jp) ) * pi/180. * R_Earth
                ELSE
                   resol_lu(ip,jp,2) =  ABS( lat_rel(ip,jp-1) - lat_rel(ip,jp+1) )/2. * pi/180. * R_Earth
                ENDIF
                !
             ENDDO
          ENDDO
          !
          !
          ! The number of maximum vegetation map points in the GCM grid is estimated.
          ! Some lmargin is taken.
          !
          IF (is_root_prc) THEN
             nix=INT(MAXVAL(resolution_g(:,1))/MAXVAL(resol_lu(:,:,1)))+2
             njx=INT(MAXVAL(resolution_g(:,2))/MAXVAL(resol_lu(:,:,2)))+2
             nbvmax = nix*njx
          ENDIF
          CALL bcast(nbvmax)
          !
          callsign="Reading monthly maps for spitfire"
          ok_interpol = .FALSE.
          DO WHILE ( .NOT. ok_interpol )
             !
             WRITE(numout,*) "Projection arrays for ",callsign," : "
             WRITE(numout,*) "nbvmax = ",nbvmax
             
             ALLOCATE(sub_index(nbpt,nbvmax,2), STAT=ALLOC_ERR)
             IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'spitfire_annual_input','Problem in allocation of variable sub_index','','')
             sub_index(:,:,:)=0
             
             ALLOCATE(sub_area(nbpt,nbvmax), STAT=ALLOC_ERR)
             IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'spitfire_annual_input','Problem in allocation of variable sub_area','','')
             sub_area(:,:)=zero
             
             CALL aggregate_p(nbpt, lalo, neighbours, resolution, contfrac, &
                  &                iml, jml, lon_rel, lat_rel, mask, callsign, &
                  &                nbvmax, sub_index, sub_area, ok_interpol)
             
             IF (.NOT. ok_interpol ) THEN
                IF (printlev_loc>=3) WRITE(numout,*) 'nbvmax will be increased from ',nbvmax,' to ', nbvmax*2
                DEALLOCATE(sub_area)
                DEALLOCATE(sub_index)
                nbvmax = nbvmax * 2
             END IF
          END DO
          !
          !
          DO ib = 1, nbpt
             !-
             !- Reinfiltration coefficient due to the slope: Calculation with parameteres maxlope_ro 
             !-
             Ninput_val = zero
             
             ! Initialize last index to the highest possible 
             idi_last=nbvmax
             DO idi=1, nbvmax
                ! Leave the do loop if all sub areas are treated, sub_area <= 0
                IF ( sub_area(ib,idi) <= zero ) THEN
                   ! Set last index to the last one used
                   idi_last=idi-1
                   ! Exit do loop
                   EXIT
                END IF
                
                ip = sub_index(ib,idi,1)
                jp = sub_index(ib,idi,2)
                Ninput_val = Ninput_val + Ninput_map(ip,jp,1) * sub_area(ib,idi)

             ENDDO
             
             IF ( idi_last >= 1 ) THEN
                output_data(ib) = Ninput_val / SUM(sub_area(ib,1:idi_last)) 
             ELSE
                CALL ipslerr_p(2,'spitfire_annual_input', '', '',&
                     &                 'No information for a point') ! Warning error
                output_data(ib) = 0.
             ENDIF
          ENDDO
          !
          
          DEALLOCATE(Ninput_map)
          DEALLOCATE(sub_index)
          DEALLOCATE(sub_area)
          DEALLOCATE(mask)
          DEALLOCATE(lon_lu)
          DEALLOCATE(lat_lu)
          DEALLOCATE(lon_rel)
          DEALLOCATE(lat_rel)
          
       END IF ! xios_interpolation

       ! Output the variables read for control only
       ! CALL xios_orchidee_send_field("interp_diag_"//TRIM(field_name),output_data)

       !
       WRITE(numout,*) 'Interpolation Done in spitfire_annual_input for ',TRIM(field_name)
       !
       !
    ELSE
       output_data(:)=zero
    ENDIF

  END SUBROUTINE spitfire_annual_input

END MODULE stomate_spitfire
