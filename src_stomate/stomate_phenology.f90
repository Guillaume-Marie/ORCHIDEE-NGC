! =================================================================================================================================
! MODULE        : stomate_phenology
!
! CONTACT       : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE       : IPSL (2006). This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        This module manages the beginning of the growing season (leaf onset).
!!      
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_stomate/stomate_phenology.f90 $ 
!! $Date: 2026-05-05 09:46:53 +0200 (mar. 05 mai 2026) $
!! $Revision: 9514 $
!! \n
!_ =================================================================================================================================

MODULE stomate_phenology

  ! modules used:
  USE xios_orchidee
  USE ioipsl_para
  USE stomate_data
  USE constantes
  USE pft_parameters
  USE dynamic_parameters
  USE sapiens_agriculture, ONLY: crop_planting
  USE function_library, ONLY: wood_to_height_eff, calculate_c0_alloc, &
       cc_to_biomass, check_vegetation_area, check_mass_balance, &
       wood_to_ba_eff, lai_to_biomass, biomass_to_lai, get_printlev
  USE time, ONLY : ts_annual_proc, julian_diff

  IMPLICIT NONE

  ! private & public routines

  PRIVATE
  PUBLIC phenology,phenology_clear

  ! first call
  LOGICAL, SAVE                           :: firstcall_all_phenology = .TRUE.
!$OMP THREADPRIVATE(firstcall_all_phenology)
  INTEGER(i_std), SAVE                    :: printlev_loc                        !! Local level of text output for current module
!$OMP THREADPRIVATE(printlev_loc)   

 
  LOGICAL, SAVE                           :: firstcall_siggdd= .TRUE.            !! flag used by siggdd for peat, in case of need
!$OMP THREADPRIVATE(firstcall_siggdd)
  
CONTAINS


!! ================================================================================================================================
!! SUBROUTINE   : phenology_clear
!!
!>\BRIEF          Flags setting   
!!
!! DESCRIPTION  : This subroutine sets flags 
!!                ::firstcall_all_phenology and therefore activates section 1.1 of each 
!!                subroutine which writes messages to the output. \n
!!                This subroutine is called at the beginning of ::stomateLpj_clear in the 
!!                ::stomate_lpj module.
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::firstcall_all_phenology
!! REFERENCE(S)  : None
!!
!! FLOWCHART     : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE phenology_clear
    firstcall_all_phenology=.TRUE.
  END SUBROUTINE phenology_clear


!! ================================================================================================================================
!! SUBROUTINE   : phenology
!!
!>\BRIEF          This subroutine controls the detection of the beginning of the growing season 
!!                (if dormancy has been long enough), leaf onset, given favourable biometeorological 
!!                conditions, and leaf growth and biomass allocation when leaf biomass is low (i.e. 
!!                at the start of the growing season.
!!
!! DESCRIPTION  : This subroutine is called by the module ::stomate_lpj and deals with the beginning of the  
!!                growing season. First it is established whether the beginning of the growing season is
!!                allowed. This occurs if the dormance period has been long enough (i.e. greater 
!!                than a minimum PFT-dependent threshold, specified by ::lowgpp_time), 
!!                AND if the last beginning of the growing season was a sufficiently long time ago 
!!                (i.e. when the growing season length is greater than a minimum threshold, specified
!!                by ::min_growthinit_time, which is defined in this module to be 300 days. \n
!!                The dormancy time-length is represented by the variable 
!!                ::time_lowgpp, which is calculated in ::stomate_season. It is increased by 
!!                the stomate time step when the weekly GPP is lower than a threshold. Otherwise
!!                it is set to zero. \n
!!                ::lowgpp_time is set for each PFT in ::stomate_data from a table of all
!!                PFT values (::lowgpp_time_tab), which is defined in ::stomate_constants. \n
!!                The growing season length is given by ::when_growthinit, which increases
!!                by the stomate time-step at each call to this phenology module, except for when
!!                leaf onset is detected, when it is set to 0. \n
!!                If these two conditions are met, leaf onset occurs if the biometeorological 
!!                conditions are also met. This is determined by the leaf onset models, which are
!!                biome-specific. Each PFT is looped over (ignoring bare soil).
!!                The onset phenology model is selected, (according to the parameter 
!!                ::pheno_model, which is initialised in stomate_data), and called. \n
!!                There are six leaf onset phenology models currently being used by ORCHIDEE. 
!!                These are: 'hum' and 'moi', which are based exclusively on moisture conditions,
!!                'humgdd' and 'moigdd', which are based on both temperature and moisture conditions,
!!                'ncdgdd', which is based on a "chilling" requirement for leaf onset, and 
!!                'ngd', which is based on the number of growing days since the temperature was 
!!                above a certain threshold, to account for the end of soil frost.
!!                Those models which are based mostly on temperature conditions are used for
!!                temperate and boreal biomes, and those which include a moisture condition are used
!!                for tropical biomes. More detail on the biometeorological conditions is provided
!!                in the sections on the individual onset models. \n
!!                The moisture conditions are based on the concept of plant "moisture availability".
!!                This is based on the soil humidity (relative soil moisture), but is moderated by
!!                the root density profile, as per the equation:
!!                \latexonly
!!                \input{phenology_vegstress_eqn1.tex}
!!                \endlatexonly
!!                \n
!!                Although some studies have shown that the length of the photoperiod is important
!!                in determining onset (and senescence) dates, this is not considered in the current
!!                versions of the onset models (Krinner et al., 2005). \n
!!                If conditions are favourable, leaf onset occurs (::plant_status is set to ibudbreak), 
!!                ::when_growthinit is set to 0.0, and the growing season has begun. \n
!!                Following the detection of leaf onset, biomass is allocated from the carbohydrate 
!!                reserves equally to the leaves and roots IF the leaf biomass is lower than a minimum
!!                threshold, which is calculated in this subroutine from the parameter
!!                ::lai_initmin, divided by the specific leaf area (both of which are
!!                PFT-dependent and set in ::stomate_constants). \n
!!                Finally, if biomass is required to be allocated from the carbohydrate reserve 
!!                because the leaf biomass is too low, the leaf age and leaf age distribution is 
!!                re-set. In this case the youngest age class fraction is set to 1 and all other   
!!                leaf age class fractions are set to 0. All leaf ages are set to 0. If there is 
!!                no biomass in the carbohydrate reserve, leaf onset will not occur and the PFT
!!                will disappear from the grid cell (Krinner et al., 2005). \n
!!                This subroutine is called in ::stomate_lpj.
!!
!! RECENT CHANGE(S):
!!
!!                   Early 2019: Nitrogen limitations can change the amount of leaves and roots
!!                               that grow during budburst.
!!
!! MAIN OUTPUT VARIABLE(S): ::biomass, 
!!                        ::when_growthinit,
!!                        ::leaf age distribution
!!                        ::leaf fraction
!!
!! REFERENCE(S) :
!! - Krinner, G., N. Viovy, N. de Noblet-Ducoudre, J. Ogee, J. Polcher, P. 
!! Friedlingstein, P. Ciais, S. Sitch and I.C. Prentice (2005), A dynamic global
!! vegetation model for studies of the coupled atmosphere-biosphere system, Global
!! Biogeochemical Cycles, 19, doi:10.1029/2003GB002199.
!!
!! FLOWCHART    : 
!! \latexonly
!! \includegraphics[scale = 1]{phenology_flowchart.png}
!! \endlatexonly
!! \n
!_ ================================================================================================================================

  SUBROUTINE phenology (npts, dt, PFTpresent, &
       veget_max, gpp_week, resp_maint_week, &
       t2m_longterm, t2m_month, t2m_week, &
       maxvegstress_lastyear, minvegstress_lastyear, &
       vegstress_month, vegstress_week, &
       gdd_m5_dormance, gdd_midwinter, ncd_dormance, ngd_minus5, &
       plant_status, time_hum_min, &
       leaf_frac, leaf_age, &
       when_growthinit, atm_to_bm, circ_class_n, &
       circ_class_biomass, KF, &
       longevity_eff_leaf, longevity_eff_sap, longevity_eff_root, age, &
       everywhere, npp_longterm, lm_lastyearmax, k_latosa_adapt, &
       cn_leaf_min_season, plant_n_uptake_daily, cn_leaf_min_2D, &
       cn_leaf_max_2D, lpft_replant, grow_season_len, doy_start_gs, &
       valid_start_gs, mean_start_gs, soil_n_min, som, &
       deepSOM_a, deepSOM_s, deepSOM_p, zf_soil, gdd_from_growthinit )
  
 !! 0. Variable and parameter declaration

   
    !! 0.1 Input variables
    !
    INTEGER(i_std), INTENT(in)                         :: npts                 !! Domain size - number of grid 
                                                                               !! cells (unitless) 
    REAL(r_std), INTENT(in)                            :: dt                   !! time step (dt_days)
    REAL(r_std), DIMENSION(:,:), INTENT(in)            :: veget_max            !! "maximal" coverage fraction of a 
                                                                               !! PFT (LAI -> infinity) on ground 
                                                                               !! (0-1, unitless)
    REAL(r_std), DIMENSION(:), INTENT(in)              :: t2m_longterm         !! "long term" 2 meter reference temperatures (K) 
    REAL(r_std), DIMENSION(:), INTENT(in)              :: t2m_month            !! "monthly" 2-meter temperatures (K) 
    REAL(r_std), DIMENSION(:), INTENT(in)              :: t2m_week             !! "weekly" 2-meter temperatures (K)
    REAL(r_std), DIMENSION(:,:), INTENT(in)            :: gpp_week             !! Weekly gross primary productivity  
    REAL(r_std), DIMENSION(:,:), INTENT(in)            :: resp_maint_week      !! Weekly maintenance respiration
    REAL(r_std), DIMENSION(:,:), INTENT(in)            :: maxvegstress_lastyear !! last year's maximum moisture 
                                                                               !! availability (0-1, unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(in)            :: minvegstress_lastyear !! last year's minimum moisture 
                                                                               !! availability (0-1, unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(in)            :: vegstress_month   !! "monthly" moisture availability 
                                                                               !! (0-1, unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(in)            :: vegstress_week    !! "weekly" moisture availability 
                                                                               !! (0-1, unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(in)            :: gdd_m5_dormance      !! growing degree days above a 
                                                                               !! threshold of -5 deg C (C) 
    REAL(r_std), DIMENSION(:,:), INTENT(in)            :: ncd_dormance         !! number of chilling days since 
                                                                               !! leaves were lost (days) 
    REAL(r_std), DIMENSION(:,:), INTENT(in)            :: longevity_eff_root   !! Effective root turnover time that accounts
                                                                               !! waterstress (days)
    REAL(r_std), DIMENSION(:,:), INTENT(in)            :: longevity_eff_sap    !! Effective sapwood turnover time that accounts
                                                                               !! waterstress (days)
    REAL(r_std), DIMENSION(:,:), INTENT(in)            :: longevity_eff_leaf   !! Effective leaf turnover time that accounts
                                                                               !! waterstress (days)
    REAL(r_std),DIMENSION(:,:), INTENT(in)             :: cn_leaf_max_2D       !! maximal leaf C/N ratio
    REAL(r_std),DIMENSION(:,:), INTENT(in)             :: cn_leaf_min_2D       !! minimal leaf C/N ratio
    LOGICAL, DIMENSION(:,:), INTENT(in)                :: lpft_replant         !! Set to true if a PFT has been clearcut
                                                                               !! and needs to be replaced by another species
    REAL(r_std),DIMENSION(:,:),INTENT(in)              :: mean_start_gs        !! mean growing season starting day for 
                                                                               !! deciduous PFTs (doy).
    REAL(r_std), DIMENSION(0:ngrnd), INTENT(in)        :: zf_soil

    !! 0.2 Ouput variables 
    REAL(r_std),DIMENSION(:,:),INTENT(out)             :: doy_start_gs         !! growing season starting day of year (DOY) for 
                                                                               !! deciduous PFTs.
    LOGICAL, DIMENSION(:,:),INTENT(out)                :: valid_start_gs       !! Doy was calculated by ORCHIDEE. It is therefore valid 
                                                                               !! and should be used to update mean_start_gs (0-1)
    
   
    !! 0.3 Modified variables

    LOGICAL, DIMENSION(:,:), INTENT(inout)             :: PFTpresent           !! PFT exists (true/false)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: time_hum_min         !! time elapsed since strongest moisture availability (days) 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: gdd_midwinter        !! growing degree days, since midwinter (C) 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: when_growthinit      !! how many days since the 
                                                                               !! beginning of the growing season (days) 
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: atm_to_bm            !! C and N taken up by carbohydrate 
                                                                               !! reserve at the beginning of the 
                                                                               !! growing season @tex ($gC m^{-2} 
                                                                               !! of total ground/day$) @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: ngd_minus5           !! number of growing days above a 
                                                                               !! threshold of -5 deg C (days)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: plant_status         !! Growth and phenological status of the plant
                                                                               !! Different stati are listed in constantes_var 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: everywhere           !! is the PFT everywhere in the grid box or 
                                                                               !! very localized (after its introduction) (?)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: circ_class_n         !! Number of individuals in each circ class
                                                                               !! @tex $(ind m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: KF                   !! Scaling factor to convert sapwood mass
                                                                               !! into leaf mass (m)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: k_latosa_adapt       !! Leaf to sapwood area adapted for long 
                                                                               !! term water stress (m)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: leaf_frac            !! fraction of leaves in leaf age 
                                                                               !! class (0-1, unitless)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: leaf_age             !! leaf age (days)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: age                  !! mean age (years)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: npp_longterm         !! "long term" net primary productivity
                                                                               !! @tex ($gC m^{-2} year^{-1}$) @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: lm_lastyearmax       !! last year's maximum leaf mass for each PFT 
                                                                               !! @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)   :: circ_class_biomass   !! Biomass components of the model tree  
                                                                               !! within a circumference class
                                                                               !! class @tex $(g C ind^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: cn_leaf_min_season   !! Minn leaf nitrogen concentration (C:N) of 
                                                                               !! the growing season
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: plant_n_uptake_daily !! Uptake of soil N by plants (gN m^(-2))
    REAL(r_std),DIMENSION(:,:),INTENT(inout)           :: grow_season_len      !! growing season length in days for deciduous PFTs.
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: soil_n_min           !! mineral nitrogen in the soil (gN/m**2)  
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: som                  !! Carbon pool: active, slow, or passive, 
                                                                               !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)                   :: n_estab_atm          !! N taken from atm at establishment [g N m-2 day-1] 
    REAL(r_std), DIMENSION(npts,nvm)                   :: n_estab_soil         !! N taken from soil at establishment [g N m-2 day-1] 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: gdd_from_growthinit  !! growing degree days, since growthinit for crops     
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: deepSOM_a            !! Soil carbon discretized with depth active (g/m**3) 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: deepSOM_s            !! Soil carbon discretized with depth slow (g/m**3) 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: deepSOM_p            !! Soil carbon discretized with depth passive (g/m**3) 

    !! 0.4 Local variables
    
    REAL(r_std)                                        :: c0_alloc             !! Root to sapwood tradeoff parameter
    LOGICAL(r_std), DIMENSION(npts)                    :: age_reset            !! does the leaf age distribution 
                                                                               !! have to be reset? (true/false)
    INTEGER(i_std)                                     :: ipts, ivm, ipar      !! indices (unitless)
    INTEGER(i_std)                                     :: imbc, iele, icir     !! indices (unitless)
    INTEGER(i_std)                                     :: iage, icarb, igrn    !! indices (unitless)
    REAL(r_std)                                        :: deficit              !! Carbon that needs to be taken from the
                                                                               !! labile and or reserve pools 
                                                                               !! @tex $(gC m^{-2})$ @endtex     
    REAL(r_std)                                        :: LF                   !! Scaling factor to convert sapwood mass
                                                                               !! into root mass (unitless)
    REAL(r_std), DIMENSION(npts,nvm)                   :: histvar              !! controls the history output 
                                                                               !! level - 0: nothing is written; 
                                                                               !! 10: everything is written 
                                                                               !! (0-10, unitless)
    REAL(r_std)                                        :: bm_wanted            !! biomass we would like to have 
                                                                               !! @tex ($gC m^{-2} of ground$) 
                                                                               !! @endtex 
    REAL(r_std)                                        :: bm_use               !! biomass we use (from 
                                                                               !! carbohydrate reserve or from 
                                                                               !! atmosphere) @tex ($gC m^{-2} of 
                                                                               !! ground$) @endtex
    REAL(r_std)                                        :: bm_wanted_n          !! biomass we would like to have 
                                                                               !! @tex ($gN m^{-2} of ground$) 
                                                                               !! @endtex 
    REAL(r_std)                                        :: bm_use_n             !! nitrogen biomass we use (from  
                                                                               !! labile or reserve pool @tex ($gN m^{-2} of  
                                                                               !! ground$) @endtex
    REAL(r_std)                                        :: cn_leaf_use          !! CN ratio used for allocationg biomass (gC/gN)
    REAL(r_std), DIMENSION(ncirc)                      :: Cl_tree              !! Intermediate variable of individual plant leaf biomass
                                                                               !! @tex $(gC. tree^{-1})$ @endtex 
    REAL(r_std), DIMENSION(ncirc)                      :: Cr_tree              !! Intermediate variable of individual plant root biomass
                                                                               !! @tex $(gC. tree^{-1})$ @endtex 
    REAL(r_std), DIMENSION(ncirc)                      :: Cs_tree              !! Individual plant, sapwood compartment
                                                                               !! @tex $(gC tree^{-1})$ @endtex   
    REAL(r_std)                                        :: Cs_grass             !! Individual plant, sapwood compartment
                                                                               !! @tex $(gC ind^{-1})$ @endtex
    REAL(r_std)                                        :: Cl_init              !! Initial leaf carbon required to start
                                                                               !! the growing season
                                                                               !! @tex $(gC tree^{-1})$ @endtex
    REAL(r_std), DIMENSION(ncirc)                      :: Cl_wanted            !! Target stand-level leaf C biomass to 
                                                                               !! sustain allometric balance.
                                                                               !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std)                                        :: Cr_init              !! Initial root carbon required to start
                                                                               !! the growing season
                                                                               !! @tex $(gC tree^{-1})$ @endtex
    REAL(r_std)                                        :: lm_min               !! minimum leaf mass @tex ($gC 
                                                                               !! m^{-2} of ground$) @endtex
    REAL(r_std), DIMENSION(ncirc)                      :: circ_class_ba_eff    !! Intermediate variable needed to calculate 
                                                                               !! effective tree height (m^2)
    REAL(r_std), DIMENSION(ncirc)                      :: height_eff           !! Effective tree height calculated from 
                                                                               !! allometric 
                                                                               !! relationships (m) but total biomass (for
                                                                               !! allocation calculations)
    REAL(r_std), DIMENSION(npts,nvm,nmbcomp,nelements) :: check_intern         !! Contains the components of the internal
                                                                               !! mass balance chech for this routine
                                                                               !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)         :: closure_intern       !! Check closure of internal mass balance
                                                                               !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)         :: pool_start           !! Start and end pool of this routine 
                                                                               !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)         :: pool_end             !! Start and end pool of this routine 
                                                                               !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(ncirc)                      :: temp_share           !! A weighting factor based on the biomass in a 
                                                                               !! circ class @tex - @endtex
    REAL(r_std), DIMENSION(ncirc)                      :: temp_share_bis       !! Temporary variable to store the share
                                                                               !! of biomass of each circumference class
                                                                               !! to the total biomass
    REAL(r_std)                                        :: bm_use_circ          !! biomass we use for a given circ class  
                                                                               !! $) @endtex
    REAL(r_std), DIMENSION(npts,nvm)                   :: veget_max_begin      !! temporary storage of veget_max to check 
                                                                               !! area conservation
    REAL(r_std)                                        :: reduction            !! Reduction factor to account for possible 
                                                                               !! nitrogen 
                                                                               !! limitation on carbon allocation (unitless)
    REAL(r_std), DIMENSION(ncirc)                      :: overspending         !! Although there is enough C, this C cannot 
                                                                               !! be allocated because there is not enough N 
                                                                               !! @tex $(gC tree^{-1})$ @endtex
    REAL(r_std), DIMENSION(nparts,nelements)           :: tmp_bm               !! Stand level biomass 
                                                                               !! tex $(gC.m^{-2})$ @endtex
    REAL(r_std)                                        :: bm_use_correction    !! Ratio between the carbon required for leaf and 
                                                                               !! root phenology and the carbon available in the
                                                                               !! reserve pools (unitless)
    REAL(r_std), DIMENSION(npts,nvm)                   :: lstress_fac          !! Light stress factor, based on total
                                                                               !! transmitted light (unitless, 0-1)
    REAL(r_std)                                        :: sla_est              !! A first estimate of sla in case its calculation is 
                                                                               !! dynamic @tex $(m^2.gC^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)                   :: tmp_xios             !! temporary variable for send xios 
    REAL(r_std), DIMENSION(npts,nvm)                   :: doy_ipresenescence   !! temporary variable for send xios 
    REAL(r_std), DIMENSION(npts,nvm)                   :: ratio_rm_gpp         !! temporary variable for send xios
    REAL(r_std), DIMENSION(npts,nvm)                   :: error_count          !! Count the number of errors in consistency checks.
    REAL(r_std), DIMENSION(npts,nvm)                   :: harvest_time         !! Prescribed harvest time adjusted for the longterm
                                                                               !! temperature at the pixel (days)
     
!_ ================================================================================================================================

    IF (firstcall_all_phenology) THEN
      !! Initialize local printlev
      ! Guillaume M. -- get_printlev calls getin_p, a COLLECTIVE bcast. Called here,
      ! inside the compute loop under a local firstcall, it hangs the run as soon as the
      ! ranks stop reaching it together. Use printlev, read once at startup by every rank.
      printlev_loc=printlev
    END IF
    IF (printlev_loc>=2) WRITE(numout,*) 'Entering phenology'
    
    !! Initialize
    doy_start_gs(:,:) = zero
    doy_ipresenescence(:,:) = zero
    valid_start_gs(:,:) = .FALSE.
    ratio_rm_gpp(:,:) = zero
    n_estab_atm(:,:) = zero
    n_estab_soil(:,:) = zero

    !! 1. Write current values to history file

    ! Current values for ::when_growthinit 
    CALL xios_orchidee_send_field("WHEN_GROWTHINIT",when_growthinit)
    CALL histwrite_p (hist_id_stomate, 'WHEN_GROWTHINIT', itime, &
         when_growthinit, npts*nvm, horipft_index)

    ! Set and write values for ::PFTpresent
    WHERE(PFTpresent)
       histvar=un
    ELSEWHERE
       histvar=zero
    ENDWHERE
    CALL xios_orchidee_send_field("PFTPRESENT",histvar)
    CALL histwrite_p (hist_id_stomate, 'PFTPRESENT', itime, &
         histvar, npts*nvm, horipft_index)

    ! Set and write values for gdd_midwinter
    WHERE(gdd_midwinter.EQ.undef)
       histvar=val_exp
    ELSEWHERE
       histvar=gdd_midwinter
    ENDWHERE
    CALL xios_orchidee_send_field("GDD_MIDWINTER",histvar)
    CALL histwrite_p (hist_id_stomate, 'GDD_MIDWINTER', itime, &
         histvar, npts*nvm, horipft_index)

    ! Set and write values for gdd_m5_dormance
    WHERE(gdd_m5_dormance.EQ.undef)
       histvar=val_exp
    ELSEWHERE
       histvar=gdd_m5_dormance
    ENDWHERE
    CALL xios_orchidee_send_field('GDD_M5_DORMANCE',histvar)
    CALL histwrite_p (hist_id_stomate, 'GDD_M5_DORMANCE', itime, &
         histvar, npts*nvm, horipft_index)

    ! Set and write values for ncd_dormance
    WHERE(ncd_dormance.EQ.undef)
       histvar=val_exp
    ELSEWHERE
       histvar=ncd_dormance
    ENDWHERE
    CALL xios_orchidee_send_field("NCD_DORMANCE",histvar)
    CALL histwrite_p (hist_id_stomate, 'NCD_DORMANCE', itime, &
         histvar, npts*nvm, horipft_index)

    !! Error checking
    error_count(:,:) = zero
    WHERE (plant_status(:,:).EQ.inone .AND. &
         veget_max(:,:).GT.min_stomate)
       error_count(:,:) = un
    END WHERE
    ! Ignore bare soil both in the sum and the loop
    IF (SUM(SUM(error_count(:,2:),2)).GT.zero) THEN
       DO ipts = 1,npts
          DO ivm = 2,nvm
             IF (error_count(ipts,ivm).GT.zero) THEN
                WRITE(numout,*) 'Pixel, PFT, ', ipts, ivm
                WRITE(numout,*) 'veget_max, plant_status, ', &
                     veget_max(ipts,ivm),plant_status(ipts,ivm)
                CALL ipslerr_p(3,'stomate_phenology',&
                     'vegetated pixel should have a plant_status', &
                     'other than inone (= 0)','')
             ENDIF
          END DO
       END DO
    END IF
    !-

    !! 2.1 Initialize check for mass balance closure
    !  The mass balance is calculated at the end of this routine
    !  in section 6
    !  Initial biomass pool
    IF (err_act.GT.1) THEN

       check_intern(:,:,:,:) = zero
       pool_start(:,:,:) = zero

       DO iele = 1,nelements

          ! atm_to_bm has as intent inout, the variable 
          ! accumulates carbon over the course of a day. 
          ! Use the difference between start and the end of 
          ! this routine
          check_intern(:,:,iatm2land,iele) = - un * &
               atm_to_bm(:,:,iele) * veget_max(:,:) * dt

          DO ipar = 1,nparts
             DO icir = 1,ncirc
                ! Initial biomass pool
                pool_start(:,:,iele) = pool_start(:,:,iele) + &
                     (circ_class_biomass(:,:,icir,ipar,iele) * &
                     circ_class_n(:,:,icir) * veget_max(:,:))
             ENDDO
          ENDDO

       ENDDO
       IF (ok_soil_carbon_discretization) THEN
          ! Soil carbon (gC m-3) * (m2 m-2)
          DO igrn = 1,ngrnd
             pool_start(:,:,initrogen) = pool_start(:,:,initrogen) + &
                  (deepSOM_a(:,igrn,:,initrogen) + deepSOM_s(:,igrn,:,initrogen) + &
                  deepSOM_p(:,igrn,:,initrogen)) * &
                  (zf_soil(igrn)-zf_soil(igrn-1)) * veget_max(:,:)
          END DO
       ELSE       
          DO icarb = 1,ncarb
             pool_start(:,:,initrogen) = pool_start(:,:,initrogen) + &
                  som(:,icarb,:,initrogen) * veget_max(:,:)
          ENDDO
       ENDIF

       pool_start(:,:,initrogen) = pool_start(:,:,initrogen) + &
            (plant_n_uptake_daily(:,:,iammonium) + &
            plant_n_uptake_daily(:,:,initrate)) * veget_max(:,:) + & 
            (soil_n_min(:,:,iammonium) + &
            soil_n_min(:,:,initrate)) * veget_max(:,:) 
       
       !! 2.2 Initialize check for surface area conservation
       !  Veget_max is a INTENT(in) variable and can therefore
       !  not be changed during the course of this subroutine
       !  Check it anyway, in case the intent get changed.
       veget_max_begin(:,:) = veget_max(:,:)

    ENDIF ! err_act.GT.1


    !! 2.4 Output messages 
    !  giving the setting of the ::always_init
    !  and ::min_growthinit_time parameters.
    IF (firstcall_all_phenology) THEN

       firstcall_all_phenology = .FALSE.

       ! Debug
       IF (printlev_loc>=3) THEN
          WRITE(numout,*) 'phenology:'
          WRITE(numout,*) '   > take carbon from atmosphere if carbohydrate' // &
               ' reserve too small (::always_init): ', always_init
          WRITE(numout,*) '   > minimum time since last beginning of a growing' // &
               ' season (d) (::min_growthinit_time): ', min_growthinit_time
       ENDIF
       !-

    ENDIF


    !! 2.5 Update labile nitrogen pool from soil nitrogen uptake
    ! Calculation of plant_uptake in nitrogen_dynamics is now borrowed 
    ! from OCN-P, which should allow us to avoid negative N_uptake values.
    DO ipts = 1,npts

       DO ivm = 1,nvm

             IF(SUM(plant_n_uptake_daily(ipts,ivm,:)) .GT. zero) THEN

             IF (SUM( (circ_class_biomass(ipts,ivm,:,ilabile,initrogen) + &
                  circ_class_biomass(ipts,ivm,:,icarbres,initrogen)) * &
                  circ_class_n(ipts,ivm,:)) .GT.zero ) THEN

                ! There are trees
                IF ( is_tree(ivm) ) THEN

                   ! Calculate the shares of ilabile
                   temp_share_bis(:) = ((circ_class_biomass(ipts,ivm,:,ilabile,initrogen) + &
                        circ_class_biomass(ipts,ivm,:,icarbres,initrogen)) * &
                        circ_class_n(ipts,ivm,:))/&
                        ((SUM(circ_class_biomass(ipts,ivm,:,ilabile,initrogen) + &
                        circ_class_biomass(ipts,ivm,:,icarbres,initrogen))* &
                        circ_class_n(ipts,ivm,:)))
                   temp_share_bis(ncirc) = zero
                   temp_share_bis(ncirc) = 1 - SUM(temp_share_bis(:))

                   ! Add the N that was taken up by the plant to its labile pool
                   DO icir = 1,ncirc

                      ! Make the distribution of plant_n_uptake_daily proportional 
                      ! to the nitrogen. plant_n_uptake_daily is in gN m-2 divide by
                      ! number of trees to get it in gN tree-1.
                      circ_class_biomass(ipts,ivm,icir,ilabile,initrogen) = &
                           circ_class_biomass(ipts,ivm,icir,ilabile,initrogen) + &
                           (plant_n_uptake_daily(ipts,ivm,iammonium) + &
                           plant_n_uptake_daily(ipts,ivm,initrate)) * &
                           temp_share_bis(icir) / circ_class_n(ipts,ivm,icir)

                   ENDDO
                   
                ELSE

                   ! Grass or cropland
                   circ_class_biomass(ipts,ivm,1,ilabile,initrogen) = &
                        circ_class_biomass(ipts,ivm,1,ilabile,initrogen) + &
                        (plant_n_uptake_daily(ipts,ivm,iammonium) + &
                        plant_n_uptake_daily(ipts,ivm,initrate)) / &
                        circ_class_n(ipts,ivm,1)

                ENDIF ! tree, grass or crop

                plant_n_uptake_daily(ipts,ivm,iammonium) = zero
                plant_n_uptake_daily(ipts,ivm,initrate) = zero                

             ELSE

                WRITE(numout,*) 'ERROR: strange. There is n_uptake but there are &
                     & no plants', plant_n_uptake_daily(ipts,ivm,initrate), &
                     plant_n_uptake_daily(ipts,ivm,iammonium), &
                     SUM(circ_class_biomass(ipts,ivm,:,ilabile,initrogen) * &
                     circ_class_n(ipts,ivm,:)),ipts,ivm

                CALL ipslerr_p(3,'stomate_phenology',&
                     'Supossedly there is n_uptake but no plants',&
                     'This should never happen - fix the problem','')

             ENDIF ! There is plant biomass

          ELSEIF ((plant_n_uptake_daily(ipts,ivm,iammonium) .LT. zero) .OR. &
               (plant_n_uptake_daily(ipts,ivm,initrate) .LT. zero)) THEN

             ! The problem needs to be solved where it 
             ! is caused. That is to say nitrogen_dynamics in 
             ! stomate_soilcarbon.f90
             WRITE(numout,*) 'ERROR: plant_n_uptake_daily has negative values',&
                  plant_n_uptake_daily(ipts,ivm,initrate),&
                  plant_n_uptake_daily(ipts,ivm,iammonium)
             CALL ipslerr_p (3,'phenology', &
                  'Negative N_uptake values.','','')

          ELSEIF ((plant_n_uptake_daily(ipts,ivm,iammonium) .EQ. zero) .OR. &
               (plant_n_uptake_daily(ipts,ivm,initrate) .EQ. zero)) THEN

             ! Nothing should be done. This test (.EQ.zero) makes sense
             ! because n_uptake may have been set to zero in nitrogen_dynamics

          ELSE

             WRITE(numout,*) 'ERROR: overlooked a condition ',ipts,ivm,&
                  plant_n_uptake_daily(ipts,ivm,initrate),&
                  plant_n_uptake_daily(ipts,ivm,iammonium)
             CALL ipslerr_p(3,'stomate_phenology',&
                  'unexpected case in IF-statement',&
                  'for n_uptake','')

          ENDIF

       ENDDO ! ivm

    ENDDO ! ipts


    !! 3. Detection of the beginning of the growing season 

    !! 3.1 allow detection of the beginning of the growing season 
    !  If dormance was long enough (i.e. when ::time_lowgpp, which 
    !  is calculated in ::stomate_season, is above a certain 
    !  PFT-dependent threshold, ::lowgpp_time, which is given in 
    !  ::stomate_constants), AND the last beginning of growing season 
    !  was a sufficiently long time ago (i.e. when ::when_growthinit, 
    !  which is calculated in this module, is greater than 
    !  ::min_growthinit_time, which is declared at the beginning of 
    !  this module). If these conditions are met, allow_initpheno is 
    !  set to TRUE. Each PFT is looped over.
    DO ivm = 2,nvm

       ! Debug
       IF(printlev_loc>=4)THEN
          DO ipts=1,npts
             IF(ipts == test_grid .AND. ivm == test_pft)THEN
                WRITE(numout,*) "Checking for senescence "
                WRITE(numout,*) "when_growthinit: ",&
                     when_growthinit(test_grid,test_pft)
                WRITE(numout,*) "min_growthinit_time: ",&
                     min_growthinit_time
                WRITE(numout,*)'plant_status', &
                     plant_status(test_grid,test_pft)
                WRITE(numout,*)'ivm', ivm
             ENDIF
          ENDDO
       ENDIF
       !-

       WHERE ( (when_growthinit(:,ivm) .GT. min_growthinit_time) .AND. &
            (.NOT.lpft_replant(:,ivm)) .AND. &
            ( (plant_status(:,ivm) .EQ. iprescribe) .OR. & 
            (plant_status(:,ivm) .EQ. idormant) ) )
          
          ! There are buds available. If the conditions further improve
          ! these buds can break and develop into leaves. As long as
          ! lpft_replant is TRUE the model should wait to replant. 
          ! lpft_replant is set to FALSE at the end of the year when
          ! species changes are accounted for.
          plant_status(:,ivm) = ibudsavail
          
       ENDWHERE

       !! 3.2 increase the ::when_growthinit counter
       WHERE (veget_max(:,ivm) .GT. min_stomate)

          !  The when_growthinit counter which gives the number of days since 
          !  the beginning of the growing season. Needed for allocation and 
          !  for the detection of the beginning of the growing season.
          !  CAREFUL: as longs as the unit of dt is day, subsequent if-statements
          !  using when_growthinit are correct. If the unit of dt would be changed
          !  to seconds the IF-statements making use of when_growthinit are no
          !  longer valid.
          when_growthinit(:,ivm) = when_growthinit(:,ivm) + dt

       ELSEWHERE

          ! There is no vegetation. Hence, days since the start of the growing
          ! season is set to zero. Only for the first day of the simulation 
          ! when_growthinit was set to undef to enable the model growing 
          ! evergreen PFTs starting from day 1.
          when_growthinit(:,ivm) = zero

       ENDWHERE

    ENDDO


     !! 4. Leaf onset

    ! The onset phenology model is selected, (according to the 
    ! parameter ::pheno_model, which is initialised in stomate_data), 
    ! and called. Each PFT is looped over (ignoring bare soil). 
    ! If conditions are favourable, plant_status is set to ibudbreak
    ! parameter used in all the differents models of phenology 
    DO ivm = 2,nvm ! Loop over # PFTs

       SELECT CASE ( pheno_model(ivm) )

       CASE ( 'hum' )

          CALL pheno_hum (npts, ivm, PFTpresent, &
               vegstress_month, vegstress_week, &
               maxvegstress_lastyear, minvegstress_lastyear, &
               plant_status, mean_start_gs, doy_start_gs, &
               valid_start_gs, grow_season_len)

       CASE ( 'moi' )

          CALL pheno_moi (npts, ivm, PFTpresent, &
               time_hum_min, &
               vegstress_month, vegstress_week, &
               plant_status, mean_start_gs, doy_start_gs, &
               valid_start_gs, grow_season_len)

       CASE ( 'ncdgdd' )

          CALL pheno_ncdgdd (npts, ivm, PFTpresent, &
               ncd_dormance, gdd_midwinter, &
               t2m_month, t2m_week, plant_status, &
               mean_start_gs, doy_start_gs, &
               valid_start_gs, grow_season_len)

       CASE ( 'ngd' )

          CALL pheno_ngd (npts, ivm, PFTpresent, ncd_dormance, &
               ngd_minus5, t2m_month, t2m_week, plant_status, &
               mean_start_gs, doy_start_gs, &
               valid_start_gs, grow_season_len)

       CASE ( 'humgdd' )

          CALL pheno_humgdd (npts, ivm, PFTpresent, &
               gdd_m5_dormance, maxvegstress_lastyear, &
               minvegstress_lastyear, &
               t2m_longterm, t2m_month, t2m_week, &
               vegstress_week, vegstress_month, &
               plant_status, mean_start_gs, doy_start_gs, &
               valid_start_gs, grow_season_len)

       CASE ( 'moigdd' )

          CALL pheno_moigdd (npts, ivm, PFTpresent, &
               gdd_m5_dormance, time_hum_min, &
               t2m_longterm, t2m_month, t2m_week, &
               vegstress_week, vegstress_month, &
               plant_status, mean_start_gs, doy_start_gs, &
               valid_start_gs, grow_season_len)

       CASE ( 'siggdd' )
          ! originally from mict-peat, with additional modifications applied according to trunk
          CALL pheno_siggdd (npts, ivm, PFTpresent, &
               gdd_m5_dormance, time_hum_min, &
               t2m_longterm, t2m_month, t2m_week, &
               vegstress_week, vegstress_month, &
               plant_status, mean_start_gs, doy_start_gs, &
               valid_start_gs, grow_season_len)

       CASE ( 'none' )

          ! no action

       CASE default

          WRITE(numout,*) 'phenology: don''t know how to treat this PFT.'
          WRITE(numout,*) '  number: (::j)',ivm
          WRITE(numout,*) '  phenology model (::pheno_model(ivm)) : ', pheno_model(ivm)
          CALL ipslerr_p(3,'stomate phenology','Cannot treat this PFT','','')

       END SELECT

    ENDDO

    ! Check whether ipresenescence happens
    WHERE (gpp_week .GT. min_stomate) 
       ratio_rm_gpp(:,:) = resp_maint_week(:,:)/gpp_week(:,:)
    ENDWHERE 

    harvest_time(:,:) = zero

    !DO ipts = 1,npts
    DO ivm = 2,nvm ! Loop over # PFTs

        IF ( .NOT. natural(ivm)) THEN
             harvest_time(:,ivm) = harvest_time_a(ivm) + &
             (t2m_longterm(:) - harvest_time_b(ivm)) * &
             ABS(t2m_longterm(:) - harvest_time_b(ivm))
        ENDIF

        DO ipts = 1,npts   

          ! The ratio_rm_gpp variable needs some time to stabilize
          ! to avoid going into presenescence at the start of the growing
          ! season because ratio_rm_gpp has not be stabilized yet, a 60
          ! day buffer was added. 
          IF ( .NOT.(pheno_model(ivm) .EQ. 'none') .AND. natural(ivm) &
               .AND. when_growthinit(ipts,ivm) .GT. 60 &
               .AND. plant_status(ipts,ivm) .EQ. icanopy &
               .AND. ratio_rm_gpp(ipts,ivm) .GT. presenescence_ratio(ivm) ) THEN
          
             plant_status(ipts,ivm) = ipresenescence
             doy_ipresenescence(ipts,ivm) = julian_diff           
        
          ! For crops to enter presenescence during the growing season and
          ! before harvest.

         ELSEIF ( ok_crops_ipresenescence .AND. &
              ( .NOT.(pheno_model(ivm) .EQ. 'none') .AND. (.NOT. natural(ivm)) &
              .AND. when_growthinit(ipts,ivm) .GT. gr_gt(ivm) &
              .AND. ( plant_status(ipts,ivm) .EQ. icanopy ) &
              .AND. ( gdd_from_growthinit(ipts,ivm) .GT. (gdd_senescence(ivm) & 
              - pre_gdd(ivm) ) .OR. ( when_growthinit(ipts,ivm) .GT. &
              ( harvest_time(ipts,ivm) - pre_harv(ivm) ) ) ) ) ) THEN

         !ELSEIF ( ok_crops_ipresenescence .AND. & 
         !     ( .NOT.(pheno_model(ivm) .EQ. 'none') .AND. (.NOT. natural(ivm)) &
         !     .AND. ( gdd_from_growthinit(ipts,ivm) .GT. (gdd_senescence(ivm) & 
         !     - pre_gdd(ivm) ) .OR. ( when_growthinit(ipts,ivm) .GT. &
         !     ( harvest_time(ipts,ivm) - pre_harv(ivm) ) ) ) ) ) THEN
                
              plant_status(ipts,ivm) = ipresenescence
              doy_ipresenescence(ipts,ivm) = julian_diff

          ENDIF
       ENDDO
    ENDDO

    ! Make a map of the pixel where budbreak occurred
    tmp_xios(:,:) = zero
    WHERE (plant_status(:,:) .EQ. ibudbreak)
       tmp_xios(:,:) = un
    END WHERE
    CALL xios_orchidee_send_field("QC_PHENO_EVENT",tmp_xios)

    ! Make a map of the pixel where budbreak was forced
    tmp_xios(:,:) = zero
    WHERE (plant_status(:,:) .EQ. ibudbreak .AND. &
           .NOT. valid_start_gs(:,:) )
       tmp_xios(:,:) = un
    ENDWHERE
    CALL xios_orchidee_send_field("QC_FORCED_PHENO_EVENT",tmp_xios)

    CALL xios_orchidee_send_field("DOY_IPRESE",doy_ipresenescence) 

  !! 5. Leaf growth and biomass allocation when leaf biomass is low.

    !! Deal with extreme cases when the evergreen leaf biomass reaches zero
    !  What we experienced with the trunk version (2.0), is that once there is 
    !  no more lai, plants do not grow anymore later. To my understanding [NVui], 
    !  there is no recruitment process. As it was a behavior we wanted to avoid 
    !  within the frame of the CMIP6 exercise, we decided to implement a patch 
    !  (see below) within the trunk. It is still a flagable option with the 
    !  ALWAYS_INIT parameter which is now PFT-dependent.
    !  I [SL] think this is not needed in CAN as we have two processes that
    !  should avoid this situation: (1) recruitment. If there is no canopy, lots 
    !  of light will reach the forest floor and so recruits will establish and
    !  should take over and (2) irrespective of whether we use recruitment we 
    !  expect no leaves -> no gpp -> labile and carbres get depleted -> 
    !  hit threshold to be killed -> kill -> replant.
    !  Nevertheless, a similar functionality was implemented in the merge. Note
    !  that the order of code was changed to deal with the circumference classes
    !  and to ensure mass balance closure.
    IF ( ts_annual_proc ) THEN 
       DO ivm = 2, nvm
          DO ipts = 1, npts
             IF (((veget_max(ipts,ivm) .GT. min_sechiba) .AND. always_init(ivm)) .AND. ((pheno_model(ivm) == "none") &
                  .AND. (SUM(circ_class_biomass(ipts,ivm,:,ileaf,icarbon)) .LT. min_sechiba)))   THEN
                
                ! Warning for the user. It is a patch after all. The code
                ! below avoids instabilities in the vegetation cover but
                ! is not rooted in ecology/physiology. Hence the classification
                ! as a patch.
                WRITE(numout,*) 'PFT number, ',ivm
                CALL ipslerr_p(2,'stomate_phenology', &
                     'The PFT seems to be dying - using the AR6 patch','','')

                ! Use minimum lai to calculate minimum leaf mass gC m-2
                lm_min = lai_to_biomass(lai_initmin(ivm),ivm)

                ! Calculate LF for allocation
                c0_alloc = calculate_c0_alloc(ivm, longevity_eff_root(ipts,ivm), &
                     longevity_eff_sap(ipts,ivm))
                LF = c0_alloc * KF(ipts,ivm)

                ! Calculate current biomass pools gC tree-1 (for grasses and trees)
                ! for grasses only the first icir will contain a value
                Cs_tree(:) = circ_class_biomass(ipts,ivm,:,isapabove,icarbon) + &
                     circ_class_biomass(ipts,ivm,:,isapbelow,icarbon)
                Cr_tree(:) = circ_class_biomass(ipts,ivm,:,iroot,icarbon)
                Cl_tree(:) = circ_class_biomass(ipts,ivm,:,ileaf,icarbon)

                ! A C/N ratio will be needed.
                cn_leaf_use=cn_leaf_min_season(ipts,ivm)
                
                IF ( is_tree(ivm)) THEN

                   !+++CHECK+++
                   ! Never had a test case for evergreen forest. It looks like
                   ! YOU found a good test case so check whether the code below
                   ! does what it is supposed to do. If so remove the ipslerr
                   ! at the bottom of this block and remove the write statements.
                   WRITE(numout,*) 'ivm, ',ivm 
                   WRITE(numout,*) 'before pheno allocation - ileaf C, ', &
                        circ_class_biomass(ipts,ivm,:,ileaf,icarbon) * veget_max(ipts,ivm)
                   WRITE(numout,*) 'before pheno allocation - iroot C, ', &
                        circ_class_biomass(ipts,ivm,:,iroot,icarbon) * veget_max(ipts,ivm)
                   WRITE(numout,*) 'before pheno allocation - ileaf N, ', &
                        circ_class_biomass(ipts,ivm,:,ileaf,initrogen) * veget_max(ipts,ivm)
                   WRITE(numout,*) 'before pheno allocation - iroot N, ', &
                        circ_class_biomass(ipts,ivm,:,iroot,initrogen) * veget_max(ipts,ivm)
                   
                   ! Basal area at the tree level (m2 tree-1)
                   circ_class_ba_eff(:) = wood_to_ba_eff(circ_class_biomass(ipts,ivm,:,:,icarbon),&
                        ivm,pipe_tune2(ipts,ivm))

                   !  Calculate tree height
                   height_eff(:) = pipe_tune2(ipts,ivm)*(4/pi*circ_class_ba_eff(:))**&
                        (pipe_tune3(ivm)/2)

                   ! There are no leaves left so the leaves that will be added are the only
                   ! leaves. reset their age
                   age_reset(ipts) = .TRUE.

                   ! Due to the IF-statement we know that have enough Cs. We will allocate
                   ! leaves and roots. We are not sure whether we really need to allocate
                   ! to the roots so check Cr_init against Cr_tree
                   ! Calculate the leaves and roots that need to be grown see 
                   ! stomate_growth_fun_alloc.f90 for more details 
                   ! Note that the units for Cl_init and Cr_init are gC m-2.
                   ! Account for the little bit of leaves we still have
                   Cl_init = MIN( lm_min, SUM(KF(ipts,ivm) * Cs_tree(:) / &
                        height_eff(:) * circ_class_n(ipts,ivm,:) ) ) - &
                        SUM(Cl_tree(:)*circ_class_n(ipts,ivm,:))

                   ! error checking
                   IF (Cl_init.LT.zero) THEN
                      WRITE(numout,*) 'ivm, Cl_init, ',ivm, icir, Cl_init
                      CALL ipslerr_p(3,'Trying to allocate a negative leaf mass',&
                           'something is very wrong in phenology','','')
                   ENDIF

                   ! Distribute the biomass over the leaves and roots of each tree. Convert from gC m-2
                   ! to gC tree-1. If Cl_init would be based on the allometric relationships this
                   ! would be straightforard but Cl_init is more likley to be determined by lm_min.
                   ! Allocate proportionaly
                   temp_share(:) = (KF(ipts,ivm) * Cs_tree(:) / &
                        height_eff(:) * circ_class_n(ipts,ivm,:) ) / &
                        SUM(KF(ipts,ivm) * Cs_tree(:) / &
                        height_eff(:) * circ_class_n(ipts,ivm,:) )
                   temp_share(ncirc) = zero
                   temp_share(ncirc) = 1 - SUM(temp_share(:))

                   ! Allocate C and N to the different diameter classes. 
                   ! Calculate atm_to_bm first. If not, circ_class_biomass has 
                   ! changed by the time it is accounted for in atm_to_bm.
                   ! circ_class_biomass is in gC tree-1. atm_to_bm is in 
                   ! gC m-2 day-1.
                   DO icir = 1,ncirc

                      ! Distribute the biomass over the leaves of the different circumference classes
                      atm_to_bm(ipts,ivm,icarbon) = atm_to_bm(ipts,ivm,icarbon) + &
                           MAX(zero, temp_share(icir) * Cl_init) / dt
                      circ_class_biomass(ipts,ivm,icir,ileaf,icarbon) = &
                           circ_class_biomass(ipts,ivm,icir,ileaf,icarbon) + &
                           MAX(zero, temp_share(icir) * Cl_init / circ_class_n(ipts,ivm,icir))

                      ! Calculate how much nitrogen is needed. Account for the available nitrogen
                      atm_to_bm(ipts,ivm,initrogen) = atm_to_bm(ipts,ivm,initrogen) + &
                           MAX(zero, (temp_share(icir) * Cl_init / cn_leaf_use ) - & 
                           circ_class_biomass(ipts,ivm,icir,ileaf,initrogen) * &
                           circ_class_n(ipts,ivm,icir)) / dt
                      circ_class_biomass(ipts,ivm,icir,ileaf,initrogen) = &
                           circ_class_biomass(ipts,ivm,icir,ileaf,initrogen) + &
                           MAX(zero, temp_share(icir) * Cl_init / circ_class_n(ipts,ivm,icir) / cn_leaf_use  - & 
                           circ_class_biomass(ipts,ivm,icir,ileaf,initrogen))
          
                      ! Distribute the biomass over the roots of the different circumference classes (gC tree-1)
                      atm_to_bm(ipts,ivm,icarbon) = atm_to_bm(ipts,ivm,icarbon) + & 
                           MAX(zero, temp_share(icir) * Cl_init / LF ) / dt
                      circ_class_biomass(ipts,ivm,icir,iroot,icarbon) = &
                           circ_class_biomass(ipts,ivm,icir,iroot,icarbon) + &
                           MAX(zero, temp_share(icir) * Cl_init / LF / circ_class_n(ipts,ivm,icir))
                      
                      ! Calculate how much nitrogen is needed. Account for the available nitrogen
                      atm_to_bm(ipts,ivm,initrogen) = atm_to_bm(ipts,ivm,initrogen) + &
                           MAX(zero,(temp_share(icir) * Cl_init / LF * fcn_root(ivm) / cn_leaf_use ) - & 
                           circ_class_biomass(ipts,ivm,icir,iroot,initrogen) * &
                           circ_class_n(ipts,ivm,icir)) / dt  
                      circ_class_biomass(ipts,ivm,icir,iroot,initrogen) = &
                           circ_class_biomass(ipts,ivm,icir,iroot,initrogen) + &
                           MAX(zero, temp_share(icir) * Cl_init / LF / circ_class_n(ipts,ivm,icir) * fcn_root(ivm) / cn_leaf_use - & 
                           circ_class_biomass(ipts,ivm,icir,iroot,initrogen))
                                  
                   ENDDO !ncirc      

                   ! Warning for the user. It is a patch after all. The code
                   ! below avoids instabilities in the vegetation cover but
                   ! is not rooted in ecology/physiology. Hence the classification
                   ! as a patch.
                   WRITE(numout,*) 'PFT number, ',ivm
                   WRITE(numout,*) 'After pheno allocation - ileaf C, ', &
                        circ_class_biomass(ipts,ivm,:,ileaf,icarbon) * veget_max(ipts,ivm)
                   WRITE(numout,*) 'After pheno allocation - iroot C, ', &
                        circ_class_biomass(ipts,ivm,:,iroot,icarbon) * veget_max(ipts,ivm)
                   WRITE(numout,*) 'After pheno allocation - ileaf N, ', &
                        circ_class_biomass(ipts,ivm,:,ileaf,initrogen) * veget_max(ipts,ivm)
                   WRITE(numout,*) 'After pheno allocation - iroot N, ', &
                        circ_class_biomass(ipts,ivm,:,iroot,initrogen) * veget_max(ipts,ivm)
                   CALL ipslerr_p(3,'stomate_phenology', &
                        'The PFT seems to be dying - using the AR6 patch',&
                        'the model was stopped because this patch has not been tested yet',&
                        'check whether the code works and remove this stop')
                   !++++++++++

                ELSE

                   !+++CHECK+++
                   ! Never had a test case for evergreen grasses and crops. 
                   ! It looks like YOU found a good test case so check whether the 
                   ! code below does what it is supposed to do. If so, remove the ipslerr
                   ! at the bottom of this block and remove the write statements.
                   WRITE(numout,*) 'ivm, ',ivm 
                   WRITE(numout,*) 'before pheno allocation - ileaf C, ', &
                        circ_class_biomass(ipts,ivm,1,ileaf,icarbon) * veget_max(ipts,ivm)
                   WRITE(numout,*) 'before pheno allocation - iroot C, ', &
                        circ_class_biomass(ipts,ivm,1,iroot,icarbon) * veget_max(ipts,ivm)
                   WRITE(numout,*) 'before pheno allocation - ileaf N, ', &
                        circ_class_biomass(ipts,ivm,1,ileaf,initrogen) * veget_max(ipts,ivm)
                   WRITE(numout,*) 'before pheno allocation - iroot N, ', &
                        circ_class_biomass(ipts,ivm,1,iroot,initrogen) * veget_max(ipts,ivm)
                   ! Grasses were simulated as evergree PFTS in CN-CAN. The code
                   ! may end up here if a grassland dies.
                   ! Calculate the available biomass in sapwood (gC m-2)
                   Cs_grass = circ_class_biomass(ipts,ivm,1,isapabove,icarbon)*&
                        circ_class_n(ipts,ivm,1)
                   
                   ! Due to the IF-statement we know that have enough Cs. We will allocate
                   ! leaves and roots. We are not sure whether we really need to allocate
                   ! to the roots so check Cr_init against Cr_tree
                   ! Calculate the leaves and roots that need to be grown see 
                   ! stomate_growth_fun_alloc.f90 for more details 
                   ! Note that the units for Cl_init and Cr_init are gC m-2.
                   ! Account for the little bit of leaves we still have
                   Cl_init = MIN( lm_min, (Cs_grass * KF(ipts,ivm))  - &
                        (Cl_tree(1)*circ_class_n(ipts,ivm,1)) )

                   ! error checking
                   IF (Cl_init.LT.zero) THEN
                      WRITE(numout,*) 'ivm, Cl_init, ',ivm, Cl_init
                      CALL ipslerr_p(3,'Trying to allocate a negative leaf mass',&
                           'something is very wrong in phenology','','')
                   ENDIF
                   
                   ! Allocate leaf carbon. Calculate atm_to_bm first. If not,
                   ! circ_class_biomass has changed by the time it is accounted
                   ! for in atm_to_bm. circ_class_biomass is in gC tree-1.
                   ! atm_to_bm is in gC m-2 day-1.
                   atm_to_bm(ipts,ivm,icarbon) = atm_to_bm(ipts,ivm,icarbon) + &
                        MAX(zero, Cl_init) / dt
                   circ_class_biomass(ipts,ivm,1,ileaf,icarbon) = circ_class_biomass(ipts,ivm,1,ileaf,icarbon) + &
                        MAX(zero, Cl_init / circ_class_n(ipts,ivm,1))    

                   ! Allocate leaf nitrogen. Account for the available nitrogen
                   atm_to_bm(ipts,ivm,initrogen) = atm_to_bm(ipts,ivm,initrogen) + &
                        MAX(zero,(Cl_init / cn_leaf_use ) - &
                        circ_class_biomass(ipts,ivm,1,ileaf,initrogen)*circ_class_n(ipts,ivm,1))/ dt
                   circ_class_biomass(ipts,ivm,1,ileaf,initrogen) = &
                        circ_class_biomass(ipts,ivm,1,ileaf,initrogen) + &
                        MAX(zero, (Cl_init / circ_class_n(ipts,ivm,1) / cn_leaf_use)  - & 
                        circ_class_biomass(ipts,ivm,1,ileaf,initrogen))
    
                   ! Allocate root carbon
                   atm_to_bm(ipts,ivm,icarbon) = atm_to_bm(ipts,ivm,icarbon) + & 
                        MAX(zero, Cl_init / LF) / dt
                   circ_class_biomass(ipts,ivm,1,iroot,icarbon) = &
                        circ_class_biomass(ipts,ivm,1,iroot,icarbon) + &
                        MAX(zero,(Cl_init / LF / circ_class_n(ipts,ivm,1)))
                  
                   ! Allocate root nitrogen. Account for the available nitrogen
                   atm_to_bm(ipts,ivm,initrogen) = atm_to_bm(ipts,ivm,initrogen) + &
                        MAX(zero,(Cl_init / LF * fcn_root(ivm) / cn_leaf_use ) - &
                        circ_class_biomass(ipts,ivm,1,iroot,initrogen)*circ_class_n(ipts,ivm,1)) / dt
                   circ_class_biomass(ipts,ivm,1,iroot,initrogen) = &
                        circ_class_biomass(ipts,ivm,1,iroot,initrogen) + &
                        MAX(zero, (Cl_init / LF / circ_class_n(ipts,ivm,1) * fcn_root(ivm) / cn_leaf_use ) - & 
                        circ_class_biomass(ipts,ivm,1,iroot,initrogen))

                   ! Warning for the user. It is a patch after all. The code
                   ! below avoids instabilities in the vegetation cover but
                   ! is not rooted in ecology/physiology. Hence the classification
                   ! as a patch.
                   WRITE(numout,*) 'PFT number, ',ivm
                   WRITE(numout,*) 'After pheno allocation - ileaf C, ', &
                        circ_class_biomass(ipts,ivm,1,ileaf,icarbon) * veget_max(ipts,ivm)
                   WRITE(numout,*) 'After pheno allocation - iroot C, ', &
                        circ_class_biomass(ipts,ivm,1,iroot,icarbon) * veget_max(ipts,ivm)
                   WRITE(numout,*) 'After pheno allocation - ileaf N, ', &
                        circ_class_biomass(ipts,ivm,1,ileaf,initrogen) * veget_max(ipts,ivm)
                   WRITE(numout,*) 'After pheno allocation - iroot N, ', &
                        circ_class_biomass(ipts,ivm,1,iroot,initrogen) * veget_max(ipts,ivm)
                   CALL ipslerr_p(3,'stomate_phenology', &
                        'The PFT seems to be dying - using the AR6 patch',&
                        'the model was stopped because this patch has not been tested yet',&
                        'check whether the code works and remove this stop')
                   !++++++++++

                ENDIF !is_tree
          
             ENDIF ! exceptional conditions
          ENDDO ! ivm
       ENDDO ! ipts
    ENDIF ! ts_annual_proc

    !  Leaves start to grow if biometeorological conditions are 
    !  favourable (plant_statuse = ibudbreak)
    DO ivm = 2,nvm ! Loop over # PFTs

       age_reset(:) = .FALSE.
       bm_use = zero
       bm_wanted = zero
 
       DO ipts = 1,npts
          
          ! Check whether the PFT is ready to grow new leaves
          IF ( (veget_max(ipts,ivm) .GT. min_sechiba) .AND. &
               (plant_status(ipts,ivm) .EQ. ibudbreak) ) THEN

             ! Calculate stand level biomass
             tmp_bm(:,:) = cc_to_biomass(npts,ivm,&
                  circ_class_biomass(ipts,ivm,:,:,:),&
                  circ_class_n(ipts,ivm,:))
          
             ! The equations assume that there is no leaf or root
             ! biomass left. Check whether this assumption is violated
             IF ( tmp_bm(ileaf,initrogen).GT.min_stomate .OR. &
                  tmp_bm(iroot,initrogen).GT.min_stomate .OR. &
                  tmp_bm(ileaf,icarbon).GT.min_stomate .OR. &
                  tmp_bm(iroot,icarbon).GT.min_stomate ) THEN

                WRITE(numout,*) 'There is leaf and/or root carbon that should '//&
                     'not be here, something could have gone wrong in senescence.' //&
                     'You may take a look in age_class_dist where plant_status is' //&
                     'merged as that is a rather subtle piece of code.'
                WRITE(numout,*) 'ipts, ivm, ',ipts, ivm
                WRITE(numout,*) 'equation assumes ileaf and iroot are zero', &
                     tmp_bm(ileaf,initrogen), tmp_bm(iroot,initrogen), &
                     tmp_bm(ileaf,icarbon), tmp_bm(iroot,icarbon)
                CALL ipslerr_p (3,'stomate_phenology',&
                     'There is leaf or root carbon that should not be here',&
                     'No carbon expected when phenology is needed',&
                     'The problem could be caused during senescence or age_class_dist')
             ENDIF

             ! NOTE. Plant_status will be updated in stomate_season
             ! If it would be updated here, the plant status would
             ! move from ibudsavail to icanopy within a single day.
             ! As a consequence the value ibudbreak would never occur
             ! in stomate_season and therefore several phenological
             ! counters would not be activated or reset. Hence we
             ! update the status from ibudsavail to ibudbreak in that
             ! routine, keep ibudbreak for the rest of this day in
             ! stomate_lpj, pass it to stomate for the next day,
             ! where it is passed to stomate_season and updated to 
             ! icanopy the day after ibudbreak.

             ! Initialize
             Cl_init = zero
             Cr_init = zero
             Cs_tree(:) = zero
             Cs_grass = zero
             Cl_wanted(:) = zero
             height_eff(:) = zero

             ! Minimum biomass is calculated 
             ! Use the following equation:
             ! \latexonly
             ! \input{phenology_lm_min_eqn2.tex}
             ! \endlatexonly
             ! \n 
             ! The allocation scheme does not need this parameter
             ! but without it larger seedlings would have to be
             ! described. The model runs with only one set of 
             ! allometric equations so a seedling has very little 
             ! leaves and needles. By setting lm_min the seedling 
             ! gets more leaves than prescribed by the allometric 
             ! relations. The first days/weeks GPP will be higher
             ! c is allocated to the wood and the seedling becomes
             ! bigger. The next time phenology is used, the seedling 
             ! already has enough wood to establish a functional
             ! canopy when the allometric relationships are 
             ! respected. In short lm_min could be removed from the
             ! code but then we would have to introduce dynamic
             ! allometric relationships that account for the 
             ! allometric differences between seedlings and mature 
             ! trees.
             lm_min = lai_to_biomass(lai_initmin(ivm),ivm)

             ! The light stress should be calculated making use of Pgap so 
             ! it accounts for LAI crown dimensions and tree distribution. 
             ! However, this would be computationally expensive so we just 
             ! use a first order estimate based on light attenuation model 
             ! by Lambert-Beer. When LAI is low, a lot of light reaches the 
             ! forest floor and so KF should increase to make use of the 
             ! available light by growing leaves. This is a not an exact 
             ! solution because the actual leaf mass could be less than
             ! lm_min but it is a good guess. A value is needed to calculate
             ! KF which in turn us needed to estimate Cl_init.
             lstress_fac = exp(lai_initmin(ivm))

             ! Calculate the sla for this amount of leaf biomass. Use a trick.
             ! use biomass_to_lai to calculate the lai with a dynamic or a 
             ! static sla calculation. Then divide the lai by the biomass to
             ! obtain the actual value for sla (m2 g-1). This is not an
             ! exact solution either but the sla of a canopy with lm_min as
             ! its biomass will be used to estimate KF (in case we use a
             ! dynamic sla).
             sla_est = biomass_to_lai(lm_min,ivm) / lm_min

             ! We might need the c0_alloc factor, so let's
             ! calculate it now.
             c0_alloc = calculate_c0_alloc(ivm, longevity_eff_root(ipts,ivm), &
                  longevity_eff_sap(ipts,ivm))

             ! Calculate KF at the beginning of the growing season                
             KF(ipts,ivm) = (k_latosa_adapt(ipts,ivm) + &
                  lstress_fac(ipts,ivm) * &
                  (k_latosa_max(ivm)-k_latosa_min(ivm))) / &
                  (sla_est * tree_ff(ivm) * pipe_density(ivm))

             ! The minimum leaf biomass is prescribed by ::lm_min 
             ! which in turn is basically prescribed through 
             ! ::lai_initmin. However, lm_min could exceed the leaf 
             ! mass that is required to respect the allometric 
             ! relationships therefore this leaf biomass should be 
             ! calculated. Calculate the allocation factors see 
             ! stomate_prsecribe.f90 and stomate_prescribe.f90 for 
             ! more details 
             LF = c0_alloc * KF(ipts,ivm)
             IF (cn_leaf_min_season(ipts,ivm).GT.zero) THEN
                cn_leaf_use=cn_leaf_min_season(ipts,ivm)
             ELSE
                ! In case it wasn't possible yet to calculate
                ! cn_leaf_min_season, use the prescribed value
                ! stored in cn_leaf_min_2D
                cn_leaf_use=cn_leaf_min_2D(ipts,ivm)
             END IF
             age_reset(ipts) = .TRUE.

             IF ( is_tree(ivm) ) THEN

                ! Calculate the stand structure (gC tree-1)
                Cs_tree(:) = &
                     ( circ_class_biomass(ipts,ivm,:,isapabove,icarbon) + &
                     circ_class_biomass(ipts,ivm,:,isapbelow,icarbon) )
                height_eff(:) = wood_to_height_eff(&
                     circ_class_biomass(ipts,ivm,:,:,icarbon), ivm, &
                     pipe_tune2(ipts,ivm))

                ! Calculate the leaves and roots that need to be grown. See 
                ! stomate_growth_fun_alloc.f90 for more details. lm_min
                ! is defined in gC m-2 so when using circ_class_biomass 
                ! (g tree-1) we need to account for the the number of trees
                ! Note that the units for Cl_init and Cr_init are gC m-2
                IF(SUM(SUM(circ_class_biomass(ipts,ivm,:,:,icarbon),1),1) .GT. min_stomate) THEN  
                   Cl_init = MAX(lm_min, SUM((KF(ipts,ivm) * Cs_tree(:) / &
                        height_eff(:)) * circ_class_n(ipts,ivm,:)))
                   Cr_init = Cl_init/LF
                   Cl_wanted(:) = Cs_tree(:) * KF(ipts,ivm)/height_eff(:)*circ_class_n(ipts,ivm,:)
                ELSE
                   Cl_init = zero
                   Cr_init = zero
                   Cl_wanted(:) = zero
                ENDIF

             ELSEIF (.NOT. is_tree(ivm) .AND. natural(ivm) ) THEN

                ! Only used when grasses are described as a mixed phenology.
                ! In ORCHIDEE-CAN grasses used to be simulated as an evergreen
                ! biome. Calculate the available biomass in sapwood (gC m-2)
                Cs_grass = circ_class_biomass(ipts,ivm,1,isapabove,icarbon)*&
                     circ_class_n(ipts,ivm,1)

                ! Calculate structure of the crop/grassland i.e. the leaves 
                ! and roots that need to be grown see 
                ! stomate_growth_fun_alloc.f90 for more details 
                ! Note that the units for Cl_init and Cr_init are gC m-2. 
                ! For grasses :: heigh_init(ivm) is the minimal height_eff 
                ! and is therefore used at the start of growing season 
                ! (it could be considered the height_eff of Cs)
                Cl_init = MAX(lm_min, Cs_grass * KF(ipts,ivm))
                Cr_init = Cl_init/LF

                ! This is the leaf mass that the grass wants to be happy,
                ! based on the amount of stemwood present and the partition
                ! factor between the sapwood and leaves.
                Cl_wanted(1) = Cs_grass * KF(ipts,ivm)

                ! Debug
                IF(Cl_wanted(1) == zero)THEN
                   WRITE(numout,*) "About to crash"
                   WRITE(numout,*) "ipts, ivm: ",ipts,ivm
                   WRITE(numout,*) "KF(ipts,ivm): ",KF(ipts,ivm)
                   WRITE(numout,*) "Cs_grass: ",Cs_grass
                   WRITE(numout,*) "circ_class_biomass icarbon: ",&
                        circ_class_biomass(ipts,ivm,1,:,icarbon)
                   WRITE(numout,*) "circ_class_biomass initrogen: ",&
                        circ_class_biomass(ipts,ivm,1,:,initrogen)
                   WRITE(numout,*) "circ_class_n: ",circ_class_n(ipts,ivm,1)
                   WRITE(numout,*) "plant_status, ",plant_status(ipts,ivm)
                   CALL ipslerr_p (3,'stomate_phenology', &
                        'We are getting ready to divide by zero.  It is possible that KF is zero', &
                        'or that Cs_grass is zero. This has been observed to happen when a PFT died last year', &
                        'but did not regrow due to trace carbon left in the labile pool.')
                ENDIF
                !-

             ELSEIF ( .NOT. natural(ivm) ) THEN

                ! Crops get planted the day that plant_status is ibudbreak. 
                ! If not, too much of there reserves are used before the 
                ! start of the growing season

                CALL crop_planting(npts, dt, ipts, ivm, &
                     veget_max, PFTpresent, c0_alloc, when_growthinit, &
                     time_hum_min, everywhere, plant_status, &
                     circ_class_n, KF, leaf_frac, age, &
                     npp_longterm, lm_lastyearmax, &
                     circ_class_biomass, atm_to_bm, k_latosa_adapt, &  
                     soil_n_min, som, deepSOM_a, deepSOM_s, deepSOM_p, &
                     zf_soil, n_estab_atm, n_estab_soil )
                
             ENDIF

             ! Debug
             IF(printlev_loc>=4)THEN
                IF(ipts == test_grid .AND. ivm == test_pft)THEN
                   WRITE(numout,*) "Phenology, getting ready to allocate"//&
                        " new leaves for trees - nonesensical for crops"
                   WRITE(numout,*) "i,j: ",ipts,ivm
                   WRITE(numout,*) "height_eff(:) ",height_eff(:)
                   WRITE(numout,*) "Cs_tree(:) ",Cs_tree(:)
                   WRITE(numout,*) "Cs_grass ",Cs_grass
                   WRITE(numout,*) "Cl_init carbon ",Cl_init
                   WRITE(numout,*) "Cr_init carbon ",Cr_init
                   WRITE(numout,*) "LF ",LF
                   WRITE(numout,*) "KF(ipts,ivm) ",KF(ipts,ivm)
                   WRITE(numout,*) "lm_min ",lm_min
                END IF
             END IF
             !-

             ! Finalize allocation for phenological growth
             IF ( natural(ivm) ) THEN

                ! Calculate C/N ratio
                ! Since, we do not have leaves during budbreak, we tried to
                ! calculate cn_leaf_use with the formulation below.
                ! However, this gave problems for the deciduous trees.  
                ! cn_leaf_use = ((LF*bm_use)/(1+LF) +
                ! (LF*bm_use*fcn_root(ivm)) / &
                !     ((1+LF)*LF)) / (0.9 * (SUM(circ_class_n(ipts,ivm,:) * &
                !     (circ_class_biomass(ipts,ivm,:,ilabile,initrogen) + &
                !     circ_class_biomass(ipts,ivm,:,icarbres,initrogen)))))

                ! Debug
                IF(printlev_loc>=4)THEN
                   IF(ipts == test_grid .AND. ivm == test_pft)THEN
                      WRITE(numout,*) 'cn_leaf_use, ', cn_leaf_use
                      WRITE(numout,*) 'labile N, ', &
                           circ_class_biomass(ipts,ivm,:,ilabile,initrogen)
                      WRITE(numout,*) 'reserve N, ', &
                           circ_class_biomass(ipts,ivm,:,icarbres,initrogen)
                   END IF
                END IF
                !-

                ! Avoid crazy cn ratios
                IF(cn_leaf_use .GT. cn_leaf_max_2D(ipts,ivm)) THEN

                   ! Use a threshold instead of the calculated value
                   cn_leaf_use = cn_leaf_max_2D(ipts,ivm)

                ENDIF

                ! Both trees and grasslands are considered to 
                ! be natural. Because grasses and trees live on during 
                ! dormance they will consume reserves and turnover may 
                ! continue. There is no guarantee that the required 
                ! amounts of C for phenology are present and can be used. 
                ! This is now being checked. Calculate how much biomass 
                ! is wanted/required. Substract leaf and root biomass
                ! to improve the mass balance closure (but the values
                ! should be less than min_stomate - see if-statement above).
                bm_wanted = MAX(Cl_init - &
                     tmp_bm(ileaf,icarbon), zero) + &
                     MAX(Cr_init - tmp_bm(iroot,icarbon), zero)
                bm_wanted_n = MAX(Cl_init/cn_leaf_use - &
                     tmp_bm(ileaf,initrogen), zero) + &
                     MAX(Cr_init*fcn_root(ivm)/cn_leaf_use - &
                     tmp_bm(iroot,initrogen), zero)

                ! Specific setting for forced phenology
                ! If the biomass in the carbohydrate reserves is less than the 
                ! required biomass take the required amount of carbon from the 
                ! atmosphere and put it into the labile pool. This only occurs 
                ! if the parameter ::always_init. Fill-up the reserve pools
                ! such that there is enough carbon and nitrogen to support
                ! phenology. Once the reserve pools are filled, continue with
                ! the phenological code.
                IF (always_init(ivm)) THEN

                   IF (tmp_bm(ilabile,icarbon) + &
                        tmp_bm(icarbres,icarbon) .LT. bm_wanted ) THEN

                      ! Force phenology - close carbon balance
                      atm_to_bm(ipts,ivm,icarbon) = atm_to_bm(ipts,ivm,icarbon) + &
                           (bm_wanted - (tmp_bm(ilabile,icarbon) + &
                           tmp_bm(icarbres,icarbon) )) / dt
                      tmp_bm(ilabile,icarbon) = &
                           tmp_bm(ilabile,icarbon) + &
                           ( bm_wanted - (tmp_bm(ilabile,icarbon) + &
                           tmp_bm(icarbres,icarbon)) )

                   ENDIF

                   IF (tmp_bm(ilabile,initrogen) + &
                        tmp_bm(icarbres,initrogen) .LT. bm_wanted_n ) THEN

                      ! Force phenology - close nitrogen balance
                      atm_to_bm(ipts,ivm,initrogen) = atm_to_bm(ipts,ivm,initrogen) + &
                           (bm_wanted_n - (tmp_bm(ilabile,initrogen) + &
                           tmp_bm(icarbres,initrogen) )) / dt
                      tmp_bm(ilabile,initrogen) = &
                           tmp_bm(ilabile,initrogen) + &
                           (bm_wanted_n - (tmp_bm(ilabile,initrogen) + &
                           tmp_bm(icarbres,initrogen) ))

                   ENDIF

                ENDIF ! always_init

                ! The biomass available to use is set to be the minimum of 
                ! 90% of the biomass of the labile pool (if carbon not taken 
                ! from the atmosphere), and the wanted biomass. 
                bm_use = MIN( 0.9 *(tmp_bm(ilabile,icarbon) + &
                     tmp_bm(icarbres,icarbon)), bm_wanted)

                ! If impose_cn=y we are sure we can fullfil the carbon
                ! demand because irrespective of how much nitrogen we 
                ! have for the moment we can just take whatever nitrogen
                ! is missing from the atmosphere
                IF (impose_cn) THEN  ! orig
                   bm_use_n = (LF*bm_use)/(1+LF)/cn_leaf_use + &
                        (LF*bm_use*fcn_root(ivm))/((1+LF)*LF*cn_leaf_use)

                   ! Check if phenology would exhaust the reserve pools
                   ! for nitrogen. Using over 90% of the reserve pools
                   ! is considered as too demanding as it would
                   ! exhausting the reserve pools (gC m-2).
                   deficit = 0.9 * (SUM(circ_class_n(ipts,ivm,:) * &
                        (circ_class_biomass(ipts,ivm,:,ilabile,initrogen) + &
                        circ_class_biomass(ipts,ivm,:,icarbres,initrogen))) - &
                        bm_use_n)

                   IF (deficit .LT. zero) THEN

                      ! Calculate the shares of ilabile
                      temp_share(:) = (circ_class_biomass(ipts,ivm,:,ilabile,initrogen)*&
                           circ_class_n(ipts,ivm,:)) / &
                           SUM(circ_class_biomass(ipts,ivm,:,ilabile,initrogen)*&
                           circ_class_n(ipts,ivm,:))
                      temp_share(ncirc) = zero
                      temp_share(ncirc) = 1 - SUM(temp_share(:))

                      ! The C/N ratio is imposed. To allocate the carbon
                      ! to new leaves we will have to take the N from
                      ! somewhere else. Increase the existing pool so it 
                      ! will be possible to satistfy the demand (gC m-2). 
                      atm_to_bm(ipts,ivm,initrogen) = atm_to_bm(ipts,ivm,initrogen) - &
                           deficit

                      ! Add the same amount of nitrogen to tmp_bm (g N m-2)
                      ! because that is the variable we will use to calculate
                      ! the allocation of the ilabile and icarbres. There is no need to
                      ! add the deficit to circ_class_biomass(ilabile) because the value
                      ! from tmp_bm(ilabile) will copies into circ_class_biomass(ilabile)
                      ! -after unit convertion- towards the end of this subroutine.
                      tmp_bm(ilabile,initrogen) = tmp_bm(ilabile,initrogen) - &
                           deficit

                   ENDIF

                ENDIF

                ! Debug
                IF(printlev_loc>=4)THEN
                   IF(ipts == test_grid .AND. ivm == test_pft)THEN
                      WRITE(numout,*) 'bm_use based on carbon, ',bm_use
                      WRITE(numout,*) 'bm_use based on nitrogen, ',&
                           0.9 * cn_leaf_use * (1+LF) * &
                           (tmp_bm(ilabile,initrogen) + tmp_bm(icarbres,initrogen)) / &
                           (LF+fcn_root(ivm)) 
                   END IF
                END IF
                !-

                ! We now have a first estimate of the carbon we can use
                ! in support of budburst (gC m-2). Calculate how much carbon we
                ! can allocate to budburst if we use 90% of the labile and
                ! reserve nitrogen. This is the final bm_use: it accounts for
                ! the allometry, the available carbon and the available nitrogen.
                ! If impose_cn=y the carbon- and nitrogen-based estimate of the
                ! carbon we can allocate to budburst will be identical.
                ! => bm_use_n = Cl_int/cn_leaf_use + Cr_init*fcn_root/cn_leaf_use
                ! => Cl_init = (LF*bm_use)/(1+LF)
                ! => Cr_init = Cl_init/LF
                ! Assume bm_use_n = 0.9 * (labile_n + carbres_n)
                ! => substitute and solve for bm_use
                ! => 0.9 * cn_leaf_use * (1+LF) * &
                !      (tmp_bm(ilabile,initrogen) + tmp_bm(icarbres,initrogen)) / &
                !      (LF+fcn_root)
                bm_use = MIN(bm_use, 0.9 * cn_leaf_use * (1+LF) * &
                     (tmp_bm(ilabile,initrogen) + tmp_bm(icarbres,initrogen)) / &
                     (LF+fcn_root(ivm)))
                bm_use_n = (LF*bm_use)/(1+LF)/cn_leaf_use + &
                     (LF*bm_use*fcn_root(ivm))/((1+LF)*LF*cn_leaf_use)

                ! Debug
                IF(printlev_loc>=4)THEN
                   IF(ipts == test_grid .AND. ivm == test_pft)THEN
                      WRITE(numout,*) 'bm_use based on carbon and nitrogen, ', &
                           ipts,ivm,bm_use, bm_use_n
                   ENDIF
                ENDIF
                !-

                ! Consistency checking. Check whether we are not overspending. 
                ! bm_use_correction should be 1 or less. If not, something 
                ! unexpected happened
                IF ((tmp_bm(ilabile,icarbon) + tmp_bm(icarbres,icarbon)) == zero) THEN
                   WRITE(numout,*) 'No labile or reserve carbon left in phenology. ', &
                        tmp_bm(ilabile,icarbon) + tmp_bm(icarbres,icarbon)
                   WRITE(numout,*) 'ipts, ivm: ',ipts,ivm
                   WRITE(numout,*) 'This means the code will crash.'
                   CALL ipslerr_p(3,'Unexpected - not enough carbon',&
                        'Zero labile and reserve carbon','','')
                ENDIF

                bm_use_correction = bm_use / (0.9* &
                     (tmp_bm(ilabile,icarbon) + tmp_bm(icarbres,icarbon)))

                IF (bm_use_correction .GT. un+min_stomate) THEN
                   WRITE(numout,*) 'About to crash', &
                        tmp_bm(ilabile,icarbon) + tmp_bm(icarbres,icarbon)
                   WRITE(numout,*) 'bm_use, bm_use_correction', &
                        ipts,ivm,bm_use,bm_use_correction
                   WRITE(numout,*) 'not enough carbon'
                   CALL ipslerr_p(3,'Unexpected - not enough carbon','','','')
                ENDIF

                bm_use_correction = bm_use_n / (0.9 * &
                     (tmp_bm(ilabile,initrogen) + tmp_bm(icarbres,initrogen)))

                IF (bm_use_correction .GT. un+min_stomate) THEN
                   WRITE(numout,*) 'bm_use, bm_use_correction', &
                        ipts,ivm,bm_use,bm_use_correction
                   WRITE(numout,*) 'not enough nitrogen'
                   CALL ipslerr_p(3,'Unexpected - not enough nitrogen','','','')
                ENDIF

                ! Debug
                IF(printlev_loc>=4)THEN
                   IF(ipts == test_grid .AND. ivm == test_pft)THEN
                      WRITE(numout,*) 'bm_use based on carbon and nitrogen, ', &
                           bm_use, bm_use_n
                   END IF
                END IF
                !-

                ! Calculate how much carbon will be left in the labile and reserve
                ! pool.
                IF (bm_use .LT. tmp_bm(ilabile,icarbon)) THEN

                   ! All carbon needed in support of the budbreak will come
                   ! from the labile pool
                   tmp_bm(ilabile,icarbon) = tmp_bm(ilabile,icarbon) - bm_use

                ELSE

                   ! There is not enough carbon in the labile pool. Use
                   ! 90% of the labile pool, leave 10% and take the rest
                   ! from the crabres
                   tmp_bm(icarbres,icarbon) = tmp_bm(icarbres,icarbon) + &
                        0.9 * tmp_bm(ilabile,icarbon) - bm_use
                   tmp_bm(ilabile,icarbon) = 0.1 * tmp_bm(ilabile,icarbon) 

                   ! Check whether we are not overspending
                   IF (tmp_bm(icarbres,icarbon) .LT. -min_stomate) THEN

                      WRITE(numout,*) 'Error: icarbres became negative'
                      WRITE(numout,*) 'tmp_bm(icarbres,icarbon) ',&
                           tmp_bm(icarbres,icarbon)
                      CALL ipslerr_p(3,'Error icrabes became negative',&
                           'something went wrong when calculating bm_use',&
                           'check stomate_phenology','')

                   ENDIF

                ENDIF

                ! Calculate how much nitrogen will be left in the labile and reserve
                ! pool.
                IF (bm_use_n .LT. tmp_bm(ilabile,initrogen)) THEN

                   ! All carbon needed in support of the budbreak will come
                   ! from the labile pool
                   tmp_bm(ilabile,initrogen) = tmp_bm(ilabile,initrogen) - bm_use_n

                ELSE

                   ! There is not enough carbon in the labile pool. Use 90%
                   ! from the labile pool, leave 10% and take the rest from 
                   ! the icarbres pool
                   tmp_bm(icarbres,initrogen) = tmp_bm(icarbres,initrogen) + &
                        0.9 * tmp_bm(ilabile,initrogen) - bm_use_n
                   tmp_bm(ilabile,initrogen) = 0.1 * tmp_bm(ilabile,initrogen)

                   ! Check whether we are not overspending
                   IF (tmp_bm(icarbres,initrogen) .LT. -min_stomate) THEN

                      WRITE(numout,*) 'Error: icarbres became negative'
                      WRITE(numout,*) 'tmp_bm(icarbres,initrogen) ',&
                           tmp_bm(icarbres,initrogen)
                      CALL ipslerr_p(3,'Error icrabes became negative',&
                           'something went wrong when calculating bm_use_n',&
                           'check stomate_phenology','')

                   ENDIF

                ENDIF
                
                ! Share of each circumference class to the total biomass
                temp_share(:) = (Cl_wanted(:) * circ_class_n(ipts,ivm,:)) / &
                     SUM(Cl_wanted(:) * circ_class_n(ipts,ivm,:))
                temp_share(ncirc) = zero
                temp_share(ncirc) = 1 - SUM(temp_share(:))

                ! Calculate the biomass components for the different circ classes  
                DO icir =1,ncirc

                   ! If there is no biomass in this circ class, we'll skip it
                   ! This should ALWAYS be the case for ncirc > 1 for grasses.
                   IF (Cl_wanted(icir) .EQ. zero) CYCLE

                   ! Biomass that will be allocated to this circumference
                   ! class. Convert biomass to gC tree-1 so the allometric 
                   ! relationships can be applied
                   bm_use_circ = bm_use * temp_share(icir) / &
                        circ_class_n(ipts,ivm,icir) 

                   ! bm_use_circ is g C / tree for this circumference class. 
                   ! These are the exact same units as circ_class_biomass.
                   ! bm_use may differ from Cr_init and Cl_init, the final 
                   ! allocation will depend on bm_use. Distribute bm_use 
                   ! over leaves and roots following allometric 
                   ! relationships. Cl_init and Cr_init take the same units 
                   ! as bm_use i.e. gC tree-1
                   ! (i) bm_use = Cl_incp + Cr_incp
                   ! (ii) Cr_incp = (Cl_incp+Cl)/LF - Cr
                   ! Substitue (ii) in (i) and solve for Cl_inc
                   ! <=> Cl_incp = (LF*(b_incp+Cr)-Cl)/(1+LF)
                   ! Because this is the start of the growing season, 
                   ! Cr = Cl = 0. Allocate carbon to the leaves and roots in 
                   ! the given circumference class 
                   Cl_init = (LF * bm_use_circ)/(un + LF) 
                   Cr_init = Cl_init/LF

                   ! The sum of Cl_init and Cr_init should equal bm_use_circ 
                   IF (bm_use_circ - Cl_init - Cr_init .LT. -min_stomate .OR. &
                        bm_use_circ - Cl_init - Cr_init .GT. min_stomate) THEN

                      WRITE(numout,*) 'over or underspending carbon '//&
                           'in phenology, ', &
                           bm_use_circ - Cl_init - Cr_init
                      CALL ipslerr_p (3,'stomate_phenology',&
                           'over or underspening carbon in phenology','','')

                   ENDIF

                   ! Distribute the biomass over the leaves and roots 
                   ! (gC tree-1) Since this whole loop is already over 
                   ! circ class, Cl_init and Cr_init are exactly what 
                   ! we need.
                   circ_class_biomass(ipts,ivm,icir,ileaf,icarbon) = Cl_init
                   circ_class_biomass(ipts,ivm,icir,ileaf,initrogen) = &
                        Cl_init / cn_leaf_use 
                   circ_class_biomass(ipts,ivm,icir,iroot,icarbon) = Cr_init
                   circ_class_biomass(ipts,ivm,icir,iroot,initrogen) = &
                        Cr_init * fcn_root(ivm) / cn_leaf_use

                   ! Use this opportunity to redistribute both carbon and nitrogen in
                   ! the labile and reserve pool. This will likely result in more
                   ! consistent cn calculations further down this subroutine
                   circ_class_biomass(ipts,ivm,icir,ilabile,:) = &
                        tmp_bm(ilabile,:) * temp_share(icir) / &
                        circ_class_n(ipts,ivm,icir)
                   circ_class_biomass(ipts,ivm,icir,icarbres,:) = &
                        tmp_bm(icarbres,:) * temp_share(icir) / &
                        circ_class_n(ipts,ivm,icir)

                   ! Debug
                   IF(printlev_loc>=4)THEN
                      IF(ipts == test_grid .AND. ivm == test_pft)THEN
                         WRITE(numout,*) "circ_class_biomass(ipts,ivm,icir,ileaf,:):", &
                              circ_class_biomass(ipts,ivm,icir,ileaf,:)
                         WRITE(numout,*) "circ_class_biomass(ipts,ivm,icir,iroot,:):",&
                              circ_class_biomass(ipts,ivm,icir,iroot,:)
                         WRITE(numout,*) "circ_class_biomass(ipts,ivm,icir,ilabile,:):",&
                              circ_class_biomass(ipts,ivm,icir,ilabile,:)
                         WRITE(numout,*) "circ_class_biomass(ipts,ivm,icir,icarbres,:):",&
                              circ_class_biomass(ipts,ivm,icir,icarbres,:)
                      END IF
                   END IF
                   !-

                ENDDO ! ncirc

                ! Set reset leaf age distribution (::age_reset) flag. 
                ! Default is TRUE.
                ! done later for better vectorization)
                age_reset(ipts) = .TRUE.

                ! Reset when_growthinit counter: start of the growing 
                ! season
                when_growthinit(ipts,ivm) = zero

                ! This was commented out, as it did not give constant CN
                ! ratios for impose_cn yes, and the carbon and nitrogen 
                ! leaf biomass moved in opposite directions 
                ! cn_leaf_min_season(ipts,ivm) = cn_leaf_max_2D(ipts,ivm)

             ELSEIF (.NOT. natural(ivm)) THEN

                ! This are crops and most that needs to
                ! be done was already done in crop_planting
                age_reset(ipts) = .TRUE.

                ! Reset when_growthinit counter: start of the growing 
                ! season
                when_growthinit(ipts,ivm) = zero

             ELSE

                ! Undefined case
                WRITE(numout,*) 'Error: stomate_phenology '//&
                     'ibudbreak is true but case is undefined'
                CALL ipslerr_p (3,'stomate_phenology', &
                     'ibudbreak is true but case is undefined','','')

             END IF ! Is tree(ivm) is natural(ivm)
             
          ENDIF ! plant_status is ibudbreak

       ENDDO ! loop over grid points

       ! Reset leaf age distribution where necessary (i.e. when 
       ! age_reset is TRUE) simply say that everything is in the 
       ! youngest age class. Set the youngest age class fraction 
       ! to 1 and all other leaf age class fractions to 0.
       WHERE ( age_reset(:) )

          leaf_frac(:,ivm,1) = un

       ENDWHERE

       DO iage = 2, nleafages

          WHERE ( age_reset(:) )

             leaf_frac(:,ivm,iage) = zero

          ENDWHERE

       ENDDO ! nleafages

       ! Ages - set all leaf ages to 0.
       DO iage = 1, nleafages

          WHERE ( age_reset(:) )

             leaf_age(:,ivm,iage) = zero

          ENDWHERE

       ENDDO ! nleafages

    ENDDO ! loop over # PFTs


 !! 6. Check numerical consistency of this routine

    IF (err_act.GT.1) THEN

       ! 6.2 Check surface area
       CALL check_vegetation_area("stomate_phenology", npts, veget_max_begin, &
            veget_max,'pft')

       ! 6.3 Mass balance closure 
       ! 6.3.1 Calculate final biomass
       pool_end(:,:,:) = zero 
       DO ipar = 1,nparts
          DO iele = 1,nelements
             DO icir = 1,ncirc
                pool_end(:,:,iele) = pool_end(:,:,iele) + &
                     (circ_class_biomass(:,:,icir,ipar,iele) * &
                     circ_class_n(:,:,icir) * veget_max(:,:))
             ENDDO
          ENDDO
       ENDDO
       IF (ok_soil_carbon_discretization) THEN
          ! Soil carbon (gC m-3) * (m2 m-2)
          DO igrn = 1,ngrnd
             pool_end(:,:,initrogen) = pool_end(:,:,initrogen) + &
                  (deepSOM_a(:,igrn,:,initrogen) + deepSOM_s(:,igrn,:,initrogen) + &
                  deepSOM_p(:,igrn,:,initrogen)) * &
                  (zf_soil(igrn)-zf_soil(igrn-1)) * veget_max(:,:)
          END DO
       ELSE
          DO icarb = 1,ncarb
             pool_end(:,:,initrogen) = pool_end(:,:,initrogen) + &
                  som(:,icarb,:,initrogen) * veget_max(:,:)
          ENDDO
       ENDIF
      pool_end(:,:,initrogen) = pool_end(:,:,initrogen) + &
          (soil_n_min(:,:,initrate) + &
          soil_n_min(:,:,iammonium)) * veget_max(:,:)  

       ! 6.3.2 Calculate mass balance       
       ! Common processes for icarbon and initrogen
      
       DO iele=1,nelements
          check_intern(:,:,iatm2land,iele) = check_intern(:,:,iatm2land,iele) + &
               atm_to_bm(:,:,iele) * veget_max(:,:) * dt
          check_intern(:,:,ipoolchange,iele) = -un * (pool_end(:,:,iele) - &
               pool_start(:,:,iele)) 
       ENDDO

       closure_intern(:,:,:) = zero
       DO imbc = 1,nmbcomp
          DO iele=1,nelements
             closure_intern(:,:,iele) = closure_intern(:,:,iele) + &
                  check_intern(:,:,imbc,iele)
          ENDDO
       ENDDO

       ! 6.3.3 Check mass balance closure
       CALL check_mass_balance("stomate_phenology", closure_intern, npts, pool_end, &
            pool_start, veget_max, 'pft')
      
    ENDIF ! err_act.GT.1

    IF (printlev_loc>=3) WRITE(numout,*) 'Leaving phenology'

 !! 7. Output growing season length
    DO ivm = 2,nvm 
       DO ipts = 1,npts
          IF ( (plant_status(ipts,ivm) .EQ. ibudbreak) .OR. &
               (plant_status(ipts,ivm) .EQ. icanopy) .OR. &
               (plant_status(ipts,ivm) .EQ. ipresenescence) .OR. &
               (plant_status(ipts,ivm) .EQ. isenescent) ) THEN
             grow_season_len(ipts,ivm) = grow_season_len(ipts,ivm)+1
          ENDIF
       ENDDO
    ENDDO 

    CALL xios_orchidee_send_field("GROW_LENGTH",grow_season_len)
    CALL xios_orchidee_send_field("DOY_START_GS",doy_start_gs)
    CALL xios_orchidee_send_field("N_ESTAB_ATM",n_estab_atm)
    CALL xios_orchidee_send_field("N_ESTAB_SOIL",n_estab_soil )

    ! Average start of the growing season (DOY)
    ! Note that mean_start_gs will only be updated in stomate_season.f90
    ! and is therefore having a one day lag on doy_start_gs.
    tmp_xios(:,:) = xios_default_val
    WHERE (mean_start_gs(:,:) .GT. zero)
       tmp_xios(:,:) = mean_start_gs(:,:)
    END WHERE
    CALL xios_orchidee_send_field("MEAN_START_GS",tmp_xios)

  END SUBROUTINE phenology


!! ================================================================================================================================
!! SUBROUTINE   : pheno_hum 
!!
!>\BRIEF          The 'hum' onset model initiate leaf onset based exclusively on moisture 
!!                availability criteria. 
!!                Currently no PFTs are assigned to this onset model.
!!
!! DESCRIPTION  : This model is for tropical biomes, where temperatures are high but moisture
!!                might be a limiting factor on growth. It is based on leaf onset model 4a in 
!!                Botta et al. (2000), which adopts the approach of Le Roux (1995). \n
!!                Leaf onset occurs if the monthly moisture availability is still quite
!!                low (i.e. lower than the weekly availability), but the weekly availability is 
!!                higher than the critical threshold ::availability_crit (as it reacts faster), 
!!                which indicates the weekly moisture availability is increasing.
!!                OR if the monthly moisture availability is high enough (i.e. above the 
!!                threshold value ::relsoilmoist_always), leaf onset is initiated if this has not 
!!                already happened. This allows vegetation in arid areas to respond to rapidly
!!                changing soil moisture conditions (Krinner et al., 2005). \n
!!                The critical weekly moisture availability threshold (::availability_crit), is
!!                calculated in this subroutine, and is a function of last year's maximum and
!!                minimum moisture availability and the PFT-dependent parameter
!!                ::hum_frac, which specifies how much of last year's available 
!!                moisture is required for leaf onset, as per the equation:
!!                \latexonly
!!                \input{phenology_moi_availcrit_eqn3.tex}
!!                \endlatexonly
!!                \n
!!                ::hum_frac is set for each PFT in ::stomate_data from a table
!!                which contains all the PFT values (::hum_frac_tab) in ::stomate_constants. \n
!!                Last year's maximum and minimum moisture availability and the monthly and 
!!                weekly moisture availability are  
!!                The ::pheno_hum subroutine is called in the subroutine ::phenology. 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::plant_status - specifies whether leaf growth can start.
!!
!! REFERENCE(S) : 
!! - Botta, A., N. Viovy, P. Ciais, P. Friedlingstein and P. Monfray (2000), 
!! A global prognostic scheme of leaf onset using satellite data,
!! Global Change Biology, 207, 337-347.
!! - Le Roux, X. (1995), Etude et modelisation des echanges d'eau et d'energie
!! sol-vegetation-atmosphere dans une savane humide, PhD Thesis, University
!! Pierre et Marie Curie, Paris, France.
!! - Krinner, G., N. Viovy, N. de Noblet-Ducoudre, J. Ogee, J. Polcher, P. 
!! Friedlingstein, P. Ciais, S. Sitch and I.C. Prentice (2005), A dynamic global
!! vegetation model for studies of the coupled atmosphere-biosphere system, Global
!! Biogeochemical Cycles, 19, doi:10.1029/2003GB002199.
!!
!! FLOWCHART    : 
!! \latexonly
!! \includegraphics[scale = 1]{pheno_hum.png}
!! \endlatexonly
!! \n             
!_ ================================================================================================================================

  SUBROUTINE pheno_hum (npts, ivm, PFTpresent, &
       vegstress_month, vegstress_week, &
       maxvegstress_lastyear, minvegstress_lastyear, &
       plant_status, mean_start_gs, doy_start_gs, &
       valid_start_gs, grow_season_len)

    !
    !! 0. Variable and parameter declarations
    !

    !
    !! 0.1 Input variables
    !
    INTEGER(i_std), INTENT(in)                                        :: npts                       !! Domain size - number of 
                                                                                                    !! grid cells (unitless) 
    INTEGER(i_std), INTENT(in)                                        :: ivm                        !! PFT index (unitless)
    LOGICAL, DIMENSION(:,:), INTENT(in)                               :: PFTpresent                 !! PFT exists (true/false)
    REAL(r_std), DIMENSION(:,:), INTENT(in)                           :: vegstress_month         !! "monthly" moisture 
                                                                                                    !! availability (0-1, unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(in)                           :: vegstress_week          !! "weekly" moisture 
                                                                                                    !! availability (0-1, unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(in)                           :: maxvegstress_lastyear   !! last year's maximum 
                                                                                                    !! moisture availability 
                                                                                                    !! (0-1, unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(in)                           :: minvegstress_lastyear   !! last year's minimum 
                                                                                                    !! moisture availability 
                                                                                                    !! (0-1, unitless)
    REAL(r_std),DIMENSION(:,:),INTENT(in)                             :: mean_start_gs              !! mean growing season starting day for 
                                                                                                    !! deciduous PFTs (doy).

    !
    !! 0.2 Output variables
    !
    REAL(r_std),DIMENSION(:,:),INTENT(out)                            :: doy_start_gs      !! growing season starting day of year (DOY) for 
                                                                                           !! deciduous PFTs.
    LOGICAL, DIMENSION(:,:),INTENT(out)                               :: valid_start_gs    !! Doy was calculated by ORCHIDEE. It is therefore valid 
                                                                                           !! and should be used to update mean_start_gs (0-1)

    !
    !! 0.3 Modified variables
    !
    REAL(r_std),DIMENSION(:,:), INTENT(inout)                         :: plant_status      !! Growth and phenological 
                                                                                           !! status of the plant
                                                                                           !! Different stati are listed 
    REAL(r_std),DIMENSION(:,:),INTENT(inout)                          :: grow_season_len   !! growing season length in days for deciduous PFTs.      
                                                                                           !! in constantes_var

    !
    !! 0.4 Local variables
    !
    REAL(r_std), DIMENSION(npts)                                      :: availability_crit     !! critical weekly moisture 
                                                                                               !! availability (0-1, unitless)
    INTEGER(i_std)                                                    :: ipts                  !! index (unitless)
    REAL(r_std)                                                       :: doy_force_pheno       !! Last day of the year that phenology should happen (doy)

!_ ================================================================================================================================

    IF (printlev>=3) WRITE(numout,*) 'Entering hum'

    !
    !! 1. Initializations
    !

    ! check the critical value ::hum_frac is defined. If not, stop.

    IF ( hum_frac(ivm) .EQ. undef ) THEN

       WRITE(numout,*) 'hum: hum_frac is undefined for PFT (::j)',ivm
       CALL ipslerr_p(3,'stomate phenology','hum_frac is undefined for this PFT','','')

    ENDIF

    !
    !! 2. Check if biometeorological conditions are favourable for leaf growth.
    !! The PFT has to be there and start of growing season must be allowed
    !

    DO ipts = 1, npts

       IF ( PFTpresent(ipts,ivm) .AND. &
            plant_status(ipts,ivm) .EQ. ibudsavail ) THEN

          !! 2.1 Calculate the critical weekly moisture availability: depends linearly on the last year 
          !! minimum and maximum moisture availabilities, and on the parameter ::hum_frac.
          availability_crit(ipts) = minvegstress_lastyear(ipts,ivm) + hum_frac(ivm) * &
               ( maxvegstress_lastyear(ipts,ivm) - minvegstress_lastyear(ipts,ivm) )

          ! Fail safe option for phenology. Calculate the last day of the year 
          ! that bud break should happen. As long as the mean_start_gs = 0
          ! the fail safe option cannot be used.
          doy_force_pheno = zero
          IF (mean_start_gs(ipts,ivm) .GT. zero) THEN
             IF (mean_start_gs(ipts,ivm) + force_pheno(ivm) .GT. 365) THEN
                doy_force_pheno = mean_start_gs(ipts,ivm) + force_pheno(ivm) - 365
             ELSE
                doy_force_pheno = mean_start_gs(ipts,ivm) + force_pheno(ivm)
             END IF
          ENDIF

          !! 2.2 Determine if growing season should start (if so, plant_status = ibudbreak).
          !!     Leaf onset occurs if the monthly moisture availability is still quite
          !!     low (i.e. lower than the weekly availability), but the weekly availability is 
          !!     already higher than the critical threshold ::availability_crit (as it reacts faster), 
          !!     which indicates the weekly moisture availability is increasing.
          !!     OR if the monthly moisture availability is high enough (i.e. above the threshold value 
          !!     ::relsoilmoist_always), leaf onset is initiated if this has not already happened.

          IF ( ( ( vegstress_week(ipts,ivm)  .GE. availability_crit(ipts) ) .AND. &
               ( vegstress_month(ipts,ivm) .LT. vegstress_week(ipts,ivm) ) ) .OR. &
               ( vegstress_month(ipts,ivm) .GE. relsoilmoist_always(ivm) ) ) THEN

             plant_status(ipts,ivm) = ibudbreak
             doy_start_gs(ipts,ivm) = julian_diff
             grow_season_len(ipts,ivm) = zero
             valid_start_gs(ipts,ivm) = .TRUE.
             
          ELSEIF (ok_force_pheno .AND. &
               doy_force_pheno .EQ. julian_diff) THEN
             
             ! Force phenology but make sure this value is not 
             ! used in the calculation of the mean doy for budbreak
             ! for this PFT
             plant_status(ipts,ivm) = ibudbreak
             doy_start_gs(ipts,ivm) = julian_diff
             grow_season_len(ipts,ivm) = zero
             valid_start_gs(ipts,ivm) = .FALSE.

          ENDIF

       ENDIF        ! PFT there and start of growing season allowed

    ENDDO ! end loop over grid points

    IF (printlev>=4) WRITE(numout,*) 'Leaving hum'

  END SUBROUTINE pheno_hum


!! ================================================================================================================================
!! SUBROUTINE   : pheno_moi
!!
!>\BRIEF          The 'moi' onset model (::pheno_moi) initiates leaf onset based exclusively 
!!                on moisture availability criteria. 
!!                It is very similar to the 'hum' onset model but instead of the weekly moisture 
!!                availability being higher than a constant threshold, the condition is that the 
!!                moisture minimum happened a sufficiently long time ago. 
!!                Currently PFT 3 (Tropical Broad-leaved Raingreen) is assigned to this model.
!!
!! DESCRIPTION  : This model is for tropical biomes, where temperatures are high but moisture
!!                might be a limiting factor on growth. It is based on leaf onset model 4b in 
!!                Botta et al. (2000).
!!                Leaf onset begins if the plant moisture availability minimum was a sufficiently  
!!                time ago, as specified by the PFT-dependent parameter ::hum_min_time 
!!                AND if the "monthly" moisture availability is lower than the "weekly"
!!                availability (indicating that soil moisture is increasing).
!!                OR if the monthly moisture availability is high enough (i.e. above the threshold 
!!                value ::relsoilmoist_always), leaf onset is initiated if this has not already 
!!                happened. \n
!!                ::hum_min_time is set for each PFT in ::stomate_data, and is 
!!                defined in the table ::hum_min_time_tab in ::stomate_constants. \n
!!                ::relsoilmoist_always is defined for both tree and grass in this subroutine 
!!                (set to 1. and 0.6 respectively). \n
!!                The ::pheno_moi subroutine is called in the subroutine ::phenology. 
!!
!! RECENT CHANGE(S): None
!!        
!! MAIN OUTPUT VARIABLE(S): ::plant_status - specifies whether leaf growth can start.
!!
!! REFERENCE(S) : 
!! - Botta, A., N. Viovy, P. Ciais, P. Friedlingstein and P. Monfray (2000), 
!! A global prognostic scheme of leaf onset using satellite data,
!! Global Change Biology, 207, 337-347.
!! - Krinner, G., N. Viovy, N. de Noblet-Ducoudre, J. Ogee, J. Polcher, P. 
!! Friedlingstein, P. Ciais, S. Sitch and I.C. Prentice (2005), A dynamic global
!! vegetation model for studies of the coupled atmosphere-biosphere system, Global
!! Biogeochemical Cycles, 19, doi:10.1029/2003GB002199.
!!
!! FLOWCHART    : 
!! \latexonly
!! \includegraphics[scale = 1]{pheno_moi.png}
!! \endlatexonly
!! \n
!_ ================================================================================================================================

  SUBROUTINE pheno_moi (npts, ivm, PFTpresent, &
       time_hum_min, &
       vegstress_month, vegstress_week, &
       plant_status, mean_start_gs, doy_start_gs, &
       valid_start_gs, grow_season_len)

    !
    !! 0. Variable and parameter declaration
    !

    !
    !! 0.1 Input variables
    !
    INTEGER(i_std), INTENT(in)                          :: npts                !! Domain size - number of grid cells (unitless)
    INTEGER(i_std), INTENT(in)                          :: ivm                 !! PFT index (unitless)
    LOGICAL, DIMENSION(:,:), INTENT(in)                 :: PFTpresent          !! PFT exists (true/false)
    REAL(r_std), DIMENSION(:,:), INTENT(in)             :: time_hum_min        !! time elapsed since strongest moisture 
                                                                               !! availability (days) 
    REAL(r_std), DIMENSION(:,:), INTENT(in)             :: vegstress_month  !! "monthly" moisture availability (0-1, unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(in)             :: vegstress_week   !! "weekly" moisture availability (0-1, unitless)
    REAL(r_std),DIMENSION(:,:),INTENT(in)               :: mean_start_gs       !! mean growing season starting day for 
                                                                               !! deciduous PFTs (doy).
    !
    !! 0.2 Output variables
    !
    REAL(r_std),DIMENSION(:,:), INTENT(inout)           :: plant_status        !! Growth and phenological status of the plant
                                                                               !! Different stati are listed in constantes_var
    REAL(r_std),DIMENSION(:,:),INTENT(out)              :: doy_start_gs        !! growing season starting day of year (DOY) for 
    LOGICAL, DIMENSION(:,:),INTENT(out)                 :: valid_start_gs      !! Doy was calculated by ORCHIDEE. It is therefore valid 
                                                                               !! and should be used to update mean_start_gs (0-1)

    !
    !! 0.3 Modified variables
    !
    REAL(r_std),DIMENSION(:,:),INTENT(inout)            :: grow_season_len     !! growing season length in days for deciduous PFTs.

    !
    !! 0.4 Local variables
    !
    INTEGER(i_std)                                           :: ipts                            !! index (unitless)
    REAL(r_std)                                              :: doy_force_pheno                 !! Last day of the year that phenology should happen (doy)

!_ ================================================================================================================================

    IF (printlev>=3) WRITE(numout,*) 'Entering moi'

    !
    !! 1. Initializations
    !

    ! check the critical value ::hum_min_time is definded. If not, stop

    IF ( hum_min_time(ivm) .EQ. undef ) THEN

       WRITE(numout,*) 'moi: hum_min_time is undefined for PFT (:,ivm) ',ivm
       CALL ipslerr_p(3,'stomate phenology','hum_min_time is undefined for this PFT','','')

    ENDIF

    !
    !! 2. Check if biometeorological conditions are favourable for leaf growth.
    !! The PFT has to be there and start of growing season must be allowed.
    !

    DO ipts = 1, npts

       IF ( PFTpresent(ipts,ivm) .AND. &
            plant_status(ipts,ivm) .EQ. ibudsavail ) THEN
          
          ! Fail safe option for phenology. Calculate the last day of the year 
          ! that bud break should happen. As long as the mean_start_gs = 0
          ! the fail safe option cannot be used.
          doy_force_pheno = zero
          IF (mean_start_gs(ipts,ivm) .GT. zero) THEN
             IF (mean_start_gs(ipts,ivm) + force_pheno(ivm) .GT. 365) THEN
                doy_force_pheno = mean_start_gs(ipts,ivm) + force_pheno(ivm) - 365
             ELSE
                doy_force_pheno = mean_start_gs(ipts,ivm) + force_pheno(ivm)
             END IF
          ENDIF

          !! 2.1 Determine if growing season should start (if so plant_status = ibudbreak).
          !!     The favorable season starts if the moisture minimum (::time_hum_min) was a sufficiently long 
          !!     time ago, i.e. greater than the threshold specified by the parameter ::hum_min_time 
          !!     and if the "monthly" moisture availability is lower than the "weekly"
          !!     availability (indicating that soil moisture is increasing).
          !!     OR if the monthly moisture availability is high enough (i.e. above the threshold value 
          !!     ::relsoilmoist_always), initiate the growing season if this has not happened yet.
          IF  ( ( ( vegstress_week(ipts,ivm) .GT. vegstress_month(ipts,ivm) ) .AND. &
               ( time_hum_min(ipts,ivm) .GT. hum_min_time(ivm) .OR. &
                time_hum_min(ipts,ivm) .EQ. undef )) .OR. &
               ( vegstress_month(ipts,ivm) .GE. relsoilmoist_always(ivm) ) ) THEN

             plant_status(ipts,ivm) = ibudbreak
             doy_start_gs(ipts,ivm) = julian_diff
             grow_season_len(ipts,ivm) = zero
             valid_start_gs(ipts,ivm) = .TRUE.
             
          ELSEIF (ok_force_pheno .AND. &
               doy_force_pheno .EQ. julian_diff) THEN
             
             ! Force phenology but make sure this value is not 
             ! used in the calculation of the mean doy for budbreak
             ! for this PFT
             plant_status(ipts,ivm) = ibudbreak
             doy_start_gs(ipts,ivm) = julian_diff
             grow_season_len(ipts,ivm) = zero
             valid_start_gs(ipts,ivm) = .FALSE.

          ENDIF

       ENDIF        ! PFT there and start of growing season allowed

    ENDDO ! end loop over grid points

    IF (printlev>=4) WRITE(numout,*) 'Leaving moi'

  END SUBROUTINE pheno_moi


!! ================================================================================================================================
!! SUBROUTINE   : pheno_humgdd
!!
!>\BRIEF          The 'humgdd' onset model initiates leaf onset based on mixed conditions of 
!!                temperature and moisture availability criteria. 
!!                Currently no PFTs are assigned to this onset model. 
!!
!! DESCRIPTION  : In this model the Growing Degree Day (GDD) model (Chuine, 2000) is combined 
!!                with the 'hum' onset model (::pheno_hum), which has previously been described,
!!                in order to account for dependence on both temperature and moisture conditions 
!!                in warmer climates. \n. 
!!                The GDD model specifies that daily temperatures above a threshold of -5  
!!                degrees C are summed, minus this threshold, giving the GDD, starting from 
!!                the beginning of the dormancy period (::time_lowgpp>0), i.e. since the leaves 
!!                were lost. \n.
!!                The dormancy time-length is represented by the variable 
!!                ::time_lowgpp, which is calculated in ::stomate_season. It is increased by 
!!                the stomate time step when the weekly GPP is lower than a threshold. Otherwise
!!                it is set to zero. \n
!!                Leaf onset begins when the a PFT-dependent GDD-threshold is reached.
!!                In addition there are temperature and moisture conditions.
!!                The temperature condition specifies that the monthly temperature has to be 
!!                higher than a constant threshold (::t_always) OR
!!                the weekly temperature is higher than the monthly temperature.
!!                There has to be at least some moisture. The moisture condition 
!!                is exactly the same as the 'hum' onset model (::pheno_hum), which has already
!!                been described. \n
!!                The GDD (::gdd_m5_dormance) is calculated in ::stomate_season. GDD is set to 
!!                undef if beginning of the growing season detected, i.e. when there is GPP 
!!                (::time_lowgpp>0).
!!                The parameter ::t_always is defined as 10 degrees C in this subroutine, 
!!                as are the parameters ::moisture_avail_tree and ::moisture_avail_grass 
!!                (set to 1 and 0.6 respectively), which are used in the moisture condition 
!!                (see ::pheno_moi onset model description). \n
!!                The PFT-dependent GDD threshold (::gdd_crit) is calculated as in the onset 
!!                model ::pheno_humgdd, using the equation:
!!                \latexonly
!!                \input{phenology_hummoigdd_gddcrit_eqn4.tex}
!!                \endlatexonly
!!                \n
!!                The three GDDcrit parameters (::gdd(ivm,*)) are set for each PFT in 
!!                ::stomate_data, and three tables defining each of the three critical GDD
!!                parameters for each PFT is given in ::gdd_crit1_tab, ::gdd_crit2_tab and 
!!                ::gdd_crit3_tab in ::stomate_constants. \n
!!                The ::pheno_humgdd subroutine is called in the subroutine ::phenology. 
!!
!! RECENT CHANGES: None
!!                
!! MAIN OUTPUT VARIABLES: plant_status - specifies whether leaf growth can start
!!
!! REFERENCE(S) : 
!! - Botta, A., N. Viovy, P. Ciais, P. Friedlingstein and P. Monfray (2000), 
!! A global prognostic scheme of leaf onset using satellite data,
!! Global Change Biology, 207, 337-347.
!! - Chuine, I (2000), A unified model for the budburst of trees, Journal of 
!! Theoretical Biology, 207, 337-347.
!! - Krinner, G., N. Viovy, N. de Noblet-Ducoudre, J. Ogee, J. Polcher, P. 
!! Friedlingstein, P. Ciais, S. Sitch and I.C. Prentice (2005), A dynamic global
!! vegetation model for studies of the coupled atmosphere-biosphere system, Global
!! Biogeochemical Cycles, 19, doi:10.1029/2003GB002199.
!!
!! FLOWCHART    : 
!! \latexonly
!! \includegraphics[scale = 1]{pheno_humgdd.png}
!! \endlatexonly
!! \n             
!_ ================================================================================================================================

  SUBROUTINE pheno_humgdd (npts, ivm, PFTpresent, gdd, &
       maxvegstress_lastyear, minvegstress_lastyear, &
       t2m_longterm, t2m_month, t2m_week, &
       vegstress_week, vegstress_month, &
       plant_status, mean_start_gs, doy_start_gs, &
       valid_start_gs, grow_season_len)

    !
    !! 0. Variable and parameter declaration
    !

    !
    !! 0.1 Input variables
    !
    INTEGER(i_std), INTENT(in)                          :: npts                    !! Domain size - number of grid cells 
                                                                                   !! (unitless) 
    INTEGER(i_std), INTENT(in)                          :: ivm                     !! PFT index (unitless)
    LOGICAL, DIMENSION(:,:), INTENT(in)                 :: PFTpresent              !! PFT exists (true/false)
    REAL(r_std), DIMENSION(:,:), INTENT(in)             :: gdd                     !! growing degree days, calculated since 
                                                                                   !! leaves have fallen (C) 
    REAL(r_std), DIMENSION(:,:), INTENT(in)             :: maxvegstress_lastyear!! last year's maximum moisture 
                                                                                   !! availability (0-1, unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(in)             :: minvegstress_lastyear!! last year's minimum moisture 
                                                                                   !! availability (0-1, unitless)
    REAL(r_std), DIMENSION(:), INTENT(in)               :: t2m_longterm            !! "long term" 2 meter temperatures (K)
    REAL(r_std), DIMENSION(:), INTENT(in)               :: t2m_month               !! "monthly" 2-meter temperatures (K)
    REAL(r_std), DIMENSION(:), INTENT(in)               :: t2m_week                !! "weekly" 2-meter temperatures (K)
    REAL(r_std), DIMENSION(:,:), INTENT(in)             :: vegstress_week       !! "weekly" moisture availability 
                                                                                   !! (0-1, unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(in)             :: vegstress_month      !! "monthly" moisture availability 
                                                                                   !! (0-1, unitless)
    REAL(r_std),DIMENSION(:,:),INTENT(in)               :: mean_start_gs           !! mean growing season starting day for 
                                                                                   !! deciduous PFTs (doy).

    !
    !! 0.2 Output variables
    !
    REAL(r_std),DIMENSION(:,:), INTENT(inout)           :: plant_status            !! Growth and phenological status of the plant
                                                                                   !! Different stati are listed in constantes_var
    REAL(r_std),DIMENSION(:,:),INTENT(out)              :: doy_start_gs            !! growing season starting day of year (DOY) for 
    LOGICAL, DIMENSION(:,:),INTENT(out)                 :: valid_start_gs          !! Doy was calculated by ORCHIDEE. It is therefore valid 
                                                                                   !! and should be used to update mean_start_gs (0-1)

    !
    !! 0.3 Modified variables
    !
    REAL(r_std),DIMENSION(:,:),INTENT(inout)            :: grow_season_len         !! growing season length in days for deciduous PFTs.

    !
    !! 0.4 Local variables
    !
    REAL(r_std), DIMENSION(npts)                             :: vegstress_crit               !! critical moisture availability 
                                                                                                !! (0-1, unitless)
    REAL(r_std), DIMENSION(npts)                             :: tl                              !! long term temperature (C)
    REAL(r_std), DIMENSION(npts)                             :: gdd_crit                        !! critical GDD (C)
    INTEGER(i_std)                                           :: ipts                            !! index (unitless)
    REAL(r_std)                                              :: doy_force_pheno                 !! Last day of the year that phenology should happen (doy)

!_ ================================================================================================================================

    IF (printlev>=3) WRITE(numout,*) 'Entering humgdd'

    !
    !! 1. Initializations
    !

    ! check the critical values ::gdd and ::pheno_crit_hum_frac are defined.
    !     If not, stop.

    IF ( ANY(pheno_gdd_crit(ivm,:) .EQ. undef) ) THEN

       WRITE(numout,*) 'humgdd: pheno_gdd_crit is undefined for PFT (::j) ',ivm
       CALL ipslerr_p(3,'stomate phenology','pheno_gdd_crit is undefined for this PFT','','')

    ENDIF

    IF ( hum_frac(ivm) .EQ. undef ) THEN

       WRITE(numout,*) 'humgdd: hum_frac is undefined for PFT (::j) ',ivm
       CALL ipslerr_p(3,'stomate phenology','hum_frac is undefined for this PFT','','')

    ENDIF

    !
    !! 2. Check if biometeorological conditions are favourable for leaf growth.
    !!   The PFT has to be there, start of growing season must be allowed, 
    !!   and GDD has to be defined.
    !

    DO ipts = 1, npts

       IF ( PFTpresent(ipts,ivm) .AND. &
            plant_status(ipts,ivm) .EQ. ibudsavail .AND. &
            ( gdd(ipts,ivm) .NE. undef ) ) THEN

          !! 2.1 Calculate the critical weekly moisture availability: depends linearly on the last year 
          !! minimum and maximum moisture availabilities, and on the parameter ::hum_frac.,
          !! (as in the ::pheno_hum model), as per the equation:
          vegstress_crit(ipts) = minvegstress_lastyear(ipts,ivm) + hum_frac(ivm) * &
               ( maxvegstress_lastyear(ipts,ivm) - minvegstress_lastyear(ipts,ivm) )

          !! 2.2 Calculate the critical GDD (::gdd_crit), which is a function of the PFT-dependent 
          !!     critical GDD and the "long term" 2 meter air temperatures.  

          tl(ipts) =  t2m_longterm(ipts) - ZeroCelsius
          gdd_crit(ipts) = pheno_gdd_crit(ivm,1) + tl(ipts)*pheno_gdd_crit(ivm,2) + &
               tl(ipts)*tl(ipts)*pheno_gdd_crit(ivm,3)

          ! Fail safe option for phenology. Calculate the last day of the year 
          ! that bud break should happen. As long as the mean_start_gs = 0
          ! the fail safe option cannot be used.
          doy_force_pheno = zero
          IF (mean_start_gs(ipts,ivm) .GT. zero) THEN
             IF (mean_start_gs(ipts,ivm) + force_pheno(ivm) .GT. 365) THEN
                doy_force_pheno = mean_start_gs(ipts,ivm) + force_pheno(ivm) - 365
             ELSE
                doy_force_pheno = mean_start_gs(ipts,ivm) + force_pheno(ivm)
             END IF
          ENDIF
          
          !! 2.3 Determine if the growing season should start (plant_status = ibudbreak).
          !!     - Has the critical gdd been reached and is the temperature increasing?
          !!     - Is there at least some humidity/moisture availability?
          !!     This occurs if the critical gdd (::gdd_crit) has been reached 
          !!     AND that is temperature increasing, which is true either if the monthly
          !!     temperature being higher than the threshold ::t_always, OR if the weekly
          !!     temperature is higher than the monthly, 
          !!     AND finally that there is sufficient moisture availability, which is 
          !!     the same condition as for the ::pheno_hum onset model.

          IF ( ( gdd(ipts,ivm) .GE. gdd_crit(ipts) ) .AND. &
               ( ( t2m_week(ipts) .GT. t2m_month(ipts) ) .OR. &
               ( t2m_month(ipts) .GT. t_always(ivm) )          ) .AND. &
               ( ( ( vegstress_week(ipts,ivm)  .GE. vegstress_crit(ipts) ) .AND. &
               ( vegstress_month(ipts,ivm) .LT. vegstress_crit(ipts) )        ) .OR. &
               ( vegstress_month(ipts,ivm) .GE. relsoilmoist_always(ivm) ) ) )  THEN

             plant_status(ipts,ivm) = ibudbreak
             doy_start_gs(ipts,ivm) = julian_diff
             grow_season_len(ipts,ivm) = zero
             valid_start_gs(ipts,ivm) = .TRUE.
             
          ELSEIF (ok_force_pheno .AND. &
               doy_force_pheno .EQ. julian_diff) THEN
             
             ! Force phenology but make sure this value is not 
             ! used in the calculation of the mean doy for budbreak
             ! for this PFT
             plant_status(ipts,ivm) = ibudbreak
             doy_start_gs(ipts,ivm) = julian_diff
             grow_season_len(ipts,ivm) = zero
             valid_start_gs(ipts,ivm) = .FALSE.

          ENDIF

       ENDIF        ! PFT there and start of growing season allowed

    ENDDO ! End loop over grid points

    IF (printlev>=4) WRITE(numout,*) 'Leaving humgdd'

  END SUBROUTINE pheno_humgdd


!! ================================================================================================================================
!! SUBROUTINE   : pheno_moigdd
!!
!>\BRIEF          The 'moigdd' onset model initiates leaf onset based on mixed temperature 
!!                and moisture availability criteria.
!!                Currently PFTs 10 - 13 (C3 and C4 grass, and C3 and C4 agriculture) 
!!                are assigned to this model. 
!!
!! DESCRIPTION  : This onset model combines the GDD model (Chuine, 2000), as described for 
!!                the 'humgdd' onset model (::pheno_humgdd), and the 'moi' model, in order 
!!                to account for dependence on both temperature and moisture conditions in
!!                warmer climates. \n
!!                Leaf onset begins when the a PFT-dependent GDD threshold is reached.
!!                In addition there are temperature and moisture conditions.
!!                The temperature condition specifies that the monthly temperature has to be 
!!                higher than a constant threshold (::t_always) OR
!!                the weekly temperature is higher than the monthly temperature.
!!                There has to be at least some moisture. The moisture condition 
!!                is exactly the same as the 'moi' onset model (::pheno_moi), which has
!!                already been described. \n
!!                GDD is set to undef if beginning of the growing season detected.
!!                As in the ::pheno_humgdd model, the parameter ::t_always is defined as 
!!                10 degrees C in this subroutine, as are the parameters ::moisture_avail_tree
!!                and ::moisture_avail_grass (set to 1 and 0.6 respectively), which are used
!!                in the moisture condition (see ::pheno_moi onset model description). \n
!!                The PFT-dependent GDD threshold (::gdd_crit) is calculated as in the onset 
!!                model ::pheno_humgdd, using the equation:
!!                \latexonly
!!                \input{phenology_hummoigdd_gddcrit_eqn4.tex}
!!                \endlatexonly
!!                \n
!!                where i and j are the grid cell and PFT respectively.
!!                The three GDDcrit parameters (::gdd(ivm,*)) are set for each PFT in 
!!                ::stomate_data, and three tables defining each of the three critical GDD
!!                parameters for each PFT is given in ::gdd_crit1_tab, ::gdd_crit2_tab and 
!!                ::gdd_crit3_tab in ::stomate_constants. \n
!!                The ::pheno_moigdd subroutine is called in the subroutine ::phenology. 
!!
!! RECENT CHANGE(S): Added temperature threshold for C4 grass (pheno_moigdd_t_crit), Dan Zhu april 2015
!!                
!! MAIN OUTPUT VARIABLE(S): ::plant_status - specifies whether leaf growth can start
!!
!! REFERENCE(S) : 
!! - Botta, A., N. Viovy, P. Ciais, P. Friedlingstein and P. Monfray (2000), 
!! A global prognostic scheme of leaf onset using satellite data,
!! Global Change Biology, 207, 337-347.
!! - Chuine, I (2000), A unified model for the budburst of trees, Journal of 
!! Theoretical Biology, 207, 337-347.
!! - Krinner, G., N. Viovy, N. de Noblet-Ducoudre, J. Ogee, J. Polcher, P. 
!! Friedlingstein, P. Ciais, S. Sitch and I.C. Prentice (2005), A dynamic global
!! vegetation model for studies of the coupled atmosphere-biosphere system, Global
!! Biogeochemical Cycles, 19, doi:10.1029/2003GB002199.
!! - Still et al., Global distribution of C3 and C4 vegetation: Carbon cycle implications, 
!! 2003, Global Biogeochemmical Cycles, DOI: 10.1029/2001GB001807. 
!!
!! FLOWCHART    : 
!! \latexonly
!! \includegraphics[scale = 1]{pheno_moigdd.png}
!! \endlatexonly
!! \n
!_ ================================================================================================================================

  SUBROUTINE pheno_moigdd (npts, ivm, PFTpresent, &
       gdd_m5_dormance, time_hum_min, &
       t2m_longterm, t2m_month, t2m_week, &
       vegstress_week, vegstress_month, &
       plant_status, mean_start_gs, doy_start_gs, &
       valid_start_gs, grow_season_len)

    !
    !! 0. Variable and parameter declaration
    !

    !
    !! 0.1 Input variables
    !
    INTEGER(i_std), INTENT(in)                          :: npts            !! Domain size - number of grid cells (unitless)
    INTEGER(i_std), INTENT(in)                          :: ivm             !! PFT index (unitless)
    LOGICAL, DIMENSION(:,:), INTENT(in)                 :: PFTpresent      !! PFT exists (true/false) 
    REAL(r_std), DIMENSION(:,:), INTENT(in)             :: gdd_m5_dormance !! growing degree days, calculated since leaves 
                                                                           !! have fallen (C) 
    REAL(r_std), DIMENSION(:,:), INTENT(in)             :: time_hum_min    !! time elapsed since strongest moisture 
                                                                           !! availability (days) 
    REAL(r_std), DIMENSION(:), INTENT(in)               :: t2m_longterm    !! "long term" 2 meter temperatures (K)
    REAL(r_std), DIMENSION(:), INTENT(in)               :: t2m_month       !! "monthly" 2-meter temperatures (K)
    REAL(r_std), DIMENSION(:), INTENT(in)               :: t2m_week        !! "weekly" 2-meter temperatures (K)
    REAL(r_std), DIMENSION(:,:), INTENT(in)             :: vegstress_week   !! "weekly" moisture availability (0-1, unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(in)             :: vegstress_month  !! "monthly" moisture availability (0-1, unitless)
    REAL(r_std),DIMENSION(:,:),INTENT(in)               :: mean_start_gs   !! mean growing season starting day for 
                                                                           !! deciduous PFTs (doy).

    !
    !! 0.2 Output variables
    !
    REAL(r_std),DIMENSION(:,:),INTENT(out)              :: doy_start_gs    !! growing season starting day of year (DOY) for 
                                                                           !! deciduous PFTs.
    LOGICAL, DIMENSION(:,:),INTENT(out)                 :: valid_start_gs  !! Doy was calculated by ORCHIDEE. It is therefore valid 
                                                                           !! and should be used to update mean_start_gs (0-1)

    !
    !! 0.3 Modified variables
    !
    REAL(r_std),DIMENSION(:,:),INTENT(inout)            :: plant_status    !! Growth and phenological status of the plant
                                                                           !! Different stati are listed in constantes_var
    REAL(r_std),DIMENSION(:,:),INTENT(inout)            :: grow_season_len !! growing season length in days for deciduous PFTs.

    !
    !! 0.4 Local variables
    !
    REAL(r_std), DIMENSION(npts)                             :: tl                              !! long term temperature (C)
    REAL(r_std), DIMENSION(npts)                             :: gdd_crit                        !! critical GDD (C)
    INTEGER(i_std)                                           :: ipts                            !! index (unitless)
    REAL(r_std)                                              :: doy_force_pheno                 !! Last day of the year that phenology should happen (doy)
!_ ================================================================================================================================

    IF (printlev>=3) WRITE(numout,*) 'Entering moigdd'

    !
    !! 1. Initializations
    !

    ! check the critical values ::gdd and ::pheno_crit_hum_min_time are defined.
    !     If not, stop.
  
    IF ( ANY(pheno_gdd_crit(ivm,:) .EQ. undef) ) THEN

       WRITE(numout,*) 'moigdd: pheno_gdd_crit is undefined for PFT',ivm
       CALL ipslerr_p(3,'stomate phenology','pheno_gdd is undefined for this PFT','','')

    ENDIF

    IF ( hum_min_time(ivm) .EQ. undef ) THEN

       WRITE(numout,*) 'moigdd: hum_min_time is undefined for PFT',ivm
       CALL ipslerr_p(3,'stomate phenology','hum_min is undefined for this PFT','','')

    ENDIF

    !
    !! 2. Check if biometeorological conditions are favourable for leaf growth.
    !!    The PFT has to be there, the start of growing season must be allowed, 
    !!    and GDD has to be defined.
    !

    DO ipts = 1, npts

       ! Debug
       IF (printlev_loc>=4)THEN
          IF(ipts == test_grid .AND. ivm == test_pft)THEN
             WRITE(numout,*) 'ivm, ',test_pft
             WRITE(numout,*) 'PFTpresent, ',PFTpresent(test_grid,test_pft)
             WRITE(numout,*) 'plant_status, ',plant_status(test_grid,test_pft)
             WRITE(numout,*) 'ipts, ivm, ',ipts, test_pft
             WRITE(numout,*) 'gdd_m5_dormance, gdd_crit', &
                  gdd_m5_dormance(ipts,test_pft), gdd_crit(ipts)
             WRITE(numout,*) 't2m_week, t2months', &
                  t2m_week(ipts),t2m_month(ipts)
             WRITE(numout,*) 't2m_month, t_always', &
                  t2m_month(ipts),t_always(ivm)
             WRITE(numout,*) 'temperature test 1', &
                  (gdd_m5_dormance(ipts,ivm) .GE. gdd_crit(ipts) .AND. &
                  t2m_week(ipts) .GT. t2m_month(ipts) )
             WRITE(numout,*) 'temperature test 2', &
                  t2m_month(ipts) .GT. t_always(ivm)
             WRITE(numout,*) 'time_hum_min, hum_min_time', &
                  time_hum_min(ipts,test_pft), &
                  hum_min_time(test_pft)
             WRITE(numout,*) 'vegstress_week, vegstress_months', &
                  vegstress_week(ipts,test_pft), &
                  vegstress_month(ipts,test_pft)
             WRITE(numout,*) 'vegstress_month, relsoilmoist_always', &
                  vegstress_month(ipts,test_pft), &
                  relsoilmoist_always(ivm)
             WRITE(numout,*) 'moisture test 1', &
                  ( time_hum_min(ipts,ivm) .GT. hum_min_time(ivm) .AND. &
                  vegstress_week(ipts,ivm) .GT. vegstress_month(ipts,ivm) )
             WRITE(numout,*) 'moisture test 2', &
                  vegstress_month(ipts,ivm) .GE. relsoilmoist_always(ivm)
             WRITE(numout,*) 'pheno_moigdd_t_crit, ', &
                  pheno_moigdd_t_crit(test_pft)
             WRITE(numout,*) 't2m_month, ', t2m_month(ipts), &
                  ZeroCelsius + pheno_moigdd_t_crit(test_pft)
             WRITE(numout,*) 'temperature test 3', &
                  ( pheno_moigdd_t_crit(ivm) == undef .OR. &
                  t2m_month(ipts) .GT. (ZeroCelsius + pheno_moigdd_t_crit(ivm)) )
          END IF
       END IF
       !-

       !ORCH2.0
       !IF ( PFTpresent(ipts,ivm) .AND. &
       !     plant_status(ipts,ivm) .EQ. ibudsavail .AND. &
       !     ( gdd_m5_dormance(ipts,ivm) .NE. undef)) THEN

       !TRUNK
       !IF ( PFTpresent(ipts,ivm) .AND. &
       !     plant_status(ipts,ivm) .EQ. ibudsavail .AND. &
       !     ( gdd_m5_dormance(ipts,ivm) .NE. undef .OR. &
       !     t2m_month(ipts) .GT. t_always(ivm) ) ) THEN

       ! Using the flag moigdd_orchidee_2_0. Combine the original ORCH2.0 and
       ! TRUNK statements (see above) in a single expression
       IF ( ( moigdd_orchidee_2_0 .AND. PFTpresent(ipts,ivm) .AND. &
            plant_status(ipts,ivm) .EQ. ibudsavail .AND. &
            ( gdd_m5_dormance(ipts,ivm) .NE. undef) ) .OR. &
            ( (.NOT.moigdd_orchidee_2_0) .AND. &
            PFTpresent(ipts,ivm) .AND. &
            plant_status(ipts,ivm) .EQ. ibudsavail .AND. &
            ( gdd_m5_dormance(ipts,ivm) .NE. undef .OR. &
            t2m_month(ipts) .GT. t_always(ivm) ) ) ) THEN

          !! 2.1 Calculate the critical GDD_M5_DORMANCE (::gdd_crit), which is a function of the PFT-dependent 
          !!     critical GDD and the "long term" 2 meter air temperatures 
          tl(ipts) = t2m_longterm(ipts) - ZeroCelsius
          gdd_crit(ipts) = pheno_gdd_crit(ivm,1) + tl(ipts)*pheno_gdd_crit(ivm,2) + &
               tl(ipts)*tl(ipts)*pheno_gdd_crit(ivm,3)

          ! Fail safe option for phenology. Calculate the last day of the year 
          ! that bud break should happen. As long as the mean_start_gs = 0
          ! the fail safe option cannot be used.
          doy_force_pheno = zero
          IF (mean_start_gs(ipts,ivm) .GT. zero) THEN
             IF (mean_start_gs(ipts,ivm) + force_pheno(ivm) .GT. 365) THEN
                doy_force_pheno = mean_start_gs(ipts,ivm) + force_pheno(ivm) - 365
             ELSE
                doy_force_pheno = mean_start_gs(ipts,ivm) + force_pheno(ivm)
             END IF
          ENDIF

          !! 2.2 Determine if the growing season should start (if so, plant_status = ibudbreak).
          ! This occurs if the critical gdd_m5_dormance (::gdd_crit) has been reached 
          ! AND that is temperature increasing, which is true either if the monthly
          ! temperature being higher than the threshold ::t_always, OR if the weekly
          ! temperature is higher than the monthly, 
          ! AND finally that there is sufficient moisture availability, which is 
          ! the same condition as for the ::pheno_moi onset model.
          ! AND when pheno_moigdd_t_crit is set(for C4 grass), if the average 
          ! temperature threshold is reached.  Note that
          ! pheno_moigdd_t_crit can be set to distinguish between tropical,
          ! temperate, and boreal PFTs.
          
          ! ORCH 2.0
          IF (moigdd_orchidee_2_0) THEN
             IF ( ( gdd_m5_dormance(ipts,ivm) .GE. gdd_crit(ipts) ) .AND. &
                  ( ( t2m_week(ipts) .GT. t2m_month(ipts) ) .OR. &
                  ( t2m_month(ipts) .GT. t_always(ivm) )  ) .AND. &
                  ( ( ( time_hum_min(ipts,ivm) .GT. hum_min_time(ivm) ) .AND. &
                  ( vegstress_week(ipts,ivm) .GT. vegstress_month(ipts,ivm) ) ) .OR. & 
                  ( vegstress_month(ipts,ivm) .GE. relsoilmoist_always(ivm) ) ) .AND. &
                  ( ( pheno_moigdd_t_crit(ivm) == undef ) .OR. &
                  (t2m_month(ipts) .GT. (ZeroCelsius + pheno_moigdd_t_crit(ivm))))) THEN
                   
                plant_status(ipts,ivm) = ibudbreak
                doy_start_gs(ipts,ivm) = julian_diff
                grow_season_len(ipts,ivm) = zero
                valid_start_gs(ipts,ivm) = .TRUE.

             ELSEIF (ok_force_pheno .AND. &
                  doy_force_pheno .EQ. julian_diff) THEN
                
                plant_status(ipts,ivm) = ibudbreak
                doy_start_gs(ipts,ivm) = julian_diff
                grow_season_len(ipts,ivm) = zero
                valid_start_gs(ipts,ivm) = .FALSE.

             END IF
          
          ! TRUNK
          ELSEIF ( ( ( gdd_m5_dormance(ipts,ivm) .GE. gdd_crit(ipts) .AND. &
               t2m_week(ipts) .GT. t2m_month(ipts) ) .OR. &
               t2m_month(ipts) .GT. t_always(ivm) ) .AND. & 
               ( ( time_hum_min(ipts,ivm) .GT. hum_min_time(ivm) .AND. &
               vegstress_week(ipts,ivm) .GT. vegstress_month(ipts,ivm) ) .OR. &
               vegstress_month(ipts,ivm) .GE. relsoilmoist_always(ivm) ) .AND. &
               ( pheno_moigdd_t_crit(ivm) == undef .OR. &
               t2m_month(ipts) .GT. (ZeroCelsius + pheno_moigdd_t_crit(ivm)) ) ) THEN
             
             plant_status(ipts,ivm) = ibudbreak
             doy_start_gs(ipts,ivm) = julian_diff
             grow_season_len(ipts,ivm) = zero
             valid_start_gs(ipts,ivm) = .TRUE.

          ELSEIF (ok_force_pheno .AND. &
               doy_force_pheno .EQ. julian_diff) THEN

             ! Force phenology but make sure this value is not 
             ! used in the calculation of the mean doy for budbreak
             ! for this PFT
             plant_status(ipts,ivm) = ibudbreak
             doy_start_gs(ipts,ivm) = julian_diff
             grow_season_len(ipts,ivm) = zero
             valid_start_gs(ipts,ivm) = .FALSE.

          END IF
     
          ! Debug
          IF (printlev_loc>=4)THEN
             IF(ipts == test_grid .AND. ivm == test_pft)THEN
                WRITE(numout,*) 'plant_status, ',plant_status(ipts,test_pft)
                WRITE(numout,*) 'doy_start_gs, ',doy_start_gs(ipts,ivm)
             END IF
          END IF
          !-  

       ENDIF ! PFT there and start of growing season allowed

    ENDDO

    IF (printlev>=4) WRITE(numout,*) 'Leaving moigdd'

  END SUBROUTINE pheno_moigdd


!! ================================================================================================================================
!! SUBROUTINE   : pheno_ncdgdd
!!
!>\BRIEF          The Number of Chilling Days - Growing Degree Day (NCD-GDD) model initiates 
!!                leaf onset if a certain relationship between the number of chilling days (NCD) 
!!                since leaves were lost, and the growing degree days (GDD) since midwinter, is 
!!                fulfilled. 
!!                Currently PFT 6 (Temperate Broad-leaved Summergreen) and PFT 8 (Boreal Broad-
!!                leaved Summergreen) are assigned to this model.
!! 
!! DESCRIPTION  : Experiments have shown that some
!!                species have a "chilling" requirement, i.e. their physiology needs cold 
!!                temperatures to trigger the mechanism that will allow the following budburst 
!!                (e.g. Orlandi et al., 2004). 
!!                An increase in chilling days, defined as a day with a daily mean air
!!                temperature below a PFT-dependent threshold, reduces a plant's GDD demand 
!!                (Cannell and Smith, 1986; Murray et al., (1989); Botta et al., 2000).
!!                The GDD threshold therefore decreases as NCD 
!!                increases, using the following empirical negative explonential law:
!!                \latexonly
!!                \input{phenology_ncdgdd_gddmin_eqn5.tex}
!!                \endlatexonly
!!                \n
!!                The constants used have been calibrated against data CHECK FOR REFERENCE OR PERSON WHO DID UPDATE.
!!                Leaf onset begins if the GDD is higher than the calculated minimum GDD
!!                (dependent upon NCD) AND if the weekly temperature is higher than the monthly 
!!                temperature. This is to ensure the temperature is increasing. \n
!!                The dormancy time-length is represented by the variable 
!!                ::time_lowgpp, which is calculated in ::stomate_season. It is increased by 
!!                the stomate time step when the weekly GPP is lower than a threshold. Otherwise
!!                it is set to zero. \n
!!                The NCD (::ncd_dormance) is calculated in ::stomate_season as  
!!                the number of days with a temperature below a PFT-dependent constant threshold
!!                (::ncdgdd_temp), starting from the beginning of the dormancy period
!!                (::time_lowgpp>0), i.e. since the leaves were lost. \n
!!                The growing degree day sum of the temperatures higher than 
!!                ::ncdgdd_temp (GDD) since midwinter (::gdd_midwinter) 
!!                is also calculated in ::stomate_season.
!!                Midwinter is detected if the monthly temperature is lower than the weekly
!!                temperature AND  the monthly temperature is lower than the long-term
!!                temperature. ::gdd_minter is therefore set to 0 at the beginning of midwinter
!!                and increased with each temperature greater than the PFT-dependent threshold.
!!                When midsummer is detected (the opposite of the above conditions), 
!!                ::gdd_midwinter is set to undef.
!!                CHECK! WHEN TO START OF DORMANCY BEEN MODIFIED FROM BOTTA- ADD IN?
!!                The ::pheno_ncdgdd subroutine is called in the subroutine ::phenology. 
!!
!! RECENT CHANGE(S): None
!!                
!! MAIN OUTPUT VARIABLE(S): ::plant_status - specifies whether leaf growth can start
!!
!! REFERENCE(S) :
!! - Botta, A., N. Viovy, P. Ciais, P. Friedlingstein and P. Monfray (2000), 
!! A global prognostic scheme of leaf onset using satellite data,
!! Global Change Biology, 207, 337-347.
!! - Cannell, M.J.R. and R.I. Smith (1986), Climatic warming, spring budburst and
!! frost damage on trees, Journal of Applied Ecology, 23, 177-191.
!! - Krinner, G., N. Viovy, N. de Noblet-Ducoudre, J. Ogee, J. Polcher, P. 
!! Friedlingstein, P. Ciais, S. Sitch and I.C. Prentice (2005), A dynamic global
!! vegetation model for studies of the coupled atmosphere-biosphere system, Global
!! Biogeochemical Cycles, 19, doi:10.1029/2003GB002199.
!! - Murray, M.B., G.R. Cannell and R.I. Smith (1989), Date of budburst of fifteen 
!! tree species in Britain following climatic warming, Journal of Applied Ecology,
!! 26, 693-700.
!! - Orlandi, F., H. Garcia-Mozo, L.V. Ezquerra, B. Romano, E. Dominquez, C. Galan,
!! and M. Fornaciari (2004), Phenological olive chilling requirements in Umbria
!! (Italy) and Andalusia (Spain), Plant Biosystems, 138, 111-116. 
!!
!! FLOWCHART    : 
!! \latexonly
!! \includegraphics[scale = 1]{pheno_ncdgdd.png}
!! \endlatexonly
!! \n
!_ ================================================================================================================================

  SUBROUTINE pheno_ncdgdd (npts, ivm, PFTpresent, &
       ncd_dormance, gdd_midwinter, &
       t2m_month, t2m_week, plant_status, &
       mean_start_gs, doy_start_gs, &
       valid_start_gs, grow_season_len)

    !
    !! 0. Variable and parameter declaration
    !

    !
    !! 0.1 Input variables
    !
    INTEGER(i_std), INTENT(in)                          :: npts            !! Domain size - number of grid cells (unitless)
    INTEGER(i_std), INTENT(in)                          :: ivm             !! PFT index (unitless)
    LOGICAL, DIMENSION(:,:), INTENT(in)                 :: PFTpresent      !! PFT exists (true/false)
    REAL(r_std), DIMENSION(:,:), INTENT(in)             :: ncd_dormance    !! number of chilling days since leaves were lost 
                                                                           !! (days) 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)          :: gdd_midwinter   !! growing degree days since midwinter (C)
    REAL(r_std), DIMENSION(:), INTENT(in)               :: t2m_month       !! "monthly" 2-meter temperatures (K)
    REAL(r_std), DIMENSION(:), INTENT(in)               :: t2m_week        !! "weekly" 2-meter temperatures (K)
    REAL(r_std),DIMENSION(:,:),INTENT(in)               :: mean_start_gs   !! mean growing season starting day for 
                                                                           !! deciduous PFTs (doy).

    !
    !! 0.2 Output variables
    !
    REAL(r_std),DIMENSION(:,:),INTENT(out)              :: doy_start_gs    !! growing season starting day of year (DOY) for 
    LOGICAL, DIMENSION(:,:),INTENT(out)                 :: valid_start_gs  !! Doy was calculated by ORCHIDEE. It is therefore valid 
                                                                           !! and should be used to update mean_start_gs (0-1)

    !
    !! 0.3 Modified variables
    !
    REAL(r_std),DIMENSION(:,:), INTENT(inout)           :: plant_status    !! Growth and phenological status of the plant
                                                                           !! Different stati are listed in constantes_var
    REAL(r_std),DIMENSION(:,:),INTENT(inout)            :: grow_season_len !! growing season length in days for deciduous PFTs.

    !
    !! 0.4 Local variables
    !
    INTEGER(i_std)                                           :: ipts            !! index (unitless)
    REAL(r_std)                                              :: gdd_min         !! critical gdd (C)
    REAL(r_std)                                              :: doy_force_pheno !! Last day of the year that phenology should happen (doy)
!_ ================================================================================================================================

    IF (printlev>=3) WRITE(numout,*) 'Entering ncdgdd'

    !
    !! 1. Initializations
    !


    !
    !! 1.2 check the critical value ::ncdgdd_temp is defined.
    !!     If not, stop.
    !

    IF ( ncdgdd_temp(ivm) .EQ. undef ) THEN

       WRITE(numout,*) 'ncdgdd: ncdgdd_temp is undefined for PFT (::j) ',ivm
       CALL ipslerr_p(3,'stomate phenology','ncdgdd_temp this PFT','','')

    ENDIF

    !
    !! 2. Check if biometeorological conditions are favourable for leaf growth.    
    !!    PFT has to be there and start of growing season must be allowed.
    !

    DO ipts = 1, npts ! loop over grid points

       ! Debug
       IF(printlev_loc>=4)THEN
          IF(ipts == test_grid .AND. ivm == test_pft)THEN
             WRITE(numout,*) 'Phenology, ',ivm
             WRITE(numout,*) 'Phenology, plant_status, ',plant_status(ipts,ivm)
             WRITE(numout,*) 'Phenology, gdd_midwinter, ',gdd_midwinter(ipts,ivm)
             WRITE(numout,*) 'Phenology, ncd_dormance, ',ncd_dormance(ipts,ivm)
             WRITE(numout,*) 'Phenology, temp, ',t2m_week(ipts),t2m_month(ipts)
          END IF
       END IF
       !-

       IF ( PFTpresent(ipts,ivm) .AND. &
            plant_status(ipts,ivm) .EQ. ibudsavail .AND. &
            ( gdd_midwinter(ipts,ivm) .NE. undef ) .AND. &
            ( ncd_dormance(ipts,ivm) .NE. undef )                  ) THEN

          !! 2.1 Calculate the critical gdd, which is related to ::ncd_dormance
          !!     using an empirical negative exponential law as described above. 
          gdd_min = ( gddncd_ref / exp(gddncd_curve*ncd_dormance(ipts,ivm)) - &
               gddncd_offset )

          ! Debug
          IF(ipts == test_grid .AND. ivm == test_pft)THEN
             IF(printlev_loc>=4) WRITE(numout,*) 'Phenology, gdd_min, ',gdd_min
          END IF
          !-

          ! Fail safe option for phenology. Calculate the last day of the year 
          ! that bud break should happen. As long as the mean_start_gs = 0
          ! the fail safe option cannot be used.
          doy_force_pheno = zero
          IF (mean_start_gs(ipts,ivm) .GT. zero) THEN
             IF (mean_start_gs(ipts,ivm) + force_pheno(ivm) .GT. 365) THEN
                doy_force_pheno = mean_start_gs(ipts,ivm) + force_pheno(ivm) - 365
             ELSE
                doy_force_pheno = mean_start_gs(ipts,ivm) + force_pheno(ivm)
             END IF
          ENDIF

          !! 2.2 Determine if the growing season should start 
          !  (if so, plant_status = ibudbreak).
          !  This occurs if the critical GDD been reached AND the 
          !  temperatures are increasing. If the growing season has 
          !  started, ::gdd_midwinter is set to "undef". 
          IF ( ( gdd_midwinter(ipts,ivm) .GE. gdd_min ) .AND. &
               ( t2m_week(ipts) .GT. t2m_month(ipts) ) ) THEN

             plant_status(ipts,ivm) = ibudbreak
             gdd_midwinter(ipts,ivm)=undef
             doy_start_gs(ipts,ivm) = julian_diff
             grow_season_len(ipts,ivm) = zero
             valid_start_gs(ipts,ivm) = .TRUE.
             
          ELSEIF (ok_force_pheno .AND. &
               doy_force_pheno .EQ. julian_diff) THEN
             
             ! Force phenology but make sure this value is not 
             ! used in the calculation of the mean doy for budbreak
             ! for this PFT
             plant_status(ipts,ivm) = ibudbreak
             gdd_midwinter(ipts,ivm)=undef
             doy_start_gs(ipts,ivm) = julian_diff
             grow_season_len(ipts,ivm) = zero
             valid_start_gs(ipts,ivm) = .FALSE.

          ENDIF

       ENDIF        ! PFT there and start of growing season allowed

    ENDDO ! end loop over grid points

    IF (printlev>=4) WRITE(numout,*) 'Leaving ncdgdd'

  END SUBROUTINE pheno_ncdgdd


!! ================================================================================================================================
!! SUBROUTINE   : pheno_ngd
!!
!>\BRIEF          The Number of Growing Days (NGD) leaf onset model initiates leaf onset if the NGD, 
!!                defined as the number of days with temperature above a constant threshold, 
!!                exceeds a critical value.
!!                Currently PFT 9 (Boreal Leedleleaf Summergreen) is assigned to this model. 
!!
!! DESCRIPTION    The NGD model is a variant of the GDD model. The model was proposed by Botta et
!!                al. (2000) for boreal and arctic biomes, and is designed to estimate 
!!                leaf onset after the end of soil frost. 
!!                The NDG (::ngd_minus5) is the number of days with a daily mean air 
!!                temperature of greater than -5 degrees C, 
!!                starting from the beginning of the dormancy period (i.e. time since the leaves 
!!                were lost/GPP below a certain threshold).
!!                Leaf onset begins if the NGD is higher than the PFT-dependent constant threshold, 
!!                ::ngd,  AND if the weekly temperature is higher than the monthly 
!!                temperature. \n
!!                The dormancy time-length is represented by the variable 
!!                ::time_lowgpp, which is calculated in ::stomate_season. It is increased by 
!!                the stomate time step when the weekly GPP is lower than a threshold. Otherwise
!!                it is set to zero. \n
!!                ::ngd_minus5 is also calculated in ::stomate_season. It is initialised at the
!!                beginning of the dormancy period (::time_lowgpp>0), and increased by the 
!!                stomate time step when the temperature > -5 degrees C. \n
!!                ::ngd is set for each PFT in ::stomate_data, and a 
!!                table defining the minimum NGD for each PFT is given in ::ngd_crit_tab
!!                in ::stomate_constants. \n
!!                The ::pheno_ngd subroutine is called in the subroutine ::phenology.      
!!
!! RECENT CHANGE(S): None
!!                
!! MAIN OUTPUT VARIABLE(S): ::plantstatus - specifies whether leaf growth can start
!!
!! REFERENCE(S) : 
  !! - Botta, A., N. Viovy, P. Ciais, P. Friedlingstein and P. Monfray (2000), 
!! A global prognostic scheme of leaf onset using satellite data,
!! Global Change Biology, 207, 337-347. 
!! - Krinner, G., N. Viovy, N. de Noblet-Ducoudre, J. Ogee, J. Polcher, P. 
!! Friedlingstein, P. Ciais, S. Sitch and I.C. Prentice (2005), A dynamic global
!! vegetation model for studies of the coupled atmosphere-biosphere system, Global
!! Biogeochemical Cycles, 19, doi:10.1029/2003GB002199.
!!
!! FLOWCHART    : 
!! \latexonly
!! \includegraphics[scale = 1]{pheno_ngd.png}
!! \endlatexonly
!! \n
!_ ================================================================================================================================

  SUBROUTINE pheno_ngd (npts, ivm, PFTpresent, ncd_dormance, &
       ngd_minus5, t2m_month, t2m_week, plant_status, &
       mean_start_gs, doy_start_gs, &
       valid_start_gs, grow_season_len)

    !
    !! 0. Variable and parameter declaration
    !

    !
    !! 0.1 Input variables
    !
    INTEGER(i_std), INTENT(in)                          :: npts            !! Domain size - number of grid cells (unitless)
    INTEGER(i_std), INTENT(in)                          :: ivm             !! PFT index (unitless)
    LOGICAL, DIMENSION(:,:), INTENT(in)                 :: PFTpresent      !! PFT exists (true/false)
    REAL(r_std), DIMENSION(:,:), INTENT(in)             :: ncd_dormance    !! number of chilling days since leaves were lost 
                                                                           !! (days) 
    REAL(r_std), DIMENSION(:), INTENT(in)                 :: t2m_month     !! "monthly" 2-meter temperatures (K)
    REAL(r_std), DIMENSION(:), INTENT(in)                 :: t2m_week      !! "weekly" 2-meter temperatures (K)
    REAL(r_std),DIMENSION(:,:),INTENT(in)               :: mean_start_gs   !! mean growing season starting day for 
                                                                           !! deciduous PFTs (doy).

    !
    !! 0.2 Output variables
    !
    REAL(r_std),DIMENSION(:,:),INTENT(out)              :: doy_start_gs    !! growing season starting day of year (DOY) for 
                                                                           !! deciduous PFTs.
    LOGICAL, DIMENSION(:,:),INTENT(out)                 :: valid_start_gs  !! Doy was calculated by ORCHIDEE. It is therefore valid 
                                                                           !! and should be used to update mean_start_gs (0-1)

    !
    !! 0.3 Modified variables
    !
    REAL(r_std),DIMENSION(:,:), INTENT(inout)           :: plant_status    !! Growth and phenological status of the plant
                                                                           !! Different stati are listed in constantes_var
    REAL(r_std), DIMENSION(:,:), INTENT(inout)          :: ngd_minus5      !! growing degree since midwinter days (C)
    REAL(r_std),DIMENSION(:,:),INTENT(inout)            :: grow_season_len !! growing season length in days for deciduous PFTs.

    !
    !! 0.4 Local variables
    !
    INTEGER(i_std)                                           :: ipts            !! index (unitless)
    REAL(r_std)                                              :: doy_force_pheno !! Last day of the year that phenology should happen (doy)

    !! =========================================================================

    IF (printlev>=3) WRITE(numout,*) 'Entering ngd'

    !
    !! 1. Initializations

    !! 1.2 check the critical value ::ngd_crit is defined.
    !!     If not, stop.
    !

    IF ( ngd_crit(ivm) .EQ. undef ) THEN

       WRITE(numout,*) 'ngd: ngd_crit is undefined for PFT (::j) ',ivm
       CALL ipslerr_p(3,'stomate phenology',&
            'ngd_crit is undefined for this PFT','','')

    ENDIF

    !
    !! 2. Check if biometeorological conditions are favourable for leaf growth.
    !!    PFT has to be there and start of growing season must be allowed.
    
    ! Debug
    IF (printlev_loc>=4) THEN
       IF(ivm == test_pft)THEN
          WRITE(numout,*) 'phenology PFTpresent, ', ivm, PFTpresent(test_grid,ivm)
          WRITE(numout,*) 'phenology plant_status, ', plant_status(test_grid,ivm)
          WRITE(numout,*) 'phenology ngd, ', ngd_minus5(test_grid,ivm)
          WRITE(numout,*) 'phenology ngd_crit, ',ngd_crit(ivm) 
          WRITE(numout,*) 'phenology ncd_dormance, ', ncd_dormance(test_grid,ivm)
          WRITE(numout,*) 'phenology ngd_min_dormance, ',ngd_min_dormance
          WRITE(numout,*) 'phenology t2m_week, ', t2m_week(test_grid)
          WRITE(numout,*) 'phenology ncdgdd_temp, ', ncdgdd_temp(ivm)+ZeroCelsius
          WRITE(numout,*) 'phenology t2m_month, ', t2m_month(test_grid)
       END IF
    END IF
    !-

    DO ipts = 1, npts

       IF ( PFTpresent(ipts,ivm) .AND. &
            plant_status(ipts,ivm)  .EQ. ibudsavail ) THEN

          ! Fail safe option for phenology. Calculate the last day of the year 
          ! that bud break should happen. As long as the mean_start_gs = 0
          ! the fail safe option cannot be used.
          doy_force_pheno = zero
          IF (mean_start_gs(ipts,ivm) .GT. zero) THEN
             IF (mean_start_gs(ipts,ivm) + force_pheno(ivm) .GT. 365) THEN
                doy_force_pheno = mean_start_gs(ipts,ivm) + force_pheno(ivm) - 365
             ELSE
                doy_force_pheno = mean_start_gs(ipts,ivm) + force_pheno(ivm)
             END IF
          ENDIF

          !! 2.1 Determine if the growing season should start 
          !  (if so, plant_status = ibudbreak). This occurs if 
          !  the critical NGD has been reached AND are temperatures 
          !  increasing.

          
          ! Trunk phenology
          IF ( ( ngd_minus5(ipts,ivm) .GE. ngd_crit(ivm)) .AND. &
               ( t2m_week(ipts) .GT. t2m_month(ipts) ) ) THEN
          
          ! ORCHIDEE-CAN phenology
!!$          IF ( ngd_minus5(ipts,ivm) .GE. ngd_crit(ivm) .AND. &
!!$               ncd_dormance(ipts,ivm) .GE. ngd_min_dormance .AND. &
!!$               t2m_week(ipts) .GT. ncdgdd_temp(ivm)+ZeroCelsius .AND. &
!!$               t2m_week(ipts) .GT. t2m_month(ipts) ) THEN

             plant_status(ipts,ivm) = ibudbreak
             ngd_minus5(ipts,ivm) = undef
             doy_start_gs(ipts,ivm) = julian_diff
             grow_season_len(ipts,ivm) = zero
             valid_start_gs(ipts,ivm) = .TRUE.
             
          ELSEIF (ok_force_pheno .AND. &
               doy_force_pheno .EQ. julian_diff) THEN
             
             ! Force phenology but make sure this value is not 
             ! used in the calculation of the mean doy for budbreak
             ! for this PFT
             plant_status(ipts,ivm) = ibudbreak
             ngd_minus5(ipts,ivm) = undef
             doy_start_gs(ipts,ivm) = julian_diff
             grow_season_len(ipts,ivm) = zero
             valid_start_gs(ipts,ivm) = .FALSE.

          ENDIF

       ENDIF        ! PFT there and start of growing season allowed

    ENDDO ! end loop over grid points

    IF (printlev>=4) WRITE(numout,*) 'Leaving ngd'

  END SUBROUTINE pheno_ngd

  !! ================================================================================================================================
!! SUBROUTINE   : pheno_siggdd
!!
!>\BRIEF          The 'siggdd' onset model based on moigdd, it is an optional model for peat PFT, 
!!                the mechanism are same to moigdd except that the gdd_crit function is sigmoid. This is based on Qiu 2019 paper. XW applied modifications necessary for the trunk.
!_ ================================================================================================================================
  SUBROUTINE pheno_siggdd (npts, ivm, PFTpresent,gdd_m5_dormance, &
       time_hum_min, &
       t2m_longterm, t2m_month, t2m_week, &
       vegstress_week, vegstress_month, &
       plant_status, mean_start_gs, doy_start_gs, &
       valid_start_gs, grow_season_len)
    !
    !! 0. Variable and parameter declaration
    !

    !
    !! 0.1 Input variables
    !
    INTEGER(i_std), INTENT(in)                               :: npts            !! Domain size - number of grid cells (unitless)
    INTEGER(i_std), INTENT(in)                               :: ivm             !! PFT index (unitless)
    LOGICAL, DIMENSION(npts,nvm), INTENT(in)                 :: PFTpresent      !! PFT exists (true/false)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)             :: gdd_m5_dormance !!! growing degree days, calculated since leaves 
                                                                                !! have fallen (C) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)             :: time_hum_min    !! time elapsed since strongest moisture 
                                                                                !! availability (days) 
    REAL(r_std), DIMENSION(npts), INTENT(in)                 :: t2m_longterm    !! "long term" 2 meter temperatures (K)
    REAL(r_std), DIMENSION(npts), INTENT(in)                 :: t2m_month       !! "monthly" 2-meter temperatures (K)
    REAL(r_std), DIMENSION(npts), INTENT(in)                 :: t2m_week        !! "weekly" 2-meter temperatures (K)
    ! modified
    REAL(r_std), DIMENSION(:,:), INTENT(in)                  :: vegstress_week   !! "weekly" moisture availability (0-1, unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(in)                  :: vegstress_month  !! "monthly" moisture availability (0-1, unitless)

    !
    !! 0.2 Output variables
    !
    !
    !! 0.3 Modified variables
    ! modified
    LOGICAL, DIMENSION(npts,nvm), INTENT(inout)              :: valid_start_gs ! original name: begin_leaves    !! signal to start putting leaves on (true/false)

    !
    !! 0.4 Local variables
    ! 
    REAL(r_std),DIMENSION(:,:),INTENT(in)               :: mean_start_gs   !! mean growing season starting day for 
                                                                           !! deciduous PFTs (doy).
    ! !
    ! !! 0.2 Output variables from trunk
    ! !
    REAL(r_std),DIMENSION(:,:),INTENT(out)              :: doy_start_gs    !! growing season starting day of year (DOY) for 
    !                                                                        !! deciduous PFTs.
    ! !
    ! !! 0.3 Two Modified variables from trunk
    ! !
    REAL(r_std),DIMENSION(:,:),INTENT(inout)            :: plant_status    !! Growth and phenological status of the plant
    !                                                                        !! Different stati are listed in constantes_var
    REAL(r_std),DIMENSION(:,:),INTENT(inout)            :: grow_season_len !! growing season length in days for deciduous PFTs.
    !
    !
    ! local variables
    REAL(r_std), DIMENSION(npts)                             :: tl                              !! long term temperature (C)
    REAL(r_std), DIMENSION(npts)                             :: gdd_crit                        !! critical GDD (C)
    INTEGER(i_std)                                           :: i, ipts !, ivm                    !! index (unitless)
    ! added by following the trunk
    REAL(r_std)                                              :: doy_force_pheno                 !! Last day of the year that phenology should happen (doy)


!_ ================================================================================================================================

    IF (printlev>=3) WRITE(numout,*) 'Entering siggdd'

    !
    !! 1. Initializations
    !

    !
    !! 1.1 first call - outputs the name of the onset model, the values of the  
    !!     moisture availability parameters for tree and grass, and the value of the 
    !!     critical monthly temperature.
    !

    IF ( firstcall_siggdd ) THEN

       WRITE(numout,*) 'pheno_siggdd:'
       WRITE(numout,*) '   > moisture availability above which moisture tendency doesn''t matter: '
       WRITE(numout,*) '         trees (::moiavail_always_tree) :', relsoilmoist_always(5)
       WRITE(numout,*) '         grasses (::moiavail_always_grass) :', relsoilmoist_always(10)
       WRITE(numout,*) '   > monthly temp. above which temp. tendency doesn''t matter (::t_always): ', &
            t_always

       firstcall_siggdd = .FALSE.
    ENDIF

    !
    !! 1.2 initialize output
    !
    ! modified trunk
    valid_start_gs(:,ivm) = .FALSE.

    !
    !! 1.3 check the critical values ::gdd and ::pheno_crit_hum_min_time are defined.
    !!     If not, stop.
    !

    IF ( ANY(pheno_gdd_crit(ivm,:) .EQ. undef) ) THEN

       WRITE(numout,*) 'siggdd: pheno_gdd_crit is undefined for PFT',ivm
       CALL ipslerr_p(3,'stomate phenology','pheno_gdd is undefined for this PFT','','')

    ENDIF

    IF ( hum_min_time(ivm) .EQ. undef ) THEN

       WRITE(numout,*) 'siggdd: hum_min_time is undefined for PFT',ivm
       CALL ipslerr_p(3,'stomate phenology','hum_min is undefined for this PFT','','')

    ENDIF
    !
    !! 2. Check if biometeorological conditions are favourable for leaf growth.
    !!    The PFT has to be there, the start of growing season must be allowed, 
    !!    and GDD has to be defined.
    !
    DO ipts = 1, npts

       IF ( PFTpresent(ipts,ivm) .AND. &
            plant_status(ipts,ivm) .EQ. ibudsavail .AND. &
            ( gdd_m5_dormance(ipts,ivm) .NE. undef .OR. &
            t2m_month(ipts) .GT. t_always(ivm) ) ) THEN

          !! 2.1 Calculate the critical GDD (::gdd_crit), which is a function of the PFT-dependent 
          !!     critical GDD and the "long term" 2 meter air temperatures 
          
          tl(ipts) = t2m_longterm(ipts) - ZeroCelsius
          gdd_crit(ipts)=1.92717865E+5/(1+EXP(-8.13142588E-2*(tl(ipts)-8.78682785E+1)))
          
          !! 2.2 Determine if the growing season should start (if so, ::begin_leaves set to TRUE).
          !!     This occurs if the critical gdd (::gdd_crit) has been reached 
          !!     AND that is temperature increasing, which is true either if the monthly
          !!     temperature being higher than the threshold ::t_always, OR if the weekly
          !!     temperature is higher than the monthly, 
          !!     AND finally that there is sufficient moisture availability, which is 
          !!     the same condition as for the ::pheno_moi onset model.
          !!     AND when pheno_moigdd_t_crit is set(for C4 grass), if the average temperature threshold is reached
          !
          ! 
          doy_force_pheno = zero 
          IF (mean_start_gs(ipts,ivm) .GT. zero) THEN
             IF (mean_start_gs(ipts,ivm) + force_pheno(ivm) .GT. 365) THEN
                doy_force_pheno = mean_start_gs(ipts,ivm) + force_pheno(ivm) - 365
             ELSE
                doy_force_pheno = mean_start_gs(ipts,ivm) + force_pheno(ivm)
             END IF
          ENDIF

          ! gdd_m5_dormance : number of days since the leaf fall: 
          IF ( ( gdd_m5_dormance(ipts,ivm) .GE. gdd_crit(ipts) ) &
               .AND. &
               ( ( t2m_week(ipts) .GT. t2m_month(ipts) ) .OR. &
                 ( t2m_month(ipts) .GT. t_always(ivm) )  ) .AND. &
               ( ( ( time_hum_min(ipts,ivm) .GT. hum_min_time(ivm) ) .AND. &
                 ( vegstress_week(ipts,ivm) .GT. vegstress_month(ipts,ivm) ) )  .OR. &
                 ( vegstress_month(ipts,ivm) .GE. relsoilmoist_always(ivm))  ) .AND. &
               ( ( pheno_moigdd_t_crit(ivm) .EQ. undef ) .OR. &
                 (t2m_month(ipts) .GT. (ZeroCelsius + pheno_moigdd_t_crit(ivm))) ) ) THEN

             plant_status(ipts,ivm) = ibudbreak
             doy_start_gs(ipts,ivm) = julian_diff
             grow_season_len(ipts,ivm) = zero
             valid_start_gs(ipts,ivm) = .TRUE.

          ELSEIF (ok_force_pheno .AND. &
               doy_force_pheno .EQ. julian_diff) THEN

             ! Force phenology but make sure this value is not 
             ! used in the calculation of the mean doy for budbreak
             ! for this PFT
             plant_status(ipts,ivm) = ibudbreak
             doy_start_gs(ipts,ivm) = julian_diff
             grow_season_len(ipts,ivm) = zero
             valid_start_gs(ipts,ivm) = .FALSE.

          ENDIF
          
          ! Debug
          IF (printlev_loc>=4)THEN
             IF(ipts == test_grid .AND. ivm == test_pft)THEN
                WRITE(numout,*) 'plant_status, ',plant_status(ipts,test_pft)
                WRITE(numout,*) 'doy_start_gs, ',doy_start_gs(ipts,ivm)
             END IF
          END IF

       ENDIF        ! PFT there and start of growing season allowed

    ENDDO

    IF (printlev>=4) WRITE(numout,*) 'Leaving siggdd'

  END SUBROUTINE pheno_siggdd


END MODULE stomate_phenology
