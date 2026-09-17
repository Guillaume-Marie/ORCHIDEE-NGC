! =================================================================================================================================
! MODULE       : sapiens_lcchange
!
! CONTACT      : orchidee-help _at_ ipsl.listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       Impact of land cover change on carbon stocks
!!
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S)	: None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_stomate/sapiens_lcchange.f90 $
!! $Date: 2026-05-14 11:18:13 +0200 (jeu. 14 mai 2026) $
!! $Revision: 9542 $
!! \n
!_ ================================================================================================================================

MODULE sapiens_lcchange

  ! modules used:
  USE ioipsl_para
  USE stomate_data
  USE pft_parameters
  USE dynamic_parameters
  USE constantes
  USE grid
  USE stomate_prescribe
  USE function_library,    ONLY: wood_to_dia, wood_to_ba, &
                                 check_vegetation_area, check_area_invariant, &
                                 check_mass_balance, sort_circ_class_biomass,&
                                 check_pixel_area, check_area_change, &
                                 check_variable_change, check_variable_swap, &
                                 check_biomass_change, get_printlev, &
                                 weibull_class_dist
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC land_cover_change_main, age_class_distr, merge_biomass_pfts, check_loss_gain, &
       check_veget, adjust_delta_veget_max, check_read_vegetmax, calculate_shares, &
       move_pft_properties

  INTEGER(i_std), SAVE       :: printlev_loc                          !! Local level of text output for this module
!$OMP THREADPRIVATE(printlev_loc)
  LOGICAL, SAVE              :: firstcall_sapiens_lcchange = .TRUE.   !! first call flag
!$OMP THREADPRIVATE(firstcall_sapiens_lcchange)
  ! Guillaume M. -- Counters for the sub-threshold classes routed to litter by
  ! merge_biomass_pfts. Reported and reset once a year by age_class_distr: without them
  ! a clean year cannot tell "the path works" from "the path was never taken".
  INTEGER(i_std), SAVE       :: orphan_litter_hits = 0                !! sub-threshold classes routed since last report
!$OMP THREADPRIVATE(orphan_litter_hits)
  REAL(r_std), SAVE          :: orphan_litter_mass = zero             !! carbon routed that way @tex $(g C m^{-2})$ @endtex
!$OMP THREADPRIVATE(orphan_litter_mass)

CONTAINS


! ================================================================================================================================
!! SUBROUTINE   : land_cover_change_main
!!
!>\BRIEF        Impact of land cover change on carbon stocks
!!
!! DESCRIPTION  : This subroutine is activated by setting VEGET_UPDATE to 1Y 
!! in the sechiba.card file. Subsequently, the impact of land cover change on 
!! in-situ carbon stocks is computed in this subroutine. Fluxes from products are
!! now dealt with in a separate subroutine CHECK. The land cover change is read from
!! driver maps containing the new "maximal" coverage fraction of a PFT (::veget_max). 
!! On the basis of this difference, the amount of 'new establishment'/'biomass export',
!! and increase/decrease of each component, are estimated. The soil carbon from the
!! PFT which experienced a decrease in ::veget_max is pooled and used as the initial
!! soil carbon for the PFTs that experienced an increase in ::veget_max.\n
!!
!! RECENT CHANGE(S) : (1) Age-classes were introduced and (2) product use was moved 
!! to a separate module
!!
!! MAIN OUTPUT VARIABLE(S) :  
!!
!! REFERENCES   : None
!!
!! FLOWCHART    : 
!! \n
!_ ================================================================================================================================

  SUBROUTINE land_cover_change_main (&
       npts,                 dt_days,              veget_max,                  veget, &
       harvest_pool,         harvest_type,         harvest_cut,                harvest_area, &
       litter,               som,                  lignin_struc,               lignin_wood, &
       PFTpresent,           everywhere,           when_growthinit,            leaf_frac, &
       circ_class_n,         circ_class_biomass,   atm_to_bm,                  forest_managed, &
       KF,                   wstress_season,       wstress_month,              plant_status, &
       npp_longterm,         croot_longterm,       age,                        lm_lastyearmax, &
       harvest_pool_bound,   bm_to_litter,         turnover_daily,             leaf_age, &
       longevity_eff_leaf,   longevity_eff_sap,    longevity_eff_root, &
       veget_max_new,        loss_gain,            age_stand,                  last_cut, &
       age_stand_bm,         age_stand_area, &
       k_latosa_adapt,       losses,               light_tran_to_floor_season, lpft_replant, &
       soil_n_min,           bact,                 species_change_map,         cn_leaf_init_2D, &
       bm_sapl_2D,           tree_bm_to_litter,    fLulccResidue,              fDeforestToProduct, &
       deepSOM_a,            deepSOM_s,            deepSOM_p,                  sugar_load, &
       frac_nobio,           frac_nobio_new,       zf_soil,                    burried_litter, &
       burried_fresh_ltr,    burried_fresh_som,    burried_bact,               burried_min_nitro, &
       burried_som,          burried_deepSOM_a,    burried_deepSOM_s,          burried_deepSOM_p, &
       burried_adj,          kill_vessels,         vessel_loss_previous,       biomass_init_drought, &
       mean_start_gs,        woody_litter_by_cut,  lignin_snag)

    IMPLICIT NONE
    
 !! 0. Variable and parameter declaration 
    
    !! 0.1 Input variables
    ! NOTE: if harvest_pool_bound is not explicitely defined 
    ! the model uses 1:ndia_harvest+2 as its dimensions. If 
    ! you like to do so, remember to change the counters in 
    ! this routine
    INTEGER, INTENT(in)                                   :: npts                   !! Domain size - number of pixels (unitless)
    INTEGER(i_std), DIMENSION (npts,nvm), INTENT(in)      :: forest_managed         !! forest management flag
    REAL(r_std), INTENT(in)                               :: dt_days                !! Time step of vegetation dynamics for stomate (days)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)          :: veget                  !! Fraction of vegetation type including 
                                                                                    !! non-biological fraction (unitless)   
    REAL(r_std), DIMENSION(0:ndia_harvest+1), INTENT(in)  :: harvest_pool_bound     !! The boundaries of the diameter classes
                                                                                    !! in the wood harvest pools @tex $(m)$ @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)          :: longevity_eff_root     !! Effective root turnover time that accounts
                                                                                    !! waterstress (days)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)          :: longevity_eff_sap      !! Effective sapwood turnover time that accounts
                                                                                    !! waterstress (days)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)          :: longevity_eff_leaf     !! Effective leaf turnover time that accounts
                                                                                    !! waterstress (days)
    REAL(r_std), DIMENSION(npts,nvm),INTENT(in)           :: veget_max_new          !! new "maximal" coverage fraction of a PFT on the  
                                                                                    !! ground.  May sum to
                                                                                    !! less than unity if the pixel has
                                                                                    !! nobio area. (unitless, 0-1)
    CHARACTER(15), INTENT(in)                              :: losses                !! Flag that determines how the losses are 
                                                                                    !! distributed. There are currently two options
                                                                                    !! 'proportional' and 'youngest'
    LOGICAL, DIMENSION(npts,nvm), INTENT(in)               :: lpft_replant          !! Set to true if a PFT has been clearcut
    INTEGER(i_std), DIMENSION (npts,nvm), INTENT(in)       :: species_change_map    !! A map which gives the PFT number that each PFT
    REAL(r_std), DIMENSION(npts,nvm),INTENT(inout)         :: light_tran_to_floor_season !! Mean seasonal fraction of light transmitted
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)           :: loss_gain             !! losses and gains due to LCC distributed over all
                                                                                    !! age classes and thus taking the age-classes into 
                                                                                    !! account (unitless, 0-1)
    REAL(r_std),DIMENSION(npts,nnobio),INTENT(in)          :: frac_nobio_new        !! New fraction of nobio per gridcell
    REAL(r_std), DIMENSION(0:ngrnd), INTENT(in)            :: zf_soil

    !! 0.2 Output variables
    
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)          :: fDeforestToProduct    !! Deforested biomass into product pool due to anthorpogeninc
                                                                                    !! land use change 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)          :: fLulccResidue         !! carbon mass flux into soil and litter due to anthropogenic land use or land cover change
    REAL(r_std), DIMENSION(npts,nelements),INTENT(out)     :: burried_adj           !! Temporary variables to simplify the calcualtions
                                                                                    !! with C and N burried under urbanized land (gC pixel-1)

                                                                              
    !! 0.3 Modified variables   
    
    LOGICAL, DIMENSION(npts,nvm), INTENT(inout)            :: PFTpresent            !! PFT present (0 or 1)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)        :: plant_status          !! Growth and phenological status of the plant
                                                                                    !! istatus = Phases defined in constantes
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)        :: KF                    !! Scaling factor to convert sapwood mass into leaf 
                                                                                    !! mass (m)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)        :: k_latosa_adapt        !! Leaf to sapwood area adapted for long 
                                                                                    !! term water stress (m)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)        :: veget_max             !! "maximal" coverage fraction of a PFT on the ground
                                                                                    !! May sum to
                                                                                    !! less than unity if the pixel has
                                                                                    !! nobio area. (unitless, 0-1)
    REAL(r_std), DIMENSION(npts,nvm,nelements), INTENT(inout) :: atm_to_bm          !! C and N taken from the atmosphere to get C and N to  
                                                                                    !! create the seedlings @tex (gC.m^{-2}dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)        :: everywhere            !! is the PFT everywhere in the grid box or 
                                                                                    !! very localized (after its introduction) (?)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)        :: when_growthinit       !! how many days ago was the beginning of 
                                                                                    !! the growing season (days)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)        :: wstress_season        !! Water stress factor, based on hum_rel_daily
                                                                                    !! (unitless, 0-1)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)        :: wstress_month         !! Water stress factor, based on hum_rel_daily
                                                                                    !! (unitless, 0-1)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)        :: npp_longterm          !! "long term" net primary productivity
                                                                                    !! @tex ($gC m^{-2} year^{-1}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)        :: croot_longterm        !! "long term" root carbon mass  
 	                                                                            !! @tex ($gC m^{-2}) @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)        :: age                   !! mean age (years)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)        :: lm_lastyearmax        !! last year's maximum leaf mass for each PFT 
                                                                                    !! @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm,nlevs), INTENT(inout)  :: lignin_struc          !! ratio Lignine/Carbon in structural litter,
                                                                                    !! above and below ground
    REAL(r_std), DIMENSION(npts,nvm,nlevs), INTENT(inout)  :: lignin_wood           !! ratio Lignine/Carbon in woody litter,
                                                                                    !! above and below ground
    REAL(r_std), DIMENSION(npts,nvm,nlevs), INTENT(inout)  :: lignin_snag           !! ratio Lignine/Carbon in snag litter,
                                                                                    !! above and below ground
    REAL(r_std), DIMENSION(npts,ncarb,nvm,nelements), INTENT(inout) :: som          !! carbon pool: active, slow, or passive 
                                                                                    !! @tex ($gC m^{-2}$) @endtex 
    REAL(r_std), DIMENSION(npts,nlitt,nvm,nlevs,nelements), INTENT(inout)  :: litter!! metabolic and structural litter, above and 
                                                                                    !! below ground @tex ($gC m^{-2}$) @endtex 
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)      :: circ_class_n               !! Number of individuals in each circ class
                                                                                    !! @tex $(ind m^{-2})$ @endtex   
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)      :: leaf_frac                  !! fraction of leaves in leaf age class (unitless;0-1)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)      :: leaf_age                   !! Leaf age (days)
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout):: bm_to_litter !! Transfer of biomass to litter
                                                                                    !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout):: tree_bm_to_litter !! conversion of biomass to litter
                                                                                         !! @tex ($gC m^{-2} day^{-1}$) @endtex
    
                                                                               
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)    ::turnover_daily              !! Turnover rates            
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)  :: circ_class_biomass         !! Biomass components of the model tree within a 
                                                                                    !! circumference class @tex $(g C ind^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)  :: harvest_pool          !! The wood and biomass that have been
                                                                               !! havested by humans @tex $(gC)$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)        :: harvest_type          !! Type of management that resulted
                                                                               !! in the harvest (unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)        :: harvest_cut           !! Type of cutting that was used for the harvest
                                                                               !! (unitless)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)      :: harvest_area          !! Harvested area (m^{2})
    INTEGER(i_std), DIMENSION(:,:), INTENT(inout)     :: age_stand             !! Age of the forest stand (years)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)        :: age_stand_bm          !! Biomass-weighted conserved mean stand age (years) - STAND_AGE
    REAL(r_std), DIMENSION(:,:), INTENT(inout)        :: age_stand_area        !! AREA-weighted conserved mean stand age (years) - STAND_AGE
    INTEGER(i_std), DIMENSION(:,:), INTENT(inout)     :: last_cut              !! Years since last thinning (years)
    REAL(r_std),DIMENSION(:,:), INTENT(in)            :: cn_leaf_init_2D       !! initial leaf C/N ratio 
    REAL(r_std),DIMENSION(:,:,:,:,:), INTENT(inout)   :: bm_sapl_2D            !! biomass of sapling
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)    :: deepSOM_a             !! Soil carbon discretized with depth active (g/m**3) 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)    :: deepSOM_s             !! Soil carbon discretized with depth slow (g/m**3) 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)    :: deepSOM_p             !! Soil carbon discretized with depth passive (g/m**3)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)      :: soil_n_min            !! mineral nitrogen in the soil (gN/m**2)  
    REAL(r_std), DIMENSION(:,:), INTENT(inout)        :: bact                  !! denitrifier biomass (gC/m**2)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)        :: sugar_load            !! Relative sugar loading of the labile pool (unitless)
    REAL(r_std),DIMENSION(:,:),INTENT(inout)          :: frac_nobio            !! Fraction of grid cell covered by lakes, land 
                                                                               !! ice, cities, ... (unitless) 
    REAL(r_std),DIMENSION(:,:,:,:),INTENT(inout)      :: burried_litter        !! Litter burried under non-biological land uses (gC or N m-2)
    REAL(r_std),DIMENSION(:,:,:),INTENT(inout)        :: burried_fresh_ltr     !! Fresh litter burried under non-biological land uses (gC or N m-2)
    REAL(r_std),DIMENSION(:,:,:),INTENT(inout)        :: burried_fresh_som     !! Fresh som burried under non-biological land uses (gC or N m-2)
    REAL(r_std),DIMENSION(:),INTENT(inout)            :: burried_bact          !! Bacteria burried under non-biological land uses (gC m-2)
    REAL(r_std),DIMENSION(:,:),INTENT(inout)          :: burried_min_nitro     !! Mineral nitrogen burried under non-biological land uses (gC or N m-2)
    REAL(r_std),DIMENSION(:,:,:),INTENT(inout)        :: burried_som           !! Som burried under non-biological land uses (gC or N m-2)
    REAL(r_std),DIMENSION(:,:,:),INTENT(inout)        :: burried_deepSOM_a     !! Som burried under non-biological land uses (gC or N m-2)
    REAL(r_std),DIMENSION(:,:,:),INTENT(inout)        :: burried_deepSOM_s     !! Som burried under non-biological land uses (gC or N m-2)
    REAL(r_std),DIMENSION(:,:,:),INTENT(inout)        :: burried_deepSOM_p     !! Som burried under non-biological land uses (gC or N m-2)
    LOGICAL, DIMENSION(:,:,:), INTENT(inout)          :: kill_vessels          !! Flag to kill vessels at the end of the day following embolism.
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)      :: vessel_loss_previous  !! Proportion of conductivity lost due to cavitation, accumulated
                                                                               !!  on the previous day (unitless).
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)  :: biomass_init_drought  !! Biomass of heartwood or sapwood before onset of drought
    REAL(r_std), DIMENSION(:,:), INTENT(inout)        :: mean_start_gs         !! mean growing season starting day for deciduous PFTs (doy).
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)      :: woody_litter_by_cut   !! Saved woody litter by icut not by iparts @tex $(gC m^{-2})$ @endtex
 
    
    !! 0.4 Local variables
    INTEGER(i_std)                                   :: ipts,ivm,imbc          !! Indices (unitless)
    INTEGER(i_std)                                   :: igroup,ifm,igrn        !! Indices (unitless)
    INTEGER(i_std)                                   :: idia,icir,ilage        !! Indices (unitless)
    INTEGER(i_std)                                   :: jcir,ilit,ilev         !! Indices (unitless)
    INTEGER(i_std)                                   :: ipar,iele,icarb        !! Indices (unitless)
    INTEGER(i_std)                                   :: ivma,iyoung,jvma       !! Indices (unitless)
    INTEGER(i_std)                                   :: nagec_temp,ispec       !! Temporary counter
    LOGICAL                                          :: lsync                  !! returns false if the biomass is not synced
    REAL(r_std)                                      :: total_losses           !! Sum of losses only in delta_veget (unitless, 0-1)
    REAL(r_std)                                      :: share_rec              !! Share of the veget_max of the existing vegetation
                                                                               !! within a PFT over the total veget_max following 
                                                                               !! expansion of that PFT (unitless, 0-1)
    REAL(r_std), DIMENSION(nlitt,nlevs,nelements)    :: litter_bank            !! Bank to store the litter that becomes available when
                                                                               !! part of a PFT is removed @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(ncarb,nelements)          :: soil_bank              !! Bank to store soil carbon that becomes available when
                                                                               !! part of a PFT is removed  @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(ngrnd,ncarb,nelements)    :: som_vertres_bank       !! Bank to store the vertically-resolved SOM dilution (gC/m**2)  

    REAL(r_std), DIMENSION(nlevs)                    :: struct_ltr_bank        !! Bank to store the lignin in structural litter that 
                                                                               !! becomes available when part of a PFT is removed 
                                                                               !! (0-1,unitless)
    REAL(r_std),DIMENSION(nlevs)                     :: woody_ltr_bank         !! Bank to store the lignin in woody litter that 
                                                                               !! becomes available when part of a PFT is removed 
                                                                               !! (0-1,unitless)
    REAL(r_std),DIMENSION(nlevs)                     :: snag_ltr_bank          !! Bank to store the lignin in snag litter that 
                                                                               !! becomes available when part of a PFT is removed 
                                                                               !! (0-1,unitless)
    REAL(r_std)                                      :: fresh_woody_ltr_bank   !! Bank to store fresh woody litter to moved to woody_litter_by_cut
                                                                               !! when part of a PFT is removed (0-1)
    REAL(r_std),DIMENSION(npts,nvm,nparts,nelements) :: fresh_litter           !! Pool of fresh litter that is being released during 
    REAL(r_std),DIMENSION(npts,nvm)                  :: fresh_litter_wood      !! Pool of fresh woody litter is being released during LCC
                                                                               !! @tex ($gC m^{-2}$) @endtex
    REAL(r_std),DIMENSION(nvm,nparts,nelements)      :: fresh_som              !! Pool of fresh som that is being released during
                                                                               !! LCC @tex ($gC m^{-2}$) @endtex 
    REAL(r_std),DIMENSION(nparts,nelements)          :: fresh_ltr_bank         !! Bank to store all the fresh litter from site
    REAL(r_std),DIMENSION(nparts,nelements)          :: fresh_som_bank         !! Bank to store all the fresh som from site 
    REAL(r_std)                                      :: bact_bank              !! Bank to store the bacteria
    REAL(r_std), DIMENSION(nnspec)                   :: soil_n_min_bank        !! Bank to store mineral soil nitrogen
                                                                               !! clearing during LCC. Weighted value of fresh_litter
    REAL(r_std), DIMENSION(npts,nvm,nelements)       :: atm_to_bm_old          !!
    
                                                                               !! @tex ($gC m^{-2}$) @endtex
    REAL(r_std)                                      :: moi_bank_season        !! bank to store vegstress_season of the cleared
                                                                               !! PFT's (unitless; 0-1)
    REAL(r_std)                                      :: moi_bank_month         !! bank to store vegstress_monthly of the cleared
                                                                               !! PFT's (unitless; 0-1)
    REAL(r_std)                                      :: trees_needed           !! Number of trees leftthat should be present in a 
                                                                               !! circumference class when merging the available and 
                                                                               !! newly established vegetation
                                                                               !! @tex $(ind m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(ncirc)                    :: diameters_temp         !! Temporary variable to store the diameters of the
                                                                               !! different diameter classes (m)
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements):: new_biomass            !! Temporary variable for biomass 
                                                                               !! @tex $(gC.m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nleafages)       :: new_leaf_frac          !! Temporary variable for leaf_frac (unitless;0-1)
    REAL(r_std), DIMENSION(ncirc)                    :: est_circ_class_n       !! Temporary variable for circ_class_n of the 
                                                                               !! established vegetation @tex $(ind m^{-2})$ @endtex
    REAL(r_std), DIMENSION(ncirc)                    :: tmp_circ_class_n       !! Temporary variable for circ_class_n of the 
                                                                               !! established vegetation @tex $(ind m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,ncirc)           :: new_circ_class_n       !! Temporary variable for circ_class_n of the
                                                                               !! @tex $(ind m^{-2})$ @endtex
    REAL(r_std), DIMENSION(ncirc,nparts,nelements)   :: est_circ_class_biomass !! Temporary variable for circ_class_biomass of the
                                                                               !! established vegetation @tex $(g C ind^{-1})$ @endtex
    REAL(r_std), DIMENSION(ncirc,nparts,nelements)   :: tmp_circ_class_biomass !! Temporary variable for circ_class_biomass of the
                                                                               !! established vegetation @tex $(g C ind^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,ncirc,nparts,nelements) :: new_circ_class_biomass !! Temporary variable for circ_class_biomass of the
                                                                                      !! new vegetation @tex $(g C ind^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,ncirc)           :: new_circ               !! Temporary variable for circ (m)
    REAL(r_std), DIMENSION(npts,nvm,nelements)       :: new_atm_to_bm          !! Temporary variable for atm_to_bm 
                                                                               !! @tex (gC.m^{-2}dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)                 :: new_age                !! Temporary variable for age (years)
    REAL(r_std), DIMENSION(npts,nvm)                 :: new_KF                 !! Temporary variable for KF
    REAL(r_std), DIMENSION(npts,nvm)                 :: new_lm_lastyearmax     !! Temporary variable for lm_last_year 
                                                                               !! @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm)                 :: new_everywhere         !! Temporary variable for everywhere (unitless, 0-1)
    REAL(r_std), DIMENSION(npts,nvm)                 :: dum_when_growthinit    !! Dummy for when_growthinit (days)
    REAL(r_std), DIMENSION(npts,nvm)                 :: dum_npp_longterm       !! Dummy for npp_longterm
                                                                               !! @tex ($gC m^{-2} year^{-1}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm)                 :: new_croot_longterm     !! Temporary variable for croot_longterm
                                                                               !! @tex ($gC m^{-2} year^{-1}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm,nmbcomp,nelements):: check_intern          !! Contains the components of the internal
                                                                               !! mass balance chech for this routine
                                                                               !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)       :: closure_intern         !! Check closure of internal mass balance
                                                                               !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)       :: pool_start             !! Start and end pool of this routine 
                                                                               !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)       :: pool_end               !! Start and end pool of this routine 
                                                                               !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    LOGICAL, DIMENSION(npts,nvm)                     :: dum_PFTpresent         !! Dummy for PFTpresent (0 or 1)
    REAL, DIMENSION(npts,nvm)                        :: dum_plant_status       !! Dummy for plant_status
    REAL(r_std)                                      :: harvest_ratio_stem     !! The fraction of the reserve pools which are
                                                                               !! located in the stem (unitless, 0-1)
    REAL(r_std)                                      :: harvest_ratio_branch   !! The fraction of the reserve pools which are
                                                                               !! located in the branches/roots (unitless, 0-1)
    REAL(r_std)                                      :: harvest_ratio_litter   !! The fraction of the reserve pools which are
                                                                               !! moved to the litter (unitless, 0-1)
    REAL(r_std)                                      :: total_nonheart_biomass !! The total amount of non-heartwood biomass
                                                                               !! in a model tree
                                                                               !! @tex $(gC tree^{-1})$ @endtex
    REAL(r_std)                                      :: total_nonheart_biomass_stem !! The total amount of non-heartwood biomass
                                                                               !! in a model tree stem
                                                                               !! @tex $(gC tree^{-1})$ @endtex
    REAL(r_std)                                      :: total_nonheart_biomass_branch !! The total amount of non-heartwood biomass
                                                                               !! in a model tree branches/roots
                                                                               !! @tex $(gC tree^{-1})$ @endtex
    REAL(r_std)                                      :: temp_sum,temp_total
    REAL(r_std), DIMENSION(npts)                     :: change_bio             !! The overall change in vegetative fraction for
                                                                               !! this pixel (0-1, unitless)
     REAL(r_std), DIMENSION(nlitt,nlevs)             :: litter_weight_rec      !! The fraction of litter on the expanded
                                                                               !! PFT.
                                                                               !! @tex $-$ @endtex
    REAL(r_std), DIMENSION(nlevs)                    :: litter_weight_struc    !! A weighting factor for the structural litter.
                                                                               !! (unitless)
    REAL(r_std), DIMENSION(nlevs)                    :: litter_weight_woody    !! A weighting factor for the woody litter.
                                                                               !! (unitless)
    REAL(r_std), DIMENSION(nlevs)                    :: litter_weight_snag     !! A weighting factor for the woody litter.
                                                                               !! (unitless)
    LOGICAL                                          :: lconverged             !! A flag for loop convergence
    INTEGER(i_std)                                   :: max_pft                !! The number of the PFT with the largest addition
                                                                               !! in land area.
    INTEGER(i_std)                                   :: min_pft                !! The number of the PFT with the largest reduction
                                                                               !! in land area (most negative value).
    INTEGER(i_std)                                   :: nsteps                 !! The number of convergence steps tried.
    INTEGER(i_std)                                   :: ncount                 !! The number of PFTs with significant changes
                                                                               !! in vegetation.
    INTEGER(i_std)                                   :: nagec_local            !! The number of age classes in this map PFT.
    INTEGER(i_std)                                   :: sv_i, sv_k             !! Sliver guard loop indices
    INTEGER(i_std)                                   :: sv_kmax                !! Slot receiving the snapped area residue
    REAL(r_std)                                      :: sv_res                 !! Area residue snapped to zero (-)
    REAL(r_std)                                      :: sv_n                   !! Trees still standing on a snapped slot
    REAL(r_std)                                      :: difference             !! How much leftover veget_max we have
    REAL(r_std), DIMENSION(npts,nvm)                 :: veget_max_begin        !! Temporary variable to check area conservation 
    REAL(r_std), DIMENSION(npts,nvm,ndia_harvest+1,nelements,nlanduse) :: harvest_pool_old  !! Temporary variable to check mass conservation
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements):: bm_to_litter_old       !! Temporary biomass to litter to check mass conservation
                                                                               !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts)                     :: change_nobio           !! Change in the non biological fraction in a pixel (0-1, unitless)
    REAL(r_std)                                      :: check_start
    REAL(r_std), DIMENSION(nelements,2)              :: test_mass              !! Temporary variable to check mass conservation 
    REAL(r_std), DIMENSION(nelements)                :: test_mass_bank         !! Temporary variable to check mass conservation
    REAL(r_std), DIMENSION(3,nelements,2)            :: test_bank              !! Temporary variable to check mass conservation
    REAL(r_std)                                      :: plant_status_init      !! Initial plant_status of a PFT
    INTEGER(i_std)                                   :: c0r_ivm, c0r_icir      !! CLASS0_REPLANT: indices creneau / classe de circonference
    REAL(r_std)                                      :: c0r_phi                !! CLASS0_REPLANT: fraction du creneau replantee (-)
    REAL(r_std)                                      :: c0r_share              !! CLASS0_REPLANT: poids du peuplement en place, = 1 - phi (-)
    REAL(r_std)                                      :: c0r_nmix               !! CLASS0_REPLANT: densite melangee d'une classe (ind m-2)
    LOGICAL                                          :: c0r_bad_in             !! SONDE C0NAN: une entree de l'accumulation est non finie
    LOGICAL                                          :: c0r_bad_out            !! SONDE C0NAN: bm_to_litter est non finie apres accumulation
    REAL(r_std), DIMENSION(npts,nvm,nelements)       :: check_fresh_litter
    REAL(r_std), DIMENSION(npts,nvm,nelements)       :: check_biomass
    REAL(r_std), DIMENSION(npts,nvm,nelements)       :: init_biomass
    REAL(r_std), DIMENSION(npts,nvm,nelements)       :: fin_biomass
    REAL(r_std), DIMENSION(npts,nvm,nelements)       :: biomass_start          !! biomass before LCC in gC or gN pixel-1
    REAL(r_std), DIMENSION(npts,nvm,nelements)       :: biomass_end            !! biomass after LCC in gC or gN pixel-1
    REAL(r_std), DIMENSION(npts,nvm,nelements)       :: atm_to_bm_loc          !! atmosphere to biomass in this subroutine (gC or gN m-2)
    REAL(r_std), DIMENSION(npts,nelements)           :: burried                !! Temporary variables to simplify the calcualtions
                                                                               !! with C and N burried under urbanized land (gC pixel-1)
    REAL(r_std), DIMENSION(npts,nvm)                 :: new_ind                !! Number of recruits grown (trees m-2 day-2). Not used in sapiens_lcchange.
                                                                               !! Added to avoid using OPTINAL arguments. 
    LOGICAL                                          :: write_debug            !! Flag to write debug statements
    REAL(r_std), DIMENSION(npts,nvm,nelements)       :: tot_res_target         !! Target for reserve carbon pool (only for grassland)
                                                                               !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm,nelements)       :: tot_lab_target         !! Target for labile carbon pool (only for grassland)
                                                                               !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm)                 :: gpp_week               !! "weekly" GPP @tex ($gC m^{-2} day^{-1}$)@endtex
    REAL(r_std), DIMENSION(npts)                     :: t2m                    !! Temperature at 2 meter (K)
                                                                               !! Added to avoid using OPTINAL arguments

    REAL(r_std)                                      :: residual_area          !! temporary variable for the residual area after a LCC or species change

!_ ================================================================================================================================

    IF (firstcall_sapiens_lcchange) THEN
      !! Initialize local printlev
      ! Guillaume M. -- get_printlev calls getin_p, a COLLECTIVE bcast. Called here,
      ! inside the compute loop under a local firstcall, it hangs the run as soon as the
      ! ranks stop reaching it together. Use printlev, read once at startup by every rank.
      printlev_loc=printlev
      firstcall_sapiens_lcchange=.FALSE.
    END IF

      ! Guillaume M. -- get_printlev calls getin_p, a COLLECTIVE bcast. Called here,
      ! inside the compute loop under a local firstcall, it hangs the run as soon as the
      ! ranks stop reaching it together. Use printlev, read once at startup by every rank.
    printlev_loc=printlev
    IF (printlev_loc.GE.2) WRITE(numout,*) 'Entering land cover change main'

    !  +++COMMENT+++
    !  This part is poorly coded as it uses the externalisation to
    !  implement the age classes. The advanatage of this approach is
    !  that is rather simple to implement. The alternative requires
    !  that a dimension for age classes is added to almost all fluxes
    !  and pool in sechiba and stomate. Also the externalisation would
    !  need to be adjusted as the parameters are now used for the
    !  different age classes.
    !
    !  The coding is considered poor because it is : (1) sensitive to 
    !  the order in which the PFTs are defined and (2) inconsistent as
    !  cohorts (within a stand got an own dimension) and age classes
    !  (across stands) are implemented through the externalisation.
    !  An example for 2 age classes and 8 tree PFTs
    !  PFT_1: bare soil
    !  PFT_2: PFT2_1 age class 1 (youngest)
    !  PFT_3: PFT2_2 age class 2 (youngest + 1)
    !  PFT_4: PFT3_1 age class 1
    !  PFT_5: PFT3_2 age class 2
    !  ...
    !  PFT_18: PFT 10 agricultural PFTs don't have age classes
    !  PFT_19: PFT 11 agricultural PFTs don't have age classes
    !  PFT_20: PFT 11 agricultural PFTs don't have age classes
    !  PFT_21: PFT 11 agricultural PFTs don't have age classes
    !  
    !  A PFT group are all age classes a a certain PFT/species
    !
    !  Some of this complexity could be avoided by preparing the maps
    !  for the exact number of age classes we want to deal with
    !  but that just shifts the complexity. It also means that we have
    !  to prescribe the age classes, which is not our goal. 
    !  For now it is assumed the maps only consider the PFTs without 
    !  age classes. So part of this code needs to conver the maps
    !  without age classes into variables in ORCHIDEE with age classes.
    !  +++++++++++

 !! 1. Ultimate checks

    IF(err_act.GT.1)THEN

       ! All initial checks should be done in slowproc right after the map
       ! is being read. If vegetation fractions or frac_nobio is adjusted
       ! afterwards, mass balance problems are unavoidable. Check whether
       ! veget_max, loss_gain and veget_max_new are still consistent.
       ! If not, something happened inbetween slowproc and the call to this
       ! subroutine. Note that the first two checks were also done at the 
       ! end of adjust_delta_veget_max.
       CALL check_area_change('start of sapiens_lcchange', npts, frac_nobio, &
            frac_nobio_new, loss_gain)
       CALL check_pixel_area('start of sapiens_lcchange', npts, veget_max_new, &
            frac_nobio_new)
       CALL check_variable_change('start of sapiens_lcchange', npts, veget_max, &
            veget_max_new, loss_gain)

    ENDIF


 !! 2. Initialization

    ! Change in non-biological fraction. Note that loss_gain contains
    ! the changes in the vegetation fractions
    change_nobio(:) = -SUM(loss_gain(:,:),2)

    ! Initialize
    fLulccResidue(:,:) = zero
    fDeforestToProduct(:,:) = zero
    burried(:,:) = zero
    burried_adj(:,:) = zero

    !! 2.2 Initialize check for mass balance closure
    IF (err_act.GT.1) THEN

       check_intern(:,:,:,:) = zero
       pool_start(:,:,:) = zero
       DO iele = 1,nelements

          ! Biomass pool + bm_to_litter
          DO ipar = 1,nparts
            DO icir=1, ncirc
                pool_start(:,:,iele) = pool_start(:,:,iele) + &
                    circ_class_biomass(:,:,icir,ipar,iele) * &
                    circ_class_n(:,:,icir) * veget_max(:,:)
            ENDDO
            pool_start(:,:,iele) = pool_start(:,:,iele) + &
                (bm_to_litter(:,:,ipar,iele) + turnover_daily(:,:,ipar,iele)) * &
                veget_max(:,:)
          ENDDO

          ! Litter pool (gC m-2) *  (m2 m-2) 
          DO ilit = 1,nlitt
             DO ilev = 1,nlevs
                pool_start(:,:,iele) = pool_start(:,:,iele) + &
                     litter(:,ilit,:,ilev,iele) * veget_max(:,:)
             ENDDO
          ENDDO

          IF (ok_soil_carbon_discretization) THEN
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

          ! The biomass harvest pool is expressed in gC pixel-1 So, it 
          ! shouldn't be multiplied by veget_max but it should be divided
          ! by area to obtain gC m-2.
          DO ipts = 1,npts
             pool_start(ipts,:,iele) = pool_start(ipts,:,iele) + &
                  SUM(SUM(harvest_pool(ipts,:,:,iele,:),3),2) / (area(ipts)*contfrac(ipts))
          ENDDO

       ENDDO
       
       !Denitrifying bacteria (only C)
       pool_start(:,:,icarbon) = pool_start(:,:,icarbon) + &
               bact(:,:) * veget_max(:,:)

       !Nitrogen only
       DO ispec = 1,nnspec
          pool_start(:,:,initrogen) = pool_start(:,:,initrogen) + &
                  soil_n_min(:,:,ispec) * veget_max(:,:)
       ENDDO

       ! Initialize area check and intermediate mass balance checks
       veget_max_begin(:,:) = veget_max(:,:)
       atm_to_bm_old(:,:,:) = atm_to_bm(:,:,:)
       pool_end(:,:,:) = zero

       ! It is difficult to find the source of mass balance problems in
       ! sapiens_lcchange because C and N are moved between PFTs and between
       ! compartments. This check tracks changes in circ_class_biomass. Further
       ! below also changes in litter, som, bact and soil_n_min and the surface
       ! area are tracked. This implies that if there is an imbalance and none
       ! of the intermediate tests failed the best place to start looking are
       ! bm_to_litter, turnover_daily and burried.  
       biomass_start(:,:,:) = zero
       biomass_end(:,:,:) = zero
       harvest_pool_old(:,:,:,:,:) = harvest_pool(:,:,:,:,:)
       bm_to_litter_old(:,:,:,:) = bm_to_litter(:,:,:,:)
       CALL check_biomass_change(npts, dt_days, area, circ_class_biomass, &
            circ_class_n, veget_max, loss_gain, atm_to_bm_old, atm_to_bm, &
            harvest_pool, harvest_pool_old, fresh_litter, bm_to_litter_old, &
            biomass_start, biomass_end, 1)
          
    ENDIF ! err_act .GT. 1

 !! 5. Calculate in-situ land cover changes

    !  NOTE: within the model framework there is large conceptual
    !  difference between management and clear cuts on the one
    !  hand and LCC on the other hand. Within the model LCC
    !  DOES NOT change the biomass, ind, ... of the PFTs that
    !  experience a loss (because we store those in gC m-2). LCC 
    !  only changes the values of veget_max for those PFTs. For
    !  this reason it was decided not to make use of lpj_kill.f90
    !  and sapiens_gap.f90 (as is done for forest management) 
    !  because this would be in conflict with the meaning of the 
    !  variable ::circ_class_kill. It is not impossible to work 
    !  around this issues but it could result in difficulties with
    !  the order of code in more complex cases were part of the PFT 
    !  experiences a surface loss and at the same time is ready to
    !  get thinned or clear cut. It appears cleaner to deal with
    !  LCC up front in a self-contained module.
    

    !! 5.1 Land cover losses within a pixel
    !  The losses need to be calculated first such that the soil
    !  carbon (and other elements) can be put in the soil&litter-bank. 
    !  Subsequently, the soil&litter-bank is used to initialize the 
    !  soil carbon of the new PFTs. The biomass is moved into a 
    !  harvest pool.
    !  NOTE: make sure that all variables for PFT1 are initialized!
    !  If not LCC won't properly work.
    DO ipts = 1, npts ! Loop over # pixels - domain size    

       !! 5.1.1 Initialize the soil&litter-bank
       !  NOTE Some of these variables have a dimension for ::nelements
       !  others don't. Make sure ::nelements is accounted for when 
       !  N is added to the model.
       soil_bank(:,:) = zero
       litter_bank(:,:,:) = zero
       struct_ltr_bank(:) = zero
       woody_ltr_bank(:) = zero
       snag_ltr_bank(:) = zero
       fresh_ltr_bank(:,:) = zero
       fresh_woody_ltr_bank = zero
       fresh_som_bank(:,:) = zero
       moi_bank_season = zero
       moi_bank_month = zero
       bact_bank = zero
       soil_n_min_bank(:) = zero
       fresh_litter(ipts,:,:,:) = bm_to_litter(ipts,:,:,:)
       fresh_litter_wood(:,:) = zero
       fresh_som(:,:,:) = turnover_daily(ipts,:,:,:)
       som_vertres_bank(:,:,:) = zero

       !! 5.5.2 Initialize additional mass balance check
       IF (err_act.GT.1) THEN
          
          !  It is difficult to check for mass balance closure in sapiens_lcchange
          !  This code checks whether mass balance is presevered for the pixels
          !  that loose surface area. It also checks mass conservation for the whole
          !  routine for litter, som, bact and soil_n_min because 
          !  these pools don't receive or donate c/n to other pools.
 
          !  Note that above changes in
          !  circ_class_biomass and surface area are tracked. This implies that
          !  if a mass balance problem occurs and none of the intermediate tests
          !  failed the best place to start looking for a solution is
          !  bm_to_litter, turnover_daily and burried because there are not yet
          !  intermediate tests for these compartments.
          CALL init_mass_balance(test_bank, test_mass,&
               veget_max, litter, som, circ_class_biomass, circ_class_n, &
               turnover_daily, bm_to_litter, bact, soil_n_min, &
               harvest_pool, deepSOM_a, deepSOM_s, deepSOM_p, ipts, zf_soil)

       ENDIF

       !! 5.1.3 Loss in PFT cover, update harvest and soil&litter-C-bank
       !  To calculate the relative change in cover fraction, PFT-specific 
       !  changes need to be normalized by the total change which is in 
       !  this case the total loss. Total losses includes losses in PFT1
       total_losses = SUM(loss_gain(ipts,:), MASK = loss_gain(ipts,:) .LT. zero)

       ! +++CHECK+++
       ! This is a powerful line of code. By adding the change_nobio to the
       ! total losses it is assumed that in case of de-urbanization or glacier
       ! retreat the new soil does not contain any carbon. This code puts some
       ! virgin soils on the C and N bank (and will dillute the C and N 
       ! concentration per m-2). This is in contradiction with the approach in 
       ! which we burry carbon in case of urbanization. If a pixel get urbanized 
       ! the C and N are burried and not decomposed but if that same pixel would 
       ! get de-urbanized the carbon would not be released but is assumed to remain 
       ! burried. To fix this issues burried should become a cumulative variable.
       IF (change_nobio(ipts).LT.zero) THEN
          total_losses = total_losses + change_nobio(ipts)
       ENDIF
       ! +++++++++++

       ! Calculate the factor to weight the lignin ratios by when merging.
       ! This depends on the amount of litter on the PFTs where vegetation
       ! is being cut.
       litter_weight_struc(:)=zero
       litter_weight_woody(:)=zero
       litter_weight_snag(:)=zero
       DO ivm = 1, nvm
          IF ( loss_gain(ipts,ivm) .LE. -min_stomate ) THEN
             DO ilev=1,nlevs
                litter_weight_struc(ilev)=litter_weight_struc(ilev)+&
                     litter(ipts,istructural,ivm,ilev,icarbon) * &
                     ABS(loss_gain(ipts,ivm))
                litter_weight_woody(ilev)=litter_weight_woody(ilev)+&
                     litter(ipts,iwoody,ivm,ilev,icarbon) * &
                     ABS(loss_gain(ipts,ivm))
                litter_weight_snag(ilev)=litter_weight_snag(ilev)+&
                     litter(ipts,isnag,ivm,ilev,icarbon) * &
                     ABS(loss_gain(ipts,ivm))
             ENDDO
          ENDIF
       ENDDO

       DO ivm = 1, nvm

          IF (loss_gain(ipts,ivm).EQ.zero) CYCLE

          ! Initialize intermediate mass balance check
          test_mass_bank(:) = zero

          ! The PFT cover decreases, move wood into the harvest pool, the 
          ! stumps and branches are supposed to be exported from the site 
          ! so they go into the harvest pool as well but a separated pool
          ! was introduced for separate treatment in the product use module.
          ! Leaves, fine roots and fruits are deposited on the soil&litter
          ! bank and will be attributed to the new PFT. This order of 
          ! calculations ensures that the litter ends up in the right PFT
          ! and that it is only decomposed once.         
          IF ( loss_gain(ipts,ivm) .LT. zero) THEN

             !! 5.1.2.1 Biomass losses of trees
             IF (is_tree(ivm)) THEN
             
                ! Store the management stategy of the pixel and PFT
                ifm = forest_managed(ipts,ivm)                

                IF(losses .EQ. 'disturbance' .AND. ifm .EQ. ifm_none) THEN
                    ! From r7635, when the natural disturbances kill more than 30% (check
                    ! replant_threshold) the vegetation fraction is moved to the youngest 
                    ! age class. In this case, we need a woody litter input to estimate
                    ! disturbance interactions better. This vegetation fraction is considered 
                    ! to be dead, so all is moved to the litter when the forest is not managed.
                    DO iele = 1,nelements
                       ! Update fresh litter for redistribution
                       ! the amount of mass will be removed is moved to fresh
                       ! litter
                       DO ipar = 1,nparts
                          fresh_litter(ipts,ivm,ipar,iele) = fresh_litter(ipts,ivm,ipar,iele) + &
                             SUM(circ_class_biomass(ipts,ivm,:,ipar,iele) * circ_class_n(ipts,ivm,:) )
                       ENDDO
                    ENDDO
                    ! Save woody litter. This will be used in pest module
                    ! (stomate_pest) for one of the susceptibility indice.
                    fresh_litter_wood(ipts,ivm) = fresh_litter_wood(ipts,ivm) + SUM( &
                            (circ_class_biomass(ipts,ivm,:,isapabove,icarbon) + &
                            circ_class_biomass(ipts,ivm,:,iheartabove,icarbon)) * &
                            circ_class_n(ipts,ivm,:) )
                ELSE
                
                    ! We will need the trunk diameters for all the circ classes in 
                    ! order to figure out which harvest pool to put them in.
                    diameters_temp(:) = &
                         wood_to_dia(circ_class_biomass(ipts,ivm,:,:,icarbon),&
                         ivm,pipe_tune2(ipts,ivm))
 
                    ! Associated harvest variables
                    ! If several cuts happen in the same year, these variables can no 
                    ! longer be correctly interpreted. 
                    harvest_type(ipts,ivm) = harvest_type(ipts,ivm) + ifm
                    harvest_cut(ipts,ivm) = harvest_cut(ipts,ivm) + icut_lcc_res
                    harvest_area(ipts,ivm,ilcc) = harvest_area(ipts,ivm,ilcc) + &
                       ABS(loss_gain(ipts,ivm) * area(ipts) * contfrac(ipts))

                    ! Now we have the stems diameters.  These are what we want to
                    ! save, and we are going to classify them by PFT type, location, 
                    ! and diameter so that we can compare to FAO statistics and do 
                    ! some more interesting post-treatment (sapiens_product_use.f90).
                    DO icir = 1,ncirc
                      DO iele = 1,nelements
                       ! We want to make sure we divide up the labile and carbohydrate
                       ! reserve pools in the same was as in other management.  We
                       ! assume that when a tree is cut down, it cannot rearrange 
                       ! these pools, and therefore the amount that is harvested
                       ! is the ratio that is stored in the harvest wood.  We assume
                       ! that the reserve pools are stored equally in all non-heart
                       ! tissue.  The harvest right now goes into three different
                       ! places: two harvest pools and litter.  We need the fractions
                       ! for the reserve pools going into each place.
                       total_nonheart_biomass = &
                           circ_class_biomass(ipts,ivm,icir,isapbelow,iele)+&
                           circ_class_biomass(ipts,ivm,icir,iroot,iele)+&
                           circ_class_biomass(ipts,ivm,icir,ileaf,iele)+&
                           circ_class_biomass(ipts,ivm,icir,ifruit,iele)+&
                           circ_class_biomass(ipts,ivm,icir,isapabove,iele)

                       ! This includes the roots and the branches, since it goes into
                       ! a different pool.
                       total_nonheart_biomass_branch=branch_ratio(ivm)*&
                           circ_class_biomass(ipts,ivm,icir,isapabove,iele)+&
                           circ_class_biomass(ipts,ivm,icir,isapbelow,iele)
                       total_nonheart_biomass_stem=(un-branch_ratio(ivm))*&
                           circ_class_biomass(ipts,ivm,icir,isapabove,iele)

                       ! These are the three ratios we need to move the reserve 
                       ! pools around. When dealing with standard land cover change 
                       ! the harvest_ratios are always defined. However, when 
                       ! dealing with a species change we try to harvest an empty 
                       ! age class, hence, the code to capture the undefined case. 
                       ! Why is the age class empty. Well in forestry we just 
                       ! harvested an older age class, in kill we replace that 
                       ! older age class by a younger age class but here we harvest 
                       ! it again before it was ever replanted.
                       IF (total_nonheart_biomass .GT. min_stomate) THEN
                          harvest_ratio_branch=total_nonheart_biomass_branch/&
                             total_nonheart_biomass
                          harvest_ratio_stem=total_nonheart_biomass_stem/&
                             total_nonheart_biomass
                          harvest_ratio_litter=un - harvest_ratio_stem - &
                             harvest_ratio_branch
                       ELSE
                          harvest_ratio_branch=zero
                          harvest_ratio_stem=zero
                          harvest_ratio_litter=un - harvest_ratio_stem - &
                             harvest_ratio_branch 
                       END IF

                       dia_loop:DO idia=1,ndia_harvest

                          ! Move the harvest to the correct diameter class and
                          ! account for the largest diameter class because
                          ! harvest_pool_bound has only ndia_harvest+1 dimensions
                          IF ( (diameters_temp(icir).LE.harvest_pool_bound(idia)).OR. &
                             ((diameters_temp(icir).GT.harvest_pool_bound(ndia_harvest)) .AND. &
                             (diameters_temp(icir).LT.harvest_pool_bound(ndia_harvest+1))) ) THEN

                             ! This is the correct diameter class
                             ! Move wood (excluding the branches) and the carbres 
                             ! and labile C contained in the wood into the harvest 
                             ! pool. Note that the unit of circ_class_biomass is gC 
                             ! tree-1, so multiply with the number of trees to obtain 
                             ! gC m-2. To avoid complexities when dealing with the 
                             ! product pool multiply with the cover change 
                             ! and the surface area of the pixel to obtain a weighted 
                             ! mean expressed in gC.
                             harvest_pool(ipts,ivm,idia,iele,ilcc) = &
                                harvest_pool(ipts,ivm,idia,iele,ilcc) + &  
                                ABS(loss_gain(ipts,ivm)) * area(ipts) * contfrac(ipts) * & 
                                ( (un - branch_ratio(ivm)) * circ_class_n(ipts,ivm,icir) * &
                                ( circ_class_biomass(ipts,ivm,icir,isapabove,iele) + &
                                circ_class_biomass(ipts,ivm,icir,iheartabove,iele) ) + &
                                ! reserves
                                harvest_ratio_stem * circ_class_n(ipts,ivm,icir) * &
                                (circ_class_biomass(ipts,ivm,icir,icarbres,iele) + &
                                circ_class_biomass(ipts,ivm,icir,ilabile,iele)) )

                             ! Move stumps and branches and the carbres and labile C 
                             ! contained in them into the harvest pool with the 
                             ! smallest dimenisons. Stumps and branches are removed 
                             ! from the site but can now be treated separately from 
                             ! stem wood in the product use module (sapiens_product_use.f90)
                             harvest_pool(ipts,ivm,1,iele,ilcc) = &
                                harvest_pool(ipts,ivm,1,iele,ilcc) + &
                                ABS(loss_gain(ipts,ivm)) * area(ipts) * contfrac(ipts) * &
                                ! branches
                                (branch_ratio(ivm) * circ_class_n(ipts,ivm,icir) * &
                                ( circ_class_biomass(ipts,ivm,icir,isapabove,iele) + &
                                circ_class_biomass(ipts,ivm,icir,iheartabove,iele) ) + &
                                ! reserves
                                harvest_ratio_branch * circ_class_n(ipts,ivm,icir) * &
                                ( circ_class_biomass(ipts,ivm,icir,icarbres,iele) + &
                                circ_class_biomass(ipts,ivm,icir,ilabile,iele) )  + &
                                ! stumps (all icarbres and ilabile is distributed
                                ! over the stems and branches)
                                circ_class_n(ipts,ivm,icir) * &
                                ( circ_class_biomass(ipts,ivm,icir,isapbelow,iele) + &
                                circ_class_biomass(ipts,ivm,icir,iheartbelow,iele) ) )
                        
                             ! Deposit the leaves, roots and fruit into a temporary 
                             ! litter pool and place it in the soil and litter bank
                             fresh_litter(ipts,ivm,ileaf,iele) = &
                                fresh_litter(ipts,ivm,ileaf,iele) + &
                                circ_class_biomass(ipts,ivm,icir,ileaf,iele) * &
                                circ_class_n(ipts,ivm,icir)
                             fresh_litter(ipts,ivm,iroot,iele) = &
                                fresh_litter(ipts,ivm,iroot,iele) + &
                                circ_class_biomass(ipts,ivm,icir,iroot,iele) * &
                                circ_class_n(ipts,ivm,icir)
                             fresh_litter(ipts,ivm,ifruit,iele) = &
                                fresh_litter(ipts,ivm,ifruit,iele) + &
                                circ_class_biomass(ipts,ivm,icir,ifruit,iele) * &
                                circ_class_n(ipts,ivm,icir)
                             fresh_litter(ipts,ivm,icarbres,iele) = &
                                fresh_litter(ipts,ivm,icarbres,iele) + &
                                harvest_ratio_litter * &
                                circ_class_biomass(ipts,ivm,icir,icarbres,iele) * &
                                circ_class_n(ipts,ivm,icir)
                             fresh_litter(ipts,ivm,ilabile,iele) = &
                                fresh_litter(ipts,ivm,ilabile,iele) + &
                                harvest_ratio_litter * &
                                circ_class_biomass(ipts,ivm,icir,ilabile,iele) * &
                                circ_class_n(ipts,ivm,icir)

                             IF (iele==icarbon) THEN
                                !+++CHECK+++
                                ! Check the definitions of these variables. There is a good chance
                                ! that this information is already contained in the prod or 
                                ! flux_prod variables which distinguish ilcc and iharvest.
                                ! AR6 output. The CN-CAN refinements of these variables were 
                                ! implemented compared to CN, i.e., the branch ratio
                                fLulccResidue(ipts,ivm) = fLulccResidue(ipts,ivm) + &
                                     ABS(loss_gain(ipts,ivm) * &
                                     ! branches
                                     ( branch_ratio(ivm) * circ_class_n(ipts,ivm,icir) * &
                                     ( circ_class_biomass(ipts,ivm,icir,isapabove,iele) + &
                                     circ_class_biomass(ipts,ivm,icir,iheartabove,iele)) + &
                                     ! reserves
                                     harvest_ratio_branch * circ_class_n(ipts,ivm,icir) * &
                                     (circ_class_biomass(ipts,ivm,icir,icarbres,iele) + &
                                     circ_class_biomass(ipts,ivm,icir,ilabile,iele) )  + &
                                     ! stumps (all icarbres and ilabile is distributed
                                     ! over the stems and branches)
                                     circ_class_n( ipts,ivm,icir) * &
                                     ( circ_class_biomass(ipts,ivm,icir,isapbelow,iele) + &
                                     circ_class_biomass(ipts,ivm,icir,iheartbelow,iele) + &
                                     circ_class_biomass(ipts,ivm,icir,iroot,iele) + &
                                     circ_class_biomass(ipts,ivm,icir,ifruit,iele) + &
                                     circ_class_biomass(ipts,ivm,icir,ileaf,iele) ) ) )
                                fDeforestToProduct(ipts,ivm) = fDeforestToProduct(ipts,ivm) + &
                                     ABS(loss_gain(ipts,ivm)) * & 
                                     ( (un - branch_ratio(ivm)) * circ_class_n(ipts,ivm,icir) * &
                                     ( circ_class_biomass(ipts,ivm,icir,isapabove,iele) + &
                                     circ_class_biomass(ipts,ivm,icir,iheartabove,iele)) +&
                                     harvest_ratio_stem * circ_class_n(ipts,ivm,icir) * &
                                     (circ_class_biomass(ipts,ivm,icir,icarbres,iele) + &
                                     circ_class_biomass(ipts,ivm,icir,ilabile,iele) ) )
                                !+++++++++++
                             END IF

                             EXIT dia_loop

                          END IF ! Tree diameter within diameter bounds of the harvest pools

                          IF(diameters_temp(icir) .GE. harvest_pool_bound(ndia_harvest+1))THEN

                             ! If we are here, it means the diameter of this trunk is larger 
                             ! than the upper boundary of our harvest bins. Considering the 
                             ! upper boundary is set to be 99999 meters in constants.f90, 
                             ! something has gone wrong!
                             WRITE(numout,*) 'ipts, ivm, icir: ',ipts, ivm, icir
                             WRITE(numout,*) 'diameters_temp(icir), harvest_pool_bound(idia)',&
                                diameters_temp(icir),harvest_pool_bound(ndia_harvest+1)
                             CALL ipslerr_p (3,'spaien_lcchange', &
                                'The tree is too big to fit into the harvest pools',&
                                'Check the diameter of the tree','')

                          END IF ! tree diameter within diameter bounds of the harvest pools

                       ENDDO dia_loop
                      ENDDO ! nelements
                   ENDDO ! ncirc                
                ENDIF ! losses disturbance

             ELSE ! is_tree

                !! 5.1.2.2 Biomass losses of grasses and crops
                !  During LCC grasses and crops are simply plowed under 
                !  and should thus be move into ::fresh_litter which is
                !  moved later on into the litter bank ::fresh_ltr_bank
                !  and finally into ::bm_to_litter. The bank is needed
                !  to take a weighted average.
                DO iele = 1,nelements
                   
                   ! there is only one circumference class for grasses 
                   !and crops
                   fresh_litter(ipts,ivm,:,iele) = fresh_litter(ipts,ivm,:,iele) + &
                            circ_class_biomass(ipts,ivm,1,:,iele) * &
                            circ_class_n(ipts,ivm,1)   
                   
                ENDDO
                
                fresh_litter_wood(ipts,ivm) = fresh_litter_wood(ipts,ivm) + SUM( &
                     (circ_class_biomass(ipts,ivm,:,isapabove,icarbon) + &
                     circ_class_biomass(ipts,ivm,:,iheartabove,icarbon)) * &
                     circ_class_n(ipts,ivm,:) )
 
             ENDIF ! is_tree(ivm)
             
             !! 5.1.2.3 Deposit soil and litter and soil&litter bank
             ! Deposit the released soil&litter pools on the soil&litter-bank. 
             ! The units of the C-pools are in gC m-2, hence, they need to be 
             ! weighted by their relative change in cover fraction to obtain 
             ! the C-pool (gC pixel-2) of the soil&litter-C-bank.
             ! Mass_removed (g C) = x1*X1 + x2*X2 + x3*X3 where x1, x2, and x3 
             ! are the change in veget_max and X1, X2 and X3 are the 
             ! C density of a given variable (g C m-2). Mass_removed will have
             ! to be redistributed over the PFTs that increased their surface
             ! area. Because we deal with net land cover changes we don't know 
             ! the bi-directional changes and so the mass_removed will be 
             ! averaged and this average will be redistributed. In other words,
             ! Mass_redistributed = y1*Y + y2*Y + y3*Y Where y1, y2 and y3 are 
             ! the increases in veget_max for PFT that gained cover fraction and
             ! Y is the average C-density. Mass conservation dictates that
             ! Mass_removed = Mass_distributed so Y (g C m-2) can de calculated
             ! as (x1*X1+x2*X2+x3*X3)/(y1+y2+y3). Note that y1+y2+y3 = total gain 
             ! and by definition total_gain = total_loss. This is the approach
             ! implemented below.
             litter_bank(:,:,:) = litter_bank(:,:,:) + &
                  loss_gain(ipts,ivm) * litter(ipts,:,ivm,:,:) / total_losses
             fresh_ltr_bank(:,:) = fresh_ltr_bank(:,:) + &
                  loss_gain(ipts,ivm) * fresh_litter(ipts,ivm,:,:) / total_losses
             fresh_som_bank(:,:) = fresh_som_bank(:,:) + &
                  loss_gain(ipts,ivm) * fresh_som(ivm,:,:) / total_losses
             bact_bank = bact_bank + &
                  loss_gain(ipts,ivm) * bact(ipts,ivm) / total_losses
             soil_n_min_bank(:) = soil_n_min_bank(:) + &
                  loss_gain(ipts,ivm) * soil_n_min(ipts,ivm,:) / total_losses
             ! Deposit fresh woody litter. This is part of fresh litter but it is
             ! separately saved for the bark beetles. Also it is different from
             ! woody_ltr_bank which is derived from litter pool.
             fresh_woody_ltr_bank = fresh_woody_ltr_bank + &
                  loss_gain(ipts,ivm) * fresh_litter_wood(ipts,ivm) / total_losses
 
             IF ( ok_soil_carbon_discretization ) THEN
                DO igrn = 1, ngrnd      
                   som_vertres_bank(igrn,iactive,:) = som_vertres_bank(igrn,iactive,:) + & 
                        loss_gain(ipts,ivm) * deepSOM_a(ipts,igrn,ivm,:) * &
                        (zf_soil(igrn)-zf_soil(igrn-1)) / total_losses 
                   som_vertres_bank(igrn,islow,:) = som_vertres_bank(igrn,islow,:) + & 
                        loss_gain(ipts,ivm) * deepSOM_s(ipts,igrn,ivm,:) * &
                        (zf_soil(igrn)-zf_soil(igrn-1)) / total_losses 
                   som_vertres_bank(igrn,ipassive,:) = som_vertres_bank(igrn,ipassive,:) + & 
                        loss_gain(ipts,ivm) * deepSOM_p(ipts,igrn,ivm,:) * &
                        (zf_soil(igrn)-zf_soil(igrn-1)) / total_losses
                ENDDO
             ELSE 
                soil_bank(:,:) =  soil_bank(:,:) + &
                     loss_gain(ipts,ivm) * som(ipts,:,ivm,:) / total_losses
             ENDIF

             DO ilev=1,nlevs

                IF(litter_weight_struc(ilev) .GE. min_stomate)THEN
           
                   struct_ltr_bank(ilev) = struct_ltr_bank(ilev) + &
                        lignin_struc(ipts,ivm,ilev) * &
                        (litter(ipts,istructural,ivm,ilev,icarbon) * &
                        ABS(loss_gain(ipts,ivm))) / &
                        litter_weight_struc(ilev)
                ENDIF

                IF(litter_weight_woody(ilev) .GE. min_stomate)THEN
                   
                   woody_ltr_bank(ilev) = woody_ltr_bank(ilev) + &
                        lignin_wood(ipts,ivm,ilev) * &
                        (litter(ipts,iwoody,ivm,ilev,icarbon) * &
                        ABS(loss_gain(ipts,ivm))) / &
                        litter_weight_woody(ilev)
                ENDIF

                IF(litter_weight_snag(ilev) .GE. min_stomate)THEN
                   
                   snag_ltr_bank(ilev) = snag_ltr_bank(ilev) + &
                        lignin_snag(ipts,ivm,ilev) * &
                        (litter(ipts,isnag,ivm,ilev,icarbon) * &
                        ABS(loss_gain(ipts,ivm))) / &
                        litter_weight_snag(ilev)
                ENDIF

             ENDDO


             ! In stomate_season.f90 this variable is only calculated if 
             ! the biomass(ileaf) > 0. As long as veget_max is 0, this 
             ! variable does not get a value. We need a value
             ! to calculate ::KF which in turn is needed to calculate the biomass
             ! of the seedlings we will prescribe. At this point in the calculation
             ! the vegstress_season from PFT1 is still zero so it is not 
             ! accounted for yet. Its exact value depends on the PFT that will
             ! be established on PFT1 and it will be accounted for later in the 
             ! code
             IF (ivm .GT. 1) THEN

                moi_bank_season = moi_bank_season + loss_gain(ipts,ivm) / &
                     total_losses * ( un - wstress_season(ipts,ivm))
                moi_bank_month = moi_bank_month + loss_gain(ipts,ivm) / &
                     total_losses * (un - wstress_month(ipts,ivm))

             ENDIF

             IF(err_act.GT.1)THEN
                
                ! It is difficult to check for mass balance closure in sapiens_lcchange
                ! because pools, fluxes and surface areas are all changing. This code 
                ! checks whether mass balance is presevered for the pixels that loose 
                ! surface area. It also checks mass conservation for the whole routine 
                ! for litter, som, bact and soil_n_min because these 
                ! pools don't receive or donate c/n to other pools.
                CALL inter_mass_balance(test_bank, test_mass,&
                     veget_max, litter, som, circ_class_biomass, circ_class_n, &
                     turnover_daily, bm_to_litter, bact, soil_n_min, &
                     harvest_pool, change_nobio, loss_gain, litter_bank, &
                     soil_bank, fresh_ltr_bank, fresh_som_bank, total_losses, &
                     test_mass_bank, bact_bank, soil_n_min_bank, &
                     pool_start, deepSOM_a, deepSOM_s, deepSOM_p, som_vertres_bank, &
                     ipts, ivm, dt_days, zf_soil)

             END IF             

             ! This subroutine is used for both land cover changes and species
             ! changes. Both options use the same variables but a slightly different
             ! approach. The best solution is probably to rewrite species change
             ! such that it follows the approach of land cover change. For the
             ! moment a small patch is propsed.
             IF (ok_change_species) THEN
                ! The species change code skips veget_max_new and direclty calculates
                ! loss_gain. Hence, loss_gain needs to be used in this statement.
                ! veget_max_new is expected to be equal to veget_max and loss_gain
                ! equals to minus veget_max_new. To work correctly, residual_area needs
                ! to be calculated as the sum of veget_amx_new and loss_gain.
               residual_area = veget_max_new(ipts,ivm) + loss_gain(ipts,ivm)
             ELSE
                ! The land cover change code reads an LCC map into veget_max_new and
                ! the calculated loss_gain as the difference between veget_max_new
                ! and veget_max. 
                residual_area = veget_max_new(ipts,ivm)
             ENDIF

             ! If the residual area after the change is zero, thus, the PFT was entirely removed
             ! we will clean the PFT. This is essential to avoid a crash on soil_n_min
             ! because we will check for this condition in stomate_litter.f90 and it
             ! is a good practice for the other variables to avoid errors later on.
             IF (residual_area.LT.10*EPSILON(un)) THEN
               
                ! Variables that were deposited on the bank
                litter(ipts,:,ivm,:,:) = zero 
                bact(ipts,ivm) = zero
                soil_n_min(ipts,ivm,:) = zero 
                IF ( ok_soil_carbon_discretization ) THEN
                   deepSOM_a(ipts,:,ivm,:) = zero
                   deepSOM_s(ipts,:,ivm,:) = zero
                   deepSOM_p(ipts,:,ivm,:) = zero
                ELSE 
                   som(ipts,:,ivm,:) = zero
                ENDIF

                ! Other variables that should be reset
                lm_lastyearmax(ipts,ivm) = zero
                croot_longterm(ipts,ivm) = zero
                age(ipts,ivm) = zero 
                leaf_frac(ipts,ivm,:) = zero 
                PFTpresent(ipts,ivm) = .FALSE.
                everywhere(ipts,ivm) = zero 
                lignin_struc(ipts,ivm,:) = zero 
                lignin_wood(ipts,ivm,:) = zero
                lignin_snag(ipts,ivm,:) = zero
                bm_to_litter(ipts,ivm,:,:) = zero
                turnover_daily(ipts,ivm,:,:) = zero 
                atm_to_bm(ipts,ivm,:) = zero 
                circ_class_biomass(ipts,ivm,:,:,:) = zero 
                circ_class_n(ipts,ivm,:) = zero 
                age_stand(ipts,ivm) = 0 
                last_cut(ipts,ivm) = 0 
                plant_status(ipts,ivm) = inone  
                sugar_load(ipts,ivm) = un
               
             END IF

          END IF ! IF (loss_gain(ivm,ipts) < zero)

       END DO ! Loop over # PFT's

       !  Deal with biological land uses being converted into non biological land use
       !  Urbanisation is increasing at the expense of vegetation cover. This is accounted
       !  for in the current generation of land cover change maps and should be dealt with.
       !  When urbanistaion occurs the biomass is harvested and so we put the C and N in 
       !  the harvest pools. It is not clear waht happens in reality with the litter C 
       !  and N. It probably gets mixed or burried during construction. The soil C and N
       !  considered in ORCHIDEE is relatively shallow compared to the depth of many
       !  modern construction works. Again it is not clear what happens with the C and N.
       !  We have two simple options: (1) we put all the carbon in the heterothropic 
       !  respiration at once and bury the N or (2) we bury both C and N. The first option
       !  will result in urbanisation being a C-source. The second option would make 
       !  urbanisation only a carbon source when biomass is harvested.Given the time of the
       !  year that land cover change is accounted for, urbanizing a cropland would be
       !  carbon neutral under the second option. In the absence of data and considering
       !  that ORCHIDEE has no urban PFTs and therefore does not account for the C uptake
       !  by urban vegetation, the second options seems to be an acceptable compromise
       !  for the moment. The other case, i.e., non_bio -> bio is discussed in 5.2 Land
       !  cover gains within a pixel
       IF (change_nobio(ipts) .GT. 10*EPSILON(un)) THEN
          
          ! Above it was tested whether loss_gains + change_nobio cancel each
          ! other out. If not, there is a problem with the LCC files. Next the
          ! data from the files are a bit cleaned/adjusted to avoid too many
          ! small changes. Small changes themselves are not a real problem but
          ! it avoids ending up with too many PFTs with very small shares.
          ! If change_nobio .GT. min_stomate we convert vegetation to no_bio.
          ! Below the model will bury all the carbon and nitrogen. ORCHIDEE treats
          ! this pool as inert, i.e., nothing is done with that variable (yet).
          ! When PFTs were removed their soil and litter carbon was deposited on the 
          ! bank. We assume that the no_bio fraction is established on all the
          ! PFTs that contributed to the bank. For this reasons we can simply use the
          ! mean soil and litter carbon on the bank (gC m-2) to initialize the soil
          ! and litter carbon of the newly established PFT's. 
          
          ! Store all the details of the different components to be able to 
          ! deal with decreases in frac_nobio. Note that we are weighting
          ! by the change_nobio. The initial farc_nobio is considered to be virgin
          ! ( = zero C and N) soil.
          burried_litter(ipts,:,:,:) = burried_litter(ipts,:,:,:) + &
               litter_bank(:,:,:) * change_nobio(ipts)
          burried_fresh_ltr(ipts,:,:) = burried_fresh_ltr(ipts,:,:) + &
               fresh_ltr_bank(:,:) * change_nobio(ipts)
          burried_fresh_som(ipts,:,:) = burried_fresh_som(ipts,:,:) + &
               fresh_som_bank(:,:) * change_nobio(ipts)
          burried_bact(ipts) = burried_bact(ipts) + &
               bact_bank * change_nobio(ipts)
          burried_min_nitro(ipts,:) = burried_min_nitro(ipts,:) + &
               soil_n_min_bank(:) * change_nobio(ipts)

          ! Temporary aggreagated variable to simplify subsequent calculations
          burried(ipts,icarbon) = (SUM(SUM(litter_bank(:,:,icarbon),1)) + &
               SUM(fresh_ltr_bank(:,icarbon)) + SUM(fresh_som_bank(:,icarbon)) + &
               bact_bank) * change_nobio(ipts)
          
          ! NOTE: until now bact has been a diagnostic carbon pool. The carbon
          ! used in this pool is never taken from the litter. Hence it should
          ! NOT be accounted in the calculation of the nbp in stomate_lpj.f90
          ! Before passing burried to stomate_lpj.f90 the carbon from the bacteria
          ! will be removed again. When bact becomes a regular pool with C and N
          ! (and thus a prognostic variable) it should be accounted for in the nbp
          ! calculation, the follow line of code should be removed, and burried
          ! should be passed rather than burried_adj.
          burried_adj(ipts,icarbon) = (SUM(SUM(litter_bank(:,:,icarbon),1)) + &
               SUM(fresh_ltr_bank(:,icarbon)) + SUM(fresh_som_bank(:,icarbon))) * &
               change_nobio(ipts)

          ! Temporary aggreagated variable to simplify subsequent calculations
          burried(ipts,initrogen) = (SUM(SUM(litter_bank(:,:,initrogen),1)) + &
               SUM(fresh_ltr_bank(:,initrogen)) + SUM(fresh_som_bank(:,initrogen)) + &
               SUM(soil_n_min_bank(:))) * change_nobio(ipts)
          
          ! Note that soil_n_min is a real n-pool that needs to
          ! be accounted for in the nbp. Hence, burried_adj and burried are 
          ! identical for nitrogen
          burried_adj(ipts,initrogen) = burried(ipts,initrogen)

          IF ( ok_soil_carbon_discretization ) THEN
    
             ! Store all the details of the different components
             burried_deepSOM_a(ipts,:,:) = burried_deepSOM_a(ipts,:,:) + &
                  som_vertres_bank(:,iactive,:) * change_nobio(ipts)
             burried_deepSOM_s(ipts,:,:) = burried_deepSOM_s(ipts,:,:) + &
                  som_vertres_bank(:,islow,:) * change_nobio(ipts)
             burried_deepSOM_p(ipts,:,:) = burried_deepSOM_p(ipts,:,:) + &
                  som_vertres_bank(:,ipassive,:) * change_nobio(ipts)

             ! Temporary aggregated variable
             burried(ipts,:) = burried(ipts,:) + &
                  (SUM(som_vertres_bank(:,iactive,:),1) + &
                  SUM(som_vertres_bank(:,islow,:),1) + &
                  SUM(som_vertres_bank(:,ipassive,:),1)) * change_nobio(ipts)
             burried_adj(ipts,:) = burried_adj(ipts,:) + &
                  (SUM(som_vertres_bank(:,iactive,:),1) + &
                  SUM(som_vertres_bank(:,islow,:),1) + &
                  SUM(som_vertres_bank(:,ipassive,:),1)) * change_nobio(ipts)

          ELSE

             ! Store all the details of the different components
             burried_som(ipts,:,:) = burried_som(ipts,:,:) + &
                  soil_bank(:,:) * change_nobio(ipts)

             ! Temporary aggregated variable
             burried(ipts,:) = burried(ipts,:) + &
                  SUM(soil_bank(:,:),1) * change_nobio(ipts)
             burried_adj(ipts,:) = burried_adj(ipts,:) + &
                  SUM(soil_bank(:,:),1) * change_nobio(ipts)

          ENDIF

       ELSE

          ! +++CHECK+++
          ! The loss of frac_nobio should come with the deposition of the burried
          ! carbon on the bank. Needs to be developed. See also earliers comments
          ! on this issue in the code.
          ! +++++++++++

       END IF

       !! 5.2 Land cover gains within a pixel
       !  The consecutive loops over nvm cannot be merged because we first
       !  need to deposit all the available pools on the bank before we
       !  can withdraw the deposits to initialize new PFTs with the weighted
       !  mean of the deposits.
       !  The following block of code deals with one vegetation type being
       !  converted into another as well as with no_bio becoming bio. De-urbanization 
       !  is an example of such a LCC but it is highly unlikely to be substantial. 
       !  More likely are the examples of glacier retreats and lake drying, 
       !  for example the Aral sea in which non biological land uses (in the ORCHIDEE 
       !  terminology) would be converted in biological land use. Current land cover 
       !  change maps are unlikley to account for such changes because their resolution
       !  appears to be modest at best. If the land cover change maps suggest
       !  such a change (= increase in the bio fraction in the pixel) anyway, 
       !  new PFTs will be established further in the code.

       !! 5.2.1 There is no land cover change or there is a gain      

       DO ivm = 1,nvm

          ! There is no land cover change
          IF (loss_gain(ipts,ivm).EQ.zero) CYCLE

          ! Debug
          IF (err_act.GT.1) THEN
             ! Initialize intermediate mass conservation checks
             test_mass(:,istart) = zero
          ENDIF
          !-

          IF ( loss_gain(ipts,ivm) .LE. -EPSILON(un) ) THEN
             
             ! There is no land cover change or there is a loss. In case there 
             ! is no land cover change nothing should be done. The case in which 
             ! there is a loss is already dealt with above - see comment 
             ! below §5.2 to learn why the case of a loss is not integrated in 
             ! this IF-statement. Basically we need to put all the lost C/N on 
             ! the bank first before we can redistribute it. 

          ! Increase in the extent of a PFT
          ELSEIF ( loss_gain(ipts,ivm) .GT. 10*EPSILON(un) ) THEN

             !! 5.2.1 Check wether this is the youngest age class
             !  If the redistribution of the LCC map in age classes is correctly 
             !  implemented a new age class should always be established in the 
             !  youngest age class which should be the first of the group. This 
             !  implies that in this case ::ivm is the first of the group and 
             !  ivm + nagec - 1 is the last age class of the group
             IF (is_tree(ivm)) THEN

                igroup = INT(agec_group(ivm))
                nagec_temp = nagec_pft(igroup)

                IF (agec_group(ivm+nagec_temp-1) .NE. igroup) THEN
                   
                   ! ::ivm and ivm+nagec-1 are in a different ::igroup 
                   ! apparently ::ivm is not the youngest age class
                   WRITE(numout,*) 'Error: problem with the redistribution of the LCC maps'
                   WRITE(numout,*) 'ipts, ivm:  ',ipts, ivm
                   WRITE(numout,*) 'veget_max, loss_gain: ',veget_max(ipts,ivm),&
                        loss_gain(ipts,ivm)
                   WRITE(numout,*) 'nagec_temp, igroup:  ',nagec_temp,igroup
                   WRITE(numout,*) 'agec_group(ivm), agec_group(ivm+nagec_temp-1):  ',&
                        agec_group(ivm), agec_group(ivm+nagec_temp-1)
                   CALL ipslerr_p (3,'land_cover_main','Cannot find proper age classes!','','')

                ENDIF

             ENDIF

             !! 5.2.2 The new PFT does not exist yet. 
             !  PFTs with a cover fraction of less than min_vegfrac should have been weeded out in 
             !  slowproc. When weeded out, the veget_max is set to zero. So if all went well in
             !  slowproc, the test of being LT min_stomate is more or less the same as testing 
             !  whether veget_max is zero. 
             IF (veget_max(ipts,ivm) .LT. min_stomate) THEN

                !! 5.2.2.1 Initialize new and dummy variables
                new_circ_class_biomass(:,:,:,:,:) = zero
                new_circ_class_n(:,:,:) = zero
                new_lm_lastyearmax(:,:) = zero
                new_age(:,:) = zero
                new_KF(:,:) = zero
                new_leaf_frac(:,:,:) = zero
                new_atm_to_bm(:,:,:) = zero
                new_everywhere(:,:) = zero 
                dum_when_growthinit(:,:) = large_value
                dum_npp_longterm(:,:) = zero
                dum_PFTpresent(:,:) = .TRUE.
                IF (loss_gain(ipts,ivm) .GT. min_stomate) THEN
                   ! veget_max is zero which means that the PFT is currently
                   ! empty, but loss_gain is > 0 which means that we have to
                   ! establish this PFT in this time step. So set the status to
                   ! iprescribe
                   plant_status(ipts,ivm) = iprescribe
                ELSE
                   ! veget_max and loss_gain are zero. This PFT does not exist
                   ! set plant_status to inone to enhance consistency within the
                   ! model.
                   plant_status(ipts,ivm) = inone
                END IF
                
                !! 5.2.2.3 Prescribe values for newly established PFTs
                ! Initialize density of individuals, biomass, allocation factors,
                ! leaf age distribution, ... to some reasonable value
                ! ::veget_max is not updated yet, therefore loss_gain is passed in 
                ! the argument list. As ::veget_max is used in the IF-THEN-statements
                ! we won't update it until the end of this routine. If we would
                ! update before that, the same PFT may be subject to conflicting actions
                ! i.e. prescribe a new PFT and then merge it with the other age-classes
                ! in the PFT-group.
             
                ! Debug
                IF (printlev_loc.GE.3) THEN
                   WRITE(numout,*) 'Calling precribe from sapiens_lcchange 1, ', ipts, ivm
                END IF
                !-

                IF (ivm .EQ. ibare_sechiba) THEN
                   
                   ! Reset forestry parameters.
                   age_stand(ipts,ivm) = zero
                   last_cut(ipts,ivm) = zero
                   
                ELSE

                   ! Set flag for writing debug statements in prescribe
                   write_debug = .FALSE.
                   IF (ipts==test_grid.AND.ivm==test_pft) write_debug=.TRUE.

                   ! Calculate the vegetation characteristics for all vegetated PFTs 
                   CALL prescribe (&
                        ivm,                                    loss_gain(ipts,ivm),           dt_days, &
                        dum_PFTpresent(ipts,ivm),               new_everywhere(ipts,ivm),      dum_when_growthinit(ipts,ivm), &
                        new_leaf_frac(ipts,ivm,:),              new_circ_class_n(ipts,ivm,:), &
                        new_circ_class_biomass(ipts,ivm,:,:,:), new_atm_to_bm(ipts,ivm,:),     forest_managed(ipts,ivm), &
                        new_KF(ipts,ivm),                       plant_status(ipts,ivm),        new_age(ipts,ivm), &
                        dum_npp_longterm(ipts,ivm),             new_lm_lastyearmax(ipts,ivm),  longevity_eff_leaf(ipts,ivm), &
                        longevity_eff_sap(ipts,ivm),            longevity_eff_root(ipts,ivm),  k_latosa_adapt(ipts,ivm), &
                        light_tran_to_floor_season(ipts,ivm),   species_change_map(ipts,ivm), &
                        cn_leaf_init_2D(ipts,ivm),              bm_sapl_2D(ipts,ivm,:,:,:),    new_ind(ipts,ivm),&
                        pipe_tune2(ipts,ivm),                   alpha_self_thinning(ipts,ivm), write_debug, &
                        tot_res_target(ipts,ivm,:),             tot_lab_target(ipts,ivm,:),    gpp_week(ipts,ivm), &
                        t2m(ipts))
                   
                   ! Update a couple of forestry parameters. I initialize these to
                   ! one instead of zero since forestry is called after this
                   ! routine, and it looks at last_cut.  This is not a real
                   ! problem since the trees we prescribe generally are the height of
                   ! saplings and not seeds.
                   age_stand(ipts,ivm) = 1
                   last_cut(ipts,ivm) = 1

                END IF

                ! Initialize the vegetation, All zero's for bare soil and
                ! non-zero values for the other PFTs.
                circ_class_biomass(ipts,ivm,:,:,:) = new_circ_class_biomass(ipts,ivm,:,:,:)
                circ_class_n(ipts,ivm,:) = new_circ_class_n(ipts,ivm,:)
                lm_lastyearmax(ipts,ivm) = new_lm_lastyearmax(ipts,ivm)
                croot_longterm(ipts,ivm) = SUM(circ_class_biomass(ipts,ivm,:,iroot,icarbon) * &
                     circ_class_n(ipts,ivm,:))
                age(ipts,ivm) = new_age(ipts,ivm)
                KF(ipts,ivm) = new_KF(ipts,ivm)
                leaf_frac(ipts,ivm,:) = new_leaf_frac(ipts,ivm,:)
                atm_to_bm(ipts,ivm,:) = new_atm_to_bm(ipts,ivm,:)
                everywhere(ipts,ivm) = new_everywhere(ipts,ivm)
                PFTpresent(ipts,ivm) = .TRUE.
                sugar_load(ipts,ivm) = un

                ! Initialize the soil and litter carbon pools. When PFTs were removed
                ! their soil and litter carbon was deposited on the bank. It is assumed
                ! that the newly established PFT are established uniformly on all the
                ! removed PFTs. Hence there is no preferential LCC sequence. For example
                ! we deforest a needle and broadleaved PFT and replace it by a crop and
                ! a grassland the crop and grassland are both established on the former
                ! needle and broadleaved PFT rather than the crop being entirely on the 
                ! former needleleaved PFT and grassland on the former broadleaved PFT.
                ! Further, the LCC-maps should guarantee that the losses equal the gains
                ! which ensures mass preservation. For these reasons we can simply use the
                ! mean soil and litter carbon on the bank (gC m-2) to initialize the soil
                ! and litter carbon of the newly established PFT's. 
                litter(ipts,:,ivm,:,:) = litter_bank(:,:,:)
                lignin_struc(ipts,ivm,:) = struct_ltr_bank(:)
                lignin_wood(ipts,ivm,:) = woody_ltr_bank(:)
                lignin_snag(ipts,ivm,:) = snag_ltr_bank(:)
                bm_to_litter(ipts,ivm,:,:) = fresh_ltr_bank(:,:)
                turnover_daily(ipts,ivm,:,:) = fresh_som_bank(:,:)
                soil_n_min(ipts,ivm,:) = soil_n_min_bank(:)
                bact(ipts,ivm) = bact_bank
                woody_litter_by_cut(ipts,ivm,icut_lcc_wood) = &
                   fresh_ltr_bank(iheartabove,icarbon) + fresh_ltr_bank(isapabove,icarbon)

                IF(printlev_loc .GE. 4 .AND. ipts == test_grid) THEN
                   WRITE(numout,*) 'move bank, ivm',ivm, woody_litter_by_cut(ipts,ivm,icut_lcc_wood)
                   WRITE(numout,*) 'move bank, previous', fresh_woody_ltr_bank
                ENDIF

                IF ( ok_soil_carbon_discretization ) THEN  
                   DO igrn = 1,ngrnd
                      deepSOM_a(ipts,igrn,ivm,:) = som_vertres_bank(igrn,iactive,:) / &
                                                   (zf_soil(igrn)-zf_soil(igrn-1))
                      deepSOM_s(ipts,igrn,ivm,:) = som_vertres_bank(igrn,islow,:) / &
                                                   (zf_soil(igrn)-zf_soil(igrn-1))   
                      deepSOM_p(ipts,igrn,ivm,:) = som_vertres_bank(igrn,ipassive,:) / &
                                                   (zf_soil(igrn)-zf_soil(igrn-1))  
                  ENDDO
                ELSE  
                   som(ipts,:,ivm,:) =  soil_bank(:,:)
                ENDIF

              ! Debug
              ! Is there enough som in soil bank?
              IF (printlev_loc .GE. 4 .AND. ipts .EQ. test_grid .AND. ivm .EQ. test_pft) THEN
                        WRITE(numout,*) 'New PFTs being prescribed. Using soil bank:'  
                        WRITE(numout,*) 'veget_max(ipts,ivm)',veget_max(ipts,ivm)
                        WRITE(numout,*) 'som(ipts,iactive,ivm,initrogen)', &
                             som(ipts,iactive,ivm,initrogen)
                        WRITE(numout,*) 'som(ipts,islow,ivm,initrogen)', &
                             som(ipts,islow,ivm,initrogen)
                        WRITE(numout,*) 'som(ipts,ipassive,ivm,initrogen)', &
                             som(ipts,ipassive,ivm,initrogen)
                        WRITE(numout,*) 'som(ipts,iactive,ivm,ic)', &
                             som(ipts,iactive,ivm,icarbon)
                        WRITE(numout,*) 'som(ipts,islow,ivm,ic)', & 
                             som(ipts,islow,ivm,icarbon)
                        WRITE(numout,*) 'som(ipts,ipassive,ivm,ic)', &
                             som(ipts,ipassive,ivm,icarbon)
                ENDIF
                
             !! 5.2.3 The new PFT does already exist
             !  See comment above. GE min_stomate should result in testing whether veget_max
             !  is greater than min_vegfrac.   
             ELSEIF (veget_max(ipts,ivm) .GE. min_stomate) THEN

                !  We want to make use of stomate_prescribe but it should be noted
                !  that this PFT already contains biomass. Therefore we will create a
                !  temporary PFT without biomass which can be used to establish the
                !  new vegetation within this PFT. For most of the INTENT(inout)
                !  variables of prescribe_prognostic we receive temporary variables 
                !  back this helps us to calculate the weighted mean of the newly 
                !  established vegetation and the vegetation that is already there.
                !  For some other variables we receive dummies as we are simply not
                !  going to use these values but will continue using the values of the
                !  vegetation that is already there.

                !! 5.2.3.1 Initialize new and dummy variables
                new_circ_class_biomass(:,:,:,:,:) = zero
                new_circ_class_n(:,:,:) = zero
                new_lm_lastyearmax(:,:) = zero
                new_age(:,:) = zero
                new_KF(:,:) = zero
                new_leaf_frac(:,:,:) = zero
                new_atm_to_bm(:,:,:) = zero
                new_everywhere(:,:) = zero 
                dum_when_growthinit(:,:) = large_value
                dum_npp_longterm(:,:) = zero
                dum_PFTpresent(:,:) = .TRUE.

                ! Store the current plant_status to re-use it after prescribe
                plant_status_init = plant_status(ipts,ivm)
                ! Set plant_status to iprescribe, else the CALL to prescribe won't work
                plant_status(ipts,ivm) = iprescribe
                
                !! 5.2.3.2 Merge the new and already available vegetation
                !  Veget_max after PFT expansion. As ::veget_max is used in the 
                !  IF-statements we won't update it until the end of this routine. If 
                !  we would update before that, the same PFT may be subject to 
                !  conflicting actions.
                CALL calculate_shares(veget_max(ipts,ivm), loss_gain(ipts,ivm), &
                        litter(ipts,:,ivm,:,:), litter_bank(:,:,:), &
                        share_rec, litter_weight_rec)

                ! Remember to add the moisture from PFT1 that will be present 
                ! if PFT1 is planted with PFT-ivm.
                IF(ABS(total_losses) .GT. min_stomate)THEN

                   ! Here we need to be careful.  loss_gain(ipts,ivm) will always be
                   ! negative, but we are using loss_gain(ipts,1) which may be negative 
                   ! or positive. So we have to take absolute values.
                   wstress_season(ipts,ivm) = un - ( (un - wstress_season(ipts,ivm)) * &
                        share_rec + (moi_bank_season  * &
                        ABS(loss_gain(ipts,1)) / ABS(total_losses) ) * (un - share_rec))
                   wstress_month(ipts,ivm) = un - ( (un - wstress_month(ipts,ivm)) * &
                        share_rec + (moi_bank_month  * &
                        ABS(loss_gain(ipts,1)) / ABS(total_losses) ) * (un - share_rec))
                
                ELSEIF(ABS(total_losses) .GT. min_stomate .AND. err_act .EQ. 2)THEN

                   ! In this case, we expect the nobio fraction to change and we're okay
                   ! with it.  So just leave the long-term water stresses at their
                   ! previous values.
                   WRITE(numout,*) 'WARNING: ipts,ivm: ',ipts,ivm
                   WRITE(numout,*) 'total_losses: ',total_losses
                   CALL ipslerr_p (2,'land_cover_main', &
                        'We have no losses on this pixel, which causes a divide by zero, but ',&
                        'you have indicated that you understand the situation so we continue.',&
                        'Look for Error in out_orchidee file for details')

                ELSEIF(ABS(total_losses) .GT. min_stomate .AND. err_act .EQ. 3)THEN
                
                   ! We are not okay with this situation.  Stop with an error.
                   WRITE(numout,*) 'ERROR: ipts,ivm: ',ipts,ivm
                   WRITE(numout,*) 'total_losses: ',total_losses
                   CALL ipslerr_p (3,'land_cover_main', &
                        'We have no losses on this pixel, which causes a divide by zero, and ',&
                        'you have not indicated that you understand the situation.',&
                        'Look for Error in out_orchidee file for details')
                   
                ENDIF


                !! 5.2.3.3 Initialize the newly established vegetation: 
                !  Initialize the newly established vegetation: density of 
                !  individuals, biomass, allocation factors, leaf age distribution, 
                !  etc to some reasonable value ::veget_max is not updated yet, 
                !  therefore loss_gain is passed in the argument list.
                write_debug = .FALSE.
                IF (ipts==test_grid.AND.ivm==test_pft) write_debug=.TRUE.
                IF (printlev_loc.GE.3) THEN
                   WRITE(numout,*) 'Calling precribe from sapiens_lcchange 2, ', ipts, ivm
                END IF

                IF (ivm == ibare_sechiba) THEN

                   ! Set all vegetation characteristics to zero for bare soil
                   plant_status(ipts,ivm) = inone
                              
                ELSE

                   ! Set flag for writing debug statements in prescribe
                   write_debug = .FALSE.
                   IF (ipts==test_grid.AND.ivm==test_pft) write_debug=.TRUE.
                   ! To avoid mass balance error if ok_dyn_grass_den is activated,
                   ! make new_circ_class_n in newly established
                   ! grassland the same as that in already existed grassland.
                   ! Under the conditions: 1) natural, 2) not tree, 3) flag
                   IF (natural(ivm) .AND. &
                      .NOT. is_tree(ivm) .AND. &
                      ok_dyn_grass_den) THEN
                      new_circ_class_n(ipts,ivm,:) = circ_class_n(ipts,ivm,:)
                   END IF


                   ! Calculate the vegetation characteristics for 
                   ! the new vegetation on expanding PFTs 
                   CALL prescribe (&
                        ivm,                                    loss_gain(ipts,ivm),           dt_days, &
                        dum_PFTpresent(ipts,ivm),               new_everywhere(ipts,ivm),      dum_when_growthinit(ipts,ivm), &
                        new_leaf_frac(ipts,ivm,:),              new_circ_class_n(ipts,ivm,:), &
                        new_circ_class_biomass(ipts,ivm,:,:,:), new_atm_to_bm(ipts,ivm,:),     forest_managed(ipts,ivm), &
                        new_KF(ipts,ivm),                       plant_status(ipts,ivm),        new_age(ipts,ivm), &
                        dum_npp_longterm(ipts,ivm),             new_lm_lastyearmax(ipts,ivm),  longevity_eff_leaf(ipts,ivm), &
                        longevity_eff_sap(ipts,ivm),            longevity_eff_root(ipts,ivm),  k_latosa_adapt(ipts,ivm), &
                        light_tran_to_floor_season(ipts,ivm),   species_change_map(ipts,ivm), &
                        cn_leaf_init_2D(ipts,ivm),              bm_sapl_2D(ipts,ivm,:,:,:),    new_ind(ipts,ivm),&
                        pipe_tune2(ipts,ivm),                   alpha_self_thinning(ipts,ivm), write_debug, &
                        tot_res_target(ipts,ivm,:),             tot_lab_target(ipts,ivm,:),    gpp_week(ipts,ivm), &
                        t2m(ipts))

                   ! Debug
                   IF ((printlev_loc>=4) .AND.(natural(ivm) .AND. &
                      .NOT. is_tree(ivm) .AND. &
                      ok_dyn_grass_den)) THEN
                      WRITE(numout,*) 'After prescribe, ipts,ivm: ',ipts,ivm
                      WRITE(numout,*) 'new density',new_circ_class_n(ipts,ivm,:)
                      WRITE(numout,*) 'current density',circ_class_n(ipts,ivm,:)
                   END IF

        
                   ! The vegetation is supposed to have the same nue as the vegetation
                   ! that was already present. Nothing should be done to nue. The call
                   ! to prescribe will set the plant status to idormant for evergreens. If the
                   ! pft is already present this is wrong. If we have leaves and roots the
                   ! plant status should be icanopy. Here we give back the initial plant_status
                   ! to the old+new pft
                   plant_status(ipts,ivm) = plant_status_init

                END IF

                IF (is_tree(ivm)) THEN

                   !! 5.2.3.3.1 Merge biomass of new and already available trees
                   ! Copy the number of individuals and the biomass of the established 
                   ! vegetation to a temporary variable. In the subroutine ::circ_class_n
                   ! and ::circ_class_biomass will be overwritten with the characteristics
                   ! of the merged vegetation
                   est_circ_class_n(:) = circ_class_n(ipts,ivm,:)
                   est_circ_class_biomass(:,:,:) = circ_class_biomass(ipts,ivm,:,:,:)
                   tmp_circ_class_n(:) = new_circ_class_n(ipts,ivm,:)
                   tmp_circ_class_biomass(:,:,:) = new_circ_class_biomass(ipts,ivm,:,:,:)

                   ! Merge the biomass. When there is an expansion in land cover
                   ! the new and established vegetation have the same PFT
                   ! number, hence the two ivm's in the argument list.        
                   CALL merge_biomass_pfts(npts, share_rec, circ_class_n, &
                        est_circ_class_n, tmp_circ_class_n, circ_class_biomass, &
                        est_circ_class_biomass, &
                        tmp_circ_class_biomass, &
                        ipts, ivm, ivm, forest_managed, bm_to_litter)
                
                ELSE

                   !! 5.2.3.3.2 Merge biomass of new and already available grass or cropland
                   ! Note that above the is_tree flag is used which implies that PFT 1 is
                   ! also dealt with in this ELSE-statement.
                   ! Dilute the existing biomass pools and number of individuals
                   IF (MIN(new_circ_class_n(ipts,ivm,1),circ_class_n(ipts,ivm,1)).GT.zero) THEN
                      ! Both circ_class_n have a value .GT.zero. Hence we can calculate
                      ! the weighted average of circ_class_n
                      circ_class_biomass(ipts,ivm,1,:,:) = share_rec * &
                           circ_class_biomass(ipts,ivm,1,:,:) + &
                           (un - share_rec) * new_circ_class_biomass(ipts,ivm,1,:,:)
                      circ_class_n(ipts,ivm,1) = share_rec * circ_class_n(ipts,ivm,1)+ &
                           (un - share_rec) * new_circ_class_n(ipts,ivm,1)
                   ELSE
                      ! One of the circ_class_n is zero so we can calculate the weighted mean of the
                      ! biomass but we have to keep the non-zero circ_class_n. If we take the mean of
                      ! circ_class_n we decrease the number of individuals we basically decrease the 
                      ! total biomass and just loose C and N.
                      circ_class_biomass(ipts,ivm,1,:,:) = share_rec * &
                           circ_class_biomass(ipts,ivm,1,:,:) + &
                           (un - share_rec) * new_circ_class_biomass(ipts,ivm,1,:,:)
                      circ_class_n(ipts,ivm,1) = MAX(circ_class_n(ipts,ivm,1),new_circ_class_n(ipts,ivm,1))
                   ENDIF

                ENDIF ! is_tree

                !! 5.2.3.4 Calculate the PFT characteristics of the merged PFT
                !  Take the weighted mean of the existing vegetation and the new 
                !  vegetation in this PFT. Note that co2_to_bm is in gC. m-2 dt-1 
                !  so we should also take the weighted mean (rather than sum if
                !  this where absolute values).
                lm_lastyearmax(ipts,ivm) = share_rec * lm_lastyearmax(ipts,ivm) + &
                     (un - share_rec) * new_lm_lastyearmax(ipts,ivm)
                croot_longterm(ipts,ivm) = share_rec * croot_longterm(ipts,ivm) + &
                     (un - share_rec) * SUM(circ_class_biomass(ipts,ivm,:,iroot,icarbon) * &
                     circ_class_n(ipts,ivm,:))
                age(ipts,ivm) = share_rec * age(ipts,ivm) + &
                     (un - share_rec) * new_age(ipts,ivm)
                KF(ipts,ivm) = share_rec * KF(ipts,ivm) + &
                     (un - share_rec) * new_KF(ipts,ivm)
                leaf_frac(ipts,ivm,:) = share_rec * leaf_frac(ipts,ivm,:) + &
                     (un - share_rec) * new_leaf_frac(ipts,ivm,:)
                atm_to_bm(ipts,ivm,:) = share_rec * atm_to_bm(ipts,ivm,:) + &
                     (un - share_rec) * new_atm_to_bm(ipts,ivm,:)
              
                ! sugar load of the new PFT is 1 (=un)
                sugar_load(ipts,ivm) = share_rec * sugar_load(ipts,ivm) + &
                     (un - share_rec) * un

                ! Everywhere deals with the migration of vegetation. Copy the
                ! status of the most migrated vegetation for the whole PFT
                everywhere(ipts,ivm) = MAX(everywhere(ipts,ivm), new_everywhere(ipts,ivm))

                ! The new soil&litter pools are the weighted mean of the newly 
                ! established vegetation for that PFT (as deposited on the bank) 
                ! and the vegetation that already exists in that PFT.
                ! NOTE that in case veget_max of PFT1 (bare soil) is increasing
                ! we will also move litter to the bare soil. stomate_litter.f90
                ! was adjusted such that decomposition also happens in PFT
                litter(ipts,:,ivm,:,:) = share_rec * litter(ipts,:,ivm,:,:) + &
                        (un - share_rec) * litter_bank(:,:,:)
                lignin_struc(ipts,ivm,:) = litter_weight_rec(istructural,:) * &
                     lignin_struc(ipts,ivm,:) + &
                     (un - litter_weight_rec(istructural,:)) * struct_ltr_bank(:)
                lignin_wood(ipts,ivm,:) = litter_weight_rec(iwoody,:) * &
                     lignin_wood(ipts,ivm,:) + &
                     (un - litter_weight_rec(iwoody,:)) * woody_ltr_bank(:)
                lignin_snag(ipts,ivm,:) = litter_weight_rec(isnag,:) * &
                     lignin_snag(ipts,ivm,:) + &
                     (un - litter_weight_rec(isnag,:)) * snag_ltr_bank(:)
                bm_to_litter(ipts,ivm,:,:) = share_rec * bm_to_litter(ipts,ivm,:,:) + & 
                     (un - share_rec) * fresh_ltr_bank(:,:)
                turnover_daily(ipts,ivm,:,:) = share_rec * turnover_daily(ipts,ivm,:,:) + &
                     (un - share_rec) * fresh_som_bank(:,:)
                soil_n_min(ipts,ivm,:) = share_rec * soil_n_min(ipts,ivm,:) + & 
                     (un - share_rec) * soil_n_min_bank(:)
                bact(ipts,ivm) = share_rec * bact(ipts,ivm) + &
                     (un - share_rec) * bact_bank
                !woody_litter_by_cut(ipts,ivm,icut_lcc_wood) = &
                !     share_rec *  woody_litter_by_cut(ipts,ivm,icut_lcc_wood) + &
                !     (un - share_rec) * fresh_woody_ltr_bank
                woody_litter_by_cut(ipts,ivm,icut_lcc_wood) = &
                     share_rec * woody_litter_by_cut(ipts,ivm,icut_lcc_wood) + &
                     (un - share_rec) * (fresh_ltr_bank(iheartabove,icarbon) + & 
                                              fresh_ltr_bank(isapabove,icarbon) )
                IF ( ok_soil_carbon_discretization ) THEN  
                   DO igrn= 1, ngrnd 
                      deepSOM_a(ipts,igrn,ivm,:) = share_rec * deepSOM_a(ipts,igrn,ivm,:) + &
                           (un - share_rec) * som_vertres_bank(igrn,iactive,:) / &
                           (zf_soil(igrn)-zf_soil(igrn-1))
                      deepSOM_s(ipts,igrn,ivm,:) = share_rec * deepSOM_s(ipts,igrn,ivm,:) + &
                           (un - share_rec) * som_vertres_bank(igrn,islow,:) / &
                           (zf_soil(igrn)-zf_soil(igrn-1))
                      deepSOM_p(ipts,igrn,ivm,:) = share_rec * deepSOM_p(ipts,igrn,ivm,:) + &
                           (un - share_rec) * som_vertres_bank(igrn,ipassive,:) / &
                           (zf_soil(igrn)-zf_soil(igrn-1))
                   ENDDO
                ELSE  
                   som(ipts,:,ivm,:) =  share_rec * som(ipts,:,ivm,:) + &
                        (un - share_rec) * soil_bank(:,:)
                ENDIF

                IF (ok_hydrol_arch .AND. ok_vessel_mortality) THEN
                   
                   CALL ipslerr_p (3,'land_cover_change',&
                        'you are the first to use vessel mortality with age classes', &
                        'some code has been proposed but has not yet been tested',&
                        'test before using!')
                   
                   ! Vessel mortality will be stopped when the biomass is moved
                   ! from one to another PFT because the variables are not suited
                   ! to take weighted means. The drought may continue at the next
                   ! time step.
                   kill_vessels(ipts,ivm,:) = .FALSE.
                   vessel_loss_previous(ipts,ivm,:) = zero
                   biomass_init_drought(ipts,ivm,:,:,:) = &
                        circ_class_biomass(ipts,ivm,:,:,:)
                                       
                END IF
                
                ! AR6 output
                ! SL did not see a difference between the calculation of bm_to_litter
                ! and tree_bm_to_litter.
                tree_bm_to_litter(ipts,ivm,:,:) = tree_bm_to_litter(ipts,ivm,:,:)+ &
                     bm_to_litter(ipts,ivm,:,:)

                ! These next two variables are integers, so we use CEILING to make sure
                ! they stay integers.  FLOOR or INT would also be fine.  The use of
                ! CEILING is arbitrary.
                age_stand(ipts,ivm)  = CEILING(share_rec * age_stand(ipts,ivm) + &
                     (un - share_rec) * 1)
                last_cut(ipts,ivm)  = CEILING(share_rec * last_cut(ipts,ivm) + &
                     (un - share_rec) * 1)

                ! Phenology is mostly climate driven so it is fair to assume that the mean
                ! start of the growing season will be rather similar between the old and
                ! merged PFT. Added this line of code to indicate that mean_start_gs should
                ! be considered in land cover change although nothing is really done here.
                mean_start_gs(ipts,ivm) = mean_start_gs(ipts,ivm) 

             ELSE
               
                WRITE(numout,*) 'Error 1: logic of the IF-statements in land_cover_change'
   
             ENDIF ! New PFT

          ELSEIF(ABS(loss_gain(ipts,ivm)) .LE. 10*EPSILON(un))THEN

             ! Do nothing.  This is okay, as far as I can tell.
             ! This seems to be the double precision limit.

          ELSE

             WRITE(numout,*) 'ipts, ivm: ',ipts,ivm
             WRITE(numout,*) 'loss_gain: ',loss_gain(ipts,ivm)
             CALL ipslerr_p(3, 'Logic of IF-statement in land_cover_change',&
                  'Second IF-statement where this could happen','','')

          ENDIF ! Land cover change?

       ENDDO  ! Loop over # PFT's

       !! 5.2.4 CLASS-0 RESTITUTION -- replant the fraction handed back to class 1
       !  Guillaume M. -- loss_gain carries ONE scalar per slot, so the conveyor nets feeding
       !  against restitution and the restitution loses its identity as a GAIN. Without that
       !  event class 1 is never restocked -- `establish` needs a slot empty of biomass and
       !  `recruite` is barred on most of the area -- and its density collapses while its
       !  QMD climbs.
       !
       !  Guillaume M. -- /!\ The area does NOT move here: the net was already applied to
       !  veget_max. Only a FRACTION phi of the slot is replanted, so NO reservoir is touched
       !  -- litter and soil are intensive (per m2). Reusing 5.2.3, which handles a real area
       !  transfer, would withdraw from the bank without a counterpart. The weight is
       !  therefore share = 1 - phi, not the calculate_shares ratio.
       IF (ok_class0_replant .AND. ALLOCATED(class0_restitute)) THEN
          DO c0r_ivm = 1,nvm
             IF (.NOT. is_tree(c0r_ivm)) CYCLE
             IF (class0_restitute(ipts,c0r_ivm) .LE. min_stomate) CYCLE
             IF (veget_max(ipts,c0r_ivm) .LE. min_vegfrac) CYCLE

             c0r_phi = MIN(un, class0_restitute(ipts,c0r_ivm) / veget_max(ipts,c0r_ivm))
             IF (c0r_phi .LE. min_stomate) CYCLE
             c0r_share = un - c0r_phi

             ! Guillaume M. -- Temporaries, exactly as in 5.2.3: prescribe refuses to
             ! establish on a slot that already carries biomass, so it is handed an EMPTY one.
             new_circ_class_biomass(ipts,c0r_ivm,:,:,:) = zero
             new_circ_class_n(ipts,c0r_ivm,:)           = zero
             new_lm_lastyearmax(ipts,c0r_ivm)           = zero
             new_age(ipts,c0r_ivm)                      = zero
             new_KF(ipts,c0r_ivm)                       = zero
             new_leaf_frac(ipts,c0r_ivm,:)              = zero
             new_atm_to_bm(ipts,c0r_ivm,:)              = zero
             new_everywhere(ipts,c0r_ivm)               = zero
             dum_when_growthinit(ipts,c0r_ivm)          = large_value
             dum_npp_longterm(ipts,c0r_ivm)             = zero
             dum_PFTpresent(ipts,c0r_ivm)               = .TRUE.

             plant_status_init = plant_status(ipts,c0r_ivm)
             plant_status(ipts,c0r_ivm) = iprescribe

             ! Guillaume M. -- /!\ The area handed to prescribe is the RESTITUTED fraction.
             ! atm_to_bm comes out per m2 of PFT (per-tree biomass * density / dt, with no
             ! area scaling), so the blend by (1 - share) below does NOT double-count.
             CALL prescribe (&
                  c0r_ivm,                       class0_restitute(ipts,c0r_ivm),     dt_days, &
                  dum_PFTpresent(ipts,c0r_ivm),  new_everywhere(ipts,c0r_ivm),       &
                  dum_when_growthinit(ipts,c0r_ivm), &
                  new_leaf_frac(ipts,c0r_ivm,:), new_circ_class_n(ipts,c0r_ivm,:), &
                  new_circ_class_biomass(ipts,c0r_ivm,:,:,:), new_atm_to_bm(ipts,c0r_ivm,:), &
                  forest_managed(ipts,c0r_ivm), &
                  new_KF(ipts,c0r_ivm),          plant_status(ipts,c0r_ivm),         &
                  new_age(ipts,c0r_ivm), &
                  dum_npp_longterm(ipts,c0r_ivm), new_lm_lastyearmax(ipts,c0r_ivm),  &
                  longevity_eff_leaf(ipts,c0r_ivm), &
                  longevity_eff_sap(ipts,c0r_ivm), longevity_eff_root(ipts,c0r_ivm), &
                  k_latosa_adapt(ipts,c0r_ivm), &
                  light_tran_to_floor_season(ipts,c0r_ivm), species_change_map(ipts,c0r_ivm), &
                  cn_leaf_init_2D(ipts,c0r_ivm), bm_sapl_2D(ipts,c0r_ivm,:,:,:),     &
                  new_ind(ipts,c0r_ivm), &
                  pipe_tune2(ipts,c0r_ivm),      alpha_self_thinning(ipts,c0r_ivm),  write_debug, &
                  tot_res_target(ipts,c0r_ivm,:), tot_lab_target(ipts,c0r_ivm,:),    &
                  gpp_week(ipts,c0r_ivm), &
                  t2m(ipts))

             plant_status(ipts,c0r_ivm) = plant_status_init

             ! Guillaume M. -- Probe: does this block CREATE the litter NaN or PROPAGATE it?
             ! isnan is blind to +-Inf, yet 0 * Inf is the only product that turns a finite
             ! per-tree biomass into a NaN in bm_to_litter below. So test NON-FINITENESS on
             ! both sides of the accumulation, per circumference class, and report density
             ! next to it: a non-finite biomass on zero density is wiped by the ELSE branch.
             c0r_bad_in = .NOT. (c0r_phi .GT. -HUGE(un) .AND. c0r_phi .LT. HUGE(un))
             DO c0r_icir = 1,ncirc
                IF ( .NOT. (circ_class_n(ipts,c0r_ivm,c0r_icir) .GT. -HUGE(un) .AND. &
                            circ_class_n(ipts,c0r_ivm,c0r_icir) .LT. HUGE(un)) .OR. &
                     ANY( .NOT. (circ_class_biomass(ipts,c0r_ivm,c0r_icir,:,:) .GT. -HUGE(un) .AND. &
                                 circ_class_biomass(ipts,c0r_ivm,c0r_icir,:,:) .LT. HUGE(un)) ) .OR. &
                     .NOT. (new_circ_class_n(ipts,c0r_ivm,c0r_icir) .GT. -HUGE(un) .AND. &
                            new_circ_class_n(ipts,c0r_ivm,c0r_icir) .LT. HUGE(un)) .OR. &
                     ANY( .NOT. (new_circ_class_biomass(ipts,c0r_ivm,c0r_icir,:,:) .GT. -HUGE(un) .AND. &
                                 new_circ_class_biomass(ipts,c0r_ivm,c0r_icir,:,:) .LT. HUGE(un)) ) ) THEN
                   c0r_bad_in = .TRUE.
                   WRITE(numout,*) '[C0NAN] ENTREE ipts',ipts,' ipft',c0r_ivm,' icir',c0r_icir, &
                        ' phi',c0r_phi,' vmax',veget_max(ipts,c0r_ivm), &
                        ' | n_old',circ_class_n(ipts,c0r_ivm,c0r_icir), &
                        ' bm_old_max',MAXVAL(circ_class_biomass(ipts,c0r_ivm,c0r_icir,:,:)), &
                        ' | n_new',new_circ_class_n(ipts,c0r_ivm,c0r_icir), &
                        ' bm_new_max',MAXVAL(new_circ_class_biomass(ipts,c0r_ivm,c0r_icir,:,:))
                ENDIF
             ENDDO
             IF (c0r_bad_in) CALL flush(numout)

             ! Guillaume M. -- /!\ THE REPLACED FRACTION MUST GO SOMEWHERE. The blend below
             ! gives the standing stand the weight (1 - phi) only, so the phi share of the
             ! old stand vanishes unless it is routed; destroying it breaks the mass balance.
             ! It goes to litter, the conservative treatment already used for orphan biomass.
             ! /!\ Part of it should physically go to harvest products; still to be settled.
             DO c0r_icir = 1,ncirc
                bm_to_litter(ipts,c0r_ivm,:,:) = bm_to_litter(ipts,c0r_ivm,:,:) + &
                     c0r_phi * circ_class_n(ipts,c0r_ivm,c0r_icir) * &
                     circ_class_biomass(ipts,c0r_ivm,c0r_icir,:,:)
             ENDDO

             ! Guillaume M. -- Verdict of the probe above. A non-finite litter reached with
             ! every input finite means THIS accumulation created it; otherwise it only
             ! carried in what it was handed. The distinction decides where the fix goes.
             c0r_bad_out = ANY( .NOT. (bm_to_litter(ipts,c0r_ivm,:,:) .GT. -HUGE(un) .AND. &
                                       bm_to_litter(ipts,c0r_ivm,:,:) .LT. HUGE(un)) )
             IF (c0r_bad_out) THEN
                WRITE(numout,*) '[C0NAN] SORTIE ipts',ipts,' ipft',c0r_ivm, &
                     ' entree_non_finie',c0r_bad_in,' verdict ', &
                     MERGE('CREE   ','PROPAGE',.NOT. c0r_bad_in)
                CALL flush(numout)
             ENDIF

             ! Guillaume M. -- Blend per circumference class. Biomass is PER TREE, so n and
             ! n*b are blended separately: n_mix*b_mix = share*n_old*b_old +
             ! (1-share)*n_new*b_new conserves carbon EXACTLY, a direct blend of b does not.
             DO c0r_icir = 1,ncirc
                c0r_nmix = c0r_share * circ_class_n(ipts,c0r_ivm,c0r_icir) + &
                     (un - c0r_share) * new_circ_class_n(ipts,c0r_ivm,c0r_icir)
                IF (c0r_nmix .GT. min_stomate) THEN
                   circ_class_biomass(ipts,c0r_ivm,c0r_icir,:,:) = &
                        ( c0r_share * circ_class_n(ipts,c0r_ivm,c0r_icir) * &
                          circ_class_biomass(ipts,c0r_ivm,c0r_icir,:,:) + &
                          (un - c0r_share) * new_circ_class_n(ipts,c0r_ivm,c0r_icir) * &
                          new_circ_class_biomass(ipts,c0r_ivm,c0r_icir,:,:) ) / c0r_nmix
                   circ_class_n(ipts,c0r_ivm,c0r_icir) = c0r_nmix
                ELSE
                   ! Guillaume M. -- Nothing on either side: never leave orphan per-tree
                   ! biomass sitting on a zero density, that path destroys mass.
                   circ_class_n(ipts,c0r_ivm,c0r_icir) = zero
                   circ_class_biomass(ipts,c0r_ivm,c0r_icir,:,:) = zero
                ENDIF
             ENDDO

             ! Guillaume M. -- Same weighting as 5.2.3 for the state variables that follow
             ! the stand.
             lm_lastyearmax(ipts,c0r_ivm) = c0r_share * lm_lastyearmax(ipts,c0r_ivm) + &
                  (un - c0r_share) * new_lm_lastyearmax(ipts,c0r_ivm)
             age(ipts,c0r_ivm) = c0r_share * age(ipts,c0r_ivm) + &
                  (un - c0r_share) * new_age(ipts,c0r_ivm)
             KF(ipts,c0r_ivm) = c0r_share * KF(ipts,c0r_ivm) + &
                  (un - c0r_share) * new_KF(ipts,c0r_ivm)
             leaf_frac(ipts,c0r_ivm,:) = c0r_share * leaf_frac(ipts,c0r_ivm,:) + &
                  (un - c0r_share) * new_leaf_frac(ipts,c0r_ivm,:)
             ! Guillaume M. -- /!\ atm_to_bm is a cumulated FLUX, not a state: it is ADDED,
             ! never blended. 5.2.3 weights it because the slot area GROWS and the flux is
             ! per m2; here the area does not move, so diluting by (1 - phi) would erase part
             ! of the flux already booked. The seedlings create phi * new_atm_to_bm per m2.
             atm_to_bm(ipts,c0r_ivm,:) = atm_to_bm(ipts,c0r_ivm,:) + &
                  c0r_phi * new_atm_to_bm(ipts,c0r_ivm,:)

             ! Guillaume M. -- STAND_AGE: the replanted share carries a freshly established
             ! stand. Blend by AREA with age_stand_estab (mirror of the progressive-harvest
             ! path): the area arrives with saplings, their biomass comes later, so a
             ! biomass weight would be inert and the indicator kept climbing without bound.
             age_stand_bm(ipts,c0r_ivm) = c0r_share * age_stand_bm(ipts,c0r_ivm) + &
                  (un - c0r_share) * age_stand_estab
             age_stand_area(ipts,c0r_ivm) = c0r_share * age_stand_area(ipts,c0r_ivm) + &
                  (un - c0r_share) * age_stand_estab

             ! Guillaume M. -- Consumed: a second read would replant the same area twice.
             class0_restitute(ipts,c0r_ivm) = zero
          ENDDO
       ENDIF

       ! Additional mass balance checks
       IF (err_act.GT.1) THEN

          ! This part of an additional mass balance check for variables
          ! that pass through the bank. Litter, som, bact and soil_n_min
          ! may be moved from one PFT to another but they are not mixed.
          ! For each of the variables their value at the start of this
          ! subroutine should equal their value at the end
          CALL fin_mass_balance(test_bank, veget_max, litter, som, &
               circ_class_biomass, circ_class_n, turnover_daily, &
               bm_to_litter, bact, soil_n_min, change_nobio, loss_gain, &
               litter_bank, soil_bank, bact_bank, &
               soil_n_min_bank, deepSOM_a, deepSOM_s, deepSOM_p, som_vertres_bank, &
               ipts, dt_days, zf_soil)

       ENDIF ! err_act

    ENDDO ! Loop over npts


 !! 6. Update the fraction areas
    
    ! Calculate new veget_max values and update variables if PFT was removed
    veget_max(:,:) = veget_max(:,:) + loss_gain(:,:)

    ! Guillaume M. -- An area below min_stomate is not an area: it is the rounding residue of a
    ! slot meant to be emptied here, and a slot with no ground but standing content is what
    ! stops stomate_prescribe on "Unexpected condition" the following year. Snap it to zero and
    ! hand the residue to the largest slot of the pixel, which is above min_vegfrac by
    ! construction, so the pixel area stays exact and no new sliver is created.
    DO sv_i = 1,npts
       sv_res = zero
       DO sv_k = 1,nvm
          IF (veget_max(sv_i,sv_k) > zero .AND. veget_max(sv_i,sv_k) < min_stomate) THEN
             ! Guillaume M. -- Only a BARE residue may be snapped. Trees standing on a slot
             ! that loses its ground would have their reservoirs orphaned: such a slot is
             ! not a rounding artefact but a real inconsistency between the two area maps,
             ! and it must be fixed where it is produced.
             sv_n = SUM(circ_class_n(sv_i,sv_k,:))
             IF (sv_n > min_stomate) THEN
                WRITE(numout,*) 'SLIVER with standing trees, ipts', sv_i, ' ipft', sv_k, &
                     ' veget_max', veget_max(sv_i,sv_k), ' n', sv_n
                CALL flush(numout)
                CALL ipslerr_p(2,'land_cover_change_main', &
                     'a slot below min_stomate still carries trees', &
                     'left untouched - the producer of that area must be fixed','')
             ELSE
                ! Guillaume M. -- Bare residue: safe to remove, but never in silence. The
                ! conveyor is meant to land emptied slots on exact zero, so reaching this
                ! branch says a transfer path is leaving one behind again.
                WRITE(numout,*) 'SLIVER snapped, ipts', sv_i, ' ipft', sv_k, &
                     ' veget_max', veget_max(sv_i,sv_k)
                CALL flush(numout)
                CALL ipslerr_p(2,'land_cover_change_main', &
                     'a bare slot below min_stomate was snapped to zero', &
                     'an area transfer is not landing on exact zero','')
                sv_res = sv_res + veget_max(sv_i,sv_k)
                veget_max(sv_i,sv_k) = zero
             ENDIF
          ENDIF
       ENDDO
       IF (sv_res > zero) THEN
          sv_kmax = MAXLOC(veget_max(sv_i,:), DIM=1)
          veget_max(sv_i,sv_kmax) = veget_max(sv_i,sv_kmax) + sv_res
       ENDIF
    ENDDO

    ! swap frac_nobio. Strictly speaking this swap could be done in slowproc
    ! because nowhere in stomate frac_nobio should be changed (because stomate
    ! deals with the vegetation not with the non-biological cover fractions).
    ! sapiens_lcchange is the exception to this and to make it explicit that
    ! this is the routine where vegetation and non-biological fractions interact
    ! this swap was implemented here. If the code is expanded with lakes and 
    ! cities, sapiens_lcchange should be taken out of stomate and should be
    ! run after lakes, cities and vegetation has been calculated.
    frac_nobio(:,:) = frac_nobio_new(:,:)

    ! Debug
    IF (printlev_loc>=4.AND.ipts==test_grid) THEN
       ! We don want to have C/N values for PFTs without veget_max. That
       ! would indicate a problem with reading and/or processing the land
       ! cover change maps.
       WRITE(numout,*) 'Checking whether biomass and veget_max are aligned'
       DO ipts = 1,npts
          WRITE(numout,*) 'ivm, veget_max, biomass'
          DO ivm = 1,nvm
             WRITE(numout,'(I2,2F20.10)') ivm, veget_max(ipts,ivm),&
                  SUM(SUM(circ_class_biomass(ipts,ivm,:,:,icarbon),2))
          ENDDO
          WRITE(numout,*) 'ivm, veget_max, bm_to_litter'  
          DO ivm = 1,nvm
             WRITE(numout,'(I2,2F20.10)') ivm, veget_max(ipts,ivm),&
                  SUM(bm_to_litter(ipts,ivm,:,icarbon))
          ENDDO
          WRITE(numout,*) 'ivm, veget_max, turnover_daily'
          DO ivm = 1,nvm
             WRITE(numout,'(I2,2F20.10)') ivm, veget_max(ipts,ivm),&
                  SUM(turnover_daily(ipts,ivm,:,icarbon))
          ENDDO
          WRITE(numout,*) 'ivm, veget_max, litter'
          DO ivm = 1,nvm
             WRITE(numout,'(I2,2F20.10)') ivm, veget_max(ipts,ivm),&
                  SUM(SUM(litter(ipts,:,ivm,:,icarbon),2))
          ENDDO
          IF (ok_soil_carbon_discretization) THEN
             WRITE(numout,*) 'ivm, veget_max, deepSOM'
             DO ivm = 1,nvm
                WRITE(numout,'(I2,2F20.10)') ivm, veget_max(ipts,ivm),&
                     SUM(deepSOM_a(ipts,:,ivm,icarbon)+&
                    deepSOM_s(ipts,:,ivm,icarbon)+deepSOM_p(ipts,:,ivm,icarbon))
             ENDDO
          ELSE
             WRITE(numout,*) 'ivm, veget_max, som'
             DO ivm = 1,nvm
                WRITE(numout,'(I2,2F20.10)') ivm, veget_max(ipts,ivm),&
                     SUM(som(ipts,:,ivm,icarbon))
             ENDDO
          END IF
       ENDDO
    ENDIF
    !-

!! 7. Mass balance closure

    IF (err_act.GT.1) THEN

       ! 7.1 Check surface area
       ! In sapiens_lcchange, both the vegetation cover and the nobio fraction
       ! can change in the pixel. The surface area check should account for 
       ! both changes.
       CALL check_vegetation_area("sapiens_lccchange", npts, veget_max_begin, &
            veget_max,'pixel',change_nobio)

       ! Quality check. At the end of adjust_delta veget_max_new
       ! was calculated as veget_max + loss_gain. Two lines above this line
       ! veget was thus set to veget_max_new. If nothing bad happened
       ! to veget_max, veget_max_new and/or loss_gain in between slowproc
       ! and sapiens_lcchange, veget_max should now equal veget_max_new
       CALL check_variable_swap("sapiens_lccchange", npts, veget_max, &
            veget_max_new)
  
       ! More detailed check at the pft-level of changes in biomass
       ! in the case of disturbance, litter has been moved to the
       !yountest age but it is not considered in the below routine
       CALL check_biomass_change(npts, dt_days, area, circ_class_biomass, &
            circ_class_n, veget_max, loss_gain, atm_to_bm_old, atm_to_bm, &
            harvest_pool, harvest_pool_old, fresh_litter, bm_to_litter_old, &
            biomass_start, biomass_end, 2)

       !! 7.2 Calculate components of the mass balance
       pool_end(:,:,:) = zero
       DO iele = 1,nelements

          ! Biomass pool + bm_to_litter
          DO ipar = 1,nparts
            DO icir=1, ncirc
                pool_end(:,:,iele) = pool_end(:,:,iele) + &
                    circ_class_biomass(:,:,icir,ipar,iele)*&
                    circ_class_n(:,:,icir) *veget_max(:,:)
            ENDDO
            pool_end(:,:,iele) = pool_end(:,:,iele) + &
                (bm_to_litter(:,:,ipar,iele) + turnover_daily(:,:,ipar,iele)) * &
                veget_max(:,:)
          ENDDO

          ! Litter pool (gC m-2) *  (m2 m-2) 
          DO ilit = 1,nlitt
             DO ilev = 1,nlevs
                pool_end(:,:,iele) = pool_end(:,:,iele) + &
                     litter(:,ilit,:,ilev,iele) * veget_max(:,:)
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

          ! C and N burried by urbanisation. Burried has no PFT
          ! dimension any more. Already multiplied with change_nobio
          ! making it easier to use the variable in the nbp calculation
          ! in stomate_lpj.f90
          pool_end(:,1,iele) = pool_end(:,1,iele) + burried(:,iele)

          ! The biomass harvest pool is expressed in gC pixel-1 So, it 
          ! shouldn't be multiplied by veget_max but it should be divided
          ! by area to obtain gC m-2. The difference between 
          ! pool_start and pool_end was accumulated over dt, therefore, 
          ! the units of ipoolchange are gC m-2 dt-1.
          DO ipts = 1,npts
             pool_end(ipts,:,iele) = pool_end(ipts,:,iele) + &
                  SUM(SUM(harvest_pool(ipts,:,:,iele,:),3),2) / &
                  (area(ipts) * contfrac(ipts))
          ENDDO

       ENDDO
       
       !Denitrifying bacteria (only C)
       pool_end(:,:,icarbon) = pool_end(:,:,icarbon) + &
               bact(:,:) * veget_max(:,:)

       ! Nitrogen only
       DO ispec = 1,nnspec
          pool_end(:,:,initrogen) = pool_end(:,:,initrogen) + &
              soil_n_min(:,:,ispec) * veget_max(:,:)
       ENDDO

       !! 7.2 Calculate mass balance
       ! Common processes
       DO iele=1,nelements
          check_intern(:,:,iatm2land,iele) = &
               check_intern(:,:,iatm2land,iele) + &
               atm_to_bm(:,:,iele) * veget_max(:,:) * dt_days
          check_intern(:,:,ipoolchange,iele) = -un * (pool_end(:,:,iele) - &
               pool_start(:,:,iele)) 
       ENDDO

       closure_intern(:,:,:) = zero
       DO imbc = 1,nmbcomp
          DO iele=1,nelements
            DO ivm=1, nvm
             ! Debug
             IF (printlev_loc>=5) THEN              
                IF(iele .EQ. 1) WRITE(numout,"(9999(G12.5,:,','))") &
                    'imbc, ivm, check_intern', test_grid , imbc, &
                    ivm, check_intern(test_grid,ivm,imbc,iele)
             ENDIF
             !-
            ENDDO
             closure_intern(:,:,iele) = closure_intern(:,:,iele) + &
                  check_intern(:,:,imbc,iele)
          ENDDO
       ENDDO

       ! 7.3 Check mass balance closure
       ! It is difficult to find the source of mass balance problems in
       ! sapiens_lcchange because C and N are moved between PFTs and between
       ! compartments. Previously changes in biomass, litter, som, bact and 
       ! soil_n_min and the surface area were checked. It implies that if 
       ! there is an imbalance and none the intermediate tests failed, the 
       ! best variables to start looking are bm_to_litter, turnover_daily and 
       ! burried. Nevertheless, there is no guarantee that the intermediate
       ! checks are 100% correct (although they already passed lots of testing).
       CALL check_mass_balance("sapiens_lccchange", closure_intern, npts, &
            pool_end, pool_start, veget_max, 'pixel')
    
    ENDIF ! err_act.GT.1

    CALL check_area_invariant('exit land_cover_change_main', npts, veget_max, circ_class_n)

    IF (printlev.GE.4) WRITE(numout,*) 'Leaving land_cover_change_main'

  END SUBROUTINE land_cover_change_main

! ================================================================================================================================
!! SUBROUTINE   : age_class_distr
!!
!>\BRIEF        Redistribute biomass, litter, soilcarbon and water across
!!              the age classes
!!
!! DESCRIPTION  : Following growth, the trees from an age class may have become
!! too big to belong to this age class. The biomass, litter, soilcarbon and
!! soil water then need to be moved from one age class to the next age class.
!! Note: also fluxes such as gpp_d, atm_to_bm, resp_maint_d, resp_growth_d, 
!! resp_hetero_d, input_daily, emissions_d and leaching_d are moved to the new
!! age class to enable nbp checks in stomate.f90. For these checks the fluxes
!! are multiplied with veget_max so every time veget_max is changing, the 
!! fluxes should be moved to the correct veget_max.
!!
!! RECENT CHANGE(S) : 
!!
!! MAIN OUTPUT VARIABLE(S) :  
!!
!! REFERENCES   : None
!!
!! FLOWCHART    : 
!! \n
!_ ================================================================================================================================

  
  SUBROUTINE age_class_distr(npts, circ_class_n, circ_class_biomass, &
       veget_max, veget_max_new, wstress_season, &
       lm_lastyearmax, lm_thisyearmax, age, leaf_frac, atm_to_bm, &
       everywhere, litter, som, lignin_struc, &
       lignin_wood, lignin_snag, bm_to_litter, tree_bm_to_litter, &
       turnover_daily,PFTpresent, when_growthinit,&
       forest_managed, KF, plant_status, &
       npp_longterm, croot_longterm, gpp_daily, gpp_year, gpp_decade, leaf_age, &
       gdd_from_growthinit, gdd_midwinter, time_hum_min, &
       hum_min_dormance, gdd_m5_dormance, ncd_dormance, &
       vegstress, humrel, season_drought_legacy, & 
       vegstress_month, vegstress_week, ngd_minus5, &
       resp_maint, resp_growth, npp_daily, &
       rue_longterm, mai, pai, &
       mai_count, previous_wood_volume, vegstress_season, &
       MatrixA, MatrixV, VectorB, VectorU, age_stand, age_stand_bm, age_stand_area, last_cut, &
       k_latosa_adapt, fm_change_map, lpft_replant, cn_leaf_min_season,&
       cn_leaf_init_2D, nstress_season, soil_n_min, p_O2, bact, &
       CN_som_litter_longterm, sugar_load, &
       deepSOM_a, deepSOM_s, deepSOM_p, kill_vessels, vessel_loss_previous, &
       biomass_init_drought, mean_start_gs, resp_hetero, &
       emission_daily, leaching_daily, n_input_daily, &
       gpp_week, resp_maint_week, light_tran_to_floor_season, &
       qsintveg, n_input, Light_Abs_Tot, Light_Tran_Tot, &
       laieff_isotrop, &
       maxvegstress_lastyear, maxvegstress_thisyear, &
       minvegstress_lastyear, minvegstress_thisyear, &
       maxgppweek_lastyear, maxgppweek_thisyear, &
       wstress_month, n_reserve_longterm, n_reserve_balance, & 
       maxfpc_lastyear, maxfpc_thisyear,&
       turnover_longterm, dead_leaves, &
       grow_season_len, us, leaf_age_crit, leaf_classes, &
       co2_fire,litterfuel,count_daylight,zf_soil, &
       gap_area_save, total_ba_init)

    IMPLICIT NONE
    
  !! 0. Variable and parameter declaration 
    
    !! 0.1 Input variables

    INTEGER, INTENT(in)                                :: npts                !! Domain size - number of pixels (unitless)
    INTEGER(i_std), DIMENSION(:,:), INTENT(in)         :: fm_change_map       !! A map which gives the desired FM strategy when
                                                                              !! the PFT will be replanted after a clearcut.
                                                                              !! (1-nvm,unitless)
    LOGICAL, DIMENSION(:,:), INTENT(in)                :: lpft_replant        !! Set to true if a PFT has been clearcut
                                                                              !! and needs to be replaced by another species
    REAL(r_std),DIMENSION(:,:), INTENT(in)             :: cn_leaf_init_2D     !! initial leaf C/N ratio
    REAL(r_std), DIMENSION(0:ngrnd), INTENT(in)        :: zf_soil

    !! 0.2 Output variables


    !! 0.3 Modified variables
    
    LOGICAL, DIMENSION(:,:), INTENT(inout)             :: PFTpresent          !! Tab indicating which PFTs are present in 
                                                                              !! each pixel
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: plant_status        !! Growth and phenological status of the plant
                                                                              !! istatus = Phases defined in constantes
    INTEGER(i_std), DIMENSION (:,:), INTENT(inout)     :: forest_managed      !! forest management flag (is the forest 
                                                                              !! being managed?)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: vegstress_month     !! "Monthly" moisture availability (0 to 1, 
                                                                              !! unitless) 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: vegstress_week      !! "Weekly" moisture availability 
                                                                              !! (0 to 1, unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: ngd_minus5          !! Number of growing days (days), threshold 
                                                                              !! -5 deg C (for phenology)   
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: resp_maint          !! Maintenance respiration  
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: resp_growth         !! Growth respiration  
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: npp_daily           !! Net primary productivity 
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: rue_longterm        !! Longterm radiation use efficiency 
                                                                              !! (??units??) 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: when_growthinit     !! How many days ago was the beginning of 
                                                                              !! the growing season (days)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: KF                  !! Scaling factor to convert sapwood mass
                                                                              !! into leaf mass (m)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: k_latosa_adapt      !! Leaf to sapwood area adapted for long 
                                                                              !! term water stress (m)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: npp_longterm        !! "Long term" mean yearly primary productivity
                                                                              !! @tex $(m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: croot_longterm      !! "long term" root carbon mass  
 	                                                                      !! @tex ($gC m^{-2}) @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: veget_max           !! "maximal" coverage fraction of a PFT on the ground
                                                                              !! May sum to
                                                                              !! less than unity if the pixel has
                                                                              !! nobio area. (unitless, 0-1)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: veget_max_new       !!"maximal" coverage fraction of a PFT on the ground
                                                                              !(with change of land cover)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: wstress_season      !! Water stress factor, based on hum_rel_daily
                                                                              !! (unitless, 0-1)   
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: lm_lastyearmax      !! last year's maximum leaf mass for each PFT 
                                                                              !! @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: lm_thisyearmax      !! this year's maximum leaf mass, for each 
                                                                              !! PFT @tex ($gC m^{-2}$) @endte
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: everywhere          !! is the PFT everywhere in the grid box or 
                                                                              !! very localized (after its introduction) (?)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: age                 !! mean age (years)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: atm_to_bm           !! CO2 and N taken from the atmosphere to get C to create  
                                                                              !! the seedlings @tex (gC.m^{-2}dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: gpp_daily           !! Daily gross primary productivity  
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: gpp_year            !! "annual" GPP @tex ($gC m^{-2} day^{-1}$)@endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: gpp_decade          !! "decadal" GPP @tex ($gC m^{-2} day^{-1}$)@endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: time_hum_min        !! Time elapsed since strongest moisture 
                                                                              !! availability (days) 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: hum_min_dormance    !! minimum moisture during dormance 
                                                                              !! (0-1, unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: gdd_midwinter       !! Growing degree days (K), since midwinter 
                                                                              !! (for phenology) - this is written to the
                                                                              !!  history files 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: gdd_from_growthinit !! growing degree days, since growthinit 
                                                                              !! for crops
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: gdd_m5_dormance     !! Growing degree days (K), threshold -5 deg 
                                                                              !! C (for phenology)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: ncd_dormance        !! Number of chilling days (days), since 
                                                                              !! leaves were lost (for phenology) 
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: circ_class_n        !! Number of individuals in each circ class
                                                                              !! @tex $(ind m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: lignin_struc        !! ratio Lignine/Carbon in structural litter,
                                                                              !! above and below ground
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: lignin_wood         !! ratio Lignine/Carbon in woody litter,
                                                                              !! above and below ground
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: lignin_snag         !! ratio Lignine/Carbon in snag litter,
                                                                              !! above and below ground
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: som                 !! carbon pool: active, slow, or passive 
                                                                              !! @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: deepSOM_a           !! Soil carbon discretized with depth active (g/m**3) 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: deepSOM_s           !! Soil carbon discretized with depth slow (g/m**3) 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: deepSOM_p           !! Soil carbon discretized with depth passive (g/m**3)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: leaf_frac           !! fraction of leaves in leaf age class (unitless;0-1)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: leaf_age            !! Leaf age (days)
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: bm_to_litter        !! Transfer of biomass to litter 
                                                                              !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: tree_bm_to_litter   !!Transfer of tree biomass to litter 
                                                                              !!@tex$(gC m^{-2} dtslow^{-1})$@endtex
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: turnover_daily      !! Transfer of litter to som 
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)   :: circ_class_biomass  !! Biomass components of the model tree within a 
                                                                              !! circumference class
                                                                              !! class @tex $(g C ind^{-1})$ @endtex  
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)   :: litter              !! metabolic and structural litter, above and 
                                                                              !! below ground @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: mai                 !! The mean annual increment
                                                                              !! @tex $(m**3 / m**2 / year)$ @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: pai                 !! The period annual increment
                                                                              !! @tex $(m**3 / m**2 / year)$ @endtex          
    INTEGER(i_std), DIMENSION(:,:),INTENT(inout)       :: mai_count           !! The number of times we've
                                                                              !! calculated the volume increment
                                                                              !! for a stand 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: previous_wood_volume!! The volume of the tree trunks
                                                                              !! in a stand for the previous year.
                                                                              !! @tex $(m**3 / m**2 )$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: vegstress_season  !! Mean growingseason moisture 
                                                                              !! availability (0 to 1, unitless)
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: MatrixA             !! Matrix containing the fluxes between the
                                                                              !! carbon pools per sechiba time step 
                                                                              !! @tex $(gC.m^2.day^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: VectorB             !! Vector containing the litter increase per
                                                                              !! sechiba time step
                                                                              !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: MatrixV             !!
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: VectorU             !!
    INTEGER(i_std), DIMENSION(:,:), INTENT(inout)      :: age_stand           !! Age of the forest stand (years)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: age_stand_area      !! AREA-weighted conserved mean stand age (years) - STAND_AGE
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: age_stand_bm        !! Biomass-weighted conserved mean stand age (years) - STAND_AGE
    INTEGER(i_std), DIMENSION(:,:), INTENT(inout)      :: last_cut            !! Years since last thinning (years)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: cn_leaf_min_season  !! Seasonal min CN ratio of leaves 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: nstress_season      !! N-related seasonal stress (used for allocation) 
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: soil_n_min          !! mineral nitrogen in the soil (gN/m**2)  
                                                                              !! (first index=npts, second index=nvm, third index=nnspec) 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: p_O2                !! partial pressure of oxygen in the soil (hPa)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: bact                !! denitrifier biomass (gC/m**2)
    REAL(r_std), DIMENSION(:,:,:),INTENT(inout)        :: CN_som_litter_longterm !! Longterm CN ratio of litter and som pools (gC/gN)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: sugar_load           !! Relative sugar loading of the labile pool (unitless)
    LOGICAL, DIMENSION(:,:,:), INTENT(inout)           :: kill_vessels         !! Flag to kill vessels at the end of the day following embolism.
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: vessel_loss_previous !! Proportion of conductivity lost due to cavitation, accumulated
                                                                               !!  on the previous day (unitless).
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)   :: biomass_init_drought !! Biomass of heartwood or sapwood before onset of drought
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: mean_start_gs        !! mean growing season starting day for deciduous PFTs (doy).
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: resp_hetero          !! Heterotrophic respiration
                                                                               !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: leaching_daily       !! mineral nitrogen leached from the soil 
                                                                               !! (gN/m**2/day)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: emission_daily       !! volatile losses of nitrogen 
                                                                               !! (gN/m**2/day)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: n_input_daily        !! nitrogen inputs into the soil  (gN/m**2/day)
                                                                               !! NH4 and NOX from the atmosphere, NH4 from BNF,
                                                                               !! agricultural fertiliser as NH4/NO3 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: gpp_week             !! Weekly gross primary productivity  
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: resp_maint_week      !! Weekly maintenance respiration
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: light_tran_to_floor_season !! Mean seasonal fraction of light transmitted (unitless, 0-1)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: qsintveg             !! Water on vegetation due to interception @tex $(kg m^{-2})$ @endtex
    REAL(r_std),DIMENSION(:,:), INTENT(inout)          :: vegstress            !! Relative soil moisture (0-1, unitless) 
    REAL(r_std),DIMENSION(:,:), INTENT(inout)          :: humrel               !! Relative humidity. Not used in stomate (needed in age_class_distr)
    REAL(r_std),DIMENSION(:,:,:), INTENT(inout)        :: season_drought_legacy!! Mean growing season moisture availability
    REAL(r_std),DIMENSION(:,:,:,:), INTENT(inout)      :: n_input              !! Nitrogen inputs into the soil (gN/m**2/timestep)
    REAL(r_std),DIMENSION (:,:,:), INTENT (inout)      :: Light_Abs_Tot        !!Absorbed radiation per level for photosynthesis
    REAL(r_std),DIMENSION (:,:,:), INTENT (inout)      :: Light_Tran_Tot       !!Transmitted radiation per level for photosynthesis
    REAL(r_std),DIMENSION(:,:,:), INTENT(inout)        :: laieff_isotrop       !! Effective LAI
    REAL(r_std),DIMENSION(:,:),INTENT(inout)           :: maxvegstress_lastyear!! last year's maximum moisture availability
    REAL(r_std),DIMENSION(:,:),INTENT(inout)           :: maxvegstress_thisyear!! this year's maximum moisture availability
    REAL(r_std),DIMENSION(:,:),INTENT(inout)           :: minvegstress_lastyear!! last year's minimum moisture availability
    REAL(r_std),DIMENSION(:,:),INTENT(inout)           :: minvegstress_thisyear!! this year's minimum moisture availability
    REAL(r_std),DIMENSION(:,:),INTENT(inout)           :: maxgppweek_lastyear  !! last year's maximum weekly GPP
    REAL(r_std),DIMENSION(:,:),INTENT(inout)           :: maxgppweek_thisyear  !! this year's maximum weekly GPP
    REAL(r_std),DIMENSION(:,:),INTENT(inout)           :: wstress_month        !! Water stress factor, based on hum_rel_daily
                                                                               !! (unitless, 0-1)
    REAL(r_std),DIMENSION(:,:),INTENT(inout)           :: n_reserve_longterm   !! "Long term" (default 3 years) actual to potential N
                                                                               !! reserve pool (0-1, unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: n_reserve_balance    !! Actual to potential N reserve pool (unitless)
    REAL(r_std),DIMENSION(:,:),INTENT(inout)           :: maxfpc_lastyear      !! Last year's maximum foliage projected
                                                                               !! coverage for each natural PFT,
                                                                               !! @tex $(m^2 m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(:,:),INTENT(inout)           :: maxfpc_thisyear      !! this year's maximum foliage projected
                                                                               !! coverage for each natural PFT,
                                                                               !! @tex $(m^2 m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(:,:,:,:),INTENT(inout)       :: turnover_longterm    !! "Long term" turnover rate
                                                                               !! @tex $(gC m^{-2} year^{-1})$ @endtex
    REAL(r_std),DIMENSION(:,:,:),INTENT(inout)         :: dead_leaves          !! Dead leaves on ground, per PFT, metabolic 
                                                                               !! and structural,  
                                                                               !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(:,:),INTENT(inout)           :: grow_season_len      !! growing season length in days for deciduous PFTs. 
    REAL(r_std),DIMENSION(:,:,:,:),INTENT(inout)       :: us                   !! Water stress index for transpiration
                                                                               !! (by soil layer and PFT) (0-1, unitless)  
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: leaf_age_crit        !! critical leaf age (days)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: leaf_classes         !! width of each leaf age class (days)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: co2_fire             !! Carbon emitted into the atmosphere by 
                                                                               !! fire (living and dead biomass)  
                                                                               !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)   :: litterfuel           !!Dead litter fuel above ground. (gC m^{-2})
    REAL(r_std),DIMENSION(:,:), INTENT(inout)          :: count_daylight       !! Number of time steps dt_radia during daylight
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: gap_area_save        !! Stand gap created by more than 30% basal area loss per year.
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: total_ba_init        !! Total basal area per pft saved at the start of the year (m^{2}/m^{-2})


    !! 0.4 Local variables
    INTEGER(i_std)                                     :: ipts,ivm,igroup     !! Indeces(unitless)
    INTEGER(i_std)                                     :: iele,ipar,ipft      !! Indeces(unitless)
    INTEGER(i_std)                                     :: iagec,imbc,icir     !! Indeces(unitless)
    INTEGER(i_std)                                     :: ilit,ilev,icarb     !! Indeces(unitless)
    INTEGER(i_std)                                     :: ivma,ispec,igrn     !! Indeces(unitless)
    INTEGER(i_std)                                     :: old_pft             !! pft number before merge
    INTEGER(i_std)                                     :: new_pft             !! pft number to merge with
    INTEGER(i_std)                                     :: ph_new_pft          !! PROGRESSIVE_HARVEST: age-class-1 slot receiving the cut area
    REAL(r_std)                                        :: ph_f                !! PROGRESSIVE_HARVEST: slot area fraction clearcut last year (-)
    REAL(r_std)                                        :: ac12_frac           !! MATURITY_CONVEYOR: fraction de la classe 1 transferee (-)
    ! Guillaume M. -- MATURITY_CONVEYOR: the "partial transfer?" test used to be replayed six
    ! times on ac12_frac, which was reset to zero in the MIDDLE of the sequence. The last
    ! reads then saw 0 and concluded "whole transfer": veget_max_new went all-or-nothing,
    ! the donor age was erased and keep_source=.FALSE. emptied the soil, litter and bact of
    ! a donor that kept most of its area. The state is now computed ONCE.
    LOGICAL                                            :: ac12_partial        !! MATURITY_CONVEYOR: transfert partiel en cours (donneur survivant)
    REAL(r_std)                                        :: ac12_r              !! MATURITY_CONVEYOR: rapport de hauteur h1/h_ref (-)
    REAL(r_std)                                        :: ac12_dref           !! MATURITY_CONVEYOR: diametre de la classe de reference (m)
    REAL(r_std)                                        :: ac12_don0           !! AC12_F_REAL: veget_max du donneur AVANT le mouvement (-)
    REAL(r_std)                                        :: mat_a_grp           !! MATURITY_TRANSFER: total area of the age-class group (-)
    REAL(r_std)                                        :: mat_r_don           !! MATURITY_TRANSFER: donor filling ratio A/(target*A_grp) (-)
    REAL(r_std)                                        :: mat_r_rec           !! MATURITY_TRANSFER: receiver filling ratio (-)
    INTEGER(i_std)                                     :: ipft_grp            !! MATURITY_TRANSFER: loop index over the group slots
    INTEGER(i_std)                                     :: iagec_rec           !! MATURITY_TRANSFER: age-class rank of the receiver
    REAL(r_std)                                        :: mat_mov_n           !! MAT_SHIFT: stem factor of the moving stand (1+s)**(1/beta) (-)
    REAL(r_std)                                        :: mat_mov_b           !! MAT_SHIFT: per-tree mass factor of the moving stand (1+s)**(2+k3) (-)
    REAL(r_std)                                        :: mat_rem_n           !! MAT_SHIFT: stem factor of the donor, from the stem identity (-)
    REAL(r_std)                                        :: mat_rem_b           !! MAT_SHIFT: per-tree mass factor of the donor, from the mass identity (-)
    REAL(r_std)                                        :: ac12_thr            !! MATURITY_CONVEYOR: seuil applique a ac12_r (-)
    REAL(r_std)                                        :: ac12_v_mov          !! MATURITY_CONVEYOR: aire reellement deplacee (-)
    REAL(r_std)                                        :: ac12_v_mov_new      !! MATURITY_CONVEYOR: idem sur veget_max_new (-)
    REAL(r_std)                                        :: ac12_v_rec, ac12_v_rem !! MATURITY_CONVEYOR: aires APRES deplacement, receveur et donneur (-)
    LOGICAL                                            :: cop_2class          !! COPPICE_2CLASS: ce creneau est conduit en taillis a deux classes
    INTEGER(i_std)                                     :: ac12_rec            !! COPPICE_2CLASS: indice du creneau receveur d'une montee
    REAL(r_std)                                        :: ac0_sum             !! INSTRUMENTATION: ecart de fermeture somme sur les classes d'age
    ! Guillaume M. -- CLASS 0. The conveyor intercepts the clearcut: instead of handing the
    ! area back to class 1 the same day, it spends N_CLASS0 years on the dominant HERBACEOUS
    ! PFT of the grid cell. The register class0_veg(ipts,group,slot) remembers the crediting
    ! group -- traceability is carried by the register, not by the PFT dimension.
    INTEGER(i_std)                                     :: c0_g                !! CLASSE 0: PFT herbace receveur de la maille (0 = aucun)
    INTEGER(i_std)                                     :: c0_k                !! CLASSE 0: indice de balayage
    REAL(r_std)                                        :: c0_best             !! CLASSE 0: veget_max du candidat herbace dominant
    LOGICAL                                            :: c0_on               !! CLASSE 0: convoyeur actif ET alloue
    REAL(r_std)                                        :: c0_due              !! CLASSE 0: aire sortant du convoyeur cette annee (-)
    REAL(r_std)                                        :: c0_sum              !! CLASSE 0: total du registre, pour la garde de coherence
    REAL(r_std)                                        :: ph_v_mov            !! PROGRESSIVE_HARVEST: veget_max moved to age class 1 (-)
    REAL(r_std)                                        :: mae_age             !! MAT_AGE_ENTRY: age of the parcel moved upward (yr)
    REAL(r_std)                                        :: mae_decay           !! MAT_AGE_ENTRY: leak of the integrator (-)
    REAL(r_std)                                        :: ph_v_mov_new        !! PROGRESSIVE_HARVEST: same for veget_max_new (-)
    REAL(r_std)                                        :: ph_v_rec            !! PROGRESSIVE_HARVEST: receiver veget_max AFTER the move (-)
    REAL(r_std)                                        :: ph_v_rem            !! PROGRESSIVE_HARVEST: donor veget_max AFTER the move (-)
    INTEGER(i_std)                                     :: ph_nskip            !! PROGRESSIVE_HARVEST: splits refused this call (min_vegfrac invariant)
    REAL(r_std)                                        :: share_rec           !! Share of the veget_max of the existing vegetation
                                                                              !! within a PFT over the total veget_max following
                                                                              !! expansion of that PFT (unitless, 0-1)
    REAL(r_std)                                        :: aw_new, aw_old      !! STAND_AGE : poids d'AIRE du receveur et de la
                                                                              !! parcelle deplacee, pour le miroir age_stand_area (-)
    REAL(r_std)                                        :: bm_w_new, bm_w_old  !! STAND_AGE : carbone SUR PIED des deux classes,
                                                                              !! capture AVANT la fusion, pour ponderer le melange
                                                                              !! de age_stand_bm par la BIOMASSE et non par l'aire
                                                                              !! (gC, aire de maille incluse)
    REAL(r_std), DIMENSION(ncirc)                      :: circ_class_ba       !! basal area of the model tree in each
                                                                              !! circ class @tex $(m^{2} m^{-2})$ @endtex 
    REAL(r_std)                                        :: dfac_ac             !! Modulation des bornes de classe par l'intensite de gestion (-)
    REAL(r_std), DIMENSION(npts,nvm)                   :: max_dia             !! Largest diamater within circ_class (m)
    REAL(r_std), DIMENSION(ncirc)                      :: est_circ_class_n    !! Temporary variable for circ_class_n of the 
                                                                              !! established vegetation 
                                                                              !! @tex $(ind m^{-2})$ @endtex
    REAL(r_std), DIMENSION(ncirc)                      :: tmp_circ_class_n    !! Temporary variable for circ_class_n of the 
                                                                              !! established vegetation 
                                                                              !! @tex $(ind m^{-2})$ @endtex
    REAL(r_std), DIMENSION(ncirc,nparts,nelements)     :: est_circ_class_biomass !! Temporary variable for circ_class_biomass 
                                                                              !! of the established vegetation 
                                                                              !! @tex $(g C ind^{-1})$ @endtex
    REAL(r_std), DIMENSION(ncirc,nparts,nelements)     :: tmp_circ_class_biomass !! Temporary variable for circ_class_biomass 
                                                                              !! of the established vegetation
                                                                              !! @tex $(g C ind^{-1})$ @endtex
    REAL(r_std), DIMENSION(ncirc)                      :: est_circ_class_dia  !! Variable to store temporary values for 
                                                                              !! circ_class (m)
    REAL(r_std), DIMENSION(ncirc)                      :: tmp_circ_class_dia  !! Variable to store temporary values for 
                                                                              !! circ_class (m)
    REAL(r_std), DIMENSION(npts,nvm,nmbcomp,nelements) :: check_intern        !! Contains the components of the internal
                                                                              !! mass balance chech for this routine
                                                                              !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)         :: closure_intern      !! Check closure of internal mass balance
                                                                              !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)         :: pool_start          !! Start and end pool of this routine 
                                                                              !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)         :: pool_end            !! Start and end pool of this routine 
                                                                              !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(nelements)                  :: temp_start          !! Start and end pool of this routine 
                                                                              !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(nelements)                  :: temp_end            !! Start and end pool of this routine 
                                                                              !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(nlitt,nlevs)                :: litter_weight_rec   !! The fraction of litter on the expanded
                                                                              !! PFT. @tex $-$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)                   :: veget_max_begin     !! Temporairy variable to check area conservation
    REAL(r_std),DIMENSION(npts,nparts,nelements)       :: tree_biomass_loss   !! tree biomass loss @tex ($gC m^{-2}$) @endtex

    ! Guillaume M. -- scaffolding [C0IN]: yearly sum of the area registered into class 0
    ! (ph_v_mov). Compared with the capped c0_new in stomate_lpj, it says how much booked
    ! area the cap silently drops. Exit criterion: the gap is explained and closed.
    REAL(r_std) :: c0in_sum
    ! Guillaume M. -- scaffolding [C0BACK]: at the EXIT of this routine, compare what class 0
    ! has booked with what the class-1 slot actually carries, per point and per group. The
    ! excess is booked area with no backing, i.e. what the feeding cap will drop later.
    REAL(r_std) :: c0b_book, c0b_carry, c0b_excess
    INTEGER(i_std) :: c0b_i, c0b_k

!_ ================================================================================================================================

    c0in_sum = zero

    ! Guillaume M. -- c0_on is loop-invariant, but it used to be assigned inside the progressive
    ! harvest block, which is conditional. A slot reaching the class-0 reservation without that
    ! block having run read an UNDEFINED logical. Computed once here, before any read.
    c0_on = (n_class0 > 0 .AND. ALLOCATED(class0_veg))

    IF (firstcall_sapiens_lcchange) THEN
      !! Initialize local printlev
      ! Guillaume M. -- get_printlev calls getin_p, a COLLECTIVE bcast. Called here,
      ! inside the compute loop under a local firstcall, it hangs the run as soon as the
      ! ranks stop reaching it together. Use printlev, read once at startup by every rank.
      printlev_loc=printlev
      firstcall_sapiens_lcchange=.FALSE.
    END IF
  
    IF (printlev_loc.GE.2) WRITE(numout,*) 'Entering age class distribution'

 !! 1. Initialize
    
    !! 1.1 Initialize variables

    !! 1.2 Initialize check for mass balance closure
    !  The mass balance is calculated at the end of this routine
    !  in section 3. Initial biomass and harvest pool all other
    !  relevant pools were just set to zero.
    IF (err_act.GT.1) THEN

       check_intern(:,:,:,:) = zero
       pool_start(:,:,:) = zero
       DO iele = 1,nelements

          ! atm_to_bm
          pool_start(:,:,iele) = pool_start(:,:,iele) + &
               atm_to_bm(:,:,iele)
          
          !  Initial biomass + bm_to_litter
          DO ipar = 1,nparts

             DO icir = 1,ncirc

                ! Initial biomass pool
                pool_start(:,:,iele) = pool_start(:,:,iele) + &
                     circ_class_biomass(:,:,icir,ipar,iele) * &
                     circ_class_n(:,:,icir) * veget_max(:,:)

             ENDDO

             ! Add bm_to_litter to the initial biomass pool
             pool_start(:,:,iele) = pool_start(:,:,iele) + &
                  ((bm_to_litter(:,:,ipar,iele)+ &
                  turnover_daily(:,:,ipar,iele)) * veget_max(:,:))

          ENDDO

          ! Litter pool (gC m-2) *  (m2 m-2) 
          DO ilit = 1,nlitt
             DO ilev = 1,nlevs
                pool_start(:,:,iele) = pool_start(:,:,iele) + &
                     litter(:,ilit,:,ilev,iele) * veget_max(:,:)
             ENDDO
          ENDDO
       
          IF (ok_soil_carbon_discretization) THEN
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
          
       ENDDO

       ! Guillaume M. -- These two terms do NOT depend on iele (bact is carbon only,
       ! soil_n_min nitrogen only). Inside the DO iele loop they were added nelements times,
       ! i.e. twice. The closure was unaffected (the duplicate cancelled between start and
       ! end) but it doubled the weight of bact in the residual. The other mass-balance
       ! blocks of this file already place them outside the loop.
       !
       ! Denitrifying bacteria (only C)
       pool_start(:,:,icarbon) = pool_start(:,:,icarbon) + &
               bact(:,:) * veget_max(:,:)

       ! Nitrogen only
       DO ispec = 1,nnspec
          pool_start(:,:,initrogen) = pool_start(:,:,initrogen) + &
               soil_n_min(:,:,ispec) * veget_max(:,:)
       ENDDO

       !! 1.3 Initialize check for area conservation
       veget_max_begin(:,:) = veget_max(:,:)

    ENDIF ! err_act.GT.1

    ! Guillaume M. -- CLASS 0: the restitution / ageing / guard blocks live in stomate_lpj,
    ! between age_class_distr and land_cover_change_main. Moving area OUT of a forest group
    ! violates the invariant check_vegetation_area('ageclass') enforces here: the group
    ! total must be conserved. The conveyor is LAND COVER CHANGE, so it now expresses
    ! itself through veget_max_new and the LCC machinery performs the move.

 !! 2. Redistribute the biomass across the age classes

    !  Following growth, the trees from an age class may have become too big to
    !  belong to this age class. The biomass, litter, soilcarbon and soil
    !  water then need to be moved from one age class to the next age class.
    ph_nskip = 0
    DO ipts = 1,npts

       ! The agec_group lets us know which PFTs belong in which age class groups.
       !  I assume that the structure is linear, that is, if PFTs 16, 17, and 18 all
       ! below to group 10, PFT 16 is the youngest, PFT 17 is the second-youngest, and
       ! PFT 18 is the oldest (largest diameter).  I also assume that all PFTs in 
       ! an age class are grouped together in the input file.  In the example above,
       ! this means that PFT 19 cannot be in group 10 unless PFT 18 is also in group 10.   

       ! This loops over all the PFTs that we have while ignoring age classes
       DO ivma=1,nvmap

          ! start_index is an array with nvmap values wich gives the pft number
          ! of the youngest age class for the nvmap groups. Say we have 28 nvmap
          ! but run with a total of 64 pfts because some of the 28 groups have
          ! age class. then start_index will have 28 dimensions and the value of
          ! start_index will vary between 1 and 64.
          ivm=start_index(ivma)

          ! If we only have a single age class for this
          ! PFT, we can skip it. nagec_pft contains for each group either 1 
          ! (no age classes) or the number of age classes.
          IF(nagec_pft(ivma) .EQ. 1)CYCLE

          ! Only forest PFTs have age classes
          IF (is_tree(ivm)) THEN

             !! 2.1 Check if the trees still belong to this age class
             !  Note that the term age class is used but that the classes used in the
             !  code are not defined on an age criterion. Instead a diameter criterion
             !  is used. The reason to do so is that growth depends on the rule of
             !  Deleuze and Dhote that in turn is based on basal area. More stable
             !  model behaviour is therefore expected when the diameter range within an
             !  age class is controlled. Also different PFTs grow at different speeds.
             !  If age would have been used to define the classes, the diameter range
             !  within the first age class would have been much larger in the 
             !  tropics compared to the boreal zone. When Deleuze and Dhote is then
             !  applied the model could start to behave differently in the tropics and
             !  boreal zone. These differences in behaviour are thought to be better 
             !  controlled by implementing a diameter criterion.
             !  As a reminder age classes (::nagec) control the diameter/height
             !  variation across stands, ::ncirc controls the diameter/height variation
             !  within a stand
             DO iagec = nagec_pft(ivma),1,-1

                ipft = ivm+iagec-1

                ! If there are no trees, there is no need to merge
                IF (veget_max(ipts,ipft).LT.min_stomate) CYCLE

                ! Guillaume M. -- PROGRESSIVE_HARVEST: split the mature slot BEFORE the
                ! diameter test below. sapiens_kill has already removed a fraction f of the
                ! individuals over the UNCHANGED area V, which is a thinning. Here it becomes
                ! an area-based clearcut: f*V of bare ground moves down to age class 1 and the
                ! surviving stand is restored to density n over the remaining (1-f)*V.
                !
                ! Guillaume M. -- The moving part carries ZERO trees, so merge_biomass_pfts
                ! only dilutes the receiver density by share_rec. Litter and SOM are INTENSIVE
                ! (per m2), so the split conserves mass exactly: hence keep_source=.TRUE. on
                ! move_pft_properties -- zeroing the source would destroy (1-f)*V*L of carbon.
                !
                ! Guillaume M. -- The one-year lag is the model's own convention, so
                ! ph_cut_frac must survive the restart. /!\ Documented approximation: the cut
                ! residues get diluted between the two parts instead of sitting entirely on
                ! the cut part; mass is conserved, only their sub-grid allocation is
                ! approximate. .NOT.lpft_replant mirrors the guard of the diameter-based move.
                IF (ok_progressive_harvest .AND. ALLOCATED(ph_cut_frac) .AND. &
                     iagec .GT. 1 .AND. .NOT.lpft_replant(ipts,ipft)) THEN
                 IF (ph_cut_frac(ipts,ipft) .GT. zero .AND. ph_cut_frac(ipts,ipft) .LT. un) THEN

                   ph_f = ph_cut_frac(ipts,ipft)

                   ! Guillaume M. -- Consume the flag straight away: the split becomes
                   ! idempotent, so a second call in the same year cannot split it twice.
                   ph_cut_frac(ipts,ipft) = zero

                   ph_new_pft = start_index(ivma)
                   ph_v_mov = ph_f * veget_max(ipts,ipft)

                   ! --- CLASS 0 (step B): route the cut through the HERBACEOUS PFT --------
                   ! Guillaume M. -- Without the conveyor the cut area rejoins class 1 the
                   ! same day: the model TELEPORTS the stand. With it, the area spends
                   ! N_CLASS0 years as herbaceous cover, which carries the ET and the albedo;
                   ! bare ground would give a wrong energy balance over those years.
                   !
                   ! Guillaume M. -- Here we only REGISTER: the cut area rejoins class 1 as
                   ! before, and the conveyor records what it is owed. Moving area OUT of a
                   ! forest group would violate the invariant of this routine, which is that
                   ! the group total is conserved. The physical move to the herbaceous PFT is
                   ! done by the LAND COVER CHANGE machinery, via veget_max_new.
                   !
                   ! Guillaume M. -- /!\ min_vegfrac INVARIANT. check_veget forbids any
                   ! fraction STRICTLY between 0 and min_vegfrac = 1e-3, while min_stomate is
                   ! 1e-8: both sides of the split can violate it. The totals are therefore
                   ! tested AFTER the move, not the moved amount -- a small gain to a class 1
                   ! that already owns area is perfectly legal.
                   !
                   ! Guillaume M. -- The donor side is not optional either: the density
                   ! restore below divides by (un - ph_f), so a slot emptied down to nothing
                   ! would have its density amplified without bound. A slot cut entirely
                   ! belongs to the clearcut path. Refusing here is safe: ph_cut_frac was
                   ! already consumed, so the cut simply stays a THINNING.
                   !
                   ! Guillaume M. -- iagec > 1 already excludes it, but the explicit
                   ! ph_new_pft .NE. ipft guard documents that source and destination differ.
                   ph_v_rec = veget_max(ipts,ph_new_pft) + ph_v_mov
                   ph_v_rem = veget_max(ipts,ipft)       - ph_v_mov
                   IF (ph_new_pft .NE. ipft .AND. ph_v_mov .GT. min_stomate .AND. &
                        ph_v_rec .GE. min_vegfrac .AND. ph_v_rem .GE. min_vegfrac) THEN

                      CALL calculate_shares(veget_max(ipts,ph_new_pft), ph_v_mov, &
                           litter(ipts,:,ph_new_pft,:,:), litter(ipts,:,ipft,:,:), &
                           share_rec, litter_weight_rec)

                      est_circ_class_n(:) = circ_class_n(ipts,ph_new_pft,:)
                      est_circ_class_biomass(:,:,:) = circ_class_biomass(ipts,ph_new_pft,:,:,:)
                      tmp_circ_class_n(:) = zero
                      tmp_circ_class_biomass(:,:,:) = zero

                      ! Guillaume M. -- The receiver is class 1, a tree PFT, so the merge is
                      ! legitimate and NEEDED. wood_to_dia is never called on an herbaceous
                      ! PFT here: the transfer to the herbaceous cover is done by the land
                      ! cover change machinery, not by this block.
                      CALL merge_biomass_pfts(npts, share_rec, circ_class_n, &
                           est_circ_class_n, tmp_circ_class_n, circ_class_biomass, &
                           est_circ_class_biomass, tmp_circ_class_biomass, ipts, &
                           ph_new_pft, ipft, forest_managed, bm_to_litter)

                      ! Guillaume M. -- Restore the surviving stand's density (see above).
                      circ_class_n(ipts,ipft,:) = circ_class_n(ipts,ipft,:) / (un - ph_f)

                      veget_max(ipts,ph_new_pft) = veget_max(ipts,ph_new_pft) + ph_v_mov
                      veget_max(ipts,ipft)       = veget_max(ipts,ipft) - ph_v_mov

                      ! Guillaume M. -- Credit the register AFTER the move actually happened,
                      ! never before: the split can still be refused by the min_vegfrac
                      ! invariant, and booking area that did not move would silently
                      ! desynchronise the register from the physics. `ivma` is the crediting
                      ! GROUP: it is the one that gets the area back in N_CLASS0 years.
                      IF (c0_on) class0_veg(ipts,ivma,1) = class0_veg(ipts,ivma,1) + ph_v_mov
                      IF (c0_on) c0in_sum = c0in_sum + ph_v_mov
                      ph_v_mov_new = ph_f * veget_max_new(ipts,ipft)
                      veget_max_new(ipts,ph_new_pft) = veget_max_new(ipts,ph_new_pft) + ph_v_mov_new
                      veget_max_new(ipts,ipft)       = veget_max_new(ipts,ipft) - ph_v_mov_new

                      ! Guillaume M. -- STAND_AGE: the moving area is BARE. It brings no
                      ! biomass today but will be planted and carry a YOUNG stand.
                      !
                      ! Guillaume M. -- Blending with the DONOR age made class 1 inherit the
                      ! age of the stand that was cut; dropping the blend left the mirror
                      ! defect, no path rejuvenating the receiver -- age_stand_bm is only
                      ! reset on a fully empty slot. Blend with the ESTABLISHED stand age
                      ! instead, weighted by AREA: the area arrives, its biomass comes later.
                      IF (veget_max(ipts,ph_new_pft) > min_stomate) THEN
                         age_stand_bm(ipts,ph_new_pft) = &
                              ( (veget_max(ipts,ph_new_pft) - ph_v_mov) &
                                * age_stand_bm(ipts,ph_new_pft) &
                              + ph_v_mov * age_stand_estab ) &
                              / veget_max(ipts,ph_new_pft)
                         age_stand_area(ipts,ph_new_pft) = &
                              ( (veget_max(ipts,ph_new_pft) - ph_v_mov) &
                                * age_stand_area(ipts,ph_new_pft) &
                              + ph_v_mov * age_stand_estab ) &
                              / veget_max(ipts,ph_new_pft)
                      ENDIF
                      ! Guillaume M. -- age_stand_bm(ipts,ipft) is left unchanged: the
                      ! surviving stand has not become younger.

                     CALL move_pft_properties(ipts, ph_new_pft, ipft, share_rec, &
                          litter_weight_rec, circ_class_biomass, fm_change_map, &
                          wstress_season(ipts,ph_new_pft), wstress_season(ipts,ipft), &
                          wstress_month(ipts,ph_new_pft), wstress_month(ipts,ipft), &
                          season_drought_legacy(ipts,ph_new_pft,:), season_drought_legacy(ipts,ipft,:), &
                          nstress_season(ipts,ph_new_pft), nstress_season(ipts,ipft), &
                          n_input(ipts,ph_new_pft,:,:), n_input(ipts,ipft,:,:), &
                          n_input_daily(ipts,ph_new_pft,:), n_input_daily(ipts,ipft,:), & 
                          vegstress_season(ipts,ph_new_pft), vegstress_season(ipts,ipft), &
                          vegstress_month(ipts,ph_new_pft), vegstress_month(ipts,ipft), &
                          vegstress_week(ipts,ph_new_pft), vegstress_week(ipts,ipft), &
                          vegstress(ipts,ph_new_pft), vegstress(ipts,ipft), &
                          humrel(ipts,ph_new_pft), humrel(ipts,ipft), &
                          everywhere(ipts,ph_new_pft), everywhere(ipts,ipft), &
                          PFTpresent(ipts,ph_new_pft), PFTpresent(ipts,ipft), &
                          sugar_load(ipts,ph_new_pft), sugar_load(ipts,ipft), &
                          cn_leaf_min_season(ipts,ph_new_pft), cn_leaf_min_season(ipts,ipft), &
                          cn_leaf_init_2D, &
                          Light_Abs_Tot(ipts,ph_new_pft,:), Light_Abs_Tot(ipts,ipft,:), &
                          Light_Tran_Tot(ipts,ph_new_pft,:), Light_Tran_Tot(ipts,ipft,:), &
                          laieff_isotrop(ipts,:,ph_new_pft), laieff_isotrop(ipts,:,ipft), &
                          lm_lastyearmax(ipts,ph_new_pft), lm_lastyearmax(ipts,ipft), &
                          lm_thisyearmax(ipts,ph_new_pft), lm_thisyearmax(ipts,ipft), &
                          age(ipts,ph_new_pft), age(ipts,ipft), &
                          leaf_frac(ipts,ph_new_pft,:), leaf_frac(ipts,ipft,:), &
                          leaf_age(ipts,ph_new_pft,:), leaf_age(ipts,ipft,:), &
                          light_tran_to_floor_season(ipts,ph_new_pft), &
                          light_tran_to_floor_season(ipts,ipft), &
                          qsintveg(ipts,ph_new_pft), qsintveg(ipts,ipft), &
                          us(ipts,ph_new_pft,:,:), us(ipts,ipft,:,:), &
                          when_growthinit(ipts,ph_new_pft), when_growthinit(ipts,ipft), &
                          gdd_from_growthinit(ipts,ph_new_pft), gdd_from_growthinit(ipts,ipft), &
                          gdd_midwinter(ipts,ph_new_pft), gdd_midwinter(ipts,ipft), &
                          time_hum_min(ipts,ph_new_pft), time_hum_min(ipts,ipft), &
                          hum_min_dormance(ipts,ph_new_pft), hum_min_dormance(ipts,ipft), &
                          gdd_m5_dormance(ipts,ph_new_pft), gdd_m5_dormance(ipts,ipft), &
                          ncd_dormance(ipts,ph_new_pft), ncd_dormance(ipts,ipft), & 
                          ngd_minus5(ipts,ph_new_pft), ngd_minus5(ipts,ipft), &
                          mean_start_gs(ipts,ph_new_pft), mean_start_gs(ipts,ipft), &
                          maxvegstress_lastyear(ipts,ph_new_pft), maxvegstress_lastyear(ipts,ipft), &
                          maxvegstress_thisyear(ipts,ph_new_pft), maxvegstress_thisyear(ipts,ipft), &
                          minvegstress_lastyear(ipts,ph_new_pft), minvegstress_lastyear(ipts,ipft), &
                          minvegstress_thisyear(ipts,ph_new_pft), minvegstress_thisyear(ipts,ipft), &
                          maxgppweek_lastyear(ipts,ph_new_pft), maxgppweek_lastyear(ipts,ipft), & 
                          maxgppweek_thisyear(ipts,ph_new_pft), maxgppweek_thisyear(ipts,ipft), &
                          n_reserve_longterm(ipts,ph_new_pft), n_reserve_longterm(ipts,ipft), &
                          n_reserve_balance(ipts,ph_new_pft), n_reserve_balance(ipts,ipft), &
                          maxfpc_lastyear(ipts,ph_new_pft), maxfpc_lastyear(ipts,ipft), &
                          maxfpc_thisyear(ipts,ph_new_pft), maxfpc_thisyear(ipts,ipft), &
                          turnover_longterm(ipts,ph_new_pft,:,:), turnover_longterm(ipts,ipft,:,:), &
                          dead_leaves(ipts,ph_new_pft,:), dead_leaves(ipts,ipft,:), &
                          grow_season_len(ipts,ph_new_pft), grow_season_len(ipts,ipft), &
                          KF(ipts,ph_new_pft), KF(ipts,ipft), &
                          atm_to_bm(ipts,ph_new_pft,:), atm_to_bm(ipts,ipft,:), &
                          npp_longterm(ipts,ph_new_pft), npp_longterm(ipts,ipft), &
                          croot_longterm(ipts,ph_new_pft), croot_longterm(ipts,ipft), &
                          gpp_daily(ipts,ph_new_pft), gpp_daily(ipts,ipft), &
                          resp_maint(ipts,ph_new_pft), resp_maint(ipts,ipft), &
                          resp_growth(ipts,ph_new_pft), resp_growth(ipts,ipft), &
                          resp_hetero(ipts,ph_new_pft), resp_hetero(ipts,ipft), &
                          npp_daily(ipts,ph_new_pft), npp_daily(ipts,ipft), &
                          rue_longterm(ipts,ph_new_pft), rue_longterm(ipts,ipft), &
                          leaching_daily(ipts,ph_new_pft,:), leaching_daily(ipts,ipft,:), &
                          emission_daily(ipts,ph_new_pft,:), emission_daily(ipts,ipft,:), &
                          gpp_week(ipts,ph_new_pft), gpp_week(ipts,ipft), &
                          resp_maint_week(ipts,ph_new_pft), resp_maint_week(ipts,ipft), &
                          gpp_year(ipts,ph_new_pft), gpp_year(ipts,ipft), &
                          gpp_decade(ipts,ph_new_pft), gpp_decade(ipts,ipft), &
                          mai(ipts,ph_new_pft), mai(ipts,ipft), &
                          pai(ipts,ph_new_pft), pai(ipts,ipft), &
                          mai_count(ipts,ph_new_pft), mai_count(ipts,ipft), &
                          previous_wood_volume(ipts,ph_new_pft), previous_wood_volume(ipts,ipft), &
                          age_stand(ipts,ph_new_pft), age_stand(ipts,ipft), &
                          last_cut(ipts,ph_new_pft), last_cut(ipts,ipft), &
                          k_latosa_adapt(ipts,ph_new_pft), k_latosa_adapt(ipts,ipft), &
                          litter(ipts,:,ph_new_pft,:,:), litter(ipts,:,ipft,:,:), &
                          bm_to_litter(ipts,ph_new_pft,:,:), bm_to_litter(ipts,ipft,:,:), &
                          tree_bm_to_litter(ipts,ph_new_pft,:,:), tree_bm_to_litter(ipts,ipft,:,:), &
                          leaf_age_crit(ipts,ph_new_pft), leaf_age_crit(ipts,ipft), &
                          leaf_classes(ipts,ph_new_pft), leaf_classes(ipts,ipft), &
                          co2_fire(ipts,ph_new_pft), co2_fire(ipts,ipft), &
                          litterfuel(ipts,:,ph_new_pft,:,:), litterfuel(ipts,:,ipft,:,:), &
                          turnover_daily(ipts,ph_new_pft,:,:), turnover_daily(ipts,ipft,:,:), &
                          soil_n_min(ipts,ph_new_pft,:), soil_n_min(ipts,ipft,:), &
                          p_O2(ipts,ph_new_pft), p_O2(ipts,ipft), &
                          bact(ipts,ph_new_pft), bact(ipts,ipft), &
                          deepSOM_a(ipts,:,ph_new_pft,:), deepSOM_a(ipts,:,ipft,:), &
                          deepSOM_s(ipts,:,ph_new_pft,:), deepSOM_s(ipts,:,ipft,:), &
                          deepSOM_p(ipts,:,ph_new_pft,:), deepSOM_p(ipts,:,ipft,:), &
                          som(ipts,:,ph_new_pft,:), som(ipts,:,ipft,:), &
                          lignin_struc(ipts,ph_new_pft,:), lignin_struc(ipts,ipft,:), &
                          lignin_wood(ipts,ph_new_pft,:), lignin_wood(ipts,ipft,:), &
                          lignin_snag(ipts,ph_new_pft,:), lignin_snag(ipts,ipft,:), &
                          forest_managed(ipts,ph_new_pft), forest_managed(ipts,ipft), &
                          plant_status(ipts,ph_new_pft), plant_status(ipts,ipft), &
                          kill_vessels(ipts,ph_new_pft,:), kill_vessels(ipts,ipft,:), &
                          vessel_loss_previous(ipts,ph_new_pft,:), vessel_loss_previous(ipts,ipft,:), &
                          biomass_init_drought(ipts,ph_new_pft,:,:,:), biomass_init_drought(ipts,ipft,:,:,:), &
                          count_daylight(ipts,ph_new_pft), count_daylight(ipts,ipft), &
                          CN_som_litter_longterm(ipts,ph_new_pft,:), CN_som_litter_longterm(ipts,ipft,:), &
                          matrixV(ipts,ph_new_pft,:,:), matrixV(ipts,ipft,:,:), &
                          matrixA(ipts,ph_new_pft,:,:), matrixA(ipts,ipft,:,:), &
                          vectorU(ipts,ph_new_pft,:), vectorU(ipts,ipft,:), &
                          vectorB(ipts,ph_new_pft,:), vectorB(ipts,ipft,:), &
                          gap_area_save(ipts,ph_new_pft,:), gap_area_save(ipts,ipft,:), &
                          total_ba_init(ipts,ph_new_pft), total_ba_init(ipts,ipft), &
                          keep_source=.TRUE.)

                   ELSE

                      ! Guillaume M. -- Split refused: it would have left a fraction strictly
                      ! between 0 and min_vegfrac, so the cut stays a THINNING. /!\ These
                      ! refusals are COUNTED: a progressive harvest silently degrading into a
                      ! thinning would look active while moving no area at all.
                      IF (ph_v_mov .GT. min_stomate) ph_nskip = ph_nskip + 1

                   ENDIF

                 ENDIF
                ENDIF

                !! 2.1.1 Basal area at the tree level (m2 tree-1)
                circ_class_ba(:) = &
                     wood_to_ba(circ_class_biomass(ipts,ipft,:,:,icarbon),&
                     ipft,pipe_tune2(ipts,ipft))

                ! Previously, quadratic mean diameter was used but that can
                ! include trees larger than age_class_bound. So for stricter
                ! distribution, maximum diameter within circ classes will be
                ! used instead.
                IF (SUM(circ_class_n(ipts,ipft,:)) .NE. zero) THEN

                   max_dia(ipts,ipft) = MAXVAL(&
                        wood_to_dia(circ_class_biomass(ipts,ipft,:,:,icarbon),&
                        ipft,pipe_tune2(ipts,ipft)))

                ELSE

                   max_dia(ipts,ipft) = zero

                ENDIF

                ! Guillaume M. -- Class bounds modulated by the management intensity: they
                ! are fractions of the cut diameter, so they follow its modulation and the
                ! whole maturity ladder of the pixel is rescaled.
                dfac_ac = un
                IF (ALLOCATED(dia_factor)) dfac_ac = dia_factor(ipts)

                !! 2.2.2 Check whether quadratic mean diameter exceeds boundaries 
                !  Check whether quadratic mean diameter exceeds boundaries of
                !  the age class.

                !Debug
                IF(printlev_loc>=4)THEN
                   WRITE(numout,*) 'Checking to merge for: '
                   WRITE(numout,*) 'ipft,ivm,iagec,ipts: ',ipft,ivm,iagec,ipts
                   WRITE(numout,*) 'nagec_pft,max_dia,age_class_bound: ',nagec_pft(ivma),&
                        max_dia(ipts,ipft),age_class_bound(iagec,ipft)*dfac_ac
                ENDIF
                !-
                
                ! If we flagged a species to be replanted the biomass
                ! pools of both the new and the old pft can be empty.
                ! Redistributing empty pools doesn't make too much sense.
                ! Replanting is taken care of in species_change.
                IF(.NOT.lpft_replant(ipts,ipft))THEN

                   ! Debug
                   IF(printlev_loc>=4) WRITE(numout,*) 'Ready to redistribute', ipts,ivm
                   !-

                   ! Guillaume M. -- /!\ MANDATORY RESET before the IF/ELSEIF chain: a value
                   ! left over from a previous slot would trigger a PARTIAL transfer on a slot
                   ! handled by another branch.
                   ac12_frac = zero
                   IF (.NOT. ALLOCATED(ac12_r_diag)) THEN
                      ALLOCATE(ac12_r_diag(npts,nvm)) ; ac12_r_diag(:,:) = zero
                      ALLOCATE(ac12_freal_diag(npts,nvm)) ; ac12_freal_diag(:,:) = zero
                   ENDIF
                   ac12_r_diag(ipts,ipft) = zero
                   ac12_freal_diag(ipts,ipft) = zero

                   ! Guillaume M. -- COPPICE_2CLASS: a coppice has two maturity classes, not
                   ! nagec. The diameter ladder is a high-forest object; a coppice keeps its
                   ! stems below the first rung for ever, so the demotion below would send
                   ! every mature slot back to class 1 and hold it there. /!\ Keyed on the
                   ! SLOT, never on the PFT group: nothing forces a group to be uniform.
                   cop_2class = ok_coppice_2class .AND. &
                        (forest_managed(ipts,ipft) .EQ. ifm_cop .OR. &
                         forest_managed(ipts,ipft) .EQ. ifm_src)
                   ! Guillaume M. -- Receiver of an upward move: the class above, or the
                   ! terminal class when the intermediate ones do not exist for this slot.
                   ac12_rec = ipft + 1
                   IF (cop_2class) ac12_rec = start_index(ivma) + nagec_pft(ivma) - 1

                   IF ( (iagec .EQ. nagec_pft(ivma)) .AND. &
                        max_dia(ipts,ipft) .GT. age_class_bound(iagec,ipft)*dfac_ac ) THEN
                      
                      ! If these conditions are satisfied our tree diameter is
                      ! very unrealistic
                      WRITE(numout,*) 'WARNING: tree diameter exceeds (m): ', &
                           age_class_bound(iagec,ipft)*dfac_ac
                      CALL ipslerr_p (3,'ERROR: age_class_distr',&
                           'tree diameter exceeds the threshold','','')

                      old_pft = ipft
                      new_pft = ipft
  
                   ELSEIF ( (iagec .NE. nagec_pft(ivma)) .AND. &
                        .NOT. cop_2class .AND. &
                        .NOT. (ok_maturity_transfer .AND. &
                               nagec_pft(ivma) > 1) .AND. &
                        max_dia(ipts,ipft) .GT. age_class_bound(iagec,ipft)*dfac_ac ) THEN

                      ! If this condition is satisfied, the trees are too
                      ! big for the current age class. They will be moved
                      ! one age class up
                      ! Guillaume M. -- This all-or-nothing bound is now the fallback used
                      ! ONLY when the progressive rule is off. Every non-terminal class is
                      ! excluded when it is on: tested first, the bound shadowed the
                      ! progressive rule and emptied the whole slot, carrying away the area
                      ! just booked into class 0.
                      old_pft = ipft
                      new_pft = ipft+1

                   ELSEIF ( cop_2class .AND. (iagec .NE. 1) .AND. &
                        (max_dia(ipts,ipft) / &
                         MAX(age_class_bound(1,ipft)*dfac_ac, min_stomate)) &
                        ** pipe_tune3(ipft) .LT. cop_r_down ) THEN

                      ! Guillaume M. -- COPPICE_2CLASS: demotion judged on the SAME height
                      ! ratio as the promotion, with a strictly lower threshold. The diameter
                      ! test further down would fire on every coppice slot and undo the
                      ! promotion the year after. /!\ Placed BEFORE the progressive branch on
                      ! purpose: that branch would swallow a razed class 2 or 3.
                      old_pft = ipft
                      new_pft = start_index(ivma)

                   ELSEIF ( (ok_maturity_transfer .OR. cop_2class) .AND. &
                        iagec .NE. nagec_pft(ivma) .AND. &
                        nagec_pft(ivma) > 1 .AND. &
                        ( cop_2class .OR. &
                          .NOT. ( iagec .NE. 1 .AND. max_dia(ipts,ipft) .LT. &
                                  age_class_bound(1,ipft)*dfac_ac &
                                  *ac12_dia_ratio_down ) ) ) THEN

                      ! Guillaume M. -- DEMOTION WINS over promotion. Without the last test
                      ! this branch would catch every non-terminal class and starve the
                      ! branch below, which is the ONLY path by which storm and beetle send
                      ! a wiped-out stand back to class 1. A slot that lost its trees must
                      ! go down, never up.

                      ! Guillaume M. -- MATURITY_CONVEYOR: PROGRESSIVE promotion, driven by a
                      ! height DIFFERENTIAL -- what creates an edge is the differential, not an
                      ! absolute height. /!\ Computed in DIAMETER: H = pipe_tune2*D^pipe_tune3
                      ! so r = (D1/D_ref)^pipe_tune3 and pipe_tune2 CANCELS OUT, which avoids
                      ! depending on that unreliable parameter.
                      !
                      ! Guillaume M. -- Reference is the CLASS BOUNDARY, not the diameter of
                      ! the first non-empty class above. That search locked the transfer out:
                      ! once class 2 empties, the reference comes from the far larger class 3
                      ! or 4, so r stalls near 0.5 and never reaches the start threshold. The
                      ! boundary makes r tend to 1 as the stand reaches its own class size.
                      !
                      ! Guillaume M. -- Applies to EVERY non-terminal class, not just class 1.
                      ! The reference is that class's own bound, which already carries a
                      ! silvicultural meaning: merchantable size for 2->3, clearcut size for
                      ! 3->4 -- the latter derived from DIA_ROTATION_TOL, not a free knob.
                      ac12_dref = age_class_bound(iagec,ipft)*dfac_ac
                      ! Guillaume M. -- Class 1 graduates on CANOPY CLOSURE, a height matter,
                      ! hence the allometric exponent. Above it the trigger is a merchantable
                      ! SIZE, judged on diameter. Dropping the exponent also makes the
                      ! threshold uniform: pipe_tune3 spans 0.42-0.52, so one height ratio
                      ! meant 0.66 to 0.72 in diameter depending on the species.
                      !
                      ! Guillaume M. -- COPPICE_2CLASS: one rung, judged in HEIGHT, whatever
                      ! the slot index. The reference is the class-1 bound because that is
                      ! the canopy-closure rung -- the only one a coppice ever meets. The
                      ! merchantable-size rungs above it describe a high forest and would
                      ! never be crossed by a stand of many thin shoots.
                      IF (cop_2class) THEN
                         ! Guillaume M. -- Reference is the regime's OWN CUT DIAMETER, not the
                         ! class-1 bound. With the bound, the terminal class was fed a
                         ! permanent stream of immature stems whose merge diluted its
                         ! diameter, so it never ripened enough to be cut and accumulated
                         ! past the self-thinning line.
                         !
                         ! Guillaume M. -- "Mature enough to enter the harvestable class" and
                         ! "big enough to be cut" must be ONE threshold, as the regular high
                         ! forest already does with its last class bound. Falls back to the
                         ! class-1 bound where no cut diameter is defined.
                         ac12_dref = age_class_bound(1,ipft)*dfac_ac
                         IF (ALLOCATED(coppice_diameter)) THEN
                            IF (coppice_diameter(ipft) > min_stomate) &
                                 ac12_dref = coppice_diameter(ipft)
                         ENDIF
                         ac12_r    = (max_dia(ipts,ipft) / MAX(ac12_dref, min_stomate)) &
                                     ** pipe_tune3(ipft)
                         ac12_thr  = cop_r_up
                      ELSEIF (iagec .EQ. 1) THEN
                         ac12_r   = (max_dia(ipts,ipft) / MAX(ac12_dref, min_stomate)) &
                                    ** pipe_tune3(ipft)
                         ac12_thr = ac12_h_ratio_start
                      ELSE
                         ac12_r   = max_dia(ipts,ipft) / MAX(ac12_dref, min_stomate)
                         ac12_thr = ac12_dia_ratio_start
                         ! Guillaume M. -- The step feeding the TERMINAL class gets its own
                         ! threshold: it is the only one whose receiver is also the harvestable
                         ! set, so the slot survives only if supply outpaces the cut. Falls back
                         ! to the shared value, hence bit-neutral until set.
                         IF (iagec .EQ. nagec_pft(ivma) - 1) ac12_thr = ac12_dia_ratio_last
                      ENDIF
                      IF (ac12_r >= ac12_thr) THEN
                         IF (ok_maturity_transfer) THEN
                            ! Guillaume M. -- MATURITY_TRANSFER: the rate is PRESCRIBED, not
                            ! derived from the size distribution. The threshold above decides
                            ! IF the stand is ripe; this rate decides HOW MUCH area moves.
                            ! What moves is an intact copy (see the assembly block below), so
                            ! no auto-limitation through max_dia exists on this path: the
                            ! only regulator is the area budget itself.
                            ac12_frac = MIN(un, MAX(zero, mat_f_ref(iagec)))
                            ! Guillaume M. -- Area regulator, dead band + cross correction.
                            ! Filling ratios are read on the CURRENT areas: the loop over
                            ! slots is descending, so the receiver has already given before
                            ! it receives and r_rec measures its real need. The correction
                            ! pushes when the donor is full AND the receiver empty; either
                            ! ratio alone would starve the terminal class or drain the
                            ! upstream ones. MAT_GAIN = 0 reproduces the prescribed rate
                            ! exactly (A/B lever).
                            mat_a_grp = zero
                            DO ipft_grp = start_index(ivma), &
                                 start_index(ivma) + nagec_pft(ivma) - 1
                               mat_a_grp = mat_a_grp + veget_max(ipts,ipft_grp)
                            ENDDO
                            IF (mat_gain > zero .AND. mat_a_grp > min_stomate) THEN
                               iagec_rec = ac12_rec - start_index(ivma) + 1
                               ! Guillaume M. -- Same effective setpoints as the harvest-feed
                               ! coupling (MAT_A3_ENTRY), so donor and receiver stay consistent.
                               mat_r_don = veget_max(ipts,ipft) / &
                                    MAX(mat_target_frac_eff(ipts, ipft, iagec) * mat_a_grp, min_stomate)
                               mat_r_rec = veget_max(ipts,ac12_rec) / &
                                    MAX(mat_target_frac_eff(ipts, ipft, iagec_rec) * mat_a_grp, min_stomate)
                               IF (mat_r_rec < mat_band_low .OR. &
                                    mat_r_rec > mat_band_high) THEN
                                  ! Guillaume M. -- An empty receiver sends the raw ratio to
                                  ! infinity; the MAT_F_MAX clamp below is the real limiter.
                                  ac12_frac = mat_f_ref(iagec) * &
                                       (mat_r_don / MAX(mat_r_rec, min_stomate))**mat_gain
                               ENDIF
                               ac12_frac = MIN(mat_f_max, MAX(mat_f_min, ac12_frac))
                            ENDIF
                         ENDIF
                         ! Guillaume M. -- Area booked into class 0 this year is a fresh
                         ! clearcut on its way to the herbaceous PFT: it is not maturing and
                         ! the conveyor must still find it in class 1. Class 1 ONLY: the
                         ! register is parked there, so subtracting it from any higher class
                         ! would shrink a transfer for nothing.
                         IF (c0_on .AND. iagec == 1 .AND. &
                              veget_max(ipts,ipft) > min_stomate) THEN
                            ac12_frac = ac12_frac * &
                                 MAX(zero, veget_max(ipts,ipft) - class0_veg(ipts,ivma,1)) &
                                 / veget_max(ipts,ipft)
                            ! Guillaume M. -- No cap below one: the prescribed rate
                            ! (MAT_F_REF, clamped by MAT_F_MIN/MAX above) is its own bound.
                            ac12_frac = MIN(un, MAX(zero, ac12_frac))
                         ENDIF
                         ! Guillaume M. -- /!\ min_vegfrac INVARIANT, same rule as the
                         ! progressive harvest. check_veget forbids any fraction STRICTLY
                         ! between 0 and min_vegfrac = 1e-3; min_stomate is far too low. Both
                         ! TOTALS are tested AFTER the move: a small gain to a class 2 that
                         ! already owns area is legal, dropping a slot below it is not.
                         ac12_v_rec = veget_max(ipts,ac12_rec) + ac12_frac*veget_max(ipts,ipft)
                         ac12_v_rem = veget_max(ipts,ipft)*(un - ac12_frac)
                         IF (ac12_v_rec < min_vegfrac .OR. ac12_v_rem < min_vegfrac) THEN
                            ! Guillaume M. -- Refusing the transfer must mean DOING NOTHING.
                            ! Zeroing ac12_frac alone left old_pft/new_pft set, so the
                            ! downstream ac12_partial test took its ELSE branch and class 1
                            ! moved up AS A BLOCK -- promoting an immature stand. old_pft =
                            ! new_pft = zero is the "no transfer" idiom of this chain.
                            ac12_frac = zero
                            ac12_r_diag(ipts,ipft) = ac12_r
                            old_pft = zero
                            new_pft = zero
                         ELSE
                            ac12_r_diag(ipts,ipft) = ac12_r
                            old_pft = ipft
                            new_pft = ac12_rec
                         ENDIF
                      ELSE
                         ! Guillaume M. -- r is filled in even below the threshold: it is what
                         ! says WHY nothing moves, and "not ripe yet" must be told apart from
                         ! "rule inactive".
                         ac12_r_diag(ipts,ipft) = ac12_r
                         old_pft = zero
                         new_pft = zero
                      ENDIF

                   ELSEIF ( .NOT. cop_2class .AND. (iagec .NE. 1) .AND. &
                        max_dia(ipts,ipft) .LT. age_class_bound(1,ipft)*dfac_ac &
                        *ac12_dia_ratio_down ) THEN

                      ! If this condition is satisfied, the trees are too
                      ! small for the current age class (after a harvest)
                      ! and they will be moved to the first age class. We
                      ! don't want this to happen if we are in the first
                      ! age class.
                      ! Guillaume M. -- HYSTERESIS: the threshold is a FRACTION of the bound,
                      ! not the bound itself. Promotion fires below the bound (height ratio
                      ! 0.85 = 0.66-0.72 in diameter), so an equal threshold demoted the same
                      ! year everything it had just promoted. The twin test in the progressive
                      ! branch above carries the same factor -- changing one alone opens a band
                      ! where a slot is neither promoted nor demoted.
                      ! Guillaume M. -- The guard used to read iagec .NE. start_index(ivma),
                      ! comparing an age-class index (1..nagec) with a PFT index (1..nvm).
                      ! It excluded the wrong class in every group whose start_index is at
                      ! most nagec: there, an oversmall class-2 stand was never demoted.
                      old_pft = ipft
                      new_pft = start_index(ivma)
                  
                      ! Debug
                      IF(printlev_loc>=4)THEN
                         WRITE(numout,*) 'age classes going down'
                         WRITE(numout,*) 'iagec, ipft, start_index, ', &
                              iagec, ipft, start_index(ivma)
                         WRITE(numout,*) 'max_dia, ', max_dia(ipts,ipft)
                         WRITE(numout,*) 'age_class_bound, ', age_class_bound(1,ipft)*dfac_ac
                         WRITE(numout,*) 'old_pft, new_pft, ',old_pft, new_pft
                      ENDIF
                      !-

                   ELSE

                      ! There is no need to change the age class. Set 
                      ! old and new pft to the same value. That way
                      ! the next IF-loop is skipped.
                      old_pft = zero
                      new_pft = zero
                      
                   ENDIF ! age_class_bounds

                ELSE

                   ! lpft_replant is TRUE. Don't change the age class. 
                   ! Set old and new pft to the same value. That way
                   ! the next IF-loop is skipped.
                   old_pft = zero
                   new_pft = zero

                ENDIF !lpft_replant


                IF (old_pft .NE. new_pft) THEN

                   ! Debug
                   IF(printlev_loc>=4 .AND. old_pft .NE. new_pft)THEN
                      WRITE(numout,*) 'Merging biomass'
                      WRITE(numout,*) 'ipts,ipft,iagec: ',ipts,ipft,iagec
                      WRITE(numout,*) 'age_class_bound(:): ',age_class_bound(:,ipft)*dfac_ac
                      WRITE(numout,*) 'max_dia: ',max_dia(ipts,ipft)
                      WRITE(numout,*) 'old_pft, new_pft, ',old_pft, new_pft
                      WRITE(numout,*) 'veget_max(new_pft), veget_max(old_pft), ', &
                           veget_max(ipts,new_pft), veget_max(ipts,old_pft)
                   ENDIF
                   !-

                   !! 2.2.3 Merge biomass
                   !  Biomass of two age classes needs to be merged. The established
                   !  vegetation is stored in new_pft, the new vegetation is stored in
                   !  old_pft
                   ! Guillaume M. -- MATURITY_CONVEYOR: share_rec = surf_rec/(surf_rec+surf_mov).
                   ! Passing the donor's TOTAL area while moving only a fraction would dilute the
                   ! receiver by 1/ac12_frac. /!\ The defect would be SILENT: merge_biomass_pfts
                   ! renormalises, so mass would be conserved and the ERR_ACT=3 check would pass;
                   ! only the receiver diameter would be wrong.
                   !
                   ! Guillaume M. -- /!\ STATE FROZEN ONCE, here, before the first read. Every
                   ! decision below (area moved, donor wipe, veget_max, veget_max_new,
                   ! age_stand_bm, keep_source) reads ac12_partial and not ac12_frac. That is
                   ! what makes the desynchronisation that broke the mass balance impossible.
                   ac12_partial = (ac12_frac > zero .AND. ac12_frac < un)
                   IF (ac12_partial) THEN
                      ac12_v_mov = ac12_frac * veget_max(ipts,old_pft)
                   ELSE
                      ac12_v_mov = veget_max(ipts,old_pft)
                   ENDIF
                   CALL calculate_shares(veget_max(ipts,new_pft), ac12_v_mov, &
                        litter(ipts,:,new_pft,:,:), litter(ipts,:,old_pft,:,:), &
                        share_rec, litter_weight_rec)

                   ! Guillaume M. -- STAND_AGE: standing carbon of both classes, captured HERE
                   ! because further down the merge has already emptied old_pft, so the weight
                   ! could no longer be computed.
                   bm_w_new = SUM(SUM(circ_class_biomass(ipts,new_pft,:,:,icarbon),2) * &
                        circ_class_n(ipts,new_pft,:)) * veget_max(ipts,new_pft)
                   bm_w_old = SUM(SUM(circ_class_biomass(ipts,old_pft,:,:,icarbon),2) * &
                        circ_class_n(ipts,old_pft,:)) * veget_max(ipts,old_pft)
                   ! Guillaume M. -- On a PARTIAL transfer only the moved parcel carries its
                   ! age into the blend: the parcel is a proportional slice, so its standing
                   ! carbon is ac12_frac times the class total. Weighting by the WHOLE donor
                   ! class dragged the receiver ~1/ac12_frac (15-30x) too hard every year and
                   ! collapsed the four class ages onto each other (1-2 yr apart).
                   IF (ac12_partial) bm_w_old = ac12_frac * bm_w_old
                   ! Guillaume M. -- AREA weights for the age_stand_area mirror: receiver
                   ! area and MOVED area (the parcel, not the whole donor class), both
                   ! captured pre-merge for the same reason as bm_w_*.
                   aw_new = veget_max(ipts,new_pft)
                   aw_old = ac12_v_mov


                   IF (trusting_hack_age_class) THEN
                      ! The calculation of the the age_class_distribution (even the multiplication with zero)
                      ! could result in a precision error. No calculations, just move all the information.
                      ! This approach only makes sense without LCC, management, disturbances, ... it is
                      ! and idealized test case that should only be used for trusting the model.
                      circ_class_biomass(ipts,new_pft,:,:,:) = circ_class_biomass(ipts,old_pft,:,:,:)
                      circ_class_biomass(ipts,old_pft,:,:,:) = zero
                      circ_class_n(ipts,new_pft,:) = circ_class_n(ipts,old_pft,:)
                      circ_class_n(ipts,old_pft,:) = zero

                   ELSE

                      ! Copy the number of individuals and the biomass of the established 
                      ! vegetation to a temporary variable. In the subroutine ::circ_class_n
                      ! and ::circ_class_biomass will be overwritten with the characteristics
                      ! of the merged vegetation
                      est_circ_class_n(:) = circ_class_n(ipts,new_pft,:)
                      est_circ_class_biomass(:,:,:) = circ_class_biomass(ipts,new_pft,:,:,:)
                      tmp_circ_class_n(:) = circ_class_n(ipts,old_pft,:)
                      tmp_circ_class_biomass(:,:,:) = circ_class_biomass(ipts,old_pft,:,:,:)
                      ! Guillaume M. -- MAT_SHIFT: the moving stand is displaced FORWARD
                      ! along the self-thinning trajectory (diameter (1+s), so stems
                      ! (1+s)**(1/beta) and per-tree mass (1+s)**(2+k3)); the donor is
                      ! rescaled backward by the two conservation identities, exact per
                      ! class for stems AND mass. This is what keeps adjacent maturity
                      ! classes apart -- the area regulator alone lets them converge.
                      IF (ok_maturity_transfer .AND. ac12_partial .AND. &
                           mat_shift > zero) THEN
                         mat_mov_n = (un + mat_shift)**(un/beta_self_thinning(old_pft))
                         mat_mov_b = (un + mat_shift)**(deux + pipe_tune3(old_pft))
                         ! Guillaume M. -- Stem identity: f*(1+nu) + (1-f)*(1-mu) = 1.
                         mat_rem_n = (un - ac12_frac*mat_mov_n) / (un - ac12_frac)
                         ! Guillaume M. -- Mass identity: f*mov_n*mov_b + (1-f)*rem_n*rem_b
                         ! = 1. The remainder leaves the trajectory by a second-order
                         ! amount in f*s: slightly younger AND slightly denser is the
                         ! price of exact conservation with a single prescribed step.
                         mat_rem_b = (un - ac12_frac*mat_mov_n*mat_mov_b) / &
                              MAX(min_stomate, (un - ac12_frac)*mat_rem_n)
                         IF (mat_rem_n <= zero .OR. mat_rem_b <= zero) THEN
                            CALL ipslerr_p(3,'age_class_distr', &
                                 'MAT_SHIFT leaves the donor with nothing', &
                                 'f*(1+s)**(1/beta) or the mass identity exceeded one', &
                                 'lower MAT_SHIFT or MAT_F_MAX')
                         ENDIF
                         tmp_circ_class_n(:) = tmp_circ_class_n(:) * mat_mov_n
                         tmp_circ_class_biomass(:,:,:) = &
                              tmp_circ_class_biomass(:,:,:) * mat_mov_b
                      ENDIF

                      ! Debug
                      ! Initiliaze mass balance of the merging procedure
                      ! MB was calculated outside the subroutine merge_biomass_pfts
                      ! because veget_max needs to be updated. Veget_max cannot be 
                      ! updated in the subroutine because the way it is called from 
                      ! the subroutine land_cover_change (where veget_max can only be
                      ! updated at the end of the whole subroutine. Therefore, this
                      ! check was included here.
                      IF (printlev_loc>=4) THEN
                         temp_start(:) = zero
                         DO ipar = 1,nparts
                            DO iele = 1,nelements
                               DO icir = 1,ncirc
                                  temp_start(iele) = temp_start(iele) + &
                                       SUM(circ_class_biomass(ipts,:,icir,ipar,iele) * &
                                       circ_class_n(ipts,:,icir) * veget_max(ipts,:))
                               ENDDO
                            ENDDO
                         ENDDO
                         ! Some output - unrelated to the mass balance
                         WRITE(numout,*) 'age distribution share_rec, ',ipts, ivm, &
                              share_rec
                         WRITE(numout,*) 'rest, ', SUM(circ_class_n(ipts,ivm,:)), &
                              SUM(est_circ_class_n), SUM(tmp_circ_class_n), &
                              SUM(SUM(circ_class_biomass(ipts,ivm,:,:,icarbon),2),1), &
                              SUM(est_circ_class_biomass), &
                              SUM(tmp_circ_class_biomass), &
                              ipts, new_pft, SUM(est_circ_class_dia), SUM(tmp_circ_class_dia), &
                              lpft_replant(ipts,ivm)
                         WRITE(numout,*) 'old_pft, new_pft, ',old_pft, new_pft
                      ENDIF
                      !-
                  
                      ! Merge the biomass of the two age classes
                      CALL merge_biomass_pfts(npts, share_rec, circ_class_n, &
                           est_circ_class_n, tmp_circ_class_n, circ_class_biomass, &
                           est_circ_class_biomass, tmp_circ_class_biomass, ipts, &
                           new_pft, old_pft, forest_managed, bm_to_litter)
                      
                      !! 2.2.4 Empty the age class that was merged and update veget_max
                      ! Guillaume M. -- /!\ WIPE THE DONOR ONLY ON A WHOLE TRANSFER. On a partial one
                      ! it keeps trees with the SAME per-m2 values, which is exactly what conserves
                      ! mass -- only its area shrank. Wiping it would leave area with NO stand.
                      ! merge_biomass_pfts writes on new_pft only, the donor enters as INTENT(in)
                      ! through tmp_circ_class_*, so conditioning this block is enough.
                      IF (.NOT. ac12_partial) THEN
                         circ_class_n(ipts,old_pft,:) = zero
                         circ_class_biomass(ipts,old_pft,:,:,:) = zero
                      ELSEIF (ok_maturity_transfer .AND. mat_shift > zero) THEN
                         ! Guillaume M. -- MAT_SHIFT counterpart: the donor steps BACKWARD by
                         ! the factors the conservation identities fixed at assembly time.
                         ! Every class keeps its rank; max_dia drops slightly at each
                         ! promotion, which restores a shred of the self-limitation the
                         ! intact copy had lost.
                         circ_class_n(ipts,old_pft,:) = &
                              circ_class_n(ipts,old_pft,:) * mat_rem_n
                         circ_class_biomass(ipts,old_pft,:,:,:) = &
                              circ_class_biomass(ipts,old_pft,:,:,:) * mat_rem_b
                      ENDIF

                   END IF ! Trusting_hack_age_class

                   ! Guillaume M. -- MATURITY_CONVEYOR: PARTIAL transfer when the promotion is
                   ! progressive. /!\ The donor now SURVIVES, so any code assuming "after a
                   ! promotion the slot is empty" is wrong. The donor keeps its PER-M2 values
                   ! unchanged and only its area shrinks, which is what conserves mass exactly,
                   ! as on the downward path.
                   IF (ac12_partial) THEN
                      ! Guillaume M. -- ac12_v_mov is computed ABOVE, before calculate_shares:
                      ! the share and the move must act on the SAME area.
                      ac12_don0 = veget_max(ipts,old_pft)
                      veget_max(ipts,new_pft) = veget_max(ipts,new_pft) + ac12_v_mov
                      veget_max(ipts,old_pft) = veget_max(ipts,old_pft) - ac12_v_mov
                      ! Guillaume M. -- Realised flux, read AFTER the write: the donor area
                      ! fraction that actually moved. The guards can clip the decided rate
                      ! to nothing; 1/AC12_F_REAL is the residence time this class realises,
                      ! the quantity checked against forest inventories.
                      IF (ac12_don0 .GT. min_stomate) THEN
                         ac12_freal_diag(ipts,old_pft) = ac12_v_mov / ac12_don0
                      ENDIF
                   ELSE
                      ! Guillaume M. -- All-or-nothing promotion: the whole donor moves, so
                      ! the realised flux is one by definition (guard: a donor already empty
                      ! moved nothing).
                      IF (veget_max(ipts,old_pft) .GT. min_stomate) THEN
                         ac12_freal_diag(ipts,old_pft) = un
                      ENDIF
                      veget_max(ipts,new_pft) = veget_max(ipts,new_pft) + &
                           veget_max(ipts,old_pft)
                      veget_max(ipts,old_pft) = zero
                   ENDIF
                   ! Guillaume M. -- /!\ DO NOT reset ac12_frac HERE: three reads still follow
                   ! (veget_max_new, age_stand_bm, keep_source) and would all conclude "whole
                   ! transfer". The reset now sits at the end of the slot, after
                   ! move_pft_properties.
                   !
                   ! Guillaume M. -- /!\ veget_max_new is a SECOND array, also emptied by the
                   ! promotion. Leaving it all-or-nothing while veget_max is fractioned would
                   ! desynchronise the two, and veget_max_new drives next year's distribution.
                   ! The downward path already does this (ph_v_mov_new).
                   IF (ac12_partial) THEN
                      ac12_v_mov_new = ac12_frac * veget_max_new(ipts,old_pft)
                      veget_max_new(ipts,new_pft) = veget_max_new(ipts,new_pft) + ac12_v_mov_new
                      veget_max_new(ipts,old_pft) = veget_max_new(ipts,old_pft) - ac12_v_mov_new
                   ELSE
                      veget_max_new(ipts,new_pft) = veget_max_new(ipts,new_pft) + &
                           veget_max_new(ipts,old_pft)
                      veget_max_new(ipts,old_pft) = zero
                   ENDIF
                   
                   ! Debug
                   IF (printlev_loc>=4) THEN
                      ! Mass balance of the merging procedure
                      temp_end(:) = zero
                      DO ipar = 1,nparts
                         DO iele = 1,nelements
                            DO icir = 1,ncirc
                               temp_end(iele) = temp_end(iele) + &
                                    SUM(circ_class_biomass(ipts,:,icir,ipar,iele) * &
                                    circ_class_n(ipts,:,icir) * veget_max(ipts,:))
                            ENDDO
                         ENDDO
                      ENDDO
                      
                      DO iele = 1,nelements
                         IF (temp_end(iele) - temp_start(iele) .LT. min_stomate .AND. &
                              temp_end(iele) - temp_start(iele) .GT. -min_stomate) THEN
                            WRITE(numout,*) 'Mass balance closure merge_biomass_pfts', iele
                         ELSE
                            WRITE(numout,*) 'Error: mass balance is not closed in merge_biomass_pfts'
                            WRITE(numout,*) '   Difference is, ',iele, temp_end(iele) - &
                                 temp_start(iele), temp_end(iele), &
                                 temp_start(iele) 
                            WRITE(numout,*) 'carbon ', &
                                 SUM(SUM(circ_class_biomass(ipts,:,:,:,icarbon),3),2)
                            WRITE(numout,*) 'nitrogen ', &
                                 SUM(SUM(circ_class_biomass(ipts,:,:,:,initrogen),3),2)
                            WRITE(numout,*) 'veget_max, ',veget_max(ipts,:)
                            IF (err_act.GT.1) THEN
                               CALL ipslerr_p(3,'ERROR:mass balance not closed',&
                                    'after merge biomass in age_class_distr','','')
                            ENDIF
                         ENDIF
                      ENDDO !nelements
                   ENDIF
                   !-

                   ! Guillaume M. -- STAND_AGE: blend of the mean age, source emptied. Done at
                   ! the caller level so as not to touch the signature of move_pft_properties
                   ! (three callers). Diagnostic only.
                   !
                   ! Guillaume M. -- The weighting is by BIOMASS, not by AREA (share_rec), as
                   ! age_stand_bm is defined. A clearcut sends class-4 area back to class 1
                   ! WITHOUT trees, and area-weighting made that bare area carry the age of the
                   ! stand that was cut. Biomass weighting settles both cases with no special
                   ! case: bare area has zero weight and the receiver keeps its age.
                   !
                   ! Guillaume M. -- /!\ The reset in stomate_lpj catches none of this: it only
                   ! acts on an EMPTY slot (veget_max <= min_stomate), and after a cut the area
                   ! is still there.
                   IF (bm_w_new + bm_w_old .GT. min_stomate) THEN
                      age_stand_bm(ipts,new_pft) = (bm_w_new * age_stand_bm(ipts,new_pft) + &
                           bm_w_old * age_stand_bm(ipts,old_pft)) / (bm_w_new + bm_w_old)
                   ELSE
                      ! Guillaume M. -- No biomass on either side: no age to carry.
                      age_stand_bm(ipts,new_pft) = zero
                   ENDIF
                   ! Guillaume M. -- AREA-weighted mirror. A BARE parcel (no standing carbon,
                   ! e.g. clearcut area sent back to class 1) carries the class-0 conveyor
                   ! residence time, n_class0 - 1 years, not the donor age: by the time that
                   ! area stands again it has sat in the conveyor for that long.
                   IF (aw_new + aw_old .GT. min_stomate) THEN
                      age_stand_area(ipts,new_pft) = (aw_new * age_stand_area(ipts,new_pft) + &
                           aw_old * MERGE(age_stand_area(ipts,old_pft), &
                                          MAX(REAL(n_class0,r_std) - un, zero), &
                                          bm_w_old .GT. min_stomate)) / (aw_new + aw_old)
                   ELSE
                      age_stand_area(ipts,new_pft) = zero
                   ENDIF
                   ! Guillaume M. -- MAT_AGE_ENTRY: on an UPWARD move the parcel's age is the age at
                   ! which area reaches the bound of the receiving class, a cohort quantity. Leaky
                   ! mean per receiving slot; a bare parcel carries the class-0 residence, as above.
                   ! Downward moves (clearcut area back to class 1) are not entries into maturity.
                   IF (ALLOCATED(mat_age_entry) .AND. new_pft .GT. old_pft .AND. &
                        aw_old .GT. min_stomate) THEN
                      mae_age = MERGE(age_stand_area(ipts,old_pft), &
                           MAX(REAL(n_class0,r_std) - un, zero), bm_w_old .GT. min_stomate)
                      mae_decay = EXP(- un / MAX(mat_age_entry_tau, un))
                      IF (mat_age_entry(ipts,new_pft) .GT. zero) THEN
                         mat_age_entry(ipts,new_pft) = mae_decay * mat_age_entry(ipts,new_pft) + &
                              (un - mae_decay) * mae_age
                      ELSE
                         mat_age_entry(ipts,new_pft) = mae_age
                      ENDIF
                   ENDIF
                   ! Guillaume M. -- /!\ Same rule as for circ_class_*: on a PARTIAL transfer the
                   ! donor survives and has not become younger, so its age is left unchanged.
                   ! The downward path already carries this reasoning; the upward one never
                   ! needed it before, having never had a survivor.
                   IF (.NOT. ac12_partial) THEN
                      age_stand_bm(ipts,old_pft) = zero
                      age_stand_area(ipts,old_pft) = zero
                   ENDIF

                   !! 2.3 Calculate the PFT characteristics of the merged PFT
                   !  Take the weighted mean of the existing vegetation and the new
                   !  vegetation in this PFT.
                   CALL move_pft_properties(ipts, new_pft, old_pft, share_rec, &
                        litter_weight_rec, circ_class_biomass, fm_change_map, &
                        wstress_season(ipts,new_pft), wstress_season(ipts,old_pft), &
                        wstress_month(ipts,new_pft), wstress_month(ipts,old_pft), &
                        season_drought_legacy(ipts,new_pft,:), season_drought_legacy(ipts,old_pft,:), &
                        nstress_season(ipts,new_pft), nstress_season(ipts,old_pft), &
                        n_input(ipts,new_pft,:,:), n_input(ipts,old_pft,:,:), &
                        n_input_daily(ipts,new_pft,:), n_input_daily(ipts,old_pft,:), & 
                        vegstress_season(ipts,new_pft), vegstress_season(ipts,old_pft), &
                        vegstress_month(ipts,new_pft), vegstress_month(ipts,old_pft), &
                        vegstress_week(ipts,new_pft), vegstress_week(ipts,old_pft), &
                        vegstress(ipts,new_pft), vegstress(ipts,old_pft), &
                        humrel(ipts,new_pft), humrel(ipts,old_pft), &
                        everywhere(ipts,new_pft), everywhere(ipts,old_pft), &
                        PFTpresent(ipts,new_pft), PFTpresent(ipts,old_pft), &
                        sugar_load(ipts,new_pft), sugar_load(ipts,old_pft), &
                        cn_leaf_min_season(ipts,new_pft), cn_leaf_min_season(ipts,old_pft), &
                        cn_leaf_init_2D, &
                        Light_Abs_Tot(ipts,new_pft,:), Light_Abs_Tot(ipts,old_pft,:), &
                        Light_Tran_Tot(ipts,new_pft,:), Light_Tran_Tot(ipts,old_pft,:), &
                        laieff_isotrop(ipts,:,new_pft), laieff_isotrop(ipts,:,old_pft), &
                        lm_lastyearmax(ipts,new_pft), lm_lastyearmax(ipts,old_pft), &
                        lm_thisyearmax(ipts,new_pft), lm_thisyearmax(ipts,old_pft), &
                        age(ipts,new_pft), age(ipts,old_pft), &
                        leaf_frac(ipts,new_pft,:), leaf_frac(ipts,old_pft,:), &
                        leaf_age(ipts,new_pft,:), leaf_age(ipts,old_pft,:), &
                        light_tran_to_floor_season(ipts,new_pft), &
                        light_tran_to_floor_season(ipts,old_pft), &
                        qsintveg(ipts,new_pft), qsintveg(ipts,old_pft), &
                        us(ipts,new_pft,:,:), us(ipts,old_pft,:,:), &
                        when_growthinit(ipts,new_pft), when_growthinit(ipts,old_pft), &
                        gdd_from_growthinit(ipts,new_pft), gdd_from_growthinit(ipts,old_pft), &
                        gdd_midwinter(ipts,new_pft), gdd_midwinter(ipts,old_pft), &
                        time_hum_min(ipts,new_pft), time_hum_min(ipts,old_pft), &
                        hum_min_dormance(ipts,new_pft), hum_min_dormance(ipts,old_pft), &
                        gdd_m5_dormance(ipts,new_pft), gdd_m5_dormance(ipts,old_pft), &
                        ncd_dormance(ipts,new_pft), ncd_dormance(ipts,old_pft), & 
                        ngd_minus5(ipts,new_pft), ngd_minus5(ipts,old_pft), &
                        mean_start_gs(ipts,new_pft), mean_start_gs(ipts,old_pft), &
                        maxvegstress_lastyear(ipts,new_pft), maxvegstress_lastyear(ipts,old_pft), &
                        maxvegstress_thisyear(ipts,new_pft), maxvegstress_thisyear(ipts,old_pft), &
                        minvegstress_lastyear(ipts,new_pft), minvegstress_lastyear(ipts,old_pft), &
                        minvegstress_thisyear(ipts,new_pft), minvegstress_thisyear(ipts,old_pft), &
                        maxgppweek_lastyear(ipts,new_pft), maxgppweek_lastyear(ipts,old_pft), & 
                        maxgppweek_thisyear(ipts,new_pft), maxgppweek_thisyear(ipts,old_pft), &
                        n_reserve_longterm(ipts,new_pft), n_reserve_longterm(ipts,old_pft), &
                        n_reserve_balance(ipts,new_pft), n_reserve_balance(ipts,old_pft), &
                        maxfpc_lastyear(ipts,new_pft), maxfpc_lastyear(ipts,old_pft), &
                        maxfpc_thisyear(ipts,new_pft), maxfpc_thisyear(ipts,old_pft), &
                        turnover_longterm(ipts,new_pft,:,:), turnover_longterm(ipts,old_pft,:,:), &
                        dead_leaves(ipts,new_pft,:), dead_leaves(ipts,old_pft,:), &
                        grow_season_len(ipts,new_pft), grow_season_len(ipts,old_pft), &
                        KF(ipts,new_pft), KF(ipts,old_pft), &
                        atm_to_bm(ipts,new_pft,:), atm_to_bm(ipts,old_pft,:), &
                        npp_longterm(ipts,new_pft), npp_longterm(ipts,old_pft), &
                        croot_longterm(ipts,new_pft), croot_longterm(ipts,old_pft), &
                        gpp_daily(ipts,new_pft), gpp_daily(ipts,old_pft), &
                        resp_maint(ipts,new_pft), resp_maint(ipts,old_pft), &
                        resp_growth(ipts,new_pft), resp_growth(ipts,old_pft), &
                        resp_hetero(ipts,new_pft), resp_hetero(ipts,old_pft), &
                        npp_daily(ipts,new_pft), npp_daily(ipts,old_pft), &
                        rue_longterm(ipts,new_pft), rue_longterm(ipts,old_pft), &
                        leaching_daily(ipts,new_pft,:), leaching_daily(ipts,old_pft,:), &
                        emission_daily(ipts,new_pft,:), emission_daily(ipts,old_pft,:), &
                        gpp_week(ipts,new_pft), gpp_week(ipts,old_pft), &
                        resp_maint_week(ipts,new_pft), resp_maint_week(ipts,old_pft), &
                        gpp_year(ipts,new_pft), gpp_year(ipts,old_pft), &
                        gpp_decade(ipts,new_pft), gpp_decade(ipts,old_pft), &
                        mai(ipts,new_pft), mai(ipts,old_pft), &
                        pai(ipts,new_pft), pai(ipts,old_pft), &
                        mai_count(ipts,new_pft), mai_count(ipts,old_pft), &
                        previous_wood_volume(ipts,new_pft), previous_wood_volume(ipts,old_pft), &
                        age_stand(ipts,new_pft), age_stand(ipts,old_pft), &
                        last_cut(ipts,new_pft), last_cut(ipts,old_pft), &
                        k_latosa_adapt(ipts,new_pft), k_latosa_adapt(ipts,old_pft), &
                        litter(ipts,:,new_pft,:,:), litter(ipts,:,old_pft,:,:), &
                        bm_to_litter(ipts,new_pft,:,:), bm_to_litter(ipts,old_pft,:,:), &
                        tree_bm_to_litter(ipts,new_pft,:,:), tree_bm_to_litter(ipts,old_pft,:,:), &
                        leaf_age_crit(ipts,new_pft), leaf_age_crit(ipts,old_pft), &
                        leaf_classes(ipts,new_pft), leaf_classes(ipts,old_pft), &
                        co2_fire(ipts,new_pft), co2_fire(ipts,old_pft), &
                        litterfuel(ipts,:,new_pft,:,:), litterfuel(ipts,:,old_pft,:,:), &
                        turnover_daily(ipts,new_pft,:,:), turnover_daily(ipts,old_pft,:,:), &
                        soil_n_min(ipts,new_pft,:), soil_n_min(ipts,old_pft,:), &
                        p_O2(ipts,new_pft), p_O2(ipts,old_pft), &
                        bact(ipts,new_pft), bact(ipts,old_pft), &
                        deepSOM_a(ipts,:,new_pft,:), deepSOM_a(ipts,:,old_pft,:), &
                        deepSOM_s(ipts,:,new_pft,:), deepSOM_s(ipts,:,old_pft,:), &
                        deepSOM_p(ipts,:,new_pft,:), deepSOM_p(ipts,:,old_pft,:), &
                        som(ipts,:,new_pft,:), som(ipts,:,old_pft,:), &
                        lignin_struc(ipts,new_pft,:), lignin_struc(ipts,old_pft,:), &
                        lignin_wood(ipts,new_pft,:), lignin_wood(ipts,old_pft,:), &
                        lignin_snag(ipts,new_pft,:), lignin_snag(ipts,old_pft,:), &
                        forest_managed(ipts,new_pft), forest_managed(ipts,old_pft), &
                        plant_status(ipts,new_pft), plant_status(ipts,old_pft), &
                        kill_vessels(ipts,new_pft,:), kill_vessels(ipts,old_pft,:), &
                        vessel_loss_previous(ipts,new_pft,:), vessel_loss_previous(ipts,old_pft,:), &
                        biomass_init_drought(ipts,new_pft,:,:,:), biomass_init_drought(ipts,old_pft,:,:,:), &
                        count_daylight(ipts,new_pft), count_daylight(ipts,old_pft), &
                        CN_som_litter_longterm(ipts,new_pft,:), CN_som_litter_longterm(ipts,old_pft,:), &
                        matrixV(ipts,new_pft,:,:), matrixV(ipts,old_pft,:,:), &
                        matrixA(ipts,new_pft,:,:), matrixA(ipts,old_pft,:,:), &
                        vectorU(ipts,new_pft,:), vectorU(ipts,old_pft,:), &
                        vectorB(ipts,new_pft,:), vectorB(ipts,old_pft,:), &
                        gap_area_save(ipts,new_pft,:), gap_area_save(ipts,old_pft,:), &
                        total_ba_init(ipts,new_pft), total_ba_init(ipts,old_pft), &
                        ! Guillaume M. -- MATURITY_CONVEYOR: on a PARTIAL transfer the source keeps
                        ! its INTENSIVE values (per m2) while its area shrinks, which is what
                        ! conserves mass when a slot is split. Without keep_source,
                        ! move_pft_properties zeroes the source and transfers litter, soil and
                        ! lignin IN FULL while only ac12_frac of the area moves.
                        !
                        ! Guillaume M. -- /!\ Scanning for '= zero' assignments CANNOT catch this:
                        ! those resets live INSIDE the called routine. Look for calls whose
                        ! default behaviour is destructive, not only for assignments.
                        keep_source=ac12_partial)

                   ! Guillaume M. -- END OF SLOT: ac12_frac has served its last read
                   ! (keep_source above). It is consumed ONLY now, so that it cannot trigger a
                   ! partial transfer on the next slot.
                   ac12_frac = zero
                   ac12_partial = .FALSE.

                ELSE

                   ! Debug
                   IF (printlev_loc>=4) WRITE(numout,*) 'old and new pft are the same'
                   !-

                ENDIF
                
             ENDDO ! loop over nagec_pft(ivma)

          ELSE ! Grasses and croplands

             ! Don't have to do anything

             ! Debug
             IF (printlev_loc>=4) WRITE(numout,*) 'Grass or cropland ',&
                  'no age classes', ivm
             !-

          ENDIF ! is_tree(ivm)

       ENDDO ! Looping over real PFTs

    ENDDO ! loop over #pixels - domain size

    ! Guillaume M. -- PROGRESSIVE_HARVEST: make the degradation into a thinning visible.
    ! Without this trace, a harvest refused on most slots (slot areas too small against
    ! min_vegfrac) would look active while moving no area at all.
    IF (ok_progressive_harvest .AND. ph_nskip .GT. 0) THEN
       WRITE(numout,*) '[PROGRESSIVE_HARVEST] decoupages refuses (invariant min_vegfrac), ', &
            'restes en eclaircie : ', ph_nskip
    ENDIF


!! 3. Check numerical consistency of this routine

    IF (err_act.GT.1) THEN

       ! 4.2 Check surface area
       CALL check_vegetation_area("age class distr", npts, veget_max_begin, &
            veget_max,'ageclass')

       ! 4.3 Mass balance closure 
       pool_end = zero
       DO iele = 1,nelements

          ! atm_to_bm
          pool_end(:,:,iele) = pool_end(:,:,iele) + &
               atm_to_bm(:,:,iele)

          ! Biomass pool + bm_to_litter
          DO ipar = 1,nparts

             DO icir = 1,ncirc

                pool_end(:,:,iele) = pool_end(:,:,iele) + &
                     circ_class_biomass(:,:,icir,ipar,iele) * &
                     circ_class_n(:,:,icir) * veget_max(:,:)

             ENDDO

             ! Add bm_to_litter to the biomass pool
             pool_end(:,:,iele) = pool_end(:,:,iele) + &
                  ((bm_to_litter(:,:,ipar,iele)+ &
                  turnover_daily(:,:,ipar,iele)) * veget_max(:,:))

          ENDDO
          
          ! Litter pool (gC m-2) *  (m2 m-2) 
          DO ilit = 1,nlitt
             DO ilev = 1,nlevs
                pool_end(:,:,iele) = pool_end(:,:,iele) + &
                     litter(:,ilit,:,ilev,iele) * veget_max(:,:)
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

       ENDDO

       ! Guillaume M. -- Exact mirror of the pool_start fix: these two terms do not depend on
       ! iele and were counted twice. The two blocks must stay symmetric, otherwise the
       ! closure would break for good.
       ! Denitrifying bacteria (only C)
       pool_end(:,:,icarbon) = pool_end(:,:,icarbon) + &
               bact(:,:) * veget_max(:,:)

       ! Nitrogen only
       DO ispec = 1,nnspec
          pool_end(:,:,initrogen) = pool_end(:,:,initrogen) + &
               soil_n_min(:,:,ispec) * veget_max(:,:)
       ENDDO
       
       ! 4.3.2 Calculate mass balance
       ! Common processes
       DO iele=1,nelements
          check_intern(:,:,ipoolchange,iele) = -un * (pool_end(:,:,iele) - &
               pool_start(:,:,iele)) 
       ENDDO

       closure_intern(:,:,:) = zero
       DO imbc = 1,nmbcomp
          DO iele=1,nelements
             ! Debug
             IF (printlev_loc>=4) &
                   WRITE(numout,*) &
                  'check_intern, ivm, imbc, iele, ', imbc, &
                  iele, check_intern(:,test_pft,imbc,iele)
             !-
             closure_intern(:,:,iele) = closure_intern(:,:,iele) + &
                  check_intern(:,:,imbc,iele)
          ENDDO
       ENDDO

       ! Guillaume M. -- INSTRUMENTATION: the standard message only gives a TOTAL per
       ! element. check_intern already holds the per-component breakdown of the balance, so
       ! it is printed rather than guessing which reservoir is drifting.
       DO ipts = 1,npts
          DO ivma = 1,nvmap
             IF (.NOT. is_tree(start_index(ivma))) CYCLE
             ac0_sum = zero
             DO iagec = 1,nagec_pft(ivma)
                ac0_sum = ac0_sum + closure_intern(ipts,start_index(ivma)+iagec-1,icarbon)
             ENDDO
             IF (ABS(ac0_sum) > min_stomate) THEN
                WRITE(numout,*) '=== FERMETURE ROMPUE, ipts,ivma,ecart :',ipts,ivma,ac0_sum
                DO iagec = 1,nagec_pft(ivma)
                   WRITE(numout,*) '   pft',start_index(ivma)+iagec-1, &
                        ' vegetmax',veget_max(ipts,start_index(ivma)+iagec-1), &
                        ' comp(atm2land,land2atm,lat2out,lat2in,poolchange)', &
                        check_intern(ipts,start_index(ivma)+iagec-1,:,icarbon)
                ENDDO
                ! Guillaume M. -- /!\ If this closure breaks again, the per-RESERVOIR breakdown
                ! is what finds the cause: six terms (atm_to_bm, standing biomass,
                ! bm_to_litter+turnover, litter, soil, bact) whose sum must equal -ac0_sum.
                ! That completeness check is what revealed the double counting of bact. It was
                ! removed once its exit criterion was met.
                CALL flush(numout)
             ENDIF
          ENDDO
       ENDDO

       ! 4.3.3 Check mass balance closure
       CALL check_mass_balance("age class distr", &
            closure_intern, npts, pool_end, pool_start, veget_max, &
            'ageclass')
      
    ENDIF ! err_act.GT.1

    IF(printlev.GE.4) WRITE(numout,*) 'Leaving age class distribution'
    
    WRITE(numout,*) '[C0IN] aire inscrite en classe 0 cette annee :', c0in_sum
    ! Guillaume M. -- [C0BACK]: is the booked area still carried by class 1 on exit?
    c0b_book = zero ; c0b_carry = zero ; c0b_excess = zero
    DO c0b_i = 1,npts
       DO c0b_k = 1,nvmap
          IF (.NOT. is_tree(start_index(c0b_k))) CYCLE
          c0b_book  = c0b_book  + class0_veg(c0b_i,c0b_k,1)
          c0b_carry = c0b_carry + veget_max(c0b_i,start_index(c0b_k))
          c0b_excess = c0b_excess + &
               MAX(zero, class0_veg(c0b_i,c0b_k,1) - veget_max(c0b_i,start_index(c0b_k)))
       ENDDO
    ENDDO
    WRITE(numout,*) '[C0BACK] inscrite', c0b_book, ' portee par classe 1', c0b_carry, &
         ' excedent sans support', c0b_excess
    CALL flush(numout)

    ! Guillaume M. -- Annual report of the sub-threshold classes routed to litter by
    ! merge_biomass_pfts, then reset. A silent zero here means the path was NOT taken,
    ! which is what tells a working fix apart from an untested one.
    WRITE(numout,*) '[ORPHLIT] classes sous seuil versees en litiere :', orphan_litter_hits, &
         ' carbone (gC m-2 cumule sur les creneaux) :', orphan_litter_mass
    CALL flush(numout)
    orphan_litter_hits = 0
    orphan_litter_mass = zero

    CALL check_area_invariant('exit age_class_distr', npts, veget_max, circ_class_n)

    !! ROTATION_GROWTH — design/MODULE_DESIGN_ROTATION_GROWTH.md
    !! Guillaume M. -- Placed HERE, and not in set_management_intensity, for two reasons:
    !! max_dia is a local of this routine, and this block sits under ts_annual_proc while
    !! set_management_intensity runs at the daily cadence. Putting it there would have
    !! integrated 365 times a year -- the trap that once made N_CLASS0 = 4 mean 4 DAYS.
    IF (rotation_growth_weight > zero .AND. ALLOCATED(rotation_growth)) THEN
       CALL update_rotation_growth(npts)
    ENDIF

  END SUBROUTINE age_class_distr


! ================================================================================================================================
!! SUBROUTINE   : update_rotation_growth
!!
!>\BRIEF        Rotation that the STAND can deliver, from the ages at which its area entered the classes.
!!
!! DESCRIPTION  : Guillaume M. -- Annual. The class bounds are diameters; mat_age_entry holds the age
!!                at which area reached each of them, so the last measured leg gives a cohort growth
!!                g = (D_{n-1} - D_{n-2}) / (a_n - a_{n-1}) and the rotation is the age at the cut
!!                diameter, R = a_n + (D_cut - D_{n-1}) / g. No diameter velocity is read on a slot:
!!                the earlier estimators (slot max_dia, stable-area filter, circumference bins) all
!!                failed because a slot is not a cohort. See design/MODULE_DESIGN_MAT_AGE_ENTRY.md.
!!
!! MAIN OUTPUT VARIABLE(S): rotation_growth
!!_ ================================================================================================================================

  SUBROUTINE update_rotation_growth(npts)

    IMPLICIT NONE

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                       :: npts            !! Number of grid cells

    !! 0.4 Local variables
    INTEGER(i_std)                                   :: ipts, ivma, iagec, ipft, n, kn, kn1
    REAL(r_std)                                      :: dfac            !! Cut-diameter modulation of the cell (-)
    REAL(r_std)                                      :: a_n, a_n1       !! Entry ages of the last two classes (yr)
    REAL(r_std)                                      :: d_n1, d_n2      !! Diameter bounds entered at those ages (m)
    REAL(r_std)                                      :: d_cut           !! Cut diameter (m)
    REAL(r_std)                                      :: g               !! Cohort diameter growth on the last leg (m/yr)
    REAL(r_std)                                      :: r_est           !! Growth-derived rotation (yr)

!_ ================================================================================================================================

    IF (nagec <= 1) RETURN
    IF (.NOT. ALLOCATED(mat_age_entry)) RETURN

    ! /!\ No global reset: rotation_growth is persisted and a group without a sample this
    ! year keeps its last valid value, otherwise the blend would flicker year to year.
    DO ipts = 1, npts
       dfac = un
       IF (ALLOCATED(dia_factor)) dfac = dia_factor(ipts)
       DO ivma = 1, nvmap
          ipft = start_index(ivma)
          IF (.NOT. is_tree(ipft)) CYCLE
          n = nagec_pft(ivma)
          IF (n <= 1) CYCLE
          IF (largest_tree_dia(ipft) <= zero) CYCLE
          kn  = ipft + n - 1
          kn1 = ipft + n - 2
          a_n = mat_age_entry(ipts,kn)
          IF (a_n <= zero) CYCLE                       ! nothing has reached the terminal class yet
          d_n1  = age_class_bound(n-1, kn) * dfac
          d_cut = largest_tree_dia(kn) * dfac * (un - dia_rotation_tol)
          ! Guillaume M. -- Growth on the last measured leg (entry into class n-1 -> entry into
          ! class n). Fallback: mean growth from age 0 when class n-1 has no sample (n = 2 or
          ! not yet fed). Both are cohort quantities: ages of moved area, not slot velocities.
          g = zero
          IF (n >= 3) THEN
             a_n1 = mat_age_entry(ipts,kn1)
             d_n2 = age_class_bound(n-2, kn1) * dfac
             IF (a_n1 > zero .AND. a_n > a_n1 .AND. d_n1 > d_n2) g = (d_n1 - d_n2) / (a_n - a_n1)
          ENDIF
          IF (g <= zero .AND. d_n1 > zero) g = d_n1 / a_n
          IF (g <= zero) CYCLE
          r_est = a_n + MAX(zero, d_cut - d_n1) / g
          r_est = MIN(rotation_growth_max, MAX(rotation_growth_min, r_est))
          DO iagec = 1, n
             rotation_growth(ipts, ipft + iagec - 1) = r_est
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE update_rotation_growth

! ================================================================================================================================
!! SUBROUTINE   : merge_biomass_pfts
!!
!>\BRIEF        Redistribute biomass, litter, soilcarbon and water across
!!              the age classes
!!
!! DESCRIPTION  : Following growth, the trees from an age class may have too
!! big to still belong to this age class. The biomass, litter, soilcarbon and
!! soil water then need to be moved from one age class to the next age class.
!!
!! RECENT CHANGE(S) : 
!!
!! MAIN OUTPUT VARIABLE(S) :  
!!
!! REFERENCES   : None
!!
!! FLOWCHART    : 
!! \n
!_ ================================================================================================================================

  SUBROUTINE merge_biomass_pfts(npts, share_rec, circ_class_n, est_circ_class_n, &
       tmp_circ_class_n, circ_class_biomass, est_circ_class_biomass, tmp_circ_class_biomass, &
       ipts, new_pft, old_pft, forest_managed, bm_to_litter)

    IMPLICIT NONE
    
  !! 0. Variable and parameter declaration 
    
    !! 0.1 Input variables
    
    INTEGER, INTENT(in)                              :: npts                  !! Domain size - number of pixels (unitless)
    REAL(r_std), INTENT(in)                          :: share_rec             !! Share of the veget_max of the existing vegetation
                                                                              !! within a PFT over the total veget_max following 
                                                                              !! expansion of that PFT (unitless, 0-1)

    REAL(r_std), DIMENSION(:,:,:), INTENT(in)        :: est_circ_class_biomass!! Temporary variable for circ_class_biomass 
                                                                              !! of the established vegetation
                                                                              !! @tex $(g C ind^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)        :: tmp_circ_class_biomass!! Temporary variable for circ_class_biomass 
                                                                              !! of the new vegetation 
                                                                              !! @tex $(g C ind^{-1})$ @endtex
    REAL(r_std), DIMENSION(:), INTENT(in)            :: est_circ_class_n      !! Temporary variable for circ_class_n
                                                                              !! of the established vegetation
                                                                              !! @tex $(ind m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:), INTENT(in)            :: tmp_circ_class_n      !! Temporary variable for circ_class_n
                                                                              !! of the new vegetation
                                                                              !! @tex $(ind m^{-2})$ @endtex
    INTEGER(i_std), INTENT(in)                       :: ipts, new_pft, old_pft!! Indices for point and PFT (unitless)
    INTEGER(i_std), DIMENSION(:,:), INTENT(in)       :: forest_managed        !! forest management flag

    !! 0.2 Output variables

    !! 0.3 Modified variables
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout) :: circ_class_biomass    !! circ_class_biomass of the merged vegetation 
                                                                              !! @tex $(g C ind^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)     :: circ_class_n          !! circ_class_n of the merged vegetation
                                                                              !! @tex $(ind m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)   :: bm_to_litter          !! Transfer of biomass to litter
                                                                              !! @tex $(g C m^{-2})$ @endtex

    !! 0.4 Local variables
    INTEGER(i_std)                                   :: ipart                 !! Indices (unitless)
    INTEGER(i_std)                                   :: icir, jcir, iele      !! Indices (unitless)
    REAL(r_std)                                      :: trees_needed          !! Number of trees leftthat should be present in a 
                                                                              !! circumference class when merging the established
                                                                              !! and newly established vegetation
                                                                              !! @tex $(ind m^{-2})$ @endtex 
    REAL(r_std)                                      :: tempa, tempb          !! Variable to store a temporary value (unit varies)
    REAL(r_std), DIMENSION(nelements)                :: orphan_to_litter      !! Mass of the sub-threshold classes routed to litter
                                                                              !! instead of being discarded @tex $(g C m^{-2})$ @endtex
    REAL(r_std), DIMENSION(nparts,nelements)         :: tempc                 !! Variable to store a temporary value of
                                                                              !! circ_class_biomass  
    REAL(r_std), DIMENSION(2*ncirc)                  :: trees_left            !! Number of trees left in a circumference class when
                                                                              !! merging the available and newly established 
                                                                              !! vegetation @tex $(ind m^{-2})$ @endtex
    REAL(r_std), DIMENSION(2*ncirc)                  :: temp_circ_class       !! Temporary variable for circ_class of the expanded
                                                                              !! PFT for sorting the availble and newly established
                                                                              !! circumference classes (m)
    REAL(r_std), DIMENSION(2*ncirc)                  :: temp_circ_class_n     !! Temporary variable for circ_class_n of the expanded
                                                                              !! PFT for sorting the availble and newly established
                                                                              !! circumference classes @tex $(ind m^{-2})$ @endtex
    REAL(r_std), DIMENSION(2*ncirc)                  :: grp_n                 !! MATURITY_TRANSFER: stems per merge group (ind m-2)
    REAL(r_std), DIMENSION(2*ncirc)                  :: grp_nd2               !! MATURITY_TRANSFER: sum of n*dia^2 per merge group
                                                                              !! (ind m-2 m2), carries the group quadratic diameter
    INTEGER(i_std)                                   :: nact                  !! MATURITY_TRANSFER: number of active merge groups
    INTEGER(i_std)                                   :: kmin                  !! MATURITY_TRANSFER: index of the closest pair
    REAL(r_std)                                      :: dmin, dist_qmd        !! MATURITY_TRANSFER: relative QMD gaps (unitless)
    REAL(r_std)                                      :: qmd_a, qmd_b          !! MATURITY_TRANSFER: group quadratic diameters (m)
    REAL(r_std), DIMENSION(2*ncirc,nparts,nelements) :: temp_circ_class_biomass !! Temporary variable for sorting circ_class_biomass
                                                                                !! @tex $(g C ind^{-1})$ @endtex
    REAL(r_std), DIMENSION(nelements)                :: total_biomass_begin     !! total biomass at the beginning of the merge (gC or gN m-2)
    REAL(r_std), DIMENSION(nelements)                :: total_biomass_merge     !! total biomass after the merge (gC or gN m-2)
    REAL(r_std)                                      :: mismatch                !! Mismatch between the merged and the initial biomass (gC or gN m-2)
    LOGICAL                                          :: fixed                   !! Keep track of whether the mass imbalance got fixed or not
    REAL(r_std), DIMENSION(ncirc)                    :: est_ba                  !! temporary variable to calculate ba of established vegetation
    REAL(r_std), DIMENSION(ncirc)                    :: tmp_ba                  !! temporary variable to calculate ba of new vegetation
    REAL(r_std)                                      :: total_ba                !! Basal area of the merged stand
    REAL(r_std)                                      :: total_n                 !! Stem number of the merged stand
    REAL(r_std)                                      :: qmd                     !! qmd of the merged stand
    REAL(r_std), DIMENSION(ncirc)                    :: est_circ_class_dia      !! Variable to store temporary values for circ_class (m)
    REAL(r_std), DIMENSION(ncirc)                    :: tmp_circ_class_dia      !! Variable to store temporary values for circ_class (m)
    REAL(r_std), DIMENSION(ncirc)                    :: exp_circ_class_n        !! Temporary variable for circ_class_n of the expanded
                                                                                !! PFT @tex $(ind m^{-2})$ @endtex

!_ ================================================================================================================================

 !! 1. Initialize

    IF (printlev >= 3) WRITE(numout,*) 'Entering merging biomass of pfts'
    
    ! Store the initial biomass for each component and each element
    total_biomass_begin(:) = zero
    total_biomass_merge(:) = zero
    DO icir = 1,ncirc
       total_biomass_begin(:) = total_biomass_begin(:) + &
            share_rec * SUM(est_circ_class_biomass(icir,:,:),1) * &
            est_circ_class_n(icir) + (1-share_rec) * &
            SUM(tmp_circ_class_biomass(icir,:,:),1) * tmp_circ_class_n(icir)
    END DO

 !! 2. Merge biomass, litter, soilcarbon and soil water

    !  We have already ::ncirc classes in this PFT and we add another 
    !  ::ncirc through calling stomate_prescribe. What we want to 
    !  do is to merge these 2*ncirc into ncirc circumference classes.     
    !  We need to populate the new distribution of trees among
    !  the circumference classes, and move the heartwood and sap masses 
    !  to the new distributions. We also need to move the other pools.  
    !  This is more tricky because the allometric relations for the 
    !  new tree sizes will NOT give the same amount of biomass for
    !  the non-woody pools. Let us first try it just distributing
    !  everything equally, and then if that doesn't work we can
    !  redistribute the non-woody pools more cleverly, even if that 
    !  means we change the total amount of biomass we have in this
    !  stand.  I want to try it this way because we conserve biomass 
    !  this way, and I hope that the stress on the trees will be
    !  small enough that it doesn't cause problems. Redistribution
    !  is not expected to de-stabilize the allocation too much if the
    !  diameter range within an age-class is not too big (or in other
    !  words if we define a sufficiently large number of age classes).
 
    ! Store circumference (m) in a temporary variable that can be
    ! sorted from small to large - note that the dimension of
    ! these variables are 2*ncirc
    est_circ_class_dia(:) = wood_to_dia(est_circ_class_biomass(:,:,icarbon),&
         new_pft,pipe_tune2(ipts,new_pft))
    tmp_circ_class_dia(:) = wood_to_dia(tmp_circ_class_biomass(:,:,icarbon),&
         old_pft,pipe_tune2(ipts,old_pft))

    ! Calculate qmd for the merged PFTs
    est_ba(:) = wood_to_ba(est_circ_class_biomass(:,:,icarbon), &
         new_pft, pipe_tune2(ipts,new_pft))
    tmp_ba(:) = wood_to_ba(tmp_circ_class_biomass(:,:,icarbon), &
         old_pft, pipe_tune2(ipts,old_pft))
    total_ba = share_rec * SUM(est_circ_class_n(:) * est_ba(:)) + &
         (un - share_rec) * SUM(tmp_circ_class_n(:) * tmp_ba(:))
    ! Guillaume M. -- Quadratic mean diameter from the stand identity BA = N * pi/4 * D**2,
    ! hence D = SQRT(4*BA/(pi*N)). wood_to_ba returns the basal area of ONE tree, so total_ba
    ! is a stand total and the stem number must divide it. Without total_n, and with pi/4
    ! where 4/pi belongs, qmd came out four times too small; it seeds weibull_class_dist, so
    ! every merge inherited the error.
    total_n = share_rec * SUM(est_circ_class_n(:)) + &
         (un - share_rec) * SUM(tmp_circ_class_n(:))
    IF (total_n > min_stomate .AND. total_ba > min_stomate) THEN
       qmd = SQRT(quatre * total_ba / (pi * total_n))
    ELSE
       qmd = zero
    ENDIF

    ! Guillaume M. -- MATURITY_TRANSFER: no law is drawn here. The target counts are
    ! built AFTER the sort (section 2.3bis) by partitioning the sorted cohorts, so the
    ! receiver structure is history, not a re-imposed Weibull. Skipping the call also
    ! avoids its fatal guard on empty classes, which a real stand may legitimately have.
    IF (.NOT. ok_maturity_transfer) THEN

       ! Calculate diameter distribution making use of qmd
       exp_circ_class_n(:) = weibull_class_dist(qmd * m_to_cm, new_pft, &
            forest_managed(ipts, new_pft))


       !! 2.1 Calculate new distribution
       !  The distribution is stored in exp_circ_class_n and determined by
       !  the Weibull distribution. Below we multiply the distribution
       !  with the total number of individual to obtain circ_class_n.
       IF (share_rec .LT. 10*EPSILON(un)) THEN

          ! The new PFT is still empty so simply copy the current tree
          ! distribution into the new PFT.
          exp_circ_class_n(:) = tmp_circ_class_n(:)

       ELSE

          exp_circ_class_n(:) = exp_circ_class_n(:) * &
               SUM(share_rec * est_circ_class_n(:) + &
               (un-share_rec) * tmp_circ_class_n(:))

       ENDIF

    ENDIF


    !! 2.2 Copy the biomass of the established and the new into a single variable
    DO icir = 1,2*ncirc

       IF (icir .LE. ncirc) THEN

          ! Circumference classes (m), number of individuals (#/m-2) and biomass
          ! of the available vegetation (gC or gN tree-1). Weight by the area share
          ! such that we don't need to do this further down. Don't weight the 
          ! circumference and biomass as those variables are expressed per tree. 
          ! This helps to simplify the equations.
          temp_circ_class(icir) = pi * est_circ_class_dia(icir)
          temp_circ_class_n(icir) = share_rec * &
               est_circ_class_n(icir)
          temp_circ_class_biomass(icir,:,:) = &
               est_circ_class_biomass(icir,:,:)  

       ELSEIF ( (icir .GT. ncirc) .AND. (icir .LE. 2*ncirc) ) THEN

          ! Circumference classes, number of individuals and biomass
          ! of the newly established vegetation (m). Weight by the area 
          ! share such that we don't need to do this further down. Don't
          ! weight the circumference and biomass as those variables are
          ! expressed per tree. This helps to simplify the equations.
          temp_circ_class(icir) = pi * tmp_circ_class_dia(icir-ncirc)
          temp_circ_class_n(icir) = (un - share_rec) * &
               tmp_circ_class_n(icir-ncirc)
          temp_circ_class_biomass(icir,:,:) = &
               tmp_circ_class_biomass(icir-ncirc,:,:)

       ELSE

          WRITE(numout,*) 'Error - problem with the dimensions of a do-loop '

       ENDIF

    ENDDO !loop over 2*ncirc

    !! 2.3 Sort the biomass and number of trees according to their 
    !  circumference 
    DO jcir = 2,2*ncirc

       tempa = temp_circ_class(jcir)
       tempb = temp_circ_class_n(jcir)
       tempc(:,:) = temp_circ_class_biomass(jcir,:,:)

       DO icir = jcir-1,1,-1

          IF (temp_circ_class(icir) .LE. tempa) EXIT
          temp_circ_class(icir+1) = temp_circ_class(icir)
          temp_circ_class_n(icir+1) = temp_circ_class_n(icir)
          temp_circ_class_biomass(icir+1,:,:) = &
               temp_circ_class_biomass(icir,:,:)

       ENDDO

       temp_circ_class(icir+1) = tempa
       temp_circ_class_n(icir+1) = tempb
       temp_circ_class_biomass(icir+1,:,:) = tempc(:,:)

    ENDDO
  
    !! 2.3bis MATURITY_TRANSFER: partition the sorted cohorts by diameter proximity
    ! Guillaume M. -- The 2*ncirc sorted cohorts are cut into ncirc CONTIGUOUS groups by
    ! repeatedly merging the two adjacent groups closest in quadratic diameter (relative
    ! gap). A wide gap is never bridged while a narrower one exists, so a bimodal stand
    ! keeps its two modes -- the information the Weibull re-draw destroyed. Empty cohorts
    ! merge first (zero distance) and are absorbed by a populated neighbour.
    IF (ok_maturity_transfer) THEN
       nact = 2*ncirc
       grp_n(1:nact) = temp_circ_class_n(1:nact)
       grp_nd2(1:nact) = temp_circ_class_n(1:nact) * (temp_circ_class(1:nact)/pi)**2
       DO WHILE (nact > ncirc)
          kmin = 1
          dmin = HUGE(un)
          DO icir = 1, nact-1
             IF (grp_n(icir) <= min_stomate .OR. grp_n(icir+1) <= min_stomate) THEN
                dist_qmd = zero
             ELSE
                qmd_a = SQRT(grp_nd2(icir) / grp_n(icir))
                qmd_b = SQRT(grp_nd2(icir+1) / grp_n(icir+1))
                ! Guillaume M. -- Sorted list + contiguous groups => qmd_b >= qmd_a.
                dist_qmd = (qmd_b - qmd_a) / MAX(qmd_a, min_stomate)
             ENDIF
             IF (dist_qmd < dmin) THEN
                dmin = dist_qmd
                kmin = icir
             ENDIF
          ENDDO
          grp_n(kmin) = grp_n(kmin) + grp_n(kmin+1)
          grp_nd2(kmin) = grp_nd2(kmin) + grp_nd2(kmin+1)
          DO icir = kmin+1, nact-1
             grp_n(icir) = grp_n(icir+1)
             grp_nd2(icir) = grp_nd2(icir+1)
          ENDDO
          nact = nact - 1
       ENDDO
       ! Guillaume M. -- The fill loop below consumes the sorted cohorts sequentially,
       ! so handing it the group sums reproduces the partition exactly (up to one
       ! rounding crumb at each group boundary, absorbed by the mass-balance rescale).
       exp_circ_class_n(1:ncirc) = grp_n(1:ncirc)
    ENDIF

    !! 2.4 Merge into a single distribution with ncirc classes
    !  Merge the 2*ncirc classes in ncirc classes. Empty one circ class
    !  after the other going from the smallest to the biggest. Here we
    !  still have 2*ncirc dimensions
    circ_class_biomass(ipts,new_pft,:,:,:) = zero
    trees_left(:) = temp_circ_class_n(:)

    DO icir=1,ncirc

       ! This is the number of trees in the ncirc circ classes thus
       ! AFTER the merge happened 
       trees_needed = exp_circ_class_n(icir)

       DO jcir = 1,2*ncirc

          IF (trees_needed .LE. zero) EXIT

          IF (trees_left(jcir) .GE. trees_needed) THEN

             ! We can get all the trees we need from this class
             ! and don't have to continue searching. Note that
             ! the units of circ_class_biomass have temporarly 
             ! changed to gC m-2
             circ_class_biomass(ipts,new_pft,icir,:,:) = &
                  circ_class_biomass(ipts,new_pft,icir,:,:) + &
                  temp_circ_class_biomass(jcir,:,:) * &
                  trees_needed
             trees_left(jcir) = trees_left(jcir) - trees_needed
             trees_needed = zero

             EXIT

          ELSE

             ! The number of trees in this class are not sufficient, so
             ! we need to move all of them to the new class. Note that
             ! the units of circ_class_biomass have temporarly changed 
             ! to gC m-2
             circ_class_biomass(ipts,new_pft,icir,:,:) = &
                  circ_class_biomass(ipts,new_pft,icir,:,:) + &
                  temp_circ_class_biomass(jcir,:,:) * &
                  trees_left(jcir)
             trees_needed = trees_needed - trees_left(jcir)
             trees_left(jcir) = zero

          ENDIF

       ENDDO ! jcir=1,ncirc   

    ENDDO ! icir=1,ncirc

    !! 2.5 Mass after the merge. Units are still in gC or gN m-2
    total_biomass_merge(:) = SUM(SUM(circ_class_biomass(ipts,new_pft,:,:,:),1),1)

    !! 2.6 Adjust to avoid small mass balance issue.
    DO iele = 1,nelements

       ! Adjust but be careful because adjusting could also hide real bugs.
       ! Checked this code several times and the mismatches were found to be 
       ! relatively small. Seems OK to adjust.
       mismatch = total_biomass_merge(iele) - total_biomass_begin(iele)
       IF (ABS(mismatch).GT.100*EPSILON(un)) THEN
          ! Only give a warning for ERR_ACT.GT.1 but always rescale
          IF(err_act.GT.1) THEN
             WRITE(numout,*) 'Imbalance, ', ipts, new_pft, iele, mismatch 
             CALL ipslerr_p(2,'merge_biomass_pfts','Small imbalances are thought to come &
                  from presicion errors', 'larger imblances i.e. 10e-7 were causing mass &
                  balance problems in sapiens_lcchange','Will adust Cs to solve the issue')
          END IF

          fixed = .FALSE.
          DO icir = ncirc,1,-1
             ! Adjust the merged biomass. Sapwood of the largest tree should be 
             ! able to easily deal with the expected mismatches (< 10e-7). Simply 
             ! subtract/add the mismatch to achieve the highest precision.
             IF(circ_class_biomass(ipts,new_pft,icir,isapbelow,iele)-mismatch.GT.zero)THEN
                circ_class_biomass(ipts,new_pft,icir,isapbelow,iele) = &
                     circ_class_biomass(ipts,new_pft,icir,isapbelow,iele) - mismatch
                fixed = .TRUE.
                EXIT
             END IF
          END DO

          ! Failed to fix the mass balance issue. This suggests that isapwood above
          ! is close to zero for all circumference classes. This is a bit strange 
          ! because it implies that there is almost no biomass to merge in the first 
          ! place.
          IF (.NOT.fixed.AND.err_act.GT.1) THEN
             WRITE(numout,*) 'circ_class_biomass after the merge, ', ipts, new_pft, iele, &
                  circ_class_biomass(ipts,new_pft,:,isapabove,iele)
             CALL ipslerr_p(3,'merge_biomass_pfts','failed to account for the mismatch', &
                  'This is unexpected - it suggest that isapabove is relly low','Check')
          END IF
       END IF ! mismatch .GT. 100*EPSILON(un)
    END DO ! iele

    !! 2.5 convert units back to the original units (gC tree-1) 
    !  Right now, circ_class_biomass gives the total biomass
    !  in each class, but it should only be for a model tree. 
    !  So let's normalize it. And, update the variables that 
    !  pass this information around
!!$    IF (SUM(exp_circ_class_n(:)).LT.min_stomate)THEN
!!$       ! This could be a suspicious case
!!$       ! or simply a coincidence. It means that the present PFT is empty
!!$       ! (just died) and that the loss_gain is less than min_stomate and
!!$       ! that for that reason we will not borther adding new biomass.    
!!$       circ_class_biomass(ipts,ivm,:,:,:) = zero
!!$          circ_class_n(ipts,ivm,:) = zero
!!$       !++++++++++++   
!!$    ELSE
       ! Guillaume M. -- /!\ GUARD RESTORED: it existed above but was COMMENTED OUT (!!$
       ! block). The division by exp_circ_class_n was unprotected, so an empty expected
       ! distribution gives 0/0 -> NaN in the receiver's circ_class_biomass. The case was
       ! IMPOSSIBLE before the PROGRESSIVE HARVEST, which moves BARE GROUND: with an empty
       ! receiver too, share_rec = 0 and both densities are zero.
       !
       ! Guillaume M. -- Mass-neutral, and bit-neutral when every density is > 0:
       ! total_biomass_begin is weighted by the two densities and vanishes when both are
       ! zero, and at class level exp_circ_class_n(icir) = 0 means loop 2.4 exited on
       ! `trees_needed .LE. zero`, leaving circ_class_biomass(icir) at its zero init.
       !
       ! Guillaume M. -- /!\ THAT REASONING HOLDS FOR EXACT ZERO ONLY. The test below is
       ! `.GT. min_stomate`: between 0 and min_stomate loop 2.4 did deposit mass, and zeroing
       ! it DESTROYED that mass. Handled below. The slot then becomes bare area without
       ! trees, which prescribe establishes afterwards.
       orphan_to_litter(:) = zero
       DO icir=1,ncirc

          IF (exp_circ_class_n(icir) .GT. min_stomate) THEN
             circ_class_biomass(ipts,new_pft,icir,:,:) = &
                  circ_class_biomass(ipts,new_pft,icir,:,:) / &
                  exp_circ_class_n(icir)
          ELSE
             ! Guillaume M. -- A class under the density threshold still carries the mass
             ! that loop 2.4 deposited into it, and zeroing it destroyed most of a slot.
             ! Send that carbon to litter instead. Biomass and bm_to_litter take the SAME
             ! area weight in the mass balance, so the transfer is exactly neutral.
             bm_to_litter(ipts,new_pft,:,:) = bm_to_litter(ipts,new_pft,:,:) + &
                  circ_class_biomass(ipts,new_pft,icir,:,:)
             orphan_to_litter(:) = orphan_to_litter(:) + &
                  SUM(circ_class_biomass(ipts,new_pft,icir,:,:),1)
             orphan_litter_hits = orphan_litter_hits + 1
             circ_class_biomass(ipts,new_pft,icir,:,:) = zero
          ENDIF
          circ_class_n(ipts,new_pft,icir) = exp_circ_class_n(icir)

       ENDDO
       orphan_litter_mass = orphan_litter_mass + orphan_to_litter(icarbon)
!!$    END IF
!!$    !+++++++++++

    !! 3. Mass balance check
    DO iele = 1,nelements

       ! If the mismatch is of 10e-10 the rescaling has little
       ! Guillaume M. -- orphan_to_litter closes the books: mass that left the biomass
       ! pool for litter is NOT missing. Without this term the test still fired on a
       ! slot that was already conserved, with numbers identical to the last digit --
       ! the deposit changes where the carbon goes, not how much of it there is.
       IF (share_rec * SUM(SUM(est_circ_class_biomass(:,:,iele),2) * &
            est_circ_class_n(:)) + (1-share_rec) * &
            SUM(SUM(tmp_circ_class_biomass(:,:,iele),2) * &
            tmp_circ_class_n(:)) - &
            SUM(SUM(circ_class_biomass(ipts,new_pft,:,:,iele),2) * &
            circ_class_n(ipts,new_pft,:)) - &
            orphan_to_litter(iele).GT.min_stomate) THEN
          WRITE(numout,*) 'Existing biomass (gC or gN m-2), ', &
               share_rec * SUM(SUM(est_circ_class_biomass(:,:,iele),2) * &
               est_circ_class_n(:))
          WRITE(numout,*) 'Newly established biomass (gC or gN m-2), ', &
               (1-share_rec) * SUM(tmp_circ_class_n(:) * &
               SUM(tmp_circ_class_biomass(:,:,iele),2))
          WRITE(numout,*) 'Biomass begin (gC or gN m-2), ', &
               share_rec * SUM(SUM(est_circ_class_biomass(:,:,iele),2) * &
               est_circ_class_n(:)) + (1-share_rec) * SUM(tmp_circ_class_n(:) * &
               SUM(tmp_circ_class_biomass(:,:,iele),2)), total_biomass_begin(iele)
          WRITE(numout,*) 'Merged biomass (gC or gN m-2), ', &
               SUM(SUM(circ_class_biomass(ipts,new_pft,:,:,iele),2) * &
               circ_class_n(ipts,new_pft,:)), total_biomass_merge(iele)
          WRITE(numout,*) 'Routed to litter (gC or gN m-2), ', orphan_to_litter(iele)
          WRITE(numout,*) 'Mismatch, ', share_rec * &
               SUM(SUM(est_circ_class_biomass(:,:,iele),2) * &
               est_circ_class_n(:)) + (1-share_rec) * &
               SUM(SUM(tmp_circ_class_biomass(:,:,iele),2) * &
               tmp_circ_class_n(:)) - &
               SUM(SUM(circ_class_biomass(ipts,new_pft,:,:,iele),2) * &
               circ_class_n(ipts,new_pft,:)) - &
               orphan_to_litter(iele)
          ! Guillaume M. -- numout is buffered: without this flush the whole diagnostic
          ! block above is lost when ipslerr_p aborts, and the run stops on a bare
          ! "mass balance problem" with no numbers.
          WRITE(numout,*) 'Context, ipts', ipts, ' new_pft', new_pft, &
               ' share_rec', share_rec, ' element', iele
          CALL flush(numout)
          CALL ipslerr_p(3,'merge_biomass_pfts','mass balance problem','','')

       END IF

    END DO

    IF (printlev.GE.4) WRITE(numout,*) 'Leaving merging biomass of pfts'
    
  END SUBROUTINE merge_biomass_pfts


! ================================================================================================================================
!! SUBROUTINE   : init_mass_balance
!!
!>\BRIEF        Initialize an extra but incomplete mass balance check
!!
!! DESCRIPTION  : It is difficult to check for mass balance closure in sapiens_lcchange
!!  because pools, fluxes and surface areas are all changing. This code checks whether 
!!  mass balance is presevered for the pixels that loose surface area. It also checks 
!!  mass conservation for the whole routine for litter, som, bact and soil_n_min because 
!! these pools don't receive or donate c/n to other pools.
!!
!! RECENT CHANGE(S) : 
!!
!! MAIN OUTPUT VARIABLE(S) : test_bank, test_mass 
!!
!! REFERENCES   : None
!!
!! FLOWCHART    : 
!! \n
!_ ================================================================================================================================

  SUBROUTINE init_mass_balance(test_bank, test_mass,&
       veget_max, litter, som, circ_class_biomass, circ_class_n, &
       turnover_daily, bm_to_litter, bact, soil_n_min, &
       harvest_pool, deepSOM_a, deepSOM_s, deepSOM_p, ipts, zf_soil)

  !! 0.1 Input variables
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(in)              :: circ_class_biomass    !! Biomass components of the model tree within a 
                                                                                        !! circumference class @tex $(g C ind^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(in)              :: harvest_pool          !! The wood and biomass that have been
                                                                                        !! havested by humans @tex $(gC)$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)                :: bm_to_litter          !! Transfer of biomass to litter
                                                                                        !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(in)              :: litter                !! metabolic and structural litter, above and 
                                                                                        !! below ground @tex ($gC m^{-2}$) @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(in)                    :: veget_max             !! "maximal" coverage fraction of a PFT on the ground
                                                                                        !! May sum to less than unity if the pixel has
                                                                                        !! nobio area. (unitless, 0-1)
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)                :: som                   !! carbon pool: active, slow, or passive 
                                                                                        !! @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)                :: deepSOM_a             !! Soil carbon discretized with depth active (g/m**3) 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)                :: deepSOM_s             !! Soil carbon discretized with depth slow (g/m**3) 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)                :: deepSOM_p             !! Soil carbon discretized with depth passive (g/m**3)
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)                :: turnover_daily        !! Transfer of litter to som
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)                  :: circ_class_n          !! Number of individuals in each circ class
                                                                                        !! @tex $(ind m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)                  :: soil_n_min            !! mineral nitrogen in the soil (gN/m**2)  
    REAL(r_std), DIMENSION(:,:), INTENT(in)                    :: bact                  !! denitrifier biomass (gC/m**2)
    INTEGER(i_std), INTENT(in)                                 :: ipts                  !! Pixel number (unitless)
    REAL(r_std), DIMENSION(0:ngrnd), INTENT(in)                :: zf_soil   
  !! 0.2 Output variables

  !! 0.3 Modified variables
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)               :: test_bank             !! Temporary variable to check mass conservation
    REAL(r_std), DIMENSION(:,:), INTENT(inout)                 :: test_mass             !! Temporary variable to check mass conservation 

  !! 0.4 Local variables
    INTEGER(i_std)                                             :: iele, ivm, ipar, icir !! Indices 
    INTEGER(i_std)                                             :: igrn                  !! Indices 
!_ ================================================================================================================================

    ! Initialize additional mass balance check
    test_bank(:,:,:) = zero
    test_mass(:,:) = zero
    
    ! Litter, som, bact and soil_n_min may get moved from 
    ! one PFT but the pools are never mixed with other pools. C/N moves
    ! between the remaining pools so those are checked as a single
    ! entity. If there is a mass balance problem, these pools can 
    ! be easily checked for mass conservation.
    DO iele = 1,nelements
       
       ! litter
       test_bank(1,iele,istart) = SUM(SUM(SUM(litter(ipts,:,:,:,iele),3),1) * &
            veget_max(ipts,:))
       ! som
       IF (ok_soil_carbon_discretization) THEN
          DO igrn = 1,ngrnd
             test_bank(2,iele,istart) = test_bank(2,iele,istart) + &
                  SUM((deepSOM_a(ipts,igrn,:,iele) + deepSOM_s(ipts,igrn,:,iele) + &
                  deepSOM_p(ipts,igrn,:,iele)) * (zf_soil(igrn)-zf_soil(igrn-1)) * &
                  veget_max(ipts,:))
          ENDDO
       ELSE
          test_bank(2,iele,istart) = SUM(SUM(som(ipts,:,:,iele),1) * &
               veget_max(ipts,:))
       END IF
       
    ENDDO
    
    ! Pure carbon or nitrogen pools: bact and soil_n_min
    test_bank(3,icarbon,istart) = SUM(bact(ipts,:) * veget_max(ipts,:))
    test_bank(3,initrogen,istart) = SUM(SUM(soil_n_min(ipts,:,:),2) * &
         veget_max(ipts,:))

  END SUBROUTINE init_mass_balance


! ================================================================================================================================
!! SUBROUTINE   : inter_mass_balance
!!
!>\BRIEF        intermediate but incomplete mass balance check
!!
!! DESCRIPTION  : It is difficult to check for mass balance closure in sapiens_lcchange
!!  because pools, fluxes and surface areas are all changing. This code checks whether 
!!  mass balance is presevered for the pixels that loose surface area. It also checks 
!!  mass conservation for the whole routine for litter, som, bact and soil_n_min because 
!! these pools don't receive or donate c/n to other pools.
!!
!! RECENT CHANGE(S) : 
!!
!! MAIN OUTPUT VARIABLE(S) : None
!!
!! REFERENCES   : None
!!
!! FLOWCHART    : 
!! \n
!_ ================================================================================================================================
  
  SUBROUTINE inter_mass_balance(test_bank, test_mass,&
       veget_max, litter, som, circ_class_biomass, circ_class_n, &
       turnover_daily, bm_to_litter, bact, soil_n_min, &
       harvest_pool, change_nobio, loss_gain, litter_bank, &
       soil_bank, fresh_ltr_bank, fresh_som_bank, total_losses, &
       test_mass_bank, bact_bank, soil_n_min_bank, &
       pool_start, deepSOM_a, deepSOM_s, deepSOM_p, som_vertres_bank, &
       ipts, ivm, dt_days, zf_soil)


!! 0.1 Input variables
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(in)              :: circ_class_biomass    !! Biomass components of the model tree within a 
                                                                                        !! circumference class @tex $(g C ind^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(in)              :: harvest_pool          !! The wood and biomass that have been
                                                                                        !! havested by humans @tex $(gC)$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)                :: bm_to_litter          !! Transfer of biomass to litter
                                                                                        !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(in)              :: litter                !! metabolic and structural litter, above and 
                                                                                        !! below ground @tex ($gC m^{-2}$) @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(in)                    :: veget_max             !! "maximal" coverage fraction of a PFT on the ground
                                                                                        !! May sum to less than unity if the pixel has
                                                                                        !! nobio area. (unitless, 0-1)
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)                :: som                   !! carbon pool: active, slow, or passive 
                                                                                        !! @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)                :: deepSOM_a             !! Soil carbon discretized with depth active (g/m**3) 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)                :: deepSOM_s             !! Soil carbon discretized with depth slow (g/m**3) 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)                :: deepSOM_p             !! Soil carbon discretized with depth passive (g/m**3)
    REAL(r_std), DIMENSION(:,:,:)                              :: som_vertres_bank      !! Bank to store the vertically-resolved SOM dilution (gC/m**2)
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)                :: turnover_daily        !! Transfer of litter to som
     REAL(r_std), DIMENSION(:,:,:), INTENT(in)                 :: circ_class_n          !! Number of individuals in each circ class
                                                                                        !! @tex $(ind m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)                  :: soil_n_min            !! mineral nitrogen in the soil (gN/m**2)
    REAL(r_std), DIMENSION(:,:), INTENT(in)                    :: bact                  !! denitrifier biomass (gC/m**2)
    REAL(r_std), DIMENSION(:), INTENT(in)                      :: change_nobio          !! Change in the non biological fraction in a pixel (0-1, unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(in)                    :: loss_gain             !! losses and gains due to LCC distributed over all
                                                                                        !! age classes and thus taking the age-classes into 
                                                                                        !! account (unitless, 0-1)
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)                  :: litter_bank           !! Bank to store the litter that becomes available when
                                                                                        !! part of a PFT is removed @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(in)                    :: soil_bank             !! Bank to store soil carbon that becomes available when
                                                                                        !! part of a PFT is removed  @tex ($gC m^{-2}$) @endtex
    REAL(r_std), INTENT(in)                                    :: bact_bank             !! Bank to store bacterial C-mass
    REAL(r_std), DIMENSION(:), INTENT(in)                      :: soil_n_min_bank       !! Bank to store mineral soil nitrogen
                                                                                        !! clearing during LCC. Weighted value of fresh_litter
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)                  :: pool_start            !! Start and end pool of this routine 
                                                                                        !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std),DIMENSION(:,:), INTENT(in)                     :: fresh_ltr_bank        !! Bank to store all the fresh litter from site
    REAL(r_std),DIMENSION(:,:), INTENT(in)                     :: fresh_som_bank        !! Bank to store all the fresh som from site 
    REAL(r_std)                                                :: total_losses          !! Sum of losses only in delta_veget (unitless, 0-1)
    INTEGER(i_std), INTENT(in)                                 :: ipts                  !! Pixel number (unitless)
    INTEGER(i_std), INTENT(in)                                 :: ivm                   !! PFT number
    REAL(r_std), INTENT(in)                                    :: dt_days               !! Time step of vegetation dynamics for stomate (days)
    REAL(r_std), DIMENSION(0:ngrnd), INTENT(in)                :: zf_soil

  !! 0.2 Output variables

  !! 0.3 Modified variables
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)               :: test_bank             !! Temporary variable to check mass conservation
    REAL(r_std), DIMENSION(:,:), INTENT(inout)                 :: test_mass             !! Temporary variable to check mass conservation 
    REAL(r_std), DIMENSION(:), INTENT(inout)                   :: test_mass_bank        !! Temporary variable to check mass conservation


  !! 0.4 Local variables
    INTEGER(i_std)                                             :: iele, par, icir !! Indices
    INTEGER(i_std)                                             :: igrn                  !! Indices  
!_ ================================================================================================================================

    ! We are only interested in the PFTs that lost part 
    ! of their share in the pixel, So we use a temporary
    ! variable for pool_start. To keep things simple
    ! the same is done for pool_end
    test_mass(:,istart) = test_mass(:,istart) + &
         pool_start(ipts,ivm,:)

    !! 6.1 Account for the carbon & nitrogen left on-site
    !  The amount of carbon/nitrogen per tree or per m-2 has
    !  not changed but the surface area covered by the PFT
    !  has decreased. Don't use the PFT dimension in pool_end
    !  because the banks no longer have PFT dimenions either
    DO iele = 1,nelements

       ! Biomass pool
       DO icir = 1,ncirc
          test_mass(iele,iend) = test_mass(iele,iend) + &
               SUM(circ_class_biomass(ipts,ivm,icir,:,iele))*&
               circ_class_n(ipts,ivm,icir) * &
               (veget_max(ipts,ivm)+loss_gain(ipts,ivm))
       ENDDO

       ! bm_to_litter and turnover_daily
       test_mass(iele,iend) = test_mass(iele,iend) + &
            SUM(bm_to_litter(ipts,ivm,:,iele) + &
            turnover_daily(ipts,ivm,:,iele)) * &
            (veget_max(ipts,ivm)+loss_gain(ipts,ivm))

       ! Litter pool (gC m-2) *  (m2 m-2) 
       test_mass(iele,iend) = test_mass(iele,iend) + &
            SUM(SUM(litter(ipts,:,ivm,:,iele),2)) * &
            (veget_max(ipts,ivm)+loss_gain(ipts,ivm))

       IF (ok_soil_carbon_discretization) THEN
          ! Soil carbon (gC m-2) *  (m2 m-2)
          DO igrn = 1,ngrnd
             test_mass(iele,iend) = test_mass(iele,iend) + &
                  (deepSOM_a(ipts,igrn,ivm,iele) + deepSOM_s(ipts,igrn,ivm,iele) + &
                  deepSOM_p(ipts,igrn,ivm,iele)) * (zf_soil(igrn)-zf_soil(igrn-1)) * &
                  (veget_max(ipts,ivm)+loss_gain(ipts,ivm))
          ENDDO
       ELSE
          ! Soil carbon (gC m-2) *  (m2 m-2)
          test_mass(iele,iend) = test_mass(iele,iend) + &
               SUM(som(ipts,:,ivm,iele)) * &
               (veget_max(ipts,ivm)+loss_gain(ipts,ivm))
       END IF

    ENDDO !nelements

    !Denitrifying bacteria (only C)
    test_mass(icarbon,iend) = test_mass(icarbon,iend) + &
         bact(ipts,ivm) * &
         (veget_max(ipts,ivm)+loss_gain(ipts,ivm))

    ! Nitrogen only
    test_mass(initrogen,iend) = test_mass(initrogen,iend) + &
         SUM(soil_n_min(ipts,ivm,:)) * & 
         (veget_max(ipts,ivm)+loss_gain(ipts,ivm))


    !! 6.2 Account for the pools where the C and N was moved to
    !  When the total surface area of a PFT was decreased, the 
    !  carbon and nitrogen that was removed was moved into several
    !  so called "banks" and the harvest pool. 
    DO iele = 1,nelements

       ! The biomass harvest pool is expressed in gC pixel-1 So, it 
       ! shouldn't be multiplied by veget_max but it should be 
       ! divided by area to obtain gC m-2.
       test_mass(iele,iend) = test_mass(iele,iend) + &
            SUM(SUM(harvest_pool(ipts,ivm,:,iele,:),2)) / &
            (area(ipts) * contfrac(ipts))
   
       test_mass_bank(iele) = test_mass_bank(iele) - &
            SUM(fresh_ltr_bank(:,iele) + fresh_som_bank(:,iele)) * &
            total_losses

       IF (ok_soil_carbon_discretization) THEN
          test_mass_bank(iele) = test_mass_bank(iele) - &
               SUM(SUM(som_vertres_bank(:,:,iele),2)) * total_losses
       ELSE
          test_mass_bank(iele) = test_mass_bank(iele) - &
               SUM(soil_bank(:,iele)) * total_losses
       ENDIF
       test_mass_bank(iele) = test_mass_bank(iele) - &
            SUM(SUM(litter_bank(:,:,iele),2)) * &
            total_losses

    ENDDO
 
    ! Pool with only carbon
    test_mass_bank(icarbon) = test_mass_bank(icarbon) - &
         bact_bank * total_losses

    ! Pools with only nitrogen
    test_mass_bank(initrogen) = test_mass_bank(initrogen) - &
         SUM(soil_n_min_bank(:)) * total_losses

    DO iele = 1,nelements
       IF (ABS(test_mass(iele,iend)-test_mass(iele,istart) + &
            test_mass_bank(iele)).GT.min_stomate) THEN
          WRITE(numout,*) 'ivm, iele, test_mass(istart), test_mass(iend), &
               & test_mass_bank, diff'
          WRITE(numout,'(2I3,4F20.10)') ivm, iele, &
               test_mass(iele,istart), test_mass(iele,iend), &
               test_mass_bank(iele), test_mass(iele,iend) - &
               test_mass(iele,istart) + test_mass_bank(iele)
          CALL ipslerr_p(3,'sapiens_lcchange','When accounting &
               & for the C and N due to decreases in PFTs', &
               'the C/N is not conserved','')
       ELSE
          ! Debug
          IF (printlev_loc>=4.AND.ipts==test_grid) THEN
             WRITE(numout,'(A,3I3)') 'OK intermediate mass balance &
                  &check for losses only, ',ipts,ivm,iele
          ENDIF
          !-
       ENDIF ! diff is zero

    ENDDO ! nelements

  END SUBROUTINE inter_mass_balance

  
! ================================================================================================================================
!! SUBROUTINE   : fin_mass_balance
!!
!>\BRIEF        finalize an extra but incomplete mass balance check
!!
!! DESCRIPTION  : It is difficult to check for mass balance closure in sapiens_lcchange
!!  because pools, fluxes and surface areas are all changing. This code checks whether 
!!  mass balance is presevered for the pixels that loose surface area. It also checks 
!!  mass conservation for the whole routine for litter, som, bact and soil_n_min because 
!! these pools don't receive or donate c/n to other pools.
!!
!! RECENT CHANGE(S) : 
!!
!! MAIN OUTPUT VARIABLE(S) : None
!!
!! REFERENCES   : None
!!
!! FLOWCHART    : 
!! \n
!_ ================================================================================================================================

  SUBROUTINE fin_mass_balance(test_bank,&
       veget_max, litter, som, circ_class_biomass, circ_class_n, &
       turnover_daily, bm_to_litter, bact, soil_n_min, &
       change_nobio, loss_gain, litter_bank, &
       soil_bank, bact_bank, soil_n_min_bank, &
       deepSOM_a, deepSOM_s, deepSOM_p, som_vertres_bank, &
       ipts, dt_days, zf_soil)


!! 0.1 Input variables
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(in)              :: circ_class_biomass    !! Biomass components of the model tree within a 
                                                                                        !! circumference class @tex $(g C ind^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)                :: bm_to_litter          !! Transfer of biomass to litter
                                                                                        !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(in)              :: litter                !! metabolic and structural litter, above and 
                                                                                        !! below ground @tex ($gC m^{-2}$) @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(in)                    :: veget_max             !! "maximal" coverage fraction of a PFT on the ground
                                                                                        !! May sum to less than unity if the pixel has
                                                                                        !! nobio area. (unitless, 0-1)
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)                :: som                   !! carbon pool: active, slow, or passive 
                                                                                        !! @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)                :: deepSOM_a             !! Soil carbon discretized with depth active (g/m**3) 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)                :: deepSOM_s             !! Soil carbon discretized with depth slow (g/m**3) 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)                :: deepSOM_p             !! Soil carbon discretized with depth passive (g/m**3)
    REAL(r_std), DIMENSION(:,:,:)                              :: som_vertres_bank      !! Bank to store the vertically-resolved SOM dilution (gC/m**2)
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)                :: turnover_daily        !! Transfer of litter to som
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)                  :: circ_class_n          !! Number of individuals in each circ class
                                                                                        !! @tex $(ind m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)                  :: soil_n_min            !! mineral nitrogen in the soil (gN/m**2)  
    REAL(r_std), DIMENSION(:,:), INTENT(in)                    :: bact                  !! denitrifier biomass (gC/m**2)
    REAL(r_std), DIMENSION(:)                                  :: change_nobio          !! Change in the non biological fraction in a pixel (0-1, unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(in)                    :: loss_gain             !! losses and gains due to LCC distributed over all
                                                                                        !! age classes and thus taking the age-classes into 
                                                                                        !! account (unitless, 0-1)
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)                  :: litter_bank           !! Bank to store the litter that becomes available when
                                                                                        !! part of a PFT is removed @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(in)                    :: soil_bank             !! Bank to store soil carbon that becomes available when
                                                                                        !! part of a PFT is removed  @tex ($gC m^{-2}$) @endtex
    REAL(r_std), INTENT(in)                                    :: bact_bank             !! Bank to store bacteria
    REAL(r_std), DIMENSION(:), INTENT(in)                      :: soil_n_min_bank       !! Bank to store mineral soil nitrogen
    INTEGER(i_std), INTENT(in)                                 :: ipts                  !! Pixel number (unitless)
    REAL(r_std), INTENT(in)                                    :: dt_days               !! Time step of vegetation dynamics for stomate (days)
    REAL(r_std), DIMENSION(0:ngrnd), INTENT(in)                :: zf_soil

  !! 0.2 Output variables

  !! 0.3 Modified variables
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)               :: test_bank              !! Temporary variable to check mass conservation

  !! 0.4 Local variables
    INTEGER(i_std)                                             :: iele, ivm, ipar, icir  !! Indices
    INTEGER(i_std)                                             :: igrn                  !! Indices  
!_ ================================================================================================================================
    
    DO iele = 1,nelements

       ! litter
       IF(change_nobio(ipts).GE.zero)THEN
          ! The litter bank has been distributed over the vegetated pixels
          ! and the new (increased) no_bio fraction.
          test_bank(1,iele,iend) = SUM(SUM(SUM(litter(ipts,:,:,:,iele),3),1) * &
               (veget_max(ipts,:)+loss_gain(ipts,:))) + &
               SUM(SUM(litter_bank(:,:,iele),1),1) * change_nobio(ipts)
       ELSE
          ! The litter bank was diluted by the virgin soil of the no bio
          ! that became vegetated. The entire litter bank was distributed over
          ! the vegetated pixels only.
          test_bank(1,iele,iend) = SUM(SUM(SUM(litter(ipts,:,:,:,iele),3),1) * &
               (veget_max(ipts,:)+loss_gain(ipts,:)))
       ENDIF

       ! som
       IF(change_nobio(ipts).GE.zero)THEN

          ! The som bank has been distributed over the vegetated pixels
          ! and the new (increased) no_bio fraction.
          IF (ok_soil_carbon_discretization) THEN
             DO igrn = 1,ngrnd     
                test_bank(2,iele,iend) = test_bank(2,iele,iend) + &
                     SUM((deepSOM_a(ipts,igrn,:,iele) + &
                     deepSOM_s(ipts,igrn,:,iele) + deepSOM_p(ipts,igrn,:,iele)) * &
                     (zf_soil(igrn)-zf_soil(igrn-1)) * &
                     (veget_max(ipts,:)+loss_gain(ipts,:))) + &
                     (SUM(som_vertres_bank(igrn,:,iele),1)) * change_nobio(ipts)
             ENDDO
          ELSE
             test_bank(2,iele,iend) = SUM(SUM(som(ipts,:,:,iele),1) * &
                  (veget_max(ipts,:)+loss_gain(ipts,:))) + &
                  SUM(soil_bank(:,iele),1) * change_nobio(ipts)
          END IF

       ELSE

          ! The som bank was diluted by the virgin soil of the no bio
          ! that became vegetated. The entire som bank was distributed over
          ! the vegetated pixels only.
          IF (ok_soil_carbon_discretization) THEN
             DO igrn = 1,ngrnd                  
                test_bank(2,iele,iend) =  test_bank(2,iele,iend) + &
                     SUM((deepSOM_a(ipts,igrn,:,iele) + &
                     deepSOM_s(ipts,igrn,:,iele) + deepSOM_p(ipts,igrn,:,iele)) * &
                     (zf_soil(igrn)-zf_soil(igrn-1)) * &
                     (veget_max(ipts,:)+loss_gain(ipts,:)))
             ENDDO
          ELSE
             test_bank(2,iele,iend) =  SUM(SUM(som(ipts,:,:,iele),1) * &
                  (veget_max(ipts,:)+loss_gain(ipts,:)))
          END IF
       ENDIF

    ENDDO

    IF(change_nobio(ipts).GE.zero)THEN
       ! The bact and soil_n_min bank has been distributed over the vegetated pixels
       ! and the new (increased) no_bio fraction.
       test_bank(3,icarbon,iend) = SUM(bact(ipts,:) * &
            (veget_max(ipts,:)+loss_gain(ipts,:))) + &
            bact_bank * change_nobio(ipts)                
       test_bank(3,initrogen,iend) = SUM(SUM(soil_n_min(ipts,:,:),2) * &
            (veget_max(ipts,:)+loss_gain(ipts,:))) + &
            SUM(soil_n_min_bank(:)) * change_nobio(ipts)

    ELSE
       ! The bact and soil_n_min bank was diluted by the virgin soil of the no bio
       ! that became vegetated. The entire bact and soil_n_min bank was distributed over
       ! the vegetated pixels only.
       test_bank(3,icarbon,iend) = SUM(bact(ipts,:) * &
            (veget_max(ipts,:)+loss_gain(ipts,:)))  
       test_bank(3,initrogen,iend) = SUM(SUM(soil_n_min(ipts,:,:),2) * &
            (veget_max(ipts,:)+loss_gain(ipts,:)))
    ENDIF

    ! Test mass conservation
    DO iele = 1,nelements

       IF (ABS(test_bank(1,iele,iend)-test_bank(1,iele,istart)).GT.min_stomate) THEN
          WRITE(numout,*) 'ipts, iele, test_bank(istart), test_bank(iend), diff'
          WRITE(numout,'(2I5,3F20.10)') ipts, iele,test_bank(1,iele,iend), test_bank(1,iele,istart),&
               test_bank(1,iele,iend)-test_bank(1,iele,istart)
          CALL ipslerr_p(plev,'sapiens_lcchange','mass balance not closed for litter','','')
       ENDIF
       IF (ABS(test_bank(2,iele,iend)-test_bank(2,iele,istart)).GT.min_stomate) THEN
          WRITE(numout,*) 'ipts, iele, test_bank(istart), test_bank(iend), diff'
          WRITE(numout,'(2I5,3F20.10)') ipts, iele,test_bank(2,iele,iend), test_bank(2,iele,istart),&
               test_bank(2,iele,iend)-test_bank(2,iele,istart)
          IF (ok_soil_carbon_discretization) THEN
             CALL ipslerr_p(plev,'sapiens_lcchange','mass balance not closed for deepSOM_x','','')
          ELSE
             CALL ipslerr_p(plev,'sapiens_lcchange','mass balance not closed for som','','')
          END IF
       ENDIF

       IF (ABS(test_bank(3,iele,iend)-test_bank(3,iele,istart)).GT.min_stomate) THEN
          WRITE(numout,*) 'ipts, iele, test_bank(istart), test_bank(iend), diff'
          WRITE(numout,'(2I5,3F20.10)') ipts, iele,test_bank(3,iele,iend), test_bank(3,iele,istart),&
               test_bank(3,iele,iend)-test_bank(3,iele,istart)
          IF (iele.EQ.1) CALL ipslerr_p(3,'sapiens_lcchange','Mass balance error for bact','','')
          IF (iele.EQ.2) CALL ipslerr_p(3,'sapiens_lcchange','Mass balance error for soil_n_min','','')
       ENDIF

       IF (ABS(test_bank(1,iele,iend)-test_bank(1,iele,istart)).GT.min_stomate) THEN
          WRITE(numout,*) 'ipts, iele, test_bank(istart), test_bank(iend), diff'
          WRITE(numout,'(2I5,3F20.10)') ipts, iele,test_bank(1,iele,iend), test_bank(1,iele,istart),&
               test_bank(1,iele,iend)-test_bank(1,iele,istart)
          CALL ipslerr_p(plev,'sapiens_lcchange','mass balance not closed for litter','','')
       ENDIF
    ENDDO ! iele

  END SUBROUTINE fin_mass_balance

!! ================================================================================================================================
!! SUBROUTINE   : check_read_vegetmax
!!
!>\BRIEF         To verify the consistency of the various fractions read from a land cover map
!!
!! DESCRIPTION  : This kind of error/precision checking should be
!!                done right after the maps are read and BEFORE any calculation
!!                is taking place. If fractions are changed later on in the code (even
!!                frac_nobio) mass balance problems are unavoidable unless changes
!!                in pools and fluxes are properly accounted for and checked which
!!                is the case in sapiens_lcchange.f90. 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): :: veget_max_new, frac_nobio_new
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

   SUBROUTINE check_read_vegetmax(kjpindex, veget_max_new, frac_nobio_new)

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                             :: kjpindex           !! Number of points for which the data needs to be interpolated

    !! 0.2 Output variables

    !! 0.3 Modified variables
    REAL(r_std),DIMENSION (:,:), INTENT(inout)             :: frac_nobio_new     !! Fraction of ice,lakes,cities, ... (unitless)
    REAL(r_std),DIMENSION (:,:), INTENT(inout)             :: veget_max_new      !! Maximum fraction of vegetation type including none 

    !! 0.4 Local variables
    INTEGER(i_std)                                         :: ipts,ivmax,ivm,jvm !! Index (unitless)
    REAL(r_std)                                            :: total              !! Total fraction returned from readvegetmax   
    INTEGER(r_std),DIMENSION(kjpindex)                     :: count              !! counter

!_ ================================================================================================================================     

    ! Start check. This kind of error/precision checking should be
    ! done right after the maps are read and BEFORE any calculation
    ! was done. If fractions are changed later on in the code (even
    ! frac_nobio) mass balance problems are unavoidable unless changes
    ! in pools and fluxes are properly accounted for and checked which
    ! is the case in sapiens_lcchange.f90.
    DO ipts=1,kjpindex

       ! The total sum of vegetation should never be greater than 1.0,
       ! since otherwise there will be a negative nobio fraction.  So
       ! check for that now, at the beginning, before doing all the other
       ! adjustments.  If there is a problem, scale all the values by a
       ! factor to set the sum exactly equal to one.
       total = SUM(veget_max_new(ipts,:)) + SUM(frac_nobio_new(ipts,:))
       IF(ABS(un-total) .GT. min_stomate)THEN

          IF(err_act.GT.1) THEN

             ! Error checking is strict and the single precision criterion
             ! has been violated. Stop the model.
             WRITE(numout,*) "ERROR: vegetation fractions exceed 1"
             WRITE(numout,*) "ipts, sum of veget_max_new: ",&
                  ipts, total
             CALL ipslerr_p (3,'check_read_vegetmax', 'The nobio fraction of &
                  this square on the new map is negative.', &
                  'Check the maps', '')

          ELSE

             ! The single precision criterion has been violated but
             ! error checking is less strict. Deal with the problem.
             ! Probably a good idea to carefully check the maps and/or
             ! calculations in slowproc.
             WRITE(numout,*) "WARNING: Scaling vegetation fractions to equal 1"
             WRITE(numout,*) "ipts, sum of veget_max_new: ",&
                  ipts, total
             WRITE(numout,*) "Write some code to monitor the problem that is &
                  hidden here and change the 3 into a 2 in ipslerr" 
             CALL ipslerr_p (3,'check_read_vegetmax', 'The nobio fraction of &
                  this square on the new map', 'is negative.  Scaling all &
                  vegetation fractions to compensate. Look for','Scaling  &
                  may cover up serious problems. Check the maps unless you &
                  are sure what you are doing.')

          END IF

       ELSE

          ! Sometimes maps are read-in which are created in single
          ! precision. Values in the code may be double precision,
          ! however min_stomate is usually set to 1e-8, which is
          ! similar to single precision. If we pass the min_stomate
          ! criterion above it is assumed that the small difference
          ! is due to a precision issue when redaing the maps. Rescale
          ! the fractions to avoid precision errors further down the code.
          ! Not accounting for this precision issue will cause mass 
          ! balance errors on the order of 1e-5 to 1e-8. The following 
          ! line of code basically converts single precision to double 
          ! precision. Given this is considered a precision issue there 
          ! is no need to notify the users. Note that the second dimension of 
          ! frac_nobio_new is nnobio.
          veget_max_new(ipts,:) = veget_max_new(ipts,:) / total
          frac_nobio_new(ipts,:) = frac_nobio_new(ipts,:) / total

       ENDIF

    END DO

    ! Quality check. It is now expected that the different vegetation 
    ! fractions in each pixel sums up to exactly one.
    CALL check_pixel_area("check_read_veget_max", kjpindex, veget_max_new, &
         frac_nobio_new)

    ! The code above is essential. If the surface area does not sum up to one
    ! the maps violate a basic assumption of ORCHIDEE. The cleaning coded below
    ! is not essential and was just introduced to avoid have many PFTs with very
    ! small cover fractions. Note that min_vegfrac is a relative threshold. It is
    ! set to 10e-6. The absolute value does differs depending on the surface area
    ! of a pixel. Although this cleaning does not relate to a basic assumption of
    ! ORCHIDEE, check_veget checks whether none of the fractions are less
    ! than min_vegfrac.
    ! At this point in the code veget_max_new is just read from the map and so
    ! it only contains values for the youngest age class of each PFT. It will be 
    ! checked whether none of the vegetation fraction read from the map are smaller
    ! than min_vegfrac. If the frac_nobio is smaller, the land surface will be put
    ! in the bare soil fraction. If the updated bare soil fraction is smaller than
    ! min_vegfrac it will be added to the PFT with the largest fraction. Given that
    ! all fraction should be positive and that only small fractions will be added
    ! to existing PFTs we don expect any issues with negative fractions.
    DO ipts = 1,kjpindex
     
       ! Check the non-biological fraction
       IF (SUM(frac_nobio_new(ipts,:)).LE.min_vegfrac) THEN
          
          ! Add the increase to the bare soil fraction. If the increase in
          ! frac_nobio_new is part of a trend, the treshold will be exceeded
          ! within a couple of years. In the meantime the land is left bare.
          ! This appears like a very reasonable assumption/approach. It may, 
          ! for example, take a couple of years between deforestation and 
          ! urbanisation. Note that in the next block of the code it will
          ! be tested whether ibare_sechiba was a suitable land use to move
          ! this residual fraction to.
          veget_max_new(ipts,ibare_sechiba) = veget_max_new(ipts,ibare_sechiba) + &
               SUM(frac_nobio_new(ipts,:))
          frac_nobio_new(ipts,:) = zero

       ENDIF

       ! Search for the largest cover fraction after land cover change to 
       ! add very small veget_max_new to. At this point in the code veget_max_new
       ! has just been read from the map and so it only contains values for the
       ! youngest age class of each PFT. Because there are no age classes
       ! yet, there will always be a solution (but the solution could be in
       ! frac_nobio_new. Note that dim=1 refers to the dimensions of the answer 
       ivmax = MAXLOC(veget_max_new(ipts,:),DIM=1)

       ! Check whether ivmax is a suitable PFT itself
       IF (veget_max_new(ipts,ivmax).LT.min_vegfrac) THEN

          ! Even the PFT with the largest vegetation fraction
          ! does not exceed the threshold. Move everything into
          ! the non-biological fraction.
          !+++CHECK+++
          ! for the moment only one non-biological land use is
          ! defined so add it in there. When more land uses will
          ! be defined, consider which land use is the most suitable
          frac_nobio_new(ipts,1) = frac_nobio_new(ipts,1) + &
               SUM(veget_max_new(ipts,:))
          veget_max_new(ipts,:) = zero
          !+++++++++++
          
       ELSE

          ! The PFT with the largest vegetation fraction is suitable
          ! to accept the residual fractions of other PFTs
          DO ivm = 1,nvm
                
             ! Check whether the fraction is large enough to be considered
             ! a liveable PFT
             IF (veget_max_new(ipts,ivm).LE.min_vegfrac) THEN

                ! Add to the cover fraction of the PFT with the largest fraction
                ! after land cover change 
                veget_max_new(ipts,ivmax) = veget_max_new(ipts,ivmax) + &
                     veget_max_new(ipts,ivm)
                veget_max_new(ipts,ivm) = zero                
             
             END IF ! veget_max_new(ivm) .LT. min_vegfrac

          END DO ! ivm

       END IF ! veget_max_new(ivmax).LT.min_vegfrac

    END DO ! ipts
    
  END SUBROUTINE check_read_vegetmax


! ================================================================================================================================
!! SUBROUTINE   : adjust_delta_veget_max
!!
!>\BRIEF        
!!
!! DESCRIPTION  :  Check whether the initially proposed changes in veget_max make sense
!! If they do so, accept this change. If the change brings the remaining veget_max below
!! the minimal threshold remove the whole PFT. Further the maps only prescribe the change
!! in nvmap and do not contain any information about the changes in age classes. Convert the
!! changes from nvmap (= age groups) to nvm (= PFT level). If a PFT looses surface area we
!! assumed that each age classes looses the same relative amount of its surface area.
!!
!! RECENT CHANGE(S) : 
!!
!! MAIN OUTPUT VARIABLE(S) :  loss_gain, veget_max_new
!!
!! REFERENCES   : None
!!
!! FLOWCHART    : 
!! \n
!_ ================================================================================================================================

  SUBROUTINE adjust_delta_veget_max(kjpindex, veget_max, frac_nobio, &
       veget_max_new, frac_nobio_new, loss_gain, losses, failed_vegfrac)

    
  !! 0. Variable and parameter declaration 
    
    !! 0.1 Input variables
    INTEGER, INTENT(in)                                  :: kjpindex                !! Domain size - number of pixels (unitless) 
    REAL(r_std), DIMENSION(:,:), INTENT(in)              :: veget_max               !! Cover fraction of a PFT (unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(in)              :: frac_nobio              !! Fraction of the mesh which is covered by ice, lakes, ...
    CHARACTER(30), INTENT(in)                            :: losses                  !! Flag that determines how the losses are 
                                                                                    !! distributed. There are currently two options
                                                                                    !! 'proportional' and 'youngest'
    
    !! 0.2 Output variables

    !! 0.3 Modified variables
    REAL(r_std), DIMENSION(:,:), INTENT(inout)           :: veget_max_new           !! Coverage fraction of a PFT (unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)           :: frac_nobio_new          !! Fraction of the mesh which is covered by ice, lakes, ...
    REAL(r_std), DIMENSION(:,:), INTENT(inout)           :: loss_gain               !! The same as delta_veget_cov but distributed of all
                                                                                    !! age classes and thus taking the age-classes into 
    LOGICAL, DIMENSION(:,:), INTENT(inout)               :: failed_vegfrac          !! Failed to find a PFT were some residual fraction could be added (true/false)
    
    !! 0.4 Local variables
    INTEGER(i_std)                                       :: ivm,ivma,ipts,ivmax     !! Indices
    INTEGER(i_std)                                       :: igroup,iyoung,innob     !! Indices
    INTEGER(i_std)                                       :: nagec_local             !! The number of age classes in this map PFT.
    REAL(r_std),DIMENSION(kjpindex,nvmap)                :: delta_veget_cov         !! changes in "maximal" coverage fraction of PFT 
    REAL(r_std), DIMENSION(kjpindex,nvmap)               :: veget_max_group         !! Total veget_max of all age clases within a PFT group 
                                                                                    !! (unitless, 0-1)
    REAL(r_std), DIMENSION(kjpindex,nvmap)               :: veget_max_new_group     !! Total veget_max_new of all age clases within a PFT group 
                                                                                    !! (unitless, 0-1)
    INTEGER(i_std), DIMENSION(kjpindex)                  :: count                   !! counter
    REAL(r_std)                                          :: total                   !! temporary variable for total area of a pixel (unitless, 0-1)
    REAL(r_std), DIMENSION(kjpindex)                     :: change_nobio            !! change in the fractions of non-biological covers (unitless, 0-1)
    
!_ ================================================================================================================================

    IF(printlev_loc.GE.2) WRITE(numout,*) 'Entering adjust_delta_veget_max'

    ! Initialize
    veget_max_group(:,:) = zero
    veget_max_new_group(:,:) = zero
    delta_veget_cov(:,:) = zero
   
    !  Remember that when age classes are being used the maps have
    !  fewer PFTs (::nvmap) than the simulation (::nvm). For several 
    !  calculations the total cover of all age classes within its group 
    !  (::veget_max_group) is required. We also need the vegetation 
    !  fraction of the PFT map arranged into groups. Notice that this 
    !  number is already given in terms of groups, but it is stored 
    !  on a grid using all nvm pfts, i.e. only nvmap values of
    !  veget_max_new(ipts,:) are non-zero.
    DO ipts = 1,kjpindex
       DO ivm=1,nvm
          igroup=agec_group(ivm)
          veget_max_group(ipts,igroup) = veget_max_group(ipts,igroup) + &
               veget_max(ipts,ivm)
          veget_max_new_group(ipts,igroup) = veget_max_new_group(ipts,igroup) + &
               veget_max_new(ipts,ivm)
       ENDDO
    ENDDO

    ! Debug
    IF (printlev_loc>=4) THEN
       DO ipts=1,kjpindex
          IF(ipts == test_grid)THEN
             WRITE(numout,*) 'Finished computing veget_max_group, ', ipts
             CALL print_lcc_error_message(kjpindex,ipts,veget_max_new,&
                  veget_max,veget_max_new_group, failed_vegfrac,&
                  veget_max_group,delta_veget_cov,loss_gain,1)
             CALL print_lcc_error_message(kjpindex,ipts,veget_max_new,&
                  veget_max,veget_max_new_group, failed_vegfrac,&
                  veget_max_group,delta_veget_cov,loss_gain,2)
             CALL print_lcc_error_message(kjpindex,ipts,veget_max_new,&
                  veget_max,veget_max_new_group, failed_vegfrac,&
                  veget_max_group,delta_veget_cov,loss_gain,5)
          ENDIF
       ENDDO
    ENDIF
    !-    
   
    ! Rather than having one big loop in which the vegetation fractions
    ! are tested PFT by PFT, a sequence of tests was set-up. The order
    ! of the tests matters. Start with calculating the change in vegetation cover
    delta_veget_cov(:,:) =  veget_max_new_group(:,:) - veget_max_group(:,:)

    ! Debug  
    IF (printlev_loc>=4) THEN
       DO ipts=1,kjpindex
          IF(ipts == test_grid)THEN
             WRITE(numout,*) 'Change in vegetation cover - before QC, ', ipts, &
                  SUM(delta_veget_cov(test_grid,:))
             CALL print_lcc_error_message(kjpindex,test_grid,veget_max_new,&
                  veget_max,veget_max_new_group, failed_vegfrac,&
                  veget_max_group,delta_veget_cov,loss_gain,3)
          ENDIF
       ENDDO
    ENDIF
    !-

    ! ::delta_veget_cov has been calculated for the sum of the age classes
    ! within a MTC. ::delta_veget_cov thus needs to be distributed over 
    ! all the age classes. This distribution requires a rule which 
    ! in the absence of better information is based on the assumption
    ! that a land cover change does not care about the age classes.
    ! Following a forestry logic, one would change an old forest into 
    ! a grassland because that would result in an income from the wood
    ! harvest. However, here we assume that LCC follows a spatialized 
    ! economic need and therefore does not respect the forestry logic.
    ! In other words, when a farmer needs more land, (s)he will not cut
    ! spatialy distributed old stands but will target a block of land
    ! bordering his/her current agricultural lands. By doing so forest
    ! of all age classes will be cut when LCC is the underlying 
    ! motivation. The reasoning applies for deforestation driven by 
    ! urbanisation. This assumption translates in that the land cover
    ! fraction of all age classes will decrease proportionally to their
    ! share when the PFT experiences a loss. When the PFT experiences an
    ! increase in the PFT cover this only happens in the youngest age
    ! class. There is one more constraint.  The new vegetation fraction
    ! at the end of the calculation must be larger than min_vegfrac. If
    ! veget_max is greater than zero but less than min_vegfrac after
    ! accounting for loss_gain, we kill the whole PFT. Note that frac_nobio
    ! is treated as the residual in these calculations.
    DO ipts = 1,kjpindex

       IF(printlev_loc.GE.4)THEN
          IF(test_grid == ipts)THEN
             WRITE(numout,*) 'Readjusting LCC fluxes for grid: ',ipts
          ENDIF
       ENDIF

       ! From here on, we are only concerned with loss_gain.  We will ignore
       ! delta_veget_cov after creating loss_gain. Note that loss_gain has 
       ! the dimensions of the simulation (::nvm) whereas delta_veget_cov has 
       ! the dimension of the maps (::nvma)
       DO ivm=1,nvm

          igroup=agec_group(ivm)
          iyoung=start_index(igroup)

          ! Only forest PFTs have age classes
          IF (is_tree(ivm)) THEN

             ! PFT cover area decreases (loss)
             IF (delta_veget_cov(ipts,igroup) .LT. zero) THEN

                SELECT CASE(losses)
                CASE('proportional')

                   ! When there is a loss in veget_max the standard code
                   ! will distribute this loss proportionaly to the populated
                   ! age classes of the species group. This reflects the
                   ! fact that a land cover change does not necessarily 
                   ! happens at the end of a forest rotation. All, young, 
                   ! medium and old forest will be deforested.
                   loss_gain(ipts,ivm) = veget_max(ipts,ivm) / &
                        veget_max_group(ipts,igroup) * &
                        delta_veget_cov(ipts,igroup)

                CASE('youngest')

                   ! This option is intended for species change only.
                   ! A species change only happens at the end of a 
                   ! rotation. So it is an old age class that will be
                   ! replaced by a young one. Note that in forestry the
                   ! old age class is marked for killing. In kill the
                   ! surface of the harvested age class is moved to
                   ! the youngest age class and here we move it 
                   ! from the youngest age class of the previous species
                   ! to the youngest age class of the new species.
                   !loss_gain(ipts,ivm) = delta_veget_cov(ipts,igroup)
                   !delta_veget_cov(ipts,igroup) = zero

                   ! Previous lines that calculate loss_gain by group level only for
                   ! the youngest age caused an error when older PFT is
                   ! harvested and species changes at the same time (ticket #667).
                   ! Because clearing out happens after age_class_distribution, the
                   ! veget_max for old PFT is not moved to the youngest PFT, and
                   ! veget_max_new is 0 for both PFTs due to species change.
                   ! This led mismatch in loss_gain for the youngest PFT.  To
                   ! prevent the error,  for the case 'youngest', loss_gain is
                   ! calculated based on PFT level (ivm), not on the age group
                   ! (igroup). 
                   loss_gain(ipts,ivm) = veget_max_new(ipts,ivm) - veget_max(ipts,ivm)

                CASE DEFAULT

                   ! Stop the simulation
                   CALL ipslerr_p (3,'adjust_delta_vegetmax', &
                        'No information on',&
                        'how the losses in veget_max should be distributed','')

                END SELECT

             ! PFT cover area increases (gain)
             ELSEIF (delta_veget_cov(ipts,igroup) .GT. zero) THEN

                ! The first age-class in a PFT group is the youngest
                ! and that is the age-class were the change in cover
                ! fraction takes place
                IF(ivm == iyoung)THEN

                   loss_gain(ipts,ivm) = delta_veget_cov(ipts,igroup)
                ELSE

                   ! The remaining age-classes in that PFT don't experience a
                   ! change in cover fraction  
                   loss_gain(ipts,ivm) = zero

                ENDIF ! looking for youngest age class

             ELSEIF (delta_veget_cov(ipts,igroup) .EQ. zero) THEN

                ! There are no changes in the age group
                loss_gain(ipts,ivm) = zero
                
             ELSE

                ! Something went wrong
                WRITE(numout,*) 'ipts, igroup, delta_veget_cov, ', &
                     ipts, igroup, delta_veget_cov(ipts,igroup)
                CALL ipslerr_p (3,'adjust_veget_max','unexpected value',&
                     'for delta_veget_cov','')
             ENDIF

          ! Grasses and croplands
          ELSE

             ! There are no age-classes for grasses and crops, hence, 
             ! delta_veget_cov should not be distributed over the age classes
             loss_gain(ipts,ivm) = delta_veget_cov(ipts,igroup)

          ENDIF ! is_tree(ivm)

       ENDDO ! looping over PFTs
       

       IF (printlev_loc.GE.4) THEN
          IF(ipts == test_grid)THEN
             WRITE(numout,*) 'First estimate of loss_gain '
             CALL print_lcc_error_message(kjpindex,ipts,veget_max_new,&
                  veget_max,veget_max_new_group, failed_vegfrac,&
                  veget_max_group,delta_veget_cov,loss_gain,4)
          ENDIF
       ENDIF

       ! frac_nobio is considered to be the most static component of the land
       ! surface. After all glaciers, lakes and cities don change too fast
       ! and if they do change they should change along a one-directional
       ! trend. For the moment the code distinguishes just a single nnobio. This code
       ! is ready for more non-biological land covers such as glacioers, lakes,
       ! and cities
       DO innob = 1,nnobio

          ! Similar to veget_max, following a land cover change the remaining
          ! cover fraction of each non-biological land cover type should exceed
          ! the threshold (::min_vegfrac)
          IF (frac_nobio_new(ipts,innob).LT.min_vegfrac) THEN 

             ! Move it to the bare soil even if this implies that the
             ! bare soil does not yet exceed the threshold of min_vegfrac
             loss_gain(ipts,ibare_sechiba) = loss_gain(ipts,ibare_sechiba) + &
                  frac_nobio_new(ipts,innob)
             
             ! Set frac_nobio_new to its new value (which is zero)
             frac_nobio_new(ipts,innob) = zero

          END IF 
          
       END DO ! innob
       
       ! After a land cover change we don't want to have fractions below
       ! min_vegfrac. Keeping small fractions would inflate the number of PFTs
       ! and thus increase the computational cost for very little added 
       ! information/results.
       DO ivm = 1,nvm

          IF (veget_max(ipts,ivm)+loss_gain(ipts,ivm).LT.-min_stomate) THEN

             ! Unexpected large mismatch. This implies that the loss exceeds
             ! the initial veget_max. Given the way loss_gain was calculated
             ! this is a very strange outcome. Something went wrong
             WRITE(numout,*) 'ipts, ivm, veget_max, loss_gain, mismatch, ', &
                  ipts, ivm, veget_max(ipts,ivm), loss_gain(ipts,ivm), &
                  veget_max(ipts,ivm)+loss_gain(ipts,ivm)
             CALL ipslerr_p (3,'adjust_veget_max','conflicting value for',&
                  'loss_gain in combination with veget_max','')

          ELSEIF ((veget_max(ipts,ivm)+loss_gain(ipts,ivm).GE.-min_stomate .AND. &
               veget_max(ipts,ivm)+loss_gain(ipts,ivm).LT.zero) .OR. &
               (veget_max(ipts,ivm)+loss_gain(ipts,ivm).GT.zero .AND. &
               veget_max(ipts,ivm)+loss_gain(ipts,ivm).LE.min_stomate)) THEN   

             ! Exclude zero from this IF-statement else all the empty PFTs will
             ! be treated by this IF-statement.
             ! Small mismatch. If it is of the order of 10e-8 it is unexpected
             ! if it is 10e-15 or smaller it is a precision issues. We'll treat
             ! both issues as a precision error for now but a warning is issued
             ! so that the users can judge for themself.
             WRITE(numout,*) 'ipts, ivm, veget_max, loss_gain, mismatch, ', &
                  ipts, ivm, veget_max(ipts,ivm), loss_gain(ipts,ivm), &
                  veget_max(ipts,ivm)+loss_gain(ipts,ivm)
             CALL ipslerr_p (2,'adjust_veget_max','small mismatch for',&
                  'loss_gain in combination with veget_max',&
                  '10e-8 is unexpected but 10e-15 is an acceptable precision issue')

             ! Try to clean these precision issues from loss_gain
             CALL check_loss_gain(ipts, ivm, veget_max, loss_gain, &
                  failed_vegfrac)

          ELSEIF (veget_max(ipts,ivm)+loss_gain(ipts,ivm).GT.zero .AND. &
               veget_max(ipts,ivm)+loss_gain(ipts,ivm).LE.min_vegfrac) THEN

             ! Try to clean too small vegetation fractions from loss_gain
             CALL check_loss_gain(ipts, ivm, veget_max, loss_gain, &
                  failed_vegfrac)
                
          ELSEIF (veget_max(ipts,ivm)+loss_gain(ipts,ivm)-un.GT.10*EPSILON(un)) THEN
             
             ! Unexpected large mismatch. This implies that new vegetation
             ! fraction will exceed 1. Given the way loss_gain was calculated
             ! something must have gone wrong.
             WRITE(numout,*) 'ipts, ivm, veget_max, loss_gain, mismatch, ', &
                  ipts, ivm, veget_max(ipts,ivm), loss_gain(ipts,ivm), &
                  veget_max(ipts,ivm)+loss_gain(ipts,ivm)
             CALL ipslerr_p (3,'adjust_veget_max','conflicting value for',&
                  'loss_gain in combination with veget_max',&
                  'new PFT fractions exceeds the entire pixel')
             
          ENDIF

       END DO ! nvm

       ! Error checking with delta_veget_cov, if delta_veget_cov is well
       ! distributed to loss_gain over PFTs.
       DO ivm=1,nvm

          igroup=agec_group(ivm)
          iyoung=start_index(igroup)
          delta_veget_cov(ipts,igroup) = delta_veget_cov(ipts,igroup) -  &
                                       loss_gain(ipts,ivm)
                                  
          IF (ivm .EQ. iyoung + nagec_pft(igroup) -1 ) THEN
              IF (delta_veget_cov(ipts,igroup) .GT. min_stomate .OR. &
                  delta_veget_cov(ipts,igroup) .LT. -min_stomate) THEN
                  CALL ipslerr_p (2,'adjust_veget_max','Not all of delta_veget_cov was',&
                      'distributed to loss_gain','')
              ENDIF
          ENDIF
       ENDDO

       ! Debug
       IF (printlev_loc.GE.4) THEN
          IF(ipts == test_grid)THEN
             WRITE(numout,*) 'Readjusting loss_gain to make sure an &
                  age class does not fall too low.'
             CALL print_lcc_error_message(kjpindex,ipts,veget_max_new,&
                  veget_max,veget_max_new_group, failed_vegfrac,&
                  veget_max_group,delta_veget_cov,loss_gain,4)
             
          END IF
       END IF
       !-
                          
    END DO ! kjpindex

    ! Recalculate veget_max_new based on realistic values for loss_gain. Note
    ! that the frac_nobio_new was updated as well. 
    veget_max_new(:,:) = veget_max(:,:) + loss_gain(:,:)

    ! Debug
    IF (printlev_loc.GE.4) THEN
       WRITE(numout,*) 'Final values for veget_max_new'
       CALL print_lcc_error_message(kjpindex,test_grid,veget_max_new,&
            veget_max,veget_max_new_group, failed_vegfrac,&
            veget_max_group,delta_veget_cov,loss_gain,4)
    END IF
    !- 
    
    ! Quality check. It is now expected that the changes cancel out each other
    ! and are exactley zero.
    CALL check_area_change("adjust_delta_vegetmax", kjpindex, frac_nobio, &
         frac_nobio_new, loss_gain)

    ! Quality check. It is still expected that the different vegetation 
    ! fractions in each pixel sums up to exactly one.
    CALL check_pixel_area("adjust_delta_vegetmax", kjpindex, &
         veget_max_new, frac_nobio_new)

    IF(printlev_loc.GT.3) WRITE(numout,*) 'Leaving adjust_delta_veget_max'
    
  END SUBROUTINE adjust_delta_veget_max


  !! ================================================================================================================================
!! SUBROUTINE   : check_veget
!!
!>\BRIEF         To verify the consistency of the various fractions defined within the grid box after having been
!!               been updated by STOMATE or the standard procedures.
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): :: none
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!
  SUBROUTINE check_veget(kjpindex, frac_nobio, veget_max, veget, tot_bare_soil, &
       soiltile, failed_vegfrac)
 
    !!  0.1 Input variables
    INTEGER(i_std), INTENT(in)                          :: kjpindex       !! Number of points for which the data needs to be interpolated
    REAL(r_std),DIMENSION (:,:), INTENT(in)             :: frac_nobio     !! Fraction of ice,lakes,cities, ... (unitless)
    REAL(r_std),DIMENSION (:,:), INTENT(in)             :: veget_max      !! Maximum fraction of vegetation type including none 
                                                                          !! biological fraction (unitless) 
    REAL(r_std),DIMENSION (:,:), INTENT(in)             :: veget          !! Vegetation fractions in the mesh (unitless)
    REAL(r_std),DIMENSION (:), INTENT(in)               :: tot_bare_soil  !! Total evaporating bare soil fraction (unitless)
    REAL(r_std),DIMENSION (:,:), INTENT(in)             :: soiltile       !! Fraction of soil tile within vegtot (0-1, unitless)
    LOGICAL, DIMENSION(:,:), INTENT(in)                 :: failed_vegfrac !! Failed to find a PFT were some residual fraction
                                                                          !! could be added (true/false)

    !!  0.2 Output variables

    !!  0.3 Modified variables
    
    !!  0.4 Local variables
    INTEGER(i_std) :: ji, jn, jv
    REAL(r_std)  :: epsilocal  !! A very small value
    REAL(r_std)  :: totfrac
    CHARACTER(len=80) :: str1, str2
 
!_ ================================================================================================================================
    
    !
    ! There is some margin added as the computing errors might bring us above EPSILON(un)
    !
    epsilocal = EPSILON(un)*100
    
    !! 1.0 Verify that none of the fractions are smaller than min_vegfrac, without beeing zero.
    !!
    DO ji=1,kjpindex
       DO jn=1,nnobio
          IF ( frac_nobio(ji,jn) > epsilocal .AND. frac_nobio(ji,jn) < min_vegfrac ) THEN
             WRITE(str1,'("Occurs on grid box", I8," and nobio type ",I3 )') ji, jn
             WRITE(str2,'("The small value obtained is ", E14.4)') frac_nobio(ji,jn)
             CALL ipslerr_p (3,'check_veget', &
                  "frac_nobio is larger than zero but smaller than min_vegfrac.", str1, str2)
          ENDIF
       ENDDO
    END DO
    
    IF (.NOT. ok_dgvm) THEN       
       DO ji=1,kjpindex
          DO jv=1,nvm
             ! If we failed to find a suitable PFT to add the residual fraction to
             ! a PFT with a large cover fraction we will have at least one PFT in this
             ! pixel with a fraction less than min_vegfrac. If we are aware about this
             ! issue it should be registered in failed_vegfrac.
             IF ( veget_max(ji,jv) > epsilocal .AND. &
                  veget_max(ji,jv) < min_vegfrac .AND. &
                  .NOT. failed_vegfrac(ji,jv) ) THEN
                WRITE(str1,'("Occurs on grid box", I8," and nvm ",I3 )') ji, jv
                WRITE(str2,'("The small value obtained is ", E14.4)') veget_max(ji,jv)
                CALL ipslerr_p (3,'check_veget', &
                     "veget_max is larger than zero but smaller than min_vegfrac.", str1, str2)
             ENDIF
          ENDDO
       ENDDO
    END IF
    
    !! 2.0 verify that with all the fractions we cover the entire grid box  
    !!
    DO ji=1,kjpindex
       totfrac = zero
       DO jn=1,nnobio
          totfrac = totfrac + frac_nobio(ji,jn)
       ENDDO
       DO jv=1,nvm
          totfrac = totfrac + veget_max(ji,jv)
       ENDDO
       IF ( ABS(totfrac - un) > epsilocal) THEN
             WRITE(str1,'("This occurs on grid box", I8)') ji
             WRITE(str2,'("The sum over all fraction and error are ", E14.4, E14.4)') totfrac, ABS(totfrac - un)
             CALL ipslerr_p (3,'check_veget', &
                   "veget_max + frac_nobio is not equal to 1.", str1, str2)
             WRITE(*,*) "EPSILON =", epsilocal 
       ENDIF
    ENDDO
    
    !! 3.0 Verify that veget is smaller or equal to veget_max
     DO ji=1,kjpindex
       DO jv=1,nvm
          IF ( jv == ibare_sechiba ) THEN
             IF ( ABS(veget(ji,jv)) > epsilocal ) THEN
                WRITE(str1,'("This occurs on grid box", I8)') ji
                WRITE(str2,'("The difference is ", E14.4)') veget(ji,jv) - veget_max(ji,jv)
                CALL ipslerr_p (3,'check_veget', &
                     "veget is not equal to zero on bare soil.", str1, str2)
             ENDIF
          ELSE
             IF ( veget(ji,jv) > veget_max(ji,jv) ) THEN
                WRITE(str1,'("This occurs on grid box", I8)') ji
                WRITE(str1,'("This occurs on PFT", I8)') jv
                WRITE(str2,'("The values for veget and veget_max :", F8.4, F8.4)') veget(ji,jv), veget_max(ji,jv)
                CALL ipslerr_p (3,'check_veget', &
                     "veget is greater than veget_max.", str1, str2)
             ENDIF
          ENDIF
       ENDDO
    ENDDO
 
   
    !! 4.0 Test tot_bare_soil in relation to the other variables
    DO ji=1,kjpindex
       IF (ok_bare_soil_new) THEN
          totfrac = veget_max(ji,ibare_sechiba)
       ELSE
          totfrac = zero
          DO jv=1,nvm
             totfrac = totfrac + (veget_max(ji,jv) - veget(ji,jv))
          ENDDO
          ! add the bare soil fraction to totfrac
          totfrac = totfrac + veget(ji,ibare_sechiba)
       ENDIF

       ! do the test
       IF ( ABS(totfrac - tot_bare_soil(ji)) > epsilocal ) THEN
          WRITE(str1,'("This occurs on grid box", I8)') ji
          WRITE(str2,'("The values for tot_bare_soil, tot frac and error :", F8.4, F8.4, E14.4)') &
               &  tot_bare_soil(ji), totfrac, ABS(totfrac - tot_bare_soil(ji))
          CALL ipslerr_p (3,'check_veget', &
               "tot_bare_soil does not correspond to the total bare soil fraction.", str1, str2)
       ENDIF
    ENDDO
    
    !! 5.0 Test that soiltile has the right sum
    DO ji=1,kjpindex
       totfrac = SUM(soiltile(ji,:))
       IF (ABS(totfrac - un) > epsilocal .AND. totfrac .NE. zero) THEN
          ! The second condition could be caused by terrestrial pixels without frac_bio such as large parts 
          ! of Greenland and some of the artic islands north of the Russian mainland. This is intended 
          ! model behavior and therefore not a problem. For all pixels with a bare soil and/or vegetation 
          ! it is expected that totfrac adds up to 1 because soiltile is normalized by the frac_bio earlier 
          ! in this module. If this not the case, there is a problem.
          WRITE(numout,*) "soiltile does not sum-up to one. This occurs on grid box", ji
          WRITE(numout,*) "The soiltiles for bare soil, short and tall vegetation are :", soiltile(ji,:)
          CALL ipslerr_p (3,'check_veget', &
               "soiltile does not sum-up to one.", "", "")
       ENDIF
    ENDDO
    
  END SUBROUTINE check_veget


  ! ================================================================================================================================
!! SUBROUTINE   : check_loss_gain
!!
!>\BRIEF        
!!
!! DESCRIPTION  :  Check whether the initially proposed changes in veget_max make sense
!! If they do so, accept this change. If the change brings the remaining veget_max below
!! the minimal threshold remove the whole PFT. Further the maps only prescribe the change
!! in nvmap and do not contain any information about the changes in age classes. Convert the
!! changes from nvmap (= age groups) to nvm (= PFT level). If a PFT looses surface area we
!! assumed that each age classes looses the same relative amount of its surface area.
!!
!! RECENT CHANGE(S) : 
!!
!! MAIN OUTPUT VARIABLE(S) :  loss_gain, veget_max_new
!!
!! REFERENCES   : None
!!
!! FLOWCHART    : 
!! \n
  !_ ================================================================================================================================

  SUBROUTINE check_loss_gain(ipts, ivm, veget_max, loss_gain, &
       failed_vegfrac)


  !! 0. Variable and parameter declaration 
    
    !! 0.1 Input variables
    INTEGER, INTENT(in)                                  :: ipts                    !! Number of pixel (unitless)
    INTEGER, INTENT(in)                                  :: ivm                     !! Number of PFT (unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(in)              :: veget_max               !! Cover fraction of a PFT (unitless)

    !! 0.2 Output variables

    !! 0.3 Modified variables
    
    REAL(r_std), DIMENSION(:,:), INTENT(inout)           :: loss_gain               !! The same as delta_veget_cov but distributed of all
                                                                                    !! age classes and thus taking the age-classes into
    LOGICAL, DIMENSION(:,:), INTENT(inout)               :: failed_vegfrac          !! Failed to find a PFT were some residual fraction could be added (true/false)

    !! 0.4 Local variables
    INTEGER(i_std)                                       :: jvm,ivmax,igroup,iyoung !! Indices
    REAL(r_std)                                          :: max_frac                !! Temporary variable to find the PFT with the maximum fraction
!================================================================================================================================

    ! Try to add the mismatch to the bare soil. The new value of 
    ! bare soil should exceed min_vegfrac.
    IF (veget_max(ipts,ibare_sechiba) + loss_gain(ipts,ibare_sechiba) + &
         veget_max(ipts,ivm) + loss_gain(ipts,ivm).GT.min_vegfrac) THEN

       ! Propagate the precision error into the bare soil PFT
       loss_gain(ipts,ibare_sechiba) = loss_gain(ipts,ibare_sechiba) + &
            veget_max(ipts,ivm) + loss_gain(ipts,ivm)

       ! Set the loss to the available veget_max such that following
       ! LCC the vegetation is completely and exactely removed.
       loss_gain(ipts,ivm) = -veget_max(ipts,ivm)

    ELSE

       ! The bare soil fraction itself was too small to add residual
       ! fraction to. Search for another PFT. New fractions can only
       ! be added to the youngest age class. Search for the youngest
       ! age class with the largest remaining fraction after the land
       ! cover change.
       max_frac = zero
       ! Second loop over nvm, note that the index is now jvm. bare
       ! soil was already checked so start at 2.
       DO jvm = 2,nvm
          igroup=agec_group(jvm)
          iyoung=start_index(igroup)
          ! Avoid repeating the test more than once for the different
          ! age classes of the same age class group.
          IF (jvm.NE.iyoung) CYCLE
          IF (veget_max(ipts,iyoung) + loss_gain(ipts,iyoung) + &
               veget_max(ipts,ivm) + loss_gain(ipts,ivm) .GT. min_vegfrac .AND. &
               veget_max(ipts,iyoung) + loss_gain(ipts,iyoung) + &
               veget_max(ipts,ivm) + loss_gain(ipts,ivm) .GT. max_frac) THEN
             ! Found a new maximum value
             max_frac = veget_max(ipts,iyoung) + loss_gain(ipts,iyoung) + &
                  veget_max(ipts,ivm) + loss_gain(ipts,ivm) 
             ivmax = iyoung
          END IF
       END DO ! jvm

       ! Check whether a suitable PFT was found 
       IF (max_frac.GT.zero) THEN

          ! Propagate the precision error into the PFT with the largest
          ! fraction follwing land cover change
          loss_gain(ipts,ivmax) = loss_gain(ipts,ivmax) + &
               veget_max(ipts,ivm) + loss_gain(ipts,ivm)

          ! Set the loss to the available veget_max such that following
          ! LCC the vegetation is completely and exactely removed.
          loss_gain(ipts,ivm) = -veget_max(ipts,ivm)

       ELSE

          ! No suitable vegetated (inc. bare soil) PFT was found. This
          ! could be because all the young age classes are empty (but
          ! there could still be vegetation in the older age classes).
          ! Will try to move the residual fraction into no_fracbio
          ! later in this subroutine
          failed_vegfrac(ipts,ivm) = .TRUE.

       ENDIF

    ENDIF ! Add residual to the bare soil

  END SUBROUTINE check_loss_gain


!================================================================================================================================
!! SUBROUTINE  : calculate_shares
!!
!>\BRIEF       
!!
!! DESCRIPTION : Calculate the shares that are used to merge the pools
!!               and state variables after a land cover change or age 
!!               class change
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): :: share_rec, litter_weight_rec
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_
!================================================================================================================================

  SUBROUTINE calculate_shares(surf_rec, surf_mov, litter_rec, litter_mov, &
       share_rec, litter_weight_rec)
       
    !! Input variables
    REAL(r_std), INTENT(in)                          :: surf_rec                !! coverage fraction of the receiving PFT
    REAL(r_std), INTENT(in)                          :: surf_mov                !! coverage fraction of the moving PFT
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)        :: litter_rec              !! receiving PFT metabolic and structural litter, above and 
                                                                                !! below ground @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)        :: litter_mov              !! moving PFT metabolic and structural litter, above and 
                                                                                !! below ground @tex ($gC m^{-2}$) @endtex

    !! Output variables
    REAL(r_std), INTENT(out)                         :: share_rec               !! Share of the veget_max of the PFT that is receiving
                                                                                !! the other PFT (unitless, 0-1)
    REAL(r_std), DIMENSION(:,:), INTENT(out)         :: litter_weight_rec       !! The fraction of litter of the PFT that is receiving
                                                                                !! the other PFT @tex $-$ @endtex (unitless, 0-1)
    
    !! Local variables
    INTEGER(i_std)                                   :: ilev, ilit              !! Indices
    
!================================================================================================================================

    ! Share of the PFT that will be added to an existing PFT
    share_rec = surf_rec / ( surf_rec + surf_mov)

    ! We also need a scaling factor which includes the litter
    DO ilev=1,nlevs
       DO ilit=1,nlitt
          
          IF(litter_mov(ilit,ilev,icarbon).GE.min_stomate)THEN
       
             ! For the calculation we only need icarbon but it is a good place
             ! to check whether we also have initrogen
             litter_weight_rec(ilit,ilev) = &
                  litter_rec(ilit,ilev,icarbon)*surf_rec/ &
                  (litter_mov(ilit,ilev,icarbon)*surf_mov + &
                  litter_rec(ilit,ilev,icarbon)*surf_rec)

          ELSE

             ! If we don have icarbon, we need to prescribe the litter weight but
             ! this is also a good place to check that we don't have initrogen
             litter_weight_rec(ilit,ilev)=zero

          END IF

       END DO !nlitt
    END DO ! nlevs

  END SUBROUTINE calculate_shares


!================================================================================================================================
!! SUBROUTINE  : move_pft_properties
!!
!>\BRIEF       
!!
!! DESCRIPTION : move pools, fluxes and state variables from one pft to another
!!               after a land cover change or age class change
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): :: long list of pools, fluxes and state variables
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_
!================================================================================================================================

  SUBROUTINE move_pft_properties(ipts, irec, imov, share_rec, &
     litter_weight_rec, circ_class_biomass, fm_change_map, &
     wstress_season_rec_loc, wstress_season_mov_loc, &
     wstress_month_rec_loc, wstress_month_mov_loc, &
     season_drought_legacy_rec, season_drought_legacy_mov, &
     nstress_season_rec_loc, nstress_season_mov_loc, &
     n_input_rec, n_input_mov, &
     n_input_daily_rec, n_input_daily_mov, &
     vegstress_season_rec_loc, vegstress_season_mov_loc, &
     vegstress_month_rec_loc, vegstress_month_mov_loc, &
     vegstress_week_rec_loc, vegstress_week_mov_loc, &
     vegstress_rec_loc, vegstress_mov_loc, &
     humrel_rec_loc, humrel_mov_loc, &
     everywhere_rec_loc, everywhere_mov_loc, &
     PFTpresent_rec_loc, PFTpresent_mov_loc, &
     sugar_load_rec_loc, sugar_load_mov_loc, &
     cn_leaf_min_season_rec_loc, cn_leaf_min_season_mov_loc, &
     cn_leaf_init_2D, &
     Light_Abs_Tot_rec, Light_Abs_Tot_mov, &
     Light_Tran_Tot_rec, Light_Tran_Tot_mov, &
     laieff_isotrop_rec, laieff_isotrop_mov, &
     lm_lastyearmax_rec_loc, lm_lastyearmax_mov_loc, &
     lm_thisyearmax_rec_loc, lm_thisyearmax_mov_loc, &
     age_rec_loc, age_mov_loc, &
     leaf_frac_rec, leaf_frac_mov, &
     leaf_age_rec, leaf_age_mov, &
     light_tran_to_floor_season_rec_loc, light_tran_to_floor_season_mov_loc, &
     qsintveg_rec_loc, qsintveg_mov_loc, &
     us_rec, us_mov, &
     when_growthinit_rec_loc, when_growthinit_mov_loc, &
     gdd_from_growthinit_rec_loc, gdd_from_growthinit_mov_loc, &
     gdd_midwinter_rec_loc, gdd_midwinter_mov_loc, &
     time_hum_min_rec_loc, time_hum_min_mov_loc, &
     hum_min_dormance_rec_loc, hum_min_dormance_mov_loc, &
     gdd_m5_dormance_rec_loc, gdd_m5_dormance_mov_loc, &
     ncd_dormance_rec_loc, ncd_dormance_mov_loc, & 
     ngd_minus5_rec_loc, ngd_minus5_mov_loc, & 
     mean_start_gs_rec_loc, mean_start_gs_mov_loc, &
     maxvegstress_lastyear_rec_loc, maxvegstress_lastyear_mov_loc, &
     maxvegstress_thisyear_rec_loc, maxvegstress_thisyear_mov_loc, &
     minvegstress_lastyear_rec_loc, minvegstress_lastyear_mov_loc, &
     minvegstress_thisyear_rec_loc, minvegstress_thisyear_mov_loc, &
     maxgppweek_lastyear_rec_loc, maxgppweek_lastyear_mov_loc, &
     maxgppweek_thisyear_rec_loc, maxgppweek_thisyear_mov_loc, &
     n_reserve_longterm_rec_loc, n_reserve_longterm_mov_loc, &
     n_reserve_balance_rec_loc, n_reserve_balance_mov_loc, &
     maxfpc_lastyear_rec_loc, maxfpc_lastyear_mov_loc, &
     maxfpc_thisyear_rec_loc, maxfpc_thisyear_mov_loc, &
     turnover_longterm_rec, turnover_longterm_mov, &
     dead_leaves_rec, dead_leaves_mov, &
     grow_season_len_rec_loc, grow_season_len_mov_loc, &
     KF_rec_loc, KF_mov_loc, &
     atm_to_bm_rec, atm_to_bm_mov, &
     npp_longterm_rec_loc, npp_longterm_mov_loc, &
     croot_longterm_rec_loc, croot_longterm_mov_loc, &
     gpp_daily_rec_loc, gpp_daily_mov_loc, &
     resp_maint_rec_loc, resp_maint_mov_loc, &
     resp_growth_rec_loc, resp_growth_mov_loc, &
     resp_hetero_rec_loc, resp_hetero_mov_loc, &
     npp_daily_rec_loc, npp_daily_mov_loc, &
     rue_longterm_rec_loc, rue_longterm_mov_loc, &
     leaching_daily_rec, leaching_daily_mov, &
     emission_daily_rec, emission_daily_mov, &
     gpp_week_rec_loc, gpp_week_mov_loc, &
     resp_maint_week_rec_loc, resp_maint_week_mov_loc, &
     gpp_year_rec_loc, gpp_year_mov_loc, &
     gpp_decade_rec_loc, gpp_decade_mov_loc, &
     mai_rec_loc, mai_mov_loc, &
     pai_rec_loc, pai_mov_loc, &
     mai_count_rec_loc, mai_count_mov_loc, &
     previous_wood_volume_rec_loc, previous_wood_volume_mov_loc, &
     age_stand_rec_loc, age_stand_mov_loc, &
     last_cut_rec_loc, last_cut_mov_loc, &
     k_latosa_adapt_rec_loc, k_latosa_adapt_mov_loc, &
     litter_rec, litter_mov, &
     bm_to_litter_rec, bm_to_litter_mov, &
     tree_bm_to_litter_rec, tree_bm_to_litter_mov, &
     leaf_age_crit_rec_loc, leaf_age_crit_mov_loc, &
     leaf_classes_rec_loc, leaf_classes_mov_loc, &
     co2_fire_rec_loc, co2_fire_mov_loc, &
     litterfuel_rec, litterfuel_mov, &
     turnover_daily_rec, turnover_daily_mov, &
     soil_n_min_rec, soil_n_min_mov, &
     p_O2_rec_loc, p_O2_mov_loc, &
     bact_rec_loc, bact_mov_loc, &
     deepSOM_a_rec, deepSOM_a_mov, &
     deepSOM_s_rec, deepSOM_s_mov, &
     deepSOM_p_rec, deepSOM_p_mov, &
     som_rec, som_mov, &
     lignin_struc_rec, lignin_struc_mov, &
     lignin_wood_rec, lignin_wood_mov, &
     lignin_snag_rec, lignin_snag_mov, &
     forest_managed_rec_loc, forest_managed_mov_loc, &
     plant_status_rec_loc, plant_status_mov_loc, &
     kill_vessels_rec, kill_vessels_mov, &
     vessel_loss_previous_rec, vessel_loss_previous_mov, &
     biomass_init_drought_rec, biomass_init_drought_mov, &
     count_daylight_rec_loc, count_daylight_mov_loc, &
     CN_som_litter_longterm_rec, CN_som_litter_longterm_mov, &
     matrixV_rec, matrixV_mov, &
     matrixA_rec, matrixA_mov, &
     vectorU_rec, vectorU_mov, &
     vectorB_rec, vectorB_mov, & 
     gap_area_save_rec, gap_area_save_mov, &
     total_ba_init_rec_loc, total_ba_init_mov_loc, keep_source)

    !! Input variables
    ! Guillaume M. -- PROGRESSIVE_HARVEST: this routine implements a WHOLE move -- the
    ! receiver gets the share_rec-weighted blend and the source is systematically reset.
    ! A partial move needs the opposite on the source side: the source keeps its INTENSIVE
    ! values (per m2) while its area shrinks. keep_source=.TRUE. skips every source reset;
    ! absent or .FALSE. reproduces the historical behaviour, hence bit-neutral.
    LOGICAL, OPTIONAL, INTENT(in)                      :: keep_source           !! .TRUE. = partial move, leave the source untouched
    INTEGER(i_std), INTENT(in)                         :: ipts                  !! pixel number
    INTEGER(i_std), INTENT(in)                         :: imov                  !! PFT number of the PFT that is moved to another PFT
    INTEGER(i_std), INTENT(in)                         :: irec                  !! PFT number of the PFT were biomass will be added
    REAL(r_std), INTENT(in)                            :: share_rec             !! Share of the veget_max of the PFT that is receiving 
                                                                                !! the other PFT (unitless, 0-1)
    REAL(r_std), DIMENSION(:,:), INTENT(in)            :: litter_weight_rec     !! The fraction of litter of the receiving PFT @tex $-$ @endtex
    REAL(r_std),DIMENSION(:,:), INTENT(in)             :: cn_leaf_init_2D       !! initial leaf C/N ratio
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(in)      :: circ_class_biomass    !! Biomass components of the model tree within a 
                                                                                !! circumference class @tex $(g C ind^{-1})$ @endtex
    INTEGER(i_std), DIMENSION(:,:), INTENT(in)         :: fm_change_map         !! A map which gives the desired FM strategy when
                                                                                !! the PFT will be replanted after a clearcut.
                                                                                !! (1-nvm,unitless)

    !! Modified variables
    !!
    REAL(r_std), INTENT(inout)                         :: wstress_season_rec_loc, &    !! Water stress factor, based on hum_rel_daily (unitless, 0-1)
                                                          wstress_season_mov_loc 
    REAL(r_std),INTENT(inout)                          :: wstress_month_rec_loc, &     !! Water stress factor, based on hum_rel_daily (unitless, 0-1)
                                                          wstress_month_mov_loc  
    REAL(r_std), DIMENSION(:), INTENT(inout)           :: season_drought_legacy_rec, &  !! mean growing season moisture availability
                                                          season_drought_legacy_mov 
    REAL(r_std), INTENT(inout)                         :: nstress_season_rec_loc, &     !! N-related seasonal stress (used for allocation)
                                                          nstress_season_mov_loc
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: n_input_rec, &                !! Nitrogen inputs into the soil (gN/m**2/timestep)
                                                          n_input_mov           
    REAL(r_std), DIMENSION(:), INTENT(inout)           :: n_input_daily_rec, &          !! nitrogen inputs into the soil  (gN/m**2/day)
                                                          n_input_daily_mov             !! NH4 and NOX from the atmosphere, NH4 from BNF,   
                                                                                        !! agricultural fertiliser as NH4/NO3 n_input_daily(ipts,imov,:) = zero   
    REAL(r_std), INTENT(inout)                         :: vegstress_season_rec_loc, &   !! Mean growingseason moisture availability
                                                          vegstress_season_mov_loc      !! (0 to 1, unitless)  
    REAL(r_std), INTENT(inout)                         :: vegstress_month_rec_loc, &    !! "Monthly" moisture availability (0 to 1,unitless) 
                                                          vegstress_month_mov_loc
    REAL(r_std), INTENT(inout)                         :: vegstress_week_rec_loc, &     !! "Weekly" moisture availability (0 to 1, unitless)
                                                          vegstress_week_mov_loc
    REAL(r_std), INTENT(inout)                         :: vegstress_rec_loc, &          !! Relative soil moisture (0-1, unitless)
                                                          vegstress_mov_loc     
    REAL(r_std), INTENT(inout)                         :: humrel_rec_loc, &             !! Relative humidity. Not used in stomate 
                                                          humrel_mov_loc                !! (needed in age_class_distr)
    REAL(r_std), INTENT(inout)                         :: everywhere_rec_loc, &         !! is the PFT everywhere in the grid box or
                                                          everywhere_mov_loc            !! very localized (after its introduction) (?)
    LOGICAL, INTENT(inout)                             :: PFTpresent_rec_loc, &         !! Tab indicating which PFTs are present in each pixel
                                                          PFTpresent_mov_loc    
    REAL(r_std), INTENT(inout)                         :: sugar_load_rec_loc, &         !! Relative sugar loading of the labile pool (unitless)
                                                          sugar_load_mov_loc   
    REAL(r_std), INTENT(inout)                         :: cn_leaf_min_season_rec_loc, & !! Seasonal min CN ratio of leaves
                                                          cn_leaf_min_season_mov_loc
    REAL(r_std), DIMENSION (:), INTENT (inout)         :: Light_Abs_Tot_rec, &          !! Absorbed radiation per level for photosynthesis
                                                          Light_Abs_Tot_mov     
    REAL(r_std), DIMENSION (:), INTENT (inout)         :: Light_Tran_Tot_rec, &         !! Transmitted radiation per level for photosynthesis
                                                          Light_Tran_Tot_mov   
    REAL(r_std), DIMENSION (:), INTENT (inout)         :: laieff_isotrop_rec, &         !! Effective LAI
                                                          laieff_isotrop_mov   
    REAL(r_std), INTENT(inout)                         :: lm_lastyearmax_rec_loc, &     !! last year's maximum leaf mass for each PFT 
                                                          lm_lastyearmax_mov_loc        !! @tex ($gC m^{-2}$) @endtex
    REAL(r_std), INTENT(inout)                         :: lm_thisyearmax_rec_loc, &     !! this year's maximum leaf mass, for each 
                                                          lm_thisyearmax_mov_loc        !! PFT @tex ($gC m^{-2}$) @endtex
    REAL(r_std), INTENT(inout)                         :: age_rec_loc, &                !! mean age (years)
                                                          age_mov_loc        
    REAL(r_std), DIMENSION(:), INTENT(inout)           :: leaf_frac_rec, &              !! fraction of leaves in leaf age class (unitless;0-1)
                                                          leaf_frac_mov         
    REAL(r_std), DIMENSION(:), INTENT(inout)           :: leaf_age_rec, &               !! Leaf age (days)
                                                          leaf_age_mov       
    REAL(r_std), INTENT(inout)                         :: light_tran_to_floor_season_rec_loc, &!! Mean seasonal fraction of light transmitted
                                                          light_tran_to_floor_season_mov_loc   !! (unitless, 0-1)
    REAL(r_std), INTENT(inout)                         :: qsintveg_rec_loc, &           !! Water on vegetation due to interception
                                                          qsintveg_mov_loc              !! @tex $(kg m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: us_rec, &                     !! Water stress index for transpiration
                                                          us_mov                        !! (by soil layer and PFT) (0-1, unitless)
    REAL(r_std), INTENT(inout)                         :: when_growthinit_rec_loc, &    !! How many days ago was the beginning of
                                                          when_growthinit_mov_loc       !! the growing season (days)
    REAL(r_std), INTENT(inout)                         :: gdd_from_growthinit_rec_loc, &!! growing degree days, since growthinit
                                                          gdd_from_growthinit_mov_loc   !! for crops
    REAL(r_std), INTENT(inout)                         :: gdd_midwinter_rec_loc, &      !! Growing degree days (K), since midwinter (for phenology)
                                                          gdd_midwinter_mov_loc         !! - this is written to the history files
    REAL(r_std), INTENT(inout)                         :: time_hum_min_rec_loc, &       !! Time elapsed since strongest moisture
                                                          time_hum_min_mov_loc          !! availability (days)
    REAL(r_std), INTENT(inout)                         :: hum_min_dormance_rec_loc, &   !! minimum moisture during dormance (0-1, unitless)
                                                          hum_min_dormance_mov_loc
    REAL(r_std), INTENT(inout)                         :: gdd_m5_dormance_rec_loc, &    !! Growing degree days (K), threshold -5 deg C
                                                          gdd_m5_dormance_mov_loc       !! (for phenology) 
    REAL(r_std), INTENT(inout)                         :: ncd_dormance_rec_loc, &       !! Number of chilling days (days), since
                                                          ncd_dormance_mov_loc          !! leaves were lost (for phenology)
    REAL(r_std), INTENT(inout)                         :: ngd_minus5_rec_loc, &         !! Number of growing days (days), threshold
                                                          ngd_minus5_mov_loc            !! -5 deg C (for phenology)
    REAL(r_std), INTENT(inout)                         :: mean_start_gs_rec_loc, &      !! mean growing season starting day for
                                                          mean_start_gs_mov_loc         !! deciduous PFTs (doy) 
    REAL(r_std), INTENT(inout)                         :: maxvegstress_lastyear_rec_loc, &!! last year's maximum moisture availability
                                                          maxvegstress_lastyear_mov_loc 
    REAL(r_std), INTENT(inout)                         :: maxvegstress_thisyear_rec_loc, &!! this year's maximum moisture availability
                                                          maxvegstress_thisyear_mov_loc 
    REAL(r_std), INTENT(inout)                         :: minvegstress_lastyear_rec_loc, &!! last year's minimum moisture availability
                                                          minvegstress_lastyear_mov_loc  
    REAL(r_std), INTENT(inout)                         :: minvegstress_thisyear_rec_loc, &!! this year's minimum moisture availability
                                                          minvegstress_thisyear_mov_loc
    REAL(r_std), INTENT(inout)                         :: maxgppweek_lastyear_rec_loc, &!! last year's maximum weekly GPP
                                                          maxgppweek_lastyear_mov_loc 
    REAL(r_std), INTENT(inout)                         :: maxgppweek_thisyear_rec_loc, &!! this year's maximum weekly GPP
                                                          maxgppweek_thisyear_mov_loc
    REAL(r_std), INTENT(inout)                         :: n_reserve_longterm_rec_loc, & !! "Long term" (default 3 years) actual to
                                                          n_reserve_longterm_mov_loc    !! potential N reserve pool (0-1, unitless)
    REAL(r_std), INTENT(inout)                         :: n_reserve_balance_rec_loc, &  !! Actual to potential N reserve pool
                                                          n_reserve_balance_mov_loc     !! (unitless)
    REAL(r_std), INTENT(inout)                         :: maxfpc_lastyear_rec_loc, &    !! Last year's maximum foliage projected
                                                          maxfpc_lastyear_mov_loc       !! coverage for each natural PFT 
                                                                                        !! @tex $(m^2 m^{-2})$ @endtex 
    REAL(r_std), INTENT(inout)                         :: maxfpc_thisyear_rec_loc, &    !! This year's maximum foliage projected
                                                          maxfpc_thisyear_mov_loc       !! coverage for each natural PFT
                                                                                        !! @tex $(m^2 m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: turnover_longterm_rec, &      !! "Long term" turnover rate
                                                          turnover_longterm_mov         !! @tex $(gC m^{-2} year^{-1})$ @endtex
    REAL(r_std), DIMENSION(:), INTENT(inout)           :: dead_leaves_rec, &            !! Dead leaves on ground, per PFT, metabolic and
                                                          dead_leaves_mov               !! structural @tex $(gC m^{-2})$ @endtex
    REAL(r_std), INTENT(inout)                         :: grow_season_len_rec_loc, &    !! growing season length in days for
                                                          grow_season_len_mov_loc       !! deciduous PFTs 
    REAL(r_std), INTENT(inout)                         :: KF_rec_loc, &                 !! Scaling factor to convert sapwood mass
                                                          KF_mov_loc                    !! into leaf mass (m)
    REAL(r_std), DIMENSION(:), INTENT(inout)           :: atm_to_bm_rec, &              !! CO2 and N taken from the atmosphere to get C
                                                          atm_to_bm_mov                 !! to create the seedlings
                                                                                        !! @tex (gC.m^{-2}dt^{-1})$ @endtex
    REAL(r_std), INTENT(inout)                         :: npp_longterm_rec_loc, &       !! "Long term" mean yearly primary productivity
                                                          npp_longterm_mov_loc          !! @tex $(m^{-2})$ @endtex
    REAL(r_std), INTENT(inout)                         :: croot_longterm_rec_loc, &     !! "long term" root carbon mass
                                                          croot_longterm_mov_loc        !! @tex ($gC m^{-2}) @endtex
    REAL(r_std), INTENT(inout)                         :: gpp_daily_rec_loc, &          !! Daily gross primary productivity
                                                          gpp_daily_mov_loc             !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), INTENT(inout)                         :: resp_maint_rec_loc, &         !! Maintenance respiration
                                                          resp_maint_mov_loc            !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), INTENT(inout)                         :: resp_growth_rec_loc, &        !! Growth respiration
                                                          resp_growth_mov_loc           !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), INTENT(inout)                         :: resp_hetero_rec_loc, &        !! Heterotrophic respiration
                                                          resp_hetero_mov_loc           !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), INTENT(inout)                         :: npp_daily_rec_loc, &          !! Net primary productivity
                                                          npp_daily_mov_loc             !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), INTENT(inout)                         :: rue_longterm_rec_loc, &       !! Longterm radiation use efficiency
                                                          rue_longterm_mov_loc          !! (??units??)
    REAL(r_std), DIMENSION(:), INTENT(inout)           :: leaching_daily_rec, &         !! mineral nitrogen leached from the soil
                                                          leaching_daily_mov            !! (gN/m**2/day)
    REAL(r_std), DIMENSION(:), INTENT(inout)           :: emission_daily_rec, &         !! volatile losses of nitrogen
                                                          emission_daily_mov            !! (gN/m**2/day)
    REAL(r_std), INTENT(inout)                         :: gpp_week_rec_loc, &           !! Weekly gross primary productivity
                                                          gpp_week_mov_loc      
    REAL(r_std), INTENT(inout)                         :: resp_maint_week_rec_loc, &    !! Weekly maintenance respiration
                                                          resp_maint_week_mov_loc 
    REAL(r_std), INTENT(inout)                         :: gpp_year_rec_loc, &           !! "annual" GPP @tex ($gC m^{-2} day^{-1}$)@endtex
                                                          gpp_year_mov_loc      
    REAL(r_std), INTENT(inout)                         :: gpp_decade_rec_loc, &         !! "decadal" GPP @tex ($gC m^{-2} day^{-1}$)@endtex
                                                          gpp_decade_mov_loc
    REAL(r_std), INTENT(inout)                         :: mai_rec_loc, &                !! The mean annual increment
                                                          mai_mov_loc                   !! @tex $(m**3 / m**2 / year)$ @endtex
    REAL(r_std), INTENT(inout)                         :: pai_rec_loc, &                !! The period annual increment
                                                          pai_mov_loc                   !! @tex $(m**3 / m**2 / year)$ @endtex
    INTEGER(i_std), INTENT(inout)                      :: mai_count_rec_loc, &          !! The number of times we've calculated the
                                                          mai_count_mov_loc             !! volume increment for a stand
    REAL(r_std), INTENT(inout)                         :: previous_wood_volume_rec_loc, &!! The volume of the tree trunks
                                                          previous_wood_volume_mov_loc   !! in a stand for the previous year.
                                                                                         !! @tex $(m**3 / m**2 )$ @endtex
    INTEGER(i_std), INTENT(inout)                      :: age_stand_rec_loc, &          !! Age of the forest stand (years)
                                                          age_stand_mov_loc      
    INTEGER(i_std), INTENT(inout)                      :: last_cut_rec_loc, &           !! Years since last thinning (years)
                                                          last_cut_mov_loc     
    REAL(r_std), INTENT(inout)                         :: k_latosa_adapt_rec_loc, &     !! Leaf to sapwood area adapted for long
                                                          k_latosa_adapt_mov_loc        !! term water stress (m)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: litter_rec, &                 !! metabolic and structural litter, above and
                                                          litter_mov                    !! below ground @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: bm_to_litter_rec, &           !! Biomass transfer to litter
                                                          bm_to_litter_mov              !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: tree_bm_to_litter_rec, &      !! Transfer of tree biomass to litter
                                                          tree_bm_to_litter_mov         !! @tex$(gC m^{-2} dtslow^{-1})$@endtex
    REAL(r_std), INTENT(inout)                         :: leaf_age_crit_rec_loc, &      !! critical leaf age (days)
                                                          leaf_age_crit_mov_loc 
    REAL(r_std), INTENT(inout)                         :: leaf_classes_rec_loc, &       !! width of each leaf age class (days)
                                                          leaf_classes_mov_loc  
    REAL(r_std), INTENT(inout)                         :: co2_fire_rec_loc, &           !! Carbon emitted into the atmosphere by
                                                          co2_fire_mov_loc              !! fire (living and dead biomass)
                                                                                        !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: litterfuel_rec, &             !! Dead litter fuel above ground. (gC m^{-2})
                                                          litterfuel_mov        
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: turnover_daily_rec, &         !! Transfer of litter to som
                                                          turnover_daily_mov   
    REAL(r_std), DIMENSION(:), INTENT(inout)           :: soil_n_min_rec, &             !! mineral nitrogen in the soil (gN/m**2)
                                                          soil_n_min_mov                !! (first index=npts, second index=nvm,
                                                                                        !! third index=nnspec) 
    REAL(r_std), INTENT(inout)                         :: p_O2_rec_loc, &               !! partial pressure of oxigen in the soil (hPa)
                                                          p_O2_mov_loc          
    REAL(r_std), INTENT(inout)                         :: bact_rec_loc, &               !! denitrifier biomass (gC/m**2)
                                                          bact_mov_loc          
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: deepSOM_a_rec, &              !! Soil carbon discretized with depth active
                                                          deepSOM_a_mov                 !! (g/m**3) 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: deepSOM_s_rec, &              !! Soil carbon discretized with depth slow
                                                          deepSOM_s_mov                 !! (g/m**3) 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: deepSOM_p_rec, &              !! Soil carbon discretized with depth passive
                                                          deepSOM_p_mov                 !! (g/m**3) 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: som_rec, &                    !! carbon pool: active, slow, or passive
                                                          som_mov                       !! @tex ($gC m^{-2}$) @endtex 
    REAL(r_std), DIMENSION(:), INTENT(inout)           :: lignin_struc_rec, &           !! ratio Lignine/Carbon in structural litter,
                                                          lignin_struc_mov              !! above and below ground
    REAL(r_std), DIMENSION(:), INTENT(inout)           :: lignin_wood_rec, &            !! ratio Lignine/Carbon in woody litter,
                                                          lignin_wood_mov               !! above and below ground
    REAL(r_std), DIMENSION(:), INTENT(inout)           :: lignin_snag_rec, &            !! ratio Lignine/Carbon in snag litter,
                                                          lignin_snag_mov               !! above and below ground
    INTEGER(i_std), INTENT(inout)                      :: forest_managed_rec_loc, &     !! forest management flag (is the forest
                                                          forest_managed_mov_loc        !! being managed?)
    REAL(r_std), INTENT(inout)                         :: plant_status_rec_loc, &       !! Growth and phenological status of the plant
                                                          plant_status_mov_loc          !! is tatus = Phases defined in constantes
    LOGICAL, DIMENSION(:), INTENT(inout)               :: kill_vessels_rec, &           !! Flag to kill vessels at the end of the day
                                                          kill_vessels_mov              !! following embolism 
    REAL(r_std), DIMENSION(:), INTENT(inout)           :: vessel_loss_previous_rec, &   !! Proportion of conductivity lost due
                                                          vessel_loss_previous_mov      !! to cavitation, accumulated on the 
                                                                                        !! previous day (unitless)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: biomass_init_drought_rec, &   !! Biomass of heartwood or sapwood before
                                                          biomass_init_drought_mov      !! onset of drought 
    REAL(r_std), INTENT(inout)                         :: count_daylight_rec_loc, &     !! Number of time steps dt_radia during daylight
                                                          count_daylight_mov_loc
    REAL(r_std), DIMENSION(:), INTENT(inout)           :: CN_som_litter_longterm_rec, & !! Longterm CN ratio of litter and
                                                          CN_som_litter_longterm_mov    !! som pools (gC/gN)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: matrixA_rec, &                !! Matrix containing the fluxes between the
                                                          matrixA_mov                   !! carbon pools per sechiba time step
                                                                                        !! @tex $(gC.m^2.day^{-1})$@endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: matrixV_rec, &                !! Matrix containing the accumulated values
                                                          matrixV_mov                   !! of matrixA
    REAL(r_std), DIMENSION(:), INTENT(inout)           :: vectorB_rec, &                !! Vector containing the litter increase per
                                                          vectorB_mov                   !! sechiba time step @tex $(gCm^{-2})$ @endtex
    REAL(r_std), DIMENSION(:), INTENT(inout)           :: vectorU_rec, &                !! Matrix containing the accumulated values
                                                          vectorU_mov                   !! of VectorB
    REAL(r_std), DIMENSION(:), INTENT(inout)           :: gap_area_save_rec, &          !! Stand gap created by more than 30% basal
                                                          gap_area_save_mov             !! area loss per year
    REAL(r_std), INTENT(inout)                         :: total_ba_init_rec_loc, &      !! Total basal area per pft saved at the start
                                                          total_ba_init_mov_loc         !! of the year (m^{2}/m^{-2})

    !! Local variables
    INTEGER(i_std)                                     :: ilev                     !! Indices
    LOGICAL                                            :: wipe_src                 !! PROGRESSIVE_HARVEST: reset the source after the move?

!================================================================================================================================

    ! Guillaume M. -- Whole move (the historical and overwhelmingly common case) unless the
    ! caller asks for a partial one.
    wipe_src = .TRUE.
    IF (PRESENT(keep_source)) wipe_src = .NOT. keep_source

    !+++CHECK+++
    ! Not sure why we merge wstress here.  wstress is recalculated
    ! everyday from vegstress_season.
    wstress_season_rec_loc = &
         share_rec * wstress_season_rec_loc + &
         wstress_season_mov_loc * (un - share_rec)
    IF (wipe_src) wstress_season_mov_loc = zero
    !-
    ! Initial value = 1
    wstress_month_rec_loc = share_rec * wstress_month_rec_loc + &
         (un - share_rec) * wstress_month_mov_loc
    IF (wipe_src) wstress_month_mov_loc = zero
    !-
    ! Initial value = 1
    nstress_season_rec_loc = &
         share_rec * nstress_season_rec_loc + &
         nstress_season_mov_loc * (un - share_rec)
    IF (wipe_src) nstress_season_mov_loc = zero
    !++++++++++++
    season_drought_legacy_rec(:) = &
         share_rec * season_drought_legacy_rec(:) + &
         season_drought_legacy_mov(:) * (un - share_rec)
    IF (wipe_src) season_drought_legacy_mov(:) = zero
    !-
    ! n_inputs. Not sure this is needed. age_class_distr is done the
    ! last day of the year and the new N-input maps are read the first
    ! day of the next year. They are in here to simplify the consistency
    ! check that can be used to confirm age_class_distr does not affect
    ! the restartability of the model
    n_input_rec(:,:) = share_rec * n_input_rec(:,:) + &
         (un - share_rec) * n_input_mov(:,:)
    IF (wipe_src) n_input_mov(:,:) = zero
    !-
    n_input_daily_rec(:)  = share_rec * n_input_daily_rec(:) + &
         (un - share_rec) * n_input_daily_mov(:)
    IF (wipe_src) n_input_daily_mov(:) = zero
    !-
    ! Initial value = 1
    vegstress_season_rec_loc = &
         share_rec * vegstress_season_rec_loc + &
         vegstress_season_mov_loc * &
         (un - share_rec)
    IF (wipe_src) vegstress_season_mov_loc = zero
    !-
    vegstress_month_rec_loc = share_rec * &
         vegstress_month_rec_loc + &
         (un - share_rec) * vegstress_month_mov_loc
    IF (wipe_src) vegstress_month_mov_loc = zero
    !-
    vegstress_week_rec_loc = share_rec * &
         vegstress_week_rec_loc + &
         (un - share_rec) * vegstress_week_mov_loc
    IF (wipe_src) vegstress_week_mov_loc = zero
    !-
    vegstress_rec_loc = share_rec * vegstress_rec_loc + &
         (un - share_rec) * vegstress_mov_loc
    IF (wipe_src) vegstress_mov_loc = zero
    !-
    humrel_rec_loc = share_rec * humrel_rec_loc + &
         (un - share_rec) * humrel_mov_loc
    IF (wipe_src) humrel_mov_loc = zero
    !-
    everywhere_rec_loc = MAX(everywhere_mov_loc, &
         everywhere_rec_loc)
    IF (wipe_src) everywhere_mov_loc = zero
    !-
    PFTpresent_rec_loc = PFTpresent_mov_loc
    IF (wipe_src) PFTpresent_mov_loc = .FALSE.
    !-
    ! Initial value = 1
    sugar_load_rec_loc = share_rec * sugar_load_rec_loc + &
         (un-share_rec) * sugar_load_mov_loc
    IF (wipe_src) sugar_load_mov_loc = zero
    !-
    ! Initial value = cn_leaf_init_2D(ipts,ivm)
    cn_leaf_min_season_rec_loc = &
         share_rec * cn_leaf_min_season_rec_loc + &
         cn_leaf_min_season_mov_loc * (un - share_rec) 
    ! Initialize. cn_leaf_min_season is only calculated the 
    ! day AFTER phenology. We need an initial value if we
    ! want to (re)grow this PFT.
    IF (wipe_src) cn_leaf_min_season_mov_loc = zero
    !-
    !-
    Light_Abs_Tot_rec(:) = Light_Abs_Tot_rec(:) * &
         share_rec + Light_Abs_Tot_mov(:) * (un - share_rec)
    IF (wipe_src) Light_Abs_Tot_mov(:) = zero
    !-
    Light_Tran_Tot_rec(:) = Light_Tran_Tot_rec(:) * &
         share_rec + Light_Tran_Tot_mov(:) * (un - share_rec)
    IF (wipe_src) Light_Tran_Tot_mov(:) = zero
    !-
    laieff_isotrop_rec(:) = laieff_isotrop_rec(:) * &
         share_rec + laieff_isotrop_mov(:) * (un - share_rec)
    IF (wipe_src) laieff_isotrop_mov(:) = zero
    !-
    lm_lastyearmax_rec_loc = share_rec * &
         lm_lastyearmax_rec_loc + &
         (un - share_rec) * lm_lastyearmax_mov_loc
    IF (wipe_src) lm_lastyearmax_mov_loc = zero
    !-
    lm_thisyearmax_rec_loc = share_rec * &
         lm_thisyearmax_rec_loc + &
         (un - share_rec) * lm_thisyearmax_mov_loc
    IF (wipe_src) lm_thisyearmax_mov_loc = zero
    !-
    age_rec_loc = share_rec * age_rec_loc + &
         (un - share_rec) * age_mov_loc
    IF (wipe_src) age_mov_loc = zero
    !-
    leaf_frac_rec(:) = share_rec * leaf_frac_rec(:) + &
         (un - share_rec) * leaf_frac_mov(:)
    IF (wipe_src) leaf_frac_mov(:) = zero
    !-
    leaf_age_rec(:) = share_rec * leaf_age_rec(:) + &
         (un - share_rec) * leaf_age_mov(:)
    IF (wipe_src) leaf_age_mov(:) = zero
    !-
    leaf_age_crit_rec_loc = share_rec * leaf_age_crit_rec_loc + &
         (un - share_rec) * leaf_age_crit_mov_loc 
    ! leaf_age_crit should not be set to zero, because if a new establishment occurs
    ! after a disturbance, a zero value can cause NaNs in `vmax`, which may then
    ! propagate into `vcmax`, NUE, GPP, and other related variables.
    !-
    leaf_classes_rec_loc = share_rec * leaf_classes_rec_loc + &
         (un - share_rec) * leaf_classes_mov_loc
    IF (wipe_src) leaf_classes_mov_loc = zero
    !- 
    co2_fire_rec_loc = share_rec * co2_fire_rec_loc + &
         (un - share_rec) * co2_fire_mov_loc
    IF (wipe_src) co2_fire_mov_loc = zero
    !-
    litterfuel_rec(:,:,:) = share_rec * &
         litterfuel_rec(:,:,:) + &
         (un - share_rec) * litterfuel_mov(:,:,:)
    IF (wipe_src) litterfuel_mov(:,:,:) = zero
    !-   
    ! light_tran_to_floor_season is only used for recruitment not for
    ! replanting. After moving the old PFT to the new PFT, the old
    ! PFT will be empty. It might be replanted at a later time step 
    ! but replanting will not make use of light_tran_to_floor_season
    ! and can therefore be set to zero. Recruitment is calculated
    ! after the call to this subroutine so light_tran_to_floor_season
    ! need to be moved to the correct age class. Also it is a seasonal
    ! variable that will be update in stomate.f90.
    ! Initial value = 1
    !-
    light_tran_to_floor_season_rec_loc = share_rec * &
         light_tran_to_floor_season_rec_loc + &
         (un - share_rec) * light_tran_to_floor_season_mov_loc
    IF (wipe_src) light_tran_to_floor_season_mov_loc = zero
    !-
    count_daylight_rec_loc = share_rec * &
         count_daylight_rec_loc + &
         (un - share_rec) * count_daylight_mov_loc
    IF (wipe_src) count_daylight_mov_loc = zero
    !-
    ! qsintveg is used in diffuco before recalculating all other
    ! biomass-related variables. Therefore it has no value after
    ! moving the old PFT to the new PFT. If there are several more
    ! variables like this in sechiba, then it'd better to re-think the
    ! order of the calculations, but for now just we move this
    ! sechiba variable in age class distribution.
    qsintveg_rec_loc = share_rec * qsintveg_rec_loc + &
         (un - share_rec) * qsintveg_mov_loc
    IF (wipe_src) qsintveg_mov_loc = zero
    !-
    us_rec(:,:) = share_rec * us_rec(:,:) + &
         (un - share_rec) * us_mov(:,:)
    IF (wipe_src) us_mov(:,:) = zero
    !-
    when_growthinit_rec_loc = share_rec * &
         when_growthinit_rec_loc + &
         (un - share_rec) * when_growthinit_mov_loc
    IF (wipe_src) when_growthinit_mov_loc = zero
    !-
    gdd_from_growthinit_rec_loc = share_rec * &
         gdd_from_growthinit_rec_loc + &
         (un - share_rec) * gdd_from_growthinit_mov_loc
    IF (wipe_src) gdd_from_growthinit_mov_loc = zero
    !-
    ! Initial value = -9999 (undef)
    gdd_midwinter_rec_loc = share_rec * &
         gdd_midwinter_rec_loc + &
         (un - share_rec) * gdd_midwinter_mov_loc
    IF (wipe_src) gdd_midwinter_mov_loc = zero
    !-
    ! Initial value = -9999 (undef)
    time_hum_min_rec_loc = share_rec * time_hum_min_rec_loc + &
         (un - share_rec) * time_hum_min_mov_loc
    IF (wipe_src) time_hum_min_mov_loc = zero
    !-
    ! Initial value = -9999 (undef)
    hum_min_dormance_rec_loc = share_rec * hum_min_dormance_rec_loc + &
         (un - share_rec) * hum_min_dormance_mov_loc
    IF (wipe_src) hum_min_dormance_mov_loc = zero
    !-
    ! Initial value = -9999 (undef)
    gdd_m5_dormance_rec_loc = share_rec * &
         gdd_m5_dormance_rec_loc + &
         (un - share_rec) * gdd_m5_dormance_mov_loc
    IF (wipe_src) gdd_m5_dormance_mov_loc = zero
    !-
    ! Initial value = -9999 (undef)
    ncd_dormance_rec_loc = share_rec * &
         ncd_dormance_rec_loc + &
         (un - share_rec) * ncd_dormance_mov_loc
    IF (wipe_src) ncd_dormance_mov_loc = zero
    !-
    ngd_minus5_rec_loc = share_rec * ngd_minus5_rec_loc + &
         (un - share_rec) * ngd_minus5_mov_loc
    IF (wipe_src) ngd_minus5_mov_loc = zero
    !-
    mean_start_gs_rec_loc = share_rec * mean_start_gs_rec_loc + &
         (un - share_rec) * mean_start_gs_mov_loc
    IF (wipe_src) mean_start_gs_mov_loc = zero
    !-
    maxvegstress_lastyear_rec_loc = share_rec * maxvegstress_lastyear_rec_loc + &
         (un - share_rec) * maxvegstress_lastyear_mov_loc
    IF (wipe_src) maxvegstress_lastyear_mov_loc = zero
    !-
    maxvegstress_thisyear_rec_loc = share_rec * maxvegstress_thisyear_rec_loc + &
         (un - share_rec) * maxvegstress_thisyear_mov_loc
    IF (wipe_src) maxvegstress_thisyear_mov_loc = zero
    !-
    ! Initial value = 1
    minvegstress_lastyear_rec_loc = share_rec * minvegstress_lastyear_rec_loc + &
         (un - share_rec) * minvegstress_lastyear_mov_loc
    IF (wipe_src) minvegstress_lastyear_mov_loc = zero
    !-
    ! Initial value = large_value
    minvegstress_thisyear_rec_loc = share_rec * minvegstress_thisyear_rec_loc + &
         (un - share_rec) * minvegstress_thisyear_mov_loc
    IF (wipe_src) minvegstress_thisyear_mov_loc = zero
    !-
    maxgppweek_lastyear_rec_loc = share_rec * maxgppweek_lastyear_rec_loc+ &
         (un - share_rec) * maxgppweek_lastyear_mov_loc
    IF (wipe_src) maxgppweek_lastyear_mov_loc = zero
    !-
    maxgppweek_thisyear_rec_loc = share_rec * maxgppweek_thisyear_rec_loc + &
         (un - share_rec) * maxgppweek_thisyear_mov_loc
    IF (wipe_src) maxgppweek_thisyear_mov_loc = zero
    !-
    ! Initial value = 1
    n_reserve_longterm_rec_loc = share_rec * n_reserve_longterm_rec_loc + &
         (un - share_rec) * n_reserve_longterm_mov_loc
    IF (wipe_src) n_reserve_longterm_mov_loc = zero
    !-
    n_reserve_balance_rec_loc = share_rec * n_reserve_balance_rec_loc + &
         (un - share_rec) * n_reserve_balance_mov_loc
    IF (wipe_src) n_reserve_balance_mov_loc = zero
    !-
    maxfpc_lastyear_rec_loc = share_rec * maxfpc_lastyear_rec_loc + &
         (un - share_rec) * maxfpc_lastyear_mov_loc
    IF (wipe_src) maxfpc_lastyear_mov_loc = zero
    !-
    maxfpc_thisyear_rec_loc = share_rec * maxfpc_thisyear_rec_loc + &
         (un - share_rec) * maxfpc_thisyear_mov_loc
    IF (wipe_src) maxfpc_thisyear_mov_loc = zero
    !-
    turnover_longterm_rec(:,:) = share_rec * turnover_longterm_rec(:,:) + &
         (un - share_rec) * turnover_longterm_mov(:,:)
    IF (wipe_src) turnover_longterm_mov(:,:) = zero
    !-
    dead_leaves_rec(:) = share_rec * dead_leaves_rec(:) + &
         (un - share_rec) * dead_leaves_mov(:)
    IF (wipe_src) dead_leaves_mov(:) = zero
    !-
    grow_season_len_rec_loc = share_rec * grow_season_len_rec_loc + &
         (un - share_rec) * grow_season_len_mov_loc
    IF (wipe_src) grow_season_len_mov_loc = zero
    !-
    KF_rec_loc = share_rec * KF_rec_loc + &
         (un - share_rec) * KF_mov_loc
    IF (wipe_src) KF_mov_loc = zero
    !-
    gap_area_save_rec(:) = share_rec * gap_area_save_rec(:) + &
         (un - share_rec) * gap_area_save_mov(:)
    IF (wipe_src) gap_area_save_mov(:) = zero
    !-
    total_ba_init_rec_loc = share_rec * total_ba_init_rec_loc + &
         (un - share_rec) * total_ba_init_mov_loc
    IF (wipe_src) total_ba_init_mov_loc = zero
    !-
    !  Note that atm_to_bm is in gC. m-2 dt-1 
    !  so we should also take the weighted mean (rather than the sum if
    !  this where absolute values).
    atm_to_bm_rec(:) = share_rec * atm_to_bm_rec(:) + &
         (un - share_rec) * atm_to_bm_mov(:)
    IF (wipe_src) atm_to_bm_mov(:) = zero
    !-
    npp_longterm_rec_loc = share_rec * npp_longterm_rec_loc + &
         (un - share_rec) * npp_longterm_mov_loc
    IF (wipe_src) npp_longterm_mov_loc = zero
    !-
    croot_longterm_rec_loc = share_rec * croot_longterm_rec_loc + &
         (un - share_rec) * croot_longterm_mov_loc
    IF (wipe_src) croot_longterm_mov_loc = zero
    !-
    gpp_daily_rec_loc = share_rec * gpp_daily_rec_loc + &
         (un - share_rec) * gpp_daily_mov_loc
    IF (wipe_src) gpp_daily_mov_loc = zero
    !-
    resp_maint_rec_loc = share_rec * resp_maint_rec_loc + &
         (un - share_rec) * resp_maint_mov_loc
    IF (wipe_src) resp_maint_mov_loc = zero
    !-
    resp_growth_rec_loc = share_rec * resp_growth_rec_loc + &
         (un - share_rec) * resp_growth_mov_loc
    IF (wipe_src) resp_growth_mov_loc = zero
    !-
    resp_hetero_rec_loc = share_rec * resp_hetero_rec_loc + &
         (un - share_rec) * resp_hetero_mov_loc
    IF (wipe_src) resp_hetero_mov_loc = zero
    !-
    npp_daily_rec_loc = share_rec * npp_daily_rec_loc + &
         (un - share_rec) * npp_daily_mov_loc
    IF (wipe_src) npp_daily_mov_loc = zero
    !-
    ! Initial value = 1
    rue_longterm_rec_loc = share_rec * rue_longterm_rec_loc + &
         (un - share_rec) * rue_longterm_mov_loc
    ! Value of one is used to indicate that this pft has no calculated 
    ! rue_longterm value yet. This matters in stomate_growth_fun_all.
    IF (wipe_src) rue_longterm_mov_loc = zero
    !-
    leaching_daily_rec(:)  = share_rec * leaching_daily_rec(:) + &
         (un - share_rec) * leaching_daily_mov(:)
    IF (wipe_src) leaching_daily_mov(:) = zero
    !-
    emission_daily_rec(:)  = share_rec * emission_daily_rec(:) + &
         (un - share_rec) * emission_daily_mov(:)
    IF (wipe_src) emission_daily_mov(:) = zero
    !-
    gpp_week_rec_loc = share_rec * gpp_week_rec_loc + &
         (un - share_rec) * gpp_week_mov_loc
    IF (wipe_src) gpp_week_mov_loc = zero
    !-
    resp_maint_week_rec_loc = share_rec * resp_maint_week_rec_loc + &
         (un - share_rec) * resp_maint_week_mov_loc
    IF (wipe_src) resp_maint_week_mov_loc = zero
    !-
    gpp_year_rec_loc = share_rec * gpp_year_rec_loc + &
         (un - share_rec) * gpp_year_mov_loc
    IF (wipe_src) gpp_year_mov_loc = zero
    !-
    gpp_decade_rec_loc = share_rec * gpp_decade_rec_loc + &
         (un - share_rec) * gpp_decade_mov_loc
    IF (wipe_src) gpp_decade_mov_loc = zero
    !-
    pre_indust_ref_gpp(ipts,irec) = share_rec * pre_indust_ref_gpp(ipts,irec) + &
         (un - share_rec) * pre_indust_ref_gpp(ipts,imov)
    IF (wipe_src) pre_indust_ref_gpp(ipts,imov) = zero
    !-
    mai_rec_loc = share_rec * mai_rec_loc + &
         (un - share_rec) * mai_mov_loc
    IF (wipe_src) mai_mov_loc = zero
    !-
    pai_rec_loc = share_rec * pai_rec_loc + &
         (un - share_rec) * pai_mov_loc
    IF (wipe_src) pai_mov_loc = zero
    !-
    ! We need an integer for mai_count.  The use of CEILING here is
    ! arbitrary.  Maybe FLOOR would work better in some cases.
    mai_count_rec_loc = CEILING(share_rec * &
         mai_count_rec_loc + &
         (un - share_rec) * mai_count_mov_loc)
    IF (wipe_src) mai_count_mov_loc = 0
    !-
    previous_wood_volume_rec_loc = share_rec * &
         previous_wood_volume_rec_loc + &
         (un - share_rec) * previous_wood_volume_mov_loc
    IF (wipe_src) previous_wood_volume_mov_loc = zero
    !-
    ! These next two variables are integers, so we use CEILING to make sure
    ! they stay integers.  FLOOR or INT would also be fine.  The use of
    ! CEILING is arbitrary.
    age_stand_rec_loc = CEILING(share_rec * &
         age_stand_rec_loc + &
         (un - share_rec) * age_stand_mov_loc)
    IF (wipe_src) age_stand_mov_loc = zero
    !-
    last_cut_rec_loc = CEILING(share_rec * last_cut_rec_loc + &
         (un - share_rec) * last_cut_mov_loc)
    IF (wipe_src) last_cut_mov_loc = zero
    !-
    ! Initial value k_latosa_min(m)
    k_latosa_adapt_rec_loc = share_rec * &
         k_latosa_adapt_rec_loc + &
         (un - share_rec) * k_latosa_adapt_mov_loc
    IF (wipe_src) k_latosa_adapt_mov_loc = zero
    !-
    litter_rec(:,:,:) =   &
         share_rec * litter_rec(:,:,:) + &
         litter_mov(:,:,:) * &
         (un - share_rec)
    IF (wipe_src) litter_mov(:,:,:) = zero
    !-
    bm_to_litter_rec(:,:) = &
         share_rec * bm_to_litter_rec(:,:) + &
         bm_to_litter_mov(:,:) * &
         (un - share_rec)
    IF (wipe_src) bm_to_litter_mov(:,:) = zero
    !-
    tree_bm_to_litter_rec(:,:) = &
         share_rec * tree_bm_to_litter_rec(:,:) + &
         tree_bm_to_litter_mov(:,:) * &
         (un - share_rec)
    IF (wipe_src) tree_bm_to_litter_mov(:,:) = zero
    !-
    
    !
    turnover_daily_rec(:,:) = share_rec * &
         turnover_daily_rec(:,:) + &
         (un - share_rec) * turnover_daily_mov(:,:)
    IF (wipe_src) turnover_daily_mov(:,:) = zero
    !-              
    soil_n_min_rec(:) = &
         share_rec * soil_n_min_rec(:) + &
         soil_n_min_mov(:) * &
         (un - share_rec)
    IF (wipe_src) soil_n_min_mov(:) = zero
    !-
    ! Initial value = 200
    p_O2_rec_loc = &
         share_rec * p_O2_rec_loc + &
         p_O2_mov_loc * (un - share_rec)
    IF (wipe_src) p_O2_mov_loc = zero
    !-
    ! Initial value = 10
    bact_rec_loc = &
         share_rec * bact_rec_loc + &
         bact_mov_loc * (un - share_rec)
    IF (wipe_src) bact_mov_loc = zero
    !-
    IF (ok_soil_carbon_discretization) THEN
       !-
       deepSOM_a_rec(:,:) =  &
            share_rec * deepSOM_a_rec(:,:) + &
            deepSOM_a_mov(:,:) * (un - share_rec)
       IF (wipe_src) deepSOM_a_mov(:,:) = zero
       !-
       deepSOM_s_rec(:,:) =  &
            share_rec * deepSOM_s_rec(:,:) + &
            deepSOM_s_mov(:,:) * (un - share_rec)
       IF (wipe_src) deepSOM_s_mov(:,:) = zero
       !-
       deepSOM_p_rec(:,:) =  &
            share_rec * deepSOM_p_rec(:,:) + &
            deepSOM_p_mov(:,:) * (un - share_rec)
       IF (wipe_src) deepSOM_p_mov(:,:) = zero
    ELSE
       !-
       som_rec(:,:) =  &
            share_rec * som_rec(:,:) + &
            som_mov(:,:) * (un - share_rec)
       IF (wipe_src) som_mov(:,:) = zero
    END IF
    !-
    ! The new soil&litter pools are the weighted mean of the newly 
    ! established vegetation for that PFT and the soil&litter pools
    ! of the vegetation that already exists in that PFT.
    ! Since it is not only the amount of vegetation present (veget_max) but also
    ! the amount of structural litter (litter) that is important, we have to 
    ! weight by both items here.
    DO ilev=1,nlevs
       lignin_struc_rec(ilev) = &
            litter_weight_rec(istructural,ilev) * &
            lignin_struc_rec(ilev) + &
            (un - litter_weight_rec(istructural,ilev)) * &
            lignin_struc_mov(ilev)
       IF (wipe_src) lignin_struc_mov(ilev) = zero
       !-
       lignin_wood_rec(ilev) = &
            litter_weight_rec(iwoody,ilev) * &
            lignin_wood_rec(ilev) + &
            (un - litter_weight_rec(iwoody,ilev) ) * &
            lignin_wood_mov(ilev)
       IF (wipe_src) lignin_wood_mov(ilev) = zero
       !-
       lignin_snag_rec(ilev) = &
            litter_weight_rec(isnag,ilev) * &
            lignin_snag_rec(ilev) + &
            (un - litter_weight_rec(isnag,ilev) ) * &
            lignin_snag_mov(ilev)
       IF (wipe_src) lignin_snag_mov(ilev) = zero
 
    ENDDO

    ! When merging age classes and accounting for management 
    ! changes, the merged pft is managed the same way as the new
    ! pft. If not, prescribed changes in forest management (as 
    ! part of the species change code) would not be consistently
    ! implemented. However, when simply merging age classes without
    ! species/management changes a conservative approach is
    ! taken and the forest management strategy of the oldest age
    ! class in the merge is used. 
    IF(ok_change_species)THEN
       ! If the old and new pft have a different management 
       ! strategy we will have to decide which strategy we 
       ! follow
       IF (forest_managed_rec_loc .NE. &
            forest_managed_mov_loc) THEN
          IF ((forest_managed_rec_loc.EQ.ifm_none) .OR. &
               (forest_managed_mov_loc).EQ.ifm_none) THEN
             ! We assume that unmanaged sites are well 
             ! protected. Actually they will even increase
             ! in surface area even if the desired fm
             ! strategy is not ifm_none.
             forest_managed_rec_loc = ifm_none 
             IF (wipe_src) forest_managed_mov_loc = ifm_none
          ELSE
             ! If the new old fm type are not equal we will 
             ! use the new prescribed forest management as prepared
             ! in fm_change_map. This should ensure that we can
             ! get the old management back.
             forest_managed_rec_loc = fm_change_map(ipts,irec) 
             forest_managed_mov_loc = forest_managed_mov_loc
             ! Debug
             IF(printlev_loc>=4)THEN
                WRITE(numout,*) 'species change - age class distr'
                WRITE(numout,*) 'fm - old PFT, ',imov, &
                     forest_managed_mov_loc
                WRITE(numout,*) 'fm - new PFT, ',irec, &
                     fm_change_map(ipts,irec)
             ENDIF
             !-
          ENDIF
       ELSE
          ! fm type for the new and old PFT are the same so nothing
          ! should be done. The code was written to make it explicit
          forest_managed_rec_loc = forest_managed_rec_loc
          forest_managed_mov_loc = forest_managed_mov_loc
       ENDIF
    ELSE ! no species change 
       ! Preserve the prescribed forest management.
       ! Do nothing but the code was written to make it explicit
       forest_managed_rec_loc = forest_managed_rec_loc
       forest_managed_mov_loc = forest_managed_mov_loc
    END IF ! species change (through replant)?

    ! Deciding which plant_status to take is a rather complex
    ! issue. We do this at the end of the year which is winter
    ! in the northern hemisphere but summer in the southern
    ! hemisphere. Furthermore part of the phenology is PFT depedent
    ! which implies that different PFTs of the same age_group may
    ! have a different plant status. Some plan_statuses require
    ! other variables to have specific values. For example, isenescent
    ! assumes that are still leaves and roots whereas ibudavailable, 
    ! ibudbreak and idormant assume that there is no root and leaf
    ! biomass. iprescribe and idead assume that there is no biomass
    ! at all. The consequence is that we have to be careful and that
    ! we cannot use a simple rule to decide on the plant_status
    ! of the merged biomass. Note that some conflicts will only 
    ! cause trouble for deciduous trees.
    IF (plant_status_rec_loc.EQ.plant_status_mov_loc)THEN
       ! Doesn't matter which plant_status is taken for the
       ! new pft.
       plant_status_rec_loc = plant_status_mov_loc
    ELSEIF (plant_status_rec_loc.EQ.inone .OR. &
         plant_status_rec_loc.EQ.iprescribe)THEN
       ! This means the irec has no biomass, so take the plant_status
       ! of the imov.
       plant_status_rec_loc = plant_status_mov_loc
    ELSEIF (plant_status_rec_loc.EQ.ibudsavail)THEN
       ! This means that the irec has no leaves/roots yet. Overwrite
       ! this plant_status with the one of the imov just in case
       ! that one already has leaves/roots.
       plant_status_rec_loc = plant_status_mov_loc
    ELSEIF (plant_status_rec_loc.EQ.ibudbreak)THEN
       ! Phenology could crash if leaf and root biomass were already present. 
       ! In earlier implementations, it was assumed that irec could not retain
       ! its plant status 
       ! when the old PFT still contained non-zero leaf or root biomass. 
       !
       ! However, since Phenology is called before any LCC modules, and leaf and root 
       ! biomass are initialized within Phenology, there is no crash risk at that stage. 
       ! Instead, a crash occurred later when the plant status changed from
       ! saplings to idormant during a merge, and then switched to ibudsavail while 
       ! leaf and root biomass were still non-zero. 
       ! Therefore, the following line is commented out to prevent this issue:
!       plant_status(ipts,irec) = plant_status(ipts,imov)
    ELSEIF (plant_status_rec_loc.EQ.icanopy)THEN
       ! This canopy cannot be simply dropped because then
       ! we may cause a conflict between the plant_status and
       ! the biomass pool. The canopy needs to be kept.
    ELSEIF (plant_status_rec_loc.EQ.ipresenescence)THEN
       ! This canopy cannot be simply dropped because then
       ! we may cause a conflict between the plant_status and
       ! the biomass pool. The canopy needs to be kept.
       plant_status_rec_loc = plant_status_rec_loc
    ELSEIF (plant_status_rec_loc.EQ.isenescent)THEN
       ! This canopy cannot be simply dropped because then
       ! we may cause a conflict between the plant_status and
       ! the biomass pool. The canopy needs to be kept.
       plant_status_rec_loc = plant_status_rec_loc
    ELSEIF (plant_status_rec_loc.EQ.idormant)THEN
       ! idormant implies it has no leaves, If the imov has leaves
       ! the plant_status should be set to that of the imov.
       plant_status_rec_loc = plant_status_mov_loc
    ELSEIF (plant_status_rec_loc.EQ.idead)THEN
       ! If the imov has some biomass it cannot be set to
       ! idead. If it doesn have biomass its plant_status is idead 
       ! or iprescribe. Set to the value of the imov. This may
       ! imply that a dead PFT is revived. It may die again in a
       ! couple of years or it may need to be replanted - impossible
       ! to say.
       plant_status_rec_loc = plant_status_mov_loc
    ELSE
       WRITE(numout,*) 'plant_status irec, ',plant_status_rec_loc
       WRITE(numout,*) 'plant_status imov, ',plant_status_mov_loc
       CALL ipslerr_p(3,'move_pft_properties','Undefined case in IF-statement', &
            'Define the missing case in the code','')
    ENDIF
    ! Remember to reset the old pft which will no longer a veget_max.
    ! Guillaume M. -- PROGRESSIVE_HARVEST: under a partial move the source DOES keep a
    ! veget_max, so the premise of this reset does not hold -- leaving inone there would
    ! break the phenology of the surviving stand.
    IF (wipe_src) plant_status_mov_loc = inone
    
    IF (ok_hydrol_arch .AND. ok_vessel_mortality) THEN
                      
       CALL ipslerr_p (3,'age_class_distr',&
            'you are the first to use vessel mortality with age classes', &
            'some code has been proposed but has not yet been tested',&
            'test before using!')
       
       ! Vessel mortality will be stopped when the biomass is moved
       ! from one to another PFT because the variables are not suited
       ! to take weighted means. The drought may continue at the next
       ! time step.
       kill_vessels_rec(:) = .FALSE.
       IF (wipe_src) kill_vessels_mov(:) = .FALSE.
       vessel_loss_previous_rec(:) = zero
       IF (wipe_src) vessel_loss_previous_mov(:) = zero
       biomass_init_drought_rec(:,:,:) = circ_class_biomass(ipts,irec,:,:,:)
       IF (wipe_src) biomass_init_drought_mov(:,:,:) = zero
       
    END IF

    ! Move dynamic parameters         
    IF (test_dynamic_alpha_self_thin) THEN
       alpha_self_thinning(ipts,irec) = share_rec * &
          alpha_self_thinning(ipts,irec) + &
          alpha_self_thinning(ipts,imov) * (un - share_rec)
       ! Guillaume M. -- PROGRESSIVE_HARVEST: the surviving stand must keep its realised
       ! alpha, not fall back on the static MTC reference. The reset used the HARD-CODED
       ! table ref_alpha_self_thin_mtc, while the live reference is the NAMELIST array
       ! ref_alpha_self_thin(ivm), the one stomate_lpj uses to recompute alpha every year.
       !
       ! Guillaume M. -- Two sources for one quantity: the gap was transient (annual
       ! recompute) but a stray alpha reads as a stand anomaly, RDI = N/Nmax(qmd, alpha).
       ! /!\ NOT bit-neutral for configurations overriding REF_ALPHA_SELF_THIN, and that is
       ! intended: the old value was the wrong one.
       IF (wipe_src) alpha_self_thinning(ipts,imov) = ref_alpha_self_thin(imov)
    ENDIF

    pipe_tune2(ipts,irec) = share_rec * pipe_tune2(ipts,irec) + &
         pipe_tune2(ipts,imov) * (un - share_rec)
    IF (wipe_src) pipe_tune2(ipts,imov) = pipe_tune2_mtc(pft_to_mtc(imov))
                    
    ! For the analytic spinup, we need longterm CN ratios
    ! calculated over the spinup period. If the tree grows
    ! and move to an older age class, this is taken care of
    ! by age_class_distr, but if PFT is killed we need to
    ! move CN back to the first age class as well.
    ! Otherwise, we will not get a CN ratio for soil and
    ! litter pools the corresponds to the spinup period.
    IF(spinup_analytic)THEN
       !-
       CN_som_litter_longterm_rec(:)=&
            share_rec * CN_som_litter_longterm_rec(:) + &
            CN_som_litter_longterm_mov(:) * (un - share_rec)
       IF (wipe_src) CN_som_litter_longterm_mov(:) = zero
       !-
       ! If we are doing an analytical spinup, we need to merge
       ! the matrices.  During our spinup, only one age class in a given
       ! PFT will be populated. This implies that the spinup should
       ! make use of LCC, and disturbances (i.e., processes that affect
       ! veget_max). Therefore, we can just move the matrix from one PFT 
       ! to the next.  If there is a chance that there is already biomass 
       ! in the new PFT and things need to be merged, you should rethink 
       ! this
       IF(share_rec.LT.min_stomate)THEN
          !-
          matrixV_rec(:,:)=matrixV_mov(:,:)
          IF (wipe_src) matrixV_mov(:,:) = zero
          !-
          matrixA_rec(:,:)=matrixA_mov(:,:)
          IF (wipe_src) matrixA_mov(:,:) = zero
          !-
          vectorU_rec(:)=vectorU_mov(:)
          IF (wipe_src) vectorU_mov(:) = zero
          !-
          vectorB_rec(:)=vectorB_mov(:)
          IF (wipe_src) vectorB_mov(:) = zero
       ELSE
          WRITE(numout,*) 'Merge problem with spinup.'
          WRITE(numout,*) 'ipts,ipft: ',ipts,imov
          WRITE(numout,*) 'share_rec: ',share_rec
          CALL ipslerr_p (3,'age_class_distr', &
               'Trying to merge spinup matrices into a PFT with existing',&
               'biomass. Stopping because we are not sure what to do.','')
       END IF ! share_rec .LT. min_stomate
    END IF ! spinup_analytic                  

  END SUBROUTINE move_pft_properties

  
  
! ================================================================================================================================
!! SUBROUTINE   : print_lcc_error_message
!!
!>\BRIEF        
!!
!! DESCRIPTION  : prints the values of the vegetation fractions and derived
!!                variables. Useful for debugging.  
!!
!! RECENT CHANGE(S) : 
!!
!! MAIN OUTPUT VARIABLE(S) :  
!!
!! REFERENCES   : None
!!
!! FLOWCHART    : 
!! \n
!_ ================================================================================================================================

  SUBROUTINE print_lcc_error_message(kjpindex, ipts, veget_max_new, veget_max, &
       veget_max_new_group, failed_vegfrac, veget_max_group, delta_veget_cov, &
       loss_gain, flag)

    IMPLICIT NONE
    
    !! 0. Variable and parameter declaration 
    
    !! 0.1 Input variables

    INTEGER, INTENT(in)                                    :: kjpindex                 !! Domain size - number of pixels (unitless, 0-1)
    INTEGER, INTENT(in)                                    :: ipts                     !! Location of current gridsquare
    INTEGER, INTENT(in)                                    :: flag                     !! Indicates what to print
    REAL(r_std), DIMENSION(:,:), INTENT(in)                :: veget_max                !! Cover fraction of a PFT (unitless, 0-1) 
    REAL(r_std), DIMENSION(:,:), INTENT(in)                :: veget_max_new            !! Cover fraction after LCC of a PFT (unitless, 0-1) 
    REAL(r_std),DIMENSION(:,:), INTENT(in)                 :: delta_veget_cov          !! changes in cover fraction of PFT (unitless, 0-1) 
    REAL(r_std), DIMENSION(:,:) , INTENT(in)               :: veget_max_new_group      !! Total veget_max_new of all age clases within a PFT group 
                                                                                       !! (unitless, 0-1)
    REAL(r_std), DIMENSION(:,:) , INTENT(in)               :: veget_max_group          !! Total veget_max of all age clases within a PFT group 
                                                                                       !! (unitless, 0-1)
    REAL(r_std), DIMENSION(:,:), INTENT(in)                :: loss_gain                !! Same as delta_veget_cov but corrected and distributed
                                                                                       !! over all age classes and thus taking the age-classes into 
                                                                                       !! account (unitless, 0-1)
    LOGICAL, DIMENSION(:,:), INTENT(inout)                 :: failed_vegfrac           !! Failed to find a PFT were some
                                                                                       !! residual fraction could be added (true/false)
    
    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables
    INTEGER(i_std)                                         :: ivm,ivma                 !! Indices
!_ ================================================================================================================================

    IF(ipts .LT. 1)THEN
       WRITE(numout,*) 'ERROR: Passing a bad point to print_lcc_error_message.'
       WRITE(numout,*) 'ipts, test_grid: ',ipts, test_grid
       RETURN
    ENDIF

    SELECT CASE(flag)
    CASE(1)

       ! Case 1
       WRITE(numout,*) 'ipts ivma veget_max_group veget_max_new_group'
       DO ivma=1,nvmap
          WRITE(numout,'(2I7,2F14.10)') ipts,ivma,veget_max_group(ipts,ivma),&
               veget_max_new_group(ipts,ivma)
       ENDDO
       WRITE(numout,*) 'SUMS: ',SUM(veget_max_group(ipts,:)),&
            SUM(veget_max_new_group(ipts,:))

    CASE(2)

       ! Case 2
       WRITE(numout,*) 'ivm veget_max veget_max_new'
       DO ivm=1,nvm
          WRITE(numout,'(2I7,2F16.12)') ipts,ivm,veget_max(ipts,ivm),&
               veget_max_new(ipts,ivm)
       ENDDO
       WRITE(numout,*) 'SUMS: ',SUM(veget_max(ipts,:)),&
            SUM(veget_max_new(ipts,:))

    CASE(3)

       ! Case 3
       WRITE(numout,*) 'ivm veget_max_group delta_veget_cov'
       DO ivma=1,nvmap
          WRITE(numout,'(2I7,3F16.12)') ipts,ivma,veget_max_group(ipts,ivma),&
               delta_veget_cov(ipts,ivma)
       ENDDO
       WRITE(numout,*) 'SUMS: ',SUM(veget_max_group(ipts,:)),&
            SUM(delta_veget_cov(ipts,:))
    CASE(4)

       !Case 4
       WRITE(numout,*) 'ipts ivm veget_max loss_gain failed_vegfrac'
       DO ivm=1,nvm
         WRITE(numout,*) ipts, ivm, veget_max(ipts,ivm),&
                loss_gain(ipts,ivm), failed_vegfrac(ipts,ivm)
       ENDDO
       WRITE(numout,*) 'SUMS: ',SUM(veget_max(ipts,:)),&
            SUM(loss_gain(ipts,:))
       CALL flush(numout)

   CASE(5)

       ! Case 5
       WRITE(numout,*) 'ivm agec_group'
       DO ivm=1,nvm
          WRITE(numout,'(2I8)') ivm,agec_group(ivm)
       ENDDO

    CASE default

       ! Default case
       WRITE(numout,*) 'ERROR: Not sure what you want me to print in LCC.'
    END SELECT
     
  END SUBROUTINE print_lcc_error_message


END MODULE sapiens_lcchange
