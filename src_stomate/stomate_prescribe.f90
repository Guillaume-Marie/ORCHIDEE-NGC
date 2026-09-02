! =================================================================================================================================
! MODULE       : stomate_prescribe
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF         Initialize and update density, height, basal area and crown area.
!!
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S)	:
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_stomate/stomate_prescribe.f90 $
!! $Date: 2026-03-20 17:55:41 +0100 (ven. 20 mars 2026) $
!! $Revision: 9442 $
!! \n
!_ ================================================================================================================================

MODULE stomate_prescribe

  ! modules used:
  USE ioipsl_para
  USE xios_orchidee  
  USE time, ONLY : ts_annual_proc
  USE pft_parameters
  USE dynamic_parameters
  USE constantes
  USE function_library,    ONLY: calculate_c0_alloc, cc_to_lai, & 
       cc_to_biomass, check_vegetation_area_point, check_mass_balance_point, &
       wood_to_qmdia, Nmax, sort_circ_class_biomass, lai_to_biomass, &
       biomass_to_lai, wood_to_ba, get_printlev, weibull_class_dist, &
       calculate_rdi_boundaries, cc_grid_max_x

  IMPLICIT NONE

  ! private & public routines

  PRIVATE
  PUBLIC prescribe, prescribe_clear

  ! first call
  LOGICAL, SAVE             :: firstcall_prescribe = .TRUE.
!$OMP THREADPRIVATE(firstcall_prescribe)
  INTEGER(i_std), SAVE       :: printlev_loc       !! Local level of text output for this module
!$OMP THREADPRIVATE(printlev_loc)

CONTAINS

! =================================================================================================================================
!! SUBROUTINE   : prescribe_clear
!!
!>\BRIEF        : Set the firstcall_prescribe flag back to .TRUE. to prepare for the next simulation.
!_=================================================================================================================================

  SUBROUTINE prescribe_clear
    firstcall_prescribe=.TRUE.
  END SUBROUTINE prescribe_clear


!! ================================================================================================================================
!! SUBROUTINE   : prescribe
!!
!>\BRIEF         Works only with static vegetation and agricultural PFTs.
!!               Initialize biomass, density, presence in the first call
!!               and update them in the following.
!!
!! DESCRIPTION (functional, design, flags): \n
!! This module works only with static vegetation and agricultural PFTs.
!! In the first call, initialize density of individuals, biomass,
!! and leaf age distribution to some reasonable value. In the following calls,
!! these variables are updated following a stand replacing disturbance or
!! if the conditions are suitable for recruitment. 
!!
!! To fulfill these purposes, pipe model are used:
!! \latexonly 
!!     \input{prescribe1.tex}
!!     \input{prescribe2.tex}
!! \endlatexonly
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLES(S): ::ind, ::biomass
!!
!! REFERENCES   :
!! - Krinner, G., N. Viovy, et al. (2005). "A dynamic global vegetation model 
!!   for studies of the coupled atmosphere-biosphere system." Global 
!!   Biogeochemical Cycles 19: GB1015, doi:1010.1029/2003GB002199.
!! - Sitch, S., B. Smith, et al. (2003), Evaluation of ecosystem dynamics,
!!   plant geography and terrestrial carbon cycling in the LPJ dynamic 
!!   global vegetation model, Global Change Biology, 9, 161-185.
!! - McDowell, N., Barnard, H., Bond, B.J., Hinckley, T., Hubbard, R.M., Ishii, 
!!   H., Köstner, B., Magnani, F. Marshall, J.D., Meinzer, F.C., Phillips, N., 
!!   Ryan, M.G., Whitehead D. 2002. The relationship between tree height and leaf 
!!   area: sapwood area ratio. Oecologia, 132:12–20.
!! - Novick, K., Oren, R., Stoy, P., Juang, F.-Y., Siqueira, M., Katul, G. 2009. 
!!   The relationship between reference canopy conductance and simplified hydraulic 
!!   architecture. Advances in water resources 32, 809-819.
!! - Ruger et al. 2009, J. of Ecol. doi: 10.1111/j.1365-2745.2009.01552.x  
!!
!! \n
!_ ================================================================================================================================

SUBROUTINE prescribe (&
     ivm,                 veget_max_loc, &
     dt,                  PFTpresent_loc, &
     everywhere_loc,      when_growthinit_loc,            leaf_frac_loc, &
     circ_class_n_loc,    circ_class_biomass_loc,         atm_to_bm_loc,         forest_managed_loc, &
     KF_loc,              plant_status_loc,               age_loc,               npp_longterm_loc, &
     lm_lastyearmax_loc,  longevity_eff_leaf_loc,         longevity_eff_sap_loc, longevity_eff_root_loc, &
     k_latosa_adapt_loc,  light_tran_to_floor_season_loc, species_change_map_loc, &
     cn_leaf_init_2D_loc, bm_sapl_2D_loc,                 new_ind_loc,           pipe_tune2_loc, &
     alpha_self_thin_loc, write_debug,                    tot_res_target_loc,    tot_lab_target_loc, &
     gpp_week,            t2m,                            qmd_init_override,     established_out)

  !! 0. Parameters and variables declaration

  !! 0.1 Input variables
    LOGICAL, INTENT(in)                                    :: write_debug            !! flag to write debug statements in the prescribe subroutine
    INTEGER(i_std), INTENT(in)                             :: ivm                    !! Specific PFT to use for the calculations (unitless)
    INTEGER(i_std), INTENT(in)                             :: forest_managed_loc     !! forest management flag
    REAL(r_std), INTENT(in)                                :: dt                     !! time step (dt_days)   
    REAL(r_std), INTENT(in)                                :: veget_max_loc          !! "maximal" coverage fraction of a PFT 
                                                                                     !! (LAI -> infinity) on ground. May sum to
                                                                                     !! less than unity if the pixel has
                                                                                     !! nobio area. (unitless; 0-1)
    REAL(r_std), INTENT(in)                                :: longevity_eff_root_loc !! Effective root turnover time that accounts
                                                                                     !! waterstress (days)
    REAL(r_std), INTENT(in)                                :: longevity_eff_sap_loc  !! Effective sapwood turnover time that accounts
                                                                                     !! waterstress (days)
    REAL(r_std), INTENT(in)                                :: longevity_eff_leaf_loc !! Effective leaf turnover time that accounts
                                                                                     !! waterstress (days)
    INTEGER(i_std), INTENT(in)                             :: species_change_map_loc !! map which gives the PFT number that each PFT
    REAL(r_std),INTENT(inout)                              :: light_tran_to_floor_season_loc  !! Seasonal fraction of light transmitted to canopy levels
    REAL(r_std), INTENT(in)                                :: cn_leaf_init_2D_loc    !! initial leaf C/N ratio
    REAL(r_std), INTENT(in)                                :: pipe_tune2_loc         !! pipe_tune2 for this ipts and ivm
    REAL(r_std), INTENT(in)                                :: alpha_self_thin_loc    !! alpha_self_thinning for this ipts and ivm
    REAL(r_std), INTENT(in)                                :: gpp_week               !! PFT gross primary productivity
    REAL(r_std), INTENT(in)                                :: t2m                    !! Temperature at 2 meter (K)
    REAL(r_std), INTENT(in), OPTIONAL                      :: qmd_init_override      !! Cold-start desync: per-(point,class) initial qmd (m)
                                                                                     !! to use instead of qmd_init(ivm) at tree establish.
                                                                                     !! >= val_exp/2 or absent => use qmd_init(ivm) (neutral).

  !! 0.2 Output variables
    REAL(r_std), INTENT(out)                               :: new_ind_loc            !! density of new individuals (trees m-2 day-1)
    LOGICAL, INTENT(out), OPTIONAL                         :: established_out        !! Cold-start desync: TRUE if a tree stand was
                                                                                     !! established (establish branch) on this call.
    REAL(r_std), DIMENSION(nelements), INTENT(out)         :: tot_res_target_loc     !! Target for reserve carbon pool (only for grassland) 
                                                                                     !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(nelements), INTENT(out)         :: tot_lab_target_loc     !! Target for labile carbon pool (only for grassland) 
                                                                                     !! @tex $(gC m^{-2})$ @endtex 


  !! 0.3 Modified variables
    LOGICAL, INTENT(inout)                                 :: PFTpresent_loc         !! PFT present (0 or 1)
    REAL(r_std), INTENT(inout)                             :: plant_status_loc       !! Growth and phenological status of the plant
                                                                                     !! Different stati are listed in in constantes_var
    REAL(r_std), INTENT(inout)                             :: everywhere_loc         !! is the PFT everywhere in the grid box or 
                                                                                     !! very localized (after its introduction) (?)
    REAL(r_std), INTENT(inout)                             :: when_growthinit_loc    !! how many days ago was the beginning of 
                                                                                     !! the growing season (days)
    REAL(r_std), INTENT(inout)                             :: npp_longterm_loc       !! "long term" net primary productivity
                                                                                     !! @tex ($gC m^{-2} year^{-1}$) @endtex
    REAL(r_std), INTENT(inout)                             :: age_loc                !! mean age (years)
    REAL(r_std), INTENT(inout)                             :: lm_lastyearmax_loc     !! last year's maximum leaf mass for each PFT 
                                                                                     !! @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(ncirc), INTENT(inout)           :: circ_class_n_loc       !! Number of individuals in each circ class
                                                                                     !! @tex $(ind m^{-2})$ @endtex
    REAL(r_std), DIMENSION(ncirc,nparts,nelements), INTENT(inout) :: circ_class_biomass_loc !! Biomass components of the model tree  
                                                                                     !! within a circumference class
                                                                                     !! class @tex $(g C ind^{-1})$ @endtex
    REAL(r_std), DIMENSION(nelements), INTENT(inout)       :: atm_to_bm_loc          !! CO2 or N taken from the atmosphere to get
                                                                                     !! the seed C and N
                                                                                     !! to create the seedlings 
                                                                                     !!  @tex (gC.m^{-2}dt^{-1})$ @endtex  
    REAL(r_std), DIMENSION(nleafages), INTENT(inout)       :: leaf_frac_loc          !! fraction of leaves in leaf age 
                                                                                     !! class (unitless;0-1)
    REAL(r_std), INTENT(inout)                             :: KF_loc                 !! Scaling factor to convert sapwood mass
                                                                                     !! into leaf mass (m) - this variable is 
                                                                                     !! passed to other routines
    REAL(r_std), INTENT(inout)                             :: k_latosa_adapt_loc     !! Leaf to sapwood area adapted for long 
                                                                                     !! term water stress (m)
    REAL(r_std), DIMENSION(ncirc,nparts,nelements), INTENT(inout)  :: bm_sapl_2D_loc !! Sapling biomass for the functional
                                                                                     !! allocation with a dimension for the 
                                                                                     !! circumference classes
                                                                                     !! @tex $(gC.ind^{-1})$ @endtex

  !! 0.4 Local variables
    REAL(r_std)                                       :: c0_alloc           !! Root to sapwood tradeoff parameter
    INTEGER(i_std)                                    :: ipar, ilev         !! index (unitless)
    INTEGER(i_std)                                    :: icir, iele         !! index (unitless)
    INTEGER(i_std)                                    :: imbc, ipart        !! index (unitless)
    INTEGER(i_std)                                    :: deb,fin, imaxt     !! index (unitless)
    REAL(r_std)                                       :: LF                 !! Scaling factor to convert sapwood mass
                                                                            !! into root mass (unitless)
    REAL(r_std)                                       :: lstress_fac        !! Light stress factor, based on total
                                                                            !! transmitted light (unitless, 0-1)
    REAL(r_std), DIMENSION(ncirc)                     :: circ_class         !! circumference of individual trees
    REAL(r_std), DIMENSION(nparts)                    :: bm_init            !! Biomass needed to initiate the next 
                                                                            !! planting @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(ncirc)                     :: nb_trees_i         !! Number of trees in each twentith 
                                                                            !! circumference quantile of the 
                                                                            !! distribution (ind) 
    REAL(r_std)                                       :: excedent           !! Number of trees after truncation to be 
                                                                            !! reallocated to smaller quantiles of the 
                                                                            !! distribution (ind) 
    REAL(r_std)                                       :: max_circ_init      !! Maximum initial circumferences of the 
                                                                            !! truncated exponential distribution (cm) 
    REAL(r_std)                                       :: dia_init_recr      !! Initial diameter of a recruit (m). Note that 
                                                                            !! this value may differ from the input argument 
                                                                            !! dia_int which is the diameter distribution of 
                                                                            !! a newly planted PFT
    REAL(r_std), DIMENSION(nmbcomp,nelements)         :: check_intern       !! Contains the components of the internal
                                                                            !! mass balance chech for this routine
                                                                            !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(nelements)                 :: closure_intern     !! Check closure of internal mass balance
                                                                            !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(nelements)                 :: pool_start         !! Start and end pool of this routine 
                                                                            !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(nelements)                 :: pool_end           !! Start and end pool of this routine 
                                                                            !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std)                                       :: sla_est            !! A first estimate of sla in case its calculation is 
                                                                            !! dynamic @tex $(m^2.gC^{-1})$ @endtex
    REAL(r_std), DIMENSION(ncirc)                     :: bypass_lm          !! leaf mass before it is set to zero
                                                                            !! for deciduous PFTs. This value is
                                                                            !! in the calculation of KF (gC tree-1)
    REAL(r_std)                                       :: veget_max_begin    !! temporary storage of veget_max_loc to check area conservation    
    REAL(r_std)                                       :: ind_init           !! Maximum number of individuals when establishing a new forest (#/m^{-2})     
    LOGICAL                                           :: establish          !! Logical flag to establish a new young vegetation
                                                                            !! after a stand replacing disturbance 
    LOGICAL                                           :: recruite           !! Logical flag to add saplings of a PFT under an
                                                                            !! existing canopy of the same PFT 
    REAL(r_std), DIMENSION(ncirc)                     :: old_circ_class_n   !! Old number of individuals in each circ class
                                                                            !! @tex $(ind m^{-2})$ @endtex
    REAL(r_std), DIMENSION(ncirc,nparts,nelements)    :: old_circ_class_biomass !! Old biomass components of the model tree  
                                                                            !! within a circumference class
                                                                            !! class @tex $(g C ind^{-1})$ @endtex
    REAL(r_std), DIMENSION(nparts,nelements)          :: old_biomass        !! Old stand level biomass 
                                                                            !! tex $(gC.m^{-2})$ @endtex
    REAL(r_std)                                       :: old_ind            !! density of old individuals                                        
    REAL(r_std)                                       :: cn_leaf,cn_root    !! CN ratio of leaves, and roots (gC/gN)
    REAL(r_std)                                       :: cn_wood            !! CN ratio of wood pool (gC/gN)
    REAL(r_std)                                       :: tmp_ind            !! temporary number of individuals per m-2
    REAL(r_std), DIMENSION(ncirc)                     :: circ_class_n_new   !! variable used to sort circ_class_n
    REAL(r_std), DIMENSION(ncirc,nparts,nelements)    :: circ_class_biomass_new !! variable used to sort circ_class_biomass
    REAL(r_std)                                       :: k_latosa_tmp       !! Temporaray variable in the calculation of 
                                                                            !! k_latosa_adapt
    REAL(r_std)                                       :: n_init             !! Initial number of individuals calculated (ind m^(-2))
    REAL(r_std)                                       :: tmp                !! Temporary variable to prepare information for xios
    REAL(r_std)                                       :: qmd, mod_qmd, max_x !! use to estimate the diameter of the tree circuference class before using make samplin 
    REAL(r_std), DIMENSION(ncirc)                     :: mod_qmd_arr        !! Multiplicateurs de diametre par classe, NORMALISES pour que le moment d ordre 2 vaille 1
    REAL(r_std)                                       :: qmd_norm           !! Facteur de normalisation du moment d ordre 2
    REAL(r_std)                                       :: weibull_resolution !! idem
    REAL(r_std), DIMENSION(ncirc)                     :: exp_circ_class_n   !! Temporary variable for expected circ class n
                                                                            !! @tex $(ind m^{-2})$ @endtex
    REAL(r_std)                                       :: real_qmd           !! Quick fix for mismatch between qmd_init and qmd from new sapplings
    INTEGER(i_std)                                    :: iage_pr            !! Age class index of this PFT slot (1..nagec), establishment RDI
    REAL(r_std)                                       :: rdi_use            !! RDI used at establishment (rdi_start, or rdi_max for cold-start desync)
    REAL(r_std), DIMENSION(ncirc)                     :: circ_class_n_grass !! temporary grass density in case of lcchange
                                                                            !! @tex $(ind m^{-2})$ @endtex
    REAL(r_std)                                       :: Cother             !! The sum of carbon from all compartments except for fruit, 
                                                                            !! reserve and labile pool (gC ind-1)
    REAL(r_std)                                       :: Trl0               !! The sum of individual reserve and labile carbon before 
                                                                            !! density adjustment (gC ind-1)
    REAL(r_std)                                       :: Nume               !! Numerator to calculate grass density
    REAL(r_std)                                       :: Deno               !! Denominator to calculate grass density
    REAL(r_std)                                       :: ex_bio             !! Sum of new reserve and labile carbon (gC ind-1)
    REAL(r_std)                                       :: ratio_res          !! Ratio of reserve carbon to the sum of reserve and labile carbon (Unitless)
    REAL(r_std)                                       :: ratio_lab          !! Ratio of labile carbon to the sum of reserve and labile carbon (Unitless)
    REAL(r_std)                                       :: ratio_den          !! Ratio of old density to new density (Unitless) 
    REAL(r_std)                                       :: lai_target         !! Target LAI @tex $(m^{2}m^{-2})$ @endtex
    REAL(r_std), DIMENSION(ncirc)                     :: Cl_target          !! Individual plant maximum leaf mass given
                                                                            !! its current sapwood mass
                                                                            !! @tex $(gC.plant^{-1})$ @endtex
    REAL(r_std)                                       :: gtemp              !! Turnover coefficient of labile C pool
    REAL(r_std)                                       :: tl                 !! Long term annual mean temperature (C)
    REAL(r_std), DIMENSION(nparts,nelements)          :: tmp_bm             !! temporary variable to indicate biomass for each PFT
                                                                            !! over per unit PFT area.
                                                                            !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std),PARAMETER                             :: scale_inc = 1      !! Scaling factor to calculate reserve and labile carbon 
                                                                            !! when increasing density (unitless)
                                                                            !! Has been used to tune the results. Could be externalized 
                                                                            !! and made PFT-specific when found to be impactful.
    REAL(r_std),PARAMETER                             :: target_c4 = 1      !! Scaling factor for reserve and labile carbon target 
                                                                            !! only in C4 grassland (unitless)
                                                                            !! Has been used to tune the results. Could be externalized 
                                                                            !! and made PFT-specific when found to be impactful.
    REAL(r_std)                                       :: temp_dens_target   !! Temporary and local target density.
!_ ================================================================================================================================

    IF (firstcall_prescribe) THEN
      !! Initialize local printlev
      ! Guillaume M. -- get_printlev calls getin_p, a COLLECTIVE bcast. Called here,
      ! inside the compute loop under a local firstcall, it hangs the run as soon as the
      ! ranks stop reaching it together. Use printlev, read once at startup by every rank.
      printlev_loc=printlev
      firstcall_prescribe=.FALSE.
    END IF

    ! Prescribe is called at the pixel level so set a higher printlev
    IF (printlev_loc>=3) WRITE(numout,*) 'Entering stomate_prescribe'

    ! Guillaume M. -- Cold-start age-class desync: default the establishment flag.
    IF (PRESENT(established_out)) established_out = .FALSE.

    !! 1. Initialize biomass at first call
    ! Write circ_class_n_grass in beginning to get the value of the
    ! current grass density, prepared for lcchange with dynamic grass
    ! density
    ! In the case of grasslands and with flag on, veget_max
    ! for grassland < min_stomate, circ_class_n will
    ! be set to zero to avoid unexpected condition.
    IF ( ( .NOT. is_tree(ivm) ) .AND. &
       ( natural(ivm) ) .AND. &
       ok_dyn_grass_den) THEN

       IF (veget_max_loc .LT. min_stomate) THEN
          circ_class_n_loc(:) = zero
       END IF

       circ_class_n_grass=circ_class_n_loc
       
       ! Debug
       IF (printlev_loc .GE. 4) THEN
          WRITE(numout,*) 'ivm, reserve carbon total', &
               ivm, circ_class_biomass_loc(1,icarbres,icarbon)*circ_class_n_loc(1)
          WRITE(numout,*) 'ivm, labile carbon total', &
               ivm, circ_class_biomass_loc(1,ilabile,icarbon)*circ_class_n_loc(1)
       END IF
       !-
       
    END IF

    !! Calculate individual reserve and labile target C only for grasslands
    ! Initialize
    tot_res_target_loc(:) = zero
    tot_lab_target_loc(:) = zero
    lstress_fac = zero
    LF = 1._r_std

    IF( ( .NOT. is_tree(ivm) ) .AND. &
       ( natural(ivm) ) .AND. &
       (ok_dyn_grass_den) .AND. &
       ( veget_max_loc .GT. min_stomate ) ) THEN
       !! Calculate gtemp and lai_target
       !! gtemp, referring to the stomate_growth_fun_all.f90
       tl = t2m - ZeroCelsius
       tmp_bm(:,:) = cc_to_biomass(1,ivm,&
                        circ_class_biomass_loc(:,:,:), circ_class_n_loc(:))
       IF (tl .GT. tmin_labile(ivm)) THEN
         gtemp = EXP((e0_labile(ivm))*(1.0/(tref_labile(ivm)-tmin_labile(ivm)) -&
            1.0/(tl-tmin_labile(ivm))))
       ELSE
         ! Too cold to grow
         gtemp = zero
       ENDIF

       ! If there is a plant, and we are either at the very start or in
       ! the growing season not during senescences, calculate labile
       ! pool use for growth.
       ! CYmark: we allow calculation of bm_alloc_tot as well for
       ! ipresenescence stage.
       IF ( SUM(circ_class_n_loc(:)) .GT. min_stomate .AND. &
          ( plant_status_loc .EQ. ibudbreak .OR. &
          plant_status_loc .EQ. icanopy .OR. &
          plant_status_loc .EQ. ipresenescence ) .AND. &
          SUM(circ_class_biomass_loc(:,ileaf,icarbon)) .GT. min_stomate) THEN

          IF ( (tmp_bm(ilabile,icarbon) .GT. min_stomate) .OR. &
               (tmp_bm(icarbres,icarbon) .GT. min_stomate) ) THEN
             ! Truncate gtemp between zero and 1. If we set the upper
             ! bound to one, we may run into numerical (precision)
             ! problems
             ! later caused by very small (10e-15) negative values. Rather
             ! than dealing with the precision issues it is easier to use
             ! 0.99 instead. 0.99 may be too high so this parameter was
             ! externalized and is pft-dependent.
             gtemp = MAX(MIN(gtemp, un-always_labile(ivm)),zero)
          ELSE
             ! There is nothing to allocate so we could as well set
             ! gtemp to zero
            gtemp = zero
          ENDIF

          ! Prioritize the use of the carbohydrate pool. Move
          ! carbohydrates to the labile pool.
          ! CYmark: Such moving reserve  to labile pool is not allowed for
          ! ipresenescence stage.
          IF ( ( plant_status_loc .EQ. ibudbreak .OR. &
                 plant_status_loc .EQ. icanopy ) .AND. &
               (tmp_bm(icarbres,icarbon) .GT. min_stomate) ) THEN
             tmp_bm(ilabile,icarbon) = &
                  tmp_bm(ilabile,icarbon) + 0.05 * &
                  tmp_bm(icarbres,icarbon)
             tmp_bm(icarbres,icarbon) = &
                  tmp_bm(icarbres,icarbon) * 0.95
          ENDIF

          ! This would be a good place to update the plant level
          ! labile and carbres pool but as it may be subject to
          ! further modifications in the next block of code, it
          ! will be done later

       ELSE
          ! The plant is absent, senescent or dead so there will be
          ! no allocation. Set gtemp to zero to keep the output files
          ! clean. Due to this line, gtemp will be zero as soon as
          ! plant_status becomes isenescent or idormant.
          gtemp = zero

       ENDIF
       
       !! lai_target, referring to the stomate_growth_fun_all.f90
       ! Calculate KF
       k_latosa_tmp = k_latosa_adapt_loc + lstress_fac * &
              (k_latosa_max(ivm)-k_latosa_min(ivm))

       IF (sla_dyn) THEN
          sla_est = &
                biomass_to_lai(circ_class_biomass_loc(1,ileaf,icarbon),ivm)/&
                circ_class_biomass_loc(1,ileaf,icarbon)
       ELSE
          sla_est = sla(ivm)
       END IF

       KF_loc = k_latosa_tmp / &
                (sla_est * tree_ff(ivm) * pipe_density(ivm))
       
       ! Calculate leaf to root area
       ! calculate c0_alloc here 
       c0_alloc = calculate_c0_alloc(ivm, &
            longevity_eff_root_loc,longevity_eff_sap_loc)
       LF = c0_alloc * KF_loc

       ! Debug
       IF (printlev_loc .GE. 4) THEN
          WRITE(numout,*) 'DebugMode,LF,c0_alloc',LF,c0_alloc
       END IF
       !-
       
       ! Calculate Cl_target
       Cl_target = MAX( circ_class_biomass_loc(1,isapabove,icarbon) * KF_loc,&
                      circ_class_biomass_loc(1,iroot,icarbon) * LF, &
                      circ_class_biomass_loc(1,ileaf,icarbon) )
       ! lai_target
       lai_target = cc_to_lai(Cl_target,circ_class_n_loc(:), ivm)

       ! Calculate areal level target for reserve
       tot_res_target_loc(icarbon) = &
                     MIN(deciduous_reserve(ivm) * (tmp_bm(iroot,icarbon) + &
                     tmp_bm(isapabove,icarbon) + &
                     tmp_bm(isapbelow,icarbon)), &
                     lai_to_biomass(lai_target,ivm) * &
                     (1.+root_reserve(ivm)/LF))

       ! Calculate areal level target for labile
       tot_lab_target_loc(icarbon) = gtemp*labile_reserve(ivm)*(tmp_bm(ileaf,initrogen)+ &
              tmp_bm(iroot,initrogen) + tmp_bm(ifruit,initrogen) + &
              tmp_bm(isapabove,initrogen) + tmp_bm(isapbelow,initrogen))
       tot_lab_target_loc(icarbon) = MAX(tot_lab_target_loc(icarbon), 10. * gpp_week )

       ! Tune reserve and labile target only for C4 grasslands
       IF (is_c4(ivm) .AND. &
            natural(ivm)) THEN
          tot_res_target_loc(icarbon)= target_c4 * tot_res_target_loc(icarbon)
          tot_lab_target_loc(icarbon) =target_c4*tot_lab_target_loc(icarbon)
       END IF

       ! Debug
       IF (printlev_loc .GE. 4) THEN
          WRITE(numout,*) 'ivm', ivm
          WRITE(numout,*) 'Density', circ_class_n_loc(1)
          WRITE(numout,*) 'KF_loc, LF ', KF_loc, LF
          WRITE(numout,*) 'gtemp', gtemp
          WRITE(numout,*) 'gpp_week', gpp_week
          WRITE(numout,*) 'Cl_target, lai_target ', Cl_target, lai_target
          WRITE(numout,*) 'tmp_bm', tmp_bm(:,icarbon)
       END IF
       !-
       
    END IF

    !! 1.0 Initialize variables that won't be calculated for all PFTs,
    !!     to avoid a NaN if the value for bare soil is ever written.
    c0_alloc=zero
    !! 1.1 Store for closing the mass balance with recruitment
    old_biomass(:,:) = cc_to_biomass(1,ivm,&
         circ_class_biomass_loc(:,:,:), & 
         circ_class_n_loc(:))
    new_ind_loc = zero

    !! 1.3 Initialize check for mass balance closure 
    IF (err_act.GT.1) THEN

       check_intern(:,:) = zero
       pool_start(:) = zero
       DO iele = 1,nelements

          ! atm_to_bm has as intent inout, the variable 
          ! accumulates carbon over the course of a day. 
          ! Use the difference between start and the end of 
          ! this routine
          check_intern(iatm2land,iele) = - un * &
               atm_to_bm_loc(iele) * veget_max_loc * dt

          DO ipar = 1,nparts
          
             DO icir = 1,ncirc

                ! Initial biomass pool. Use circ_class_biomass
                ! because that variable is the prognostic variable
                ! biomass is syncronized to circ_class_biomass. If 
                ! many diameter classes are used, each diameter class
                ! separatly passes all criteria but when combined
                ! a mass balance problem occurs.
                pool_start(iele) = pool_start(iele) + &
                   circ_class_biomass_loc(icir,ipar,iele) * &
                   circ_class_n_loc(icir) * veget_max_loc

             ENDDO
             
          ENDDO

       ENDDO

       !! 1.3 Initialize check for area conservation
       veget_max_begin = veget_max_loc

    ENDIF ! err_act.GT.1

  !! 2. Prescribe carbon and nitrogen biomass

    ! We want the vegetation to regrow after a stand replacing disturbance
    ! If there is supposed to be vegetation here (according to veget_max_loc)
    ! but there isn't (according to ind), then establish a new young
    ! vegetation in that PFT.
    
    IF (is_tree(ivm)) THEN
       ! Select the target density as a function of forest
       ! management
       temp_dens_target = dens_target(ivm,forest_managed_loc)
    ELSE
       ! Prescribe a value that enables passing the loops
       ! below. For grasses and croplands, recruitment is
       ! not used. 
       temp_dens_target = 1.5
    END IF
    
    ! Debug
    IF(printlev_loc.GE.4 .AND. write_debug)THEN
       WRITE(numout,*) 'plant_status ', plant_status_loc
       WRITE(numout,*) 'veget_max ', veget_max_loc
       WRITE(numout,*) 'recruitment_pft ',recruitment_pft
       WRITE(numout,*) 'ok_dgvm ',  ok_dgvm
       WRITE(numout,*) 'biomass - C, ', &
            SUM(SUM(circ_class_biomass_loc(:,:,icarbon),1))
       WRITE(numout,*) 'biomass - N, ', &
            SUM(SUM(circ_class_biomass_loc(:,:,initrogen),1))
       WRITE(numout,*) 'ts_annual_proc ', ts_annual_proc
       WRITE(numout,*) 'circ_class_n, ', SUM(circ_class_n_loc(:))
       WRITE(numout,*) 'dens_target, ', temp_dens_target*ha_to_m2
       WRITE(numout,*) 'do_now_recruit, ',do_now_recruit
    ENDIF
    !
    
    ! Debug check for FM and PFT values in dens_target
    IF (printlev_loc.GE.4 ) THEN
       WRITE(numout,*) 'ivm and forest_managed_loc used in dens_target : ivm = ', ivm
       WRITE(numout,*) 'forest_managed_loc = ', forest_managed_loc
       WRITE(numout,*) 'dens_target = ', temp_dens_target
       CALL flush(numout) 
    ENDIF

    ! Decide whether we have to establish new vegetation (::establish)
    ! allow recruitment (::recruite). These two local variables were
    ! introduced to allow maximal flexibility so within a single
    ! run, PFTs with and without recruitment can co-exist.
    IF( ( (plant_status_loc .EQ. iprescribe) .OR. & 
         (plant_status_loc .EQ. idead) ) .AND. &
         (SUM(SUM(circ_class_biomass_loc(:,:,icarbon),1)) .LE. min_stomate) .AND. &
         (veget_max_loc .GT. min_stomate) .AND. &
         (.NOT. ok_dgvm) ) THEN

       ! The vegetation just died from replacing disturbances
       ! (two first conditions). The new vegetation could be 
       ! of the same or a different PFT than the previous 
       ! vegetation. 
       establish = .TRUE.
       recruite = .FALSE.
       ! Debug
       IF (printlev_loc>=4) THEN
          WRITE(numout,*) 'We are calling stomate_prescribe.f90 &
               & to initialize PFT  ',ivm
       ENDIF
       !-

    ELSEIF( ( (plant_status_loc .NE. iprescribe) .OR. &
         (plant_status_loc .NE. idead) ) .AND. &
         (SUM(SUM(circ_class_biomass_loc(:,:,icarbon),1)) .GT. min_stomate) .AND. &
         (veget_max_loc .GT. min_stomate) .AND. &
         ((forest_managed_loc .EQ. ifm_none) .OR. (forest_managed_loc .EQ. ifm_uneven)) .AND. &
         (.NOT. ok_dgvm) .AND. &
         ts_annual_proc .AND. &
         (SUM(circ_class_n_loc(:)) .GT. temp_dens_target*ha_to_m2) .AND. &
         do_now_recruit) THEN

       ! There is already vegetation (first three conditions)
       ! and this specific PFT has recruitment (condition 5)
       ! Only use recruitment when the recruitment flag is TRUE
       ! (condition 7)
       ! no recruit when the stand density is below dens_target because 
       ! the model will slowly  but surely kill the forest before a new
       ! establishment
       establish = .FALSE.
       recruite = .TRUE. 
       
       ! Debug
       IF (printlev_loc>=4) THEN
          WRITE(numout,*) 'We are calling stomate_prescribe.f90 &
               & for sapling recruitment in PFT  ',ivm
          WRITE(numout,*) 'Density of stems', SUM(circ_class_n_loc(:))
       ENDIF
       !-
   
    ELSEIF( ok_dyn_grass_den .AND. &
         (plant_status_loc .EQ. ibudbreak .OR. &
         plant_status_loc .EQ. icanopy ) .AND. &
         .NOT. is_tree(ivm) .AND. &
         natural(ivm) .AND. &
         veget_max_loc .GT. min_stomate .AND. &
         .NOT. ok_dgvm .AND. &
         (circ_class_n_loc(1)*(circ_class_biomass_loc(1,icarbres,icarbon) + &
         circ_class_biomass_loc(1,ilabile,icarbon)) .GT. &
         tot_res_target_loc(icarbon)+tot_lab_target_loc(icarbon)) .AND. &
         (SUM(circ_class_biomass_loc(1,:,icarbon)) .GT. min_stomate) .AND. &
         (circ_class_n_loc(1) .LT. 1) .AND. &
         (tot_res_target_loc(icarbon)+tot_lab_target_loc(icarbon) .GT. min_stomate)) THEN

       ! For grasslands to increase density
       ! The conditions include:
       ! (1) the flag ok_dyn_grass_den for the calculation of dynamic grass
       ! density is switched on; (2) during the growing season (ibudbreak, icanopy); 
       ! (3) This PFT belongs to grasslands category (.NOT. is_tree(ivm) .AND. natural(ivm)); 
       ! (4) The value of Vfra for each grassland is greater than 10-8; 
       ! (5) Not in the ok_dgvm (Dynamic Global Vegetation Model) mode; 
       ! (6) Current reserve and labile carbon are greater than the target of reserve and
       ! labile carbon; 
       ! (7) The sum of carbon biomass is greater than 10-8; 
       ! (8) The current density is lower than 1;
       ! (9) the target of reserve and labile carbon is greater than 10-8.
       
       establish = .FALSE.
       recruite = .TRUE.
       
    ELSEIF( (veget_max_loc .LT. min_stomate) .AND. &
         (SUM(SUM(circ_class_biomass_loc(:,:,icarbon),1)) .LT. min_stomate) .AND. &
         (SUM(circ_class_n_loc(:)) .LT. min_stomate) )THEN
       
       ! The PFT is not defined. Initialize the variables
       ! to ensure that all goes well when moving and merging 
       ! different age classes.
       circ_class_n_loc(:) = zero
       circ_class_biomass_loc(:,:,icarbon) = zero 

       ! Conditions are not fullfilled for establishing
       ! new vegetation
       establish = .FALSE.
       recruite = .FALSE.
       
    ELSEIF (veget_max_loc .GT. min_stomate) THEN

       ! Conditions are not fullfilled for either establishing 
       ! a new young vegetation or having recruitment under an
       ! exisisting vegetation of the same PFT.
       establish = .FALSE.
       recruite = .FALSE.
       
       ! Debug
       IF (printlev_loc>=4) THEN
          WRITE(numout,*) 'We are calling stomate_prescribe.f90 &
               & but nothing should be done for PFT ',ivm
       ENDIF
       !-
             
    ELSE

       ! Unexpected condition
       WRITE(numout,*) 'plant_status ', plant_status_loc
       WRITE(numout,*) 'veget_max ', veget_max_loc
       WRITE(numout,*) 'recruitment_pft ',recruitment_pft
       WRITE(numout,*) 'ok_dgvm ',  ok_dgvm
       WRITE(numout,*) 'biomass - C, ', &
            SUM(SUM(circ_class_biomass_loc(:,:,icarbon),1))
       WRITE(numout,*) 'biomass - N, ', &
            SUM(SUM(circ_class_biomass_loc(:,:,initrogen),1))
       WRITE(numout,*) 'ts_annual_proc ', ts_annual_proc
       WRITE(numout,*) 'circ_class_n, ', SUM(circ_class_n_loc(:))
       WRITE(numout,*) 'dens_target, ', temp_dens_target*ha_to_m2
       WRITE(numout,*) 'do_now_recruit, ',do_now_recruit
       CALL flush(numout)
       
       CALL ipslerr_p (3,'stomate_prescribe.f90',&
            'Unexpected condition - time to rethink the', &
            'IF that choses between establish and recruite','')
      
    ENDIF

    !! 2.2 Calculate variables that will determine the size of the saplings
    IF( establish .OR. recruite )THEN

       !! 2.2.1 Calculate c0_alloc
       ! The calculation starts with c0_alloc, it only depends on PFT and 
       ! the effective longiveties.
       c0_alloc = calculate_c0_alloc(ivm, &
            longevity_eff_root_loc,longevity_eff_sap_loc)
             
       !! 2.2.2 Calculate k_latosa for newly established PFTs
       !  If a new vegetation has to be established we recalculate
       !  k_latosa from the parameter file (see note on physiological
       !  memory of k_latosa below). If we add recruits we use the KF of
       !  the existing vegetation.
       IF (establish) THEN

          ! Initialize tree level variables
          circ_class_n_loc(:) = zero
          circ_class_biomass_loc(:,:,:) = zero
          
          ! Assume that there is no physiological memory for k_latosa
          ! between different generations (this is probably true when
          ! trees are planted but seems less likely for natural
          ! regeneration. The first test that implemented memory has
          ! caused problems with the LAI from the second generation
          ! onwards. Nevertheless, we still use k_latosa_adapt. In case
          ! someone want to try again in the future, the structure of
          ! the code is ready to account for memory
          k_latosa_adapt_loc = k_latosa_min(ivm)

          ! Debugging
          IF (printlev_loc>=4) THEN
             WRITE(numout,*) 'Establish a PFT (stomate_prescribe.f90):'
             WRITE(numout,*) ' Imposing initial biomass for prescribed &
                  & trees, initial reserve mass for prescribed grasses.'
             WRITE(numout,*) ' Declaring prescribed PFTs present.'
             WRITE(numout,*) ' Establish ',establish,' Recruite ',recruite 
          ENDIF
          !-
                
       ENDIF ! establish

       !! 2.2.3 Calculate Lightstress
       !  Note that light and water stress have an opposite effect on KF. 
       !  Waterstress will decrease KF because more C should be allocated 
       !  to the roots. Light stress will increase KF because more C should 
       !  be allocated to the leaves.
       !  Lightstress varies from 0 to 1 and is calculated from the canopy 
       !  structure (veget).
       IF (establish) THEN

          ! Given that there is no vegetation at this point, the
          ! lightstress cannot be calculated and is set to 1.
          ! However, as soon as we grow a canopy it will experience
          ! light stress and so KF should be adjusted.
          lstress_fac = un

       ELSEIF (recruite) THEN 

          ! We already have vegetation, so we should do recruitment if 
          ! we want it for this PFT. First lets compute light stress based 
          ! on the fraction of light reaching the ground (level 1),
          ! averaged over all daytime time steps (sinang>0) of the vegetative 
          ! period (gpp>0). This is based on Pgap and accounts for canopy
          ! structure. 
          lstress_fac = light_tran_to_floor_season_loc

          ! This value is accumulated in stomate.f90. It is accumulated
          ! over an entire year and used only the last day of the year.
          ! After it has been used, it can be reset. Note that 
          ! age_class_dist takes place before recruitment so it has be
          ! moved from one PFT to another first.
          light_tran_to_floor_season_loc = zero

          ! Debug
          IF(printlev_loc>=4)THEN
             WRITE(numout,*) 'lstress_fac in stomate_prescribe is ', &
                  lstress_fac*100, ' (%)'
          ENDIF
          !-
          
       ENDIF ! establish/recruite for light

    ENDIF ! establish .OR. recruite for initialization

          
    !! 2.3 Initial vegetation has to be prescribed only when the vegetation is 
    ! static
    IF (establish .OR. recruite) THEN

       !! 2.3.1 Initialize woody PFT's 
       !  Use veget_max_loc to check whether the PFT is present. If the PFT 
       !  is present but it has no biomass, prescribe its biomass (gC m-{2}).
       IF ( is_tree(ivm) ) THEN

          ! prescribe/recruitment branching for actual process
          IF (establish) THEN
             
             ! Allocation factors
             ! Sapwood to root ratio
             ! Following Magnani et al. 2000 "In order to decreases hydraulic 
             ! resistance, the investment of carbon in fine roots or sapwood 
             ! yields to the plant very different returns, both because of 
             ! different hydraulic conductivities and because of the strong 
             ! impact of plant height on shoot resistance. On the other hand, 
             ! fine roots and sapwood have markedly different longevities and 
             ! the cost of production, discounted for turnover, will differ 
             ! accordingly. Optimal growth under hydraulic constraints requires 
             ! that the ratio of marginal hydraulic returns to marginal annual 
             ! cost for carbon investment in either roots or sapwood be the 
             ! same (Bloom, Chapin & Mooney 1985; Case & Fair 1989). This
             ! is formalized in equation (13) and further derived to obtain 
             ! equation (17) in Magnani et al 2000. The latter is implemented 
             ! here.Pipe_density is given in gC/m-3, convert to kg/m-3. And 
             ! apply equation (17) in Magnani et al 2000. Note that c0_alloc 
             ! was calculated at the start of this routine. The calculation 
             ! itself is done in function_library
             
             ! Calculate leaf area to sapwood area
             ! To be consistent with the hydraulic limitations and pipe theory,
             ! k_latosa is calculated from equation (18) in Magnani et al.
             ! To do so, total hydraulic resistance and tree height need to 
             ! known. This poses a problem as the resistance depends on the leaf 
             ! area and the leaf area on the resistance. There is no independent 
             ! equation and equations 12 and 18 depend on each other and 
             ! substitution would be circular. Hence prescribed k_latosa values 
             ! were obtained from observational records and are given in 
             ! mtc_parameters.f90. 
             
             ! The relationship between height and k_latosa as reported in 
             ! McDowell et al 2002 and Novick et al 2009 could be implemented 
             ! to adjust k_latosa for the height of the stand. This did NOT 
             ! result in a realistic model behavior 
             ! k_latosa(ipts,ivm) = wstress_fac(ipts,ivm) * &
             !    (k_latosa_max(ivm) - (latosa_height(ivm) * &
             !    (SUM( nb_trees_i(:) *ave_tree_height(:) ) / &
             !       SUM( nb_trees_i(:) ))))
             
             ! Also k_latosa has been reported to be a function of CO2 
             ! concentration (Atwell et al. 2003, Tree Physiology, 23:13-21 
             ! and Pakati et al. 2000, Global Change Biology, 6:889-897). 
             ! This effect is not accounted for in the current code
             
             ! Finally, k_latosa is also reported to be a function of 
             ! diameter (i.e. stand thinning, Simonin et al 2006, Tree 
             ! Physiology, 26:493-503). Here the relationship with thinning was 
             ! interpreted as a realtionship with light stress. Note that light 
             ! stress cannot be calculated at this time in the model because 
             ! there is no canopy (that's why we are in prescribe!) and
             ! so there is no lightstress. lstress was therefore set to one 
             ! (see above). We prefered to keep this redundancy in the code 
             ! because it makes it clear that k_latosa is calculated in the 
             ! same way in prescribe.f90 and growth_fun_all.f90
             ! k_latosa(ipts,ivm) = (k_latosa_adapt(ipts,ivm) + &
             !     (lstress_fac(ipts,ivm) * &
             !     (k_latosa_max(ivm)-k_latosa_min(ivm))))
             
             ! +++CHECK+++
             ! Ideally waterstress and lightstress should both
             ! affect k_latosa. One of the following approaches could
             ! be used as a starting point.
             !  k_latosa(ipts,ivm) = k_latosa_min(ivm) + &
             !     (wstress_fac(ipts,ivm) * lstress_fac(ipts,ivm) * &
             !     (k_latosa_max(ivm)-k_latosa_min(ivm)))
             !  k_latosa(ipts,ivm) = wstress_fac(ipts,ivm) * (k_latosa_min(ivm) + &
             !      (lstress_fac(ipts,ivm) * &
             !      (k_latosa_max(ivm)-k_latosa_min(ivm))))
             ! ++++++++++++
             
             ! Calculate conversion coefficient for sapwood area to leaf 
             ! area 
             ! (1) The scaling parameter between leaf and sapwood mass 
             ! is derived from LA_ind = k_latosa * SA_ind, where 
             ! LA_ind = leaf area of an individual, SA_ind is the sapwood 
             ! area of an individual and k_latosa a pipe-model parameter
             ! (2) LA_ind = Cl * sla
             ! (3) Cs = SA_ind * height * wooddensity * tree_ff
             ! Substitute (2) and (3) in (1)
             ! Cl = Cs * k1 / (wooddensity * sla * tree_ff * height)
             ! Cl = Cs*KF/height, where KF is in (m)
             ! KF is passed to the allocation routine and it is saved in the 
             ! restart file.
             !+++CHECK+++
             ! At one point it looked like a good idea to take the max of two
             ! options but by doing so we cannot recalculate KF in phenology.
             ! It has not been confirmed that this is really a problem. Use the
             ! simpelest approach but leave the alternative in the code as a
             ! suggestion of a possible solution in case something goes wrong.
!!$             k_latosa_tmp = MAX(k_latosa_min(ivm),k_latosa_adapt(ipts,ivm) + &
!!$                  (lstress_fac(ipts,ivm) * &
!!$                  (k_latosa_max(ivm)-k_latosa_min(ivm))))
             !+++++++++++
             k_latosa_tmp = k_latosa_adapt_loc + (lstress_fac * &
                  (k_latosa_max(ivm)-k_latosa_min(ivm)))
             
             IF (sla_dyn) THEN 
                KF_loc = k_latosa_tmp / &
                     (slainit(ivm) * pipe_density(ivm) * tree_ff(ivm))
             ELSE
                KF_loc = k_latosa_tmp / &
                     (sla(ivm) * pipe_density(ivm) * tree_ff(ivm))
             END IF


             ! quadratic diameter is now a pft parameter and it is use the guess
             ! the weibull circumference class distribution.
             qmd = qmd_init(ivm)

             ! Guillaume M. -- Cold-start age-class desync: when the caller supplies a
             ! per-(point,age-class) initial diameter (OK_DIA_STAGGER), use it instead of
             ! the uniform qmd_init(ivm). A value >= val_exp/2, or an absent argument,
             ! means "no override" and keeps trunk behaviour bit-for-bit.
             IF (PRESENT(qmd_init_override)) THEN
                IF (qmd_init_override < val_exp/2._r_std) qmd = qmd_init_override
             ENDIF
             IF (PRESENT(established_out)) established_out = .TRUE.

             ! We know the quadratic mean diameter qmd but the use make_sapling
             ! we have to estimate the diameter of each circufenrence class.
             ! weibull_resolution is based on the definition and calculation 
             ! made in the weibull_dist function in order to be consistante. 
             exp_circ_class_n(:) = weibull_class_dist(qmd * m_to_cm, ivm, &
                  forest_managed_loc)

             ! Guillaume M. -- Single source with weibull_class_dist: this block used to
             ! carry its own copy of the same SELECT CASE. The class multipliers below and
             ! the class WEIGHTS come from the same grid, so two copies drifting apart
             ! would build a stand whose shape and whose diameters disagree.
             max_x = cc_grid_max_x(forest_managed_loc)

             weibull_resolution = max_x / ncirc

             ! Guillaume M. -- Second-moment normalisation: the class multipliers are
             ! rescaled so that SUM(w_i * mod_qmd_i**2) = 1, hence the constructed stand
             ! has quadratic mean diameter qmd by construction. Without it n_init was
             ! computed on a diameter the stand did not have. The distribution shape is
             ! untouched; the weights are those of weibull_class_dist (exp_circ_class_n).
             DO icir = 1,ncirc
                mod_qmd_arr(icir) = un/ncirc + REAL(icir-1,r_std)*weibull_resolution
             ENDDO
             qmd_norm = SUM(exp_circ_class_n(:))
             IF (qmd_norm > min_stomate) THEN
                qmd_norm = SQRT(SUM(exp_circ_class_n(:)*mod_qmd_arr(:)**2) / qmd_norm)
                IF (qmd_norm > min_stomate) mod_qmd_arr(:) = mod_qmd_arr(:) / qmd_norm
             ENDIF

             DO icir = 1,ncirc
                ! Debug
                IF (printlev_loc.GE.4 .AND. write_debug) THEN
                   WRITE(numout,*) ' make_saplings icir:',icir
                END IF
                !-
                ! Make use of qmd_init (calculated in stomate.f90 as a parameter)
                ! do calculate the wood mass of the saplings of a newly
                ! planted PFT.

                CALL make_saplings(bm_sapl_2D_loc(icir,:,:), qmd * mod_qmd_arr(icir), &
                   ivm, KF_loc, &
                   c0_alloc, plant_status_loc, establish, cn_leaf_init_2D_loc, &
                   bypass_lm(icir), pipe_tune2_loc, write_debug)
                ! Guillaume M. -- mod_qmd_arr is built and normalised before this loop:
                ! normalisation is impossible while stepping through the classes.
             ENDDO

             ! Prescribe the newly planted population.
             ! qmd_init is the diameter of the tree
             ! in meters.  This is the diameter that we will use
             ! to calculate the new distribution of trees among
             ! the circumference classes.
             IF (forest_managed_loc .NE. ifm_cop) THEN
                ! Currently, the actual qmd after make_sapling is smaller than qmd_init. 
                ! n_init is calculated based on qmd_init (which is larger than the
                ! actual qmd) and rdi_start. Consequently, the resulting rdi is lower than 
                ! rdi_start (e.g., ~0.2), since n_init is determined to
                ! match qmd_init while the actual qmd is smaller. 
                ! As a result, rdi_start has no practical effect during parameterization. 
                ! The line below is a temporary fix to ensure the correct rdi
                ! at initialization. (when qmd is replaced to real_qmd in n_init
                ! calculation)
                real_qmd = wood_to_qmdia(bm_sapl_2D_loc(:,:,icarbon), &
                 exp_circ_class_n(:), ivm, pipe_tune2_loc)
                ! Guillaume M. -- A staggered older age class is a mid-rotation managed
                ! stand, not a freshly-planted one. At rdi_start a large qmd yields so few
                ! stems that density falls under dens_target and the stand is wiped by the
                ! density-collapse clearcut of sapiens_forestry. Desync classes therefore
                ! start at the management ceiling rdi_max, keeping density above target.
                iage_pr = -1     ! sentinelle : le bloc PRESCRIBE_RDI_FRAC n'a pas tourne
                rdi_use = rdi_start(ivm)
                IF (PRESENT(qmd_init_override)) THEN
                   IF (qmd_init_override < val_exp/2._r_std) rdi_use = rdi_max(ivm, forest_managed_loc)
                ENDIF
                ! Guillaume M. -- Establishment RDI is a fraction of rdi_max, set per age
                ! class: class 1 stays far below the self-thinning ceiling (a true gap),
                ! later classes stay close to it. This subsumes the desync case above,
                ! which is the fraction 1.0 special case. A fraction <= 0 leaves the value
                ! above untouched, so the block is bit-neutral when unset.
                IF (ALLOCATED(prescribe_rdi_frac) .AND. nagec > 1) THEN
                   iage_pr = ivm - start_index(agec_group(ivm)) + 1
                   IF (iage_pr >= 1 .AND. iage_pr <= nagec) THEN
                      IF (prescribe_rdi_frac(iage_pr) > zero) &
                           rdi_use = prescribe_rdi_frac(iage_pr) &
                                     * rdi_max(ivm, forest_managed_loc)
                   ENDIF
                ENDIF
                ! Guillaume M. -- Nmax must be evaluated on real_qmd, the diameter the
                ! stand actually gets from make_saplings, not on the larger target qmd:
                ! real_qmd was computed but never used, which pinned established RDI near
                ! 0.2 and made rdi_start and PRESCRIBE_RDI_FRAC inert. NOT bit-neutral:
                ! density rises threefold. rdi_max stays below inventory RDI -- other job.
                IF (real_qmd > min_stomate) THEN
                   n_init = MIN(Nmax(real_qmd * m_to_cm, &
                      alpha_self_thin_loc, ivm, forest_managed_loc) * rdi_use, &
                      nmaxplants(ivm))
                ELSE
                   ! Guillaume M. -- Degenerate real_qmd (zero sapling mass): fall back
                   ! on qmd rather than divide by a null diameter.
                   n_init = MIN(Nmax(qmd * m_to_cm, &
                      alpha_self_thin_loc, ivm, forest_managed_loc) * rdi_use, &
                      nmaxplants(ivm))
                ENDIF

             ELSE
                n_init = MIN(Nmax(coppice_diameter(ivm)*m_to_cm, &
                   alpha_self_thin_loc, ivm, forest_managed_loc) * rdi_max(ivm, forest_managed_loc), &
                   nmaxplants(ivm))
             END IF


             n_init = n_init*ha_to_m2
             circ_class_n_loc(:) = exp_circ_class_n(:) / SUM(exp_circ_class_n(:)) * n_init

             circ_class_biomass_loc(:,:,:) = bm_sapl_2D_loc(:,:,:)

             ! Debug
             IF (printlev_loc.GE.4 .AND. write_debug) THEN
                DO icir=1,ncirc
                   WRITE(numout,*)'icir, circ_class_n, ', icir, &
                        circ_class_n_loc(icir)
                   WRITE(numout,*)'Initial biomass distribution - carbon'
                   WRITE(numout,*) SUM(circ_class_biomass_loc(icir,:,icarbon))
                   WRITE(numout,*) SUM(circ_class_biomass_loc(icir,:,initrogen))
                ENDDO
                WRITE(numout,*)' qmd', qmd
                WRITE(numout,*)' n_init(ivm)', n_init
                WRITE(numout,*)' Maximun number of saplings at RDI=1 ', &
                     Nmax(qmd_init(ivm)*m_to_cm, alpha_self_thin_loc, ivm, forest_managed_loc)
             END IF
             !-

          ELSEIF (recruite) THEN

             !! Calculate the number of succesful recruits
             ! We want to avoid adding a lot of recruits which are then killed a bit 
             ! later in the simulation. Such an approach would also be in conflict
             ! with the approach in which we keep a constant number of diameter
             ! and thus have to merge the biomass of the recruits with the biomass
             ! in the smallest diameter class. For that reason we calculate the
             ! the number of succesful recruits that will survive and will make it
             ! into the next diameter class. This number is expected to be rather low
             ! and cannot be compared against observations of actual recruitment. 
             ! It should be compared against observations of in growth.
             
             ! Recruitment is based on the light stress factor [0 1] computed
             ! above. It makes use of the functional response of recruitment
             ! to light described by Ruger et al (2009, J. of Ecol.,
             ! doi: 10.1111/j.1365-2745.2009.01552.x) in the 50-ha plot of Barro
             ! Colorado Island (BCI) in Panama: 
             ! Log10(Nrecruits)= A + B * ( Log10(lstress_fac) - MU )
             !
             ! Where Nrecruits is the number of recruits for an area of 25 m2.
             ! The parameter A is the intercept (mean log10 recruits expected), 
             ! B the slope (functional response), and MU=Log10(0.02)=-1.7 
             ! Here the recruitment response to light can be negative (B<0),
             ! decelerating (0<B<1), accelerating (B>1) or linear (B=1). 
             ! MU is needed here to express the fraction of light between 0 and 1
             ! centred on mean light 2% at BCI. The linearized (log-log) power
             ! model is implemented here for each PFT following the original setup
             ! of Ruger et al 2009 to ensure realistic values. IT COULD EASILY BE
             ! SIMPLIFIED. The predicted number of recruits is for an area of
             ! 25 m2 so we have to compute the antilogarithm and divide by 25 to
             ! get the estimate by 1 m2 which is the unit used in ORCHIDEE.
             IF (recruitment_alpha(ivm).LT.min_stomate .AND. &
                  recruitment_beta(ivm).LT.min_stomate) THEN
                
                ! Sanity check. If the users sets both values to zero
                ! new_ind_loc is set to zero. The simulation should be 
                ! identical to a simulation without recruitment.
                new_ind_loc = zero
                
             ELSE
                      
                ! Calculate the actual number of trees before recruitment
                tmp_ind = SUM(circ_class_n_loc(:))
                
                ! Calculate the number of recruits m2 day-1
                ! New approach in which stand density (tmp_ind) is the main driver which
                ! is further modulated by light. This approach simulates fewer recruits
                ! prior to canopy closure. Als the response to a thinning depends on 
                ! the stand structure following the thinning. This approach seems a bit
                ! as double counting because the light should be a consequence of the stand
                ! density.
                new_ind_loc = recruitment_alpha(ivm)**(LOG10(recruitment_beta(ivm)*tmp_ind)) * &
                     SQRT(LOG10(MAX(100.*lstress_fac,un)))
                
             ENDIF

             ! Update the number of individuals
             old_ind = tmp_ind
             tmp_ind = tmp_ind + new_ind_loc

             ! Debug
             IF(printlev_loc.GE.4 .AND. write_debug)THEN
                WRITE(numout,*) 'canopy cover (LAI), ',&
                     cc_to_lai(circ_class_biomass_loc(:,ileaf,icarbon),&
                     circ_class_n_loc(:),ivm) 
                WRITE(numout,*) 'used canopy_gap (%) is, ', lstress_fac*100 
                WRITE(numout,*) 'Original ind is, ', (tmp_ind - new_ind_loc)*m2_to_ha
                WRITE(numout,*) 'new_ind is, ', new_ind_loc*m2_to_ha
                WRITE(numout,*) 'ind+new_ind is, ', tmp_ind*m2_to_ha
                WRITE(numout,*) 'Dg (m) ', qmd_init(ivm)
                ind_init = Nmax(qmd_init(ivm)*m_to_cm, alpha_self_thin_loc, ivm, forest_managed_loc)
                WRITE(numout,*) 'Nmax:, ',ind_init
             ENDIF
             !-

             !! Calculate the dimensions of the recruits    
             ! The height of the saplings is prescribed and determines the reserves
             ! which are especially important for deciduous species, which need to
             ! survive on their reserves for the first couple of months to a whole
             ! year because the phenology scheme requires annual mean values to get
             ! started. Calculate sapling diameter and biomass for prescribed sapling
             ! height (above). Note that the dimensions of a newly planted sapling and 
             ! a recruit are not the same.
             dia_init_recr = ( recruitment_height(ivm) / pipe_tune2_loc ) ** &
                  ( 1. / pipe_tune3(ivm) )

             ! Use the variables from the pipe model and the allometric
             ! relationships to calculate the different biomass components
             ! of the saplings in each diameter class.
             CALL make_saplings(bm_sapl_2D_loc(1,:,:), dia_init_recr, ivm, KF_loc, &
                  c0_alloc, plant_status_loc, establish, cn_leaf_init_2D_loc, &
                  bypass_lm(1), pipe_tune2_loc, write_debug)
             
             ! Debug
             IF (printlev_loc.GE.4 .AND. write_debug) THEN
                WRITE(numout,*)'Prescribed recruitment height (m): ',recruitment_height(ivm)
                WRITE(numout,*)'Recruitment diameter (m) for prescribed height: ', dia_init_recr
                WRITE(numout,*)'bm_sapl after calling make_saplings: ',bm_sapl_2D_loc(1,:,:)
             ENDIF
             !-

             ! Update the number of individuals in the first size class. All
             ! recruitment occurs in the circ_class with the smallest diameter
             old_circ_class_n(:) = circ_class_n_loc(:)
             circ_class_n_loc(1) = circ_class_n_loc(1) + new_ind_loc                        

             !! Determine the biomass in the first circumference class.
             ! This is a dillution approach were an average tree is made by making use
             ! of the biomass in the available trees and the biomass in the recruits.
             ! This is a simple mean weighted by the number of trees. The model needs
             ! the biomass in an individual tree (gC tree-1)
             DO iele = 1,nelements
                old_circ_class_biomass(:,:,iele) = &
                     circ_class_biomass_loc(:,:,iele) 
                circ_class_biomass_loc(1,:,iele) = &
                     (circ_class_biomass_loc(1,:,iele) * &
                     old_circ_class_n(1) + &
                     bm_sapl_2D_loc(1,:,iele) * new_ind_loc ) / &
                     circ_class_n_loc(1)
             ENDDO

             ! Debug
             IF (printlev_loc.GE.4 .AND. write_debug) THEN
                WRITE(numout,*) "Canopy gap (%)", lstress_fac*100
                WRITE(numout,*) "Original density (N/ha): ",&
                     old_circ_class_n*m2_to_ha
                WRITE(numout,*) "New density (N/ha): ",circ_class_n_loc*m2_to_ha
                WRITE(numout,*) "Recruited saplings per ha in first size class ", &
                     new_ind_loc * m2_to_ha
                WRITE(numout,*) "Original biomass per individual", &
                     SUM(old_circ_class_biomass(:,:,icarbon),2),&
                     SUM(old_circ_class_biomass(:,:,initrogen),2)
                WRITE(numout,*) "New biomass per individual", &
                     SUM(circ_class_biomass_loc(:,:,icarbon),2),&
                     SUM(circ_class_biomass_loc(:,:,initrogen),2)
                WRITE(numout,*) "Biomass of new saplings gC/sapling", &
                     bm_sapl_2D_loc(1,:,icarbon)
                WRITE(numout,*) "Biomass of new saplings gN/sapling", &
                     bm_sapl_2D_loc(1,:,initrogen)
                WRITE(numout,*) "KF ", KF_loc
             ENDIF
             !-

          ENDIF ! end prescribe/recruitment block

          IF (establish) THEN

             ! Set leaf age classes, all leaves are current year leaves
             leaf_frac_loc(:) = zero
             leaf_frac_loc(1) = un
             
             ! Set time since last beginning of growing season but only
             ! for the first day of the whole simulation. When the model
             ! is initialized when_growthinit is set to undef. In subsequent
             ! time steps it should have a value. For trees without phenology
             ! the growing season starts at the moment the PFT is prescribed
             IF (when_growthinit_loc .EQ. undef) THEN
                when_growthinit_loc = 200
             ENDIF
             
          ENDIF

          ! The carbon and nitrogen to build the saplings is taken from the 
          ! atmosphere, keep track of amount to calculate the C-balance 
          ! closure
          IF (establish) THEN

             ! The whole biomass of this PFT was taken from the
             ! atmosphere
             atm_to_bm_loc(:) = atm_to_bm_loc(:) + &
                (SUM(cc_to_biomass(1, ivm, &
                circ_class_biomass_loc(:,:,:), &
                circ_class_n_loc(:)),1) / dt )
         
          ELSEIF (recruite) THEN

             ! Only the increase in biomass of this PFT was taken
             ! from the atmosphere
             atm_to_bm_loc(:) = atm_to_bm_loc(:) + &
                  (SUM(cc_to_biomass(1, ivm, &
                  circ_class_biomass_loc(:,:,:), &
                  circ_class_n_loc(:)) - &
                  old_biomass(:,:),1)  / dt )
             
          ENDIF

          ! Check whether the order was preserved. We only want to
          ! preserve the order for the carbon biomass. Nitrogen follows
          ! carbon. We expect that the biomass of an individual tree 
          ! increases with an increasing circ class. This is probably 
          ! not important for the correct functioning of the model. So,
          ! if this turns out to be computationally expensive it is 
          ! worth trying without. We try to preserve the order as
          ! an additional mean to control the quality of the simulations.
          ! A strict order also helps to interpret the output. Hence, 
          ! if the order was not preserved we will sort the biomasses. 
          ! Note that this subroutine first checks whether the order 
          ! was preseverd.
          circ_class_biomass_new = circ_class_biomass_loc(:,:,:)
          circ_class_n_new = circ_class_n_loc(:)
          CALL sort_circ_class_biomass(circ_class_biomass_new, &
               circ_class_n_new)
          circ_class_biomass_loc(:,:,:) = circ_class_biomass_new
          circ_class_n_loc(:) = circ_class_n_new

       ELSEIF( ( .NOT. is_tree(ivm) ) .AND. &
            ( natural(ivm) ) .AND. &
            ( veget_max_loc .GT. min_stomate ) ) THEN

          IF (establish) THEN

             !! 2.3.2 Initialize grassy PFTs
             !  Use veget_max_loc to check whether the PFT is present. If the PFT is 
             !  present but it has no biomass, prescribe its biomass (gC m-{2}). 
             !  It is assumed that at day 1 grasses have all their biomass in the 
             !  reserve pool. The criteria exclude crops. Crops are no longer 
             !  prescribed but planted the day that begin_leaves is true.
             
             ! Debug
             IF(printlev_loc.GE.4 .AND. write_debug)THEN
                WRITE(numout,*) 'We will prescribe a new grass  '//&
                     'vegetation, the old one died'
             ENDIF
             ! -

             ! For grasses, we assume that newly established density is
             !  based on nmaxplants in pft_parameters.f90.
             tmp_ind = nmaxplants(ivm) * ha_to_m2

             ! With flag activated, if the established grassland is new
             ! circ_class_n_grass=0 so it shouldn't go
             ! through the following calculation, in this case, newly
             ! established density is still based on nmaxplants
             ! Only in case of lcchange, to avoid crash,
             ! the newly established grass density should be identical to that 
             ! in already existed grassland, where circ_class_n_grass should
             ! have a non-zero value.
             IF(ok_dyn_grass_den) THEN
                IF(circ_class_n_grass(1) .GT. zero ) THEN
                   tmp_ind = circ_class_n_grass(1)
                END IF
             END IF

             ! initialize everything to make sure there are not random values 
             ! floating around for ncirc = 1
             bm_sapl_2D_loc(1,:,:) = val_exp

             ! Similar as for trees, the initial height of the vegetation was 
             ! defined. By using lai_to_biomass dynamic sla is accounted for. 
             ! Root and sapwood mass calculated below.
             ! Leaf biomass calculated first to calculate KF and LF.
             bm_sapl_2D_loc(1,ileaf,icarbon) = lai_to_biomass(&
                height_init(ivm)/lai_to_height(ivm),ivm)

             !+++CHECK+++
             ! We do not deal with circumference classes for grasses
             ! and crops, but we want to keep the arrays the same as for the trees 
             ! so we put the young grass information into the first circumference 
             ! class. Calculate the sapwood to leaf mass in a similar way as has 
             ! been done for trees. For trees this approach had been justified by 
             ! observations. For grasses such justification is not supported by 
             ! observations but we didn't try to find it. Needs more work by 
             ! someone interested in grasses. There might be a more elegant
             ! solution making use of a well observed parameter. 
             ! The mass of the structural carbon relates to the mass of the leaves 
             ! through a prescribed parameter ::k_latosa
             !+++CHECK+++
             ! At one point it looked like a good idea to take the max of two
             ! options but by doing so we cannot recalculate KF in phenology.
             ! It has not been confirmed that this is really a problem. Use the
             ! simpelest approach but leave the alternative in the code as a
             ! suggestion of a possible solution in case something goes wrong.
!!$             k_latosa_tmp = MAX(k_latosa_min(ivm),k_latosa_adapt_loc + &
!!$                  (lstress_fac * &
!!$                  (k_latosa_max(ivm)-k_latosa_min(ivm))))
             !+++++++++++
             k_latosa_tmp = k_latosa_adapt_loc + lstress_fac * &
                  (k_latosa_max(ivm)-k_latosa_min(ivm))

             IF (sla_dyn) THEN 
                sla_est = biomass_to_lai(bm_sapl_2D_loc(1,ileaf,icarbon),ivm)/&
                        bm_sapl_2D_loc(1,ileaf,icarbon)
             ELSE
                sla_est = sla(ivm)
             END IF

             KF_loc = k_latosa_tmp / &
                (sla_est * tree_ff(ivm) * pipe_density(ivm))

             ! Calculate leaf to root area    
             LF = c0_alloc * KF_loc

             ! Debug
             IF(printlev_loc.GE.4 .AND. write_debug)THEN
                WRITE(numout,*) 'KF, ', KF_loc
                WRITE(numout,*) 'LF, c0_alloc, ', c0_alloc * KF_loc, &
                     c0_alloc
             ENDIF

             ! Similar as for trees, the initial height of the vegetation was 
             ! defined. By using lai_to_biomass dynamic sla is accounted for.
             IF (test_biomass_init) THEN 
                bm_sapl_2D_loc(1,ileaf,icarbon) = biomass_init(ivm)
             ELSE
                bm_sapl_2D_loc(1,ileaf,icarbon) = lai_to_biomass(&
                     height_init(ivm)/lai_to_height(ivm),ivm)
             ENDIF

             ! Use allometric relationships to define the root mass based on 
             ! leaf mass. Some 'sapwood mass' is needed to store the reserves. 
             ! An arbitrairy fraction of 5% was used 
             bm_sapl_2D_loc(1,iroot,icarbon) = &
                  bm_sapl_2D_loc(1,ileaf,icarbon) / LF
             bm_sapl_2D_loc(1,isapabove,icarbon) = &
                  bm_sapl_2D_loc(1,ileaf,icarbon) / KF_loc

             ! Some of the biomass components that exist for trees are undefined 
             ! for grasses
             bm_sapl_2D_loc(1,isapbelow,icarbon) = zero
             bm_sapl_2D_loc(1,iheartabove,icarbon) = zero
             bm_sapl_2D_loc(1,iheartbelow,icarbon) = zero
             bm_sapl_2D_loc(1,ifruit,icarbon) = zero
             bm_sapl_2D_loc(1,icarbres,icarbon) = zero

             ! Pools that are defined in the same way for trees and grasses
             bm_sapl_2D_loc(1,ilabile,icarbon) = labile_to_total * &
                  (bm_sapl_2D_loc(1,ileaf,icarbon) + &
                  fcn_root(ivm) * bm_sapl_2D_loc(1,iroot,icarbon) + fcn_wood(ivm) * &
                  (bm_sapl_2D_loc(1,isapabove,icarbon) + &
                  bm_sapl_2D_loc(1,isapbelow,icarbon) + &
                  bm_sapl_2D_loc(1,icarbres,icarbon)))

             ! Allocate the nitrogen
             cn_leaf=cn_leaf_init_2D_loc 
             cn_wood=cn_leaf/fcn_wood(ivm) 
             cn_root=cn_leaf/fcn_root(ivm) 

             !+++CHECK+++
             ! See comments at a similmar block of code for the
             ! tree saplings
             ! Use the C/N ratio to calculate the N content for
             ! all the biomass components.
             bm_sapl_2D_loc(1,:,initrogen) = &
                  bm_sapl_2D_loc(1,:,icarbon) / cn_root

             ! Overwrite the N content for the leaves and wood
             ! components. Note that only sapabove is defined the
             ! other wood components are zero
             bm_sapl_2D_loc(1,ileaf,initrogen) = &
                  bm_sapl_2D_loc(1,ileaf,icarbon) / cn_leaf
             bm_sapl_2D_loc(1,isapabove,initrogen) = &
                  bm_sapl_2D_loc(1,isapabove,icarbon) / cn_wood
             !++++++++++++
             
             ! Avoid deciduous PFTs to have leaves out at establishment
             ! Whether the saplings have leaves or don't have leaves the first 
             ! year doesn't really matter. Either the approach is correct in the 
             ! northern hemisphere or in the southern hemisphere. Anyhow, a 
             ! spin-up is recommended to avoid issues with the initial conditions
             IF ( pheno_type(ivm) .NE. 1 ) THEN

                ! Not evergreen. Deciduous PFTs now need to survive half a year 
                ! before bud burst. To ensure survival there are several options: 
                ! (a) either the height of the initial vegetation is increased 
                ! (this results in more reserves) or (b) the reserves could be
                ! increased. The second option may result in result in numerical 
                ! issues further down the code as the optimal reserve level is 
                ! calculated from the other biomass pools
                
                ! Notice that the sapwood is not put into reserves.  Sapwood is
                ! important in phenology when the amount of leaves and roots
                ! are determined, so it needs to be already present.
                bm_sapl_2D_loc(1,icarbres,:) = bm_sapl_2D_loc(1,icarbres,:) + &
                     bm_sapl_2D_loc(1,ileaf,:) + &
                     bm_sapl_2D_loc(1,iroot,:)
                bm_sapl_2D_loc(1,ileaf,:) = zero
                bm_sapl_2D_loc(1,iroot,:) = zero
                
                ! When there are no leaves, the crop/grass is dormant
                plant_status_loc = idormant
                
             ELSE

                ! Initilize carbohydrate reserves for evergreen PFTs
                bm_sapl_2D_loc(1,icarbres,:) = zero 

                ! Evergreen plants always have a canopy
                plant_status_loc = icanopy

             ENDIF

             ! Debug
             IF (printlev_loc.GE.4 .AND. write_debug) THEN
                DO iele=1,nelements
                   WRITE(numout,*) ' young grass biomass (gC):',ivm
                   WRITE(numout,*) '  bm_sapl_2D(ileaf) ',&
                        bm_sapl_2D_loc(1,ileaf,iele)
                   WRITE(numout,*) '  bm_sapl_2D(ispabove) ',&
                        bm_sapl_2D_loc(1,isapabove,iele)
                   WRITE(numout,*) '  bm_sapl_2D(isapbelow) ',&
                        bm_sapl_2D_loc(1,isapbelow,iele)
                   WRITE(numout,*) '  bm_sapl_2D(iheartabove) ',&
                        bm_sapl_2D_loc(1,iheartabove,iele)
                   WRITE(numout,*) '  bm_sapl_2D(iheartbelow) ',&
                        bm_sapl_2D_loc(1,iheartbelow,iele)
                   WRITE(numout,*) '  bm_sapl_2D(iroot) ',&
                        bm_sapl_2D_loc(1,iroot,iele)
                   WRITE(numout,*) '  bm_sapl_2D(ifruit,)', &
                        bm_sapl_2D_loc(1,ifruit,iele)
                   WRITE(numout,*) '  bm_sapl_2D(icarbres,)', &
                        bm_sapl_2D_loc(1,icarbres,iele)
                   WRITE(numout,*) '  bm_sapl_2D(ilabile,)', &
                        bm_sapl_2D_loc(1,ilabile,iele)
                ENDDO
             ENDIF
             ! -

             ! circ_class_biomass (gC tree-1). This 
             ! allows to more easily check the mass balance closure
             circ_class_biomass_loc(1,:,:) = bm_sapl_2D_loc(1,:,:)
             circ_class_n_loc(1) = tmp_ind

             ! Set leaf age classes -> all leaves will be current year leaves
             leaf_frac_loc(:) = zero
             leaf_frac_loc(1) = un

             ! Set time since last beginning of growing season but only
             ! for the first day of the whole simulation. When the model
             ! is initialized when_growthinit is set to undef. In subsequent
             ! time steps it should have a value.
             IF (when_growthinit_loc .EQ. undef) THEN
                when_growthinit_loc = 200
             ENDIF
                   
             ! The carbon and nitrogen biomass to build the saplings is taken 
             ! from the atmosphere, keep track of amount to calculate the mass 
             ! balance closure
             DO iele = 1,nelements
                atm_to_bm_loc(iele) = atm_to_bm_loc(iele) + &
                     ((SUM(circ_class_biomass_loc(1,:,iele))*&
                     circ_class_n_loc(1)) / dt)
             ENDDO

             ! Debug
             IF (printlev_loc.GE.4 .AND. write_debug) THEN
                DO iele=1,nelements
                   WRITE(numout,*) ' root to sapwood tradeoff (LF) : ', c0_alloc
                   WRITE(numout,*) ' grass sapling biomass: ',SUM(bm_sapl_2D_loc(1,:,iele))
                ENDDO
             ENDIF
             ! -
           
          ELSEIF (recruite) THEN
             
             ! Debug
             IF (printlev_loc .GE. 4) THEN
                WRITE(numout,*) 'Before recruite, ivm', ivm
                WRITE(numout,*) 'Density before recruite', circ_class_n_loc(1)
                WRITE(numout,*) 'biomass carbon before recruite', &
                                circ_class_biomass_loc(1,:,icarbon)
             END IF
             !-
             
             ! The theory to increase density through asexual and vegetative
             ! reproduction is based on sufficient fruit, reserves, and
             ! labile carbon. We assumed that after (vegetative) reproduction
             ! the areal reserve carbon and labile carbon should equal their
             ! respective targets (gC m-2), and the fruit carbon pool is zero because
             ! the seeds have been used.
             ! 1. Mass balance should be preserved before and after the density
             ! change, so we have:
             ! N0*Call = N1*Cother + N1*Trl1
             ! Call is carbon from all compartment (gC ind-1) and Trl1 is sum of new
             ! reserve+labile carbon:
             ! Trl1=(scale_inc*N0/N1)*Trltarget_area/N0,
             ! combining the two, new density is computed as:
             ! N1=(N0*Call-scale_inc*Ttarget_area)/(Cother)
             ! Trl0 = individual reserve+labile
             Trl0=circ_class_biomass_loc(1,ilabile,icarbon) + &
                  circ_class_biomass_loc(1,icarbres,icarbon)

             ! The individual biomass from other compartments should be the
             ! same before and after dynamic grass density, so we calculate
             ! it as: sum of individual biomass C from all compartments
             ! – individual labile C – individual reserve C  - individual
             ! fruit C, represented in Cother, with unit of gC ind-1
             Cother= SUM(circ_class_biomass_loc(1,:,icarbon)) &
                  -circ_class_biomass_loc(1,ifruit,icarbon)-&
                  circ_class_biomass_loc(1,ilabile,icarbon) &
                  -circ_class_biomass_loc(1,icarbres,icarbon)

             ! Keep old density and biomass
             old_circ_class_n(:) = circ_class_n_loc(:)
             old_circ_class_biomass(:,:,:) = circ_class_biomass_loc(:,:,:)

             ! calculate new density using Nume/Deno
             Nume=circ_class_n_loc(1)*SUM(circ_class_biomass_loc(1,:,icarbon))- &
                 scale_inc*(tot_res_target_loc(icarbon)+tot_lab_target_loc(icarbon))
             Deno=Cother
             IF ((Nume .GT. min_stomate) .AND. &
                (Deno .GT. min_stomate)) THEN
                circ_class_n_loc(1) = Nume/Deno
             END IF

             ! Update the individual biomass for reserve and labile carbon, and
             ! judge whether density exceed 1.
             IF (circ_class_n_loc(1) .GT. 1) THEN

                ! If density exceeds the maximum which is 1, keep it at 1.
                circ_class_n_loc(1) = 1

                ! in this case, to keep mass balance, thus
                ! ex_bio=N0*Call/N1-Cother 
                ex_bio= MAX(min_stomate,((SUM(old_circ_class_biomass(1,:,icarbon))* &
                     old_circ_class_n(1))/circ_class_n_loc(1) - Cother))

                ! Calculate ratio_res and ratio_lab
                ratio_res = old_circ_class_biomass(1,icarbres,icarbon) / &
                     (old_circ_class_biomass(1,icarbres,icarbon) + &
                     old_circ_class_biomass(1,ilabile,icarbon))
                ratio_lab = un - ratio_res 

                ! Calculate individual reserve or labile carbon based on ex_bio
                ! and ratio
                circ_class_biomass_loc(1,icarbres,icarbon) = ratio_res * &
                     ex_bio    
                circ_class_biomass_loc(1,ilabile,icarbon) = ratio_lab * &
                     ex_bio
                
             ELSE
                
                ! In the case of density<=1, Trl1=scale_inc*Ttarget_area/N1
                circ_class_biomass_loc(1,icarbres,icarbon) = tot_res_target_loc(1)* &
                     scale_inc /circ_class_n_loc(1)
                circ_class_biomass_loc(1,ilabile,icarbon) = tot_lab_target_loc(1)* &
                     scale_inc /circ_class_n_loc(1)

             END IF ! density > 1

             ! Update circ_class_biomass of fruit to 0.
             circ_class_biomass_loc(1,ifruit,icarbon)=zero

             ! Dilute the circ_class_biomass for nitrogen for mass balance,
             ! according to the ratio of old density to new density
             ratio_den = old_circ_class_n(1)/circ_class_n_loc(1)
             DO ipart = 1,nparts
                circ_class_biomass_loc(1,ipart,initrogen) = &
                     circ_class_biomass_loc(1,ipart,initrogen)*ratio_den
             ENDDO ! ipart = 1,nparts

             ! Debug
             IF (printlev_loc .GE. 4) THEN
                WRITE(numout,*) 'After recruite, ivm', ivm
                WRITE(numout,*) 'Density after recruite', circ_class_n_loc(1)
                WRITE(numout,*) 'biomass carbon after recruite', &
                     circ_class_biomass_loc(1,:,icarbon)
                WRITE(numout,*) 'reserve target', &
                     tot_res_target_loc
                WRITE(numout,*) 'labile target', &
                     tot_lab_target_loc
             END IF
             !-
             
           END IF ! establish or recruite

       ELSEIF (.NOT. natural(ivm)) THEN

          IF (establish) THEN

             
             ! Initialize croplands - else leaves_begin will never become true
             ! Set time since last beginning of growing season but only
             ! for the first day of the whole simulation. When the model
             ! is initialized when_growthinit is set to undef. In subsequent
             ! time steps it should have a value.
             IF (when_growthinit_loc .EQ. undef) THEN
                when_growthinit_loc = 200
             ENDIF
             
             ! Debug
             IF (printlev_loc.GE.4 .AND. write_debug) WRITE(numout,*) 'Initialize crops',ivm
             !-

          ELSEIF (recruite) THEN

             ! No recruitment in croplands
             
          ENDIF ! establish or recruite
          
       ENDIF ! forest, grassland or cropland
       
       !! 2.3.3 Declare PFT present
       !  Now that the PFT has biomass it should be declared 'present'
       !  everywhere in that grid box. Assign some additional properties
       IF (establish) THEN
          PFTpresent_loc = .TRUE.
          everywhere_loc = un
          age_loc = zero
          npp_longterm_loc = npp_longterm_init
          lm_lastyearmax_loc = SUM(bm_sapl_2D_loc(:,ileaf,icarbon) * &
               circ_class_n_loc(:))
       ENDIF

    ENDIF ! establish or recruite

 !! 3. Check numerical consistency of this routine

    IF (err_act.GT.1) THEN

       ! 3.1 Check surface area
       CALL check_vegetation_area_point("stomate_prescribe", veget_max_begin, &
            veget_max_loc)

       ! 3.2 Mass balance closure 
       ! 3.2.1 Calculate final biomass
       pool_end(:) = zero 
       DO ipar = 1,nparts
          DO iele = 1,nelements
             DO icir = 1,ncirc
                pool_end(iele) = pool_end(iele) + &
                     (circ_class_biomass_loc(icir,ipar,iele) * &
                     circ_class_n_loc(icir) * veget_max_loc)
             ENDDO
          ENDDO
       ENDDO

       ! 3.2.2 Calculate mass balance
       ! Common processes
       DO iele=1,nelements
          check_intern(iatm2land,iele) = check_intern(iatm2land,iele) + &
               atm_to_bm_loc(iele) * veget_max_loc * dt
          check_intern(ipoolchange,iele) = -un * (pool_end(iele) - &
               pool_start(iele)) 
       ENDDO

       closure_intern(:) = zero
       DO imbc = 1,nmbcomp
          DO iele=1,nelements
             closure_intern(iele) = closure_intern(iele) + &
                  check_intern(imbc,iele)
          ENDDO
       ENDDO

       ! 3.3 Check mass balance closure
       CALL check_mass_balance_point("stomate_prescribe", closure_intern, &
            pool_end, pool_start, veget_max_loc)
      
    ENDIF ! err_act .GT. 1

    IF(printlev_loc>=3) WRITE(numout,*) 'Leaving stomate_prescribe.f90'

  END SUBROUTINE prescribe



!! ================================================================================================================================
!! SUBROUTINE   : make_sapling
!!
!>\BRIEF Given the dimensions and the allocation factors, calculate the biomass of a sapling        
!!
!! DESCRIPTION (functional, design, flags):
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLES(S): ::bm_sapl, ::KF
!!
!! REFERENCES   :
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE make_saplings(bm_sapl_2D_loc, tree_dia, ivm, KF_loc, &
       c0_alloc, plant_status_loc, establish, cn_leaf_init_2D_loc, &
       bypass_lm, pipe_tune2_loc, write_debug)

    !! 0. Parameters and variables declaration

    !! 0.1 Input variables
    LOGICAL, INTENT(in)                                :: write_debug        !! flag to write debug statements in the prescribe subroutine
    REAL(r_std), INTENT(in)                            :: tree_dia           !! Tree diameter (m)
    INTEGER(i_std), INTENT(in)                         :: ivm                !! index (unitless)
    REAL(r_std), INTENT(in)                            :: c0_alloc           !! Root to sapwood tradeoff parameter
    LOGICAL, INTENT(in)                                :: establish          !! yes: establish a new young vegetation. 
                                                                             !! No: add saplings under an existing vegetation
    REAL(r_std), INTENT(in)                            :: cn_leaf_init_2D_loc!! initial leaf C/N ratio 
    REAL(r_std), INTENT(in)                            :: KF_loc             !! Scaling factor to convert sapwood mass
                                                                             !! into leaf mass (m) - this variable is 
                                                                             !! passed to other routines
    REAL(r_std), INTENT(in)                            :: pipe_tune2_loc     !! pipe_tune2
    
    !! 0.2 Output variables
    REAL(r_std), INTENT(out)                           :: bypass_lm          !! leaf mass before it is set to zero
                                                                             !! for deciduous PFTs. This value is
                                                                             !! in the calculation of KF (gC tree-1)

    !! 0.3 Modified variables
    REAL(r_std), DIMENSION(:,:),INTENT(inout)          :: bm_sapl_2D_loc     !! Sapling biomass for the functional 
    REAL, INTENT(inout)                                :: plant_status_loc   !! Growth and phenological status of the plant 
                                                                             !! Different stati defined in constants
    !! 0.4  Local variables
    INTEGER(i_std)                                    :: iele                !! index (unitless) 
    REAL(r_std)                                       :: cn_leaf,cn_root     !! CN ratio of leaves and roots (gC/gN)
    REAL(r_std)                                       :: cn_wood             !! CN ratio of wood pool (gC/gN)        
    REAL(r_std)                                       :: tree_height         !! Tree height (m)
    REAL(r_std)                                       :: wood_above_loc      !! Total aboveground woody mass to be split sap/heart (gC tree-1)
    REAL(r_std)                                       :: f_heart_loc         !! Heartwood mass fraction of the woody cross-section (-)
    REAL(r_std), PARAMETER                            :: sap_depth = 0.05_r_std         !! Sapwood radial depth (m), ~5 cm temperate average
    REAL(r_std), PARAMETER                            :: dia_heartwood_min = 0.05_r_std !! Diameter (m) above which heartwood is partitioned explicitly (cold-start desync)

!_ ================================================================================================================================

      ! The woody biomass is contained in four components. Thus, 
      ! wood_init = isapabove + isapbelow + iheartabove + iheartbelow. 
      ! Given that isapbelow = isapbelow and iheartabove = 
      ! iheartbelow =  bm_sapl_heartabove*isapabove. If 
      ! bm_sapl_heartbelow = bm_sapl_heartabove = 0.2, then 
           
      IF (qmd_init(ivm) .GT. 0.05) THEN

        CALL ipslerr_p(3,'Initial diameter of sapling was set to more than 5cm',&
        'The current estimate of a sapling biomass does not account &
        for heartwood','If you want to use larger diameters add a &
        heartwood estimate in the code','')

      ELSE

         ! Inialize biomass of a sapling. Taking the same biomass
         ! above and belowground is a bit arbitrairy. It will result
         ! in slightly taller recruits than precribed
         bm_sapl_2D_loc(isapabove,icarbon) = tree_dia**(2+pipe_tune3(ivm))*pi/4* &
            tree_ff(ivm)*pipe_density(ivm)*pipe_tune2_loc
         ! We need something in the heartwood
         bm_sapl_2D_loc(iheartabove,icarbon) =  0.01 * bm_sapl_2D_loc(isapabove,icarbon) 
         bm_sapl_2D_loc(iheartbelow,icarbon) =  bm_sapl_2D_loc(iheartbelow,icarbon)
         ! To maintain consistency between mass and volume, subtract 
         ! heartwood mass.
         bm_sapl_2D_loc(isapabove,icarbon) = tree_dia**(2+pipe_tune3(ivm))*pi/4* &
            tree_ff(ivm)*pipe_density(ivm)*pipe_tune2_loc - &
            bm_sapl_2D_loc(iheartabove,icarbon)
         bm_sapl_2D_loc(isapbelow,icarbon) = bm_sapl_2D_loc(isapabove,icarbon)

         ! Guillaume M. -- Re-partition sap/heartwood. The split above leaves a fixed 1%
         ! heartwood, unphysical at the large establishment diameters: a 40 cm tree would
         ! be 99% sapwood. Sapwood is taken as the outer ring of radial depth sap_depth
         ! and heartwood as the inner core, f_heart = ((D - 2*sap_depth)/D)**2, floored at
         ! 0.01. Gated on ok_dia_stagger, so every other call keeps the legacy split.
         IF (ok_dia_stagger .AND. tree_dia .GT. dia_heartwood_min) THEN
            wood_above_loc = bm_sapl_2D_loc(isapabove,icarbon) + &
                 bm_sapl_2D_loc(iheartabove,icarbon)
            IF (tree_dia .GT. 2._r_std*sap_depth) THEN
               f_heart_loc = ((tree_dia - 2._r_std*sap_depth)/tree_dia)**2
            ELSE
               f_heart_loc = zero
            ENDIF
            f_heart_loc = MAX(f_heart_loc, 0.01_r_std)
            bm_sapl_2D_loc(iheartabove,icarbon) = f_heart_loc * wood_above_loc
            bm_sapl_2D_loc(iheartbelow,icarbon) = bm_sapl_2D_loc(iheartabove,icarbon)
            bm_sapl_2D_loc(isapabove,icarbon)   = (un - f_heart_loc) * wood_above_loc
            bm_sapl_2D_loc(isapbelow,icarbon)   = bm_sapl_2D_loc(isapabove,icarbon)
         ENDIF

         ! Error checking
         IF (bm_sapl_2D_loc(isapabove,icarbon) .LT. zero .OR. &
              bm_sapl_2D_loc(isapbelow,icarbon) .LT. zero .OR. &
              bm_sapl_2D_loc(iheartabove,icarbon) .LT. zero .OR. &
              bm_sapl_2D_loc(iheartbelow,icarbon) .LT. zero) THEN
            ! If one of the biomass pools is negative, problems will
            ! occur later in the code.
            WRITE(numout,*) 'PFT, ',ivm
            WRITE(numout,*) '  bm_sapl_2D_loc(ispabove) ', &
                 bm_sapl_2D_loc(isapabove,icarbon)
            WRITE(numout,*) '  bm_sapl_2D_loc(isapbelow) ', &
                 bm_sapl_2D_loc(isapbelow,icarbon)
            WRITE(numout,*) '  bm_sapl_2D_loc(iheartabove) ', &
                 bm_sapl_2D_loc(iheartabove,icarbon)
            WRITE(numout,*) '  bm_sapl_2D_loc(iheartbelow) ', &
                 bm_sapl_2D_loc(iheartbelow,icarbon)
            CALL ipslerr_p(3, 'One of the biomass pools of a sapling is negative', &
                 'This will cause problems later', &
                 'fix it now in make_sapling','')
         ELSEIF (bm_sapl_2D_loc(isapabove,icarbon) .LT. min_stomate .OR. &
              bm_sapl_2D_loc(isapbelow,icarbon) .LT. min_stomate .OR. &
              bm_sapl_2D_loc(iheartabove,icarbon) .LT. min_stomate .OR. &
              bm_sapl_2D_loc(iheartbelow,icarbon) .LT. min_stomate) THEN
            ! One of the biomass pools is very small, recruitment will have
            ! little impact
            WRITE(numout,*) 'PFT, ',ivm
            WRITE(numout,*) '  bm_sapl_2D_loc(ispabove) ', &
                 bm_sapl_2D_loc(isapabove,icarbon)
            WRITE(numout,*) '  bm_sapl_2D_loc(isapbelow) ', &
                 bm_sapl_2D_loc(isapbelow,icarbon)
            WRITE(numout,*) '  bm_sapl_2D_loc(iheartabove) ', &
                 bm_sapl_2D_loc(iheartabove,icarbon)
            WRITE(numout,*) '  bm_sapl_2D_loc(iheartbelow) ', &
                 bm_sapl_2D_loc(iheartbelow,icarbon)
            CALL ipslerr_p(2, 'One of the biomass pools of a sapling is very small', &
                 'recruitment will have little impact', &
                 'check the parameter values','')
         END IF

      ENDIF  

      ! Use the allometric relationships to calculate initial leaf 
      ! and root mass. Note that wstress should be accounted for
      ! through c0_alloc. There was an inconsistency in the
      ! original calculation (OCN) - most pools were in gN but 
      ! leaves were in gC. A correction is proposed. Given that 
      ! icarbes was just set to zero the equation for ilabile is 
      ! a bit overkill and icarbes could be omitted.
      tree_height =  ((((bm_sapl_2D_loc(isapabove,icarbon)+ &
           bm_sapl_2D_loc(isapbelow,icarbon)+bm_sapl_2D_loc(iheartabove,icarbon)+ &
           bm_sapl_2D_loc(iheartbelow,icarbon))/(tree_ff(ivm)* &
           pipe_density(ivm))*4/pi)**(pipe_tune3(ivm)/2)) &
               *pipe_tune2_loc)**(1/(pipe_tune3(ivm)/2+1))

      bm_sapl_2D_loc(ileaf,icarbon) = ( bm_sapl_2D_loc(isapabove,icarbon) + &
           bm_sapl_2D_loc(isapbelow,icarbon) ) * KF_loc /  tree_height 
      bm_sapl_2D_loc(iroot,icarbon) = bm_sapl_2D_loc(ileaf,icarbon) / &
           ( KF_loc * c0_alloc )
      bm_sapl_2D_loc(icarbres,icarbon)=zero
      bm_sapl_2D_loc(ifruit,icarbon) = zero
      bm_sapl_2D_loc(ilabile,icarbon) =  labile_to_total * &
           (bm_sapl_2D_loc(ileaf,icarbon)  + &
           fcn_root(ivm) * bm_sapl_2D_loc(iroot,icarbon) + &
           fcn_wood(ivm) * (bm_sapl_2D_loc(isapabove,icarbon) + &
           bm_sapl_2D_loc(isapbelow,icarbon) + &
           bm_sapl_2D_loc(icarbres,icarbon)))      

      ! Debug
      IF (printlev_loc.GE.4 .AND. write_debug) THEN 
         WRITE(numout,*) ' +++CHECK INSIDE make_sapling+++'
         WRITE(numout,*) ' root to sapwood tradeoff p :', c0_alloc
         WRITE(numout,*) ' sapling height, sapling diameter ', &
              tree_height , ( tree_height/pipe_tune2_loc ) ** &
              ( 1. / pipe_tune3(ivm) )
         WRITE(numout,*) 'pipe_density, ',pipe_density(ivm)
         WRITE(numout,*) '++++++++++++++++++++++++++++++++'
      ENDIF
      ! -

      ! Allocate the nitrogen
      cn_leaf=cn_leaf_init_2D_loc
      cn_wood=cn_leaf/fcn_wood(ivm)
      cn_root=cn_leaf/fcn_root(ivm)

      ! Use the C/N ratio to calculate the N content for
      ! all the biomass components.
      bm_sapl_2D_loc(:,initrogen) = &
           bm_sapl_2D_loc(:,icarbon) / cn_root

      ! Overwrite the N content for the leaves and wood
      ! components.
      bm_sapl_2D_loc(ileaf,initrogen) = &
           bm_sapl_2D_loc(ileaf,icarbon) / cn_leaf
      bm_sapl_2D_loc(isapabove,initrogen) = &
           bm_sapl_2D_loc(isapabove,icarbon) / cn_wood
      bm_sapl_2D_loc(isapbelow,initrogen) = &
           bm_sapl_2D_loc(isapbelow,icarbon) / cn_wood
      bm_sapl_2D_loc(iheartabove,initrogen) = &
           bm_sapl_2D_loc(iheartabove,icarbon) / cn_wood
      bm_sapl_2D_loc(iheartbelow,initrogen) = &
           bm_sapl_2D_loc(iheartbelow,icarbon) / cn_wood

      ! Store the leaf mass so it can be used in the calculation of KF
      ! when using a dynamic sla.
      bypass_lm = bm_sapl_2D_loc(ileaf,icarbon)

      ! Avoid deciduous PFTs to have leaves out at establishment
      ! Whether the saplings have leaves or don't have leaves the
      ! first year doesn't really matter. Either the approach is
      ! correct in the northern hemisphere or in the southern
      ! hemisphere. Anyhow, a spin-up is needed to avoid issues
      ! with the initial conditions
      IF ( pheno_type(ivm) .NE. 1 ) THEN
         
         ! Not evergreen. Deciduous PFTs now need to survive an extra
         ! year before bud burst. To ensure survival there are several
         ! options: (a) either the height of the initial vegetation is
         ! increased (this results in more reserves) or (b) the reserves
         ! could be increased. The second option may result in numerical
         ! issues further down the code as the optimal reserve level is
         ! calculated from the other biomass pools. Also some initial
         ! tests showed that higher reserves simply resulted in more
         ! respiration.
         ! Tip: use taller trees to start with if they die between
         ! prescribe and leaf onset
         bm_sapl_2D_loc(icarbres,:) = (bm_sapl_2D_loc(icarbres,:) + &
              bm_sapl_2D_loc(ileaf,:) + bm_sapl_2D_loc(iroot,:))
         bm_sapl_2D_loc(ileaf,:) = zero
         bm_sapl_2D_loc(iroot,:) = zero
      
         ! When deciduous trees have no leaves they are dormant
         ! If the code is used for recruitment use the actual 
         ! plant_status of the overstorey trees for the saplings
         IF (establish) THEN
            plant_status_loc = idormant
         ENDIF
         
      ELSE
         
         ! Initilize carbohydrate reserves for evergreen PFTs
         bm_sapl_2D_loc(icarbres,:) = zero

         ! Evergreen trees always have a canopy. If the code
         ! is used for recruitment use the actual plant_status 
         ! of the overstorey trees for the saplings
         IF (establish) THEN
            plant_status_loc = icanopy
         ENDIF
    
      ENDIF

      ! Debug
      IF (printlev_loc>=4) THEN
         DO iele=1,nelements
            WRITE(numout,*) '  bm_sapl_2D_loc(ileaf) ',&
                 bm_sapl_2D_loc(ileaf,iele)
            WRITE(numout,*) '  bm_sapl_2D_loc(ispabove) ',&
                 bm_sapl_2D_loc(isapabove,iele)
            WRITE(numout,*) '  bm_sapl_2D_loc(isapbelow) ',&
                 bm_sapl_2D_loc(isapbelow,iele)
            WRITE(numout,*) '  bm_sapl_2D_loc(iheartabove) ',&
                 bm_sapl_2D_loc(iheartabove,iele)
            WRITE(numout,*) '  bm_sapl_2D_loc(iheartbelow) ',&
                 bm_sapl_2D_loc(iheartbelow,iele)
            WRITE(numout,*) '  bm_sapl_2D_loc(iroot) ',&
                 bm_sapl_2D_loc(iroot,iele)
            WRITE(numout,*) '  bm_sapl_2D_loc(ifruit)', &
                 bm_sapl_2D_loc(ifruit,iele)
            WRITE(numout,*) '  bm_sapl_2D_loc(icarbres)', &
                 bm_sapl_2D_loc(icarbres,iele)
            WRITE(numout,*) '  bm_sapl_2D_loc(ilabile)', &
                 bm_sapl_2D_loc(ilabile,iele)
         ENDDO
      ENDIF
      
   END SUBROUTINE make_saplings
  
END MODULE stomate_prescribe
