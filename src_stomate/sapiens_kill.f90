! =================================================================================================================================
! MODULE       : sapiens_kill
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       Kill all vegetation scheduled for killing by human activities.
!!
!!\n DESCRIPTION : Actually kills biomass.  The biomass is scheduled for killing
!!       by the variable circ_class_kill in other sapiens routines.  Here the
!!       biomass is moved to the litter and harvest pools, as appropriate.
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S) : 
!! - Sitch, S., B. Smith, et al. (2003), Evaluation of ecosystem dynamics,
!!         plant geography and terrestrial carbon cycling in the LPJ dynamic 
!!         global vegetation model, Global Change Biology, 9, 161-185.\n
!! - Waring, R. H. (1983). "Estimating forest growth and efficiency in relation 
!!         to canopy leaf area." Advances in Ecological Research 13: 327-354.\n
!! 
!! SVN          :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_stomate/sapiens_kill.f90 $
!! $Date: 2026-04-02 12:01:32 +0200 (jeu. 02 avril 2026) $
!! $Revision: 9457 $
!! \n
!_ ================================================================================================================================

MODULE sapiens_kill

  ! modules used:

  USE stomate_data
  USE grid
  USE pft_parameters
  USE ioipsl_para
  USE stomate_mark_kill,      ONLY : distribute_mortality_biomass
  USE function_library,       ONLY : nmax,  wood_to_dia, check_vegetation_area, &
                                     check_mass_balance, cc_to_biomass, get_printlev,&
                                     weibull_class_dist, wood_to_qmdia

  USE constantes

  IMPLICIT NONE

  ! private & public routines

  PRIVATE
  PUBLIC anthropogenic_mortality, anthropogenic_mortality_clear 
  
  ! Variable declaration

  LOGICAL, SAVE             :: firstcall_sapiens_kill = .TRUE.     !! first call flag
!$OMP THREADPRIVATE(firstcall_sapiens_kill)
   INTEGER(i_std), SAVE       :: printlev_loc                      !! Local level of text output for this module
!$OMP THREADPRIVATE(printlev_loc)


CONTAINS

!! ================================================================================================================================
!! SUBROUTINE   : anthropogenic_mortality_clear
!!
!>\BRIEF        Set the firstcall flag back to .TRUE. to prepare for the next simulation.
!_ ================================================================================================================================
  
  SUBROUTINE anthropogenic_mortality_clear
     firstcall_sapiens_kill = .TRUE.
  END SUBROUTINE anthropogenic_mortality_clear


!! ================================================================================================================================
!! SUBROUTINE   : anthropogenic_mortality
!!
!>\BRIEF        Transfer dead biomass to litter and update stand density as a result
!!              of human activities.
!!
!! DESCRIPTION  : This routine has a simple purpose: kill individuals by human activities.  
!!   The variable circ_class_kill indicates how many individuals need to be killed by the 
!!   various processes in the code, and is calculated in other modules.  anthropogenic_mortality 
!!   takes this variable and decides on the correct order of which to actually kill the
!!   trees.  As these are all human activities, some of the biomass may be moved to
!!   harvest pools, while other parts are moved to litter.
!!
!! RECENT CHANGE(S): 
!!
!! MAIN OUTPUT VARIABLE(S): ::circ_class_biomass; ::circ_class_n
!!
!! REFERENCE(S)   :
!! - Sitch, S., B. Smith, et al. (2003), Evaluation of ecosystem dynamics,
!!         plant geography and terrestrial carbon cycling in the LPJ dynamic 
!!         global vegetation model, Global Change Biology, 9, 161-185.
!! - Waring, R. H. (1983). "Estimating forest growth and efficiency in relation 
!!         to canopy leaf area." Advances in Ecological Research 13: 327-354.
!!
!! FLOWCHART    : None
!!\n
!_ ================================================================================================================================
 
  SUBROUTINE anthropogenic_mortality (npts,  &
       bm_to_litter, woody_litter_by_cut, circ_class_biomass, circ_class_kill, &
       circ_class_n, harvest_pool_bound, harvest_pool, harvest_type, &
       harvest_cut, harvest_area, &
       veget_max, age_stand, last_cut, mai_count, &
       coppice_dens, plant_status, biomass_cut)

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    ! NOTE: if harvest_pool_bound is not explicitely defined 
    ! the model uses 1:ndia_harvest+2 as its dimensions. If 
    ! you like to do so remember to change the counters in
    ! this routine
    INTEGER(i_std), INTENT(in)                                   :: npts               !! Domain size (-)
    REAL(r_std), DIMENSION(0:), INTENT(in)                       :: harvest_pool_bound !! The boundaries of the diameter classes

                                                                                       !! in the wood harvest pools
                                                                                       !! @tex $(m)$ @endtex
    REAL(r_std),DIMENSION(:,:),INTENT(in)                        :: veget_max          !! Maximum fraction of vegetation type 
    REAL(r_std), DIMENSION(:,:),INTENT(in)                       :: coppice_dens       !! The density of a coppice at the first
                                                                                       !! cutting.
                                                                                       !! @tex $( 1 / m**2 )$ @endtex 

    !! 0.2 Output variables
    REAL(r_std), DIMENSION(:,:,:,:,:,:), INTENT(out)             :: biomass_cut        !! Biomass change due to  ncut

    !! 0.3 Modified variables   
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)             :: circ_class_biomass !! Biomass of the components of the model  
                                                                                       !! tree within a circumference
                                                                                       !! class @tex $(gC ind^{-1})$ @endtex  
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)                 :: circ_class_n       !! Number of individuals in each circ class
                                                                                       !! @tex $(m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)             :: circ_class_kill    !! Number of trees within a circ that needs
                                                                                       !! to be killed @tex $(ind m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)               :: bm_to_litter       !! Biomass transfer to litter 
                                                                                       !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)                 :: woody_litter_by_cut!! Saved woody litter by icut not by iparts
                                                                                       !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)             :: harvest_pool       !! The wood and biomass that have been
                                                                                       !! havested by humans @tex $(gC)$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)                   :: harvest_type       !! Type of management that resulted
                                                                                       !! in the harvest (unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)                   :: harvest_cut        !! Type of cutting that was used for the harvest
                                                                                       !! (unitless)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)                 :: harvest_area       !! Harvested area (m^{2})
    INTEGER(i_std), DIMENSION(:,:),INTENT(inout)                 :: last_cut           !! Years since last thinning or clearcut (years)
    INTEGER(i_std), DIMENSION(:,:),INTENT(inout)                 :: age_stand          !! Age of stand (years)
    INTEGER(i_std), DIMENSION(:,:),INTENT(inout)                 :: mai_count          !! The number of times we've
                                                                                       !! calculated the volume increment
                                                                                       !! for a stand
    REAL(r_std), DIMENSION(:,:), INTENT(inout)                   :: plant_status       !! Growth and phenological status of the plant    

    !! 0.4 Local variables
    INTEGER(i_std)                                               :: ipts,ivm           !! Indices
    INTEGER(i_std)                                               :: iele, icir         !! Indices
    INTEGER(i_std)                                               :: ipar, imbc         !! Indices
    
    REAL(r_std), DIMENSION(npts,nvm,nmbcomp,nelements)           :: check_intern       !! Contains the components of the internal
                                                                                       !! mass balance chech for this routine
                                                                                       !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)                   :: closure_intern     !! Check closure of internal mass balance
                                                                                       !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)                   :: pool_start         !! Start pool of this routine 
                                                                                       !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)                   :: pool_end           !! End pool of this routine 
                                                                                       !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)                             :: veget_max_begin    !! temporary storage of veget_max to check area conservation
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements)            :: init_biomass       !! Temporary array of the initial biomass
    REAL(r_std), DIMENSION(nparts,nelements)                     :: actu_biomass       !! Temporary array of the actual biomass
!_ ================================================================================================================================

    IF (firstcall_sapiens_kill) THEN
       !! Initialize local printlev
      ! Guillaume M. -- get_printlev calls getin_p, a COLLECTIVE bcast. Called here,
      ! inside the compute loop under a local firstcall, it hangs the run as soon as the
      ! ranks stop reaching it together. Use printlev, read once at startup by every rank.
       printlev_loc=printlev
       firstcall_sapiens_kill=.FALSE.
    END IF

    IF (printlev_loc.GE.2) WRITE(numout,*) 'Entering anthropogenic_mortality'

 !! 1. Initialize variables

    ! 1.2 Initialize check for mass balance closure
    IF (err_act.GT.1) THEN

       check_intern(:,:,:,:) = zero
       pool_start = zero
       DO iele = 1,nelements
          DO ipar = 1,nparts
             DO icir = 1,ncirc
                !  Initial litter and biomass
                pool_start(:,:,iele) = pool_start(:,:,iele) + &
                     (circ_class_biomass(:,:,icir,ipar,iele) * &
                     circ_class_n(:,:,icir) * veget_max(:,:))
             ENDDO
             pool_start(:,:,iele) = pool_start(:,:,iele) + &
                  (bm_to_litter(:,:,ipar,iele) * veget_max(:,:))
          ENDDO

          ! The biomass harvest pool shouldn't be multiplied by 
          ! veget_max its units are already in gC pixel-1. Right 
          ! now the harvest pool is in gC.  We will divide by the
          ! area to make it per pixel. 
          DO ipts = 1, npts
             pool_start(ipts,:,iele) = pool_start(ipts,:,iele) + &
                  SUM(SUM(harvest_pool(ipts,:,:,iele,:),3),2) / &
                  (area(ipts) * contfrac(ipts))
          ENDDO
          
       ENDDO 

       !! 1.3 Initialize check for area conservation
       veget_max_begin(:,:) = veget_max(:,:)
       
    ENDIF ! err_act.GT.1

    ! The total biomass at the start of the routine
    init_biomass = cc_to_biomass(npts,nvm,circ_class_biomass(:,:,:,:,:),&
         circ_class_n(:,:,:))


 !! 3. Kill plants

    ! All trees marked to die from forest management or salavage 
    ! logging will be killed here (trees that died under a 
    ! conservation strategy are dealt with natural_mortality.  
    DO  ipts = 1,npts ! loop over land_points

       DO ivm = 2,nvm  ! loop over #PFT

          IF(veget_max(ipts,ivm) == zero)THEN
              ! this vegetation type is not present, so no reason to do the 
              ! calculation
              CYCLE
          ENDIF

          ! Clearcuts, salvage logging and thinning can only happen on trees.
          IF(is_tree(ivm))THEN
             
             ! Salvage logging after fire in thin/high stand
             ! (managed forest)
             CALL clearcut_harvest_aboveground_dead_roots(npts, ipts, &
                  ivm, ifm_thin, icut_fire, bm_to_litter, &
                  woody_litter_by_cut, circ_class_biomass, circ_class_kill, &
                  circ_class_n, harvest_pool, harvest_type, harvest_cut, &
                  harvest_area, harvest_pool_bound, veget_max, &
                  last_cut, age_stand, mai_count, &
                  max(zero, salvage_dist(ivm,ifm_thin) - &
                   branch_harvest(ivm,ifm_thin) * branch_ratio(ivm)))

             ! Store the biomass change of the stand due to fire
             ! Update initial biomass in case there will be more killing
             actu_biomass = cc_to_biomass(npts,nvm,&
                  circ_class_biomass(ipts,ivm,:,:,:),&
                  circ_class_n(ipts,ivm,:))
             biomass_cut(ipts,ivm,:,:,ifm_thin,icut_fire) = &
                  init_biomass(ipts,ivm,:,:) - actu_biomass(:,:)
             init_biomass(ipts,ivm,:,:) =  actu_biomass(:,:)
             
             ! Salvage logging after fire in uneven stand
             ! (managed forest)
             CALL clearcut_harvest_aboveground_dead_roots(npts, ipts, &
                  ivm, ifm_uneven, icut_fire, bm_to_litter, &
                  woody_litter_by_cut, circ_class_biomass, circ_class_kill, &
                  circ_class_n, harvest_pool, harvest_type, harvest_cut, &
                  harvest_area, harvest_pool_bound, veget_max, &
                  last_cut, age_stand, mai_count, &
                  max(zero, salvage_dist(ivm,ifm_uneven) - &
                   branch_harvest(ivm,ifm_uneven) * branch_ratio(ivm)))
             
             ! Store the biomass change of the stand due to fire
             ! Update initial biomass in case there will be more killing
             actu_biomass = cc_to_biomass(npts,nvm,&
                  circ_class_biomass(ipts,ivm,:,:,:),&
                  circ_class_n(ipts,ivm,:))
             biomass_cut(ipts,ivm,:,:,ifm_uneven,icut_fire) = &
                  init_biomass(ipts,ivm,:,:) - actu_biomass(:,:)
             init_biomass(ipts,ivm,:,:) =  actu_biomass(:,:)

             ! Salvage logging after fire in coppice stand (managed
             ! forest )                
             CALL clearcut_harvest_aboveground_dead_roots(npts, ipts, &
                  ivm, ifm_cop, icut_fire, bm_to_litter, &
                  woody_litter_by_cut, circ_class_biomass, circ_class_kill, &
                  circ_class_n, harvest_pool, harvest_type, harvest_cut, &
                  harvest_area, harvest_pool_bound, veget_max, &
                  last_cut, age_stand, mai_count,&
                  max(zero, salvage_dist(ivm,ifm_cop) - &
                   branch_harvest(ivm,ifm_cop) * branch_ratio(ivm)))

             ! Store the biomass change of the stand due to fire
             ! Update initial biomass in case there will be more killing
             actu_biomass = cc_to_biomass(npts,nvm,&
                  circ_class_biomass(ipts,ivm,:,:,:),&
                  circ_class_n(ipts,ivm,:))
             biomass_cut(ipts,ivm,:,:,ifm_cop,icut_fire) = &
                  init_biomass(ipts,ivm,:,:) - actu_biomass(:,:)
             init_biomass(ipts,ivm,:,:) =  actu_biomass(:,:)

             ! Salvage logging after fire in SRC (Short Rotation Coppice) 
             ! stand (managed forest )
             CALL clearcut_harvest_aboveground_dead_roots(npts, ipts, &
                  ivm, ifm_src, icut_fire, bm_to_litter, &
                  woody_litter_by_cut, circ_class_biomass, circ_class_kill, &
                  circ_class_n, harvest_pool, harvest_type, harvest_cut, &
                  harvest_area, harvest_pool_bound, veget_max, &
                  last_cut, age_stand, mai_count,&
                  max(zero, salvage_dist(ivm,ifm_src) - &
                   branch_harvest(ivm,ifm_src) * branch_ratio(ivm)))

             ! Store the biomass change of the stand due to fire
             ! Update initial biomass in case there will be more killing
             actu_biomass = cc_to_biomass(npts,nvm,&
                  circ_class_biomass(ipts,ivm,:,:,:),&
                  circ_class_n(ipts,ivm,:))
             biomass_cut(ipts,ivm,:,:,ifm_src,icut_fire) = &
                  init_biomass(ipts,ivm,:,:) - actu_biomass(:,:)
             init_biomass(ipts,ivm,:,:) =  actu_biomass(:,:)


             ! Salvage logging after bark beetle attack in thin/high stand (managed
             ! forest)
             CALL clearcut_harvest_aboveground_dead_roots(npts, ipts, &
                  ivm, ifm_thin, icut_beetle, bm_to_litter, &
                  woody_litter_by_cut, circ_class_biomass, circ_class_kill, &
                  circ_class_n, harvest_pool, harvest_type, harvest_cut, &
                  harvest_area, harvest_pool_bound, veget_max, &
                  last_cut, age_stand, mai_count, &
                  max(zero, salvage_dist(ivm,ifm_thin) - &
                   branch_harvest(ivm,ifm_thin) * branch_ratio(ivm)))

             ! Store the biomass change of the stand due to bark beetle
             ! Update initial biomass in case there will be more killing
             actu_biomass = cc_to_biomass(npts,nvm,&
                  circ_class_biomass(ipts,ivm,:,:,:),&
                  circ_class_n(ipts,ivm,:))
             biomass_cut(ipts,ivm,:,:,ifm_thin,icut_beetle) = &
                  init_biomass(ipts,ivm,:,:) - actu_biomass(:,:)
             init_biomass(ipts,ivm,:,:) =  actu_biomass(:,:)

             CALL clearcut_harvest_aboveground_dead_roots(npts, ipts, &
                  ivm, ifm_uneven, icut_beetle, bm_to_litter, &
                  woody_litter_by_cut, circ_class_biomass, circ_class_kill, &
                  circ_class_n, harvest_pool, harvest_type, harvest_cut, &
                  harvest_area, harvest_pool_bound, veget_max, &
                  last_cut, age_stand, mai_count, &
                  max(zero, salvage_dist(ivm,ifm_uneven) - &
                   branch_harvest(ivm,ifm_uneven) * branch_ratio(ivm)))

             ! Store the biomass change of the stand due to bark beetle
             ! Update initial biomass in case there will be more killing
             actu_biomass = cc_to_biomass(npts,nvm,&
                  circ_class_biomass(ipts,ivm,:,:,:),&
                  circ_class_n(ipts,ivm,:))
             biomass_cut(ipts,ivm,:,:,ifm_uneven,icut_beetle) = &
                  init_biomass(ipts,ivm,:,:) - actu_biomass(:,:)
             init_biomass(ipts,ivm,:,:) =  actu_biomass(:,:)

             ! Salvage logging after windthrow  in thin/high stand (managed
             ! forest)
             CALL clearcut_harvest_aboveground_dead_roots(npts, ipts, &
                  ivm, ifm_thin, icut_storm_break, bm_to_litter, &
                  woody_litter_by_cut, circ_class_biomass, circ_class_kill, &
                  circ_class_n, harvest_pool, harvest_type, harvest_cut, &
                  harvest_area, harvest_pool_bound, veget_max, &
                  last_cut, age_stand, mai_count, &
                  max(zero, salvage_dist(ivm,ifm_thin) - &
                   branch_harvest(ivm,ifm_thin) * branch_ratio(ivm)))

             ! Store the biomass change of the stand due to storm damage
             ! Update initial biomass in case there will be more killing
             actu_biomass = cc_to_biomass(npts,nvm,&
                  circ_class_biomass(ipts,ivm,:,:,:),&
                  circ_class_n(ipts,ivm,:))
             biomass_cut(ipts,ivm,:,:,ifm_thin,icut_storm_break) = &
                  init_biomass(ipts,ivm,:,:) - actu_biomass(:,:)
             init_biomass(ipts,ivm,:,:) =  actu_biomass(:,:)
             
             ! Salvage logging after windthrow  in thin/high stand (managed
             ! forest)
             CALL clearcut_harvest_aboveground_dead_roots(npts, ipts, &
                  ivm, ifm_uneven, icut_storm_break, bm_to_litter, &
                  woody_litter_by_cut, circ_class_biomass, circ_class_kill, &
                  circ_class_n, harvest_pool, harvest_type, harvest_cut, &
                  harvest_area, harvest_pool_bound, veget_max, &
                  last_cut, age_stand, mai_count, &
                  max(zero, salvage_dist(ivm,ifm_uneven) - &
                   branch_harvest(ivm,ifm_uneven) * branch_ratio(ivm)))

             ! Store the biomass change of the stand due to storm damage
             ! Update initial biomass in case there will be more killing
             actu_biomass = cc_to_biomass(npts,nvm,&
                  circ_class_biomass(ipts,ivm,:,:,:),&
                  circ_class_n(ipts,ivm,:))
             biomass_cut(ipts,ivm,:,:,ifm_uneven,icut_storm_break) = &
                  init_biomass(ipts,ivm,:,:) - actu_biomass(:,:)
             init_biomass(ipts,ivm,:,:) =  actu_biomass(:,:)


             ! Salvage logging after windthrow in coppice stand (managed
             ! forest )
             CALL clearcut_harvest_aboveground_dead_roots(npts, ipts, &
                  ivm, ifm_cop, icut_storm_break, bm_to_litter, &
                  woody_litter_by_cut, circ_class_biomass, circ_class_kill, &
                  circ_class_n, harvest_pool, harvest_type, harvest_cut, &
                  harvest_area, harvest_pool_bound, veget_max, &
                  last_cut, age_stand, mai_count,&
                  max(zero, salvage_dist(ivm,ifm_cop) - &
                   branch_harvest(ivm,ifm_cop) * branch_ratio(ivm)))

             ! Store the biomass change of the stand due to storm damage
             ! Update initial biomass in case there will be more killing
             actu_biomass = cc_to_biomass(npts,nvm,&
                  circ_class_biomass(ipts,ivm,:,:,:),&
                  circ_class_n(ipts,ivm,:))
             biomass_cut(ipts,ivm,:,:,ifm_cop,icut_storm_break) = &
                  init_biomass(ipts,ivm,:,:) - actu_biomass(:,:)
             init_biomass(ipts,ivm,:,:) =  actu_biomass(:,:)
             
             ! Salvage logging after windthrow in SRC (Short Rotation Coppice) stand (managed
             ! forest )
             CALL clearcut_harvest_aboveground_dead_roots(npts, ipts, &
                  ivm, ifm_src, icut_storm_break, bm_to_litter, &
                  woody_litter_by_cut, circ_class_biomass, circ_class_kill, &
                  circ_class_n, harvest_pool, harvest_type, harvest_cut, &
                  harvest_area, harvest_pool_bound, veget_max, &
                  last_cut, age_stand, mai_count,&
                  max(zero, salvage_dist(ivm,ifm_src) - &
                   branch_harvest(ivm,ifm_src) * branch_ratio(ivm)))

             ! Store the biomass change of the stand due to storm damage
             ! Update initial biomass in case there will be more killing
             actu_biomass = cc_to_biomass(npts,nvm,&
                  circ_class_biomass(ipts,ivm,:,:,:),&
                  circ_class_n(ipts,ivm,:))
             biomass_cut(ipts,ivm,:,:,ifm_src,icut_storm_break) = &
                  init_biomass(ipts,ivm,:,:) - actu_biomass(:,:)
             init_biomass(ipts,ivm,:,:) =  actu_biomass(:,:)

             ! Salvage logging after windthrow  in thin/high stand (managed
             ! forest) 
             CALL clearcut_harvest_aboveground_dead_roots(npts, ipts, &
                  ivm, ifm_thin, icut_storm_uproot, bm_to_litter, &
                  woody_litter_by_cut, circ_class_biomass, circ_class_kill, &
                  circ_class_n, harvest_pool, harvest_type, harvest_cut, &
                  harvest_area, harvest_pool_bound, veget_max, &
                  last_cut, age_stand, mai_count, &
                  max(zero, salvage_dist(ivm,ifm_thin) - &
                   branch_harvest(ivm,ifm_thin) * branch_ratio(ivm)))

             ! Store the biomass change of the stand due to storm damage
             ! Update initial biomass in case there will be more killing
             actu_biomass = cc_to_biomass(npts,nvm,&
                  circ_class_biomass(ipts,ivm,:,:,:),&
                  circ_class_n(ipts,ivm,:))
             biomass_cut(ipts,ivm,:,:,ifm_thin,icut_storm_uproot) = &
                  init_biomass(ipts,ivm,:,:) - actu_biomass(:,:)
             init_biomass(ipts,ivm,:,:) =  actu_biomass(:,:)
             
             ! Salvage logging after windthrow  in thin/high stand (managed
             ! forest) 
             CALL clearcut_harvest_aboveground_dead_roots(npts, ipts, &
                  ivm, ifm_uneven, icut_storm_uproot, bm_to_litter, &
                  woody_litter_by_cut, circ_class_biomass, circ_class_kill, &
                  circ_class_n, harvest_pool, harvest_type, harvest_cut, &
                  harvest_area, harvest_pool_bound, veget_max, &
                  last_cut, age_stand, mai_count, &
                  max(zero, salvage_dist(ivm,ifm_uneven) - &
                   branch_harvest(ivm,ifm_uneven) * branch_ratio(ivm)))

             ! Store the biomass change of the stand due to storm damage
             ! Update initial biomass in case there will be more killing
             actu_biomass = cc_to_biomass(npts,nvm,&
                  circ_class_biomass(ipts,ivm,:,:,:),&
                  circ_class_n(ipts,ivm,:))
             biomass_cut(ipts,ivm,:,:,ifm_uneven,icut_storm_uproot) = &
                  init_biomass(ipts,ivm,:,:) - actu_biomass(:,:)
             init_biomass(ipts,ivm,:,:) =  actu_biomass(:,:)

             ! Salvage logging after windthrow in coppice stand (managed
             ! forest )                
             CALL clearcut_harvest_aboveground_dead_roots(npts, ipts, &
                  ivm, ifm_cop, icut_storm_uproot, bm_to_litter, &
                  woody_litter_by_cut, circ_class_biomass, circ_class_kill, &
                  circ_class_n, harvest_pool, harvest_type, harvest_cut, &
                  harvest_area, harvest_pool_bound, veget_max, &
                  last_cut, age_stand, mai_count,&
                  max(zero, salvage_dist(ivm,ifm_cop) - &
                   branch_harvest(ivm,ifm_cop) * branch_ratio(ivm)))

             ! Store the biomass change of the stand due to storm damage
             ! Update initial biomass in case there will be more killing
             actu_biomass = cc_to_biomass(npts,nvm,&
                  circ_class_biomass(ipts,ivm,:,:,:),&
                  circ_class_n(ipts,ivm,:))
             biomass_cut(ipts,ivm,:,:,ifm_cop,icut_storm_uproot) = &
                  init_biomass(ipts,ivm,:,:) - actu_biomass(:,:)
             init_biomass(ipts,ivm,:,:) =  actu_biomass(:,:)

             ! Salvage logging after windthrow in SRC (Short Rotation Coppice) 
             ! stand (managed forest )
             CALL clearcut_harvest_aboveground_dead_roots(npts, ipts, &
                  ivm, ifm_src, icut_storm_uproot, bm_to_litter, &
                  woody_litter_by_cut, circ_class_biomass, circ_class_kill, &
                  circ_class_n, harvest_pool, harvest_type, harvest_cut, &
                  harvest_area, harvest_pool_bound, veget_max, &
                  last_cut, age_stand, mai_count,&
                  max(zero, salvage_dist(ivm,ifm_src) - &
                   branch_harvest(ivm,ifm_src) * branch_ratio(ivm)))

             ! Store the biomass change of the stand due to storm damage
             ! Update initial biomass in case there will be more killing
             actu_biomass = cc_to_biomass(npts,nvm,&
                  circ_class_biomass(ipts,ivm,:,:,:),&
                  circ_class_n(ipts,ivm,:))
             biomass_cut(ipts,ivm,:,:,ifm_src,icut_storm_uproot) = &
                  init_biomass(ipts,ivm,:,:) - actu_biomass(:,:)
             init_biomass(ipts,ivm,:,:) =  actu_biomass(:,:)
           
             ! Now see if there are any plants to be killed by a clearcut.
             ! Only the stems are harvested. This should be done before
             ! a thinning so that we don't overcount.
             CALL clearcut_harvest_aboveground_dead_roots(npts, ipts, &
                  ivm, ifm_thin, icut_clear, bm_to_litter, &
                  woody_litter_by_cut, circ_class_biomass, circ_class_kill, &
                  circ_class_n, harvest_pool, harvest_type, harvest_cut, &
                  harvest_area, harvest_pool_bound, veget_max,&
                  last_cut, age_stand, mai_count,&
                  max(zero, clearcut_residue(ivm,ifm_thin) - &
                   branch_harvest(ivm,ifm_thin) * branch_ratio(ivm)))

             ! Store the biomass change of the stand due to storm damage
             ! Update initial biomass in case there will be more killing
             actu_biomass = cc_to_biomass(npts,nvm,&
                  circ_class_biomass(ipts,ivm,:,:,:),&
                  circ_class_n(ipts,ivm,:))
             biomass_cut(ipts,ivm,:,:,ifm_thin,icut_clear) = &
                  init_biomass(ipts,ivm,:,:) - actu_biomass(:,:)
             init_biomass(ipts,ivm,:,:) =  actu_biomass(:,:)

             ! Not only do we harvest stems during thinning, we have to 
             ! adjust the mortality.  We'll do all this in a subroutine.
             CALL thinning_harvest_aboveground_dead_roots(npts, ipts, &
                  ivm, ifm_thin, icut_thin, bm_to_litter, &
                  woody_litter_by_cut, circ_class_biomass, circ_class_kill, &
                  circ_class_n, harvest_pool, harvest_type, harvest_cut, &
                  harvest_area, harvest_pool_bound, veget_max, &
                  last_cut, max(zero, clearcut_residue(ivm,ifm_thin) - &
                   branch_harvest(ivm,ifm_thin) * branch_ratio(ivm)))

             ! Store the biomass change of the stand due to storm damage
             ! Update initial biomass in case there will be more killing
             actu_biomass = cc_to_biomass(npts,nvm,&
                  circ_class_biomass(ipts,ivm,:,:,:),&
                  circ_class_n(ipts,ivm,:))
             biomass_cut(ipts,ivm,:,:,ifm_thin,icut_thin) = &
                  init_biomass(ipts,ivm,:,:) - actu_biomass(:,:)
             init_biomass(ipts,ivm,:,:) =  actu_biomass(:,:)

             ! Now see if there are any plants to be killed by a clearcut.
             ! Only the stems are harvested. This should be done before
             ! a thinning so that we don't overcount. 
             ! We should never clearcut in an uneven forest. Remove?
             CALL clearcut_harvest_aboveground_dead_roots(npts, ipts, &
                  ivm, ifm_uneven, icut_clear, bm_to_litter, &
                  woody_litter_by_cut, circ_class_biomass, circ_class_kill, &
                  circ_class_n, harvest_pool, harvest_type, harvest_cut, &
                  harvest_area, harvest_pool_bound, veget_max,&
                  last_cut, age_stand, mai_count,&
                  max(zero, clearcut_residue(ivm,ifm_uneven) - &
                   branch_harvest(ivm,ifm_uneven) * branch_ratio(ivm)))

             ! Store the biomass change of the stand due to storm damage
             ! Update initial biomass in case there will be more killing
             actu_biomass = cc_to_biomass(npts,nvm,&
                  circ_class_biomass(ipts,ivm,:,:,:),&
                  circ_class_n(ipts,ivm,:))
             biomass_cut(ipts,ivm,:,:,ifm_uneven,icut_clear) = &
                  init_biomass(ipts,ivm,:,:) - actu_biomass(:,:)
             init_biomass(ipts,ivm,:,:) =  actu_biomass(:,:)

             ! Not only do we harvest stems during thinning, we have to 
             ! adjust the mortality.  We'll do all this in a subroutine.
             CALL thinning_harvest_aboveground_dead_roots(npts, ipts, &
                  ivm, ifm_uneven, icut_thin, bm_to_litter, &
                  woody_litter_by_cut, circ_class_biomass, circ_class_kill, &
                  circ_class_n, harvest_pool, harvest_type, harvest_cut, &
                  harvest_area, harvest_pool_bound, veget_max, &
                  last_cut, max(zero, clearcut_residue(ivm,ifm_uneven) - &
                    branch_harvest(ivm,ifm_uneven) * branch_ratio(ivm)))

             ! Store the biomass change of the stand due to storm damage
             ! Update initial biomass in case there will be more killing
             actu_biomass = cc_to_biomass(npts,nvm,&
                  circ_class_biomass(ipts,ivm,:,:,:),&
                  circ_class_n(ipts,ivm,:))
             biomass_cut(ipts,ivm,:,:,ifm_uneven,icut_thin) = &
                  init_biomass(ipts,ivm,:,:) - actu_biomass(:,:)
             init_biomass(ipts,ivm,:,:) =  actu_biomass(:,:)


             ! And the final cut of an SRC, where the underground wood 
             ! dies. This is exactly the same as the clearcut above, 
             ! except that we take all of the aboveground biomass.
             CALL clearcut_harvest_aboveground_dead_roots(npts, ipts, &
                  ivm, ifm_src, icut_cop3, bm_to_litter, &
                  woody_litter_by_cut, circ_class_biomass, circ_class_kill, &
                  circ_class_n, harvest_pool, harvest_type, harvest_cut, &
                  harvest_area, harvest_pool_bound, veget_max,&
                  last_cut, age_stand, mai_count, &
                  max(zero, clearcut_residue(ivm,ifm_src) - &
                   branch_harvest(ivm,ifm_src) * branch_ratio(ivm)))

             ! Store the biomass change of the stand due to storm damage
             ! Update initial biomass in case there will be more killing
             actu_biomass = cc_to_biomass(npts,nvm,&
                  circ_class_biomass(ipts,ivm,:,:,:),&
                  circ_class_n(ipts,ivm,:))
             biomass_cut(ipts,ivm,:,:,ifm_src,icut_cop3) = &
                  init_biomass(ipts,ivm,:,:) - actu_biomass(:,:)
             init_biomass(ipts,ivm,:,:) =  actu_biomass(:,:)

             ! Have we cut a coppice system for the first time?
             CALL first_coppice_cut(npts, ipts, ivm, ifm_cop, &
                  icut_cop1, harvest_pool_bound, veget_max,&
                  circ_class_biomass, circ_class_n, &
                  circ_class_kill, bm_to_litter, harvest_pool, &
                  harvest_type, harvest_cut, harvest_area, last_cut)

             ! Store the biomass change of the stand due to storm damage
             ! Update initial biomass in case there will be more killing
             actu_biomass = cc_to_biomass(npts,nvm,&
                  circ_class_biomass(ipts,ivm,:,:,:),&
                  circ_class_n(ipts,ivm,:))
             biomass_cut(ipts,ivm,:,:,ifm_cop,icut_cop1) = &
                  init_biomass(ipts,ivm,:,:) - actu_biomass(:,:)
             init_biomass(ipts,ivm,:,:) =  actu_biomass(:,:)

             ! Have we cut an SRC system for the first time?
             CALL first_coppice_cut(npts, ipts, ivm, ifm_src, &
                  icut_cop1, harvest_pool_bound, veget_max,&
                  circ_class_biomass, circ_class_n, &
                  circ_class_kill, bm_to_litter, harvest_pool, &
                  harvest_type, harvest_cut, harvest_area, last_cut)

             ! Store the biomass change of the stand due to storm damage
             ! Update initial biomass in case there will be more killing
             actu_biomass = cc_to_biomass(npts,nvm,&
                  circ_class_biomass(ipts,ivm,:,:,:),&
                  circ_class_n(ipts,ivm,:))
             biomass_cut(ipts,ivm,:,:,ifm_src,icut_cop1) = &
                  init_biomass(ipts,ivm,:,:) - actu_biomass(:,:)
             init_biomass(ipts,ivm,:,:) =  actu_biomass(:,:)

             ! What if we coppice for the second (and subsequent) cuts?
             CALL second_coppice_cut(npts, ipts, ivm, ifm_cop, &
                  icut_cop2, harvest_pool_bound, veget_max,&
                  circ_class_biomass, circ_class_n, &
                  circ_class_kill, bm_to_litter, harvest_pool, &
                  harvest_type, harvest_cut, harvest_area, &
                  last_cut, coppice_dens)
             
             ! Store the biomass change of the stand due to storm damage
             ! Update initial biomass in case there will be more killing
             actu_biomass = cc_to_biomass(npts,nvm,&
                  circ_class_biomass(ipts,ivm,:,:,:),&
                  circ_class_n(ipts,ivm,:))
             biomass_cut(ipts,ivm,:,:,ifm_cop,icut_cop2) = &
                  init_biomass(ipts,ivm,:,:) - actu_biomass(:,:)
             init_biomass(ipts,ivm,:,:) =  actu_biomass(:,:)

             ! For short rotation coppices
             CALL second_coppice_cut(npts, ipts, ivm, ifm_src, &
                  icut_cop2, harvest_pool_bound, veget_max,&
                  circ_class_biomass, circ_class_n, &
                  circ_class_kill, bm_to_litter, harvest_pool, &
                  harvest_type, harvest_cut, harvest_area, &
                  last_cut, coppice_dens)

             ! Store the biomass change of the stand due to storm damage
             ! Update initial biomass in case there will be more killing
             actu_biomass = cc_to_biomass(npts,nvm,&
                  circ_class_biomass(ipts,ivm,:,:,:),&
                  circ_class_n(ipts,ivm,:))
             biomass_cut(ipts,ivm,:,:,ifm_src,icut_cop2) = &
                  init_biomass(ipts,ivm,:,:) - actu_biomass(:,:)
             init_biomass(ipts,ivm,:,:) =  actu_biomass(:,:)
             
          END IF ! is_tree

          ! Here some circ classes may be empty.  We do not want to do
          ! anything about that yet, since we have not killed the natural
          ! mortality biomass.  If we redistribute the biomass now we will
          ! change the circ_class_biomass, but the old values for this
          ! were used to determine circ_class_kill(inatural) in 
          ! stomate_mark_kill.  If we change those values we will kill 
          ! the wrong amount of biomass in stomate_kill.
          !
          ! It is possible that all the biomass was killed. The only way this
          ! can been seen is because circ_class_n is zero. circ_class_biomass is
          ! the biomass of an individual tree and is not affected by the above
          ! mortality unless all trees are killed. In the latter case plant_status 
          ! should then be set to idead. In mortality_clean all the necessary 
          ! actions will be taken to make sure that the PFT will be prescribed 
          ! at the start of the year.
          IF (SUM(circ_class_n(ipts,ivm,:)).LT.min_stomate)THEN

             ! Kill the PFT
             plant_status(ipts,ivm) = idead
  
            ! Move the remaining biomass into the proper
            ! litter pools. Note that these lines help to ensure
            ! mass balance closure because circ_class_n in this
            ! multiplication is less than min_stomate!
            DO ipar = 1,nparts
               DO iele = 1,nelements
                     bm_to_litter(ipts,ivm,ipar,iele) = &
                          bm_to_litter(ipts,ivm,ipar,iele) + &
                          SUM(circ_class_biomass(ipts,ivm,:,ipar,iele) * &
                          circ_class_n(ipts,ivm,:))
               ENDDO  
            ENDDO

            ! Set the biomass and the number of individuals to zero
            circ_class_biomass(ipts,ivm,:,:,:) = zero
            circ_class_n(ipts,ivm,:) = zero

          ENDIF ! 

       ENDDO  ! loop over pfts

    END DO ! loop over land points

 !! 4. Check numerical consistency of this routine

    IF (err_act.GT.1) THEN

       ! 4.2 Check surface area
       CALL check_vegetation_area("stomate_anthropogenic_mortality", npts, &
            veget_max_begin, veget_max,'pft')

       ! 4.3 Mass balance closure 
       pool_end = zero
       DO iele = 1,nelements
          DO ipar = 1,nparts
             DO icir = 1,ncirc
                pool_end(:,:,iele) = pool_end(:,:,iele) + &
                     (circ_class_biomass(:,:,icir,ipar,iele) * &
                     circ_class_n(:,:,icir) * veget_max(:,:)) 
             ENDDO
             pool_end(:,:,iele) = pool_end(:,:,iele) + &
                     bm_to_litter(:,:,ipar,iele) * veget_max(:,:)
          ENDDO

          ! The biomass harvest pool shouldn't be multiplied by veget_max
          ! its units are already in gC pixel-1. Right now the harvest 
          ! pool is in gC.  We will divide by the area to make it 
          ! per pixel.
          DO ipts = 1, npts
             pool_end(ipts,:,iele) = pool_end(ipts,:,iele) + &
                  SUM(SUM(harvest_pool(ipts,:,:,iele,:),3),2) / &
                  (area(ipts) * contfrac(ipts))
          ENDDO
       ENDDO ! nelements
   
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
             IF (printlev_loc>=4) WRITE(numout,*) &
                  'check_intern, ivm, imbc, iele, ', imbc, &
                  iele, check_intern(:,test_pft,imbc,iele)
             !-
             closure_intern(:,:,iele) = closure_intern(:,:,iele) + &
                  check_intern(:,:,imbc,iele)
          ENDDO
       ENDDO

       ! 4.3.3 Check mass balance closure
       CALL check_mass_balance("stomate_anthropogenic_mortality", closure_intern, &
            npts, pool_end, pool_start, veget_max, 'pft')
      
    ENDIF ! err_act.GT.1

    IF (printlev.GE.4) WRITE(numout,*) 'Leaving anthropogenic_mortality'

  END SUBROUTINE anthropogenic_mortality


!! ================================================================================================================================
!! SUBROUTINE   : thinning_harvest_aboveground_dead_roots
!!
!>\BRIEF        Transfer dead biomass to litter and update stand density due
!!              to thinning operations
!!
!! DESCRIPTION  : Moves dead biomass to the litter pool and the harvest pool
!!                where appropriate.  Since not all the trees are killed,
!!                fewer counters have to be reset here than in the clearcut
!!                case.  In addition, thinning will change the fraction of
!!                trees killed by mortality, so we account for that as well.
!!
!!                NOTE: If harvest_fraction is 1, all the aboveground biomass is
!!                      harvested.  If it's 1-branch_ratio, only stems are 
!!                      harvested.  Belowground wood, roots, fruits, and
!!                      leaves all become litter.  A portion of aboveground
!!                      wood is put in the harvest pool.  Reserves and labile
!!                      carbon are assumed to be distributed evenly across
!!                      all non-heartwood components.
!!
!! RECENT CHANGE(S): 
!!
!! MAIN OUTPUT VARIABLE(S): ::circ_class_biomass; ::circ_class_n; ::circ_class_kill
!!
!! REFERENCE(S)   :
!!
!! FLOWCHART    : None
!!\n
!_ ================================================================================================================================
 
  SUBROUTINE thinning_harvest_aboveground_dead_roots (npts, ipts, &
       ivm, ifm, icut, bm_to_litter, woody_litter_by_cut, &
       circ_class_biomass, circ_class_kill, circ_class_n, &
       harvest_pool, harvest_type, harvest_cut, harvest_area, &
       harvest_pool_bound, veget_max, last_cut, &
       harvest_fraction)


    !! 0. Variable and parameter declaration

    !! 0.1 Input variables 
    ! NOTE: if harvest_pool_bound is not explicitely defined 
    ! the model uses 1:ndia_harvest+2 as its dimensions. If 
    ! you like to do so remember to change the counters in
    ! this routine
    INTEGER(i_std), INTENT(in)                                   :: npts               !! Domain size (-)
    INTEGER(i_std), INTENT(in)                                   :: ipts               !! current grid square
    INTEGER(i_std), INTENT(in)                                   :: ivm                !! current PFT
    INTEGER(i_std), INTENT(in)                                   :: ifm                !! current fm strategy
    INTEGER(i_std), INTENT(in)                                   :: icut               !! current cut time
    REAL(r_std), DIMENSION(0:), INTENT(in)                       :: harvest_pool_bound !! The boundaries of the diameter classes
                                                                                       !! in the wood harvest pools
                                                                                       !! @tex $(m)$ @endtex
    REAL(r_std),DIMENSION(:,:),INTENT(in)                        :: veget_max          !! Maximum fraction of vegetation type 
    REAL(r_std),INTENT(in)                                       :: harvest_fraction   !! The fraction of the aboveground
                                                                                       !! biomass that is harvested.  This is
                                                                                       !! generally 1-branch_ratio if you only
                                                                                       !! want to harvest the stems, and 1 if
                                                                                       !! you want to harvest everything.

    !! 0.2 Output variables


    !! 0.3 Modified variables   
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)             :: circ_class_biomass !! Biomass of the components of the model  
                                                                                       !! tree within a circumference
                                                                                       !! class @tex $(gC ind^{-1})$ @endtex  
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)                 :: circ_class_n       !! Number of individuals in each circ class
                                                                                       !! @tex $(m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)             :: circ_class_kill    !! Number of trees within a circ that needs
                                                                                       !! to be killed @tex $(ind m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)               :: bm_to_litter       !! Biomass transfer to litter 
                                                                                       !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)                 :: woody_litter_by_cut!! Saved woody litter by icut not by iparts
                                                                                       !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)             :: harvest_pool       !! The wood and biomass that have been
                                                                                       !! havested by humans @tex $(gC)$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)                   :: harvest_type       !! Type of management that resulted
                                                                                       !! in the harvest (unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)                   :: harvest_cut        !! Type of cutting that was used for the harvest
                                                                                       !! (unitless)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)                 :: harvest_area       !! Harvested area (m^{2})
    INTEGER(i_std), DIMENSION(:,:),INTENT(inout)                 :: last_cut           !! Years since last thinning or clearcut (years)

    !! 0.4 Local variables

    INTEGER(i_std)                                               :: icir, imbc, iele   !! Indices
    INTEGER(i_std)                                               :: ipar               !! Indices
    REAL(r_std)                                                  :: total_bm
    REAL(r_std)                                                  :: thinned_bm
    REAL(r_std)                                                  :: dead_bm
    REAL(r_std)                                                  :: bm_difference
    REAL(r_std), DIMENSION(npts,nvm,nmbcomp,nelements)           :: check_intern       !! Contains the components of the internal
                                                                                       !! mass balance chech for this routine
                                                                                       !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)                   :: closure_intern     !! Check closure of internal mass balance
                                                                                       !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)                   :: pool_start         !! Start pool of this routine 
                                                                                       !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)                   :: pool_end           !! End pool of this routine 
                                                                                       !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    !_ ================================================================================================================================

    IF (printlev.GE.3) WRITE(numout,*) &
         'Entering thinning_harvest_aboveground_dead_roots sapiens_kill.f90'

    ! 1.2 Initialize check for mass balance closure
    IF (err_act.GT.1) THEN
       
       check_intern(:,:,:,:) = zero
       pool_start = zero
       DO iele = 1,nelements
          DO ipar = 1,nparts
             ! Initial litter and biomass. This routine is called
             ! within an ipts-LOOP and an ivm-LOOP so we test one 
             ! pixel and one pft at a time. Store in the first position.
             DO icir = 1,ncirc
                pool_start(1,1,iele) = pool_start(1,1,iele) + &
                     (circ_class_biomass(ipts,ivm,icir,ipar,iele) * &
                     circ_class_n(ipts,ivm,icir) * veget_max(ipts,ivm))
             ENDDO
             pool_start(1,1,iele) = pool_start(1,1,iele) + &
                  (bm_to_litter(ipts,ivm,ipar,iele) * veget_max(ipts,ivm))
          ENDDO

          ! The biomass harvest pool shouldn't be multiplied by 
          ! veget_max its units are already in gC pixel-1. Right 
          ! now the harvest pool is in gC.  We will divide by the
          ! area to make it per pixel. 
          pool_start(1,1,iele) = pool_start(1,1,iele) + &
               SUM(SUM(harvest_pool(ipts,ivm,:,iele,:),2)) / &
               (area(ipts) * contfrac(ipts))
          
       ENDDO 
       
    ENDIF ! err_act.GT.1


    IF( SUM(circ_class_kill(ipts,ivm,:,ifm,icut)) .GT. min_stomate ) THEN

       CALL harvest_aboveground_dead_roots(npts, ipts, ivm, &
            ifm, icut, bm_to_litter, woody_litter_by_cut, &
            circ_class_biomass, circ_class_kill, circ_class_n, harvest_pool, &
            harvest_type, harvest_cut, harvest_area, &
            harvest_pool_bound, veget_max, harvest_fraction)

       ! thinning trees generally avoids natural mortality and self-thinning, 
       ! so let's remove all thinned trees from the pools.  This is tricky 
       ! because mortality depends on a fraction of the total biomass, which 
       ! may be spread out over different circ classes than what were killed 
       ! by thinning.
       thinned_bm = zero
       dead_bm = zero
       DO icir=1,ncirc
          total_bm=SUM(circ_class_biomass(ipts,ivm,icir,:,icarbon))
          thinned_bm=thinned_bm+circ_class_kill(ipts,ivm,icir,ifm,icut)*&
               total_bm
          dead_bm=dead_bm+circ_class_kill(ipts,ivm,icir,ifm_none,icut_thin)*&
               total_bm
       ENDDO

       bm_difference = dead_bm - thinned_bm
       IF(bm_difference .LE. min_stomate)THEN

          ! if the thinned biomass is more than what would have died from 
          ! natural causes, we can just assume that thinning killed all 
          ! these trees. Zero out natural mortality.
          circ_class_kill(ipts,ivm,:,ifm_none,icut_thin)=zero

       ELSE

          ! if not, we have to reduce the amount of biomass that is 
          ! killed naturally by what was removed by thinning, 
          ! redistributing the remainder properly across all circ 
          ! classes as we did in stomate_mark_kill.
          circ_class_kill(ipts,ivm,:,ifm_none,icut_thin) = zero
          CALL distribute_mortality_biomass( bm_difference, &
               death_distribution_factor(ivm), &
               circ_class_n(ipts,ivm,:), &
               circ_class_biomass(ipts,ivm,:,:,icarbon), &
               circ_class_kill(ipts,ivm,:,ifm_none,icut_thin) )
       ENDIF ! bm_difference .LE. min_stomate

       ! make sure nothing has gone negative
       WHERE ( circ_class_kill(ipts,ivm,:,ifm_none,icut_thin) .LT. zero )
          circ_class_kill(ipts,ivm,:,ifm_none,icut_thin)=zero
       END WHERE

       ! now zero out this circ_class_kill, since we should never try
       ! to kill these same trees again.
       circ_class_kill(ipts,ivm,:,ifm,icut)=zero

       ! Reset the values for the time since the last operation
       last_cut(ipts,ivm)=0
       
    ENDIF ! checking to see if we have biomass to thin


 !! Check numerical consistency of this routine

    IF (err_act.GT.1) THEN

       ! 4.3 Mass balance closure 
       pool_end(:,:,:) = zero
       DO iele = 1,nelements
          DO ipar = 1,nparts
             DO icir = 1,ncirc
                pool_end(1,1,iele) = pool_end(1,1,iele) + &
                     (circ_class_biomass(ipts,ivm,icir,ipar,iele) * &
                     circ_class_n(ipts,ivm,icir) * veget_max(ipts,ivm))
             ENDDO
             pool_end(1,1,iele) = pool_end(1,1,iele) + &
                  (bm_to_litter(ipts,ivm,ipar,iele) * veget_max(ipts,ivm))
          ENDDO

          ! The biomass harvest pool shouldn't be multiplied by veget_max
          ! its units are already in gC pixel-1. Right now the harvest 
          ! pool is in gC.  We will divide by the area to make it 
          ! per pixel.
          pool_end(1,1,iele) = pool_end(1,1,iele) + &
               SUM(SUM(harvest_pool(ipts,ivm,:,iele,:),2)) / &
               (area(ipts) * contfrac(ipts))
       ENDDO ! nelements
   
       ! 4.3.2 Calculate mass balance
       ! Common processes
       DO iele=1,nelements
          check_intern(1,1,ipoolchange,iele) = -un * (pool_end(1,1,iele) - &
               pool_start(1,1,iele)) 
       ENDDO

       closure_intern(:,:,:) = zero
       DO imbc = 1,nmbcomp
          DO iele=1,nelements
             ! Debug
             IF (printlev_loc>=4) WRITE(numout,*) &
                  'check_intern, ivm, imbc, iele, ', imbc, &
                  iele, check_intern(1,1,imbc,iele)
             !-
             closure_intern(1,1,iele) = closure_intern(1,1,iele) + &
                  check_intern(1,1,imbc,iele)
          ENDDO
       ENDDO

       ! 4.3.3 Check mass balance closure
       CALL check_mass_balance("thinning_harvest_aboveground_dead_roots", &
            closure_intern, npts, pool_end, pool_start, veget_max, 'ipts')
      
    ENDIF ! err_act.GT.1

    IF (printlev.GE.4) WRITE(numout,*) &
         'Leaving thinning_harvest_aboveground_dead_roots'

  END SUBROUTINE thinning_harvest_aboveground_dead_roots

!! ================================================================================================================================
!! SUBROUTINE   : harvest_aboveground_dead_roots
!!
!>\BRIEF        Transfer dead biomass to litter and wood to harvest pool.
!!
!! DESCRIPTION  : Moves dead biomass to the litter pool and the harvest pool
!!                where appropriate.  
!!
!!                NOTE: If harvest_fraction is 1, all the aboveground biomass is
!!                      harvested.  If it's 1-branch_ratio, only stems are 
!!                      harvested.  Belowground wood, roots, fruits, and
!!                      leaves all become litter.  A portion of aboveground
!!                      wood is put in the harvest pool.  Reserves and labile
!!                      carbon are assumed to be distributed evenly across
!!                      all non-heartwood components.
!!
!! RECENT CHANGE(S): 
!!
!! MAIN OUTPUT VARIABLE(S): ::bm_to_litter, ::harvest_pool
!!
!! REFERENCE(S)   :
!!
!! FLOWCHART    : None
!!\n
!_ ================================================================================================================================
 
  SUBROUTINE harvest_aboveground_dead_roots (npts, ipts, ivm, &
       ifm, icut, bm_to_litter, woody_litter_by_cut, &
       circ_class_biomass, circ_class_kill, circ_class_n, harvest_pool, &
       harvest_type, harvest_cut, harvest_area, &
       harvest_pool_bound, veget_max, harvest_fraction)


    !! 0. Variable and parameter declaration

    !! 0.1 Input variables 
    ! NOTE: if harvest_pool_bound is not explicitely defined 
    ! the model uses 1:ndia_harvest+2 as its dimensions. If 
    ! you like to do so remember to change the counters in
    ! this routine
    INTEGER(i_std), INTENT(in)                                   :: npts               !! Domain size (-)
    INTEGER(i_std), INTENT(in)                                   :: ipts               !! current grid square
    INTEGER(i_std), INTENT(in)                                   :: ivm                !! current PFT
    INTEGER(i_std), INTENT(in)                                   :: ifm                !! current fm strategy
    INTEGER(i_std), INTENT(in)                                   :: icut               !! current cut time
    REAL(r_std), DIMENSION(0:), INTENT(in)                       :: harvest_pool_bound !! The boundaries of the diameter classes
                                                                                       !! in the wood harvest pools
                                                                                       !! @tex $(m)$ @endtex
    REAL(r_std),DIMENSION(:,:),INTENT(in)                        :: veget_max          !! Maximum fraction of vegetation type 
    REAL(r_std),INTENT(in)                                       :: harvest_fraction   !! The fraction of the aboveground
                                                                                       !! biomass that is harvested.  This is
                                                                                       !! generally 1-branch_ratio if you only
                                                                                       !! want to harvest the stems, and 1 if
                                                                                       !! you want to harvest everything.

    !! 0.2 Output variables


    !! 0.3 Modified variables   

    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)             :: circ_class_biomass !! Biomass of the components of the model  
                                                                                       !! tree within a circumference
                                                                                       !! class @tex $(gC ind^{-1})$ @endtex  
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)                 :: circ_class_n       !! Number of individuals in each circ class
                                                                                       !! @tex $(m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)             :: circ_class_kill    !! Number of trees within a circ that needs
                                                                                       !! to be killed @tex $(ind m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)               :: bm_to_litter       !! Biomass transfer to litter 
                                                                                       !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)             :: harvest_pool       !! The wood and biomass that have been
                                                                                       !! havested by humans @tex $(gC)$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)                   :: harvest_type       !! Type of management that resulted
                                                                                       !! in the harvest (unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)                   :: harvest_cut        !! Type of cutting that was used for the harvest
                                                                                       !! (unitless)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)                 :: harvest_area       !! Harvested area (m^{2})
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)                 :: woody_litter_by_cut!! Saved woody litter by icut not by iparts
                                                                                       !! @tex $(gC m^{-2})$ @endtex


    !! 0.4 Local variables
    INTEGER(i_std)                                               :: idia, icir, imbc   !! Indices
    INTEGER(i_std)                                               :: iele, ipar         !! Indices
    LOGICAL                                                      :: ph_area_scaled     !! PROGRESSIVE_HARVEST: credit only the fraction actually cut
    REAL(r_std), DIMENSION(nelements)                            :: total_nonheart_biomass_harvest
    REAL(r_std), DIMENSION(nelements)                            :: total_nonheart_biomass_kill
    REAL(r_std), DIMENSION(nelements)                            :: harvest_ratio
    REAL(r_std), DIMENSION(ncirc)                                :: diameters_temp
    REAL(r_std), DIMENSION(npts,nvm,nmbcomp,nelements)           :: check_intern       !! Contains the components of the internal
                                                                                       !! mass balance chech for this routine
                                                                                       !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)                   :: closure_intern     !! Check closure of internal mass balance
                                                                                       !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)                   :: pool_start         !! Start pool of this routine 
                                                                                       !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)                   :: pool_end           !! End pool of this routine 
                                                                                       !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex

    !_ ================================================================================================================================

    IF (printlev.GE.3) WRITE(numout,*) &
         'Entering harvest_aboveground_dead_roots sapiens_kill.f90'

    !! 1.1 Initialize local printlev if not already done
    IF (firstcall_sapiens_kill) THEN
       !! Initialize local printlev
      ! Guillaume M. -- get_printlev calls getin_p, a COLLECTIVE bcast. Called here,
      ! inside the compute loop under a local firstcall, it hangs the run as soon as the
      ! ranks stop reaching it together. Use printlev, read once at startup by every rank.
       printlev_loc=printlev
       firstcall_sapiens_kill=.FALSE.
    END IF

    ! 1.2 Initialize check for mass balance closure
    IF (err_act.GT.1) THEN
       
       check_intern(:,:,:,:) = zero
       pool_start = zero
       DO iele = 1,nelements
          DO ipar = 1,nparts
             ! Initial litter and biomass. This routine is called
             ! within an ipts-LOOP and an ivm-LOOP so we test one 
             ! pixel and one pft at a time. Store in the first position.
             DO icir = 1,ncirc
                pool_start(1,1,iele) = pool_start(1,1,iele) + &
                     (circ_class_biomass(ipts,ivm,icir,ipar,iele) * &
                     circ_class_n(ipts,ivm,icir) * veget_max(ipts,ivm))
             ENDDO
             pool_start(1,1,iele) = pool_start(1,1,iele) + &
                  (bm_to_litter(ipts,ivm,ipar,iele) * veget_max(ipts,ivm))
          ENDDO

          ! The biomass harvest pool shouldn't be multiplied by 
          ! veget_max its units are already in gC pixel-1. Right 
          ! now the harvest pool is in gC.  We will divide by the
          ! area to make it per pixel. 
          pool_start(1,1,iele) = pool_start(1,1,iele) + &
               SUM(SUM(harvest_pool(ipts,ivm,:,iele,:),2)) / &
               (area(ipts)*contfrac(ipts))
          
       ENDDO 
       
    ENDIF ! err_act.GT.1
 
    ! We will need the trunk diameters for all the circ classes 
    ! in order to figure out which harvest pool to put them in.
    diameters_temp(:)=wood_to_dia(circ_class_biomass(ipts,ivm,:,:,icarbon),&
         ivm,pipe_tune2(ipts,ivm)) 

    DO icir=1,ncirc

       ! It's possible that we only kill individuals in certain circ classes.
       ! In this case, we need to check for each circ class to make sure
       ! we don't divide by zero.
       IF(circ_class_kill(ipts,ivm,icir,ifm,icut) .LT. min_stomate) CYCLE
       ! By any chance, the total anthropogenic mortality can exceed
       ! circ_class_n, which causes negative circ_class_n.
       IF(circ_class_kill(ipts,ivm,icir,ifm,icut) .GT. &
                circ_class_n(ipts,ivm,icir)) THEN
          circ_class_kill(ipts,ivm,icir,ifm,icut) = circ_class_n(ipts,ivm,icir)
       ENDIF

       ! Debug
       IF(printlev_loc>=4 .AND. ipts == test_grid .AND. ivm == test_pft)THEN
          WRITE(numout,*) 'Killing trees in harvest_aboveground_dead_roots'
          WRITE(numout,*) 'ipts,ivm,icir,ifm,icut: ',ipts,ivm,icir,ifm,icut
          WRITE(numout,*) 'circ_class_kill,circ_class_n: ',&
               circ_class_kill(ipts,ivm,icir,ifm,icut),&
               circ_class_n(ipts,ivm,icir)
       ENDIF
       !-

       ! The decisions of what to do with the carboyhydrate reserves and the
       ! labile carbon pool are tricky. In turnover, they are not touched
       ! when leaves/roots/fruits die, which implies that they are stored
       ! primarily in the above and belowground sapwood (heartwood is dead). 
       ! However, here a tree does not know thinning is coming, so it cannot
       ! move this carbon around.  We will therefore determine where is goes
       ! based on the total amount of non-heartwood pools and where they go.
       total_nonheart_biomass_kill(:)=circ_class_kill(ipts,ivm,icir,ifm,icut)*&
            (circ_class_biomass(ipts,ivm,icir,isapbelow,:)+&
            circ_class_biomass(ipts,ivm,icir,iroot,:)+&
            circ_class_biomass(ipts,ivm,icir,ileaf,:)+&
            circ_class_biomass(ipts,ivm,icir,ifruit,:)+&
            circ_class_biomass(ipts,ivm,icir,isapabove,:))

       total_nonheart_biomass_harvest(:)=harvest_fraction*&
            circ_class_kill(ipts,ivm,icir,ifm,icut)*&
            circ_class_biomass(ipts,ivm,icir,isapabove,:)

       ! This is now the fraction of the nonheartwood biomass that will 
       ! go into the harvest pool.  We will use this to partition the 
       ! carbohydrate reserve and labile pools into litter and harvest.
       DO iele = 1,nelements
          IF (total_nonheart_biomass_kill(iele).LE.zero) THEN
             harvest_ratio(iele) = un
          ELSE
             harvest_ratio(iele)=total_nonheart_biomass_harvest(iele) / &
                  total_nonheart_biomass_kill(iele)
          END IF
       END DO

       ! Guillaume M. -- Probe: which FACTOR is non finite? The writes below are products
       ! biomass x kill with no division, so a NaN can only come from an already non finite
       ! factor, or from Inf x 0 -- the latter poisons 9 compartments x 2 elements at once.
       ! /!\ Test isnan AND abs > 1e30: a diluted 1e33 sentinel gives an Inf, not a NaN, and
       ! a probe asking isnan alone misses it.
       IF ( ANY(isnan(circ_class_biomass(ipts,ivm,icir,:,:))) .OR. &
            ANY(ABS(circ_class_biomass(ipts,ivm,icir,:,:)) > 1.e30) .OR. &
            isnan(circ_class_kill(ipts,ivm,icir,ifm,icut)) .OR. &
            ABS(circ_class_kill(ipts,ivm,icir,ifm,icut)) > 1.e30 ) THEN
          WRITE(numout,*) '[KILLNF] ipts',ipts,' ipft',ivm,' icir',icir, &
               ' ifm',ifm,' icut',icut, &
               ' | bm nan',ANY(isnan(circ_class_biomass(ipts,ivm,icir,:,:))), &
               ' bm inf',ANY(ABS(circ_class_biomass(ipts,ivm,icir,:,:)) > 1.e30), &
               ' | kill',circ_class_kill(ipts,ivm,icir,ifm,icut), &
               ' | bm max',MAXVAL(ABS(circ_class_biomass(ipts,ivm,icir,:,:)))
          CALL flush(numout)
       ENDIF

       ! all of the following biomass becomes litter
       bm_to_litter(ipts,ivm,isapbelow,:) = &
            bm_to_litter(ipts,ivm,isapbelow,:) + &
            circ_class_biomass(ipts,ivm,icir,isapbelow,:)*&
            circ_class_kill(ipts,ivm,icir,ifm,icut)
       bm_to_litter(ipts,ivm,iheartbelow,:) = &
            bm_to_litter(ipts,ivm,iheartbelow,:) + &
            circ_class_biomass(ipts,ivm,icir,iheartbelow,:)*&
            circ_class_kill(ipts,ivm,icir,ifm,icut)
       bm_to_litter(ipts,ivm,iroot,:) = bm_to_litter(ipts,ivm,iroot,:) + &
            circ_class_biomass(ipts,ivm,icir,iroot,:)*&
            circ_class_kill(ipts,ivm,icir,ifm,icut)
       bm_to_litter(ipts,ivm,ileaf,:) = bm_to_litter(ipts,ivm,ileaf,:) + &
            circ_class_biomass(ipts,ivm,icir,ileaf,:)*&
            circ_class_kill(ipts,ivm,icir,ifm,icut)
       bm_to_litter(ipts,ivm,ifruit,:) = &
            bm_to_litter(ipts,ivm,ifruit,:) + &
            circ_class_biomass(ipts,ivm,icir,ifruit,:)*&
            circ_class_kill(ipts,ivm,icir,ifm,icut)

       ! Include the branches as litter, if necessary.
       bm_to_litter(ipts,ivm,isapabove,:) = &
            bm_to_litter(ipts,ivm,isapabove,:) + &
            (un - harvest_fraction)*&
            circ_class_biomass(ipts,ivm,icir,isapabove,:)*&
            circ_class_kill(ipts,ivm,icir,ifm,icut)
       bm_to_litter(ipts,ivm,iheartabove,:) = &
            bm_to_litter(ipts,ivm,iheartabove,:) + &
            (un - harvest_fraction)*&
            circ_class_biomass(ipts,ivm,icir,iheartabove,:)*&
            circ_class_kill(ipts,ivm,icir,ifm,icut)

       ! Save woody litter by icut. This will be used in pest module
       ! (stomate_pest) for one of the susceptibility indice.
       woody_litter_by_cut(ipts,ivm,icut) = woody_litter_by_cut(ipts,ivm,icut) + &
             (un - harvest_fraction) * &
             ( circ_class_biomass(ipts,ivm,icir,isapabove,icarbon) + &
                circ_class_biomass(ipts,ivm,icir,iheartabove,icarbon) )* &
             circ_class_kill(ipts,ivm,icir,ifm,icut)

       ! Don't forget to include the proper proportion of the labile and
       ! carbres pools.
       bm_to_litter(ipts,ivm,icarbres,:) = &
            bm_to_litter(ipts,ivm,icarbres,:) + &
            (un - harvest_ratio(:))*&
            circ_class_biomass(ipts,ivm,icir,icarbres,:)*&
            circ_class_kill(ipts,ivm,icir,ifm,icut)
       bm_to_litter(ipts,ivm,ilabile,:) = &
            bm_to_litter(ipts,ivm,ilabile,:) + &
            (un - harvest_ratio(:))*&
            circ_class_biomass(ipts,ivm,icir,ilabile,:)*&
            circ_class_kill(ipts,ivm,icir,ifm,icut)

       ! now we harvest the stems.  These are what we want to save, and we are
       ! going to classify them by PFT type, location, and diameter so that
       ! we can compare to FAO statistics and do some more interesting
       ! post-treatment.

       ! UNITS: The units of the biomass pools are per meter squared.
       ! It seems to make more sense to store the harvest pool in absolute
       ! units, which means we have to multiple the result by the fraction of
       ! this grid cell covered by this PFT (veget_max) and the size of the 
       ! grid cell.  It seems that the fraction of nobio land does not need 
       ! to be used here, since frac_nobio and veget_max are summed together 
       ! in slowproc and they must equal one.  This gives
       ! us the grams of carbon harvest, which is a big number but that's
       ! why people invented scientific notation.
       dia_loop:DO idia=1,ndia_harvest

          ! Move the harvest to the correct diameter class and
          ! account for the largest diameter class because
          ! harvest_pool_bound has only ndia_harvest+1 dimensions
          IF ( (diameters_temp(icir).LE.harvest_pool_bound(idia)).OR.&
               ((diameters_temp(icir).GT.harvest_pool_bound(ndia_harvest)).AND.&
               (diameters_temp(icir).LT.harvest_pool_bound(ndia_harvest+1))))& 
               THEN

             ! this is the correct diameter class 
             harvest_pool(ipts,ivm,idia,:,iharvest) = &
                  harvest_pool(ipts,ivm,idia,:,iharvest) + &
                  harvest_fraction*circ_class_kill(ipts,ivm,icir,ifm,icut)*&
                  (circ_class_biomass(ipts,ivm,icir,isapabove,:)+&
                  circ_class_biomass(ipts,ivm,icir,iheartabove,:)) * &
                  veget_max(ipts,ivm) * &
                  area(ipts) * contfrac(ipts)

             ! also add the proper fractions of the labile and carbres pools
             harvest_pool(ipts,ivm,idia,:,iharvest) = &
                  harvest_pool(ipts,ivm,idia,:,iharvest)+ &
                  harvest_ratio(:) * circ_class_kill(ipts,ivm,icir,ifm,icut)*&
                  ( circ_class_biomass(ipts,ivm,icir,ilabile,:) + &
                  circ_class_biomass(ipts,ivm,icir,icarbres,:) ) * &
                  veget_max(ipts,ivm) * &
                  area(ipts) * contfrac(ipts)

             ! The harvest of tree diameter icir is send to a single idia
             ! once this idia has been found, there is no need to loop
             ! over the remaining dia classes.  
             EXIT dia_loop

          ENDIF ! diameters_temp(icir) .LE. harvest_pool_bound(idia)

          ! If we are here, it means the diameter of this trunk is larger
          ! than the upper boundary of our harvest bins.  Considering the
          ! upper boundary is set to be 99999 meters in constants.f90,
          ! something has gone wrong!
          IF (diameters_temp(icir).GE.harvest_pool_bound(ndia_harvest+1)) THEN

             WRITE(numout,*) 'ERROR: Our tree is too big to fit '&
                  // 'into the harvest pools.'
             WRITE(numout,*) 'ipts, ivm, icir: ',ipts, ivm, icir
             WRITE(numout,*) 'diameters_temp(icir), harvest_pool_bound(idia)',&
                  diameters_temp(icir),harvest_pool_bound(ndia_harvest+1)
             CALL ipslerr_p (3,'harvest_aboveground_dead_roots', &
                  'Tree too big to fit into the harvest pool.',&
                  'Look into the output file for details.','')
          ENDIF

       ENDDO dia_loop

       ! remove the number of individuals that were harvested
       circ_class_n(ipts,ivm,icir) = circ_class_n(ipts,ivm,icir) - &
            circ_class_kill(ipts,ivm,icir,ifm,icut)

    ENDDO ! loop over circ classes

    ! Associated harvest variables
    ! If several cuts happen in the same year, these variables can no 
    ! longer be correctly interpreted.
    harvest_type(ipts,ivm) = harvest_type(ipts,ivm) + ifm
    harvest_cut(ipts,ivm) = harvest_cut(ipts,ivm) + icut
    ! Guillaume M. -- Crediting the WHOLE slot area whatever the number of trees killed
    ! overstates HARVEST_AREA well above what the rotation allows. Under
    ! ok_progressive_harvest a clearcut takes only ph_cut_frac of the slot, so only that
    ! share counts; thinnings and coppice cuts keep whole-slot semantics. /!\ The
    ! ph_cut_frac > 0 test is required: other clearcut routes carry 0 and would credit ZERO.
    ph_area_scaled = .FALSE.
    IF (ok_progressive_harvest .AND. icut .EQ. icut_clear .AND. ALLOCATED(ph_cut_frac)) THEN
       IF (ph_cut_frac(ipts,ivm) .GT. zero) ph_area_scaled = .TRUE.
    ENDIF
    IF (ph_area_scaled) THEN
       harvest_area(ipts,ivm,iharvest) = harvest_area(ipts,ivm,iharvest) + &
            ph_cut_frac(ipts,ivm) * veget_max(ipts,ivm) * area(ipts) * contfrac(ipts)
    ELSE
       harvest_area(ipts,ivm,iharvest) = harvest_area(ipts,ivm,iharvest) + &
            veget_max(ipts,ivm) * area(ipts) * contfrac(ipts)
    ENDIF
    
    !! Check numerical consistency of this routine

    IF (err_act.GT.1) THEN

       ! 4.3 Mass balance closure 
       pool_end(:,:,:) = zero
       DO iele = 1,nelements
          DO ipar = 1,nparts
             DO icir = 1,ncirc
                pool_end(1,1,iele) = pool_end(1,1,iele) + &
                     (circ_class_biomass(ipts,ivm,icir,ipar,iele) * &
                     circ_class_n(ipts,ivm,icir) * veget_max(ipts,ivm))
             ENDDO
             pool_end(1,1,iele) = pool_end(1,1,iele) + &
                  (bm_to_litter(ipts,ivm,ipar,iele) * veget_max(ipts,ivm))
          ENDDO

          ! The biomass harvest pool shouldn't be multiplied by veget_max
          ! its units are already in gC pixel-1. Right now the harvest 
          ! pool is in gC.  We will divide by the area to make it 
          ! per pixel.
          pool_end(1,1,iele) = pool_end(1,1,iele) + &
               SUM(SUM(harvest_pool(ipts,ivm,:,iele,:),2)) / &
               (area(ipts) * contfrac(ipts))
       ENDDO ! nelements
   
       ! 4.3.2 Calculate mass balance
       ! Common processes
       DO iele=1,nelements
          check_intern(1,1,ipoolchange,iele) = -un * (pool_end(1,1,iele) - &
               pool_start(1,1,iele)) 
       ENDDO

       closure_intern(:,:,:) = zero
       DO imbc = 1,nmbcomp
          DO iele=1,nelements
             ! Debug
             IF (printlev_loc>=4) WRITE(numout,*) &
                  'check_intern, ivm, imbc, iele, ', imbc, &
                  iele, check_intern(1,1,imbc,iele)
             !-
             closure_intern(1,1,iele) = closure_intern(1,1,iele) + &
                  check_intern(1,1,imbc,iele)
          ENDDO
       ENDDO

       ! 4.3.3 Check mass balance closure
       CALL check_mass_balance("harvest_aboveground_dead_roots", closure_intern, &
            npts, pool_end, pool_start, veget_max, 'ipts')
      
    ENDIF ! err_act.GT.1

    IF (printlev.GE.4) WRITE(numout,*) &
         'Leaving harvest_aboveground_dead_roots'

  END SUBROUTINE harvest_aboveground_dead_roots

!! ================================================================================================================================
!! SUBROUTINE   : harvest_aboveground_wood_living_roots
!!
!>\BRIEF        Transfer dead biomass to litter and wood to harvest pool.
!!
!! DESCRIPTION  : Moves dead biomass to the litter pool and the harvest pool
!!                where appropriate.  
!!
!!                NOTE: All branches and stems are harvested.  Basically, this means
!!                      all aboveground woody biomass.  The fine roots
!!                      are killed and left to litter.  The sapwood and heartwood
!!                      belowground is left alive.  The leaves are killed, and the
!!                      labile and carbres pools are left untouched in the living
!!                      biomass.
!!
!! RECENT CHANGE(S): 
!!
!! MAIN OUTPUT VARIABLE(S): ::bm_to_litter, ::harvest_pool
!!
!! REFERENCE(S)   :
!! - Bellassen, V., Le Maire, G., Dhote, J.F., Viovy, N., Ciais, P., 2010. 
!! Modeling forest management within a global vegetation model – Part 1: 
!! model structure and general behaviour. Ecological Modelling 221, 2458–2474.
!!
!! FLOWCHART    : None
!!\n
!_ ================================================================================================================================
 
  SUBROUTINE harvest_aboveground_wood_living_roots (npts, ipts, ivm, &
       ifm, icut, bm_to_litter, circ_class_n, &
       circ_class_biomass, circ_class_kill, harvest_pool, &
       harvest_type, harvest_cut, harvest_area, &
       harvest_pool_bound, veget_max)


    !! 0. Variable and parameter declaration

    !! 0.1 Input variables 
    ! NOTE: if harvest_pool_bound is not explicitely defined 
    ! the model uses 1:ndia_harvest+2 as its dimensions. If 
    ! you like to do so remember to change the counters in
    ! this routine
    INTEGER(i_std), INTENT(in)                                   :: npts               !! Domain size (-)
    INTEGER(i_std), INTENT(in)                                   :: ipts               !! current grid square
    INTEGER(i_std), INTENT(in)                                   :: ivm                !! current PFT
    INTEGER(i_std), INTENT(in)                                   :: ifm                !! current fm strategy
    INTEGER(i_std), INTENT(in)                                   :: icut               !! current cut time
    REAL(r_std), DIMENSION(0:), INTENT(in)                       :: harvest_pool_bound !! The boundaries of the diameter classes
                                                                                       !! in the wood harvest pools
                                                                                       !! @tex $(m)$ @endtex
    REAL(r_std),DIMENSION(:,:),INTENT(in)                        :: veget_max          !! Maximum fraction of vegetation type 


    !! 0.2 Output variables


    !! 0.3 Modified variables
   
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)                 :: circ_class_n       !! Number of individuals in each circ class
                                                                                       !! @tex $(m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)             :: circ_class_biomass !! Biomass of the components of the model  
                                                                                       !! tree within a circumference
                                                                                       !! class @tex $(gC ind^{-1})$ @endtex  
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)             :: circ_class_kill    !! Number of trees within a circ that needs
                                                                                       !! to be killed @tex $(ind m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)               :: bm_to_litter       !! Biomass transfer to litter 
                                                                                       !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)             :: harvest_pool       !! The wood and biomass that have been
                                                                                       !! havested by humans @tex $(gC)$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)                   :: harvest_type       !! Type of management that resulted
                                                                                       !! in the harvest (unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)                   :: harvest_cut        !! Type of cutting that was used for the harvest
                                                                                       !! (unitless)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)                 :: harvest_area       !! Harvested area (m^{2})

    !! 0.4 Local variables

    INTEGER(i_std)                                               :: idia, icir,imbc    !! Indices
    INTEGER(i_std)                                               :: iele, ipar         !! Indices
    LOGICAL                                                      :: ph_area_scaled     !! PROGRESSIVE_HARVEST: credit only the fraction actually cut
    REAL(r_std), DIMENSION(nelements)                            :: total_nonheart_biomass_harvest
    REAL(r_std), DIMENSION(nelements)                            :: total_nonheart_biomass_living
    REAL(r_std), DIMENSION(nelements)                            :: total_nonheart_biomass_kill
    REAL(r_std), DIMENSION(nelements)                            :: harvest_ratio
    REAL(r_std), DIMENSION(nelements)                            :: living_ratio
    REAL(r_std), DIMENSION(nelements)                            :: litter_ratio
    REAL(r_std), DIMENSION(ncirc)                                :: diameters_temp
    REAL(r_std), DIMENSION(nelements)                            :: dead_heartwood_be
    REAL(r_std), DIMENSION(nelements)                            :: dead_sapwood_be
    REAL(r_std)                                                  :: cc_taken           !! Trees actually taken from this
                                                                                       !! class, bounded by the population
                                                                                       !! @tex $(ind m^{-2})$ @endtex
    REAL(r_std)                                                  :: cop_frac           !! Share of the class taken (-)
    REAL(r_std), DIMENSION(npts,nvm,nmbcomp,nelements)           :: check_intern       !! Contains the components of the internal
                                                                                       !! mass balance chech for this routine
                                                                                       !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)                   :: closure_intern     !! Check closure of internal mass balance
                                                                                       !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)                   :: pool_start         !! Start pool of this routine 
                                                                                       !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)                   :: pool_end           !! End pool of this routine 
                                                                                       !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
 
    !_ ================================================================================================================================

    IF (printlev.GE.3) WRITE(numout,*) &
         'Entering harvest_aboveground_wood_living_roots sapiens_kill.f90'

        ! 1.2 Initialize check for mass balance closure
    IF (err_act.GT.1) THEN
       
       check_intern(:,:,:,:) = zero
       pool_start = zero
       DO iele = 1,nelements
          DO ipar = 1,nparts
             ! Initial litter and biomass. This routine is called
             ! within an ipts-LOOP and an ivm-LOOP so we test one 
             ! pixel and one pft at a time. Store in the first position.
             DO icir = 1,ncirc
                pool_start(1,1,iele) = pool_start(1,1,iele) + &
                     (circ_class_biomass(ipts,ivm,icir,ipar,iele) * &
                     circ_class_n(ipts,ivm,icir) * veget_max(ipts,ivm))
             ENDDO
             pool_start(1,1,iele) = pool_start(1,1,iele) + &
                  (bm_to_litter(ipts,ivm,ipar,iele) * veget_max(ipts,ivm))
          ENDDO

          ! The biomass harvest pool shouldn't be multiplied by 
          ! veget_max its units are already in gC pixel-1. Right 
          ! now the harvest pool is in gC.  We will divide by the
          ! area to make it per pixel. 
          pool_start(1,1,iele) = pool_start(1,1,iele) + &
               SUM(SUM(harvest_pool(ipts,ivm,:,iele,:),2)) / &
               (area(ipts) * contfrac(ipts))
          
       ENDDO 
       
    ENDIF ! err_act.GT.1

    ! we will need the trunk diameters for all the circ classes in 
    ! order to figure out which harvest pool to put them in.
    diameters_temp(:)=wood_to_dia(circ_class_biomass(ipts,ivm,:,:,icarbon),&
         ivm,pipe_tune2(ipts,ivm)) 

    DO icir=1,ncirc

       ! It's possible that we only kill individuals in certain circ classes.
       ! In this case, we need to check for each circ class to make sure
       ! we don't divide by zero.
       IF(circ_class_kill(ipts,ivm,icir,ifm,icut) .LT. min_stomate) CYCLE

       ! Guillaume M. -- circ_class_kill is set equal to circ_class_n when the cut is decided,
       ! but earlier cuts in this same dispatcher shrink circ_class_n before we get here, so
       ! kill can exceed the surviving population. Book what is really there and remove exactly
       ! that share: the routine used to blank the whole class while crediting kill only, which
       ! created mass. Nominal case kill == n is unchanged (cop_frac = 1).
       cc_taken = MIN(circ_class_kill(ipts,ivm,icir,ifm,icut), circ_class_n(ipts,ivm,icir))
       IF (circ_class_n(ipts,ivm,icir) .GT. min_stomate) THEN
          cop_frac = cc_taken / circ_class_n(ipts,ivm,icir)
       ELSE
          cop_frac = zero
       ENDIF
       IF (cc_taken .LT. min_stomate) CYCLE

       ! Debug
       IF(printlev_loc>=4 .AND. ipts == test_grid .AND. ivm == test_pft)THEN
          WRITE(numout,*) 'Killing trees in '&
               // 'harvest_aboveground_wood_living_roots'
          WRITE(numout,*) 'ipts,ivm,icir,ifm,icut: ',ipts,ivm,icir,ifm,icut
          WRITE(numout,*) 'circ_class_kill: ',&
               circ_class_kill(ipts,ivm,icir,ifm,icut)
       ENDIF
       !-

       !+++CHECK+++
       ! Ratios are calculated but never used. The ratios are typically used
       ! to decide how to distribute the labile and reserve pools. is there a 
       ! good reason why it is not done here?

       ! The decision of what to do with the carboyhydrate reserves and the
       ! labile carbon pool are tricky.  In turnover, they are not touched
       ! when leaves/roots/fruits die, which implies that they are stored
       ! primarily in the above and belowground sapwood (heartwood is dead). 
       ! However, here a tree does not know thinning is coming, so it cannot
       ! move this carbon around.  In theory, we should remove part of the
       ! labile and reserve pools to the different product pools, but
       ! these pools are very sensitive for the allocation.  Therefore
       ! we don't touch them, which is equivalent to saying the reserve
       ! and labile carbon are all found in the roots.
       total_nonheart_biomass_kill(:)=cc_taken*&
            (circ_class_biomass(ipts,ivm,icir,isapbelow,:)+&
            circ_class_biomass(ipts,ivm,icir,iroot,:)+&
            circ_class_biomass(ipts,ivm,icir,ileaf,:)+&
            circ_class_biomass(ipts,ivm,icir,ifruit,:)+&
            circ_class_biomass(ipts,ivm,icir,isapabove,:))

       total_nonheart_biomass_harvest(:)=&
            cc_taken*&
            (circ_class_biomass(ipts,ivm,icir,isapabove,:))

       total_nonheart_biomass_living(:)=&
            cc_taken*&
            (circ_class_biomass(ipts,ivm,icir,isapbelow,:))

       ! This is now the fraction of the nonheartwood biomass that will 
       ! go into the harvest pool.  We will use this to partition the 
       ! carbohydrate reserve and labile pools into litter and harvest.
       harvest_ratio(:)=total_nonheart_biomass_harvest(:)/&
            total_nonheart_biomass_kill(:)
       living_ratio(:)=total_nonheart_biomass_living(:)/&
            total_nonheart_biomass_kill(:)
       litter_ratio(:)=un-harvest_ratio(:)-living_ratio(:)
       IF(litter_ratio(icarbon) .LT. -min_stomate .OR. &
            litter_ratio(initrogen) .LT. -min_stomate)THEN

          WRITE(numout,*) 'We have a problem with coppicing!'
          WRITE(numout,*) 'ipts,ivm,icir,ifm,icut ',ipts,ivm,icir,ifm,icut
          WRITE(numout,*) 'harvest_ratio, living_ratio, litter_ratio ',&
               harvest_ratio(:), living_ratio(:), litter_ratio(:)
          CALL ipslerr_p (3,'harvest_aboveground_wood_living_roots', &
               'litter_ratio is negative.',&
               'Look into the output file for details.','') 
       ENDIF
       !+++++++++++

       ! all of the following biomass becomes litter
       bm_to_litter(ipts,ivm,iroot,:) = bm_to_litter(ipts,ivm,iroot,:) + &
            circ_class_biomass(ipts,ivm,icir,iroot,:)*&
            cc_taken
       bm_to_litter(ipts,ivm,ileaf,:) = bm_to_litter(ipts,ivm,ileaf,:) + &
            circ_class_biomass(ipts,ivm,icir,ileaf,:)*&
            cc_taken
       bm_to_litter(ipts,ivm,ifruit,:) = bm_to_litter(ipts,ivm,ifruit,:) + &
            circ_class_biomass(ipts,ivm,icir,ifruit,:)*&
            cc_taken

       ! We also turn part of the belowground wood into litter. The reason
       ! for this is that some coppiced forests were not being able to 
       ! regenerate due to the amount of belowground wood, which kept the 
       ! trees out of balance for the whole next year.  The leaves grew, 
       ! but never enough to put wood into the stems, and the leaves died 
       ! at the end of the year. Then everything started over exactly the 
       ! same the following years. The biomass stayed constant (the 
       ! belowground mass was two orders of magnitude larger than any
       ! leaf mass which grew) and it never died.  Killing part of the 
       ! belowground mass should allow the tree to allometrically allocate 
       ! faster. At the same time, it's not entirely unrealistic.  Parts of 
       ! the roots stay alive during coppicing, but do they all?  Unknown.
       ! This is probably only an issue for full coppicing, not for short 
       ! rotation coppices, since those stands are much younger.
       dead_sapwood_be(:)=circ_class_biomass(ipts,ivm,icir,isapbelow,:)*&
            coppice_kill_be_wood(ivm)
       bm_to_litter(ipts,ivm,isapbelow,:) = &
            bm_to_litter(ipts,ivm,isapbelow,:) + dead_sapwood_be*&
            cc_taken
       dead_heartwood_be(:)=circ_class_biomass(ipts,ivm,icir,iheartbelow,:)*&
            coppice_kill_be_wood(ivm)
       bm_to_litter(ipts,ivm,iheartbelow,:) = &
            bm_to_litter(ipts,ivm,iheartbelow,:) + dead_heartwood_be*&
            cc_taken

       ! now we have the stems.  These are what we want to save, and we are
       ! going to classify them by PFT type, location, and diameter so that
       ! we can compare to FAO statistics and do some more interesting
       ! post-treatment.

       ! UNITS: The units of the biomass pools are per meter squared.
       ! It seems to make more sense to store the harvest pool in absolute
       ! units, which means we have to multiple the result by the fraction of
       ! this grid cell covered by this PFT (veget_max) and the size of the 
       ! grid cell.  It seems that the fraction of nobio land does not need 
       ! to be used here, since frac_nobio and veget_max are summed together 
       ! in slowproc and they must equal one.  This gives
       ! us the grams of carbon harvest, which is a big number but that's
       ! why people invented scientific notation.

       dia_loop:DO idia=1,ndia_harvest

          ! Move the harvest to the correct diameter class and
          ! account for the largest diameter class because
          ! harvest_pool_bound has only ndia_harvest+1 dimensions
          IF ( (diameters_temp(icir).LE.harvest_pool_bound(idia)) .OR. &
               ((diameters_temp(icir).GT.harvest_pool_bound(ndia_harvest)).AND.&
               (diameters_temp(icir).LT.harvest_pool_bound(ndia_harvest+1))) )&
               THEN

             ! this is the correct diameter class 
             harvest_pool(ipts,ivm,idia,:,iharvest)=&
                  harvest_pool(ipts,ivm,idia,:,iharvest)+ &
                  cc_taken*&
                  (circ_class_biomass(ipts,ivm,icir,isapabove,:)+&
                  circ_class_biomass(ipts,ivm,icir,iheartabove,:)) * &
                  veget_max(ipts,ivm) * &
                  area(ipts) * contfrac(ipts)

             ! Debug
             IF(printlev_loc>=4 .AND. ipts == test_grid .AND. &
                  ivm == test_pft)THEN
                WRITE(numout,*) 'harvest_aboveground_living_roots 1 ',&
                     ipts,ivm,icir
                WRITE(numout,*) 'harvest_area,harvest_type,harvest_cut: ',&
                     harvest_area(ipts,ivm,iharvest),harvest_type(ipts,ivm),&
                     harvest_cut(ipts,ivm)
                WRITE(numout,*) 'Woody biomass - C: ',&
                     circ_class_kill(ipts,ivm,icir,ifm,icut)*&
                     (circ_class_biomass(ipts,ivm,icir,isapabove,icarbon)+&
                     circ_class_biomass(ipts,ivm,icir,iheartabove,icarbon))
                WRITE(numout,*) 'Woody biomass - N: ',&
                     circ_class_kill(ipts,ivm,icir,ifm,icut)*&
                     (circ_class_biomass(ipts,ivm,icir,isapabove,initrogen)+&
                     circ_class_biomass(ipts,ivm,icir,iheartabove,initrogen))
             ENDIF
             !-

             EXIT dia_loop

          ENDIF ! diameters_temp(icir) .LE. harvest_pool_bound(idia)

          ! If we are here, it means the diameter of this trunk is larger
          ! than the upper boundary of our harvest bins.  Considering the
          ! upper boundary is set to be 99999 meters in constants.f90,
          ! something has gone wrong!
          IF (diameters_temp(icir).GE.harvest_pool_bound(ndia_harvest+1)) THEN

             WRITE(numout,*) 'Our tree is too big to fit into the harvest pools.'
             WRITE(numout,*) 'ipts, ivm, icir: ',ipts, ivm, icir
             WRITE(numout,*) 'diameters_temp(icir), harvest_pool_bound(idia)',&
                  diameters_temp(icir),harvest_pool_bound(ndia_harvest+1)
             CALL ipslerr_p (3,'harvest_aboveground_wood_living_roots', &
                  'Tree is too big for the harvest pool.',&
                  'Look into the output file for details.','')

          ENDIF

       ENDDO dia_loop

       ! No individuals die here since the roots are still alive
       ! Remove the dead biomass.  This includes everything except
       ! the belowground biomass.  
       circ_class_biomass(ipts,ivm,icir,ileaf,:)=&
            circ_class_biomass(ipts,ivm,icir,ileaf,:)*(un-cop_frac)
       circ_class_biomass(ipts,ivm,icir,isapabove,:)=&
            circ_class_biomass(ipts,ivm,icir,isapabove,:)*(un-cop_frac)
       circ_class_biomass(ipts,ivm,icir,iheartabove,:)=&
            circ_class_biomass(ipts,ivm,icir,iheartabove,:)*(un-cop_frac)
       circ_class_biomass(ipts,ivm,icir,iroot,:)=&
            circ_class_biomass(ipts,ivm,icir,iroot,:)*(un-cop_frac)
       circ_class_biomass(ipts,ivm,icir,ifruit,:)=&
            circ_class_biomass(ipts,ivm,icir,ifruit,:)*(un-cop_frac)

       ! Reduce the root biomass.
       circ_class_biomass(ipts,ivm,icir,isapbelow,:)=&
            circ_class_biomass(ipts,ivm,icir,isapbelow,:)*&
            (un-coppice_kill_be_wood(ivm)*cop_frac)
       circ_class_biomass(ipts,ivm,icir,iheartbelow,:)=&
            circ_class_biomass(ipts,ivm,icir,iheartbelow,:)*&
            (un-coppice_kill_be_wood(ivm)*cop_frac)

    ENDDO ! loop over circ classes

    ! Associated harvest variables
    ! If several cuts happen in the same year, these variables can no 
    ! longer be correctly interpreted.
    harvest_type(ipts,ivm) = harvest_type(ipts,ivm) + ifm
    harvest_cut(ipts,ivm) = harvest_cut(ipts,ivm) + icut
    ! Guillaume M. -- Crediting the WHOLE slot area whatever the number of trees killed
    ! overstates HARVEST_AREA well above what the rotation allows. Under
    ! ok_progressive_harvest a clearcut takes only ph_cut_frac of the slot, so only that
    ! share counts; thinnings and coppice cuts keep whole-slot semantics. /!\ The
    ! ph_cut_frac > 0 test is required: other clearcut routes carry 0 and would credit ZERO.
    ph_area_scaled = .FALSE.
    IF (ok_progressive_harvest .AND. icut .EQ. icut_clear .AND. ALLOCATED(ph_cut_frac)) THEN
       IF (ph_cut_frac(ipts,ivm) .GT. zero) ph_area_scaled = .TRUE.
    ENDIF
    IF (ph_area_scaled) THEN
       harvest_area(ipts,ivm,iharvest) = harvest_area(ipts,ivm,iharvest) + &
            ph_cut_frac(ipts,ivm) * veget_max(ipts,ivm) * area(ipts) * contfrac(ipts)
    ELSE
       harvest_area(ipts,ivm,iharvest) = harvest_area(ipts,ivm,iharvest) + &
            veget_max(ipts,ivm) * area(ipts) * contfrac(ipts)
    ENDIF
       
    !! Check numerical consistency of this routine

    IF (err_act.GT.1) THEN

       ! 4.3 Mass balance closure 
       pool_end(:,:,:) = zero
       DO iele = 1,nelements
          DO ipar = 1,nparts
             DO icir = 1,ncirc
                pool_end(1,1,iele) = pool_end(1,1,iele) + &
                     (circ_class_biomass(ipts,ivm,icir,ipar,iele) * &
                     circ_class_n(ipts,ivm,icir) * veget_max(ipts,ivm))
             ENDDO
             pool_end(1,1,iele) = pool_end(1,1,iele) + &
                  (bm_to_litter(ipts,ivm,ipar,iele) * veget_max(ipts,ivm))
          ENDDO

          ! The biomass harvest pool shouldn't be multiplied by veget_max
          ! its units are already in gC pixel-1. Right now the harvest 
          ! pool is in gC.  We will divide by the area to make it 
          ! per pixel.
          pool_end(1,1,iele) = pool_end(1,1,iele) + &
               SUM(SUM(harvest_pool(ipts,ivm,:,iele,:),2)) / &
               (area(ipts) * contfrac(ipts))
       ENDDO ! nelements
   
       ! 4.3.2 Calculate mass balance
       ! Common processes
       DO iele=1,nelements
          check_intern(1,1,ipoolchange,iele) = -un * (pool_end(1,1,iele) - &
               pool_start(1,1,iele)) 
       ENDDO

       closure_intern(:,:,:) = zero
       DO imbc = 1,nmbcomp
          DO iele=1,nelements
             ! Debug
             IF (printlev_loc>=4) WRITE(numout,*) &
                  'check_intern, ivm, imbc, iele, ', imbc, &
                  iele, check_intern(1,1,imbc,iele)
             !-
             closure_intern(1,1,iele) = closure_intern(1,1,iele) + &
                  check_intern(1,1,imbc,iele)
          ENDDO
       ENDDO

       ! 4.3.3 Check mass balance closure
       CALL check_mass_balance("harvest_aboveground_wood_living_roots", &
            closure_intern, npts, pool_end, pool_start, veget_max, 'ipts')
      
    ENDIF ! err_act.GT.1

    IF (printlev.GE.4) WRITE(numout,*) &
         'Leaving harvest_aboveground_wood_living_roots'

  END SUBROUTINE harvest_aboveground_wood_living_roots

!! ================================================================================================================================
!! SUBROUTINE   : second_coppice_cut
!!
!>\BRIEF        Simulates cutting a coppice between the first cut and the
!!              final cut.
!!
!! DESCRIPTION  : Moves dead biomass to the litter pool and the harvest pool
!!                where appropriate.  A few counters are reset.  The leftover
!!                biomass (the woody root system) is patitioned among enough
!!                individuals to give the same density as is found in coppice_dens,
!!                which is the density of stumps left after the first cut.
!!
!!                NOTE: All branches and stems are harvested.  Basically, this means
!!                      all aboveground woody biomass.  The fine roots
!!                      are killed and left to litter.  The sapwood and heartwood
!!                      belowground is left alive.  The leaves are killed, and the
!!                      labile and carbres pools are left untouched in the living
!!                      biomass.
!!
!! RECENT CHANGE(S): 
!!
!! MAIN OUTPUT VARIABLE(S): ::bm_to_litter, ::harvest_pool
!!
!! REFERENCE(S)   :
!! - Bellassen, V., Le Maire, G., Dhote, J.F., Viovy, N., Ciais, P., 2010. 
!! Modeling forest management within a global vegetation model – Part 1: 
!! model structure and general behaviour. Ecological Modelling 221, 2458–2474.
!!
!! FLOWCHART    : None
!!\n
!_ ================================================================================================================================
 
  SUBROUTINE second_coppice_cut (npts, ipts, ivm, ifm, &
       icut, harvest_pool_bound, veget_max, &
       circ_class_biomass, circ_class_n, circ_class_kill, &
       bm_to_litter, harvest_pool, harvest_type, harvest_cut, &
       harvest_area, last_cut, coppice_dens)


    !! 0. Variable and parameter declaration

    !! 0.1 Input variables 
    ! NOTE: if harvest_pool_bound is not explicitely defined 
    ! the model uses 1:ndia_harvest+2 as its dimensions. If 
    ! you like to do so remember to change the counters in
    ! this routine
    INTEGER(i_std), INTENT(in)                                   :: npts               !! Domain size (-)
    INTEGER(i_std), INTENT(in)                                   :: ipts               !! current grid square
    INTEGER(i_std), INTENT(in)                                   :: ivm                !! current PFT
    INTEGER(i_std), INTENT(in)                                   :: ifm                !! current fm strategy
    INTEGER(i_std), INTENT(in)                                   :: icut               !! current cut time
    REAL(r_std), DIMENSION(0:), INTENT(in)                       :: harvest_pool_bound !! The boundaries of the diameter classes
                                                                                       !! in the wood harvest pools
                                                                                       !! @tex $(m)$ @endtex
    REAL(r_std),DIMENSION(:,:),INTENT(in)                        :: veget_max          !! Maximum fraction of vegetation type 
    REAL(r_std), DIMENSION(:,:),INTENT(in)                       :: coppice_dens       !! The density of a coppice at the first
                                                                                       !! cutting.
                                                                                       !! @tex $( 1 / m**2 )$ @endtex 

    !! 0.2 Output variables


    !! 0.3 Modified variables   
    
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)             :: circ_class_biomass !! Biomass of the components of the model  
                                                                                       !! tree within a circumference
                                                                                       !! class @tex $(gC ind^{-1})$ @endtex  
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)                 :: circ_class_n       !! Number of individuals in each circ class
                                                                                       !! @tex $(m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)             :: circ_class_kill    !! Number of trees within a circ that needs
                                                                                       !! to be killed @tex $(ind m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)               :: bm_to_litter       !! Biomass transfer to litter 
                                                                                       !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)             :: harvest_pool       !! The wood and biomass that have been
                                                                                       !! havested by humans @tex $(gC)$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)                   :: harvest_type       !! Type of management that resulted
                                                                                       !! in the harvest (unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)                   :: harvest_cut        !! Type of cutting that was used for the harvest
                                                                                       !! (unitless)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)                 :: harvest_area       !! Harvested area (m^{2})
    INTEGER(i_std), DIMENSION(:,:),INTENT(inout)                 :: last_cut           !! Years since last thinning or clearcut (years)


    !! 0.4 Local variables

    REAL(r_std), DIMENSION(ncirc,nparts,nelements)               :: circ_class_biomass_new !! Biomass of the components of the model  
                                                                                       !! tree within a circumference
                                                                                       !! class @tex $(gC ind^{-1})$ @endtex  
    REAL(r_std), DIMENSION(ncirc)                                :: circ_class_n_new   !! Number of individuals in each circ class
                                                                                       !! @tex $(m^{-2})$ @endtex
    REAL(r_std), DIMENSION(ncirc)                                :: old_trees_left
    REAL(r_std)                                                  :: trees_needed
    REAL(r_std)                                                  :: qm_dia
    REAL(r_std)                                                  :: scale_factor
    INTEGER                                                      :: icir, jcir, imbc
    INTEGER(i_std)                                               :: iele, ipar         !! Indices
    REAL(r_std), DIMENSION(npts,nvm,nmbcomp,nelements)           :: check_intern       !! Contains the components of the internal
                                                                                       !! mass balance chech for this routine
                                                                                       !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)                   :: closure_intern     !! Check closure of internal mass balance
                                                                                       !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)                   :: pool_start         !! Start pool of this routine 
                                                                                       !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)                   :: pool_end           !! End pool of this routine 
                                                                                       !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(ncirc)                                :: exp_circ_class_n   !! temporary variable for expected circ class n
                                                                                       !! @tex $(ind m^{-2})$ @endtex 
    !_ ================================================================================================================================

    IF (printlev.GE.3) WRITE(numout,*) &
         'Entering second_coppice_cut sapiens_kill.f90'

        ! 1.2 Initialize check for mass balance closure
    IF (err_act.GT.1) THEN
       
       check_intern(:,:,:,:) = zero
       pool_start = zero
       DO iele = 1,nelements
          DO ipar = 1,nparts
             ! Initial litter and biomass. This routine is called
             ! within an ipts-LOOP and an ivm-LOOP so we test one 
             ! pixel and one pft at a time. Store in the first position.
             DO icir = 1,ncirc
                pool_start(1,1,iele) = pool_start(1,1,iele) + &
                     (circ_class_biomass(ipts,ivm,icir,ipar,iele) * &
                     circ_class_n(ipts,ivm,icir) * veget_max(ipts,ivm))
             ENDDO
             pool_start(1,1,iele) = pool_start(1,1,iele) + &
                  (bm_to_litter(ipts,ivm,ipar,iele) * veget_max(ipts,ivm))
          ENDDO

          ! The biomass harvest pool shouldn't be multiplied by 
          ! veget_max its units are already in gC pixel-1. Right 
          ! now the harvest pool is in gC.  We will divide by the
          ! area to make it per pixel. 
          pool_start(1,1,iele) = pool_start(1,1,iele) + &
               SUM(SUM(harvest_pool(ipts,ivm,:,iele,:),2)) / &
               (area(ipts) * contfrac(ipts)) 
          
       ENDDO 
       
    ENDIF ! err_act.GT.1

    ! Any stands that need to be coppiced for the second time?
    IF( SUM(circ_class_kill(ipts,ivm,:,ifm, icut)).GT. min_stomate ) THEN

       ! Debug
       IF(printlev_loc>=4)THEN
          WRITE(numout,*) 'Doing a second coppice cut'
          WRITE(numout,*) 'ipts,ivm,ifm,icut ',ipts,ivm,ifm,icut
       ENDIF
       !-

       ! Again, we harvest all aboveground biomass and leave the roots
       CALL harvest_aboveground_wood_living_roots(npts, ipts, ivm, &
            ifm, icut, bm_to_litter, circ_class_n, &
            circ_class_biomass, circ_class_kill, harvest_pool, &
            harvest_type, harvest_cut, harvest_area, &
            harvest_pool_bound, veget_max)

       ! reset the time since the last cut
       last_cut(ipts,ivm)=0

       ! We also reset circ_class_kill for this grid square, since
       ! we clearly don't need to do mortality if we have no more
       ! aboveground biomass left.
       circ_class_kill(ipts,ivm,:,:,:)=zero

       qm_dia = wood_to_qmdia(&
            circ_class_biomass(ipts,ivm,:,:,icarbon), &
            circ_class_n(ipts,ivm,:),ivm,pipe_tune2(ipts,ivm))

       ! qm_dia is the diameter at the base of the tree
       ! in meters.  This is the diameter that we will use
       ! to calculate the new distribution of trees among
       ! the circumference classes.
       exp_circ_class_n(:) = weibull_class_dist(qm_dia*m_to_cm, ivm, ifm)
       ! Now we have a certain number density and a certain amount of
       ! biomass.  We want to have the same number of trees that we had
       ! during the first cut, since the stumps from the first part aren't
       ! going to change (maybe a couple will die, but not a large difference).
       circ_class_n_new(:)=exp_circ_class_n(:)*SUM(circ_class_n(ipts,ivm,:))
       circ_class_biomass_new(:,:,:)=zero

       ! now we populate the new classes
       old_trees_left(:)=circ_class_n(ipts,ivm,:)
       DO icir=1,ncirc

          trees_needed=circ_class_n_new(icir)
          DO jcir=1,ncirc

             IF(trees_needed .LE. zero)EXIT

             IF(old_trees_left(jcir) .GE. trees_needed )THEN
                ! we can get all the trees we need from this class
                ! and don't have to continue searching
                circ_class_biomass_new(icir,:,:)=&
                     circ_class_biomass_new(icir,:,:)+&
                     circ_class_biomass(ipts,ivm,jcir,:,:)*trees_needed
                old_trees_left(jcir)=old_trees_left(jcir)-trees_needed
                trees_needed = zero
                EXIT
             ELSE

                ! the trees in this class are not sufficient, so we
                ! need to move all of them to the new class.
                circ_class_biomass_new(icir,:,:)=&
                     circ_class_biomass_new(icir,:,:)+&
                     circ_class_biomass(ipts,ivm,jcir,:,:)*&
                     old_trees_left(jcir)
                trees_needed = trees_needed -old_trees_left(jcir)
                old_trees_left(jcir)=zero

             ENDIF ! checking if we have enough trees

          ENDDO ! jcir=1,ncirc

       ENDDO ! icir=1,ncirc

       ! right now, circ_class_biomass_new gives the total biomass
       ! in each class, but it should only be for a model tree. 
       ! So let's normalize it.
       ! Guillaume M. -- A Weibull class can legitimately expect no tree, and the numerator
       ! is then zero as well: 0/0 yields NaN over that class's 9 parts x 2 elements.
       ! /!\ Divide whenever the count is strictly positive, and NEVER write a zero count:
       ! growth_fun_all divides b_inc_tot by circ_class_n class by class, guarded only on
       ! the SUM, so an exactly empty class returns NaN for the whole PFT the next day.
       DO icir=1,ncirc
          IF (circ_class_n_new(icir) .GT. zero) THEN
             circ_class_biomass_new(icir,:,:)=circ_class_biomass_new(icir,:,:)/&
                  circ_class_n_new(icir)
          ELSE
             circ_class_biomass_new(icir,:,:)=zero
          ENDIF
       ENDDO

       ! now update the variables that pass this information around
       DO icir=1,ncirc
          circ_class_biomass(ipts,ivm,icir,:,:)=&
               circ_class_biomass_new(icir,:,:)
          circ_class_n(ipts,ivm,icir)=circ_class_n_new(icir)
       ENDDO

       ! reset the time since the last cut
       last_cut(ipts,ivm)=0

       ! We also reset circ_class_kill for this grid square, since
       ! we clearly don't need to do mortality if we have no more
       ! aboveground biomass left.
       circ_class_kill(ipts,ivm,:,:,:)=zero

       ! After a first coppicing, we will regrow a certain number
       ! of shoots per stool.  For a second coppice cut, our starting
       ! density is exactly the same as the first cut.
       ! Guillaume M. -- Both divisions are unguarded: a stand left with no stem gives
       ! 0/0 here, and a zero factor turns every biomass into Inf below. The rescaling
       ! is mass-neutral by construction (n x biomass is invariant), so skipping it when
       ! there is nothing to rescale conserves mass exactly.
       IF (SUM(circ_class_n(ipts,ivm,:)) .GT. min_stomate) THEN
          scale_factor=coppice_dens(ipts,ivm)/SUM(circ_class_n(ipts,ivm,:))*&
               REAL(shoots_per_stool(ivm),r_std)
       ELSE
          scale_factor=zero
       ENDIF
       IF (scale_factor .GT. min_stomate) THEN
          circ_class_n(ipts,ivm,:)=circ_class_n(ipts,ivm,:)*&
               scale_factor

          ! in order not to change the total biomass, we need to
          ! make all our model trees smaller.
          circ_class_biomass(ipts,ivm,:,:,:)=circ_class_biomass(ipts,ivm,:,:,:)/&
               scale_factor
       ENDIF

    ENDIF ! any individuals to second cut coppice?

    !! Check numerical consistency of this routine
    IF (err_act.GT.1) THEN

       ! 4.3 Mass balance closure 
       pool_end(:,:,:) = zero
       DO iele = 1,nelements
          DO ipar = 1,nparts
             DO icir = 1,ncirc
                pool_end(1,1,iele) = pool_end(1,1,iele) + &
                     (circ_class_biomass(ipts,ivm,icir,ipar,iele) * &
                     circ_class_n(ipts,ivm,icir) * veget_max(ipts,ivm))
             ENDDO
             pool_end(1,1,iele) = pool_end(1,1,iele) + &
                  (bm_to_litter(ipts,ivm,ipar,iele) * veget_max(ipts,ivm))
          ENDDO

          ! The biomass harvest pool shouldn't be multiplied by veget_max
          ! its units are already in gC pixel-1. Right now the harvest 
          ! pool is in gC.  We will divide by the area to make it 
          ! per pixel.
          pool_end(1,1,iele) = pool_end(1,1,iele) + &
               SUM(SUM(harvest_pool(ipts,ivm,:,iele,:),2)) / &
               (area(ipts) * contfrac(ipts))
       ENDDO ! nelements
   
       ! 4.3.2 Calculate mass balance
       ! Common processes
       DO iele=1,nelements
          check_intern(1,1,ipoolchange,iele) = -un * (pool_end(1,1,iele) - &
               pool_start(1,1,iele)) 
       ENDDO

       closure_intern(:,:,:) = zero
       DO imbc = 1,nmbcomp
          DO iele=1,nelements
             ! Debug
             IF (printlev_loc>=4) WRITE(numout,*) &
                  'check_intern, ivm, imbc, iele, ', imbc, &
                  iele, check_intern(1,1,imbc,iele)
             !-
             closure_intern(1,1,iele) = closure_intern(1,1,iele) + &
                  check_intern(1,1,imbc,iele)
          ENDDO
       ENDDO

       ! 4.3.3 Check mass balance closure
       CALL check_mass_balance("second_coppice_cut", closure_intern, &
            npts, pool_end, pool_start, veget_max, 'ipts')
      
    ENDIF ! err_act.GT.1

    IF (printlev.GE.4) WRITE(numout,*) &
         'Leaving second_coppice_cut'

  END SUBROUTINE second_coppice_cut

!! ================================================================================================================================
!! SUBROUTINE   : first_coppice_cut
!!
!>\BRIEF        Transfer dead biomass to litter and wood to harvest pool.
!!
!! DESCRIPTION  : This is a total harvest of the aboveground biomass of a
!!                stand in such a way that the belowground wood continues to
!!                live.  This results in a sprouting in the following year of
!!                several shoots per stump.  Biomass is turned to litter, harvested,
!!                and left alive, depending on the pool.  The number of individuals
!!                goes up and the size of the model tree goes down.
!!
!!                NOTE: All branches and stems are harvested.  Basically, 
!!                      this means all aboveground woody biomass.  The fine roots
!!                      are killed and left to litter.  The sapwood and heartwood
!!                      belowground is left alive.  The leaves are killed, and the
!!                      labile and carbres pools are left untouched in the living
!!                      biomass.
!!
!! RECENT CHANGE(S): 
!!
!! MAIN OUTPUT VARIABLE(S): ::bm_to_litter, ::harvest_pool
!!
!! REFERENCE(S)   :
!!
!! FLOWCHART    : None
!!\n
!_ ================================================================================================================================
 
  SUBROUTINE first_coppice_cut (npts, ipts, ivm, ifm, &
       icut, harvest_pool_bound, veget_max, &
       circ_class_biomass, circ_class_n, circ_class_kill, &
       bm_to_litter, harvest_pool, harvest_type, harvest_cut, &
       harvest_area, last_cut)


    !! 0. Variable and parameter declaration

    !! 0.1 Input variables 
    ! NOTE: if harvest_pool_bound is not explicitely defined 
    ! the model uses 1:ndia_harvest+2 as its dimensions. If 
    ! you like to do so remember to change the counters in
    ! this routine
    INTEGER(i_std), INTENT(in)                                   :: npts               !! Domain size (-)
    INTEGER(i_std), INTENT(in)                                   :: ipts               !! current grid square
    INTEGER(i_std), INTENT(in)                                   :: ivm                !! current PFT
    INTEGER(i_std), INTENT(in)                                   :: ifm                !! current fm strategy
    INTEGER(i_std), INTENT(in)                                   :: icut               !! current cut time
    REAL(r_std), DIMENSION(0:), INTENT(in)                       :: harvest_pool_bound !! The boundaries of the diameter classes
                                                                                       !! in the wood harvest pools
                                                                                       !! @tex $(m)$ @endtex
    REAL(r_std),DIMENSION(:,:),INTENT(in)                        :: veget_max          !! Maximum fraction of vegetation type 


    !! 0.2 Output variables


    !! 0.3 Modified variables   

    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)             :: circ_class_biomass !! Biomass of the components of the model  
                                                                                       !! tree within a circumference
                                                                                       !! class @tex $(gC ind^{-1})$ @endtex  
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)                 :: circ_class_n       !! Number of individuals in each circ class
                                                                                       !! @tex $(m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)             :: circ_class_kill    !! Number of trees within a circ that needs
                                                                                       !! to be killed @tex $(ind m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)               :: bm_to_litter       !! Biomass transfer to litter 
                                                                                       !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)             :: harvest_pool       !! The wood and biomass that have been
                                                                                       !! havested by humans @tex $(gC)$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)                   :: harvest_type       !! Type of management that resulted
                                                                                       !! in the harvest (unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)                   :: harvest_cut        !! Type of cutting that was used for the harvest
                                                                                       !! (unitless)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)                 :: harvest_area       !! Harvested area (m^{2})
    INTEGER(i_std), DIMENSION(:,:),INTENT(inout)                 :: last_cut           !! Years since last thinning or clearcut (years)

    !! 0.4 Local variables
    INTEGER(i_std)                                               :: imbc, iele         !! Indices
    INTEGER(i_std)                                               :: ipar, icir         !! Indices
    REAL(r_std), DIMENSION(npts,nvm,nmbcomp,nelements)           :: check_intern       !! Contains the components of the internal
                                                                                       !! mass balance chech for this routine
                                                                                       !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)                   :: closure_intern     !! Check closure of internal mass balance
                                                                                       !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)                   :: pool_start         !! Start pool of this routine 
                                                                                       !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)                   :: pool_end           !! End pool of this routine 
                                                                                       !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex

    !_ ================================================================================================================================

    IF (printlev.GE.3) WRITE(numout,*) &
         'Entering first_coppice_cut sapiens_kill.f90'

    ! 1.2 Initialize check for mass balance closure
    IF (err_act.GT.1) THEN
       
       check_intern(:,:,:,:) = zero
       pool_start = zero
       DO iele = 1,nelements
          DO ipar = 1,nparts
             ! Initial litter and biomass. This routine is called
             ! within an ipts-LOOP and an ivm-LOOP so we test one 
             ! pixel and one pft at a time. Store in the first position.
             DO icir = 1,ncirc
                pool_start(1,1,iele) = pool_start(1,1,iele) + &
                     (circ_class_biomass(ipts,ivm,icir,ipar,iele) * &
                     circ_class_n(ipts,ivm,icir) * veget_max(ipts,ivm))
             ENDDO
             pool_start(1,1,iele) = pool_start(1,1,iele) + &
                  (bm_to_litter(ipts,ivm,ipar,iele) * veget_max(ipts,ivm))
          ENDDO

          ! The biomass harvest pool shouldn't be multiplied by 
          ! veget_max its units are already in gC pixel-1. Right 
          ! now the harvest pool is in gC.  We will divide by the
          ! area to make it per pixel. 
          pool_start(1,1,iele) = pool_start(1,1,iele) + &
               SUM(SUM(harvest_pool(ipts,ivm,:,iele,:),2)) / &
               (area(ipts) * contfrac(ipts))
          
       ENDDO 
       
    ENDIF ! err_act.GT.1

    ! What if we coppice for the first cut?
    IF( SUM(circ_class_kill(ipts,ivm,:,ifm, icut))&
         .GT. min_stomate ) THEN
       
       ! Debug
       IF(printlev>=4)THEN
          WRITE(numout,*) 'Doing a first coppice cut'
          WRITE(numout,*) 'ipts,ivm,ifm,icut ',ipts,ivm,ifm,icut
       ENDIF
       !-

       CALL harvest_aboveground_wood_living_roots(npts, ipts, ivm, &
            ifm, icut, bm_to_litter, circ_class_n, &
            circ_class_biomass, circ_class_kill, harvest_pool, &
            harvest_type, harvest_cut, harvest_area, &
            harvest_pool_bound, veget_max)

       ! reset the time since the last cut
       last_cut(ipts,ivm)=0

       ! We also reset circ_class_kill for this grid square, since
       ! we clearly don't need to do mortality if we have no more
       ! aboveground biomass left.
       circ_class_kill(ipts,ivm,:,:,:)=zero

       ! After a first coppicing, we will regrow a certain number
       ! of shoots per stool.  Therefore, our density starting
       ! next year will be much higher.
       circ_class_n(ipts,ivm,:)=circ_class_n(ipts,ivm,:)*&
            REAL(shoots_per_stool(ivm),r_std)

       ! in order not to change the total biomass, we need to
       ! make all our model trees smaller.
       circ_class_biomass(ipts,ivm,:,:,:)=circ_class_biomass(ipts,ivm,:,:,:)/&
            REAL(shoots_per_stool(ivm),r_std)

    ENDIF ! any trees to first cut coppice?

    !! Check numerical consistency of this routine

    IF (err_act.GT.1) THEN

       ! 4.3 Mass balance closure 
       pool_end(:,:,:) = zero
       DO iele = 1,nelements
          DO ipar = 1,nparts
             DO icir = 1,ncirc
                pool_end(1,1,iele) = pool_end(1,1,iele) + &
                     (circ_class_biomass(ipts,ivm,icir,ipar,iele) * &
                     circ_class_n(ipts,ivm,icir) * veget_max(ipts,ivm))
             ENDDO
             pool_end(1,1,iele) = pool_end(1,1,iele) + &
                  (bm_to_litter(ipts,ivm,ipar,iele) * veget_max(ipts,ivm))
          ENDDO

          ! The biomass harvest pool shouldn't be multiplied by veget_max
          ! its units are already in gC pixel-1. Right now the harvest 
          ! pool is in gC.  We will divide by the area to make it 
          ! per pixel.
          pool_end(1,1,iele) = pool_end(1,1,iele) + &
               SUM(SUM(harvest_pool(ipts,ivm,:,iele,:),2)) / &
               (area(ipts) * contfrac(ipts))
       ENDDO ! nelements
   
       ! 4.3.2 Calculate mass balance
       ! Common processes
       DO iele=1,nelements
          check_intern(1,1,ipoolchange,iele) = -un * (pool_end(1,1,iele) - &
               pool_start(1,1,iele)) 
       ENDDO

       closure_intern(:,:,:) = zero
       DO imbc = 1,nmbcomp
          DO iele=1,nelements
             ! Debug
             IF (printlev_loc>=4) WRITE(numout,*) &
                  'check_intern, ivm, imbc, iele, ', imbc, &
                  iele, check_intern(1,1,imbc,iele)
             !-
             closure_intern(1,1,iele) = closure_intern(1,1,iele) + &
                  check_intern(1,1,imbc,iele)
          ENDDO
       ENDDO

       ! 4.3.3 Check mass balance closure
       CALL check_mass_balance("first_coppice_cut", closure_intern, &
            npts, pool_end, pool_start, veget_max, 'ipts')
      
    ENDIF ! err_act.GT.1

    IF (printlev.GE.4) WRITE(numout,*) &
         'Leaving first_coppice_cut'

  END SUBROUTINE first_coppice_cut

!! ================================================================================================================================
!! SUBROUTINE   : clearcut_harvest_aboveground_dead_roots
!!
!>\BRIEF        Transfer dead biomass to litter and wood to harvest pool 
!!              for a whole stand.
!!
!! DESCRIPTION  : Moves dead biomass to the litter pool and the harvest pool
!!                where appropriate, assuming the whole stand is harvested.  Various
!!                counters are reset.
!!
!!                NOTE: If harvest_fraction is 1, all the aboveground biomass is
!!                      harvested.  If it's 1-branch_ratio, only stems are 
!!                      harvested.  Belowground wood, roots, fruits, and
!!                      leaves all become litter.  A portion of aboveground
!!                      wood is put in the harvest pool.  Reserves and labile
!!                      carbon are assumed to be distributed evenly across
!!                      all non-heartwood components.
!!
!! RECENT CHANGE(S): 
!!
!! MAIN OUTPUT VARIABLE(S): ::bm_to_litter, ::harvest_pool
!!
!! REFERENCE(S)   :
!!
!! FLOWCHART    : None
!!\n
!_ ================================================================================================================================
 
  SUBROUTINE clearcut_harvest_aboveground_dead_roots (npts, ipts, ivm, &
       ifm, icut, bm_to_litter, woody_litter_by_cut, &
       circ_class_biomass, circ_class_kill, circ_class_n, harvest_pool, &
       harvest_type, harvest_cut, harvest_area, &
       harvest_pool_bound, veget_max, last_cut, age_stand, &
       mai_count, harvest_fraction)


    !! 0. Variable and parameter declaration

    !! 0.1 Input variables 
    ! NOTE: if harvest_pool_bound is not explicitely defined 
    ! the model uses 1:ndia_harvest+2 as its dimensions. If 
    ! you like to do so remember to change the counters in
    ! this routine
    INTEGER(i_std), INTENT(in)                                   :: npts               !! Domain size (-)
    INTEGER(i_std), INTENT(in)                                   :: ipts               !! current grid square
    INTEGER(i_std), INTENT(in)                                   :: ivm                !! current PFT
    INTEGER(i_std), INTENT(in)                                   :: ifm                !! current fm strategy
    INTEGER(i_std), INTENT(in)                                   :: icut               !! current cut time
    REAL(r_std), DIMENSION(0:), INTENT(in)                       :: harvest_pool_bound !! The boundaries of the diameter classes
                                                                                       !! in the wood harvest pools
                                                                                       !! @tex $(m)$ @endtex
    REAL(r_std),DIMENSION(:,:),INTENT(in)                        :: veget_max          !! Maximum fraction of vegetation type 
    REAL(r_std),INTENT(in)                                       :: harvest_fraction   !! The fraction of the aboveground
                                                                                       !! biomass that is harvested.  This is
                                                                                       !! generally 1-branch_ratio if you only
                                                                                       !! want to harvest the stems, and 1 if
                                                                                       !! you want to harvest everything.


    !! 0.2 Output variables


    !! 0.3 Modified variables   
   
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)             :: circ_class_biomass !! Biomass of the components of the model  
                                                                                       !! tree within a circumference
                                                                                       !! class @tex $(gC ind^{-1})$ @endtex  
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)                 :: circ_class_n       !! Number of individuals in each circ class
                                                                                       !! @tex $(m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)             :: circ_class_kill    !! Number of trees within a circ that needs
                                                                                       !! to be killed @tex $(ind m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)               :: bm_to_litter       !! Biomass transfer to litter 
                                                                                       !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)                 :: woody_litter_by_cut!! Saved woody litter by icut not by iparts

    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)             :: harvest_pool       !! The wood and biomass that have been
                                                                                       !! havested by humans @tex $(gC)$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)                   :: harvest_type       !! Type of management that reulted
                                                                                       !! in the harvest (unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)                   :: harvest_cut        !! Type of cutting that was used for the harvest
                                                                                       !! (unitless)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)                 :: harvest_area       !! Harvested area (m^{2})
    INTEGER(i_std), DIMENSION(:,:),INTENT(inout)                 :: last_cut           !! Years since last thinning or clearcut (years)
    INTEGER(i_std), DIMENSION(:,:),INTENT(inout)                 :: age_stand          !! Age of stand (years)
    INTEGER(i_std), DIMENSION(:,:),INTENT(inout)                 :: mai_count          !! The number of times we've
                                                                                       !! calculated the volume increment
                                                                                       !! for a stand
    !! 0.4 Local variables
    INTEGER(i_std)                                               :: imbc, iele         !! Indices
    INTEGER(i_std)                                               :: ipar, icir         !! Indices
    REAL(r_std), DIMENSION(npts,nvm,nmbcomp,nelements)           :: check_intern       !! Contains the components of the internal
                                                                                       !! mass balance chech for this routine
                                                                                       !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)                   :: closure_intern     !! Check closure of internal mass balance
                                                                                       !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)                   :: pool_start         !! Start pool of this routine 
                                                                                       !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)                   :: pool_end           !! End pool of this routine 
                                                                                       !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    
    !_ ================================================================================================================================

    IF (printlev.GE.3) WRITE(numout,*) &
         'Entering clearcut_harvest_aboveground_dead_roots sapiens_kill.f90'

        ! 1.2 Initialize check for mass balance closure
    IF (err_act.GT.1) THEN
       
       check_intern(:,:,:,:) = zero
       pool_start = zero
       DO iele = 1,nelements
          DO ipar = 1,nparts
             ! Initial litter and biomass. This routine is called
             ! within an ipts-LOOP and an ivm-LOOP so we test one 
             ! pixel and one pft at a time. Store in the first position.
             DO icir = 1,ncirc
                pool_start(1,1,iele) = pool_start(1,1,iele) + &
                     (circ_class_biomass(ipts,ivm,icir,ipar,iele) * &
                     circ_class_n(ipts,ivm,icir) * veget_max(ipts,ivm))
             ENDDO
             pool_start(1,1,iele) = pool_start(1,1,iele) + &
                  (bm_to_litter(ipts,ivm,ipar,iele) * veget_max(ipts,ivm))
          ENDDO

          ! The biomass harvest pool shouldn't be multiplied by 
          ! veget_max its units are already in gC pixel-1. Right 
          ! now the harvest pool is in gC.  We will divide by the
          ! area to make it per pixel. 
          pool_start(1,1,iele) = pool_start(1,1,iele) + &
               SUM(SUM(harvest_pool(ipts,ivm,:,iele,:),2)) / &
               (area(ipts) * contfrac(ipts))
          
       ENDDO 
       
    ENDIF ! err_act.GT.1

    IF( SUM(circ_class_kill(ipts,ivm,:,ifm, icut)) .GT. min_stomate ) THEN

       CALL harvest_aboveground_dead_roots(npts, ipts, ivm, &
            ifm, icut, bm_to_litter, woody_litter_by_cut, &
            circ_class_biomass, circ_class_kill, circ_class_n, harvest_pool, &
            harvest_type, harvest_cut, harvest_area, &
            harvest_pool_bound, veget_max, &
            harvest_fraction)

       ! We have to reset the age of the stand and the time since the
       ! last human intervention.
       age_stand(ipts,ivm)=0
       last_cut(ipts,ivm)=0
       mai_count(ipts,ivm)=0

       ! We also reset circ_class_kill for this grid square, since
       ! we clearly don't need to do mortality if we have no more
       ! trees left.
       circ_class_kill(ipts,ivm,:,ifm,icut)=zero

    ENDIF

    !! Check numerical consistency of this routine

    IF (err_act.GT.1) THEN

       ! 4.3 Mass balance closure 
       pool_end(:,:,:) = zero
       DO iele = 1,nelements
          DO ipar = 1,nparts
             DO icir = 1,ncirc
                pool_end(1,1,iele) = pool_end(1,1,iele) + &
                     (circ_class_biomass(ipts,ivm,icir,ipar,iele) * &
                     circ_class_n(ipts,ivm,icir) * veget_max(ipts,ivm))
             ENDDO
             pool_end(1,1,iele) = pool_end(1,1,iele) + &
                  (bm_to_litter(ipts,ivm,ipar,iele) * veget_max(ipts,ivm))
          ENDDO

          ! The biomass harvest pool shouldn't be multiplied by veget_max
          ! its units are already in gC pixel-1. Right now the harvest 
          ! pool is in gC.  We will divide by the area to make it 
          ! per pixel.
          pool_end(1,1,iele) = pool_end(1,1,iele) + &
               SUM(SUM(harvest_pool(ipts,ivm,:,iele,:),2)) / &
               (area(ipts) * contfrac(ipts))
       ENDDO ! nelements
   
       ! 4.3.2 Calculate mass balance
       ! Common processes
       DO iele=1,nelements
          check_intern(1,1,ipoolchange,iele) = -un * (pool_end(1,1,iele) - &
               pool_start(1,1,iele)) 
       ENDDO

       closure_intern(:,:,:) = zero
       DO imbc = 1,nmbcomp
          DO iele=1,nelements
             ! Debug
             IF (printlev_loc>=4) WRITE(numout,*) &
                  'check_intern, ivm, imbc, iele, ', imbc, &
                  iele, check_intern(1,1,imbc,iele)
             !-
             closure_intern(1,1,iele) = closure_intern(1,1,iele) + &
                  check_intern(1,1,imbc,iele)
          ENDDO
       ENDDO

       ! 4.3.3 Check mass balance closure
       CALL check_mass_balance("clearcut_harvest_aboveground_dead_roots", &
            closure_intern, npts, pool_end, pool_start, veget_max, 'ipts')
      
    ENDIF ! err_act.GT.1

    IF (printlev.GE.4) WRITE(numout,*) &
         'Leaving clearcut_harvest_aboveground_dead_roots'

  END SUBROUTINE clearcut_harvest_aboveground_dead_roots

END MODULE sapiens_kill
