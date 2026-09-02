! =================================================================================================================================
! MODULE       : stomate_kill
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       Simulate mortality of individuals and update biomass, litter and 
!! stand density of the PFT
!!
!!\n DESCRIPTION : Simulate mortality of individuals and update biomass, litter and 
!! stand density of the PFT. This module replaces lpj_gap/kill.  Biomass actually
!! dies here and is moved to the litter pools.  This does not work with the DGVM
!! at the moment.
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
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_stomate/stomate_kill.f90 $
!! $Date: 2026-05-14 11:18:13 +0200 (jeu. 14 mai 2026) $
!! $Revision: 9542 $
!! \n
!_ ================================================================================================================================

MODULE stomate_kill

  ! modules used:

  USE stomate_data
  USE pft_parameters
  USE constantes
  USE ioipsl_para 
  USE stomate_prescribe
  USE stomate_mark_kill,      ONLY : distribute_mortality_biomass
  USE function_library,       ONLY : nmax, wood_to_dia, wood_to_qmdia, check_vegetation_area, &
                                     check_mass_balance, sort_circ_class_biomass, &
                                     weibull_class_dist, &
                                     cc_to_biomass, get_printlev
  USE sapiens_lcchange,       ONLY : merge_biomass_pfts, calculate_shares, move_pft_properties

  IMPLICIT NONE

  ! private & public routines

  PRIVATE
  PUBLIC natural_mortality, natural_mortality_clear, mortality_clean, disturbance_kill
  
  ! Variable declaration 

  LOGICAL, SAVE             :: firstcall_stomate_kill = .TRUE.     !! first call flag
!$OMP THREADPRIVATE(firstcall_stomate_kill)
  INTEGER(i_std), SAVE      :: printlev_loc                        !! Local level of text output for current module
!$OMP THREADPRIVATE(printlev_loc)
CONTAINS


!! ================================================================================================================================
!! SUBROUTINE   : natural_mortality_clear
!!
!>\BRIEF        Set the firstcall flag back to .TRUE. to prepare for the next simulation.
!_ ================================================================================================================================
  
  SUBROUTINE natural_mortality_clear
     firstcall_stomate_kill = .TRUE.
  END SUBROUTINE natural_mortality_clear


!! ================================================================================================================================
!! SUBROUTINE   : natural_mortality
!!
!>\BRIEF        Transfer dead biomass to litter and update stand density for trees
!!              which die of natural causes.
!!
!! DESCRIPTION  : This routine has a simple purpose: kill individuals by natural causes.  
!!   The variable circ_class_kill indicates how many individuals need to be killed by the various
!!   processes in the code, and is calculated in other modules. natural_mortality takes
!!   this variable and decides on the correct order of which to actually kill the
!!   trees, making sure that we don't double count.  
!!
!! RECENT CHANGE(S): 
!!
!! MAIN OUTPUT VARIABLE(S): ::circ_class_biomass; ::circ_class_n density of individuals, 
!! ::mortality mortality (fraction of trees that is dying per time step)
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
 
  SUBROUTINE natural_mortality (npts,  bm_to_litter, woody_litter_by_cut, &
       circ_class_biomass, circ_class_kill, circ_class_n,emissions_fire, &
       veget_max, biomass_cut, plant_status, tot_res_target, &
       tot_lab_target, forest_managed)

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables 
    INTEGER(i_std), INTENT(in)                                   :: npts               !! Domain size (-)
    REAL(r_std),DIMENSION(:,:),INTENT(in)                        :: veget_max          !! Maximum fraction of vegetation type 
    INTEGER(i_std),DIMENSION(:,:),INTENT(in)                     :: forest_managed     !! forest management flag



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
    REAL(r_std), DIMENSION(:,:,:,:,:,:), INTENT(inout)           :: biomass_cut        !! Biomass change due to  ncut
    REAL(r_std), DIMENSION(:,:), INTENT(inout)                   :: plant_status       !! Growth and phenological status of the plant 
                                                                                       !! see constants voor values
    REAL(r_std), DIMENSION(:,:,:),INTENT(inout)                  :: tot_res_target     !! Target for reserve carbon pool (only for grassland) 
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(:,:,:),INTENT(inout)                  :: tot_lab_target     !! Target for labile carbon pool (only for grassland) 
                                                                                       !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:),INTENT(inout)                :: emissions_fire     !! Emissions from fire. (gC m^{-2})


    !! 0.4 Local variables
    INTEGER(i_std)                                               :: ipts, ivm, ipar    !! Indices
    INTEGER(i_std)                                               :: iele, icir, imbc   !! Indices
    INTEGER(i_std)                                               :: ifm,icut,ipart     !! Indices
    INTEGER,DIMENSION(nvm)                                       :: nkilled            !! the number of grid points at which a given
                                                                                       !! given PFT's biomass has been reduced to zero
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
    REAL(r_std), DIMENSION(nparts,nelements)                     :: vartmp_bm          !! Temporary variable
    REAL(r_std)                                                  :: ndying             !! Number of trees that die in a single day (trees/m2)
    REAL(r_std), DIMENSION(npts,nvm,ncirc)                       :: old_circ_class_n   !! Grass density before it is adjusted to the
                                                                                       !! available resevre and labile carbon (ind m-2)
    REAL(r_std)                                                  :: Trl0               !! The sum of individual reserve and labile carbon before 
                                                                                       !! density adjustment (gC ind-1)
    REAL(r_std)                                                  :: Ttarget_area       !! The sum of areal target for reserve and labile carbon (gC m-2)
    REAL(r_std)                                                  :: Cother             !! The sum of carbon from all compartments except for reserve and
                                                                                       !! labile pool (gC ind-1)
    REAL(r_std), DIMENSION(npts,nvm,ncirc,nparts,nelements)      :: old_circ_class_biomass !! biomass components before the density adjustment (gC ind-1)
    REAL(r_std)                                                  :: ratio_den          !! Ratio of old density to new density (Unitless)
    REAL(r_std)                                                  :: ratio_res          !! Ratio of reserve carbon to the sum of reserve and labile carbon
                                                                                       !! (Unitless)
    REAL(r_std)                                                  :: ratio_lab          !! Ratio of labile carbon to the sum of reserve and labile carbon
                                                                                       !! (Unitless)
    LOGICAL, DIMENSION(npts,nvm)                                 :: has_subtracted     !! for litter pool
    REAL(r_std),PARAMETER                                        :: scale_dec = 0      !! Scaling factor to calculate reserve and labile carbon
                                                                                       !! when decreasing grassland density (unitless)
                                                                                       !! Has been used to tune the results. Could be externalized 
                                                                                       !! and made PFT-specific when found to be impactful.
    REAL(r_std),PARAMETER                                        :: scale_carbon = 0   !! Scaling factor to control the percentage of carbon into litter
                                                                                       !! pool when decreasing grassland density (unitless)
                                                                                       !! Has been used to tune the results. Could be externalized 
                                                                                       !! and made PFT-specific when found to be impactful.
    REAL(r_std)                                                  :: temp_dens_target   !! Temporary and local target density.
!_ ================================================================================================================================

    !! Initialize variables

    !! 1.1 Set firstcall flag
    IF ( firstcall_stomate_kill ) THEN
       !! Initialize local printlev
      ! Guillaume M. -- get_printlev calls getin_p, a COLLECTIVE bcast. Called here,
      ! inside the compute loop under a local firstcall, it hangs the run as soon as the
      ! ranks stop reaching it together. Use printlev, read once at startup by every rank.
       printlev_loc=printlev
       firstcall_stomate_kill = .FALSE.
    ENDIF

    IF (printlev_loc.GE.2) WRITE(numout,*) 'Entering natural mortality. Use constant mortality = ', &
         ok_constant_mortality

    !! 1.2 Initialize variables
    nkilled(:)=0
    has_subtracted(:,:)= .FALSE.

    ! 1.3 Initialize check for mass balance closure
    IF (err_act.GT.1) THEN

       check_intern(:,:,:,:) = zero
       pool_start(:,:,:) = zero
       DO ipar = 1,nparts
          DO iele = 1,nelements        
             !  Initial litter and biomass
             DO icir = 1,ncirc
                pool_start(:,:,iele) = pool_start(:,:,iele) + &
                     (circ_class_biomass(:,:,icir,ipar,iele) * &
                     circ_class_n(:,:,icir) * veget_max(:,:))

             ENDDO
             pool_start(:,:,iele) = pool_start(:,:,iele) + &
                  bm_to_litter(:,:,ipar,iele) * veget_max(:,:)
          ENDDO
       ENDDO

       !! 1.3 Initialize check for area conservation
       veget_max_begin(:,:) = veget_max(:,:)

    ENDIF ! err_act.GT.1

    ! The total biomass at the start of the routine
    init_biomass = cc_to_biomass(npts,nvm,circ_class_biomass(:,:,:,:,:),&
         circ_class_n(:,:,:))


    DO  ipts = 1,npts ! loop over land_points

       DO ivm = 2,nvm  ! loop over #PFT

          IF(veget_max(ipts,ivm) == zero)THEN
             ! this vegetation type is not present, so no reason 
             ! to do the calculation
             CYCLE
          ENDIF

          IF (is_tree(ivm)) THEN
             ! Select the target density as a function of forest
             ! management
             temp_dens_target = dens_target(ivm,forest_managed(ipts,ivm))
          ELSE
             ! Prescribe a value that enables passing the loops
             ! below. For grasses and croplands, recruitment is
             ! not used. 
             temp_dens_target = 1.5
          END IF
          
          !! Make unmanaged forest die slowly
          ! This IF statement is used to overcome a huge wood litter input 
          ! when a forest has to be killed when it reach dens_target(ivm). 
          ! Instead, the forest will slowly die at a rate define by 365*15 
          ! (days. Note that 15 is externalized as ndying_year) in order 
          ! to distribute the wood litter over multiple years. During this 
          ! state, trees recruitment is disabled.  
          ! Guillaume M. -- Slow death of understocked stands, on a DEDICATED flag, distinct
          ! from ok_min_density_reset: this block is not a reset but a progressive mortality
          ! sink (circ_class_n -= ndying, below). Keeping the two flags separate is what makes
          ! the class-1 RDI response attributable to this sink rather than to establishment.
          IF(ok_slow_death_density .AND. &
             SUM(circ_class_n(ipts,ivm,:)) .LT. temp_dens_target*ha_to_m2) THEN

             ! If we force the clear cut, we should not make use of the slow 
             ! mortality approach
             IF (.NOT. forced_clear_cut) THEN

                DO icir=1,ncirc

                   ! Number of tree to be removed every days since dens_target
                   ! is exceeded. The first time ndying is calculated, the 
                   ! value will result in circ_calss_n = 0 after ndying years. 
                   ! But ndying_years is long and trees might get killed for 
                   ! other reasons as well so circ_class_n might reach zero 
                   ! in less than ndying_year. Make ndying never exceeds
                   ! circ_class_n(:,:,icir).
                   ndying = MIN(circ_class_n(ipts,ivm,icir), &
                           (temp_dens_target/(365*ndying_year(ivm)))*ha_to_m2)

                   ! Move the dead trees to the bm_to_litter pool
                   bm_to_litter(ipts,ivm,:,:) = bm_to_litter(ipts,ivm,:,:) + &
                        circ_class_biomass(ipts,ivm,icir,:,:)*ndying

                   ! Save woody litter by icut. This will be used in pest module
                   ! (stomate_pest) for one of the susceptibility indice.
                   ! Since this process clears the stand it will be saved to icut_clear
                   woody_litter_by_cut(ipts,ivm,icut_clear) = &
                        woody_litter_by_cut(ipts,ivm,icut_clear) + &
                        ( circ_class_biomass(ipts,ivm,icir,isapabove,icarbon) + &
                        circ_class_biomass(ipts,ivm,icir,iheartabove,icarbon) )* ndying

                   ! remove the number of individuals that died
                   circ_class_n(ipts,ivm,icir) = circ_class_n(ipts,ivm,icir)-&
                        ndying

                END DO

             END IF ! .NOT. forced_clear_cut

             ! Debug
             IF (printlev_loc>=4) THEN
                WRITE(numout,*) 'Slowly kill PFT, ',ipts,ivm
                WRITE(numout,*) 'ndying, ', ndying
                WRITE(numout,*) 'circ_class_n, ', SUM(circ_class_n(ipts,ivm,:))
                WRITE(numout,*) 'nb strem thresold,', temp_dens_target*0.01*ha_to_m2
             ENDIF
             !-

             ! When only 1% of the dens_target tree remain we kill the 
             ! rest of the stand. If we force the clear cut we do it in a single
             ! time step
             IF(forced_clear_cut .OR. & 
                     SUM(circ_class_n(ipts,ivm,:)) .LE. temp_dens_target*0.01*ha_to_m2)THEN

                DO icir = 1,ncirc

                   ! Move the dead trees to the bm_to_litter pool
                   bm_to_litter(ipts,ivm,:,:) = bm_to_litter(ipts,ivm,:,:)+ &
                        circ_class_biomass(ipts,ivm,icir,:,:)*&
                        circ_class_n(ipts,ivm,icir)

                   ! Save woody litter by icut. This will be used in the pest module
                   ! (stomate_pest) for one of the susceptibility indice.
                   ! Since this process clears the stand it will be saved to
                   ! icut_clear
                   woody_litter_by_cut(ipts,ivm,icut_clear) = &
                        woody_litter_by_cut(ipts,ivm,icut_clear) + &
                        ( circ_class_biomass(ipts,ivm,icir,isapabove,icarbon) + &
                        circ_class_biomass(ipts,ivm,icir,iheartabove,icarbon) ) * &
                        circ_class_n(ipts,ivm,icir)
 
                   ! zero the number of individuals in this circ class
                   circ_class_n(ipts,ivm,icir) = zero
                   circ_class_biomass(ipts,ivm,icir,:,:) = zero

                END DO

                ! No plants left. Set plant_status to idead
                plant_status(ipts,ivm) = idead

                ! The PFT is dead so if we already assigned some trees
                ! to be killed in the rest of this module, this is no
                ! longer possible. Set circ_class_kill to zero.
                circ_class_kill(ipts,ivm,:,:,:) = zero

                ! Debug
                IF (printlev_loc>=4) THEN
                   WRITE(numout,*) 'End of slow killing, ', ipts,ivm
                   WRITE(numout,*) 'circ_class_n, ', SUM(circ_class_n(ipts,ivm,:))
                   WRITE(numout,*) 'Plant_status, ',plant_status(ipts,ivm)
                ENDIF
                !-

             ENDIF

             ! Actual biomass after mortality has been accounted for
             actu_biomass = cc_to_biomass(npts,nvm,&
                  circ_class_biomass(ipts,ivm,:,:,:),&
                  circ_class_n(ipts,ivm,:))

             biomass_cut(ipts,ivm,:,:,ifm_none,icut_clear) = &
                  biomass_cut(ipts,ivm,:,:,ifm_none,icut_clear) + &
                  init_biomass(ipts,ivm,:,:) - actu_biomass(:,:)

             ! Update the initial biomass to avoid double counting
             init_biomass(ipts,ivm,:,:) = actu_biomass(:,:)

          ENDIF

          !! 3. Remove individuals that were send to circ_class_kill

          ! First, all the plants that are supposed to be killed by 
          ! fire are killed. We do not yet know how to handle this, 
          ! so this will have to be taken into account when SPITFIRE 
          ! is coupled.

          ! All of the natural death is grouped 
          ! together in one pool since all of the biomass will be left on 
          ! site (moved to the litter pools).
          DO ifm=1,nfm_types

             DO icut=1,ncut_times

                ! Since we gave circ_class_kill the same dimensions as the
                ! harvest pool, we need to explicitly say which combinations
                ! of ifm and icut are killed in which ways.  Currently, I am
                ! putting all natural death into ifm_none, even if it happens
                ! in another management type (for example, self-thinning will 
                ! happen with ifm_thin, but it's really a natural death).
                IF(ifm .NE. ifm_none)CYCLE
                IF((icut .NE. icut_thin) .AND. &
                     (icut .NE. icut_clear) .AND. &
                     (icut .NE. icut_beetle) .AND. &
                     (icut .NE. icut_storm_break) .AND. &
                     (icut .NE. icut_fire) .AND. &
                     (icut .NE. icut_storm_uproot)) CYCLE

                ! Killing is done a little differently for trees and 
                ! non-trees, since circ_classes are not defined for non-trees.
                IF(is_tree(ivm))THEN

                   DO icir=1,ncirc

                      IF (circ_class_kill(ipts,ivm,icir,ifm,icut).EQ.zero) CYCLE
                      
                      ! Debug
                      IF(test_pft == ivm .AND. test_grid == ipts .AND. &
                           printlev_loc>=4) THEN
                         WRITE(numout,*) 'killing before: ',ipts,ivm,icir,ifm,icut,&
                              circ_class_kill(ipts,ivm,icir,ifm,icut),&
                              circ_class_n(ipts,ivm,icir)
                      ENDIF
                      !-
                      
                      ! Consistency checking and subsequent correction
                      IF(circ_class_kill(ipts,ivm,icir,ifm,icut) .GT. &
                                circ_class_n(ipts,ivm,icir)) THEN
                           circ_class_kill(ipts,ivm,icir,ifm,icut) = circ_class_n(ipts,ivm,icir)
                      ENDIF                      

                      ! move the dead biomass to the respective litter pool
                      vartmp_bm(:,:) = circ_class_biomass(ipts,ivm,icir,:,:) * &
                                       circ_class_kill(ipts,ivm,icir,ifm,icut)
                      IF (icut .NE. icut_fire) THEN
                         bm_to_litter(ipts,ivm,:,:) = bm_to_litter(ipts,ivm,:,:) + vartmp_bm(:,:)
                      ELSE
                         ! part of dead trees are assumed being burned through live crown burning.
                         ! Combusted biomass components include all leaf and fruit and biomass being
                         ! considered as 1-hour fuel (i.e., small branches). A more refined
                         ! approach could be introduced later on.
                         emissions_fire(ipts,ivm,:,ifirecrown) = emissions_fire(ipts,ivm,:,ifirecrown) + &
                             vartmp_bm(ileaf,:) + vartmp_bm(ifruit,:) + &
                             ( vartmp_bm(isapabove,:) + vartmp_bm(iheartabove,:) + &
                               vartmp_bm(icarbres,:) + vartmp_bm(ilabile,:))*alloc_firefuel(ihour1)

                         bm_to_litter(ipts,ivm,isapabove,:) = bm_to_litter(ipts,ivm,isapabove,:) + &
                                    vartmp_bm(isapabove,:)*(1-alloc_firefuel(ihour1))
                         bm_to_litter(ipts,ivm,iheartabove,:) = bm_to_litter(ipts,ivm,iheartabove,:) + &
                                    vartmp_bm(iheartabove,:)*(1-alloc_firefuel(ihour1))
                         bm_to_litter(ipts,ivm,icarbres,:) = bm_to_litter(ipts,ivm,icarbres,:) + &
                                    vartmp_bm(icarbres,:)*(1-alloc_firefuel(ihour1))
                         bm_to_litter(ipts,ivm,ilabile,:) = bm_to_litter(ipts,ivm,ilabile,:) + &
                                    vartmp_bm(ilabile,:)*(1-alloc_firefuel(ihour1))

                         bm_to_litter(ipts,ivm,isapbelow,:) = bm_to_litter(ipts,ivm,isapbelow,:) + &
                                    vartmp_bm(isapbelow,:)
                         bm_to_litter(ipts,ivm,iheartbelow,:) = bm_to_litter(ipts,ivm,iheartbelow,:) + &
                                    vartmp_bm(iheartbelow,:)
                         bm_to_litter(ipts,ivm,iroot,:) = bm_to_litter(ipts,ivm,iroot,:) + &
                                    vartmp_bm(iroot,:)
                      ENDIF

                      ! Save woody litter by icut. This will be used in pest
                      ! module
                      ! (stomate_pest) for one of the susceptibility indice.
                      IF (icut .NE. icut_fire) THEN
                        woody_litter_by_cut(ipts,ivm,icut) = &
                            woody_litter_by_cut(ipts,ivm,icut) + &
                            (circ_class_biomass(ipts,ivm,icir,isapabove,icarbon) + &
                            circ_class_biomass(ipts,ivm,icir,iheartabove,icarbon)) * &
                            circ_class_kill(ipts,ivm,icir,ifm,icut)
                      ELSE
                        woody_litter_by_cut(ipts,ivm,icut) = &
                            woody_litter_by_cut(ipts,ivm,icut) + &
                            (circ_class_biomass(ipts,ivm,icir,isapabove,icarbon) + &
                            circ_class_biomass(ipts,ivm,icir,iheartabove,icarbon)) * &
                            circ_class_kill(ipts,ivm,icir,ifm,icut) * (1-alloc_firefuel(ihour1))

                      ENDIF

                      ! remove the number of individuals that died
                      circ_class_n(ipts,ivm,icir) = circ_class_n(ipts,ivm,icir) - &
                           circ_class_kill(ipts,ivm,icir,ifm,icut)

                      ! Here, it's possible that the number of individuals left
                      ! will be a very, very small amount.  If it's very low,
                      ! we will just move all that biomass to the dead pool
                      ! and set the number of individuals equal to zero,
                      ! effectively killing the circ class.
                      IF(circ_class_n(ipts,ivm,icir) .LE. min_stomate)THEN

                         bm_to_litter(ipts,ivm,:,:) = bm_to_litter(ipts,ivm,:,:) + &
                              circ_class_biomass(ipts,ivm,icir,:,:) * &
                              circ_class_n(ipts,ivm,icir)

                         ! Save woody litter by icut. This will be used in the pest module
                         ! (stomate_pest) for one of the susceptibility indice.
                         ! Since this process clears the stand it will be saved to
                         ! icut_clear
                         woody_litter_by_cut(ipts,ivm,icut_clear) = &
                              woody_litter_by_cut(ipts,ivm,icut_clear) + &
                              ( circ_class_biomass(ipts,ivm,icir,isapabove,icarbon) + &
                              circ_class_biomass(ipts,ivm,icir,iheartabove,icarbon) ) * &
                              circ_class_n(ipts,ivm,icir)

                         ! zero the number of individuals in this circ class
                         circ_class_n(ipts,ivm,icir) = circ_class_n(ipts,ivm,icir) - &
                              circ_class_n(ipts,ivm,icir)

                         ! If there are no individuals left then the biomass
                         ! in that circ should be set to zero. If not some of 
                         ! the IF-loops will fail because circ_class_n = 0 and
                         ! circ_class_biomass = 0 is used as a logic test later
                         ! in the code
                         circ_class_biomass(ipts,ivm,icir,:,:) = zero

                      ENDIF

                      ! Debug
                      IF(test_pft == ivm .AND. test_grid == ipts .AND. &
                           printlev_loc>=4) THEN
                         WRITE(numout,*) 'killing after: ',ipts,ivm,icir,ifm,icut,&
                              circ_class_kill(ipts,ivm,icir,ifm,icut),&
                              circ_class_n(ipts,ivm,icir)
                      ENDIF
                      !-

                   ENDDO ! loop over circ classes

                   IF (SUM(circ_class_n(ipts,ivm,:)).LT.min_stomate) THEN

                      IF (SUM(SUM(SUM(circ_class_biomass(ipts,ivm,:,:,:),1),1),1).LT.&
                           min_stomate) THEN

                         ! The PFT is dead. Set plant_status to idead to ensure
                         ! plant_status will be set to iprescribe in 
                         ! mortality_clean.
                         plant_status(ipts,ivm) = idead

                      ELSE

                         WRITE(numout,*) 'sum of circ_class_n, ', &
                              SUM(circ_class_n(ipts,ivm,:))
                         WRITE(numout,*) 'sum of circ_biomass_n - C, ', &
                              SUM(SUM(circ_class_biomass(ipts,ivm,:,:,icarbon),1),1)
                         WRITE(numout,*) 'sum of circ_biomass_n - N, ', &
                              SUM(SUM(circ_class_biomass(ipts,ivm,:,:,initrogen),1),1)
                         CALL ipslerr_p(3,'inconsistency in stomate_kill',&
                              'Both indiviudals and biomass should be zero', &
                              'or both should be above min_stomate','')
                      END IF
                   END IF

                   ! Actual biomass after mortality has been accounted for
                   actu_biomass = cc_to_biomass(npts,nvm,&
                        circ_class_biomass(ipts,ivm,:,:,:),&
                        circ_class_n(ipts,ivm,:))

                   ! Calculate the biomass change due to ncuts. Add to the values
                   ! that were already accumulated in anthropogenic_mortality in 
                   ! sapiens_kill.f90
                   biomass_cut(ipts,ivm,:,:,ifm,icut) = &
                        biomass_cut(ipts,ivm,:,:,ifm,icut) + &
                        init_biomass(ipts,ivm,:,:) - actu_biomass(:,:)

                   ! Update the initial biomass to avoid double counting
                   init_biomass(ipts,ivm,:,:) = actu_biomass(:,:)


                ELSE ! Grasses

                   IF(natural(ivm)) THEN
                    
                      IF(circ_class_n(ipts,ivm,1) .GT. min_stomate) THEN
                         IF ((ifm .EQ. ifm_none) .AND. (icut .EQ. icut_thin)) THEN
                            ! Conditions to decrease density include:
                            ! (1)The flag ok_dyn_grass_den for the calculation of dynamic grass density is
                            ! switched on; (2) veget_max for each grassland is greater than 10-8;
                            ! (3) The plant_status is at ipresenescence, isenescent or idormant;
                            ! (4) Current reserve and labile carbon are lower than the target of reserve
                            ! and labile carbon;(5) circ_class_kill >0, which
                            ! means density will be decreased; (6) The target of reserve and labile carbon 
                            ! is greater than 10-8;(7)The sum of carbon biomass is greater than 10-8; 
                            ! (8)The reserve and labile carbon is greater than 10-8; 
                            ! (9)Grass density is greater than the minimal value 0.05.
                            IF (( ok_dyn_grass_den) .AND. &
                               (veget_max(ipts,ivm) .GT. min_stomate) .AND. &
                               ((plant_status(ipts,ivm) .EQ. ipresenescence) .OR. &
                               (plant_status(ipts,ivm) .EQ. isenescent) .OR. &
                               (plant_status(ipts,ivm) .EQ. idormant)) .AND. &
                               (circ_class_n(ipts,ivm,1)*(circ_class_biomass(ipts,ivm,1,icarbres,icarbon) + &
                               circ_class_biomass(ipts,ivm,1,ilabile,icarbon)) .LT.&
                               (tot_res_target(ipts,ivm,icarbon) + &
                               tot_lab_target(ipts,ivm,icarbon))) .AND. &
                               (circ_class_kill(ipts, ivm, 1,ifm_none,icut_thin) .GT. 0) .AND. &
                               ((tot_res_target(ipts,ivm,icarbon) + &
                               tot_lab_target(ipts,ivm,icarbon)) .GT.min_stomate) .AND. &
                               (SUM(circ_class_biomass(ipts,ivm,1,:,icarbon)).GT. min_stomate) .AND. &
                               (circ_class_biomass(ipts,ivm,1,icarbres,icarbon) .GE. min_stomate) .AND. &
                               (circ_class_biomass(ipts,ivm,1,ilabile,icarbon) .GE. min_stomate) .AND. &
                               (circ_class_n(ipts,ivm,1) .GT. 0.05)) THEN
                               
                               ! Debug
                               IF (printlev_loc>=4) THEN
                                 WRITE(numout,*) 'Before density reduction,ipts,ivm',ipts, ivm
                                 WRITE(numout,*) 'Density, before',circ_class_n(ipts,ivm,1)
                                 WRITE(numout,*) 'current reserve',circ_class_biomass(ipts,ivm,1,icarbres,icarbon)
                                 WRITE(numout,*) 'total reserve target',tot_res_target(ipts,ivm,icarbon)
                                 WRITE(numout,*) 'current labile',circ_class_biomass(ipts,ivm,1,ilabile,icarbon)
                                 WRITE(numout,*) 'total labile target',tot_lab_target(ipts,ivm,icarbon)
                                 WRITE(numout,*) 'biomass', &
                                       circ_class_biomass(ipts,ivm,1,:,icarbon)
                                 WRITE(numout,*) 'bm_to_litter', bm_to_litter(ipts, ivm, :, icarbon)

                               END IF

                               ! Update density and biomass 
                               ! To decrease grassland density,the calculation should follow the two principles:
                               ! 1. Mass balance should be preserved before and after the density change, so:
                               ! eq.1: N0*(Cother+Trl0) = N1*(Cother + Trl1) +scale_carbon*(N0-N1)*Cother
                               ! 2. Cother will remain constant, and individual reserve and labile
                               ! carbon will reach the target:((N0+scale_dec*N1)/N1)*(Ttarget_area)/N0
                               ! Combining 1 & 2, we can calculate new density
                               ! Keep old density and biomass before decreasing
                               ! density
                               old_circ_class_n(ipts,ivm,1) = circ_class_n(ipts,ivm,1)
                               old_circ_class_biomass(:,:,:,:,:) =circ_class_biomass(:,:,:,:,:)

                               ! Calculate the sum of reserve and labile carbon
                               ! Trl0 (gC ind-1)
                               Trl0 =old_circ_class_biomass(ipts,ivm,1,ilabile,icarbon) +&
                                     old_circ_class_biomass(ipts,ivm,1,icarbres,icarbon)

                               ! Calculate the sum of carbon from all compartments except for reserve and
                               ! labile pool, Cother (gC ind-1)
                               Cother = SUM(old_circ_class_biomass(ipts,ivm,1,:,icarbon))-&
                                     Trl0

                               ! Calculate the sum of areal target for reserve
                               ! and labile carbon, Ttarget_area (gC m-2) 
                               Ttarget_area = tot_res_target(ipts,ivm,icarbon) + &
                                  tot_lab_target(ipts,ivm,icarbon)

                               ! New density can be calculated using
                               ! circ_class_kill passed from
                               ! stomate_mark_kill.f90
                               circ_class_n(ipts,ivm,1) = circ_class_n(ipts,ivm,1) - &
                                   circ_class_kill(ipts,ivm,1,ifm_none,icut_thin)

                               ! Calculate ratio of old density to new density
                               ratio_den = old_circ_class_n(ipts,ivm,1)/circ_class_n(ipts,ivm,1)

                               ! Update biomass for reserve and labile carbon,
                               ! when density > 0.05, use
                               ! (scale_dec+N0/N1)*Ttarget_area/N0
                               IF (circ_class_n(ipts,ivm,1) .GT. 0.050001) THEN
                                  circ_class_biomass(ipts,ivm,1,icarbres,icarbon) = & 
                                     (scale_dec+old_circ_class_n(ipts,ivm,1)/circ_class_n(ipts,ivm,1))*&
                                     tot_res_target(ipts,ivm,icarbon)/ old_circ_class_n(ipts,ivm,1)
                                  circ_class_biomass(ipts,ivm,1,ilabile,icarbon) = &
                                     (scale_dec+old_circ_class_n(ipts,ivm,1)/circ_class_n(ipts,ivm,1))*&
                                     tot_lab_target(ipts,ivm,icarbon)/ old_circ_class_n(ipts,ivm,1)

                               ELSEIF (circ_class_n(ipts,ivm,1) .LE. 0.050001) THEN
                                  
                                  ! When density is 0.05, which means it has reached the minimal limit
                                  ! and constrained to 0.05, so  we will keep mass balance in this case.
                                  ! To make eq.1 valid, Trl1 (sum of new reserve
                                  ! and labile carbon,gC ind-1) is calculated:
                                  ! Trl1 = (N0*((1-scale_carbon)*Cother+Trl0))/N1-(1-scale_carbon)*Cother

                                  ! Individual reserve and labile carbon is
                                  ! calculated based on ratio and Trl1
                                  ratio_res = old_circ_class_biomass(ipts,ivm,1,icarbres,icarbon) / Trl0
                                  ratio_lab = un - ratio_res                       
                                  circ_class_biomass(ipts,ivm,1,icarbres,icarbon)= ratio_res*&
                                      ((old_circ_class_n(ipts,ivm,1)*((1-scale_carbon)*Cother+Trl0))/ &
                                      circ_class_n(ipts,ivm,1)-(1-scale_carbon)*Cother)
                                  circ_class_biomass(ipts,ivm,1,ilabile,icarbon)= ratio_lab*&
                                      ((old_circ_class_n(ipts,ivm,1)*((1-scale_carbon)*Cother+Trl0))/ &
                                      circ_class_n(ipts,ivm,1)-(1-scale_carbon)*Cother)

                               END IF ! new density > or < 0.05
                               
                               ! Calculate the circ_class_biomass for nitrogen for
                               ! mass balance, by multiplying by the ratio of
                               ! old density to new density   
                               DO ipar = 1,nparts
                                 circ_class_biomass(ipts,ivm,1,ipar,initrogen) = &
                                 circ_class_biomass(ipts,ivm,1,ipar,initrogen)*ratio_den
                               ENDDO ! ipar = 1,nparts

                               ! move the dead biomass to the respective litter
                               ! pool, except for reserve and labile pool. The
                               ! proportion is controled by scale_carbon
                               DO ipar = 1, nparts
                                  IF (ipar .NE. icarbres .AND. ipar .NE. ilabile) THEN
                                     bm_to_litter(ipts, ivm, ipar, icarbon) = bm_to_litter(ipts, ivm, ipar, icarbon) + &
                                        scale_carbon*circ_class_biomass(ipts, ivm, 1, ipar, icarbon) * &
                                        circ_class_kill(ipts, ivm, 1, ifm_none, icut_thin)
                                  END IF ! not reserve and labile carbon

                               ENDDO
 
                               ! Debug
                               IF (printlev_loc >= 4) THEN
                                  WRITE(numout,*) 'After density reduction,ipts,ivm', ipts,ivm
                                  WRITE(numout,*) 'Density, after',circ_class_n(ipts,ivm,1)
                                  WRITE(numout,*) 'current reserve',circ_class_biomass(ipts,ivm,1,icarbres,icarbon)
                                  WRITE(numout,*) 'total reserve target',tot_res_target(ipts,ivm,icarbon)
                                  WRITE(numout,*) 'current labile',circ_class_biomass(ipts,ivm,1,ilabile,icarbon)
                                  WRITE(numout,*) 'total labile target',tot_lab_target(ipts,ivm,icarbon)
                                  WRITE(numout,*) 'biomass', &
                                    circ_class_biomass(ipts,ivm,1,:,icarbon)
                                  WRITE(numout,*) 'killed density',circ_class_kill(ipts,ivm,1,ifm_none,icut_thin)
                                  WRITE(numout,*) 'bm_to_litter', bm_to_litter(ipts, ivm, :, icarbon)
                               END IF !Debug
   
                            ELSE
                               
                               ! This is not very clean.  We don't use circ classes
                               ! for grasses and crops.  All mortality is done based on biomass.
                               ! However, it's cleaner here to just pass circ_class_kill from
                               ! stomate_mark_kill, and it is in line with what is done for the forests.
                               ! So we take the total grass/crop biomass that is scheduled for
                               ! killing in stomate_mark_kill and define
                               ! circ_class_kill(:,:,inatural,1) and
                               ! circ_class_biomass(:,:,1,:,:) so that the product of
                               ! these two matches the amount of biomass scheduled for killing.
                               ! move the dead biomass to the respective litter pool
                               bm_to_litter(ipts,ivm,:,:) = bm_to_litter(ipts,ivm,:,:) +&
                                    circ_class_biomass(ipts,ivm,1,:,:)*&
                                    circ_class_kill(ipts,ivm,1,ifm_none,icut_thin)

                               ! circ_class_biomass has to change such that the biomass
                               ! and circ_class_biomass*circ_class_n are in sync.  So
                               ! we  will just sync it here.  This has to be done this way
                               ! since circ_class_n never changes for grass/crops.
                               ! Below we combine the calculation of the new biomass
                               ! biomass = circ_class_biomass * circ_class_n
                               ! biomass_killed = circ_class_biomass * circ_class_kill
                               ! Given that we keep circ_class_n constant for grass
                               ! and crops we recalculate circ_class_biomass as
                               ! (biomass-biomass_killed)/circ_class_n
                               circ_class_biomass(ipts,ivm,1,:,:) = &
                                    circ_class_biomass(ipts,ivm,1,:,:) * &
                                    (circ_class_n(ipts,ivm,1) - &
                                    circ_class_kill(ipts,ivm,1,ifm_none,icut_thin)) / &
                                    circ_class_n(ipts,ivm,1)                        
                            END IF ! ok_dyn_grass_den
                         END IF !ifm, icut

                         ! Debug
                         IF (printlev_loc>=4) THEN
                            WRITE(numout,*) 'Before judge to idead,ipts,ivm',ipts, ivm
                            WRITE(numout,*) 'Biomass',&
                                 circ_class_biomass(ipts,ivm,1,:,icarbon)
                            WRITE(numout,*) 'Density',&
                                 circ_class_n(ipts,ivm,1)
                         END IF
                         !-

                         ! For grasslands, when reserve and labile carbon is
                         ! zero, it should indicate idead. 
                         IF (circ_class_biomass(ipts,ivm,1,ilabile,icarbon) + &
                              circ_class_biomass(ipts,ivm,1,icarbres,icarbon) .LE. min_stomate) THEN
                            
                            ! The plant is dead. Set plant_status to idead to ensure
                            ! plant_status will be set to iprescribe in 
                            ! mortality_clean.
                            plant_status(ipts,ivm) = idead
                            bm_to_litter(ipts,ivm,:,:) = bm_to_litter(ipts,ivm,:,:) + &
                                 circ_class_biomass(ipts,ivm,1,:,:)*&
                                 circ_class_n(ipts,ivm,1)

                            !Debug
                            IF (printlev_loc>=4) THEN
                               WRITE(numout,*) 'Yes!idead,ipts,ivm',ipts, ivm
                            END IF
                            !-

                            ! zero the number of individuals in this circ class
                            circ_class_n(ipts,ivm,1) = zero

                            ! If there are no individuals left then the biomass
                            ! in that circ should be set to zero. If not some of 
                            ! the IF-loops will fail because circ_class_n = 0 and
                            ! circ_class_biomass = 0 is used as a logic test
                            ! laterin the code
                            circ_class_biomass(ipts,ivm,1,:,:) = zero

                         ENDIF

                      ENDIF

                   ENDIF ! is.natural

                ENDIF ! checking for a tree

                ! This PFT should never be killed again.  
                ! To make sure of that, reset all the
                ! variables to zero.
                circ_class_kill(ipts,ivm,:,ifm,icut) = zero

             ENDDO ! reason of mortality

          ENDDO ! management types
       
       ENDDO  ! loop over pfts

    END DO ! loop over land points



    !! 4. Check numerical consistency of this routine

    IF (err_act.GT.1) THEN


       ! 4.2 Check surface area
       CALL check_vegetation_area("stomate_natural_mortality", npts, &
            veget_max_begin, veget_max,'pft')

       ! 4.3 Mass balance closure 
       pool_end = zero
       DO ipar = 1,nparts
          DO iele = 1,nelements
             DO icir = 1,ncirc
                pool_end(:,:,iele) = pool_end(:,:,iele) + &
                     (circ_class_biomass(:,:,icir,ipar,iele) * &
                     circ_class_n(:,:,icir) * veget_max(:,:))
             ENDDO
             pool_end(:,:,iele) = pool_end(:,:,iele) + &
                  bm_to_litter(:,:,ipar,iele) * veget_max(:,:)
          ENDDO
       ENDDO

       ! 4.3.2 Calculate mass balance
       ! Common processes
       DO iele=1,nelements
          check_intern(:,:,ipoolchange,iele) = -un * (pool_end(:,:,iele) - &
               pool_start(:,:,iele))
          check_intern(:,:,iland2atm,iele) = -un * emissions_fire(:,:,iele,ifirecrown) * veget_max(:,:) 
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
       CALL check_mass_balance("stomate_natural_mortality", closure_intern, &
            npts, pool_end, pool_start, veget_max, 'pft')

    ENDIF ! err_act.GT.1

    IF (printlev.GE.4) WRITE(numout,*) 'Leaving natural_mortality'

  END SUBROUTINE natural_mortality


!! ================================================================================================================================
!! SUBROUTINE   : mortality_clean
!!
!>\BRIEF        After biomass has been killed, there are some clean-up operations
!!              to do.
!!
!! DESCRIPTION  : If all the biomass has been killed for a PFT, we want to reset
!!                some counters.  There is also the chance that we have
!!                killed all the biomass in a circ class of trees.  If this is
!!                the case, we want to redistribute the biomass in the remaining
!!                classes so that all circ classes have some biomass in them.
!!                A circ class with no biomass causes allocation to crash.
!!
!! RECENT CHANGE(S): 
!!
!! MAIN OUTPUT VARIABLE(S): ::circ_class_biomass
!!
!! REFERENCE(S)   :
!!
!! FLOWCHART    : None
!!\n
!_ ================================================================================================================================
 
  SUBROUTINE mortality_clean (npts, circ_class_biomass, circ_class_n, dt_days, &
       veget_max, species_change_map, fm_change_map, &
       longevity_eff_leaf, longevity_eff_sap, longevity_eff_root, &
       lpft_replant, bm_sapl_2D, wstress_season, wstress_month, &
       nstress_season, n_input, n_input_daily, vegstress_season, &
       vegstress_month, vegstress_week, vegstress, humrel, season_drought_legacy, &
       everywhere, PFTpresent, sugar_load, cn_leaf_min_season, &
       cn_leaf_init_2D, Light_Abs_Tot, Light_Tran_Tot, &
       laieff_isotrop, lm_lastyearmax, &
       lm_thisyearmax, age, &
       leaf_frac, leaf_age, light_tran_to_floor_season, &
       qsintveg, us, when_growthinit, gdd_from_growthinit, &
       gdd_midwinter, time_hum_min, hum_min_dormance, gdd_m5_dormance, &
       ncd_dormance, ngd_minus5, mean_start_gs, &
       maxvegstress_lastyear, maxvegstress_thisyear, &
       minvegstress_lastyear, minvegstress_thisyear, &
       maxgppweek_lastyear, maxgppweek_thisyear, &
       n_reserve_longterm, n_reserve_balance, &
       maxfpc_lastyear, maxfpc_thisyear, &
       turnover_longterm, dead_leaves, grow_season_len, &
       KF, atm_to_bm, npp_longterm, croot_longterm, &
       gpp_daily, gpp_year, gpp_decade, resp_maint, resp_growth, resp_hetero, &
       npp_daily, rue_longterm, leaching_daily, &
       emission_daily, gpp_week, resp_maint_week, mai, &
       pai, mai_count, previous_wood_volume, age_stand, &
       last_cut, k_latosa_adapt, litter, bm_to_litter, &
       tree_bm_to_litter, leaf_age_crit, leaf_classes, co2_fire, litterfuel, &
       turnover_daily, soil_n_min, p_O2, bact, &
       deepSOM_a, deepSOM_s, deepSOM_p, som, &
       lignin_struc, lignin_wood, lignin_snag, forest_managed, plant_status, &
       kill_vessels, vessel_loss_previous, biomass_init_drought, &
       count_daylight, CN_som_litter_longterm, zf_soil, &
       MatrixV, MatrixA, VectorU, VectorB, gap_area_save, total_ba_init)
    
    !! 0. Variable and parameter declaration
    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                         :: npts                  !! Domain size (-)
    REAL(r_std), DIMENSION(:,:), INTENT(in)            :: longevity_eff_root    !! Effective root turnover time that accounts
                                                                                !! waterstress (days)
    REAL(r_std), DIMENSION(:,:), INTENT(in)            :: longevity_eff_sap     !! Effective sapwood turnover time that accounts
                                                                                !! waterstress (days)
    REAL(r_std), DIMENSION(:,:), INTENT(in)            :: longevity_eff_leaf    !! Effective leaf turnover time that accounts
                                                                                !! waterstress (days)  
    REAL(r_std), INTENT(in)                            :: dt_days               !! Time step of vegetation dynamics for stomate (days)
    LOGICAL, DIMENSION(:,:), INTENT(in)                :: lpft_replant          !! Set to true if a PFT has been clearcut
                                                                                !! and needs to be replaced by another species
    INTEGER(i_std), DIMENSION (:,:), INTENT(in)        :: species_change_map    !! A map which gives the PFT number that each
    INTEGER(i_std), DIMENSION(:,:), INTENT(in)         :: fm_change_map         !! A map which gives the desired FM strategy when
                                                                                !! the PFT will be replanted after a clearcut.
                                                                                !! (1-nvm,unitless)
    REAL(r_std),DIMENSION(:,:), INTENT(in)             :: cn_leaf_init_2D       !! initial leaf C/N ratio
    REAL(r_std), DIMENSION(0:ngrnd), INTENT(in)        :: zf_soil    

    !! 0.3 Modified variables
    REAL(r_std),DIMENSION(:,:,:,:,:), INTENT(inout)    :: bm_sapl_2D            !! biomass of sapling
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)   :: circ_class_biomass    !! Biomass of the components of the model  
                                                                                !! tree within a circumference
                                                                                !! class @tex $(gC ind^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: circ_class_n          !! Number of individuals in each circ class
                                                                                !! @tex $(ind m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: veget_max             !! "maximal" coverage fraction of a PFT on the ground
                                                                                !! (unitless, 0-1)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: wstress_season        !! Water stress factor, based on hum_rel_daily
                                                                                !! (unitless, 0-1)  
    REAL(r_std),DIMENSION(:,:),INTENT(inout)           :: wstress_month         !! Water stress factor, based on hum_rel_daily
                                                                                !! (unitless, 0-1)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: nstress_season        !! N-related seasonal stress (used for allocation) 
    REAL(r_std),DIMENSION(:,:,:,:), INTENT(inout)      :: n_input               !! Nitrogen inputs into the soil (gN/m**2/timestep)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: n_input_daily         !! nitrogen inputs into the soil  (gN/m**2/day)
                                                                                !! NH4 and NOX from the atmosphere, NH4 from BNF,
                                                                                !! agricultural fertiliser as NH4/NO3 n_input_daily(ipts,imov,:) = zero
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: vegstress_season      !! Mean growingseason moisture availability (0 to 1, unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: vegstress_month       !! "Monthly" moisture availability (0 to 1, 
                                                                                !! unitless) 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: vegstress_week        !! "Weekly" moisture availability 
                                                                                !! (0 to 1, unitless)
    REAL(r_std),DIMENSION(:,:), INTENT(inout)          :: vegstress             !! Relative soil moisture (0-1, unitless) 
    REAL(r_std),DIMENSION(:,:), INTENT(inout)          :: humrel                !! Relative humidity. Not used in stomate (needed in age_class_distr)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: season_drought_legacy !! Mean growing season moisture availability
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: everywhere            !! is the PFT everywhere in the grid box or 
                                                                                !! very localized (after its introduction) (?)
    LOGICAL, DIMENSION(:,:), INTENT(inout)             :: PFTpresent            !! Tab indicating which PFTs are present in 
                                                                                !! each pixel
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: sugar_load            !! Relative sugar loading of the labile pool (unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: cn_leaf_min_season    !! Seasonal min CN ratio of leaves 
    REAL(r_std),DIMENSION (:,:,:), INTENT (inout)      :: Light_Abs_Tot         !! Absorbed radiation per level for photosynthesis
    REAL(r_std),DIMENSION (:,:,:), INTENT (inout)      :: Light_Tran_Tot        !! Transmitted radiation per level for photosynthesis
    REAL(r_std),DIMENSION(:,:,:), INTENT(inout)        :: laieff_isotrop        !! Effective LAI
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: lm_lastyearmax        !! last year's maximum leaf mass for each PFT 
                                                                                !! @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: lm_thisyearmax        !! this year's maximum leaf mass, for each 
                                                                                !! PFT @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: age                   !! mean age (years)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: leaf_frac             !! fraction of leaves in leaf age class (unitless;0-1)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: leaf_age              !! Leaf age (days)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: light_tran_to_floor_season !! Mean seasonal fraction of light transmitted (unitless, 0-1)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: qsintveg              !! Water on vegetation due to interception @tex $(kg m^{-2})$ @endtex
    REAL(r_std),DIMENSION(:,:,:,:),INTENT(inout)       :: us                    !! Water stress index for transpiration
                                                                                !! (by soil layer and PFT) (0-1, unitless)   
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: when_growthinit       !! How many days ago was the beginning of 
                                                                                !! the growing season (days)

    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: gdd_from_growthinit   !! growing degree days, since growthinit 
                                                                                !! for crops
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: gdd_midwinter         !! Growing degree days (K), since midwinter 
                                                                                !! (for phenology) - this is written to the
                                                                                !!  history files 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: time_hum_min          !! Time elapsed since strongest moisture 
                                                                                !! availability (days) 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: hum_min_dormance      !! minimum moisture during dormance 
                                                                                !! (0-1, unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: gdd_m5_dormance       !! Growing degree days (K), threshold -5 deg 
                                                                                !! C (for phenology)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: ncd_dormance          !! Number of chilling days (days), since 
                                                                                !! leaves were lost (for phenology)
 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: ngd_minus5            !! Number of growing days (days), threshold 
                                                                                !! -5 deg C (for phenology) 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: mean_start_gs         !! mean growing season starting day for deciduous PFTs (doy).
    REAL(r_std),DIMENSION(:,:),INTENT(inout)           :: maxvegstress_lastyear !! last year's maximum moisture availability
    REAL(r_std),DIMENSION(:,:),INTENT(inout)           :: maxvegstress_thisyear !! this year's maximum moisture availability
    REAL(r_std),DIMENSION(:,:),INTENT(inout)           :: minvegstress_lastyear !! last year's minimum moisture availability
    REAL(r_std),DIMENSION(:,:),INTENT(inout)           :: minvegstress_thisyear !! this year's minimum moisture availability
    REAL(r_std),DIMENSION(:,:),INTENT(inout)           :: maxgppweek_lastyear   !! last year's maximum weekly GPP
    REAL(r_std),DIMENSION(:,:),INTENT(inout)           :: maxgppweek_thisyear   !! this year's maximum weekly GPP
    REAL(r_std),DIMENSION(:,:),INTENT(inout)           :: n_reserve_longterm    !! "Long term" (default 3 years) actual to potential N
                                                                                !! reserve pool (0-1, unitless)
    REAL(r_std),DIMENSION(:,:), INTENT(inout)          :: n_reserve_balance     !! Actual to potential N reserve pool (unitless)
    REAL(r_std),DIMENSION(:,:),INTENT(inout)           :: maxfpc_lastyear       !! Last year's maximum foliage projected
                                                                                !! coverage for each natural PFT,
                                                                                !! @tex $(m^2 m^{-2})$ @endtex
    REAL(r_std),DIMENSION(:,:),INTENT(inout)           :: maxfpc_thisyear       !! this year's maximum foliage projected
                                                                                !! coverage for each natural PFT,
                                                                                !! @tex $(m^2 m^{-2})$ @endtex
    REAL(r_std),DIMENSION(:,:,:,:),INTENT(inout)       :: turnover_longterm     !! "Long term" turnover rate
                                                                                !! @tex $(gC m^{-2} year^{-1})$ @endtex
    REAL(r_std),DIMENSION(:,:,:),INTENT(inout)         :: dead_leaves           !! Dead leaves on ground, per PFT, metabolic 
                                                                                !! and structural,  
                                                                                !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(:,:),INTENT(inout)           :: grow_season_len       !! growing season length in days for deciduous PFTs.
    

    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: KF                    !! Scaling factor to convert sapwood mass
                                                                                !! into leaf mass (m)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: atm_to_bm             !! CO2 and N taken from the atmosphere to get C to create  
                                                                                !! the seedlings @tex (gC.m^{-2}dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: npp_longterm          !! "Long term" mean yearly primary productivity
                                                                                !! @tex $(m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: croot_longterm        !! "long term" root carbon mass  
 	                                                                        !! @tex ($gC m^{-2}) @endtex
  
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: gpp_daily             !! Daily gross primary productivity  
                                                                                !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: gpp_year              !! "annual" GPP @tex ($gC m^{-2} day^{-1}$)@endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: gpp_decade            !! "decadal" GPP @tex ($gC m^{-2} day^{-1}$)@endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: resp_maint            !! Maintenance respiration  
                                                                                !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: resp_growth           !! Growth respiration  
                                                                                !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: resp_hetero           !! Heterotrophic respiration
                                                                                !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: npp_daily             !! Net primary productivity 
                                                                                !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: rue_longterm          !! Longterm radiation use efficiency 
                                                                                !! (??units??) 
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: leaching_daily        !! mineral nitrogen leached from the soil 
                                                                                !! (gN/m**2/day)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: emission_daily        !! volatile losses of nitrogen 
                                                                                !! (gN/m**2/day)

    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: gpp_week              !! Weekly gross primary productivit
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: resp_maint_week       !! Weekly maintenance respiration
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: mai                   !! The mean annual increment
                                                                                !! @tex $(m**3 / m**2 / year)$ @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: pai                   !! The period annual increment
                                                                                !! @tex $(m**3 / m**2 / year)$ @endtex          
    INTEGER(i_std), DIMENSION(:,:),INTENT(inout)       :: mai_count             !! The number of times we've
                                                                                !! calculated the volume increment
                                                                                !! for a stand 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: previous_wood_volume  !! The volume of the tree trunks
                                                                                !! in a stand for the previous year.
                                                                                !! @tex $(m**3 / m**2 )$ @endtex
    INTEGER(i_std), DIMENSION(:,:), INTENT(inout)      :: age_stand             !! Age of the forest stand (years)  
    INTEGER(i_std), DIMENSION(:,:), INTENT(inout)      :: last_cut              !! Years since last thinning (years)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: k_latosa_adapt        !! Leaf to sapwood area adapted for long 
                                                                                !! term water stress (m)
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)   :: litter                !! metabolic and structural litter, above and 
                                                                                !! below ground @tex ($gC m^{-2}$) @endtex
    
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: bm_to_litter          !! Biomass transfer to litter 
                                                                                !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: tree_bm_to_litter     !! Transfer of tree biomass to litter
                                                                                !! @tex$(gC m^{-2} dtslow^{-1})$@endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: leaf_age_crit         !! critical leaf age (days)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: leaf_classes          !! width of each leaf age class (days) 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: co2_fire              !! Carbon emitted into the atmosphere by 
                                                                                !! fire (living and dead biomass) 
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)   :: litterfuel            !! Dead litter fuel above ground. (gC m^{-2}) 
                                                                                !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: turnover_daily        !! Transfer of litter to som
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: soil_n_min            !! mineral nitrogen in the soil (gN/m**2)  
                                                                                !! (first index=npts, second index=nvm, third index=nnspec) 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: p_O2                  !! partial pressure of oxigen in the soil (hPa)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: bact                  !! denitrifier biomass (gC/m**2)
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: deepSOM_a             !! Soil carbon discretized with depth active (g/m**3) 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: deepSOM_s             !! Soil carbon discretized with depth slow (g/m**3) 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: deepSOM_p             !! Soil carbon discretized with depth passive (g/m**3)
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: som                   !! carbon pool: active, slow, or passive 
                                                                                !! @tex ($gC m^{-2}$) @endtex 
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: lignin_struc          !! ratio Lignine/Carbon in structural litter,
                                                                                !! above and below ground
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: lignin_wood           !! ratio Lignine/Carbon in woody litter,
                                                                                !! above and below ground
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: lignin_snag           !! ratio Lignine/Carbon in snag litter,
                                                                                !! above and below ground

    INTEGER(i_std), DIMENSION (:,:), INTENT(inout)     :: forest_managed        !! forest management flag (is the forest 
                                                                                !! being managed?)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: plant_status          !! Growth and phenological status of the plant
                                                                                !! istatus = Phases defined in constantes
    LOGICAL, DIMENSION(:,:,:), INTENT(inout)           :: kill_vessels          !! Flag to kill vessels at the end of the day following embolism.
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: vessel_loss_previous  !! Proportion of conductivity lost due to cavitation, accumulated
                                                                                !!  on the previous day (unitless).
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)   :: biomass_init_drought  !! Biomass of heartwood or sapwood before onset of drought
    REAL(r_std),DIMENSION(:,:), INTENT(inout)          :: count_daylight        !! Number of time steps dt_radia during daylight
    REAL(r_std), DIMENSION(:,:,:),INTENT(inout)        :: CN_som_litter_longterm!! Longterm CN ratio of litter and som pools (gC/gN)
    
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: matrixA               !! Matrix containing the fluxes between the carbon
                                                                                !! pools per sechiba time step @tex $(gC.m^2.day^{-1})$
                                                                                !! @endtex
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: matrixV               !! Matrix containing the accumulated values of matrixA
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: vectorB               !! Vector containing the litter increase per sechiba
                                                                                !! time step @tex $(gCm^{-2})$ @endtex
    
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: vectorU               !! Matrix containing the accumulated values of VectorB   
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: gap_area_save         !! Stand gap created by more than 30% basal area loss per year.
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: total_ba_init         !! Total basal area per pft saved at the start of the year (m^{2}/m^{-2})

 
    !! 0.4 Local variables
    REAL(r_std)                                        :: share_rec             !! Share of the veget_max of the youngest age class
                                                                                !! (unitless, 0-1)
    INTEGER(i_std)                                     :: ipts,ivm,ilage        !! Indices
    INTEGER(i_std)                                     :: icir,jcir,ilev        !! Indices
    INTEGER(i_std)                                     :: ipar,iele,imbc        !! Indices
    INTEGER(i_std)                                     :: ilit,icarb,igrn       !! Indices
    INTEGER(i_std)                                     :: inc, ispec            !! Indices
    LOGICAL                                            :: lredistribute         !! Flag if we need to redistribute individuals
                                                                                !! among the circ classes
    REAL(r_std), DIMENSION(npts,nvm,nmbcomp,nelements) :: check_intern          !! Contains the components of the internal
                                                                                !! mass balance chech for this routine
                                                                                !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)         :: closure_intern        !! Check closure of internal mass balance
                                                                                !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)         :: pool_start            !! Start pool of this routine 
                                                                                !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)         :: pool_end              !! End pool of this routine 
                                                                                !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(ncirc)                      :: old_trees_left
    REAL(r_std), DIMENSION(ncirc)                      :: circ_class_n_new
    REAL(r_std), DIMENSION(ncirc,nparts,nelements)     :: circ_class_biomass_new
    REAL(r_std)                                        :: trees_needed
    LOGICAL                                            :: lyoungest
    REAL(r_std), DIMENSION(npts,nvm,nleafages)         :: new_leaf_frac         !! Temporary variable for leaf_frac (unitless;0-1)
    REAL(r_std), DIMENSION(ncirc)                      :: est_circ_class_n      !! Temporary variable for circ_class_n of the 
                                                                                !! established vegetation @tex $(ind m^{-2})$ @endtex
    REAL(r_std), DIMENSION(ncirc)                      :: tmp_circ_class_n      !! Temporary variable for circ_class_n of the 
                                                                                !! established vegetation @tex $(ind m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,ncirc)             :: new_circ_class_n      !! Temporary variable for circ_class_n of the
                                                                                !! new vegetation
                                                                                !! @tex $(ind m^{-2})$ @endtex
    REAL(r_std), DIMENSION(ncirc,nparts,nelements)     :: est_circ_class_biomass!! Temporary variable for circ_class_biomass of the
                                                                                !! established vegetation @tex $(g C ind^{-1})$ @endtex
    REAL(r_std), DIMENSION(ncirc,nparts,nelements)     :: tmp_circ_class_biomass!! Temporary variable for circ_class_biomass of the
                                                                                !! established vegetation @tex $(g C ind^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,ncirc,nparts,nelements) &
                                                     :: new_circ_class_biomass !! Temporary variable for circ_class_biomass of the
                                                                               !! new vegetation @tex $(g C ind^{-1})$ @endtex
    REAL(r_std), DIMENSION(ncirc)                    :: est_circ_class_dia     !! Variable to store temporary values for 
                                                                               !! circ_class (m)
    REAL(r_std), DIMENSION(ncirc)                    :: tmp_circ_class_dia     !! Variable to store temporary values for 
                                                                               !! circ_class (m)
    REAL(r_std), DIMENSION(npts,nvm,ncirc)           :: new_circ               !! Temporary variable for circ (m)
    REAL(r_std), DIMENSION(npts,nvm,nelements)       :: new_atm_to_bm          !! Temporary variable for atm_to_bm 
                                                                               !! @tex (gC.m^{-2}dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)                 :: new_age                !! Temporary variable for age (years)
    REAL(r_std), DIMENSION(npts,nvm)                 :: new_lm_lastyearmax     !! Temporary variable for lm_last_year 
                                                                               !! @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm)                 :: new_everywhere         !! Temporary variable for everywhere (unitless, 0-1)
    REAL(r_std), DIMENSION(npts,nvm)                 :: dum_when_growthinit    !! Dummy for when_growthinit (days)
    REAL(r_std), DIMENSION(npts,nvm)                 :: dum_npp_longterm       !! Dummy for npp_longterm
                                                                               !! @tex ($gC m^{-2} year^{-1}$) @endtex
    INTEGER(i_std)                                   :: ivma
    LOGICAL, DIMENSION(npts,nvm)                     :: dum_PFTpresent         !! Dummy for PFTpresent (0 or 1)
    REAL(r_std), DIMENSION(npts,nvm)                 :: dum_plant_status       !! Dummy for plant_status (true/false)
    REAL(r_std)                                      :: moi_bank               !! bank to store vegstress_season of the cleared
                                                                               !! PFT's (unitless; 0-1)
    REAL(r_std), DIMENSION(npts,nvm)                 :: loss_gain              !! The same as delta_veget but distributed of all
                                                                               !! age classes and thus taking the age-classes into 
                                                                               !! account (unitless, 0-1)
    REAL(r_std)                                      :: total_losses           !! Sum of losses only in delta_veget (unitless, 0-1)
    REAL(r_std), DIMENSION(nlitt,nlevs,nelements)    :: litter_bank            !! Bank to store the litter that becomes available when
                                                                               !! part of a PFT is removed @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(ncarb)                    :: soil_bank              !! Bank to store soil carbon that becomes available when
                                                                               !! part of a PFT is removed  @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(nlevs)                    :: struct_ltr_bank        !! Bank to store the lignin in structural litter that 
                                                                               !! becomes available when part of a PFT is removed 
                                                                               !! (0-1,unitless)
    REAL(r_std),DIMENSION(nlevs)                     :: woody_ltr_bank         !! Bank to store the lignin in woody litter that 
                                                                               !! becomes available when part of a PFT is removed 
                                                                               !! (0-1,unitless)
    REAL(r_std),DIMENSION(nvm,nparts,nelements)      :: fresh_litter           !! Pool of fresh litter that is being released during 
                                                                               !! LCC @tex ($gC m^{-2}$) @endtex 
    REAL(r_std),DIMENSION(nparts,nelements)          :: fresh_ltr_bank         !! Bank to store all the fresh litter from site 
                                                                               !! clearing during LCC. Weighted value of fresh_litter
                                                                               !! @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(nlitt,nlevs)              :: litter_weight_rec      !! The fraction of litter of the receiving PFT @tex $-$ @endtex
    INTEGER(i_std)                                   :: iyoung                 !! index of the youngest age class of that species
    REAL(r_std), DIMENSION(npts,nvm,norphans,nelements) &
                                               :: orphan_flux_local            !! Storage for fluxes of PFTs that no longer exist. 
                                                                               !! This is only for co2bm at the moment, since that
                                                                               !! is the only one used in this routine. Following the 
                                                                               !! total destruction of a PFT by death
                                                                               !! (veget_max_new = 0), a flux from before death needs  
                                                                               !! storage @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), DIMENSION(nlitt,nlevs)                :: litter_weight_young  !! The fraction of litter on the young
                                                                               !! PFT.
                                                                               !! @tex $-$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)                   :: veget_max_begin      !! Temporary variable to check area conservation 
    REAL(r_std)                                        :: dummy                !! Dummy variable for prescribe
    REAL(r_std), DIMENSION(npts,nvm)                   :: new_ind              !! Number of recruits grown (trees m-2 day-2). Not used in sapiens_lcchange.
                                                                               !! Added to avoid using OPTINAL arguments. 
    REAL(r_std)                                        :: qm_dia               !! Quadratic Diameter (m)
    REAL(r_std)                                        :: ccn_before           !! Stand stem count before the daily
                                                                               !! redistribution (ind m-2)
    REAL(r_std)                                        :: ccn_after            !! Stand stem count after (ind m-2)
    REAL(r_std)                                        :: cbm_before           !! Stand carbon before the daily
                                                                               !! redistribution (gC m-2)
    REAL(r_std)                                        :: cbm_after            !! Stand carbon after (gC m-2)
    REAL(r_std),DIMENSION(ncirc)                       :: exp_circ_class_n     !! temporary variable for expected circ class n
                                                                               !! @tex $(ind m^{-2}) @endtex
    LOGICAL                                            :: write_debug          !! Flag to write debug statements in prescribe
    REAL(r_std), DIMENSION(npts,nvm,nelements)         :: tot_res_target       !! Target for reserve carbon pool (only for grassland)
                                                                               !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)         :: tot_lab_target       !! Target for labile carbon pool (only for grassland) 
                                                                               !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts)                       :: t2m                  !! Temperature at 2 meter (K)

!_ ================================================================================================================================

    IF (printlev.GE.2) WRITE(numout,*) 'Entering mortality_clean.'

 !! 1. Initialize check for mass balance closure
  
    !! 1.2 Initialize check for mass balance closure 
    IF (err_act.GT.1) THEN

       check_intern(:,:,:,:) = zero
       pool_start(:,:,:) = zero
       DO iele = 1,nelements
          DO ipar = 1,nparts
             DO icir = 1,ncirc
                !  Initial biomass
                pool_start(:,:,iele) = pool_start(:,:,iele) + &
                     (circ_class_biomass(:,:,icir,ipar,iele) * &
                     circ_class_n(:,:,icir) * veget_max(:,:)) 
             ENDDO
             ! bm_to_litter
             pool_start(:,:,iele) = pool_start(:,:,iele) + &
                  bm_to_litter(:,:,ipar,iele) * veget_max(:,:)
          ENDDO

          ! The litter and soil carbon pools can be moved around
          ! if age classes are changed.
       
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
                     deepSOM_p(:,igrn,:,iele))  * &
                    (zf_soil(igrn)-zf_soil(igrn-1)) * veget_max(:,:)
             END DO
          ELSE
             ! Soil carbon (gC m-2) *  (m2 m-2)
             DO icarb = 1,ncarb
                pool_start(:,:,iele) = pool_start(:,:,iele) + &
                     som(:,icarb,:,iele) * veget_max(:,:)
             ENDDO
          ENDIF

          check_intern(:,:,iatm2land,iele) = -un * &
               atm_to_bm(:,:,iele) * veget_max(:,:)

       ENDDO

       ! Denitrifying bacteria (only C)
       pool_start(:,:,icarbon) = pool_start(:,:,icarbon) + &
            bact(:,:) * veget_max(:,:)

       ! Nitrogen only
       DO ispec = 1,nnspec
          pool_start(:,:,initrogen) = pool_start(:,:,initrogen) + &
               soil_n_min(:,:,ispec) * veget_max(:,:)
       ENDDO

       ! Specific fluxes
       ! The fluxes should not be accounted for here because this
       ! code deals with mortality. Following mortality the veget_max
       ! has to be moved to the youngest age class but the fluxes should 
       ! stay in the age class for which they were generated. Clearly, 
       ! closing the mass balance with the fluxes, requires moving
       ! the fluxes. If done so, the NBP consistency check will fail
       ! because those fluxes are double counted: once after iage
       ! in the old age class and once after imor in the youngest
       ! age class. NThe mass balance is still checked for the pools

       !! 1.3 Initialize check for area conservation
       veget_max_begin(:,:) = veget_max(:,:)

    ENDIF ! err_act.GT.1


 !! 2. Redistribute biomass
  
    DO  ipts = 1,npts ! loop over land_points

       DO ivm = 2,nvm  ! loop over #PFT

          IF(veget_max(ipts,ivm) == zero)THEN
             ! this vegetation type is not present, so no reason to do the 
             ! calculation.
             CYCLE
          ENDIF

          ! First we check to see if one of the circ classes is empty.  
          ! We only need to do this if there is some biomass, since we 
          ! don't care if all of the circ classes are empty.
          IF(SUM(SUM(circ_class_biomass(ipts,ivm,:,:,icarbon),1)) .GT. min_stomate)THEN
             
             IF(is_tree(ivm))THEN

                ! By setting this to .TRUE. redistribution
                ! will be taken care of every day. Initially this flag was 
                ! only true at the end of the year but that resulted in 
                ! spikes and pulses.
                lredistribute=.TRUE.

                IF(lredistribute) THEN

                   ! Guillaume M. -- Stem-count conservation check. The mass balance blocks
                   ! only run under err_act > 1 and only aggregate MASSES: a redistribution
                   ! that leaks individuals (fill loop remainder, overwritten class) was
                   ! invisible in production. Read before, compared after the write-back.
                   ccn_before = SUM(circ_class_n(ipts,ivm,:))
                   cbm_before = zero
                   DO icir = 1,ncirc
                      cbm_before = cbm_before + circ_class_n(ipts,ivm,icir) * &
                           SUM(circ_class_biomass(ipts,ivm,icir,:,icarbon))
                   ENDDO

                   ! Debug
                   IF(printlev_loc>=4 .AND. ivm==test_pft)THEN
                      WRITE(numout,*) 'Start of redistributing biomass in mortality_clean'
                      WRITE(numout,*) 'ipts,ivm',ipts,ivm
                      WRITE(numout,*) 'circ_class_n, ',circ_class_n(ipts,ivm,:)
                      WRITE(numout,*) 'circ_class_biomass, ', &
                           SUM(circ_class_biomass(ipts,ivm,:,:,icarbon),2)
                   ENDIF
                   !-

                   ! What we want to do is recreate the circumference class distribution,
                   ! since if we have hit this point we have at least one class that
                   ! is empty and this will cause the allocation to crash.  

                   ! The first step is to decide of the distribution of individuals 
                   ! that we want in each class, for example, a uniform or exponential 
                   ! distribution. This was made with an input parameter. Note
                   ! that this is a decisison with very high consequences. The diamter
                   ! range may vary but the distribution is always the same.

                   ! Now we need to populate the new distribution of trees among
                   ! the circumference classes, and move the heartwood and sap masses 
                   ! to the new distributions. We also need to move the other pools.  
                   ! This is tricky because the allometric relations for the new
                   ! tree sizes will NOT give the same amount of biomass for
                   ! the non-woody pools.  Let us first try it just distributing
                   ! everything equally, and then if that doesn't work we can
                   ! redistribute the non-woody pools more cleverly, even if that 
                   ! means we change the total amount of biomass we have in this
                   ! stand. I want to try it this way because we conserve biomass 
                   ! this way, and I hope that the stress on the trees will be
                   ! small enough that it doesn't cause problems.  This
                   ! redistribution should only happen rarely, so the trees
                   ! should have a chance to equilibrate.
                   old_trees_left(:)=circ_class_n(ipts,ivm,:)

                   qm_dia = wood_to_qmdia(&
                   circ_class_biomass(ipts,ivm,:,:,icarbon), &
                   circ_class_n(ipts,ivm,:),ivm,pipe_tune2(ipts,ivm))

                    ! qm_dia is the diameter at the base of the tree
                    ! in meters.  This is the diameter that we will use
                    ! to calculate the new distribution of trees among
                    ! the circumference classes. 
                   exp_circ_class_n(:) = weibull_class_dist(qm_dia*m_to_cm, ivm, &
                      forest_managed(ipts,ivm))
                     ! now we populate the new classes
                   circ_class_n_new(:) = exp_circ_class_n(:)*SUM(circ_class_n(ipts,ivm,:))

                   circ_class_biomass_new(:,:,:) = zero

                   ! The subsequent calculation nicely closes the carbon
                   ! balance but within the precision 10-16. If there are
                   ! enough trees (first case in the IF-loop), we introduce
                   ! a precision error in the number of trees. This implies
                   ! that when the number of trees is very small, the
                   ! precision issue can become a numerical issue. This is 
                   ! especially true for the largest circumference classes 
                   ! because they contain the fewest individuals. 
                   ! We observed the largest circumference classes 
                   ! containing 10-4 to 10-5 gC less than the previous 
                   ! class. The rest of the code shouldn't care too much 
                   ! whether the circumferences classes are in order or not 
                   ! but to reduce the chances this problem occurs
                   ! we inverted the DO-loops. That way the precision error
                   ! is carried to the smaller diameter classes 
                   ! which have more individuals and so the precision error
                   ! stays a precision issue. Note that this routine is
                   ! called every day so the problem is likely to be 
                   ! corrected the next day.
                   DO icir=ncirc,1,-1

                      trees_needed=circ_class_n_new(icir)

                      DO jcir=ncirc,1,-1

                         IF(trees_needed .LE. zero)EXIT

                         ! Debug
                         IF(printlev_loc>=4 .AND. ivm==test_pft)THEN
                            WRITE(numout,*) 'start of loop'
                            WRITE(numout,*) 'circ_class_biomass_new, C',&
                                 icir,SUM(circ_class_biomass_new(icir,:,icarbon))
                             WRITE(numout,*) 'circ_class_biomass_new, N',&
                                 icir,SUM(circ_class_biomass_new(icir,:,initrogen))
                            WRITE(numout,*) 'old_trees_left, ', &
                                 jcir, old_trees_left(jcir)
                            WRITE(numout,*) 'trees needed, ',trees_needed
                         ENDIF
                         !-

                         IF(old_trees_left(jcir) .GE. trees_needed )THEN

                            ! we can get all the trees we need from this class
                            ! and don't have to continue searching
                            circ_class_biomass_new(icir,:,:)=circ_class_biomass_new(icir,:,:)+&
                                 circ_class_biomass(ipts,ivm,jcir,:,:)*trees_needed
                            old_trees_left(jcir)=old_trees_left(jcir)-trees_needed
                            trees_needed = zero

                            ! Debug
                            IF(printlev_loc>=4 .AND. ivm==test_pft)THEN
                               WRITE(numout,*) 'enough trees'
                            ENDIF
                            !-
   
                            EXIT

                         ELSE

                            ! the trees in this class are not sufficient, so we
                            ! need to move all of them to the new class.
                            circ_class_biomass_new(icir,:,:)=circ_class_biomass_new(icir,:,:)+&
                                 circ_class_biomass(ipts,ivm,jcir,:,:)*&
                                 old_trees_left(jcir)
                            trees_needed = trees_needed-old_trees_left(jcir)
                            old_trees_left(jcir)=zero

                            ! Debug
                            IF(printlev_loc>=4 .AND. ivm==test_pft)THEN
                               WRITE(numout,*) 'not enough trees'
                            ENDIF
                            !-

                         ENDIF
                         
                         ! Debug
                         IF(printlev_loc>=4 .AND. ivm==test_pft)THEN
                            WRITE(numout,*) 'start of loop'
                            WRITE(numout,*) 'circ_class_biomass_new, C',&
                                 icir,SUM(circ_class_biomass_new(icir,:,icarbon))
                             WRITE(numout,*) 'circ_class_biomass_new, N',&
                                 icir,SUM(circ_class_biomass_new(icir,:,initrogen))
                            WRITE(numout,*) 'old_trees_left, ', &
                                 jcir, old_trees_left(jcir)
                            WRITE(numout,*) 'trees needed, ',trees_needed
                         ENDIF
                         !-

                      ENDDO ! jcir=1,ncirc

                   ENDDO ! icir=1,ncirc

                   ! right now, circ_class_biomass_new gives the total biomass
                   ! in each class, but it should only be for a model tree. 
                   ! So let's normalize it.
                   !
                   ! right now, circ_class_biomass_new gives the total biomass
                   ! in each class, but it should only be for a model tree. 
                   ! So let's normalize it.
                   DO icir=1,ncirc
                      IF (circ_class_n_new(icir) .GT. zero) THEN
                         circ_class_biomass_new(icir,:,:) = &
                              circ_class_biomass_new(icir,:,:)/&
                              circ_class_n_new(icir)
                      ENDIF
                   ENDDO

                   ! Guillaume M. -- A class holding almost no stem carried a runaway per-tree
                   ! biomass: the divisor above walks to zero, and `[MESURE]` a class with
                   ! 0,02 stem/ha reached 94 cm growing 1,9 cm/yr where the stand grows 0,5.
                   ! Such a class now ADOPTS the size of its neighbour instead of keeping a
                   ! diameter its handful of trees cannot justify.
                   ! /!\ Its stems are LEFT IN PLACE. Emptying the class was tried and breaks
                   ! the invariant this very routine exists to hold -- growth_fun_all divides
                   ! b_inc_tot by circ_class_n CLASS BY CLASS (stomate_growth_fun_all:2431,
                   ! vector expression), guarded only on the SUM, so an exactly empty class
                   ! yields Inf x 0 = NaN over the whole PFT the next day.
                   IF (cc_n_floor .GT. zero) THEN
                      DO icir=ncirc,2,-1
                         IF (circ_class_n_new(icir) .LE. cc_n_floor .AND. &
                              circ_class_n_new(icir-1) .GT. cc_n_floor) THEN
                            circ_class_biomass_new(icir,:,:) = &
                                 circ_class_biomass_new(icir-1,:,:)
                         ENDIF
                      ENDDO
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
                   CALL sort_circ_class_biomass(circ_class_biomass_new,&
                           circ_class_n_new)

                   ! Now update the variables that pass this information around
                   DO icir=1,ncirc
                      circ_class_biomass(ipts,ivm,icir,:,:) = &
                           circ_class_biomass_new(icir,:,:)
                      circ_class_n(ipts,ivm,icir) = circ_class_n_new(icir)
                   ENDDO

                   ! Guillaume M. -- Conservation verdict. Purely diagnostic (reads + warning),
                   ! strictly bit-neutral, same pattern as check_area_invariant. The expected
                   ! error is machine rounding (~1e-15 relative); a real leak is orders of
                   ! magnitude larger, so the threshold sits between the two.
                   ccn_after = SUM(circ_class_n(ipts,ivm,:))
                   cbm_after = zero
                   DO icir = 1,ncirc
                      cbm_after = cbm_after + circ_class_n(ipts,ivm,icir) * &
                           SUM(circ_class_biomass(ipts,ivm,icir,:,icarbon))
                   ENDDO
                   IF (ccn_before .GT. min_stomate .AND. &
                        ABS(ccn_after-ccn_before) .GT. 1.E-10_r_std*ccn_before) THEN
                      WRITE(numout,*) '[STEMCHECK] ipts',ipts,' ivm',ivm, &
                           ' stems before',ccn_before,' after',ccn_after, &
                           ' rel',ABS(ccn_after-ccn_before)/ccn_before
                      CALL ipslerr_p(2,'mortality_clean', &
                           'stem count not conserved by the daily redistribution',&
                           'see [STEMCHECK] in the log for the point and PFT','')
                   ENDIF
                   IF (cbm_before .GT. min_stomate .AND. &
                        ABS(cbm_after-cbm_before) .GT. 1.E-10_r_std*cbm_before) THEN
                      WRITE(numout,*) '[STEMCHECK] ipts',ipts,' ivm',ivm, &
                           ' stand C before',cbm_before,' after',cbm_after, &
                           ' rel',ABS(cbm_after-cbm_before)/cbm_before
                      CALL ipslerr_p(2,'mortality_clean', &
                           'stand carbon not conserved by the daily redistribution',&
                           'see [STEMCHECK] in the log for the point and PFT','')
                   ENDIF

                ENDIF ! redistribute
                
                ! Debug
                IF(printlev_loc>=4 .AND. ivm==test_pft)THEN
                   WRITE(numout,*) 'End of redistributing biomass in mortality_clean'
                   WRITE(numout,*) 'ipts,ivm',ipts,ivm
                   WRITE(numout,*) 'circ_class_n, ',circ_class_n(ipts,ivm,:)
                   WRITE(numout,*) 'circ_class_biomass - C, ', &
                        SUM(circ_class_biomass(ipts,ivm,:,:,icarbon),2)
                    WRITE(numout,*) 'circ_class_biomass - N, ', &
                        SUM(circ_class_biomass(ipts,ivm,:,:,initrogen),2)
                ENDIF
                !-

             ENDIF ! is_tree
          
          ELSEIF ((SUM(SUM(circ_class_biomass(ipts,ivm,:,:,icarbon),1),1) .LT. min_stomate .AND. &
               plant_status(ipts,ivm) .EQ. idead) .OR. &
               ((circ_class_biomass(ipts,ivm,1,ilabile,icarbon) + &
               circ_class_biomass(ipts,ivm,1,icarbres,icarbon) .LE. min_stomate) .AND. &
               (plant_status(ipts,ivm) .EQ. idead) .AND. &
               .NOT. is_tree(ivm) .AND. natural(ivm)) ) THEN 

             ! The PFT just died. The biomass was already removed in 
             ! stomate_kill so there is no biomass to redistribute. The 
             ! cleaning that should be done when the same PFT is replanted 
             ! and there is only 1 age class is minimal and is done in this
             ! IF-statemet. The moving and cleaning that needs to be done 
             ! when there is more than 1 age class is done in one of the
             ! following IF-statements.

             ! All the biomass was already killed for this site. Some flags to
             ! be reset in this case. It is OK to kill the PFT even if 
             ! we are using species change but we can not replant yet
             ! because we want to replant with a different species this is
             ! basically the same as a land cover change and so it is best 
             ! taken care of with the sapiens_lcchange code. To end up
             ! here we need a veget_max for this PFT. Another check is 
             ! through using ::PFTpresent.
             IF(PFTpresent(ipts,ivm))THEN
                
                !! Reinitialize vegetation characteristics in STOMATE
                IF (veget_max(ipts,ivm).GT.min_stomate) THEN
                   ! The same PFT needs to be replanted
                   plant_status(ipts,ivm) = iprescribe

                ELSE
                   ! The text above suggest that PFTpresent can only be true
                   ! if veget_max is more then zero. That seems to make sense.
                   ! If true, we could reinitialize the PFT right after killing
                   ! it. Double check whether this is true.
                   plant_status(ipts,ivm) = inone

                   ! For the moment make the model crash and have a closer look
                   ! at how it is possible to have PFTpresent = TRUE without
                   ! veget_max. If there is no veget_max, plant status should be
                   ! none.
                   CALL ipslerr_p(3,'mortality_clean','wrong assumption',&
                        'look in the code for more comments','')
                END IF
                age(ipts,ivm) = zero

                ! We could keep the previous sugar_load. That could
                ! be considered local adaptation. By resetting it
                ! as done here, we give the replanted vegetation a bit
                ! of a boost because it will take a bit of time before
                ! gpp is reduced (if at all) by sugar_load. Note that
                ! this is done before move_pft_properties so this new
                ! value will used in that subroutine
                sugar_load(ipts,ivm) = un

                ! Specific variables for forest stands
                IF(is_tree(ivm))THEN
                   !+++CHECK+++
                   ! We could keep the previous when_growthinit. It would
                   ! make the phenology a bit more stable. The problem
                   ! of an unstable phenology has not been documented yet.
                   ! This is more a theoretical consideration.
                   when_growthinit(ipts,ivm) = large_value
                   !+++++++++++
                   last_cut(ipts,ivm) = zero
                   mai_count(ipts,ivm) = zero
                   age_stand(ipts,ivm) = zero
                ENDIF

                !! 1.6 Update leaf ages
                DO ilage = 1, nleafages
                   leaf_age(ipts,ivm,ilage) = zero
                   leaf_frac(ipts,ivm,ilage) = zero
                ENDDO

                ! If we have grasses,
                ! we reset the long term GPP.  This value as 500 in the
                ! old code.  We externalized the variable, but we are
                ! leaving the default at 500.  Many PFTs should have
                ! a long term value of less than 500, so someone could
                ! do some work on this in the future.
                IF (.NOT. is_tree(ivm) .AND. .NOT. ok_constant_mortality) THEN
                   npp_longterm(ipts,ivm) = npp_reset_value(ivm)
                ENDIF
                
             ENDIF ! is PFT present?

          ENDIF ! check if biomass is equal to zero

          ! The part is above is doing the cleaning that is common to
          ! configurations with and without age classes. In the 
          ! following code, the added complexity of the age classes is
          ! being dealt with. 
          IF(is_tree(ivm))THEN
             
             ! This is essentially the same code that is found in land cover change,
             ! since the same process happens there.  Notice that here we do not
             ! need to use the _bank variables, since we assume that if a stand dies
             ! due to natural causes it will always be replaced by the same species.
             ! We also assume that if a stand is clearcut, the model will replace it
             ! by the same species.  To do so otherwise would require something more
             ! like a DGVM. That is what is being done in the species change code.
             ! If the species change code is used, the PFT will be replanted at the
             ! end of the year. Not in the middle of the year as being done here.

             ! If we are using age classes, we need to reset a lot of counters
             ! and change veget_max, moving veget_max from this age class to
             ! the lowest age class.  If veget_max is greater than zero and there
             ! is no biomass, prescribe will attempt to grow trees here in the
             ! next timestep.  This is fine if we have no age classes, but it doesn't
             ! make any sense to prescribe trees that are 50 years old.  We should
             ! only prescribe trees which are 0 years old.
             IF( SUM(SUM(circ_class_biomass(ipts,ivm,:,:,icarbon),1),1) .LT. min_stomate .AND. &
                  (agec_group(ivm)==species_change_map(ipts,ivm) .OR. &
                  .NOT.lpft_replant(ipts,ivm))) THEN 

                ! Debug
                IF(printlev_loc>=4)THEN
                   WRITE(numout,*) 'Getting reset in mortality_clean'
                   WRITE(numout,*) 'ivm,ipts',ivm,ipts
                ENDIF
                !-
                
                IF(nagec .GT. 1)THEN

                   ! First, we need to find the youngest age class of this PFT.
                   iyoung=start_index(agec_group(ivm))
                   IF(ivm == iyoung)THEN
                      lyoungest=.TRUE.
                   ELSE
                      lyoungest=.FALSE.
                   ENDIF

                   ! If we are in the youngest age class, things will be handled properly with
                   ! prescribe in the next step. Nothing needs to be done. Only deal with the
                   ! case where we are not in the youngest age class.
                   IF(.NOT.lyoungest)THEN
                      
                      ! All of our counters for this stand are zero.  We need to merge that
                      ! with the youngest age class.  There are a couple possibilities:
                      ! (1) nothing exists right now in the youngest age class, in 
                      ! which case the veget_max of the old age class becomes that of the 
                      ! youngest and we prescribe .  (2) something 
                      ! exist in the youngest age class, in which case we prescribe 
                      ! to the current age class and then merge biomass with the youngest.  
                      ! (3) both the youngest age class and the current age 
                      ! class have died.  In that case we can do the same thing as if there 
                      ! is no vegetation in the young age class, we just have to make sure 
                      ! to add the veget_max and not replace it.
                      ! As ::veget_max is used in the 
                      ! IF-statements we won't update it until the end of this routine. If 
                      ! we would update before that, the same PFT may be subject to 
                      ! conflicting actions.
                      CALL calculate_shares(veget_max(ipts,iyoung), veget_max(ipts,ivm), &
                           litter(ipts,:,iyoung,:,:), litter(ipts,:,ivm,:,:), &
                           share_rec, litter_weight_rec)

                      IF( SUM(SUM(circ_class_biomass(ipts,iyoung,:,:,icarbon),1),1) .LT. min_stomate)THEN

                         ! Example: an old age class died. In a future time step we will replant
                         ! but we will replant in the youngest age class of this species group.
                         ! This youngest age has no biomass but it may already contain pools and 
                         ! state variables from other age class changes. Therefore pools and state
                         ! variables will be merged with the receiving (young) age class. Replanting
                         ! is simple so we will do it in stomate_prescribe

                         ! Debug
                         IF(printlev_loc>=4)THEN
                            WRITE(numout,*) 'Merging our biomass to the youngest age class.'
                            WRITE(numout,*) 'No current biomass in the youngest age class.'
                            WRITE(numout,*) 'ipts,iyoung,ivm: ',ipts,iyoung,ivm
                            WRITE(numout,*) 'share_young: ',share_rec
                            WRITE(numout,*) 'veget_max(young and current): ',&
                                 veget_max(ipts,iyoung),veget_max(ipts,ivm)
                         ENDIF
                         !-
                
                         !+++CHECK+++
                         ! There was a probably a good reason why not 
                         ! the PFT is not replanted here but it is not
                         ! documented and it is not obvious. Replanting
                         ! here seems like a good thing to do and it
                         ! would make this part of the code more 
                         ! consistent with other parts were replanting
                         ! is done right away. As this functionality 
                         ! is working correctly for the moment it will
                         ! not be touched until there is a good reason
                         ! to change it.
                         !+++++++++++

                         ! Note that sugar_load and when_growthinit 
                         ! are already touched above. This is not a
                         ! really a duplication. In the above code the 
                         ! value of the new pft is given but the new
                         ! pft is still thge same as the old pft. Here
                         ! the value of the old pft will be merged with 
                         ! the value of the receiving pft. What is done
                         ! above affects the final values but does not
                         ! not replace the moving done below. 
                         CALL move_pft_properties(ipts, iyoung, ivm, share_rec, &
                              litter_weight_rec, circ_class_biomass, fm_change_map, &
                              wstress_season(ipts,iyoung), wstress_season(ipts,ivm), &
                              wstress_month(ipts,iyoung), wstress_month(ipts,ivm), &
                              season_drought_legacy(ipts,iyoung,:), &
                              season_drought_legacy(ipts,ivm,:), &
                              nstress_season(ipts,iyoung), nstress_season(ipts,ivm), &
                              n_input(ipts,iyoung,:,:), n_input(ipts,ivm,:,:), &
                              n_input_daily(ipts,iyoung,:), n_input_daily(ipts,ivm,:), &
                              vegstress_season(ipts,iyoung), vegstress_season(ipts,ivm), &
                              vegstress_month(ipts,iyoung), vegstress_month(ipts,ivm), &
                              vegstress_week(ipts,iyoung), vegstress_week(ipts,ivm), &
                              vegstress(ipts,iyoung), vegstress(ipts,ivm), &
                              humrel(ipts,iyoung), humrel(ipts,ivm), &
                              everywhere(ipts,iyoung), everywhere(ipts,ivm), &
                              PFTpresent(ipts,iyoung), PFTpresent(ipts,ivm), &
                              sugar_load(ipts,iyoung), sugar_load(ipts,ivm), &
                              cn_leaf_min_season(ipts,iyoung), cn_leaf_min_season(ipts,ivm), &
                              cn_leaf_init_2D, &
                              Light_Abs_Tot(ipts,iyoung,:), Light_Abs_Tot(ipts,ivm,:), &
                              Light_Tran_Tot(ipts,iyoung,:), Light_Tran_Tot(ipts,ivm,:), &
                              laieff_isotrop(ipts,:,iyoung), laieff_isotrop(ipts,:,ivm), &
                              lm_lastyearmax(ipts,iyoung), lm_lastyearmax(ipts,ivm), &
                              lm_thisyearmax(ipts,iyoung), lm_thisyearmax(ipts,ivm), &
                              age(ipts,iyoung), age(ipts,ivm), &
                              leaf_frac(ipts,iyoung,:), leaf_frac(ipts,ivm,:), &
                              leaf_age(ipts,iyoung,:), leaf_age(ipts,ivm,:), &
                              light_tran_to_floor_season(ipts,iyoung), &
                              light_tran_to_floor_season(ipts,ivm), &
                              qsintveg(ipts,iyoung), qsintveg(ipts,ivm), &
                              us(ipts,iyoung,:,:), us(ipts,ivm,:,:), &
                              when_growthinit(ipts,iyoung), when_growthinit(ipts,ivm), & 
                              gdd_from_growthinit(ipts,iyoung), gdd_from_growthinit(ipts,ivm), &
                              gdd_midwinter(ipts,iyoung), gdd_midwinter(ipts,ivm), &
                              time_hum_min(ipts,iyoung), time_hum_min(ipts,ivm), &
                              hum_min_dormance(ipts,iyoung), hum_min_dormance(ipts,ivm), & 
                              gdd_m5_dormance(ipts,iyoung), gdd_m5_dormance(ipts,ivm), &
                              ncd_dormance(ipts,iyoung), ncd_dormance(ipts,ivm), &
                              ngd_minus5(ipts,iyoung), ngd_minus5(ipts,ivm), &
                              mean_start_gs(ipts,iyoung), mean_start_gs(ipts,ivm), &
                              maxvegstress_lastyear(ipts,iyoung), maxvegstress_lastyear(ipts,ivm), &
                              maxvegstress_thisyear(ipts,iyoung), maxvegstress_thisyear(ipts,ivm), &
                              minvegstress_lastyear(ipts,iyoung), minvegstress_lastyear(ipts,ivm), &
                              minvegstress_thisyear(ipts,iyoung), minvegstress_thisyear(ipts,ivm), &
                              maxgppweek_lastyear(ipts,iyoung), maxgppweek_lastyear(ipts,ivm), &
                              maxgppweek_thisyear(ipts,iyoung), maxgppweek_thisyear(ipts,ivm), &
                              n_reserve_longterm(ipts,iyoung), n_reserve_longterm(ipts,ivm), &
                              n_reserve_balance(ipts,iyoung), n_reserve_balance(ipts,ivm), &
                              maxfpc_lastyear(ipts,iyoung), maxfpc_lastyear(ipts,ivm), &
                              maxfpc_thisyear(ipts,iyoung), maxfpc_thisyear(ipts,ivm), &
                              turnover_longterm(ipts,iyoung,:,:), turnover_longterm(ipts,ivm,:,:), &
                              dead_leaves(ipts,iyoung,:), dead_leaves(ipts,ivm,:), &
                              grow_season_len(ipts,iyoung), grow_season_len(ipts,ivm), &
                              KF(ipts,iyoung), KF(ipts,ivm), &
                              atm_to_bm(ipts,iyoung,:), atm_to_bm(ipts,ivm,:), &
                              npp_longterm(ipts,iyoung), npp_longterm(ipts,ivm), &
                              croot_longterm(ipts,iyoung), croot_longterm(ipts,ivm), &
                              gpp_daily(ipts,iyoung), gpp_daily(ipts,ivm), &
                              resp_maint(ipts,iyoung), resp_maint(ipts,ivm), &
                              resp_growth(ipts,iyoung), resp_growth(ipts,ivm), &
                              resp_hetero(ipts,iyoung), resp_hetero(ipts,ivm), &
                              npp_daily(ipts,iyoung), npp_daily(ipts,ivm), &
                              rue_longterm(ipts,iyoung), rue_longterm(ipts,ivm), &
                              leaching_daily(ipts,iyoung,:), leaching_daily(ipts,ivm,:), &
                              emission_daily(ipts,iyoung,:), emission_daily(ipts,ivm,:), &
                              gpp_week(ipts,iyoung), gpp_week(ipts,ivm), &
                              resp_maint_week(ipts,iyoung), resp_maint_week(ipts,ivm), &
                              gpp_year(ipts,iyoung), gpp_year(ipts,ivm), &
                              gpp_decade(ipts,iyoung), gpp_decade(ipts,ivm), &
                              mai(ipts,iyoung), mai(ipts,ivm), &
                              pai(ipts,iyoung), pai(ipts,ivm), &
                              mai_count(ipts,iyoung), mai_count(ipts,ivm), &
                              previous_wood_volume(ipts,iyoung), previous_wood_volume(ipts,ivm), &
                              age_stand(ipts,iyoung), age_stand(ipts,ivm), &
                              last_cut(ipts,iyoung), last_cut(ipts,ivm), &
                              k_latosa_adapt(ipts,iyoung), k_latosa_adapt(ipts,ivm), &
                              litter(ipts,:,iyoung,:,:), litter(ipts,:,ivm,:,:), &
                              bm_to_litter(ipts,iyoung,:,:), bm_to_litter(ipts,ivm,:,:), &
                              tree_bm_to_litter(ipts,iyoung,:,:), tree_bm_to_litter(ipts,ivm,:,:), &
                              leaf_age_crit(ipts,iyoung), leaf_age_crit(ipts,ivm), &
                              leaf_classes(ipts,iyoung), leaf_classes(ipts,ivm), &
                              co2_fire(ipts,iyoung), co2_fire(ipts,ivm), &
                              litterfuel(ipts,:,iyoung,:,:), litterfuel(ipts,:,ivm,:,:), &
                              turnover_daily(ipts,iyoung,:,:), turnover_daily(ipts,ivm,:,:), &
                              soil_n_min(ipts,iyoung,:), soil_n_min(ipts,ivm,:), &
                              p_O2(ipts,iyoung), p_O2(ipts,ivm), &
                              bact(ipts,iyoung), bact(ipts,ivm), &
                              deepSOM_a(ipts,:,iyoung,:), deepSOM_a(ipts,:,ivm,:), &
                              deepSOM_s(ipts,:,iyoung,:), deepSOM_s(ipts,:,ivm,:), &
                              deepSOM_p(ipts,:,iyoung,:), deepSOM_p(ipts,:,ivm,:), &
                              som(ipts,:,iyoung,:), som(ipts,:,ivm,:), &
                              lignin_struc(ipts,iyoung,:), lignin_struc(ipts,ivm,:), &
                              lignin_wood(ipts,iyoung,:), lignin_wood(ipts,ivm,:), &
                              lignin_snag(ipts,iyoung,:), lignin_snag(ipts,ivm,:), &
                              forest_managed(ipts,iyoung), forest_managed(ipts,ivm), &
                              plant_status(ipts,iyoung), plant_status(ipts,ivm), &
                              kill_vessels(ipts,iyoung,:), kill_vessels(ipts,ivm,:), &
                              vessel_loss_previous(ipts,iyoung,:), vessel_loss_previous(ipts,ivm,:), &
                              biomass_init_drought(ipts,iyoung,:,:,:), biomass_init_drought(ipts,ivm,:,:,:), &
                              count_daylight(ipts,iyoung), count_daylight(ipts,ivm), &
                              CN_som_litter_longterm(ipts,iyoung,:), CN_som_litter_longterm(ipts,ivm,:), &
                              matrixV(ipts,iyoung,:,:), matrixV(ipts,ivm,:,:), &
                              matrixA(ipts,iyoung,:,:), matrixA(ipts,ivm,:,:), &
                              vectorU(ipts,iyoung,:), vectorU(ipts,ivm,:), &
                              vectorB(ipts,iyoung,:), vectorB(ipts,ivm,:), &
                              gap_area_save(ipts,iyoung,:), gap_area_save(ipts,ivm,:), &
                              total_ba_init(ipts,iyoung), total_ba_init(ipts,ivm))
                         !+++++++++++

                         ! Case specific variables
                         plant_status(ipts,iyoung) = iprescribe
                         veget_max(ipts,iyoung) = veget_max(ipts,iyoung)+veget_max(ipts,ivm)
                         
                         !+++CHECK+++
                         ! could be taken care of in the subroutine that cleans the pft that died
                         plant_status(ipts,ivm) = inone
                         veget_max(ipts,ivm) = zero
                         !+++++++++++
                         
                      ELSE

                         ! Example: the youngest age class already contains some
                         ! biomass. We will increase the veget_max with the veget_max
                         ! of another age class that died.

                         ! There are two important things to do here.  First,
                         ! we need to prescribe biomass to this PFT, since part of it
                         ! is currently empty.  Then we need to merge this biomass
                         ! with what is already in the youngest age class of this PFT.

                         ! We want to make use of prescribe but it should 
                         ! be noted that the youngest PFT already contains biomass. 
                         ! Therefore we will create a temporary PFT without biomass which 
                         ! can be used to establish the  new vegetation within the 
                         ! youngest. For most of the INTENT(inout)  variables of 
                         ! prescribe we receive temporary variables back this 
                         ! helps us to calculate the weighted mean of the newly established 
                         ! vegetation and the vegetation that is already there. For some 
                         ! other variables we receive dummies as we are simply not
                         ! going to use these values but will continue using the values of the
                         ! vegetation that is already there.

                         ! Debug
                         IF(printlev_loc>=4)THEN
                            WRITE(numout,*) 'Merging our biomass to the youngest age class.'
                            WRITE(numout,*) 'Biomass already in the youngest age class.'
                            WRITE(numout,*) 'ipts,iyoung,ivm: ',ipts,iyoung,ivm
                            WRITE(numout,*) 'share_young: ',share_rec
                            WRITE(numout,*) 'veget_max(young and current): ',&
                                 veget_max(ipts,iyoung),veget_max(ipts,ivm)
                         ENDIF
                         !-

                         !! 5.2.3.1 Initialize new and dummy variables
                         new_circ_class_biomass(:,:,:,:,:) = zero
                         new_circ_class_n(:,:,:) = zero
                         new_lm_lastyearmax(:,:) = zero
                         new_age(:,:) = zero
                         new_leaf_frac(:,:,:) = zero
                         new_atm_to_bm(:,:,:) = zero
                         new_everywhere(:,:) = zero 
                         dum_when_growthinit(:,:) = large_value
                         dum_npp_longterm(:,:) = zero
                         dum_PFTpresent(:,:) = .TRUE.
                         dum_plant_status(:,:) = iprescribe
                
                         !! 5.2.3.3 Initialize the newly established vegetation: 
                         !  Initialize the newly established vegetation: density of 
                         !  individuals, biomass, allocation factors, leaf age distribution, 
                         !  etc to some reasonable value ::veget_max is not updated yet, 
                         !  therefore loss_gain is passed in the argument list.
                         !  For loss_gain, we only want to replant the current PFT.
                         !  If lpft_replat is TRUE we should never be here in the 
                         !  first place but to be safe lpft_replant is passed into
                         !  prescribe. This is an optional variable that will prevent
                         !  prescribing biomass if species change is used.

                         ! Debug
                         IF(printlev_loc>=4)THEN
                            IF(lpft_replant(ipts,ivm))THEN
                               WRITE(numout,*) 'ERROR - mortality_clean should never be here'
                               IF (err_act.GT.1) CALL ipslerr_p(3, &
                                    'mortality_clean in stomate_kill.f90',&
                                    'species change was activated, replanting needed',&
                                    'do not replant in mortality_clean, wait for lcchange','')
                            ENDIF
                         ENDIF
                         !-

                         ! In mortality_clean stomate_prescribe is only used to establish new
                         ! vegetation. The variables required for recruitment are not 
                         ! passed all the way into this routine but are made dummy
                         ! variables. A more elegant solution is to use OPTIONAL variables
                         ! but the code became a bit confused because lpft_replant is
                         ! already an OPTIONAL variable but for a completly different 
                         ! functionality.
                         loss_gain(:,:)=zero
                         loss_gain(ipts,ivm)=veget_max(ipts,ivm)
                         ! Dummy for light_tran_tot_season
                         dummy = un
                         write_debug = .FALSE.
                         IF (ipts==test_grid.AND.ivm==test_pft) write_debug=.TRUE.
                         IF (printlev_loc.GE.3) THEN
                            WRITE(numout,*) 'Calling precribe from stomate_kill, ', ipts, ivm
                         END IF
                         CALL prescribe (&
                              ivm,                                    loss_gain(ipts,ivm),           dt_days, &
                              dum_PFTpresent(ipts,ivm),               new_everywhere(ipts,ivm),      dum_when_growthinit(ipts,ivm), &
                              new_leaf_frac(ipts,ivm,:),              new_circ_class_n(ipts,ivm,:), &
                              new_circ_class_biomass(ipts,ivm,:,:,:), new_atm_to_bm(ipts,ivm,:),     forest_managed(ipts,ivm), &
                              KF(ipts,ivm),                           plant_status(ipts,ivm),        new_age(ipts,ivm), &
                              dum_npp_longterm(ipts,ivm),             new_lm_lastyearmax(ipts,ivm),  longevity_eff_leaf(ipts,ivm), &
                              longevity_eff_sap(ipts,ivm),            longevity_eff_root(ipts,ivm),  k_latosa_adapt(ipts,ivm), &
                              dummy,                                  species_change_map(ipts,ivm), &
                              cn_leaf_init_2D(ipts,ivm),              bm_sapl_2D(ipts,ivm,:,:,:),    new_ind(ipts,ivm),&
                              pipe_tune2(ipts,ivm),                   alpha_self_thinning(ipts,ivm), write_debug, &
                              tot_res_target(ipts,ivm,:),             tot_lab_target(ipts,ivm,:),    gpp_week(ipts,ivm), &
                              t2m(ipts))



                         !! 5.2.3.3.1 Merge biomass of new and already available trees
                         ! Unlike LCC, we are only at this point if we are dealing
                         ! with forests, so we don't have to check for that here.

                         ! Copy the number of individuals and the biomass of the established 
                         ! vegetation to a temporary variable. In the subroutine ::circ_class_n
                         ! and ::circ_class_biomass will be overwritten with the characteristics
                         ! of the merged vegetation
                         ! The term "established" here refers to the youngest age class, while
                         ! "tmp" is the prescribed vegetation that we just created.
                         est_circ_class_n(:) = circ_class_n(ipts,iyoung,:)
                         est_circ_class_biomass(:,:,:) = circ_class_biomass(ipts,iyoung,:,:,:)
                         tmp_circ_class_n(:) = new_circ_class_n(ipts,ivm,:)
                         tmp_circ_class_biomass(:,:,:) = new_circ_class_biomass(ipts,ivm,:,:,:)
                         
                         ! Debug
                         IF(printlev_loc>=4)THEN
                            WRITE(numout,*) 'rest, ', SUM(circ_class_n(ipts,ivm,:)), &
                              SUM(est_circ_class_n), SUM(tmp_circ_class_n), &
                              SUM(SUM(circ_class_biomass(ipts,ivm,:,:,icarbon),2),1), &
                              SUM(est_circ_class_biomass), &
                              SUM(tmp_circ_class_biomass), &
                              ipts, iyoung, SUM(est_circ_class_dia), SUM(tmp_circ_class_dia), &
                              lpft_replant(ipts,ivm)
                         ENDIF
                         !-

                         ! Merge the biomass: new pft into the established pft
                         ! (ivm not iyoung)
                         CALL merge_biomass_pfts(npts, share_rec, circ_class_n, &
                              est_circ_class_n, tmp_circ_class_n, circ_class_biomass, &
                              est_circ_class_biomass, &
                              tmp_circ_class_biomass, &
                              ipts, iyoung, ivm, forest_managed, bm_to_litter)


                         ! As saplings are planted in the current PFT, atm_to_bm should be considered
                         ! when moving the PFT to the youngest PFT. 
                         ! Currently, there is a line that sets atm_to_bm to zero in stomate_lpj after 
                         ! LCC 'proportional', so below adding is not necessary for now. 
                         ! However, this line is added to address potential errors in the future.
                         atm_to_bm(ipts,ivm,:) = atm_to_bm(ipts,ivm,:) + new_atm_to_bm(ipts,ivm,:)

                         !+++CHECK+++
                         ! sugar_load may be duplicated -> see CHECK above
                         ! when_growthinit might get duplicxated as well -> see CHECK above
                         ! Note that we do not use the
                         ! "bank" concept present in LCC because we are replanting
                         ! the same PFT which already existed, just in a younger class.
                         CALL move_pft_properties(ipts, iyoung, ivm, share_rec, &
                              litter_weight_rec, circ_class_biomass, fm_change_map, &
                              wstress_season(ipts,iyoung), wstress_season(ipts,ivm), &
                              wstress_month(ipts,iyoung), wstress_month(ipts,ivm), &
                              season_drought_legacy(ipts,iyoung,:), &
                              season_drought_legacy(ipts,ivm,:), &
                              nstress_season(ipts,iyoung), nstress_season(ipts,ivm), &
                              n_input(ipts,iyoung,:,:), n_input(ipts,ivm,:,:), &
                              n_input_daily(ipts,iyoung,:), n_input_daily(ipts,ivm,:), &
                              vegstress_season(ipts,iyoung), vegstress_season(ipts,ivm), &
                              vegstress_month(ipts,iyoung), vegstress_month(ipts,ivm), &
                              vegstress_week(ipts,iyoung), vegstress_week(ipts,ivm), &
                              vegstress(ipts,iyoung), vegstress(ipts,ivm), &
                              humrel(ipts,iyoung), humrel(ipts,ivm), &
                              everywhere(ipts,iyoung), everywhere(ipts,ivm), &
                              PFTpresent(ipts,iyoung), PFTpresent(ipts,ivm), &
                              sugar_load(ipts,iyoung), sugar_load(ipts,ivm), &
                              cn_leaf_min_season(ipts,iyoung), cn_leaf_min_season(ipts,ivm), &
                              cn_leaf_init_2D, &
                              Light_Abs_Tot(ipts,iyoung,:), Light_Abs_Tot(ipts,ivm,:), &
                              Light_Tran_Tot(ipts,iyoung,:), Light_Tran_Tot(ipts,ivm,:), &
                              laieff_isotrop(ipts,:,iyoung), laieff_isotrop(ipts,:,ivm), &
                              lm_lastyearmax(ipts,iyoung), lm_lastyearmax(ipts,ivm), &
                              lm_thisyearmax(ipts,iyoung), lm_thisyearmax(ipts,ivm), &
                              age(ipts,iyoung), age(ipts,ivm), &
                              leaf_frac(ipts,iyoung,:), leaf_frac(ipts,ivm,:), &
                              leaf_age(ipts,iyoung,:), leaf_age(ipts,ivm,:), &
                              light_tran_to_floor_season(ipts,iyoung), &
                              light_tran_to_floor_season(ipts,ivm), &
                              qsintveg(ipts,iyoung), qsintveg(ipts,ivm), &
                              us(ipts,iyoung,:,:), us(ipts,ivm,:,:), &
                              when_growthinit(ipts,iyoung), when_growthinit(ipts,ivm), &
                              gdd_from_growthinit(ipts,iyoung), gdd_from_growthinit(ipts,ivm), &
                              gdd_midwinter(ipts,iyoung), gdd_midwinter(ipts,ivm), &
                              time_hum_min(ipts,iyoung), time_hum_min(ipts,ivm), &
                              hum_min_dormance(ipts,iyoung), hum_min_dormance(ipts,ivm), &
                              gdd_m5_dormance(ipts,iyoung), gdd_m5_dormance(ipts,ivm), &
                              ncd_dormance(ipts,iyoung), ncd_dormance(ipts,ivm), &
                              ngd_minus5(ipts,iyoung), ngd_minus5(ipts,ivm), &
                              mean_start_gs(ipts,iyoung), mean_start_gs(ipts,ivm), &
                              maxvegstress_lastyear(ipts,iyoung), maxvegstress_lastyear(ipts,ivm), &
                              maxvegstress_thisyear(ipts,iyoung), maxvegstress_thisyear(ipts,ivm), &
                              minvegstress_lastyear(ipts,iyoung), minvegstress_lastyear(ipts,ivm), &
                              minvegstress_thisyear(ipts,iyoung), minvegstress_thisyear(ipts,ivm), &
                              maxgppweek_lastyear(ipts,iyoung), maxgppweek_lastyear(ipts,ivm), &
                              maxgppweek_thisyear(ipts,iyoung), maxgppweek_thisyear(ipts,ivm), &
                              n_reserve_longterm(ipts,iyoung), n_reserve_longterm(ipts,ivm), &
                              n_reserve_balance(ipts,iyoung), n_reserve_balance(ipts,ivm), &
                              maxfpc_lastyear(ipts,iyoung), maxfpc_lastyear(ipts,ivm), &
                              maxfpc_thisyear(ipts,iyoung), maxfpc_thisyear(ipts,ivm), &
                              turnover_longterm(ipts,iyoung,:,:), turnover_longterm(ipts,ivm,:,:), &
                              dead_leaves(ipts,iyoung,:), dead_leaves(ipts,ivm,:), &
                              grow_season_len(ipts,iyoung), grow_season_len(ipts,ivm), &
                              KF(ipts,iyoung), KF(ipts,ivm), &
                              atm_to_bm(ipts,iyoung,:), atm_to_bm(ipts,ivm,:), &
                              npp_longterm(ipts,iyoung), npp_longterm(ipts,ivm), &
                              croot_longterm(ipts,iyoung), croot_longterm(ipts,ivm), &
                              gpp_daily(ipts,iyoung), gpp_daily(ipts,ivm), &
                              resp_maint(ipts,iyoung), resp_maint(ipts,ivm), &
                              resp_growth(ipts,iyoung), resp_growth(ipts,ivm), &
                              resp_hetero(ipts,iyoung), resp_hetero(ipts,ivm), &
                              npp_daily(ipts,iyoung), npp_daily(ipts,ivm), &
                              rue_longterm(ipts,iyoung), rue_longterm(ipts,ivm), &
                              leaching_daily(ipts,iyoung,:), leaching_daily(ipts,ivm,:), &
                              emission_daily(ipts,iyoung,:), emission_daily(ipts,ivm,:), &
                              gpp_week(ipts,iyoung), gpp_week(ipts,ivm), &
                              resp_maint_week(ipts,iyoung), resp_maint_week(ipts,ivm), &
                              gpp_year(ipts,iyoung), gpp_year(ipts,ivm), &
                              gpp_decade(ipts,iyoung), gpp_decade(ipts,ivm), &
                              mai(ipts,iyoung), mai(ipts,ivm), &
                              pai(ipts,iyoung), pai(ipts,ivm), &
                              mai_count(ipts,iyoung), mai_count(ipts,ivm), &
                              previous_wood_volume(ipts,iyoung), previous_wood_volume(ipts,ivm), &
                              age_stand(ipts,iyoung), age_stand(ipts,ivm), &
                              last_cut(ipts,iyoung), last_cut(ipts,ivm), &
                              k_latosa_adapt(ipts,iyoung), k_latosa_adapt(ipts,ivm), &
                              litter(ipts,:,iyoung,:,:), litter(ipts,:,ivm,:,:), &
                              bm_to_litter(ipts,iyoung,:,:), bm_to_litter(ipts,ivm,:,:), &
                              tree_bm_to_litter(ipts,iyoung,:,:), tree_bm_to_litter(ipts,ivm,:,:), &
                              leaf_age_crit(ipts,iyoung), leaf_age_crit(ipts,ivm), &
                              leaf_classes(ipts,iyoung), leaf_classes(ipts,ivm), &
                              co2_fire(ipts,iyoung), co2_fire(ipts,ivm), &
                              litterfuel(ipts,:,iyoung,:,:), litterfuel(ipts,:,ivm,:,:), &
                              turnover_daily(ipts,iyoung,:,:), turnover_daily(ipts,ivm,:,:), &
                              soil_n_min(ipts,iyoung,:), soil_n_min(ipts,ivm,:), &
                              p_O2(ipts,iyoung), p_O2(ipts,ivm), &
                              bact(ipts,iyoung), bact(ipts,ivm), &
                              deepSOM_a(ipts,:,iyoung,:), deepSOM_a(ipts,:,ivm,:), &
                              deepSOM_s(ipts,:,iyoung,:), deepSOM_s(ipts,:,ivm,:), &
                              deepSOM_p(ipts,:,iyoung,:), deepSOM_p(ipts,:,ivm,:), &
                              som(ipts,:,iyoung,:), som(ipts,:,ivm,:), &
                              lignin_struc(ipts,iyoung,:), lignin_struc(ipts,ivm,:), &
                              lignin_wood(ipts,iyoung,:), lignin_wood(ipts,ivm,:), &
                              lignin_snag(ipts,iyoung,:), lignin_snag(ipts,ivm,:), &
                              forest_managed(ipts,iyoung), forest_managed(ipts,ivm), &
                              plant_status(ipts,iyoung), plant_status(ipts,ivm), &
                              kill_vessels(ipts,iyoung,:), kill_vessels(ipts,ivm,:), &
                              vessel_loss_previous(ipts,iyoung,:), vessel_loss_previous(ipts,ivm,:), &
                              biomass_init_drought(ipts,iyoung,:,:,:), biomass_init_drought(ipts,ivm,:,:,:), &
                              count_daylight(ipts,iyoung), count_daylight(ipts,ivm), &
                              CN_som_litter_longterm(ipts,iyoung,:), CN_som_litter_longterm(ipts,ivm,:), &
                              matrixV(ipts,iyoung,:,:), matrixV(ipts,ivm,:,:), &
                              matrixA(ipts,iyoung,:,:), matrixA(ipts,ivm,:,:), &
                              vectorU(ipts,iyoung,:), vectorU(ipts,ivm,:), &
                              vectorB(ipts,iyoung,:), vectorB(ipts,ivm,:), &
                              gap_area_save(ipts,iyoung,:), gap_area_save(ipts,ivm,:), &
                              total_ba_init(ipts,iyoung), total_ba_init(ipts,ivm))
                         !+++++++++++
                         
                         ! Now move the veget_max from the current class to the young one.
                         veget_max(ipts,iyoung)=veget_max(ipts,iyoung)+veget_max(ipts,ivm)
                         veget_max(ipts,ivm) = zero

                      ENDIF
                
                      !+++CHECK+++
                      ! Should this be done in another place?
                      circ_class_biomass(ipts,ivm,:,:,:) = zero
                      circ_class_n(ipts,ivm,:) = zero
                      !+++++++++++

                   ENDIF ! If this is the youngest age class

                ENDIF ! If we are using age classes

             ENDIF ! All biomass in ivm and plant the same species

          ENDIF ! is_tree

       ENDDO  ! loop over pfts

    END DO ! loop over land points


 !! 3. Check numerical consistency of this routine

    IF (err_act.GT.1) THEN

       ! 3.2 Check surface area
       CALL check_vegetation_area("mortality_clean", npts, veget_max_begin, &
            veget_max,'ageclass')
       
       ! 3.3 Mass balance closure 
       ! 3.3.1 Calculate final biomass
       pool_end(:,:,:) = zero
       
       DO iele = 1,nelements
          DO ipar = 1,nparts
             DO icir = 1,ncirc
                pool_end(:,:,iele) = pool_end(:,:,iele) + &
                     (circ_class_biomass(:,:,icir,ipar,iele) * &
                     circ_class_n(:,:,icir) * veget_max(:,:)) 
             ENDDO
             ! bm_to_litter
             pool_end(:,:,iele) = pool_end(:,:,iele) + &
                  bm_to_litter(:,:,ipar,iele) * veget_max(:,:)
          ENDDO

          ! Litter pool (gC m-2) *  (m2 m-2) 
          DO ilit = 1,nlitt
             DO ilev = 1,nlevs
                pool_end(:,:,iele) = pool_end(:,:,iele) + &
                     litter(:,ilit,:,ilev,iele) * veget_max(:,:)
             ENDDO
          ENDDO

          IF (ok_soil_carbon_discretization) THEN
             ! Soil carbon (gC m-3) * (m2 m-2)
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

       ! Denitrifying bacteria (only C)
       pool_end(:,:,icarbon) = pool_end(:,:,icarbon) + &
            bact(:,:) * veget_max(:,:)

       ! Nitrogen only
       DO ispec = 1,nnspec
          pool_end(:,:,initrogen) = pool_end(:,:,initrogen) + &
               soil_n_min(:,:,ispec) * veget_max(:,:)
       ENDDO
    
       ! 3.3.2 Calculate mass balance
       ! Common processes
       DO iele=1,nelements
          check_intern(:,:,iatm2land,iele) = &
               check_intern(:,:,iatm2land,iele) + &
               atm_to_bm(:,:,iele) * veget_max(:,:) 
          check_intern(:,:,ipoolchange,iele) = -un * (pool_end(:,:,iele) - &
               pool_start(:,:,iele)) 
       ENDDO
   
       ! Specific fluxes
       ! The fluxes should not be accounted for here because this
       ! code deals with mortality. Following mortality the veget_max
       ! has to be moved to the youngest age class but the fluxes should 
       ! stay in the age class for which they were generated. Clearly, 
       ! closing the mass balance with the fluxes, requires moving
       ! the fluxes. If done so, the NBP consistency check will fail
       ! because those fluxes are double counted: once after iage
       ! in the old age class and once after imor in the youngest
       ! age class. NThe mass balance is still checked for the pools
       closure_intern(:,:,:) = zero
       DO imbc = 1,nmbcomp
          DO iele=1,nelements
             ! Debug
             IF (printlev_loc>=4) WRITE(numout,*) &
                  'check_intern, ivm, imbc, iele, ',test_pft, imbc, &
                  iele, check_intern(:,test_pft,imbc,iele)
             !-
             closure_intern(:,:,iele) = closure_intern(:,:,iele) + &
                  check_intern(:,:,imbc,iele)
          ENDDO
       ENDDO

       ! 4.3.3 Check mass balance closure
       CALL check_mass_balance("mortality_clean", closure_intern, npts, &
            pool_end, pool_start, veget_max, 'ageclass')
       
    ENDIF ! err_act.GT.1
    
    IF (printlev.GE.4) WRITE(numout,*) 'Leaving mortality_clean'
    
  END SUBROUTINE mortality_clean

!!
!================================================================================================================================
!! SUBROUTINE  : disturbance_kill
!!
!>\BRIEF       
!!
!! DESCRIPTION : Calculates total mortality from the disturbances and decides
! to replace the stand or not
!!              If the total mortality is higher than the threshold then the
!veget_max_disturb is calculated 
!!              to call land_cover_change.
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): :: veget_max_disturb,loss_gain
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_
!================================================================================================================================

  SUBROUTINE disturbance_kill(npts, circ_class_biomass, circ_class_n, veget_max, loss_gain, &
                              veget_max_disturb,circ_class_kill)
       
 !! 0. Variable and parameter description

    !! 0.1 Input variables
    INTEGER, INTENT(in)                             :: npts                 !! Domain size - number of pixels (unitless,0-1) 
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(in)   :: circ_class_biomass   !! Biomass of the componets of the model tree
                                                                            !! within a circumference
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)       :: circ_class_n         !! Number of trees within each circumference 
                                                                            !! @tex $(indm^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(in)         :: veget_max            !! Cover fraction of a PFT (unitless, 0-1)

    !! 0.2 Output variables
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)   :: loss_gain            !! The same as delta_veget but distributed of all
                                                                            !! age classes and thus taking the age-classes into 
                                                                            !! account (unitless, 0-1)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)   :: veget_max_disturb    !! New cover fraction of a PFT after stand-replacing
                                                                            !! disturbances (unitless, 0-1)

    !! 0.3 Modified variables
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout):: circ_class_kill      !! Number of trees within a circ that needs
                                                                            !! to be killed @tex $(indm^{-2})$ @endtex

    !! 0.4 Local variables
    INTEGER(i_std)                                  :: ipts, ivm            !! Indice
    INTEGER(i_std)                                  :: icut, igroup, iyoung !! Indice
    REAL(r_std), DIMENSION(npts,nvmap)              :: loss_gain_group      !! Total veget_max of all age clases within a PFT group 
                                                                            !! (unitless, 0-1)
    REAL(r_std), DIMENSION(npts,nvm)                :: total_damage_disturb !! The total rate of damage by natural disturbances
    REAL(r_std)                                     :: damage_rate_wind     !! 
    REAL(r_std)                                     :: damage_rate_beetle   !!
    REAL(r_std)                                     :: damage_rate_fire
    REAL(r_std)                                     :: temp_dens_target     !! Temporary and local target density.
!_
!================================================================================================================================

    IF (printlev >= 2) WRITE(numout,*) 'Entering disturbance_kill'

    !! 1.1 Set firstcall flag
    IF ( firstcall_stomate_kill ) THEN
       !! Initialize local printlev
      ! Guillaume M. -- get_printlev calls getin_p, a COLLECTIVE bcast. Called here,
      ! inside the compute loop under a local firstcall, it hangs the run as soon as the
      ! ranks stop reaching it together. Use printlev, read once at startup by every rank.
       printlev_loc=printlev
       firstcall_stomate_kill = .FALSE.
    ENDIF

    !! 1.2  Initialize
    veget_max_disturb = veget_max
    loss_gain = zero
    loss_gain_group = zero
    total_damage_disturb = zero
  
    !! 2. Calculate the total mortality from the natural disturbances

    ! 2.1 total kill cannot exceed the number of individuals 
    DO ipts = 1, npts  
      DO ivm = 1, nvm
        !! 3. Replacing the stand or not
        ! If the total mortality is less than the threshold, it is considered
        ! to create small gaps. 
        ! But if it is over the threshold, then the dead stand will be replanted.
        ! This is not applied when the pft is the youngest.
        
        damage_rate_wind = zero
        damage_rate_beetle = zero
        damage_rate_fire = zero
        ! Currently, the disturbances are calculated with trees
        IF(.NOT. is_tree(ivm)) CYCLE
    
        igroup=agec_group(ivm)
        iyoung=start_index(igroup)
    
        IF ((iyoung + nagec_pft(igroup) -1) .EQ. iyoung) THEN
        ! There is only one age class. This should be dealt again later.
        ELSE IF (ivm .NE. iyoung) THEN
            ! In modeled heterogeneous stand with size structures, the ratio of biomass
            ! loss is a better indicator than the number of individuals since we are 
            ! interested more in a few dominant big trees than many small understory trees  

            IF(SUM(circ_class_n(ipts,ivm,:)*SUM(circ_class_biomass(ipts,ivm,:,:,icarbon),2)) .GT. &
                        min_stomate) THEN

                damage_rate_wind = &  
                   MIN( SUM( SUM( (circ_class_kill(ipts,ivm,:,:,icut_storm_break) + &
                                  circ_class_kill(ipts,ivm,:,:,icut_storm_uproot)),2 ) * &
                        SUM(circ_class_biomass(ipts,ivm,:,:,icarbon),2) )/ &
                        SUM(circ_class_n(ipts,ivm,:)*SUM(circ_class_biomass(ipts,ivm,:,:,icarbon),2) ), un)
                damage_rate_beetle = &
                   MIN( SUM( SUM(circ_class_kill(ipts,ivm,:,:,icut_beetle),2) * &
                        SUM(circ_class_biomass(ipts,ivm,:,:,icarbon),2) )/ &
                        SUM(circ_class_n(ipts,ivm,:)*SUM(circ_class_biomass(ipts,ivm,:,:,icarbon),2) ), un)
                damage_rate_fire = &
                   MIN( SUM( SUM(circ_class_kill(ipts,ivm,:,:,icut_fire),2) * &
                        SUM(circ_class_biomass(ipts,ivm,:,:,icarbon),2) )/ &
                        SUM(circ_class_n(ipts,ivm,:)*SUM(circ_class_biomass(ipts,ivm,:,:,icarbon),2) ), un)
            ENDIF

           ! This code combines wind and bark beetle damage rates into a new cohort,even if only one disturbance 
           ! exceeds the threshold. While this approach simplifies calculations significantly, it does not
           ! differentiate between mass damage from a single disturbance and smaller, gap-creating damage  
           ! from other disturbances. 
           IF ( damage_rate_wind .GE. replace_threshold_wind .OR. &
                 damage_rate_beetle .GE. replace_threshold_beetle .OR. &
                 damage_rate_fire .GE. replace_threshold_fire) THEN
               total_damage_disturb(ipts,ivm) = MIN(damage_rate_wind + damage_rate_beetle + damage_rate_fire, un)

               ! Guillaume M. -- Route the stand-replacing area through CLASS 0 instead of class 1:
               ! otherwise only the progressive harvest feeds the register and fire, storm and
               ! beetle skip the regeneration phase. The three agents are booked TOGETHER on the
               ! summed rate, as the condition above treats them. /!\ This runs after
               ! age_class_distr and the conveyor, so the credit is used on the NEXT annual pass.
               IF (ok_disturb_via_class0 .AND. n_class0 > 0 .AND. ALLOCATED(class0_veg)) THEN
                  class0_veg(ipts,agec_group(ivm),1) = class0_veg(ipts,agec_group(ivm),1) + &
                       total_damage_disturb(ipts,ivm) * veget_max(ipts,ivm)
               ENDIF
           ENDIF

           ! If the individuals died more than dens target then let the pft die
           ! through the killing process
           !IF ( total_damage_disturb .GE. dens_target(ivm)/SUM(circ_class_n(ipts,ivm,:)) THEN
           !     total_damage_disturb(ipts,ivm) = zero
           !ENDIF

           IF(printlev_loc .GE. 4 .AND. ipts .EQ. test_grid .AND. SUM(total_damage_disturb(ipts,:)) .GT. zero) THEN
               WRITE(numout,*) 'ivm', ivm
               WRITE(numout,*) 'total_damage_disturb',total_damage_disturb(ipts,ivm)
               WRITE(numout,*) 'damage_rate_beetle',damage_rate_beetle 
               WRITE(numout,*) 'damage_rate_wind',damage_rate_wind
               WRITE(numout,*) 'damage_rate_fire',damage_rate_fire
               WRITE(numout,*) '3',SUM( SUM((circ_class_kill(ipts,ivm,:,:,icut_storm_break) + &
                                  circ_class_kill(ipts,ivm,:,:,icut_storm_uproot)),2) * &
                        SUM(circ_class_biomass(ipts,ivm,:,:,icarbon),2) )
               WRITE(numout,*) '4', SUM(circ_class_n(ipts,ivm,:)*SUM(circ_class_biomass(ipts,ivm,:,:,icarbon),2) )
           ENDIF

           IF(total_damage_disturb(ipts,ivm) .GT. min_stomate) THEN
           
                veget_max_disturb(ipts,ivm) = veget_max(ipts,ivm) - veget_max(ipts,ivm) * &
                    total_damage_disturb(ipts,ivm)  

                IF(veget_max_disturb(ipts,ivm) .LT. min_vegfrac .OR. &
                    veget_max(ipts,ivm) - veget_max_disturb(ipts,ivm) .LT. min_vegfrac ) THEN
                    veget_max_disturb(ipts,ivm) = veget_max(ipts,ivm)
                ELSE
                   ! Not to apply it again, empty the circ_class_kill
                   circ_class_kill(ipts,ivm,:,:,icut_storm_uproot) = zero
                   circ_class_kill(ipts,ivm,:,:,icut_storm_break) = zero
                   circ_class_kill(ipts,ivm,:,:,icut_beetle) = zero
                   circ_class_kill(ipts,ivm,:,:,icut_fire) = zero
                ENDIF

                ! If remained veget frac is less than min_vegfrac then move all 
                IF(veget_max_disturb(ipts,ivm) .LE. min_vegfrac) THEN
                ENDIF
                ! Update veget_max_disturb for the youngest age class
                veget_max_disturb(ipts,iyoung) =  veget_max_disturb(ipts,iyoung) + &
                    veget_max(ipts,ivm) - veget_max_disturb(ipts,ivm)  
                ! Calculate loss gain
                loss_gain(ipts,ivm) = veget_max_disturb(ipts,ivm) - veget_max(ipts,ivm) 
                loss_gain(ipts,iyoung) = veget_max_disturb(ipts,iyoung) - veget_max(ipts,iyoung)
                loss_gain_group(ipts,igroup) = loss_gain_group(ipts,igroup) + &
                    loss_gain(ipts,ivm)
                
            ENDIF

            ! The calculation below should be performed even if there is no biomass in the  
            ! oldest age class to verify the accuracy of the loss_gain calculation.  
            ! Please suggest improvements if a better approach exists.  
            IF(ivm .EQ. (iyoung + nagec_pft(igroup) -1)) THEN ! oldest
                loss_gain_group(ipts,igroup) = loss_gain_group(ipts,igroup) + loss_gain(ipts,iyoung)
             ENDIF
             
        ENDIF !IF iyoung
       
        IF(printlev_loc .GE. 4 .AND. ipts .EQ. test_grid .AND. loss_gain(ipts,ivm) .GT. zero) THEN
            WRITE(numout,*) 'ipts, ivm',ipts, ivm
            WRITE(numout,*) 'igroup,iyoung',igroup,iyoung
            WRITE(numout,*) 'loss_gain',loss_gain(ipts,ivm)
            WRITE(numout,*) 'loss_gain_y',loss_gain(ipts,iyoung)
            WRITE(numout,*) 'rate_wind',damage_rate_wind
            WRITE(numout,*) 'rate_beete',damage_rate_beetle
            WRITE(numout,*) 'rate_fire',damage_rate_fire
            WRITE(numout,*) 'total_rate',total_damage_disturb(ipts,ivm)
            WRITE(numout,*) 'sum_ccn',SUM(circ_class_n(ipts,ivm,:))
            WRITE(numout,*) 'veget_max',veget_max(ipts,ivm)
            WRITE(numout,*) 'veget_max_disturb',veget_max_disturb(ipts,ivm)
        ENDIF                           
      ENDDO ! ivm    
    ENDDO ! ipts

    !! 4. error check
    DO ipts = 1,npts 
      IF(SUM(loss_gain_group(ipts,:)) .GT. 10*EPSILON(un) .OR. &
         SUM(loss_gain_group(ipts,:)) .LT. -10*EPSILON(un)  ) THEN
        ! Since veget_max removed added to the youngest PFT, the total loss_gain
        ! should be zero
        ! Error checking is strict and the single precision criterion
        ! has been violated. Stop the model.
        WRITE(numout,*) "ERROR: sum of loss_gain_group is not zero"
        WRITE(numout,*) 'ipts, loss_gain sum', ipts, SUM(loss_gain_group(ipts,:))
        CALL ipslerr_p (3,'disturbance_kill: Check_loss_gain', '', &
               '', '')
      ENDIF
    ENDDO   
  END SUBROUTINE
END MODULE stomate_kill
