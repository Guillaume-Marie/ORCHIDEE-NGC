! =================================================================================================================================
! MODULE       : stomate_mark_kill
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF         Marks vegetation for killing by natural causes.  Replaces
!                lpj_kill.  Does not work with the DGVM.
!!
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S) :
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_stomate/stomate_mark_kill.f90 $
!! $Date: 2026-01-29 14:53:46 +0100 (jeu. 29 janv. 2026) $
!! $Revision: 9352 $
!! \n
!_ ================================================================================================================================


MODULE stomate_mark_kill

  ! modules used:

  USE ioipsl_para
  USE xios_orchidee
  USE stomate_data
  USE pft_parameters
  USE constantes
  USE function_library,       ONLY : nmax, wood_to_qmdia, check_vegetation_area, &
                                     check_mass_balance, biomass_to_lai, &
                                     get_printlev, weibull_class_dist, &
                                     wood_to_qmdia_up_half, wood_to_dia, &
                                     calculate_rdi_boundaries

  IMPLICIT NONE

  ! private & public routines

  PRIVATE
  PUBLIC mark_to_kill, distribute_mortality_biomass

  INTEGER(i_std), SAVE       :: printlev_loc       !! Local level of text output for this module
!$OMP THREADPRIVATE(printlev_loc)
  LOGICAL, SAVE              :: firstcall_mark_kill = .TRUE.   !! first call flag
!$OMP THREADPRIVATE(firstcall_mark_kill)
CONTAINS

!! ================================================================================================================================
!! SUBROUTINE  : mark_to_kill
!!
!>\BRIEF       Kills natural pfts that have low biomass or low number of individuals.
!!             Does not change the biomass or litter pools.
!!
!! DESCRIPTION : Kills natural PFTS.  Forests can die through self-thinning or
!!             environmental conditions.  Crops can die through environmental
!!             conditions, and grasses don't die at all.  The total number
!!             of individuals to be killed in each circ class are outputted
!!             in circ_class_kill, and the biomass is actually killed in
!!             gap_prognostic.
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::circ_class_kill
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE mark_to_kill (npts, whichroutine, lm_lastyearmax, &
       PFTpresent, npp_longterm, circ_class_biomass, &
       circ_class_n, circ_class_kill, rue_longterm, turnover_longterm, &
       dt, veget_max, forest_managed, tot_res_target, &
       tot_lab_target, plant_status)

 !! 0. Variable and parameter description

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                                :: npts                !! Domain size (unitless)
    CHARACTER(LEN=10), INTENT(in)                             :: whichroutine        !! Message (unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(in)                   :: veget_max           !! "maximal" coverage fraction of a PFT on
                                                                                     !! the ground (unitless, 0-1)
    REAL(r_std), DIMENSION(:,:), INTENT(in)                   :: lm_lastyearmax      !! Last year's maximum leaf mass, for each
                                                                                     !! PFT @tex $(gC.m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(in)                   :: rue_longterm        !! Longterm radiation use efficiency
                                                                                     !! (??units??)
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)               :: turnover_longterm   !! "Long term" (default 3-year) turnover
                                                                                     !! rate @tex $(gC m^{-2} year^{-1})$ @endtex
    REAL(r_std), INTENT(in)                                   :: dt                  !! Time step (days)
    REAL(r_std), DIMENSION(:,:), INTENT(in)                   :: npp_longterm        !! "Long term" (default = 3-year) net primary
                                                                                     !! productivity
                                                                                     !! @tex $(gC.m^{-2} year^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:,:),INTENT(in)                  :: tot_res_target      !! Target for reserve carbon pool (only for grassland) 
                                                                                     !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:),INTENT(in)                  :: tot_lab_target      !! Target for labile carbon pool (only for grassland)
                                                                                     !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(in)                   :: plant_status        !! Growth and phenological status of the plant

    !! 0.2 Output variables


    !! 0.3 Modified variables

    LOGICAL, DIMENSION(:,:), INTENT(inout)                    :: PFTpresent          !! Is pft there (true/false)
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)          :: circ_class_biomass  !! Biomass of the componets of the model 
                                                                                     !! tree within a circumference
                                                                                     !! class @tex $(gC ind^{-1})$ @endtex 
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)              :: circ_class_n        !! Number of individuals in each circ class
                                                                                     !! @tex $(m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)          :: circ_class_kill     !! Number of trees within a circ that needs
                                                                                     !! to be killed @tex $(ind m^{-2})$ @endtex
    INTEGER(i_std), DIMENSION (:,:), INTENT(in)               :: forest_managed      !! forest management flag (is the forest
                                                                                     !! being managed?)

    !! 0.4 Local variables

    INTEGER(i_std)                                            :: ipts, ivm,icir,ipar !! Indices       (unitless)
    INTEGER(i_std)                                            :: iele, imbc, istep   !! Indices       (unitless)
    INTEGER(i_std)                                            :: ifm, iout, istart   !! Indices       (unitless)
    LOGICAL, DIMENSION(npts)                                  :: was_killed          !! Bookkeeping   (true/false)
    REAL(r_std), DIMENSION(ncirc)                             :: diameters_temp      !! Trunk diameter calculated from allometry
    REAL(r_std)                                               :: Dg                  !! Quadratic mean diameter of the stand (cm)
    REAL(r_std)                                               :: Dg_up_half          !! Target quadratic mean diameter of the stand (cm)
    REAL(r_std)                                               :: bm_difference       !! the difference between self-thinning biomass
                                                                                     !! and that due to environmental mortality
    REAL(r_std)                                               :: difference          !! the difference between the self-thinning
                                                                                     !! density
                                                                                     !! and the true density (individuals per m**2)
    REAL(r_std), DIMENSION(npts,nvm)                          :: rdi                 !! Relative density index (unitless, 0-2)
    REAL(r_std), DIMENSION(npts,nvm)                          :: rdi_target_upper    !! Upper limit of RDI. When reached the stand
                                                                                     !! needs to be thinned (managed) / kille (unmanaged) (unitless)   
    REAL(r_std), DIMENSION(npts,nvm)                          :: rdi_target_lower    !! Lower limit of RDI. When thin/kill, thin/kill to this             
                                                                                     !! stand density (unitless)
    REAL(r_std), DIMENSION(npts,nvm)                          :: st_dens             !! the self-thinning density of the stand
                                                                                     !! (individuals per m**2)
    REAL(r_std), DIMENSION(npts,nvm)                          :: st_mortality        !! biomass killed due to self-thinning density
    REAL(r_std), DIMENSION(npts,nvm)                          :: d_mortality         !! biomass killed due to the environment
    REAL(r_std)                                               :: delta_biomass       !! Net biomass increase for the previous
                                                                                     !! year @tex $(gC m^{-2} year^{-1})$ @endtex
    REAL(r_std)                                               :: vigour              !! Growth efficiency, an indicator of tree
                                                                                     !! vitality, used to calculate mortality
    REAL(r_std)                                               :: mortality_greff     !! Mortality rate derived by growth
                                                                                     !! efficiency @tex $(year^{-1})$ @endtex
    REAL(r_std)                                               :: mortality           !! Mortality (fraction of trees that is
                                                                                     !! dying per time step)
    REAL(r_std), DIMENSION(npts,nvm,nmbcomp,nelements)        :: check_intern        !! Contains the components of the internal
                                                                                     !! mass balance chech for this routine
                                                                                     !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)                :: closure_intern      !! Check closure of internal mass balance
                                                                                     !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)                :: pool_start, pool_end!! Start and end pool of this routine
                                                                                     !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    LOGICAL                                                   :: lconverged          !! Checking for loop convergence
    REAL(r_std), DIMENSION(ncirc)                             :: target_st_kill      !! How many trees from each circumference
                                                                                     !! class we would like to kill.
                                                                                     !! @tex $(trees m^{-2})$ @endtex
    REAL(r_std), DIMENSION(ncirc)                             :: true_st_kill        !! How many trees from each circumference
                                                                                     !! class that we actually killed.
                                                                                     !! @tex $(trees m^{-2})$ @endtex
    REAL(r_std)                                               :: excess_ind          !! THe difference between the amount of
                                                                                     !! trees we want to kill and how many we
                                                                                     !! can actually kill based on our
                                                                                     !! circ_class distribution.
                                                                                     !! @tex $(trees m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)                          :: veget_max_begin     !! temporary storage of veget_max to check area conservation
    REAL(r_std), DIMENSION(npts,nvm,ncirc)                    :: tmp                 !! Temporary variable to prepare information for xios
    REAL(r_std)                                               :: n_init              !! Initial number of individuals calculated
                                                                                     !! based on gradient calculation (ind m^(-2))
    REAL(r_std), DIMENSION(npts,nvm,noutdiaclass)             :: temp1, temp2        !! Temporary variable for writing to history files
    REAL(r_std), DIMENSION(ncirc)                             :: circ_dia            !! temporary variable for vegetation diameter
    REAL(r_std), DIMENSION(ncirc)                             :: st_dist             !! distribution of the cutting trees in unmanaged forest
    REAL(r_std)                                               :: Trl0                !! The sum of individual reserve and labile carbon before
                                                                                     !! density adjustment (gC ind-1)
    REAL(r_std)                                               :: Ttarget_area        !! The sum of areal target for reserve and labile carbon (gC m-2)
    REAL(r_std)                                               :: Cother              !! The sum of carbon from all compartments except for reserve and
                                                                                     !! labile pool (gC ind-1)
    REAL(r_std), DIMENSION(npts,nvm,ncirc)                    :: old_circ_class_n    !! old grass density (ind m-2)
    REAL(r_std), DIMENSION(npts,nvm,ncirc)                    :: new_circ_class_n    !! new grass density (ind m-2)
    REAL(r_std),PARAMETER                                     :: scale_dec = 0       !! Scaling factor to calculate reserve and labile carbon 
                                                                                     !! when decreasing grassland density (unitless)
                                                                                     !! Has been used to tune the results. Could be externalized 
                                                                                     !! and made PFT-specific when found to be impactful.
    REAL(r_std),PARAMETER                                     :: scale_carbon = 0    !! Scaling factor to control the percentage of carbon into litter
                                                                                     !! pool when decreasing grassland density (unitless)
                                                                                     !! Has been used to tune the results. Could be externalized 
                                                                                     !! and made PFT-specific when found to be impactful.    

!_ ================================================================================================================================
   
   IF (firstcall_mark_kill) THEN
      !! Initialize local printlev
      ! Guillaume M. -- get_printlev calls getin_p, a COLLECTIVE bcast. Called here,
      ! inside the compute loop under a local firstcall, it hangs the run as soon as the
      ! ranks stop reaching it together. Use printlev, read once at startup by every rank.
      printlev_loc=printlev
      firstcall_mark_kill=.FALSE.
   END IF

   IF (printlev_loc >= 2) WRITE(numout,*) 'Entering mark_to_kill'

 !! 1. Initialize

   !! 1.2 Initialize check for mass balance closure
   IF (err_act.GT.1) THEN

      check_intern = zero
      pool_start = zero
      DO ipar = 1,nparts
         DO iele = 1,nelements
            DO icir = 1,ncirc
             
               !  Initial biomass pool
               pool_start(:,:,iele) = pool_start(:,:,iele) + &
                    circ_class_biomass(:,:,icir,ipar,iele) * &
                    circ_class_n(:,:,icir) * veget_max(:,:)

            ENDDO
         ENDDO
      ENDDO
   
      !! 1.3 Initialize check for area conservation
       veget_max_begin(:,:) = veget_max(:,:)

    ENDIF ! err_act.GT.1

    !! 1.4 Initialize variables
    d_mortality(:,:)=zero
    st_mortality(:,:)=zero
   
 !! 2. Kill PFTs

   ! Kill plants if number of individuals or last year's leaf mass is close to zero.
   ! the "was_killed" business is necessary for a more efficient code on the VPP
   ! processor
   DO ivm = 2,nvm ! loop plant functional types

      was_killed(:) = .FALSE.
      IF ( natural(ivm) ) THEN
         
         !! 2.2 Kill natural PFTs when running STOMATE without DGVM

         !! 2.2.1 Kill PFTs
         !  Kill PFTs but not in the first year because one year
         !  is needed to fill the long-term and seasonal variables.
         !  Therefore no leaves are grown in the first year. As
         !  long as the plant has never grown leaves, ::rue_longterm
         !  keeps its initial value of 1. Note that the logic of
         !  this statement differs from the original verison that
         !  was used for the resource limited allocation where the
         !  carbohydrate reserves pool was treated differently than
         !  in the functional allocation.
         WHERE ( PFTpresent(:,ivm) .AND. &
              (rue_longterm(:,ivm) .NE. un) .AND. & 
              (SUM(circ_class_biomass(:,ivm,:,icarbres,icarbon)*circ_class_n(:,ivm,:),2) .LE. min_stomate .AND. & 
              SUM(circ_class_biomass(:,ivm,:,ilabile,icarbon)*circ_class_n(:,ivm,:),2) .LE. min_stomate .AND. &
              SUM(circ_class_biomass(:,ivm,:,ileaf,icarbon)*circ_class_n(:,ivm,:),2) .LE. min_stomate ) .AND. & 
              SUM(circ_class_n(:,ivm,:),2) .GT. zero)

            was_killed(:) = .TRUE.
           
         ENDWHERE         
 
         !+++CHECK+++
         ! If the density of a forest goes below a certain level,
         ! the forest dies off. The density limit is referred to
         ! in Bellasen et al 2010. Because the dying was very abrupt
         ! resulting in a huge litter input in a single time step
         ! this cause of mortality has been moved in stomate_kill
         ! where the trees are made to die slowly (thus over the 
         ! course of several years) resulting in more realistic 
         ! litter inputs. The benefit of having this code here is
         ! that it makes use circ_cass_kill and therefore enables
         ! dealing with several mortality causes in a single time
         ! step. Ideally mortality caused by exceeding a density
         ! treshold should be taken out of stomate_kill and moved 
         ! backed in stomate_mark_kill.
         !+++++++++++

         !+++CHECK+++
         ! When running the model with a poor parameter set for PFT9, the
         ! wood turnover exceeded wood growth for a pixel in Greenland.
         ! This resulted in a Cs of zero which was causing mass balance
         ! problems in stomate_phenology. This should not happen with
         ! reasonable parameters but we check for it anyway. Ideally only
         ! the circumference class with a Cs = zero should be killed and
         ! mortality_clean should take care of it. As we are working towards
         ! a dealine and this is a rare case (happened after 355 years of a
         ! global simulation with 15 PFTs) we implemented a more quick and
         ! dirty solution.
         IF(is_tree(ivm))THEN

            DO icir = 1,ncirc
               
               WHERE ( circ_class_n(:,ivm,icir).GT.zero .AND. &
                    (circ_class_biomass(:,ivm,icir,isapabove,icarbon) + &
                    circ_class_biomass(:,ivm,icir,isapbelow,icarbon)) .LE. min_stomate)

                  was_killed(:) = .TRUE.
                  
               ENDWHERE

            END DO
            
         END IF
         !+++++++++++
         
         ! Debug
         IF(printlev_loc>4)THEN
            IF(is_tree(ivm))THEN
               DO ipts=1,npts
                  ifm=forest_managed(ipts,ivm)
                  IF( PFTpresent(ipts,ivm) .AND. &
                          (SUM(circ_class_n(ipts,ivm,:)) .LE. dens_target(ivm,ifm)*ha_to_m2))THEN
                     WRITE(numout,*) 'Killing PFT, ipts: ',ivm,ipts
                     WRITE(numout,*) 'Forest density too low'
                     WRITE(numout,*) 'circ_class_n(ipts,ivm,:)',circ_class_n(ipts,ivm,:)
                     WRITE(numout,*) 'dens_target(ivm,ifm)*ha_to_m2',dens_target(ivm,ifm)*ha_to_m2
                  ENDIF
               ENDDO
            ENDIF

            DO ipts=1,npts
               IF ( PFTpresent(ipts,ivm) .AND. &
                    (rue_longterm(ipts,ivm) .NE. un) .AND. & 
                    (SUM(circ_class_biomass(ipts,ivm,:,icarbres,icarbon)) .LE. min_stomate .AND. & 
                    SUM(circ_class_biomass(ipts,ivm,:,ilabile,icarbon)) .LE. min_stomate .AND. &
                    SUM(circ_class_biomass(ipts,ivm,:,ileaf,icarbon)) .LE. min_stomate ) .AND. & 
                    SUM(circ_class_n(ipts,ivm,:)) .GT. zero)THEN
                  WRITE(numout,*) 'Killing PFT, ipts: ',ivm,ipts
                  WRITE(numout,*) 'No roots, leaves, labile pool.'
                  WRITE(numout,*) 'circ_class_n(ipts,ivm,:)',SUM(circ_class_n(ipts,ivm,:))
                  WRITE(numout,*) 'Reserve pool',SUM(circ_class_biomass(ipts,ivm,:,icarbres,icarbon))
                  WRITE(numout,*) 'Labile pool',SUM(circ_class_biomass(ipts,ivm,:,ilabile,icarbon))
                  WRITE(numout,*) 'Leaves',SUM(circ_class_biomass(ipts,ivm,:,ileaf,icarbon))
               ENDIF
            ENDDO
         ENDIF
         !-

         !! 2.3 Bookkeeping for PFTs that were killed
         !  For PFTs that were killed, return biomass pools to litter
         !  and update biomass pools to zero         
         DO ipts = 1,npts

            IF ( was_killed(ipts) ) THEN

               ! If a PFT is killed, this just means that we mark all the
               ! individuals for a natural death (all debris stays on site).
               ! The killing actually takes place in lpj_gap.f90.
               ! We have to do this differently for trees and non-trees, since
               ! circ_class_n is not defined for non-trees.
               IF(is_tree(ivm))THEN

                  ! In reality, if the forest is managed we will not lose this
                  ! wood. The forest would be harvested. So, let's put these into
                  ! different pools depending on the management type.
                  ifm=forest_managed(ipts,ivm)

                  DO icir=1,ncirc
                     
                     circ_class_kill(ipts,ivm,icir,ifm,icut_clear)=&
                          circ_class_n(ipts,ivm,icir)

                     ! Debug
                     IF ((printlev_loc.GT.4) .AND. &
                             (circ_class_kill(ipts,ivm,icir,ifm,icut_clear) .GT. zero) &
                             .AND. (ipts .EQ. test_grid) .AND. (ivm .EQ. test_pft)) THEN

                             WRITE(numout,*) 'found missing mortality event', circ_class_kill(ipts,ivm,icir,ifm,icut_clear)
                     ENDIF
                     !
                  ENDDO

               ELSE

                  ! Since stomate_mortality only knows about circ_class_kill and circ_class_biomass,
                  ! we need to give values to icir=1 for these variables for grasses
                  ! and crops so that the correct amount of biomass is killed.
                  circ_class_kill(ipts,ivm,1,ifm_none,icut_clear) = &
                       SUM(circ_class_n(ipts,ivm,:))

               ENDIF ! is_tree(ivm)

               !! 2.7 Print sub routine messages
               IF (printlev_loc>=4) THEN

                  WRITE(numout,*) 'kill: eliminated ',PFT_name(ivm)
                  WRITE(numout,*) '  pft number: ',ivm
                  WRITE(numout,*) '  at grid: ',ipts

               ENDIF

            ENDIF ! was_killed

         ENDDO ! loop over points

      ENDIF ! PFT is natural

      !! We need to enforce two different types of killing at the moment,
      !! self-thinning and environmental mortality. Self-thinning is due
      !! to resouce competition in dense stands, and environmental
      !! mortality is due to environmental stress. Self-thinning is
      !! only applicable to forests. 
      DO ipts=1,npts

         ! Don't waste time if the whole PFT was already scheduled to die,
         ! or if the PFT is not here.
         IF(was_killed(ipts))CYCLE
         IF(.NOT.PFTpresent(ipts,ivm))CYCLE
 
         !! let's calculate the self-thinning limit first
         ! Number of indiviudals for trees and grasses are calculated 
         ! differently
         IF ( is_tree(ivm) ) THEN

            ! Calculate the self-thinning density.  This is the target density
            ! for which small trees are killed until the selfthinning relationship
            ! is satisfied. There are unit conversions here because the forestry 
            ! routines prefer to work with trees per hectare. Convert everything to
            ! individuals per m**2 when possible. Self-thinning relationship
            ! is fitted for the aboveground biomass/diameter so wood_to_qmdia 
            ! should be used here (wood_to_dia_eff accounts for both the above
            ! and belowground biomass).
            Dg = wood_to_qmdia(circ_class_biomass(ipts,ivm,:,:,icarbon), &
                 circ_class_n(ipts,ivm,:),ivm,pipe_tune2(ipts,ivm))

            Dg_up_half = wood_to_qmdia_up_half(circ_class_biomass(ipts,ivm,:,:,icarbon), &
                 circ_class_n(ipts,ivm,:),ivm,pipe_tune2(ipts,ivm))            

            ! After a coppice at the end of the year, it's possible that all
            ! the aboveground biomass is killed, but we still have individuals
            ! because the roots are alive. This means we don't have
            ! a diameter, though, and therefore we will divide by zero
            ! in the next block of code. So we skip it, which is okay because
            ! there will be no self-thinning in this case anyways.
            IF(Dg .LT. min_stomate) CYCLE

            ! Focus on unmanaged forests. Managed forests are mostly dealt 
            ! with in sapiens_forestry.f90 ecept for overstocking (see below).
            IF( forest_managed(ipts,ivm) .EQ. ifm_none) THEN

               ! Although the underlying ecology might be very different
               ! i.e. thinning (managed) vs gap dynamics (unmanaged) both
               ! processes can be described by an RDI approach but with
               ! different parameters. One of the strengths of the RDI 
               ! approach is that it can reduce the number of indiviudals
               ! rather quickly.
               ! It is a bit strange to start the model with a high 
               ! number of indiviudal (n_init, based on qmd_init 
               ! calculated as a parameter in stomate.f90) and then kill
               ! many of them in a relatively short time period. The high 
               ! initial number is needed because we plant small trees
               ! but require a reasonable total LAI to get enough GPP
               ! to actually grow. This transition phase is also needed
               ! to get OK tree ring width at the start of a simulation.
               ifm=forest_managed(ipts,ivm)
               ! The self-thinning relationships are for 
               ! diameters expressed in cm. rdi itself is unitless
               rdi(ipts,ivm) = (SUM(circ_class_n(ipts,ivm,:))*m2_to_ha)/ &
                       Nmax(Dg*m_to_cm, alpha_self_thinning(ipts,ivm), ivm, ifm)

               ! Calculate the rdi boundaries for unmanaged forests
               ! Polynomial approach
                CALL calculate_rdi_boundaries(Dg_up_half, &
                                             ivm, ifm, &
                                             rdi_target_upper(ipts,ivm), &
                                             rdi_target_lower(ipts,ivm))
               
               ! We initiate some gap-dynamics and self-thinning if rdi 
               ! of our stand is outside rdi range thought to be 
               ! acceptable for unmanaged forests
               IF(rdi(ipts,ivm) .GT. rdi_target_upper(ipts,ivm)) THEN
                  
                  ! Reduce the stand density to the lower rdi
                  st_dens(ipts,ivm) = &
                          Nmax(Dg*m_to_cm, alpha_self_thinning(ipts,ivm), ivm, forest_managed(ipts,ivm)) * & 
                       ha_to_m2 * rdi_target_lower(ipts,ivm)
               ELSE
                    
                  ! Don't apply gap dynamics and/or self thinning
                  st_dens(ipts,ivm) = SUM(circ_class_n(ipts,ivm,:))
                    
               ENDIF
 
            ELSE

               ! For all management strategies except unmanaged forests
               st_dens(ipts,ivm) = SUM(circ_class_n(ipts,ivm,:))
               
            END IF
               
            ! Debug
            IF(printlev_loc>=4 .AND. ivm == test_pft .AND. &
                 ipts == test_grid)THEN
               WRITE(numout,*) 'Does self thinning happen?'
               WRITE(numout,*) 'ipts,ivm: ',ipts,ivm
               WRITE(numout,*) 'circ_class_n(ipts,ivm,:),st_dens(ipts,ivm): ',&
                    SUM(circ_class_n(ipts,ivm,:)),st_dens(ipts,ivm)
            ENDIF
            !-

            ! Actual density exceeds the expected density so trees will be removed
            IF(SUM(circ_class_n(ipts,ivm,:)) .GT. st_dens(ipts,ivm))THEN

               IF(printlev_loc>=4 .AND. ivm == test_pft)THEN
                  WRITE(numout,*) 'Self thinning is taking place.'
                  WRITE(numout,*) 'ipts,ivm: ',ipts,ivm
                  WRITE(numout,*) 'circ_class_n(ipts,ivm,:),st_dens(ipts,ivm): ',&
                       SUM(circ_class_n(ipts,ivm,:)),st_dens(ipts,ivm)
               ENDIF

               ! Here is where we actually do the self thinning.  The number of
               ! trees that we want to kill is the differnce between the
               ! self thinning density and the true density.  Originally, we
               ! killed only small trees until the densities were equal, the
               ! reason being that self-thinning represents competition
               ! driven mortality, and small trees are not able to compete for
               ! resources (e.g. light) as well as larger trees. However, this
               ! seemed to cause strange results for both temperate and tropical
               ! PFTs, so now we kill trees based on a user-defined distribution.

               ! The self-thinning relationship is a theoretical maximum for a 
               ! relative large region (because it is based on inventory data)
               ! We can use this observed relationship for managed forests because
               ! the RDI will keep us below this maximum. It has been observed that
               ! unmanaged forests are also below the maximum (Luyssaert et al 2011).
               ! This is accounted for here. Being below the self-thinning, especially
               ! at early stages when canopy closure has not been reached yet, should
               ! help to make unmanaged forest grow faster (reaching more realistic
               ! heights).
               ! This is what we want to kill.

        
                ! This is the diameter that we will use
                ! to update the distribution of dead trees among
                ! the circumference classes. If for now st_dist =
                ! circ_class_dist this should change in the future. This is the
                ! reason why keep to distinct variable for stand ciconference class
                ! distribution (circ_class_dist) and dead tree circonference
                ! class distribution (st_dist).
               difference=SUM(circ_class_n(ipts,ivm,:)) - st_dens(ipts,ivm)
               st_dist(:) = weibull_class_dist(Dg*m_to_cm, ivm, &
                 forest_managed(ipts,ivm))
               target_st_kill(:)=st_dist(:)*difference
            
               ! However, we may not have a sufficient number of trees in each class,
               ! so we might have to adjust things a bit.
               true_st_kill(:)=target_st_kill(:)
               istep=1

               DO

                  lconverged=.FALSE.

                  ! Find out how many indivuals we cannot kill with the current
                  ! distribution.
                  excess_ind=zero
                  DO icir=1,ncirc
                     IF(circ_class_n(ipts,ivm,icir) .LT. true_st_kill(icir))THEN
                        excess_ind=excess_ind+true_st_kill(icir)-circ_class_n(ipts,ivm,icir)
                        true_st_kill(icir)=circ_class_n(ipts,ivm,icir)
                     ENDIF
                  ENDDO

                  IF(excess_ind .GT. zero)THEN
                     ! We have something to redistribute.
                     ! Put all of the excess in the first class which is not already completely
                     ! killed.  If we add too much, we'll catch it in the next loop.
                     DO icir=1,ncirc
                        IF(circ_class_n(ipts,ivm,icir) .GT. true_st_kill(icir))THEN
                           true_st_kill(icir)=true_st_kill(icir)+excess_ind
                           excess_ind=zero
                        ENDIF
                     ENDDO
                  ELSE
                     ! Done.
                     lconverged=.TRUE.
                  ENDIF

                  IF(lconverged)EXIT

                  istep=istep+1
                  IF(istep .GT. 100)THEN

                     ! This is an arbitrary exit here, just to make sure
                     ! we don't get stuck in an infinite loop.  We should
                     ! exit before we ever get here.  Otherwise, you can increase
                     ! the number in the if statement and try again.
                     WRITE (numout,*) 'ipts, ivm : ',ipts,ivm
                     CALL ipslerr_p (3,'stomate_mark_kill', &
                          'Was not able to distribute self-thinning individuals!', &
                          'See the comments at this part of the code for a solution','')

                  END IF
               ENDDO

               ! It is now known how many trees in each diameter class need to be killed
               circ_class_kill(ipts,ivm,:,ifm_none,icut_thin)=true_st_kill(:)

               ! Debug
               IF(printlev_loc>=4 .AND. ivm == test_pft)THEN
                  WRITE(numout,*) 'After self thinning.'
                  WRITE(numout,*) 'circ_class_kill: ',&
                       circ_class_kill(ipts,ivm,:,ifm_none,icut_thin)
               ENDIF
               !-

            ENDIF ! SUM(circ_class_n(ipts,ivm,:)) .GT. st_dens(ipts,ivm)

            !+++CEHCK+++
            ! now let's look at environmental mortality.  The calculation of this
            ! rate depends on if it is fixed for the whole simulation or if it's
            ! allowed to respond to things like NPP. The more the disturbances and
            ! mortality are developed in ORCHIDEE the less clear the difference between
            ! a non-stand replacing disturbance, self-thinning and background mortality
            ! becomes. This environmental background moratlity is left in the code as
            ! a sort of safety valve. It ensures that the number of individuals 
            ! changes a little bit every time step which prevents ORCHIDEE from being
            ! stuck. For the moment the environmental mortality is so small that
            ! it is too little to result in realistic results with (self)-thinning.
            ! As environmental mortality depends on a statistical description, i.e.,
            ! residence time, this approach should (one day) be replaced by a more
            ! mechanistic approach that simulates different casuses of mortality:
            ! drought, storms, pests, fire, competition, ...
            ! We tried running ORCHIDEE r7529 with this mortality. Running without
            ! environmenatl mortality resulted in pixels where the growth increased
            ! as well as pixel where the growth decreased. For most PFTs changes in
            ! growth were reslatively small. At this time it was decided to keep this
            ! code as a safety valve.

            !! 3.1 Use growth efficiency or constant mortality?
            IF ( .NOT. ok_constant_mortality  ) THEN

               ! Was the PFT alive last year?
               IF ( ( lm_lastyearmax(ipts,ivm) .GT. min_stomate ) ) THEN

                  !! 3.1.1 Estimate net biomass increment
                  ! Estimate net biomass increment by subtracting turnover from
                  ! last year's NPP.
                  ! Note that npp_longterm is the longterm growth efficiency
                  ! (NPP/LAI) to be fair to deciduous trees
                  delta_biomass = MAX( npp_longterm(ipts,ivm) - &
                       ( turnover_longterm(ipts,ivm,ileaf,icarbon) + &
                       turnover_longterm(ipts,ivm,iroot,icarbon) + &
                       turnover_longterm(ipts,ivm,ifruit,icarbon) + & 
                       turnover_longterm(ipts,ivm,isapabove,icarbon) + &
                       turnover_longterm(ipts,ivm,isapbelow,icarbon) ), zero)

                  !! 3.1.2 Calculate growth efficiency
                  !  Calculate growth efficiency by dividing net biomass increment
                  !  by last year's maximum LAI. (corresponding actually to the
                  !  maximum LAI of the previous year). By using the function 
                  !  biomass_to_lai dynamic sla is accounted for.
                  vigour = delta_biomass / biomass_to_lai(lm_lastyearmax(ipts,ivm),ivm)

               ELSE
                 
                  ! PFT does not exist or is dead
                  vigour = zero

               END IF ! ( lm_lastyearmax(ipts,ivm) .GT. min_stomate)

               !! 3.1.3 Calculate growth efficiency mortality rate
               mortality_greff = mortality_max(ivm) / &
                    ( un + ref_mortality(ivm) * vigour )

               ! Scale mortality rate (fraction, 0-1) by timesteps per year
               mortality = MAX(mortality_min(ivm), mortality_greff) * dt/&
                    one_year 

            ELSE ! .NOT. ok_constant_mortality

               !! 3.2 Use constant mortality accounting for the residence time
               !! of each tree PFT

               !  Scale mortality rate (fraction, 0-1) by timesteps per year
               !  NOTE: the meaning of residence_time is very different between
               !  the DOFOCO branch and the trunk. In the trunk biomass has no age
               !  thus the residence time accounts for all forest dynamics including
               !  self-thinning, pests, diseases and windthrow. In the DOFOCO branch
               !  biomass does have an age and self-thinning is explicitly accounted
               !  for, hence, the residence time should be much higher as it only
               !  accounts for pest, diseases and windthrow. Even the latter is not
               !  exact because as long as those disturbances are small scale they
               !  are probably accounted for in the parametrization of self-thinning.
               mortality = dt/(residence_time(ivm)*one_year)

            ENDIF ! .NOT.  ok_constant_mortality
            !+++++++++++

            ! Debug
            IF (printlev_loc>=4 .AND. ivm .EQ. test_pft) THEN
               WRITE(numout,*) 'Mortality rate (0-1), ', mortality
            ENDIF
            !-

            ! Now the mortality is just a percentage of the total biomass
            ! on the site.  Since the main variable we output here is
            ! circ_class_kill, we need to convert that biomass into the
            ! number of individuals that are killed, and from which
            ! circumference classes.  Mortality in excess of the self-thinning
            ! mortality will preferenctially kill large trees, since they are
            ! assumed to be less resistant to climatic change.

            ! Aggregate over the different biomass components
            DO ipar = 1,nparts

               ! Aggregate over the different circumference classes
               DO icir = 1,ncirc

                  ! Mortality from self-thinning
                  ! Note that the unit of st_mortaility (g C m-2 dt-1) differs
                  ! from that of mortality (rate dt-1)
                  st_mortality(ipts,ivm) = st_mortality(ipts,ivm) + &
                       circ_class_kill(ipts,ivm,icir,ifm_none,icut_thin) * &
                       circ_class_biomass(ipts,ivm,icir,ipar,icarbon)

               ENDDO

               ! Ageing and environmental mortality
               ! Note that the unit of d_mortaility (g C m-2 dt-1) differs from
               ! that of mortality (rate dt-1)
               d_mortality(ipts,ivm) =  d_mortality(ipts,ivm) + &
                    mortality * SUM(circ_class_biomass(ipts,ivm,:,ipar,icarbon)*&
                    circ_class_n(ipts,ivm,:))
       
            ENDDO ! ipar = 1,nparts

            IF(printlev_loc>=4 .AND. test_pft .EQ. ivm .AND. &
                 test_grid == ipts)THEN
               WRITE(numout,*) 'Does mortality happen?'
               WRITE(numout,*) 'ipts,ivm: ',ipts,ivm
               WRITE(numout,*) 'd_mortality(ipts,ivm): ',d_mortality(ipts,ivm)
               WRITE(numout,*) 'st_mortality(ipts,ivm): ',st_mortality(ipts,ivm)
            ENDIF


            !+++CHECK+++
            ! There are different ways to combine these two sources of 
            ! mortality. First one could argue that self-thinning
            ! includes the background mortality. If background mortality
            ! exceeds self-thinning this is due to the simplified way
            ! background mortality is calculated. 
            IF (st_mortality(ipts,ivm) .GT. min_stomate) THEN
               
               ! There is self-thinning mortality, no need to account
               ! for environmental mortality
               bm_difference = zero

            ELSE

               ! There is no self_thinning mortality, use the 
               ! environmental mortality
               bm_difference = d_mortality(ipts,ivm)

            ENDIF

            ! The second way of looking at this is that the calculation of
            ! background mortality is reliable and should therefore always 
            ! be taken into account. In this reasoning self-thinning is
            ! considered part of the background mortality. Account for both
            ! but avoid double counting. Only one approach can be used.
            !bm_difference=d_mortality(ipts,ivm) - st_mortality(ipts,ivm)
            !+++++++++++
      
            !+++CHECK +++
            ! The question still stands of which trees we kill with the
            ! environmental mortality.  If we only kill the biggest, we
            ! lose a circ class within the first year sometimes.
            IF(bm_difference .GT. zero )THEN

               ! We know the total biomass that we want to kill.  Let's partition
               ! this biomass among all our circumference classes based on an
               ! exponential distribution which takes more biomass from
               ! the largest class than the smallest.  There is a chance that we
               ! will not have enough biomass to kill in some classes, in which
               ! case we need to rearrange it a bit.

               ! Debug
               IF(printlev_loc>=4 .AND. ivm == test_pft .AND. &
                    ipts == test_grid)THEN
                  WRITE(numout,*) 'Mortality is taking place.'
                  WRITE(numout,*) 'ipts,ivm: ',ipts,ivm
                  WRITE(numout,*) 'bm_difference: ',bm_difference
                  WRITE(numout,*) 'circ_class_kill init: ', &
                       circ_class_kill(ipts,ivm,:,ifm_none,icut_thin)
                  WRITE(numout,*) 'circ_class_n init: ', &
                        circ_class_n(ipts,ivm,:)
               ENDIF
               !-

               CALL distribute_mortality_biomass( bm_difference, &
                    death_distribution_factor(ivm), circ_class_n(ipts,ivm,:), &
                    circ_class_biomass(ipts,ivm,:,:,icarbon), &
                    circ_class_kill(ipts,ivm,:,ifm_none,icut_thin) )

               ! Debug
               IF(printlev_loc>=4 .AND. ivm == test_pft .AND. &
                    ipts == test_grid)THEN
                  WRITE(numout,*) 'after circ_class_kill: ', &
                       circ_class_kill(ipts,ivm,:,ifm_none,icut_thin)
                  WRITE(numout,*) 'after circ_class_n : ', &
                       circ_class_n(ipts,ivm,:)
               ENDIF
               !-

            ENDIF ! bm_difference .GT. zero
            
         ELSE ! IF(is_tree(ivm))

            ! Calculate how many individuals should be killed due
            ! to low reserve + labile carbon (reserve+labile<target)
            IF(natural(ivm)) THEN
               ! Conditions to decrease density include:
               ! (1) The flag ok_dyn_grass_den for the calculation of dynamic 
               ! grass density is switched on; 
               ! (2) The veget_max for each grassland is greater than 10-8;
               ! (3) The plant_status is at ipresenescence, isenescent or idormant; (4)Current reserve and labile 
               ! carbon are lower than the target of reserve and labile carbon;
               ! (5)The reserve and labile carbon is greater than 10-8; 
               ! (6)Grass density is greater than the minimal value 0.05.
               IF ( (ok_dyn_grass_den) .AND. &
                    (veget_max(ipts,ivm) .GT. min_stomate) .AND. &
                    ((plant_status(ipts,ivm) .EQ. ipresenescence) .OR. &
                    (plant_status(ipts,ivm) .EQ. isenescent) .OR. &
                    (plant_status(ipts,ivm) .EQ. idormant)) .AND. &
                    (circ_class_n(ipts,ivm,1)*(circ_class_biomass(ipts,ivm,1,icarbres,icarbon)+ &
                    circ_class_biomass(ipts,ivm,1,ilabile,icarbon)) .LT.&
                    (tot_res_target(ipts,ivm,icarbon) + &
                    tot_lab_target(ipts,ivm,icarbon))) .AND. &
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
                     WRITE(numout,*) 'total biomass', &
                          circ_class_biomass(ipts,ivm,1,:,icarbon)*circ_class_n(ipts,ivm,1)
                  END IF

                  ! To decrease grassland density,the calculation should follow the two principles:
                  ! 1. Mass balance should be preserved before and after the density change, so:
                  ! eq.1: N0*(Cother+Trl0) = N1*(Cother + Trl1) +scale_carbon*(N0-N1)*Cother
                  ! 2. Cother will remain constant, and individual reserve and labile
                  ! carbon will reach the target:((N0+scale_dec*N1)/N1)*(Ttarget_area)/N0
                  ! Combining 1 & 2, we can calculate new density as:
                  ! N1=(N0*((1-scale_carbon)*Cother+Trl0)-Ttarget_area) / 
                  ! ((1-scale_carbon)*Cother+scale_dec*Ttarget_area/N0)

                  ! Calculate the sum of reserve and labile carbon, Trl0 (gC ind-1) 
                  Trl0=circ_class_biomass(ipts,ivm,1,icarbres,icarbon)+ &
                       circ_class_biomass(ipts,ivm,1,ilabile,icarbon)

                  ! Calculate the sum of carbon from all compartments except for
                  ! reserve and labile pool, Cother (gC ind-1)
                  Cother = SUM(circ_class_biomass(ipts,ivm,1,:,icarbon))-Trl0

                  ! Calculate the sum of areal target for reserve and labile
                  ! carbon, Ttarget_area (gC m-2)
                  Ttarget_area = tot_res_target(ipts,ivm,icarbon)+ &
                       tot_lab_target(ipts,ivm,icarbon)

                  ! Calculate new density
                  new_circ_class_n(ipts,ivm,1)=(circ_class_n(ipts,ivm,1)*((1-scale_carbon)*Cother+Trl0)-Ttarget_area)/&
                       ((1-scale_carbon)*Cother+scale_dec*Ttarget_area/circ_class_n(ipts,ivm,1))

                  ! The minimum density is 0.05 to avoid numerical error/crash,
                  ! constrain to 0.05 if new density<0.05
                  new_circ_class_n(ipts,ivm,1) = MAX(new_circ_class_n(ipts,ivm,1), 0.05)
                  old_circ_class_n(ipts,ivm,1) = circ_class_n(ipts,ivm,1)
   
                  !! Calculate the number of individuals to be killed,
                  !circ_class_kill, ind m-2
                  circ_class_kill(ipts,ivm,1,ifm_none,icut_thin) = &
                       old_circ_class_n(ipts,ivm,1)-new_circ_class_n(ipts,ivm,1)

                  !! To prevent precision error, if circ_class_kill is lower
                  !than 10000*EPSILON(un), we keep it at zero
                  IF (circ_class_kill(ipts,ivm,1,ifm_none,icut_thin) .LT. &
                       10000*EPSILON(un)) THEN
                     circ_class_kill(ipts,ivm,1,ifm_none,icut_thin) = zero
                  END IF !

                  ! Debug
                  IF (printlev_loc >= 4) THEN
                      WRITE(numout,*) 'After density reduction,ipts,ivm', ipts, ivm
                      WRITE(numout,*) 'Density, after',circ_class_n(ipts,ivm,1)
                      WRITE(numout,*) 'current reserve',circ_class_biomass(ipts,ivm,1,icarbres,icarbon)
                      WRITE(numout,*) 'total reserve target',tot_res_target(ipts,ivm,icarbon)
                      WRITE(numout,*) 'current labile',circ_class_biomass(ipts,ivm,1,ilabile,icarbon)
                      WRITE(numout,*) 'total labile target',tot_lab_target(ipts,ivm,icarbon)
                      WRITE(numout,*) 'total biomass', &
                          circ_class_biomass(ipts,ivm,1,:,icarbon)*circ_class_n(ipts,ivm,1)
                      WRITE(numout,*) 'killed density',circ_class_kill(ipts,ivm,1,ifm_none,icut_thin)
                  END IF !Debug
               ELSE
                  circ_class_kill(ipts,ivm,:,ifm_none,icut_thin)=zero
               END IF !Condition to reduce density
            END IF ! natural

            ! This is the original code for the mortality.  Notice that it doesn't actually do
            ! anything, it is left as documentation for when the dgvm is restored.
!!$             IF ( ok_dgvm .OR. .NOT. ok_constant_mortality) THEN
!!$
!!$                ! +++CHECK+++
!!$                ! For grasses, if last year's NPP is very small (less than 10 gCm^{-2}year{-1})
!!$                ! the grasses completely die. Grasses and crops are deciduous so in the first
!!$                ! year they have no leaves. Consequently npp_longterm is less then
!!$                ! ::npp_longterm_init. At the beginning of the second year when rue_longterm is
!!$                ! no longer 1, the PFT is killed immediatly. The IF statement needs to be refined
!!$                ! to properly allow for this type of PFT killing.
!!$                ! The more fundamental and conceptual question is whether we want to kill grasses?
!!$                ! Unless we simulate desertification, it is not clear why we would kill the grass
!!$                ! most likely to replace it by ... grass. Because grass has hardly any biomass it
!!$                ! seem unlikely that killing the grass will affect the C-fluxes.
!!$                ! For both reasons i.e. the technical and conceptual, it was decided NOT to kill
!!$                ! the grass PFTs
!!$                IF ( PFTpresent(ipts,ivm) .AND. ( npp_longterm(ipts,ivm) .LE. npp_longterm_init ) ) THEN
!!$
!!$                   ! ORIGINAL code - see above to justify the change
!!$                   !!$mortality = un
!!$                   mortality = zero
!!$
!!$                END IF
!!$                ! +++++++++++
!!$
!!$                ! Update biomass and litter pools
!!$                DO ielem = 1,nelements
!!$
!!$                   DO ipart = 1, nparts
!!$
!!$                      IF ( PFTpresent(ipts,ivm) ) THEN
!!$
!!$                         d_mortality(ipts,ivm) =  mortality * biomass(ipts,ivm,ipart,ielem)
!!$                         bm_to_litter(ipts,ivm,ipart,ielem) = bm_to_litter(ipts,ivm,ipart,ielem) + &
!!$                              d_mortality(ipts,ivm)
!!$                         biomass(ipts,ivm,ipart,ielem) = biomass(ipts,ivm,ipart,ielem) - &
!!$                              d_mortality(ipts,ivm)
!!$
!!$                      END IF ! PFTpresent(ipts,ivm)
!!$
!!$                   ENDDO ! nparts
!!$
!!$                END DO ! nelements
!!$
!!$             ENDIF ! .NOT.ok_dgvm .AND. .NOT.ok_constant_mortality

         ENDIF ! end check for trees
         
      ENDDO ! loop over grid points

   ENDDO ! loop over PFTs

   ! Write output as stored in ORCHIDEE -> NCIRC
   ! Useful to understand what happens in the model
   CALL xios_orchidee_send_field("CCMORTALITY",SUM(SUM(circ_class_kill(:,:,:,:,:),5),4))
   ! Use same definition for aboveground and belowground components as for the recruits
   tmp(:,:,:) = (circ_class_biomass(:,:,:,isapabove,icarbon) + &
        circ_class_biomass(:,:,:,iheartabove,icarbon)) * &
        SUM(SUM(circ_class_kill(:,:,:,:,:),5),4)
   CALL xios_orchidee_send_field("CCMORTALITY_M_AB_c",tmp(:,:,:))
   tmp(:,:,:) = (circ_class_biomass(:,:,:,isapbelow,icarbon) + &
        circ_class_biomass(:,:,:,iheartbelow,icarbon)) * &
        SUM(SUM(circ_class_kill(:,:,:,:,:),5),4)
   CALL xios_orchidee_send_field("CCMORTALITY_M_BE_c",tmp(:,:,:))


   ! Write output file with fixed diamter classes -> NOUTDIACLASS
   ! Useful for model intercomparison.
    ! convert biomass into diameters
    temp1(:,:,:) = zero
    temp2(:,:,:) = zero
    DO ipts = 1,npts
       DO ivm = 1,nvm
          IF(is_tree(ivm))THEN
             circ_dia(:) = wood_to_dia(circ_class_biomass(ipts,ivm,:,:,icarbon),&
                  ivm,pipe_tune2(ipts,ivm))
             istart = 2
             DO icir = 1,ncirc
                DO iout = istart,noutdiaclass+1
                   IF (circ_dia(icir).GE.out_dia_class(iout-1) .AND. & 
                        circ_dia(icir).LT.out_dia_class(iout)) THEN
                      ! Number of dead trees
                      temp1(ipts,ivm,iout-1) = temp1(ipts,ivm,iout-1) + &
                           SUM(SUM(circ_class_kill(ipts,ivm,icir,:,:),2),1)
                      ! Aboveground woody biomass of dying trees
                      temp2(ipts,ivm,iout-1) = temp2(ipts,ivm,iout-1) + &
                           (circ_class_biomass(ipts,ivm,icir,isapabove,icarbon) + &
                           circ_class_biomass(ipts,ivm,icir,iheartabove,icarbon)) * &
                           SUM(SUM(circ_class_kill(ipts,ivm,icir,:,:),2),1)
                      istart = iout
                   END IF
                END DO
             END DO ! icirc
          END IF ! is_tree
       END DO ! ivm
    END DO ! ipts
        
    ! Write to xios
    CALL xios_orchidee_send_field("stemmort_dben",temp1(:,:,:)*m2_to_ha)
    CALL xios_orchidee_send_field("cmort_size_dben",temp2(:,:,:)/1e3)


 !! 4. Check numerical consistency of this routine

    IF (err_act.GT.1) THEN
       
       ! 4.2 Check surface area
       CALL check_vegetation_area("stomate_mark_to_kill", npts, veget_max_begin, &
            veget_max,'pft')
       
       ! 4.3 Mass balance closure
       ! 4.3.1 Calculate final biomass
       pool_end = zero
       DO ipar = 1,nparts
          DO iele = 1,nelements
             DO icir = 1,ncirc
                pool_end(:,:,iele) = pool_end(:,:,iele) + &
                     circ_class_biomass(:,:,icir,ipar,iele) * &
                     circ_class_n(:,:,icir) * veget_max(:,:)
             ENDDO
          ENDDO
       ENDDO

       !! 4.2 Calculate components of the mass balance
       ! 4.3.2 Calculate mass balance
       ! Common processes
       DO iele=1,nelements
          check_intern(:,:,ipoolchange,iele) = -un * (pool_end(:,:,iele) - &
               pool_start(:,:,iele)) 
       ENDDO
   
       closure_intern = zero
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
       CALL check_mass_balance("stomate_mark_to_kill", closure_intern, npts, &
            pool_end, pool_start, veget_max, 'pft')

    ENDIF ! err_act.GT.1 
   
    IF (printlev>=4) WRITE(numout,*) 'Leaving mark_to_kill'

 END SUBROUTINE mark_to_kill

!! ================================================================================================================================
!! SUBROUTINE  : distribute_mortality_biomass
!!
!>\BRIEF       Distributes biomass that is going to be killed by natural
!!             causes (not self thinning) over circ classes.
!!
!! DESCRIPTION : Mortality is going to kill a certain amount of biomass
!!               in forests every day.  Since we have circumference classes
!!               now for our forests, we need to determine which classes
!!               of trees will suffer from this environmental mortality.
!!               Right now we are taking an exponential distribution.
!!               Notice that this is NOT the same as redistributing biomass
!!               after one of the circ classes becomes empty.
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::circ_class_kill
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE distribute_mortality_biomass ( bm_difference, ddf_temp, circ_class_n_temp, &
       circ_class_biomass_temp, circ_class_kill_temp )

 !! 0. Variable and parameter description

    !! 0.1 Input variables
    REAL(r_std),INTENT(in)                       :: bm_difference           !! the biomass to distribute
    REAL(r_std),INTENT(in)                       :: ddf_temp                !! the death_distribution_factor for this pft
    REAL(r_std),DIMENSION(:),INTENT(in)          :: circ_class_n_temp       !! circ_class_n for this point/PFT
    REAL(r_std),DIMENSION(:,:),INTENT(in)        :: circ_class_biomass_temp !! circ_class_biomass for this point/PFT
 
   !! 0.2 Output variables

    !! 0.3 Modified variables
    REAL(r_std), DIMENSION(:), INTENT(out)      :: circ_class_kill_temp     !! Number of trees within a circ that needs
                                                                            !! to be killed @tex $(ind m^{-2})$ @endtex

    !! 0.4 Local variables
    REAL(r_std), DIMENSION(ncirc)               :: biomass_desired          !! The biomass that dies naturally
    REAL(r_std), DIMENSION(ncirc)               :: death_distribution       !! The fraction of biomass taken from
                                                                            !! each circ class for mortality.
    LOGICAL                                     :: ldone                    !! Flag to exit a loop
    REAL(r_std)                                 :: scale_factor             !!
    REAL(r_std)                                 :: sum_total                !!
    REAL(r_std)                                 :: leftover_bm              !! excess biomass that we need to kill
    REAL(r_std)                                 :: living_biomass           !! summed biomass
    REAL(r_std)                                 :: living_trees             !! summed biomass
    INTEGER                                     :: icir

!_ ================================================================================================================================

    IF(ncirc == 1)THEN

       ! This is the easy case.  All of our biomass will be taken from the only
       ! circumference class that we have.
       circ_class_kill_temp(1)=circ_class_kill_temp(1)+&
            bm_difference/SUM(circ_class_biomass_temp(1,:))

       RETURN
    ENDIF

    ! Here we assume an exponential distribution, arranged so that
    ! ddf_temp times more biomass is taken from the largest circ class
    ! compared to the smallest.
    biomass_desired(:)=zero
    death_distribution(:)=un
    scale_factor=ddf_temp**(un/REAL(ncirc-1))
    DO icir=2,ncirc
       death_distribution(icir)=death_distribution(icir-1)*scale_factor
    ENDDO
    ! Normalize it
    sum_total=SUM(death_distribution(:))
    death_distribution(:)=death_distribution(:)/sum_total

    ! Ideally, how much is killed from each class?  Be careful to include
    ! what was killed in self-thinning here!  If we don't, we may try to kill
    ! more biomass than is available in the loop below.
    DO icir=1,ncirc
       biomass_desired(icir)=death_distribution(icir)*bm_difference+&
            circ_class_kill_temp(icir)*SUM(circ_class_biomass_temp(icir,:))
    ENDDO

    ! Right now, we know how much biomass we want to kill in each
    ! class (biomass_desired).  What we will do now is to loop through
    ! all the circ_classes and see if we have this much biomass
    ! in each class still alive.  The total amount of vegetation still
    ! alive is circ_class_n_temp(icir)-circ_class_kill_temp(icir).  If
    ! this number of individuals cannot give us the total biomass that
    ! we need, we keep track of a residual quantity, leftover_bm, which
    ! we try to take from other circ_classes.  It's possible that we
    ! will have to loop several times, which makes things more complicated.

    ldone=.FALSE.
    leftover_bm=zero
    DO
       DO icir=ncirc,1,-1

          living_trees=circ_class_n_temp(icir)-circ_class_kill_temp(icir)
          living_biomass=&
               SUM(circ_class_biomass_temp(icir,:))*living_trees
          biomass_desired(icir)=biomass_desired(icir)+leftover_bm

          IF(living_biomass .LE. biomass_desired(icir))THEN

             ! We can't get everything from this class, so we kill whatever is left.
             leftover_bm=biomass_desired(icir)-living_biomass
             biomass_desired(icir)=zero
             circ_class_kill_temp(icir)=circ_class_kill_temp(icir)+&
                  living_trees

          ELSE

             ! We have enough in this class.
             circ_class_kill_temp(icir)=circ_class_kill_temp(icir)+&
                  biomass_desired(icir)/SUM(circ_class_biomass_temp(icir,:))
             biomass_desired(icir)=zero
             leftover_bm=zero

          ENDIF ! living_biomass .LE. biomass_desired(icir)

       ENDDO ! loop over circ classes

       IF(leftover_bm .LE. min_stomate) EXIT

       ! it's possible that we don't have enough biomass left to kill what needs to be
       ! killed, so everything just dies.  I cannot think of a case where this
       ! would happen, though, since mortality should always be a percentage of the
       ! total biomass.
       ldone=.TRUE.
       DO icir=1,ncirc

          IF( circ_class_kill_temp(icir) .LT. circ_class_n_temp(icir) ) &
               ldone=.FALSE. 

       ENDDO

       IF(ldone)EXIT ! All our biomass is dead, and we still want to kill more!

    ENDDO

  END SUBROUTINE distribute_mortality_biomass

END MODULE stomate_mark_kill
