! =================================================================================================================================
! MODULE       : stomate_turnover.f90
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        This module manages the end of the growing season and calculates herbivory 
!               and turnover of leaves, fruits, fine roots.
!! %
!!\n DESCRIPTION: This subroutine calculates leaf senescence due to climatic conditions or 
!! as a function of leaf age and new LAI, and subsequent turnover of the different plant 
!! biomass compartments (sections 1 to 6), herbivory (section 7), fruit turnover for trees 
!! (section 8) and sapwood conversion (section 9). 
!!
!! RECENT CHANGE(S): None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_stomate/stomate_turnover.f90 $
!! $Date: 2025-12-03 11:32:58 +0100 (mer. 03 déc. 2025) $
!! $Revision: 9224 $
!_ ================================================================================================================================

MODULE stomate_turnover

  ! modules used:
  USE xios_orchidee
  USE ioipsl_para
  USE stomate_data
  USE constantes
  USE grid
  USE pft_parameters
  USE sapiens_agriculture, ONLY: crop_harvest
  USE function_library, ONLY : biomass_to_lai, &
       check_vegetation_area, check_mass_balance, lai_to_biomass, get_printlev
  USE time, ONLY : julian_diff

  IMPLICIT NONE

  ! private & public routines

  PRIVATE
  PUBLIC turn_over, turnover_clear, drought_mortality

  LOGICAL, SAVE        :: firstcall_turnover = .TRUE.      !! first call (true/false)
!$OMP THREADPRIVATE(firstcall_turnover)
  INTEGER(i_std), SAVE :: printlev_loc                     !! Local level of text output for current module
!$OMP THREADPRIVATE(printlev_loc)
CONTAINS


!! ================================================================================================================================
!! SUBROUTINE   : turnover_clear
!!
!>\BRIEF        Set flag ::firstcall_turnover to .TRUE., and therefore activate section 1 
!!              of subroutine turn which writes a message to the output.
!!                
!_ ================================================================================================================================

  SUBROUTINE turnover_clear
    firstcall_turnover=.TRUE.
  END SUBROUTINE turnover_clear


!! ================================================================================================================================
!! SUBROUTINE    : turn_over
!!
!>\BRIEF         Calculate turnover of leaves, roots, fruits and sapwood
!!               due to aging or climatic induced senescence. Calculate
!!               herbivory.
!!
!! DESCRIPTION : This subroutine determines the turnover of leaves and fine
!! roots (and stems for grasses) and simulates following processes:
!! 1. Mean leaf age is calculated from leaf ages of separate leaf age
!!    classes. Should actually be recalculated at the end of the routine,
!!    but it does not change too fast. The mean leaf age is calculated
!!    using the following equation:
!!    \latexonly
!!    \input{turnover_lma_update_eqn1.tex}
!!    \endlatexonly
!!    \n 
!! 2. Meteorological senescence: the detection of the end of the growing
!!    season and shedding of leaves, fruits and fine roots due to
!!    unfavourable meteorological conditions. The model distinguishes
!!    three different types of "climatic" leaf senescence, that do not 
!!    change the age structure: sensitivity to cold temperatures, to lack
!!    of water, or both. If meteorological conditions are fulfilled, a
!!    flag ::senescence is set to TRUE. Note that evergreen species do not
!!    experience climatic senescence. Climatic senescence is triggered by
!!    sensitivity to cold temperatures where the critical temperature for
!!    senescence is calculated using the following equation:
!!    \latexonly
!!    \input{turnover_temp_crit_eqn2.tex}
!!    \endlatexonly
!!    \n
!!    Climatic senescence is triggered by sensitivity to lack of water
!!    availability where the moisture availability critical level is
!!    calculated using the following equation:
!!    \latexonly
!!    \input{turnover_moist_crit_eqn3.tex}
!!    \endlatexonly
!!    \n
!!    Climatic senescence is triggered by sensitivity to temperature or to
!!    lack of water where critical temperature and moisture availability
!!    are calculated as above.\n
!!    Trees in climatic senescence lose their fine roots at the same rate
!!    as they lose their leaves. The rate of biomass loss of both fine
!!    roots and leaves is presribed through the equation:
!!    \latexonly
!!    \input{turnover_clim_senes_biomass_eqn4.tex}
!!    \endlatexonly
!!    \n
!!    with ::leaffall(j) a PFT-dependent time constant which is given in 
!!    ::stomate_constants. In grasses, leaf senescence is extended to
!!    the whole plant (all carbon pools) except to its carbohydrate reserve.    
!! 3. Senescence due to aging: the loss of leaves, fruits and  biomass due
!!    to aging At a certain age, leaves fall off, even if the climate
!!     would allow a green plant all year round. Even if the meteorological
!!    conditions are favorable for leaf maintenance, plants, and in
!!    particular, evergreen trees, have to renew their leaves simply because
!!    the old leaves become inefficient. Roots, fruits (and stems for
!!    grasses) follow leaves. The ??senescence?? rate varies with leaf
!!    age. Note that plant is not declared senescent in this case (wchich
!!    is important for allocation: if the plant loses leaves because of 
!!    their age, it can renew them). The leaf turnover rate due to aging
!!    of leaves is calculated using the following equation:
!!    \latexonly
!!    \input{turnover_age_senes_biomass_eqn5.tex}
!!    \endlatexonly
!!    \n
!!    Drop all leaves if there is a very low leaf mass during senescence.
!!    After this, the biomass of different carbon pools both for trees
!!    and grasses is set to zero and the mean leaf age is reset to zero.
!!    Finally, the leaf fraction and leaf age of the different leaf age
!!    classes is set to zero. For deciduous trees: next to leaves, also
!!    fruits and fine roots are dropped.
!!    For grasses: all aboveground carbon pools, except the carbohydrate
!!    reserves are affected:
!! 4. Update the leaf biomass, leaf age class fraction and the LAI
!!    Older leaves will fall more frequently than younger leaves and
!!    therefore the leaf age distribution needs to be recalculated after
!!    turnover. The fraction of biomass in each leaf class is updated using
!!    the following equation:
!!    \latexonly
!!    \input{turnover_update_LeafAgeDistribution_eqn6.tex}
!!    \endlatexonly
!!    \n
!! 5. Simulate herbivory activity and update leaf and fruits biomass.
!!    Herbivore activity affects the biomass of leaves and fruits as well
!!    as stalks (only for grasses). However, herbivores do not modify leaf
!!    age structure.
!! 6. Calculates fruit turnover for trees. Trees simply lose their fruits
!!    with a time constant ::longevity_fruit(j), that is set to 90 days for all
!!    PFTs in ::stomate_constants 
!! 7. Convert sapwood to heartwood for trees and update heart and softwood
!!    above and belowground biomass. Sapwood biomass is converted into
!!    heartwood biomass with a time constant tau ::longevity_sap(j) of 1 year.
!!    Note that this biomass conversion is not added to "turnover" as the
!!    biomass is not lost. For the updated heartwood, the sum of new
!!    heartwood above and new heartwood below after converting sapwood to 
!!    heartwood, is saved as ::hw_new(:). Creation of new heartwood
!!    decreases the age of the plant ??carbon?? with a factor that is
!!    determined by: old heartwood ::hw_old(:) divided by the new
!!    heartwood ::hw_new(:)
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLES: ::Biomass of leaves, fruits, fine roots and
!!  sapwood above (latter for grasses only), ::Update leaf age distribution
!!  with new leaf age class fraction 
!!
!! REFERENCE(S) : 
!! - Krinner, G., N. Viovy, N. de Noblet-Ducoudre, J. Ogee, J. Polcher, P. 
!! Friedlingstein, P. Ciais, S. Sitch and I.C. Prentice (2005), A dynamic global
!! vegetation model for studies of the coupled atmosphere-biosphere system,
!! Global Biogeochemical Cycles, 19, doi:10.1029/2003GB002199. 
!! - McNaughton, S. J., M. Oesterheld, D. A. Frank and K. J. Williams (1989), 
!! Ecosystem-level patterns of primary productivity and herbivory in
!! terrestrial habitats, Nature, 341, 142-144, 1989. 
!! - Sitch, S., C. Huntingford, N. Gedney, P. E. Levy, M. Lomas, S. L. Piao,
!! Betts, R., Ciais, P., Cox, P., Friedlingstein, P., Jones, C. D.,
!! Prentice, I. C. and F. I. Woodward : Evaluation of the terrestrial carbon  
!! cycle, future plant geography and climate-carbon cycle feedbacks using 5
!! dynamic global vegetation models (dgvms), Global Change Biology, 14(9),
!! 2015 –2039, 2008. 
!!
!! FLOWCHART    : 
!! \latexonly
!! \includegraphics[scale=0.5]{turnover_flowchart_1.png}
!! \includegraphics[scale=0.5]{turnover_flowchart_2.png}
!! \endlatexonly
!! \n
!_ ================================================================================================================================

  SUBROUTINE turn_over (npts, dt, PFTpresent, herbivores, &
       gpp_week, resp_maint_week, &
       maxvegstress_lastyear, minvegstress_lastyear, vegstress_week, &
       vegstress_month, t2m_longterm, t2m_month, t2m_week, &
       veget_max, gdd_from_growthinit, leaf_age, leaf_frac, age, &
       turnover, plant_status, turnover_time, &
       circ_class_biomass, circ_class_n, &
       when_growthinit, longevity_eff_leaf, longevity_eff_sap, longevity_eff_root, &
       harvest_pool, harvest_type, harvest_cut, harvest_area, &
       wstress_month, leaf_age_crit, doy_end_gs)

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables 

    INTEGER(i_std), INTENT(in)                       :: npts                 !! Domain size - number of grid cells 
                                                                             !! (unitless) 
    REAL(r_std), INTENT(in)                          :: dt                   !! time step (dt_days)
    LOGICAL, DIMENSION(:,:), INTENT(in)              :: PFTpresent           !! PFT exists (true/false)
    REAL(r_std), DIMENSION(:,:), INTENT(in)          :: herbivores           !! time constant of probability of a leaf to 
                                                                             !! be eaten by a herbivore (days) 
    REAL(r_std), DIMENSION(:,:), INTENT(in)          :: maxvegstress_lastyear !! last year's maximum moisture availability 
                                                                             !! (0-1, unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(in)          :: minvegstress_lastyear !! last year's minimum moisture availability 
                                                                             !! (0-1, unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(in)          :: vegstress_week        !! "weekly" moisture availability 
                                                                             !! (0-1, unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(in)          :: vegstress_month       !! "monthly" moisture availability 
                                                                             !! (0-1, unitless)
    REAL(r_std), DIMENSION(:), INTENT(in)            :: t2m_longterm         !! "longterm" 2-meter temperatures (K)
    REAL(r_std), DIMENSION(:), INTENT(in)            :: t2m_month            !! "monthly" 2-meter temperatures (K)
    REAL(r_std), DIMENSION(:), INTENT(in)            :: t2m_week             !! "weekly" 2 meter temperatures (K)
    REAL(r_std), DIMENSION(:,:), INTENT(in)          :: veget_max            !! "maximal" coverage fraction of a PFT (LAI 
                                                                             !! -> infinity) on ground (unitless) 
    REAL(r_std), DIMENSION(:,:), INTENT(in)          :: gpp_week             !! PFT gross primary productivity 
    REAL(r_std), DIMENSION(:,:), INTENT(in)          :: resp_maint_week      !! Weekly maintenance respiration
    REAL(r_std), DIMENSION(:,:), INTENT(in)          :: gdd_from_growthinit  !! gdd senescence for crop
    REAL(r_std), DIMENSION(:,:), INTENT(in)          :: longevity_eff_root   !! Effective root turnover time that accounts
                                                                             !! waterstress (days)
    REAL(r_std), DIMENSION(:,:), INTENT(in)          :: longevity_eff_sap    !! Effective sapwood turnover time that accounts
                                                                             !! waterstress (days)
    REAL(r_std), DIMENSION(:,:), INTENT(in)          :: longevity_eff_leaf   !! Effective leaf turnover time that accounts
                                                                             !! waterstress (days)
    REAL(r_std), DIMENSION(:,:), INTENT(in)          :: wstress_month        !! Water stress factor, based on hum_rel_daily
                                                                             !! (unitless, 0-1)
    REAL(r_std), DIMENSION(:,:), INTENT(in)          :: leaf_age_crit        !! critical leaf age (days)

    !! 0.2 Output variables
    REAL(r_std),DIMENSION(:,:),INTENT(out)           :: doy_end_gs           !! growing season end day of year (DOY) for 
                                                                             !! deciduous PFTs.
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(out)     :: turnover             !! Turnover @tex ($gC m^{-2}$) @endtex

    !! 0.3 Modified variables

    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)     :: leaf_age             !! age of the leaves (days)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)     :: leaf_frac            !! fraction of leaves in leaf age class 
                                                                             !! (0-1, unitless)
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout) :: harvest_pool         !! The wood and biomass that have been
                                                                             !! havested by humans @tex $(gC)$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)       :: harvest_type         !! Type of management that resulted
                                                                             !! in the harvest (unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)       :: harvest_cut          !! Type of cutting that was used for the harvest
                                                                             !! (unitless)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)     :: harvest_area         !! Harvested area (m^{2})
    REAL(r_std), DIMENSION(:,:), INTENT(inout)       :: age                  !! age (years)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)     :: turnover_time        !! turnover_time of grasses (days) 
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout) :: circ_class_biomass   !! Biomass of the componets of the model  
                                                                             !! tree within a circumference
                                                                             !! class @tex $(gC ind^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)       :: plant_status         !! Growth and phenological status of the plant
                                                                             !! Different stati are defined in constantes
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)     :: circ_class_n         !! Number of individuals in each circ class
                                                                             !! @tex $(m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)       :: when_growthinit      !! How many days ago was the beginning of 
                                                                             !! the growing season (days)


    !! 0.4 Local  variables

    LOGICAL, DIMENSION(npts)                         :: shed_rest            !! shed the remaining leaves? (true/false)
    INTEGER(i_std)                                   :: ivm,iele,ilage,ipts  !! Index (unitless)
    INTEGER(i_std)                                   :: ipar,icir,imbc,ij    !! Index (unitless)
    INTEGER(i_std), SAVE                             :: turn_count           !! Counter 
!$OMP THREADPRIVATE(turn_count)
    REAL(r_std), DIMENSION(npts,nvm)                 :: lai                  !! leaf area index @tex ($m^2 m^{-2}$)
    REAL(r_std), DIMENSION(npts,nvm)                 :: leaf_meanage         !! mean age of the leaves (days)
    REAL(r_std), DIMENSION(npts,ncirc)               :: dturnover            !! Intermediate variable for turnover ??
                                                                             !! @tex ($gC m^{-2}$) @endtex 
    REAL(r_std), DIMENSION(npts)                     :: vegstress_crit        !! critical moisture availability, function 
                                                                             !! of last year's moisture availability 
                                                                             !! (0-1, unitless)
    REAL(r_std), DIMENSION(npts)                     :: tl                   !! long term annual mean temperature, (C)
    REAL(r_std), DIMENSION(npts)                     :: t_crit               !! critical senescence temperature, function 
                                                                             !! of long term annual temperature (K) 
    REAL(r_std), DIMENSION(npts)                     :: sapconv              !! Sapwood conversion @tex ($gC m^{-2}$) 
                                                                             !! @endtex 
    REAL(r_std), DIMENSION(npts)                     :: hw_old               !! old heartwood mass @tex ($gC m^{-2}$) 
                                                                             !! @endtex 
    REAL(r_std), DIMENSION(npts)                     :: hw_new               !! new heartwood mass @tex ($gC m^{-2}$) 
                                                                             !! @endtex 
    REAL(r_std), DIMENSION(npts)                     :: lm_old               !! old leaf mass @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(npts)                     :: init_biomass         !! Biomass before turnover. This variable is used
                                                                             !! in IF-statements to ensure that the same 
                                                                             !! initial biomass is used for N and C 
                                                                             !! @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(npts,nleafages)           :: delta_lm             !! leaf mass change for each age class @tex 
                                                                             !! ($gC m^{-2}$) @endtex 
    REAL(r_std), DIMENSION(npts)                     :: turnover_rate        !! turnover rate (unitless) 
    REAL(r_std), DIMENSION(npts,nvm)                 :: new_turnover_time    !! instantaneous turnover time (days)
    REAL(r_std), DIMENSION(npts,nvm)                 :: harvest_time         !! Prescribed harvest time adjusted for the
                                                                             !! longterm temperature at the pixel (days)
    REAL(r_std)                                      :: branch_turn          !! turnover of branches (how much goes 
                                                                             !! to litter each day) (cf article CO2fix)
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements):: bm_old               !! old root mass
    REAL(r_std), DIMENSION(npts,nvm)                 :: sapwood_age          !! sapwood age (days)
    REAL(r_std), DIMENSION(npts)                     :: sw_old               !! old sapwood mass 
                                                                             !! (gC/(m**2 of nat/agri ground))
    REAL(r_std), DIMENSION(npts)                     :: sw_new               !! new sapwood mass 
                                                                             !! (gC/(m**2 of nat/agri ground)
    REAL(r_std), DIMENSION(npts,nvm,nmbcomp,nelements)&
                                                     :: check_intern         !! Contains the components of the internal
                                                                             !! mass balance chech for this routine
                                                                             !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)       :: closure_intern       !! Check closure of internal mass balance
                                                                             !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)       :: pool_start, pool_end !! Start and end pool of this routine 
                                                                             !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)                 :: veget_max_begin      !! temporary storage of veget_max to check
                                                                             !! area conservation (unitless, 0-1)
    REAL(r_std), DIMENSION(npts,nvm)                 :: leaf_turn_ageing     !! temporary variable for preparing output variable.
                                                                             !! Leaf turnover due to leaf ageing (excl. senescence)
                                                                             !! @tex $(gC m^{-2} dt^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm)                 :: last_plant_status    !! Plant status at the start of this routine
    REAL(r_std), DIMENSION(npts,nvm)                 :: doy_isenescent       !! date of year for isenescence
    REAL(r_std)                                      :: recycle
    REAL(r_std),DIMENSION(npts)                      :: temp_w
    INTEGER(i_std), DIMENSION(4)                     :: parts
    REAL(r_std), DIMENSION(npts,nvm)                 :: ratio_rm_gpp         !! temporary variable for send xios 
!_ ================================================================================================================================
                                                                             
    IF (firstcall_turnover) THEN
      !! Initialize local printlev
      ! Guillaume M. -- get_printlev calls getin_p, a COLLECTIVE bcast. Called here,
      ! inside the compute loop under a local firstcall, it hangs the run as soon as the
      ! ranks stop reaching it together. Use printlev, read once at startup by every rank.
      printlev_loc=printlev
    END IF

    IF (printlev_loc>=2) WRITE(numout,*) 'Entering turnover'

    !! Initialize
    ratio_rm_gpp(:,:) = zero
    WHERE (gpp_week .GT. min_stomate) 
      ratio_rm_gpp(:,:) = resp_maint_week(:,:)/gpp_week(:,:)
    ENDWHERE 

    ! Debug
    IF (printlev_loc>=4) THEN
       DO ipar = 1,nparts
          DO iele = 1,nelements
            DO icir =1, ncirc
               IF(icir == 1 .AND. iele == 1)&
               WRITE(numout,*) 'Biomass check 01: ',&
                    circ_class_biomass(test_grid,test_pft,1,ipar,iele) * &
                    circ_class_n(test_grid,test_pft,1)
            ENDDO
         ENDDO
      ENDDO
    ENDIF

!! 1. first call - output messages

    ! Initialize first call
    IF ( firstcall_turnover ) THEN

       IF(printlev_loc>=4) THEN
          WRITE(numout,*) 'turnover:'
          WRITE(numout,*) ' > minimum mean leaf age for senescence ',&
               '(days) (::min_leaf_age_for_senescence) : ',&
               min_leaf_age_for_senescence
       ENDIF

       turn_count = 0
       firstcall_turnover = .FALSE.

    ENDIF
    parts(1)=iroot
    parts(2)=ileaf
    parts(3)=ifruit
    parts(4)=isapabove

 !! 2. Initializations 

    !! 2.1 set output to zero
    turnover(:,:,:,:) = zero
    new_turnover_time(:,:) = zero
    harvest_time(:,:) = zero 
    doy_end_gs(:,:) = zero
    doy_isenescent(:,:) = zero
    leaf_turn_ageing(:,:) = zero

    ! Archive plant_status at the start of this routine
    last_plant_status(:,:) = plant_status(:,:)

    !! 2.2 Initialize check for mass balance closure
    !  The mass balance is calculated at the end of this routine
    !  in section 10.
    ! turnover is always equal to zero, so maybe that term doesn't need
    ! to be included here.
    IF (err_act.GT.1) THEN

       pool_start = zero

       DO iele = 1,nelements
          ! Biomass pool + turnover
          DO ipar = 1,nparts
             DO icir = 1,ncirc
                ! Initial biomass pool
                pool_start(:,:,iele) = pool_start(:,:,iele) + &
                     circ_class_biomass(:,:,icir,ipar,iele) * &
                     circ_class_n(:,:,icir) * veget_max(:,:)
             ENDDO
             ! Add turnover to the initial biomass pool
             pool_start(:,:,iele) = pool_start(:,:,iele) + &
                  (turnover(:,:,ipar,iele) * veget_max(:,:))
                
          ENDDO
       ENDDO

       ! The biomass harvest pool is expressed in gC pixel-1 So, it 
       ! shouldn't be multiplied by veget_max but it should be divided
       ! by area to obtain gC m-2.
       DO ivm = 1,nvm
          DO iele = 1,nelements
             pool_start(:,ivm,iele) = pool_start(:,ivm,iele) + &
                  SUM(SUM(harvest_pool(:,ivm,:,iele,:),3),2) / &
                  (area(:) * contfrac(:))
          ENDDO
       ENDDO

       !! 1.3 Initialize check for area conservation
       veget_max_begin(:,:) = veget_max(:,:)

    ENDIF ! err_act.GT.1

    !+++CHECK+++
    ! No branch turnover in CAN but it is not clear 
    ! why that would be justified. Could be accounted
    ! for somewhere else (seems unlikely).
!!$    ! Compute the turnover of branches 
!!$    ! (2.5% per year : Orchidee standard, Masera 2003, Lehtonen 2004,
!!$    ! DeAngelis 1981)
!!$    branch_turn = 0.025/(one_year/dt)*ss_branch_turn
    !++++++++++++

    ! Calculate lai
    DO ivm = 2,nvm
       DO ipts=1,npts
          lai(ipts,ivm) = biomass_to_lai(SUM(circ_class_biomass(ipts,ivm,:,ileaf,icarbon)),ivm)
       ENDDO
    ENDDO
    lai(:,1) = zero

    !! 2.3 Initialize check for surface area conservation
    !  Veget_max is a INTENT(in) variable and can therefore
    !  not be changed during the course of this subroutine
    !  No need to check whether the subroutine preserves the
    !  total surface area of the pixel.

    !! 2.4 Recalculate mean leaf age
    !  Mean leaf age is recalculated from leaf ages of separate leaf age
    !  classes. Leaf age is used as a criterion in several of the
    !  following IF/WHERE loops. The mean leaf age is calculated using
    !  the following equation:
    !  \latexonly
    !  \input{turnover_lma_update_eqn1.tex}
    !  \endlatexonly
    !  \n
    leaf_meanage(:,:) = zero

    DO ilage = 1, nleafages
       leaf_meanage(:,:) = leaf_meanage(:,:) + &
            leaf_age(:,:,ilage) * leaf_frac(:,:,ilage)
    ENDDO

    ! Debug
    IF (printlev_loc>=4) THEN
       DO ipar = 1,nparts
          DO iele = 1,nelements
            DO icir =1, ncirc
               IF(icir == 1 .AND. iele == 1)&
               WRITE(numout,*) 'Biomass/turnover check 02: ',&
                    circ_class_biomass(test_grid,test_pft,1,ipar,iele) * &
                    circ_class_n(test_grid,test_pft,1),&
                    turnover(test_grid,test_pft,ipar,iele)
            ENDDO
         ENDDO
      ENDDO
    ENDIF

!! 3. Climatic senescence

    ! Three different types of "climatic" leaf senescence,
    ! that do not change the age structure. 

    DO ivm = 2,nvm ! Loop over # PFTs
       !! 3.1 Determine if there is climatic senescence. 
       !  The climatic senescence can be of three types: sensitivity to
       !  cold temperatures, to lack of water, or both. If meteorological
       !  conditions are fulfilled, a flag senescence is set to TRUE.
       !  Evergreen species do not experience climatic senescence.
       SELECT CASE ( senescence_type(ivm) )

       CASE ('crop' )

          !+++CHECK+++
          ! This is the original code. CAN took a different approach (see below)
          ! but that different approach need to be carefuly revisited
          ! Crop senescence is based on a GDD criterium as in crop models.
          ! The problem with this approach is the use of a global gdd_senescence. For 
          ! C3 crops this value is 2500 but is in the temperate zone not reached by 
          ! the end of summer. As such crops keep their leaves through winter and
          ! only exceed this threshold in the next spring. If the harvest finally 
          ! takes place, the conditions are good to replant. The fallow season is 
          ! therefore reduced to a single day. An alternative approach is now used.
!!$          WHERE ( (SUM(circ_class_biomass(:,ivm,:,ileaf,icarbon),2) .GT. zero ) .AND. &
!!$               ( leaf_meanage(:,ivm) .GT. min_leaf_age_for_senescence(ivm) )&
!!$               .AND.( gdd_from_growthinit(:,ivm) .GT. gdd_senescence(ivm)))
!!$             
!!$              plant_status(:,ivm) = isenescent
!!$              doy_isenescent(:,ivm) = julian_diff
!!$
!!$          ENDWHERE
          
          
          ! This is pure plumbing under time pressure (says the one who
          ! introduced it). A growing season should
          ! last about 4 months  where the longterm temperature is 282 K
          ! (crude observation from N. France - harvest_time_b) before 
          ! senescence can set in unless the growing degrees are reached. 
          ! In colder and warmer places the growing season is shortened 
          ! or lengthened. Because in warmer places the growing season 
          ! lasts longer, crop can be planted that required more light. 
          ! The very long growing season in warm places does not account 
          ! for wet/dry seasons but should be overruled by the 
          ! gdd_senescence criterium.
          harvest_time(:,ivm) = harvest_time_a(ivm) + &
               (t2m_longterm(:) - harvest_time_b(ivm)) * &
               ABS(t2m_longterm(:) - harvest_time_b(ivm))
   
          ! Ecophysiological approach
          WHERE ( ( plant_status(:,ivm) .EQ. ipresenescence .OR. &
               plant_status(:,ivm) .EQ. icanopy)  .AND. &
               leaf_meanage(:,ivm) .GT. min_leaf_age_for_senescence(ivm) .AND. &
               ( gdd_from_growthinit(:,ivm) .GT. gdd_senescence(ivm) .OR. &
               when_growthinit(:,ivm) .GT. harvest_time(:,ivm) ) )
             
             plant_status(:,ivm) = isenescent
             doy_isenescent(:,ivm) = julian_diff
             
          ENDWHERE
          
       CASE ( 'cold' )

          !! 3.1.2 Summergreen species
          !  Climatic senescence is triggered by sensitivity to cold.
          !  Critical temperature for senescence may depend on long term
          !  annual mean temperature
          tl(:) = t2m_longterm(:) - ZeroCelsius
          t_crit(:) = ZeroCelsius + senescence_temp(ivm,1) + &
               tl(:) * senescence_temp(ivm,2) + &
               tl(:)*tl(:) * senescence_temp(ivm,3)

          ! Mixed approach
          WHERE ( (plant_status(:,ivm) .EQ. ipresenescence .OR. &
               plant_status(:,ivm) .EQ. icanopy) .AND. &
               leaf_meanage(:,ivm) .GT. min_leaf_age_for_senescence(ivm) .AND. &
               (ratio_rm_gpp(:,ivm) .GT. senescence_ratio(ivm) .OR. &
               (t2m_month(:) .LT. t_crit(:) .AND. &
               t2m_week(:) .LT. t2m_month(:))) )
             
             plant_status(:,ivm) = isenescent
             doy_isenescent(:,ivm) = julian_diff
             
          END WHERE

          ! Fail safe option. If something goes wrong with t_crit, phenology 
          ! will get screwed up, leaf age is getting too high, vcmax is 
          ! decreasing, etc. We added an option that should prevent this from 
          ! happening. 1.2 is an arbitrary threshold.
          WHERE ( ( plant_status(:,ivm).EQ.ipresenescence .OR. &
                    plant_status(:,ivm).EQ.icanopy).AND. &
               when_growthinit(:,ivm).GT.1.2*leaf_age_crit(:,ivm))

             plant_status(:,ivm) = isenescent
             doy_isenescent(:,ivm) = julian_diff

          ENDWHERE
                
       CASE ( 'dry' )

          !! 3.1.3 Raingreen species
          !  Climatic senescence is triggered by sensitivity to lack of
          !  water availability
          vegstress_crit(:) = &
               MIN( MAX( minvegstress_lastyear(:,ivm) + hum_frac(ivm) * &
               ( maxvegstress_lastyear(:,ivm) - minvegstress_lastyear(:,ivm) ), &
               senescence_hum(ivm) ), nosenescence_hum(ivm) )

          WHERE ( (plant_status(:,ivm) .EQ. ipresenescence .OR. &
               plant_status(:,ivm) .EQ. icanopy) .AND. &
               (leaf_meanage(:,ivm).GT.min_leaf_age_for_senescence(ivm)).AND.&
               (ratio_rm_gpp(:,ivm) .GT. senescence_ratio(ivm) .OR. &
               vegstress_week(:,ivm) .LT. vegstress_crit(:)) )
             
             plant_status(:,ivm) = isenescent
             doy_isenescent(:,ivm) = julian_diff
             
          ENDWHERE

          ! Fail safe option. In some wet years (or years with very low LAI) the
          ! vegstress_week never drops below vegstress_crit at some isolated pixels
          ! If that is the case phenology is getting screwed up, leaf age is 
          ! getting too high, vcmax is decreaseing, etc. We added an option that 
          ! should prevent this to happen. 1.2 is an arbitrary threshold.
          WHERE ((plant_status(:,ivm).EQ.ipresenescence .OR. &
                  plant_status(:,ivm).EQ.icanopy ) .AND. &
               when_growthinit(:,ivm).GT.1.2*leaf_age_crit(:,ivm))
             
             plant_status(:,ivm) = isenescent
             doy_isenescent(:,ivm) = julian_diff

          ENDWHERE

       CASE ( 'mixed' )

          !! 3.1.4 Mixed criterion: Climatic senescence is triggered
          !!       by sensitivity to temperature or to lack of water  
          vegstress_crit(:) = MIN( MAX( minvegstress_lastyear(:,ivm) + &
               hum_frac(ivm) * (maxvegstress_lastyear(:,ivm) - &
               minvegstress_lastyear(:,ivm) ), senescence_hum(ivm) ), &
               nosenescence_hum(ivm) )
          tl(:) = t2m_longterm(:) - ZeroCelsius
          t_crit(:) = ZeroCelsius + senescence_temp(ivm,1) + &
               tl(:) * senescence_temp(ivm,2) + &
               tl(:)*tl(:) * senescence_temp(ivm,3)

          IF ( is_tree(ivm) ) THEN

             !+++CHECK+++
             CALL ipslerr_p(3,'You are the first the use a mixed', &
                  'senescence model for a tree PFT', &
                  'Check the conditions contained in the WHERE in', &
                  'stomate_turnover.f90')
             
             WHERE ( ( plant_status(:,ivm) .EQ. ipresenescence .OR. &
                       plant_status(:,ivm) .EQ. icanopy) .AND. &
                  (leaf_meanage(:,ivm).GT.min_leaf_age_for_senescence(ivm)).AND.&
                  ( ( vegstress_week(:,ivm) .LT. vegstress_crit(:) ) .OR. &
                  ( ( t2m_month(:) .LT. t_crit(:) ) .AND. &
                  ( t2m_week(:) .LT. t2m_month(:) ) ) ) )

                plant_status(:,ivm) = isenescent
                doy_isenescent(:,ivm) = julian_diff

             ENDWHERE
             !+++++++++++

             
          ELSE

             IF (natural(ivm)) THEN
                ! Grasses

                !++++CHECK++++
                ! This is still a patch in the trunck (Tag 3.0). The variable
                ! have been changed since ORCHIDEE-CN-CAn so we tried to follow
                ! the idea of the patch rather than its identic code. The idea is 
                ! that, if the leaf age gets too high (i.e., too high of fraction 
                ! of leaves are in the oldest leaf age class), senescence is 
                ! triggered.  In ORC2.2, this is done with lgrassleafage=TRUE. 
                ! for sensesence type "mixed", and it triggers the senescence flag.  
                ! In the trunk, this is the plant_status variable. In ORC2.2
                ! the patch increases turnover_time. In the trunk, we also increased 
                ! the turnover_time, but in a different way.
               
                ! Physiological approach
                WHERE ( (plant_status(:,ivm) .EQ. ipresenescence .OR. &
                     plant_status(:,ivm) .EQ. icanopy) .AND. &
                     ratio_rm_gpp(:,ivm) .GT. senescence_ratio(ivm) )
                   
                   ! Shed leaves, roots and fruits at a high rate = senescence
                   turnover_time(:,ivm,ileaf) = leaffall(ivm)
                   turnover_time(:,ivm,iroot) = leaffall(ivm)
                   turnover_time(:,ivm,ifruit) = leaffall(ivm)
                   
                   ! Slowly turnover the structural carbon this allows to store 
                   ! reserves. This structural pool is essential for the
                   ! functioning of the labile and reserve pools
                   turnover_time(:,ivm,isapabove) = longevity_eff_sap(:,ivm)
                   plant_status(:,ivm) = isenescent                     
                   doy_isenescent(:,ivm) = julian_diff
                   
                ENDWHERE               
            
             ELSE

                ! croplands
                CALL ipslerr_p (3,'Select a different senescence model',&
                     'croplands now have their own senescence',&
                     'It is called crop','')

             END IF ! is_natural

          END IF ! is_tree

       CASE ( 'none' )

          !! 3.1.5 Evergreen species 
          !  Evergreen species do not experience climatic senescence
          
       CASE default

          !! 3.1.6 Other cases
          !  In case no climatic senescence type is recognized.
          WRITE(numout,*) '  turnover: don''t know how to treat this PFT.'
          WRITE(numout,*) '  number (::ivm) : ',ivm
          WRITE(numout,*) '  senescence type (::senescence_type(ivm)) : ',&
               senescence_type(ivm)
          CALL ipslerr_p (3,'stomate_turnover',&
               'turnover: don''t know how to treat this PFT','case 1','')

       END SELECT

       ! Debug
       IF (printlev_loc >= 4 .AND. ivm .EQ. test_pft) THEN
          WRITE(numout,*) 'phenology, plant_status, ',plant_status(:,ivm)
       ENDIF
       !-

       !! 3.2 Drop leaves and roots, plus stems and fruits for grasses
       IF ( is_tree(ivm) ) THEN
          !! 3.2.1 Trees in climatic senescence lose their fine roots at
          !!       the same rate as they lose their leaves. 
          !  The rate of biomass loss of both fine roots and leaves
          !  is presribed through the equation:
          !  \latexonly
          !  \input{turnover_clim_senes_biomass_eqn4.tex}
          !  \endlatexonly
          !  \n
          !  with ::leaffall(ivm) a PFT-dependent time constant which is
          !  given in ::stomate_constants),
          ! Calculate stand-level turnover for
          ! both carbon and nitrogen (gC tree-1)).
          ! For nitrogen, we first calculate the turnover that
          ! is lost from the respective biomass pools, i.e.,
          ! all the nitrogen is lost from the leaf pool. Then
          ! recycle some of the leaf nitrogen into the labile pool
          ! and finaly correct the turnover pool for the recycled
          ! nitrogen.

          DO iele =1, nelements   
            DO ij = 1, 2 ! iroot and ileaf
                ipar=parts(ij)
                IF (ipar .EQ. ileaf .AND. iele .EQ. initrogen) THEN
                    recycle=recycle_leaf(ivm)
                ELSEIF (ipar .EQ. iroot .AND. iele .EQ. initrogen) THEN 
                    recycle=recycle_root(ivm)
                ELSE
                    recycle=zero
                ENDIF
                dturnover(:,:) = zero
                DO icir=1, ncirc
                    WHERE ( plant_status(:,ivm) .EQ. isenescent )
                        dturnover(:,icir) = dt/ leaffall(ivm)* &
                            circ_class_biomass(:,ivm,icir,ipar,iele)
                    ENDWHERE
                ENDDO
                circ_class_biomass(:,ivm,:,ilabile,iele)= &
                    circ_class_biomass(:,ivm,:,ilabile,iele) + &
                    dturnover(:,:)*recycle
                turnover(:,ivm,ipar,iele) = turnover(:,ivm,ipar,iele)+ &
                    SUM(dturnover(:,:) * circ_class_n(:,ivm,:),2) * &
                     ( un - recycle )
                circ_class_biomass(:,ivm,:,ipar,iele) = &
                    circ_class_biomass(:,ivm,:,ipar,iele) - &
                    dturnover(:,:)
            ENDDO
          ENDDO 

       ELSE    
          !! 3.2.2.1 crops with 'crop' phenological model
          IF (senescence_type(ivm) .EQ. 'crop') THEN
             DO ipts = 1,npts
                IF (plant_status(ipts,ivm) .EQ. isenescent) THEN
                   ! Crops are planted every year. So make sure to remove
                   ! everything at the end of senescence which for crops
                   ! is actually harvest. Next stomate_prescribe should
                   ! prescribe a new crop. If sapwood is not removed, some
                   ! is left on site and the model will try to grow a crop.
                   ! However, there are not enough reserves left so the 
                   ! crop will die. This results in one year of normal crop
                   ! growth (following the prescribe) and then one year of
                   ! almost no growth (the year starting from the too low
                   ! reserves). This issue has been fixed.
                   CALL crop_harvest(npts, ipts, ivm, dt, &
                        veget_max, turnover, harvest_pool, harvest_type, &
                        harvest_cut, harvest_area, &
                        circ_class_biomass,circ_class_n,leaf_meanage)
                   
                   ! Update plant_status. Crops were harvested. No living
                   ! biomass was left on-site thus we will have to 
                   ! prescribe a new vegetation during the next growing
                   ! season.
                   plant_status(ipts,ivm) = idead

                ENDIF
             ENDDO

          !! 3.2.2.2 grass based on 'mixed' phenological model
          ELSEIF (senescence_type(ivm) .EQ. 'mixed') THEN
            DO iele = 1,nelements ! carbon and nitrogen (gC m-2)
                DO ij = 1, 4 ! ileaf, iroot, ifruit, isapabove
                    ipar=parts(ij)
                    IF (ipar .EQ. ileaf .AND. iele .EQ. initrogen) THEN
                        recycle=recycle_leaf(ivm)
                    ELSEIF (ipar .EQ. iroot .AND. iele .EQ. initrogen) THEN
                        recycle=recycle_root(ivm)
                    ELSE
                        recycle=zero
                    ENDIF
                    ! First calculate the total turnover. Part of the
                    ! turnover will be recycled and stored in the labile
                    ! pool. The remaining mass is added to the turnover
                    ! pool that will eventually become litter. Once the
                    ! turnover is known the biomass pool can be updated.
                    dturnover(:,:) = zero
                    WHERE (plant_status(:,ivm) .EQ. isenescent)
                        dturnover(:,1) = circ_class_biomass(:,ivm,1,ipar,iele)*&
                            dt/turnover_time(:,ivm,ipar)
                    ENDWHERE
                    circ_class_biomass(:,ivm,1,ilabile,iele) = &
                        circ_class_biomass(:,ivm,1,ilabile,iele) + &
                        recycle * dturnover(:,1)
                    turnover(:,ivm,ipar,iele) = turnover(:,ivm,ipar,iele)+ &
                        dturnover(:,1) * circ_class_n(:,ivm,1)* (un-recycle)
                    circ_class_biomass(:,ivm,1,ipar,iele) = &
                        circ_class_biomass(:,ivm,1,ipar,iele) - &
                        dturnover(:,1)
                    ! Debug
                    IF (printlev_loc>4) THEN
                       IF(ivm == test_pft)THEN
                          WRITE(numout,*) 'total turnover, ', test_grid,test_pft,ipar,iele
                          WRITE(numout,*) 'dturnover, ', dturnover(test_grid,1)
                          WRITE(numout,*) 'recycle, ', recycle
                          WRITE(numout,*) 'turnover, ', turnover(test_grid,test_pft,ipar,iele)
                          WRITE(numout,*) 'circ_class_biomass, ', &
                               circ_class_biomass(test_grid,test_pft,1,ipar,iele)
                       ENDIF
                    ENDIF
                    
                ENDDO   
             ENDDO
 
             ! Debug
             IF (printlev_loc>4) THEN
                WRITE(numout,*) 'test - recycle, ', ivm, recycle_root(ivm), &
                     recycle_leaf(ivm)
             ENDIF
             !-
          ELSEIF  (senescence_type(ivm) .EQ. 'none') THEN
             !! 3.2.2.3 Grasses modeled as an evergreen biome
             turnover(:,ivm,:,:) = zero
          ELSE
             !! 3.2.2.4 Conceptual problem
             WRITE(numout,*) 'ERROR: senenescence type not known - 2'
             CALL ipslerr_p(3,'stomate_turnover',&
                  'turnover: senescence type not known','case 2','')
          END IF
       ENDIF ! tree/grass
    ENDDO ! loop over PFTs

    ! Debug
    IF (printlev_loc>=4) THEN
       DO ipar = 1,nparts
          DO iele = 1,nelements
            DO icir =1, ncirc
               IF(icir == 1 .AND. iele == 1)&
               WRITE(numout,*) 'Biomass/turnover check 03: ',&
                    circ_class_biomass(test_grid,test_pft,1,ipar,iele) * &
                    circ_class_n(test_grid,test_pft,1),&
                    turnover(test_grid,test_pft,ipar,iele)
            ENDDO
         ENDDO
      ENDDO
    ENDIF

!! 4. Leaf fall
    !  At a certain age, leaves fall off, even if the climate would allow a 
    !  green plant all year round. Even if the meteorological conditions are 
    !  favorable for leaf maintenance, plants, and in particular, evergreen 
    !  trees, have to renew their leaves simply because the old leaves become 
    !  inefficient. Another reason for leaf fall may be waterstress. When 
    !  the plant experiences water stress some of the leaves will die.
    !  Wstress is calculated from the ratio between stressed and unstressed
    !  GPP. If the wstress is less than 1, this implies that we have to maintain
    !  a complete canopy (and therefore have the Ra cost) but that only part 
    !  of the canopy will contribute to GPP. Although this is exactely what 
    !  happens on a warm summer day, what we are trying to account for here is
    !  the impact of a lasting drought (days to weeks). Adaptation to the
    !  growing conditions are accounted for through KF (see allocation)
    !  Roots, fruits (and stems) follow leaves. The decay rate varies with 
    !  leaf age. Note that the plant is not declared senescent in this case 
    !  (wchich is important for allocation: if the plant looses leaves because 
    !  of their age, it can renew them). 
    !  Notice that we do not change the reserve pools here (ilabile and 
    !  icarbres). We assume that the tree will move these pools out of old 
    !  leaves and therefore the carbon will not be lost.
    !  The leaf turnover rate due to aging of leaves is calculated using 
    !  the following equation:
    !  \latexonly
    !  \input{turnover_age_senes_biomass_eqn5.tex}
    !  \endlatexonly
    !  \n
    DO ivm = 2,nvm 

       !! 4.1 Calculate critical leaf age       
       !  save old leaf mass
       lm_old(:) = SUM(circ_class_biomass(:,ivm,:,ileaf,icarbon)* &
            circ_class_n(:,ivm,:),2)

       !! initialize leaf mass change in age class
       delta_lm(:,:) = zero
            
       !! 4.2 Turnover for different leaf age classes
       DO ilage = 1,nleafages

          turnover_rate(:) = zero

          ! Age-driven turnover can only occur when there is a canopy
          ! and if the leaves exceeds a certain age. The plant_status 
          ! of a crop is set for only one day to isenescent. The crop 
          ! is harvested on that day. For the deciudous trees turn-over during senesce
          ! is accounted for above. For all these reasons, the turnover
          ! calculated here was limited to the plant_status = icanopy.
          WHERE ( (plant_status(:,ivm) .EQ. icanopy .OR. &
                  plant_status(:,ivm) .EQ. ipresenescence ).AND. &
               leaf_age(:,ivm,ilage) .GT. leaf_age_crit(:,ivm) / deux)

             ! Not clear where this equation comes from or what it is
             ! based on. 
             turnover_rate(:) =  &
                  MIN( 0.99_r_std, dt / ( leaf_age_crit(:,ivm) * &
                  ( leaf_age_crit(:,ivm) / leaf_age(:,ivm,ilage) )**quatre ) )

          ENDWHERE
          
          IF ( is_tree(ivm) ) THEN

             ! Stand level turnover (gC m-2) for leaves
             ! Water stress is increasing (note that ::wstress is thus 
             ! decreasing) so the LAI will be lowered if the decrease is 
             ! small, this will result in a cosmetic change in LAI. 
             ! If water stress is increasing over a the course of a week, 
             ! this decrease may become important. The correction is 
             ! cummulative, for example, ::biomass(ileaf) is decreased 
             ! by 10% in day one and the remaining fraction of 
             ! ::biomass(ileaf) can be decreased by another 10% the next 
             ! day, etc ... The choice to multiply with ::wstress_month
             ! is arbitrary, we could have well used ::wstress_week. By
             ! using the monthly value we hope a smoother decrease 
             ! compared to using ::wstress_week or ::wstress_day.
             ! Update biomass and turnover. During turnover all carbon
             ! in leaves, roots and fruits will be lost. This is not the
             ! case for N. Plants are more careful with teir nitrogen
             ! part of the nitrogen will be resorbed prior to leaf and
             ! root turnover. The fraction of nitrogen that is resorbed
             ! is given by ::recycle_leaf. The resorbed nitrogen is
             ! moved into the labile pool the rest will be lost in
             ! turnover. Note that no resorprion happens for fruits so
             ! all nitrogen in a fruit is lost at turnover.
             ! It is not clear why roots follow the leaves. A turnover
             ! rate for roots based on longevity_eff_root can be easily calculated
             ! See crops and grasses where we made use of longevity_eff_root
             ! Note that this turnover is a bit arbitrairy there is
             ! no reason why the same fraction of roots and leaves
             ! should die. This formulation will disturb the allometric
             ! allocation and will need to be restored by phenological
             ! growth in the allocation routine.
             !+++CHECK+++
!!$             WHERE ( wstress_month(:,ivm) .GE. un-min_stomate)
!!$                temp_w(:) = turnover_rate(:) / &
!!$                     ((min_stomate)**(1/tau_hum_month))
!!$             ELSEWHERE( wstress_month(:,ivm) .LT. un-min_stomate)
!!$                temp_w(:) = turnover_rate(:) / &
!!$                     ((un-wstress_month(:,ivm))**(1/tau_hum_month))
!!$             END WHERE
!!$             
!!$             ! temp_w can be larger than 1. If so there is a chance that we
!!$             ! don have enough leaves to shed. Make sure this does npt happen
!!$             ! and truncate temp_w to 1.
!!$             WHERE ( temp_w(:) .GE. un)
!!$                ! No longer account for wstress
!!$                temp_w(:) = turnover_rate(:)
!!$             ELSE
!!$                temp_w(:) = temp_w(:) * turnover_rate(:)
!!$             ENDWHERE
             
             ! This line replaces all lines above. These lines are
             ! not present in the trunk. Not clear who added them
             ! and for what reason
             temp_w(:) = turnover_rate(:)
             !+++++++++++

             ! Note that delta_lm is in gC m-2 (dturnover is in gC tree-1)
             delta_lm(:,ilage) = - temp_w(:) * leaf_frac(:,ivm,ilage) * &
                  SUM(circ_class_biomass(:,ivm,:,ileaf,icarbon) * & 
                  circ_class_n(:,ivm,:),2)
             ! Accumulate the leaf turnover due to ageing in gC m-2
             leaf_turn_ageing(:,ivm) = leaf_turn_ageing(:,ivm) + delta_lm(:,ilage)

             DO iele=1, nelements
                ! Turnover (gC tree-1): leaves, roots and fruit
                ! Sapwood turnover is accounted for separatly later
                ! in this routine
                DO ij =1, 3 !iroot, ileaf, ifruit
                    ipar=parts(ij)
                    IF (ipar .EQ. ileaf .AND. iele .EQ. initrogen) THEN
                        recycle=recycle_leaf(ivm)
                    ELSEIF (ipar .EQ. iroot .AND. iele .EQ. initrogen) THEN
                        recycle=recycle_root(ivm)
                    ELSE
                        recycle= zero
                    ENDIF
                    ! Note that all of this code is still within an loop over te
                    ! leaf age classes. Avoid double counting by accounting for 
                    ! the fractions of the leaf age classes when calculating root,
                    ! sapwood and fruit turnover as well
                    dturnover(:,:) = zero
                    DO icir=1, ncirc
                       IF(ipar .EQ. iroot) THEN
                          dturnover(:,icir) = leaf_frac(:,ivm,ilage) * dt / &
                               longevity_eff_root(:,ivm) * &
                               circ_class_biomass(:,ivm,icir,ipar,iele)
                       ELSEIF (ipar .EQ. ileaf) THEN 
                          dturnover(:,icir) = leaf_frac(:,ivm,ilage) * temp_w(:) * &
                               circ_class_biomass(:,ivm,icir,ipar,iele)
                       ELSEIF(ipar .EQ. ifruit) THEN
                          dturnover(:,icir) = leaf_frac(:,ivm,ilage) * dt / &
                               longevity_fruit(ivm) * &
                               circ_class_biomass(:,ivm,icir,ipar,iele)
                       ENDIF
                    ENDDO
                    ! Nitrogen resorption
                    circ_class_biomass(:,ivm,:,ilabile,iele) = &
                         circ_class_biomass(:,ivm,:,ilabile,iele) + &
                         recycle * dturnover(:,:)
                    !the stand level turnover (gC m-2):
                    turnover(:,ivm,ipar,iele) = turnover(:,ivm,ipar,iele) + &
                         SUM(dturnover(:,:)* circ_class_n(:,ivm,:),2) * (un - recycle)
                    ! Update circ_class_biomass
                    circ_class_biomass(:,ivm,:,ipar,iele) = &
                         circ_class_biomass(:,ivm,:,ipar,iele) - &
                         dturnover(:,:)
                 ENDDO
              ENDDO
             
           ELSEIF ( .NOT. is_tree(ivm) .AND. .NOT. natural(ivm) ) THEN

              ! For crops the sapwood acts as the structure to supports
              ! leaves and roots and store the reserves and labile.
              ! Because the growing season is short for crops, the 
              ! turnover of tissues should be kept to a minimum or the GPP
              ! and NPP will be too low.
              ! Note that dturnover was redefined so that it can be used
              ! for both N and C
              !+++CHECK+++
!!$              WHERE ( wstress_month(:,ivm) .GE. un-min_stomate)
!!$                 temp_w(:) = turnover_rate(:) / &
!!$                      ((min_stomate)**(1/tau_hum_month))
!!$              ELSEWHERE( wstress_month(:,ivm) .LT. un-min_stomate)
!!$                 temp_w(:) = turnover_rate(:) / &
!!$                      ((un-wstress_month(:,ivm))**(1/tau_hum_month))
!!$              END WHERE
!!$              ! Do we have enough biomass to satisfy the turnover?
!!$              WHERE ( temp_w(:) .GE. un)
!!$                 ! No longer account for wstress
!!$                 temp_w(:) = turnover_rate(:)
!!$              ELSEWHERE
!!$                 temp_w(:) = temp_w(:)*turnover_rate(:)
!!$              ENDWHERE              
              ! Most basic approach
              temp_w(:) = turnover_rate(:)
              !+++++++++++

              ! Note that delta_lm is in gC m-2 (dturnover is in gC tree-1)
              delta_lm(:,ilage) = - temp_w(:) * leaf_frac(:,ivm,ilage) * &
                   circ_class_biomass(:,ivm,1,ileaf,icarbon)*&
                   circ_class_n(:,ivm,1)
              
              ! Accumulate the leaf turnover due to ageing
              leaf_turn_ageing(:,ivm) = leaf_turn_ageing(:,ivm) + delta_lm(:,ilage)

              DO iele=1, nelements
                 ! Turnover (gC plant-1):
                DO ij =1,4 ! ileaf, iroot, ifruit, isapabove
                   ipar=parts(ij)
                   IF (ipar .EQ. ileaf .AND. iele .EQ. initrogen) THEN
                      recycle=recycle_leaf(ivm)
                   ELSEIF (ipar .EQ. iroot .AND. iele .EQ. initrogen) THEN
                      recycle=recycle_root(ivm)
                   ELSE
                      recycle= zero
                   ENDIF
                   ! Note that all of this code is still within an loop over te
                   ! leaf age classes. Avoid double counting by accounting for 
                   ! the fractions of the leaf age classes when calculating root,
                   ! sapwood and fruit turnover as well.
                   dturnover(:,:) = zero
                   IF(ipar .EQ. iroot) THEN
                      dturnover(:,1) = circ_class_biomass(:,ivm,1,ipar,iele)*&
                           dt / longevity_eff_root(:,ivm) * leaf_frac(:,ivm,ilage)
                   ELSEIF (ipar .EQ. isapabove) THEN
                      dturnover(:,1) = circ_class_biomass(:,ivm,1,ipar,iele)*&
                           dt / longevity_eff_sap(:,ivm) * leaf_frac(:,ivm,ilage)
                   ELSEIF (ipar .EQ. ileaf) THEN
                      dturnover(:,1) = circ_class_biomass(:,ivm,1,ipar,iele)*&
                           temp_w(:) * leaf_frac(:,ivm,ilage)
                   ENDIF
                   ! Resorb some of the nitrogen
                   circ_class_biomass(:,ivm,1,ilabile,iele) = &
                        circ_class_biomass(:,ivm,1,ilabile,iele) + &
                        recycle * dturnover(:,1)
                   turnover(:,ivm,ipar,iele) = turnover(:,ivm,ipar,iele)+ &
                        dturnover(:,1)*circ_class_n(:,ivm,1)*(un - recycle)
                   ! Update circ_class_biomass
                   circ_class_biomass(:,ivm,1,ipar,iele) = &
                        circ_class_biomass(:,ivm,1,ipar,iele) - &
                        dturnover(:,1)
                ENDDO
             ENDDO
             
          ELSEIF ( .NOT. is_tree(ivm) .AND. natural(ivm) ) THEN

             ! In CAN grasses were considered managed ecosystem and the turnover
             ! went into a harvest pool. Here we follow the previous trunk
             ! definitions and thus assume grasslands are unamanged. The turnover 
             ! should end up in the litter and the within-season turnover should 
             ! be relatively low 10-20% compared to a manage systsem where the
             ! turnover could be much higher (100%) due to mowing and grazing.
             !+++CHECK+++
!!$             WHERE ( wstress_month(:,ivm) .GE. un-min_stomate)
!!$                temp_w(:) = leaf_frac(:,ivm,ilage) * turnover_rate(:) / &
!!$                     ((min_stomate)**(1/tau_hum_month))
!!$             ELSEWHERE( wstress_month(:,ivm) .LT. un-min_stomate)
!!$                temp_w(:) = leaf_frac(:,ivm,ilage) * turnover_rate(:) / &
!!$                     ((un-wstress_month(:,ivm))**(1/tau_hum_month))
!!$             END WHERE
!!$             ! Do we have enough biomass to satisfy the turnover?
!!$             WHERE ( temp_w(:) .GE. un)
!!$                ! No longer account for wstress
!!$                 temp_w(:) = turnover_rate(:)
!!$              ELSEWHERE
!!$                 temp_w(:) = temp_w(:)*turnover_rate(:)
!!$              ENDWHERE              
              ! Most basic approach
              temp_w(:) = turnover_rate(:)
              !+++++++++++

              ! Note that delta_lm is in gC m-2 (dturnover is in gC tree-1)
              delta_lm(:,ilage) = - temp_w(:) * leaf_frac(:,ivm,ilage) * &
                   circ_class_biomass(:,ivm,1,ileaf,icarbon)*&
                   circ_class_n(:,ivm,1)

              ! Accumulate the leaf turnover due to ageing
              leaf_turn_ageing(:,ivm) = leaf_turn_ageing(:,ivm) + delta_lm(:,ilage)

              DO iele=1, nelements
                 ! Turnover (gC tree-1):
                 DO ij =1,4 ! ileaf, iroot, ifruit, isapabove
                    ipar=parts(ij)
                    IF (ipar .EQ. ileaf .AND. iele .EQ. initrogen) THEN
                       recycle=recycle_leaf(ivm)
                    ELSEIF (ipar .EQ. iroot .AND. iele .EQ. initrogen) THEN
                       recycle=recycle_root(ivm)
                    ELSE
                       recycle= zero
                    ENDIF
                    ! Note that all of this code is still within an loop over te
                    ! leaf age classes. Avoid double counting by accounting for 
                    ! the fractions of the leaf age classes when calculating root,
                    ! sapwood and fruit turnover as well.
                    dturnover(:,:) = zero
                    IF(ipar .EQ. iroot) THEN
                       dturnover(:,1) = circ_class_biomass(:,ivm,1,ipar,iele)*&
                            dt / longevity_eff_root(:,ivm) * leaf_frac(:,ivm,ilage)
                    ELSEIF (ipar .EQ.isapabove) THEN
                       dturnover(:,1) =circ_class_biomass(:,ivm,1,ipar,iele)*&
                            dt / longevity_eff_sap(:,ivm) * leaf_frac(:,ivm,ilage)
                    ELSEIF (ipar .EQ.ileaf) THEN
                       dturnover(:,1) = temp_w(:) * leaf_frac(:,ivm,ilage) * &
                            circ_class_biomass(:,ivm,1,ipar,iele)
                    ELSEIF (ipar .EQ.ifruit) THEN
                       dturnover(:,1) = temp_w(:) * leaf_frac(:,ivm,ilage) * &
                            circ_class_biomass(:,ivm,1,ipar,iele)
                    ENDIF
                    ! Nitrogen resorption
                    circ_class_biomass(:,ivm,1,ilabile,iele) = &
                         circ_class_biomass(:,ivm,1,ilabile,iele) + &
                         recycle  * dturnover(:,1)
                    ! Stand level turnover (gC m-2) 
                    turnover(:,ivm,ipar,iele) = turnover(:,ivm,ipar,iele)+ &
                         dturnover(:,1)*circ_class_n(:,ivm,1) * & 
                         (un - recycle)
                    ! Update circ_class_biomass
                    circ_class_biomass(:,ivm,1,ipar,iele) = &
                         circ_class_biomass(:,ivm,1,ipar,iele) - &
                         dturnover(:,1)
                    
                 ENDDO
              ENDDO
              
           ELSE
              WRITE(numout,*) 'ERROR: turnover has not been defined - 3'
              CALL ipslerr_p(3,'stomate_turnover',&
                   'turnover: senescence type not known','case 3','')
              
           ENDIF ! is_tree(ivm)
           
       ENDDO ! # leaf ages

       !! 4.3 Recalculate the fraction of leaf biomass in each leaf age class.
       !  Older leaves will be dropped more than younger leaves and therefore 
       !  the leaf age distribution needs to be recalculated after turnover. 
       !  The fraction of biomass in each leaf class is updated using the
       !  following equation:
       !  \latexonly
       !  \input{turnover_update_LeafAgeDistribution_eqn6.tex}
       !  \endlatexonly
       !  \n
       !
       !  new fraction = new leaf mass of that fraction / new total leaf mass
       !               = (old fraction*old total leaf mass ::lm_old(:) + &
       !                  biomass change of that fraction ::delta_lm(:,ipar))/&
       !                  new total leaf mass ::biomass(:,ivm,ileaf)
       DO ilage = 1, nleafages
          WHERE (SUM(circ_class_biomass(:,ivm,:,ileaf,icarbon)*circ_class_n(:,ivm,:),2) .GT. min_stomate )
             leaf_frac(:,ivm,ilage) = (leaf_frac(:,ivm,ilage)*lm_old(:)+ & 
                delta_lm(:,ilage)) / SUM(circ_class_biomass(:,ivm,:,ileaf,icarbon)* &
                circ_class_n(:,ivm,:),2)
          ELSEWHERE
             leaf_frac(:,ivm,ilage) = zero
          ENDWHERE
       ENDDO
    ENDDO         ! loop over PFTs
      
   ! Debug
    IF (printlev_loc>=4) THEN
       DO ipar = 1,nparts
          DO iele = 1,nelements
            DO icir =1, ncirc
               IF(icir == 1 .AND. iele == 1)&
               WRITE(numout,*) 'Biomass/turnover check 04: ',&
                    circ_class_biomass(test_grid,test_pft,1,ipar,iele) * &
                    circ_class_n(test_grid,test_pft,1),&
                    turnover(test_grid,test_pft,ipar,iele)
            ENDDO
         ENDDO
      ENDDO
    ENDIF

    !! 5. Drop all leaves if there is a very low leaf mass during senescence
    
    !  Both for deciduous trees and grasses same conditions are checked:
    !  If biomass is large enough (i.e. when it is greater than zero), 
    !  AND when senescence is set to true
    !  AND the leaf biomass drops below a critical minimum biomass level
    !  (i.e. when it is lower than half the minimum initial LAI
    !  ::lai_initmin(ivm) divided by the specific leaf area ::slainit(ivm),
    !  ::lai_initmin(ivm) is set to 0.3 in stomate_data.f90 and slainit is a
    !  PFT-specific constant, If these conditions are met,
    !  the flag ::shed_rest(:) is set to TRUE.
    !
    !  After this, the biomass of different carbon pools both for trees and
    !  grasses is set to zero and the mean leaf age is reset to zero. Finally,
    !  the leaf fraction and leaf age of the different leaf age classes is set
    !  to zero.
    DO ivm = 2,nvm ! Loop over # PFTs
       ! Debug
       IF (ivm .EQ. test_pft .AND. printlev_loc >= 4) THEN
          WRITE(numout,*) 'low leaf mass, '
       ENDIF
       shed_rest(:) = .FALSE.
       !! 5.1 For deciduous trees: leaves, fruits and fine roots are dropped 
       !  For deciduous trees: next to leaves, also fruits and fine roots
       !  are dropped: fruit ::biomass(:,ivm,ifruit) and fine root
       !  ::biomass(:,ivm,iroot) carbon pools are set to zero.
       IF ( is_tree(ivm) .AND. ( senescence_type(ivm) .NE. 'none' ) ) THEN
          ! Check whether we shed the remaining leaves. The condition
          ! depends on biomass(ileaf) so first calculate the sheding
          ! at the tree level. because biomass(ileaf) is not changed
          ! at the tree level the same statement can then be used at
          ! the stand level.
          ! slainit can be used because a threshold is calculated. The
          ! leaf mass calculated by lai_initmin/2.)/slainit -see below
          ! is not used in any other calculations.
          init_biomass(:) = SUM(circ_class_biomass(:,ivm,:,ileaf,icarbon)*&
               circ_class_n(:,ivm,:),2)

          DO iele = 1,nelements
            DO ij =1,3 !ileaf, iroot,ifruit
                ipar=parts(ij)
                IF (ipar .EQ. ileaf .AND. iele .EQ. initrogen)THEN
                    recycle=recycle_leaf(ivm)
                ELSEIF (ipar .EQ. iroot .AND. iele .EQ. initrogen ) THEN
                    recycle=recycle_root(ivm)
                ELSE
                    recycle=zero
                ENDIF
                dturnover(:,:) = zero
                DO icir=1, ncirc
                   WHERE ( ( init_biomass(:) .GT. zero ) .AND.&
                        plant_status(:,ivm) .EQ. isenescent .AND. &
                        ( init_biomass(:) .LT. lai_to_biomass(lai_initmin(ivm),ivm)/2))
                      dturnover(:,icir) = circ_class_biomass(:,ivm,icir,ipar,iele)
                      circ_class_biomass(:,ivm,icir,ipar,iele)= zero
                      ! Account for resorption at the tree level (gC tree-1)
                      ! The last day of senescence, the recycled nitrogen is moved into 
                      ! the reserve pool rather than the labile pool (which is done all 
                      ! the other days during senescence). The reason is that the first 
                      ! day of dormancy, there are no leaves left and then the labile 
                      ! pool is moved into the reserve pool. Deciduous species showed an
                      ! unrealistic one-day spike in the labile pool. By moving the N directly
                      ! into the reserve pool, this spike is avoided.
                      circ_class_biomass(:,ivm,icir,icarbres,iele) = &
                           circ_class_biomass(:,ivm,icir,icarbres,iele)+ &
                           recycle * dturnover(:,icir)
                      turnover(:,ivm,ipar,iele) = turnover(:,ivm,ipar,iele)+&
                           dturnover(:,icir)*circ_class_n(:,ivm,icir)*(un-recycle)
                   ENDWHERE
                ENDDO
             ENDDO
          ENDDO

          ! slainit can be used because a threshold is calculated. The
          ! leaf mass calculated by lai_initmin/2.)/slainit -see below
          ! is not used in any other calculations.
          WHERE(( init_biomass(:) .GT. zero ) .AND. &
               plant_status(:,ivm) .EQ. isenescent .AND. &
               ( init_biomass(:) .LT. &
               lai_to_biomass(lai_initmin(ivm),ivm)/2) )
             plant_status(:,ivm) = idormant
             leaf_meanage(:,ivm) = zero
             shed_rest(:) = .TRUE.
          ENDWHERE
          
       ELSEIF (.NOT. is_tree(ivm)) THEN

          IF (senescence_type(ivm) .EQ. 'crop') THEN
          
             ! Nothing should be done because we now
             ! have a harvest module for crops.
             
          ELSEIF (senescence_type(ivm) .EQ. 'mixed') THEN  

             !! 6.2 Drop leaves, roots, fruit and sapwood for grasses
             !  For grasses all aboveground carbon pools, except the
             !  carbohydrate reserves are affected: fruit
             !  ::biomass(:,ivm,ifruit,icarbon), fine root
             !  ::biomass(:,ivm,iroot,icarbon) and sapwood above 
             !  ::biomass(:,ivm,isapabove,icarbon) carbon pools are set
             !  to zero.
             !  Shed the remaining leaves if LAI very low.
             init_biomass(:) = circ_class_biomass(:,ivm,1,ileaf,icarbon)*&
                  circ_class_n(:,ivm,1)
             DO iele = 1,nelements
                ! Grasses should live on after senescence. Do not
                ! shed the sapwood because then the whole 
                ! system collapse because the reserves are
                ! calculated based on the sapwood. no sapwood = no 
                ! reserves = no phenology = no growth.
                !+++ CHECK +++
                ! The trunk seems to take part of the sapwood for grasses
                ! in senescence, but we take none here.
                DO ij =1,3 !ileaf, iroot, ifruit
                    ipar=parts(ij)
                    IF (ipar .EQ. ileaf .AND. iele .EQ. initrogen) THEN
                        recycle=recycle_leaf(ivm)
                    ELSEIF (ipar .EQ. iroot .AND. iele .EQ. initrogen) THEN
                        recycle=recycle_root(ivm)
                    ELSE
                        recycle=zero
                    ENDIF
                    dturnover(:,:) = zero
                    ! slainit can be used because a threshold is calculated. The
                    ! leaf mass calculated by lai_initmin/2.)/slainit -see below
                    ! is not used in any other calculations.
                    WHERE ( ( init_biomass(:) .GT. min_stomate ) .AND.&
                         plant_status(:,ivm) .EQ. isenescent .AND. &
                         ( init_biomass(:) .LT. &
                         lai_to_biomass(lai_initmin(ivm),ivm)/2) )
                       ! Account for resorption at the tree level (gC tree-1)
                       ! This should be out of the iele loop else labile is
                       ! calculated twice and the mass balance is no longer preserved.
                       dturnover(:,1) = circ_class_biomass(:,ivm,1,ipar,iele)
                       circ_class_biomass(:,ivm,1,ipar,iele)= zero
                       circ_class_biomass(:,ivm,1,ilabile,iele) = &
                            circ_class_biomass(:,ivm,1,ilabile,iele)+ &
                            recycle * dturnover(:,1)
                       turnover(:,ivm,ipar,iele) = turnover(:,ivm,ipar,iele)+ &
                            dturnover(:,1)*circ_class_n(:,ivm,1)*(un-recycle)
                    ENDWHERE
                 ENDDO
            ENDDO
            ! slainit can be used because a threshold is calculated. The
            ! leaf mass calculated by lai_initmin/2.)/slainit -see below
            ! is not used in any other calculations.
            WHERE( ( init_biomass(:) .GT. min_stomate ) .AND.&
                 plant_status(:,ivm) .EQ. isenescent .AND. &
                 ( init_biomass(:) .LT.&
                 lai_to_biomass(lai_initmin(ivm),ivm)/2) )
               shed_rest(:) = .TRUE.
               plant_status(:,ivm) = idormant
               leaf_meanage(:,ivm) = zero
            ENDWHERE

          ELSEIF (senescence_type(ivm) .EQ. 'none') THEN
             ! Nothing should be done
             
          ELSE
             WRITE(numout,*) 'ERROR: senescence type is not known - 4'
             CALL ipslerr_p(3,'stomate_turnover',&
                  'turnover: senescence type not known','case 4','')
          
          ENDIF
          
          ! Debug
          IF (ivm .EQ. test_pft .AND. printlev_loc >= 4) THEN
             WRITE(numout,*) 'phenology, shed_rest, ',shed_rest(:)
          ENDIF
          !-
       
       ENDIF ! is_tree(ivm)

 
       !! Kill PFTs that have no chance to survive their dormancy
       ! If the plant goes into dormancy the leaves and roots will be
       ! shed. If the labile and reserve C are empty, the plant won't be 
       ! able to grow again so it should be killed. In most cases this
       ! will be dealt with in stomate_mark_to_kill but we had a couple of 
       ! pixels that went to idormant at time t and back to ibudavail at 
       ! t+1. Hence the model never passed mortality and crashed in phenology
       ! because there were no reserve and/or labile carbon pools to 
       ! grow the initial leaves. If a PFT is killed we will set the
       ! when_growthinit to zero so it is no longer possible to go to 
       ! ibudavail in phenology.
       DO ipts = 1,npts

          IF(SUM((circ_class_biomass(ipts,ivm,:,ileaf,icarbon) + &
               circ_class_biomass(ipts,ivm,:,ilabile,icarbon) + & 
               circ_class_biomass(ipts,ivm,:,icarbres,icarbon)) * & 
               circ_class_n(ipts,ivm,:)).LT.min_stomate .AND. &
               plant_status(ipts,ivm).EQ.idormant) THEN

             ! Make sure the plant won't go into budavail or budbreak 
             ! the next day. This PFT should get killed.
             when_growthinit(ipts,ivm) = zero

          END IF

       END DO

       !! 5.3 Reset the leaf age structure
       !  The leaf fraction and leaf age of the different leaf
       !  age classes is set to zero.
       DO ilage = 1, nleafages
          WHERE ( shed_rest(:) )
             leaf_age(:,ivm,ilage) = zero
             leaf_frac(:,ivm,ilage) = zero
          ENDWHERE
       ENDDO
    ENDDO  ! loop over PFTs
    
    ! Debug
    IF (printlev_loc>=4) THEN
       DO ipar = 1,nparts
          DO iele = 1,nelements
             DO icir =1, ncirc
                IF(icir == 1 .AND. iele == 1)&
                     WRITE(numout,*) 'Biomass/turnover check  05: ',&
                     circ_class_biomass(test_grid,test_pft,1,ipar,iele) * &
                     circ_class_n(test_grid,test_pft,1),&
                     turnover(test_grid,test_pft,ipar,iele)
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    
    !! 6. Herbivore activity: elephants, cows, gazelles but no lions.
    
    !  Herbivore activity affects the biomass of leaves and fruits as well 
    !  as stalks (only for grasses). Herbivore activity does not modify leaf 
    !  age structure. Herbivores ::herbivores(:,ivm) is the time constant of 
    !  probability of a leaf to be eaten by a herbivore, and is calculated in 
    !  ::stomate_season. following Mc Naughton et al. [1989].
    IF ( ok_herbivores ) THEN

       ! If the herbivore activity is allowed (if ::ok_herbivores is true,
       ! which is set in run.def), remove the amount of biomass consumed
       ! by herbivory from the leaf biomass ::biomass(:,ivm,ileaf,icarbon) and 
       ! the fruit biomass ::biomass(:,ivm,ifruit,icarbon). The daily amount
       ! consumed equals the biomass multiplied by 1 day divided by the time
       ! constant ::herbivores(:,ivm).
       DO ivm = 2,nvm ! Loop over # PFTs
         IF ( is_tree(ivm) ) THEN
             ! Debug
             IF (ivm .EQ. test_pft .AND. printlev_loc >= 4) THEN
                WRITE(numout,*) 'herbivores, ',herbivores
             ENDIF
             !-
             ! For trees: only the leaves and fruit carbon pools are affected
             ! No time to re-allocate nitrogen during browsing
             !+++CHECK+++
             ! Why is root herbivory not accounted for. Could be difficult
             ! to parameterize but if leaf browsing is accounted for, there
             ! is no reason to ignore root browsing
             init_biomass(:) = SUM(circ_class_biomass(:,ivm,:,ileaf,icarbon)*&
                circ_class_n(:,ivm,:),2)
             DO iele=1, nelements
                ! Turnover (gC tree-1)
                DO ij =2,3 ! ileaf and ifruit
                    ipar=parts(ij)
                    dturnover(:,:) = zero
                    DO icir=1, ncirc
                        WHERE (init_biomass(:) .GT. zero .AND.herbivores(:,ivm) .GT. min_stomate)
                            ! Tree Circ level turnover (gC tree-1):
                            dturnover(:,icir) = circ_class_biomass(:,ivm,icir,ipar,iele) * &
                                dt/herbivores(:,ivm)
                        ENDWHERE
                    ENDDO
                    !the stand level turnover (gC m-2):
                    turnover(:,ivm,ipar,iele) = turnover(:,ivm,ipar,iele)+ &
                        SUM(dturnover(:,:)*circ_class_n(:,ivm,:),2)
                    ! Update circ_class_biomass
                    circ_class_biomass(:,ivm,:,ipar,iele) = &
                        circ_class_biomass(:,ivm,:,ipar,iele)-dturnover(:,:)
                ENDDO
             ENDDO

          ELSEIF (.NOT. is_tree(ivm)) THEN 

             init_biomass(:) = circ_class_biomass(:,ivm,1,ileaf,icarbon)*&
                circ_class_n(:,ivm,1)
             DO iele=1, nelements
                ! Turnover (gC tree-1)
                DO ij =2,4 !leaf, ifruit, isapwood
                    ipar=parts(ij)
                    dturnover(:,:) = zero
                    WHERE (init_biomass(:) .GT. zero .AND. herbivores(:,ivm) .GT. min_stomate)
                        ! Tree Circ level turnover (gC tree-1):
                        dturnover(:,1) = circ_class_biomass(:,ivm,1,ipar,iele)* &
                            dt/herbivores(:,ivm)
                    ENDWHERE
                    !the stand level turnover (gC m-2):
                    turnover(:,ivm,ipar,iele) = turnover(:,ivm,ipar,iele)+&
                        SUM(dturnover(:,:)*circ_class_n(:,ivm,:),2)
                    ! Update circ_class_biomass
                    circ_class_biomass(:,ivm,1,ipar,iele) = &
                        circ_class_biomass(:,ivm,1,ipar,iele)-dturnover(:,1)
                ENDDO
             ENDDO
             
          ELSE
             WRITE(numout,*) 'ERROR: vegetation type is not known'
             CALL ipslerr_p(3,'stomate_turnover',&
                  'vegetation type is not known','','')
             
          ENDIF  ! tree/grass?
       ENDDO    ! loop over PFT
    ENDIF ! end herbivores

    ! Debug
    IF (printlev_loc>=4) THEN
       DO ipar = 1,nparts
          DO iele = 1,nelements
             DO icir =1, ncirc
               IF(icir == 1 .AND. iele == 1)&
               WRITE(numout,*) 'Biomass/turnover check 06: ',&
                    circ_class_biomass(test_grid,test_pft,1,ipar,iele) * &
                    circ_class_n(test_grid,test_pft,1),&
                    turnover(test_grid,test_pft,ipar,iele)
            ENDDO
         ENDDO
      ENDDO
    ENDIF
    
!! 8. Conversion of sapwood to heartwood

    !  Conversion of sapwood to heartwood both for aboveground and
    !  belowground sapwood and heartwood. Following LPJ (Sitch et al., 2003),
    !  sapwood biomass is converted into heartwood biomass with a time
    !  constant tau ::longevity_sap(ivm) in years. Note that this biomass conversion
    !  is not added to "turnover" as the biomass is not lost!
    DO ivm = 2,nvm ! Loop over # PFTsi
       IF ( is_tree(ivm) ) THEN
          ! For the recalculation of age in 9.2 (in case the vegetation is
          ! not dynamic ie. ::ok_dgvm is FALSE), the heartwood above and
          ! below is stored in ::hw_old(:).
          IF ( .NOT. ok_dgvm ) THEN
             hw_old(:) = SUM(circ_class_biomass(:,ivm,:,iheartabove,icarbon),2) + &
                  SUM(circ_class_biomass(:,ivm,:,iheartbelow,icarbon),2)
             sw_old(:) = SUM(circ_class_biomass(:,ivm,:,isapabove,icarbon),2) + &
                  SUM(circ_class_biomass(:,ivm,:,isapbelow,icarbon),2) 
          ENDIF

          !! 8.1 Calculate the rate of sapwood to heartwood conversion 
          !  Calculate the rate of sapwood to heartwood conversion with
          !  the time constant ::longevity_sap(ivm) and update aboveground and
          !  belowground sapwood ::biomass(:,ivm,isapabove) and
          !  ::biomass(:,ivm,isapbelow) and heartwood
          !  ::biomass(:,ivm,iheartabove) and ::biomass(:,ivm,iheartbelow).
            ! Tree level above and belowground (gC tree-1)
          DO iele = 1,nelements
            DO icir = 1,ncirc
                circ_class_biomass(:,ivm,icir,iheartabove,iele) =  &
                    circ_class_biomass(:,ivm,icir,iheartabove,iele) + &
                    circ_class_biomass(:,ivm,icir,isapabove,iele) * dt / &
                    longevity_eff_sap(:,ivm)
                circ_class_biomass(:,ivm, icir,isapabove,iele) = &
                    circ_class_biomass(:,ivm, icir,isapabove,iele) * &
                    (un - dt / longevity_eff_sap(:,ivm)) 
                circ_class_biomass(:,ivm,icir,iheartbelow,iele) = &
                    circ_class_biomass(:,ivm, icir,iheartbelow,iele) + &
                    circ_class_biomass(:,ivm, icir,isapbelow,iele) * dt / &
                    longevity_eff_sap(:,ivm)
                circ_class_biomass(:,ivm, icir,isapbelow,iele) = &
                    circ_class_biomass(:,ivm, icir,isapbelow,iele) * &
                    (un - dt / longevity_eff_sap(:,ivm))
            ENDDO
          ENDDO
          !! 8.2 If the vegetation is not dynamic, the age of the plant
          !!     is decreased. 
          !  The updated heartwood, the sum of new heartwood above and new
          !  heartwood below after converting sapwood to heartwood, is saved
          !  as ::hw_new(:). Creation of new heartwood decreases the age of
          !  the plant with a factor that is determined by: old heartwood
          !  ::hw_old(:) divided by the new heartwood ::hw_new(:)
          IF ( .NOT. ok_dgvm ) THEN
             hw_new(:) = SUM(circ_class_biomass(:,ivm,:,iheartabove,icarbon),2) + &
                  SUM(circ_class_biomass(:,ivm,:,iheartbelow,icarbon),2)
             sw_new(:) = SUM(circ_class_biomass(:,ivm,:,isapabove,icarbon),2) + &
                  SUM(circ_class_biomass(:,ivm,:,isapbelow,icarbon),2)

             WHERE ( hw_new(:) .GT. min_stomate )
                age(:,ivm) = age(:,ivm) * hw_old(:)/hw_new(:)
             ENDWHERE
          ENDIF
       ENDIF
    ENDDO       ! loop over PFTs
       
    
    ! Write to output file
    CALL xios_orchidee_send_field("HERBIVORES",herbivores)
    
    ! the xios operator for leaf_age and leaf_age_crit are maximum so 
    ! the zero values are correctly accounted for in the history file
    CALL xios_orchidee_send_field("LEAF_AGE",leaf_meanage)
    CALL xios_orchidee_send_field("LEAF_AGE_CRIT",leaf_age_crit)
    
    ! By dividing the LEAF_M_MAX_c by LEAF_TURN_AGEING_c one gets
    ! the annual turnover of leaves during prior to senescence. Based on
    ! litter traps this value should be 10 to 20%. 
    CALL xios_orchidee_send_field("LEAF_TURN_AGEING_c",leaf_turn_ageing)

    ! End of growing season. For phenology ibudbreak indicates the start
    ! of the growing season. There is no such single-day setting for senescence
    ! and all stati lasts several days. Use a change in plant_status
    ! to identify the single-day at which the leaf season has ended.
    DO ivm = 2,nvm
       DO ipts = 1,npts
          ! Deciduous trees and grasses should go through isenescent. By
          ! using .GE.idormant also idead is included.
          IF(natural(ivm) .AND. &
               last_plant_status(ipts,ivm) .EQ. isenescent .AND. &
               plant_status(ipts,ivm) .GE. idormant) THEN       
             doy_end_gs(ipts,ivm) = julian_diff
          ELSEIF (.NOT. natural(ivm) .AND. &
               (last_plant_status(ipts,ivm) .EQ. icanopy .OR. &
                last_plant_status(ipts,ivm) .EQ. ipresenescence ) .AND. &
                plant_status(ipts,ivm) .GE. idormant) THEN       
             doy_end_gs(ipts,ivm) = julian_diff
          ENDIF
       ENDDO
    ENDDO ! loop of nvm

    CALL xios_orchidee_send_field("DOY_END_GS",doy_end_gs)
    CALL xios_orchidee_send_field("DOY_ISENE",doy_isenescent)

    CALL histwrite_p (hist_id_stomate, 'LEAF_AGE', itime, &
         leaf_meanage, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'HERBIVORES', itime, &
         herbivores, npts*nvm, horipft_index)

    

 !! 9. Check numerical consistency of this routine

    ! Debug
    IF (printlev_loc>=4) THEN
       DO ipar = 1,nparts
          DO iele = 1,nelements
            DO icir =1, ncirc
               IF(icir == 1 .AND. iele == 1)&
               WRITE(numout,*) 'Biomass/turnover check  08: ',&
                    circ_class_biomass(test_grid,test_pft,1,ipar,iele) * &
                    circ_class_n(test_grid,test_pft,1),&
                    turnover(test_grid,test_pft,ipar,iele)
            ENDDO
         ENDDO
      ENDDO
    ENDIF

    IF (err_act.GT.1) THEN

       ! 9.1 Check surface area
       CALL check_vegetation_area("stomate_turnover", npts, veget_max_begin, &
            veget_max,'pft')

       ! 9.2 Calculate final biomass
       pool_end(:,:,:) = zero
       DO ipar = 1,nparts
          DO iele = 1,nelements
            DO icir =1, ncirc
                pool_end(:,:,iele) = pool_end(:,:,iele) + &
                    (circ_class_biomass(:,:,icir,ipar,iele) * &
                    circ_class_n(:,:,icir)* veget_max(:,:))
            ENDDO
            pool_end(:,:,iele) = pool_end(:,:,iele) + &
                (turnover(:,:,ipar,iele) * veget_max(:,:))

          ENDDO
       ENDDO

       ! The biomass harvest pool is expressed in gC pixel-1 So, it 
       ! shouldn't be multiplied by veget_max but it should be divided
       ! by area to obtain gC m-2.
       DO ivm = 1,nvm
          DO iele = 1,nelements
             pool_end(:,ivm,iele) = pool_end(:,ivm,iele) + &
                  SUM(SUM(harvest_pool(:,ivm,:,iele,:),3),2) / &
                  (area(:) * contfrac(:))
          ENDDO
       ENDDO

       !! 9.3 Calculate components of the mass balance
       check_intern(:,:,iatm2land,:) = zero
       check_intern(:,:,iland2atm,:) = -un * zero
       check_intern(:,:,ilat2out,:) = zero
       check_intern(:,:,ilat2in,:) = -un * zero
       check_intern(:,:,ipoolchange,:) = &
            un * (pool_end(:,:,:) - pool_start(:,:,:))
       closure_intern(:,:,:) = zero
       DO imbc = 1,nmbcomp
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
       
       ! 9.4 Check mass balance closure
       CALL check_mass_balance("stomate_turnover", closure_intern, npts, &
            pool_end, pool_start, veget_max,'pft')
      
    ENDIF ! err_act.GT.1
    
    IF(printlev>=3) WRITE(numout,*) 'Leaving turnover'

  END SUBROUTINE turn_over


!================================================================================================================================
!! SUBROUTINE    : drought_mortality
!!
!>\BRIEF         Calculate turnover of sapwood to heartwood induced by drought.
!!
!! DESCRIPTION : This subroutine determines the turnover of sapwood into
!! heartwood induced by drought.
!! RECENT CHANGE(S) : None.
!!
!! MAIN OUTPUT VARIABLES: ::circ_class_biomass
!!
!! REFERENCE(S) :
!!
!! FLOWCHART    :
!_
!================================================================================================================================

 SUBROUTINE drought_mortality (npts, kill_vessels, vessel_mortality_daily, &
      veget_max, biomass_init_drought, bm_to_litter, circ_class_biomass, &
      circ_class_n)

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                       :: npts                           !! Domain size (number of grid cells).
    LOGICAL, DIMENSION(:,:,:),   INTENT(in)          :: kill_vessels                   !! Flag to kill vessels at the end of the day following embolism.
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)        :: vessel_mortality_daily         !! Proportion of daily vessel mortality due to cavitation in the xylem. Equal to 
                                                                                       !! a fraction of vl_diff.
    REAL(r_std), DIMENSION(:,:), INTENT(in)          :: veget_max                      !! 'Maximal' coverage fraction of a PFT (LAI -> Infinity) on ground (unitless).

    !! 0.2 Output variables

    !! 0.3 Modified variables
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout) :: biomass_init_drought           !! Biomass of heartwood or sapwood before onset of drought. Used to compute
                                                                                       !! turnover on same reference biomass in stomate_turnover.f90. Should be the same
                                                                                       !! along entire drought episode.
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)   :: bm_to_litter                   !! Background mortality of biomass (not senescence-driven).

    !! 0.3 Modified variables
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout) :: circ_class_biomass             !! Biomass of the components of the model tree within a circumference class
                                                                                       !! (gC/ind.).
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)     :: circ_class_n                   !! Number of individuals in each circ class.

    !! 0.4 Local variables
    INTEGER(i_std)                                   :: ivm, iele, ilage, ipts         !! Index (unitless).
    INTEGER(i_std)                                   :: ipar, icir, imbc, ij           !! Index (unitless).
    REAL(r_std)                                      :: total_sap_biomass              !! Total sapwood biomass of model tree.
    REAL(r_std), DIMENSION(nparts)                   :: vessel_turnover                !! Turnover for sapwood of model tree.
    REAL(r_std), DIMENSION(npts,nvm,ncirc)           :: total_cc_biomass               !! Total wood biomass of model tree.
    REAL(r_std), DIMENSION(npts)                     :: hw_old                         !! Old heartwood mass (gC/m2).
    REAL(r_std), DIMENSION(npts)                     :: hw_new                         !! New heartwood mass (gC/m2).
    REAL(r_std), DIMENSION(npts)                     :: sw_new                         !! New sapwood mass (gC/m2 of nat./agri. ground).
    REAL(r_std), DIMENSION(npts,nvm,nmbcomp,nelements) &
                                                     :: check_intern                   !! Contains the components of the internal mass balance check for this routine
                                                                                       !!(gC/pix./dt).
    REAL(r_std), DIMENSION(npts,nvm,nelements)       :: closure_intern                 !! Check closure of internal mass balance (gC/pix./dt).
    REAL(r_std), DIMENSION(npts,nvm,nelements)       :: pool_start, pool_end           !! Start and end pool of this routine gC/pix./dt).
    REAL(r_std), DIMENSION(npts,nvm)                 :: veget_max_begin                !! Temporary storage of veget_max to check area vegetation (unitless, [0-1]).

!================================================================================================================================

    !! 1.1. Initialize check for mass balance closure

    ! The mass balance is calculated at the end of this routine
    ! in section 3.
    IF(err_act .GT. 1) THEN
       pool_start = zero
       DO iele = 1,nelements
          ! Biomass pool + bm_to_litter.
          DO ipar = 1,nparts
             DO icir = 1,ncirc
                ! Initial biomass pool.
                pool_start(:,:,iele) = pool_start(:,:,iele) + &
                     circ_class_biomass(:,:,icir,ipar,iele) * &
                     circ_class_n(:,:,icir) * veget_max(:,:)
             ENDDO
             ! Add turnover to the initial biomass pool.
             pool_start(:,:,iele) = pool_start(:,:,iele) + &
                  (bm_to_litter(:,:,ipar,iele) * veget_max(:,:))
          ENDDO
       ENDDO
       !! 1.2. Initialize check for area conservation
       veget_max_begin(:,:) = veget_max(:,:)

    ENDIF ! IF(err_act.GT.1)

    ! 2. Calculate the rate of sapwood to heartwood conversion with
    !    the rate of daily mortality ::vessel_mortality_daily(:,ivm).
    DO ivm = 2,nvm

       IF(is_tree(ivm)) THEN

          ! Calculation of effect of embolism on heartwood and sapwood
          ! biomass. Embolism turns sapwood into heartwood. Every time,
          ! fraction of sapwood biomass is subtracted from sapwood biomass
          ! and added to heartwood biomass. Distinction between aboveground
          ! and belowground biomass.
          DO ipts = 1,npts
             DO iele = 1,nelements
                DO icir = 1,ncirc
                   
                   ! Total sapwood in gC tree-1
                   total_sap_biomass = &
                        circ_class_biomass(ipts,ivm,icir,isapabove,iele) + &
                        circ_class_biomass(ipts,ivm,icir,isapbelow,iele)
                      
                   ! Part for drought: If the flag ::kill_vessels(:,ivm) is TRUE,
                   ! and there is sapwood remaining, which means the tree is still
                   ! alive, we compute the turnover induced by drought.
                   IF(kill_vessels(ipts,ivm,icir) .AND. &
                        total_sap_biomass .GT. min_stomate) THEN
                      
                      ! First, compute the mass of sapwood that will be
                      ! subtracted from the model tree. Note that vessel_mortality_daily
                      ! contains the increase in mortality (thus not the entire mortality)
                      ! It is multiplied with the biomass at the start of the drought
                      ! to avoid double counting.
                      vessel_turnover(:) = biomass_init_drought(ipts,ivm,icir,:,iele) * &
                           vessel_mortality_daily(ipts,ivm,icir)

                      ! Compare ::vessel_turnover(ipar) to the current sapwood
                      ! biomass of the model tree to make sure we do not end up
                      ! with a negative sapwood biomass, which is not realistic. If
                      ! ::vessel_turnover(ipar) is greater than current sapwood
                      ! biomass, we attribute to ::vessel_turnover(ipar) the value
                      ! of the current sapwood biomass, so we will get a null
                      ! sapwood biomass at the end of the turnover calculation.
                      ! Aboveground wood biomass.
                      IF(vessel_turnover(isapabove) .GE. &
                           circ_class_biomass(ipts,ivm,icir,isapabove,iele)) THEN
                         vessel_turnover(isapabove) = &
                              circ_class_biomass(ipts,ivm,icir,isapabove,iele)
                      ENDIF

                      ! Use biomass_init_drought in the calculation of
                      ! turnover induced by drought so mortality is not
                      ! overestimated, since ::vessel_mortality_daily(:,ivm) is
                      ! calculated as a proportion of sapwood biomass.
                      ! Biomass of aboveground heartwood. We add dead sapwood to
                      ! heartwood biomass.
                      circ_class_biomass(ipts,ivm,icir,iheartabove,iele) = &
                           circ_class_biomass(ipts,ivm,icir,iheartabove,iele) + &
                           vessel_turnover(isapabove)

                      ! Biomass of aboveground sapwood. We subtract dead sapwood to
                      ! sapwood biomass.
                      circ_class_biomass(ipts,ivm, icir,isapabove,iele) = &
                           circ_class_biomass(ipts,ivm,icir,isapabove,iele) - &
                           vessel_turnover(isapabove)

                      ! Belowground wood biomass.
                      IF(vessel_turnover(isapbelow) .GE. &
                           circ_class_biomass(ipts,ivm,icir,isapbelow,iele)) THEN
                         vessel_turnover(isapbelow) = &
                              circ_class_biomass(ipts,ivm,icir,isapbelow,iele)
                      ENDIF

                      ! Biomass of belowground heartwood. We add dead sapwood to
                      ! heartwood biomass.
                      circ_class_biomass(ipts,ivm,icir,iheartbelow,iele) = &
                           circ_class_biomass(ipts,ivm, icir,iheartbelow,iele) + &
                           vessel_turnover(isapbelow)

                      ! Biomass of belowground sapwood. We subtract dead sapwood to
                      ! sapwood biomass.
                      circ_class_biomass(ipts,ivm, icir,isapbelow,iele) = &
                           circ_class_biomass(ipts,ivm,icir,isapbelow,iele) - &
                           vessel_turnover(isapbelow)

                   ENDIF ! IF(kill_vessels)

                   ! Update total sapwood biomass of model tree after vessel
                   ! moratlity has been dealt with
                   total_sap_biomass = &
                        circ_class_biomass(ipts,ivm,icir,isapabove,iele) + &
                        circ_class_biomass(ipts,ivm,icir,isapbelow,iele)

                   ! Where sapwood biomass becomes zero, model must process
                   ! tree death. The model is still in a DO-loop over iele
                   ! hence the code below will kill the whole tree is either
                   ! icarbon or initrogen is empty.
                   IF(total_sap_biomass .LE. min_stomate) THEN

                      ! Current biomass of model tree is moved into litter.
                      ! Carbon biomass. We are still in a do-loop over iele
                      ! so we cannot have another do-loop over iele. It was
                      ! written explicitly.
                      bm_to_litter(ipts,ivm,:,icarbon) = &
                           bm_to_litter(ipts,ivm,:,icarbon) + &
                           circ_class_biomass(ipts,ivm,icir,:,icarbon) * &
                           circ_class_n(ipts,ivm,icir)

                      ! Nitrogen biomass.
                      bm_to_litter(ipts,ivm,:,initrogen) = &
                           bm_to_litter(ipts,ivm,:,initrogen) + &
                           circ_class_biomass(ipts,ivm,icir,:,initrogen) * &
                           circ_class_n(ipts,ivm,icir)

                      ! To get a closed mass balance, we set biomass and number of
                      ! individuals of dead circumference class to zero since the
                      ! biomass has been moved into litter. We also set the initial 
                      ! biomass to zero so that the tree does not revive.
                      circ_class_biomass(ipts,ivm,icir,:,:) = zero
                      circ_class_n(ipts,ivm,icir) = zero
                      biomass_init_drought(ipts,ivm,icir,:,:) = zero

                   ENDIF ! IF(total_sap_biomass .LE. min_stomate)

                ENDDO ! DO icir = 1,ncirc

             ENDDO ! DO iele = 1,nelements

          ENDDO ! DO ipts = 1,kjpindex

       ENDIF ! is_tree

    ENDDO ! DO ivm = 2,nvm
    

    !! 3. Check numerical consistency of this routine
    IF (err_act.GT.1) THEN

       ! 3.1. Check surface area
       CALL check_vegetation_area("drought_mortality", npts, veget_max_begin, &
            veget_max,'pft')

       ! 3.2. Calculate final biomass
       pool_end(:,:,:) = zero
       DO ipar = 1,nparts
          DO iele = 1,nelements
             DO icir =1, ncirc
                pool_end(:,:,iele) = pool_end(:,:,iele) + &
                     (circ_class_biomass(:,:,icir,ipar,iele) * &
                     circ_class_n(:,:,icir)* veget_max(:,:))
             ENDDO
             ! Added ::bm_to_litter(:,ivm,ipar,iele) here.
             pool_end(:,:,iele) = pool_end(:,:,iele) + &
                  bm_to_litter(:,:,ipar,iele) * veget_max(:,:)
          ENDDO
       ENDDO

       ! 3.3. Calculate components of the mass balance
       check_intern(:,:,iatm2land,:) = zero
       check_intern(:,:,iland2atm,:) = -un * zero
       check_intern(:,:,ilat2out,:) = zero
       check_intern(:,:,ilat2in,:) = -un * zero
       check_intern(:,:,ipoolchange,:) = &
            un * (pool_end(:,:,:) - pool_start(:,:,:))
       closure_intern(:,:,:) = zero
       DO imbc = 1,nmbcomp
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
       CALL check_mass_balance("drought_mortality", closure_intern, npts, &
            pool_end, pool_start, veget_max,'pft')

    ENDIF ! IF(err_act .GT. 1)

    IF(printlev>=3) WRITE(numout,*) 'Leaving drought_mortality.'

  END SUBROUTINE drought_mortality

END MODULE stomate_turnover
