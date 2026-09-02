! ============================================================================================================================
! MODULE        : stomate_pest
!
! CONTACT       : orchidee-help _at_ ipsl.listes.ipsl.fr
!
! LICENCE       : IPSL (2006).
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        : Calculates tree damage due to bark-beetle outbreaks 
!!
!!\n DESCRIPTION: See subroutine bark_beetle_damage for a more elaborate description of the
!!              principles.
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCES: - C. TEMPERLI, H. BUGMANN, AND C. ELKIN. Cross-scale interactions
!!             among bark beetles, climate change, and wind disturbances: a
!!             landscape modeling approach. Ecological Monographs, 83(3), 2013, pp. 383–402.
!!             - Marie, G., Jeong, J., Jactel, H., Petter, G., Cailleret, M., McGrath, M., 
!!             Bastrikov, V., Ghattas, J., Guenet, B., Lansø, A.-S., Naudts, K., Valade, A., 
!!             Yue, C., and Luyssaert, S.: Simulating Bark Beetle Outbreak Dynamics and their 
!!             Influence on Carbon Balance Estimates with ORCHIDEE r7791. GMD, 2024. 
!!
!! SVN           :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-DOFOCO/ORCHIDEE/src_stomate/stomate_pest.f90 
!! $Date: 2019-01-22 13:41:13 +0100 (mar., 06 aout. 2019)
!! $Revision: 9352 $
!! \n
!_ ===========================================================================================================================

MODULE stomate_pest

  ! modules used:
  USE ioipsl_para
  USE xios_orchidee
  USE stomate_data
  USE constantes
  USE constantes_soil
  USE pft_parameters
  USE pft_parameters_var
  USE function_library, ONLY: wood_to_ba, wood_to_qmdia, Nmax, wood_to_dia, cc_kill_to_area
  IMPLICIT NONE

  ! private & public routines

  PRIVATE
  PUBLIC pest_clear, bark_beetle_damage, pest_write

  LOGICAL, SAVE                   :: firstcall_bark_beetle_damage = .TRUE.     !! First call 
!$OMP THREADPRIVATE(firstcall_bark_beetle_damage)
  INTEGER(i_std), SAVE            :: printlev_loc                       !! Local level of text output
!$OMP THREADPRIVATE(printlev_loc)

CONTAINS

!!
!==================================================================================================================================
!!  SUBROUTINE   : pest_clear 
!!
!>\BRIEF        Set the firstcall flag to .TRUE. and activate initialization
!! 
!_
!================================================================================================================================

  SUBROUTINE pest_clear
    firstcall_bark_beetle_damage = .TRUE.
  END SUBROUTINE pest_clear


!================================================================================================================================
!! SUBROUTINE    : bark_beetle_damage
!!
!>\BRIEF         Calculates tree damage due to bark-beetle outbreaks 
!!
!!\n DESCRIPTION : This subroutine is used to calculate the annual amount of trees killed by
!!              bark beetles. This module is strongly inspired by the work of Temperli
!!              et al. 2013 in Ecological Monographs. The bark beetle module developped 
!!              for LANDCLIM landscape model was the basis of this work. The module was
!!              adapted to fit the constraints of ORCHIDEE in regard to its spatial and
!!              temporal resolution. In short, this module calculate a suceptibility of
!!              of a forest to be attack by bark beetles at the pixel/PFT level. At 
!!              present the model has been parameterized for Picea abies. The risk
!!              probability of a pixel to be attacked depends on pixel suceptibility and
!!              bark beetle pressure where pressure is calculated as an index that
!!              indicates the activity of the beetles in this pixel. Then, a mortality
!!              rate for each Picea abies PFT is calculated based on PFT suceptibility
!!              and beetle pressure. Finally, the wood biomass kill by bark beetle of a
!!              given PFT is the product of the pixel risk * mortality rate * the PFT
!!              wood biomass. For each PFT we first kill the small trees and if more
!!              biomass needs to be killed, ORCHIDEE will continue with killing the
!!              larger diameter classes.
!!              The most important difference beetween this module and the original one
!!              in LANDCLIM is that ORCHIDEE lost all information about the spatialization
!!              of the attack. In other words the attack is no longer spatialy explicit.
!!  
!! RECENT CHANGE(S): Added in December 2019
!!
!! MAIN OUTPUT VARIABLE(S) : circ_class_kill
!!
!! REFERENCES: CHRISTIAN TEMPERLI, HARALD BUGMANN, AND CHE ELKIN. Cross-scale
!!             interactions among bark beetles, climate change, and wind
!!             disturbances: a landscape modeling approach. Ecological Monographs,
!!             83(3), 2013, pp. 383–402.
!!
!! FLOWCHART     : 
!! \n
!_================================================================================================================================

  SUBROUTINE bark_beetle_damage( &
       npts,  circ_class_biomass, i_beetles_generation,&
       veget_max, circ_class_n, &
       age_stand, season_drought_legacy, resolution, &
       circ_class_kill, forest_managed, &
       wood_leftover_legacy, &
       i_beetles_activity_legacy, P_beetles_attacked_legacy, & 
       B_beetles_kill_legacy)

    IMPLICIT NONE

    !! 0. Variable declarations

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                         :: npts                     !! number of pixels
    INTEGER(i_std), DIMENSION (:, :), INTENT(in)       :: forest_managed           !! forest management flag 
    !! (0 : unmanaged, >0 : managed)
    REAL(r_std), DIMENSION(:, :, :, :, :), INTENT(in)  :: circ_class_biomass       !! Biomass of an individual tree (gC.m-2) 
    !! in each circ class
    REAL(r_std), DIMENSION(:, :, :), INTENT(in)        :: circ_class_n             !! Number of trees(tree.m-2) 
    !! within each circumference
    REAL(r_std), DIMENSION(:, :), INTENT(in)           :: veget_max                !! "maximal" coverage fraction of 
    !! a PFT on the ground
    REAL(r_std), DIMENSION(:, :), INTENT(in)           :: resolution               !! The length in m of each grid-box in 
    !! X and Y axis
    INTEGER(i_std), DIMENSION(:, :), INTENT(in)        :: age_stand                !! Age of stand (years)
    REAL(r_std), DIMENSION(:, :, :), INTENT(in)        :: season_drought_legacy    !! mean growing season moisture availability 
    REAL(r_std), DIMENSION(:, :, :), INTENT(in)        :: wood_leftover_legacy     !! Coarse woody bedris, snags and logs
    !! (gC.m-2) left on the PFT 
    !! during previous timestep
    REAL(r_std), DIMENSION(:, :), INTENT(inout)        :: P_beetles_attacked_legacy!! risk index stored every year, n time 


    !! 0.2 Output variables

    REAL(r_std), DIMENSION(npts, nvm)     :: beetle_damage_wood_vol   !! Wood volume loss from bark beetle @tex $(m^3 m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts, nvm)     :: DR_beetles               !! rate of biomass loss from bark beetle per pft
    REAL(r_std), DIMENSION(npts, nvm)     :: area_damaged_pest        !! Forest area killed by bark beetle, cc-kill basal-area route (m2)

    !! 0.3 Modified variables

    REAL(r_std), DIMENSION(:, :, :, :, :), INTENT(inout):: circ_class_kill         !! Number of trees within a circ that needs
    REAL(r_std), DIMENSION(:, :, :),INTENT(inout)      :: i_beetles_activity_legacy!! Index (0-1) Biomass of tree from the 
    !! same PFT that was
    REAL(r_std), DIMENSION(:, :, :), INTENT(inout)     :: B_beetles_kill_legacy    !! temporary variable(gC.m-2)
    !! infected during the previous timestep  
    REAL(r_std), DIMENSION(npts, nvm)                  :: i_hosts_dead_monitor     !! variable used to monitor windthrow susceptibility (unitless)
    REAL(r_std), DIMENSION(npts, nvm)                  :: i_hosts_attractivity_monitor!! variable used to monitor overall susceptibility index (unitless)
    REAL(r_std), DIMENSION(npts, nvm)                  :: i_hosts_competition_monitor !! variable used to monitor rdi susceptibility (unitless)
    REAL(r_std), DIMENSION(npts, nvm)                  :: i_hosts_alive_monitor
    !! variable used to monitor rdi susceptibility (unitless)
    REAL(r_std), DIMENSION(npts, nvm)        :: i_beetles_massattack_monitor
    REAL(r_std), DIMENSION(npts, nvm)        :: i_beetles_survival_monitor
    !! variable used to monitor rdi susceptibility (unitless)
    REAL(r_std), DIMENSION(npts, nvm)        :: i_beetles_activity_monitor
    REAL(r_std), DIMENSION(:, :, :), INTENT(inout)     :: i_beetles_generation

    !! variable used to monitor rdi susceptibility (unitless)
    REAL(r_std), DIMENSION(npts, nvm)        :: i_hosts_defence_monitor      !! variable used to monitor ability of the trees to defend themselves against bark beetle attacks (unitless)
    REAL(r_std), DIMENSION(npts, nvm)        :: i_beetles_generation_monitor !! variable used to monitor generation index susceptibility

    !! 0.4 Local variables
    INTEGER(i_std)                                    :: ipts, ileg, ivm, icir, &
         ivmap, iagec, iyear, xvm, &
         ncirc_loc
    INTEGER(i_std), DIMENSION(npts, nvmap, nagec)     :: nc                        !! number of PFT counter
    !! within a species group that are
    !! infected by bark beetles (unitless) 
    REAL(r_std), DIMENSION(npts, nvm)                 :: PWS                       !! average drought index on this PFT 
    !! during previous time steps.
    REAL(r_std), DIMENSION(npts, nvm)                 :: N_wood                    !! average biomass (gC m-2)
    !! of dead wood left on this PFT
    !! during previous time steps.
    REAL(r_std), DIMENSION(npts, nvmap)               :: i_hosts_dead_gp           !! Species group-level susceptibility 
    !! of bark beetle survival due to the 
    !! leftover of windthrow debris (unitless)
    REAL(r_std), DIMENSION(npts, nvmap)               :: i_hosts_defence_gp        !! Species group-level susceptibility 
    !! to bark beetle attacks due to weakening 
    !! from droughts (unitless)
    REAL(r_std), DIMENSION(npts, nvmap)               :: i_hosts_competition_gp    !! Species group-level susceptibility 
    !! to relative density of stem
    REAL(r_std), DIMENSION(npts, nvmap)               :: i_hosts_alive_gp          !! Species group-level susceptibility 
    REAL(r_std), DIMENSION(npts, nvmap)               :: i_beetles_massattack_gp

    !! of beetle survival to host availability
    !! due to the stand age (unitless)
    REAL(r_std), DIMENSION(npts, nvmap)               :: i_hosts_share_gp          !! Species group-level susceptibility 
    !! to bark beetle attacks due to
    !! the landscape heterogeneity (unitless)
    REAL(r_std), DIMENSION(npts, nvmap)               :: i_hosts_attractivity_gp   !! Species group-level susceptibility 
    !! to bark beetle attacks due to 
    !! the combined effect of windthrow, 
    !! stand age and


    REAL(r_std), DIMENSION(npts, nvmap)               :: i_beetles_survival_gp
    REAL(r_std), DIMENSION(npts, nvmap)               :: P_beetles_attacked_gp     !! risk index averaged for all
    !! age classes of the same PFT
    REAL(r_std), DIMENSION(npts, nvm)                 :: i_beetles_pressure_monitor!! variable used to monitor beetle pressure index (unitles)

    REAL(r_std), DIMENSION(npts, nvmap)               :: B_beetles_attacked_gp     !! Total biomass of a species group 
    !! that is infected by bark
    !! beetles (gC or gN m-2). 
    REAL(r_std), DIMENSION(npts, nvmap)               :: vegetmax_gp               !! share of all PFTs of the same species 
    !! (called species group) within
    REAL(r_std), DIMENSION(npts)                      :: vegetmax_deciduous        !! share of non-beetle trees ("smoke screen")
    !! within a pixel (0-1, unitless)
    REAL(r_std), DIMENSION(npts)                      :: vegetmax_beetle           !! share of bark beetle-host PFTs within
    !! a pixel (0-1, unitless)
    REAL(r_std), DIMENSION(npts, nvm)                 :: i_hosts_defence           !! Pft-level susceptibility
    !! due to weakening from droughts 
    !! (0-1, unitless)
    REAL(r_std), DIMENSION(npts, nvm)                 :: i_hosts_competition
    REAL(r_std), DIMENSION(npts, nvm)                 :: i_hosts_health            !! Pft-level susceptibility of tree 
    !! from bark beetle infestation 
    REAL(r_std), DIMENSION(npts, nvm)                 :: P_success                 !! probability that an infected tree die

    REAL(r_std), DIMENSION(npts, nvmap)               :: i_beetles_pressure
    REAL(r_std), DIMENSION(npts, nvm)                 :: B_beetles_kill            !! bark beetle damage in gC m-2 year-1
    REAL(r_std), DIMENSION(npts, nvmap)               :: i_beetles_generation_gp   !! Average bark beetle pressure during 
    !! the three recent years
    REAL(r_std), DIMENSION(npts, nvmap)               :: i_hosts_susceptibility_gp

    REAL(r_std)                                       :: qm_dia                    !! quadratic diameter
    REAL(r_std), DIMENSION(ncirc)                     :: qdia                      !! quadratic diameter per circumference class
    REAL(r_std), DIMENSION(npts, nvm)                 :: rdi                       !! relative density index
    REAL(r_std), DIMENSION(npts, nvmap)               :: rdi_gp                    !! rdi for all age-classe of a same PFT
    !! (0-1, unitless)
    REAL(r_std), DIMENSION(npts, nvmap)               :: B_total                   !! Actual biomass of the PFT (gC or gN m-2)
    REAL(r_std), DIMENSION(npts, nvmap)               :: B_wood                    !! Actual aboveground woody biomass of the PFT 
    !! (gC or gN m-2)
    REAL(r_std)                                       :: share                     !! Deciduous/coniferous share in a pixel
    REAL(r_std), DIMENSION(npts, nvm, legacy_years)   :: beetle_temp1              !! temporary variable
    REAL(r_std), DIMENSION(npts, nvm, beetle_legacy)  :: beetle_temp2
    LOGICAL, DIMENSION(npts, nvmap)                   :: beetle_gp                 !! Flag that indicates whether or not 
    !! bark beetle damage should
    !! be calculated for that species group. 
    REAL(r_std), DIMENSION(npts,nvm)                  :: wood_damage_all           !! Wood mass loss from wind throw  @tex $(gC m^{-2})$ @endtex


    !================================================================================================================================

    IF (firstcall_bark_beetle_damage) THEN
       !! Initialize local printlev
      ! Guillaume M. -- get_printlev calls getin_p, a COLLECTIVE bcast. Called here,
      ! inside the compute loop under a local firstcall, it hangs the run as soon as the
      ! ranks stop reaching it together. Use printlev, read once at startup by every rank.
       printlev_loc=printlev
       firstcall_bark_beetle_damage = .FALSE.
    END IF

    IF (printlev_loc .GE. 2) WRITE(numout, *) 'Entering stomate pest damage'
    ! 1.1 Initialize variables
    vegetmax_gp(:, :) = zero
    i_hosts_dead_gp(:, :) = zero
    i_hosts_defence_gp(:, :) = zero
    i_hosts_share_gp(:, :) = zero
    i_hosts_alive_gp(:, :) = zero
    i_hosts_competition_gp(:, :) = zero
    i_beetles_massattack_gp(:, :) = zero
    i_hosts_susceptibility_gp(:, :) = zero
    i_hosts_dead_monitor(:, :) = zero
    i_hosts_defence_monitor(:, :) = zero
    i_hosts_alive_monitor(:, :) = zero
    i_hosts_competition_monitor(:, :) = zero
    i_beetles_massattack_monitor(:, :) = zero
    i_beetles_generation_monitor(:, :) = zero
    i_beetles_survival_monitor(:, :) = zero
    i_beetles_activity_monitor(:, :) = zero
    i_hosts_attractivity_monitor(:, :) = zero
    i_beetles_pressure_monitor(:, :) = zero
    rdi(:, :) = zero
    rdi_gp(:, :) = zero
    i_hosts_competition(:, :) = zero
    B_total(:, :) = zero
    B_wood(:, :) = zero
    PWS(:, :) = zero
    N_wood(:, :) = zero
    B_beetles_kill(:, :) = zero
    P_success(:, :) = zero
    P_beetles_attacked_gp(:, :) = zero
    i_beetles_generation_gp(:, :) = zero
    nc(:, :, :) = zero
    vegetmax_deciduous = zero
    vegetmax_beetle = zero
    i_hosts_health = zero
    i_beetles_pressure = zero
    B_beetles_attacked_gp = zero
    wood_damage_all = zero
    beetle_damage_wood_vol = zero
    DR_beetles = zero

    ! Assume that we are not interested in calculating bark beetle
    ! damage for this species group.
    beetle_gp(:,:) = .FALSE.

    ! This is the main loop in wich the beetle damage will be calculated
    DO ipts = 1, npts

       ! 2. Preliminary PFT loop
       DO ivmap = 1, nvmap

          ! 2.1 Loop for wood and stand composition       
          ! Preliminary PFT loop is for calculating the variables shared within the  same
          ! species for all age classes. But it calculates twice; FIRST for wood  and the
          ! proportion of deciduous trees for ALL age classes; and SECOND for the rdi,
          ! biomass, and defence index for all ages EXCEPT the youngest age class.
          ! Basically, the module does not calculate bark beetle damage for trees with
          ! smaller than 7.5 centimeter of diameter and all trees in the youngest age class
          ! are smaller than that. But we think woody litter from the youngest stand can
          ! affect other ages by colonizing beetles and flying away, and deciduous trees
          ! from the youngest stand can affect by disturbing the bark beetle spreading.
          ! Therefore woody litter and veget_max of deciduous trees from the youngest class
          ! will be included in the first loop and PFT counter will be reset in the second
          ! loop.
          DO iagec = 1, nagec_pft(ivmap)
             ! Determine the PFT index for the age classes within the species
             ! group under consideration
             ivm = start_index(ivmap) + (iagec - 1)
             !-
             IF (.NOT. beetle_pft(ivm) .OR. (veget_max(ipts, ivm) .LE. min_stomate) ) THEN
                ! We don't want to calculate bark beetle infections for this PFT
                ! or there is no vegetation present in the PFT. nc is a counter.
                ! Because at least one PFT of the age class group has no beetles in it,
                ! the counter is decreased. It will be used later in other calculations.
                nc(ipts, ivmap, iagec) = 0
                ! Guillaume M. -- Non-host trees act as a smoke screen: their area is
                ! accumulated here, in the non-host branch.
                IF (is_tree(ivm)) THEN
                   vegetmax_deciduous(ipts) = vegetmax_deciduous(ipts) + veget_max(ipts, ivm)
                ENDIF
             ELSE
                nc(ipts, ivmap, iagec) = 1
                ! Calculate the total area of this species group across all age
                ! classes. Indicate that there bark beetles damage needs to be
                ! calculated in this species group.
                vegetmax_gp(ipts, ivmap) = vegetmax_gp(ipts, ivmap) + veget_max(ipts, ivm)
                beetle_gp(ipts, ivmap) = .TRUE.
                vegetmax_beetle(ipts) = vegetmax_beetle(ipts) + veget_max(ipts, ivm)

                ! Deciduous tree prevent picea trees to be infested by bark beetle
                ! by acting as "smoke screen". We include that effect in the code
                ! later through share susceptibility.
                ! Guillaume M. -- The smoke screen is now every non-host tree, not only
                ! the deciduous ones, so the area is accumulated in the branch above.
             ENDIF
          ENDDO ! iagec

          IF(printlev_loc >= 4)THEN
             WRITE(numout,*) 'ivmap, ipts', ivmap, ipts
             WRITE(numout,*) 'veget_max_gp for wood', vegetmax_gp(ipts, ivmap)
          END IF

          IF (beetle_gp(ipts, ivmap) .AND. vegetmax_gp(ipts, ivmap) .GT. min_stomate) THEN
             DO iagec = 1, nagec_pft(ivmap)
                IF (nc(ipts, ivmap, iagec) .EQ. 0) THEN
                   CYCLE
                ENDIF
                ivm = start_index(ivmap) + (iagec - 1)

                B_wood(ipts, ivmap) = B_wood(ipts, ivmap) + &
                     SUM((circ_class_biomass(ipts, ivm, :, isapabove,icarbon) + &
                     circ_class_biomass(ipts, ivm, :, iheartabove, icarbon)) * &
                     circ_class_n(ipts, ivm, :))
                ! Update the average biomass (gC or gN m-2)
                ! of dead wood left on this PFT during previous time steps.
                ! Only biomass added during the legacy period is accounted for.
                N_wood(ipts, ivmap) = N_wood(ipts,ivmap) + &
                     SUM(wood_leftover_legacy(ipts, ivm, :)) / legacy_years_wood

                IF(N_wood(ipts, ivmap) .LT. zero) THEN
                   ! needs a small error message
                   N_wood(ipts,ivmap) = zero
                ENDIF

             ENDDO
          ENDIF

          ! Guillaume M. -- The smoke-screen share is no longer computed here but in the
          ! main PFT loop below, where vegetmax_beetle (host area) is available.

          ! vegetmax_gp will be recalculated for other variables.
          ! In this case, the youngest age class will be excluded since anyways
          ! no damage is calculated with the youngest age class.
          vegetmax_gp(ipts,ivmap) = zero

          ! 2.2 Loop for rdi, biomass, and defence and reset nc
          DO iagec = 1, nagec_pft(ivmap)

             ! Determine the PFT index for the age classes within the species
             ! group under consideration
             ivm = start_index(ivmap) + (iagec - 1)

             IF(printlev_loc >= 4)THEN
                WRITE(numout,*) 'PFT, ivmap, ipts', ivm, ivmap, ipts
                WRITE(numout,*) 'beetle_pft', beetle_pft(ivm)
                WRITE(numout,*) 'beetle_group', beetle_gp(ipts, ivmap)
                WRITE(numout,*) 'veget_max', veget_max(ipts, ivm)
             END IF
             !-

             IF (.NOT. beetle_pft(ivm) .OR. (veget_max(ipts, ivm) .LE. min_stomate) &
                  .OR. (nagec_pft(ivmap) .GT. un .AND. ivm .EQ. start_index(ivmap)) &
                  .OR. (ALL(circ_class_n(ipts, ivm, :) .LE. min_stomate)) ) THEN
                ! We don't want to calculate bark beetle infections for this PFT
                ! or there is no vegetation present in the PFT. nc is a counter.
                ! Because at least one PFT of the age class group has no beetles
                ! in it,
                ! the counter is decreased. It will be used later in other
                ! calculations.
                nc(ipts, ivmap, iagec) = 0
             ELSE
                nc(ipts, ivmap, iagec) = 1
                !Relative density index (RDI) is needed here in order to
                !calculate susceptibility to RDI in the next circ_class loop 
                qm_dia = wood_to_qmdia( &
                     circ_class_biomass(ipts, ivm, :, :, icarbon), &
                     circ_class_n(ipts, ivm, :), ivm, pipe_tune2(ipts,ivm))

                ! Error checking
                IF (qm_dia == 0.0 ) THEN
                   WRITE(numout,*) 'ZERO VALUE detected before Nmax computation:', &
                        ' Dg = qm_dia =', qm_dia, &
                        ' biomass_temp = circ_class_biomass =', &
                        circ_class_biomass(ipts, ivm, :, :, icarbon), &
                        ' ind = circ_class_n =', circ_class_n(ipts, ivm, :), &
                        ' pft = ivm =', ivm, &
                        ' maxheight = pipe_tune2 =', pipe_tune2(ipts,ivm), &
                        ' ipts =', ipts
                   CALL flush(numout)
                   CALL ipslerr_p(3,'stomate_pest.f90', 'dia equals zero before calculating Nmax', &
                        'this should not happen', 'Check the code')
                END IF
                !-

                rdi(ipts, ivm) = (SUM(circ_class_n(ipts, ivm, :)) * m2_to_ha) / &
                     Nmax(qm_dia * m_to_cm, alpha_self_thinning(ipts,ivm), &
                     ivm, forest_managed(ipts, ivm))

                ! Calculate the total area of this species group across all age
                ! classes. Indicate that there bark beetles damage needs to be
                ! calculated in this species group.
                vegetmax_gp(ipts, ivmap) = vegetmax_gp(ipts, ivmap) + veget_max(ipts, ivm)
                beetle_gp(ipts, ivmap) = .TRUE.
                
             ENDIF

          ENDDO ! iagec

          IF (beetle_gp(ipts, ivmap) .AND. vegetmax_gp(ipts, ivmap) .GT. min_stomate) THEN
             DO iagec = 1, nagec_pft(ivmap)
                IF (nc(ipts, ivmap, iagec) .EQ. 0) THEN
                   CYCLE
                ENDIF
                ivm = start_index(ivmap) + (iagec - 1)
                ! All picea pft class are conrcerned by suceptibility to RDI
                rdi_gp(ipts, ivmap) = rdi_gp(ipts, ivmap) + rdi(ipts, ivm) * un/ SUM(nc(ipts, ivmap, :))

                ! Actual biomass of the PFT (gC or gN m-2)
                B_total(ipts, ivmap) = B_total(ipts, ivmap) + &
                     SUM(SUM(circ_class_biomass(ipts, ivm, :, :, icarbon), 2) * &
                     circ_class_n(ipts, ivm, :))
                PWS(ipts, ivm) = &
                     SUM(season_drought_legacy(ipts, ivm, :)) / legacy_years

             ENDDO
          ENDIF

       ENDDO ! DO ivmap

       ! 3. The main PFT loop 
       ! Bark beetle damage is calculated across all age classes of the
       ! same species group. Hence the loop over nvmap instead of nvm.
       DO ivmap = 1, nvmap

          ! The fist age class PFT use when PFT parameter is needed in
          ! equation outside iagec loop  
          xvm = start_index(ivmap)

          ! In the original model "share" represent the purity of the stand
          ! (monospecific vs mixed). This index reflects the fact that beetles are
          ! sensitive to deciduous (AGE_GROUP 6 and 8 for 15 PFT) species 
          ! which act as natural "smoke screen".

          !Debug
          IF(printlev_loc >= 4)THEN
             WRITE(numout,*) 'ivmap',ivmap
             WRITE(numout,*) 'vegetmax_group', vegetmax_gp(ipts, ivmap)
             WRITE(numout,*) 'rdi_group', rdi_gp(ipts, ivmap)
             WRITE(numout,*) 'vegetmax_deciduous', vegetmax_deciduous(ipts)
             WRITE(numout,*) 'beetle_group', beetle_gp(ipts, ivmap)
             WRITE(numout,*) 'nc', nc(ipts, ivmap, :)
          END IF
          !-

          IF (beetle_gp(ipts, ivmap) .AND. vegetmax_gp(ipts, ivmap) .GT. min_stomate) THEN
             ! Bark beetle damage needs to be calculated for this species group and
             ! there is vegetation present.
             DO iagec = 1, nagec_pft(ivmap)
                IF (nc(ipts, ivmap, iagec) .EQ. 0) THEN
                   CYCLE
                ENDIF
                ivm = start_index(ivmap) + (iagec - 1)

             ENDDO ! iagec

             ! Guillaume M. -- Smoke-screen share: non-host trees dilute host attractivity.
             ! The share is normalised by the host area (vegetmax_beetle), not by the
             ! species-group area.
             IF(vegetmax_deciduous(ipts) .GT. zero) THEN
                DO iagec = 1, nagec_pft(ivmap)
                   IF (nc(ipts, ivmap, iagec) .EQ. 0) THEN
                      CYCLE
                   ENDIF
                   ivm = start_index(ivmap) + (iagec - 1)
                   share = vegetmax_deciduous(ipts) / &
                      (vegetmax_beetle(ipts) + vegetmax_deciduous(ipts))
                   i_hosts_share_gp(ipts, ivmap) = &
                      un / (un + exp(S_share(ivm) * &
                      (share - SH_limit(ivm))))
                ENDDO
             ELSE
                i_hosts_share_gp(ipts, ivmap) = un
             END IF

             DO iagec = 1, nagec_pft(ivmap)
                IF (nc(ipts, ivmap, iagec) .EQ. 0) THEN
                   CYCLE
                ENDIF
                ivm = start_index(ivmap) + (iagec - 1)

                ! 3.2 Calculate susceptibility index
                i_beetles_generation_gp(ipts, ivmap) = i_beetles_generation_gp(ipts, ivmap) + &
                     i_beetles_generation(ipts, ivm, 1)/ &
                     (veget_max(ipts,ivm) / vegetmax_gp(ipts, ivmap))

                !+++JINA Do we still need this condition? nc should filter this?
                IF (.NOT. beetle_pft(ivm) .OR. veget_max(ipts, ivm) .LT. min_stomate) THEN
                   ! No reason to calculate bark beetle damage
                   CYCLE
                ELSE
                   IF(printlev_loc >= 4) THEN
                      WRITE(numout,*) 'ivmap, ivm',ivmap,ivm
                      WRITE(numout,*) 'wood_leftover', wood_leftover_legacy(ipts, ivm, :)
                      WRITE(numout,*) 'B_beeltes_kill_legacy', B_beetles_kill_legacy(ipts, ivm, :)
                      WRITE(numout,*) 'PWS', PWS(ipts, ivm)
                      WRITE(numout,*) 'B_total', B_total(ipts,ivmap)
                      WRITE(numout,*) 'B_wood', B_wood(ipts, ivmap)
                      WRITE(numout,*) 'RDI group', rdi_gp(ipts, ivmap)
                      WRITE(numout,*) 'RDI', rdi(ipts, ivm)
                   ENDIF

                   i_hosts_competition_gp(ipts, ivmap) = &
                        (un / (un + exp(S_competition(xvm) * &
                        (rdi_gp(ipts, ivmap) - RDi_limit(xvm)))))

                   i_hosts_competition(ipts, ivm) = &
                        (un / (un + exp(S_competition(xvm) * &
                        (rdi(ipts, ivm) - RDi_limit(xvm)))))

                   ! Long term drought susceptibility is an indirect measurement
                   ! tree defence to beetle colonisation. Multiple drought events
                   ! will decrease the amount of repelent produced by the trees.
                   i_hosts_defence(ipts,ivm) = un / (un + exp(S_defence(ivm) * &
                        ((un - season_drought_legacy(ipts, ivm, 1)) - PWS_limit(ivm))))
                   i_hosts_defence_gp(ipts, ivmap) = i_hosts_defence_gp(ipts, ivmap) + &
                        i_hosts_defence(ipts, ivm) * veget_max(ipts, ivm) / vegetmax_gp(ipts, ivmap)

                   i_hosts_health(ipts, ivm) =  (i_hosts_defence(ipts, ivm) + &
                        i_hosts_competition(ipts, ivm))/deux
                ENDIF
             ENDDO ! end loop nagec
             i_beetles_massattack_gp(ipts, ivmap) = un  / ( un + exp(S_massattack(xvm)* &
                  (P_beetles_attacked_legacy(ipts, xvm) - BP_limit(xvm))))

             i_hosts_susceptibility_gp(ipts, ivmap) = (un / (un + exp(S_susceptibility(xvm) * &
                  (rdi_gp(ipts, ivmap) - i_rd_susceptibility(xvm)))))     
         
             i_hosts_alive_gp(ipts, ivmap) = i_hosts_susceptibility_gp(ipts, ivmap) * &
                  i_beetles_massattack_gp(ipts, ivmap)

             ! The calculation of the windthrow suceptibility differs from the
             ! original paper of Temperli et al. 2013. First, the wood biomass
             ! lost by windthrow is replace by the wood leftover. This change
             ! implied that windthrow suceptibility is only representing beetle x
             ! windthrow interaction anymore. Of course, this index is still
             ! increasing drasticaly after a windthrow event. Second, the
             ! equation is using the actual woody biomass instead of the maximal
             ! potential woody biomass, a fix value of 300 t/ha in the paper, in
             ! ORCHIDEE we never reach this value and it will
             ! underestimate the contribution of windthrows.
             ! The amount of dead wood left on-site by previous wind storms is
             ! smaller than the living aboveground biomass. This will be accounted
             ! for to decrease the pft-level susceptibility to bark beetle attacks
             ! due to the left-over of windthrow debris. When 30% of
             ! the current wooody biomass is on the floor
             ! susceptibility reach its maximum(1)(unitless)
             i_hosts_dead_gp(ipts, ivmap) = MIN(un, (N_wood(ipts, ivmap) / &
                  B_wood(ipts, ivmap)) / max_Nwood(xvm))

             ! Calculate pixel attractivity at the level of the whole species group. It is
             ! assumed that half of the vulnerability comes from the relative density index. 
             ! But here RDI is not used as a proxy for competition but rather
             ! a proxy for bettle spead rate wich is stimulatated when tree
             ! canopies are close to each other. Long term drought (proxy of plant health) 
             ! and pixel deciduous share also account in the calculation. 
             ! wood leftover is also used to trigger the epidemic phase by increasing
             ! the risk index right after a storm but don't have any impact
             ! during the epidemic phase.
             ! Guillaume M. -- Third survival channel: drought-weakened host defence scaled
             ! by voltinism (i_beetles_generation). Without it survival only follows the
             ! alive/dead wood channels, so acute drought alone cannot start an epidemic.
             i_beetles_survival_gp(ipts, ivmap) = &
                  MAX(i_hosts_alive_gp(ipts, ivmap), &
                  i_hosts_dead_gp(ipts, ivmap), &
                  i_hosts_defence_gp(ipts,ivmap) * &
                        (1 + 0.4 * (i_beetles_generation_gp(ipts,ivmap) - 1)) )

             i_hosts_attractivity_gp(ipts, ivmap) = &
                  MAX(i_hosts_defence_gp(ipts, ivmap) , &
                  i_hosts_competition_gp(ipts, ivmap)) * &
                  i_hosts_share_gp(ipts,ivmap) 

             ! Assign the values calculated at the PFT-level to each individual age class
             ! in ORCHIDEE.
             DO iagec = 1, nagec_pft(ivmap)
                IF (nc(ipts, ivmap, iagec) .EQ. 0) THEN
                   CYCLE
                ENDIF
                ivm = start_index(ivmap) + (iagec - 1)
                i_hosts_attractivity_monitor(ipts, ivm) = i_hosts_attractivity_gp(ipts, ivmap)
                i_beetles_survival_monitor(ipts, ivm) = i_beetles_survival_gp(ipts,ivmap)
                i_hosts_dead_monitor(ipts, ivm) = i_hosts_dead_gp(ipts,ivmap)
                i_hosts_alive_monitor(ipts, ivm) = i_hosts_alive_gp(ipts, ivmap)  
                i_hosts_competition_monitor(ipts, ivm) = i_hosts_competition_gp(ipts, ivmap)
                i_beetles_activity_monitor(ipts, ivm) = i_beetles_activity_legacy(ipts, xvm, 1)
                i_beetles_massattack_monitor(ipts, ivm) = i_beetles_massattack_gp(ipts, ivmap)
                i_hosts_defence_monitor(ipts, ivm) = i_hosts_defence_gp(ipts, ivmap)
                i_beetles_generation_monitor(ipts, ivm) = i_beetles_generation_gp(ipts, ivmap) 
             END DO


             !Debug
             IF(printlev_loc >= 4) THEN
                WRITE(numout,*) 'i_hosts_defence', i_hosts_defence_gp(ipts, ivmap)
                WRITE(numout,*) 'i_beetles_massattack', i_beetles_massattack_gp(ipts,ivmap)
                WRITE(numout,*) 'i_hosts_alive',i_hosts_alive_gp(ipts, ivmap)
                WRITE(numout,*) 'i_hosts_competition', i_hosts_competition_gp(ipts, ivmap)
                WRITE(numout,*) 'i_hosts_share', i_hosts_share_gp(ipts, ivmap)
                WRITE(numout,*) 'i_hosts_dead', i_hosts_dead_gp(ipts, ivmap)
             END IF
             !-

             ! Calculate as the real pressure of the beetles. Avegare of the average bark
             ! beetle of the past three years and the beetle legacy (which is a ratio between
             ! affected dead wood (breeding grounds) and fresh biomass (food)). This average
             ! is an index for the available beetles when multiplied with the susceptibility
             ! of the species group we get the beetle pressure (which could be seen as a proxy
             ! for beetle biomass). 
             i_beetles_pressure(ipts, ivmap) = i_beetles_survival_gp(ipts, ivmap) * &
                  (i_beetles_generation_gp(ipts, ivmap) + &
                  i_beetles_activity_legacy(ipts, xvm, 1)) / deux

             ! Estimate the risk index at the pixel level.
             P_beetles_attacked_gp(ipts, ivmap) =  i_hosts_attractivity_gp(ipts, ivmap) * &
                  i_beetles_pressure(ipts, ivmap)

             ! This loop is use for monitoring intermediate result.
             ! ivmap ---> ivm, because there is no ivmap dimension in the
             ! output files 
             DO iagec = 1, nagec_pft(ivmap)
                IF (nc(ipts, ivmap, iagec) .EQ. 0) THEN
                   CYCLE
                ENDIF
                ivm = start_index(ivmap) + (iagec - 1)
                i_beetles_pressure_monitor(ipts, ivm) = i_beetles_pressure(ipts, ivmap)
                P_beetles_attacked_legacy(ipts, ivm) = i_beetles_pressure(ipts, ivmap)
             END DO

             !Debug
             IF(printlev_loc >= 4) THEN
                WRITE(numout,*) 'i_beetles_generation_gp', i_beetles_generation_gp(ipts, ivmap)
                WRITE(numout,*) "i_beetles_survival", i_beetles_survival_gp(ipts, ivmap)
                WRITE(numout,*) 'i_beetles_pressure', i_beetles_pressure(ipts, ivmap)
                WRITE(numout,*) 'i_beetles_activity_legacy', i_beetles_activity_legacy(ipts, xvm, 1)
                WRITE(numout,*) 'P_beetles_attacked_gp', P_beetles_attacked_gp(ipts, ivmap)
                WRITE(numout,*) 'i_hosts_attractivity_gp', i_hosts_attractivity_gp(ipts, ivmap)
             END IF
             !-

             ! 3.4 Calculate Bark beetle damage                
             ! Estimate the total biomass of spruce attacked by bark beetle. The
             ! calculation diverge from the original equation which used a fix
             ! value corresponding the half of the maximun potential biomass (also fix).
             ! In our case maximun potential biomass is not a fix value and depend
             ! climate forcing, soil properties, interspecies water competition. Seems
             ! easier to use the actual biomass and recalibrate suceptibility parameters
             ! in order to decrease the overall risk index values.
             B_beetles_attacked_gp(ipts, ivmap) = B_total(ipts, ivmap) * P_beetles_attacked_gp(ipts, ivmap)

             ! Distribute the total infected biomass over the different age classes
             ! of the same species
             DO iagec = 1, nagec_pft(ivmap)
                ! Determine the PFT index for the age classes within the species
                ! group under consideration
                ivm = start_index(ivmap) + (iagec - 1)

                ! Estimate mortality rate at the PFT level (unitless).
                P_success(ipts, ivm) = (i_hosts_health(ipts, ivm) + &
                     i_beetles_pressure(ipts, ivmap)) / deux

                ! Estimate living biomass killed by bark beetles (gC or gN m-2)
                B_beetles_kill(ipts, ivm) = P_success(ipts, ivm) * &
                     B_beetles_attacked_gp(ipts, ivmap) * &
                     veget_max(ipts, ivm) / vegetmax_gp(ipts, ivmap)

                beetle_temp2(ipts, ivm, :) = zero

                DO iyear = 1, beetle_legacy - 1
                   beetle_temp2(ipts, ivm, iyear + 1) = &
                        B_beetles_kill_legacy(ipts, ivm, iyear)
                ENDDO
                beetle_temp2(ipts, ivm, 1) = B_beetles_kill(ipts, ivm)
                B_beetles_kill_legacy(ipts, ivm, :) = beetle_temp2(ipts, ivm, :)

                IF (nc(ipts, ivmap, iagec) .EQ. 0) THEN
                   CYCLE
                ENDIF

                !Debug
                IF(printlev_loc >= 4) THEN
                   WRITE(numout,*) 'B_total', B_total(ipts,ivmap)
                   WRITE(numout,*) 'B_beetles_attacked_gp', B_beetles_attacked_gp(ipts, ivmap)
                   WRITE(numout,*) 'B_beetles_kill', B_beetles_kill(ipts, ivm)
                   WRITE(numout,*) 'DR_beetles', B_beetles_kill(ipts, ivm) / &
                        B_total(ipts,ivmap) * 100
                   WRITE(numout,*) 'P_success', P_success(ipts, ivm)
                   WRITE(numout,*) 'vegetmax_group', vegetmax_gp(ipts, ivmap)
                   WRITE(numout,*) 'veget_max', veget_max(ipts, ivm)
                END IF
                !-
                ! Select which circ_class would be damaged based on the thier
                ! diameters
                ncirc_loc = 0
                qdia(:) = wood_to_dia(circ_class_biomass(ipts, ivm, :, :, icarbon), &
                     ivm, pipe_tune2(ipts,ivm))
                DO icir = 1, ncirc
                   IF (qdia(icir) > 0.075) THEN 
                      !filter based on D (rather than age in Temperli)
                      ncirc_loc = ncirc_loc + 1
                   ENDIF
                ENDDO

                ! By using an decreasing index, I made the assumption that the
                ! trees which die first from bark beetle infection are the bigest one.
                DO icir = ncirc_loc, 1, -1
                   IF(SUM(circ_class_biomass(ipts, ivm, icir, :, icarbon)) * &
                        circ_class_n(ipts, ivm, icir) .LT. B_beetles_kill(ipts, ivm)) THEN

                      ! All individuals in the circumference class need to be killed to satisfy
                      ! the total amount of biomass that needs to be killed. Note that all of 
                      ! the above code is calculated all the time (thus irrespective of whether
                      ! pest are accounted for or not). The damage in the code below will only
                      ! be accounted for when ok_pest = .TRUE.
                      IF(ok_pest) THEN 
                         circ_class_kill(ipts, ivm, icir, forest_managed(ipts, ivm), icut_beetle) = &
                              circ_class_n(ipts, ivm, icir)
                      ENDIF
                      B_beetles_kill(ipts, ivm) = B_beetles_kill(ipts, ivm) - &
                           SUM(circ_class_biomass(ipts,ivm, icir, :, icarbon)) * circ_class_n(ipts, ivm, icir)
                   ELSE
                      IF(B_beetles_kill(ipts, ivm) .GT. zero) THEN
                         ! Only part of the trees in this circumference class should be killed
                         IF(ok_pest) THEN 
                            circ_class_kill(ipts, ivm, icir, forest_managed(ipts, ivm), icut_beetle) = &
                                 circ_class_n(ipts, ivm, icir) * B_beetles_kill(ipts, ivm) / &
                                 (SUM(circ_class_biomass(ipts, ivm, icir, :, icarbon)) * circ_class_n(ipts, ivm, icir))
                         ENDIF
                         B_beetles_kill(ipts, ivm) = zero
                      ELSE
                         ! The total demand of mortality has been satisfied. No more killing is
                         ! required.
                         circ_class_kill(ipts, ivm, icir, forest_managed(ipts, ivm), icut_beetle) = zero
                      ENDIF
                   ENDIF

                   wood_damage_all(ipts,ivm) = wood_damage_all(ipts,ivm) + &
                        (circ_class_biomass(ipts, ivm, icir, isapabove, icarbon) + &
                        circ_class_biomass(ipts, ivm, icir, iheartabove, icarbon)) * (un - branch_ratio(ivm)) *  &
                        circ_class_kill(ipts, ivm, icir, forest_managed(ipts, ivm), icut_beetle)

                ENDDO ! End loop ncirc
                beetle_damage_wood_vol(ipts, ivm) = wood_damage_all(ipts, ivm) / (pipe_density(ivm))
                DR_beetles(ipts, ivm) = B_beetles_kill_legacy(ipts, ivm, 1) / &
                     SUM(SUM(circ_class_biomass(ipts, ivm, :, :,icarbon), 2) * circ_class_n(ipts, ivm, :))

                ! Calculate Beetle pressure index. The actual biomass(biomass_current) is
                ! used instead of the maximal potential biomass, i.e., 300 t/ha, as
                ! was done in temperli et al. 2013. The
                ! ratio between dead wood from previous beetle outbreaks (breeding
                ! grounds) and uninfected biomass (fresh food) is seen as a driver of new
                ! infections. Make sure the legacy stays inbetween zero and one.
                beetle_temp1(ipts, ivm, :) = zero

                DO iyear = 1, legacy_years - 1
                   beetle_temp1(ipts, ivm, iyear + 1) = i_beetles_activity_legacy(ipts, ivm, iyear)
                ENDDO

                beetle_temp1(ipts, ivm, 1) = un / (un + exp(S_activity(ivm)* &
                     (DR_beetles(ipts, ivm) - act_limit(ivm))))

                i_beetles_activity_legacy(ipts, ivm, :) = beetle_temp1(ipts, ivm, :)

             ENDDO   ! End loop nagec

          ELSE
             ! No mortality from bark beetles. Set the mortality to zero to avoid problems
             ! with uninitialized values
             DO iagec = 1, nagec_pft(ivmap) 
                IF (nc(ipts, ivmap, iagec) .EQ. 0) THEN !++JINA Do we need this line?
                   CYCLE
                ENDIF
                ! Determine the PFT index for the age classes within the species
                ! group under consideration
                ivm = start_index(ivmap) + (iagec - 1)
                ! No trees are to be killed
                circ_class_kill(ipts, ivm, :, forest_managed(ipts, ivm), icut_beetle) = zero
             ENDDO ! nagec
          ENDIF ! beetle_group .EQV. .TRUE. .AND. vegetmax_group(ipts,ivmap) .GT. min_stomate

       ENDDO ! End loop nvmap     
    ENDDO ! End loop npts

    CALL xios_orchidee_send_field("i_hosts_dead",i_hosts_dead_monitor)
    CALL xios_orchidee_send_field("i_hosts_competition",i_hosts_competition_monitor)
    CALL xios_orchidee_send_field("i_hosts_alive",i_hosts_alive_monitor)
    CALL xios_orchidee_send_field("i_beetles_massattack",i_beetles_massattack_monitor)
    CALL xios_orchidee_send_field("i_hosts_defence",i_hosts_defence_monitor)
    CALL xios_orchidee_send_field("i_beetles_generation",i_beetles_generation_monitor)
    CALL xios_orchidee_send_field("i_beetles_pressure",i_beetles_pressure_monitor)
    CALL xios_orchidee_send_field("i_beetles_survival",i_beetles_survival_monitor)
    CALL xios_orchidee_send_field("i_hosts_attractivity",i_hosts_attractivity_monitor)
    CALL xios_orchidee_send_field("i_beetles_activity",i_beetles_activity_monitor)
    CALL xios_orchidee_send_field("B_beetles_kill",B_beetles_kill_legacy(:,:,1))
    CALL xios_orchidee_send_field("Beetle_damage_woodvol",beetle_damage_wood_vol)
    CALL xios_orchidee_send_field("DR_beetles", DR_beetles)

    ! Guillaume M. -- Beetle-affected area fed to the edge_length feedback: killed
    ! biomass (gC m-2 yr-1) is converted through pest_biomass_ref, as
    ! B_beetles_kill / pest_biomass_ref * veget_max * cell area, summed over tree PFTs.
    ! Cell area is resolution(:,1)*resolution(:,2) (m2).
    IF (ok_aed_feedback .AND. ALLOCATED(dA_pest)) THEN
       IF (pest_biomass_ref > min_stomate) THEN
          DO ivm = 1, nvm
             IF (is_tree(ivm)) THEN
                ! Guillaume M. -- Gate on the epidemic phase: only above
                ! aed_beetle_epidemic_threshold do beetles overcome host defences and kill
                ! at stand level, opening a gap. Below it the kills are endemic background
                ! mortality of isolated trees, which opens no forest/non-forest edge and
                ! must stay out of the edge budget.
                WHERE (i_beetles_massattack_monitor(:, ivm) > aed_beetle_epidemic_threshold)
                   dA_pest(:) = dA_pest(:) +                                    &
                        B_beetles_kill(:, ivm) / pest_biomass_ref *            &
                        veget_max(:, ivm) * resolution(:, 1) * resolution(:, 2)
                END WHERE
             ENDIF
          ENDDO
       ENDIF
    ENDIF

    ! Guillaume M. -- Perturbed forest area from beetle kills, through the shared
    ! cc_kill_to_area basal-area route: same metric and unit (m2) as
    ! AREA_DAMAGED_{FIRE,STORM,HARVEST}. Built from circ_class_kill(icut_beetle), so it
    ! follows the mortality actually applied to the stand. Annual field, emitted once a
    ! year here and zero-filled by pest_write on the other time steps.
    area_damaged_pest(:,:) = zero
    DO ivm = 1, nvm
       IF (is_tree(ivm)) THEN
          CALL cc_kill_to_area(npts, ivm, SUM(circ_class_kill(:,ivm,:,:,icut_beetle), DIM=3), &
               circ_class_biomass(:,ivm,:,:,:), circ_class_n(:,ivm,:), veget_max(:,ivm),      &
               resolution(:,1)*resolution(:,2), pipe_tune2(:,ivm), area_damaged_pest(:,ivm))
       ENDIF
    ENDDO
    CALL xios_orchidee_send_field("AREA_DAMAGED_PEST", area_damaged_pest)

  END SUBROUTINE bark_beetle_damage

!!
!==================================================================================================================================
!!  SUBROUTINE   : pest_write
!!
!>\BRIEF        Wirte xios outputs when .NOT. ts_annual_proc to avoid that
!values are written only once.
!! 
!_
!================================================================================================================================

  SUBROUTINE pest_write (npts, B_beetles_kill_legacy)

    !! 0. Variable declarations

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                  :: npts                       !! number of pixels

    REAL(r_std), DIMENSION(:, :, :), INTENT(in) :: B_beetles_kill_legacy      !! temporary variable(gC.m-2)
    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables
    REAL(r_std), DIMENSION(npts, nvm)          :: beetle_damage_wood_vol      !! Wood volume loss from bark beetle @tex $(m^3 m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts, nvm)          :: DR_beetles                  !! rate of biomass loss from bark beetle per pft
    REAL(r_std), DIMENSION(npts, nvm)          :: area_damaged_pest           !! Forest area killed by bark beetle, cc-kill basal-area route (m2)

    !! infected during the previous timestep  
    REAL(r_std), DIMENSION(npts, nvm)          :: i_hosts_dead_monitor        !! variable used to monitor windthrow susceptibility (unitless)
    REAL(r_std), DIMENSION(npts, nvm)          :: i_hosts_attractivity_monitor!! variable used to monitor overall susceptibility index (unitless)
    REAL(r_std), DIMENSION(npts, nvm)          :: i_hosts_competition_monitor !! variable used to monitor rdi susceptibility (unitless)
    REAL(r_std), DIMENSION(npts, nvm)          :: i_hosts_alive_monitor
    REAL(r_std), DIMENSION(npts, nvm)          :: i_beetles_pressure_monitor  !! variable used to monitor beetle pressure index (unitles)
    !! variable used to monitor rdi susceptibility (unitless)
    REAL(r_std), DIMENSION(npts, nvm)          :: i_beetles_massattack_monitor
    REAL(r_std), DIMENSION(npts, nvm)          :: i_beetles_survival_monitor
    !! variable used to monitor rdi susceptibility (unitless)
    REAL(r_std), DIMENSION(npts, nvm)          :: i_beetles_activity_monitor
    !! variable used to monitor rdi susceptibility (unitless)
    REAL(r_std), DIMENSION(npts, nvm)          :: i_hosts_defence_monitor     !! variable used to monitor ability of the trees to defend themselves against
                                                                              !! bark beetle attacks (unitless)
    REAL(r_std), DIMENSION(npts, nvm)          :: i_beetles_generation_monitor!! variable used to monitor generation index susceptibility
    
!================================================================================================================================

    IF (printlev.GE.2) WRITE(numout,*) 'Entering pest_write'

    i_hosts_dead_monitor(:,:) = zero
    i_hosts_competition_monitor(:,:) = zero
    i_hosts_alive_monitor(:,:) = zero
    i_beetles_massattack_monitor(:,:) = zero
    i_hosts_defence_monitor(:,:) = zero
    i_beetles_generation_monitor(:,:) = zero
    i_beetles_pressure_monitor(:,:) = zero
    i_beetles_survival_monitor(:,:) = zero
    i_hosts_attractivity_monitor(:,:) = zero
    i_beetles_activity_monitor(:,:) = zero
    beetle_damage_wood_vol(:,:) = zero
    DR_beetles(:,:) = zero
    area_damaged_pest(:,:) = zero

    CALL xios_orchidee_send_field("i_hosts_dead",i_hosts_dead_monitor)
    CALL xios_orchidee_send_field("i_hosts_competition",i_hosts_competition_monitor)
    CALL xios_orchidee_send_field("i_hosts_alive",i_hosts_alive_monitor)
    CALL xios_orchidee_send_field("i_beetles_massattack",i_beetles_massattack_monitor)
    CALL xios_orchidee_send_field("i_hosts_defence",i_hosts_defence_monitor)
    CALL xios_orchidee_send_field("i_beetles_generation",i_beetles_generation_monitor)
    CALL xios_orchidee_send_field("i_beetles_pressure",i_beetles_pressure_monitor)
    CALL xios_orchidee_send_field("i_beetles_survival",i_beetles_survival_monitor)
    CALL xios_orchidee_send_field("i_hosts_attractivity",i_hosts_attractivity_monitor)
    CALL xios_orchidee_send_field("i_beetles_activity",i_beetles_activity_monitor)
    CALL xios_orchidee_send_field("B_beetles_kill",B_beetles_kill_legacy(:,:,1))
    CALL xios_orchidee_send_field("Beetle_damage_woodvol",beetle_damage_wood_vol)
    CALL xios_orchidee_send_field("DR_beetles", DR_beetles)
    CALL xios_orchidee_send_field("AREA_DAMAGED_PEST", area_damaged_pest)

  END SUBROUTINE pest_write

  
END MODULE stomate_pest

