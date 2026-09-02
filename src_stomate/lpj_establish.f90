! =================================================================================================================================
! MODULE       : lpj_establish
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        Establish pft's
!!
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S)	: 
!! - Sitch, S., B. Smith, et al. (2003), Evaluation of ecosystem dynamics,
!!        plant geography and terrestrial carbon cycling in the LPJ dynamic 
!!        global vegetation model, Global Change Biology, 9, 161-185.\n
!! - Haxeltine, A. and I. C. Prentice (1996), BIOME3: An equilibrium
!!        terrestrial biosphere model based on ecophysiological constraints, 
!!        resource availability, and competition among plant functional types,
!!        Global Biogeochemical Cycles, 10(4), 693-709.\n
!! - Smith, B., I. C. Prentice, et al. (2001), Representation of vegetation
!!        dynamics in the modelling of terrestrial ecosystems: comparing two
!!        contrasting approaches within European climate space,
!!        Global Ecology and Biogeography, 10, 621-637.\n
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_stomate/lpj_establish.f90 $
!! $Date: 2023-11-23 17:13:54 +0100 (jeu. 23 nov. 2023) $
!! $Revision: 8318 $
!! \n
!_ ================================================================================================================================

MODULE lpj_establish

  ! modules used:
  USE xios_orchidee
  USE ioipsl_para
  USE stomate_data
  USE constantes
  USE dynamic_parameters
  USE grid
  USE function_library,    ONLY: biomass_to_lai

  IMPLICIT NONE

  ! private & public routines
  PRIVATE
  PUBLIC establish,establish_clear

  LOGICAL, SAVE                          :: firstcall_establish = .TRUE.           !! first call
!$OMP THREADPRIVATE(firstcall_establish)
CONTAINS


!! ================================================================================================================================
!! SUBROUTINE   : fire_clear 
!!
!>\BRIEF       Set the firstcall_establish flag to .TRUE. and activate initialization
!_ ================================================================================================================================

  SUBROUTINE establish_clear
    firstcall_establish = .TRUE.
  END SUBROUTINE establish_clear


! =================================================================================================================================
! SUBROUTINE   : establish 
!
!>\BRIEF       Calculate sstablishment of new woody PFT and herbaceous PFTs
!!
!! DESCRIPTION : Establishments of new woody and herbaceous PFT are simulated. 
!! Maximum establishment rate (0.12) declines due to competition for light (space).
!! There are two establishment estimates: one for the for DGVM and one for the 
!! static cases.\n
!! In the case of DGVM, competitive process of establishment for the area of 
!! available space is represented using more detailed description compared with static 
!! one. Biomass and distribution of plant age are updated on the basis of changes 
!! in number of individuals. Finally excess sapwood of is converted to heartwood.
!!
!! \latexonly
!! \input{equation_lpj_establish.tex}
!! \endlatexonly
!! \n
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S)    :
!! Smith, B., I. C. Prentice, et al. (2001), Representation of vegetation
!!    dynamics in the modelling of terrestrial ecosystems: comparing two
!!    contrasting approaches within European climate space,
!!    Global Ecology and Biogeography, 10, 621-637.
!!
!! FLOWCHART       : 
!! \latexonly
!! \includegraphics[scale = 0.7]{establish.png}
!! \endlatexonly
!! \n
!_ ================================================================================================================================
 
  SUBROUTINE establish (npts, dt, PFTpresent, regenerate, &
       neighbours, resolution, need_adjacent, herbivores, &
       precip_annual, gdd0, lm_lastyearmax, &
       cn_ind, lai, avail_tree, avail_grass,  npp_longterm, &
       leaf_age, leaf_frac, &
       ind, biomass, age, everywhere, co2_to_bm, &
       soil_n_min, n_uptake, nstress_season, veget_cov_max, woodmass_ind, &
       mortality, bm_sapl_2D, bm_to_litter)

    !! 0. Variable and parameter declaration
    
    !! 0.1 Input variables
    
    INTEGER(i_std), INTENT(in)                                :: npts            !! Domain size - number of pixels (dimensionless)    
    REAL(r_std), INTENT(in)                                   :: dt              !! Time step of vegetation dynamics for stomate 
                                                                                 !! (days)            
    LOGICAL, DIMENSION(npts,nvm), INTENT(in)                  :: PFTpresent      !! Is pft there (unitless)    
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: regenerate      !! Winter sufficiently cold (unitless)   
    INTEGER(i_std), DIMENSION(npts,NbNeighb), INTENT(in)      :: neighbours      !! indices of the neighbours of each grid point 
                                                                                 !! (unitless);  
                                                                                 !! [1=North and then clockwise]  
    REAL(r_std), DIMENSION(npts,2), INTENT(in)                :: resolution      !! resolution at each grid point (m); 1=E-W, 2=N-S     
    LOGICAL, DIMENSION(npts,nvm), INTENT(in)                  :: need_adjacent   !! in order for this PFT to be introduced, does it
                                                                                 !! have to be present in an adjacent grid box?  
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: herbivores      !! time constant of probability of a leaf to 
                                                                                 !! be eaten by a herbivore (days)     
    REAL(r_std), DIMENSION(npts), INTENT(in)                  :: precip_annual   !! annual precipitation (mm year^{-1}) 
    REAL(r_std), DIMENSION(npts), INTENT(in)                  :: gdd0            !! growing degree days (degree C)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: lm_lastyearmax  !! last year's maximum leaf mass for each PFT 
                                                                                 !! (gC m^{-2 })
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: cn_ind          !! crown area of individuals (m^2)        
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: lai             !! leaf area index OF an individual plant 
                                                                                 !! (m^2 m^{-2})           
    REAL(r_std), DIMENSION(npts), INTENT(in)                  :: avail_tree      !! space availability for trees (unitless)     
    REAL(r_std), DIMENSION(npts), INTENT(in)                  :: avail_grass     !! space availability for grasses (unitless)  
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: npp_longterm    !! longterm NPP, for each PFT (gC m^{-2})   
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: veget_cov_max   !! "maximal" coverage fraction of a PFT 
                                                                                 !! (LAI -> infinity) on ground (unitless)
    REAL(r_std), DIMENSION(npts,nvm),INTENT(in)               :: mortality       !! Fraction of individual dying this time 
                                                                                 !! step (0 to 1, unitless) 
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements),INTENT(in)  :: bm_sapl_2D

    !! 0.2 Output variables
    
    !! 0.3 Modified variables

    REAL(r_std), DIMENSION(npts,nvm,nleafages), INTENT(inout) :: leaf_age        !! leaf age (days)     
    REAL(r_std), DIMENSION(npts,nvm,nleafages), INTENT(inout) :: leaf_frac       !! fraction of leaves in leaf age class (unitless)    
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: ind             !! Number of individuals (individuals m^{-2})          
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout):: biomass   !! biomass (gC m^{-2 })     
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: age             !! mean age (years)       
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: everywhere      !! is the PFT everywhere in the grid box or very 
                                                                                 !! localized (unitless)    
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: co2_to_bm       !! biomass up take for establishment i.e. 
                                                                                 !! pseudo-photosynthesis (gC m^{-2} day^{-1}) 
    REAL(r_std), DIMENSION(npts,nvm,nnspec), INTENT(inout)    :: soil_n_min      !! Nitrogen in mineral soil,    
 !! (gN/(m**2 of PFT)) 
 !! to correct for nitrogen creation associated with establishmnet   
    REAL(r_std), DIMENSION(npts,nvm,nionspec), INTENT(inout)  :: n_uptake        !! Plant uptake of ammonium and nitrate  
 !! (gN/m**2/day) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: nstress_season  !! seasonal nitrogen stress  
 !! (unitless) 
 !! (get's initialised at first establishment) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: woodmass_ind    !! woodmass of the individual, needed to calculate
                                                                                 !! crownarea in lpj_crownarea (gC m^{-2 })  
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout) :: bm_to_litter !!Biomass transfer to litter

    !! 0.4 Local variables

    REAL(r_std)                                               :: tau_eatup       !! time during which a sapling can be entirely 
                                                                                 !! eaten by herbivores (days) 
    REAL(r_std), DIMENSION(npts,nvm)                          :: fpc_nat         !! new fpc, foliage projective cover: fractional
                                                                                 !! coverage (unitless)       
    REAL(r_std), DIMENSION(npts)                              :: estab_rate_max_climate_tree  !! maximum tree establishment rate, 
                                                                                              !! based on climate only (unitless)  
    REAL(r_std), DIMENSION(npts)                              :: estab_rate_max_climate_grass !! maximum grass establishment rate,
                                                                                              !! based on climate only (unitless)
    REAL(r_std), DIMENSION(npts)                              :: estab_rate_max_tree          !! maximum tree establishment rate, 
                                                                                              !! based on climate and fpc 
                                                                                              !! (unitless) 
    REAL(r_std), DIMENSION(npts)                              :: estab_rate_max_grass         !! maximum grass establishment rate,
                                                                                              !! based on climate and fpc 
                                                                                              !! (unitless) 
    REAL(r_std), DIMENSION(npts)                              :: sumfpc          !! total natural fpc (unitless)
    REAL(r_std), DIMENSION(npts)                              :: fracnat         !! total fraction occupied by natural 
                                                                                 !! vegetation (unitless)  
    REAL(r_std), DIMENSION(npts)                              :: sumfpc_wood     !! total woody fpc (unitless)    
    REAL(r_std), DIMENSION(npts)                              :: spacefight_grass!! for grasses, measures the total concurrence 
                                                                                 !! for available space (unitless)
    REAL(r_std), DIMENSION(npts,nvm)                          :: d_ind           !! change in number of individuals per time step 
                                                                                 !! (individuals m^{-2} day{-1})         
    REAL(r_std), DIMENSION(npts,nelements)                    :: bm_new          !! biomass increase (gC m^{-2 })        
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements)         :: biomass_old     !! Save the original biomass passed into the subroutine 
    REAL(r_std), DIMENSION(npts,nelements)                    :: bm_non          !! Non-effective establishment: the "virtual" saplings that die instantly
    REAL(r_std), DIMENSION(npts,nelements)                    :: bm_eff          !! Effective (or real) establishment  
    REAL(r_std), DIMENSION(npts)                              :: dia             !! stem diameter (m)    
    REAL(r_std), DIMENSION(npts)                              :: b1              !! temporary variable           
    REAL(r_std), DIMENSION(npts)                              :: woodmass        !! woodmass of an individual (gC m^{-2})  
    REAL(r_std), DIMENSION(npts)                              :: leaf_mass_young !! carbon mass in youngest leaf age class 
                                                                                 !! (gC m^{-2})
    REAL(r_std), DIMENSION(npts)                              :: factor          !! reduction factor for establishment if many 
                                                                                 !! trees or grasses are present (unitless)  
    REAL(r_std), DIMENSION(npts,nelements)                    :: total_bm        !! Total carbon (or nitrogen) mass for all pools (gC m^{-2})    
    REAL(r_std), DIMENSION(npts,nelements)                    :: total_bm_sapl   !! Total sappling biomass for all pools 
                                                                                 !! (gC(orN) m^{-2})  
    REAL(r_std), DIMENSION(npts,nelements)                    :: total_bm_sapl_non !! total non-effective sapling biomass 

    INTEGER(i_std)                                            :: nfrontx         !! from how many sides is the grid box invaded
                                                                                 !! (unitless?)   
    INTEGER(i_std)                                            :: nfronty         !! from how many sides is the grid box invaded
                                                                                 !! (unitless?)   
    REAL(r_std), DIMENSION(npts)                              :: vn              !! flow due to new individuals veget_cov_max after 
                                                                                 !! establishment, to get a proper estimate of 
                                                                                 !! carbon and nitrogen 
    REAL(r_std), DIMENSION(npts)                              :: lai_ind         !! lai on each PFT surface (m^2 m^{-2})   
    REAL(r_std), DIMENSION(npts)                              :: veget_cov_max_tree  !! Sum of veget_cov_max for the pft's which are trees
    INTEGER(i_std)                                            :: nbtree          !! Number of PFT's which are trees
    INTEGER(i_std)                                            :: i,j,k,m         !! indices (unitless)       
    REAL(r_std), DIMENSION(npts)                              :: soil_n_constraint !! Soil inorganic N pool (gN m-2) 
    REAL(r_std), DIMENSION(npts)                              :: ntake           !! total N associated to establishment (gN) 
    REAL(r_std), DIMENSION(npts)                              :: ntake_nitrate   !! nitrate associated to establishment (gN) 
!_ ================================================================================================================================

    IF (printlev>=3) WRITE(numout,*) 'Entering establish'

  !! 1. messages and initialization

    ! time during which young plants can be completely eaten by herbivores after germination
    ! (and then individual die) assume to be half year   
    ! No reference
    tau_eatup = one_year/2.
    
    ! Calculate the sum of the vegetation over the tree pft's and the number of pft's which are trees
    veget_cov_max_tree(:) = 0.0
    nbtree=0
    DO j = 1, nvm 
       IF (is_tree(j)) THEN
          veget_cov_max_tree(:) = veget_cov_max_tree(:) + veget_cov_max(:,j)
          nbtree = nbtree + 1
       END IF
    END DO
    ! Set nbtree=1 to avoid zero division later if there are no tree PFT's.
    ! For that case veget_cov_max_tree=0 so there will not be any impact.
    IF (nbtree == 0) nbtree=1

    !! 1.1 First call only
    IF ( firstcall_establish ) THEN

       WRITE(numout,*) 'establish:'

       WRITE(numout,*) '   > time during which a sapling can be entirely eaten by herbivores (d): ', &
            tau_eatup

       firstcall_establish = .FALSE.

    ENDIF

  !! 2. recalculate fpc

    IF (ok_dgvm) THEN
       fracnat(:) = un

       !! 2.1 Only natural part of the grid cell
       do j = 2,nvm ! Loop over # PFTs
          
          IF ( .NOT. natural(j) ) THEN
             fracnat(:) = fracnat(:) - veget_cov_max(:,j)
          ENDIF
       ENDDO ! Loop over # PFTs
       
       sumfpc(:) = zero

       !! 2.2 Total natural fpc on grid
       !      The overall fractional coverage of a PFT in a grid is calculated here.
       !      FPC is related to mean individual leaf area index by the Lambert-Beer law.
       !      See Eq. (1) in tex file.\n
       DO j = 2,nvm ! Loop over # PFTs
          IF ( natural(j) ) THEN
             WHERE(fracnat(:).GT.min_stomate)
                WHERE (lai(:,j) == val_exp) 
                   fpc_nat(:,j) = cn_ind(:,j) * ind(:,j) / fracnat(:)
                ELSEWHERE
                   fpc_nat(:,j) = cn_ind(:,j) * ind(:,j) / fracnat(:) & 
                        * ( un - exp( - biomass_to_lai(lm_lastyearmax(:,j),npts,j)* ext_coeff(j) ) )
                ENDWHERE
             ENDWHERE

             WHERE ( PFTpresent(:,j) )
                sumfpc(:) = sumfpc(:) + fpc_nat(:,j)
             ENDWHERE
          ELSE

             fpc_nat(:,j) = zero

          ENDIF

       ENDDO ! Loop over # PFTs
       
       !! 2.3 Total woody fpc on grid and number of regenerative tree pfts
       !      Total woody FPC increases by adding new FPC.
       !      Under the condition that temperature in last winter is higher than a threshold, 
       !      woody plants is exposed in higher competitive environment.
       sumfpc_wood(:) = zero

       DO j = 2,nvm ! Loop over # PFTs

          IF ( is_tree(j) .AND. natural(j) ) THEN

             ! total woody fpc
             WHERE ( PFTpresent(:,j) )
                sumfpc_wood(:) = sumfpc_wood(:) + fpc_nat(:,j)
             ENDWHERE

          ENDIF

       ENDDO ! Loop over # PFTs

       !! 2.4 Total number of natural grasses on grid\n
       !     Grass increment equals 'everywhere'\n
       spacefight_grass(:) = zero

       DO j = 2,nvm ! Loop over # PFTs

          IF ( .NOT. is_tree(j) .AND. natural(j) ) THEN

             ! Count a PFT fully only if it is present on a grid.
             WHERE ( PFTpresent(:,j) )
                spacefight_grass(:) = spacefight_grass(:) + everywhere(:,j)
             ENDWHERE

          ENDIF

       ENDDO ! Loop over # PFTs

       !! 2.5 Maximum establishment rate, based on climate only\n
       WHERE ( ( precip_annual(:) .GE. precip_crit ) .AND. ( gdd0(:) .GE. gdd_crit_estab ) )

          estab_rate_max_climate_tree(:) = estab_max_tree ! 'estab_max_*'; see 'stomate_constants.f90'
          estab_rate_max_climate_grass(:) = estab_max_grass

       ELSEWHERE

          estab_rate_max_climate_tree(:) = zero
          estab_rate_max_climate_grass(:) = zero

       ENDWHERE

       !! 2.6 Reduce maximum tree establishment rate if many trees are present.
       !      In the original DGVM, this is done using a step function which yields a
       !      reduction by factor 4 if sumfpc_wood(i) .GT.  fpc_crit - 0.05.
       !      This can lead to small oscillations (without consequences however).
       !      Here, a steady linear transition is used between fpc_crit-0.075 and
       !      fpc_crit-0.025.
       !      factor(:) = 1. - 15. * ( sumfpc_wood(:) - (fpc_crit-.075))
       !      factor(:) = MAX( 0.25_r_std, MIN( 1._r_std, factor(:)))
       !      S. Zaehle modified according to Smith et al. 2001
       !      See Eq. (2) in header
       factor(:)=(un - exp(- establish_scal_fact * (un - sumfpc_wood(:))))*(un - sumfpc_wood(:))
       estab_rate_max_tree(:) = estab_rate_max_climate_tree(:) * factor(:)

       !! 2.7 Modulate grass establishment rate.
       !      If canopy is not closed (fpc < fpc_crit-0.05), normal establishment.
       !      If canopy is closed, establishment is reduced by a factor 4.
       !      Factor is linear between these two bounds.
       !      This is different from the original DGVM where a step function is
       !      used at fpc_crit-0.05 (This can lead to small oscillations,
       !      without consequences however).
       !      factor(:) = 1. - 15. * ( sumfpc(:) - (fpc_crit-.05))
       !      factor(:) = MAX( 0.25_r_std, MIN( 1._r_std, factor(:)))
       !      estab_rate_max_grass(:) = estab_rate_max_climate_grass(:) * factor(:)
       !      S. Zaehle modified to true LPJ formulation, grasses are only allowed in the
       !      fpc fraction not occupied by trees..., 080806
       !      estab_rate_max_grass(:)=MAX(0.98-sumfpc(:),zero)
       !      See Eq. (3) in header
       estab_rate_max_grass(:) = MAX(MIN(estab_rate_max_climate_grass(:), max_tree_coverage - sumfpc(:)),zero)

       !! 2.8 Longterm grass NPP for competition between C4 and C3 grasses
       !      to avoid equal veget_cov_max, the idea is that more reestablishment
       !      is possible for the more productive PFT
       factor(:) = min_stomate

       DO j = 2,nvm ! Loop over # PFTs
          IF ( natural(j) .AND. .NOT.is_tree(j)) & 
               factor(:) = factor(:) + npp_longterm(:,j) * &
               biomass_to_lai(lm_lastyearmax(:,j),npts,j)
       ENDDO ! Loop over # PFTs

       !! 2.9 Establish natural PFTs
       d_ind(:,:) = zero

       IF ( NbNeighb /= 8 ) THEN
          CALL ipslerr_p(3, "establish", "This routine needs to be adapted to non rectengular grids", "Talk to Jan Polcher", " ")
       ENDIF

       DO j = 2,nvm ! Loop over # PFTs

          IF ( natural(j) ) THEN ! only for natural PFTs

             !! 2.9.1 PFT expansion across the grid box. Not to be confused with areal coverage.
             IF ( treat_expansion ) THEN

                ! only treat plants that are regenerative and present and still can expand
                DO i = 1, npts ! Loop over # pixels - domain size

                   IF ( PFTpresent(i,j) .AND. &
                        ( everywhere(i,j) .LT. un ) .AND. &
                        ( regenerate(i,j) .GT. regenerate_crit ) ) THEN

                      ! from how many sides is the grid box invaded (separate x and y directions
                      ! because resolution may be strongly anisotropic)
                      ! For the moment we only look into 4 direction but that can be expanded (JP) 
                      nfrontx = 0
                      IF ( neighbours(i,3) .GT. 0 ) THEN
                         IF ( everywhere(neighbours(i,3),j) .GT. 1.-min_stomate ) nfrontx = nfrontx+1
                      ENDIF
                      IF ( neighbours(i,7) .GT. 0 ) THEN
                         IF ( everywhere(neighbours(i,7),j) .GT. 1.-min_stomate ) nfrontx = nfrontx+1
                      ENDIF

                      nfronty = 0
                      IF ( neighbours(i,1) .GT. 0 ) THEN
                         IF ( everywhere(neighbours(i,1),j) .GT. 1.-min_stomate ) nfronty = nfronty+1
                      ENDIF
                      IF ( neighbours(i,5) .GT. 0 ) THEN
                         IF ( everywhere(neighbours(i,5),j) .GT. 1.-min_stomate ) nfronty = nfronty+1
                      ENDIF
                      
                      everywhere(i,j) = &
                           everywhere(i,j) + migrate(j) * dt/one_year * &
                           ( nfrontx / resolution(i,1) + nfronty / resolution(i,2) )
                      
                      IF ( .NOT. need_adjacent(i,j) ) THEN
                         
                         ! in that case, we also assume that the PFT expands from places within
                         ! the grid box (e.g., oasis).
                         ! What is this equation? No reference.
                         everywhere(i,j) = &
                              everywhere(i,j) + migrate(j) * dt/one_year * &
                              2. * SQRT( pi*everywhere(i,j)/(area(i)) )

                      ENDIF

                      everywhere(i,j) = MIN( everywhere(i,j), un )

                   ENDIF

                ENDDO ! Loop over # pixels - domain size

             ENDIF ! treat expansion?

             !! 2.9.2 Establishment rate
             !      - Is lower if the PFT is only present in a small part of the grid box
             !        (after its introduction), therefore multiplied by "everywhere".
             !      - Is divided by the number of PFTs that compete ("spacefight").
             !      - Is modulated by space availability (avail_tree, avail_grass).

             !! 2.9.2.1 present and regenerative trees
             IF ( is_tree(j) ) THEN

                WHERE ( PFTpresent(:,j) .AND. ( regenerate(:,j) .GT. regenerate_crit ) )
                   d_ind(:,j) = estab_rate_max_tree(:)*everywhere(:,j) * &
                        avail_tree(:) * dt/one_year
                ENDWHERE

             !! 2.9.2.2 present and regenerative grasses
             ELSE
                WHERE ( PFTpresent(:,j) .AND. ( regenerate(:,j) .GT. regenerate_crit )  & 
                     .AND.factor(:).GT.min_stomate .AND. spacefight_grass(:).GT. min_stomate) 
                   
                   d_ind(:,j) = estab_rate_max_grass(:)*everywhere(:,j)/spacefight_grass(:) * &
                        MAX(min_stomate,npp_longterm(:,j)*biomass_to_lai(lm_lastyearmax(:,j),npts,j)/factor(:)) * &
                        fracnat(:) * dt/one_year
                ENDWHERE
                
             ENDIF  ! tree/grass
             
          ENDIF ! if natural
       ENDDO  ! Loop over # PFTs
       
  !! 3. Lpj establishment in static case 

    !     Lpj establishment in static case, S. Zaehle 080806, account for real LPJ dynamics in
    !     prescribed vegetation, i.e. population dynamics within a given area of the grid cell.
    ELSE 

       d_ind(:,:) = zero

       DO j = 2,nvm ! Loop over # PFTs
          WHERE(ind(:,j)*cn_ind(:,j).GT.min_stomate)
             lai_ind(:) = biomass_to_lai(lm_lastyearmax(:,j),npts,j)/(ind(:,j)*cn_ind(:,j))
          ELSEWHERE
             lai_ind(:) = zero
          ENDWHERE

          !! 3.1 For natural woody PFTs
          IF ( natural(j) .AND. is_tree(j)) THEN 

             ! See Eq. (4) in tex file.            
             fpc_nat(:,j) =  MIN(un, cn_ind(:,j) * ind(:,j) * & 
                  MAX( ( un - exp( - ext_coeff(j) * lai_ind(:) ) ), min_cover ) )


             WHERE (veget_cov_max(:,j).GT.min_stomate.AND.ind(:,j).LE.2.)

                !! 3.1.1 Only establish into growing stands 
                !        Only establish into growing stands, ind can become very
                !        large in the static mode because LAI is very low in poor 
                !        growing conditions, favouring continuous establishment. 
                !        To avoid this a maximum IND is set. BLARPP: This should be
                !        replaced by a better stand density criteria.
                factor(:)=(un - exp(-establish_scal_fact * (un - fpc_nat(:,j))))*(un - fpc_nat(:,j))

                estab_rate_max_tree(:) = estab_max_tree * factor(:) * MAX(0.11 , nstress_season(:,j) )

                !! 3.1.2 do establishment for natural PFTs\n
                d_ind(:,j) = MAX( zero, estab_rate_max_tree(:) * dt/one_year)

             ENDWHERE

             !S. Zaehle: quickfix: to simulate even aged stand, uncomment the following lines...
             !where (ind(:,j) .LE. min_stomate)
             !d_ind(:,j) = 0.1 !MAX( 0.0, estab_rate_max_tree(:) * dt/one_year)
             WHERE (veget_cov_max(:,j).GT.min_stomate .AND. ind(:,j).EQ.zero)
                d_ind(:,j) = ind_0_estab
             ENDWHERE

          !! 3.2 For natural grass PFTs
          ELSEIF ( natural(j) .AND. .NOT.is_tree(j)) THEN 

             WHERE (veget_cov_max(:,j).GT.min_stomate)

                fpc_nat(:,j) =  cn_ind(:,j) * ind(:,j) * & 
                     MAX( ( un - exp( - ext_coeff(j) * lai_ind(:) ) ), min_cover )

                d_ind(:,j) = MAX(zero , (un - fpc_nat(:,j)) * dt/one_year )

             ENDWHERE

             WHERE (veget_cov_max(:,j).GT.min_stomate .AND. ind(:,j).EQ. zero)
                d_ind(:,j) = ind_0_estab 
             ENDWHERE

          ENDIF

       ENDDO ! Loop over # PFTs

    ENDIF ! DGVM OR NOT

  !! 4. Biomass calculation

    DO j = 2,nvm ! Loop over # PFTs

       IF ( natural(j) ) THEN ! only for natural PFTs

          !! 4.1 Herbivores reduce establishment rate
          !      We suppose that saplings are vulnerable during a given time after establishment.
          !      This is taken into account by preventively reducing the establishment rate.
          IF ( ok_herbivores ) THEN

             d_ind(:,j) = d_ind(:,j) * EXP( - tau_eatup/herbivores(:,j) )

          ENDIF

          !! 4.2 Total biomass.
          !      Add biomass only if d_ind, over one year, is of the order of ind.
          !      save old leaf mass to calculate leaf age
          leaf_mass_young(:) = leaf_frac(:,j,1) * biomass(:,j,ileaf,icarbon)

          ! total biomass of existing PFT to limit biomass added from establishment
          total_bm(:,:) = zero

          DO k = 1, nparts
             total_bm(:,:) = total_bm(:,:) + biomass(:,j,k,:)
          ENDDO
          ! soil nitrogen available to feed plant establishment.  
          ! This needs to be taken from soil N. Reason: plants must not produce nitrogen by establishment 
          ! (otherwise this creates N) 
          soil_n_constraint(:)=soil_n_min(:,j,iammonium)+soil_n_min(:,j,initrate) 


          IF(ok_dgvm) THEN
             vn(:) = veget_cov_max(:,j)
          ELSE
             vn(:) = un
          ENDIF

          !! 4.3 Woodmass calculation

          ! +++CHECK+++
          ! These calculations do not make use of circ_class. Most likley
          ! the results are not very precise if not simply wrong
          ! +++++++++++
          
          !! 4.3.1 with DGVM
          IF(ok_dgvm) THEN

             ! S. Zaehle calculate new woodmass_ind and veget_cov_max after establishment (needed for correct scaling!)
             ! essential correction for MERGE!
             IF(is_tree(j))THEN
                DO i=1,npts ! Loop over # pixels - domain size
                   IF((d_ind(i,j)+ind(i,j)).GT.min_stomate) THEN

                      IF((total_bm(i,icarbon).LE.min_stomate) .OR. (veget_cov_max(i,j) .LE. min_stomate)) THEN

                         ! new wood mass of PFT
                         woodmass_ind(i,j) = &
                              (((biomass(i,j,isapabove,icarbon) + biomass(i,j,isapbelow,icarbon) &
                              + biomass(i,j,iheartabove,icarbon) + biomass(i,j,iheartbelow,icarbon))*veget_cov_max(i,j)) &
                              + (bm_sapl_2D(i,j,isapabove,icarbon) + bm_sapl_2D(i,j,isapbelow,icarbon) &
                              + bm_sapl_2D(i,j,iheartabove,icarbon) + &
                                bm_sapl_2D(i,j,iheartbelow,icarbon))*d_ind(i,j))/(ind(i,j) + d_ind(i,j))

                      ELSE
 
                         ! new biomass is added to the labile pool, hence there is no change 
                         ! in CA associated with establishment
                         woodmass_ind(i,j) = &
                              & (biomass(i,j,isapabove,icarbon) + biomass(i,j,isapbelow,icarbon) &
                              & +biomass(i,j,iheartabove,icarbon) + biomass(i,j,iheartbelow,icarbon))*veget_cov_max(i,j) &
                              & /(ind(i,j) + d_ind(i,j))

                      ENDIF

                      ! new diameter of PFT
                      ! +++CHECK+++
                      ! Use functions from function_library to be consistent
                      dia(i) = (woodmass_ind(i,j)/(pipe_density(j)*pi/4.*pipe_tune2(i,j)))**(1./(2.+pipe_tune3(j)))
                      ! +++++++++++
                      ! Calculation of crown area has been changed
                      vn(i) = (ind(i,j) + d_ind(i,j))*crown_vertohor_dia(j)*crown_to_height(j)*pipe_tune2(i,j)*(dia(i)**pipe_tune3(j))

                   ENDIF
                ENDDO ! Loop over # pixels - domain size
             ELSE ! for grasses, cnd=1, so the above calculation cancels
                vn(:) = ind(:,j) + d_ind(:,j)
             ENDIF

          !! 4.3.2 without DGVM (static)\n
          ELSE 
             DO i=1,npts ! Loop over # pixels - domain size
                IF(is_tree(j).AND.(d_ind(i,j)+ind(i,j)).GT.min_stomate) THEN
                   IF(total_bm(i,icarbon).LE.min_stomate) THEN

                      ! new wood mass of PFT
                      woodmass_ind(i,j) = &
                           & (((biomass(i,j,isapabove,icarbon) + biomass(i,j,isapbelow,icarbon) &
                           & + biomass(i,j,iheartabove,icarbon) + biomass(i,j,iheartbelow,icarbon))) &
                           & + (bm_sapl_2D(i,j,isapabove,icarbon) + bm_sapl_2D(i,j,isapbelow,icarbon) &
                           & + bm_sapl_2D(i,j,iheartabove,icarbon) + &
                           &   bm_sapl_2D(i,j,iheartbelow,icarbon))*d_ind(i,j))/(ind(i,j)+d_ind(i,j))

                   ELSE
 
                      ! new biomass is added to the labile pool, hence there is no change 
                      ! in CA associated with establishment
                      woodmass_ind(i,j) = &
                           & (biomass(i,j,isapabove,icarbon) + biomass(i,j,isapbelow,icarbon) &
                           & + biomass(i,j,iheartabove,icarbon) + biomass(i,j,iheartbelow,icarbon)) &
                           & /(ind(i,j) + d_ind(i,j))

                   ENDIF
                ENDIF
             ENDDO ! Loop over # pixels - domain size

             vn(:) = un ! cannot change in static!, and veget_cov_max implicit in d_ind

          ENDIF

          !! 4.4 total biomass of PFT added by establishment defined over veget_cov_max ...

          total_bm_sapl(:,:) = zero
          total_bm_sapl_non(:,:) = zero
          biomass_old(:,j,:,:)=biomass(:,j,:,:)
          DO k = 1, nparts ! Loop over # litter tissues (nparts=8); see 'stomate_constants.f90'
             WHERE(d_ind(:,j).GT.min_stomate.AND.total_bm(:,icarbon).GT.min_stomate.AND.vn(:).GT.min_stomate)

                total_bm_sapl(:,icarbon) = total_bm_sapl(:,icarbon) + bm_sapl_2D(:,j,k,icarbon) * d_ind(:,j) / vn(:)
                total_bm_sapl(:,initrogen) = total_bm_sapl(:,initrogen) + bm_sapl_2D(:,j,k,initrogen) * d_ind(:,j)  / vn(:) 

                ! non-effective establishment
                total_bm_sapl_non(:,icarbon) = total_bm_sapl_non(:,icarbon) + &
                     bm_sapl_2D(:,j,k,icarbon) * (ind(:,j)+d_ind(:,j))*mortality(:,j) / vn(:)

                total_bm_sapl_non(:,initrogen) = total_bm_sapl_non(:,initrogen) + &
                     bm_sapl_2D(:,j,k,initrogen) * (ind(:,j)+d_ind(:,j))*mortality(:,j) / vn(:)
                

             ENDWHERE
          ENDDO ! Loop over # litter tissues


          total_bm_sapl(:,initrogen)=MIN(total_bm_sapl(:,initrogen),soil_n_constraint(:)) 
          ntake(:)=zero 
!Dan Zhu modification: there is a problem here, if DGVM is activated, co2_to_bm will never reach
!0 due to establishment, where d_ind is still large at equilibrium (=ind*mortality). So we
!need to subtract it from litter (not biomass, because the
!corresponding biomass has been lost in lpj_gap).

          !! 4.5 Update biomass at each component
          DO k = 1, nparts ! Loop over # litter tissues

             bm_new(:,:) = zero
             bm_non(:,:) = zero
             bm_eff(:,:) = zero

             ! first ever establishment, C flows
             WHERE( d_ind(:,j).GT.min_stomate .AND. &
                  total_bm(:,icarbon).LE.min_stomate .AND. &
                  veget_cov_max(:,j).GT.min_stomate)

                bm_new(:,icarbon) = d_ind(:,j) * bm_sapl_2D(:,j,k,icarbon) / vn(:)
                biomass(:,j,k,icarbon) = biomass(:,j,k,icarbon) + bm_new(:,icarbon)

                bm_new(:,initrogen) = d_ind(:,j) * (bm_sapl_2D(:,j,k,initrogen)) / vn(:)
                biomass(:,j,k,initrogen) = biomass(:,j,k,initrogen) + bm_new(:,initrogen)

                ! bm_to_litter minus the 'non-effective' establishment (mortality), but cannot be less than 0
                WHERE((veget_cov_max_tree(:) .GT. 0.1) .AND. (veget_cov_max(:,j) .LT. veget_cov_max_tree(:)/nbtree) )

                   bm_non(:,icarbon) = MIN( biomass(:,j,k,icarbon)+bm_to_litter(:,j,k,icarbon), &
                        (ind(:,j)+d_ind(:,j))*mortality(:,j) * bm_sapl_2D(:,j,k,icarbon)/vn(:) )
                   bm_eff(:,icarbon) = MIN( npp_longterm(:,j)/one_year, bm_new(:,icarbon)-bm_non(:,icarbon) )
                   bm_non(:,icarbon) = MIN( biomass(:,j,k,icarbon)+bm_to_litter(:,j,k,icarbon), &
                        bm_new(:,icarbon)-bm_eff(:,icarbon) )

                   bm_non(:,initrogen) = MIN( biomass(:,j,k,initrogen)+bm_to_litter(:,j,k,initrogen), &
                        (ind(:,j)+d_ind(:,j))*mortality(:,j) * bm_sapl_2D(:,j,k,initrogen)/vn(:) )
                   bm_eff(:,initrogen) = bm_new(:,initrogen)-bm_non(:,initrogen) 
                   bm_non(:,initrogen) = MIN( biomass(:,j,k,initrogen)+bm_to_litter(:,j,k,initrogen), &
                        bm_new(:,initrogen)-bm_eff(:,initrogen) )
                   ! I think there was an error 
                   ! co2_to_bm(:,j)=co2_to_bm(:,j) + bm_new(:) - bm_non(:)
                   co2_to_bm(:,j)=co2_to_bm(:,j) + bm_new(:,icarbon)/dt - bm_non(:,icarbon)/dt
                   ntake(:) = ntake(:) + bm_new(:,initrogen) - bm_non(:,initrogen)

                   WHERE( bm_to_litter(:,j,k,icarbon) .LT. bm_non(:,icarbon) )
                      biomass(:,j,k,icarbon) = biomass(:,j,k,icarbon) - &
                           ( bm_non(:,icarbon) - bm_to_litter(:,j,k,icarbon) )
                   ENDWHERE
                   bm_to_litter(:,j,k,icarbon) = bm_to_litter(:,j,k,icarbon) - &
                        MIN(bm_to_litter(:,j,k,icarbon), bm_non(:,icarbon) )

                   WHERE( bm_to_litter(:,j,k,initrogen) .LT. bm_non(:,initrogen) )
                      biomass(:,j,k,initrogen) = biomass(:,j,k,initrogen) - &
                           ( bm_non(:,initrogen) - bm_to_litter(:,j,k,initrogen) )
                   ENDWHERE
                   bm_to_litter(:,j,k,initrogen) = bm_to_litter(:,j,k,initrogen) - &
                        MIN(bm_to_litter(:,j,k,initrogen), bm_non(:,initrogen) )

                ELSEWHERE

                   bm_non(:,icarbon) = MIN( bm_to_litter(:,j,k,icarbon), &
                        (ind(:,j)+d_ind(:,j))*mortality(:,j) * bm_sapl_2D(:,j,k,icarbon)/vn(:) )
                   co2_to_bm(:,j)=co2_to_bm(:,j) + bm_new(:,icarbon)/dt - bm_non(:,icarbon)/dt
                   bm_to_litter(:,j,k,icarbon)=bm_to_litter(:,j,k,icarbon)- bm_non(:,icarbon)

                   bm_non(:,initrogen) = MIN( bm_to_litter(:,j,k,initrogen), &
                        (ind(:,j)+d_ind(:,j))*mortality(:,j) * bm_sapl_2D(:,j,k,initrogen)/vn(:) )
                   ntake(:)= ntake(:) + bm_new(:,initrogen) - bm_non(:,initrogen)
                   bm_to_litter(:,j,k,initrogen)=bm_to_litter(:,j,k,initrogen)- bm_non(:,initrogen)

                ENDWHERE

             ENDWHERE

             ! establishment into existing population, C flows
             WHERE(d_ind(:,j).GT.min_stomate.AND.total_bm(:,icarbon).GT.min_stomate)

                bm_new(:,icarbon) = total_bm_sapl(:,icarbon) * biomass_old(:,j,k,icarbon) / total_bm(:,icarbon)
                biomass(:,j,k,icarbon) = biomass(:,j,k,icarbon) + bm_new(:,icarbon)

                bm_new(:,initrogen) = total_bm_sapl(:,initrogen) * biomass_old(:,j,k,initrogen) / total_bm(:,initrogen)
                biomass(:,j,k,initrogen) = biomass(:,j,k,initrogen) + bm_new(:,initrogen)

                WHERE((veget_cov_max_tree(:) .GT. 0.1) .AND. (veget_cov_max(:,j) .LT. veget_cov_max_tree(:)/nbtree) )

                   bm_non(:,icarbon) = MIN( biomass(:,j,k,icarbon)+bm_to_litter(:,j,k,icarbon), &
                        total_bm_sapl_non(:,icarbon) *biomass_old(:,j,k,icarbon)/total_bm(:,icarbon) )
                   bm_eff(:,icarbon) = MIN( npp_longterm(:,j)/one_year, bm_new(:,icarbon)-bm_non(:,icarbon) )
                   bm_non(:,icarbon) = MAX( zero, MIN( biomass(:,j,k,icarbon)+bm_to_litter(:,j,k,icarbon)-min_stomate, &
                        bm_new(:,icarbon)-bm_eff(:,icarbon) ) )

                   bm_non(:,initrogen) = MIN( biomass(:,j,k,initrogen)+bm_to_litter(:,j,k,initrogen), &
                        total_bm_sapl_non(:,initrogen) *biomass_old(:,j,k,initrogen)/total_bm(:,initrogen) )
                   bm_eff(:,initrogen) = bm_new(:,initrogen)-bm_non(:,initrogen)
                   bm_non(:,initrogen) = MAX( zero, MIN( biomass(:,j,k,initrogen)+bm_to_litter(:,j,k,initrogen)-min_stomate, &
                        bm_new(:,initrogen)-bm_eff(:,initrogen) ) )

                   co2_to_bm(:,j)=co2_to_bm(:,j) + bm_new(:,icarbon) - bm_non(:,icarbon)
                   ntake(:) = ntake(:) + bm_new(:,initrogen) - bm_non(:,initrogen)

                   WHERE( bm_to_litter(:,j,k,icarbon) .LT. bm_non(:,icarbon) )
                      biomass(:,j,k,icarbon) = biomass(:,j,k,icarbon) - ( bm_non(:,icarbon) - bm_to_litter(:,j,k,icarbon) )
                   ENDWHERE
                   bm_to_litter(:,j,k,icarbon) = bm_to_litter(:,j,k,icarbon) - MIN(bm_to_litter(:,j,k,icarbon), bm_non(:,icarbon) )

                   WHERE( bm_to_litter(:,j,k,initrogen) .LT. bm_non(:,initrogen) )
                      biomass(:,j,k,initrogen) = biomass(:,j,k,initrogen) - ( bm_non(:,initrogen) - bm_to_litter(:,j,k,initrogen) )
                   ENDWHERE
                   bm_to_litter(:,j,k,initrogen) = bm_to_litter(:,j,k,initrogen) - &
                        MIN(bm_to_litter(:,j,k,initrogen), bm_non(:,initrogen) )

                ELSEWHERE

                   bm_non(:,icarbon) = MIN( bm_to_litter(:,j,k,icarbon), &
                        total_bm_sapl_non(:,icarbon) *biomass_old(:,j,k,icarbon)/total_bm(:,icarbon) )
                   co2_to_bm(:,j) = co2_to_bm(:,j) + bm_new(:,icarbon)/dt - bm_non(:,icarbon)/dt
                   bm_to_litter(:,j,k,icarbon)=bm_to_litter(:,j,k,icarbon)- bm_non(:,icarbon)

                   bm_non(:,initrogen) = MIN( bm_to_litter(:,j,k,initrogen), &
                        total_bm_sapl_non(:,initrogen) *biomass_old(:,j,k,initrogen)/total_bm(:,initrogen) )
                   ntake(:) = ntake(:) + bm_new(:,initrogen) - bm_non(:,initrogen)
                   bm_to_litter(:,j,k,initrogen)=bm_to_litter(:,j,k,initrogen)- bm_non(:,initrogen)
                ENDWHERE

             ENDWHERE

          ENDDO ! Loop over # litter tissues

          ntake_nitrate(:)=MIN(ntake(:),MAX(soil_n_min(:,j,iammonium)-min_stomate,0.0)) 
          soil_n_min(:,j,iammonium)=soil_n_min(:,j,iammonium)-ntake_nitrate(:) 
          n_uptake(:,j,iammonium)=n_uptake(:,j,iammonium) + ntake_nitrate(:)/dt 
  
          ntake(:)=MIN(ntake(:)-ntake_nitrate(:),MAX(soil_n_min(:,j,initrate)-min_stomate,0.0)) 
          soil_n_min(:,j,initrate)=soil_n_min(:,j,initrate)-ntake(:) 
          n_uptake(:,j,initrate)=n_uptake(:,j,initrate) + ntake(:)/dt 

          IF (ANY( bm_to_litter(:,j,:,icarbon) .LT. 0.0 ) .OR. ANY( biomass(:,j,:,icarbon) .LT. 0.0 ) ) THEN
             CALL ipslerr_p(3,'establish','something wrong in establish/gap.','','')
          ENDIF 

          !! 4.6 Decrease leaf age in youngest class if new leaf biomass is higher than old one.
          WHERE ( d_ind(:,j) * bm_sapl_2D(:,j,ileaf,icarbon) .GT. min_stomate )
 
             ! reset leaf ages. Should do a real calculation like in the npp routine, 
             ! but this case is rare and not worth messing around.
             ! S. Zaehle 080806, added real calculation now, because otherwise leaf_age/leaf_frac
             ! are not initialised for the calculation of vmax, and hence no growth at all.
             ! logic follows that of stomate_npp.f90, just that it's been adjusted for the code here
             leaf_age(:,j,1) = leaf_age(:,j,1) * leaf_mass_young(:) / &
                  ( leaf_mass_young(:) + d_ind(:,j) * bm_sapl_2D(:,j,ileaf,icarbon) )

          ENDWHERE

          leaf_mass_young(:) = leaf_mass_young(:) + d_ind(:,j) * bm_sapl_2D(:,j,ileaf,icarbon)   

          !! 4.7 Youngest class: new mass in youngest class divided by total new mass
          WHERE ( biomass(:,j,ileaf,icarbon) .GT. min_stomate )
             ! new age class fractions (fraction in youngest class increases)
             leaf_frac(:,j,1) = leaf_mass_young(:) / biomass(:,j,ileaf,icarbon)

          ENDWHERE

          !! 4.8 Other classes: old mass in leaf age class divided by new mass
          DO m = 2, nleafages

             WHERE ( biomass(:,j,ileaf,icarbon) .GT. min_stomate )

                leaf_frac(:,j,m) = leaf_frac(:,j,m) * & 
                     ( biomass(:,j,ileaf,icarbon) + d_ind(:,j) * bm_sapl_2D(:,j,ileaf,icarbon) ) /  biomass(:,j,ileaf,icarbon)

             ENDWHERE

          ENDDO

          !! 4.9 Update age and number of individuals
          WHERE ( d_ind(:,j) .GT. min_stomate )

             age(:,j) = age(:,j) * ind(:,j) / ( ind(:,j) + d_ind(:,j) )

             ind(:,j) = ind(:,j) + d_ind(:,j)

          ENDWHERE

          !! 4.10 Convert excess sapwood to heartwood
          !!      No longer done : supressed by S. Zaehle given that the LPJ logic of carbon allocation was 
          !!      contradictory to SLAVE allocation. See CVS tag 1_5 for initial formulation. 


       ENDIF ! natural

    ENDDO ! Loop over # PFTs

  !! 5. history

    d_ind = d_ind / dt

    CALL xios_orchidee_send_field("IND_ESTAB",d_ind)
    CALL xios_orchidee_send_field("ESTABTREE",estab_rate_max_tree)
    CALL xios_orchidee_send_field("ESTABGRASS",estab_rate_max_grass)

    CALL histwrite_p (hist_id_stomate, 'IND_ESTAB', itime, d_ind, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'ESTABTREE', itime, estab_rate_max_tree, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'ESTABGRASS', itime, estab_rate_max_grass, npts, hori_index)

    IF (printlev>=4) WRITE(numout,*) 'Leaving establish'

  END SUBROUTINE establish

END MODULE lpj_establish
