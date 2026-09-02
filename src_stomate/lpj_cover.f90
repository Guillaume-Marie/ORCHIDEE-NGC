! =================================================================================================================================
! MODULE       : lpj_cover
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
!                This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        Recalculate vegetation cover and LAI
!!
!!\n DESCRIPTION : None
!!
!! RECENT CHANGE(S) : None
!!
!! REFERENCE(S) : 
!!        Sitch, S., B. Smith, et al. (2003), Evaluation of ecosystem dynamics,
!!        plant geography and terrestrial carbon cycling in the LPJ dynamic 
!!        global vegetation model, Global Change Biology, 9, 161-185.\n
!!        Smith, B., I. C. Prentice, et al. (2001), Representation of vegetation
!!        dynamics in the modelling of terrestrial ecosystems: comparing two
!!        contrasting approaches within European climate space,
!!        Global Ecology and Biogeography, 10, 621-637.\n
!!
!! SVN :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_stomate/lpj_cover.f90 $
!! $Date: 2024-02-28 07:37:25 +0100 (mer. 28 févr. 2024) $
!! $Revision: 8446 $
!! \n
!_ ================================================================================================================================

MODULE lpj_cover

  ! modules used:

  USE ioipsl_para
  USE stomate_data
  USE pft_parameters

  IMPLICIT NONE

  ! private & public routines

  PRIVATE
  PUBLIC cover

CONTAINS

!! ================================================================================================================================
!! SUBROUTINE     : cover
!!
!>\BRIEF          Recalculate vegetation cover and LAI
!!
!!\n DESCRIPTION : Veget_cov_max is first renewed here according to newly calculated foliage biomass in this calculation step 
!! Then, litter, soil carbon, and biomass are also recalcuted with taking into account the changes in Veget_cov_max (i.e. delta_veg)
!! Grid-scale fpc (foliage projected coverage) is calculated to obtain the shadede ground area by leaf's light capture
!! Finally, grid-scale fpc is adjusted not to exceed 1.0
!!
!! RECENT CHANGE(S) : 
!!
!!                    Early 2019.  Added soil nitrogen pools.
!!
!! MAIN OUTPUT VARIABLE(S) : ::lai (leaf area index, @tex $(m^2 m^{-2})$ @endtex), 
!! :: veget (fractional vegetation cover, unitless)
!!
!! REFERENCE(S)   : None
!! 
!! FLOWCHART :
!! \latexonly 
!!     \includegraphics[scale=0.5]{lpj_cover_flowchart.png}
!! \endlatexonly
!! \n
!_ ================================================================================================================================

  SUBROUTINE cover (npts, cn_ind, ind, biomass, &
       veget_cov_max, veget_cov_max_old, lai, litter, som, turnover_daily, bm_to_litter, &
       co2_to_bm, co2_fire, resp_hetero, resp_hetero_litter, resp_hetero_soil, resp_maint, resp_growth, resp_excess, gpp_daily,&
       deepSOM_a, deepSOM_s, deepSOM_p, &
       lignin_struc, lignin_wood, lignin_snag, soil_n_min)

!! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                                  :: npts             !! Domain size (unitless)  
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)                :: cn_ind           !! Crown area 
                                                                                    !! @tex $(m^2)$ @endtex per individual
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)                :: ind              !! Number of individuals 
                                                                                    !! @tex $(m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)                :: veget_cov_max_old!! "Maximal" coverage fraction of a PFT (LAI-> 
                                                                                    !! infinity) on ground at beginning of time 

    !! 0.2 Output variables

    !! 0.3 Modified variables

    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)             :: lai                      !! Leaf area index OF AN INDIVIDUAL PLANT 
                                                                                            !! @tex $(m^2 m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nlitt,nvm,nlevs,nelements), INTENT(inout) :: litter         !! Metabolic and structural litter, above and 
                                                                                            !! below ground @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,ncarb,nvm,nelements), INTENT(inout) :: som                  !! Carbon pool: active, slow, or passive @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout) :: biomass             !! Biomass @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)                  :: veget_cov_max       !! "Maximal" coverage fraction of a PFT (LAI->
                                                                                            !! infinity) on ground (unitless)
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout) :: turnover_daily      !! Turnover rates (gC m^{-2} day^{-1})
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout) :: bm_to_litter        !! Conversion of biomass to litter 
                                                                                            !! @tex $(gC m^{-2} day^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)                  :: co2_to_bm           !! biomass up take for establishment           
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)                  :: co2_fire            !! Carbon emitted to the atmosphere by fire(living
                                                                                            !! and dead biomass)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)                  :: resp_hetero         !! Heterotrophic respiration
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)                  :: resp_hetero_litter  !! Heterotrophic respiration from litter
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)                  :: resp_hetero_soil    !! Heterotrophic respiration from soil
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)                  :: resp_maint          !! Maintenance respiration 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)                  :: resp_growth         !! Growth respiration
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)                  :: resp_excess         !! Excess respiration
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)                  :: gpp_daily           !! Daily gross primary productivity
    REAL(r_std), DIMENSION(npts,ngrnd,nvm,nelements), INTENT(inout)  :: deepSOM_a           !! Soil carbon discretized with depth (g/m**3) active
    REAL(r_std), DIMENSION(npts,ngrnd,nvm,nelements), INTENT(inout)  :: deepSOM_s           !! Soil carbon discretized with depth (g/m**3) slow
    REAL(r_std), DIMENSION(npts,ngrnd,nvm,nelements), INTENT(inout)  :: deepSOM_p           !! Soil carbon discretized with depth (g/m**3) passive
    REAL(r_std), DIMENSION(npts,nvm,nlevs), INTENT(inout)            :: lignin_struc        !! ratio Lignine/Carbon in structural litter,
                                                                                            !! above and below ground
    REAL(r_std), DIMENSION(npts,nvm,nlevs), INTENT(inout)            :: lignin_wood         !! ratio Lignine/Carbon in woody litter,
                                                                                            !! above and below ground
    REAL(r_std), DIMENSION(npts,nvm,nlevs), INTENT(inout)            :: lignin_snag         !! ratio Lignine/Carbon in snag litter,
                                                                                            !! above and below ground
    REAL(r_std),DIMENSION(npts,nvm,nnspec), INTENT(inout)            :: soil_n_min          !! mineral nitrogen in the soil (gN/m**2)  
    !! 0.4 Local variables

    INTEGER(i_std)                                              :: i,j,k,m                  !! Index (unitless)
    REAL(r_std), DIMENSION(npts,nlitt,nlevs,nelements)          :: dilu_lit                 !! Dilution for carbon variables
    REAL(r_std), DIMENSION(npts,ncarb,nelements)                :: dilu_som                 !! Dilution for soil Organic Matter
                                                                                            !! @tex $(gC(or N) m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,ngrnd,nelements)                :: dilu_deepSOM_a           !! Dilution for active soil Organic Matter when ok_soil_carbon_discretization
                                                                                            !! @tex $(gC(or N) m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,ngrnd,nelements)                :: dilu_deepSOM_s           !! Dilution for slow soil Organic Matter when ok_soil_carbon_discretization
                                                                                            !! @tex $(gC(or N) m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,ngrnd,nelements)                :: dilu_deepSOM_p           !! Dilution for passive soil Organic Matter when ok_soil_carbon_discretization
                                                                                            !! @tex $(gC(or N) m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nnspec)                         :: dilu_sin                 !! Dilution for soil Inorganic Nitrogen
                                                                                            !! @tex $(gN m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nlevs)                          :: dilu_lf_struc            !! Dilution for fraction of structural litter that is lignin
    REAL(r_std), DIMENSION(npts,nlevs)                          :: dilu_lf_wood             !! Dilution for fraction of woody litter that is lignin
    REAL(r_std), DIMENSION(npts,nlevs)                          :: dilu_lf_snag             !! Dilution for fraction of snag litter that is lignin
    REAL(r_std), DIMENSION(npts,nparts,nelements)               :: dilu_bio                 !! Dilution for carbon variables
    REAL(r_std), DIMENSION(npts)                                :: dilu_TCarbon             !! Dilution for carbon variables
    REAL(r_std), DIMENSION(npts,nparts,nelements)               :: dilu_turnover_daily      !! Dilution for carbon variables
    REAL(r_std), DIMENSION(npts,nparts,nelements)               :: dilu_bm_to_litter        !! Dilution for carbon variables
    REAL(r_std), DIMENSION(npts)                                :: dilu_co2flux_new         !! Dilution for carbon variables
    REAL(r_std), DIMENSION(npts)                                :: dilu_gpp_daily           !! Dilution for carbon variables
    REAL(r_std), DIMENSION(npts)                                :: dilu_resp_growth         !! Dilution for carbon variables
    REAL(r_std), DIMENSION(npts)                                :: dilu_resp_excess         !! Dilution for carbon variables
    REAL(r_std), DIMENSION(npts)                                :: dilu_resp_maint          !! Dilution for carbon variables
    REAL(r_std), DIMENSION(npts)                                :: dilu_resp_hetero         !! Dilution for carbon variables
    REAL(r_std), DIMENSION(npts)                                :: dilu_resp_hetero_litter  !! Dilution for carbon variables
    REAL(r_std), DIMENSION(npts)                                :: dilu_resp_hetero_soil    !! Dilution for carbon variables
    REAL(r_std), DIMENSION(npts)                                :: dilu_co2_to_bm           !! Dilution for carbon variables
    REAL(r_std), DIMENSION(npts)                                :: dilu_co2_fire            !! Dilution for carbon variables
    REAL(r_std), DIMENSION(npts,nvm)                            :: TCarbon                  !! Total carbon
    REAL(r_std), DIMENSION(npts,ncarb,nvm,nelements)            :: som_old                  !! Carbon pool stored before update: active, slow, or passive @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,ngrnd,nvm,nelements)            :: deepSOM_a_old            !! Soil carbon discretized with depth (g/m**3) active stored before update
    REAL(r_std), DIMENSION(npts,ngrnd,nvm,nelements)            :: deepSOM_s_old            !! Soil carbon discretized with depth (g/m**3) slow stored before update
    REAL(r_std), DIMENSION(npts,ngrnd,nvm,nelements)            :: deepSOM_p_old            !! Soil carbon discretized with depth (g/m**3) passive stored before update
    REAL(r_std), DIMENSION(npts,nvm)                            :: co2flux_new              !! NBP after re-calculation in order to conserve carbon 
    REAL(r_std), DIMENSION(npts,nvm)                            :: co2flux_old              !! NBP before re-calculation 
    REAL(r_std), DIMENSION(nvm)                                 :: delta_veg                !! Conversion factors (unitless)
    REAL(r_std), DIMENSION(nvm)                                 :: reduct                   !! Conversion factors (unitless)
    REAL(r_std)                                                 :: delta_veg_sum            !! Conversion factors (unitless)
    REAL(r_std)                                                 :: diff                     !! Conversion factors (unitless)
    REAL(r_std)                                                 :: sr                       !! Conversion factors (unitless)
    REAL(r_std), DIMENSION(npts)                                :: frac_nat                 !! Conversion factors (unitless)
    REAL(r_std), DIMENSION(npts)                                :: sum_vegettree            !! Conversion factors (unitless)
    REAL(r_std), DIMENSION(npts)                                :: sum_vegetgrass           !! Conversion factors (unitless) 
    REAL(r_std), DIMENSION(npts)                                :: sum_veget_natveg         !! Conversion factors (unitless)
    REAL(r_std), DIMENSION(npts)                                :: vartmp                   !! Temporary variable used to add history

!_ ================================================================================================================================

 !! 1. If the vegetation is dynamic, calculate new maximum vegetation cover for natural plants
  
    IF ( ok_dgvm ) THEN

       !! 1.1  Calculate initial values of vegetation cover
       frac_nat(:) = un
       sum_veget_natveg(:) = zero
       veget_cov_max(:,ibare_sechiba) = un
       som_old(:,:,:,:) = som(:,:,:,:)
       deepSOM_a_old(:,:,:,:) = deepSOM_a(:,:,:,:) 
       deepSOM_s_old(:,:,:,:) = deepSOM_s(:,:,:,:)
       deepSOM_p_old(:,:,:,:) = deepSOM_p(:,:,:,:)

       DO j = 2,nvm ! loop over PFTs

          IF ( natural(j) ) THEN
	     
             ! Summation of individual tree crown area to get total foliar projected coverage
             veget_cov_max(:,j) = ind(:,j) * cn_ind(:,j)
             sum_veget_natveg(:) = sum_veget_natveg(:) + veget_cov_max(:,j)

          ELSE
             
             !fraction occupied by agriculture needs to be substracted for the DGVM
             !this is used below to constrain veget for natural vegetation, see below
             frac_nat(:) = frac_nat(:) - veget_cov_max(:,j)

          ENDIF

       ENDDO ! loop over PFTs

       DO i = 1, npts ! loop over grid points
	  
          ! Recalculation of vegetation projected coverage when ::frac_nat was below ::sum_veget_natveg
          ! It means that non-natural vegetation will recover ::veget_cov_max as natural vegetation
          IF (sum_veget_natveg(i) .GT. frac_nat(i) .AND. frac_nat(i) .GT. min_stomate) THEN

             DO j = 2,nvm ! loop over PFTs
                IF( natural(j) ) THEN
                   veget_cov_max(i,j) =  veget_cov_max(i,j) * frac_nat(i) / sum_veget_natveg(i)
                ENDIF
             ENDDO ! loop over PFTs

          ENDIF
       ENDDO ! loop over grid points
	
       ! Renew veget_cov_max of bare soil as 0 to difference of veget_cov_max (ibare_sechiba) 
       ! to current veget_cov_max
       DO j = 2,nvm ! loop over PFTs
          veget_cov_max(:,ibare_sechiba) = veget_cov_max(:,ibare_sechiba) - veget_cov_max(:,j)
       ENDDO ! loop over PFTs
       veget_cov_max(:,ibare_sechiba) = MAX( veget_cov_max(:,ibare_sechiba), zero )

       !! 1.2 Calculate carbon fluxes between PFTs to maintain mass balance
       !! Assure carbon closure when veget_cov_max changes(delta_veg): if veget_cov_max of some PFTs decrease, we use "dilu" to 
       !! record the corresponding lost in carbon (biomass, litter, soil carbon, gpp, respiration etc.) for 
       !! these PFTs, and re-allocate "dilu" to those PFTs with increasing veget_cov_max.
       DO i = 1, npts ! loop over grid points
          
          ! Calculate the change in veget_cov_max between previous time step and current time step
          delta_veg(:) = veget_cov_max(i,:)-veget_cov_max_old(i,:)
          delta_veg_sum = SUM(delta_veg,MASK=delta_veg.LT.zero)

          dilu_lit(i,:,:,:) = zero
          dilu_som(i,:,:) = zero
          dilu_deepSOM_a(i,:,:) = zero
          dilu_deepSOM_s(i,:,:) = zero
          dilu_deepSOM_p(i,:,:) = zero
          dilu_sin(i,:) = zero
          dilu_lf_struc(i,:) = zero
          dilu_lf_wood(i,:) = zero
          dilu_lf_snag(i,:) = zero
          dilu_bio(i,:,:) = zero
          dilu_TCarbon(i)=zero
          dilu_turnover_daily(i,:,:)=zero
          dilu_bm_to_litter(i,:,:)=zero
          dilu_co2flux_new(i)=zero
          dilu_gpp_daily(i)=zero
          dilu_resp_growth(i)=zero
          dilu_resp_excess(i)=zero
          dilu_resp_maint(i)=zero
          dilu_resp_hetero(i)=zero
          dilu_resp_hetero_litter(i)=zero
          dilu_resp_hetero_soil(i)=zero
          dilu_co2_to_bm(i)=zero
          dilu_co2_fire(i)=zero

          ! Calculate TCarbon: total carbon including biomass, litter and soil carbon, as well as "today's" turnover and 
          ! bm_to_litter due to mortality, because today's turnover and bm_to_litter are not yet added into "litter" until tomorrow. 
          DO j=1, nvm
                TCarbon(i,j)=SUM(biomass(i,j,:,icarbon))+SUM(som(i,:,j,icarbon))+SUM(litter(i,:,j,:,icarbon))+&
                     SUM(turnover_daily(i,j,:,icarbon))+SUM(bm_to_litter(i,j,:,icarbon))
                co2flux_old(i,j)=resp_maint(i,j)+resp_growth(i,j)+resp_excess(i,j)+resp_hetero(i,j) + &
                     co2_fire(i,j)-co2_to_bm(i,j)-gpp_daily(i,j)
                co2flux_new(i,j)=resp_maint(i,j)+resp_growth(i,j)+resp_excess(i,j)+resp_hetero(i,j) + &
                     co2_fire(i,j)-co2_to_bm(i,j)-gpp_daily(i,j)
          ENDDO

          DO j=1, nvm ! loop over PFTs
             IF ( delta_veg(j) < -min_stomate ) THEN 
                dilu_lit(i,:,:,:) =  dilu_lit(i,:,:,:) + delta_veg(j) * litter(i,:,j,:,:) / delta_veg_sum
                dilu_som(i,:,:) =  dilu_som(i,:,:) + delta_veg(j) * som(i,:,j,:) / delta_veg_sum
                dilu_deepSOM_a(i,:,:) =  dilu_deepSOM_a(i,:,:) + delta_veg(j) * deepSOM_a(i,:,j,:) / delta_veg_sum
                dilu_deepSOM_s(i,:,:) =  dilu_deepSOM_s(i,:,:) + delta_veg(j) * deepSOM_s(i,:,j,:) / delta_veg_sum
                dilu_deepSOM_p(i,:,:) =  dilu_deepSOM_p(i,:,:) + delta_veg(j) * deepSOM_p(i,:,j,:) / delta_veg_sum
                dilu_sin(i,:) =  dilu_sin(i,:) + delta_veg(j) * soil_n_min(i,j,:) / delta_veg_sum 
                dilu_lf_struc(i,:) = dilu_lf_struc(i,:) + & 
                     delta_veg(j) * lignin_struc(i,j,:) * litter(i,istructural,j,:,icarbon)/ delta_veg_sum
                dilu_lf_wood(i,:) = dilu_lf_wood(i,:) + &
                     delta_veg(j) * lignin_wood(i,j,:)*litter(i,iwoody,j,:,icarbon)  / delta_veg_sum
                dilu_lf_snag(i,:) = dilu_lf_snag(i,:) + &
                     delta_veg(j) * lignin_snag(i,j,:)*litter(i,isnag,j,:,icarbon)  / delta_veg_sum
                dilu_TCarbon(i)= dilu_TCarbon(i) + delta_veg(j) * TCarbon(i,j) / delta_veg_sum
                dilu_turnover_daily(i,:,:)=dilu_turnover_daily(i,:,:)+delta_veg(j)*turnover_daily(i,j,:,:)/delta_veg_sum
                dilu_bm_to_litter(i,:,:)=dilu_bm_to_litter(i,:,:)+delta_veg(j)*bm_to_litter(i,j,:,:)/delta_veg_sum
                dilu_co2flux_new(i)=dilu_co2flux_new(i)+delta_veg(j)*co2flux_old(i,j)/delta_veg_sum
                dilu_gpp_daily(i)=dilu_gpp_daily(i)+delta_veg(j)*gpp_daily(i,j)/delta_veg_sum
                dilu_resp_growth(i)=dilu_resp_growth(i)+delta_veg(j)*resp_growth(i,j)/delta_veg_sum
                dilu_resp_excess(i)=dilu_resp_excess(i)+delta_veg(j)*resp_excess(i,j)/delta_veg_sum
                dilu_resp_maint(i)=dilu_resp_maint(i)+delta_veg(j)*resp_maint(i,j)/delta_veg_sum
                dilu_resp_hetero(i)=dilu_resp_hetero(i)+delta_veg(j)*resp_hetero(i,j)/delta_veg_sum
                dilu_resp_hetero_litter(i)=dilu_resp_hetero_litter(i)+delta_veg(j)*resp_hetero_litter(i,j)/delta_veg_sum
                dilu_resp_hetero_soil(i)=dilu_resp_hetero_soil(i)+delta_veg(j)*resp_hetero_soil(i,j)/delta_veg_sum
                dilu_co2_to_bm(i)=dilu_co2_to_bm(i)+delta_veg(j)*co2_to_bm(i,j)/delta_veg_sum
                dilu_co2_fire(i)=dilu_co2_fire(i)+delta_veg(j)*co2_fire(i,j)/delta_veg_sum
             ENDIF
          ENDDO ! loop over PFTs

          DO j=1, nvm ! loop over PFTs
             IF ( delta_veg(j) > min_stomate) THEN

                ! Dilution of reservoirs
                ! Recalculate the litter and soil carbon with taking into accout the change in 
                ! veget_cov_max (delta_veg)
                ! Litter
                ! Lignin fraction of structural litter
                lignin_struc(i,j,:)=(lignin_struc(i,j,:) * veget_cov_max_old(i,j)* litter(i,istructural,j,:,icarbon) + & 
                     dilu_lf_struc(i,:) * delta_veg(j)) / veget_cov_max(i,j) 

                ! Lignin fraction of snag litter
                lignin_snag(i,j,:)=(lignin_snag(i,j,:) * veget_cov_max_old(i,j)* litter(i,isnag,j,:,icarbon) + & 
                     dilu_lf_snag(i,:) * delta_veg(j)) / veget_cov_max(i,j)

                ! Lignin fraction of woody litter
                lignin_wood(i,j,:)=(lignin_wood(i,j,:) * veget_cov_max_old(i,j)* litter(i,iwoody,j,:,icarbon) + & 
                     dilu_lf_wood(i,:) * delta_veg(j)) / veget_cov_max(i,j)

                litter(i,:,j,:,:)=(litter(i,:,j,:,:) * veget_cov_max_old(i,j) + dilu_lit(i,:,:,:) * delta_veg(j)) &
                                  / veget_cov_max(i,j)

                WHERE ( litter(i,istructural,j,:,icarbon) > min_stomate )
                   lignin_struc(i,j,:) = lignin_struc(i,j,:)/litter(i,istructural,j,:,icarbon)
                ELSEWHERE
                   lignin_struc(i,j,:) = LC_leaf(j)
                ENDWHERE

                WHERE ( litter(i,iwoody,j,:,icarbon) > min_stomate )
                   lignin_wood(i,j,:) = lignin_struc(i,j,:)/litter(i,iwoody,j,:,icarbon)
                ELSEWHERE
                   lignin_wood(i,j,:) = LC_heartabove(j)
                ENDWHERE

                WHERE ( litter(i,isnag,j,:,icarbon) > min_stomate )
                   lignin_snag(i,j,:) = lignin_struc(i,j,:)/litter(i,isnag,j,:,icarbon)
                ELSEWHERE
                   lignin_snag(i,j,:) = LC_heartabove(j)
                ENDWHERE

                ! Soil carbon
                som(i,:,j,:)=(som(i,:,j,:) * veget_cov_max_old(i,j) + dilu_som(i,:,:) * delta_veg(j)) / veget_cov_max(i,j)
                ! BG: Here only done for carbon but I keep as much as possible the possibility to do it with N
                IF ( ok_soil_carbon_discretization ) THEN
                  deepSOM_a(i,:,j,:)=(deepSOM_a(i,:,j,:) * veget_cov_max_old(i,j) + &
                       dilu_deepSOM_a(i,:,:) * delta_veg(j)) / veget_cov_max(i,j)
                  deepSOM_s(i,:,j,:)=(deepSOM_s(i,:,j,:) * veget_cov_max_old(i,j) + &
                       dilu_deepSOM_s(i,:,:) * delta_veg(j)) / veget_cov_max(i,j)
                  deepSOM_p(i,:,j,:)=(deepSOM_p(i,:,j,:) * veget_cov_max_old(i,j) + &
                       dilu_deepSOM_p(i,:,:) * delta_veg(j)) / veget_cov_max(i,j)
                ENDIF

                ! Soil inorganic nitrogen  
                soil_n_min(i,j,:)=(soil_n_min(i,j,:) * veget_cov_max_old(i,j) + dilu_sin(i,:) &  
                     * delta_veg(j)) / veget_cov_max(i,j) 
                TCarbon(i,j)=(TCarbon(i,j) * veget_cov_max_old(i,j) + dilu_TCarbon(i) * delta_veg(j)) / veget_cov_max(i,j)
                turnover_daily(i,j,:,:)=(turnover_daily(i,j,:,:)*veget_cov_max_old(i,j)+&
                     dilu_turnover_daily(i,:,:)*delta_veg(j))/veget_cov_max(i,j)
                bm_to_litter(i,j,:,:)=(bm_to_litter(i,j,:,:)*veget_cov_max_old(i,j)+&
                     dilu_bm_to_litter(i,:,:)*delta_veg(j))/veget_cov_max(i,j)
                co2flux_new(i,j)=(co2flux_old(i,j)*veget_cov_max_old(i,j)+dilu_co2flux_new(i)*delta_veg(j))/veget_cov_max(i,j)
                gpp_daily(i,j)=(gpp_daily(i,j)*veget_cov_max_old(i,j)+dilu_gpp_daily(i)*delta_veg(j))/veget_cov_max(i,j)
                resp_growth(i,j)=(resp_growth(i,j)*veget_cov_max_old(i,j)+dilu_resp_growth(i)*delta_veg(j))/veget_cov_max(i,j)
                resp_excess(i,j)=(resp_excess(i,j)*veget_cov_max_old(i,j)+dilu_resp_excess(i)*delta_veg(j))/veget_cov_max(i,j)
                resp_maint(i,j)=(resp_maint(i,j)*veget_cov_max_old(i,j)+dilu_resp_maint(i)*delta_veg(j))/veget_cov_max(i,j)
                resp_hetero(i,j)=(resp_hetero(i,j)*veget_cov_max_old(i,j)+dilu_resp_hetero(i)*delta_veg(j))/veget_cov_max(i,j)
                resp_hetero_litter(i,j)=(resp_hetero_litter(i,j)*veget_cov_max_old(i,j)+ &
                              dilu_resp_hetero_litter(i)*delta_veg(j))/veget_cov_max(i,j)
                resp_hetero_soil(i,j)=(resp_hetero_soil(i,j)*veget_cov_max_old(i,j)+ &
                              dilu_resp_hetero_soil(i)*delta_veg(j))/veget_cov_max(i,j)
                co2_to_bm(i,j)=(co2_to_bm(i,j)*veget_cov_max_old(i,j)+dilu_co2_to_bm(i)*delta_veg(j))/veget_cov_max(i,j)
                co2_fire(i,j)=(co2_fire(i,j)*veget_cov_max_old(i,j)+dilu_co2_fire(i)*delta_veg(j))/veget_cov_max(i,j)
             ENDIF

             IF(veget_cov_max(i,j).GT.min_stomate) THEN

                ! Correct biomass densities to conserve mass
                ! since it's defined on veget_cov_max
                biomass(i,j,:,:) = biomass(i,j,:,:) * veget_cov_max_old(i,j) / veget_cov_max(i,j)

             ENDIF

          ENDDO ! loop over PFTs
       ENDDO ! loop over grid points

       vartmp(:)=SUM(gpp_daily*veget_cov_max,dim=2)
       CALL histwrite_p (hist_id_stomate, "tGPP", itime, vartmp, npts, hori_index)
       vartmp(:)=SUM(resp_growth*veget_cov_max,dim=2)
       CALL histwrite_p (hist_id_stomate, "tRESP_GROWTH", itime, vartmp, npts, hori_index)
       vartmp(:)=SUM(resp_excess*veget_cov_max,dim=2)
       CALL histwrite_p (hist_id_stomate, "tRESP_EXCESS", itime, vartmp, npts, hori_index)
       vartmp(:)=SUM(resp_maint*veget_cov_max,dim=2)
       CALL histwrite_p (hist_id_stomate, "tRESP_MAINT", itime, vartmp, npts, hori_index)
       vartmp(:)=SUM(resp_hetero*veget_cov_max,dim=2)
       CALL histwrite_p (hist_id_stomate, "tRESP_HETERO", itime, vartmp, npts, hori_index)
       vartmp(:)=SUM(co2_to_bm*veget_cov_max,dim=2)
       CALL histwrite_p (hist_id_stomate, "tCO2_TAKEN", itime, vartmp, npts, hori_index)
       vartmp(:)=SUM(co2_fire*veget_cov_max,dim=2)
       CALL histwrite_p (hist_id_stomate, "tCO2_FIRE", itime, vartmp, npts, hori_index)
       vartmp(:)=SUM(co2flux_new*veget_cov_max,dim=2)
       CALL histwrite_p (hist_id_stomate, "tCO2FLUX", itime, vartmp, npts, hori_index)
       vartmp(:)=SUM(co2flux_old*veget_cov_max_old,dim=2)
       CALL histwrite_p (hist_id_stomate, "tCO2FLUX_OLD", itime, vartmp, npts, hori_index)
       vartmp(:)=SUM(TCarbon*veget_cov_max,dim=2)
       CALL histwrite_p (hist_id_stomate, "tCARBON", itime, vartmp, npts, hori_index)
       vartmp(:)=SUM(SUM(biomass(:,:,:,icarbon),dim=3)*veget_cov_max,dim=2)
       CALL histwrite_p (hist_id_stomate, "tBIOMASS", itime, vartmp, npts, hori_index)
       vartmp(:)=SUM(SUM(SUM(litter(:,:,:,:,icarbon),dim=4),dim=2)*veget_cov_max,dim=2)
       CALL histwrite_p (hist_id_stomate, "tLITTER", itime, vartmp, npts, hori_index)
       vartmp(:)=SUM(SUM(som(:,:,:,icarbon),dim=2)*veget_cov_max,dim=2)
       CALL histwrite_p (hist_id_stomate, "tSOILC", itime, vartmp, npts, hori_index)

      IF ( ok_soil_carbon_discretization ) THEN
        vartmp(:)=SUM(SUM(deepSOM_a(:,:,:,icarbon),dim=2)*veget_cov_max,dim=2)
        CALL histwrite_p (hist_id_stomate, "tDEEPCa", itime, vartmp, npts, hori_index)
        vartmp(:)=SUM(SUM(deepSOM_s(:,:,:,icarbon),dim=2)*veget_cov_max,dim=2)
        CALL histwrite_p (hist_id_stomate, "tDEEPCs", itime, vartmp, npts, hori_index)
        vartmp(:)=SUM(SUM(deepSOM_p(:,:,:,icarbon),dim=2)*veget_cov_max,dim=2)
        CALL histwrite_p (hist_id_stomate, "tDEEPCp", itime, vartmp, npts, hori_index)
        vartmp(:)=SUM(SUM(deepSOM_a(:,:,:,initrogen),dim=2)*veget_cov_max,dim=2)
        CALL histwrite_p (hist_id_stomate, "tDEEPNa", itime, vartmp, npts, hori_index)
        vartmp(:)=SUM(SUM(deepSOM_s(:,:,:,initrogen),dim=2)*veget_cov_max,dim=2)
        CALL histwrite_p (hist_id_stomate, "tDEEPNs", itime, vartmp, npts, hori_index)
        vartmp(:)=SUM(SUM(deepSOM_p(:,:,:,initrogen),dim=2)*veget_cov_max,dim=2)
        CALL histwrite_p (hist_id_stomate, "tDEEPNp", itime, vartmp, npts, hori_index)
      ENDIF

    ENDIF

  END SUBROUTINE cover

END MODULE lpj_cover
