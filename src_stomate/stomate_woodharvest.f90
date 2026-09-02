! =================================================================================================================================
! MODULE       : stomate_woodharvest
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
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
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_stomate/stomate_woodharvest.f90 $
!! $Date: 2026-02-11 18:11:31 +0100 (mer. 11 févr. 2026) $
!! $Revision: 9376 $
!! \n
!_ ================================================================================================================================


MODULE stomate_woodharvest

  USE ioipsl_para
  USE stomate_data
  USE pft_parameters
  USE constantes
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC woodharvest_main
  
CONTAINS


!! ================================================================================================================================
!! SUBROUTINE   : woodharvest_main
!!
!>\BRIEF        Impact of wood harvest on carbon stocks
!!
!! DESCRIPTION  : This subroutine is always activate if DO_WOOD_HARVEST=y in the configuration file.
!! woodharvest_main is called from stomateLpj every year (at the beginning)
!! The impact of wood harvest on carbon stocks is computed in this subroutine. The land cover change is written
!! by the difference of current and previous "maximal" coverage fraction of a PFT. 
!! On the basis of this difference, the amount of 'new establishment'/'biomass export',
!! and increase/decrease of each component, are estimated.\n
!!
!! Main structure of woodharvest is:
!! 1. Initialization
!! 2. Calculation of changes in carbon stocks and biomass by wood harvest
!! 3. Update 10 year- and 100 year-turnover pool contents
!! 4. History
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : ::prod10, ::prod100, ::flux10, ::flux100,
!!   :: cflux_prod10 and :: cflux_prod100 
!!
!! REFERENCES   : None
!!
!! FLOWCHART    : 
!! \n
!_ ================================================================================================================================

  
  SUBROUTINE woodharvest_main ( npts, dt_days, veget_max_old, &
       biomass, &        
       flux10_harvest,flux100_harvest, prod10_harvest,prod100_harvest,&
       convflux_harvest,cflux_prod10_harvest,cflux_prod100_harvest,&
       harvestwood,harvestpft,fHarvestToProduct,nflux_prod_total_harvest)
    
    IMPLICIT NONE
    
  !! 0. Variable and parameter declaration 
    
    !! 0.1 Input variables
    
    INTEGER, INTENT(in)                            :: npts             !! Domain size - number of pixels (unitless)
    REAL(r_std), INTENT(in)                        :: dt_days          !! Time step of vegetation dynamics for stomate (days)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)   :: veget_max_old    !! Current "maximal" coverage fraction of a PFT (LAI
                                                                       !! -> infinity) on ground
    REAL(r_std), DIMENSION(npts), INTENT(in)       :: harvestwood      !! Harvested wood biomass (gC m-2 yr-1)
 
    !! 0.2 Output variables

    REAL(r_std), DIMENSION(npts), INTENT(out)      :: convflux_harvest        !! flux release during first year following wood harvest
    REAL(r_std), DIMENSION(npts), INTENT(out)      :: cflux_prod10_harvest    !! total annual release from the 10 year-turnover
                                                                              !! pool @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(npts), INTENT(out)      :: cflux_prod100_harvest   !! total annual release from the 100 year-
                                                                              !! turnover pool @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)  :: harvestpft              !! wood harvested for each PFT
                                                                              !! @tex ($gC m^{-2} day^{-1}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)  :: fHarvestToProduct       !!  Harvested biomass into product pool due to anthorpogenic 
                                                                              !! land use (gC m-2 )
    REAL(r_std), DIMENSION(npts), INTENT(out)      :: nflux_prod_total_harvest !! release of N associated to wood harvest  @tex ($gN m^{-2}$) @endtex   
!! 0.3 Modified variables   
    
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout) :: biomass    !! biomass @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(npts,0:10), INTENT(inout)  :: prod10_harvest  !! products remaining in the 10/100 year-turnover pool 
                                                                         !! after the annual release for each compartment
                                                                         !! (10 or 100 + 1 : input from year of land cover change)
    REAL(r_std), DIMENSION(npts,0:100), INTENT(inout) :: prod100_harvest !! see prod10_havest
    REAL(r_std), DIMENSION(npts,10), INTENT(inout)    :: flux10_harvest  !! annual release from the 10 year-turnover pool compartments

    REAL(r_std), DIMENSION(npts,100), INTENT(inout)   :: flux100_harvest !! annual release from the 100 year-turnover pool compartments
                                                
    !! 0.4 Local variables

    INTEGER(i_std)            :: i, j, l, m    !! indices (unitless)
    REAL(r_std)               :: harvest       !! wood harvest renomalized by forest fraction
    REAL(r_std)               :: relharvest    !! relative harvest fraction of each pft 
    REAL(r_std)               :: relbiomass    !! relative biomass of each PFT 
    REAL(r_std)               :: vegettree     !! vegfrac of trees
    REAL(r_std)               :: abovemean     !! mean aboveground biomass @tex ($gC m^{-2}$) @endtex
    REAL(r_std)               :: remove_biomass_heartabove_carbon  !! Removed carbon biomass from heartabove pool
    REAL(r_std)               :: remove_biomass_sapabove_carbon    !! Removed carbon biomass from sapabove pool
    REAL(r_std)               :: remove_biomass_heartabove_nitrogen!! Removed nitrogen biomass from heartabove pool
    REAL(r_std)               :: remove_biomass_sapabove_nitrogen  !! Removed nitrogen biomass from sapabove pool
    
!_ ================================================================================================================================

    IF (printlev>=3) WRITE(numout,*) 'Entering woodharvest_main'
    
  !! 1. initialization
    
    !! 2. calculation of changes in carbon stocks and biomass by wood harvest \n
    !! 2.1 initialization of carbon stocks\n
    prod10_harvest(:,0)           = zero
    prod100_harvest(:,0)          = zero   
    convflux_harvest(:)           = zero
    cflux_prod10_harvest(:)       = zero
    cflux_prod100_harvest(:)      = zero
    harvestpft(:,:)               = zero

    fHarvestToProduct(:,:)        = zero
    nflux_prod_total_harvest(:)   = zero

    DO i = 1, npts  
       abovemean=0          
       vegettree=0
       !! 2.3 calculation of total forest fraction and mean biomass
       DO j=1, nvm
          if (is_tree(j)) THEN
             abovemean=abovemean+((biomass(i,j,iheartabove,icarbon)+biomass(i,j,isapabove,icarbon)))*veget_max_old(i,j)
             vegettree=vegettree+veget_max_old(i,j)
          ENDIF
       ENDDO
       IF (vegettree > min_stomate) THEN
          
          harvest=harvestwood(i) 
          abovemean=abovemean/vegettree

          !! 2.4 wood harvest for each PFT is estimated proportionally to the relative biomass.
          DO j=1,nvm
             relharvest=0
             IF ((is_tree(j)) .AND. (abovemean > min_stomate)  .AND. &
                  (biomass(i,j,iheartabove,icarbon)+biomass(i,j,isapabove,icarbon)>min_stomate))  THEN
                relbiomass=((biomass(i,j,iheartabove,icarbon)+biomass(i,j,isapabove,icarbon))/abovemean) 
                relharvest=harvest*relbiomass

                !! 2.5 the harvest could not be higher than 20% of the biomasse to avoid rapid killing of vegetation
                IF ((biomass(i,j,iheartabove,icarbon)+biomass(i,j,isapabove,icarbon))*0.2<relharvest) &
                     relharvest=(biomass(i,j,iheartabove,icarbon)+biomass(i,j,isapabove,icarbon))*0.2; 
                
                remove_biomass_heartabove_carbon = relharvest*(biomass(i,j,iheartabove,icarbon) / &
                     (biomass(i,j,iheartabove,icarbon)+biomass(i,j,isapabove,icarbon)))
                remove_biomass_heartabove_nitrogen = remove_biomass_heartabove_carbon * biomass(i,j,iheartabove,initrogen) / &
                     biomass(i,j,iheartabove,icarbon)
                remove_biomass_sapabove_carbon = relharvest*(biomass(i,j,isapabove,icarbon) / &
                     (biomass(i,j,iheartabove,icarbon) + biomass(i,j,isapabove,icarbon)))
                remove_biomass_sapabove_nitrogen = remove_biomass_sapabove_carbon * biomass(i,j,isapabove,initrogen) / &
                     biomass(i,j,isapabove,icarbon)
                
                
                biomass(i,j,iheartabove,icarbon)= biomass(i,j,iheartabove,icarbon) - remove_biomass_heartabove_carbon
                biomass(i,j,isapabove,icarbon)= biomass(i,j,isapabove,icarbon) - remove_biomass_sapabove_carbon
                biomass(i,j,iheartabove,initrogen)= biomass(i,j,iheartabove,initrogen) - remove_biomass_heartabove_nitrogen
                biomass(i,j,isapabove,initrogen)= biomass(i,j,isapabove,initrogen) - remove_biomass_sapabove_nitrogen

                fHarvestToProduct(i,j)=relharvest*veget_max_old(i,j)
                ! 2.6 calculation of different flux like for LU
                convflux_harvest(i)  = convflux_harvest(i) +coeff_lcchange_1(j)*relharvest*veget_max_old(i,j)
                harvestpft(i,j)=relharvest*dt_days/one_year
                prod10_harvest(i,0)  = prod10_harvest(i,0)  +coeff_lcchange_10(j)*relharvest*veget_max_old(i,j)
                prod100_harvest(i,0) = prod100_harvest(i,0) +coeff_lcchange_100(j)*relharvest*veget_max_old(i,j)


                nflux_prod_total_harvest(i) = nflux_prod_total_harvest(i) + &
                     (remove_biomass_heartabove_nitrogen+remove_biomass_sapabove_nitrogen)*veget_max_old(i,j)
             ENDIF
          ENDDO
       ENDIF

       !! 3.4 update 10 year-turnover pool content following flux emission
       !!     (linear decay (10%) of the initial carbon input)
       DO  l = 0, 8
          m = 10 - l
          cflux_prod10_harvest(i) =  cflux_prod10_harvest(i) + flux10_harvest(i,m)
          prod10_harvest(i,m)   =  prod10_harvest(i,m-1)   - flux10_harvest(i,m-1)
          flux10_harvest(i,m)   =  flux10_harvest(i,m-1)
          IF (prod10_harvest(i,m) .LT. 1.0) prod10_harvest(i,m) = zero
       ENDDO
       
       cflux_prod10_harvest(i) = cflux_prod10_harvest(i) + flux10_harvest(i,1) 
       flux10_harvest(i,1)     = 0.1 * prod10_harvest(i,0)
       prod10_harvest(i,1)     = prod10_harvest(i,0)
       
       !! 3.5 update 100 year-turnover pool content following flux emission\n
       DO   l = 0, 98
          m = 100 - l
          cflux_prod100_harvest(i)  =  cflux_prod100_harvest(i) + flux100_harvest(i,m)
          prod100_harvest(i,m)      =  prod100_harvest(i,m-1)   - flux100_harvest(i,m-1)
          flux100_harvest(i,m)      =  flux100_harvest(i,m-1)
          IF (prod100_harvest(i,m).LT.1.0) prod100_harvest(i,m) = zero
       ENDDO
       
       cflux_prod100_harvest(i)  = cflux_prod100_harvest(i) + flux100_harvest(i,1) 
       flux100_harvest(i,1)      = 0.01 * prod100_harvest(i,0)
       prod100_harvest(i,1)      = prod100_harvest(i,0)
       prod10_harvest(i,0)        = 0.0
       prod100_harvest(i,0)       = 0.0 
       
    ENDDO ! Loop over # pixels - domain size
    
  !! 4. history
    convflux_harvest        = convflux_harvest/one_year*dt_days
    cflux_prod10_harvest    = cflux_prod10_harvest/one_year*dt_days
    cflux_prod100_harvest   = cflux_prod100_harvest/one_year*dt_days
    fHarvestToProduct       = fHarvestToProduct/one_year*dt_days    
    nflux_prod_total_harvest = nflux_prod_total_harvest/one_year*dt_days

    IF (printlev>=4) WRITE(numout,*) 'Leaving woodharvest_main'
    
  END SUBROUTINE woodharvest_main
  
END MODULE stomate_woodharvest
