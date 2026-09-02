! =================================================================================================================================
! MODULE       : sapiens_product_use
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       Following lateral transport, calculate C-stored in the product pools
!!             and the CO2-release from product decomposition  
!!
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S)	: None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_stomate/sapiens_product_use.f90 $
!! $Date: 2025-09-08 14:34:13 +0200 (lun. 08 sept. 2025) $
!! $Revision: 9087 $
!! \n
!_ ================================================================================================================================

MODULE sapiens_product_use

! modules used:
  
  USE stomate_data
  USE pft_parameters
  USE constantes
  USE grid
  USE time, ONLY : ts_annual_proc, one_year
  USE function_library,    ONLY: check_mass_balance, get_printlev
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC product_add, product_decomp, product_init

  INTEGER(i_std), SAVE       :: printlev_loc                             !! Local level of text output for this module
!$OMP THREADPRIVATE(printlev_loc)
  LOGICAL, SAVE              :: firstcall_sapiens_product_add = .TRUE.   !! first call flag
!$OMP THREADPRIVATE(firstcall_sapiens_product_add)
  LOGICAL, SAVE              :: firstcall_sapiens_product_decomp = .TRUE.   !! first call flag
!$OMP THREADPRIVATE(firstcall_sapiens_product_decomp)

CONTAINS 


!! ================================================================================================================================
!! SUBROUTINE   : product_add
!!
!>\BRIEF        Calculate the C stored in biomass-based products.
!!
!! DESCRIPTION  : The basic approach distinguished 3 pools (as introduced by
!! Shilong Piao based on Houghton. Wood harvest with small dimensions is being 
!! used as fire wood the remaining wood goes into the short, medium and long 
!! product pools. This implies that when the dimensions of the harvest change, 
!! the product pools will change as well.\n
!! The basic approach distinguished 3 pools (as introduced by
!! Shilong Piao based on Houghton. This module does NOT make use of the 
!! dimensions of the woody harvest or of the management and cutting flags.
!! All wood is distributed across the short, medium and long product pools. 
!! This basically means that wood use in ..., 1600, 1700, 1800, 1900 and 2000 
!! was identical.\n
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : prod_s, prod_m, prod_l 
!!
!! REFERENCES   : None
!!
!! FLOWCHART    : None
!!
!_ ================================================================================================================================

 SUBROUTINE product_add(npts, dt_days, harvest_pool_bound, harvest_pool_acc, &
            harvest_type, harvest_cut, harvest_area, &
            prod_s, prod_m, prod_l, prod_s_total, prod_m_total, &
            prod_l_total, flux_s, flux_m, flux_l, veget_max, flux_s_pft)

   IMPLICIT NONE


 !! 0. Variable and parameter declaration 
    
    !! 0.1 Input variables
   
    INTEGER, INTENT(in)                                  :: npts              !! Domain size - number of pixels (unitless)
    REAL(r_std), INTENT(in)                              :: dt_days           !! Time step of vegetation dynamics for stomate
                                                                              !! (days)
    REAL(r_std), DIMENSION(0:ndia_harvest+1), INTENT(in) :: harvest_pool_bound!! The boundaries of the diameter classes
                                                                              !! in the wood harvest pools @tex $(m)$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(in)              :: veget_max         !! Passed so we can check_mass_balance 
                                                                              !! Making it OPTIONAL would have been nicer.


 !! 0.2 Output variables

    REAL(r_std), DIMENSION(:,:,:,:), INTENT(out)         :: prod_s_total      !! Total products remaining in the short-lived  
                                                                              !! pool after the annual decomposition 
                                                                              !! @tex $(gC)$ @endtex 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(out)         :: prod_m_total      !! Total products remaining in the medium-lived
                                                                              !! pool after the annual decomposition 
                                                                              !! @tex $(gC)$ @endtex 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(out)         :: prod_l_total      !! Total products remaining in the long-lived
                                                                              !! pool after the annual release 
                                                                              !! @tex $(gC)$ @endtex 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(out)         :: flux_s_pft        !! Total flux from decomposition the same year of the lcc and harvest
                                                                              !! @tex $(gC m^{-1} dt^{-1})$ @endtex

    !! 0.3 Modified variables

    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)     :: prod_s            !! Short-lived product pool after the annual 
                                                                              !! release of each compartment (short + 1 : 
                                                                              !! input from year of land cover change) 
                                                                              !! @tex ($gC pixel^{-1}$) @endtex    
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)     :: prod_m            !! Medium-lived product pool after the annual 
                                                                              !! release of each compartment (medium + 1 : 
                                                                              !! input from year of land cover change) 
                                                                              !! @tex ($gC pixel^{-1}$) @endtex 
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)     :: prod_l            !! long-lived product pool after the annual 
                                                                              !! release of each compartment (long + 1 : 
                                                                              !! input from year of land cover change) 
                                                                              !! @tex ($gC pixel^{-1}$) @endtex
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)     :: harvest_pool_acc  !! The wood and biomass that have been
                                                                              !! havested by humans @tex $(gC) pixel-1$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)           :: harvest_type      !! Type of management that resulted
                                                                              !! in the harvest (unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)           :: harvest_cut       !! Type of cutting that was used for the harvest
                                                                              !! (unitless)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)         :: harvest_area      !! Harvested area (m^{2})
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)     :: flux_s            !! Annual release from the short-lived product 
                                                                              !! pool @tex $(gC) pixel^{-1} year^{-1}$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)     :: flux_m            !! Annual release from the medium-lived product  
                                                                              !! pool @tex $(gC) pixel^{-1} year^{-1}$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)     :: flux_l            !! Annual release from the long-lived product
                                                                              !! pool @tex $(gC) pixel^{-1} year^{-1}$ @endtex

  
    
    !! 0.4 Local variables
 
    INTEGER(i_std)                                       :: ipts, ivm, iage   !! Indices (unitless)
    INTEGER(i_std)                                       :: imbc, idia, iele  !! Indices (unitless)
    INTEGER(i_std)                                       :: ilan, ilct        !! Indices (unitless)
    REAL(r_std), DIMENSION(nelements,nlctypes)           :: harvest_acc_tmp   !! Accumulated harvest per land cover type (gC pixel-1) 
    REAL(r_std), DIMENSION(nelements,nlanduse)           :: exported_biomass  !! Total biomass exported from the PFT (gC)
    REAL(r_std), DIMENSION(nelements,nlanduse)           :: fuelwood_biomass  !! Harvest with small dimensions that will
                                                                              !! be used as fire wood (gC)
    REAL(r_std), DIMENSION(nelements,nlanduse)           :: timber_biomass    !! Harvest with larger dimensions that will
                                                                              !! be used for all kinds of uses (gC)
    REAL(r_std), DIMENSION(nelements)                    :: residual          !! Residual @tex $(gC m^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm,nmbcomp,nelements)   :: check_intern      !! Contains the components of the internal
                                                                              !! mass balance chech for this routine
                                                                              !! @tex $(gC m^{-2} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)           :: closure_intern    !! Check closure of internal mass balance
                                                                              !! @tex $(gC m^{-2} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)           :: pool_start        !! Start pool of this routine 
                                                                              !! @tex $(gC m^{-2} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)           :: pool_end          !! End pool of this routine 
                                                                              !! @tex $(gC m^{-2} dt^{-1})$ @endtex

!_ ================================================================================================================================    

    !! 1. initialization

    !! 1.1 Initialize local printlev
    ! Only first time one subroutine in the module is called
      ! Guillaume M. -- get_printlev calls getin_p, a COLLECTIVE bcast. Called here,
      ! inside the compute loop under a local firstcall, it hangs the run as soon as the
      ! ranks stop reaching it together. Use printlev, read once at startup by every rank.
    IF (firstcall_sapiens_product_add) THEN
      printlev_loc=printlev
      firstcall_sapiens_product_add=.FALSE.
    END IF

    IF (printlev_loc.GE.2) WRITE(numout,*) 'Entering product add'
    
    !! 1.2 Fluxes and pools
    prod_s(:,1,:,:,:) = zero
    prod_m(:,1,:,:,:) = zero
    prod_l(:,1,:,:,:) = zero
    flux_s_pft(:,:,:,:) = zero
    
    !! 1.4 Initialize check for mass balance closure 
    IF (err_act.GT.1) THEN

       ! The biomass harvest pool shouldn't be multiplied by veget_max
       ! its units are already in gC pixel-1. The total amount for the 
       ! pixel was stored in harvest_pool_acc. Account for the harvest and
       ! wood product pools. There are no longer PFTs
       check_intern(:,:,:,:) = zero
       pool_start(:,:,:) = zero
       DO iele = 1,nelements
          pool_start(:,1,iele) = &
               ( SUM(SUM(SUM(harvest_pool_acc(:,:,:,iele,:),4),3),2) + &
               SUM(SUM(SUM(prod_l(:,:,iele,:,:),2),2),2) + &
               SUM(SUM(SUM(prod_m(:,:,iele,:,:),2),2),2) + &
               SUM(SUM(SUM(prod_s(:,:,iele,:,:),2),2),2) ) / &
               (area(:) * contfrac(:))
       ENDDO
       
    ENDIF ! err_act.GT.1

    !! 1.5 Logical check of coeff_lcchange
    !  Within a PFT the coeff_lcchange should add up to one. Else the 
    !  carbon balance will not be closed! Do not include PFT 1 in this 
    !  check because its coeff_lcchange were set to undef (-9999).
    DO ivm = 2,nvm

       IF ( ABS(coeff_lcchange_s(ivm) + coeff_lcchange_m(ivm) + &
            coeff_lcchange_l(ivm)-un) .GT. min_stomate ) THEN

          WRITE(numout,*) 'Error: coeff_lcchange in sapiens_product_add.f90 ' //&
               'do not add up to 1 in PFT - cannot close the mass balance, ',ivm
          CALL ipslerr_p (3,'sapiens_product_add', &
               'product_add','coeff_lcchange do not add up to 1 ','')
         
       ENDIF

    ENDDO


    !! 2. Partition harvest in a short, medium and long-lived product pool
    
    DO ipts = 1,npts

       harvest_acc_tmp(:,:) = zero

       !! 2.1 Calculate the initial mass of the product pool (gC)
       DO ivm = 2,nvm

          ! Determine the land cover type
          IF (is_tree(ivm)) THEN
             ilct = iforest
          ELSEIF (.NOT. is_tree(ivm) .AND. natural(ivm)) THEN
             ilct = igrass
          ELSEIF (.NOT. is_tree(ivm) .AND. .NOT. natural(ivm)) THEN
             ilct = icrop
          ELSE
             CALL ipslerr_p (3,'sapiens_product_add.f90',&
                  'not clear to which land cover type the PFT belongs',&
                  'check the code', '')
          END IF
 
          DO iele = 1,nelements
             harvest_acc_tmp(iele,ilct) = harvest_acc_tmp(iele,ilct) + &
                  SUM(SUM(harvest_pool_acc(ipts,ivm,:,iele,:),2))
          END DO
          
          IF (ok_dimensional_product_use) THEN
       
             ! Wood harvest with small dimensions is being used as fire wood
             ! the remaining wood goes into the short, medium and long product
             ! pools. This implies that when the dimensions of the harvest 
             ! change, the product pools will change as well

             ! Product use is calculated separatly for the biomass export due
             ! to lcc and harvest. Be aware that for the 
             ! moment these values can be overwritten within a simulation.
             ! Although the sequence of the calls in stomate_lpj is believed
             ! to protect against this, it is not enforced by the code. For
             ! example, if we first thin and then harvest, wood of both
             ! actions will be stored in the harvest pool but it will be 
             ! recored as just harvest (thinning will be overwritten by harvest)
             fuelwood_biomass(:,:) = zero
             timber_biomass(:,:) = zero
             exported_biomass(:,:) = zero
             
             IF (fuelwood_diameter(ivm) .LT. min_stomate) THEN

                ! Code for grasslands and croplands.
                ! All biomass exported.  No distinction between
                ! fuelwood and timber.
                DO iele = 1,nelements
                   exported_biomass(iele,ilcc) = &
                        SUM(harvest_pool_acc(ipts,ivm,:,iele,ilcc))
                   exported_biomass(iele,iharvest) = &
                        SUM(harvest_pool_acc(ipts,ivm,:,iele,iharvest))
                END DO
                
                ! Attribute the harvested biomass to pools with a different 
                ! turnover time. The turnover times are PFT depend which 
                ! accounts for the different uses different PFT's are put to.
                ! Grasses and crops do not provide fuelwood. They do have
                ! a very short turnover time. Medium and long-lived are 
                ! nevertheless considered in case grasses and crops become
                ! to be used in longer-lived products (e.g. bio-plastics)
                ! It also makes the code more consistent. The values stored in 
                ! harvest_pool_acc are absolute numbers gC (for that pixel), hence, 
                ! the spatial extent of the LCC change, harvest or thinning has 
                ! already been accounted for.
                prod_s(ipts,1,:,:,ilct) = prod_s(ipts,1,:,:,ilct) + &
                     (coeff_lcchange_s(ivm) * exported_biomass(:,:))
                prod_m(ipts,1,:,:,ilct) = prod_m(ipts,1,:,:,ilct) + &
                     (coeff_lcchange_m(ivm) * exported_biomass(:,:))
                prod_l(ipts,1,:,:,ilct) = prod_l(ipts,1,:,:,ilct) + &
                     (coeff_lcchange_l(ivm) * exported_biomass(:,:))

                ! Calculate how much C and N is released back into the
                ! atmosphere the same year (needed for AR6 output)
                ! Note: when nshort=1, this makes completely sense. For nshort>1,
                ! this only accounts for the carbon release happening within the
                ! the year of harvest.
                DO iele = 1, nelements
                   flux_s_pft(ipts,ivm,iele,ilct) = &
                        coeff_lcchange_s(ivm) * SUM(exported_biomass(iele,:)) / nshort
                END DO
                
             ELSE

                ! Use of forest biomass.
                ! All biomass is assigned to timber. Later the fuelwood will
                ! be subtracted. This ensure mass balance closure.
                DO iele = 1,nelements
                   timber_biomass(iele,ilcc) = &
                        SUM(harvest_pool_acc(ipts,ivm,:,iele,ilcc))
                   timber_biomass(iele,iharvest) = &
                        SUM(harvest_pool_acc(ipts,ivm,:,iele,iharvest))
                END DO

                ! Use the dimensions of the stems to move the biomass into
                ! the right product pool
                DO idia = 1,ndia_harvest

                   IF (harvest_pool_bound(idia) .LE. fuelwood_diameter(ivm)) THEN

                      ! If the diameter of the harvest pool is below the treshold
                      ! the wood is used as fuelwood and will end up in the short
                      ! lived pool
                      DO iele = 1,nelements
                         fuelwood_biomass(iele,ilcc) = fuelwood_biomass(iele,ilcc) + & 
                              harvest_pool_acc(ipts,ivm,idia,iele,ilcc)
                         fuelwood_biomass(iele,iharvest) = &
                              fuelwood_biomass(iele,iharvest) + & 
                              harvest_pool_acc(ipts,ivm,idia,iele,iharvest)
                      END DO

                   ENDIF

                ENDDO ! harvest diameter

                ! Calculate the long lived pools as the residual of all
                ! wood mines the fuel wood.
                timber_biomass(:,:) = timber_biomass(:,:) - fuelwood_biomass(:,:)

                ! Attribute the harvested biomass to pools with a different 
                ! turnover time. The turnover times are PFT depend which 
                ! accounts for the different uses different PFT's are put to i.e. 
                ! wood from forests vs grass or crops. The values stored in 
                ! harvest_pool_acc are absolute numbers gC (for that pixel), hence, 
                ! the spatial extent of the LCC change, harvest or thinning has 
                ! already been accounted for.
                prod_s(ipts,1,:,:,ilct) = prod_s(ipts,1,:,:,ilct) + fuelwood_biomass(:,:)

                ! Calculate how much C and N is released back into the
                ! atmosphere the same year (needed for AR6 output)
                DO iele = 1, nelements
                   flux_s_pft(ipts,ivm,iele,ilct) = (fuelwood_biomass(iele,ilcc) + &
                        fuelwood_biomass(iele,iharvest)) / nshort
                END DO

                ! Allocated timber to medium and long-lived product pool
                IF  (coeff_lcchange_m(ivm) + coeff_lcchange_l(ivm) .GT. &
                     min_stomate) THEN

                   ! If the condition is satisfied, calculate the fractions
                   ! if it is not satisfied there is no medium and long lived
                   ! pool for the products. Nothing should then be done. Avoid
                   ! divide by zero
                   prod_m(ipts,1,:,:,ilct) = prod_m(ipts,1,:,:,ilct) + &
                        (coeff_lcchange_m(ivm)/ &
                        (coeff_lcchange_m(ivm)+coeff_lcchange_l(ivm)) * &
                        timber_biomass(:,:))
                   prod_l(ipts,1,:,:,ilct) = prod_l(ipts,1,:,:,ilct) + &
                        (coeff_lcchange_l(ivm)/ &
                        (coeff_lcchange_m(ivm)+coeff_lcchange_l(ivm)) * &
                        timber_biomass(:,:))

                ENDIF ! Check coeff_lcchange

             ENDIF ! grasses/crops vs forests

          ELSE

             ! The scheme does not account for the harvest dimensions. All
             ! wood is distributed across the short, medium and long product
             ! pools. This basically means that wood use in 1600, 1700, 1800,
             ! 1900 and 2000 was identical.

             ! Product use is calculated separatly for the biomass export due
             ! to lcc and harvest. Be aware that for the 
             ! moment these values can be overwritten within a simulation.
             ! Although the sequence of the calls in stomate_lpj is believed
             ! to protect against this, it is not enforced by the code. For
             ! example, if we first thin and then harvest, wood of both
             ! actions will be stored in the harvest pool but it will be 
             ! recored as just harvest (thinning will be overwritten by harvest)
             exported_biomass(:,:) = zero

             ! All biomass exported. No distinction
             ! between fuelwood and timber.
             DO iele = 1,nelements
                exported_biomass(iele,ilcc) = &
                     SUM(harvest_pool_acc(ipts,ivm,:,iele,ilcc))
                exported_biomass(iele,iharvest) = &
                     SUM(harvest_pool_acc(ipts,ivm,:,iele,iharvest))
             END DO

             !! 2.2 Process the exported biomass
             ! Attribute the harvested biomass to pools with a different 
             ! turnover time. The turnover times are PFT depend which 
             ! accounts for the different uses different PFT's are put to i.e. 
             ! wood from forests vs grass or crops. The values stored in 
             ! harvest_pool_acc are absolute numbers gC (for that pixel), hence, 
             ! the spatial extent of the LCC change, harvest or thinning has 
             ! already been accounted for.
             prod_s(ipts,1,:,:,ilct) = prod_s(ipts,1,:,:,ilct) + &
                  (coeff_lcchange_s(ivm) * exported_biomass(:,:))
             prod_m(ipts,1,:,:,ilct) = prod_m(ipts,1,:,:,ilct) + &
                  (coeff_lcchange_m(ivm) * exported_biomass(:,:))
             prod_l(ipts,1,:,:,ilct) = prod_l(ipts,1,:,:,ilct) + &
                  (coeff_lcchange_l(ivm) * exported_biomass(:,:))

             ! Calculate how much C and N is released back into the
             ! atmosphere the same year (needed for AR6 output)
             ! Note: when nshort=1, this makes completely sense. For nshort>1,
             ! this only accounts for the carbon release happening within the
             ! the year of harvest.
             DO iele = 1, nelements
                flux_s_pft(ipts,ivm,iele,ilct) = &
                     coeff_lcchange_s(ivm) * SUM(exported_biomass(iele,:)) / nshort
             END DO

          END IF ! ok_dimensional_product_use   
             
       ENDDO ! loop over PFTs 

       ! Error checking and ensuring mass balance closure
       DO ilct = 1,nlctypes

          ! Check whether all the harvest was used as is to be expected
          ! The precision we are aiming for is 10-8 per m-2. The harvest pool
          ! is expressed for the whole pixel
          residual(:) = harvest_acc_tmp(:,ilct) - &
               SUM(prod_s(ipts,1,:,:,ilct),2) - &
               SUM(prod_m(ipts,1,:,:,ilct),2) - &
               SUM(prod_l(ipts,1,:,:,ilct),2)
    
          DO iele = 1,nelements

             IF (ABS(residual(iele)/(area(ipts)*contfrac(ipts))) .GT. min_stomate) THEN
                
                WRITE (numout,*) 'Error: some harvest in '//&
                     'sapiens_product_add has not been attributed'
                WRITE (numout,*) 'Residual in pixel, iele, ', &
                     iele, ipts, residual
                CALL ipslerr_p (3,'sapiens_product_add', &
                     'product_add','some harvest has not been attributed','')

             ELSE

                ! We could still have residuals within the precision at the 
                ! m2 level but outside the precision at the pixel level. 
                ! Move the residual to the short lived pool.
                IF ( prod_s(ipts,1,iele,iharvest,ilct).GT.zero) THEN
                   
                   ! If there was harvest put the residual in the product pool of
                   ! wood harvest
                   prod_s(ipts,1,iele,iharvest,ilct) = &
                        prod_s(ipts,1,iele,iharvest,ilct) + &
                        residual(iele)
                   
                ELSE
                
                   ! There was no harvest but there is some residual so the residual
                   ! should come from land cover change
                   prod_s(ipts,1,iele,ilcc,ilct) = prod_s(ipts,1,iele,ilcc,ilct) + &
                        residual(iele)

                ENDIF
             
             ENDIF
             
          END DO ! iele

       END DO ! land cover types

       ! Calculate the annual decomposition - this is used in product_decomp
       flux_s(ipts,1,:,:,:) = (un/nshort) * prod_s(ipts,1,:,:,:)
       flux_m(ipts,1,:,:,:) = (un/nmedium) * prod_m(ipts,1,:,:,:)
       flux_l(ipts,1,:,:,:) = (un/nlong) * prod_l(ipts,1,:,:,:)
       
       ! All the carbon in the harvest pool was allocated to the product 
       ! pools the harvest pool can be set to zero so we can us it with 
       ! confidence in the mass balance calculation
       harvest_pool_acc(ipts,:,:,:,:) = zero
       harvest_type(ipts,:) = zero
       harvest_cut(ipts,:) = zero
       harvest_area(ipts,:,:) = zero

    END DO ! #pixels - spatial domain

 !! 3. Check numerical consistency of this routine

    IF (err_act.GT.1) THEN

       ! 3.1 Mass balance closure 
       ! 3.1.1 Calculate final biomass
       ! The biomass harvest pool shouldn't be multiplied by veget_max
       ! its in gC pixel-1 so divide by area to get a value in
       ! gC m-2. Account for the harvest and wood product pools. The harvest 
       ! pool should be empty, it was added for completeness.
       pool_end = zero
       DO iele=1,nelements
          pool_end(:,1,iele) = &
               ( SUM(SUM(SUM(harvest_pool_acc(:,:,:,iele,:),4),3),2) + &
               SUM(SUM(SUM(prod_l(:,:,iele,:,:),2),2),2) + &
               SUM(SUM(SUM(prod_m(:,:,iele,:,:),2),2),2) + &
               SUM(SUM(SUM(prod_s(:,:,iele,:,:),2),2),2) ) / &
               (area(:) * contfrac(:))
       ENDDO

       !! 3.1.2 Calculate mass balance
       ! Common processes
       DO iele=1,nelements
          check_intern(:,1,ipoolchange,iele) = &
               -un * (pool_end(:,1,iele) - &
               pool_start(:,1,iele))
       ENDDO

       closure_intern = zero
       DO imbc = 1,nmbcomp
          DO iele=1,nelements
             closure_intern(:,1,iele) = closure_intern(:,1,iele) + &
                  check_intern(:,1,imbc,iele)
          ENDDO
       ENDDO
       
       ! 3.1.3 Check mass balance closure
       CALL check_mass_balance("product_add", closure_intern, npts, &
            pool_end, pool_start, veget_max, 'products')
      
    ENDIF ! err_act.GT.1

 !! 4. Convert the pools into fluxes and aggregate

    !  Convert pools into fluxes
    !  The pools are given in the absolute amount of carbon for the 
    !  whole pixel and the harvest was accumulated over the year. 
    !  The unit is thus gc pixel-1 year-1. These fluxes should be 
    !  expressed as gC m-2 dt-1. Therefore, divide by the area 
    !  of the pixel and the number of time steps within a year
    DO iele = 1,nelements
       
       DO ilan = 1,nlanduse

          DO ilct = 1,nlctypes
          
             ! Convert prod_x_total in gC m-2. prod_x_total is not used in any mass
             ! balance check so there is no need to use a residual term to achieve
             ! higher precision.
             prod_s_total(:,iele,ilan,ilct) = SUM(prod_s(:,:,iele,ilan,ilct),2) 
             prod_m_total(:,iele,ilan,ilct) = SUM(prod_m(:,:,iele,ilan,ilct),2) 
             prod_l_total(:,iele,ilan,ilct) = SUM(prod_l(:,:,iele,ilan,ilct),2) 
             
          END DO
          
       END DO

       ! AR6 output. Convert a pool into a flux.
       flux_s_pft(:,:,iele,:) =  flux_s_pft(:,:,iele,:) / one_year * dt_days 

    END DO

    IF (printlev.GE.4) WRITE(numout,*) 'Leaving product add'


  END SUBROUTINE product_add

  
!! ================================================================================================================================
!! SUBROUTINE   : product_decomp
!!
!>\BRIEF        Calculate the CO2 release through product decomposition. 
!!
!! DESCRIPTION  : 
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : prod_s, prod_m, prod_l, flux_s, flux_m, flux_l,
!!   flux_prod_s flux_prod_m and flux_prod_l 
!!
!! REFERENCES   : None
!!
!! FLOWCHART    : None
!!
!_ ================================================================================================================================

 SUBROUTINE product_decomp(npts, dt_days, &
            flux_prod_s, flux_prod_m, flux_prod_l, prod_s, &
            prod_m, prod_l, prod_s_total, prod_m_total, &
            prod_l_total, flux_prod_total, flux_s, flux_m, &
            flux_l, veget_max)

   IMPLICIT NONE


 !! 0. Variable and parameter declaration 
    
    !! 0.1 Input variables
   
    INTEGER, INTENT(in)                                  :: npts              !! Domain size - number of pixels (unitless)
    REAL(r_std), INTENT(in)                              :: dt_days           !! Time step of vegetation dynamics for stomate
                                                                              !! (days)
    REAL(r_std), DIMENSION(:,:), INTENT(in)              :: veget_max         !! Passed so we can check_mass_balance 
                                                                              !! Making it OPTIONAL would have been nicer.


 !! 0.2 Output variables

    REAL(r_std), DIMENSION(:,:,:,:), INTENT(out)         :: flux_prod_s       !! C-released during first years (short term) 
                                                                              !! following land cover change 
                                                                              !! @tex ($gC m^{-1} dt^{-1}$) @endtex 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(out)         :: flux_prod_m       !! Total annual release from decomposition of 
                                                                              !! the medium-lived product pool 
                                                                              !! @tex ($gC m^{-1} dt^{-1}$) @endtex
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(out)         :: flux_prod_l       !! Total annual release from decomposition of 
                                                                              !! the long-lived product pool 
                                                                              !! @tex ($gC m^{-1} dt^{-1}$) @endtex
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(out)         :: prod_s_total      !! Total products remaining in the short-lived  
                                                                              !! pool after the annual decomposition 
                                                                              !! @tex $(gC)$ @endtex 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(out)         :: prod_m_total      !! Total products remaining in the medium-lived
                                                                              !! pool after the annual decomposition 
                                                                              !! @tex $(gC)$ @endtex 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(out)         :: prod_l_total      !! Total products remaining in the long-lived
                                                                              !! pool after the annual release 
                                                                              !! @tex $(gC)$ @endtex 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(out)         :: flux_prod_total   !! Total flux from decomposition of the short,
                                                                              !! medium and long lived product pools 
                                                                              !! @tex $(gC m^{-1} dt^{-1})$ @endtex

    
    !! 0.3 Modified variables

    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)     :: prod_s            !! Short-lived product pool after the annual 
                                                                              !! release of each compartment (short + 1 : 
                                                                              !! input from year of land cover change) 
                                                                              !! @tex ($gC pixel^{-1}$) @endtex    
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)     :: prod_m            !! Medium-lived product pool after the annual 
                                                                              !! release of each compartment (medium + 1 : 
                                                                              !! input from year of land cover change) 
                                                                              !! @tex ($gC pixel^{-1}$) @endtex 
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)     :: prod_l            !! long-lived product pool after the annual 
                                                                              !! release of each compartment (long + 1 : 
                                                                              !! input from year of land cover change) 
                                                                              !! @tex ($gC pixel^{-1}$) @endtex
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)     :: flux_s            !! Annual release from the short-lived product 
                                                                              !! pool @tex $(gC) pixel^{-1} year^{-1}$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)     :: flux_m            !! Annual release from the medium-lived product  
                                                                              !! pool @tex $(gC) pixel^{-1} year^{-1}$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)     :: flux_l            !! Annual release from the long-lived product
                                                                              !! pool @tex $(gC) pixel^{-1} year^{-1}$ @endtex


    !! 0.4 Local variables
 
    INTEGER(i_std)                                       :: ipts, ivm, iage   !! Indices (unitless)
    INTEGER(i_std)                                       :: imbc, idia, iele  !! Indices (unitless)
    INTEGER(i_std)                                       :: ilan, ilct        !! Indices (unitless)
    REAL(r_std), DIMENSION(nelements,nlctypes)           :: harvest_acc_tmp   !! Accumulated harvest per land cover type (gC pixel-1) 
    REAL(r_std), DIMENSION(nelements,nlanduse)           :: exported_biomass  !! Total biomass exported from the PFT (gC)
    REAL(r_std), DIMENSION(nelements,nlanduse)           :: fuelwood_biomass  !! Harvest with small dimensions that will
                                                                              !! be used as fire wood (gC)
    REAL(r_std), DIMENSION(nelements,nlanduse)           :: timber_biomass    !! Harvest with larger dimensions that will
                                                                              !! be used for all kinds of uses (gC)
    REAL(r_std), DIMENSION(nelements)                    :: residual          !! Residual @tex $(gC m^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm,nmbcomp,nelements)   :: check_intern      !! Contains the components of the internal
                                                                              !! mass balance chech for this routine
                                                                              !! @tex $(gC m^{-2} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)           :: closure_intern    !! Check closure of internal mass balance
                                                                              !! @tex $(gC m^{-2} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)           :: pool_start        !! Start pool of this routine 
                                                                              !! @tex $(gC m^{-2} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)           :: pool_end          !! End pool of this routine 
                                                                              !! @tex $(gC m^{-2} dt^{-1})$ @endtex
    REAL(r_std)                                          :: daily_share       !! Daily share in decomposition
!_ ================================================================================================================================    


!! 1. initialization

    !! 1.1 Initialize local printlev only first time one subroutine in the module is called
      ! Guillaume M. -- get_printlev calls getin_p, a COLLECTIVE bcast. Called here,
      ! inside the compute loop under a local firstcall, it hangs the run as soon as the
      ! ranks stop reaching it together. Use printlev, read once at startup by every rank.
    IF (firstcall_sapiens_product_decomp) THEN
      printlev_loc=printlev
      firstcall_sapiens_product_decomp=.FALSE.
    END IF

    IF (printlev_loc.GE.2) WRITE(numout,*) 'Entering product decomposition'

    !! 1.2 Fluxes and pools
    prod_s_total(:,:,:,:) = zero
    prod_m_total(:,:,:,:) = zero
    prod_l_total(:,:,:,:) = zero   
    flux_prod_s(:,:,:,:) = zero
    flux_prod_m(:,:,:,:) = zero
    flux_prod_l(:,:,:,:) = zero
    flux_prod_total(:,:,:,:) = zero

    ! Take 366 days for leap years. When no leap years are used the last
    ! day decomposes twice as much as other days. This may helps to prevent
    ! mass balance problems. This is a preventive idea, not sure whether
    ! it is needed. The cost the small: ~2/365th of emissions the last day
    ! of the year instead of 1/365th.
    daily_share = 1./REAL(one_year+1)
    
    !! 1.3 Initialize check for mass balance closure 
    IF (err_act.GT.1) THEN

       ! The biomass harvest pool shouldn't be multiplied by veget_max
       ! its units are already in gC pixel-1. The total amount for the 
       ! pixel was stored in harvest_pool_acc. Account for the harvest and
       ! wood product pools. There are no longer PFTs
       check_intern(:,:,:,:) = zero
       pool_start(:,:,:) = zero
       DO iele = 1,nelements
          pool_start(:,1,iele) = &
               ( SUM(SUM(SUM(prod_l(:,:,iele,:,:),2),2),2) + &
               SUM(SUM(SUM(prod_m(:,:,iele,:,:),2),2),2) + &
               SUM(SUM(SUM(prod_s(:,:,iele,:,:),2),2),2) ) / &
               (area(:) * contfrac(:))
       ENDDO
       
    ENDIF ! err_act.GT.1
    
    !! 2. Update the product pools and decomposition fluxes
    IF (.NOT.ts_annual_proc) THEN

       ! Repeat this for all days except the last day of the year
       ! Decompose 1/365th of the annual flux. Update the product
       ! pools and decomposition fluxes
       DO iele = 1,nelements
          DO ilan = 1,nlanduse
             DO ilct = 1,nlctypes
                flux_prod_l(:,iele,ilan,ilct) = flux_prod_l(:,iele,ilan,ilct) + &
                     daily_share*SUM(flux_l(:,:,iele,ilan,ilct),2)
                prod_l(:,:,iele,ilan,ilct) = prod_l(:,:,iele,ilan,ilct) - &
                     daily_share*flux_l(:,:,iele,ilan,ilct)
                flux_prod_m(:,iele,ilan,ilct) = flux_prod_m(:,iele,ilan,ilct) + &
                     daily_share*SUM(flux_m(:,:,iele,ilan,ilct),2)
                prod_m(:,:,iele,ilan,ilct) = prod_m(:,:,iele,ilan,ilct) - &
                     daily_share*flux_m(:,:,iele,ilan,ilct)
                flux_prod_s(:,iele,ilan,ilct) = flux_prod_s(:,iele,ilan,ilct) + &
                     daily_share*SUM(flux_s(:,:,iele,ilan,ilct),2)
                prod_s(:,:,iele,ilan,ilct) = prod_s(:,:,iele,ilan,ilct) - &
                     daily_share*flux_s(:,:,iele,ilan,ilct)
             END DO ! nlctypes
          END DO ! nlanduse
       END DO ! nelements
       
    ELSE
       
       ! Update the pools and fluxes the last day of the year. Move
       ! the cohorts to the next age class, so that the harvest of
       ! the previous 365 days can be added to the product pools.
       DO iele = 1,nelements
          DO ilan = 1,nlanduse
             DO ilct = 1,nlctypes

                ! Empty the last age class.
                ! Last age class should become zero. Treat separatly to close the
                ! mass balance.We take 364 times 1/365th of the flux. For the last
                ! time step we could take 1/365 of flux_x again but this could
                ! cause a mass balance problem. Instead we take prod_x as this is
                ! the residual of 364*1/365*flux_x. This should close the mass
                ! balance problem.
                flux_prod_l(:,iele,ilan,ilct) = flux_prod_l(:,iele,ilan,ilct) + &
                     prod_l(:,nlong,iele,ilan,ilct)
                flux_l(:,nlong,iele,ilan,ilct) = zero
                prod_l(:,nlong,iele,ilan,ilct) = zero
                flux_prod_m(:,iele,ilan,ilct) = flux_prod_m(:,iele,ilan,ilct) + &
                     prod_m(:,nmedium,iele,ilan,ilct)
                flux_m(:,nmedium,iele,ilan,ilct) = zero
                prod_m(:,nmedium,iele,ilan,ilct) = zero
                flux_prod_s(:,iele,ilan,ilct) = flux_prod_s(:,iele,ilan,ilct) + &
                     prod_s(:,nshort,iele,ilan,ilct)
                flux_s(:,nshort,iele,ilan,ilct) = zero
                prod_s(:,nshort,iele,ilan,ilct) = zero

                END DO ! nlctypes
          END DO ! nlanduse
       END DO ! nelements

       ! Update the other age classes/cohorts
       DO iele = 1,nelements
          DO ilan = 1,nlanduse
             DO ilct = 1,nlctypes

                ! Long-lived pools
                ! Update all age classes at once
                flux_prod_l(:,iele,ilan,ilct) = flux_prod_l(:,iele,ilan,ilct) + &
                     daily_share*SUM(flux_l(:,:,iele,ilan,ilct),2)
                prod_l(:,:,iele,ilan,ilct) = prod_l(:,:,iele,ilan,ilct) - &
                     daily_share*flux_l(:,:,iele,ilan,ilct)
                
                ! Move all age classes one up. The last age class was
                ! emptied in the previously
                DO iage = nlong+1,2,-1       
                   prod_l(:,iage,iele,ilan,ilct) = prod_l(:,iage-1,iele,ilan,ilct) 
                   flux_l(:,iage,iele,ilan,ilct) = flux_l(:,iage-1,iele,ilan,ilct)
                ENDDO

                ! Set zeroth age class to zero so there
                ! is space to store the harvest of the past year
                prod_l(:,1,iele,ilan,ilct) = zero

                ! Medium-lived pools
                ! Update all age classes at once
                flux_prod_m(:,iele,ilan,ilct) = flux_prod_m(:,iele,ilan,ilct) + &
                     daily_share*SUM(flux_m(:,:,iele,ilan,ilct),2)
                prod_m(:,:,iele,ilan,ilct) = prod_m(:,:,iele,ilan,ilct) - &
                     daily_share*flux_m(:,:,iele,ilan,ilct)

                ! Move all age classes one up. The last age class was
                ! emptied in the previously
                DO iage = nmedium+1,2,-1
                   prod_m(:,iage,iele,ilan,ilct) = prod_m(:,iage-1,iele,ilan,ilct)
                   flux_m(:,iage,iele,ilan,ilct) = flux_m(:,iage-1,iele,ilan,ilct)
                ENDDO

                ! Set zeroth age class to zero so there
                ! is space to store the harvest of the past year
                prod_m(:,1,iele,ilan,ilct) = zero

                ! Short-lived pools
                IF (nshort .GE. 2) THEN

                   ! Update all age classes at once
                   flux_prod_s(:,iele,ilan,ilct) = flux_prod_s(:,iele,ilan,ilct) + &
                        daily_share*SUM(flux_s(:,:,iele,ilan,ilct),2)
                   prod_s(:,:,iele,ilan,ilct) = prod_s(:,:,iele,ilan,ilct) - &
                        daily_share*flux_s(:,:,iele,ilan,ilct)

                   ! Move all age classes one up. The last age class was
                   ! emptied in the previously
                   DO iage = nshort+1,2,-1
                      prod_s(:,iage,iele,ilan,ilct) = prod_s(:,iage-1,iele,ilan,ilct)
                      flux_s(:,iage,iele,ilan,ilct) = flux_s(:,iage-1,iele,ilan,ilct)
                   ENDDO

                ELSE

                   ! With only one age class, age class 1 and zero
                   ! should be emptied to make sure the mass balance
                   ! is preserved.
                   flux_prod_s(:,iele,ilan,ilct) = flux_prod_s(:,iele,ilan,ilct) + &
                        prod_s(:,1,iele,ilan,ilct)
                   
                END IF

                ! Set zeroth age class to zero so there
                ! is space to store the harvest of the past year
                prod_s(:,1,iele,ilan,ilct) = zero

             END DO ! nlctypes
          END DO ! nlanduse
       END DO ! nelements
       
    END IF ! ts_annual_proc
    

 !! 3. Check numerical consistency of this routine

    IF (err_act.GT.1) THEN

       ! 3.1 Mass balance closure 
       ! 3.1.1 Calculate final biomass
       ! The biomass harvest pool shouldn't be multiplied by veget_max
       ! its in gC pixel-1 so divide by area to get a value in
       ! gC m-2. Account for the harvest and wood product pools. The harvest 
       ! pool should be empty, it was added for completeness.
       pool_end = zero
       DO iele=1,nelements
          pool_end(:,1,iele) = &
               ( SUM(SUM(SUM(prod_l(:,:,iele,:,:),2),2),2) + &
               SUM(SUM(SUM(prod_m(:,:,iele,:,:),2),2),2) + &
               SUM(SUM(SUM(prod_s(:,:,iele,:,:),2),2),2) ) / &
               (area(:) * contfrac(:))
       ENDDO

       !! 3.1.2 Calculate mass balance
       ! Common processes
       DO iele=1,nelements
          check_intern(:,1,iland2atm,iele) = &
               -un * ( SUM(SUM(flux_prod_s(:,iele,:,:) + &
               flux_prod_m(:,iele,:,:) + &
               flux_prod_l(:,iele,:,:),2),2) ) / &
               (area(:) * contfrac(:))
          check_intern(:,1,ipoolchange,iele) = &
               -un * (pool_end(:,1,iele) - &
               pool_start(:,1,iele))
       ENDDO

       closure_intern = zero
       DO imbc = 1,nmbcomp
          DO iele=1,nelements
             closure_intern(:,1,iele) = closure_intern(:,1,iele) + &
                  check_intern(:,1,imbc,iele)
          ENDDO
       ENDDO
       
       ! 3.1.3 Check mass balance closure
       CALL check_mass_balance("product_decomp", closure_intern, npts, &
            pool_end, pool_start, veget_max, 'products')
      
    ENDIF ! err_act.GT.1

    !! 4. Convert the pools into fluxes and aggregate
 
    !  Convert pools into fluxes
    !  The pools are given in the absolute amount of carbon for the 
    !  whole pixel and the harvest was accumulated over the year. 
    !  The unit is thus gc pixel-1 year-1. These fluxes should be 
    !  expressed as gC m-2 dt-1. Therefore, divide by the area 
    !  of the pixel and the number of time steps within a year
    DO iele = 1,nelements
       DO ilan = 1,nlanduse
          DO ilct = 1,nlctypes
          
             ! Convert flux_prod_total and flux_prod_x to gC m-2 dt-1. Calculate
             ! flux_prod_l as the residual for better mass conservation.
             flux_prod_total(:,iele,ilan,ilct) = (flux_prod_s(:,iele,ilan,ilct) + &
                  flux_prod_m(:,iele,ilan,ilct) + flux_prod_l(:,iele,ilan,ilct)) / &
                  (area(:) * contfrac(:))
             flux_prod_s(:,iele,ilan,ilct) = flux_prod_s(:,iele,ilan,ilct) / &
                  (area(:) * contfrac(:)) 
             flux_prod_m(:,iele,ilan,ilct) = flux_prod_m(:,iele,ilan,ilct) / &
                  (area(:) * contfrac(:)) 
             flux_prod_l(:,iele,ilan,ilct) = flux_prod_total(:,iele,ilan,ilct) - &
                  flux_prod_s(:,iele,ilan,ilct) - flux_prod_m(:,iele,ilan,ilct)

             ! Convert prod_x_total in gC m-2. prod_x_total is not used in any mass
             ! balance check so there is no need to use a residual term to achieve
             ! higher precision.
             prod_s_total(:,iele,ilan,ilct) = SUM(prod_s(:,:,iele,ilan,ilct),2) 
             prod_m_total(:,iele,ilan,ilct) = SUM(prod_m(:,:,iele,ilan,ilct),2) 
             prod_l_total(:,iele,ilan,ilct) = SUM(prod_l(:,:,iele,ilan,ilct),2) 
             
          END DO ! iclt
       END DO ! ilan
    END DO ! iele

    IF (printlev.GE.4) WRITE(numout,*) 'Leaving product decomposition'


  END SUBROUTINE product_decomp

  
!! ================================================================================================================================
!! SUBROUTINE   : product_init
!!
!>\BRIEF        Initialize the product use variables when product use is not
!!              calculated 
!!
!! DESCRIPTION  : Initialize the product use variables when product use is not
!! calculated
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : prod_s, prod_m, prod_l, flux_s, flux_m, flux_l,
!!   flux_prod_s flux_prod_m and flux_prod_l 
!!
!! REFERENCES   : None
!!
!! FLOWCHART    : None
!!
!_ ================================================================================================================================

  SUBROUTINE product_init(flux_prod_s, flux_prod_m, flux_prod_l, &
               prod_s_total, prod_m_total, prod_l_total, flux_s, &
               flux_m, flux_l, flux_prod_total, flux_s_pft, &
               prod_s, prod_m, prod_l)

   IMPLICIT NONE


 !! 0. Variable and parameter declaration 
    
    !! 0.1 Input variables

    !! 0.2 Output variables

    REAL(r_std), DIMENSION(:,:,:,:), INTENT(out)         :: flux_prod_s       !! C-released during first years (short term) 
                                                                              !! following land cover change 
                                                                              !! @tex ($gC m^{-1} dt^{-1}$) @endtex 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(out)         :: flux_prod_m       !! Total annual release from decomposition of 
                                                                              !! the medium-lived product pool 
                                                                              !! @tex ($gC m^{-1} dt^{-1}$) @endtex
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(out)         :: flux_prod_l       !! Total annual release from decomposition of 
                                                                              !! the long-lived product pool 
                                                                              !! @tex ($gC m^{-1} dt^{-1}$) @endtex
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(out)         :: prod_s_total      !! Total products remaining in the short-lived  
                                                                              !! pool after the annual decomposition 
                                                                              !! @tex $(gC)$ @endtex 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(out)         :: prod_m_total      !! Total products remaining in the medium-lived
                                                                              !! pool after the annual decomposition 
                                                                              !! @tex $(gC)$ @endtex 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(out)         :: prod_l_total      !! Total products remaining in the long-lived
                                                                              !! pool after the annual release 
                                                                              !! @tex $(gC)$ @endtex 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(out)         :: flux_prod_total   !! Total flux from decomposition of the short,
                                                                              !! medium and long lived product pools 
                                                                              !! @tex $(gC m^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(out)         :: flux_s_pft        !! Total flux from decomposition the same year of the lcc and harvest
                                                                              !! @tex $(gC m^{-1} dt^{-1})$ @endtex

    !! 0.3 Modified variables

    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)     :: flux_s            !! Annual release from the short-lived product 
                                                                              !! pool @tex $(gC) pixel^{-1} year^{-1}$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)     :: flux_m            !! Annual release from the medium-lived product  
                                                                              !! pool @tex $(gC) pixel^{-1} year^{-1}$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)     :: flux_l            !! Annual release from the long-lived product
                                                                              !! pool @tex $(gC) pixel^{-1} year^{-1}$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)     :: prod_s            !! Short-lived product pool after the annual 
                                                                              !! release of each compartment (short + 1 : 
                                                                              !! input from year of land cover change) 
                                                                              !! @tex ($gC$) @endtex    
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)     :: prod_m            !! Medium-lived product pool after the annual
                                                                              !! release of each compartment (medium + 1 : 
                                                                              !! input from year of land cover change) 
                                                                              !! @tex ($gC$) @endtex 
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)     :: prod_l            !! long-lived product pool after the annual 
                                                                              !! release of each compartment (long + 1 : 
                                                                              !! input from year of land cover change) 
                                                                              !! @tex ($gC$) @endtex

    !! 0.4 Local variables

!_ ================================================================================================================================    

    IF (printlev.GE.2) WRITE(numout,*) 'Entering product init'
 
    !! 1. initialization
    prod_s(:,1,:,:,:) = zero
    prod_m(:,1,:,:,:) = zero
    prod_l(:,1,:,:,:) = zero
    flux_prod_s(:,:,:,:) = zero
    flux_prod_m(:,:,:,:) = zero
    flux_prod_l(:,:,:,:) = zero
    prod_s_total(:,:,:,:) = SUM(prod_s(:,:,:,:,:),2)
    prod_m_total(:,:,:,:) = SUM(prod_m(:,:,:,:,:),2)
    prod_l_total(:,:,:,:) = SUM(prod_l(:,:,:,:,:),2)
    flux_s_pft(:,:,:,:) = zero
    flux_prod_total(:,:,:,:) = zero

    IF (printlev.GE.4) WRITE(numout,*) 'Leaving product init'

  END SUBROUTINE product_init

END MODULE sapiens_product_use



