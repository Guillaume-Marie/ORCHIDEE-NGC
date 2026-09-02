! =================================================================================================================================
! MODULE 	: stomate_vmax
!
! CONTACT	: orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      	: IPSL (2006). This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        calculates the leaf efficiency.
!!	
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): None
!!
!! SVN		:
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_stomate/stomate_vmax.f90 $ 
!! $Date: 2025-08-14 17:10:02 +0200 (jeu. 14 août 2025) $
!! $Revision: 9070 $
!! \n
!_ =================================================================================================================================

MODULE stomate_vmax

  ! modules used:

  USE ioipsl_para
  USE xios_orchidee
  USE stomate_data
  USE constantes
  USE pft_parameters

  IMPLICIT NONE

  ! private & public routines

  PRIVATE
  PUBLIC vmax, vmax_clear

  ! first call
  LOGICAL, SAVE                                              :: firstcall_vmax = .TRUE.
!$OMP THREADPRIVATE(firstcall_vmax)

CONTAINS

!! ================================================================================================================================
!! SUBROUTINE	: vmax_clear
!!
!>\BRIEF	  Flag setting 
!!
!!\n DESCRIPTION: This subroutine sets flags ::firstcall_vmax, to .TRUE., and therefore activates   
!!		  section 1.1 of the ::vmax subroutine which writes messages to the output. \n
!!		  This subroutine is called at the end of the subroutine ::stomate_clear, in the 
!!		  module ::stomate.
!!
!! RECENT CHANGE(S):None
!!
!! MAIN OUTPUT VARIABLE(S): ::firstcall_vmax
!!
!! REFERENCE(S)  : None 
!!
!! FLOWCHART     : None
!! \n		  
!_ =================================================================================================================================

  SUBROUTINE vmax_clear
    firstcall_vmax=.TRUE.
  END SUBROUTINE vmax_clear



!! ================================================================================================================================
!! SUBROUTINE    : vmax
!!
!>\BRIEF         This subroutine computes vcmax photosynthesis parameters 
!! given optimal vcmax parameter values and a leaf age-related efficiency.
!!
!! DESCRIPTION (functional, design, flags): 
!! Leaf age classes are introduced to take into account the fact that photosynthetic activity depends on leaf age
!! (Ishida et al., 1999). There are \f$nleafages\f$ classes (constant defined in stomate_constants.f90).
!! This subroutine first calculates the new age of each leaf age-class based on fraction of leaf 
!! that goes from one to another class.                                              
!! Then calculation of the new fraction of leaf in each class is performed.      
!! Last, leaf efficiency is calculated for each PFT and for each leaf age class.
!! vcmax is defined as vcmax25 and vjmax_opt weighted by a mean leaf
!! efficiency. vcmax25 is PFT-dependent constants defined in constants_mtc.f90.
!!
!! This routine is called once at the beginning by stomate_var_init and then at each stomate time step by stomateLpj.
!!
!! RECENT CHANGE(S): 
!!       1) In r5639, calculation of nitrogen use efficiency was added, which will constrain vcmax.
!!          The nue efficiency takes into account varying efficiencies of different leaf age classes.
!!
!! MAIN OUTPUT VARIABLE(S): vcmax, nue
!!
!! REFERENCE(S)	: 
!! - Ishida, A., A. Uemura, N. Koike, Y. Matsumoto, and A. Lai Hoe (1999),
!! Interactive effects of leaf age and self-shading on leaf structure, photosynthetic
!! capacity and chlorophyll fluorescence in the rain forest tree,
!! dryobalanops aromatica, Tree Physiol., 19, 741-747
!!
!! FLOWCHART    : None
!!
!! REVISION(S)	: None
!! \n
!_ ================================================================================================================================

  SUBROUTINE vmax (npts, dt, leaf_age, leaf_frac, assim_param, &
       circ_class_biomass, circ_class_n, sugar_load, leaf_age_crit, &
       leaf_classes, veget_max)

    !! 0. Variable and parameter declaration
 
    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                                 :: npts                    !! Domain size (unitless)
    REAL(r_std), INTENT(in)                                    :: dt                      !! time step of stomate (days)
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(in)              :: circ_class_biomass      !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:),INTENT(in)                   :: circ_class_n            !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(in)                    :: sugar_load              !! Function reflecting how 
                                                                                          !! reserve+labile pools are filled
    REAL(r_std), DIMENSION(:,:), INTENT(in)                    :: veget_max               !! "Maximal" coverage fraction of a PFT (LAI 
                                                                                          !! -> infinity) on ground
   
    !! 0.2 Output variables 
    REAL(r_std), DIMENSION(:,:,:), INTENT(out)                 :: assim_param             !! vmax, nue and leaf N for photosynthesis
                                                                                          !! @tex $(\mu mol m^{-2}s^{-1})$ @endtex
   
    !! 0.3 Modified variables
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)               :: leaf_age                !! Leaf age (days)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)               :: leaf_frac               !! fraction of leaves in leaf age 
                                                                                          !! classes 
                                                                                          !! (unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)                 :: leaf_classes            !! width of each leaf age class (days)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)                 :: leaf_age_crit           !! critical leaf age (days)
    
    !! 0.4 Local variables
    REAL(r_std), DIMENSION(npts)                               :: leaf_efficiency         !! leaf efficiency (vcmax/vcmax25)
                                                                                          !! (unitless)
    REAL(r_std), DIMENSION(npts,nvm,nleafages)                 :: d_leaf_frac             !! turnover between age classes
                                                                                          !! (unitless)
    REAL(r_std), DIMENSION(npts,nleafages)                     :: leaf_age_new            !! new leaf age (days)
    REAL(r_std), DIMENSION(npts)                               :: sumfrac                 !! sum of leaf age fractions, 
                                                                                          !! for normalization
                                                                                          !! (unitless)
    REAL(r_std), DIMENSION(npts)                               :: rel_age                 !! relative leaf age (age/critical age)
                                                                                          !! (unitless)
    INTEGER(i_std)                                             :: j,m                     !! indices (unitless)
    REAL(r_std), DIMENSION(npts,nvm)                           :: vcmax                   !! Maximum rate of carboxylation 
                                                                                          !! @tex ($\mu mol m^{-2} s^{-1}$) @endtex
    REAL(r_std),DIMENSION (npts,nvm)                           :: nue                     !! Nitrogen use Efficiency with impact of leaf age (umol CO2 (gN)-1 s-1) 
!_ ================================================================================================================================

    IF (printlev>=2) WRITE(numout,*) 'Entering vmax'

    !! 1 Initialization

    !! 1.1 first call: info about flags and parameters.
    IF ( firstcall_vmax ) THEN
       
       IF (printlev >= 2) THEN
          WRITE(numout,*) 'vmax:'
          WRITE(numout,*) '   > offset (minimum vcmax/vmax_opt):' , vmax_offset
          WRITE(numout,*) '   > relative leaf age at which vmax reaches vcmax_opt:', leafage_firstmax 
          WRITE(numout,*) '   > relative leaf age at which vmax falls below vcmax_opt:', leafage_lastmax
          WRITE(numout,*) '   > relative leaf age at which vmax reaches its minimum:', leafage_old
       END IF

      firstcall_vmax = .FALSE.

    ENDIF

    ! Calculation how many days each leaf age class contains
    leaf_classes(:,:) = zero
    DO j = 2,nvm 
       WHERE (veget_max(:,j).GT.min_stomate)
          leaf_classes(:,j) = leaf_age_crit(:,j) / REAL(nleafages,r_std)
       END WHERE 
    END DO

    !! 2 leaf age: general increase and turnover between age classes.
   
    !! 2.1 increase leaf age
    !  The age of the leaves in each leaf-age-class increases by 1 time step.
    DO m = 1, nleafages ! Loop over # leaf age classes 
       DO j = 2,nvm ! Loop over # PFTs
          WHERE ( leaf_frac(:,j,m) .GT. min_stomate .AND. &
               veget_max(:,j) .GT. min_stomate)
             ! Increase leaf age
             leaf_age(:,j,m) = leaf_age(:,j,m) + dt            
          ENDWHERE
       ENDDO ! Loop over # PFTs 

    ENDDO ! Loop over # leaf age classes

    !! 2.2 turnover between leaf age classes
    !     d_leaf_frac(:,:,m) = what leaves m-1 and goes into m
    DO j = 2,nvm ! Loop over # PFTs

       !! 2.2.1 fluxes

       !! nothing goes into first age class
       d_leaf_frac(:,j,:) = zero

       !! for others age classes (what goes from m-1 to m)
       DO m = 2, nleafages 

          WHERE (veget_max(:,j).GT.min_stomate)
             ! leaf_classes is calculated in stomate_season.f90 as the quotient of the critical leaf age per the number of age classes.
             ! The critical leaf age is a PFT-dependent and calculated in stomate_season as a function of leaf_longevity and t2m_longterm, 
             ! representing the leaf life span.
             ! leaf_classes determines the turnover between the nleafages different leaf age classes
             ! (see section [118] in Krinner et al. (2005)).
             d_leaf_frac(:,j,m) = leaf_frac(:,j,m-1) * dt/leaf_classes(:,j)
          END WHERE

       ENDDO       

       !! 2.2.2 new leaf age in class
       !       new age = ( old age * (old fraction - fractional loss) + fractional increase * age of the source class ) / new fraction
       !       The leaf age of the youngest class (m=1) is updated into stomate_alloc          
       leaf_age_new(:,:) = zero

       DO m = 2, nleafages-1       ! Loop over age classes
          ! For all age classes except first and last 
          WHERE ( d_leaf_frac(:,j,m) .GT. min_stomate )

             leaf_age_new(:,m) = ( ( (leaf_frac(:,j,m)- d_leaf_frac(:,j,m+1)) * leaf_age(:,j,m) )  + &
                  ( d_leaf_frac(:,j,m) * leaf_age(:,j,m-1) ) ) / &
                  ( leaf_frac(:,j,m) + d_leaf_frac(:,j,m)- d_leaf_frac(:,j,m+1) )

          ENDWHERE

       ENDDO       ! Loop over age classes

       ! For last age class, there is no leaf fraction leaving the class. 
       WHERE ( d_leaf_frac(:,j,nleafages) .GT. min_stomate )

          leaf_age_new(:,nleafages) = ( ( leaf_frac(:,j,nleafages) * leaf_age(:,j,nleafages) )  + &
               ( d_leaf_frac(:,j,nleafages) * leaf_age(:,j,nleafages-1) ) ) / &
               ( leaf_frac(:,j,nleafages) + d_leaf_frac(:,j,nleafages) )

       ENDWHERE

       DO m = 2, nleafages       ! Loop over age classes

          WHERE ( d_leaf_frac(:,j,m) .GT. min_stomate )

             leaf_age(:,j,m) = leaf_age_new(:,m)

          ENDWHERE

       ENDDO       ! Loop over age classes

       !! 2.2.3 calculate new fraction
       DO m = 2, nleafages       ! Loop over age classes

          ! where the change comes from
          leaf_frac(:,j,m-1) = leaf_frac(:,j,m-1) - d_leaf_frac(:,j,m)

          ! where it goes to
          leaf_frac(:,j,m) = leaf_frac(:,j,m) + d_leaf_frac(:,j,m)

       ENDDO       ! Loop over age classes

       !! 2.2.4 renormalize fractions in order to prevent accumulation 
       !       of numerical errors

       ! correct small negative values
       DO m = 1, nleafages
          leaf_frac(:,j,m) = MAX( zero, leaf_frac(:,j,m) )
       ENDDO

       ! total of fractions, should be very close to one where there is leaf mass
       sumfrac(:) = zero

       DO m = 1, nleafages ! Loop over age classes

          sumfrac(:) = sumfrac(:) + leaf_frac(:,j,m)

       ENDDO ! Loop over age classes

       ! normalize
       DO m = 1, nleafages       ! Loop over age classes

          WHERE ( sumfrac(:) .GT. min_stomate )

             leaf_frac(:,j,m) = leaf_frac(:,j,m) / sumfrac(:) 

          ELSEWHERE

             leaf_frac(:,j,m) = zero

          ENDWHERE

       ENDDO       ! Loop over age classes

    ENDDO         ! Loop over PFTs

    !
    !! 3 calculate vmax as a function of the age
    !
   
    ! Define parameters for sugar-loading
    vcmax(:,:) = zero
    nue(:,:) = zero

    DO j = 2,nvm 

       ! sum up over the different age classes
       IF (ok_dgvm .AND. pheno_type(j)==1 .AND. leaf_tab(j)==2) THEN
          ! pheno_typ=evergreen and leaf_tab=needleleaf
          vcmax(:,j) = Vcmax25(j)

       ELSE
          
          DO m = 1, nleafages

             WHERE(veget_max(:,j).GT.min_stomate)
                !! 3.1 efficiency in each of the age classes
                !     it varies from vmax_offset to 1 
                !     linearly increases from vmax_offset to 1 for 0 < rel_age < leafage_firstmax
                !     is 1 when leafage_firstmax < rel_age < leafage_lastmax
                !     linearly decreases from 1 to vmax_offset for leafage_lastmax < rel_age < leafage_firstmax
                !     vmax_offset for rel_age >= leafage_old
                !     (Ishida et al., 1999)
                rel_age(:) = leaf_age(:,j,m) / leaf_age_crit(:,j)
                leaf_efficiency(:) = MAX( vmax_offset, MIN( un, &
                     vmax_offset + (un - vmax_offset) * rel_age(:) / leafage_firstmax, &
                     un - (un - vmax_offset) * ( rel_age(:) - leafage_lastmax ) / &
                     ( leafage_old - leafage_lastmax ) ) )


                !
                !! 3.2 add to mean vmax
                ! 
                ! vcmax calculated here is no longer used
                vcmax(:,j) = vcmax(:,j) + Vcmax25(j) * leaf_efficiency(:) * leaf_frac(:,j,m)

                ! Calculate nue. Nue is used in diffuco_trans_co2 to calculate vcmax
                ! from nue and Nleaf. So nue is the varaibable of interest
                nue(:,j) = nue(:,j) + nue_opt(j) * sugar_load(:,j) * &
                     leaf_efficiency(:) * leaf_frac(:,j,m)
                
             END WHERE
                
          ENDDO     ! loop over age classes

       ENDIF

    ENDDO       ! loop over PFTs

    !! 4. Photosynthesis parameters
    assim_param(:,:,ivcmax) = zero
    assim_param(:,:,inue)   = zero
    assim_param(:,:,ileafN) = zero
 
    ! CNP approach
    DO j = 2,nvm
       assim_param(:,j,ivcmax) = vcmax(:,j)
       assim_param(:,j,inue)   = nue(:,j)
       ! This approach distinguished between the structural N (0.4%) and the active
       ! nitrogen in the leaves. Note here leaf nitrogen per ground area is used.
       ! As our photosynthesis is calculated over different layers of LAI, we will
       ! need vcmax per leaf area per layer, thus requiring leaf nitrogen per leaf
       ! area per layer as well. This is done in the module where we calculate
       ! photosynthesis.
       assim_param(:,j,ileafn) = SUM(circ_class_biomass(:,j,:,ileaf,initrogen) * &
            circ_class_n(:,j,:),2) * (un - snc) 
    ENDDO

    CALL xios_orchidee_send_field("SUGAR_LOAD",sugar_load(:,:))
    CALL xios_orchidee_send_field("NUE",nue(:,:))

    IF (printlev>=4) WRITE(numout,*) 'Leaving vmax'

  END SUBROUTINE vmax

END MODULE stomate_vmax
