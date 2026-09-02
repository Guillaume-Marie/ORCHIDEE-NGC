! =================================================================================================================================
! MODULE       : stomate_somdynamics
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       Calculate soil dynamics largely following the Century model
!!	
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): None
!!
!! SVN		:
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_stomate/stomate_soilcarbon.f90 $ 
!! $Date: 2026-05-11 16:12:42 +0200 (lun. 11 mai 2026) $
!! $Revision: 9529 $
!! \n
!_ ================================================================================================================================

MODULE stomate_som_dynamics

  ! modules used:

  USE ioipsl_para
  USE stomate_data
  USE constantes
  USE constantes_soil
  USE xios_orchidee
  USE function_library, ONLY : check_mass_balance, get_printlev

  IMPLICIT NONE

  ! private & public routines

  PRIVATE
  PUBLIC som_dynamics,som_dynamics_clear,nitrogen_dynamics,nitrogen_dynamics_clear

  ! Variables shared by all subroutines in this module
  
  LOGICAL, SAVE                                           :: firstcall_soilcarbon = .TRUE.   !! First time a subroutine in current module is called
!$OMP THREADPRIVATE(firstcall_soilcarbon)
  LOGICAL, SAVE                                           :: firstcall_som = .TRUE.          !! First time subroutine som_dynamics is called
!$OMP THREADPRIVATE(firstcall_som)
  LOGICAL, SAVE                                           :: firstcall_nitrogen = .TRUE.     !! First time subroutine nitrogen_dynamics is called
!$OMP THREADPRIVATE(firstcall_nitrogen)
  INTEGER(i_std), SAVE                                    :: printlev_loc                    !! Local level of text output for current module
!$OMP THREADPRIVATE(printlev_loc)

  ! flux fractions within carbon pools
  REAL(r_std),ALLOCATABLE,SAVE, DIMENSION(:,:,:)          :: frac_carb        !! Flux fractions between carbon pools 
                                                                              !! (second index=origin, third index=destination) 
                                                                              !! (unitless, 0-1)
!$OMP THREADPRIVATE(frac_carb)
  REAL(r_std), ALLOCATABLE,SAVE, DIMENSION(:,:)           :: frac_resp        !! Flux fractions from carbon pools to the atmosphere (respiration) (unitless, 0-1)
!$OMP THREADPRIVATE(frac_resp)



CONTAINS


!! ================================================================================================================================
!!  SUBROUTINE   : som_dynamics_clear
!!
!>\BRIEF        Set the flag ::firstcall_som to .TRUE. and as such activate sections 1.1.2 and 1.2 of the subroutine som_dynamics 
!! (see below).
!! 
!_ ================================================================================================================================
  
  SUBROUTINE som_dynamics_clear
    firstcall_som=.TRUE.
  END SUBROUTINE som_dynamics_clear


!! ================================================================================================================================
!!  SUBROUTINE   : som_dynamics
!!
!>\BRIEF        Computes the soil respiration and nutrient stocks, essentially 
!! following Parton et al. (1987).  Additional dynamics for the nitrogen pools are
!! described in the nitrogen_dynamics subroutine.
!!
!! DESCRIPTION	: The soil is divided into 3 pools, with different 
!! characteristic turnover times : active (1-5 years), slow (20-40 years) 
!! and passive (200-1500 years).\n
!! There are three types of nutrient (carbon, nitrogen) transferred into the soil:\n
!! - input to active and slow pools from litter decomposition,\n
!! - nutrient fluxes between the three pools,\n
!! - nurtient losses from the pools to the atmosphere, i.e., soil respiration.\n
!!
!! The subroutine performs the following tasks:\n
!!
!! Section 1.\n
!! The flux fractions (f) between carbon pools are defined based on Parton et 
!! al. (1987). The fractions are constants, except for the flux fraction from
!! the active pool to the slow pool, which depends on the clay content,\n
!! \latexonly
!! \input{soilcarbon_eq1.tex}
!! \endlatexonly\n
!! In addition, to each pool is assigned a constant turnover time.  No nitrogen
!! is considered in this section.\n
!!
!! Section 2.\n
!! The carbon and nitrogen inputs, calculated in the stomate_litter module, are 
!! added to the carbon and nitrogen stocks of the different pools.\n
!!
!! Section 3.\n
!! First, the outgoing flux of each pool is calculated. It is 
!! proportional to the product of the carbon stock and the ratio between the 
!! iteration time step and the residence time:\n
!! \latexonly
!! \input{soilcarbon_eq2.tex}
!! \endlatexonly
!! ,\n
!! Note that in the case of crops, the additional multiplicative factor 
!! integrates the faster decomposition due to tillage (following Gervois et 
!! al. (2008)).
!! In addition, the flux from the active pool depends on the clay content:\n
!! \latexonly
!! \input{soilcarbon_eq3.tex}
!! \endlatexonly
!! ,\n
!! Each pool is then cut from the carbon amount corresponding to each outgoing
!! flux:\n
!! \latexonly
!! \input{soilcarbon_eq4.tex}
!! \endlatexonly\n
!! Note that the total flux for both carbon and nitrogen out of the pools is
!! calculated first and grouped together.  This total material flux is then
!! partitioned between the elements and the pools based on a target CN ratio.\n
!! Second, the flux fractions lost to the atmosphere is calculated in each pool
!! by subtracting from 1 the pool-to-pool flux fractions. The soil respiration 
!! is then the summed contribution of all the pools,\n
!! \latexonly
!! \input{soilcarbon_eq5.tex}
!! \endlatexonly\n
!! Note that soil respiration only happens for the carbon pools.\n
!! Finally, each soil nutrient pool accumulates the contribution of the other pools:
!! \latexonly
!! \input{soilcarbon_eq6.tex}
!! \endlatexonly
!! The lefotver nitrogen flux that doesn't go into one of the four pools is 
!! mineralized.\n
!! Section 4.\n
!! If the flag SPINUP_ANALYTIC is set to true, the matrix A is updated following
!! Lardy (2011).
!!
!! RECENT CHANGE(S): Merge with the Nitrogen cycle, 2018.
!! 
!! MAIN OUTPUTS VARIABLE(S): carbon, resp_hetero_soil
!!
!! REFERENCE(S)   :
!! - Parton, W.J., D.S. Schimel, C.V. Cole, and D.S. Ojima. 1987. Analysis of 
!!   factors controlling soil organic matter levels in Great Plains grasslands. 
!!   Soil Sci. Soc. Am. J., 51, 1173-1179.
!! - Gervois, S., P. Ciais, N. de Noblet-Ducoudre, N. Brisson, N. Vuichard, 
!!   and N. Viovy (2008), Carbon and water balance of European croplands 
!!   throughout the 20th century, Global Biogeochem. Cycles, 22, GB2022, 
!!   doi:10.1029/2007GB003018.
!! - Lardy, R, et al., A new method to determine soil organic carbon equilibrium,
!!   Environmental Modelling & Software (2011), doi:10.1016|j.envsoft.2011.05.016
!! - S. Zaehle and A. D. Friend (2010), Carbon and nitrogen cycle dynamics in the
!!   O-CN land surface model: 1. Model description, site-scale evaluation, and 
!!   sensitivity to parameter estimates. Global Biogeochem. Cycles, 24, GB1005, 
!!   doi:10.1029/2009GB003521.
!!
!! FLOWCHART    :
!! \latexonly
!! \includegraphics[scale=0.5]{soilcarbon_flowchart.jpg}
!! \endlatexonly
!! \n
!_ ================================================================================================================================

  SUBROUTINE som_dynamics (npts, lalo, clay, silt, veget_max, & 
       som_input, control_temp, control_moist, drainage_pft,&
       CN_target, som, soil_n_min, resp_hetero_soil, &
       MatrixA, n_mineralisation, CN_som_litter_longterm, tau_CN_longterm, &
       height_acro,height_cato,carbon_acro,carbon_cato,tcarbon_acro,tcarbon_cato,resp_acro_oxic, &
       resp_acro_anoxic,resp_cato,acro_to_cato,litter_to_acro, wtp_peat)

!! 0. Variable and parameter declaration

    !! 0.1 Input variables
    
    INTEGER(i_std), INTENT(in)                            :: npts             !! Domain size (unitless)
    REAL(r_std), DIMENSION(:,:),INTENT(in)                :: lalo             !! Geogr. coordinates (latitude,longitude) (degrees)
    REAL(r_std), DIMENSION(:), INTENT(in)                 :: clay             !! Clay fraction (unitless, 0-1) 
    REAL(r_std), DIMENSION(:), INTENT(in)                 :: silt             !! Silt fraction (unitless, 0-1)
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)           :: som_input        !! Amount of Organic Matter going into the SOM pools from litter decomposition \f$(gC m^{-2} day^{-1})$\f
    REAL(r_std), DIMENSION(:,:), INTENT(in)               :: control_temp     !! Temperature control of heterotrophic respiration (unitless: 0->1)
    REAL(r_std), DIMENSION(:,:), INTENT(in)               :: control_moist    !! Moisture control of heterotrophic respiration (unitless: 0.25->1)
    REAL(r_std), DIMENSION(:,:), INTENT(in)               :: veget_max        !! Fractional coverage: maximum share of the pixel taken by a pft 
    REAL(r_std), DIMENSION(:,:), INTENT(in)               :: drainage_pft     !! Drainage per PFT (mm/m2 /dt_sechiba) (fraction of water content)   
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)             :: CN_target        !! C to N ratio of SOM flux from one pool to another (gN m-2 dt-1)
    REAL(r_std), DIMENSION(:,:,:),INTENT(in)              :: soil_n_min       !! mineral nitrogen in the soil (gN/m**2) 
    REAL(r_std), DIMENSION(npts),INTENT(in)               :: wtp_peat         !! water table position when peat is activated

    !! 0.2 Output variables
    
    REAL(r_std), DIMENSION(:,:), INTENT(out)              :: resp_hetero_soil !! Soil heterotrophic respiration \f$(gC m^{-2} (dt_sechiba one_day^{-1})^{-1})$\f

    ! peatland: they will be deleted later, and comment not fully added.
    REAL(r_std), DIMENSION(npts),INTENT(out)               :: tcarbon_acro
    REAL(r_std), DIMENSION(npts),INTENT(out)               :: tcarbon_cato
    REAL(r_std), DIMENSION(npts,nvm,nelements),INTENT(out) :: resp_acro_oxic   !!respiration of acrotelm( oxic )
    REAL(r_std), DIMENSION(npts,nvm,nelements),INTENT(out) :: resp_acro_anoxic !!respiration of acrotelm( anoxic )
    REAL(r_std), DIMENSION(npts,nvm,nelements),INTENT(out) :: resp_cato        !!respiration of catotelm
    REAL(r_std), DIMENSION(npts,nvm,nelements),INTENT(out) :: litter_to_acro   !!total carbon input from litter to acrotelm
    REAL(r_std), DIMENSION(npts,nvm,nelements),INTENT(out) :: acro_to_cato     !!carbon flux from acrotelm to catotelm

    REAL(r_std), DIMENSION(npts), INTENT(out)              :: height_cato


    !! 0.3 Modified variables
    
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)        :: som              !! SOM pools: active, slow, or passive, \f$(gC m^{2})$\f
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)        :: MatrixA          !! Matrix containing the fluxes between the carbon pools
                                                                              !! per sechiba time step 
                                                                              !! @tex $(gC.m^2.day^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)            :: n_mineralisation !! net nitrogen mineralisation of decomposing SOM,
                                                                              !! (gN/m**2/day), assumed to be NH4
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)          :: CN_som_litter_longterm !! Longterm CN ratio of litter and som pools (gC/gN)
    REAL(r_std), INTENT(inout)                            :: tau_CN_longterm  !! Counter used for calculating the longterm CN ratio of SOM and litter pools (seconds)

    
    REAL(r_std), DIMENSION(npts,nvm,nelements), INTENT(inout)              ::carbon_acro
    REAL(r_std), DIMENSION(npts,nvm,nelements), INTENT(inout)              ::carbon_cato
    REAL(r_std), DIMENSION(npts), INTENT(inout)                   ::height_acro

    !! 0.4 Local variables
    REAL(r_std)                                           :: dt               !! Time step \f$(dt_sechiba one_day^{-1})$\f
    LOGICAL                                               :: l_error          !! Diagnostic boolean for error allocation (true/false)
    INTEGER(i_std)                                        :: ier              !! Check errors in netcdf call
    REAL(r_std), DIMENSION(ncarb,nvm)                     :: som_turn         !! Residence time in SOM pools (days)
    REAL(r_std), DIMENSION(npts,nvm,ncarb,nelements)      :: fluxtot          !! Total flux out of carbon pools \f$(gC m^{2})$\f
    REAL(r_std), DIMENSION(npts,nvm,ncarb,ncarb,nelements):: flux             !! Fluxes between carbon pools \f$(gC m^{2})$\f    
    REAL(r_std), DIMENSION(npts,ncarb,ncarb,nelements)    :: distributed      !!
    REAL(r_std), DIMENSION(npts,ncarb,nelements)          :: residual         !!
    
    REAL(r_std), DIMENSION(npts,nvm,nelements)            :: total_som_n_flux !! Total flux between all soil pools. Used 
                                                                              !! to adjust som if overspending 
                                                                              !! @tex $(N gN m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)            :: old_som_n_flux   !! Total flux between all soil pools. Used 
                                                                              !! to adjust som if overspending 
                                                                              !! @tex $(N gN m^{-2})$ @endtex
    REAL(r_std)                                           :: delta_n_min      !! Amount of nitrogen that was overspend 
                                                                              !! @tex $(N gN m^{-2})$ @endtex
    REAL(r_std)                                           :: reduction_factor !! Factor to reduce the som pools in case of 
                                                                              !! overspending (unitless)
    REAL(r_std), DIMENSION(ncarb,nelements)               :: share_flux       !! share of the diiferent som components in case of
                                                                              !! overspending (unitless)
    REAL(r_std), DIMENSION(nelements)                     :: delta_som        !! Absolute reduction gN m-2 or gC m-2 in case of overspending
    CHARACTER(LEN=7), DIMENSION(ncarb)                    :: soilpools_str    !! Name of the soil pools for informative outputs (unitless)
    REAL(r_std), DIMENSION(npts,nvm,ncarb)                :: decomp_rate_soilcarbon !! Decomposition rate of the soil carbon pools (s) 
    REAL(r_std), DIMENSION(npts,ncarb)                    :: tsoilpools       !! Diagnostic for soil carbon turnover rate by pool (1/s)
    INTEGER(i_std)                                        :: iicarb           !! Indices (unitless)
    INTEGER(i_std)                                        :: icarb, imbc      !! Indices (unitless)
    INTEGER(i_std)                                        :: ipts, ivm,iele   !! Indices (unitless)
    INTEGER(i_std)                                        :: inbpools         !! Indices (unitless)
    REAL(r_std), DIMENSION(npts,nvm,nmbcomp,nelements)    :: check_intern     !! Contains the components of the internal
                                                                              !! mass balance check for this routine
                                                                              !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)            :: closure_intern   !! Check closure of internal mass balance
                                                                              !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)            :: pool_start       !! Start and end pool of this routine 
                                                                              !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)            :: pool_end         !! Start and end pool of this routine 
                                                                              !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    INTEGER(i_std)                                        :: ix               !! Row/column number of the maximal value
    REAL(r_std), DIMENSION(npts,ncarb,nvm,nelements)      :: som_start        !! som at the start of the subroutine
    REAL(r_std), DIMENSION(npts,nvm)                      :: n_mineralisation_start !! n_mineralisation at the start of the subroutine
    REAL(r_std)                                           :: mismatch         !! Carbon and nitrogen lost/gained during the calculations (precision error)
    REAL(r_std), DIMENSION(npts)                          :: error_count     !! Count the number of errors

    ! peatland 
    REAL(r_std), DIMENSION(npts)                          ::B                !! acrotelm oxic respiration 
    REAL(r_std), DIMENSION(npts)                          ::KA               !! acrotelm decomposition rate
    REAL(r_std), DIMENSION(npts)                          ::KP               !! catotelm formation rate
    REAL(r_std), DIMENSION(npts)                          ::KC               !! catotelm decomposition rate
    INTEGER(i_std)                                        ::pft_peat
    REAL(r_std)                                           ::wtd_min          !! interface between acrotelm and catotelm
    REAL(r_std), DIMENSION(npts)                          ::wtd              !! water table depth
!_ ================================================================================================================================

    !! Initialize printlev_loc only first time one of the subroutines 
    IF ( firstcall_soilcarbon ) THEN
      ! Guillaume M. -- get_printlev calls getin_p, a COLLECTIVE bcast. Called here,
      ! inside the compute loop under a local firstcall, it hangs the run as soon as the
      ! ranks stop reaching it together. Use printlev, read once at startup by every rank.
       printlev_loc=printlev
       firstcall_soilcarbon = .FALSE.
    ENDIF

    IF (printlev_loc>=3) WRITE(numout,*) 'Entering som_dynamics' 

    dt = dt_sechiba/one_day

    !! 0. Initialize check for mass balance closure
    !  The mass balance is calculated at the end of this routine
    !  in section 14. Initial biomass and harvest pool and all other
    !  relevant pools were set to zero before calling this routine.
    pool_start(:,:,:) = zero
    DO iele = 1,nelements
       DO icarb = 1,ncarb
          pool_start(:,:,iele) = pool_start(:,:,iele) + &
               ( som(:,icarb,:,iele) + (som_input(:,icarb,:,iele) * dt) ) * &
               veget_max(:,:)
       ENDDO
    ENDDO

    pool_start(:,:,initrogen) = pool_start(:,:,initrogen) + &
         ( n_mineralisation(:,:) ) * veget_max(:,:)
    

    !! 1. Initializations

    IF ( firstcall_som ) THEN
    !! 1. Initializations
        l_error = .FALSE.
        ALLOCATE (frac_carb(npts,ncarb,ncarb), stat=ier)
        l_error = l_error .OR. (ier.NE.0)
        ALLOCATE (frac_resp(npts,ncarb), stat=ier)
        l_error = l_error .OR. (ier.NE.0)
        IF (l_error) THEN
           STOP 'stomate_som_dynamics: error in memory allocation'
        ENDIF

        frac_carb(:,:,:) = 0.0
        !! 1.1 Get soil "constants"
        !! 1.1.1 Flux fractions between carbon pools: depend on clay content, recalculated each time
        ! From active pool: depends on clay content
        frac_carb(:,iactive,ipassive) = active_to_pass_ref_frac + active_to_pass_clay_frac*clay(:)
        frac_carb(:,iactive,islow) = un - frac_carb(:,iactive,ipassive) - (active_to_co2_ref_frac - &
             active_to_co2_clay_silt_frac*(clay(:)+silt(:)))
    
        ! From slow pool    
        frac_carb(:,islow,ipassive) = slow_to_pass_ref_frac + slow_to_pass_clay_frac*clay(:)
        ! OCN doesn't use Parton 1993 formulation for frac_carb(:,islow,ipassive) 
        ! but the one of 1987 : ie = 0.03 .... 
        frac_carb(:,islow,iactive) = un - frac_carb(:,islow,ipassive) - slow_to_co2_ref_frac
    
        ! From passive pool
        frac_carb(:,ipassive,iactive) = pass_to_active_ref_frac
        frac_carb(:,ipassive,islow) = pass_to_slow_ref_frac

        ! From surface pool
        frac_carb(:,isurface,islow) = surf_to_slow_ref_frac

        !! 1.1.2 Determine the respiration fraction : what's left after
        ! subtracting all the 'pool-to-pool' flux fractions
        ! Diagonal elements of frac_carb are zero
        frac_resp(:,:) = un - frac_carb(:,:,isurface) - frac_carb(:,:,iactive) - frac_carb(:,:,islow) - &
             frac_carb(:,:,ipassive) 
        
        !! 1.2 Messages : display the residence times  
        soilpools_str(iactive) = 'active'
        soilpools_str(islow) = 'slow'
        soilpools_str(ipassive) = 'passive'
        soilpools_str(isurface) = 'surface'

       IF(printlev_loc.GT.2)THEN
          WRITE(numout,*) 'som_dynamics:'        
          WRITE(numout,*) '   > minimal SOM residence time in soil pools (d):'
          DO icarb = 1, ncarb ! Loop over soil pools
             WRITE(numout,*) '(1, ::soilpools_str(icarb)):', &
                  soilpools_str(icarb)
             WRITE(numout,*) 'FRACCARB icarb=',icarb,' ', &
                  frac_carb(1,icarb,isurface),frac_carb(1,icarb,iactive),&
                  frac_carb(1,icarb,islow),frac_carb(1,icarb,ipassive)
          ENDDO
          WRITE(numout,*) '   > flux fractions between soil pools: depend on clay content'
       ENDIF



        firstcall_som = .FALSE.
        
    ENDIF

    ! som_turn(ipassive,ivm) now uses PFT-dependent parameters. The others, iactive,
    ! islow and isurface have a PFT-dimension but the parameter value is still the
    ! same for all PFTs
    DO ivm = 1,nvm
       !! 1.1.3 Turnover in SOM pools (in days)
       !! som_turn_ipool are the turnover (in year)
       !! It is weighted by Temp and Humidity function later
       som_turn(iactive,ivm) = som_turn_iactive / one_year    
       som_turn(islow,ivm) = som_turn_islow / one_year                
       som_turn(isurface,ivm) = som_turn_isurface / one_year
       ! Note that som_turn_ipassive is defined at the PFT level as an
       ! easy way to tune the total soil carbon without disrupting the
       ! N-cycle too much. This is considered a hack, thus not a
       ! permanent solution.
       som_turn(ipassive,ivm) = som_turn_ipassive(ivm) / one_year
    END DO
    
    
    !! 1.3 Initialise
    resp_hetero_soil(:,:) = zero
    decomp_rate_soilcarbon(:,:,:) = zero
    total_som_n_flux(:,:,:) = zero
    old_som_n_flux(:,:,:) = zero

!! 2. Update the SOM stocks with the different soil carbon and nitrogen input

    IF(printlev_loc.GT.4)THEN
       WRITE(numout,*) 'som, ',som(:,:,:,:)
       WRITE(numout,*) 'som_input, ',som_input(:,:,:,:)
       WRITE(numout,*) 'som_turn, ',som_turn(:,:)
       WRITE(numout,*) 'decomp_factor, ',decomp_factor(:)
    ENDIF

    som(:,:,:,:) = som(:,:,:,:) + som_input(:,:,:,:) * dt 

    ! Store value for high precision mass conservation 
    som_start(:,:,:,:) = som(:,:,:,:)
    n_mineralisation_start(:,:) = n_mineralisation(:,:)

!! 3. Fluxes between nutrient reservoirs, and to the atmosphere (respiration) \n

    !! 3.2. Calculate fluxes
    DO ivm = 1,nvm ! Loop over # PFTs

       !! 3.2.1. Flux out of pools
       ! Loop over carbon pools from which the flux comes
       !
       IF (.NOT. is_peat(ivm)) THEN  ! this line of condition was introduced by Chunjing. 
          DO icarb = 1, ncarb 
             
             DO iele = 1, nelements 
                
                ! Determine total flux out of pool [gC/N m-2]
                fluxtot(:,ivm,icarb,iele) = dt*som_turn(icarb,ivm) * &
                     som(:,icarb,ivm,iele) * control_moist(:,ibelow) * &
                     control_temp(:,ibelow) * decomp_factor(ivm)
                
                ! Flux from active pools depends on clay content [gC/N m-2]
                IF ( icarb .EQ. iactive ) THEN
                   fluxtot(:,ivm,icarb,iele) = fluxtot(:,ivm,icarb,iele) * &
                        ( un - som_turn_iactive_clay_frac * clay(:) )
                ENDIF
                
                ! Update the loss in each carbon pool [gC/N m-2]. Because som
                ! is calculated as the residual of som and fluxtot mass should be
                ! preserved with a very high precision.
                som(:,icarb,ivm,iele) = som(:,icarb,ivm,iele) - &
                     fluxtot(:,ivm,icarb,iele)
                
                ! Calculate diagnostic for decomposition rate of the soil carbon pools
                IF ( iele == icarbon ) THEN
                   decomp_rate_soilcarbon(:,ivm,icarb) = dt*som_turn(icarb,ivm) * &
                        control_moist(:,ibelow) * control_temp(:,ibelow) * decomp_factor(ivm)
                   IF ( icarb .EQ. iactive ) THEN
                      decomp_rate_soilcarbon(:,ivm,icarb) = decomp_rate_soilcarbon(:,ivm,icarb) * &
                           ( un - som_turn_iactive_clay_frac * clay(:) )
                   ENDIF
                END IF
                
             ENDDO !nelements
             
             ! Fluxes towards the other pools (icarb -> iicarb)
             ! Loop over the SOM pools where the flux goes
             DO iicarb = 1, ncarb 
                
                ! Carbon flux [gN m-2]
                flux(:,ivm,icarb,iicarb,icarbon) = frac_carb(:,icarb,iicarb) * &
                     fluxtot(:,ivm,icarb,icarbon)
                
                ! Nitrogen flux - Function of the C stock of the 'departure' pool 
                ! and of the C to N target ratio of the 'arrival' pool [gN m-2]
                flux(:,ivm,icarb,iicarb,initrogen) = frac_carb(:,icarb,iicarb) * &
                     fluxtot(:,ivm,icarb,icarbon) / & 
                     CN_target(:,ivm,iicarb)
                
             ENDDO ! receiving SOM pools
             
          ENDDO ! donating SOM pools

          !! 3.2.2 respiration
          !  Fluxtot is in [gC m-2], resp_hetero_soil is in [gC m-2 dt-1] 
          !  divide by dt.
          resp_hetero_soil(:,ivm) = &
               ( frac_resp(:,iactive) * fluxtot(:,ivm,iactive,icarbon) + &
               frac_resp(:,islow) * fluxtot(:,ivm,islow,icarbon) + &
               frac_resp(:,ipassive) * fluxtot(:,ivm,ipassive,icarbon) + &
               frac_resp(:,isurface) * fluxtot(:,ivm,isurface,icarbon) ) / dt
          
          !! 3.2.3 add fluxes to active, slow, and passive pools
          ! Loop over SOM pools
          DO icarb = 1, ncarb 
             
             ! Calculate the components of the som fluxes [gN m-2]
             ! or [gC m-2]
             som(:,icarb,ivm,:) = som(:,icarb,ivm,:) + &
                  flux(:,ivm,iactive,icarb,:) + &
                  flux(:,ivm,ipassive,icarb,:) + &
                  flux(:,ivm,islow,icarb,:) + &
                  flux(:,ivm,isurface,icarb,:)
             
             ! Calculate the total flux [gN m-2]
             total_som_n_flux(:,ivm,:) = total_som_n_flux(:,ivm,:) + &
                  flux(:,ivm,iactive,icarb,:) + &
                  flux(:,ivm,ipassive,icarb,:) + &
                  flux(:,ivm,islow,icarb,:) + &
                  flux(:,ivm,isurface,icarb,:)
             
          ENDDO ! Loop over SOM pools

          !peatland++
          ! added one dimension of nelements, xw
          IF (ok_peat_NoDiscretisation) THEN      
             carbon_acro(:,ivm,:)=zero
             carbon_cato(:,ivm,:)=zero
             acro_to_cato(:,ivm,:)=zero
             litter_to_acro(:,ivm,:)=zero
             resp_cato(:,ivm,:)=zero
             resp_acro_anoxic(:,ivm,:)=zero
             resp_acro_oxic(:,ivm,:)=zero
          ENDIF
       ENDIF  ! is_peat 
       
    ENDDO ! End loop over PFTs

    ! In the ORC3, n_mineralisation is calculated in the loop above as:
    ! n_mineralisation(:,ivm) = n_mineralisation(:,ivm) + fluxtot(:,k,initrogen) - &
    !   (flux(:,k,iactive,initrogen)+flux(:,k,ipassive,initrogen)+&
    !   flux(:,k,islow,initrogen)+flux(:,k,isurface,initrogen)). In CN-CAN we
    ! really want to avoid negative n_mineralisation so the code below is to ensure
    ! n_mineralisation will remain positive.
    
    ! Store to calculate absolute overspending
    old_som_n_flux(:,:,:) = total_som_n_flux(:,:,:)

    ! Choose between two approaches to deal with mass balance problems: close
    ! the mass balance at the end of the calculations or close it in every
    ! subroutine seperatly.
    IF (.NOT.test_adjust_rhetero) THEN

       DO ivm = 1,nvm
          
          DO ipts = 1,npts

             ! Follow the ORC3 approach
             ! Calculate the required mineralisation, possible mass issues
             ! will be dealt with in nitrogen_dynamics.
             n_mineralisation(ipts,ivm) = n_mineralisation(ipts,ivm) + &
                  SUM(fluxtot(ipts,ivm,:,initrogen)) - &
                  total_som_n_flux(ipts,ivm,initrogen)
             
          END DO
    
       END DO

    ELSE
       
       ! If the proposed decomposition consumes more nitrogen
       ! then currently available in the soil_n_min pool, the
       ! decomposition needs to be adjusted/reduced.
       DO ivm = 1,nvm

          DO ipts = 1,npts

             IF ((soil_n_min(ipts,ivm,initrate)+soil_n_min(ipts,ivm,iammonium)) + &
                  n_mineralisation(ipts,ivm) + SUM(fluxtot(ipts,ivm,:,initrogen)) - &
                  total_som_n_flux(ipts,ivm,initrogen) .GE. zero) THEN

                ! There is enough nitrogen in the mineralised pool
                ! to allow for the immobiolization (if any) in 
                ! soil_carbon [gN m-2]. This is the same as what is
                ! done in the trunk.
                n_mineralisation(ipts,ivm) = n_mineralisation(ipts,ivm) + &
                     SUM(fluxtot(ipts,ivm,:,initrogen)) - &
                     total_som_n_flux(ipts,ivm,initrogen)

             ELSE

                ! There is not enough soil nitrogen to immobilize later on. 
                ! Need to decrease the som decomposition. Calculate the 
                ! amount of nitrogen [gN m-2] that was overspend [gN m-2]
                delta_n_min = SUM(fluxtot(ipts,ivm,:,initrogen)) - &
                     total_som_n_flux(ipts,ivm,initrogen) + &
                     ((soil_n_min(ipts,ivm,initrate)+soil_n_min(ipts,ivm,iammonium)) + &
                     n_mineralisation(ipts,ivm))

                ! Because total_som_n_flux is subtracted, the most 
                ! effective method is adjust this variable. However,
                ! if we decrease total_som_n_flux, we should also decrease
                ! fluxtot proportionally. 
                reduction_factor = (SUM(fluxtot(ipts,ivm,:,initrogen)) - &
                     total_som_n_flux(ipts,ivm,initrogen)-delta_n_min)/&
                     (SUM(fluxtot(ipts,ivm,:,initrogen)) - &
                     total_som_n_flux(ipts,ivm,initrogen))

                ! Return the nitrogen and carbon that was decomposed in excess
                ! of the available nitrogen pools.
                delta_som(:) = (SUM(fluxtot(ipts,ivm,:,:),1) - &
                     total_som_n_flux(ipts,ivm,:))*(1-reduction_factor)
                fluxtot(ipts,ivm,:,:) = fluxtot(ipts,ivm,:,:)*reduction_factor
                total_som_n_flux(ipts,ivm,initrogen) = reduction_factor * &
                     total_som_n_flux(ipts,ivm,initrogen)

                ! Recalculate n_mineralisation with the reduced fluxes
                n_mineralisation(ipts,ivm) = n_mineralisation(ipts,ivm) + &
                     SUM(fluxtot(ipts,ivm,:,initrogen)) - &
                     total_som_n_flux(ipts,ivm,initrogen)

                ! Error checking
                IF (ABS(delta_som(initrogen)-delta_n_min).GT.min_stomate)THEN

                   WRITE(numout,*) 'ERROR: inconsistency with nitrogen reduction',&
                        delta_som(initrogen),delta_n_min
                   CALL ipslerr_p(3,'som_dynamics',&
                        'inconsistency in calculating nitrogen overspending','','')

                ENDIF

                ! Error checking
                IF ((soil_n_min(ipts,ivm,initrate)+soil_n_min(ipts,ivm,iammonium)) + &
                     n_mineralisation(ipts,ivm) .LT. SUM(fluxtot(ipts,ivm,:,initrogen)) - &
                     total_som_n_flux(ipts,ivm,initrogen)) THEN

                   WRITE(numout,*) 'ERROR: failed to solve nitrogen overspending',&
                        (soil_n_min(ipts,ivm,initrate)+soil_n_min(ipts,ivm,iammonium))+&
                        n_mineralisation(ipts,ivm), SUM(fluxtot(ipts,ivm,:,initrogen)) - &
                        total_som_n_flux(ipts,ivm,initrogen)
                   CALL ipslerr_p(3,'som_dynamics','Still overspending',&
                        'Something wrong with the calculations of the reduction','')

                ENDIF

                ! Proportionaly distribue delta_som over the flux components
                share_flux(:,:) = zero
                DO iele = 1,nelements

                   DO icarb = 1,ncarb

                      share_flux(icarb,iele) = SUM(flux(ipts,ivm,:,icarb,iele)) / &
                           SUM(SUM(flux(ipts,ivm,:,:,iele),1))

                   ENDDO

                ENDDO

                ! Avoid introducing mass balance problems. Calculate the last share 
                ! as the residual of the other shares to make sure that SUM(share)=1
                share_flux(ncarb,:) = 1
                DO icarb = 1,ncarb-1

                   share_flux(ncarb,:) = share_flux(ncarb,:) - share_flux(icarb,:)

                ENDDO

                ! Error checking
                DO iele = 1,nelements

                   IF (ABS(SUM(share_flux(:,iele))-1).GT.min_stomate) THEN

                      WRITE(numout,*) 'ERROR: sum of shares should equal 1',&
                           iele, SUM(share_flux(:,iele))
                      CALL ipslerr_p(3,'som_dynamics','sum of shares should equal','','')

                   ENDIF

                ENDDO

                ! Restore the som pools. The implicit dimension denotes the carbon pools.
                DO iele = 1,nelements

                   som(ipts,:,ivm,iele) = som(ipts,:,ivm,iele) + &
                        delta_som(iele) * share_flux(:,iele)
                   
                ENDDO

                ! Some of the nitrogen that was initially mineralized is no longer 
                ! mineralized. A proportional amount of C has been restored in the 
                ! soil pools. This implies that there was less decomposition than 
                ! initially estimated so heterotrophic respiration also needs to be
                ! adjusted.

                IF (resp_hetero_soil(ipts,ivm)-delta_som(icarbon)/dt.GT.-min_stomate)THEN

                   ! Adjust heterotrophic respiration
                   resp_hetero_soil(ipts,ivm) = resp_hetero_soil(ipts,ivm) - &
                        delta_som(icarbon)/dt

                ELSE
                   
                   ! The decomposition was adjusted but should never become negative
                   WRITE(numout,*) 'ERROR: adjusted heterotrophic respiration &
                        & became negative',resp_hetero_soil(ipts,ivm), &
                        delta_som(icarbon)/dt,resp_hetero_soil(ipts,ivm)-delta_som(icarbon)/dt
                   WRITE(numout,*) "ipts, lon, lat: ",ipts, lalo(ipts,2), lalo(ipts,1)
                   WRITE(numout,*) "ivm: ",ivm
                   WRITE(numout,*) "If you got this error message during a spinup in coupled mode, "
                   WRITE(numout,*) "then you can perturbate the atmosphere slightly and relaunch again"
                   CALL ipslerr_p(3,'stomate_soilcarbon',&
                        'adjusted heterotrophic respiration became negative','','')
                   
                ENDIF ! resp_hetero_litter .GT. SUM(delta_som_input)

             ENDIF ! overspending nitrogen?
          ENDDO !nvm

       ENDDO !npts

    END IF ! test_adjust_rhetero

    !! peatland ++ with question:
    ! I put the following added part before the mass balance check in the trunk. Hope that it is fine.
    ! ok_peat_NoDiscretisation !!Flag to activite peatland carbon accumumodule, from Kleinen et al,2012
    ! this flag is False par default. xw 202506
    IF (ok_peat_NoDiscretisation) THEN
       !!! decomposition rates and catotelm formation rates are constrained by temperature
       KA(:)=control_temp(:,ibelow)*KA_ini*dt/one_year !control_temp(:,ibelow): integrated temperature of of 11layers
       KP(:)=control_temp(:,ibelow)*KP_ini*dt/one_year
       KC(:)=KC_ini*dt/one_year   !*control_temp(:,ibelow)?

       height_cato(:)=zero
       tcarbon_acro(:) = zero
       tcarbon_cato(:) = zero
       resp_acro_oxic(:,:,:)=zero
       resp_acro_anoxic(:,:,:)=zero
       resp_cato(:,:,:)=zero
       litter_to_acro(:,:,:)=zero
       acro_to_cato(:,:,:)=zero
       wtd_min=300.  !in mm
       CALL getin_p('wtd_min', wtd_min)

       DO ipts=1,npts
           wtd(ipts)=wtd_min-wtp_peat(ipts) !in mm
           wtd(ipts)=wtd(ipts)*0.001  ! in m 
           IF (wtd(ipts) .LE. 0.0) THEN
              B(ipts)=1.0
           ELSE IF (wtd(ipts) .LT. height_acro(ipts)) THEN
              B(ipts)=(height_acro(ipts)-wtd(ipts))/height_acro(ipts)
           ELSE IF (wtd(ipts) .GE. height_acro(ipts)) THEN
              B(ipts)=0.0
           ENDIF
       ENDDO
       
       ! Question : add a new dimension of nelements in some variables below. need to make sure that they are fine. XW.
       ! nelements = 2    !! Number of isotopes considered
       DO ipts=1,npts
          DO ivm=2,nvm
             DO iele = 1,nelements
                IF (is_peat(ivm)) THEN
                   
                   ! modofied to som for trunk, xw, by adding one more dimension
                   ! som has dimension of npts, ncarb, ivm, nelement
                   som(ipts,:,ivm, iele) =zero
               
                   ! modified for trunk xw
                   resp_acro_oxic(ipts,ivm,iele)=resp_acro_oxic(ipts,ivm,iele)+B(ipts)*KA(ipts)*carbon_acro(ipts,ivm,iele)
                   resp_acro_anoxic(ipts,ivm,iele)=resp_acro_anoxic(ipts,ivm,iele)+(1.-B(ipts))*v_ratio*KA(ipts)*carbon_acro(ipts,ivm,iele)
                   acro_to_cato(ipts,ivm,iele)=acro_to_cato(ipts,ivm,iele)+KP(ipts)*carbon_acro(ipts,ivm,iele)
                   resp_cato(ipts,ivm,iele)=resp_cato(ipts,ivm,iele)+KC(ipts)*carbon_cato(ipts,ivm,iele)

                   ! modified for trunk xw
                   resp_hetero_soil(ipts,ivm) = (resp_acro_oxic(ipts,ivm,iele)+resp_acro_anoxic(ipts,ivm,iele)+resp_cato(ipts,ivm,iele))/dt

                   !modified xw for trunk
                   litter_to_acro(ipts,ivm,iele)=litter_to_acro(ipts,ivm,iele)+som_input(ipts,icarb,ivm,iele)*dt+ &
                        &   som_input(ipts,icarb,ivm,iele)*dt
                   carbon_acro(ipts,ivm,iele)=carbon_acro(ipts,ivm,iele)+litter_to_acro(ipts,ivm,iele)- &
                        & acro_to_cato(ipts,ivm,iele)-resp_acro_oxic(ipts,ivm,iele)-resp_acro_anoxic(ipts,ivm,iele)
                   carbon_cato(ipts,ivm,iele)=carbon_cato(ipts,ivm,iele)+acro_to_cato(ipts,ivm,iele)-resp_cato(ipts,ivm,iele)
                   
                   tcarbon_acro(ipts)=tcarbon_acro(ipts)+carbon_acro(ipts,ivm,iele)
                   tcarbon_cato(ipts)=tcarbon_cato(ipts)+carbon_cato(ipts,ivm,iele)
                   ! P_A: acrotelm density (g/m3)
                   ! P_C: catotelm density (g/m3)
                   ! cf_A: carbon fraction in acrotelm peat
                   ! cf_C: carbon fraction in catotelm peat
                   height_acro(ipts)=tcarbon_acro(ipts)/(p_A*cf_A)
                   height_cato(ipts)=tcarbon_cato(ipts)/(p_C*cf_C) 
                ENDIF
            ENDDO
          ENDDO
       ENDDO
    ENDIF
    ! peatland++

    !! 4. Check mass balance closure
    !! 4.1 For Carbon 
    !! 4.1.1 Calculate components of the mass balance
    pool_end(:,:,:) = zero 
    DO iele = 1,nelements
       DO icarb = 1,ncarb
          pool_end(:,:,iele) = pool_end(:,:,iele) + &
               ( som(:,icarb,:,iele) * veget_max(:,:) )
       ENDDO
    ENDDO

    pool_end(:,:,initrogen) = pool_end(:,:,initrogen) + &
         (n_mineralisation(:,:) * veget_max(:,:) ) 

    !! 4.1.2 Calculate mass balance
    !  Note that soilcarbon is transfered to other pools but that the pool
    !  was not updated. We should not account for it in ::pool_end
    check_intern(:,:,:,:) = zero
    check_intern(:,:,iatm2land,icarbon) = zero
    check_intern(:,:,iland2atm,icarbon) = -un * resp_hetero_soil(:,:) * veget_max(:,:) * dt 
    check_intern(:,:,ilat2out,icarbon) = -un *zero
    check_intern(:,:,ilat2in,icarbon) =  zero
    check_intern(:,:,ipoolchange,icarbon) = -un * (pool_end(:,:,icarbon) - &
         pool_start(:,:,icarbon))
    closure_intern = zero
    DO imbc = 1,nmbcomp
       closure_intern(:,:,icarbon) = closure_intern(:,:,icarbon) + &
            check_intern(:,:,imbc,icarbon)
       IF(printlev_loc.GT.4)THEN
          WRITE(numout,*) "check_intern(:,:,imbc,icarbon)",&
               check_intern(:,:,imbc,icarbon)
       ENDIF
    ENDDO

    !! 4.2.2 Calculate mass balance
    !  Note that soilcarbon is transfered to other pools but that the pool
    !  was not updated. We should not account for it in ::pool_end
    check_intern(:,:,iatm2land,initrogen) = zero
    check_intern(:,:,iland2atm,initrogen) = -un * zero
    check_intern(:,:,ilat2out,initrogen) = -un* zero
    check_intern(:,:,ilat2in,initrogen) =  zero
    check_intern(:,:,ipoolchange,initrogen) = -un * &
         (pool_end(:,:,initrogen) - pool_start(:,:,initrogen))

    DO imbc = 1,nmbcomp
       closure_intern(:,:,initrogen) = closure_intern(:,:,initrogen) + &
            check_intern(:,:,imbc,initrogen)
       IF(printlev_loc.GT.4) THEN
          WRITE(numout,*) "check_intern(:,:,imbc,initrogen)",&
               check_intern(:,:,imbc,initrogen)
       ENDIF
    ENDDO
    
    CALL check_mass_balance("som_dynamics", closure_intern, npts, pool_end, &
         pool_start, veget_max, 'pft')

    ! The above mass balance check, checks for imbalances of 10-8. If the
    ! test fails there is a serious problem. If the test succeeds we
    ! assume that the imbalance is due to precision errors. We will close
    ! the mass balance here with a high precision by redistributing the 
    ! residual.
    DO ipts = 1,npts
       DO ivm = 1,nvm
          IF (veget_max(ipts,ivm) .EQ. zero) CYCLE
          DO iele = 1,nelements

             ! Calculate the mismatch
             IF (iele.EQ.icarbon) THEN
                mismatch = SUM(som_start(ipts,:,ivm,iele)) - &
                     SUM(som(ipts,:,ivm,iele)) - resp_hetero_soil(ipts,ivm) * dt
             ELSE
                mismatch = SUM(som_start(ipts,:,ivm,iele)) - &
                     SUM(som(ipts,:,ivm,iele)) + n_mineralisation_start(ipts,ivm) - &
                     n_mineralisation(ipts,ivm)
             END IF

             ! The mismatch is likely a precision error. Fix the issue rather than 
             ! stopping the model
             IF(ABS(mismatch).GT.1000*EPSILON(un)) THEN
                
                ! Close the mass balance by adding the residual to passive pool. Note
                ! that the difference between som (10e3) and mismatch (<10e-13) may be
                ! so large that the correction has no effect due to rounding errors. It
                ! is not considered a big problem but it could frustrating when looking 
                ! into the details of the calculations at a single pixel.
                ix = MAXLOC(som(ipts,:,ivm,iele),1)
                som(ipts,ix,ivm,iele) = som(ipts,ix,ivm,iele) - mismatch

                ! Debug
                IF (printlev_loc.GE.4) THEN
                   IF (ipts.EQ.test_grid.AND.ivm.EQ.test_pft) THEN
                      WRITE(numout,*) 'Imbalance, ', iele, mismatch
                      IF (iele.EQ.icarbon) THEN
                         WRITE(numout,*) 'Adjusted imbalance C, ', &
                              EPSILON(som(ipts,ix,ivm,iele)), EPSILON(mismatch),&
                              pool_start(ipts,ivm,iele)/veget_max(ipts,ivm) - &
                              SUM(som(ipts,:,ivm,iele)) - &
                              resp_hetero_soil(ipts,ivm) * dt 
                      ELSE
                         WRITE(numout,*) 'Adjusted imbalance N, ', &
                              EPSILON(som(ipts,ix,ivm,iele)), EPSILON(mismatch), &
                              pool_start(ipts,ivm,iele)/veget_max(ipts,ivm) - &
                              SUM(som(ipts,:,ivm,iele)) - &
                              n_mineralisation(ipts,ivm)                     
                      END IF
                   END IF
                END IF
                !-
                IF (som(ipts,ix,ivm,iele).LT.zero) THEN
                   ! If this really was a precision error, it should be 
                   ! a small correction compared to the actual value of
                   ! flux. If flux becomes zero, the correction
                   ! wasn that small after all. Check.
                   WRITE(numout,*) 'ipts, ivm, ielement, icarb, ', ipts, ivm, iele, ix
                   WRITE(numout,*) 'som, ', som(ipts,ix,ivm,iele) 
                   WRITE(numout,*) 'closure_intern, ', closure_intern(ipts,ivm,iele) / &
                        veget_max(ipts,ivm)
                   CALL ipslerr_p(3,'soilcarbon','Trying to correct what was', &
                        'thought to be a precision error',&
                        'The negative values suggests there is a real problem')
                END IF
             END IF ! closure_intern
          END DO ! iele
       END DO ! ivm
    END DO ! ipts
    
 !! 5. (Quasi-)Analytical Spin-up
    
    !! 5.2 Finish to fill MatrixA with fluxes between soil pools
    
    IF (spinup_analytic) THEN

       ! During the spinup PFT1 (bare soil) should still be empty
       ! for the moment the only way to get carbon/nitrogen in PFT1
       ! is after a LCC. LCC should not be used during spinup.
       DO ivm = 2,nvm 

          ! flux leaving the active pool
          MatrixA(:,ivm,iactive_pool,iactive_pool) = moins_un * &
               dt*som_turn(iactive,ivm) * &
               control_moist(:,ibelow) * control_temp(:,ibelow) * &
               ( 1. - som_turn_iactive_clay_frac * clay(:)) * decomp_factor(ivm)
          
          ! flux received by the active pool from the slow pool
          MatrixA(:,ivm,iactive_pool,islow_pool) =  frac_carb(:,islow,iactive)*dt*som_turn(islow,ivm) * &
               control_moist(:,ibelow) * control_temp(:,ibelow) * decomp_factor(ivm)

          ! flux received by the active pool from the passive pool
          MatrixA(:,ivm,iactive_pool,ipassive_pool) =  frac_carb(:,ipassive,iactive)*dt*som_turn(ipassive,ivm) * &
               control_moist(:,ibelow) * control_temp(:,ibelow) * decomp_factor(ivm)

          ! flux leaving the slow pool
          MatrixA(:,ivm,islow_pool,islow_pool) = moins_un * &
               dt*som_turn(islow,ivm) * &
               control_moist(:,ibelow) * control_temp(:,ibelow) * decomp_factor(ivm)
          
          ! flux received by the slow pool from the active pool
          MatrixA(:,ivm,islow_pool,iactive_pool) =  frac_carb(:,iactive,islow) *&
               dt*som_turn(iactive,ivm) * &
               control_moist(:,ibelow) * control_temp(:,ibelow) * &
               ( 1. - som_turn_iactive_clay_frac * clay(:) ) * decomp_factor(ivm)
          
          ! flux received by the slow pool from the surface pool
          MatrixA(:,ivm,islow_pool,isurface_pool) =  frac_carb(:,isurface,islow) *&
               dt*som_turn(isurface,ivm) * &
               control_moist(:,ibelow) * control_temp(:,ibelow) * decomp_factor(ivm)
          
          ! flux leaving the passive pool
          MatrixA(:,ivm,ipassive_pool,ipassive_pool) =  moins_un * &
               dt*som_turn(ipassive,ivm) * &
               control_moist(:,ibelow) * control_temp(:,ibelow) * decomp_factor(ivm)      
          
          ! flux received by the passive pool from the active pool
          MatrixA(:,ivm,ipassive_pool,iactive_pool) =  &
               frac_carb(:,iactive,ipassive)* &
               dt*som_turn(iactive,ivm) * &
               control_moist(:,ibelow) * control_temp(:,ibelow) *&
               ( 1. - som_turn_iactive_clay_frac * clay(:) ) * decomp_factor(ivm)
          
          ! flux received by the passive pool from the slow pool
          MatrixA(:,ivm,ipassive_pool,islow_pool) =  &
               frac_carb(:,islow,ipassive) * &
               dt*som_turn(islow,ivm) * &
               control_moist(:,ibelow) * control_temp(:,ibelow) * decomp_factor(ivm)

          ! flux leaving the surface pool
          MatrixA(:,ivm,isurface_pool,isurface_pool) =  moins_un * &
               dt*som_turn(isurface,ivm) * &
               control_moist(:,ibelow) * control_temp(:,ibelow) * decomp_factor(ivm)

          IF (printlev_loc>=4) WRITE(numout,*)'Finish to fill MatrixA'

!!$          ! Error checking
!!$          ! If there is no veget_max, the litter pools should be empty
!!$          error_count(:) = zero
!!$          WHERE ( SUM(ABS(som(:,:,ivm,initrogen)),2) .GT. min_stomate .AND. &
!!$               veget_max(:,ivm) .LT. min_stomate)
!!$             error_count(:) = un
!!$          END WHERE
!!$          IF (SUM(error_count(:)).GT.zero) THEN
!!$             DO ipts = 1,npts
!!$                IF (error_count(ipts).GT.zero) THEN
!!$                   WRITE(numout,*) 'There is litter but no veget_max in pixel and PFT, ', ipts, ivm
!!$                   WRITE(numout,*) 'veget_max, ', veget_max(ipts,ivm)
!!$                   DO icarb = 1,ncarb
!!$                      WRITE(numout,*) 'som C, ', ipts, icarb, ivm, som(ipts,icarb,ivm,icarbon)
!!$                      WRITE(numout,*) 'som N, ', ipts, icarb, ivm, som(ipts,icarb,ivm,initrogen)
!!$                   END DO
!!$                   CALL ipslerr_p(3, 'stomate_soilcarbon','Problem when preparing for the analytical spinup',&
!!$                        'veget_max indicates that the PFT is not available','however there is soilcarbon')
!!$                END IF
!!$             END DO
!!$          END IF
!!$          !-
          
          WHERE (som(:,isurface,ivm,initrogen) .GT. min_stomate)
             CN_som_litter_longterm(:,ivm,isurface_pool) = ( CN_som_litter_longterm(:,ivm,isurface_pool) * (tau_CN_longterm-dt) &
                  + som(:,isurface,ivm,icarbon)/som(:,isurface,ivm,initrogen) * dt)/ (tau_CN_longterm)
          ENDWHERE

          WHERE (som(:,iactive,ivm,initrogen) .GT. min_stomate)
             CN_som_litter_longterm(:,ivm,iactive_pool)  = ( CN_som_litter_longterm(:,ivm,iactive_pool) * (tau_CN_longterm-dt) &
                  + som(:,iactive,ivm,icarbon)/som(:,iactive,ivm,initrogen) * dt)/ (tau_CN_longterm)
          ENDWHERE
          
          WHERE(som(:,islow,ivm,initrogen) .GT. min_stomate)
             CN_som_litter_longterm(:,ivm,islow_pool)    = ( CN_som_litter_longterm(:,ivm,islow_pool) * (tau_CN_longterm-dt) &
                  + som(:,islow,ivm,icarbon)/som(:,islow,ivm,initrogen) * dt)/ (tau_CN_longterm)
          ENDWHERE

          WHERE(som(:,ipassive,ivm,initrogen) .GT. min_stomate)
             CN_som_litter_longterm(:,ivm,ipassive_pool) =  ( CN_som_litter_longterm(:,ivm,ipassive_pool) * (tau_CN_longterm-dt) &
                  + som(:,ipassive,ivm,icarbon)/som(:,ipassive,ivm,initrogen) * dt)/ (tau_CN_longterm)
          ENDWHERE
      
       ENDDO ! Loop over # PFTS

       ! Add Identity for each submatrix(7,7) 
       DO inbpools = 1,nbpools
          MatrixA(:,:,inbpools,inbpools) = MatrixA(:,:,inbpools,inbpools) + un 
       END DO
    ENDIF ! (spinup_analytic)

    ! Output diagnostics
    DO icarb = 1, ncarb ! Loop over carbon pools
       DO ipts = 1, npts
          IF (SUM(decomp_rate_soilcarbon(ipts,:,icarb)*veget_max(ipts,:)) > min_sechiba) THEN
             tsoilpools(ipts,icarb) = 1./(SUM(decomp_rate_soilcarbon(ipts,:,icarb)*veget_max(ipts,:))/dt_sechiba)
          ELSE
             tsoilpools(ipts,icarb) = xios_default_val
          END IF
       END DO
    END DO
    CALL xios_orchidee_send_field("tSoilPools",tsoilpools)


    IF (printlev_loc>=4) WRITE(numout,*) 'Leaving som_dynamics'
 
  END SUBROUTINE som_dynamics


!! ================================================================================================================================
!!  SUBROUTINE   : nitrogen_dynamics_clear
!!
!>\BRIEF        Set the flag ::firstcall to .TRUE. 
!! 
!! 
!_ ================================================================================================================================
  
  SUBROUTINE nitrogen_dynamics_clear
    firstcall_nitrogen=.TRUE.
  END SUBROUTINE nitrogen_dynamics_clear



!! ================================================================================================================================
!!  SUBROUTINE   : nitrogen_dynamics
!!
!>\BRIEF        : Describes mineralization dynamics of nitrogen species in the
!!  soil.  Inspired by the DNDC model of Li et al (1992,2000), but very simplified 
!!  to avoid having to calculate microbe growth for the moment. Builds on the 
!!  physico-chemical reactions with some fixed assumptions.
!!
!! DESCRIPTION	: Plants can only uptake nitrogen in the form of NH4 and NO3.
!!   In the soil, nitrogen coverts between NH3, NH4, NO3, NOX, N2O, and N2,
!!   depending on things like the concentrations of each species, the pH of
!!   the soil, the temperature, the amount of oxygen present, and the soil
!!   moisture content.  This subroutine describes the dynamics between
!!   all these species, influencing the amount of nitrogen available for
!!   both plant uptake and litter decomposition.
!!
!!   Ammonium (NH4) may be the most important species, as all litter
!!   N is assumed to decompose to NH4.  NH4 can then be lost via
!!   adsorpotion to soil clays, transformed to NH3 or NO3, or uptaken
!!   by plants.
!!
!! The subroutine performs the following tasks:\n
!!
!! Section 2.\n
!! Immobilizes nitrogen, reducing the amount of nitrogen available in first
!! the ammonium (NH4) pool and then the nitrate (NO3) pool.  This amount is
!! later added back into the mineralized soil NH4 pool.  In essence, it
!! saves this quantity to be used in the decomposition routines.  Without
!! this step, all of the nitrogen could be used up in this subroutine, thus
!! preventing litter decomposition from happening. \n
!!
!! Section 3.\n
!! The goal of this section is to calculate the volumetric fraction of
!! aneorobic microsites in the soil, which gives an idea of how much
!! aneorobic bacteria can be transforming soil nitrogen.  It depends
!! on the amount of oxygen in the soil, which depends on the concentration
!! difference of oxygen between the atmosphere and the soil, the diffusion
!! rate of oxygen in the soil, and how much oxygen is consumed by 
!! respiration in the soil.  For ideal gases, "concentration" is given by
!! the partial pressure.\n
!!
!! Section 4.\n
!! This section takes into account the losses of ammonium in the soil
!! due to adsorption onto clays, as well as the leaching of ammonium
!! and nitrate out of the system.\n
!!
!! Section 5.\n
!! Nitrification, that is, the change from ammonium (NH4+) to nitrate,
!! including effects of soil temperature, moisture, and pH, and
!! losses from NH4+ conversion to N2O and NO. \n
!!
!! Section 6.\n
!! Denitrification, that is the process of nitrate (NO3-) losing
!! oxygen atoms to eventually form N2, via the intermediate steps 
!! of NO2-, NO, and N2O.  All of this happens in the soil, and
!! depends on bacterial biomass, soil temperature, pH, and
!! soil water content.\n
!!
!! Section 7.\n
!! Calculates how much ammonium and nitrate is taken up by the plants.\n
!!
!! Section 8.\n
!! Calculates how much nitrogen is lost through emission into the
!! atmosphere through the surface of the soil.  Depends on the rate of
!! diffusion through the soil.  Nitrate emission is set to zero, but
!! the rest can pass from the soil to the atmosphere.  Emission is given
!! the lowest priority of any loss pathway, in the sense that if the 
!! amount of the species after changes due to leeching, nitrification, and
!! denitrification is less than what is calculated for the emissions, the
!! amount of emission is reduced to match.\n
!!
!! Section 9.\n
!! Update the values of the soil nitrogen pools, taking into account
!! all the losses and gains computed above.  After this, add any
!! external nitrogen inputs, like fertilizer and biological nitrogen
!! fixation (BNF).\n
!!
!! Section 10.\n
!! Write output variables.\n
!!
!! RECENT CHANGE(S): 
!! 
!! MAIN OUTPUTS VARIABLE(S): soil_n_min, mineralisation
!!
!! REFERENCE(S)   :
!! - S. Zaehle and A. D. Friend (2010), Carbon and nitrogen cycle dynamics in the
!!   O-CN land surface model: 1. Model description, site-scale evaluation, and 
!!   sensitivity to parameter estimates. Global Biogeochem. Cycles, 24, GB1005, 
!!   doi:10.1029/2009GB003521.
!! - Li C. S., J. Aber, F. Stange, K. Butterbach-Bahl, and H. Papen (2000), 
!!   A process-oriented model of N2O and NO emissions from forest soils: 
!!   1. Model development, Journal of Geophysical Research-Atmospheres, 
!!   105, 4369-4384.
!! - Li, C., S. Frolking, and T. A. Frolking (1992), A model of nitrous 
!!   oxide evolution from soil driven by rainfall events: 1. Model 
!!   structure and sensitivity, J. Geophys. Res., 97(D9), 9759–9776, 
!!   doi:10.1029/92JD00509. 
!! - Saxton, K.E., Rawls, W.J., Romberger, J.S., Papendick, R.I., 1986
!!   Estimating generalized soil-water characteristics from texture. 
!!   Soil Sci. Soc. Am. J. 50, 1031-1036
!! - Kesik, M., Ambus, P., Baritz, R., Brüggemann, N., et al.
!!   Inventories of N2O and NO emissions from European forest soils, 
!!   Biogeosciences, 2, 353-375, https://doi.org/10.5194/bg-2-353-2005, 2005. 
!! - Schmid, M., Neftel, A., Riedo, M. et al. 
!!   Process-based modelling of nitrous oxide emissions from different nitrogen sources in mown grassland
!!   Nutrient Cycling in Agroecosystems 60: 177. https://doi.org/10.1023/A:1012694218748, 2001.
!! - Garnett, T., Plett, D., Conn, V., Conn, S., Rabie, H., Rafalski, J. A., ... & Kaiser, 
!!   B. N. (2015). Variation for N uptake system in maize: genotypic response to N supply. 
!!   Frontiers in plant science, 6, 936.
!! - Kulmatiski, A., Adler, P. B., Stark, J. M., & Tredennick, A. T. (2017). Water and nitrogen 
!!   uptake are better associated with resource availability than root biomass. Ecosphere, 8(3), e01738.
!! - Min X., M. Y. Siddiqi, R. D. Guy, A. D. M. Glass, and H. J. Kronzucker (1998), Induction of 
!!   nitrate uptake and nitrate reductase activity in trembling aspen and lodgepole pine, Plant, 
!!  Cell and Environment, 21, 1039-1046. 
!! 
!! FLOWCHART    :
!!
!_ ================================================================================================================================
  SUBROUTINE nitrogen_dynamics(npts, njsc, clay, sand, &
       tsoil_decomp, tmc_pft, drainage_pft, runoff_pft, swc_pft, veget_max, resp_sol, &
       som, n_input, month, pH, &
       n_mineralisation, pb, plant_n_uptake, Bd, soil_n_min, &
       p_O2, bact, atm_to_immob, leaching, emission, ld_redistribute, &
       circ_class_biomass, circ_class_n, cn_leaf_min_2D, cn_leaf_max_2D, cn_leaf_init_2D, &
       mcs_hydrol, mcfc_hydrol, croot_longterm, n_reserve_longterm, sugar_load)
    
    ! 0 declarations
    
    ! 0.1 input
    INTEGER(i_std), INTENT(in)                             :: npts            !! Domain size
    INTEGER(i_std), DIMENSION(:), INTENT(in)               :: njsc            !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless)
    REAL(r_std), DIMENSION(:), INTENT(in)                  :: clay            !! clay fraction (between 0 and 1)
    REAL(r_std), DIMENSION(:), INTENT(in)                  :: sand            !! sand fraction (between 0 and 1)
    REAL(r_std), DIMENSION(:), INTENT(in)                  :: tsoil_decomp    !! soil temperature (degC)
    REAL(r_std), DIMENSION(:,:), INTENT(in)                :: tmc_pft         !! Total soil water per PFT (mm/m2)
    REAL(r_std), DIMENSION(:,:), INTENT(in)                :: drainage_pft    !! Drainage per PFT (mm/m2)
    REAL(r_std), DIMENSION(:,:), INTENT(in)                :: runoff_pft      !! Runoff per PFT (mm/m2 /dt_sechiba) 
    REAL(r_std), DIMENSION(:,:), INTENT(in)                :: swc_pft         !! Relative Soil water content [tmcr:tmcs] per pft (-)
    REAL(r_std), DIMENSION(:,:), INTENT(in)                :: veget_max       !! fraction of a vegetation (0-1)
    REAL(r_std), DIMENSION(:,:), INTENT(in)                :: resp_sol        !! carbon respired from below ground 
                                                                              !! (hetero+autotrophic) (gC/m**2/day)
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)            :: som             !! SOM (gC(or N)/m**2)
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)            :: n_input         !! nitrogen inputs into the soil  (gN/m**2/day)
                                                                              !! NH4 and NOX from the atmosphere, NH4 from BNF,
                                                                              !! agricultural fertiliser as NH4/NO3
    INTEGER(i_std), INTENT(in)                             :: month           !! month number required for n_input (1-12)
    REAL(r_std), DIMENSION(:), INTENT(in)                  :: pH              !! soil pH
    REAL(r_std), DIMENSION(:,:,:,:,:),INTENT(in)           :: circ_class_biomass !! Biomass components of the model tree  
                                                                              !! within a circumference class
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)              :: circ_class_n    !! Number of trees in each circumference class
    REAL(r_std), DIMENSION(:,:), INTENT(inout)             :: n_mineralisation!! net nitrogen mineralisation of decomposing SOM
                                                                              !! (gN/m**2/day), supposed to be NH4
    REAL(r_std), DIMENSION(:), INTENT(in)                  :: pb              !! Air pressure (hPa)
    REAL(r_std), DIMENSION(:), INTENT(in)                  :: Bd              !! Bulk density (kg/m**3)
    REAL(r_std),DIMENSION(:,:), INTENT(in)                 :: cn_leaf_min_2D  !! minimal leaf C/N ratio 
    REAL(r_std),DIMENSION(:,:), INTENT(in)                 :: cn_leaf_max_2D  !! maximal leaf C/N ratio 
    REAL(r_std),DIMENSION(:,:), INTENT(in)                 :: cn_leaf_init_2D !! initial leaf C/N ratio
    REAL(r_std),DIMENSION(:,:), INTENT(in)                 :: croot_longterm  !! "Long term" (default 3 years) root carbon mass
                                                                              !! per ground area
                                                                              !! @tex $(gC m^{-2} year^{-1})$ @endtex
    REAL(r_std),DIMENSION(:), INTENT(in)                   :: mcs_hydrol      !! Saturated volumetric water content output to be used in stomate_soilcarbon
    REAL(r_std),DIMENSION(:), INTENT(in)                   :: mcfc_hydrol     !! Volumetric water content at field capacity output to be used in stomate_soilcarbon
    REAL(r_std), DIMENSION(:), INTENT(in)                  :: ld_redistribute !! logical set to redistribute som and litter
    REAL(r_std), DIMENSION(:,:), INTENT(in)                :: n_reserve_longterm !! 'longterm' balance of actual N reserve to 
                                                                              !! potential N reserve
    REAL(r_std), DIMENSION(:,:), INTENT(in)                :: sugar_load      !! Relative sugar loading of the labile pool (unitless)
    REAL(r_std), DIMENSION(:,:,:),INTENT(inout)            :: soil_n_min      !! mineral nitrogen in the soil (gN/m**2) 
    REAL(r_std), DIMENSION(:,:),INTENT(inout)              :: p_O2            !! partial pressure of oxygen in the soil (hPa)
    REAL(r_std), DIMENSION(:,:),INTENT(inout)              :: bact            !! denitrifier biomass (gC/m**2)

    ! 0.2 output
    REAL(r_std), DIMENSION(:,:,:), INTENT(out)             :: plant_n_uptake  !! Uptake of soil N by plants 
                                                                              !! (gN/m**2/timestep)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)             :: atm_to_immob    !! Nitrogen taken from the atmosphere to 
                                                                              !! support immobilisation (gN m-2 dt-1)

    ! 0.3 local
    REAL(r_std)                                                        :: dt              !! Time step \f$(dt_sechiba one_day^{-1})$\f
    LOGICAL                                                            :: l_error         !! Diagnostic boolean for error allocation (true/false)
    INTEGER(i_std)                                                     :: ier             !! Check errors in netcdf call
    REAL(r_std), DIMENSION(npts,nvm,nionspec)                          :: leaching        !! mineral nitrogen leached from the soil 
                                                                                          !! (gN/m**2/timestep)
    REAL(r_std), DIMENSION(npts,nvm)                                   :: immob           !! N immobilized (gN/m**2/day)
    REAL(r_std), DIMENSION(npts)                                       :: afps_max        !! maximum pore-volume of the soil
                                                                                          !! Table 2, Li et al., 2000
                                                                                          !! (fraction)                       
    REAL(r_std), DIMENSION(npts,nvm)                                   :: afps            !! pore-volume of the soil
                                                                                          !! Table 2, Li et al., 2000
                                                                                          !! (fraction)   
    REAL(r_std), DIMENSION(npts,nvm)                                   :: D_s             !! Oxygen diffusion in soil (m**2/day)
    REAL(r_std), DIMENSION(npts,nvm)                                   :: mol_O2_resp     !! Density of Moles of O2 related to respiration term 
                                                                                          !! ((molesO2 m-3) (gC m-2)-1)
    REAL(r_std), DIMENSION(npts,nvm)                                   :: p_O2_resp       !! O2 partial pressure related to the respiration term 
                                                                                          !! ((hPa) (gC m-2)-1)
    REAL(r_std), DIMENSION(npts)                                       :: p_O2air         !! Oxygen partial pressure in air (hPa)
    REAL(r_std), DIMENSION(npts)                                       :: g_O2            !! gradient of O2 partial pressure (hPa m-1)
    REAL(r_std), DIMENSION(npts)                                       :: d_O2            !! change in O2 partial pressure (hPa day-1)
    REAL(r_std), DIMENSION(npts,nvm)                                   :: anvf            !! Volumetric fraction of anaerobic microsites 
                                                                                          !! (fraction of pore-volume)
    REAL(r_std), DIMENSION(npts)                                       :: FixNH4          !! Fraction of adsorbed NH4+ (-)
    REAL(r_std), DIMENSION(npts,nvm)                                   :: n_adsorbed      !! Ammonium adsorpted (gN/m**2)
                                                                                          !! based on Li et al. 1992, JGR, Table 4
    REAL(r_std), DIMENSION(npts,nvm)                                   :: fwnit           !! Effect of soil moisture on nitrification (-)
                                                                                          !! Zhang et al. 2002, Ecological Modelling, appendix A, page 101
    REAL(r_std), DIMENSION(npts,nvm)                                   :: fwdenit         !! Effect of soil moisture on nitrification (-)
                                                                                          !! Zhang et al. 2002, Ecological Modelling, appendix A, page 101
    REAL(r_std), DIMENSION(npts)                                       :: var_temp_sol    !! Temperature function used for calc. ft_nit (-)
                                                                                          !! Zhang et al. 2002, Ecological Modelling, appendix A, page 101
    REAL(r_std), DIMENSION(npts)                                       :: ft_nit          !! Effect of temperature on nitrification (-)
                                                                                          !! Zhang et al. 2002, Ecological Modelling, appendix A, page 101
    REAL(r_std), DIMENSION(npts)                                       :: fph             !! Effect of pH on nitrification (-)
                                                                                          !! Zhang et al. 2002, Ecological Modelling, appendix A, page 101
    REAL(r_std), DIMENSION(npts)                                       :: ftv             !! Effect of temperature on NO2 or NO production 
                                                                                          !! during nitrification (-)
                                                                                          !! Zhang et al. 2002, Ecological Modelling, appendix A, page 102
    REAL(r_std), DIMENSION(npts,nvm,3)                                 :: nitrification   !! N-compounds production (NO3, N2O, NO) related 
                                                                                          !! to nitrification process (gN/m**2/tstep)
    REAL(r_std), DIMENSION(npts)                                       :: ft_denit        !! Temperature response of relative growth rate of 
                                                                                          !! total denitrifiers (-)
                                                                                          !! Eq. 2 Table 4 of Li et al., 2000
    REAL(r_std), DIMENSION(npts)                                       :: fph_no3         !! Soil pH response of relative growth rate of 
                                                                                          !! NO3 denitrifiers (-)
                                                                                          !! Eq. 2 Table 4 of Li et al., 2000
    REAL(r_std), DIMENSION(npts)                                       :: fph_no          !! Soil pH response of relative growth rate of 
                                                                                          !! NO denitrifiers (-)
                                                                                          !! Eq. 2 Table 4 of Li et al., 2000
    REAL(r_std), DIMENSION(npts)                                       :: fph_n2o         !! Soil pH response of relative growth rate of 
                                                                                          !! N2O denitrifiers (-)
                                                                                          !! Eq. 2 Table 4 of Li et al., 2000
    REAL(r_std), DIMENSION(npts,nvm)                                   :: Kn_conv         !! Conversion from kgN/m3 to gN/m2
    REAL(r_std), DIMENSION(npts)                                       :: mu_NO3          !! Relative growth rate of NO3 denitrifiers (hour**-1)
                                                                                          !! Eq.1 Table 4 Li et al., 2000 
    REAL(r_std), DIMENSION(npts)                                       :: mu_N2O          !! Relative growth rate of N2O denitrifiers (hour**-1)
                                                                                          !! Eq.1 Table 4 Li et al., 2000
    REAL(r_std), DIMENSION(npts)                                       :: mu_NO           !! Relative growth rate of NO denitrifiers (hour**-1)
                                                                                          !! Eq.1 Table 4 Li et al., 2000  
    REAL(r_std), DIMENSION(npts)                                       :: sum_n           !! sum of all N species in the soil (gN/m**2)
    REAL(r_std), DIMENSION(npts,nvm,3)                                 :: denitrification !! N-compounds consumption (NO3, N2O, NO) related 
                                                                                          !! to denitrificaiton process (gN/m**2/tstep)
    REAL(r_std), DIMENSION(npts)                                       :: dn_bact         !! denitrifier biomass change (kgC/m**3/timestep)
    REAL(r_std), DIMENSION(npts)                                       :: ft_uptake       !! Temperature response of N uptake by plants (-)
    REAL(r_std)                                                        :: conv_fac_vmax   !! Conversion from (umol (gDW)-1 h-1) to (gN (gC)-1 timestep-1)
    REAL(r_std), DIMENSION(npts,nvm)                                   :: nc_leaf_min     !! Minimal NC ratio of leaf (gN / gC)
    REAL(r_std), DIMENSION(npts,nvm)                                   :: nc_leaf_max     !! Maximal NC ratio of leaf (gN / gC)
    REAL(r_std), DIMENSION(npts)                                       :: lab_n           !! Labile nitrogen in plants (gN/m**2) 
    REAL(r_std), DIMENSION(npts)                                       :: lab_c           !! Labile carbon in plants (gC/m**2) 
    REAL(r_std), DIMENSION(npts)                                       :: NCplant         !! NC ratio of the plant (gN / gC)
                                                                                          !! Eq. (9) p. 3 of SM of Zaehle & Friend, 2010
    REAL(r_std), DIMENSION(npts,nvm)                                   :: f_NCplant       !! Response of Nitrogen uptake by plants 
                                                                                          !! to N/C ratio of the labile pool
    REAL(r_std)                                                        :: conv_fac_concent!! Conversion factor from (umol per litter) to (gN m-2) 
    REAL(r_std), DIMENSION(npts)                                       :: frac_nh3        !! dissociation of [NH3] to [NH4+] (-)
    REAL(r_std), DIMENSION(npts)                                       :: F_clay          !! Response of N-emissions to clay fraction (-)
    REAL(r_std), DIMENSION(npts,nvm,nnspec)                            :: emission        !! volatile losses of nitrogen 
                                                                                          !! (gN/m**2/timestep)
    INTEGER(i_std)                                                     :: ipts,ivm,inspec,ininput,iion,imbc,ipar,iele !! Index 
    REAL(r_std), DIMENSION(npts,nvm)                                   :: f_drain         !! Fraction of tmc which has been drained 
                                                                                          !! in the last time step
    REAL(r_std), DIMENSION(npts,nvm,nmbcomp,nelements)                 :: check_intern    !! Contains the components of the internal
                                                                                          !! mass balance chech for this routine
                                                                                          !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)                         :: closure_intern  !! Check closure of internal mass balance
                                                                                          !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)                         :: pool_start      !! Start and end pool of this routine 
                                                                                          !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)                         :: pool_end        !! Start and end pool of this routine 
                                                                                          !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)                                   :: err             !! Array to store errors (dimensionless)
    REAL(r_std), DIMENSION(npts)                                       :: temp_sol        !! soil temperature (C) 
    REAL(r_std), DIMENSION(npts)                                       :: active_biomass  !! temporary variable of active root biomass (gC m-2)
 !_ ================================================================================================================================

    !! Initialize printlev_loc only first time one of the subroutines 
    IF ( firstcall_soilcarbon ) THEN
      ! Guillaume M. -- get_printlev calls getin_p, a COLLECTIVE bcast. Called here,
      ! inside the compute loop under a local firstcall, it hangs the run as soon as the
      ! ranks stop reaching it together. Use printlev, read once at startup by every rank.
       printlev_loc=printlev
       firstcall_soilcarbon = .FALSE.
    ENDIF

    IF (printlev_loc>=3) WRITE(numout,*) 'Entering nitrogen_dynamics' 

    !! 1. Initializations
    dt = dt_sechiba/one_day

    IF ( firstcall_nitrogen ) THEN                      
       firstcall_nitrogen = .FALSE. 
    ENDIF     
   
    ! Transform tsoil_decomp into degree C 
    temp_sol(:) = tsoil_decomp(:) - tp_00 
   
    ! Debug
    IF(printlev_loc.GE.4)THEN
       WRITE(numout,*) 'dt, ',dt
       WRITE(numout,*) 'test whether we have negative soil_n_min at the start'
       DO ivm = 1,nvm
          WRITE(numout,*) "n_input_dt, ",n_input(test_grid,ivm,month,:) * dt
          IF (n_mineralisation(test_grid,ivm).LT.0) THEN
             WRITE(numout,*)'Negative mineralisation,',ivm,&
                  n_mineralisation(test_grid,ivm)
          ENDIF
          IF (soil_n_min(test_grid,ivm,initrate).LT.0 .OR. &
               soil_n_min(test_grid,ivm,iammonium).LT.0) THEN
             WRITE(numout,*)'Hit exception', dt
             WRITE(numout,*)'PFT, ', ivm 
             WRITE(numout,*)'soil_n_min_initrate:', &
                  soil_n_min(test_grid,ivm,initrate)
             WRITE(numout,*)'soil_n_min_initrate:', &
                  soil_n_min(test_grid,ivm,iammonium)   
          ENDIF
       ENDDO
    ENDIF
    !-   

    !! 1.2 Initialize check for mass balance closure 
    IF (err_act.GT.1) THEN

       check_intern(:,:,:,:) = zero
       pool_start(:,:,:) = zero

       ! atm_to_immob has as intent inout, uptake may happen to
       ! enable immobilisation
       check_intern(:,:,iatm2land,initrogen) = - un * &
            atm_to_immob(:,:) * veget_max(:,:) * dt

       DO inspec = 1,nnspec
          pool_start(:,:,initrogen) = pool_start(:,:,initrogen) + &
               soil_n_min(:,:,inspec) * veget_max(:,:)
       ENDDO    

    ENDIF
    
    !! 2. Preparation for decomposition
    ! 
    ! 2.1 conservation of mass from decomposition
    ! immobilisation has absolute priority to avoid mass conservation problems
    ! the code in litter and soilcarbon has to make sure that soil ammonium can never be 
    ! more than exhausted completely by immobilisation!!!
    ! As a reminder, mineralisation is the amount of litter destined to become soil NH4
    immob(:,:) = zero
    leaching(:,:,:) = zero
    f_drain(:,:) = zero 
    emission(:,:,:) = zero
    plant_n_uptake(:,:,:) = zero
    nitrification(:,:,:) = zero
    denitrification(:,:,:) = zero
    n_adsorbed(:,:) = zero
    f_NCplant(:,:) = zero

    ! 1.4 Conservation of mass from decomposition
    ! Immobilisation has absolute priority. To avoid mass 
    ! conservation problems the code in litter and soilcarbon 
    ! has to make sure that soil ammonium can never be 
    ! more than exhausted completely by immobilisation!
    !
    ! We could try to correct the issues when it happens in 
    ! stomate_litter or we could try to fix it here. Because
    ! less pools are accounted for here, it seems more easy
    ! to deal with the negative mineralization here.
    DO ipts = 1,npts

       DO ivm = 1,nvm 

          ! Choose between two approaches to deal with mass balance
          ! problems: close the mass balance at the end of the calculations
          ! or close it in every subroutine seperatly.
          IF (.NOT.test_adjust_rhetero) THEN
             
             ! There is immobilisation
             IF(veget_max(ipts,ivm).GT.min_stomate .AND. & 
                  n_mineralisation(ipts,ivm).LT.zero) THEN

                ! Move mineralisation to immobilisation to have
                ! a temporary variable that will be used to update
                ! the soil pools
                immob(ipts,ivm) = - n_mineralisation(ipts,ivm)

                ! ORC3 approach: take the missing N from the atmosphere.
                soil_n_min(ipts,ivm,initrate) = soil_n_min(ipts,ivm,initrate) - &
                     MAX(zero,(immob(ipts,ivm)-soil_n_min(ipts,ivm,iammonium)))
                soil_n_min(ipts,ivm,iammonium) = soil_n_min(ipts,ivm,iammonium) - &
                     MIN(immob(ipts,ivm),soil_n_min(ipts,ivm,iammonium))
                
                IF (soil_n_min(ipts,ivm,initrate) .LT. 0.) THEN

                   ! Take the missing nitrogen from the atmosphere. Leave a little
                   ! bit of mineral nitrate in the soil
                   atm_to_immob(ipts,ivm) = atm_to_immob(ipts,ivm) - &
                        soil_n_min(ipts,ivm,initrate)/dt + 2 * min_n/dt
                   soil_n_min(ipts,ivm,initrate) = min_n
                   soil_n_min(ipts,ivm,iammonium) = min_n
                   
                END IF

             ELSEIF(veget_max(ipts,ivm).GT.min_stomate .AND. & 
                  n_mineralisation(ipts,ivm).GE.zero) THEN
                
                ! Add the mineral nitrogen from litter decomposition
                soil_n_min(ipts,ivm,iammonium) = soil_n_min(ipts,ivm,iammonium) + &
                     n_mineralisation(ipts,ivm)
                
             ELSEIF(veget_max(ipts,ivm).LE.min_stomate) THEN
                
                soil_n_min(ipts,ivm,initrate) = zero
                soil_n_min(ipts,ivm,iammonium) = zero
                
             ELSE
                
                WRITE(numout,*)'at ipts veget_max for ivm: ',ipts, ivm, veget_max(ipts,ivm)
                WRITE(numout,*)'n_mineralisation', n_mineralisation(ipts,ivm)                     
                CALL ipslerr_p (3,'nitrogen_dynamics', 'Unforeseen case. for mineralization','','')
                
             ENDIF ! n_mineralisation .LT. zero
                
          ELSE

             ! There is immobilisation
             IF(veget_max(ipts,ivm).GT.min_stomate .AND. & 
                  n_mineralisation(ipts,ivm).LT.zero) THEN

                ! Move mineralisation to immobilisation to have
                ! a temporary variable that will be used to update
                ! the soil pools
                immob(ipts,ivm) = - n_mineralisation(ipts,ivm)

                ! There is quite some code in stomate_litter.f90 and
                ! som_dynamics (in this module) to ensure that there is
                ! enough soil nitrogen to sustain the immobilisation
                ! Test whether the previous correction were succesful. 
                ! If not, we will to go back to the aforementioned subroutines
                ! and fix the problems there.
                IF (immob(ipts,ivm)-(soil_n_min(ipts,ivm,initrate)+ &
                     soil_n_min(ipts,ivm,iammonium)).LT.min_stomate)THEN

                   ! In case the N related to immobilisation is higher that the N
                   ! available in the [NH4+] pool, we take the remaining N from the 
                   ! [NO3-] pool. This code is pretty condensed and replaces a set of 
                   ! IF-statements.
                   soil_n_min(ipts,ivm,initrate) = soil_n_min(ipts,ivm,initrate) - &
                        MAX(zero,(immob(ipts,ivm)-soil_n_min(ipts,ivm,iammonium)))
                   soil_n_min(ipts,ivm,iammonium) = soil_n_min(ipts,ivm,iammonium) - &
                        MIN(immob(ipts,ivm),soil_n_min(ipts,ivm,iammonium))

                ELSE

                   WRITE(numout,*) 'ERROR: immob exceeds the available soil nitrogen &
                        & pool',ipts,ivm,immob(ipts,ivm),(soil_n_min(ipts,ivm,initrate)+&
                        soil_n_min(ipts,ivm,iammonium))
                   CALL ipslerr_p(3,'nitrogen_dynamics',&
                        'immob exceeds the available soil nitrogen',&
                        'this problem should be fixed in stomate_litter',&
                        'and/or in som_dynamics')
                END IF

             ELSEIF(veget_max(ipts,ivm).GT.min_stomate .AND. & 
                  n_mineralisation(ipts,ivm).GE.zero) THEN

                ! Add the mineral nitrogen from litter decomposition
                soil_n_min(ipts,ivm,iammonium) = soil_n_min(ipts,ivm,iammonium) + &
                     n_mineralisation(ipts,ivm)

             ELSEIF(veget_max(ipts,ivm).LE.min_stomate) THEN

                soil_n_min(ipts,ivm,initrate) = zero
                soil_n_min(ipts,ivm,iammonium) = zero
                
             ELSE

                WRITE(numout,*)'at ipts veget_max for ivm: ',ipts, ivm, veget_max(ipts,ivm)
                WRITE(numout,*)'n_mineralisation', n_mineralisation(ipts,ivm)                     
                CALL ipslerr_p (3,'nitrogen_dynamics', 'Unforeseen case. for mineralization','','')

             ENDIF ! n_mineralisation .LT. zero

          END IF ! test_adjust_rhetero

       ENDDO ! ivm
       
    ENDDO ! ipts

    ! Debug
    err(:,:) = zero
    WHERE(soil_n_min(:,:,initrate) .LT. -min_stomate)
       err(:,:) = 1.
    ENDWHERE

    ! Check nummerical consitency 
    WHERE(soil_n_min(:,:,iammonium) .LT. -min_stomate)
       err(:,:) = err(:,:) + 1
    ENDWHERE

    ! Check whether errors occured    
    IF (SUM(SUM(err(:,:),2)).GT.zero) THEN
       WRITE(numout,*) 'Number of errors found: ',SUM(SUM(err(:,:),2))
       WRITE(numout,*) 'Negative soil_n_min_pool, initrate',soil_n_min(:,:,initrate)
       WRITE(numout,*) 'Negative soil_n_min_pool, iammonium',soil_n_min(:,:,iammonium)
       WRITE(numout,*) 'That does not make any sense'
       WRITE(numout,*) 'Have a closer look at nitrogen_dynamics'
       CALL ipslerr_p(3,'Negative soil_n_min in nitrogen_dynamics',&
           'Compensated for by taking N from the atmosphere',&
           'This should not happen','')
    END IF

  
    !! 3. ANVF 
    !
    ! The goal of this section is to calculate the volumetric fraction of
    ! aneorobic microsites (ANVF) in the soil, which gives an idea of how much
    ! aneorobic bacteria can be transforming soil nitrogen.  It depends
    ! on the oxygen diffusion and partial pressure in the soil.  These euqations
    ! are taken from Table 2 of Li et al. (2000) Table 2, though some things 
    ! remain unclear or undefinied.  afps, for example, is never explicitly given.
    ! S. Zaehle defined it with the following equation - but did not use it.
    ! In OCN, diffusion is not of afps/afps_max ratio but accounts for soilhum and 
    ! soil temperature according to Monteith & Unsworth, 1990 
    ! d_ox(:) = 1.73664 * ( 0.15 * (exp(-(soilhum_av(:)**3.)/0.44)-exp(-1./0.44))) * & 
    !          (1.+0.007*tsoil_av(:))
  
    ! The equation for afpsmax, the maximum pore volume of the soil, 
    ! is taken from Table 2 of Saxton (1986)
    afps_max(:) = h_saxton + j_saxton * (sand(:)*100.) + &
         k_saxton * log10(clay(:)*100.)
  
    DO ivm = 1,nvm
       ! pore volume of the soil.  Unclear where this comes from.
       afps(:,ivm) = MAX(0.1,( un - swc_pft(:,ivm) ) * afps_max(:))

       ! diffusion through the soil. Table 2 of Li et al. (2000)
       D_s(:,ivm) = D_air * ( afps(:,ivm)**diffusionO2_power_1 ) / &
         ( afps_max(:)**diffusionO2_power_2 )
  
       ! Account for the impact of frost on diffusion. Table 2 of Li et al. (2000)
       WHERE ( temp_sol(:) .GT. zero )
          D_s(:,ivm) = D_s(:,ivm) * F_nofrost
       ELSEWHERE
          D_s(:,ivm) = D_s(:,ivm) * F_frost
       ENDWHERE
  
       ! Equation (3) - Oxygen partial pressure
       ! Written in an odd differential form in Table 2 of Li et al (2000), with the
       ! time derivative of the partial pressure being related to the depth derivative
       ! of both the partial pressure and the diffusion rate.
   
       ! So we do it this way, calculating the depth gradient afterwards.
       ! Oxygen loss from respiration has to be expressed as a partial pressure
       ! using the ideal gas PV=nRRT relationship with P in Pa, V in m-3 and T in K
       ! RR is the ideal gas constat (J mol-1 K-1) and n the number of moles
       ! Respiration, one mole of C requires two moles of oxygen (CO2)
       ! Resp_below expressed in gC m-2 day-1
       ! mol_O2_resp : Density of Moles of O2 related to respiration 
       ! term (molesO2 m-3) (gC m-2)-1
       ! In OCN, afps_max is used instead of afps. We keep here the 
       ! original formulation
       mol_O2_resp(:,ivm) = (un / C_molar_mass * 2. ) / (zmaxh * afps(:,ivm))

       ! O2 partial pressure related to the respiration term (hPa) (gC m-2)-1
       p_O2_resp(:,ivm)   = mol_O2_resp(:,ivm) * RR * ( temp_sol(:) + tp_00 ) * Pa_to_hPa
    ENDDO
  
    ! Change in oxygen partial pressure d_O2 - Ds x gradient
    ! Not sure that we should use z_decomp (could be a fraction of zmaxh, too)
    ! oxygen partial pressure in air - p_O2air (hPa)
    p_O2air(:)= V_O2 * pb(:)


    DO ivm = 1, nvm

       ! assume the partial pressure of oxygen in the soil is uniform over the
       ! whole depth that soil decomposers are active, in order to calculate
       ! the change in partial pressure of oxygen in the soil due to the
       ! pressure different between the air and soil and the diffusion rate
       ! of oxygen through the soil
       g_O2(:) = ( p_O2air(:) - p_O2(:,ivm) ) / z_decomp
       d_O2(:) = D_s(:,ivm) / z_decomp * g_O2(:)
       
       ! D_s / z_decomp has the unit of a conductivity (m day-1)
       ! much oxgen has entered or left the soil due to the pressure difference
       ! between the atmosphere and the soil, in addition to any losses from
       ! soil respiration consuming oxygen to create CO2
       p_O2(:,ivm) = p_O2(:,ivm) + d_O2(:)*dt - p_O2_resp(:,ivm)*resp_sol(:,ivm)*dt

       ! Equation (4) Volumetric fraction of anaerobic microsites (ANVF)
       ! a and b constants are not specified in Li et al., 2000
       ! S. Zaehle used a=0.85 and b=1 without mention to any publication
       ! a_anvf=0.85
       ! b_anvf=1.
       anvf(:,ivm) = a_anvf * ( 1 - b_anvf * p_O2(:,ivm) / p_O2air(:) )
       anvf(:,ivm) = MAX(zero, anvf(:,ivm))

    ENDDO
    
    ! To test the restartability of the model and as a general
    ! QC of the restart files, p_O2 should be zero if the PFT
    ! is not defined in the simulation.
    WHERE (veget_max(:,:).LT.min_stomate)
       p_O2(:,:) = zero
    END WHERE

    ! Debug
    IF(printlev_loc>=4 .AND. ivm==test_pft)THEN
       WRITE(numout,*) 'CHECK values after nitrification nh4 to no3'
       WRITE(numout,*) 'PFT=',test_pft
       WRITE(numout,*) 'anvf(:,ivm)=',anvf(test_grid,test_pft)
       WRITE(numout,*) 'p_O2(:,ivm)=',p_O2(test_grid,test_pft)
       WRITE(numout,*) 'p_O2air(:)=',p_O2air(test_grid)
       WRITE(numout,*) 'pb(:)=',pb(test_grid)
       WRITE(numout,*) 'd_O2(:)=',d_O2(test_grid)
       WRITE(numout,*) 'p_O2_resp(:)=',p_O2_resp(test_grid,:)
       WRITE(numout,*) 'mol_O2_resp(:)=',mol_O2_resp(test_grid,:)
       WRITE(numout,*) 'resp_soil(:,ivm)=',resp_sol(test_grid,test_pft)
       WRITE(numout,*) 'afps(:)=',afps(test_grid,:)        
    ENDIF
    !-
  
    !! 4. Physical removal of NH4 and NO3
    !
    ! 4.1 Adsorption of NH4 into clay
    !
    ! Reduce soil ammonium concentrations according to the equations in
    ! Table 4 of Li et al. (1992).  This represents adsorption onto soil clays.
    ! FixNH4 : Fraction of adsorbed NH4+
    ! FixNH4=[0.41-0.47 log(NH4) clay/clay_max]
    ! NH4+ concentration in the soil liquid, gN kg-1 soil 
    ! (p. 9774 of Li et al., 1992)
    ! but Zhang et al. (2002) Appendix A/B define NH4 as
    ! NH4+ in a soil layer in kgN ha-1
    ! In OCN, NH4 seems defined as kgN kg-1 or m-3 of water 
    ! Comment of N. Vuichard: I don't know which definition is the good one ...
    DO ivm = 1, nvm

       WHERE(soil_n_min(:,ivm,iammonium).GT.min_stomate)
          FixNH4(:) = MAX(0., (a_FixNH4 + b_FixNH4 * &
               MAX(0.,log10(soil_n_min(:,ivm,iammonium)*10 )))) &
               * MIN(clay(:),clay_max) / clay_max
       ELSEWHERE
          FixNH4(:) = 0.0   
       ENDWHERE

       ! In OCN, we do not multiply by FixNH4 but by FixNH4/(1+FixNH4)
       ! Comment of N. Vuichard: It is not clear if we should keep the 
       ! formulation of OCN or not
       ! So far, we keep the original formulation
       n_adsorbed(:,ivm) = MIN(soil_n_min(:,ivm,iammonium),&
            soil_n_min(:,ivm,iammonium) * FixNH4(:))
    ENDDO

    ! Pool update after ammonium adsorption. The previous line
    ! of code protects against significant negative values. The
    ! MAX(zero,x) below protect against insignificant negative 
    ! values that may cause some tests in stomate_litter to fail.
    ! Set soil_n_min to zero and use that value in the subsequent
    ! calculations. Note that towards the end n_absorbed is added
    ! back into soil_n_min
    WHERE(soil_n_min(:,:,iammonium) - n_adsorbed(:,:) .LT. zero)
       n_adsorbed(:,:) =  soil_n_min(:,:,iammonium)
       soil_n_min(:,:,iammonium) = zero
    ELSEWHERE
       soil_n_min(:,:,iammonium) = soil_n_min(:,:,iammonium) - n_adsorbed(:,:)
    ENDWHERE

    IF(printlev_loc>=4)THEN
       WRITE(numout,*) 'CHECK values after adsorption'
       WRITE(numout,*) 'soil_n_min ',soil_n_min(test_grid,test_pft,:)
       WRITE(numout,*) "n_adsorbed(:,:)",n_adsorbed(test_grid,test_pft)
    ENDIF
    !-

    ! 4.2 Leaching of NH4 and NO3

    WHERE((tmc_pft(:,:)+runoff_pft(:,:)) .NE. 0) 
    	f_drain(:,:) = (fracn_drainage*drainage_pft(:,:) + &
          fracn_runoff*runoff_pft)/(tmc_pft(:,:)+runoff_pft(:,:))
    ENDWHERE
    f_drain(:,:) = MAX(zero, MIN(f_drain(:,:),un))

    DO ivm = 1, nvm
       leaching(:,ivm,iammonium) =  leaching(:,ivm,iammonium) + &
            MIN(soil_n_min(:,ivm,iammonium) * f_drain(:,ivm), &
            soil_n_min(:,ivm,iammonium))
       leaching(:,ivm,initrate)  = leaching(:,ivm,initrate) + &
            MIN(soil_n_min(:,ivm,initrate)  * f_drain(:,ivm), &
            soil_n_min(:,ivm,initrate))
    ENDDO

  
    !! 5. Nitrification of NH4 in the oxygenated part of the soil
    !
    ! Nitrification refers to ammonium (NH4+) being oxidized to 
    ! nitrate (NO3-) under aeroboic (i.e., in the prescence of
    ! oxygen) conditions.  This sections thus reduces the
    ! concentration of ammonium in the soil.  These equations
    ! are taken from Zhang et al. (2002), Appendix A
    ! 
    ! 5.1 Effect of environmental factors (Temp, Humidity, pH)
    ! 5.1.1 Effect of soil moisture on nitrification
    !
    ! 202602 xw: the following variable is npts x nvm
    ! if swc_pft is 0, then fwnit is fw_nit_0 = -0.0243 firstly
    ! then becomses 0
    fwnit(:,:) = fw_nit_0 + fw_nit_1 * swc_pft(:,:) + fw_nit_2 * swc_pft(:,:)**2 + &
         fw_nit_3 * swc_pft(:,:)**3 + fw_nit_4 * swc_pft(:,:)**4
    fwnit(:,:) = MAX (0., MIN( un, fwnit(:,:) ) )
  
    ! 5.1.2 Effect of temperature on nitrification
    ! Note that Zhang et al 2002 fold in the factor 0.1 with the parameter
    ! for the term linear in temperature (ft_nit_1), while we group it 
    ! with the soil temperature as per the rest of the terms.  Checking 
    ! the parameter values, this leads to a first impression that they 
    ! differ by a factor of 10, when in reality the same overall result 
    ! is calculated.
    ! ft_nit_0 = -0.0233, note that temp_sol can be < 0, 0 or > 0
    ! ft_nit minimal value will be 0.001
    var_temp_sol(:) = temp_sol(:) * 0.1
    ft_nit(:) = ft_nit_0 + ft_nit_1 * var_temp_sol(:) + &
         ft_nit_2 * var_temp_sol(:)**2 + ft_nit_3 * var_temp_sol(:)**3 + &
         ft_nit_4 * var_temp_sol(:)**4
    ft_nit(:) = MAX (0.001, MIN( un, ft_nit(:) ) )
  
    ! 5.1.3 Effect of pH on nitrification
    ! fph_0 = -1.2314
    fph(:) = fph_0 + fph_1 * pH(:) + fph_2 * ph(:)**2
    fph(:) = MAX(0.0, fph(:) )
  
    ! 5.1.4 Effect of temperature on NO2 or NO production during nitrification
    ! I am not sure why the factor of 0.0025*NO3,N seems to be missing from
    ! the equation reported in the literature.
    !ftv_0 = 2.72 , ftv_2 = 9615.
    ftv(:) = ftv_0 **(ftv_1 - ftv_2 /(temp_sol(:) + tp_00) )
    ftv(:) = MAX(0.0, ftv(:) )
  
    DO ivm = 1, nvm

       ! 5.2 Actual nitrification rate per PFT NH4 soil pool - NH4 to NO3
       ! 
       ! Equation for the nitrification rate probably from 
       ! Schmid et al., 2001, Nutr. Cycl. Agro (eq.1)
       ! but the environmental factors used are from Zhang (2002) who 
       ! used an other equation
       ! I don't know how this can be mixed together ?
       ! In addition, the formulation mixed the one from Schmid with 
       ! the use of the anaerobic balloon defined by Li et al., 2000. 
       ! I don't know how this can be mixed together ?
       ! Last, in OCN, the default fraction of N-NH4 which is converted
       ! to N-NO3 appears equal to 2 day-1 (a factor "2" in the equation)
       ! - In Schmid et al., the nitrification rate at 20 ◦C and field 
       ! capacity (knitrif,20) was set to 0.2 d−1 (Ten time less...)
       ! k_nitrif = 0.2
       ! nitrification(:,ivm,i_nh4_to_no3) = MIN(fw(:) * fph(:) * &
       ! ft_nit(:) * k_nitrif * dt * & 
       ! soil_n_min(:,ivm,iammonium) * (1.0 - anvf(:,ivm)) , &
       ! soil_n_min(:,ivm,iammonium)- leaching(:,ivm,iammonium))
       nitrification(:,ivm,i_nh4_to_no3) = MIN(fwnit(:,ivm) * fph(:) * &
            ft_nit(:) * k_nitrif * dt * soil_n_min(:,ivm,iammonium) * &
            MAX(zero,(1.0 - anvf(:,ivm))), &
            soil_n_min(:,ivm,iammonium)-leaching(:,ivm,iammonium))
       
       ! Debug
       IF(printlev_loc>=4 .AND. ivm == test_pft)THEN
          WRITE(numout,*) 'CHECK values after nitrification nh4 to no3'
          WRITE(numout,*) 'PFT=',ivm
          WRITE(numout,*) 'nitrification(:,ivm,i_nh4_to_no3) ',&
               nitrification(test_grid,ivm,i_nh4_to_no3)
          WRITE(numout,*) 'fwnit=',fwnit(test_grid,ivm)
          WRITE(numout,*) 'fph=',fph(test_grid)
          WRITE(numout,*) 'ft_nit=',ft_nit(test_grid)
          WRITE(numout,*) 'anvf=',anvf(test_grid,ivm)
       ENDIF
       
       !
       ! 5.3 Emission of N2O during nitrification - NH4 to N2O
       ! 
       nitrification(:,ivm,i_nh4_to_n2o) = ftv(:) * swc_pft(:,ivm) * &
            n2o_nitrif_p * nitrification(:,ivm,i_nh4_to_no3)
       
       ! Debug
       IF(printlev_loc>=4 .AND. ivm==test_pft)THEN
          WRITE(numout,*) 'CHECK values after nitrification nh4 to n2o'
          WRITE(numout,*) 'nitrification(:,ivm,i_nh4_to_n2o) ',&
               nitrification(test_grid,ivm,i_nh4_to_n2o)
          WRITE(numout,*) 'ftv=',ftv(test_grid)
          WRITE(numout,*) 'swc_pft=',swc_pft(test_grid,:)
       ENDIF
       
       nitrification(:,ivm,i_nh4_to_n2o) = &
            MIN(nitrification(:,ivm,i_nh4_to_no3), &
            nitrification(:,ivm,i_nh4_to_n2o))       
     
       ! 5.4 Production of NO during nitrification - NH4 to NO
       nitrification(:,ivm,i_nh4_to_no) = ftv(:) * no_nitrif_p * &
            nitrification(:,ivm,i_nh4_to_no3)

       IF(printlev_loc>=4 .AND. ivm == test_pft)THEN
          WRITE(numout,*) 'CHECK values after nitrification nh4 to no'          
          WRITE(numout,*) 'nitrification(:,ivm,i_nh4_to_no) ',&
               nitrification(test_grid,ivm,i_nh4_to_no)
       ENDIF

       ! 5.5 Production of NO from chemodenitrification
       ! Chemodenitrification is the conversion of NO2 to NO through a purely
       ! chemical process, dependent on temperature, soil pH, and concentration of
       ! NO2.  This is the final step of the process NH4+ -> NO3- -> NO2- -> NO
       ! Based on Chem_NO in Kesik et al.(2005)
       ! BUT Kesik et al. used NO2 concentration and not NO3 production 
       ! as it is done in OCN and modification of a multiplicative constant
       ! from 300 (Kesik) to 30 (OCN) without clear motivation - I don't &
       ! how this is reliable 
       ! [Matt also questions this equation...as it is written below, it is the amount
       ! of NH4 being converted to NO3 multiplied by the decomposition of NO2 to NO, but
       ! there is no NO3 to NO2 conversion included, which implies that all NO3 is
       ! instantaneously converted to NO2...which is nonsense based on the fact that an NO3
       ! pool exists]
       nitrification(:,ivm,i_nh4_to_no) = nitrification(:,ivm,i_nh4_to_no) + &
            ( chemo_0 * chemo_1 * exp(chemo_ph0 * pH(:) ) * & 
            exp(chemo_t0/((temp_sol(:)+tp_00)*RR))) * &
            nitrification(:,ivm,i_nh4_to_no3)
       
       nitrification(:,ivm,i_nh4_to_no) = &
            MIN(nitrification(:,ivm,i_nh4_to_no3)-&
            nitrification(:,ivm,i_nh4_to_n2o), &
            nitrification(:,ivm,i_nh4_to_no))
       
       
       IF(printlev_loc>=4 .AND. ivm == test_pft)THEN
          WRITE(numout,*) 'CHECK values after nitrification nh4 to no chemodenitrification'
          WRITE(numout,*) 'nitrification(:,ivm,i_nh4_to_no) ',nitrification(test_grid,ivm,i_nh4_to_no)
          WRITE(numout,*) 'pH(:)=',pH(test_grid)
          WRITE(numout,*) 'temp_sol(:)=',temp_sol(test_grid)
       ENDIF
       
       ! In OCN, NO production and N2O production is deducted from NO3 production due to NH4
       nitrification(:,ivm,i_nh4_to_no3) = MAX(zero,nitrification(:,ivm,i_nh4_to_no3) - &
            (nitrification(:,ivm,i_nh4_to_no) + nitrification(:,ivm,i_nh4_to_n2o)))
    ENDDO
     
  
    !! 6. Denitrification processes
    !
    ! Denitrification is the conversion of any oxidized N comound (NO2-, NO, N2O)
    ! to the final product of gaseous nitrogen (N2).
    !
    ! Denitrifiers are organisisms responsible for denitrification (bacteria?).
    !
    ! Li et al, 2000, JGR Table 4 eq 1, 2 and 4, ignoring NO2 (similar to NO)
    ! Comment of S. Zaehle: "includes treatment of bacteria dynamics, but I do not have any confidence in this;
    ! at least it provides rates that appear resonable." 
  
    ! 6.1 Temperature response of relative growth rate of total denitrifiers
    ft_denit(:) = 2.**((temp_sol(:)-22.5)/10.)
  
    ! 6.2 Soil pH response of relative growth rate of total denitrifiers
    ! See also comment about parenthesis' position in Thesis of Vincent Prieur (page 50)
    fph_no3(:) = 1.0 - 1.0 / ( 1.0 + EXP((pH(:)-fph_no3_0)/fph_no3_1)) 
    fph_no(:) = 1.0 - 1.0 / ( 1.0 + EXP((pH(:)-fph_no_0)/fph_no_1)) 
    fph_n2o(:) = 1.0 - 1.0 / ( 1.0 + EXP((pH(:)-fph_n2o_0)/fph_n2o_1))
  
    ! Half saturation value of N oxides (kgN/m3), Kn
    ! Unit conversion factor from kgN/m3 to gN/m2
    ! OCN used the value 0.087, but 0.083 is listed in Li et al (2000).
    ! In OCN, we account for max_eau_eau. Not clear if this is correct or not.
    Kn_conv(:,:) = Kn * tmc_pft(:,:)    
    
    ! Relative growth rate of Nox denitrifiers
    ! Eq.1 Table 4 Li et al., 2000 - but there is an error because Eq. 1 it does not 
    ! account for Temp and ph responses.  
    ! In OCN it does not account for [DOC], though Li et al. (2000) does.  We keep
    ! the formulation from OCN.  In addition, in OCN the term mu_nox(max) is missing 
    ! in the equation that defines the relative growth rate of total denitrifiers 
    ! (dn_bact) - we add the term back in here.
    !
    ! 6.3  Maximum Relative growth rate of Nox denitrifiers (hour-1): mu_no3, mu_no2, mu_no
    ! kfwdenit = -5. in constants.f99
    DO ivm = 1, nvm
       WHERE((tmc_pft(:,ivm)/(zmaxh*1000.)) .LT. mcfc_hydrol(:))
          fwdenit(:,ivm) = fwdenitfc * &
               exp(kfwdenit * (mcfc_hydrol(:)-tmc_pft(:,ivm)/(zmaxh*1000.))/mcfc_hydrol(:) )
       ELSEWHERE
          fwdenit(:,ivm) = fwdenitfc + (1 - fwdenitfc) * &
               (tmc_pft(:,ivm)/(zmaxh*1000.) - mcfc_hydrol(:)) / (mcs_hydrol(:)-mcfc_hydrol(:))
       ENDWHERE

       WHERE((soil_n_min(:,ivm,initrate) + Kn_conv(:,ivm)) .GT. min_stomate)
          mu_no3(:) = mu_no3_max * soil_n_min(:,ivm,initrate) / (soil_n_min(:,ivm,initrate) + Kn_conv(:,ivm))
       ELSEWHERE
          mu_no3(:) = zero
       ENDWHERE
       WHERE((soil_n_min(:,ivm,inox) + Kn_conv(:,ivm)*fact_kn_no) .GT. min_stomate)
          mu_no(:)  = mu_no_max  * soil_n_min(:,ivm,inox) / (soil_n_min(:,ivm,inox) + Kn_conv(:,ivm)*fact_kn_no)
       ELSEWHERE
          mu_no(:) = zero 
       ENDWHERE
       
       WHERE((soil_n_min(:,ivm,initrous) + Kn_conv(:,ivm)*fact_kn_n2o) .GT. min_stomate)
          mu_n2o(:) = mu_n2o_max * soil_n_min(:,ivm,initrous) / (soil_n_min(:,ivm,initrous) + Kn_conv(:,ivm)*fact_kn_n2o)
       ELSEWHERE
          mu_n2o(:) = zero
       ENDWHERE
  
       ! stocks of all N species (NO3, NO, N2O) (gN/m**2)
       sum_n(:) = soil_n_min(:,ivm,inox) + soil_n_min(:,ivm,initrous) + & 
            soil_n_min(:,ivm,initrate)
    
       WHERE(sum_n(:).GT.min_stomate)

	   ! Maint_c maintenance coefficient of carbon (kgC/kgC/h)
           ! Maint_c = 0.0076
           ! Yc maximum growth yield on soluble carbon (kgC/kgC)
           ! Yc = 0.503
           ! denitrifier bacterial population change
           ! Eq. 3 (Rg-Rd) of Table 4 Li et al., 2000
           ! In OCN, dn_bact is ranged between ,0.0005* &
           ! som(:,iactive,ivm,icarbon) and 
           ! 0.25 * som(:,iactive,ivm,icarbon)) - no clear justification.
           ! This is not kept in the current version
           ! In addition, in OCN, Only dn_bact is defined, not bact. &
           ! I don't know how it may work
           !+++CHECK+++
           ! Nitrogen in bacteria is calculated, but we do not specify where this
           ! nitrogen is taken from. Because bacteria are not included in the
           ! mass balance check, this is not causing a problem but in these
           ! lines we create bacteria compose from pure N. If bact is ever used, 
           ! their carbon use must be accounted for. In revision 3460 bact is further 
           ! determined by the max-min constraint.
           ! Original code
           !dn_bact(:) = ( anvf(:,ivm) * ( mu_no3 + mu_no + mu_n2o ) - &
           !     Maint_c * Yc ) * bact(:,ivm) * 24. * dt
           !bact(:,ivm) = bact(:,ivm) + dn_bact(:)
           !bact(:,ivm) = MAX(MIN(bact(:,ivm),0.25*som(:,iactive,ivm,icarbon)),&
           !              0.0005*som(:,iactive,ivm,icarbon))
          
          ! This appears to be an assumption that the denitrifier biomass is
          ! 0.005% of the active carbon soil organic matter.  Unclear where
          ! this comes from.
          bact(:,ivm) = cte_bact*som(:,iactive,ivm,icarbon)
          !+++++++++++


          ! 6.4 Consumption rate of N oxides
          ! Table 4 Li et al.(2000)
          ! Based on the maximum growth rate on N oxides (Y_nox), maintainance
          ! coefficient on N oxides (M_nox), denitrifier biomass (bact), 
          ! relative growth rate of NOX denitrifiers (mu_nox), concentration
          ! of NOX in the soil (soil_n_min), total concentration of all nitrogen
          ! species in the soil (sum_n), the amount of NOX leached out of the
          ! soil (leaching).  The fwdenit, ft_denit, and fph_nox do not come
          ! from Table 4.  Conversion from (per hour) to (per timestep) is 
          ! handled by 24*dt at the end.
          !
          ! The reactions are NO3 -> (NO + NO2) -> N2O -> N2
          !
          ! NO3 consumption
          ! In OCN, multiplication by 0.1 of the NO3 consumption - no justification
          ! This is not kept in the current version
          denitrification(:,ivm,i_no3_to_nox) = MAX(zero, MIN(soil_n_min(:,ivm,initrate)-leaching(:,ivm,initrate), &
               fwdenit(:,ivm) * ft_denit(:) * fph_no3(:) * ( mu_no3(:) / Y_no3 + M_no3 * soil_n_min(:,ivm,initrate) / sum_n(:)) * &
               bact(:,ivm) * 24. * dt ))
          !
          ! NO consumption
          denitrification(:,ivm,i_nox_to_n2o) = MIN(soil_n_min(:,ivm,inox), &
               fwdenit(:,ivm) * ft_denit(:) * fph_no(:) * ( mu_no(:) / Y_no + M_no * soil_n_min(:,ivm,inox) / sum_n(:)) * &
               bact(:,ivm) * 24. * dt )        
          !
          ! N2O consumption
          denitrification(:,ivm,i_n2o_to_n2) = MIN(soil_n_min(:,ivm,initrous), &
               fwdenit(:,ivm) * ft_denit(:) * fph_n2o(:) * ( mu_n2o(:) / Y_n2o + M_n2o * soil_n_min(:,ivm,initrous) / sum_n(:)) &
               * bact(:,ivm) * 24. * dt )
  
          ! Dynamcis on denitrifier bacterial population change is not used, due to
          ! lack of documention in OCN and lack of clarity in Li et al (2000).
          ! In addition, OCN used dn_bact, but bact is used here.  Not clear what
          ! the consequences of this are.

        ENDWHERE
  
    ENDDO

    ! 7. Plant uptake of ionic nitrogen
    !
    ! This section comes from Section 3 of the supporting material of Zaehle & Friend (2010).
    ! Only ammonium and nitrate are uptaken by plants.
    ! 
    ! Comment from OCN : Temperature function of uptake is similar to SOM decomposition
    ! to avoid N accumulation at low temperatures 
    ! In addition in the SM of Zaehle & Friend, 2010 P. 3 it is mentioned
    ! "The uptake rate is observed to be sensitive to root temperature 
    ! which is included as f(T), 
    ! thereby following the temperature sensitivity of net N mineralization"
    IF (test_Q10_on_plant_uptake) THEN
       ft_uptake(:) = control_temp_func (npts, temp_sol(:)+tp_00)
    ELSE  
       ft_uptake(:) = 1.0  
    ENDIF
    ! Comment of N. Vuichard: Plenty of parameter values used in the formulation of the plant uptake are 
    ! not traceable to anything. I tend to rely much of the param values to the values reported 
    ! in the reference publications
    
    ! Vmax of nitrogen uptake (in umol (g DryWeight_root)-1 h-1)
    ! 
    ! In OCN, the same values are used both for uptake of NH4+ and NO3-
    ! See p. 3 of the SM of Zaehle & Friend, 2010:
    ! "As a first approximation average values for vmax, kNmin and KNmin
    ! are assumed for all PFT (Table S1), and for both ammonium and nitrate"
    ! However based on the two papers of Kronzucker et al. (1995, 1996), it
    ! seems that vmax should be much higher for NH4+ than for NO3- . 
    ! The value of Vmax reported in Table S1 of Zaehle & Friend, 2010 is 5.14 ugN (g-1C) d-1
    ! Comment of N. Vuichard : I can not relate the value of vmax in OCN expressed in 
    ! (umol (g DryWeight_root)-1 h-1) (vmax=3) to the one reported in Table S1 :  5.14 ugN (g-1C) d-1
    ! The use of the conv_fac conversion factor used in OCN should help to convert (umol (g DryWeight_root)-1 h-1)
    ! in (gN (g-1C) tstep-1). conv_fac is defined as 24. * dt / 2. / 1000000. * 14. 
    ! - 24 conversion of hour to day
    ! - dt conversion of day to dt 
    ! - 1/2 conversion of g(DryWeight) to gC
    ! - 1000000. conversion of ug to g
    ! - 14 conversion of umol to ugN
    ! So using conv_fac, vmax expressed in ugN (g-1C) d-1 should be equal to 
    ! 3*24./2.*14. = 504 ugN (gC)-1 d-1 not 5.14
    ! In addition, I think there is an error in the conversion from (gDW)-1 to (gC)-1
    ! When expressed per gC of root, the N uptake should be twice more 
    ! than the one expressed per gDW of root, not twice less. So one should
    ! multiply by 2., not divide by 2. (so vmax = 2016 ug (gC)-1 d-1 )
    !
    ! Vmax of nitrogen uptake (in umol (g DryWeight_root)-1 h-1)  
    ! Conversion from (umol (gDW)-1 h-1) to (gN (gC)-1 timestep-1)
    conv_fac_vmax= 24. * dt * 2. / 1000000. * 14. 
  
    ! minimal and maximal NC ratio of leaf (gN / gC)
    ! used in the response of N uptake to NC ratio
    ! see eq. 10 , p. 3 of SM of Zaehle & Friend, 2010
    nc_leaf_max(:,1)     = 0.
    nc_leaf_max(:,2:nvm) = 1. / cn_leaf_min_2D(:,2:nvm)
  
    nc_leaf_min(:,1)     = -1./min_stomate
    nc_leaf_min(:,2:nvm) = 1. / cn_leaf_max_2D(:,2:nvm)
  
    !+++CHECK+++
    ! In OCN, lab_n and lab_c are not calculated for bare soil. Why do we do it
    ! here?
    DO ivm = 1, nvm

       ! NCplant (gN / gC)
       ! See equation (9) p. 3 of SM of Zaehle & Friend, 2010. Zaehle and Friend
       ! did not account for the reserve carbon and nitrogen. Since then we
       ! implemented an optimal distribution between the labile and reserve pool
       ! (in stomate_growth_fun_all.f90) so there can be a lot of N in the
       ! reserve pool. The reserve pool is now included in the calculation of the
       ! C/N ratio of the plant.
       lab_n(:) = SUM(circ_class_biomass(:,ivm,:,ilabile,initrogen)*circ_class_n(:,ivm,:),2) + &
            SUM(circ_class_biomass(:,ivm,:,iroot,initrogen)*circ_class_n(:,ivm,:),2) + & 
            SUM(circ_class_biomass(:,ivm,:,ileaf,initrogen)*circ_class_n(:,ivm,:),2) + &
            SUM(circ_class_biomass(:,ivm,:,icarbres,initrogen)*circ_class_n(:,ivm,:),2)
       lab_c(:) = SUM(circ_class_biomass(:,ivm,:,ilabile,icarbon)*circ_class_n(:,ivm,:),2) + &
            SUM(circ_class_biomass(:,ivm,:,iroot,icarbon)*circ_class_n(:,ivm,:),2) + &  
            SUM(circ_class_biomass(:,ivm,:,ileaf,icarbon)*circ_class_n(:,ivm,:),2) + &
            SUM(circ_class_biomass(:,ivm,:,icarbres,icarbon)*circ_class_n(:,ivm,:),2)
       
       WHERE(lab_c(:).GT.min_stomate)
          NCplant(:) = lab_n(:) / lab_c(:)
       ELSEWHERE
          NCplant(:) = fcn_root(ivm)/cn_leaf_init_2D(:,ivm) 
       ENDWHERE
       
       ! Nitrogen demand responds to N/C ratio of the labile pools for 
       ! each PFT zero uptake at a NC in the labile pool corresponding 
       ! to maximal leaf NC. If the N-stress increases, f_NCplant increases
       ! up to a maximum of 1 for an NC corresponding to a minimal leaf NC
       ! See equation (10) p. 3 of SM of Zaehle & Friend, 2010  
       f_NCplant(:,ivm)=MIN(MAX( ( nc_leaf_max(:,ivm) - NCplant(:) ) / & 
           ( nc_leaf_max(:,ivm) - nc_leaf_min(:,ivm) ), 0. ), 1.) 
       
       ! Plant N uptake (gN m-2 tstep-1)
       ! See eq. (8) p. 3 of SM of Zaehle & Friend, 2010
       ! In OCN, there was a multiplicative factor "2" that I can not explain
       ! It may partly compensate the error mentioned above about conv_fac_vmax
       !
       ! Coefficient K_N_min (umol per litter)
       !
       ! It corresponds to the [NH4+] (resp. [NO3-]) for which the Nuptake 
       ! equals vmax/2. Kronzucker, 1995 reports values that range between
       ! 20 and 40 umol for NH4+ uptake OCN seems to use 30 umol for both
       ! NH4+ and NO3-. The value used is 0.84 expressed in (gN m-2)
       ! In Kronzucker, 1995 and 1996, the concentrations are always 
       ! expressed in umol. No clear reference about the standard volume it
       ! is related to (litter ?).
       ! Coefficient K_N_min (umol per litter)
       ! K_N_min(:) = 30
       ! Conversion factor from (umol per litter) to (gN m-2) 
       conv_fac_concent = 14.0 * 1e3 * 1e-6 * 2.0

       ! 14   : molar mass for N (gN mol-1)
       ! 10^3 : conversion factor (dm3 to m3)
       ! 10-6 : conversion factor (ug to g)
       ! 2    : m3 per m2 of soil
       !
       ! Comment of N. Vuichard : I wonder why it assumes 2 m3 per m2 of 
       ! soil. The depth of the soil is 2 meter. But it doesn't mean that
       ! all the column contains only water, there is also soil, no ?
       ! The question is : what is the volume that is considered for 
       ! the concentration of NH4+ and NO3-. Is it a volume of soil or of
       ! solution ? If it is a volume of solution, I would suggest to
       ! multiply by a factor corresponding to the relative volume 
       ! of water within the soil. - Not done yet. 
       ! 
       ! K_N_min * conv_fac_concent = 0.84
       !
       ! Coefficient low_K_N_min ((gN m-2)-1)
       !
       ! See eq. 8 of SM of Zaehle et al. (2010) and Table S1, In table S1, it
       ! is defined as "Rate of N uptake not associated with Michaelis- Menten
       ! Kinetics" with a value of 0.05. It is mentioned also (unitless) but to
       ! my opinion (N. Vuichard) it should have the same unit that 1/K_N_min or
       ! 1/N_min, so ((gN m-2)-1)
       ! +++CECK+++
       ! Comment of N. Vuichard
       ! So far, I cannot relate the value of low_K_N_min (0.05) to any reference
       ! especially Kronzucker (1996). If I refer to Figure 4 of Kronzucker showing
       ! the NH4+ influx as a function of NH4+ concentration, the slope of the
       ! relationship could be used to define Vmax*low_K_N_min. For a concentration
       ! of NH4+ of 50 mmol, the influx equals 35 umol g-1 h-1. For a concentration
       ! of NH4+ of 20 mmol, the influx equals 17 umol g-1 h-1.
       ! slope = 10-3 * (35 - 17) / (50 - 20) = 10-3 * 18 / 30 = 0.0006 g-1 h-1 
       ! low_K_N_min = slope / Vmax = 0.0006 / 3 = 0.0002 (umol)-1 
       ! using the Conversion factor from (umol per litter) to (gN m-2) 
       ! conv_fac_concent = 14 * 1e3 * 1e-6 * 2, one should get 
       ! low_K_N_min = 0.0002 / ( 14 * 1e3 * 1e-6 * 2 ) = 0.007 ((gN m-2)-1)
       ! The value of 0.007 does not match with the one in OCN (0.05) - 
       ! This needs to be clarified
       ! ++++++++++
       ! low_K_N_min ((umol)-1)
  
       ! In eq. 8, I (N. Vuichard) think there is an error. It should not be
       ! (Nmin X KNmin) but (Nmin + KNmin). The source code of OCN is correct
       
       ! Plant uptake of ammonium and nitrate, To allow trees to take up 
       ! some N during winter uptake was made a function of longterm root 
       ! biomass (last 3 years). croot_longterm now has a value even if there
       ! is no vegetation (i.e. crops), plant_n_uptake, however, requires
       ! the availability of plants
       DO iion = 1,nionspec
          ! Run the model without the additional control of sugar_load
          ! sugar_load was removed - See comments in above
          WHERE (SUM(circ_class_n(:,ivm,:),2).GT.min_stomate)
             plant_n_uptake(:,ivm,iion) = vmax_uptake(ivm,iion) * &
                  conv_fac_vmax * croot_longterm(:,ivm) * ft_uptake(:) * &
                  soil_n_min(:,ivm,iion) * ( low_K_N_min(iion) / &
                  conv_fac_concent + 1. / ( soil_n_min(:,ivm,iion) + &
                  K_N_min(iion) * conv_fac_concent ) ) * &
                  f_NCplant(:,ivm)
          ELSEWHERE
             plant_n_uptake(:,ivm,iion) = zero
          ENDWHERE
       ENDDO

       ! Debug
       IF ((printlev_loc>=4).AND.(ivm==test_pft)) THEN
          WRITE(numout,*)'in nitrogen_dyn after plant_uptake, ivm', ivm  
          WRITE(numout,*) 'plant_n_uptake ', &
               plant_n_uptake(test_grid,test_pft,:)
          WRITE(numout,*) 'min plant_n_uptake ammonium', soil_n_min(test_grid,ivm,iammonium) - &
            (leaching(test_grid,ivm,iammonium) + nitrification(test_grid,ivm,i_nh4_to_no3) + &
            nitrification(test_grid:,ivm,i_nh4_to_no) + nitrification(test_grid,ivm,i_nh4_to_n2o))   
          WRITE(numout,*) 'min plant_n_uptake nitrate',soil_n_min(test_grid,ivm,initrate) - &
            leaching(test_grid,ivm,initrate) - denitrification(test_grid,ivm,i_no3_to_nox) + &
            nitrification(test_grid,ivm,i_nh4_to_no3)
          WRITE(numout,*) ' ccb isapbelow', &
               circ_class_biomass(test_grid,ivm,:,isapbelow,:)
          WRITE(numout,*) ' ccb iroot', &
               circ_class_biomass(test_grid,ivm,:,iroot,:)
          WRITE(numout,*) 'ft_uptake(test_grid)',ft_uptake(test_grid)
          WRITE(numout,*) 'soil_n_min ',soil_n_min(test_grid,test_pft,:)
          WRITE(numout,*) 'low_K_N_min', low_K_N_min(:)
          WRITE(numout,*) 'K_N_min(iammonium)', K_N_min(:)
          WRITE(numout,*) 'nc_leaf_max(ivm)',nc_leaf_max(test_grid,ivm)
          WRITE(numout,*) 'nc_leaf_min(ivm)',nc_leaf_min(test_grid,ivm)
          WRITE(numout,*) 'leaching ',leaching(test_grid,test_pft,:)
          WRITE(numout,*) 'nitrification ',nitrification(test_grid,test_pft,:)
          WRITE(numout,*) 'bact',bact(test_grid,ivm) 
          WRITE(numout,*) 'denitrification ',denitrification(test_grid,test_pft,:)

          WRITE(numout,*) 'tmc_pft',tmc_pft(test_grid,test_pft)
       ENDIF
       !-

       ! Calculate plant uptake
       ! This code should prevent plant_n_uptake from being negative. This is
       ! essential to prevent errors in stomate_growth_fun_all
       plant_n_uptake(:,ivm,iammonium) = MAX(MIN(soil_n_min(:,ivm,iammonium) - & 
            (leaching(:,ivm,iammonium) + nitrification(:,ivm,i_nh4_to_no3) + &
            nitrification(:,ivm,i_nh4_to_no) + nitrification(:,ivm,i_nh4_to_n2o)), &
            plant_n_uptake(:,ivm,iammonium)), zero)

       plant_n_uptake(:,ivm,initrate) = MAX(MIN(soil_n_min(:,ivm,initrate) - &
            leaching(:,ivm,initrate) - denitrification(:,ivm,i_no3_to_nox) + &
            nitrification(:,ivm,i_nh4_to_no3), plant_n_uptake(:,ivm,initrate)), zero)

    END DO

    ! Loss of ionic forms of N through drainage and gaseous forms through
    ! volatilisation (the emission part should eventually be calculated in 
    ! Sechiba's diffuco routines
    DO ivm = 1,nvm

       ! 8. Loss of N through drainage and gaseous emission
       ! (the emission part should eventually be calculated in Sechiba's diffuco routines)
  
       ! 8.1 Loss of NH4 due to evaporation of NH3
       ! NH4+ is first converted to NH3, and then it evaporates from the soil.
       ! These equations come from Li et al. (1992), Table 4, in
       ! addition to Appendix A of Zhang et al. (2002)
  
       ! Current dissociation of [NH3] to [NH4+]
       ! See Table 4 of Li et al. 1992 and Appendix A of Zhang et al. 2002
       ! log(K_NH4) - log(K_H20) = log(NH4/NH3) + pH
       ! See also formula in "DISSOCIATION CONSTANTS OF INORGANIC ACIDS AND 
       ! BASES" pdf file
       ! pK_H2O = -log(K_H2O) = 14
       ! pk_NH4 = -log(K_NH4) = 9.25
       ! pK_H2O - pK_NH4 = log([NH4+]/[NH3]) + pH
       ! [NH4+]/[NH3] = 1O^(pK_H2O - pK_NH4 - pH) = 10^(4.75-pH)
       
       ! In OCN, one makes use of frac_nh3. That should be the NH3/NH4 
       ! ratio. It is defined as:
       ! frac_nh3(:) = 10.0**(4.25-pH(:)) / (1. + 10.0**(4.25-pH(:)))
       
       ! CommentS of N. Vuichard : I have several interrogations about 
       ! the equation and value used
       ! in OCN. But I have also concerns about the formulation of Li et al. ... 
       ! 1/ 
       ! The formulation of Li et al. doesn't match with the formulas in 
       ! "DISSOCIATION CONSTANTS OF INORGANIC ACIDS AND BASES" or with 
       ! the formulas
       ! of http://www.onlinebiochemistry.com/obj-512/Chap4-StudNotes.html
       ! To my opinion, one should replace log(K_NH4+) by log(K_NH3) in the 
       ! formulation of Li et al.
       ! with the relationship pKa + pKb = pKwater (where a and b are acid
       ! and base)
       ! this leads to [NH4+]/[NH3] = 10^(pK_NH4 - pH) = 10^(9.25 - pH)
       ! 2/
       ! Wether I'm right or not about 1/, I don't understand the value 
       ! used in OCN (4.25). 
       ! It should be either 4.75 or 9.25 but 4.25 looks strange
       ! 3/ 
       ! OCN used a formulation for [NH3]/[NH4+] of the type: X/(1+X) 
       ! with X=[NH3]/[NH4+]
       ! This leads to X/(1+X)=([NH3]/[NH4+])/(([NH4+]/[NH4+])+([NH3]/[NH4+]))
       ! or X/(1+X) = [NH3] / ( [NH3] + [NH4+] )
       ! This means that the value stored in soil_n_min(:,:,iammonium) 
       ! corresponds to the total
       ! N of both [NH4+] and [NH3]. This makes sense to my opinion. 
       ! But I wonder if one should not
       ! account for this partitioning between [NH4+] and [NH3] in other 
       ! processes. To check.
       ! 4/ 
       ! The X value should relate to [NH3]/[NH4+] but to my opinion 
       ! the X value used 
       ! in the equation in OCN corresponds to [NH4+]/[NH3]. Is this a bug ? 
   
       ! In conclusion, I propose the formulation
       frac_nh3(:) = 10.0**(0.15*(pH(:)-pk_NH4)) / (1. + 10.0**(0.15*(pH(:)-pk_NH4)))
  
       ! This seems a patch added to OCN in order to (as mentioned in OCN) 
       ! reduced emissions at low concentration (problem of only one soil layer)
       ! high conentrations are usually associated with fertiliser 
       ! events -> top layer
       ! and thus increased emission
       ! NOT ACTIVATE HERE 
       ! emm_fac(:) = 0.01
       ! WHERE(soil_n_min(:,ivm,iammonium).GT.0.01.AND.&
       !     soil_n_min(:,ivm,iammonium).LE.4.)
       !      emm_fac(:) = MAX(0.01,1-exp(-(soil_n_min(:,ivm,iammonium)/1.75)**8)) 
       ! ENDWHERE       
       ! WHERE(soil_n_min(:,ivm,iammonium).GT.4.)
       !      emm_fac(:)=1.0
       ! ENDWHERE
       
  
       ! 8.2 Volatilisation of gasous species, Table 4, Li et al. 2000 I,
       !     using diffusivity of oxigen in air as a surrogate (from OCN)
       !     assumes no effect of air concentration on diffusion
       !     takes standard depth as reference
       
       ! Clay limitation
       ! F_clay_0 = 0.13
       ! F_clay_1 = -0.079
       F_clay(:) = F_clay_0 + F_clay_1 * clay(:)
  
       ! In OCN, one used formulation of the type
       ! emission(:,ivm,inox-1) = & 
       !       MIN( d_ox(:) * soil_n_min(:,ivm,inox) * (0.13-0.079*clay(:)) * 
       !       dt / z_decomp
       ! It should have the unit d_ox * soil_n_min * dt / z_decomp
       !                         m2 day-1 * gN m-2 * day / m
       !                         gN / m which is not homogeneous with 
       ! the unit expected (gn m-2)
       ! To my opinion, one should not use soil_n_min (gN m-2) but a 
       ! volumetric concentration (gN m-3)
       ! But I'm not clear what is the appropriate volume to consider 
       ! (volume of soil ?)
  
       ! NH4 emission (gN m-2 per time step)
       emission(:,ivm,iammonium) = D_s(:,ivm) * emm_fac * frac_nh3(:) * &
            soil_n_min(:,ivm,iammonium) / zmaxh &
            * F_clay(:) / z_decomp * dt
                   
       ! NO3 emission 
       emission(:,ivm,initrate) = zero
         
       ! NO2 emission (gN m-2 per time step)
       emission(:,ivm,inox) = D_s(:,ivm) * soil_n_min(:,ivm,inox) / &
            zmaxh * F_clay(:) / z_decomp * dt
         
       ! N2O emission (gN m-2 per time step)
       emission(:,ivm,initrous) = D_s(:,ivm) * soil_n_min(:,ivm,initrous) /&
            zmaxh * F_clay(:) / z_decomp * dt
                         
       ! N2 emission (gN m-2 per time step)
       emission(:,ivm,idinitro) = D_s(:,ivm) * soil_n_min(:,ivm,idinitro) /&
            zmaxh * F_clay(:) / z_decomp * dt
              
    ENDDO

    ! Calculate the emissions of the different N-species
    emission(:,:,iammonium) = MAX(zero, MIN(emission(:,:,iammonium), &
         soil_n_min(:,:,iammonium) - nitrification(:,:,i_nh4_to_no3) &
         - nitrification(:,:,i_nh4_to_no) - nitrification(:,:,i_nh4_to_n2o) &
         - leaching(:,:,iammonium) - plant_n_uptake(:,:,iammonium)))  
    emission(:,:,inox) = MAX(zero, MIN(emission(:,:,inox), &
         soil_n_min(:,:,inox) + nitrification(:,:,i_nh4_to_no) & 
         + denitrification(:,:,i_no3_to_nox) - denitrification(:,:,i_nox_to_n2o)))
    emission(:,:,initrous) = MAX(zero, MIN(emission(:,:,initrous), &
         soil_n_min(:,:,initrous) + nitrification(:,:,i_nh4_to_n2o) & 
         + denitrification(:,:,i_nox_to_n2o) - denitrification(:,:,i_n2o_to_n2)))
    emission(:,:,idinitro) = MAX(zero,MIN(emission(:,:,idinitro), &
         soil_n_min(:,:,idinitro)  + denitrification(:,:,i_n2o_to_n2)))

    !
    ! 9. Update pools
    !

    !+++CHECK+++
    ! To my [WHO?] opinion, for a better consistency, I would recommand
    ! to consider the leaching separately when calculating the ammonium
    ! and nitrate budget. Leaching is calculated from sechiba. Might be
    ! better to remove leaching at the top of the routine especially due
    ! to the later calculation of NH4+ and NO3- concentration that will
    ! vary with the soil water content. Let's imagine that from one time
    ! step to another, the change in soil water content is only due to
    ! leaching. We don't want that the NH4+ and NO3- concentration vary
    ! from one time step to the other. The best way to avoid this is to
    ! remove first the leaching from the NH4+ and NO3- pools
    ! THIS IS NOT DONE YET. 
    ! Note that deposition is added to leaching if there is no active soil
    ! carbon pool. If leaching is dealt with up front, this issues should
    ! be considered when checking the mass balance
    !+++++++++++
  
    ! 9.1 update pools of nitrogen in the soil from nitrification and 
    ! denitrification, plant uptake, leaching and volatile emissions,
    ! desorption from clay, and net mineralisation
    !
    ! In my opinion, for a better consistency, I would recommand to consider the leaching separately 
    ! when calculating the ammonium and nitrate budget. Leaching is calculated from sechiba
    ! Might be better to remove leaching at the top of the routine especially due to the later calculation of
    ! NH4+ and NO3- concentration that will vary with the soil water content
    ! Let's imagine that from one time step to another, the change in soil water content is only due to leaching
    ! We don't want that the NH4+ and NO3- concentration vary from one time step to the other. The best way to avoid
    ! this is to remove first the leaching from the NH4+ and NO3- pools
    ! THIS IS NOT DONE YET

    IF(printlev_loc>=4)THEN
       WRITE(numout,*) 'CHECK values before update'
       WRITE(numout,*) 'nitrification ',nitrification(test_grid,test_pft,:)
       WRITE(numout,*) 'denitrification ',denitrification(test_grid,test_pft,:)
       WRITE(numout,*) 'leaching ',leaching(test_grid,test_pft,:)
       WRITE(numout,*) 'emission ',emission(test_grid,test_pft,:)
       WRITE(numout,*) 'plant_n_uptake ',plant_n_uptake(test_grid,test_pft,:)
       WRITE(numout,*) 'mineralisation ',n_mineralisation(test_grid,test_pft)
       WRITE(numout,*) 'immob ',immob(test_grid,test_pft)
       WRITE(numout,*) 'soil_n_min',soil_n_min(test_grid,test_pft,:)
    ENDIF

    ! Nitrogen species: iammonium
    ! Note that n_absorbed was removed at the start of this subroutine to make sure
    ! that the fluxes correctly accounted for the ammonium that is adsorbed and cannot
    ! take part in part of the processes. Here it is added again.
    soil_n_min(:,:,iammonium) = soil_n_min(:,:,iammonium) + n_adsorbed(:,:)

    ! First check whether we can substract the fluxes. If not truncate some of these
    ! fluxes. if they are truncated for the subtraction they should also be truncated
    ! for the addition, else the mass balance will not be closed. The order is arbitrary
    ! but determines which fluxes will get truncated. One key assumption is that all
    ! the calculated fluxes are positive, only then are we sure that these tests make
    ! sense. The code above was adhusted to ensure that this is indeed the case.
    ! Account for nitrification(:,:,i_nh4_to_no3)
    WHERE (soil_n_min(:,:,iammonium)-nitrification(:,:,i_nh4_to_no3) .LT.zero)
       nitrification(:,:,i_nh4_to_no3) = soil_n_min(:,:,iammonium)
       soil_n_min(:,:,iammonium) = zero
    ELSEWHERE
       soil_n_min(:,:,iammonium) = soil_n_min(:,:,iammonium)-nitrification(:,:,i_nh4_to_no3)
    ENDWHERE

    ! Account for nitrification(:,:,i_nh4_to_no)
    WHERE (soil_n_min(:,:,iammonium)-nitrification(:,:,i_nh4_to_no) .LT.zero)
       nitrification(:,:,i_nh4_to_no) = soil_n_min(:,:,iammonium)
       soil_n_min(:,:,iammonium) = zero
    ELSEWHERE
       soil_n_min(:,:,iammonium) = soil_n_min(:,:,iammonium)-nitrification(:,:,i_nh4_to_no)
    ENDWHERE

    ! Account for nitrification(:,:,i_nh4_to_n2o)
    WHERE (soil_n_min(:,:,iammonium)-nitrification(:,:,i_nh4_to_n2o) .LT.zero)
       nitrification(:,:,i_nh4_to_n2o) = soil_n_min(:,:,iammonium)
       soil_n_min(:,:,iammonium) = zero
    ELSEWHERE
       soil_n_min(:,:,iammonium) = soil_n_min(:,:,iammonium)-nitrification(:,:,i_nh4_to_n2o)
    ENDWHERE

    ! Account for leaching(:,:,iammonium)
    WHERE (soil_n_min(:,:,iammonium)-leaching(:,:,iammonium) .LT.zero)
       leaching(:,:,iammonium) = soil_n_min(:,:,iammonium)
       soil_n_min(:,:,iammonium) = zero
    ELSEWHERE
       soil_n_min(:,:,iammonium) = soil_n_min(:,:,iammonium)-leaching(:,:,iammonium)
    ENDWHERE

    ! Account for emission(:,:,iammonium)
    WHERE (soil_n_min(:,:,iammonium)-emission(:,:,iammonium) .LT.zero)
       emission(:,:,iammonium) = soil_n_min(:,:,iammonium)
       soil_n_min(:,:,iammonium) = zero
    ELSEWHERE
       soil_n_min(:,:,iammonium) = soil_n_min(:,:,iammonium)-emission(:,:,iammonium)
    ENDWHERE

    ! Account for plant_n_uptake(:,:,iammonium)
    WHERE (soil_n_min(:,:,iammonium)-plant_n_uptake(:,:,iammonium) .LT.zero)
       plant_n_uptake(:,:,iammonium) = soil_n_min(:,:,iammonium)
       soil_n_min(:,:,iammonium) = zero
    ELSEWHERE
       soil_n_min(:,:,iammonium) = soil_n_min(:,:,iammonium)-plant_n_uptake(:,:,iammonium)
    ENDWHERE

    ! Nitrogen species: initrate
    soil_n_min(:,:,initrate) = soil_n_min(:,:,initrate) + nitrification(:,:,i_nh4_to_no3)

    ! Account for denitrification(:,:,i_no3_to_nox)
    WHERE (soil_n_min(:,:,initrate)-denitrification(:,:,i_no3_to_nox) .LT.zero)
       denitrification(:,:,i_no3_to_nox) = soil_n_min(:,:,initrate)
       soil_n_min(:,:,initrate) = zero
    ELSEWHERE
       soil_n_min(:,:,initrate) = soil_n_min(:,:,initrate)-denitrification(:,:,i_no3_to_nox)
    ENDWHERE

    ! Account for leaching(:,:,initrate)
    WHERE (soil_n_min(:,:,initrate)-leaching(:,:,initrate) .LT.zero)
       leaching(:,:,initrate) = soil_n_min(:,:,initrate)
       soil_n_min(:,:,initrate) = zero
    ELSEWHERE
       soil_n_min(:,:,initrate) = soil_n_min(:,:,initrate)-leaching(:,:,initrate)
    ENDWHERE

    ! Account for plant_n_uptake(:,:,initrate)
    WHERE (soil_n_min(:,:,initrate)-plant_n_uptake(:,:,initrate) .LT.zero)
       plant_n_uptake(:,:,initrate) = soil_n_min(:,:,initrate)
       soil_n_min(:,:,initrate) = zero
    ELSEWHERE
       soil_n_min(:,:,initrate) = soil_n_min(:,:,initrate)-plant_n_uptake(:,:,initrate)
    ENDWHERE

    ! Nitrogen species: inox
    soil_n_min(:,:,inox) =  soil_n_min(:,:,inox) + nitrification(:,:,i_nh4_to_no) & 
         + denitrification(:,:,i_no3_to_nox) 

    ! Account for denitrification(:,:,i_nox_to_n2o)
    WHERE (soil_n_min(:,:,inox)-denitrification(:,:,i_nox_to_n2o) .LT.zero)
       denitrification(:,:,i_nox_to_n2o) = soil_n_min(:,:,inox)
       soil_n_min(:,:,inox) = zero
    ELSEWHERE
       soil_n_min(:,:,inox) = soil_n_min(:,:,inox)-denitrification(:,:,i_nox_to_n2o)
    ENDWHERE

    ! Account for emission(:,:,inox)
    WHERE (soil_n_min(:,:,inox)-emission(:,:,inox) .LT.zero)
       emission(:,:,inox) = soil_n_min(:,:,inox)
       soil_n_min(:,:,inox) = zero
    ELSEWHERE
       soil_n_min(:,:,inox) = soil_n_min(:,:,inox)-emission(:,:,inox)
    ENDWHERE

    ! Nitrogen species: initrous
    soil_n_min(:,:,initrous) = soil_n_min(:,:,initrous) + nitrification(:,:,i_nh4_to_n2o) & 
         + denitrification(:,:,i_nox_to_n2o)

    ! Account for denitrification(:,:,i_n2o_to_n2)
    WHERE (soil_n_min(:,:,initrous)-denitrification(:,:,i_n2o_to_n2) .LT.zero)
       denitrification(:,:,i_n2o_to_n2) = soil_n_min(:,:,initrous)
       soil_n_min(:,:,initrous) = zero
    ELSEWHERE
       soil_n_min(:,:,initrous) = soil_n_min(:,:,initrous)-denitrification(:,:,i_n2o_to_n2)
    ENDWHERE

    ! Account for emission(:,:,initrous)
    WHERE (soil_n_min(:,:,initrous)-emission(:,:,initrous) .LT.zero)
       emission(:,:,initrous) = soil_n_min(:,:,initrous)
       soil_n_min(:,:,initrous) = zero
    ELSEWHERE
       soil_n_min(:,:,initrous) = soil_n_min(:,:,initrous)-emission(:,:,initrous)
    ENDWHERE
   
    ! Nitrogen species: idinitro
    soil_n_min(:,:,idinitro) = soil_n_min(:,:,idinitro)  + denitrification(:,:,i_n2o_to_n2)

    ! Account for emission(:,:,idinitro)
    WHERE (soil_n_min(:,:,idinitro)-emission(:,:,idinitro) .LT.zero)
       emission(:,:,idinitro) = soil_n_min(:,:,idinitro)
       soil_n_min(:,:,idinitro) = zero
    ELSEWHERE
       soil_n_min(:,:,idinitro) = soil_n_min(:,:,idinitro)-emission(:,:,idinitro)
    ENDWHERE

    ! 9.2 add nitrogen either from deposition (read from fields or prescribed in run.def) as
    ! well as BNF, needs to be revised...
    ! iatm_ammo=1, iatm_nitr=2, ibnf=3, imanure=4, ifert_ammo=5, ifert_nitr=6
    DO ipts = 1,npts
       DO ivm = 1,nvm

          ! Deposition of NHx and NOy
          IF(veget_max(ipts,ivm).GT.min_stomate .AND. &
               som(ipts,iactive,ivm,icarbon).GT.min_stomate) THEN
   
             ! If there is an active carbon pool, the n-deposition
             ! is added to the soil pool
             soil_n_min(ipts,ivm,iammonium) = soil_n_min(ipts,ivm,iammonium) + &
                  n_input(ipts,ivm,month,iatm_ammo) * dt
             soil_n_min(ipts,ivm,initrate) = soil_n_min(ipts,ivm,initrate) + &
                  n_input(ipts,ivm,month,iatm_nitr) * dt
   
          ELSEIF(veget_max(ipts,ivm).GT.min_stomate .AND. & 
               som(ipts,iactive,ivm,icarbon).LE.min_stomate) THEN

             ! If there is no active soil carbon pool, the N-deposition
             ! is leached.
             leaching(ipts,ivm,iammonium)=leaching(ipts,ivm,iammonium) + &
                  n_input(ipts,ivm,month,iatm_ammo) * dt
             leaching(ipts,ivm,initrate)=leaching(ipts,ivm,initrate) + &
                 n_input(ipts,ivm,month,iatm_nitr) * dt
   
          ELSEIF(veget_max(ipts,ivm).LT.min_stomate) THEN

             soil_n_min(ipts,ivm,iammonium) = zero
             soil_n_min(ipts,ivm,initrate) = zero
             n_input(ipts,ivm,month,iatm_nitr) = zero
             n_input(ipts,ivm,month,iatm_ammo) = zero

          ELSE

             WRITE(numout,*)'at ipts veget_max for ivm: ', ipts, ivm, veget_max(ipts,ivm)
             WRITE(numout,*)'som', som(ipts,iactive,ivm,icarbon)
             CALL ipslerr_p (3,'nitrogen_dynamics', 'Unforeseen case. for som','','')    
  
          ENDIF 

          ! BNF
          IF(veget_max(ipts,ivm).GT.min_stomate) THEN
                
             IF ( natural(ivm) ) THEN

                IF (ivm.EQ.ibare_sechiba) THEN

                   !  The deposition file contains the same values for all PFTs. If there is
                   !  no veget_max, the n_input should be set to zero. Likewise for the
                   !  natural PFT BNF and fertilization should not be accounted for.
                   n_input(ipts,ivm,month,ibnf) = zero

                ELSE

                   ! Consider biological N-fixation in natural PFTs. Note
                   ! that bare soil is considered a natural PFT but that due
                   ! to its construct bare soil will be treate differently
                   soil_n_min(ipts,ivm,iammonium) = soil_n_min(ipts,ivm,iammonium) + &
                        n_input(ipts,ivm,month,ibnf) * dt
                   
                   ! Fertiliser use for agriculture, no BNF, that's already 
                   ! accounted for in fertil using the average global ratio of 
                   ! ammonium to nitrate to separate the two species
                   ! ratio_nh4_fert = 7./8.
                   soil_n_min(ipts,ivm,iammonium) = soil_n_min(ipts,ivm,iammonium) + &
                        n_input(ipts,ivm,month,ifert_ammo)*dt
                   
                   IF (soil_n_min(ipts,ivm,initrate) + &
                        n_input(ipts,ivm,month,ifert_nitr)*dt .LE. 0.0) THEN

                      !IF (printlev_loc >=4)  WRITE(numout,*) 'untruncated 6, ',soil_n_min(ipts,ivm,initrate) + & 
                      !  n_input(ipts,ivm,month,ifert_nitr)*dt

                      soil_n_min(ipts,ivm,initrate) = MAX(zero,soil_n_min(ipts,ivm,initrate) + & 
                           n_input(ipts,ivm,month,ifert_nitr)*dt)
                   
                   ELSE
                      soil_n_min(ipts,ivm,initrate) = soil_n_min(ipts,ivm,initrate) +  &
                                 n_input(ipts,ivm,month,ifert_nitr)*dt
                   ENDIF

                ENDIF
   
             ELSE
   
                ! Fertiliser use for agriculture, no BNF, that's already 
                ! accounted for in fertil using the average global ratio of 
                ! ammonium to nitrate to separate the two species
                ! ratio_nh4_fert = 7./8.
                soil_n_min(ipts,ivm,iammonium) = soil_n_min(ipts,ivm,iammonium) + &
                     n_input(ipts,ivm,month,ifert_ammo)*dt

                IF (soil_n_min(ipts,ivm,initrate) + &
                     n_input(ipts,ivm,month,ifert_nitr)*dt .LE. 0.0) THEN

                   !IF (printlev_loc >=4)  WRITE(numout,*) 'untruncated 6, ',soil_n_min(ipts,ivm,initrate) + & 
                   !  n_input(ipts,ivm,month,ifert_nitr)*dt

                   soil_n_min(ipts,ivm,initrate) = MAX(zero,soil_n_min(ipts,ivm,initrate) + & 
                     n_input(ipts,ivm,month,ifert_nitr)*dt)

                ELSE
                    soil_n_min(ipts,ivm,initrate) = soil_n_min(ipts,ivm,initrate) +  &
                                 n_input(ipts,ivm,month,ifert_nitr)*dt
                ENDIF

                !  The deposition file contains the same values for all PFTs. If there is
                !  no veget_max, the n_input should be set to zero. Likewise for the
                !  natural PFT fertilization should not be accounted for and for the crops
                !  BNF should not be accounted for
                n_input(ipts,ivm,month,ibnf) = zero

             ENDIF

          ELSE ! veget_max.GT.min_stomate

             !  The deposition file contains the same values for all PFTs. If there is
             !  no veget_max, the n_input should be set to zero. Likewise for the
             !  natural PFT fertilization should not be accounted for and for the crops
             !  BNF should not be accounted for
             n_input(ipts,ivm,month,:) = zero

          ENDIF

       ENDDO
    ENDDO

    ! Error checking
    IF (err_act.GT.1) THEN

       !! 6. Mass balance closure
       !! 6.1 Calculate components of the mass balance
       !  Note that soilcarbon is transfered to other pools but that the pool
       !  itself was not updated. We can just use ::pool_end
       pool_end(:,:,:) = zero
       DO inspec = 1,nnspec
          pool_end(:,:,initrogen) = pool_end(:,:,initrogen) + &
               soil_n_min(:,:,inspec) * veget_max(:,:)
       ENDDO

       
       !! Calculate mass balance
       check_intern(:,:,iatm2land,initrogen) = &
            check_intern(:,:,iatm2land,initrogen) + &
            atm_to_immob(:,:) * veget_max(:,:) * dt
       
       DO ininput = 1,ninput
          IF (ininput .NE. imanure) THEN
            check_intern(:,:,iatm2land,initrogen) = &
                 check_intern(:,:,iatm2land,initrogen) + &
                 n_input(:,:,month,ininput) * dt * veget_max(:,:)
           END IF
       ENDDO
      
       DO inspec= 1, nnspec
          check_intern(:,:,iland2atm,initrogen) = &
               check_intern(:,:,iland2atm,initrogen) &
               -un * (emission(:,:,inspec) * veget_max(:,:)) 
       ENDDO

       DO iion = 1, nionspec
          check_intern(:,:,ilat2out,initrogen) = &
               check_intern(:,:,ilat2out,initrogen) &
               -un * ( plant_n_uptake(:,:,iion) + &
               leaching(:,:,iion) ) * veget_max(:,:)
       ENDDO

       check_intern(:,:,ilat2in,initrogen) = &
            check_intern(:,:,ilat2in,initrogen) + &
            n_mineralisation(:,:) * veget_max(:,:)

       check_intern(:,:,ipoolchange,initrogen) = -un * &
            (pool_end(:,:,initrogen) - pool_start(:,:,initrogen))

       closure_intern = zero
       DO imbc = 1,nmbcomp
          closure_intern(:,:,initrogen) = closure_intern(:,:,initrogen) + &
               check_intern(:,:,imbc,initrogen)
       ENDDO

       CALL check_mass_balance("nitrogen_dynamics", closure_intern, &
          npts, pool_end, pool_start, veget_max,'pft')

    ENDIF ! (ERR_ACT.GT.1) 


    ! 10.  Write output values
    CALL histwrite_p (hist_id_stomate, 'N_UPTAKE_NH4', itime, &
         plant_n_uptake(:,:,iammonium)/dt, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'N_UPTAKE_NO3', itime, &
         plant_n_uptake(:,:,initrate)/dt, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'N_MINERALISATION', itime, &
         n_mineralisation(:,:)/dt, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'SOIL_NH4', itime, &
         soil_n_min(:,:,iammonium), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'SOIL_NO3', itime, &
         soil_n_min(:,:,initrate), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'SOIL_NOX', itime, &
         soil_n_min(:,:,inox), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'SOIL_N2O', itime, &
         soil_n_min(:,:,initrous), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'SOIL_N2', itime, &
         soil_n_min(:,:,idinitro), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'SOIL_P_OX', itime, &
         p_O2(:,:), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'BACT', itime, &
         bact(:,:), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'NH3_EMISSION', itime, &
         emission(:,:,iammonium)/dt, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'NOX_EMISSION', itime, &
         emission(:,:,inox)/dt, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'N2O_EMISSION', itime, &
         emission(:,:,initrous)/dt, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'N2_EMISSION', itime, &
         emission(:,:,idinitro)/dt, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'NH4_LEACHING', itime, &
         leaching(:,:,iammonium)/dt, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'NO3_LEACHING', itime, &
         leaching(:,:,initrate)/dt, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'NITRIFICATION', itime, &
         nitrification(:,:,i_nh4_to_no3), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'DENITRIFICATION', itime, &
         denitrification(:,:,i_n2o_to_n2), npts*nvm, horipft_index)  
    CALL histwrite_p (hist_id_stomate, 'NHX_DEPOSITION', itime, &
         n_input(:,nvm,month,iatm_ammo), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'NOX_DEPOSITION', itime, &
         n_input(:,nvm,month,iatm_nitr), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'BNF', itime, &
         n_input(:,nvm,month,ibnf), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'N_FERTILISER_AMMO', itime, &
         n_input(:,nvm,month,ifert_ammo), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'N_FERTILISER_NITR', itime, &
         n_input(:,nvm,month,ifert_nitr), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'N_MANURE', itime, &
         n_input(:,:,month,imanure), npts*nvm, horipft_index)  
 
    CALL xios_orchidee_send_field("f_NCPLANT",f_NCplant(:,:)) 
    CALL xios_orchidee_send_field("N_UPTAKE_NH4",plant_n_uptake(:,:,iammonium)/dt)
    CALL xios_orchidee_send_field("N_UPTAKE_NO3",plant_n_uptake(:,:,initrate)/dt)
    CALL xios_orchidee_send_field("N_MINERALISATION",n_mineralisation(:,:)/dt)
    CALL xios_orchidee_send_field("SOIL_NH4",soil_n_min(:,:,iammonium))
    CALL xios_orchidee_send_field("SOIL_NO3",soil_n_min(:,:,initrate))
    CALL xios_orchidee_send_field("SOIL_NOX",soil_n_min(:,:,inox))
    CALL xios_orchidee_send_field("SOIL_N2O",soil_n_min(:,:,initrous))
    CALL xios_orchidee_send_field("SOIL_N2",soil_n_min(:,:,idinitro))
    CALL xios_orchidee_send_field("SOIL_P_OX",p_O2(:,:))
    CALL xios_orchidee_send_field("BACT",bact(:,:))
    CALL xios_orchidee_send_field("NH3_EMISSION",emission(:,:,iammonium)/dt)
    CALL xios_orchidee_send_field("NOX_EMISSION",emission(:,:,inox)/dt)
    CALL xios_orchidee_send_field("N2O_EMISSION",emission(:,:,initrous)/dt)
    CALL xios_orchidee_send_field("N2_EMISSION",emission(:,:,idinitro)/dt)
    CALL xios_orchidee_send_field("NH4_LEACHING",leaching(:,:,iammonium)/dt)
    CALL xios_orchidee_send_field("NO3_LEACHING",leaching(:,:,initrate)/dt)
    CALL xios_orchidee_send_field("NITRIFICATION",nitrification(:,:,i_nh4_to_no3))
    CALL xios_orchidee_send_field("DENITRIFICATION",denitrification(:,:,i_n2o_to_n2))
    CALL xios_orchidee_send_field("NHX_DEPOSITION",n_input(:,:,month,iatm_ammo))
    CALL xios_orchidee_send_field("NOX_DEPOSITION",n_input(:,:,month,iatm_nitr))
    CALL xios_orchidee_send_field("BNF",n_input(:,:,month,ibnf))
    CALL xios_orchidee_send_field("N_FERTILISER_AMMO",n_input(:,:,month,ifert_ammo))  
    CALL xios_orchidee_send_field("N_FERTILISER_NITR",n_input(:,:,month,ifert_nitr))  
    CALL xios_orchidee_send_field("N_MANURE",n_input(:,:,month,imanure))

    CALL xios_orchidee_send_field("fBNF",SUM(n_input(:,:,month,ibnf)*veget_max,dim=2)/1e3/one_day)

    CALL xios_orchidee_send_field("fNdep",SUM((n_input(:,:,month,iatm_ammo)+n_input(:,:,month,iatm_nitr))*veget_max,dim=2)/1e3/one_day)
    CALL xios_orchidee_send_field("fNfert",SUM((n_input(:,:,month,ifert_ammo)+n_input(:,:,month,ifert_nitr))*veget_max,dim=2)/1e3/one_day)
    CALL xios_orchidee_send_field("fNgas",SUM((emission(:,:,iammonium)+emission(:,:,inox)+emission(:,:,initrous)+emission(:,:,idinitro))*veget_max,dim=2)/dt/1e3/one_day)
    CALL xios_orchidee_send_field("fN2O",SUM(emission(:,:,initrous)*veget_max,dim=2)/dt/1e3/one_day)
    CALL xios_orchidee_send_field("fNOx",SUM(emission(:,:,inox)*veget_max,dim=2)/dt/1e3/one_day)
    CALL xios_orchidee_send_field("fNleach",SUM((leaching(:,:,iammonium)+leaching(:,:,initrate))*veget_max,dim=2)/dt/1e3/one_day)
    CALL xios_orchidee_send_field("fNnetmin",SUM((n_mineralisation(:,:)+immob(:,:))*veget_max,dim=2)/dt/1e3/one_day)
    CALL xios_orchidee_send_field("fNup",SUM((plant_n_uptake(:,:,iammonium)+plant_n_uptake(:,:,initrate))*veget_max,dim=2)/dt/1e3/one_day)
    CALL xios_orchidee_send_field("nMineral",SUM((soil_n_min(:,:,iammonium)+soil_n_min(:,:,initrate)+soil_n_min(:,:,inox)+&
         soil_n_min(:,:,initrous)+soil_n_min(:,:,idinitro))*veget_max,dim=2)/1e3)
    CALL xios_orchidee_send_field("nMineralNH4",SUM(soil_n_min(:,:,iammonium)*veget_max,dim=2)/1e3)
    CALL xios_orchidee_send_field("nMineralNO3",SUM(soil_n_min(:,:,initrate)*veget_max,dim=2)/1e3)


  END SUBROUTINE nitrogen_dynamics



!! ================================================================================================================================
!!  FUNCTION   : control_temp_func
!!
!>\BRIEF        : Unclear.
!!
!! DESCRIPTION	: Unable to find where this comes from.
!!   Referenced by Zaehle and Friend (2010) in the Appendix,
!!   which points to Krinner et al (2005), but the closest I found
!!   there was Eq. A32, which is simillar but does not
!!   mention Q10 at all.
!!
!! RECENT CHANGE(S): 
!! 
!! MAIN OUTPUTS VARIABLE(S): 
!!
!! REFERENCE(S)   :
!! - S. Zaehle and A. D. Friend (2010), Carbon and nitrogen cycle dynamics in the
!!   O-CN land surface model: 1. Model description, site-scale evaluation, and 
!!   sensitivity to parameter estimates. Global Biogeochem. Cycles, 24, GB1005, 
!!   doi:10.1029/2009GB003521.
!! - Krinner G., N. Viovy, N. de Noblet-Ducoudre, J. Ogee, J. Polcher, P. Friedlingstein, 
!!   P. Ciais, S. Sitch, and I. C. Prentice (2005), A dynamic global vegetation model 
!!   for studies of the coupled atmosphere-biosphere system, Global Biogeochemical 
!!   Cycles, 19, doi.:10.1029/2003/GB002199.
!! 
!! FLOWCHART    :
!!
!_ ================================================================================================================================
  FUNCTION control_temp_func (npts, temp_in) RESULT (tempfunc_result)

  !! 0. Variable and parameter declaration
    
    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                 :: npts            !! Domain size - number of land pixels (unitless)
    REAL(r_std), DIMENSION(:), INTENT(in)   :: temp_in         !! Temperature (K)

    !! 0.2 Output variables
    REAL(r_std), DIMENSION(npts)               :: tempfunc_result !! Temperature control factor (0-1, unitless)

    !! 0.3 Modified variables

    !! 0.4 Local variables

!_ ================================================================================================================================

    tempfunc_result(:) = exp( soil_Q10_uptake * ( temp_in(:) - (ZeroCelsius+tsoil_ref)) / Q10 )
    tempfunc_result(:) = MIN( un, tempfunc_result(:) )

  END FUNCTION control_temp_func

  ! The following subroutines soilcarbon_leak, altcalc_DOC, cryoturbate_doc_POC are not included in the trunk.

END MODULE stomate_som_dynamics
