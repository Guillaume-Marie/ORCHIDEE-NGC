! =================================================================================================================================
! MODULE           : stomate_resp
!
! CONTACT          : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE          : IPSL (2006)
!                  This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF           Calculates maintenance respiration for different plant components
!!
!!\n DESCRIPTION   : None
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S)	   :
!!- McCree KJ. An equation for the respiration of white clover plants grown under controlled conditions. 
!! In: Setlik I, editor. Prediction and measurement of photosynthetic productivity. Wageningen, The Netherlands: 
!! Pudoc; 1970. p. 221-229.
!! - Krinner G, Viovy N, de Noblet-Ducoudre N, Ogee J, Polcher J, Friedlingstein P,
!! Ciais P, Sitch S, Prentice I C (2005) A dynamic global vegetation model for studies
!! of the coupled atmosphere-biosphere system. Global Biogeochemical Cycles, 19, GB1015,
!! doi: 10.1029/2003GB002199.\n
!! Ruimy A., Dedieu G., Saugier B. (1996), TURC: A diagnostic model
!! of continental gross primary productivity and net primary productivity,
!! Global Biogeochemical Cycles, 10, 269-285.\n

!! SVN :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_stomate/stomate_resp.f90 $
!! $Date: 2025-03-24 17:53:05 +0100 (lun. 24 mars 2025) $
!! $Revision: 8879 $
!! \n
!_ ================================================================================================================================
 
MODULE stomate_resp

  ! modules used:
  USE stomate_data
  USE pft_parameters
  USE constantes  
  USE constantes_soil
  USE xios_orchidee

  IMPLICIT NONE

  ! private & public routines
  PRIVATE
  PUBLIC maint_respiration,maint_respiration_clear

  LOGICAL, SAVE                                              :: firstcall_resp = .TRUE.                 !! first call
!$OMP THREADPRIVATE(firstcall_resp)

  CONTAINS


!! ================================================================================================================================
!! SUBROUTINE 	: maint_respiration_clear
!!
!>\BRIEF        : Set the flag ::firstcall_resp to .TRUE. and as such activate section 
!!                1.1 of the subroutine maint_respiration (see below).
!_ ================================================================================================================================

  SUBROUTINE maint_respiration_clear
    firstcall_resp=.TRUE.
  END SUBROUTINE maint_respiration_clear


!! ================================================================================================================================
!! SUBROUTINE 	: maint_respiration
!!
!>\BRIEF         Calculate PFT maintenance respiration of each living plant part by 
!! multiplying the biomass of plant part by maintenance respiration coefficient which
!! depends on long term mean annual temperature. PFT maintenance respiration is carbon flux 
!! with the units @tex $(gC.m^{-2}dt_sechiba^{-1})$ @endtex, and the convention is from plants to the 
!! atmosphere.
!!
!! DESCRIPTION : The maintenance respiration of each plant part for each PFT is the biomass of the plant 
!! part multiplied by maintenance respiration coefficient. The biomass allocation to different 
!! plant parts is done in routine stomate_alloc.f90. The maintenance respiration coefficient is 
!! calculated in this routine.\n
!!
!! The maintenance respiration coefficient is the fraction of biomass that is lost during 
!! each time step, which increases linearly with temperature (2-meter air temperature for aboveground plant
!! tissues; root-zone temperature for below-ground tissues). Air temperature is an input forcing variable. 
!! Root-zone temperature is a convolution of root and soil temperature profiles and also calculated 
!! in this routine.\n
!!
!! The calculation of maintenance respiration coefficient (fraction of biomass respired) depends linearly
!! on temperature:
!! - the relevant temperature for different plant parts (air temperature or root-zone temperature)\n
!! - intercept: prescribed maintenance respiration coefficients at 0 Degree Celsius for 
!!   different plant parts for each PFT in routine stomate_constants.f90\n
!! - slope: calculated with a quadratic polynomial with the multi-annual mean air temperature 
!! (the constants are in routine stomate_constants.f90) as follows\n 
!!    \latexonly
!!      \input{resp3.tex} 
!!    \endlatexonly
!!   Where, maint_resp_slope1, maint_resp_slope2, maint_resp_slope3 are constant in stomate_constants.f90.
!!   Then coeff_maint is calculated as follows:\n
!!    \latexonly
!!      \input{resp4.tex} 
!!    \endlatexonly  
!! If the calculation result is negative, maintenance respiration coefficient will take the value 0.
!! Therefore the maintenance respiration will also be 0.\n
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): PFT maintenance respiration of different plant parts (::resp_maint_part_radia)
!!
!! REFERENCE(S)	:
!! McCree KJ. An equation for the respiration of white clover plants grown under controlled conditions. In: 
!! Setlik I, editor. Prediction and measurement of photosynthetic productivity. Wageningen, 
!! The Netherlands: Pudoc; 1970. p. 221-229.
!! Krinner G, Viovy N, de Noblet-Ducoudre N, Ogee J, Polcher J, Friedlingstein P,
!! Ciais P, Sitch S, Prentice I C (2005) A dynamic global vegetation model for studies
!! of the coupled atmosphere-biosphere system. Global Biogeochemical Cycles, 19, GB1015,
!! doi: 10.1029/2003GB002199.\n
!! Ruimy A., Dedieu G., Saugier B. (1996), TURC: A diagnostic model
!! of continental gross primary productivity and net primary productivity,
!! Global Biogeochemical Cycles, 10, 269-285.\n
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE maint_respiration ( npts, t2m, t2m_longterm, stempdiag, &
       root_profile, circ_class_n, circ_class_biomass,resp_maint_part_radia, cn_leaf_init_2D)

!! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                         :: npts                !! Domain size - number of grid cells (unitless)
    REAL(r_std), DIMENSION(:), INTENT(in)              :: t2m                 !! 2 meter air temperature - forcing variable (K)
    REAL(r_std), DIMENSION(:), INTENT(in)              :: t2m_longterm        !! Long term annual mean 2 meter reference air temperatures 
                                                                              !! calculated in stomate_season.f90 (K)
    REAL(r_std), DIMENSION(:,:), INTENT (in)           :: stempdiag           !! Soil temperature of each soil layer (K)
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)        :: root_profile        !! Normalized root mass/length fraction in each soil layer 
                                                                              !! (0-1, unitless)
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)          :: circ_class_n        !! Number of individuals in each circ class
                                                                              !! @tex $(ind m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(in)      :: circ_class_biomass  !! Biomass components of the model tree  
                                                                              !! within a circumference class
                                                                              !! class @tex $(g C ind^{-1})$ @endtex
    REAL(r_std),DIMENSION(:,:), INTENT(in)             :: cn_leaf_init_2D     !! initial leaf C/N ratio

    !! 0.2 Output variables

    REAL(r_std), DIMENSION(npts,nvm,nparts),INTENT(out):: resp_maint_part_radia !! PFT maintenance respiration of different
                                                                              !! plant parts @tex $(gC.m^{-2}dt^{-1} )$ @endtex

    !! 0.3 Modified variables
   
 
    !! 0.4 Local variables

    INTEGER(i_std)                                    :: ipts,ivm,ipar,islm   !! Indeces (unitless)
    REAL(r_std), DIMENSION(npts,2:nvm)                :: t_root               !! PFT root temperature (convolution of root and soil 
                                                                              !! temperature profiles) (K)
    REAL(r_std)                                       :: coeff_maint          !! PFT maintenance respiration coefficients of different 
                                                                              !! plant compartments at 0 deg C 
    REAL(r_std), DIMENSION(npts,nparts)               :: t_maint_radia        !! Temperature which is pertinent for maintenance respiration, 
                                                                              !! which is air/root temperature for above/below-ground 
                                                                              !! compartments (K)     
    REAL(r_std), DIMENSION(npts,nvm,nparts)           :: gtemp                !! Temperature response of respiration in the
                                                                              !! Lloyd-Taylor Model (-)
    REAL(r_std), DIMENSION(npts,nvm,nparts)           :: cn                   !! CN ratio of a biomass pool ((gC)(gN)-1) 
    REAL(r_std)                                       :: ref_cn               !! Prescribed reference C/N ratio
    REAL(r_std), DIMENSION(npts,nvm,nparts)           :: adjust_resp          !! C/N-based modulator of respiration 
    REAL(r_std)                                       :: tl                   !! Long term reference temperature in degrees Celcius 
                                                                              !! (= t2m_longterm - 273.15) (C)
    REAL(r_std), DIMENSION(npts,nvm,nparts)           :: slope                !! slope of the temperature dependence of maintenance 
                                                                              !! respiration coefficient (1/K)
    REAL(r_std), DIMENSION(npts,nvm,nparts)           :: temp                 !! temporary variable to write to XIOS
    REAL(r_std), DIMENSION(npts,nvm,nparts)           :: temp2                !! temporary variable to write to XIOS
    REAL(r_std)                                       :: err_resp_control    !! 
    REAL(r_std)                                       :: tmps_1
!_ ================================================================================================================================
        
    IF (printlev>=3) WRITE(numout,*) 'Entering maintenance respiration'
    err_resp_control = zero
 !! 1. Initializations
    !! 1.2. Calculate root temperature
    !  Calculate root temperature as the convolution of root and soil temperature profiles
    DO ivm = 2,nvm
       DO ipts = 1,npts

          ! Calculate root temperature
          ! Use the root profile temperature (K) to weight the soil layers 
          ! (number of layers = nslm) at different depths. Note that the vertical axis
          ! of root_profile and stempdiag are centered around the nodes. For root_profile
          ! the center of each layer us given by znh (see vertical_soil.f90). The top
          ! and bottom of the layer are calculated in hydrol_root_profile. If the
          ! naming is correct, the discretisation of the variable stempdiag should follow
          ! diaglev. diaglev is defined in control.f90 making use of znt (see 
          ! vertical_soil.f90) 

          tmps_1 = zero
          DO islm = 1, nslm ! Loop over # soil layers
             tmps_1 = tmps_1 + stempdiag(ipts,islm) * root_profile(ipts,ivm,islm,istruc)
          ENDDO
          t_root(ipts,ivm) = tmps_1

       ENDDO 
    ENDDO 

    ! Initialise
    slope(:,:,:) = zero
    gtemp(:,:,:) = zero
    cn(:,:,:) = zero
    adjust_resp(:,:,:) = zero
    resp_maint_part_radia(:,:,:) = zero
    temp2(:,:,:) = zero

 !! 2. Define maintenance respiration coefficients
    DO ipts = 1,npts

       !! 2.1 Temperature for maintenanace respiration
       !  Temperature which is used to calculate maintenance respiration for different 
       !  plant compartments (above- and belowground):
       !  - for aboveground parts, we use 2-meter air temperature, t2m
       !  - for belowground parts, we use root temperature calculated in section 1.2 of 
       !    this subroutine

       ! 2.1.1 Aboveground biomass

       t_maint_radia(ipts,ileaf) = t2m(ipts)
       t_maint_radia(ipts,isapabove) = t2m(ipts)
       t_maint_radia(ipts,ifruit) = t2m(ipts)
       t_maint_radia(ipts,ilabile) = t2m(ipts)
       t_maint_radia(ipts,iheartabove) = t2m(ipts)

       DO ivm = 2,nvm ! Loop over # PFTs

          ! 2.1.2 Belowground biomass
          t_maint_radia(ipts,isapbelow) = t_root(ipts,ivm)
          t_maint_radia(ipts,iroot) = t_root(ipts,ivm)
          t_maint_radia(ipts,iheartbelow) = t_root(ipts,ivm)

          ! 2.1.3 Depending on the PFT
          IF ( is_tree(ivm) ) THEN
             t_maint_radia(ipts,icarbres) = t2m(ipts)
          ELSE
             t_maint_radia(ipts,icarbres) = t_root(ipts,ivm)
          ENDIF

          !! 2.2 Calculate maitenance respiration coefficients (coeff_maint)
          ! The calculation of the maintenance repiration has been a topic of long and unresolved
          ! debate. This is reflected in the different approaches that can be found in the CMIP5 
          ! trunk, ORCHIDEE-CAN, ORCHIDEE-CNP, O-CN, and the CMIP6 trunk. The approach for CMIP5
          ! was inspired on Krinner et al 2005. The later approaches were inspired on Sitch et al 
          ! 2003 and the equations were consistent with that paper. The parameter setting for 
          ! coeff_maint is in the range of 0.066 to 0.011 as reported in the paper but exact values 
          ! are not given. Although, the principle of a climate correction for coeff_maint is
          ! mentioned in Sitch et al 2003, the reduction factors themselves were not given. As it 
          ! appears now this block of code pretends much more knowledge then we actually have. Rather 
          ! than using a baseline coeff_maint that is later corrected for the climate region, the 
          ! parameter values for coeff_maint could be simply prescribed and made pft-specific.

          ! There are however a couple of problems with that approach: 
          ! (1) a PFT specific respiration coefficient is used to 
          ! address the observation that plants that grow in warmer regions repsire for a given 
          ! air temperature less than plants growing in a colder region. In Sitch et al, the PFT
          ! specific maint_coeff thus compensate for the temperature effect calculated by gtemp
          ! (see below). (2) Maintenance respiration increase if the N pool increases resulting in
          ! an apparent decrease in NPP (but the absolute value of GPP and NPP should still increase). 
          ! Following increased N-availability, NPP/GPP in ORCHIDEE does not changes a  lot. 
          ! This is the opposite of what has been observed in over a century of fertilization experiments 
          ! and in more recent meta-analyses such as Vicca et al 2012. Vicca et al 2012 suggest that 
          ! an increase in NPP/GPP following fertilization is due to the fact that the C loss to 
          ! myccorgizae are decreasing. In ORCHIDEE this loss is accounted for the maintenance 
          ! respiration, it is hidden in COEFF_MAINT_RESP, and (3) estimates of Ra_maint also 
          ! include the C-fluxes to leaching, BVOCs mycorrhizae, etc. and will therefore be higher
          ! than the observations. This issue should be addressed by adding mycorrhizae in the 
          ! in the model.

          ! Given Vicca et al 2012 (Ecology Letters) an NPP/GPP ratio of 0.5 is 'universal' for 
          ! forests given a sufficient nutrient supply and strictly defining NPP as solely its 
          ! biomass components (thus excluding VOC, exudation and subsidies to myccorrhizae as is 
          ! the case in ORCHIDEE). Unless these currently missing fluxes are added to ORCHIDEE and
          ! observation based values become available for coeff_maint_init, this parameter will be
          ! adjusted to obtain an NPP/GPP of 0.5 in the absence of nutrient limitations.

          SELECT CASE (maint_resp_control)

             CASE ('nitrogen')

                ! This approach follows the idea that maintenance respiration is driven by the 
                ! nitrogen pools. A PFT specific coeff_maint_init was kept as was the temperature
                ! correction. Following Ali et al 2016 a small correction was made for structural
                ! nitrogen (which is considered to be 4% of the nitrogen in the respiring pools. 
                ! We no longer use the parameters values for coeff_maint_init as presented in 
                ! Sitch et al 2003 as a plausible range. Values have been adjusted to obtain 
                ! reasonable NPP/GPP values which is more important for the rest of the 
                ! simulation.
                DO ipar = 1, nparts
                   IF ( ipar.EQ.ileaf .OR. ipar.EQ.iroot .OR. ipar.EQ.ifruit .OR. &
                           ipar.EQ.isapabove .OR. ipar.EQ.isapbelow) THEN

                      ! Plant parts that respire - the values have been optimized to reproduce 
                      ! observed NPP/GPP ratios
                      coeff_maint = coeff_maint_init(ivm) * dt_sechiba/one_day

                   ELSE

                      ! None respiring plant parts: heartwood, reserve pool 
                      coeff_maint = zero

                   ENDIF

                   !! Calculate maintenance respiration coefficients
                   ! LPJ respiration factors based on Sitch et al. 2003 - second part of the calculation
                   ! Temperature response, LLoyd and Taylor, 1994. E0 = 308.56 comes from the paper of
                   ! Lloyd and Taylor but was fitted for soil respiration which only partly consists of
                   ! authotrophic (root) respiration. With E0 = 308.56 the temperature sensitivity of
                   ! resp_maint is way too high for PFTs that occur along a substantial temperature gradient
                   ! such as the C3 grasslands. Too high means that the simulated NPP/GPP ratio (0.1 to 0.9) 
                   ! exceeds the observed NPP/GPP across PFTs (0.3 to 0.7). There are two ways to tackle this
                   ! issue: (1) reduce the range of a single PFT (this has been done when moving from 13 to 15  
                   ! PFTs) and (2) reduce the temperature sensitivity. The latter is what is being done by 
                   ! introducing the slope_ra parameter. Slope_ra is arbitrary so reducing the spatial extent 
                   ! of a single PFT is probably the more scientific way forward.

                   IF (t_maint_radia(ipts,ipar)-ZeroCelsius-tmin_maint_resp(ivm).GT.min_stomate) THEN
                      gtemp(ipts,ivm,ipar) = EXP((e0_maint_resp(ivm)/slope_ra)*&
                              (1.0/(tref_maint_resp(ivm)-tmin_maint_resp(ivm))-1.0 / &
                              (t_maint_radia(ipts,ipar)-ZeroCelsius-tmin_maint_resp(ivm))))
                   ELSE
                      ! No gtemp below -46.01 degrees Celsius
                      gtemp(ipts,ivm,ipar) = 0.0
                   ENDIF
             
                   ! maint_resp seems very low for deciduous species with this formulation. 
                   ! Following Ali et al. (2016) we account for structural N which does not 
                   ! contribute to respiration
                   resp_maint_part_radia(ipts,ivm,ipar) = coeff_maint * gtemp(ipts,ivm,ipar) * & 
                           (SUM(circ_class_biomass(ipts,ivm,:,ipar,initrogen)*circ_class_n(ipts,ivm,:)) - &
                           snc*SUM(circ_class_biomass(ipts,ivm,:,ipar,icarbon)*circ_class_n(ipts,ivm,:)))

                ENDDO ! Loop over # plant parts

             CASE ('cn')
                
                ! ORCHIDEE 3 approach
                ! This approach refines Sitch et al 2003 by constraining the respiration with C/N ratios. 
                ! The C/N ratios were reset to the values presented in Sitch et al 2003 but still seem on the low side. 
                ! If pft-specific values are to be used, changes in respiration could be compensated for by changing
                ! coeff_maint_init. Given Vicca et al 2012 (Ecology Letters) an NPP/GPP ratio of 0.5 is 'universal'
                ! for forests given a sufficient nutrient supply and strictly defining NPP as solely its biomass
                ! components (thus excluding VOC, exudation and subsidies to myccorrhizae as is the case in ORCHIDEE).
                ! Unless observation based values are available for coeff_maint_init, these values could be adjusted
                ! within the range of 0.066 to 0.011 to obtain an NPP/GPP of 0.5 in the absence of nutrient
                ! limitations.
                DO ipar = 1, nparts ! Loop over # plant parts
                   ! LPJ respiration factors based on Sitch et al. 2003 - first part of the calculation
                   ! in the orc4 case ilabile is excluded but icarbres is included
                   IF ( ipar.EQ.iheartabove .OR. ipar.EQ.iheartbelow .OR. ipar.EQ.ilabile) THEN
                      coeff_maint = zero
                   ELSE
                      ! Use a  PFT-specific value - Values from OCN are used
                      coeff_maint = coeff_maint_init(ivm)*dt_sechiba/one_day
                   ENDIF

                   ! Fall back on Krinner et al 2005 to calculate the temperature
                   ! dependency of each PFT. Note that the calculations of the slope is tricky. In 
                   ! ORCHIDEE the calculation of Ra includes Rm (which is likely to be temperature 
                   ! dependent), Rg (which depends on the growth) as well as C-subsidies to mycorrhizae 
                   ! and leaching. C-subsisdies to mycorrhizae are nutrient-dependent. We are using 
                   ! a strict temperature dependency to simulate Rm and to account for the fact 
                   ! that we don calculate C-subsidies (Ra will be too high but NPP should still be 
                   ! correct in ORCHIDEE). NPP/GPP should be highest in the temperate zone and lowest
                   ! in the boreal and tropics. If we allow tl to become negative (which is the case
                   ! in the high artic) the slope calculations becomes extremely difficult to control.
                   ! Hence, tl is trucated at zero.
                   tl = MAX(zero, t2m_longterm(ipts) - ZeroCelsius)
                   slope(ipts,ivm,ipar) = maint_resp_slope(ivm,1) + tl * maint_resp_slope(ivm,2) + &
                   tl*tl * maint_resp_slope(ivm,3)

                   ! When, in the original formulation, the temperature dropped below zero, 
                   ! gtemp dropped below one but it was prevented that it became negative.
                   ! Such and approached implied that we believed that respiration decreases
                   ! at sub zero temperature until it becomes zero somewhere between -5 and -10 
                   ! depending the value of slope (which depends on the PFT). In other words at 
                   ! sub zero temperatures Rm was believed to stop. In the new formulation
                   ! we think that at zero degrees Rm reaches its minimum but can no longer
                   ! decrease. The difference between the new and the old approach is either
                   ! zero (old) or a one (un; new) in the MAX statements
                   ! in the orc4 case the maximum and the calculation and 1 is taken where 1
                   ! denotes the minimum respiration.
                   gtemp(ipts,ivm,ipar) = MAX(( un + slope(ipts,ivm,ipar) * &
                           (t_maint_radia(ipts,ipar)-ZeroCelsius) ), zero)

                   ! Calculate the actual C/N ratio
                   IF (SUM(circ_class_biomass(ipts,ivm,:,ipar,initrogen)*circ_class_n(ipts,ivm,:)).GT.min_stomate) THEN
                      cn(ipts,ivm,ipar) = SUM(circ_class_biomass(ipts,ivm,:,ipar,icarbon)*circ_class_n(ipts,ivm,:)) / &
                              SUM(circ_class_biomass(ipts,ivm,:,ipar,initrogen)*circ_class_n(ipts,ivm,:))
                   ELSE
                      cn(ipts,ivm,ipar) = zero
                   ENDIF

                   ! Calculate the limiting C/N ratio
                   ! orc3 calculates the variable limit_cn but the variable is no longer used.
                   ! That code has not been copied.
                   IF ( ipar.EQ.ileaf ) THEN
                      ref_cn=45.
                   ELSEIF ( ipar.EQ.iroot .OR. ipar.EQ.ifruit ) THEN
                      ref_cn=45./fcn_root(ivm)
                   ELSEIF ( ipar.EQ.isapabove .OR. ipar.EQ.isapbelow .OR. ipar.EQ.icarbres) THEN
                      ref_cn=45./fcn_wood(ivm)
                   ELSE
                      ref_cn=45.
                   ENDIF

                   ! Use the ref_cn to calculate a reduction factor. Several of the values
                   ! in this equation were tuned for the ORCHIDEE 3.0
                   adjust_resp(ipts,ivm,ipar)=cn(ipts,ivm,ipar)/ref_cn * &
                           MAX(MIN( ( 1. + ( 1. - cn(ipts,ivm,ipar) / ref_cn ) * 3./10. ), 1.2),&
                           0.8)
                   IF(SUM(circ_class_biomass(ipts,ivm,:,ipar,initrogen)*circ_class_n(ipts,ivm,:)) .GT. min_stomate)THEN
                      resp_maint_part_radia(ipts,ivm,ipar) = coeff_maint * gtemp(ipts,ivm,ipar) * &
                           SUM(circ_class_biomass(ipts,ivm,:,ipar,initrogen)*circ_class_n(ipts,ivm,:)) * &
                           adjust_resp(ipts,ivm,ipar)
                   ELSE
                      resp_maint_part_radia(ipts,ivm,ipar) = zero
                   ENDIF

                ENDDO ! Loop over # plant parts

             CASE ('trunk')
       
                ! This approach refines Sitch et al 2003 by constraining the respiration with C/N ratios. 
                ! The C/N ratios were reset to the values presented in Sitch et al 2003 but still seem on the low side. 
                ! If pft-specific values are to be used, changes in respiration could be compensated for by changing
                ! coeff_maint_init. Given Vicca et al 2012 (Ecology Letters) an NPP/GPP ratio of 0.5 is 'universal'
                ! for forests given a sufficient nutrient supply and strictly defining NPP as solely its biomass
                ! components (thus excluding VOC, exudation and subsidies to myccorrhizae as is the case in ORCHIDEE).
                ! Unless observation based values are available for coeff_maint_init, these values could be adjusted
                ! within the range of 0.066 to 0.011 to obtain an NPP/GPP of 0.5 in the absence of nutrient
                ! limitations.
                DO ipar = 1, nparts ! Loop over # plant parts
                   ! LPJ respiration factors based on Sitch et al. 2003 - first part of the calculation
                   IF ( ipar.EQ.iheartabove .OR. ipar.EQ.iheartbelow .OR. ipar.EQ.icarbres) THEN
                      coeff_maint = zero
                   ELSE
                      ! Use a  PFT-specific value - Values from OCN are used
                      coeff_maint = coeff_maint_init(ivm)*dt_sechiba/one_day
                   ENDIF

                   ! Fall back on Krinner et al 2005 to calculate the temperature
                   ! dependency of each PFT. Note that the calculations of the slope is tricky. In 
                   ! ORCHIDEE the calculation of Ra includes Rm (which is likely to be temperature 
                   ! dependent), Rg (which depends on the growth) as well as C-subsidies to mycorrhizae 
                   ! and leaching. C-subsisdies to mycorrhizae are nutrient-dependent. We are using 
                   ! a strict temperature dependency to simulate Rm and to account for the fact 
                   ! that we don calculate C-subsidies (Ra will be too high but NPP should still be 
                   ! correct in ORCHIDEE). NPP/GPP should be highest in the temperate zone and lowest
                   ! in the boreal and tropics. If we allow tl to become negative (which is the case
                   ! in the high artic) the slope calculations becomes extremely difficult to control.
                   ! Hence, tl is trucated at zero.
                   tl = MAX(zero, t2m_longterm(ipts) - ZeroCelsius)
                   slope(ipts,ivm,ipar) = maint_resp_slope(ivm,1) + tl * maint_resp_slope(ivm,2) + &
                   tl*tl * maint_resp_slope(ivm,3)

                   ! When, in the original formulation, the temperature dropped below zero, 
                   ! gtemp dropped below one but it was prevented that it became negative.
                   ! Such and approached implied that we believed that respiration decreases
                   ! at sub zero temperature until it becomes zero somewhere between -5 and -10 
                   ! depending the value of slope (which depends on the PFT). In other words at 
                   ! sub zero temperatures Rm was believed to stop. In the new formulation
                   ! we think that at zero degrees Rm reaches its minimum but can no longer
                   ! decrease. The difference between the new and the old approach is either
                   ! zero (old) or a one (un; new) in the MAX statements
                   gtemp(ipts,ivm,ipar) = MAX(( un + slope(ipts,ivm,ipar) * &
                           (t_maint_radia(ipts,ipar)-ZeroCelsius) ), un)

                   ! Calculate the actual C/N ratio
                   IF (SUM(circ_class_biomass(ipts,ivm,:,ipar,initrogen)*circ_class_n(ipts,ivm,:)).GT.min_stomate) THEN
                      cn(ipts,ivm,ipar) = SUM(circ_class_biomass(ipts,ivm,:,ipar,icarbon)*circ_class_n(ipts,ivm,:)) / &
                              SUM(circ_class_biomass(ipts,ivm,:,ipar,initrogen)*circ_class_n(ipts,ivm,:))
                   ELSE
                      cn(ipts,ivm,ipar) = zero
                   ENDIF

                   ! The model does not control the C/N ratio of the labile and carbohydrate 
                   ! pools. The C/N ratio is truncated here. This has probably a relative small 
                   ! impact because a very high c/n happens when the nitrogen in this pool is 
                   ! very low. When nitrogen is very low it has little impact on the calculation 
                   ! of resp_maint_part_radia because that is based on the nitrogen pool. The
                   ! threshold of 200 is a bit arbitrairy but it should ensure that it stays
                   ! with the physiological boundaries. 
                   IF (cn(ipts,ivm,ipar) .GT. 200) THEN
                      cn(ipts,ivm,ipar) = 200
                   ENDIF

                   ! Calculate the limiting C/N ratio
                   IF ( ipar.EQ.ileaf ) THEN
                      ref_cn=45.
                   ELSEIF ( ipar.EQ.iroot .OR. ipar.EQ.ifruit ) THEN
                      ref_cn=45./fcn_root(ivm)
                   ELSEIF ( ipar.EQ.isapabove .OR. ipar.EQ.isapbelow .OR. ipar.EQ.icarbres) THEN
                      ref_cn=45./fcn_wood(ivm)
                   ELSE
                      ref_cn=45.
                   ENDIF

                   ! Use the ref_cn to calculate a reduction factor. Several of the values
                   ! in this equation were tuned for the ORCHIDEE 3.0
                   adjust_resp(ipts,ivm,ipar)=cn(ipts,ivm,ipar)/ref_cn * &
                           MAX(MIN( ( 1. + ( 1. - cn(ipts,ivm,ipar) / ref_cn ) * 3./10. ), 1.2),&
                           0.8)
                   IF(SUM(circ_class_biomass(ipts,ivm,:,ipar,initrogen)*circ_class_n(ipts,ivm,:)) .GT. min_stomate)THEN
                      resp_maint_part_radia(ipts,ivm,ipar) = coeff_maint * gtemp(ipts,ivm,ipar) * &
                              SUM(circ_class_biomass(ipts,ivm,:,ipar,initrogen)*circ_class_n(ipts,ivm,:)) * adjust_resp(ipts,ivm,ipar)
                      ! Debug - understanding the impact of the nitrogen pool
                      temp2(ipts,ivm,ipar) = adjust_resp(ipts,ivm,ipar) * &
                              SUM(circ_class_biomass(ipts,ivm,:,ipar,initrogen)*circ_class_n(ipts,ivm,:))
                   ELSE
                      resp_maint_part_radia(ipts,ivm,ipar) = zero
                      ! Debug - understanding the impact of the nitrogen pool 
                      temp2(ipts,ivm,ipar) = zero
                   ENDIF

                ENDDO ! Loop over # plant parts
                
             CASE default
                ! All other cases
                err_resp_control = err_resp_control + un

          END SELECT
       ENDDO ! Loop over # PFTs  
    ENDDO

    IF (err_resp_control.GT.zero) THEN
       WRITE(numout,*) 'MAINT_RESP_CONTROL was set to: ',maint_resp_control
       CALL ipslerr_p (3,'stomate_resp', 'don''t know how to calculate maint_resp', &
               'Check orchidee.def', '')
    ENDIF

    ! Write details for debugging and tuning 
    WHERE (resp_maint_part_radia(:,:,:) .EQ. zero)
       temp(:,:,:) = xios_default_val
    ELSEWHERE
       temp(:,:,:) = resp_maint_part_radia(:,:,:)
    ENDWHERE
    CALL xios_orchidee_send_field("RESP_MAINT_PART",temp)
    
    WHERE (slope(:,:,:) .EQ. zero)
       temp(:,:,:) = xios_default_val
    ELSEWHERE
       temp(:,:,:) = slope(:,:,:)
    ENDWHERE
    CALL xios_orchidee_send_field("SLOPE_MAINT_PART",temp)
    
    WHERE (gtemp(:,:,:) .EQ. zero)
       temp(:,:,:) = xios_default_val
    ELSEWHERE
       temp(:,:,:) = gtemp(:,:,:)
    ENDWHERE
    CALL xios_orchidee_send_field("GTEMP_MAINT_PART",temp)
    
    WHERE (cn(:,:,:) .EQ. zero)
       temp(:,:,:) = xios_default_val
    ELSEWHERE
       temp(:,:,:) = cn(:,:,:)
    ENDWHERE
    CALL xios_orchidee_send_field("CN_MAINT_PART",temp)
    
    WHERE (adjust_resp(:,:,:) .EQ. zero)
       temp(:,:,:) = xios_default_val
    ELSEWHERE
       temp(:,:,:) = adjust_resp(:,:,:)
    ENDWHERE
    CALL xios_orchidee_send_field("ADJUST_MAINT_PART",temp)

    WHERE (temp2(:,:,:) .EQ. zero)
       temp2(:,:,:) = xios_default_val
    ENDWHERE
    CALL xios_orchidee_send_field("NPOOL_MAINT_PART",temp2)

 !! 4. Check consistency of this routine

    !  This routine only calculates respiration factors but respiration itself
    !  is not accounted for through pools and fluxes. Hence, there is no need to
    !  CALL check_veget_max and CALL check_mass_balance 

    IF (printlev>=3) WRITE(numout,*) 'Leaving maintenance respiration'

  END SUBROUTINE maint_respiration

END MODULE stomate_resp
