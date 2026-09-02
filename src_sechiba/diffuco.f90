! ================================================================================================================================
!  MODULE       : diffuco
!
!  CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
!  LICENCE      : IPSL (2006)
!  This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF   This module calculates the limiting coefficients, both aerodynamic
!! and hydrological, for the turbulent heat fluxes.
!!
!!\n DESCRIPTION: The aerodynamic resistance R_a is used to limit
!! the transport of fluxes from the surface layer of vegetation to the point in the atmosphere at which
!! interaction with the LMDZ atmospheric circulation model takes place. The aerodynamic resistance is
!! calculated either within the module r_aerod (if the surface drag coefficient is provided by the LMDZ, and 
!! if the flag 'ldq_cdrag_from_gcm' is set to TRUE) or r_aero (if the surface drag coefficient must be calculated).\n
!!
!! Within ORCHIDEE, evapotranspiration is a function of the Evaporation Potential, but is modulated by a
!! series of resistances (canopy and aerodynamic) of the surface layer, here represented by beta.\n
!!
!! DESCRIPTION	:
!!
!! This module calculates the beta for several different scenarios: 
!! - diffuco_snow calculates the beta coefficient for sublimation by snow, 
!! - diffuco_inter calculates the beta coefficient for interception loss by each type of vegetation, 
!! - diffuco_bare calculates the beta coefficient for bare soil, 
!! - diffuco_trans_co2 calculates the beta coefficient for transpiration for each type of vegetation, using Farqhar's formula
!! - chemistry_bvoc calculates the beta coefficient for emissions of biogenic compounds \n
!!
!! Finally, the module diffuco_comb computes the combined $\alpha$ and $\beta$ coefficients for use 
!! elsewhere in the module. \n

!! RECENT CHANGE(S): Nathalie le 28 mars 2006 - sur proposition de Fred Hourdin, ajout
!! d'un potentiometre pour regler la resistance de la vegetation (rveg is now in pft_parameters)
!! October 2018: Removed diffuco_trans using Jarvis formula for calculation of beta coefficient
!!
!! REFERENCE(S)	: None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_sechiba/diffuco.f90 $
!! $Date: 2026-03-03 14:02:01 +0100 (mar. 03 mars 2026) $
!! $Revision: 9417 $
!! \n
!_ ================================================================================================================================

MODULE diffuco

  ! modules used :
  USE constantes
  USE constantes_soil
  USE qsat_moisture
  USE sechiba_io_p
  USE ioipsl
  USE pft_parameters
  USE grid
  USE time, ONLY : one_day, dt_sechiba, day_end
  USE ioipsl_para 
  USE xios_orchidee
  USE chemistry, ONLY : chemistry_initialize, chemistry_bvoc, chemistry_clear
  USE function_library, ONLY: Arrhenius, Arrhenius_modified, get_printlev
  IMPLICIT NONE

  ! public routines :
  PRIVATE
  PUBLIC :: diffuco_main, diffuco_initialize, diffuco_finalize, diffuco_clear

  !
  ! variables used inside diffuco module : declaration and initialisation
  !
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: wind                      !! Wind module (m s^{-1})
!$OMP THREADPRIVATE(wind)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: test_overPET_intensity    !! Test variable to check how much evap is greater than evapot 
!$OMP THREADPRIVATE(test_overPET_intensity)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: test_overPET_freq         !! Number of occurances for test 
!$OMP THREADPRIVATE(test_overPET_freq)
  INTEGER(i_std), SAVE, PRIVATE                      :: printlev_loc              !! Local level of text output for current module
!$OMP THREADPRIVATE(printlev_loc)


CONTAINS


!!  =============================================================================================================================
!! SUBROUTINE:    diffuco_initialize
!!
!>\BRIEF	  Allocate module variables, read from restart file or initialize with default values
!!
!! DESCRIPTION:	  Allocate module variables, read from restart file or initialize with default values.
!!                Call chemistry_initialize for initialization of variables needed for the calculations of BVOCs.
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S): None
!! 
!! FLOWCHART: None
!! \n
!_ ==============================================================================================================================
  SUBROUTINE diffuco_initialize (kjit,    kjpindex, index,                  &
                                 rest_id, lalo,     neighbours, resolution, &
                                 rstruct, q_cdrag)
    
    !! 0. Variable and parameter declaration
    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                         :: kjit             !! Time step number (-) 
    INTEGER(i_std), INTENT(in)                         :: kjpindex         !! Domain size (-)
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)   :: index            !! Indeces of the points on the map (-)
    INTEGER(i_std),INTENT (in)                         :: rest_id          !! _Restart_ file identifier (-)
    REAL(r_std),DIMENSION (kjpindex,2),   INTENT (in)  :: lalo             !! Geographical coordinates
    INTEGER(i_std),DIMENSION (kjpindex,NbNeighb),INTENT (in):: neighbours  !! Vector of neighbours for each 
    REAL(r_std),DIMENSION (kjpindex,2), INTENT(in)     :: resolution       !! The size in km of each grid-box in X and Y
    
    !! 0.2 Output variables
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: rstruct          !! Structural resistance for the vegetation
    
    !! 0.3 Modified variables
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: q_cdrag          !! Surface drag coefficient  (-)
    
    !! 0.4 Local variables
    INTEGER                                            :: ilai
    INTEGER                                            :: jv
    INTEGER                                            :: ier
    CHARACTER(LEN=4)                                   :: laistring
    CHARACTER(LEN=80)                                  :: var_name        
    !_ ================================================================================================================================
   
    !! Initialize local printlev variable
    printlev_loc=get_printlev('diffuco')
    
 
    !! 1. Define flag ldq_cdrag_from_gcm. This flag determines if the cdrag should be taken from the GCM or be calculated. 
    !!    The default value is true if the q_cdrag variables was already initialized. This is the case when coupled to the LMDZ.

    !Config Key   = CDRAG_FROM_GCM
    !Config Desc  = Keep cdrag coefficient from gcm.
    !Config If    = OK_SECHIBA
    !Config Def   = y
    !Config Help  = Set to .TRUE. if you want q_cdrag coming from GCM (if q_cdrag on initialization is non zero).
    !Config         Keep cdrag coefficient from gcm for latent and sensible heat fluxes.
    !Config Units = [FLAG]
    IF ( ABS(MAXVAL(q_cdrag)) .LE. EPSILON(q_cdrag)) THEN
       ldq_cdrag_from_gcm = .FALSE.
    ELSE
       ldq_cdrag_from_gcm = .TRUE.
    ENDIF
    CALL getin_p('CDRAG_from_GCM', ldq_cdrag_from_gcm)
    IF (printlev_loc>=2) WRITE(numout,*) "ldq_cdrag_from_gcm = ",ldq_cdrag_from_gcm

    !! 2. Allocate module variables
    ALLOCATE (wind(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'diffuco_initialize','Problem in allocate of variable wind','','')

    !! 3. Read variables from restart file
    IF (printlev>=3) WRITE (numout,*) 'Read DIFFUCO variables from restart file'

    CALL ioconf_setatt_p('UNITS', 's/m')
    CALL ioconf_setatt_p('LONG_NAME','Structural resistance')
    CALL restget_p (rest_id, 'rstruct', nbp_glo, nvm, 1, kjit, .TRUE., rstruct, "gather", nbp_glo, index_g)
    IF ( ALL(rstruct(:,:) == val_exp) ) THEN
       DO jv = 1, nvm
          rstruct(:,jv) = rstruct_const(jv)
       ENDDO
    ENDIF
    
    !! 4. Initialize chemistry module
    IF (printlev>=3) WRITE(numout,*) "ok_bvoc:",ok_bvoc
    IF ( ok_bvoc ) CALL chemistry_initialize(kjpindex, lalo, neighbours, resolution)

    ALLOCATE (test_overPET_freq(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'diffuco_initialize','Problem in allocate of variable test_overPET_freq','','')

    ALLOCATE (test_overPET_intensity(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'diffuco_initialize','Problem in allocate of variable test_overPET_intensity','','')

    ! Note that test_overPET_freq is not restartable but it is only a secondary diagnostic varibale
    test_overPET_intensity(:)=xios_default_val
    test_overPET_freq(:)=zero


  END SUBROUTINE diffuco_initialize



!! ================================================================================================================================
!! SUBROUTINE    : diffuco_main
!!
!>\BRIEF	 The root subroutine for the module, which calls all other required
!! subroutines.
!! 
!! DESCRIPTION   : 

!! This is the main subroutine for the module. 
!! First it calculates the surface drag coefficient (via a call to diffuco_aero), using available parameters to determine
!! the stability of air in the surface layer by calculating the Richardson Nubmber. If a value for the 
!! surface drag coefficient is passed down from the atmospheric model and and if the flag 'ldq_cdrag_from_gcm' 
!! is set to TRUE, then the subroutine diffuco_aerod is called instead. This calculates the aerodynamic coefficient. \n
!!
!! Following this, an estimation of the saturated humidity at the surface is made (via a call
!! to qsatcalc in the module qsat_moisture). Following this the beta coefficients for sublimation (via 
!! diffuco_snow), interception (diffuco_inter), bare soil (diffuco_bare), and transpiration (via 
!! diffuco_trans_co2) are calculated in sequence. Finally 
!! the alpha and beta coefficients are combined (diffuco_comb). \n
!!
!! The surface drag coefficient is calculated for use within the module enerbil. It is required to to
!! calculate the aerodynamic coefficient for every flux. \n
!!
!! The various beta coefficients are used within the module enerbil for modifying the rate of evaporation, 
!! as appropriate for the surface. As explained in Chapter 2 of Guimberteau (2010), that module (enerbil) 
!! calculates the rate of evaporation essentially according to the expression $E = /beta E_{pot}$, where
!! E is the total evaporation and $E_{pot}$ the evaporation potential. If $\beta = 1$, there would be
!! essentially no resistance to evaporation, whereas should $\beta = 0$, there would be no evaporation and
!! the surface layer would be subject to some very stong hydrological stress. \n
!!
!! The following processes are calculated:
!! - call diffuco_aero for aerodynamic transfer coeficient
!! - call diffuco_snow for partial beta coefficient: sublimation
!! - call diffuco_inter for partial beta coefficient: interception for each type of vegetation
!! - call diffuco_bare for partial beta coefficient: bare soil
!! - call diffuco_trans_co2 for partial beta coefficient: transpiration for each type of vegetation, using Farqhar's formula
!! - call diffuco_comb for alpha and beta coefficient
!! - call chemistry_bvoc for alpha and beta coefficients for biogenic emissions
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): humrel, q_cdrag, vbeta, vbeta1, vbeta4,
!! vbeta2, vbeta3, rveget, cimean   
!!
!! REFERENCE(S) :				        
!! - de Noblet-Ducoudré, N, Laval, K & Perrier, A, 1993. SECHIBA, a new set of parameterisations
!! of the hydrologic exchanges at the land-atmosphere interface within the LMD Atmospheric General
!! Circulation Model. Journal of Climate, 6, pp.248-273.
!! - de Rosnay, P, 1999. Représentation des interactions sol-plante-atmosphère dans le modèle de circulation générale
!! du LMD, 1999. PhD Thesis, Université Paris 6, available (25/01/12): 
!! http://www.ecmwf.int/staff/patricia_de_rosnay/publications.html#8
!! - Ducharne, A, 1997. Le cycle de l'eau: modélisation de l'hydrologie continentale, étude de ses interactions avec 
!! le climat, PhD Thesis, Université Paris 6
!! - Guimberteau, M, 2010. Modélisation de l'hydrologie continentale et influences de l'irrigation
!! sur le cycle de l'eau, PhD Thesis, available (25/01/12):
!! http://www.sisyphe.upmc.fr/~guimberteau/docs/manuscrit_these.pdf
!! - Lathière, J, 2005. Evolution des émissions de composés organiques et azotés par la biosphère continentale dans le 
!! modèle LMDz-INCA-ORCHIDEE, Université Paris 6
!!
!! FLOWCHART	:
!! \latexonly 
!!     \includegraphics[scale=0.5]{diffuco_main_flowchart.png}
!! \endlatexonly
!! \n
!_ ================================================================================================================================

  SUBROUTINE diffuco_main (kjit, kjpindex, index, indexveg, indexlai, u, v, &
     & zlev, z0m, z0h, roughheight, temp_sol, temp_air, temp_growth, rau, q_cdrag, qsurf, qair, pb, &
     & evap_bare_lim, evap_bare_lim_ns, evapot, evapot_corr, snow, flood_frac, flood_res, &
     & frac_nobio, totfrac_nobio, &
     & swnet, swdown, coszang, ccanopy, humrel, veget, veget_max, lai, &
     & qsintveg, qsintmax, assim_param, &
     & resist_equiv_evap, resist_equiv_evap_veg, resist_equiv_tot, frac_sub_snow, frac_evap_intercept, resist_intercept, &
     & frac_evap_transp, resist_transp, vbeta3pot, frac_evap_bare, frac_evap_flood, raero, &
     & gsmean, rveget, rstruct, cimean, gpp, cresist, &
     & lalo, neighbours, resolution, ptnlev1, precip_rain, frac_age, tot_bare_soil, frac_snow_veg, frac_snow_nobio, &
     & hist_id, hist2_id, vbeta2sum, vbeta3sum, Light_Abs_Tot, &
     & Light_Tran_Tot, lai_per_level, vbeta23, leaf_ci, &
     & gs_distribution, gs_diffuco_output, gstot_component, gstot_frac,&
     & VPD, air_relhum, f_gco2_out, warnings, u_speed, profile_vbeta3, profile_rveget, &
     & delta_c13_assim, leaf_ci_out, info_limitphoto, JJ_out, assimi_lev, &
     & psi_leaf_next, assimtot, Rdtot_pft)

  !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                         :: kjit               !! Time step number (-) 
    INTEGER(i_std), INTENT(in)                         :: kjpindex           !! Domain size (-)
    INTEGER(i_std),INTENT (in)                         :: hist_id            !! _History_ file identifier (-)
    INTEGER(i_std),INTENT (in)                         :: hist2_id           !! _History_ file 2 identifier (-)
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)     :: index            !! Indeces of the points on the map (-)
    INTEGER(i_std),DIMENSION (kjpindex*(nlevels_tot+1)), INTENT (in) :: indexlai  !! Indeces of the points on the 3D map
    INTEGER(i_std),DIMENSION (kjpindex*nvm), INTENT (in) :: indexveg         !! Indeces of the points on the 3D map (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: u                  !! Eastward Lowest level wind speed (m s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: v                  !! Northward Lowest level wind speed (m s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: zlev               !! Height of first layer (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: z0m                !! Surface roughness Length for momentum (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: z0h                !! Surface roughness Length for heat (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: roughheight        !! Effective height for roughness (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: temp_sol           !! Skin temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: temp_air           !! Lowest level temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: temp_growth        !! Growth temperature (°C) - Is equal to t2m_month
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: rau                !! Air Density (kg m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: qsurf              !! Near surface air specific humidity (kg kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: qair               !! Lowest level air specific humidity (kg kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: snow               !! Snow mass (kg)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: flood_frac         !! Fraction of floodplains
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: flood_res          !! Reservoir in floodplains (estimation to avoid over-evaporation)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: pb                 !! Surface level pressure (hPa)
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)      :: evap_bare_lim   !! Limit to the bare soil evaporation when the 
                                                                             !! 11-layer hydrology is used (-)
    REAL(r_std),DIMENSION (kjpindex,nstm), INTENT (inout) :: evap_bare_lim_ns!! Limit to the bare soil evaporation when the 
                                                                             !! 11-layer hydrology is used (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: evapot             !! Soil Potential Evaporation (mm day^{-1}) 
                                                                             !! NdN This variable does not seem to be used at 
                                                                             !! all in diffuco
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: evapot_corr        !! Soil Potential Evaporation
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (in) :: frac_nobio       !! Fraction of ice,lakes,cities,... (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: totfrac_nobio      !! Total fraction of ice+lakes+cities+... (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: swnet              !! Net surface short-wave flux (W m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: swdown             !! Down-welling surface short-wave flux (W m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: coszang            !! Cosine of the solar zenith angle (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: ccanopy            !! CO2 concentration inside the canopy (ppm)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: veget              !! Fraction of vegetation type (-)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: veget_max          !! Max. fraction of vegetation type (LAI->infty)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(in)   :: lai                !! PFT leaf area index (m^{2} m^{-2})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: qsintveg           !! Water on vegetation due to interception (kg m^{-2})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: qsintmax           !! Maximum water on vegetation for interception 
                                                                             !! (kg m^{-2})
    REAL(r_std),DIMENSION (kjpindex,nvm,npco2), INTENT (in) :: assim_param   !! referenced vcmax, nue, leaf N for photosynthesis
                                                                             !! Nitrogen use Efficiency with impact of leaf age (umol CO2 (gN)-1 s-1) 
                                                                             !! Total canopy (eg leaf) N (gN m-2[ground]) 
    REAL(r_std),DIMENSION (kjpindex,2),   INTENT (in)  :: lalo               !! Geographical coordinates
    INTEGER(i_std),DIMENSION (kjpindex,NbNeighb),INTENT (in):: neighbours    !! Vector of neighbours for each 
                                                                             !! grid point (1=N, 2=E, 3=S, 4=W)
    REAL(r_std),DIMENSION (kjpindex,2), INTENT(in)     :: resolution         !! The size in m of each grid-box in X and Y
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: ptnlev1            !! 1st level of soil temperature (Kelvin)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: precip_rain        !! Rain precipitation expressed in mm/tstep
    REAL(r_std),DIMENSION (kjpindex,nvm,nleafages), INTENT (in)  :: frac_age !! Age efficiency for isoprene emissions (from STOMATE)
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)       :: tot_bare_soil      !! Total evaporating bare soil fraction 
    REAL(r_std),DIMENSION (:,:,:), INTENT (in)         :: Light_Abs_Tot      !! Absorbed radiation per level for photosynthesis
    REAL(r_std),DIMENSION (:,:,:), INTENT (in)         :: Light_Tran_Tot     !! Transmitted radiation per level for photosynthesis
    REAL(r_std), DIMENSION(:,:,:), INTENT(IN)          :: lai_per_level      !! This is the LAI per vertical level
                                                                             !! @tex $(m^{2} m^{-2})$
    REAL(r_std), DIMENSION (:,:), INTENT (in)            :: u_speed          !! wind speed at specified canopy level (m/s)
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)       :: frac_snow_veg     !! Snow cover fraction on vegeted area
    REAL(r_std), DIMENSION (kjpindex,nnobio), INTENT(in):: frac_snow_nobio   !! Snow cover fraction on non-vegeted area
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(in)   :: psi_leaf_next     !! Approximated Leaf water potential at time step n+1 (MPa) (=psi_leaf when no stress)


    !! 0.2 Output variables

    REAL(r_std),DIMENSION (kjpindex), INTENT (out)     :: resist_equiv_evap     !! Equivalent resistance to the fluxes involved in the evaporation flux (interception loss + transpiration + bare soil evaporation + floodplains evaporation) (s m^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)     :: resist_equiv_evap_veg !! Equivalent resistance to the fluxes involved in the evaporation flux from the vegetation (interception loss + transpiration + bare soil evaporation) (s m^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)     :: resist_equiv_tot      !! Equivalent resistance to the fluxes involved in the total evaporative flux (all evaporations + snow sublimation) (s m^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)     :: frac_sub_snow         !! Fraction of gridcell where snow sublimation occurs (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)     :: frac_evap_bare        !! Fraction of gridcell where bare soil evaporation occurs (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)     :: frac_evap_flood       !! Fraction of gridcell where floodplains evaporation occurs
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: gsmean                !! Mean stomatal conductance to CO2 (mol m-2 s-1) 
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: frac_evap_intercept   !! Fraction of gridcell where intercepted water evaporation occurs (-)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: resist_intercept      !! Sum of the resistances applied to intercepted water evaporation (s m^{-1})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: frac_evap_transp      !! Fraction of gridcell where transpiration occurs (-)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: resist_transp         !! Sum of the resistances applied to transpiration (s m^{-1})
    REAL(r_std),DIMENSION (kjpindex),     INTENT (out) :: raero                 !! Aerodynamic resistances (s m^{-1})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: vbeta3pot        !! Beta for potential transpiration
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: rveget           !! Stomatal resistance for the whole canopy (s m^{-1})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout) :: rstruct        !! Structural resistance for the vegetation
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: cimean           !! Mean leaf Ci (ppm)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: gpp              !! Assimilation ((gC m^{-2} dt_sechiba^{-1}), total area)        
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: cresist          !! coefficient for resistances (??)

    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: vbeta23          !! Beta for fraction of wetted foliage that will
                                                                           !! transpire once intercepted water has evaporated (-)
    REAL(r_std),DIMENSION (kjpindex,nvm,nlevels_tot), INTENT (out)       :: leaf_ci    !! intercellular CO2 concentration (ppm)
    REAL(r_Std),DIMENSION (kjpindex,nvm,nlevels_tot), INTENT (out)       :: gs_distribution
    REAL(r_Std),DIMENSION (kjpindex,nvm,nlevels_tot), INTENT (out)       :: gs_diffuco_output
    REAL(r_Std),DIMENSION (kjpindex,nvm), INTENT (out)                   :: f_gco2_out

    REAL(r_Std),DIMENSION (kjpindex,nvm,nlevels_tot), INTENT (out)       :: gstot_component
    REAL(r_Std),DIMENSION (kjpindex,nvm,nlevels_tot), INTENT (out)       :: gstot_frac
    REAL(r_std),DIMENSION (kjpindex,nvm,nlevels_tot), INTENT (out)       :: profile_vbeta3
    REAL(r_std),DIMENSION (kjpindex,nvm,nlevels_tot), INTENT (out)       :: profile_rveget

    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)                   :: delta_c13_assim !! C13 concentration in delta notation
                                                                                            !! @tex $ permille $ @endtex (per thousand)   
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)                   :: leaf_ci_out     !! Ci    
    REAL(r_std),DIMENSION (kjpindex,nvm,nlevels_tot+1),INTENT(out)       :: info_limitphoto
    REAL(r_std),DIMENSION (kjpindex,nvm,nlevels_tot),INTENT(out)         :: JJ_out 
    REAL(r_std),DIMENSION (kjpindex,nvm,nlevels_tot),INTENT(out)         :: assimi_lev  !! assimilation (at all LAI levels)
                                                                                        !! ($\mumolm^{-2}s^{-1}$)
    REAL(r_std), DIMENSION(kjpindex,nvm),INTENT(out)                    :: assimtot          !! total assimilation 
                                                                                 !! @tex ($\mu mol m^{-2} s^{-1}$) @endtex
    REAL(r_std), DIMENSION(kjpindex,nvm),INTENT(out)                    :: Rdtot_pft          !! total respiration 
                                                                                 !! @tex ($\mu mol m^{-2} s^{-1}$) @endtex

    !! 0.3 Modified variables

    REAL(r_std),DIMENSION (kjpindex, nvm), INTENT (inout)                :: humrel      !! Soil moisture stress (within range 0 to 1)
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)                     :: air_relhum  !! Air relative humidity at 2m
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)                     :: VPD         !! Vapor Pressure Deficit (kPa)
 
    REAL(r_std),DIMENSION(kjpindex,nvm,nwarns), INTENT(inout) &
                                                          :: warnings      !! A warning counter
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)      :: q_cdrag       !! Surface drag coefficient  (-)


    !! 0.4 Local variables
    REAL(r_std),DIMENSION(kjpindex)          :: vbeta2sum, vbeta3sum
    INTEGER(i_std)                           :: ilai
    CHARACTER(LEN=4)                         :: laistring
    CHARACTER(LEN=80)                        :: var_name           !! To store variables names for I/O
    REAL(r_std),DIMENSION(kjpindex)          :: qsatt              !! Surface saturated humidity (kg kg^{-1})
    REAL(r_std),DIMENSION (kjpindex,nvm)     :: cim                !! Intercellular CO2 over nlevels_tot

    ! Isotope_c13
    INTEGER(i_std), DIMENSION(kjpindex)                     :: index_assi        !! indices (unitless)
    INTEGER(i_std)                                          :: nia,inia,nina,iainia   !! counter/indices (unitless)
    INTEGER(i_std)                                          :: ji,jl,jv          !! Indices (-)    
    REAL(r_std), DIMENSION(nlevels_tot+1)                   :: laitab            !! LAI per layer (m^2.m^{-2})

!_ ================================================================================================================================
    

    gs_distribution = 0.0d0
    wind(:) = SQRT (u(:)*u(:) + v(:)*v(:))
  
  !! 1. Calculate the different coefficients

    IF (.NOT.ldq_cdrag_from_gcm) THEN
        ! Case 1a)
       CALL diffuco_aero (kjpindex, kjit, u, v, zlev, z0h, z0m, roughheight, temp_sol, temp_air, &
                          qsurf, qair, snow, q_cdrag)
    ENDIF

    ! Case 1b)
    CALL diffuco_raerod (kjpindex, u, v, q_cdrag, raero)

  !! 2. Make an estimation of the saturated humidity at the surface

    CALL qsatcalc (kjpindex, temp_sol, pb, qsatt)

  !! 3. Calculate the beta coefficient for sublimation
  
    CALL diffuco_snow (kjpindex, qair, qsatt, rau, u, v, q_cdrag, &
         snow, frac_nobio, totfrac_nobio, frac_snow_veg, frac_snow_nobio, &
         raero, frac_sub_snow)

    CALL diffuco_flood (kjpindex, qair, qsatt, rau, u, v, q_cdrag, evapot, evapot_corr, &
         & flood_frac, flood_res, raero, frac_evap_flood)

  !! 4. Calculate the beta coefficient for interception

    CALL diffuco_inter (kjpindex, qair, qsatt, rau, u, v, q_cdrag, humrel, veget, &
       & qsintveg, qsintmax, rstruct, raero, frac_evap_intercept, resist_intercept, vbeta23) 

  !! 5. Calculate the beta coefficient for transpiration
    CALL diffuco_trans_co2 (kjpindex, swdown, pb, qair, temp_air, &
         temp_growth, rau, q_cdrag, humrel, &
         assim_param, ccanopy, lai, &
         veget_max, qsintveg, qsintmax, frac_evap_transp, resist_transp, vbeta3pot, rveget, &
         rstruct, cimean, gsmean, gpp, cresist, vbeta23, Light_Abs_Tot,&
         Light_Tran_Tot, lai_per_level, gs_distribution, gs_diffuco_output, &
         VPD, air_relhum, f_gco2_out, gstot_component, gstot_frac, veget, warnings, u_speed, profile_vbeta3, &
         profile_rveget, index_assi, leaf_ci, assimtot, assimi_lev, hist_id, &
         indexveg, indexlai, index, kjit, cim, info_limitphoto, JJ_out,&
         psi_leaf_next, raero, Rdtot_pft)

    ! calculate C13 isotopic signature
    IF (ok_c13) THEN
       CALL isotope_c13(kjpindex, index_assi, leaf_ci, ccanopy, lai_per_level, &
            assimtot, delta_c13_assim, leaf_ci_out, nvm, &
            nlevels_tot, assimi_lev)
    ENDIF

    !
    !biogenic emissions
    !

    IF ( ok_bvoc ) THEN
       CALL chemistry_bvoc (kjpindex, swdown, coszang, temp_air, &
            temp_sol, ptnlev1, precip_rain, humrel, veget_max, &
            lai, frac_age, lalo, ccanopy, cim, wind, snow, &
            veget, hist_id, hist2_id, kjit, index, &
            indexlai, indexveg)
    ENDIF
    !
    ! combination of coefficient : alpha and beta coefficient
    ! beta coefficient for bare soil
    !

    CALL diffuco_bare (kjpindex, tot_bare_soil, veget_max, frac_evap_intercept, resist_intercept, frac_evap_transp, resist_transp, raero, &
                          evap_bare_lim, evap_bare_lim_ns, frac_evap_bare)

  !! 6. Combine the alpha and beta coefficients

    CALL diffuco_comb (kjpindex, humrel, rau, u, v, q_cdrag, pb, qair, temp_sol, temp_air, &
        snow, veget, lai,  &    
        tot_bare_soil, frac_evap_flood, raero, frac_sub_snow, frac_evap_intercept, resist_intercept, &
        frac_evap_transp, resist_transp, frac_evap_bare, &
        evap_bare_lim, evap_bare_lim_ns, veget_max, resist_equiv_evap_veg, resist_equiv_evap, resist_equiv_tot, &
        qsintmax, vbeta2sum, vbeta3sum)


    CALL xios_orchidee_send_field("q_cdrag",q_cdrag)
    CALL xios_orchidee_send_field("raero",raero)
    CALL xios_orchidee_send_field("wind",wind)
    CALL xios_orchidee_send_field("qsatt",qsatt)
    CALL xios_orchidee_send_field("coszang",coszang)
    CALL xios_orchidee_send_field('cim', cim)
    CALL xios_orchidee_send_field("test_overPET_intensity",test_overPET_intensity)
    CALL xios_orchidee_send_field("test_overPET_freq",test_overPET_freq)

    IF ( .NOT. almaoutput ) THEN
       CALL histwrite_p(hist_id, 'raero', kjit, raero, kjpindex, index)
       CALL histwrite_p(hist_id, 'cdrag', kjit, q_cdrag, kjpindex, index)
       CALL histwrite_p(hist_id, 'Wind', kjit, wind, kjpindex, index)
       CALL histwrite_p(hist_id, 'qsatt', kjit, qsatt, kjpindex, index)
       CALL histwrite_p(hist_id, 'cim', kjit, cim, kjpindex*nvm, indexveg)

       IF ( hist2_id > 0 ) THEN
          CALL histwrite_p(hist2_id, 'raero', kjit, raero, kjpindex, index)
          CALL histwrite_p(hist2_id, 'cdrag', kjit, q_cdrag, kjpindex, index)
          CALL histwrite_p(hist2_id, 'Wind', kjit, wind, kjpindex, index)
          CALL histwrite_p(hist2_id, 'qsatt', kjit, qsatt, kjpindex, index)
       ENDIF
    ELSE
       CALL histwrite_p(hist_id, 'cim', kjit, cim, kjpindex*nvm, indexveg)
    ENDIF

    IF (printlev>=3) WRITE (numout,*) ' diffuco_main done '

  END SUBROUTINE diffuco_main

!!  =============================================================================================================================
!! SUBROUTINE: diffuco_finalize
!!
!>\BRIEF          Write to restart file
!!
!! DESCRIPTION:   This subroutine writes the module variables and variables calculated in diffuco
!!                to restart file
!!
!! RECENT CHANGE(S): None
!! REFERENCE(S): None
!! FLOWCHART: None
!! \n
!_ ==============================================================================================================================
  SUBROUTINE diffuco_finalize (kjit, kjpindex, rest_id, rstruct )

    !! 0. Variable and parameter declaration
    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                         :: kjit             !! Time step number (-) 
    INTEGER(i_std), INTENT(in)                         :: kjpindex         !! Domain size (-)
    INTEGER(i_std),INTENT (in)                         :: rest_id          !! _Restart_ file identifier (-)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: rstruct          !! Structural resistance for the vegetation

    !! 0.4 Local variables
    INTEGER                                            :: ilai
    CHARACTER(LEN=4)                                   :: laistring
    CHARACTER(LEN=80)                                  :: var_name        

!_ ================================================================================================================================
    
  !! 1. Prepare the restart file for the next simulation
    IF (printlev>=3) WRITE (numout,*) 'Complete restart file with DIFFUCO variables '
    
    CALL restput_p (rest_id, 'rstruct', nbp_glo, nvm, 1, kjit, rstruct, 'scatter',  nbp_glo, index_g)

  END SUBROUTINE diffuco_finalize


!! ================================================================================================================================
!! SUBROUTINE		 			: diffuco_clear
!!
!>\BRIEF					Housekeeping module to deallocate the variables
!! rstruct and raero
!!
!! DESCRIPTION				        : Housekeeping module to deallocate the variables
!! rstruct and raero
!!
!! RECENT CHANGE(S)                             : None
!!
!! MAIN OUTPUT VARIABLE(S)	                : None
!!
!! REFERENCE(S)				        : None
!!
!! FLOWCHART                                    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE diffuco_clear()

    ! Deallocate and reset variables in chemistry module
    CALL chemistry_clear
    IF (ALLOCATED  (test_overPET_intensity)) DEALLOCATE (test_overPET_intensity)
    IF (ALLOCATED  (test_overPET_freq)) DEALLOCATE (test_overPET_freq)

  END SUBROUTINE diffuco_clear


!! ================================================================================================================================
!! SUBROUTINE	: diffuco_aero
!!
!>\BRIEF	This module first calculates the surface drag 
!! coefficient, for cases in which the surface drag coefficient is NOT provided by the coupled 
!! atmospheric model LMDZ or when the flag ldq_cdrag_from_gcm is set to FALSE 
!!
!! DESCRIPTION	: Computes the surface drag coefficient, for cases 
!! in which it is NOT provided by the coupled atmospheric model LMDZ. The module first uses the 
!! meteorolgical input to calculate the Richardson Number, which is an indicator of atmospheric 
!! stability in the surface layer. The formulation used to find this surface drag coefficient is 
!! dependent on the stability determined. \n
!!
!! Designation of wind speed
!! \latexonly 
!!     \input{diffucoaero1.tex}
!! \endlatexonly
!!
!! Calculation of geopotential. This is the definition of Geopotential height (e.g. Jacobson 
!! eqn.4.47, 2005). (required for calculation of the Richardson Number)
!! \latexonly 
!!     \input{diffucoaero2.tex}
!! \endlatexonly
!! 
!! \latexonly 
!!     \input{diffucoaero3.tex}
!! \endlatexonly
!!
!! Calculation of the virtual air temperature at the surface (required for calculation
!! of the Richardson Number)
!! \latexonly 
!!     \input{diffucoaero4.tex}
!! \endlatexonly
!!
!! Calculation of the virtual surface temperature (required for calculation of th
!! Richardson Number)
!! \latexonly 
!!     \input{diffucoaero5.tex}
!! \endlatexonly
!!
!! Calculation of the squared wind shear (required for calculation of the Richardson
!! Number)
!! \latexonly 
!!     \input{diffucoaero6.tex}
!! \endlatexonly
!! 
!! Calculation of the Richardson Number. The Richardson Number is defined as the ratio 
!! of potential to kinetic energy, or, in the context of atmospheric science, of the
!! generation of energy by wind shear against consumption
!! by static stability and is an indicator of flow stability (i.e. for when laminar flow 
!! becomes turbulent and vise versa). It is approximated using the expression below:
!! \latexonly 
!!     \input{diffucoaero7.tex}
!! \endlatexonly
!!
!! The Richardson Number hence calculated is subject to a minimum value:
!! \latexonly 
!!     \input{diffucoaero8.tex}
!! \endlatexonly
!! 
!! Computing the drag coefficient. We add the add the height of the vegetation to the 
!! level height to take into account that the level 'seen' by the vegetation is actually 
!! the top of the vegetation. Then we we can subtract the displacement height.
!! \latexonly 
!!     \input{diffucoaero9.tex}
!! \endlatexonly
!! 
!! For the stable case (i.e $R_i$ $\geq$ 0)
!! \latexonly 
!!     \input{diffucoaero10.tex}
!! \endlatexonly
!!
!! \latexonly 
!!     \input{diffucoaero11.tex}
!! \endlatexonly
!!          
!! For the unstable case (i.e. $R_i$ < 0)
!! \latexonly 
!!     \input{diffucoaero12.tex}
!! \endlatexonly
!!
!! \latexonly 
!!     \input{diffucoaero13.tex}
!! \endlatexonly
!!               
!! If the Drag Coefficient becomes too small than the surface may uncouple from the atmosphere.
!! To prevent this, a minimum limit to the drag coefficient is defined as:
!!
!! \latexonly 
!!     \input{diffucoaero14.tex}
!! \endlatexonly
!! 
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): q_cdrag
!!
!! REFERENCE(S)	: 
!! - de Noblet-Ducoudré, N, Laval, K & Perrier, A, 1993. SECHIBA, a new set of parameterisations
!! of the hydrologic exchanges at the land-atmosphere interface within the LMD Atmospheric General
!! Circulation Model. Journal of Climate, 6, pp.248-273
!! - Guimberteau, M, 2010. Modélisation de l'hydrologie continentale et influences de l'irrigation
!! sur le cycle de l'eau, PhD Thesis, available from:
!! http://www.sisyphe.upmc.fr/~guimberteau/docs/manuscrit_these.pdf
!! - Jacobson M.Z., Fundamentals of Atmospheric Modeling (2nd Edition), published Cambridge 
!! University Press, ISBN 0-521-54865-9
!!
!! FLOWCHART	:
!! \latexonly 
!!     \includegraphics[scale=0.5]{diffuco_aero_flowchart.png}
!! \endlatexonly
!! \n
!_ ================================================================================================================================

  SUBROUTINE diffuco_aero (kjpindex, kjit, u, v, zlev, z0h, z0m, roughheight, temp_sol, temp_air, &
                           qsurf, qair, snow, q_cdrag)

  !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                          :: kjpindex, kjit   !! Domain size
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: u                !! Eastward Lowest level wind speed (m s^{-1}) 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: v                !! Northward Lowest level wind speed (m s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: zlev             !! Height of first atmospheric layer (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: z0h               !! Surface roughness Length for heat (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: z0m               !! Surface roughness Length for momentum (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: roughheight      !! Effective roughness height (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: temp_sol         !! Ground temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: temp_air         !! Lowest level temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: qsurf            !! near surface specific air humidity (kg kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: qair             !! Lowest level specific air humidity (kg kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: snow             !! Snow mass (kg)

    !! 0.2 Output variables
    
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)      :: q_cdrag          !! Surface drag coefficient  (-)

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                                      :: ji, jv
    REAL(r_std)                                         :: speed, zg, zdphi, ztvd, ztvs, zdu2
    REAL(r_std)                                         :: zri, cd_neut, zscf, cd_tmp
!_ ================================================================================================================================

  !! 1. Initialisation

    ! test if we have to work with q_cdrag or to calcul it
    DO ji=1,kjpindex
       
       !! 1a).1 Designation of wind speed
       !! \latexonly 
       !!     \input{diffucoaero1.tex}
       !! \endlatexonly
       speed = wind(ji)
    
       !! 1a).2 Calculation of geopotentiel
       !! This is the definition of Geopotential height (e.g. Jacobson eqn.4.47, 2005). (required
       !! for calculation of the Richardson Number)
       !! \latexonly 
       !!     \input{diffucoaero2.tex}
       !! \endlatexonly
       zg = zlev(ji) * cte_grav
      
       !! \latexonly 
       !!     \input{diffucoaero3.tex}
       !! \endlatexonly
       zdphi = zg/cp_air
       
       !! 1a).3 Calculation of the virtual air temperature at the surface 
       !! required for calculation of the Richardson Number
       !! \latexonly 
       !!     \input{diffucoaero4.tex}
       !! \endlatexonly
       ztvd = (temp_air(ji) + zdphi / (un + rvtmp2 * qair(ji))) * (un + retv * qair(ji)) 
       
       !! 1a).4 Calculation of the virtual surface temperature 
       !! required for calculation of the Richardson Number
       !! \latexonly 
       !!     \input{diffucoaero5.tex}
       !! \endlatexonly
       ztvs = temp_sol(ji) * (un + retv * qsurf(ji))
     
       !! 1a).5 Calculation of the squared wind shear 
       !! required for calculation of the Richardson Number
       !! \latexonly 
       !!     \input{diffucoaero6.tex}
       !! \endlatexonly
       zdu2 = MAX(cepdu2,speed**2)
       
       !! 1a).6 Calculation of the Richardson Number
       !!  The Richardson Number is defined as the ratio of potential to kinetic energy, or, in the 
       !!  context of atmospheric science, of the generation of energy by wind shear against consumption
       !!  by static stability and is an indicator of flow stability (i.e. for when laminar flow 
       !!  becomes turbulent and vise versa).\n
       !!  It is approximated using the expression below:
       !!  \latexonly 
       !!     \input{diffucoaero7.tex}
       !! \endlatexonly
       zri = zg * (ztvd - ztvs) / (zdu2 * ztvd)
      
       !! The Richardson Number hence calculated is subject to a minimum value:
       !! \latexonly 
       !!     \input{diffucoaero8.tex}
       !! \endlatexonly       
       zri = MAX(MIN(zri,5.),-5.)
       
       !! 1a).7 Computing the drag coefficient
       !!  We add the add the height of the vegetation to the level height to take into account
       !!  that the level 'seen' by the vegetation is actually the top of the vegetation. Then we 
       !!  we can subtract the displacement height.
       !! \latexonly 
       !!     \input{diffucoaero9.tex}
       !! \endlatexonly

       !! 7.0 Snow smoothering
       !! Snow induces low levels of turbulence.
       !! Sensible heat fluxes can therefore be reduced of ~1/3. Pomeroy et al., 1998
       cd_neut = ct_karman ** 2. / ( LOG( (zlev(ji) + roughheight(ji)) / z0m(ji) ) * LOG( (zlev(ji) + roughheight(ji)) / z0h(ji) ) )
       
       !! 1a).7.1 - for the stable case (i.e $R_i$ $\geq$ 0)
       IF (zri .GE. zero) THEN
          
          !! \latexonly 
          !!     \input{diffucoaero10.tex}
          !! \endlatexonly
          zscf = SQRT(un + cd * ABS(zri))
         
          !! \latexonly 
          !!     \input{diffucoaero11.tex}
          !! \endlatexonly          
          cd_tmp=cd_neut/(un + trois * cb * zri * zscf)
       ELSE
          
          !! 1a).7.2 - for the unstable case (i.e. $R_i$ < 0)
          !! \latexonly 
          !!     \input{diffucoaero12.tex}
          !! \endlatexonly
          zscf = un / (un + trois * cb * cc * cd_neut * SQRT(ABS(zri) * &
               & ((zlev(ji) + roughheight(ji)) / z0m(ji))))

          !! \latexonly 
          !!     \input{diffucoaero13.tex}
          !! \endlatexonly               
          cd_tmp=cd_neut * (un - trois * cb * zri * zscf)
       ENDIF
       
       !! If the Drag Coefficient becomes too small than the surface may uncouple from the atmosphere.
       !! To prevent this, a minimum limit to the drag coefficient is defined as:
       
       !! \latexonly 
       !!     \input{diffucoaero14.tex}
       !! \endlatexonly
       !!
       q_cdrag(ji) = MAX(cd_tmp, min_qc/MAX(speed,min_wind))

       ! In some situations it might be useful to give an upper limit on the cdrag as well. 
       ! The line here should then be uncommented.
      !q_cdrag(ji) = MIN(q_cdrag(ji), 0.5/MAX(speed,min_wind))

    END DO

    IF (printlev_loc>=3) WRITE (numout,*) ' not ldqcdrag_from_gcm : diffuco_aero done '

  END SUBROUTINE diffuco_aero


!! ================================================================================================================================
!! SUBROUTINE    : diffuco_snow
!!
!>\BRIEF         This subroutine computes the beta coefficient for snow sublimation.
!!
!! DESCRIPTION   : This routine computes beta coefficient for snow sublimation, which
!! integrates the snow on both vegetation and other surface types (e.g. ice, lakes,
!! cities etc.) \n
!!
!! A critical depth of snow (snowcri) is defined to calculate the fraction of each grid-cell
!! that is covered with snow (snow/snowcri) while the remaining part is snow-free.
!! We also carry out a first calculation of sublimation (subtest) to lower down the beta
!! coefficient if necessary (if subtest > snow). This is a predictor-corrector test. 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): :: frac_sub_snow
!!
!! REFERENCE(S) :
!! - de Noblet-Ducoudré, N, Laval, K & Perrier, A, 1993. SECHIBA, a new set of parameterisations
!! of the hydrologic exchanges at the land-atmosphere interface within the LMD Atmospheric General
!! Circulation Model. Journal of Climate, 6, pp. 248-273
!! - Guimberteau, M, 2010. Modélisation de l'hydrologie continentale et influences de l'irrigation
!! sur le cycle de l'eau, PhD Thesis, available from:
!! http://www.sisyphe.upmc.fr/~guimberteau/docs/manuscrit_these.pdf
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  
SUBROUTINE diffuco_snow (kjpindex, qair, qsatt, rau, u, v,q_cdrag, &
       & snow, frac_nobio, totfrac_nobio, frac_snow_veg, frac_snow_nobio, &
       & raero, frac_sub_snow)

  !! 0. Variable and parameter declaration
    
    !! 0.1 Input variables
 
    INTEGER(i_std), INTENT(in)                           :: kjpindex       !! Domain size (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: qair           !! Lowest level specific air humidity (kg kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: qsatt          !! Surface saturated humidity (kg kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: rau            !! Air density (kg m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: u              !! Eastward Lowest level wind speed (m s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: v              !! Northward Lowest level wind speed (m s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: q_cdrag        !! Surface drag coefficient  (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: snow           !! Snow mass (kg m^{-2})
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (in) :: frac_nobio     !! Fraction of ice, lakes, cities etc. (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: totfrac_nobio  !! Total fraction of ice, lakes, cities etc. (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)         :: frac_snow_veg  !! Snow cover fraction on vegeted area
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT(in)  :: frac_snow_nobio!! Snow cover fraction on non-vegeted area
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)         :: raero          !! Snow cover fraction on vegeted area
    
    !! 0.2 Output variables

    REAL(r_std),DIMENSION (kjpindex), INTENT (out)       :: frac_sub_snow  !! Fraction of the gridcell where snow sublimation occurs (-)
    
    !! 0.3 Modified variables

    !! 0.4 Local variables

    REAL(r_std),DIMENSION (kjpindex)                     :: frac_sub_snow_veg_predict    !! Prediction of the fraction of the vegetated area where snow sublimation occurs (-)
    REAL(r_std),DIMENSION (kjpindex)                     :: frac_sub_snow_veg_corr       !! Correction of the fraction of the vegetated area where snow sublimation occurs (-)
    REAL(r_std)                                          :: frac_sub_snow_nobio_predict  !! Prediction of the fraction of the non-vegetated area where snow sublimation occurs (-) 
    REAL(r_std)                                          :: frac_sub_snow_nobio_corr     !! Correction of the fraction of the non-vegetated area where snow sublimation occurs (-)
    REAL(r_std),DIMENSION (kjpindex)                     :: frac_sub_snow_veg            !! Fraction of the vegetated area where snow sublimation occurs (-)
    REAL(r_std)                                          :: frac_sub_snow_nobio          !! Fraction of the non-vegetated area where snow sublimation occurs (-)
    REAL(r_std),DIMENSION (kjpindex)                     :: frac_sub_snow_nobio_tot      !! Total fraction of the non-vegetated area where snow sublimation occurs (-) (useless?)
    REAL(r_std)                                          :: subtest        !! Sublimation for test (kg m^{-2})
    REAL(r_std)                                          :: zrapp          !! Modified factor (ratio)
    REAL(r_std)                                          :: speed          !! Wind speed (m s^{-1})
    REAL(r_std)                                          :: vbeta1_add     !! Beta for sublimation (ratio)
    INTEGER(i_std)                                       :: ji, jv         !! Indices (-)
!_ ================================================================================================================================

  !! 1. Calculate beta coefficient for snow sublimation on the vegetation\n

    DO ji=1,kjpindex  ! Loop over # pixels - domain size

       ! Fraction of mesh that can sublimate snow
       frac_sub_snow_veg_predict(ji) = (un - totfrac_nobio(ji)) * frac_snow_veg(ji)

       ! Limitation of sublimation in case of snow amounts smaller 
       ! than the atmospheric demand. 
       speed = MAX(min_wind, wind(ji))

       subtest = frac_sub_snow_veg_predict(ji) * dt_sechiba * rau(ji) * &
               & ( qsatt(ji) - qair(ji) ) / raero(ji)

       frac_sub_snow_veg(ji) = frac_sub_snow_veg_predict(ji)

       IF ( subtest .GT. min_sechiba ) THEN
          zrapp = snow(ji) / subtest
          IF ( zrapp .LT. un ) THEN
             frac_sub_snow_veg_corr(ji) = frac_sub_snow_veg_predict(ji) * zrapp
             frac_sub_snow_veg(ji) = frac_sub_snow_veg_corr(ji) 
          ENDIF
       ENDIF

    END DO ! Loop over # pixels - domain size

  !! 2. Add the beta coefficients calculated from other surfaces types (snow on ice,lakes, cities...)

    frac_sub_snow_nobio_tot(:) = zero

    DO jv = iice, iice ! Loop over # other surface types
       DO ji=1,kjpindex ! Loop over # pixels - domain size

           frac_sub_snow_nobio_predict = frac_nobio(ji,jv) * frac_snow_nobio(ji, jv)

          ! Limitation of sublimation in case of snow amounts smaller than
          ! the atmospheric demand. 
          speed = MAX(min_wind, wind(ji))

          !!     Limitation of sublimation by the snow accumulated on the ground 
          !!     A first approximation is obtained with the old values of
          !!     qair and qsol_sat: function of temp-sol and pb. (see call of qsatcalc)
          subtest = frac_sub_snow_nobio_predict * dt_sechiba * rau(ji) * &
                & ( qsatt(ji) - qair(ji) ) / raero(ji)

           frac_sub_snow_nobio = frac_sub_snow_nobio_predict

           IF ( subtest .GT. min_sechiba ) THEN
              zrapp = snow(ji) / subtest
              IF ( zrapp .LT. un ) THEN
                 frac_sub_snow_nobio_corr = frac_sub_snow_nobio_predict*zrapp
                 frac_sub_snow_nobio = frac_sub_snow_nobio_corr
              ENDIF
           ENDIF

           frac_sub_snow_nobio_tot(ji) = frac_sub_snow_nobio_tot(ji) + frac_sub_snow_nobio

       ENDDO ! Loop over # pixels - domain size
    ENDDO ! Loop over # other surface types

    frac_sub_snow(:) = frac_sub_snow_veg(:) + frac_sub_snow_nobio_tot(:)

    IF (printlev>=3) WRITE (numout,*) ' diffuco_snow done '

  END SUBROUTINE diffuco_snow


!! ================================================================================================================================
!! SUBROUTINE		 			: diffuco_flood 
!!
!>\BRIEF				       	This routine computes partial beta coefficient : floodplains
!!
!! DESCRIPTION				        : 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S)	                : frac_evap_flood
!!
!! REFERENCE(S)				        : None
!!
!! FLOWCHART                                    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE diffuco_flood (kjpindex, qair, qsatt, rau, u, v, q_cdrag, evapot, evapot_corr, &
       & flood_frac, flood_res, raero, frac_evap_flood)

    ! interface description
    ! input scalar 
    INTEGER(i_std), INTENT(in)                               :: kjpindex   !! Domain size
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: qair       !! Lowest level specific humidity
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: qsatt      !! Surface saturated humidity
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: rau        !! Density
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: u          !! Lowest level wind speed 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: v          !! Lowest level wind speed
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: q_cdrag    !! Surface drag coefficient  (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: flood_res  !! water mass in flood reservoir
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: flood_frac !! fraction of floodplains
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: evapot     !! Potential evaporation
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: evapot_corr!! Potential evaporation2
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: raero      !! water mass in flood reservoir
    ! output fields
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: frac_evap_flood          !! Fraction of the gridcell where floodplains evaporatio occurs (-)

    ! local declaration
    REAL(r_std),DIMENSION (kjpindex)                         :: frac_evap_flood_predict  !! Prediction of the fraction of the gridcell where floodplains evaporatio occurs (-)
    REAL(r_std),DIMENSION (kjpindex)                         :: frac_evap_flood_corr     !! Correction of the fraction of the gridcell where floodplains evaporatio occurs (-)
    REAL(r_std)                                              :: subtest, zrapp, speed
    INTEGER(i_std)                                           :: ji, jv

!_ ================================================================================================================================
    !
    ! beta coefficient for sublimation for floodplains
    !
    DO ji=1,kjpindex
       !

       IF (evapot(ji) .GT. min_sechiba) THEN
          frac_evap_flood_predict(ji) = flood_frac(ji) *evapot_corr(ji)/evapot(ji)
       ELSE
          frac_evap_flood_predict(ji) = flood_frac(ji)
       ENDIF
       !
       ! -- Limitation of evaporation in case of water amounts smaller than
       !    the atmospheric demand. 
       
       !
       speed = MAX(min_wind, SQRT (u(ji)*u(ji) + v(ji)*v(ji)))
       !
       subtest = frac_evap_flood_predict(ji) * dt_sechiba * rau(ji) * &
               & ( qsatt(ji) - qair(ji) ) / raero(ji)
       !  
       frac_evap_flood(ji) = frac_evap_flood_predict(ji)
       IF ( subtest .GT. min_sechiba ) THEN
          zrapp = flood_res(ji) / subtest
          IF ( zrapp .LT. un ) THEN
             frac_evap_flood_corr(ji) = frac_evap_flood_predict(ji) * zrapp
             frac_evap_flood(ji)=frac_evap_flood_corr(ji)
          ENDIF
       ENDIF
       !
    END DO

    IF (printlev>=3) WRITE (numout,*) ' diffuco_flood done '

  END SUBROUTINE diffuco_flood


!! ================================================================================================================================
!! SUBROUTINE    : diffuco_inter
!!
!>\BRIEF	 This routine computes the partial beta coefficient
!! for the interception for each type of vegetation
!!
!! DESCRIPTION   : We first calculate the dry and wet parts of each PFT (wet part = qsintveg/qsintmax).
!! The former is submitted to transpiration only (vbeta3 coefficient, calculated in 
!! diffuco_trans_co2), while the latter is first submitted to interception loss 
!! (vbeta2 coefficient) and then to transpiration once all the intercepted water has been evaporated 
!! (vbeta23 coefficient). Interception loss is also submitted to a predictor-corrector test, 
!! as for snow sublimation. \n
!!
!! \latexonly 
!!     \input{diffucointer1.tex}
!! \endlatexonly
!! Calculate the wet fraction of vegetation as  the ration between the intercepted water and the maximum water 
!! on the vegetation. This ratio defines the wet portion of vegetation that will be submitted to interception loss.
!!
!! \latexonly 
!!     \input{diffucointer2.tex}
!! \endlatexonly
!!
!! Calculation of $\beta_3$, the canopy transpiration resistance
!! \latexonly 
!!     \input{diffucointer3.tex}
!! \endlatexonly            
!! 
!! We here determine the limitation of interception loss by the water stored on the leaf. 
!! A first approximation of interception loss is obtained using the old values of
!! qair and qsol_sat, which are functions of temp-sol and pb. (see call of 'qsatcalc')
!! \latexonly 
!!     \input{diffucointer4.tex}
!! \endlatexonly
!!
!! \latexonly
!!     \input{diffucointer5.tex}
!! \endlatexonly
!!
!! \latexonly 
!!     \input{diffucointer6.tex}
!! \endlatexonly
!!
!! Once the whole water stored on foliage has evaporated, transpiration can take place on the fraction
!! 'zqsvegrap'.
!! \latexonly 
!!     \input{diffucointer7.tex}
!! \endlatexonly
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::frac_evap_intercept ::resist_intercept ::vbeta23
!!
!! REFERENCE(S) :
!! - de Noblet-Ducoudré, N, Laval, K & Perrier, A, 1993. SECHIBA, a new set of parameterisations
!! of the hydrologic exchanges at the land-atmosphere interface within the LMD Atmospheric General
!! Circulation Model. Journal of Climate, 6, pp. 248-273
!! - Guimberteau, M, 2010. Modélisation de l'hydrologie continentale et influences de l'irrigation
!! sur le cycle de l'eau, PhD Thesis, available from:
!! http://www.sisyphe.upmc.fr/~guimberteau/docs/manuscrit_these.pdf
!! - Perrier, A, 1975. Etude physique de l'évaporation dans les conditions naturelles. Annales 
!! Agronomiques, 26(1-18): pp. 105-123, pp. 229-243
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE diffuco_inter (kjpindex, qair, qsatt, rau, u, v, q_cdrag, humrel, veget, &
     & qsintveg, qsintmax, rstruct, raero, frac_evap_intercept, resist_intercept, vbeta23)
   
  !! 0 Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                           :: kjpindex   !! Domain size (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: qair       !! Lowest level specific air humidity (kg kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: qsatt      !! Surface saturated humidity (kg kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: rau        !! Air Density (kg m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: u          !! Eastward Lowest level wind speed (m s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: v          !! Northward Lowest level wind speed (m s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: q_cdrag    !! Surface drag coefficient  (-)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)    :: humrel     !! Soil moisture stress (within range 0 to 1)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)    :: veget      !! vegetation fraction for each type (fraction)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)    :: qsintveg   !! Water on vegetation due to interception (kg m^{-2})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)    :: qsintmax   !! Maximum water on vegetation (kg m^{-2})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)    :: rstruct    !! architectural resistance (s m^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: raero      !! Aerodynamic resistance (s m^{-1})
    
    !! 0.2 Output variables
    
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)   :: frac_evap_intercept    !! Fraction of the gridcell where intercted water evaporation occurs (-)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)   :: resist_intercept       !! Sum of the resistances applied to intercepted water evaporation (s m^{-1}) 
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)   :: vbeta23                !! Beta for fraction of wetted foliage that will 
                                                                                   !! transpire (-)

    !! 0.4 Local variables

    REAL(r_std),DIMENSION (kjpindex,nvm)                 :: frac_evap_intercept_predict    !! Prediction of the fraction of the gridcell where intercted water evaporation occurs (-)
    REAL(r_std),DIMENSION (kjpindex,nvm)                 :: frac_evap_intercept_corr       !! Correction of the fraction of the gridcell where intercted water evaporation occurs (-)
    REAL(r_std),DIMENSION (kjpindex,nvm)                 :: vbeta_intercept                !! Beta for interception loss (-) (for easier calculation of beta23)
    INTEGER(i_std)                                       :: ji, jv                               !! (-), (-)
    REAL(r_std)                                          :: zqsvegrap, ziltest, zrapp, speed     !!
!_ ================================================================================================================================

  !! 1. Initialize

    vbeta_intercept(:,:) = zero
    frac_evap_intercept(:,:) = zero
    resist_intercept(:,:) = un
    vbeta23(:,:) = zero
   
  !! 2. The beta coefficient for interception by vegetation. 
    
    DO jv = 2,nvm

      DO ji=1,kjpindex

         IF (veget(ji,jv) .GT. min_sechiba .AND. qsintveg(ji,jv) .GT. zero ) THEN

            zqsvegrap = zero
            IF (qsintmax(ji,jv) .GT. min_sechiba ) THEN

            !! \latexonly 
            !!     \input{diffucointer1.tex}
            !! \endlatexonly
            !!
            !! We calculate the wet fraction of vegetation as  the ration between the intercepted water and the maximum water 
            !! on the vegetation. This ratio defines the wet portion of vegetation that will be submitted to interception loss.
            !!
                zqsvegrap = MAX(zero, qsintveg(ji,jv) / qsintmax(ji,jv))
            END IF

            !! \latexonly 
            !!     \input{diffucointer2.tex}
            !! \endlatexonly
            speed = MAX(min_wind, wind(ji))

            !! Calculation of $\beta_3$, the canopy transpiration resistance
            !! \latexonly 
            !!     \input{diffucointer3.tex}
            !! \endlatexonly
            frac_evap_intercept_predict(ji,jv) = veget(ji,jv) * zqsvegrap
            resist_intercept(ji,jv) = raero(ji) + rstruct(ji,jv) 
            vbeta_intercept(ji,jv) = veget(ji,jv) * zqsvegrap * (un / (un + speed * q_cdrag(ji) * rstruct(ji,jv)))
            
            !! We here determine the limitation of interception loss by the water stored on the leaf. 
            !! A first approximation of interception loss is obtained using the old values of
            !! qair and qsol_sat, which are functions of temp-sol and pb. (see call of 'qsatcalc')
            !! \latexonly 
            !!     \input{diffucointer4.tex}
            !! \endlatexonly
            ziltest = frac_evap_intercept_predict(ji,jv) * dt_sechiba * rau(ji) * (qsatt(ji) - qair(ji)) / resist_intercept(ji,jv)

            frac_evap_intercept(ji,jv) = frac_evap_intercept_predict(ji,jv)

            IF ( ziltest .GT. min_sechiba ) THEN

                !! \latexonly 
                !!     \input{diffucointer5.tex}
                !! \endlatexonly
                zrapp = qsintveg(ji,jv) / ziltest
                IF ( zrapp .LT. un ) THEN
                   
                    !! \latexonly 
                    !!     \input{diffucointer6.tex}
                    !! \endlatexonly
                    !!
		    !! Once the whole water stored on foliage has evaporated, transpiration can take place on the fraction
                    !! 'zqsvegrap'.
                   IF ( humrel(ji,jv) >= min_sechiba ) THEN
                      vbeta23(ji,jv) = MAX(vbeta_intercept(ji,jv) - vbeta_intercept(ji,jv) * zrapp, zero)
                   ELSE
                      ! We don't want transpiration when the soil cannot deliver it
                      vbeta23(ji,jv) = zero
                   ENDIF
                    
                    !! \latexonly 
                    !!     \input{diffucointer7.tex}
                    !! \endlatexonly
                    vbeta_intercept(ji,jv) = vbeta_intercept(ji,jv) * zrapp
                    frac_evap_intercept_corr(ji,jv) = frac_evap_intercept_predict(ji,jv)*zrapp
                    frac_evap_intercept(ji,jv) = frac_evap_intercept_corr(ji,jv)
                ENDIF
            ENDIF
        END IF
!        ! Autre formulation possible pour l'evaporation permettant une transpiration sur tout le feuillage
!        !commenter si formulation Nathalie sinon Tristan
!        speed = MAX(min_wind, wind(ji))
!        
!        vbeta23(ji,jv) = MAX(zero, veget(ji,jv) * (un / (un + speed * q_cdrag(ji) * rstruct(ji,jv))) - vbeta2(ji,jv))

      END DO

    END DO

    IF (printlev>=3) WRITE (numout,*) ' diffuco_inter done '

  END SUBROUTINE diffuco_inter


!! ==============================================================================================================================
!! SUBROUTINE      : diffuco_bare
!!
!>\BRIEF	   This routine computes the partial beta coefficient corresponding to
!! bare soil
!!
!! DESCRIPTION	   : Bare soil evaporation is submitted to a maximum possible flow (evap_bare_lim)
!! 
!! Calculation of wind speed
!! \latexonly 
!!     \input{diffucobare1.tex}
!! \endlatexonly
!!             
!! The calculation of $\beta_4$
!! \latexonly 
!!     \input{diffucobare2.tex}
!! \endlatexonly
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::frac_evap_bare
!!
!! REFERENCE(S)	 :
!! - de Noblet-Ducoudré, N, Laval, K & Perrier, A, 1993. SECHIBA, a new set of parameterisations
!! of the hydrologic exchanges at the land-atmosphere interface within the LMD Atmospheric General
!! Circulation Model. Journal of Climate, 6, pp.248-273
!! - Guimberteau, M, 2010. Modélisation de l'hydrologie continentale et influences de l'irrigation
!! sur le cycle de l'eau, PhD Thesis, available from:
!! http://www.sisyphe.upmc.fr/~guimberteau/docs/manuscrit_these.pdf
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE diffuco_bare (kjpindex, tot_bare_soil, veget_max, frac_evap_intercept, resist_intercept, frac_evap_transp, resist_transp, raero, &
                           evap_bare_lim, evap_bare_lim_ns, frac_evap_bare)

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                         :: kjpindex            !! Domain size (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)       :: tot_bare_soil       !! Total evaporating bare soil fraction
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: veget_max           !! Max. fraction of vegetation type (LAI->infty)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout)  :: frac_evap_intercept !! Fraction of gridcell where intercepted water evaporation occurs (-)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: resist_intercept    !! Sum of resistances applied to interceped water evaporation flux (s m^{-1})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout)  :: frac_evap_transp !! Fraction of gridcell where transpiration occurs (-)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: resist_transp       !! Sum of resistances applied to the transpiration flux (s m^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: raero               !! Aerodynamic resistance (s m^{-1})

    !! 0.2 Output variables
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)     :: frac_evap_bare      !! Fraction of gridcell where bare soil evaporation occurs (-)
    
    !! 0.3 Modified variables 
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)   :: evap_bare_lim           !! limiting factor for bare soil evaporation
                                                                                  !! when the 11-layer hydrology is used (-)
    REAL(r_std),DIMENSION (kjpindex,nstm), INTENT (inout) :: evap_bare_lim_ns     !! limiting factor for bare soil evaporation 
                                                                                  !! when the 11-layer hydrology is used (-)       
    !! 0.4 Local variables
    REAL(r_std),DIMENSION (kjpindex,nvm)               :: vbeta_intercept         !! Beta for Interception (for easier checking that the sum of the evaporative fluxes is not exceeding the total potential evaporation)
    REAL(r_std),DIMENSION (kjpindex,nvm)               :: vbeta_transp            !! Beta for Transpiration (for easier checking that the sum of the evaporative fluxes is not exceeding the total potential evaporation)
    REAL(r_std)                                        :: test                    !! Test for checking that sum of transpiration, IL and bare soil evaporation does not exceed potential evaporation
    INTEGER(i_std)                                     :: ji
    REAL(r_std), DIMENSION(kjpindex)                   :: vegtot
!_ ================================================================================================================================

    IF (printlev>=3) WRITE (numout,*) 'Entering diffuco_bare'

    !! 1. Calculation of the soil resistance and the beta (beta_4) for bare soil

    ! To use the new version of hydrol_split_soil and ensure the water conservation, we must have throughout sechiba at all time:
    !        a) evap_bare_lim(ji) = SUM(evap_bare_lim_ns(ji,:)*soiltile(ji,:)*vegtot(ji))
    !        b) all the terms (vbeta4, evap_bare_lim, evap_bare_lim_ns) =0 if evap_bare_lim(ji) LE min_sechiba
    ! This must also be kept true in diffuco_comb
    

    DO ji = 1, kjpindex
       vbeta_transp(ji,:) = frac_evap_transp(ji,:)*raero(ji)/resist_transp(ji,:)
       vbeta_intercept(ji,:) = frac_evap_intercept(ji,:)*raero(ji)/resist_intercept(ji,:)
       

       test = SUM(vbeta_intercept(ji,:)+vbeta_transp(ji,:))
       IF ( test .GT. un) THEN

          IF (printlev .GE. 3) WRITE(numout,*) "WARNING : BETA2 + BETA3 greater than 1, each reduced by: ", (test-1)*100, " % for pixel", ji

          !! Sending result in %
          test_overPET_intensity(ji)=(test-1)*100
          test_overPET_freq(ji)=test_overPET_freq(ji)+1

          !! Decreasing both fluxes (transpiration and IL by test in order to not exceed PET)
          frac_evap_transp(ji,:) = frac_evap_transp(ji,:)/test
          frac_evap_intercept(ji,:) = frac_evap_intercept(ji,:)/test

          !! Redefining the local vbetas (used a bit after in the routine)
          vbeta_transp(ji,:) = frac_evap_transp(ji,:)*raero(ji)/resist_transp(ji,:)
          vbeta_intercept(ji,:) = frac_evap_intercept(ji,:)*raero(ji)/resist_intercept(ji,:)
       ENDIF


       ! The limitation by 1-beta2-beta3 is due to the fact that evaporation under vegetation is possible
       !! \latexonly 
       !!     \input{diffucobare3.tex}
       !! \endlatexonly
       
       IF ( (evap_bare_lim(ji) .GT. min_sechiba) .AND. &  
            ! in this case we can't have vegtot LE min_sechina, cf hydrol_soil
            (un - SUM(vbeta_intercept(ji,:)+vbeta_transp(ji,:)) .GT. min_sechiba) ) THEN 
          ! eventually, none of the left-hand term is close to zero

          vegtot(ji) = SUM(veget_max(ji,:))
          IF (evap_bare_lim(ji) < (un - SUM(vbeta_intercept(ji,:)+vbeta_transp(ji,:)))) THEN
             ! Standard case

             frac_evap_bare(ji) = evap_bare_lim(ji)
          ELSE

             frac_evap_bare(ji) = un - SUM(vbeta_intercept(ji,:)+vbeta_transp(ji,:))

             ! We now have to redefine evap_bare_lim & evap_bare_lim_ns
             IF (evap_bare_lim(ji) .GT. min_sechiba) THEN

                evap_bare_lim_ns(ji,:) = evap_bare_lim_ns(ji,:) * frac_evap_bare(ji) / evap_bare_lim(ji)

             ELSE ! we must re-invent evap_bare_lim_ns => uniform across soiltiles             
                evap_bare_lim_ns(ji,:) = tot_bare_soil(ji)/vegtot(ji)
             ENDIF

             evap_bare_lim(ji) = frac_evap_bare(ji)
             ! consistent with evap_bare_lim(ji) =
             ! SUM(evap_bare_lim_ns(ji,:)*soiltile(ji,:)*vegtot(ji))
             ! as SUM(soiltile(ji,:)) = 1
          END IF
          
       ELSE ! instead of having very small vbeta4, we set everything to zero

          frac_evap_bare(ji) = zero
          evap_bare_lim(ji) = zero
          evap_bare_lim_ns(ji,:) = zero
       ENDIF
       
    END DO
    
    IF (printlev>=3) WRITE (numout,*) ' diffuco_bare done '
    
  END SUBROUTINE diffuco_bare


!! ==============================================================================================================================
!! SUBROUTINE   : diffuco_trans_co2
!!
!>\BRIEF        This subroutine computes carbon assimilation and stomatal 
!! conductance, following respectively Farqhuar et al. (1980) and Ball et al. (1987) 
!! and Yin X and Struik PC. 2009.
!!
!! DESCRIPTION  :\n
!! *** General:\n 
!! The equations are different depending on the photosynthesis mode (C3 versus C4).
!! Based on the work of Farquhar, von Caemmerer and Berry who published in 1980 a
!! biochemical model for C3 photosynthetic rates (the FvCB model), Yin and
!! Struik 2009 propose an analytical algorithm that incorporates a gs model and uses gm as a
!! temperature-dependent parameter for calculating assimilation based on the FvCB
!! model for various environmental scenarios. Yin and Struik 2009 also propose a
!! C4-equivalent version of the FvCB model and an associated analytical algorithm
!! for this model. They believe that applying FvCB-type models to both C3 and C4
!! crops is recommended to accurately predict the response of crop photosynthesis
!! to multiple, interactive environmental variables.  \n
!!
!!  Besides an analytical solution for photosynthesis the scheme includes a
!! modified Arrhenius function for the temperature dependence that accounts for a
!! decrease of Vcmax and Jmax at high temperatures and a temperature dependent
!! J max/Vcmax ratio (Kattge & Knorr 2007). The temperature response of Vcmax and
!! J max was parametrized with values from reanalysed data in literature (Kattge &
!! Knorr 2007), whereas Vcmax and Jmax at a reference temperature of 25° C were
!! derived from observed species-specific values in the TRY database (Kattge et
!! al. 2011).
!!   
!! In this version here we replace the old trunk assimilation and conductance 
!! which were computed over 20 levels of LAI and then integrated at the canopy level. 
!! These artificial LAI levels  were introduced in order to allow for a saturation of
!! photosynthesis with LAI, but they have no physical meaning. The approach of
!  the two stream radiation transfer model was extended for a multi-layer canopy.
!! The user can choose the numbers of layers in the run.def and the available light at
!! each layers decides about the numbers of layers used for the photosynthesis
!! calculation.
!!
!! NOTE: The variable Light_Abs_Tot is the amount of light abosrbed by the canopy,
!!       according to the two stream albedo routine.  Notice that this is a simple mean
!!       between the the absorption profile for light coming from an isotropic source 
!!       (e.g. clouds, and what light is scattered by aerosols on a sunny day) and the
!!       incoming radiation directly from the sun.  If the fraction of light from each
!!       source is ever provided from the atmospheric model or the forcing file, that
!!       will have to be changed in the albedo routines.
!!
!!       Also note that I feel sunlight and shaded leaf absorption cannot be calculated
!!       by using Collim_Abs_Tot and Isotrop_Abs_Tot from the albedo routines, since
!!       direct sunlight will scatter light through the canopy.  This light can hit
!!       shaded leaves.  Therefore, part of Collim_Abs_Tot will be used for sunlight
!!       leaves and part for shaded leaves, while all of Isotrop_Abs_Tot will be for
!!       shaded leaves.
!!
!! RECENT CHANGE(S): 
!!
!!   2018/12 added a nitrogen use efficiency constraint on maximum rate of carboxylation
!!
!!   2018/11 replaced t2m and q2m by temp_air and qair; prognostic variables on first atmospheric model layer
!!
!! MAIN OUTPUT VARIABLE(S): beta coefficients, resistances, CO2 intercellular 
!! concentration
!!
!! REFERENCE(S) :
!!
!!     Yin X, Struik PC. C3 and C4 photosynthesis models: an overview from the
!! perspective of crop modelling. NJAS-Wageningen Journal of Life Sciences 2009.
!!     C.J. Bernacchi, A.R. Portis, H. Nakano, S. von Caemmerer, S.P. Long,
!! Temperature response of mesophyll conductance. Implication for the determination
!! of Rubisco enzyme kinetics and for limitations to photosynthesis in vivo, Plant
!! Physiology 130 (2002) 1992–1998.
!!     C.J. Bernacchi, E.L. Singsaas, C. Pimentel, A.R. Portis Jr., S.P. Long,
!! Improved temperature response functions for models of Rubisco-limited
!! photosynthesis, Plant, Cell and Environment 24 (2001) 253–259.
!!     B.E. Medlyn, E. Dreyer, D. Ellsworth, M. Forstreuter, P.C. Harley, M.U.F.
!! Kirschbaum, X. Le Roux, P. Montpied, J. Strassemeyer, A. Walcroft, K. Wang, D.
!! Loustau, Temperature response of parameters of a biochemically based model of
!! photosynthesis. II. A review of experimental data, Plant, Cell and Environ- ment
!! 25 (2002) 1167–1179.
!!     Kattge, J. and Knorr, W.: Temperature acclimation in a biochemical model of
!! photosynthesis: a reanalysis of data from 36 species, Plant Cell Environ., 30
!! (2007) 1176–1190.
!!     Keenan, T., et al., Soil water stress and coupled photosynthesis–conductance
!! models: Bridging the gap between conflicting reports on the relative roles of
!! stomatal, mesophyll conductance and biochemical limitations to photosynthesis.
!! Agric. Forest Meteorol. (2010) 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

SUBROUTINE diffuco_trans_co2 (kjpindex, swdown, pb, qair, temp_air, &
                                temp_growth, rau, q_cdrag, humrel, &
                                assim_param, Ca, lai, &
                                veget_max, qsintveg, qsintmax, frac_evap_transp, resist_transp, vbeta3pot, rveget, rstruct, &
                                cimean, gsmean, gpp, cresist, vbeta23, Light_Abs_Tot,&
                                Light_Tran_Tot, lai_per_level, gs_distribution, gs_diffuco_output, &
                                VPD, air_relhum, f_gco2_out, gstot_component, gstot_frac, veget, warnings, u_speed, profile_vbeta3, &
                                profile_rveget, index_assi, leaf_ci, assimtot, assimi_lev, &
                                hist_id, indexveg, indexlai, index, kjit, cim, info_limitphoto, &
                                JJ_out,psi_leaf_next, raero, Rdtot_pft)


    !
    !! 0. Variable and parameter declaration
    !

    !
    !! 0.1 Input variables
    !
    INTEGER(i_std), INTENT(in)                               :: kjpindex         !! Domain size (unitless)
    REAL(r_std),DIMENSION (:), INTENT (in)                   :: swdown           !! Downwelling short wave flux 
                                                                                 !! @tex ($W m^{-2}$) @endtex 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: pb               !! Lowest level pressure (hPa)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: qair             !! Specific humidity at first atmospheric model layer
                                                                                 !! @tex ($kg kg^{-1}$) @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: temp_air         !! Air temperature at first atmospheric model layer (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: temp_growth      !! Growth temperature (°C) - Is equal to t2m_month
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: rau              !! air density @tex ($kg m^{-3}$) @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: q_cdrag          !! Surface drag coefficient (-)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)        :: humrel           !! Soil moisture stress (0-1,unitless)
    REAL(r_std),DIMENSION (:,:,:), INTENT (in)               :: assim_param      !! min+max+opt temps (K), vcmax, vjmax for 
                                                                                 !! photosynthesis 
                                                                                 !! @tex ($\mu mol m^{-2} s^{-1}$) @endtex 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: Ca               !! CO2 concentration inside the canopy

    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(in)         :: lai 
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)        :: veget_max        !! Maximum vegetation fraction of each PFT inside 
                                                                                 !! the grid box (0-1, unitless) 
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)        :: qsintveg         !! Water on vegetation due to interception 
                                                                                 !! @tex ($kg m^{-2}$) @endte
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)        :: qsintmax         !! Maximum water on vegetation
                                                                                 !! @tex ($kg m^{-2}$) @endtex
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)        :: vbeta23          !! Beta for fraction of wetted foliage that will 
                                                                                 !! transpire (unitless)
    REAL(r_std),DIMENSION (:,:,:),INTENT (in)                :: Light_Abs_Tot    !! Absorbed radiation per level for photosynthesis
    REAL(r_std),DIMENSION (:,:,:),INTENT (in)                :: Light_Tran_Tot   !! Transmitted radiation per level for photosynthesis
    REAL(r_std),DIMENSION(:,:,:),INTENT (in)                 :: lai_per_level    !! This is the LAI per vertical level
                                                                                 !! @tex $(m^{2} m^{-2})$
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT (in)        :: veget
    REAL(r_std), DIMENSION(:,:), INTENT (in)                   :: u_speed        !! wind speed at specified canopy level (m/s)   
    
    INTEGER(i_std),INTENT(in)                                :: hist_id          !! History_ file identifier (-)
    INTEGER(i_std),DIMENSION (kjpindex*nvm),INTENT(in)       :: indexveg         !! Indices of the points on the 3D map (-)    
    INTEGER(i_std),DIMENSION (kjpindex*(nlevels_tot+1)),INTENT(in)  :: indexlai  !! Indices of the points on the 3D map 
    INTEGER(i_std),DIMENSION (kjpindex),INTENT(in)           :: index            !! Indices of the points on the map (-) 
    INTEGER(i_std),INTENT(in)                                :: kjit             !! Time step number (-) 

    REAL(r_std),DIMENSION(kjpindex,nvm), INTENT(in)          :: psi_leaf_next    !! Leaf water potential at time step n-1 (MPa)
    REAL(r_std),DIMENSION(kjpindex), INTENT(in)              :: raero    !! Leaf water potential at time step n-1 (MPa)

    !
    !! 0.2 Output variables
    !
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)       :: frac_evap_transp !! Fraction of gridcell where transpiration occurs (-)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)       :: resist_transp    !! Sum of the resistances applied to the transpiration flux (s m^{-1})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)       :: vbeta3pot        !! Beta for Potential Transpiration
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)       :: rveget           !! stomatal resistance of vegetation 
                                                                                 !! @tex ($s m^{-1}$) @endtex
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)       :: rstruct          !! structural resistance @tex ($s m^{-1}$) @endtex
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)       :: cimean           !! mean intercellular CO2 concentration 
                                                                                 !! @tex ($\mu mol mol^{-1}$) @endtex
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)       :: gsmean           !! mean stomatal conductance to CO2 (umol m-2 s-1)
    REAL(r_Std),DIMENSION (kjpindex,nvm), INTENT (out)       :: gpp              !! Assimilation ((gC m^{-2} dt_sechiba^{-1}), total area)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)       :: cresist          !! coefficient for resistances (??)

    REAL(r_Std),DIMENSION (kjpindex,nvm,nlevels_tot), INTENT (out) &
                                                             :: gs_distribution    
    REAL(r_Std),DIMENSION (kjpindex,nvm,nlevels_tot), INTENT (out) &
                                                             :: gs_diffuco_output
    REAL(r_Std),DIMENSION (kjpindex,nvm), INTENT (out) &
                                                             :: f_gco2_out
    REAL(r_Std),DIMENSION (kjpindex,nvm,nlevels_tot), INTENT (out) &       
                                                             :: gstot_component
    REAL(r_Std),DIMENSION (kjpindex,nvm,nlevels_tot), INTENT (out) &       
                                                             :: gstot_frac
    REAL(r_std),DIMENSION (kjpindex,nvm,nlevels_tot), INTENT (out) &       
                                                             :: profile_vbeta3
    REAL(r_std),DIMENSION (kjpindex,nvm,nlevels_tot), INTENT (out) &
                                                             :: profile_rveget
    INTEGER(i_Std),DIMENSION (kjpindex), INTENT (out)        :: index_assi       !! pixels for a given PFT where photosynthesis takes place
    REAL(r_std),DIMENSION (kjpindex,nvm,nlevels_tot), INTENT (out)  :: leaf_ci          !! intercellular CO2 concentration (ppm)
    REAL(r_std),DIMENSION(kjpindex,nvm), INTENT (out)        :: assimtot         !! total assimilation 
    REAL(r_std),DIMENSION(kjpindex,nvm), INTENT (out)        :: Rdtot_pft        !! total respiration 
    REAL(r_std),DIMENSION(kjpindex,nvm,nlevels_tot), INTENT (out) &
                                                             :: assimi_lev       !! assimilation (at all LAI levels)
                                                                                 !! @tex ($\mu mol m^{-2} s^{-1}$) @endtex

    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)       :: cim              !! Intercellular CO2 over nlevels_tot 
    REAL(r_std),DIMENSION (kjpindex,nvm,nlevels_tot+1),INTENT(out) &
                                                             :: info_limitphoto  !! Store information of limitation for photosynthesis

    REAL(r_std),DIMENSION (kjpindex,nvm,nlevels_tot),INTENT(out) &
                                                             :: JJ_out           !! Store JJ
    !
    !! 0.3 Modified variables
    REAL(r_std),DIMENSION(kjpindex,nvm,nwarns), INTENT(inout) &
                                                             :: warnings         !! A warning counter

    REAL(r_std), DIMENSION(kjpindex),INTENT(inout)           :: air_relhum       !! air relative humidity at 2m
                                                                                 !! @tex ($kg kg^{-1}$) @endtex
    REAL(r_std), DIMENSION(kjpindex),INTENT(inout)           :: VPD              !! Vapor Pressure Deficit (kPa)
 
    !
    !! 0.4 Local variables
    !
    
    REAL(r_std), DIMENSION(kjpindex,nvm,nlevels_tot+1)       :: assimi           !! assimilation (at a specific LAI level)
                                                                                 !! @tex ($\mu mol m^{-2} s^{-1}$) @endtex
    REAL(r_std), DIMENSION(kjpindex,nvm)                     :: coupled_lai
    !
    REAL(r_Std),DIMENSION (kjpindex,nvm,nlevels_tot)         :: gs_store
    REAL(r_Std),DIMENSION (kjpindex,nvm)                     :: gs_column_total
    REAL(r_Std),DIMENSION (kjpindex,nvm)                     :: gstot_store


    REAL(r_std),DIMENSION (kjpindex,nvm,nlevels_tot)         :: profile_gs_store
    REAL(r_std),DIMENSION (kjpindex,nvm,nlevels_tot)         :: profile_cimean
    REAL(r_std),DIMENSION (kjpindex,nlevels_tot)             :: profile_assimtot

    REAL(r_std),DIMENSION (kjpindex,nlevels_tot)             :: profile_Rdtot
    REAL(r_std),DIMENSION (kjpindex,nlevels_tot)             :: profile_gstot
    REAL(r_std),DIMENSION (kjpindex,nlevels_tot)             :: profile_gstop

    REAL(r_std),DIMENSION (kjpindex,nvm,nlevels_tot)         :: profile_gsmean

    REAL(r_std),DIMENSION (kjpindex,nvm,nlevels_tot)         :: profile_gpp

    
    REAL(r_std),DIMENSION (kjpindex,nvm,nlevels_tot)         :: profile_rstruct
    REAL(r_std),DIMENSION (nlevels_tot)                      :: profile_speed
    REAL(r_std),DIMENSION (nlevels_tot)                      :: profile_cresist
    REAL(r_std),DIMENSION(kjpindex,nvm)                      :: resist_stom    !! Leaf water potential at time step n-1 (MPa)
    
    !
    !
    REAL(r_std),DIMENSION(kjpindex,nvm,nlevels_tot)          :: Abs_light        !! fraction of absorbed light within each layer (just to check)
    REAL(r_std),DIMENSION(kjpindex,nvm)                      :: Abs_light_can    !! fraction of absorbed light within the canopy (just to check)  
    REAL(r_std),DIMENSION (kjpindex,nvm)                     :: vcmax            !! maximum rate of carboxylation Notice that,
                                                                                 !! as it is currently written, this is only
                                                                                 !! a function of nue and leafN_top
                                                                                 !! @tex ($\mu mol CO2 m^{-2} s^{-1}$) @endtex
    REAL(r_std)                           :: N_profile                           !! N-content reduction factor for a specific canopy depth (unitless) 
    REAL(r_std)                           :: N_normalisation                     !! Normalisation factor for the N-profile, 
                                                                                 !! to get the correct total N content (unitless) 
    REAL(r_std),DIMENSION (kjpindex)      :: leafN_top                           !! Leaf N content - top of canopy (gN m-2[leaf]) 
    INTEGER(i_std)                        :: ji, jv, ilev, jl, limit_photo       !! indices (unitless)
    REAL(r_std), DIMENSION(kjpindex)      :: leaf_ci_lowest                      !! intercellular CO2 concentration at the lowest 
                                                                                 !! LAI level
                                                                                 !! @tex ($\mu mol mol^{-1}$) @endtex
    REAL(r_std), DIMENSION(kjpindex)                         :: zqsvegrap        !! relative water quantity in the water 
                                                                                 !! interception reservoir (0-1,unitless) 
    REAL(r_std)                                              :: speed            !! wind speed @tex ($m s^{-1}$) @endtex
    ! Assimilation
    LOGICAL, DIMENSION(kjpindex)                             :: assimilate       !! where assimilation is to be calculated 
                                                                                 !! (unitless) 
    INTEGER(i_std)                                           :: nia,inia,nina    !! counter/indices (unitless)
    INTEGER(i_std)                                           :: inina,iainia     !! counter/indices (unitless)
    INTEGER(i_std), DIMENSION(kjpindex)                      :: index_non_assi   !! indices (unitless)
    REAL(r_std), DIMENSION(kjpindex, nlevels_tot+1)          :: vc2              !! rate of carboxylation (at a specific LAI level) 
                                                                                 !! @tex ($\mu mol CO2 m^{-2} s^{-1}$) @endtex 

    REAL(r_std), DIMENSION(kjpindex,nlevels_tot+1)           :: vj2              !! rate of Rubisco regeneration (at a specific LAI 
                                                                                 !! level) @tex ($\mu mol e- m^{-2} s^{-1}$) @endtex 
    REAL(r_std), DIMENSION(kjpindex)                         :: gstop            !! stomatal conductance to H2O at topmost level 
                                                                                 !! @tex ($m s^{-1}$) @endtex

    REAL(r_std), DIMENSION(kjpindex,nlevels_tot+1)           :: gs               !! stomatal conductance to CO2 
                                                                                 !! @tex ($\mol m^{-2} s^{-1}$) @endtex
    REAL(r_std), DIMENSION(kjpindex)                         :: gamma_star       !! CO2 compensation point (ppm)
                                                                                 !! @tex ($\mu mol mol^{-1}$) @endtex
    REAL(r_std), DIMENSION(kjpindex)                         :: water_lim        !! water limitation factor (0-1,unitless)

    REAL(r_std), DIMENSION(kjpindex)                         :: gstot            !! total stomatal conductance to H2O
                                                                                 !! Final unit is
                                                                                 !! @tex ($m s^{-1}$) @endtex
    REAL(r_std), DIMENSION(kjpindex)                         :: Rdtot            !! Total Day respiration (respiratory CO2 release other than by photorespiration) (mumol CO2 m−2 s−1) 
    REAL(r_std), DIMENSION(kjpindex)                         :: gs_top           !! leaf stomatal conductance to H2O across all top levels
                                                                                 !! @tex ($\mol H2O m^{-2} s^{-1}$) @endtex
    REAL(r_std), DIMENSION(kjpindex)                         :: qsatt                               !! surface saturated humidity at 2m (??) 
                                                                                 !! @tex ($g g^{-1}$) @endtex
                                                                                 !! levels (0-1,unitless)
    REAL(r_std), DIMENSION(kjpindex)                         :: T_Vcmax          !! Temperature dependance of Vcmax (unitless)
    REAL(r_std), DIMENSION(kjpindex)                         :: S_Vcmax_acclim_temp                 !! Entropy term for Vcmax 
                                                                                 !! accounting for acclimation to temperature (J K-1 mol-1)
    REAL(r_std), DIMENSION(kjpindex)                         :: T_Jmax           !! Temperature dependance of Jmax
    REAL(r_std), DIMENSION(kjpindex)                         :: S_Jmax_acclim_temp                  !! Entropy term for Jmax 
                                                                                 !! accounting for acclimation toxs temperature (J K-1 mol-1)
    REAL(r_std), DIMENSION(kjpindex)                         :: T_gm             !! Temperature dependance of gmw
    REAL(r_std), DIMENSION(kjpindex)                         :: T_Rd             !! Temperature dependance of Rd (unitless)
    REAL(r_std), DIMENSION(kjpindex)                         :: T_Kmc            !! Temperature dependance of KmC (unitless)
    REAL(r_std), DIMENSION(kjpindex)                         :: T_KmO            !! Temperature dependance of KmO (unitless)
    REAL(r_std), DIMENSION(kjpindex)                         :: T_Sco            !! Temperature dependance of Sco
    REAL(r_std), DIMENSION(kjpindex)                         :: T_gamma_star     !! Temperature dependance of gamma_star (unitless)    
    REAL(r_std), DIMENSION(kjpindex)                         :: vc               !! Maximum rate of Rubisco activity-limited carboxylation (mumol CO2 m−2 s−1)
    REAL(r_std), DIMENSION(kjpindex)                         :: vj               !! Maximum rate of e- transport under saturated light (mumol CO2 m−2 s−1)
    REAL(r_std), DIMENSION(kjpindex)                         :: gm               !! Mesophyll diffusion conductance (molCO2 mâ2 sâ1 barâ1) 
    REAL(r_std), DIMENSION(kjpindex)                         :: g0var 
    REAL(r_std), DIMENSION(kjpindex,nlevels_tot+1)           :: Rd               !! Day respiration (respiratory CO2 release other than by photorespiration)(mumol CO2 mâ2 sâ1)
    REAL(r_std), DIMENSION(kjpindex)                         :: Kmc              !! Michaelis–Menten constant of Rubisco for CO2 (mubar)
    REAL(r_std), DIMENSION(kjpindex)                         :: KmO              !! Michaelis–Menten constant of Rubisco for O2 (mubar)
    REAL(r_std), DIMENSION(kjpindex)                         :: Sco              !! Relative CO2 /O2 specificity factor for Rubisco (bar bar-1) 
    REAL(r_std), DIMENSION(kjpindex)                         :: gb_co2           !! Boundary-layer conductance (molCO2 mâ2 sâ1 barâ1) 
    REAL(r_std), DIMENSION(kjpindex)                         :: gb_h2o           !! Boundary-layer conductance (molH2O mâ2 sâ1 barâ1)
    REAL(r_std), DIMENSION(kjpindex)                         :: fvpd             !! Factor for describing the effect of leaf-to-air vapour difference on gs (-)
    REAL(r_std), DIMENSION(kjpindex,nvm)                     :: f_psi            !! Sentivity of stomata to leaf water potential
    REAL(r_std), DIMENSION(kjpindex)                         :: f_gco2           !! Sentivity function that controls the stomatal conductance (either fvpd or f_psi according to the resolution wanted)
    REAL(r_std), DIMENSION(kjpindex,nvm)                     :: psi_leaf_approx  !! Smooth value of Leaf water potential on timepsi time-steps (MPa)
    REAL(r_std), DIMENSION(kjpindex)                         :: low_gamma_star   !! Half of the reciprocal of Sc/o (bar bar-1)
    REAL(r_std)                                              :: fcyc             !! Fraction of electrons at PSI that follow cyclic transport around PSI (-)
    REAL(r_std)                                              :: z                !! A lumped parameter (see Yin et al. 2009) ( mol mol-1)                          
    REAL(r_std)                                              :: Rm               !! Day respiration in the mesophyll (umol CO2 m−2 s−1)
    REAL(r_std)                                              :: Cs_star          !! Cs -based CO2 compensation point in the absence of Rd (ubar)
    REAL(r_std), DIMENSION(kjpindex)                         :: Iabs             !! Photon flux density absorbed by leaf photosynthetic pigments (umol photon m−2 s−1)
    REAL(r_std), DIMENSION(kjpindex)                         :: Jmax             !! Maximum value of J under saturated light (umol e− m−2 s−1)
    REAL(r_std), DIMENSION(kjpindex,nlevels_tot+1)           :: JJ               !! Rate of e− transport (umol e− m−2 s−1)
    REAL(r_std)                                              :: J2               !! Rate of all e− transport through PSII (umol e− m−2 s−1)
    REAL(r_std)                                              :: VpJ2             !! e− transport-limited PEP carboxylation rate (umol CO2 m−2 s−1)
    REAL(r_std)                                              :: A_1, A_3         !! Lowest First and third roots of the analytical solution for a general 
                                                                                 !! cubic equation (see Appendix A of Yin et al. 2009) (umol CO2 m−2 s−1)
    REAL(r_std)                                              :: A_1_tmp, A_3_tmp            !! Temporary First and third roots of the analytical solution 
                                                                                            !! for a general cubic equation (see Appendix A of Yin et al. 2009) (umol CO2 m−2 s−1)
    REAL(r_std)                                              :: Obs                         !! Bundle-sheath oxygen partial pressure (ubar)
    REAL(r_std), DIMENSION(kjpindex,nlevels_tot+1)           :: Cc_loc                      !! Chloroplast CO2 partial pressure (ubar)
                                                                                            !! Notice we had _loc to the name to avoid
                                                                                            !! conflict with a variable already declared
                                                                                            !! in constantes_var
    REAL(r_std)                                              :: ci_star                     !! Ci -based CO2 compensation point in the absence of Rd (ubar)        
    REAL(r_std)                                              :: a,b,c,d,m,f,j,g,h,i,l,p,q,r !! Variables used for solving the cubic equation (see Yin et al. (2009))
    REAL(r_std)                                              :: QQ,UU,PSI,x1,x2,x3          !! Variables used for solving the cubic equation (see Yin et al. (2009))
    rEAL(r_std)                                              :: lai_min                     !! mininum lai to calculate photosynthesis       
    REAL(r_std), DIMENSION(kjpindex)                         :: laisum                      !! when calculating cim over nlevels_tot

    LOGICAL,DIMENSION(kjpindex,nlevels_tot)                  :: top_level                   !! A flag to see if a given level is considered
                                                                                            !! a "top" level or not.  Must be done this way
                                                                                            !! because there could be several "top" levels.
    REAL(r_std)                                              :: totlai_top                  !! The sum of the LAI in all the levels that we
                                                                                            !! declare as the "top".
    INTEGER(i_std)                                           :: lev_count                   !! How many levels we have already taken for
                                                                                            !! the "top"

! @defgroup Photosynthesis Photosynthesis
! @{   
    ! 1. Preliminary calculations\n
    REAL(r_std),DIMENSION (kjpindex,nvm,nlevels_tot)         :: light_tran_to_level         !! Cumulative amount of light transmitted to a given level
    REAL(r_std), DIMENSION(kjpindex)                         :: gab_h2o                     !! Boundary-layer + aerodynamic conductance (molH2O m-2 s-1 bar-1)
    LOGICAL,DIMENSION (kjpindex)                             :: lskip                       !! A flag to indicate we've hit an unusual
                                                                                            !! case were we may have numerical instability.
    INTEGER                                                  :: ic, jc, kc


    INTEGER                                                  :: ipts
    LOGICAL                                                  :: llight
    LOGICAL                                                  :: llai
    INTEGER,SAVE                                             :: istep=0
!$OMP THREADPRIVATE(istep)

!_ ================================================================================================================================

    leaf_ci(:,:,:) = zero

   IF (printlev>=4) WRITE (numout,*) ' diffuco_trans_co2'
   IF (printlev_loc>=5) WRITE(numout,*) 'diffuco_trans_co2: nlevels_tot,lai_per_level',nlevels_tot,lai_per_level(:,:,:)

    profile_gs_store = 0.0d0
    profile_cimean = 0.0d0
    profile_assimtot = 0.0d0

    profile_Rdtot = 0.0d0
    profile_gstot = 0.0d0
    profile_gstop = 0.0d0

    profile_gsmean = 0.0d0
    profile_gpp = 0.0d0
    profile_rveget = 2.5d20
    profile_rstruct = 0.0d0
    profile_speed = 0.0d0
    profile_cresist = 0.0d0
    profile_vbeta3 = 0.0d0

    JJ_out = 0.0d0
    
    ! Debug 
    IF(printlev_loc>=4)THEN 
       WRITE(numout,*) 'Absorbed light in diffuco. ',test_grid,test_pft 
       DO ilev=nlevels_tot,1,-1 
          WRITE(numout,'(I5,10F20.10)') & 
               ilev, Light_Abs_Tot(test_grid,test_pft,ilev),& 
               Light_Tran_Tot(test_grid,test_pft,ilev),& 
               lai_per_level(test_grid,test_pft,ilev) 
       ENDDO
    ENDIF
    !-

    gs_column_total = 0.0d0
    gs_distribution = 0.0d0
    gs_diffuco_output = 0.0d0
    gstot_frac = 0.0d0
    f_gco2_out = 0.0d0
    gstot_component = 0.d0
    gstot_store = 0.0d0
    gs_store = 0.0d0
    !
    ! @addtogroup Photosynthesis
    ! @{   
    ! 1.3 Estimate relative humidity of air (for calculation of the stomatal conductance).\n
    !! \latexonly
    !! \input{diffuco_trans_co2_1.3.tex}
    !! \endlatexonly
    ! @}
    !
    CALL qsatcalc (kjpindex, temp_air, pb, qsatt)
    air_relhum(:) = &
      ( qair(:) * pb(:) / (Tetens_1+qair(:)* Tetens_2) ) / &
      ( qsatt(:)*pb(:) / (Tetens_1+qsatt(:)*Tetens_2 ) )

    ! WRITE(numout, *) '050315 air_relhum is: ', air_relhum

    VPD(:) = ( qsatt(:)*pb(:) / (Tetens_1+qsatt(:)*Tetens_2 ) ) &
         - ( qair(:) * pb(:) / (Tetens_1+qair(:)* Tetens_2) )
    ! We limit the impact of VPD in the range of [0:6] kPa
    ! This seems to be the typical range of values
    ! that the vapor pressure difference can take, although
    ! there is a paper by Mott and Parkhurst, Plant, Cell & Environment,
    ! Vol. 14, pp. 509--519 (1991) which gives a value of 15 kPa.
    !!$ VPD(:) = MAX(zero,MIN(VPD(:)/10.,6.))
    ! Convert VPD from hPa to kPa
    VPD(:) = VPD(:)/10.
    

    !
    ! 2. beta coefficient for vegetation transpiration
    !
    rstruct(:,1) = rstruct_const(1)
    !+++CHECK+++
    ! Setting rveget makes the model crash when running in debug mode becaus 
    ! the default xios settings cannot write undefined values.
!!$    rveget(:,:) = undef_sechiba
    rveget(:,:) = un
    !+++++++++++
    frac_evap_transp(:,:) = zero
    resist_transp(:,:) = un
    vbeta3pot(:,:) = zero
    gsmean(:,:) = zero
    gpp(:,:) = zero
    light_tran_to_level(:,:,:)=un
    assimi_lev(:,:,:) = zero
    cim(:,:) = zero
    cimean(:,1) = Ca(:)
    info_limitphoto(:,:,:) = zero 
    f_psi(:,:) = zero
    !
    ! 2.1 Initializations
    ! Loop over vegetation types
    ! 

    DO jv = 2,nvm
       gamma_star(:) = zero 
       Kmo(:) = zero 
       Kmc(:) = zero 
       gm(:) = zero 
       g0var(:) = zero 
       
       Cc_loc(:,:) = zero 
       vc2(:,:) = zero
       vj2(:,:) = zero 
       JJ(:,:) = zero 
       gs(:,:) = zero 
       assimi(:,:,:) = zero 
       Rd(:,:) = zero
         
      !
      ! beta coefficient for vegetation transpiration
      !
      rstruct(:,jv) = rstruct_const(jv)
      cimean(:,jv) = Ca(:)

    ! ---TEMP---
    ! calculation of vcmax before N is included  
    !  IF (downregulation_co2) THEN
    !     vcmax(:,jv) = assim_param(:,jv,ivcmax)*&
    !          (un-downregulation_co2_coeff(jv)*&
    !          log(Ca(:)/downregulation_co2_baselevel))
    !  ELSE
    !     vcmax(:,jv) = assim_param(:,jv,ivcmax)
    !  ENDIF
    !-----------

      WHERE (lai(:,jv) .GT. min_sechiba)
         !The normalisation factor is now computed later, based on
         !the profile of the transmitted light.
         !leafN_top(:) = assim_param(:,jv,ileafN) * ext_coeff_N(jv) / &
         !     ( 1. -exp(-ext_coeff_N(jv) * lai(:,jv)) ) 
         leafN_top(:) = assim_param(:,jv,ileafN) 

      ELSEWHERE

         leafN_top(:) = zero

      ENDWHERE

      ! Very import line of code. Here vcmas is CALCULATED from nue and ileafN
      ! which are both stored in assim_param. Note that the calculation of
      ! vcmax in stomate_vmax is no longer used in the model.
      ! Note here leaf nitrogen per ground area is used (by integrating all leaves),
      ! thus vcmax is calculated based on ground area as well. The vcmax will be
      ! converted to be based per leaf area per LAI layer later on.
      vcmax(:,jv) = assim_param(:,jv,inue)*leafN_top(:)

      !WRITE(numout,*) '210515 assim_param(grid,pft,npco2)','grid:'grid,'pft:',pft,'npco2:',npco2
      !WRITE(numout,*)  assim_param(:,:,:)      
      !! mask that contains points where there is photosynthesis
      !! For the sake of vectorisation [DISPENSABLE], computations are done only for convenient points.
      !! nia is the number of points where the assimilation is calculated and nina the number of points where photosynthesis is not
      !! calculated (based on criteria on minimum or maximum values on LAI, vegetation fraction, shortwave incoming radiation, 
      !! temperature and relative humidity).
      !! For the points where assimilation is not calculated, variables are initialized to specific values. 
      !! The assimilate(kjpindex) array contains the logical value (TRUE/FALSE) relative to this photosynthesis calculation.
      !! The index_assi(kjpindex) array indexes the nia points with assimilation, whereas the index_no_assi(kjpindex) array indexes
      !! the nina points with no assimilation.
      !          
      ! set minimum lai value to test if photosynthesis should be calculated
      ! 0.0001 is abitrary
      ! could be externalised 
      lai_min = 0.0001
      !
      nia=0
      nina=0
      index_assi=0
      index_non_assi=0
      !
      DO ji=1,kjpindex     

         IF(veget_max(ji,jv) == zero)THEN
            ! this vegetation type is not present, so no reason to do the 
            ! calculation. However, WHERE loops later on depend on the
            ! value of assimilate, so we set it here.
            assimilate(ji) = .FALSE.
            CYCLE
         ENDIF

         !
         ! This is tricky.  We need a certain amount of LAI to have photosynthesis.
         ! However, we divide our canopy up into certain levels, so even if we
         ! have a total amount that is sufficient, there might not be enough
         ! in any particular layer.  So we check each layer to see if there
         ! is one that has enough to give a non-zero photosynthesis.  In
         ! the worst-case scenario (small LAI in every other level), this becomes 
         ! our top level for the gs calculation.  If gs_top_level is zero below,
         ! we have a divide by zero issue.
         llai=.FALSE.
         DO ilev=1,nlevels_tot
            IF(lai_per_level(ji,jv,ilev) .GT. lai_min) THEN 
             llai=.TRUE.
            ELSE
            ENDIF
         ENDDO


         IF ( llai .AND. &
              ( veget_max(ji,jv) .GT. min_sechiba ) ) THEN
            
            ! Here we need to do something different than in the previous model.
            ! We explicitly need to check for absorbed light.  This occurs
            ! because our albedo in the two stream model is a function
            ! of the solar angle.  If the sun has set, we have no
            ! light, and therefore we cannot have photosynthesis.
            llight=.FALSE.
            DO ilev=1,nlevels_tot
               IF(Light_Abs_Tot(ji,jv,ilev) .GT. min_sechiba) llight=.TRUE.
            ENDDO
            IF ( ( swdown(ji) .GT. min_sechiba )   .AND. &
                 ( humrel(ji,jv) .GT. min_sechiba) .AND. &
                 ( temp_growth(ji) .GT. tphoto_min(jv) ) .AND. &
                 llight .AND. &
                 ( temp_growth(ji) .LT. tphoto_max(jv) )) THEN
               !
               assimilate(ji) = .TRUE.
               nia=nia+1
               index_assi(nia)=ji
               !
            ELSE
               !
               assimilate(ji) = .FALSE.
               nina=nina+1
               index_non_assi(nina)=ji
               !
            ENDIF
         ELSE
            !
            assimilate(ji) = .FALSE.
            nina=nina+1
            index_non_assi(nina)=ji
            !
         ENDIF
        ! 
      ENDDO 
      !


      IF(printlev_loc>=4 .AND. jv == test_pft)THEN
         WRITE(numout,*) 'jv: ',jv
         WRITE(numout,*) 'test_grid: ',test_grid
         WRITE(numout,*) 'llai, veget_max: ',llai,veget_max(:,jv)
         WRITE(numout,*) 'temp_growth: ',temp_growth(test_grid)
         WRITE(numout,*) 'tphoto_min: ',tphoto_min(jv)
         WRITE(numout,*) 'tphoto_max: ',tphoto_max(jv)
         WRITE(numout,*) 'swdown: ',swdown(test_grid)
         WRITE(numout,*) 'index_non_assi(:): ',index_non_assi(test_grid)
         WRITE(numout,*) 'index_assi(:): ',index_assi(test_grid)
         WRITE(numout,*) 'nia: ',nia
         WRITE(numout,*) 'nina',nina
      ENDIF

      top_level(:,:)=.FALSE.
      gstot(:) = zero
      gstop(:) = zero
      assimtot(:,jv) = zero
      Rdtot_pft(:,jv) = zero
      Rdtot(:)=zero
      gs_top(:) = zero
      lskip(:)=.FALSE.
      !
      zqsvegrap(:) = zero
      WHERE (qsintmax(:,jv) .GT. min_sechiba)
         !! relative water quantity in the water interception reservoir
         zqsvegrap(:) = MAX(zero, qsintveg(:,jv) / qsintmax(:,jv))
      ENDWHERE

      !! Calculate the water limitation factor that restricts Vcmax when hydrol_arch is turned off.
      !!  If hydrol_architecture is used, VCmax is not controled but instead there is a control on
      !!  stomatal conductance
      IF (.NOT. ok_hydrol_arch) THEN
         WHERE ( assimilate(:) )
            water_lim(:) =  humrel(:,jv)
         ELSEWHERE
            water_lim(:) = zero
         ENDWHERE
      END IF

      ! give a default value of ci for all pixel that do not assimilate
      DO ilev=1,nlevels_tot
         DO inina=1,nina
            leaf_ci(index_non_assi(inina),jv,ilev) = Ca(index_non_assi(inina))
         ENDDO
      ENDDO 

      !
      !
      ! Here is the calculation of assimilation and stomatal conductance
      ! based on the work of Farquahr, von Caemmerer and Berry (FvCB model) 
      ! as described in Yin et al. 2009
      ! Yin et al. developed a extended version of the FvCB model for C4 plants
      ! and proposed an analytical solution for both photosynthesis pathways (C3 and C4)
      ! Photosynthetic parameters used are those reported in Yin et al. 
      ! Except For Vcmax25, relationships between Vcmax25 and Jmax25 for which we use 
      ! Medlyn et al. (2002) and Kattge & Knorr (2007)
      ! Because these 2 references do not consider mesophyll conductance, we neglect this term
      ! in the formulations developed by Yin et al. 
      ! Consequently, gm (the mesophyll conductance) tends to the infinite
      ! This is of importance because as stated by Kattge & Knorr and Medlyn et al.,
      ! values of Vcmax and Jmax derived with different model parametrizations are not 
      ! directly comparable and the published values of Vcmax and Jmax had to be standardized
      ! to one consistent formulation and parametrization

      ! See eq. 6 of Yin et al. (2009)
      ! Parametrization of Medlyn et al. (2002) - from Bernacchi et al. (2001)
      T_KmC(:)        = Arrhenius(kjpindex,temp_air,298.,E_KmC(jv))
      T_KmO(:)        = Arrhenius(kjpindex,temp_air,298.,E_KmO(jv))
      T_Sco(:)        = Arrhenius(kjpindex,temp_air,298.,E_Sco(jv))
      T_gamma_star(:) = Arrhenius(kjpindex,temp_air,298.,E_gamma_star(jv))


      ! Parametrization of Yin et al. (2009) - from Bernacchi et al. (2001)
      T_Rd(:)         = Arrhenius(kjpindex,temp_air,298.,E_Rd(jv))


      ! For C3 plants, we assume that the Entropy term for Vcmax and Jmax 
      ! acclimates to temperature as shown by Kattge & Knorr (2007) - Eq. 9 and 10
      ! and that Jmax and Vcmax respond to temperature following a modified Arrhenius function
      ! (with a decrease of these parameters for high temperature) as in Medlyn et al. (2002) 
      ! and Kattge & Knorr (2007).
      ! In Yin et al. (2009), temperature dependance to Vcmax is based only on a Arrhenius function
      ! Concerning this apparent unconsistency, have a look to the section 'Limitation of 
      ! Photosynthesis by gm' of Bernacchi (2002) that may provide an explanation
      
      ! Growth temperature tested by Kattge & Knorr range from 11 to 35°C
      ! So, we limit the relationship between these lower and upper limits
      S_Jmax_acclim_temp(:) = aSJ(jv) + bSJ(jv) * MAX(11., MIN(temp_growth(:),35.))      
      T_Jmax(:)  = Arrhenius_modified(kjpindex,temp_air,298.,E_Jmax(jv),D_Jmax(jv),S_Jmax_acclim_temp)

      S_Vcmax_acclim_temp(:) = aSV(jv) + bSV(jv) * MAX(11., MIN(temp_growth(:),35.))   
      T_Vcmax(:) = Arrhenius_modified(kjpindex,temp_air,298.,E_Vcmax(jv),D_Vcmax(jv),S_Vcmax_acclim_temp)

      vc(:) = vcmax(:,jv) * T_Vcmax(:)
  

      ! As shown by Kattge & Knorr (2007), we make use
      ! of Jmax25/Vcmax25 ratio (rJV) that acclimates to temperature for C3 plants
      ! rJV is written as a function of the growth temperature
      ! rJV = arJV + brJV * T_month 
      ! See eq. 10 of Kattge & Knorr (2007)
      ! and Table 3 for Values of arJV anf brJV 
      ! Growth temperature is monthly temperature (expressed in °C) - See first paragraph of
      ! section Methods/Data of Kattge & Knorr
      vj(:) = ( arJV(jv) + brJV(jv) *  MAX(11., MIN(temp_growth(:),35.)) ) * vcmax(:,jv) * T_Jmax(:)

      T_gm(:)  = Arrhenius_modified(kjpindex,temp_air,298.,E_gm(jv),D_gm(jv),S_gm(jv))

      IF (ok_hydrol_arch) THEN
         !! If the hydraulic architecture, by Tuzet et al. (2017), is activated, the water_lim factor is 
         !! calculated and applied later.
         gm(:) = gm25(jv) * T_gm(:)
      ELSE
         gm(:) = gm25(jv) * T_gm(:) * MAX(1-stress_gm(jv), water_lim(:))
      ENDIF

      ! @endcodeinc
      !
      KmC(:)=KmC25(jv)*T_KmC(:)
      KmO(:)=KmO25(jv)*T_KmO(:)
      Sco(:)=Sco25(jv)*T_sco(:)
      gamma_star(:) = gamma_star25(jv)*T_gamma_star(:)
      ! low_gamma_star is defined by Yin et al. (2009)
      ! as the half of the reciprocal of Sco - See Table 2
      !!$ We derive its value from Equation 7 of Yin et al. (2009)
      !!$ assuming a constant O2 concentration (Oi)
      !!$ low_gamma_star(:) = gamma_star(:)/Oi
      low_gamma_star(:) = 0.5 / Sco(:)

      ! VPD expressed in kPa
      ! Note : MIN(1.-min_sechiba,MAX(min_sechiba,(a1(jv) - b1(jv) * VPD(:))))
      ! is always between 0-1 not including 0 and 1


      !! According to the way the hydraulic architecture is solved, the stomatal conductance is either linked to the VPD or the
      !! leaf water status (thanks to the leaf water potential). In this configuration, as the system is fully coupled, it 
      !! resolution can lead to instabilities when a stress begins to be felt. To avoid this instability, a relaxation function 
      !! is implemented. The aim of this relaxation is to smooth the oscillations resulting from the instability.

     
      IF (ok_hydrol_arch) THEN

         !! If the hydraulic architecture is activated, the calculation of f_gco2 now depends on the 
         !! leaf water potential as in Tuzet et al. (2003) (Eq.1 and 2)
         
         !! Then, the f_psi function from Tuzet et al. (2003) is computed and multiplied by the a_fpsi coefficient. The value 
         !! value is given to the variable f_gco2.

         IF ( (printlev_loc >=4).AND.(jv==test_pft) ) THEN   
            WRITE(numout,*) 'psi_leaf_next', psi_leaf_next(test_grid,jv)
            CALL flush(numout)   
         ENDIF
         f_psi(:,jv) = MAX(min_sechiba, a_fpsi(jv) * (1. + EXP(sf(jv) * psi_ref_g0(jv) ) )/ &
                           (1. + EXP(sf(jv) * (psi_ref_g0(jv) - psi_leaf_next(:,jv)))))
         
         f_gco2(:) = f_psi(:,jv)

         !! For this hydraulic architecture, the water limitation factor is the same that controls stomatal conductance.
         ! Non Stomatal Limitations which proxy is the function of VPD below is
         ! applied to water_lim which will constrain mesophyll conductance and
         ! maximum carboxylation rate

         IF (ok_nsl_vpd) THEN
            WHERE ( assimilate(:) )
               water_lim(:) = f_gco2(:)/(1+VPD(:)/d_l)/a_fpsi(jv)
            ELSEWHERE 
               water_lim(:) = zero
            ENDWHERE
         ELSE
            WHERE ( assimilate(:) )
               water_lim(:) = f_gco2(:)/a_fpsi(jv)
            ELSEWHERE 
               water_lim(:) = zero
            ENDWHERE
         ENDIF

         gm(:) = gm(:) * MAX(1-stress_gm(jv), water_lim(:))

      ELSE

         !! If no hydraulic architecture is called, the previous definition of fvpd is computed.

         fvpd(:) = 1. / ( 1. / MIN(1.-min_sechiba,MAX(min_sechiba,(a1(jv) - b1(jv) * VPD(:)))) - 1. ) & 
              * MAX(1-stress_gs(jv), water_lim(:))

         f_gco2(:) = fvpd(:)

      ENDIF

      !! To avoid problems of two strong residual transpiration, g0 is limitated by the f_gco2 function. This "patch" is not really 
      !! appropriate and maybe highlights an oversestimation of g0.
      IF (no_g0) THEN 
         g0var(:) = g0(jv) * MAX(1-stress_gm(jv), water_lim(:))
      ELSE 
         g0var(:) = g0(jv) 
      ENDIF

      ! We do not convert the boundary layer conductance from m/s to
      ! molH2O/m2/s anymore as we now convert the conductance equivalent
      ! to the aerodynamic and boundary layer conductances later
      gb_h2o(:) = gb_ref

      ! Addition of a aerodynamic conductance in the computation of GPP
      ! Conversion from a conductance in (m s-1) to (mol H2O m-2 s-1)
      ! (44.6 * (tp_00/temp_air(:) * (pb(:)/pb_std)) from Peary et al. 1991     

      IF (rough_dyn) THEN
         ! If ROUGH_DYN is activated, then the drag coefficient
         ! includes the aerodynamic and boundary layer conductances.
         ! Then we do not consider the boundary layer conductance
         ! computed from Pearcy et al. 1991.
         gab_h2o(:) = MAX(min_wind, wind(:)) * q_cdrag(:) *&
              44.6 * (tp_00/temp_air(:)) * (pb(:)/pb_std)
      ELSE
         ! If ROUGH_DYN is not activated, then the drag coefficient
         ! includes only the aerodynamic conductance. Then we add
         ! the boundary layer conductance from Pearcy et al 1991 !
         ! to the aerodynamic conductance.
         gab_h2o(:) = (gb_h2o(:) * MAX(min_wind, wind(:)) * q_cdrag(:) *&
              44.6 * (tp_00/temp_air(:)) * (pb(:)/pb_std) ) / &
              (gb_h2o(:) + (MAX(min_wind, wind(:)) * q_cdrag(:) ) )

      END IF

      ! conversion from (mol H2O m-2 s-1) to (mol CO2 m-2 s-1)
      gb_co2(:) = gab_h2o(:) / ratio_H2O_to_CO2

      !
      ! 2.4 Loop over LAI levels to estimate assimilation and conductance        
      !
      JJ(:,:)=zero
      vc2(:,:)=zero
      vj2(:,:)=zero
      Cc_loc(:,:)=zero
      gs(:,:)=zero
      Rd(:,:)=zero


      IF (nia .GT. 0) THEN

         assim_loop:DO inia=1,nia

            !
            iainia=index_assi(inia)

            !
            ! Find the top layer that still contains lai.  We have also
            ! added a check in here to make sure this layer has absorbed
            ! light.  If there is no absorbed light, we will not calculate
            ! stomatal conductance for the top level, which means gstop
            ! will be zero later and we'll crash when we divide by it.
            ! Include a warning so we know how often this happens.
            ! The old system of chosing just one top level caused some bad results
            ! in the transpir, evidently due to the calculation of vbeta3 being
            ! too low due to too low of LAI in the top level.  Therefore, we are
            ! now going to take several top levels in an effort to make sure
            ! we have enough LAI.  We aim for an arbitrary amount of 0.1 or
            ! the top three levels.  This has not been rigoursly tested, but it
            ! is theoretically justified by the fact that if the LAI of the top
            ! level is very low, other leaves will be getting a lot of sun and
            ! well-exposed to the atmosphere, and therefore contribute equal
            ! to the resistence of the vegetation.
            totlai_top=zero
            lev_count=0
            !   
            DO ilev=nlevels_tot, 1, -1
               ! lai_min is 0.0001. If a layer with so little lai occurs at the
               ! bottom of the canopy, the absorbed light can be below 10e-8.
               ! This implies that we should test wether Light_Abs_Tot exceeds
               ! zero rather than 10-8 (min_sechiba). It is difficult to set a
               ! strict threshold above which Light_Abs_Tot should be if
               ! the lai_per_level is above lai_min. Setting the threshold to
               ! zero is save but not very useful as it will only track the 
               ! worst possible errors.
               IF(lai_per_level(iainia,jv,ilev) .GT. lai_min)THEN
                  IF(Light_Abs_Tot(iainia,jv,ilev) .GE. zero)THEN
                     IF(totlai_top .LT. lai_top(jv) .AND. lev_count .LT. nlev_top)THEN
                        lev_count=lev_count+1
                        top_level(iainia,ilev)=.TRUE.
                        totlai_top=totlai_top+lai_per_level(iainia,jv,ilev)
                     ELSE
                        ! We have enough levels now.
                        !WRITE(numout,*) 'Top levels ',totlai_top,lev_count
                        EXIT
                     ENDIF
                  ELSE
                     WRITE(numout,*) 'ERROR: Found a level with LAI, ' // &
                          'but no absorbed light.'
                     WRITE(numout,*) 'ERROR: iainia,jv,ilev, ',iainia,jv,ilev
                     WRITE(numout,*) 'ERROR: Light_Abs_Tot(iainia,jv,ilev),' //&
                          'lai_per_level(iainia,jv,ilev),',iainia,jv, &
                          ilev,Light_Abs_Tot(iainia,jv,ilev), &
                          lai_per_level(iainia,jv,ilev)
                     CALL ipslerr_p(3,'diffuco_trans_co2', &
                          'Found a level with LAI but no absorbed light','','')
                  ENDIF
               ENDIF
            ENDDO
            IF(totlai_top .LT. min_sechiba)THEN
               ! There seem to be some fringe cases where this happens, usually
               ! at sunrise or sunset with low LAI levels.  Instead of stopping
               ! here, I'm going to issue a warning.  If this case happens,
               ! we produce no GPP, but we're not going to crash.  We need to
               ! be careful in one of the loops below, since a stomatal 
               ! conductance of zero will cause a crash.
                IF(err_act.GT.1)THEN
                   WRITE(numout,*) 'WARNING: Did not find a level with LAI in diffuco!'
                   WRITE(numout,*) 'WARNING: iainia,ivm: ',iainia,jv
                ENDIF
               assimtot(iainia,jv) = zero
               Rdtot_pft(iainia,jv) = zero
               Rdtot(iainia) = zero
               gstot(iainia) = zero
               assimi_lev = zero
               lskip(iainia)=.TRUE.
               CYCLE assim_loop
            ENDIF

            ! We need to know how much is transmitted to this level
            ! from the top of the canopy, since this absolute number
            ! is what is actually used in the photosyntheis.  Right now
            ! Light_Tran_Tot and Light_Abs_Tot are the 
            ! relative quantities transmitted through and absorbed
            ! by the layer.
            DO ilev = nlevels_tot-1,1,-1

               light_tran_to_level(iainia,jv,ilev) = &
                    light_tran_to_level(iainia,jv,ilev+1) * &
                    Light_Tran_Tot(iainia,jv,ilev+1)

            ENDDO

            !We sum N_profile over the levels where it will be used in the next loop below.
            !The total will give the normalisation factor.
            N_normalisation=0
            DO ilev = 1, nlevels_tot
               IF (Light_Abs_Tot(iainia,jv,ilev) .GT. min_sechiba) THEN
                  IF(lai_per_level(iainia,jv,ilev) .LE. min_sechiba)THEN
                     CYCLE
                  ENDIF
                  N_profile = ( un - .7_r_std * &
                       ( un - light_tran_to_level(iainia,jv,ilev) ) )
                  N_normalisation = N_normalisation+N_profile*lai_per_level(iainia,jv,ilev)
               ENDIF
            ENDDO

            ! Debug
            IF(printlev_loc>=4 .AND. jv == test_pft .AND. iainia == test_grid)THEN
               WRITE(numout,*) 'Absolute transmitted light to each level',&
                    jv,iainia
               DO ilev = nlevels_tot,1,-1
                  WRITE(numout,*) ilev,light_tran_to_level(iainia,jv,ilev)

               ENDDO
            ENDIF
            !-

            DO ilev = 1, nlevels_tot
 
               IF (Light_Abs_Tot(iainia,jv,ilev) .GT. min_sechiba) THEN

                  IF(lai_per_level(iainia,jv,ilev) .LE. min_sechiba)THEN

                     ! This is not necessarily a problem.  LAI is updated in 
                     ! slowproc, while absorbed light is updated in condveg, 
                     ! which is called after diffuco. So we might be one 
                     ! timestep off.  It generally should be at night,
                     ! though, so there is no absorbed light. If this 
                     ! persists, that is a real problem.
                     IF(err_act.GT.1)THEN
                        WRITE(numout,*) 'WARNING: We have absorbed light but no LAI ' // &
                             'in this level. Skipping.'
                        WRITE(numout,'(A,3I8)') 'WARNING: iainia,jv,ilev: ',iainia,jv,ilev
                        WRITE(numout,'(A,2F20.10)') 'WARNING: Light_Abs_Tot(iainia,jv,ilev), ' // &
                             'lai_per_level(iainia,jv,ilev): ',&
                             Light_Abs_Tot(iainia,jv,ilev), &
                             lai_per_level(iainia,jv,ilev)
                     ENDIF
                     CYCLE
                  ENDIF
                  
                  ! nitrogen is scaled according to the transmitted light
                  ! before the caclulation was: N_profile = 
                  ! exp( -ext_coeff_N(jv)*laitab(jl) ). Where laitab is the 
                  ! lai contained in a fixed lai layer. ext_coeff_N is 
                  ! prescribed so this was a prescribed profile. We replaced
                  ! it by the transmitted light calculated by the 2stream
                  ! model, maybe ( un - Light_Tran_Tot(iainia,jv,ilev) )
                  ! could be replaced by the directly calculated absorbed
                  ! light Light_Abs_Tot(iainia,jv,ilev).
                  ! Notice that light in the old code was the cumulative 
                  ! light transmitted to the layer, while Light_Tran_Tot 
                  ! is the relative light transmitted through the layer.
                  ! Note the trick to convert vcmax from per ground area to 
                  ! per leaf area is made here and below. N_profile calculated here 
                  ! and used below, is implicitly divided by the leaf area index 
                  ! (leaf area per ground area). c.f. how N_normalisation is calculated
                  ! above.
                  N_profile = ( un - .7_r_std * &
                       ( un - light_tran_to_level(iainia,jv,ilev) ) )/N_normalisation
                 
                  ! the vcmax and photosynthesis is only restricted by 
                  ! stomatal closure when transpiration is higher than the 
                  ! water supply to the leaves as calculated in hydrol_arch. 
                  ! This more consistent with the substantial consensus in 
                  ! plant physiology that the main restriction of 
                  ! photosynthesis during drought occurs through diffusion 
                  ! limitation rather than photochemical or biochemical failure
                  ! (Flexas et al. 2006). Effects of drought stress on 
                  ! respiration are less well known: increases, decreases and 
                  ! no change have been found. For now there is no effect of 
                  ! drought stress on dark respiration when hydrol_arch is 
                  ! used.
                  vc2(iainia,ilev) = vc(iainia) * N_profile * MAX(1-stress_vcmax(jv), water_lim(iainia))
                  vj2(iainia,ilev) = vj(iainia) * N_profile * MAX(1-stress_vcmax(jv), water_lim(iainia))
                     
                  ! see Comment in legend of Fig. 6 of Yin et al. (2009)
                  ! Rd25 is assumed to equal 0.01 Vcmax25 
                  Rd(iainia,ilev) = vcmax(iainia,jv) * N_profile * 0.01 * &
                       T_Rd(iainia) * MAX(1-stress_vcmax(jv), water_lim(iainia))
                  ! As Light_Abs_Tot is the absorbed radiation per layer, 
                  ! expressed as the fraction of overall light hitting the 
                  ! canopy, there is no need for an additional multiplication 
                  ! by light_tran_to_level, which is is the cumulative 
                  ! fraction of transmitted light (otherwise we are counting
                  ! twice the extinction by the above layers).
                  !Iabs(iainia) =swdown(iainia) * W_to_mol * RG_to_PAR *&
                  !     Light_Abs_Tot(iainia,jv,ilev) * &
                  !     light_tran_to_level(iainia,jv,ilev) / &
                  !     lai_per_level(iainia,jv,ilev)
                  Iabs(iainia) =swdown(iainia) * W_to_mol * RG_to_PAR *&
                       Light_Abs_Tot(iainia,jv,ilev) / &
                       lai_per_level(iainia,jv,ilev)


                  ! If we have a very dense canopy, the lower levels might
                  ! get almost no light.  If the light is close to zero, it
                  ! causes a numerical problem below, as x1 will be zero. Of
                  ! course, if the absorbed light is really small 
                  ! photosynthesis will not happen, so we can just skip this 
                  ! level. This number is fairly arbitrary, but it should be 
                  ! as small as possible.
                  IF(Iabs(iainia) .LT. 1e-12)CYCLE

                  ! eq. 4 of Yin et al (2009)
                  ! Note that Iabs is very sensitive to the number of layers 
                  ! use in the photosynthesis module and to the canopy 
                  ! structure. The first sensitivity is related to 
                  ! none-linearities to the light conditions and is an 
                  ! numerical issue. The sensitivity to the canopy structure
                  ! is a wanted behavior (that's why canopy structure was 
                  ! implemented) but makes a substantial difference for the 
                  ! total GPP in tropical compared to temperate forests. An 
                  ! easy way to prescribe a different canopy structure is by 
                  ! changing the range of the prescribed height. The parameter
                  ! ::alpha_LL describes the relationship between Jmax and
                  ! light saturation. It is thought to depend on the amount of 
                  ! chlorophyl in the leaves which in turn depends on the 
                  ! N-content. For the moment we use a prescribed value. After 
                  ! the introduction of the structured canopy we noticed an 
                  ! important sensitivity of GPP to this parameter which may 
                  ! indicate that a substantial part of our canopy is under
                  ! light saturation. 
                  Jmax(iainia)     = vj2(iainia,ilev)
                  JJ(iainia,ilev)  = ( alpha_LL(jv) * Iabs(iainia) + Jmax(iainia) - &
                       sqrt((alpha_LL(jv) * Iabs(iainia) + Jmax(iainia) )**2. &
                       - 4 * theta(jv) * Jmax(iainia) * alpha_LL(jv) * &
                       Iabs(iainia)) ) &
                       / ( 2 * theta(jv))

                  ! Store JJ for recalculation in mleb_diffuco_co2. To avoid expensive
                  ! calculations with analytic solutions and loops, JJ will be out.
                  JJ_out(iainia,jv,ilev) = JJ(iainia,ilev)                 
 
                 IF (printlev_loc>=4 .AND. test_pft == jv) THEN
                     WRITE(numout,*) '*************************'
                     WRITE(numout,*) 'jv,ilev: ',jv,ilev
                     WRITE(numout,*) 'JJ(iainia,ilev),Jmax(iainia),Iabs(iainia),Rd(iainia,ilev): ',&
                          JJ(iainia,ilev),Jmax(iainia),Iabs(iainia),Rd(iainia,ilev)
                     WRITE(numout,*) 'vj2(iainia,ilev),vc2(iainia,ilev), N_profile:',&
                          vj2(iainia,ilev),vc2(iainia,ilev),N_profile
                     WRITE(numout,*) 'Light_Abs_Tot(iainia,jv,ilev),Light_Tran_Tot(iainia,jv,ilev):',&
                          Light_Abs_Tot(iainia,jv,ilev), &
                          Light_Tran_Tot(iainia,jv,ilev)
                     
                     WRITE(numout,*) 'water_lim(iainia),T_Rd(iainia):',&
                          water_lim(iainia),T_Rd(iainia)
                     WRITE(numout,*) 'vc(iainia),vj(iainia), vcmax(iainia,jv):',&
                          vc(iainia),vj(iainia),vcmax(iainia,jv)
                     WRITE(numout,*) '*************************'
                 ENDIF

                 IF ( is_c4(jv) )  THEN 
                    !
                    ! @addtogroup Photosynthesis
                    ! @{   
                    !
                    ! 2.4.2 Assimilation for C4 plants (Collatz et al., 1992)\n
                    !! \latexonly
                    !! \input{diffuco_trans_co2_2.4.2.tex}
                    !! \endlatexonly
                    ! @}           
                    !
                    ! Analytical resolution of the Assimilation based 
                    ! Yin et al. (2009)
                    ! Eq. 28 of Yin et al. (2009)
                    fcyc= 1. - ( 4.*(1.-fpsir(jv))*(1.+fQ(jv)) + &
                         3.*h_protons(jv)*fpseudo(jv) ) / &
                         ( 3.*h_protons(jv) - 4.*(1.-fpsir(jv)))

                    ! See paragraph after eq. (20b) of Yin et al.
                    Rm=Rd(iainia,ilev)/2.

                    ! We assume that cs_star equals ci_star 
                    ! (see Comment in legend of Fig. 6 of Yin et al. (2009)
                    ! Equation 26 of Yin et al. (2009)
                    Cs_star = (gbs(jv) * low_gamma_star(iainia) * &
                         Oi - ( 1. + low_gamma_star(iainia) * alpha(jv) / &
                         0.047) * Rd(iainia,ilev) + Rm ) / ( gbs(jv) + kp(jv) ) 

                    ! eq. 11 of Yin et al (2009)
                    J2 = JJ(iainia,ilev) / ( 1. - fpseudo(jv) / ( 1. - fcyc ) )
                    
                    ! Equation right after eq. (20d) of Yin et al. (2009)
                    z = ( 2. + fQ(jv) - fcyc ) / ( h_protons(jv) * (1. - fcyc ))
                    
                    VpJ2 = fpsir(jv) * J2 * z / 2.
                    
                    A_3=9999.

                    ! See eq. right after eq. 18 of Yin et al. (2009)
                    ! The loop over limit_photo calculates the assimultion
                    ! rate of the Rubisco-limited CO2 assimilation (Ac) and
                    ! the rate of electron transport-limited CO2 assimultion
                    ! (Aj).  The overall assimilation will be the slower
                    ! of these two processes.
                    DO limit_photo=1,2
                       ! Is Vc limiting the Assimilation
                       IF ( limit_photo .EQ. 1 ) THEN
                          a = 1. + kp(jv) / gbs(jv)
                          b = 0.
                          x1 = vc2(iainia,ilev)
                          x2 = KmC(iainia)/KmO(iainia)
                          x3 = KmC(iainia)
                          ! Is J limiting the Assimilation
                       ELSE
                          a = 1.
                          b = VpJ2
                          x1 = (1.- fpsir(jv)) * J2 * z / 3.
                          x2 = 7. * low_gamma_star(iainia) / 3.
                          x3 = 0.
                       ENDIF

                       m=f_gco2(iainia)-g0var(iainia)/gb_co2(iainia)
                       d=g0var(iainia)*(Ca(iainia)-Cs_star) + f_gco2(iainia)*Rd(iainia,ilev)
                       f=(b-Rm-low_gamma_star(iainia)*Oi*gbs(jv))*x1*d + a*gbs(jv)*x1*Ca(iainia)*d
                       j=(b-Rm+gbs(jv)*x3 + x2*gbs(jv)*Oi)*m + (alpha(jv)*x2/0.047-1.)*d &
                            + a*gbs(jv)*(Ca(iainia)*m - d/gb_co2(iainia) - (Ca(iainia) - Cs_star ))

                       g=(b-Rm-low_gamma_star(iainia)*Oi*gbs(jv))*x1*m - (alpha(jv)*low_gamma_star(iainia)/0.047+1.)*x1*d &
                            + a*gbs(jv)*x1*(Ca(iainia)*m - d/gb_co2(iainia) - (Ca(iainia)-Cs_star ))

                       h=-((alpha(jv)*low_gamma_star(iainia)/0.047+1.)*x1*m + (a*gbs(jv)*x1*(m-1.))/gb_co2(iainia) )
                       i= ( b-Rm + gbs(jv)*x3 + x2*gbs(jv)*Oi )*d + a*gbs(jv)*Ca(iainia)*d
                       l= ( alpha(jv)*x2/0.047 - 1.)*m - (a*gbs(jv)*(m-1.))/gb_co2(iainia)
                       
                       p = (j-(h-l*Rd(iainia,ilev))) / l
                       q = (i+j*Rd(iainia,ilev)-g) / l
                       r = -(f-i*Rd(iainia,ilev)) / l 
                       
                       ! See Yin et al. (2009) and  Baldocchi (1994)
                       QQ = ( (p**2._r_std) - 3._r_std * q) / 9._r_std
                       UU = ( 2._r_std* (p**3._r_std) - 9._r_std *p*q + 27._r_std *r) /54._r_std

                       
                       ! +++ CHECK +++
                       ! The line below looked better however for quick trusting
                       ! passing we'll overlook it for now
                       ! IF  (QQ .GT. EPSILON(1._r_std)) THEN 
                       IF  (QQ .GE. 0._r_std) THEN
                          IF (ABS(UU/(QQ**1.5_r_std) ) .LE. 1._r_std)  THEN
                             PSI = ACOS(UU/(QQ**1.5_r_std))
                             A_3_tmp = -2._r_std * SQRT(QQ) * COS(( PSI + 4._r_std * PI)/3._r_std ) - p / 3._r_std
                             IF (( A_3_tmp .LT. A_3 )) THEN
                                A_3 = A_3_tmp
                                info_limitphoto(iainia,jv,ilev)=2.
                             ELSE
                                ! In case, J is not limiting the assimilation
                                ! we have to re-initialise a, b, x1, x2 and x3 values
                                ! in agreement with a Vc-limited assimilation 
                                a = 1. + kp(jv) / gbs(jv)
                                b = 0.
                                x1 = vc2(iainia,ilev)
                                x2 = KmC(iainia)/KmO(iainia)
                                x3 = KmC(iainia)
                                info_limitphoto(iainia,jv,ilev)=1.
                             ENDIF
                          ENDIF
                          ! Used in post-proc for counting the proportion of errors
                          IF ((printlev_loc >= 4).AND.(jv==test_pft).AND.(iainia==test_grid)) THEN
                             WRITE(numout,*) 'No pb with QQ (diffuco.f90)'
                          ENDIF
                       ELSE
                          IF ((printlev_loc >= 4).AND.(jv==test_pft).AND.(iainia==test_grid)) THEN
                             WRITE(numout,*) 'Problem with QQ (diffuco.f90)'
                          ENDIF
                          WRITE(numout,*) 'WARNING: unexpected value for QQ (diffuco.f90)'
                       ENDIF

                       IF ( ( A_3 .EQ. 9999. ) .OR. ( A_3 .LT. (-Rd(iainia,ilev)) ) ) THEN
                          IF ( printlev>=4 ) THEN
                             WRITE(numout,*) 'We have a problem in diffuco_trans_co2'
                             WRITE(numout,*) 'no real positive solution found for pft:',jv
                             WRITE(numout,*) 'temp_air:',temp_air(iainia)
                             WRITE(numout,*) 'vpd:',VPD(iainia)
                             WRITE(numout,*) 'Rd:',Rd(iainia,ilev)
                          END IF
                          A_3 = -Rd(iainia,ilev)
                       ENDIF
                 
                       assimi(iainia,jv,ilev) = A_3

                       IF(printlev_loc>=5 .AND. jv == test_pft .AND. iainia == test_pft)THEN
                          WRITE(numout,*) '( x1 - ( assimi(iainia,jv,ilev) + Rd(iainia,ilev) ) ', &
                               ( x1 - (assimi(iainia,jv,ilev) + Rd(iainia,ilev) ))
                          WRITE(numout,*) 'x1,assimi(iainia,jv,ilev),Rd(iainia,ilev): ', &
                               x1,assimi(iainia,jv,ilev),Rd(iainia,ilev)
                          WRITE(numout,*) 'jv,ilev,iainia: ', jv,ilev,iainia
                       ENDIF

                       !++++ TEMP +++++
                       IF(( x1 - (assimi(iainia,jv,ilev) + Rd(iainia,ilev) )) == zero)THEN
                          ! This causes a problem, since we really should not
                          ! be dividing by zero.  However, it seems to happen
                          ! fairly rarely, so all we'll do is cycle to the
                          ! next loop and keep track of this warning so we
                          ! can look at it afterwards.
                          warnings(iainia,jv,iwphoto)=warnings(iainia,jv,iwphoto)+un

                          IF(err_act.GT.1)THEN
                             WRITE(numout,*) 'Diffuco, Getting ready to divide by zero! C4'
                             WRITE(numout,*) '( x1 - ( assimi(iainia,jv,ilev) + Rd(iainia,ilev) ) ', &
                                  ( x1 - (assimi(iainia,jv,ilev) + Rd(iainia,ilev) ))
                             WRITE(numout,*) 'x1,assimi(iainia,jv,ilev),Rd(iainia,ilev): ', &
                                  x1,assimi(iainia,jv,ilev),Rd(iainia,ilev)
                             WRITE(numout,*) 'jv,ilev,iainia: ', jv,ilev,iainia
                          ENDIF

                          CYCLE

                       ENDIF

                       IF ( ABS( assimi(iainia,jv,ilev) + Rd(iainia,ilev) ) .LT. min_sechiba ) THEN
                          gs(iainia,ilev) = g0var(iainia)
                          !leaf_ci keeps its initial value (Ca)
                       ELSE
                          ! Eq. 24 of Yin et al. (2009) 
                          Obs = ( alpha(jv) * assimi(iainia,jv,ilev) ) / ( 0.047 * gbs(jv) ) + Oi
                          ! Eq. 23 of Yin et al. (2009)
                          Cc_loc(iainia,ilev) = ( ( assimi(iainia,jv,ilev) + Rd(iainia,ilev) ) * &
                               ( x2 * Obs + x3 ) + low_gamma_star(iainia) * Obs * x1 ) &
                               / MAX(min_sechiba, x1 - ( assimi(iainia,jv,ilev) + Rd(iainia,ilev) ))
                          ! Eq. 22 of Yin et al. (2009)
                          leaf_ci(iainia,jv,ilev) = ( Cc_loc(iainia,ilev) - ( b - assimi(iainia,jv,ilev) - Rm ) / gbs(jv) ) / a
                                                    
                          !This is eq. 25 of Yin et al. (2009), adjusted
                          !followin Duursma 2019, but we want to
                          !use Ca instead of Cs. Therefore we need to take into
                          !account the CO2 transfer along the path from Ca to Cs,
                          !thus accounting for the boundary layer conductance (gb_co2).
                          IF (ok_hydrol_arch) THEN 
                             gs(iainia,ilev) = MAX(g0var(iainia),(assimi(iainia,jv,ilev) + Rd(iainia,ilev) ) / &
                                  ((Ca(iainia)-assimi(iainia,jv,ilev))/gb_co2(iainia)-Cs_star ) * f_gco2(iainia))
                          ELSE
                             gs(iainia,ilev) = g0var(iainia) + (assimi(iainia,jv,ilev) + Rd(iainia,ilev) ) / &
                                  ((Ca(iainia)-assimi(iainia,jv,ilev))/gb_co2(iainia)-Cs_star ) * f_gco2(iainia)                               
                          ENDIF
                       ENDIF
                       !                       
                    ENDDO !! ENDDO about limit photo
                    !
                 ELSE  
                    !
                    ! @addtogroup Photosynthesis
                    ! @{   
                    !
                    ! 2.4.3 Assimilation for C3 plants (Farqhuar et al., 1980)\n
                    !! \latexonly
                    !! \input{diffuco_trans_co2_2.4.3.tex}
                    !! \endlatexonly
                    ! @}           
                    !
                    !
                    !
                    A_1=9999.

                    ! See eq. right after eq. 18 of Yin et al. (2009)
                    DO limit_photo=1,2
                       ! Is Vc limiting the Assimilation
                       IF ( limit_photo .EQ. 1 ) THEN

		        IF (printlev_loc>=4) THEN
	                    WRITE(numout,*) 'IF - A_1, Do in limit_photo limit_photo .EG. 1 '
		        ENDIF
                          x1 = vc2(iainia,ilev)
                          ! It should be O not Oi (comment from Vuichard)
                          x2 = KmC(iainia) * ( 1. + 2*gamma_star(iainia)*Sco(iainia) / KmO(iainia) )
                          ! Is J limiting the Assimilation
                       ELSE
                          x1 = JJ(iainia,ilev)/4.
                          x2 = 2. * gamma_star(iainia)
                       ENDIF


                       ! See Appendix B of Yin et al. (2009)
                       a = g0var(iainia) * ( x2 + gamma_star(iainia) ) + &
                            ( g0var(iainia) / gm(iainia) + f_gco2(iainia) ) * ( x1 - Rd(iainia,ilev) )
                       b = Ca(iainia) * ( x1 - Rd(iainia,ilev) ) - gamma_star(iainia) * x1 - Rd(iainia,ilev) * x2
                       c = Ca(iainia) + x2 + ( 1./gm(iainia) + 1./gb_co2(iainia) ) * ( x1 - Rd(iainia,ilev) ) 
                       d = x2 + gamma_star(iainia) + ( x1 - Rd(iainia,ilev) ) / gm(iainia)
                       m = 1./gm(iainia) + ( g0var(iainia)/gm(iainia) + f_gco2(iainia) ) * ( 1./gm(iainia) + 1./gb_co2(iainia) )

                       p = - ( d + (x1 - Rd(iainia,ilev) ) / gm(iainia) + a * (1./gm(iainia) + 1./gb_co2(iainia) ) + &
                            ( g0var(iainia)/gm(iainia) + f_gco2(iainia)) * c) / m

                       q = ( d * ( x1 - Rd(iainia,ilev) ) + a*c + ( g0var(iainia)/gm(iainia) + f_gco2(iainia) ) * b ) / m
                       r = - a * b / m

                       ! See Yin et al. (2009) 
                       QQ = ( (p**2._r_std) - 3._r_std * q) / 9._r_std
                       UU = ( 2._r_std* (p**3._r_std) - 9._r_std *p*q + 27._r_std *r) /54._r_std

                       ! IF  (QQ .GT. EPSILON(1._r_std)) THEN
                       IF  (QQ .GE. 0._r_std) THEN
                          IF (ABS(UU/(QQ**1.5_r_std) ) .LE. 1._r_std)  THEN
                          PSI = ACOS(UU/(QQ**1.5_r_std))
                          A_1_tmp = -2._r_std * SQRT(QQ) * COS( PSI / 3._r_std ) - p / 3._r_std

                             IF (( A_1_tmp .LT. A_1 )) THEN
                                A_1 = A_1_tmp
                                info_limitphoto(iainia,jv,ilev)=2.
                             ELSE
                                ! In case, J is not limiting the assimilation
                                ! we have to re-initialise x1 and x2 values
                                ! in agreement with a Vc-limited assimilation 
                                x1 = vc2(iainia,ilev)
                                ! It should be O not Oi (comment from Vuichard)
                                x2 = KmC(iainia) * ( 1. + 2*gamma_star(iainia)*Sco(iainia) / KmO(iainia) )
                                info_limitphoto(iainia,jv,ilev)=1.
                             ENDIF
                          ENDIF 
                          ! Used in post-proc for counting the proportions of
                          ! errors
                          IF ((printlev_loc >= 4).AND.(jv==test_pft).AND.(iainia==test_grid)) THEN
                             WRITE(numout,*) 'No pb with QQ (diffuco.f90)'
                          ENDIF
                       ELSE
                          IF ((printlev_loc >= 4).AND.(jv==test_pft).AND.(iainia==test_grid)) THEN
                             WRITE(numout,*) 'Problem with QQ (diffuco.f90)'
                          ENDIF
                          WRITE(numout,*) 'WARNING: unexpected value for QQ (diffuco.f90)'
                       ENDIF 
                    ENDDO

                    IF ( (A_1 .EQ. 9999.) .OR. ( A_1 .LT. (-Rd(iainia,ilev)) ) ) THEN
                       IF ( printlev>=4 ) THEN
                          WRITE(numout,*) 'We have a problem in diffuco_trans_co2'
                          WRITE(numout,*) 'no real positive solution found for pft:',jv
                          WRITE(numout,*) 'temp_air:',temp_air(iainia)
                          WRITE(numout,*) 'vpd:',VPD(iainia)
                          WRITE(numout,*) 'Setting the solution to -Rd(iainia,ilev):',-Rd(iainia,ilev)
                       END IF
                       A_1 = -Rd(iainia,ilev)
                    ENDIF
                    assimi(iainia,jv,ilev) = A_1

                    ! Eq. 18 of Yin et al. (2009)

                    IF(printlev_loc>=5 .AND. jv == test_pft .AND. iainia == test_pft)THEN
                       WRITE(numout,*) '( x1 - ( assimi(iainia,jv,ilev) + Rd(iainia,ilev) ) ', &
                            ( x1 - (assimi(iainia,jv,ilev) + Rd(iainia,ilev) ))
                       WRITE(numout,*) 'x1,assimi(iainia,jv,ilev),Rd(iainia,ilev): ', &
                            x1,assimi(iainia,jv,ilev),Rd(iainia,ilev)
                       WRITE(numout,*) 'jv,ilev,iainia: ', jv,ilev,iainia
                    ENDIF
                    
                    !++++ TEMP +++++
                    IF(( x1 - (assimi(iainia,jv,ilev) + Rd(iainia,ilev) )) == zero)THEN

                       ! This causes a problem, since we really should not
                       ! be dividing by zero.  However, it seems to happen
                       ! fairly rarely, so all we'll do is cycle to the
                       ! next loop and keep track of this warning so we
                       ! can look at it afterwards.
                       warnings(iainia,jv,iwphoto)=warnings(iainia,jv,iwphoto)+un

                       IF(err_act.GT.1)THEN
                          WRITE(numout,*) 'Diffuco, Getting ready to divide by zero! C3'
                          WRITE(numout,*) '( x1 - ( assimi(iainia,jv,ilev) + Rd(iainia,ilev) ) ', &
                               ( x1 - (assimi(iainia,jv,ilev) + Rd(iainia,ilev) ))
                          WRITE(numout,*) 'x1,assimi(iainia,jv,ilev),Rd(iainia,ilev): ', &
                               x1,assimi(iainia,jv,ilev),Rd(iainia,ilev)
                          WRITE(numout,*) 'jv,ilev,iainia: ', jv,ilev,iainia
                       ENDIF

                       CYCLE

                    ENDIF
                    !++++++++++++++++++++++++++++++

                    IF ( ABS( assimi(iainia,jv,ilev) + Rd(iainia,ilev) ) .LT. min_sechiba ) THEN
                       gs(iainia,ilev) = g0var(iainia)
                    ELSE
                       ! Eq. 18 of Yin et al. (2009)
                       Cc_loc(iainia,ilev) = ( gamma_star(iainia) * x1 + ( assimi(iainia,jv,ilev) + Rd(iainia,ilev) ) * x2 )  &
                            / MAX( min_sechiba, x1 - ( assimi(iainia,jv,ilev) + Rd(iainia,ilev) ) )
                       ! Eq. 17 of Yin et al. (2009)
                       leaf_ci(iainia,jv,ilev) = Cc_loc(iainia,ilev) + assimi(iainia,jv,ilev) / gm(iainia) 
                       ! See eq. right after eq. 15 of Yin et al. (2009)
                       ci_star = gamma_star(iainia) - Rd(iainia,ilev) / gm(iainia)
                       ! 
                       ! Eq. 15 of Yin et al. (2009)
                        IF (ok_hydrol_arch) THEN 
                           gs(iainia,ilev) = MAX(g0var(iainia),( assimi(iainia,jv,ilev) + Rd(iainia,ilev) ) / ( leaf_ci(iainia,jv,ilev) &
                               - ci_star ) * f_gco2(iainia) )
                        ELSE
                           gs(iainia,ilev) = g0var(iainia) + ( assimi(iainia,jv,ilev) + Rd(iainia,ilev) ) / ( leaf_ci(iainia,jv,ilev) &
                               - ci_star ) * f_gco2(iainia) 
                        ENDIF
                    ENDIF


                 ENDIF !C3 vs C4

                 !
                 !
                 ! this is only for output to netcdf
                 ! gs_diffuco_output(iainia,jv,ilev) = gs(iainia,ilev)
                 f_gco2_out(iainia,jv) = f_gco2_out(iainia,jv) + f_gco2(iainia)*lai_per_level(iainia,jv,ilev)
                 
                 !
                 !
                 ! @addtogroup Photosynthesis
                 ! @{   
                 !
                 !! 2.4.4 Estimatation of the stomatal conductance (Ball et al., 1987).\n
                 !! \latexonly
                 !! \input{diffuco_trans_co2_2.4.4.tex}
                 !! \endlatexonly
                 ! @}           
                 !
                 !
                 ! keep stomatal conductance of the levels which we decide
                 ! are the top levels.  Notice that we multiply by the LAI
                 ! in the level here, so this is no longer a per leaf quantity,
                 ! but the whole quantity.
                 IF ( top_level(iainia,ilev) ) THEN
                    gs_top(iainia) = gs_top(iainia)+gs(iainia,ilev)*lai_per_level(iainia,jv,ilev)
                    !
                 ENDIF
                 !
                 ! @addtogroup Photosynthesis
                 ! @{   
                 !
                 !! 2.4.5 Integration at the canopy level\n
                 !! \latexonly
                 !! \input{diffuco_trans_co2_2.4.5.tex}
                 !! \endlatexonly
                 ! @}           
                 ! total assimilation and conductance
                 ! 
                 !++++++++++++++++++++
                 ! The equations give the amount of carbon produced per m**2 if
                 ! leaf surface.  Therefore, since we want GPP in units of m**2 of land surface,
                 ! we have to multiply by the LAI.

                 assimtot(iainia,jv) = assimtot(iainia,jv) + &
                      assimi(iainia,jv,ilev) * lai_per_level(iainia,jv,ilev)
                 Rdtot_pft(iainia,jv) = Rdtot_pft(iainia,jv) + &
                      Rd(iainia,ilev) * lai_per_level(iainia,jv,ilev)

                 IF (ok_mleb .AND. .NOT. hack_can_structure ) THEN 
                     profile_assimtot(iainia,ilev) = profile_assimtot(iainia,ilev) + &
                                                   & assimi(iainia,jv,ilev) * lai_per_level(iainia,jv,ilev)
                     profile_Rdtot(iainia,ilev) = profile_Rdtot(iainia,ilev) + &
                                                   & Rd(iainia,ilev) * lai_per_level(iainia,jv,ilev)
                     profile_gstot(iainia,ilev) = profile_gstot(iainia,ilev) + &
                                                   & gs(iainia,ilev) * lai_per_level(iainia,jv,ilev)
                 END IF ! (hack_can_structure ) THEN 

!                      assimi(iainia) 

                 Rdtot(iainia) = Rdtot(iainia) + &
                      Rd(iainia,ilev) * lai_per_level(iainia,jv,ilev)
                 gstot(iainia) = gstot(iainia) + &
                      gs(iainia,ilev) * lai_per_level(iainia,jv,ilev)
                 assimi_lev(iainia,jv,ilev) = assimi(iainia,jv,ilev)

                 IF(test_pft == jv .AND. printlev_loc>=4)THEN
                    WRITE(numout,'(A,2I5,10ES20.8)') 'assimi middle',ilev,&
                         iainia,assimi(iainia,jv,ilev),assimtot(iainia,jv),gs(iainia,ilev),&
                         gstot(iainia),lai_per_level(iainia,jv,ilev)
                 ENDIF

                 !! (JR 100714) we need to retain the gs (stomatal conductance) per level and per PFT for use
                 !! outside of this routine, in the case that the supply term constraint is applied and the 
                 !! stomatal conductance profile needs to be re-calculated.
                 gs_store(iainia, jv, ilev) = gs(iainia,ilev)
                 !
                 
                 !! (JR 100714) we also retain the gstot (stomatal conductance over all canopy levels, constrained
                 !! by LAI, again for the recalculation of the stomatal conductance profile.
                 gstot_store(iainia, jv) = gstot_store(iainia, jv) + (gs(iainia,ilev) * lai_per_level(iainia,jv,ilev))
                 !
              ENDIF ! if there is absorbed light
              
           ENDDO  ! loop over nlevels_tot

           !! Calculated intercellular CO2 over nlai needed for the chemistry module
           cim(iainia,jv) = zero
           laisum(iainia) = zero
           DO ilev = 1,nlevels_tot
              IF(lai_per_level(iainia,jv,ilev).GT.zero .AND. leaf_ci(iainia,jv,ilev).GT.zero)THEN     
                 cim(iainia,jv)= cim(iainia,jv)+leaf_ci(iainia,jv,ilev)*lai_per_level(iainia,jv,ilev)
                 laisum(iainia)=laisum(iainia) + lai_per_level(iainia,jv,ilev)
              ENDIF
           ENDDO
           IF (laisum(iainia).GT.zero) THEN
              cim(iainia,jv)= cim(iainia,jv)/laisum(iainia)
           ELSE
              cim(iainia,jv)=zero
           ENDIF

        ENDDO assim_loop  ! loop over assimulation points
        
!        IF(jv==test_pft) THEN 
!           CALL histwrite_p(hist_id, 'Cc_loc', kjit, Cc_loc, kjpindex*(nlevels_tot+1), indexlai) 
!           CALL histwrite_p(hist_id, 'Vc', kjit, Vc2, kjpindex*(nlevels_tot+1), indexlai) 
!           CALL histwrite_p(hist_id, 'Vj', kjit, JJ, kjpindex*(nlevels_tot+1), indexlai) 
!           CALL histwrite_p(hist_id, 'limitphoto', kjit, info_limitphoto, kjpindex*nvm*(nlevels_tot+1), indexlai) 
!           CALL histwrite_p(hist_id, 'gammastar', kjit, gamma_star, kjpindex, index) 
!           CALL histwrite_p(hist_id, 'Kmo', kjit, Kmo, kjpindex, index) 
!           CALL histwrite_p(hist_id, 'Kmc', kjit, Kmc, kjpindex, index) 
!           CALL histwrite_p(hist_id, 'gm', kjit, gm, kjpindex, index) 
!           CALL histwrite_p(hist_id, 'gs', kjit, gs, kjpindex*(nlevels_tot+1), indexlai) 
!           CALL histwrite_p(hist_id, 'assimi', kjit, assimi, kjpindex*nvm*(nlevels_tot+1), indexlai) 
!           CALL histwrite_p(hist_id, 'Rd', kjit, Rd, kjpindex*(nlevels_tot+1), indexlai) 
!        ENDIF
        

      !
      !! 2.5 Calculate resistances

        DO inia=1,nia

           !
           iainia=index_assi(inia)
           !

           ! We hit a fringe case above where we could not
           ! calculate the GPP.  Skip this loop so we don't
           ! divide by zero.
           IF(lskip(iainia)) CYCLE


           !! Mean stomatal conductance for CO2 (mol m-2 s-1)
           gsmean(iainia,jv) = gstot(iainia)
           
           !
           ! cimean is the "mean ci" calculated in such a way that assimilation 
           ! calculated in enerbil is equivalent to assimtot
           !

           IF (ok_mleb .AND. .NOT. hack_can_structure ) THEN 
               DO ilev = 1, nlevels_tot
                   profile_gsmean(iainia, jv, ilev) = profile_gstot(iainia, ilev)
                    IF(printlev_loc >= 4) WRITE(numout,*) '110215 profile_gsmean by level: ', &
                         ilev, ' is: ', profile_gsmean(iainia, jv, ilev)
               ENDDO ! ilev = 1, nlevels_tot
           END IF ! (hack_can_structure ) THEN 


           ! cimean is completely diagnostic and not
           ! used for anything.  If the analytical photosynthesis is unable to find
           ! a solution for the A_1 coefficients in any of the levels, assimtot will
           ! be equal to -Rd, which means we will be dividing by zero.  This loop
           ! is to prevent that, and since cimean is diagnostic is doesn't matter
           ! what the value is.
           IF ( ABS(gsmean(iainia,jv)-g0var(iainia)*laisum(iainia)) .GT. min_sechiba) THEN
              cimean(iainia,jv) = (f_gco2(iainia)*(assimtot(iainia,jv)+Rdtot(iainia))) /&
                   (gsmean(iainia,jv)-g0var(iainia)*laisum(iainia)) + gamma_star(iainia) 
           ELSE
              cimean(iainia,jv) = gamma_star(iainia) 
           ENDIF            


           IF (ok_mleb .AND. .NOT. hack_can_structure ) THEN 
              DO ilev = 1, nlevels_tot
                 IF((profile_assimtot(iainia,ilev)+profile_Rdtot(iainia,ilev)) .GT. min_sechiba)THEN
                    profile_cimean(iainia,jv,ilev) = (profile_gs_store(iainia,jv,ilev)-g0(jv)) &
                         / (f_gco2(iainia)*(profile_assimtot(iainia,ilev)+profile_Rdtot(iainia,ilev))) &
                         + gamma_star(iainia)
                    IF(printlev_loc >= 4) THEN
                         WRITE(numout,*) '210515 f_gco2:', f_gco2(iainia)
                         WRITE(numout,*) 'profile_assimtot:', profile_assimtot(iainia,ilev)
                         WRITE(numout,*) 'profile Rdtot:', profile_Rdtot(iainia,ilev)  
                    ENDIF ! IF(printlev >= 4)
                 ELSE
                    ! According to Nicolas, this variable is completely diagnostic and not
                    ! used for anything.  If the analytical photosynthesis is unable to find
                    ! a solution for the A_1 coefficients in any of the levels, assimtot will
                    ! be equal to -Rd, which means we will be dividing by zero.  This loop
                    ! is to prevent that, and since cimean is diagnostic is doesn't matter
                    ! what the value is.
                    profile_cimean(iainia,jv,ilev) = val_exp
                 ENDIF
                 ! WRITE(numout, *) '100215b profile_cimean by level: ', ilev, ' is: ', profile_cimean(iainia, jv, ilev)
              ENDDO ! ilev = 1, nlevels_tot
           END IF ! (hack_can_structure ) THEN 


           ! conversion from umol m-2 (PFT) s-1 to gC m-2 (mesh area) tstep-1
           gpp(iainia,jv) = assimtot(iainia,jv)*12e-6*veget_max(iainia,jv)*dt_sechiba
           istep=istep+1
           
           IF (printlev_loc>=4 .AND. jv==test_pft) THEN
              WRITE(numout,*) 'GPP DIFFUCO: ',istep,gpp(iainia,jv),iainia,jv,&
                   assimtot(iainia,jv),veget_max(iainia,jv),dt_sechiba
           ENDIF

           IF (ok_mleb .AND. .NOT. hack_can_structure ) THEN 
              DO ilev = 1, nlevels_tot
                 profile_gpp(iainia,jv,ilev) = profile_assimtot(iainia, ilev)*12e-6*veget_max(iainia,jv)*dt_sechiba
                 IF(printlev_loc >= 4)  WRITE(numout, *) '100215b profile_gpp by level: ', &
                      ilev, ' is: ', profile_gpp(iainia, jv, ilev)
              END DO ! ilev = 1, nlevels_tot
           END IF ! (hack_can_structure ) THEN 

           ! As in Pearcy, Schulze and Zimmermann
           ! Measurement of transpiration and leaf conductance
           ! Chapter 8 of Plant Physiological Ecology
           ! Field methods and instrumentation, 1991
           ! Editors:
           !
           !    Robert W. Pearcy,
           !    James R. Ehleringer,
           !    Harold A. Mooney,
           !    Philip W. Rundel
           !
           ! ISBN: 978-0-412-40730-7 (Print) 978-94-010-9013-1 (Online)
           !
           !
           !
           !++++++++ CHECK +++++++++
           ! Do we need to normalise gs_top here by the LAI in the level?
           ! In the trunk, it seems like they use an LAI of 0.1, so renormalise
           ! to that.  The default value of lai_top is 0.1.
           !
           !! (JR 100714)
           !! gs_top is normalised in the trunk over the entire canopy LAI, as it
           !! meant to provide backwards compatability with the 'old' method for calculating
           !! the stomatal conductance (in diffuco_trans) that calculates gs according to 
           !! incident SW radiation (rather than CO2 concentration, as here in diffuco_trans_co2)
           !!
           !! rstruct in diffuco_trans is an extra resistance term added to reflect the fact
           !! that less light reaches the lower canopy, and the total resistance used to calculate
           !! vbeta3, and hence the transpiration, is (rstruct + rveget)
           !!
           !! The rstruct term makes the code more difficult to follow, and is no longer has
           !! a physical requirement. We are getting small values of gstop as calculated below,
           !! which was resulting in big resistances when rstruct = min_sechiba and 
           !! rveget = un/gstop and total resistance = (rstruct + rveget) - see commented out
           !! code below.
           !!
           !! The most transparent structure is to set total resistance = rveget = (1 / gstot)
           !! This is analogous to the way in which assimtot, and hence canopy GPP, is calculated
           !! above
           gs_top(iainia)=gs_top(iainia)/totlai_top*lai_top(jv)
           !+++++++++++++++++++++++++

           gstot(iainia) =  mol_to_m_1 *(temp_air(iainia)/tp_00)*&
                (pb_std/pb(iainia))*gstot(iainia)*ratio_H2O_to_CO2
           gstop(iainia) =  mol_to_m_1 * (temp_air(iainia)/tp_00)*&
                (pb_std/pb(iainia))*gs_top(iainia)*ratio_H2O_to_CO2

           IF (ok_mleb .AND. .NOT. hack_can_structure ) THEN 
              DO ilev = 1, nlevels_tot

                  profile_gsmean(iainia,jv,ilev) = profile_gsmean(iainia,jv,ilev) * 1e-6

                  profile_gstot(iainia,ilev) =  mol_to_m_1 *(temp_air(iainia)/tp_00)*&
                    (pb_std/pb(iainia))*profile_gstot(iainia,ilev)*ratio_H2O_to_CO2

                  profile_gstop(iainia,ilev) =  mol_to_m_1 * (temp_air(iainia)/tp_00)*&
                    (pb_std/pb(iainia))*gs_top(iainia)*ratio_H2O_to_CO2*&
                    lai_per_level(iainia,jv,ilev)

              END DO ! ilev = 1, nlevels_tot
           END IF ! (hack_can_structure ) THEN                  

   
           !IF (jv == test_pft .AND. test_grid == 1) THEN 
           !         WRITE(numout, '(A20)'        ) 'Without water stress:'
           !         WRITE(numout, '(A20, I5)'    ) 'pft type is:    ', jc
           !        WRITE(numout, '(A20, ES15.5)') 'Old_Gs_total:', gstot(iainia)            
           !ENDIF   
            
           !++++++ CHECK ++++++
           ! We introduce ::coupled_lai to calculate the :: vbeta3 
           ! In an closed canopy not all the leaves fully interact with the
           ! atmosphere because leaves can shelter each other. In a more open
           ! canopy most leaves can interact with the atmosphere. We use the
           ! ratio of lai over ::Pgap as a proxy for the amount of canopy that
           ! can interact with the atmosphere. Clearly that amount of lai that
           ! interacts should never exceed the:: lai.

           ! Use the lai also as coupled lai
           coupled_lai(iainia, jv) = lai(iainia, jv)

           !---TEMP---
!!$           IF (printlev_loc>=4) THEN
!!$              WRITE(numout,*) 'coupled LAI, ',  jv, coupled_lai(iainia, jv)
!!$              WRITE(numout,*) 'total LAI, ',  jv, lai(iainia, jv)
!!$              WRITE(numout,*) 'veget, ', jv, veget(iainia, jv)
!!$              WRITE(numout,*) 'biomass(ileaf), ',biomass(iainia, jv, ileaf, icarbon) 
!!$           ENDIF
           !----------
      
           !+++++CHECK++++
           !Here, we try to use coupled_lai to reset/refine the total stomata conductance 
           !which interacted with the atmosphere for the transpiration for each pft.
           IF (lai(iainia,jv) > min_sechiba) THEN
               gstot(iainia) = gstot(iainia)* (coupled_lai(iainia,jv)/lai(iainia,jv))
           ELSE
               gstot(iainia) = 0
           END IF
           !++++++++++++++
 
           !+++++++++CHECK++++++++
           ! To get the resonable transpiration we had to set :: lai_top to
           ! ::lai
           ! so we were using the entire canopy. Therefore we no longer use
           ! ::lai_top and simply make use of the entire canopy to calculate the
           ! vegetation resistance from inverse of ::gstot so, we calculate the
           ! :: rveget using :: 1/gstot. This case where gstot(iana)=zero hints
           ! at another problem: we have vegetation but no photosynthesis.
           IF (gstot(iainia).GT.zero)THEN
              rveget(iainia, jv) = un / gstot(iainia)
           ELSE
              rveget(iainia, jv) = undef_sechiba
              CALL ipslerr(2,'diffuco_trans_co2','gstot is zero! This hints at another problem',&
                   'do we have canopy layers with LAI but without absorbed light?','')
           ENDIF
           !+++++++++++++++++++++++++++++++++++++++           

           !rveget(iainia,jv) = un/gstop(iainia)
!!$           rveget_min(iainia,jv) = (defc_plus / kzero(jv)) * &
!!$                (un / SUM(lai_per_level(iainia,jv,:)))
           !
           ! rstruct is the difference between rtot (=1./gstot) and rveget
           !
           ! Correction Nathalie - le 27 Mars 2006 - Interdire a rstruct d'etre negatif
           !rstruct(iainia,jv) = un/gstot(iainia) - &
           !     rveget(iainia,jv)
           
           !++++++++ CHECK +++++++
           !It seems that the calcualtion of rstruct is duplicated with the calculation of rvegt
           !So, we set the ::rstruct as zero 
           !rstruct(iainia,jv) = MAX( un/gstot(iainia) - &
           !     rveget(iainia,jv), min_sechiba)
           !
           !
           !! wind is a global variable of the diffuco module.
           speed = MAX(min_wind, wind(iainia))
           !
           ! beta for transpiration
           !
           ! Corrections Nathalie - 28 March 2006 - on advices of Fred Hourdin
           !! Introduction of a potentiometer rveg_pft to settle the rveg+rstruct sum problem in the coupled mode.
           !! rveg_pft=1 in the offline mode. rveg_pft is a global variable declared in the diffuco module.
           !vbeta3(iainia,jv) = veget_max(iainia,jv) * &
           !  (un - zqsvegrap(iainia)) * &
           !  (un / (un + speed * q_cdrag(iainia) * (rveget(iainia,jv) + &
           !   rstruct(iainia,jv))))
           !! Global resistance of the canopy to evaporation
           !cresist=(un / (un + speed * q_cdrag(iainia) * &
           !     veget(iainia,jv)/veget_max(iainia,jv) * &
           !     (rveg_pft(jv)*(rveget(iainia,jv) + rstruct(iainia,jv)))))

           ! (JR 100714)
           ! previously the resistance term in this expression was (rveget + rstruct), but we now
           ! use rveget alone  
           
           resist_stom(iainia,jv) = veget(iainia,jv)/veget_max(iainia,jv) * (rveg_pft(jv) * (rveget(iainia,jv)))

           cresist(iainia,jv)=(un / (un + speed * q_cdrag(iainia) * resist_stom(iainia,jv)))

           frac_evap_transp(iainia,jv) = veget(iainia,jv)
           
           resist_transp(iainia,jv) = raero(iainia) + resist_stom(iainia,jv)
           

!           cresist=un / (un + speed * q_cdrag(iainia) * &
!                veget(iainia,jv)/veget_max(iainia,jv) * &
!                (rveg_pft(jv)*rveget(iainia,jv)))

           !+++CHECK+++
           !For the interception we assume -as earlier- that the surface area
           !that intercepts is well estimated by ::veget. This is wrong because
           !every leav can contribute to the interception. Depending on weather
           !the parameters for interception are bulk values or single leave. Our
           !approach is correct or needd to be improved.  
           !
           IF ( humrel(iainia,jv) <= min_sechiba ) THEN
              ! Because of a minimum conductance g0, vbeta3 cannot be zero even if humrel=0
              ! in the above equation.
              ! Here, we force transpiration to be zero when the soil cannot deliver it
              frac_evap_transp(iainia,jv) = zero
           END IF

           IF (ok_mleb .AND. .NOT. hack_can_structure ) THEN 
               DO ilev = 1, nlevels_tot
                   
                   ! Fix the gstot equal zero condition to aviod the issue of divided by zero
                   IF ( profile_gstot(iainia,ilev) .LE. min_sechiba ) THEN
                       profile_gstot(iainia,ilev) = 1.0d-6 
                   ENDIF

                   profile_rveget(iainia,jv,ilev) = un/profile_gstot(iainia,ilev)

                   ! need to resolve the rstruct issue
                   profile_rstruct(iainia,jv,ilev) = MAX (un/profile_gstot(iainia,ilev) - &
                                                  & profile_rveget(iainia,jv,ilev), min_sechiba)

                   profile_speed(ilev) = MAX(min_wind, u_speed(iainia,ilev))  ! n.b. need to replace this with another wind profile parameterisation

                   profile_cresist(ilev) = (un / (un + u_speed(iainia,ilev) * q_cdrag(iainia) * & 
                                        & (rveg_pft(jv) * (profile_rveget(iainia,jv,ilev) + profile_rstruct(iainia,jv,ilev)) ) ) )

                   profile_vbeta3(iainia,jv,ilev) = veget_max(iainia,jv) * (un - zqsvegrap(iainia)) * profile_cresist(ilev) + & 
                                        & MIN( vbeta23(iainia,jv), &
                                        veget_max(iainia,jv) * zqsvegrap(iainia) * profile_cresist(ilev) )

               END DO ! ilev = 1, nlevels_tot
           END IF ! (hack_can_structure ) THEN 
           
           IF (printlev_loc>=4) THEN
               DO ilev = 1, nlevels_tot
                   WRITE(numout, *) 'profile_rveget(iainia, jv, ilev) is: ', profile_rveget(iainia, jv, ilev)
               END DO

               DO ilev = 1, nlevels_tot
                   WRITE(numout, *) 'profile_vbeta3(iainia, jv, ilev) is: ', profile_vbeta3(iainia, jv, ilev)
               END DO

               DO ilev = 1, nlevels_tot            
                   WRITE(numout, *) 'profile_cresist(ilev) is: ', (profile_cresist(ilev))
               END DO
           END IF !(printlev_loc) 

           ! vbeta3pot for computation of potential transpiration (needed for irrigation)
           vbeta3pot(iainia,jv) = MAX(zero, veget(iainia,jv) * cresist(iainia,jv))
           !
           !
            
        ENDDO !loop over grid points

        !
     ENDIF ! if nia LT zero
       
  ENDDO  ! loop over vegetation types


  DO ic = 1, kjpindex
     DO jc = 1, nvm 
     !++++++TEMP++++++++++++++++++++++++++++++++++++
     !  IF (printlev_loc>=4 .AND. jc == test_pft) THEN
     !     WRITE(714,   '(A20, ES20.4)') 'Rcanopy,', (un / (un + speed *q_cdrag(ic) * &
     !                        (rveg_pft(jc)*(rveget(ic,jc) + rstruct(ic,jc))))) 
                                                   
     !     WRITE(numout,'(A20, ES20.4)') 'Rcanopy',(un / (un + speed *q_cdrag(ic) * &
     !                        (rveg_pft(jc)*(rveget(ic,jc) + rstruct(ic,jc))))) 
     !  ENDIF
     !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     

     
        IF(veget_max(ic,jc) == zero)THEN
            ! this vegetation type is not present, so no reason to do the 
            ! calculation
            CYCLE
        ENDIF

        DO kc = 1, nlevels_tot
            gs_column_total(ic,jc) = gs_column_total(ic,jc) + gs_store(ic,jc,kc) * lai_per_level(ic,jc,kc) 

            gstot_component(ic,jc,kc) = gs_store(ic,jc,kc) * lai_per_level(ic,jc,kc) 

            IF (gstot_store(ic,jc) .eq. 0.0d0) THEN
                gstot_frac(ic,jc,kc) = 0.0d0
            ELSE
                gstot_frac(ic,jc,kc) = gstot_component(ic,jc,kc) / gstot_store(ic,jc)
            END IF
        END DO ! kc = 1, nlevels_tot 

     END DO ! jc = 1, nvm
  END DO ! ic = 1, kjpindex

    
   Abs_light(:,:,:) = zero
   Abs_light_can(:,:) = zero
   
  DO ic = 1, kjpindex
     DO jc = 1, nvm 
        IF(veget_max(ic,jc) == zero)THEN
            ! this vegetation type is not present, so no reason to do the 
            ! calculation
            CYCLE
        ENDIF

        DO kc = 1, nlevels_tot
           IF (gs_column_total(ic,jc) .ne. zero) THEN
              gs_distribution(ic,jc,kc) = gs_store(ic,jc,kc)* lai_per_level(ic,jc,kc) / gs_column_total(ic,jc)
           ELSE
              gs_distribution(ic,jc,kc) = zero
           END IF
              Abs_light(ic,jc,kc)=Light_Abs_Tot(ic,jc,kc)*light_tran_to_level(ic,jc,kc)
           END DO ! kc = 1, nlevels_tot 
           IF (jv==test_pft .AND. printlev_loc>=4) THEN  
              Abs_light_can(:,:)=SUM(Abs_light(:,:,:),3)
              WRITE(numout,*) 'Absorbed light in the canopy', Abs_light_can(:,jc)
           ENDIF
     END DO ! jc = 1, nvm
  END DO ! ic = 1, kjpindex   

!!$   !++++++ WRITING OUT SOME STUFF
!!$   counter=counter+1
!!$   WRITE(123,'(I8,20F20.10)') counter,gs_distribution(test_grid,test_pft,:)
!!$   WRITE(124,'(I8,20F20.10)') counter,gstot_frac(test_grid,test_pft,:)
!!$   WRITE(125,'(I8,20F20.10)') counter,gstot_store(test_grid,test_pft)
!!$   WRITE(126,'(I8,20F20.10)') counter,gs_store(test_grid,test_pft,:)
!!$   WRITE(127,'(I8,20F20.10)') counter,lai_per_level(test_grid,test_pft,:)
!!$   WRITE(128,'(I8,20F20.10)') counter,gstot_component(test_grid,test_pft,:)
!!$   WRITE(129,'(I8,20F20.10)') counter,gs_column_total(test_grid,test_pft)
!!$   WRITE(130,'(I8,20F20.10)') counter,SUM(lai_per_level(test_grid,test_pft,:))
!!$   IF(rveget(test_grid,test_pft) .LT. 1e19) &
!!$        WRITE(131,'(I8,20ES20.10)') counter,rveget(test_grid,test_pft)
!!$   IF(rveget_min(test_grid,test_pft) .LT. 1e19) &
!!$        WRITE(132,'(I8,20ES20.10)') counter,rveget_min(test_grid,test_pft)
!!$   WRITE(133,'(I8,20ES20.10)') counter,vbeta3(test_grid,test_pft)
!!$   WRITE(134,'(I8,20ES20.10)') counter,vbeta3pot(test_grid,test_pft)
   
  ! Write some output for Machine learning trials
  CALL xios_orchidee_send_field("temp_growth",temp_growth)
  CALL xios_orchidee_send_field("rau",rau)
  CALL xios_orchidee_send_field("Ca",Ca)
  CALL xios_orchidee_send_field("vbeta3pot",vbeta3pot)
  CALL xios_orchidee_send_field("LIGHT_TRAN_TO_LEVEL",Light_Tran_Tot)
  CALL xios_orchidee_send_field("LIGHT_ABS_TO_LEVEL",Light_Abs_Tot)

  IF (printlev_loc>=3) WRITE (numout,*) ' diffuco_trans_co2 done ' 

END SUBROUTINE diffuco_trans_co2


!! ================================================================================================================================
!! SUBROUTINE	   : diffuco_comb
!!
!>\BRIEF           This routine combines the previous partial beta 
!! coefficients and calculates the total alpha and complete beta coefficients.
!!
!! DESCRIPTION	   : Those integrated coefficients are used to calculate (in enerbil.f90) the total evapotranspiration 
!! from the grid-cell. \n
!!
!! In the case that air is more humid than surface, dew deposition can occur (negative latent heat flux). 
!! In this instance, for temperature above zero, all of the beta coefficients are set to 0, except for 
!! interception (vbeta2) and bare soil (vbeta4 with zero soil resistance). The amount of water that is 
!! intercepted by leaves is calculated based on the value of LAI of the surface. In the case of freezing 
!! temperatures, water is added to the snow reservoir, and so vbeta4 and vbeta2 are set to 0, and the 
!! total vbeta is set to 1.\n
!!
!! \latexonly 
!!     \input{diffucocomb1.tex}
!! \endlatexonly
!!
!! The beta and alpha coefficients are initially set to 1.
!! \latexonly 
!!     \input{diffucocomb2.tex}
!! \endlatexonly
!!
!! If snow is lower than the critical value:
!! \latexonly 
!!     \input{diffucocomb3.tex}
!! \endlatexonly
!! If in the presence of dew:
!! \latexonly 
!!     \input{diffucocomb4.tex}
!! \endlatexonly
!!
!! Determine where the water goes (soil, vegetation, or snow)
!! when air moisture exceeds saturation.
!! \latexonly 
!!     \input{diffucocomb5.tex}
!! \endlatexonly
!!
!! If it is not freezing dew is put into the interception reservoir and onto the bare soil. If it is freezing, 
!! water is put into the snow reservoir. 
!! Now modify vbetas where necessary: for soil and snow
!! \latexonly 
!!     \input{diffucocomb6.tex}
!! \endlatexonly
!!
!! and for vegetation
!! \latexonly 
!!     \input{diffucocomb7.tex}
!! \endlatexonly
!!
!! Then compute part of dew that can be intercepted by leafs.
!!
!! There will be no transpiration when air moisture is too high, under any circumstance
!! \latexonly 
!!     \input{diffucocomb8.tex}
!! \endlatexonly
!!
!! There will also be no interception loss on bare soil, under any circumstance.
!! \latexonly 
!!     \input{diffucocomb9.tex}
!! \endlatexonly
!!
!! The flowchart details the 'decision tree' which underlies the module. 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): frac_sub_snow, frac_evap_bare, frac_evap_intercept, frac_evap_transp, resist_equiv_evap, resist_equiv_tot, humrel
!!
!! REFERENCE(S) :
!! - de Noblet-Ducoudré, N, Laval, K & Perrier, A, 1993. SECHIBA, a new set of parameterisations
!! of the hydrologic exchanges at the land-atmosphere interface within the LMD Atmospheric General
!! Circulation Model. Journal of Climate, 6, pp.248-273
!! - Guimberteau, M, 2010. Modélisation de l'hydrologie continentale et influences de l'irrigation
!! sur le cycle de l'eau, PhD Thesis, available from:
!! http://www.sisyphe.upmc.fr/~guimberteau/docs/manuscrit_these.pdf
!!
!! FLOWCHART    :
!! \latexonly 
!!     \includegraphics[scale=0.25]{diffuco_comb_flowchart.png}
!! \endlatexonly
!! \n
!_ ================================================================================================================================

  SUBROUTINE diffuco_comb (kjpindex, humrel, rau, u, v, q_cdrag, pb, qair, temp_sol, temp_air, &
        snow, veget, lai,  &    
        tot_bare_soil, frac_evap_flood, raero, frac_sub_snow, frac_evap_intercept, resist_intercept, &
        frac_evap_transp, resist_transp, frac_evap_bare, &
        evap_bare_lim, evap_bare_lim_ns, veget_max, resist_equiv_evap_veg, resist_equiv_evap, resist_equiv_tot, &
        qsintmax, vbeta2sum, vbeta3sum)

    
    ! Ajout qsintmax dans les arguments de la routine Nathalie / le 13-03-2006

  !! 0. Variable and parameter declaration
    
    !! 0.1 Input variables
    
    INTEGER(i_std), INTENT (in)                          :: kjpindex   !! Domain size (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: rau        !! Air Density (kg m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: u          !! Eastward Lowest level wind speed (m s^{-1}) 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: v          !! Nortward Lowest level wind speed (m s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: q_cdrag    !! Surface drag coefficient  (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: pb         !! Lowest level pressure (hPa)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: qair       !! Lowest level specific air humidity (kg kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: temp_sol   !! Skin temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: temp_air   !! Lower air temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: snow       !! Snow mass (kg)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)    :: veget      !! Fraction of vegetation type (fraction)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)    :: qsintmax   !! Maximum water on vegetation (kg m^{-2})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)    :: lai        !! Leaf area index (m^2 m^{-2})
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)        :: tot_bare_soil       !! Total evaporating bare soil fraction 
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)    :: veget_max           !! Max. fraction of vegetation type (LAI->infty)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: frac_evap_flood     !! Fraction of gricell where floodplains evaporation occurs (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: raero               !! Aerodynamic resistance

    !! 0.2 Output variables
    
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)       :: resist_equiv_evap_veg      !! Equivalent resistance to the fluxes involved in the evaporation flux from the vegetation (interception loss + transpiration + bare soil evaporation) (s m^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)       :: resist_equiv_evap          !! Equivalent resistance to the fluxes involved in the evaporation flux (interception loss + transpiration + bare soil evaporation + floodplains evaporation) (s m^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)       :: resist_equiv_tot           !! Equivalent resistance to the fluxes involved in the total evaporative flux (all evaporations + snow sublimation) (s m^{-1})

    !! 0.3 Modified variables 
    
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)     :: frac_sub_snow       !! Fraction of gridcell where snow sublimation occurs (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)     :: frac_evap_bare      !! Fraction of gridcell where bare soil evaporation occurs (-)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout) :: humrel              !! Soil moisture stress (within range 0 to 1)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout) :: frac_evap_intercept !! Fraction of gridcell where intercepted water evaporation occurs (-)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout) :: frac_evap_transp    !! Fraction of gridcell where transpiration occurs (-)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout) :: resist_intercept    !! Sum of resistances applied to intercepted water evaporation flux (s m^{-1})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout) :: resist_transp       !! Sum of resistances applied to transpiration flux (s m^{-1})
    REAL(r_std), DIMENSION(kjpindex), INTENT (inout)     :: vbeta2sum, vbeta3sum
    
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)     :: evap_bare_lim       !! limiting factor for bare soil evaporation  
                                                                                !! when the 11-layer hydrology is used (-)
    REAL(r_std),DIMENSION (kjpindex,nstm), INTENT (inout):: evap_bare_lim_ns    !! limiting factor for bare soil evaporation 
                                                                                !! when the 11-layer hydrology is used (-)    
    !! 0.4 Local variables
    
    INTEGER(i_std)                                       :: ji, jv
    REAL(r_std)                                          :: zevtest, zsoil_moist, zrapp
    REAL(r_std), DIMENSION(kjpindex)                     :: qsatt
    LOGICAL, DIMENSION(kjpindex)                         :: toveg, tosnow
    REAL(r_std)                                          :: coeff_dew_veg
    REAL(r_std), DIMENSION(kjpindex)                     :: vbeta2sum_temp,vbeta3sum_temp
    REAL(r_std), DIMENSION(kjpindex)                     :: vegtot
!_ ================================================================================================================================
    
    !! \latexonly 
    !!     \input{diffucocomb1.tex}
    !! \endlatexonly

    !! 0. Initialisation of vbeta2sum and vbeta3sum. They  are needed in the 
    !! calculations for the multi-layer  energy budget.
    vbeta2sum(:) = zero
    vbeta3sum(:) = zero
    
    !! 1. If dew is present:
     
    CALL qsatcalc (kjpindex, temp_sol, pb, qsatt)

    
    !! 1.1 Determine where the water goes 
    !! Determine where the water goes (soil, vegetation, or snow)
    !! when air moisture exceeds saturation.
    !! \latexonly 
    !!     \input{diffucocomb5.tex}
    !! \endlatexonly
    toveg(:) = .FALSE.
    tosnow(:) = .FALSE.
    DO ji = 1, kjpindex
      IF ( qsatt(ji) .LT. qair(ji) ) THEN
          IF (temp_air(ji) .GT. tp_00) THEN

              !! If it is not freezing dew is put into the 
              !! interception reservoir and onto the bare soil.
              toveg(ji) = .TRUE.
          ELSE

              !! If it is freezing water is put into the 
              !! snow reservoir.
              tosnow(ji) = .TRUE.
          ENDIF
      ENDIF
    END DO

    !! 1.2 Now modify vbetas where necessary.

    !! 1.2.1 Soil and snow 

    !! \latexonly 
    !!     \input{diffucocomb6.tex}
    !! \endlatexonly

    ! We need to keep consistency between evap_bare_lim, evap_bare_lim_ns and vbeta4 (thus vevapnu)
    ! or we have a water conservation issue in hydrol_split_soil

    DO ji = 1, kjpindex

       IF ( toveg(ji) ) THEN

          frac_sub_snow(ji) = zero

          vegtot(ji) = SUM(veget_max(ji,:))

          IF ( (tot_bare_soil(ji) .GT. min_sechiba) .AND. (vegtot(ji).GT. min_sechiba) ) THEN
             
             frac_evap_bare(ji) = tot_bare_soil(ji)
             
             evap_bare_lim_ns(ji,:) = tot_bare_soil(ji)/vegtot(ji)               

             evap_bare_lim(ji) = frac_evap_bare(ji)
             ! consistent with evap_bare_lim(ji) = SUM(evap_bare_lim_ns(ji,:)*soiltile(ji,:)*vegtot(ji))
             ! as SUM(soiltile(ji,:)) = 1

          ELSE          
             frac_evap_bare(ji) = zero
             evap_bare_lim_ns(ji,:) = zero
             evap_bare_lim(ji) = zero
          ENDIF
       ENDIF

       IF ( tosnow(ji) ) THEN
          frac_sub_snow(ji) = un
          frac_evap_bare(ji) = zero
          evap_bare_lim_ns(ji,:) = zero
          evap_bare_lim(ji) = zero
       ENDIF

    ENDDO

    !! 1.2.2 Vegetation and interception loss
    !! \latexonly 
    !!     \input{diffucocomb7.tex}
    !! \endlatexonly
    DO jv = 1, nvm
      
      DO ji = 1, kjpindex
        IF ( toveg(ji) ) THEN
           IF (qsintmax(ji,jv) .GT. min_sechiba) THEN
              
              ! Compute part of dew that can be intercepted by leafs.
              IF ( lai(ji,jv) .GT. min_sechiba) THEN
                IF (lai(ji,jv) .GT. 1.5) THEN
                   coeff_dew_veg= &
                         &   dew_veg_poly_coeff(6)*lai(ji,jv)**5 &
                         & - dew_veg_poly_coeff(5)*lai(ji,jv)**4 &
                         & + dew_veg_poly_coeff(4)*lai(ji,jv)**3 &
                         & - dew_veg_poly_coeff(3)*lai(ji,jv)**2 &
                         & + dew_veg_poly_coeff(2)*lai(ji,jv) &
                         & + dew_veg_poly_coeff(1)
                 ELSE
                    coeff_dew_veg=un
                 ENDIF
              ELSE
                 coeff_dew_veg=zero
              ENDIF
              IF (jv .EQ. 1) THEN
                ! This line may not work with CWRR when frac_bare is distributed
                ! among three soiltiles. Fortunately, qsintmax(ji,1)=0 (LAI=0 in
                ! PFT1) so we never pass here
                 frac_evap_intercept(ji,jv) = coeff_dew_veg*tot_bare_soil(ji)
              ELSE
                 frac_evap_intercept(ji,jv) = coeff_dew_veg*veget(ji,jv)
              ENDIF
           ELSE
             frac_evap_intercept(ji,jv) = zero ! if qsintmax=0, vbeta2=0
           ENDIF
        ENDIF
        IF ( tosnow(ji) ) frac_evap_intercept(ji,jv) = zero
        
      ENDDO
      
    ENDDO

    !! 1.2.3 Vegetation and transpiration  
    !! There will be no transpiration when air moisture is too high, under any circumstance
    !! \latexonly 
    !!     \input{diffucocomb8.tex}
    !! \endlatexonly
    DO jv = 1, nvm
      DO ji = 1, kjpindex
        IF ( qsatt(ji) .LT. qair(ji) ) THEN
          frac_evap_transp(ji,jv) = zero
          humrel(ji,jv) = zero
       ENDIF
      ENDDO
    ENDDO
    
    
    !! 1.2.4 Overrules 1.2.2
    !! There will also be no interception loss on bare soil, under any circumstance.
    !! \latexonly 
    !!     \input{diffucocomb9.tex}
    !! \endlatexonly
    DO ji = 1, kjpindex
       IF ( qsatt(ji) .LT. qair(ji) ) THEN
          frac_evap_intercept(ji,1) = zero
       ENDIF
    ENDDO

    !! 2. Now calculate vbeta in all cases (the equality needs to hold for enerbil to be consistent)
    DO ji = 1, kjpindex



          vbeta2sum(ji) = SUM(frac_evap_intercept(ji,:)/resist_intercept(ji,:))
          vbeta3sum(ji) = SUM(frac_evap_transp(ji,:)/resist_transp(ji,:))
          IF ((frac_evap_bare(ji)/raero(ji) + SUM(frac_evap_intercept(ji,:)/resist_intercept(ji,:) &
               + frac_evap_transp(ji,:)/resist_transp(ji,:))) .LT. min_sechiba) THEN

             frac_evap_bare(ji) = zero
             frac_evap_intercept(ji,:)= zero
             frac_evap_transp(ji,:)= zero

             evap_bare_lim_ns(ji,:) = zero
             evap_bare_lim(ji) = zero
          END IF

          !! Calculation of the equivalent resistances :
          !! In the current calculations, all evaporative fluxes are considered in parallel between the surface and atmosphere.
          !! When trying to calculate the total fluxes in enerbil, one needs to know the equivalent resistance related to the parallel fluxes involved.
          !! To calculate them, we sum the invert of the resistances (conductances) multiplied by their fraction of gridcell to get the invert of the 
          !! equivalent resistance (equivalent conductance). See the technical note about this calculation. In order to avoid division by zero, we 
          !! limitate the value to 10**(-20).

          resist_equiv_evap(ji) = un / MAX(10.0**(-20),( (un - frac_sub_snow(ji)) * (un - frac_evap_flood(ji)) * &
               ( frac_evap_bare(ji) /raero(ji) + SUM(frac_evap_intercept(ji,:)/resist_intercept(ji,:)) + &
               SUM(frac_evap_transp(ji,:)/resist_transp(ji,:)) ) + frac_evap_flood(ji)/raero(ji) ))

          resist_equiv_evap_veg(ji) = un / MAX(10.0**(-20),( (un - frac_sub_snow(ji)) * (un - frac_evap_flood(ji)) * &
               ( frac_evap_bare(ji) /raero(ji) + SUM(frac_evap_intercept(ji,:)/resist_intercept(ji,:)) + &
               SUM(frac_evap_transp(ji,:)/resist_transp(ji,:)) ) ))

          resist_equiv_tot(ji) = un /  MAX(10.0**(-20),( frac_sub_snow(ji) * (un - frac_evap_flood(ji)) / raero(ji) + &
                (un - frac_sub_snow(ji)) * (un - frac_evap_flood(ji)) * ( frac_evap_bare(ji) /raero(ji) + &
                SUM(frac_evap_intercept(ji,:)/resist_intercept(ji,:)) + SUM(frac_evap_transp(ji,:)/resist_transp(ji,:)) ) + &
                frac_evap_flood(ji)/raero(ji)))
    ENDDO

    CALL xios_orchidee_send_field("evap_bare_lim",evap_bare_lim) 
    CALL xios_orchidee_send_field("evap_bare_lim_ns",evap_bare_lim_ns)

    IF (printlev>=3) WRITE (numout,*) ' diffuco_comb done '

  END SUBROUTINE diffuco_comb


!! ================================================================================================================================
!! SUBROUTINE	: diffuco_raerod
!!
!>\BRIEF	Computes the aerodynamic resistance, for cases in which the
!! surface drag coefficient is provided by the coupled atmospheric model LMDZ and  when the flag
!! 'ldq_cdrag_from_gcm' is set to TRUE
!!
!! DESCRIPTION	: Simply computes the aerodynamic resistance, for cases in which the
!! surface drag coefficient is provided by the coupled atmospheric model LMDZ. If the surface drag coefficient
!! is not provided by the LMDZ or signalled by the flag 'ldq_cdrag_from_gcm' set to FALSE, then the subroutine
!! diffuco_aero is called instead of this one.
!!
!! Calculation of the aerodynamic resistance, for diganostic purposes. First calculate wind speed:
!! \latexonly 
!!     \input{diffucoaerod1.tex}
!! \endlatexonly       
!!
!! next calculate ::raero
!! \latexonly 
!!     \input{diffucoaerod2.tex}
!! \endlatexonly
!! 
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::raero
!!
!! REFERENCE(S)	:
!! - de Noblet-Ducoudré, N, Laval, K & Perrier, A, 1993. SECHIBA, a new set of parameterisations
!! of the hydrologic exchanges at the land-atmosphere interface within the LMD Atmospheric General
!! Circulation Model. Journal of Climate, 6, pp.248-273
!! - Guimberteau, M, 2010. Modélisation de l'hydrologie continentale et influence de l'irrigation
!! sur le cycle de l'eau, PhD Thesis, available from:
!! http://www.sisyphe.upmc.fr/~guimberteau/docs/manuscrit_these.pdf
!!
!! FLOWCHART    :  None
!! \n
!_ ================================================================================================================================

  SUBROUTINE diffuco_raerod (kjpindex, u, v, q_cdrag, raero)
    
    IMPLICIT NONE
    
  !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                     :: kjpindex     !! Domain size (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)  :: u            !! Eastward Lowest level wind velocity (m s^{-1}) 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)  :: v            !! Northward Lowest level wind velocity (m s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)  :: q_cdrag      !! Surface drag coefficient  (-)
    
    !! 0.2 Output variables 
    
    REAL(r_std),DIMENSION (kjpindex), INTENT (out) :: raero        !! Aerodynamic resistance (s m^{-1})
     
    !! 0.3 Modified variables

    !! 0.4 Local variables
    
    INTEGER(i_std)                                 :: ji           !! (-)
    REAL(r_std)                                    :: speed        !! (m s^{-1})
!_ ================================================================================================================================
   
  !! 1. Simple calculation of the aerodynamic resistance, for diganostic purposes.

    DO ji=1,kjpindex

       !! \latexonly 
       !!     \input{diffucoaerod1.tex}
       !! \endlatexonly       
       speed = MAX(min_wind, wind(ji))

       !! \latexonly 
       !!     \input{diffucoaerod2.tex}
       !! \endlatexonly
       raero(ji) = un / (q_cdrag(ji)*speed)
       
    ENDDO
  
  END SUBROUTINE diffuco_raerod


!! ================================================================================================================================
!! SUBROUTINE   : isotope_c13
!!
!>\BRIEF         Calculate the fractionation of C13 during photosynthesis\n
!!
!! DESCRIPTION  :   
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): :: delta_c13_assim
!!
!! REFERENCE(S) : 
!!	- Farquhar, G. D., Ehleringer, J. R., & Hubick, K. T. (1989). CARBON ISOTOPE DISCRIMINATION AND PHOTOSYNTHESIS. 
!!	  Annual Review of Plant Physiology and Plant Molecular Biology, 40, 503-537. doi: 10.1146/annurev.arplant.40.1.50
!!
!! FLOWCHART    : None
!_ ================================================================================================================================

  SUBROUTINE isotope_c13(kjpindex, index_assi, leaf_ci, ccanopy, lai_per_level, &
       assimtot, delta_c13_assim, leaf_ci_out, nvm, &
       nlevels_tot, assimi_lev)
    
    IMPLICIT NONE

    !! 0. Variables and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT (in)                            :: kjpindex          !! number of land pixels
    INTEGER(i_std), INTENT (in)                            :: nvm               !! number of pfts
    INTEGER(i_std), INTENT (in)                            :: nlevels_tot       !! number of lai layers in canopy
    INTEGER(i_std), DIMENSION(:), INTENT (in)              :: index_assi        !! indices (unitless)   
    REAL(r_std), DIMENSION (:,:,:), INTENT (in)            :: leaf_ci           !! intercellular CO2 concentration (ppm)
    REAL(r_std),DIMENSION (:), INTENT (in)                 :: ccanopy           !! CO2 concentration inside the canopy (ppm)   
    REAL(r_std), DIMENSION(:,:,:), INTENT (in)             :: lai_per_level     !! This is the LAI per vertical level
                                                                                !! @tex $(m^{2} m^{-2})$ @endtex 
    !+++CHECK+++
    ! Looks like assimi is a local variable in trans_co2 and that
    ! its values are overwritten for every new layer. If this is 
    ! indeed the case, it should not be used outside trans_co2    
!    REAL(r_std), DIMENSION(:,:), INTENT (in)               :: assimi            !! assimilation (at a specific LAI level)
!                                                                                !! @tex ($\mu mol m^{-2} s^{-1}$) @endtex
    !+++++++++++    
    REAL(r_std), DIMENSION(:,:), INTENT (in)               :: assimtot          !! total assimilation
                                                                                !! @tex ($\mu mol m^{-2} s^{-1}$) @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT (in)             :: assimi_lev        !! assimilation (at all LAI levels)
                                                                                !! @tex ($\mu mol m^{-2} s^{-1}$) @endtex

    !! 0.2 Output variables
    REAL(r_std),DIMENSION (:,:), INTENT (out)              :: delta_c13_assim   !! C13 concentration in delta notation 
                                                                                !! @tex $ permille $ @endtex (per thousand)   
    REAL(r_std),DIMENSION (:,:), INTENT (out)              :: leaf_ci_out       !! Ci @tex $ ??unit?? $ @endtex      
    
    !! 0.3 Modified variables

    !! 0.4 Local variables 
    INTEGER(i_std)                                          :: inia,iainia      !! counter/indices (unitless)
    INTEGER(i_std)                                          :: ji,jv,jl         !! Indices (-)    
    

!_ ================================================================================================================================    
    
    !! Initialize
    delta_c13_assim(:,:) = zero
    leaf_ci_out(:,:) = zero

    !! Calculate delta_C13
    DO ji=1,kjpindex

       ! Skip the whole pixel of none of the PFTs did photosynthesis
       IF (SUM(SUM(assimi_lev(ji,:,:),2),1) .LE. zero) THEN
            ! There was no photsynthesis at this pixel, no reason to do the 
            ! calculation
          CYCLE
       ENDIF

       DO jv = 1,nvm

          ! Skip the PFT is it didn't do photosynthesis
          IF (SUM(assimi_lev(ji,jv,:)) .LE. zero) THEN
             ! There was no photsynthesis at this pixel, no reason to do the 
             ! calculation
             CYCLE
          ENDIF

          ! There was photosynthesis at this PFT, calculate the C13 fractionation
          DO jl = 1,nlevels_tot
  
             IF (lai_per_level(ji,jv,jl) .GT. zero) THEN

                IF (printlev_loc>=4) THEN
                   WRITE(numout,*) 'leaf_ci, ',leaf_ci(ji,jv,jl)
                   WRITE(numout,*) 'assimi_lev, ',assimi_lev(ji,jv,jl)
                   WRITE(numout,*) 'assimtot - isotope - before calc, ',assimtot(ji,jv)
                   WRITE(numout,*) 'lai_per_level, ',lai_per_level(ji,jv,jl)
                   WRITE(numout,*) 'day, ',day_end
                   WRITE(numout,*) 'pft, ',jv
                ENDIF

                ! Simple version of Farquhar model
                ! delta_c13_assim(:,jv) = (-8 - c13_a - (c13_b - c13_a) * cimean(:,jv) / ccanopy(:)) 
                delta_c13_assim(ji,jv) = delta_c13_assim(ji,jv) + (- c13_a - (c13_b - c13_a)* &
                     leaf_ci(ji,jv,jl)/ccanopy(ji)) * &
                     assimi_lev(ji,jv,jl) * lai_per_level(ji,jv,jl)
                leaf_ci_out(ji,jv) = leaf_ci_out(ji,jv) +  &
                     leaf_ci(ji,jv,jl) * &
                     assimi_lev(ji,jv,jl) * lai_per_level(ji,jv,jl)
                
             ENDIF
             
             ! Simple version of Farquhar model with photorespiratioin, gi and gamma*
             ! This code copied from the work by Thomas Eglin and Thomas Launois where it
             ! was already commented out. When moving it to ORCHIDEE-CAN is has not been 
             ! tested with ORCHIDEE-CAN. It is not clear why it was commented out.
             ! was commented out
!!$             IF (assimi(ji,jv) .GT. 0) THEN 
!!$                delta_c13_assim(ji,jv) = delta_c13_assim(ji,jv) + (&
!!$                     - c13_ab * (ccanopy(ji)-leaf_cs(ji,jv,jl))/(ccanopy(ji)) &
!!$                     - c13_a * (leaf_cs(ji,jv,jl)-leaf_ci(ji,jv,jl))/(ccanopy(ji)) &
!!$                     - (C13_es + c13_al) * (leaf_ci(ji,jv,jl)-leaf_cc(ji,jv,jl))/(ccanopy(ji)) &
!!$                     - (c13_b+2.0) * leaf_cc(ji,jv,jl)/(ccanopy(ji)) &
!!$                     - (C13_e * (Rd(ji)/(Rd(ji) + assimi(ji,jv))) + C13_f * CP(ji))/(ccanopy(ji))&
!!$                     ) * assimi(ji,jv) * (laitab(jl+1)-laitab(jl))
!!$             ENDIF

             
!!$                  (Rd(ji) + assimi(ji,jv)), Rd(ji),assimi(ji,jv)
             
             ! This code copied from the work by Thomas Eglin and Thomas Launois where it
             ! was already commented out. When moving it to ORCHIDEE-CAN is has not been 
             ! tested with ORCHIDEE-CAN. It is not clear why it was commented out.
             ! was commented out
             ! Modification M; Bordigoni 07/2008 : New calculation of c3-plants fractionation
             ! according to Lloyd and Farquhar_Oecologia_1994
!!$             delta_c13_assim(ji,jv) =  delta_c13_assim(ji,jv) +(- (&
!!$                  4.4*(1-(leaf_ci(ji,jv,jl) / ccanopy(ji) )+0.025)+ &
!!$                  0.075*(1.1+0.7)+ &
!!$                  27.5*(( leaf_ci(ji,jv,jl) / ccanopy(ji))-0.1)- &
!!$                  (11*0.1+8*(1.54*1.05*((temp_air(ji)-273.16)+2.5)))/ ccanopy(ji)) )* &
!!$                  assimi(ji,jv) * (laitab(jl+1)-laitab(jl))

          ENDDO ! photosynthesis layers

          IF (assimtot(ji,jv) .GT. threshold_c13_assim) THEN

             delta_c13_assim(ji,jv) = delta_c13_assim(ji,jv)/assimtot(ji,jv)
             leaf_ci_out(ji,jv) = leaf_ci_out(ji,jv)/assimtot(ji,jv)

          ELSE

             delta_c13_assim(ji,jv) = zero
             leaf_ci_out(ji,jv) = zero

          ENDIF

       ENDDO ! # PFTs
          
    ENDDO ! #kjpindex

  END SUBROUTINE isotope_c13


END MODULE diffuco
