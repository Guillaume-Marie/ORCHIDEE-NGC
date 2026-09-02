! ===================================================================================================\n
! MODULE        : hydrol
!
! CONTACT       : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE       : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        This module computes the soil moisture processes on continental points. 
!!
!!\n DESCRIPTION : contains hydrol_main, hydrol_initialize, hydrol_finalise, hydrol_init, 
!!                 hydrol_var_init, hydrol_waterbal, hydrol_alma,
!!                 hydrol_vegupd, hydrol_canop, hydrol_flood, hydrol_soil, hydrol_root_profile.
!!                 The assumption in this module is that very high vertical resolution is
!!                 needed in order to properly resolve the vertical diffusion of water in
!!                 the soils. Furthermore we have taken into account the sub-grid variability
!!                 of soil properties and vegetation cover by allowing the co-existence of
!!                 different soil moisture columns in the same grid box.
!!                 This routine was originaly developed by Patricia deRosnay.
!!
!! RECENT CHANGE(S) : November 2020: It is possible to define soil hydraulic parameters from maps,
!!                    as needed for the SP-MIP project (Tafasca Salma and Ducharne Agnes).
!!                    Here, it leads to change dimensions and indices. 
!!                    We can also impose kfact_root=1 in all soil layers to cancel the effect of
!!                    roots on ks profile (keyword KFACT_ROOT_CONST).
!!                    Octobre 2023: New irrigation scheme. Here new irrigation demand based in
!!                    soil moisture deficit, and irrigation application.
!!                 
!! REFERENCE(S) :
!! - de Rosnay, P., J. Polcher, M. Bruen, and K. Laval, Impact of a physically based soil
!! water flow and soil-plant interaction representation for modeling large-scale land surface
!! processes, J. Geophys. Res, 107 (10.1029), 2002. \n
!! - de Rosnay, P. and Polcher J. (1998) Modeling root water uptake in a complex land surface scheme coupled 
!! to a GCM. Hydrology and Earth System Sciences, 2(2-3):239-256. \n
!! - de Rosnay, P., M. Bruen, and J. Polcher, Sensitivity of surface fluxes to the number of layers in the soil 
!! model used in GCMs, Geophysical research letters, 27 (20), 3329 - 3332, 2000. \n
!! - d’Orgeval, T., J. Polcher, and P. De Rosnay, Sensitivity of the West African hydrological
!! cycle in ORCHIDEE to infiltration processes, Hydrol. Earth Syst. Sci. Discuss, 5, 2251 - 2292, 2008. \n
!! - Carsel, R., and R. Parrish, Developing joint probability distributions of soil water retention
!! characteristics, Water Resources Research, 24 (5), 755 - 769, 1988. \n
!! - Mualem, Y., A new model for predicting the hydraulic conductivity of unsaturated porous
!! media, Water Resources Research, 12 (3), 513 - 522, 1976. \n
!! - Van Genuchten, M., A closed-form equation for predicting the hydraulic conductivity of
!! unsaturated soils, Soil Science Society of America Journal, 44 (5), 892 - 898, 1980. \n
!! - Campoy, A., Ducharne, A., Cheruy, F., Hourdin, F., Polcher, J., and Dupont, J.-C., Response 
!! of land surface fluxes and precipitation to different soil bottom hydrological conditions in a
!! general circulation model,  J. Geophys. Res, in press, 2013. \n
!! - Gouttevin, I., Krinner, G., Ciais, P., Polcher, J., and Legout, C. , 2012. Multi-scale validation
!! of a new soil freezing scheme for a land-surface model with physically-based hydrology.
!! The Cryosphere, 6, 407-430, doi: 10.5194/tc-6-407-2012. \n
!! - Tafasca S. (2020). Evaluation de l’impact des propriétés du sol sur l’hydrologie simulee dans le
!! modèle ORCHIDEE, PhD thesis, Sorbonne Universite. \n
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_sechiba/hydrol.f90 $
!! $Date: 2026-05-12 21:40:35 +0200 (mar. 12 mai 2026) $
!! $Revision: 9536 $
!! \n
!_ ===============================================================================================\n
MODULE hydrol

  USE ioipsl
  USE xios_orchidee
  USE constantes
  USE time, ONLY : one_day, dt_sechiba, julian_diff, one_year, FirstTsYear
  USE constantes_soil
  USE pft_parameters
  USE dynamic_parameters
  USE sechiba_io_p
  USE grid
  USE explicitsnow
  USE function_library, ONLY : get_printlev, Arrhenius, wood_to_dia, wood_to_cv, wood_to_height

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: hydrol_main, hydrol_initialize, hydrol_finalize, hydrol_clear

  !
  ! variables used inside hydrol module : declaration and initialisation
  !
  LOGICAL, SAVE                                   :: doponds=.FALSE.           !! Reinfiltration flag (true/false)
!$OMP THREADPRIVATE(doponds)
  REAL(r_std), SAVE                               :: froz_frac_corr            !! Coefficient for water frozen fraction correction
!$OMP THREADPRIVATE(froz_frac_corr)
  REAL(r_std), SAVE                               :: max_froz_hydro            !! Coefficient for water frozen fraction correction
!$OMP THREADPRIVATE(max_froz_hydro)
  REAL(r_std), SAVE                               :: smtot_corr                !! Coefficient for water frozen fraction correction
!$OMP THREADPRIVATE(smtot_corr)
  LOGICAL, SAVE                                   :: do_rsoil                  !! Flag to calculate rsoil for bare soile evap (true/false)
!$OMP THREADPRIVATE(do_rsoil)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: rsoil_scale_tile          !! Factor to scale rsoil for bare soil evap per soil tile
!$OMP THREADPRIVATE(rsoil_scale_tile)
  LOGICAL, SAVE                                   :: kfact_root_const          !! Control kfact_root calculation, set constant kfact_root=1 if kfact_root_const=true
!$OMP THREADPRIVATE(kfact_root_const)
  CHARACTER(LEN=80) , SAVE                        :: var_name                  !! To store variables names for I/O
!$OMP THREADPRIVATE(var_name)

  !
  REAL(r_std), PARAMETER                          :: allowed_err =  2.0E-8_r_std
  REAL(r_std), PARAMETER                          :: EPS1 = EPSILON(un)      !! A small number
  ! one dimension array allocated, computed, saved and got in hydrol module
  ! Values per soil type

  REAL(r_std), SAVE                               :: n0                    !! Variable for parameter CWRR_NKS_N0 (unitless)
!$OMP THREADPRIVATE(n0)
  REAL(r_std), SAVE                               :: nk_rel                !! Variable for parameter CWRR_NKS_POWER (unitless)
!$OMP THREADPRIVATE(nk_rel)

  INTEGER(i_std), SAVE :: printlev_loc                                    !! Local level of text output for current module
!$OMP THREADPRIVATE(printlev_loc)
  !! Campbell Parametrisation sometimes used in the muff

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: b_muff               !! Campbell coeficients b (unitless)
!$OMP THREADPRIVATE(b_muff)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: psi_air_entry        !! Air entry water potential (MPa)
!$OMP THREADPRIVATE(psi_air_entry)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: mcr_sup              !! Residual soil water content in the superior soil layer (m3/m3) (tuzet)
!$OMP THREADPRIVATE(mcr_sup)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: mcr_inf              !! Residual soil water content in the inferior soil layer (m3/m3) (tuzet)
!$OMP THREADPRIVATE(mcr_inf)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: mcs_sup              !! Saturated soil water content in the superior soil layer (m3/m3) (tuzet)
!$OMP THREADPRIVATE(mcs_sup)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: mcs_inf              !! Saturated soil water content in the inferior soil layer (m3/m3) (tuzet)
!$OMP THREADPRIVATE(mcs_inf)
  LOGICAL, SAVE                                   :: is_vg                !! Flag to control the hydraulic conductivity and water potential
                                                                          !! calculations in the new hydrualic architecture
!$OMP THREADPRIVATE(is_vg)


  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: pcent               !! Fraction of saturated volumetric soil moisture above 
                                                                         !! which transpir is max (0-1, unitless)
!$OMP THREADPRIVATE(pcent)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: mc_awet             !! Vol. wat. cont. above which albedo is cst
                                                                         !!  @tex $(m^{3} m^{-3})$ @endtex 
!$OMP THREADPRIVATE(mc_awet)                                             
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: mc_adry             !! Vol. wat. cont. below which albedo is cst
                                                                         !!  @tex $(m^{3} m^{-3})$ @endtex 
!$OMP THREADPRIVATE(mc_adry)                                             
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: tot_watveg_beg   !! Total amount of water on vegetation at start of time 
                                                                         !! step @tex $(kg m^{-2})$ @endtex 
!$OMP THREADPRIVATE(tot_watveg_beg)                                      
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: tot_watveg_end   !! Total amount of water on vegetation at end of time step
                                                                         !!  @tex $(kg m^{-2})$ @endtex 
!$OMP THREADPRIVATE(tot_watveg_end)                                      
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: tot_watsoil_beg  !! Total amount of water in the soil at start of time step
                                                                         !!  @tex $(kg m^{-2})$ @endtex 
!$OMP THREADPRIVATE(tot_watsoil_beg)                                     
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: tot_watsoil_end  !! Total amount of water in the soil at end of time step
                                                                         !!  @tex $(kg m^{-2})$ @endtex 
!$OMP THREADPRIVATE(tot_watsoil_end)                                     
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: snow_beg         !! Total amount of snow at start of time step
                                                                         !!  @tex $(kg m^{-2})$ @endtex 
!$OMP THREADPRIVATE(snow_beg)                                            
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: snow_end         !! Total amount of snow at end of time step
                                                                         !!  @tex $(kg m^{-2})$ @endtex 
!$OMP THREADPRIVATE(snow_end)                                           
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: delsoilmoist     !! Change in soil moisture @tex $(kg m^{-2})$ @endtex 
!$OMP THREADPRIVATE(delsoilmoist)                                         
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: delintercept     !! Change in interception storage
                                                                         !!  @tex $(kg m^{-2})$ @endtex 
!$OMP THREADPRIVATE(delintercept)                                        
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: delswe           !! Change in SWE @tex $(kg m^{-2})$ @endtex 
!$OMP THREADPRIVATE(delswe)
  REAL(r_std),ALLOCATABLE, SAVE, DIMENSION (:)       :: undermcr         !! Nb of tiles under mcr for a given time step
!$OMP THREADPRIVATE(undermcr) 
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: mask_veget       !! zero/one when veget fraction is zero/higher (1)
!$OMP THREADPRIVATE(mask_veget)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: mask_soiltile    !! zero/one where soil tile is zero/higher (1)
!$OMP THREADPRIVATE(mask_soiltile)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: humrelv          !! Water stress index for transpiration
                                                                         !! for each soiltile x PFT couple (0-1, unitless)
!$OMP THREADPRIVATE(humrelv)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: vegstressv       !! Water stress index for vegetation growth
                                                                         !! for each soiltile x PFT couple (0-1, unitless)
!$OMP THREADPRIVATE(vegstressv)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: precisol         !! Throughfall+Totmelt per PFT 
                                                                         !!  @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(precisol)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: throughfall      !! Throughfall per PFT 
                                                                         !!  @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(throughfall)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: precisol_ns      !! Throughfall per soiltile 
                                                                         !!  @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(precisol_ns)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: ae_ns            !! Bare soil evaporation per soiltile
                                                                         !!  @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(ae_ns)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: free_drain_coef  !! Coefficient for free drainage at bottom
                                                                         !!  (0-1, unitless) 
!$OMP THREADPRIVATE(free_drain_coef) 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: zwt_force        !! Prescribed water table depth (m)
!$OMP THREADPRIVATE(zwt_force) 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: frac_bare_ns     !! Evaporating bare soil fraction per soiltile
                                                                         !!  (0-1, unitless)
!$OMP THREADPRIVATE(frac_bare_ns)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: rootsink         !! Transpiration sink by soil layer and soiltile
                                                                         !! @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(rootsink)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: subsnowveg       !! Sublimation of snow on vegetation 
                                                                         !!  @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(subsnowveg)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: subsnownobio     !! Sublimation of snow on other surface types  
                                                                         !! (ice, lakes,...) @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(subsnownobio)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: icemelt          !! Ice melt @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(icemelt)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: subsinksoil      !! Excess of sublimation as a sink for the soil
                                                                         !! @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(subsinksoil)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: vegtot           !! Total Total fraction of grid-cell covered by PFTs
                                                                         !! (bare soil + vegetation) (1; 1)
!$OMP THREADPRIVATE(vegtot)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: resdist          !! Soiltile values from previous time-step (1; 1)
!$OMP THREADPRIVATE(resdist)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: vegtot_old       !! Total Total fraction of grid-cell covered by PFTs
                                                                         !! from previous time-step (1; 1)
!$OMP THREADPRIVATE(vegtot_old)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: mx_eau_var       !! Maximum water content of the soil @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(mx_eau_var)

  ! arrays used by cwrr scheme
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: kfact_root       !! Factor to increase Ks towards the surface
                                                                         !! (unitless)
                                                                         !! DIM = kjpindex * nslm * nstm
!$OMP THREADPRIVATE(kfact_root)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: kfact            !! Factor to reduce Ks with depth (unitless)
                                                                         !! DIM = nslm * kjpindex
!$OMP THREADPRIVATE(kfact)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: zz               !! Depth of nodes [znh in vertical_soil] transformed into (mm) 
!$OMP THREADPRIVATE(zz)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: dz               !! Internode thickness [dnh in vertical_soil] transformed into (mm)
!$OMP THREADPRIVATE(dz)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: dh               !! Layer thickness [dlh in vertical_soil] transformed into (mm)
!$OMP THREADPRIVATE(dh)
  INTEGER(i_std), SAVE                               :: itopmax          !! Number of layers where the node is above 0.1m depth
!$OMP THREADPRIVATE(itopmax)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: mc_lin   !! 50 Vol. Wat. Contents to linearize K and D, for each texture 
                                                                 !!  @tex $(m^{3} m^{-3})$ @endtex
                                                                 !! DIM = imin:imax * kjpindex
!$OMP THREADPRIVATE(mc_lin)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: k_lin    !! 50 values of unsaturated K, for each soil layer and texture
                                                                 !!  @tex $(mm d^{-1})$ @endtex
                                                                 !! DIM = imin:imax * nslm * kjpindex
!$OMP THREADPRIVATE(k_lin)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: d_lin    !! 50 values of diffusivity D, for each soil layer and texture
                                                                 !!  @tex $(mm^2 d^{-1})$ @endtex
                                                                 !! DIM = imin:imax * nslm * kjpindex
!$OMP THREADPRIVATE(d_lin)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: a_lin    !! 50 values of the slope in K=a*mc+b, for each soil layer and texture
                                                                 !!  @tex $(mm d^{-1})$ @endtex
                                                                 !! DIM = imin:imax * nslm * kjpindex
!$OMP THREADPRIVATE(a_lin)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: b_lin    !! 50 values of y-intercept in K=a*mc+b, for each soil layer and texture
                                                                 !!  @tex $(m^{3} m^{-3})$ @endtex
                                                                 !! DIM = imin:imax * nslm * kjpindex
!$OMP THREADPRIVATE(b_lin)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: humtot   !! Total Soil Moisture @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(humtot)
  LOGICAL, ALLOCATABLE, SAVE, DIMENSION (:)          :: resolv   !! Mask of land points where to solve the diffusion equation
                                                                 !! (true/false)
!$OMP THREADPRIVATE(resolv)

!! for output
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: kk_moy   !! Mean hydraulic conductivity over soiltiles (mm/d)
!$OMP THREADPRIVATE(kk_moy)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: kk       !! Hydraulic conductivity for each soiltiles (mm/d)
!$OMP THREADPRIVATE(kk)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: avan_mod_tab  !! VG parameter a modified from  exponantial profile
                                                                      !! @tex $(mm^{-1})$ @endtex !! DIMENSION (nslm,kjpindex)
!$OMP THREADPRIVATE(avan_mod_tab)  
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: nvan_mod_tab  !! VG parameter n  modified from  exponantial profile
                                                                      !! (unitless) !! DIMENSION (nslm,kjpindex)  
!$OMP THREADPRIVATE(nvan_mod_tab)
  
!! linarization coefficients of hydraulic conductivity K (hydrol_soil_coef)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: k        !! Hydraulic conductivity K for each soil layer
                                                                 !!  @tex $(mm d^{-1})$ @endtex
                                                                 !! DIM = (:,nslm)
!$OMP THREADPRIVATE(k)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: a        !! Slope in K=a*mc+b(:,nslm)
                                                                 !!  @tex $(mm d^{-1})$ @endtex
                                                                 !! DIM = (:,nslm)
!$OMP THREADPRIVATE(a)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: b        !! y-intercept in K=a*mc+b
                                                                 !!  @tex $(m^{3} m^{-3})$ @endtex
                                                                 !! DIM = (:,nslm)
!$OMP THREADPRIVATE(b)
!! linarization coefficients of hydraulic diffusivity D (hydrol_soil_coef)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: d        !! Diffusivity D for each soil layer
                                                                 !!  @tex $(mm^2 d^{-1})$ @endtex
                                                                 !! DIM = (:,nslm)
!$OMP THREADPRIVATE(d)
!! matrix coefficients (hydrol_soil_tridiag and hydrol_soil_setup), see De Rosnay (1999), p155-157
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: e        !! Left-hand tridiagonal matrix coefficients
!$OMP THREADPRIVATE(e)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: f        !! Left-hand tridiagonal matrix coefficients
!$OMP THREADPRIVATE(f)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: g1       !! Left-hand tridiagonal matrix coefficients
!$OMP THREADPRIVATE(g1)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: ep       !! Right-hand matrix coefficients
!$OMP THREADPRIVATE(ep)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: fp       !! Right-hand atrix coefficients
!$OMP THREADPRIVATE(fp)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: gp       !! Right-hand atrix coefficients
!$OMP THREADPRIVATE(gp)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: rhs      !! Right-hand system
!$OMP THREADPRIVATE(rhs)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: srhs     !! Temporarily stored rhs
!$OMP THREADPRIVATE(srhs)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: tmat             !! Left-hand tridiagonal matrix
!$OMP THREADPRIVATE(tmat)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: stmat            !! Temporarily stored tmat
  !$OMP THREADPRIVATE(stmat)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: water2infilt     !! Water to be infiltrated
                                                                         !! @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(water2infilt)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tmc              !! Total moisture content per soiltile 
                                                                         !!  @tex $(kg m^{-2})$ @endtex 
!$OMP THREADPRIVATE(tmc)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tmcr             !! Total moisture content at residual per soiltile 
                                                                         !!  @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(tmcr)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tmcs             !! Total moisture content at saturation per soiltile 
                                                                         !!  @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(tmcs)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tmcfc            !! Total moisture content at field capacity per soiltile 
                                                                         !!  @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(tmcfc)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tmcw             !! Total moisture content at wilting point per soiltile 
                                                                         !!  @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(tmcw)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tmc_litter       !! Total moisture in the litter per soiltile
                                                                         !!  @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(tmc_litter)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: tmc_litt_mea     !! Total moisture in the litter over the grid
                                                                         !!  @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(tmc_litt_mea)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tmc_litter_wilt  !! Total moisture of litter at wilt point per soiltile 
                                                                         !!  @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(tmc_litter_wilt)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tmc_litter_field !! Total moisture of litter at field cap. per soiltile 
                                                                         !!  @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(tmc_litter_field)
!!! A CHANGER DANS TOUT HYDROL: tmc_litter_res et sat ne devraient pas dependre de ji - tdo
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tmc_litter_res   !! Total moisture of litter at residual moisture per soiltile 
                                                                         !!  @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(tmc_litter_res)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tmc_litter_sat   !! Total moisture of litter at saturation per soiltile 
                                                                         !!  @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(tmc_litter_sat)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tmc_litter_awet  !! Total moisture of litter at mc_awet per soiltile 
                                                                         !!  @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(tmc_litter_awet)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tmc_litter_adry  !! Total moisture of litter at mc_adry per soiltile 
                                                                         !!  @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(tmc_litter_adry)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: tmc_litt_wet_mea !! Total moisture in the litter over the grid below which 
                                                                         !! albedo is fixed constant
                                                                         !!  @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(tmc_litt_wet_mea)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: tmc_litt_dry_mea !! Total moisture in the litter over the grid above which 
                                                                         !! albedo is constant
                                                                         !!  @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(tmc_litt_dry_mea)
  LOGICAL, SAVE                                      :: tmc_init_updated = .FALSE. !! Flag allowing to determine if tmc is initialized.
!$OMP THREADPRIVATE(tmc_init_updated)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: v1               !! Temporary variable (:)
!$OMP THREADPRIVATE(v1)

  !! par type de sol :
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: ru_ns            !! Surface runoff per soiltile
                                                                         !!  @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(ru_ns)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: dr_ns            !! Drainage per soiltile
                                                                         !!  @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(dr_ns)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tr_ns            !! Transpiration per soiltile
!$OMP THREADPRIVATE(tr_ns)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: vegetmax_soil    !! (:,nvm,nstm) percentage of each veg. type on each soil 
                                                                         !! of each grid point 
!$OMP THREADPRIVATE(vegetmax_soil)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: mc               !! Total volumetric water content at the calculation nodes
                                                                         !! (eg : liquid + frozen)
                                                                         !!  @tex $(m^{3} m^{-3})$ @endtex
!$OMP THREADPRIVATE(mc)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: root_mc_fc         !! Max Field Capacity moisture content, for layers with roots, in soil tile of irrig_st
                                                                         !!  @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(root_mc_fc)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION (:) :: nslm_root          !! max. layers defining the root zone
                                                                         !!  @tex $(layer)$ @endtex
!$OMP THREADPRIVATE(nslm_root)
  LOGICAL, ALLOCATABLE, SAVE, DIMENSION (:)        :: is_irrigated       !! is the PFT in an irrigated soil tile ?
                                                                         !! i.e pref_soil_veg(jv) == irrig_st
!$OMP THREADPRIVATE(is_irrigated)

   REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:) :: mc_read_prev       !! Soil moisture from file at previous timestep in the file
!$OMP THREADPRIVATE(mc_read_prev)
   REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:) :: mc_read_next       !! Soil moisture from file at next time step in the file
!$OMP THREADPRIVATE(mc_read_next)
   REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:) :: mc_read_current    !! For nudging, linear time interpolation bewteen mc_read_prev and mc_read_next
!$OMP THREADPRIVATE(mc_read_current)
   REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:) :: mask_mc_interp     !! Mask of valid data in soil moisture nudging file
!$OMP THREADPRIVATE(mask_mc_interp)
   REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)   :: tmc_aux            !! Temporary variable needed for the calculation of diag nudgincsm for nudging
!$OMP THREADPRIVATE(tmc_aux)
   REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)   :: snowdz_read_prev   !! snowdz read from file at previous timestep in the file [m]
!$OMP THREADPRIVATE(snowdz_read_prev)
   REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)   :: snowdz_read_next   !! snowdz read from file at next time step in the file [m]
!$OMP THREADPRIVATE(snowdz_read_next)
   REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)   :: snowrho_read_prev  !! snowrho read from file at previous timestep in the file (Kg/m^3)
!$OMP THREADPRIVATE(snowrho_read_prev)
   REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)   :: snowrho_read_next  !! snowrho read from file at next time step in the file (Kg/m^3)
!$OMP THREADPRIVATE(snowrho_read_next)
   REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)   :: snowtemp_read_prev !! snowtemp read from file at previous timestep in the file
!$OMP THREADPRIVATE(snowtemp_read_prev)
   REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)   :: snowtemp_read_next !! snowtemp read from file at next time step in the file
!$OMP THREADPRIVATE(snowtemp_read_next)
   REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)   :: mask_snow_interp   !! Mask of valid data in snow nudging file
!$OMP THREADPRIVATE(mask_snow_interp)

   REAL(r_std),ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: mcl              !! Liquid water content
                                                                         !!  @tex $(m^{3} m^{-3})$ @endtex
!$OMP THREADPRIVATE(mcl)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: soilmoist        !! (:,nslm) Mean of each soil layer's moisture
                                                                         !! across soiltiles
                                                                         !!  @tex $(kg m^{-2})$ @endtex 
!$OMP THREADPRIVATE(soilmoist)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: soilmoist_s      !! (:,nslm) Mean of each soil layer's moisture
                                                                         !! per soiltiles
                                                                         !!  @tex $(kg m^{-2})$ @endtex 
!$OMP THREADPRIVATE(soilmoist_s)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: soilmoist_liquid !! (:,nslm) Mean of each soil layer's liquid moisture
                                                                         !! across soiltiles
                                                                         !!  @tex $(kg m^{-2})$ @endtex 
!$OMP THREADPRIVATE(soilmoist_liquid)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: soil_wet_ns      !! Soil wetness above mcw (0-1, unitless)
!$OMP THREADPRIVATE(soil_wet_ns)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: soil_wet_litter  !! Soil wetness aove mvw in the litter (0-1, unitless)
!$OMP THREADPRIVATE(soil_wet_litter)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: qflux_ns         !! Diffusive water fluxes between soil layers
                                                                         !! (at lower interface)
!$OMP THREADPRIVATE(qflux_ns)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: check_top_ns     !! Diagnostic calculated in hydrol_diag_soil_flux
                                                                         !! (water balance residu of top soil layer)
!$OMP THREADPRIVATE(check_top_ns)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: profil_froz_hydro !! Frozen fraction for each hydrological soil layer
!$OMP THREADPRIVATE(profil_froz_hydro)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: profil_froz_hydro_ns  !! As  profil_froz_hydro per soiltile
!$OMP THREADPRIVATE(profil_froz_hydro_ns)

  ! useful flags or parameters when the natual peatland is activated via ok_peat_hydro
  LOGICAL,SAVE                                       :: ok_run2peat       !! Flag to activate the moving the runoff of other tiles to peat tile 
!$OMP THREADPRIVATE(ok_run2peat)
  LOGICAL,SAVE                                       :: ok_wt_ab         !! When flags on, there are above surface resevoir exist 
!$OMP THREADPRIVATE(ok_wt_ab)
  REAL(r_std), SAVE                                 :: max_wt_ab        !! The max height of the above surface sevoir if peat exist
!$OMP THREADPRIVATE(max_wt_ab)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: wt_ab            !! Above surface resevoir of soiltile4 for boreal natual peatland 
!$OMP THREADPRIVATE(wt_ab)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: run2peat         !! Sum of runoff from soiltile1-3, which were added into soiltile4 (boreal natual peatland)
!$OMP THREADPRIVATE(run2peat)
  
  ! useful flag or parameters of tides by Victoria Naipal (Mangro)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: wt_ab_tide       !! Above surface resevoir in case of tides
!$OMP THREADPRIVATE(wt_ab_tide)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: run2man          !! Sum of runoff from soiltile1-3, which were added into soiltile6 when tides is activated 
!$OMP THREADPRIVATE(run2man)
  
  ! flag or parameters useful when agri_peat is activated
  LOGICAL, SAVE                                      :: agri_drain=.FALSE.!! Flag to allow computing the drainage differently in case of agri_peat
!$OMP THREADPRIVATE(agri_drain)
  REAL(r_std), SAVE                                  :: drain_factor      !! New drain factor if agri_drain=True
!$OMP THREADPRIVATE(drain_factor)
  
  ! hydrological parameters for peatland
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)     :: kfact_peat        !! Factor to reduce Ks with depth (unitless) for peatland
!$OMP THREADPRIVATE(kfact_peat)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)     :: mc_lin_peat       !! 50 Vol. Wat. Contents to linearize K and D, for peatland
                                                                         !!  @tex $(m^{3} m^{-3})$ @endtex
!$OMP THREADPRIVATE(mc_lin_peat)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)   :: k_lin_peat        !! 50 values of unsaturated K, for each soil layer of peat
                                                                         !!  @tex $(mm d^{-1})$ @endtex
!$OMP THREADPRIVATE(k_lin_peat)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)   :: d_lin_peat        !! 50 values of unsaturated D, for each soil layer of peat
                                                                         !!  @tex $(mm d^{-1})$ @endtex
!$OMP THREADPRIVATE(d_lin_peat)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)   :: a_lin_peat        !! 50 values of unsaturated A, for each soil layer of peat
                                                                         !!  @tex $(mm d^{-1})$ @endtex
!$OMP THREADPRIVATE(a_lin_peat)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)   :: b_lin_peat        !! 50 values of unsaturated B, for each soil layer of peat
                                                                         !!  @tex $(mm d^{-1})$ @endtex
!$OMP THREADPRIVATE(b_lin_peat)
  
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)     :: mcs_grid          !! Soil moisture in saturation per grid, dim kjpindex
!$OMP THREADPRIVATE(mcs_grid)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)     :: mcw_grid          !! Soil moisture in wilting point per grid, dim kjpindex  
!$OMP THREADPRIVATE(mcw_grid)

CONTAINS

!! ================================================================================================================================
!! SUBROUTINE 	: hydrol_initialize
!!
!>\BRIEF         Allocate module variables, read from restart file or initialize with default values
!!
!! DESCRIPTION :
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE hydrol_initialize ( ks,             nvan,      avan,          mcr,              &
                                 mcs,            mcfc,      mcw,           kjit,             &
                                 kjpindex,       index,     rest_id,       njsc,             &
                                 soiltile,       veget,     veget_max,     frac_nobio,altmax,&
                                 humrel,         vegstress, drysoil_frac,  shumdiag_perma,   &
                                 qsintveg,       evap_bare_lim,            evap_bare_lim_ns, &
                                 snow,           snow_age,  snow_nobio_age,                  &
                                 snowrho,        snowtemp,  snowgrain,     snowdz,           &
                                 snowheat,       mc_layh,   mcl_layh,      soilmoist_out,    &
                                 mc_layh_s,      mcl_layh_s,soilmoist_out_s,mc_out,          &
                                 ksoil,          root_profile, us,         icetemp,          &
                                 icedz,          ice_sheet_mask,                             &
                                 psi_leaf,       psi_leaf_next,            psi_leaf_midday,  &
                                 psi_sto_leaf_save, psi_sto_wood_save,     psi_root_sup,     &
                                 psi_root_inf,   psi_soil_sup,             psi_soil_inf,     &
                                 psi_xylem_trunk,           psi_xylem_leaf,                  &
                                 psi_xylem_collar,  psi_sto_wood,          psi_sto_leaf,     &
                                 mc_i_sup,       mc_i_inf,  F_absorption, &
                                 wtp,liqwt_ratio,shumdiag_peat,shumdiag_croppeat, shumdiag_man)
                                


    !! 0. Variable and parameter declaration
    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                         :: kjit             !! Time step number 
    INTEGER(i_std), INTENT(in)                         :: kjpindex         !! Domain size
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)   :: index            !! Indeces of the points on the map
    INTEGER(i_std),INTENT (in)                         :: rest_id          !! Restart file identifier
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)   :: njsc             !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: ks               !! Hydraulic conductivity at saturation (mm {-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: nvan             !! Van Genuchten coeficients n (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: avan             !! Van Genuchten coeficients a (mm-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: mcr              !! Residual volumetric water content (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: mcs              !! Saturated volumetric water content (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: mcfc             !! Volumetric water content at field capacity (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: mcw              !! Volumetric water content at wilting point (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex,nstm), INTENT (in) :: soiltile         !! Fraction of each soil tile within vegtot (0-1, unitless)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: veget            !! Fraction of vegetation type           
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: veget_max        !! Max. fraction of vegetation type (LAI -> infty)
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (in)  :: frac_nobio    !! Fraction of ice, lakes, ...
    REAL(r_std),DIMENSION (kjpindex,nvm),INTENT(in)    :: altmax           !! Maximul active layer thickness (m). Be careful, here active means non frozen.
                                                                           !! Not related with the active soil carbon pool.
    !! 0.2 Output variables

    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)                :: humrel           !! Relative humidity
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)                :: vegstress        !! Veg. moisture stress (only for vegetation growth)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)                    :: drysoil_frac     !! function of litter wetness
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (out)               :: shumdiag_perma   !! Percent of porosity filled with water (mc/mcs) used for the thermal computations
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)                :: qsintveg         !! Water on vegetation due to interception
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)                     :: evap_bare_lim    !! Limitation factor for bare soil evaporation
    REAL(r_std),DIMENSION (kjpindex,nstm), INTENT(out)                :: evap_bare_lim_ns !! Limitation factor for bare soil evaporation
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)                    :: snow             !! Snow mass [Kg/m^2]
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)                    :: snow_age         !! Snow age
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (out)             :: snow_nobio_age   !! Snow age on ice, lakes, ...
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(out)              :: snowrho          !! Snow density (Kg/m^3)
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(out)              :: snowtemp         !! Snow temperature
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(out)              :: snowgrain        !! Snow grainsize 
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(out)              :: snowdz           !! Snow layer thickness [m]
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(out)              :: snowheat         !! Snow heat content
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (out)               :: mc_layh          !! Volumetric moisture content for each layer in hydrol (liquid+ice) m3/m3
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (out)               :: mcl_layh         !! Volumetric moisture content for each layer in hydrol (liquid) m3/m3
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (out)               :: soilmoist_out    !! Total soil moisture content for each layer in hydrol (liquid+ice), mm
    REAL(r_std),DIMENSION (kjpindex,nslm,nstm), INTENT (out)          :: mc_layh_s        !! Volumetric moisture content for each layer in hydrol (liquid+ice) m3/m3
    REAL(r_std),DIMENSION (kjpindex,nslm,nstm), INTENT (out)          :: mcl_layh_s       !! Volumetric moisture content for each layer in hydrol (liquid) m3/m3
    REAL(r_std),DIMENSION (kjpindex,nslm,nstm), INTENT (out)          :: soilmoist_out_s  !! Total soil moisture content for each layer in hydrol (liquid+ice), mm
    REAL(r_std),DIMENSION (kjpindex,nslm,nstm), INTENT (out)          :: mc_out
    REAL(r_std),DIMENSION (kjpindex,nslm,nstm), INTENT (out)          :: ksoil
    REAL(r_std),DIMENSION (kjpindex,nvm,nslm,nroot_prof), INTENT(out) :: root_profile     !! Normalized root mass/length fraction in each soil layer 
                                                                                          !! (0-1, unitless)
    REAL(r_std),DIMENSION (kjpindex,nvm,nstm,nslm), INTENT (out)      :: us               !! Water stress index for transpiration
                                                                                          !! (by soil layer and PFT) (0-1, unitless)
    REAL(r_std), DIMENSION (kjpindex,nice), INTENT(out)               :: icetemp          !! Ice temperature
    REAL(r_std), DIMENSION (kjpindex,nice), INTENT(out)               :: icedz            !! Ice layer thickness [m]
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (out)                 :: ice_sheet_mask   !! Ice sheet mask 

    REAL(r_std),DIMENSION (kjpindex,nvm,ncirc),       INTENT(out)     :: psi_leaf              !! Leaf Water potential (MPa)
    REAL(r_std),DIMENSION (kjpindex,nvm),       INTENT(out)           :: psi_leaf_next         !! Approximated Leaf Water potential at time step n+1 (MPa) (= psi_leaf when no stress)
    REAL(r_std),DIMENSION (kjpindex,nvm,ncirc),       INTENT(out)     :: psi_leaf_midday       !! Leaf Water potential at midday (MPa)
    REAL(r_std),DIMENSION (kjpindex,nvm,ncirc),       INTENT(out)     :: psi_sto_leaf_save     !! Leaf storage Water potential at time step n-1 (MPa)
    REAL(r_std),DIMENSION (kjpindex,nvm,ncirc),       INTENT(out)     :: psi_sto_wood_save     !! Wood storage Water potential at time step n-1 (MPa)
    REAL(r_std),DIMENSION(kjpindex,nvm,ncirc),        INTENT(out)     :: psi_root_sup          !! Superficial roots Water Potential (MPa)
    REAL(r_std),DIMENSION(kjpindex,nvm,ncirc),        INTENT(out)     :: psi_root_inf          !! Inferior roots Water Potential (MPa)
    REAL(r_std),DIMENSION(kjpindex,nstm),             INTENT(out)     :: psi_soil_sup          !! Superficial soil Water Potential (MPa)
    REAL(r_std),DIMENSION(kjpindex,nstm),             INTENT(out)     :: psi_soil_inf          !! Inferior soil Water Potential (MPa)

    REAL(r_std),DIMENSION(kjpindex,nvm,ncirc),        INTENT(out)     :: psi_xylem_trunk       !! Xylem (trunk level) Water Potential (MPa)
    REAL(r_std),DIMENSION(kjpindex,nvm,ncirc),        INTENT(out)     :: psi_xylem_leaf        !! Xylem (leaf level) Water Potential (MPa)
    REAL(r_std),DIMENSION(kjpindex,nvm,ncirc),        INTENT(out)     :: psi_xylem_collar      !! Xylem (collar level) Water Potential (MPa)
    REAL(r_std),DIMENSION(kjpindex,nvm,ncirc),        INTENT(out)     :: psi_sto_wood          !! Wood storage Water Potential (MPa)
    REAL(r_std),DIMENSION(kjpindex,nvm,ncirc),        INTENT(out)     :: psi_sto_leaf          !! Leaf storage Water Potential (MPa)
    REAL(r_std),DIMENSION(kjpindex,nvm,ncirc,nrp),    INTENT(out)     :: mc_i_sup              !! Water content at each node of the absorption muff of the superficial soil layer (m^3.m^-3)
    REAL(r_std),DIMENSION(kjpindex,nvm,ncirc,nrp),    INTENT(out)     :: mc_i_inf              !! Water content at each node of the absorption muff of the inferior soil layer (m^3.m^-3)
    REAL(r_std),DIMENSION(kjpindex,nvm),        INTENT(out)           :: F_absorption          !! Total root absorption flux (m^3/s)
    
    REAL(r_std), DIMENSION (kjpindex),INTENT(out)                     :: wtp                   !! water table position, computed per soil tile, when peat is activated
    REAL(r_std), DIMENSION (kjpindex),INTENT(out)                     :: liqwt_ratio           !! liquid water ratio, compared to the saturated one
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (out)               :: shumdiag_peat         !! soil humidity for peat, when peat is activated
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (out)               :: shumdiag_croppeat     !! soil humidity for crop peat, when agri_peat is activated
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (out)               :: shumdiag_man          !! soil humidity for mangrive, when tides is activated

    !! 0.4 Local variables
    INTEGER(i_std)                                                    :: jsl
    REAL(r_std),DIMENSION (kjpindex)                                  :: soilwetdummy          !! Temporary variable never used
!_ ================================================================================================================================
    !
    printlev_loc = get_printlev('hydrol')
    !
    CALL hydrol_init (ks, nvan, avan, mcr, mcs, mcfc, mcw, &
         kjit, kjpindex, index, rest_id, veget_max, frac_nobio, soiltile, &
         humrel, vegstress, snow, snow_age, &
         snow_nobio_age, qsintveg, snowdz, snowgrain, &
         snowrho, snowtemp, snowheat, &
         drysoil_frac, evap_bare_lim, evap_bare_lim_ns, mc_out, ksoil, &
         root_profile, us, icetemp, icedz, ice_sheet_mask, &
         psi_leaf,          psi_leaf_next,  psi_leaf_midday, psi_sto_leaf_save,  psi_sto_wood_save, &
         psi_root_sup,      psi_root_inf,   psi_soil_sup,       psi_soil_inf, &
         psi_xylem_trunk,    psi_xylem_leaf, &
         psi_xylem_collar,  psi_sto_wood,   psi_sto_leaf, &
         mc_i_sup,          mc_i_inf,       F_absorption, & 
         wtp,liqwt_ratio)

    CALL hydrol_var_init (ks, nvan, avan, mcr, mcs, mcfc, mcw, &
                          kjpindex, veget, veget_max, soiltile, njsc, altmax, &
         mx_eau_var, shumdiag_perma,                                          &
         drysoil_frac, qsintveg, mc_layh, mcl_layh, mc_layh_s, mcl_layh_s, &
         shumdiag_peat,shumdiag_croppeat, shumdiag_man) 

    !! Initialize hydrol_alma routine if the variables were not found in the restart file. This is done in the end of 
    !! hydrol_initialize so that all variables(humtot,..) that will be used are initialized.
    IF (ALL(tot_watveg_beg(:)==val_exp) .OR.  ALL(tot_watsoil_beg(:)==val_exp) .OR. ALL(snow_beg(:)==val_exp)) THEN
       ! The output variable soilwetdummy is not calculated at first call to hydrol_alma.
       CALL hydrol_alma(kjpindex, index, .TRUE., qsintveg, snow, soilwetdummy)
    END IF
    
    !! Calculate itopmax indicating the number of layers where the node is above 0.1m depth
    itopmax=1
    DO jsl = 1, nslm
       ! znh : depth of nodes
       IF (znh(jsl) <= 0.1) THEN
          itopmax=jsl
       END IF
    END DO
    IF (printlev>=3) WRITE(numout,*) "Number of layers where the node is above 0.1m depth: itopmax=",itopmax

    ! Copy soilmoist into a local variable to be sent to thermosoil
    soilmoist_out(:,:) = soilmoist(:,:)
    soilmoist_out_s(:,:,:) = soilmoist_s(:,:,:)

  END SUBROUTINE hydrol_initialize


!! ================================================================================================================================
!! SUBROUTINE 	: hydrol_main
!!
!>\BRIEF         
!!
!! DESCRIPTION :
!! - called every time step
!! - initialization and finalization part are not done in here
!!
!! - 1 computes snow  ==> explicitsnow
!! - 2 computes vegetations reservoirs  ==> hydrol_vegupd
!! - 3 computes canopy  ==> hydrol_canop
!! - 4 computes surface reservoir  ==> hydrol_flood
!! - 5 computes soil hydrology ==> hydrol_soil
!!
!! IMPORTANT NOTICE : The water fluxes are used in their integrated form, over the time step 
!! dt_sechiba, with a unit of kg m^{-2}.
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE hydrol_main (ks, nvan, avan, mcr, mcs, mcfc, mcw, kjit, kjpindex, &
       & index, indexveg, indexsoil, indexlayer, indexnslm, &
       & temp_sol_new, floodout, runoff, drainage, frac_nobio, totfrac_nobio, frac_snow_nobio, vevapwet, veget, veget_max, njsc, &
       & qsintmax, qsintveg, vevapnu, vevapsno, vevapflo, snow, snow_age, snow_nobio_age,  &
       & tot_melt, transpir, precip_rain, precip_snow, returnflow, reinfiltration, irrigation, &
       & humrel, vegstress, drysoil_frac, evapot, evapot_penm, evap_bare_lim, evap_bare_lim_ns, &
       & flood_frac, flood_res, z0m, &
       & shumdiag,shumdiag_perma, k_litt, litterhumdiag, soilcap, soiltile, fraclut, reinf_slope_soil, rest_id, hist_id, hist2_id,&
       & contfrac, stempdiag, &
       & temp_air, pb, u, v, tq_cdrag, swnet, pgflux, &
       & snowrho, snowtemp, snowgrain, snowdz, snowheat, snowliq, &
       & grndflux,gtemp,tot_bare_soil, &
       & lambda_snow,cgrnd_snow,dgrnd_snow,frac_snow_veg,temp_sol_add, &
       & lambda_ice,cgrnd_ice,dgrnd_ice, ice_sheet_mask, icetemp, icedz, snowmelt, zrainfall,&
       & mc_layh, mcl_layh, tmc_pft, drainage_pft, runoff_pft, swc_pft, soilmoist_out, mc_layh_s, mcl_layh_s, soilmoist_out_s,&
       & mc_out, e_frac,ksoil, mcs_hydrol, mcfc_hydrol, altmax, root_profile, &
       & root_depth, root_deficit, &
       & circ_class_biomass,  us, run_off_lic, run_off_lic_frac, circ_class_n, gsmean,&
       & psi_leaf, psi_leaf_next, psi_leaf_midday, psi_sto_leaf_save, psi_sto_wood_save, &
       & psi_root_sup, psi_root_inf, psi_soil_sup, psi_soil_inf, psi_xylem_trunk, psi_xylem_leaf, psi_xylem_collar, &
       & psi_sto_wood, psi_sto_leaf, &
       & mc_i_sup, mc_i_inf, F_absorption, lalo, z_array_out, &
       & qair, qsol_sat_new, rau, frac_evap_flood, frac_sub_snow, frac_evap_transp, cimean, assimtot, Rdtot_pft, vessel_loss,&
       & wtp,mc_peat_above,liqwt_ratio,shumdiag_peat,shumdiag_croppeat,mc_croppeat_above, &
       & shumdiag_man,mc_man_above ) 

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables
  
    INTEGER(i_std), INTENT(in)                         :: kjit             !! Time step number 
    INTEGER(i_std), INTENT(in)                         :: kjpindex         !! Domain size
    INTEGER(i_std),INTENT (in)                         :: rest_id,hist_id  !! _Restart_ file and _history_ file identifier
    INTEGER(i_std),INTENT (in)                         :: hist2_id         !! _history_ file 2 identifier
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)   :: index            !! Indeces of the points on the map
    INTEGER(i_std),DIMENSION (kjpindex*nvm), INTENT (in):: indexveg        !! Indeces of the points on the 3D map for veg
    INTEGER(i_std),DIMENSION (kjpindex*nstm), INTENT (in):: indexsoil      !! Indeces of the points on the 3D map for soil
    INTEGER(i_std),DIMENSION (kjpindex*nslm), INTENT (in):: indexlayer     !! Indeces of the points on the 3D map for soil layers
    INTEGER(i_std),DIMENSION (kjpindex*nslm), INTENT (in):: indexnslm      !! Indeces of the points on the 3D map for of diagnostic soil layers

    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: precip_rain      !! Rain precipitation
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: precip_snow      !! Snow precipitation
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: returnflow       !! Routed water which comes back into the soil (from the 
                                                                           !! bottom) 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: reinfiltration   !! Routed water which comes back into the soil (at the 
                                                                           !! top) 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: irrigation       !! Water from irrigation returning to soil moisture  
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: temp_sol_new     !! New soil temperature

    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)   :: njsc             !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless)
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (in) :: frac_nobio     !! Fraction of ice, lakes, ...
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: totfrac_nobio    !! Total fraction of ice+lakes+...
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT(in):: frac_snow_nobio  !! Snow cover fraction on non-vegeted area
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)       :: frac_snow_veg    !! Snow cover fraction on vegetation    
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: soilcap          !! Soil capacity
    REAL(r_std),DIMENSION (kjpindex,nstm), INTENT (in) :: soiltile         !! Fraction of each soil tile within vegtot (0-1, unitless)
    REAL(r_std),DIMENSION (kjpindex,nlut), INTENT (in) :: fraclut          !! Fraction of each landuse tile (0-1, unitless)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: vevapwet         !! Interception loss
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: veget            !! Fraction of vegetation type           
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: veget_max        !! Max. fraction of vegetation type (LAI -> infty)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: qsintmax         !! Maximum water on vegetation for interception
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: transpir         !! Transpiration (mm day-1)
    REAL(r_std),DIMENSION (kjpindex,nstm), INTENT (in) :: reinf_slope_soil !! Slope coef per soil tile
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: ks               !! Hydraulic conductivity at saturation (mm {-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: nvan             !! Van Genuchten coeficients n (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: avan             !! Van Genuchten coeficients a (mm-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: mcr              !! Residual volumetric water content (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: mcs              !! Saturated volumetric water content (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: mcfc             !! Volumetric water content at field capacity (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: mcw              !! Volumetric water content at wilting point (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: evapot           !! Soil Potential Evaporation
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: evapot_penm      !! Soil Potential Evaporation Correction
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: flood_frac       !! flood fraction
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: contfrac         !! Fraction of continent in the grid
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (in) :: stempdiag        !! Diagnostic temp profile from thermosoil
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)       :: temp_air         !! Air temperature
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)       :: u,v              !! Horizontal wind speed
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)       :: tq_cdrag         !! Surface drag coefficient (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)       :: pb               !! Surface pressure
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)       :: swnet            !! Net shortwave radiation
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)       :: pgflux           !! Net energy into snowpack
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)       :: gtemp            !! First soil layer temperature
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)       :: tot_bare_soil    !! Total evaporating bare soil fraction 
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)       :: lambda_snow      !! Coefficient of the linear extrapolation of surface temperature 
    REAL(r_std),DIMENSION (kjpindex,nsnow), INTENT (in):: cgrnd_snow       !! Integration coefficient for snow numerical scheme
    REAL(r_std),DIMENSION (kjpindex,nsnow), INTENT (in):: dgrnd_snow       !! Integration coefficient for snow numerical scheme
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)       :: lambda_ice       !! Coefficient of the linear extrapolation of surface ice temperature 
    REAL(r_std),DIMENSION (kjpindex,nice), INTENT (in) :: cgrnd_ice        !! Integration coefficient for ice numerical scheme
    REAL(r_std),DIMENSION (kjpindex,nice), INTENT (in) :: dgrnd_ice        !! Integration coefficient for ice numerical scheme
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)   :: ice_sheet_mask   !! Ice sheet mask
    REAL(r_std), DIMENSION (kjpindex,nvm),INTENT(in)   :: altmax           !! Maximul active layer thickness (m). Be careful, here active means non frozen.
                                                                           !! Not related with the active soil carbon pool.
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(in)      :: circ_class_biomass !! Biomass components of the model tree  
                                                                           !! within a circumference class
                                                                           !! class @tex $(g C ind^{-1})$ @endtex
    REAL(r_std), DIMENSION(kjpindex,nvm,ncirc),INTENT(in) :: circ_class_n
    REAL(r_std), DIMENSION(kjpindex), INTENT(in)       :: z0m              !! Surface roughness for momentum (m)
    REAL(r_std), DIMENSION (kjpindex),INTENT(in)       :: qair             !! Lowest level atmospheric air humidity (kg/kg)
    REAL(r_std), DIMENSION (kjpindex),INTENT(in)       :: qsol_sat_new     !! Surface saturated humidity (kg/kg)
    REAL(r_std), DIMENSION (kjpindex),INTENT(in)       :: rau              !! Air density (kg/m3)
    REAL(r_std), DIMENSION (kjpindex),INTENT(in)       :: frac_evap_flood  !! Fraction of gridcell that evaporates floodplains (-)
    REAL(r_std), DIMENSION (kjpindex),INTENT(in)       :: frac_sub_snow    !! Fraction of gridcell that sublimates snow (-)
    REAL(r_std), DIMENSION (kjpindex,nvm),INTENT(in)   :: frac_evap_transp !! Fraction of gridcell that transpires (-)
    REAL(r_std), DIMENSION (kjpindex,nvm),INTENT(in)   :: cimean           !! Mean leaf ci (ppm)
    REAL(r_std), DIMENSION (kjpindex,nvm),INTENT(in)   :: assimtot         !! Total assimilation per pft (mu mol m-2 s-1)
    REAL(r_std), DIMENSION (kjpindex,nvm),INTENT(in)   :: Rdtot_pft        !! Total respiration per pft (mu mol m-2 s-1) 
    REAL(r_std), DIMENSION(kjpindex,nvm,ncirc,nlevels_tot),INTENT(in)   :: z_array_out  


    !! 0.2 Output variables
    
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)                :: vegstress        !! Veg. moisture stress (only for vegetation growth)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)                    :: drysoil_frac     !! function of litter wetness
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (out)               :: shumdiag         !! Relative soil moisture in each soil layer 
                                                                                          !! with respect to (mcfc-mcw)
                                                                                          !! (unitless; can be out of 0-1)
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (out)               :: shumdiag_perma   !! Percent of porosity filled with water (mc/mcs) used for the thermal computations 
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)                    :: k_litt           !! litter approximate conductivity
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)                    :: litterhumdiag    !! litter humidity
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)                    :: tot_melt         !! Total melt    
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)                    :: floodout         !! Flux out of floodplains
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)                :: tmc_pft          !! Total soil water per PFT (mm/m2) 
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)                :: drainage_pft     !! Drainage per PFT (mm/m2)   
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)                :: runoff_pft       !! Runoff per PFT (mm/m2) 
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)                :: swc_pft          !! Relative Soil water content [tmcr:tmcs] per pft (-)
    REAL(r_std),DIMENSION (kjpindex,nslm,nstm), INTENT(out)           :: mc_out           !! Soil water content (copy of mc), which is need for 
                                                                                          !! the hydraulic architecture
    REAL(r_std),DIMENSION (kjpindex,nvm,nslm,nroot_prof), INTENT(out) :: root_profile     !! Normalized root mass/length fraction in each soil layer 
                                                                                          !! (0-1, unitless)
    REAL(r_std), DIMENSION (kjpindex,nvm,ndepths), INTENT(out)        :: root_depth       !! Node and interface numbers at which the deepest roots
                                                                                          !! occur (1 to nslm, unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)                    :: root_deficit     !! water deficit to reach field capacity of soil

    ! peatland, tides, without leak
    REAL(r_std), DIMENSION (kjpindex),INTENT(out)                     :: wtp              !! water table position
    REAL(r_std), DIMENSION (kjpindex),INTENT(out)                     :: liqwt_ratio      !! Liquid water ratio in peat, compared to the saturated one
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)                    :: mc_peat_above    !! liquid soil moisture content within the first 4 layers when peat is activated
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)                    :: mc_croppeat_above !! liquid soil moisture content within the first 4 layers when agri_peat is activated
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)                    :: mc_man_above      !! liquid soil moisture content within the first 4 layers when tides is activated
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (out)               :: shumdiag_peat     !! soil moisture content useful for the decomposition over peat
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (out)               :: shumdiag_croppeat !! soil moisture content useful for the decomposition when agri_peat is activated 
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (out)               :: shumdiag_man      !! soil moisture content useful for the decomposition when tides is activated 
    

    !! 0.3 Modified variables

    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(inout):: qsintveg                         !! Water on vegetation due to interception
    REAL(r_std),DIMENSION (kjpindex), INTENT(inout)    :: evap_bare_lim                    !! Limitation factor (beta) for bare soil evaporation 
    REAL(r_std),DIMENSION (kjpindex,nstm),INTENT(inout):: evap_bare_lim_ns                 !! Limitation factor (beta) for bare soil evaporation    
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(inout):: humrel                           !! Relative humidity
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)   :: vevapnu                          !! Bare soil evaporation
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)   :: vevapsno                         !! Snow evaporation
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)   :: vevapflo                         !! Floodplain evaporation
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)   :: flood_res                        !! flood reservoir estimate
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)   :: snow                             !! Snow mass [kg/m^2]
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)   :: snow_age                         !! Snow age
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (inout) :: snow_nobio_age              !! Snow age on ice, lakes, ...
    !! We consider that any water on the ice is snow and we only peforme a water balance to have consistency.
    !! The water balance is limite to + or - 10^6 so that accumulation is not endless

    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: runoff       !! Complete surface runoff
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: drainage     !! Drainage
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(inout) :: snowrho      !! Snow density (Kg/m^3)
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(inout) :: snowtemp     !! Snow temperature (K)
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(inout) :: snowgrain    !! Snow grainsize    
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(inout) :: snowdz       !! Snow layer thickness [m]
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(inout) :: snowheat     !! Snow heat content
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(out)   :: snowliq      !! Snow liquid content (m)
    REAL(r_std), DIMENSION (kjpindex), INTENT(out)         :: grndflux     !! Net flux into soil W/m2
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT(out)     :: mc_layh      !! Volumetric moisture content for each layer in hydrol(liquid + ice) [m3/m3)]
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT(out)     :: mcl_layh     !! Volumetric moisture content for each layer in hydrol(liquid) [m3/m3]
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT(out)     :: soilmoist_out!! Total soil moisture content for each layer in hydrol(liquid + ice) [mm]
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)       :: temp_sol_add !! additional surface temperature due to the melt of first layer
                                                                           !! at the present time-step @tex ($K$) @endtex
    REAL(r_std), DIMENSION (kjpindex), INTENT(out)         :: snowmelt     !! Snow melt [mm/dt_sechiba]
    REAL(r_std), DIMENSION (kjpindex), INTENT(out)         :: zrainfall    !! Rain precipitation on snow 
                                                                           !! @tex $(kg m^{-2})$ @endtex 
    REAL(r_std), DIMENSION (kjpindex,nice), INTENT(inout)  :: icetemp      !! Ice temperature
    REAL(r_std), DIMENSION (kjpindex,nice), INTENT(inout)  :: icedz        !! Ice layer thickness
    REAL(r_std),DIMENSION (kjpindex,nslm,nstm), INTENT(out)  :: mc_layh_s  !! Volumetric moisture content for each layer in hydrol(liquid + ice) [m3/m3)]
    REAL(r_std),DIMENSION (kjpindex,nslm,nstm), INTENT(out)  :: mcl_layh_s !! Volumetric moisture content for each layer in hydrol(liquid) [m3/m3]/
    REAL(r_std), DIMENSION (kjpindex,nslm,nstm), INTENT(inout) :: ksoil    !! Soil conductivity (a copy of k for each soil type) (mm/d)
    REAL(r_std),DIMENSION (kjpindex,nslm,nstm), INTENT(out)  :: soilmoist_out_s  !! Total soil moisture content for each layer in hydrol(liquid + ice) [mm]
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)            :: mcs_hydrol !! Saturated volumetric water content output to be used in stomate_soilcarbon
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)            :: mcfc_hydrol!! Volumetric water content at field capacity output to be used in stomate_soilcarbon
    REAL(r_std),DIMENSION (kjpindex,nvm,nstm,nslm), INTENT(inout) :: us    !! Water stress index for transpiration
                                                                           !! (by soil layer and PFT) (0-1, unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: run_off_lic       !! Contains calving, melting and liquid precipitation on continental ice
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: run_off_lic_frac  !! Contains cell fraction corresponding to run_off_lic
    REAL(r_std),DIMENSION (:,:,:), INTENT (out)              :: vessel_loss       !! Cavitation during timestep
    REAL(r_std), DIMENSION(kjpindex,nvm),            INTENT(inout)    :: gsmean                  !! Stomatal conductance 
    REAL(r_std), DIMENSION(kjpindex,nvm,ncirc),      INTENT(inout)    :: psi_leaf                !! Leaf Water Potential (MPa)
    REAL(r_std), DIMENSION(kjpindex,nvm),            INTENT(inout)    :: psi_leaf_next           !! Approximated Leaf Water Potential at time step n+1 (MPa) (=psi_leaf when no stress)
    REAL(r_std), DIMENSION(kjpindex,nvm,ncirc),      INTENT(inout)    :: psi_leaf_midday         !! Leaf Water Potential at midday (MPa)
    REAL(r_std), DIMENSION(kjpindex,nvm,ncirc),      INTENT(inout)    :: psi_sto_leaf_save       !! Leaf storage Water Potential at time step n-1 (MPa)
    REAL(r_std), DIMENSION(kjpindex,nvm,ncirc),      INTENT(inout)    :: psi_sto_wood_save       !! Wood storage Water Potential at time step n-1 (MPa)
    REAL(r_std), DIMENSION(kjpindex,nvm,nslm,nstm),  INTENT(inout)    :: e_frac                  !! Fraction of water transpired supplied by individual layers (no units) 
    REAL(r_std), DIMENSION(kjpindex,nvm,ncirc),      INTENT(inout)    :: psi_root_sup            !! Superficial roots Water Potential (MPa) 
    REAL(r_std), DIMENSION(kjpindex,nvm,ncirc),      INTENT(inout)    :: psi_root_inf            !! Inferior roots Water Potential (MPa)
    REAL(r_std), DIMENSION(kjpindex,nstm),           INTENT(inout)    :: psi_soil_sup            !! Superficial soil Water Potential (MPa)
    REAL(r_std), DIMENSION(kjpindex,nstm),           INTENT(inout)    :: psi_soil_inf            !! Inferior soil Water Potential (MPa)
    REAL(r_std), DIMENSION(kjpindex,nvm,ncirc),      INTENT(inout)    :: psi_xylem_trunk         !! Xylem (Trunk level) Water Potential (MPa)
    REAL(r_std), DIMENSION(kjpindex,nvm,ncirc),      INTENT(inout)    :: psi_xylem_leaf          !! Xylem (Leaf level) Water Potential (MPa)
    REAL(r_std), DIMENSION(kjpindex,nvm,ncirc),      INTENT(inout)    :: psi_xylem_collar        !! Xylem (Collar level) Water Potential (MPa)
    REAL(r_std), DIMENSION(kjpindex,nvm,ncirc),      INTENT(inout)    :: psi_sto_wood            !! Wood storage Water Potential (MPa)
    REAL(r_std), DIMENSION(kjpindex,nvm,ncirc),      INTENT(inout)    :: psi_sto_leaf            !! Leaf storage Water Potential (MPa)
    REAL(r_std), DIMENSION(kjpindex,nvm,ncirc,nrp),  INTENT(inout)    :: mc_i_sup                !! Water content at each node of the absorption muff of the superficial soil layer (m^3.m^-3)
    REAL(r_std), DIMENSION(kjpindex,nvm,ncirc,nrp),  INTENT(inout)    :: mc_i_inf                !! Water content at each node of the absorption  muff of the inferior soil layer (m^3.m^-3)
    REAL(r_std), DIMENSION(kjpindex,nvm),            INTENT(inout)    :: F_absorption            !! Total root absorption flux (m^3/s)
    REAL(r_std), DIMENSION(kjpindex,2),              INTENT(in)       :: lalo                    !! Lat/Lon


     
    !! 0.4 Local variables

    INTEGER(i_std)                                     :: jst              !! Index of soil tiles (unitless, 1-3)
    INTEGER(i_std)                                     :: jsl              !! Index of soil layers (unitless)
    INTEGER(i_std)                                     :: ji, jv
    LOGICAL(i_std), SAVE                               :: first_hydrol_main=.TRUE.
    REAL(r_std),DIMENSION (kjpindex)                   :: soilwet          !! A temporary diagnostic of soil wetness
    REAL(r_std),DIMENSION (kjpindex)                   :: snowdepth_diag   !! Depth of snow layer containing default values, only for diagnostics
    REAL(r_std),DIMENSION (kjpindex, nsnow)            :: snowdz_diag      !! Depth of snow layer on all layers containing default values, 
                                                                           !! only for diagnostics [m]
    REAL(r_std),DIMENSION (kjpindex)                   :: njsc_tmp         !! Temporary REAL value for njsc to write it
    REAL(r_std), DIMENSION (kjpindex,nstm)             :: tmc_top          !! Moisture content in the itopmax upper layers, per tile
    REAL(r_std), DIMENSION (kjpindex)                  :: humtot_top       !! Moisture content in the itopmax upper layers, for diagnistics
    REAL(r_std), DIMENSION(kjpindex)                   :: histvar          !! Temporary variable when computations are needed
    REAL(r_std), DIMENSION (kjpindex,nvm)              :: frac_bare        !! Fraction(of veget_max) of bare soil in each vegetation type
    INTEGER(i_std), DIMENSION(kjpindex*imax)           :: mc_lin_axis_index
    REAL(r_std), DIMENSION(kjpindex)                   :: twbr             !! Grid-cell mean of TWBR Total Water Budget Residu[kg/m2/dt]
    REAL(r_std), DIMENSION(kjpindex,nvm)               :: F_storage        !! Water supplied by storage [kg/m2/dt]
    REAL(r_std),DIMENSION (kjpindex,nslm)              :: land_root_profile!! To ouput the grid-cell mean of nroot
    REAL(r_std),DIMENSION (kjpindex,nslm)              :: land_dlh         !! To ouput the soil layer thickness on all grid points [m]
    REAL(r_std),DIMENSION (kjpindex,nslm)              :: land_mcs         !! To ouput the mean of mcs
    REAL(r_std),DIMENSION (kjpindex,nslm)              :: land_mcfc        !! To ouput the mean of mcfc 
    REAL(r_std),DIMENSION (kjpindex,nslm)              :: land_mcw         !! To ouput the mean of mcw 
    REAL(r_std),DIMENSION (kjpindex,nslm)              :: land_mcr         !! To ouput the mean of mcr 
    REAL(r_std),DIMENSION (kjpindex)                   :: land_tmcs        !! To ouput the grid-cell mean of tmcs 
    REAL(r_std),DIMENSION (kjpindex)                   :: land_tmcfc       !! To ouput the grid-cell mean of tmcfc
    REAL(r_std),DIMENSION (kjpindex)                   :: drain_upd        !! Change in drainage due to decrease in vegtot
                                                                           !! on mc [kg/m2/dt]
    REAL(r_std),DIMENSION (kjpindex)                   :: runoff_upd       !! Change in runoff due to decrease in vegtot
                                                                           !! on water2infilt[kg/m2/dt]
    REAL(r_std),DIMENSION (kjpindex)                   :: mrsow            !! Soil wetness above wilting point for CMIP6 (humtot-WP)/(SAT-WP)
    REAL(r_std), DIMENSION (kjpindex,nlut)             :: humtot_lut       !! Moisture content on landuse tiles, for diagnostics
    REAL(r_std), DIMENSION (kjpindex,nlut)             :: humtot_top_lut   !! Moisture content in upper layers on landuse tiles, for diagnostics
    REAL(r_std), DIMENSION (kjpindex,nlut)             :: mrro_lut         !! Total runoff from landuse tiles, for diagnostics
    REAL(r_std),DIMENSION (kjpindex)                   :: mcw_tmp, mcs_tmp !! combining according to soil tile
    REAL(r_std),DIMENSION (nstm)                       :: veg_tmp          !! Temporary variable for sum of vegetation per soil tile

!_ ================================================================================================================================
    !! 1. Update vegtot_old and recalculate vegtot
    vegtot_old(:) = vegtot(:)

    DO ji = 1, kjpindex
       vegtot(ji) = SUM(veget_max(ji,:))
    ENDDO

    !! Calculate rsoil_scale_tile: Average of rsoil_scale for each tile depending on veget_max
    !!
    !! To save computing time, this is done only when veget_max could have changed:
    !! only at start, once per year or if DGVM is activated
    IF ((first_hydrol_main .OR. FirstTsYear) .OR. OK_DGVM) THEN
       rsoil_scale_tile(:,:)=0.0
       DO ji = 1, kjpindex
          veg_tmp(:) = 0.0
          DO jv=1,nvm
             jst = pref_soil_veg(jv) 
             rsoil_scale_tile(ji,jst) = rsoil_scale_tile(ji,jst) + veget_max(ji,jv)*rsoil_scale(jv)
             veg_tmp(jst) = veg_tmp(jst)+veget_max(ji,jv)
          END DO
          DO jst=1,nstm
             IF (veg_tmp(jst) > min_sechiba) THEN
                rsoil_scale_tile(ji,jst) = rsoil_scale_tile(ji,jst)/veg_tmp(jst)
             ELSE
                ! Too small vegetaion on this soil tile, set variable to 0
                rsoil_scale_tile(ji,jst) = 0.0
             END IF
          END DO
       END DO

       ! Check for too high values
       DO ji = 1, kjpindex
          IF (MAXVAL(rsoil_scale_tile(ji,:)) .GT. MAXVAL(rsoil_scale(:))) THEN
             WRITE(numout,*) 'rsoil_scale_tile(ji,:)=', rsoil_scale_tile(ji,:)
             WRITE(numout,*) 'rsoil_scale(:)=', rsoil_scale(:)
             CALL ipslerr_p(3, 'hydrol_main',&
                  'At least one value for rsoil_scale_tile is higher than the biggest value in RSOIL_SCALE','','')
          END IF
       END DO
       CALL xios_orchidee_send_field("rsoil_scale_tile",rsoil_scale_tile)
       first_hydrol_main=.FALSE.
    END IF
    
    !! 2. Applay nudging for soil moisture and/or snow variables

    ! For soil moisture, here only read and interpolate the soil moisture from file to current time step. 
    ! The values will be applayed in hydrol_soil after the soil moisture has been updated. 
    IF (ok_nudge_mc) THEN
       CALL hydrol_nudge_mc_read(kjit)
    END IF

    ! Read, interpolate and applay nudging of snow variables
    IF ( ok_nudge_snow) THEN
     CALL hydrol_nudge_snow(kjit, kjpindex, snowdz, snowrho, snowtemp )
    END IF
    !
    !! 3. Shared time step
    IF ((printlev>=3).OR.(printlev_loc>=4)) THEN
       WRITE (numout,*) 'hydrol pas de temps = ',dt_sechiba
       CALL flush(numout)
    ENDIF
    !
    ! Loop on soiltiles to compute the variables (ji,jst)
    DO jv=1,nvm
       DO ji = 1, kjpindex
          IF (soiltile(ji,pref_soil_veg(jv)) .LE. min_sechiba) THEN
             tmc_pft(ji,jv) = 0.0
             swc_pft(ji,jv) = 0.0
          ELSE
             tmc_pft(ji,jv)      = MAX(tmc(ji,pref_soil_veg(jv)),tmcr(ji,pref_soil_veg(jv)))
             swc_pft(ji,jv)      = MIN(un, (MAX(zero, (( tmc(ji,pref_soil_veg(jv)) - tmcr(ji,pref_soil_veg(jv)) ) &
                  / ( tmcs(ji,pref_soil_veg(jv)) - tmcr(ji,pref_soil_veg(jv)))))))
          ENDIF
       ENDDO
    ENDDO
    !
    !! Calculate hydraulic architecture, by Tuzet et al. 2017
    IF (ok_hydrol_arch) THEN
       CALL hydrol_hydraulic_arch_tuzet_calc(kjit, kjpindex, ks, nvan, avan, transpir, mc_out, mcw, root_profile, veget, veget_max, &
            njsc, soiltile, circ_class_n, circ_class_biomass, u, v, tq_cdrag, gsmean, pb, temp_air, qair, qsol_sat_new, rau,    & 
            frac_evap_flood, frac_sub_snow, frac_evap_transp, cimean, assimtot, Rdtot_pft, lalo, z_array_out, &
            psi_leaf, psi_leaf_next, psi_leaf_midday, psi_sto_leaf_save, psi_sto_wood_save, e_frac, psi_root_sup, psi_root_inf, &
            psi_soil_sup, psi_soil_inf, psi_xylem_trunk, psi_xylem_leaf, psi_xylem_collar, psi_sto_wood, psi_sto_leaf, &
            mc_i_sup, mc_i_inf, F_absorption, vessel_loss)
    ENDIF

    ! 
    !! 3.1 Calculate snow processes with explicit snow model
    CALL explicitsnow_main(kjpindex,    precip_rain,  precip_snow,   temp_air,    pb,       &
         u,           v,            temp_sol_new,  soilcap,     pgflux,   &
         frac_nobio,  totfrac_nobio,frac_snow_nobio,gtemp,                &
         lambda_snow, cgrnd_snow,   dgrnd_snow,    contfrac,    z0m,      & 
         lambda_ice,  cgrnd_ice,    dgrnd_ice,     ice_sheet_mask,        &
         vevapsno,    snow_age,     snow_nobio_age, snowrho,              & !inout
         snowgrain,   snowdz,                                             & !inout
         snowtemp,     snowheat,    snow,                                 & ! inout
         temp_sol_add,icetemp,      icedz,                                &
         snowliq,     subsnownobio, grndflux,      snowmelt,    tot_melt, &
         subsinksoil, zrainfall,    frac_snow_veg, veget,       veget_max,&
         run_off_lic, run_off_lic_frac)            

    !We normalize per vegtot for grid boxes that are not 100% nobio(ice).       
    WHERE(vegtot(:) .NE. 0.)
       subsinksoil(:)=subsinksoil(:)/vegtot(:)
    ENDWHERE       
    !
    !! 3.2 computes vegetations reservoirs  ==>hydrol_vegupd
    CALL hydrol_vegupd(kjpindex, veget, veget_max, soiltile, qsintveg, frac_bare, drain_upd, runoff_upd)
    !
    !! 3.3 computes canopy  ==>hydrol_canop
    ! additional variables could be added below for the leak case, but not added here
    CALL hydrol_canop(kjpindex, precip_rain, vevapwet, veget_max, veget, qsintmax, qsintveg,precisol,tot_melt, frac_snow_veg)

    !
    !! 3.4 computes surface reservoir  ==>hydrol_flood
    CALL hydrol_flood(kjpindex,  vevapflo, flood_frac, flood_res, floodout)

    !
    !! 3.5 computes soil hydrology ==>hydrol_soil

    CALL hydrol_soil(ks, nvan, avan, mcr, mcs, mcfc, mcw, &
         kjpindex, veget, veget_max, soiltile, njsc, reinf_slope_soil,  &
         transpir, vevapnu, evapot, evapot_penm, runoff, drainage, & 
         returnflow, reinfiltration, irrigation, &
         tot_melt,evap_bare_lim,evap_bare_lim_ns, shumdiag, shumdiag_perma, &
         k_litt, litterhumdiag, humrel, vegstress, drysoil_frac,&
         stempdiag,snow,snowdz, tot_bare_soil,  u, v, tq_cdrag, &
         mc_layh, mcl_layh, mc_layh_s, mcl_layh_s, e_frac, ksoil, &
         altmax, root_profile, root_depth, root_deficit, circ_class_biomass, &
         us, precip_rain, totfrac_nobio, frac_snow_nobio, F_absorption, &
         wtp,mc_peat_above,liqwt_ratio,shumdiag_peat,shumdiag_croppeat,mc_croppeat_above,&
         shumdiag_man,mc_man_above)
 
    ! The update fluxes come from hydrol_vegupd
    ! drainage and runoff: for the whole grid, by combining different soiltiles in the grid
    ! draindage: computed based ed dr_ns, 
    ! runoff: computed based on ru_ns
    drainage(:) =  drainage(:) +  drain_upd(:)
    runoff(:) =  runoff(:) +  runoff_upd(:)

    DO jv=1,nvm
       ! dr_ns: per soil tile
       ! ru_ns: per soil tile
       drainage_pft(:,jv) = dr_ns(:,pref_soil_veg(jv))
       runoff_pft(:,jv) = ru_ns(:,pref_soil_veg(jv))
    ENDDO

    !! 12.10 Copy mc to mc_out to return from hydrol_main. It's needed in hydrolaulic_arch, when we do not use
    !! soil to root resistance
    mc_out(:,:,:) = mc(:,:,:)

    !! 4 write out file  ==> hydrol_alma/histwrite(*)
    !
    ! If we use the ALMA standards
    CALL hydrol_alma(kjpindex, index, .FALSE., qsintveg, snow, soilwet)
    
    ! Calculate the moisture in the upper itopmax layers corresponding to 0.1m (humtot_top): 
    ! For ORCHIDEE with nslm=11 and zmaxh=2, itopmax=6. 
        ! We compute tmc_top as tmc but only for the first itopmax layers. Then we compute a humtot with this variable.
    DO jst=1,nstm
       DO ji=1,kjpindex
          tmc_top(ji,jst) = dz(2) * ( trois*mc(ji,1,jst) + mc(ji,2,jst) )/huit
          DO jsl = 2, itopmax
             tmc_top(ji,jst) = tmc_top(ji,jst) + dz(jsl) * (trois*mc(ji,jsl,jst)+mc(ji,jsl-1,jst))/huit &
                  + dz(jsl+1) * (trois*mc(ji,jsl,jst)+mc(ji,jsl+1,jst))/huit
          ENDDO
       ENDDO
    ENDDO
 
    ! We average the values of each soiltile and multiply by vegtot to transform to a grid-cell mean
    humtot_top(:) = zero
    DO jst=1,nstm
       DO ji=1,kjpindex
          humtot_top(ji) = humtot_top(ji) + soiltile(ji,jst) * tmc_top(ji,jst) * vegtot(ji)
       ENDDO
    ENDDO

    ! Calculate the Total Water Budget Residu (in kg/m2 over dt_sechiba)
    ! All the delstocks and fluxes below are averaged over the mesh
    ! Does not include the routing reservoirs, although the flux to/from routing are integrated
    IF (ok_hydrol_arch) THEN 
       DO ji=1,kjpindex
          twbr(ji) = (delsoilmoist(ji) + delintercept(ji) + delswe(ji)) &
               - ( precip_rain(ji) + precip_snow(ji) + irrigation(ji) + floodout(ji) &
               + returnflow(ji) + reinfiltration(ji) ) &
               + ( runoff(ji) + drainage(ji) + SUM(vevapwet(ji,:)) &
               + SUM(F_absorption(ji,:)*(dt_sechiba * kilo_to_unit* soiltile(ji,pref_soil_veg(:)) * vegtot(ji))) &
               + vevapnu(ji) + vevapsno(ji) + vevapflo(ji) )
               ! + SUM(transpir(ji,:)) + vevapnu(ji) + vevapsno(ji) + vevapflo(ji) )
       ENDDO
    ELSE 
       DO ji=1,kjpindex
          twbr(ji) = (delsoilmoist(ji) + delintercept(ji) + delswe(ji)) &
               - ( precip_rain(ji) + precip_snow(ji) + irrigation(ji) + floodout(ji) &
               + returnflow(ji) + reinfiltration(ji) ) &
               + ( runoff(ji) + drainage(ji) + SUM(vevapwet(ji,:)) &
               + SUM(transpir(ji,:)) + vevapnu(ji) + vevapsno(ji) + vevapflo(ji) )
       ENDDO
    ENDIF

    ! Calculating flux from storage
    F_storage(:,:) = zero
    IF (ok_hydrol_arch) THEN
       DO ji=1,kjpindex
         DO jv =1,nvm
             F_storage(ji,jv) = - F_absorption(ji,jv) * dt_sechiba * kilo_to_unit + &
                         transpir(ji,jv)
             !F_storage(ji,jv) = - F_absorption(ji,jv)*(dt_sechiba * kilo_to_unit* soiltile(ji,pref_soil_veg(jv)) * vegtot(ji)) + &
             !            transpir(ji,jv)
         ENDDO
       ENDDO
    ENDIF

    ! Transform unit from kg/m2/dt to kg/m2/s (or mm/s)
    CALL xios_orchidee_send_field("twbr",twbr/dt_sechiba)
    IF (ok_hydrol_arch) CALL xios_orchidee_send_field("F_storage",F_storage/dt_sechiba)
    CALL xios_orchidee_send_field("undermcr",undermcr) ! nb of tiles undermcr at end of timestep

    ! Calculate land_root_profile : grid-cell mean of the structural root_profile 
    ! Do not treat PFT1 because it has no roots
    land_root_profile(:,:) = zero
    DO jsl=1,nslm
       DO jv=2,nvm
          DO ji=1,kjpindex
               IF ( vegtot(ji) > min_sechiba ) THEN 
               land_root_profile(ji,jsl) = land_root_profile(ji,jsl) + veget_max(ji,jv) * root_profile(ji,jv,jsl,istruc) / vegtot(ji) 
            END IF
          END DO
       ENDDO
    ENDDO

    CALL xios_orchidee_send_field("land_root_profile",land_root_profile)   

    DO jsl=1,nslm
       land_dlh(:,jsl)=dlh(jsl)
    ENDDO
    CALL xios_orchidee_send_field("dlh",land_dlh)

    ! Particular soil moisture values, spatially averaged over the grid-cell
    ! (a) total SM in kg/m2
    !     we average the total values of each soiltile and multiply by vegtot to transform to a grid-cell mean (over total land)
    !
    ! peatland
    ! tmcs is initialized differently for peat and nonpeat in hydrol_var_init, by using mcs and mcs_peat
    ! tmcs is not modified any more after the initialization
    ! land_tmcs is good for the case that the peatland is included in the grid
    land_tmcs(:) = zero
    land_tmcfc(:) = zero
    DO jst=1,nstm
       DO ji=1,kjpindex
          land_tmcs(ji) = land_tmcs(ji) + soiltile(ji,jst) * tmcs(ji,jst) * vegtot(ji)
          land_tmcfc(ji) = land_tmcfc(ji) + soiltile(ji,jst) * tmcfc(ji,jst) * vegtot(ji)
       ENDDO
    ENDDO
    CALL xios_orchidee_send_field("tmcs",land_tmcs) ! in kg/m2
    CALL xios_orchidee_send_field("tmcfc",land_tmcfc) ! in kg/m2

    ! (b) volumetric moisture content by layers in m3/m3
    !     mcs etc are identical in all layers (no normalization by vegtot to be comparable to mc)
    ! keep the trunk original, because land_mcs is not really useful in the computation
    ! 
    ! land_mcs etc contains only non peatland 
    ! if we want to include peatland in the land_mcs, 
    ! we have to use soiltie to combine both
    DO jsl=1,nslm
       land_mcs(:,jsl) = mcs(:)
       land_mcfc(:,jsl) = mcfc(:)
       land_mcw(:,jsl) = mcw(:)
       land_mcr(:,jsl) = mcr(:)
    ENDDO
    CALL xios_orchidee_send_field("mcs",land_mcs) ! in m3/m3
    CALL xios_orchidee_send_field("mcfc",land_mcfc) ! in m3/m3 
    CALL xios_orchidee_send_field("mcw",land_mcw) ! in m3/m3 
    CALL xios_orchidee_send_field("mcr",land_mcr) ! in m3/m3 
    CALL xios_orchidee_send_field("mc_layh",mc_layh)
    CALL xios_orchidee_send_field("mcl_layh",mcl_layh) 

    IF (ok_peat_hydro) THEN
       CALL xios_orchidee_send_field("run2peat",run2peat/dt_sechiba)         
       CALL xios_orchidee_send_field("wt_ab",wt_ab)
       CALL xios_orchidee_send_field("wtp",wtp)
       CALL xios_orchidee_send_field("mcl",mcl)
    ENDIF

    IF (tides) THEN
       CALL xios_orchidee_send_field("wt_ab_tide",wt_ab_tide)
       CALL xios_orchidee_send_field("run2man",run2man/dt_sechiba) 
    ENDIF
      
    CALL xios_orchidee_send_field("water2infilt",water2infilt)   
    CALL xios_orchidee_send_field("mc",mc)
    CALL xios_orchidee_send_field("kfact_root",kfact_root)
    CALL xios_orchidee_send_field("rootsink",rootsink)
    CALL xios_orchidee_send_field("vegetmax_soil",vegetmax_soil)
    CALL xios_orchidee_send_field("evapnu_soil",ae_ns/dt_sechiba)
    CALL xios_orchidee_send_field("drainage_soil",dr_ns/dt_sechiba)
    CALL xios_orchidee_send_field("transpir_soil",tr_ns/dt_sechiba)
    CALL xios_orchidee_send_field("runoff_soil",ru_ns/dt_sechiba)
    CALL xios_orchidee_send_field("humrel",humrel)     
    CALL xios_orchidee_send_field("drainage",drainage/dt_sechiba) ! [kg m-2 s-1]
    CALL xios_orchidee_send_field("runoff",runoff/dt_sechiba) ! [kg m-2 s-1]
    CALL xios_orchidee_send_field("precisol",precisol/dt_sechiba)
    CALL xios_orchidee_send_field("throughfall",throughfall/dt_sechiba)
    CALL xios_orchidee_send_field("precip_rain",precip_rain/dt_sechiba)
    CALL xios_orchidee_send_field("precip_snow",precip_snow/dt_sechiba)
    CALL xios_orchidee_send_field("qsintmax",qsintmax)
    CALL xios_orchidee_send_field("qsintveg",qsintveg)
    CALL xios_orchidee_send_field("qsintveg_tot",SUM(qsintveg(:,:),dim=2))
    histvar(:)=(precip_rain(:)-SUM(throughfall(:,:),dim=2))
    CALL xios_orchidee_send_field("prveg",histvar/dt_sechiba)

    IF ( do_floodplains ) THEN
       CALL xios_orchidee_send_field("floodout",floodout/dt_sechiba)
    END IF

    CALL xios_orchidee_send_field("snowmelt",snowmelt/dt_sechiba)
    CALL xios_orchidee_send_field("tot_melt",tot_melt/dt_sechiba)

    CALL xios_orchidee_send_field("soilmoist",soilmoist)
    CALL xios_orchidee_send_field("soilmoist_liquid",soilmoist_liquid)
    CALL xios_orchidee_send_field("shumdiag_perma",shumdiag_perma)
    CALL xios_orchidee_send_field("humtot_frozen",SUM(soilmoist(:,:),2)-SUM(soilmoist_liquid(:,:),2))

    CALL xios_orchidee_send_field("tmc",tmc)
    CALL xios_orchidee_send_field("humtot",humtot)
    CALL xios_orchidee_send_field("humtot_top",humtot_top)

    ! For the soil wetness above wilting point for CMIP6 (mrsow)
    ! modified by taking into account the soiltiles
    ! But mrsow is not really useful in the current computation in trunk
    DO ji=1,kjpindex
       !
       IF( .NOT. ok_peat_hydro .AND. (nstm == 3)) THEN
          mcw_tmp(ji) = mcw(ji)
          mcs_tmp(ji) = mcs(ji)
       ELSE IF(ok_peat_hydro .AND. (nstm == 4)) THEN
          mcw_tmp(ji) = mcw(ji)*(un-soiltile(ji,4))+mcw_peat*soiltile(ji,4)
          mcs_tmp(ji) = mcs(ji)*(un-soiltile(ji,4))+mcs_peat*soiltile(ji,4)
       ELSE IF(ok_peat_hydro .AND. (nstm == 5)) THEN
          mcw_tmp(ji) = mcw(ji)*(un-soiltile(ji,4)-soiltile(ji,5))+mcw_peat*(soiltile(ji,5)+ soiltile(ji,5))
          mcs_tmp(ji) = mcs(ji)*(un-soiltile(ji,4)-soiltile(ji,5))+mcs_peat*(soiltile(ji,5)+ soiltile(ji,5))
        ELSE IF(ok_peat_hydro .AND. (nstm == 6)) THEN
          mcw_tmp(ji) = mcw(ji)*(un-soiltile(ji,4)-soiltile(ji,5)-soiltile(ji,6))+mcw_peat*(soiltile(ji,5)+ soiltile(ji,5)-soiltile(ji,6))
          mcs_tmp(ji) = mcs(ji)*(un-soiltile(ji,4)-soiltile(ji,5)-soiltile(ji,6))+mcs_peat*(soiltile(ji,5)+ soiltile(ji,5)-soiltile(ji,6))
       ENDIF
    ENDDO
    mrsow(:) = MAX( zero,humtot(:) - zmaxh*mille*mcw_tmp(:) ) &
         / ( zmaxh*mille*( mcs_tmp(:) - mcw_tmp(:) ) )
    CALL xios_orchidee_send_field("mrsow",mrsow)

   ! Output irrigation related variables
    CALL xios_orchidee_send_field("root_deficit", root_deficit)
    CALL xios_orchidee_send_field("root_mc_fc", root_mc_fc)
    
    ! Prepare diagnostic snow variables
    !  Add XIOS default value where no snow
    DO ji=1,kjpindex
       IF (snow(ji) > 0) THEN
          snowdz_diag(ji,:) = snowdz(ji,:)
          snowdepth_diag(ji) = SUM(snowdz(ji,:))*(1-totfrac_nobio(ji))*frac_snow_veg(ji)
       ELSE
          snowdz_diag(ji,:) = xios_default_val
          snowdepth_diag(ji) = xios_default_val             
       END IF
    END DO
    CALL xios_orchidee_send_field("snowdz",snowdz_diag)
    CALL xios_orchidee_send_field("snowdepth",snowdepth_diag)

    CALL xios_orchidee_send_field("frac_bare",frac_bare)
    CALL xios_orchidee_send_field("soilwet",soilwet)
    CALL xios_orchidee_send_field("delsoilmoist",delsoilmoist)
    CALL xios_orchidee_send_field("delswe",delswe)
    CALL xios_orchidee_send_field("delintercept",delintercept)  

    IF (ok_freeze_cwrr) THEN
       CALL xios_orchidee_send_field("profil_froz_hydro",profil_froz_hydro)
    END IF
    CALL xios_orchidee_send_field("profil_froz_hydro_ns", profil_froz_hydro_ns)
    CALL xios_orchidee_send_field("kk_moy",kk_moy) ! in mm/d

    !! Calculate diagnostic variables on Landuse tiles for LUMIP/CMIP6
    humtot_lut(:,:)=0
    humtot_top_lut(:,:)=0
    mrro_lut(:,:)=0
    DO jv=1,nvm
       jst=pref_soil_veg(jv) ! soil tile index
       IF (natural(jv)) THEN
          humtot_lut(:,id_psl) = humtot_lut(:,id_psl) + tmc(:,jst)*veget_max(:,jv)
          humtot_top_lut(:,id_psl) = humtot_top_lut(:,id_psl) + tmc_top(:,jst)*veget_max(:,jv)
          mrro_lut(:,id_psl) = mrro_lut(:,id_psl) + (dr_ns(:,jst)+ru_ns(:,jst))*veget_max(:,jv)
       ELSE
          humtot_lut(:,id_crp) = humtot_lut(:,id_crp) + tmc(:,jst)*veget_max(:,jv)
          humtot_top_lut(:,id_crp) = humtot_top_lut(:,id_crp) + tmc_top(:,jst)*veget_max(:,jv)
          mrro_lut(:,id_crp) = mrro_lut(:,id_crp) + (dr_ns(:,jst)+ru_ns(:,jst))*veget_max(:,jv)
       ENDIF
    END DO

    WHERE (fraclut(:,id_psl)>min_sechiba)
       humtot_lut(:,id_psl) = humtot_lut(:,id_psl)/fraclut(:,id_psl)
       humtot_top_lut(:,id_psl) = humtot_top_lut(:,id_psl)/fraclut(:,id_psl)
       mrro_lut(:,id_psl) = mrro_lut(:,id_psl)/fraclut(:,id_psl)/dt_sechiba
    ELSEWHERE
       humtot_lut(:,id_psl) = val_exp
       humtot_top_lut(:,id_psl) = val_exp
       mrro_lut(:,id_psl) = val_exp
    END WHERE
    WHERE (fraclut(:,id_crp)>min_sechiba)
       humtot_lut(:,id_crp) = humtot_lut(:,id_crp)/fraclut(:,id_crp)
       humtot_top_lut(:,id_crp) = humtot_top_lut(:,id_crp)/fraclut(:,id_crp)
       mrro_lut(:,id_crp) = mrro_lut(:,id_crp)/fraclut(:,id_crp)/dt_sechiba
    ELSEWHERE
       humtot_lut(:,id_crp) = val_exp
       humtot_top_lut(:,id_crp) = val_exp
       mrro_lut(:,id_crp) = val_exp
    END WHERE

    humtot_lut(:,id_pst) = val_exp
    humtot_lut(:,id_urb) = val_exp
    humtot_top_lut(:,id_pst) = val_exp
    humtot_top_lut(:,id_urb) = val_exp
    mrro_lut(:,id_pst) = val_exp
    mrro_lut(:,id_urb) = val_exp

    CALL xios_orchidee_send_field("humtot_lut",humtot_lut)
    CALL xios_orchidee_send_field("humtot_top_lut",humtot_top_lut)
    CALL xios_orchidee_send_field("mrro_lut",mrro_lut)

    ! Write diagnistic for soil moisture nudging
    IF (ok_nudge_mc) CALL hydrol_nudge_mc_diag(kjpindex, soiltile)


    IF ( .NOT. almaoutput ) THEN
       CALL histwrite_p(hist_id, 'frac_bare', kjit, frac_bare, kjpindex*nvm, indexveg)

       CALL histwrite_p(hist_id, 'moistc', kjit,mc, kjpindex*nslm*nstm, indexlayer)
       CALL histwrite_p(hist_id, 'kfactroot', kjit, kfact_root, kjpindex*nslm*nstm, indexlayer)
       CALL histwrite_p(hist_id, 'vegetsoil', kjit,vegetmax_soil, kjpindex*nvm*nstm, indexveg)
       CALL histwrite_p(hist_id, 'evapnu_soil', kjit, ae_ns, kjpindex*nstm, indexsoil)
       CALL histwrite_p(hist_id, 'drainage_soil', kjit, dr_ns, kjpindex*nstm, indexsoil)
       CALL histwrite_p(hist_id, 'transpir_soil', kjit, tr_ns, kjpindex*nstm, indexsoil)
       CALL histwrite_p(hist_id, 'runoff_soil', kjit, ru_ns, kjpindex*nstm, indexsoil)
       CALL histwrite_p(hist_id, 'humtot_soil', kjit, tmc, kjpindex*nstm, indexsoil)
       ! mrso is a perfect duplicate of humtot
       CALL histwrite_p(hist_id, 'humtot', kjit, humtot, kjpindex, index)
       CALL histwrite_p(hist_id, 'mrso', kjit, humtot, kjpindex, index)
       CALL histwrite_p(hist_id, 'mrsos', kjit, humtot_top, kjpindex, index)
       njsc_tmp(:)=njsc(:)
       CALL histwrite_p(hist_id, 'soilindex', kjit, njsc_tmp, kjpindex, index)
       CALL histwrite_p(hist_id, 'humrel',   kjit, humrel,   kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'drainage', kjit, drainage, kjpindex, index)
       ! NB! According to histdef in intersurf, the variables 'runoff' and 'mrros' have different units
       CALL histwrite_p(hist_id, 'runoff', kjit, runoff, kjpindex, index)
       CALL histwrite_p(hist_id, 'mrros', kjit, runoff, kjpindex, index)
       histvar(:)=(runoff(:)+drainage(:))
       CALL histwrite_p(hist_id, 'mrro', kjit, histvar, kjpindex, index)
       CALL histwrite_p(hist_id, 'precisol', kjit, precisol, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'rain', kjit, precip_rain, kjpindex, index)

       histvar(:)=(precip_rain(:)-SUM(throughfall(:,:),dim=2))
       CALL histwrite_p(hist_id, 'prveg', kjit, histvar, kjpindex, index)

       CALL histwrite_p(hist_id, 'snowf', kjit, precip_snow, kjpindex, index)
       CALL histwrite_p(hist_id, 'qsintmax', kjit, qsintmax, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'qsintveg', kjit, qsintveg, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'snowmelt',kjit,snowmelt,kjpindex,index)
       CALL histwrite_p(hist_id, 'shumdiag_perma',kjit,shumdiag_perma,kjpindex*nslm,indexnslm)

       IF ( do_floodplains ) THEN
          CALL histwrite_p(hist_id, 'floodout', kjit, floodout, kjpindex, index)
       ENDIF
       !
       IF ( hist2_id > 0 ) THEN
          CALL histwrite_p(hist2_id, 'moistc', kjit,mc, kjpindex*nslm*nstm, indexlayer)
          CALL histwrite_p(hist2_id, 'kfactroot', kjit, kfact_root, kjpindex*nslm*nstm, indexlayer)
          CALL histwrite_p(hist2_id, 'vegetsoil', kjit,vegetmax_soil, kjpindex*nvm*nstm, indexveg)
          CALL histwrite_p(hist2_id, 'evapnu_soil', kjit, ae_ns, kjpindex*nstm, indexsoil)
          CALL histwrite_p(hist2_id, 'drainage_soil', kjit, dr_ns, kjpindex*nstm, indexsoil)
          CALL histwrite_p(hist2_id, 'transpir_soil', kjit, tr_ns, kjpindex*nstm, indexsoil)
          CALL histwrite_p(hist2_id, 'runoff_soil', kjit, ru_ns, kjpindex*nstm, indexsoil)
          CALL histwrite_p(hist2_id, 'humtot_soil', kjit, tmc, kjpindex*nstm, indexsoil)
          ! mrso is a perfect duplicate of humtot
          CALL histwrite_p(hist2_id, 'humtot', kjit, humtot, kjpindex, index)
          CALL histwrite_p(hist2_id, 'mrso', kjit, humtot, kjpindex, index)
          CALL histwrite_p(hist2_id, 'mrsos', kjit, humtot_top, kjpindex, index)
          njsc_tmp(:)=njsc(:)
          CALL histwrite_p(hist2_id, 'soilindex', kjit, njsc_tmp, kjpindex, index)
          CALL histwrite_p(hist2_id, 'humrel',   kjit, humrel,   kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'drainage', kjit, drainage, kjpindex, index)
          ! NB! According to histdef in intersurf, the variables 'runoff' and 'mrros' have different units
          CALL histwrite_p(hist2_id, 'runoff', kjit, runoff, kjpindex, index)
          CALL histwrite_p(hist2_id, 'mrros', kjit, runoff, kjpindex, index)
          histvar(:)=(runoff(:)+drainage(:))
          CALL histwrite_p(hist2_id, 'mrro', kjit, histvar, kjpindex, index)

          IF ( do_floodplains ) THEN
             CALL histwrite_p(hist2_id, 'floodout', kjit, floodout, kjpindex, index)
          ENDIF
          CALL histwrite_p(hist2_id, 'precisol', kjit, precisol, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'rain', kjit, precip_rain, kjpindex, index)
          CALL histwrite_p(hist2_id, 'snowf', kjit, precip_snow, kjpindex, index)
          CALL histwrite_p(hist2_id, 'snowmelt',kjit,snowmelt,kjpindex,index)
          CALL histwrite_p(hist2_id, 'qsintmax', kjit, qsintmax, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'qsintveg', kjit, qsintveg, kjpindex*nvm, indexveg)
       ENDIF
    ELSE
       CALL histwrite_p(hist_id, 'Snowf', kjit, precip_snow, kjpindex, index)
       CALL histwrite_p(hist_id, 'Rainf', kjit, precip_rain, kjpindex, index)
       CALL histwrite_p(hist_id, 'Qs', kjit, runoff, kjpindex, index)
       CALL histwrite_p(hist_id, 'Qsb', kjit, drainage, kjpindex, index)
       CALL histwrite_p(hist_id, 'Qsm', kjit, snowmelt, kjpindex, index)
       CALL histwrite_p(hist_id, 'DelSoilMoist', kjit, delsoilmoist, kjpindex, index)
       CALL histwrite_p(hist_id, 'DelSWE', kjit, delswe, kjpindex, index)
       CALL histwrite_p(hist_id, 'DelIntercept', kjit, delintercept, kjpindex, index)
       !
       CALL histwrite_p(hist_id, 'SoilMoist', kjit, soilmoist, kjpindex*nslm, indexlayer)
       CALL histwrite_p(hist_id, 'SoilWet', kjit, soilwet, kjpindex, index)
       !
       CALL histwrite_p(hist_id, 'RootMoist', kjit, tot_watsoil_end, kjpindex, index)
       CALL histwrite_p(hist_id, 'SubSnow', kjit, vevapsno, kjpindex, index)

       IF ( hist2_id > 0 ) THEN
          CALL histwrite_p(hist2_id, 'Snowf', kjit, precip_snow, kjpindex, index)
          CALL histwrite_p(hist2_id, 'Rainf', kjit, precip_rain, kjpindex, index)
          CALL histwrite_p(hist2_id, 'Qs', kjit, runoff, kjpindex, index)
          CALL histwrite_p(hist2_id, 'Qsb', kjit, drainage, kjpindex, index)
          CALL histwrite_p(hist2_id, 'Qsm', kjit, snowmelt, kjpindex, index)
          CALL histwrite_p(hist2_id, 'DelSoilMoist', kjit, delsoilmoist, kjpindex, index)
          CALL histwrite_p(hist2_id, 'DelSWE', kjit, delswe, kjpindex, index)
          CALL histwrite_p(hist2_id, 'DelIntercept', kjit, delintercept, kjpindex, index)
          !
          CALL histwrite_p(hist2_id, 'SoilMoist', kjit, soilmoist, kjpindex*nslm, indexlayer)
          CALL histwrite_p(hist2_id, 'SoilWet', kjit, soilwet, kjpindex, index)
          !
          CALL histwrite_p(hist2_id, 'RootMoist', kjit, tot_watsoil_end, kjpindex, index)
          CALL histwrite_p(hist2_id, 'SubSnow', kjit, vevapsno, kjpindex, index)
       ENDIF
    ENDIF

    IF (ok_freeze_cwrr) THEN
       CALL histwrite_p(hist_id, 'profil_froz_hydro', kjit,profil_froz_hydro , kjpindex*nslm, indexlayer)
    ENDIF
    !CALL histwrite_p(hist_id, 'kk_moy', kjit, kk_moy,kjpindex*nslm, indexlayer) ! averaged over soiltiles
    CALL histwrite_p(hist_id, 'profil_froz_hydro', kjit, profil_froz_hydro_ns, kjpindex*nslm*nstm, indexlayer)

    ! Copy soilmoist into a local variable to be sent to thermosoil
    soilmoist_out(:,:) = soilmoist(:,:)
    soilmoist_out_s(:,:,:) = soilmoist_s(:,:,:)

    ! Copy mcs and mcfc into local variables to be sent to stomate_soilcarbon.
    ! note: these two variabels are usesful in nitrogen_dynamics subroutine of stomate_soilcarbon
    ! but they were defined initially without considering peatland.
    ! peat is now considered also.
    mcs_hydrol(:) =  mcs(:)
    mcfc_hydrol(:) =  mcfc(:)

    IF (ok_peat_hydro) THEN
       DO ji=1, kjpindex
          mcs_hydrol(ji) = mcs(ji)*(un-soiltile(ji,4))+mcs_peat*soiltile(ji,4)
          mcfc_hydrol(ji) = mcfc(ji)*(un-soiltile(ji,4))+mcf_peat*soiltile(ji,4)
       ENDDO
     
    ELSE IF (agri_peat) THEN
        DO ji=1, kjpindex
           mcs_hydrol(ji) = mcs(ji)*(un-soiltile(ji,4)-soiltile(ji,5))+mcs_peat*(soiltile(ji,4)+ soiltile(ji,5))
           mcfc_hydrol(ji) = mcfc(ji)*(un-soiltile(ji,4)-soiltile(ji,5))+mcf_peat*(soiltile(ji,4)+ soiltile(ji,5))
        ENDDO
        
    ELSE IF (agri_peat .AND. tides) THEN
       DO ji=1, kjpindex
          mcs_hydrol(ji) = mcs(ji)*(un-soiltile(ji,4)-soiltile(ji,5)-soiltile(ji,5))+mcs_peat*(soiltile(ji,4)+ soiltile(ji,5)+ soiltile(ji,5))
          mcfc_hydrol(ji) = mcfc(ji)*(un-soiltile(ji,4)-soiltile(ji,5)-soiltile(ji,5))+mcf_peat*(soiltile(ji,4)+ soiltile(ji,5)+ soiltile(ji,5))
       ENDDO
    ENDIF
    !     
    IF ((printlev>=3).OR.(printlev_loc>=4)) THEN
       WRITE (numout,*) ' hydrol_main Done '
       CALL flush(numout)
    ENDIF

  END SUBROUTINE hydrol_main


!! ================================================================================================================================
!! SUBROUTINE 	: hydrol_finalize
!!
!>\BRIEF         
!!
!! DESCRIPTION : This subroutine writes the module variables and variables calculated in hydrol to restart file
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE hydrol_finalize( kjit,           kjpindex,   rest_id,  vegstress,  &
                              qsintveg,       humrel,     snow,     snow_age,   &
                              snow_nobio_age, snowrho,    snowtemp, snowdz,     &
                              snowheat,       snowgrain,                        &
                              drysoil_frac, evap_bare_lim, evap_bare_lim_ns,    &
                              mc_out, ksoil, root_profile, us,                  &
                              icetemp, psi_leaf, psi_leaf_next, psi_leaf_midday,&
                              psi_sto_leaf_save,          psi_sto_wood_save,    &
                              psi_root_sup,   psi_root_inf,        psi_soil_sup,&
                              psi_soil_inf,   psi_xylem_trunk,   psi_xylem_leaf,&
                              psi_xylem_collar,           psi_sto_wood,         &
                              psi_sto_leaf,   mc_i_sup,   mc_i_inf, F_absorption,&
                              wtp,liqwt_ratio)

    !! 0. Variable and parameter declaration
    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                                         :: kjit           !! Time step number 
    INTEGER(i_std), INTENT(in)                                         :: kjpindex       !! Domain size
    INTEGER(i_std),INTENT (in)                                         :: rest_id        !! Restart file identifier
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)                  :: vegstress      !! Veg. moisture stress (only for vegetation growth)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)                  :: qsintveg       !! Water on vegetation due to interception
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)                  :: humrel
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)                      :: snow           !! Snow mass [Kg/m^2]
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)                      :: snow_age       !! Snow age
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (in)               :: snow_nobio_age !! Snow age on ice, lakes, ...
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(in)                :: snowrho        !! Snow density (Kg/m^3)
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(in)                :: snowtemp       !! Snow temperature (K)
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(in)                :: snowdz         !! Snow layer thickness [m]
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(in)                :: snowheat       !! Snow heat content
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(in)                :: snowgrain      !! Snow grainsize
    REAL(r_std),DIMENSION (kjpindex),INTENT(in)                        :: drysoil_frac   !! function of litter wetness
    REAL(r_std),DIMENSION (kjpindex),INTENT(in)                        :: evap_bare_lim
    REAL(r_std),DIMENSION (kjpindex,nstm),INTENT(in)                   :: evap_bare_lim_ns
    REAL(r_std),DIMENSION (kjpindex,nslm,nstm), INTENT (in)            :: mc_out
    REAL(r_std),DIMENSION (kjpindex,nslm,nstm), INTENT (in)            :: ksoil
    REAL(r_std),DIMENSION (kjpindex,nvm,nslm,nroot_prof), INTENT(in)   :: root_profile    !! Normalized root mass/length fraction in each soil layer
                                                                                          !! (0-1, unitless)
    REAL(r_std),DIMENSION (kjpindex,nvm,nstm,nslm), INTENT(in)         :: us              !! Water stress index for transpiration
                                                                                          !! (by soil layer and PFT) (0-1, unitless)
    REAL(r_std), DIMENSION (kjpindex,nice), INTENT(in)                 :: icetemp         !! Ice temperature (K)    
    REAL(r_std), DIMENSION (kjpindex,nvm,ncirc),      INTENT(in)       :: psi_leaf               !! Leaf Water Potential (MPa)
    REAL(r_std), DIMENSION (kjpindex,nvm),            INTENT(in)       :: psi_leaf_next          !! Approximated Leaf water potential at time step n+1 (MPa) (=psi_leaf when no stress)
    REAL(r_std), DIMENSION (kjpindex,nvm,ncirc),      INTENT(in)       :: psi_leaf_midday        !! Leaf Water Potential at midday (MPa)
    REAL(r_std), DIMENSION (kjpindex,nvm,ncirc),      INTENT(in)       :: psi_sto_leaf_save      !! Leaf storage water potential at time step n-1 (MPa)
    REAL(r_std), DIMENSION (kjpindex,nvm,ncirc),      INTENT(in)       :: psi_sto_wood_save      !! Wood storage water potential at time step n-1 (MPa)
    REAL(r_std), DIMENSION (kjpindex,nvm,ncirc),      INTENT(in)       :: psi_root_sup           !! Superficial root water potential (MPa)
    REAL(r_std), DIMENSION (kjpindex,nvm,ncirc),      INTENT(in)       :: psi_root_inf           !! Inferior root water potential (MPa)
    REAL(r_std), DIMENSION (kjpindex,nstm),     INTENT(in)             :: psi_soil_sup           !! Superficial soil water potential (MPa)
    REAL(r_std), DIMENSION (kjpindex,nstm),     INTENT(in)             :: psi_soil_inf           !! Inferior soil water potential (MPa)
    REAL(r_std), DIMENSION (kjpindex,nvm,ncirc),      INTENT(in)       :: psi_xylem_trunk        !! Xylem (trunk level) water potential (MPa)
    REAL(r_std), DIMENSION (kjpindex,nvm,ncirc),      INTENT(in)       :: psi_xylem_leaf         !! Xylem (leaf level) water potential (MPa)
    REAL(r_std), DIMENSION (kjpindex,nvm,ncirc),      INTENT(in)       :: psi_xylem_collar       !! Xylem (collar level) water potential (MPa)
    REAL(r_std), DIMENSION (kjpindex,nvm,ncirc),      INTENT(in)       :: psi_sto_wood           !! Wood storage water potential (MPa)
    REAL(r_std), DIMENSION (kjpindex,nvm,ncirc),      INTENT(in)       :: psi_sto_leaf           !! Leaf storage water potential (MPa)
    REAL(r_std), DIMENSION (kjpindex,nvm,ncirc,nrp),  INTENT(in)       :: mc_i_sup               !! Water content at each node of the absorption muff in the superficial soil layer (m^3/m^3)
    REAL(r_std), DIMENSION (kjpindex,nvm,ncirc,nrp),  INTENT(in)       :: mc_i_inf               !! Water content at each node of the absorption muff in the inferior soil layer (m^3/m^3)  
    REAL(r_std), DIMENSION (kjpindex,nvm),      INTENT(in)             :: F_absorption           !! Total root absorption flux (m^3/s)
    
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)                     :: wtp                    !! water table position
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)                     :: liqwt_ratio            !! liquid water ratio in peat, compared to the satured one
    
    !! 0.4 Local variables
    INTEGER(i_std)                                                     :: jst, jsl
   
!_ ================================================================================================================================


    IF (printlev>=3) WRITE (numout,*) 'Write restart file with HYDROLOGIC variables '

    CALL restput_p(rest_id, 'moistc', nbp_glo,  nslm, nstm, kjit, mc, 'scatter',  nbp_glo, index_g)
    !-
    CALL restput_p(rest_id, 'moistcl', nbp_glo,  nslm, nstm, kjit, mcl, 'scatter',  nbp_glo, index_g)
    !-
    IF (ok_nudge_mc) THEN
       CALL restput_p(rest_id, 'mc_read_next', nbp_glo,  nslm, nstm, kjit, mc_read_next, 'scatter',  nbp_glo, index_g)
    END IF

    IF (ok_nudge_snow) THEN
       CALL restput_p(rest_id, 'snowdz_read_next', nbp_glo,  nsnow, 1, kjit, snowdz_read_next(:,:), &
            'scatter',  nbp_glo, index_g)
       !-
       CALL restput_p(rest_id, 'snowrho_read_next', nbp_glo,  nsnow, 1, kjit, snowrho_read_next(:,:), &
            'scatter',  nbp_glo, index_g)
       !-
       CALL restput_p(rest_id, 'snowtemp_read_next', nbp_glo,  nsnow, 1, kjit, snowtemp_read_next(:,:), &
            'scatter',  nbp_glo, index_g)
    END IF

    CALL restput_p(rest_id, 'us', nbp_glo, nvm, nstm, nslm, kjit, us,'scatter',nbp_glo,index_g)
    !-
    CALL restput_p(rest_id, 'free_drain_coef', nbp_glo,   nstm, 1, kjit,  free_drain_coef, 'scatter',  nbp_glo, index_g)
    !-
    CALL restput_p(rest_id, 'zwt_force', nbp_glo,   nstm, 1, kjit,  zwt_force, 'scatter',  nbp_glo, index_g)
    !-
    CALL restput_p(rest_id, 'water2infilt', nbp_glo,   nstm, 1, kjit,  water2infilt, 'scatter',  nbp_glo, index_g)
    !-
    CALL restput_p(rest_id, 'ae_ns', nbp_glo,   nstm, 1, kjit,  ae_ns, 'scatter',  nbp_glo, index_g)
    !-
    CALL restput_p(rest_id, 'vegstress', nbp_glo,   nvm, 1, kjit,  vegstress, 'scatter',  nbp_glo, index_g)
    !-
    CALL restput_p(rest_id, 'snow', nbp_glo,   1, 1, kjit,  snow, 'scatter',  nbp_glo, index_g)
    !-
    CALL restput_p(rest_id, 'snow_age', nbp_glo,   1, 1, kjit,  snow_age, 'scatter',  nbp_glo, index_g)
    !-
    CALL restput_p(rest_id, 'snow_nobio_age', nbp_glo,   nnobio, 1, kjit,  snow_nobio_age, 'scatter', nbp_glo, index_g)
    !-
    CALL restput_p(rest_id, 'qsintveg', nbp_glo, nvm, 1, kjit,  qsintveg, 'scatter',  nbp_glo, index_g)
    !-
    CALL restput_p(rest_id, 'evap_bare_lim_ns', nbp_glo, nstm, 1, kjit,  evap_bare_lim_ns, 'scatter',  nbp_glo, index_g)
    !-
    CALL restput_p(rest_id, 'evap_bare_lim', nbp_glo, 1, 1, kjit,  evap_bare_lim, 'scatter',  nbp_glo, index_g)
    !-
    CALL restput_p(rest_id, 'root_profile_struc', nbp_glo, nvm, nslm, kjit,  root_profile(:,:,:,istruc), 'scatter',  nbp_glo, index_g)
    !-
    CALL restput_p(rest_id, 'root_profile_func', nbp_glo, nvm, nslm, kjit,  root_profile(:,:,:,ifunc), 'scatter',  nbp_glo, index_g)
    !-
    CALL restput_p(rest_id, 'resdist', nbp_glo, nstm, 1, kjit,  resdist, 'scatter',  nbp_glo, index_g)
    !-
    CALL restput_p(rest_id, 'vegtot_old', nbp_glo, 1, 1, kjit,  vegtot_old, 'scatter',  nbp_glo, index_g)            
    !-
    CALL restput_p(rest_id, 'drysoil_frac', nbp_glo,   1, 1, kjit, drysoil_frac, 'scatter', nbp_glo, index_g)
    !-
    CALL restput_p(rest_id, 'humrel', nbp_glo,   nvm, 1, kjit,  humrel, 'scatter',  nbp_glo, index_g)
    !-
    CALL restput_p(rest_id, 'tot_watveg_beg', nbp_glo,  1, 1, kjit,  tot_watveg_beg, 'scatter',  nbp_glo, index_g)
    !-
    CALL restput_p(rest_id, 'tot_watsoil_beg', nbp_glo, 1, 1, kjit,  tot_watsoil_beg, 'scatter',  nbp_glo, index_g)
    !-
    CALL restput_p(rest_id, 'snow_beg', nbp_glo,        1, 1, kjit,  snow_beg, 'scatter',  nbp_glo, index_g)
    !-
    IF(ok_hydrol_arch)THEN
       CALL restput_p(rest_id, 'mc_out', nbp_glo, nslm,  nstm, kjit, mc_out, 'scatter',  nbp_glo, index_g)
       CALL restput_p(rest_id, 'ksoil', nbp_glo, nslm,  nstm, kjit, ksoil, 'scatter',  nbp_glo, index_g) 
       CALL restput_p(rest_id, 'psi_leaf', nbp_glo,   nvm, ncirc, kjit,  psi_leaf, 'scatter',  nbp_glo, index_g)
       CALL restput_p(rest_id, 'psi_leaf_next', nbp_glo,   nvm, 1, kjit,  psi_leaf_next, 'scatter',  nbp_glo, index_g)
       CALL restput_p(rest_id, 'psi_leaf_midday', nbp_glo, nvm, ncirc, kjit, psi_leaf_midday, 'scatter',  nbp_glo, index_g)
       CALL restput_p(rest_id, 'psi_sto_leaf_save', nbp_glo,   nvm, ncirc, kjit,  psi_sto_leaf_save, 'scatter',  nbp_glo, index_g)
       CALL restput_p(rest_id, 'psi_sto_wood_save', nbp_glo,   nvm, ncirc, kjit,  psi_sto_wood_save, 'scatter',  nbp_glo, index_g)
       CALL restput_p(rest_id, 'psi_root_sup', nbp_glo,   nvm, ncirc, kjit,  psi_root_sup, 'scatter',  nbp_glo, index_g)
       CALL restput_p(rest_id, 'psi_root_inf', nbp_glo,   nvm, ncirc, kjit,  psi_root_inf, 'scatter',  nbp_glo, index_g)
       CALL restput_p(rest_id, 'psi_soil_sup', nbp_glo,   nstm, 1, kjit, psi_soil_sup, 'scatter',  nbp_glo, index_g)
       CALL restput_p(rest_id, 'psi_soil_inf', nbp_glo,   nstm, 1, kjit, psi_soil_inf, 'scatter',  nbp_glo, index_g) 
       CALL restput_p(rest_id, 'psi_xylem_trunk', nbp_glo,   nvm, ncirc, kjit,  psi_xylem_trunk, 'scatter',  nbp_glo, index_g)
       CALL restput_p(rest_id, 'psi_xylem_leaf', nbp_glo,   nvm, ncirc, kjit,  psi_xylem_leaf, 'scatter',  nbp_glo, index_g)
       CALL restput_p(rest_id, 'psi_xylem_collar', nbp_glo,   nvm, ncirc, kjit,  psi_xylem_collar, 'scatter',  nbp_glo, index_g)
       CALL restput_p(rest_id, 'psi_sto_wood', nbp_glo,   nvm, ncirc, kjit,  psi_sto_wood, 'scatter',  nbp_glo, index_g)
       CALL restput_p(rest_id, 'psi_sto_leaf', nbp_glo,   nvm, ncirc, kjit,  psi_sto_leaf, 'scatter',  nbp_glo, index_g)
       CALL restput_p(rest_id, 'F_absorption', nbp_glo,   nvm, 1, kjit,  F_absorption, 'scatter',  nbp_glo, index_g)
       CALL restput_p(rest_id, 'mc_i_sup', nbp_glo,   nvm, ncirc, nrp, kjit,  mc_i_sup, 'scatter',  nbp_glo, index_g)
       CALL restput_p(rest_id, 'mc_i_inf', nbp_glo,   nvm, ncirc, nrp, kjit,  mc_i_inf, 'scatter',  nbp_glo, index_g)
    END IF
    
    ! Write variables for explictsnow module to restart file
    CALL explicitsnow_finalize ( kjit,     kjpindex, rest_id,    snowrho,   &
         snowtemp, snowdz,   snowheat,   snowgrain, icetemp)

    ! peatland and tides
    IF ( ok_peat_hydro ) THEN
       var_name= 'run2peat'
       CALL restput_p(rest_id, var_name, nbp_glo,   1, 1, kjit,  run2peat, 'scatter', nbp_glo, index_g)
       
       var_name= 'wt_ab'
       CALL restput_p(rest_id, var_name, nbp_glo,   1, 1, kjit, wt_ab,'scatter', nbp_glo, index_g)
       
       var_name= 'wtp'
       CALL restput_p(rest_id, var_name, nbp_glo,   1, 1, kjit,  wtp, 'scatter', nbp_glo, index_g)
       
       var_name= 'liqwt_ratio'
       CALL restput_p(rest_id, var_name, nbp_glo,   1, 1, kjit,  liqwt_ratio, 'scatter', nbp_glo, index_g)
    ENDIF

    IF (tides) THEN
       var_name= 'wt_ab_tide'
       CALL restput_p(rest_id, var_name, nbp_glo,   1, 1, kjit, wt_ab_tide,'scatter', nbp_glo, index_g)
       
       var_name= 'run2man'
       CALL restput_p(rest_id, var_name, nbp_glo,   1, 1, kjit,  run2man,'scatter', nbp_glo, index_g)
    ENDIF

  END SUBROUTINE hydrol_finalize


!! ================================================================================================================================
!! SUBROUTINE   : hydrol_init
!!
!>\BRIEF        Initializations and memory allocation   
!!
!! DESCRIPTION  :
!! - 1 Some initializations
!! - 2 make dynamic allocation with good dimension
!! - 2.1 array allocation for soil textur
!! - 2.2 Soil texture choice
!! - 3 Other array allocation
!! - 4 Open restart input file and read data for HYDROLOGIC process
!! - 5 get restart values if none were found in the restart file
!! - 6 Vegetation array      
!! - 7 set humrelv from us
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!!_ hydrol_init

  SUBROUTINE hydrol_init(ks, nvan, avan, mcr, mcs, mcfc, mcw, &
       kjit, kjpindex, index, rest_id, veget_max, frac_nobio, soiltile, &
       humrel,  vegstress, snow,       snow_age,   snow_nobio_age, qsintveg, &
       snowdz,  snowgrain, snowrho,    snowtemp,   snowheat, &
       drysoil_frac,  evap_bare_lim,   evap_bare_lim_ns,     &
       mc_out,        ksoil,           root_profile,      us, &
       icetemp,       icedz,           ice_sheet_mask,  &
       psi_leaf,      psi_leaf_next,   psi_leaf_midday,   psi_sto_leaf_save, psi_sto_wood_save, &
       psi_root_sup,  psi_root_inf,    psi_soil_sup, psi_soil_inf, psi_xylem_trunk,   psi_xylem_leaf, &
       psi_xylem_collar, psi_sto_wood, psi_sto_leaf, &
       mc_i_sup,         mc_i_inf,     F_absorption, &
       wtp,liqwt_ratio)
    

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT (in)                         :: kjit               !! Time step number 
    INTEGER(i_std), INTENT (in)                         :: kjpindex           !! Domain size
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)    :: index              !! Indeces of the points on the map
    INTEGER(i_std), INTENT (in)                         :: rest_id            !! _Restart_ file identifier 
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)   :: veget_max          !! Carte de vegetation max
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (in):: frac_nobio         !! Fraction of ice, lakes, ...
    REAL(r_std),DIMENSION (kjpindex,nstm), INTENT (in)  :: soiltile           !! Fraction of each soil tile within vegtot (0-1, unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: ks                 !! Hydraulic conductivity at saturation (mm {-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: nvan               !! Van Genuchten coeficients n (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: avan               !! Van Genuchten coeficients a (mm-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: mcr                !! Residual volumetric water content (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: mcs                !! Saturated volumetric water content (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: mcfc               !! Volumetric water content at field capacity (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: mcw                !! Volumetric water content at wilting point (m^{3} m^{-3})

    !! 0.2 Output variables

    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)                  :: humrel           !! Stress hydrique, relative humidity
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)                  :: vegstress        !! Veg. moisture stress (only for vegetation growth)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)                      :: snow             !! Snow mass [Kg/m^2]
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)                      :: snow_age         !! Snow age
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (out)               :: snow_nobio_age   !! Snow age on ice, lakes, ...
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)                  :: qsintveg         !! Water on vegetation due to interception
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT(out)                  :: snowdz           !! Snow depth [m]
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT(out)                  :: snowgrain        !! Snow grain size  
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT(out)                  :: snowheat         !! Snow heat content
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT(out)                  :: snowtemp         !! Snow temperature (K)
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT(out)                  :: snowrho          !! Snow density (Kg/m^3)
    REAL(r_std),DIMENSION (kjpindex),INTENT(out)                        :: drysoil_frac     !! function of litter wetness
    REAL(r_std),DIMENSION (kjpindex),INTENT(out)                        :: evap_bare_lim
    REAL(r_std),DIMENSION (kjpindex,nstm),INTENT(out)                   :: evap_bare_lim_ns
    REAL(r_std),DIMENSION (kjpindex,nslm,nstm), INTENT (out)            :: mc_out
    REAL(r_std),DIMENSION (kjpindex,nslm,nstm), INTENT (out)            :: ksoil
    REAL(r_std),DIMENSION (kjpindex,nvm,nslm,nroot_prof), INTENT(out)   :: root_profile       !! Normalized root mass/length fraction in each soil layer 
                                                                                              !! (0-1, unitless)
    REAL(r_std),DIMENSION (kjpindex,nvm,nstm,nslm), INTENT (out)        :: us                 !! Water stress index for transpiration
                                                                                              !! (by soil layer and PFT) (0-1, unitless)
    REAL(r_std), DIMENSION (kjpindex,nvm,ncirc),      INTENT(out)       :: psi_leaf           !! Leaf water potential (MPa)
    REAL(r_std), DIMENSION (kjpindex,nvm),      INTENT(out)             :: psi_leaf_next      !! Approximated Leaf water potential at time step n+1 (MPa) (=psi_leaf when no stress)
    REAL(r_std), DIMENSION (kjpindex,nvm,ncirc),      INTENT(out)       :: psi_leaf_midday    !! Leaf water potential at midday (MPa)
    REAL(r_std), DIMENSION (kjpindex,nvm,ncirc),      INTENT(out)       :: psi_sto_leaf_save  !! Leaf storage water potential at time step n-1 (MPa)
    REAL(r_std), DIMENSION (kjpindex,nvm,ncirc),      INTENT(out)       :: psi_sto_wood_save  !! Wood storage water potential at time step n-1 (MPa)
    REAL(r_std), DIMENSION (kjpindex,nvm,ncirc),      INTENT(out)       :: psi_root_sup       !! Superficial root water potential (MPa)
    REAL(r_std), DIMENSION (kjpindex,nvm,ncirc),      INTENT(out)       :: psi_root_inf       !! Inferior root water potential (MPa)
    REAL(r_std), DIMENSION (kjpindex,nstm),     INTENT(out)             :: psi_soil_sup       !! Superficial soil water potential (MPa)
    REAL(r_std), DIMENSION (kjpindex,nstm),     INTENT(out)             :: psi_soil_inf       !! Inferior soil water potential (MPa)
    REAL(r_std), DIMENSION (kjpindex,nvm,ncirc),      INTENT(out)       :: psi_xylem_trunk    !! Xylem (trunk level) water potential (MPa)
    REAL(r_std), DIMENSION (kjpindex,nvm,ncirc),      INTENT(out)       :: psi_xylem_leaf     !! Xylem (leaf level) water potential (MPa)
    REAL(r_std), DIMENSION (kjpindex,nvm,ncirc),      INTENT(out)       :: psi_xylem_collar   !! Xylem (collar level) water potential (MPa)
    REAL(r_std), DIMENSION (kjpindex,nvm,ncirc),      INTENT(out)       :: psi_sto_wood       !! Wood storage water potential (MPa)
    REAL(r_std), DIMENSION (kjpindex,nvm,ncirc),      INTENT(out)       :: psi_sto_leaf       !! Leaf storage water potential (MPa)
    REAL(r_std), DIMENSION (kjpindex,nvm,ncirc,nrp),  INTENT(out)       :: mc_i_sup           !! Water content at each node of the absorption muff in the superficial soil layer (m^3/m^3)
    REAL(r_std), DIMENSION (kjpindex,nvm,ncirc,nrp),  INTENT(out)       :: mc_i_inf           !! Water content at each node of the absorption muff in the inferior soil layer (m^3/m^3)
    REAL(r_std), DIMENSION (kjpindex,nvm),      INTENT(out)             :: F_absorption       !! Total root absorption flux (m^3/s)

    REAL(r_std),DIMENSION (kjpindex,nice),INTENT(out)                   :: icetemp          !! Ice temperature
    REAL(r_std), DIMENSION (kjpindex,nice),INTENT(out)                  :: icedz            !! Ice layer thickness [m]
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (out)                   :: ice_sheet_mask   !! Ice sheet mask
    
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)                       :: wtp              !! water table position
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)                       :: liqwt_ratio      !! Liquid water ratio in the peat, compared the the satured one


    !! 0.4 Local variables

    INTEGER(i_std)                                     :: ier                   !! Error code
    INTEGER(i_std)                                     :: ji                    !! Index of land grid cells (1)
    INTEGER(i_std)                                     :: jv                    !! Index of PFTs (1)
    INTEGER(i_std)                                     :: jst                   !! Index of soil tiles (1)
    INTEGER(i_std)                                     :: jsl                   !! Index of soil layers (1)
    INTEGER(i_std)                                     :: jsc                   !! Index of soil texture (1)
    INTEGER(i_std), PARAMETER                          :: error_level = 3       !! Error level for consistency check
                                                                                !! Switch to 2 tu turn fatal errors into warnings  
    REAL(r_std), ALLOCATABLE, DIMENSION (:)            :: free_drain_max        !! Temporary var for initialization of free_drain_coef 
    REAL(r_std), ALLOCATABLE, DIMENSION (:)            :: zwt_default           !! Temporary variable for initialization of zwt_force
    LOGICAL                                            :: zforce                !! To test if we force the WT in any of the soiltiles

!_ ================================================================================================================================

    !! 1 Some initializations
    !
    !Config Key   = DO_PONDS
    !Config Desc  = Should we include ponds 
    !Config Def   = n
    !Config If    = 
    !Config Help  = This parameters allows the user to ask the model
    !Config         to take into account the ponds and return 
    !Config         the water into the soil moisture. If this is 
    !Config         activated, then there is no reinfiltration 
    !Config         computed inside the hydrol module.
    !Config Units = [FLAG]
    !
    doponds = .FALSE.
    CALL getin_p('DO_PONDS', doponds)

    !Config Key   = FROZ_FRAC_CORR 
    !Config Desc  = Coefficient for the frozen fraction correction
    !Config Def   = 1.0
    !Config If    = OK_FREEZE
    !Config Help  =
    !Config Units = [-]
    froz_frac_corr = 1.0
    CALL getin_p("FROZ_FRAC_CORR", froz_frac_corr)

    !Config Key   = MAX_FROZ_HYDRO
    !Config Desc  = Coefficient for the frozen fraction correction
    !Config Def   = 1.0
    !Config If    = OK_FREEZE
    !Config Help  =
    !Config Units = [-]
    max_froz_hydro = 1.0
    CALL getin_p("MAX_FROZ_HYDRO", max_froz_hydro)

    !Config Key   = SMTOT_CORR
    !Config Desc  = Coefficient for the frozen fraction correction
    !Config Def   = 2.0
    !Config If    = OK_FREEZE
    !Config Help  =
    !Config Units = [-]
    smtot_corr = 2.0
    CALL getin_p("SMTOT_CORR", smtot_corr)

       
    !Config Key   = DO_RSOIL
    !Config Desc  = Should we reduce soil evaporation with a soil resistance
    !Config Def   = y if RSOIL_SCALE>0
    !Config If    = 
    !Config Help  = This parameters allows the user to ask the model
    !Config         to calculate a soil resistance to reduce the soil evaporation
    !Config Units = [FLAG]

    ! Note that DO_RSOIL could be removed but is kept to keep retro-compatibility with older orchidee.def parameter files
    IF (ALL(rsoil_scale(:)==zero)) THEN
       do_rsoil = .FALSE.
    ELSE
       do_rsoil = .TRUE.
    END IF
    CALL getin_p('DO_RSOIL', do_rsoil) 


    IF (ALL(rsoil_scale(:)==zero) .AND. do_rsoil) THEN
       ! Incoherence between rsoil_scale and do_rsoil
       CALL ipslerr_p(3,'hydrol_init','Incoherence between RSOIL_SCALE and DO_RSOIL',&
            'Remove DO_RSOIL from run.def or orchidee.def if using RSOIL_SCALE','')
    END IF

    !!! peatland
    ! 
    ! peatland-related definitions useful for hydrol
    !
    IF(ok_peat_hydro .OR. tides) THEN       
       !Config Key   = OK_RUN2PEAT
       !Config Desc  = Flag to allow adding the runoff of other soil tiles to the peat soil tile  
       !Config Def   = n
       !Config If    = OK_PEAT_HYDRO or PEAT_TIDES
       !Config Help  = 
       !Config Units = [FLAG]
       ok_run2peat= .TRUE.
       CALL getin_p('OK_RUN2PEAT', ok_run2peat)
       !
       !Config Key   = OK_WT_AB
       !Config Desc  = Flag to allow water table above the surface 
       !Config Def   = n
       !Config If    = OK_PEAT_HYDRO or PEAT_TIDES
       !Config Help  = 
       !Config Units = [FLAG]
       ok_wt_ab= .TRUE.
       CALL getin_p('OK_WT_AB', ok_wt_ab)
       !
       !Config Key   = MAX_WT_AB
       !Config Desc  = Flag to allow water table being above the surface 
       !Config Def   = 100.0
       !Config If    = OK_PEAT_HYDRO or PEAT_TIDES
       !Config Help  = 
       !Config Units = [mm]
       max_wt_ab= 100.
       CALL getin_p('MAX_WT_AB', max_wt_ab)
    ELSE
       ok_run2peat = .FALSE.
       ok_wt_ab = .FALSE.
       max_wt_ab = 100.
    ENDIF
 
    
    
    !Config Key   = AGRI_DRAIN
    !Config Desc  = Should we consider drainage over agriculture and peat 
    !Config Def   = n
    !Config If    = OK_PEAT_HYDRO AND AGRI_PEAT
    !Config Help  = 
    !Config Units = [FLAG]
    agri_drain= .FALSE. 
    CALL getin_p('AGRI_DRAIN',agri_drain)
    !
    !Config Key   = DRAIN_FACTOR 
    !Config Desc  = Coefficient for drainage over agriculture and peat
    !Config Def   = 1.0
    !Config If    = AGRI_DRAIN
    !Config Help  =
    !Config Units = [-]
    drain_factor = 1.0
    CALL getin_p('DRAIN_FACTOR',drain_factor)

    IF (ok_peat_hydro .AND. ok_thermodynamical_freezing) THEN
       ! Write a warning for this case which needs future improvements
       CALL ipslerr_p(2,'hydrol_init','ok_thermodynamical_freezing is True and peatland exits:',&
            'hydrol_soil_froz is not really adapted for peatland computation',&
            'because mcs_peat and mcr_peat does not represent the soil satuation in the whole grid.')
    END IF
    
    !! 2 make dynamic allocation with good dimension

    !! 2.1 array allocation for soil texture


    !! Campbell parametrisation sometimes used in the muff (Parameters for the new hydraulic architecture)

    ALLOCATE (b_muff(nscm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable b_muff','','')
    ALLOCATE (psi_air_entry(nscm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable psi_air_entry','','')
    ALLOCATE (mcr_sup(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable mcr_sup','','')
    ALLOCATE (mcr_inf(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable mcr_inf','','')
    ALLOCATE (mcs_sup(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable mcs_sup','','')
    ALLOCATE (mcs_inf(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable mcs_inf','','')
     

    ALLOCATE (pcent(nscm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable pcent','','')

    ALLOCATE (mc_awet(nscm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable mc_awet','','')

    ALLOCATE (mc_adry(nscm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable mc_adry','','')

    ALLOCATE (rsoil_scale_tile(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable rsoil_scale_tile','','')

    ! peat
    ALLOCATE (kfact_peat(nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable kfact_peat','','') 
    ALLOCATE (mc_lin_peat(imin:imax),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable mc_lin_peat','','')
    ALLOCATE (k_lin_peat(imin:imax, nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable k_lin_peat','','')
    ALLOCATE (d_lin_peat(imin:imax, nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable d_lin_peat','','')
    ALLOCATE (a_lin_peat(imin:imax, nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable a_lin_peat','','')
    ALLOCATE (b_lin_peat(imin:imax, nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable b_lin_peat','','')
    ALLOCATE (mcw_grid(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable mcw_grid','','')
    ALLOCATE (mcs_grid(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable mcs_grid','','')
    ALLOCATE (run2peat(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable run2peat','','')
    ALLOCATE (wt_ab(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable wt_ab','','')
    ALLOCATE (wt_ab_tide(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable wt_ab_tide','','')
    ALLOCATE (run2man(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable run2man','','')

       
    !! 2.3 Read in the run.def the parameters values defined by the user

    !! Campbell parameterisation sometimes used in the muff

    !Config Key   = B_MUFF
    !Config Desc  = Campbell coefficient b
    !Config If    = 
    !Config Def   = 2.708, 2.708, 2.708
    !Config Help  = This parameter will be constant over the entire 
    !Config         simulated domain, thus independent from soil
    !Config         texture.   
    !Config Units = [-]
    b_muff(:) = b_muff_param(:)  
    CALL getin_p("B_MUFF",b_muff)

    !! Check parameter value (correct range)
    IF ( ANY(b_muff(:) <= zero) ) THEN
       CALL ipslerr_p(error_level, "hydrol_init.", &
            &     "Wrong parameter value for B_MUFF.", &
            &     "This parameter should be positive. ", &
            &     "Please, check parameter value in run.def. ")
    END IF

    !Config Key   = PSI_AIR_ENTRY
    !Config Desc  = Air entry water potential
    !Config If    = 
    !Config Def   = -0.015, -0.015, -0.015
    !Config Help  = This parameter will be constant over the entire 
    !Config         simulated domain, thus independent from soil
    !Config         texture.   
    !Config Units = [-]
    psi_air_entry(:) = psi_air_entry_param(:)
    CALL getin_p("PSI_AIR_ENTRY",psi_air_entry)

    !! Check parameter value (correct range)
    IF ( ANY(psi_air_entry(:) >= zero) ) THEN
       CALL ipslerr_p(error_level, "hydrol_init.", &
            &     "Wrong parameter value for PSI_AIR_ENTRY.", &
            &     "This parameter should be positive. ", &
            &     "Please, check parameter value in run.def. ")
    END IF

    !Config Key   = IS_VG
    !Config Desc  = Flag to control the calculation of the hydraulic conductivity and water potential
    !Config If    = 
    !Config Def   = TRUE
    !Config Help  = This parameter will be constant over the entire 
    !Config         simulated domain, thus independent from soil
    !Config         texture.   
    !Config Units = [FLAG]
    is_vg = .FALSE.
    CALL getin_p("IS_VG",is_vg)

    IF (split_soil_properties) THEN
       mcr_sup(:) = mcr_sup_param
       mcr_inf(:) = mcr_inf_param
       mcs_sup(:) = mcs_sup_param
       mcs_inf(:) = mcs_inf_param
    ELSE       
       mcr_sup(:) = mcr(:)
       mcr_inf(:) = mcr(:)
       mcs_sup(:) = mcs(:)
       mcs_inf(:) = mcs(:)
    ENDIF



    !Config Key   = WETNESS_TRANSPIR_MAX
    !Config Desc  = Soil moisture above which transpir is max, for each soil texture class 
    !Config If    = 
    !Config Def   = 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8
    !Config Help  = This parameter is independent from soil texture for
    !Config         the time being.
    !Config Units = [-]
    pcent(:) = pcent_usda(:)  
    CALL getin_p("WETNESS_TRANSPIR_MAX",pcent)

    !! Check parameter value (correct range)
    IF ( ANY(pcent(:) <= zero) .OR. ANY(pcent(:) > 1.) ) THEN
       CALL ipslerr_p(error_level, "hydrol_init.", &
            &     "Wrong parameter value for WETNESS_TRANSPIR_MAX.", &
            &     "This parameter should be positive and less or equals than 1. ", &
            &     "Please, check parameter value in run.def. ")
    END IF


    !Config Key   = VWC_MIN_FOR_WET_ALB
    !Config Desc  = Vol. wat. cont. above which albedo is cst
    !Config If    = 
    !Config Def   = 0.25, 0.25, 0.25
    !Config Help  = This parameter is independent from soil texture for
    !Config         the time being.
    !Config Units = [m3/m3]
    mc_awet(:) = mc_awet_usda(:)
    CALL getin_p("VWC_MIN_FOR_WET_ALB",mc_awet)

    !! Check parameter value (correct range)
    IF ( ANY(mc_awet(:) < 0) ) THEN
       CALL ipslerr_p(error_level, "hydrol_init.", &
            &     "Wrong parameter value for VWC_MIN_FOR_WET_ALB.", &
            &     "This parameter should be positive. ", &
            &     "Please, check parameter value in run.def. ")
    END IF


    !Config Key   = VWC_MAX_FOR_DRY_ALB
    !Config Desc  = Vol. wat. cont. below which albedo is cst
    !Config If    = 
    !Config Def   = 0.1, 0.1, 0.1
    !Config Help  = This parameter is independent from soil texture for
    !Config         the time being.
    !Config Units = [m3/m3]
    mc_adry(:) = mc_adry_usda(:)
    CALL getin_p("VWC_MAX_FOR_DRY_ALB",mc_adry)

    !! Check parameter value (correct range)
    IF ( ANY(mc_adry(:) < 0) .OR. ANY(mc_adry(:) > mc_awet(:)) ) THEN
       CALL ipslerr_p(error_level, "hydrol_init.", &
            &     "Wrong parameter value for VWC_MAX_FOR_DRY_ALB.", &
            &     "This parameter should be positive and not greater than VWC_MIN_FOR_WET_ALB.", &
            &     "Please, check parameter value in run.def. ")
    END IF


    !! 3 Other array allocation
    ALLOCATE (mask_veget(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable mask_veget','','')

    ALLOCATE (mask_soiltile(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable mask_soiltile','','')

    ALLOCATE (humrelv(kjpindex,nvm,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable humrelv','','')

    ALLOCATE (vegstressv(kjpindex,nvm,nstm),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable vegstressv','','')

    ALLOCATE (precisol(kjpindex,nvm),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable precisol','','')

    ALLOCATE (throughfall(kjpindex,nvm),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable throughfall','','')

    ALLOCATE (precisol_ns(kjpindex,nstm),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable precisol_nc','','')

    ALLOCATE (free_drain_coef(kjpindex,nstm),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable free_drain_coef','','')

    ALLOCATE (zwt_force(kjpindex,nstm),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable zwt_force','','')

    ALLOCATE (frac_bare_ns(kjpindex,nstm),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable frac_bare_ns','','')

    ALLOCATE (water2infilt(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable water2infilt','','')

    ALLOCATE (ae_ns(kjpindex,nstm),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable ae_ns','','')

    ALLOCATE (rootsink(kjpindex,nslm,nstm),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable rootsink','','')

    ALLOCATE (subsnowveg(kjpindex),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable subsnowveg','','')

    ALLOCATE (subsnownobio(kjpindex,nnobio),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable subsnownobio','','')

    ALLOCATE (icemelt(kjpindex),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable icemelt','','')

    ALLOCATE (subsinksoil(kjpindex),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable subsinksoil','','')

    ALLOCATE (mx_eau_var(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable mx_eau_var','','')

    ALLOCATE (vegtot(kjpindex),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable vegtot','','')

    ALLOCATE (vegtot_old(kjpindex),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable vegtot_old','','')

    ALLOCATE (resdist(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable resdist','','')

    ALLOCATE (humtot(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable humtot','','')

    ALLOCATE (resolv(kjpindex),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable resolv','','')

    ALLOCATE (k(kjpindex,nslm),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable k','','')

    ALLOCATE (kk_moy(kjpindex,nslm),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable kk_moy','','')
    kk_moy(:,:) = 276.48
    
    ALLOCATE (kk(kjpindex,nslm,nstm),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable kk','','')
    kk(:,:,:) = 276.48
    
    ALLOCATE (avan_mod_tab(nslm,kjpindex),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable avan_mod_tab','','')
    
    ALLOCATE (nvan_mod_tab(nslm,kjpindex),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable nvan_mod_tab','','')

    ALLOCATE (a(kjpindex,nslm),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable a','','')

    ALLOCATE (b(kjpindex,nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable b','','')

    ALLOCATE (d(kjpindex,nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable d','','')

    ALLOCATE (e(kjpindex,nslm),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable e','','')

    ALLOCATE (f(kjpindex,nslm),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable f','','')

    ALLOCATE (g1(kjpindex,nslm),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable g1','','')

    ALLOCATE (ep(kjpindex,nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable ep','','')

    ALLOCATE (fp(kjpindex,nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable fp','','')

    ALLOCATE (gp(kjpindex,nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable gp','','')

    ALLOCATE (rhs(kjpindex,nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable rhs','','')

    ALLOCATE (srhs(kjpindex,nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable srhs','','')

    ALLOCATE (tmc(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable tmc','','')

    ALLOCATE (tmcs(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable tmcs','','')

    ALLOCATE (tmcr(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable tmcr','','')

    ALLOCATE (tmcfc(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable tmcfc','','')

    ALLOCATE (tmcw(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable tmcw','','')

    ALLOCATE (tmc_litter(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable tmc_litter','','')

    ALLOCATE (tmc_litt_mea(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable tmc_litt_mea','','')

    ALLOCATE (tmc_litter_res(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable tmc_litter_res','','')

    ALLOCATE (tmc_litter_wilt(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable tmc_litter_wilt','','')

    ALLOCATE (tmc_litter_field(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable tmc_litter_field','','')

    ALLOCATE (tmc_litter_sat(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable tmc_litter_sat','','')

    ALLOCATE (tmc_litter_awet(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable tmc_litter_awet','','')

    ALLOCATE (tmc_litter_adry(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable tmc_litter_adry','','')

    ALLOCATE (tmc_litt_wet_mea(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable tmc_litt_wet_mea','','')

    ALLOCATE (tmc_litt_dry_mea(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable tmc_litt_dry_mea','','')

    ALLOCATE (v1(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable v1','','')

    ALLOCATE (ru_ns(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable ru_ns','','')
    ru_ns(:,:) = zero

    ALLOCATE (dr_ns(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable dr_ns','','')
    dr_ns(:,:) = zero

    ALLOCATE (tr_ns(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable tr_ns','','')

    ALLOCATE (vegetmax_soil(kjpindex,nvm,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable vegetmax_soil','','')

    ALLOCATE (mc(kjpindex,nslm,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable mc','','')

    ALLOCATE (root_mc_fc(kjpindex),stat=ier) !
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable root_mc_fc','','') !

    ALLOCATE (nslm_root(kjpindex),stat=ier) !
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable nslm_root','','') !

    ALLOCATE (is_irrigated(nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable is_irrigated','','')


    ! Variables for nudging of soil moisture
    IF (ok_nudge_mc) THEN
       ALLOCATE (mc_read_prev(kjpindex,nslm,nstm),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable mc_read_prev','','')
       ALLOCATE (mc_read_next(kjpindex,nslm,nstm),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable mc_read_next','','')
       ALLOCATE (mc_read_current(kjpindex,nslm,nstm),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable mc_read_current','','')
       ALLOCATE (mask_mc_interp(kjpindex,nslm,nstm),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable mask_mc_interp','','')
       ALLOCATE (tmc_aux(kjpindex,nstm),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable tmc_aux','','')
    END IF

    ! Variables for nudging of snow variables
    IF (ok_nudge_snow) THEN
       ALLOCATE (snowdz_read_prev(kjpindex,nsnow),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable snowdz_read_prev','','')
       ALLOCATE (snowdz_read_next(kjpindex,nsnow),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable snowdz_read_next','','')
       
       ALLOCATE (snowrho_read_prev(kjpindex,nsnow),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable snowrho_read_prev','','')
       ALLOCATE (snowrho_read_next(kjpindex,nsnow),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable snowrho_read_next','','')
       
       ALLOCATE (snowtemp_read_prev(kjpindex,nsnow),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable snowtemp_read_prev','','')
       ALLOCATE (snowtemp_read_next(kjpindex,nsnow),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable snowtemp_read_next','','')
       
       ALLOCATE (mask_snow_interp(kjpindex,nsnow),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable mask_snow_interp','','')
    END IF

    ALLOCATE (mcl(kjpindex, nslm, nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable mcl','','')

    IF (ok_freeze_cwrr) THEN
       ALLOCATE (profil_froz_hydro(kjpindex, nslm),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable profil_froz_hydrol','','')
       profil_froz_hydro(:,:) = zero
    ENDIF
    
    ALLOCATE (profil_froz_hydro_ns(kjpindex, nslm, nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable profil_froz_hydro_ns','','')
    profil_froz_hydro_ns(:,:,:) = zero
    
    ALLOCATE (soilmoist(kjpindex,nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable soilmoist','','')

    ALLOCATE (soilmoist_s(kjpindex,nslm,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable soilmoist_s','','')

    ALLOCATE (soilmoist_liquid(kjpindex,nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable soilmoist_liquid','','')

    ALLOCATE (soil_wet_ns(kjpindex,nslm,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable soil_wet_ns','','')

    ALLOCATE (soil_wet_litter(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable soil_wet_litter','','')

    ALLOCATE (qflux_ns(kjpindex,nslm,nstm),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable qflux_ns','','')

    ALLOCATE (check_top_ns(kjpindex,nstm),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable check_top_ns','','')

    ALLOCATE (tmat(kjpindex,nslm,3),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable tmat','','')

    ALLOCATE (stmat(kjpindex,nslm,3),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable stmat','','')

    ALLOCATE (kfact_root(kjpindex, nslm, nstm), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable kfact_root','','')

    ALLOCATE (kfact(nslm, kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable kfact','','')

    ALLOCATE (zz(nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable zz','','')

    ALLOCATE (dz(nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable dz','','')
    
    ALLOCATE (dh(nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable dh','','')

    ALLOCATE (mc_lin(imin:imax, kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable mc_lin','','')

    ALLOCATE (k_lin(imin:imax, nslm, kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable k_lin','','')

    ALLOCATE (d_lin(imin:imax, nslm, kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable d_lin','','')

    ALLOCATE (a_lin(imin:imax, nslm, kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable a_lin','','')

    ALLOCATE (b_lin(imin:imax, nslm, kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable b_lin','','')

    ALLOCATE (undermcr(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable undermcr','','')

    ALLOCATE (tot_watveg_beg(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable tot_watveg_beg','','')
    
    ALLOCATE (tot_watveg_end(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable tot_watvag_end','','')
    
    ALLOCATE (tot_watsoil_beg(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable tot_watsoil_beg','','')
    
    ALLOCATE (tot_watsoil_end(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable tot_watsoil_end','','')
    
    ALLOCATE (delsoilmoist(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable delsoilmoist','','')
    
    ALLOCATE (delintercept(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable delintercept','','')
    
    ALLOCATE (delswe(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable delswe','','')
    
    ALLOCATE (snow_beg(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable snow_beg','','')
    
    ALLOCATE (snow_end(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable snow_end','','')
    
    !! 4 Open restart input file and read data for HYDROLOGIC process
    IF (printlev>=3) WRITE (numout,*) ' we have to read a restart file for HYDROLOGIC variables'

    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME', 'moistc')
    CALL restget_p (rest_id, 'moistc', nbp_glo, nslm , nstm, kjit, .TRUE., mc, "gather", nbp_glo, index_g)

    IF (ok_nudge_mc) THEN
       CALL ioconf_setatt_p('LONG_NAME','Soil moisture read from nudging file')
       CALL restget_p (rest_id, 'mc_read_next', nbp_glo, nslm , nstm, kjit, .TRUE., mc_read_next, &
            "gather", nbp_glo, index_g)
    END IF

    IF (ok_nudge_snow) THEN
       CALL ioconf_setatt_p('UNITS', 'm')
       CALL ioconf_setatt_p('LONG_NAME','Snow layer thickness read from nudging file')
       CALL restget_p (rest_id, 'snowdz_read_next', nbp_glo, nsnow, 1, kjit, .TRUE., snowdz_read_next, &
            "gather", nbp_glo, index_g)

       CALL ioconf_setatt_p('UNITS', 'kg/m^3')
       CALL ioconf_setatt_p('LONG_NAME','Snow density profile read from nudging file')
       CALL restget_p (rest_id, 'snowrho_read_next', nbp_glo, nsnow, 1, kjit, .TRUE., snowrho_read_next, &
            "gather", nbp_glo, index_g)

       CALL ioconf_setatt_p('UNITS', 'K')
       CALL ioconf_setatt_p('LONG_NAME','Snow temperature read from nudging file')
       CALL restget_p (rest_id, 'snowtemp_read_next', nbp_glo, nsnow, 1, kjit, .TRUE., snowtemp_read_next, &
            "gather", nbp_glo, index_g)
    END IF

    CALL restget_p (rest_id, 'moistcl', nbp_glo, nslm , nstm, kjit, .TRUE., mcl, "gather", nbp_glo, index_g)
    !
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','us')
    CALL restget_p (rest_id, 'us', nbp_glo, nvm, nstm, nslm, kjit, .TRUE., us, "gather", nbp_glo, index_g)
    !
    var_name= 'free_drain_coef'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Coefficient for free drainage at bottom of soil')
    CALL restget_p (rest_id, var_name, nbp_glo, nstm, 1, kjit, .TRUE., free_drain_coef, "gather", nbp_glo, index_g)
    !
    var_name= 'zwt_force'
    CALL ioconf_setatt_p('UNITS', 'm')
    CALL ioconf_setatt_p('LONG_NAME','Prescribed water table depth')
    CALL restget_p (rest_id, var_name, nbp_glo, nstm, 1, kjit, .TRUE., zwt_force, "gather", nbp_glo, index_g)
    !
    var_name= 'water2infilt'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Remaining water to be infiltrated on top of the soil')
    CALL restget_p (rest_id, var_name, nbp_glo, nstm, 1, kjit, .TRUE., water2infilt, "gather", nbp_glo, index_g)
    !
    var_name= 'ae_ns'
    CALL ioconf_setatt_p('UNITS', 'kg/m^2')
    CALL ioconf_setatt_p('LONG_NAME','Bare soil evap on each soil type')
    CALL restget_p (rest_id, var_name, nbp_glo, nstm, 1, kjit, .TRUE., ae_ns, "gather", nbp_glo, index_g)
    !
    var_name= 'snow'        
    CALL ioconf_setatt_p('UNITS', 'kg/m^2')
    CALL ioconf_setatt_p('LONG_NAME','Snow mass')
    CALL restget_p (rest_id, var_name, nbp_glo, 1  , 1, kjit, .TRUE., snow, "gather", nbp_glo, index_g)
    !
    var_name= 'snow_age'
    CALL ioconf_setatt_p('UNITS', 'd')
    CALL ioconf_setatt_p('LONG_NAME','Snow age')
    CALL restget_p (rest_id, var_name, nbp_glo, 1  , 1, kjit, .TRUE., snow_age, "gather", nbp_glo, index_g)
    !
    var_name= 'snow_nobio_age'
    CALL ioconf_setatt_p('UNITS', 'd')
    CALL ioconf_setatt_p('LONG_NAME','Snow age on other surface types')
    CALL restget_p (rest_id, var_name, nbp_glo, nnobio  , 1, kjit, .TRUE., snow_nobio_age, "gather", nbp_glo, index_g)
    !
    var_name= 'qsintveg'
    CALL ioconf_setatt_p('UNITS', 'kg/m^2')
    CALL ioconf_setatt_p('LONG_NAME','Intercepted moisture')
    CALL restget_p (rest_id, var_name, nbp_glo, nvm, 1, kjit, .TRUE., qsintveg, "gather", nbp_glo, index_g)

    var_name= 'evap_bare_lim_ns'
    CALL ioconf_setatt_p('UNITS', '?')
    CALL ioconf_setatt_p('LONG_NAME','?')
    CALL restget_p (rest_id, var_name, nbp_glo, nstm, 1, kjit, .TRUE., evap_bare_lim_ns, "gather", nbp_glo, index_g)
    CALL setvar_p (evap_bare_lim_ns, val_exp, 'NO_KEYWORD', 0.0)

    var_name= 'resdist'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','soiltile values from previous time-step')
    CALL restget_p (rest_id, var_name, nbp_glo, nstm, 1, kjit, .TRUE., resdist, "gather", nbp_glo, index_g)

    var_name= 'vegtot_old'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','vegtot from previous time-step')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., vegtot_old, "gather", nbp_glo, index_g)       

    ! Read drysoil_frac. It will be initalized later in hydrol_var_init if the varaible is not find in restart file.
    CALL ioconf_setatt_p('UNITS', '')
    CALL ioconf_setatt_p('LONG_NAME','Function of litter wetness')
    CALL restget_p (rest_id, 'drysoil_frac', nbp_glo, 1  , 1, kjit, .TRUE., drysoil_frac, "gather", nbp_glo, index_g)

    IF (ok_peat_hydro) THEN
       var_name= 'run2peat'
       CALL ioconf_setatt('UNITS', 'mm/d')
       CALL ioconf_setatt('LONG_NAME','runoff from soiltile1-3')
       CALL restget_p (rest_id, var_name, nbp_glo, 1  , 1, kjit, .TRUE.,run2peat , "gather", nbp_glo, index_g)
       
       var_name= 'wt_ab'
       CALL ioconf_setatt('UNITS', 'mm')
       CALL ioconf_setatt('LONG_NAME','water above surface')
       CALL restget_p (rest_id, var_name, nbp_glo, 1  , 1, kjit, .TRUE.,wt_ab, "gather", nbp_glo, index_g)
       
       var_name= 'wtp'
       CALL ioconf_setatt_p('UNITS', 'mm')
       CALL ioconf_setatt_p('LONG_NAME','mean water table position per soiltile')
       CALL restget_p (rest_id, var_name, nbp_glo, 1  , 1, kjit, .TRUE.,wtp , "gather", nbp_glo, index_g)

       var_name= 'liqwt_ratio'
       CALL ioconf_setatt('UNITS', '-')
       CALL ioconf_setatt('LONG_NAME','Liquid water for peat veg growth')
       CALL restget_p (rest_id, var_name, nbp_glo, 1  , 1, kjit, .TRUE.,liqwt_ratio, "gather", nbp_glo, index_g)
   ENDIF
   
   IF (tides) THEN
       var_name= 'wt_ab_tide'
       CALL ioconf_setatt('UNITS', 'mm')
       CALL ioconf_setatt('LONG_NAME','water above surface tide')
       CALL restget_p (rest_id, var_name, nbp_glo, 1  , 1, kjit, .TRUE.,wt_ab_tide, "gather", nbp_glo, index_g)
       
       var_name= 'run2man'
       CALL ioconf_setatt('UNITS', 'mm/d')
       CALL ioconf_setatt('LONG_NAME','runoff from soiltile1-3 to mangroves')
       CALL restget_p (rest_id, var_name, nbp_glo, 1  , 1, kjit, .TRUE.,run2man, "gather", nbp_glo, index_g)
       
    ENDIF


    IF (ok_hydrol_arch) THEN
       var_name='mc_out' 
       CALL restget_p (rest_id, var_name, nbp_glo, nslm, nstm,  kjit, .TRUE., mc_out, "gather", nbp_glo, index_g)
       IF (ALL(mc_out(:,:,:) == val_exp)) mc_out(:,:,:) = zero 

       var_name='ksoil'
       CALL restget_p (rest_id, var_name, nbp_glo, nslm, nstm,  kjit, .TRUE., ksoil, "gather", nbp_glo, index_g)
       IF (ALL(ksoil(:,:,:) == val_exp)) ksoil(:,:,:) = min_sechiba
    ENDIF

    !! 5 get restart values if none were found in the restart file
    !
    !Config Key   = HYDROL_MOISTURE_CONTENT
    !Config Desc  = Soil moisture on each soil tile and levels
    !Config If    = 
    !Config Def   = 0.3
    !Config Help  = The initial value of mc if its value is not found
    !Config         in the restart file. This should only be used if the model is 
    !Config         started without a restart file.
    !Config Units = [m3/m3]
    !
    CALL setvar_p (mc, val_exp, 'HYDROL_MOISTURE_CONTENT', 0.3_r_std)
    ! peatland 
    ! Initialization of 0.5 for mc peatland.
    ! This value can be further discussed later
    ! mc: kjpindex, nslm, nstm
    IF (ok_peat_hydro .AND. ALL(mc(:,:,:)== 0.3_r_std) ) THEN
       IF (nstm >= 4 ) THEN
          DO jst = 4, nstm
             DO ji=1,kjpindex
                IF (soiltile(ji,jst) .GT. min_sechiba) THEN
                   mc(ji,:,jst)  = 0.5_r_std
                ENDIF
             ENDDO
          ENDDO
       ENDIF
    ENDIF

    ! Initialize mcl as mc if it is not found in the restart file
    IF ( ALL(mcl(:,:,:)==val_exp) ) THEN
       mcl(:,:,:) = mc(:,:,:)
    END IF

    !Config Key   = US_INIT
    !Config Desc  = US_NVM_NSTM_NSLM
    !Config If    = 
    !Config Def   = 0.0
    !Config Help  = The initial value of us (relative moisture) if its value is not found
    !Config         in the restart file. This should only be used if the model is 
    !Config         started without a restart file.
    !Config Units = [-]
    !
    DO jsl=1,nslm
       CALL setvar_p (us(:,:,:,jsl), val_exp, 'US_INIT', zero)
    ENDDO
    !
    !Config Key   = ZWT_FORCE
    !Config Desc  = Prescribed water depth, dimension nstm
    !Config If    = 
    !Config Def   = undef undef undef
    !Config Help  = The initial value of zwt_force if its value is not found
    !Config         in the restart file. undef corresponds to a case whith no forced WT. 
    !Config         This should only be used if the model is started without a restart file.
    !Config Units = [m]

    ALLOCATE (zwt_default(nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable zwt_default','','')
    zwt_default(:) = undef_sechiba
    CALL setvar_p (zwt_force, val_exp, 'ZWT_FORCE', zwt_default )

    zforce = .FALSE.
    DO jst=1,nstm
       IF (zwt_force(1,jst) <= zmaxh) zforce = .TRUE. ! AD16*** check if OK with vertical_soil
    ENDDO
    !
    !Config Key   = FREE_DRAIN_COEF
    !Config Desc  = Coefficient for free drainage at bottom, dimension nstm
    !Config If    = 
    !Config Def   = 1.0 1.0 1.0
    !Config Help  = The initial value of free drainage coefficient if its value is not found
    !Config         in the restart file. This should only be used if the model is 
    !Config         started without a restart file.
    !Config Units = [-]

    ALLOCATE (free_drain_max(nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable free_drain_max','','')
    !free_drain_max(:)=1.0 ! trunk orig
    !
    ! peatland
    IF (ok_peat_hydro ) THEN
       free_drain_max(:)=1.0
       free_drain_max(4)= zero
       
       IF (tides) THEN
          free_drain_max(6)= zero
       ENDIF
    ELSE
       free_drain_max(:)=1.0
    ENDIF

    IF (agri_peat .AND. agri_drain) THEN
       free_drain_max(5)=drain_factor
    ENDIF
    
    CALL setvar_p (free_drain_coef, val_exp, 'FREE_DRAIN_COEF', free_drain_max)

    IF (agri_peat .AND. agri_drain) THEN
       free_drain_coef(:,5)=free_drain_max(5)
    ENDIF

    IF (printlev>=2) WRITE (numout,*) ' hydrol_init => free_drain_coef = ',free_drain_coef(1,:)
    DEALLOCATE(free_drain_max)

    !
    !Config Key   = WATER_TO_INFILT
    !Config Desc  = Water to be infiltrated on top of the soil
    !Config If    = 
    !Config Def   = 0.0
    !Config Help  = The initial value of free drainage if its value is not found
    !Config         in the restart file. This should only be used if the model is 
    !Config         started without a restart file.
    !Config Units = [mm]
    !
    CALL setvar_p (water2infilt, val_exp, 'WATER_TO_INFILT', zero)
    !
    !Config Key   = EVAPNU_SOIL
    !Config Desc  = Bare soil evap on each soil if not found in restart
    !Config If    = 
    !Config Def   = 0.0
    !Config Help  = The initial value of bare soils evap if its value is not found
    !Config         in the restart file. This should only be used if the model is 
    !Config         started without a restart file.
    !Config Units = [mm]
    !
    CALL setvar_p (ae_ns, val_exp, 'EVAPNU_SOIL', zero)
    !
    !Config Key  = HYDROL_SNOW
    !Config Desc  = Initial snow mass if not found in restart
    !Config If    = OK_SECHIBA
    !Config Def   = 0.0
    !Config Help  = The initial value of snow mass if its value is not found
    !Config         in the restart file. This should only be used if the model is 
    !Config         started without a restart file.
    !Config Units =
    !
    CALL setvar_p (snow, val_exp, 'HYDROL_SNOW', zero)
    !
    !Config Key   = HYDROL_SNOWAGE
    !Config Desc  = Initial snow age if not found in restart
    !Config If    = OK_SECHIBA
    !Config Def   = 0.0
    !Config Help  = The initial value of snow age if its value is not found
    !Config         in the restart file. This should only be used if the model is 
    !Config         started without a restart file.
    !Config Units = ***
    !
    CALL setvar_p (snow_age, val_exp, 'HYDROL_SNOWAGE', zero)
    !
    !Config Key   = HYDROL_SNOW_NOBIO_AGE
    !Config Desc  = Initial snow age on ice, lakes, etc. if not found in restart
    !Config If    = OK_SECHIBA
    !Config Def   = 0.0
    !Config Help  = The initial value of snow age if its value is not found
    !Config         in the restart file. This should only be used if the model is 
    !Config         started without a restart file.
    !Config Units = ***
    !
    CALL setvar_p (snow_nobio_age, val_exp, 'HYDROL_SNOW_NOBIO_AGE', zero)
    !
    !Config Key   = HYDROL_QSV
    !Config Desc  = Initial water on canopy if not found in restart
    !Config If    = OK_SECHIBA
    !Config Def   = 0.0
    !Config Help  = The initial value of moisture on canopy if its value 
    !Config         is not found in the restart file. This should only be used if
    !Config         the model is started without a restart file. 
    !Config Units = [mm]
    !
    CALL setvar_p (qsintveg, val_exp, 'HYDROL_QSV', zero)

    ! peatland
    IF (ok_peat_hydro) THEN
       !
       !Config Key   = RUN2PEAT
       !Config Desc  = Runoff from other soil tiles to peat
       !Config If    = OK_PEAT_HYDRO (or OK_RUN2PEAT)
       !Config Def   = 0.0
       !Config Help  = The initial value of run2peat if its value 
       !Config         is not found in the restart file. This should only be used if
       !Config         the model is started without a restart file. 
       !Config Units = [mm]
       CALL setvar_p (run2peat, val_exp,'HYDROL_run2peat', zero)
       !
       !Config Key   = WT_AB
       !Config Desc  = Height of water reservatory above the surface
       !Config If    = OK_PEAT_HYDRO (or ok_wt_ab)
       !Config Def   = 0.0
       !Config Help  = The initial value of wt_ab if its value 
       !Config         is not found in the restart file. This should only be used if
       !Config         the model is started without a restart file. 
       !Config Units = [mm]
       CALL setvar_p (wt_ab, val_exp,'HYDROL_wt_ab', zero)
       !
       !Config Key   = WTP
       !Config Desc  = Water table position
       !Config If    = OK_PEAT_HYDRO 
       !Config Def   = 0.0
       !Config Help  = The initial value of wtp if its value 
       !Config         is not found in the restart file. This should only be used if
       !Config         the model is started without a restart file. 
       !Config Units = [mm]
       CALL setvar_p (wtp, val_exp,'HYDROL_wtp', zero)
       !Config Key   = liqwt_ratio
       !Config Desc  = Liquid water ratio in peat
       !Config If    = OK_PEAT_HYDRO 
       !Config Def   = 0.0
       !Config Help  = The initial value of liqwt_ratio if its value 
       !Config         is not found in the restart file. This should only be used if
       !Config         the model is started without a restart file. 
       !Config Units = [mm]
       CALL setvar_p (liqwt_ratio, val_exp,'HYDROL_liqwt_ratio', zero)   
       !Config Key   = wt_ab_tide
       !Config Desc  = Height of water reservatory above the surface when flag tides is activated
       !Config If    = TIDES
       !Config Def   = 0.0
       !Config Help  = The initial value of wt_ab_tide if its value 
       !Config         is not found in the restart file. This should only be used if
       !Config         the model is started without a restart file. 
       !Config Units = [mm]
    ENDIF
    IF (tides) THEN
       CALL setvar_p (wt_ab_tide, val_exp,'HYDROL_wt_ab_tide', zero)
       !Config Key   = run2man
       !Config Desc  = Runoff from other soil tiles to mangrove when the tides is activated
       !Config If    = TIDES
       !Config Def   = 0.0
       !Config Help  = The initial value of wt_ab_tide if its value 
       !Config         is not found in the restart file. This should only be used if
       !Config         the model is started without a restart file. 
       !Config Units = [mm]
       CALL setvar_p (run2man, val_exp,'HYDROL_run2man', zero)
    ENDIF


    !! 6 Vegetation array      
    !
    ! If resdist is not in restart file, initialize with soiltile
    IF ( MINVAL(resdist) .EQ.  MAXVAL(resdist) .AND. MINVAL(resdist) .EQ. val_exp) THEN
       resdist(:,:) = soiltile(:,:)
    ENDIF

    !
    !  Remember that it is only frac_nobio + SUM(veget_max(,:)) that is equal to 1. Thus we need vegtot
    !
    IF ( ALL(vegtot_old(:) == val_exp) ) THEN
       ! vegtot_old was not found in restart file
       DO ji = 1, kjpindex
          vegtot_old(ji) = SUM(veget_max(ji,:))
       ENDDO
    ENDIF

    ! In the initialization phase, vegtot must take the value from previous time-step. 
    ! This is because hydrol_main is done before veget_max is updated in the end of the time step. 
    vegtot(:) = vegtot_old(:)
    !
    !
    ! compute the masks for veget

    mask_veget(:,:) = 0
    mask_soiltile(:,:) = 0

    DO jst=1,nstm
       DO ji = 1, kjpindex
          IF(soiltile(ji,jst) .GT. min_sechiba) THEN
             mask_soiltile(ji,jst) = 1
          ENDIF
       END DO
    ENDDO

    DO jv = 1, nvm
       DO ji = 1, kjpindex
          IF(veget_max(ji,jv) .GT. min_sechiba) THEN
             mask_veget(ji,jv) = 1
          ENDIF
       END DO
    END DO

    humrelv(:,:,:) = SUM(us,dim=4)


    !! 7a. Set vegstress

    var_name= 'vegstress'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Vegetation growth moisture stress')
    CALL restget_p (rest_id, var_name, nbp_glo, nvm, 1, kjit, .TRUE., vegstress, "gather", nbp_glo, index_g)

    vegstressv(:,:,:) = humrelv(:,:,:)
    ! Calculate vegstress if it is not found in restart file
    IF (ALL(vegstress(:,:)==val_exp)) THEN
       DO jv=1,nvm
          DO ji=1,kjpindex
             vegstress(ji,jv)=vegstress(ji,jv) + vegstressv(ji,jv,pref_soil_veg(jv))
          END DO
       END DO
    END IF
    !! 7b. Set humrel   
    ! Read humrel from restart file
    var_name= 'humrel'
    IF (is_root_prc) THEN
       CALL ioconf_setatt_p('UNITS', '')
       CALL ioconf_setatt_p('LONG_NAME','Relative humidity')
    ENDIF
    CALL restget_p (rest_id, var_name, nbp_glo, nvm, 1, kjit, .TRUE., humrel, "gather", nbp_glo, index_g)

    ! Calculate humrel if it is not found in restart file
    IF (ALL(humrel(:,:)==val_exp)) THEN
       ! set humrel from humrelv, assuming equi-repartition for the first time step
       humrel(:,:) = zero
       DO jv=1,nvm
          DO ji=1,kjpindex
             humrel(ji,jv)=humrel(ji,jv) + humrelv(ji,jv,pref_soil_veg(jv))      
          END DO
       END DO
    END IF

    ! Read evap_bare_lim from restart file
    var_name= 'evap_bare_lim'
    IF (is_root_prc) THEN
       CALL ioconf_setatt_p('UNITS', '')
       CALL ioconf_setatt_p('LONG_NAME','Limitation factor for bare soil evaporation')
    ENDIF
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., evap_bare_lim, "gather", nbp_glo, index_g)

    ! Calculate evap_bare_lim if it was not found in the restart file.
    IF ( ALL(evap_bare_lim(:) == val_exp) ) THEN
       DO ji = 1, kjpindex
          evap_bare_lim(ji) =  SUM(evap_bare_lim_ns(ji,:)*vegtot(ji)*soiltile(ji,:))
       ENDDO
    END IF

    ! Read root profile from restart file (it is used in hydraul_arch before it is 
    ! calculated in hydrol.f90. Putting it in the restarts avoids zero values the
    ! first day of each year.
    var_name= 'root_profile_struc'
    IF (is_root_prc) THEN
       CALL ioconf_setatt_p('UNITS', '')
       CALL ioconf_setatt_p('LONG_NAME','Structural root profile')
    ENDIF
    CALL restget_p (rest_id, var_name, nbp_glo, nvm, nslm, kjit, .TRUE., root_profile(:,:,:,istruc), "gather", nbp_glo, index_g)
    IF (ALL(root_profile(:,:,:,istruc) == val_exp)) root_profile(:,:,:,istruc) = zero
    
    var_name= 'root_profile_func'
    IF (is_root_prc) THEN
       CALL ioconf_setatt_p('UNITS', '')
       CALL ioconf_setatt_p('LONG_NAME','Functional root profile')
    ENDIF
    CALL restget_p (rest_id, var_name, nbp_glo, nvm, nslm, kjit, .TRUE., root_profile(:,:,:,ifunc), "gather", nbp_glo, index_g)
    IF (ALL(root_profile(:,:,:,ifunc) == val_exp)) root_profile(:,:,:,ifunc) = zero
       
    IF (ok_hydrol_arch) THEN
       var_name= 'psi_leaf'
       IF (is_root_prc) THEN
          CALL ioconf_setatt_p('UNITS', '')
          CALL ioconf_setatt_p('LONG_NAME','Leaf water potential')
       ENDIF
       CALL restget_p (rest_id, var_name, nbp_glo, nvm, ncirc, kjit, .TRUE., psi_leaf, "gather", nbp_glo, index_g)
       
       ! Calculate psi_leaf if it is not found in restart file
       IF (ALL(psi_leaf(:,:,:)==val_exp)) THEN
          psi_leaf(:,:,:) = zero
       END IF
       
       var_name= 'psi_leaf_next'
       IF (is_root_prc) THEN
          CALL ioconf_setatt_p('UNITS', '')
          CALL ioconf_setatt_p('LONG_NAME','Approximated Leaf water potential at time step n+1')
       ENDIF
       CALL restget_p (rest_id, var_name, nbp_glo, nvm, 1, kjit, .TRUE., psi_leaf_next, "gather", nbp_glo, index_g)
       
       ! Calculate psi_leaf_next if it is not found in restart file
       IF (ALL(psi_leaf_next(:,:)==val_exp)) THEN
          psi_leaf_next(:,:) = zero
       END IF

       var_name= 'psi_leaf_midday'
       IF (is_root_prc) THEN
          CALL ioconf_setatt_p('UNITS', '')
          CALL ioconf_setatt_p('LONG_NAME','Midday leaf water potential')
       ENDIF
       CALL restget_p (rest_id, var_name, nbp_glo, nvm, ncirc, kjit, .TRUE., psi_leaf_midday, "gather", nbp_glo, index_g)

       ! Calculate psi_leaf_midday if it is not found in restart file
       IF (ALL(psi_leaf_midday(:,:,:)==val_exp)) THEN
          psi_leaf_midday(:,:,:) = zero
       END IF

       var_name= 'psi_sto_leaf_save'
       IF (is_root_prc) THEN
          CALL ioconf_setatt_p('UNITS', '')
          CALL ioconf_setatt_p('LONG_NAME','Leaf storage water potential at time step n-1')
       ENDIF
       CALL restget_p (rest_id, var_name, nbp_glo, nvm, ncirc, kjit, .TRUE., psi_sto_leaf_save, "gather", nbp_glo, index_g)
       
       ! Calculate psi_sto_leaf_save if it is not found in restart file
       IF (ALL(psi_sto_leaf_save(:,:,:)==val_exp)) THEN
          psi_sto_leaf_save(:,:,:) = zero
       END IF
       
       var_name= 'psi_sto_wood_save'
       IF (is_root_prc) THEN
          CALL ioconf_setatt_p('UNITS', '')
          CALL ioconf_setatt_p('LONG_NAME','Wood storage water potential at time step n-1')
       ENDIF
       CALL restget_p (rest_id, var_name, nbp_glo, nvm, ncirc, kjit, .TRUE., psi_sto_wood_save, "gather", nbp_glo, index_g)
       
       ! Calculate psi_sto_wood_save if it is not found in restart file
       IF (ALL(psi_sto_wood_save(:,:,:)==val_exp)) THEN
          psi_sto_wood_save(:,:,:) = zero
       END IF
       
       var_name= 'psi_root_sup'
       IF (is_root_prc) THEN
          CALL ioconf_setatt_p('UNITS', '')
          CALL ioconf_setatt_p('LONG_NAME','Superficial roots water potential')
       ENDIF
       CALL restget_p (rest_id, var_name, nbp_glo, nvm, ncirc, kjit, .TRUE., psi_root_sup, "gather", nbp_glo, index_g)
       
       ! Calculate psi_root_sup if it is not found in restart file
       IF (ALL(psi_root_sup(:,:,:)==val_exp)) THEN
          psi_root_sup(:,:,:) = zero
       END IF
       
       var_name= 'psi_root_inf'
       IF (is_root_prc) THEN
          CALL ioconf_setatt_p('UNITS', '')
          CALL ioconf_setatt_p('LONG_NAME','Inferior roots water potential')
       ENDIF
       CALL restget_p (rest_id, var_name, nbp_glo, nvm, ncirc, kjit, .TRUE., psi_root_inf, "gather", nbp_glo, index_g)
       
       ! Calculate psi_root_inf if it is not found in restart file
       IF (ALL(psi_root_inf(:,:,:)==val_exp)) THEN
          psi_root_inf(:,:,:) = zero
       END IF

       var_name= 'psi_soil_sup'
       IF (is_root_prc) THEN
          CALL ioconf_setatt_p('UNITS', '')
          CALL ioconf_setatt_p('LONG_NAME','Superficial soil water potential')
       ENDIF
       CALL restget_p (rest_id, var_name, nbp_glo, nstm, 1, kjit, .TRUE., psi_soil_sup, "gather", nbp_glo, index_g)

       ! Calculate psi_soil_sup if it is not found in restart file
       IF (ALL(psi_soil_sup(:,:)==val_exp)) THEN
          psi_soil_sup(:,:) = zero
       END IF

       var_name= 'psi_soil_inf'
       IF (is_root_prc) THEN
          CALL ioconf_setatt_p('UNITS', '')
          CALL ioconf_setatt_p('LONG_NAME','Inferior soil water potential')
       ENDIF
       CALL restget_p (rest_id, var_name, nbp_glo, nstm, 1, kjit, .TRUE., psi_soil_inf, "gather", nbp_glo, index_g)

       ! Calculate psi_soil_inf if it is not found in restart file
       IF (ALL(psi_soil_inf(:,:)==val_exp)) THEN
          psi_soil_inf(:,:) = zero
       END IF
 
       var_name= 'psi_xylem_trunk'
       IF (is_root_prc) THEN
          CALL ioconf_setatt_p('UNITS', '')
          CALL ioconf_setatt_p('LONG_NAME','Xylem (trunk level) water potential')
       ENDIF
       
       CALL restget_p (rest_id, var_name, nbp_glo, nvm, ncirc, kjit, .TRUE., psi_xylem_trunk, "gather", nbp_glo, index_g)
       
       ! Calculate psi_xylem_trunk if it is not found in restart file
       IF (ALL(psi_xylem_trunk(:,:,:)==val_exp)) THEN
          psi_xylem_trunk(:,:,:) = zero
       END IF
       
       var_name= 'psi_xylem_leaf'
       IF (is_root_prc) THEN
          CALL ioconf_setatt_p('UNITS', '')
          CALL ioconf_setatt_p('LONG_NAME','Xylem (leaf level) water potential')
       ENDIF
       CALL restget_p (rest_id, var_name, nbp_glo, nvm, ncirc, kjit, .TRUE., psi_xylem_leaf, "gather", nbp_glo, index_g)
       
       ! Calculate psi_xylem_leaf if it is not found in restart file
       IF (ALL(psi_xylem_leaf(:,:,:)==val_exp)) THEN
          psi_xylem_leaf(:,:,:) = zero
       END IF
       
       var_name= 'psi_xylem_collar'
       IF (is_root_prc) THEN
          CALL ioconf_setatt_p('UNITS', '')
          CALL ioconf_setatt_p('LONG_NAME','Xylem (collar level) water potential')
       ENDIF
       CALL restget_p (rest_id, var_name, nbp_glo, nvm, ncirc, kjit, .TRUE., psi_xylem_collar, "gather", nbp_glo, index_g)
       
       ! Calculate psi_xylem_collar if it is not found in restart file
       IF (ALL(psi_xylem_collar(:,:,:)==val_exp)) THEN
          psi_xylem_collar(:,:,:) = zero
       END IF
       
       var_name= 'psi_sto_wood'
       IF (is_root_prc) THEN
          CALL ioconf_setatt_p('UNITS', '')
          CALL ioconf_setatt_p('LONG_NAME','Storage wood water potential')
       ENDIF
       CALL restget_p (rest_id, var_name, nbp_glo, nvm, ncirc, kjit, .TRUE., psi_sto_wood, "gather", nbp_glo, index_g)
       
       ! Calculate psi_sto_wood if it is not found in restart file
       IF (ALL(psi_sto_wood(:,:,:)==val_exp)) THEN
          psi_sto_wood(:,:,:) = zero
       END IF
       
       var_name= 'psi_sto_leaf'
       IF (is_root_prc) THEN
          CALL ioconf_setatt_p('UNITS', '')
          CALL ioconf_setatt_p('LONG_NAME','Leaf storage water potential')
       ENDIF
       CALL restget_p (rest_id, var_name, nbp_glo, nvm, ncirc, kjit, .TRUE., psi_sto_leaf, "gather", nbp_glo, index_g)
       
       ! Calculate psi_sto_leaf if it is not found in restart file
       IF (ALL(psi_sto_leaf(:,:,:)==val_exp)) THEN
          psi_sto_leaf(:,:,:) = zero
       END IF
       
       var_name= 'F_absorption'
       IF (is_root_prc) THEN
          CALL ioconf_setatt_p('UNITS', '')
          CALL ioconf_setatt_p('LONG_NAME','Total root absorption')
       ENDIF
       CALL restget_p (rest_id, var_name, nbp_glo, nvm, 1, kjit, .TRUE., F_absorption, "gather", nbp_glo, index_g)
       
       ! Calculate F_absorption if it is not found in restart file
       IF (ALL(F_absorption(:,:)==val_exp)) THEN
          F_absorption(:,:) = zero
       END IF
       
       var_name= 'mc_i_sup'
       IF (is_root_prc) THEN
          CALL ioconf_setatt_p('UNITS', '')
          CALL ioconf_setatt_p('LONG_NAME','Water content at each node of the absorption muff of the superficial soil layer')
       ENDIF
       CALL restget_p (rest_id, var_name, nbp_glo, nvm, ncirc, nrp, kjit, .TRUE., mc_i_sup, "gather", nbp_glo, index_g)
       
       ! Calculate mc_i_sup if it is not found in restart file
       IF (ALL(mc_i_sup(:,:,:,:)==val_exp)) THEN
          mc_i_sup(:,:,:,:) = un
       END IF
       
       var_name= 'mc_i_inf'
       IF (is_root_prc) THEN
          CALL ioconf_setatt_p('UNITS', '')
          CALL ioconf_setatt_p('LONG_NAME','Water content at each node of the absorption muff of the inferior soil layer')
       ENDIF
       CALL restget_p (rest_id, var_name, nbp_glo, nvm, ncirc, nrp, kjit, .TRUE., mc_i_inf, "gather", nbp_glo, index_g)
       
       ! Calculate mc_i_inf if it is not found in restart file
       IF (ALL(mc_i_inf(:,:,:,:)==val_exp)) THEN
          mc_i_inf(:,:,:,:) = un
       END IF
    END IF ! ok_hydrol_arch


    ! Read from restart file       
    ! The variables tot_watsoil_beg, tot_watsoil_beg and snwo_beg will be initialized in the end of 
    ! hydrol_initialize if they were not found in the restart file.

    var_name= 'tot_watveg_beg'
    IF (is_root_prc) THEN
       CALL ioconf_setatt_p('UNITS', '?')
       CALL ioconf_setatt_p('LONG_NAME','?')
    ENDIF
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., tot_watveg_beg, "gather", nbp_glo, index_g)

    var_name= 'tot_watsoil_beg'
    IF (is_root_prc) THEN
       CALL ioconf_setatt_p('UNITS', '?')
       CALL ioconf_setatt_p('LONG_NAME','?')
    ENDIF
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., tot_watsoil_beg, "gather", nbp_glo, index_g)

    var_name= 'snow_beg'
    IF (is_root_prc) THEN
       CALL ioconf_setatt_p('UNITS', '?')
       CALL ioconf_setatt_p('LONG_NAME','?')
    ENDIF
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., snow_beg, "gather", nbp_glo, index_g)


    ! Initialize variables for explictsnow module by reading restart file
    CALL explicitsnow_initialize( kjit,  kjpindex,   rest_id,   frac_nobio,  &
         snowrho,   snowtemp,  snowdz,   snowheat,   snowgrain,              &
         icetemp,   icedz,     ice_sheet_mask)


    ! Initialize soil moisture for nudging if not found in restart file
    IF (ok_nudge_mc) THEN
       IF ( ALL(mc_read_next(:,:,:)==val_exp) ) mc_read_next(:,:,:) = mc(:,:,:)
    END IF

    ! Initialize snow variables for nudging if not found in restart file
    IF (ok_nudge_snow) THEN
       IF ( ALL(snowdz_read_next(:,:)==val_exp) ) snowdz_read_next(:,:) = snowdz(:,:)
       IF ( ALL(snowrho_read_next(:,:)==val_exp) ) snowrho_read_next(:,:) = snowrho(:,:)
       IF ( ALL(snowtemp_read_next(:,:)==val_exp) ) snowtemp_read_next(:,:) = snowtemp(:,:)
    END IF


    IF (printlev>=3) WRITE (numout,*) ' hydrol_init done '

  END SUBROUTINE hydrol_init


!! ================================================================================================================================
!! SUBROUTINE 	: hydrol_clear
!!
!>\BRIEF        Deallocate arrays 
!!
!_ ================================================================================================================================
!_ hydrol_clear

  SUBROUTINE hydrol_clear()
   
    ! Allocation for soiltile related parameters
    IF ( ALLOCATED (pcent)) DEALLOCATE (pcent)
    IF ( ALLOCATED (mc_awet)) DEALLOCATE (mc_awet)
    IF ( ALLOCATED (mc_adry)) DEALLOCATE (mc_adry)
    ! Other arrays
    IF (ALLOCATED (mask_veget)) DEALLOCATE (mask_veget)
    IF (ALLOCATED (mask_soiltile)) DEALLOCATE (mask_soiltile)
    IF (ALLOCATED (humrelv)) DEALLOCATE (humrelv)
    IF (ALLOCATED (vegstressv)) DEALLOCATE (vegstressv)
    IF (ALLOCATED  (precisol)) DEALLOCATE (precisol)
    IF (ALLOCATED  (throughfall)) DEALLOCATE (throughfall)
    IF (ALLOCATED  (precisol_ns)) DEALLOCATE (precisol_ns)
    IF (ALLOCATED  (free_drain_coef)) DEALLOCATE (free_drain_coef)
    IF (ALLOCATED  (frac_bare_ns)) DEALLOCATE (frac_bare_ns)
    IF (ALLOCATED  (water2infilt)) DEALLOCATE (water2infilt)
    IF (ALLOCATED  (ae_ns)) DEALLOCATE (ae_ns)
    IF (ALLOCATED  (rootsink)) DEALLOCATE (rootsink)
    IF (ALLOCATED  (subsnowveg)) DEALLOCATE (subsnowveg)
    IF (ALLOCATED  (subsnownobio)) DEALLOCATE (subsnownobio)
    IF (ALLOCATED  (icemelt)) DEALLOCATE (icemelt)
    IF (ALLOCATED  (subsinksoil)) DEALLOCATE (subsinksoil)
    IF (ALLOCATED  (mx_eau_var)) DEALLOCATE (mx_eau_var)
    IF (ALLOCATED  (vegtot)) DEALLOCATE (vegtot)
    IF (ALLOCATED  (vegtot_old)) DEALLOCATE (vegtot_old)
    IF (ALLOCATED  (resdist)) DEALLOCATE (resdist)
    IF (ALLOCATED  (tot_watveg_beg)) DEALLOCATE (tot_watveg_beg)
    IF (ALLOCATED  (tot_watveg_end)) DEALLOCATE (tot_watveg_end)
    IF (ALLOCATED  (tot_watsoil_beg)) DEALLOCATE (tot_watsoil_beg)
    IF (ALLOCATED  (tot_watsoil_end)) DEALLOCATE (tot_watsoil_end)
    IF (ALLOCATED  (delsoilmoist)) DEALLOCATE (delsoilmoist)
    IF (ALLOCATED  (delintercept)) DEALLOCATE (delintercept)
    IF (ALLOCATED  (snow_beg)) DEALLOCATE (snow_beg)
    IF (ALLOCATED  (snow_end)) DEALLOCATE (snow_end)
    IF (ALLOCATED  (delswe)) DEALLOCATE (delswe)
    IF (ALLOCATED  (undermcr)) DEALLOCATE (undermcr)
    IF (ALLOCATED  (v1)) DEALLOCATE (v1)
    IF (ALLOCATED  (humtot)) DEALLOCATE (humtot)
    IF (ALLOCATED  (resolv)) DEALLOCATE (resolv)
    IF (ALLOCATED  (k)) DEALLOCATE (k)
    IF (ALLOCATED  (kk)) DEALLOCATE (kk)
    IF (ALLOCATED  (kk_moy)) DEALLOCATE (kk_moy)
    IF (ALLOCATED  (avan_mod_tab)) DEALLOCATE (avan_mod_tab)
    IF (ALLOCATED  (nvan_mod_tab)) DEALLOCATE (nvan_mod_tab)
    IF (ALLOCATED  (a)) DEALLOCATE (a)
    IF (ALLOCATED  (b)) DEALLOCATE (b)
    IF (ALLOCATED  (d)) DEALLOCATE (d)
    IF (ALLOCATED  (e)) DEALLOCATE (e)
    IF (ALLOCATED  (f)) DEALLOCATE (f)
    IF (ALLOCATED  (g1)) DEALLOCATE (g1)
    IF (ALLOCATED  (ep)) DEALLOCATE (ep)
    IF (ALLOCATED  (fp)) DEALLOCATE (fp)
    IF (ALLOCATED  (gp)) DEALLOCATE (gp)
    IF (ALLOCATED  (rhs)) DEALLOCATE (rhs)
    IF (ALLOCATED  (srhs)) DEALLOCATE (srhs)
    IF (ALLOCATED  (tmc)) DEALLOCATE (tmc)
    IF (ALLOCATED  (tmcs)) DEALLOCATE (tmcs)
    IF (ALLOCATED  (tmcr)) DEALLOCATE (tmcr)
    IF (ALLOCATED  (tmcfc)) DEALLOCATE (tmcfc)
    IF (ALLOCATED  (tmcw)) DEALLOCATE (tmcw)
    IF (ALLOCATED  (tmc_litter)) DEALLOCATE (tmc_litter)
    IF (ALLOCATED  (tmc_litt_mea)) DEALLOCATE (tmc_litt_mea)
    IF (ALLOCATED  (tmc_litter_res)) DEALLOCATE (tmc_litter_res)
    IF (ALLOCATED  (tmc_litter_wilt)) DEALLOCATE (tmc_litter_wilt)
    IF (ALLOCATED  (tmc_litter_field)) DEALLOCATE (tmc_litter_field)
    IF (ALLOCATED  (tmc_litter_sat)) DEALLOCATE (tmc_litter_sat)
    IF (ALLOCATED  (tmc_litter_awet)) DEALLOCATE (tmc_litter_awet)
    IF (ALLOCATED  (tmc_litter_adry)) DEALLOCATE (tmc_litter_adry)
    IF (ALLOCATED  (tmc_litt_wet_mea)) DEALLOCATE (tmc_litt_wet_mea)
    IF (ALLOCATED  (tmc_litt_dry_mea)) DEALLOCATE (tmc_litt_dry_mea)
    IF (ALLOCATED  (ru_ns)) DEALLOCATE (ru_ns)
    IF (ALLOCATED  (dr_ns)) DEALLOCATE (dr_ns)
    IF (ALLOCATED  (tr_ns)) DEALLOCATE (tr_ns)
    IF (ALLOCATED  (vegetmax_soil)) DEALLOCATE (vegetmax_soil)
    IF (ALLOCATED  (mc)) DEALLOCATE (mc)
    IF (ALLOCATED  (root_mc_fc)) DEALLOCATE (root_mc_fc)
    IF (ALLOCATED  (nslm_root)) DEALLOCATE (nslm_root)
    IF (ALLOCATED  (is_irrigated)) DEALLOCATE (is_irrigated)
    IF (ALLOCATED  (soilmoist)) DEALLOCATE (soilmoist)
    IF (ALLOCATED  (soilmoist_s)) DEALLOCATE (soilmoist_s)
    IF (ALLOCATED  (soilmoist_liquid)) DEALLOCATE (soilmoist_liquid)
    IF (ALLOCATED  (soil_wet_ns)) DEALLOCATE (soil_wet_ns)
    IF (ALLOCATED  (soil_wet_litter)) DEALLOCATE (soil_wet_litter)
    IF (ALLOCATED  (qflux_ns)) DEALLOCATE (qflux_ns)
    IF (ALLOCATED  (tmat)) DEALLOCATE (tmat)
    IF (ALLOCATED  (stmat)) DEALLOCATE (stmat)
    IF (ALLOCATED  (kfact_root)) DEALLOCATE (kfact_root)
    IF (ALLOCATED  (kfact)) DEALLOCATE (kfact)
    IF (ALLOCATED  (zz)) DEALLOCATE (zz)
    IF (ALLOCATED  (dz)) DEALLOCATE (dz)
    IF (ALLOCATED  (dh)) DEALLOCATE (dh)
    IF (ALLOCATED  (mc_lin)) DEALLOCATE (mc_lin)
    IF (ALLOCATED  (k_lin)) DEALLOCATE (k_lin)
    IF (ALLOCATED  (d_lin)) DEALLOCATE (d_lin)
    IF (ALLOCATED  (a_lin)) DEALLOCATE (a_lin)
    IF (ALLOCATED  (b_lin)) DEALLOCATE (b_lin)

    !! peatland
    IF ( ALLOCATED (mcs_grid)) DEALLOCATE (mcs_grid)
    IF ( ALLOCATED (mcw_grid)) DEALLOCATE (mcw_grid)
    IF ( ALLOCATED (run2peat)) DEALLOCATE (run2peat)
    IF ( ALLOCATED (wt_ab)) DEALLOCATE (wt_ab)
    IF ( ALLOCATED (wt_ab_tide)) DEALLOCATE (wt_ab_tide)
    IF ( ALLOCATED (run2man)) DEALLOCATE (run2man)
    IF (ALLOCATED  (kfact_peat)) DEALLOCATE (kfact_peat)
    IF (ALLOCATED  (mc_lin_peat)) DEALLOCATE (mc_lin_peat)
    IF (ALLOCATED  (k_lin_peat)) DEALLOCATE (k_lin_peat)
    IF (ALLOCATED  (d_lin_peat)) DEALLOCATE (d_lin_peat)
    IF (ALLOCATED  (a_lin_peat)) DEALLOCATE (a_lin_peat)
    IF (ALLOCATED  (b_lin_peat)) DEALLOCATE (b_lin_peat)

  END SUBROUTINE hydrol_clear

!! ================================================================================================================================
!! SUBROUTINE 	: hydrol_tmc_update
!!
!>\BRIEF        This routine updates the soil moisture profiles when the vegetation fraction have changed. 
!!
!! DESCRIPTION  :
!! 
!!    This routine update tmc and mc with variation of veget_max (LAND_USE or DGVM activated)
!! 
!!
!!
!!
!! RECENT CHANGE(S) : Adaptation to excluding nobio from soiltile(1)
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE hydrol_tmc_update ( kjpindex, veget_max, soiltile, qsintveg, drain_upd, runoff_upd)

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                            :: kjpindex         !! domain size
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)     :: veget_max        !! max fraction of vegetation type
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT (in)   :: soiltile         !! Fraction of each soil tile (0-1, unitless)
    
    !! 0.2 Output variables
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)         :: drain_upd        !! Change in drainage due to decrease in vegtot
                                                                              !! on mc [kg/m2/dt]
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)         :: runoff_upd       !! Change in runoff due to decrease in vegtot
                                                                              !! on water2infilt[kg/m2/dt]
    
    !! 0.3 Modified variables
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout)  :: qsintveg         !! Amount of water in the canopy interception 

    !! 0.4 Local variables
    INTEGER(i_std)                             :: ji, jv, jst,jsl, index      !! Indices
    LOGICAL                                    :: soil_upd                    !! True if soiltile changed since last time step
    LOGICAL                                    :: vegtot_upd                  !! True if vegtot changed since last time step
    REAL(r_std), DIMENSION(kjpindex,nstm)      :: vmr                         !! Change in soiltile (within vegtot)
    REAL(r_std), DIMENSION(kjpindex)           :: vmr_sum
    REAL(r_std), DIMENSION(kjpindex)           :: delvegtot    
    REAL(r_std), DIMENSION(kjpindex,nslm)      :: mc_dilu                     !! Total loss of moisture content
    REAL(r_std), DIMENSION(kjpindex)           :: infil_dilu                  !! Total loss for water2infilt
    REAL(r_std), DIMENSION(kjpindex,nstm)      :: tmc_old                     !! tmc before calculations
    REAL(r_std), DIMENSION(kjpindex,nstm)      :: water2infilt_old            !! water2infilt before calculations
    REAL(r_std), DIMENSION (kjpindex,nvm)      :: qsintveg_old                !! qsintveg before calculations
    REAL(r_std), DIMENSION(kjpindex)           :: test
    REAL(r_std), DIMENSION(kjpindex,nslm,nstm) :: mcaux                       !! serves to hold the chnage in mc when vegtot decreases

    
    !! 1. Update canopy interception following a land cover change
    !     If a PFT has disapperead as result from a veget_max change, 
    !     the intercepted water will have been lost during the removal of the vegetation. 
    !     The water previously stored on the canopy will now be added to surface water.
    !     Other adaptations of qsintveg are delt by the normal functioning of hydrol_canop
    DO ji=1,kjpindex
       IF (vegtot_old(ji) .GT.min_sechiba) THEN
          DO jv=1,nvm
             IF ((veget_max(ji,jv).LT.min_sechiba).AND.(qsintveg(ji,jv).GT.0.)) THEN

                ! The PFT has been removed but there is still some water on the canopy a solution need to be
                ! found for this water. If it is a forest PFT that was removed we will just add the water to
                ! soil water column of the tall vegetation. Note that it is also possible that last forest
                ! was removed. In that case there is no longer a tall vegetation water column. In that case
                ! we need to find a different water column to add the canopy water to. Ideally that would be
                ! to water column to which the new PFT belongs. For example if the last forest became a cropland
                ! the water previously stored in the forest canopy should be added to the soil water column
                ! of the short vegetation. Because the current land cover change functionality only deals
                ! with net land cover changes we don know the exact changes. An approximation will be used.

                ! Search for a suitable soil tile index to move the canopy water into
                jst=pref_soil_veg(jv)
                IF(resdist(ji,jst).GT.zero)THEN
                   index = jst   
                ELSE
                   ! Note that dim=1 refers to the dimensions of the answer
                   index = MAXLOC(resdist(ji,:),DIM=1)
                   IF(resdist(ji,index).LE.zero)THEN
                      WRITE(numout,*) 'ipts, index, resdist, ', ji, index, resdist(ji,:)
                      CALL ipslerr_p(3,'hydrol_tmc_update','if all resdist - see above- are zero',&
                           'the last vegetation may have been replaced by a non biological land cover',&
                           'This transfer has not yet been implemented in the code')
                   END IF
                END IF

                ! Move the canopy water into the surface water
                water2infilt(ji,index) = water2infilt(ji,index) + qsintveg(ji,jv)/(resdist(ji,index)*vegtot_old(ji))
                qsintveg(ji,jv) = zero
                   
             END IF
          END DO
       END IF
    END DO
   
    !! 2. We now deal with the changes of soiltile and corresponding soil moistures
    !!    Because sum(soiltile)=1 whatever vegtot, we need to distinguish two cases:
    !!    - when vegtot changes (meaning that the nobio fraction changes too),
    !!    - and when vegtot does not changes (a priori the most frequent case)

    vegtot_upd = SUM(ABS((vegtot(:)-vegtot_old(:)))) .GT. zero ! True if at least one land point with a vegtot change
    runoff_upd(:) = zero
    drain_upd(:) = zero
    IF (vegtot_upd) THEN

       ! We find here the processing specific to the chnages of nobio fraction and vegtot
       delvegtot(:) = vegtot(:) - vegtot_old(:)

       DO jst=1,nstm
          DO ji=1,kjpindex

             IF (delvegtot(ji) .GT. min_sechiba) THEN

                !! 2.1. If vegtot increases (nobio decreases), then the mc in each soiltile is decreased
                !!      assuming the same proportions for each soiltile, and each soil layer
                
                mc(ji,:,jst) = mc(ji,:,jst) * vegtot_old(ji)/vegtot(ji) ! vegtot cannot be zero as > vegtot_old
                water2infilt(ji,jst) = water2infilt(ji,jst) * vegtot_old(ji)/vegtot(ji)

             ELSE

                !! 2.2 If vegtot decreases (nobio increases), then the mc in each soiltile should increase,
                !!     but should not exceed mcs
                !!     For simplicity, we choose to send the corresponding water volume to drainage
                !!     We do the same for water2infilt but send the excess to surface runoff

                IF (vegtot(ji) .GT.min_sechiba) THEN
                   mcaux(ji,:,jst) =  mc(ji,:,jst) * (vegtot_old(ji)-vegtot(ji))/vegtot(ji) ! mcaux is the delta mc
                ELSE ! we just have nobio in the grid-cell
                   mcaux(ji,:,jst) =  mc(ji,:,jst)
                ENDIF
                
                drain_upd(ji) = drain_upd(ji) + dz(2) * ( trois*mcaux(ji,1,jst) + mcaux(ji,2,jst) )/huit
                DO jsl = 2,nslm-1
                   drain_upd(ji) = drain_upd(ji) + dz(jsl) * (trois*mcaux(ji,jsl,jst)+mcaux(ji,jsl-1,jst))/huit &
                        + dz(jsl+1) * (trois*mcaux(ji,jsl,jst)+mcaux(ji,jsl+1,jst))/huit
                ENDDO
                drain_upd(ji) = drain_upd(ji) + dz(nslm) * (trois*mcaux(ji,nslm,jst) + mcaux(ji,nslm-1,jst))/huit

                IF (vegtot(ji) .GT.min_sechiba) THEN
                   runoff_upd(ji) = runoff_upd(ji) + water2infilt(ji,jst) * (vegtot_old(ji)-vegtot(ji))/vegtot(ji)
                ELSE ! we just have nobio in the grid-cell
                   runoff_upd(ji) = runoff_upd(ji) + water2infilt(ji,jst)
                ENDIF

             ENDIF
             
          ENDDO
       ENDDO
       
    ENDIF
    
    !! 3. At the end of step 2, we are back to a case where vegtot changes are treated, so we can use soiltile
    !!    as a fraction of vegtot to process the mc transfers between soil tiles due to the changes of vegetation map
   
    !! 3.1 Check if soiltiles changed since last time step
    soil_upd=SUM(ABS(soiltile(:,:)-resdist(:,:))) .GT. zero
    IF (printlev>=3) WRITE (numout,*) 'soil_upd ', soil_upd
        
    IF (soil_upd) THEN
     
       !! 3.2 Define the change in soiltile
       vmr(:,:) = soiltile(:,:) - resdist(:,:)  ! resdist is the previous values of soiltiles, previous timestep, so before new map

       ! Total area loss by the three soil tiles
       DO ji=1,kjpindex
          vmr_sum(ji)=SUM(vmr(ji,:),MASK=vmr(ji,:).LT.zero)
       ENDDO

       !! 3.3 Shrinking soil tiles
       !! 3.3.1 Total loss of moisture content from the shrinking soil tiles, expressed by soil layer
       mc_dilu(:,:)=zero
       DO jst=1,nstm
          DO jsl = 1, nslm
             DO ji=1,kjpindex
                IF ( vmr(ji,jst) < -min_sechiba ) THEN
                   mc_dilu(ji,jsl) = mc_dilu(ji,jsl) + mc(ji,jsl,jst) * vmr(ji,jst) / vmr_sum(ji)
                ENDIF
             ENDDO
          ENDDO
       ENDDO

       !! 3.3.2 Total loss of water2inft from the shrinking soil tiles
       infil_dilu(:)=zero
       DO jst=1,nstm
          DO ji=1,kjpindex
             IF ( vmr(ji,jst) < -min_sechiba ) THEN
                infil_dilu(ji) = infil_dilu(ji) + water2infilt(ji,jst) * vmr(ji,jst) / vmr_sum(ji)
             ENDIF
          ENDDO
       ENDDO

       !! 3.4 Each gaining soil tile gets moisture proportionally to both the total loss and its areal increase 

       ! As the original mc from each soil tile are in [mcr,mcs] and we do weighted avrage, the new mc are in [mcr,mcs]
       ! The case where the soiltile is created (soiltile_old=0) works as the other cases

       ! 3.4.1 Update mc(kjpindex,nslm,nstm) !m3/m3
       DO jst=1,nstm
          DO jsl = 1, nslm
             DO ji=1,kjpindex
                IF ( vmr(ji,jst) > min_sechiba ) THEN
                   mc(ji,jsl,jst) = ( mc(ji,jsl,jst) * resdist(ji,jst) + mc_dilu(ji,jsl) * vmr(ji,jst) ) / soiltile(ji,jst)
                   ! NB : soiltile can not be zero for case vmr > zero, see slowproc_veget
                ENDIF
             ENDDO
          ENDDO
       ENDDO
       
       ! 3.4.2 Update water2inft
       DO jst=1,nstm
          DO ji=1,kjpindex
             IF ( vmr(ji,jst) > min_sechiba ) THEN !donc soiltile>0     
                water2infilt(ji,jst) = ( water2infilt(ji,jst) * resdist(ji,jst) + infil_dilu(ji) * vmr(ji,jst) ) / soiltile(ji,jst)
             ENDIF !donc resdist>0
          ENDDO
       ENDDO

       ! 3.4.3 Case where soiltile < min_sechiba 
       DO jst=1,nstm
          DO ji=1,kjpindex
             IF ( soiltile(ji,jst) .LT. min_sechiba ) THEN
                water2infilt(ji,jst) = zero
                mc(ji,:,jst) = zero
             ENDIF
          ENDDO
       ENDDO

    ENDIF ! soil_upd

    !! 4. Update tmc and humtot
    
    DO jst=1,nstm
       DO ji=1,kjpindex
             tmc(ji,jst) = dz(2) * ( trois*mc(ji,1,jst) + mc(ji,2,jst) )/huit
             DO jsl = 2,nslm-1
                tmc(ji,jst) = tmc(ji,jst) + dz(jsl) * (trois*mc(ji,jsl,jst)+mc(ji,jsl-1,jst))/huit &
                     + dz(jsl+1) * (trois*mc(ji,jsl,jst)+mc(ji,jsl+1,jst))/huit
             ENDDO
             tmc(ji,jst) = tmc(ji,jst) + dz(nslm) * (trois*mc(ji,nslm,jst) + mc(ji,nslm-1,jst))/huit
             tmc(ji,jst) = tmc(ji,jst) + water2infilt(ji,jst)
             ! WARNING tmc is increased by includes water2infilt(ji,jst)
       ENDDO
    ENDDO

    humtot(:) = zero
    DO jst=1,nstm
       DO ji=1,kjpindex
          humtot(ji) = humtot(ji) + vegtot(ji) * soiltile(ji,jst) * tmc(ji,jst) ! average over grid-cell (i.e. total land)
       ENDDO
    ENDDO


    !! Now that the work is done, update resdist
    resdist(:,:) = soiltile(:,:)

    IF (printlev>=3) WRITE (numout,*) ' hydrol_tmc_update done '

  END SUBROUTINE hydrol_tmc_update

!! ================================================================================================================================
!! SUBROUTINE 	: hydrol_var_init
!!
!>\BRIEF        This routine initializes hydrologic parameters to define K and D, and diagnostic hydrologic variables.  
!!
!! DESCRIPTION  :
!! - 1 compute the depths
!! - 2 compute the profile for roots
!! - 3 compute the profile for a and n Van Genuchten parameter
!! - 4 compute the linearized values of k, a, b and d for the resolution of Fokker Planck equation
!! - 5 water reservoirs initialisation
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_var_init

  SUBROUTINE hydrol_var_init (ks, nvan, avan, mcr, mcs, mcfc, mcw, & 
       kjpindex, veget, veget_max, soiltile, njsc, altmax, &
       mx_eau_var, shumdiag_perma, &
       drysoil_frac, qsintveg, mc_layh, mcl_layh, mc_layh_s, mcl_layh_s, &
       shumdiag_peat,shumdiag_croppeat, shumdiag_man) 

    ! interface description

    !! 0. Variable and parameter declaration

    ! input scalar 
    INTEGER(i_std), INTENT(in)                          :: kjpindex      !! Domain size (number of grid cells) (1)
    ! input fields
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)   :: veget_max     !! PFT fractions within grid-cells (1; 1)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)   :: veget         !! Effective fraction of vegetation by PFT (1; 1)
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)    :: njsc          !! Index of the dominant soil textural class 
                                                                         !! in the grid cell (1-nscm, unitless) 
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT (in) :: soiltile      !! Fraction of each soil tile within vegtot (0-1, unitless)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)   :: altmax        !! Maximul active layer thickness (m). Be careful, here active means non frozen.
                                                                         !! Not related with the active soil carbon pool.
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: ks               !! Hydraulic conductivity at saturation (mm {-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: nvan             !! Van Genuchten coeficients n (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: avan             !! Van Genuchten coeficients a (mm-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: mcr              !! Residual volumetric water content (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: mcs              !! Saturated volumetric water content (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: mcfc             !! Volumetric water content at field capacity (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: mcw              !! Volumetric water content at wilting point (m^{3} m^{-3})

    !! 0.2 Output variables

    REAL(r_std),DIMENSION (kjpindex), INTENT (out)      :: mx_eau_var    !! Maximum water content of the soil 
                                                                         !! @tex $(kg m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (out) :: shumdiag_perma!! Percent of porosity filled with water (mc/mcs)
                                                                         !! used for the thermal computations
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)    :: drysoil_frac  !! function of litter humidity
    REAL(r_std), DIMENSION (kjpindex,nslm), INTENT (out):: mc_layh       !! Volumetric soil moisture content for each layer in hydrol(liquid+ice) [m3/m3]
    REAL(r_std), DIMENSION(kjpindex,nslm,nstm), INTENT (out):: mc_layh_s !! Volumetric soil moisture content for each layer in hydrol(liquid+ice) [m3/m3]
    REAL(r_std), DIMENSION(kjpindex,nslm), INTENT (out):: mcl_layh       !! Volumetric soil moisture content for each layer in hydrol(liquid) [m3/m3]
    REAL(r_std), DIMENSION(kjpindex,nslm,nstm), INTENT (out):: mcl_layh_s!! Volumetric soil moisture content for each layer in hydrol(liquid) [m3/m3]

    ! peatland
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (out) :: shumdiag_peat !! soil moisture content for the decomposition when peat is activated
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (out) :: shumdiag_croppeat !! soil moisture content for the decomposition when agri_peat is activated
    ! tides
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (out) :: shumdiag_man  !! !! soil moisture content for the decomposition when tides is activated

    !! 0.3 Modified variables
    REAL(r_std),DIMENSION(kjpindex,nvm), INTENT (inout)    :: qsintveg   !! Water on vegetation due to interception
                                                                         !! @tex $(kg m^{-2})$ @endtex
    !! peatland
    REAL(r_std),DIMENSION(nslm)                         :: afact_peat    !! Multiplicative factor for decay of a with depth, for peat
    REAL(r_std),DIMENSION(nslm)                         :: nfact_peat    !! Multiplicative factor for decay of n with depth, for peat
    REAL(r_std)                                         :: m_peat        !! m=1-1/n (unitless) for peat
    REAL(r_std)                                         :: frac_peat     !! Relative linearized VWC (unitless for peat
    REAL(r_std)                                         :: avan_mod_peat !!  VG parameter a  modified from  exponantial profile, for peat
    REAL(r_std)                                         :: nvan_mod_peat !! VG parameter n  modified from  exponantial profile, for peat
    
    !! 0.4 Local variables
    INTEGER(i_std)                                      :: ji, jv        !! Grid-cell and PFT indices (1)
    INTEGER(i_std)                                      :: jst, jsc, jsl !! Soiltile, Soil Texture, and Soil layer indices (1)
    INTEGER(i_std)                                      :: i             !! Index (1)
    REAL(r_std)                                         :: m             !! m=1-1/n (unitless)
    REAL(r_std)                                         :: frac          !! Relative linearized VWC (unitless)
    REAL(r_std)                                         :: avan_mod      !! VG parameter a modified from  exponantial profile
                                                                         !! @tex $(mm^{-1})$ @endtex
    REAL(r_std)                                         :: nvan_mod      !! VG parameter n  modified from  exponantial profile
                                                                         !! (unitless)
    REAL(r_std), DIMENSION(nslm,kjpindex)               :: afact, nfact  !! Multiplicative factor for decay of a and n with depth
                                                                         !! (unitless)
    ! parameters for "soil densification" with depth
    REAL(r_std)                                         :: dp_comp       !! Depth at which the 'compacted' value of ksat
                                                                         !! is reached (m)
    REAL(r_std)                                         :: f_ks          !! Exponential factor for decay of ksat with depth 
                                                                         !! @tex $(m^{-1})$ @endtex 
    REAL(r_std)                                         :: a0            !! fitted value for relation log((a-a0)/(a_ref-a0)) = 
                                                                         !! ak_rel * log(k/k_ref) 
                                                                         !! @tex $(mm^{-1})$ @endtex
    REAL(r_std)                                         :: ak_rel        !! fitted value for relation log((a-a0)/(a_ref-a0)) = 
                                                                         !! ak_rel * log(k/k_ref)
                                                                         !! (unitless) 
    REAL(r_std)                                         :: kfact_max     !! Maximum factor for Ks decay with depth (unitless)
    REAL(r_std)                                         :: k_tmp, tmc_litter_ratio
    INTEGER(i_std), PARAMETER                           :: error_level = 3 !! Error level for consistency check
                                                                           !! Switch to 2 tu turn fatal errors into warnings
    REAL(r_std), DIMENSION (kjpindex,nslm)              :: alphavg         !! VG param a modified with depth at each node
                                                                           !! @tex $(mm^{-1})$ @endtexe
    REAL(r_std), DIMENSION (kjpindex,nslm)              :: nvg             !! VG param n modified with depth at each node
                                                                           !! (unitless)
                                                                           !! need special treatment
    INTEGER(i_std)                                      :: ii
    INTEGER(i_std)                                      :: iiref           !! To identify the mc_lins where k_lin and d_lin
                                                                           !! need special treatment
    REAL(r_std)                                         :: root_profile_tmp
    REAL(r_std)                                         :: cum_dh          !! Depth to bottom layer (m)
    INTEGER(i_std)                                      :: nslm_root_tmp   !! Temporal, deeper root zone soil layer 
    REAL(r_std), DIMENSION (kjpindex,nslm)              :: smf             !! Soil moisture of each layer at field capacity
                                                                           !!  @tex $(kg m^{-2})$ @endtex
    
    INTEGER(i_std)                                      :: jiref           !! To identify the mc_lins where k_lin and d_lin 
                                                                           !! need special treatment
!_ ================================================================================================================================

    !Config Key   = CWRR_NKS_N0 
    !Config Desc  = fitted value for relation log((n-n0)/(n_ref-n0)) = nk_rel * log(k/k_ref)
    !Config Def   = 0.0
    !Config If    = 
    !Config Help  =
    !Config Units = [-]
    n0 = 0.0
    CALL getin_p("CWRR_NKS_N0",n0)

    !! Check parameter value (correct range)
    IF ( n0 < zero ) THEN
       CALL ipslerr_p(error_level, "hydrol_var_init.", &
            &     "Wrong parameter value for CWRR_NKS_N0.", &
            &     "This parameter should be non-negative. ", &
            &     "Please, check parameter value in run.def. ")
    END IF


    !Config Key   = CWRR_NKS_POWER
    !Config Desc  = fitted value for relation log((n-n0)/(n_ref-n0)) = nk_rel * log(k/k_ref)
    !Config Def   = 0.0
    !Config If    = 
    !Config Help  =
    !Config Units = [-]
    nk_rel = 0.0
    CALL getin_p("CWRR_NKS_POWER",nk_rel)

    !! Check parameter value (correct range)
    IF ( nk_rel < zero ) THEN
       CALL ipslerr_p(error_level, "hydrol_var_init.", &
            &     "Wrong parameter value for CWRR_NKS_POWER.", &
            &     "This parameter should be non-negative. ", &
            &     "Please, check parameter value in run.def. ")
    END IF


    !Config Key   = CWRR_AKS_A0 
    !Config Desc  = fitted value for relation log((a-a0)/(a_ref-a0)) = ak_rel * log(k/k_ref)
    !Config Def   = 0.0
    !Config If    = 
    !Config Help  =
    !Config Units = [1/mm]
    a0 = 0.0
    CALL getin_p("CWRR_AKS_A0",a0)

    !! Check parameter value (correct range)
    IF ( a0 < zero ) THEN
       CALL ipslerr_p(error_level, "hydrol_var_init.", &
            &     "Wrong parameter value for CWRR_AKS_A0.", &
            &     "This parameter should be non-negative. ", &
            &     "Please, check parameter value in run.def. ")
    END IF


    !Config Key   = CWRR_AKS_POWER
    !Config Desc  = fitted value for relation log((a-a0)/(a_ref-a0)) = ak_rel * log(k/k_ref)
    !Config Def   = 0.0
    !Config If    = 
    !Config Help  =
    !Config Units = [-]
    ak_rel = 0.0
    CALL getin_p("CWRR_AKS_POWER",ak_rel)

    !! Check parameter value (correct range)
    IF ( nk_rel < zero ) THEN
       CALL ipslerr_p(error_level, "hydrol_var_init.", &
            &     "Wrong parameter value for CWRR_AKS_POWER.", &
            &     "This parameter should be non-negative. ", &
            &     "Please, check parameter value in run.def. ")
    END IF


    !Config Key   = KFACT_DECAY_RATE
    !Config Desc  = Factor for Ks decay with depth
    !Config Def   = 2.0
    !Config If    = 
    !Config Help  =  
    !Config Units = [1/m]
    f_ks = 2.0
    CALL getin_p ("KFACT_DECAY_RATE", f_ks)

    !! Check parameter value (correct range)
    IF ( f_ks < zero ) THEN
       CALL ipslerr_p(error_level, "hydrol_var_init.", &
            &     "Wrong parameter value for KFACT_DECAY_RATE.", &
            &     "This parameter should be positive. ", &
            &     "Please, check parameter value in run.def. ")
    END IF


    !Config Key   = KFACT_STARTING_DEPTH
    !Config Desc  = Depth for compacted value of Ks 
    !Config Def   = 0.3
    !Config If    = 
    !Config Help  =  
    !Config Units = [m]
    dp_comp = 0.3
    CALL getin_p ("KFACT_STARTING_DEPTH", dp_comp)

    !! Check parameter value (correct range)
    IF ( dp_comp <= zero ) THEN
       CALL ipslerr_p(error_level, "hydrol_var_init.", &
            &     "Wrong parameter value for KFACT_STARTING_DEPTH.", &
            &     "This parameter should be positive. ", &
            &     "Please, check parameter value in run.def. ")
    END IF


    !Config Key   = KFACT_MAX
    !Config Desc  = Maximum Factor for Ks increase due to vegetation
    !Config Def   = 10.0
    !Config If    = 
    !Config Help  =
    !Config Units = [-]
    kfact_max = 10.0
    CALL getin_p ("KFACT_MAX", kfact_max)

    !! Check parameter value (correct range)
    IF ( kfact_max < 10. ) THEN
       CALL ipslerr_p(error_level, "hydrol_var_init.", &
            &     "Wrong parameter value for KFACT_MAX.", &
            &     "This parameter should be greater than 10. ", &
            &     "Please, check parameter value in run.def. ")
    END IF



    !Config Key   = KFACT_ROOT_CONST
    !Config Desc  = Set constant kfact_root in every soil layer. Otherwise kfact_root increase over soil depth in the rootzone. 
    !Config If    = 
    !Config Def   = n
    !Config Help  = Use KFACT_ROOT_CONST=true to impose kfact_root=1 in every soil layer. Otherwise kfact_root increase over soil depth in the rootzone. 
    !Config Units = [y/n]
    kfact_root_const = .FALSE.
    CALL getin_p("KFACT_ROOT_CONST",kfact_root_const)

    
    ! Set is_irrigated logical for each PFT
    DO jv=1,nvm
       IF (pref_soil_veg(jv).eq.irrig_st) THEN
          is_irrigated(jv) = .TRUE.
       ELSE ! PFT is not from irrigated soil tile
          is_irrigated(jv) = .FALSE.
       ENDIF   
    ENDDO
    
    !-
    !! 1 Create local variables in mm for the vertical depths
    !!   Vertical depth variables (znh, dnh, dlh) are stored in module vertical_soil_var in m.
    DO jsl=1,nslm
       zz(jsl) = znh(jsl)*kilo_to_unit
       dz(jsl) = dnh(jsl)*kilo_to_unit
       dh(jsl) = dlh(jsl)*kilo_to_unit
    ENDDO

    !! 2 Compute the maximum layers defining the root zone
    cum_dh = zero
    nslm_root(:) =  nslm
    nslm_root_tmp =  nslm
    DO jsl =1, nslm 
       IF (( cum_dh ) < cum_dh_thr*mille) THEN  
          cum_dh = cum_dh + dh(jsl) 
          nslm_root_tmp = jsl 
       ENDIF
    END DO

    DO ji=1,kjpindex 
       IF (SUM(veget_max(ji, : ), MASK= is_irrigated(:) ) > min_sechiba) THEN 
          nslm_root(ji) = nslm_root_tmp 
       ENDIF
    ENDDO

    ! Calculates field capacity soil moisture per soil layers
    ! then calculate field capacity soil moisture over root zone
    smf(:,:) = zero
    root_mc_fc(:) = zero

    smf(:,1) = dz(2) * (quatre*mcfc(:))/huit 
    DO jsl = 2,nslm-1
       smf(:,jsl) = dz(jsl) * ( quatre*mcfc(:) )/huit &
            + dz(jsl+1) * ( quatre*mcfc(:) )/huit
    ENDDO
    smf(:,nslm) = dz(nslm) * (quatre*mcfc(:))/huit
    
    DO ji = 1,kjpindex
       root_mc_fc(ji) = SUM(smf(ji,1:nslm_root(ji) ))
    ENDDO
    
    !! 3 Compute the profile for a and n

    ! For every soil texture
    DO ji = 1, kjpindex 
       DO jsl=1,nslm
          ! PhD thesis of d'Orgeval, 2006, p81, Eq. 4.38; d'Orgeval et al. 2008, Eq. 2
          ! Calibrated against Hapex-Sahel measurements
          kfact(jsl,ji) = MIN(MAX(EXP(- f_ks * (zz(jsl)/mille - dp_comp)), un/kfact_max),un)
          ! PhD thesis of d'Orgeval, 2006, p81, Eqs. 4.39; 4.42, and Fig 4.14

          nfact(jsl,ji) = ( kfact(jsl,ji) )**nk_rel
          afact(jsl,ji) = ( kfact(jsl,ji) )**ak_rel
       ENDDO
    ENDDO

    ! peatland
    ! kfact_peat is used to compue d_lin_peat and then to compue a, b, k
    ! soiltile condition is added later to deal with pft peat is 0, by using kfact_peat
    ! kfact_peat is defined the same as kfact
    IF (ok_peat_hydro) THEN
       DO jsl=1,nslm
          kfact_peat(jsl) = MIN(MAX(EXP(- f_ks * (zz(jsl)/mille - dp_comp)),un/kfact_max),un)
          nfact_peat(jsl)=( kfact_peat(jsl) )**nk_rel
          afact_peat(jsl) = ( kfact_peat(jsl) )**ak_rel
       ENDDO
    ENDIF

    ! For every grid cell
     DO ji = 1, kjpindex
       !-
       !! 4 Compute the linearized values of k, a, b and d
       !!   The effect of kfact_root on ks thus on k, a, n and d, is taken into account further in the code,
       !!   in hydrol_soil_coef.
       !-
       ! Calculate the matrix coef for Dublin model (de Rosnay, 1999; p149)
       ! piece-wise linearised hydraulic conductivity k_lin=alin * mc_lin + b_lin
       ! and diffusivity d_lin in each interval of mc, called mc_lin,
       ! between imin, for residual mcr, and imax for saturation mcs.

       ! We define 51 bounds for 50 bins of mc between mcr and mcs
       mc_lin(imin,ji)=mcr(ji)
       mc_lin(imax,ji)=mcs(ji)
       DO ii= imin+1, imax-1 ! ii=2,50
          mc_lin(ii,ji) = mcr(ji) + (ii-imin)*(mcs(ji)-mcr(ji))/(imax-imin)
       ENDDO

       DO jsl = 1, nslm
          ! From PhD thesis of d'Orgeval, 2006, p81, Eq. 4.42
          nvan_mod = n0 + (nvan(ji)-n0) * nfact(jsl,ji)
          avan_mod = a0 + (avan(ji)-a0) * afact(jsl,ji)
          m = un - un / nvan_mod
          ! Creation of arrays for SP-MIP output by landpoint
          nvan_mod_tab(jsl,ji) = nvan_mod
          avan_mod_tab(jsl,ji) = avan_mod
          ! We apply Van Genuchten equation for K(theta) based on Ks(z)=ks(ji) * kfact(jsl,ji)
          DO ii = imax,imin,-1
             frac=MIN(un,(mc_lin(ii,ji)-mcr(ji))/(mcs(ji)-mcr(ji)))
             k_lin(ii,jsl,ji) = ks(ji) * kfact(jsl,ji) * (frac**0.5) * ( un - ( un - frac ** (un/m)) ** m )**2
          ENDDO

          ! k_lin should not be zero, nor too small
          ! We track iiref, the bin under which mc is too small and we may get zero k_lin
          !salma: ji replaced with ii and jiref replaced with iiref and jsc with ji
          ii=imax-1
          DO WHILE ((k_lin(ii,jsl,ji) > 1.e-32) .and. (ii>0))
             iiref=ii
             ii=ii-1
          ENDDO
          DO ii=iiref-1,imin,-1
             k_lin(ii,jsl,ji)=k_lin(ii+1,jsl,ji)/10.
          ENDDO

          DO ii = imin,imax-1 ! ii=1,50
             ! We deduce a_lin and b_lin based on continuity between segments k_lin = a_lin*mc-lin+b_lin
             a_lin(ii,jsl,ji) = (k_lin(ii+1,jsl,ji)-k_lin(ii,jsl,ji)) / (mc_lin(ii+1,ji)-mc_lin(ii,ji))
             b_lin(ii,jsl,ji)  = k_lin(ii,jsl,ji) - a_lin(ii,jsl,ji)*mc_lin(ii,ji)

             ! We calculate the d_lin for each mc bin, from Van Genuchten equation for D(theta)
             ! d_lin is constant and taken as the arithmetic mean between the values at the bounds of each bin
             IF (ii.NE.imin .AND. ii.NE.imax-1) THEN
                frac=MIN(un,(mc_lin(ii,ji)-mcr(ji))/(mcs(ji)-mcr(ji)))
                d_lin(ii,jsl,ji) =(k_lin(ii,jsl,ji) / (avan_mod*m*nvan_mod)) *  &
                     ( (frac**(-un/m))/(mc_lin(ii,ji)-mcr(ji)) ) * &
                     (  frac**(-un/m) -un ) ** (-m)
                frac=MIN(un,(mc_lin(ii+1,ji)-mcr(ji))/(mcs(ji)-mcr(ji)))
                d_lin(ii+1,jsl,ji) =(k_lin(ii+1,jsl,ji) / (avan_mod*m*nvan_mod))*&
                     ( (frac**(-un/m))/(mc_lin(ii+1,ji)-mcr(ji)) ) * &
                     (  frac**(-un/m) -un ) ** (-m)
                d_lin(ii,jsl,ji) = undemi * (d_lin(ii,jsl,ji)+d_lin(ii+1,jsl,ji))
             ELSE IF(ii.EQ.imax-1) THEN
                d_lin(ii,jsl,ji) =(k_lin(ii,jsl,ji) / (avan_mod*m*nvan_mod)) * &
                     ( (frac**(-un/m))/(mc_lin(ii,ji)-mcr(ji)) ) *  &
                     (  frac**(-un/m) -un ) ** (-m)
             ENDIF
          ENDDO !Salma end loop over landpoints

          ! Special case for ii=imin
          d_lin(imin,jsl,ji) = d_lin(imin+1,jsl,ji)/1000.

          ! We adjust d_lin where k_lin was previously adjusted otherwise we might get non-monotonous variations
          ! We don't want d_lin = zero
          DO ii=iiref-1,imin,-1
             d_lin(ii,jsl,ji)=d_lin(ii+1,jsl,ji)/10.
          ENDDO

       ENDDO
    ENDDO

    ! Output of alphavg and nvg at each node for SP-MIP
    ! note here alphavg and nvg are only for nonpead land grids
    ! if we want to include peatland in these output variables, 
    ! then we need to use soiltile to combine avant_mod_tab and avant_mod_tab_peat
    DO jsl = 1, nslm
       alphavg(:,jsl) = avan_mod_tab(jsl,:)*1000. ! from mm-1 to m-1
       nvg(:,jsl) = nvan_mod_tab(jsl,:)
    ENDDO
    CALL xios_orchidee_send_field("alphavg",alphavg) ! in m-1
    CALL xios_orchidee_send_field("nvg",nvg) ! unitless

    
    ! peatland
    ! We define 51 bounds for 50 bins of mc between mcr and mcs
    ! no need to modify soiltile, because some change is applied in later code which uses these variables
    IF (ok_peat_hydro) THEN
       mc_lin_peat(imin)=mcr_peat
       mc_lin_peat(imax)=mcs_peat
       DO ji= imin+1, imax-1 
          mc_lin_peat(ji) = mcr_peat + (ji-imin)*(mcs_peat-mcr_peat)/(imax-imin)
       ENDDO
       
       DO jsl = 1, nslm
          nvan_mod_peat = n0 + (nvan_peat-n0) * nfact_peat(jsl)
          avan_mod_peat = a0 + (avan_peat-a0) * afact_peat(jsl)
          m_peat = un - un / nvan_mod_peat
          !
          ! We apply Van Genuchten equation for K(theta) based on Ks(z)=ks(jsc) * kfact(jsl,jsc)
          DO ji = imax,imin,-1
             frac_peat=MIN(un,(mc_lin_peat(ji)-mcr_peat)/(mcs_peat-mcr_peat))
             k_lin_peat(ji,jsl) = ks_peat * kfact_peat(jsl) * (frac_peat**0.5) * ( un - ( un - frac_peat ** (un/m_peat)) ** m_peat )**2
          ENDDO
          ! k_lin should not be zero, nor too small
          ! We track jiref, the bin under which mc is too small and we may get zero k_lin     
          ji=imax-1
          DO WHILE ((k_lin_peat(ji,jsl) > 1.e-32) .and. (ji>0))
             jiref=ji
             ji=ji-1
          ENDDO
          DO ji=jiref-1,imin,-1
             k_lin_peat(ji,jsl)=k_lin_peat(ji+1,jsl)/10.
          ENDDO
          DO ji = imin,imax-1
             ! We deduce a_lin and b_lin based on continuity between segments k_lin = a_lin*mc-lin+b_lin
             a_lin_peat(ji,jsl) = (k_lin_peat(ji+1,jsl)-k_lin_peat(ji,jsl)) / (mc_lin_peat(ji+1)-mc_lin_peat(ji))
             b_lin_peat(ji,jsl)  = k_lin_peat(ji,jsl) - a_lin_peat(ji,jsl)*mc_lin_peat(ji)
             
             ! We calculate the d_lin for each mc bin, from Van Genuchten equation for D(theta)
             ! d_lin is constant and taken as the arithmetic mean between the values at the bounds of each bin 
             IF (ji.NE.imin .AND. ji.NE.imax-1) THEN
                frac_peat=MIN(un,(mc_lin_peat(ji)-mcr_peat)/(mcs_peat-mcr_peat))
                d_lin_peat(ji,jsl) =(k_lin_peat(ji,jsl) / (avan_mod_peat*m_peat*nvan_mod_peat)) *  &
                     ( (frac_peat**(-un/m_peat))/(mc_lin_peat(ji)-mcr_peat) ) * &
                     (  frac_peat**(-un/m_peat) -un ) ** (-m_peat)
                frac_peat=MIN(un,(mc_lin_peat(ji+1)-mcr_peat)/(mcs_peat-mcr_peat))
                d_lin_peat(ji+1,jsl) =(k_lin_peat(ji+1,jsl) / (avan_mod_peat*m_peat*nvan_mod_peat))*&
                     ( (frac_peat**(-un/m_peat))/(mc_lin_peat(ji+1)-mcr_peat) ) * &
                     (  frac_peat**(-un/m_peat) -un ) ** (-m_peat)
                d_lin_peat(ji,jsl) = undemi * (d_lin_peat(ji,jsl)+d_lin_peat(ji+1,jsl))
             ELSE IF(ji.EQ.imax-1) THEN
                d_lin_peat(ji,jsl) =(k_lin_peat(ji,jsl) / (avan_mod_peat*m_peat*nvan_mod_peat)) * &
                     ( (frac_peat**(-un/m_peat))/(mc_lin_peat(ji)-mcr_peat) ) *  &
                     (  frac_peat**(-un/m_peat) -un ) ** (-m_peat)
             ENDIF
          ENDDO
          
          ! Special case for ji=imin
          d_lin_peat(imin,jsl) = d_lin_peat(imin+1,jsl)/1000.
          
          ! We adjust d_lin where k_lin was previously adjusted otherwise we might get non-monotonous variations
          ! We don't want d_lin = zero
          DO ji=jiref-1,imin,-1
             d_lin_peat(ji,jsl)=d_lin_peat(ji+1,jsl)/10.
          ENDDO
          
       ENDDO
    ENDIF

    !! 5 Water reservoir initialisation
    !
!!$    DO jst = 1,nstm
!!$       DO ji = 1, kjpindex
!!$          mx_eau_var(ji) = mx_eau_var(ji) + soiltile(ji,jst)*&
!!$               &   zmaxh*mille*mcs(njsc(ji))
!!$       END DO
!!$    END DO

    mx_eau_var(:) = zero

    ! peatland 
    IF (ok_peat_hydro) THEN
       IF (agri_peat) THEN
          mx_eau_var(:)=zmaxh*mille*((mcs(:)*(un-soiltile(:,4)-soiltile(:,5)-soiltile(:,6)))+mcs_peat*(soiltile(:,4)+soiltile(:,5)+soiltile(:,6)))
       ELSE
          IF (tides) THEN 
             mx_eau_var(:)=zmaxh*mille*((mcs(:)*(un-soiltile(:,4)-soiltile(:,6)))+mcs_peat*(soiltile(:,4)+soiltile(:,6)))
          ELSE
             ! no need to modify for the case of pft peat is 0, because soiltile4 is 0 in this equation
             mx_eau_var(:)=zmaxh*mille*((mcs(:)*(un-soiltile(:,4)))+mcs_peat*soiltile(:,4))
          ENDIF
       ENDIF
    ELSE
       mx_eau_var(:) = zmaxh*mille*mcs(:) 
    ENDIF

    ! Compute the litter humidity, shumdiag and fry
    shumdiag_perma(:,:) = zero
    humtot(:) = zero
    tmc(:,:) = zero

    shumdiag_peat(:,:) = zero
    shumdiag_croppeat(:,:) = zero
    shumdiag_man(:,:) = zero

    ! Loop on soiltiles to compute the variables (ji,jst)
    DO jst=1,nstm 
       DO ji = 1, kjpindex
          ! peatland
          ! tmcs tmcr are useful to compute swc_pft, land_tmcs etc
          IF ( (ok_peat_hydro) .AND. (jst .GE. 4) .AND. (jst .LE. 6)) Then           
             If (soiltile(ji,jst) .GT. min_sechiba) THEN
                ! here the variables are computed per soil tile, ie, peat here
                tmcs(ji,jst)=zmaxh* mille*mcs_peat
                tmcr(ji,jst)=zmaxh* mille*mcr_peat
                tmcfc(ji,jst)=zmaxh* mille*mcf_peat
                tmcw(ji,jst)=zmaxh* mille*mcw_peat
             ELSE
                tmcs(ji,jst) = 0.0
                tmcr(ji,jst) = 0.0
                tmcs(ji,jst) = 0.0
                tmcfc(ji,jst) = 0.0
                tmcw(ji,jst) = 0.0
             ENDIF
          ELSE
             tmcs(ji,jst)=zmaxh* mille*mcs(ji)
             tmcr(ji,jst)=zmaxh* mille*mcr(ji)
             tmcfc(ji,jst)=zmaxh* mille*mcfc(ji)
             tmcw(ji,jst)=zmaxh* mille*mcw(ji)
          ENDIF
       ENDDO
    ENDDO
       
    ! The total soil moisture for each soiltile:
    DO jst=1,nstm
       DO ji=1,kjpindex
          tmc(ji,jst)= dz(2) * ( trois*mc(ji,1,jst)+ mc(ji,2,jst))/huit
       END DO
    ENDDO

    DO jst=1,nstm 
       DO jsl=2,nslm-1
          DO ji=1,kjpindex
             tmc(ji,jst) = tmc(ji,jst) + dz(jsl) * ( trois*mc(ji,jsl,jst) + mc(ji,jsl-1,jst))/huit &
                  & + dz(jsl+1)*(trois*mc(ji,jsl,jst) + mc(ji,jsl+1,jst))/huit
          END DO
       END DO
    ENDDO

    DO jst=1,nstm 
       DO ji=1,kjpindex
          tmc(ji,jst) = tmc(ji,jst) +  dz(nslm) * (trois * mc(ji,nslm,jst) + mc(ji,nslm-1,jst))/huit
          tmc(ji,jst) = tmc(ji,jst) + water2infilt(ji,jst)
       ENDDO
    END DO

    ! Initialize humtot such that twbr is also closed at the first time step
    humtot(:) = zero
    DO jst=1,nstm
       DO ji=1,kjpindex
          !average over grid-cell (i.e. total land)
          humtot(ji) = humtot(ji) + vegtot(ji) * resdist(ji,jst) * tmc(ji,jst)
       ENDDO
    ENDDO

!JG: hydrol_tmc_update should not be called in the initialization phase. Call of hydrol_tmc_update makes the model restart differenlty.    
!    ! If veget has been updated before restart (with LAND USE or DGVM),
!    ! tmc and mc must be modified with respect to humtot conservation.
!   CALL hydrol_tmc_update ( kjpindex, veget_max, soiltile, qsintveg)

    ! The litter variables:
    ! level 1
    DO jst=1,nstm 
       DO ji=1,kjpindex
          IF ( (ok_peat_hydro) .AND. (jst .GE. 4) .AND. (jst .LE. 6)) Then
             IF (soiltile(ji,jst) .GT. min_sechiba) THEN
                tmc_litter(ji,jst) = dz(2) * (trois*mcl(ji,1,jst)+mcl(ji,2,jst))/huit
                tmc_litter_wilt(ji,jst) = dz(2) * mcw_peat / deux
                tmc_litter_res(ji,jst) = dz(2) * mcr_peat / deux
                tmc_litter_field(ji,jst) = dz(2) * mcf_peat / deux
                tmc_litter_sat(ji,jst) = dz(2) * mcs_peat / deux
                tmc_litter_awet(ji,jst) = dz(2) * mc_awet_peat / deux
                tmc_litter_adry(ji,jst) = dz(2) * mc_adry_peat / deux
             ELSE
                tmc_litter(ji,jst) = 0.0
                tmc_litter_wilt(ji,jst) = 0.0
                tmc_litter_res(ji,jst) = 0.0
                tmc_litter_field(ji,jst) = 0.0
                tmc_litter_sat(ji,jst) = 0.0
                tmc_litter_awet(ji,jst) = 0.0
                tmc_litter_adry(ji,jst) = 0.0
             ENDIF
          ELSE
             tmc_litter(ji,jst) = dz(2) * (trois*mcl(ji,1,jst)+mcl(ji,2,jst))/huit
             tmc_litter_wilt(ji,jst) = dz(2) * mcw(ji) / deux
             tmc_litter_res(ji,jst) = dz(2) * mcr(ji) / deux
             tmc_litter_field(ji,jst) = dz(2) * mcfc(ji) / deux
             tmc_litter_sat(ji,jst) = dz(2) * mcs(ji) / deux
             tmc_litter_awet(ji,jst) = dz(2) * mc_awet(njsc(ji)) / deux
             tmc_litter_adry(ji,jst) = dz(2) * mc_adry(njsc(ji)) / deux
          ENDIF
       ENDDO
    END DO
    ! sum from level 2 to 4
    DO jst=1,nstm 
       DO jsl=2,4
          DO ji=1,kjpindex
             IF ( (ok_peat_hydro) .AND. (jst .GE. 4) .AND. (jst .LE. 6)) THEN
                IF (soiltile(ji,jst) .GT. min_sechiba) THEN
                   tmc_litter(ji,jst) = tmc_litter(ji,jst) + dz(jsl) * &
                        & ( trois*mcl(ji,jsl,jst) + mcl(ji,jsl-1,jst))/huit &
                        & + dz(jsl+1)*(trois*mcl(ji,jsl,jst) + mcl(ji,jsl+1,jst))/huit
                   tmc_litter_wilt(ji,jst) = tmc_litter_wilt(ji,jst) + &
                        &(dz(jsl)+ dz(jsl+1))* mcw_peat/deux
                   tmc_litter_res(ji,jst) = tmc_litter_res(ji,jst) + &
                        &(dz(jsl)+ dz(jsl+1))* mcr_peat/deux
                   tmc_litter_sat(ji,jst) = tmc_litter_sat(ji,jst) + &
                        &(dz(jsl)+ dz(jsl+1))*  mcs_peat/deux
                   tmc_litter_field(ji,jst) = tmc_litter_field(ji,jst) + &
                        & (dz(jsl)+ dz(jsl+1))* mcf_peat/deux
                   tmc_litter_awet(ji,jst) = tmc_litter_awet(ji,jst) + &
                        &(dz(jsl)+ dz(jsl+1))* &
                        & mc_awet_peat/deux
                   tmc_litter_adry(ji,jst) = tmc_litter_adry(ji,jst) + &
                        & (dz(jsl)+ dz(jsl+1))* &
                        & mc_adry_peat/deux  
                   ELSE
                      tmc_litter(ji,jst) = 0.0
                      tmc_litter_wilt(ji,jst) = 0.0
                      tmc_litter_res(ji,jst) = 0.0
                      tmc_litter_field(ji,jst) = 0.0
                      tmc_litter_sat(ji,jst) = 0.0
                      tmc_litter_awet(ji,jst) = 0.0
                      tmc_litter_adry(ji,jst) = 0.0
                ENDIF
             ELSE
                tmc_litter(ji,jst) = tmc_litter(ji,jst) + dz(jsl) * &
                     & ( trois*mcl(ji,jsl,jst) + mcl(ji,jsl-1,jst))/huit &
                     & + dz(jsl+1)*(trois*mcl(ji,jsl,jst) + mcl(ji,jsl+1,jst))/huit
                tmc_litter_wilt(ji,jst) = tmc_litter_wilt(ji,jst) + &
                     &(dz(jsl)+ dz(jsl+1))*&
                     & mcw(ji)/deux
                tmc_litter_res(ji,jst) = tmc_litter_res(ji,jst) + &
                     &(dz(jsl)+ dz(jsl+1))*&
                     & mcr(ji)/deux
                tmc_litter_sat(ji,jst) = tmc_litter_sat(ji,jst) + &
                     &(dz(jsl)+ dz(jsl+1))* &
                     & mcs(ji)/deux
                tmc_litter_field(ji,jst) = tmc_litter_field(ji,jst) + &
                     & (dz(jsl)+ dz(jsl+1))* &
                     & mcfc(ji)/deux
                tmc_litter_awet(ji,jst) = tmc_litter_awet(ji,jst) + &
                     &(dz(jsl)+ dz(jsl+1))* &
                     & mc_awet(njsc(ji))/deux
                tmc_litter_adry(ji,jst) = tmc_litter_adry(ji,jst) + &
                     & (dz(jsl)+ dz(jsl+1))* &
                     & mc_adry(njsc(ji))/deux
             ENDIF
          END DO
       END DO
    END DO


    DO jst=1,nstm 
       DO ji=1,kjpindex
          ! here we set that humrelv=0 in PFT1
          humrelv(ji,1,jst) = zero
       ENDDO
    END DO


    ! Calculate shumdiag_perma for thermosoil
    ! Use resdist instead of soiltile because we here need to have 
    ! shumdiag_perma at the value from previous time step.
    ! Here, soilmoist is only used as a temporary variable to calculate shumdiag_perma
    ! (based on resdist=soiltile from previous timestep, but normally equal to soiltile)
    ! For consistency with hydrol_soil, we want to calculate a grid-cell average
    soilmoist(:,:) = zero
    DO jst=1,nstm
       DO ji=1,kjpindex
          soilmoist(ji,1) = soilmoist(ji,1) + resdist(ji,jst) * &
               dz(2) * ( trois*mc(ji,1,jst) + mc(ji,2,jst) )/huit
          DO jsl = 2,nslm-1
             soilmoist(ji,jsl) = soilmoist(ji,jsl) + resdist(ji,jst) * &
                  ( dz(jsl) * (trois*mc(ji,jsl,jst)+mc(ji,jsl-1,jst))/huit &
                  + dz(jsl+1) * (trois*mc(ji,jsl,jst)+mc(ji,jsl+1,jst))/huit )
          END DO
          soilmoist(ji,nslm) = soilmoist(ji,nslm) + resdist(ji,jst) * &
               dz(nslm) * (trois*mc(ji,nslm,jst) + mc(ji,nslm-1,jst))/huit
       ENDDO
    ENDDO
    DO ji=1,kjpindex
        soilmoist(ji,:) = soilmoist(ji,:) * vegtot_old(ji) ! grid cell average 
    ENDDO

    soilmoist_s(:,:,:) = zero
    DO jst=1,nstm
       DO ji=1,kjpindex
          soilmoist_s(ji,1,jst) = soilmoist_s(ji,1,jst) + resdist(ji,jst) * &
               dz(2) * ( trois*mc(ji,1,jst) + mc(ji,2,jst) )/huit
          DO jsl = 2,nslm-1
             soilmoist_s(ji,jsl,jst) = soilmoist_s(ji,jsl,jst) + resdist(ji,jst) * &
                  ( dz(jsl) * (trois*mc(ji,jsl,jst)+mc(ji,jsl-1,jst))/huit &
                  + dz(jsl+1) * (trois*mc(ji,jsl,jst)+mc(ji,jsl+1,jst))/huit )
          END DO
          soilmoist_s(ji,nslm,jst) = soilmoist_s(ji,nslm,jst) + resdist(ji,jst) * &
               dz(nslm) * (trois*mc(ji,nslm,jst) + mc(ji,nslm-1,jst))/huit
       ENDDO
    ENDDO
    DO ji=1,kjpindex
        soilmoist_s(ji,:,:) = soilmoist_s(ji,:,:) * vegtot_old(ji) ! grid cell average 
    ENDDO
    
    ! peatland
    IF (ok_peat_hydro) THEN
       shumdiag_peat(:,:)=zero
       DO ji=1,kjpindex
          IF ( soiltile(ji,4) .GT. min_sechiba ) THEN
             shumdiag_peat(ji,1)= dz(2) * ( trois*mc(ji,1,4) + mc(ji,2,4) )/huit 
             DO jsl = 2,nslm-1
                shumdiag_peat(ji,jsl)=dz(jsl) * (trois*mc(ji,jsl,4)+mc(ji,jsl-1,4))/huit &
                     + dz(jsl+1) * (trois*mc(ji,jsl,4)+mc(ji,jsl+1,4))/huit 
             ENDDO
             shumdiag_peat(ji,nslm)= dz(nslm) * (trois*mc(ji,nslm,4) + mc(ji,nslm-1,4))/huit
          ENDIF
       ENDDO
    ENDIF
    
    IF (tides) THEN
       shumdiag_man(:,:)=zero
       DO ji=1,kjpindex
          shumdiag_man(ji,1)= dz(2) * ( trois*mc(ji,1,6) + mc(ji,2,6) )/huit
          DO jsl = 2,nslm-1
             shumdiag_man(ji,jsl)=dz(jsl) *(trois*mc(ji,jsl,6)+mc(ji,jsl-1,6))/huit &
                  + dz(jsl+1) * (trois*mc(ji,jsl,6)+mc(ji,jsl+1,6))/huit
          ENDDO
          shumdiag_man(ji,nslm)= dz(nslm) * (trois*mc(ji,nslm,6) + mc(ji,nslm-1,6))/huit
       ENDDO
    ENDIF
    
    IF (ok_peat_hydro) THEN
       DO jsl=1,nslm
          DO ji=1,kjpindex
             IF ( soiltile(ji,4) .GT. min_sechiba ) THEN
                ! stomate_litter need this variable 
                shumdiag_peat(ji,jsl)=shumdiag_peat(ji,jsl)/(dh(jsl)) 
                shumdiag_peat(ji,jsl)=MAX(MIN(shumdiag_peat(ji,jsl), un), zero)
             ENDIF
          ENDDO
       ENDDO
    ENDIF
    
    IF (tides) THEN
       DO jsl=1,nslm
          DO ji=1,kjpindex
             shumdiag_man(ji,jsl)=shumdiag_man(ji,jsl)/(dh(jsl))
             shumdiag_man(ji,jsl)=MAX(MIN(shumdiag_man(ji,jsl), un), zero)
          ENDDO
       ENDDO
    ENDIF
    
    IF (agri_peat) THEN
       shumdiag_croppeat(:,:)=zero
       DO ji=1,kjpindex
          shumdiag_croppeat(ji,1)= dz(2) * ( trois*mc(ji,1,5) + mc(ji,2,5) )/huit
          DO jsl = 2,nslm-1
             shumdiag_croppeat(ji,jsl)=dz(jsl) * (trois*mc(ji,jsl,5)+mc(ji,jsl-1,5))/huit &
                  + dz(jsl+1) * (trois*mc(ji,jsl,5)+mc(ji,jsl+1,5))/huit
          ENDDO
          shumdiag_croppeat(ji,nslm)= dz(nslm) * (trois*mc(ji,nslm,5) + mc(ji,nslm-1,5))/huit
       ENDDO
    ENDIF
    IF (agri_peat) THEN
       DO jsl=1,nslm
          DO ji=1,kjpindex
             shumdiag_croppeat(ji,jsl)=shumdiag_croppeat(ji,jsl)/(dh(jsl))
             shumdiag_croppeat(ji,jsl)=MAX(MIN(shumdiag_croppeat(ji,jsl), un), zero)
          ENDDO
       ENDDO
    ENDIF

    ! -- shumdiag_perma for restart
    ! For consistency with hydrol_soil, we want to calculate a grid-cell average
    ! shumdiag_perma is initialized to 0 already; not necessary to redefine its value if peat pft is 0.
    ! it is relative humidity.
    DO jsl = 1, nslm
       DO ji=1,kjpindex  
          !dh: soil layer thickness
          IF ( (ok_peat_hydro)) THEN
             IF (soiltile (ji,4) .GT. min_sechiba) THEN  
                ! consider the average of peat and non-peat for this parameter
                shumdiag_perma(ji,jsl) = soilmoist(ji,jsl)/(dh(jsl)*(mcs(ji)*(un-soiltile(ji,4))+mcs_peat*soiltile(ji,4))) 
             ENDIF
          ELSE
             shumdiag_perma(ji,jsl) = soilmoist(ji,jsl) / (dh(jsl)*mcs(ji))
          ENDIF
          shumdiag_perma(ji,jsl) = MAX(MIN(shumdiag_perma(ji,jsl), un), zero) 
       ENDDO
    ENDDO
               
    ! Calculate drysoil_frac if it was not found in the restart file
    ! For simplicity, we set drysoil_frac to 0.5 in this case
    IF (ALL(drysoil_frac(:) == val_exp)) THEN
       DO ji=1,kjpindex
          drysoil_frac(ji) = 0.5
       END DO
    END IF

    !! Calculate the volumetric soil moisture content (mc_layh and mcl_layh) needed in 
    !! thermosoil for the thermal conductivity. 
    ! These values are only used in thermosoil_init in absence of a restart file

    mc_layh(:,:) = zero
    mcl_layh(:,:) = zero
    mc_layh_s = mc
    mcl_layh_s = mc

    DO jst=1,nstm
       DO jsl=1,nslm
          DO ji=1,kjpindex
            mc_layh(ji,jsl) = mc_layh(ji,jsl) + mc(ji,jsl,jst) * resdist(ji,jst)  * vegtot_old(ji)
            mcl_layh(ji,jsl) = mcl_layh(ji,jsl) + mcl(ji,jsl,jst) * resdist(ji,jst) * vegtot_old(ji)
         ENDDO
      END DO
    END DO

    IF (printlev>=3) WRITE (numout,*) ' hydrol_var_init done '

  END SUBROUTINE hydrol_var_init



   
!! ================================================================================================================================
!! SUBROUTINE 	: hydrol_canop
!!
!>\BRIEF        This routine computes canopy processes.
!!
!! DESCRIPTION  :
!! - 1 evaporation off the continents
!! - 1.1 The interception loss is take off the canopy. 
!! - 1.2 precip_rain is shared for each vegetation type
!! - 1.3 Limits the effect and sum what receives soil
!! - 1.4 swap qsintveg to the new value
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_canop

  SUBROUTINE hydrol_canop (kjpindex, precip_rain, vevapwet, veget_max, veget, qsintmax, &
       & qsintveg,precisol,tot_melt, frac_snow_veg)

    ! 
    ! interface description
    !

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                               :: kjpindex    !! Domain size
    ! input fields
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: precip_rain !! Rain precipitation
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(in)        :: veget_max   !! max fraction of vegetation type
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(in)        :: veget       !! Fraction of vegetation type 
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(in)        :: qsintmax    !! Maximum water on vegetation for interception
    REAL(r_std), DIMENSION  (kjpindex), INTENT (in)          :: tot_melt    !! Total melt
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(in)        :: vevapwet    !! Interception loss
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: frac_snow_veg !! Snow cover fraction on vegetation    

    !! 0.2 Output variables
    
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(out)       :: precisol    !! Water fallen onto the ground (throughfall+Totmelt)

    !! 0.3 Modified variables
    
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(inout)     :: qsintveg    !! Water on vegetation due to interception

    !! 0.4 Local variables

    INTEGER(i_std)                                           :: ji, jv
    REAL(r_std), DIMENSION (kjpindex,nvm)                    :: zqsintvegnew

!_ ================================================================================================================================

    ! boucle sur les points continentaux
    ! calcul de qsintveg au pas de temps suivant
    ! par ajout du flux interception loss
    ! calcule par enerbil en fonction
    ! des calculs faits dans diffuco
    ! calcul de ce qui tombe sur le sol
    ! avec accumulation dans precisol
    ! essayer d'harmoniser le traitement du sol nu
    ! avec celui des differents types de vegetation
    ! fait si on impose qsintmax ( ,1) = 0.0
    
    !! 1 evaporation off the continents
    
    !! 1.1 Precipitation on bare soil
    !  Precip_rain (mm) needs to be distributed over the different PFTs. Bare soil will also
    !  receive precipitation but because there is no canopy on bare soil, there is no precipitation
    !  accumulated on the leaves. 
    qsintveg(:,1) = zero
    precisol(:,1) = veget_max(:,1)*(1-frac_snow_veg(:))*precip_rain(:)
    
    !! 1.2 Interception loss
    !  Interception loss is taken off the water that is stored on the leaves of the canopy.
    !  qsintveg has been observed to take on small negative values (-10-e5 to -10e-11). This was
    !  assumed to be be a consequence of the implicit coupling (see ticket 201). The negative value 
    !  should hovere be small (not clear what small means in this context), At the next time step
    !  vbeta2 should be zero. In diffuco there are some efforts to avoid this situation but it seems
    !  that those efforts are not 100%.
    DO jv=2,nvm
       qsintveg(:,jv) = qsintveg(:,jv) - vevapwet(:,jv)
    END DO

    !! 1.3 Calculate the water stored on the leaves
    !  It is raining: precip_rain is shared over the different PFTs. Because the time step
    !  is rather long (30 minutes) it is unrealistic to assume that all the precipitation
    !  falling during the time step will be stored on the leaves. If this assumption is not made
    !  the interception loss will likely be too high. ORCHIDEE overcomes this issue by
    !  assuming that part of the preciption that intercats with the canopy will be stored
    !  on the leaves. The leaves will, however, become too heavy and tip. This tipping water will
    !  become throughfall before it can be intercepted and evaporated. This approach overcomes the
    !  need to explicitly calculate leaf tipping due to precipitation accumulation. The share of the
    !  intercepted water that will contribute to this type of throughfall is given by the parameter
    !  throughfall_by_pft.
    DO jv=2,nvm
       qsintveg(:,jv) = qsintveg(:,jv) + veget(:,jv) * ((1-throughfall_by_pft(jv))*precip_rain(:))
    END DO


    !! 1.4 Limits the effect and sum what receives soil
    !  Calculate the precipitation that is stored on the leaves (zqsintvegnew). Precipitation that
    !  passes through a canopy gap will not interact with the canopy (veget_max - veget) and will
    !  therefore contribute directly to throughfall. veget is calculated as the projected leaf area.
    !  By definition all the precipitation that falls over veget will interact with the canopy
    !  (as there are no gaps left in veget) but part of this precipitation (described by throughfall_by_pft)
    !  is moved directly to the throughfall.
    precisol(:,1)=veget_max(:,1)*(1-frac_snow_veg(:))*precip_rain(:)
    DO jv=2,nvm
       DO ji = 1, kjpindex
          ! Calculate the water that remains on the leaf as a water film
          zqsintvegnew(ji,jv) = MIN (qsintveg(ji,jv),qsintmax(ji,jv))
          ! Throughfall is composed by a the precipitation that passes through the gaps without interaction
          ! and the fraction that interacts. A share of the fraction that intercats with the canopy
          ! is expected to drip of the leaves due to leaf tipping. This fraction also contributes
          ! to the throughfall.
          precisol(ji,jv) = (veget(ji,jv)*throughfall_by_pft(jv)*(1-frac_snow_veg(ji))*precip_rain(ji)) + &
               qsintveg(ji,jv) - zqsintvegnew (ji,jv) + &
               (veget_max(ji,jv) - veget(ji,jv))*(1-frac_snow_veg(ji))*precip_rain(ji)
       ENDDO
    END DO
        
    ! Precisol is currently the same as throughfall, save it for diagnostics
    throughfall(:,:) = precisol(:,:)

    !! 1.5 Account for the contribution of snowmelt to throughfall 
    DO jv=1,nvm
       DO ji = 1, kjpindex
          IF (vegtot(ji).GT.min_sechiba) THEN
             precisol(ji,jv) = precisol(ji,jv)+tot_melt(ji)*veget_max(ji,jv)/vegtot(ji)
          ENDIF
       END DO
    END DO
    
    !! 1.6 swap qsintveg to the new value
    DO jv=2,nvm
       qsintveg(:,jv) = zqsintvegnew (:,jv)
    END DO

    IF (printlev>=3) WRITE (numout,*) ' hydrol_canop done '

  END SUBROUTINE hydrol_canop


!! ================================================================================================================================
!! SUBROUTINE 	: hydrol_vegupd
!!
!>\BRIEF        Vegetation update   
!!
!! DESCRIPTION  :
!!   The vegetation cover has changed and we need to adapt the reservoir distribution 
!!   and the distribution of plants on different soil types.
!!   You may note that this occurs after evaporation and so on have been computed. It is
!!   not a problem as a new vegetation fraction will start with humrel=0 and thus will have no
!!   evaporation. If this is not the case it should have been caught above.
!!
!! - 1 Update of vegetation is it needed?
!! - 2 calculate water mass that we have to redistribute
!! - 3 put it into reservoir of plant whose surface area has grown
!! - 4 Soil tile gestion
!! - 5 update the corresponding masks
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_vegupd

  SUBROUTINE hydrol_vegupd(kjpindex, veget, veget_max, soiltile, qsintveg, frac_bare, drain_upd, runoff_upd)


    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    ! input scalar 
    INTEGER(i_std), INTENT(in)                            :: kjpindex 
    ! input fields
    REAL(r_std), DIMENSION (kjpindex, nvm), INTENT(in)    :: veget            !! New vegetation map
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)     :: veget_max        !! Max. fraction of vegetation type
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT (in)   :: soiltile         !! Fraction of each soil tile within vegtot (0-1, unitless)

    !! 0.2 Output variables
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(out)    :: frac_bare        !! Fraction(of veget_max) of bare soil
                                                                              !! in each vegetation type
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)         :: drain_upd        !! Change in drainage due to decrease in vegtot
                                                                              !! on mc [kg/m2/dt]
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)         :: runoff_upd       !! Change in runoff due to decrease in vegtot
                                                                              !! on water2infilt[kg/m2/dt]
    

    !! 0.3 Modified variables

    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout)  :: qsintveg         !! Water on old vegetation 

    !! 0.4 Local variables

    INTEGER(i_std)                                 :: ji,jv,jst

!_ ================================================================================================================================

    !! 1 If veget has been updated at last time step (with LAND USE or DGVM),
    !! tmc and mc must be modified with respect to humtot conservation.
    CALL hydrol_tmc_update ( kjpindex, veget_max, soiltile, qsintveg, drain_upd, runoff_upd)


    ! Compute the masks for veget
    
    mask_veget(:,:) = 0
    mask_soiltile(:,:) = 0
    
    DO jst=1,nstm
       DO ji = 1, kjpindex
          IF(soiltile(ji,jst) .GT. min_sechiba) THEN
             mask_soiltile(ji,jst) = 1
          ENDIF
       END DO
    ENDDO
          
    DO jv = 1, nvm
       DO ji = 1, kjpindex
          IF(veget_max(ji,jv) .GT. min_sechiba) THEN
             mask_veget(ji,jv) = 1
          ENDIF
       END DO
    END DO

    ! Compute vegetmax_soil 
    vegetmax_soil(:,:,:) = zero
    DO jv = 1, nvm
       jst = pref_soil_veg(jv)
       DO ji=1,kjpindex
          ! for veget distribution used in sechiba via humrel
          IF (mask_soiltile(ji,jst).GT.0 .AND. vegtot(ji) > min_sechiba) THEN
             vegetmax_soil(ji,jv,jst)=veget_max(ji,jv)/soiltile(ji,jst)
          ENDIF
       ENDDO
    ENDDO

    ! Calculate frac_bare (previosly done in slowproc_veget)
    DO ji =1, kjpindex
       IF( veget_max(ji,1) .GT. min_sechiba ) THEN
          frac_bare(ji,1) = un
       ELSE
          frac_bare(ji,1) = zero
       ENDIF
    ENDDO

    IF (ok_bare_soil_new) THEN
       ! Since the flag ok_bare_soil_new no longer treats the gaps in the canopy as
       ! bare soil, frac_bare for other PFTs than 1 will be zero.
       !  Note that the same thing is done in slowproc for tot_bare_soil.
       frac_bare(:,2:nvm) = zero

    ELSE

       DO jv = 2, nvm
          DO ji =1, kjpindex
             IF( veget_max(ji,jv) .GT. min_sechiba ) THEN
                frac_bare(ji,jv) = un - veget(ji,jv)/veget_max(ji,jv)
             ELSE
                frac_bare(ji,jv) = zero
             ENDIF
          ENDDO
       ENDDO
    ENDIF

    ! Tout dans cette routine est maintenant certainement obsolete (veget_max etant constant) en dehors des lignes 
    ! suivantes et le calcul de frac_bare:
    frac_bare_ns(:,:) = zero
    DO jst = 1, nstm
       DO jv = 1, nvm
          DO ji =1, kjpindex
             IF(vegtot(ji) .GT. min_sechiba) THEN
                frac_bare_ns(ji,jst) = frac_bare_ns(ji,jst) + vegetmax_soil(ji,jv,jst) * frac_bare(ji,jv) / vegtot(ji)
             ENDIF
          END DO
       ENDDO
    END DO
    
    IF (printlev>=3) WRITE (numout,*) ' hydrol_vegupd done '

  END SUBROUTINE hydrol_vegupd


!! ================================================================================================================================
!! SUBROUTINE 	: hydrol_flood
!!
!>\BRIEF        This routine computes the evolution of the surface reservoir (floodplain).  
!!
!! DESCRIPTION  :
!! - 1 Take out vevapflo from the reservoir and transfer the remaining to subsinksoil
!! - 2 Compute the total flux from floodplain floodout (transfered to routing)
!! - 3 Discriminate between precip over land and over floodplain
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_flood

  SUBROUTINE hydrol_flood (kjpindex, vevapflo, flood_frac, flood_res, floodout)

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    ! input scalar 
    INTEGER(i_std), INTENT(in)                               :: kjpindex         !!
    ! input fields
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: flood_frac       !! Fraction of floodplains in grid box

    !! 0.2 Output variables

    REAL(r_std), DIMENSION (kjpindex), INTENT(out)           :: floodout         !! Flux to take out from floodplains

    !! 0.3 Modified variables

    REAL(r_std), DIMENSION (kjpindex), INTENT(inout)         :: flood_res        !! Floodplains reservoir estimate
    REAL(r_std), DIMENSION (kjpindex), INTENT(inout)         :: vevapflo         !! Evaporation over floodplains

    !! 0.4 Local variables

    INTEGER(i_std)                                           :: ji, jv           !! Indices
    REAL(r_std), DIMENSION (kjpindex)                        :: temp             !! 

!_ ================================================================================================================================
    !- 
    !! 1 Take out vevapflo from the reservoir and transfer the remaining to subsinksoil 
    !-
    DO ji = 1,kjpindex
       temp(ji) = MIN(flood_res(ji), vevapflo(ji))
    ENDDO
    DO ji = 1,kjpindex
       flood_res(ji) = flood_res(ji) - temp(ji)
       subsinksoil(ji) = subsinksoil(ji) + vevapflo(ji) - temp(ji)
       vevapflo(ji) = temp(ji)
    ENDDO

    !- 
    !! 2 Compute the total flux from floodplain floodout (transfered to routing) 
    !-
    DO ji = 1,kjpindex
       floodout(ji) = vevapflo(ji) - flood_frac(ji) * SUM(precisol(ji,:))
    ENDDO

    !-
    !! 3 Discriminate between precip over land and over floodplain
    !-
    DO jv=1, nvm
       DO ji = 1,kjpindex
          precisol(ji,jv) = precisol(ji,jv) * (1 - flood_frac(ji))
       ENDDO
    ENDDO 

    IF (printlev>=3) WRITE (numout,*) ' hydrol_flood done'

  END SUBROUTINE hydrol_flood


!! ================================================================================================================================
!! SUBROUTINE 	: hydrol_soil
!!
!>\BRIEF        This routine computes soil processes with CWRR scheme (Richards equation solved by finite differences).
!! Note that the water fluxes are in kg/m2/dt_sechiba.
!!
!! DESCRIPTION  :
!! 0. Initialisation, and split 2d variables to 3d variables, per soil tile
!! -- START MAIN LOOP (prognostic loop to update mc and mcl) OVER SOILTILES
!! 1. FIRSTLY, WE CHANGE MC BASED ON EXTERNAL FLUXES, ALL APPLIED AT THE SOIL SURFACE
!! 1.1 Reduces water2infilt and water2extract to their difference 
!! 1.2 To remove water2extract (including bare soilevaporation) from top layer
!! 1.3 Infiltration
!! 1.4 Reinfiltration of surface runoff : compute temporary surface water and extract from runoff
!! 2. SECONDLY, WE UPDATE MC FROM DIFFUSION, INCLUDING DRAINAGE AND ROOTSINK
!!    This will act on mcl (liquid water content) only
!! 2.1 K and D are recomputed after infiltration
!! 2.2 Set the tridiagonal matrix coefficients for the diffusion/redistribution scheme
!! 2.3 We define mcl (liquid water content) based on mc and profil_froz_hydro_ns
!! 2.4 We calculate the total SM at the beginning of the routine tridiag for water conservation check
!! 2.5 Defining where diffusion is solved : everywhere
!! 2.6 We define the system of linear equations for mcl redistribution
!! 2.7 Solves diffusion equations
!! 2.8 Computes drainage = bottom boundary condition, consistent with rhs(ji,jsl=nslm)
!! 2.9 For water conservation check during redistribution, we calculate the total liquid SM
!!     at the end of the routine tridiag, and we compare the difference with the flux...
!! 3. AFTER DIFFUSION/REDISTRIBUTION
!! 3.1 Updating mc, as all the following checks against saturation will compare mc to mcs
!! 3.2 Correct here the possible over-saturation values (subroutine hydrol_soil_smooth_over_mcs2 acts on mc)
!!     Here hydrol_soil_smooth_over_mcs2 discard all excess as ru_corr_ns, oriented to either ru_ns or dr_ns
!! 3.3 Negative runoff is reported to drainage
!! 3.4 Optional block to force saturation below zwt_force
!! 3.5 Diagnosing the effective water table depth
!! 3.6 Diagnose under_mcr to adapt water stress calculation below
!! 4. At the end of the prognostic calculations, we recompute important moisture variables
!! 4.1 Total soil moisture content (water2infilt added below)
!! 4.2 mcl is a module variable; we update it here for calculating bare soil evaporation,
!! 5. Optional check of the water balance of soil column (if check_cwrr)
!! 5.1 Computation of the vertical water fluxes
!! 6. SM DIAGNOSTICS FOR OTHER ROUTINES, MODULES, OR NEXT STEP
!! 6.1 Total soil moisture, soil moisture at litter levels, soil wetness, us, humrelv, vesgtressv
!! 6.2 We need to turn off evaporation when is_under_mcr
!! 6.3 Calculate the volumetric soil moisture content (mc_layh and mcl_layh) needed in thermosoil
!! 6.4 The hydraulic conductivities exported here are the ones used in the diffusion/redistribution
!! -- ENDING THE MAIN LOOP ON SOILTILES
!! 7. Summing 3d variables into 2d variables
!! 8. XIOS export of local variables, including water conservation checks
!! 9. COMPUTING EVAP_BARE_LIM_NS FOR NEXT TIME STEP, WITH A LOOP ON SOILTILES
!!    The principle is to run a dummy integration of the water redistribution scheme
!!    to check if the SM profile can sustain a potential evaporation.
!!    If not, the dummy integration is redone from the SM profile of the end of the normal integration,
!!    with a boundary condition leading to a very severe water limitation: mc(1)=mcr
!! 10. evap_bar_lim is the grid-cell scale beta
!!
!! RECENT CHANGE(S) : 2016 by A. Ducharne
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) :
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE hydrol_soil (ks, nvan, avan, mcr, mcs, mcfc, mcw, & 
       & kjpindex, veget, veget_max, soiltile, njsc, reinf_slope_soil, &
       & transpir, vevapnu, evapot, evapot_penm, runoff, drainage, &
       & returnflow, reinfiltration, irrigation, &
       & tot_melt, evap_bare_lim, evap_bare_lim_ns, shumdiag, shumdiag_perma,&
       & k_litt, litterhumdiag, humrel,vegstress, drysoil_frac, &
       & stempdiag,snow, &
       & snowdz, tot_bare_soil, u, v, tq_cdrag, mc_layh, mcl_layh, mc_layh_s, &
       & mcl_layh_s, e_frac, ksoil, altmax, root_profile, root_depth, root_deficit, &
       & circ_class_biomass, us, precip_rain, totfrac_nobio, frac_snow_nobio, F_absorption, &
       & wtp,mc_peat_above,liqwt_ratio,shumdiag_peat,shumdiag_croppeat,mc_croppeat_above,&
       & shumdiag_man,mc_man_above)

    ! 
    ! interface description

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                               :: kjpindex
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)        :: veget            !! Fraction of vegetation type
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT (in)       :: veget_max        !! Map of max vegetation types [-]
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)         :: njsc             !! Index of the dominant soil textural class 
                                                                                 !! in the grid cell (1-nscm, unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: ks               !! Hydraulic conductivity at saturation (mm {-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: nvan             !! Van Genuchten coeficients n (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: avan             !! Van Genuchten coeficients a (mm-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: mcr              !! Residual volumetric water content (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: mcs              !! Saturated volumetric water content (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: mcfc             !! Volumetric water content at field capacity (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: mcw              !! Volumetric water content at wilting point (m^{3} m^{-3})
    REAL(r_std), DIMENSION (kjpindex, nstm), INTENT (in)     :: soiltile         !! Fraction of each soil tile within vegtot (0-1, unitless)
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(in)        :: transpir         !! Transpiration  
                                                                                 !!  @tex $(kg m^{-2} dt\_sechiba^{-1})$ @endtex
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(in)        :: F_absorption     !! Total root absorption (ok_hydrol_arch = .TRUE.)  
                                                                                 !!  @tex $(m^3 s^{-1})$ @endtex
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT (in)      :: reinf_slope_soil !! Fraction of surface runoff that reinfiltrates per soil tile
                                                                                 !!  (unitless, [0-1])
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: returnflow       !! Water returning to the soil from the bottom
                                                                                 !!  @tex $(kg m^{-2} dt\_sechiba^{-1})$ @endtex
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: reinfiltration   !! Water returning to the top of the soil
                                                                                 !!  @tex $(kg m^{-2} dt\_sechiba^{-1})$ @endtex
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: irrigation       !! Irrigation
                                                                                 !!  @tex $(kg m^{-2} dt\_sechiba^{-1})$ @endtex
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: evapot           !! Potential evaporation
                                                                                 !!  @tex $(kg m^{-2} dt\_sechiba^{-1})$ @endtex
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: evapot_penm      !! Potential evaporation "Penman" (Milly's correction)
                                                                                 !!  @tex $(kg m^{-2} dt\_sechiba^{-1})$ @endtex
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: tot_melt         !! Total melt from snow and ice
                                                                                 !!  @tex $(kg m^{-2} dt\_sechiba^{-1})$ @endtex
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (in)       :: stempdiag        !! Diagnostic temp profile from thermosoil
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: snow             !! Snow mass
                                                                                 !!  @tex $(kg m^{-2})$ @endtex
    REAL(r_std), DIMENSION (kjpindex,nsnow),INTENT(in)       :: snowdz           !! Snow depth (m)
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: tot_bare_soil    !! Total evaporating bare soil fraction 
                                                                                 !!  (unitless, [0-1])
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)             :: u,v              !! Horizontal wind speed
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)             :: tq_cdrag         !! Surface drag coefficient 
    REAL(r_std), DIMENSION(kjpindex,nvm,nslm,nstm),INTENT(in):: e_frac           !! Fraction of water transpired supplied by individual layers (no units)
    REAL(r_std),DIMENSION (kjpindex,nvm),INTENT(in)          :: altmax           !! Maximul active layer thickness (m). Be careful, here active means non frozen.
                                                                                 !! Not related with the active soil carbon pool.
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(in)            :: circ_class_biomass !! Biomass components of the model tree  
                                                                                 !! within a circumference class
                                                                                 !! class @tex $(g C ind^{-1})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: precip_rain      !! Rain precipitation				
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: totfrac_nobio    !! Total fraction of continental ice+lakes+...	
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT(in)      :: frac_snow_nobio  !! Snow cover fraction on non-vegeted area	


    !! 0.2 Output variables

    REAL(r_std), DIMENSION (kjpindex), INTENT(out)                   :: runoff            !! Surface runoff
                                                                                          !!  @tex $(kg m^{-2} dt\_sechiba^{-1})$ @endtex
    REAL(r_std), DIMENSION (kjpindex), INTENT(out)                   :: drainage          !! Drainage
                                                                                          !!  @tex $(kg m^{-2} dt\_sechiba^{-1})$ @endtex
    REAL(r_std), DIMENSION (kjpindex), INTENT(out)                   :: evap_bare_lim     !! Limitation factor (beta) for bare soil evaporation  
                                                                                          !! on each soil column (unitless, [0-1])
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT(out)              :: evap_bare_lim_ns  !! Limitation factor (beta) for bare soil evaporation  
                                                                                          !! on each soil column (unitless, [0-1])
    REAL(r_std), DIMENSION (kjpindex,nslm), INTENT (out)             :: shumdiag          !! Relative soil moisture in each diag soil layer 
                                                                                          !! with respect to (mcfc-mcw) (unitless, [0-1])
    REAL(r_std), DIMENSION (kjpindex,nslm), INTENT (out)             :: shumdiag_perma    !! Percent of porosity filled with water (mc/mcs)
                                                                                          !! in each diag soil layer (for the thermal computations)
                                                                                          !! (unitless, [0-1])
    REAL(r_std), DIMENSION (kjpindex), INTENT (out)                  :: k_litt            !! Litter approximated hydraulic conductivity
                                                                                          !!  @tex $(mm d^{-1})$ @endtex
    REAL(r_std), DIMENSION (kjpindex), INTENT (out)                  :: litterhumdiag     !! Mean of soil_wet_litter across soil tiles
                                                                                          !! (unitless, [0-1])
    REAL(r_std), DIMENSION (kjpindex, nvm), INTENT(out)              :: vegstress         !! Veg. moisture stress (only for vegetation 
                                                                                          !! growth) (unitless, [0-1])
    REAL(r_std), DIMENSION (kjpindex), INTENT (out)                  :: drysoil_frac      !! Function of the litter humidity
    REAL(r_std), DIMENSION (kjpindex,nslm), INTENT (out)             :: mc_layh           !! Volumetric water content (liquid + ice) for each soil layer
                                                                                          !! averaged over the mesh (for thermosoil)
                                                                                          !!  @tex $(m^{3} m^{-3})$ @endtex
    REAL(r_std), DIMENSION (kjpindex,nslm), INTENT (out)             :: mcl_layh          !! Volumetric liquid water content for each soil layer 
                                                                                          !! averaged over the mesh (for thermosoil)
                                                                                          !!  @tex $(m^{3} m^{-3})$ @endtex 
    REAL(r_std), DIMENSION (kjpindex,nslm,nstm), INTENT (out)         :: mc_layh_s        !! Volumetric soil moisture content for each layer in hydrol(liquid + ice) [m3/m3]
    REAL(r_std), DIMENSION (kjpindex,nslm,nstm), INTENT (out)         :: mcl_layh_s       !! Volumetric soil moisture content for each layer in hydrol(liquid) [m3/m3]
    REAL(r_std),DIMENSION (kjpindex,nvm,nslm,nroot_prof), INTENT(out) :: root_profile     !! Normalized root mass/length fraction in each soil layer
                                                                                          !! (0-1, unitless)
    REAL(r_std), DIMENSION (kjpindex,nvm,ndepths), INTENT(out)        :: root_depth       !! Node and interface numbers at which the deepest roots
                                                                                          !! occur (1 to nslm, unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)                     :: root_deficit     !! water deficit to reach SM target of soil column, for irrigation demand
    
    ! peatland
    REAL(r_std), DIMENSION (kjpindex), INTENT(out)                    :: wtp              !! water table position, when peat is activated
    REAL(r_std), DIMENSION (kjpindex), INTENT(out)                    :: liqwt_ratio      !! liquid water ratio in peat.
    REAL(r_std), DIMENSION (kjpindex,nslm)                            :: h_eau            !! relative moisture useful for water table computation
    REAL(r_std), DIMENSION (kjpindex)                                 :: tmcl_peat        !! total liquid water content over peat tile when the peat is activated
    REAL(r_std)                                                       :: tmcs_peat        !! satuated water content over peat tile when the peat is activated
    REAL(r_std), DIMENSION (kjpindex), INTENT (out)                   :: mc_peat_above    !! liquid soil moisture content within the first 4 layers when peat is activated
    REAL(r_std), DIMENSION (kjpindex), INTENT (out)                   :: mc_croppeat_above!! liquid soil moisture content within the first 4 layers when agri_peat is activated
    REAL(r_std), DIMENSION (kjpindex), INTENT (out)                   :: mc_man_above     !! liquid soil moisture content within the first 4 layers when tides is activated
    REAL(r_std), DIMENSION (kjpindex,nslm), INTENT (out)              :: shumdiag_peat    !! soil moisture content for the decomposition when peat is activated
    REAL(r_std), DIMENSION (kjpindex,nslm), INTENT (out)              :: shumdiag_croppeat!! soil moisture content for the decomposition when agri_peat is activated
    REAL(r_std), DIMENSION (kjpindex,nslm), INTENT (out)              :: shumdiag_man     !! soil moisture content for the decomposition when tides is activated
    REAL(r_std), DIMENSION (kjpindex,nstm)                            :: tmc_soil         !! total water content without including water2infilt if peat is activated

    !! 0.3 Modified variables

    REAL(r_std), DIMENSION (kjpindex), INTENT(inout)         :: vevapnu                   !! Bare soil evaporation
                                                                                          !!  @tex $(kg m^{-2} dt\_sechiba^{-1})$ @endtex
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT (inout)    :: humrel                    !! Relative humidity (0-1, dimensionless)
    REAL(r_std), DIMENSION (kjpindex,nslm,nstm), INTENT(inout) :: ksoil                   !! Soil conductivity (a copy of k for each soil type)
    REAL(r_std),DIMENSION (kjpindex,nvm,nstm,nslm), INTENT(inout) :: us                   !! Water stress index for transpiration
                                                                                          !! (by soil layer and PFT) (0-1, unitless)

    

    !! 0.4 Local variables

    INTEGER(i_std)                                 :: ji, jv, jsl, jst           !! Indice
    INTEGER(i_std)                                 :: jst_kfact_root             !! Indice for kfact_root calculation
    REAL(r_std), PARAMETER                         :: frac_mcs = 0.66            !! Temporary depth
    REAL(r_std), DIMENSION(kjpindex)               :: temp                       !! Temporary value for fluxes
    REAL(r_std), DIMENSION(kjpindex)               :: tmcold                     !! Total SM at beginning of hydrol_soil (kg/m2)
    REAL(r_std), DIMENSION(kjpindex)               :: tmcint                     !! Ancillary total SM (kg/m2)
    REAL(r_std), DIMENSION(kjpindex,nslm)          :: mcint                      !! To save mc values for future use
    REAL(r_std), DIMENSION(kjpindex,nslm)          :: mclint                     !! To save mcl values for future use
    LOGICAL, DIMENSION(kjpindex,nstm)              :: is_under_mcr               !! Identifies under residual soil moisture points
    REAL(r_std), DIMENSION(kjpindex)               :: deltahum,diff              !!
    LOGICAL(r_std), DIMENSION(kjpindex)            :: test                       !!
    REAL(r_std), DIMENSION(kjpindex)               :: water2extract              !! Water flux to be extracted at the soil surface
                                                                                 !!  @tex $(kg m^{-2} dt\_sechiba^{-1})$ @endtex
    REAL(r_std), DIMENSION(kjpindex)               :: returnflow_soil            !! Water from the routing back to the bottom of 
                                                                                 !! the soil @tex $(kg m^{-2} dt\_sechiba^{-1})$ @endtex
    REAL(r_std), DIMENSION(kjpindex)               :: reinfiltration_soil        !! Water from the routing back to the top of the 
                                                                                 !! soil @tex $(kg m^{-2} dt\_sechiba^{-1})$ @endtex
    REAL(r_std), DIMENSION(kjpindex,nstm)          :: irrigation_soil            !! Water from irrigation returning to soil moisture per soil tile
                                                                                 !!  @tex $(kg m^{-2} dt\_sechiba^{-1})$ @endtex
    REAL(r_std), DIMENSION(kjpindex)               :: flux_infilt                !! Water to infiltrate
                                                                                 !!  @tex $(kg m^{-2})$ @endtex
    REAL(r_std), DIMENSION(kjpindex)               :: flux_bottom                !! Flux at bottom of the soil column
                                                                                 !!  @tex $(kg m^{-2})$ @endtex
    REAL(r_std), DIMENSION(kjpindex)               :: flux_top                   !! Flux at top of the soil column (for bare soil evap)
                                                                                 !!  @tex $(kg m^{-2})$ @endtex
    REAL(r_std), DIMENSION (kjpindex,nstm)         :: qinfilt_ns                 !! Effective infiltration flux per soil tile
                                                                                 !!  @tex $(kg m^{-2})$ @endtex   
    REAL(r_std), DIMENSION (kjpindex)              :: qinfilt                    !! Effective infiltration flux  
                                                                                 !!  @tex $(kg m^{-2})$ @endtex 
    REAL(r_std), DIMENSION (kjpindex,nstm)         :: ru_infilt_ns               !! Surface runoff from hydrol_soil_infilt per soil tile
                                                                                 !!  @tex $(kg m^{-2})$ @endtex   
    REAL(r_std), DIMENSION (kjpindex)              :: ru_infilt                  !! Surface runoff from hydrol_soil_infilt 
                                                                                 !!  @tex $(kg m^{-2})$ @endtex 
    REAL(r_std), DIMENSION (kjpindex,nstm)         :: ru_corr_ns                 !! Surface runoff produced to correct excess per soil tile
                                                                                 !!  @tex $(kg m^{-2})$ @endtex
    REAL(r_std), DIMENSION (kjpindex)              :: ru_corr                    !! Surface runoff produced to correct excess 
                                                                                 !!  @tex $(kg m^{-2} dt\_sechiba^{-1})$ @endtex  
    REAL(r_std), DIMENSION (kjpindex,nstm)         :: ru_corr2_ns                !! Correction of negative surface runoff per soil tile
                                                                                 !!  @tex $(kg m^{-2})$ @endtex
    REAL(r_std), DIMENSION (kjpindex)              :: ru_corr2                   !! Correction of negative surface runoff
                                                                                 !!  @tex $(kg m^{-2})$ @endtex 
    REAL(r_std), DIMENSION (kjpindex,nstm)         :: dr_corr_ns                 !! Drainage produced to correct excess
                                                                                 !!  @tex $(kg m^{-2})$ @endtex
    REAL(r_std), DIMENSION (kjpindex,nstm)         :: dr_corrnum_ns              !! Drainage produced to correct numerical errors in tridiag
                                                                                 !!  @tex $(kg m^{-2})$ @endtex    
    REAL(r_std), DIMENSION (kjpindex)              :: dr_corr                    !! Drainage produced to correct excess 
                                                                                 !!  @tex $(kg m^{-2} dt\_sechiba^{-1})$ @endtex
    REAL(r_std), DIMENSION (kjpindex)              :: dr_corrnum                 !! Drainage produced to correct numerical errors in tridiag
                                                                                 !!  @tex $(kg m^{-2} dt\_sechiba^{-1})$ @endtex
    REAL(r_std), DIMENSION (kjpindex,nslm)         :: dmc                        !! Delta mc when forcing saturation (zwt_force) 
                                                                                 !!  @tex $(m^{3} m^{-3})$ @endtex 
    REAL(r_std), DIMENSION (kjpindex,nstm)         :: dr_force_ns                !! Delta drainage when forcing saturation (zwt_force) 
                                                                                 !!  per soil tile  @tex $(kg m^{-2})$ @endtex
    REAL(r_std), DIMENSION (kjpindex)              :: dr_force                   !! Delta drainage when forcing saturation (zwt_force)
                                                                                 !!  @tex $(kg m^{-2})$ @endtex  
    REAL(r_std), DIMENSION (kjpindex,nstm)         :: wtd_ns                     !! Effective water table depth (m)

    REAL(r_std), DIMENSION (kjpindex)              :: wtd                        !! Mean water table depth in the grid-cell (m)

    ! For the calculation of soil_wet_ns and us/humrel/vegstress 
    
    REAL(r_std), DIMENSION (kjpindex,nstm)         :: humtot_soil_liquid         !! Soil water content per soiltile 
    
    REAL(r_std), DIMENSION (kjpindex,nslm,nstm)    :: sm                         !! Soil moisture of each layer per soil tile (liquid phase)
                                                                                 !!  @tex $(kg m^{-2})$ @endtex
    REAL(r_std), DIMENSION (kjpindex,nslm,nstm)    :: smt                        !! Soil moisture of each layer per soil tile (liquid+solid phase)
                                                                                 !!  @tex $(kg m^{-2})$ @endtex
    REAL(r_std), DIMENSION (kjpindex,nslm)         :: smw                        !! Soil moisture of each layer at wilting point
                                                                                 !!  @tex $(kg m^{-2})$ @endtex 
    REAL(r_std), DIMENSION (kjpindex,nslm)         :: smf                        !! Soil moisture of each layer at field capacity
                                                                                 !!  @tex $(kg m^{-2})$ @endtex   
    REAL(r_std), DIMENSION (kjpindex,nslm)         :: sms                        !! Soil moisture of each layer at saturation
                                                                                 !!  @tex $(kg m^{-2})$ @endtex 
    REAL(r_std), DIMENSION (kjpindex,nslm)         :: sm_nostress                !! Soil moisture of each layer at which us reaches 1
                                                                                 !!  @tex $(kg m^{-2})$ @endtex 
    ! For water conservation checks (in mm/dtstep unless otherwise mentioned)
    REAL(r_std), DIMENSION (kjpindex,nstm)         :: check_infilt_ns             !! Water conservation diagnostic at routine scale
    REAL(r_std), DIMENSION (kjpindex,nstm)         :: check1_ns                   !! Water conservation diagnostic at routine scale
    REAL(r_std), DIMENSION (kjpindex,nstm)         :: check_tr_ns                 !! Water conservation diagnostic at routine scale
    REAL(r_std), DIMENSION (kjpindex,nstm)         :: check_over_ns               !! Water conservation diagnostic at routine scale
    REAL(r_std), DIMENSION (kjpindex,nstm)         :: check_under_ns              !! Water conservation diagnostic at routine scale
    REAL(r_std), DIMENSION(kjpindex)               :: tmci                        !! Total soil moisture at beginning of routine (kg/m2)
    REAL(r_std), DIMENSION(kjpindex)               :: tmcf                        !! Total soil moisture at end of routine (kg/m2)
    REAL(r_std), DIMENSION(kjpindex)               :: diag_tr                     !! Transpiration flux
    REAL(r_std), DIMENSION (kjpindex)              :: check_infilt                !! Water conservation diagnostic at routine scale 
    REAL(r_std), DIMENSION (kjpindex)              :: check1                      !! Water conservation diagnostic at routine scale 
    REAL(r_std), DIMENSION (kjpindex)              :: check_tr                    !! Water conservation diagnostic at routine scale 
    REAL(r_std), DIMENSION (kjpindex)              :: check_over                  !! Water conservation diagnostic at routine scale 
    REAL(r_std), DIMENSION (kjpindex)              :: check_under                 !! Water conservation diagnostic at routine scale
    ! For irrigation triggering
    INTEGER(i_std), DIMENSION (kjpindex)           :: lai_irrig_trig              !! Number of PFT per cell with LAI> LAI_IRRIG_MIN -
    ! Diagnostic of the vertical soil water fluxes  
    REAL(r_std), DIMENSION (kjpindex,nslm)         :: qflux                       !! Local upward flux into soil layer 
                                                                                  !! from lower interface
                                                                                  !!  @tex $(kg m^{-2})$ @endtex 
    REAL(r_std), DIMENSION (kjpindex)              :: check_top                   !! Water budget residu in top soil layer
                                                                                  !!  @tex $(kg m^{-2})$ @endtex 

    ! Variables for calculation of a soil resistance, option do_rsoil (following the formulation of Sellers et al 1992, implemented in Oleson et al. 2008)
    REAL(r_std)                     		   :: speed			 !! magnitude of wind speed required for Aerodynamic resistance 
    REAL(r_std)                                    :: ra			 !! diagnosed aerodynamic resistance
    REAL(r_std), DIMENSION(kjpindex)               :: mc_rel                     !! first layer relative soil moisture, required for rsoil
    REAL(r_std), DIMENSION(kjpindex)               :: evap_soil                  !! soil evaporation from Oleson et al 2008
    REAL(r_std), DIMENSION(kjpindex,nstm)          :: r_soil_ns	                 !! soil resistance from Oleson et al 2008
    REAL(r_std), DIMENSION(kjpindex)               :: r_soil	                 !! soil resistance from Oleson et al 2008
    REAL(r_std), DIMENSION(kjpindex)               :: tmcs_litter                !! Saturated soil moisture in the 4 "litter" soil layers
    REAL(r_std), DIMENSION(nslm)                   :: root_profile_tmp           !! Temporary variable to calculate the root_profile

    ! For CMIP6 and SP-MIP : ksat and matric pressure head psi(theta)
    REAL(r_std)                                    :: mc_ratio, mvg, avg, mc_ratio_peat, mvg_peat, avg_peat
    REAL(r_std)                                    :: psi,  psi_peat             !! Matric head (per soil layer and soil tile) [mm=kg/m2]
    REAL(r_std), DIMENSION (kjpindex,nslm)         :: psi_moy                    !! Mean matric head per soil layer [mm=kg/m2]  
    REAL(r_std), DIMENSION (kjpindex,nslm)         :: ksat                       !! Saturated hydraulic conductivity at each node (mm/d)  
    REAL(r_std), DIMENSION (kjpindex,nvm,nslm,nroot_prof) :: tmp                 !! temporary variable for writing the root profiles to XIOS

    !
    ! peat and tides:
    REAL(r_std)					   :: wtp_tide			 !! water table tide forcing if tides is true
    REAL(r_std)					   :: dwtp_tide  		 !! change in water table tide forcing if tides is true (wtp_tide[t]-wtp_tide[t-1])
    REAL(r_std),DIMENSION(nslm)                    :: nfact_peat                 !! Multiplicative factor for decay of n with depth, for peatland
    REAL(r_std)                                    :: nvan_mod_peat              !! VG parameter n  modified from  exponantial profile,for peat
    !
!_ ================================================================================================================================

    !! 0.1 Arrays with DIMENSION(kjpindex)
    
    returnflow_soil(:) = zero
    reinfiltration_soil(:) = zero
    irrigation_soil(:,:) = zero
    qflux_ns(:,:,:) = zero
    mc_layh(:,:) = zero ! for thermosoil
    mcl_layh(:,:) = zero ! for thermosoil
    kk(:,:,:) = zero
    kk_moy(:,:) = zero
    undermcr(:) = zero ! needs to be initialized outside from jst loop
    ksat(:,:) = zero
    psi_moy(:,:) = zero
    humtot_soil_liquid(:,:) = zero

    !! Calculate kfact_root
    IF ( kfact_root_const ) THEN
       kfact_root(:,:,:) = un
    ELSE
       !! An exponential factor is used to increase ks near the surface depending on the amount of roots in the soil 
       !! through a geometric average over the vegets
       !! This comes from the PhD thesis of d'Orgeval, 2006, p82; d'Orgeval et al. 2008, Eqs. 3-4
       !! (Calibrated against Hapex-Sahel measurements)
       !! Since rev 2916: veget_max/2 is used instead of veget
       kfact_root(:,:,:) = un
       DO jsl = 1, nslm
          DO jv = 2, nvm
             jst_kfact_root = pref_soil_veg(jv)
             DO ji = 1, kjpindex
                IF (soiltile(ji,jst_kfact_root) .GT. min_sechiba) THEN
                   ! peatland : 
                   ! no need to modify for the case of peat pft 0, because of the soiltile condition just above
                   ! also because the kfact_root is initialized to 1 already.
                   ! Note: ks_peat = 2120, ks = ks_usda if we do not use sp_mip soil texteure
                   ! ks_peat is defined between usda loamy sand and sandy loam, with quite high value of ks.
                   ! Besides, vegetmax is considered in the computation of kfact_root. so it computes only the peatland according to its proportion. 
                   IF ( (ok_peat_hydro) .AND. ((jst_kfact_root  .EQ. 4 ) .OR. (jst_kfact_root  .EQ. 5))) THEN
                      IF (soiltile(ji,jst_kfact_root ) .GT. min_sechiba) THEN
                         kfact_root(ji,jsl,jst_kfact_root) = kfact_root(ji,jsl,jst_kfact_root ) *&
                              & MAX((MAXVAL(ks_usda)/ks_peat)**(- vegetmax_soil(ji,jv,jst_kfact_root )/2 * (humcste(jv)*zz(jsl)/mille - un)/deux), un)
                      ENDIF
                   ELSE
                      kfact_root(ji,jsl,jst_kfact_root) = kfact_root(ji,jsl,jst_kfact_root) * &
                           MAX((MAXVAL(ks_usda)/ks(ji))**(- vegetmax_soil(ji,jv,jst_kfact_root)/2 * (humcste(jv)*zz(jsl)/mille - un)/deux), un)
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    END IF



    IF (ok_freeze_cwrr) THEN
       
       ! 0.1 Calculate the temperature and fozen fraction at the hydrological levels
       ! Calculates profil_froz_hydro_ns as a function of stempdiag and mc if ok_thermodynamical_freezing
       ! These values will be kept till the end of the prognostic loop
       DO jst=1,nstm 
          CALL hydrol_soil_froz(nvan, avan, mcr, mcs, kjpindex,jst,stempdiag, soiltile)
       ENDDO

    ELSE
  
       profil_froz_hydro_ns(:,:,:) = zero
             
    ENDIF
    
    !! 0.2 Split 2d variables to 3d variables, per soil tile
    !  Here, the evaporative fluxes are distributed over the soiltiles as a function of the 
    !    corresponding control factors; they are normalized to vegtot
    !  At step 7, the reverse transformation is used for the fluxes produced in hydrol_soil
    !    flux_cell(ji)=sum(flux_ns(ji,jst)*soiltile(ji,jst)*vegtot(ji))


    CALL hydrol_split_soil (kjpindex, veget_max, soiltile, vevapnu, transpir, humrel, &
         evap_bare_lim, evap_bare_lim_ns, tot_bare_soil, us, e_frac, F_absorption)

    
    !! 0.3 Common variables related to routing, with all return flow applied to the soil surface
    ! The fluxes coming from the routing are uniformly splitted into the soiltiles,
    !    but are normalized to vegtot like the above fluxes: 
    !            flux_ns(ji,jst)=flux_cell(ji)/vegtot(ji)
    ! It is the case for : irrigation_soil(ji) and reinfiltration_soil(ji) cf below
    ! It is also the case for subsinksoil(ji), which is divided by (1-tot_frac_nobio) at creation in hydrol_snow
    ! AD16*** The transformation in 0.2 and 0.3 is likely to induce conservation problems
    !         when tot_frac_nobio NE 0, since sum(soiltile) NE vegtot in this case
    IF (.NOT. old_irrig_scheme) THEN !
       IF (.NOT. irrigated_soiltile) THEN
          DO ji=1,kjpindex
             IF(vegtot(ji).GT.min_sechiba ) THEN
                returnflow_soil(ji) = zero
                reinfiltration_soil(ji) = (returnflow(ji) + reinfiltration(ji))/vegtot(ji)
                IF(soiltile(ji, irrig_st).GT.min_sechiba) THEN
                   !irrigation_soil(ji, 1:2) = 0, if irrig_st = 3. Not put because Values
                   !are already zero, due to initialization
                   irrigation_soil(ji, irrig_st) = irrigation(ji) / (soiltile(ji, irrig_st) * vegtot(ji) )
                   !Irrigation is kg/m2 of grid cell. Here, all that water is put on
                   !irrig_st (irrigated soil tile), by default = 3, for the others
                   !soil tiles irrigation = zero
                ENDIF
             ENDIF
          ENDDO
       ENDIF
    ELSE  ! old_irrig_scheme is used
       DO ji=1,kjpindex
          IF(vegtot(ji).GT.min_sechiba) THEN
             ! returnflow_soil is assumed to enter from the bottom, but it is not possible with CWRR
             returnflow_soil(ji) = zero
             reinfiltration_soil(ji) = (returnflow(ji) + reinfiltration(ji))/vegtot(ji)
             irrigation_soil(ji,:) = irrigation(ji)/vegtot(ji)
             ! irrigation_soil(ji) = irrigation(ji)/vegtot(ji)
             ! Computed for all the grid cell. New way is equivalent, and coherent
             ! with irrigation_soil new dimensions (cells, soil tiles)
             ! Irrigation is kg/m2 of grid cell. For the old irrig. scheme,
             ! irrigation soil is the same for every soil tile
             ! Next lines are in tag 2.0, deleted because values are already init to zero
             ! ELSE
             ! returnflow_soil(ji) = zero
             ! reinfiltration_soil(ji) = zero
             ! irrigation_soil(ji) = zero
             ! ENDIF
          ENDIF
       ENDDO
    ENDIF
    
    !! -- START MAIN LOOP (prognostic loop to update mc and mcl) OVER SOILTILES
    !!    The called subroutines work on arrays with DIMENSION(kjpindex),
    !!    recursively used for each soiltile jst
    
    DO jst = 1,nstm

       is_under_mcr(:,jst) = .FALSE.
       
       !! 0.4. Keep initial values for future check-up
       
       ! Total moisture content (including water2infilt) is saved for balance checks at the end
       ! In hydrol_tmc_update, tmc is increased by water2infilt(ji,jst), but mc is not modified !
       tmcold(:) = tmc(:,jst)
       
       ! The value of mc is kept in mcint (nstm dimension removed), in case needed for water balance checks 
       DO jsl = 1, nslm
          DO ji = 1, kjpindex
             mcint(ji,jsl) = mask_soiltile(ji,jst) * mc(ji,jsl,jst)
          ENDDO
       ENDDO
       !
       ! Initial total moisture content : tmcint does not include water2infilt, contrarily to tmcold
       DO ji = 1, kjpindex
          tmcint(ji) = dz(2) * ( trois*mcint(ji,1) + mcint(ji,2) )/huit 
       ENDDO
       DO jsl = 2,nslm-1
          DO ji = 1, kjpindex
             tmcint(ji) = tmcint(ji) + dz(jsl) &
                  & * (trois*mcint(ji,jsl)+mcint(ji,jsl-1))/huit &
                  & + dz(jsl+1) * (trois*mcint(ji,jsl)+mcint(ji,jsl+1))/huit
          ENDDO
       ENDDO
       DO ji = 1, kjpindex
          tmcint(ji) = tmcint(ji) + dz(nslm) &
               & * (trois * mcint(ji,nslm) + mcint(ji,nslm-1))/huit
       ENDDO

       !! 1. FIRSTLY, WE CHANGE MC BASED ON EXTERNAL FLUXES, ALL APPLIED AT THE SOIL SURFACE
       !!   Input = water2infilt(ji,jst) + irrigation_soil(ji) + reinfiltration_soil(ji) + precisol_ns(ji,jst)
       !!      - negative evaporation fluxes (MIN(ae_ns(ji,jst),zero)+ MIN(subsinksoil(ji),zero))
       !!   Output = MAX(ae_ns(ji,jst),zero) + subsinksoil(ji) = positive evaporation flux = water2extract
       ! In practice, negative subsinksoil(ji) is not possible 

       !! 1.1 Reduces water2infilt and water2extract to their difference 

       ! Compares water2infilt and water2extract to keep only difference
       ! Here, temp is used as a temporary variable to store the min of water to infiltrate vs evaporate
       DO ji = 1, kjpindex
          temp(ji) = MIN(water2infilt(ji,jst) + irrigation_soil(ji,jst) + reinfiltration_soil(ji) &
                         - MIN(ae_ns(ji,jst),zero) - MIN(subsinksoil(ji),zero) + precisol_ns(ji,jst), &
                           MAX(ae_ns(ji,jst),zero) + MAX(subsinksoil(ji),zero) )
       ENDDO

       ! The water to infiltrate at the soil surface is either 0, or the difference to what has to be evaporated
       !   - the initial water2infilt (right hand side) results from qsintveg changes with vegetation updates
       !   - irrigation_soil is the input flux to the soil surface from irrigation 
       !   - reinfiltration_soil is the input flux to the soil surface from routing 'including returnflow)
       !   - eventually, water2infilt holds all fluxes to the soil surface except precisol (reduced by water2extract)
       DO ji = 1, kjpindex
         !Note that in tag 2.0, irrigation_soil(ji), changed to be coherent with new variable dimension
          water2infilt(ji,jst) = water2infilt(ji,jst) + irrigation_soil(ji, jst) + reinfiltration_soil(ji) &
                - MIN(ae_ns(ji,jst),zero) - MIN(subsinksoil(ji),zero) + precisol_ns(ji,jst) &
                - temp(ji) 
       ENDDO       
             
       ! The water to evaporate from the sol surface is either the difference to what has to be infiltrated, or 0
       !   - subsinksoil is the residual from sublimation is the snowpack is not sufficient
       !   - how are the negative values of ae_ns taken into account ??? 
       DO ji = 1, kjpindex
          water2extract(ji) =  MAX(ae_ns(ji,jst),zero) + MAX(subsinksoil(ji),zero) - temp(ji) 
       ENDDO

       ! Here we acknowledge that subsinksoil is part of ae_ns, but ae_ns is not used further
       ae_ns(:,jst) = ae_ns(:,jst) + subsinksoil(:)  

       !! 1.2 To remove water2extract (including bare soil) from top layer
       flux_top(:) = water2extract(:)

       !! 1.3 Infiltration

       !! Definition of flux_infilt
       DO ji = 1, kjpindex
          ! Initialise the flux to be infiltrated  
          flux_infilt(ji) = water2infilt(ji,jst) 
       ENDDO
       
       !! K and D are computed for the profile of mc before infiltration
       !! They depend on the fraction of soil ice, given by profil_froz_hydro_ns
       CALL hydrol_soil_coef(mcr, mcs, kjpindex,jst, soiltile)

       !! Infiltration and surface runoff are computed
       !! Infiltration stems from comparing liquid water2infilt to initial total mc (liquid+ice)
       !! The conductivity comes from hydrol_soil_coef and relates to the liquid phase only
       !  This seems consistent with ok_freeze
       CALL hydrol_soil_infilt(ks, mcr, mcs, mcfc, mcw, kjpindex, jst, flux_infilt, stempdiag, &
            soiltile, qinfilt_ns, ru_infilt_ns, check_infilt_ns)
       ru_ns(:,jst) = ru_infilt_ns(:,jst)

       !! 1.4 Reinfiltration of surface runoff : compute temporary surface water and extract from runoff
       ! Evrything here is liquid
       ! RK: water2infilt is both a volume for future reinfiltration (in mm) and a correction term for surface runoff (in mm/dt_sechiba) 
       IF ( .NOT. doponds ) THEN ! this is the general case...
          DO ji = 1, kjpindex
             water2infilt(ji,jst) = reinf_slope_soil(ji,jst) * ru_ns(ji,jst)
          ENDDO
       ELSE
          DO ji = 1, kjpindex           
             water2infilt(ji,jst) = zero
          ENDDO
       ENDIF
       !
       DO ji = 1, kjpindex           
          ru_ns(ji,jst) = ru_ns(ji,jst) - water2infilt(ji,jst)
       END DO

       !! 2. SECONDLY, WE UPDATE MC FROM DIFFUSION, INCLUDING DRAINAGE AND ROOTSINK
       !!    This will act on mcl only
       
       !! 2.1 K and D are recomputed after infiltration
       !! They depend on the fraction of soil ice, still given by profil_froz_hydro_ns 
       CALL hydrol_soil_coef(mcr, mcs, kjpindex,jst, soiltile)

 
       !! 2.2 Set the tridiagonal matrix coefficients for the diffusion/redistribution scheme
       !! This process will further act on mcl only, based on a, b, d from hydrol_soil_coef
       CALL hydrol_soil_setup(kjpindex,jst)

       !! 2.3 We define mcl (liquid water content) based on mc and profil_froz_hydro_ns
       DO jsl = 1, nslm
          DO ji =1, kjpindex
             IF ( (ok_peat_hydro)  .AND. (jst .GE. 4) .AND. (jst .LE. 6)) THEN
                IF (soiltile (ji,jst) .GT. min_sechiba) THEN
                   mcl(ji,jsl,jst)= MIN( mc(ji,jsl,jst), mcr_peat + &
                        (un-profil_froz_hydro_ns(ji,jsl,jst))*(mc(ji,jsl,jst)-mcr_peat) )
                ELSE
                   mcl(ji,jsl,jst)= 0.0
                ENDIF
             ELSE 
                mcl(ji,jsl,jst)= MIN( mc(ji,jsl,jst), mcr(ji) + &
                     (un-profil_froz_hydro_ns(ji,jsl,jst))*(mc(ji,jsl,jst)-mcr(ji)) )
             ENDIF
             ! we always have mcl<=mc
             ! if mc>mcr, then mcl>mcr; if mc=mcr,mcl=mcr; if mc<mcr, then mcl<mcr
             ! if profil_froz_hydro_ns=0 (including NOT ok_freeze_cwrr) we keep mcl=mc
          ENDDO
       ENDDO

       ! The value of mcl is kept in mclint (nstm dimension removed), used in the flux computation after diffusion
       DO jsl = 1, nslm
          DO ji = 1, kjpindex
             mclint(ji,jsl) = mask_soiltile(ji,jst) * mcl(ji,jsl,jst)
          ENDDO
       ENDDO

       !! 2.3bis Diagnostic of the matric potential used for redistribution by Richards/tridiag (in m)
       !  We use VG relationship giving psi as a function of mc (mcl in our case)
       !  With patches against numerical pbs when (mc_ratio - un) becomes very slightly negative (gives NaN) 
       !  or if psi become too strongly negative (pbs with xios output)
       !
       DO jsl=1, nslm
          DO ji = 1, kjpindex
             IF (soiltile(ji,jst) .GT. zero) THEN                
                ! 
                mvg = un - un / nvan_mod_tab(jsl,ji)
                avg = avan_mod_tab(jsl,ji)*1000. ! to convert in m-1
                mc_ratio = MAX( 10.**(-14*mvg), (mcl(ji,jsl,jst) - mcr(ji))/(mcs(ji) - mcr(ji)) )**(-un/mvg)
                psi = - MAX(zero,(mc_ratio - un))**(un/nvan_mod_tab(jsl,ji)) / avg ! in m 
                psi_moy(ji,jsl) = psi_moy(ji,jsl) + soiltile(ji,jst) * psi ! average across soil tiles

                IF (ok_peat_hydro .AND. (jst .GE. 4) .AND. (jst .LE. 6)) THEN
                   mvg_peat = un - un / nvan_peat
                   avg_peat = avan_peat*1000. ! to convert in m-1
                   nfact_peat(jsl)=( kfact_peat(jsl) )**nk_rel
                   nvan_mod_peat = n0 + (nvan_peat-n0) * nfact_peat(jsl)
                   
                   mc_ratio_peat = MAX( 10.**(-14*mvg_peat), (mcl(ji,jsl,jst) - mcr_peat)/(mcs_peat - mcr_peat) )**(-un/mvg_peat)
                   psi_peat = - MAX(zero,(mc_ratio_peat - un))**(un/nvan_mod_peat) / avg_peat ! in m 
                ELSE
                   mc_ratio = MAX( 10.**(-14*mvg), (mcl(ji,jsl,jst) - mcr(ji))/(mcs(ji) - mcr(ji)) )**(-un/mvg)
                   psi_peat = - MAX(zero,(mc_ratio - un))**(un/nvan_mod_tab(jsl,ji)) / avg
                ENDIF
                psi_moy(ji,jsl) = psi_moy(ji,jsl) + soiltile(ji,jst) * psi_peat ! average across soil tiles
             ENDIF
          ENDDO
       ENDDO

       !! 2.4 We calculate the total SM at the beginning of the routine tridiag for water conservation check
       !  (on mcl only, since the diffusion only modifies mcl)
       tmci(:) = dz(2) * ( trois*mcl(:,1,jst) + mcl(:,2,jst) )/huit
       DO jsl = 2,nslm-1
          tmci(:) = tmci(:) + dz(jsl) * (trois*mcl(:,jsl,jst)+mcl(:,jsl-1,jst))/huit &
               + dz(jsl+1) * (trois*mcl(:,jsl,jst)+mcl(:,jsl+1,jst))/huit
       ENDDO
       tmci(:) = tmci(:) + dz(nslm) * (trois*mcl(:,nslm,jst) + mcl(:,nslm-1,jst))/huit

       ! peatland 
       IF (ok_freeze_cwrr) THEN
          DO ji =1, kjpindex
             DO jsl = 1, nslm
                IF ( (ok_peat_hydro) .AND. (jst .GE. 4) .AND. (jst .LE. 6) ) THEN
                   IF (soiltile (ji,jst) .GT. min_sechiba) THEN
                      mcl(ji,jsl,jst)= MIN(mc(ji,jsl,jst),mcr_peat+(1-profil_froz_hydro_ns(ji,jsl,jst))*(mc(ji,jsl,jst)-mcr_peat))
                   ELSE
                      mcl(ji,jsl,jst)= 0.0
                   ENDIF
                ELSE 
                   mcl(ji,jsl,jst)= MIN(mc(ji,jsl,jst),mcr(ji)+(1-profil_froz_hydro_ns(ji,jsl,jst))*(mc(ji,jsl,jst)-mcr(ji)))
                ENDIF
             ENDDO
          ENDDO
       ELSE
          mcl(:,:,jst)=mc(:,:,jst)
       ENDIF
       

       !! 2.5 Defining where diffusion is solved : everywhere
       !! Since mc>mcs is not possible after infiltration, and we accept that mc<mcr
       !! (corrected later by shutting off all evaporative fluxes in this case)
       !  Nothing is done if resolv=F
       !
       ! Given a soiltile jst,resolv is a 1D array with logical values:
       ! given the interested soil tile jst and a grid, if the soiltile of this jst > 0ru
       ! the resol of this grid is True.
       resolv(:) = (mask_soiltile(:,jst) .GT. 0)

       !! 2.6 We define the system of linear equations for mcl redistribution,
       !! based on the matrix coefficients from hydrol_soil_setup
       !! following the PhD thesis of de Rosnay (1999), p155-157
       !! The bare soil evaporation (subtracted from infiltration) is used directly as flux_top
       ! rhs for right-hand side term; fp for f'; gp for g'; ep for e'; with flux=0 !
      
       !- First layer
       DO ji = 1, kjpindex
          tmat(ji,1,1) = zero
          tmat(ji,1,2) = f(ji,1)
          tmat(ji,1,3) = g1(ji,1)
          rhs(ji,1)    = fp(ji,1) * mcl(ji,1,jst) + gp(ji,1)*mcl(ji,2,jst) &
               &  - flux_top(ji) - (b(ji,1)+b(ji,2))/deux *(dt_sechiba/one_day) - rootsink(ji,1,jst)
       ENDDO
       !- soil body
       DO jsl=2, nslm-1
          DO ji = 1, kjpindex
             tmat(ji,jsl,1) = e(ji,jsl)
             tmat(ji,jsl,2) = f(ji,jsl)
             tmat(ji,jsl,3) = g1(ji,jsl)
             rhs(ji,jsl) = ep(ji,jsl)*mcl(ji,jsl-1,jst) + fp(ji,jsl)*mcl(ji,jsl,jst) &
                  & +  gp(ji,jsl) * mcl(ji,jsl+1,jst) & 
                  & + (b(ji,jsl-1) - b(ji,jsl+1)) * (dt_sechiba/one_day) / deux & 
                  & - rootsink(ji,jsl,jst) 
          ENDDO
       ENDDO       
       !- Last layer, including drainage
       DO ji = 1, kjpindex
          jsl=nslm
          tmat(ji,jsl,1) = e(ji,jsl)
          tmat(ji,jsl,2) = f(ji,jsl)
          tmat(ji,jsl,3) = zero
          rhs(ji,jsl) = ep(ji,jsl)*mcl(ji,jsl-1,jst) + fp(ji,jsl)*mcl(ji,jsl,jst) &
               & + (b(ji,jsl-1) + b(ji,jsl)*(un-deux*free_drain_coef(ji,jst))) * (dt_sechiba/one_day) / deux &
               & - rootsink(ji,jsl,jst)
       ENDDO
       !- Store the equations in case needed again
       DO jsl=1,nslm
          DO ji = 1, kjpindex
             srhs(ji,jsl) = rhs(ji,jsl)
             stmat(ji,jsl,1) = tmat(ji,jsl,1)
             stmat(ji,jsl,2) = tmat(ji,jsl,2)
             stmat(ji,jsl,3) = tmat(ji,jsl,3) 
          ENDDO
       ENDDO
       
       !! 2.7 Solves diffusion equations, but only in grid-cells where resolv is true, i.e. everywhere (cf 2.2)
       !!     The result is an updated mcl profile

       CALL hydrol_soil_tridiag(kjpindex,jst)

       !! 2.8 Computes drainage = bottom boundary condition, consistent with rhs(ji,jsl=nslm)
       ! dr_ns in mm/dt_sechiba, from k in mm/d
       ! This should be done where resolv=T, like tridiag (drainage is part of the linear system !)
       DO ji = 1, kjpindex
          IF (resolv(ji)) THEN
             ! if a soiltile fraction is > 0 at jst, mask_soiltile is 1, resolv for this grid is T
             ! dr_ns is computed according to soil conductivite and free_drain_coef.
             ! free_drain_coef is imporant.
             dr_ns(ji,jst) = mask_soiltile(ji,jst)*k(ji,nslm)*free_drain_coef(ji,jst) * (dt_sechiba/one_day)
          ELSE
             !if a soiltile fraction is 0 at jst, mask_soiltile is 0, resolv for this grid is F
             ! dr_ns for this jst is 0. 
             dr_ns(ji,jst) = zero
          ENDIF
       ENDDO

       !! 2.9 For water conservation check during redistribution AND CORRECTION, 
       !!     we calculate the total liquid SM at the end of the routine tridiag
       tmcf(:) = dz(2) * ( trois*mcl(:,1,jst) + mcl(:,2,jst) )/huit
       DO jsl = 2,nslm-1
          tmcf(:) = tmcf(:) + dz(jsl) * (trois*mcl(:,jsl,jst)+mcl(:,jsl-1,jst))/huit &
               + dz(jsl+1) * (trois*mcl(:,jsl,jst)+mcl(:,jsl+1,jst))/huit
       ENDDO
       tmcf(:) = tmcf(:) + dz(nslm) * (trois*mcl(:,nslm,jst) + mcl(:,nslm-1,jst))/huit
          
       !! And we compare the difference with the flux...
       ! Normally, tcmf=tmci-flux_top(ji)-transpir-dr_ns
       DO ji=1,kjpindex
          diag_tr(ji)=SUM(rootsink(ji,:,jst))
       ENDDO
       ! Here, check_tr_ns holds the inaccuracy during the redistribution phase
       check_tr_ns(:,jst) = tmcf(:)-(tmci(:)-flux_top(:)-dr_ns(:,jst)-diag_tr(:))

       !! We solve here the numerical errors that happen when the soil is close to saturation
       !! and drainage very high, and which lead to negative check_tr_ns: the soil dries more 
       !! than what is demanded by the fluxes, so we need to increase the fluxes.
       !! This is done by increasing the drainage.
       !! There are also instances of positive check_tr_ns, larger when the drainage is high
       !! They are similarly corrected by a decrease of dr_ns, in the limit of keeping a positive drainage.
       DO ji=1,kjpindex
          IF ( check_tr_ns(ji,jst) .LT. zero ) THEN
              dr_corrnum_ns(ji,jst) = -check_tr_ns(ji,jst)
          ELSE
              dr_corrnum_ns(ji,jst) = -MIN(dr_ns(ji,jst),check_tr_ns(ji,jst))             
          ENDIF
          dr_ns(ji,jst) = dr_ns(ji,jst) + dr_corrnum_ns(ji,jst) ! dr_ns increases/decrease if check_tr negative/positive
       ENDDO
       !! For water conservation check during redistribution
       IF (check_cwrr) THEN          
          check_tr_ns(:,jst) = tmcf(:)-(tmci(:)-flux_top(:)-dr_ns(:,jst)-diag_tr(:)) 
       ENDIF

       !! 3. AFTER DIFFUSION/REDISTRIBUTION

       !! 3.1 Updating mc, as all the following checks against saturation will compare mc to mcs
       !      The frozen fraction is constant, so that any water flux to/from a layer changes
       !      both mcl and the ice amount. The assumption behind this is that water entering/leaving
       !      a soil layer immediately freezes/melts with the proportion profil_froz_hydro_ns/(1-profil_...)
       DO jsl = 1, nslm
          DO ji =1, kjpindex
             IF ( (ok_peat_hydro)  .AND. (jst .GE. 4) .AND. (jst .LE. 6)) THEN
                IF (soiltile (ji,jst) .GT. min_sechiba) THEN 
                   mc(ji,jsl,jst) = MAX( mcl(ji,jsl,jst), mcl(ji,jsl,jst) + &
                        profil_froz_hydro_ns(ji,jsl,jst)*(mc(ji,jsl,jst)-mcr_peat) )
                ELSE
                   mc(ji,jsl,jst) = 0.0
                ENDIF
             ELSE
                mc(ji,jsl,jst) = MAX( mcl(ji,jsl,jst), mcl(ji,jsl,jst) + &
                     profil_froz_hydro_ns(ji,jsl,jst)*(mc(ji,jsl,jst)-mcr(ji)) )
             ENDIF
             ! if profil_froz_hydro_ns=0 (including NOT ok_freeze_cwrr) we get mc=mcl
          ENDDO
       ENDDO

       !! 3.2 Correct here the possible over-saturation values (subroutine hydrol_soil_smooth_over_mcs2 acts on mc)
       !    Oversaturation results from numerical inaccuracies and can be frequent if free_drain_coef=0
       !    Here hydrol_soil_smooth_over_mcs2 discard all excess as ru_corr_ns, oriented to either ru_ns or dr_ns
       dr_corr_ns(:,jst) = zero
       ru_corr_ns(:,jst) = zero
       call hydrol_soil_smooth_over_mcs2(mcs, kjpindex, jst, soiltile, ru_corr_ns, check_over_ns)
       
       ! In absence of freezing, if F is large enough, the correction of oversaturation is sent to drainage       
       DO ji = 1, kjpindex
          IF ((free_drain_coef(ji,jst) .GE. 0.5) .AND. (.NOT. ok_freeze_cwrr) ) THEN
             dr_corr_ns(ji,jst) = ru_corr_ns(ji,jst) 
             ru_corr_ns(ji,jst) = zero
          ENDIF
       ENDDO
       dr_ns(:,jst) = dr_ns(:,jst) + dr_corr_ns(:,jst)
       ru_ns(:,jst) = ru_ns(:,jst) + ru_corr_ns(:,jst)

       !! 3.3 Negative runoff is reported to drainage
       !  Since we computed ru_ns directly from hydrol_soil_infilt, ru_ns should not be negative
       ru_corr2_ns(:,jst) = zero
       DO ji = 1, kjpindex
          IF (ru_ns(ji,jst) .LT. zero) THEN
             IF (printlev>=3)  WRITE (numout,*) 'NEGATIVE RU_NS: runoff and drainage before correction',&
                  ru_ns(ji,jst),dr_ns(ji,jst)
             dr_ns(ji,jst)=dr_ns(ji,jst)+ru_ns(ji,jst)
             ru_corr2_ns(ji,jst) = -ru_ns(ji,jst)
             ru_ns(ji,jst)= 0.
          END IF          
       ENDDO

       !! 3.4.1 Optional nudging for soil moisture 
       IF (ok_nudge_mc) THEN
          CALL hydrol_nudge_mc(kjpindex, jst, mc)
       END IF


       !! 3.4.2 Optional block to force saturation below zwt_force
       ! This block is not compatible with freezing; in this case, mcl must be corrected too
       ! We test if zwt_force(1,jst) <= zmaxh, to avoid steps 1 and 2 if unnecessary
       
       IF (zwt_force(1,jst) <= zmaxh) THEN

          !! We force the nodes below zwt_force to be saturated
          !  As above, we compare mc to mcs 
          DO jsl = 1,nslm
             DO ji = 1, kjpindex
                dmc(ji,jsl) = zero
                IF ( ( zz(jsl) >= zwt_force(ji,jst)*mille ) ) THEN
                   !! note that we are still in the loop of nstm (soil tiles).
                   IF ( (ok_peat_hydro) .AND. (jst .GE. 4) .AND. (jst .LE. 6)) THEN
                      IF (soiltile (ji,jst) .GT. min_sechiba) THEN
                         dmc(ji,jsl) = mcs_peat - mc(ji,jsl,jst)
                         mc(ji,jsl,jst) = mcs_peat
                      ELSE
                         dmc(ji,jsl) = mcs(ji) - mc(ji,jsl,jst) ! addition to reach mcs (m3/m3) = positive value
                         mc(ji,jsl,jst) = mcs(ji) ! should not set to 0 
                      ENDIF
                   ELSE
                      ! orig trunk
                      dmc(ji,jsl) = mcs(ji) - mc(ji,jsl,jst) ! addition to reach mcs (m3/m3) = positive value
                      mc(ji,jsl,jst) = mcs(ji)
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
          
          !! To ensure conservation, this needs to be balanced by a negative change in drainage (in kg/m2/dt)
          DO ji = 1, kjpindex
             dr_force_ns(ji,jst) = dz(2) * ( trois*dmc(ji,1) + dmc(ji,2) )/huit ! top layer = initialization
          ENDDO
          DO jsl = 2,nslm-1 ! intermediate layers
             DO ji = 1, kjpindex
                dr_force_ns(ji,jst) = dr_force_ns(ji,jst) + dz(jsl) &
                     & * (trois*dmc(ji,jsl)+dmc(ji,jsl-1))/huit &
                     & + dz(jsl+1) * (trois*dmc(ji,jsl)+dmc(ji,jsl+1))/huit
             ENDDO
          ENDDO
          DO ji = 1, kjpindex
             dr_force_ns(ji,jst) = dr_force_ns(ji,jst) + dz(nslm) & ! bottom layer
                  & * (trois * dmc(ji,nslm) + dmc(ji,nslm-1))/huit
             dr_ns(ji,jst) = dr_ns(ji,jst) - dr_force_ns(ji,jst) ! dr_force_ns is positive and dr_ns must be reduced
          END DO

       ELSE          

          dr_force_ns(:,jst) = zero 

       ENDIF

       !! 3.5 Diagnosing the effective water table depth:
       !!     Defined as as the smallest jsl value when mc(jsl) is no more at saturation (mcs), starting from the bottom
       !      If there is a part of the soil which is saturated but underlain with unsaturated nodes,
       !      this is not considered as a water table
       DO ji = 1, kjpindex
          wtd_ns(ji,jst) = undef_sechiba ! effective water table depth in meters
          jsl=nslm
          IF ( (ok_peat_hydro) .AND. (jst .EQ. 4) .AND. (jst .LE. 6)) THEN
             IF (soiltile (ji,jst) .GT. min_sechiba) THEN
                DO WHILE ( (mc(ji,jsl,jst) .EQ. mcs_peat) .AND. (jsl > 1) )
                   wtd_ns(ji,jst) = zz(jsl)/mille ! in meters 
                   jsl=jsl-1
                ENDDO

             ELSE
                DO WHILE ( (mc(ji,jsl,jst) .EQ. mcs(ji)) .AND. (jsl > 1) )
                   wtd_ns(ji,jst) = zz(jsl)/mille ! in meters
                   jsl=jsl-1
                ENDDO
             ENDIF
          ELSE
             DO WHILE ( (mc(ji,jsl,jst) .EQ. mcs(ji)) .AND. (jsl > 1) )
                wtd_ns(ji,jst) = zz(jsl)/mille ! in meters
                jsl=jsl-1
             ENDDO
          ENDIF
       ENDDO

       !! 3.6 Diagnose under_mcr to adapt water stress calculation below
       !      This routine does not change tmc but decides where we should turn off ET to prevent further mc decrease
       !      Like above, the tests are made on total mc, compared to mcr
       CALL hydrol_soil_smooth_under_mcr(mcr, kjpindex, jst, soiltile, is_under_mcr, check_under_ns)

       ! The tides was made by VN, and is included in the orchidee trunk by xw, which is not tested yet. 
       ! agri_peat is included, and not tested yet.
       ! mangrove: low relief, runoff from soiltile1-3 are leaded into soiltile 6 according to soil tile area
       IF (tides) THEN 
          DO ji = 1, kjpindex
             ! ok_run2peat is true if ok_peat_hydro or tides are activated.
             IF ( (ok_peat_hydro) .AND. (ok_run2peat) .AND. ((jst .GE. 4) .AND. (jst .LE. 6))) THEN
                IF ((jst .EQ. 4) .AND. (soiltile(ji,4) .GT. min_sechiba) ) THEN
                   run2peat(ji)=zero
                   ! change water added into pealtand according to the peatland area
                   run2peat(ji)= ( soiltile(ji,4) / (soiltile(ji,4) + soiltile(ji,6)) ) & 
                        * ( ru_ns(ji,1)*soiltile(ji,1) & 
                        + ru_ns(ji,2)*soiltile(ji,2) & 
                        + ru_ns(ji,3)*soiltile(ji,3) )
                   IF ( (soiltile(ji,6) .LT. min_sechiba) ) THEN
                      ru_ns(ji,1)=zero
                      ru_ns(ji,2)=zero
                      ru_ns(ji,3)=zero
                   ENDIF
                ELSEIF ((jst .EQ. 6) .AND. (soiltile(ji,6) .GT. min_sechiba) ) THEN
                   run2man(ji)=zero
                   run2man(ji)=(soiltile(ji,6)/(soiltile(ji,4)+soiltile(ji,6))) &
                        *(ru_ns(ji,1)*soiltile(ji,1)+ru_ns(ji,2)*soiltile(ji,2)+ru_ns(ji,3)*soiltile(ji,3))
                   ru_ns(ji,1)=zero
                   ru_ns(ji,2)=zero
                   ru_ns(ji,3)=zero
                ENDIF
             ENDIF
          ENDDO
       ELSE IF ((.NOT. tides) .AND. (ok_peat_hydro) .AND. (ok_run2peat)) THEN
          ! if peatland activated, and tides not activiated: low relief, runoff from soiltile1-3 are leaded into peatland
          ! and will be infiltrated into the soil in the next time step
          ! the ru_ns of non-peatland is set to 0.
          ! thus this will change the results of runoff_pft(ji, jvm) (which correspond to ru_ns(ji, jst)), esp for non-peat pft
          ! for runoff_pft(ji,jvm) of peatland, it will be also modified in the next time step due to infiltration
          DO ji = 1, kjpindex
             IF ((jst .EQ. 4)) THEN
                IF ((soiltile(ji,jst) .GT. min_sechiba)) THEN
                   ! this is the orginal code from mict-peat
                   ! the condition soiltile here would cancel the possible effects of free_drain_coef when pft peat is 0 ??
                   run2peat(ji)=zero
                   run2peat(ji)= ru_ns(ji,1)*soiltile(ji,1)+ru_ns(ji,2)*soiltile(ji,2)+ru_ns(ji,3)*soiltile(ji,3)
                   ru_ns(ji,1)=zero
                   ru_ns(ji,2)=zero
                   ru_ns(ji,3)=zero
                ELSE
                   run2peat(ji)=zero
                   ! do not modify ru_ns of other soil tiles in this case
                ENDIF
             ENDIF
          ENDDO
       ENDIF
   
       ! peatland: water table can sustain above soil surface, part of runoff
       !      of soiltile4 is conserved by above surface resevoir wt_ab 
       ! ok_wt_ab is true if ok_peat_hydro or tides are activiated
       DO ji = 1, kjpindex
          ! ok_wt_ab is true when ok_peat_hydro and ok_wt_ab are T
          IF ( (ok_peat_hydro) .AND. (ok_wt_ab) )THEN
             ! the soiltile 6 for tides is not considered in the condition below
             ! only the soiltile 4.
             IF ((jst .EQ. 4) .AND. (soiltile(ji,4) .GT. min_sechiba) ) THEN
                wt_ab(ji)=zero
                ! the runoff of soiltile 4 is used to compute wt_ab
                ! without adding run2peat. Is this normal ? 
                wt_ab(ji)=wt_ab(ji) + ru_ns(ji,jst)
                ru_ns(ji,jst)=MAX(wt_ab(ji) - max_wt_ab,zero)
                wt_ab(ji)=MIN(max_wt_ab, wt_ab(ji))
                !! water of above surface reservoir will be infiltrated into the
                !   soil in the next time step
                !         water2infilt(ji,jst)=water2infilt(ji,jst)+ wt_ab(ji)
             ELSE
                wt_ab(ji)=zero
             ENDIF
          ENDIF
       ENDDO

       ! tides: water table can sustain above soil surface, tidal water is conserved by above-surface reservoir wt_ab_tide   
       !        If wt_ab_tide decreases compared to previous timestep, you have a saltwater runoff for soil tile 4
       !        If wt_ab_tide increases there is no runoff for soil tile 4; if wt_ab is 0, runoff = runoff (no runoff_tide)
       IF ( tides ) THEN
          DO ji = 1, kjpindex
    	     IF ( jst .EQ. 4 ) THEN
                IF ( (wtp_tide .GT. zero) .AND. (dwtp_tide .GT. zero) ) THEN !! tide is positive (above-ground) and increasing 
                   wt_ab_tide(ji)=wtp_tide+ru_ns(ji,jst)
                   ru_ns(ji,jst)=zero
                   water2infilt(ji,jst)=water2infilt(ji,jst)+ wt_ab_tide(ji)
                ELSE IF ( (wtp_tide .LE. zero) .AND. (dwtp_tide .LE. zero) ) THEN !! no tide, wt_ab_tide is empty
                   wt_ab_tide(ji) = zero 
                   ru_ns(ji,jst)=ru_ns(ji,jst) ! runoff is non-infiltrated excess precip water
                ELSE IF ( (wtp_tide .LE. zero) .AND. (dwtp_tide .LT. zero) ) THEN !! tide decreasing, reaching zero at this timestep
                   ru_ns(ji,jst) = wt_ab_tide(ji) ! runoff is equal to non-infiltrated water + tidal retreat (what is left over in wt_ab_tide reservoir)
                   wt_ab_tide(ji) = zero 	
                ELSE IF ( (wtp_tide .GT. zero) .AND. (dwtp_tide .LT. zero) ) THEN !! tide decreasing
                   ru_ns(ji,jst)=ru_ns(ji,jst)-dwtp_tide  !!! runoff is equal to non-infiltrated water + tidal retreat
                   wt_ab_tide(ji) = wtp_tide
                   water2infilt(ji,jst)=water2infilt(ji,jst)+ wt_ab_tide(ji)	
                ENDIF
             ENDIF
          ENDDO
       ENDIF
      

       !! 4. At the end of the prognostic calculations, we recompute important moisture variables

       !! 4.1 Total soil moisture content (water2infilt added below)
       DO ji = 1, kjpindex
          tmc(ji,jst) = dz(2) * ( trois*mc(ji,1,jst) + mc(ji,2,jst) )/huit
       ENDDO
       DO jsl = 2,nslm-1
          DO ji = 1, kjpindex
             tmc(ji,jst) = tmc(ji,jst) + dz(jsl) &
                  & * (trois*mc(ji,jsl,jst)+mc(ji,jsl-1,jst))/huit &
                  & + dz(jsl+1) * (trois*mc(ji,jsl,jst)+mc(ji,jsl+1,jst))/huit
          ENDDO
       ENDDO
       DO ji = 1, kjpindex
          tmc(ji,jst) = tmc(ji,jst) + dz(nslm) &
               & * (trois * mc(ji,nslm,jst) + mc(ji,nslm-1,jst))/huit
       END DO

       !! 4.2 mcl is a module variable; we update it here for calculating bare soil evaporation,
       !!     and in case we would like to export it (xios)
       DO jsl = 1, nslm
          DO ji =1, kjpindex
             IF ( (ok_peat_hydro) .AND. (jst .GE. 4) .AND. (jst .LE. 6)) THEN
                IF (soiltile(ji,jst) .GT. min_sechiba) THEN
                   mcl(ji,jsl,jst)= MIN( mc(ji,jsl,jst), mcr_peat + &
                        (un-profil_froz_hydro_ns(ji,jsl,jst))*(mc(ji,jsl,jst)-mcr_peat) )
                ELSE
                   mcl(ji,jsl,jst)= 0.0
                ENDIF
             ELSE
                mcl(ji,jsl,jst)= MIN( mc(ji,jsl,jst), mcr(ji) + &
                     (un-profil_froz_hydro_ns(ji,jsl,jst))*(mc(ji,jsl,jst)-mcr(ji)) )
                ! if profil_froz_hydro_ns=0 (including NOT ok_freeze_cwrr) we keep mcl=mc
             ENDIF
          ENDDO
       ENDDO
       
       !! 5. Optional check of the water balance of soil column (if check_cwrr)

       IF (check_cwrr) THEN

          !! 5.1 Computation of the vertical water fluxes and water balance of the top layer
          CALL hydrol_diag_soil_flux(kjpindex,jst,mclint,flux_top)

       ENDIF

       !! 6. SM DIAGNOSTICS FOR OTHER ROUTINES, MODULES, OR NEXT STEP
       !    Starting here, mc and mcl should not change anymore 
       
       !! 6.1 Total soil moisture, soil moisture at litter levels, soil wetness, us, humrelv, vesgtressv
       !!     (based on mc)

       !! In classical output, tmc includes water2infilt(ji,jst)
       ! If peat is activated:
       !    tmc computed here will be given to define tmc_soil in the later part of code
       !    tmc will be then updated by adding water2infilt which also takes into account the runoff
       ! So if peatland is activated, tmc will be larger without peatland, due to runoff from other tiles included
       ! Adding the ok_peat_hydro condition is good to avoid having water2infilt for twice.
       !
       ! IF peatland is not activated:
       !     tmc computed here will represent tmc plus water2infilt
       ! 
       IF (.NOT. ok_peat_hydro) THEN
          DO ji=1,kjpindex
             tmc(ji,jst) = tmc(ji,jst) + water2infilt(ji,jst)
          END DO
       ENDIF
              
       ! The litter is the 4 top levels of the soil
       ! Compute various field of soil moisture for the litter (used for stomate and for albedo)
       ! We exclude the frozen water from the calculation
       DO ji=1,kjpindex
          tmc_litter(ji,jst) = dz(2) * ( trois*mcl(ji,1,jst)+ mcl(ji,2,jst))/huit
       END DO
       ! sum from level 1 to 4
       DO jsl=2,4
          DO ji=1,kjpindex
             tmc_litter(ji,jst) = tmc_litter(ji,jst) + dz(jsl) * & 
                  & ( trois*mcl(ji,jsl,jst) + mcl(ji,jsl-1,jst))/huit &
                  & + dz(jsl+1)*(trois*mcl(ji,jsl,jst) + mcl(ji,jsl+1,jst))/huit
          END DO
       END DO

       !! peatland, tides
       ! According to the comments of Chunjing:
       ! mc_peat_above is used in littercalc, use water content instead of relative hum for peat.
       ! use the first 4 layers only.
       ! But this comment is not correct, because mc_peat_above is now relative humidity
       !
       mc_man_above(:)=0
       mc_peat_above(:)=0
       IF (perma_peat) THEN
          DO ji=1,kjpindex
             IF ((soiltile (ji,4) .GT. min_sechiba)) THEN 
                mc_peat_above(ji)= tmc_litter(ji,4)/(dh(1)+dh(2)+dh(3)+dh(4))
                mc_croppeat_above(ji)=0
                IF (tides) THEN
                   mc_man_above(ji)= tmc_litter(ji,6)/(dh(1)+dh(2)+dh(3)+dh(4))
                ENDIF
                IF (agri_peat) THEN
                   mc_croppeat_above(ji)= tmc_litter(ji,5)/(dh(1)+dh(2)+dh(3)+dh(4))
                ENDIF
             ENDIF
          ENDDO
       ENDIF
       !! 202603 mc_peat_above is relative humidity 0-1
       !! but we did not apply min max as other variables !!
       !! to be verified

       ! Subsequent calculation of soil_wet_litter (tmc-tmcw)/(tmcfc-tmcw)
       ! Based on liquid water content
       DO ji=1,kjpindex
          IF (soiltile (ji,jst) .GT. min_sechiba) THEN 
             soil_wet_litter(ji,jst) = MIN(un, MAX(zero,&
                  & (tmc_litter(ji,jst)-tmc_litter_wilt(ji,jst)) / &
                  & (tmc_litter_field(ji,jst)-tmc_litter_wilt(ji,jst)) ))
          ELSE
             soil_wet_litter(ji,jst) = 0.0
          ENDIF
       END DO
       
       ! Preliminary calculation of various soil moistures (for each layer, in kg/m2)
       sm(:,1,jst)  = dz(2) * (trois*mcl(:,1,jst) + mcl(:,2,jst))/huit
       smt(:,1,jst) = dz(2) * (trois*mc(:,1,jst) + mc(:,2,jst))/huit
       ! smw, smf, sms: dimension kjpindex,nslm
       smw(:,1) = dz(2) * (quatre*mcw(:))/huit
       smf(:,1) = dz(2) * (quatre*mcfc(:))/huit
       sms(:,1) = dz(2) * (quatre*mcs(:))/huit

       IF ( (ok_peat_hydro) .AND. (jst .GE. 4) .AND. (jst .LE. 6) ) THEN
          DO ji = 1, kjpindex
             IF (soiltile(ji,jst) .GT. min_sechiba) THEN
                smw(ji,1) = dz(2) * (quatre*mcw_peat)/huit
                smf(ji,1) = dz(2) * (quatre*mcf_peat)/huit
                sms(ji,1) = dz(2) * (quatre*mcs_peat)/huit
             ENDIF
          ENDDO
       ELSE
          ! ! original ORC trunk
          smw(:,1) = dz(2) * (quatre*mcw(:))/huit
          smf(:,1) = dz(2) * (quatre*mcfc(:))/huit
          sms(:,1) = dz(2) * (quatre*mcs(:))/huit
       ENDIF

       DO jsl = 2,nslm-1
          sm(:,jsl,jst)  = dz(jsl) * (trois*mcl(:,jsl,jst)+mcl(:,jsl-1,jst))/huit &
               + dz(jsl+1) * (trois*mcl(:,jsl,jst)+mcl(:,jsl+1,jst))/huit
          smt(:,jsl,jst) = dz(jsl) * (trois*mc(:,jsl,jst)+mc(:,jsl-1,jst))/huit &
               + dz(jsl+1) * (trois*mc(:,jsl,jst)+mc(:,jsl+1,jst))/huit
          smw(:,jsl) = dz(jsl) * ( quatre*mcw(:) )/huit &
               + dz(jsl+1) * ( quatre*mcw(:) )/huit
          smf(:,jsl) = dz(jsl) * ( quatre*mcfc(:) )/huit &
               + dz(jsl+1) * ( quatre*mcfc(:) )/huit
          sms(:,jsl) = dz(jsl) * ( quatre*mcs(:) )/huit &
               + dz(jsl+1) * ( quatre*mcs(:) )/huit
          IF ( (ok_peat_hydro) .AND. (jst .GE. 4) .AND. (jst .LE. 6)) THEN
             DO ji = 1, kjpindex
                IF (soiltile(ji,jst) .GT. min_sechiba) THEN
                   smw(ji,jsl) = dz(jsl) * ( quatre*mcw_peat )/huit &
                        + dz(jsl+1) * ( quatre*mcw_peat)/huit
                   smf(ji,jsl) = dz(jsl) * ( quatre*mcf_peat )/huit &
                        + dz(jsl+1) * ( quatre*mcf_peat )/huit
                   sms(ji,jsl) = dz(jsl) * ( quatre*mcs_peat )/huit &
                        + dz(jsl+1) * ( quatre*mcs_peat )/huit
                ENDIF
             ENDDO
          ELSE
          !    ! orig orc trunk
             smw(:,jsl) = dz(jsl) * ( quatre*mcw(:) )/huit &
                  + dz(jsl+1) * ( quatre*mcw(:) )/huit
             smf(:,jsl) = dz(jsl) * ( quatre*mcfc(:) )/huit &
                  + dz(jsl+1) * ( quatre*mcfc(:) )/huit
             sms(:,jsl) = dz(jsl) * ( quatre*mcs(:) )/huit &
                  + dz(jsl+1) * ( quatre*mcs(:) )/huit
          ENDIF
       ENDDO

       sm(:,nslm,jst)  = dz(nslm) * (trois*mcl(:,nslm,jst) + mcl(:,nslm-1,jst))/huit
       smt(:,nslm,jst) = dz(nslm) * (trois*mc(:,nslm,jst) + mc(:,nslm-1,jst))/huit

       !!!peatland
       smw(:,nslm) = dz(nslm) * (quatre*mcw(:))/huit
       smf(:,nslm) = dz(nslm) * (quatre*mcfc(:))/huit
       sms(:,nslm) = dz(nslm) * (quatre*mcs(:))/huit
       IF ( (ok_peat_hydro) .AND. (jst .GE. 4) .AND. (jst .LE. 6)  ) THEN
          DO ji = 1, kjpindex
             IF (soiltile(ji,jst) .GT. min_sechiba) THEN
                smw(ji,nslm) = dz(nslm) * (quatre*mcw_peat)/huit
                smf(ji,nslm) = dz(nslm) * (quatre*mcf_peat)/huit
                sms(ji,nslm) = dz(nslm) * (quatre*mcs_peat)/huit
             ENDIF
          ENDDO 
       ELSE
          smw(:,nslm) = dz(nslm) * (quatre*mcw(:))/huit
          smf(:,nslm) = dz(nslm) * (quatre*mcfc(:))/huit
          sms(:,nslm) = dz(nslm) * (quatre*mcs(:))/huit
       ENDIF
       !
       ! sm_nostress = soil moisture of each layer at which us reaches 1, here at the middle of [smw,smf]
       DO jsl = 1,nslm
          sm_nostress(:,jsl) = smw(:,jsl) + pcent(njsc(:)) * (smf(:,jsl)-smw(:,jsl))
          IF ( (ok_peat_hydro) .AND. (jst .GE. 4) .AND. (jst .LE. 6)  ) THEN
             DO ji = 1, kjpindex
                IF (soiltile(ji,jst) .GT. min_sechiba) THEN
                   sm_nostress(ji,jsl) = smw(ji,jsl) + pcent_peat * (smf(ji,jsl)-smw(ji,jsl))
                ENDIF
             ENDDO
          ELSE
             sm_nostress(:,jsl) = smw(:,jsl) + pcent(njsc(:)) * (smf(:,jsl)-smw(:,jsl))
          ENDIF
       ENDDO
       
       ! sm_nostress = soil moisture of each layer at which us reaches 1, here at the middle of [smw,smf]
       ! pcent_peat is equal to pcent in the current version, but since smw and swf are computed different for peatland and
       ! non peatland, so there is sm_nostress would be different as well.
   
       ! Saturated litter soil moisture for rsoil
       tmcs_litter(:) = zero
       DO jsl = 1,4
          tmcs_litter(:) = tmcs_litter(:) + sms(:,jsl)
       END DO

       ! Here we compute root zone deficit, to have an estimate of water demand in irrigated soil column (i.e. crop and grass)
       IF(jst .EQ. irrig_st ) THEN
           !It computes water deficit only on the root zone, and only on the layers where
           !there is actually a deficit. If there is not deficit, it does not take into account that layer
            DO ji = 1,kjpindex

                root_deficit(ji) = SUM( MAX(zero, beta_irrig*smf(ji,1:nslm_root(ji) ) &
                - sm(ji,1:nslm_root(ji),jst  ) )) - water2infilt(ji,jst)

                root_deficit(ji) = MAX( root_deficit(ji) ,  zero)

            ENDDO
            !It COUNTS the number of pft with LAI > lai_irrig_min, inside the soiltile
            !It compares veget, but it is the same as they are related by a function
            lai_irrig_trig(:) = 0

            DO jv = 1, nvm
              IF( is_irrigated(jv) ) THEN
                DO ji = 1,kjpindex

                  IF(veget(ji, jv) > veget_max(ji, jv) * ( un - &
                  exp( - lai_irrig_min * ext_coeff_vegetfrac(jv) ) ) ) THEN

                    lai_irrig_trig(ji) = lai_irrig_trig(ji) + 1

                  ENDIF

                ENDDO
              ENDIF

            ENDDO
            !If none of the PFT inside the soil tile have LAI >  lai_irrig_min (I.E. lai_irrig_trig(ji) = 0 )
            !The root deficit is set to zero, and irrigation is not triggered
            DO ji = 1,kjpindex

              IF( lai_irrig_trig(ji) < 1 ) THEN
                root_deficit(ji) = zero
              ENDIF

            ENDDO
       ENDIF

       ! Soil wetness profiles (W-Ww)/(Ws-Ww)
       ! soil_wet_ns is the ratio of available soil moisture to max available soil moisture
       ! (ie soil moisture at saturation minus soil moisture at wilting point).
       ! soil wet is a water stress for stomate, to control C decomposition
       ! Based on liquid water content
       ! 
       ! 202603:
       ! In mict peat:
       !soil_wet_ns(ji,jsl,jst) = MIN(un, MAX(zero, &
       !           (sm(ji,jsl)-smw_tmp(ji,jsl))*(sms(ji,jsl)-smw(ji,jsl))/(sms_tmp(ji,jsl)-smw_tmp(ji,jsl))**2 ))
       DO jsl=1,nslm
          DO ji=1,kjpindex
             IF (soiltile(ji,jst) .GT. min_sechiba) THEN
                soil_wet_ns(ji,jsl,jst) = MIN(un, MAX(zero, &
                     (sm(ji,jsl,jst)-smw(ji,jsl))/(sms(ji,jsl)-smw(ji,jsl)) ))
             ELSE
                soil_wet_ns(ji,jsl,jst) = zero
             ENDIF
          END DO
       END DO
      
       DO ji=1,kjpindex
          humtot_soil_liquid(ji,jst) = SUM(sm(ji,:,jst))
       END DO

       ! Compute us and the new humrelv to use in sechiba (with loops on the vegetation types)
       ! This is the water stress for transpiration (diffuco) and photosynthesis (diffuco)
       ! humrel is never used in stomate
       ! Based on liquid water content

       ! -- PFT1
       humrelv(:,1,jst) = zero        
       ! -- Top layer
       DO jv = 2,nvm
          DO ji=1,kjpindex
             !- Here we make the assumption that roots do not take water from the 1st layer. 
             us(ji,jv,jst,1) = zero
             humrelv(ji,jv,jst) = zero ! initialisation of the sum
          END DO
       ENDDO

    ENDDO ! Interrupt MAIN LOOP over soil tiles, corresponds to the first DO loop over nstm
     
    ! There are two different ways of looking at a root profile in the code. It 
    ! could reflect "structure" or "function". The code uses a different 
    ! root profile depending on what it is used for. A structural and functional 
    ! root profile are calculated below.
    CALL hydrol_root_profile(kjpindex, altmax, sm, smw, root_profile, root_depth)
       
    ! Make root_dens XIOS proof. Use NAN instead of zero to obtain the correct mean
    ! value for the period that roots are present.
    tmp(:,:,:,:) = root_profile(:,:,:,:)
    DO jv = 1,nvm
       DO jsl=1,nslm
          WHERE (SUM(circ_class_biomass(:,jv,:,iroot,icarbon),dim=2).LT.min_stomate)
             tmp(:,jv,jsl,istruc) = xios_default_val
             tmp(:,jv,jsl,ifunc) = xios_default_val
          END WHERE
       END DO
    END DO
    CALL xios_orchidee_send_field("ROOT_PROF_STRUC",tmp(:,:,:,istruc))
    CALL xios_orchidee_send_field("ROOT_PROF_FUNC",tmp(:,:,:,ifunc))

    !
    !Restart MAIN LOOP over soil tiles, the second loop over nstm, 
    !in order to to compute us, humrelv etc.
    !
    DO jst = 1,nstm 
       ! Intermediate and bottom layers
       DO jsl = 2,nslm
          DO jv = 2, nvm
             DO ji=1,kjpindex
                IF (soiltile(ji, jst) .GT. min_sechiba) THEN
                   ! AD16*** Although plants can only withdraw liquid water, we compute here the water stress
                   ! based on mc and the corresponding thresholds mcs, pcent, or potentially mcw and mcfc
                   ! This is consistent with assuming that ice is uniformly distributed within the poral space
                   ! In such a case, freezing makes mcl and the "liquid" porosity smaller than the "total" values
                   ! And it is the same for all the moisture thresholds, which are proportional to porosity.
                   ! Since the stress is based on relative moisture, it could thus independent from the porosity
                   ! at first order, thus independent from freezing.   
                   ! 26-07-2017: us and humrel now based on liquid soil moisture, so the stress is stronger
                   IF(new_watstress) THEN
                      IF((sm(ji,jsl,jst)-smw(ji,jsl)) .GT. min_sechiba) THEN
                         us(ji,jv,jst,jsl) = MIN(un, MAX(zero, &
                              (EXP(- alpha_watstress_pft(jv) * &
                              ( (smf(ji,jsl) - smw(ji,jsl)) / ( sm_nostress(ji,jsl) - smw(ji,jsl)) ) * &
                              ( (sm_nostress(ji,jsl) - sm(ji,jsl,jst)) / ( sm(ji,jsl,jst) - smw(ji,jsl)) ) ) ) ))&
                              * root_profile(ji,jv,jsl,ifunc)
                      ELSE
                         us(ji,jv,jst,jsl) = 0.
                      ENDIF
                   ELSE
                      ! orig trunk
                      us(ji,jv,jst,jsl) = MIN(un, MAX(zero, &
                           (sm(ji,jsl,jst)-smw(ji,jsl))/(sm_nostress(ji,jsl)-smw(ji,jsl)) )) * root_profile(ji,jv,jsl,ifunc)
                   ENDIF
                ELSE
                   us(ji,jv,jst,jsl) = 0.
                ENDIF
                humrelv(ji,jv,jst) = humrelv(ji,jv,jst) + us(ji,jv,jst,jsl)
             END DO
          END DO
       ENDDO

       !! vegstressv is the water stress for phenology in stomate
       !! It varies linearly from zero at wilting point to 1 at field capacity
       vegstressv(:,:,jst) = zero
       DO jv = 2, nvm
          DO ji=1,kjpindex
             DO jsl=1,nslm
                IF (soiltile(ji, jst) .GT. min_sechiba) THEN
                   vegstressv(ji,jv,jst) = vegstressv(ji,jv,jst) + &
                        MIN(un, MAX(zero, (sm(ji,jsl,jst)-smw(ji,jsl))/(smf(ji,jsl)-smw(ji,jsl)) )) &
                        * root_profile(ji,jv,jsl,ifunc)
                ELSE
                   vegstressv(ji,jv,jst) = 0.
                ENDIF
             END DO
          END DO
       END DO


       ! -- If the PFT is absent, the corresponding humrelv and vegstressv = 0
       DO jv = 2, nvm
          DO ji = 1, kjpindex
             IF (vegetmax_soil(ji,jv,jst) .LT. min_sechiba) THEN
                humrelv(ji,jv,jst) = zero
                vegstressv(ji,jv,jst) = zero
                us(ji,jv,jst,:) = zero
             ENDIF
          END DO
       END DO

       !! 6.2 We need to turn off evaporation when is_under_mcr
       !!     We set us, humrelv and vegstressv to zero in this case
       !!     WARNING: It's different from having locally us=0 in the soil layers(s) where mc<mcr
       !!              This part is crucial to preserve water conservation
       DO jsl = 1,nslm
          DO jv = 2, nvm
             WHERE (is_under_mcr(:,jst))
                us(:,jv,jst,jsl) = zero
             ENDWHERE
          ENDDO
       ENDDO
       DO jv = 2, nvm
          WHERE (is_under_mcr(:,jst))
             humrelv(:,jv,jst) = zero
          ENDWHERE
       ENDDO
       !rwilt and soil_wet_ns to zero in this case. 
       ! They are used later for shumdiag and shumdiag_perma
       DO jsl = 1,nslm
          WHERE (is_under_mcr(:,jst))
             soil_wet_ns(:,jsl,jst) = zero
          ENDWHERE
       ENDDO

       ! Counting the nb of under_mcr occurences in each grid-cell 
       WHERE (is_under_mcr(:,jst))
          undermcr = undermcr + un
       ENDWHERE

       !! 6.3 Calculate the volumetric soil moisture content (mc_layh and mcl_layh) needed in 
       !!     thermosoil for the thermal conductivity. 
       !! The multiplication by vegtot creates grid-cell average values
       ! *** To be checked for consistency with the use of nobio properties in thermosoil
       mc_layh_s = mc
       mcl_layh_s = mc
       DO jsl=1,nslm
          DO ji=1,kjpindex
             mc_layh(ji,jsl) = mc_layh(ji,jsl) + mc(ji,jsl,jst) * soiltile(ji,jst) * vegtot(ji) 
             mcl_layh(ji,jsl) = mcl_layh(ji,jsl) + mcl(ji,jsl,jst) * soiltile(ji,jst) * vegtot(ji)
          ENDDO
       END DO

       !! 6.4 The hydraulic conductivities exported here are the ones used in the diffusion/redistribution
       ! (no call of hydrol_soil_coef since 2.1)
       ! We average the values of each soiltile and keep the specific value (no multiplication by vegtot)
       !
       ! This step hydrol_soil_coef is added here, not to update a, b, c, k values that were computed in the previous iteration of nstm
       ! but only to get a, b, c ,k values corresponding to the right soil tile used below.
       ! thus, not need to call hydrol_setup here.
       CALL hydrol_soil_coef(mcr, mcs, kjpindex,jst,soiltile)
       DO ji = 1, kjpindex
          ! k: ji x jsl (grids x vertical levels)
          kk_moy(ji,:) = kk_moy(ji,:) + soiltile(ji,jst) * k(ji,:) 
          kk(ji,:,jst) = k(ji,:)
       ENDDO
       
       
       !! 6.5 We also want to export ksat at each node for CMIP6
       !  (In the output, done only once according to field_def_orchidee.xml; same averaging as for kk)
       ! No need to modify it for peat, because this one is only for the purpose of outputing the variable
       DO jsl = 1, nslm
          ksat(:,jsl) = ksat(:,jsl) + soiltile(:,jst) * &
               ( ks(:) * kfact(jsl,:) * kfact_root(:,jsl,jst) ) 
       ENDDO

      IF (printlev>=3) WRITE (numout,*) ' prognostic/diagnostic part of hydrol_soil done for jst =', jst          

    END DO  ! end of loop on soiltile 

    !! -- ENDING THE MAIN LOOP ON SOILTILES

    ! here we compute tmc again with the condition of ok_peat_hydro
    ! need to update tmc for peatland by taking into account run2peat or run2man
    tmc_soil(:,:) = zero
    IF (ok_peat_hydro) THEN
       DO ji=1,kjpindex
          IF (agri_peat .AND. ok_run2peat) THEN
             run2peat(ji)= run2peat(ji) + ru_ns(ji,5)*soiltile(ji,5)
             ru_ns(ji,5) = zero
          ENDIF
          DO jst = 1, nstm
             IF (jst .EQ. 4) THEN
                IF (ok_run2peat .AND. (soiltile(ji,4) .GT. min_sechiba)) THEN
                   water2infilt(ji,jst)=water2infilt(ji,jst)+ run2peat(ji)/soiltile(ji,4)
                   ! no need to modify the code here if peat pft is 0:
                   ! in that case, we do not modify the value of water2infilter simply
                ENDIF
                water2infilt(ji,jst)=water2infilt(ji,jst)+ wt_ab(ji)
             ENDIF
             
             IF (tides) THEN
                IF (jst .EQ. 6) THEN
                   IF (ok_run2peat .AND. (soiltile(ji,6) .GT. min_sechiba)) THEN
                      water2infilt(ji,jst)=water2infilt(ji,jst)+ run2man(ji)/soiltile(ji,6)
                   ENDIF
                   water2infilt(ji,jst)=water2infilt(ji,jst)+ wt_ab(ji)
                ENDIF
             ENDIF
          ENDDO
          DO jst = 1, nstm
             ! when peatland exits, the water2infilt for differen soiltile would be different than usual
             ! tmc_soil: the classical tmc 
             ! tmc: the modified tmc by considering water2infilt (with runoff and above-ground watertable included)
             ! tides not considered here.
             IF ((jst .GE. 4) .AND. (jst .LE. 5) .AND. (soiltile(ji,jst) .GT. min_sechiba)) THEN
                tmc_soil(ji,jst)= tmc(ji,jst)
                tmc(ji,jst) = tmc(ji,jst) + water2infilt(ji,jst)
             ENDIF
          ENDDO
       ENDDO
    ENDIF

    !! 7. Summing 3d variables into 2d variables 
    CALL hydrol_diag_soil (ks, nvan, avan, mcr, mcs, mcfc, mcw, kjpindex, veget_max, soiltile, runoff, drainage, &
         & evapot, vevapnu, returnflow, reinfiltration, irrigation, &
         & shumdiag,shumdiag_perma, k_litt, litterhumdiag, humrel, &
         & vegstress, drysoil_frac, tot_melt, us, &
         & precip_rain, totfrac_nobio, frac_snow_nobio, &
         & shumdiag_peat, shumdiag_croppeat, shumdiag_man)

    ! Means of wtd, runoff and drainage corrections, across soiltiles    
    wtd(:) = zero 
    ru_corr(:) = zero
    ru_corr2(:) = zero
    dr_corr(:) = zero
    dr_corrnum(:) = zero
    dr_force(:) = zero
    DO jst = 1, nstm
       DO ji = 1, kjpindex 
          wtd(ji) = wtd(ji) + soiltile(ji,jst) * wtd_ns(ji,jst) ! average over vegtot only
          IF (vegtot(ji) .GT. min_sechiba) THEN ! to mimic hydrol_diag_soil
             ! We average the values of each soiltile and multiply by vegtot to transform to a grid-cell mean
             ru_corr(ji) = ru_corr(ji) + vegtot(ji) * soiltile(ji,jst) * ru_corr_ns(ji,jst) 
             ru_corr2(ji) = ru_corr2(ji) + vegtot(ji) * soiltile(ji,jst) * ru_corr2_ns(ji,jst) 
             dr_corr(ji) = dr_corr(ji) + vegtot(ji) * soiltile(ji,jst) * dr_corr_ns(ji,jst) 
             dr_corrnum(ji) = dr_corrnum(ji) + vegtot(ji) * soiltile(ji,jst) * dr_corrnum_ns(ji,jst)
             dr_force(ji) = dr_force(ji) - vegtot(ji) * soiltile(ji,jst) * dr_force_ns(ji,jst)
                                       ! the sign is OK to get a negative drainage flux
          ENDIF
       ENDDO
    ENDDO

    ! Means local variables, including water conservation checks
    ru_infilt(:)=0.
    qinfilt(:)=0.
    check_infilt(:)=0.
    check_tr(:)=0.
    check_over(:)=0.
    check_under(:)=0.
    qflux(:,:)=0.
    check_top(:)=0.
    DO jst = 1, nstm
       DO ji = 1, kjpindex 
          IF (vegtot(ji) .GT. min_sechiba) THEN ! to mimic hydrol_diag_soil
             ! We average the values of each soiltile and multiply by vegtot to transform to a grid-cell mean
             ru_infilt(ji) = ru_infilt(ji) + vegtot(ji) * soiltile(ji,jst) * ru_infilt_ns(ji,jst)
             qinfilt(ji) = qinfilt(ji) + vegtot(ji) * soiltile(ji,jst) * qinfilt_ns(ji,jst)
          ENDIF
       ENDDO
    ENDDO
 
    IF (check_cwrr) THEN
       DO jst = 1, nstm
          DO ji = 1, kjpindex 
             IF (vegtot(ji) .GT. min_sechiba) THEN ! to mimic hydrol_diag_soil
                ! We average the values of each soiltile and multiply by vegtot to transform to a grid-cell mean
                check_infilt(ji) = check_infilt(ji) + vegtot(ji) * soiltile(ji,jst) * check_infilt_ns(ji,jst)
                check_tr(ji) = check_tr(ji) + vegtot(ji) * soiltile(ji,jst) * check_tr_ns(ji,jst)
                check_over(ji) = check_over(ji) + vegtot(ji) * soiltile(ji,jst) * check_over_ns(ji,jst)
                check_under(ji) =  check_under(ji) + vegtot(ji) * soiltile(ji,jst) * check_under_ns(ji,jst)
                !
                qflux(ji,:) = qflux(ji,:) + vegtot(ji) * soiltile(ji,jst) * qflux_ns(ji,:,jst)
                check_top(ji) =  check_top(ji) + vegtot(ji) * soiltile(ji,jst) * check_top_ns(ji,jst)
             ENDIF
          ENDDO
       ENDDO
    END IF

    !! 8. COMPUTING EVAP_BARE_LIM_NS FOR NEXT TIME STEP, WITH A LOOP ON SOILTILES
    !!    The principle is to run a dummy integration of the water redistribution scheme
    !!    to check if the SM profile can sustain a potential evaporation.
    !!    If not, the dummy integration is redone from the SM profile of the end of the normal integration,
    !!    with a boundary condition leading to a very severe water limitation: mc(1)=mcr

    ! evap_bare_lim = beta factor for bare soil evaporation
    evap_bare_lim(:) = zero
    evap_bare_lim_ns(:,:) = zero

    ! Loop on soil tiles  
    DO jst = 1,nstm

       !! 8.1 Save actual mc, mcl, and tmc for restoring at the end of the time step
       !!      and calculate tmcint corresponding to mc without water2infilt
       DO jsl = 1, nslm
          DO ji = 1, kjpindex
             mcint(ji,jsl) = mask_soiltile(ji,jst) * mc(ji,jsl,jst)
             mclint(ji,jsl) = mask_soiltile(ji,jst) * mcl(ji,jsl,jst)
          ENDDO
       ENDDO

       DO ji = 1, kjpindex
          temp(ji) = tmc(ji,jst)
          tmcint(ji) = temp(ji) - water2infilt(ji,jst) ! to estimate bare soil evap based on water budget
       ENDDO

       !! 8.2 Since we estimate bare soile evap for the next time step, we update profil_froz_hydro and mcl
       !     (effect of mc only, the change in stempdiag is neglected)
       IF ( ok_freeze_cwrr ) CALL hydrol_soil_froz(nvan, avan, mcr, mcs,kjpindex,jst,stempdiag, soiltile)

        DO jsl = 1, nslm
          DO ji =1, kjpindex
             IF ( (ok_peat_hydro) .AND. (jst .GE. 4) .AND. (jst .LE. 6)  ) THEN
                IF ((soiltile(ji,jst) .GT. min_sechiba)) THEN
                   mcl(ji,jsl,jst)= MIN( mc(ji,jsl,jst), mcr_peat + &
                        (un-profil_froz_hydro_ns(ji,jsl,jst))*(mc(ji,jsl,jst)-mcr_peat) )
                ELSE
                    mcl(ji,jsl,jst) = 0.0
                ENDIF
             ELSE
                mcl(ji,jsl,jst)= MIN( mc(ji,jsl,jst), mcr(ji) + &
                     (un-profil_froz_hydro_ns(ji,jsl,jst))*(mc(ji,jsl,jst)-mcr(ji)) )
                ! if profil_froz_hydro_ns=0 (including NOT ok_freeze_cwrr) we keep mcl=mc
             ENDIF
          ENDDO
       ENDDO          

       !! 8.3 K and D are recomputed for the updated profile of mc/mcl
       CALL hydrol_soil_coef(mcr, mcs, kjpindex,jst, soiltile)

       !! for the hydraulic architecture we need to pass the hydraulic
       !  conductivity. We save this variable in ksoil 
       ksoil(:,:,jst) = k

       !! 8.4 Set the tridiagonal matrix coefficients for the diffusion/redistribution scheme
       CALL hydrol_soil_setup(kjpindex,jst)
       resolv(:) = (mask_soiltile(:,jst) .GT. 0) 

       !! 8.5 We define the system of linear equations, based on matrix coefficients, 

       !- Impose potential evaporation as flux_top in mm/step, assuming the water is available
       ! Note that this should lead to never have evapnu>evapot_penm(ji)

       DO ji = 1, kjpindex
          IF (vegtot(ji).GT.min_sechiba) THEN
             
             ! We calculate a reduced demand, by means of a soil resistance (Sellers et al., 1992)
             ! It is based on the liquid SM only, like for us and humrel
             IF (do_rsoil) THEN
                mc_rel(ji) = tmc_litter(ji,jst)/tmcs_litter(ji) ! tmc_litter based on mcl
                ! based on SM in the top 4 soil layers (litter) to smooth variability
                r_soil_ns(ji,jst) = rsoil_scale_tile(ji,jst)*exp(8.206 - 4.255 * mc_rel(ji))
             ELSE
                r_soil_ns(ji,jst) = zero
             ENDIF

             ! Aerodynamic resistance
             speed = MAX(min_wind, SQRT (u(ji)*u(ji) + v(ji)*v(ji)))
             IF (speed * tq_cdrag(ji) .GT. min_sechiba) THEN
                ra = un / (speed * tq_cdrag(ji))
                evap_soil(ji) = evapot_penm(ji) / (un + r_soil_ns(ji,jst)/ra)
             ELSE
                evap_soil(ji) = evapot_penm(ji)
             ENDIF
                          
             flux_top(ji) = evap_soil(ji) * &
                  AINT(frac_bare_ns(ji,jst)+un-min_sechiba)
          ELSE

             flux_top(ji) = zero
             ! r_soil_ns needs a value to support the calculation in 
             ! section "evap_bar_lim is the grid-cell scale beta"
             r_soil_ns(ji,jst) = zero
             
          ENDIF
       ENDDO

       ! We don't use rootsinks, because we assume there is no transpiration in the bare soil fraction (??)
       !- First layer
       DO ji = 1, kjpindex
          tmat(ji,1,1) = zero
          tmat(ji,1,2) = f(ji,1)
          tmat(ji,1,3) = g1(ji,1)
          rhs(ji,1)    = fp(ji,1) * mcl(ji,1,jst) + gp(ji,1)*mcl(ji,2,jst) &
               - flux_top(ji) - (b(ji,1)+b(ji,2))/deux *(dt_sechiba/one_day) - rootsink(ji,1,jst)
       ENDDO
       !- soil body
       DO jsl=2, nslm-1
          DO ji = 1, kjpindex
             tmat(ji,jsl,1) = e(ji,jsl)
             tmat(ji,jsl,2) = f(ji,jsl)
             tmat(ji,jsl,3) = g1(ji,jsl)
             rhs(ji,jsl) = ep(ji,jsl)*mcl(ji,jsl-1,jst) + fp(ji,jsl)*mcl(ji,jsl,jst) &
                  +  gp(ji,jsl) * mcl(ji,jsl+1,jst) &
                  + (b(ji,jsl-1) - b(ji,jsl+1)) * (dt_sechiba/one_day) / deux &
                  - rootsink(ji,jsl,jst)
          ENDDO
       ENDDO
       !- Last layer
       DO ji = 1, kjpindex
          jsl=nslm
          tmat(ji,jsl,1) = e(ji,jsl)
          tmat(ji,jsl,2) = f(ji,jsl)
          tmat(ji,jsl,3) = zero
          rhs(ji,jsl) = ep(ji,jsl)*mcl(ji,jsl-1,jst) + fp(ji,jsl)*mcl(ji,jsl,jst) &
               + (b(ji,jsl-1) + b(ji,jsl)*(un-deux*free_drain_coef(ji,jst))) * (dt_sechiba/one_day) / deux &
               - rootsink(ji,jsl,jst)
       ENDDO
       !- Store the equations for later use (9.6)
       DO jsl=1,nslm
          DO ji = 1, kjpindex
             srhs(ji,jsl) = rhs(ji,jsl)
             stmat(ji,jsl,1) = tmat(ji,jsl,1)
             stmat(ji,jsl,2) = tmat(ji,jsl,2)
             stmat(ji,jsl,3) = tmat(ji,jsl,3)
          ENDDO
       ENDDO

       !! 8.6 Solve the diffusion equation, assuming that flux_top=evapot_penm (updates mcl)
       CALL hydrol_soil_tridiag(kjpindex,jst)

       !! 8.7 Alternative solution with mc(1)=mcr in points where the above solution leads to mcl<mcr 
       ! hydrol_soil_tridiag calculates mc recursively from the top as a fonction of rhs and tmat 
       ! We re-use these the above values, but for mc(1)=mcr and the related tmat
       
       DO ji = 1, kjpindex
          ! by construction, mc and mcl are always on the same side of mcr, so we can use mcl here
          IF ( (ok_peat_hydro) .AND. (jst .GE. 4) .AND. (jst .LE. 6)) THEN
             IF ((soiltile(ji,jst) .GT. min_sechiba) ) THEN
                resolv(ji) = (mcl(ji,1,jst).LT.(mcr_peat).AND.flux_top(ji).GT.min_sechiba)
             ELSE
                resolv(ji) = (mcl(ji,1,jst).LT.(mcr(ji)).AND.flux_top(ji).GT.min_sechiba)
             ENDIF
                    
          ELSE
              resolv(ji) = (mcl(ji,1,jst).LT.(mcr(ji)).AND.flux_top(ji).GT.min_sechiba)
           ENDIF
       ENDDO
       !! Reset the coefficient for diffusion (tridiag is only solved if resolv(ji) = .TRUE.)O
       DO jsl=1,nslm
          !- The new condition is to put the upper layer at residual soil moisture
          DO ji = 1, kjpindex
             rhs(ji,jsl) = srhs(ji,jsl)
             tmat(ji,jsl,1) = stmat(ji,jsl,1)
             tmat(ji,jsl,2) = stmat(ji,jsl,2)
             tmat(ji,jsl,3) = stmat(ji,jsl,3)
          END DO
       END DO
       
       DO ji = 1, kjpindex
          tmat(ji,1,2) = un
          tmat(ji,1,3) = zero
          
          rhs(ji,1) = mcr(ji)
          IF ( (ok_peat_hydro) .AND. (jst .GE. 4) .AND. (jst .LE. 6) ) THEN
             IF ((soiltile(ji,jst) .GT. min_sechiba)) THEN
                rhs(ji,1) = mcr_peat
             ENDIF
          ENDIF
       ENDDO
       
       ! Solves the diffusion equation with new surface bc where resolv=T 
       CALL hydrol_soil_tridiag(kjpindex,jst)

       !! 8.8 In both case, we have drainage to be consistent with rhs
       DO ji = 1, kjpindex
          flux_bottom(ji) = mask_soiltile(ji,jst)*k(ji,nslm)*free_drain_coef(ji,jst) * (dt_sechiba/one_day)
       ENDDO
       
       !! 8.9 Water budget to assess the top flux = soil evaporation
       !      Where resolv=F at the 2nd step (9.6), it should simply be the potential evaporation

       ! Total soil moisture content for water budget

       DO jsl = 1, nslm
          DO ji =1, kjpindex
              IF ( (ok_peat_hydro) .AND. (jst .GE. 4) .AND. (jst .LE. 6) ) THEN
                IF (soiltile(ji,jst) .GT. min_sechiba ) THEN
                   mc(ji,jsl,jst) = MAX( mcl(ji,jsl,jst), mcl(ji,jsl,jst) + &
                        profil_froz_hydro_ns(ji,jsl,jst)*(mc(ji,jsl,jst)-mcr_peat))
                ELSE
                   mc(ji,jsl,jst) = 0.0 
                ENDIF
             ELSE 
                mc(ji,jsl,jst) = MAX( mcl(ji,jsl,jst), mcl(ji,jsl,jst) + &
                     profil_froz_hydro_ns(ji,jsl,jst)*(mc(ji,jsl,jst)-mcr(ji)) )
                ! if profil_froz_hydro_ns=0 (including NOT ok_freeze_cwrr) we get mc=mcl
             ENDIF
          ENDDO
       ENDDO
       
       DO ji = 1, kjpindex
          tmc(ji,jst) = dz(2) * ( trois*mc(ji,1,jst) + mc(ji,2,jst) )/huit
       ENDDO       
       DO jsl = 2,nslm-1
          DO ji = 1, kjpindex
             tmc(ji,jst) = tmc(ji,jst) + dz(jsl) &
                  * (trois*mc(ji,jsl,jst)+mc(ji,jsl-1,jst))/huit &
                  + dz(jsl+1) * (trois*mc(ji,jsl,jst)+mc(ji,jsl+1,jst))/huit
          ENDDO
       ENDDO
       DO ji = 1, kjpindex
          tmc(ji,jst) = tmc(ji,jst) + dz(nslm) &
               * (trois * mc(ji,nslm,jst) + mc(ji,nslm-1,jst))/huit
       END DO
    
       ! Deduce upper flux from soil moisture variation and bottom flux
       ! TMCi-D-BSE=TMC (BSE=bare soil evap=TMCi-TMC-D)
       ! The numerical errors of tridiag close to saturation cannot be simply solved here,
       ! we can only hope they are not too large because we don't add water at this stage... 
       DO ji = 1, kjpindex
          evap_bare_lim_ns(ji,jst) = mask_soiltile(ji,jst) * &
               (tmcint(ji)-tmc(ji,jst)-flux_bottom(ji)-SUM(rootsink(ji,:,jst)))
       END DO

       !! 8.10 evap_bare_lim_ns is turned from an evaporation rate to a beta
       DO ji = 1, kjpindex
          ! Here we weight evap_bare_lim_ns by the fraction of bare evaporating soil. 
          ! This is given by frac_bare_ns, taking into account bare soil under vegetation
          IF(vegtot(ji) .GT. min_sechiba) THEN
             evap_bare_lim_ns(ji,jst) = evap_bare_lim_ns(ji,jst) * frac_bare_ns(ji,jst)
          ELSE
             evap_bare_lim_ns(ji,jst) = 0.
          ENDIF
       END DO

       ! We divide by evapot, which is consistent with diffuco (evap_bare_lim_ns < evapot_penm/evapot)
       ! Further decrease if tmc_litter is below the wilting point

       IF (do_rsoil) THEN
          DO ji=1,kjpindex
             IF (evapot(ji).GT.min_sechiba) THEN
                evap_bare_lim_ns(ji,jst) = evap_bare_lim_ns(ji,jst) / evapot(ji)
             ELSE
                evap_bare_lim_ns(ji,jst) = zero ! not redundant with the is_under_mcr case below
                                                ! but not necessarily useful
             END IF
             evap_bare_lim_ns(ji,jst)=MAX(MIN(evap_bare_lim_ns(ji,jst),1.),0.)
          END DO
       ELSE
          DO ji=1,kjpindex
             IF ((evapot(ji).GT.min_sechiba) .AND. &
                  (tmc_litter(ji,jst).GT.(tmc_litter_wilt(ji,jst)))) THEN
                evap_bare_lim_ns(ji,jst) = evap_bare_lim_ns(ji,jst) / evapot(ji)
             ELSEIF((evapot(ji).GT.min_sechiba).AND. &
                  (tmc_litter(ji,jst).GT.(tmc_litter_res(ji,jst)))) THEN
                evap_bare_lim_ns(ji,jst) =  (un/deux) * evap_bare_lim_ns(ji,jst) / evapot(ji)
                ! This is very arbitrary, with no justification from the literature
             ELSE
                evap_bare_lim_ns(ji,jst) = zero
             END IF
             evap_bare_lim_ns(ji,jst)=MAX(MIN(evap_bare_lim_ns(ji,jst),1.),0.)
          END DO
       ENDIF

       !! 8.11 Set evap_bare_lim_ns to zero if is_under_mcr at the end of the prognostic loop
       !!      (cf us, humrelv, vegstressv in 5.2)
       WHERE (is_under_mcr(:,jst))
          evap_bare_lim_ns(:,jst) = zero
       ENDWHERE

       !! 8.12 Restores mc, mcl, and tmc, to erase the effect of the dummy integrations
       !!      on these prognostic variables 
       DO jsl = 1, nslm
          DO ji = 1, kjpindex
             mc(ji,jsl,jst) = mask_soiltile(ji,jst) * mcint(ji,jsl)
             mcl(ji,jsl,jst) = mask_soiltile(ji,jst) * mclint(ji,jsl)
          ENDDO
       ENDDO
       DO ji = 1, kjpindex
          tmc(ji,jst) = temp(ji)
       ENDDO

    ENDDO !end loop on tiles for dummy integration

    ! peatland
    ! modified the following code (of mict-peat) to be consistent with the trunk
    IF (ok_peat_hydro) THEN       
       tmcs_peat=dz(2) * ( trois*mcs_peat + mcs_peat )/huit
       DO jsl = 2,nslm-1
          tmcs_peat=tmcs_peat+dz(jsl)*(trois*mcs_peat+mcs_peat)/huit&
               + dz(jsl+1) * (trois*mcs_peat+mcs_peat)/huit
       ENDDO
       tmcs_peat = tmcs_peat + dz(nslm) &
               * (trois * mcs_peat + mcs_peat)/huit

       ! compute tmcl_peat: 
       ! mcl of peatland (4) is used directly in the computation below.
       ! We suppose that when ok_peat_hydro is activated,soiltile 4 is peat boreal
       ! soiltile 6 is tropical,soiltile 5 is agriPeat
       ! only soiltile 4 is used in the computation below.
       DO ji = 1, kjpindex
          tmcl_peat(ji) = dz(2) * ( trois*mcl(ji,1,4) + mcl(ji,2,4) )/huit
       ENDDO
       DO jsl = 2,nslm-1
          DO ji = 1, kjpindex
             tmcl_peat(ji) = tmcl_peat(ji) + dz(jsl) &
                  * (trois*mcl(ji,jsl,4)+mcl(ji,jsl-1,4))/huit &
                  + dz(jsl+1) * (trois*mcl(ji,jsl,4)+mcl(ji,jsl+1,4))/huit
          ENDDO
       ENDDO
     
       DO ji = 1, kjpindex
          tmcl_peat(ji) = tmcl_peat(ji) + dz(nslm) &
               * (trois * mcl(ji,nslm,4) + mcl(ji,nslm-1,4))/huit
       END DO
       ! liqwt_ratio is useful in stomate_main. liquid water ratio in peat.
       ! used in sechiba and slow proc
       DO ji=1,kjpindex
          IF (soiltile(ji,4) .GT. min_sechiba) THEN
             liqwt_ratio(ji)=tmcl_peat(ji)/tmcs_peat
          ENDIF
       ENDDO
    ENDIF

    ! peatland  peatland water table: water table of soiltile4 only
    ! We do not consider agri_peat or tides here
    IF (ok_peat_hydro) THEN
       DO ji=1,kjpindex  
          IF (soiltile (ji,4) .GT. min_sechiba) THEN    
             wtp(ji)= zero
             DO jsl=1,nslm
                h_eau(ji,jsl)=MIN(MAX(zero,(mc(ji,jsl,4)-mcr_peat)/(mcs_peat-mcr_peat)),un)
                wtp(ji)=wtp(ji) + h_eau(ji,jsl)*dz(jsl)
             ENDDO
             IF (wt_ab(ji) .GT. zero) THEN
                wtp(ji)=wtp(ji)+wt_ab(ji)
             ENDIF
             wtp(ji)= deux*mille - wtp(ji)     !!in mm  
          ELSE
             wtp(ji)= zero
          ENDIF
       ENDDO
    ENDIF

    !!! tides, mangrove water table: water table of soiltile6
    IF (tides) THEN
       DO ji=1,kjpindex  
          IF (soiltile (ji,6) .GT. min_sechiba) THEN    
             wtp(ji)= zero
             DO jsl=1,nslm
                h_eau(ji,jsl)=MIN(MAX(zero,(mc(ji,jsl,6)-mcr_peat)/(mcs_peat-mcr_peat)),un)
                wtp(ji)=wtp(ji) + h_eau(ji,jsl)*dz(jsl)
             ENDDO
             IF (wt_ab(ji) .GT. zero) THEN
                wtp(ji)=wtp(ji)+wt_ab(ji)
             ENDIF
             wtp(ji)= deux*mille - wtp(ji)     !!in mm  
          ELSE
             wtp(ji)= zero
          ENDIF
       ENDDO
    ENDIF

    !! 9. evap_bar_lim is the grid-cell scale beta
    DO ji = 1, kjpindex
       evap_bare_lim(ji) =  SUM(evap_bare_lim_ns(ji,:)*vegtot(ji)*soiltile(ji,:))
       r_soil(ji) =  SUM(r_soil_ns(ji,:)*vegtot(ji)*soiltile(ji,:))
    ENDDO
    ! si vegtot LE min_sechiba, evap_bare_lim_ns et evap_bare_lim valent zero


    !! 10. XIOS export of local variables, including water conservation checks
    CALL xios_orchidee_send_field("ksat",ksat) ! mm/d (for CMIP6, once)
    CALL xios_orchidee_send_field("psi_moy",psi_moy) ! mm (for SP-MIP)
    CALL xios_orchidee_send_field("wtd",wtd) ! in m
    CALL xios_orchidee_send_field("ru_corr",ru_corr/dt_sechiba)   ! adjustment flux added to surface runoff (included in runoff)
    CALL xios_orchidee_send_field("ru_corr2",ru_corr2/dt_sechiba)
    CALL xios_orchidee_send_field("dr_corr",dr_corr/dt_sechiba)   ! adjustment flux added to drainage (included in drainage)
    CALL xios_orchidee_send_field("dr_corrnum",dr_corrnum/dt_sechiba) 
    CALL xios_orchidee_send_field("dr_force",dr_force/dt_sechiba) ! adjustement flux added to drainage to sustain a forced wtd
    CALL xios_orchidee_send_field("qinfilt",qinfilt/dt_sechiba)
    CALL xios_orchidee_send_field("ru_infilt",ru_infilt/dt_sechiba)
    CALL xios_orchidee_send_field("r_soil",r_soil) ! s/m
    CALL xios_orchidee_send_field("humtot_soiliq",humtot_soil_liquid) 

    IF (check_cwrr) THEN
       CALL xios_orchidee_send_field("check_infilt",check_infilt/dt_sechiba)
       CALL xios_orchidee_send_field("check_tr",check_tr/dt_sechiba)
       CALL xios_orchidee_send_field("check_over",check_over/dt_sechiba)
       CALL xios_orchidee_send_field("check_under",check_under/dt_sechiba) 
       ! Variables calculated in hydrol_diag_soil_flux
       CALL xios_orchidee_send_field("qflux",qflux/dt_sechiba) ! upward water flux at the low interface of each layer
       CALL xios_orchidee_send_field("check_top",check_top/dt_sechiba) !water budget residu in top layer 
    END IF

    IF (ok_peat_hydro) THEN
       CALL xios_orchidee_send_field("tmc_soil",tmc_soil)
    ENDIF


  END SUBROUTINE hydrol_soil


!! ================================================================================================================================
!! SUBROUTINE   : hydrol_soil_infilt
!!
!>\BRIEF        Infiltration
!!
!! DESCRIPTION  :
!! 1. We calculate the total SM at the beginning of the routine
!! 2. Infiltration process
!! 2.1 Initialization of time counter and infiltration rate
!! 2.2 Infiltration layer by layer, accounting for an exponential law for subgrid variability
!! 2.3 Resulting infiltration and surface runoff
!! 3. For water conservation check, we calculate the total SM at the beginning of the routine,
!!    and export the difference with the flux
!! 5. Local verification 
!!
!! RECENT CHANGE(S) : 2016 by A. Ducharne
!! Adding checks and interactions variables with hydrol_soil, but the processes are unchanged
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_soil_infilt

  SUBROUTINE hydrol_soil_infilt(ks, mcr, mcs, mcfc, mcw, kjpindex, ins, flux_infilt, stempdiag, soiltile, &
                                qinfilt_ns, ru_infilt, check)

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    ! GLOBAL (in or inout)
    INTEGER(i_std), INTENT(in)                         :: kjpindex         !! Domain size
    INTEGER(i_std), INTENT(in)                         :: ins
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: ks               !! Hydraulic conductivity at saturation (mm {-1})
        
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: mcr              !! Residual volumetric water content (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: mcs              !! Saturated volumetric water content (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: mcfc             !! Volumetric water content at field capacity (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: mcw              !! Volumetric water content at wilting point (m^{3} m^{-3})
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)     :: flux_infilt      !! Water to infiltrate
                                                                           !!  @tex $(kg m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (in) :: stempdiag        !! Diagnostic temp profile from thermosoil
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT (in):: soiltile         !! Fraction of each soil tile within vegtot (0-1, unitless)


    !! 0.2 Output variables
    REAL(r_std), DIMENSION(kjpindex,nstm), INTENT(out) :: check       !! delta SM - flux (mm/dt_sechiba)
    REAL(r_std), DIMENSION(kjpindex,nstm), INTENT(out) :: ru_infilt   !! Surface runoff from soil_infilt (mm/dt_sechiba)
    REAL(r_std), DIMENSION(kjpindex,nstm), INTENT(out) :: qinfilt_ns  !! Effective infiltration flux (mm/dt_sechiba)
    

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                                :: ji, jsl      !! Indices
    REAL(r_std), DIMENSION (kjpindex)             :: wat_inf_pot  !! infiltrable water in the layer
    REAL(r_std), DIMENSION (kjpindex)             :: wat_inf      !! infiltrated water in the layer
    REAL(r_std), DIMENSION (kjpindex)             :: dt_tmp       !! time remaining before the end of the time step
    REAL(r_std), DIMENSION (kjpindex)             :: dt_inf       !! the time it takes to complete the infiltration in the 
                                                                  !! layer 
    REAL(r_std)                                   :: k_m          !! the mean conductivity used for the saturated front
    REAL(r_std), DIMENSION (kjpindex)             :: infilt_tmp   !! infiltration rate for the considered layer 
    REAL(r_std), DIMENSION (kjpindex)             :: infilt_tot   !! total infiltration 
    REAL(r_std), DIMENSION (kjpindex)             :: flux_tmp     !! rate at which precip hits the ground

    REAL(r_std), DIMENSION(kjpindex)              :: tmci         !! total SM at beginning of routine (kg/m2)
    REAL(r_std), DIMENSION(kjpindex)              :: tmcf         !! total SM at end of routine (kg/m2)
    
    REAL(r_std)                                   :: mcs_tmp      !! variable similar to mcs when peat is activated

!_ ================================================================================================================================

    ! If data (or coupling with GCM) was available, a parameterization for subgrid rainfall could be performed

    !! 1. We calculate the total SM at the beginning of the routine
    IF (check_cwrr) THEN
       tmci(:) = dz(2) * ( trois*mc(:,1,ins) + mc(:,2,ins) )/huit
       DO jsl = 2,nslm-1
          tmci(:) = tmci(:) + dz(jsl) * (trois*mc(:,jsl,ins)+mc(:,jsl-1,ins))/huit &
               + dz(jsl+1) * (trois*mc(:,jsl,ins)+mc(:,jsl+1,ins))/huit
       ENDDO
       tmci(:) = tmci(:) + dz(nslm) * (trois*mc(:,nslm,ins) + mc(:,nslm-1,ins))/huit
    ENDIF

    !! 2. Infiltration process

    !! 2.1 Initialization

    DO ji = 1, kjpindex
       !! First we fill up the first layer (about 1mm) without any resistance and quasi-immediately
       !
       ! soiltile 4: boreal natural peat, 
       ! soiltile 6: tropical peat
       ! soiltile 5: agriPeat 
       !
       ! given jst being peatland and its fraction is 0, then the computation below use the hydrol param of the whole grid,
       ! which is the same as for other soil tiles in the grid
       mcs_tmp = mcs(ji)
       IF ((ok_peat_hydro) .AND. ((ins .EQ. 4) .OR. (ins .EQ. 5))) THEN
          IF ((soiltile(ji,ins) .GT. min_sechiba)) THEN
             mcs_tmp = mcs_peat
          ENDIF
       ENDIF   
       wat_inf_pot(ji) = MAX((mcs_tmp-mc(ji,1,ins)) * dz(2) / deux, zero)
       wat_inf(ji) = MIN(wat_inf_pot(ji), flux_infilt(ji))
       mc(ji,1,ins) = mc(ji,1,ins) + wat_inf(ji) * deux / dz(2)
       !
    ENDDO

    !! Initialize a countdown for infiltration during the time-step and the value of potential runoff
    dt_tmp(:) = dt_sechiba / one_day
    infilt_tot(:) = wat_inf(:)
    !! Compute the rate at which water will try to infiltrate each layer
    ! flux_temp is converted here to the same unit as k_m
    flux_tmp(:) = (flux_infilt(:)-wat_inf(:)) / dt_tmp(:)

    !! 2.2 Infiltration layer by layer 
    DO jsl = 2, nslm-1
       DO ji = 1, kjpindex
          !! Infiltrability of each layer if under a saturated one
          ! This is computed by an simple arithmetic average because 
          ! the time step (30min) is not appropriate for a geometric average (advised by Haverkamp and Vauclin)
          !
          IF ((ok_peat_hydro) .AND. (ins .GE. 4) .AND. (ins .LE. 6)) THEN
             ! kfact_root takes into account the peatland and its fraction in the previous computation
             ! if soiltile is 0, kfact_root is 1
             IF ((soiltile(ji,ins) .GT. min_sechiba)) THEN
                k_m = (k(ji,jsl) + ks_peat*kfact_peat(jsl-1)*kfact_root(ji,jsl,ins)) / deux
             ELSE
                k_m = (k(ji,jsl) + ks(ji)*kfact(jsl-1,ji)*kfact_root(ji,jsl,ins)) / deux 
             ENDIF
          ELSE
             k_m = (k(ji,jsl) + ks(ji)*kfact(jsl-1,ji)*kfact_root(ji,jsl,ins)) / deux 
          ENDIF

          IF (ok_freeze_cwrr) THEN
             IF (stempdiag(ji, jsl) .LT. ZeroCelsius) THEN
                k_m = k(ji,jsl)
             ENDIF
          ENDIF

          !! We compute the mean rate at which water actually infiltrate:
          ! Subgrid: Exponential distribution of k around k_m, but average p directly used 
          ! See d'Orgeval 2006, p 78, but it's not fully clear to me (AD16***)
          infilt_tmp(ji) = k_m * (un - EXP(- flux_tmp(ji) / k_m)) 

          !! From which we deduce the time it takes to fill up the layer or to end the time step...
          !
          IF ((ok_peat_hydro) .AND. (ins .GE. 4) .AND. (ins .LE. 6)) THEN
             IF (soiltile(ji,ins) .GT. min_sechiba ) THEN
                wat_inf_pot(ji) =  MAX((mcs_peat-mc(ji,jsl,ins)) * (dz(jsl) + dz(jsl+1)) / deux, zero)
             ELSE
                wat_inf_pot(ji) =  MAX((mcs(ji)-mc(ji,jsl,ins)) * (dz(jsl) + dz(jsl+1)) / deux, zero)
             ENDIF
          ELSE
             wat_inf_pot(ji) =  MAX((mcs(ji)-mc(ji,jsl,ins)) * (dz(jsl) + dz(jsl+1)) / deux, zero)
          ENDIF
          !
          
          !! From which we deduce the time it takes to fill up the layer or to end the time step...
          IF ( infilt_tmp(ji) > min_sechiba) THEN
             dt_inf(ji) =  MIN(wat_inf_pot(ji)/infilt_tmp(ji), dt_tmp(ji))
             ! The water infiltration TIME has to limited by what is still available for infiltration.
             IF ( dt_inf(ji) * infilt_tmp(ji) > flux_infilt(ji)-infilt_tot(ji) ) THEN
                dt_inf(ji) = MAX(flux_infilt(ji)-infilt_tot(ji),zero)/infilt_tmp(ji)
             ENDIF
          ELSE
             dt_inf(ji) = dt_tmp(ji)
          ENDIF

          !! The water enters in the layer
          wat_inf(ji) = dt_inf(ji) * infilt_tmp(ji)
          ! bviously the moisture content
          mc(ji,jsl,ins) = mc(ji,jsl,ins) + &
               & wat_inf(ji) * deux / (dz(jsl) + dz(jsl+1))
          ! the time remaining before the next time step
          dt_tmp(ji) = dt_tmp(ji) - dt_inf(ji)
          ! and finally the infilt_tot (which is just used to check if there is a problem, below) 
          infilt_tot(ji) = infilt_tot(ji) + infilt_tmp(ji) * dt_inf(ji)
       ENDDO
    ENDDO

    !! 2.3 Resulting infiltration and surface runoff
    ru_infilt(:,ins) = flux_infilt(:) - infilt_tot(:)
    qinfilt_ns(:,ins) = infilt_tot(:)

    !! 3. For water conservation check: we calculate the total SM at the beginning of the routine
    !!    and export the difference with the flux
    IF (check_cwrr) THEN
       tmcf(:) = dz(2) * ( trois*mc(:,1,ins) + mc(:,2,ins) )/huit
       DO jsl = 2,nslm-1
          tmcf(:) = tmcf(:) + dz(jsl) * (trois*mc(:,jsl,ins)+mc(:,jsl-1,ins))/huit &
               + dz(jsl+1) * (trois*mc(:,jsl,ins)+mc(:,jsl+1,ins))/huit
       ENDDO
       tmcf(:) = tmcf(:) + dz(nslm) * (trois*mc(:,nslm,ins) + mc(:,nslm-1,ins))/huit
       ! Normally, tcmf=tmci+infilt_tot
       check(:,ins) = tmcf(:)-(tmci(:)+infilt_tot(:))
    ENDIF
    
    !! 5. Local verification
    DO ji = 1, kjpindex
       IF (infilt_tot(ji) .LT. -min_sechiba .OR. infilt_tot(ji) .GT. flux_infilt(ji) + min_sechiba) THEN
          WRITE (numout,*) 'Error in the calculation of infilt tot', infilt_tot(ji)
          WRITE (numout,*) 'k, ji, jst, mc', k(ji,1:2), ji, ins, mc(ji,1,ins)
          CALL ipslerr_p(3, 'hydrol_soil_infilt', 'We will STOP now.','Error in calculation of infilt tot','')
       ENDIF
    ENDDO

  END SUBROUTINE hydrol_soil_infilt


!! ================================================================================================================================
!! SUBROUTINE   : hydrol_soil_smooth_under_mcr
!!
!>\BRIEF        : Modifies the soil moisture profile to avoid under-residual values, 
!!                then diagnoses the points where such "excess" values remain. 
!!
!! DESCRIPTION  :
!! The "excesses" under-residual are corrected from top to bottom, by transfer of excesses 
!! to the lower layers. The reverse transfer is performed to smooth any remaining "excess" in the bottom layer.
!! If some "excess" remain afterwards, the entire soil profile is at the threshold value (mcs or mcr), 
!! and the remaining "excess" is necessarily concentrated in the top layer. 
!! This allowing diagnosing the flag is_under_mcr. 
!! Eventually, the remaining "excess" is split over the entire profile
!! 1. We calculate the total SM at the beginning of the routine
!! 2. Smoothes the profile to avoid negative values of punctual soil moisture
!! Note that we check that mc > min_sechiba in hydrol_soil 
!! 3. For water conservation check, We calculate the total SM at the beginning of the routine,
!!    and export the difference with the flux
!!
!! RECENT CHANGE(S) : 2016 by A. Ducharne
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_soil_smooth_under_mcr

  SUBROUTINE hydrol_soil_smooth_under_mcr(mcr, kjpindex, ins, soiltile, is_under_mcr, check)

    !- arguments

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                         :: kjpindex     !! Domain size
    INTEGER(i_std), INTENT(in)                         :: ins          !! Soiltile index (1-nstm, unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: mcr          !! Residual volumetric water content (m^{3} m^{-3})
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT (in)          :: soiltile         !! Fraction of each soil tile within vegtot (0-1, unitless)
    
    !! 0.2 Output variables

    LOGICAL, DIMENSION(kjpindex,nstm), INTENT(out)     :: is_under_mcr !! Flag diagnosing under residual soil moisture
    REAL(r_std), DIMENSION(kjpindex,nstm), INTENT(out) :: check        !! delta SM - flux
    

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                       :: ji,jsl
    REAL(r_std)                          :: excess
    REAL(r_std), DIMENSION(kjpindex)     :: excessji
    REAL(r_std), DIMENSION(kjpindex)     :: tmci      !! total SM at beginning of routine
    REAL(r_std), DIMENSION(kjpindex)     :: tmcf      !! total SM at end of routine

!_ ================================================================================================================================       

    !! 1. We calculate the total SM at the beginning of the routine
    IF (check_cwrr) THEN
       tmci(:) = dz(2) * ( trois*mc(:,1,ins) + mc(:,2,ins) )/huit
       DO jsl = 2,nslm-1
          tmci(:) = tmci(:) + dz(jsl) * (trois*mc(:,jsl,ins)+mc(:,jsl-1,ins))/huit &
               + dz(jsl+1) * (trois*mc(:,jsl,ins)+mc(:,jsl+1,ins))/huit
       ENDDO
       tmci(:) = tmci(:) + dz(nslm) * (trois*mc(:,nslm,ins) + mc(:,nslm-1,ins))/huit
    ENDIF

    !! 2. Smoothes the profile to avoid negative values of punctual soil moisture

    ! 2.1 smoothing from top to bottom
    DO jsl = 1,nslm-2
       DO ji=1, kjpindex
          IF ((ok_peat_hydro) .AND. (ins .GE. 4) .AND. (ins .LE. 6)) THEN
             IF ((soiltile(ji,ins) .GT. min_sechiba)) THEN
                excess = MAX(mcr_peat-mc(ji,jsl,ins),zero)
             ELSE
                excess = MAX(mcr(ji)-mc(ji,jsl,ins),zero)
             ENDIF
          ELSE
             excess = MAX(mcr(ji)-mc(ji,jsl,ins),zero)
          ENDIF
          mc(ji,jsl,ins) = mc(ji,jsl,ins) + excess
          mc(ji,jsl+1,ins) = mc(ji,jsl+1,ins) - excess * &
               &  (dz(jsl)+dz(jsl+1))/(dz(jsl+1)+dz(jsl+2))
       ENDDO
    ENDDO

    jsl = nslm-1
    DO ji=1, kjpindex
       !
       IF ((ok_peat_hydro) .AND. (ins .GE. 4) .AND. (ins .LE. 6)) THEN
          IF (soiltile(ji,ins) .GT. min_sechiba) THEN
             excess = MAX(mcr_peat-mc(ji,jsl,ins),zero)
          ELSE
             excess = MAX(mcr(ji)-mc(ji,jsl,ins),zero)
          ENDIF
       ELSE
          excess = MAX(mcr(ji)-mc(ji,jsl,ins),zero)
       ENDIF
       mc(ji,jsl,ins) = mc(ji,jsl,ins) + excess
       mc(ji,jsl+1,ins) = mc(ji,jsl+1,ins) - excess * &
            &  (dz(jsl)+dz(jsl+1))/dz(jsl+1)
    ENDDO

    ! at the bottom level, use mcr_peat or mcr for soiltile 4 and 5
    ! to be consistent with mict-peat for this case
    jsl = nslm
    DO ji=1, kjpindex
       IF ( (ok_peat_hydro) .AND. ( (ins .EQ. 4) .OR. (ins .EQ. 5)) ) THEN
          IF (soiltile(ji,ins) .GT. min_sechiba) THEN
             excess = MAX(mcr_peat-mc(ji,jsl,ins),zero)
          ELSE
             excess = MAX(mcr(ji)-mc(ji,jsl,ins),zero)
          ENDIF
       ELSE
          excess = MAX(mcr(ji)-mc(ji,jsl,ins),zero)
       ENDIF
       mc(ji,jsl,ins) = mc(ji,jsl,ins) + excess
       mc(ji,jsl-1,ins) = mc(ji,jsl-1,ins) - excess * &
               &  dz(jsl)/(dz(jsl-1)+dz(jsl))
    ENDDO

    ! 2.2 smoothing from bottom to top
    DO jsl = nslm-1,2,-1
       DO ji=1, kjpindex
          IF ((ok_peat_hydro) .AND. (ins .GE. 4) .AND. (ins .LE. 6)) THEN
             IF (soiltile(ji,ins) .GT. min_sechiba) THEN
                excess = MAX(mcr_peat-mc(ji,jsl,ins),zero)
             ELSE
                excess = MAX(mcr(ji)-mc(ji,jsl,ins),zero)
             ENDIF
          ELSE
             excess = MAX(mcr(ji)-mc(ji,jsl,ins),zero)
          ENDIF
          mc(ji,jsl,ins) = mc(ji,jsl,ins) + excess
          mc(ji,jsl-1,ins) = mc(ji,jsl-1,ins) - excess * &
               &  (dz(jsl)+dz(jsl+1))/(dz(jsl-1)+dz(jsl))
       ENDDO
    ENDDO

    ! 2.3 diagnoses is_under_mcr(ji), and updates the entire profile 
    ! excess > 0
    DO ji=1, kjpindex
       ! allowing agripeat also in this case
       IF ( (ok_peat_hydro) .AND. ( (ins .EQ. 4) .OR. (ins .EQ. 5)) ) THEN
          IF (soiltile(ji,ins) .GT. min_sechiba) THEN
             excessji(ji) = mask_soiltile(ji,ins) * MAX(mcr_peat-mc(ji,1,ins),zero)
          ELSE
             excessji(ji) = mask_soiltile(ji,ins) * MAX(mcr(ji)-mc(ji,1,ins),zero)
          ENDIF
       ELSE
          excessji(ji) = mask_soiltile(ji,ins) * MAX(mcr(ji)-mc(ji,1,ins),zero)
       ENDIF
    ENDDO
    DO ji=1, kjpindex
       mc(ji,1,ins) = mc(ji,1,ins) + excessji(ji) ! then mc(1)=mcr
       is_under_mcr(ji,ins) = (excessji(ji) .GT. min_sechiba)
    ENDDO

    ! 2.4 The amount of water corresponding to excess in the top soil layer is redistributed in all soil layers
      ! -excess(ji) * dz(2) / deux donne le deficit total, negatif, en mm
      ! diviser par la profondeur totale en mm donne des delta_mc identiques en chaque couche, en mm
      ! retransformes en delta_mm par couche selon les bonnes eqs (eqs_hydrol.pdf, Eqs 13-15), puis sommes
      ! retourne bien le deficit total en mm
    DO jsl = 1, nslm
       DO ji=1, kjpindex
          mc(ji,jsl,ins) = mc(ji,jsl,ins) - excessji(ji) * dz(2) / (deux * zmaxh*mille)
       ENDDO
    ENDDO
    ! This can lead to mc(jsl) < mcr depending on the value of excess,
    ! but this is no major pb for the diffusion
    ! Yet, we need to prevent evaporation if is_under_mcr
    
    !! Note that we check that mc > min_sechiba in hydrol_soil 

    ! We just make sure that mc remains at 0 where soiltile=0
    DO jsl = 1, nslm
       DO ji=1, kjpindex
          mc(ji,jsl,ins) = mask_soiltile(ji,ins) * mc(ji,jsl,ins)
       ENDDO
    ENDDO

    !! 3. For water conservation check, We calculate the total SM at the beginning of the routine,
    !!    and export the difference with the flux
    IF (check_cwrr) THEN
       tmcf(:) = dz(2) * ( trois*mc(:,1,ins) + mc(:,2,ins) )/huit
       DO jsl = 2,nslm-1
          tmcf(:) = tmcf(:) + dz(jsl) * (trois*mc(:,jsl,ins)+mc(:,jsl-1,ins))/huit &
               + dz(jsl+1) * (trois*mc(:,jsl,ins)+mc(:,jsl+1,ins))/huit
       ENDDO
       tmcf(:) = tmcf(:) + dz(nslm) * (trois*mc(:,nslm,ins) + mc(:,nslm-1,ins))/huit
       ! Normally, tcmf=tmci since we just redistribute the deficit 
       check(:,ins) = tmcf(:)-tmci(:)
    ENDIF
       
  END SUBROUTINE hydrol_soil_smooth_under_mcr

!! ================================================================================================================================
!! SUBROUTINE   : hydrol_soil_smooth_over_mcs2
!!
!>\BRIEF        : Modifies the soil moisture profile to avoid over-saturation values, 
!!                by putting the excess in ru_ns
!!                Thus, no point remain where such "excess" values remain 
!!
!! DESCRIPTION  :
!! The "excesses" over-saturation are corrected, by directly discarding the excess as rudr_corr,
!! to be added to ru_ns or dr_nsrunoff (via rudr_corr).
!! Therefore, there is no more smoothing, and this helps preventing the saturation of too many layers,
!! which leads to numerical errors with tridiag.
!! 1. We calculate the total SM at the beginning of the routine
!! 2. In case of over-saturation, we directly eliminate the excess via rudr_corr
!!    The calculation of the adjustement flux needs to account for nodes n-1 and n+1.
!! 3. For water conservation checks, we calculate the total SM at the beginning of the routine,
!!    and export the difference with the flux   
!!
!! RECENT CHANGE(S) : 2016 by A. Ducharne
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_soil_smooth_over_mcs2

  SUBROUTINE hydrol_soil_smooth_over_mcs2(mcs, kjpindex, ins, soiltile, rudr_corr, check)

    !- arguments

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                           :: kjpindex        !! Domain size
    INTEGER(i_std), INTENT(in)                           :: ins             !! Soiltile index (1-nstm, unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: mcs             !! Saturated volumetric water content (m^{3} m^{-3})
    REAL(r_std), DIMENSION(kjpindex,nstm), INTENT (in)   :: soiltile    !! Fraction of each soil tile within vegtot (0-1, unitless)
    
    !! 0.2 Output variables

    REAL(r_std), DIMENSION(kjpindex,nstm), INTENT(out)   :: check           !! delta SM - flux
    
    !! 0.3 Modified variables   
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT(inout):: rudr_corr         !! Surface runoff produced to correct excess (mm/dtstep)

    !! 0.4 Local variables

    INTEGER(i_std)                        :: ji,jsl
    REAL(r_std), DIMENSION(kjpindex,nslm) :: excess
    REAL(r_std), DIMENSION(kjpindex)      :: tmci    !! total SM at beginning of routine
    REAL(r_std), DIMENSION(kjpindex)      :: tmcf    !! total SM at end of routine
    REAL(r_std)                           :: mcs_tmp

!_ ================================================================================================================================       
    !-

    !! 1. We calculate the total SM at the beginning of the routine
    IF (check_cwrr) THEN
       tmci(:) = dz(2) * ( trois*mc(:,1,ins) + mc(:,2,ins) )/huit
       DO jsl = 2,nslm-1
          tmci(:) = tmci(:) + dz(jsl) * (trois*mc(:,jsl,ins)+mc(:,jsl-1,ins))/huit &
               + dz(jsl+1) * (trois*mc(:,jsl,ins)+mc(:,jsl+1,ins))/huit
       ENDDO
       tmci(:) = tmci(:) + dz(nslm) * (trois*mc(:,nslm,ins) + mc(:,nslm-1,ins))/huit
    ENDIF  

    !! 2. In case of over-saturation, we don't do any smoothing,
    !! but directly eliminate the excess as runoff (via rudr_corr)
    !    we correct the calculation of the adjustement flux, which needs to account for nodes n-1 and n+1  
    !    for the calculation to remain simple and accurate, we directly drain all the oversaturated mc,
    !    without transfering to lower layers        

    !! 2.1 thresholding from top to bottom, with excess defined along jsl
    DO jsl = 1, nslm
       DO ji=1, kjpindex
          ! 
          ! add agri_peat also. Tides not considered here at this moment.
          mcs_tmp = mcs(ji)
          IF ( (ok_peat_hydro) .AND. ( (ins .EQ. 4) .OR. (ins .EQ. 5)) ) THEN
              IF (soiltile(ji,ins) .GT. min_sechiba) THEN
                 mcs_tmp = mcs_peat
              ELSE
                 mcs_tmp = mcs(ji)
              ENDIF
          ENDIF
          excess(ji,jsl) = MAX(mc(ji,jsl,ins)-mcs_tmp,zero)
          mc(ji,jsl,ins) = mc(ji,jsl,ins) - excess(ji,jsl) 
       ENDDO
    ENDDO

    !! 2.2 To ensure conservation, this needs to be balanced by additional drainage (in kg/m2/dt)                        
    DO ji = 1, kjpindex
       rudr_corr(ji,ins) = dz(2) * ( trois*excess(ji,1) + excess(ji,2) )/huit ! top layer = initialisation  
    ENDDO
    DO jsl = 2,nslm-1 ! intermediate layers      
       DO ji = 1, kjpindex
          rudr_corr(ji,ins) = rudr_corr(ji,ins) + dz(jsl) &
               & * (trois*excess(ji,jsl)+excess(ji,jsl-1))/huit &
               & + dz(jsl+1) * (trois*excess(ji,jsl)+excess(ji,jsl+1))/huit
       ENDDO
    ENDDO
    DO ji = 1, kjpindex
       rudr_corr(ji,ins) = rudr_corr(ji,ins) + dz(nslm) &    ! bottom layer 
            & * (trois * excess(ji,nslm) + excess(ji,nslm-1))/huit
    END DO

    !! 3. For water conservation checks, we calculate the total SM at the beginning of the routine,
    !!    and export the difference with the flux

    IF (check_cwrr) THEN
       tmcf(:) = dz(2) * ( trois*mc(:,1,ins) + mc(:,2,ins) )/huit
       DO jsl = 2,nslm-1
          tmcf(:) = tmcf(:) + dz(jsl) * (trois*mc(:,jsl,ins)+mc(:,jsl-1,ins))/huit &
               + dz(jsl+1) * (trois*mc(:,jsl,ins)+mc(:,jsl+1,ins))/huit
       ENDDO
       tmcf(:) = tmcf(:) + dz(nslm) * (trois*mc(:,nslm,ins) + mc(:,nslm-1,ins))/huit
       ! Normally, tcmf=tmci-rudr_corr
       check(:,ins) = tmcf(:)-(tmci(:)-rudr_corr(:,ins))
    ENDIF
    
  END SUBROUTINE hydrol_soil_smooth_over_mcs2


!! ================================================================================================================================
!! SUBROUTINE   : hydrol_diag_soil_flux
!!
!>\BRIEF        : This subroutine diagnoses the vertical liquid water fluxes between the 
!!                different soil layers, based on each layer water budget. It also checks the
!!                corresponding water conservation (during redistribution).
!!
!! DESCRIPTION  :
!! 1. Initialize qflux_ns from the bottom, with dr_ns
!! 2. Between layer nslm and nslm-1, by means of water budget knowing mc changes and flux at the lowest interface
!! 3. We go up, and deduct qflux_ns(1:nslm-2), still by means of water budget
!! 4. Water balance verification: pursuing upward water budget, the flux at the surface should equal -flux_top  
!!
!! RECENT CHANGE(S) : 2016 by A. Ducharne to fit hydrol_soil
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE hydrol_diag_soil_flux(kjpindex,ins,mclint,flux_top)
    !
    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                         :: kjpindex        !! Domain size
    INTEGER(i_std), INTENT(in)                         :: ins             !! index of soil type
    REAL(r_std), DIMENSION (kjpindex,nslm), INTENT(in) :: mclint          !! mc values at the beginning of the time step
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)      :: flux_top        !! Exfiltration (bare soil evaporation minus infiltration)
    
    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables
    REAL(r_std), DIMENSION (kjpindex)                  :: check_temp      !! Diagnosed flux at soil surface, should equal -flux_top
    INTEGER(i_std)                                     :: jsl,ji

    !_ ================================================================================================================================

    !- Compute the diffusion flux at every level from bottom to top (using mcl,mclint, and sink values)
    DO ji = 1, kjpindex

       !! 1. Initialize qflux_ns from the bottom, with dr_ns
       jsl = nslm
       qflux_ns(ji,jsl,ins) = dr_ns(ji,ins)
       !! 2. Between layer nslm and nslm-1, by means of water budget 
       !!    knowing mc changes and flux at the lowest interface
       !     qflux_ns is downward
       jsl = nslm-1
       qflux_ns(ji,jsl,ins) = qflux_ns(ji,jsl+1,ins) & 
            &  + (mcl(ji,jsl,ins)-mclint(ji,jsl) &
            &  + trois*mcl(ji,jsl+1,ins) - trois*mclint(ji,jsl+1)) &
            &  * (dz(jsl+1)/huit) &
            &  + rootsink(ji,jsl+1,ins) 
    ENDDO

    !! 3. We go up, and deduct qflux_ns(1:nslm-2), still by means of water budget
    ! Here, qflux_ns(ji,1,ins) is the downward flux between the top soil layer and the 2nd one
    DO jsl = nslm-2,1,-1
       DO ji = 1, kjpindex
          qflux_ns(ji,jsl,ins) = qflux_ns(ji,jsl+1,ins) & 
               &  + (mcl(ji,jsl,ins)-mclint(ji,jsl) &
               &  + trois*mcl(ji,jsl+1,ins) - trois*mclint(ji,jsl+1)) &
               &  * (dz(jsl+1)/huit) &
               &  + rootsink(ji,jsl+1,ins) &
               &  + (dz(jsl+2)/huit) &
               &  * (trois*mcl(ji,jsl+1,ins) - trois*mclint(ji,jsl+1) &
               &  + mcl(ji,jsl+2,ins)-mclint(ji,jsl+2)) 
       END DO
    ENDDO
    
    !! 4. Water balance verification: pursuing upward water budget, the flux at the surface (check_temp) 
    !! should equal -flux_top
    DO ji = 1, kjpindex

       check_temp(ji) =  qflux_ns(ji,1,ins) + (dz(2)/huit) &
            &  * (trois* (mcl(ji,1,ins)-mclint(ji,1)) + (mcl(ji,2,ins)-mclint(ji,2))) &
            &  + rootsink(ji,1,ins)    
       ! flux_top is positive when upward, while check_temp is positive when downward
       check_top_ns(ji,ins) = flux_top(ji)+check_temp(ji)

       IF (ABS(check_top_ns(ji,ins))/dt_sechiba .GT. min_sechiba) THEN
          ! Diagnosed (check_temp) and imposed (flux_top) differ by more than 1.e-8 mm/s 
          WRITE(numout,*) 'Problem in the water balance, qflux_ns computation, surface fluxes', flux_top(ji),check_temp(ji)
          WRITE(numout,*) 'Diagnosed and imposed fluxes differ by more than 1.e-8 mm/s: ', check_top_ns(ji,ins)
          WRITE(numout,*) 'ji', ji, 'jsl',jsl,'ins',ins
          WRITE(numout,*) 'mclint', mclint(ji,:)
          WRITE(numout,*) 'mcl', mcl(ji,:,ins)
          WRITE (numout,*) 'rootsink', rootsink(ji,1,ins)
          CALL ipslerr_p(1, 'hydrol_diag_soil_flux', 'NOTE:',&
               & 'Problem in the water balance, qflux_ns computation','')
       ENDIF
    ENDDO

  END SUBROUTINE hydrol_diag_soil_flux


!! ================================================================================================================================
!! SUBROUTINE   : hydrol_soil_tridiag
!!
!>\BRIEF        This subroutine solves a set of linear equations which has a tridiagonal coefficient matrix. 
!!
!! DESCRIPTION  : It is only applied in the grid-cells where resolv(ji)=TRUE
!! 
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : mcl (global module variable)
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_soil_tridiag 

  SUBROUTINE hydrol_soil_tridiag(kjpindex,ins)

    !- arguments

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                         :: kjpindex        !! Domain size
    INTEGER(i_std), INTENT(in)                         :: ins             !! number of soil type

    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                                     :: ji,jsl
    REAL(r_std), DIMENSION(kjpindex)                   :: bet
    REAL(r_std), DIMENSION(kjpindex,nslm)              :: gam

!_ ================================================================================================================================
    DO ji = 1, kjpindex

       IF (resolv(ji)) THEN
          bet(ji) = tmat(ji,1,2)
          mcl(ji,1,ins) = rhs(ji,1)/bet(ji)
       ENDIF
    ENDDO

    DO jsl = 2,nslm
       DO ji = 1, kjpindex
          
          IF (resolv(ji)) THEN

             gam(ji,jsl) = tmat(ji,jsl-1,3)/bet(ji)
             bet(ji) = tmat(ji,jsl,2) - tmat(ji,jsl,1)*gam(ji,jsl)
             mcl(ji,jsl,ins) = (rhs(ji,jsl)-tmat(ji,jsl,1)*mcl(ji,jsl-1,ins))/bet(ji)
          ENDIF

       ENDDO
    ENDDO

    DO ji = 1, kjpindex
       IF (resolv(ji)) THEN
          DO jsl = nslm-1,1,-1
             mcl(ji,jsl,ins) = mcl(ji,jsl,ins) - gam(ji,jsl+1)*mcl(ji,jsl+1,ins)
          ENDDO
       ENDIF
    ENDDO

  END SUBROUTINE hydrol_soil_tridiag


!! ================================================================================================================================
!! SUBROUTINE   : hydrol_soil_coef
!!
!>\BRIEF        Computes coef for the linearised hydraulic conductivity 
!! k_lin=a_lin mc_lin+b_lin and the linearised diffusivity d_lin. 
!!
!! DESCRIPTION  :
!! First, we identify the interval i in which the current value of mc is located.
!! Then, we give the values of the linearized parameters to compute 
!! conductivity and diffusivity as K=a*mc+b and d.
!!
!! RECENT CHANGE(S) : Addition of the dependence to profil_froz_hydro_ns 
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_soil_coef
 
  SUBROUTINE hydrol_soil_coef(mcr, mcs, kjpindex,ins, soiltile)

    IMPLICIT NONE
    !
    !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                        :: kjpindex         !! Domain size
    INTEGER(i_std), INTENT(in)                        :: ins              !! Index of soil type
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)     :: mcr              !! Residual volumetric water content (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)     :: mcs              !! Saturated volumetric water content (m^{3} m^{-3})
    REAL(r_std), DIMENSION(kjpindex,nstm), INTENT (in) :: soiltile        !! Fraction of each soil tile within vegtot (0-1, unitless)

    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                                    :: jsl,ji,i
    REAL(r_std)                                       :: mc_ratio
    REAL(r_std)                                       :: mc_used    !! Used liquid water content
    REAL(r_std)                                       :: x,m
    
!_ ================================================================================================================================

    IF (ok_freeze_cwrr) THEN
    
       ! Calculation of liquid and frozen saturation degrees with respect to residual
       ! x=liquid saturation degree/residual=(mcl-mcr)/(mcs-mcr)
       ! 1-x=frozen saturation degree/residual=(mcfc-mcr)/(mcs-mcr) (=profil_froz_hydro)
       
       DO jsl=1,nslm
          DO ji=1,kjpindex
             
             x = 1._r_std - profil_froz_hydro_ns(ji, jsl,ins)
             ! 
             IF ((ok_peat_hydro)) THEN
                IF ((ins .GE. 4) .AND. (ins .LE. 6) .AND. (soiltile(ji,ins) .GT. min_sechiba) ) THEN
                   mc_used = mcr_peat+x*MAX((mc(ji,jsl, ins)-mcr_peat),zero)
                   i= MAX(imin, MIN(imax-1,INT(imin+(imax-imin)*(mc_used-mcr_peat)/(mcs_peat-mcr_peat))))
                   a(ji,jsl) = a_lin_peat(i,jsl) * kfact_root(ji,jsl,ins)
                   b(ji,jsl) = b_lin_peat(i,jsl) * kfact_root(ji,jsl,ins)
                   d(ji,jsl) = d_lin_peat(i,jsl) * kfact_root(ji,jsl,ins)
                   k(ji,jsl) = MAX(k_lin_peat(imin+1,jsl), &
                        a_lin_peat(i,jsl) * mc_used + b_lin_peat(i,jsl))
                 ELSE 
                    mc_used = mcr(ji)+x*MAX((mc(ji,jsl, ins)-mcr(ji)),zero)
                    !
                    ! calcul de k based on mc_liq
                    ! k = kfrac * (a*theta + b)
                    i= MAX(imin, MIN(imax-1, INT(imin +(imax-imin)*(mc_used-mcr(ji))/(mcs(ji)-mcr(ji)))))
                    a(ji,jsl) = a_lin(i,jsl,ji) * kfact_root(ji,jsl,ins) ! in mm/d
                    b(ji,jsl) = b_lin(i,jsl,ji) * kfact_root(ji,jsl,ins) ! in mm/d
                    d(ji,jsl) = d_lin(i,jsl,ji) * kfact_root(ji,jsl,ins) ! in mm^2/d
                    k(ji,jsl) = kfact_root(ji,jsl,ins) * MAX(k_lin(imin+1,jsl,ji), &
                         a_lin(i,jsl,ji) * mc_used + b_lin(i,jsl,ji)) ! in mm/d
                ENDIF
             ELSE  
                ! mc_used is used in the calculation of hydrological properties
                ! It corresponds to a liquid mc, but the expression is different from mcl in hydrol_soil,
                ! to ensure that we get the a, b, d of the first bin when mcl<mcr
                mc_used = mcr(ji)+x*MAX((mc(ji,jsl, ins)-mcr(ji)),zero)
             
                ! calcul de k based on mc_liq
                ! k = kfrac * (a*theta + b)
                i= MAX(imin, MIN(imax-1, INT(imin +(imax-imin)*(mc_used-mcr(ji))/(mcs(ji)-mcr(ji)))))
                a(ji,jsl) = a_lin(i,jsl,ji) * kfact_root(ji,jsl,ins) ! in mm/d
                b(ji,jsl) = b_lin(i,jsl,ji) * kfact_root(ji,jsl,ins) ! in mm/d
                d(ji,jsl) = d_lin(i,jsl,ji) * kfact_root(ji,jsl,ins) ! in mm^2/d
                k(ji,jsl) = kfact_root(ji,jsl,ins) * MAX(k_lin(imin+1,jsl,ji), &
                     a_lin(i,jsl,ji) * mc_used + b_lin(i,jsl,ji)) ! in mm/d
             ENDIF
          ENDDO ! loop on grid
       ENDDO
       
             
    ELSE
       ! .NOT. ok_freeze_cwrr
       DO jsl=1,nslm
          DO ji=1,kjpindex 
             IF ((ok_peat_hydro) .AND. ((ins .GE. 4) .AND. (ins .LE. 6))) THEN
                IF (soiltile(ji,ins) .GT. min_sechiba) THEN
                   mc_ratio = MAX(mc(ji,jsl,ins)-mcr_peat,zero)/(mcs_peat-mcr_peat)
                   i= MAX(MIN(INT((imax-imin)*mc_ratio)+imin , imax-1), imin)
                   a(ji,jsl) = a_lin_peat(i,jsl) * kfact_root(ji,jsl,ins) 
                   b(ji,jsl) = b_lin_peat(i,jsl) * kfact_root(ji,jsl,ins)
                   d(ji,jsl) = d_lin_peat(i,jsl) * kfact_root(ji,jsl,ins)
                   k(ji,jsl) = MAX(k_lin_peat(imin+1,jsl), &
                        a_lin_peat(i,jsl) *  mc(ji,jsl,ins)+ b_lin_peat(i,jsl))
                ELSE
                   mc_ratio = MAX(mc(ji,jsl,ins)-mcr(ji), zero)/(mcs(ji)-mcr(ji))
                
                   i= MAX(MIN(INT((imax-imin)*mc_ratio)+imin , imax-1), imin)
                   a(ji,jsl) = a_lin(i,jsl,ji) * kfact_root(ji,jsl,ins) ! in mm/d
                   b(ji,jsl) = b_lin(i,jsl,ji) * kfact_root(ji,jsl,ins) ! in mm/d
                   d(ji,jsl) = d_lin(i,jsl,ji) * kfact_root(ji,jsl,ins) ! in mm^2/d
                   k(ji,jsl) = kfact_root(ji,jsl,ins) * MAX(k_lin(imin+1,jsl,ji), &
                        a_lin(i,jsl,ji) * mc(ji,jsl,ins) + b_lin(i,jsl,ji))  ! in mm/d
                ENDIF
             ELSE  
                ! it is impossible to consider a mc<mcr for the binning
                mc_ratio = MAX(mc(ji,jsl,ins)-mcr(ji), zero)/(mcs(ji)-mcr(ji))
                
                i= MAX(MIN(INT((imax-imin)*mc_ratio)+imin , imax-1), imin)
                a(ji,jsl) = a_lin(i,jsl,ji) * kfact_root(ji,jsl,ins) ! in mm/d
                b(ji,jsl) = b_lin(i,jsl,ji) * kfact_root(ji,jsl,ins) ! in mm/d
                d(ji,jsl) = d_lin(i,jsl,ji) * kfact_root(ji,jsl,ins) ! in mm^2/d
                k(ji,jsl) = kfact_root(ji,jsl,ins) * MAX(k_lin(imin+1,jsl,ji), &
                     a_lin(i,jsl,ji) * mc(ji,jsl,ins) + b_lin(i,jsl,ji))  ! in mm/d
             ENDIF
          END DO 
       END DO
    ENDIF
    
  END SUBROUTINE hydrol_soil_coef

!! ================================================================================================================================
!! SUBROUTINE   : hydrol_soil_froz
!!
!>\BRIEF        Computes profil_froz_hydro_ns, the fraction of frozen water in the soil layers. 
!!
!! DESCRIPTION  :
!!
!! RECENT CHANGE(S) : Created by A. Ducharne in 2016.
!!
!! MAIN OUTPUT VARIABLE(S) : profil_froz_hydro_ns
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_soil_froz
 
  SUBROUTINE hydrol_soil_froz(nvan, avan, mcr, mcs, kjpindex,ins,stempdiag, soiltile)

    IMPLICIT NONE
    !
    !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                        :: kjpindex         !! Domain size
    INTEGER(i_std), INTENT(in)                        :: ins              !! Index of soil type
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)     :: nvan             !! Van Genuchten coeficients n (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)     :: avan             !! Van Genuchten coeficients a (mm-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)     :: mcr              !! Residual volumetric water content (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)     :: mcs              !! Saturated volumetric water content (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (in):: stempdiag        !! Diagnostic temp profile from thermosoil
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT(in):: soiltile         !! Fraction of each soil tile within vegtot (0-1, unitless)


    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                                    :: jsl,ji,i
    REAL(r_std)                                       :: x,m
    REAL(r_std)                                       :: denom
    REAL(r_std),DIMENSION (kjpindex)                  :: froz_frac_moy
    REAL(r_std),DIMENSION (kjpindex)                  :: smtot_moy
    REAL(r_std),DIMENSION (kjpindex,nslm)             :: mc_ns
    REAL(r_std)                                       :: mcs_tmp, mcr_tmp
    
!_ ================================================================================================================================

!    ONLY FOR THE (ok_freeze_cwrr) CASE
    
       ! Calculation of liquid and frozen saturation degrees above residual moisture
       !   x=liquid saturation degree/residual=(mcl-mcr)/(mcs-mcr)
       !   1-x=frozen saturation degree/residual=(mcfc-mcr)/(mcs-mcr) (=profil_froz_hydro)
       ! It's important for the good work of the water diffusion scheme (tridiag) that the total
       ! liquid water also includes mcr, so mcl > 0 even when x=0
       
       DO jsl=1,nslm
          DO ji=1,kjpindex
             ! Van Genuchten parameter for thermodynamical calculation
             ! peatland
             mcs_tmp = mcs(ji) 
             mcr_tmp = mcr(ji) 
             m = 1. -1./nvan(ji)
             IF ((ok_peat_hydro) .AND. ((ins .GE. 4) .AND. (ins .LE. 6))) THEN
                IF (soiltile(ji,ins) .GT. min_sechiba) THEN
                   m = 1. -1./nvan_peat
                   mcs_tmp = mcs_peat
                   mcr_tmp = mcr_peat
                ENDIF
             ENDIF
                
             IF ((.NOT. ok_thermodynamical_freezing).OR.(mc(ji,jsl, ins).LT.(mcr_tmp+min_sechiba))) THEN
                ! Linear soil freezing or soil moisture below residual
                IF (stempdiag(ji, jsl).GE.(fr_center+fr_dT/2.)) THEN
                   x=1._r_std
                ELSE IF ( (stempdiag(ji,jsl) .GE. (fr_center-fr_dT/2.)) .AND. &
                     (stempdiag(ji,jsl) .LT. (fr_center+fr_dT/2.)) ) THEN 
                   x=(stempdiag(ji, jsl)-(fr_center-fr_dT/2.))/fr_dT
                ELSE 
                   x=0._r_std
                ENDIF
             ELSE IF (ok_thermodynamical_freezing) THEN                
                IF (stempdiag(ji, jsl).GE.(fr_center+fr_dT/2.)) THEN
                   x=1._r_std
                ELSE IF ( (stempdiag(ji,jsl) .GE. (fr_center-fr_dT/2.)) .AND. &
                     (stempdiag(ji,jsl) .LT. (fr_center+fr_dT/2.)) ) THEN
                   ! Factor 2.2 from the PhD of Isabelle Gouttevin
                   x=MIN(((mcs_tmp-mcr_tmp) &
                        *((2.2*1000.*avan(ji)*(fr_center+fr_dT/2.-stempdiag(ji, jsl)) &
                        *lhf/ZeroCelsius/10.)**nvan(ji)+1.)**(-m)) / &
                        (mc(ji,jsl, ins)-mcr_tmp),1._r_std)
                ELSE
                   x=0._r_std 
                ENDIF
             ENDIF
             
             profil_froz_hydro_ns(ji,jsl,ins) = 1._r_std-x
             mc_ns(ji,jsl)=mc(ji,jsl,ins)/mcs_tmp
          ENDDO ! loop on grid
       ENDDO
    
       ! Applay correction on the frozen fraction
       ! Depends on two external parameters: froz_frac_corr and smtot_corr
       froz_frac_moy(:)=zero
       denom=zero
       DO jsl=1,nslm
          froz_frac_moy(:)=froz_frac_moy(:)+dh(jsl)*profil_froz_hydro_ns(:,jsl,ins)
          denom=denom+dh(jsl)
       ENDDO
       froz_frac_moy(:)=froz_frac_moy(:)/denom

       smtot_moy(:)=zero
       denom=zero
       DO jsl=1,nslm-1
          smtot_moy(:)=smtot_moy(:)+dh(jsl)*mc_ns(:,jsl)
          denom=denom+dh(jsl)
       ENDDO
       smtot_moy(:)=smtot_moy(:)/denom

       DO jsl=1,nslm
          profil_froz_hydro_ns(:,jsl,ins)=MIN(profil_froz_hydro_ns(:,jsl,ins)* &
                                              (froz_frac_moy(:)**froz_frac_corr)*(smtot_moy(:)**smtot_corr), max_froz_hydro)
       ENDDO

     END SUBROUTINE hydrol_soil_froz
     

!! ================================================================================================================================
!! SUBROUTINE   : hydrol_soil_setup
!!
!>\BRIEF        This subroutine computes the matrix coef.  
!!
!! DESCRIPTION  : None 
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : matrix coef
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE hydrol_soil_setup(kjpindex,ins)


    IMPLICIT NONE
    !
    !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                        :: kjpindex          !! Domain size
    INTEGER(i_std), INTENT(in)                        :: ins               !! index of soil type

    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std) :: jsl,ji
    REAL(r_std)                        :: temp3, temp4

!_ ================================================================================================================================
    !-we compute tridiag matrix coefficients (LEFT and RIGHT) 
    ! of the system to solve [LEFT]*mc_{t+1}=[RIGHT]*mc{t}+[add terms]: 
    ! e(nslm),f(nslm),g1(nslm) for the [left] vector
    ! and ep(nslm),fp(nslm),gp(nslm) for the [right] vector

    ! w_time=1 (in constantes_soil) indicates implicit computation for diffusion 
    temp3 = w_time*(dt_sechiba/one_day)/deux
    temp4 = (un-w_time)*(dt_sechiba/one_day)/deux

    ! Passage to arithmetic means for layer averages also in this subroutine : Aurelien 11/05/10

    !- coefficient for first layer
    DO ji = 1, kjpindex
       e(ji,1) = zero
       f(ji,1) = trois * dz(2)/huit  + temp3 &
            & * ((d(ji,1)+d(ji,2))/(dz(2))+a(ji,1))
       g1(ji,1) = dz(2)/(huit)       - temp3 &
            & * ((d(ji,1)+d(ji,2))/(dz(2))-a(ji,2))
       ep(ji,1) = zero
       fp(ji,1) = trois * dz(2)/huit - temp4 &
            & * ((d(ji,1)+d(ji,2))/(dz(2))+a(ji,1))
       gp(ji,1) = dz(2)/(huit)       + temp4 &
            & * ((d(ji,1)+d(ji,2))/(dz(2))-a(ji,2))
    ENDDO

    !- coefficient for medium layers

    DO jsl = 2, nslm-1
       DO ji = 1, kjpindex
          e(ji,jsl) = dz(jsl)/(huit)                        - temp3 &
               & * ((d(ji,jsl)+d(ji,jsl-1))/(dz(jsl))+a(ji,jsl-1))

          f(ji,jsl) = trois * (dz(jsl)+dz(jsl+1))/huit  + temp3 &
               & * ((d(ji,jsl)+d(ji,jsl-1))/(dz(jsl)) + &
               & (d(ji,jsl)+d(ji,jsl+1))/(dz(jsl+1)) )

          g1(ji,jsl) = dz(jsl+1)/(huit)                     - temp3 &
               & * ((d(ji,jsl)+d(ji,jsl+1))/(dz(jsl+1))-a(ji,jsl+1))

          ep(ji,jsl) = dz(jsl)/(huit)                       + temp4 &
               & * ((d(ji,jsl)+d(ji,jsl-1))/(dz(jsl))+a(ji,jsl-1))

          fp(ji,jsl) = trois * (dz(jsl)+dz(jsl+1))/huit - temp4 &
               & * ( (d(ji,jsl)+d(ji,jsl-1))/(dz(jsl)) + &
               & (d(ji,jsl)+d(ji,jsl+1))/(dz(jsl+1)) )

          gp(ji,jsl) = dz(jsl+1)/(huit)                     + temp4 &
               & *((d(ji,jsl)+d(ji,jsl+1))/(dz(jsl+1))-a(ji,jsl+1))
       ENDDO
    ENDDO

    !- coefficient for last layer
    DO ji = 1, kjpindex
       e(ji,nslm) = dz(nslm)/(huit)        - temp3 &
            & * ((d(ji,nslm)+d(ji,nslm-1)) /(dz(nslm))+a(ji,nslm-1))
       f(ji,nslm) = trois * dz(nslm)/huit  + temp3 &
            & * ((d(ji,nslm)+d(ji,nslm-1)) / (dz(nslm)) &
            & -a(ji,nslm)*(un-deux*free_drain_coef(ji,ins)))
       g1(ji,nslm) = zero
       ep(ji,nslm) = dz(nslm)/(huit)       + temp4 &
            & * ((d(ji,nslm)+d(ji,nslm-1)) /(dz(nslm))+a(ji,nslm-1))
       fp(ji,nslm) = trois * dz(nslm)/huit - temp4 &
            & * ((d(ji,nslm)+d(ji,nslm-1)) /(dz(nslm)) &
            & -a(ji,nslm)*(un-deux*free_drain_coef(ji,ins)))
       gp(ji,nslm) = zero
    ENDDO

  END SUBROUTINE hydrol_soil_setup

  
!! ================================================================================================================================
!! SUBROUTINE   : hydrol_split_soil
!!
!>\BRIEF        Splits 2d variables into 3d variables, per soiltile (_ns suffix), at the beginning of hydrol
!!              At this stage, the forcing fluxes to hydrol are transformed from grid-cell averages 
!!              to mean fluxes over vegtot=sum(soiltile)  
!!
!! DESCRIPTION  :
!! 1. Split 2d variables into 3d variables, per soiltile
!! 1.1 Throughfall
!! 1.2 Bare soil evaporation
!! 1.2.2 ae_ns new
!! 1.3 transpiration
!! 1.4 root sink
!! 2. Verification: Check if the deconvolution is correct and conserves the fluxes
!! 2.1 precisol 
!! 2.2 ae_ns and evapnu
!! 2.3 transpiration
!! 2.4 root sink
!!
!! RECENT CHANGE(S) : 2016 by A. Ducharne to match the simplification of hydrol_soil
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================


  SUBROUTINE hydrol_split_soil (kjpindex, veget_max, soiltile, vevapnu, transpir, humrel, &
       evap_bare_lim, evap_bare_lim_ns, tot_bare_soil, us, e_frac, F_absorption)

    ! 
    ! interface description

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                                    :: kjpindex
    REAL(r_std), DIMENSION (kjpindex, nvm), INTENT(in)            :: veget_max        !! max Vegetation map 
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT (in)           :: soiltile         !! Fraction of each soiltile within vegtot (0-1, unitless)
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)                :: vevapnu          !! Bare soil evaporation
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT (in)            :: transpir         !! Transpiration
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT (in)            :: humrel           !! Relative humidity
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)                 :: evap_bare_lim    !!   
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT(in)            :: evap_bare_lim_ns !!   
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)                 :: tot_bare_soil    !! Total evaporating bare soil fraction
    REAL(r_std),DIMENSION (kjpindex,nvm,nstm,nslm), INTENT(inout) :: us               !! Water stress index for transpiration
                                                                                      !! (by soil layer and PFT) (0-1, unitless)
    REAL(r_std), DIMENSION (kjpindex,nvm,nslm,nstm), INTENT(in)   :: e_frac           !! Relative humidity per layer
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(in)             :: F_absorption     !! Total root absorption (ok_hydraulic_arch = .TRUE.)

    !! 0.4 Local variables

    INTEGER(i_std)                                :: ji, jv, jsl, jst
    REAL(r_std), DIMENSION (kjpindex)             :: tmp_check1
    REAL(r_std), DIMENSION (kjpindex)             :: tmp_check2
    REAL(r_std), DIMENSION (kjpindex,nstm)        :: tmp_check3
    LOGICAL                                       :: error
!_ ================================================================================================================================
    
    !! 1. Split 2d variables into 3d variables, per soiltile
    
    ! Reminders:
    !  corr_veg_soil(:,nvm,nstm) = PFT fraction per soiltile in each grid-cell
    !      corr_veg_soil(ji,jv,jst)=veget_max(ji,jv)/soiltile(ji,jst) 
    !  soiltile(:,nstm) = fraction of vegtot covered by each soiltile (0-1, unitless) 
    !  vegtot(:) = total fraction of grid-cell covered by PFTs (fraction with bare soil + vegetation)
    !  veget_max(:,nvm) = PFT fractions of vegtot+frac_nobio 
    !  veget(:,nvm) =  fractions (of vegtot+frac_nobio) covered by vegetation in each PFT 
    !       BUT veget(:,1)=veget_max(:,1) 
    !  frac_bare(:,nvm) = fraction (of veget_max) with bare soil in each PFT
    !  tot_bare_soil(:) = fraction of grid mesh covered by all bare soil (=SUM(frac_bare*veget_max))
    !  frac_bare_ns(:,nstm) = evaporating bare soil fraction (of vegtot) per soiltile (defined in hydrol_vegupd)
    
    !! 1.1 Throughfall
    ! Transformation from precisol (flux from PFT jv in m2 of grid-mesh)
    ! to  precisol_ns (flux from contributing PFTs with another unit, in m2 of soiltile)
    precisol_ns(:,:)=zero
    DO jv=1,nvm
       DO ji=1,kjpindex
          jst=pref_soil_veg(jv)
          IF((veget_max(ji,jv).GT.min_sechiba) .AND. ((soiltile(ji,jst)*vegtot(ji)) .GT. min_sechiba)) THEN
             precisol_ns(ji,jst) = precisol_ns(ji,jst) + &
                     precisol(ji,jv) / (soiltile(ji,jst)*vegtot(ji))                
          ENDIF
       END DO
    END DO
    

    !! 1.2 Bare soil evaporation and ae_ns
    ae_ns(:,:)=zero
    DO jst=1,nstm
       DO ji=1,kjpindex
          IF(evap_bare_lim(ji).GT.min_sechiba) THEN       
             ae_ns(ji,jst) = vevapnu(ji) * evap_bare_lim_ns(ji,jst)/evap_bare_lim(ji)
          ENDIF
       ENDDO
    ENDDO

    !! 1.3 transpiration
    ! Transformation from transpir (flux from PFT jv in m2 of grid-mesh)
    ! to tr_ns (flux from contributing PFTs with another unit, in m2 of soiltile)
    ! To do next: simplify the use of humrelv(ji,jv,jst) /humrel(ji,jv), since both are equal
    tr_ns(:,:)=zero
    DO jv=1,nvm
       jst=pref_soil_veg(jv)
       DO ji=1,kjpindex
          IF ((humrel(ji,jv).GT.min_sechiba) .AND. ((soiltile(ji,jst)*vegtot(ji)) .GT.min_sechiba))THEN 
             tr_ns(ji,jst)= tr_ns(ji,jst) &
                  + transpir(ji,jv) * (humrelv(ji,jv,jst) / humrel(ji,jv)) &
                  / (soiltile(ji,jst)*vegtot(ji))
                     
             ENDIF
       END DO
    END DO

    !! 1.4 root sink
    ! Transformation from transpir (flux from PFT jv in m2 of grid-mesh)
    ! to root_sink (flux from contributing PFTs and soil layer with another unit, in m2 of soiltile)
    rootsink(:,:,:)=zero

    IF (ok_hydrol_arch)THEN

       DO jv = 1,nvm
          jst=pref_soil_veg(jv)  ! OBS jst = 1,nstm
          DO jsl=1,nslm
             DO ji=1,kjpindex
                IF (humrel(ji,jv).GT.min_sechiba .AND. ((soiltile(ji,jst)*vegtot(ji)) .GT. min_sechiba)) THEN
                   rootsink(ji,jsl,jst) = rootsink(ji,jsl,jst) &
                        + (F_absorption(ji,jv)*e_frac(ji,jv,jsl,jst)*dt_sechiba * kilo_to_unit)
                END IF
             END DO
          END DO
       END DO
      
    ELSE

       DO jv=1,nvm
          jst=pref_soil_veg(jv)
          DO jsl=1,nslm
             DO ji=1,kjpindex
                IF ((humrel(ji,jv).GT.min_sechiba) .AND. ((soiltile(ji,jst)*vegtot(ji)) .GT.min_sechiba)) THEN 
                   rootsink(ji,jsl,jst) = rootsink(ji,jsl,jst) &
                        + transpir(ji,jv) * (us(ji,jv,jst,jsl) / humrel(ji,jv)) &
                        / (soiltile(ji,jst)*vegtot(ji))                     
                   ! rootsink(ji,1,jst)=0 as us(ji,jv,jst,1)=0
                END IF
             END DO
          END DO
       END DO

    END IF ! ok_hydrol_arch

    !! 2. Verification: Check if the deconvolution is correct and conserves the fluxes (grid-cell average)

    IF (check_cwrr) THEN

       error=.FALSE.

       !! 2.1 precisol 

       tmp_check1(:)=zero
       DO jst=1,nstm
          DO ji=1,kjpindex
             tmp_check1(ji)=tmp_check1(ji) + precisol_ns(ji,jst)*soiltile(ji,jst)*vegtot(ji)
          END DO
       END DO
       
       tmp_check2(:)=zero  
       DO jv=1,nvm
          DO ji=1,kjpindex
             tmp_check2(ji)=tmp_check2(ji) + precisol(ji,jv)
          END DO
       END DO

       DO ji=1,kjpindex   
          IF(ABS(tmp_check1(ji) - tmp_check2(ji)).GT.allowed_err) THEN
             WRITE(numout,*) 'PRECISOL SPLIT FALSE:ji=',ji,tmp_check1(ji),tmp_check2(ji)
             WRITE(numout,*) 'err',ABS(tmp_check1(ji)- tmp_check2(ji))
             WRITE(numout,*) 'vegtot',vegtot(ji)
             DO jv=1,nvm
                WRITE(numout,'(a,i2.2,"|",F13.4,"|",F13.4,"|",3(F9.6))') &
                     'jv,veget_max, precisol, vegetmax_soil ', &
                     jv,veget_max(ji,jv),precisol(ji,jv),vegetmax_soil(ji,jv,:)
             END DO
             DO jst=1,nstm
                WRITE(numout,*) 'jst,precisol_ns',jst,precisol_ns(ji,jst)
                WRITE(numout,*) 'soiltile', soiltile(ji,jst)
             END DO
             error=.TRUE.
             CALL ipslerr_p(2, 'hydrol_split_soil', 'We will STOP in the end of this subroutine.',&
                  & 'check_CWRR','PRECISOL SPLIT FALSE')
          ENDIF
       END DO
       
       !! 2.2 ae_ns and evapnu

       tmp_check1(:)=zero
       DO jst=1,nstm
          DO ji=1,kjpindex
             tmp_check1(ji)=tmp_check1(ji) + ae_ns(ji,jst)*soiltile(ji,jst)*vegtot(ji)
          END DO
       END DO

       DO ji=1,kjpindex   

          IF(ABS(tmp_check1(ji) - vevapnu(ji)).GT.allowed_err) THEN
             WRITE(numout,*) 'VEVAPNU SPLIT FALSE:ji, Sum(ae_ns), vevapnu =',ji,tmp_check1(ji),vevapnu(ji)
             WRITE(numout,*) 'err',ABS(tmp_check1(ji)- vevapnu(ji))
             WRITE(numout,*) 'ae_ns',ae_ns(ji,:)
             WRITE(numout,*) 'vegtot',vegtot(ji)
             WRITE(numout,*) 'evap_bare_lim, evap_bare_lim_ns',evap_bare_lim(ji), evap_bare_lim_ns(ji,:)
             DO jst=1,nstm
                WRITE(numout,*) 'jst,ae_ns',jst,ae_ns(ji,jst)
                WRITE(numout,*) 'soiltile', soiltile(ji,jst)
             END DO
             error=.TRUE.
             CALL ipslerr_p(2, 'hydrol_split_soil', 'We will STOP in the end of this subroutine.',&
                  & 'check_CWRR','VEVAPNU SPLIT FALSE')
          ENDIF
       ENDDO

    !! 2.3 transpiration

       tmp_check1(:)=zero
       DO jst=1,nstm
          DO ji=1,kjpindex
             tmp_check1(ji)=tmp_check1(ji) + tr_ns(ji,jst)*soiltile(ji,jst)*vegtot(ji)
          END DO
       END DO
       
       tmp_check2(:)=zero  
       DO jv=1,nvm
          DO ji=1,kjpindex
             tmp_check2(ji)=tmp_check2(ji) + transpir(ji,jv)
          END DO
       END DO

       DO ji=1,kjpindex   
          IF(ABS(tmp_check1(ji)- tmp_check2(ji)).GT.allowed_err) THEN
             WRITE(numout,*) 'TRANSPIR SPLIT FALSE:ji=',ji,tmp_check1(ji),tmp_check2(ji)
             WRITE(numout,*) 'err',ABS(tmp_check1(ji)- tmp_check2(ji))
             WRITE(numout,*) 'vegtot',vegtot(ji)
             DO jv=1,nvm
                WRITE(numout,*) 'jv,veget_max, transpir',jv,veget_max(ji,jv),transpir(ji,jv)
                DO jst=1,nstm
                   WRITE(numout,*) 'vegetmax_soil:ji,jv,jst',ji,jv,jst,vegetmax_soil(ji,jv,jst)
                END DO
             END DO
             DO jst=1,nstm
                WRITE(numout,*) 'jst,tr_ns',jst,tr_ns(ji,jst)
                WRITE(numout,*) 'soiltile', soiltile(ji,jst)
             END DO
             error=.TRUE.
             CALL ipslerr_p(2, 'hydrol_split_soil', 'We will STOP in the end of this subroutine.',&
                  & 'check_CWRR','TRANSPIR SPLIT FALSE')
          ENDIF

       END DO

    !! 2.4 root sink

       tmp_check3(:,:)=zero
       DO jst=1,nstm
          DO jsl=1,nslm
             DO ji=1,kjpindex
                tmp_check3(ji,jst)=tmp_check3(ji,jst) + rootsink(ji,jsl,jst)
             END DO
          END DO
       ENDDO

       DO jst=1,nstm
          DO ji=1,kjpindex
             IF(ABS(tmp_check3(ji,jst) - tr_ns(ji,jst)).GT.allowed_err) THEN
                WRITE(numout,*) 'ROOTSINK SPLIT FALSE:ji,jst=', ji,jst,&
                     & tmp_check3(ji,jst),tr_ns(ji,jst)
                WRITE(numout,*) 'err',ABS(tmp_check3(ji,jst)- tr_ns(ji,jst))
                WRITE(numout,*) 'HUMREL(jv=1:13)',humrel(ji,:)
                WRITE(numout,*) 'TRANSPIR',transpir(ji,:)
                DO jv=1,nvm 
                   WRITE(numout,*) 'jv=',jv,'us=',us(ji,jv,jst,:)
                ENDDO
                error=.TRUE.
                CALL ipslerr_p(2, 'hydrol_split_soil', 'We will STOP in the end of this subroutine.',&
                  & 'check_CWRR','ROOTSINK SPLIT FALSE')
             ENDIF
          END DO
       END DO


       !! Exit if error was found previously in this subroutine
       IF ( error ) THEN
          WRITE(numout,*) 'One or more errors have been detected in hydrol_split_soil. Model stops.'
          CALL ipslerr_p(3, 'hydrol_split_soil', 'We will STOP now.',&
               & 'One or several fatal errors were found previously.','')
       END IF

    ENDIF ! end of check_cwrr


  END SUBROUTINE hydrol_split_soil
  

!! ================================================================================================================================
!! SUBROUTINE   : hydrol_diag_soil
!!
!>\BRIEF        Calculates diagnostic variables at the grid-cell scale
!!
!! DESCRIPTION  :
!! - 1. Apply mask_soiltile
!! - 2. Sum 3d variables in 2d variables with fraction of vegetation per soil type
!!
!! RECENT CHANGE(S) : 2016 by A. Ducharne for the claculation of shumdiag_perma
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_diag_soil

  SUBROUTINE hydrol_diag_soil (ks, nvan, avan, mcr, mcs, mcfc, mcw, kjpindex, veget_max, soiltile, runoff, drainage, &
       & evapot, vevapnu, returnflow, reinfiltration, irrigation, &
       & shumdiag,shumdiag_perma, k_litt, litterhumdiag, humrel, vegstress, drysoil_frac, tot_melt, us, &
       & precip_rain, totfrac_nobio, frac_snow_nobio, &
       & shumdiag_peat, shumdiag_croppeat, shumdiag_man)
    ! 
    ! interface description

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    ! input scalar 
    INTEGER(i_std), INTENT(in)                               :: kjpindex 
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)        :: veget_max       !! Max. vegetation type
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT (in)      :: soiltile        !! Fraction of each soil tile within vegtot (0-1, unitless)
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: evapot          !! 
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: returnflow      !! Water returning to the deep reservoir
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: reinfiltration  !! Water returning to the top of the soil
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: irrigation      !! Water from irrigation
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: tot_melt        !!
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: ks               !! Hydraulic conductivity at saturation (mm {-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: nvan             !! Van Genuchten coeficients n (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: avan             !! Van Genuchten coeficients a (mm-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: mcr              !! Residual volumetric water content (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: mcs              !! Saturated volumetric water content (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: mcfc             !! Volumetric water content at field capacity (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: mcw              !! Volumetric water content at wilting point (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: precip_rain      !! Rain precipitation				
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: totfrac_nobio    !! Total fraction of continental ice+lakes+...		
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT(in)      :: frac_snow_nobio  !! Snow cover fraction on non-vegeted area
    ! peatland
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (out)      :: shumdiag_peat    !! soil moisture content useful for the decomposition when peat is activated
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (out)      :: shumdiag_croppeat!! soil moisture content useful for the decomposition when agri_peat is activated
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (out)      :: shumdiag_man     !! soil moisture content useful for the decomposition when tides is activated

    !! 0.2 Output variables

    REAL(r_std), DIMENSION (kjpindex), INTENT (out)          :: drysoil_frac    !! Function of litter wetness
    REAL(r_std), DIMENSION (kjpindex), INTENT(out)           :: runoff          !! complete runoff
    REAL(r_std), DIMENSION (kjpindex), INTENT(out)           :: drainage        !! Drainage
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (out)      :: shumdiag        !! relative soil moisture
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (out)      :: shumdiag_perma  !! Percent of porosity filled with water (mc/mcs) used for the thermal computations
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: k_litt          !! litter cond.
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: litterhumdiag   !! litter humidity
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)       :: humrel          !! Relative humidity
    REAL(r_std), DIMENSION (kjpindex, nvm), INTENT(out)      :: vegstress       !! Veg. moisture stress (only for vegetation growth)
 
    !! 0.3 Modified variables 

    REAL(r_std), DIMENSION (kjpindex), INTENT(inout)         :: vevapnu         !!
    REAL(r_std),DIMENSION (kjpindex,nvm,nstm,nslm), INTENT(inout) :: us         !! Water stress index for transpiration
                                                                                !! (by soil layer and PFT) (0-1, unitless)
    !! 0.4 Local variables

    INTEGER(i_std)                                           :: ji, jv, jsl, jst, i
    REAL(r_std), DIMENSION (kjpindex)                        :: mask_vegtot
    REAL(r_std)                                              :: k_tmp, tmc_litter_ratio

!_ ================================================================================================================================
    !
    ! Put the prognostics variables of soil to zero if soiltype is zero

    !! 1. Apply mask_soiltile
    
    DO jst=1,nstm 
       DO ji=1,kjpindex

             ae_ns(ji,jst) = ae_ns(ji,jst) * mask_soiltile(ji,jst)
             dr_ns(ji,jst) = dr_ns(ji,jst) * mask_soiltile(ji,jst)
             ru_ns(ji,jst) = ru_ns(ji,jst) * mask_soiltile(ji,jst)
             tmc(ji,jst) =  tmc(ji,jst) * mask_soiltile(ji,jst)

             DO jv=1,nvm
                humrelv(ji,jv,jst) = humrelv(ji,jv,jst) * mask_soiltile(ji,jst)
                DO jsl=1,nslm
                   us(ji,jv,jst,jsl) = us(ji,jv,jst,jsl)  * mask_soiltile(ji,jst)
                END DO
             END DO

             DO jsl=1,nslm          
                mc(ji,jsl,jst) = mc(ji,jsl,jst)  * mask_soiltile(ji,jst)
             END DO

       END DO
    END DO

    runoff(:) = zero
    drainage(:) = zero
    humtot(:) = zero
    shumdiag(:,:)= zero
    shumdiag_perma(:,:)=zero
    k_litt(:) = zero
    litterhumdiag(:) = zero
    tmc_litt_dry_mea(:) = zero
    tmc_litt_wet_mea(:) = zero
    tmc_litt_mea(:) = zero
    humrel(:,:) = zero
    vegstress(:,:) = zero
    shumdiag_peat(:,:)=zero   
    shumdiag_croppeat(:,:)=zero
    shumdiag_man(:,:)=zero  

    IF (ok_freeze_cwrr) THEN
       profil_froz_hydro(:,:)=zero ! initialisation for the mean of profil_froz_hydro_ns
    ENDIF
    
    !! 2. Sum 3d variables in 2d variables with fraction of vegetation per soil type

    DO ji = 1, kjpindex
       mask_vegtot(ji) = 0
       IF(vegtot(ji) .GT. min_sechiba) THEN
          mask_vegtot(ji) = 1
       ENDIF
    END DO
    
    DO ji = 1, kjpindex 
       ! Here we weight ae_ns by the fraction of bare evaporating soil. 
       ! This is given by frac_bare_ns, taking into account bare soil under vegetation
       ae_ns(ji,:) = mask_vegtot(ji) * ae_ns(ji,:) * frac_bare_ns(ji,:)
    END DO

    ! We average the values of each soiltile and multiply by vegtot to transform to a grid-cell mean
    DO jst = 1, nstm
       DO ji = 1, kjpindex 
          drainage(ji) = mask_vegtot(ji) * (drainage(ji) + vegtot(ji)*soiltile(ji,jst) * dr_ns(ji,jst))
          runoff(ji) = mask_vegtot(ji) *  (runoff(ji) +   vegtot(ji)*soiltile(ji,jst) * ru_ns(ji,jst)) &
               &   + (1 - mask_vegtot(ji)) * (tot_melt(ji) + irrigation(ji) + returnflow(ji) + reinfiltration(ji))
          humtot(ji) = mask_vegtot(ji) * (humtot(ji) + vegtot(ji)*soiltile(ji,jst) * tmc(ji,jst)) 
          IF (ok_freeze_cwrr) THEN 
             !  profil_froz_hydro_ns comes from hydrol_soil, to remain the same as in the prognotic loop
             profil_froz_hydro(ji,:)=mask_vegtot(ji) * &
                  (profil_froz_hydro(ji,:) + vegtot(ji)*soiltile(ji,jst) * profil_froz_hydro_ns(ji,:, jst))
          ENDIF
       END DO
    END DO

    ! we add the excess of snow sublimation to vevapnu
    ! - because vevapsno is modified in hydrol_snow if subsinksoil
    ! - it is multiplied by vegtot because it is devided by 1-tot_frac_nobio at creation in hydrol_snow

    DO ji = 1,kjpindex
       IF (vegtot(ji) .NE. 0.) THEN
          vevapnu(ji) = vevapnu (ji) + subsinksoil(ji)*vegtot(ji)
       ELSE
	  vevapnu(ji) = vevapnu (ji) + subsinksoil(ji)
       ENDIF
       runoff(ji) = runoff(ji) + precip_rain(ji)*totfrac_nobio(ji)*(1-frac_snow_nobio(ji,iice))
    END DO

    DO jst=1,nstm
       DO jv=1,nvm
          DO ji=1,kjpindex
             IF(veget_max(ji,jv).GT.min_sechiba) THEN
                vegstress(ji,jv)=vegstress(ji,jv)+vegstressv(ji,jv,jst)
                vegstress(ji,jv)= MAX(vegstress(ji,jv),zero)
             ENDIF
          END DO
       END DO
    END DO
    DO jst=1,nstm
       DO jv=1,nvm
          DO ji=1,kjpindex
             humrel(ji,jv)=humrel(ji,jv)+humrelv(ji,jv,jst)
             humrel(ji,jv)=MAX(humrel(ji,jv),zero)
          END DO
       END DO
    END DO

    !! Litter... the goal is to calculate drysoil_frac, to calculate the albedo in condveg
    ! In condveg, drysoil_frac serve to calculate the albedo of drysoil, excluding the nobio contribution which is further added
    ! In conclusion, we calculate drysoil_frac based on moisture averages restricted to the soiltile (no multiplication by vegtot)
    ! BUT THIS IS NOT USED ANYMORE WITH THE NEW BACKGROUNG ALBEDO
    !! k_litt is calculated here as a grid-cell average (for consistency with drainage)   
    !! litterhumdiag, like shumdiag, is averaged over the soiltiles for transmission to stomate
    !
    !
    tmc_litter_ratio  = 0.0
    DO jst=1,nstm       
       DO ji=1,kjpindex
          ! We compute here a mean k for the 'litter' used for reinfiltration from floodplains of ponds        
          IF ( tmc_litter(ji,jst) < tmc_litter_res(ji,jst)) THEN  
             i = imin
          ELSE IF ((tmc_litter(ji,jst) .GE. tmc_litter_res(ji,jst)) .AND. (soiltile(ji,jst) .GT. min_sechiba)) THEN
             tmc_litter_ratio = (tmc_litter(ji,jst)-tmc_litter_res(ji,jst)) / &
                  & (tmc_litter_sat(ji,jst)-tmc_litter_res(ji,jst))
             i= MAX(MIN(INT((imax-imin)*tmc_litter_ratio)+imin, imax-1), imin)
          ELSE IF (soiltile(ji,jst) .LE. min_sechiba) THEN
             tmc_litter_ratio = 0.0
             i = imin
          ENDIF
        
          IF ((ok_peat_hydro) .AND. ((jst .GE. 4) .AND. (jst .LE. 6))) THEN
             IF (soiltile(ji,jst) .GT. min_sechiba) THEN
                k_tmp = MAX(k_lin_peat(i,1)*ks_peat, zero)
             ELSE
                k_tmp = MAX(k_lin(i,1,ji)*ks(ji), zero)
             ENDIF
          ELSE
             k_tmp = MAX(k_lin(i,1,ji)*ks(ji), zero)
          ENDIF
          k_litt(ji) = k_litt(ji) + vegtot(ji)*soiltile(ji,jst) * SQRT(k_tmp) ! grid-cell average
       ENDDO      
       DO ji=1,kjpindex
          litterhumdiag(ji) = litterhumdiag(ji) + &
               & soil_wet_litter(ji,jst) * soiltile(ji,jst)

          tmc_litt_wet_mea(ji) =  tmc_litt_wet_mea(ji) + & 
               & tmc_litter_awet(ji,jst)* soiltile(ji,jst)

          tmc_litt_dry_mea(ji) = tmc_litt_dry_mea(ji) + &
               & tmc_litter_adry(ji,jst) * soiltile(ji,jst) 

          tmc_litt_mea(ji) = tmc_litt_mea(ji) + &
               & tmc_litter(ji,jst) * soiltile(ji,jst) 
       ENDDO
    ENDDO
    
    DO ji=1,kjpindex
       IF ( tmc_litt_wet_mea(ji) - tmc_litt_dry_mea(ji) > zero ) THEN
          drysoil_frac(ji) = un + MAX( MIN( (tmc_litt_dry_mea(ji) - tmc_litt_mea(ji)) / &
               & (tmc_litt_wet_mea(ji) - tmc_litt_dry_mea(ji)), zero), - un)
       ELSE
          drysoil_frac(ji) = zero
       ENDIF
    END DO
    
    ! Calculate soilmoist, as a function of total water content (mc)
    ! We average the values of each soiltile and multiply by vegtot to transform to a grid-cell mean
    soilmoist(:,:) = zero
    DO jst=1,nstm
       DO ji=1,kjpindex
             soilmoist(ji,1) = soilmoist(ji,1) + soiltile(ji,jst) * &
                  dz(2) * ( trois*mc(ji,1,jst) + mc(ji,2,jst) )/huit
             DO jsl = 2,nslm-1
                soilmoist(ji,jsl) = soilmoist(ji,jsl) + soiltile(ji,jst) * &
                     ( dz(jsl) * (trois*mc(ji,jsl,jst)+mc(ji,jsl-1,jst))/huit &
                     + dz(jsl+1) * (trois*mc(ji,jsl,jst)+mc(ji,jsl+1,jst))/huit )
             END DO
             soilmoist(ji,nslm) = soilmoist(ji,nslm) + soiltile(ji,jst) * &
                  dz(nslm) * (trois*mc(ji,nslm,jst) + mc(ji,nslm-1,jst))/huit
       END DO
    END DO
    DO ji=1,kjpindex
       soilmoist(ji,:) = soilmoist(ji,:) * vegtot(ji) ! conversion to grid-cell average
    ENDDO

    soilmoist_s(:,:,:) = zero
    DO jst=1,nstm
       DO ji=1,kjpindex
             soilmoist_s(ji,1,nstm) = soilmoist_s(ji,1,nstm) + soiltile(ji,jst) * &
                  dz(2) * ( trois*mc(ji,1,jst) + mc(ji,2,jst) )/huit
             DO jsl = 2,nslm-1
                soilmoist_s(ji,jsl,nstm) = soilmoist_s(ji,jsl,nstm) + soiltile(ji,jst) * &
                     ( dz(jsl) * (trois*mc(ji,jsl,jst)+mc(ji,jsl-1,jst))/huit &
                     + dz(jsl+1) * (trois*mc(ji,jsl,jst)+mc(ji,jsl+1,jst))/huit )
             END DO
             soilmoist_s(ji,nslm,nstm) = soilmoist_s(ji,nslm,nstm) + soiltile(ji,jst) * &
                  dz(nslm) * (trois*mc(ji,nslm,jst) + mc(ji,nslm-1,jst))/huit
       END DO
    END DO
    DO ji=1,kjpindex
       soilmoist_s(ji,:,:) = soilmoist_s(ji,:,:) * vegtot(ji) ! conversion to grid-cell average
    ENDDO

    soilmoist_liquid(:,:) = zero
    DO jst=1,nstm
       DO ji=1,kjpindex
          soilmoist_liquid(ji,1) = soilmoist_liquid(ji,1) + soiltile(ji,jst) * &
               dz(2) * ( trois*mcl(ji,1,jst) + mcl(ji,2,jst) )/huit
          DO jsl = 2,nslm-1
             soilmoist_liquid(ji,jsl) = soilmoist_liquid(ji,jsl) + soiltile(ji,jst) * &
                  ( dz(jsl) * (trois*mcl(ji,jsl,jst)+mcl(ji,jsl-1,jst))/huit &
                  + dz(jsl+1) * (trois*mcl(ji,jsl,jst)+mcl(ji,jsl+1,jst))/huit )
          END DO
          soilmoist_liquid(ji,nslm) = soilmoist_liquid(ji,nslm) + soiltile(ji,jst) * &
               dz(nslm) * (trois*mcl(ji,nslm,jst) + mcl(ji,nslm-1,jst))/huit
       ENDDO
    ENDDO
    DO ji=1,kjpindex
        soilmoist_liquid(ji,:) = soilmoist_liquid(ji,:) * vegtot_old(ji) ! grid cell average 
    ENDDO
    
    ! peat and tides
    !Used in stomate_litter for computing litter decomposition for peatland. non relative humidity.
    IF (ok_peat_hydro) THEN
       shumdiag_peat(:,:)=zero
       DO ji=1,kjpindex
          IF (soiltile(ji,4) .GT. min_sechiba) THEN
             shumdiag_peat(ji,1)= dz(2) * ( trois*mc(ji,1,4) + mc(ji,2,4) )/huit
             DO jsl = 2,nslm-1
                shumdiag_peat(ji,jsl)= dz(jsl) * (trois*mc(ji,jsl,4)+mc(ji,jsl-1,4))/huit &
                     + dz(jsl+1) * (trois*mc(ji,jsl,4)+mc(ji,jsl+1,4))/huit
             ENDDO
             shumdiag_peat(ji,nslm)= dz(nslm) * (trois*mc(ji,nslm,4) + mc(ji,nslm-1,4))/huit
          ENDIF
       ENDDO
    ENDIF

    IF (ok_peat_hydro) THEN
       DO jsl=1,nslm
          DO ji=1,kjpindex
             IF ( soiltile(ji,4) .GT. min_sechiba ) THEN
                shumdiag_peat(ji,jsl)=shumdiag_peat(ji,jsl)/(dh(jsl)) 
                shumdiag_peat(ji,jsl)=MAX(MIN(shumdiag_peat(ji,jsl), un), zero)
             ENDIF
          ENDDO
       ENDDO
    ENDIF
    
    IF (tides) THEN
       shumdiag_man(:,:)=zero
       DO ji=1,kjpindex
          shumdiag_man(ji,1)= dz(2) * ( trois*mc(ji,1,6) + mc(ji,2,6) )/huit
          DO jsl = 2,nslm-1
             shumdiag_man(ji,jsl)= dz(jsl) * (trois*mc(ji,jsl,6)+mc(ji,jsl-1,6))/huit &
                  + dz(jsl+1) * (trois*mc(ji,jsl,6)+mc(ji,jsl+1,6))/huit
          ENDDO
          shumdiag_man(ji,nslm)= dz(nslm) * (trois*mc(ji,nslm,6) + mc(ji,nslm-1,6))/huit
       ENDDO
    ENDIF
        
    IF (agri_peat) THEN
       shumdiag_croppeat(:,:)=zero
       DO ji=1,kjpindex
          shumdiag_croppeat(ji,1)= dz(2) * ( trois*mc(ji,1,5) + mc(ji,2,5) )/huit
          DO jsl = 2,nslm-1
             shumdiag_croppeat(ji,jsl)=dz(jsl) * (trois*mc(ji,jsl,5)+mc(ji,jsl-1,5))/huit &
                  + dz(jsl+1) * (trois*mc(ji,jsl,5)+mc(ji,jsl+1,5))/huit
          ENDDO
          shumdiag_croppeat(ji,nslm)= dz(nslm) * (trois*mc(ji,nslm,5) + mc(ji,nslm-1,5))/huit
       ENDDO
    ENDIF

    IF (tides) THEN
       DO jsl=1,nslm
          DO ji=1,kjpindex
             IF ( soiltile(ji,6) .GT. min_sechiba ) THEN
                shumdiag_man(ji,jsl)=shumdiag_man(ji,jsl)/(dh(jsl)) 
                shumdiag_man(ji,jsl)=MAX(MIN(shumdiag_man(ji,jsl), un), zero)
             ENDIF
          ENDDO
       ENDDO
    ENDIF

    IF (agri_peat) THEN
       DO jsl=1,nslm
          DO ji=1,kjpindex
             IF ( soiltile(ji,5) .GT. min_sechiba ) THEN
                shumdiag_croppeat(ji,jsl)=shumdiag_croppeat(ji,jsl)/(dh(jsl)) 
                shumdiag_croppeat(ji,jsl)=MAX(MIN(shumdiag_croppeat(ji,jsl), un), zero)
             ENDIF
          ENDDO
       ENDDO
    ENDIF
        
    ! Shumdiag: we start from soil_wet_ns, but then we change the range via ((mcs_peat-mcw_peat)/(mcf_peat-mcw_peat))
    ! then do a spatial average, excluding the nobio fraction on which stomate doesn't act
    ! Used in stomate_litter for computing litter decomposition for nonpeat 
    ! 
    DO jst=1,nstm      
       DO jsl=1,nslm
          DO ji=1,kjpindex
             IF ( (ok_peat_hydro) .AND. ((jst .GE. 4) .AND. (jst .LE. 6))) THEN
                IF (soiltile (ji,jst) .GT. min_sechiba)  THEN 
                   ! average between the peat and non peat
                   shumdiag(ji,jsl) = shumdiag(ji,jsl) + soil_wet_ns(ji,jsl,jst) * &
                        soiltile(ji,jst)*((mcs_peat-mcw_peat)/(mcf_peat-mcw_peat))
                ELSE
                   shumdiag(ji,jsl) = shumdiag(ji,jsl) + soil_wet_ns(ji,jsl,jst) * soiltile(ji,jst) * &
                               ((mcs(ji)-mcw(ji))/(mcfc(ji)-mcw(ji)))
                ENDIF
             ELSE
                shumdiag(ji,jsl) = shumdiag(ji,jsl) + soil_wet_ns(ji,jsl,jst) * soiltile(ji,jst) * &
                               ((mcs(ji)-mcw(ji))/(mcfc(ji)-mcw(ji)))
             ENDIF
             !
             shumdiag(ji,jsl) = MAX(MIN(shumdiag(ji,jsl), un), zero) 
          ENDDO
       ENDDO
    ENDDO
    
    ! Shumdiag_perma is based on soilmoist / moisture at saturation in the layer
    ! Her we start from grid averages by hydrol soil layer and transform it to the diag levels
    ! We keep a grid-cell average, like for all variables transmitted to ok_freeze
    ! Finally Shumdiag_perma is relative humidity , used in thermosol code to compute paca.
    DO jsl=1,nslm              
       DO ji=1,kjpindex
          IF ( (ok_peat_hydro) .AND. (jst .EQ. 4)) THEN
             IF (soiltile (ji,4) .GT. min_sechiba) THEN 
                ! average between the peat and non peat
                ! soilmoist is water content scaled by the soil tile, and it is not relative humidty
                shumdiag_perma(ji,jsl) = soilmoist(ji,jsl)/(dh(jsl) *(mcs(ji)*(un-soiltile(ji,4))+mcs_peat*soiltile(ji,4)))
             ELSE
                shumdiag_perma(ji,jsl) = soilmoist(ji,jsl) / (dh(jsl)*mcs(ji))
             ENDIF
           ELSE
             shumdiag_perma(ji,jsl) = soilmoist(ji,jsl) / (dh(jsl)*mcs(ji))
          ENDIF
          shumdiag_perma(ji,jsl) = MAX(MIN(shumdiag_perma(ji,jsl), un), zero) 
       ENDDO
    ENDDO

  END SUBROUTINE hydrol_diag_soil  


!! ================================================================================================================================
!! SUBROUTINE   : hydrol_alma 
!!
!>\BRIEF        This routine computes the changes in soil moisture and interception storage for the ALMA outputs.  
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_alma

  SUBROUTINE hydrol_alma (kjpindex, index, lstep_init, qsintveg, snow, soilwet)
    !
    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT (in)                        :: kjpindex     !! Domain size
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)   :: index        !! Indeces of the points on the map
    LOGICAL, INTENT (in)                               :: lstep_init   !! At which time is this routine called ?
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: qsintveg     !! Water on vegetation due to interception
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: snow         !! Snow water equivalent

    !! 0.2 Output variables

    REAL(r_std),DIMENSION (kjpindex), INTENT (out)      :: soilwet     !! Soil wetness

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std) :: ji
    REAL(r_std) :: watveg

!_ ================================================================================================================================
    !
    !
    IF ( lstep_init ) THEN
       ! Initialize variables if they were not found in the restart file

       DO ji = 1, kjpindex
          watveg = SUM(qsintveg(ji,:))
          tot_watveg_beg(ji) = watveg
          tot_watsoil_beg(ji) = humtot(ji)
          snow_beg(ji)        = snow(ji)

       ENDDO

       RETURN

    ENDIF
    !
    ! Calculate the values for the end of the time step
    !
    DO ji = 1, kjpindex
       watveg = SUM(qsintveg(ji,:)) ! average within the mesh
       tot_watveg_end(ji) = watveg
       tot_watsoil_end(ji) = humtot(ji) ! average within the mesh
       snow_end(ji) = snow(ji)          ! average within the mesh
       delintercept(ji) = tot_watveg_end(ji) - tot_watveg_beg(ji) ! average within the mesh
       delsoilmoist(ji) = tot_watsoil_end(ji) - tot_watsoil_beg(ji)
       delswe(ji)       = snow_end(ji) - snow_beg(ji) ! average within the mesh
       
    ENDDO
    !
    !
    ! Transfer the total water amount at the end of the current timestep top the begining of the next one.
    !
    tot_watveg_beg = tot_watveg_end
    tot_watsoil_beg = tot_watsoil_end
    snow_beg(:) = snow_end(:)
    !
    DO ji = 1,kjpindex
       IF ( mx_eau_var(ji) > 0 ) THEN
          soilwet(ji) = tot_watsoil_end(ji) / mx_eau_var(ji)
       ELSE
          soilwet(ji) = zero
       ENDIF
    ENDDO
    !
  END SUBROUTINE hydrol_alma


!! ================================================================================================================================
!! SUBROUTINE   : hydrol_nudge_mc_read
!!
!>\BRIEF         Read soil moisture from file and interpolate to the current time step
!!
!! DESCRIPTION  : Nudging of soil moisture and/or snow variables is done if OK_NUDGE_MC=y and/or OK_NUDGE_SNOW=y in run.def. 
!!                This subroutine reads and interpolates spatialy if necessary and temporary the soil moisture from file. 
!!                The values for the soil moisture will be applaied later using hydrol_nudge_mc
!!
!! RECENT CHANGE(S) : None
!!
!! \n
!_ ================================================================================================================================

  SUBROUTINE hydrol_nudge_mc_read(kjit)

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                         :: kjit         !! Timestep number

    !! 0.3 Locals variables
    REAL(r_std)                                :: tau                   !! Position between to values in nudge mc file
    REAL(r_std), DIMENSION(iim_g,jjm_g,nslm,1) :: mc_read_glo2D_1       !! mc from file at global 2D(lat,lon) grid per soiltile
    REAL(r_std), DIMENSION(iim_g,jjm_g,nslm,1) :: mc_read_glo2D_2       !! mc from file at global 2D(lat,lon) grid per soiltile
    REAL(r_std), DIMENSION(iim_g,jjm_g,nslm,1) :: mc_read_glo2D_3       !! mc from file at global 2D(lat,lon) grid per soiltile
    REAL(r_std), DIMENSION(nbp_glo,nslm,nstm)  :: mc_read_glo1D         !! mc_read_glo2D on land-only vector form, in global
    INTEGER(i_std), SAVE                       :: istart_mc !! start index to read from input file
!$OMP THREADPRIVATE(istart_mc)
    INTEGER(i_std)                             :: iend                  !! end index to read from input file
    INTEGER(i_std)                             :: i, j, ji, jg, jst, jsl!! loop index
    INTEGER(i_std)                             :: iim_file, jjm_file, llm_file !! Dimensions in input file
    INTEGER(i_std), SAVE                       :: ttm_mc                !! Time dimensions in input file
!$OMP THREADPRIVATE(ttm_mc)
    INTEGER(i_std), SAVE                       :: mc_id                 !! index for netcdf files
!$OMP THREADPRIVATE(mc_id)
    LOGICAL, SAVE                              :: firsttime_mc=.TRUE.
!$OMP THREADPRIVATE(firsttime_mc)

 
    !! 1. Nudging of soil moisture

       !! 1.2 Read mc from file, once a day only
       !!     The forcing file must contain daily frequency variable for the full year of the simulation
       IF (MOD(kjit,INT(one_day/dt_sechiba)) == 1) THEN 
          ! Save mc read from file from previous day
          mc_read_prev = mc_read_next

          IF (nudge_interpol_with_xios) THEN
             ! Read mc from input file. XIOS interpolates it to the model grid before it is received here.
             CALL xios_orchidee_recv_field("moistc_interp", mc_read_next)

             ! Read and interpolation the mask for variable mc from input file. 
             ! This is only done to be able to output the mask it later for validation purpose.
             ! The mask corresponds to the fraction of the input source file which was underlaying the model grid cell. 
             ! If the msask is 0 for a model grid cell, then the default value 0.2 set in field_def_orchidee.xml, is used for that grid cell. 
             CALL xios_orchidee_recv_field("mask_moistc_interp", mask_mc_interp)

          ELSE

             ! Only read fields from the file. We here suppose that no interpolation is needed.
             IF (is_root_prc) THEN 
                IF (firsttime_mc) THEN
                   ! Open and read dimenions in file
                   CALL flininfo('nudge_moistc.nc',  iim_file, jjm_file, llm_file, ttm_mc, mc_id)
                   
                   ! Coherence test between dimension in the file and in the model run
                   IF ((iim_file /= iim_g) .OR. (jjm_file /= jjm_g)) THEN
                      WRITE(numout,*) 'hydrol_nudge: iim_file, jjm_file, llm_file, ttm_mc=', &
                           iim_file, jjm_file, llm_file, ttm_mc
                      WRITE(numout,*) 'hydrol_nudge: iim_g, jjm_g=', iim_g, jjm_g
                      CALL ipslerr_p(3,'hydrol_nudge','Problem in coherence between dimensions in nudge_moistc.nc file and model',&
                           'No interpolation is done on this file','This input file must be on the same horizontal resolution as the model.')
                   END IF
                   
                   firsttime_mc=.FALSE.
                   istart_mc=julian_diff-1 ! initialize time counter to read
                   IF (printlev>=2) WRITE(numout,*) "Start read nudge_moistc.nc file at time step: ", istart_mc+1
                END IF

                istart_mc=istart_mc+1  ! read next time step in the file
                iend=istart_mc         ! only read 1 time step
                
                ! Read mc from file, one variable per soiltile
                IF (printlev>=3) WRITE(numout,*) &
                     "Read variables moistc_1, moistc_2 and moistc_3 from nudge_moistc.nc at time step: ", istart_mc
                CALL flinget (mc_id, 'moistc_1', iim_g, jjm_g, nslm, ttm_mc, istart_mc, iend, mc_read_glo2D_1)
                CALL flinget (mc_id, 'moistc_2', iim_g, jjm_g, nslm, ttm_mc, istart_mc, iend, mc_read_glo2D_2)
                CALL flinget (mc_id, 'moistc_3', iim_g, jjm_g, nslm, ttm_mc, istart_mc, iend, mc_read_glo2D_3)

                ! Transform from global 2D(iim_g, jjm_g) into into land-only global 1D(nbp_glo)
                ! Put the variables on the 3 soiltiles in the same file
                DO ji = 1, nbp_glo
                   j = ((index_g(ji)-1)/iim_g) + 1
                   i = (index_g(ji) - (j-1)*iim_g)
                   mc_read_glo1D(ji,:,1) = mc_read_glo2D_1(i,j,:,1)
                   mc_read_glo1D(ji,:,2) = mc_read_glo2D_2(i,j,:,1)
                   mc_read_glo1D(ji,:,3) = mc_read_glo2D_3(i,j,:,1)
                END DO
             END IF

             ! Distribute the fields on all processors
             CALL scatter(mc_read_glo1D, mc_read_next)

             ! No interpolation is done, set the mask to 1
             mask_mc_interp(:,:,:) = 1

          END IF ! nudge_interpol_with_xios
       END IF ! MOD(kjit,INT(one_day/dt_sechiba)) == 1
       
      
       !! 1.3 Linear time interpolation between daily fields to the current time step
       tau   = (kjit-1)*dt_sechiba/one_day - AINT((kjit-1)*dt_sechiba/one_day)
       mc_read_current(:,:,:) = (1.-tau)*mc_read_prev(:,:,:) + tau*mc_read_next(:,:,:)

       !! 1.4 Output daily fields and time interpolated fields only for debugging and validation purpose
       CALL xios_orchidee_send_field("mc_read_next", mc_read_next)
       CALL xios_orchidee_send_field("mc_read_current", mc_read_current)
       CALL xios_orchidee_send_field("mc_read_prev", mc_read_prev)
       CALL xios_orchidee_send_field("mask_mc_interp_out", mask_mc_interp)


  END SUBROUTINE hydrol_nudge_mc_read

!! ================================================================================================================================
!! SUBROUTINE   : hydrol_nudge_mc
!!
!>\BRIEF         Applay nuding for soil moisture
!!
!! DESCRIPTION  : Applay nudging for soil moisture. The nuding values were previously read and interpolated using 
!!                the subroutine hydrol_nudge_mc_read
!!                This subroutine is called from a loop over all soil tiles.
!!
!! RECENT CHANGE(S) : None
!!
!! \n
!_ ================================================================================================================================
  SUBROUTINE hydrol_nudge_mc(kjpindex, jst, mc_loc)

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                                :: kjpindex    !! Domain size
    INTEGER(i_std), INTENT(in)                                :: jst         !! Index for current soil tile
       
    !! 0.2 Modified variables
    REAL(r_std), DIMENSION(kjpindex,nslm,nstm), INTENT(inout) :: mc_loc      !! Soil moisture
    
    !! 0.2 Locals variables
    REAL(r_std), DIMENSION(kjpindex,nslm,nstm) :: mc_aux                     !! Temorary variable for calculation of nudgincsm
    INTEGER(i_std)                             :: ji, jsl                    !! loop index    
    
    
    !! 1.5 Applay nudging of soil moisture using alpha_nudge_mc at each model sechiba time step.
    !!     alpha_mc_nudge calculated using the parameter for relaxation time NUDGE_TAU_MC set in module constantes.
    !!     alpha_nudge_mc is between 0-1
    !!     If alpha_nudge_mc=1, the new mc will be replaced by the one read from file 
    mc_loc(:,:,jst) = (1-alpha_nudge_mc)*mc_loc(:,:,jst) + alpha_nudge_mc * mc_read_current(:,:,jst)
    
    
    !! 1.6 Calculate diagnostic for nudging increment of water in soil moisture
    !!     Here calculate tmc_aux for the current soil tile. Later in hydrol_nudge_mc_diag, this will be used to calculate nudgincsm
    mc_aux(:,:,jst)  = alpha_nudge_mc * ( mc_read_current(:,:,jst) - mc_loc(:,:,jst))
    DO ji=1,kjpindex
       tmc_aux(ji,jst) = dz(2) * ( trois*mc_aux(ji,1,jst) + mc_aux(ji,2,jst) )/huit
       DO jsl = 2,nslm-1
          tmc_aux(ji,jst) = tmc_aux(ji,jst) + dz(jsl) *  (trois*mc_aux(ji,jsl,jst)+mc_aux(ji,jsl-1,jst))/huit &
               + dz(jsl+1) * (trois*mc_aux(ji,jsl,jst)+mc_aux(ji,jsl+1,jst))/huit
       ENDDO
       tmc_aux(ji,jst) = tmc_aux(ji,jst) + dz(nslm) * (trois*mc_aux(ji,nslm,jst) + mc_aux(ji,nslm-1,jst))/huit
    ENDDO
       

  END SUBROUTINE hydrol_nudge_mc


  SUBROUTINE hydrol_nudge_mc_diag(kjpindex, soiltile)
    !! 0.1 Input variables    
    INTEGER(i_std), INTENT(in)                         :: kjpindex    !! Domain size
    REAL(r_std), DIMENSION(kjpindex,nstm), INTENT (in) :: soiltile    !! Fraction of each soil tile within vegtot (0-1, unitless)
    
    !! 0.2 Locals variables
    REAL(r_std), DIMENSION(kjpindex)           :: nudgincsm             !! Nudging increment of water in soil moisture
    INTEGER(i_std)                             :: ji, jst               !! loop index


    ! Average over grid-cell
    nudgincsm(:) = zero
    DO jst=1,nstm
       DO ji=1,kjpindex
          nudgincsm(ji) = nudgincsm(ji) + vegtot(ji) * soiltile(ji,jst) * tmc_aux(ji,jst)
       ENDDO
    ENDDO
    
    CALL xios_orchidee_send_field("nudgincsm", nudgincsm)

  END SUBROUTINE hydrol_nudge_mc_diag


  !! ================================================================================================================================
  !! SUBROUTINE   : hydrol_nudge_snow
  !!
  !>\BRIEF         Read, interpolate and applay nudging snow variables
  !!
  !! DESCRIPTION  : Nudging of snow variables is done if OK_NUDGE_SNOW=y is set in run.def
  !!
  !! RECENT CHANGE(S) : None
  !!
  !! MAIN IN-OUTPUT VARIABLE(S) : snowdz, snowrho, snowtemp
  !!
  !! REFERENCE(S) : 
  !!
  !! \n
  !_ ================================================================================================================================


  SUBROUTINE hydrol_nudge_snow(kjit,   kjpindex, snowdz, snowrho, snowtemp )

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                         :: kjit        !! Timestep number
    INTEGER(i_std), INTENT(in)                         :: kjpindex    !! Domain size

    !! 0.2 Modified variables
    REAL(r_std), DIMENSION(kjpindex,nsnow), INTENT(inout)     :: snowdz      !! Snow layer thickness [m]
    REAL(r_std), DIMENSION(kjpindex,nsnow), INTENT(inout)     :: snowrho     !! Snow density (Kg/m^3)
    REAL(r_std), DIMENSION(kjpindex,nsnow), INTENT(inout)     :: snowtemp    !! Snow temperature (K)



    !! 0.3 Locals variables
    REAL(r_std)                                :: tau                   !! Position between to values in nudge mc file
    REAL(r_std), DIMENSION(kjpindex,nsnow)     :: snowdz_read_current   !! snowdz from file interpolated to current timestep [m]
    REAL(r_std), DIMENSION(kjpindex,nsnow)     :: snowrho_read_current  !! snowrho from file interpolated to current timestep (Kg/m^3)
    REAL(r_std), DIMENSION(kjpindex,nsnow)     :: snowtemp_read_current !! snowtemp from file interpolated to current timestep (K)
    REAL(r_std), DIMENSION(kjpindex)           :: nudgincswe            !! Nudging increment of water in snow
    REAL(r_std), DIMENSION(iim_g,jjm_g,nsnow,1):: snowdz_read_glo2D     !! snowdz from file at global 2D(lat,lon) grid [m]
    REAL(r_std), DIMENSION(iim_g,jjm_g,nsnow,1):: snowrho_read_glo2D    !! snowrho from file at global 2D(lat,lon) grid (Kg/m^3)
    REAL(r_std), DIMENSION(iim_g,jjm_g,nsnow,1):: snowtemp_read_glo2D   !! snowrho from file at global 2D(lat,lon) grid (K)
    REAL(r_std), DIMENSION(nbp_glo,nsnow)      :: snowdz_read_glo1D     !! snowdz_read_glo2D on land-only vector form, in global (m)
    REAL(r_std), DIMENSION(nbp_glo,nsnow)      :: snowrho_read_glo1D    !! snowdz_read_glo2D on land-only vector form, in global (Kg/m^3)
    REAL(r_std), DIMENSION(nbp_glo,nsnow)      :: snowtemp_read_glo1D   !! snowdz_read_glo2D on land-only vector form, in global (K)
    INTEGER(i_std), SAVE                       ::  istart_snow!! start index to read from input file
!$OMP THREADPRIVATE(istart_snow)
    INTEGER(i_std)                             :: iend                  !! end index to read from input file
    INTEGER(i_std)                             :: i, j, ji, jg, jst, jsl!! loop index
    INTEGER(i_std)                             :: iim_file, jjm_file, llm_file !! Dimensions in input file
    INTEGER(i_std), SAVE                       :: ttm_snow      !! Time dimensions in input file
!$OMP THREADPRIVATE(ttm_snow)
    INTEGER(i_std), SAVE                       :: snow_id        !! index for netcdf files
!$OMP THREADPRIVATE(snow_id)
    LOGICAL, SAVE                              :: firsttime_snow=.TRUE.
!$OMP THREADPRIVATE(firsttime_snow)

 
    !! 2. Nudging of snow variables
    IF (ok_nudge_snow) THEN

       !! 2.1 Read snow variables from file, once a day only
       !!     The forcing file must contain daily frequency values for the full year of the simulation
       IF (MOD(kjit,INT(one_day/dt_sechiba)) == 1) THEN 
          ! Save variables from previous day
          snowdz_read_prev   = snowdz_read_next
          snowrho_read_prev  = snowrho_read_next
          snowtemp_read_prev = snowtemp_read_next
          
          IF (nudge_interpol_with_xios) THEN
             ! Read and interpolation snow variables and the mask from input file
             CALL xios_orchidee_recv_field("snowdz_interp", snowdz_read_next)
             CALL xios_orchidee_recv_field("snowrho_interp", snowrho_read_next)
             CALL xios_orchidee_recv_field("snowtemp_interp", snowtemp_read_next)
             CALL xios_orchidee_recv_field("mask_snow_interp", mask_snow_interp)

          ELSE
             ! Only read fields from the file. We here suppose that no interpolation is needed.
             IF (is_root_prc) THEN 
                IF (firsttime_snow) THEN
                   ! Open and read dimenions in file
                   CALL flininfo('nudge_snow.nc',  iim_file, jjm_file, llm_file, ttm_snow, snow_id)
                   
                   ! Coherence test between dimension in the file and in the model run
                   IF ((iim_file /= iim_g) .OR. (jjm_file /= jjm_g)) THEN
                      WRITE(numout,*) 'hydrol_nudge: iim_file, jjm_file, llm_file, ttm_snow=', &
                           iim_file, jjm_file, llm_file, ttm_snow
                      WRITE(numout,*) 'hydrol_nudge: iim_g, jjm_g=', iim_g, jjm_g
                      CALL ipslerr_p(3,'hydrol_nudge','Problem in coherence between dimensions in nudge_snow.nc file and model',&
                           'iim_file should be equal to iim_g','jjm_file should be equal to jjm_g')
                   END IF
                                         
                   firsttime_snow=.FALSE.
                   istart_snow=julian_diff-1  ! initialize time counter to read
                   IF (printlev>=2) WRITE(numout,*) "Start read nudge_snow.nc file at time step: ", istart_snow+1
                END IF

                istart_snow=istart_snow+1  ! read next time step in the file
                iend=istart_snow      ! only read 1 time step
                
                ! Read snowdz, snowrho and snowtemp from file
                IF (printlev>=2) WRITE(numout,*) &
                  "Read variables snowdz, snowrho and snowtemp from nudge_snow.nc at time step: ", istart_snow,ttm_snow
                CALL flinget (snow_id, 'snowdz', iim_g, jjm_g, nsnow, ttm_snow, istart_snow, iend, snowdz_read_glo2D)
                CALL flinget (snow_id, 'snowrho', iim_g, jjm_g, nsnow, ttm_snow, istart_snow, iend, snowrho_read_glo2D)
                CALL flinget (snow_id, 'snowtemp', iim_g, jjm_g, nsnow, ttm_snow, istart_snow, iend, snowtemp_read_glo2D)


                ! Transform from global 2D(iim_g, jjm_g) variables into into land-only global 1D variables (nbp_glo)
                DO ji = 1, nbp_glo
                   j = ((index_g(ji)-1)/iim_g) + 1
                   i = (index_g(ji) - (j-1)*iim_g)
                   snowdz_read_glo1D(ji,:) = snowdz_read_glo2D(i,j,:,1)
                   snowrho_read_glo1D(ji,:) = snowrho_read_glo2D(i,j,:,1)
                   snowtemp_read_glo1D(ji,:) = snowtemp_read_glo2D(i,j,:,1)
                END DO
             END IF

             ! Distribute the fields on all processors
             CALL scatter(snowdz_read_glo1D, snowdz_read_next)
             CALL scatter(snowrho_read_glo1D, snowrho_read_next)
             CALL scatter(snowtemp_read_glo1D, snowtemp_read_next)

             ! No interpolation is done, set the mask to 1
             mask_snow_interp=1

          END IF ! nudge_interpol_with_xios

          
          ! Test if the values for depth of snow is in a valid range when read from the file, 
          ! else set as no snow cover
          DO ji=1,kjpindex
             IF ((SUM(snowdz_read_next(ji,:)) .LE. 0.0) .OR. (SUM(snowdz_read_next(ji,:)) .GT. 100)) THEN
                ! Snowdz has no valide values in the file, set here as no snow
                snowdz_read_next(ji,:)   = 0
                snowrho_read_next(ji,:)  = 50.0
                snowtemp_read_next(ji,:) = tp_00
             END IF
          END DO

       END IF ! MOD(kjit,INT(one_day/dt_sechiba)) == 1
       
      
       !! 2.2 Linear time interpolation between daily fields for current time step
       tau   = (kjit-1)*dt_sechiba/one_day - AINT((kjit-1)*dt_sechiba/one_day)
       snowdz_read_current(:,:) = (1.-tau)*snowdz_read_prev(:,:) + tau*snowdz_read_next(:,:)
       snowrho_read_current(:,:) = (1.-tau)*snowrho_read_prev(:,:) + tau*snowrho_read_next(:,:)
       snowtemp_read_current(:,:) = (1.-tau)*snowtemp_read_prev(:,:) + tau*snowtemp_read_next(:,:)

       !! 2.3 Output daily fields and time interpolated fields only for debugging and validation purpose
       CALL xios_orchidee_send_field("snowdz_read_next", snowdz_read_next)
       CALL xios_orchidee_send_field("snowdz_read_current", snowdz_read_current)
       CALL xios_orchidee_send_field("snowdz_read_prev", snowdz_read_prev)
       CALL xios_orchidee_send_field("snowrho_read_next", snowrho_read_next)
       CALL xios_orchidee_send_field("snowrho_read_current", snowrho_read_current)
       CALL xios_orchidee_send_field("snowrho_read_prev", snowrho_read_prev)
       CALL xios_orchidee_send_field("snowtemp_read_next", snowtemp_read_next)
       CALL xios_orchidee_send_field("snowtemp_read_current", snowtemp_read_current)
       CALL xios_orchidee_send_field("snowtemp_read_prev", snowtemp_read_prev)
       CALL xios_orchidee_send_field("mask_snow_interp_out", mask_snow_interp)

       !! 2.4 Applay nudging of snow variables using alpha_nudge_snow at each model sechiba time step.
       !!     alpha_snow_nudge calculated using the parameter for relaxation time NUDGE_TAU_SNOW set in module constantes.
       !!     alpha_nudge_snow is between 0-1
       !!     If alpha_nudge_snow=1, the new snow variables will be replaced by the ones read from file.
       snowdz(:,:) = (1-alpha_nudge_snow)*snowdz(:,:) + alpha_nudge_snow * snowdz_read_current(:,:)
       snowrho(:,:) = (1-alpha_nudge_snow)*snowrho(:,:) + alpha_nudge_snow * snowrho_read_current(:,:)
       snowtemp(:,:) = (1-alpha_nudge_snow)*snowtemp(:,:) + alpha_nudge_snow * snowtemp_read_current(:,:)

       !! 2.5 Calculate diagnostic for the nudging increment of water in snow
       nudgincswe=0.
       DO jg = 1, nsnow 
          nudgincswe(:) = nudgincswe(:) +  &
               alpha_nudge_snow*(snowdz_read_current(:,jg)*snowrho_read_current(:,jg)-snowdz(:,jg)*snowrho(:,jg))
       END DO
       CALL xios_orchidee_send_field("nudgincswe", nudgincswe)
       
    END IF

  END SUBROUTINE hydrol_nudge_snow


!! ================================================================================================================================
!! SUBROUTINE 	: hydrol_hydraulic_arch_tuzet_calc
!!
!>\BRIEF       Calculates the water status of the vegetation and the water transport inside the plants
!!             accounting for the effects of hydrualic architecture following (Tuzet et al. (2003/2008/2017))
!!
!!\n DESCRIPTION : 
!!   
!!   The module calculate the water potential inside the architecture of the plant. The scheme of resolution is the following:
!!      1) Calculation of the fluxes at each stage of the vegetation; 
!!      2) Calculation of the root water absorption from the soil: two methods here, classical resistance scheme are resolution of 
!!         Richard's equation in cylindrical coordinates in an "absorption muff". The outpur of the method is the root water
!!         potential in each soil layer (see Tuzet et al. (2003))
!!      3) Calculation of the water potential at each stage of the vegetation
!!      
!! RECENT CHANGE(S): Added by Julien Alléon (December 2022)
!!
!! MAIN OUTPUT VARIABLE(S): :: psi_leaf  
!!
!! REFERENCE(S)	: Tuzet et al. 2017
!!                Tuzet et al. 2008
!!                Tuzet et al. 2003
!!                Bonan et al. 2014
!!
!! FLOWCHART    : 
!! 
!! ================================================================================================================================


  SUBROUTINE hydrol_hydraulic_arch_tuzet_calc (kjit, kjpindex, ks, nvan, avan, transpir, mc_out, mcw, root_profile, veget, veget_max,         &
       njsc, soiltile, circ_class_n, circ_class_biomass, u, v, tq_cdrag, gsmean, pb, temp_air, qair, qsol_sat_new, rau,                       & 
       frac_evap_flood, frac_sub_snow, frac_evap_transp, cimean, assimtot, Rdtot_pft, lalo, z_array_out, psi_leaf, psi_leaf_next, psi_leaf_midday,&
       psi_sto_leaf_save, psi_sto_wood_save, e_frac, psi_root_sup, psi_root_inf, psi_soil_sup, psi_soil_inf, psi_xylem_trunk, psi_xylem_leaf, &
       psi_xylem_collar, psi_sto_wood, psi_sto_leaf, mc_i_sup, mc_i_inf, F_absorption, vessel_loss)

    !! Variable declaration

    !! Input Variables

    INTEGER(i_std),                                   INTENT(in)      :: kjit                     !! Time step number (unitless)
    INTEGER(i_std),                                   INTENT(in)      :: kjpindex                 !! Domain size, terrestrial pixels only (unitless)
    REAL(r_std),        DIMENSION(:),                 INTENT(in)      :: ks                       !! Hydraulic conductivity at saturation (mm {-1})
    REAL(r_std),        DIMENSION(:),                 INTENT(in)      :: nvan                     !! Van Genuchten coeficients n (unitless)
    REAL(r_std),        DIMENSION(:),                 INTENT(in)      :: avan                     !! Van Genuchten coeficients a (mm-1})
    REAL(r_std),        DIMENSION(:,:),               INTENT(in)      :: transpir                 !! Unstressed transpiration calculated by the energy module (mm/d)
    REAL(r_std),        DIMENSION(:,:,:),             INTENT(in)      :: mc_out                   !! Soil water content per soil layer per soil tile (m^3/m^3)
    REAL(r_std),        DIMENSION(:),                 INTENT(in)      :: mcw                      !! Volumetric water content at wilting point (m^{3} m^{-3})
    REAL(r_std),        DIMENSION(:,:,:,:),           INTENT(in)      :: root_profile             !! Normalized root mass/length fraction in each soil layer 
    REAL(r_std),        DIMENSION(:,:),               INTENT(in)      :: veget                    !! Fraction of vegetation type (unitless)
    REAL(r_std),        DIMENSION(:,:),               INTENT(in)      :: veget_max                !! Maximum fraction of vegetation type (unitless)
    INTEGER(i_std),     DIMENSION(:),                 INTENT(in)      :: njsc                     !! Index of the dominant soil textural class (unitless)
    REAL(r_std),        DIMENSION(:,:),               INTENT(in)      :: soiltile                 !! Fraction of each soiltile within vegtot (0-1, unitless) 
    REAL(r_std),        DIMENSION(:,:,:),             INTENT(in)      :: circ_class_n             !! Number of trees within a circumference class (tree m-2)
    REAL(r_std),        DIMENSION(:,:,:,:,:),         INTENT(in)      :: circ_class_biomass       !! Biomass components of the model tree  
                                                                                                  !! within a circumference class
                                                                                                  !! class @tex $(g C ind^{-1})$ @endtex  
    REAL(r_std),        DIMENSION (:),                INTENT(in)      :: u,v                      !! Horizontal wind speed (m/s)
    REAL(r_std),        DIMENSION (:),                INTENT(in)      :: tq_cdrag                 !! Surface drag coefficient (-)
    REAL(r_std),        DIMENSION (:,:),              INTENT(in)      :: gsmean                   !! Stomatal conductance (mol m-2 s-1) 
    REAL(r_std),        DIMENSION (:),                INTENT(in)      :: pb                       !! Lowest level atmospheric air pressure (hPa)
    REAL(r_std),        DIMENSION (:),                INTENT(in)      :: temp_air                 !! Lowest level atmospheric air temperature (K)
    REAL(r_std),        DIMENSION (:),                INTENT(in)      :: qair                     !! Lowest level atmospheric air humidity (kg/kg)
    REAL(r_std),        DIMENSION (:),                INTENT(in)      :: qsol_sat_new             !! Surface saturated humidity (kg/kg)
    REAL(r_std),        DIMENSION (:),                INTENT(in)      :: rau                      !! Air density (kg/m3)
    REAL(r_std),        DIMENSION (:),                INTENT(in)      :: frac_evap_flood          !! Fraction of gridcell that evaporates floodplains (-)
    REAL(r_std),        DIMENSION (:),                INTENT(in)      :: frac_sub_snow            !! Fraction of gridcell that sublimates snow (-)
    REAL(r_std),        DIMENSION (:,:),              INTENT(in)      :: frac_evap_transp         !! Fraction of gridcell that transpires (-)
    REAL(r_std),        DIMENSION (:,:),              INTENT(in)      :: cimean                   !! Mean leaf ci (ppm)
    REAL(r_std),        DIMENSION (:,:),              INTENT(in)      :: assimtot                 !! Total assimilation per pft (mu mol m-2 s-1)
    REAL(r_std),        DIMENSION (:,:),              INTENT(in)      :: Rdtot_pft                !! Total respiration per pft (mu mol m-2 s-1)
    REAL(r_std),        DIMENSION (:,:),              INTENT(in)      :: lalo                     !! Latitude and Longitude (°)
    REAL(r_std),        DIMENSION (:,:,:,:),          INTENT(in)      :: z_array_out              !! Heights above soil of the different PGap nodes (m)


    !! Output Variables

    REAL(r_std),        DIMENSION(:,:),               INTENT(out)     :: psi_leaf_next            !! Approximated Leaf Water Potential at time step n+1 (MPa) (= psi_leaf when no stress)
    REAL(r_std),        DIMENSION(:,:,:),             INTENT(out)     :: psi_leaf                 !! Leaf Water Potential (MPa)
    REAL(r_std),        DIMENSION(:,:,:),             INTENT(out)     :: psi_leaf_midday          !! Leaf Water Potential at midday (MPa)

    REAL(r_std),        DIMENSION(:,:,:,:),           INTENT(out)     :: e_frac                   !! Fraction of water transpired supplied by individual layers (no units) 
    REAL(r_std),        DIMENSION(:,:,:),             INTENT(out)     :: vessel_loss              !! Cavitation during timestep


    !! Modified Variables

    ! Water potentials
    REAL(r_std),        DIMENSION(:,:),               INTENT(inout)   :: psi_soil_sup             !! Water potential of the superficial soil level (MPa)
    REAL(r_std),        DIMENSION(:,:),               INTENT(inout)   :: psi_soil_inf             !! Water potential of the inferior soil level (MPa)
    REAL(r_std),        DIMENSION(:,:,:),             INTENT(inout)   :: psi_root_sup             !! Water potential of the roots inside the superficial soil level (MPa)
    REAL(r_std),        DIMENSION(:,:,:),             INTENT(inout)   :: psi_root_inf             !! Water potential of the roots inside the inferior soil level (MPa)
    REAL(r_std),        DIMENSION(:,:,:),             INTENT(inout)   :: psi_xylem_trunk          !! Xylem (trunk level) Water potential (MPa)
    REAL(r_std),        DIMENSION(:,:,:),             INTENT(inout)   :: psi_xylem_leaf           !! Xylem (leaf level) Water potential (MPa)
    REAL(r_std),        DIMENSION(:,:,:),             INTENT(inout)   :: psi_xylem_collar         !! Xylem (collar level) Water potential (MPa)
    REAL(r_std),        DIMENSION(:,:,:),             INTENT(inout)   :: psi_sto_wood             !! Wood storage Water potential (MPa)
    REAL(r_std),        DIMENSION(:,:,:),             INTENT(inout)   :: psi_sto_leaf             !! Leaf Water potential at time step n-1 (MPa)
    REAL(r_std),        DIMENSION(:,:,:),             INTENT(inout)   :: psi_sto_wood_save        !! Wood storage Water potential at time step n-1 (MPa)
    REAL(r_std),        DIMENSION(:,:,:),             INTENT(inout)   :: psi_sto_leaf_save        !! Leaf storage Water potential at time step n-1 (MPa)

    ! Soil water contents
    REAL(r_std),        DIMENSION(:,:,:,:),           INTENT(inout)   :: mc_i_sup                 !! Water content at each node of the absorption muff in the superficial soil level (m^3/m^3)
    REAL(r_std),        DIMENSION(:,:,:,:),           INTENT(inout)   :: mc_i_inf                 !! Water content at each node of the absorption muff in the inferior soil level (m^3/m^3)

    ! Water fluxes
    REAL(r_std),        DIMENSION(:,:),               INTENT(inout)   :: F_absorption             !! Total root absorption flux (m^3/s)



    !! Local Variables

    ! Soil related variables

    REAL(r_std),        DIMENSION(kjpindex,nstm)                      :: mc_sup                   !! Water content in the superficial soil level (m^3/m^3)
    REAL(r_std),        DIMENSION(kjpindex,nstm)                      :: mc_inf                   !! Water potential of the inferior soil level (m^3/m^3)
    REAL(r_std),        DIMENSION(kjpindex,nstm)                      :: mc_sup_save              !! Water content in the superficial soil level at time step n-1 (m^3/m^3)
    REAL(r_std),        DIMENSION(kjpindex,nstm)                      :: mc_inf_save              !! Water potential of the inferior soil level at time step n-1 (m^3/m^3)

    ! Fluxes    
    REAL(r_std),        DIMENSION(kjpindex,nvm,nslm,nstm)             :: Fi                       !! Root absorption flux absorbed from each soil layers (m^3/s)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: Fxylem_up                !! Water flux in the xylem from the trunk stage to the leaf one (m^3/s)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: Fxylem_low               !! Water flux in the xylem from the collar stage to the trunk one (m^3/s)
    ! Below are storage fluxes which when positive, mean water is being supplied to the water storage organs
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: Fsto_leaf                !! Water supply flux to the leaf water storage (m^3/s)
    REAL(r_std),        DIMENSION(kjpindex,nvm,ncirc)                 :: Fsto_leaf_save           !! Saved water supply flux to the leaf water storage (kg/m^2/s)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: Fsto_wood                !! Water supply flux to the sapwood water storage (m^3/s)
    REAL(r_std),        DIMENSION(kjpindex,nvm,ncirc)                 :: Fsto_wood_save           !! Saved water supply flux to the sapwood water storage (kg/m^2/s)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: Fsup                     !! Superficial roots absorption flux (m^3/s)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: Finf                     !! Inferior roots absorption flux (m^3/s)

    ! Storage related variables
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: V_wood                   !! Current amount of water in the sapwood water storage (m^3)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: V_wood_save              !! Amount of water in the sapwood water storage at time step n-1 (m^3)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: V_leaf                   !! Current amount of water in the leaf water storage (m^3)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: V_leaf_save              !! Amount of water in the leaf water storage at time step n-1 (m^3)
    REAL(r_std),        DIMENSION(kjpindex,nvm,ncirc)                 :: V_leaf_out               !! Current normalized amount of water in the leaf water storage for output (m^3)
    REAL(r_std),        DIMENSION(kjpindex,nvm,ncirc)                 :: V_wood_out               !! Current normalized amount of water in the wood water storage for output (m^3)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: C_wood                   !! Capacitance of the sapwood water storage (m^3/MPa)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: C_leaf                   !! Capacitance of the leaf water storage (m^3/MPa)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: tau_wood                 !! Time constant of the sapwood water storage differential equation (s)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: tau_leaf                 !! Time constant of the leaf water storage differential equation (s)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: Vr_wood                  !! Residual amount of water in the sapwood water storage (m^3/s)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: Vr_leaf                  !! Residual amount of water in the leaf water storage (m^3/s)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: Vmax_wood                !! Maximum amount of water in the sapwood water storage (m^3)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: Vmax_leaf                !! Maximum amount of water in the leaf water storage (m^3)

    ! Water resistances
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: Res_sto_leaf             !! Resistance between leaf storage and xylem at the leaf level (MPa.s/m^3)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: Res_sto_wood             !! Resistance between sapwood storage and xylem at the trunk level (MPa.s/m^3)
    REAL(r_std),        DIMENSION(kjpindex,nvm,ncirc)                 :: Res_xylem_up             !! Xylem resistance between the trunk level and the leaf one (MPa.s/m^3)
    REAL(r_std),        DIMENSION(kjpindex,nvm,ncirc)                 :: Res_xylem_low            !! Xylem resistance between the collar level and the trunk one (MPa.s/m^3)
    REAL(r_std),        DIMENSION(kjpindex,nvm,ncirc)                 :: Res_mesophyll            !! Resistance between the xylem at the leaf level and the stomatal cavities (MPa.s/m^3)
    REAL(r_std),        DIMENSION(kjpindex,nvm,ncirc)                 :: Res_root_sup             !! Resistance between the superficial soil/root interface and the xylem at the collar level (MPa.s/m^3)
    REAL(r_std),        DIMENSION(kjpindex,nvm,ncirc)                 :: Res_root_inf             !! Resistance between the superficial soil/root interface and the xylem at the collar level (MPa.s/m^3)
    REAL(r_std),        DIMENSION(kjpindex,nvm,ncirc)                 :: R_root                   !! Equivalent resistance for the root system (MPa.s/m^3)
    REAL(r_std),        DIMENSION(kjpindex,nvm,ncirc)                 :: K_plant                  !! Whole plant conductance (kg/MPa/s)
    REAL(r_std),        DIMENSION(kjpindex,nvm,ncirc)                 :: delta_psi                !! Soil to leaf potential drop (MPa)
    REAL(r_std),        DIMENSION(kjpindex,nvm,ncirc)                 :: k_xylem_up_max           !! Maximum xylem up conductance  (m^3/s/MPa)
    REAL(r_std),        DIMENSION(kjpindex,nvm,ncirc)                 :: k_xylem_low_max          !! Maximum xylem low conductance  (m^3/s/MPa)

    ! Estimations/predictions related variables

    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: psi_sto_wood_predict     !! Prediction of the wood storage water potential at time step n+1 (MPa)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: psi_sto_leaf_predict     !! Prediction of the leaf storage water potential at time step n+1 (MPa)

    ! Temporary variables (/IF\ (no water stress or current tstep) : = normal variables /ELSE\ = next tstep values)

    REAL(r_std),        DIMENSION(kjpindex,nstm)                      :: psi_soil_sup_temp        !! Superficial soil layer Water Potential (MPa)
    REAL(r_std),        DIMENSION(kjpindex,nstm)                      :: psi_soil_inf_temp        !! Inferior soil layer Water Potential (MPa)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: psi_xylem_leaf_temp      !! Temporary value of the xylem (leaf level) water potential (MPa)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: psi_xylem_trunk_temp     !! Temporary value of the xylem (trunk level) water potential (MPa)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: psi_xylem_collar_temp    !! Temporary value of the xylem (collar level) water potential (MPa)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: psi_root_sup_temp        !! Temporary value of the superficial root water potential (MPa)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: psi_root_inf_temp        !! Temporary value of the inferior root water potential (MPa)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: psi_leaf_temp            !! Temporary value of leaf water potential (MPa)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: psi_leaf_temp_save       !! Saved value of the temporary value of the leaf water potential (MPa)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: psi_sto_wood_temp        !! Temporary value of the wood storage water potential (MPa)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: psi_sto_leaf_temp        !! Temporary value of the leaf storage water potential (MPa)
    REAL(r_std),        DIMENSION(kjpindex,nvm,ncirc)                 :: psi_leaf_previous        !! Leaf water potential at time step n-1 (MPa)
    REAL(r_std),        DIMENSION(kjpindex,nvm,ncirc)                 :: psi_leaf_mean_gs         !! Seasonal mean of psi_leaf (MPa)
    REAL(r_std)                                                       :: grav_grad                !! gravity gradient difference in hydraulic potential

    REAL(r_std),        DIMENSION(kjpindex,nvm,nslm,nstm)             :: Fi_temp                  !! Temporary value of the Root absorption flux absorbed from each soil layers (m^3/s)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: Fsup_temp                !! Temporary value of the superficial roots absorption flux (m^3/s)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: Finf_temp                !! Temporary value of the inferior roots absorption flux (m^3/s)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: F_absorption_temp        !! Temporary value of the total root absorption flux (m^3/s)

    REAL(r_std),        DIMENSION(kjpindex,nslm,nstm)                 :: mc_out_temp              !! Temporary value of the evolution of mc_out when the prediction of the next time step values is launched(m^3/m^-3)
    REAL(r_std),        DIMENSION(kjpindex,nstm)                      :: mc_sup_temp              !! Temporary value of the water content in the superficial soil level (m^3/m^3)
    REAL(r_std),        DIMENSION(kjpindex,nstm)                      :: mc_inf_temp              !! Temporary value of the water content in the inferior soil level (m^3/m^3)
    REAL(r_std),        DIMENSION(kjpindex,nvm,nrp)                   :: mc_i_sup_temp            !! Temporary value of the water content at each node of the absorption muff in the superficial soil level (m^3/m^3)
    REAL(r_std),        DIMENSION(kjpindex,nvm,nrp)                   :: mc_i_inf_temp            !! Temporary value of the water content at each node of the absorption muff in the inferior soil level (m^3/m^3)


    ! Variables used in the predictor/corrector scheme

    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: func_predict             !! Function used for the Predictor/Corrector scheme of Adams-Moulton 3 (water storages)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: func_n                   !! Function used for the Predictor/Corrector scheme of Adams-Moulton 3 (water storages)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: func_n_1                 !! Function used for the Predictor/Corrector scheme of Adams-Moulton 3 (water storages)


    ! Variables used when the simili-implicit scheme is launched

    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: f_psi_1                  !! Sensitivity of the stomatal conductance to psi_leaf at the previous time step (-)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: f_psi_2                  !! Sensitivity of the stomatal conductance to psi_leaf at the next time step (-)
    REAL(r_std)                                                       :: stress_log               !! Term inside a LOG used when reverting a stress equation (-)

    REAL(r_std),        DIMENSION(kjpindex)                           :: speed                    !! Wind speed (m/s)
    INTEGER(i_std),     DIMENSION(kjpindex,nvm,ncirc)                 :: nsub_step                !! Number of iterations in the simili-implicit scheme
    LOGICAL,            DIMENSION(kjpindex,nvm,ncirc)                 :: launch_next_calc         !! Flag to launch the approximation of the next psi_leaf with the simili-implicit scheme
    LOGICAL,            DIMENSION(kjpindex,nvm,ncirc)                 :: this_tstep_calc          !! Flag to differ this time step calculation and the future approximation
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: transpir_test            !! Estimation of the new time-step transpiration (simili-implicit scheme)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: g_stom                   !! Stomatal conductance at the previous time step (s/m)  
    REAL(r_std),        DIMENSION(kjpindex)                           :: delta                    !! Ratio of sine between two time steps (used in the simili-impliit scheme)
    REAL(r_std),        DIMENSION(kjpindex)                           :: h0                       !! Difference in ° between noon and sunrise
    REAL(r_std),        DIMENSION(kjpindex)                           :: sin_t1                   !! Sine value at the previous time step (-)
    REAL(r_std),        DIMENSION(kjpindex)                           :: sin_t2                   !! Sine value at the new time step (-)  
    INTEGER(i_std)                                                    :: kjit_temp                !! This time-step in current day
    INTEGER(i_std),     DIMENSION(kjpindex)                           :: kjit_sunrise             !! Time-step for sunrise
    INTEGER(i_std),     DIMENSION(kjpindex)                           :: kjit_delta               !! Number of time-step between sunrise and sunset
    REAL(r_std)                                                       :: gamma                    !! Solar Day angle   
    REAL(r_std)                                                       :: dec                      !! Sun declination

    ! Variables used for the transitions between PFT and circumference classes
    REAL(r_std),        DIMENSION(ncirc)                              :: pft_to_indiv             !! Factor to convert PFT to circ_class
    REAL(r_std),        DIMENSION(ncirc)                              :: circ_class_to_indiv      !! Factor to convert circ class to individual tree
    REAL(r_std)                                                       :: dia_median               !! Median of the diameters of the classes used in Deleuze split equation
    REAL(r_std)                                                       :: dia_mean                 !! Mean of the diameters of the classes used in Deleuze split equation
    REAL(r_std)                                                       :: sigma                    !! Sigma constant in Deleuze split equation
    REAL(r_std)                                                       :: deleuze_power            !! Param in Deleuze split equation
    REAL(r_std),        DIMENSION(ncirc)                              :: deleuze_func             !! Deleuze split equation
    REAL(r_std),        DIMENSION(ncirc)                              :: circ_dia                 !! Tree diameters of each circumference class (m)
    REAL(r_std),        DIMENSION(ncirc)                              :: circ_height              !! Tree height of each circumference class (m)
    REAL(r_std),        DIMENSION(ncirc)                              :: crown_volume             !! Crown volume of each circumference class (m^3)
    REAL(r_std),        DIMENSION(ncirc)                              :: crown_width              !! Crown width of each circumference class (m)
    REAL(r_std),        DIMENSION(ncirc)                              :: transpir_circ            !! Transpiration of each circumference class after Deleuze split
    REAL(r_std),        DIMENSION(ncirc)                              :: gs_circ                  !! Stomatal conductance calculated for each circumference class after Deleuze split
    REAL(r_std),        DIMENSION(ncirc)                              :: assim_circ               !! Assimilation of each circumference class after Deleuze split
    REAL(r_std),        DIMENSION(ncirc)                              :: resp_circ                !! Respiration of each circumference class after Deleuze split
    REAL(r_std),        DIMENSION(ncirc)                              :: lai_circ                 !! LAI of each circumference class after Deleuze split
    REAL(r_std),        DIMENSION(kjpindex)                           :: gamma_star               !! CO2 compensation point (ppm)
    REAL(r_std),        DIMENSION(kjpindex,nvm,ncirc)                 :: dummy_transpir           !! Transpiration of each circ class at end of the calculation

    ! Indices
    INTEGER(i_std)                                                    :: ipts,ivm,isl             !! Indices (respectively: grid-cells, PFTs, soil layers)
    INTEGER(i_std)                                                    :: ist,jst,i,icir          !! Indices (respectively: soil textures, "pref_soil_veg", implicit scheme iterations, circumference classes)
    ! Counter
    INTEGER(i_std)                                                    :: counter                  !! Counts if dynamic resistances are curtailed



!! =================================================================================================================================================
    

    !! The hydraulic architecture calculation differs according to the impact of the water stress. 

    !! Calculation setup:

    !! In most of the cases, the calculation will follow the following process:
    !!   0) Hydraulic resistances set up;
    !!   1) Determination of the soil water content in both layers;
    !!   2) Calculation of the water fluxes inside the plant;
    !!   3) Call to the root absorption module;
    !!   4) Link with the hydrol.f90 module;
    !!   5) Calculation of the water potential (especially the leaf water potential which will be used in diffuco.f90).
    
    !! In period of drought, the strong sensibility of stomatal conductance to leaf water potential engenders the need
    !! to solve the vegetation hydraulic budget implicitly. If the budget is solved explicitly, the leaf water potential
    !! and the transpiration flux will oscillate. However, implementing a fully implicit scheme is impossible in our case
    !! as this would mean iterate over at least 3 modules (diffuco.f90 to calculate the stomatal conductance and vbeta_3,
    !! enerbil.f90 to calculate the transpiration flux and hydrol.f90 to calculate the hydraulic architecture). A simili-
    !! implicit scheme thus needs to be implemented. To do so, it has been decided to iterate only on the hydraulic 
    !! architecture routines. The 2 first modules which permits to calculates transpiration are replaced by an "estimation"
    !! of the transpiration flux based on a calculation presented inside the technical note. 

    !! Transpiration basically varies because of photosynthesis, air humidity and leaf water potential. Estimating
    !! transpiration at a next time step can be done by estimatng that photosynthesis and air humidity (mainly driven by
    !! the energy captured by vegetation) are following a kind of sine curve during the day. Thus, using the values of the
    !! previous time step for photosynthesis (via gsmean) and air humidity (via transpir) and multiplying them by a ratio of 
    !! sines between two time steps permit to have a good approximation of the transpiration. This good estimation is then
    !! used in an iterative process that permits to converge toward an equilibrium value of psi_leaf for the next time
    !! step that will be sent to diffuco.f90 in order to compute the stomatal conductance.

    !! This iterative process is not launched at all time steps because in most of the cases, the explicit scheme is enough.
    !! Consequently, according to the value of the leaf water potential at the present time step, the estimation is launched
    !! and flagged by the two booleans "launch_next_calc" and "this_tstep_calc".

    !! When the iterative process is needed, the number of iterations is set to 25 (at maximum, in reality 4 to 5 iterations
    !! are needed). The first iteration (this_tstep_calc = TRUE) calculates the values of the present time step. The 
    !! following iterations are estimating the next value of the leaf water potential. The main input variable (transpir) is
    !! estimated following the process described above, a little soil water budget is calculated and the hydraulic 
    !! architecture is launched at each time step. At the end of each iteration a leaf water potential is calculated. The 
    !! mean value between the previous and current iterations is sent to the model in order to accelerate the convergence.

    !! The function is refined to circ class level to account for inter circ_class variability (height, circumference, biomass)
    !! which all influence the hydraulic behavior. It is being even refined to individual tree level 
    !! (which is not more demanding that circ class in terms of computer time, as we still run once for each circ_class). So it is a nice (because "cheap") 
    !! opportunity to check how well we can capture emerging phenomena (like hydraulic behaviors including total plant resistance)

    !! Set up of the iterative process.

    !! Initialisation:
    !! At the beginning, the number of iterations is set to 1, only one iteration that permits to calculate the values of
    !! this time step.

    nsub_step(:,:,:) = 1                !! Only this time step calculation
    launch_next_calc(:,:,:) = .FALSE.   !! The calculation of the next time step value is not launched
    this_tstep_calc(:,:,:) = .TRUE.     !! Only the first (and only) iteration is launched
    delta(:) = 1.

    ! Debug
    IF (ANY(transpir(:,:)<zero)) THEN
       WRITE(numout,*) 'WARNING : transpir entering tuzet_calc sometimes negative :', transpir(:,:) 
       DO ipts = 1, kjpindex
          DO ivm = 1, nvm
             IF (transpir(ipts,ivm) .LT. zero) THEN
                WRITE(numout,*) 'transpir negative for ipts =', ipts 
                WRITE(numout,*) 'ivm =', ivm 
                EXIT
             ENDIF  
          ENDDO
       ENDDO
    ENDIF 
    !-

    !! If the value of psi_leaf is below a threshold, the iterative process
    !! is launched.
    
    !! DEBUG:
    !! Sometimes, at the beginning of the time step, psi_leaf is not below this threshold whereas at the end, it is.
    !! This leads to phase shift: iterations are launched too late and the instabilities are arriving. Maybe another
    !! condition should be added at the end of the code to counter this problem.
    !! This problem appears during 3 days in 2003, I consider it minor for now.
    
    DO ipts = 1, kjpindex
       
       DO ivm = 1, nvm  

          DO icir = 1,ncirc

             IF (transpir(ipts,ivm) .GT. min_sechiba) THEN

                IF ((psi_leaf(ipts,ivm,icir) .LT. (psi_ref_g0(ivm)-(psi_ref_g0(ivm)/2.))) .AND. (gsmean(ipts,ivm) .NE. 0)) THEN

                   nsub_step(ipts,ivm,icir) = n_iter_tuzet
                   launch_next_calc(ipts,ivm,icir) = .TRUE.
                   psi_leaf_previous(ipts,ivm,icir) = psi_leaf(ipts,ivm,icir)
                   ! ELSE IF (printlev_loc >= 4) THEN
                   !    WRITE(numout,*) 'no iter does happen!!; psi_leaf ='
                   !    WRITE(numout,*) psi_leaf(ipts,ivm)
                   !    CALL flush(numout)
                ENDIF
             ENDIF
          ENDDO
       ENDDO 
       ! Bare soil must show a psi_leaf of 0 
       psi_leaf_previous(ipts,1,:) = xios_default_val
    ENDDO
    CALL xios_orchidee_send_field("transpir_unstressed",transpir)


    dummy_transpir(:,:,:) = zero 
    Finf(:,:) = zero
    Fsup(:,:) = zero
    !! Start of the loop over the grid cells and the PFTs.

    DO ipts = 1, kjpindex
 
       DO ivm = 1, nvm 

          ! Trick to separate the PFT level transpiration/assimilation and respiration into circumference class level variables 
          ! (using Deleuze et al. (2004) formulation) similarly as in allocation module
          
          ! ++++ CHECK ++++++
          ! This function also used in allocation uses a deleuze_b parametrized as 0
          ! thus deleuze_power is always 2. While this parametrization continues
          ! the calculation of deleuze_power and dia_median are useless
          
          IF (ivm .NE. 1) THEN

             IF (is_tree(ivm)) THEN 

                circ_dia(:) = wood_to_dia(circ_class_biomass(ipts,ivm,:,:,icarbon),ivm,pipe_tune2(ipts,ivm))
                circ_height(:) = wood_to_height(circ_class_biomass(ipts,ivm,:,:,icarbon),ivm,pipe_tune2(ipts,ivm))
                crown_volume(:) = wood_to_cv(circ_class_biomass(ipts,ivm,:,:,icarbon),circ_class_n(ipts,ivm,:),ivm,pipe_tune2(ipts,ivm))
                crown_width(:) = ((crown_volume(:)*6)/pi)**(1.0/3.0)

                IF (MOD(ncirc,2) == 0) THEN
                   dia_median = (circ_dia(ncirc/2+1)+circ_dia(ncirc/2))/2.0
                ELSE
                   dia_median = circ_dia(ncirc/2+1)
                ENDIF
                sigma = deleuze_a(ivm) + deleuze_b(ivm) * dia_median
                
                ! +++ CHECK +++
                ! Should not even need this check. How to not enter if no trees ? 
                dia_mean=0
                IF (SUM(circ_class_n(ipts,ivm,:)) .GT. min_sechiba) THEN 
                   DO icir = 1, ncirc
                      dia_mean = dia_mean + circ_dia(icir)*circ_class_n(ipts,ivm,icir)
                   ENDDO
                   dia_mean = dia_mean/SUM(circ_class_n(ipts,ivm,:))
                ENDIF 

                ! everything hard-coded can be found in Deleuze & Dhoté 2004
                deleuze_power = MAX(MIN(1.8 + deleuze_power_a(ivm)*dia_mean,3.5),2.0) 
                 
                deleuze_func(:)=circ_dia(:)-m_dv(ivm)*sigma+((m_dv(ivm)*sigma+circ_dia(:))**2.0-4.0*circ_dia(:)*sigma)**(1/deleuze_power)
                
                ! +++ CHECK ++++++
                ! stomate variables should have also veget_max as scaling factor ? 
                ! This conversion factor unit is m^2(grid cell)/tree
                IF (SUM(circ_class_n(ipts,ivm,:)) .GT. min_sechiba) THEN
		   pft_to_indiv(:) = deleuze_func(:)/SUM(deleuze_func(:))/SUM(circ_class_n(ipts,ivm,:)) 
		ELSE
		   pft_to_indiv(:) = zero
		ENDIF
                ! ++++++++++++++++++++

             ELSE ! grass PFTs

	        !! If there is no diameter, we consider similar fluxes for each circumference class
                pft_to_indiv(:) = 1/ncirc

             ENDIF

             IF (veget_max(ipts,ivm).GT.min_sechiba) THEN
		transpir_circ(:) = transpir(ipts,ivm)*pft_to_indiv(:)/veget_max(ipts,ivm) ! now in /m2(PFT)
	     ELSE
		transpir_circ(:) = zero
	     ENDIF
             assim_circ(:)    = assimtot(ipts,ivm)*pft_to_indiv(:)      ! already in m2PFT
             resp_circ(:)     = Rdtot_pft(ipts,ivm)*pft_to_indiv(:)     ! already in m2PFT
             lai_circ(:)      = circ_class_biomass(ipts,ivm,:,ileaf,icarbon)*sla(ivm)*circ_class_n(ipts,ivm,:)

             IF (.NOT. ok_tuzet_indiv) THEN
		! In that case, it is what each circ. class transpires per m2(grid cell)
		! Otherwise it is the fluxes in each individual tree belonging to each circ. class
                transpir_circ(:) = transpir_circ(:) * circ_class_n(ipts,ivm,:)
                assim_circ(:)    = assim_circ(:) * circ_class_n(ipts,ivm,:)
                resp_circ(:)     = resp_circ(:) * circ_class_n(ipts,ivm,:)
                lai_circ(:)      = lai_circ(:) * circ_class_n(ipts,ivm,:) 
             ENDIF 

          ELSE ! bare soil
             transpir_circ(:) = zero
             assim_circ(:)    = zero
             resp_circ(:)     = zero
             lai_circ(:)      = zero
          ENDIF

          gamma_star(:) = Arrhenius(kjpindex,temp_air,298.,E_gamma_star(ivm))*gamma_star25(ivm)

          IF (ivm==test_pft .AND. ipts==test_grid .AND. printlev_loc >= 3) THEN
             WRITE(numout,*) 'circ_dia =', circ_dia(:) 
             WRITE(numout,*) 'dia_median =', dia_median 
             WRITE(numout,*) 'sigma =', sigma 
             WRITE(numout,*) 'circ_height(:) =', circ_height(:) 
             WRITE(numout,*) 'crown_volume(:) =', crown_volume(:) 
             WRITE(numout,*) 'crown_width(:) =', crown_width(:) 
             WRITE(numout,*) 'deleuze_func =', deleuze_func(:) 
             WRITE(numout,*) 'transpir_cc =', transpir_circ(:) 
             !IF (.NOT. ANY(circ_class_n(ipts,ivm,:).LT.min_sechiba)) THEN
             !   WRITE(numout,*) 'transpir_indiv =',transpir_circ(:)*veget_max(ipts,ivm)/circ_class_n(ipts,ivm,:)
             !ENDIF
          ENDIF
                
          ! Debug         
          IF (printlev_loc.GE.3) THEN
             WRITE(numout,*) 'transpir =', transpir(ipts,ivm) 
             WRITE(numout,*) 'transpir_cc =', transpir_circ(3) 
          ENDIF
          !-

          Fi(ipts,ivm,:,:) = zero   
          F_absorption(ipts,ivm) = zero   

          DO icir = 1, ncirc
         

             !Fi(ipts,ivm,:,:) = zero   
             Fi_temp(ipts,ivm,:,:) = zero 
             !F_absorption(ipts,ivm) = zero  
             mc_i_sup_temp(ipts,ivm,:) = mc_i_sup(ipts,ivm,icir,:)
             mc_i_inf_temp(ipts,ivm,:) = mc_i_inf(ipts,ivm,icir,:)
             mc_out_temp(ipts,:,:) = mc_out(ipts,:,:)
             Fsto_leaf_save(ipts,ivm,icir) = zero
             Fsto_wood_save(ipts,ivm,icir) = zero
             V_leaf_out(ipts,ivm,icir) = zero
             V_wood_out(ipts,ivm,icir) = zero
             vessel_loss(ipts,ivm,icir) = zero
             K_plant(ipts,ivm,icir) = zero

             jst=pref_soil_veg(ivm)

             IF ((veget(ipts,ivm).LT.min_sechiba) .OR. (ivm==1)) THEN
                ! initialization when cycling when pft is bare soil
                ! or absent
                ! ++++ CHECK +++++++++
                ! Will this interfere with the LCC ?
                ! ++++++++++++++++++++ 

                ! Temporary var.
                psi_root_sup_temp(ipts,ivm) = psi_root_sup(ipts,ivm,icir)
                psi_root_inf_temp(ipts,ivm) = psi_root_inf(ipts,ivm,icir)   
                psi_soil_sup_temp(ipts,jst) = psi_soil_sup(ipts,jst)
                psi_soil_inf_temp(ipts,jst) = psi_soil_inf(ipts,jst)
                psi_xylem_leaf_temp(ipts,ivm) = psi_xylem_leaf(ipts,ivm,icir)
                psi_xylem_trunk_temp(ipts,ivm) = psi_xylem_trunk(ipts,ivm,icir)
                psi_xylem_collar_temp(ipts,ivm) = psi_xylem_collar(ipts,ivm,icir)
                psi_sto_leaf_temp(ipts,ivm) = psi_sto_leaf(ipts,ivm,icir)
                psi_sto_wood_temp(ipts,ivm) = psi_sto_wood(ipts,ivm,icir)
                psi_leaf_temp(ipts,ivm) = psi_leaf(ipts,ivm,icir)
                ! The others         
                psi_leaf_mean_gs(ipts,ivm,icir) = zero
                Fxylem_up(ipts,ivm) = zero
                Fxylem_low(ipts,ivm) = zero
                Finf(ipts,ivm) = zero
                Fsup(ipts,ivm) = zero
                mc_sup(ipts,jst) = SUM(mc_out(ipts,:,jst))/nslm ! average (no water uptake)
                mc_inf(ipts,jst) = SUM(mc_out(ipts,:,jst))/nslm
                mc_sup_temp(ipts,jst) = mc_sup(ipts,jst) 
                mc_inf_temp(ipts,jst) = mc_inf(ipts,jst)
                f_psi_1(ipts,ivm) = un ! more consistent than 0, but unused anyway
                f_psi_2(ipts,ivm) = un
                transpir_test(ipts,ivm) = zero
                g_stom(ipts,ivm) = zero
                Res_sto_leaf(ipts,ivm)= zero
                Res_sto_wood(ipts,ivm)= zero
                Res_mesophyll(ipts,ivm,icir)= zero
                Res_xylem_up(ipts,ivm,icir)= zero
                Res_xylem_low(ipts,ivm,icir)= zero
                Res_root_sup(ipts,ivm,icir)= zero
                Res_root_inf(ipts,ivm,icir)= zero
                R_root(ipts,ivm,icir)= zero
                Fsup_temp(ipts,ivm) = zero
                CYCLE 
             ENDIF

             !! Loop over the iterations
             !! In most of the cases, the number of iteration is 1 and all the temporary variables are equal to the real ones.

             next_calc_loop: DO i = 1, nsub_step(ipts,ivm,icir) 

                !! Few initialisations
                !! If the number of iteration is more than 1, the values "mc_i_sup/inf", "psi_root_sup/inf" can differ.

                Fxylem_up(ipts,ivm) = 0.                
                Fxylem_low(ipts,ivm) = 0.        
                mc_sup_temp(ipts,:) = 0.
                mc_inf_temp(ipts,:) = 0.    
                mc_i_sup_temp(ipts,ivm,:) = mc_i_sup(ipts,ivm,icir,:)
                mc_i_inf_temp(ipts,ivm,:) = mc_i_inf(ipts,ivm,icir,:)   
                mc_out_temp(ipts,:,:) = mc_out(ipts,:,:)


                IF(this_tstep_calc(ipts,ivm,icir)) THEN
                   psi_xylem_leaf_temp(ipts,ivm) = psi_xylem_leaf(ipts,ivm,icir)
                   psi_xylem_trunk_temp(ipts,ivm) = psi_xylem_trunk(ipts,ivm,icir)
                   psi_xylem_collar_temp(ipts,ivm) = psi_xylem_collar(ipts,ivm,icir)
                   psi_sto_leaf_temp(ipts,ivm) = psi_sto_leaf(ipts,ivm,icir)
                   psi_sto_wood_temp(ipts,ivm) = psi_sto_wood(ipts,ivm,icir)
                   psi_leaf_temp(ipts,ivm) = psi_leaf(ipts,ivm,icir) 
                   psi_root_sup_temp(ipts,ivm) = psi_root_sup(ipts,ivm,icir)
                   psi_root_inf_temp(ipts,ivm) = psi_root_inf(ipts,ivm,icir)   
                ENDIF

                !! If this is the first iteration, the transpiration corresponds to the one calculated at this time step 
                !! (converted into m^3/s)

                IF (this_tstep_calc(ipts,ivm,icir)) THEN
                   IF ((soiltile(ipts,pref_soil_veg(ivm)) * vegtot(ipts)) .GT. min_sechiba) THEN
                      transpir_test(ipts,ivm) = transpir_circ(icir)/&
                           (dt_sechiba * kilo_to_unit * soiltile(ipts,pref_soil_veg(ivm)) * vegtot(ipts))
                      !! BUG QUAND LA TRANSPIRATION EST NULLE CAR ON CREE UNE TRANSPIRATION EGALE A MIN_SECHIBA²
                   ELSE 

                      transpir_test(ipts,ivm) = zero!min_sechiba*min_sechiba ! min_sec was too big 

                   ENDIF

                   !! If this is not the first iteration, the estimation of the next value of transpiration is launched.

                ELSE

                   IF ((soiltile(ipts,pref_soil_veg(ivm)) * vegtot(ipts)) .GT. min_sechiba) THEN

                      !! The estimation mainly relies on the estimation of the sine values over two time steps.
                      !! The sine is center over the zenith hour at the time step. The sine should be equal to 
                      !! 0 at sunrise and sunset.

                      !! The next calculation focus on the calculation of sunrise and sunset. This relies on the
                      !! declination angle of the sun. This declination angle follows the calculation in src_global/solar.f90.

                      !! Calculation of the declination angle:

                      gamma = 2. * pi * julian_diff/one_year
                      dec =  ( 0.006918-0.399912*COS(gamma)+0.070257*SIN(gamma) &
                           &       -0.006758*COS(2*gamma)+0.000907*SIN(2*gamma)      &
                           &       -0.002697*COS(3*gamma)+0.00148*SIN(3*gamma))


                      !! Calculation of the delta (in °) between zenith and sunrise (converted into hours) and then into 
                      !! the number of the time step at sunrise (kjit_sunrise).
                      ! ++++ CHECK ++++ 
                      ! could be a problem in summer (mistmatch compared to formula
                      ! above)
!!$                   h0(ipts)= ACOS(-TAN(lalo(ipts,1)*pi/180)*TAN(dec))*180/pi
!!$
!!$                   kjit_sunrise(ipts) = INT((24+INT(lalo(ipts,2)/30))-INT(h0(ipts)/30))

                      !! However, this calculation leads to problems (sunrise is most of the time equal to 8:30 am which
                      !! is not a good approximation in summer in Europe). Consequently, for now, the sunrise is set to 
                      !! 6:30 am (kjit_sunrise = 13) and sunset is set to 9 pm (kjit= 42).
                      kjit_sunrise(ipts)=INT(6.5*3600/dt_sechiba)
                      kjit_delta(ipts)=INT(21*3600/dt_sechiba) - kjit_sunrise(ipts)  ! sunset-sunrise

                      kjit_temp = MOD(kjit,INT(one_day/dt_sechiba))

                      !! If the current time step is between sunrise and sunset, the ratio of the sines is calculated.

                      IF ((kjit_temp - kjit_sunrise(ipts).LT. kjit_delta(ipts)) .AND. &!(INT(2*h0(ipts)/30))) .AND.&
                           (kjit_temp+1 - kjit_sunrise(ipts).LT.  kjit_delta(ipts)) .AND. &!(INT(2*h0(ipts)/30))) .AND.&
                           (kjit_temp .GT. kjit_sunrise(ipts)) .AND.&
                           (kjit_temp+1 .GT. kjit_sunrise(ipts))) THEN

                         sin_t1(ipts) = ABS(SIN((kjit_temp - kjit_sunrise(ipts))*pi/kjit_delta(ipts)))!(2*INT(h0(ipts)/30))))
                         sin_t2(ipts) = ABS(SIN((kjit_temp+1 - kjit_sunrise(ipts))*pi/kjit_delta(ipts)))!(2*INT(h0(ipts)/30))))

                         delta(ipts) = sin_t2(ipts)/sin_t1(ipts)
                      ELSE
                         delta(ipts) = 1
                      ENDIF


                      !! The calculation of the next value of transpiration is started.


                      !! The function f_psi at the previous and next time step is computed.
                      IF (printlev_loc>=3) THEN
                         WRITE(numout,*) 'psi_leaf_prev', psi_leaf_previous(ipts,ivm,icir)
                         WRITE(numout,*) 'psi_leaf_temp', psi_leaf_temp(ipts,ivm)
                      ENDIF
                      
                      ! Minimum stress set
                      f_psi_1(ipts,ivm) = MAX(min_sechiba,(1. + EXP(sf(ivm) * psi_ref_g0(ivm))) /&
                                          (1. + EXP(sf(ivm) * (psi_ref_g0(ivm) - psi_leaf_previous(ipts,ivm,icir)))))
                      f_psi_2(ipts,ivm) = MAX(min_sechiba,(1. + EXP(sf(ivm) * psi_ref_g0(ivm))) /&
                                          (1. + EXP(sf(ivm) * (psi_ref_g0(ivm) - psi_leaf_temp(ipts,ivm)))))

                      !! Stomatal conductance is converted into m/s 

                      IF ((cimean(ipts,ivm) - gamma_star(ipts)) .GT. min_sechiba) THEN
                         gs_circ(icir)= MAX(g0(ivm)*lai_circ(icir),(assim_circ(icir)+resp_circ(icir)) /(cimean(ipts,ivm) - gamma_star(ipts)))*f_psi_1(ipts,ivm)

                      ELSE
                         gs_circ(icir)= g0(ivm)*lai_circ(icir)
                      ENDIF
                        
                      ! Debug
                      IF (printlev_loc.GE.3) THEN
                         WRITE(numout,*) 'lai_circ(icir)', lai_circ(icir) 
                         WRITE(numout,*) 'gs_circ(icir)', gs_circ(icir) 
                         WRITE(numout,*) 'g0(ivm)*lai_circ(icir)',g0(ivm)*lai_circ(icir)
                      ENDIF
                      !-

                      ! Putting a minimum to stomatal conductance 
                      gs_circ(icir) = MAX(min_sechiba,gs_circ(icir))
                      

                      g_stom(ipts,ivm) = mol_to_m_1 *(temp_air(ipts)/tp_00)*&
                           (pb_std/pb(ipts))*gs_circ(icir)*ratio_H2O_to_CO2



                      !! Wind speed is computed (used for the aerodynamic resistance, considered constant over 2 time steps for this estimation)

                      speed(ipts)=MIN(min_wind,SQRT(u(ipts)*u(ipts) + v(ipts)*v(ipts)))

                      !! The next value of the calculation is computed (see technical note for more details).
                      IF (printlev_loc >= 4) THEN 
                         WRITE(numout,*) "nextcalc i =",i         
                         WRITE(numout,*) "ipts =",ipts
                         WRITE(numout,*) "ivm =",ivm
                         WRITE(numout,*) "time step=",kjit
                         CALL flush(numout)                
                      ENDIF

                      transpir_test(ipts,ivm) = transpir_circ(icir)/&!MAX(min_sechiba*min_sechiba, transpir_circ(icir)/&    
                           (dt_sechiba * kilo_to_unit* soiltile(ipts,pref_soil_veg(ivm)) * vegtot(ipts)) * & 
                           delta(ipts)*((1+(veget(ipts,ivm)/veget_max(ipts,ivm))*(tq_cdrag(ipts)*speed(ipts))/g_stom(ipts,ivm))/&
                           (1+(veget(ipts,ivm)/veget_max(ipts,ivm))*(tq_cdrag(ipts)*speed(ipts)*&
                           f_psi_1(ipts,ivm))/(delta(ipts)*g_stom(ipts,ivm)*f_psi_2(ipts,ivm)))) !)

                      ! Debug
                      IF (printlev_loc >= 4) THEN
                         WRITE(numout,*) "transpir_test(ipts,ivm)", transpir_test(ipts,ivm)
                         WRITE(numout,*) "transpir(ipts,ivm)", transpir(ipts,ivm)
                         WRITE(numout,*) "g_stom(ipts,ivm)",g_stom(ipts,ivm) 
                         WRITE(numout,*) "f_psi_1(ipts,ivm)", f_psi_1(ipts,ivm) 
                         WRITE(numout,*) "f_psi_2(ipts,ivm)",f_psi_2(ipts,ivm) 
                         WRITE(numout,*) "transpir_test estimate factor", &
                           delta(ipts)*((1+(veget(ipts,ivm)/veget_max(ipts,ivm))* &
                           (tq_cdrag(ipts)*speed(ipts))/g_stom(ipts,ivm)) / &
                           (1+(veget(ipts,ivm)/veget_max(ipts,ivm))*(tq_cdrag(ipts)*speed(ipts) * &
                           f_psi_1(ipts,ivm))/(delta(ipts)*g_stom(ipts,ivm)*f_psi_2(ipts,ivm)))) 
                      ENDIF
                      !-

                   ELSE 
                      transpir_test(ipts,ivm) = zero!min_sechiba*min_sechiba
                   ENDIF

                ENDIF


                !! As the value of the water contents will differ between the first iteration calculation and the one for the future estimation,
                !! a quick water budget is implemented. In the first iteration, Fsup_temp and Finf_temp are equal to 0 so mc_out_temp = mc_out.
                !! In the other iterations, the transpiration is removed from the layers.


                DO isl= 1, nslm

                   IF (SUM(Fi_temp(ipts,ivm,:,jst)) .GT. min_sechiba) THEN

                      IF (isl .LE. lim_layer) THEN

                         mc_out_temp(ipts,isl,jst) = mc_out(ipts,isl,jst)-Fsup_temp(ipts,ivm)* &
                              Fi_temp(ipts,ivm,isl,jst) / SUM(Fi_temp(ipts,ivm,:,jst))*&
                              dt_sechiba/SUM(dlh(1:lim_layer))

                      ELSE

                         mc_out_temp(ipts,isl,jst) = mc_out(ipts,isl,jst)-Finf_temp(ipts,ivm)* &
                              Fi_temp(ipts,ivm,isl,jst) / SUM(Fi_temp(ipts,ivm,:,jst))*&
                              dt_sechiba/SUM(dlh(lim_layer+1:))

                      ENDIF

                   ELSE

                      mc_out_temp(ipts,isl,:) = mc_out(ipts,isl,:)

                   ENDIF

                ENDDO

                IF ((printlev_loc .GT. 3).AND.(ipts.EQ.test_grid).AND. &
                     ANY(mc_out_temp(ipts,:,jst)-mcr_sup(ipts).LT.zero)) THEN
                   WRITE(numout,*) 'mc_out_temp-mcr < 0'
                   CALL flush(numout)
                ENDIF


                !! 0. Calculation of the dynamic resistances:

                !! The resistances are either dynamic, relying on a sigmoid function as in Pammenter and Vander Willingen 1998, or fixed. When the resistances are dynamic, 
                !! their values depend on the value of the water potential at the stage below. If this value is below a precise
                !! value (50% decay value), the value of the resistance increases in order to limitate the water transport.
                !! If the resistances are fixed, the values are relying on Tuzet et al. (2017).


                IF (ok_hydrol_arch_dyn_res) THEN
                   ! ++++ CHECK ++++
                   ! quickfix for resistances and potential going to too
                   ! extreme values -> max set to the resistances.
                   ! Does not apply to root and storage resistances which seem to
                   ! cause less problems
                   ! Temporary. Fortunately, when mortality is implemented, it
                   ! could be removed
                   
                   ! +++ CHECK +++
                   ! LGC: Use pipe_density instead ? 

                   IF (ok_sap_feedback) THEN                   
                     
                     IF (ipts==test_grid.AND.ivm==test_pft.AND.printlev_loc>=3) THEN
                        WRITE(numout,*) 'ks_max(test_pft)', ks_max(ivm)
                        WRITE(numout,*) 'cc_bm(sapabove)', circ_class_biomass(ipts,ivm,icir,isapabove,icarbon)
                        WRITE(numout,*) 'cc_height(icir)', circ_height(icir)
                        CALL flush(numout)
                     ENDIF 
           
                     ! Special case for xylem up and low in which we simulate embolism
                     ! We consider all pipes/path lengths to leaves (on which we
                     ! integrate conductivity) equal
                     ! Another underlying hypothesis is that sapwood area is conserved
                     ! along the tree (Da Vinci's rule).

                     ! +++ CHECK +++
                     ! However in calculating equivalent
                     ! conductance in the below commented out way (* circ_class_n)i we overlook the capacitance  
                     ! Then we should implement an equivalent resistance formula that
                     ! depends on system's frequency.. --> Hence, keep to tree scale
                     ! If we upscaled to stand level
                     !k_xylem_up_max(ipts,ivm,icir) = (ks_max(ivm)*circ_class_biomass(ipts,ivm,icir,isapabove,icarbon)*crown_width/circ_height(icir)) / &
                     !                          ((circ_height(icir) ** 2) * sapwood_density(ivm)) * circ_class_n(ipts,ivm,icir)
                     ! +++++++++++++++++++
                     ! Doing on tree scale here
                     k_xylem_up_max(ipts,ivm,icir)  = (ks_max(ivm)*circ_class_biomass(ipts,ivm,icir,isapabove,icarbon)* &
                         crown_width(icir)/circ_height(icir)) / ((crown_width(icir) ** 2) * sapwood_density(ivm))
                     k_xylem_low_max(ipts,ivm,icir) = (ks_max(ivm)*circ_class_biomass(ipts,ivm,icir,isapabove,icarbon)* &
                        (circ_height(icir)-crown_width(icir))/circ_height(icir)) / &
                        (((circ_height(icir)-crown_width(icir)) ** 2) * sapwood_density(ivm))
                     
                     IF (.NOT. ok_tuzet_indiv) THEN
                        k_xylem_up_max(ipts,ivm,icir)  = k_xylem_up_max(ipts,ivm,icir) * circ_class_n(ipts,ivm,icir)
                        k_xylem_low_max(ipts,ivm,icir) = k_xylem_low_max(ipts,ivm,icir) * circ_class_n(ipts,ivm,icir)
                     ENDIF

                     ! +++ CHECK +++
                     ! We have a first "curtailement" of Res_xylem_up/low and Res_mesophyll because
                     ! they showed numerical problems (when potentials go too
                     ! low in the EXP). We keep the equations consistent. We
                     ! don't count what results from a minimum applied to the
                     ! conductance as a resistance curtailment  

                     ! xylem up
                     IF (k_xylem_up_max(ipts,ivm,icir) .GT. 1/undef_sechiba) THEN
                        Res_xylem_up(ipts,ivm,icir) = MIN(undef_sechiba, &
                                                     (1+EXP(a_wood(ivm)*(psi_xylem_trunk_temp(ipts,ivm)-psi_50_res_wood(ivm)))) &
                                                     / k_xylem_up_max(ipts,ivm,icir) )
                     ELSE
                        Res_xylem_up(ipts,ivm,icir) = undef_sechiba
                        k_xylem_up_max(ipts,ivm,icir) = 1/undef_sechiba
                        WRITE(numout,*) 'WARNING Xylem up conductance is 0; sapwood probably dead'
                        IF (ipts==test_grid.AND.ivm==test_pft) WRITE(numout,*) 'cc_bm(isapabove,icarbon)', &
                           circ_class_biomass(ipts,ivm,icir,isapabove,icarbon) 
                     ENDIF

                     ! xylem low
                     IF (k_xylem_low_max(ipts,ivm,icir) .GT. 1/undef_sechiba) THEN
                        Res_xylem_low(ipts,ivm,icir) = MIN(undef_sechiba, &
                                                     (1+EXP(a_wood(ivm)*(psi_xylem_collar_temp(ipts,ivm)-psi_50_res_wood(ivm)))) &
                                                     / k_xylem_low_max(ipts,ivm,icir) )
                     ELSE
                        Res_xylem_low(ipts,ivm,icir) = undef_sechiba
                        k_xylem_low_max(ipts,ivm,icir) = 1/undef_sechiba
                        WRITE(numout,*) 'WARNING Xylem low conductance is 0; sapwood probably dead'
                        IF (ipts==test_grid.AND.ivm==test_pft) WRITE(numout,*) 'cc_bm(isapbelow,icarbon)', &
                           circ_class_biomass(ipts,ivm,icir,isapbelow,icarbon) 
                     ENDIF
                   
                   ELSE 
                     Res_xylem_up(ipts,ivm,icir) = MIN(undef_sechiba, &
                        (1+EXP(a_wood(ivm)*(psi_xylem_trunk_temp(ipts,ivm)-psi_50_res_wood(ivm))))/k_wood_max(ivm))
                     Res_xylem_low(ipts,ivm,icir) = MIN(undef_sechiba, &
                        (1+EXP(a_wood(ivm)*(psi_xylem_collar_temp(ipts,ivm)-psi_50_res_wood(ivm))))/k_wood_max(ivm))
                   ENDIF

                   Res_mesophyll(ipts,ivm,icir) = MIN(undef_sechiba, &
                      (1+EXP(a_leaf(ivm)*(psi_xylem_leaf_temp(ipts,ivm)-psi_50_res_leaf(ivm))))/k_leaf_max(ivm))
                   Res_sto_leaf(ipts,ivm) = (1+EXP(a_leaf(ivm)*(psi_sto_leaf_temp(ipts,ivm)-psi_50_res_leaf(ivm))))/k_leaf_max(ivm)
                   Res_sto_wood(ipts,ivm) = (1+EXP(a_wood(ivm)*(psi_sto_wood_temp(ipts,ivm)-psi_50_res_wood(ivm))))/k_wood_max(ivm)
                   Res_root_sup(ipts,ivm,icir) = (1+EXP(a_root(ivm)*(psi_root_sup_temp(ipts,ivm)-psi_50_res_root(ivm))))/k_root_max(ivm)
                   Res_root_inf(ipts,ivm,icir) = (1+EXP(a_root(ivm)*(psi_root_inf_temp(ipts,ivm)-psi_50_res_root(ivm))))/k_root_max(ivm)
                   ! ++++++++++++++++++++

                   IF ((printlev_loc >= 4).AND.(ivm==test_pft).AND.(ipts==test_grid)) THEN
                      WRITE(numout,*) 'Res_mesophyll', Res_mesophyll(ipts,ivm,icir)
                      WRITE(numout,*) 'psi_xylem_leaf', psi_xylem_leaf(ipts,ivm,:)
                      !WRITE(numout,*) 'Res_sto_leaf', Res_sto_leaf(ipts,ivm)
                      !WRITE(numout,*) 'psi_sto_leaf_temp', psi_sto_leaf_temp(ipts,ivm)
                      !WRITE(numout,*) 'Res_sto_wood', Res_sto_wood(ipts,ivm)
                      !WRITE(numout,*) 'psi_sto_wood_temp', psi_sto_wood_temp(ipts,ivm)
                      WRITE(numout,*) 'Res_xylem_up', Res_xylem_up(ipts,ivm,icir)
                      WRITE(numout,*) 'psi_xylem_trunk', psi_xylem_trunk(ipts,ivm,:)
                      WRITE(numout,*) 'Res_xylem_low', Res_xylem_low(ipts,ivm,icir)
                      WRITE(numout,*) 'psi_xylem_collar', psi_xylem_collar(ipts,ivm,:)
                      !WRITE(numout,*) 'Res_root_sup', Res_root_sup(ipts,ivm,icir)
                      !WRITE(numout,*) 'psi_root_sup_temp', psi_root_sup_temp(ipts,ivm)
                      !WRITE(numout,*) 'Res_root_inf', Res_root_inf(ipts,ivm,icir)
                      !WRITE(numout,*) 'psi_root_inf_temp', psi_root_inf_temp(ipts,ivm)
                      CALL flush(numout)
                   ENDIF

                   IF ((printlev_loc .EQ. 4).AND.(veget(ipts,ivm).GT.min_sechiba)) THEN 
                      ! Here we do writing for a basic counter in out_orchidee
                      ! files of the proportion of resistances calculated that are
                      ! 'curtailed'
                      counter=un
                      IF ((1+EXP(a_wood(ivm)*(psi_xylem_trunk_temp(ipts,ivm)-psi_50_res_wood(ivm))))&
                         /k_wood_max(ivm).GT.undef_sechiba) THEN 
                         WRITE(numout,*) 'resistance curtailed : xylem trunk'
                         counter=zero
                      ENDIF
                      IF ((1+EXP(a_leaf(ivm)*(psi_xylem_leaf_temp(ipts,ivm)-psi_50_res_leaf(ivm))))&
                         /k_leaf_max(ivm) .GT. undef_sechiba) THEN 
                         WRITE(numout,*) 'resistance curtailed : leaf'
                         counter=zero
                      ENDIF
                      IF ((1+EXP(a_wood(ivm)*(psi_xylem_collar_temp(ipts,ivm)-psi_50_res_wood(ivm))))&
                         /k_wood_max(ivm) .GT. undef_sechiba) THEN 
                         WRITE(numout,*) 'resistance curtailed : xylem collar'
                         counter=zero
                      ENDIF
                      IF (counter==un) THEN
                         WRITE(numout,*) 'resistance not curtailed'
                      ENDIF

                   ENDIF

                ELSE

                   Res_sto_leaf(ipts,ivm)=Rsto_leaf(ivm)
                   Res_sto_wood(ipts,ivm)=Rsto_wood(ivm)
                   Res_mesophyll(ipts,ivm,icir)=Rleaf(ivm)
                   Res_xylem_up(ipts,ivm,icir)=Rxylem(ivm)
                   Res_xylem_low(ipts,ivm,icir)=Rxylem(ivm)
                   Res_root_sup(ipts,ivm,icir)=Rroot_sup(ivm)
                   Res_root_inf(ipts,ivm,icir)=Rroot_inf(ivm)

                ENDIF
               
                ! We make the non dynamic resistances compatible with individual
                ! level Tuzet 
                IF ( ok_tuzet_indiv .AND. (circ_class_n(ipts,ivm,icir).GT.min_sechiba) ) THEN
                   Res_sto_leaf(ipts,ivm)=Res_sto_leaf(ipts,ivm)/circ_class_n(ipts,ivm,icir)
                   Res_sto_wood(ipts,ivm)=Res_sto_wood(ipts,ivm)/circ_class_n(ipts,ivm,icir)
                   Res_mesophyll(ipts,ivm,icir)=Res_mesophyll(ipts,ivm,icir)/circ_class_n(ipts,ivm,icir)
                   Res_root_sup(ipts,ivm,icir)=Res_root_sup(ipts,ivm,icir)/circ_class_n(ipts,ivm,icir)
                   Res_root_inf(ipts,ivm,icir)=Res_root_inf(ipts,ivm,icir)/circ_class_n(ipts,ivm,icir)

                   IF (.NOT. ok_hydrol_arch_dyn_res .OR. .NOT. ok_sap_feedback) THEN
                      Res_xylem_up(ipts,ivm,icir)=Res_xylem_up(ipts,ivm,icir)/circ_class_n(ipts,ivm,icir)
                      Res_xylem_low(ipts,ivm,icir)=Res_xylem_low(ipts,ivm,icir)/circ_class_n(ipts,ivm,icir)
                   ENDIF                
                ELSEIF ( ok_tuzet_indiv .AND. circ_class_n(ipts,ivm,icir).LT.min_sechiba) THEN! if no trees
                   Res_sto_leaf(ipts,ivm)=undef_sechiba
                   Res_sto_wood(ipts,ivm)=undef_sechiba
                   Res_mesophyll(ipts,ivm,icir)=undef_sechiba
                   Res_root_sup(ipts,ivm,icir)=undef_sechiba
                   Res_root_inf(ipts,ivm,icir)=undef_sechiba
                ENDIF
 
                ! Calculating percentage loss of conductivity (PLC), xylem low
                ! and up being in series
                IF ( ok_vessel_mortality ) THEN
                   vessel_loss(ipts,ivm,icir) = zero
                   IF ( ok_hydrol_arch_dyn_res ) THEN
                      IF (.NOT. Res_xylem_up(ipts,ivm,icir)+Res_xylem_low(ipts,ivm,icir).LT.min_stomate) THEN
                         IF (ok_sap_feedback) THEN
                            vessel_loss(ipts,ivm,icir) = 1 - (1/k_xylem_low_max(ipts,ivm,icir)+1/k_xylem_up_max(ipts,ivm,icir)) &
                               / (Res_xylem_up(ipts,ivm,icir)+Res_xylem_low(ipts,ivm,icir))
                         ELSE
                            vessel_loss(ipts,ivm,icir) = 1 - (2/k_wood_max(ivm)) / (Res_xylem_up(ipts,ivm,icir)+Res_xylem_low(ipts,ivm,icir))
                         ENDIF
                      ELSE
                         WRITE(numout,*) 'WARNING (Res_xylem_up + Res_xylem_low) is too low'
                         IF (ipts==test_grid.AND.ivm==test_pft) WRITE(numout,*) 'Res_xylem_up(ipts,ivm,icir)', Res_xylem_up(ipts,ivm,icir)
                      ENDIF
                   ENDIF
                ENDIF

                !! 1. Determination of the soil water content in both soil layer, superficial and inferior

                !! The first step of the resolution is to define two layers of soil thanks to the 11 ones present in hydrol. To do so,
                !! the soil water content of both soil layers is defined thanks to the ones of the 11 soil layers weighted by the 
                !! height of each layer.

                !! The 2-layers scheme for the soil resolution follows the scheme of Tuzet et al. (2017). This scheme can be changed.
                !! However, adding soil layers adds a lot of complexity for the resolution of the root absorption. Moreover, it is
                !! maybe not relevant to define thin soil layers (~1-5 mm) to solve the root absorption as the roots are thicker. 
                !! Consequently, it has been decided for now to limitate the number of layers to only 2.

                DO isl = 1, nslm      ! Loop over the soil layers

                   DO ist = 1, nstm   ! Loop over the soil tiles

                      IF (isl .LE. lim_layer) THEN  ! lim_layer is a parameter that define the node which delimits the superficial and inferior layers

                         mc_sup_temp(ipts,ist) = mc_sup_temp(ipts,ist) + mc_out_temp(ipts,isl,ist) * dh(isl)/SUM(dh(1:lim_layer))

                      ELSE

                         mc_inf_temp(ipts,ist) = mc_inf_temp(ipts,ist) + mc_out_temp(ipts,isl,ist) * dh(isl)/SUM(dh(lim_layer+1:))

                      ENDIF

                   ENDDO

                ENDDO

                IF ((printlev_loc .GT. 4) .AND. ((mc_sup_temp(ipts,jst) - mcr_sup(ipts) .LT. zero) .OR. &
                     (mc_inf_temp(ipts,jst) - mcr_inf(ipts) .LT. zero))) THEN
                   WRITE(numout,*) 'WARNING: mc less than residual'
                ENDIF

                !! 2. Calculation of the fluxes inside the vegetation:


                !! The first step of the calculation is to determine the water transport fluxes inside the vegetation.
                !! Those fluxes depend on the hydraulic architecture implemented (with or without water storage). 

                !! Without water storage, the split of the transpiration only occurs at the collar interface. At this
                !! stage, the transpiration flux has to be split between the superficial and inferior root absorptions.
                !! The split is done thanks to a Kirchhoff's current law at the collar level. The split depends on the 
                !! value of the root water potentials at the previous time step.

                !! With water storage, the differential equations that control the water dynamic inside the water storages 
                !! have to be solved first. As this resolution brings a lot of instability, two predictor/corrector 
                !! schemes have been implemented for each water storage: Runge-Kutta 2 or Adams-Moulton order 3. However,
                !! Runge-Kutta 2 scheme seems to be sufficiently effective. We kept Adams-Moulton order 3 scheme in case of
                !! need. The aim of the scheme is to predict the value of the water potential at the new time step in order
                !! to mimic an implicit scheme and avoid instabilities. This permits to calculate the amount of water
                !! extracted from the storage compartment and, thus, calculate the flux from the storage. Once all the 
                !! fluxes are calculated, the split at the collar level is the same as without storage.

                IF (ok_hydrol_arch_storage) THEN

                   Vmax_wood(ipts,ivm) = zero 
                   Vmax_leaf(ipts,ivm) = zero
                   Vr_wood(ipts,ivm)   = zero 
                   Vr_leaf(ipts,ivm)   = zero

                   !! The water flux supplied by the water storage compartments are calculated by solving the two differential equations.
                   IF (soiltile(ipts,jst) * vegtot(ipts) .GT. min_sechiba) THEN
                      IF (circ_class_biomass(ipts,ivm,icir,ileaf,icarbon)*circ_class_n(ipts,ivm,icir) .GT. min_sechiba) THEN
                         !! If there is no leaf biomass, the storage is disabled.
                         !! The previous value of the storage water potential is saved for the prediction/correction scheme

                         !! psi_sto_leaf_save(ipts,ivm) = psi_sto_leaf(ipts,ivm)

                         !! The leaf biomass permits to calculate the maximum amount of water that can be stored in the leaves.
                         !! We consider that 80% of the leaf biomass permits to store water.

                         !! ++++ CHECK ++++ 
                         !! The factor leaf_wat_store_fac permits to get similar values to Tuzet et al. (2017) model (in the case of PFT6), this
                         !! can maybe outline a biomass underestimation. The factor should be removed.
                         Vmax_leaf(ipts,ivm) = deux*circ_class_biomass(ipts,ivm,icir,ileaf,icarbon)*0.8/(kilo_to_unit*ph2o)* &
                                               leaf_wat_store_fac(ivm)
                         
                         IF (.NOT. ok_tuzet_indiv) THEN
                           Vmax_leaf(ipts,ivm) = Vmax_leaf(ipts,ivm)*circ_class_n(ipts,ivm,icir)
                         ENDIF

                         !! +++++++++++++++

                         !! As discussed with Andrée Tuzet, the residual water storage is empirically set to 40% of the maximum water storage
                         Vr_leaf(ipts,ivm) = 0.4*Vmax_leaf(ipts,ivm)

                         !! Previous amount of water stored in the compartment, compartment capacitance and time constant are calculated
                         !! as in Tuzet et al. (2017) (Eq. 7 and 8). They will permit to calculate the amount of water removed from or added to the compartment.
                         V_leaf_save(ipts,ivm) = (1+exp(lambda_leaf(ivm)*psi_ref_sto_leaf(ivm)))&
                              /(1+exp(lambda_leaf(ivm)*(-psi_sto_leaf(ipts,ivm,icir)+psi_ref_sto_leaf(ivm))))&
                              *(Vmax_leaf(ipts,ivm)-Vr_leaf(ipts,ivm))+Vr_leaf(ipts,ivm) 

                         C_leaf(ipts,ivm) = lambda_leaf(ivm)*(((1+exp(lambda_leaf(ivm)*psi_ref_sto_leaf(ivm)))*Vmax_leaf(ipts,ivm)&
                              +Vr_leaf(ipts,ivm)*exp(lambda_leaf(ivm)*psi_ref_sto_leaf(ivm)))-V_leaf_save(ipts,ivm))&
                              /(1+exp(lambda_leaf(ivm)*psi_ref_sto_leaf(ivm)))*&
                              ((1+exp(lambda_leaf(ivm)*psi_ref_sto_leaf(ivm)))&
                              /(1+exp(lambda_leaf(ivm)*(-psi_sto_leaf(ipts,ivm,icir)+psi_ref_sto_leaf(ivm)))))

                         tau_leaf(ipts,ivm) = C_leaf(ipts,ivm) * (Res_xylem_up(ipts,ivm,icir) + Res_sto_leaf(ipts,ivm))

                         !! Here comes the resolution for the new time step. the prediction/correction scheme permits to mimic an implicit 
                         !! resolution of the differential equation (Tuzet et al. (2017) Eq.5).
                         !! The values of the variables that are not yet calculated are relaxed thanks to the approximations calculated at the beginning
                         !! of the subroutine.

                         !! The Adams-Moulton Order 3 method can be used but does not really improve the results:

                         !! Adams-Moulton Order 3 method:

                         psi_sto_leaf_predict(ipts,ivm) = psi_sto_leaf(ipts,ivm,icir)+dt_sechiba*(1/tau_leaf(ipts,ivm))*&
                              & (psi_xylem_trunk(ipts,ivm,icir)-Res_xylem_up(ipts,ivm,icir)*transpir_test(ipts,ivm)-psi_sto_leaf(ipts,ivm,icir))

                         func_n(ipts,ivm) = (1/tau_leaf(ipts,ivm))*&
                              & (psi_xylem_trunk(ipts,ivm,icir)-Res_xylem_up(ipts,ivm,icir)*transpir_test(ipts,ivm)-psi_sto_leaf(ipts,ivm,icir))

                         func_predict(ipts,ivm) = (1/tau_leaf(ipts,ivm))*&
                              & (psi_xylem_trunk(ipts,ivm,icir)-Res_xylem_up(ipts,ivm,icir)*transpir_test(ipts,ivm)-psi_sto_leaf_predict(ipts,ivm))

                         func_n_1(ipts,ivm) = (1/tau_leaf(ipts,ivm))*&
                              & (psi_xylem_trunk(ipts,ivm,icir)-Res_xylem_up(ipts,ivm,icir)*transpir_test(ipts,ivm)-psi_sto_leaf_save(ipts,ivm,icir))

                         psi_sto_leaf_save(ipts,ivm,icir) = psi_sto_leaf(ipts,ivm,icir)

                         psi_sto_leaf_temp(ipts,ivm) = MAX(psi_xylem_leaf(ipts,ivm,icir),&
                              MIN(0.,psi_sto_leaf_save(ipts,ivm,icir) + dt_sechiba/12.*(5.*func_predict(ipts,ivm)+8.*func_n(ipts,ivm)-func_n_1(ipts,ivm))))


                         !! The other method implemented (less voracious in terms of memory) is the method of Runge-Kutta 2

                         !! Runge-Kutta 2 Method:

!!$                   psi_sto_leaf_predict(ipts,ivm) = psi_sto_leaf_save(ipts,ivm)+dt_sechiba*(1/tau_leaf(ipts,ivm))*&
!!$                        (psi_xylem_trunk(ipts,ivm)-Res_xylem_up(ipts,ivm,icir)*transpir(ipts,ivm)&
!!$                        /(dt_sechiba * kilo_to_unit* soiltile(ipts,jst) * vegtot(ipts))-psi_sto_leaf_save(ipts,ivm))
!!$
!!$                   psi_sto_leaf(ipts,ivm) = MIN(0.,psi_sto_leaf_save(ipts,ivm)+dt_sechiba/deux*&
!!$                        ((1/tau_leaf(ipts,ivm))*(psi_xylem_trunk(ipts,ivm)-Res_xylem_up(ipts,ivm,icir)*transpir(ipts,ivm)&
!!$                        /(dt_sechiba * kilo_to_unit* soiltile(ipts,jst) * vegtot(ipts))-psi_sto_leaf_save(ipts,ivm))+&
!!$                        (1/tau_leaf(ipts,ivm))*(psi_xylem_trunk(ipts,ivm)-Res_xylem_up(ipts,ivm,icir)*transpir(ipts,ivm)&
!!$                        /(dt_sechiba * kilo_to_unit* soiltile(ipts,jst) * vegtot(ipts))-psi_sto_leaf_predict(ipts,ivm))))


                         !! Explicit method is given as example but leads to strong instabilities when water stress appears.

                         !! Explicit method:

!!$                   psi_sto_leaf(ipts,ivm) = psi_sto_leaf(ipts,ivm) * exp(-dt_sechiba/tau_leaf(ipts,ivm)) &
!!$                        + (psi_xylem_trunk(ipts,ivm)-Res_xylem_up(ipts,ivm,icir)*transpir(ipts,ivm)&
!!$                        /(dt_sechiba * kilo_to_unit* soiltile(ipts,jst) * vegtot(ipts)))*(1-exp(-dt_sechiba/tau_leaf(ipts,ivm)))


                         !! The value of the storage water potential obtained permits to calculate the new amount of water
                         !! inside the storage compartment, and then the flux removed from or added to the storage compartment.

                         V_leaf(ipts,ivm) = MIN(Vmax_leaf(ipts,ivm),(1+exp(lambda_leaf(ivm)*psi_ref_sto_leaf(ivm)))&
                              /(1+exp(lambda_leaf(ivm)*(-psi_sto_leaf_temp(ipts,ivm)+psi_ref_sto_leaf(ivm))))&
                              *(Vmax_leaf(ipts,ivm)-Vr_leaf(ipts,ivm))+Vr_leaf(ipts,ivm) )

                         Fsto_leaf(ipts,ivm) = (V_leaf(ipts,ivm) - V_leaf_save(ipts,ivm))/dt_sechiba

                         IF ((Fsto_leaf(ipts,ivm) .GT. 0.) .AND. (transpir_circ(icir) .GT. min_sechiba)) THEN

                            Fsto_leaf(ipts,ivm)=0.
                            psi_sto_leaf_save(ipts,ivm,icir) = psi_sto_leaf(ipts,ivm,icir)
                            psi_sto_leaf_temp(ipts,ivm)=psi_sto_leaf(ipts,ivm,icir)

                         ENDIF

                         IF ((Fsto_leaf(ipts,ivm) .LT. 0.) .AND. (transpir_circ(icir) .LT. min_sechiba)) THEN

                            Fsto_leaf(ipts,ivm)=0.
                            psi_sto_leaf_save(ipts,ivm,icir) = psi_sto_leaf(ipts,ivm,icir)
                            psi_sto_leaf_temp(ipts,ivm)=psi_sto_leaf(ipts,ivm,icir)

                         ENDIF

                      ELSE ! If there is no leaf biomass
                         !! Storage water potential is set to the xylem water potential

                         psi_sto_leaf_save(ipts,ivm,icir) = psi_sto_leaf(ipts,ivm,icir)
                         psi_sto_leaf_temp(ipts,ivm)=psi_xylem_leaf(ipts,ivm,icir)
                         Fsto_leaf(ipts,ivm) = 0.

                      ENDIF

                   ELSE ! If there is no vegetation
                      !! Same as if there is no biomass

                      psi_sto_leaf_save(ipts,ivm,icir) = psi_sto_leaf(ipts,ivm,icir)
                      psi_sto_leaf_temp(ipts,ivm)=psi_xylem_leaf(ipts,ivm,icir)
                      Fsto_leaf(ipts,ivm) = 0.

                   ENDIF

                   !! If vegetation is not present, there is no transpiration. If there is no biomass, there is no storage flux
                   IF ((soiltile(ipts,jst) * vegtot(ipts) .LT. min_sechiba).AND.transpir_test(ipts,ivm).GT.min_sechiba) THEN
                      WRITE(numout,*) 'WARNING :  no vegetation should mean no transpir_test, transpir_test =', &
                        transpir_test(ipts,ivm)
                   ENDIF
                   IF ((circ_class_biomass(ipts,ivm,icir,ileaf,icarbon)*circ_class_n(ipts,ivm,icir) .LT. min_sechiba) &
                       .AND. Fsto_leaf(ipts,ivm).GT.min_sechiba ) THEN
                      WRITE(numout,*) 'WARNING :  no biomass should mean no Fsto_leaf, Fsto_leaf =', Fsto_leaf(ipts,ivm)
                   ENDIF

                   ! If more water is supplied from storage than is transpired,
                   ! apply the following correction such that all water needed
                   ! is supplied by storage (anyway, this case seems to happen
                   ! only at very low flux values)
                   IF ( (transpir_test(ipts,ivm).LT.-Fsto_leaf(ipts,ivm)) .AND. (Fsto_leaf(ipts,ivm).LT.zero) ) THEN
                       Fsto_leaf(ipts,ivm) = -transpir_test(ipts,ivm)
                       V_leaf(ipts,ivm) = Fsto_leaf(ipts,ivm)*dt_sechiba + V_leaf_save(ipts,ivm)
                   ENDIF

                   !! The flux below is calculated thanks to the two fluxes above.
                   Fxylem_up(ipts,ivm) = transpir_test(ipts,ivm) + Fsto_leaf(ipts,ivm)


                   IF ((circ_class_biomass(ipts,ivm,icir,isapabove,icarbon)&
                        +circ_class_biomass(ipts,ivm,icir,isapbelow,icarbon))*circ_class_n(ipts,ivm,icir) .GT. min_sechiba) THEN
                      !! If there is no sapwood biomass, the storage compartment is disabled (should not be the case)

                      !! Previous value of the wood water storage is saved for the future prediction/correction scheme

!!$                psi_sto_wood_save(ipts,ivm) = psi_sto_wood(ipts,ivm)

                      !! As for the leaves, the maximum wood water storage is calculated thanks to the sapwwod biomass
                      !! We consider that 80% of the sapwood biomass permits to store water.

                      ! ++++ CHECK ++++++++
                      ! -- 3.5 is a patch factor to come back to Tuzet et al. (2017) values of maximum storage (PFT6).
                      ! This factor should be removed. It maybe highlights a biomass underestimation.
                      ! -- 0.8 comes from the hypothesis that 80% of the sapwood biomass permits to store water
                      ! --'deux' originates from the hypothesis that 1 gC stores 2g of water

                      Vmax_wood(ipts,ivm) = deux*(circ_class_biomass(ipts,ivm,icir,isapabove,icarbon)&
                           +circ_class_biomass(ipts,ivm,icir,isapbelow,icarbon))&
                           *0.8/(kilo_to_unit*ph2o)* wood_wat_store_fac(ivm)
                      IF (.NOT. ok_tuzet_indiv) THEN
                         Vmax_wood(ipts,ivm) = Vmax_wood(ipts,ivm)*circ_class_n(ipts,ivm,icir) 
                      ENDIF
                      ! +++++++++++++++++++
  
                      !! Residual water storage is set to 40% of the maximum water storage as discussed with Andrée Tuzet.

                      Vr_wood(ipts,ivm) = 0.4*Vmax_wood(ipts,ivm)

                      !! Previous values of the amount of stored water, compartment capacitance and time constant are calculated
                      !! following Tuzet et al. (2017) (Eq. 7 and 8)

!!$                V_wood_save(ipts,ivm) = (1+exp(lambda_wood(ivm)*psi_ref_sto_wood(ivm)))&
!!$                     /(1+exp(lambda_wood(ivm)*(-psi_sto_wood_save(ipts,ivm)+psi_ref_sto_wood(ivm))))&
!!$                     *(Vmax_wood(ipts,ivm)-Vr_wood(ipts,ivm))+Vr_wood(ipts,ivm)
!!$
!!$                C_wood(ipts,ivm) = lambda_wood(ivm)*(((1+exp(lambda_wood(ivm)*psi_ref_sto_wood(ivm)))*Vmax_wood(ipts,ivm)&
!!$                     +Vr_wood(ipts,ivm)*exp(lambda_wood(ivm)*psi_ref_sto_wood(ivm)))-V_wood_save(ipts,ivm))&
!!$                     /(1+exp(lambda_wood(ivm)*psi_ref_sto_wood(ivm)))*&
!!$                     ((1+exp(lambda_wood(ivm)*psi_ref_sto_wood(ivm)))&
!!$                     /(1+exp(lambda_wood(ivm)*(-psi_sto_wood_save(ipts,ivm)+psi_ref_sto_wood(ivm)))))
!!$
!!$                tau_wood(ipts,ivm) =  C_wood(ipts,ivm) * (Res_xylem_low(ipts,ivm,icir) + Res_sto_wood(ipts,ivm))

!!$                   IF (this_tstep_calc(ipts,ivm)) THEN
!!$                      psi_sto_wood_temp(ipts,ivm) = psi_sto_wood(ipts,ivm)
!!$                   ENDIF

                      V_wood_save(ipts,ivm) = (1+exp(lambda_wood(ivm)*psi_ref_sto_wood(ivm)))&
                           /(1+exp(lambda_wood(ivm)*(-psi_sto_wood(ipts,ivm,icir)+psi_ref_sto_wood(ivm))))&
                           *(Vmax_wood(ipts,ivm)-Vr_wood(ipts,ivm))+Vr_wood(ipts,ivm)
                      ! +++ CHECK +++ (cleaning) : two lines seem redundant
                      C_wood(ipts,ivm) = lambda_wood(ivm)*(((1+exp(lambda_wood(ivm)*psi_ref_sto_wood(ivm)))*Vmax_wood(ipts,ivm)&
                           +Vr_wood(ipts,ivm)*exp(lambda_wood(ivm)*psi_ref_sto_wood(ivm)))-V_wood_save(ipts,ivm))&
                           /(1+exp(lambda_wood(ivm)*psi_ref_sto_wood(ivm)))*&
                           ((1+exp(lambda_wood(ivm)*psi_ref_sto_wood(ivm)))&
                           /(1+exp(lambda_wood(ivm)*(-psi_sto_wood(ipts,ivm,icir)+psi_ref_sto_wood(ivm)))))
                      ! ++++++++++++++++++++++++
                      tau_wood(ipts,ivm) =  C_wood(ipts,ivm) * (Res_xylem_low(ipts,ivm,icir) + Res_sto_wood(ipts,ivm))

                      !! Here comes the resolution for the new time step. the prediction/correction scheme permits to mimic an implicit 
                      !! resolution of the differential equation (Tuzet et al. (2017) Eq.5)

                      !! The Adams-Moulton Order 3 method can be used but does not really improve the results:

                      !! Adams-Moulton Order 3 method:

                      psi_sto_wood_predict(ipts,ivm) =  psi_sto_wood(ipts,ivm,icir)+dt_sechiba*(1/tau_wood(ipts,ivm))*&
                           & (psi_xylem_collar(ipts,ivm,icir)-Res_xylem_low(ipts,ivm,icir)*Fxylem_up(ipts,ivm)-psi_sto_wood(ipts,ivm,icir))

                      func_n(ipts,ivm) = (1/tau_wood(ipts,ivm))*&
                           & (psi_xylem_collar(ipts,ivm,icir)-Res_xylem_low(ipts,ivm,icir)*Fxylem_up(ipts,ivm)-psi_sto_wood(ipts,ivm,icir))

                      func_predict(ipts,ivm) = (1/tau_wood(ipts,ivm))*&
                           & (psi_xylem_collar(ipts,ivm,icir)-Res_xylem_low(ipts,ivm,icir)*Fxylem_up(ipts,ivm)-psi_sto_wood_predict(ipts,ivm))

                      func_n_1(ipts,ivm) = (1/tau_wood(ipts,ivm))*&
                           & (psi_xylem_collar(ipts,ivm,icir)-Res_xylem_low(ipts,ivm,icir)*Fxylem_up(ipts,ivm)-psi_sto_wood_save(ipts,ivm,icir))

                      psi_sto_wood_save(ipts,ivm,icir) = psi_sto_wood(ipts,ivm,icir)

                      psi_sto_wood_temp(ipts,ivm) = MAX(psi_xylem_trunk(ipts,ivm,icir),&
                           MIN(0.,psi_sto_wood_save(ipts,ivm,icir) + dt_sechiba/12.*(5.*func_predict(ipts,ivm)+8.*func_n(ipts,ivm)-func_n_1(ipts,ivm))))



                      !! The other method implemented (less voracious in terms of memory) is the method of Runge-Kutta 2

                      !! Runge-Kutta 2 Method:

!!$                psi_sto_wood_predict(ipts,ivm) = psi_sto_wood_save(ipts,ivm)+dt_sechiba*(1/tau_wood(ipts,ivm))*&
!!$                     (psi_xylem_collar(ipts,ivm)-Res_xylem_low(ipts,ivm,icir)*Fxylem_up(ipts,ivm)-psi_sto_wood_save(ipts,ivm))
!!$
!!$                psi_sto_wood(ipts,ivm)=MIN(0.,psi_sto_wood_save(ipts,ivm)+dt_sechiba/deux*&
!!$                     ((1/tau_wood(ipts,ivm))*(psi_xylem_collar(ipts,ivm)-Res_xylem_low(ipts,ivm,icir)*Fxylem_up(ipts,ivm)-psi_sto_wood_save(ipts,ivm))+&
!!$                     (1/tau_wood(ipts,ivm))*(psi_xylem_collar(ipts,ivm)-Res_xylem_low(ipts,ivm,icir)*Fxylem_up(ipts,ivm)-psi_sto_wood_predict(ipts,ivm))))


                      !! Explicit method is given as example but leads to strong instabilities when water stress appears.

                      !! Explicit method:

!!$                psi_sto_wood(ipts,ivm) = psi_sto_wood(ipts,ivm) * exp(-dt_sechiba/tau_wood(ipts,ivm)) &
!!$                     + (psi_xylem_collar(ipts,ivm)-Res_xylem_low(ipts,ivm,icir)*Fxylem_up(ipts,ivm))*(1-exp(-dt_sechiba/tau_wood(ipts,ivm)))


                      !! The new value of the amount of stored water is calculated thanks to the prediction/correction value and permits
                      !! to calculate the flux removed from or added to the storage compartment.

                      V_wood(ipts,ivm) = MIN(Vmax_wood(ipts,ivm),(1+exp(lambda_wood(ivm)*psi_ref_sto_wood(ivm)))&
                           /(1+exp(lambda_wood(ivm)*(-psi_sto_wood_temp(ipts,ivm)+psi_ref_sto_wood(ivm))))&
                           *(Vmax_wood(ipts,ivm)-Vr_wood(ipts,ivm))+Vr_wood(ipts,ivm))

                      Fsto_wood(ipts,ivm) = (V_wood(ipts,ivm) - V_wood_save(ipts,ivm))/dt_sechiba


                      IF ((Fsto_wood(ipts,ivm) .GT. 0.) .AND. (transpir_circ(icir) .GT. min_sechiba)) THEN

                         Fsto_wood(ipts,ivm)=0.
                         psi_sto_wood_save(ipts,ivm,icir) = psi_sto_wood(ipts,ivm,icir)
                         psi_sto_wood_temp(ipts,ivm)=psi_sto_wood(ipts,ivm,icir)

                      ENDIF


                      IF ((Fsto_wood(ipts,ivm) .LT. 0.) .AND. (transpir_circ(icir) .LT. min_sechiba)) THEN

                         Fsto_wood(ipts,ivm)=0.
                         psi_sto_wood_save(ipts,ivm,icir) = psi_sto_wood(ipts,ivm,icir)
                         psi_sto_wood_temp(ipts,ivm)=psi_sto_wood(ipts,ivm,icir)

                      ENDIF

                   ELSE !! If there is no sapwood biomass

                      psi_sto_wood_save(ipts,ivm,icir)=psi_sto_wood(ipts,ivm,icir)
                      psi_sto_wood_temp(ipts,ivm)=psi_xylem_trunk(ipts,ivm,icir)
                      Fsto_wood(ipts,ivm) = 0.

                   ENDIF

                   IF ((circ_class_biomass(ipts,ivm,icir,isapabove,icarbon)+circ_class_biomass(ipts,ivm,icir,isapbelow,icarbon))&
                      *circ_class_n(ipts,ivm,icir) .LT. min_sechiba .AND. Fsto_wood(ipts,ivm).GT.min_sechiba) THEN
                      WRITE(numout,*) 'WARNING :  no biomass should mean no Fsto_wood, Fsto_wood =', Fsto_wood(ipts,ivm)
                   ENDIF
                      
                   ! If more water is supplied from storage than is transpired,
                   ! apply the following correction such that all water needed
                   ! is supplied by storage (anyway, this case seems to happen
                   ! only at very low flux values)
                   IF ( (Fxylem_up(ipts,ivm).LT.-Fsto_wood(ipts,ivm)) .AND. (Fsto_wood(ipts,ivm).LT.zero) ) THEN
                       Fsto_wood(ipts,ivm) = -Fxylem_up(ipts,ivm)
                       V_wood(ipts,ivm) = Fsto_wood(ipts,ivm)*dt_sechiba + V_wood_save(ipts,ivm)
                   ENDIF

                   Fxylem_low(ipts,ivm) = Fxylem_up(ipts,ivm) + Fsto_wood(ipts,ivm)

                   !! Finally, the total root absorption is equal to the water flux from the xylem at the collar level
                   !! and the xylem at the trunk level. This F_absorption ill be used for the resolution of the root
                   !! absorption

                   F_absorption_temp(ipts,ivm) = Fxylem_low(ipts,ivm)

                ELSE
                   !! If there is no storage, the flux of absorption corresponds to the flux of transpiration (with the unit conversion)

                   !! If vegetation is not present, there is no transpiration. 
                   IF ((soiltile(ipts,jst) * vegtot(ipts) .LT. min_sechiba).AND.transpir_test(ipts,ivm).GT.min_sechiba) THEN
                      WRITE(numout,*) 'WARNING :  no vegetation should mean no transpir_test, transpir_test =', transpir_test(ipts,ivm)
                   ENDIF

                   !! As there is no storage in the scheme, all the fluxes are the same at each stage of the architecture

                   Fxylem_up(ipts,ivm) = transpir_test(ipts,ivm)

                   Fxylem_low(ipts,ivm) = Fxylem_up(ipts,ivm)

                   F_absorption_temp(ipts,ivm) = Fxylem_low(ipts,ivm)

                ENDIF

                ! Debug
                IF ( (Fxylem_low(ipts,ivm) .LT. zero) .OR. ((printlev_loc>=4).AND.(ipts==test_grid).AND.(ivm==test_pft)) ) THEN 
                   IF (Fxylem_low(ipts,ivm) .LT. zero) WRITE(numout,*) 'WARNING : Fxylem_low < 0'
                   WRITE(numout,*) 'Fluxes : transpir_test', transpir_test(ipts,ivm)
                   WRITE(numout,*) '1transpir_cc', transpir_circ(:)
                   WRITE(numout,*) '1transpir', transpir(ipts,ivm)
                   WRITE(numout,*) 'Fxylem_up', Fxylem_up(ipts,ivm)
                   WRITE(numout,*) 'Fsto_wood(ipts,ivm)',Fsto_wood(ipts,ivm) 
                   WRITE(numout,*) 'Fxylem_low', Fxylem_low(ipts,ivm)
                ENDIF
                !-

                !! 3. Call to the root absorption resolution method:

                !! This part calls the routine that will be used to solve the root absorption.

                !! If the Asborption muffs are called, a radial resolution of Richards equation around a fictive root is launched
                !! The aim of the resolution is to model the gradient of water content that appears around the roots when they 
                !! are absorbing water. This water content gradient controls the amount of water absorbed by the roots. As the 
                !! soil water content is decreasing, the soil hydraulic conductivity is decreasing and it becomes more difficult
                !! for the plant to absorb water. Consequently, the root water potential will decrease in order to permit the 
                !! absorption. the output of the subroutine is this root water potential. The resolution of Richard's equation
                !! relies on the same method as for the vertical resolution of the hydrol.f90 module. The muff model is taken from 
                !! Tuzet et al. (2003).

                !! If the absorption muffs are not calculated, the root absorption is modelled thanks to a dynamic resistance
                !! which increases when the soil water content decreases. The model of the resistance is taken from Bonan et al. (2014).

                IF (ok_hydrol_arch_muff) THEN

                   CALL hydrol_hydraulic_arch_tuzet_muff(kjit, kjpindex, ipts, ivm, icir, soiltile, veget_max,  &
                        njsc, ks, nvan, avan, F_absorption_temp, circ_class_biomass, circ_class_n,             &
                        Res_root_sup, Res_root_inf, root_profile, Fsup_temp, Finf_temp, psi_soil_sup_temp, psi_soil_inf_temp, &
                        psi_root_sup_temp, psi_root_inf_temp, mc_sup_temp, mc_inf_temp, mc_i_sup_temp, mc_i_inf_temp)
                ELSE 

                   CALL hydrol_hydraulic_arch_tuzet_resist(kjpindex, ipts, ivm, icir, soiltile, njsc, ks,    &
                        nvan, avan, F_absorption_temp, circ_class_biomass, circ_class_n, mc_sup_temp, mc_inf_temp,    &
                        Res_root_sup, Res_root_inf, root_profile, Fsup_temp, Finf_temp, psi_root_sup_temp,   &
                        psi_root_inf_temp, psi_soil_sup_temp, psi_soil_inf_temp)

                ENDIF


                !! 4. Link with the hydrol.f90 module:
 
                !! The outputs of the previous subroutine call are the root water potential in each layer and the absorption fluxes in each layers.
                !! Before using the root water potentials and solve the hydraulic architecture, we have to make the link with the hydrol.f90 module
                !! by recalculating the root absorption but not inside our 2 layers but in the nslm layers of the hydrol.f90 scheme.
                
  		!! To link the absorption fluxes determined in the two previously defined layers with the nslm layers of hydrol.f90,
                !! a weighting of those two fluxes is made in order to dispatch the flux inside the layers.
                !! The weighting is based on the amount of water present in each layer compared to the mean one. 
                !! This weighting is purely arbitrary, we consider here that even if there is a few roots, the trees can absorb
                !! water from the layer. The main driving variable is the amount of water available in the layer.
                !! As the absorption is sometimes divided into nsub_step, the "sub-fluxes" are added together in order to get the
                !! total one at the end of the loop.
                !! e_frac is calculated after the nsub_step loop.

                IF (test_functional_roots_tuzet) THEN                
        
                   DO isl = 1 , nslm  !! Loop over soil layers
                      
                      IF (isl .LE. lim_layer) THEN
                         ! Will test this functional root profile based on wilting point
                         Fi_temp(ipts,ivm,isl,jst) = Fsup_temp(ipts,ivm) &
                              * (mc_out_temp(ipts,isl,jst)-mcw(ipts))*dh(isl)          &
                              /(MAX(0.001,(mc_sup_temp(ipts,jst)-mcw(ipts)))*SUM(dh(:lim_layer)))
                      ELSE

                         Fi_temp(ipts,ivm,isl,jst) = Finf_temp(ipts,ivm) &
                              * (mc_out_temp(ipts,isl,jst)-mcw(ipts))*dh(isl)          &
                              /(MAX(0.001,(mc_inf_temp(ipts,jst)-mcw(ipts)))*SUM(dh(lim_layer+1:))) 
                      ENDIF

                   ENDDO

                ELSE 

                   DO isl = 1 , nslm  !! Loop over soil layers
                      
                      IF (isl .LE. lim_layer) THEN
                         Fi_temp(ipts,ivm,isl,jst) = Fsup_temp(ipts,ivm) &
                              * (mc_out_temp(ipts,isl,jst)-mcr_sup(ipts))*dh(isl)          &
                              /(MAX(0.001,(mc_sup_temp(ipts,jst)-mcr_sup(ipts)))*SUM(dh(:lim_layer)))
                      ELSE

                         Fi_temp(ipts,ivm,isl,jst) = Finf_temp(ipts,ivm) &
                              * (mc_out_temp(ipts,isl,jst)-mcr_inf(ipts))*dh(isl)          &
                              /(MAX(0.001,(mc_inf_temp(ipts,jst)-mcr_inf(ipts)))*SUM(dh(lim_layer+1:))) 
                      ENDIF

                   ENDDO
                ENDIF

		!! 5. Calculation of the potentials
                !! determine the water potential at each level of the vegetation. Those water potential are calculated thanks to a 
                !! classical flux/potential function

                !! The previous collar water potential is saved for future estimation and the new one is calculated thanks to the 
                !! superficial root absorption and the superficial root water potential.
                IF (printlev>=4.AND.ipts==test_grid.AND.ivm==test_pft) THEN
                   WRITE(numout,*) 'Res_root_sup', Res_root_sup(ipts,ivm,icir)
                   WRITE(numout,*) 'psi_root_sup_temp', psi_root_sup_temp(ipts,ivm)
                   CALL flush(numout)
                ENDIF

                ! +++ CHECK +++
                ! Was gravity taken into account when calculating Fsup/inf_temp ?
                IF (.NOT. tuzet_inf_layer) THEN
                   psi_xylem_collar_temp(ipts,ivm) = MAX(psi_min,psi_root_sup_temp(ipts,ivm) &
                                                         - Fsup_temp(ipts,ivm) * Res_root_sup(ipts,ivm,icir) &
                                                         - rho_h2o*cte_grav*SUM(dlh(1:lim_layer))/kilo_to_unit)
                ELSE
                   psi_xylem_collar_temp(ipts,ivm) = MAX(psi_min,psi_root_inf_temp(ipts,ivm) &
                                                         - Finf_temp(ipts,ivm) * Res_root_inf(ipts,ivm,icir) &
                                                         - rho_h2o*cte_grav*SUM(dlh(lim_layer+1:))/kilo_to_unit)
                ENDIF 

                IF (printlev_loc>=4.AND.ipts==test_grid.AND.ivm==test_pft)  WRITE(numout,*) 'psi_xylem_collar_temp', psi_xylem_collar_temp(ipts,ivm)
                
                psi_leaf_temp_save(ipts,ivm) = psi_leaf_temp(ipts,ivm)
                
		IF (ok_hydrol_arch_storage) THEN
                   !! If the storage is enabled, all the hydraulic architecture has to be calculated.

                   !! Xylem water potentials at the trunk and leaf levels are saved and recalculated thanks to the corresponding fluxes,
                   !! resistances and water potentials 
                   ! See below for Res_mesophyll the explanations. In short,
                   ! resistances are sometimes set at maximum values, potential
                   ! need to remain below a minimum value
                   
                   grav_grad = rho_h2o*cte_grav * (circ_height(icir)-crown_width(icir))/kilo_to_unit

                   CALL resistances_curtail(ipts,ivm,icir,Fxylem_low,Res_xylem_low,psi_xylem_collar_temp,grav_grad)

                   ! +++++ CHECK +++++++++++++
                   ! Not using z_array_out which even if indexed with cc
                   ! is set for all cc on biggest circ class height
                   ! using crown width instead
                   psi_xylem_trunk_temp(ipts,ivm) = -Fxylem_low(ipts,ivm) * Res_xylem_low(ipts,ivm,icir) &
                      + psi_xylem_collar_temp(ipts,ivm) - grav_grad

                   ! Debug
                   IF ( (psi_xylem_trunk_temp(ipts,ivm).GT.zero) .OR. &
                           (psi_xylem_trunk_temp(ipts,ivm).LT.psi_min-min_sechiba) .AND. &
                           printlev_loc>=4 ) THEN 
                      WRITE(numout,*) 'WARNING : psi_xylem_trunk_temp not in [psi_min, 0] '
                      WRITE(numout,*) 'transpir :', transpir(ipts,ivm)
                      WRITE(numout,*) 'scale factor of transpir :', 1/(dt_sechiba * kilo_to_unit * soiltile(ipts,pref_soil_veg(ivm))&
                        * vegtot(ipts)) 
                      WRITE(numout,*) "this_tstep_calc(ipts,ivm,icir)", this_tstep_calc(ipts,ivm,icir)
                      WRITE(numout,*) 'transpir_circ :', transpir_circ(:)
                      WRITE(numout,*) 'transpir_test :', transpir_test(ipts,ivm)
                      IF (ok_hydrol_arch_storage) WRITE(numout,*) 'Fsto_wood :', Fsto_wood(ipts,ivm)
                      WRITE(numout,*) 'Fxylem_low :', Fxylem_low(ipts,ivm)
                      WRITE(numout,*) 'Res_xylem_low :', Res_xylem_low(ipts,ivm,icir)
                      WRITE(numout,*) 'psi_xylem_collar_temp :', psi_xylem_collar_temp(ipts,ivm)
                      WRITE(numout,*) 'psi_xylem_trunk_temp :', psi_xylem_trunk_temp(ipts,ivm)
                      WRITE(numout,*) 'psi_leaf_temp :', psi_leaf_temp(ipts,ivm)
                      IF ((soiltile(ipts,pref_soil_veg(ivm)) * vegtot(ipts)) .GT. min_sechiba) THEN
                         WRITE(numout,*) "transpir_test estimate factor", &
                            1/(dt_sechiba * kilo_to_unit* soiltile(ipts,pref_soil_veg(ivm)) * vegtot(ipts)) * & 
                           delta(ipts)*((1+(veget(ipts,ivm)/veget_max(ipts,ivm))* &
                           (tq_cdrag(ipts)*speed(ipts))/g_stom(ipts,ivm)) / &
                           (1+(veget(ipts,ivm)/veget_max(ipts,ivm))*(tq_cdrag(ipts)*speed(ipts) * &
                           f_psi_1(ipts,ivm))/(delta(ipts)*g_stom(ipts,ivm)*f_psi_2(ipts,ivm)))) 
                      ENDIF
                   ENDIF
                   !-

                   grav_grad = rho_h2o*cte_grav * (circ_height(icir)-crown_width(icir)/2)/kilo_to_unit

                   CALL resistances_curtail(ipts,ivm,icir,Fxylem_up,Res_xylem_up,psi_xylem_trunk_temp,grav_grad)

                   psi_xylem_leaf_temp(ipts,ivm) = -Fxylem_up(ipts,ivm) * Res_xylem_up(ipts,ivm,icir) + psi_xylem_trunk_temp(ipts,ivm) - grav_grad


                   !! Finally, the leaf water potential is calculated

                   IF ((soiltile(ipts,jst) * vegtot(ipts)) .GT. min_sechiba) THEN
                      ! This verification is to make sure that the resistance is
                      ! not too high compared to transpir test, as it occurs sometimes, leading to a
                      ! psi_leaf_temp that is way too low, which makes the model
                      ! crash. We set a minimum water potential and when the
                      ! resistance has been set to a maximum earlier, or if with
                      ! its current value it is going to lead to potentials that
                      ! are too low, it is here set such that the potential will be
                      ! equal to the minimum potential
                      grav_grad=zero
                      CALL resistances_curtail(ipts,ivm,icir,transpir_test,Res_mesophyll,psi_xylem_leaf_temp,grav_grad)

                      psi_leaf_temp(ipts,ivm) = -transpir_test(ipts,ivm) * Res_mesophyll(ipts,ivm,icir) + psi_xylem_leaf_temp(ipts,ivm) - grav_grad
                   ELSE
                      psi_leaf_temp(ipts,ivm) = psi_xylem_leaf_temp(ipts,ivm)
                   ENDIF


                ELSE !! If storage is disabled, we calculates only the new leaf water potential

                   IF (circ_class_biomass(ipts,ivm,icir,iroot,icarbon) .GE. min_sechiba) THEN !! If there is root biomass 

                      grav_grad = rho_h2o*cte_grav * (circ_height(icir)-crown_width(icir))/kilo_to_unit
                      CALL resistances_curtail(ipts,ivm,icir,Fxylem_low,Res_xylem_low,psi_xylem_collar_temp,grav_grad)
                      psi_xylem_trunk_temp(ipts,ivm) = -Fxylem_low(ipts,ivm) * Res_xylem_low(ipts,ivm,icir) + psi_xylem_collar_temp(ipts,ivm) - grav_grad 
        
                      grav_grad = rho_h2o*cte_grav* (circ_height(icir)-crown_width(icir)/2)/kilo_to_unit
                      CALL resistances_curtail(ipts,ivm,icir,Fxylem_up,Res_xylem_up,psi_xylem_trunk_temp,grav_grad)
                      psi_xylem_leaf_temp(ipts,ivm) = -Fxylem_up(ipts,ivm) * Res_xylem_up(ipts,ivm,icir) + psi_xylem_trunk_temp(ipts,ivm) - grav_grad

                      grav_grad = zero
                      CALL resistances_curtail(ipts,ivm,icir,transpir_test,Res_mesophyll,psi_xylem_leaf_temp,grav_grad)
                      psi_leaf_temp(ipts,ivm) = -transpir_test(ipts, ivm) * Res_mesophyll(ipts,ivm,icir) + psi_xylem_leaf_temp(ipts,ivm) - grav_grad

                   ELSE !! If there is no leaf biomass

                      psi_xylem_trunk_temp(ipts,ivm) = psi_xylem_collar_temp(ipts,ivm)

                      psi_xylem_leaf_temp(ipts,ivm) = psi_xylem_trunk_temp(ipts,ivm)

                      psi_leaf_temp(ipts,ivm) = psi_xylem_leaf_temp(ipts,ivm)

                   ENDIF

                ENDIF

                !! To avoid any problem in case of really low transpiration, potentials are set to zero when positive 
                !! (positive potential isn't physically possible)
                
                IF ( psi_root_sup_temp(ipts,ivm).GT.zero .OR. &
                psi_root_inf_temp(ipts,ivm).GT.zero .OR.      &
                psi_soil_sup_temp(ipts,jst).GT.zero .OR.      &
                psi_soil_inf_temp(ipts,jst).GT.zero .OR.      &
                psi_xylem_collar_temp(ipts,ivm).GT.zero .OR.  &
                psi_leaf_temp(ipts,ivm).GT.zero ) THEN
                   WRITE(numout,*) 'WARNING : root, soil, collar or leaf potential >0 !'
                   WRITE(numout,*) 'psi_xylem_collar_temp(ipts,ivm) =', psi_xylem_collar_temp(ipts,ivm)
                   WRITE(numout,*) 'psi_leaf_temp(ipts,ivm) =', psi_leaf_temp(ipts,ivm)
                   WRITE(numout,*) 'psi_root_inf_temp(ipts,ivm) =', psi_root_inf_temp(ipts,ivm)
                   WRITE(numout,*) 'psi_soil_sup_temp(ipts,ivm) =', psi_soil_sup_temp(ipts,jst)
                ENDIF

                psi_root_sup_temp(ipts,ivm) = MIN(psi_root_sup_temp(ipts,ivm), 0.)
                psi_root_inf_temp(ipts,ivm) = MIN(psi_root_inf_temp(ipts,ivm), 0.)
                psi_soil_sup_temp(ipts,jst) = MIN(psi_soil_sup_temp(ipts,jst), 0.)
                psi_soil_inf_temp(ipts,jst) = MIN(psi_soil_inf_temp(ipts,jst), 0.)
                psi_xylem_collar_temp(ipts,ivm) = MIN(psi_xylem_collar_temp(ipts,ivm), 0.)
                psi_leaf_temp(ipts,ivm) = MIN(psi_leaf_temp(ipts,ivm), 0.)



                !! As everyhting is now calculated, a distinction has to be done between the first iteration and the other ones.
                !! For the first iteration, all the values calculated are sent to the real values (the values of the current
                !! time step calculation). Then, if launch_next_calc = TRUE, the loop over the iterations is launched.
                !! If the current iteration is not the first in the loop, the new value of psi_leaf_temp is set and the loop
                !! continues.

                IF (this_tstep_calc(ipts,ivm,icir)) THEN

                   psi_leaf(ipts,ivm,icir) = psi_leaf_temp(ipts,ivm)
                   psi_xylem_leaf(ipts,ivm,icir) = psi_xylem_leaf_temp(ipts,ivm)
                   psi_xylem_trunk(ipts,ivm,icir) = psi_xylem_trunk_temp(ipts,ivm)
                   psi_xylem_collar(ipts,ivm,icir) = psi_xylem_collar_temp(ipts,ivm)
                   psi_sto_leaf(ipts,ivm,icir) = psi_sto_leaf_temp(ipts,ivm)
                   psi_sto_wood(ipts,ivm,icir) = psi_sto_wood_temp(ipts,ivm)
                   psi_root_sup(ipts,ivm,icir) = psi_root_sup_temp(ipts,ivm)
                   psi_root_inf(ipts,ivm,icir) = psi_root_inf_temp(ipts,ivm)
                   psi_soil_sup(ipts,jst) = psi_soil_sup_temp(ipts,jst)
                   psi_soil_inf(ipts,jst) = psi_soil_inf_temp(ipts,jst)
                   Finf(ipts,ivm) = Finf_temp(ipts,ivm)
                   Fsup(ipts,ivm) = Fsup_temp(ipts,ivm)
                   mc_sup(ipts,pref_soil_veg(ivm)) = mc_sup_temp(ipts,pref_soil_veg(ivm))
                   mc_inf(ipts,pref_soil_veg(ivm)) = mc_inf_temp(ipts,pref_soil_veg(ivm))
                   mc_i_sup(ipts,ivm,icir,:) = mc_i_sup_temp(ipts,ivm,:)
                   mc_i_inf(ipts,ivm,icir,:) = mc_i_inf_temp(ipts,ivm,:)
                   Fi(ipts,ivm,:,:) = Fi(ipts,ivm,:,:) + Fi_temp(ipts,ivm,:,:) !To add the fluxes from all the circ_classes
                   F_absorption(ipts,ivm) = F_absorption(ipts,ivm) + F_absorption_temp(ipts,ivm)
                
                   IF ((printlev_loc .GT. 4) .AND. ((mc_sup_temp(ipts,jst) - mcr_sup(ipts) .LT. zero) .OR. &
                      (mc_inf_temp(ipts,jst) - mcr_inf(ipts) .LT. zero))) THEN
                      WRITE(numout,*) 'WARNING: mc below residual 2'
                   ENDIF
                   
                   IF (ok_hydrol_arch_storage) THEN 
                      Fsto_leaf_save(ipts,ivm,icir) = Fsto_leaf(ipts,ivm)*soiltile(ipts,pref_soil_veg(ivm))*vegtot(ipts)*kilo_to_unit
                      Fsto_wood_save(ipts,ivm,icir) = Fsto_wood(ipts,ivm)*soiltile(ipts,pref_soil_veg(ivm))*vegtot(ipts)*kilo_to_unit

                      IF (.NOT. Vmax_leaf(ipts,ivm)-Vr_leaf(ipts,ivm) .LT. min_sechiba) THEN 
                         V_leaf_out(ipts,ivm,icir) = (V_leaf(ipts,ivm) -Vr_leaf(ipts,ivm))/(Vmax_leaf(ipts,ivm) - Vr_leaf(ipts,ivm))
                      ENDIF
                      IF (.NOT. Vmax_wood(ipts,ivm)-Vr_wood(ipts,ivm) .LT. min_sechiba) THEN 
                         V_wood_out(ipts,ivm,icir) = (V_wood(ipts,ivm) -Vr_wood(ipts,ivm))/(Vmax_wood(ipts,ivm) - Vr_wood(ipts,ivm))
                      ENDIF

                   ENDIF

                   IF (launch_next_calc(ipts,ivm,icir)) THEN

                      this_tstep_calc(ipts,ivm,icir) = .FALSE.
                      dummy_transpir(ipts,ivm,icir) =  transpir_test(ipts,ivm) * &
                        (dt_sechiba*kilo_to_unit*soiltile(ipts,pref_soil_veg(ivm))*vegtot(ipts))
                      CYCLE next_calc_loop
                   ELSE
                      dummy_transpir(ipts,ivm,icir) =  transpir_test(ipts,ivm) * &
                        (dt_sechiba*kilo_to_unit*soiltile(ipts,pref_soil_veg(ivm))*vegtot(ipts))
                   ENDIF


                ENDIF

                IF (.NOT. (this_tstep_calc(ipts,ivm,icir))) THEN

                   !! If the difference between 2 time steps is less than a threshold, the loop is exited.
                   IF (printlev .GT. 3) THEN
                      WRITE(numout,*) "This Time step values"
                      WRITE(numout,*) "psi_leaf_temp = ", psi_leaf_temp(ipts,ivm)
                      WRITE(numout,*) "transpir_test = ", transpir_test(ipts,ivm)
                      WRITE(numout,*) "transpir_fsupinf =", transpir(test_grid,test_pft)  / &
                           (dt_sechiba * kilo_to_unit* soiltile(test_grid,pref_soil_veg(test_pft)) * vegtot(test_grid) )
                      WRITE(numout,*) "Fxylem_up =", Fxylem_up(test_grid,test_pft)
                      WRITE(numout,*) "Fsto_leaf =", Fsto_leaf(test_grid,test_pft)
                      WRITE(numout,*) "Fxylem_low =", Fxylem_low(test_grid,test_pft)
                      WRITE(numout,*) "Fsto_wood =", Fsto_wood(test_grid,test_pft)
                      WRITE(numout,*) "F_absorption_temp =", F_absorption_temp(test_grid,test_pft)
                      WRITE(numout,*) "Fsup_temp =", Fsup_temp(test_grid,test_pft)
                      WRITE(numout,*) "Finf_temp =", Finf_temp(test_grid,test_pft)
                      WRITE(numout,*) "psi_root_sup_temp =", psi_root_sup_temp(test_grid,test_pft)
                      WRITE(numout,*) "psi_root_inf_temp =", psi_root_inf_temp(test_grid,test_pft)
                      WRITE(numout,*) "psi_soil_sup_temp =", psi_soil_sup_temp(test_grid,pref_soil_veg(test_pft))
                      WRITE(numout,*) "psi_soil_inf_temp =", psi_soil_inf_temp(test_grid,pref_soil_veg(test_pft))

                   ENDIF

                   IF (ABS(psi_leaf_temp(ipts,ivm) - psi_leaf_temp_save(ipts,ivm)) .GT. 0.005) THEN 


                      psi_leaf_temp(ipts,ivm) = (psi_leaf_temp_save(ipts,ivm) + psi_leaf_temp(ipts,ivm))/2.
                      dummy_transpir(ipts,ivm,icir) =  transpir_test(ipts,ivm) * &
                        (dt_sechiba*kilo_to_unit*soiltile(ipts,pref_soil_veg(ivm))*vegtot(ipts))

                   ELSE
                      
                      dummy_transpir(ipts,ivm,icir) =  transpir_test(ipts,ivm) * &
                        (dt_sechiba*kilo_to_unit*soiltile(ipts,pref_soil_veg(ivm))*vegtot(ipts))

                      EXIT next_calc_loop

                   ENDIF
                ENDIF
 
                dummy_transpir(ipts,ivm,icir) =  transpir_test(ipts,ivm) * &
                  (dt_sechiba*kilo_to_unit*soiltile(ipts,pref_soil_veg(ivm))*vegtot(ipts))

                IF (printlev_loc .GT. 3) WRITE(numout,*) 'tuzet iteration nb ', i

             ENDDO next_calc_loop

             IF (circ_class_biomass(ipts,ivm,icir,ileaf,icarbon)*circ_class_n(ipts,ivm,icir) .GT. min_sechiba) THEN  
                psi_leaf_mean_gs(ipts,ivm,icir) = psi_leaf(ipts,ivm,icir)
             ELSE
                psi_leaf_mean_gs(ipts,ivm,icir) = 0.
             ENDIF
            
             ! In summer astronomical midday is rather around 13:00 - CHECK - to
             ! refine
             IF (MOD(kjit,INT(one_day/dt_sechiba)).EQ.(13+1)*2) THEN
                psi_leaf_midday(ipts,ivm,icir) = psi_leaf(ipts,ivm,icir)
             ELSE
                psi_leaf_midday(ipts,ivm,icir) = 0.
             ENDIF

             !+++CHECK+++
             ! LGC: CHECK UNITS !!
             ! LGC: PSI or PSI_TEMP ?? 
             ! transpir .GT. min_sechiba ? Trop restrictif non ? A comparer aux
             ! valeurs usuelles
             
             IF (tuzet_inf_layer) THEN
                delta_psi(ipts,ivm,icir) = psi_leaf(ipts,ivm,icir)-psi_soil_inf(ipts,jst)
             ELSE
                delta_psi(ipts,ivm,icir) = psi_leaf(ipts,ivm,icir)-psi_soil_sup(ipts,jst)
             ENDIF

             IF (transpir_test(ipts,ivm) .GT. min_sechiba) THEN  
                IF (ok_tuzet_indiv) THEN
                   K_plant(ipts,ivm,icir) = delta_psi(ipts,ivm,icir)/transpir_test(ipts,ivm)
                ELSE 
                   K_plant(ipts,ivm,icir) = delta_psi(ipts,ivm,icir)/transpir_test(ipts,ivm)*circ_class_n(ipts,ivm,icir)
                ENDIF
             ELSE
                K_plant(ipts,ivm,icir) = min_sechiba*min_sechiba
             ENDIF
                
             ! Calculation of an equivalent resistance of the root system (to output)
             IF (Res_root_inf(ipts,ivm,icir).GT.min_sechiba .AND. Res_root_sup(ipts,ivm,icir).GT.min_sechiba) THEN
                R_root(ipts,ivm,icir) = 1 / (1/Res_root_inf(ipts,ivm,icir) + 1/Res_root_sup(ipts,ivm,icir))
             ELSE
                R_root(ipts,ivm,icir) = min_sechiba
             ENDIF

          ENDDO !ncirc
          
          !! The value of the leaf water potential sent to diffuco.f90 is either psi_leaf (when there is only one iteration)
          !! or the last value of psi_leaf_temp at the end of the iterative process.

          speed(ipts)=MIN(min_wind,SQRT(u(ipts)*u(ipts) + v(ipts)*v(ipts)))
          
          ! The inversion of the transpiration estimation equation to retrieve
          ! the leaf water potential is as below
          ! psi_leaf_next(ipts,ivm) = MAX(psi_min, psi_ref_g0(ivm) -&
          ! (LOG(((un - frac_evap_flood(ipts)) * (un - frac_sub_snow(ipts)) * frac_evap_transp(ipts,ivm) *&
          ! (rau(ipts)*dt_sechiba*delta(ipts)*(qsol_sat_new(ipts) - qair(ipts)))*&
          ! (dt_sechiba*kilo_to_unit*soiltile(ipts,pref_soil_veg(ivm))*vegtot(ipts))/SUM(dummy_transpir(ipts,ivm,:))-&
          ! un/(speed(ipts)*tq_cdrag(ipts)))*&
          ! delta(ipts)*(MAX(g0(ivm),(assimtot(ipts,ivm)+Rdtot_pft(ipts,ivm))/(cimean(ipts,ivm)-gamma_star(ipts))))*&
          ! (1+EXP(sf(ivm)*psi_ref_g0(iivm))) - un)/sf(ivm)))

          ! +++ CHECK +++
          ! If cimean-gamma_star<0 OR qair > qsat, do we want leaf water
          ! potential set at the previous organ's value ? 
          IF ((SUM(dummy_transpir(ipts,ivm,:)) .GT. zero) .AND. (cimean(ipts,ivm)-gamma_star(ipts).GT. zero)) THEN
             stress_log = ((un - frac_evap_flood(ipts)) * (un - frac_sub_snow(ipts)) * frac_evap_transp(ipts,ivm) * &
                     (rau(ipts)*dt_sechiba*delta(ipts)*(qsol_sat_new(ipts) - qair(ipts))) * &
                     (dt_sechiba*kilo_to_unit*soiltile(ipts,pref_soil_veg(ivm))*vegtot(ipts))/SUM(dummy_transpir(ipts,ivm,:)) - &
                     un/(speed(ipts)*tq_cdrag(ipts))) * &
                     delta(ipts)*(MAX(g0(ivm),(assimtot(ipts,ivm)+Rdtot_pft(ipts,ivm))/(cimean(ipts,ivm)-gamma_star(ipts)))) * &
                     (1+EXP(sf(ivm)*psi_ref_g0(ivm))) - un
 
             ! Quickfix if qair > qsat -> or other cases where stress_log <0 ? 
             IF (stress_log.GE.1/undef_sechiba)  THEN
                psi_leaf_next(ipts,ivm) = MAX(psi_min, psi_ref_g0(ivm) - LOG(stress_log)/sf(ivm))
             ELSE
                psi_leaf_next(ipts,ivm) = MINVAL(psi_xylem_leaf(ipts,ivm,:))
                IF (stress_log .LT. zero) WRITE(numout,*) 'WARNING : stress_log negative, qsat-qair=', qsol_sat_new(ipts) - qair(ipts) 
             ENDIF
          ELSE
            psi_leaf_next(ipts,ivm) = MINVAL(psi_xylem_leaf(ipts,ivm,:)) 
          ENDIF
          ! ++++++++++++

          IF (printlev_loc.GE.3) THEN
             WRITE(numout,*) 'dummy_transpir(ipts,ivm,:) = ', dummy_transpir(ipts,ivm,:) 
             WRITE(numout,*) 'transpir_test(ipts,ivm) = ', transpir_test(ipts,ivm) 
             WRITE(numout,*) 'transpir_circ(:) = ', transpir_circ(:) 
             WRITE(numout,*) 'transpir(ipts,ivm) = ', transpir(ipts,ivm) 
          ENDIF 

         ! psi_leaf_next(ipts,ivm) = psi_leaf_temp(ipts,ivm,icir)

          !! Finally, the link with the "rootsink" term in hydrol.f90 resolution is done thanks to the variable e_frac
          !! which represents the fraction of the F_absorption flux that is absorbed in each layer (same as the previous 
          !! hydraulic_arch.f90 module).           

          DO isl = 1, nslm

             IF (SUM(Fi(ipts,ivm,:,jst)) .GT. 0.) THEN

                e_frac(ipts,ivm,isl,jst) = Fi(ipts,ivm,isl,jst) / SUM(Fi(ipts,ivm,:,jst))

             ELSE 

                e_frac(ipts,ivm,isl,jst)=0.

             ENDIF

          ENDDO
       
       ENDDO  ! nvm

    ENDDO  ! kjpindex

    IF (printlev_loc .GT. 3) THEN
       counter=zero
       DO ipts=1,kjpindex
          DO jst=2,nstm
             IF ((mc_sup(ipts,jst) - mcr_sup(ipts) .LT. zero) .OR. &
                 (mc_inf(ipts,jst) - mcr_inf(ipts).LT. zero)) THEN
                counter=counter+1
             ENDIF
          ENDDO
       ENDDO
       IF (counter .NE. zero) WRITE(numout,*) 'WARNING : There are ', counter, &
          ' soiltilesXpixels with mc_sup/inf_temp below mcr'
    ENDIF

    IF (printlev .GT. 3) THEN
       WRITE(numout,*) "Next time step values"
       WRITE(numout,*) "dh =", dh(:)
       WRITE(numout,*) "SUM dh =", SUM(dh(:))
       WRITE(numout,*) "zmaxh =", zmaxh
       WRITE(numout,*) "veget_max", veget_max(:,:)
       WRITE(numout,*) "psi_sto_leaf_save =", psi_sto_leaf_save(test_grid,test_pft,:)
       WRITE(numout,*) "psi_xylem_trunk =", psi_xylem_trunk(test_grid,test_pft,:)
       WRITE(numout,*) "Res_sto_leaf =", Res_sto_leaf(test_grid,test_pft)
       WRITE(numout,*) "Res_xylem_up =", Res_xylem_up(test_grid,test_pft,:)
       WRITE(numout,*) "transpir_fsupinf =", transpir(test_grid,test_pft)  / &
            (dt_sechiba * kilo_to_unit* soiltile(test_grid,pref_soil_veg(test_pft)) * vegtot(test_grid) )
       WRITE(numout,*) "psi_sto_wood_save =", psi_sto_wood_save(test_grid,test_pft,:)
       WRITE(numout,*) "psi_xylem_collar =", psi_xylem_collar(test_grid,test_pft,:)
       WRITE(numout,*) "Res_xylem_low =", Res_xylem_low(test_grid,test_pft,:)
       WRITE(numout,*) "Fxylem_up =", Fxylem_up(test_grid,test_pft)
       WRITE(numout,*) "test_grid =", test_grid
       WRITE(numout,*) "test_pft =", test_pft
       WRITE(numout,*) "transpir =", transpir(test_grid,test_pft)
       WRITE(numout,*) "Fsto_leaf =", Fsto_leaf(test_grid,test_pft)
       WRITE(numout,*) "Fxylem_low =", Fxylem_low(test_grid,test_pft)
       WRITE(numout,*) "Fsto_wood =", Fsto_wood(test_grid,test_pft)
       WRITE(numout,*) "F_absorption =", F_absorption(test_grid,test_pft)
       WRITE(numout,*) "mc_out =", mc_out(test_grid,:,pref_soil_veg(test_pft))
       WRITE(numout,*) "mc_sup =", mc_sup(test_grid,pref_soil_veg(test_pft))
       WRITE(numout,*) "mc_inf =", mc_inf(test_grid,pref_soil_veg(test_pft))
       WRITE(numout,*) "Fsup =", Fsup(test_grid,test_pft)
       WRITE(numout,*) "Finf =", Finf(test_grid,test_pft)
       WRITE(numout,*) "Fi =", Fi(test_grid,test_pft,:,pref_soil_veg(test_pft))
       WRITE(numout,*) "SUM(Fi) =", SUM(Fi(test_grid,test_pft,:,pref_soil_veg(test_pft)))
       WRITE(numout,*) "e_frac =", e_frac(test_grid,test_pft,:,pref_soil_veg(test_pft))
       WRITE(numout,*) "SUM(e_frac) =", SUM(e_frac(test_grid,test_pft,:,pref_soil_veg(test_pft)))
       WRITE(numout,*) "psi_root_sup =", psi_root_sup(test_grid,test_pft,:)
       WRITE(numout,*) "psi_root_inf =", psi_root_inf(test_grid,test_pft,:)
       WRITE(numout,*) "psi_soil_sup =", psi_soil_sup(test_grid,pref_soil_veg(test_pft))
       WRITE(numout,*) "psi_soil_inf =", psi_soil_inf(test_grid,pref_soil_veg(test_pft))
       WRITE(numout,*) "psi_xylem_collar =", psi_xylem_collar(test_grid,test_pft,:)
       WRITE(numout,*) "psi_xylem_trunk =", psi_xylem_trunk(test_grid,test_pft,:)
       WRITE(numout,*) "psi_xylem_leaf =", psi_xylem_leaf(test_grid,test_pft,:)
       WRITE(numout,*) "psi_sto_wood =", psi_sto_wood(test_grid,test_pft,:)
       WRITE(numout,*) "psi_sto_leaf =", psi_sto_leaf(test_grid,test_pft,:)
       WRITE(numout,*) "psi_leaf =", psi_leaf(test_grid,test_pft,:)
       WRITE(numout,*) "psi_leaf_next =", psi_leaf_next(test_grid,test_pft)
       WRITE(numout,*) "psi_leaf_midday =", psi_leaf_midday(test_grid,test_pft,:)
       WRITE(numout,*) "Vr_leaf =", Vr_leaf(test_grid,test_pft)
       WRITE(numout,*) "Vmax_leaf =", Vmax_leaf(test_grid,test_pft)
       WRITE(numout,*) "V_leaf =", V_leaf(test_grid,test_pft)
       WRITE(numout,*) "C_leaf =", C_leaf(test_grid,test_pft)
       WRITE(numout,*) "tau_leaf =", tau_leaf(test_grid,test_pft)
       WRITE(numout,*) "Vr_wood =", Vr_wood(test_grid,test_pft)
       WRITE(numout,*) "Vmax_wood =", Vmax_wood(test_grid,test_pft)
       WRITE(numout,*) "V_wood =", V_wood(test_grid,test_pft)
       WRITE(numout,*) "C_wood =", C_wood(test_grid,test_pft)
       WRITE(numout,*) "tau_wood =", tau_wood(test_grid,test_pft)
    ENDIF

    !! Send variables to xios
    
    CALL xios_orchidee_send_field("Fsup",Fsup*one_day*kilo_to_unit) 
    CALL xios_orchidee_send_field("Finf",Finf*one_day*kilo_to_unit) 
    CALL xios_orchidee_send_field("mc_sup",mc_sup) 
    CALL xios_orchidee_send_field("mc_inf",mc_inf) 
    CALL xios_orchidee_send_field("psi_leaf",psi_leaf) 
    CALL xios_orchidee_send_field("psi_leaf_mean_gs",psi_leaf_mean_gs) 
    CALL xios_orchidee_send_field("psi_leaf_min",psi_leaf_mean_gs) 
    CALL xios_orchidee_send_field("psi_leaf_midday",psi_leaf_midday)
    CALL xios_orchidee_send_field("psi_xylem_trunk",psi_xylem_trunk) 
    CALL xios_orchidee_send_field("psi_xylem_leaf",psi_xylem_leaf) 
    CALL xios_orchidee_send_field("psi_xylem_collar",psi_xylem_collar) 
    CALL xios_orchidee_send_field("psi_sto_wood",psi_sto_wood) 
    CALL xios_orchidee_send_field("psi_sto_leaf",psi_sto_leaf) 
    CALL xios_orchidee_send_field("psi_root_sup",psi_root_sup)
    CALL xios_orchidee_send_field("psi_root_inf",psi_root_inf) 
    CALL xios_orchidee_send_field("psi_soil_sup",psi_soil_sup)
    CALL xios_orchidee_send_field("psi_soil_inf",psi_soil_inf)
    CALL xios_orchidee_send_field("Fsto_leaf",-Fsto_leaf_save) 
    CALL xios_orchidee_send_field("Fsto_wood",-Fsto_wood_save) 
    CALL xios_orchidee_send_field("leaf_stor_norm",V_leaf_out)
    CALL xios_orchidee_send_field("wood_stor_norm",V_wood_out)
    CALL xios_orchidee_send_field("K_plant",K_plant)
    ! Converting resistances from MPA.s/m3 to MPA.s/kg so it compares to K_plant
    CALL xios_orchidee_send_field("R_xylem",(Res_xylem_up+Res_xylem_low)/kilo_to_unit)
    CALL xios_orchidee_send_field("R_leaf",Res_mesophyll/kilo_to_unit)
    CALL xios_orchidee_send_field("R_root",R_root/kilo_to_unit)

  END SUBROUTINE hydrol_hydraulic_arch_tuzet_calc





!! ================================================================================================================================
!! SUBROUTINE 	: hydrol_hydraulic_arch_tuzet_resist
!!
!>\BRIEF           Calculates root absorption from the soil thanks to a dynamic resistance function of root biomass and soil
!!                 hydraulic conductivity. Calculates the root water potential.
!!
!!\n DESCRIPTION : 
!!   
!!       1. The subroutine starts to calculate the soil to root resistance thanks to a formula taken from Bonan et al. (2014).
!!       2. Once the resistance is determined in each layer (superficial and inferior), the routine calculates the flux split
!!          between the superficial and inferior layers.
!!       3. Finally, the root water potential is calculated thanks to the soil water potential and the previously calculated flux.
!!      
!! RECENT CHANGE(S): Added by Julien Alléon (December 2022)
!!
!! MAIN OUTPUT VARIABLE(S): psi_root_sup, psi_root_inf, Fsup, Finf  
!!
!! REFERENCE(S)	: Bonan et al. (2014)
!!
!! FLOWCHART    : 
!! 
!! ================================================================================================================================


  SUBROUTINE hydrol_hydraulic_arch_tuzet_resist (kjpindex, igrid, ipft, icir, soiltile, njsc, ks,              &
       nvan, avan, F_abs, circ_class_biomass, circ_class_n, mc_sup_temp, mc_inf_temp,           &
       Res_root_sup, Res_root_inf, root_profile, Fsup_temp, Finf_temp, psi_root_sup_temp,       &
       psi_root_inf_temp, psi_soil_sup_temp, psi_soil_inf_temp)

    !! Variable declaration

    !! Input Variables

    INTEGER(i_std),                                   INTENT(in)      :: kjpindex             !! Domain size, terrestrial pixels only (unitless)
    INTEGER(i_std),                                   INTENT(in)      :: igrid                !! Index of the pixel considered on the grid (unitless)
    INTEGER(i_std),                                   INTENT(in)      :: ipft                 !! Index of the pft considered (unitless)
    INTEGER(i_std),                                   INTENT(in)      :: icir                 !! Index of the considered circumference class (unitless)
    REAL(r_std),        DIMENSION(:,:),               INTENT(in)      :: soiltile             !! Fraction of each soiltile within vegtot (0-1, unitless) 
    INTEGER(i_std),     DIMENSION(:),                 INTENT(in)      :: njsc                 !! Index of the dominant soil textural class (unitless)  
    REAL(r_std),        DIMENSION(:),                 INTENT(in)      :: ks                   !! Hydraulic conductivity at saturation (mm {-1})
    REAL(r_std),        DIMENSION(:),                 INTENT(in)      :: nvan                 !! Van Genuchten coeficients n (unitless)
    REAL(r_std),        DIMENSION(:),                 INTENT(in)      :: avan                 !! Van Genuchten coeficients a (mm-1})
    REAL(r_std),        DIMENSION(:,:),               INTENT(in)      :: F_abs                !! Root absorption flux (m^3/s)
    REAL(r_std),        DIMENSION(:,:,:,:,:),         INTENT(in)      :: circ_class_biomass   !! Biomass per individual in each circumference class (gC/ind)
    REAL(r_std),        DIMENSION(:,:,:),             INTENT(in)      :: circ_class_n         !! Number of individual per circumference class (ind/m^2)
    REAL(r_std),        DIMENSION(:,:),               INTENT(in)      :: mc_sup_temp          !! Temporary value of the soil water potential in the superficial soil layer (MPa)
    REAL(r_std),        DIMENSION(:,:),               INTENT(in)      :: mc_inf_temp          !! Temporary value of the soil water potential in the inferior soil layer (MPa)
    REAL(r_std),        DIMENSION(:,:,:),             INTENT(in)      :: Res_root_sup         !! Root resistance in the superficial soil layer (MPa.s/m^3)
    REAL(r_std),        DIMENSION(:,:,:),             INTENT(in)      :: Res_root_inf         !! Root resistance in the inferior soil layer (MPa.s/m^3)
    REAL(r_std),        DIMENSION(:,:,:,:),           INTENT(in)      :: root_profile         !! Normalized root mass/length fraction in each soil layer 


    !! Output Variables

    REAL(r_std),        DIMENSION(:,:),               INTENT(out)     :: Fsup_temp            !! Temporary value of the water flux inside the superficial soil level (m^3/s)
    REAL(r_std),        DIMENSION(:,:),               INTENT(out)     :: Finf_temp            !! Temporary value of the water flux inside the inferior soil level (m^3/s)
    REAL(r_std),        DIMENSION(:,:),               INTENT(out)     :: psi_root_sup_temp    !! Temporary value of the superficial roots Water Potential (MPa)
    REAL(r_std),        DIMENSION(:,:),               INTENT(out)     :: psi_root_inf_temp    !! Temporary value of the inferior roots Water Potential (MPa)
    REAL(r_std),        DIMENSION(kjpindex,nstm),     INTENT(out)     :: psi_soil_sup_temp    !! Superficial soil layer Water Potential (MPa)
    REAL(r_std),        DIMENSION(kjpindex,nstm),     INTENT(out)     :: psi_soil_inf_temp    !! Inferior soil layer Water Potential (MPa)

    !! Modified Variables


    !! Local Variables

    INTEGER(i_std)                                                    :: jst                  !! Indices (respectively: grid-cells, PFTs, pref_soil_veg, soiltiles) 
    INTEGER(i_std)                                                    :: jsl                  !! Soil layer index
    REAL(r_std),        DIMENSION(nslm+1)                             :: z_soil               !! Depth of each node in the soil column (m)
    REAL(r_std),        DIMENSION(nvm)                                :: rpc                  !! Scaling factor for the calculation of the fraction of root biomass (unitless)
    REAL(r_std),        DIMENSION(nvm)                                :: root_dens_sup        !! Fraction of root biomass in the superficial soil column (unitless)
    REAL(r_std),        DIMENSION(nvm)                                :: root_dens_inf        !! Fraction of root biomass in the inferior soil column (unitless)
    REAL(r_std)                                                       :: lr_sup               !! Total root length in the superficial soil column per unit of soil volumne (m/m^3)
    REAL(r_std)                                                       :: lr_inf               !! Total root length in the inferior soil column per unit of soil volumne (m/m^3)
    REAL(r_std)                                                       :: rs_sup               !! Half distance between roots in the superficial soil column (m)
    REAL(r_std)                                                       :: rs_inf               !! Half distance between roots in the inferior soil column (m)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: Res_soil_sup         !! Soil to root resistance in the superficial soil column (MPa.s/m^3)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: Res_soil_inf         !! Soil to root resistance in the inferior soil column (MPa.s/m^3)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: ksoil_sup            !! Soil hydraulic conductivity in the superficial soil layer (mm/d)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: ksoil_inf            !! Soil hydraulic conductivity in the inferior soil layer (mm/d)
    REAL(r_std)                                                       :: m                    !! VG coefficient (unitless)
    REAL(r_std)                                                       :: frac                 !! Fraction to calculate VG conductivity (unitless)
    LOGICAL                                                           :: is_root              !! to deal with index for water potential calculation

!! ==========================================================================================================================================================================



    !! The subroutine calculation relies on 4 steps:
    !!   - The calculation of the soil water potential and soil hydraulic conductivity;
    !!   - The calculation of the soil root resistance thanks to the soil hydraulic conductivity and the root biomass;
    !!   - The calculation of the flux split between the superficial and inferior layers;
    !!   - The calculation of the root water potentials.
       
    !! Initialisations
    jst = pref_soil_veg(ipft)  !! Soil type corresponding to the PFT
    
    ksoil_sup(igrid,ipft) = 0.             
    ksoil_inf(igrid,ipft) = 0.    
    psi_soil_sup_temp(:,:) = 0. 
    psi_soil_inf_temp(:,:) = 0.

    !! As the soil water content in both layer is defined, the soil water potential can be calculated thanks to the
    !! Van Genuchten relationship (Van Genuchten (1980)). 
    IF ((printlev_loc .GT. 3).AND.(igrid==test_grid).AND. &
       (mc_sup_temp(igrid,jst)-mcr_sup(igrid).LT.zero)) THEN
       WRITE(numout,*) 'mc_sup-mcr_sup no muff < 0' 
       CALL flush(numout)
    ENDIF
       
    is_root = .FALSE.
    CALL calc_soil_potentials(igrid, ipft, jst, njsc, avan, nvan, mc_sup_temp, mc_inf_temp,&
                              psi_soil_sup_temp, psi_soil_inf_temp, is_root)

    ! Soil water debugging
    IF (printlev_loc .GT. 3) THEN 
       WRITE(numout,*), "jst : ", jst
       WRITE(numout,*), "psi_soil_sup_temp(testgrid,pref(ipft)) : ", psi_soil_sup_temp(test_grid,jst)
       WRITE(numout,*), "psi_soil_inf_temp(testgrid,pref(ipft)) : ", psi_soil_inf_temp(test_grid,jst)
       WRITE(numout,*), "mc_inf_temp(testgrid,pref(ipft)) : ", mc_inf_temp(test_grid,jst)
       WRITE(numout,*), "mc_sup_temp(testgrid,pref(ipft)) : ", mc_sup_temp(test_grid,jst)
       CALL flush(numout)
    ENDIF

    m=(1.-1./nvan(igrid))

    frac=MAX(0.001,MIN(0.999,(mc_sup_temp(igrid,jst)-mcr_sup(igrid))/(mcs_sup(igrid)-mcr_sup(igrid))))
    ksoil_sup(igrid,jst) = ks(igrid) * (frac**0.5) * ( un - ( un - frac ** (un/m)) ** m )**2

    frac=MAX(0.001,MIN(0.999,(mc_inf_temp(igrid,jst)-mcr_inf(igrid))/(mcs_inf(igrid)-mcr_inf(igrid))))
    ksoil_inf(igrid,jst) = ks(igrid) * (frac**0.5) * ( un - ( un - frac ** (un/m)) ** m )**2

    !! Because the soil water content cannot be lower than mcr, the soil water potential admit a minimum values. 
    !! The values below -5MPa are set to this lower limit

    psi_soil_sup_temp(igrid,jst) = MAX(psi_soil_sup_temp(igrid,jst),psi_root_min(ipft))
    psi_soil_inf_temp(igrid,jst) = MAX(psi_soil_inf_temp(igrid,jst),psi_root_min(ipft))

    !! The first calculation consists in determining the root biomass used in the calculation of the soil to root resistance.

    z_soil(1)= 0.
    z_soil(2:nslm+1) = zlt(1:nslm)
    root_dens_sup(ipft) = zero
    root_dens_inf(ipft) = zero

    !! We define a scaling factor in order to reduce the fraction of root biomass in each layer ("root_dens") to a unitless
    !! factor between 0 and 1.

    IF (circ_class_biomass(igrid,ipft,icir,iroot,icarbon) .GE. min_sechiba) THEN 
       !! There is no calculation when there is no roots

       !! Caluclation of the fraction of root biomass in each soil layer. The calculation relies on the exponential decrease of
       !! root density in the soil column. The difference of exponentials corresponds to the fraction of biomass between the 
       !! two corresponding layer (layer 1 to lim_layer for the superficial soil layer and layer lim_layer+1 to nslm for the 
       !! inferior one). The difference of exponentials is then weigthed by the rpc factor which permits to calculate the fraction
       !! of total biomass in each layer.

       ! +++ CHECK +++
       ! CHECK other place in the file where this should be duplicated 
       ! At the first time step, as root_profile has not been calculated, we
       ! calculate a strucural root profile (which does not need the soil water
       ! variables calculated only late, in hydrol_soil)
       ! +++++++++++++ 
       IF ( SUM(root_profile(igrid,ipft,:,ifunc)).EQ.zero .OR. test_functional_roots_tuzet) THEN
          !! We define a scaling factor in order to reduce the fraction of root biomass in each layer ("root_dens") to a unitless
          !! factor between 0 and 1.
          rpc(ipft) = un / (un - EXP(-z_soil(nslm+1) * humcste(ipft)))
          root_dens_sup(ipft) = rpc(ipft) * (EXP(-z_soil(1)*humcste(ipft)) - &
               EXP( -z_soil(lim_layer+1)*humcste(ipft)))

          root_dens_inf(ipft) = rpc(ipft) * (EXP(-z_soil(lim_layer + 1)*humcste(ipft)) - &
               EXP( -z_soil(nslm+1)*humcste(ipft)))
       ELSE
       ! Weighted mean of functional root profile in sup and inf layers
          DO jsl = 1, lim_layer  
             root_dens_sup(ipft) = root_dens_sup(ipft) + root_profile(igrid,ipft,jsl,ifunc) * dlh(jsl)/SUM(dlh(1:lim_layer)) ! replace
          ENDDO
          DO jsl = lim_layer+1, nslm
             root_dens_inf(ipft) = root_dens_inf(ipft) + root_profile(igrid,ipft,jsl,ifunc) * dlh(jsl)/SUM(dlh(lim_layer+1:))
          ENDDO
       ENDIF

       !! The root length is calculated thanks to the fraction of root biomass (unitless) multipled by the root biomass (g/m^2)
       !! and a factor which links the root biomass to the root length (srl: PFT parameter). The total is then divided by the 
       !! total soil height in order to have a length of root per unit volume of soil.

       lr_sup = root_dens_sup(ipft)*deux*circ_class_biomass(igrid,ipft,icir,iroot,icarbon)*circ_class_n(igrid,ipft,icir)*srl(ipft)&
            /SUM(dlh(:lim_layer))

       !! The half distance between roots and the soil to root resistance are calculated thanks to the root length and the soil hydraulic
       !! conductivity(Bonan et al. (2014) Eq.A22). 

       rs_sup = (un/(pi*lr_sup))**(undemi)

       Res_soil_sup(igrid,ipft)=(LOG(rs_sup/r_froot(ipft)))/(deux*pi*lr_sup*SUM(dlh(:lim_layer))&
            *ksoil_sup(igrid,jst)*mega_to_unit/(cte_grav*ph2o*mm_m * one_day))

       !! Same calculation for the inferior layer.

       lr_inf = root_dens_inf(ipft)*deux*circ_class_biomass(igrid,ipft,icir,iroot,icarbon)*circ_class_n(igrid,ipft,icir)*srl(ipft)&
            /SUM(dlh(lim_layer+1:))

       rs_inf = (un/(pi*lr_inf))**(undemi)

       Res_soil_inf(igrid,ipft)=(LOG(rs_inf/r_froot(ipft)))/(deux*pi*lr_inf*SUM(dlh(lim_layer+1:))&
            *ksoil_inf(igrid,jst)*mega_to_unit/(cte_grav*ph2o*mm_m * one_day))


       !! The next step of the resolution is to determine the absorption fluxes at the soil/root interface. Those fluxes are 
       !! determined thanks to the transpiration at the previous timestep. The calculation results from a Kirchhoff's current law.
       !! The fluxes are determined thanks to the difference of water potentials at the previous time-step. 


       IF (F_abs(igrid,ipft) .GT. 0.) THEN

          IF ((soiltile(igrid,jst) * vegtot(igrid)) .GT. min_sechiba) THEN

             Finf_temp(igrid,ipft) =  MIN(F_abs(igrid,ipft), &
                  MAX(0.,(F_abs(igrid,ipft)&
                  + (psi_soil_inf_temp(igrid,jst) - psi_soil_sup_temp(igrid,jst)) / (Res_soil_sup(igrid,ipft) + Res_root_sup(igrid,ipft,icir))) / &
                  ( 1. + (Res_soil_inf(igrid,ipft) + Res_root_inf(igrid,ipft,icir)) / (Res_soil_sup(igrid,ipft) + Res_root_sup(igrid,ipft,icir)))))


             Fsup_temp(igrid,ipft) = F_abs(igrid,ipft) - Finf_temp(igrid,ipft)

          ENDIF

       ELSE

          ! To avoid problems in case of null transpiration, the fluxes are set to zero when transpiration is null

          Fsup_temp(igrid,ipft) = 0.
          Finf_temp(igrid,ipft) = 0.

       ENDIF

    ELSE ! If there is no root biomass

       Res_soil_sup(igrid,ipft) = 0.
       Res_soil_inf(igrid,ipft) = 0.
       Fsup_temp(igrid,ipft) = 0.
       Finf_temp(igrid,ipft) = 0.

    ENDIF

    !! Finally, the roots water potentials are calculated thanks to the soil water potentials and the fluxes previously calculated.

    psi_root_sup_temp(igrid,ipft) = -Fsup_temp(igrid,ipft) * Res_soil_sup(igrid,ipft) + psi_soil_sup_temp(igrid,jst) 
    psi_root_inf_temp(igrid,ipft) = -Finf_temp(igrid,ipft) * Res_soil_inf(igrid,ipft) + psi_soil_inf_temp(igrid,jst) 

    IF ( printlev_loc >= 4 ) THEN
       WRITE(numout,*) "psi_soil_sup =", psi_soil_sup_temp(igrid,jst)
       WRITE(numout,*) "psi_soil_inf =", psi_soil_inf_temp(igrid,jst)
    ENDIF
    
  END SUBROUTINE hydrol_hydraulic_arch_tuzet_resist




!! ================================================================================================================================
!! SUBROUTINE 	: hydrol_hydraulic_arch_tuzet_muff
!!
!>\BRIEF        Calculates the root absorption with a radial resolution of Richard's equation around a fictive root. Subroutine 
!!              called thanks to the flag "ok_hydrol_arch_muff".   
!!
!!\n DESCRIPTION : 
!!   
!!              This subroutine calculates the root absorption by resolving Richard's equation radially around a fictive root.
!!              The steps of resolution are the following:
!!                 1. Calculation of the root length;
!!                 2. Calculation of the muff raidus;
!!                 3. Resolution of Richard's equation;
!!                 4. Calculation of the root water potentials.
!!      
!! RECENT CHANGE(S): Added by Julien Alléon (December 2022)
!!
!! MAIN OUTPUT VARIABLE(S): Fsup, Finf, psi_root_sup, psi_root_inf, mc_i_sup, mc_i_inf.  
!!
!! REFERENCE(S)	: Tuzet et al. 2017
!!                Tuzet et al. 2003
!!
!! FLOWCHART    : 
!! 
!! ================================================================================================================================


  SUBROUTINE hydrol_hydraulic_arch_tuzet_muff (kjit, kjpindex, igrid, ipft, icir, soiltile, veget_max, njsc, ks,    &
       nvan, avan, F_abs, circ_class_biomass, circ_class_n, Res_root_sup, Res_root_inf, root_profile, Fsup_temp, Finf_temp, &
       psi_soil_sup_temp, psi_soil_inf_temp, psi_root_sup_temp, psi_root_inf_temp, mc_sup_temp, mc_inf_temp, mc_i_sup_temp, mc_i_inf_temp)

!! Variable declaration

    !! Input Variables

    INTEGER(i_std),                                   INTENT(in)      :: kjit                 !! Time step number (unitless)
    INTEGER(i_std),                                   INTENT(in)      :: kjpindex             !! Domain size, terrestrial pixels only (unitless)
    INTEGER(i_std),                                   INTENT(in)      :: igrid                !! Index of the pixel considered on the grid (unitless)
    INTEGER(i_std),                                   INTENT(in)      :: ipft                 !! Index of the pft considered (unitless)
    INTEGER(i_std),                                   INTENT(in)      :: icir                !! Index of the considered considered circumference class (unitless)
    REAL(r_std),        DIMENSION(:,:),               INTENT(in)      :: soiltile             !! Fraction of each soiltile within vegtot (0-1, unitless)  
    REAL(r_std),        DIMENSION(:,:),               INTENT(in)      :: veget_max            !! Maximum fraction of vegetation type (unitless) 
    INTEGER(i_std),     DIMENSION(:),                 INTENT(in)      :: njsc                 !! Index of the dominant soil textural class (unitless)
    REAL(r_std),        DIMENSION(:),                 INTENT(in)      :: ks                   !! Hydraulic conductivity at saturation (mm {-1})
    REAL(r_std),        DIMENSION(:),                 INTENT(in)      :: nvan                 !! Van Genuchten coeficients n (unitless)
    REAL(r_std),        DIMENSION(:),                 INTENT(in)      :: avan                 !! Van Genuchten coeficients a (mm-1})
    REAL(r_std),        DIMENSION(:,:),               INTENT(in)      :: F_abs                !! Root absorption flux (m^3/s)
    REAL(r_std),        DIMENSION(:,:,:,:,:),         INTENT(in)      :: circ_class_biomass   !! Biomass per individual for each circumference class (g/ind)
    REAL(r_std),        DIMENSION(:,:,:),             INTENT(in)      :: circ_class_n         !! Number of individuals per circumference class (ind/m^2)
    REAL(r_std),        DIMENSION(:,:,:),             INTENT(in)      :: Res_root_sup         !! Root resistance in the superficial soil layer (MPa.s/m^3)
    REAL(r_std),        DIMENSION(:,:,:),             INTENT(in)      :: Res_root_inf         !! Root resistance in the inferior soil layer (MPa.s/m^3)
    REAL(r_std),        DIMENSION(:,:,:,:),           INTENT(in)      :: root_profile         !! Normalized root mass/length fraction in each soil layer 

    !! Output Variables

    REAL(r_std),        DIMENSION(:,:),               INTENT(out)     :: Fsup_temp                 !! Temporary variable of the root absorption flux in the superficial soil layer (m^3/s)
    REAL(r_std),        DIMENSION(:,:),               INTENT(out)     :: Finf_temp                 !! Temporary variable of the root absorption flux in the superficial soil layer (m^3/s)

    !! Modified Variables

    !! Water potentials
    REAL(r_std),        DIMENSION(:,:),               INTENT(inout)   :: psi_soil_sup_temp         !! Temporary variable of the superficial soil Water Potential (MPa)
    REAL(r_std),        DIMENSION(:,:),               INTENT(inout)   :: psi_soil_inf_temp         !! Temporary variable of the inferior soil Water Potential (MPa)
    REAL(r_std),        DIMENSION(:,:),               INTENT(inout)   :: psi_root_sup_temp         !! Temporary variable of the superficial roots Water Potential (MPa)
    REAL(r_std),        DIMENSION(:,:),               INTENT(inout)   :: psi_root_inf_temp         !! Temporary variable of the inferior roots Water Potential (MPa)

    !! Water contents
    REAL(r_std),        DIMENSION(:,:),               INTENT(inout)   :: mc_sup_temp               !! Temporary variable of the water content of the superficial layer (m^3/m^3)
    REAL(r_std),        DIMENSION(:,:),               INTENT(inout)   :: mc_inf_temp               !! Temporary variable of the water content of the inferior layer (m^3/m^3)
    REAL(r_std),        DIMENSION(:,:,:),             INTENT(inout)   :: mc_i_sup_temp             !! Temporary variable of the discretized values of the soil water content inside the superficial cylinder (m^3/m^3)
    REAL(r_std),        DIMENSION(:,:,:),             INTENT(inout)   :: mc_i_inf_temp             !! Temporary variable of the discretized values of the soil water content inside the inferior cylinder (m^3/m^3)

    !! Local Variables

    !! Soil depth
    REAL(r_std),        DIMENSION(nslm+1)                             :: z_soil               !! Depth of each node in the soil column (m)

    !! Root properties
    REAL(r_std),        DIMENSION(nvm)                                :: rpc                  !! Scaling factor for the calculation of the fraction of root biomass (unitless)
    REAL(r_std),        DIMENSION(nvm)                                :: root_dens_sup        !! Fraction of root biomass in the superficial soil column (unitless)
    REAL(r_std),        DIMENSION(nvm)                                :: root_dens_inf        !! Fraction of root biomass in the inferior soil column (unitless)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: lr_muff_sup          !! Length of the fine roots inside the superficial soil layer (m/m^3)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: lr_muff_inf          !! Length of the fine roots inside the inferior soil layer (m/m^3) 
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: root_radius          !! Fine root radius in the muff
    REAL(r_std),        DIMENSION(ncirc)                              :: circ_root_share      !! Share of root mass per circ_class

    !! Soil tile fluxes
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: F_sup_st             !! Water flux inside the superficial soil level (m^3/s)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: F_inf_st             !! Water flux inside the inferior soil level (m^3/s)

    !! Radius
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: Rad_sup              !! Radius of the superficial muff (m)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: Rad_inf              !! Radius of the infrerior muff (m)
    REAL(r_std),        DIMENSION(kjpindex,nvm,nrp)                   :: R_i_sup              !! Radius discretization of the superficial muff (m)
    REAL(r_std),        DIMENSION(kjpindex,nvm,nrp)                   :: R_i_inf              !! Radius discretization of the infrerior muff (m)
    REAL(r_std),        DIMENSION(kjpindex,nvm,nrp)                   :: dri_sup              !! Thickness of each radial layer in the discretization of the superficial muff (m)
    REAL(r_std),        DIMENSION(kjpindex,nvm,nrp)                   :: dri_inf              !! Thickness of each radial layer in the discretization of the infrerior muff (m)
    
    !! Radial soil properties
    REAL(r_std),        DIMENSION(kjpindex,nvm,nrp)                   :: d_lin_rad_sup        !! Radial diffusivity inside the superficial muff (mm^2/d)
    REAL(r_std),        DIMENSION(kjpindex,nvm,nrp)                   :: d_lin_rad_inf        !! Radial diffusivity inside the inferior muff (mm^2/d)
    REAL(r_std),        DIMENSION(kjpindex,nvm,nrp)                   :: ksoil_rad_sup        !! Radial conductivity inside the superficial muff (mm/d)
    REAL(r_std),        DIMENSION(kjpindex,nvm,nrp)                   :: ksoil_rad_inf        !! Radial conductivity inside the inferior muff (mm/d)
    REAL(r_std)                                                       :: m                    !! VG coefficient  (unitless)
    REAL(r_std)                                                       :: frac                 !! Fraction to calculate VG conductivity (unitless)

    !! Water contents
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: delta_mc_sup         !! Difference ratio of the superficial soil water content between two time steps (unitless)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: delta_mc_inf         !! Difference ratio of the inferior soil water content between two time steps (unitless)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: mc_pft_sup           !! Water content at the superficial root level of each pft (m^3/m^3)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: mc_pft_inf           !! Water content at the inferior root level of each pft (m^3/m^3)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: mc_i_sup_mean        !! Water content at the superficial root level of each pft (m^3/m^3)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                       :: mc_i_inf_mean        !! Water content at the inferior root level of each pft (m^3/m^3)
    INTEGER(i_std)                                                    :: mc_ratio             !! Ratio used in the calculation of the diffusivity (unitless)

    !! Richard's equation resolution variables
    REAL(r_std),        DIMENSION(kjpindex,nrp,3)                     :: tmat_rad             !! Left hand matrix for the radial resolution
    REAL(r_std),        DIMENSION(kjpindex,nrp)                       :: rhs_rad              !! Right hand matrix for the radial resolution

    !! Indices
    INTEGER(i_std)                                                    :: i                    !! Indice used in the calculation of the diffusivity
    INTEGER(i_std)                                                    :: jrp                  !! Indice for discretization inside the muff
    INTEGER(i_std)                                                    :: ipts,ivm,isl,ist,jst !! Indices (respectively: grid-cells, PFTs, soil layers, soil tiles, pref_soil_veg)

    LOGICAL                                                           :: is_sup               !! Flag to check if the resolution occurs in the superficial or inferior soil layer
    LOGICAL                                                           :: is_root              !! to deal with index for water potential calculation

!! ==========================================================================================================================================================================
    
    !! Initializations
    
    lr_muff_inf(igrid,ipft) = 0.                
    lr_muff_sup(igrid,ipft) = 0.
    F_inf_st(igrid,ipft) = 0.                
    F_sup_st(igrid,ipft) = 0.
    delta_mc_sup(igrid,ipft) = 1.                
    delta_mc_inf(igrid,ipft) = 1.
    d_lin_rad_sup(igrid,ipft,:) = 0.                
    d_lin_rad_inf(igrid,ipft,:) = 0.
    ksoil_rad_sup(igrid,ipft,:) = 0.                
    ksoil_rad_inf(igrid,ipft,:) = 0.


    !! The subroutine relies on several steps:
    !!   1. The definition of the muffs
    !!   2. The radial resolution of Richard's Equation
    !!   3. The calculation of the root water potentials.


    !! 1. Definition of the muffs:

    !! The principle of the muffs is the following. An "absorption muff" is a cylinder of soil around a fictive root of length lr_sup or lr_inf. 
    !! The root has a radius (root_radius) and, around this root, the voulme of soil in each soil horizon is distributed in order to form a 
    !! cylinder. In this cylinder (the so called "muff"), of radius Rad_sup or Rad_inf, the radial resolution of Richard's equation will occur.
    !! To do so, the muff is discretized in nrp nodes on which we will determine the water content. The water gradient between all those nodes
    !! is saved at each time step (mc_i_sup and mc_i_inf). At the next time step, this saved gradient is multiplied by the ratio between the mean
    !! water content in each layer at the present and previosu time step. The idea is to conserve the gradient and remove or add water in the muff 
    !! according to the resolution of hydrol. The aim of the multiplication is to conserve water.


    !! Calculation of the root radius inside the muff at the soil tile level. Mean value of all the fine roots radiuses 
    !! according to the soil tile and the number of PFTs it regroups.

    !! DEBUG !!
    !! This part should be better coded as the number of PFTs differs from one configuration to the other.

    root_radius(igrid,ipft)=0.

    jst=pref_soil_veg(ipft)
    root_radius(igrid,ipft)=r_froot(ipft)

    !! Soil tile 1 corresponds to bare soil and has no roots inside it
    
    is_root = .FALSE.
    CALL calc_soil_potentials(igrid, ipft, jst, njsc, avan, nvan, mc_sup_temp, mc_inf_temp,&
                              psi_soil_sup_temp, psi_soil_inf_temp, is_root)


    !! Initialisation of the muff water content gradient at the first time step in order to avoid problems 

    IF (kjit==1) THEN
       DO jrp = 1, nrp
          mc_i_sup_temp(igrid,ipft,jrp)=mc_sup_temp(igrid,jst)
          mc_i_inf_temp(igrid,ipft,jrp)=mc_inf_temp(igrid,jst)
       ENDDO
    ENDIF


    z_soil(1)= 0.
    z_soil(2:nslm+1) = zlt(1:nslm)
    root_dens_sup(ipft) = zero
    root_dens_inf(ipft) = zero

    IF (circ_class_biomass(igrid,ipft,icir,iroot,icarbon)*circ_class_n(igrid,ipft,icir) .GE. min_sechiba) THEN 
       !! There is no calculation when there is no roots

       !! Calculation of the fraction of root biomass in each soil layer. The calculation relies on the exponential decrease of
       !! root density in the soil column. The difference of exponentials corresponds to the fraction of biomass between the 
       !! two corresponding layer (layer 1 to lim_layer for the superficial soil layer and layer lim_layer+1 to nslm for the 
       !! inferior one). The difference of exponentials is then weighted by the rpc factor which permits to calculate the fraction
       !! of total biomass in each layer.

       ! +++ CHECK +++
       ! CHECK other place in the file where this should be duplicated 
       ! At the first time step, as root_profile has not been calculated, we
       ! calculate a strucural root profile (which does not need the soil water
       ! variables calculated only late, in hydrol_soil)
       ! +++++++++++++ 
       IF ( SUM(root_profile(igrid,ipft,:,ifunc)).EQ.zero .OR. test_functional_roots_tuzet) THEN
          !! We define a scaling factor in order to reduce the fraction of root biomass in each layer ("root_dens") to a unitless
          !! factor between 0 and 1.
          rpc(ipft) = un / (un - EXP(-z_soil(nslm+1) * humcste(ipft)))
          root_dens_sup(ipft) = rpc(ipft) * (EXP(-z_soil(1)*humcste(ipft)) - &
               EXP( -z_soil(lim_layer+1)*humcste(ipft)))

          root_dens_inf(ipft) = rpc(ipft) * (EXP(-z_soil(lim_layer + 1)*humcste(ipft)) - &
               EXP( -z_soil(nslm+1)*humcste(ipft)))
       ELSE
       ! Weighted mean of functional root profile in sup and inf layers
          DO isl = 1, lim_layer  
             root_dens_sup(ipft) = root_dens_sup(ipft) + root_profile(igrid,ipft,isl,ifunc) * dlh(isl)/SUM(dlh(1:lim_layer)) 
          ENDDO
          DO isl = lim_layer+1, nslm
             root_dens_inf(ipft) = root_dens_inf(ipft) + root_profile(igrid,ipft,isl,ifunc) * dlh(isl)/SUM(dlh(lim_layer+1:))
          ENDDO
       ENDIF

       !! The root length is calculated thanks to the fraction of root biomass (unitless) multipled by the root biomass (g/m^2)
       !! and a factor which links the root biomass to the root length (srl: PFT parameter). The total is then divided by the 
       !! total soil height in order to have a length of root per unit volume of soil.

       circ_root_share(icir) = circ_class_biomass(igrid,ipft,icir,iroot,icarbon)*circ_class_n(igrid,ipft,icir)/&
            SUM(circ_class_biomass(igrid,ipft,:,iroot,icarbon)*circ_class_n(igrid,ipft,:))

       lr_muff_sup(igrid,ipft) = lr_muff_sup(igrid,ipft) + 2.*root_dens_sup(ipft)*&
            circ_class_biomass(igrid,ipft,icir,iroot,icarbon)*circ_class_n(igrid,ipft,icir)*srl(ipft)&
            /SUM(dlh(:lim_layer))


       lr_muff_inf(igrid,ipft) = lr_muff_inf(igrid,ipft) + 2.*root_dens_inf(ipft)*&
            circ_class_biomass(igrid,ipft,icir,iroot,icarbon)*circ_class_n(igrid,ipft,icir)*srl(ipft)&
            /SUM(dlh(lim_layer + 1:))

       !! The fluxes in the superficial and inferior soil layers are determined thanks to the resolution of a Kirchhoff's law 
       !! of currents. The fluxes are determined thanks to the difference of root water potentials in the two soil layers at 
       !! the previous time step.

       IF (F_abs(igrid,ipft) .GT. 0.) THEN
          IF ((soiltile(igrid,jst) * vegtot(igrid)) .GT. min_sechiba) THEN
             !! If ther is no vegetation, there is no flux
             IF (ipft .EQ. test_pft) THEN 
             ENDIF
             Finf_temp(igrid,ipft) =  MIN(F_abs(igrid,ipft), MAX(0.,(F_abs(igrid,ipft)&
                  + (psi_root_inf_temp(igrid,ipft) - psi_root_sup_temp(igrid,ipft)) / (Res_root_sup(igrid,ipft,icir))) / &
                  ( 1. + Res_root_inf(igrid,ipft,icir) / Res_root_sup(igrid,ipft,icir))))

             Fsup_temp(igrid,ipft) = F_abs(igrid,ipft) - Finf_temp(igrid,ipft)
          ENDIF
       ELSE
          ! To avoid problems in case of null transpiration, the fluxes are set to zero when transpiration is null
          Fsup_temp(igrid,ipft) = 0.
          Finf_temp(igrid,ipft) = 0.
       ENDIF
    ELSE ! If there is no root biomass
       Fsup_temp(igrid,ipft) = 0.
       Finf_temp(igrid,ipft) = 0.
    ENDIF

    !! Calculation of the muffs' radiuses. According to a discussion with Andrée Tuzet, the radius is calculated as follows thanks to 
    !! the fine roots length in each layer. We, then, substract the root radius in order to consider only a muff of soil.

    IF (lr_muff_sup(igrid,ipft) .GT. min_sechiba) THEN
       !! In order to consider only PFTs where there is actually vegetation
       Rad_sup(igrid,ipft)=(pi*lr_muff_sup(igrid,ipft)/circ_root_share(icir))**(-1./2.)-root_radius(igrid,ipft)
    ELSE
       Rad_sup(igrid,ipft)=0.
    ENDIF


    IF (lr_muff_inf(igrid,ipft) .GT. min_sechiba) THEN
       !! In order to consider only PFTs where there is actually vegetation
       Rad_inf(igrid,ipft)=(pi*lr_muff_inf(igrid,ipft)/circ_root_share(icir))**(-1./2.)-root_radius(igrid,ipft)
    ELSE
       Rad_inf(igrid,ipft)=0.
    ENDIF


    !! Discretisation of the muffs once the radius is known.

    DO jrp = 1, nrp ! Loop over muff nodes
       !! The muff is discretised in nrp nodes following an order 3 equation
       !! The aim is to have more nodes near the root (jrp=0), where the water
       !! gradient is strong.
       R_i_sup(igrid,ipft,jrp) = Rad_sup(igrid,ipft) * (REAL(jrp)/REAL(nrp))**3         
       R_i_inf(igrid,ipft,jrp) = Rad_inf(igrid,ipft) * (REAL(jrp)/REAL(nrp))**3        
    ENDDO

    !! Distance between two nodes in the muff discretisation

    dri_sup(igrid,ipft,1) = 0.                    
    dri_inf(igrid,ipft,1) = 0.                  

    DO jrp = 2, nrp
       dri_sup(igrid,ipft,jrp) = (R_i_sup(igrid,ipft,jrp)-R_i_sup(igrid,ipft,jrp-1))*kilo_to_unit           
       dri_inf(igrid,ipft,jrp) = (R_i_inf(igrid,ipft,jrp)-R_i_inf(igrid,ipft,jrp-1))*kilo_to_unit       
    ENDDO

    !! We apply here the difference in mean water content from one time step to another by multiplying the saved
    !! gradient by a ratio between the mean water content in each layer from the previous time step to the new one.

    !! Caluclation of the mean water content inside the muff weighted by the distance between each node.
    
    !! Superficial muff :
    
    IF (SUM(dri_sup(igrid,ipft,:)) .GT. min_sechiba) THEN
       mc_i_sup_mean(igrid,ipft)=SUM(mc_i_sup_temp(igrid,ipft,:)*dri_sup(igrid,ipft,:))/SUM(dri_sup(igrid,ipft,:))
    ELSE
       mc_i_sup_mean(igrid,ipft)  = mc_sup_temp(igrid,jst)
    ENDIF

    !! Inferior muff :

    IF (SUM(dri_inf(igrid,ipft,:)) .GT. min_sechiba) THEN
       mc_i_inf_mean(igrid,ipft)=SUM(mc_i_inf_temp(igrid,ipft,:)*dri_inf(igrid,ipft,:))/SUM(dri_inf(igrid,ipft,:))
    ELSE
       mc_i_inf_mean(igrid,ipft)  = mc_inf_temp(igrid,jst)
    ENDIF

    !! Calculation of the ratio between the mean value in the muff (previous time step) and the mean value in the soil horizon 
    !! (new time step, calculated in the subroutine hydrol_hydraulic_arch_tuzet_calc)

    IF ((mc_sup_temp(igrid,jst) .GT. min_sechiba) .AND. (mc_inf_temp(igrid,jst) .GT. min_sechiba)) THEN
       delta_mc_sup(igrid,ipft) = mc_i_sup_mean(igrid,ipft) / mc_sup_temp(igrid,jst) 
       delta_mc_inf(igrid,ipft) = mc_i_inf_mean(igrid,ipft) / mc_inf_temp(igrid,jst) 
    ELSE 
       delta_mc_sup(igrid,ipft) = 0.
       delta_mc_inf(igrid,ipft) = 0.
    ENDIF

    !! We apply the delta between the previous and the new time step to the water content at each node in the muff

    DO jrp = 1, nrp   ! Loop over the muff nodes
       IF (lr_muff_sup(igrid,ipft) .GT. min_sechiba) THEN
          IF (delta_mc_sup(igrid,ipft) .GT. min_sechiba) THEN
             mc_i_sup_temp(igrid,ipft,jrp) = MAX(mcr_sup(igrid)+0.001,MIN(mcs_sup(igrid),mc_i_sup_temp(igrid,ipft,jrp)/delta_mc_sup(igrid,ipft)))
          ENDIF
       ELSE
          mc_i_sup_temp(igrid,ipft,jrp) = mc_sup_temp(igrid,jst)
       ENDIF

       IF (lr_muff_inf(igrid,ipft) .GT. min_sechiba) THEN
          IF (delta_mc_inf(igrid,ipft) .GT. min_sechiba) THEN
             mc_i_inf_temp(igrid,ipft,jrp) = MAX(mcr_inf(igrid)+0.001,MIN(mcs_inf(igrid),mc_i_inf_temp(igrid,ipft,jrp)/delta_mc_inf(igrid,ipft)))
          ENDIF
       ELSE
          mc_i_inf_temp(igrid,ipft,jrp) = mc_inf_temp(igrid,jst)
       ENDIF
    ENDDO

    !! As the muffs are defined at the soil tile level instead of the PFT level, the total root absorption flux at the soil tile level is then
    !! the sum of the root absorption fluxes of each PFTs that composes it

    F_sup_st(igrid,ipft) = Fsup_temp(igrid,ipft)
    F_inf_st(igrid,ipft) = Finf_temp(igrid,ipft)

    !! Finally, the soil tile root absorption flux is divided by the section area of the muff and the root length in order to limit the radial 
    !! resolution of Richard's equation to a 1D resolution along the radial axis. 

    IF (lr_muff_sup(igrid,ipft) .GT. min_sechiba) THEN
       F_sup_st(igrid,ipft) = F_sup_st(igrid,ipft) * dt_sechiba * kilo_to_unit/(pi*Rad_sup(igrid,ipft)&
            *SUM(dlh(:lim_layer))*lr_muff_sup(igrid,ipft))
    ELSE
       F_sup_st(igrid,ipft) = 0.
    ENDIF

    IF (lr_muff_inf(igrid,ipft) .GT. min_sechiba) THEN
       F_inf_st(igrid,ipft) = F_inf_st(igrid,ipft) * dt_sechiba * kilo_to_unit/(pi*Rad_inf(igrid,ipft)&
            *SUM(dlh(lim_layer+1:))*lr_muff_inf(igrid,ipft))
    ELSE
       F_inf_st(igrid,ipft) = 0.
    ENDIF


    !! 2. Resolution of Richard's equation :

    !! The radial resolution is an adaptation of Patricia de Rosnay's vertical resolution of Richard's equation. The idea is exactly the same, the only
    !! modifications made on the subtoutines are to adapt the resolution to the radial axis. We call only one subroutine for each soil layer. The flag is_sup
    !! permits to make the subroutine understand that we are considering the superficial or inferior soil layer.


    is_sup= .TRUE.
    CALL hydrol_muff_radial_coef_setup(kjpindex, igrid, ipft, njsc, ks, nvan, avan, lr_muff_sup, Rad_sup, &
         dri_sup, mc_i_sup_temp, F_sup_st, is_sup, tmat_rad, rhs_rad)

    is_sup= .FALSE.
    CALL hydrol_muff_radial_coef_setup(kjpindex, igrid, ipft, njsc, ks, nvan, avan, lr_muff_inf, Rad_inf, &
         dri_inf, mc_i_inf_temp, F_inf_st, is_sup, tmat_rad, rhs_rad)

    IF ( (printlev .GE. 4) .OR. (printlev_loc .GE. 4) ) THEN
       IF (ipft.EQ.test_pft) THEN
          !WRITE(numout,*) "mc_i_sup_temp =", mc_i_sup_temp(test_grid,test_pft,:)
          !WRITE(numout,*) "mc_i_inf_temp =", mc_i_inf_temp(test_grid,test_pft,:)
          WRITE(numout,*) "lr_muff_inf =", lr_muff_inf(test_grid,test_pft)
          WRITE(numout,*) "lr_muff_sup =", lr_muff_sup(test_grid,test_pft)
          WRITE(numout,*) "F_sup_st =", F_sup_st(test_grid,test_pft)
          WRITE(numout,*) "F_inf_st =", F_inf_st(test_grid,test_pft)
       ENDIF
    ENDIF

    !! 3. Determination of the root water potential:

    !! The previous subroutines permit to calculate the soil water content in each radial layer of the muffs. At jrp=0, we assume a continuity of the water
    !! potential between the soil and the root. Consequently, the soil water content is converted to a water potential following the empiric equations of Van
    !! Genuchten (1980) or Campbell (1985) (parameterisation used in Tuzet et al. (2003-2017)).

    !! As the water content inside the muff is defined at the soil tile level, we have to determine a water content at the soil/root interface of each
    !! PFT that composes the soil tile. We decided to weight the water content at the soil tile level with the fraction of the grid occupied by the PFT.

    mc_pft_sup(igrid,ipft) = mc_i_sup_temp(igrid,ipft,1) 
    mc_pft_inf(igrid,ipft) = mc_i_inf_temp(igrid,ipft,1) 

    !! Finally, this "PFT water content" is converted into a root water potential.
    is_root = .TRUE.
    CALL calc_soil_potentials(igrid, ipft, jst, njsc, avan, nvan, mc_pft_sup, mc_pft_inf,&
                              psi_root_sup_temp, psi_root_inf_temp, is_root)

    !! To avoid too low values

    !psi_root_sup_temp(igrid,ipft) = MAX(psi_root_sup_temp(igrid,ipft),psi_root_min(ipft))
    !psi_root_inf_temp(igrid,ipft) = MAX(psi_root_inf_temp(igrid,ipft),psi_root_min(ipft))

    IF ( (printlev .GT. 4).OR.(printlev_loc>=4) ) THEN
       WRITE(numout,*) "psi_root_sup_temp calcul =", psi_root_sup_temp(igrid,ipft)
       WRITE(numout,*) "psi_root_inf_temp calcul =", psi_root_inf_temp(igrid,ipft)
       WRITE(numout,*) "mc_pft_sup =", mc_pft_sup(igrid,ipft)
       WRITE(numout,*) "mc_pft_inf =", mc_pft_inf(igrid,ipft)
       CALL flush(numout)
    ENDIF

  END SUBROUTINE hydrol_hydraulic_arch_tuzet_muff




!! ================================================================================================================================
!! SUBROUTINE 	: hydrol_muff_radial_coef_setup
!!
!>\BRIEF         This routine follows the resolution of de Rosnay's thesis. The matrix coefficients are defined before solving the 
!!               system in the next subroutine. The soil properties are also calculated
!!
!!\n DESCRIPTION : 
!!               The subroutine calculates:
!!                  1. The soil hydraulic properties according to the parameterisation wanted;
!!                  2. The matrix coefficients to solve the system.
!!
!!      
!! RECENT CHANGE(S): Added by Julien Alléon (December 2022)
!!
!! MAIN OUTPUT VARIABLE(S): tmat_rad, rhs_rad
!!
!! REFERENCE(S)	: Patricia de Rosnay's thesis
!!                hydrol.f90 technical documentation
!!
!! FLOWCHART    : 
!! 
!! ================================================================================================================================


  SUBROUTINE hydrol_muff_radial_coef_setup (kjpindex, igrid, ipft, njsc, ks, nvan, avan, lr_muff, Rad, &
       drad, mc_i, F_st, is_sup, tmat_rad, rhs_rad)

    !! Variable declaration

    !! Input Variables

    INTEGER(i_std),                               INTENT(in)          :: kjpindex             !! Domain size, terrestrial pixels only (unitless)
    INTEGER(i_std),                               INTENT(in)          :: igrid                !! Index of the grid-cell considered (unitless)
    INTEGER(i_std),                               INTENT(in)          :: ipft                 !! Index of the pft considered (unitless)
    INTEGER(i_std),     DIMENSION(:),             INTENT(in)          :: njsc                 !! Index of the dominant soil textural class (unitless)
    REAL(r_std),        DIMENSION(:),             INTENT(in)          :: ks                   !! Hydraulic conductivity at saturation (mm {-1})
    REAL(r_std),        DIMENSION(:),             INTENT(in)          :: nvan                 !! Van Genuchten coeficients n (unitless)
    REAL(r_std),        DIMENSION(:),             INTENT(in)          :: avan                 !! Van Genuchten coeficients a (mm-1})
    REAL(r_std),        DIMENSION(:,:),           INTENT(in)          :: lr_muff              !! Root length in the soil horizon (m) 
    REAL(r_std),        DIMENSION(:,:),           INTENT(in)          :: Rad                  !! Radius of the muff (m)
    REAL(r_std),        DIMENSION(:,:,:),         INTENT(in)          :: drad                 !! Radial thickness of the discretization of the muff (m)
    REAL(r_std),        DIMENSION(:,:,:),         INTENT(inout)       :: mc_i                 !! Water content at each node of the muff (m^3/m^3)
    REAL(r_std),        DIMENSION(:,:),           INTENT(inout)       :: F_st                 !! Root absorption flux at the soiltil level (sum of the flux from each PFT) (m^3/s)
    LOGICAL,                                      INTENT(in)          :: is_sup               !! Flag to determine if the resolution occurs in the superficial or inferior muff


    !! Output Variables

    REAL(r_std),        DIMENSION(:,:,:),         INTENT(out)         :: tmat_rad             !! Left hand matrix for the radial resolution
    REAL(r_std),        DIMENSION(:,:),           INTENT(out)         :: rhs_rad              !! Right hand matrix for the radial resolution

    !! Modified Variables


    !! Local Variables

    !! Soil properties
    REAL(r_std),        DIMENSION(kjpindex,nvm,nrp)                   :: d_lin_rad            !! Radial diffusivity inside the considered muff (mm^2/d)
    REAL(r_std),        DIMENSION(kjpindex,nvm,nrp)                   :: ksoil_rad            !! Radial conductivity inside the considered muff (mm/d)
    REAL(r_std)                                                       :: m                    !! VG coefficient (unitless)
    REAL(r_std)                                                       :: frac                 !! Fraction to calculate VG conductivity (unitless)

    !! Matrix coefficients
    REAL(r_std),        DIMENSION(kjpindex,nrp)                       :: er                   !! Matrix coefficient
    REAL(r_std),        DIMENSION(kjpindex,nrp)                       :: fr                   !! Matrix coefficient
    REAL(r_std),        DIMENSION(kjpindex,nrp)                       :: gr1                  !! Matrix coefficient
    REAL(r_std),        DIMENSION(kjpindex,nrp)                       :: erp                  !! Matrix coefficient
    REAL(r_std),        DIMENSION(kjpindex,nrp)                       :: frp                  !! Matrix coefficient
    REAL(r_std),        DIMENSION(kjpindex,nrp)                       :: grp                  !! Matrix coefficient
    REAL(r_std),        DIMENSION(kjpindex,nrp)                       :: dr                   !! Radial value of the diffusivity inside the muff (mm^2/d)
    REAL(r_std)                                                       :: temp3, temp4         !! Time constants used in the solving


    !! Indices
    INTEGER(i_std)                                                    :: jrp,jst,ipts,i         !! Indices (respectively: muff discretisation, 
                                                                                                !! pref_soil_veg, grid-cell, number of iterations)


!! ================================================================================================================================================

    
    !! Comments from hydrol subroutine:
    
    !-we compute tridiag matrix coefficients (LEFT and RIGHT) 
    ! of the system to solve [LEFT]*mc_{t+1}=[RIGHT]*mc{t}+[add terms]: 
    ! e(nslm),f(nslm),g1(nslm) for the [left] vector
    ! and ep(nslm),fp(nslm),gp(nslm) for the [right] vector
    
    
    !! Comments from hydrol subroutine:
    ! w_time=1 (in constantes_soil) indicates implicit computation for diffusion 
    temp3 = w_time*(dt_sechiba/(one_day))/deux
    temp4 = (un-w_time)*(dt_sechiba/(one_day))/deux


    m=(1.-1./nvan(igrid))

    DO jrp=1,nrp      ! Loop over muff discretisation

       IF (lr_muff(igrid,ipft) .GT. min_sechiba) THEN

          !! Calculation of the soil hydraulic properties according to the parameterisation wanted (Van Genuchten (1980) or Campbell (1985), used in
          !! Tuzet et al. (2003 - 2017)). We calculate the soil hydraulic conductivity and diffusivity (needed for the resolution).

          IF (is_vg) THEN

             !! Van Genuchten (1980) parameterisation

             IF (is_sup) THEN !! To distinguish the soil properties from one layer to another (to adapt)

                frac=MAX(0.001,MIN(0.999,(mc_i(igrid,ipft,jrp)-mcr_sup(igrid))/(mcs_sup(igrid)-mcr_sup(igrid))))

                ksoil_rad(igrid,ipft,jrp) =  ks(igrid)* (frac**0.5) * ( un - ( un - frac ** (un/m)) ** m )**2

                d_lin_rad(igrid,ipft,jrp) = ((ksoil_rad(igrid,ipft,jrp) / (avan(igrid)*m*nvan(igrid))) *  &
                     ( (frac**(-un/m))/(MAX(0.001,mc_i(igrid,ipft,jrp)-mcr_sup(igrid))) ) * &
                     (  frac**(-un/m) -un ) ** (-m))

             ELSE

                frac=MAX(0.001,MIN(0.999,(mc_i(igrid,ipft,jrp)-mcr_inf(igrid))/(mcs_inf(igrid)-mcr_inf(igrid))))

                ksoil_rad(igrid,ipft,jrp) =  ks(igrid)* (frac**0.5) * ( un - ( un - frac ** (un/m)) ** m )**2

                d_lin_rad(igrid,ipft,jrp) = ((ksoil_rad(igrid,ipft,jrp) / (avan(igrid)*m*nvan(igrid))) *  &
                     ( (frac**(-un/m))/(MAX(0.001,mc_i(igrid,ipft,jrp)-mcr_inf(igrid))) ) * &
                     (  frac**(-un/m) -un ) ** (-m))
             ENDIF

          ELSE

             !! Campbell (1985) parameterisation (used by Tuzet et al. (2003 - 2017))

             IF (is_sup) THEN !! To distinguish the soil properties from one layer to another (to adapt)

                d_lin_rad(igrid,ipft,jrp) = - ks(igrid) * b_muff(njsc(igrid)) * psi_air_entry(njsc(igrid)) &
                     * mega_to_unit/(cte_grav*rho_h2o*kilo_to_unit)/(mcs_sup(igrid)-mcr_sup(igrid))*&
                     ((mc_i(igrid,ipft,jrp)-mcr_sup(igrid))/(mcs_sup(igrid)-mcr_sup(igrid)))**(2.+b_muff(njsc(igrid)))

             ELSE

                d_lin_rad(igrid,ipft,jrp) = - ks(igrid) * b_muff(njsc(igrid)) * psi_air_entry(njsc(igrid)) &
                     * mega_to_unit/(cte_grav*rho_h2o*kilo_to_unit)/(mcs_inf(igrid)-mcr_inf(igrid))*&
                     ((mc_i(igrid,ipft,jrp)-mcr_inf(igrid))/(mcs_inf(igrid)-mcr_inf(igrid)))**(2.+b_muff(njsc(igrid)))

             ENDIF


          ENDIF

       ELSE

          d_lin_rad(igrid,ipft,jrp) = 0.

       ENDIF

    ENDDO

    !! As Van Genuchten (1980) parameterisation considers the conductivity to be infinitely low when approaching the residual moisture 
    !! content, if the conductivity is below a threshold, we limit it and adapt all the other values inside the muff in order to keep a 
    !! good gradient of conductivity and avoid possible counter-gradients (bugs) --> Done in hydrol.f90 too.

    IF ((lr_muff(igrid,ipft) .GT. min_sechiba) .AND. (is_vg)) THEN

       DO jrp=nrp-1,1,-1

          IF (ksoil_rad(igrid,ipft,jrp) .LT. 1e-32) THEN

             IF (is_sup) THEN

                frac=MAX(0.001,MIN(0.999,(mc_i(igrid,ipft,jrp)-mcr_sup(igrid))/(mcs_sup(igrid)-mcr_sup(igrid))))
                ksoil_rad(igrid,ipft,jrp)=ksoil_rad(igrid,ipft,jrp+1)/10.

                d_lin_rad(igrid,ipft,jrp) = ((ksoil_rad(igrid,ipft,jrp) / (avan(igrid)*m*nvan(igrid))) *  &
                     ( (frac**(-un/m))/(MAX(0.001,mc_i(igrid,ipft,jrp)-mcr_sup(igrid))) ) * &
                     (  frac**(-un/m) -un ) ** (-m))
             ELSE

                frac=MAX(0.001,MIN(0.999,(mc_i(igrid,ipft,jrp)-mcr_inf(igrid))/(mcs_inf(igrid)-mcr_inf(igrid))))
                ksoil_rad(igrid,ipft,jrp)=ksoil_rad(igrid,ipft,jrp+1)/10.

                d_lin_rad(igrid,ipft,jrp) = ((ksoil_rad(igrid,ipft,jrp) / (avan(igrid)*m*nvan(igrid))) *  &
                     ( (frac**(-un/m))/(MAX(0.001,mc_i(igrid,ipft,jrp)-mcr_inf(igrid))) ) * &
                     (  frac**(-un/m) -un ) ** (-m))

             ENDIF

          ENDIF

       ENDDO
       
    ENDIF

    !! To facilitate the reading of the calculation, we rename the diffusivity
    
    dr(igrid,:) = d_lin_rad(igrid,ipft,:)


    !! Calculation of the matrix coefficients:

    !! Coefficient for first layer

    IF (drad(igrid,ipft,2) .GT. min_sechiba) THEN

       er(igrid,1) = zero
       fr(igrid,1) = trois * drad(igrid,ipft,2)/huit  + temp3 &
            & * (dr(igrid,1)+dr(igrid,2))/drad(igrid,ipft,2)
       gr1(igrid,1) = drad(igrid,ipft,2)/(huit)       - temp3 &
            & * (dr(igrid,1)+dr(igrid,2))/drad(igrid,ipft,2)
       erp(igrid,1) = zero
       frp(igrid,1) = trois * drad(igrid,ipft,2)/huit - temp4 &
            & * (dr(igrid,1)+dr(igrid,2))/drad(igrid,ipft,2)
       grp(igrid,1) = drad(igrid,ipft,2)/(huit)       + temp4 &
            & * (dr(igrid,1)+dr(igrid,2))/drad(igrid,ipft,2)
    ELSE 

       er(igrid,1)=0
       fr(igrid,1)=0
       gr1(igrid,1)=0
       erp(igrid,1)=0
       frp(igrid,1)=0
       grp(igrid,1)=0

    ENDIF

    !! Coefficient for medium layers

    DO jrp = 2, nrp-1

       IF ((drad(igrid,ipft,jrp) .GT. min_sechiba) .AND. (drad(igrid,ipft,jrp+1) .GT. min_sechiba))  THEN

          er(igrid,jrp) = drad(igrid,ipft,jrp)/(huit)                        - temp3 &
               & * (dr(igrid,jrp)+dr(igrid,jrp-1))/drad(igrid,ipft,jrp)

          fr(igrid,jrp) = trois * (drad(igrid,ipft,jrp)+drad(igrid,ipft,jrp+1))/huit  + temp3 &
               & * ((dr(igrid,jrp)+dr(igrid,jrp-1))/(drad(igrid,ipft,jrp)) + &
               & (dr(igrid,jrp)+dr(igrid,jrp+1))/(drad(igrid,ipft,jrp+1)) )

          gr1(igrid,jrp) = drad(igrid,ipft,jrp+1)/(huit)                     - temp3 &
               & * (dr(igrid,jrp)+dr(igrid,jrp+1))/drad(igrid,ipft,jrp+1)

          erp(igrid,jrp) = drad(igrid,ipft,jrp)/(huit)                       + temp4 &
               & * (dr(igrid,jrp)+dr(igrid,jrp-1))/drad(igrid,ipft,jrp)

          frp(igrid,jrp) = trois * (drad(igrid,ipft,jrp)+drad(igrid,ipft,jrp+1))/huit - temp4 &
               & * ( (dr(igrid,jrp)+dr(igrid,jrp-1))/(drad(igrid,ipft,jrp)) + &
               & (dr(igrid,jrp)+dr(igrid,jrp+1))/(drad(igrid,ipft,jrp+1)) )

          grp(igrid,jrp) = drad(igrid,ipft,jrp+1)/(huit)                     + temp4 &
               & *(dr(igrid,jrp)+dr(igrid,jrp+1))/drad(igrid,ipft,jrp+1)             
       ELSE 

          er(igrid,jrp)=0
          fr(igrid,jrp)=0
          gr1(igrid,jrp)=0
          erp(igrid,jrp)=0
          frp(igrid,jrp)=0
          grp(igrid,jrp)=0

       ENDIF

    ENDDO

    !! Coefficient for last layer 

    IF (drad(igrid,ipft,nrp) .GT. min_sechiba) THEN

       er(igrid,nrp) = drad(igrid,ipft,nrp)/(huit)        - temp3 &
            & * (dr(igrid,nrp)+dr(igrid,nrp-1)) /drad(igrid,ipft,nrp)
       fr(igrid,nrp) = trois * drad(igrid,ipft,nrp)/huit  + temp3 &
            & * (dr(igrid,nrp)+dr(igrid,nrp-1)) / drad(igrid,ipft,nrp)
       gr1(igrid,nrp) = zero
       erp(igrid,nrp) = drad(igrid,ipft,nrp)/(huit)       + temp4 &
            & * (dr(igrid,nrp)+dr(igrid,nrp-1)) /drad(igrid,ipft,nrp)
       frp(igrid,nrp) = trois * drad(igrid,ipft,nrp)/huit - temp4 &
            & * (dr(igrid,nrp)+dr(igrid,nrp-1)) /drad(igrid,ipft,nrp)
       grp(igrid,nrp) = zero

    ELSE 

       er(igrid,nrp)=0
       fr(igrid,nrp)=0
       gr1(igrid,nrp)=0
       erp(igrid,nrp)=0
       frp(igrid,nrp)=0
       grp(igrid,nrp)=0

    ENDIF

    !! Calculation of the matrix

    !! First layer, this is where the root absorption takes place

    tmat_rad(igrid,1,1) = zero
    tmat_rad(igrid,1,2) = fr(igrid,1)
    tmat_rad(igrid,1,3) = gr1(igrid,1)
    rhs_rad(igrid,1)    = frp(igrid,1) * mc_i(igrid,ipft,1) + grp(igrid,1)*mc_i(igrid,ipft,2) &
         &  - F_st(igrid,ipft)

    !! Muff body, pure diffusion

    DO jrp=2, nrp-1
       tmat_rad(igrid,jrp,1) = er(igrid,jrp)
       tmat_rad(igrid,jrp,2) = fr(igrid,jrp)
       tmat_rad(igrid,jrp,3) = gr1(igrid,jrp)
       rhs_rad(igrid,jrp) = erp(igrid,jrp)*mc_i(igrid,ipft,jrp-1) + frp(igrid,jrp)*mc_i(igrid,ipft,jrp) &
            & +  grp(igrid,jrp) * mc_i(igrid,ipft,jrp+1)
    ENDDO

    !! Last layer, we consider no exchanges between te edge of the muff and outside it

    tmat_rad(igrid,nrp,1) = er(igrid,nrp)
    tmat_rad(igrid,nrp,2) = fr(igrid,nrp)
    tmat_rad(igrid,nrp,3) = zero
    rhs_rad(igrid,nrp) = erp(igrid,nrp)*mc_i(igrid,ipft,nrp-1) + frp(igrid,nrp)*mc_i(igrid,ipft,nrp) 


    CALL hydrol_muff_radial_resolution(kjpindex, igrid, ipft, tmat_rad, rhs_rad,is_sup, mc_i)


  END SUBROUTINE hydrol_muff_radial_coef_setup


!! ================================================================================================================================
!! SUBROUTINE 	: hydrol_muff_radial_resolution
!!
!>\BRIEF        Solves the tridiagonal system resulting from Richard's equation discretisation.   
!!
!!\n DESCRIPTION : 
!!           
!!              This subroutine solves the tridiagonal system resulting from the discretisation of the radial Richard's equation.
!!              The resolution is the same as in Patricia de Rosnay's thesis. It is only adapted to the radial resolution.
!!
!!      
!! RECENT CHANGE(S): Added by Julien Alléon (December 2022)
!!
!! MAIN OUTPUT VARIABLE(S): mc_i
!!
!! REFERENCE(S)	: Patricia de Rosnay's thesis
!!                hydrol.f90 technical note
!!
!! FLOWCHART    : 
!! 
!! ================================================================================================================================


  SUBROUTINE hydrol_muff_radial_resolution (kjpindex, igrid, ipft, tmat_rad, rhs_rad,is_sup, mc_i)

    !! Variable declaration

    !! Input Variables

    INTEGER(i_std),                               INTENT(in)          :: kjpindex             !! Domain size, terrestrial pixels only (unitless)
    INTEGER(i_std),                               INTENT(in)          :: igrid                !! Index of the grid-cell considered (unitless)
    INTEGER(i_std),                               INTENT(in)          :: ipft                 !! Index of the pft considered (unitless)
    REAL(r_std),        DIMENSION(:,:,:),         INTENT(in)          :: tmat_rad             !! Left hand matrix for the radial resolution
    REAL(r_std),        DIMENSION(:,:),           INTENT(in)          :: rhs_rad              !! Right hand matrix for the radial resolution
    LOGICAL,                                      INTENT(in)          :: is_sup               !! Flag to control if the resolution occurs in the superficial or inferior muff

    !! Output Variables



    !! Modified Variables

    REAL(r_std),        DIMENSION(:,:,:),         INTENT(inout)       :: mc_i                 !! Water content at each node of the muff (m^3/m^3)

    !! Local Variables

    !! Matrix
    REAL(r_std),        DIMENSION(kjpindex)                           :: bet_rad              !! Resolution term
    REAL(r_std),        DIMENSION(kjpindex,nrp)                       :: gam_rad              !! Resolution term

    !! Indices
    INTEGER(i_std)                                                    :: ipts,jrp             !! Indices (respectively: grid-cells, muff discretisation)


!! ================================================================================================================================================


    bet_rad(igrid) = tmat_rad(igrid,1,2)

    IF (bet_rad(igrid) .GT. min_sechiba) THEN

       mc_i(igrid,ipft,1) = rhs_rad(igrid,1)/bet_rad(igrid)

    ENDIF

    gam_rad(igrid,1) = 0.

    DO jrp = 2,nrp

       IF (bet_rad(igrid) .GT. min_sechiba) THEN

          gam_rad(igrid,jrp) = tmat_rad(igrid,jrp-1,3)/bet_rad(igrid)

       ELSE

          gam_rad(igrid,jrp) = 0.

       ENDIF

       bet_rad(igrid) = tmat_rad(igrid,jrp,2) - tmat_rad(igrid,jrp,1)*gam_rad(igrid,jrp)

       IF (bet_rad(igrid) .GT. min_sechiba) THEN

          mc_i(igrid,ipft,jrp) = (rhs_rad(igrid,jrp)-tmat_rad(igrid,jrp,1)*mc_i(igrid,ipft,jrp-1))/bet_rad(igrid)

       ENDIF

    ENDDO

    DO jrp = nrp-1,1,-1
       IF (is_sup) THEN
          mc_i(igrid,ipft,jrp) = MAX(mcr_sup(igrid)+0.0001,mc_i(igrid,ipft,jrp) - gam_rad(igrid,jrp+1)*mc_i(igrid,ipft,jrp+1))
       ELSE
          mc_i(igrid,ipft,jrp) = MAX(mcr_inf(igrid)+0.0001,mc_i(igrid,ipft,jrp) - gam_rad(igrid,jrp+1)*mc_i(igrid,ipft,jrp+1))
       ENDIF
    ENDDO

  END SUBROUTINE hydrol_muff_radial_resolution



!! ================================================================================================================================
  !! SUBROUTINE   : hydrol_root_profile
  !!
  !>\BRIEF         Calculates the share of the root biomass in each soil layer based on
  !!               structural and functional approach.Calculate the root profile
  !!
  !! DESCRIPTION  : Root structure is probably how most 
  !! of us think about roots (i.e. digging a whole and observing where the roots 
  !! are). When thinking at root structure, the profile should be relatively 
  !! constant over time. A logic time integrator to determine this constancy is 
  !! longevity_root as the profile cannot grow faster then the roots grow and die.
  !! Currently the structural root profile is simply fixed it over time. This could
  !! be changed by, for example, making humcste a function of the tree diameter.
  
  !! In ORCHIDEE root structure is used in the calculation of kfact_root which is water 
  !! infiltration along roots (accounted for in hydrol.f90) and the input of soil 
  !! carbon and nitrogen at depth due to the turnover of roots which is accounted 
  !! for stomate_soil_carbon_discretization.f90. Furthermore, it is used to calculate 
  !! the root temperature in stomate_resp.f90. When thinking about root function it 
  !! is not so important where the roots are located but it is more important at 
  !! which depth the roots will be active. The function approach could be used to
  !! calculate from which soil layers the plants take most the soil water for their
  !! transpiration. This way of looking at the roots is similar to how we look at 
  !! the canopy where we have a lot of leaves at places in the canopy where little
  !! light can penetrate and where a large part of the photosynthesis is taken care 
  !! of by the leaves in the top layers of the canopy.
  !!
  !! NOTE: for the moment root structure and root function are only coupled through
  !! the depth of the soil. In the absence of roots, the root functions cannot be 
  !! fullfilled. This is the most minimalistic coupling. It basically implies that
  !! a very small fraction of the roots, e.g., < 1% could take up all the water 
  !! required for transpiration. A thighter coupling between structure and function
  !! is to be expected but this needs to be checked in the literature.
  !!
  !! RECENT CHANGE(S) : None
  !!
  !! MAIN OUTPUT VARIABLE(S) : root_profile
  !!
  !! REFERENCE(S) : 
  !!
  !! \n
  !_ ================================================================================================================================

  SUBROUTINE hydrol_root_profile(kjpindex, altmax, sm, smw, root_profile, root_depth)

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                         :: kjpindex           !! Domain size
    REAL(r_std), DIMENSION(:,:), INTENT(in)            :: altmax             !! Maximul active layer thickness (m). Be careful, here active means not frozen.
                                                                             !! Not related with the active soil carbon pool.
    REAL(r_std), DIMENSION (:,:,:), INTENT(in)         :: sm                 !! Soil moisture of each layer for every soil tile (liquid phase)
                                                                             !!  @tex $(kg m^{-2})$ @endtex
    REAL(r_std), DIMENSION (:,:), INTENT(in)           :: smw                !! Soil moisture of each layer at wilting point
                                                                             !!  @tex $(kg m^{-2})$ @endtex 

    !! 0.2 Output variables
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(out)       :: root_profile       !! Normalized root mass/length fraction in each soil layer 
                                                                             !! (0-1, unitless)
                                                                             !! DIM = kjpindex * nvm * nslm
    REAL(r_std), DIMENSION (:,:,:), INTENT(out)        :: root_depth         !! Node and interface numbers at which the deepest roots
                                                                             !! occur (1 to nslm, unitless)

    
    !! 0.3 Modified variables
    
    !! 0.4 Local variables
    INTEGER(i_std)                                     :: ji,jv,jsl,jst      !! Indices  
    REAL(r_std)                                        :: rpc                !! Integration constant for vertical decomposer 
    REAL(r_std)                                        :: z_top, z_bottom    !! top and bottom node in between which to integrate the root profile
    REAL(r_std), DIMENSION(kjpindex)                   :: count              !! Count the number of errors
    REAL(r_std), DIMENSION(nslm)                       :: root_profile_tmp   !! Temporary variable
    REAL(r_std)                                        :: root_depth_tmp     !! Temporary variable
    
!_ ================================================================================================================================

  !! 1.Rooting depth

    ! Calculate the maximum depth roots occur at. Two constraints are already accouned
    ! for: (1) crop can root to 0.8 m, grasslands are assumed to root no deeper than 1 m. 
    ! Trees can root down to 2 m. (2) Roots do not extend into frozen soils.
    ! ORCHIDEE assumes that plant roots go down to the depth of the soil water profile
    ! right from day 1. In other words, even a very young tree sapling, crop or grasslands
    ! has roots that extend down to their :: max_root_depth.
    ! Hence, plant age/height does NOT constraint rooting depth for now. A future
    ! development could make rooting depth a function of plant height. This approach
    ! would correctly increase the water stress of young vegetation but may wrongly
    ! increase the water stress of young plants growing in semi-arid regions. In such
    ! regions vegetative reproduction is likley to be common as it allows the offspring
    ! to use water from the parent plant as long as the offsprong is too small to reach
    ! the deeper water layers.
    root_depth(:,:,:) = zero
    DO ji = 1,kjpindex
       DO jv = 1,nvm
           
          ! Plants can root up to 2 m (zdr(jsl), the prescribed max_root_depth
          ! or the first permafrost layer
          root_depth_tmp = MIN(MIN(altmax(ji,jv),maxaltmax), &
               MIN(zdr(nslm),max_root_depth(jv)))
             
          ! Find the index of the node which is the closest to the actual 
          ! root_depth. This layer is used to truncate the root profile.
          root_depth(ji,jv,inode) = MINLOC(ABS(root_depth_tmp-znh(:)),DIM=1)

          ! zdr is defined as a 0:nslm matrix. MINLOC does not know that 
          ! and gives the indices ad 1:nslm+1. Convert the indices back to 
          ! 0:nslm by subtracting 1.
          root_depth(ji,jv,iinterface) = MINLOC(ABS(root_depth_tmp-zdr(:)),DIM=1)-1

       END DO ! jv
    END DO ! ji
   
    ! Prescribe a solution for very shallow root systems
    WHERE (root_depth(:,:,inode).LE.2)
       root_depth(:,:,inode) = 2
    ENDWHERE

    WHERE (root_depth(:,:,iinterface).LE.2)
       root_depth(:,:,iinterface) = 2
    ENDWHERE
    

  !! 2.Structural root profile
 
    ! The structural root profile is calculated as an exponentially decreasing 
    ! root mass with depth which. ORCHIDEE uses the structural root profile to
    ! calculate rain infiltration along the roots, som inputs at depth and the
    ! root temperature for autotrophic respiration. 
    
    ! NOTE: the shape of the profile is determined by its depth and the parameter 
    ! humcste. For the moment humcste is PFT-dependent but constant over time. The 
    ! depth is PFT-dependent (see above) and is a function of the active layer
    ! thickness when the soil discretization is used. The depth of the root profile
    ! could/should be made a function of tree diameter or plant biomass. 
    root_profile(:,:,:,istruc) = zero
    DO jv = 1,nvm
       DO ji = 1,kjpindex

          ! Note that the integration will start at the top of the second layer (the 
          ! top of the first layer is zdr(0) hence the zdr(1) is the top of the 
          ! second layer) and will continue until the bottom of the profile. The
          ! bottom is the profile is calculated above but never extends deeper than 
          ! the depth used in hydrol.f90 (in thermosoil.f90 the profile extends much
          ! deeper. The first layer is excluded because the profile will be used to
          ! extract water from the soil. The first layer is very thin and by extracting
          ! water it could dry to quickly. Also in reality not too many roots are found
          ! in the top mm of the soil. zdr describes the nodes and interfaces of the soil
          ! layers as proposed by de Rosnay 1999 (PhD thesis).
          z_top = zdr(1)
          z_bottom = zdr(root_depth(ji,jv,iinterface))

          ! Calculate the total surface area under an exponential curve between
          ! zdr(1) and the zdr(nslm)
          rpc = un / ( EXP(-z_top * humcste(jv)) - EXP(-z_bottom*humcste(jv)) )

          DO jsl = 2, root_depth(ji,jv,iinterface)

             ! Calculate the share of the total roots for layers which are "centered"
             ! at the nodes of the hydrology scheme. Centered was written in quotes
             ! because the layer is not symmetric. Using the nodes as the center of the
             ! layers follows De Rosnay (PhD, figure C.2 page 156). The root profile
             ! starts at the top of the second layer, ends at the bottom of the 11th 
             ! layer and is calculated for the nodes in between.
             ! The following equation are derived (but rewritten) from the integrals 
             ! of equations C9 to C11 of De Rosnay's (1999) PhD thesis (page 158).
             root_profile(ji,jv,jsl,istruc) = rpc * &
                  ( EXP(-zdr(jsl-1)*humcste(jv)) - EXP(-zdr(jsl)*humcste(jv)) )

          END DO ! root_depth

          ! Top layer does not contain structural roots (see z_top)
          ! This line is not needed but was added as a reminder.
          root_profile(ji,jv,1,istruc) = zero

          ! Error checking. Each root profile should add up to 1.
          IF (err_act.GT.1) THEN
             IF (ABS(SUM(root_profile(ji,jv,:,istruc))-un).GT.100*EPSILON(un)) THEN
                WRITE(numout,*) 'pixel, PFT, sum of structural root_profile, ', ji, jv, &
                     SUM(root_profile(ji,jv,:,istruc))
                WRITE(numout,*) 'hydrol_root_profile, structural root_profile, ', &
                     root_profile(ji,jv,:,istruc)
                CALL ipslerr_p(plev,'hydrol.f90', &
                     'structural root profile does not add up to 1', &
                     'Check its calculation', '')
             ENDIF
          END IF

       END DO ! kjpindex
    ENDDO ! nvm


  !! 3. Functional root profile

    ! Calculates the share of the root biomass in each soil layer based on
    ! the soil water content in each layer. The roots now follow the water.
    ! this results in a very dynamic root profile that may change every half
    ! hour (i.e. the time step for hydrol). ORCHIDEE uses the root function
    ! to calculate plant water stress (hydrol_hydraulic_arch_tuzet_calc)
    ! and the soil layers from which the transpiration is taken (hydrol.f90).
    
    ! NOTE: yet another root function is nutrient uptake. A separate root 
    ! profile could be calculated to be used in the calculation of N uptake
    ! is stomate_soilcarbon.f90 (nitrogen_dynamics).

    
    root_profile(:,:,:,ifunc) = zero
    DO jv = 2, nvm
       jst = pref_soil_veg(jv) ! soil tile index
       DO ji = 1, kjpindex

          ! Plant available soil water per layer. 
          ! The calculations could make use of sm and smw. These variables has 
          ! been labelled soil moisture in kg/m2 and soil moisture at each layer 
          ! at wilting point also in kg/m2. As an alternative swc and mcr could
          ! be used. swc is calculated in hydrol.f90 based on mc (m3/m3). mc is 
          ! used to calculate smt (total soil water thus liquid + ice) and mcl is 
          ! used to calculate sm (liquid only). sm denotes only liquid water, swc
          ! denotes liquid and frozen water. Given the focus on root function, a
          ! variable describing the liquid water seems the better choice.
          
          root_profile_tmp(:) = zero
          DO jsl = 2, root_depth(ji,jv,iinterface)
             root_profile_tmp(jsl) = MAX(zero,sm(ji,jsl,jst)-smw(ji,jsl) )
          ENDDO
          
          ! Normalize to obtain a root profile (fraction between 0 and 1)
          IF (SUM(root_profile_tmp(:)) .GT. min_sechiba ) THEN
             root_profile(ji,jv,:,ifunc) = root_profile_tmp(:) / &
                  SUM(root_profile_tmp(:))
          ELSE
             root_profile(ji,jv,:,ifunc) = zero
          END IF

          ! Top layer is not used for water uptake
          ! This line is not needed but was added as a reminder.
          root_profile(ji,jv,1,ifunc) = zero

          ! The functional profile should also equal to one if there is soil moisture.
          IF (err_act.GT.1) THEN
             IF (ABS(SUM(root_profile(ji,jv,:,ifunc))-un).GT.100*EPSILON(un) .AND. &
                  SUM(root_profile_tmp(:)) .GT. min_sechiba) THEN
                WRITE(numout,*) 'pixel, PFT, sum of functional root_profile, ', ji, jv, &
                     SUM(root_profile(ji,jv,:,ifunc))
                WRITE(numout,*) 'hydrol_root_profile, functional root_profile, ', &
                     root_profile(ji,jv,:,ifunc)
                CALL ipslerr_p(plev,'hydrol.f90', &
                     'functional root profile does not add up to 1', &
                     'Check its calculation', '')
             ENDIF
          ENDIF

       END DO ! kjpindex
    ENDDO ! ivm 

  END SUBROUTINE hydrol_root_profile


!! ================================================================================================================================
!! SUBROUTINE 	 : calc_soil_potentials
!!
!>\BRIEF        Uses Campbell's or Van Genuchten to calculate soil water potentials 
!!
!!\n DESCRIPTION : 
!!           
!!      
!! RECENT CHANGE :
!! MAIN OUTPUT VARIABLE(S): psi_soil sup and psi_soil_inf
!!
!! REFERENCE(S)	:
!!              
!!
!! FLOWCHART    : 
!! 
!! ================================================================================================================================


  SUBROUTINE calc_soil_potentials (igrid, ipft, jst, njsc, avan, nvan, mc_sup, mc_inf, psi_sup, psi_inf, is_root)

    !! Variable declaration

    !! Input Variables

    LOGICAL,                                      INTENT(in)          :: is_root               !! Flag to check if we're calculating root or soil potential
    INTEGER(i_std),                               INTENT(in)          :: igrid                 !! Index of the grid-cell considered (unitless)
    INTEGER(i_std),                               INTENT(in)          :: ipft                  !! Index of the pft considered (unitless)
    INTEGER(i_std),                               INTENT(in)          :: jst                   !! Index of soil tiles (unitless, 1-3)
    INTEGER(i_std),     DIMENSION(:),             INTENT(in)          :: njsc                  !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless)
    REAL(r_std),        DIMENSION(:),             INTENT(in)          :: avan                  !! Van Genuchten coeficients a (mm-1})
    REAL(r_std),        DIMENSION(:),             INTENT(in)          :: nvan             !! Van Genuchten coeficients n (unitless)
    REAL(r_std),        DIMENSION(:,:),           INTENT(in)          :: mc_sup                !! Temporary value of the water content in the superficial soil level (m^3/m^3)
    REAL(r_std),        DIMENSION(:,:),           INTENT(in)          :: mc_inf                !! Temporary value of the water content in the inferior soil level (m^3/m^3)

    !! Output Variables

    !! Modified Variables

    REAL(r_std),       DIMENSION(:,:),            INTENT(inout)       :: psi_sup               !! Superficial soil Water Potential (MPa)
    REAL(r_std),       DIMENSION(:,:),            INTENT(inout)       :: psi_inf               !! Inferior soil Water Potential (MPa)

    !! Local Variables

    !! Matrix

    !! Indices
    INTEGER(i_std)                                                    :: ipsi                 !! Used to distinguish the root case
    

!! ================================================================================================================================================

    IF (is_root) THEN
       ipsi = ipft
    ELSE
       ipsi = jst
    ENDIF

    IF (mc_sup(igrid,ipsi) .GT. min_stomate .AND. &
         (mc_sup(igrid,ipsi) - mcr_sup(igrid)) .GT. min_stomate) THEN

       IF (is_vg) THEN

          !! Parameterisation of Van Genuchten (1980)

          psi_sup(igrid,ipsi) =  &
               - ( cte_grav * rho_h2o * kilo_to_unit * & 
               ( 1./(avan(igrid)*kilo_to_unit) ) * &
               ( ( (mc_sup(igrid,ipsi) - mcr_sup(igrid) ) / &    
               ( mcs_sup(igrid) - mcr_sup(igrid) ) ) ** &
               ( -1. / (1.-1./nvan(igrid)) ) - 1.) ** &
               (1/nvan(igrid) ) ) / mega_to_unit

       ELSE

          !! Parameterisation of Campbell (1985)

          psi_sup(igrid,ipsi)=psi_air_entry(njsc(igrid))*&
               & ((mcs_sup(igrid)-mcr_sup(igrid))/(mc_sup(igrid,ipsi)-mcr_sup(igrid)))**b_muff(njsc(igrid))

       ENDIF

    ELSE !! To avoid bugs when mc_sup appears to be lesser than mcr

       psi_sup(igrid,ipsi) = psi_root_min(ipft)
       WRITE(numout,*) 'Superior layer below mcr'
 
    ENDIF

    IF (mc_inf(igrid,ipsi) .GT. min_stomate .AND. &
         (mc_inf(igrid,ipsi) - mcr_inf(igrid)) .GT. min_stomate) THEN

       IF (is_vg) THEN


          !! Parameterisation of Van Genuchten (1980)

          psi_inf(igrid,ipsi) =  &
               - ( cte_grav * rho_h2o * kilo_to_unit * & 
               ( 1./(avan(igrid)*kilo_to_unit) ) * &
               ( ( (mc_inf(igrid,ipsi) - mcr_inf(igrid) ) / &    
               ( mcs_inf(igrid) - mcr_inf(igrid) ) ) ** &
               ( -1. / (1.-1./nvan(igrid)) ) - 1.) ** &
               (1/nvan(igrid) ) ) / mega_to_unit

       ELSE


          !! Parameterisation of Campbell (1985)

          psi_inf(igrid,ipsi)=psi_air_entry(njsc(igrid))*&
               & ((mcs_inf(igrid)-mcr_inf(igrid))/(mc_inf(igrid,ipsi)-mcr_inf(igrid)))**b_muff(njsc(igrid))

       ENDIF

    ELSE !! To avoid bugs when mc_inf appears to be lesser than mcr

       psi_inf(igrid,ipsi) = psi_root_min(ipft)
       WRITE(numout,*) 'Inferior layer below mcr'
    
    ENDIF

    psi_inf(igrid,ipsi) = MAX(psi_root_min(ipft),psi_inf(igrid,ipsi))
    psi_sup(igrid,ipsi) = MAX(psi_root_min(ipft),psi_sup(igrid,ipsi))
 
  END SUBROUTINE calc_soil_potentials

!! ================================================================================================================================
!! SUBROUTINE 	 : resistances_curtail
!!
!>\BRIEF        Curtailement of the hydraulic resistances. It has been checked that it happens sufficiently rarely to be a pure 
!!              numerical fix and does not often interfere with the scheme  
!!
!!\n DESCRIPTION : 
!!           
!!      
!! RECENT CHANGE :
!! MAIN OUTPUT VARIABLE(S): Res
!!
!! REFERENCE(S)	:
!!              
!!
!! FLOWCHART    : 
!! 
!! ================================================================================================================================


  SUBROUTINE resistances_curtail (ipts, ivm, icir, Flux, Res, psi_temp, grav_grad)

    !! Variable declaration

    !! Input Variables

    INTEGER(i_std),                               INTENT(in)          :: ipts                  !! Index of the grid-cell considered (unitless)
    INTEGER(i_std),                               INTENT(in)          :: ivm                   !! Index of the pft considered (unitless)
    INTEGER(i_std),                               INTENT(in)          :: icir                  !! Index of the circ class considered (unitless)
    REAL(r_std),                                  INTENT(in)          :: grav_grad             !! gravity gradient difference in hydraulic potential
    REAL(r_std),        DIMENSION(:,:),           INTENT(in)          :: Flux                  !! water flux through the organ considered
    REAL(r_std),        DIMENSION(:,:),           INTENT(in)          :: psi_temp              !! water potential at the organ considered 

    !! Output Variables

    !! Modified Variables

    REAL(r_std),       DIMENSION(:,:,:),          INTENT(inout)       :: Res                   !! Hydraulic resistance of the organ considered

    !! Local Variables
    !! Matrix
    !! Indices
    
    IF ( (-Flux(ipts,ivm) * Res(ipts,ivm,icir) + psi_temp(ipts,ivm) - grav_grad < psi_min).AND.(Flux(ipts,ivm).GT.zero) &
        .OR. ((Res(ipts,ivm,icir) .EQ. undef_sechiba).AND.(Flux(ipts,ivm) .GT. min_sechiba)) ) THEN 
        Res(ipts,ivm,icir) = (psi_temp(ipts,ivm) - grav_grad - psi_min) / Flux(ipts,ivm)
      IF (printlev_loc>=4) WRITE(numout,*) "resistance curtailed"
    ELSE
      IF (printlev_loc>=4) WRITE(numout,*) "resistance not curtailed"
    ENDIF

  END SUBROUTINE resistances_curtail

!! ================================================================================================================================================



END MODULE hydrol
