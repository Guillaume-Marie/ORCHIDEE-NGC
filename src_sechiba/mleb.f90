!  ==============================================================================================================================
!  MODULE		 			: mleb
!
!  CONTACT		 			: orchidee-help _at_ listes.ipsl.fr
!
!  LICENCE	 			        : IPSL (2006)
!  This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF					This module computes the energy balance on 
!! continental surfaces either with the multi-layer energy scheme, which collapses to the
!! original enerbil scheme, when jnlvls is equal to 1.
!! The module contains the following subroutines: mleb_initialize, mleb_main, mleb_begin,
!! mleb_profile, mleb_lwdrad, mleb_boundarylayer_resistance,
!! mleb_stomatal_resistance, mleb_alpha_beta_coeff, mleb_calc_column,
!! mleb_pottemp, mleb_flux, mleb_evapveg
!!
!!\n DESCRIPTION                                : 
!! \n
!! \latexonly 
!!     \input{enerbil_intro2.tex}
!! \endlatexonly
!!
!! IMPORTANT NOTE: The coefficients A and B are defined differently than those in the referenced 
!! literature and from those in the code and documentation of the atmospheric model LMDZ. For the
!! avoidance of doubt, the coefficients as described here always refer to the ORCHIDEE coefficients, 
!! and are denoted as such (with the marker: ORC). The re-definition of the coefficients takes place 
!! within LMDZ before they are passed to ORCHIDEE. The following sequence of expressions is to be 
!! found within the LMDZ module 'surf_land_orchidee':\n
!!
!! \latexonly 
!!     \input{surflandLMDZ1.tex}
!!     \input{surflandLMDZ2.tex}
!!     \input{surflandLMDZ3.tex}
!!     \input{surflandLMDZ4.tex}
!! \endlatexonly
!! \n
!!
!! \latexonly 
!!     \input{enerbil_symbols.tex}
!! \endlatexonly
!!
!! RECENT CHANGE(S)                             : The multi-layer energy budget
!code was moved to its own module in r6346. 
!!
!! REFERENCE(S)	                                : Ryder et al., 2016, GMD, doi:10.5194/gmd-9-223-2016 
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_sechiba/mleb.f90 $
!! $Date: 2025-11-06 17:58:38 +0100 (jeu. 06 nov. 2025) $
!! $Revision: 9155 $
!! \n
!_ ================================================================================================================================

MODULE mleb

  ! routines called : restput, restget
  !
  ! modules used
  USE ioipsl
  USE xios_orchidee
  USE ioipsl_para 
  USE constantes
  USE time, ONLY : one_day, dt_sechiba
  USE pft_parameters
  USE qsat_moisture
  USE sechiba_io_p
  USE constantes_soil
  USE explicitsnow

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: mleb_write, mleb_initialize, mleb_finalize, mleb_main, mleb_clear

  ! variables used inside enerbil module : declaration and initialisation

  INTEGER(i_std), SAVE                          :: printlev_loc  !! Local printlev for this module
!$OMP THREADPRIVATE(printlev_loc)

  ! one dimension array allocated, computed and used in enerbil module exclusively

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:) :: psold         !! Old surface dry static energy (J kg^{-1})
!$OMP THREADPRIVATE(psold)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:) :: qsol_sat      !! Saturated specific humudity for old temperature (kg kg^{-1})
!$OMP THREADPRIVATE(qsol_sat)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:) :: pdqsold       !! Deriv. of saturated specific humidity at old temp
  !! (kg (kg s)^{-1})
!$OMP THREADPRIVATE(pdqsold)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:) :: psnew         !! New surface static energy (J kg^{-1})
!$OMP THREADPRIVATE(psnew)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:) :: qsol_sat_new  !! New saturated surface air moisture (kg kg^{-1})
!$OMP THREADPRIVATE(qsol_sat_new)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:) :: lwabs         !! LW radiation absorbed by the surface (W m^{-2})
!$OMP THREADPRIVATE(lwabs)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:) :: lwup          !! Long-wave up radiation (W m^{-2})
!$OMP THREADPRIVATE(lwup)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:) :: lwnet         !! Net Long-wave radiation (W m^{-2})
!$OMP THREADPRIVATE(lwnet)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:) :: fluxsubli     !! Energy of sublimation (mm day^{-1})
!$OMP THREADPRIVATE(fluxsubli)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:) :: qsat_air      !! Air saturated specific humidity (kg kg^{-1})
!$OMP THREADPRIVATE(qsat_air)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:) :: tair          !! Air temperature (K)
!$OMP THREADPRIVATE(tair)

  !
  ! Multi-layer energy budget variables that are used throughout the program / module
  !

  REAL(r_std),SAVE                                :: canopy_height = 35.0d0       !! canopy height (m)
!$OMP THREADPRIVATE(canopy_height)
  REAL(r_std),SAVE                                :: rho_veg = 1000.0d0           !! density of vegetation
!$OMP THREADPRIVATE(rho_veg)
  REAL(r_std),SAVE                                :: leaf_tks = 0.001d0           !! leaf thickness
!$OMP THREADPRIVATE(leaf_tks)

  !
  ! Various parameterisation  factors
  !
  REAL(r_std), SAVE                               :: r_sto_v_fac
!$OMP THREADPRIVATE(r_sto_v_fac)
  REAL(r_std), SAVE                               :: big_r_h2o_fac
!$OMP THREADPRIVATE(big_r_h2o_fac)

  REAL(r_std), SAVE                               :: big_k_lw       !! LW extinction coefficient
!$OMP THREADPRIVATE(big_k_lw)
  REAL(r_std), SAVE                               :: big_k_sw       !! SW extinction coefficient 
!$OMP THREADPRIVATE(big_k_sw)

  REAL(r_std), SAVE                               :: k_ext_fac
!$OMP THREADPRIVATE(k_ext_fac)

  REAL(r_std), SAVE                               :: k_eddy_slope   !! coefficients for S-shape fuction for K_eddy
!$OMP THREADPRIVATE(k_eddy_slope)
  REAL(r_std), SAVE                               :: k_eddy_ustar
!$OMP THREADPRIVATE(k_eddy_ustar)

  REAL(r_std), SAVE                               :: ks_tune        !! coefficients for S-shape fuction for KS (K_eddy only for surface)
!$OMP THREADPRIVATE(ks_tune)
  REAL(r_std), SAVE                               :: ks_slope
!$OMP THREADPRIVATE(ks_slope)
  REAL(r_std), SAVE                               :: ks_veget
!$OMP THREADPRIVATE(ks_veget)

  REAL(r_std), SAVE                               :: sr_fac         !! stomatal resistance factor 
!$OMP THREADPRIVATE(sr_fac)
  REAL(r_std), SAVE                               :: br_fac         !! boundary resistance factor
!$OMP THREADPRIVATE(br_fac)

  REAL(r_std), SAVE                               :: surf_lev
!$OMP THREADPRIVATE(surf_lev)

  LOGICAL, SAVE                                   :: WC02    !! logical paramter for applying (Wohlfahrt & Cernusca 2002) Dynamic drag 
!$OMP THREADPRIVATE(WC02)
  REAL(r_std), SAVE                               :: a_1     !! a1 to a5: parameters for the sub-canopy wind profile  
  REAL(r_std), SAVE                               :: a_2     !! the default values were set as the 
  REAL(r_std), SAVE                               :: a_3     !! a1= 0.065; a2=0.001; a3=0.434; a4= -0.751; a5=0.071 
  REAL(r_std), SAVE                               :: a_4
  REAL(r_std), SAVE                               :: a_5
!$OMP THREADPRIVATE(a_1)
!$OMP THREADPRIVATE(a_2)
!$OMP THREADPRIVATE(a_3)
!$OMP THREADPRIVATE(a_4)
!$OMP THREADPRIVATE(a_5)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: jtheta                   !! specific heat capacity of air (J / (kg K) )
!$OMP THREADPRIVATE(jtheta)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: pad_deltah               !! Plant Area Density as a function of height above the ground (m^2/m^3)
!$OMP THREADPRIVATE(pad_deltah)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: box_height               !! height of each box edge, starting at the surface layer (m)
!$OMP THREADPRIVATE(box_height)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: box_height_old           !! divisionof canopy levels from the previous time step
!$OMP THREADPRIVATE(box_height_old)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: temp_atmos_pres_grid     !! Atmospheric temperature down the column (K) PRESENT STEP
!$OMP THREADPRIVATE(temp_atmos_pres_grid)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: q_atmos_pres_grid        !! Atmospheric specific humididy down the column (kg/kg) PRESENT STEP
!$OMP THREADPRIVATE(q_atmos_pres_grid)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: temp_leaf_pres_grid      !! Leaf temperature down the column (K) PRESENT STEP
!$OMP THREADPRIVATE(temp_leaf_pres_grid)
  REAL(r_std),ALLOCATABLE, SAVE, DIMENSION (:)     :: temp_surf_pres           !! temperature at surface (K) PRESENT STEP
!$OMP THREADPRIVATE(temp_surf_pres)

  !
  ! Settings for the speciffic multi-layer energy budget uotput file  
  !
  CHARACTER(LEN=256), SAVE                        :: mleb_cdfname             !!To store the name of the multilayer netCDF output file 
  !! (seperate to the rest)
!$OMP THREADPRIVATE(mleb_cdfname)
  LOGICAL,SAVE                                    :: first_call_mleb_histwrite!! detemine when to open and write headers
  !! in the mleb_netcdf file
!$OMP THREADPRIVATE(first_call_mleb_histwrite)



CONTAINS

  ! --------------------------------------------------------------------------------------

  ! insert function for the multilayer LW radiation scheme at the start of module

  REAL(r_std) FUNCTION gu (x_in)

    REAL(r_std) :: k_ext = 0.9 
    REAL(r_std) :: x_in, k_ext_fac

    !+++ TEMP +++
    k_ext_fac = 1.0
    k_ext=k_ext*k_ext_fac 
    !++++++++++++
    gu = exp(-k_ext*x_in)    
  END FUNCTION gu ! gu(x_in)


  ! --------------------------------------------------------------------------------------




  !!  =============================================================================================================================
  !! SUBROUTINE		 			: mleb_initialize
  !!
  !> BRIEF					    : Initisalisation of the multi-layer energy budget model
  !!
  !! DESCRIPTION				: Initisalisation of the multi-layer energy budget model, allocate module variables,
  !!                              read from the run.def, set some of the parameterisation variables and read from the restart files.
  !!
  !! MAIN OUTPUT VARIABLE(S)	:
  !!
  !! REFERENCE(S)				:
  !!
  !! FLOWCHART  : None
  !! \n
  !_ ==============================================================================================================================

  SUBROUTINE mleb_initialize(kjit, kjpindex, rest_id,                    &
       temp_air, qair,                             &
       temp_sol, temp_sol_new, tsol_rad,           &
       evapot,   evapot_corr,  qsurf,    fluxsens, &
       fluxlat,  vevapp,                           & 
       u_speed, z_array_out,  &
       max_height_store )

    IMPLICIT NONE

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                            :: kjit             !! Time step number (unitless) 
    INTEGER(i_std), INTENT(in)                            :: kjpindex         !! Domain size (unitless)
    INTEGER(i_std),INTENT (in)                            :: rest_id          !! Restart_ file and history file identifier (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)         :: temp_air         !! Air temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)         :: qair             !! Lowest level specific humidity
    !! @tex $(kg kg^{-1})$ @endtex
    !! 0.2 Output variables

    !! 0.2 Output variables
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)        :: temp_sol         !! Surface temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)        :: temp_sol_new     !! Surface temperature at the new time step(K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)        :: tsol_rad         !! Radiative surface temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)        :: evapot           !! Surface Potential Evaporation (mm/tstep)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)        :: evapot_corr      !! Surface Potential Evaporation corrected with Milly's Correction (mm/tstep)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)        :: qsurf            !! Surface specific humidity (kg kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)        :: fluxsens         !! Sensible heat flux (W m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)        :: fluxlat          !! Latent heat flux (W m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)        :: vevapp           !! Total evaporation flux (mm day^{-1})
    REAL(r_std), DIMENSION(kjpindex,jnlvls),INTENT(out)   :: u_speed          !! Canopy wind speed profile - see comment in mleb_finalize

    REAL(r_std), DIMENSION (kjpindex,nvm,ncirc,nlevels_tot), INTENT(inout)  :: z_array_out         !! Heights of levels inside the canopy (from stomate) (m)
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(out)                      :: max_height_store    !! Maximum height of the canopy (m)


    !! 0.4 Local variables
    INTEGER(i_std)                                     :: ier                 !! Allocation variable (unitless)
    INTEGER(i_std)                                     :: ji, ilev            !! Indice (respectively: grid-cells, canopy levels) (unitless)
    CHARACTER(LEN=80) :: var_name


    !_ ================================================================================================================================

    IF (printlev>=3) WRITE (numout,*) 'mleb_initialize : Start initillization'

    ! Initialize local printlev
!!$    printlev_loc=get_printlev('mleb')
    !printlev_loc=4

    !! 1. Allocate module variables
    ALLOCATE (psold(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'mleb_initialize','Problem in allocation of variable ','','')

    ALLOCATE (qsol_sat(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'mleb_initialize','Problem in allocation of variable ','','')

    ALLOCATE (pdqsold(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'mleb_initialize','Problem in allocation of variable ','','')

    ALLOCATE (psnew(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'mleb_initialize','Problem in allocation of variable ','','')

    ALLOCATE (qsol_sat_new(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'mleb_initialize','Problem in allocation of variable qsol_sat_new','','')

    ALLOCATE (lwabs(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'mleb_initialize','Problem in allocation of variable lwabs','','')

    ALLOCATE (lwup(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'mleb_initialize','Problem in allocation of variable lwup','','')

    ALLOCATE (lwnet(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'mleb_initialize','Problem in allocation of variable lwnet','','')

    ALLOCATE (fluxsubli(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'mleb_initialize','Problem in allocation of variable fluxsubli','','')

    ALLOCATE (qsat_air(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'mleb_initialize','Problem in allocation of variable qsat_air','','')

    ALLOCATE (tair(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'mleb_initialize','Problem in allocation of variable tair','','')

    ALLOCATE (temp_atmos_pres_grid(kjpindex,jnlvls),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'mleb_initialize','Problem in allocation of variable temp_atmos_pres_grid','','')

    ALLOCATE (q_atmos_pres_grid(kjpindex,jnlvls),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'mleb_initialize','Problem in allocation of variable q_atoms_pres_grid','','')

    ALLOCATE (temp_leaf_pres_grid(kjpindex,jnlvls),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'mleb_initialize','Problem in allocation of variable temp_leaf_pres_grid','','')

    ALLOCATE (temp_surf_pres(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'mleb_initialize','Problem in allocation of variable temp_surf_pres','','')



    !! 2. Initialize variables from restart file or by default values
    !! The variables read are: temp_sol (surface temperature), qsurf (near surface specific humidity),
    !! evapot (soil potential evaporation), evapot_corr (corrected soil potential evaporation), tsolrad
    !! (radiative surface temperature), evapora (evaporation), fluxlat (latent heat flux), fluxsens
    !! (sensible heat flux) and temp_sol_new (new surface temperature).
    IF (printlev>=3) WRITE (numout,*) 'Read a restart file for ENERBIL variables'

    CALL ioconf_setatt_p('UNITS', 'K')
    CALL ioconf_setatt_p('LONG_NAME','Surface temperature')
    CALL restget_p (rest_id, 'temp_sol', nbp_glo, 1, 1, kjit, .TRUE., temp_sol, "gather", nbp_glo, index_g)

    !Config Key   = ENERBIL_TSURF
    !Config Desc  = Initial temperature if not found in restart
    !Config If    = OK_SECHIBA
    !Config Def   = 280.
    !Config Help  = The initial value of surface temperature if its value is not found
    !Config         in the restart file. This should only be used if the model is 
    !Config         started without a restart file.
    !Config Units = Kelvin [K]
    CALL setvar_p (temp_sol, val_exp,'ENERBIL_TSURF', 280._r_std)

    ! Initialize temp_sol_new with temp_sol. These variables are always equal in the beginning of a new time step.
    temp_sol_new(:) = temp_sol(:)

    CALL ioconf_setatt_p('UNITS', 'g/g')
    CALL ioconf_setatt_p('LONG_NAME','near surface specific humidity')
    CALL restget_p (rest_id, 'qsurf', nbp_glo, 1, 1, kjit, .TRUE., qsurf, "gather", nbp_glo, index_g)
    IF ( ALL( qsurf(:) .EQ. val_exp ) ) THEN 
       qsurf(:) = qair(:)
    ENDIF

    CALL ioconf_setatt_p('UNITS', 'mm day^{-1}')
    CALL ioconf_setatt_p('LONG_NAME','Soil Potential Evaporation')
    CALL restget_p (rest_id, 'evapot', nbp_glo, 1, 1, kjit, .TRUE., evapot, "gather", nbp_glo, index_g)
    CALL setvar_p (evapot, val_exp, 'ENERBIL_EVAPOT', zero)

    CALL ioconf_setatt_p('UNITS', 'mm day^{-1}')
    CALL ioconf_setatt_p('LONG_NAME','Corrected Soil Potential Evaporation')
    CALL restget_p (rest_id, 'evapot_corr', nbp_glo, 1, 1, kjit, .TRUE., evapot_corr, "gather", nbp_glo, index_g)
    !Config Key   = ENERBIL_EVAPOT
    !Config Desc  = Initial Soil Potential Evaporation
    !Config If    = OK_SECHIBA       
    !Config Def   = 0.0
    !Config Help  = The initial value of soil potential evaporation if its value 
    !Config         is not found in the restart file. This should only be used if
    !Config         the model is started without a restart file. 
    !Config Units = 
    CALL setvar_p (evapot_corr, val_exp, 'ENERBIL_EVAPOT', zero)

    CALL ioconf_setatt_p('UNITS', 'K')
    CALL ioconf_setatt_p('LONG_NAME','Radiative surface temperature')
    CALL restget_p (rest_id, 'tsolrad', nbp_glo, 1, 1, kjit, .TRUE., tsol_rad, "gather", nbp_glo, index_g)
    IF ( ALL( tsol_rad(:) .EQ. val_exp ) ) THEN 
       tsol_rad(:) = temp_sol(:)
    ENDIF

    !! 1.3 Set the fluxes so that we have something reasonable and not NaN on some machines
    CALL ioconf_setatt_p('UNITS', 'Kg/m^2/dt')
    CALL ioconf_setatt_p('LONG_NAME','Evaporation')
    CALL restget_p (rest_id, 'evapora', nbp_glo, 1, 1, kjit, .TRUE., vevapp, "gather", nbp_glo, index_g)
    IF ( ALL( vevapp(:) .EQ. val_exp ) ) THEN 
       vevapp(:) = zero
    ENDIF

    CALL ioconf_setatt_p('UNITS', 'W/m^2')
    CALL ioconf_setatt_p('LONG_NAME','Latent heat flux')
    CALL restget_p (rest_id, 'fluxlat', nbp_glo, 1, 1, kjit, .TRUE., fluxlat, "gather", nbp_glo, index_g)
    IF ( ALL( fluxlat(:) .EQ. val_exp ) ) THEN 
       fluxlat(:) = zero
    ENDIF

    CALL ioconf_setatt_p('UNITS', 'W/m^2')
    CALL ioconf_setatt_p('LONG_NAME','Sensible heat flux')
    CALL restget_p (rest_id, 'fluxsens', nbp_glo, 1, 1, kjit, .TRUE., fluxsens, "gather", nbp_glo, index_g)
    IF ( ALL( fluxsens(:) .EQ. val_exp ) ) THEN 
       fluxsens(:) = zero
    ENDIF

    !! 1. Initialize arrays
    ALLOCATE (pad_deltah(0:jnlvls),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in pad_deltah allocation. '
    END IF

    ALLOCATE (box_height(jnlvls+1),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in box_height allocation.'
    END IF

    ALLOCATE (box_height_old(jnlvls+1),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in box_height_old allocation.'
    END IF

    ALLOCATE (jtheta(jnlvls),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in jtheta  allocation.'
    END IF

    !! 2. Read parameters from the run.def

    k_ext_fac = 1.0
    CALL getin_p('K_EXT_FAC',k_ext_fac)

    IF (jnlvls .eq. 1) THEN
       ! no read-in profile for single layer case

       ! set a default value of box_height
       box_height(:) = 1.0d0

       ! set a default value of pad_deltah
       pad_deltah(:) = 1.0d0

       k_eddy_slope = 1.0_r_std
       k_eddy_ustar = 0.001_r_std

    ELSE

       !Config Key   = K_EDDY_SLOPE
       !Config Desc  = coefficients for S-shape fuction for K_eddy
       !Config If    = OK_MLEB
       !Config Def   = 10.0
       !Config Help  = 
       !Config Units = [-]
       ! set a default value of k_eddy_slope
       k_eddy_slope = 20.0_r_std
       CALL getin_p('K_EDDY_SLOPE', k_eddy_slope)

       !Config Key   = K_EDDY_USTAR
       !Config Desc  = Parameter for k_eddy function
       !Config If    = OK_MLEB
       !Config Def   = 0.1
       !Config Help  = 
       !Config Units = [-]
       ! set a default value of k_eddy_ustar
       k_eddy_ustar = 0.4_r_std
       CALL getin_p('K_EDDY_USTAR', k_eddy_ustar)
    ENDIF

    !Config Key   = mleb_cdfname
    !Config Desc  = Name of the output CDF file
    !Config If    = OK_MLEB
    !Config Def   = col_cdf.nc
    !Config Help  = This values gives the name of the output CDF file. 
    !Config Units = [-]
    ! set a default value of mleb_cdfname
    mleb_cdfname = 'col_cdf.nc'
    CALL getin_p('mleb_cdfname', mleb_cdfname)

    !Config Key   = LW_EXTINCTION_COEF
    !Config If    = OK_MLEB
    !Config Desc  = Value of the extinction coeffecient
    !Config Def   = 0.8
    !Config Help  = This values gives the value of the LW extinction coefficient.
    !Config Units = [-]
    ! set a default value of big_k_lw (if there is nothing provided in  run.def)
    big_k_lw = 0.75d0 ! first test the result with out this, as it has
    ! not been included for many revesions
    CALL getin_p('LW_EXTINCTION_COEF', big_k_lw)

    !Config Key   = SW_EXTINCTION_COEF
    !Config If    = OK_MLEB
    !Config Desc  = Value of the extinction coeffecient
    !Config Def   = 0.45 
    !Config Help  = This values gives the value of the LW extinction coefficient.
    !Config Units = [-]
    ! set a default value of big_k_sw (if there is nothing provided in run.def)
    big_k_sw = 0.45d0  ! first test the result with out this, as it
    ! has  not been included for many revesions
    CALL getin_p('SW_EXTINCTION_COEF', big_k_sw)

    !Config Key   = KS_TUNE
    !Config Desc  = Coefficient for tuning the forest floor evapotranspiration
    !Config If    = OK_MLEB
    !Config Def   = 1.0
    !Config Help  = 
    !Config Units = [-]
    ! set a default value of ks_tune
    ks_tune = 1.0_r_std
    CALL getin_p('KS_TUNE', ks_tune)

    !Config Key   = KS_SLOPE
    !Config Desc  = Coefficient for surface evapotranspiration
    !Config If    = OK_MLEB
    !Config Def   = 20.0
    !Config Help  = 
    !Config Units = [-]
    ! set a default value of 
    ks_slope = 17.0_r_std
    CALL getin_p('KS_SLOPE', ks_slope)

    !Config Key   = KS_VEGET
    !Config Desc  = Coefficient for surface evapotranspiration
    !Config If    = OK_MLEB
    !Config Def   = 0.35
    !Config Help  = 
    !Config Units = [-]
    ! set a default value of ks_veget
    ks_veget = 0.9_r_std
    CALL getin_p('KS_VEGET', ks_veget)

    !Config Key   = SR_FAC
    !Config If    = OK_MLEB
    !Config Desc  = Value of the stomatal resistance factor
    !Config Def   = 1.0
    !Config Help  = 
    !Config Units = [-]
    ! set a default value of sr_fac
    sr_fac = 2.426_r_std
    CALL getin_p('SR_FAC', sr_fac)

    !Config Key   = BR_FAC
    !Config If    = OK_MLEB
    !Config Desc  = Value of the boundary resistance factor
    !Config Def   = 1.0
    !Config Help  = 
    !Config Units = [-]
    ! set a default value of br_fac
    br_fac = 0.857_r_std
    CALL getin_p('BR_FAC', br_fac)

    !Config Key   = WC02
    !Config Desc  = Logic variable for using vertical dynamic drag coefficient
    !Config If    = OK_MLEB
    !Config Def   = FALSE
    !Config Help  = Logic variable for using vertical dynamic drag coefficient 
    !               based on Wohlfahrt and Cernusca, 2002 
    !Config Units = [FLAG]
    WC02 = .TRUE.
    CALL getin_p('WC02', WC02)

    !Config Key  = A_1 to A_5 
    !Config Desc = Sub-canopy wind spped parametrization
    !Config IF   = OK_MLEB
    !Config Def  = A_1:6.140, A_2:0.00, A_3:0.434, A_4:-0.751, A_5:0.071
    !Config Help = 
    !Config Units = [-]
    a_1 =  6.140d0
    a_2 =  0.001d0
    a_3 =  0.360d0
    a_4 =  -0.081d0
    a_5 =  0.028d0
    CALL getin_p('A_1',a_1)
    CALL getin_p('A_2',a_2)
    CALL getin_p('A_3',a_3)
    CALL getin_p('A_4',a_4)
    CALL getin_p('A_5',a_5)


    !! 3.1 Set values for the mleb calculations

    rho_veg = 1000.0d0              ! assume that density of vegetation is as for pure water for now (1000 kg/m^3)
    leaf_tks = 0.00014d0            ! leaf thickness (derived from TRYdatabase), values tried by JR,  0.00014d0,
    ! 0.001d, 0.0017d0 

    jtheta = 4.18d3                 ! specific heat capacity of vegetation (assumed equal to water J/(kg K))
    !jtheta = jtheta * 0.5          ! see reference of
    ! http://link.springer.com/article/10.1007%2Fs10765-010-0877-7
    ! M. S. Jayalakshmy and  J. Philip  
    ! Thermophysical Properties of Plant Leaves and Their Influence 
    ! on the Environment Temperature, International Journal of 
    ! Thermophysics, 2010, Volume 31, Issue 11-12, pp 2295-2304
    r_sto_v_fac = 1.0d0
    big_r_h2o_fac = 1.0d0


    surf_lev = 0.0d0
    pad_deltah(0) = 0.001d0


    ! 3.2 Set the  mleb_cdfname
    mleb_cdfname = '/home/surface7/jalleon/'//mleb_cdfname

    first_call_mleb_histwrite = .TRUE.

    ! 3.2 Read from the restart files

    ! The leaf temperature profile for the first time step (for the
    ! multi-level energy budget) must have an initial value, which we
    ! also set to the temperature of air (just for this initial step)

    var_name = 'temp_atmos_pres_grid'
    CALL restget_p (rest_id, var_name, nbp_glo, jnlvls, 1, kjit, &
         &        .TRUE., temp_atmos_pres_grid(:,:), 'gather', nbp_glo, index_g)
    IF ( ALL(temp_atmos_pres_grid(:,:) == val_exp) ) THEN
       DO ji = 1,kjpindex
          temp_atmos_pres_grid(ji,:) = temp_air(ji)
       ENDDO
    ENDIF

    q_atmos_pres_grid(:,:) = val_exp
    var_name = 'q_atmos_pres_grid'
    CALL restget_p (rest_id, var_name, nbp_glo, jnlvls , 1, kjit, &
         &        .TRUE., q_atmos_pres_grid(:,:), 'gather', nbp_glo, index_g)
    DO ji = 1,kjpindex
       IF ( ALL(q_atmos_pres_grid(ji,:) == val_exp) ) q_atmos_pres_grid(ji,:) = qair(ji)
    ENDDO

    temp_leaf_pres_grid(:,:) = val_exp
    var_name = 'temp_leaf_pres_grid'
    CALL restget_p (rest_id, var_name, nbp_glo, jnlvls , 1, kjit, &
         &        .TRUE., temp_leaf_pres_grid(:,:), 'gather', nbp_glo, index_g)
    DO ji = 1,kjpindex
       IF ( ALL(temp_leaf_pres_grid(ji,:) == val_exp) ) temp_leaf_pres_grid(ji,:) = temp_air(ji)
    ENDDO

    lwup(:) = val_exp
    var_name = 'lwup'
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, &
         &        .TRUE., lwup(:), 'gather', nbp_glo, index_g)
    IF ( ALL(lwup(:) == val_exp) ) lwup(:) = 0.97 * c_stefan * temp_air(:)**4


    temp_surf_pres(:) = val_exp
    var_name = 'temp_surf_pres'
    CALL restget_p (rest_id, var_name, nbp_glo,1,1, kjit, &
         &        .TRUE., temp_surf_pres, 'gather', nbp_glo, index_g )
    ! WARNING temp_surf_pres has no dimensions, whihc is okay for now since
    ! used for either jnlvls greater than 1 or in loops where where just
    ! overwritten at each kjpindex. In cases that we run for jnlvls greater
    ! than 1, it will only be for one pixel. But when the multi-layer energy
    ! budget possibly will be further refined, dimensions should be added here.
    IF ( ALL(temp_surf_pres(:) == val_exp) )  temp_surf_pres(:) = temp_air(:)


    ! --- CHECK ---!
    ! I (ASL) am doubting whether below code section is needed, because we
    ! already get z_array_out from slowproc_initialize and max_height_store.
    ! This part of code  was originally positioned in sechiba_initialize and
    ! bound be the flag ok_new_enerbil_nextstep. No recent tests have been made
    ! with hack_can_structure. Once this will be done check whether or not
    ! we need the below code.
    ! --- CHECK ---!

    ! The 'hack_can_structure' flag activates the sections of code which directly
    ! link the energy budget scheme to the the size and LAI profile of the canopy for
    ! the respective PFT and age class that is calculated in stomate, for the albedo.
    ! If 'hack_can_structure' is TRUE, and jnlvls > 1, then the model takes LAI
    ! profile information and canopy level heights from the run.def.
    ! If 'hack_can_structure' is FALSE, and jnlvls > 1, then the profile infomation
    ! and canopy levels heights comes from the PGap-based processes for calculation of
    ! stand profile information in stomate.
    ! If jnlvls=1, then the code behaves as for the single layer model.

    IF (.NOT. hack_can_structure) THEN

       ! for enerbil nextstep initialisation

       ! we have a problem in that z_array_out has a value of zero on the first call 
       ! of slowproc_main (above), so we define z_array_out with temporary values.

       u_speed(:,:) = val_exp
       var_name = 'u_speed'
       CALL restget_p (rest_id, var_name, nbp_glo, jnlvls, 1, kjit, &
            &              .TRUE., u_speed, 'gather', nbp_glo, index_g)
       IF (ALL(u_speed(:,:) == val_exp)) u_speed(:,:) = zero
       !--- CHECK ---!
       ! Determine the dimensions of u_speed. In diffuco it is over nlevels_tot,
       ! while jnlvls in mleb and sechiba.f90. Currently, we are not affected by
       ! this because our standard setting is hack_can_structure=true.
       !--- CHECK ---! 

       z_array_out(:,:,:,:) = val_exp
       var_name = 'z_array_out'
       CALL restget_p (rest_id, var_name, nbp_glo, nvm, ncirc, nlevels_tot, kjit, &
            &   .TRUE., z_array_out, 'gather', nbp_glo, index_g)
       IF (ALL(z_array_out(:,:,:,:) ==val_exp)) THEN

          !Config Key   = MAX_HEIGHT_STORE
          !Config Desc  = The initial maximum canopy height
          !Config If    = OK_MLEB
          !Config Def   = 0.0
          !Config Help  = 
          !Config Units = [-]
          max_height_store(:,:) = 20.0d0
          !+++CHECK+++
          ! getin_p can not put a single value into an array
          ! This will not work in parallel. See getin_p forest_managed
          ! for a solution
          CALL getin_p('MAX_HEIGHT_STORE', max_height_store(:,:))
          !++++++++++++

          z_array_out(:,:,:,:) = zero
          DO ilev = 1, nlevels_tot
             z_array_out(:,:,:,ilev) = (REAL(ilev)) * (max_height_store(1,test_pft)/REAL(nlevels_tot))
          END DO ! i = 1, nlevels_tot
       END IF

    END IF ! IF (hack_can_structure)

  END SUBROUTINE mleb_initialize

  !!  =============================================================================================================================
  !! SUBROUTINE		 		    : mleb_clear
  !!
  !>\BRIEF				    Routine deallocates clear output variables if already allocated.
  !!
  !! DESCRIPTION			    : This is a 'housekeeping' routine that deallocates the key output
  !! variables, if they have already been allocated. The variables that are deallocated are psold,
  !! qsol_sat, pdqsold, psnew, qsol_sat_new, lwup, qsat_air, tair
  !!
  !! RECENT CHANGE(S)			    : None
  !!
  !! MAIN OUTPUT VARIABLE(S)	            : None
  !!
  !! REFERENCES				    : None
  !! 
  !! FLOWCHART                              : None
  !! \n
  !_ ==============================================================================================================================

  SUBROUTINE mleb_clear ()
    IF ( ALLOCATED (psold)) DEALLOCATE (psold)
    IF ( ALLOCATED (qsol_sat)) DEALLOCATE (qsol_sat)
    IF ( ALLOCATED (pdqsold)) DEALLOCATE (pdqsold)
    IF ( ALLOCATED (psnew)) DEALLOCATE (psnew)
    IF ( ALLOCATED (qsol_sat_new)) DEALLOCATE (qsol_sat_new)
    IF ( ALLOCATED (lwabs)) DEALLOCATE (lwabs)
    IF ( ALLOCATED (lwup)) DEALLOCATE (lwup)
    IF ( ALLOCATED (lwnet)) DEALLOCATE (lwnet)
    IF ( ALLOCATED (fluxsubli)) DEALLOCATE (fluxsubli)
    IF ( ALLOCATED (qsat_air)) DEALLOCATE (qsat_air)
    IF ( ALLOCATED (tair)) DEALLOCATE (tair)
    IF ( ALLOCATED (temp_atmos_pres_grid)) DEALLOCATE (temp_atmos_pres_grid)
    IF ( ALLOCATED (q_atmos_pres_grid)) DEALLOCATE (q_atmos_pres_grid)
    IF ( ALLOCATED (temp_surf_pres)) DEALLOCATE (temp_surf_pres)
    IF ( ALLOCATED (temp_leaf_pres_grid)) DEALLOCATE (temp_leaf_pres_grid)
    IF ( ALLOCATED (pad_deltah)) DEALLOCATE (pad_deltah)
    IF ( ALLOCATED (box_height)) DEALLOCATE (box_height)
    IF ( ALLOCATED (box_height_old)) DEALLOCATE (box_height_old)
    IF ( ALLOCATED (jtheta)) DEALLOCATE (jtheta)

  END SUBROUTINE mleb_clear


  !! ===========================================================================================================================
  !! SUBROUTINE: mleb_main
  !!
  !!
  !>\BRIEF: calls the main processes for the multi-layer version of the energy budget
  !!          which may also be run in single-layer mode, when 'jnlvls' is set equal
  !!          to 1
  !!
  !! DESCRIPTION: Mleb_main contains a series of subroutine designed to
  !! calculate the energy budget through the canopy. Currently, the full
  !! functionality is only applicable at the pixel scale with 1 PFT in the given
  !! pixel. Thus, as a default the model is currently run with jnlvls = 1, where
  !! jnlvls is the number of levels within the canopy, to be able to run on
  !! larger spatial scales. 
  !! 
  !!
  !! RECENT CHANGE(S): None
  !!
  !! TO DO: In general the module and subroutines could benefit from an overall
  !! cleaning, which includes: removal of commented variables in the declaritions,
  !! check dimensions of variables with vertical dimensions, check variables
  !! with no given dimensions, remove code which is substituded by the
  !! mleb_box_height subroutine, explain all declared variables.
  !!
  !! MAIN OUTPUT VARIABLE(S): vevapnu, vevapsno, vevapflo, transpir, 
  !!                          transpot, vevapwet, 
  !!                          temp_surf_pres, flux_ground_h, 
  !!                          flux_ground_le, 
  !!                          qsol_sat_new_out, qair_new, 
  !!                          Light_Abs_Tot_mean, Light_Alb_Tot_mean
  !!                          vbeta3, evapot, evapot_corr, temp_sol, qsurf
  !!                          fluxsens, fluxlat, netrad, tsol_rad, vevapp
  !!                          temp_sol_new, vbeta2sum, vbeta3sum
  !!
  !! REFERENCE(S): Ryder et al., 2016, GMD, doi:10.5194/gmd-9-223-2016 
  !!
  !! FLOWCHART:
  !_ ==============================================================================================================================

  !+++CHECK+++
  ! Variable dimensions xxx_Tot_mean are for a single pixel. Causing 1+1 problems
  ! when running a larger domain. Needs to be corrected when implementing
  ! a global use of the multi-layer energy budget.

  ! Secondly since the old hydr arch removal transpir_supply is not calculated
  ! anymore even though it is used as E_supply below
  SUBROUTINE mleb_main (kjit, kjpindex, ldrestart_read, ldrestart_write, &
       & index, indexveg, lwdown, swnet, swdown, &
       & epot_air, temp_air, u, v, petAcoef, petBcoef, &
       & qair, peqAcoef, peqBcoef, pb, rau, &
       & vbeta, vbeta1, vbeta2, vbeta3, &
       & vbeta3pot, vbeta4, vbeta5, &
       & emis, soilflx, soilcap, q_cdrag, &
       & humrel, fluxsens, fluxlat, vevapp, transpir, &
       & transpot, vevapnu, vevapwet, vevapsno, &
       & vevapflo, temp_sol, tsol_rad, temp_sol_new, &
       & qsurf, evapot, evapot_corr, rest_id, hist_id, & 
       & hist2_id, &
       & ok_mleb_history_file, &
       & E_supply, hydrol_flag, &
       & hydrol_flag2, hydrol_flag3, vbeta2sum, vbeta3sum, veget_max, &
       & qsol_sat_new_out, qair_new, veget, &
       & Light_Abs_Tot_mean, Light_Alb_Tot_mean, &
       & laieff_isotrop, lai_per_level, z_array_out, &
       & transpir_supply_column, u_speed, profile_vbeta3, profile_rveget, &
       & max_height_store, &
       & precip_rain, pgflux, snowdz, temp_sol_add)
    !+++++++++++

    IMPLICIT NONE



    !! 0 Variable and parameter description 

    !! 0.1 Input variables

    LOGICAL,                                                   INTENT(in)     :: ldrestart_read           !! Logical for _restart_ file to read (-)
    LOGICAL,                                                   INTENT(in)     :: ldrestart_write          !! Logical for _restart_ file to write (-)
    LOGICAL,                                                   INTENT(in)     :: ok_mleb_history_file     !! Flag that controls the netCDF output routine
    LOGICAL,                                                   INTENT(in)     :: hydrol_flag              !! Flag that recalculates the energy budget for the new supply term (true/false)    
    LOGICAL,        DIMENSION(kjpindex),                       INTENT(in)     :: hydrol_flag2             !! flag that 'trips' the alternative energy budget for each grid square
    LOGICAL,        DIMENSION(kjpindex,nvm),                   INTENT(in)     :: hydrol_flag3             !! flag that 'trips' the alternative energy budget for each grid square and PFT

    INTEGER(i_std),                                            INTENT(in)     :: kjit                     !! Time step number (-)
    INTEGER(i_std),                                            INTENT(in)     :: kjpindex                 !! Domain size (-) 

    INTEGER(i_std),                                            INTENT(in)     :: rest_id,hist_id          !! _Restart_ file and _history_ file identifier (-)
    INTEGER(i_std),                                            INTENT(in)     :: hist2_id                 !! _history_ file 2 identifier (-)

    INTEGER(i_std), DIMENSION(kjpindex),                       INTENT(in)     :: index                    !! Indeces of the points on the map (-)
    INTEGER(i_std), DIMENSION(kjpindex*nvm),                   INTENT(in)     :: indexveg                 !! Indeces of the points on the 3D map

    REAL(r_std),    DIMENSION(kjpindex,nvm),                   INTENT(in)     :: E_supply                 !! Supply of water for transpiration (mm timestep^{-1})

    REAL(r_std),    DIMENSION(kjpindex),                       INTENT(in)     :: lwdown                   !! Down-welling long-wave flux (W m^{-2})
    REAL(r_std),    DIMENSION(kjpindex),                       INTENT(in)     :: swnet                    !! Net surface short-wave flux (W m^{-2})
    REAL(r_std),    DIMENSION(kjpindex),                       INTENT(in)     :: swdown                   !! Down-welling short-wave flux (W/m^2)
    REAL(r_std),    DIMENSION(kjpindex),                       INTENT(in)     :: epot_air                 !! Air potential energy (=cp_air*t_air (J))

    REAL(r_std),    DIMENSION(kjpindex),                       INTENT(in)     :: temp_air                 !! Air temperature (K)
    REAL(r_std),    DIMENSION(kjpindex),                       INTENT(in)     :: u                        !! Eastward Lowest level wind speed  (m s^{-1})
    REAL(r_std),    DIMENSION(kjpindex),                       INTENT(in)     :: v                        !! Northward Lowest level wind speed (m s^{-1})
    REAL(r_std),    DIMENSION(kjpindex),                       INTENT(in)     :: petAcoef                 !! Coefficient used for the implicit calculation of the temperature (see technical note)
    REAL(r_std),    DIMENSION(kjpindex),                       INTENT(in)     :: petBcoef                 !! Coefficient used for the implicit calculation of the temperature (see technical note)
    REAL(r_std),    DIMENSION(kjpindex),                       INTENT(in)     :: peqAcoef                 !! Coefficient used for the implicit calculation of the humidity (see technical note)
    REAL(r_std),    DIMENSION(kjpindex),                       INTENT(in)     :: peqBcoef                 !! Coefficient used for the implicit calculation of the humidity (see technical note)
    REAL(r_std),    DIMENSION(kjpindex),                       INTENT(in)     :: pb                       !! Lowest level pressure (hPa)
    REAL(r_std),    DIMENSION(kjpindex),                       INTENT(in)     :: rau                      !! Air density (kg m^{-3})

    REAL(r_std),    DIMENSION(kjpindex),                       INTENT(inout)  :: vbeta                    !! Resistance coefficient to the total evaporation flux (-)
    REAL(r_std),    DIMENSION(kjpindex),                       INTENT(in)     :: vbeta1                   !! Resistance to the snow sublimation flux (-) 
    REAL(r_std),    DIMENSION(kjpindex),                       INTENT(in)     :: vbeta4                   !! Resistance to the bare soil evaporation flux (-)
    REAL(r_std),    DIMENSION(kjpindex),                       INTENT(in)     :: vbeta5                   !! Resistance to th Floodplains evaporation flux (-)
    REAL(r_std),    DIMENSION(kjpindex),                       INTENT(in)     :: emis                     !! Surface emissivity (-)

    REAL(r_std),    DIMENSION(kjpindex),                       INTENT(in)     :: soilflx                  !! Soil heat flux (W m^{-2})
    REAL(r_std),    DIMENSION(kjpindex),                       INTENT(in)     :: soilcap                  !! Soil calorific capacity (J K^{-1])
    REAL(r_std),    DIMENSION(kjpindex),                       INTENT(in)     :: q_cdrag                  !! Surface drag coefficient (-)
!!$    REAL(r_std),    DIMENSION(kjpindex),                       INTENT(in)     :: ccanopy                  !! CO_2 concentration in the canopy (ppm)
    REAL(r_std),    DIMENSION(kjpindex,nvm),                   INTENT(in)     :: humrel                   !! Soil moisture stress coefficient (within range 0 to 1)
    REAL(r_std),    DIMENSION(kjpindex,nvm),                   INTENT(in)     :: vbeta2                   !! Resistance to the intercepted water evaporation (-)
    REAL(r_std),    DIMENSION(kjpindex,nvm),                   INTENT(in)     :: vbeta3pot                !! Resistance to the potential transpiration flux
!!$    REAL(r_std),    DIMENSION(kjpindex,nvm),                   INTENT(in)     :: vbetaco2                 !! Vegetation resistance to CO2 (-)
!!$    REAL(r_std),    DIMENSION(kjpindex,nvm),                   INTENT(in)     :: cimean                   !! mean intercellular ci (from STOMATE) 
                                                                                                             !! (\mu mole m^{-2} s^{-1})

    REAL(r_std),    DIMENSION(kjpindex),                       INTENT(in)     :: qair                     !! Atmospheric lowest level specific humidity (kg kg^{-1})
    REAL(r_std),    DIMENSION(kjpindex,nvm),                   INTENT(in)     :: veget_max                !! Max. fraction of vegetation type (between 0 and 1, LAI -> infty, unitless)
    REAL(r_std),    DIMENSION(kjpindex,nvm),                   INTENT(in)     :: veget                    !! Vegetation fraction for each pft (0-1)
    REAL(r_std),    DIMENSION(kjpindex,nvm,ncirc,nlevels_tot), INTENT(inout)  :: z_array_out              !! Height of the different nodes in the canopy (An output of h_array, to use in sechiba)                             
    REAL(r_std),    DIMENSION(nlevels_tot, kjpindex, nvm),     INTENT(in)     :: transpir_supply_column   !! Supply of water for transpiration at each level in the canopy (mm/tstep) 
    REAL(r_std),    DIMENSION(kjpindex,nvm),                   INTENT(inout)  :: max_height_store         !! Maximum height of the canopy (m)
    REAL(r_std),    DIMENSION(kjpindex,nlevels_tot,nvm),       INTENT(in)     :: laieff_isotrop           !! Effective LAI
    REAL(r_std),    DIMENSION(kjpindex,nvm,nlevels_tot),       INTENT(in)     :: lai_per_level            !! LAI per level
    REAL(r_std),    DIMENSION(kjpindex,nvm,nlevels_tot),       INTENT(in)     :: profile_vbeta3           !! Profile of vbeta3 inside the canopy
    REAL(r_std),    DIMENSION(kjpindex,nvm,nlevels_tot),       INTENT(in)     :: profile_rveget           !! Profile of stomatal resistance in the canopy
    REAL(r_std),    DIMENSION(kjpindex),                       INTENT(in)     :: precip_rain              !! Rainfall
    REAL(r_std),    DIMENSION(kjpindex,nsnow),                 INTENT(in)     :: snowdz                   !! Snow depth at each snow layer [m]


    !! 0.2 Output variables

    REAL(r_std),    DIMENSION(kjpindex),                       INTENT(out)    :: vevapnu                  !! Bare soil evaporation flux (mm day^{-1})            
    REAL(r_std),    DIMENSION(kjpindex),                       INTENT(out)    :: vevapsno                 !! Snow sublimation flux (mm day^{-1})               
    REAL(r_std),    DIMENSION(kjpindex),                       INTENT(out)    :: vevapflo                 !! Floodplains evaporation flux (mm/d)                       
    REAL(r_std),    DIMENSION(kjpindex,nvm),                   INTENT(out)    :: transpir                 !! Transpiration flux (mm day^{-1})                  
    REAL(r_std),    DIMENSION(kjpindex,nvm),                   INTENT(out)    :: transpot                 !! Potential transpiration flux (mm/d)                  
    REAL(r_std),    DIMENSION(kjpindex,nvm),                   INTENT(out)    :: vevapwet                 !! Intercepted water evaporation flux (mm day^{-1})                    
    REAL(r_std),    DIMENSION(kjpindex),                       INTENT(out)    :: qsol_sat_new_out         !! Saturated surface air moisture at the new time step (kg kg^{-1})
    REAL(r_std),    DIMENSION(kjpindex),                       INTENT(out)    :: qair_new                 !! Air moisture at the new time step (kg.kg^(-1))
    REAL(r_std),    DIMENSION(kjpindex),                       INTENT(out)    :: temp_sol_add             !! Additional energy to melt snow for snow ablation case (K) 

    !! 0.3 Modified variables
    !+++CHECK+++
    ! Variable dimensions are for a single pixel. Causing 1+1 problems
    ! when running a larger domain. Needs to be corrected when implementing
    ! a global use of the multi-layer energy budget.
    REAL(r_std),    DIMENSION(nlevels_tot),                    INTENT(in)     :: Light_Abs_Tot_mean       !! Total light absorption for a given canopy level (unitless, 0-1)
    REAL(r_std),    DIMENSION(nlevels_tot),                    INTENT(in)     :: Light_Alb_Tot_mean       !! Total albedo for a given level (unitless, 0-1)
    !+++++++++++
    REAL(r_std),    DIMENSION(kjpindex,nvm),                   INTENT(inout)  :: vbeta3                   !! Resistance to the transpiration flux (-)
    REAL(r_std),    DIMENSION(kjpindex),                       INTENT(inout)  :: evapot                   !! Soil Potential Evaporation flux (mm/tstep)
    REAL(r_std),    DIMENSION(kjpindex),                       INTENT(inout)  :: evapot_corr              !! Soil Potential Evaporation Corrected by Milly's correction (mm/tstep)
    REAL(r_std),    DIMENSION(kjpindex),                       INTENT(inout)  :: temp_sol                 !! Surface temperature (K)
    REAL(r_std),    DIMENSION(kjpindex),                       INTENT(inout)  :: qsurf                    !! Surface specific humidity (kg kg^{-1})
    REAL(r_std),    DIMENSION(kjpindex),                       INTENT(out)    :: fluxsens                 !! Top level sensible heat flux (W m^{-2})
    REAL(r_std),    DIMENSION(kjpindex),                       INTENT(out)    :: fluxlat                  !! Top level latent heat flux (W m^{-2})
    REAL(r_std),    DIMENSION(kjpindex),                       INTENT(out)    :: tsol_rad                 !! Radiative temperature equivalent (W m^{-2})
    REAL(r_std),    DIMENSION(kjpindex),                       INTENT(out)    :: vevapp                   !! Total of evaporation at the top level (mm day^{-1})

    REAL(r_std),    DIMENSION(kjpindex),                       INTENT(out)    :: temp_sol_new             !! Surface temperature at the new time step (K)
    REAL(r_std),    DIMENSION(kjpindex),                       INTENT(inout)  :: vbeta2sum                !! Sum of all vbeta2 according to PFTs (-)
    REAL(r_std),    DIMENSION(kjpindex),                       INTENT(inout)  :: vbeta3sum                !! Sum of all vbeta3 according to PFTs (-)

    REAL(r_std),    DIMENSION(kjpindex),                       INTENT(inout)  :: pgflux                   !! Net energy injected into snowpack(W/m^2)
    REAL(r_std),    DIMENSION(kjpindex,jnlvls),                INTENT(out)    :: u_speed                  !! Wind speed at specified canopy level (m/s)

    !! 0.4 Local variables
    LOGICAL,        DIMENSION(kjpindex)                                       :: warning_correction       !! Warning if Milly's correction denominator is zero over a pixel

    INTEGER(i_std)                           	                              :: i, j                     !! Indexes (levels from soil surface to top of the canopy)
    INTEGER(i_std)                    	                                      :: ji, jv                   !! Indexes (grid-cell pixels, PFTs)
    INTEGER(i_std)                                                            :: ier                      !! Allocation variable (useless)

    REAL(r_std),    DIMENSION(jnlvls)                                         :: horiz_flux_t             !! Flux of sensible heat from leaves to air inside the canopy (W.m^(-2))
    REAL(r_std),    DIMENSION(jnlvls)                                         :: horiz_flux_q             !! Flux of latent heat from leaves to air inside the canopy (W.m^(-2))
    REAL(r_std)                                                               :: u_star                   !! Friction velocity (m/s)
    REAL(r_std),    DIMENSION(0:jnlvls)                                       :: jrpad                    !! Plant Area Density, replacement for values above (as definitions between different
                                                                                                          !!    parts of the model was becoming confusing)
                                                                                                          !! 'jrpad' is the total one-sided leaf area per unit of layer volume
    REAL(r_std)                                                               :: temp_surf_next           !! Surface temperature at the new time step (K) (same as temp_soil_next?)

    REAL(r_std),    DIMENSION(kjpindex)                                       :: netrad                   !! Net radiation flux (W m^{-2})
    REAL(r_std)                                                               :: flux_ground_h            !! Ground level sensible heat flux exchanged between surface and first canopy air layer (W/m^2)
    REAL(r_std)                                                               :: flux_ground_le           !! Ground level latent heat flux exchanged between surface and first canopy air layer (W/m^2)

    REAL(r_std),    DIMENSION(kjpindex,jnlvls)                                :: flux_h_grid              !! Sensible heat flux between each level in the column (W/m^2)
    REAL(r_std),    DIMENSION(kjpindex,jnlvls)                                :: flux_le_grid             !! latent heat flux between each level in the column (W/m^2)
    REAL(r_std),    DIMENSION(kjpindex,jnlvls)                                :: q_a_next_grid            !! Atmospheric specific humidity profile at the new timestep
    REAL(r_std),    DIMENSION(kjpindex,jnlvls)                                :: t_a_next_grid            !! Atmospheric temperature profile at the new timestep
    REAL(r_std),    DIMENSION(jnlvls)                                         :: h_flux, le_flux          !! Sensible and latent heat fluxes at the top of the canopy (W/m^2)
    REAL(r_std),    DIMENSION(jnlvls)                                         :: jalpha                   !! 
    REAL(r_std),    DIMENSION(jnlvls)                                         :: jbeta                    !! 
    REAL(r_std),    DIMENSION(jnlvls)                                         :: k_eddy                   !! Eddy transport coefficient (m^2 s-1)
    REAL(r_std),    DIMENSION(jnlvls)                                         :: down_sw_tab              !! 
    REAL(r_std),    DIMENSION(jnlvls)                                         :: down_lw_tab              !! 
    REAL(r_std),    DIMENSION(jnlvls)                                         :: big_r                    !! 
    REAL(r_std),    DIMENSION(jnlvls)                                         :: big_r_prime              !! 
    REAL(r_std)                                                               :: layer_height             !! 
    REAL(r_std),    DIMENSION(jnlvls)                                         :: delta_h                  !! Height of each level (m)
    REAL(r_std),    DIMENSION(jnlvls)                                         :: delta_z                  !! Vertical difference between the middle of each box (m)
    !   and the middle of the first layer (m)
    REAL(r_std),    DIMENSION(kjpindex)                                       :: jalpha_surf              !!
    REAL(r_std),    DIMENSION(kjpindex)                                       :: jbeta_surf               !! 
    REAL(r_std),    DIMENSION(jnlvls)                                         :: big_e, big_f, big_g      !!
    REAL(r_std),    DIMENSION(jnlvls)                                         :: a_t, b_t, c_t, d_t       !!
    REAL(r_std),    DIMENSION(jnlvls)                                         :: a_q, b_q, c_q, d_q       !!
    REAL(r_std),    DIMENSION(jnlvls)                                         :: jemissivity              !!
    REAL(r_std),    DIMENSION(jnlvls)                                         :: big_pad                  !! plant area density at each level (m^2/m^3)
    REAL(r_std),    DIMENSION(jnlvls)                                         :: eta_1, eta_2             !!
    REAL(r_std),    DIMENSION(jnlvls)                                         :: eta_3, eta_4             !!
    REAL(r_std),    DIMENSION(jnlvls)                                         :: jomega                   !!
    REAL(r_std)                                                               :: eta_1_surf               !!
    REAL(r_std)                                                               :: eta_2_surf               !!
    REAL(r_std)                                                               :: eta_3_surf               !!
    REAL(r_std)                                                               :: eta_4_surf               !!
    REAL(r_std)                                                               :: eta_1_above              !!
    REAL(r_std)                                                               :: eta_2_above              !!
    REAL(r_std)                                                               :: eta_3_above              !!
    REAL(r_std)                                                               :: eta_4_above              !!

    REAL(r_std),    DIMENSION(kjpindex,jnlvls)                                :: temp_leaf_next_grid      !! leaf temperature down the column (K) PRESENT STEP

    REAL(r_std)                                                               :: delta_energy_surf        !!
    REAL(r_std),    DIMENSION(kjpindex,jnlvls)                                :: delta_energy             !! 
    REAL(r_std),    DIMENSION(jnlvls)                                         :: jomega_sw                !!
    REAL(r_std)                                                               :: jomega_surf              !! 
    REAL(r_std)                                                               :: jemissivity_surf         !!
    REAL(r_std)                                                               :: jomega_sw_surf           !! emissivity of surface layer
    REAL(r_std),    DIMENSION(jnlvls+1)                                       :: top_lev                  !!
    REAL(r_std),    DIMENSION(kjpindex)                                       :: q_sat_surf               !!
    REAL(r_std)                                                               :: zikt, zikq, speed        !!
    REAL(r_std)                                                               :: delta_temp               !! change in temperature over the time step (Kelvins)
    REAL(r_std),    DIMENSION(kjpindex)                                       :: jpdqsold                 !!
    REAL(r_std),    DIMENSION(kjpindex)                                       :: epot_air_new             !!
    CHARACTER(LEN=256)                                                        :: var_name                 !! To store variables names for I/O
    REAL(r_std),    DIMENSION(jnlvls)                                         :: t_l_pres, t_l_next       !!
    REAL(r_std),    DIMENSION(0:jnlvls+1, 0:jnlvls+1)                         :: ngu_alpha                !! 
    REAL(r_std),    DIMENSION(kjpindex)                                       :: dev_qsol                 !!
    REAL(r_std),    DIMENSION(kjpindex)                                       :: gs_effective             !! 'effective' stomatal conductance
    REAL(r_std),    DIMENSION(kjpindex)                                       :: vbeta3_effective         !!

    REAL(r_std),    DIMENSION(kjpindex)                                       :: grad_qsat                !!
    REAL(r_std)                                                               :: correction               !!
    REAL(r_std)                                                               :: qc                       !! Surface drag coefficient (??)

    REAL(r_std)                                                               :: light_surf_absorbed
    REAL(r_std),    DIMENSION(jnlvls)                                         :: light_layer_absorbed     !!

    REAL(r_std)                                                               :: pft_running_total        !!
    REAL(r_std),    DIMENSION(0:jnlvls)                                       :: pad_deltah_running_total !!
    REAL(r_std),    DIMENSION(jnlvls+1)                                       :: box_height_running_total !!


    ! each potential input variable is listed separately - I can't see a cleaner way to 
    !   achieve that this using Fortran. This places a limit of 50 on the number of canopy 
    !   levels in the energy budget profile, which is sufficient for the moment.

    REAL(r_std)                                                               :: pad_lev_1, pad_lev_2, pad_lev_3, pad_lev_4 
    REAL(r_std)                                                               :: pad_lev_5, pad_lev_6, pad_lev_7, pad_lev_8 
    REAL(r_std)                                                               :: pad_lev_9, pad_lev_10, pad_lev_11, pad_lev_12 
    REAL(r_std)                                                               :: pad_lev_13, pad_lev_14, pad_lev_15, pad_lev_16 
    REAL(r_std)                                                               :: pad_lev_17, pad_lev_18, pad_lev_19, pad_lev_20
    REAL(r_std)                                                               :: pad_lev_21, pad_lev_22, pad_lev_23, pad_lev_24 
    REAL(r_std)                                                               :: pad_lev_25, pad_lev_26, pad_lev_27, pad_lev_28 
    REAL(r_std)                                                               :: pad_lev_29, pad_lev_30, pad_lev_31, pad_lev_32 
    REAL(r_std)                                                               :: pad_lev_33, pad_lev_34, pad_lev_35, pad_lev_36 
    REAL(r_std)                                                               :: pad_lev_37, pad_lev_38, pad_lev_39, pad_lev_40
    REAL(r_std)                                                               :: pad_lev_41, pad_lev_42, pad_lev_43, pad_lev_44 
    REAL(r_std)                                                               :: pad_lev_45, pad_lev_46, pad_lev_47, pad_lev_48 
    REAL(r_std)                                                               :: pad_lev_49, pad_lev_50    

    REAL(r_std)                                                               :: height_lev_1, height_lev_2, height_lev_3, height_lev_4 
    REAL(r_std)                                                               :: height_lev_5, height_lev_6, height_lev_7, height_lev_8 
    REAL(r_std)                                                               :: height_lev_9, height_lev_10, height_lev_11, height_lev_12 
    REAL(r_std)                                                               :: height_lev_13, height_lev_14, height_lev_15, height_lev_16 
    REAL(r_std)                                                               :: height_lev_17, height_lev_18, height_lev_19, height_lev_20
    REAL(r_std)                                                               :: height_lev_21, height_lev_22, height_lev_23, height_lev_24 
    REAL(r_std)                                                               :: height_lev_25, height_lev_26, height_lev_27, height_lev_28 
    REAL(r_std)                                                               :: height_lev_29, height_lev_30, height_lev_31, height_lev_32 
    REAL(r_std)                                                               :: height_lev_33, height_lev_34, height_lev_35, height_lev_36 
    REAL(r_std)                                                               :: height_lev_37, height_lev_38, height_lev_39, height_lev_40
    REAL(r_std)                                                               :: height_lev_41, height_lev_42, height_lev_43, height_lev_44 
    REAL(r_std)                                                               :: height_lev_45, height_lev_46, height_lev_47, height_lev_48 
    REAL(r_std)                                                               :: height_lev_49, height_lev_50

    CHARACTER(len=256)                                                        :: jname
    CHARACTER(len=256)                                                        :: jname2
    CHARACTER(len=256)                                                        :: format_string


    REAL(r_std),    DIMENSION(0:jnlvls+1)                                     :: ngu_bign
    REAL(r_std),    DIMENSION(0:jnlvls+1)                                     :: ngu_alpha_insert


    REAL(r_std),    DIMENSION(jnlvls)                                         :: flux_netrad       ! net radiation flux for each layer
    REAL(r_std),    DIMENSION(jnlvls)                                         :: flux_netsw        ! net radiation contributed from short-wave
    REAL(r_std),    DIMENSION(jnlvls)                                         :: flux_netlw        ! net radiation contributed from long-wave

    LOGICAL                                                                   :: new_sub=.FALSE. !temporarry for development of new sub routine

    !_ =======================================================================================


    !! This main routine calculated the multi-layer energy budget. The steps are the following:
    !!     - Define the discretisation all along the canopy (two options: mleb_boxheight (to test) or the following code);
    !!     - Define the lai profile in each layer of the canopy with the routine mleb_profile;
    !!     - Solve the Long-Wave Radiation budget in each layer with mleb_lwrad;
    !!     - Calculate the boundary layer and stomatal resistances in each layer with mleb_boundarylayer_resistance
    !!       and mleb_stomatal_resistance;
    !!     - Calculate the alpha and beta coefficients necessary for the Short Wave Radiation budget with mleb_alpha_beta_coef;
    !!     - Solve the energy budget in each layer thanks to mleb_calc_column;
    !!     - Calculate the fluxes emanating from the forest;
    !!     - Split the water flux in its components with mleb_evapveg.
    !!
    !! The different layers are defined as follows:
    !!     - Soil surface:                                                 Layer:   1                                  ;
    !!     - From soil surface to first layer of canopy:                   Layers:  1                to  jnlvls_under  ;
    !!     - Inside the canopy:                                            Layers:  jnlvls_under+1   to  jnlvls_canopy ;
    !!     - From the top of the canopy to the first layer of atmospher:   Layers:  jnlvls_canopy+1  to  jnlvls        .


    !! Debug:

    IF (printlev_loc>=4) THEN
       WRITE (numout,*) '240215 start of enerbil, tair is: ', tair
       WRITE (numout,*) '061116 in mleb_main, epot_air is: ', epot_air 
       WRITE (numout,*) '20150327 lwdown in mleb_main:', lwdown(:)
       WRITE (numout,*) 'start of mleb_main vbeta',vbeta(:)
       WRITE (numout,*) 'hydrol_flag', hydrol_flag
    END IF ! (printlev_loc>=4) THEN


    WRITE (numout,*) 'MLEB'
    WRITE (numout,*) 'hydrol_flag', hydrol_flag


    ! for cases in which we exit the loop if hydrol_flag2 is false
    epot_air_new = epot_air
    temp_sol_new = temp_sol

    ! DEBUG, check laieff_isotrop before pad and bx_height calc.
    IF (printlev_loc>=4) THEN
       DO ji = 1,1
          !DO jv = 4,4
          DO jv=1,nvm
             DO i = (jnlvls_under + 1) , (jnlvls_under + jnlvls_canopy)
                WRITE (numout, *) 'before ns laieff_isotrop:',ji, i-jnlvls_under+1, jv, &
                     &'is', laieff_isotrop(ji, i-jnlvls_under, jv)
             END DO ! i 
          END DO ! jv
       END DO ! ji
       WRITE(numout,*)'kjit',kjit 
    END IF ! (printlev_loc>=4)
    !




    !! TO DO: To test the new subroutine vs the original code (the idea is to move the previous code 
    !! and have a proper code
    
    IF(new_sub) THEN

       ! call the new subroutine - make this more precise
       CALL mleb_boxheight(kjpindex,laieff_isotrop, z_array_out,max_height_store, &
            & rau,temp_atmos_pres_grid,q_atmos_pres_grid,temp_leaf_pres_grid, &
            & box_height,top_lev, pad_deltah, jrpad)


    ELSE
       ! Use the original code
       ! ------------------------------------------------------------------------------------
       !
       ! there are two ways to arrive at the MLA (mass per unit leaf area, kg/m^2), 
       ! which we then multiply by the LAD (leaf area density, m^2/m^3) to arrive at the 
       ! vegetation matter per m in height

       ! 1) as in Nobel, Park S. (2005) 'Physiochemical and Environmental Plant Physiology' 
       ! Edition 3, published Elsevier, chapter 7.1, p. 309
       !    - assume leaf has a specific heat of water (4.19 kJ kg^-1 degC^-1 at 20 degC)
       !    - assume leaf has thickness 300micrometres
       !    - assume leaf has density of 700kg/m^3 (water is 1000kg/m^3 but leaves are often 
       !        30% air)
       !    - this gives (300 x 10^-6 m)(700 kg m^-3) = 0.21 kg m^-2

       ! 2) Using LMA (leaf per unit mass area, m^2/kg) from the TRY database
       !    - for Eucalyptus delegatensis we have 155 cm^2/g ---> 0.155 m^2/g. 
       !      This is dried weight of the leaf.
       !    - 100 gC/m^2 soil
       !    - 100 * 2 * 10 = 2000 g C / m^2 soil
       !    Dry weight / Fresh weight = 8 / SLA ---> g FW / m^2 (here we are assuming a dry 
       !                                               weight to fresh weight ratio of 1:9)
       !                              = 8 / 0.155
       !                              = 516 g FW / m^2
       !                              = 0.516 kg / m^2
       ! 
       ! ------------------------------------------------------------------------------------

       !! Debug

       IF (printlev_loc>=4) THEN  
          WRITE (numout, *) 'SUM(laieff_isotrop(1:jnlvls_canopy, test_grid,test_pft)) is: ', &
               &  SUM(laieff_isotrop(test_grid,1:jnlvls_canopy,test_pft)) 
       END IF !(printlev_loc>=4)

       IF (jnlvls .ne. 1) THEN
          IF (jnlvls .ne. (jnlvls_under + jnlvls_canopy + jnlvls_over)) THEN
             WRITE (numout, *) 'number of levels in canopy column imbalance'
             STOP
          END IF
       END IF

       IF (jnlvls .ne. 1) THEN
          IF (nlai .ne. (jnlvls_canopy + 1)) THEN
             WRITE (numout, *) 'number of levels in canopy column imbalance - nlai should be equal to jnlvls_canopy+1'
             STOP
          END IF
       END IF

       !! The definition of the LAI imposes the height of each layer and the LAI inside each layer.
       !! The code presented here only concerns the case where LAI is not imposed


       IF (.NOT. hack_can_structure) THEN 

          pft_running_total = 0.0d0
          box_height_running_total = 0.0d0
          pad_deltah_running_total = 0.0d0

          !! TO DO: Do we keep the loop over the PFTs or we let mleb at the PFT level?

          DO ji = 1,1
             ! Currently, the mleb only runs for one PFT, and therefore
             ! jv is set to the PFT that is used instead of 1,nvm. 
             DO jv = test_pft,test_pft 
                !DO jv = 1,nvm

                !! TO DO: Clean the definition of max_height_store and z_array_out and see if the IF conditions
                !!        are still mandatory.

                !! Those IF conditions should be used to better intialise the values. max_height_store and z_array_out
                !! should not be changed here. This would lead to inconsistencies.

!!$                IF (veget_max(ji,jv) .GT. zero .AND. veget(ji,jv) .GT. zero) THEN 
!!$
!!$                   IF ((SUM(laieff_isotrop(ji, 1:jnlvls_canopy, jv)) .gt. min_sechiba) .AND. &
!!$                        & (SUM(z_array_out(ji, jv, 1, :)) .gt. min_sechiba)) THEN
!!$
!!$                      IF (SUM(z_array_out(ji, jv, 1, :)) .gt. min_sechiba) THEN

                !! Debug ID1.
                !! Quick patch to make the code run. Should be initialized properly

                !! If there are no leaves, the layers are initialized thanks to the max height of the forest.
                !! This should be changed by considering local variables for max_height_store and z_array_out.
 
                IF (ANY(z_array_out(ji,jv,1,2:) .LE. min_sechiba)) THEN
                   IF (ANY(max_height_store(:,test_pft) .LE. min_sechiba)) THEN
                      max_height_store(:,:)=20.
                   ENDIF
                   DO i = 1, nlevels_tot
                      z_array_out(:,:,:,i) = (REAL(i)) * (max_height_store(1,test_pft)/REAL(nlevels_tot))
                   END DO ! i = 1, nlevels_tot

                END IF
                
                !! A faire plus proprement

                IF(read_lai .EQ. 1) THEN
                   max_height_store(test_grid,test_pft) = max_canopy_height_prescribed
                ENDIF
                !! End ID1
                

                !! The plant area density in each layer of canopy (from jnlvls_under+1 to jnlvls_canopy)
                !! is updated considering the effective LAI from stomate

                pad_deltah(:) = 0.001d0

                !! NB: The first level of lai_eff_isotrop represents the understorey.
                !!     We set this equal to the min.

                DO i = (jnlvls_under + 1) , (jnlvls_under + jnlvls_canopy+1)

                   pad_deltah(i) = lai_per_level(ji,jv,i-jnlvls_under)
!!$                   pad_deltah(i) = laieff_isotrop(ji, i-jnlvls_under+1, jv) * cos(pi/3.0d0)

                   IF (pad_deltah(i) .lt. 1.0d-3) THEN

                      pad_deltah(i) = 1.0d-3

                   END IF

                   !! Debug

                   IF (printlev_loc>=4) THEN
                      WRITE (numout, *) '251016 laieff_isotrop:',ji, i-jnlvls_under+1, jv, &
                           'is', laieff_isotrop(ji, i-jnlvls_under, jv)
                      WRITE (numout, *) '160315 pad_deltah level: ', i,' is: ', pad_deltah(i)
                   END IF ! (printlev_loc>=4)


                END DO ! i = (jnlvls_under + 1) , (jnlvls_under + jnlvls_canopy)





                !! Debug

                IF (printlev_loc>=4) THEN
                   WRITE (numout, *) 'in calculate box_height section,pft',jv
                   WRITE (numout, *) 'max_height_store(ji,jv) is: ', max_height_store(ji, jv)
                   WRITE (numout, *) 'z_array_out(ji,jv,1,:) is: ', z_array_out(ji, jv, 1, :)
                   WRITE(numout,*) "pad_deltah", pad_deltah(:)
                END IF ! (printlev_loc>=4) THEN


                !! Definition of the box_heights from the soil surface to the first layer of the atmosphere.

                !! Under the canopy:
                !! Z_array_out corresponds to the height of the layers inside the canopy. The first value corresponds
                !! to the soil surface. The difference between the two first values corresponds thus to the under-canopy
                !! space. This space is devided into jnlvls_under boxes of similar width

                DO i = 1, jnlvls_under - 1

                   box_height(i) = (i) * (z_array_out(ji, jv, 1, 1)) / (jnlvls_under)

                END DO ! 1, jnlvls_under - 1 

                !! The height of the under-canopy interface corresponds to the second value of z_array_out.
                
                box_height(jnlvls_under) = z_array_out(ji, jv, 1, 1)


                !! Debug ID2.
                !! Review indexation. Should be reviewed a bit to be proper.

                !! Inside the canopy:
                !! Inside the canopy, the height of each box corresponds precisely to the height of each canopy layer 
                !! represented by z_array_out.

                DO i = (jnlvls_under + 1), (jnlvls_under + jnlvls_canopy - 1)

                   box_height(i) = z_array_out(ji, jv, 1, i - jnlvls_under + 1)

                END DO ! i = (jnlvls_under + 1), (jnlvls_under + jnlvls_canopy - 1)

                !! At the interface between the cnaopy and the atmosphere, the height corresponds to the top of the 
                !! canopy. The MAX value is here to avoid bugs when there are nos leaves. +0.1 is here to avoid having
                !! a last canopy layer with a 0 width. This should be double-checked.

                box_height(jnlvls_under + jnlvls_canopy) = MAX(box_height(jnlvls_under + jnlvls_canopy-1)+0.1, max_height_store(ji, jv))


                !! Above the canopy:
                !! The first layer of the atmosphere is set at a height of 50m. The difference between 50m and the 
                !! height of the canopy is then divided into jnlvls_over boxes of similar width.

                DO i = (jnlvls_under + jnlvls_canopy + 1), (jnlvls_under + jnlvls_canopy + jnlvls_over - 1)

                   box_height(i) = MAX(box_height(jnlvls_under + jnlvls_canopy), max_height_store(ji, jv)) + &
                        & (( i - (jnlvls_under + jnlvls_canopy) ) * &
                        & (49.0d0 - MAX(box_height(jnlvls_under + jnlvls_canopy), max_height_store(ji, jv)) )/ jnlvls_over )

                END DO ! i = (jnlvls_under + jnlvls_canopy + 1), (jnlvls_under + jnlvls_canopy + jnlvls_over)

                !! End ID2

                !! The two last layers are fixed

                box_height(jnlvls_under + jnlvls_canopy + jnlvls_over) = 49.0d0
                box_height(jnlvls_under + jnlvls_canopy + jnlvls_over+1) = 50.0d0

                jrpad(:) = pad_deltah(:)

                !! Debug

                IF (printlev_loc>=4) THEN
                   DO i = 1, jnlvls
                      WRITE(numout, *) 'lev: ', i, ' b_h: ', box_height(i), ' p_d: ', pad_deltah(i)
                   END DO ! i = 1, jnlvls
                END IF ! (printlev_loc>=4) THEN

                ! END IF ! (SUM(laieff_isotrop(ji, 1:jnlvls_canopy, jv)) .ne. 0.0)


                !! TO DO: Check if this is still needed

                ! IF (SUM(laieff_isotrop(ji, 1:jnlvls_canopy, jv)) .gt. min_sechiba) THEN
                !  pft_running_total = pft_running_total + 1
                ! END IF

                !  box_height_running_total = box_height_running_total + box_height ! NOT needed for one PFT
                !  pad_deltah_running_total = pad_deltah_running_total + pad_deltah ! NOT needed for one PFT 

                ! END IF !(veget_max(ji,jv) .GT. zero .AND. veget(ji,jv) .GT. zero) THEN 

             END DO ! jv = 2, nvm

                IF (printlev_loc>=4) THEN
                   WRITE(numout,*) "max_height", max_height_store(ji,:)
                   WRITE(numout,*) "box_height", box_height(:)
                ENDIF
             !! Debug

             IF (printlev_loc>=4) THEN
                WRITE(numout, *) '100515 pft_running_total is: ', pft_running_total
                WRITE(numout, *) '100515 box_height_running_total(:) is: ', box_height_running_total(:)
                WRITE(numout, *) '100515 pad_deltah_running_total(:) is: ', pad_deltah_running_total(:)
             END IF ! (printlev_loc>=4) THEN

             
             !! TO DO: Check if this is still needed

             ! IF (pft_running_total .gt. min_sechiba) THEN
             !    box_height = box_height_running_total / pft_running_total
             !    pad_deltah = pad_deltah_running_total / pft_running_total
             ! ELSE
             !    box_height = 0.0d0
             !    pad_deltah = 0.0d0
             ! ENDIF

             !  box_height = box_height_running_total / 1.0d0  ! NOT needed for one PFT
             !  pad_deltah = pad_deltah_running_total / 1.0d0  ! NOT needed for one PFT

             !! TO DO: To remove

             box_height_old = box_height


             !! Debug

             IF(printlev_loc >= 4)THEN
                WRITE(numout,*) 'box_height_old', box_height_old(:)
                WRITE(numout,*) 'box_height', box_height(:)
                WRITE(numout,*) 'temp_atmos_pres_grid',temp_atmos_pres_grid(ji,:)
                WRITE(numout,*) 'q_atmos_pres_grid',q_atmos_pres_grid(ji,:)
                WRITE(numout,*) 'temp_leaf_pres_grid',temp_leaf_pres_grid(ji,:)
             ENDIF ! (printlev_loc >= 4) 


             !! TO DO: To remove

             box_height_old = box_height

          END DO ! ji = 1, kjpindex

       END IF ! (hack_can_structure) THEN 

       top_lev(:) = box_height(:)

    END IF ! new_sub

    ! end the new subruotine here. do not know whether top_lev should included
    ! or not, or if we need if at all. Maybe just substitute top_lev with
    ! box_height (aslanso).

    ! calculate wind speed based on input velocities
    speed = MAX(min_wind, SQRT (u(1)*u(1) + v(1)*v(1)))! something odd here,missing dimension

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>


    !! Here comes the resolution of the multi-layer energy budget.


    !! TO DO: Surely useless, to remove.

    ! as per the original scheme
    CALL mleb_begin (kjpindex, temp_sol, lwdown, swnet, pb, &
         & psold, qsol_sat, pdqsold, netrad, emis, &
         & dev_qsol, lwabs)


    !! Calculates the transport profile inside the canopy (LAI and turbulences)

    IF (jnlvls .eq. 1) THEN
       !! No transport profile for single layer case
       delta_h(:) = 1.0d0
       delta_z(:) = 1.0d0
       u_speed(:,:) = zero
    ELSE
       CALL mleb_profile (pad_deltah, layer_height, delta_z, &
            & delta_h, box_height, k_eddy, &
            & u_star, speed, u_speed, big_pad, &
            & k_eddy_slope,k_eddy_ustar, surf_lev, top_lev, swnet(1), jrpad, & 
            & a_1, a_2, a_3, a_4, a_5, WC02)
    END IF !(jnlvls .eq. 1)    

    !! Calculates the Long Wave Radiations profile inside the canopy following Gu et al. (1999)

    IF (jnlvls .eq. 1) THEN
       !! No radiation profile for single layer case
       eta_1(:)=zero
       eta_2(:)=zero
       eta_3(:)=zero
       eta_4(:)=zero
       jomega(:)=zero
       jomega_sw(:)=zero
       jemissivity(:)=zero
       jemissivity_surf=zero
       eta_1_surf=zero
       eta_2_surf=zero
       eta_3_surf=zero
       eta_4_surf=zero
       ngu_alpha(:,:)=zero
       jomega_surf=zero
       jomega_sw_surf=zero
       eta_1_above=zero
       eta_2_above=zero
       eta_3_above=zero
       eta_4_above=zero
    ELSE
       CALL mleb_lwrad (jomega, jemissivity, temp_leaf_pres_grid, temp_surf_pres, &
            & pad_deltah, jemissivity_surf, jomega_surf, jomega_sw_surf, &
            & eta_3_above, eta_1_above, ngu_alpha, ngu_alpha_insert, kjpindex, &
            & big_k_lw, big_k_sw) 
    END IF !(jnlvls .eq. 1)    


    !! Calculation of the boundary layer resistance to heat transport in each layer. This resistance will be applied
    !! to the transport of sensible heat flux between the leaves and the atmosphere inside the canopy.

    IF (jnlvls .eq. 1) THEN
       !! Use assigned beta coefficients for single layer case
    ELSE
       CALL mleb_boundarylayer_resistance (big_r, u_speed, temp_atmos_pres_grid, &
            & kjpindex, br_fac, q_cdrag, pad_deltah, &
            & jrpad)
    END IF !(jnlvls .eq. 1) 

    !! Calculation of the total resistance to water transport in each layer (boundary layer+stomatal resistances). 
    !! This resistance will be applied to the transport of latent heat flux between the leaves and the atmosphere
    !! inside the canopy.

    IF (jnlvls .eq. 1) THEN
       !! Use assigned beta coefficients for single layer case
    ELSE
       CALL mleb_stomatal_resistance (pad_deltah, swnet, lwdown, pb, & 
            & down_sw_tab, down_lw_tab, big_r_prime, &
            & temp_atmos_pres_grid, q_atmos_pres_grid, big_r, &
            & u_speed, r_sto_v_fac, big_r_h2o_fac, rau, kjpindex, &
            & sr_fac, q_cdrag, &
            & delta_h, jrpad, profile_rveget)
    END IF !(jnlvls .eq. 1)        


    !! Calculation of the alpha and beta coefficients needed for the solving of the Short Wave Radiations
    !! scheme in each layers.

    CALL mleb_alpha_beta_coeff (temp_atmos_pres_grid, temp_leaf_pres_grid, &
         & jalpha, jbeta, jalpha_surf, jbeta_surf, & 
         & temp_surf_pres, q_sat_surf, &
         & kjpindex, pb, temp_sol, &
         & jpdqsold)


    !! Solving of the multi-layer energy budget.

    !! Debug
    IF(printlev_loc>=4) WRITE(numout,*)'before mleb_calc_column, epot_air',epot_air

    IF(printlev_loc>=4) WRITE(numout,*) "laieff_isotrop", laieff_isotrop(1,:,test_pft)

    !+++CHECK+++
    ! Variable dimensions xxx_Tot_mean are for a single pixel. Causing 1+1 problems
    ! when running a larger domain. Needs to be corrected when implementing
    ! a global use of the multi-layer energy budget.
    CALL mleb_calc_column (pad_deltah, layer_height, big_r, big_r_prime, delta_z, &
         & k_eddy, swnet, swdown, lwdown, &
         & temp_atmos_pres_grid, q_atmos_pres_grid, & 
         & a_t, b_t, c_t, d_t, a_q, b_q, c_q, d_q, & 
         & big_e, big_f, big_g, jalpha, jbeta, temp_leaf_pres_grid, delta_h, &
         & eta_1, eta_2, eta_3, eta_4, jtheta, &
         & jalpha_surf, jbeta_surf, rau, &
         & soilcap, soilflx, speed, q_cdrag, kjpindex, &
         & flux_ground_h, flux_ground_le, flux_h_grid, flux_le_grid, &
         & u_speed, t_a_next_grid, q_a_next_grid, & 
         & temp_surf_next, eta_1_surf, eta_2_surf, eta_3_surf, eta_4_surf, &
         & temp_leaf_next_grid, epot_air_new, delta_temp, &
         & horiz_flux_t, horiz_flux_q, &
         & ngu_bign, ngu_alpha, ngu_alpha_insert, qair_new, qair, &
         & emis, temp_sol_new, vbeta1, vbeta4, vbeta5, vbeta2, &
         & vbeta3, vbeta3pot, vbeta, fluxsens, fluxlat, netrad, &
         & petAcoef, petBcoef, peqAcoef, peqBcoef, &
         & u, v, E_supply, hydrol_flag, hydrol_flag2, &
         & hydrol_flag3, vbeta2sum, vbeta3sum, veget_max, veget, &
         & Light_Abs_Tot_mean, Light_Alb_Tot_mean, &
         & light_surf_absorbed, light_layer_absorbed, &
         & transpir_supply_column, pb, profile_vbeta3, profile_rveget, &
         & flux_netrad, flux_netsw, flux_netlw, &
         & ks_tune, ks_slope, ks_veget,qsol_sat_new, transpir)
    !+++++++++++

    !! Debug

    IF (printlev_loc>= 4) THEN
       WRITE(numout,*)'before mleb_flux, epot_air',epot_air
       WRITE(numout,*)'before mleb_flux, epot_air',epot_air_new
    ENDIF !printlev_loc >= 4

    
    !! Calculation of the total fluxes needed in other parts of the code (same as enerbil)

    CALL mleb_flux (kjpindex, emis, temp_sol, rau, u, v, &
         & q_cdrag, vbeta, vbeta1, vbeta5, &
         & qair_new, epot_air_new, psnew, qsurf, &
         & fluxsens , fluxlat , fluxsubli, vevapp, &
         & temp_sol_new, lwdown, swnet, &
         & lwup, lwnet, pb, tsol_rad, netrad, &
         & evapot, evapot_corr, &
         & precip_rain, snowdz, temp_air, pgflux, soilcap, temp_sol_add,qsol_sat_new)


    !! Diagnoses the values for evaporation and transpiration: 
    !!    - vevapsno (snow evaporation);
    !!    - vevapnu (bare soil evaporation); 
    !!    - transpir (transpiration); 
    !!    - gpp (assimilation); 
    !!    - vevapwet (interception).

    CALL mleb_evapveg (kjpindex,  vbeta1, vbeta2, vbeta3, &
         & vbeta3pot, vbeta4, vbeta5, &
         & rau, u, v, q_cdrag, qair_new, humrel, &
         & vevapsno, vevapnu, vevapflo, vevapwet, &
         & transpir, transpot, evapot,qsol_sat_new)


    !! TO DO: to remove
    ! surface layer check
    !ASLA delta_energy_surf = (temp_surf_next - temp_surf_pres) * soilcap(1)

    
    !! Routine to have an output file with all the variables of the multi-layer energy budget.

    !! TO DO: we should use the xios_4dim process and get rid of this.

    IF (jnlvls .eq. 1) THEN
       ! rely on history file output single layer case
    ELSE
       ! only open and write to the mleb history file, if the
       ! ok_mleb_history_file flag is true
       !IF(ok_mleb_history_file) THEN
       IF(.NOT. hydrol_flag) THEN
          IF (mleb_cdfname .NE. '/home/surface7/jalleon/col_cdf.nc') THEN
             CALL mleb_netcdf(first_call_mleb_histwrite, t_a_next_grid(1,:), q_a_next_grid(1,:), &
                  & temp_leaf_next_grid(1,:), &
                  & flux_h_grid(1,:), flux_le_grid(1,:), flux_ground_h, flux_ground_le, &
                  & k_eddy, u_speed, u_star, temp_surf_next, &
                  & jalpha, jbeta, jalpha_surf(1), jbeta_surf(1), & 
                  & eta_1, eta_2, eta_3, eta_4, jomega, &
                  & eta_1_surf, eta_2_surf, eta_3_surf, eta_4_surf, &
                  & lwdown(1), swdown(1), &
                  & temp_atmos_pres_grid(1,:), q_atmos_pres_grid(1,:), &
                  & big_e, big_f, big_g, a_t, b_t, c_t, d_t, a_q, b_q, c_q, d_q, soilflx(1), &
                  & jomega_surf, jemissivity_surf, jomega_sw_surf, jemissivity, &
                  & jomega_sw, temp_leaf_pres_grid(1,:), temp_surf_pres, &
                  & big_r, big_r_prime, horiz_flux_t,  horiz_flux_q, &
                  & soilcap(1), eta_3_above, eta_1_above, eta_4_above, &
                  & eta_2_above,delta_z, delta_h, pad_deltah, mleb_cdfname, &
                  & ngu_bign, ngu_alpha(1,:), ngu_alpha_insert, qsol_sat_new(1), kjit, &
                  & light_surf_absorbed, light_layer_absorbed, box_height, profile_rveget(test_grid,test_pft,:)) 
          ENDIF
       ENDIF
       !ELSE
       ! do nothing
       !ENDIF ! (ok_mleb_history_file)
    END IF !(jnlvls .eq. 1)        


    !! TO DO: to remove
    ! Write to history files
!!$     CALL histwrite_p(hist_id, 'evapot', kjit, evapot, kjpindex, index)
!!$     CALL histwrite_p(hist_id, 'evapot_corr', kjit, evapot_corr, kjpindex, index)
!!$     CALL histwrite_p(hist2_id, 'transpot', kjit, transpot, kjpindex, index)


    !! Saving of the needed variables for the next time step

    temp_atmos_pres_grid(:,:) = t_a_next_grid(:,:)
    q_atmos_pres_grid(:,:) = q_a_next_grid(:,:)
    temp_leaf_pres_grid(:,:) = temp_leaf_next_grid(:,:)

    temp_surf_pres(:) = temp_surf_next

    qsol_sat_new_out(:) = qsol_sat_new

  END SUBROUTINE mleb_main



  !!  ===========================================================================================================================
  !! SUBROUTINE: mleb_netcdf
  !!
  !!
  !>\BRIEF: produces an output netCDF file for the multilayer scheme 
  !!           (though it significantly slows the runtime, so this 
  !!           subroutine is not for general use) 
  !!
  !! DESCRIPTION:
  !! 
  !!
  !! RECENT CHANGE(S): None
  !!
  !! MAIN OUTPUT VARIABLE(S): none
  !!
  !! REFERENCE(S):
  !!
  !! FLOWCHART:
  !_ ==============================================================================================================================


  SUBROUTINE mleb_netcdf(first_call_mleb_histwrite, t_a_next_grid, q_a_next_grid, &
       & temp_leaf_next_grid, &
       & flux_h, flux_le, flux_ground_h, flux_ground_le, &
       & k_eddy, u_speed, u_star, temp_surf_next, &
       & jalpha, jbeta, jalpha_surf, jbeta_surf, & 
       & eta_1, eta_2, eta_3, eta_4, jomega, &
       & eta_1_surf, eta_2_surf, eta_3_surf, eta_4_surf, &
       & down_lw_ac, down_sw_ac, &
       & temp_atmos_pres, q_atmos_pres,big_e, big_f, big_g, a_t, b_t, c_t, d_t, a_q, b_q, c_q, d_q, heat_soil, &
       & jomega_surf, jemissivity_surf, jomega_sw_surf, jemissivity, & 
       & jomega_sw, temp_leaf_pres, temp_surf_pres, &
       & big_r, big_r_prime, horiz_flux_h, horiz_flux_le, &
       & theta_zero, eta_3_above, eta_1_above, eta_4_above, &
       & eta_2_above, delta_z, delta_h, pad_deltah, mleb_cdfname, &
       & ngu_bign, ngu_alpha, ngu_alpha_insert, qsol_sat_new, kjit, &
       & light_surf_absorbed, light_layer_absorbed, box_height, profile_rveget)  

    USE netcdf

    IMPLICIT NONE

    !! 0.1 Input variables
    LOGICAL                                              :: first_call_mleb_histwrite ! detemine whether we open and write headers in
    ! the mleb history file

    CHARACTER(LEN=80), INTENT(in)                        :: mleb_cdfname          

    REAL(r_std), DIMENSION (jnlvls), INTENT(in)          :: jemissivity       
    REAL(r_std), DIMENSION (jnlvls), INTENT(in)          :: temp_atmos_pres   
    REAL(r_std), DIMENSION(jnlvls), INTENT(in)          :: big_e, big_f, big_g
    REAL(r_std), DIMENSION(jnlvls), INTENT(in)          :: a_t, b_t, c_t, d_t
    REAL(r_std), DIMENSION(jnlvls), INTENT(in)          :: a_q, b_q, c_q, d_q
    REAL(r_std), DIMENSION (jnlvls), INTENT(in)          :: q_atmos_pres     
    REAL(r_std), DIMENSION (jnlvls), INTENT(in)          :: big_r, big_r_prime  
    REAL(r_std), DIMENSION (jnlvls), INTENT(in)          :: jalpha, jbeta  
    REAL(r_std),                     INTENT(in)          :: jalpha_surf, jbeta_surf  
    REAL(r_std), DIMENSION (jnlvls), INTENT(in)          :: t_a_next_grid     
    REAL(r_std), DIMENSION (jnlvls), INTENT(in)          :: temp_leaf_next_grid, temp_leaf_pres         ! leaf temperature down the column (K) PRESENT STEP
    REAL(r_std), DIMENSION (jnlvls), INTENT(in)          :: q_a_next_grid        ! atmospheric specific humididy down the column (kg/kg) PRESENT STEP
    REAL(r_std), DIMENSION (jnlvls), INTENT(in)          :: flux_h, flux_le     ! leaf temperature down the column (K) PRESENT STEP
    REAL(r_std), DIMENSION (jnlvls) , INTENT(in)         :: horiz_flux_h,  horiz_flux_le        
    REAL(r_std), DIMENSION (jnlvls), INTENT(in)          :: k_eddy              ! eddy transport coefficient (m^2 s-1)
    REAL(r_std), DIMENSION (:,:), INTENT(in)          :: u_speed             ! (m/s)
    REAL(r_std), INTENT(in)                           :: u_star             ! (m/s)
    REAL(r_std), DIMENSION (jnlvls), INTENT(in)          :: jomega_sw
    REAL(r_std), DIMENSION (jnlvls), INTENT(in)          :: eta_1, eta_2, eta_3, eta_4, jomega

    REAL(r_std), DIMENSION (jnlvls+1), INTENT(in)        :: box_height

    REAL, DIMENSION (0:jnlvls+1), INTENT(in)             :: ngu_alpha
    REAL(r_std), DIMENSION (0:jnlvls+1), INTENT(in)      :: ngu_alpha_insert

    REAL(r_std), DIMENSION (jnlvls), INTENT(in)        :: delta_z         ! box altitude 
    REAL(r_std), DIMENSION (jnlvls), INTENT(in)        :: delta_h         ! box height 
    REAL(r_std), DIMENSION (0:jnlvls), INTENT(in)        :: pad_deltah         ! eddy transport coefficient (m^2 s-1)    
    REAL(r_std), DIMENSION (0:jnlvls+1), INTENT(in)      :: ngu_bign

    REAL(r_std), INTENT(in)                              :: flux_ground_h, flux_ground_le
    REAL(r_std), INTENT(in)                              :: jomega_surf   
    REAL(r_std), INTENT(in)                              :: jemissivity_surf    
    REAL(r_std), INTENT(in)                              :: jomega_sw_surf    
    REAL(r_std), INTENT(in)                              :: temp_surf_next, qsol_sat_new    
    REAL(r_std), DIMENSION (:),INTENT(in)                :: temp_surf_pres
    REAL(r_std), INTENT(in)                              :: eta_1_surf, eta_2_surf
    REAL(r_std), INTENT(in)                              :: eta_3_surf, eta_4_surf, eta_3_above
    REAL(r_std), INTENT(in)                              :: eta_1_above, eta_4_above, eta_2_above   
    REAL(r_std), INTENT(in)                              :: down_lw_ac             ! input downwelling lw radiation above the canopy (watts/m^2)
    REAL(r_std), INTENT(in)                              :: down_sw_ac             ! input downwelling sw radiation above the canopy (watts/m^2) 
    REAL(r_std), INTENT(in)                              :: heat_soil              ! heat flux from or to deeper soil (layered scheme in ORCHIDEE) (W m^-2)
    REAL(r_std), INTENT(in)                              :: theta_zero             ! heat capacity of the infinitesimal surface layer (J/ (K m^2))

    INTEGER(i_std), INTENT(in)                           :: kjit

    REAL(r_std), INTENT(in)                              :: light_surf_absorbed
    REAL(r_std), DIMENSION(:), INTENT(in)                :: light_layer_absorbed
    REAL(r_std),DIMENSION (nlevels_tot), INTENT(in)  :: profile_rveget

    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std), SAVE                     :: tap_id, qap_id, tlp_id, qlp_id              ! ID integers
!$OMP THREADPRIVATE(tap_id, qap_id, tlp_id, qlp_id)
    INTEGER(i_std), SAVE                     :: tan_id, qan_id, tln_id, qln_id              ! ID integers
!$OMP THREADPRIVATE(tan_id, qan_id, tln_id, qln_id)
    INTEGER(i_std), SAVE                     :: bige_id, bigf_id, bigg_id                   ! ID integers
!$OMP THREADPRIVATE(bige_id, bigf_id, bigg_id)
    INTEGER(i_std), SAVE                     :: at_id, bt_id, ct_id, dt_id                  ! ID integers
!$OMP THREADPRIVATE(at_id, bt_id, ct_id, dt_id)
    INTEGER(i_std), SAVE                     :: aq_id, bq_id, cq_id, dq_id                  ! ID integers
!$OMP THREADPRIVATE(aq_id, bq_id, cq_id, dq_id)
    INTEGER(i_std), SAVE                     :: bigr_id, bigrprime_id                       ! ID integers
!$OMP THREADPRIVATE(bigr_id, bigrprime_id)
    INTEGER(i_std), SAVE                     :: horiz_f_h_id, horiz_f_le_id                 ! ID integers
!$OMP THREADPRIVATE(horiz_f_h_id, horiz_f_le_id)
    INTEGER(i_std), SAVE                     :: hflux_id, leflux_id, keddy_id, uspeed_id    ! ID integers
!$OMP THREADPRIVATE(hflux_id, leflux_id, keddy_id, uspeed_id)
    INTEGER(i_std), SAVE                     :: tsn_id, tsp_id, qsn_id                              ! ID integers
!$OMP THREADPRIVATE(tsn_id, tsp_id, qsn_id)
    INTEGER(i_std), SAVE                     :: e1_id, e2_id, e3_id, e4_id                  ! ID integers
!$OMP THREADPRIVATE(e1_id, e2_id, e3_id, e4_id)
    INTEGER(i_std), SAVE                     :: e1s_id, e2s_id, e3s_id, e4s_id              ! ID integers
!$OMP THREADPRIVATE(e1s_id, e2s_id, e3s_id, e4s_id)
    INTEGER(i_std), SAVE                     :: e3a_id, e1a_id, e2a_id, e4a_id              ! ID integers
!$OMP THREADPRIVATE(e3a_id, e1a_id, e2a_id, e4a_id)
    INTEGER(i_std), SAVE                     :: lw_id, sw_id, hs_id                         ! ID integers
!$OMP THREADPRIVATE(lw_id, sw_id, hs_id)
    INTEGER(i_std), SAVE                     :: jos_id, jes_id, josws_id                    ! ID integers
!$OMP THREADPRIVATE(jos_id, jes_id, josws_id)
    INTEGER(i_std), SAVE                     :: jas_id, jbs_id, ja_id, jb_id                ! ID integers
!$OMP THREADPRIVATE(jas_id, jbs_id, ja_id, jb_id)
    INTEGER(i_std), SAVE                     :: qleaf_id, qsurf_id                          ! ID integers
!$OMP THREADPRIVATE(qleaf_id, qsurf_id)
    INTEGER(i_std), SAVE                     :: rveg_id                                       ! ID integers
!$OMP THREADPRIVATE(rveg_id)
    INTEGER(i_std), SAVE                     :: je_id, jo_id, josw_id                       ! ID integers
!$OMP THREADPRIVATE(je_id, jo_id, josw_id)
    INTEGER(i_std), SAVE                     :: tz_id, dz_id, dh_id, pad_id                 ! ID integers  
!$OMP THREADPRIVATE(tz_id, pad_id)
    INTEGER(i_std), SAVE                     :: ngu_bign_id, ngu_bign_s_id, ngu_bign_a_id   ! ID integers
!$OMP THREADPRIVATE(ngu_bign_id, ngu_bign_s_id, ngu_bign_a_id)
    INTEGER(i_std), SAVE                     :: ngu_alpha_id, ngu_alpha_s_id, ngu_alpha_a_id   ! ID integers
!$OMP THREADPRIVATE(ngu_alpha_id, ngu_alpha_s_id, ngu_alpha_a_id)
    INTEGER(i_std), SAVE                     :: ngu_alphai_id, ngu_alphai_s_id, ngu_alphai_a_id   ! ID integers
!$OMP THREADPRIVATE(ngu_alphai_id, ngu_alphai_s_id, ngu_alphai_a_id)
    INTEGER(i_std), SAVE                     :: msa_id, mla_id                            ! ID integers
!$OMP THREADPRIVATE(msa_id, mla_id)
    INTEGER(i_std), SAVE                     :: mleb_step_count    
!$OMP THREADPRIVATE(mleb_step_count)
    INTEGER(i_std), SAVE                     :: iret, iret2, fid, fid2, fid3
!$OMP THREADPRIVATE(iret, iret2, fid, fid2, fid3)
    INTEGER(i_std), SAVE                     :: ibh_id
!$OMP THREADPRIVATE(ibh_id)
    INTEGER(i_std), SAVE                     :: ustar_id
!$OMP THREADPRIVATE(ustar_id)

    INTEGER(i_std)                           :: ii, jj, nlonid, nlatid, varid
    INTEGER(i_std)                           :: TimeDimID, jtid, ts_id, HeightDimID, i
    INTEGER(i_std)                           :: height_dim_id, height_dim_id2, time_dim_id
    INTEGER(i_std)                           :: jstartnc_b, jcountnc_b
    INTEGER(i_std)                           :: dimids_b
    INTEGER(i_std)                           :: dimids(2) 
    INTEGER(i_std)                           :: dimids2(2) 
    INTEGER, DIMENSION(2)                    :: jstartnc, jcountnc    
    LOGICAL                                  :: ok_mleb_history_file    

    REAL(r_std), DIMENSION(jnlvls)           :: flux_h_combo, flux_le_combo
    REAL(r_std), DIMENSION (jnlvls)          :: array_jomega_surf
    REAL(r_std), DIMENSION (jnlvls)          :: array_jemissivity_surf 
    REAL(r_std), DIMENSION (jnlvls)          :: array_jomega_sw_surf 
    REAL(r_std), DIMENSION (jnlvls)          :: array_temp_surf_pres
    REAL(r_std), DIMENSION (jnlvls)          :: array_temp_surf_next
    REAL(r_std), DIMENSION (jnlvls)          :: array_qsol_sat_new
    REAL(r_std), DIMENSION (jnlvls)          :: array_theta_zero
    REAL(r_std), DIMENSION (jnlvls)          :: array_u_star
    REAL(r_std), DIMENSION (jnlvls)          :: array_eta_1_surf, array_eta_2_surf
    REAL(r_std), DIMENSION (jnlvls)          :: array_eta_3_surf, array_eta_4_surf 
    REAL(r_std), DIMENSION (jnlvls)          :: array_eta_3_above, array_eta_1_above
    REAL(r_std), DIMENSION (jnlvls)          :: array_eta_2_above, array_eta_4_above
    REAL(r_std), DIMENSION (jnlvls)          :: array_down_lw_ac, array_down_sw_ac
    REAL(r_std), DIMENSION (jnlvls)          :: array_heat_soil
    REAL(r_std), DIMENSION (jnlvls)          :: array_ngu_bign_surf
    REAL(r_std), DIMENSION (jnlvls)          :: array_ngu_bign_above
    REAL(r_std), DIMENSION (jnlvls)          :: array_ngu_bign
    REAL(r_std), DIMENSION (jnlvls)          :: array_ngu_alpha_surf
    REAL(r_std), DIMENSION (jnlvls)          :: array_ngu_alpha_above
    REAL(r_std), DIMENSION (jnlvls)          :: array_ngu_alpha

    REAL(r_std), DIMENSION (jnlvls)          :: array_ngu_alphai_surf
    REAL(r_std), DIMENSION (jnlvls)          :: array_ngu_alphai_above
    REAL(r_std), DIMENSION (jnlvls)          :: array_ngu_alphai

    REAL(r_std), DIMENSION (jnlvls)          :: array_light_surf_absorbed
    REAL(r_std), DIMENSION (jnlvls)          :: array_light_layer_absorbed
    REAL(r_std), DIMENSION (jnlvls)          :: array_q_surf
    REAL(r_std), DIMENSION (jnlvls)          :: array_q_leaf
    REAL(r_std), DIMENSION (jnlvls)          :: array_jalpha_surf
    REAL(r_std), DIMENSION (jnlvls)          :: array_jbeta_surf
    REAL(r_std), DIMENSION (jnlvls)          :: array_rveget
    !_ ==============================================================================================================================

    REAL(r_std), DIMENSION (jnlvls)          :: intermediate_box_height

    ! parameters, dimensions and variable set-ups (only runs once)

    IF (first_call_mleb_histwrite) THEN

       first_call_mleb_histwrite = .FALSE.

       ! … define the parameters here

       iret = nf90_create(mleb_cdfname, nf90_write, fid)


       ! define the dimensions

       iret = nf90_def_dim(fid, "height", jnlvls, height_dim_id)
       iret = nf90_def_dim(fid, "time", nf90_unlimited, time_dim_id)

       dimids = (/ height_dim_id, time_dim_id /)


       ! define the variables
       iret = nf90_def_var(fid,"temp_atmos_pres",nf90_real4, dimids, tap_id)
       iret = nf90_def_var(fid,"q_atmos_pres",nf90_real4, dimids, qap_id)

       iret = nf90_def_var(fid,"t_a_next_grid",nf90_real4, dimids, tan_id)
       iret = nf90_def_var(fid,"q_a_next_grid",nf90_real4, dimids, qan_id)

       iret = nf90_def_var(fid,"temp_leaf_pres",nf90_real4, dimids, tlp_id)

       iret = nf90_def_var(fid,"temp_leaf_next_grid",nf90_real4, dimids, tln_id)
       iret = nf90_def_var(fid,"q_leaf",nf90_real4, dimids, qleaf_id)
       iret = nf90_def_var(fid,"q_surf",nf90_real4, dimids, qsurf_id)
       iret = nf90_def_var(fid,"jalpha_surf",nf90_real4, dimids, jas_id)
       iret = nf90_def_var(fid,"jbeta_surf",nf90_real4, dimids,jbs_id)
       iret = nf90_def_var(fid,"jalpha",nf90_real4, dimids, ja_id)
       iret = nf90_def_var(fid,"jbeta",nf90_real4, dimids, jb_id)
       iret = nf90_def_var(fid,"rveget",nf90_real4, dimids, rveg_id)

       iret = nf90_def_var(fid,"big_r",nf90_real4, dimids, bigr_id)
       iret = nf90_def_var(fid,"big_r_prime",nf90_real4, dimids, bigrprime_id)

       iret = nf90_def_var(fid,"h_flux",nf90_real4, dimids, hflux_id)
       iret = nf90_def_var(fid,"le_flux",nf90_real4, dimids, leflux_id)

       iret = nf90_def_var(fid,"horiz_flux_h",nf90_real4, dimids, horiz_f_h_id)
       iret = nf90_def_var(fid,"horiz_flux_le",nf90_real4, dimids, horiz_f_le_id)

       iret = nf90_def_var(fid,"k_eddy",nf90_real4, dimids, keddy_id)

       iret = nf90_def_var(fid,"delta_z",nf90_real4, dimids, dz_id)
       iret = nf90_def_var(fid,"delta_h",nf90_real4, dimids, dh_id)
       iret = nf90_def_var(fid,"pad_deltah",nf90_real4, dimids, pad_id)

       iret = nf90_def_var(fid,"u_speed",nf90_real4, dimids, uspeed_id)
       iret = nf90_def_var(fid,"u_star",nf90_real4, dimids, ustar_id)
!!$
       iret = nf90_def_var(fid,"big_e",nf90_real4, dimids, bige_id)
       iret = nf90_def_var(fid,"big_f",nf90_real4, dimids, bigf_id)
       iret = nf90_def_var(fid,"big_g",nf90_real4, dimids, bigg_id)
       iret = nf90_def_var(fid,"a_t",nf90_real4, dimids, at_id)
       iret = nf90_def_var(fid,"b_t",nf90_real4, dimids, bt_id)
       iret = nf90_def_var(fid,"c_t",nf90_real4, dimids, ct_id)
       iret = nf90_def_var(fid,"d_t",nf90_real4, dimids, dt_id)
       iret = nf90_def_var(fid,"a_q",nf90_real4, dimids, aq_id)
       iret = nf90_def_var(fid,"b_q",nf90_real4, dimids, bq_id)
       iret = nf90_def_var(fid,"c_q",nf90_real4, dimids, cq_id)
       iret = nf90_def_var(fid,"d_q",nf90_real4, dimids, dq_id)

       iret = nf90_def_var(fid,"temp_surf_pres",nf90_real4, dimids, tsp_id)
       iret = nf90_def_var(fid,"temp_surf_next",nf90_real4, dimids, tsn_id)
!!$       iret = nf90_def_var(fid,"qsol_sat_new",nf90_real4, dimids, qsn_id)

       iret = nf90_def_var(fid,"qsol_sat_new",nf90_real4, dimids, qsn_id)
       iret = nf90_def_var(fid,"theta_zero",nf90_real4, dimids, tz_id)
!!$
       iret = nf90_def_var(fid,"eta_1",nf90_real4, dimids, e1_id)
       iret = nf90_def_var(fid,"eta_2",nf90_real4, dimids, e2_id)
       iret = nf90_def_var(fid,"eta_3",nf90_real4, dimids, e3_id)
       iret = nf90_def_var(fid,"eta_4",nf90_real4, dimids, e4_id)

       iret = nf90_def_var(fid,"eta_1_surf",nf90_real4, dimids, e1s_id)
       iret = nf90_def_var(fid,"eta_2_surf",nf90_real4, dimids, e2s_id)
       iret = nf90_def_var(fid,"eta_3_surf",nf90_real4, dimids, e3s_id)
       iret = nf90_def_var(fid,"eta_4_surf",nf90_real4, dimids, e4s_id)

       iret = nf90_def_var(fid,"eta_1_above",nf90_real4, dimids, e1a_id)
       iret = nf90_def_var(fid,"eta_2_above",nf90_real4, dimids, e2a_id)
       iret = nf90_def_var(fid,"eta_3_above",nf90_real4, dimids, e3a_id)
       iret = nf90_def_var(fid,"eta_4_above",nf90_real4, dimids, e4a_id)

       iret = nf90_def_var(fid,"lw_down_ac",nf90_real4, dimids, lw_id)
       iret = nf90_def_var(fid,"sw_down_ac",nf90_real4, dimids, sw_id)

       iret = nf90_def_var(fid,"heat_soil",nf90_real4, dimids, hs_id)
!!$
!!$       iret = nf90_def_var(fid,"jomega",nf90_real4, dimids, jo_id)
       !iret = nf90_def_var(fid,"jomega_sw",nf90_real4, dimids, josw_id)
!!$       iret = nf90_def_var(fid,"jomega_surf",nf90_real4, dimids, jos_id)
!!$       iret = nf90_def_var(fid,"jemissivity",nf90_real4, dimids, je_id)
!!$       iret = nf90_def_var(fid,"jemissivity_surf",nf90_real4, dimids, jes_id)
!!$       iret = nf90_def_var(fid,"jomega_sw_surf",nf90_real4, dimids, josws_id)
!!$
!!$       iret = nf90_def_var(fid,"ngu_bign",nf90_real4, dimids, ngu_bign_id)
!!$       iret = nf90_def_var(fid,"ngu_bign_surf",nf90_real4, dimids, ngu_bign_s_id)
!!$       iret = nf90_def_var(fid,"ngu_bign_above",nf90_real4, dimids, ngu_bign_a_id)
!!$       iret = nf90_def_var(fid,"ngu_alpha",nf90_real4, dimids, ngu_alpha_id)
!!$       iret = nf90_def_var(fid,"ngu_alpha_surf",nf90_real4, dimids, ngu_alpha_s_id)
!!$       iret = nf90_def_var(fid,"ngu_alpha_above",nf90_real4, dimids, ngu_alpha_a_id)
!!$
!!$       iret = nf90_def_var(fid,"ngu_alphai",nf90_real4, dimids, ngu_alphai_id)
!!$       iret = nf90_def_var(fid,"ngu_alphai_surf",nf90_real4, dimids, ngu_alphai_s_id)
!!$       iret = nf90_def_var(fid,"ngu_alphai_above",nf90_real4, dimids, ngu_alphai_a_id)

       iret = nf90_def_var(fid,"light_surf_absorbed",nf90_real4, dimids, msa_id)
       iret = nf90_def_var(fid,"light_layer_absorbed",nf90_real4, dimids, mla_id)

       iret = nf90_def_var(fid,"intermediate_box_height",nf90_real4, dimids, ibh_id)

       iret = nf90_enddef (fid)

    END IF !(first_call_mleb_histwrite)

    intermediate_box_height(1:jnlvls) = box_height(1:jnlvls)

    ! open the netCDF file at each step
    iret = nf90_open(mleb_cdfname, nf90_Write, fid) 

    jstartnc = (/ 1, 1 /)  ! vector of integers specifying the index in the variable from 
    !   which the first of the data values will be read
    jcountnc = (/ jnlvls, 1 /)  ! vector of integers specifying the number of indices 
    !   selected along each dimension

    mleb_step_count = mleb_step_count + 1    ! this is the time step counter

    jstartnc(2) = mleb_step_count   ! (2) refers to the time dimension at present (it 
    !   must always come last)

    flux_h_combo(1) = flux_ground_h
    flux_h_combo(2:jnlvls) = flux_h(2:jnlvls)

    flux_le_combo(1) = flux_ground_le
    flux_le_combo(2:jnlvls) = flux_le(2:jnlvls) 

    array_temp_surf_pres(1) = temp_surf_pres(1) ! OBS on dimensions here, if
    !ever running wiht jnlvls greater than 1 and 
    ! kjpindex greater than 1, fix this!
    array_temp_surf_pres(2:jnlvls) = 0.0d0

    array_temp_surf_next(1) = temp_surf_next
    array_temp_surf_next(2:jnlvls) = 0.0d0

    array_q_leaf(:) =  jalpha(:)*temp_leaf_next_grid(:) + jbeta(:)

    array_q_surf(1) = jalpha_surf*temp_surf_next + jbeta_surf
    array_q_surf(2:jnlvls) = 0.0d0

    array_jbeta_surf(1) = jbeta_surf
    array_jbeta_surf(2:jnlvls) = 0.0d0

    array_jalpha_surf(1) = jalpha_surf
    array_jalpha_surf(2:jnlvls) = 0.0d0

    array_rveget(1:10) = 0.0d0
    array_rveget(11:20) = profile_rveget(1:10)
    array_rveget(20:29)=0.0d0

    array_qsol_sat_new(1) = qsol_sat_new
    array_qsol_sat_new(2:jnlvls) = 0.0d0

    array_eta_1_surf(1) = eta_1_surf
    array_eta_1_surf(2:jnlvls) = 0.0d0

    array_eta_2_surf(1) = eta_2_surf
    array_eta_2_surf(2:jnlvls) = 0.0d0    

    array_eta_3_surf(1) = eta_3_surf
    array_eta_3_surf(2:jnlvls) = 0.0d0    

    array_eta_1_above(1) = 0.0d0!eta_1_above
    array_eta_1_above(2:jnlvls) = 0.0d0  

    array_eta_2_above(1) = 0.0d0!eta_2_above
    array_eta_2_above(2:jnlvls) = 0.0d0  

    array_eta_3_above(1) = 0.0d0!eta_3_above
    array_eta_3_above(2:jnlvls) = 0.0d0  

    array_eta_4_above(1) = 0.0d0!eta_4_above
    array_eta_4_above(2:jnlvls) = 0.0d0  

    array_eta_4_surf(1) = eta_4_surf
    array_eta_4_surf(2:jnlvls) = 0.0d0

    array_down_lw_ac(1) = eta_2_surf * temp_surf_next - eta_3_surf!down_lw_ac
    array_down_lw_ac(2:jnlvls) = 0.0d0 
    array_down_sw_ac(1) = eta_4_surf * down_sw_ac
    array_down_sw_ac(2:jnlvls) = 0.0d0

    array_heat_soil(1) = heat_soil
    array_heat_soil(2:jnlvls) = 0.0d0
!!$
!!$    array_jomega_surf(1) = jomega_surf
!!$    array_jomega_surf(2:jnlvls) = 0.0d0
!!$
    array_theta_zero(1) = theta_zero
    array_theta_zero(2:jnlvls) = 0.0d0
!!$
    array_u_star(1) = u_star
    array_u_star(2:jnlvls) = 0.0d0

!!$    array_jemissivity_surf(1) = jemissivity_surf
!!$    array_jemissivity_surf(2:jnlvls) = 0.0d0
!!$
!!$    array_jomega_sw_surf(1) = jomega_sw_surf    
!!$    array_jomega_sw_surf(2:jnlvls) = 0.0d0     
!!$
!!$    array_ngu_bign_surf(1) = ngu_bign(0)    
!!$    array_ngu_bign_surf(2:jnlvls) = 0.0d0            
!!$
!!$    array_ngu_bign_above(1) = ngu_bign(jnlvls+1)
!!$    array_ngu_bign_above(2:jnlvls) = 0.0d0 
!!$
!!$    array_ngu_bign(1:jnlvls) = ngu_bign(1:jnlvls)    
!!$
!!$
!!$
!!$    array_ngu_alpha_surf(1) = ngu_alpha(0)    
!!$    array_ngu_alpha_surf(2:jnlvls) = 0.0d0            
!!$
!!$    array_ngu_alpha_above(1) = ngu_alpha(jnlvls+1)
!!$    array_ngu_alpha_above(2:jnlvls) = 0.0d0 
!!$
!!$    array_ngu_alpha(1:jnlvls) = ngu_alpha(1:jnlvls) 
!!$
!!$
!!$
!!$    array_ngu_alphai_surf(1) = ngu_alpha_insert(0)    
!!$    array_ngu_alphai_surf(2:jnlvls) = 0.0d0            
!!$
!!$    array_ngu_alphai_above(1) = ngu_alpha_insert(jnlvls+1)
!!$    array_ngu_alphai_above(2:jnlvls) = 0.0d0 
!!$
!!$    array_ngu_alphai(1:jnlvls) = ngu_alpha_insert(1:jnlvls) 

    array_light_surf_absorbed(1) = light_surf_absorbed
    array_light_surf_absorbed(2:jnlvls) = 0.0d0 

    array_light_layer_absorbed(1:jnlvls) = light_layer_absorbed(1:jnlvls) 

    IF (printlev_loc>=4) THEN
       WRITE(numout, *) '240215 in mleb_netcdf, after column allocation'
    END IF ! (printlev_loc>=4) THEN


    iret = nf90_put_var(fid, bige_id, big_e(:), start = jstartnc, count = jcountnc) 
    iret = nf90_put_var(fid, bigf_id, big_f(:), start = jstartnc, count = jcountnc) 
    iret = nf90_put_var(fid, bigg_id, big_g(:), start = jstartnc, count = jcountnc) 
    iret = nf90_put_var(fid, at_id, a_t(:), start = jstartnc, count = jcountnc) 
    iret = nf90_put_var(fid, bt_id, b_t(:), start = jstartnc, count = jcountnc) 
    iret = nf90_put_var(fid, ct_id, c_t(:), start = jstartnc, count = jcountnc) 
    iret = nf90_put_var(fid, dt_id, d_t(:), start = jstartnc, count = jcountnc) 
    iret = nf90_put_var(fid, aq_id, a_q(:), start = jstartnc, count = jcountnc) 
    iret = nf90_put_var(fid, bq_id, b_q(:), start = jstartnc, count = jcountnc) 
    iret = nf90_put_var(fid, cq_id, c_q(:), start = jstartnc, count = jcountnc) 
    iret = nf90_put_var(fid, dq_id, d_q(:), start = jstartnc, count = jcountnc)

    iret = nf90_put_var(fid, tap_id, temp_atmos_pres(:), start = jstartnc, count = jcountnc) 
    iret = nf90_put_var(fid, qap_id, q_atmos_pres(:), start = jstartnc, count = jcountnc)

    iret = nf90_put_var(fid, tan_id, t_a_next_grid(:), start = jstartnc, count = jcountnc) 
    iret = nf90_put_var(fid, qan_id, q_a_next_grid(:), start = jstartnc, count = jcountnc)
    iret = nf90_put_var(fid, qleaf_id, array_q_leaf(:), start = jstartnc, count = jcountnc)  
    iret = nf90_put_var(fid, qsurf_id, array_q_surf(:), start = jstartnc, count = jcountnc)     
    iret = nf90_put_var(fid, ja_id, jalpha(:), start = jstartnc, count = jcountnc)     
    iret = nf90_put_var(fid, jb_id, jbeta(:), start = jstartnc, count = jcountnc)     
    iret = nf90_put_var(fid, jas_id, array_jalpha_surf(:), start = jstartnc, count = jcountnc)     
    iret = nf90_put_var(fid, jbs_id, array_jbeta_surf(:), start = jstartnc, count = jcountnc)      
    iret = nf90_put_var(fid, rveg_id, array_rveget(:), start = jstartnc, count = jcountnc)         

    iret = nf90_put_var(fid, qsn_id, array_qsol_sat_new(:), start = jstartnc, count = jcountnc)         

    iret = nf90_put_var(fid, tlp_id, temp_leaf_pres(:), start = jstartnc, count = jcountnc)         

    iret = nf90_put_var(fid, tln_id, temp_leaf_next_grid(:), start = jstartnc, count = jcountnc)         

    iret = nf90_put_var(fid, bigr_id, big_r(:), start = jstartnc, count = jcountnc)         
    iret = nf90_put_var(fid, bigrprime_id, big_r_prime(:), start = jstartnc, count = jcountnc)

    iret = nf90_put_var(fid, hflux_id, flux_h_combo(:), start = jstartnc, count = jcountnc)         
    iret = nf90_put_var(fid, leflux_id, flux_le_combo(:), start = jstartnc, count = jcountnc)         

    iret = nf90_put_var(fid, horiz_f_h_id, horiz_flux_h(:), start = jstartnc, count = jcountnc)         
    iret = nf90_put_var(fid, horiz_f_le_id, horiz_flux_le(:), start = jstartnc, count = jcountnc)

    iret = nf90_put_var(fid, keddy_id, k_eddy(:), start = jstartnc, count = jcountnc)              

    iret = nf90_put_var(fid, dz_id, delta_z(:), start = jstartnc, count = jcountnc) 
    iret = nf90_put_var(fid, dh_id, delta_h(:), start = jstartnc, count = jcountnc) 
    iret = nf90_put_var(fid, pad_id, pad_deltah(1:jnlvls), start = jstartnc, count = jcountnc)              

    iret = nf90_put_var(fid, uspeed_id, u_speed(:,:), start = jstartnc, count = jcountnc) 
    iret = nf90_put_var(fid, ustar_id, array_u_star(:), start = jstartnc, count = jcountnc)              

    iret = nf90_put_var(fid, e1_id, eta_1(:), start = jstartnc, count = jcountnc)              
    iret = nf90_put_var(fid, e2_id, eta_2(:), start = jstartnc, count = jcountnc)              
    iret = nf90_put_var(fid, e3_id, eta_3(:), start = jstartnc, count = jcountnc)              
    iret = nf90_put_var(fid, e4_id, eta_4(:), start = jstartnc, count = jcountnc)              

    iret = nf90_put_var(fid, tsp_id, array_temp_surf_pres(:), start = jstartnc, count = jcountnc)              

    iret = nf90_put_var(fid, tsn_id, array_temp_surf_next(:), start = jstartnc, count = jcountnc)              

    iret = nf90_put_var(fid, e1s_id, array_eta_1_surf(:), start = jstartnc, count = jcountnc)              
    iret = nf90_put_var(fid, e2s_id, array_eta_2_surf(:), start = jstartnc, count = jcountnc)              
    iret = nf90_put_var(fid, e3s_id, array_eta_3_surf(:), start = jstartnc, count = jcountnc)              
    iret = nf90_put_var(fid, e4s_id, array_eta_4_surf(:), start = jstartnc, count = jcountnc)              

    iret = nf90_put_var(fid, e1a_id, array_eta_1_above(:), start = jstartnc, count = jcountnc)              
    iret = nf90_put_var(fid, e2a_id, array_eta_2_above(:), start = jstartnc, count = jcountnc)              
    iret = nf90_put_var(fid, e3a_id, array_eta_3_above(:), start = jstartnc, count = jcountnc)              
    iret = nf90_put_var(fid, e4a_id, array_eta_4_above(:), start = jstartnc, count = jcountnc)              
!!$
!!$
    iret = nf90_put_var(fid, lw_id, array_down_lw_ac(:), start = jstartnc, count = jcountnc)              
    iret = nf90_put_var(fid, sw_id, array_down_sw_ac(:), start = jstartnc, count = jcountnc) 

    iret = nf90_put_var(fid, hs_id, array_heat_soil(:), start = jstartnc, count = jcountnc) 
!!$
    iret = nf90_put_var(fid, tz_id, array_theta_zero(:), start = jstartnc, count = jcountnc) 
!!$
!!$    iret = nf90_put_var(fid, jo_id, jomega(:), start = jstartnc, count = jcountnc) 
!!$    iret = nf90_put_var(fid, josw_id, jomega_sw(:), start = jstartnc, count = jcountnc) 
!!$    iret = nf90_put_var(fid, jos_id, array_jomega_surf(:), start = jstartnc, count = jcountnc) 
!!$    iret = nf90_put_var(fid, je_id, jemissivity(:), start = jstartnc, count = jcountnc) 
!!$    iret = nf90_put_var(fid, jes_id, array_jemissivity_surf(:), start = jstartnc, count = jcountnc) 
!!$    iret = nf90_put_var(fid, josws_id, array_jomega_sw_surf(:), start = jstartnc, count = jcountnc) 
!!$
!!$    iret = nf90_put_var(fid, ngu_bign_id, array_ngu_bign(:), start = jstartnc, count = jcountnc)              
!!$    iret = nf90_put_var(fid, ngu_bign_s_id, array_ngu_bign_surf(:), start = jstartnc, count = jcountnc)              
!!$    iret = nf90_put_var(fid, ngu_bign_a_id, array_ngu_bign_above(:), start = jstartnc, count = jcountnc)     
!!$
!!$
!!$    iret = nf90_put_var(fid, ngu_alpha_id, array_ngu_alpha(:), start = jstartnc, count = jcountnc)              
!!$    iret = nf90_put_var(fid, ngu_alpha_s_id, array_ngu_alpha_surf(:), start = jstartnc, count = jcountnc)              
!!$    iret = nf90_put_var(fid, ngu_alpha_a_id, array_ngu_alpha_above(:), start = jstartnc, count = jcountnc)   
!!$
!!$
!!$    iret = nf90_put_var(fid, ngu_alphai_id, array_ngu_alphai(:), start = jstartnc, count = jcountnc)              
!!$    iret = nf90_put_var(fid, ngu_alphai_s_id, array_ngu_alphai_surf(:), start = jstartnc, count = jcountnc)              
!!$    iret = nf90_put_var(fid, ngu_alphai_a_id, array_ngu_alphai_above(:), start = jstartnc, count = jcountnc)   

    iret = nf90_put_var(fid, msa_id, array_light_surf_absorbed(:), start = jstartnc, count = jcountnc)              
    iret = nf90_put_var(fid, mla_id, array_light_layer_absorbed(:), start = jstartnc, count = jcountnc) 

    iret = nf90_put_var(fid, ibh_id, intermediate_box_height(:), start = jstartnc, count = jcountnc) 

    ! close the netCDF file and write to disk                   
    iret = nf90_close(fid)


  END SUBROUTINE mleb_netcdf



  !! ===========================================================================================================================
  !! SUBROUTINE: mleb_lwrad
  !!
  !!
  !>\BRIEF: calculates the radiation balance with a multilayer canopy, including the radiation matrix, 
  !!           as outlined in reference below  
  !!
  !! DESCRIPTION:
  !! 
  !!
  !! RECENT CHANGE(S): None
  !!
  !! MAIN OUTPUT VARIABLE(S): eta_1, eta_2, eta_3, eta_4, jomega, jemissivity, 
  !!                            temp_leaf_pres_grid, temp_surf_pres,
  !!                            pad_deltah, jomega_sw, eta_1_surf, 
  !!                            eta_2_surf, eta_3_surf, eta_4_surf,
  !!                            jemissivity_surf, jomega_surf, jomega_sw_surf, 
  !!                            eta_3_above, eta_1_above, eta_4_above, 
  !!                            eta_2_above, ngu_alpha, ngu_alpha_insert, 
  !!                            kjpindex, big_k_lw, big_k_sw
  !!
  !! REFERENCE(S): Gu, L. et al., 1999. Micrometeorology, biophysical exchanges and NEE 
  !!                  decomposition in a two-storey boreal forest - development and test 
  !!                  of an integrated model. Agricultural and Forest Meteorology, 94, pp.123–148.
  !!
  !! FLOWCHART:
  !_ ==============================================================================================================================

  SUBROUTINE mleb_lwrad (jomega, jemissivity, temp_leaf_pres_grid, temp_surf_pres, &
       & pad_deltah, jemissivity_surf, jomega_surf, jomega_sw_surf, &
       & eta_3_above, eta_1_above, ngu_alpha, ngu_alpha_insert, kjpindex, &
       & big_k_lw, big_k_sw) 

    IMPLICIT NONE


    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                             :: kjpindex               !! Domain size (-)

    REAL(r_std), DIMENSION (kjpindex,jnlvls), INTENT(in)   :: temp_leaf_pres_grid    ! leaf temperature down the column (K) PRESENT STEP
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)          :: temp_surf_pres
    REAL(r_std), DIMENSION (0:jnlvls), INTENT(in)          :: pad_deltah
    REAL(r_std), INTENT(in)                                :: big_k_lw
    REAL(r_std), INTENT(in)                                :: big_k_sw

    !! 0.2 Output variables
    REAL(r_std), INTENT(out)                               :: jemissivity_surf                     ! emissivity of surface layer
    REAL(r_std), INTENT(out)                               :: jomega_surf
    REAL(r_std), INTENT(out)                               :: jomega_sw_surf
    REAL(r_std), INTENT(out)                               :: eta_3_above, eta_1_above

    REAL, DIMENSION (0:jnlvls+1, 0:jnlvls+1), INTENT(out)  :: ngu_alpha
    REAL(r_std), DIMENSION(0:jnlvls+1), INTENT(out)        :: ngu_alpha_insert

    REAL(r_std), DIMENSION(jnlvls), INTENT(out)            :: jomega

    REAL(r_std), DIMENSION(jnlvls), INTENT(out)            :: jemissivity


    !! 0.3 Modified variables


    !! 0.4 Local variables
    INTEGER(i_std)                                         :: i, j, k, m, ji, jk

    REAL(r_std)                                            :: rho_albedo_surf
    REAL(r_std)                                            :: jcounter, jcounter2
    REAL(r_std)                                            :: series_sub1_surf
    REAL(r_std)                                            :: series_sub2_surf, series_sub1_sw_surf
    REAL(r_std)                                            :: temp_surf
    REAL(r_std)                                            :: jfactor

    REAL(r_std), DIMENSION(jnlvls)                         :: temp_leaf_pres    ! leaf temperature down the column (K) PRESENT STEP
    REAL(r_std), DIMENSION(jnlvls)                         :: series_sub1, series_sub2
    REAL(r_std), DIMENSION(jnlvls)                         :: series_sub1_sw, series_sub2_sw  
    REAL(r_std), DIMENSION(jnlvls)                         :: rho_albedo                           ! albedo of each vegetation layer
    REAL(r_std), DIMENSION (jnlvls)                        :: jomega_sw

    REAL(r_std), DIMENSION (0:jnlvls)                      :: jpad    
    REAL(r_std), DIMENSION (0:jnlvls)                      :: acc_pad        
    REAL(r_std)                                            :: total_ngu_bign 

    !! ==============================================================================================================================

    !! TO DO : Clean what is not used (Basically all the first part)

    DO jk = 1, kjpindex

       !! Initialisations:

       temp_leaf_pres(:) = temp_leaf_pres_grid(jk,:)  
       jemissivity(:) = 0.96d0 != 0.96d0
       jomega_sw_surf = 1.0d0
       jomega_surf = 1.0d0         
       jemissivity_surf = 0.96d0 != 0.96d0
       temp_surf =  0.0d0  !(not used)             
       series_sub1(:) = 1.0d0
       series_sub2(:) = 1.0d0 
       series_sub1_sw(:) = 1.0d0
       series_sub2_sw(:) = 1.0d0  
       rho_albedo(:) = 0.4d0
       rho_albedo_surf = 0.4d0           
       rho_albedo(:) = 0.4d0
       rho_albedo_surf = 0.4d0    
       jomega_surf = 1.0d0
       jomega_sw_surf = 1.0d0
       jemissivity_surf = 0.96d0
       jemissivity(:) = 0.96d0
       jcounter = 1.0d0
       jfactor = 1.0d0
       series_sub1_surf = 1.0d0

       ! calculate omega coefficient for each level
       DO i=1, jnlvls
          jomega(i) = 1 - (exp(-1.0d0 * big_k_lw * pad_deltah(i) ) )
       END DO ! i = 1, jnlvls

       DO i=1, jnlvls
          jomega_sw(i) = 1 - (exp(-1.0d0 * big_k_sw * pad_deltah(i) ) ) 
       END DO ! i = 1, jnlvls 



       eta_1_above = 0.0d0     

       eta_1_above = (series_sub1(1) * jomega(1) * (1.0d0 - jemissivity(1))) 

       DO i = 2, jnlvls-1
          eta_1_above = eta_1_above + &
               & (series_sub1(i) * jomega(i) * (1.0d0 - jemissivity(i))) 
       END DO ! i = 2, jnlvls-2

       eta_1_above = eta_1_above + &
            & series_sub1(jnlvls) * jomega(jnlvls) * (1.0d0 - jemissivity(jnlvls))

       eta_1_above = eta_1_above + &
            & series_sub1_surf * jomega_surf * (1.0d0 - jemissivity(jnlvls))


       eta_3_above = 0.0d0

       DO j = jnlvls, 1
          jfactor = 1.0d0
          DO k = jnlvls, j
             jfactor = jfactor * (1.0d0 - jomega(k))
          END DO ! k = jnlvls, j
          eta_3_above = eta_3_above + & 
               & (jomega(j) * jfactor * c_stefan * temp_leaf_pres(j)**4.0d0)
       END DO ! j = jnlvls, 1

       eta_3_above = eta_3_above + (jomega_surf * jfactor * c_stefan * temp_surf_pres(jk)**4.0d0)


       ! ------------------------------------------------------------------------------------

       !! END TO DO

       ! we here calculate the radiation matrix

       m = jnlvls

       DO i = 0, jnlvls
          jpad(i) = pad_deltah(i)
       END DO ! i = 0, jnlvls


       ! define acc_pad

       acc_pad(jnlvls) = 0.0d0

       DO i = m-1, 0, -1
          acc_pad(i) = jpad(i+1) + acc_pad(i+1)
       END DO ! i = m-1, 0, -1

       IF (printlev_loc>=4) THEN
          WRITE (numout,*) 'in jradiation matrix  acc_pad: ', acc_pad(:)
       END IF ! (printlev_loc>=4) THEN


       DO i = 0, m+1   
          DO j = 0, m+1

             IF ((i .eq. 0) .and. (j .eq. 0)) THEN
                ngu_alpha(i,j) = -1.0d0
                ngu_alpha_insert(i) = -1.0d0 

             ELSE IF ((i .eq. 0) .and. ((j .ge. 1) .and. (j .le. m))) THEN
                ngu_alpha(i,j) = gu(acc_pad(0) - acc_pad(j-1)) - gu(acc_pad(0)-acc_pad(j))

             ELSE IF ((i .eq. 0) .and. (j .eq. (m+1)))  THEN
                ngu_alpha(i,j) = gu(acc_pad(0))

             ELSE IF (((i .ge. 1) .and. (i .le. m)) &
                  &  .and. ((j .ge. 1) .and. (j .le. (i-1)))) THEN
                ngu_alpha(i,j) = gu(acc_pad(j) - acc_pad(i-1)) - gu(acc_pad(j-1)-acc_pad(i-1)) &
                     & - gu(acc_pad(j) - acc_pad(i)) + gu(acc_pad(j-1)-acc_pad(i))

             ELSE IF (((i .ge. 1) .and. (i .le. m)) .and. (j .eq. i))  THEN

                ! this is a departure from the matrix scheme, because I am treating 
                !   the temperature of the infinitesimal soil surface layer, and each 
                !   vegetation layer, as implicit when (i = j). Radiation from each of the 
                !   other levels is treated explicitly - i.e. by using the Gu scheme

                ! following line is from the full scheme, for completeness
                ! ngu_alpha(i,j) = 2.0d0 * gu(acc_pad(i-1) - acc_pad(i)) - 2.0d0 
                ngu_alpha_insert(i) = gu(acc_pad(i-1) - acc_pad(i)) 
                ! test to see influence of absorbed/emitted from the same layer in the cycle
                ngu_alpha(i,j) = 0.0d0 

             ELSE IF (((i .ge. 1).and.(i .le. m)).and.(((j .ge. (i+1)) .and. (j .le. m))))  THEN
                ngu_alpha(i,j) = gu(acc_pad(i) - acc_pad(j-1)) - gu(acc_pad(i)-acc_pad(j))  &
                     & - gu(acc_pad(i-1) - acc_pad(j-1)) + gu(acc_pad(i-1)-acc_pad(j))

             ELSE IF ((i .eq. (m+1)) .and. (j .eq. 0))  THEN
                ngu_alpha(i,j) = gu(acc_pad(0))

             ELSE IF ((i .eq. (m+1)) .and. ((j .ge. 1) .and. (j .le. m)))  THEN
                ngu_alpha(i,j) = gu(acc_pad(j)) - gu(acc_pad(j-1))   

             ELSE IF ((i .eq. (m+1)) .and. (j .eq. (m+1)))  THEN
                ngu_alpha(i,j) = -1.0d0
                ngu_alpha_insert(i) = -1.0d0   

             ELSE IF (((i .ge. 1) .and. (i .le. m)) .and. (j .eq. (0)))  THEN
                ngu_alpha(i,j) = gu(acc_pad(0) - acc_pad(i-1)) - gu(acc_pad(0) - acc_pad(i))

             ELSE IF (((i .ge. 1) .and. (i .le. m)) .and. (j .eq. (m+1)))  THEN
                ngu_alpha(i,j) = gu(acc_pad(i)) - gu(acc_pad(i-1))  

             ELSE 

             END IF

          END DO ! j = 0, m+1

       END DO ! i = 0, m+1

    END DO ! jk = 1, kjpindex

    ! end of radiation matrix calculate

  END SUBROUTINE mleb_lwrad



  !!  ===========================================================================================================================
  !! SUBROUTINE: mleb_calc_column
  !!
  !!
  !>\BRIEF: calculates the energy, for either the multilayer or the single layer case 
  !!
  !! DESCRIPTION: This module calculates the coefficients that are required to determine the energy 
  !!              budget at each time step. See also sec. 3.5 in Ryder et al., 2016.
  !! 
  !!
  !! RECENT CHANGE(S): None
  !!
  !! MAIN OUTPUT VARIABLE(S): big_e, big_f, big_g, flux_h_grid, flux_le_grid, 
  !!                               a_t, b_t, c_t, d_t, a_q, b_q, c_q, d_q, 
  !!                               t_a_next_grid, q_a_next_grid, horiz_flux_t, horiz_flux_q, 
  !!                               temp_surf_next, flux_ground_h, flux_ground_le, delta_temp, 
  !!                               temp_leaf_next_grid, temp_sol_new, epot_air_new, 
  !!                               big_r, big_r_prime, speed, vbeta3, eta_1, eta_2, eta_3, eta_4, 
  !!                               eta_1_surf, eta_2_surf, eta_3_surf, eta_4_surf, qair_new, 
  !!                               qair, fluxsens, fluxlat, vbeta2sum, vbeta3sum
  !!
  !! REFERENCE(S): Best, M.J. et al., 2004. A proposed structure for coupling tiled surfaces 
  !!                 with the planetary boundary layer. Journal of Hydrometeorology, 
  !!                 5, pp.1271–1278.
  !!
  !!               Polcher, J. et al., 1998. A proposal for a general interface between land 
  !!                  surface schemes and general circulation models. Global and Planetary 
  !!                  Change, 19, pp.261–276
  !!    
  !!               Ryder et al., 2016, GMD, doi:10.5194/gmd-9-223-2016 
  !!
  !! FLOWCHART: 
  !_ ==============================================================================================================================
  !+++CHECK+++
  ! Variable dimensions xxx_Tot_mean are for a single pixel. Causing 1+1 problems
  ! when running a larger domain. Needs to be corrected when implementing
  ! a global use of the multi-layer energy budget.
  SUBROUTINE mleb_calc_column (pad_deltah, layer_height, big_r, big_r_prime, delta_z,      &
       & k_eddy, swnet, swdown, lwdown,                                    &
       & temp_atmos_pres_grid, q_atmos_pres_grid,                          &
       & a_t, b_t, c_t, d_t,                                               &
       & a_q, b_q, c_q, d_q,                                               &
       & big_e, big_f, big_g,                                              &
       & jalpha, jbeta, temp_leaf_pres_grid, delta_h,                      &
       & eta_1, eta_2, eta_3, eta_4,                                       &
       & jtheta,                                                           &
       & jalpha_surf, jbeta_surf, rau,                                     &
       & soilcap, soilflx, speed,                                          &
       & q_cdrag,                                                          &
       & kjpindex,                                                         &
       & flux_ground_h, flux_ground_le, flux_h_grid, flux_le_grid,         &
       & u_speed, t_a_next_grid, q_a_next_grid,                            &
       & temp_surf_next,                                                   &
       & eta_1_surf, eta_2_surf, eta_3_surf, eta_4_surf,                   &
       & temp_leaf_next_grid, epot_air_new, delta_temp,                    &
       & horiz_flux_t, horiz_flux_q,                                       &
       & ngu_bign, ngu_alpha, ngu_alpha_insert, qair_new, qair,            &
       & emis, temp_sol_new, vbeta1, vbeta4, vbeta5, vbeta2,               &
       & vbeta3, vbeta3pot, vbeta, fluxsens, fluxlat, netrad,              &
       & petAcoef, petBcoef, peqAcoef, peqBcoef,                           &
       & u, v, E_supply,                                                   &
       & hydrol_flag, hydrol_flag2,                                        &
       & hydrol_flag3, vbeta2sum, vbeta3sum, veget_max, veget,             &
       & Light_Abs_Tot_mean, Light_Alb_Tot_mean, &
       & light_surf_absorbed,      &
       & light_layer_absorbed, transpir_supply_column, pb, profile_vbeta3, &
       & profile_rveget,                                                   &
       & flux_netrad, flux_netsw, flux_netlw,                              &
       & ks_tune, ks_slope, ks_veget,qsol_sat_new, transpir) 
    !+++++++++++


    ! 

    IMPLICIT NONE


    ! 0.1 input variables


    INTEGER(i_std), INTENT(in)                                    :: kjpindex                   !! Domain size (-)
    REAL(r_std), INTENT (in)                                      :: layer_height               !! layer height (m)

    REAL(r_std), DIMENSION (kjpindex, jnlvls), INTENT(in)         :: temp_leaf_pres_grid        !! leaf temperature down the column (K) PRESENT STEP
    REAL(r_std), DIMENSION (kjpindex, jnlvls), INTENT(in)         :: temp_atmos_pres_grid       !! atmospheric temperature at the present timestep (K)
    REAL(r_std), DIMENSION (kjpindex, jnlvls), INTENT(in)         :: q_atmos_pres_grid          !! specific humidity at the present timestep (kg/kg)


    REAL(r_std), DIMENSION (kjpindex), INTENT (in)                :: q_cdrag                    !! Surface drag coefficient (-)
    REAL(r_std), DIMENSION (jnlvls), INTENT(in)                   :: jtheta                     !! specific heat capacity of air (J / (kg K) )

    REAL(r_std), DIMENSION (jnlvls), INTENT (inout)               :: k_eddy                     !! eddy diffusivity between each level (m^2 / s)
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)                :: lwdown                     !! Down-welling long-wave flux (W m^{-2})
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)                :: swnet                      !! Net surface short-wave flux (W m^{-2})
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)                :: swdown                     !! Down-welling short-wave flux (W m^2)
    REAL(r_std), DIMENSION (0:jnlvls), INTENT(in)                 :: pad_deltah                 !! leaf area density at each level (m^2 / m^3)
    REAL(r_std), DIMENSION (jnlvls), INTENT(in)                   :: jalpha, jbeta
    REAL(r_std), DIMENSION(kjpindex), INTENT(in)                  :: jalpha_surf, jbeta_surf
    REAL(r_std), DIMENSION (jnlvls), INTENT(in)                   :: delta_z                    !! separation of levels (m)
    REAL(r_std), DIMENSION(jnlvls), INTENT(in)                    :: delta_h                    !! height of each level (m)


    REAL(r_std),DIMENSION (kjpindex), INTENT (in)                 :: soilcap                    !! Soil calorific capacity (J K^{-1])
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)                :: soilflx                    !! Soil flux (W m^{-2})

    REAL(r_std), DIMENSION (kjpindex,jnlvls), INTENT(in)                   :: u_speed                    !! wind speed at each level (m/s)

    REAL(r_std), DIMENSION (kjpindex), INTENT (in)                :: emis                       !! Emissivity (-)
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)                 :: rau                        !! Air density (kg m^{-3})
    REAL(r_std), DIMENSION (0:jnlvls+1, 0:jnlvls+1), INTENT(in)   :: ngu_alpha


    REAL(r_std), DIMENSION (kjpindex), INTENT (in)                :: u                          !! Eastward Lowest level wind speed  (m s^{-1})
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)                :: v                          !! Northward Lowest level wind speed (m s^{-1})
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)                :: petAcoef                   !! PetAcoef (see note)
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)                :: petBcoef                   !! PetBcoef (see note)
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)                :: peqAcoef                   !! PeqAcoef (see note)
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)                :: peqBcoef                   !! PeqBcoef (see note) 
    REAL(r_std), DIMENSION (kjpindex), INTENT (inout)             :: vbeta                      !! Resistance coefficient (-)
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)                :: vbeta1                     !! Snow resistance (-) 
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)                :: vbeta4                     !! Bare soil resistance (-)
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)                :: vbeta5                     !! Floodplains resistance
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT (in)            :: vbeta2                     !! Interception resistance (-)
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT (in)            :: vbeta3pot                  !! Vegetation resistance for potential transpiration
    !    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT (in)            :: vbetaco2                   !! Vegetation resistance to CO2 (-)
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT (in)            :: E_supply                   !! Supply of water for transpiration @tex $(m.s^{-1})$ @endtex
    REAL(r_std),  DIMENSION (kjpindex, nvm), INTENT(in)           :: veget_max                  !! Max. fraction of vegetation type (LAI -> infty, unitless)
    REAL(r_std),  DIMENSION (kjpindex, nvm), INTENT(in)           :: veget                      !! Vegetation fraction for each pft 
    REAL(r_std), DIMENSION(nlevels_tot,kjpindex,nvm), INTENT(in)  :: transpir_supply_column     !! Supply of water for transpiration 

    REAL(r_std),DIMENSION (kjpindex,nvm,nlevels_tot), INTENT(in)  :: profile_vbeta3
    REAL(r_std),DIMENSION (kjpindex,nvm,nlevels_tot), INTENT(in)  :: profile_rveget
    REAL(r_std), INTENT(in)                                       :: ks_tune                    !! coefficient for tuning the forest floor evapotranspiration
    REAL(r_std), INTENT(in)                                       :: ks_slope                   !! coefficient for S-shape function 
    REAL(r_std), INTENT(in)                                       :: ks_veget                   !! coefficient for S-shape function

    ! 0.2 output variables
    REAL(r_std), DIMENSION (jnlvls), INTENT(out)                  :: big_e, big_f, big_g
    REAL(r_std), DIMENSION (kjpindex, jnlvls), INTENT(out)        :: flux_h_grid, flux_le_grid  !! sensible and latent heat flux respectively at each level (both W/m^2)
    REAL(r_std), DIMENSION (jnlvls), INTENT(out)                  :: a_t, b_t, c_t, d_t
    REAL(r_std), DIMENSION (jnlvls), INTENT(out)                  :: a_q, b_q, c_q, d_q 
    REAL(r_std), DIMENSION (kjpindex, jnlvls), INTENT(out)        :: t_a_next_grid
    REAL(r_std), DIMENSION (kjpindex, jnlvls), INTENT(out)        :: q_a_next_grid
    REAL(r_std), DIMENSION (jnlvls), INTENT(out)                  :: horiz_flux_t,  horiz_flux_q
    REAL(r_std), INTENT(out)                                      :: temp_surf_next             !! temperature of surface at respectively present and next timestep (K)
    REAL(r_std), INTENT(out)                                      :: flux_ground_h, flux_ground_le    !! sensible and latent heat flux respectively from the ground (both W/m^2)
    REAL(r_std), INTENT(out)                                      :: delta_temp                 !! change in temperature over the time step (Kelvins)
    REAL(r_std), DIMENSION (kjpindex, jnlvls), INTENT(out)        :: temp_leaf_next_grid
    REAL(r_std), INTENT(out)                                      :: light_surf_absorbed
    REAL(r_std), DIMENSION(jnlvls), INTENT(out)                   :: light_layer_absorbed

    REAL(r_std), DIMENSION (jnlvls), INTENT(out)                  :: flux_netrad                !! net radiation flux for each layer  
    REAL(r_std), DIMENSION (jnlvls), INTENT(out)                  :: flux_netsw                 !! net radiation contributed from short-wave
    REAL(r_std), DIMENSION (jnlvls), INTENT(out)                  :: flux_netlw                 !! net radiation contributed from long-wave

    ! 0.3 modified variables
    REAL(r_std), DIMENSION (kjpindex), INTENT(inout)              :: temp_sol_new               !! New surface temperature (K)
    REAL(r_std), DIMENSION (kjpindex), INTENT(inout)              :: epot_air_new
    REAL(r_std), DIMENSION (jnlvls), INTENT(inout)                :: big_r, big_r_prime         !! stomatal resistance for sensible and latent heat flux at each level (s/m)    
    REAL(r_std), INTENT(inout)                                    :: speed                      !! wind speed above the canopy (m/s), as defined by old enerbil
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT (inout)         :: vbeta3                     !! Vegetation resistance (-)
    REAL(r_std), DIMENSION (jnlvls), INTENT (inout)               :: eta_1, eta_2
    REAL(r_std), DIMENSION (jnlvls), INTENT (inout)               :: eta_3, eta_4
    REAL(r_std), INTENT (inout)                                   :: eta_1_surf, eta_2_surf
    REAL(r_std), INTENT (inout)                                   :: eta_3_surf, eta_4_surf
    REAL(r_std), DIMENSION (kjpindex), INTENT (inout)             :: qair_new
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)                :: qair
    REAL(r_std), DIMENSION (kjpindex), INTENT(inout)              :: fluxsens                   !! Sensible heat flux (W m^{-2})
    REAL(r_std), DIMENSION (kjpindex), INTENT(inout)              :: fluxlat                    !! Latent heat flux (W m^{-2})     
    REAL(r_std), DIMENSION (kjpindex), INTENT (inout)             :: vbeta2sum, vbeta3sum
    REAL(r_std), DIMENSION (kjpindex), INTENT(inout)              :: netrad                     !! Net radiation flux (W m^{-2})   added by YC 

    !! Debug ID3.
    !! Transpiration considered as the sum of all the layers' transpiration
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT (inout)         :: transpir                     !! Transpiration
    !! End ID3
    

    !+++CHECK+++
    ! Variable dimensions xxx_Tot_mean are for a single pixel. Causing 1+1 problems
    ! when running a larger domain. Needs to be corrected when implementing
    ! a global use of the multi-layer energy budget.
    REAL(r_std),DIMENSION(nlevels_tot), INTENT(in)                :: Light_Abs_Tot_mean         !! total light absorption for a given canopy level
    REAL(r_std),DIMENSION(nlevels_tot), INTENT(in)                :: Light_Alb_Tot_mean         !! total albedo for a given level
    !+++++++++++

    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)   :: qsol_sat_new         !! Air saturated humidity (kg.kg-1)

    ! 0.4 local variables
    INTEGER(i_std)                                                :: i, j, m, n, ji, jv, jk, jm
    REAL(r_std), DIMENSION (jnlvls)                               :: x_1, x_2, x_3, x_4, x_5, x_6
    REAL(r_std), DIMENSION (jnlvls)                               :: y_1, y_2, y_3, y_4, y_5, y_6 
    REAL(r_std), DIMENSION (jnlvls)                               :: t_a_pres
    REAL(r_std), DIMENSION (jnlvls)                               :: q_a_pres
    REAL(r_std), DIMENSION (jnlvls)                               :: t_l_pres, t_l_next
    REAL(r_std), DIMENSION (jnlvls)                               :: q_l_pres, q_l_next
    REAL(r_std), DIMENSION (jnlvls)                               :: big_e_top
    REAL(r_std), DIMENSION (jnlvls)                               :: big_f_top
    REAL(r_std), DIMENSION (jnlvls)                               :: big_g_top
    REAL(r_std), DIMENSION (jnlvls)                               :: big_efg_bottom
    REAL(r_std), DIMENSION (jnlvls)                               :: source_sink_t, source_sink_q
    REAL(r_std), DIMENSION (jnlvls)                               :: source_sink_sum_t, source_sink_sum_q
    REAL(r_std), DIMENSION (jnlvls)                               :: source_sink_t_l
    REAL(r_std), DIMENSION (jnlvls)                               :: temp_leaf_pres           !! leaf temperature down the column (K) PRESENT STEP
    REAL(r_std), DIMENSION (jnlvls)                               :: temp_atmos_pres          !! atmospheric temperature at the present timestep (K)
    REAL(r_std), DIMENSION (jnlvls)                               :: q_atmos_pres             !! specific humidity at the present timestep (kg/kg)
    REAL(r_std), DIMENSION (jnlvls)                               :: flux_h                   !! sensible heat flux at each level (both W/m^2)
    REAL(r_std), DIMENSION (jnlvls)                               :: flux_le                  !! latent heat flux at each level (both W/m^2)
    REAL(r_std), DIMENSION (jnlvls)                               :: t_a_next
    REAL(r_std), DIMENSION (jnlvls)                               :: q_a_next

    REAL(r_std)                                                   :: big_omega_1, big_omega_2
    REAL(r_std)                                                   :: big_omega_3, big_omega_4
    REAL(r_std)                                                   :: big_omega_5, big_omega_6
    REAL(r_std)                                                   :: big_omega_7, big_omega_8
    REAL(r_std)                                                   :: big_omega_5_sub, big_omega_6_sub
    REAL(r_std)                                                   :: big_omega_7_sub, big_omega_8_sub
    REAL(r_std)                                                   :: xi_1, xi_2, xi_3, xi_4

    !    REAL(r_std)                                                   :: source_sink_sum_t
    !    REAL(r_std)                                                   :: source_sink_sum_q

    REAL(r_std)                                                   :: fevap, turb_factor
    REAL(r_std)                                                   :: M_1_h, N_1_h, M_1_q, N_1_q

    REAL(r_std)                                                   :: jtchalev0, jtcp_air, jtrau
    REAL(r_std), DIMENSION (kjpindex)                             :: beta_diffuco                !! the 'beta coefficient' from the existing scheme, for the surface layer (-) 
    REAL(r_std), DIMENSION (kjpindex)                             :: beta_diffuco_evap           !! the 'beta coefficient' from the existing scheme, for the surface layer (-)
    REAL(r_std), DIMENSION (kjpindex)                             :: beta_diffuco_sub            !! the 'beta coefficient' from the existing scheme, for the surface layer (-)
    REAL(r_std), DIMENSION (kjpindex)                             :: beta_diffuco_evap_notrans
    REAL(r_std), DIMENSION (kjpindex)                             :: beta_diffuco_evap_trans
    REAL(r_std), DIMENSION (kjpindex,nvm)                             :: delta_temp_canopy

    REAL(r_std), DIMENSION (kjpindex)                             :: evap_fraction               !! a evapotranspiration fraction based on vbeta3, vbeta4 and other paramters  

    REAL(r_std), DIMENSION (jnlvls)                               :: delta_t_a, sum_flux_a
    REAL(r_std), DIMENSION (jnlvls)                               :: delta_q_a, sum_flux_q_a

    REAL(r_std), DIMENSION (jnlvls)                               :: delta_temp_leaf
    REAL(r_std), DIMENSION (jnlvls)                               :: sum_flux_leaf

    REAL(r_std), DIMENSION (0:jnlvls+1)                           :: gu_alpha_sum
    REAL(r_std), DIMENSION (0:jnlvls+1)                           :: ngu_temp_solid
    REAL(r_std), DIMENSION (0:jnlvls+1), INTENT(out)              :: ngu_bign
    REAL(r_std), DIMENSION (0:jnlvls+1), INTENT(in)               :: ngu_alpha_insert    
    REAL(r_std), DIMENSION (0:jnlvls+1)                           :: ngu_alpha_insert2    

    LOGICAL                                                       :: hydrol_flag                  !! Flag that recalculates the energy budget for the new supply term (true/false)    
    LOGICAL, DIMENSION (kjpindex), INTENT(in)                     :: hydrol_flag2
    LOGICAL, DIMENSION (kjpindex, nvm), INTENT(in)                :: hydrol_flag3                 !! flag that 'trips' the alternative energy budget for each grid square

    REAL(r_std), DIMENSION (0:jnlvls)                             :: pad_deltah_temp              !! leaf area density at each level (m^2 / m^3)

    INTEGER                                                       :: convergence_count

    REAL(r_std), DIMENSION (jnlvls)                               :: phi_trans


    ! 090215 new column variables

    REAL(r_std), DIMENSION (jnlvls)                               :: column_air_relhum
    REAL(r_std), DIMENSION (jnlvls)                               :: column_q_surf_leaf
    REAL(r_std), DIMENSION (jnlvls)                               :: column_r_sv_leaf
    REAL(r_std), DIMENSION (jnlvls)                               :: column_r_sto_v_leaf
    REAL(r_std), DIMENSION (jnlvls)                               :: column_evap
    REAL(r_std), DIMENSION (jnlvls)                               :: column_qsatt_leaf

    REAL(r_std), DIMENSION (kjpindex), INTENT (in)                :: pb                           !! Lowest level pressure (Pa)

    REAL(r_std), DIMENSION (kjpindex)                             :: ks_weight_1                  !! for the calculation of evap_fraction 
    REAL(r_std), DIMENSION (kjpindex)                             :: ks_weight_2                  !! weight_1 for vbeta3; weight_2 for vbeta4 
    REAL(r_std)                                                   :: temp_sol_err                 !! the absolute error of temp_sol surface temperatrure 
    REAL(r_std)                                                   :: epsilon_err                  !! mechain error + numerical error set to 1.0E-15  
    REAL(r_std)                                                   :: temp_sol_dummy1
    REAL(r_std)                                                   :: temp_sol_dummy2
    REAL(r_std)                                                   :: vbeta3_pixel                 !! temporary variable to calculate vbeta3 at the pixel level
    !! ==============================================================================================================================
    !      
    epsilon_err = 1.0E-15

    !! Initialisation of several variables:

    phi_trans(:) = 0.0d0 !1.0d4

    flux_h(:) = 0.0d0; flux_le(:) = 0.0d0

    source_sink_t(:) = 0.0d0; source_sink_q(:) = 0.0d0

    temp_surf_next = 0.0d0

    big_e = 0.0d0; big_f = 0.0d0; big_g = 0.0d0

    i = 0; ji = 0

    x_1 = 0.0d0; x_2 = 0.0d0; x_3 = 0.0d0; 
    x_4 = 0.0d0; x_5 = 0.0d0; x_6 = 0.0d0
    y_1 = 0.0d0; y_2 = 0.0d0; y_3 = 0.0d0; 
    y_4 = 0.0d0; y_5 = 0.0d0; y_6 = 0.0d0

    ! t_a_pres = temp_atmos_pres
    t_a_next = 0.0d0
    ! q_a_pres = q_atmos_pres
    q_a_next = 0.0d0

    t_l_next = 0.0d0; q_l_pres = 0.0d0; q_l_next = 0.0d0

    !! Set resistances of uppermost layer to <approaching infinity>
    big_r(jnlvls) = 1.0d20; big_r_prime(jnlvls) = 1.0d20

    light_surf_absorbed = 0.0d0
    light_layer_absorbed = 0.0d0


    !! Debug

    ! 121013 we need to re-calculate beta_diffuco_evap_trans, specifically vbeta3sum

    IF (printlev_loc>=4) THEN
       WRITE(numout, *) 'in mleb_calc_column, big_r_prime is: ', big_r_prime
       WRITE(numout, *) 'in mleb_calc_column, profile_rveget is: ', profile_rveget(test_grid,test_pft, :)
    END IF ! (printlev_loc>=4)


    ! hydrol_flag = .FALSE. ! AL: this is added at somepoint after DOFOCO
    ! rev2566. I have asked JR whether this is needed or not.


    gridloop: DO ji=1, kjpindex  ! -----------------------------------------------------------------

       !! TO DO: Think about what to do with this convergence loop:
       !!            - Do we keep it with in order to keep compatibility with previous hydraulic architecture?
       !!            - Do we remove everything?

       ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  start convergence loop here  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

       temp_sol_err = 0.0_r_std

       convergence_loop: DO jm = 1, 20

          !++++++++++copy the solved temperature for next time step i
          ! from previous covergence loop step result   
          temp_sol_dummy1 = temp_surf_next 
          !++++++++++++++++++++++++++++++


          ! calculate the revised vbeta3sum, based on the revised vbeta3

          IF (hydrol_flag) THEN  ! overall hydrol flag


             vbeta3sum(ji) = 0.0d0
             DO jv = 1, nvm
!!$                   !++++++++++
!!$                   IF(ji == test_grid .AND. jv == test_pft)THEN
!!$                      WRITE(numout,*) 'vbeta3 vnkldsvd 1',vbeta3(ji,jv)
!!$                   ENDIF
!!$                   !+++++++++++
                vbeta3sum(ji) = vbeta3sum(ji) + vbeta3(ji,jv)
             ENDDO ! jv = 1, nvm    

             vbeta(ji) = vbeta4(ji) + vbeta2sum(ji) + vbeta3sum(ji)

          END IF ! IF (hydrol_flag .eq. .TRUE.) 

          temp_leaf_pres(:) = temp_leaf_pres_grid(ji,:)
          t_a_pres(:) = temp_atmos_pres_grid(ji,:) 
          q_a_pres(:) = q_atmos_pres_grid(ji,:)
          t_l_pres(:) = temp_leaf_pres_grid(ji,:)


          !+++++ CHECK ++++++
          ! This loop causes a bug where we have incorrect results for large scale simulations.
          ! Apparently, if we exit here some variables are not getting properly calculated for the pixels
          ! which have no water stress.  This will cause a performance problem, so best to continue
          ! looking into this.  One variable where it was very noticable was temp_sol_new. For now, we 
          ! skip it using HACK_ENERBIL_HYDROL=y.
          IF(.NOT. HACK_ENERBIL_HYDROL) THEN
             IF (hydrol_flag) THEN
                ! are we on the second call of enerbil, with the corrections to transpiration
                IF (.NOT. hydrol_flag2(ji)) THEN
                   ! no need to run mleb_calc_column again if (E_supply .gt. transpir_mod) for all PFTs at this
                   !   particular grid point
                   ! WRITE (*,*) '091013 hydrol_flag2 is false so we are exiting the gridloop'
                   t_a_next_grid(ji,:) = temp_atmos_pres_grid(ji,:)
                   q_a_next_grid(ji,:) = q_atmos_pres_grid(ji,:)
                   temp_leaf_next_grid(ji,:) = temp_leaf_pres_grid(ji,:)

                   temp_sol_new(ji) = temp_sol_new(ji)
                   IF(test_grid == ji) WRITE(numout,*) "ruieozenerbil mleb_calc_column 1",temp_sol_new(ji)
                   temp_surf_next = temp_surf_next
                   qsol_sat_new(ji) = qsol_sat_new(ji)
                   qair_new(ji) = qair(ji)

                   fluxsens(ji) = fluxsens(ji)
                   fluxlat(ji) = fluxlat(ji)
                   vbeta2sum(ji) = vbeta2sum(ji)
                   vbeta3sum(ji) = vbeta3sum(ji)
                   vbeta3(ji,:) = vbeta3(ji,:)

                   CYCLE gridloop

                ELSE

                   ! hydrol_flag is not set, so carry on with the energy budget calculation

                END IF ! IF (hydrol_flag2(ji) .eq. .FALSE.)
             ELSE
                ! we are on the first call of enerbil
             END IF ! IF (hydrol_flag .eq. .TRUE.) 
          ENDIF
          !+++++++++++++++++++++++++++++++++++++++++++

          ! we are inside the loop over grid squares - speed is calculated from the horizontal
          !   velocities u and v, which are the input parameters
          speed = MAX(min_wind, SQRT (u(ji)*u(ji) + v(ji)*v(ji)))


          ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

          !! A note on the A and B coefficients:
          !! The values passed down from the coupled model are petAcoef, petBcoef, peqAcoef and
          !!   peqBcoef, allocated per gridsquare.
          !!
          !! NB: For the uncoupled mode, only b_t and b_q (at the uppermost level) have a value. 
          !!     All of the other coefficients are zero.
          !!
          !! However, for the coupled mode, a_t, b_t, a_q and b_q should all have values.
          !! c_t, d_t, c_q and d_q should not have values as they are derived internal to enerbil
          !! (they represent the interaction between the atmosphere in the canopy and the vegetation)
          !! and so are all initialised here.
          !!
          !! The interaction between the coefficients within the module 'surf_land_orchidee' (within LMDZ) 
          !! is as follows:
          !! 
          !! petAcoef (orchidee) = petBcoef (lmdz) * delta_t
          !!
          !! petBcoef (orchidee) = petAcoef (lmdz)
          !!
          !! peqAcoef (orchidee) = peqBcoef (lmdz) * delta_t
          !!
          !! peqBcoef (orchidee) = peqAcoef (lmdz)

          !! enerbil_equiv ========================================

          
          !! Debug:
          
          IF(printlev_loc >= 4) THEN
             WRITE (numout,*) 'start mleb_calc_column'
             WRITE (numout,*) 'a_t',a_t
             WRITE (numout,*) 'b_t',b_t
          ENDIF

          !! Definition of the coefficients for temperature:

          a_t(jnlvls) = petAcoef(ji)        ! NB: kjpindex = 1 corresponds to the bare soil PFT
          b_t(jnlvls) = petBcoef(ji)/cp_air ! NB: kjpindex = 1 corresponds to the bare soil PFT 
          
          !! NB: petBcoef is in units of surface static energy - (cp_air * temp) as pressure 
          !!     difference between the surface and the lowest layer of the atmosphere is
          !!     assumed negligable


          !! Debug:

          IF(printlev_loc >= 4) THEN
             WRITE (numout,*) 'a_t(jnlvls)',a_t(jnlvls)
             WRITE (numout,*) 'b_t(jnlvls)',b_t(jnlvls)
          ENDIF


          
          c_t(jnlvls) = 0.0d0
          d_t(jnlvls) = 0.0d0



          !! Definition of the coefficients for humidity:

          a_q(jnlvls) = peqAcoef(ji) ! NB: kjpindex = 1 corresponds to the bare soil PFT
          b_q(jnlvls) = peqBcoef(ji) ! NB: kjpindex = 1 corresponds to the bare soil PFT 
          !! NB: peqBcoef is already in units of (kg/kg)
          c_q(jnlvls) = 0.0d0
          d_q(jnlvls) = 0.0d0


          !! Definition of the beta coefficients (simplifies the code later):

          beta_diffuco(ji) =  1.0d0 * vbeta1(ji) * (un - vbeta5(ji)) + &
               & (un - vbeta1(ji)) * (un - vbeta5(ji)) * vbeta(ji) + vbeta5(ji)


          !! Seperation the beta_diffuco term above into expressions for sublimination (snow) and
          !! evaporation

          beta_diffuco_sub(ji) = vbeta1(ji) * (1.0d0 - vbeta5(ji))
          beta_diffuco_evap(ji) = (1.0d0 - vbeta1(ji)) * &
               & (1.0d0 -vbeta5(ji)) * vbeta(ji) + vbeta5(ji)


          !! Same definition of the version of the evaporation expression above for which the transpiration
          !! term is removed

          beta_diffuco_evap_notrans(ji) = (1.0d0 - vbeta1(ji)) * &
               & (1.0d0 - vbeta5(ji)) * (vbeta4(ji) + vbeta2sum(ji)) + vbeta5(ji)

          !! Same definition of the version of the evaporation expression above for which the transpiration
          !! term is included

          beta_diffuco_evap_trans(ji) = (1.0d0 - vbeta1(ji)) * &
               & (1.0d0 - vbeta5(ji)) * vbeta3sum(ji)





          !! Use of the results from mleb_lwrad to determine the eta coefficients in each layer.
          !! eta_4 is calculated thanks to variables coming from the module albedo_surface.
          !! mleb_lwrad solves Gu's Long Wave Radiation scheme (Gu et al. (1999)).

          IF (jnlvls .eq. 1) THEN

             eta_1_surf = 1.0d0 - (1.0d0 - emis(ji))
             eta_4_surf = 1.0d0

             eta_2_surf = - c_stefan * 4.0d0 * emis(ji) * (psold(ji)/cp_air)**3
             eta_3_surf = - 3.0d0 * c_stefan * emis(ji) * ((psold(ji)/cp_air)**4)

             eta_1(1) = 0.0d0 
             eta_2(1) = c_stefan * 4.0d0 * temp_leaf_pres(1)**3
             eta_3(1) = c_stefan * 3.0d0 * temp_leaf_pres(1)**4 
             eta_4(1) = 0.0d0

          ELSE

             !! Last calculations related to Gu's radiation scheme:

             eta_4 = 0.0d0
             eta_4_surf = 0.1d0

             ngu_temp_solid(0) = c_stefan * ((psold(ji)/cp_air)**4)

             DO i = 1, jnlvls
                ngu_temp_solid(i) = c_stefan * (temp_leaf_pres(i)**4) 
             END DO ! i = 1, jnlvls

             ngu_temp_solid(jnlvls+1) = lwdown(ji)

             ngu_bign(:) = 0.0d0

             ngu_alpha_insert2(0) = - 1.0d0 * ngu_temp_solid(0)
             DO i = 1, jnlvls
                ngu_alpha_insert2(i) = (2.0d0 * ngu_alpha_insert(i) - 2.0d0) * ngu_temp_solid(i)
             END DO ! i = 1, jnlvls 
             ngu_alpha_insert2(jnlvls+1) = - 1.0d0 * ngu_temp_solid(jnlvls+1)

             DO i = 0, jnlvls + 1 
                DO j = 0, jnlvls + 1
                   IF (i .ne. j) THEN
                      ngu_bign(i) = ngu_bign(i) + ( ngu_alpha(i,j) * ngu_temp_solid(j) )
                   ELSE
                      ngu_bign(i) = ngu_bign(i) + 0.0d0
                   END IF ! (i .ne. j)            
                END DO ! j = 0, jnlvls + 1
             END DO ! i = 0, jnlvls + 1



             !! Definition of all the eta coefficients for the LW radiation scheme:

             DO i = 1, jnlvls

                eta_1(i) = 0.0d0 

                eta_2(i) = (2.0d0*ngu_alpha_insert(i)-2.0d0) &
                     &  * c_stefan * 4.0d0 * temp_leaf_pres(i)**3

                eta_3(i) = -(2.0d0*ngu_alpha_insert(i) - 2.0d0) &
                     &  * c_stefan * 3.0d0 * temp_leaf_pres(i)**4 + ngu_bign(i)

             END DO ! i = 1, jnlvls

             !! The Gu method assumes black body emitters at the surface for simplicity:

             eta_1_surf = 0.0d0
             
             !! Debug : eta scheme to be the same as in the article

!!$             eta_2_surf = - 4.0d0 * c_stefan * (((1.0d0/cp_air) * psold(ji))**3)
             eta_2_surf = - 4.0d0 * c_stefan * (((1.0d0/cp_air) * psold(ji))**3)      

!!$             eta_3_surf = - (c_stefan * (3.0d0 * (psold(ji)/cp_air)**4) + ngu_bign(0))
             eta_3_surf = (ngu_bign(0) - (-1.0 *c_stefan * (3.0d0 * (psold(ji)/cp_air)**4)))

             !! End Debug

             !! Definition of the eta_4 coefficients thanks to the results of albedo_surface:

             !+++CHECK+++
             ! Variable dimensions xxx_Tot_mean are for a single pixel. Causing 1+1 problems
             ! when running a larger domain. Needs to be corrected when implementing
             ! a global use of the multi-layer energy budget.
             IF ((.NOT. hack_can_structure).AND. (swdown(ji) .GT. min_sechiba)) THEN !(SUM(Light_Abs_Tot_mean(:)) > min_sechiba)) THEN

                !! At the surface layer:

                light_surf_absorbed = 0.0d0

                !! Debug ID4.
                !! The light absorbed at the bottom layer is the light absorbed by the first layer in the albedo surface scheme
                !! Discussed with Matthew.

!!$                light_surf_absorbed = Light_Abs_Tot_mean(1)

                DO i=1,jnlvls_canopy+1
                   light_surf_absorbed = light_surf_absorbed + Light_Abs_Tot_mean(i)! + Light_Alb_Tot_mean(i)
                ENDDO
                light_surf_absorbed = MAX(0.,1 - (light_surf_absorbed + Light_Alb_Tot_mean(nlevels_tot)))

                !! End ID4
                
                !! Debug

                IF (printlev_loc>=4) THEN
                   WRITE (numout, *) 'light_surf_absorbed is: ', light_surf_absorbed 
                   WRITE (numout, *) 'Light_Alb_Tot_mean is: ', Light_Alb_Tot_mean
                   WRITE (numout, *) 'Light_Abs_Tot_mean is: ', Light_Abs_Tot_mean
                END IF

                !! Along the canopy layers:

                light_layer_absorbed = 0.0d0

                DO i = 1, jnlvls_canopy+1
                   light_layer_absorbed(i + jnlvls_under) = Light_Abs_Tot_mean(i)
                END DO
                
                !! Debug:

                IF (printlev_loc>=4) THEN
                   WRITE (numout, *) 'light_layer_absorbed is: ', light_layer_absorbed
                END IF ! (printlev_loc>=4)



                eta_4(:) = light_layer_absorbed 
                eta_4_surf = light_surf_absorbed


                !! Debug:

                IF (printlev_loc>=4) THEN
                   WRITE (numout, *) 'leaf level absorption eta_4:', eta_4(:)
                   WRITE (numout, *) 'soil surface absorption eta_4_surf:', eta_4_surf
                END IF ! (printlev_loc>=4) THEN


             END IF ! (hack_can_structure)
             !+++++++++++

          END IF ! (jnlvls doesn't equal 1)    



          !! Beginning of the calculation of all the coefficients all along the canopy:
          !! In order to summarize the principle, here are some explanations:
          !!    - For each layer, the temperature and humidity are calculated as an expression of:
          !!          - The temperature/humidity of the air of the layer below at the same time step;
          !!          - The temperature/humidity of the leaves at the same time step;
          !!          - A set of coefficients (A_t, B_t, C_t and D_t for temperature and
          !!            A_q, B_q, C_q and D_q for humidity)
          !!    - The coefficients are calculated thanks to the values of the coefficients of the
          !!      layer above;
          !!    - The leaf temperature is also an expression of air temperature/humidity of the layer below
          !!      and a set of coefficients based on their values in the upper layer.
          !!
          !! --> The technical note about energy budgets in ORCHIDEE explains this in more details.
          !!
          !! Consequently, as the coefficients in the first layer of the atmosphere are known (petAcoef...),
          !! the coefficients are all expressed one by one going from the top of the canopy until the surface.
          !! As soon as all those coefficients are calculated, the temperatures and humidities are calculated
          !! from the bottom to the top.


          !! Definition of the coefficients :

          IF (jnlvls .eq. 1) THEN !! Case where there is only one layer.

             x_1(1) = 0.0d0; x_2(1) = 0.0d0; x_3(1) = 0.0d0; x_4(1) = 0.0d0; x_5(1) = 0.0d0

             y_1(1) = 0.0d0; y_2(1) = 0.0d0; y_3(1) = 0.0d0; y_4(1) = 0.0d0; y_5(1) = 0.0d0

             !! TO DO: remove this useless loop

             DO i = jnlvls, jnlvls

                big_e_top(i) = (-(dt_sechiba*a_q(i)*chalev0*rau(ji))/                                &
                     &  (rho_veg*leaf_tks*pad_deltah(i)*delta_h(i)*                                 &
                     &  big_r_prime(i)*jtheta(i)))                                                  &
                     &  - ((dt_sechiba* d_t(i) * cp_air * rau(ji)) /                                   & 
                     & (rho_veg*leaf_tks * pad_deltah(i) * delta_h(i) * big_r(i) * jtheta(i)))   

                big_f_top(i) = (-(dt_sechiba*a_t(i)*cp_air*rau(ji))/                                 &
                     & (rho_veg*leaf_tks*pad_deltah(i)*delta_h(i)*                                  &
                     & big_r(i)*jtheta(i)))                                                         &
                     & -((dt_sechiba* d_q(i) * chalev0 * rau(ji)) /                                    &
                     & (rho_veg*leaf_tks * pad_deltah(i)* delta_h(i) * big_r_prime(i) * jtheta(i)))

                big_g_top(i) = (t_l_pres(i) + dt_sechiba* eta_3(i)/                                  &
                     & (rho_veg*leaf_tks * pad_deltah(i) * delta_h(i) * jtheta(i))                  &  
                     &  + (eta_1(i) * lwdown(ji) * dt_sechiba)/                                        &
                     & (rho_veg*leaf_tks * pad_deltah(i) * delta_h(i) * jtheta(i))                  &
                     &  + (eta_4(i) * swdown(ji) * dt_sechiba)/                                        &
                     & (rho_veg*leaf_tks * pad_deltah(i) * delta_h(i) * jtheta(i))                  &
                     &  + (chalev0 * rau(ji) * dt_sechiba* jbeta(i))/                                  &
                     & (rho_veg*leaf_tks * pad_deltah(i) *delta_h(i) *  big_r_prime(i) * jtheta(i)) &
                     &  - (dt_sechiba* cp_air * rau(ji) * b_t(i))/                                     &
                     & (rho_veg*leaf_tks * pad_deltah(i) * delta_h(i) * big_r(i) * jtheta(i))       &
                     &  - (dt_sechiba* chalev0 * rau(ji) * b_q(i))/                                    &
                     & (rho_veg*leaf_tks * pad_deltah(i) * delta_h(i) * big_r_prime(i) * jtheta(i)))

                big_efg_bottom(i) =                                                               &
                     & ( 1.0d0 - (dt_sechiba* jalpha(i) * chalev0 * rau(ji))/                          &
                     &  (rho_veg*leaf_tks * pad_deltah(i) * delta_h(i) * big_r_prime(i)*jtheta(i))  &
                     &  - (dt_sechiba* cp_air * rau(ji)) /                                             &
                     & (rho_veg*leaf_tks * pad_deltah(i) * delta_h(i) * jtheta(i) * big_r(i))       &
                     &  - (eta_2(i) * dt_sechiba) /                                                    & 
                     & (rho_veg *leaf_tks * pad_deltah(i) * delta_h(i) * jtheta(i))                 &
                     &  + (dt_sechiba* chalev0 * rau(ji) * cp_air * c_t(i)) /                          &
                     & (rho_veg *leaf_tks * pad_deltah(i) * delta_h(i) *                            &
                     &  jtheta(i) * big_r(i)) )

                big_e(i) = big_e_top(i) / big_efg_bottom(i)

                big_f(i) = big_f_top(i) / big_efg_bottom(i)

                big_g(i) = big_g_top(i) / big_efg_bottom(i)

             END DO ! i = jnlvls, jnlvls

          ELSE ! Case with several layers

             !! Calculation of the coefficients in the highest layer.

             !! TO DO: Remove this useless loop

             DO i = jnlvls, jnlvls

                !! TO DO: Check if this IF condition is useless or not (no matter if we keep the flag or not)
                !!        --> What leads to the distinction between hydrol_flag=TRUE or FALSE?

                IF (hydrol_flag) THEN

                   big_e_top(i) =                                                                     &
                        &   ((dt_sechiba* d_t(i) * cp_air * rau(ji)) /                                    & 
                        & (rho_veg*leaf_tks * pad_deltah(i) * delta_h(i) * jtheta(i) * big_r(i)))   

                   big_f_top(i) = ( + (dt_sechiba * a_t(i) * cp_air * rau(ji))/                          &
                        & (rho_veg*leaf_tks * pad_deltah(i)* delta_h(i)* jtheta(i) * big_r(i)) )    

                   big_g_top(i) = (t_l_pres(i) + dt_sechiba* eta_3(i)/                                   &
                        & (rho_veg*leaf_tks * pad_deltah(i) * delta_h(i) * jtheta(i))                  & 
                        &  + (eta_1(i) * lwdown(ji) * dt_sechiba)/                                        &
                        & (rho_veg*leaf_tks * pad_deltah(i) * delta_h(i) * jtheta(i))                  &
                        &  + (eta_4(i) * swdown(ji) * dt_sechiba)/                                        & 
                        & (rho_veg*leaf_tks * pad_deltah(i) * delta_h(i) * jtheta(i))                  &
                        &  + (dt_sechiba * cp_air * rau(ji) * b_t(i))/                                    &
                        & (rho_veg*leaf_tks * pad_deltah(i) * delta_h(i) * jtheta(i) * big_r(i) )      &
                        &  - (dt_sechiba * phi_trans(i)) /                                                &
                        & (rho_veg*leaf_tks * pad_deltah(i) * delta_h(i) * jtheta(i)))

                   big_efg_bottom(i) =                                                                 &
                        & ( 1.0d0                                                                       &
                        &  + (dt_sechiba * cp_air * rau(ji)) /                                             &
                        & (rho_veg*leaf_tks * pad_deltah(i) * delta_h(i) * jtheta(i) * big_r(i))        &
                        &  - (eta_2(i) * dt_sechiba) /                                                     & 
                        & (rho_veg*leaf_tks * pad_deltah(i) * delta_h(i) * jtheta(i))                   &
                        &  - (dt_sechiba * chalev0 * rau(ji) * cp_air * c_t(i)) /                          &
                        & (rho_veg*leaf_tks * pad_deltah(i) * delta_h(i) * jtheta(i) * big_r(i)) )

                   big_e(i) = big_e_top(i) / big_efg_bottom(i)

                   big_f(i) = big_f_top(i) / big_efg_bottom(i)

                   big_g(i) = big_g_top(i) / big_efg_bottom(i)

                ELSE ! following if not (hydrol_flag)

                   big_e_top(i) = (+(dt_sechiba*a_q(i)*chalev0*rau(ji))/                                 &
                        &  (rho_veg*leaf_tks*pad_deltah(i)*delta_h(i)*                                 &
                        &  big_r_prime(i)*jtheta(i)))                                                  &
                        &  + ((dt_sechiba* d_t(i) * cp_air * rau(ji)) /                                   & 
                        & (rho_veg*leaf_tks * pad_deltah(i) * delta_h(i) * big_r(i) * jtheta(i)))   

                   big_f_top(i) = (+(dt_sechiba*a_t(i)*cp_air*rau(ji))/                                  &
                        & (rho_veg*leaf_tks*pad_deltah(i)*delta_h(i)*                                  &
                        & big_r(i)*jtheta(i)))                                                         &
                        & +((dt_sechiba* d_q(i) * chalev0 * rau(ji)) /                                    & 
                        & (rho_veg*leaf_tks * pad_deltah(i)* delta_h(i) * big_r_prime(i) * jtheta(i)))

                   big_g_top(i) = (t_l_pres(i) + dt_sechiba* eta_3(i)/                                   &
                        & (rho_veg*leaf_tks * pad_deltah(i) * delta_h(i) * jtheta(i))                  & 
                        &  + (eta_1(i) * lwdown(ji) * dt_sechiba)/                                        &
                        & (rho_veg*leaf_tks * pad_deltah(i) * delta_h(i) * jtheta(i))                  &
                        &  + (eta_4(i) * swdown(ji) * dt_sechiba)/                                        & 
                        & (rho_veg*leaf_tks * pad_deltah(i) * delta_h(i) * jtheta(i))                  &
                        &  - (chalev0 * rau(ji) * dt_sechiba* jbeta(i))/                                  &
                        & (rho_veg*leaf_tks * pad_deltah(i) *delta_h(i) *  big_r_prime(i) * jtheta(i)) &
                        &  + (dt_sechiba* cp_air * rau(ji) * b_t(i))/                                     &
                        & (rho_veg*leaf_tks * pad_deltah(i) * delta_h(i) * big_r(i) * jtheta(i))       &
                        &  + (dt_sechiba* chalev0 * rau(ji) * b_q(i))/                                    & 
                        & (rho_veg*leaf_tks * pad_deltah(i) * delta_h(i) * big_r_prime(i) * jtheta(i)))

                   big_efg_bottom(i) =                                                                &
                        & ( 1.0d0 + (dt_sechiba* jalpha(i) * chalev0 * rau(ji))/                          &
                        &  (rho_veg*leaf_tks * pad_deltah(i) * delta_h(i) * big_r_prime(i)*jtheta(i))  &
                        &  + (dt_sechiba* cp_air * rau(ji)) /                                             &
                        & (rho_veg*leaf_tks * pad_deltah(i) * delta_h(i) * jtheta(i) * big_r(i))       &
                        &  - (eta_2(i) * dt_sechiba) /                                                    & 
                        & (rho_veg *leaf_tks * pad_deltah(i) * delta_h(i) * jtheta(i))                 &
                        &  - (dt_sechiba* chalev0 * rau(ji) * cp_air * c_t(i)) /                          &
                        & (rho_veg *leaf_tks * pad_deltah(i) * delta_h(i) *                            &
                        &  jtheta(i) * big_r(i)) )

                   big_e(i) = big_e_top(i) / big_efg_bottom(i)

                   big_f(i) = big_f_top(i) / big_efg_bottom(i)

                   big_g(i) = big_g_top(i) / big_efg_bottom(i)

                END IF ! (hydrol_flag)

             END DO ! i = jnlvls, jnlvls

             
             !! Calculation of the coefficients from the top the bottom of the surface:

             DO i=jnlvls-1, 2, -1

                ! --------------------------------------------------------------------

                !! TO DO: Check if the flag is needed or not

                IF (hydrol_flag) THEN


                   ! -------------------------------------------------------------------

                   x_1(i) = 1 - ( dt_sechiba* ( (a_t(i+1) * (k_eddy(i)/(delta_z(i)*delta_h(i))))  &
                        &           - (k_eddy(i)/(delta_z(i)*delta_h(i)))                  &
                        &           - (k_eddy(i-1)/(delta_z(i-1)*delta_h(i)))              &
                        &           - ((1.0d0)/(big_r(i)*delta_h(i)))))                    &
                        &           - ((k_eddy(i) * c_t(i+1)) /                            &
                        &             (delta_z(i)*delta_h(i))) * big_f(i+1) * dt_sechiba

                   x_2(i) = t_a_pres(i) + ((b_t(i+1)*k_eddy(i)*dt_sechiba)/(delta_z(i) * delta_h(i)))     &
                        & + ((k_eddy(i) * c_t(i+1)) / (delta_z(i) * delta_h(i))) * big_g(i+1) * dt_sechiba

                   x_3(i) = (k_eddy(i-1)/(delta_z(i-1)*delta_h(i))) * dt_sechiba

                   x_4(i) = (k_eddy(i) * d_t(i+1) / (delta_z(i) * delta_h(i))) * dt_sechiba               &
                        &  + ((k_eddy(i) * c_t(i+1)) / (delta_z(i) * delta_h(i))) * big_e(i+1) * dt_sechiba

                   x_5(i) = + (dt_sechiba) / (delta_h(i) * big_r(i))

                   ! -------------------------------------------------------------------

                   y_1(i) = 1 - ( dt_sechiba* ( (a_q(i+1) * k_eddy(i)/(delta_z(i)*delta_h(i)))       &
                        &  - (k_eddy(i)/(delta_z(i)*delta_h(i)))                             &
                        &  - (k_eddy(i-1)/(delta_z(i-1)*delta_h(i)))                         &
                        &  ))  

                   y_2(i) = q_a_pres(i) +                                                                 &
                        & ((b_q(i+1) * k_eddy(i) * dt_sechiba)/(delta_z(i) * delta_h(i)))                     &
                        & + ((k_eddy(i) * c_q(i+1)) / (delta_z(i) * delta_h(i))) * big_g(i+1) * dt_sechiba    &
                        & + (dt_sechiba * phi_trans(i)) / (delta_h(i) * chalev0 * rau(ji))

                   y_3(i) = (k_eddy(i-1)/(delta_z(i-1)*delta_h(i))) * dt_sechiba



                   y_4(i) = (k_eddy(i) * d_q(i+1) / (delta_z(i) * delta_h(i))) * dt_sechiba          

                   y_5(i) = 0.0d0 ! + (jalpha(i) * dt_sechiba) / (delta_h(i) * big_r_prime(i))

                   ! -------------------------------------------------------------------
                   a_t(i) = x_3(i) / (x_1(i) - (x_4(i) * (y_4(i) / y_1(i))))
                   b_t(i) = ( x_2(i) + (x_4(i) * (y_2(i)/y_1(i))) )  /        &
                        & (x_1(i) - (x_4(i) *(y_4(i)/y_1(i))))
                   c_t(i) = ( (x_4(i) * (y_5(i)/y_1(i))) + x_5(i)) /          &
                        & (x_1(i) - (x_4(i)*(y_4(i)/y_1(i))))
                   d_t(i) = ( x_4(i) * (y_3(i)/y_1(i))  )   /                 &
                        & (x_1(i) - (x_4(i) * (y_4(i) / y_1(i))))

                   !! Debug:

                   IF(printlev_loc >= 4) THEN
                      WRITE(numout,*)'hydrol_flag TRUE'      
                      WRITE(numout,*)'a_t(i)', a_t(i)
                      WRITE(numout,*)'b_t(i)', b_t(i)
                   ENDIF

                   ! -------------------------------------------------------------------

                   a_q(i) = y_3(i) / (y_1(i) - (y_4(i) * (x_4(i) / x_1(i))))
                   b_q(i) = ( y_2(i) + (y_4(i) * (x_2(i)/x_1(i))) )  /        &
                        & (y_1(i) - (y_4(i) *(x_4(i)/x_1(i))))
                   c_q(i) = ( (y_4(i) * (x_5(i)/x_1(i))) + y_5(i)) /          &
                        & (y_1(i) - (y_4(i)*(x_4(i)/x_1(i))))
                   d_q(i) = ( y_4(i) * (x_3(i)/x_1(i))  )   /                 &
                        & (y_1(i) - (y_4(i) * (x_4(i) / x_1(i))))

                   ! -------------------------------------------------------------------


                   big_e_top(i) =                                                                        &
                        &   ((dt_sechiba* d_t(i) * cp_air * rau(ji)) /                                    & 
                        & (rho_veg*leaf_tks * pad_deltah(i) * delta_h(i)* jtheta(i) * big_r(i)))   

                   big_f_top(i) = ( + (dt_sechiba* a_t(i) * cp_air * rau(ji))/                              &
                        & (rho_veg*leaf_tks * pad_deltah(i)* delta_h(i) * jtheta(i)* big_r(i)) )    


                   big_g_top(i) = (t_l_pres(i) + dt_sechiba* eta_3(i)/                                      &
                        & (rho_veg*leaf_tks*pad_deltah(i)*delta_h(i)*jtheta(i))                        & 
                        &  + (eta_1(i) * lwdown(ji) * dt_sechiba)/                                        &
                        & (rho_veg*leaf_tks*pad_deltah(i)*delta_h(i)*jtheta(i))                        &
                        &  + (eta_4(i) * swdown(ji) * dt_sechiba)/                                        & 
                        & (rho_veg*leaf_tks*pad_deltah(i)*delta_h(i)*jtheta(i))                        &
                        &  + (dt_sechiba * cp_air * rau(ji) * b_t(i))/                                    &
                        & (rho_veg*leaf_tks*pad_deltah(i)*delta_h(i)* jtheta(i) * big_r(i) )           &
                                !  &  - (dt_sechiba * phi_trans) / (delta_h(i) * jtheta(i) * rho_veg))
                        &  - (dt_sechiba * phi_trans(i)) /                                                &
                        & (rho_veg*leaf_tks * pad_deltah(i) * delta_h(i) * jtheta(i)))


                   big_efg_bottom(i) =                                                                  &
                        & ( 1.0d0                                                                      &
                        &  + (dt_sechiba * cp_air * rau(ji)) /                                            &
                        & (rho_veg*leaf_tks * pad_deltah(i) * delta_h(i) * jtheta(i)  * big_r(i))      &
                        &  - (eta_2(i) * dt_sechiba) /                                                    & 
                        & (rho_veg *leaf_tks * pad_deltah(i) * delta_h(i) * jtheta(i))                 &
                        &  - (dt_sechiba * cp_air * rau(ji) * c_t(i)) /                                   &
                        & (rho_veg *leaf_tks * pad_deltah(i) * delta_h(i) *jtheta(i) * big_r(i)) )


                   big_e(i) = big_e_top(i) / big_efg_bottom(i)

                   big_f(i) = big_f_top(i) / big_efg_bottom(i)

                   big_g(i) = big_g_top(i) / big_efg_bottom(i)


                ELSE ! following if not (hydrol_flag)
                   
                   !! Debug:

                   IF(printlev_loc >= 4 ) THEN
                      WRITE(numout,*) 'dt_sechiba',dt_sechiba
                      WRITE(numout,*) 'a_t(i+1)',a_t(i+1)
                      WRITE(numout,*) 'k_eddy(i)', k_eddy(i)
                      WRITE(numout,*) 'delta_z(i)',delta_z(i)
                      WRITE(numout,*) 'delta_h(i)',delta_h(i)
                      WRITE(numout,*) 'k_eddy(i-1)',k_eddy(i-1)
                      WRITE(numout,*) 'delta_z(i-1)', delta_z(i-1)
                      WRITE(numout,*) 'big_r(i)', big_r(i)
                      WRITE(numout,*) 'c_t(i+1)', c_t(i+1)
                      WRITE(numout,*) 'big_f(i+1)', big_f(i+1)
                   ENDIF ! (printlev_lov >=  4)

                   ! --------------------------------------------------------------------

                   x_1(i) = 1 - ( dt_sechiba* ( (a_t(i+1) * (k_eddy(i)/(delta_z(i)*delta_h(i))))       &
                        &           - (k_eddy(i)/(delta_z(i)*delta_h(i)))                      &
                        &           - (k_eddy(i-1)/(delta_z(i-1)*delta_h(i)))                  &
                        &           - ((1.0d0)/(big_r(i)*delta_h(i)))))                        &
                        &           - ((k_eddy(i) * c_t(i+1)) /                                &
                        &             (delta_z(i)*delta_h(i))) * big_f(i+1) * dt_sechiba

                   x_2(i) = t_a_pres(i) + ((b_t(i+1)*k_eddy(i)*dt_sechiba)/(delta_z(i) * delta_h(i)))  &
                        & + ((k_eddy(i) * c_t(i+1)) / (delta_z(i) * delta_h(i))) * big_g(i+1) * dt_sechiba

                   x_3(i) = (k_eddy(i-1)/(delta_z(i-1)*delta_h(i))) * dt_sechiba

                   x_4(i) = (k_eddy(i) * d_t(i+1) / (delta_z(i) * delta_h(i))) * dt_sechiba            &
                        &  + ((k_eddy(i) * c_t(i+1)) / (delta_z(i) * delta_h(i))) * big_e(i+1) * dt_sechiba

                   x_5(i) = + (dt_sechiba) / (delta_h(i) * big_r(i))

                   ! -------------------------------------------------------------------

                   y_1(i) = 1.0d0 - ( dt_sechiba* ( (a_q(i+1) * (k_eddy(i)/(delta_z(i)*delta_h(i))))   &
                        &  - (k_eddy(i)/(delta_z(i)*delta_h(i)))                                      &
                        &  - (k_eddy(i-1)/(delta_z(i-1)*delta_h(i)))                                  &
                        &  - ((1.0d0)/(big_r_prime(i)*delta_h(i)))))                                  &
                        &  - ((k_eddy(i) * c_q(i+1)) / (delta_z(i)*delta_h(i))) * big_e(i+1) * dt_sechiba

                   y_2(i) = q_a_pres(i) +                                                           &
                        & ((b_q(i+1) * k_eddy(i) * dt_sechiba)/(delta_z(i) * delta_h(i)))                &
                        & + (jbeta(i) * dt_sechiba)  / (big_r_prime(i)*delta_h(i))                       &
                        & + ((k_eddy(i) * c_q(i+1)) / (delta_z(i) * delta_h(i))) * big_g(i+1) * dt_sechiba

                   y_3(i) = (k_eddy(i-1)/(delta_z(i-1)*delta_h(i))) * dt_sechiba

                   y_4(i) = (k_eddy(i) * d_q(i+1) / (delta_z(i) * delta_h(i))) * dt_sechiba            &
                        &  + ((k_eddy(i) * c_q(i+1)) / (delta_z(i) * delta_h(i))) * big_f(i+1) * dt_sechiba

                   y_5(i) = + (jalpha(i) * dt_sechiba) / (delta_h(i) * big_r_prime(i))

                   ! -------------------------------------------------------------------

                   a_t(i) = x_3(i) / (x_1(i) - (x_4(i) * (y_4(i) / y_1(i))))
                   b_t(i) = ( x_2(i) + (x_4(i) * (y_2(i)/y_1(i))) )  /        &
                        & (x_1(i) - (x_4(i) *(y_4(i)/y_1(i))))
                   c_t(i) = ( (x_4(i) * (y_5(i)/y_1(i))) + x_5(i)) /          &
                        & (x_1(i) - (x_4(i)*(y_4(i)/y_1(i))))
                   d_t(i) = ( x_4(i) * (y_3(i)/y_1(i))  )   /                 &
                        & (x_1(i) - (x_4(i) * (y_4(i) / y_1(i))))

                   !! Debug:

                   IF(printlev_loc >= 4) THEN 
                      WRITE(numout,*)'hydrol_flag FALSE'
                      WRITE(numout,*)'a_t(i)', a_t(i)
                      WRITE(numout,*)'b_t(i)', b_t(i)
                   ENDIF
                   ! -------------------------------------------------------------------

                   a_q(i) = y_3(i) / (y_1(i) - (y_4(i) * (x_4(i) / x_1(i))))
                   b_q(i) = ( y_2(i) + (y_4(i) * (x_2(i)/x_1(i))) )  /        &
                        & (y_1(i) - (y_4(i) *(x_4(i)/x_1(i))))
                   c_q(i) = ( (y_4(i) * (x_5(i)/x_1(i))) + y_5(i)) /          &
                        & (y_1(i) - (y_4(i)*(x_4(i)/x_1(i))))
                   d_q(i) = ( y_4(i) * (x_3(i)/x_1(i))  )   /                 &
                        & (y_1(i) - (y_4(i) * (x_4(i) / x_1(i))))

                   ! -------------------------------------------------------------------

                   big_e_top(i) = (+(dt_sechiba*a_q(i)*chalev0*rau(ji))/                             &
                        & (rho_veg*leaf_tks * pad_deltah(i)*delta_h(i)*big_r_prime(i)*jtheta(i)))  &
                        & + ((dt_sechiba* d_t(i) * cp_air * rau(ji)) /                                &
                        & (rho_veg*leaf_tks * pad_deltah(i) * delta_h(i) * big_r(i) * jtheta(i)))   

                   big_f_top(i) = ( + (dt_sechiba* a_t(i) * cp_air * rau(ji))/                       &
                        & (rho_veg*leaf_tks * pad_deltah(i)* delta_h(i)                            &
                        & * big_r(i) * jtheta(i)) )                                                &
                        & + ((dt_sechiba* d_q(i) * chalev0 * rau(ji)) /                               &
                        & (rho_veg*leaf_tks * pad_deltah(i) * delta_h(i)                           &
                        & * big_r_prime(i) * jtheta(i)))

                   big_g_top(i) = (t_l_pres(i)+dt_sechiba*eta_3(i)/                                  &
                        & (rho_veg*leaf_tks*pad_deltah(i)*delta_h(i)*jtheta(i))                    & 
                        & + (eta_1(i)*lwdown(ji)*dt_sechiba)/                                         &
                        & (rho_veg*leaf_tks*pad_deltah(i)*delta_h(i)*jtheta(i))                    &
                        & + (eta_4(i)*swdown(ji)*dt_sechiba)/                                         &
                        & (rho_veg*leaf_tks*pad_deltah(i)*delta_h(i)*jtheta(i))                    &
                        & - (chalev0*rau(ji)*dt_sechiba*jbeta(i))/                                    &
                        & (rho_veg*leaf_tks*pad_deltah(i)*delta_h(i) *                             &
                        &  big_r_prime(i)*jtheta(i))                                               &
                        & + (dt_sechiba*cp_air*rau(ji)*b_t(i))/                                       &
                        & (rho_veg*leaf_tks*pad_deltah(i)*delta_h(i)*big_r(i)* jtheta(i))          &
                        & + (dt_sechiba*chalev0*rau(ji)*b_q(i))/                                      &
                        & (rho_veg*leaf_tks*pad_deltah(i)*delta_h(i)*                              &
                        & big_r_prime(i)*jtheta(i)) )

                   big_efg_bottom(i) =                                                            &
                        & ( 1.0d0 + (dt_sechiba* jalpha(i) * chalev0 * rau(ji)) /                     &
                        & (rho_veg *leaf_tks * pad_deltah(i) * delta_h(i) *                        &
                        & big_r_prime(i) * jtheta(i))                                              &
                        & + (dt_sechiba* cp_air * rau(ji)) /                                          &
                        & (rho_veg *leaf_tks * pad_deltah(i) * delta_h(i) * jtheta(i) * big_r(i))  &
                        & - (eta_2(i) * dt_sechiba) /                                                 &
                        & (rho_veg*leaf_tks * pad_deltah(i)  * delta_h(i) * jtheta(i))             &
                        & - (dt_sechiba* chalev0 * rau(ji) * c_q(i)) /                                & 
                        & (rho_veg *leaf_tks* pad_deltah(i) * delta_h(i) * jtheta(i)               &
                        & * big_r_prime(i))                                                        &
                        & - (dt_sechiba* cp_air * rau(ji) * c_t(i))/                                  &
                        & (rho_veg *leaf_tks * pad_deltah(i) * delta_h(i) * jtheta(i) * big_r(i)))

                   big_e(i) = big_e_top(i) / big_efg_bottom(i)

                   big_f(i) = big_f_top(i) / big_efg_bottom(i)

                   big_g(i) = big_g_top(i) / big_efg_bottom(i)

                END IF ! (hydrol_flag)               
                ! -------------------------------------------------------------------

             END DO !i = jnlvls-1, 2, -1

             

             !! Calculation of the coefficients in the first layer (surface layer):
             
             !! TO DO: Check if the IF condition is needed or not.

             IF (hydrol_flag) THEN

                ! -------------------------------------------------------------------

                x_1(1) = 1.0d0 - dt_sechiba* ((a_t(2) * k_eddy(1)) / (delta_z(1) * delta_h(1))     &
                     &  +  (k_eddy(1) * c_t(2) * big_f(2)) / (delta_z(1) * delta_h(1))       &
                     &  -  (k_eddy(1)/ (delta_z(1) * delta_h(1)))                            &
                     &  -  ( 1.0d0 / (delta_h(1) * big_r(1))) )

                x_2(1) = t_a_pres(1)+dt_sechiba*(  (k_eddy(1) * b_t(2))/(delta_z(1) * delta_h(1))  &
                     &   +  (k_eddy(1) * c_t(2) * big_g(2))/(delta_z(1) * delta_h(1)) )

                x_3(1) = - dt_sechiba* (1.0d0 / delta_h(1)) * ( 1.0d0 / (rau(ji) * cp_air) )    

                x_4(1) = dt_sechiba* ( (k_eddy(1) * c_t(2) * big_e(2)) + (k_eddy(1) * d_t(2)) ) /  &
                     & (delta_z(1) * delta_h(1)) 

                x_5(1) = + dt_sechiba* ( 1.0d0 / (delta_h(1) * big_r(1)) )

                ! -------------------------------------------------------------------

                y_1(1) = 1.0d0 - dt_sechiba* (    (a_q(2) * k_eddy(1)) / (delta_z(1)*delta_h(1))   &
                     &  +  (k_eddy(1) * c_q(2) * big_e(2)) / (delta_z(1) * delta_h(1))       &
                     &  -  (k_eddy(1) / (delta_z(1) * delta_h(1)))                           &
                     &  -  ( 1.0d0 / (delta_h(1) * big_r_prime(1))) )

                y_2(1) = q_a_pres(1)+dt_sechiba*(  (k_eddy(1) * b_q(2))/(delta_z(1) * delta_h(1))  &
                     &   +  (k_eddy(1) * c_q(2) * big_g(2))/(delta_z(1) * delta_h(1))        &
                     &   + jbeta(1) / (delta_h(1) * big_r_prime(1)) )

                y_3(1) = - dt_sechiba* (1.0d0 / delta_h(1)) * ( 1.0d0 / (rau(ji) * chalev0) )

                y_4(1) = dt_sechiba* ( (k_eddy(1)* c_q(2) * big_f(2)) + (k_eddy(1) * d_q(2)) ) /  &
                     &   (delta_z(1) * delta_h(1)) 

                y_5(1) = + dt_sechiba* ( jalpha(1) / (delta_h(1) * big_r_prime(1)) )

                ! -------------------------------------------------------------------

                a_t(1) = x_3(1) / (x_1(1) - (x_4(1) * (y_4(1) / y_1(1))))
                b_t(1) = ( x_2(1) + (x_4(1) * (y_2(1)/y_1(1))) )  /         &
                     & (x_1(1) - (x_4(1) *(y_4(1)/y_1(1))))
                c_t(1) = ( (x_4(1) * (y_5(1)/y_1(1))) + x_5(1)) /           &
                     & (x_1(1) - (x_4(1)*(y_4(1)/y_1(1))))
                d_t(1) = ( x_4(1) * (y_3(1)/y_1(1))  )   /                  &
                     & (x_1(1) - (x_4(1) * (y_4(1) / y_1(1))))

                !! Debug:

                IF(printlev_loc >= 4) THEN
                   WRITE(numout,*)'hydrol_flag TRUE'
                   WRITE(numout,*)'a_t(i)', a_t(i)
                   WRITE(numout,*)'b_t(i)', b_t(i)
                ENDIF
                ! -------------------------------------------------------------------

                a_q(1) = y_3(1) / (y_1(1) - (y_4(1) * (x_4(1) / x_1(1))))
                b_q(1) = ( y_2(1) + (y_4(1) * (x_2(1)/x_1(1))) )  /         &
                     & (y_1(1) - (y_4(1) *(x_4(1)/x_1(1))))
                c_q(1) = ( (y_4(1) * (x_5(1)/x_1(1))) + y_5(1)) /           & 
                     & (y_1(1) - (y_4(1)*(x_4(1)/x_1(1))))
                d_q(1) = ( y_4(1) * (x_3(1)/x_1(1))  )   /                  &
                     & (y_1(1) - (y_4(1) * (x_4(1) / x_1(1))))

             ELSE ! following if not (hydrol_flag)

                x_1(1) = 1.0d0 - dt_sechiba* ((a_t(2) * k_eddy(1)) / (delta_z(1)*delta_h(1))       &
                     &  +  (k_eddy(1) * c_t(2) * big_f(2)) / (delta_z(1) * delta_h(1))       &
                     &  -  (k_eddy(1) / (delta_z(1) * delta_h(1)))                           &
                     &  -  ( 1.0d0 / (delta_h(1) * big_r(1))) )

                x_2(1) = t_a_pres(1)+dt_sechiba*(  (k_eddy(1) * b_t(2))/(delta_z(1) * delta_h(1))  &
                     &   +  (k_eddy(1) * c_t(2) * big_g(2))/(delta_z(1) * delta_h(1)) )

                x_3(1) = - dt_sechiba* (1.0d0 / delta_h(1)) * ( 1.0d0 / (rau(ji) * cp_air) )    

                x_4(1) = dt_sechiba* ( (k_eddy(1) * c_t(2) * big_e(2)) + (k_eddy(1) * d_t(2)) ) /  &
                     & (delta_z(1) * delta_h(1)) 

                x_5(1) = + dt_sechiba* ( 1.0d0 / (delta_h(1) * big_r(1)) )

                ! -------------------------------------------------------------------

                y_1(1) = 1.0d0 - dt_sechiba* (    (a_q(2) * k_eddy(1)) / (delta_z(1)*delta_h(1))   &
                     &  +  (k_eddy(1) * c_q(2) * big_e(2)) / (delta_z(1) * delta_h(1))       &
                     &  -  (k_eddy(1) / (delta_z(1) * delta_h(1)))                           &
                     &  -  ( 1.0d0 / (delta_h(1) * big_r_prime(1))) )

                y_2(1) = q_a_pres(1)+dt_sechiba*(  (k_eddy(1) * b_q(2))/(delta_z(1) * delta_h(1))  &
                     &   +  (k_eddy(1) * c_q(2) * big_g(2))/(delta_z(1) * delta_h(1))        &
                     &   + jbeta(1) / (delta_h(1) * big_r_prime(1)) )

                y_3(1) =  - dt_sechiba* (1.0d0 / delta_h(1)) * ( 1.0d0 / (rau(ji) * chalev0) )

                y_4(1) = dt_sechiba* ( (k_eddy(1) * c_q(2) * big_f(2)) + (k_eddy(1) * d_q(2)) ) /  &
                     &   (delta_z(1) * delta_h(1)) 

                y_5(1) = + dt_sechiba* ( jalpha(1) / (delta_h(1) * big_r_prime(1)) )

                ! -------------------------------------------------------------------

                a_t(1) = x_3(1) / (x_1(1) - (x_4(1) * (y_4(1) / y_1(1))))
                b_t(1) = ( x_2(1) + (x_4(1) * (y_2(1)/y_1(1))) )  /         &
                     & (x_1(1) - (x_4(1) *(y_4(1)/y_1(1))))
                c_t(1) = ( (x_4(1) * (y_5(1)/y_1(1))) + x_5(1)) /           &
                     & (x_1(1) - (x_4(1)*(y_4(1)/y_1(1))))
                d_t(1) = ( x_4(1) * (y_3(1)/y_1(1))  )   /                  &
                     & (x_1(1) - (x_4(1) * (y_4(1) / y_1(1))))

                !! Debug:

                IF(printlev_loc >= 4) THEN
                   WRITE(numout,*)'hydrol_flag FALSE'
                   WRITE(numout,*)'a_t(i)', a_t(i)
                   WRITE(numout,*)'b_t(i)', b_t(i)
                ENDIF
                ! -------------------------------------------------------------------

                a_q(1) = y_3(1) / (y_1(1) - (y_4(1) * (x_4(1) / x_1(1))))
                b_q(1) = ( y_2(1) + (y_4(1) * (x_2(1)/x_1(1))) )  /         &
                     & (y_1(1) - (y_4(1) *(x_4(1)/x_1(1))))
                c_q(1) = ( (y_4(1) * (x_5(1)/x_1(1))) + y_5(1)) /           & 
                     & (y_1(1) - (y_4(1)*(x_4(1)/x_1(1))))
                d_q(1) = ( y_4(1) * (x_3(1)/x_1(1))  )   /                  &
                     & (y_1(1) - (y_4(1) * (x_4(1) / x_1(1))))

                ! -------------------------------------------------------------------

             END IF ! (hydrol_flag)

          END IF ! (jnlvls=1)

          IF (printlev_loc>=4) THEN
             WRITE(numout,*) 'k_eddy', k_eddy(:)
             WRITE(numout,*) 'a_t', a_t(:)
             WRITE(numout,*) 'a_q', a_q(:)
             WRITE(numout,*) 'b_t', b_t(:)
             WRITE(numout,*) 'b_q', b_q(:)
             WRITE(numout,*) 'c_t', c_t(:)
             WRITE(numout,*) 'c_q', c_q(:)
             WRITE(numout,*) 'd_t', d_t(:)
             WRITE(numout,*) 'd_q', d_q(:)
             WRITE(numout,*) 'big_e', big_e(:)
             WRITE(numout,*) 'big_f', big_f(:)
             WRITE(numout,*) 'delta_z', delta_z(:)
             WRITE(numout,*) 'delta_h', delta_h(:)
             WRITE(numout,*) 'j_alpha', jalpha(:)
             WRITE(numout,*) 'big_r', big_r(:)
             WRITE(numout,*) 'big_r_prime', big_r_prime(:)

             WRITE(numout,*) 'x_1', x_1(:)
             WRITE(numout,*) 'x_3', x_3(:)
             WRITE(numout,*) 'x_4', x_4(:)
             WRITE(numout,*) 'y_1', y_1(:)
             WRITE(numout,*) 'y_3', y_3(:)
             WRITE(numout,*) 'y_4', y_4(:)
          ENDIF

          !! Calculation of the last coefficients at the surface layer:

          big_e_top(1) = (+(dt_sechiba* a_q(1) * chalev0 * rau(ji)) /                              &
               & (rho_veg *leaf_tks * pad_deltah(1) * delta_h(1) * jtheta(1) * big_r_prime(1)))   &
               & +((dt_sechiba* d_t(1) * cp_air * rau(ji)) /                                         &
               & (rho_veg*leaf_tks * pad_deltah(1) * delta_h(1) * jtheta(1) * big_r(1)) )   


          big_f_top(1) = (+(dt_sechiba* a_t(1) * cp_air * rau(ji)) /                              &
               & (rho_veg*leaf_tks * pad_deltah(1) * delta_h(1)* jtheta(1) * big_r(1)))          &
               & + ((dt_sechiba* d_q(1) * chalev0 * rau(ji)) /                                      &
               & (rho_veg*leaf_tks * pad_deltah(1) * delta_h(1) * jtheta(1) * big_r_prime(1)))


          big_g_top(1) = (t_l_pres(1) + dt_sechiba* eta_3(1) /                                    &
               & (rho_veg*leaf_tks * pad_deltah(1) * delta_h(1) * jtheta(1))                     & 
               & + (eta_1(1) * lwdown(ji) * dt_sechiba)/                                            &
               & (rho_veg *leaf_tks* pad_deltah(1) * delta_h(1) * jtheta(1))                     &
               &  + (eta_4(1) * swdown(ji) * dt_sechiba)/                                           &
               & (rho_veg *leaf_tks* pad_deltah(1) * delta_h(1) * jtheta(1))                     &
               &  - (chalev0 * rau(ji) * dt_sechiba* jbeta(1))/                                     & 
               & (rho_veg *leaf_tks* pad_deltah(1) * delta_h(1) * jtheta(1) * big_r_prime(1))    &
               &  + (dt_sechiba* cp_air * rau(ji) * b_t(1)) /                                       &
               & (rho_veg *leaf_tks* pad_deltah(1) * delta_h(1) * jtheta(1) * big_r(1))          &
               &  + (dt_sechiba * chalev0 * rau(ji) * b_q(1)) /                                     &
               & (rho_veg *leaf_tks * pad_deltah(1) * delta_h(1) * jtheta(1) * big_r_prime(1)))


          big_efg_bottom(1) = ( 1.0d0 + (dt_sechiba *jalpha(1) * chalev0 * rau(ji))/              &
               & (rho_veg *leaf_tks * pad_deltah(1) * delta_h(1) * jtheta(1) * big_r_prime(1) )  &
               & + (dt_sechiba * cp_air * rau(ji)) /                                                &
               & (rho_veg *leaf_tks * pad_deltah(1) * delta_h(1) * jtheta(1) * big_r(1))         &
               & - (eta_2(1) * dt_sechiba) /                                                        &
                                ! & - (0.0d0 * dt_sechiba) /                                                         &
               & (rho_veg*leaf_tks * pad_deltah(1)  * delta_h(1) * jtheta(1))                    &
               & - (dt_sechiba * chalev0 * rau(ji) * c_q(1)) /                                      &
               & (rho_veg *leaf_tks * pad_deltah(1) * delta_h(1) * jtheta(1) * big_r_prime(1))   &
               & - (dt_sechiba * cp_air * rau(ji) * c_t(1))/                                        &
               & (rho_veg *leaf_tks * pad_deltah(1) * delta_h(1) * jtheta(1) * big_r(1)))


          big_e(1) = big_e_top(1) / big_efg_bottom(1)

          big_f(1) = big_f_top(1) / big_efg_bottom(1)

          big_g(1) = big_g_top(1) / big_efg_bottom(1)    


          ! --------------------------------------------------------------------------------------------------------------


          !! Solving of the energy budget equation at the surface layer:
          !! The aim is to determine the temperature of the surface layer in order to calculate all the 
          !! temperatures from the bottom of the column to the top.


          IF (jnlvls .eq. 1) THEN
             turb_factor = speed * q_cdrag(ji)                    !! Surface turbulence calculated as in existing ORCHIDEE
          ELSE
             turb_factor = u_speed(ji,1) * q_cdrag(ji) * ks_tune  !! ks_tune=1 
          END IF


          !! TEST JULIEN : Add the beta coefficients before in order to close the scheme, See after for the previous version

          ks_weight_1(ji) = 1.0_r_std / (1.0_r_std + EXP (-ks_slope * (SUM(veget(ji,:)) - ks_veget) )) 
          ks_weight_2(ji) = 1.0_r_std - ks_weight_1(ji) 

          vbeta3_pixel = zero
          DO jv = 1,nvm
             IF (vbeta3(ji,jv).GT.min_sechiba) THEN
                vbeta3_pixel = vbeta3_pixel + 1/(vbeta3(ji,jv))
             ENDIF
          ENDDO


          IF (printlev_loc>=4) THEN
             WRITE(numout,*) 'ks_weight_1', ks_weight_1(1)
             WRITE(numout,*) 'ks_weight_2', ks_weight_2(1)
             WRITE(numout,*) 'veget', veget(1,test_pft)
             WRITE(numout,*) 'vbeta3_mleb', vbeta3(1,test_pft)
             WRITE(numout,*) 'vbeta3_pixel', vbeta3_pixel
             WRITE(numout,*) 'vbeta4', vbeta4
          ENDIF

          ! Calculate the total weighting coefficient both from vbeta4 and vbeta3
          IF (vbeta3_pixel.GT.min_sechiba) THEN
             evap_fraction(ji) =  (  ks_weight_1(ji) * vbeta4(ji) & 
                  &  + ks_weight_2(ji) * (1/vbeta3_pixel) )
          ELSE
             evap_fraction(ji) =  MAX(0.001,ks_weight_1(ji) * vbeta4(ji))
          ENDIF
          !+++++++++++++

          ! evap_fraction(ji) =  vbeta3(ji,jv) 
          ! set up a upper limited for the potential evapotranspiration rate 
          ! thus if the veget ~1.0, there is little forest floor evapotranspiration   
          IF  ( evap_fraction(ji) .GT. 1.0_r_std ) THEN
             evap_fraction(ji) =  1.0_r_std
          ENDIF






          !! Calculation of the components for surface sensible heat flux: 

          big_omega_1 = 1.0d0 - (rau(ji) * 1.0d0) * turb_factor * &
               & (a_t(1) + c_t(1) * big_f(1)) 

          big_omega_2 = + (rau(ji) * cp_air) * turb_factor * &
               & (b_t(1) + c_t(1) * big_g(1))

          big_omega_3 = - rau(ji) * cp_air * turb_factor   

          big_omega_4 = + rau(ji) * cp_air * turb_factor * &
               & (c_t(1) * big_e(1) + d_t(1))


          IF (printlev_loc>=4) THEN
             WRITE(numout,*) 'c_t(1)', c_t(1)
             WRITE(numout,*) 'c_q(1)', c_q(1)
             WRITE(numout,*) 'd_t(1)', d_t(1)
             WRITE(numout,*) 'd_q(1)', d_q(1)
             WRITE(numout,*) 'big_e(1)', big_e(1)
             WRITE(numout,*) 'big_f(1)', big_f(1)
             WRITE(numout,*) 'turb_factor', turb_factor
             WRITE(numout,*) 'rau', rau(1)
             WRITE(numout,*) 'cp_air', cp_air
          ENDIF

          !! Calcualtion of the components for surface latent heat flux (evaporation):


          !! TEST JULIEN : add the beta coefficients

          big_omega_5 = 1.0d0 - evap_fraction(ji)*(rau(ji) * 1.0d0) * turb_factor  * & 
               & (a_q(1) + c_q(1) * big_e(1)) 

          big_omega_6 = +  evap_fraction(ji)*(rau(ji) * chalev0) * turb_factor  * &
               & (b_q(1) + c_q(1) * big_g(1) - jbeta_surf(ji) )

          big_omega_7 = - jalpha_surf(ji) * evap_fraction(ji)* rau(ji) * chalev0 * turb_factor 

          big_omega_8 = +  evap_fraction(ji)*rau(ji) * chalev0 * turb_factor  * &
               & (c_q(1) * big_f(1) + d_q(1))

          !! Calculation of the components for surface latent heat flux (sublimination):


          big_omega_5_sub = 1.0d0 - beta_diffuco_sub(ji) *(rau(ji) * 1.0d0) * turb_factor* & 
               & (a_q(1) + c_q(1) * big_e(1)) 

          big_omega_6_sub = + beta_diffuco_sub(ji) *(rau(ji) * chalsu0) * turb_factor * &
               & (b_q(1) + c_q(1) * big_g(1) - jbeta_surf(ji) )

          big_omega_7_sub = - jalpha_surf(ji) *beta_diffuco_sub(ji) * rau(ji) * chalsu0 * turb_factor  

          big_omega_8_sub = +beta_diffuco_sub(ji) * rau(ji) * chalsu0 * turb_factor *  &
               & (c_q(1) * big_f(1) + d_q(1))

          !! Calculation of the two first components of the surface fluxes:

          xi_1 = ( big_omega_2 + (big_omega_4/big_omega_5) * big_omega_6 )  / &
               &  ( big_omega_1 - (big_omega_4/big_omega_5) * big_omega_8 )

          xi_2 = ( big_omega_3 + (big_omega_4/big_omega_5) * big_omega_7 )  / &
               &  ( big_omega_1 - (big_omega_4/big_omega_5) * big_omega_8 )



          !! TEST JULIEN SEE BEFORE
!!$
!!$          !! Parameterisation of the surface evaporation (see Chen et al. (2016))
!!$
!!$          
!!$          !! Here we try to paramterize the surface evaporation  
!!$          !! with the combination of the vebta3, vbeta4, ks_slope, ks_veget
!!$          !! and other variables such as pgap(=1-veget), day of growth etc  
!!$          !! air_temperature as an indicator to account 
!!$          !! seasonal grass leaf-fall at forest floor     
!!$          
!!$
!!$          !+++CHECK+++
!!$          ! veget and vbeta3 have a PFT-dimension. How should this dimension be accounted for?
!!$          ! Set up a S-shape function (from 0 to 1) based on the site dependent vegetation cover fraction  
!!$          ! for the calculation of the weighting for vebta3 or vbeta4  
!!$          ! paramteriziation the surface evapotranspiration coefficient for the boundary condition for the 
!!$          ! multilayer scheme
!!$          ks_weight_1(ji) = 1.0_r_std / (1.0_r_std + EXP (-ks_slope * (SUM(veget(ji,:)) - ks_veget) )) 
!!$          ks_weight_2(ji) = 1.0_r_std - ks_weight_1(ji) 
!!$
!!$          !++++++++++++
!!$
!!$          ! Add a conditional rule to consider the temperature constrain the growing season
!!$          ! of forest floor when the vegetation cover below certain threshood or critiria             
!!$          ! "ks_veget" and longterm average temperature below certain value (15oC) 
!!$          ! the grass was assume to be reduced having little contribution on total evapotraspiration.  
!!$          !IF ( veget(ji,jv) .LT. ks_veget) THENN
!!$          !	 ks_weight_1(ji) = 1.0_r_std / (((1.0_r_std + abs( t2m_month_out(ji) - 288.15)/15.) ) +  & 
!!$          !				 &               EXP (-ks_slope *(veget(ji,jv) - ks_veget) ))
!!$          !	 ks_weight_2(ji) = 1.0_r_std - ks_weight_1(ji)  
!!$
!!$          !ENDIF     
!!$
!!$          !+++CHECK++++
!!$          ! veget and vbeta3 have a PFT-dimension. How should this dimension be accounted for? vbeta3 is the 
!!$          ! vbeta for each PFT. Use a the analogy of a parallel electric circuit to calculate the total vbeta3 
!!$          ! for the pixel. 
!!$          vbeta3_pixel = zero
!!$          DO jv = 1,nvm
!!$             IF (vbeta3(ji,jv).GT.min_sechiba) THEN
!!$                vbeta3_pixel = vbeta3_pixel + 1/vbeta3(ji,jv)
!!$             ENDIF
!!$          ENDDO
!!$
!!$          ! Calculate the total weighting coefficient both from vbeta4 and vbeta3
!!$          IF (vbeta3_pixel.GT.min_sechiba) THEN
!!$             evap_fraction(ji) =  (  ks_weight_1(ji) * vbeta4(ji) & 
!!$                  &  + ks_weight_2(ji) * (1/vbeta3_pixel) )
!!$          ELSE
!!$             evap_fraction(ji) =  ks_weight_1(ji) * vbeta4(ji)
!!$          ENDIF
!!$          !+++++++++++++
!!$
!!$          ! evap_fraction(ji) =  vbeta3(ji,jv) 
!!$          ! set up a upper limited for the potential evapotranspiration rate 
!!$          ! thus if the veget ~1.0, there is little forest floor evapotranspiration   
!!$          IF  ( evap_fraction(ji) .GT. 1.0_r_std ) THEN
!!$             evap_fraction(ji) =  1.0_r_std
!!$          ENDIF
!!$
!!$
!!$          ! set a lowerbound/miinum of evapofraction, as we tune this parameter 
!!$          ! IF (evap_fraction(ji) .LT.  0.01_r_std) THEN
!!$          !     evap_fraction(ji) = 0.01_r_std
!!$          ! ENDIF
!!$
!!$          ! We still need to add some parmeters from the  phenology for the
!!$          ! decedious forest type.
!!$          ! IF the is decidious tree.....
!!$          ! evapo_fraction(ji) = f(day start to growth, vb3, vb4)
!!$          !
!!$
!!$          IF (printlev_loc>=4) THEN
!!$
!!$             WRITE (numout, '(A40,F10.5)') 'veget_fraction:', veget(test_grid, test_pft) 
!!$             WRITE (numout, '(A40,F10.5)') 'beta3:', vbeta3(test_grid, test_pft)
!!$             WRITE (numout, '(A40,F10.5)') 'beta4:', vbeta4(ji)
!!$             WRITE (numout, '(A40,F10.5)') 'ks_weight_1:', ks_weight_1(ji)
!!$             WRITE (numout, '(A40,F10.5)') 'ks_weight_2:', ks_weight_2(ji)
!!$             WRITE (numout, '(A40,F10.5)') 'ks_tune:', ks_tune
!!$             WRITE (numout, '(A40,F10.5)') 'forest_floor_evapo_fraction:', evap_fraction(ji)
!!$             !+++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$
!!$          END IF ! (printlev_loc>=4) THEN
!!$


          !! Calculation of the last two components of the surface fluxes:

          !! TO DO: Check why there is a distinction here.

          IF (printlev_loc>=4) THEN
             WRITE(numout,*) 'ks_tune', ks_tune
             WRITE(numout,*) 'qcdrag', q_cdrag(1)
             WRITE(numout,*) 'u_speed', u_speed(1,1)
             WRITE(numout,*) 'omega 1', big_omega_1
             WRITE(numout,*) 'omega 2', big_omega_2
             WRITE(numout,*) 'omega 3', big_omega_3
             WRITE(numout,*) 'omega 4', big_omega_4
             WRITE(numout,*) 'omega 5', big_omega_5
             WRITE(numout,*) 'omega 6', big_omega_6
             WRITE(numout,*) 'omega 7', big_omega_7
             WRITE(numout,*) 'omega 8', big_omega_8
          ENDIF

          IF (.NOT. hack_can_structure) THEN


             xi_3 =  ( big_omega_6 + (big_omega_8/big_omega_1) * big_omega_2 ) /         &
                  & ( big_omega_5 - (big_omega_8/big_omega_1) * big_omega_4 )  +        &
                  & ( big_omega_6_sub + (big_omega_8_sub/big_omega_1) * big_omega_2 ) / &
                  & ( big_omega_5_sub - (big_omega_8_sub/big_omega_1) * big_omega_4 ) 

             xi_4 =  ( big_omega_7 + (big_omega_8/big_omega_1) * big_omega_3 ) /         &
                  &  ( big_omega_5 - (big_omega_8/big_omega_1) * big_omega_4 ) +        &
                  & ( big_omega_7_sub + (big_omega_8_sub/big_omega_1) * big_omega_3 ) / &
                  &  ( big_omega_5_sub - (big_omega_8_sub/big_omega_1) * big_omega_4 ) 

          ELSE

             xi_3 = beta_diffuco_evap(ji)  *                                             &
                  & ( big_omega_6 + (big_omega_8/big_omega_1) * big_omega_2 ) /         &
                  & ( big_omega_5 - (big_omega_8/big_omega_1) * big_omega_4 )  +        &
                  & beta_diffuco_sub(ji) *                                              &
                  & ( big_omega_6_sub + (big_omega_8_sub/big_omega_1) * big_omega_2 ) / &
                  & ( big_omega_5_sub - (big_omega_8_sub/big_omega_1) * big_omega_4 ) 

             xi_4 = beta_diffuco_evap(ji) *                                              &
                  & ( big_omega_7 + (big_omega_8/big_omega_1) * big_omega_3 ) /         &
                  &  ( big_omega_5 - (big_omega_8/big_omega_1) * big_omega_4 ) +        &
                  & beta_diffuco_sub(ji) *                                              &
                  & ( big_omega_7_sub + (big_omega_8_sub/big_omega_1) * big_omega_3 ) / &
                  &  ( big_omega_5_sub - (big_omega_8_sub/big_omega_1) * big_omega_4 )     

          END IF ! (hack_can_structure)


          IF (printlev_loc>=4) THEN
             WRITE (numout, *) 'dt_sechiba:', dt_sechiba
             WRITE (numout, *) 'soilcap(ji):', soilcap(1)
             WRITE (numout, *) 'xi_1:', xi_1
             WRITE (numout, *) 'xi_2:', xi_2
             WRITE (numout, *) 'xi_3:', xi_3
             WRITE (numout, *) 'xi_4:', xi_4
             WRITE (numout, *) 'evap_fraction', evap_fraction(1)
             WRITE (numout, *) 'turb_factor', turb_factor
          ENDIF

          !! Debug:

          IF (printlev_loc>=4) THEN

             WRITE (numout, *) 'prior to temp_surf_next'
             WRITE (numout, *) 'psold(ji)/cp_air:', psold(ji)/cp_air
             WRITE (numout, *) 'dt_sechiba:', dt_sechiba
             WRITE (numout, *) '(dt_sechiba/soilcap(ji)):', (dt_sechiba/soilcap(ji))
             WRITE (numout, *) 'xi_1:', xi_1
             WRITE (numout, *) 'xi_2:', xi_2
             WRITE (numout, *) 'xi_3:', xi_3
             WRITE (numout, *) 'xi_4:', xi_4

             WRITE (numout, *) 'eta_1_surf*lwdown(ji):', eta_1_surf*lwdown(ji)
             WRITE (numout, *) 'eta_2_surf:', eta_2_surf
             WRITE (numout, *) 'eta_3_surf:', eta_3_surf
             WRITE (numout, *) 'eta_4_surf:', eta_4_surf
             WRITE (numout, *) 'swnet(ji):', swnet(ji)
             WRITE (numout, *) 'swdown(ji):', swdown(ji)
             WRITE (numout, *) 'eta_4_surf*swnet(ji):', eta_4_surf*swdown(ji)

             WRITE (numout, *) 'soilflx(ji):', soilflx(ji)

             WRITE (numout, *) 'x_1(1), x_2(1), x_3(1), x_4(1), x_5(1):', x_1(1), x_2(1), x_3(1), x_4(1), x_5(1)
             WRITE (numout, *) 'y_1(1), y_2(1), y_3(1), y_4(1), y_5(1):', y_1(1), y_2(1), y_3(1), y_4(1), y_5(1)

             WRITE (numout, *) 'a_t(1), b_t(1), c_t(1), d_t(1):', a_t(1), b_t(1), c_t(1), d_t(1)
             WRITE (numout, *) 'a_q(1), b_q(1), c_q(1), d_q(1):', a_q(1), b_q(1), c_q(1), d_q(1)

             WRITE (numout, *) 'a_t(jnlvls), b_t(jnlvls), c_t(jnlvls), d_t(jnlvls):', a_t(1), b_t(1), c_t(1), d_t(1)
             WRITE (numout, *) 'a_q(jnlvls), b_q(jnlvls), c_q(jnlvls), d_q(jnlvls):', a_q(1), b_q(1), c_q(1), d_q(1)

             WRITE (numout, *) 'big_r_prime(1): ', big_r_prime(1)

             WRITE (numout, *) 'big_omega_1, big_omega_2: ', big_omega_1, big_omega_2
             WRITE (numout, *) 'big_omega_3, big_omega_4: ', big_omega_3, big_omega_4
             WRITE (numout, *) 'big_omega_5, big_omega_6: ', big_omega_5, big_omega_6
             WRITE (numout, *) 'big_omega_7, big_omega_8: ', big_omega_7, big_omega_8

             WRITE (numout, *) 'big_e(1), big_f(1), big_g(1): ', big_e(1), big_f(1), big_g(1)

          END IF ! (printlev_loc>=4)




          !! Calculation of the surface temperature and other surface variables:

          temp_surf_next = ( psold(ji)/cp_air  + (dt_sechiba/(soilcap(ji))) * ( xi_1 + xi_3   &
               & + eta_1_surf*lwdown(ji)                                  &
               & + eta_4_surf*swdown(ji)                                  &
               & + eta_3_surf + soilflx(ji)) ) /  &
               & ( 1.0d0 - (xi_2 + xi_4 + eta_2_surf) * (dt_sechiba/(soilcap(ji))))

          IF (printlev_loc>=4) THEN
             WRITE(numout,*) 'psold', psold(1)
             WRITE(numout,*) 'soilcap', soilcap(1)
             WRITE(numout,*) 'xi_1', xi_1
             WRITE(numout,*) 'xi_2', xi_2
             WRITE(numout,*) 'xi_3', xi_3
             WRITE(numout,*) 'xi_4', xi_4
             WRITE(numout,*) 'eta_1_surf', eta_1_surf
             WRITE(numout,*) 'eta_2_surf', eta_2_surf
             WRITE(numout,*) 'eta_3_surf', eta_3_surf
             WRITE(numout,*) 'eta_4_surf', eta_4_surf
             WRITE(numout,*) 'lw_down', lwdown(1)
             WRITE(numout,*) 'sw_down', swdown(1)
             WRITE(numout,*) 'soilflw', soilflx(1)
             WRITE(numout,*) 'temp_surf_next', temp_surf_next
          ENDIF


          temp_sol_new(ji) = temp_surf_next


          psnew(ji) = temp_surf_next * cp_air

          qsol_sat_new(ji) = qsol_sat(ji) + &
               &  (1.0d0/(cp_air)) * pdqsold(ji) * &
               &  (psnew(ji) - psold(ji))


          delta_temp = temp_surf_next - psold(ji)/cp_air

          epot_air_new(ji) = (1.0d0/( rau(ji) * speed * q_cdrag(ji))) *                 &
               & (cp_air * b_t(1) - psold(ji) - cp_air * delta_temp) /          &
               & ((1.0d0/( rau(ji) * speed * q_cdrag(ji))) - a_t(1) ) +         &
               & cp_air * temp_surf_next



          IF (printlev_loc>=4) THEN
             WRITE(numout,*) 'qsol_sat_new', qsol_sat_new
          ENDIF

          !! Debug:

          IF (printlev_loc>=4) THEN
             WRITE (numout,*) 'temp_surf_next:', temp_surf_next
             WRITE (numout,*) 'psold and psnew:', psold,psnew
             WRITE (numout,*) 'q_cdrag:', q_cdrag
             WRITE (numout,*) 'speed:', speed
             WRITE (numout,*) 'delta_temp:', delta_temp
             WRITE (numout,*) 'rau:', rau
             WRITE (numout,*) 'cp_air:', cp_air
          ENDIF ! (printlev_loc=>4)   
          ! section modified from enerbil_flux

          !! \latexonly 
          !!     \input{enerbilsurftemp17.tex}
          !! \endlatexonly




          fevap = (chalev0 * beta_diffuco_evap(ji) +                                         &
               & chalsu0 * beta_diffuco_sub(ji)) * ((b_q(1) - qsol_sat(ji)) /             &
               & ((1.0d0/(rau(ji)*speed*q_cdrag(ji))) - a_q(1) ) )                        &
               & - (chalev0 * beta_diffuco_evap(ji) + chalsu0 * beta_diffuco_sub(ji))     &
               & * ((pdqsold(ji)/cp_air) /                                                &
               & ((1.0d0/(rau(ji)*speed*q_cdrag(ji))) - a_q(1) )  )  *                    &
               & (cp_air*temp_surf_next - psold(ji))




          IF (ABS(fevap) < EPSILON(un)) THEN            
             !! \latexonly 
             !!     \input{enerbilsurftemp18.tex}
             !! \endlatexonly

             qair_new(ji) = qair(ji)                                                                   
          ELSE
             !! \latexonly 
             !!     \input{enerbilsurftemp19.tex}
             !! \endlatexonly         
             qair_new(ji) = (1.0d0/( rau(ji) * speed * q_cdrag(ji))) * un / &
                  & ( chalsu0 *  vbeta1(ji) * (un - vbeta5(ji)) + &
                  & chalev0 * ((un - vbeta1(ji))*(un - vbeta5(ji)) * vbeta(ji) * 1.0d0 + vbeta5(ji)) ) &
                  & * fevap + qsol_sat_new(ji)  
          ENDIF ! (ABS(fevap) < EPSILON(un))


          !! Debug:

          IF(printlev_loc >= 4)THEN
             WRITE(numout,*)'Why is qair_new the same as qsol_sat_new?'
             WRITE(numout,*)'first term', 1.0d0/( rau(ji) * speed * q_cdrag(ji))  
             WRITE(numout,*)'vbeta',vbeta(ji)
             WRITE(numout,*)'vbeta1',vbeta1(ji)
             WRITE(numout,*)'vbeta5',vbeta5(ji)
             WRITE(numout,*)'chalev0',chalev0
             WRITE(numout,*)'chalsu0',chalsu0
             WRITE(numout,*)'fevap',fevap
             WRITE(numout,*)'qair_new(ji)',qair_new(ji)
          ENDIF




          !! In case of water stress and previous hydraulic architecture, the solution has
          !! to converge toward a precise value of transpiration. Consequently, vbeta3 is
          !! recalculated.

          IF (hydrol_flag) THEN  ! overall hydrol flag
             ! vbeta2sum(:) = zero
             IF (hydrol_flag2(ji)) THEN
                DO jv = 1, nvm
                   IF (.NOT. hydrol_flag3(ji, jv)) THEN 
                      ! this is for cases when the supply term is higher than the 
                      !     transpiration) - i.e. the normal state of affairs
                      ! WRITE(numout,*) '201013 vbeta3(ji,jv) is unchanged as: ', vbeta3(ji,jv)
                   ELSE
                      ! WRITE(numout,*) '261113 vbeta3(ji,jv) is originally: ', vbeta3(ji,jv)
                      IF(printlev_loc >= 4) THEN
                         WRITE(numout,*)'dt_sechiba',dt_sechiba
                         WRITE(numout,*)'vbeta1',vbeta1(ji)
                         WRITE(numout,*)'qsol_sat_new',qsol_sat_new(ji)
                         WRITE(numout,*)'qair_new',qair_new(ji)
                         WRITE(numout,*)'rau',rau(ji)
                         WRITE(numout,*)'speed',speed
                         WRITE(numout,*)'q_cdrag',q_cdrag(ji)
                      ENDIF
                      !+++CHECK+++
                      ! The if-statement was added to prevent crashes
                      ! Check whether it makes sense that vbeta3 = 0 when 
                      ! vbeta1 = 1
                      IF ( (vbeta1(ji) .NE. 1.0d0) .AND. &
                           & (qsol_sat_new(ji) .NE. qair_new(ji)) ) THEN
                         vbeta3(ji,jv) = ( veget_max(ji,jv)* E_supply(ji,jv) ) / &
                              & ( dt_sechiba * (1.0d0 - vbeta1(ji)) * (qsol_sat_new(ji) - qair_new(ji)) * &
                              & ( rau(ji) * speed * q_cdrag(ji)) )
                      ELSE
                         vbeta3(ji,jv) = zero
                      ENDIF
                      !+++++++++++
                      ! WRITE(numout,*) '261113 new vbeta3(ji,jv) is: ', vbeta3(ji,jv)    
                   END IF ! (hydrol_flag3(ji, jv) .eq. .FALSE.)
                ENDDO ! jv = 1, nvm
             ELSE
                ! do nothing
             END IF ! (hydrol_flag2(ji) .eq. .TRUE.)        
          END IF ! IF (hydrol_flag .eq. .TRUE.)       




          !! If hydrol_flag==TRUE, the solution has not converged yet so we do not exit the loop

          IF (hydrol_flag) THEN
             ! We aim to calculate a new delta_temp using the new vbeta3
          ELSE
             ! in non-hydrol interface case (the first call to the energy budget in the sechiba cycle,
             !     we can exit this loop straightaway)
             EXIT convergence_loop 
          END IF


          !! Values for the convergence loop

          !++++ copy the current step in convergence loop for checking the error ++++++ 
          temp_sol_dummy2 = temp_surf_next
          !+++++++++++++++++++++++

          temp_sol_err = abs(temp_sol_dummy2 - temp_sol_dummy1) 
          ! print out the absolute error
!!$                 WRITE(numout, '(A40,ES18.6)') "120615 Absolute Error:", temp_sol_err

          !+++++++ set a critiria then we exit the loop +++++
          !IF ( temp_sol_error .LT. epsilon_err ) THEN
          !    EXIT convergence_loop 
          !ENDIF    
          !++++++++++++++++++++++++++++++++++++++++++++++++++



       END DO convergence_loop


       ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  end convergence loop here  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



       !! The surface temperature is now calculated, all the coefficients in all layers and the 
       !! surface temperature are known. We can calculate all the temperature in all layers 
       !! starting from the bottom toward the top layer.

       !! Calculation of the latent and sensible heat fluxes at the surface:

       flux_ground_h =   (xi_1 + (xi_2 * temp_surf_next))
       flux_ground_le =   (xi_3 + (xi_4 * temp_surf_next))  




       !! Calculation of all the temperatures and humidities (air+leaves): 


       t_l_next(1) = big_e(1) * (flux_ground_le) + big_f(1) * (flux_ground_h) + big_g(1)

       t_a_next(1) = (a_t(1) * (flux_ground_h)) + b_t(1) + &
            & (c_t(1) * t_l_next(1)) + (d_t(1) * (flux_ground_le)) 

       q_a_next(1) = (a_q(1) * (flux_ground_le)) + b_q(1) + &
            & (c_q(1) * t_l_next(1)) + (d_q(1) * (flux_ground_h)) 

       IF (jnlvls .gt. 1) THEN    
          DO i = 2, jnlvls
             t_l_next(i) = (big_e(i) * q_a_next(i-1)) + &
                  & (big_f(i) * t_a_next(i-1)) + big_g(i)

             t_a_next(i) = (a_t(i) * t_a_next(i-1)) + b_t(i) + &
                  & (c_t(i) * t_l_next(i)) + (d_t(i) * q_a_next(i-1))

             q_a_next(i) = (a_q(i) * q_a_next(i-1)) + b_q(i) + &
                  & (c_q(i) * t_l_next(i)) + (d_q(i) * t_a_next(i-1))
          END DO ! i = 1, jnlvls
       END IF ! (jnlvls .gt. 1) THEN

       t_a_next_grid(ji, :) = t_a_next(:)

       q_a_next_grid(ji, :) = q_a_next(:)

       temp_leaf_next_grid(ji, :) = t_l_next(:)

       column_qsatt_leaf(:) = t_l_next(:) * jalpha(:) + jbeta(:)


       delta_temp_canopy(ji,:) = t_a_next_grid(ji,jnlvls_canopy) - t_a_next_grid(ji,1)

       !! Calculation of the latent and sensible heat fluxes in all layers:

       !! The lowest level flux is always the flux from the ground
       !! NB: Sign convention: positive flux is the flow of heat from high temperature to low
       !! temperature or from high humidity to low humidity respectably

       flux_h(1) = flux_ground_h
       flux_le(1) = flux_ground_le


       IF (jnlvls .gt. 1) THEN
          DO i = 2, jnlvls
             flux_h(i) = - cp_air * rau(ji) * &
                  & (t_a_next(i) - t_a_next(i-1)) * k_eddy(i-1)/delta_z(i-1)
             flux_le(i) = - chalev0 * rau(ji) * &
                  & (q_a_next(i) - q_a_next(i-1)) * k_eddy(i-1)/delta_z(i-1)
          END DO ! i = 2, jnlvls
       END IF !  (jnlvls .eq. 1) THEN   


       flux_h_grid(ji,  :) = flux_h(:)
       flux_le_grid(ji, :) = flux_le(:)
       ! here, we save the wind profile for each grid.  


       IF (printlev_loc>=4) THEN
          WRITE(numout,*) 'flux_h', flux_h(:)
          WRITE(numout,*) 'flux_le', flux_le(:)
          WRITE(numout,*) 't_a_next', t_a_next(:)
          WRITE(numout,*) 'q_a_next', q_a_next(:)
       ENDIF

       !! TO DO: To remove

       IF (.NOT. hack_can_structure) THEN 

          DO i = 11, 20 ! 090215 calculate provisional column values of q_surf_leaf
             !  column_q_surf_leaf(i) = vbeta1(ji) * ( ( 1.0d0 - vbeta5(ji)) + vbeta5(ji) ) * column_qsatt_leaf(i) + &
             !                         &  (1.0d0 - vbeta1(ji)) * (1.0d0 - vbeta5(ji)) * (vbeta4(ji) + vbeta2sum(ji) + &
             !                         &  vbeta3sum(ji)) * column_qsatt_leaf(i)
             !aslacolumn_q_surf_leaf(i) = profile_vbeta3(ji,jv,(i-10)) * column_qsatt_leaf(i)
             !aslacolumn_q_surf_leaf(i) = MAX(column_q_surf_leaf(i), q_a_next(i))
          END DO ! i = 1, jnlvls

       END IF ! (hack_can_structure) THEN 



       !! Calculation of the net shortwave & longwave radiations in each canopy layer:
       !! Short-wave radiations term can be calculated from eta4 
       !! Long-wave radiations tern can be calculated from eta1 to eta3   

       DO i =1, jnlvls

          IF (i .EQ. 1) THEN
             !! IF (i .eq.1), the surface layer net radiation term should be equals 
             !! (LE + H + G) + (h_horizontal + le_horizontal)    
             !! Note the sign of eta_3 is not the same for surface and leave layer balance calculation
             !! Here, the longwave and shortwave were separated according to eta1, eta2, eta3 and eta4 

             flux_netsw(i) = eta_4_surf*swdown(ji) 
             flux_netsw(i) = flux_netsw(i) + (eta_4(i) * swdown(ji))

             flux_netlw(i) = eta_1_surf*lwdown(ji) + (eta_2_surf * temp_surf_next) - eta_3_surf 
             flux_netlw(i) = flux_netlw(i) + (eta_1(i) * lwdown(ji)) + (eta_2(i)*t_l_next(i)) +  eta_3(i) 

          ELSE IF (i .EQ. jnlvls) THEN

             flux_netsw(i) =  swnet(ji) 

             flux_netlw(i) =  lwdown(ji) - lwup(ji) 

          ELSE

             !! For other layers, the net rdiation can be calculated as the summation of leave layer 
             !! energy balance (hrizontal terms) and previous layer net radiation term (vertical term).  
             !! Here, the longwave and shortwave were separated according to eta1, eta2, eta3 and eta4 

             flux_netsw(i) =  (eta_4(i) * swdown(ji)) + flux_netsw(i-1)                         
             flux_netlw(i) =  (eta_1(i) * lwdown(ji)) + (eta_2(i)*t_l_next(i)) +  eta_3(i)  +  flux_netlw(i-1) 

          ENDIF

          flux_netrad(i) = flux_netsw(i) + flux_netlw(i) - soilflx(ji)                                 

       ENDDO

       !! TO DO: Already written, to remove

       flux_h_grid(ji, :)  = flux_h(:)
       flux_le_grid(ji, :) = flux_le(:)


       !! The fluxes at the top of the column are the fluxes of the gridcell:

       fluxsens(ji) = flux_h(jnlvls)
       fluxlat(ji) = flux_le(jnlvls)  

       
       !! TO DO: Check if it is still necessary:

!!$       IF (.NOT. hack_can_structure) THEN 
!!$          DO i = 1, jnlvls ! 090215 calculate relhum of air
!!$             column_air_relhum(i) = (  (column_q_surf_leaf(i) * pb(ji)) / (tetens_1 + column_q_surf_leaf(i) * tetens_2)   )  / &
!!$                  & ( (column_qsatt_leaf(i) * pb(ji)) / (tetens_1 + column_qsatt_leaf(i) * tetens_2) )
!!$          END DO ! i = 1, jnlvls
!!$       END IF ! (hack_can_structure) THEN 


       !! Debug:

       IF (printlev_loc>=4) THEN

          WRITE (numout,*) 'balance test delta temp: ', temp_surf_next - (psold(ji)/cp_air)

          WRITE (numout,*) 'balance test fluxes sum: ', (-flux_h(1) - flux_le(1) + &
               & eta_1_surf*lwdown(ji) + eta_4_surf*swdown(ji) + eta_3_surf + soilflx(ji)   &
               & + eta_2_surf * temp_surf_next) &
               & * (dt_sechiba/soilcap(ji))

          WRITE (numout,*) 'balance test delta temp / (dt_sechiba/soilcap(ji): ', &
               & (temp_surf_next - (psold(ji)/cp_air)) / (dt_sechiba/soilcap(ji))

          WRITE (numout,*) 'balance test fluxes sum: ', -flux_h(1) - flux_le(1) + &
               & eta_1_surf*lwdown(ji) + eta_4_surf*swdown(ji) + soilflx(ji)   &
               & + (eta_2_surf * temp_surf_next) - eta_3_surf

          WRITE (numout,*) 'total balance test: ', &
               & (temp_surf_next - (psold(ji)/cp_air)) / (dt_sechiba/soilcap(ji)) - &
               & (-flux_h(1) - flux_le(1) + &
               & eta_1_surf*lwdown(ji) + eta_4_surf*swdown(ji) + soilflx(ji)   &
               & + (eta_2_surf * temp_surf_next) - eta_3_surf)

          WRITE (numout,*) 'balance test flux_h(1): ', flux_h(1)
          WRITE (numout,*) 'balance test flux_le(1): ', flux_le(1)
          WRITE (numout,*) 'balance test eta_1_surf*lwdown(ji): ', eta_1_surf*lwdown(ji)
          WRITE (numout,*) 'balance test eta_1_surf', eta_1_surf
          WRITE (numout,*) 'balance test lwdown', lwdown(ji)
          WRITE (numout,*) 'balance test eta_2_surf * temp_surf_next: ', eta_2_surf * temp_surf_next
          WRITE (numout,*) 'balance test eta_3_surf: ', eta_3_surf
          WRITE (numout,*) 'balance test eta_4_surf*swdown(ji): ', eta_4_surf*swdown(ji)
          WRITE (numout,*) 'balance test soilflx(ji): ', soilflx(ji)

       END IF !(printlev_loc>=4) 


       !! Calculation of the horizontal fluxes in each layer 
       !! (fluxes between the leaves and the atmosphere):

       source_sink_t(1) = cp_air * rau(ji) * &
            & ((t_a_next(1) - t_a_pres(1)) * delta_h(1)) / dt_sechiba

       source_sink_q(1) = chalev0 * rau(ji) * &
            & ((q_a_next(1) - q_a_pres(1)) * delta_h(1)) / dt_sechiba

       source_sink_t_l(1) = ((t_l_next(1) - t_l_pres(1)) * delta_h(1)) / dt_sechiba

       horiz_flux_t(1) = (t_l_next(1) - t_a_next(1))/big_r(1) * cp_air * rau(ji) 

       horiz_flux_q(1) = (((jalpha(1) * t_l_next(1) + jbeta(1)) - q_a_next(1))/ &
            & big_r_prime(1)) * chalev0 * rau(ji)

       IF (jnlvls .gt. 1) THEN
          DO i = 2, jnlvls

             source_sink_t(i) = cp_air * rau(ji) * &
                  & ((t_a_next(i) - t_a_pres(i)) * delta_h(i)) / dt_sechiba
             source_sink_q(i) = chalev0 * rau(ji) * &
                  & ((q_a_next(i) - q_a_pres(i)) * delta_h(i)) / dt_sechiba

             source_sink_t_l(i) = ((t_l_next(i) - t_l_pres(i)) * delta_h(i)) / dt_sechiba

             horiz_flux_t(i) = (t_l_next(i) - t_a_next(i))/big_r(i) * cp_air * rau(ji) 


             !! TO DO: Chack and remove what is commented:

             !    IF (.NOT. hack_can_structure) THEN
             !        horiz_flux_q(i) = phi_trans(i)
             !    ELSE
             horiz_flux_q(i) = (((jalpha(i) * t_l_next(i) + jbeta(i)) - q_a_next(i))/ &
                  & big_r_prime(i)) * chalev0 * rau(ji)
             !    END IF 

          END DO ! i = 1, jnlvls
       END IF ! (jnlvls .gt. 1) THEN



       !! Debug ID3.
       !! Transpir was not equal to the sum of the transpiration in all the layers
       transpir(ji,:) = 0
       DO jv=1,nvm
          transpir(ji,jv) = SUM(horiz_flux_q(:))/chalev0*dt_sechiba
       ENDDO
       !hydrol_flag = .FALSE.

       !! End ID3


       !! Check if the energy balance is closed in each layer:

       IF (jnlvls .gt. 1) THEN  

          delta_t_a(1) = cp_air * rau(ji) * ((t_a_next(1) - t_a_pres(1)) &
               & * delta_h(1)) / dt_sechiba
          sum_flux_a(1) = - flux_h(1+1) + ( - flux_h(1) + horiz_flux_t(1) )

          delta_q_a(1) = chalev0 * rau(ji) * &
               & ((q_a_next(1) - q_a_pres(1)) * delta_h(1)) / &
               & dt_sechiba

          sum_flux_q_a(1) = - flux_le(1+1) + ( - flux_le(1) + horiz_flux_q(1) )

          DO i = 2, jnlvls-1
             delta_t_a(i) = cp_air * rau(ji) * ((t_a_next(i) - t_a_pres(i)) &
                  & * delta_h(i)) / dt_sechiba
             sum_flux_a(i) = - flux_h(i+1) + ( flux_h(i) + horiz_flux_t(i) )
          END DO ! i = 1, jnlvls-1

          DO i = 2, jnlvls-1
             delta_q_a(i) = chalev0 * rau(ji) * &
                  & ((q_a_next(i) - q_a_pres(i)) * delta_h(i)) / &
                  & dt_sechiba
             sum_flux_q_a(i) = - flux_le(i+1) + ( flux_le(i) + horiz_flux_q(i) )
          END DO ! i = 1, jnlvls-1


          DO i = 1, jnlvls-1
             delta_temp_leaf(i) = rho_veg*leaf_tks * pad_deltah(i) * jtheta(i) * &
                  & ((t_l_next(i) - t_l_pres(i)) * delta_h(i)) / &
                  & dt_sechiba
             sum_flux_leaf(i) = (eta_1(i) * lwdown(ji)) + (eta_2(i)*t_l_next(i)) + &
                  & eta_3(i) + (eta_4(i) * swdown(ji)) - &
                  & ( horiz_flux_t(i) + horiz_flux_q(i))
          END DO ! i = 1, jnlvls-1

          DO i = 1, jnlvls-1
             IF (ABS(source_sink_t(i) - sum_flux_a(i)) .gt. 1.0d-5) THEN
                WRITE (numout,*) 'level: ', i, ' does not balance at atmos temperature'
             ELSE
                ! do nothing
             END IF

             IF (ABS(source_sink_q(i) - sum_flux_q_a(i)) .gt. 1.0d-5) THEN
                WRITE (numout,*) 'level: ', i, ' does not balance at atmos specific humidity'
             ELSE
                ! do nothing
             END IF

             IF (ABS(delta_temp_leaf(i) - sum_flux_leaf(i)) .gt. 1.0d-5) THEN
                WRITE (numout,*) 'level: ', i, ' does not balance at leaf temperature'
             ELSE
                ! do nothing
             END IF

          END DO ! i = 1, jnlvls-1

!!$          TO DO : To check

!!$          DO i = 1, jnlvls-1
!!$             source_sink_sum_t = (source_sink_sum_t + source_sink_t(i) + horiz_flux_t(i))
!!$             source_sink_sum_q = (source_sink_sum_q + source_sink_q(i) + horiz_flux_q(i))
!!$          END DO ! i = 1, jnlvls-1

       END IF ! (jnlvls .gt. 1) THEN


       !! End of the calculation

    END DO gridloop ! ji=1,kjpindex 

    CALL xios_orchidee_send_field("delta_temp_canopy",delta_temp_canopy)

  END SUBROUTINE mleb_calc_column





  !!  ===========================================================================================================================
  !! SUBROUTINE: mleb_alpha_beta_coeff 
  !!
  !!
  !>\BRIEF: calculates the rate of change of saturated surface humidity of the leaf surface 
  !!
  !! DESCRIPTION: Calculates the alpha and beta coefficients for the specific
  !! leaf humidity (see Ryder et al., 2016, section 3.1). 
  !! 
  !!
  !! RECENT CHANGE(S): None
  !!
  !! MAIN OUTPUT VARIABLE(S) : jalpha, jbeta, jalpha_surf, jbeta_surf, q_sat_surf, temp_surf_pres
  !!
  !! REFERENCE(S): Monteith & Unsworth (2008 edition)
  !!
  !! FLOWCHART: 
  !_ ==============================================================================================================================


  SUBROUTINE mleb_alpha_beta_coeff (temp_atmos_pres_grid, temp_leaf_pres_grid, jalpha, jbeta, &
       & jalpha_surf, jbeta_surf, temp_surf_pres, q_sat_surf, kjpindex, &
       & pb, temp_sol, jpdqsold)

    ! we use the approach outlined in Monteith & Unsworth (2008 edition)

    IMPLICIT NONE

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                               :: kjpindex                       !! Domain size (-)
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)           :: pb                             !! Lowest level pressure (hPa)
    REAL(r_std), DIMENSION (kjpindex, jnlvls), INTENT (in)   :: temp_atmos_pres_grid           !! atmospheric temperature down the column (K)
    REAL(r_std), DIMENSION (kjpindex, jnlvls), INTENT (in)   :: temp_leaf_pres_grid            !! leaf temperature down the column (K)
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: temp_sol                       !! Surface temperature (K)

    !! 0.2 Output variables

    REAL(r_std), DIMENSION (jnlvls), INTENT (out)            :: jalpha, jbeta
    REAL(r_std), DIMENSION (kjpindex), INTENT (out)          :: jalpha_surf, jbeta_surf
    REAL(r_std), DIMENSION (kjpindex), INTENT (out)          :: q_sat_surf
    REAL(r_std), DIMENSION (kjpindex), INTENT(out)           :: jpdqsold

    !! 0.3 Modified variables
    REAL(r_std), DIMENSION (kjpindex),INTENT (inout)         :: temp_surf_pres                 !! surface temperature at PRESENT time step (Kelvins)

    !! 0.4 Local variables

    REAL(r_std), DIMENSION (jnlvls)                          :: temp_atmos                     !! atmospheric temperature down the column (K)
    REAL(r_std), DIMENSION (jnlvls)                          :: temp_leaf                      !! leaf temperature down the column (K)

    INTEGER(i_std)                                           :: ji                             !! grid point loop
    INTEGER(i_std)                                           :: i                              !! levels loop

    REAL(r_std)                                              :: lambda_const = 2480.0d0        !! (J g^-1) at 10degC n.b. this changes by 2.4 per degree
    REAL(r_std)                                              :: bigm_w = 18.0d0                !! (g mol^-1) approximate molar mass of water
    REAL(r_std)                                              :: bigr_const = 8.314d0           !! (J mol^-1 K^-1) the molar gas constant
    REAL(r_std)                                              :: eta_const = 0.622d0            !! ratio of molecular mass of water to molecular mass of air
    REAL(r_std)                                              :: biga_const = 17.27d0
    REAL(r_std)                                              :: tstar_const = 273.0d0          !! (Kelvins)
    REAL(r_std)                                              :: tprime_const = 36.0d0          !! (Kelvins)
    REAL(r_std)                                              :: e_sat_tstar = 0.611d0          !! (kPa)
    REAL(r_std)                                              :: zrapp,zcorr,ztemperature,zqsat !! Temporary vector variables      
    REAL(r_std), DIMENSION (jnlvls)                          :: q_sat                          !! (kg/kg)
    REAL(r_std), DIMENSION (jnlvls)                          :: big_delta                      !! delta q_sat / delta temp at leaf temperature
    REAL(r_std), DIMENSION (jnlvls)                          :: leaf_qsol_sat                    !! Soil temperature (K)
    REAL(r_std), DIMENSION (jnlvls)                          :: leaf_dev_qsol                    !! Soil temperature (K)
    REAL(r_std), DIMENSION (jnlvls)                          :: leaf_pdqsold                     !! Soil temperature (K)
    REAL(r_std), DIMENSION (kjpindex)                        :: big_delta_surf
    REAL(r_std), DIMENSION (kjpindex)                        :: j_qsol_sat, j_dev_qsol
    REAL(r_std), DIMENSION (kjpindex)                        :: j_dev_q, j_q_sat


    !   =from qsat_moisture=================================================================================================


    !! This routine calculates the alpha and beta coefficients needed for the calculation of the saturated humidity.

    !! TO DO: Clean a bit the calculations.

    !! 1. Initialisation
    zrapp = msmlr_h2o/msmlr_air
    zcorr = 0.00320991_r_std

    jalpha(:) = 0.0d0
    jbeta(:) = 0.0d0
    jalpha_surf(:) = 0.0d0 
    jbeta_surf(:) = 0.0d0

    q_sat(:) = 0.0d0
    big_delta(:) = 0.0d0
    leaf_qsol_sat(:) = 0.0d0
    leaf_dev_qsol(:) = 0.0d0
    leaf_pdqsold(:) = 0.0d0
    big_delta_surf(:) = 0.0d0 
    j_qsol_sat(:) = 0.0d0
    j_dev_qsol(:) = 0.0d0
    jpdqsold(:) = 0.0d0

    DO ji = 1, kjpindex

       !! 2. Computes saturated humidity one time and store in qsfrict local array

       ztemperature = temp_sol(ji)
       j_qsol_sat(:) = 0.0d0
       j_dev_qsol(:) = 0.0d0

       CALL qsatcalc (kjpindex, temp_sol, pb, j_qsol_sat)

       CALL dev_qsatcalc (kjpindex, temp_sol, pb, j_dev_qsol)

       temp_surf_pres(ji) = temp_sol(ji)

       !+++ According to r2093 in the trunk, the line below is modified to avoid possible instability in the coupling.

       ! jpdqsold(ji) = j_dev_qsol(ji) * ( pb(ji)**kappa ) / cp_air

       jpdqsold(ji) = j_dev_qsol(ji)

       q_sat(:) = 0.0d0
       big_delta(:) = 0.0d0

       ! calculate jalpha_surf and jbeta_surf first of all:

       ! Using Tetens' method for close approximation of saturated specific humidity 
       !   (e.g. 2.27 of Monteith & Unsworth (2008)) and the specific humidity approximation 
       !   in terms of vapour pressure (e.g. 2.32 of Monteith & Unsworth (2008))

       q_sat_surf(ji) = (eta_const / 100000.0d0)  * (e_sat_tstar*1000.0d0) *  &
            & EXP((biga_const*(temp_surf_pres(ji) - tstar_const)/(temp_surf_pres(ji) - tprime_const)))


       !! Debug ID5.
       !! Crash floating overflow problem. To re-test
!!$       big_delta_surf(ji) = ( lambda_const * bigm_w * q_sat_surf(ji) ) / &
!!$                      & ( bigr_const * (temp_surf_pres(ji))**2 )
       !! End ID5


       IF (printlev_loc>=4) THEN
          WRITE (numout,*) 'in calc_jalpha_jbeta, big_delta_surf is: ', big_delta_surf(ji)    
       END IF !(printlev_loc>=4)

       ! jalpha_surf(ji) = j_dev_qsol(ji)
       jalpha_surf(ji) = jpdqsold(ji) !* ( pb(ji)**kappa ) / cp_air
       jbeta_surf(ji) = j_qsol_sat(ji) - temp_sol(ji) * jpdqsold(ji) !* ( pb(ji)**kappa ) / cp_air
       q_sat_surf(ji) = j_qsol_sat(ji)

       IF (jnlvls .gt. 1) THEN

          DO i=1, jnlvls

             temp_atmos(i) = temp_atmos_pres_grid(ji,i)
             temp_leaf(i) = temp_leaf_pres_grid(ji,i)

             ! Using Tetens' method for close approximation of saturated specific humidity 
             !   (see 2.27 of Monteith & Unsworth (2008)) and the specific humidity 
             !   approximation in terms of vapour pressure (see 2.32 of Monteith & Unsworth (2008))

             q_sat(i) = (eta_const / 100000.0d0)  * (e_sat_tstar*1000.0d0) *  &
                  & EXP((biga_const*(temp_leaf(i) - tstar_const)/(temp_leaf(i) - tprime_const)))

             big_delta(i) = ( lambda_const * bigm_w * q_sat(i) ) / &
                  & ( bigr_const * temp_leaf(i)**2 )

             jalpha(i) = big_delta(i) 
             jbeta(i) = q_sat(i) - (temp_leaf(i) * jalpha(i))

             CALL qsatcalc (kjpindex, temp_leaf_pres_grid(:,i), pb, j_q_sat)

             !   CALL dev_qsatcalc (kjpindex, temp_sol, pb, dev_qsol) jtest3
             CALL dev_qsatcalc (kjpindex, temp_leaf_pres_grid(:,i), pb, j_dev_q)

             jalpha(i) = j_dev_q(ji) * ( pb(ji)**kappa ) / cp_air

             jbeta(i) = j_q_sat(ji) - temp_leaf(i) * j_dev_q(ji) * ( pb(ji)**kappa ) / cp_air

          END DO !i = 1, jnlvls

       END IF ! IF (jnlvls .gt. 1)

    END DO ! ji = 1, kjpindex

  END SUBROUTINE mleb_alpha_beta_coeff



  !!  ===========================================================================================================================
  !! SUBROUTINE                  : mleb_profile
  !!
  !!
  !>\BRIEF                       : calculates the rate of change of saturated surface humidity of the leaf surface 
  !!
  !! DESCRIPTION                 : This routine will allocate a LAD to each level of the canopy using the approach 
  !!                                           of Yang et al. (1999) and the parameters they allocate to an aspen stand (if not 
  !!                                           already known through read-in data). Based on this LAD profile, canopy turbulence 
  !!                                           statistics will be allocated using the method of Massman & Weil (1999).
  !! 
  !!
  !! RECENT CHANGE(S): None
  !!
  !! MAIN OUTPUT VARIABLE(S)     : jalpha, jbeta, jalpha_surf, jbeta_surf, q_sat_surf, temp_surf_pres
  !!
  !! REFERENCE(S)                : Massman, W.J. & Weil, J.C., 1999. An analytical one-dimensional second-order 
  !!                               closure model of turbu lence statistics and the lagrangian time scale within and 
  !!                               above plant canopies of arbitrary structure. Boundary-Layer Meteorology, 91, 
  !!                               pp.81–107.
  !!
  !!                               Raupach, M.R., Antonia R.A. and Rajagpalan, S., 1991. 'Rough-wall turbulent 
  !!                               boundary layers', Applied Mechanics Review, 44, pp. 1-25
  !!
  !!                               Yang, X., Witcosky, J.J. & Miller, D.R., 1999. Vertical overstory canopy 
  !!                               architecture of temperate deciduous in the Eastern United States. Forest Science, 
  !!                               45(3), pp.349–358.
  !!
  !!                               Wohlfahrt G. and A. Cernusca 2002, Monmentum tansfer by a mountain
  !!                               meadow canopy: a simulation analysis based on Massman's (1997) model,
  !!                               Boundary Layer Meteorology, 103, 391-407  
  !!
  !! FLOWCHART           : 
  !_ ==============================================================================================================================


  SUBROUTINE mleb_profile (pad_deltah, layer_height, delta_z, &
       & delta_h, box_height, k_eddy, &
       & u_star, speed, u_speed, big_pad, & 
       & k_eddy_slope,k_eddy_ustar, surf_lev, top_lev, swnet, jrpad, &  
       & a_1, a_2, a_3, a_4, a_5, WC02)

    IMPLICIT NONE

    !! 0 Variable and parameter description

    !! 0.1 Input variables
    REAL(r_std), DIMENSION (0:jnlvls), INTENT(in)    :: pad_deltah             !! Plant Area Density as a function of height above the ground (m^2/m^3)
    REAL(r_std), DIMENSION(jnlvls+1), INTENT(in)     :: box_height             !! height of each box edge, starting at the surface layer (m)
    REAL(r_std), INTENT(in)                          :: surf_lev
    REAL(r_std), DIMENSION(jnlvls+1), INTENT(in)     :: top_lev
    REAL(r_std), INTENT(in)                          :: swnet                  !! Net surface short-wave flux (W m^{-2})
    REAL(r_std), INTENT(in)                          :: k_eddy_slope           !! for k_eddy parameterization  
    REAL(r_std), INTENT(in)                          :: k_eddy_ustar           !! for k_eddy parameterization 
    REAL(r_std), INTENT(in)                          :: speed

    !! 0.2 Output variables
    REAL(r_std), INTENT(out)                         :: u_star                 !! friction velocity (m/s)
    REAL(r_std), INTENT(out)                         :: layer_height           !! the height of each layer in the canopy (assumed equal for now) (m)
    REAL(r_std), DIMENSION (jnlvls)                  :: cum_pad_deltah         !! Cumulative Plant Area Density multiplied by layer height (m^2/m^3)
    REAL(r_std), DIMENSION (jnlvls), INTENT(out)     :: delta_z                !! change in height (m)
    REAL(r_std), DIMENSION (:,:), INTENT(out)     :: u_speed                !! vertical wind speed (m/s)
    REAL(r_std), DIMENSION (jnlvls), INTENT(out)     :: big_pad 
    REAL(r_std), DIMENSION (jnlvls), INTENT(out)     :: k_eddy                 !! calculated eddy diffusivity coefficient (m^2/s)
    REAL(r_std), DIMENSION (jnlvls), INTENT(out)     :: delta_h                !! height of each level (m)
    REAL(r_std), DIMENSION (0:jnlvls), INTENT(out)   :: jrpad                  !! Plant Area Density, replacement for values above (as definitions between different
    !!    parts of the model was becoming confusing)
    !!    'jrpad' is the total one-sided leaf area per unit of layer volume

    !! 0.3 Modified variables

    !! 0.4 Local variables
    REAL(r_std), DIMENSION (jnlvls)                   :: cum_jrpad             !! Cumulative Plant Area Density (m^2/m^3)

    !! 0.4 Local variables
    REAL(r_std)                                       :: param_b=0.426d0
    REAL(r_std)                                       :: param_c=1.856d0
    REAL(r_std)                                       :: c_d=0.2d0
    REAL(r_std)                                       :: p_m=1.0d0
    REAL(r_std)                                       :: jlambda
    REAL(r_std)                                       :: jalpha=0.05d0
    REAL(r_std)                                       :: small_n
    INTEGER(i_std)                                    :: i
    INTEGER(i_std)                                    :: j, k, m, n
    REAL(r_std)                                       :: c_1, c_2, c_3         !! coefficients as defined in Massman & Weil (1999) p. 83
    REAL(r_std), DIMENSION (jnlvls+1)                 :: jepsilon              !! epsilon coefficient as in Massman & Weil (1999), 
    !!   equation (6) and just after
    REAL(r_std), DIMENSION (jnlvls)                   :: eta_z                 !! eta coefficient as in Massman & Weil (1999), equation (3)
    REAL(r_std), DIMENSION (jnlvls)                   :: sigma_exp_1
    REAL(r_std), DIMENSION (jnlvls)                   :: sigma_exp_2
    REAL(r_std), DIMENSION (jnlvls)                   :: sigma_exp_3
    REAL(r_std), DIMENSION (jnlvls)                   :: sigma_e_over_u_star   !! sigma_e normalised by u_star (-)
    REAL(r_std), DIMENSION (jnlvls)                   :: sigma_u_over_u_star   !! sigma_u normalised by u_star (-)
    REAL(r_std), DIMENSION (jnlvls)                   :: sigma_v_over_u_star   !! sigma_v normalised by u_star (-)
    REAL(r_std), DIMENSION (jnlvls)                   :: sigma_w_over_u_star   !! sigma_w normalised by u_star (-)
    REAL(r_std), DIMENSION (jnlvls)                   :: sigma_u               !! standard deviation in u-stream (horizontal) velocity (m/s); Massman & Weil eqn. (11)
    REAL(r_std), DIMENSION (jnlvls)                   :: sigma_v               !! standard deviation in v-stream (horizontal) velocity (m/s); Massman & Weil eqn. (11)
    REAL(r_std), DIMENSION (jnlvls)                   :: sigma_w               !! standard deviation in w-stream (vertical) velocity (m/s); Massman & Weil eqn. (11)

    REAL(r_std)                                       :: gamma_1=2.4d0
    REAL(r_std)                                       :: gamma_2=1.9d0
    REAL(r_std)                                       :: gamma_3=1.25d0
    REAL(r_std)                                       :: nu_1, nu_2, nu_3
    REAL(r_std)                                       :: bigb1
    REAL(r_std)                                       :: rp_a_1 = 1.25d0       !! Raupach et al. (1991) co-efficients
    REAL(r_std)                                       :: rp_a_0 = 0.25d0       !! Raupach et al. (1991) co-efficients
    REAL(r_std)                                       :: rp_c_0 = 0.3d0        !! Raupach et al. (1991) co-efficients 
    REAL(r_std), DIMENSION (0:jnlvls)                 :: gu_l
    REAL(r_std)                                       :: gu_ext = 0.7d0
    REAL(r_std), DIMENSION (jnlvls)                   :: middle_layer_height   !! the height in the middle of each layer in the canopy (m)
    REAL(r_std), DIMENSION (jnlvls+1)                 :: z_h                   !! height above ground at the base of each layer (m)
    REAL(r_std), DIMENSION (jnlvls+1)                 :: jeta
    REAL(r_std), DIMENSION (jnlvls)                   :: t_lagrangian          !! Lagrangian timescale (s)

    ! Parameter for dynamic drag coefficient cd_z
    ! Ref:  Wohlfahrt G. and A. Cernusca 2002, Monmentum tansfer by a mountain
    ! meadow canopy: a simulation analysis based on Massman's (1997) model,
    ! Boundary Layer Meteorology, 103, 391-407 
    ! a negtive sign was adeed on a_4 (this could be the typo in origial paper) 
    LOGICAL, INTENT(in)                               :: WC02                     ! flage for the drag_coefficient is depend on LAI profile
    REAL(r_std), DIMENSION(jnlvls)                    :: c_dz
    REAL(r_std), INTENT(in)                           :: a_1 
    REAL(r_std), INTENT(in)                           :: a_2 
    REAL(r_std), INTENT(in)                           :: a_3 
    REAL(r_std), INTENT(in)                           :: a_4 
    REAL(r_std), INTENT(in)                           :: a_5 
    !_ ==============================================================================================================================

    
    !! This routine calculates the eddy diffusivity inside each layer of the canopy following the model developed by
    !! Massman & Weil (1999).


    !! Initialisation of some variables:

    layer_height = 0.0_r_std
    middle_layer_height = 0.0_r_std
    z_h = 0.0_r_std
    jepsilon = 0.0_r_std
    cum_pad_deltah = 0.0_r_std
    small_n = 0.0_r_std
    bigb1 = 0.0_r_std
    sigma_exp_1 = 0.0_r_std
    sigma_exp_2 = 0.0_r_std
    sigma_exp_3 = 0.0_r_std

    nu_1 = (gamma_1**2.0_r_std+gamma_2**2.0_r_std+gamma_3**2.0_r_std)**(-0.5_r_std)
    nu_3 = (gamma_1**2.0_r_std+gamma_2**2.0_r_std+gamma_3**2.0_r_std)**(3.0_r_std/2.0_r_std)   
    nu_2 = nu_3/6.0_r_std-(gamma_3**2.0_r_std/(2.0_r_std*nu_1))


    !! Calculation of the central node of each layer and there width:

    !! TO DO: Clean a bit the calculations of delta_h and jrpad overthere.

    DO i=1, jnlvls
       delta_h(i) = box_height(i+1) - box_height(i)

       IF(printlev_loc>= 4) THEN
          WRITE(numout,*)'191016 in enerbil  box_height(i+1),box_height(i)',box_height(i+1),box_height(i)
          WRITE(numout,*)'191016 in enerbil, delta_h(i)',delta_h(i) 
       ENDIF
    END DO !i=1, jnlvls


    DO i=1, jnlvls-1
       delta_z(i) = (box_height(i+1) + (delta_h(i+1)/2.0_r_std)) - &
            &  (box_height(i) + (delta_h(i)/2.0_r_std))
    END DO !i=1, jnlvls

    DO i=1, jnlvls+1
       z_h(i) = box_height(i)
       jeta(i) = z_h(i) / canopy_height
       ! this is from Yang et al. (1999), equation 1
       !++++TEMP++++
       ! for single layer case without STOMATE
       !jepsilon(i) = 1.0_r_std - exp ((-1.0_r_std)*(( (1-jeta(i)) /param_b ) ** (param_c)))
       jepsilon(i) = 1.0_r_std 
       !+++++++++++++
    END DO !i=1, jnlvls+1

    DO i=1, jnlvls
       middle_layer_height(i) = (z_h(i) + z_h(i+1)) / 2.0_r_std
    END DO !i=1, jnlvls


    delta_h(1) = top_lev(1) - surf_lev

    DO i = 2, jnlvls
       delta_h(i) = top_lev(i) - top_lev(i-1)
    END DO ! i = 2, jnlvls 

    DO i = 1, jnlvls-1
       delta_z(i) = (delta_h(i) + delta_h(i+1))/2.0_r_std
    END DO ! i = 1, jnlvls-1

    delta_z(jnlvls) = box_height(jnlvls) - box_height(jnlvls-1)



    !! Coefficients c_1, c_2 and c_3 as defined in Massman & Weil (1999) p. 83

    c_1 = 0.320_r_std
    c_2 = 0.264_r_std
    c_3 = 15.1_r_std


    !! Calculation of the cumulative plant area density inside each layer:

    cum_pad_deltah(jnlvls) = pad_deltah(jnlvls)
    DO i=jnlvls-1, 1, -1
       cum_pad_deltah(i) = cum_pad_deltah(i+1) + pad_deltah(i)
    END DO !i=1, jnlvls

    DO i=1, jnlvls
       big_pad(i) = pad_deltah(i) / delta_h(i)
    END DO ! i=1, jnlvls

    ! here we define jrpad correctly - the total one-sided leaf area per unit of layer volume
    DO i=1, jnlvls
       jrpad(i) = pad_deltah(i) / delta_h(i)
    END DO !i=1, jnlvls



    !! Calculation of the Massman & Weil (1999) eddy diffusivity model: 


    !+++++ CHECK +++++++
    !Calculate the eta_z with Dynamic Drag Coeffiient (c_dz) or Static Drag Coefficient Constance(c_d)
    IF ( WC02 ) THEN
       !Dyanaic Drage Coefficient vertical
       c_dz(:) = c_d !set default value c_dz to c_d
       DO i = 1, jnlvls
          c_dz(i) =  a_1**(-1.0_r_std*cum_pad_deltah(i)/a_2) + a_3**(-1.0_r_std*cum_pad_deltah(i)/a_4) + a_5   
          !WRITE(numout, *) 'C_Deff(z):', c_dz(i), 'PAIcum:',cum_pad_deltah(i)
       ENDDO
       eta_z(1) = (( c_dz(1) * pad_deltah(1)  )  )
       DO i=2, jnlvls
          eta_z(i) = (( c_dz(i) * pad_deltah(i)  )  ) + eta_z(i-1)
          !   WRITE(numout, *) 'eta_z:', eta_z(i)
       END DO

    ELSE
       !Constant Drag Coefficient
       eta_z(1) = (( c_d * pad_deltah(1) ) / p_m )
       DO i=2, jnlvls
          eta_z(i) = (( c_d * pad_deltah(i) ) / p_m ) + eta_z(i-1)
       END DO

    ENDIF
    !+++++++++++++++++++++++

    u_star = speed * (c_1 - c_2 * exp (- c_3 * eta_z(jnlvls) ) )

    jlambda = SQRT( 7.0_r_std / (3.0_r_std * jalpha**2.0_r_std * nu_1 * nu_3 ) + & 
         & (1.0_r_std / 3.0_r_std - (gamma_3**2.0_r_std * nu_1**2.0_r_std ))/  &
         & (3.0_r_std * jalpha**2.0_r_std * nu_1 * nu_2))

    small_n = eta_z(jnlvls)/((2.0_r_std*(u_star**2.0_r_std))/(speed**2.0_r_std))


    bigb1 = -(9.0_r_std* u_star / speed)/(2.0_r_std* jalpha * nu_1 * & 
         & (9.0_r_std/4.0_r_std - ((jlambda**2.0_r_std*u_star**4.0_r_std)/speed**4.0_r_std)))

    DO i=1, jnlvls
       sigma_exp_1(i) = -jlambda * eta_z(jnlvls) * ( 1.0_r_std - eta_z(i) / eta_z(jnlvls))
       sigma_exp_2(i) = -3.0_r_std * small_n * ( 1.0_r_std - ( eta_z(i) / eta_z(jnlvls)))
       sigma_exp_3(i) = -jlambda * eta_z(jnlvls) * (1.0_r_std - eta_z(i) / eta_z(jnlvls))
       sigma_e_over_u_star(i) = ( nu_3 * &
            & EXP (sigma_exp_1(i)) + bigb1 *( EXP( sigma_exp_2(i)) - &
            & EXP( sigma_exp_3(i))))**(1.0_r_std/3.0_r_std)
       sigma_u_over_u_star(i) = gamma_1 * nu_1 * sigma_e_over_u_star(i)
       sigma_v_over_u_star(i) = gamma_2 * nu_1 * sigma_e_over_u_star(i)
       sigma_w_over_u_star(i) = gamma_3 * nu_1 * sigma_e_over_u_star(i)
       sigma_u(i) = sigma_u_over_u_star(i) * u_star
       sigma_v(i) = sigma_v_over_u_star(i) * u_star
       sigma_w(i) = sigma_w_over_u_star(i) * u_star
       u_speed(:,i) = speed * exp (-small_n * (1-(eta_z(i)/eta_z(jnlvls))))

       !++++SET a MINIMUM wind speed  +++++
       u_speed(:,i) = MAX (0.05_r_std, u_speed(:,i) ) 
       !++++++++++++++++++++++++++++++++++++   
    END DO !i=1, jnlvls


    IF (jnlvls .eq. 1) THEN
       u_speed(:,1) = speed
    END IF

    ! approximate the Lagrangian timescale using the parameterisation of Raupach (1989)
    ! where 0.4_r_std is the von Karmen constant

    DO i=1, jnlvls
       t_lagrangian(i) = (canopy_height/u_star) *                                      &
            &    (max (rp_c_0, 0.4_r_std * ( z_h(i) - (0.63_r_std*canopy_height) )   &
            &     / (rp_a_1 * canopy_height)))
    END DO !i=1, jnlvls

    ! this is the eddy diffusivity approximation for the far-field (also Raupach, 1989)

    DO i=1, jnlvls

       k_eddy(i) = (sigma_w(i))**2 * t_lagrangian(i)!/3.

    END DO !i=1, jnlvls

    IF (printlev_loc>=4) THEN
       WRITE(numout,*) 'k_eddy before', k_eddy(:)
       WRITE(numout,*) 'sigma_w', sigma_w(:)
       WRITE(numout,*) 't_lagrangian', t_lagrangian(:)
    ENDIF
    !+++ FOR NIGHTTIME K_EDDY was parameterized as a S-shape funtion based on threshood ustar value ++++
    ! Here we use a paramterization method for scaling the K_eddy based on u*
    ! 
    ! u_star < 0.3 cut 50% 
    ! u_star > 0.3 increase from 50 %
    ! K_EDDY_SLOPE determine the increasing/decreasing rate 
!!$    !     
!!$    k_eddy = k_eddy * (1.0 / (1.0 + EXP (-k_eddy_slope * (u_star- k_eddy_ustar))) ) 
    !++++++++++++++++++++++++++++++++++++++

  END SUBROUTINE mleb_profile


  !!  ===========================================================================================================================
  !! SUBROUTINE                  : mleb_stomatal_resistance
  !!
  !!
  !>\BRIEF                       : Calculates the resistance to the latent heat flux. 
  !!
  !! DESCRIPTION                 : The resistance to the latent heat flux is calucalated as the sum of the boundary-layer 
  !! resistance and the leaf stomatla resistance (see Ryder et al., 2016 section 3.2). Albeit the equation for big_r_prime 
  !! is a tiny bit different in the model, than the paper - an extra term is present in the code below, likely related to
  !! the boundary layer resistant.
  !! 
  !!
  !! RECENT CHANGE(S)            : None
  !!
  !! TO BE DONE                  : Much code is repeated below the impose_can_structure true or false. This should be reduced.
  !! Moreover, big_r_prime are possible overwritten, although the equations are almost identical, the differences and correctness
  !! should be verified. Should change down_sw_tab and down_lw_tab from out variables to local.
  !!
  !! MAIN OUTPUT VARIABLE(S)     :     down_sw_tab (but not used else where), down_lw_tab (but not used else where), big_r_prime
  !!
  !! REFERENCE(S)                :     Gao, W., Wesely, M.L. & Doskey, P.V., 1993. Numerical modeling of the turbulent 
  !!                                   diffusion and chemistry of NOx, O3, Isoprene, and other reactive trace gases in 
  !!                                   and above a forest canopy. Journal of Geophysical Research, 98(D10), pp.18339–18353
  !!
  !!                                   Guimberteau, M., 2010. Modélisation de l’hydrologie continentale et influences de 
  !!                                   l'irrigation sur le cycle de l'eau.
  !!
  !!                                   Monteith J., Unsworth M., 2008. Principles of environmental physics. Published 
  !!                                   Academic Press (Elsevier), ISBN: 978-0-12-505103-3
  !!
  !! FLOWCHART           : 
  !_ ==============================================================================================================================


  ! -------------------------------------------------------------------------------------------------

  SUBROUTINE mleb_stomatal_resistance (pad_deltah, swnet, lwdown, pb, &
       & down_sw_tab, down_lw_tab, big_r_prime, temp_atmos_pres_grid, &
       & q_atmos_pres_grid, big_r, u_speed, r_sto_v_fac, big_r_h2o_fac, &
       & rau, kjpindex, sr_fac, q_cdrag, delta_h, jrpad, &
       & profile_rveget)




    IMPLICIT NONE


    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                                   :: kjpindex         !! Domain size (-)
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)                :: rau              !! Air density (kg m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)                 :: swnet            !! Net surface short-wave flux (W m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)                 :: lwdown           !! Down-welling long-wave flux (W m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)                 :: pb               !! Lowest level pressure (hPa)
    REAL(r_std), DIMENSION (0:jnlvls), INTENT(in)                :: pad_deltah       !! Plant Area Density multiplied by layer height
    REAL(r_std), DIMENSION (kjpindex, jnlvls), INTENT(in)        :: temp_atmos_pres_grid
    REAL(r_std), DIMENSION (kjpindex, jnlvls), INTENT(in)        :: q_atmos_pres_grid
    REAL(r_std), DIMENSION(jnlvls), INTENT(in)                   :: big_r    
    REAL(r_std), DIMENSION(kjpindex,jnlvls), INTENT(in)                   :: u_speed
    REAL(r_std), INTENT(in)                                      :: r_sto_v_fac, big_r_h2o_fac    
    REAL(r_std), INTENT(in)                                      :: sr_fac           !! parameterisation factor   
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)               :: q_cdrag          !! Surface drag coefficient (-)
    REAL(r_std), DIMENSION(jnlvls), INTENT(in)                   :: delta_h          !! height of each level (m)
    REAL(r_std),DIMENSION (kjpindex,nvm,nlevels_tot), INTENT(in) :: profile_rveget
    REAL(r_std), DIMENSION (0:jnlvls), INTENT(in)                :: jrpad            !! Plant Area Density, replacement for values above
    !!   (as definitions between different
    !! parts of the model was becoming confusing)
    !!   'jrpad' is the total one-sided leaf area per unit of layer volume

    !! 0.2 Output variables
    REAL(r_std), DIMENSION (jnlvls), INTENT (out)                :: down_sw_tab, down_lw_tab
    REAL(r_std), DIMENSION (jnlvls), INTENT(out)                 :: big_r_prime


    !! 0.3 Modified variables


    !! 0.4 Local variables
    INTEGER(i_std)                                               :: i, j, jk             !! loop counters

    REAL(r_std), DIMENSION (jnlvls)                              :: temp_atmos_pres    
    REAL(r_std), DIMENSION (jnlvls)                              :: q_atmos_pres    
    ! constants from Tetens' expression for saturated vapor pressure (see Monteith &
    !   Unsworth 2.1)  bigA (-), tstar (K), tprime(K)
    REAL(r_std)                                                  :: biga_const = 17.27_r_std
    REAL(r_std)                                                  :: tstar_const = 273.0_r_std
    REAL(r_std)                                                  :: tprime_const = 36.0_r_std 
    REAL(r_std)                                                  :: e_sat_tstar = 0.611_r_std  !vapor pressure at T* (kPa)
    !! laitab is the accumulated LAI from the top of the canopy (see Gao et al., 1993)
    REAL(r_std), DIMENSION (jnlvls)                              :: laitab
    REAL(r_std), DIMENSION (jnlvls)                              :: e_sat_tair 
    !! extinction coefficient used is that from from src_parameters of 0.5_r_std for all PFTs               
    REAL(r_std)                                                  :: ext_coef_sw = 0.5_r_std
    REAL(r_std)                                                  :: ext_coef_lw = 0.8_r_std
    ! this is the assumed minimum average stomatal resistance (s/m)
    REAL(r_std)                                                  :: k_0 = 25.0d-2
    ! PFT dependent number that is characterised by vegetation (value for PFT8, 
    !   from Guimberteau, Table 2.2) (-)
    REAL(r_std)                                                  :: bigr_0 = 125.0_r_std ! solar radiation constant (Wm^-2)
    REAL(r_std)                                                  :: lambda = 1.5d-3
    ! coefficient for the increase of stomatal closure under hydrological deficit 
    !   (m^2 s kg^(-1))
    REAL(r_std), DIMENSION (jnlvls)                              :: delta_c
    ! deficit of water in the atmosphere at the level in question (kg m^-3)
    REAL(r_std), DIMENSION (jnlvls)                              :: q_sat_ta, r_sto_v
    REAL(r_std)                                                  :: big_l =2.46d6           !(J/kg)
    REAL(r_std)                                                  :: r_b = 30.0_r_std            ! leaf boundary layer resistance
    ! (see Guimberteau section 4.2)
    REAL(r_std)                                                  :: r_sv = 25.0_r_std           ! structural resistance for PFT8 (s/m) 
    REAL(r_std)                                                  :: r_a = 0.0_r_std             ! aerodynamic resistance represented by the eddy
    REAL(r_std)                                                  :: bigu_sv = 0.8_r_std         ! water uptake function
    REAL(r_std)                                                  :: small_d                 ! characteristic leaf length (m)
    REAL(r_std), DIMENSION(jnlvls)                               :: bigd_h2o                ! molecular diffusivity of water
    REAL(r_std), DIMENSION(jnlvls)                               :: nu_air                  ! kinematic viscosity of air
    REAL(r_std), DIMENSION(jnlvls)                               :: schmidt                 ! schmidt number
    REAL(r_std), DIMENSION(jnlvls)                               :: sherwood                ! sherwood number
    REAL(r_std), DIMENSION(jnlvls)                               :: reynolds                ! reynolds number
    REAL(r_std), DIMENSION(jnlvls)                               :: big_r_h2o               ! h2o leaf boundary layer resistance
    REAL(r_std)                                                  :: smalla = 23.0_r_std        ! this is the assumed minimum average stomatal resistance (s/m)

  !! ==============================================================================================================================
    
    
    !! This routine calculates the resistance to latent heat flux between the leaves and the atmosphere inside the canopy.
    !! The first resistance (big_r_h2o) corresponds to the boundary layer resistance around the leaves and the second one 
    !! to the stomatal conductance (r_sto_v). The ressistance to latent heat flux big_r_prime is the sum of both.

    !! Initialisation:

    r_sto_v=zero
    big_r_h2o=zero


    !---initialized the big_r_prime
    big_r_prime(:)= 1.0d3
    !----


    !!J.Alléon: This part of the code should be cleaned up

    !! TO DO: Check if the flag has to be removed or not.

    IF (.NOT. hack_can_structure) THEN 

       DO jk = 1, kjpindex

          temp_atmos_pres(:) = temp_atmos_pres_grid(jk,:)
          q_atmos_pres(:) = q_atmos_pres_grid(jk,:)

          small_d = 0.05_r_std                            ! set the leaf length to 5cm for now
          ! physical dimensions for aspen leaf length: www.ag.ndsu.edu/trees/handbook/th-3-103.pdf

          ! n.b. for Eucalyptus delegatensis:
          ! http://davisla.wordpress.com/2012/04/07/plant-of-the-week-eucalyptus-delegatensis/
          ! 'juvenile leaf being up to 25cm long and 10cm broad the mature 18cm long and 3cm broad'


          ! calculate steps of LAI (using the same formulation as in diffuco.f90 for co2)
          laitab = 0.0_r_std

          ! calculate light fraction that comes through at a given LAI for each vegetation type
          ! laitab(jnlvls) = 0

          IF (jnlvls .eq. 1) THEN
             ! omit the layered calculation for single layer case
          ELSE
             
             !! TO DO: The part between "###" should be removed.

             !! ################

             DO i=jnlvls-1, 1, -1

                laitab(i) = (jrpad(i+1) * delta_h(i+1) ) + laitab(i+1)             

             END DO !i=jnlvls, 1, -1

             DO i=1, jnlvls
                down_sw_tab(i) = swnet(1) * exp (-ext_coef_sw * laitab(i))
                down_lw_tab(i) = lwdown(1) * exp (-ext_coef_lw * laitab(i))
             END DO !i=1, jnlvls


             DO i=1, jnlvls
                q_sat_ta(i) = (0.622_r_std/ pb(1))  *e_sat_tstar* 1000.0_r_std * &
                     & EXP((biga_const*(temp_atmos_pres(i) - tstar_const) / &
                     & (temp_atmos_pres(i) - tprime_const)))

                delta_c(i) = rau(1) * MAX (q_sat_ta(i) - q_atmos_pres(i), zero) /100.0_r_std
             END DO !i=1, jnlvls

             r_sto_v(:) = 1.0d20

             DO i=1, nlevels_tot
                IF (down_sw_tab(i) .ne. 0.0d0) THEN
                   r_sto_v(i+jnlvls_under) = (smalla/k_0) * (1/pad_deltah(i+jnlvls_under)) * &
                        & ((down_sw_tab(i+jnlvls_under)+bigr_0)/down_sw_tab(i+jnlvls_under))  &
                        & * (un + lambda * (delta_c(i+jnlvls_under)/smalla) )
                ELSE
                   r_sto_v(i) = 1.0d20
                END IF
             END DO !i=1,jnlvls


             IF(printlev_loc >= 4) THEN
                WRITE(numout, *) '100316 nlevels_tot is: ', nlevels_tot
                WRITE(numout, *) '100316 (jnlvls_under + nlevels_tot) is: ', (jnlvls_under + nlevels_tot)
             ENDIF

             !! Debug ID6.
             !! The expression linking the stomatal resistance calculated in diffuco and the one used in mleb is re-established.
             !! The calculation above should be removed.

             IF (.NOT. hack_can_structure) THEN 

             !! ###############
                
                !! Initialisation: except inside the canopy, the stomatal conductance should be set to infinity.

!!$                r_sto_v(:) = 1.0d20

                !! TO DO: To remove.  

                ! r_sto_v(10:19) = profile_rveget(ji,jv,1:nlevels_tot) 
                ! r_sto_v(jnlvls_under:(jnlvls_under + nlevels_tot)) = profile_rveget(ji,jv,1:nlevels_tot) !asla 

                !! Calculation of the stomatal resistance:

                !! The stomatal conductance in each layer is calculated in diffuco.f90.

!!$                DO i=1,nlevels_tot
!!$                   IF (profile_rveget(jk,test_pft,i) .GT. zero) THEN
!!$                      r_sto_v(i+jnlvls_under) = profile_rveget(jk,test_pft,i) 
!!$                      !! Changer 6 et faire une boucle pft 
!!$                   ELSE
!!$                      r_sto_v(i+jnlvls_under) = 1.0d20
!!$                   ENDIF
!!$                ENDDO
                

             END IF ! (hack_can_structure) THEN 

             !! End ID6

             !! Debug.

             IF(printlev_loc >=4 ) THEN
                DO i = 1, jnlvls
                   WRITE (numout,*) 'r_sto_v(i) is: ', r_sto_v(i)
                END DO
             ENDIF

             !! TO DO: Remove this comment.

             ! We're using the expression from ORCHIDEE as in Guimberteau (2010)


             !! Calculation of the boundary layer resistance:
             !! /!\ This boundary layer resistance is different from the resistance to sensible heat flux as it 
             !!     set for the latent heat flux.

             DO i = 1, jnlvls

                !! NB: The transition point from laminar to turbulent flow occurs when the 
                !!     Reynolds number passes 8,000 to 25,000 (Baldochhi, 1988)

                !! Debug ID7.
                !! Previous expressions of air diffusivity and viscosity where causing problems.
                !! Changed for other expressions (see GGSheet)

!!$                bigd_h2o(i) = 2.26d-5 + (1.51d-7 * (temp_atmos_pres(i) - 273.15_r_std))
!!$                nu_air(i) = 1.35d-5 + (1.0d-7 * (temp_atmos_pres(i) - 273.15_r_std))

                bigd_h2o(i) = MAX(1d-6, (1.0_r_std + 7d-3 * (temp_atmos_pres(i) - 273.15_r_std)) * 21.2d-6)
                nu_air(i) = MAX(1d-6, -1.363528d-14 * temp_atmos_pres(i)**3 + 1.00881778d-10 * temp_atmos_pres(i)**2 &
                     & + 3.452139d-8 * temp_atmos_pres(i) - 3.400747d-6)

                !! End ID7

                schmidt(i) = nu_air(i) / bigd_h2o(i)

                reynolds(i) = (small_d * u_speed(jk,i))/nu_air(i)

                ! is the flow laminar or turbulent?
                IF (reynolds(i) .le. 8000.0_r_std) THEN
                   sherwood(i) = 0.66_r_std * (reynolds(i))**(0.5_r_std) + (schmidt(i))**(0.33_r_std)
                ELSE IF (reynolds(i) .gt. 8000.0_r_std) THEN
                   sherwood(i) = 0.03_r_std * (reynolds(i))**(0.8_r_std) + (schmidt(i))**(0.33_r_std)
                END IF

                big_r_h2o(i) = small_d / (bigd_h2o(i) * sherwood(i))

                !! TO DO: To remove

                !! Debug ID8.
                !! Line commented in order to be the same as in the article. To comment out if we want to have the same as
                !! in diffuco. We should look at which one is better.

                   ! +++ TEMP +++ use the same expression of leaf boundary+++++++++++++++++++++++++++++++++++++++
                   ! resistance in the subroutineL diffuco_trans_co2 
!!$                       big_r_h2o(i) = (1.0_r_std /25.0_r_std)/(22.4_r_std*temp_atmos_pres(i)/273._r_std/1000._r_std) 
!!$                       !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$                       big_r_prime(i) = ((r_sto_v(i) * r_sto_v_fac) + (big_r_h2o(i) * big_r_h2o_fac)) !/ pad_deltah(i)


                !! End ID8

             END DO ! i=1, jnlvls

          END IF ! (jnlvls .eq. 1)

          !! Debug ID9.
          !! big_r_h2o already takes into account the boundary layer resistance for water. We should not add big_r as it was before.

          !! The final resistance is the sum of the stomatal resistance and the boundary layer one.

          big_r_prime = (r_sto_v * r_sto_v_fac) + (big_r_h2o * big_r_h2o_fac)

          !! End ID9

          
          big_r_prime = big_r_prime * sr_fac

          big_r_prime(jnlvls) = 1.0d20   
          !! NB: The top level resistance should be approaching 
          !!     infinity (see diagram) need to implement this 
          !!     more elegantly
          !! TO DO: To remove as the resistance is set to infinity at the beginning of the routine.

          IF (printlev_loc>=4) WRITE(numout, *) 'big_r_prime is: ', big_r_prime


       END DO ! jk = 1, kjpindex



    ELSE ! (.NOT. hack_can_structure)

       !! TO DO: Check if this is still needed.

       DO jk = 1, kjpindex

          temp_atmos_pres(:) = temp_atmos_pres_grid(jk,:)
          q_atmos_pres(:) = q_atmos_pres_grid(jk,:)

          small_d = 0.05_r_std                            ! set the leaf length to 5cm for now
          ! physical dimensions for aspen leaf length: www.ag.ndsu.edu/trees/handbook/th-3-103.pdf

          ! n.b. for Eucalyptus delegatensis:
          ! http://davisla.wordpress.com/2012/04/07/plant-of-the-week-eucalyptus-delegatensis/
          ! 'juvenile leaf being up to 25cm long and 10cm broad the mature 18cm long and 3cm broad'


          ! calculate steps of LAI (using the same formulation as in diffuco.f90 for co2)
          laitab = 0.0_r_std

          ! calculate light fraction that comes through at a given LAI for each vegetation type
          ! laitab(jnlvls) = 0

          IF (jnlvls .eq. 1) THEN
             ! omit the layered calculation for single layer case
          ELSE

             DO i=jnlvls-1, 1, -1

                ! laitab(i) = (pad_deltah(i+1) * delta_h(i+1)) + laitab(i+1)   
                ! laitab(i) = (pad_deltah(i+1) ) + laitab(i+1)             
                laitab(i) = (jrpad(i+1) * delta_h(i+1) ) + laitab(i+1)             

             END DO !i=jnlvls, 1, -1

             DO i=1, jnlvls
                down_sw_tab(i) = swnet(1) * exp (-ext_coef_sw * laitab(i))
                down_lw_tab(i) = lwdown(1) * exp (-ext_coef_lw * laitab(i))
             END DO !i=1, jnlvls

             IF (printlev_loc>=4) THEN

                DO i=jnlvls-1, 1, -1
                   WRITE (numout,*) 'laitab, level: ', i, ' is: ', laitab(i)
                END DO !i=1, jnlvls

                DO i=jnlvls-1, 1, -1
                   WRITE (numout,*) 'down_sw_tab, level: ', i, ' is: ', down_sw_tab(i)
                END DO !i=1, jnlvls

                DO i=jnlvls-1, 1, -1
                   WRITE (numout,*) 'down_lw_tab, level: ', i, ' is: ', down_lw_tab(i)
                END DO !i=1, jnlvls

             END IF ! (printlev_loc>=4)       

             DO i=1, jnlvls
                q_sat_ta(i) = (0.622_r_std/ pb(1))  *e_sat_tstar* 1000.0_r_std * &
                     & EXP((biga_const*(temp_atmos_pres(i) - tstar_const) / &
                     & (temp_atmos_pres(i) - tprime_const)))

                delta_c(i) = rau(1) * MAX (q_sat_ta(i) - q_atmos_pres(i), zero) /100.0_r_std
             END DO !i=1, jnlvls

             DO i=1, jnlvls
                IF (down_sw_tab(i) .ne. 0.0d0) THEN
                   r_sto_v(i) = (smalla/k_0) * (1/pad_deltah(i)) * &
                        & ((down_sw_tab(i)+bigr_0)/down_sw_tab(i))  &
                        & * (un + lambda * (delta_c(i)/smalla) )
                ELSE
                   r_sto_v(i) = 1.0d20
                END IF
             END DO !i=1,jnlvls


             ! We're using the expression from ORCHIDEE as in Guimberteau (2010)

             DO i = 1, jnlvls

                ! n.b. the transition point from laminar to turbulent flow occurs when the 
                !   Reynolds number passes 8,000 to 25,000 (Baldochhi, 1988)


                !! Debug ID7.
                !! Previous expressions of air diffusivity and viscosity where causing problems.
                !! Changed for other expressions (see GGSheet)

!!$                bigd_h2o(i) = 2.26d-5 + (1.51d-7 * (temp_atmos_pres(i) - 273.15_r_std))
!!$                nu_air(i) = 1.35d-5 + (1.0d-7 * (temp_atmos_pres(i) - 273.15_r_std))

                bigd_h2o(i) = MAX(1d-6, (1.0_r_std + 7d-3 * (temp_atmos_pres(i) - 273.15_r_std)) * 21.2d-6)
                nu_air(i) = MAX(1d-6, -1.363528d-14 * temp_atmos_pres(i)**3 + 1.00881778d-10 * temp_atmos_pres(i)**2 &
                     & + 3.452139d-8 * temp_atmos_pres(i) - 3.400747d-6)

                !! End ID7

                schmidt(i) = nu_air(i) / bigd_h2o(i)

                reynolds(i) = (small_d * u_speed(jk,i))/nu_air(i)

                ! is the flow laminar or turbulent?
                IF (reynolds(i) .le. 8000.0_r_std) THEN
                   sherwood(i) = 0.66_r_std * (reynolds(i))**(0.5_r_std) + (schmidt(i))**(0.33_r_std)
                ELSE IF (reynolds(i) .gt. 8000.0_r_std) THEN
                   sherwood(i) = 0.03_r_std * (reynolds(i))**(0.8_r_std) + (schmidt(i))**(0.33_r_std)
                END IF

                big_r_h2o(i) = small_d / (bigd_h2o(i) * sherwood(i))

                big_r_prime(i) = ((r_sto_v(i) * r_sto_v_fac) + (big_r_h2o(i) * big_r_h2o_fac)) !/ pad_deltah(i)

             END DO ! i=1, jnlvls

          END IF ! (jnlvls .eq. 1)

          !! Debug ID9.
          !! big_r_h2o already takes into account the boundary layer resistance for water. We should not add big_r as it was before.

          big_r_prime = ((r_sto_v * r_sto_v_fac) + (big_r_h2o * big_r_h2o_fac))

          !! End ID9
          
          big_r_prime = big_r_prime * sr_fac

          big_r_prime(jnlvls) = 1.0d20   ! the top level resistance should be approaching 
          !   infinity (see diagram) need to implement this 
          !   more elegantly

          IF (printlev_loc>=4) WRITE (numout, *) 'big_r_prime is: ', big_r_prime

       END DO ! jk = 1, kjpindex

    END IF ! (hack_can_structure)

  END SUBROUTINE mleb_stomatal_resistance


  ! -------------------------------------------------------------------------------------------------

  !!  ===========================================================================================================================
  !! SUBROUTINE                  : mleb_boundarylayer_resistance
  !!
  !!
  !>\BRIEF                       : Calculates the boundary layer resistance. 
  !!
  !! DESCRIPTION                 : The resistance to the sensible heat flux is equal to the boundary-layer 
  !! resistance (Rb,i) at the leaf surface, thus this subroutine essentially gives us the resistance to the sensible heaf flux, 
  !! which is needed in the  mleb_stomatal_resistance subroutine to calculated the resistance to the latent heaf flux (see Ryder et al.,
  !! 2016, sec. 3.2 eqn. 5 to 8).  
  !! 
  !!
  !! RECENT CHANGE(S): None
  !!
  !! MAIN OUTPUT VARIABLE(S)     : big_r
  !!
  !! REFERENCE(S)                :     Baldocchi, D.D., 1988. A multi-layer model for estimating sulfur dioxide deposition 
  !!                                      to a deciduous oak forest canopy. Atmospheric Environment (1967), 22(5), 
  !!                                      pp.869–884
  !!
  !! FLOWCHART                   : 
  !_ ==============================================================================================================================



  SUBROUTINE mleb_boundarylayer_resistance (big_r, u_speed, temp_atmos_pres_grid, &
       & kjpindex, br_fac, q_cdrag, &
       & pad_deltah, jrpad)


    IMPLICIT NONE

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                                     :: kjpindex             !! Domain size (-)
    REAL(r_std), DIMENSION (kjpindex, jnlvls), INTENT(in)          :: temp_atmos_pres_grid !! atmospheric temperature down the column (K) PRESENT STEP
    REAL(r_std), DIMENSION (kjpindex,jnlvls), INTENT(in)                    :: u_speed
    REAL(r_std), INTENT(in)                                        :: br_fac
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)                 :: q_cdrag              !! Surface drag coefficient (-)
    REAL(r_std), DIMENSION (0:jnlvls), INTENT(in)                  :: pad_deltah           !! Plant Area Density multiplied by layer height
    REAL(r_std), DIMENSION (0:jnlvls), INTENT(in)                  :: jrpad                !! Plant Area Density, replacement for values above (as definitions between different
    !!    parts of the model was becoming confusing)
    !!   'jrpad' is the total one-sided leaf area per unit of layer volume


    !! 0.2 Output variables
    REAL(r_std), DIMENSION(jnlvls), INTENT(out)                    :: big_r                !! Boundary layer resistance (s m-1)   

    !! 0.3 Modified variables


    !! 0.4 Local variables

    REAL(r_std), DIMENSION (jnlvls)                                :: temp_atmos_pres      !! atmospheric temperature down the column (K) PRESENT STEP
    REAL(r_std)                                                    :: r_a = 0.0_r_std      !! aerodynamic resistance is represented by the eddy
    ! diffusivity coefficients so I have set this to zero 
    REAL(r_std)                                                    :: r_b = 30.0_r_std     !! leaf boundary layer resistance
    INTEGER(i_std)                                                 :: i, jk
    REAL(r_std)                                                    :: small_d              !! characteristic leaf length (m)
    REAL(r_std)                                                    :: big_l                !! characteristic leaf dimension (m)
    REAL(r_std)                                                    :: big_d                !! molecular diffusivity of air
    REAL(r_std)                                                    :: sherwood             !! Sherwood number
    REAL(r_std)                                                    :: prandtl              !! Prandtl number 
    REAL(r_std), DIMENSION (jnlvls)                                :: nusselt              !! Nusselt number
    REAL(r_std), DIMENSION (jnlvls)                                :: bigd_h               !! heat diffusivity of air (m^2 s^-1)
    REAL(r_std), DIMENSION (jnlvls)                                :: nu_air               !! the kinematic viscosity of air
    REAL(r_std), DIMENSION (jnlvls)                                :: reynolds             !! reynolds number
    REAL(r_std), DIMENSION(jnlvls)                                 :: big_r_heat           !! Resistance to the sensible heat flux (s m-1) 

  ! ==============================================================================================================================


    !! This routine calculates the resistance to sensible heat flux between the leaves and the atmosphere in the canopy. 
    !! The model uses an equation developed by Baldocchi et al. (1988).


    DO jk = 1, kjpindex

       temp_atmos_pres(:) = temp_atmos_pres_grid(jk,:)

       IF (jnlvls .eq. 1) THEN

          big_r(:)=un

          ! omit the layered calculation for single layer case
       ELSE

          !! TO DO: Clean this commented section.

          small_d = 0.01_r_std                                ! set the leaf length to 1cm for now
          big_l = 0.01_r_std                                  ! set the leaf length to 1cm for now
          prandtl = 1.0_r_std                                 ! set the Prandtl number to 1 for now
          !            DO i = 1, jnlvls
          !                big_r(i) = 200.0_r_std * (0.05_r_std / u_speed (i))
          !            END DO ! i = 1, jnlvls

          !! Debug ID10.
          !! big_r is set to 10^20 everywhere except inside the canopy, the loop over all layers is reduced to a loop inside the canopy.

          big_r_heat(:)=1.0d20

          DO i = 1, nlevels_tot

             !! TO DO: Check if the previous version of the variables is working.

             !! Debug ID11 (same ID7)

!!$             bigd_h(i+jnlvls_under) = 1.9d-5 + (1.26d-5 * (temp_atmos_pres(i+jnlvls_under) - 273.15_r_std))
!!$             nu_air(i+jnlvls_under) = 1.35d-5 + (1.0d-7 * (temp_atmos_pres(i+jnlvls_under) - 273.15_r_std))


             bigd_h(i+jnlvls_under) = MAX(1d-6, (1.0_r_std + 7d-3 * (temp_atmos_pres(i+jnlvls_under) - 273.15_r_std)) * 18.9d-6)
             nu_air(i+jnlvls_under) = MAX(1d-6, -1.363528d-14 * temp_atmos_pres(i+jnlvls_under)**3 + 1.00881778d-10 * temp_atmos_pres(i+jnlvls_under)**2 &
                  & + 3.452139d-8 * temp_atmos_pres(i+jnlvls_under) - 3.400747d-6)

             !! End ID11

             reynolds(i+jnlvls_under) = (big_l * u_speed(jk,i+jnlvls_under)) / nu_air(i+jnlvls_under)

             !! Is the flow laminar or turbulent?
             IF (prandtl .lt. 0.6_r_std) THEN
                nusselt(i+jnlvls_under) = 0.332_r_std * (reynolds(i+jnlvls_under))**(1.0_r_std/2.0_r_std) + (prandtl)**(1.0_r_std/3.0_r_std)
             ELSE IF (prandtl .ge. 0.6_r_std) THEN
                nusselt(i+jnlvls_under) = 0.0296_r_std * (reynolds(i+jnlvls_under))**(4.0_r_std/5.0_r_std) + (prandtl)**(1.0_r_std/3.0_r_std)
             END IF

             big_r_heat(i+jnlvls_under) = small_d / (bigd_h(i+jnlvls_under) * nusselt(i+jnlvls_under))

          END DO ! i=1, jnlvls

          !! End ID10

       END IF ! (jnlvls .eq. 1)

       !! The variable big_r is the one used in the column energy budget calculation.

       big_r = big_r_heat

       big_r = big_r * br_fac

       big_r(jnlvls) = 1.0d20     
    
       !! NB: The top level resistance should be approaching 
       !!     infinity (see diagram) need to implement this more 
       !!     elegantly.
       !! TO DO: This is already set at the beginning of the routine, should be removed.

    END DO ! jk = 1, kjpindex

  END SUBROUTINE mleb_boundarylayer_resistance


  !!  ====================================================================================
  !! SUBROUTINE                 : mleb_begin
  !!
  !>\BRIEF                  Preliminary variables required for the 
  !! calculation of the energy budget are derived.
  !!
  !! DESCRIPTION                : This routines computes preliminary 
  !! variables required for the calculation of the energy budget: the old surface static 
  !! energy (psold), the surface saturation humidity (qsol_sat), the derivative of satured 
  !! specific humidity at the old temperature (pdqsold) and the net radiation (netrad).
  !!
  !! MAIN OUTPUT VARIABLE(S)                    : psold, qsol_sat, pdqsold, netrad
  !!
  !! REFERENCE(S)               :
  !! - Best, MJ, Beljaars, A, Polcher, J & Viterbo, P, 2004. A proposed structure for 
  !! coupling tiled surfaces with the planetary boundary layer. Journal of 
  !! Hydrometeorology, 5, pp.1271-1278
  !! - de Noblet-Ducoudré, N, Laval, K & Perrier, A, 1993. SECHIBA, a new set of 
  !! parameterisations of the hydrologic exchanges at the land-atmosphere interface within
  !! the LMD Atmospheric General Circulation Model. Journal of Climate, 6, pp.248-273
  !! - Dufresne, J-L & Ghattas, J, 2009. Description du schéma de la couche limite 
  !! turbulente et la interface avec la surface planetaire dans LMDZ, Technical note, 
  !! available (22/12/11):
  !! http://lmdz.lmd.jussieu.fr/developpeurs/notes-techniques/ressources/pbl_surface.pdf
  !! - Polcher, J. McAvaney, B, Viterbo, P, Gaertner, MA, Hahmann, A, Mahfouf, J-F, 
  !! Noilhan, J Phillips, TJ, Pitman, AJ, Schlosser, CA, Schulz, J-P, Timbal, B, Verseghy,
  !! D Xue, Y, 1998. A proposal for a general interface between land surface schemes and
  !! general circulation models. Global and Planetary Change, 19, pp.261-276
  !! - Richtmeyer, RD, Morton, KW, 1967. Difference Methods for Initial-Value Problems.
  !! Interscience Publishers\n
  !! - Schulz, Jan-Peter, Lydia Dümenil, Jan Polcher, 2001: On the Land Surface–Atmosphere 
  !! Coupling and Its Impact in a Single-Column Atmospheric Model. J. Appl. Meteor., 
  !! 40, 642–663.
  !!
  !! FLOWCHART  : None                     
  !! \n
  !_ =====================================================================================

  SUBROUTINE mleb_begin (kjpindex, temp_sol, lwdown, swnet, pb, psold, qsol_sat, &
       & pdqsold, netrad, emis, dev_qsol,lwabs)

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                         :: kjpindex         !! Domain size (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: temp_sol         !! Surface temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: lwdown           !! Down-welling long-wave flux (W m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: swnet            !! Net surface short-wave flux (W m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: pb               !! Lowest level pressure (hPa)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: emis             !! Emissivity (-)


    !! 0.2 Output variables

    REAL(r_std),DIMENSION (kjpindex), INTENT (out)     :: psold            !! Old surface dry static energy (J kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)     :: qsol_sat         !! Saturated specific humudity for old temperature 
    !! (kg kg^{-1})    
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)     :: pdqsold          !! Derivative of satured specific humidity at the old 
    !! temperature (kg (kg s)^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)     :: netrad           !! Net radiation (W m^{-2})
    REAL(r_std), DIMENSION(kjpindex), INTENT (out)     :: dev_qsol
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)      :: lwabs           !! Absorbed long-wave flux (W m^{-2})

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                                     :: ji

    REAL(r_std), PARAMETER                             :: missing = 999998.
    !_ ================================================================================================================================


    !! TO DO: This routine is surely useless now. Should be removed.


    !  write(numout, *) 'mleb_begin lwdown:', lwdown(:)
    !! 1. Computes psold (the surface static energy for the old timestep)

    !! We here define the surface static energy for the 'old' timestep, in terms of the surface
    !! temperature and heat capacity.
    !! \latexonly 
    !!     \input{enerbilbegin1.tex}
    !! \endlatexonly
    psold(:) = temp_sol(:)*cp_air

    !! 2. Computes qsol_sat (the surface saturated humidity).

    !! We call the routine 'qsatcalc' from within the module 'src_parameters/constantes_veg'
    IF(printlev_loc >= 4) WRITE(numout,*)'271016 mleb_begin, temp_sol', temp_sol
    CALL qsatcalc (kjpindex, temp_sol, pb, qsol_sat)

    IF ( diag_qsat ) THEN
       IF ( ANY(ABS(qsol_sat(:)) .GT. missing) ) THEN
          DO ji = 1, kjpindex
             IF ( ABS(qsol_sat(ji)) .GT. missing) THEN
                WRITE(numout,*) 'ERROR on ji = ', ji
                WRITE(numout,*) 'temp_sol(ji),  pb(ji) :', temp_sol(ji),  pb(ji)
                CALL ipslerr (3,'mleb_begin', &
                     &           'qsol too high ','','')
             ENDIF
          ENDDO
       ENDIF
    ENDIF

    !! 3. Computes pdqsold 

    !! Computes pdqsold (the derivative of the saturated humidity with respect to temperature
    !! analysed at the surface temperature at the 'old' timestep.
    !! We call the routine 'dev_qsatcalc' from within the module 'src_parameters/constantes_veg'.
    CALL dev_qsatcalc (kjpindex, temp_sol, pb, dev_qsol)

    IF (printlev_loc>=4) THEN
!!$       WRITE(numout,*) '***** 010213 *****'
!!$       WRITE(numout,*) ' '
!!$       WRITE(numout,*) 'after CALL dev_qsatcalc'
!!$       WRITE(numout,*) 'temp_sol(1) is: ', temp_sol(1)
!!$       WRITE(numout,*) 'pb(1) is: ', pb(1)
!!$       WRITE(numout,*) 'dev_qsol(1) is: ', dev_qsol(1)        
!!$       WRITE(numout,*) 'EFFECTIVE jalpha is: ', dev_qsol(1)
!!$       WRITE(numout,*) 'EFFECTIVE jbeta is: ', qsol_sat(1) - temp_sol(1) * dev_qsol(1)
    ENDIF  ! (printlev_loc>=4)


    !+++ According to r2093 in the trunk, this causes and instabilty in the coupling
    !    DO ji = 1, kjpindex
!!$    
!!$    !! \latexonly 
!!$    !!     \input{enerbilbegin2.tex}
!!$    !! \endlatexonly
    !      pdqsold(ji) = dev_qsol(ji) * ( pb(ji)**kappa ) / cp_air
    !    ENDDO

    pdqsold(:) = dev_qsol(:)
    !+++++++++

    IF ( diag_qsat ) THEN
       IF ( ANY(ABS( pdqsold(:)) .GT. missing) ) THEN
          DO ji = 1, kjpindex
             IF ( ABS( pdqsold(ji)) .GT. missing ) THEN
                WRITE(numout,*) 'ERROR on ji = ', ji
                WRITE(numout,*) 'temp_sol(ji),  pb(ji) :', temp_sol(ji),  pb(ji)
                CALL ipslerr (3,'enerbil_begin', &
                     &           'pdqsold too high ','','')
             ENDIF
          ENDDO
       ENDIF
    ENDIF


    !! 4. Computes the net radiation and the absorbed LW radiation absorbed at the surface. 

    !! Long wave radiation absorbed by the surface is the product of emissivity and downwelling LW radiation    
    !! \latexonly 
    !!     \input{enerbilbegin3.tex}
    !! \endlatexonly
    lwabs(:) = emis(:) * lwdown(:)

    !! Net radiation is calculated as:
    !! \latexonly 
    !!     \input{enerbilbegin4.tex}
    !! \endlatexonly
    netrad(:) = lwdown(:) + swnet (:) - (emis(:) * c_stefan * temp_sol(:)**4 + &
         & (un - emis(:)) * lwdown(:)) 
    !
    IF (printlev>=3) WRITE (numout,*) ' mleb_begin done '


    !+++ TEMP FOR NIGHTTIME K_EDDY_FAC ++++
    !DO ji =1, kjpindex
    !      IF (swnet(ji) <= 0 ) THEN
    !          k_eddy_fac =0.01
    !      ELSE
    !          CALL getin_p('K_EDDY_FAC', k_eddy_fac) 
    !      ENDIF
    !ENDDO
    !+++++++++++++++++++++++++++++++++++++


  END SUBROUTINE mleb_begin

  !===========================================================================================================================
  !! SUBROUTINE                   : mleb_evapveg
  !!
  !>\BRIEF                        : Diagnoses the values for evaporation and transpiration: 
  !! vevapsno (snow evaporation), vevapnu (bare soil evaporation), transpir (transpiration) and vevapwet
  !! (interception)
  !!
  !! DESCRIPTION                    : Based on the estimation of the fluxes in mleb_calc_column , 
  !! this routine splits the total evaporation into transpiration interception
  !! loss from vegetation, bare soil evaporation and snow sublimation. 
  !!
  !! MAIN OUTPUT VARIABLE(S)        : vevapsno, vevapnu, transpir,vevapwet, transpot, vevapflo
  !!
  !! REFERENCE(S)                   :
  !! - Best, MJ, Beljaars, A, Polcher, J & Viterbo, P, 2004. A proposed structure for coupling tiled
  !! surfaces with the planetary boundary layer. Journal of Hydrometeorology, 5, pp.1271-1278
  !! - de Noblet-Ducoudré, N, Laval, K & Perrier, A, 1993. SECHIBA, a new set of parameterisations
  !! of the hydrologic exchanges at the land-atmosphere interface within the LMD Atmospheric General
  !! Circulation Model. Journal of Climate, 6, pp.248-273
  !! - Dufresne, J-L & Ghattas, J, 2009. Description du schéma de la couche limite turbulente et la
  !! interface avec la surface planetaire dans LMDZ, Technical note, available (22/12/11):
  !!
  !http://lmdz.lmd.jussieu.fr/developpeurs/notes-techniques/ressources/pbl_surface.pdf
  !! - Guimberteau, M, 2010. Modélisation de l'hydrologie continentale et influences de l'irrigation
  !! sur le cycle de l'eau, PhD Thesis, available (22/12/11): 
  !! http://www.sisyphe.upmc.fr/~guimberteau/docs/manuscrit_these.pdf
  !! - Polcher, J. McAvaney, B, Viterbo, P, Gaertner, MA, Hahmann, A, Mahfouf, J-F, Noilhan, J
  !! Phillips, TJ, Pitman, AJ, Schlosser, CA, Schulz, J-P, Timbal, B, Verseghy, D &
  !! Xue, Y, 1998. A proposal for a general interface between land surface schemes and
  !! general circulation models. Global and Planetary Change, 19, pp.261-276
  !! - Richtmeyer, RD, Morton, KW, 1967. Difference Methods for Initial-Value Problems.
  !! Interscience Publishers
  !! - Schulz, Jan-Peter, Lydia Dümenil, Jan Polcher, 2001: On the Land Surface–Atmosphere 
  !! Coupling and Its Impact in a Single-Column Atmospheric Model. J. Appl. Meteor., 40, 642–663. 
  !!
  !! FLOWCHART   : None
  !! \n
  !_
  !==============================================================================================================================



  SUBROUTINE mleb_evapveg (kjpindex, vbeta1, vbeta2, vbeta3, vbeta3pot, &
       & vbeta4, vbeta5, rau, u, v, q_cdrag, qair, humrel, &
!!$     & vbeta4, vbeta5, cimean, ccanopy, rau, u, v, q_cdrag, qair, humrel, &
       & vevapsno, vevapnu , vevapflo, vevapwet, transpir, transpot, evapot,qsol_sat_new)

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                          :: kjpindex          !! Domain size (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: vbeta1            !! Snow resistance (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: vbeta4            !! Bare soil resistance (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: vbeta5            !! Floodplains resistance
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: rau               !! Density (kg m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: u, v              !! Wind velocity in directions u and v (m s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: q_cdrag           !! Surface drag coefficient (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: qair              !! Lowest level specific humidity (kg kg^{-1})
!!$    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: ccanopy           !! CO2 concentration in the canopy (ppm)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: evapot            !! Soil Potential Evaporation (mm/tstep)
    REAL(r_std),DIMENSION (kjpindex, nvm), INTENT (in)  :: humrel            !! Soil moisture stress (within range 0 to 1)
!!$ DS 15022011 humrel was used in a previous version of Orchidee, developped by Nathalie. Need to be discussed if it should be introduce again             
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)   :: vbeta2            !! Interception resistance (-)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)   :: vbeta3            !! Vegetation resistance (-)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)   :: vbeta3pot         !! Vegetation resistance for potential transpiration
!!$    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)   :: cimean            !! STOMATE: mean intercellular ci \mu mole m^{-2} s^{-1}
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)   :: qsol_sat_new         !! Air saturated humidity (kg.kg-1)

    !! 0.2 Output variables

    REAL(r_std),DIMENSION (kjpindex), INTENT (out)      :: vevapsno          !! Snow evaporation (mm day^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)      :: vevapnu           !! Bare soil evaporation (mm day^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)      :: vevapflo         !! Floodplains evaporation
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)  :: transpir          !! Transpiration (mm day^{-1})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)  :: transpot          !! Potential Transpiration
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)  :: vevapwet          !! Interception (mm day^{-1})

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                                      :: ji, jv
    REAL(r_std), DIMENSION(kjpindex)                    :: xx
    REAL(r_std), DIMENSION(kjpindex)                    :: vbeta2sum, vbeta3sum
    REAL(r_std)                                     :: speed         !! Speed (m s^{-1})
    !_ ==============================================================================================================================


    !! This routine aims at dividing the total evaporation flux inside its components according to the beta coefficients calculated in 
    !! diffuco.f90.
    !! This routine should maybe be reviewed if we consider the multi-layer energy budget only applicable to single canopies. 

    !! All the fluxes are calculated thanks to the corresponding beta (which takes into account fraction of grid cell and resistances to 
    !! evaporation flux) and the difference of humidity between the surface and the atmosphere.



    !! Initialisation: usefull to calculate the evaporation of floodplains where there is vegetation

    vbeta2sum(:) = 0.
    vbeta3sum(:) = 0.
    DO jv = 1, nvm
       vbeta2sum(:) = vbeta2sum(:) + vbeta2(:,jv)
       vbeta3sum(:) = vbeta3sum(:) + vbeta3(:,jv)
    ENDDO

    IF (printlev_loc>=4) WRITE (numout,*) '041113 vbeta3sum is: ', vbeta3sum   

    !! 1. Compute vevapsno (evaporation from snow) and vevapnu (bare soil evaporation)

    DO ji=1,kjpindex 

       !! \latexonly 
       !!     \input{enerbilevapveg1.tex}
       !! \endlatexonly
       speed = MAX(min_wind, SQRT (u(ji)*u(ji) + v(ji)*v(ji)))


       !! 1.1 Snow sublimation
       !! \latexonly 
       !!     \input{enerbilevapveg2.tex}
       !! \endlatexonly
       vevapsno(ji) = (un - vbeta5(ji)) * vbeta1(ji) * dt_sechiba * rau(ji) * speed * &
            & q_cdrag(ji) * (qsol_sat_new(ji) - qair(ji))

       !! 1.2 Bare soil evaporation
       !! \latexonly 
       !!     \input{enerbilevapveg3.tex}
       !! \endlatexonly
       vevapnu(ji) = (un - vbeta1(ji)) * (un-vbeta5(ji)) * vbeta4(ji) * dt_sechiba * &
            & rau(ji) * speed * q_cdrag(ji) * (qsol_sat_new(ji) - qair(ji))

       !
       ! 1.3 floodplains evaporation - transpiration et interception prioritaires dans les floodplains
       !
       vevapflo(ji) = vbeta5(ji) &
            & * dt_sechiba * rau(ji) * speed * q_cdrag(ji) * (qsol_sat_new(ji) - qair(ji))

    END DO

    !! 2. Compute transpir (transpiration) and vevapwet (interception)

    !! Preliminaries
    DO ji = 1, kjpindex

       !! \latexonly 
       !!     \input{enerbilevapveg4.tex}
       !! \endlatexonly
       speed = MAX(min_wind, SQRT (u(ji)*u(ji) + v(ji)*v(ji)))

       !! \latexonly 
       !!     \input{enerbilevapveg5.tex}
       !! \endlatexonly

       !! TO DO: To clean
       
       !+++TEMP+++
!!$        IF (printlev_loc>=4) THEN
!!$           WRITE(numout,*) 'Checking transpir'
!!$           WRITE(numout,*) 'vbeta1(ji), ',vbeta1(ji)
!!$           WRITE(numout,*) 'qsol_sat_new(ji)', qsol_sat_new(ji)
!!$           WRITE(numout,*) 'qair(ji), ', qair(ji)
!!$           WRITE(numout,*) 'rau(ji), ', rau(ji)
!!$           WRITE(numout,*) 'speed, ', speed
!!$           WRITE(numout,*) 'q_cdrag(ji)', q_cdrag(ji)
!!$           WRITE(numout,*) 'vbeta3(ji,jv)', vbeta3(ji,:)
!!$        ENDIF



       xx(ji) = dt_sechiba * (un-vbeta1(ji)) * (un-vbeta5(ji)) * (qsol_sat_new(ji)-qair(ji)) &
            & * rau(ji) * speed * q_cdrag(ji)
       !++++++++++

    ENDDO

    DO jv=1,nvm 
       DO ji=1,kjpindex 

          !! 2.1 Calculate interception loss
          !! \latexonly 
          !!     \input{enerbilevapveg6.tex}
          !! \endlatexonly
          vevapwet(ji,jv) = xx(ji) * vbeta2(ji,jv)
          ! 
          !! 2.2 Calculate transpiration
          !! \latexonly 
          !!     \input{enerbilevapveg7.tex}
          !! \endlatexonly


          !! TO DO: Think again about this as it is strange.

          !! Debug ID3.
          !! Transpiration should be the sum of all layers' one

          transpir(ji,jv) = transpir(ji,jv)

          !! End ID3

          !! TO DO: To clean.

          !IF (jv .eq. 7) THEN
          !     WRITE(numout,*) '261113b in first calc PFT4, vbeta3(ji,jv) is: ', vbeta3(ji,jv)
          !     WRITE(numout,*) '261113b in first calc PFT4, transpir (ji,jv) is: ', transpir (ji,jv)
          !     WRITE(numout,*) '261113b in first calc PFT4, (1.0_r_std - vbeta1(ji)) is: ', (1.0_r_std - vbeta1(ji))
          !     WRITE(numout,*) '261113b in first calc PFT4, qsol_sat_new(ji) is: ', qsol_sat_new(ji)
          !     WRITE(numout,*) '261113b in first calc PFT4, qair(ji) is: ', qair(ji)
          !     WRITE(numout,*) '261113b in first calc PFT4, ( rau(ji) * speed * q_cdrag(ji)) is: ', ( rau(ji) * speed * q_cdrag(ji))
          !END IF

          transpot(ji,jv) = xx(ji) * vbeta3pot(ji,jv)

       END DO
    END DO

    IF (printlev>=4) WRITE (numout,*) ' enerbil_evapveg done '

  END SUBROUTINE mleb_evapveg


  !!  =============================================================================================================================
  !! SUBROUTINE                 : mleb_flux
  !!
  !>\BRIEF                  Computes the new soil temperature, net radiation and
  !! latent and sensible heat flux for the new time step.
  !!
  !! DESCRIPTION                : This routine diagnoses, based on the new soil temperature, 
  !! the net radiation, the total evaporation, the latent heat flux, the sensible heat flux and the 
  !! sublimination flux. It also diagnoses the potential evaporation used for the fluxes (evapot) and the potential
  !! as defined by Penman & Monteith (Monteith, 1965) based on the correction term developed by Chris 
  !! Milly (1992). This Penman-Monteith formulation is required for the estimation of bare soil evaporation 
  !! when the 11 layer CWRR moisture scheme is used for the soil hydrology.
  !!
  !! MAIN OUTPUT VARIABLE(S)    : qsurf, fluxsens, fluxlat, fluxsubli, vevapp, lwup, lwnet,
  !! tsol_rad, netrad, evapot, evapot_corr
  !!
  !! REFERENCE(S)                   :
  !! - Best, MJ, Beljaars, A, Polcher, J & Viterbo, P, 2004. A proposed structure for coupling tiled
  !! surfaces with the planetary boundary layer. Journal of Hydrometeorology, 5, pp.1271-1278
  !! - de Noblet-Ducoudré, N, Laval, K & Perrier, A, 1993. SECHIBA, a new set of parameterisations
  !! of the hydrologic exchanges at the land-atmosphere interface within the LMD Atmospheric General
  !! Circulation Model. Journal of Climate, 6, pp.248-273
  !! - Dufresne, J-L & Ghattas, J, 2009. Description du schéma de la couche limite turbulente et la
  !! interface avec la surface planetaire dans LMDZ, Technical note, available (22/12/11):
  !! http://lmdz.lmd.jussieu.fr/developpeurs/notes-techniques/ressources/pbl_surface.pdf
  !! - Guimberteau, M, 2010. Modélisation de l'hydrologie continentale et influences de l'irrigation
  !! sur le cycle de l'eau, PhD Thesis, available (22/12/11):
  !! http://www.sisyphe.upmc.fr/~guimberteau/docs/manuscrit_these.pdf
  !! - Monteith, JL, 1965. Evaporation and Environment, paper presented at Symposium of the Society
  !! for Experimental Biology
  !! - Monteith & Unsworth, 2008. Principles of Environmental Physics (third edition), published Elsevier
  !! ISBN 978-0-12-505103-3
  !! - Milly, P. C. D., 1992: Potential Evaporation and Soil Moisture in General Circulation Models. 
  !! Journal of Climate, 5, pp. 209–226.
  !! - Polcher, J. McAvaney, B, Viterbo, P, Gaertner, MA, Hahmann, A, Mahfouf, J-F, Noilhan, J
  !! Phillips, TJ, Pitman, AJ, Schlosser, CA, Schulz, J-P, Timbal, B, Verseghy, D &
  !! Xue, Y, 1998. A proposal for a general interface between land surface schemes and
  !! general circulation models. Global and Planetary Change, 19, pp.261-276
  !! - Richtmeyer, RD, Morton, KW, 1967. Difference Methods for Initial-Value Problems.
  !! Interscience Publishers
  !! - Schulz, Jan-Peter, Lydia Dümenil, Jan Polcher, 2001: On the Land Surface–Atmosphere 
  !! Coupling and Its Impact in a Single-Column Atmospheric Model. J. Appl. Meteor., 40, 642–663.
  !!
  !! FLOWCHART                  : None
  !! \n
  !_ ==============================================================================================================================

  SUBROUTINE mleb_flux (kjpindex, emis, temp_sol, rau, u, v, q_cdrag, vbeta, vbeta1, vbeta5, &
       & qair, epot_air, psnew, qsurf, fluxsens, fluxlat, fluxsubli, vevapp, temp_sol_new, &
       & lwdown, swnet, lwup, lwnet, pb, tsol_rad, netrad, evapot, evapot_corr, &
       & precip_rain, snowdz,temp_air, pgflux, soilcap, temp_sol_add, qsol_sat_new)

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                               :: kjpindex      !! Domain size (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: emis          !! Emissivity (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: temp_sol      !! Surface temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: rau           !! Density (kg m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: u,v           !! Wind velocity in components u and v (m s^{-1}) 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: q_cdrag       !! Surface drag coefficient (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: vbeta         !! Resistance coefficient (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: vbeta1        !! Snow resistance  (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: vbeta5        !! Flood resistance 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: qair          !! Lowest level specific humidity (kg kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: epot_air      !! Air potential energy (J)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: psnew         !! New surface static energy (J kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: temp_sol_new  !! New surface temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: pb            !! Lowest level pressure (hPa)
    REAL(r_std),DIMENSION (kjpindex,nsnow), INTENT(in)       :: snowdz        !! Snow depth [m]
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: precip_rain   !! Rainfall
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: temp_air      !! Air temperature
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: lwdown        !! Downward long wave radiation (W m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: swnet         !! Net short wave radiation (W m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: soilcap       !! Soil calorific capacity including snow and soil (J K^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: fluxsens      !! Sensible heat flux (W m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: fluxlat       !! Latent heat flux (W m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: qsol_sat_new         !! Air saturated humidity (kg.kg-1)

    !! 0.2 Output variables

    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: qsurf         !! Surface specific humidity (kg kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: fluxsubli     !! Energy of sublimation (W m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: vevapp        !! Total of evaporation (mm day^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: lwup          !! Long-wave up radiation (W m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: lwnet         !! Long-wave net radiation (W m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: tsol_rad      !! Radiative surface temperature (W m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)         :: netrad        !! Net radiation (W m^{-2})
    REAL(r_std), DIMENSION (kjpindex),INTENT(out)            :: temp_sol_add  !! Additional energy to melt snow for snow ablation case (K)

    !! 0.3 Modified variables
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)         :: evapot        !! Soil Potential Evaporation (mm/tstep)
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)         :: evapot_corr   !! Soil Potential Evaporation Correction (mm/tstep)
    REAL(r_std), DIMENSION (kjpindex),INTENT(inout)          :: pgflux        !! Net energy into snowpack(W/m^2)

    !! 0.4 Local variables
    INTEGER(i_std)                                           :: ji
    REAL(r_std),DIMENSION (kjpindex)                         :: grad_qsat
    REAL(r_std)                                              :: correction
    REAL(r_std)                                              :: speed              !! Speed (m s^{-1}), 
    REAL(r_std)                                              :: qc          !! Surface drag coefficient (??)
    REAL(r_std), DIMENSION (kjpindex)                        :: zgflux
    REAL(r_std), DIMENSION (kjpindex)                        :: qsol_sat_tmp
    REAL(r_std), DIMENSION (kjpindex)                        :: zerocelcius
    REAL(r_std), DIMENSION (kjpindex)                        :: PHPSNOW
    REAL(r_std), DIMENSION (kjpindex)                        :: lwup_tmp
    REAL(r_std), DIMENSION (kjpindex)                        :: netrad_tmp
    REAL(r_std), DIMENSION (kjpindex)                        :: fluxsens_tmp
    REAL(r_std), DIMENSION (kjpindex)                        :: fluxlat_tmp
    LOGICAL,DIMENSION (kjpindex)                             :: warning_correction
    !_ ================================================================================================================================


    
    !! This routine calculates the fluxes emanating from the pixel (the forest) based on the new surface temperature.
    !! This routine is more or less the same as the one in enerbil.f90 (except for the latent and sensible heat fluxes, 
    !! calculated thanks to the transport inside the canopy.

    !! TO DO: As we consider the module as applicable only for one pixel, we need to see if this routine has to be reviewed
    !!        or not!


    zerocelcius(:) = tp_00
    CALL qsatcalc (kjpindex, zerocelcius, pb, qsol_sat_tmp)

    DO ji=1,kjpindex

       !! 1. Determination of 'housekeeping' variables

       !! The horizontal wind speed is calculated from the velocities along 
       !! each axis using Pythagorus' theorem. 
       !! \latexonly 
       !!     \input{enerbilflux1.tex}
       !! \endlatexonly
       speed = MAX(min_wind, SQRT (u(ji)*u(ji) + v(ji)*v(ji)))

       !! From definition, the surface drag coefficient is defined:
       !! \latexonly 
       !!     \input{enerbilflux2.tex}
       !! \endlatexonly
       qc = speed * q_cdrag(ji)

       !! 2. Calculation of the net upward flux of longwave radiation
       !! We first of all calculate the radiation as a result of the Stefan-Boltzmann equation,
       !! which is the sum of the calculated values at the surface temperature  at the 'old' 
       !! temperature and the value that corresponds to the difference between the old temperature
       !! and the temperature at the 'new' timestep.
       !! \latexonly 
       !!     \input{enerbilflux3.tex}
       !! \endlatexonly


       lwup(ji) = MAX(0.,emis(ji) * c_stefan * temp_sol(ji)**4 + &
            &     quatre * emis(ji) * c_stefan * temp_sol(ji)**3 * &
            &     (temp_sol_new(ji) - temp_sol(ji)) )

       !! We then add to this value that from the reflected longwave radiation:
       !! \latexonly 
       !!     \input{enerbilflux4.tex}
       !! \endlatexonly
       lwup(ji) = lwup(ji)  +  (un - emis(ji)) * lwdown(ji)

       !! The radiative surface temperature is calculated according to the Stefan-Boltzmann relation:
       !! \latexonly 
       !!     \input{enerbilflux5.tex}
       !! \endlatexonly      
       ! tsol_rad(ji) = emis(ji) * c_stefan * temp_sol(ji)**4 + lwup(ji)

       !! The implicit solution computes an emitted long-wave flux which is a limited Taylor expansion 
       !! of the future surface temperature around the current values. Thus the long-wave flux does not  
       !! correspond to the new surface temperature but some intermediate value. So we need to deduce the 
       !! radiative surface temperature to which this flux corresponds. 
       !! 
       tsol_rad(ji) = (lwup(ji)/ (emis(ji) * c_stefan)) **(1./quatre) 

       !! qsurf (the surface specific humidity) is a simple diagnostic which will be used by the 
       !! GCM to compute the dependence of of the surface layer stability on moisture.
       !! \latexonly 
       !!     \input{enerbilflux6.tex}
       !! \endlatexonly
       qsurf(ji) = (vbeta1(ji) * (un - vbeta5(ji)) + vbeta5(ji)) * qsol_sat_new(ji) + &
            & (un - vbeta1(ji))*(un - vbeta5(ji)) * vbeta(ji) * qsol_sat_new(ji)

       qsurf(ji)=qsurf(ji) + (1 - (vbeta1(ji) * (un - vbeta5(ji)) + vbeta5(ji)))*qair(ji) + &                          
            & (1- (un - vbeta1(ji))*(un - vbeta5(ji)) * vbeta(ji) )*qair(ji)

       ! When no turbulence and no wind (qc small) the near surface air moisture is equal to qair.
       ! This is related to problem with interpolation in stdlevvar in LMDZ. 
       IF ( qc .LE. min_qc ) qsurf(ji) = qair(ji)

       IF (ji .eq. 1 .AND. (printlev_loc>=4)) THEN
          WRITE(numout,*) "qsurf(1) is: ", qsurf(1) 
          WRITE(numout,*) "qair(1) is: ", qair(1) 
       END IF

       !! \latexonly 
       !!     \input{enerbilflux7.tex}
       !! \endlatexonly
       qsurf(ji) = MAX(qsurf(ji), qair(ji))

       !! Net downward radiation is the sum of the down-welling less the up-welling long wave flux, plus
       !! the short wave radiation.
       !! \latexonly 
       !!     \emissivity absorbed lw radiationinput{enerbilflux8.tex}
       !! \endlatexonly
       netrad(ji) = lwdown(ji) + swnet(ji) - lwup(ji) 

       !! 'vevapp' is the sum of the total evaporative processes (snow plus non-snow processes).
       !! \latexonly 
       !!     \input{enerbilflux9.tex}
       !! \endlatexonly
       vevapp(ji) = dt_sechiba * rau(ji) * qc * (vbeta1(ji) * (un - vbeta5(ji)) + vbeta5(ji)) * &
            & (qsol_sat_new(ji) - qair(ji)) + &
            &  dt_sechiba * rau(ji) * qc * (un - vbeta1(ji))*(un-vbeta5(ji)) * vbeta(ji) * &
            & (qsol_sat_new(ji) - qair(ji))

       !! The total latent heat flux is the sum of the snow plus non-snow processes.
       !! \latexonly 
       !!     \input{enerbilflux10.tex}
       !! \endlatexonly
       ! commented out CN-CAN because it is calculate in mleb_calc_column but it
       ! now is an input variable to calculate the energy in melting the snowpack 
!!$       fluxlat(ji) = chalsu0 * rau(ji) * qc * vbeta1(ji) * (un - vbeta5(ji)) * &
!!$             & (qsol_sat_new(ji) - qair(ji)) + &
!!$             &  chalev0 * rau(ji) * qc * vbeta5(ji) *&
!!$             & (qsol_sat_new(ji) - qair(ji)) + &
!!$             &  chalev0 * rau(ji) * qc * (un - vbeta1(ji)) * (un - vbeta5(ji)) * vbeta(ji) * &
!!$             & (valpha(ji) * qsol_sat_new(ji) - qair(ji))

       !! The sublimination flux concerns is calculated using vbeta1, the snow resistance.
       !! \latexonly 
       !!     \input{enerbilflux11.tex}
       !! \endlatexonly
       fluxsubli(ji) = chalsu0 * rau(ji) * qc * vbeta1(ji) * (un - vbeta5(ji)) * &
            & (qsol_sat_new(ji) - qair(ji)) 

       !! The sensible heat flux is a factor of the difference between the new surface static energy
       !! and the potential energy of air.
       !! \latexonly 
       !!     \input{enerbilflux12.tex}
       !! \endlatexonly
       ! commented out CN-CAN because it is calculate in mleb_calc_column but it
       ! now is an input variable to calculate the energy in melting the snowpack 
!!$       fluxsens(ji) =  rau(ji) * qc * (psnew(ji) - epot_air(ji))

       !! This is the net longwave downwards radiation.
       !! \latexonly 
       !!     \input{enerbilflux13.tex}
       !! \endlatexonly
       lwnet(ji) = lwdown(ji) - lwup(ji)

       !! Diagnoses the potential evaporation used for the fluxes (evapot)
       !! \latexonly 
       !!     \input{enerbilflux14.tex}
       !! \endlatexonly  
       evapot(ji) =  MAX(zero, dt_sechiba * rau(ji) * qc * (qsol_sat_new(ji) - qair(ji)))

       !! From definition we can say:
       !! \latexonly 
       !!     \input{enerbilflux15.tex}
       !! \endlatexonly 
       IF(printlev_loc >= 4) THEN
          WRITE(numout,*)'071116 mleb_flux tair before',tair(ji)
          WRITE(numout,*)'071116 mleb_flux epot_air before',epot_air(ji)
       ENDIF
       tair(ji)  =  epot_air(ji) / cp_air
       IF(printlev_loc >= 4) WRITE(numout,*)'071116 mleb_flux tair after',tair(ji)

       !! To calculate net energy flux into the snowpack
       PHPSNOW(ji) = precip_rain(ji)*(4.218E+3)*(MAX(tp_00,temp_air(ji))-tp_00)/dt_sechiba ! (w/m2)
       pgflux(ji)  = netrad(ji) - fluxsens(ji) - fluxlat(ji) + PHPSNOW(ji)

       !! To get the extra energy used to melt the snowpack
       IF (temp_sol_new (ji) > tp_00 .AND. &
            SUM(snowdz(ji,:)) .GT. zero .AND. soilcap(ji) .GT. min_sechiba) THEN

          lwup_tmp(ji) = emis(ji) * c_stefan * temp_sol(ji)**4 + &
               quatre * emis(ji) * c_stefan * temp_sol(ji)**3 * (tp_00 - temp_sol(ji))

          lwup_tmp(ji) = lwup_tmp(ji)  +  (un - emis(ji)) * lwdown(ji)
          netrad_tmp(ji) = lwdown(ji) + swnet(ji) - lwup_tmp(ji)
          fluxsens_tmp(ji) =  rau(ji) * qc * cp_air * (tp_00 - epot_air(ji)/cp_air)
          fluxlat_tmp(ji) = chalsu0 * rau(ji) * qc * vbeta1(ji) * (un-vbeta5(ji)) * &
               & (qsol_sat_tmp(ji) - qair(ji)) + &
               &  chalev0 * rau(ji) * qc * vbeta5(ji) *&
               & (qsol_sat_tmp(ji) - qair(ji)) + &
               &  chalev0 * rau(ji) * qc * (un - vbeta1(ji)) * (un-vbeta5(ji))* vbeta(ji) * &
               & (qsol_sat_tmp(ji) - qair(ji))

          zgflux(ji)  = netrad_tmp(ji) - fluxsens_tmp(ji) - fluxlat_tmp(ji)+PHPSNOW(ji)

          temp_sol_add(ji) = -(pgflux(ji) - zgflux(ji))*dt_sechiba/soilcap(ji)

          pgflux(ji) = zgflux(ji)

       ELSE

          temp_sol_add(ji) = zero

       ENDIF

    ENDDO ! DO ji=1,kjpindex

    !! 3. Define qsat_air with the subroutine src_parameter:

    ! WRITE(*,*) '181113a test, pre call to qsatcalc enerbil line 6210; tair is: ', tair
    CALL qsatcalc(kjpindex, tair, pb, qsat_air)
    CALL dev_qsatcalc(kjpindex, tair, pb, grad_qsat)

    ! grad_qsat(:)= (qsol_sat_new(:)- qsat_air(:)) / ((psnew(:) - epot_air(:)) / cp_air) ! * dt_sechiba

    warning_correction(:)=.FALSE.
    DO ji=1,kjpindex

       !! \latexonly 
       !!     \input{enerbilflux16.tex}
       !! \endlatexonly 
       speed = MAX(min_wind, SQRT (u(ji)*u(ji) + v(ji)*v(ji)))

       !! \latexonly 
       !!     \input{enerbilflux17.tex}
       !! \endlatexonly 
       qc = speed * q_cdrag(ji)

       !! Derive the potential as defined by Penman & Monteith (Monteith, 1965) based on the correction 
       !! term developed by Chris Milly (1992).
       IF ((evapot(ji) .GT. min_sechiba) .AND. ((psnew(ji) - epot_air(ji)) .NE. zero )) THEN

          !! \latexonly 
          !!     \input{enerbilflux18.tex}
          !! \endlatexonly 
          correction =  (quatre * emis(ji) * c_stefan * tair(ji)**3 + rau(ji) * qc * cp_air + &
               &                  chalev0 * rau(ji) * qc * grad_qsat(ji) * vevapp(ji) / evapot(ji) )

          !! \latexonly 
          !!     \input{enerbilflux19.tex}
          !! \endlatexonly 
          IF (ABS(correction) .GT. min_sechiba) THEN
             correction = chalev0 * rau(ji) * qc * grad_qsat(ji) * (un - vevapp(ji)/evapot(ji)) / correction
          ELSE
             warning_correction(ji)=.TRUE.
          ENDIF

       ELSE

          correction = zero

       ENDIF

       correction = MAX (zero, correction)

       !! \latexonly 
       !!     \input{enerbilflux20.tex}
       !! \endlatexonly 
       evapot_corr(ji) = evapot(ji) / (un + correction)

    ENDDO

    IF ( ANY(warning_correction) ) THEN
       DO ji=1,kjpindex
          IF ( warning_correction(ji) ) THEN
             WRITE(numout,*) ji,"Denumerator of the Milly Correction is zero, so no correction has been applied"
          ENDIF
       ENDDO
    ENDIF

    IF (printlev>=4) WRITE (numout,*) ' enerbil_flux done '

  END SUBROUTINE mleb_flux


  !!  ===========================================================================================================================
  !! SUBROUTINE                 : mleb_write
  !!
  !!
  !>\BRIEF                      Write a bunch of variables to the history file
  !!
  !! DESCRIPTION                : 
  !! Originally the routine enerbil_main was used for this, but that seemed to create
  !! problems with uninitialized variables when it was called after mleb_main.
  !! \n
  !!
  !! RECENT CHANGE(S): None
  !!
  !! MAIN OUTPUT VARIABLE(S)    : None
  !!
  !! REFERENCE(S)       :
  !!
  !! FLOWCHART          :
  !! \n
  !_ ==============================================================================================================================

  SUBROUTINE mleb_write (kjit, kjpindex, &
       & index, lwdown, temp_sol, &
!!$       & index, lwdown, ccanopy, temp_sol, &
       & temp_sol_new, evapot, evapot_corr, hist_id, hist2_id, &
       & vevapp, vevapwet, transpir, vevapnu, vevapsno, vevapflo)

    !! 0 Variable and parameter description

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                           :: kjit          !! Time step number (-)
    INTEGER(i_std), INTENT(in)                           :: kjpindex      !! Domain size (-)
    INTEGER(i_std),INTENT (in)                           :: hist_id       !! _Restart_ file and _history_ file identifier (-)
    INTEGER(i_std),INTENT (in)                           :: hist2_id      !! _history_ file 2 identifier (-)

    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)     :: index         !! Indeces of the points on the map (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: lwdown        !! Down-welling long-wave flux (W m^{-2})
!!$    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: ccanopy       !! CO_2 concentration in the canopy (ppm)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: evapot        !! Soil Potential Evaporation (mm/tstep)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: evapot_corr   !! Soil Potential Evaporation Correction (mm/tstep)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: temp_sol      !! Surface temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: temp_sol_new  !! New surface temperature (K)
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)        :: vevapp        !! Total of evaporation (mm day^{-1})
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(in)    :: vevapwet      !! Interception (mm day^{-1})
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(in)    :: transpir      !! Transpiration (mm day^{-1})
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)        :: vevapnu       !! Bare soil evaporation (mm day^{-1})
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)        :: vevapsno      !! Snow evaporation (mm day^{-1})               
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)        :: vevapflo      !! Floodplains evaporation      


    !! 0.2 Output variables

    !! 0.3 Modified variables


    !! 0.4 Local variables
    REAL(r_std),DIMENSION (kjpindex)             :: diffevap         !! Difference betwence vevapp and composing fluxes (Kg/m^2/s)
    INTEGER(i_std)                               :: ji               !! Local index    


    !_ ================================================================================================================================



    !CALL xios_orchidee_send_field("netrad",netrad)
    CALL xios_orchidee_send_field("evapot",evapot/dt_sechiba)
    CALL xios_orchidee_send_field("evapot_corr",evapot_corr/dt_sechiba)
    CALL xios_orchidee_send_field("lwdown",lwabs)
    CALL xios_orchidee_send_field("lwnet",lwnet)
    CALL xios_orchidee_send_field("Qv",fluxsubli)

    DO ji=1,kjpindex 
       diffevap(ji) = vevapp(ji) - ( SUM(vevapwet(ji,:)) + & 
            SUM(transpir(ji,:)) + vevapnu(ji) + vevapsno(ji) + vevapflo(ji) )  
    ENDDO
    CALL xios_orchidee_send_field("diffevap",diffevap/dt_sechiba) ! mm/s 


    IF ( .NOT. almaoutput ) THEN
       CALL histwrite_p(hist_id, 'evapot', kjit, evapot, kjpindex, index)
       CALL histwrite_p(hist_id, 'evapot_corr', kjit, evapot_corr, kjpindex, index)
       CALL histwrite_p(hist_id, 'lwdown', kjit, lwabs,  kjpindex, index)
       CALL histwrite_p(hist_id, 'lwnet',  kjit, lwnet,  kjpindex, index)
       IF ( hist2_id > 0 ) THEN
          CALL histwrite_p(hist2_id, 'evapot', kjit, evapot, kjpindex, index)
          CALL histwrite_p(hist2_id, 'evapot_corr', kjit, evapot_corr, kjpindex, index)
          CALL histwrite_p(hist2_id, 'lwdown', kjit, lwabs,  kjpindex, index)
          CALL histwrite_p(hist2_id, 'lwnet',  kjit, lwnet,  kjpindex, index)
       ENDIF
    ELSE
       CALL histwrite_p(hist_id, 'LWnet', kjit, lwnet, kjpindex, index)
       CALL histwrite_p(hist_id, 'Qv', kjit, fluxsubli, kjpindex, index)
       CALL histwrite_p(hist_id, 'PotEvap', kjit, evapot_corr, kjpindex, index)
       CALL histwrite_p(hist_id, 'PotEvapOld', kjit, evapot, kjpindex, index)
       IF ( hist2_id > 0 ) THEN
          CALL histwrite_p(hist2_id, 'LWnet', kjit, lwnet, kjpindex, index)
          CALL histwrite_p(hist2_id, 'Qv', kjit, fluxsubli, kjpindex, index)
          CALL histwrite_p(hist2_id, 'PotEvap', kjit, evapot_corr, kjpindex, index)
       ENDIF
    ENDIF

    IF (printlev>=4) WRITE (numout,*) ' enerbil_main Done '

  END SUBROUTINE mleb_write


  !!
  !===========================================================================================================================
  !! SUBROUTINE: mleb_boxheight
  !!
  !!
  !>\BRIEF: calculates box heights and plant area density 
  !!
  !! DESCRIPTION: The subroutine calculates box heights and plant area density as
  !! a function of height. It was made to move a lot of code from the begining of
  !! mleb_main (which is still there), so needs to veryfy that this subroutine
  !! gives identical results.  
  !! 
  !!
  !! RECENT CHANGE(S): None
  !!
  !! TO BE DONE: Verify that this subroutine gives in identical results as the
  !! code currently in the begining of mleb main. Delete the call to
  !! mleb_reallocate_column, and remove the variables temp_atmos_pres_grid,
  !! q_atmos_pres_grid,temp_leaf_pres_grid  
  !!
  !! MAIN OUTPUT VARIABLE(S): box_height, pad_deltah 
  !!
  !! REFERENCE(S):
  !!
  !! FLOWCHART:
  !_
  !==============================================================================================================================


  SUBROUTINE  mleb_boxheight(kjpindex,laieff_isotrop, z_array_out, max_height_store, &
       & rau,temp_atmos_pres_grid, q_atmos_pres_grid,temp_leaf_pres_grid, &
       & box_height,top_lev, pad_deltah, jrpad)

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                                            :: kjpindex                !! Domain size (-)
    REAL(r_std), DIMENSION(kjpindex,nlevels_tot,nvm), INTENT(in)          :: laieff_isotrop          !! Effective LAI
    REAL(r_std),DIMENSION(kjpindex,nvm,ncirc,nlevels_tot), INTENT(in)     :: z_array_out             !! An output of h_array, to use in sechiba 
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(in)                     :: max_height_store
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)                        :: rau                     !! Air density (kg m^{-3})


    !! 0.2 Output variables
    REAL(r_std), DIMENSION(jnlvls+1), INTENT(out)            :: top_lev
    REAL(r_std), DIMENSION (0:jnlvls), INTENT(out)            :: jrpad                  !! PlantArea Density, replacement for values above 
    !!(as definitions between different)
    REAL(r_std), DIMENSION (kjpindex, jnlvls), INTENT(in)     :: q_atmos_pres_grid


    !! 0.3 Modified variables
    REAL(r_std), DIMENSION(jnlvls+1), INTENT(out)               :: box_height    !!height of each box edge, starting at the surface layer (m)
    REAL(r_std), DIMENSION(0:jnlvls), INTENT(out)               :: pad_deltah    !! Plant Area Density as a function 
    !! of height above the ground (m^2/m^3)
!!$    REAL(r_std), DIMENSION(jnlvls+1), INTENT(out)               :: pad_deltah    !! Plant Area Density as a function 
!!$                                                                                 !! of height above the ground (m^2/m^3)
    REAL(r_std), DIMENSION (kjpindex, jnlvls), INTENT(in)      :: temp_leaf_pres_grid
    REAL(r_std), DIMENSION (kjpindex, jnlvls), INTENT(in)      :: temp_atmos_pres_grid


    !! 0.4 Local variables
    INTEGER(i_std)                               :: iv,i,ji,kjit,jv ! fix the indexing
    INTEGER(i_std)                               :: ij    

    REAL(r_std)                                  :: pft_running_total
    REAL(r_std), DIMENSION(0:jnlvls)             :: pad_deltah_running_total
    REAL(r_std), DIMENSION(jnlvls+1)             :: box_height_running_total
    REAL(r_std), DIMENSION(jnlvls+1)             :: box_height_old

    ! ------------------------------------------------------------------------------------
    !
    ! there are two ways to arrive at the MLA (mass per unit leaf area, kg/m^2), 
    ! which we then multiply by the LAD (leaf area density, m^2/m^3) to arrive
    ! at the 
    ! vegetation matter per m in height

    ! 1) as in Nobel, Park S. (2005) 'Physiochemical and Environmental Plant
    ! Physiology' 
    ! Edition 3, published Elsevier, chapter 7.1, p. 309
    !    - assume leaf has a specific heat of water (4.19 kJ kg^-1 degC^-1 at 20
    !    degC)
    !    - assume leaf has thickness 300micrometres
    !    - assume leaf has density of 700kg/m^3 (water is 1000kg/m^3 but leaves
    !    are often 
    !        30% air)
    !    - this gives (300 x 10^-6 m)(700 kg m^-3) = 0.21 kg m^-2

    ! 2) Using LMA (leaf per unit mass area, m^2/kg) from the TRY database
    !    - for Eucalyptus delegatensis we have 155 cm^2/g ---> 0.155 m^2/g. 
    !      This is dried weight of the leaf.
    !    - 100 gC/m^2 soil
    !    - 100 * 2 * 10 = 2000 g C / m^2 soil
    !    Dry weight / Fresh weight = 8 / SLA ---> g FW / m^2 (here we are
    !    assuming a dry 
    !                                               weight to fresh weight ratio
    !                                               of 1:9)
    !                              = 8 / 0.155
    !                              = 516 g FW / m^2
    !                              = 0.516 kg / m^2
    ! 
    ! ------------------------------------------------------------------------------------

    ! asla WRITE(numout, *) 'kjpindex', kjpindex
    ! asla  WRITE(numout, *) 'rau',rau
    ! asla WRITE(numout, *) 'jtheta',jtheta
    ! asla WRITE(numout, *) 'rho_veg', rho_veg
    ! asla WRITE(numout, *) 'leaf_tks',leaf_tks
    ! asla WRITE(numout, *) 'jnlvls', jnlvls
    IF (printlev_loc>=4) THEN
       WRITE (numout, *) 'SUM(laieff_isotrop(test_grid,1:jnlvls_canopy,test_pft)) is: ', &
            & SUM(laieff_isotrop(test_grid,1:jnlvls_canopy,test_pft))
    END IF !(printlev_loc>=4)

    IF (jnlvls .ne. 1) THEN
       IF (jnlvls .ne. (jnlvls_under + jnlvls_canopy + jnlvls_over)) THEN
          WRITE (numout, *) 'number of levels in canopy column imbalance'
          STOP
       END IF
    END IF

    IF (.NOT. hack_can_structure) THEN

       pft_running_total = 0.0d0
       box_height_running_total = 0.0d0
       pad_deltah_running_total = 0.0d0

       DO ji = 1,1
          ! Currently, the mleb only runs for one PFT, and therefore
          ! jv is set to the PFT that is used instead of 1,nvm. 
          ! DO jv = 1,nvm 
          DO jv = test_pft,test_pft

             ! IF (veget_max(ji,jv) .GT. zero .AND. veget(ji,jv) .GT. zero)
             ! THEN 

             ! IF ((SUM(laieff_isotrop(ji, 1:jnlvls_canopy, jv)) .gt.
             ! min_sechiba) .AND. &
             !      & (SUM(z_array_out(ji, jv, 1, :)) .gt. min_sechiba))
             !      T/HEN

             !IF (SUM(z_array_out(ji, jv, 1, :)) .gt. min_sechiba) THEN

             ! update plant area density
             pad_deltah(:) = 0.001d0
             ! n.b. the first level of lai_eff_isotrop represents the
             ! understorey.
             ! We set this equal to the min!
             DO i = (jnlvls_under + 1) , (jnlvls_under + jnlvls_canopy)
                pad_deltah(i) = laieff_isotrop(ji, i-jnlvls_under+1, jv) * cos(pi/3.0d0)
                IF (pad_deltah(i) .lt. 1.0d-3) THEN
                   pad_deltah(i) = 1.0d-3
                END IF
                IF (printlev_loc>=4) THEN
                   WRITE (numout, *) '251016 laieff_isotrop:',ji, i-jnlvls_under+1, jv, 'is', &
                        laieff_isotrop(ji, i-jnlvls_under, jv)
                   WRITE (numout, *) '160315 pad_deltah level: ',i,' is: ', pad_deltah(i)
                END IF ! (printlev_loc>=4)
             END DO ! i = (jnlvls_under + 1) , (jnlvls_under +jnlvls_canopy)

             IF (printlev_loc>=4) THEN
                WRITE (numout, *) 'in calculate box_heightsection,pft',jv
                WRITE (numout, *) 'max_height_store(ji,jv) is: ',max_height_store(ji, jv)
                WRITE (numout, *) 'z_array_out(ji,jv,1,:) is: ',z_array_out(ji, jv, 1, :)
             END IF ! (printlev_loc>=4) THEN

             DO i = 1, jnlvls_under - 1
                box_height(i) = (i) * (z_array_out(ji, jv, 1, 2) -z_array_out(ji, jv, 1, 1)) / (jnlvls_under)
             END DO ! 1, jnlvls_under - 1 

             box_height(jnlvls_under) = z_array_out(ji, jv, 1, 2)



             DO i = (jnlvls_under + 1), (jnlvls_under + jnlvls_canopy- 1)
                box_height(i) = z_array_out(ji, jv, 1, i -jnlvls_under + 2)
             END DO ! i = (jnlvls_under + 1), (jnlvls_under +jnlvls_canopy - 1)


             box_height(jnlvls_under + jnlvls_canopy) = max_height_store(ji, jv)

             DO i = (jnlvls_under + jnlvls_canopy + 1), (jnlvls_under+ jnlvls_canopy + jnlvls_over - 1)
                box_height(i) = max_height_store(ji, jv) + &
                     & (( i - (jnlvls_under +jnlvls_canopy) ) * &
                     & (49.0d0 - max_height_store(ji, jv))/ jnlvls_over )
             END DO ! i = (jnlvls_under + jnlvls_canopy + 1),(jnlvls_under + jnlvls_canopy + jnlvls_over)


             box_height(jnlvls_under + jnlvls_canopy + jnlvls_over) =49.0d0
             box_height(jnlvls_under + jnlvls_canopy + jnlvls_over+1)= 50.0d0

             jrpad(:) = pad_deltah(:)

             IF (printlev_loc>=4) THEN
                DO i = 1, jnlvls
                   WRITE(numout, *) 'lev: ', i, ' b_h: ',box_height(i), ' p_d: ', pad_deltah(i)
                END DO ! i = 1, jnlvls
             END IF ! (printlev_loc>=4) THEN

             ! END IF ! (SUM(laieff_isotrop(ji, 1:jnlvls_canopy, jv)) .ne.
             ! 0.0)



             ! IF (SUM(laieff_isotrop(ji, 1:jnlvls_canopy, jv)) .gt.
             ! min_sechiba) THEN
             !  pft_running_total = pft_running_total + 1
             ! END IF

             !  box_height_running_total = box_height_running_total +
             !  box_height ! NOT needed for one PFT
             !  pad_deltah_running_total = pad_deltah_running_total +
             !  pad_deltah ! NOT needed for one PFT 

             ! END IF !(veget_max(ji,jv) .GT. zero .AND. veget(ji,jv) .GT.
             ! zero) THEN 

          END DO ! jv = 2, nvm
          IF (printlev_loc>=4) THEN
             WRITE(numout, *) '100515 pft_running_total is: ',pft_running_total
             WRITE(numout, *) '100515 box_height_running_total(:) is: ',box_height_running_total(:)
             WRITE(numout, *) '100515 pad_deltah_running_total(:) is: ',pad_deltah_running_total(:)
          END IF ! (printlev_loc>=4) THEN

          ! IF (pft_running_total .gt. min_sechiba) THEN
          !    box_height = box_height_running_total / pft_running_total
          !    pad_deltah = pad_deltah_running_total / pft_running_total
          ! ELSE
          !    box_height = 0.0d0
          !    pad_deltah = 0.0d0
          ! ENDIF

          !  box_height = box_height_running_total / 1.0d0  ! NOT needed for
          !  one PFT
          !  pad_deltah = pad_deltah_running_total / 1.0d0  ! NOT needed for
          !  one PFT


          box_height_old = box_height


          IF(printlev_loc >= 4)THEN
             WRITE(numout,*) 'box_height_old', box_height_old(:)
             WRITE(numout,*) 'box_height', box_height(:)
             WRITE(numout,*) 'temp_atmos_pres_grid',temp_atmos_pres_grid(ji,:)
             WRITE(numout,*) 'q_atmos_pres_grid',q_atmos_pres_grid(ji,:) 
             WRITE(numout,*) 'temp_leaf_pres_grid',temp_leaf_pres_grid(ji,:)
          ENDIF ! (printlev_loc >= 4)
          box_height_old = box_height


       END DO ! ji = 1, kjpindex
    ELSE
       ! in cases hack_can_structure is true mleb does not work


    END IF ! (hack_can_structure) THEN 

    top_lev(:) = box_height(:)


  END SUBROUTINE mleb_boxheight

  ! --------------------------------------------------------------------------------------------------------

  !!
  !=============================================================================================================================
  !! SUBROUTINE:               mleb_finalize
  !!
  !>\BRIEF                     Write to restart file
  !!
  !! DESCRIPTION:              This subroutine writes the module variables and
  !!                           variables calculated in mleb to the restart file.
  !!
  !! RECENT CHANGE(S): None
  !!
  !! REFERENCE(S): None
  !! 
  !! FLOWCHART: None
  !! \n
  !_
  !==============================================================================================================================
  SUBROUTINE mleb_finalize (kjit,    kjpindex,    rest_id,                &
       evapot,  evapot_corr, temp_sol, tsol_rad,     &
       qsurf,   fluxsens,    fluxlat,  vevapp,       &
       u_speed,  z_array_out,  &
       max_height_store )

    IMPLICIT NONE

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                            :: kjit             !! Time step number (unitless) 
    INTEGER(i_std), INTENT(in)                            ::kjpindex                !! Domain size (-)
    INTEGER(i_std),INTENT (in)                            :: rest_id          !! Restart_ file and history file identifier 


    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: evapot              !! Soil Potential Evaporation (mm/tstep)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: evapot_corr         !! Soil Potential Evaporation Correction (mm/tstep)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: temp_sol            !! Surface temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: tsol_rad            !! Tsol_rad (W m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: qsurf               !! Surface specific humidity (kg kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: fluxsens            !! Sensible heat flux (W m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: fluxlat             !! Latent heat flux (W m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: vevapp              !! Total of evaporation (mm day^{-1})
    REAL(r_std), DIMENSION(kjpindex,jnlvls),INTENT(in)          :: u_speed             !! Canopy wind speed profile (the dimension of u_speed only make sense, 
    !! when running site scale simulations, thus .NOT. hack_can_structure
    !! will only work at the site scale??
    REAL(r_std), DIMENSION(kjpindex,nvm,ncirc,nlevels_tot), INTENT(in):: z_array_out  !! heights of tree levels from stomate
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(in)     :: max_height_store


    !! 0.4 Local variables
    INTEGER(i_std)                                     :: ier
    INTEGER(i_std)                                     :: ji, ilev
    CHARACTER(LEN=80) :: var_name



    !_ ================================================================================================================================

    !! 1. Write variables to restart file to be used for the next simulation
    IF (printlev_loc>=3) WRITE (numout,*) 'Write restart file with mleb variables'

    CALL restput_p(rest_id, 'temp_sol', nbp_glo, 1, 1, kjit,  temp_sol, 'scatter',  nbp_glo, index_g)
    CALL restput_p(rest_id, 'qsurf', nbp_glo, 1, 1, kjit,  qsurf, 'scatter',  nbp_glo, index_g)
    CALL restput_p(rest_id, 'evapot', nbp_glo, 1, 1, kjit,  evapot, 'scatter',  nbp_glo, index_g)
    CALL restput_p(rest_id, 'evapot_corr', nbp_glo, 1, 1, kjit,  evapot_corr, 'scatter',  nbp_glo, index_g)
    CALL restput_p(rest_id, 'tsolrad', nbp_glo, 1, 1, kjit,  tsol_rad, 'scatter',  nbp_glo, index_g)
    CALL restput_p(rest_id, 'evapora', nbp_glo, 1, 1, kjit,  vevapp, 'scatter',  nbp_glo, index_g)
    CALL restput_p(rest_id, 'fluxlat', nbp_glo, 1, 1, kjit,  fluxlat, 'scatter',  nbp_glo, index_g)
    CALL restput_p(rest_id, 'fluxsens', nbp_glo, 1, 1, kjit,  fluxsens, 'scatter',  nbp_glo, index_g)


    CALL restput_p (rest_id, 'temp_atmos_pres_grid', nbp_glo, jnlvls, 1, kjit, &
         temp_atmos_pres_grid, 'scatter', nbp_glo, index_g)

    CALL restput_p (rest_id, 'q_atmos_pres_grid', nbp_glo, jnlvls, 1, kjit, &
         q_atmos_pres_grid, 'scatter', nbp_glo, index_g)

    CALL restput_p (rest_id, 'temp_leaf_pres_grid', nbp_glo, jnlvls, 1, kjit, &
         temp_leaf_pres_grid, 'scatter', nbp_glo, index_g)

    CALL restput_p (rest_id, 'temp_surf_pres', nbp_glo, 1, 1, kjit, &
         temp_surf_pres, 'scatter', nbp_glo, index_g)

    IF (.NOT. hack_can_structure) THEN
       CALL restput_p(rest_id, 'u_speed', nbp_glo, jnlvls, 1, kjit,  u_speed, &
            'scatter',  nbp_glo, index_g)

       CALL restput_p(rest_id, 'z_array_out', nbp_glo, nvm, ncirc, nlevels_tot, kjit,  z_array_out, &
            'scatter',  nbp_glo, index_g)
    END IF

  END SUBROUTINE mleb_finalize

END MODULE mleb
