!================================================================================================================================
!! MODULE   : initsnowpack
!!
!>\BRIEF        Initialize snowpack from parameters calculated with observations (Arthern et al)
!!                
!! DESCRIPTION  : Initialize snowpack from parameters calculated with observations (Arthern et al)
!!
!! \n
!_
!================================================================================================================================



MODULE initsnowpack
    USE ioipsl_para
    USE constantes_soil
    USE constantes
    USE time, ONLY : one_day, dt_sechiba
    USE pft_PARAMETERs 
    USE qsat_moisture 
    USE sechiba_io_p
    USE xios_orchidee
    USE grid
 
    ! TODO : BIND C -> rename routines
    IMPLICIT NONE
    PRIVATE
    PUBLIC initsnowpack_snowtemperature_antarctica, initsnowpack_snowtemperature_greenland, &
         initsnowpack_snowdz_greenland, initsnowpack_logaccu_antarctica, &
         initsnowpack_surfacedensity_antarctica, initsnowpack_snowdensity, &
         initsnowpack_logaccu_greenland, initsnowpack_surfacedensity_greenland

    ! constants for density profiles
    ! ==============================
    ! C0, C1_init : constants for activation energy, 0.07 for rho <= 550 kg m-3, 0.03 for rho >= 550 kg m-3
    REAL(r_std), PARAMETER :: C0 = 0.07
    REAL(r_std), PARAMETER :: C1_init = 0.03
    ! alpha : model-to-observation ratio (MO) snow density correction from Ligtenberg et al 2011
    ! define minimal value for alpha
    REAL(r_std), PARAMETER :: alpha_min = 0.25
    ! alpha0 = beta0 - gamma0 * logaccu : alpha for rho <= 550 kg m-3 (Ligtenberg et al 2011)
    REAL(r_std), PARAMETER :: beta0 = 1.435
    REAL(r_std), PARAMETER :: gamma0 = 0.151
    ! alpha1 = beta1 - gamma1 * logaccu : alpha for rho >= 550 kg m-3 (Ligtenberg et al 2011)
    ! alpha1 = beta1 * beta1_mult - gamma1 * gamma1_mult * logaccu -> adaptation by C. Agosta to fit observed z830
    REAL(r_std), PARAMETER :: beta1 = 2.366
    REAL(r_std), PARAMETER :: gamma1 = 0.293
    REAL(r_std), PARAMETER :: beta1_mult = 0.9
    REAL(r_std), PARAMETER :: gamma1_mult = 0.9
    ! rho_ice_init : ice density [kg m-3]
    REAL(r_std), PARAMETER :: rho_ice_init = 917.
    ! g : gravity [m s-2]
    REAL(r_std), PARAMETER :: gravit = 9.81
    ! E_c : activation energy [J mol-1]
    REAL(r_std), PARAMETER :: E_c = 60000
    ! E_g : activation energy [J mol-1]
    REAL(r_std), PARAMETER :: E_g = 42400
    ! Ra : gas constant [J mol-1 K-1]
    REAL(r_std), PARAMETER :: Ra = 8.3144621
    ! TzeroC2K : Temperature of 0°C in °K
    REAL(r_std), PARAMETER :: TzeroC2K = 273.15
    ! rho550 = 550 kg m-3
    REAL(r_std), PARAMETER :: rho550 = 550.
    ! fill value for data initialization
    REAL(r_std) :: nan_value
    
    TYPE PARAMETERs_snowtemperature
        ! ========================================================
        ! snowtemperature = lat_scale * (lat - lat_ref) - sh_scale * (sh - sh_ref) / 1000. + st_ref
        ! * sh_scale is the dependency of mean annual temperature to the altitude (in °C/km)
        ! * lat_scale is the dependency of mean annual temperature to the latitude (in °C/°N)
        ! * lat_scale and sh_scale are linear between sh_bnd_min->val_scale_min and sh_bnd_max->val_scale_max,
        !   with fixed values bellow sh_bnd_min=val_scale_min and above sh_bnd_max=val_scale_max
        ! ========================================================
        ! st_ref: fixed (known) mean surface temperature (in °C) at the latitude lat_ref and altitude sh_ref
        REAL(r_std) st_ref
        ! lat_ref: latitude of known mean surface temperature (in °N)
        REAL(r_std) lat_ref
        ! sh_ref: altitude of known mean surface temperature (in m asl)
        REAL(r_std) sh_ref
        ! sh_bnd_min: lower bound of altitude for lat_scale and sh_scale linear variations (in m asl)
        REAL(r_std) sh_bnd_min
        ! sh_bnd_max: upper bound of altitude for lat_scale and sh_scale linear variations (in m asl)
        REAL(r_std) sh_bnd_max
        ! lat_scale_min: value of lat_scale at sh_bnd_min (in °C / °N)
        REAL(r_std) lat_scale_min
        ! lat_scale_max: value of lat_scale at sh_bnd_max (in °C / °N)
        REAL(r_std) lat_scale_max
        ! sh_scale_min: value of sh_scale at sh_bnd_min (in °C / km)
        REAL(r_std) sh_scale_min
        ! sh_scale_max: value of sh_scale at sh_bnd_max (in °C / km)
        REAL(r_std) sh_scale_max
    END TYPE PARAMETERs_snowtemperature
    
    TYPE PARAMETERs_logaccu
        ! logaccu_ref: fixed (known) logaccu (in ln(kg m-2 yr-1)) at the mean annual temperature st_ref
        REAL(r_std) logaccu_ref
        ! st_ref: altitude of known mean surface temperature (in m asl)
        REAL(r_std) st_ref
        ! st_scale: logaccu increase with mean annual temperature (in ln(kg m-2 yr-1) / °C)
        REAL(r_std) st_scale
        ! logaccu_min: lower bound for logaccu (in ln(kg m-2 yr-1)) / negative for no lower bound
        REAL(r_std) logaccu_min
        ! logaccu_max: upper bound for logaccu (in ln(kg m-2 yr-1)) / negative for no upper bound
        REAL(r_std) logaccu_max
    END TYPE PARAMETERs_logaccu
    
    TYPE PARAMETERs_surfacedensity
        ! st_ref_min : temperature at which surfacedensity = rhos_min (in °C)
        REAL(r_std) st_ref_min
        ! st_ref_max : temperature at which surfacedensity = rhos_max (in °C)
        REAL(r_std) st_ref_max
        ! rhos_ref_min : surface firn density at temperature st_ref_min (in kg m-3)
        REAL(r_std) rhos_ref_min
        ! rhos_ref_max : surface firn density at temperature st_ref_max (in kg m-3)
        REAL(r_std) rhos_ref_max
        ! rhos_min : lower bound for surfacedensity (in kg m-3) / negative for no lower bound
        REAL(r_std) rhos_min
        ! rhos_max : upper bound for surfacedensity (in kg m-3) / negative for no upper bound
        REAL(r_std) rhos_max
    END TYPE PARAMETERs_surfacedensity

    TYPE PARAMETERs_snowdz
        ! st_ref_min : temperature at which snowdz =snowdz_min (in °C)
        REAL(r_std) st_ref_min
        ! st_ref_max : temperature at which snowdz = snowdz_max (in °C)
        REAL(r_std) st_ref_max
        ! sdz_ref_min : snowdz at temperature st_ref_min (in kg m-3)
        REAL(r_std) sdz_ref_min
        ! sdz_ref_max : snowdz at temperature st_ref_max (in kg m-3)
        REAL(r_std) sdz_ref_max
        ! sdz_min : lower bound for snowdz (in m) 
        REAL(r_std) sdz_min
        ! sdz_max : upper bound for snowdz (in m) 
        REAL(r_std) sdz_max
    END TYPE PARAMETERs_snowdz
     
    
 CONTAINS
    
     SUBROUTINE default_param_snowtemperature_antarctica(param_snowtemperature)
        ! PARAMETERization mean temperature = function(latitude, altitude) for Antarctica
        ! ===============================================================================
        ! From observations and MAR, we find 2 regimes in Antarctica:
        ! * Bellow 500 m a.s.l. : from 1.3 °C/°lat and 6.5 °C/km
        ! * Above 500 m a.s.l. : 0.6 °C/°lat and 9.8 °C/km
        ! param_snowtemperature : TYPE containing default values of param_snowtemperature for Antarctica
        IMPLICIT NONE
        TYPE(PARAMETERs_snowtemperature), INTENT(out) :: param_snowtemperature
        param_snowtemperature%st_ref = 0.
        param_snowtemperature%lat_ref = -57.04
        param_snowtemperature%sh_ref = 0.
        param_snowtemperature%sh_bnd_min = 0.
        param_snowtemperature%sh_bnd_max = 500.
        param_snowtemperature%lat_scale_min = 1.141
        param_snowtemperature%lat_scale_max = 0.734
        param_snowtemperature%sh_scale_min = 6.5
        param_snowtemperature%sh_scale_max = 9.8
    END SUBROUTINE default_param_snowtemperature_antarctica

    SUBROUTINE default_param_snowtemperature_greenland(param_snowtemperature)
        ! PARAMETERization mean temperature = function(latitude, altitude) for Greenland
        ! ===============================================================================
        ! From observations and Hanna et al (2005), we find 2 regimes in Greenland:
        ! * Bellow 1000 m a.s.l. : 5.9 °C/km
        ! * Above 1000 m a.s.l. : 8.2 °C/km
        ! param_snowtemperature : TYPE containing default values of param_snowtemperature for Greenland
        IMPLICIT NONE
        TYPE(PARAMETERs_snowtemperature), INTENT(out) :: param_snowtemperature
        param_snowtemperature%st_ref = 1.108
        param_snowtemperature%lat_ref = 60
        param_snowtemperature%sh_ref = 1000.24
        param_snowtemperature%sh_bnd_min = 0.
        param_snowtemperature%sh_bnd_max = 1000.
        param_snowtemperature%lat_scale_min = -0.658
        param_snowtemperature%lat_scale_max = -0.916
        param_snowtemperature%sh_scale_min = 5.9
        param_snowtemperature%sh_scale_max = 8.2
    END SUBROUTINE default_param_snowtemperature_greenland
    
    SUBROUTINE init_snowtemperature(kjpindex, ice_sheet_mask, topography, lalo, param_snowtemperature, snowtemperature)
        ! ================================================================================
        ! Parametrisation of mean annual temperature as a function of (altitude, latitude)
        ! ================================================================================
        ! -----------------------------------------------------------
        ! Theoritical temperature gradient with altitude and latitude
        ! -----------------------------------------------------------
        ! from Feulner et al., 2013 (DOI: 10.1175/JCLI-D-12-00636.1)
        ! Fig. 3 and 5 : the lapse rate vs. latitude at high latitude is about 0.55 °C °lat-1
        ! with a moist-adiabatic lapse rate of 5 °C km-1 everywhere except for Antarctica,
        ! for Antarctica, a dry-adiabatic lapse rate of 9.8 °C km-1 is assumed.
        ! ------------------------------------------------------------
        ! Parametrized temperature gradient with altitude and latitude
        ! ------------------------------------------------------------
        ! Deduced from observations, MAR outputs and physical considerations
        ! snowtemperature = lat_scale * (lat - lat_ref) - sh_scale * (sh - sh_ref) / 1000. + st_ref
        ! * st_ref is the mean annual surface temperature (in °C) at the latitude lat_ref (in °N) and altitude sh_ref (in m)
        ! * lat_scale is the dependency of mean annual temperature to the latitude (in °C/°N)
        ! * sh_scale is the dependency of mean annual temperature to the altitude (in °C/km)
        ! lat_scale and sh_scale are linear between sh_bnd_min->val_scale_min and sh_bnd_max->val_scale_max,
        ! with fixed values bellow sh_bnd_min=val_scale_min and above sh_bnd_max=val_scale_max
        IMPLICIT NONE
        ! inputs
        ! ======
        INTEGER(i_std), INTENT(in)                           :: kjpindex       !! Domain size
        INTEGER(i_std), DIMENSION (kjpindex), INTENT(in)     :: ice_sheet_mask !! Ice sheet mask 0=no ice, 1=ice
        REAL(r_std), DIMENSION (kjpindex), INTENT(in)        :: topography     !! Topography [m]
        REAL(r_std), DIMENSION(kjpindex,2), INTENT(in)       :: lalo
        TYPE(PARAMETERs_snowtemperature), INTENT(in) :: param_snowtemperature  !! param_snowtemperature: TYPE containing PARAMETERs for snowtemperature
        ! output
        ! ======
        
        REAL(r_std), DIMENSION (kjpindex), INTENT(inout)       :: snowtemperature !! snowtemperature : parametrised mean annual surface temperature (in °C)
        ! local
        ! =====
        INTEGER(i_std)                                       :: ji
        REAL(r_std) lat_scale, sh_scale
        REAL(r_std) st_ref, lat_ref, sh_ref, sh_bnd_min, sh_bnd_max, lat_scale_min, lat_scale_max, sh_scale_min, sh_scale_max
        
        ! load PARAMETERs values
        st_ref = param_snowtemperature%st_ref
        lat_ref = param_snowtemperature%lat_ref
        sh_ref = param_snowtemperature%sh_ref
        sh_bnd_min = param_snowtemperature%sh_bnd_min
        sh_bnd_max = param_snowtemperature%sh_bnd_max
        lat_scale_min = param_snowtemperature%lat_scale_min
        lat_scale_max = param_snowtemperature%lat_scale_max
        sh_scale_min = param_snowtemperature%sh_scale_min
        sh_scale_max = param_snowtemperature%sh_scale_max
        
        lat_scale = 0.
        sh_scale = 0.
        DO ji=1, kjpindex
          IF (ice_sheet_mask(ji) .GT. 0.0) THEN
             ! lat_scale : st dependency to latitude (°C/°N)
             lat_scale = (lat_scale_max - lat_scale_min) / (sh_bnd_max - sh_bnd_min) * &
                   (topography(ji) - sh_bnd_min) + lat_scale_min
             lat_scale = min(lat_scale, lat_scale_min)
             lat_scale = max(lat_scale, lat_scale_max)
             ! sh_scale : st dependency to elevation (°C/km)
             sh_scale = (sh_scale_max - sh_scale_min) / (sh_bnd_max - sh_bnd_min) * &
                   (topography(ji) - sh_bnd_min) + sh_scale_min
             sh_scale = min(sh_scale, sh_scale_max)
             sh_scale = max(sh_scale, sh_scale_min)
             ! snowtemperature : parametrised mean annual surface temperature (in °C)
             snowtemperature(ji) = lat_scale * (lalo(ji,1) - lat_ref) - &
                   sh_scale * (topography(ji) - sh_ref) / 1000. + st_ref
             IF (snowtemperature(ji) .GT. 0) THEN
                snowtemperature(ji) = 0
             END IF

           END IF
        END DO
    
    END SUBROUTINE init_snowtemperature
    
    SUBROUTINE initsnowpack_snowtemperature_antarctica(kjpindex, ice_sheet_mask, topography, lalo, snowtemperature)
        ! ===========================================================
        ! init_snowtemperature with default PARAMETERs for Antarctica
        ! ===========================================================
        IMPLICIT NONE
        ! inputs
        ! ======
        INTEGER(i_std), INTENT(in)                           :: kjpindex       !! Domain size
        INTEGER(i_std), DIMENSION (kjpindex), INTENT(in)     :: ice_sheet_mask !! Ice sheet mask 0=no ice, 1=ice
        REAL(r_std), DIMENSION (kjpindex), INTENT(in)        :: topography     !! Topography [m]
        REAL(r_std), DIMENSION(kjpindex,2), INTENT(in)       :: lalo
        ! output
        ! ======
        REAL(r_std), DIMENSION (kjpindex), INTENT(out)       :: snowtemperature(kjpindex) !! snowtemperature : parametrised mean annual surface temperature (in °C)
        ! local
        ! =====
        ! param_snowtemperature: TYPE containing PARAMETERs for snowtemperature
        TYPE(PARAMETERs_snowtemperature) :: param_snowtemperature
        
        CALL default_param_snowtemperature_antarctica(param_snowtemperature)
        CALL init_snowtemperature(kjpindex, ice_sheet_mask, topography, lalo, param_snowtemperature, snowtemperature)
    END SUBROUTINE initsnowpack_snowtemperature_antarctica
 
    SUBROUTINE initsnowpack_snowtemperature_greenland(kjpindex, ice_sheet_mask, topography, lalo, snowtemperature)
        ! ===========================================================
        ! init_snowtemperature with default PARAMETERs for Greenland
        ! ===========================================================
        IMPLICIT NONE
        ! inputs
        ! ======
        INTEGER(i_std), INTENT(in)                           :: kjpindex       !! Domain size
        INTEGER(i_std), DIMENSION (kjpindex), INTENT(in)     :: ice_sheet_mask !! Ice sheet mask 0=no ice, 1=ice
        REAL(r_std), DIMENSION (kjpindex), INTENT(in)        :: topography     !! Topography [m]
        REAL(r_std), DIMENSION(kjpindex,2), INTENT(in)       :: lalo
        ! output
        ! ======
        ! snowtemperature : parametrised mean annual surface temperature (in °C)
        REAL(r_std), DIMENSION (kjpindex), INTENT(out)       :: snowtemperature !! snowtemperature : parametrised mean annual surface temperature (in °C)
        ! local
        ! =====
        ! param_snowtemperature: TYPE containing PARAMETERs for snowtemperature
        TYPE(PARAMETERs_snowtemperature) :: param_snowtemperature
        
        CALL default_param_snowtemperature_greenland(param_snowtemperature)
        CALL init_snowtemperature(kjpindex, ice_sheet_mask, topography, lalo, param_snowtemperature, snowtemperature)
    END SUBROUTINE initsnowpack_snowtemperature_greenland
    
    
    SUBROUTINE default_param_logaccu_antarctica(param_logaccu)
        ! # PARAMETERization ln(accumulation) = function(temperature) for Antarctica
        ! ==========================================================================
        ! logaccu = st_scale * (st - st_ref) + logaccu_ref
        IMPLICIT NONE
        REAL(r_std), PARAMETER :: accu_min = 10.        ! accu_min : minimal accumulation considered (in kg m-2 yr-1), must be > 0 for computing ln(accumulation)
        ! output
        ! ======
        ! TYPE containing default values of param_logaccu for Antarctica
        TYPE(PARAMETERs_logaccu), INTENT(out) :: param_logaccu
        
        param_logaccu%st_ref = -59.
        param_logaccu%logaccu_ref = 3.434
        param_logaccu%st_scale = 0.056
        param_logaccu%logaccu_min = 3.434 ! negative means no upper bound
        param_logaccu%logaccu_max = -1 ! negative means no upper bound
    END SUBROUTINE default_param_logaccu_antarctica


    SUBROUTINE default_param_logaccu_greenland(param_logaccu)
        ! # PARAMETERization ln(accumulation) = function(temperature) for Greenland
        ! ==========================================================================
        ! logaccu = st_scale * (st - st_ref) + logaccu_ref
        IMPLICIT NONE
        ! accu_min : minimal accumulation considered (in kg m-2 yr-1), must be > 0 for computing ln(accumulation)
        REAL(r_std), PARAMETER :: accu_min = 10.
        ! output
        ! ======
        ! TYPE containing default values of param_logaccu for Antarctica
        TYPE(PARAMETERs_logaccu), INTENT(out) :: param_logaccu
        
        param_logaccu%st_ref = -43
        param_logaccu%logaccu_ref = 3.594
        param_logaccu%st_scale = 0.064
        param_logaccu%logaccu_min = 5.095 ! negative means no upper bound
        param_logaccu%logaccu_max = -1 ! negative means no upper bound
    END SUBROUTINE default_param_logaccu_greenland
 
    SUBROUTINE init_logaccu(kjpindex, ice_sheet_mask, st, param_logaccu, logaccu)
        ! ====================================================================================
        ! # Parametrisation of snow accumulation ln(accumulation) as a function of temperature
        ! ====================================================================================
        ! ---------------------------------------------------------
        ! ## Theoritical ln(accumulation) gradient with temperature
        ! ---------------------------------------------------------
        ! ln(accumulation) is linear with temperature when humidity reach saturation, which is usual in polar regions.
        ! -----------------------------------------------------------
        ! ## PARAMETERized ln(accumulation) gradient with temperature
        ! -----------------------------------------------------------
        ! Deduced from MAR output and physical considerations
        ! logaccu = st_scale * (st - st_ref) + logaccu_ref
        ! * logaccu_ref is the mean annual ln(accumulation) (in ln(kg m-2 yr-1)) at the mean annual temperature st_ref (in °C)
        ! * st_scale is the dependency of ln(accumulation) to mean annual temperature  (in ln(kg m-2 yr-1)/°C)
        ! * logaccu_min and logaccu_max are the minimal and maximal bounds for logaccu
        IMPLICIT NONE
        ! inputs
        ! ======
        INTEGER(i_std), INTENT(in)                           :: kjpindex       !! Domain size
        INTEGER(i_std), DIMENSION (kjpindex), INTENT(in)     :: ice_sheet_mask !! Ice sheet mask 0=no ice, 1=ice
        REAL(r_std), DIMENSION(kjpindex), INTENT(in)         :: st             !! st: mean annual temperature (in °C)
        ! param_logaccu: TYPE containing PARAMETERs for logaccu
        TYPE(PARAMETERs_logaccu), INTENT(in) :: param_logaccu
        ! output
        ! ======
        REAL(r_std), DIMENSION(kjpindex), INTENT(inout)        :: logaccu        !! logaccu: parametrised ln(accumulation) (in ln(kg m-2 yr-1))
        ! local
        ! =====
        INTEGER(i_std)                                       :: ji
        ! logaccu PARAMETERs
        REAL(r_std) logaccu_ref, st_ref, st_scale, logaccu_min, logaccu_max
        
        ! load PARAMETERs
        logaccu_ref = param_logaccu%logaccu_ref
        st_ref = param_logaccu%st_ref
        st_scale = param_logaccu%st_scale
        logaccu_min = param_logaccu%logaccu_min
        logaccu_max = param_logaccu%logaccu_max
        
        ! logaccu : parametrised ln(accumulation)
        DO ji=1, kjpindex
          IF (ice_sheet_mask(ji) .GT. 0.0) THEN
                logaccu(ji) = st_scale * (st(ji) - st_ref) + logaccu_ref
                IF (logaccu_min > 0.) THEN
                   logaccu(ji) = max(logaccu(ji), logaccu_min)
                END IF
                IF (logaccu_max > 0.) THEN
                   logaccu(ji) = min(logaccu(ji), logaccu_max)
                END IF
          END IF
        END DO
    END SUBROUTINE init_logaccu
    
    
    SUBROUTINE initsnowpack_logaccu_antarctica(kjpindex, ice_sheet_mask, st, logaccu)
        ! ===================================================
        ! init_logaccu with default PARAMETERs for Antarctica
        ! ===================================================
        IMPLICIT NONE
        ! inputs
        ! ======
        INTEGER(i_std), INTENT(in)                           :: kjpindex       !! Domain size
        INTEGER(i_std), DIMENSION (kjpindex), INTENT(in)     :: ice_sheet_mask !! Ice sheet mask 0=no ice, 1=ice
        REAL(r_std), DIMENSION(kjpindex), INTENT(in)         :: st             !! st: mean annual temperature (in °C)
        ! output
        ! ======
        REAL(r_std), DIMENSION(kjpindex), INTENT(inout)        :: logaccu        !! logaccu: parametrised ln(accumulation) (in ln(kg m-2 yr-1))
        ! local
        ! =====
        ! param_logaccu: TYPE containing PARAMETERs for logaccu
        TYPE(PARAMETERs_logaccu) :: param_logaccu
        
        CALL default_param_logaccu_antarctica(param_logaccu)
        CALL init_logaccu(kjpindex, ice_sheet_mask, st, param_logaccu, logaccu)
    
    END SUBROUTINE initsnowpack_logaccu_antarctica
 
    SUBROUTINE initsnowpack_logaccu_greenland(kjpindex, ice_sheet_mask, st, logaccu)
        ! ===================================================
        ! init_logaccu with default PARAMETERs for Greenland
        ! ===================================================
        IMPLICIT NONE
        ! inputs
        ! ======
        INTEGER(i_std), INTENT(in)                           :: kjpindex       !! Domain size
        INTEGER(i_std), DIMENSION (kjpindex), INTENT(in)     :: ice_sheet_mask !! Ice sheet mask 0=no ice, 1=ice
        REAL(r_std), DIMENSION(kjpindex), INTENT(in)         :: st             !! st: mean annual temperature (in °C)
        ! output
        ! ======
        REAL(r_std), DIMENSION(kjpindex), INTENT(inout)        :: logaccu        !! logaccu: parametrised ln(accumulation) (in ln(kg m-2 yr-1))
 
        ! local
        ! =====
        ! param_logaccu: TYPE containing PARAMETERs for logaccu
        TYPE(PARAMETERs_logaccu) :: param_logaccu
        
        CALL default_param_logaccu_greenland(param_logaccu)
        CALL init_logaccu(kjpindex, ice_sheet_mask, st, param_logaccu, logaccu)
    
    END SUBROUTINE initsnowpack_logaccu_greenland
    
    SUBROUTINE default_param_surfacedensity_antarctica(param_surfacedensity)
        ! # PARAMETERization of surface firn density = function(temperature) for Antarctica
        ! =================================================================================
        IMPLICIT NONE
        TYPE(PARAMETERs_surfacedensity), INTENT(out) :: param_surfacedensity
        param_surfacedensity%st_ref_min = -46.511
        param_surfacedensity%st_ref_max = -3.863
        param_surfacedensity%rhos_ref_min = 340.614
        param_surfacedensity%rhos_ref_max = 446.89
        param_surfacedensity%rhos_min = 345.38
        param_surfacedensity%rhos_max = 450.0
    END SUBROUTINE default_param_surfacedensity_antarctica 

    SUBROUTINE default_param_surfacedensity_greenland(param_surfacedensity)
        ! # PARAMETERization of surface firn density = function(temperature) for Greenland
        ! =================================================================================
        IMPLICIT NONE
        TYPE(PARAMETERs_surfacedensity), INTENT(out) :: param_surfacedensity
        param_surfacedensity%st_ref_min = -28.98
        param_surfacedensity%st_ref_max = 12.01
        param_surfacedensity%rhos_ref_min = 310.38
        param_surfacedensity%rhos_ref_max = 415.73
        param_surfacedensity%rhos_min = 301.543
        param_surfacedensity%rhos_max = 398.34
    END SUBROUTINE default_param_surfacedensity_greenland
    
    SUBROUTINE init_surfacedensity(kjpindex, ice_sheet_mask, st, param_surfacedensity, surfacedensity)
        ! ======================================================================
        ! # Parametrisation of surface snow density as a function of temperature
        ! ======================================================================
        ! -----------------------------------------------------------------
        ! ## Theoritical surface firn density and gradient with temperature
        ! -----------------------------------------------------------------
        ! Surface firn density is defined as the mean first meter density of firn.
        ! At Dome C, mean annual temperature is -45 °C and surface firn density is 320. kg m-3
        ! Near Dumont d'Urville, mean annual temperature is -15 °C and surface firn density is 450. kg m-3
        ! ------------------------------------
        ! ## Parametrized surface firn density
        ! ------------------------------------
        ! surfacedensity = (rhos_ref_max - rhos_ref_min) / (st_ref_max - st_ref_min)
        ! * surfacedensity is linear with temperature, between (st_ref_min, rhos_ref_min) and (st_ref_max, rhos_ref_max)
        ! * surfacedensity minimum value is rhos_min and surfacedensity maximum value is rhos_max
        IMPLICIT NONE
 
        ! inputs
        ! ======
        INTEGER(i_std), INTENT(in)                           :: kjpindex       !! Domain size
        INTEGER(i_std), DIMENSION (kjpindex), INTENT(in)     :: ice_sheet_mask !! Ice sheet mask 0=no ice, 1=ice
        REAL(r_std), DIMENSION(kjpindex), INTENT(in)         :: st             !! st: mean annual temperature (in °C)
        ! param_surfacedensity : TYPE containing PARAMETERs for surface density
        TYPE(PARAMETERs_surfacedensity), INTENT(in) :: param_surfacedensity
        ! output
        ! ======
        REAL(r_std), DIMENSION(kjpindex), INTENT(inout)         :: surfacedensity!! surfacedensity : parametrised surface firn density (in kg m-3)
        ! local
        ! =====
        INTEGER(i_std)                                       :: ji
        ! surfacedensity PARAMETERs
        REAL(r_std) st_ref_min, st_ref_max, rhos_ref_min, rhos_ref_max, rhos_min, rhos_max
        
        ! load PARAMETERs values
        st_ref_min = param_surfacedensity%st_ref_min
        st_ref_max = param_surfacedensity%st_ref_max
        rhos_ref_min = param_surfacedensity%rhos_ref_min
        rhos_ref_max = param_surfacedensity%rhos_ref_max
        rhos_min = param_surfacedensity%rhos_min
        rhos_max = param_surfacedensity%rhos_max
        
        DO ji=1, kjpindex
          IF (ice_sheet_mask(ji) .GT. 0.0) THEN
             ! surfacedensity : parametrised surface firn density
             surfacedensity(ji) = (rhos_ref_max - rhos_ref_min) / (st_ref_max - st_ref_min) * &
                   (st(ji) - st_ref_min) + rhos_ref_min
             IF (rhos_min > 0.) THEN
                surfacedensity(ji) = max(surfacedensity(ji), rhos_min)
             END IF
             IF (rhos_max > 0.) THEN
                surfacedensity(ji) = min(surfacedensity(ji), rhos_max)
             END IF
          END IF
        END DO
    END SUBROUTINE init_surfacedensity
    
    
    SUBROUTINE initsnowpack_surfacedensity_antarctica(kjpindex, ice_sheet_mask, st, surfacedensity)
        ! ==========================================================
        ! init_surfacedensity with default PARAMETERs for Antarctica
        ! ==========================================================
        IMPLICIT NONE
 
        INTEGER(i_std), INTENT(in)                           :: kjpindex       !! Domain size
        INTEGER(i_std), DIMENSION (kjpindex), INTENT(in)     :: ice_sheet_mask !! Ice sheet mask 0=no ice, 1=ice
        REAL(r_std), DIMENSION(kjpindex), INTENT(in)         :: st             !! st: mean annual temperature (in °C)
        ! output
        ! ======
        REAL(r_std), DIMENSION(kjpindex), INTENT(inout)         :: surfacedensity!! surfacedensity : parametrised surface firn density (in kg m-3)
        ! local
        ! =====
        ! param_surfacedensity : TYPE containing PARAMETERs for surface density
        TYPE(PARAMETERs_surfacedensity) :: param_surfacedensity
        
        CALL default_param_surfacedensity_antarctica(param_surfacedensity)
        CALL init_surfacedensity(kjpindex, ice_sheet_mask, st, param_surfacedensity, surfacedensity)
    
    END SUBROUTINE initsnowpack_surfacedensity_antarctica
 
    SUBROUTINE initsnowpack_surfacedensity_greenland(kjpindex, ice_sheet_mask, st, surfacedensity)
        ! ==========================================================
        ! init_surfacedensity with default PARAMETERs for Greenland
        ! ==========================================================
        IMPLICIT NONE
        INTEGER(i_std), INTENT(in)                           :: kjpindex       !! Domain size
        INTEGER(i_std), DIMENSION (kjpindex), INTENT(in)     :: ice_sheet_mask !! Ice sheet mask 0=no ice, 1=ice
        REAL(r_std), DIMENSION(kjpindex), INTENT(in)         :: st             !! st: mean annual temperature (in °C)
        ! output
        ! ======
        REAL(r_std), DIMENSION(kjpindex), INTENT(inout)         :: surfacedensity!! surfacedensity : parametrised surface firn density (in kg m-3)
        ! local
        ! =====
        ! param_surfacedensity : TYPE containing PARAMETERs for surface density
        TYPE(PARAMETERs_surfacedensity) :: param_surfacedensity
        
        CALL default_param_surfacedensity_greenland(param_surfacedensity)
        CALL init_surfacedensity(kjpindex, ice_sheet_mask, st, param_surfacedensity, surfacedensity)
    
    END SUBROUTINE initsnowpack_surfacedensity_greenland

        SUBROUTINE default_param_snowdz_greenland(param_snowdz)
        ! # PARAMETERization of surface firn density = function(temperature) for Greenland
        ! =================================================================================
        IMPLICIT NONE
        TYPE(PARAMETERs_snowdz), INTENT(out) :: param_snowdz
        param_snowdz%st_ref_min = -14.7
        param_snowdz%st_ref_max = -8.7
        param_snowdz%sdz_ref_min = 0.5
        param_snowdz%sdz_ref_max = 10.
        param_snowdz%sdz_min = 0.5
        param_snowdz%sdz_max = 10.
    END SUBROUTINE default_param_snowdz_greenland

    SUBROUTINE init_snowdz(kjpindex, ice_sheet_mask, st, param_snowdz, snowdz)
        ! ===================================================================================
        ! # Parametrisation of snowdz as a function of mean surface temperature
        ! ===================================================================================
        ! ------------------------------------
        ! ## Parametrized snowpack thickness
        ! ------------------------------------
        ! snowdz = sdz_ref_max - (sdz_ref_max - sdz_ref_min) / (st_ref_max - st_ref_min) * (st - st_ref_min)
        IMPLICIT NONE

        ! inputs
        ! ======
        INTEGER(i_std), INTENT(in)                           :: kjpindex       !! Domain size
        INTEGER(i_std), DIMENSION (kjpindex), INTENT(in)     :: ice_sheet_mask !! Ice sheet mask 0=no ice, 1=ice
        REAL(r_std), DIMENSION(kjpindex), INTENT(in)         :: st             !! surface temperature
        ! param_surfacedensity : TYPE containing PARAMETERs for surface density
        TYPE(PARAMETERs_snowdz), INTENT(in) :: param_snowdz
        ! output
        ! ======
        REAL(r_std), DIMENSION(kjpindex), INTENT(inout)         :: snowdz      !! snowdz : parametrised snow thickness (in m)
        ! local
        ! =====
        INTEGER(i_std)                                       :: ji
        ! snowdz PARAMETERs
        REAL(r_std) st_ref_min, st_ref_max, sdz_ref_min, sdz_ref_max, sdz_min, sdz_max

        ! load PARAMETERs values
        st_ref_min = param_snowdz%st_ref_min
        st_ref_max = param_snowdz%st_ref_max
        sdz_ref_min = param_snowdz%sdz_ref_min
        sdz_ref_max = param_snowdz%sdz_ref_max
        sdz_min = param_snowdz%sdz_min
        sdz_max = param_snowdz%sdz_max

        DO ji=1, kjpindex
          IF (ice_sheet_mask(ji) .GT. 0.0) THEN
             ! sdz : parametrised snow dz
             snowdz(ji) = sdz_ref_max - (sdz_ref_max - sdz_ref_min) / (st_ref_max - st_ref_min) * (st(ji) - st_ref_min)
             IF (sdz_min > 0.) THEN
                snowdz(ji) = max(snowdz(ji), sdz_min)
             END IF
             IF (sdz_max > 0.) THEN
                snowdz(ji) = min(snowdz(ji), sdz_max)
             END IF
          END IF
        END DO

    END SUBROUTINE init_snowdz

    SUBROUTINE initsnowpack_snowdz_greenland(kjpindex, ice_sheet_mask, st, snowdz)
        ! ==========================================================
        ! init_surfacedensity with default PARAMETERs for Greenland
        ! ==========================================================
        IMPLICIT NONE
        INTEGER(i_std), INTENT(in)                           :: kjpindex       !! Domain size
        INTEGER(i_std), DIMENSION (kjpindex), INTENT(in)     :: ice_sheet_mask !! Ice sheet mask 0=no ice, 1=ice
        REAL(r_std), DIMENSION(kjpindex), INTENT(in)         :: st             !! st: mean annual temperature (in °C)
        ! output
        ! ======
        REAL(r_std), DIMENSION(kjpindex), INTENT(inout)         :: snowdz !! parametrised snowpack thickness (in kg m-3)
        ! local
        ! =====
        ! param_surfacedensity : TYPE containing PARAMETERs for surface density
        TYPE(PARAMETERs_snowdz) :: param_snowdz

        CALL default_param_snowdz_greenland(param_snowdz)
        CALL init_snowdz(kjpindex, ice_sheet_mask, st, param_snowdz, snowdz)

    END SUBROUTINE initsnowpack_snowdz_greenland
    
    SUBROUTINE E_formula(kjpindex, ice_sheet_mask, st, Const, alpha, eformula)
        ! E_formula used for integrated snow density
        IMPLICIT NONE
        ! input
        ! =====
        INTEGER(i_std), INTENT(in)                           :: kjpindex       !! Domain size
        INTEGER(i_std), DIMENSION (kjpindex), INTENT(in)     :: ice_sheet_mask !! Ice sheet mask 0=no ice, 1=ice
        REAL(r_std), DIMENSION(kjpindex), INTENT(in)         :: st             !! st: mean annual temperature (in °C)
        
        REAL(r_std), INTENT(in)                              :: Const          !! Const : Constant for activation energy
        REAL(r_std), DIMENSION (kjpindex), INTENT(in)        :: alpha  !! alpha : model-to-observation ratio (MO) correction from Ligtenberg et al 2011, depend of logaccu
 
        ! output
        ! ======
        REAL(r_std), DIMENSION (kjpindex), INTENT(inout)     :: eformula
        ! local
        ! =====
        INTEGER(i_std)                                       :: ji
        
        DO ji=1, kjpindex
          IF (ice_sheet_mask(ji) .GT. 0.0) THEN
             eformula(ji) = Const * gravit * exp((E_g - E_c) / (Ra * (st(ji) + TzeroC2K))) * &
                    rho_ice_init * alpha(ji)
          END IF  
        END DO
    END SUBROUTINE E_formula
    
    
    SUBROUTINE rho_formula(kjpindex, ice_sheet_mask, z_profile, st, surfacedensity, Const, alpha, rho)
        ! E_formula used for integrated snow density
        IMPLICIT NONE
        ! input
        ! =====
        INTEGER(i_std), INTENT(in)                           :: kjpindex       !! Domain size
        INTEGER(i_std), DIMENSION (kjpindex), INTENT(in)     :: ice_sheet_mask !! Ice sheet mask 0=no ice, 1=ice
        REAL(r_std), DIMENSION (kjpindex,nsnow)              :: z_profile      !! z_profile : depth profile on which to compute the density (in m)
        REAL(r_std), DIMENSION(kjpindex), INTENT(in)         :: st             !! st: mean annual temperature (in °C)
        REAL(r_std), DIMENSION(kjpindex), INTENT(in)         :: surfacedensity !! surfacedensity : parametrised surface firn density (in kg m-3)
        REAL(r_std), INTENT(in)                              :: Const          !! Const : Constant for activation energy
        REAL(r_std), DIMENSION(kjpindex), INTENT(in)         :: alpha          !! alpha : model-to-observation ratio (MO) correction from Ligtenberg et al 2011, depend of logaccu
        ! output
        ! ======
        REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(out) :: rho        !! Snow density (Kg/m^3)
        ! local
        ! =====
        REAL(r_std), DIMENSION (kjpindex)                    :: eformula
        INTEGER(i_std)                                       :: ji,jj
        
        CALL E_formula(kjpindex, ice_sheet_mask, st, Const, alpha, eformula)
        
        DO ji=1,kjpindex
          DO jj=1,nsnow
             IF ((ice_sheet_mask(ji) .GT. 0.0) .and. (z_profile(ji, jj) .GE. 0.)) THEN
                rho(ji, jj) = exp(eformula(ji) * z_profile(ji, jj)) / &
                         (rho_ice_init / surfacedensity(ji) - 1 + exp(eformula(ji) * z_profile(ji, jj))) * rho_ice_init
             END IF
          END DO
        END DO
 
    END SUBROUTINE rho_formula
    
    SUBROUTINE z550(kjpindex, ice_sheet_mask, st, surfacedensity, alpha0, z550_loc)
        ! depth of density 550 km m-3
        IMPLICIT NONE
        INTEGER(i_std), INTENT(in)                           :: kjpindex       !! Domain size
        INTEGER(i_std), DIMENSION (kjpindex), INTENT(in)     :: ice_sheet_mask !! Ice sheet mask 0=no ice, 1=ice
        REAL(r_std), DIMENSION(kjpindex), INTENT(in)         :: st             !! st: mean annual temperature (in °C)
        REAL(r_std), DIMENSION(kjpindex), INTENT(in)         :: surfacedensity !! surfacedensity : parametrised surface firn density (in kg m-3)
        ! alpha0 : MO correction from Ligtenberg et al 2011 for rho <= 550 kg m-3
        REAL(r_std), DIMENSION(kjpindex), INTENT(in)         :: alpha0         !! alpha0 : MO correction from Ligtenberg et al 2011 for rho <= 550 kg m-3
 
        ! surfacedensity
        ! output
        ! ======
        REAL(r_std), DIMENSION (kjpindex), INTENT(inout)     :: z550_loc       !! z550_loc : depth to which the density of 550 kg m-3 is reached (in m)
        ! local
        ! =====
        REAL(r_std), DIMENSION (kjpindex)                    :: eformula
        INTEGER(i_std)                                       :: ji
        
        ! alpha0 : MO correction from Ligtenberg et al 2011 for rho <= 550 kg m-3
        CALL E_formula(kjpindex, ice_sheet_mask, st, C0, alpha0, eformula)
        
        DO ji=1, kjpindex
          IF (ice_sheet_mask(ji) .GT. 0.0) THEN
                    z550_loc(ji) = log((rho_ice_init / surfacedensity(ji) - 1.) / (rho_ice_init / rho550 - 1.)) / eformula(ji)
          END IF
        END DO
    END SUBROUTINE z550
    
    SUBROUTINE alphaMO(kjpindex, ice_sheet_mask, logaccu, beta, gamma, alpha)
        ! alpha(logaccu) model-to-obs ratio from Ligtenberg11
        IMPLICIT NONE
        ! input
        ! =====
        INTEGER(i_std), INTENT(in)                           :: kjpindex       !! Domain size
        INTEGER(i_std), DIMENSION (kjpindex), INTENT(in)     :: ice_sheet_mask !! Ice sheet mask 0=no ice, 1=ice
        ! logaccu: ln(accumulation) (in ln(kg m-2 yr-1))
        REAL(r_std), DIMENSION(kjpindex), INTENT(in)         :: logaccu        !! logaccu: ln(accumulation) (in ln(kg m-2 yr-1))
        ! beta, gamma : factors for logaccu dependency
        REAL(r_std), INTENT(in) :: beta
        REAL(r_std), INTENT(in) :: gamma
        ! output
        ! ======
        ! alpha : alpha(logaccu) model-to-obs ratio
        REAL(r_std), DIMENSION(kjpindex), INTENT(inout)         :: alpha          !! alpha : model-to-observation ratio (MO) correction from Ligtenberg et al 2011, depend of logaccu
        INTEGER(i_std)                                       :: ji
        
        DO ji=1, kjpindex
          IF (ice_sheet_mask(ji) .GT. 0.0) THEN
                ! corrective factor for density
                alpha(ji) = max(beta - gamma * logaccu(ji), alpha_min)
          END IF
       END DO
    END SUBROUTINE alphaMO
    
    SUBROUTINE initsnowpack_snowdensity(kjpindex, ice_sheet_mask, z_profile, st, logaccu, surfacedensity, rho)
        IMPLICIT NONE
        ! input
        ! =====
        INTEGER(i_std), INTENT(in)                           :: kjpindex       !! Domain size
        INTEGER(i_std), DIMENSION (kjpindex), INTENT(in)     :: ice_sheet_mask !! Ice sheet mask 0=no ice, 1=ice
        REAL(r_std), DIMENSION (kjpindex,nsnow)              :: z_profile      !! z_profile : depth profile on which to compute the density (in m)
        REAL(r_std), DIMENSION(kjpindex), INTENT(in)         :: st             !! st: mean annual temperature (in °C)
        REAL(r_std), DIMENSION(kjpindex), INTENT(in)         :: surfacedensity !! surfacedensity : parametrised surface firn density (in kg m-3)
        REAL(r_std), DIMENSION(kjpindex), INTENT(in)         :: logaccu        !! logaccu: ln(accumulation) (in ln(kg m-2 yr-1))
        ! output
        ! ======
        REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(out) :: rho        !! Snow density (Kg/m^3)
        ! local
        ! =====
        REAL(r_std), DIMENSION (kjpindex)                    :: z550_loc       !! z550_loc : depth to which the density of 550 kg m-3 is reached (in m)
        REAL(r_std), DIMENSION(kjpindex)                     :: alpha0         !! alpha0 : MO correction from Ligtenberg et al 2011 for rho <= 550 kg m-3
        REAL(r_std), DIMENSION (kjpindex,nsnow)              :: rho0        !! rho0 : snow density for density <= 550 kg m-3
        REAL(r_std), DIMENSION(kjpindex)                     :: alpha1         !! alpha1 : alpha(logaccu) model-to-obs ratio for density >= 550 kg m-3
        REAL(r_std), DIMENSION (kjpindex,nsnow)              :: dz550_profile      !! dz_profile : z bellow z550 = z - z550
        REAL(r_std), DIMENSION (kjpindex,nsnow)              :: rho1        !! rho1 : snow density for density >= 550 kg m-3
        REAL(r_std), DIMENSION (kjpindex)                    :: rho550_2d        !! rho550_2d = 550 kg m-3
 
        ! loops
        INTEGER(i_std)                                       :: ji,jj
        
        nan_value = -9999.
        
        ! corrective factor Ligtenberg11
        ! ------------------------------
        ! corrective factor for density <= 550 kg m-3
        CALL alphaMO(kjpindex, ice_sheet_mask, logaccu, beta0, gamma0, alpha0)
        ! snow density for z <= 550 kg m-3
        rho0 = nan_value
        CALL rho_formula(kjpindex, ice_sheet_mask, z_profile, st, surfacedensity, C0, alpha0, rho0)
        
        ! define z550_loc to compute the density bellow this depth
        CALL z550(kjpindex, ice_sheet_mask, st, surfacedensity, alpha0, z550_loc)
        
        ! initialize dz550_profile with negative value
        dz550_profile = -1.
        DO ji=1,kjpindex
          DO jj=1,nsnow
             IF (ice_sheet_mask(ji) .GT. 0.0) THEN
                dz550_profile(ji, jj) = z_profile(ji, jj) - z550_loc(ji)
             END IF
          END DO
        END DO
 
        
        ! corrective factor for density >= 550 kg m-3
        CALL alphaMO(kjpindex, ice_sheet_mask, logaccu, beta1 * beta1_mult, gamma1 * gamma1_mult, alpha1)
        ! snow density for z >= 550 kg m-3
        rho550_2d = rho550
        rho1 = nan_value
        CALL rho_formula(kjpindex, ice_sheet_mask, dz550_profile, st, rho550_2d, C1_init, alpha1, rho1)
        
        ! compute rho
        DO ji=1,kjpindex
          DO jj=1,nsnow
             IF (ice_sheet_mask(ji) .GT. 0.0) THEN
                IF (z_profile(ji, jj) <= z550_loc(ji)) THEN
                   rho(ji, jj) = rho0(ji, jj)
                else
                   rho(ji, jj) = rho1(ji, jj)
                END IF
             END IF
          END DO
        END DO
 
    END SUBROUTINE initsnowpack_snowdensity
 
 
 END module initsnowpack
