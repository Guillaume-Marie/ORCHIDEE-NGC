! ================================================================================================================================
!  MODULE       : explicitsnow
!
!  CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
!  LICENCE      : IPSL (2006)
!  This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF  Computes hydrologic snow processes on continental points.
!!
!!\n DESCRIPTION: This module computes hydrologic snow processes on continental points.
!!
!! RECENT CHANGES: 
!!
!! REFERENCES	: None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_sechiba/explicitsnow.f90 $
!! $Date: 2026-01-14 14:32:41 +0100 (mer. 14 janv. 2026) $
!! $Revision: 9306 $
!! \n
!_ ================================================================================================================================

MODULE explicitsnow
  USE ioipsl_para
  USE constantes_soil
  USE constantes
  USE time, ONLY : one_day, dt_sechiba
  USE pft_parameters 
  USE qsat_moisture 
  USE sechiba_io_p
  USE xios_orchidee
  USE grid
  USE initsnowpack, ONLY : initsnowpack_snowtemperature_antarctica, initsnowpack_snowtemperature_greenland, &
       initsnowpack_snowdz_greenland, initsnowpack_logaccu_antarctica, &
       initsnowpack_surfacedensity_antarctica, initsnowpack_snowdensity, &
       initsnowpack_logaccu_greenland, initsnowpack_surfacedensity_greenland

  IMPLICIT NONE

  ! Public routines :
  PRIVATE
  PUBLIC :: explicitsnow_main, explicitsnow_initialize, explicitsnow_xios_initialize, explicitsnow_finalize

  !! Possible values for COMPACT_SNOW_METHOD are COMPACT_SNOW_METHOD, _DEFAULTCOMPACT_SNOW_DECHARME16 and COMPACT_SNOW_KAMPENHOUT
  INTEGER(i_std), PARAMETER                    :: COMPACT_SNOW_METHOD_DEFAULT = 1 !! Index for original method
  INTEGER(i_std), PARAMETER                    :: COMPACT_SNOW_DECHARME16 = 2     !! Index to follow the method by Decharme et al 2016
  INTEGER(i_std), PARAMETER                    :: COMPACT_SNOW_KAMPENHOUT = 3     !! Index to follow the method by Kampenhout 2017
  !! Possible values for VISCOSITY_METHOD are VISCOSITY_METHOD_DEFAULT and VISCOSITY_METHOD_VIONNET
  INTEGER(i_std), PARAMETER                    :: VISCOSITY_METHOD_DEFAULT = 1    !! Index for original method Anderson
  INTEGER(i_std), PARAMETER                    :: VISCOSITY_METHOD_VIONNET = 2    !! Index to follow the method by Vionnet et al 2012 Crocus

  
CONTAINS

!================================================================================================================================
!! SUBROUTINE   : explicitsnow_initialize
!!
!>\BRIEF        Read variables for explictsnow module from restart file
!!                
!! DESCRIPTION  : Read variables for explictsnow module from restart file
!!                Initialisation of the snowpack following Arthern et al, 2010 and Ligtenberg et al, 2011 (P. Conesa 2024)
!!
!! \n
!_
!================================================================================================================================
  SUBROUTINE explicitsnow_initialize( kjit,     kjpindex, rest_id,  frac_nobio, &
                                      snowrho,  snowtemp, snowdz,   snowheat, snowgrain, &
                                      icetemp,  icedz,    ice_sheet_mask)

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                           :: kjit           !! Time step number 
    INTEGER(i_std), INTENT(in)                           :: kjpindex       !! Domain size
    INTEGER(i_std),INTENT (in)                           :: rest_id        !! Restart file identifier
    REAL(r_std), DIMENSION (kjpindex,nnobio), INTENT(in) :: frac_nobio     !! Fraction of continental ice, lakes, ...

    !! 0.2 Output variables
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(out)  :: snowrho        !! Snow density (Kg/m^3)
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(out)  :: snowtemp       !! Snow temperature
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(out)  :: snowdz         !! Snow layer thickness [m]
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(out)  :: snowheat       !! Snow heat content (J/m^2)
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(out)  :: snowgrain      !! Snow grainsize (m)
    REAL(r_std), DIMENSION (kjpindex,nice), INTENT(out)   :: icetemp        !! Ice temperature
    REAL(r_std), DIMENSION (kjpindex,nice), INTENT(out)   :: icedz          !! Ice layer thickness [m]
    INTEGER(i_std), DIMENSION (kjpindex), INTENT(out)     :: ice_sheet_mask !! Ice sheet mask 0=no ice, 1=ice
    
    !! local variables
    INTEGER(i_std)                                        :: ji,jj
    REAL(r_std), DIMENSION (kjpindex,nsnow)               :: snowliq        !! Snow liquid water content [m]
    REAL(r_std), DIMENSION (kjpindex)                     :: snow_depth_tmp !! Snow depth [m]
    REAL(r_std), DIMENSION (kjpindex)                     :: topography     !! Topography [m]
    LOGICAL                                               :: read_topography!! Read file topography
    REAL(r_std), DIMENSION (kjpindex)                     :: logaccu        !! Logaccu 
    REAL(r_std), DIMENSION (kjpindex)                     :: surfacetemp    !! Surface Temp [°C]
    REAL(r_std), DIMENSION (kjpindex,nsnow)               :: snowdz_cumsum  !! Snowpack thickness [m]
    REAL(r_std), DIMENSION (kjpindex)                     :: ice_sheet_mask_diag  !! ice_sheet_mask output
    REAL(r_std), DIMENSION (kjpindex)                     :: surfacedensity !! Surface Density [kg.m-3]
    

    !! 1. Read from restart file
    CALL restget_p (rest_id, 'snowrho', nbp_glo, nsnow, 1, kjit,.TRUE.,snowrho, "gather", nbp_glo, index_g)
    CALL setvar_p (snowrho, val_exp, 'Snow Density profile', xrhosmin)
    
    CALL restget_p (rest_id, 'snowtemp', nbp_glo,  nsnow, 1, kjit,.TRUE.,snowtemp, "gather", nbp_glo,  index_g)
    CALL setvar_p (snowtemp, val_exp, 'Snow Temperature profile', tp_00)
    
    CALL restget_p (rest_id, 'snowdz', nbp_glo,  nsnow, 1, kjit,.TRUE.,snowdz, "gather", nbp_glo,  index_g)
    CALL setvar_p (snowdz, val_exp, 'Snow depth profile', 0.0)
    
    CALL restget_p (rest_id, 'snowheat', nbp_glo,  nsnow, 1, kjit,.TRUE.,snowheat, "gather", nbp_glo,  index_g)
    CALL setvar_p (snowheat, val_exp, 'Snow Heat profile', 0.0)
    
    CALL restget_p (rest_id, 'snowgrain', nbp_glo,  nsnow, 1, kjit,.TRUE.,snowgrain, "gather", nbp_glo,  index_g)
    CALL setvar_p (snowgrain, val_exp, 'Snow grain profile', 0.0)
    
    !initialize
    surfacedensity(:) = 0.
    ice_sheet_mask(:) = 0
    
    IF (ok_ice_sheet) THEN
       CALL restget_p (rest_id, 'icetemp', nbp_glo,  nice, 1, kjit,.TRUE.,icetemp, "gather", nbp_glo,  index_g)
    
       IF (ALL(icetemp(:,:)==val_exp)) THEN
          ! icetemp was not found in restart file
          IF (read_reftempice) THEN
             ! Read variable ptn from file
             CALL explicitsnow_read_reftempicefile(kjpindex,lalo,icetemp)
          ELSE
             ! Initialize icetemp with a constant value which can be set in run.def

             !Config Key   = EXPLICITSNOW_TICEPRO
             !Config Desc  = Initial ice temperature profile if not found in restart
             !Config Def   = 273.
             !Config If    = OK_SECHIBA
             !Config Help  = The initial value of the temperature profile in the soil if 
             !Config         its value is not found in the restart file. Here
             !Config         we only require one value as we will assume a constant 
             !Config         throughout the column.
             !Config Units = Kelvin [K]
             CALL setvar_p (icetemp, val_exp,'EXPLICITSNOW_TICEPRO',273._r_std)
          END IF
       END IF
       DO ji=1,kjpindex
          icedz(ji,:) = ZSICOEF1(:)
       ENDDO
      
  
       ! Read ice sheet mask from file or determine it using frac_nobio
       IF (read_ice_sheet_mask) THEN
          ! Read ice sheet mask from file
          IF (printlev>=2) WRITE(numout,*) 'Read ice_sheet_mask from file'
          CALL explicitsnow_read_ice_sheet_mask(kjpindex,lalo,ice_sheet_mask)
       ELSE
          ! Determine the ice sheet mask using the fraction of nobio (ice_sheet_mask will be only 0 or 1)
          IF (printlev>=2) WRITE(numout,*) 'Determine ice_sheet_mask from frac_nobio using frac_min_ice=',frac_min_ice
          DO ji=1, kjpindex
             IF (frac_nobio(ji,iice) >= frac_min_ice) THEN
                ice_sheet_mask(ji) = 1
             ELSE
                ice_sheet_mask(ji) = 0
             END IF
          END DO
       END IF
      
       IF (ALL(snowtemp(:,:)==273.15)) THEN
          IF (INIT_SNOWPACK) THEN
             ! Read relief from file in order to initialize snowpack density
             IF (printlev>=2) WRITE(numout,*) 'Read topography from file'
             CALL explicitsnow_read_relief(kjpindex, topography)

             ! Definition of snow depth in ice sheet points
             DO ji=1, kjpindex
                IF (ice_sheet_mask(ji) .EQ. 2) THEN ! For Antarctica : same snowpack thickness everywhere
                   snow_depth_tmp(ji) = 10.0
                ELSE IF (ice_sheet_mask(ji) .EQ. 1) THEN  !For Greenland : snowpack thickness parameterization 
                   !"Compute snow temperature"
                   call initsnowpack_snowtemperature_greenland(kjpindex, ice_sheet_mask, topography, lalo, surfacetemp) 
                   !"Compute snow depth tmp"
                   call initsnowpack_snowdz_greenland(kjpindex, ice_sheet_mask, surfacetemp, snow_depth_tmp)
                ELSE
                   snow_depth_tmp(ji) = 0.0 
                ENDIF
             END DO

             ! Discretization of snowpack
             CALL explicitsnow_levels(kjpindex,snow_depth_tmp, snowdz)
             IF (printlev>=2) WRITE(numout,*) 'snowpack discretization'

             !Compute snowpack thickness at each layer
             DO ji=1, kjpindex
                IF (ice_sheet_mask(ji) .GT. 0) THEN
                   snowdz_cumsum(ji,1) = snowdz(ji,1)
                   DO jj=2, nsnow
                      snowdz_cumsum(ji,jj) = snowdz_cumsum(ji,jj-1) + snowdz(ji,jj)
                   END DO
                END IF
             END DO
          
             IF (ANY(ice_sheet_mask .EQ. 2)) THEN
                !"Compute snow temperature"
                call initsnowpack_snowtemperature_antarctica(kjpindex, ice_sheet_mask, topography, lalo, surfacetemp) ! add snowtemp
                !"Compute logaccu"
                call initsnowpack_logaccu_antarctica(kjpindex, ice_sheet_mask, surfacetemp, logaccu)
                !"Compute surface density"
                call initsnowpack_surfacedensity_antarctica(kjpindex, ice_sheet_mask, surfacetemp, surfacedensity)
                ! "Compute snow density"
                call initsnowpack_snowdensity(kjpindex, ice_sheet_mask, snowdz_cumsum, surfacetemp, logaccu, surfacedensity, snowrho) ! add snowrho
                IF (printlev>=2) WRITE(numout,*) 'init snowrho'
                ! snowtemperature in the subroutine is initialized in [°C], we need to convert it in K
                IF (printlev>=2) WRITE(numout,*) 'Valeur maximale de snowrho : ',MAXVAL(snowrho)
                IF (printlev>=2) WRITE(numout,*) 'Valeur minimale de snowrho : ',MINVAL(snowrho)
             END IF
             IF (ANY(ice_sheet_mask .EQ. 1)) THEN
                !"Compute snow temperature"
                call initsnowpack_snowtemperature_greenland(kjpindex, ice_sheet_mask, topography, lalo, surfacetemp) ! add snowtemp
                !"Compute logaccu"
                call initsnowpack_logaccu_greenland(kjpindex, ice_sheet_mask, surfacetemp, logaccu)
                !"Compute surface density"
                call initsnowpack_surfacedensity_greenland(kjpindex, ice_sheet_mask, surfacetemp, surfacedensity)
                ! "Compute snow density"
                call initsnowpack_snowdensity(kjpindex, ice_sheet_mask, snowdz_cumsum, surfacetemp, logaccu, surfacedensity, snowrho) ! add snowrho
                IF (printlev>=2) WRITE(numout,*) 'init snowrho gr'
                ! snowtemperature in the subroutine is initialized in [°C], we need to convert it in K
                IF (printlev>=2) WRITE(numout,*) 'Valeur maximale de surfacedensity : ',MAXVAL(surfacedensity)
                IF (printlev>=2) WRITE(numout,*) 'Valeur minimale de surfacedensity : ',MINVAL(surfacedensity)
             END IF
             ! Update snow heat and snowliq
             DO ji = 1, kjpindex
                IF (ice_sheet_mask(ji) .GT. 0) THEN 
                   snowtemp(ji,:) = surfacetemp(ji) + tp_00
                   icetemp(ji,:) = surfacetemp(ji) + tp_00
                   snowliq(ji,:) = 0.0 ! no liquid water
                   snowheat(ji,:) = snow3lheat_1d(snowliq(ji,:),snowrho(ji,:),snowdz(ji,:),snowtemp(ji,:))
                END IF
             END DO
          END IF
       END IF
       ice_sheet_mask_diag(:) = ice_sheet_mask(:) ! to have a real for output
       CALL xios_orchidee_send_field("ice_sheet_mask",ice_sheet_mask_diag)
       CALL xios_orchidee_send_field("surfacedensity",surfacedensity)
    END IF
    

  END SUBROUTINE explicitsnow_initialize

!!  =============================================================================================================================
!! SUBROUTINE:    explicitsnow_xios_initialize
!!
!>\BRIEF	  Initialize xios dependant definition before closing context defintion
!!
!! DESCRIPTION:	  Initialize xios dependant defintion before closing context defintion
!!
!! RECENT CHANGE(S): Initialization of XIOS to read ice_sheet_mask and
!!                   relief init_snowpack maps
!! \n
!_ ==============================================================================================================================

  SUBROUTINE explicitsnow_xios_initialize

    CHARACTER(LEN=255) :: filename, name

    IF (printlev>=3) WRITE(numout,*) 'In explicitsnow_xios_initialize'  

      !! 12. Prepare for reading of relief_initsnowpack file
    filename = 'relief_initsnowpack.nc'
    name = filename(1:LEN_TRIM(FILENAME)-3)
    IF (INIT_SNOWPACK) THEN
       CALL xios_orchidee_set_file_attr("relief_file",name=name)
    ELSE
       ! Deactivate relief reading
       IF (printlev>=2) WRITE(numout,*) 'The relief file will not be read by XIOS'
       CALL xios_orchidee_set_file_attr("relief_file",name=name,enabled=.FALSE.)
       CALL xios_orchidee_set_field_attr("relief",name=name,enabled=.FALSE.)
    END IF

    !!
    !! 13. Prepare for reading of ice_sheet_mask file
    !!
    filename = 'ice_sheet_mask.nc'
    name = filename(1:LEN_TRIM(FILENAME)-3)
    IF (read_ice_sheet_mask) THEN
       CALL xios_orchidee_set_file_attr("ice_sheet_mask_file",name=name)
    ELSE
       ! Deactivate ice_sheet_mask reading
       IF (printlev>=2) WRITE(numout,*) 'The ice_sheet_mask file will not be read by XIOS'
       CALL xios_orchidee_set_file_attr("ice_sheet_mask_file",name=name,enabled=.FALSE.)
       CALL xios_orchidee_set_field_attr("icemask",name=name,enabled=.FALSE.)
    END IF
    
    IF (printlev>=3) WRITE(numout,*) 'End explicitsnow_xios_initialize'
   
  END SUBROUTINE explicitsnow_xios_initialize
  

!================================================================================================================================
!! SUBROUTINE   : explicitsnow_main
!!
!>\BRIEF        Main call for snow calculations
!!                
!! DESCRIPTION  : Main routine for calculation of the snow processes with explictsnow module.
!!
!! RECENT CHANGE(S) : None 
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_
!================================================================================================================================
  SUBROUTINE explicitsnow_main(kjpindex,    precip_rain,  precip_snow,   temp_air,    pb,       & ! in
                               u,           v,            temp_sol_new,  soilcap,     pgflux,   & ! in
                               frac_nobio,  totfrac_nobio,frac_snow_nobio,gtemp,                & ! in
                               lambda_snow, cgrnd_snow,   dgrnd_snow,     contfrac,   z0m,      & ! in
                               lambda_ice,  cgrnd_ice,    dgrnd_ice,      ice_sheet_mask,       & ! in  
                               vevapsno,    snow_age,     snow_nobio_age, snowrho,              & !inout
                               snowgrain,   snowdz,                                             & !inout
                               snowtemp,     snowheat,    snow,                                 & ! inout
                               temp_sol_add, icetemp,     icedz,                                & ! inout
                               snowliq,     subsnownobio, grndflux,      snowmelt,    tot_melt, & ! output
                               subsinksoil, zrainfall,    frac_snow_veg, veget,      veget_max, & ! output
                               run_off_lic, run_off_lic_frac )                                    ! output
             

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                               :: kjpindex         !! Domain size
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: precip_rain      !! Rainfall
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: precip_snow      !! Snowfall
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: temp_air         !! Air temperature
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: pb               !! Surface pressure
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: u,v              !! Horizontal wind speed
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: temp_sol_new     !! Surface temperature
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: soilcap          !! Soil capacity
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: pgflux           !! Net energy into snowpack
    REAL(r_std), DIMENSION (kjpindex,nnobio), INTENT(in)     :: frac_nobio       !! Fraction of continental ice, lakes, ...
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: totfrac_nobio    !! Total fraction of continental ice+lakes+ ...
    REAL(r_std), DIMENSION (kjpindex,nnobio), INTENT(in)     :: frac_snow_nobio  !! Snow cover fraction on non-vegeted area
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: frac_snow_veg    !! Snow cover fraction on vegetation
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: gtemp            !! First soil layer temperature 
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: lambda_snow      !! Coefficient of the linear extrapolation of surface temperature 
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT (in)     :: cgrnd_snow       !! Integration coefficient for snow numerical scheme
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT (in)     :: dgrnd_snow       !! Integration coefficient for snow numerical scheme
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: lambda_ice       !! Coefficient of the linear extrapolation of ice surface temperature
    REAL(r_std), DIMENSION (kjpindex,nice), INTENT (in)      :: cgrnd_ice        !! Integration coefficient for ice numerical scheme
    REAL(r_std), DIMENSION (kjpindex,nice), INTENT (in)      :: dgrnd_ice        !! Integration coefficient for ice numerical scheme
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)        :: veget            !! Fraction of vegetation type
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)        :: veget_max        !! Max. fraction of vegetation type
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)         :: ice_sheet_mask   !! Ice sheet mask
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)           :: contfrac         !! Fraction of continent in the grid
    REAL(r_std), DIMENSION(kjpindex), INTENT(in)             :: z0m               !! Surface roughness for momentum (m)


    !! 0.2 Output fields
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(out)     :: snowliq          !! Snow liquid content (m)
    REAL(r_std), DIMENSION (kjpindex,nnobio), INTENT(out)    :: subsnownobio     !! Sublimation of snow on other surface types (ice, lakes, ...)
    REAL(r_std), DIMENSION (kjpindex), INTENT(out)           :: grndflux         !! Net flux into soil [W/m2]
    REAL(r_std), DIMENSION (kjpindex), INTENT(out)           :: snowmelt         !! Snow melt [mm/dt_sechiba]
    REAL(r_std), DIMENSION (kjpindex), INTENT(out)           :: tot_melt         !! Total melt from ice and snow
    REAL(r_std), DIMENSION (kjpindex), INTENT(out)           :: subsinksoil      !! Excess of sublimation as a sink for the soil
    REAL(r_std), DIMENSION (kjpindex), INTENT(out)           :: zrainfall        !! Rain precipitation on snow [kg/m2] 
    REAL(r_std), DIMENSION (kjpindex), INTENT(out)           :: run_off_lic      !! Contains calving, melting and liquid precipitation on continental ice
    REAL(r_std), DIMENSION (kjpindex), INTENT(out)           :: run_off_lic_frac !! Contains cell fraction corresponding to run_off_lic


    !! 0.3 Modified fields
    REAL(r_std), DIMENSION (kjpindex), INTENT(inout)         :: vevapsno         !! Snow evaporation  @tex ($kg m^{-2}$) @endtex 
    REAL(r_std), DIMENSION (kjpindex), INTENT(inout)         :: snow_age         !! Snow age
    REAL(r_std), DIMENSION (kjpindex,nnobio), INTENT(inout)  :: snow_nobio_age   !! Snow age on ice, lakes, ...
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(inout)   :: snowrho          !! Snow density (Kg/m^3)
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(inout)   :: snowgrain        !! Snow grainsize (m)
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(inout)   :: snowdz           !! Snow layer thickness [m]
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(inout)   :: snowtemp         !! Snow temperature
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(inout)   :: snowheat         !! Snow heat content (J/m^2)
    REAL(r_std), DIMENSION (kjpindex,nice), INTENT(inout)    :: icetemp          !! Ice temperature
    REAL(r_std), DIMENSION (kjpindex,nice), INTENT(inout)    :: icedz            !! Ice layer thickness [m] 
    REAL(r_std), DIMENSION (kjpindex), INTENT(inout)         :: snow             !! Snow mass [Kg/m^2]
    REAL(r_std), DIMENSION (kjpindex), INTENT(inout)         :: temp_sol_add     !! Additional energy to melt snow for snow ablation case (K)

 
    !! 0.4 Local declaration
    INTEGER(i_std)                                           :: ji, iv, jj,m,jv
    REAL(r_std),DIMENSION  (kjpindex)                        :: snow_depth_tmp
    REAL(r_std),DIMENSION  (kjpindex)                        :: snowmelt_from_maxmass
    REAL(r_std)                                              :: snowdzm1
    REAL(r_std), DIMENSION (kjpindex)                        :: thrufal          !! Water leaving snow pack [kg/m2/s]
    REAL(r_std), DIMENSION (kjpindex)                        :: snowmelt_tmp,temp_sol_new_old
    REAL(r_std), DIMENSION (kjpindex,nsnow)                  :: snowdz_old
    REAL(r_std), DIMENSION (kjpindex)                        :: ZLIQHEATXS
    REAL(r_std)                                              :: grndflux_tmp
    REAL(r_std), DIMENSION (nsnow)                           :: snowtemp_tmp
    REAL(r_std)                                              :: s2flux_tmp,fromsoilflux
    REAL(r_std), DIMENSION (kjpindex,nsnow)                  :: pcapa_snow 
    REAL(r_std), DIMENSION (kjpindex)                        :: psnowhmass
    REAL(r_std), PARAMETER                                   :: XP00 = 1.E5
    REAL(r_std), DIMENSION (kjpindex)                        :: ice_sheet_melt   !! Ice melt [mm/dt_sechiba]
    REAL(r_std), DIMENSION (kjpindex)                        :: snowheattot_begin,snowheattot_end,del_snowheattot !! [J/m2]
    REAL(r_std), DIMENSION (kjpindex,nsnow)                  :: snowvisco        !! Snow viscosity [Pa/s]
    REAL(r_std), DIMENSION (kjpindex,nsnow)                  :: snowliq_diag     !! Snow liquid content per snow layer [kg/m2]
    REAL(r_std), DIMENSION (kjpindex)                        :: snowliqtot_diag  !! Snow liquid content integrated over snow depth [kg/m2]
    REAL(r_std), DIMENSION (kjpindex,nsnow)                  :: snowrho_diag
    REAL(r_std), DIMENSION (kjpindex,nsnow)                  :: snowvisco_diag
    REAL(r_std), DIMENSION (kjpindex,nsnow)                  :: snowheat_diag
    REAL(r_std), DIMENSION (kjpindex,nsnow)                  :: snowgrain_diag
    


    !! 1. Initialization
    
    temp_sol_new_old = temp_sol_new
    tot_melt(:) = zero
    snowmelt(:) = zero
    snowmelt_from_maxmass(:) = zero
    ice_sheet_melt(:) = zero 
    
    !! 2. on Vegetation
    ! 2.1 Snow fall
    snowheattot_begin(:)=SUM(snowheat(:,:),2)
    CALL explicitsnow_fall(kjpindex,precip_snow,temp_air,u,v,snowrho,snowdz,ice_sheet_mask,&
         & z0m,snowheat,snowgrain,snowtemp,psnowhmass)
    
    ! 2.2 calculate the new snow discretization
    snow_depth_tmp(:) = SUM(snowdz(:,:),2)
    
    snowdz_old = snowdz
    
   ! update snowdz
    CALL explicitsnow_levels(kjpindex,snow_depth_tmp, snowdz)

    ! 2.3 Snow heat redistribution
    CALL explicitsnow_transf(kjpindex,snowdz_old,snowdz,snowrho,snowheat,snowgrain)
    
    ! 2.4 Diagonize water portion of the snow from snow heat content:
    DO ji=1, kjpindex
       IF (SUM(snowdz(ji,:)) .GT. 0.0) THEN
          snowtemp(ji,:) = snow3ltemp_1d(snowheat(ji,:),snowrho(ji,:),snowdz(ji,:))
          snowliq(ji,:) = snow3lliq_1d(snowheat(ji,:),snowrho(ji,:),snowdz(ji,:),snowtemp(ji,:))
       ELSE
          snowliq(ji,:) = zero
          snowtemp(ji,:) = tp_00
       ENDIF
    END DO
    
    ! 2.5 snow compaction
    SELECTCASE (compact_snow_method)
    CASE (COMPACT_SNOW_METHOD_DEFAULT)
       ! original scheme
       CALL explicitsnow_compactn(kjpindex,snowtemp,snowrho,snowdz,snowvisco,snowliq)
    CASE (COMPACT_SNOW_DECHARME16)
       ! snow compaction resulting from changes in snow viscosity : compaction Decharme 2016
       CALL explicitsnow_compactn_up(kjpindex,snowtemp,snowrho,snowdz,snowliq)
       ! snow compaction due to snowdrift (used with explicitsnow_compactn_up)
       CALL explicitsnow_drift(kjpindex,u,v,snowrho,snowdz)
    CASE (COMPACT_SNOW_KAMPENHOUT)
       ! snow compaction due to snowdrift Vionnet et al 2012 + Kampenhout 2017  
       CALL explicitsnow_compactn(kjpindex,snowtemp,snowrho,snowdz,snowvisco,snowliq)    
       CALL WindDriftCompaction(kjpindex,u,v,snowrho,snowdz) 
    CASE DEFAULT
       CALL ipslerr_p(3,'explicitsnow_main','Unsupported COMPACT_SNOW_METHOD', &
               'Currently supported methods are ','default(1), Decharme16(2) or Kampenhout(3)')
    ENDSELECT
    
    ! Add XIOS default value for snowvisco_diag where no snow
    DO ji=1,kjpindex
       IF (snow(ji) .GT. zero) THEN
          snowvisco_diag(ji,:) = snowvisco(ji,:)
       ELSE
          snowvisco_diag(ji,:) = xios_default_val
       END IF
    END DO
    
    ! Update snow heat 
    DO ji = 1, kjpindex
       snowheat(ji,:) = snow3lheat_1d(snowliq(ji,:),snowrho(ji,:),snowdz(ji,:),snowtemp(ji,:))
    ENDDO

    ! 2.6 Calculate the snow temperature profile based on heat diffusion
    CALL explicitsnow_profile (kjpindex,cgrnd_snow,dgrnd_snow,lambda_snow,temp_sol_new, snowtemp,snowdz,temp_sol_add)
    
    ! 2.7 Test whether snow is existed on the ground or not
    grndflux(:)=0.0
    CALL explicitsnow_gone(kjpindex,pgflux,&
         & snowheat,snowtemp,snowdz,snowrho,snowliq,grndflux,snowmelt)
    
    ! 2.8 Calculate snow melt processes
    CALL explicitsnow_melt(kjpindex,pcapa_snow,precip_rain,totfrac_nobio,frac_snow_nobio, &
     frac_snow_veg,veget,veget_max,zrainfall,snowtemp,snowdz,snowrho,snowliq,snowmelt,grndflux)
    
    ! 2.9 Calculate snow refreezing processes
    CALL explicitsnow_refrz(kjpindex,zrainfall,snowtemp,snowdz,snowrho,snowliq,snowmelt)
    
    IF (ok_ice_sheet) THEN
       ! Calculate ice temperature   
       CALL explicitsnow_iceprofile(kjpindex, cgrnd_ice,dgrnd_ice,lambda_ice,temp_sol_new,icedz, &
            & ice_sheet_mask, snowdz, cgrnd_snow, dgrnd_snow, snowtemp, temp_sol_add, icetemp)
         
       ! Calculate ice melt
       CALL explicitsnow_icemelt(kjpindex,ice_sheet_mask,snowdz,precip_snow, zrainfall, snowmelt, vevapsno, icetemp,icedz, &
            & grndflux,ice_sheet_melt)

       ! Update icetemp and icedz
       CALL explicitsnow_icelevels(kjpindex,ice_sheet_mask,ice_sheet_melt,icetemp,icedz)
    ENDIF  
    
    
    ! 2.10 Snow sublimation changing snow thickness
    CALL explicitsnow_subli(kjpindex,snowrho,frac_nobio,snowdz,vevapsno, &
            snowliq,snowtemp,snow,subsnownobio,subsinksoil)
    
    ! 2.11 Snowmelt_from_maxmass 
    CALL explicitsnow_maxmass(kjpindex,snowrho,soilcap,snow,snowdz,snowmelt_from_maxmass)

        
    ! 2.12 Calculate snow grain size using the updated thermal gradient
    !      This subroutine is not used in current version
    !    CALL explicitsnow_grain(kjpindex,snowliq,snowdz,gtemp,snowtemp,pb,snowgrain)
    
    ! 2.13  Update snow heat
    ! Update the heat content (variable stored each time step)
    ! using current snow temperature and liquid water content:
    !
    ! First, make check to make sure heat content not too large
    ! (this can result due to signifigant heating of thin snowpacks):
    ! add any excess heat to ground flux:
    !
    DO ji=1,kjpindex
       DO jj=1,nsnow
          ZLIQHEATXS(ji)  = MAX(0.0, snowliq(ji,jj)*ph2o - 0.10*snowdz(ji,jj)*snowrho(ji,jj))*chalfu0/dt_sechiba
          snowliq(ji,jj) = snowliq(ji,jj) - ZLIQHEATXS(ji)*dt_sechiba/(ph2o*chalfu0)
          snowliq(ji,jj) = MAX(0.0, snowliq(ji,jj))
          grndflux(ji)   = grndflux(ji)   + ZLIQHEATXS(ji)
       ENDDO
    ENDDO
    
    snow(:) = 0.0
    DO ji=1,kjpindex !domain size
       snow(ji) = SUM(snowrho(ji,:) * snowdz(ji,:))
    ENDDO
        
    DO ji = 1, kjpindex
       snowheat(ji,:) = snow3lheat_1d(snowliq(ji,:),snowrho(ji,:),snowdz(ji,:),snowtemp(ji,:))
    ENDDO
    snowheattot_end(:)=SUM(snowheat(:,:),2)
    del_snowheattot(:)=snowheattot_end(:)-snowheattot_begin(:)
    
    
    !! 4. On other surface types - not done yet
    
    IF ( nnobio .GT. 1 ) THEN
       WRITE(numout,*) 'WE HAVE',nnobio-1,' SURFACE TYPES I DO NOT KNOW'
       WRITE(numout,*) 'CANNOT TREAT SNOW ON THESE SURFACE TYPES'
       CALL ipslerr_p(3,'explicitsnow_main','Unsupported surface type','','')
    ENDIF
    
    !! 5. Computes snow age on land and land ice (for albedo)
    call explicitsnow_age(kjpindex,snow,precip_snow,precip_rain,frac_snow_nobio, &
            temp_sol_new,snow_age,snow_nobio_age)
    
  
    !! 6. Check the snow on land 
    DO ji=1,kjpindex
       IF (snow(ji) .EQ. 0) THEN
          snowrho(ji,:)=50.0
          snowgrain(ji,:)=0.0
          snowdz(ji,:)=0.0
          snowliq(ji,:)=0.0
       ENDIF
    ENDDO
    
    
    ! Snow melt only if there is more than a given mass : maxmass_snow
    ! Here I suggest to remove the snow based on a certain threshold of snow
    ! depth instead of snow mass because it is quite difficult for
    ! explictsnow module to calculate other snow properties following the
    ! removal of snow mass
    ! to define the threshold of snow depth based on old snow density (330
    ! kg/m3)
    !       maxmass_snowdepth=maxmass_snow/sn_dens 
!    snowmelt_from_maxmass(:) = 0.0
    
    !! 7. compute total melt 
    
    DO ji=1,kjpindex 
       tot_melt(ji) = snowmelt(ji) +  snowmelt_from_maxmass(ji)
    ENDDO


    !! Calculate variables for coupling with LMDZ. 
    run_off_lic(:) =  ice_sheet_melt(:) + (snowmelt(:) + snowmelt_from_maxmass(:))*frac_nobio(:,iice)
    run_off_lic_frac(:) = frac_nobio(:,iice)


    !! 8. Recalculate the new snow discretization
    !     This is done to make sure that all levels have snow
    snow_depth_tmp(:) = SUM(snowdz(:,:),2)
    snowdz_old = snowdz
    CALL explicitsnow_levels(kjpindex,snow_depth_tmp, snowdz)
    CALL explicitsnow_transf(kjpindex,snowdz_old,snowdz,snowrho,snowheat,snowgrain)

    ! Add XIOS default value where no snow
    DO ji=1,kjpindex
       IF (snow(ji) .GT. zero) THEN
          ! Output for snowliq and snowliqtot change unit from m to kg/m2 and are multiplied by contfrac 
          ! to take into acount the portion form land fraction only
          snowliq_diag(ji,:) = snowliq(ji,:) * 1000 * contfrac(ji)
          snowliqtot_diag(ji) = SUM(snowliq_diag(ji,:)) * contfrac(ji)
          snowrho_diag(ji,:) = snowrho(ji,:)
          snowheat_diag(ji,:) = snowheat(ji,:)
          snowgrain_diag(ji,:) = snowgrain(ji,:)
       ELSE
          snowliq_diag(ji,:) = xios_default_val
          snowliqtot_diag(ji) = xios_default_val
          snowrho_diag(ji,:) = xios_default_val
          snowheat_diag(ji,:) = xios_default_val
          snowgrain_diag(ji,:) = xios_default_val
       END IF
    END DO

    CALL xios_orchidee_send_field("snowliq",snowliq_diag)
    CALL xios_orchidee_send_field("snowliqtot", snowliqtot_diag)
    CALL xios_orchidee_send_field("snowrho",snowrho_diag)
    CALL xios_orchidee_send_field("snow_viscosity",snowvisco_diag)
    CALL xios_orchidee_send_field("snowheat",snowheat_diag)
    CALL xios_orchidee_send_field("snowgrain",snowgrain_diag)
    CALL xios_orchidee_send_field("snowmelt_from_maxmass",snowmelt_from_maxmass/dt_sechiba)
    CALL xios_orchidee_send_field("soilcap",soilcap)
    CALL xios_orchidee_send_field("del_snowheattot",del_snowheattot)
    CALL xios_orchidee_send_field("run_off_lic",run_off_lic/dt_sechiba)
    CALL xios_orchidee_send_field("run_off_lic_frac",run_off_lic_frac)

    IF (printlev>=3) WRITE(numout,*) 'explicitsnow_main done'

  END SUBROUTINE explicitsnow_main
  

!================================================================================================================================
!! SUBROUTINE   : explicitsnow_finalize
!!
!>\BRIEF        Write variables for explictsnow module to restart file
!!                
!! DESCRIPTION  : Write variables for explictsnow module to restart file
!!
!! \n
!_
!================================================================================================================================
  SUBROUTINE explicitsnow_finalize ( kjit,     kjpindex, rest_id,    snowrho,   &
       snowtemp, snowdz, snowheat, snowgrain, icetemp)
    
    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                           :: kjit           !! Time step number 
    INTEGER(i_std), INTENT(in)                           :: kjpindex       !! Domain size
    INTEGER(i_std),INTENT (in)                           :: rest_id        !! Restart file identifier
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(in)  :: snowrho        !! Snow density (Kg/m^3)
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(in)  :: snowtemp       !! Snow temperature (K)
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(in)  :: snowdz         !! Snow layer thickness [m]
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(in)  :: snowheat       !! Snow heat content (J/m^2)
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(in)  :: snowgrain      !! Snow grainsize (m)
    REAL(r_std), DIMENSION (kjpindex,nice), INTENT(in)   :: icetemp        !! Ice temperature (K)
    
    !! 1. Write to restart file
    CALL restput_p(rest_id, 'snowrho', nbp_glo, nsnow, 1, kjit, snowrho, 'scatter', nbp_glo, index_g)
    !-
    CALL restput_p(rest_id, 'snowtemp', nbp_glo, nsnow, 1, kjit, snowtemp, 'scatter', nbp_glo, index_g)
    !-
    CALL restput_p(rest_id, 'snowdz', nbp_glo, nsnow, 1, kjit, snowdz, 'scatter', nbp_glo, index_g)
    !-
    CALL restput_p(rest_id, 'snowheat', nbp_glo, nsnow, 1, kjit, snowheat, 'scatter', nbp_glo, index_g)
    !-
    CALL restput_p(rest_id, 'snowgrain', nbp_glo, nsnow, 1, kjit, snowgrain, 'scatter', nbp_glo, index_g)

    IF ( ok_ice_sheet) CALL restput_p(rest_id, 'icetemp', nbp_glo, nice, 1, kjit, icetemp, 'scatter', nbp_glo, index_g)
    
  END SUBROUTINE explicitsnow_finalize
  

!================================================================================================================================
!! SUBROUTINE   : explicitsnow_grain
!!
!>\BRIEF        Compute evolution of snow grain size
!!                
!! DESCRIPTION  : 
!!
!! RECENT CHANGE(S) : None 
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : R. Jordan (1991)
!!
!! FLOWCHART    : None
!! \n
!_
!================================================================================================================================


  SUBROUTINE explicitsnow_grain(kjpindex,snowliq,snowdz,gtemp,snowtemp,pb,snowgrain)

    !! 0.1 Input variables
    INTEGER(i_std),INTENT(in)                                  :: kjpindex         !! Domain size
    REAL(r_std),DIMENSION(kjpindex,nsnow),INTENT(in)           :: snowliq          !! Liquid water content
    REAL(r_std),DIMENSION(kjpindex,nsnow),INTENT(in)           :: snowdz           !! Snow depth (m)
    REAL(r_std),DIMENSION(kjpindex),INTENT(in)                 :: gtemp            !! First soil layer temperature
    REAL(r_std),DIMENSION(kjpindex,nsnow),INTENT(in)           :: snowtemp         !! Snow temperature (K)
    REAL(r_std),DIMENSION (kjpindex),INTENT(in)                :: pb               !! Surface pressure (hpa) 

    !! 0.2 Output variables

    !! 0.3 Modified variables
    REAL(r_std),DIMENSION(kjpindex,nsnow),INTENT(inout)        :: snowgrain        !! Snow grain size (m)

    !! 0.4 Local variables
    REAL(r_std),DIMENSION(kjpindex,nsnow)                      :: zsnowdz,zdz,ztheta
    REAL(r_std),DIMENSION(kjpindex,0:nsnow)                    :: ztemp,zdiff,ztgrad,zthetaa,zfrac,&
         zexpo,zckt_liq,zckt_ice,zckt
    REAL(r_std),DIMENSION(kjpindex,nsnow)                      :: zrhomin,zgrainmin
    INTEGER(i_std) :: ji,jj

    !! 0.5 Local parameters
    REAL(r_std), PARAMETER                                     :: ztheta_crit = 0.02     !! m3 m-3
    REAL(r_std), PARAMETER                                     :: zc1_ice     = 8.047E+9 !! kg m-3 K
    REAL(r_std), PARAMETER                                     :: zc1_liq     = 5.726E+8 !! kg m-3 K
    REAL(r_std), PARAMETER                                     :: zdeos       = 0.92E-4  !! effective diffusion
    !! coef for water vapor in snow
    !! at 0C and 1000 mb (m2 s-1)
    REAL(r_std), PARAMETER                                     :: zg1         = 5.0E-7   !! m4 kg-1
    REAL(r_std), PARAMETER                                     :: zg2         = 4.0E-12  !! m2 s-1
    REAL(r_std), PARAMETER                                     :: ztheta_w      = 0.05   !! m3 m-3
    REAL(r_std), PARAMETER                                     :: ztheta_crit_w = 0.14   !! m3 m-3
    REAL(r_std), PARAMETER                                     :: zdzmin        = 0.01   !! m : minimum thickness
    !! for thermal gradient evaluation:
    !! to prevent excessive gradients
    !! for vanishingly thin snowpacks.
    REAL(r_std), PARAMETER                                     :: xp00=1.E5

    !! 1. initialize

    DO ji=1,kjpindex 


       zsnowdz(ji,:)  = MAX(xsnowdmin/nsnow, snowdz(ji,:))

       DO jj=1,nsnow-1
          zdz(ji,jj)      = zsnowdz(ji,jj) + zsnowdz(ji,jj+1)
       ENDDO
       zdz(ji,nsnow)     = zsnowdz(ji,nsnow)

       ! compute interface average volumetric water content (m3 m-3):
       ! first, layer avg VWC:
       !
       ztheta(ji,:) = snowliq(ji,:)/MAX(xsnowdmin, zsnowdz(ji,:))

       ! at interfaces:
       zthetaa(ji,0)      = ztheta(ji,1)
       DO jj=1,nsnow-1
          zthetaa(ji,jj)  = (zsnowdz(ji,jj)  *ztheta(ji,jj)   +             &
               zsnowdz(ji,jj+1)*ztheta(ji,jj+1))/zdz(ji,jj)
       ENDDO
       zthetaa(ji,nsnow) = ztheta(ji,nsnow)
       ! compute interface average temperatures (K):
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
       !
       ztemp(ji,0)      = snowtemp(ji,1)
       DO jj=1,nsnow-1
          ztemp(ji,jj)  = (zsnowdz(ji,jj)  *snowtemp(ji,jj)   +             &
               zsnowdz(ji,jj+1)*snowtemp(ji,jj+1))/zdz(ji,jj)
       ENDDO
       ztemp(ji,nsnow) = snowtemp(ji,nsnow)
       !
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
       ! compute variation of saturation vapor pressure with temperature
       ! for solid and liquid phases:
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
       zexpo(ji,:)    = chalsu0/(xrv*ztemp(ji,:))
       zckt_ice(ji,:) = (zc1_ice/ztemp(ji,:)**2)*(zexpo(ji,:) - 1.0)*EXP(-zexpo(ji,:))
       !
       zexpo(ji,:)    = chalev0/(xrv*ztemp(ji,:))
       zckt_liq(ji,:) = (zc1_liq/ztemp(ji,:)**2)*(zexpo(ji,:) - 1.0)*EXP(-zexpo(ji,:))
       !
       ! compute the weighted ice/liquid total variation (N m-2 K):
       !
       zfrac(ji,:)    = MIN(1.0, zthetaa(ji,:)/ztheta_crit)
       zckt(ji,:)     = zfrac(ji,:)*zckt_liq(ji,:) + (1.0 - zfrac(ji,:))*zckt_ice(ji,:)

       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
       ! Compute effective diffusion coefficient (m2 s-1):
       ! -diffudivity relative to a reference diffusion at 1000 mb and freezing point
       !  multiplied by phase energy coefficient
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
       !
       DO jj=0,nsnow
          zdiff(ji,jj) = zdeos*(xp00/(pb(ji)*100.))*((ztemp(ji,jj)/tp_00)**6)*zckt(ji,jj)
       ENDDO

       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
       ! Temperature gradient (K m-1):

       ztgrad(ji,0)      = 0.0 ! uppermost layer-mean and surface T's are assumed to be equal
       DO jj=1,nsnow-1
          ztgrad(ji,jj)  = 2*(snowtemp(ji,jj)-snowtemp(ji,jj+1))/MAX(zdzmin, zdz(ji,jj))
       ENDDO
       !
       ! assume at base of snow, temperature is in equilibrium with soil
       ! (but obviously must be at or below freezing point!)
       !
       ztgrad(ji,nsnow) = 2*(snowtemp(ji,nsnow) - MIN(tp_00, gtemp(ji)))/MAX(zdzmin, zdz(ji,nsnow))
       ! prognostic grain size (m) equation:
       !-------------------------------------------------------------------
       ! first compute the minimum grain size (m):
       !
       zrhomin(ji,:)     = xrhosmin
       zgrainmin(ji,:)   = snow3lgrain_1d(zrhomin(ji,:))

       ! dry snow:
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
       !
       DO jj=1,nsnow

          IF(ztheta(ji,jj) == 0.0) THEN

             ! only allow growth due to vapor flux INTO layer: 
             ! aab add sublimation(only condensation) as upper BC...?

             snowgrain(ji,jj)  = snowgrain(ji,jj) +                                      &
                  (dt_sechiba*zg1/MAX(zgrainmin(ji,jj),snowgrain(ji,jj)))*      &
                  ( zdiff(ji,jj-1)*MAX(0.0,ztgrad(ji,jj-1)) -                &
                  zdiff(ji,jj)  *MIN(0.0,ztgrad(ji,jj)) )
          ELSE

             ! wet snow
             ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
             !
             snowgrain(ji,jj)  = snowgrain(ji,jj) +                                      &
                  (dt_sechiba*zg2/MAX(zgrainmin(ji,jj),snowgrain(ji,jj)))*      &
                  MIN(ztheta_crit_w, ztheta(ji,jj) + ztheta_w)
          END IF

       ENDDO


    ENDDO


  END SUBROUTINE explicitsnow_grain

!================================================================================================================================
!! SUBROUTINE   : explicitsnow_compactn
!!
!>\BRIEF        Compute Compaction/Settling
!!                
!! DESCRIPTION  : 
!!     Snow compaction due to overburden and settling.
!!     Mass is unchanged: layer thickness is reduced
!!     in proportion to density increases. Method
!!     of Anderson (1976): see Loth and Graf, 1993,
!!     J. of Geophys. Res., 98, 10,451-10,464.
!!
!! RECENT CHANGE(S) : None 
!!
!! MAIN OUTPUT VARIABLE(S): snowrho, snowdz
!!
!! REFERENCE(S) : Loth and Graf (1993), Mellor (1964) and Anderson (1976)
!!
!! FLOWCHART    : None
!! \n
!_
!================================================================================================================================


  SUBROUTINE explicitsnow_compactn(kjpindex,snowtemp,snowrho,snowdz,snowvisco,snowliq)

    !! 0.1 Input variables
    INTEGER(i_std),INTENT(in)                                 :: kjpindex         !! Domain size
    REAL(r_std),DIMENSION(kjpindex,nsnow),INTENT(in)          :: snowtemp         !! Snow temperature
    REAL(r_std),DIMENSION(kjpindex,nsnow),INTENT(in)          :: snowliq        !! Liquid water content

    !! 0.2 Output variables
    REAL(r_std),DIMENSION(kjpindex,nsnow),INTENT(out)         :: snowvisco        !! Snow viscosity (Pa/s)

    !! 0.3 Modified variables
    REAL(r_std),DIMENSION(kjpindex,nsnow),INTENT(inout)       :: snowrho          !! Snow density (Kg/m^3)
    REAL(r_std),DIMENSION(kjpindex,nsnow),INTENT(inout)       :: snowdz           !! Snow depth [m]

    !! 0.4 Local variables
    REAL(r_std),DIMENSION(kjpindex,nsnow)                     :: zwsnowdz,zsmass,zsnowrho2,zviscocity,zsettle 
    REAL(r_std),DIMENSION(kjpindex)                           :: zsmassc          !! cummulative snow mass (kg/m2)
    REAL(r_std),DIMENSION(kjpindex)                           :: snowdepth_crit
    INTEGER(i_std)                                            :: ji,jj

    !! 1. initialize

    zsnowrho2  = snowrho
    zsettle(:,:)    = ZSNOWCMPCT_ACM
    zviscocity(:,:) = ZSNOWCMPCT_V0

    !! 2. Calculating Cumulative snow mass (kg/m2):
    
    DO ji=1, kjpindex

       
       IF (SUM(snowdz(ji,:)) .GT. 0.0) THEN 

          zwsnowdz(ji,:)= snowdz(ji,:)*snowrho(ji,:)

          zsmassc (ji)= 0.0

          DO jj=1,nsnow
             zsmass(ji,jj)   = zsmassc(ji) + zwsnowdz(ji,jj)
             zsmassc(ji)      = zsmassc(ji) + zwsnowdz(ji,jj)
          ENDDO


          !! 3. Computing compaction/Settling
          ! ----------------------
          ! Compaction/settling if density below upper limit
          ! (compaction is generally quite small above ~ 500 kg m-3):
          !
          DO jj=1,nsnow
             IF (snowrho(ji,jj) .LT. xrhosmax) THEN

                ! First calculate settling due to freshly fallen snow: (NOTE:bug here for the snow temperature profile)
                !

                zsettle(ji,jj)     = ZSNOWCMPCT_ACM*EXP(                                      &
                     -ZSNOWCMPCT_BCM*(tp_00-MIN(tp_00,snowtemp(ji,jj)))                      &
                     -ZSNOWCMPCT_CCM*MAX(0.0,                                  &
                     snowrho(ji,jj)-ZSNOWCMPCT_RHOD))
                !
                ! Snow viscocity:
                !              

                !!Vicosity formulation :
                !! VISCOSITY_METHOD=1 Anderson et al., 1979 formulation
                !! VISCOSITY_METHOD =2 Vionnet et al., 2012 formulation
                SELECTCASE (VISCOSITY_METHOD)
                CASE (VISCOSITY_METHOD_DEFAULT) !! Anderson et al., 1979
                   zviscocity(ji,jj)   = ZSNOWCMPCT_V0*EXP( ZSNOWCMPCT_VT*(tp_00-MIN(tp_00,snowtemp(ji,jj))) + &
                        ZSNOWCMPCT_VR*snowrho(ji,jj) )

                   ! Calculate snow density: compaction from weight/over-burden
                   ! Anderson 1976 method:
                   zsnowrho2(ji,jj)    = snowrho(ji,jj) + snowrho(ji,jj)*dt_sechiba*  &
                        ((cte_grav*zsmass(ji,jj)/zviscocity(ji,jj))                   &
                        + zsettle(ji,jj) )
                   
                CASE (VISCOSITY_METHOD_VIONNET) !! Vionnet (2012) viscosity
                   zviscocity(ji,jj) =  (1./(1.+60.*(snowliq(ji,jj)/1000*snowdz(ji,jj))))* &
                        MIN(4.,exp(MIN(0.4,gs_vionnet-0.2)/0.1))*eta_0*(snowrho(ji,jj)/c_eta)* &
                        exp(a_eta*(tp_00-snowtemp(ji,jj))+b_eta*snowrho(ji,jj))
                        
                   zsnowrho2(ji,jj)    = snowrho(ji,jj) + snowrho(ji,jj)*dt_sechiba*       &
                        ((cte_grav*zsmass(ji,jj)/zviscocity(ji,jj)))
                   
                CASE DEFAULT
                   CALL ipslerr_p(3,'explicitsnow_compactn','Unsupported VISCOSITY_METHOD', &
                        'Currently supported methods are ','default(1) or vionnet(2)')
                ENDSELECT
                
                ! Conserve mass by decreasing grid thicknesses in response
                ! to density increases
                !
                snowdz(ji,jj)  = snowdz(ji,jj)*(snowrho(ji,jj)/zsnowrho2(ji,jj))

             ENDIF
          ENDDO

          ! Update density (kg m-3):
          snowrho(ji,:) = zsnowrho2(ji,:)
          snowvisco(ji,:) = zviscocity(ji,:)
          
       ENDIF
       
    ENDDO
    
  END SUBROUTINE explicitsnow_compactn

!================================================================================================================================
!! SUBROUTINE   : explicitsnow_compactn_up
!!
!>\BRIEF        Compute Compaction/Settling
!!                
!! DESCRIPTION  : 
!!     Snow compaction due to overburden and settling.
!!     Mass is unchanged: see Decharme et al. (2016)
!!
!! RECENT CHANGE(S) : 12/07/2021 
!!
!! MAIN OUTPUT VARIABLE(S): snowrho, snowdz
!!
!! REFERENCE(S) : Decharme et al. (2016) 
!!
!! FLOWCHART    : None
!! \n
!_
!================================================================================================================================


  SUBROUTINE explicitsnow_compactn_up(kjpindex,snowtemp,snowrho,snowdz,snowliq)

    !! 0.1 Input variables
    INTEGER(i_std),INTENT(in)                                 :: kjpindex
    REAL(r_std),DIMENSION(kjpindex,nsnow),INTENT(in)          :: snowtemp       !! Snow temperature
    REAL(r_std),DIMENSION(kjpindex,nsnow),INTENT(in)          :: snowliq        !! Liquid water content

    !! 0.2 Output variables

    !! 0.3 Modified variables
    REAL(r_std),DIMENSION(kjpindex,nsnow),INTENT(inout)       :: snowrho        !! Snow density (Kg/m^3)
    REAL(r_std),DIMENSION(kjpindex,nsnow),INTENT(inout)       :: snowdz         !! Snow depth

    !! 0.4 Local variables
    REAL(r_std),DIMENSION(kjpindex,nsnow)                     :: zsmass, zsnowrho2, zviscocity
    REAL(r_std),DIMENSION(kjpindex,nsnow)                     :: ztemp          !! temperature dependance for snow viscosity
    REAL(r_std),DIMENSION(kjpindex,nsnow)                     :: zwholdmax      !! max water liquid content
    REAL(r_std),DIMENSION(kjpindex,nsnow)                     :: Fwliq          !! describes the decrease of viscosity in presence of liquid water
    INTEGER(i_std)                                            :: ji,jj

    !! 1. initialize
    zsnowrho2(:,:) = snowrho(:,:)
    zsmass(:,:) = 0.0

    !! 2. Calculating Cumulative snow mass (kg/m2):

    DO ji=1, kjpindex
       IF (SUM(snowdz(ji,:)) .GT. 0.0) THEN 

          ! 1. Cumulative snow mass (kg/m2):
          DO jj=2,nsnow
             zsmass(ji,jj) = zsmass(ji,jj-1) + snowdz(ji,jj) * snowrho(ji,jj-1)
          ENDDO
          ! 1st layer : half the mass of the uppermost layer applied to itself
          zsmass(ji,1)= 0.5 * snowdz(ji,1) * snowrho(ji,1)

          ! 2. Compaction
          ! Liquid water effect
          zwholdmax(ji,:) = snow3lhold_1d(snowrho(ji,:),snowdz(ji,:))
          zwholdmax(ji,:) = max(1.e-10, zwholdmax(ji,:))
          Fwliq(ji,:) = 1./ (1. + 10.*MIN(1.0,snowliq(ji,:)/zwholdmax(ji,:)))

          ! Snow viscocity, density and grid thicknesses
          DO jj=1,nsnow
             IF (snowrho(ji,jj) .LT. xrhosmax) THEN
                ! Temperature dependence limited to 5K: Schleef et al. (2014)
                ztemp(ji,jj) = ZSNOWCMPCT_XT * MIN(DTref,tp_00-MIN(tp_00,snowtemp(ji,jj)))

                ! Calculate snow viscocity: Brun et al. (1989), Vionnet et al. (2012)
                zviscocity(ji,jj)   = ZSNOWCMPCT_X0 * Fwliq(ji,jj) * EXP(ztemp(ji,jj) &
                     + ZSNOWCMPCT_XR * snowrho(ji,jj)) * snowrho(ji,jj)/xrhosref
                ! Calculate new snow density:                  
                zsnowrho2(ji,jj) = snowrho(ji,jj) + dt_sechiba &
                     *(snowrho(ji,jj)*cte_grav*zsmass(ji,jj)/zviscocity(ji,jj))
                ! Conserve mass by decreasing grid thicknesses in response to density increases
                snowdz(ji,jj)  = snowdz(ji,jj)*(snowrho(ji,jj)/zsnowrho2(ji,jj))
             ENDIF
          ENDDO
          ! Update density (kg m-3):
          snowrho(ji,:) = zsnowrho2(ji,:)
       ENDIF
    ENDDO

  END SUBROUTINE explicitsnow_compactn_up
      
!!
!================================================================================================================================
!! SUBROUTINE   : explicitsnow_drift
!!
!>\BRIEF        Compute Compaction due to snow drift : wind induced densification of near-surface snow layers
!!                
!! DESCRIPTION  : 
!!     Snow compaction due to snowdrift
!!     Mass is unchanged: see Decharme et al. (2016)
!!
!! RECENT CHANGE(S) : 04/10/2021 
!!
!! MAIN OUTPUT VARIABLE(S): snowrho, snowdz
!!
!! REFERENCE(S) : Decharme et al. (2016)
!!
!! FLOWCHART    : None
!! \n
!_
!================================================================================================================================


  SUBROUTINE explicitsnow_drift(kjpindex,u,v,snowrho,snowdz)

    !! 0.1 Input variables
    INTEGER(i_std),INTENT(in)                                 :: kjpindex  
    REAL(r_std),DIMENSION(kjpindex),INTENT(in)                :: u,v

    !! 0.2 Output variables

    !! 0.3 Modified variables
    REAL(r_std),DIMENSION(kjpindex,nsnow),INTENT(inout)       :: snowrho          !! Snow density (Kg/m^3)
    REAL(r_std),DIMENSION(kjpindex,nsnow),INTENT(inout)       :: snowdz           !! Snow depth

    !! 0.4 Local variables
    REAL(r_std),DIMENSION(kjpindex,nsnow)                     :: zsnowrho2, zsettle
    REAL(r_std),DIMENSION(kjpindex)                           :: zsmobc      !! cumulative mobility index
    REAL(r_std),DIMENSION(kjpindex,nsnow)                     :: zsmob       !! mobility index
    REAL(r_std),DIMENSION(kjpindex,nsnow)                     :: zwmob       !! for mobility index
    REAL(r_std),DIMENSION(kjpindex,nsnow)                     :: Zmob        !! Mobility index
    REAL(r_std),DIMENSION(kjpindex,nsnow)                     :: Zwind       !! Wind-driven compaction index
    REAL(r_std),DIMENSION(kjpindex,nsnow)                     :: Ztau        !! Function used for the calculation of teh compaction rate
    REAL(r_std),DIMENSION(kjpindex,nsnow)                     :: tau_cmpct   !! Compaction rate (s)
    REAL(r_std)                                               :: speed       !! Wind speed
    INTEGER(i_std)                                            :: ji,jj
    LOGICAL                                                   :: Zwindup     !! True if Zwind(jj-1) > 0

    !! 1. initialize
    zsnowrho2(:,:) = snowrho(:,:)
    zsettle(:,:) = 0.
    zsmobc(:) = 0.0

    !! 2. Calculating Cumulative snow mass (kg/m2):

    DO ji=1, kjpindex
       Zwindup = .true.
       IF (SUM(snowdz(ji,:)) .GT. 0.0) THEN 

          ! Computation of zsmob : Analogous to zsmass (used in Equation B3 Decharme 2016)
          DO jj=1,nsnow
             ! Mobility index : Equation B1 dans Appendix B Decharme 2016
             Zmob(ji,jj) = Amob*(1-MAX(0.,(snowrho(ji,jj)-xrhosmin)/xrhosmob))  ! Amob= 1.25, xrhosmin= 50 kg.m-3, xrhosmob= 295 kg.m-3
             zwmob(ji,jj) = snowdz(ji,jj)*(b_tau - Zmob(ji,jj))
             zsmob(ji,jj) = zsmobc(ji) + zwmob(ji,jj)
             zsmobc(ji) = zsmobc(ji) + zwmob(ji,jj)
          ENDDO

          ! Snow viscocity, density and grid thicknesses
          DO jj=1,nsnow
             IF (snowrho(ji,jj) .LT. xrhosmax) THEN
                ! Calculate wind densification : Brun 1997
                ! Mobility index : Equation B1 dans Appendix B Decharme 2016
                Zmob(ji,jj) = Amob*(1-MAX(0.,(snowrho(ji,jj)-xrhosmin)/xrhosmob))  ! Amob= 1.25, xrhosmin= 50 kg.m-3, xrhosmob= 295 kg.m-3
                ! Computation of zsmob : Analogous to zsmass (used in Equation B3 Decharme 2016)
                zwmob(ji,jj) = snowdz(ji,jj)*(b_tau - Zmob(ji,jj))
                zsmob(ji,jj) = zsmobc(ji) + zwmob(ji,jj)
                zsmobc(ji) = zsmobc(ji) + zwmob(ji,jj)

                speed=SQRT(u(ji)*u(ji)+v(ji)*v(ji)) ! wind speed

                ! Wind driven compaction index : Equation B2 in Appendix B Decharme 2016                     
                Zwind(ji,jj) = 1 - ACMPCT*EXP(-BCMPCT*KCMPCT*speed) + Zmob(ji,jj) ! ACMPCT=2.868, BCMPCT=0.085 s.m-1, KCMPCT=1.25

                IF ((Zwind(ji,jj).GT.0.).AND.Zwindup) THEN ! snowdrift occurs only if Zwind>0 
                   Ztau(ji,jj) = MAX(0., Zwind(ji,jj))*EXP(a_tau*zsmob(ji,jj))    

                   ! Wind compaction rate (Equation B3 in Appendix B)
                   tau_cmpct(ji,jj) = 2*KCMPCT*one_day/Ztau(ji,jj)

                   ! zsettle : Second term in the right-hand side of Equation 13 in Decharme et al. 2016
                   zsettle(ji,jj) = MAX(0.,(xrhoswind-snowrho(ji,jj))/tau_cmpct(ji,jj)) !xrhoswind=350
                ELSE
                   Zwindup = .false.
                   Ztau(ji,jj) = 0.
                   tau_cmpct(ji,jj) = 0.
                   zsettle(ji,jj) = 0.
                ENDIF

                ! Calculate new snow density:                  
                ! settling by wind transport only in case of not too dense snow
                IF (snowrho(ji,jj).LT.xrhoswind) THEN
                   ! settling by wind cannot lead to densities above xrhoswind                     
                   zsnowrho2(ji,jj) = min(xrhoswind, snowrho(ji,jj) + dt_sechiba * zsettle(ji,jj))
                   ! Conserve mass by decreasing grid thicknesses in response to density increases
                   snowdz(ji,jj) = snowdz(ji,jj)*(snowrho(ji,jj)/zsnowrho2(ji,jj))
                ENDIF
             ENDIF
          ENDDO
          ! Update density (kg m-3):
          snowrho(ji,:) = zsnowrho2(ji,:)
       ENDIF
    ENDDO

  END SUBROUTINE explicitsnow_drift
  
!================================================================================================================================
!>\BRIEF        Compute Compaction due to snow drift : wind induced densification of near-surface snow layers
!!                
!! DESCRIPTION  : 
!!     Snow compaction due to snowdrift
!!     Mass is unchanged: see Decharme et al. (2016)
!!
!! RECENT CHANGE(S) : 04/10/2021 
!!
!! MAIN OUTPUT VARIABLE(S): snowrho, snowdz
!!
!! REFERENCE(S) : Decharme et al. (2016)
!!
!! FLOWCHART    : None
!! \n
!_
!================================================================================================================================


  

subroutine WindDriftCompaction(kjpindex,u,v,snowrho,snowdz)
    !
    ! !DESCRIPTION:
    !
    ! Compute wind drift compaction for a single column and level.
    !
    ! Also updates zpseudo and mobile for this column. However, zpseudo remains unchanged
    ! if mobile is already false or becomes false within this subroutine.
    !
    ! The structure of the updates done here for zpseudo and mobile requires that this
    ! subroutine be called first for the top layer of snow, then for the 2nd layer down,
    ! etc. - and finally for the bottom layer. Before beginning the loops over layers,
    ! mobile should be initialized to .true. and zpseudo should be initialized to 0.
    !
    ! !USES:
    !

    !! 0.1 Input variables
    INTEGER(i_std),INTENT(in)                                 :: kjpindex
    REAL(r_std),DIMENSION(kjpindex),INTENT(in)                :: u,v

    !! 0.2 Output variables

    !! 0.3 Modified variables
    REAL(r_std),DIMENSION(kjpindex,nsnow),INTENT(inout)       :: snowrho          !! Snow density (Kg/m^3)
    REAL(r_std),DIMENSION(kjpindex,nsnow),INTENT(inout)       :: snowdz           !! Snow depth

    !! 0.4 Local variables
    REAL(r_std),DIMENSION(kjpindex,nsnow)                     :: zsnowrho2, compaction_rate
    REAL(r_std),DIMENSION(kjpindex)                           :: zpseudoc      ! wind drift compaction / pseudo depth cumulative
    REAL(r_std),DIMENSION(kjpindex,nsnow)                     :: zpseudo      ! wind drift compaction / pseudo depth for this column at this layer
    REAL(r_std),DIMENSION(kjpindex,nsnow)                     :: Frho       ! Mobility density factor [-]
    REAL(r_std),DIMENSION(kjpindex,nsnow)                     :: MO          ! Mobility index [-]
    REAL(r_std),DIMENSION(kjpindex,nsnow)                     :: SI         ! Driftability index [-]
    REAL(r_std),DIMENSION(kjpindex,nsnow)                     :: gamma_drift        ! Scaling factor for wind drift time scale [-]
    REAL(r_std),DIMENSION(kjpindex,nsnow)                     :: tau_inverse   ! Inverse of the effective time scale [1/s]
    REAL(r_std)                                               :: speed       !! Wind speed
    INTEGER(i_std)                                            :: ji,jj
    LOGICAL                                                   :: mobile      !! mobile boolean

    REAL(r_std)                                               :: rho_min = 50.0      ! wind drift compaction / minimum density [kg/m3]
    REAL(r_std)                                               :: rho_max = 350.0     ! wind drift compaction / maximum density [kg/m3]
    REAL(r_std)                                               :: tau_ref = 48.0 * 3600.0  ! wind drift compaction / reference time [s]


    !-----------------------------------------------------------------------

    !! 1. initialize
    zsnowrho2(:,:) = snowrho(:,:)
    compaction_rate(:,:) = 0.
    zpseudoc(:) = 0.0

    !! 2. Calculating Wind Drif Compaction:

    DO ji=1, kjpindex
       IF (SUM(snowdz(ji,:)) .GT. 0.0) THEN
          mobile=.TRUE.
          DO jj=1,nsnow
          IF (snowrho(ji,jj) .LT. xrhosmax) THEN
             IF (mobile) THEN
                Frho(ji,jj) = 1.25 - 0.0042*(max(rho_min, snowrho(ji,jj))-rho_min)
                ! assuming dendricity = 0, sphericity = 1, grain size = 0.35 mm Non-dendritic snow
                MO(ji,jj) = -0.069 + 0.66*Frho(ji,jj)
                speed=SQRT(u(ji)*u(ji)+v(ji)*v(ji)) ! wind speed
                SI(ji,jj) = -2.868 * exp(-0.085*speed) + 1.0 + MO(ji,jj)

                     IF (SI(ji,jj) > 0.0) THEN
                           SI(ji,jj) = min(SI(ji,jj), 3.25)
                          ! Increase zpseudo (wind drift / pseudo depth) to the middle of
                        ! the pseudo-node for the sake of the following calculation
                        zpseudo(ji,jj) = zpseudoc(ji) + 0.5 * snowdz(ji,jj) * (3.25 - SI(ji,jj))
                        zpseudoc(ji) = zpseudoc(ji) + zpseudo(ji,jj)

                        gamma_drift(ji,jj) = SI(ji,jj)*exp(-zpseudo(ji,jj)/0.1)
                        tau_inverse(ji,jj) = gamma_drift(ji,jj) / tau_ref
                        compaction_rate(ji,jj) = max(0.0, rho_max-snowrho(ji,jj)) * tau_inverse(ji,jj)
                      ELSE  ! SI <= 0
                        mobile = .false.
                        compaction_rate(ji,jj) = 0.0
                     END IF
            END IF
            IF (snowrho(ji,jj).LT.xrhoswind) THEN
            ! Calculate new snow density:                 
            ! settling by wind cannot lead to densities above xrhoswind                     
            zsnowrho2(ji,jj) = snowrho(ji,jj) + dt_sechiba * compaction_rate(ji,jj)
            ! Conserve mass by decreasing grid thicknesses in response to density increases
            snowdz(ji,jj) = snowdz(ji,jj)*(snowrho(ji,jj)/zsnowrho2(ji,jj))
            END IF
            END IF
         END DO
         ! Update density (kg m-3):
         snowrho(ji,:) = zsnowrho2(ji,:)
      END IF
             
   END DO

  end subroutine WindDriftCompaction






!!
!================================================================================================================================
!! SUBROUTINE   : explicitsnow_transf
!!
!>\BRIEF        Computing snow mass and heat redistribution due to grid thickness configuration resetting
!!                
!! DESCRIPTION  : Snow mass and heat redistibution due to grid thickness
!!                configuration resetting. Total mass and heat content
!!                of the overall snowpack unchanged/conserved within this routine. 
!! RECENT CHANGE(S) : None 
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_
!================================================================================================================================

  SUBROUTINE explicitsnow_transf(kjpindex,snowdz_old,snowdz,snowrho,snowheat,snowgrain)

    !! 0.1 Input variables
    INTEGER(i_std),INTENT(in)                                           :: kjpindex         !! Domain size
    REAL(r_std),DIMENSION(kjpindex,nsnow),INTENT(in)                    :: snowdz_old       !! Snow depth at the previous time step 

    !! 0.2 Output variables

    !! 0.3 Modified variables
    REAL(r_std),DIMENSION(kjpindex,nsnow),INTENT(inout)                 :: snowrho          !! Snow density (Kg/m^3)
    REAL(r_std),DIMENSION(kjpindex,nsnow),INTENT(inout)                 :: snowgrain        !! Snow grain size (m)
    REAL(r_std),DIMENSION(kjpindex,nsnow),INTENT(inout)                 :: snowdz           !! Snow depth (m)
    REAL(r_std),DIMENSION(kjpindex,nsnow),INTENT(inout)                 :: snowheat         !! Snow heat content/enthalpy (J/m2)
                                                                        
    !! 0.4 Local varibles                                              
    REAL(r_std),DIMENSION(kjpindex,nsnow)                               :: zsnowrhon        !! snow density (new)
    REAL(r_std),DIMENSION(kjpindex,nsnow)                               :: zsnowgrainn      !! snow grain size (new)
    REAL(r_std),DIMENSION(kjpindex,nsnow)                               :: zsnowgraino      !! snow grain size (old)
    REAL(r_std),DIMENSION(kjpindex)                                     :: zsumheat         !! snow heat (z sum)
    REAL(r_std),DIMENSION(kjpindex)                                     :: zsumswe          !! snow water equivalent (z sum)
    REAL(r_std),DIMENSION(kjpindex)                                     :: zsumgrain        !! snow grain size (z sum)
    REAL(r_std),DIMENSION(kjpindex)                                     :: psnow_new        !! total snow height (new)
    REAL(r_std),DIMENSION(kjpindex)                                     :: psnow_old        !! total snow height (old)
    REAL(r_std),DIMENSION(kjpindex,nsnow)                               :: zsnowrhoo        !! snow density (old)
    REAL(r_std),DIMENSION(kjpindex,nsnow)                               :: zsnowheato       !! snow heat (old)
    REAL(r_std),DIMENSION(kjpindex,nsnow)                               :: zsnowheatn       !! snow heat (new)
    REAL(r_std),DIMENSION(kjpindex,nsnow)                               :: zmassdzo         !! snow mass (old) intercepting new layers
    REAL(r_std),DIMENSION(kjpindex,nsnow)                               :: zsnowztop_old    !! vertical top grid limits (old)
    REAL(r_std),DIMENSION(kjpindex,nsnow)                               :: zsnowztop_new    !! vertical top grid limits (new)
    REAL(r_std),DIMENSION(kjpindex,nsnow)                               :: zsnowzbot_old    !! vertical bottom grid limits (old)
    REAL(r_std),DIMENSION(kjpindex,nsnow)                               :: zsnowzbot_new    !! vertical bottom grid limits (new)
    REAL(r_std),DIMENSION(kjpindex,nsnow)                               :: zmasstotn        !! snow mass (new)
    REAL(r_std)                                                         :: ZPROPOR          !! proportional coef
    INTEGER(i_std)                                                      :: ji,jj,jjo
    

    ! Initialization
    zsumheat(:)       = 0.0
    zsumswe(:)        = 0.0
    zsumgrain(:)      = 0.0
    zsnowheatn(:,:)   = 0.0
    zsnowrhon(:,:)    = 0.0
    zmasstotn(:,:)    = 0.0
    zsnowgrainn(:,:)  = 0.0
    zsnowgraino(:,:)  = 0.0
    
    
    DO ji=1, kjpindex
       psnow_new(ji) = SUM(snowdz(ji,:))
       IF (psnow_new(ji) > 0.) THEN
          IF (psnow_new(ji) .GE. xsnowcritd .AND. ALL(snowdz_old(ji,:).NE.0.)) THEN
             ! total snow height of old snow
             psnow_old(ji) = SUM(snowdz_old(ji,:))

             ! initialization of variables describing the initial snowpack 
             zsnowrhoo(ji,:) = snowrho(ji,:)
             zsnowheato(ji,:) = snowheat(ji,:)
             zmassdzo(ji,:) = undef_sechiba

             ! 1. Calculate vertical grid limits (m):
             ! --------------------------------------
             zsnowztop_old(ji,1) = psnow_old(ji)
             zsnowztop_new(ji,1) = psnow_new(ji)
             zsnowzbot_old(ji,1) = zsnowztop_old(ji,1) - snowdz_old(ji,1)
             zsnowzbot_new(ji,1) = zsnowztop_new(ji,1) - snowdz(ji,1)

             ! calculate vertical position of every snow layer limits
             DO jj=2,nsnow
                zsnowztop_old(ji,jj) = zsnowzbot_old(ji,jj-1)
                zsnowztop_new(ji,jj) = zsnowzbot_new(ji,jj-1)
                zsnowzbot_old(ji,jj) = zsnowztop_old(ji,jj) - snowdz_old(ji,jj)
                zsnowzbot_new(ji,jj) = zsnowztop_new(ji,jj) - snowdz(ji,jj)
             ENDDO

             zsnowzbot_old(ji,nsnow) = 0.0
             zsnowzbot_new(ji,nsnow) = 0.0

             ! 3. Calculate mass, heat, charcateristics mixing due to vertical grid resizing:
             ! --------------------------------------------------------------------
             !
             ! loop over the new snow layers
             ! Summ or avergage of the constituting quantities of the old snow layers
             ! which are totally or partially inserted in the new snow layer
             ! For snow age, mass weighted average is used.

             DO jj=1,nsnow
                DO jjo=1,nsnow
                   IF((zsnowztop_old(ji,jjo)>zsnowzbot_new(ji,jj)).AND.(zsnowzbot_old(ji,jjo)<zsnowztop_new(ji,jj)))THEN                
                      ZPROPOR = (MIN(zsnowztop_old(ji,jjo), zsnowztop_new(ji,jj)) &
                           -  MAX(zsnowzbot_old(ji,jjo), zsnowzbot_new(ji,jj))) &
                           / snowdz_old(ji,jjo)
                      zmassdzo(ji,jjo) = zsnowrhoo(ji,jjo)*snowdz_old(ji,jjo)*ZPROPOR
                      zmasstotn(ji,jj) = zmasstotn(ji,jj) + zmassdzo(ji,jjo)
                      zsnowheatn(ji,jj) = zsnowheatn(ji,jj) + zsnowheato(ji,jjo)*ZPROPOR
                      zsnowgraino(ji,jjo) = snowgrain(ji,jjo)*snowdz_old(ji,jjo)*ZPROPOR
                      zsnowgrainn(ji,jj) = zsnowgrainn(ji,jj) + zsnowgraino(ji,jjo)
                   ENDIF
                ENDDO
             ENDDO

             ! the new layer inherits from the weighted average properties of the old ones
             ! heat and mass
             !   zsnowheatn(ji,:) = ZSNOWHEAN(ji,:)  ! remplace par le calcul direct de zsnowheatn !!!!
             zsnowrhon(ji,:) = zmasstotn(ji,:)/snowdz(ji,:)
             snowrho(ji,:)   = zsnowrhon(ji,:)
             snowheat(ji,:)  = zsnowheatn(ji,:)
             snowgrain(ji,:)  = zsnowgrainn(ji,:)/snowdz(ji,:)

          ENDIF

          ! 4. Vanishing or very thin snowpack check:
          ! -----------------------------------------
          !
          ! NOTE: ONLY for very shallow snowpacks, mix properties (homogeneous):
          ! this avoids problems related to heat and mass exchange for
          ! thin layers during heavy snowfall or signifigant melt: one
          ! new/old layer can exceed the thickness of several old/new layers.
          ! Therefore, mix (conservative):
          IF(psnow_new(ji) < xsnowcritd .OR. ANY(snowdz_old(ji,:).EQ. 0.)) THEN
             zsumheat(ji) = SUM(snowheat(ji,:))
             zsumswe(ji)  = SUM(snowrho(ji,:)*snowdz_old(ji,:))
             zsumgrain(ji)= SUM(snowgrain(ji,:)*snowdz_old(ji,:))
             DO jj=1,nsnow
                zsnowheatn(ji,jj)  = zsumheat(ji)/nsnow
                snowdz(ji,jj)      = psnow_new(ji)/nsnow
                zsnowrhon(ji,jj)   = zsumswe(ji)/psnow_new(ji)
                zsnowgrainn(ji,jj) = zsumgrain(ji)/psnow_new(ji)
             ENDDO
             snowrho(ji,:)   = zsnowrhon(ji,:)
             snowheat(ji,:)  = zsnowheatn(ji,:)
             snowgrain(ji,:) = zsnowgrainn(ji,:)
          ENDIF
       ENDIF
    ENDDO
   

  END SUBROUTINE explicitsnow_transf

  
!!
!================================================================================================================================
!! SUBROUTINE   : explicitsnow_fall
!!
!>\BRIEF    Computes snowfall    
!!                
!! DESCRIPTION  : Computes snowfall    
!routine. 
!! RECENT CHANGE(S) : None 
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_
!================================================================================================================================
 
  SUBROUTINE explicitsnow_fall(kjpindex,precip_snow,temp_air,u,v,snowrho,snowdz,ice_sheet_mask,&
                        & z0m,snowheat,snowgrain,snowtemp,psnowhmass)

    !! 0.1 Input variables
    INTEGER(i_std),INTENT(in)                              :: kjpindex            !! Domain size
    REAL(r_std),DIMENSION(kjpindex),INTENT(in)             :: precip_snow         !! Snow rate (SWE) (kg/m2 per dt_sechiba)
    REAL(r_std),DIMENSION(kjpindex),INTENT(in)             :: temp_air            !! Air temperature
    REAL(r_std),DIMENSION(kjpindex),INTENT(in)             :: u,v                 !! Horizontal wind speed
    REAL(r_std),DIMENSION(kjpindex,nsnow),INTENT(in)       :: snowtemp            !! Snow temperature
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)       :: ice_sheet_mask      !! Ice sheet mask
    REAL(r_std), DIMENSION(kjpindex), INTENT(in)           :: z0m               !! Surface roughness for momentum (m)
    
    !! 0.2 Output variables
    REAL(r_std), DIMENSION(kjpindex),INTENT(out)           :: psnowhmass          !! Heat content of snowfall (J/m2)

    !! 0.3 Modified variables
    REAL(r_std),DIMENSION(kjpindex,nsnow),INTENT(inout)    :: snowrho             !! Snow density profile (kg/m3)
    REAL(r_std),DIMENSION(kjpindex,nsnow),INTENT(inout)    :: snowdz              !! Snow layer thickness profile (m)
    REAL(r_std),DIMENSION(kjpindex,nsnow),INTENT(inout)    :: snowheat            !! Snow heat content/enthalpy (J/m2)
    REAL(r_std),DIMENSION(kjpindex,nsnow),INTENT(inout)    :: snowgrain           !! Snow grain size (m)

    !! 0.4 Local variables
    REAL(r_std), DIMENSION(kjpindex)                       :: rhosnew             !! Snowfall density
    REAL(r_std), DIMENSION(kjpindex)                       :: dsnowfall           !! Snowfall thickness (m)
    REAL(r_std), DIMENSION(kjpindex,nsnow)                 :: snowdz_old
    REAL(r_std), DIMENSION(kjpindex)                       :: snow_depth
    REAL(r_std), DIMENSION(kjpindex)                       :: snow_depth_old
    REAL(r_std), DIMENSION(kjpindex)                       :: newgrain
    REAL(r_std)                                            :: snowfall_delta
    REAL(r_std)                                            :: speed
    INTEGER(i_std)                                         :: ji,jj
    REAL(r_std)                                            :: rhos_max
    REAL(r_std)                                            :: ws_cut
    REAL(r_std)                                            :: intervalle_u
    REAL(r_std), DIMENSION(kjpindex)                       :: speed_10m
    REAL(r_std)                                            :: rhos_min_tanh
    REAL(r_std)                                            :: rhos_max_tanh
    REAL(r_std), DIMENSION(kjpindex)                       :: z                   !! height in m of the forced wind
    REAL(r_std), DIMENSION(kjpindex)                       :: u_10m, v_10m
    REAL(r_std), PARAMETER                                 :: SIGMA = 0.999747931957245 !! ORCHIDEE wind vertical level in sigma coordinate (2m)
    
    !! 1. initialize the variables
    rhosnew(:) = 0.0
    rhos_max = xrhosmax
    snowdz_old = snowdz
    DO ji=1,kjpindex
      snow_depth(ji) = SUM(snowdz(ji,:))
    ENDDO

    snow_depth_old = snow_depth
    snowfall_delta = 0.0 

    !! 2. incorporate snowfall into snowpack 
    DO ji = 1, kjpindex
       
      z(ji) = ((1. - SIGMA)/9.81)*287.058*temp_air(ji)    
      u_10m(ji)  = u(ji) * (LOG(10.0/z0m(ji))/LOG(z(ji)/z0m(ji)))
      v_10m(ji)  = v(ji) * (LOG(10.0/z0m(ji))/LOG(z(ji)/z0m(ji)))
      speed_10m(ji) = MAX(min_wind, SQRT (u_10m(ji)*u_10m(ji) + v_10m(ji)*v_10m(ji)))
      
      speed = MAX(min_wind, SQRT (u(ji)*u(ji) + v(ji)*v(ji)))
        
      ! new snow fall on snowpack 
      ! NOTE: when the surface temperature is zero, it means that snowfall has no
      ! heat content but it can be used for increasing the thickness and changing the density (maybe it is a bug)
      psnowhmass(ji) = 0.0
      IF ( (precip_snow(ji) .GT. 0.0) ) THEN

         ! calculate
         psnowhmass(ji) = precip_snow(ji)*(xci*(temp_air(ji)-tp_00)-chalfu0)
                
         IF (SURFDENSITY_XRHOSMIN) THEN
            ! Snowfall density with a fixed value
            rhosnew(ji) = xrhosmin
         ELSE IF (RHOS_WS) THEN
            ! Snowfall density computed as surfacedensity from initialisation
            rhos_min_tanh = 260.
            rhos_max_tanh = 450.
            ws_cut = 8.5
            intervalle_u = 4.
            IF (ice_sheet_mask(ji).EQ.1) rhos_max=400. ! Greenland
            
            rhosnew(ji) = rhos_min_tanh + ( rhos_max_tanh - rhos_min_tanh ) * &
                 ( TANH( (speed_10m(ji) - ws_cut ) / intervalle_u ) + 1. ) / 2.
            
            rhosnew(ji) = min(rhosnew(ji), rhos_max)
         ELSE      
            ! Snowfall density: Following CROCUS (Pahaut 1976)
            rhosnew(ji) = MAX(xrhosmin, snowfall_a_sn + snowfall_b_sn*(temp_air(ji)-tp_00) + &
                 snowfall_c_sn*SQRT(speed))
         END IF

         ! Augment total pack depth:
         dsnowfall(ji) = precip_snow(ji)/rhosnew(ji) !snowfall thickness (m)
         snow_depth(ji) = snow_depth(ji) + dsnowfall(ji) 

         ! Fresh snowfall changes the snowpack density and liquid content in uppermost layer 
         IF (dsnowfall(ji) .GT. zero) THEN
            snowrho(ji,1) = (snowdz(ji,1)*snowrho(ji,1) + dsnowfall(ji)*rhosnew(ji))/ &
                 (snowdz(ji,1)+dsnowfall(ji))
            
            snowdz(ji,1) = snowdz(ji,1) + dsnowfall(ji)

            ! Add energy of snowfall to snowpack:
            ! Update heat content (J/m2) (therefore the snow temperature
            ! and liquid content):
            snowheat(ji,1)  = snowheat(ji,1) + psnowhmass(ji)
            
            ! Incorporate snowfall grain size:
            newgrain(ji)    = MIN(dgrain_new_max, snow3lgrain_0d(rhosnew(ji)))
            snowgrain(ji,1) = (snowdz_old(ji,1)*snowgrain(ji,1) + dsnowfall(ji)*newgrain(ji))/ &
                 snowdz(ji,1)
         ENDIF
      ELSE
        dsnowfall(ji) = 0.
      ENDIF

      ! new snow fall on snow free surface. 
      ! we use the linearization for the new snow fall on snow-free ground 
      IF ( (dsnowfall(ji) .GT. zero) .AND. (snow_depth_old(ji) .EQ. zero) ) THEN
        snowfall_delta = 1.0
        DO jj=1,nsnow
          snowdz(ji,jj) = snowfall_delta*(dsnowfall(ji)/nsnow) + &
                          (1.0-snowfall_delta)*snowdz(ji,jj)

          snowheat(ji,jj) = snowfall_delta*(psnowhmass(ji)/nsnow) + &
                            (1.0-snowfall_delta)*snowheat(ji,jj)

          snowrho(ji,jj) = snowfall_delta*rhosnew(ji) + &
                          (1.0-snowfall_delta)*snowrho(ji,jj)

          snowgrain(ji,jj) = snowfall_delta*newgrain(ji) + &
                            (1.0-snowfall_delta)*snowgrain(ji,jj)
        ENDDO
      ENDIF
    ENDDO 

  END SUBROUTINE explicitsnow_fall

!!
!================================================================================================================================
!! SUBROUTINE   : explicitsnow_gone
!!
!>\BRIEF        Check whether snow is gone 
!!                
!! DESCRIPTION  : If so, set thickness (and therefore mass and heat) and liquid
!!                content to zero, and adjust fluxes of water, evaporation and
!!                heat into underlying surface. 
!! RECENT CHANGE(S) : None 
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_
!================================================================================================================================

  SUBROUTINE explicitsnow_gone(kjpindex,pgflux,&
       snowheat,snowtemp,snowdz,snowrho,snowliq,grndflux,snowmelt)

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                                 :: kjpindex     !! Domain size
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)              :: pgflux       !! Net energy into snow pack(w/m2)
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT(in)          :: snowheat     !! Snow heat content (J/m^2)

    !! 0.2 Output variables

    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(inout)     :: snowtemp     !! Snow temperature
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(inout)     :: snowdz       !! Snow depth [m]
    REAL(r_std),DIMENSION (kjpindex,nsnow), INTENT(inout)      :: snowrho      !! Snow density (Kg/m^3)
    REAL(r_std),DIMENSION (kjpindex,nsnow), INTENT(inout)      :: snowliq      !! Liquid water content
    REAL(r_std),DIMENSION(kjpindex), INTENT(inout)             :: grndflux     !! Soil/snow interface heat flux (W/m2)
    REAL(r_std),DIMENSION(kjpindex),INTENT(inout)              :: snowmelt     !! Snow melt
    REAL(r_std),DIMENSION(kjpindex)                            :: thrufal      !! Water leaving snowpack(kg/m2/s)

    !! 0.4 Local variables

    INTEGER(i_std)                                             :: ji,jj
    REAL(r_std),DIMENSION(kjpindex)                            :: snowgone_delta
    REAL(r_std),DIMENSION (kjpindex)                           :: totsnowheat  !!snow heat content at each layer 
    REAL(r_std),DIMENSION(kjpindex)                            :: snowdepth_crit

    ! first caculate total snowpack snow heat content
    snowgone_delta(:) = un
    thrufal(:)=0.0
    snowmelt(:)=0
    totsnowheat(:)  = SUM(snowheat(:,:),2) 

    DO ji = 1, kjpindex

       IF ( pgflux(ji) >= (-totsnowheat(ji)/dt_sechiba) ) THEN
          ! all the snow melts
          grndflux(ji) = pgflux(ji) + (totsnowheat(ji)/dt_sechiba)
          thrufal(ji)=SUM(snowrho(ji,:)*snowdz(ji,:))
          snowgone_delta(ji) = 0.0       
          snowmelt(ji) = snowmelt(ji)+thrufal(ji)
       ENDIF

       ! update of snow state (either still present or not)
       DO jj=1,nsnow
          snowdz(ji,jj)  =   snowdz(ji,jj) *snowgone_delta(ji)
          snowliq(ji,jj)   =   snowliq(ji,jj) *snowgone_delta(ji)
          snowtemp(ji,jj) = (1.0-snowgone_delta(ji))*tp_00 + snowtemp(ji,jj)*snowgone_delta(ji)
       ENDDO
    ENDDO

  END SUBROUTINE explicitsnow_gone


!================================================================================================================================
!! SUBROUTINE   : explicitsnow_melt
!!
!>\BRIEF        Computes snow melt processes within snowpack
!!                
!! DESCRIPTION  : 
!! RECENT CHANGE(S) : None 
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_
!================================================================================================================================

  SUBROUTINE explicitsnow_melt(kjpindex,pcapa_snow,precip_rain,totfrac_nobio,frac_snow_nobio, &
       frac_snow_veg,veget,veget_max,zrainfall,snowtemp,snowdz,snowrho,snowliq,snowmelt,grndflux)

    !! 0.1 Input variables
    INTEGER(i_std), INTENT (in)                             :: kjpindex        !! Domain size
    REAL(r_std),DIMENSION (kjpindex,nsnow)                  :: pcapa_snow      !! Heat capacity for snow
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)            :: precip_rain     !! Rainfall
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)           :: totfrac_nobio   !! Total fraction of continental ice+lakes+ ...
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT(in)     :: frac_snow_nobio !! Snow cover fraction on non-vegeted area
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)            :: frac_snow_veg   !! Snow cover fraction on vegetation
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)       :: veget           !! Fraction of vegetation type
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)       :: veget_max       !! Max. fraction of vegetation type
   
    !! 0.2 Output variables
    REAL(r_std), DIMENSION (kjpindex), INTENT(out)          :: zrainfall       !! Rain precipitation on snow (mm)
    
    !! 0.3 Modified variables
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(inout)  :: snowtemp     !! Snow temperature 
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(inout)  :: snowdz       !! Snow depth [m]
    REAL(r_std),DIMENSION (kjpindex,nsnow), INTENT(inout)   :: snowrho      !! Snow layer density (Kg/m^3)
    REAL(r_std),DIMENSION (kjpindex,nsnow), INTENT(inout)   :: snowliq      !! Liquid water content
    REAL(r_std),DIMENSION (kjpindex), INTENT(inout)         :: snowmelt     !! Snowmelt
    REAL(r_std),DIMENSION(kjpindex),INTENT(inout)           :: grndflux     !! Net energy input to soil 

    !! 0.4 Local variables
    REAL(r_std),DIMENSION (nsnow)                           :: zmeltxs      !! Excess melt
    REAL(r_std),DIMENSION (kjpindex)                        :: meltxs       !! Residual snowmelt energy applied to underlying soil
    REAL(r_std),DIMENSION (nsnow)                           :: zsnowlwe     !! total liquid equivalent water content of each snow layer (m)
    REAL(r_std),DIMENSION (kjpindex)                        :: snowmass     !! total mass of snow
    REAL(r_std),DIMENSION (nsnow)                           :: zphase       !! Phase change (from ice to water) (J/m2)
    REAL(r_std),DIMENSION (nsnow)                           :: zsnowmelt    !! Snow melt (liquid water) (m)
    REAL(r_std),DIMENSION (nsnow)                           :: zsnowtemp

    REAL(r_std),DIMENSION (nsnow)                           :: zwholdmax    !! Maximum liquid water holding (m)
    REAL(r_std),DIMENSION (nsnow)                           :: zcmprsfact   !! Compression factor due to densification from melting
    REAL(r_std),DIMENSION(kjpindex,nsnow)                   :: melt_snow    !! Total snowmelt (before refreezing) (m)
    REAL(r_std),DIMENSION (kjpindex)                        :: frac_rain_veg 
    INTEGER(i_std)                                          :: ji,jv


    ! Initialize
    melt_snow(:,:) = 0.
    zrainfall(:) = 0.
    
    DO ji = 1, kjpindex
       ! rain fraction in snow computed as in precisol
       frac_rain_veg(ji) = veget_max(ji,1)*frac_snow_veg(ji)
       DO jv = 2,nvm
          frac_rain_veg(ji) = frac_snow_veg(ji)*throughfall_by_pft(jv)*veget(ji,jv) &
               + frac_snow_veg(ji)*(veget_max(ji,jv) - veget(ji,jv))  &
               + frac_rain_veg(ji)
       ENDDO
       
       snowmass(ji) = SUM(snowrho(ji,:) * snowdz(ji,:))
       IF ((snowmass(ji) .GT. min_sechiba)) THEN

          !! 1 snow melting due to positive snowpack snow temperature

          !! 1.0 total liquid equivalent water content of each snow layer

          zsnowlwe(:) = snowrho(ji,:) * snowdz(ji,:)/ ph2o 

          !! 1.1 phase change (J/m2)

          pcapa_snow(ji,:) = snowrho(ji,:)*xci

          zphase(:)  = MIN(pcapa_snow(ji,:)*MAX(0.0, snowtemp(ji,:)-tp_00)*      &
               snowdz(ji,:),                                       &
               MAX(0.0,zsnowlwe(:)-snowliq(ji,:))*chalfu0*ph2o)

          !! 1.2 update snow liq water and temperature if melting

          zsnowmelt(:) = zphase(:)/(chalfu0*ph2o)

          !! 1.3 cool off the snow layer temperature due to melt

          zsnowtemp(:) = snowtemp(ji,:) - zphase(:)/(pcapa_snow(ji,:)* snowdz(ji,:))

          snowtemp(ji,:) = MIN(tp_00, zsnowtemp(:))

          zmeltxs(:)   = (zsnowtemp(:)-snowtemp(ji,:))*pcapa_snow(ji,:)*snowdz(ji,:)

          !! 1.4 loss of snowpack depth and liquid equivalent water

          zwholdmax(:) = snow3lhold_1d(snowrho(ji,:),snowdz(ji,:)) ! 1 dimension

          zcmprsfact(:) = (zsnowlwe(:)-MIN(snowliq(ji,:)+zsnowmelt(:),zwholdmax(:)))/ &
               (zsnowlwe(:)-MIN(snowliq(ji,:),zwholdmax(:)))

          snowdz(ji,:)    = snowdz(ji,:)*zcmprsfact(:)

          snowrho(ji,:)     = zsnowlwe(:)*ph2o/snowdz(ji,:)

          snowliq(ji,:)   = snowliq(ji,:) + zsnowmelt(:)

          melt_snow(ji,:) = zsnowmelt(:)

          ! excess heat from melting, using it to warm underlying ground to conserve energy
          meltxs(ji) = SUM(zmeltxs(:))/dt_sechiba  ! (W/m2)

          ! energy flux into the soil 
          grndflux(ji) = grndflux(ji) + meltxs(ji)
          
          ! Rain entering snow (m):
          zrainfall(ji) = (frac_snow_nobio(ji,iice)*totfrac_nobio(ji) &
               + frac_rain_veg(ji) )* precip_rain(ji)

       ELSE
          snowdz(ji,:)=0.
          snowliq(ji,:)=0.
          snowmelt(ji)=snowmelt(ji)+SUM(snowrho(ji,:)*snowdz(ji,:))
          !This addition is to get the precipitation that falls on the snow in case the snow has just totally disppeared in explicitsnow_gone.
          snowmelt(ji) = snowmelt(ji) + (frac_snow_nobio(ji,iice)*totfrac_nobio(ji) &
               + frac_rain_veg(ji) )*precip_rain(ji)
       ENDIF

    ENDDO
    
    CALL xios_orchidee_send_field("melt_snow",melt_snow)
    
  END SUBROUTINE explicitsnow_melt



!================================================================================================================================
!! SUBROUTINE   : explicitsnow_refrz
!!
!>\BRIEF        Computes snow refreezing processes within snowpack
!!                
!! DESCRIPTION  : 
!! RECENT CHANGE(S) : None 
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_
!================================================================================================================================

  SUBROUTINE explicitsnow_refrz(kjpindex,zrainfall,snowtemp,snowdz,snowrho,snowliq,snowmelt)

    !! 0.1 Input variables
    INTEGER(i_std), INTENT (in)                             :: kjpindex        !! Domain size
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)           :: zrainfall       !! Rain precipitation on snow (mm)

    !! 0.2 Output variables

    !! @tex $(kg m^{-2})$ @endtex
    !! 0.3 Modified variables
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(inout)  :: snowtemp     !! Snow temperature 
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(inout)  :: snowdz       !! Snow depth [m]
    REAL(r_std),DIMENSION (kjpindex,nsnow), INTENT(inout)   :: snowrho      !! Snow layer density (Kg/m^3)
    REAL(r_std),DIMENSION (kjpindex,nsnow), INTENT(inout)   :: snowliq      !! Liquid water content
    REAL(r_std),DIMENSION (kjpindex), INTENT(inout)         :: snowmelt     !! Snowmelt

    !! 0.4 Local variables
    REAL(r_std),DIMENSION (nsnow)                           :: flowliq
    REAL(r_std),DIMENSION (kjpindex)                        :: snowmass
    REAL(r_std),DIMENSION (nsnow)                           :: zphase       !! Phase change (from water to ice) (J/m2)
    REAL(r_std),DIMENSION (nsnow)                           :: zsnowdz      !! Snow layer depth [m]
    REAL(r_std),DIMENSION (nsnow)                           :: zwholdmax    !! Maximum liquid water holding (m)
    REAL(r_std),DIMENSION (nsnow)                           :: zscap        !! Snow heat capacity (J/m3 K)
    REAL(r_std),DIMENSION (nsnow)                           :: zsnowliq     !! (m)
    REAL(r_std),DIMENSION (nsnow)                           :: snowtemp_old
    REAL(r_std),DIMENSION (0:nsnow)                         :: zflowliqt    !!(m)
    REAL(r_std)                                             :: zpcpxs
    REAL(r_std)                                             :: ztotwcap
    INTEGER(i_std)                                          :: jj,ji
    REAL(r_std),DIMENSION(kjpindex,nsnow)                   :: zrefrz       !! (mm)
    REAL(r_std),DIMENSION(nsnow)                            :: zsnowrho
    REAL(r_std)                                             :: znumer

    ! Initialize
    zrefrz(:,:) = 0. 

    DO ji = 1, kjpindex

       snowmass(ji) = SUM(snowrho(ji,:) * snowdz(ji,:))
       IF ((snowmass(ji) .GT. min_sechiba)) THEN

          ! Refreezing is limited by available energy as in CROCUS
          !
          ! 1. Increases Liquid Water of top layer from rain
          !    ---------------------------------------------
          !
          !  Rainfall (m) initialises the liquid flow which feeds the top layer         
          zflowliqt(0) = zrainfall(ji)/ph2o ! zflowliqt on top layer = zrainfall
          snowtemp_old(:) = snowtemp(ji,:)  ! keep in memory snowtemp before refreezing
          zwholdmax(:) = snow3lhold_1d(snowrho(ji,:),snowdz(ji,:)) ! 1 dimension
          zscap(:) = snow3lscap_1d(snowrho(ji,:)) !     = snowrho(ji,:)*xci
          DO jj=1,nsnow ! loop on snow layers
             !
             ! 2. Increases Liquid Water from the upper layers flow (or rain for top layer)
             !    -----------------------------
             zsnowliq(jj) = snowliq(ji,jj) + zflowliqt(jj-1) ! add liquid water coming from top layer
             !
             ! Calculate the maximum possible refreezing
             zphase(jj) = min(zscap(jj)* max(0.0, tp_00 - snowtemp_old(jj))*snowdz(ji,jj),  &
               zsnowliq(jj)*chalfu0*ph2o)
             !
             ! warm layer and reduce liquid if freezing occurs
             zsnowdz(jj) = MAX(xsnowdmin/nsnow, snowdz(ji,jj))
             snowtemp(ji,jj) = snowtemp(ji,jj) + zphase(jj)/(zscap(jj)*zsnowdz(jj))
             !
             ! Reduce liquid portion if freezing occurs:
             zsnowliq(jj) = zsnowliq(jj) - ( (snowtemp(ji,jj)-snowtemp_old(jj))*       &
                  zscap(jj)*zsnowdz(jj)/(chalfu0*ph2o) )
             !! Refreezing:
             zrefrz(ji,jj) = (snowtemp(ji,jj)-snowtemp_old(jj))*zscap(jj)*zsnowdz(jj)/chalfu0
             zsnowliq(jj) = MAX(zsnowliq(jj), 0.0)
             !
             ! 4. Calculate flow from the excess of holding capacity
             !    --------------------------------------------------------------
             !
             ! Any water in excess of the maximum holding space for liquid water
             ! amount is drained into next layer down.
             zflowliqt(jj) = max(0.0,zsnowliq(jj)-zwholdmax(jj))
             zsnowliq(jj)  = zsnowliq(jj) - zflowliqt(jj) ! zsnowliq limited to zwholdmax
             !
             ! 5. Density is adjusted to conserve the mass
             !    --------------------------------------------------------------
             flowliq(jj) = zflowliqt(jj) - zflowliqt(jj-1)
             IF (flowliq(jj).GT.0.) THEN ! mass loss
                ! update snowdz but avoid 0 values
                zsnowdz(jj) = MAX(xsnowdmin/nsnow, snowdz(ji,jj) - flowliq(jj)*ph2o/snowrho(ji,jj))
                zsnowrho(jj) = (snowrho(ji,jj) * snowdz(ji,jj) - flowliq(jj)*ph2o) / zsnowdz(jj)
             ELSE ! mass gain
                znumer =  snowrho(ji,jj) * snowdz(ji,jj) - flowliq(jj) * ph2o
                zsnowrho(jj) = znumer / snowdz(ji,jj)
                zsnowdz(jj) = snowdz(ji,jj)
             ENDIF
             snowliq(ji,jj) = zsnowliq(jj)
             snowrho(ji,jj) = zsnowrho(jj)
             snowdz(ji,jj)  = zsnowdz(jj)
             ! keeps snow density below ice density
             IF ( snowrho(ji,jj) > xrhosmax ) THEN
                snowdz(ji,jj) = snowdz(ji,jj) * snowrho(ji,jj) / xrhosmax
                snowrho(ji,jj) = xrhosmax
             ENDIF         
          ENDDO
          snowmelt(ji)  = snowmelt(ji) + zflowliqt(nsnow) * ph2o
       ENDIF
       
    ENDDO
   
    CALL xios_orchidee_send_field("zrainfall",zrainfall/dt_sechiba)
    CALL xios_orchidee_send_field("zrefrz",zrefrz/dt_sechiba)

  END SUBROUTINE explicitsnow_refrz
  

!================================================================================================================================
!! SUBROUTINE   : explicitsnow_icemelt
!!
!>\BRIEF        Computes ice melt on ice sheet area (no refreezing process)
!!                
!! DESCRIPTION  : 
!! RECENT CHANGE(S) : None 
!!
!! MAIN OUTPUT VARIABLE(S): ice_sheet_melt, icetemp, icedz, grndflux
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_
!================================================================================================================================

  SUBROUTINE explicitsnow_icemelt(kjpindex,ice_sheet_mask,snowdz, precip_snow, zrainfall, snowmelt, vevapsno, icetemp,icedz,&
       grndflux, ice_sheet_melt)

    !! 0.1 Input variables
    INTEGER(i_std), INTENT (in)                             :: kjpindex        !! Domain size
    INTEGER(i_std),DIMENSION (kjpindex), INTENT(in)         :: ice_sheet_mask  !! Ice sheet mask
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(in)     :: snowdz          !! Snow depth
    REAL(r_std),DIMENSION(kjpindex),INTENT(in)              :: precip_snow     !! Snow rate (SWE) (kg/m2 per dt_sechiba)
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)           :: zrainfall       !! Rain precipitation on snow 
    !! @tex $(kg m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)            :: snowmelt        !! Snowmelt
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)           :: vevapsno        !! Snow evaporation  @tex ($kg m^{-2}$) @endtex 

    !! 0.2 Output variables
    REAL(r_std),DIMENSION(kjpindex),INTENT(out)             :: ice_sheet_melt  !! Ice melt [mm/dt_sechiba]

    !! 0.3 Modified variables
    REAL(r_std), DIMENSION (kjpindex,nice), INTENT(inout)   :: icetemp         !! Snow temperature 
    REAL(r_std), DIMENSION (kjpindex,nice), INTENT(inout)   :: icedz           !! Snow depth
    REAL(r_std),DIMENSION(kjpindex),INTENT(inout)           :: grndflux        !! Net energy input to soil

    !! 0.4 Local variables
    REAL(r_std),DIMENSION (kjpindex)                        :: meltxs          !! Residual snowmelt energy applied to underlying soil
    REAL(r_std),DIMENSION (nice)                            :: zphase_ice      !! Phase change (from ice to water) (J/m2)
    REAL(r_std),DIMENSION (nice)                            :: zicemelt        !! Ice melt (liquid water) (m)
    REAL(r_std),DIMENSION (nice)                            :: zicetemp
    REAL(r_std),DIMENSION (nice)                            :: zmeltxs         !! Excess melt
    INTEGER(i_std)                                          :: ji
    REAL(r_std),DIMENSION (kjpindex,nice)                   :: pcapa_ice       !! Heat capacity for ice
    REAL(r_std), DIMENSION (kjpindex)                       :: surf_massbal    !! surface mass balance  [mm/dt_sechiba]


    !initialize
    ice_sheet_melt(:) = 0.

    DO ji = 1, kjpindex

       IF (ice_sheet_mask(ji).GT.0) THEN ! ice grid point

          !! 1 Ice melting due to positive ice temperature

          !! 1.1 phase change (J/m2)

          pcapa_ice(ji,:) = rho_ice * (ZICETHRMHEAT1 + ZICETHRMHEAT2*(icetemp(ji,:)-tp_00)) ! formulation as in GRISLI prop_th_icetemp.f90

          zphase_ice(:)  = pcapa_ice(ji,:)*MAX(0.0, icetemp(ji,:)-tp_00) * icedz(ji,:)

          !! 1.2 ice melt     J/m2/(J/kg)=kg/m2 = mm

          zicemelt(:) = zphase_ice(:)/chalfu0

          !! 1.3 cool off the snow layer temperature due to melt

          zicetemp(:) = icetemp(ji,:) - zphase_ice(:)/(pcapa_ice(ji,:)* icedz(ji,:))

          icetemp(ji,:) = MIN(tp_00, zicetemp(:))

          zmeltxs(:)   = (zicetemp(:)-icetemp(ji,:))*pcapa_ice(ji,:)*icedz(ji,:)

          icedz(ji,:) = max(0.,icedz(ji,:) - zicemelt(:)/ rho_ice)

          ice_sheet_melt(ji)  = sum(zicemelt(:))

          ! excess heat from melting, using it to warm underlying ground to conserve energy
          meltxs(ji) = SUM(zmeltxs(:))/dt_sechiba  ! (W/m2)

          ! energy flux into the soil 
          grndflux(ji) = grndflux(ji) + meltxs(ji)  

       ELSE
          ice_sheet_melt(ji) = 0.
       ENDIF
    ENDDO
    surf_massbal(:) = precip_snow(:) + zrainfall(:) - snowmelt(:) - ice_sheet_melt(:) - vevapsno(:)

    CALL xios_orchidee_send_field("surf_massbal",surf_massbal/dt_sechiba)
    CALL xios_orchidee_send_field("ice_sheet_melt",ice_sheet_melt/dt_sechiba)

  END SUBROUTINE explicitsnow_icemelt
    


!================================================================================================================================
!! SUBROUTINE   : explicitsnow_icelevels
!!
!>\BRIEF        Update icelevels (constant) and icetemp
!!                
!! DESCRIPTION  : 
!! RECENT CHANGE(S) : None 
!!
!! MAIN OUTPUT VARIABLE(S): icetemp, icedz
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_
!================================================================================================================================

  SUBROUTINE explicitsnow_icelevels(kjpindex,ice_sheet_mask,ice_sheet_melt,icetemp,icedz)

    !! 0.1 Input variables

    INTEGER(i_std), INTENT (in)                             :: kjpindex        !! Domain size
    INTEGER(i_std),DIMENSION (kjpindex), INTENT(in)         :: ice_sheet_mask  !! Ice sheet mask
    REAL(r_std),DIMENSION(kjpindex),INTENT(in)              :: ice_sheet_melt  !! Ice melt [mm/dt_sechiba]

    !! 0.2 Output variables                                                                            

    !! 0.3 Modified variables
    REAL(r_std), DIMENSION (kjpindex,nice), INTENT(inout)  :: icetemp     !! Snow temperature 
    REAL(r_std), DIMENSION (kjpindex,nice), INTENT(inout)  :: icedz       !! Snow depth

    !! 0.4 Local variables           
    INTEGER(i_std)                                          :: ji,jj,xx,locxx
    REAL(r_std),DIMENSION(nice)                             :: znt_tmp              !! Ice grid layer boundary after melt
    REAL(r_std),DIMENSION (kjpindex,nice)                   :: icetemp_new          !! Update ice temperature
    REAL(r_std)                                             :: totdz                !! melt cumulated on all ice layers

    DO ji = 1, kjpindex
       ! icetemp is computed only on ice sheet point
       IF (ice_sheet_mask(ji).GT.0) THEN
          IF (ice_sheet_melt(ji).GT.min_sechiba) THEN
             totdz = 0.
             DO jj=1,nice-1
                icetemp_new(ji,jj) =  icetemp(ji,jj)*icedz(ji,jj)/ZSICOEF1(jj) + icetemp(ji,jj+1)*(ZSICOEF1(jj)-icedz(ji,jj))/ZSICOEF1(jj)
                totdz = totdz + ZSICOEF1(jj)-icedz(ji,jj) ! melt cumulated on all ice layers
             ENDDO
             icetemp_new(ji,nice) = icetemp(ji,nice)*(ZSICOEF1(nice) - totdz)/ZSICOEF1(nice) + tp_00 * totdz/ZSICOEF1(nice)
             icetemp(ji,:)=icetemp_new(ji,:)
             icedz(ji,:)=ZSICOEF1(:)
          ENDIF
       ENDIF
    ENDDO
  END SUBROUTINE explicitsnow_icelevels

!================================================================================================================================
!! SUBROUTINE   : explicitsnow_levels
!!
!>\BRIEF        Computes snow discretization based on given total snow depth
!!                
!! DESCRIPTION  : 
!! RECENT CHANGE(S) : compatible with 3 and 12 layers 
!!
!! MAIN OUTPUT VARIABLE(S): snowdz
!!
!! REFERENCE(S) :
!!
!! FLOWCHART    : None
!! \n
!_
!================================================================================================================================ 
  SUBROUTINE explicitsnow_levels(kjpindex,snow_thick, snowdz)

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                                     :: kjpindex                         !! Domain size
    REAL(r_std),DIMENSION (kjpindex),   INTENT (in)                :: snow_thick                       !! Total snow depth

    !! 0.2 Output variables

    !! 0.3 Modified variables
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(inout)         :: snowdz                           !! Snow depth [m]

    !! 0.4 Local variables
    INTEGER(i_std)                                                 :: ji, jj

    !! parameters for 3 layers snow model
    REAL(r_std), PARAMETER, DIMENSION(3)                           :: ZSGCOEF1  = (/0.25, 0.50, 0.25/) !! Snow grid parameters 
    REAL(r_std), PARAMETER, DIMENSION(2)                           :: ZSGCOEF2  = (/0.05, 0.34/)       !! Snow grid parameters 
    REAL(r_std), PARAMETER                                         :: ZSNOWTRANS = 0.20                !! Minimum total snow depth at which surface layer thickness is constant (m)
    REAL(r_std), PARAMETER                                         :: XSNOWCRITD = 0.03                !! (m)

    !! parameters for 12 layers snow model
    REAL(r_std), PARAMETER, DIMENSION(3)                           :: ZSGCOEF11  = (/0.3, 0.4, 0.3/)   !! Snow grid parameters : distribution of snow on layers 6, 7 and 8
    REAL(r_std), PARAMETER, DIMENSION(3)                           :: ZSGCOEF22  = (/0.3, 0.3, 0.3/)
    LOGICAL, PARAMETER, DIMENSION(12)                              :: mask_d_r = (/.true.,.true.,.true.,.true.,.true.,.false.,.false.,.false.,.true.,.true.,.true.,.true./)  ! True for nsnow=1,5 and nsnow=9,12 
    REAL(r_std), DIMENSION(kjpindex)                               :: d_r ! sum of the snow depth of layers 6, 7 and 8.


    IF (nsnow .eq. 3) THEN
       DO ji=1,kjpindex
          IF ( snow_thick(ji) .LE. (XSNOWCRITD+0.01)) THEN
             snowdz(ji,1) = MIN(0.01, snow_thick(ji)/nsnow)
             snowdz(ji,3) = MIN(0.01, snow_thick(ji)/nsnow)
             snowdz(ji,2) = snow_thick(ji) - snowdz(ji,1) - snowdz(ji,3)
          ENDIF
       ENDDO

       WHERE ( snow_thick(:) .LE. ZSNOWTRANS .AND. &
            snow_thick(:) .GT. (XSNOWCRITD+0.01) )
          snowdz(:,1) = snow_thick(:)*ZSGCOEF1(1)
          snowdz(:,2) = snow_thick(:)*ZSGCOEF1(2) 
          snowdz(:,3) = snow_thick(:)*ZSGCOEF1(3)
       END WHERE

       DO ji = 1,kjpindex
          IF (snow_thick(ji) .GT. ZSNOWTRANS) THEN
             snowdz(ji,1) = ZSGCOEF2(1)
             snowdz(ji,2) = (snow_thick(ji)-ZSGCOEF2(1))*ZSGCOEF2(2) + ZSGCOEF2(1)
             ! When using simple finite differences, limit the thickness
             ! factor between the top and 2nd layers to at most 10
             snowdz(ji,2)  = MIN(10*ZSGCOEF2(1),  snowdz(ji,2) )
             snowdz(ji,3)  = snow_thick(ji) - snowdz(ji,2) - snowdz(ji,1)
          ENDIF
       ENDDO

    ELSE IF (nsnow .eq. 12) THEN ! snow layering as in Decharme 2016 3.1.1      
       DO ji=1,kjpindex
          !! test if snowdz must be updated
          IF ((snowdz(ji,1).LT.(0.5*min(Zmax(1),snow_thick(ji)/12)) &
               .OR. (snowdz(ji,1).GT.(1.5*min(Zmax(1),snow_thick(ji)/12)))) &
               .OR. (snowdz(ji,2).LT.(0.5*min(Zmax(2),snow_thick(ji)/12))   &
               .OR. (snowdz(ji,2).GT.(1.5*min(Zmax(2),snow_thick(ji)/12)))) &
               .OR. (snowdz(ji,12).LT.(0.5*min(Zmax(12),snow_thick(ji)/12)) &
               .OR. (snowdz(ji,12).GT.(1.5*min(Zmax(12),snow_thick(ji)/12))))) THEN

             DO jj=1,5
                snowdz(ji,jj) = min(Zmax(jj),snow_thick(ji)/12)
             ENDDO
             DO jj=9,nsnow
                snowdz(ji,jj) = min(Zmax(jj),snow_thick(ji)/12)
             ENDDO
             d_r(ji) = snow_thick(ji)
             do jj=1,nsnow
                if (mask_d_r(jj)) then
                   d_r(ji) = d_r(ji) - snowdz(ji,jj)
                endif
             enddo
             snowdz(ji,6) = ZSGCOEF11(1)*d_r(ji) - min(0., ZSGCOEF22(1)*d_r(ji) - snowdz(ji,5))
             snowdz(ji,7) = ZSGCOEF11(2)*d_r(ji) + min(0., ZSGCOEF22(2)*d_r(ji) - snowdz(ji,5)) &
                  + min(0.,ZSGCOEF22(2)*d_r(ji) - snowdz(ji,9))
             snowdz(ji,8) = ZSGCOEF11(3)*d_r(ji) - min(0., ZSGCOEF22(3)*d_r(ji) - snowdz(ji,9))
          ENDIF
       ENDDO
    ENDIF

  END SUBROUTINE explicitsnow_levels


!!
!================================================================================================================================
!! SUBROUTINE   : explicitsnow_profile
!!
!>\BRIEF        
!!
!! DESCRIPTION  : In this routine solves the numerical snow thermal scheme, ie calculates the new snow temperature profile.
!!
!! RECENT CHANGE(S) : None
!! 
!! MAIN OUTPUT VARIABLE(S): snowtemp, temp_sol_add
!!
!! REFERENCE(S) :
!!
!! FLOWCHART    : None 
!! \n 
!_
!================================================================================================================================
  SUBROUTINE explicitsnow_profile (kjpindex, cgrnd_snow,dgrnd_snow,lambda_snow,temp_sol_new, snowtemp,snowdz,temp_sol_add)

    !! 0. Variables and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                               :: kjpindex     !! Domain size (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: temp_sol_new !! skin temperature
    REAL(r_std), DIMENSION (kjpindex,nsnow),INTENT (in)      :: cgrnd_snow   !! Integration coefficient for snow numerical scheme
    REAL(r_std), DIMENSION (kjpindex,nsnow),INTENT (in)      :: dgrnd_snow   !! Integration coefficient for snow numerical scheme
    REAL(r_std), DIMENSION (kjpindex),INTENT(in)             :: lambda_snow  !! Coefficient of the linear extrapolation of surface temperature 
    REAL(r_std), DIMENSION (kjpindex,nsnow),INTENT(in)       :: snowdz       !! Snow layer thickness [m]

    !! 0.2 Output variables

    !! 0.3 Modified variables
    REAL(r_std),DIMENSION (kjpindex,nsnow), INTENT (inout)   :: snowtemp
    REAL(r_std), DIMENSION (kjpindex),INTENT(inout)          :: temp_sol_add !! Additional energy to melt snow for snow ablation case (K)

    !! 0.4 Local variables
    INTEGER(i_std)                                           :: ji, jg
    !_
    !================================================================================================================================
    !! 1. Computes the snow temperatures
    DO ji = 1,kjpindex
       IF (SUM(snowdz(ji,:)) .GT. 0) THEN

          snowtemp(ji,1) = (lambda_snow(ji) * cgrnd_snow(ji,1) + (temp_sol_new(ji)+temp_sol_add(ji))) / &
               (lambda_snow(ji) * (un - dgrnd_snow(ji,1)) + un)
          temp_sol_add(ji) = zero 
          DO jg = 1,nsnow-1
             snowtemp(ji,jg+1) = cgrnd_snow(ji,jg) + dgrnd_snow(ji,jg) * snowtemp(ji,jg)
          ENDDO

       ENDIF
    ENDDO
    IF (printlev>=3) WRITE (numout,*) ' explicitsnow_profile done '

  END SUBROUTINE explicitsnow_profile
  
!!
!================================================================================================================================
!! SUBROUTINE   : explicitsnow_iceprofile
!!
!>\BRIEF        
!!
!! DESCRIPTION  : In this routine solves the numerical ice thermal scheme, ie calculates the new ice temperature profile.
!!
!! RECENT CHANGE(S) : None
!! 
!! MAIN OUTPUT VARIABLE(S): icetemp, temp_sol_add
!!
!! REFERENCE(S) :
!!
!! FLOWCHART    : None 
!! \n 
!_
!================================================================================================================================
  SUBROUTINE explicitsnow_iceprofile (kjpindex, cgrnd_ice,dgrnd_ice,lambda_ice,temp_sol_new,icedz,&
       ice_sheet_mask, snowdz, cgrnd_snow, dgrnd_snow, snowtemp, temp_sol_add, icetemp)

    !! 0. Variables and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                               :: kjpindex     !! Domain size (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)             :: temp_sol_new !! skin temperature
    REAL(r_std), DIMENSION (kjpindex,nice),INTENT(in)        :: cgrnd_ice    !! Integration coefficient for ice numerical scheme
    REAL(r_std), DIMENSION (kjpindex,nice),INTENT(in)        :: dgrnd_ice    !! Integration coefficient for ice numerical scheme
    REAL(r_std), DIMENSION (kjpindex),INTENT(in)             :: lambda_ice   !! Coefficient of the linear extrapolation of surface temperature 
    REAL(r_std), DIMENSION (kjpindex,nice),INTENT(in)        :: icedz        !! ice layer thickness
    INTEGER(i_std),DIMENSION (kjpindex), INTENT(in)          :: ice_sheet_mask !! Ice sheet mask
    REAL(r_std), DIMENSION (kjpindex,nsnow),INTENT(in)       :: snowdz       !! Snow layer thickness
    REAL(r_std), DIMENSION (kjpindex,nsnow),INTENT(in)       :: cgrnd_snow   !! Integration coefficient for snow numerical scheme
    REAL(r_std), DIMENSION (kjpindex,nsnow),INTENT(in)       :: dgrnd_snow   !! Integration coefficient for snow numerical scheme
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(in)      :: snowtemp     !! Snow temperature

    !! 0.2 Output variables

    !! 0.3 Modified variables
    REAL(r_std), DIMENSION (kjpindex),INTENT(inout)          :: temp_sol_add !! Additional energy to melt ice for ice ablation case (K)
    REAL(r_std), DIMENSION (kjpindex,nice), INTENT(inout)    :: icetemp      !! Ice temperature

    !! 0.4 Local variables
    INTEGER(i_std)                                           :: ji, jg
    !_
    !================================================================================================================================
    !! 1. Computes the ice temperatures
    DO ji = 1,kjpindex
       IF (ice_sheet_mask(ji).GT.0) THEN ! ice grid point
          IF (sum(snowdz(ji,:)).GT.0.) THEN ! grid point with snow
             icetemp(ji,1) = cgrnd_snow(ji,nsnow) + dgrnd_snow(ji,nsnow) * snowtemp(ji,nsnow)
             temp_sol_add(ji) = zero
             DO jg = 1,nice-1
                icetemp(ji,jg+1) = cgrnd_ice(ji,jg) + dgrnd_ice(ji,jg) * icetemp(ji,jg)
             ENDDO
          ELSE  ! no snow
             icetemp(ji,1) = (lambda_ice(ji) * cgrnd_ice(ji,1) + (temp_sol_new(ji)+temp_sol_add(ji))) / &
                  (lambda_ice(ji) * (un - dgrnd_ice(ji,1)) + un)
             temp_sol_add(ji) = zero 
             DO jg = 1,nice-1
                icetemp(ji,jg+1) = cgrnd_ice(ji,jg) + dgrnd_ice(ji,jg) * icetemp(ji,jg)
             ENDDO
          ENDIF
       ELSE
          icetemp(ji,:) = tp_00   ! ice temperature on non ice-sheet area (could be an other value, to be tested)         
       ENDIF
    ENDDO
    IF (printlev>=3) WRITE (numout,*) ' explicitsnow_iceprofile done '

  END SUBROUTINE explicitsnow_iceprofile
  
!! ================================================================================================================================
!! SUBROUTINE   : explicitsnow_read_reftempicefile
!!
!>\BRIEF          
!!
!! DESCRIPTION	: Read file with longterm ice temperature
!!                
!!
!! RECENT CHANGE(S) : None
!! 
!! MAIN OUTPUT VARIABLE(S): reftempice : Reference temperature for ice
!!                          
!! REFERENCE(S) :
!!
!! FLOWCHART    : None 
!! \n 
!_ ================================================================================================================================
  SUBROUTINE explicitsnow_read_reftempicefile(kjpindex,lalo,reftempice)

    USE interpweight

    IMPLICIT NONE

    !! 0. Variables and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in) :: kjpindex
    REAL(r_std), DIMENSION(kjpindex,2), INTENT(in) :: lalo

    !! 0.2 Output variables
    REAL(r_std), DIMENSION(kjpindex, nice), INTENT(out) :: reftempice

    !! 0.3 Local variables
    INTEGER(i_std) :: ib
    CHARACTER(LEN=80) :: filename
    REAL(r_std),DIMENSION(kjpindex) :: reftempice_file                       !! Horizontal temperature field interpolated from file [K]
    INTEGER(i_std),DIMENSION(kjpindex,8) :: neighbours
    REAL(r_std)                                          :: vmin, vmax       !! min/max values to use for the
    !!   renormalization
    REAL(r_std), DIMENSION(kjpindex)                     :: areftemp         !! Availability of data for the interpolation
    CHARACTER(LEN=80)                                    :: variablename     !! Variable to interpolate the file
    CHARACTER(LEN=80)                                    :: lonname, latname !! lon, lat names in input file
    REAL(r_std), DIMENSION(:), ALLOCATABLE               :: variabletypevals !! Values for all the types of the variable
    !!   (variabletypevals(1) = -un, not used)
    CHARACTER(LEN=50)                                    :: fractype         !! method of calculation of fraction
    !!   'XYKindTime': Input values are kinds 
    !!     of something with a temporal 
    !!     evolution on the dx*dy matrix'
    LOGICAL                                              :: nonegative       !! whether negative values should be removed
    CHARACTER(LEN=50)                                    :: maskingtype      !! Type of masking
    !!   'nomask': no-mask is applied
    !!   'mbelow': take values below maskvals(1)
    !!   'mabove': take values above maskvals(1)
    !!   'msumrange': take values within 2 ranges;
    !!      maskvals(2) <= SUM(vals(k)) <= maskvals(1)
    !!      maskvals(1) < SUM(vals(k)) <= maskvals(3)
    !!       (normalized by maskvals(3))
    !!   'var': mask values are taken from a 
    !!     variable inside the file (>0)
    REAL(r_std), DIMENSION(3)                            :: maskvals         !! values to use to mask (according to 
    !!   `maskingtype') 
    CHARACTER(LEN=250)                                   :: namemaskvar      !! name of the variable to use to mask 
    REAL(r_std)                                          :: reftemp_norefinf
    REAL(r_std)                                          :: reftemp_default  !! Default value


    !Config Key   = SOIL_REFTEMP_FILE
    !Config Desc  = File with climatological soil temperature
    !Config If    = READ_REFTEMPICE
    !Config Def   = reftempice.nc
    !Config Help  = 
    !Config Units = [FILE]
    filename = 'reftempice.nc'
    CALL getin_p('REFTEMPICE_FILE',filename)
    variablename = 'icetemp'

    IF (printlev >= 3) WRITE(numout,*) " in thermosoil_read_reftempicefile filename '" // TRIM(filename) // &
         "' variable name: '" //TRIM(variablename) // "'"

    ! For this case there are not types/categories. We have 'only' a continuos field
    ! Assigning values to vmin, vmax

    vmin = 0.
    vmax = 9999.

    !   For this file we do not need neightbours!
    neighbours = 0

    !! Variables for interpweight
    ! Type of calculation of cell fractions
    fractype = 'default'
    ! Name of the longitude and latitude in the input file
    lonname = 'nav_lon'
    latname = 'nav_lat'
    ! Default value when no value is get from input file
    reftemp_default = 1.
    ! Reference value when no value is get from input file
    reftemp_norefinf = 1.
    ! Should negative values be set to zero from input file?
    nonegative = .FALSE.
    ! Type of mask to apply to the input data (see header for more details)
    maskingtype = 'nomask'
    ! Values to use for the masking (here not used)
    maskvals = (/ undef_sechiba, undef_sechiba, undef_sechiba /)
    ! Name of the variable with the values for the mask in the input file (only if maskkingtype='var') (here not used)
    namemaskvar = ''

    CALL interpweight_2Dcont(kjpindex, 0, 0, lalo, resolution, neighbours,                            &
         contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype,        &
         maskvals, namemaskvar, -1, fractype, reftemp_default, reftemp_norefinf,                         &
         reftempice_file, areftemp)
    IF (printlev >= 5) WRITE(numout,*)'  thermosoil_read_reftempicefile after interpweight_2Dcont'

    ! Copy reftempice_file temperature to all ground levels
    DO ib=1, kjpindex
       reftempice(ib, :) = reftempice_file(ib)
    END DO

  END SUBROUTINE explicitsnow_read_reftempicefile


  !! ================================================================================================================================
  !! SUBROUTINE   : explicitsnow_read_ice_sheet_mask
  !!
  !>\BRIEF          
  !!
  !! DESCRIPTION	: Read file with ice sheet mask
  !!                
  !!
  !! RECENT CHANGE(S) : None
  !! 
  !! MAIN OUTPUT VARIABLE(S): ice_sheet_mask : ice sheet mask where there is ice between snow and soil
  !!                          
  !! REFERENCE(S) :
  !!
  !! FLOWCHART    : None 
  !! \n 
  !_ ================================================================================================================================
  SUBROUTINE explicitsnow_read_ice_sheet_mask(kjpindex,lalo,ice_sheet_mask)

    USE interpweight

    IMPLICIT NONE

    !! 0. Variables and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in) :: kjpindex
    REAL(r_std), DIMENSION(kjpindex,2), INTENT(in) :: lalo !! lon/lat

    !! 0.2 Output variables
    INTEGER(i_std), DIMENSION(kjpindex), INTENT(out) :: ice_sheet_mask

    !! 0.3 Local variables
    INTEGER(i_std) :: ji
    REAL(r_std),DIMENSION(kjpindex) :: ice_sheet_mask_real                   !! Ice sheet mask read in file 0 no ice, 1 ice sheet
    
    ! read ice sheet mask with xios (need for Antarctica)
    IF (xios_interpolation) THEN
       CALL xios_orchidee_recv_field('icemask',ice_sheet_mask_real)
       WRITE(numout,*)'explicitsnow_read_ice_sheet_mask after xios interpolation'
    ELSE
       call ipslerr_p(3,'explicitsnow_read_ice_sheet_mask', &
            'explicitsnow_read_ice_sheet_mask works only with XIOS_INTERPOLATION = yes','','')
    ENDIF

    ! Value must be an integer : a point could not be partialy ice covered
    DO ji=1, kjpindex
       IF (lalo(ji,1).LT.-60.) THEN ! Antarctic
          IF (ice_sheet_mask_real(ji).GT.0.) ice_sheet_mask(ji)=2 ! To avoid value between 0 and 2 on the coast by interpolation
       ELSE
          ice_sheet_mask(ji) = nint(ice_sheet_mask_real(ji))
       ENDIF
    END DO

    CALL xios_orchidee_send_field("interp_diag_ice_sheet_mask",real(ice_sheet_mask, r_std))

  END SUBROUTINE explicitsnow_read_ice_sheet_mask

  !! ================================================================================================================================
  !! SUBROUTINE   : explicitsnow_read_relief
  !!
  !>\BRIEF          
  !!
  !! DESCRIPTION	: Read topography file necessary for the initialization of density snowpack
  !!                
  !!
  !! RECENT CHANGE(S) : None
  !! 
  !! MAIN OUTPUT VARIABLE(S): topography from Relief.nc
  !!                          
  !! REFERENCE(S) :
  !!
  !! FLOWCHART    : None 
  !! \n 
  !_ ================================================================================================================================
  SUBROUTINE explicitsnow_read_relief(kjpindex, topography)

   IMPLICIT NONE

   !! 0. Variables and PARAMETER declaration

   !! 0.1 Input variables
   INTEGER(i_std), INTENT(in) :: kjpindex
   !! 0.2 Output variables
   REAL(r_std), DIMENSION(kjpindex), INTENT(out) :: topography              !! Topography field interpolated from file [m]

   !! 0.3 Local variables
   
   IF (xios_interpolation) THEN
      CALL xios_orchidee_recv_field('relief',topography)
      WRITE(numout,*)'  explicitsnow_read_relief after xios interpolation'
   END IF

 END SUBROUTINE explicitsnow_read_relief


!! ================================================================================================================================
!! SUBROUTINE   : explicitsnow_maxmass
!!
!>\BRIEF          
!!
!! DESCRIPTION	: limits the snow mass below a threshold (maxmass_snow)
!!                
!!
!! RECENT CHANGE(S) : adapted for nsnow=3 and nsnow=12
!! 
!! MAIN OUTPUT VARIABLE(S): snow, snowdz & snowmelt_from_maxmass
!!                          
!! REFERENCE(S) :
!!
!! FLOWCHART    : None 
!! \n 
!_ ================================================================================================================================  
  SUBROUTINE explicitsnow_maxmass(kjpindex,snowrho,soilcap,snow,snowdz,snowmelt_from_maxmass)

    IMPLICIT NONE
    ! variables pas necessaires a transmettre ?  :    maxmass_snow, chalfu0     

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                               :: kjpindex
    REAL(r_std),DIMENSION(kjpindex,nsnow),INTENT(in)         :: snowrho                !! Snow density (Kg/m^3)
    REAL(r_std),DIMENSION (kjpindex),INTENT(in)              :: soilcap                !! Soil heat capacity

    !! 0.2 Output variables
    REAL(r_std), DIMENSION (kjpindex), INTENT(out)           :: snow                   !! Snow mass [Kg/m^2]
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(inout)   :: snowdz                 !! Snow depth
    REAL(r_std), DIMENSION(kjpindex), INTENT(out)            :: snowmelt_from_maxmass  !! snowmelt to limit snow mass under maxmass_snow

    !! 0.3 Local variables
    INTEGER(i_std)                                           :: ji, jj
    INTEGER(i_std)                                           :: nsnowstart, nsnowstop  !! start and stop index layer 
    INTEGER(i_std)                                           :: locjj
    REAL(r_std)                                              :: snow_d1k
    REAL(r_std),DIMENSION(kjpindex,nsnow)                    :: WSNOWDZ !! snow mass (kg/m2)
    REAL(r_std),DIMENSION(kjpindex)                          :: SMASSC !! snow mass cumulated 
    REAL(r_std),DIMENSION(kjpindex,nsnow)                    :: SMASS  !! snow mass cumulated starting from bottom layer 
    REAL(r_std),DIMENSION(kjpindex)                          :: ZSNOWMELT_XS !! thickness of the layer that has been partially removed
    REAL(r_std),DIMENSION(kjpindex)                          :: ZSNOWDZ !! new thickness of the snow layer
    REAL(r_std),DIMENSION(kjpindex)                          :: snow_remove !! snow mass of layers that completely disapears


    ! initialisation
    !    snowmelt_from_maxmass(:) = zero

    IF (nsnow.EQ.3) THEN      ! snowmelt_from_maxmass can be applied to layers 3 to 2
       nsnowstart=3
       nsnowstop=2
    ELSEIF (nsnow.EQ.12) THEN ! snowmelt_from_maxmass can be applied to layers 9 to 6 because 10, 11 and 12 are thin
       nsnowstart=9
       nsnowstop=6
    ELSE
       CALL ipslerr_p(3,'explicitsnow_maxmass','Unsupported nsnow value', &
               'Currently supported nsnow values are 3 or 12','')
    ENDIF

    DO ji=1,kjpindex !domain size
       snow(ji) = SUM(snowrho(ji,:) * snowdz(ji,:))
       IF(snow(ji).GT.maxmass_snow)THEN
          IF (snow(ji).GT.(2.5*maxmass_snow)) THEN 
             ! melt much faster for points where accumulation is very very high
             snow_d1k = 10. * soilcap(ji) / chalfu0
          ELSE IF (snow(ji).GT.(1.5*maxmass_snow)) THEN 
             ! melt much faster for points where accumulation is very high
             snow_d1k = six * soilcap(ji) / chalfu0
          ELSE IF (snow(ji).GT.(1.2*maxmass_snow)) THEN 
             ! melt faster for points where accumulation is too high
             snow_d1k = trois * soilcap(ji) / chalfu0
          ELSE
             ! standard melting
             snow_d1k = un * soilcap(ji) / chalfu0
          ENDIF
          snowmelt_from_maxmass(ji) = MIN((snow(ji) - maxmass_snow),snow_d1k)
          ! Calculating the snow accumulation  
          WSNOWDZ(ji,:)= snowdz(ji,:)*snowrho(ji,:) ! WSNOWDZ in kg/m2
          SMASSC(ji)= 0.0
          DO jj=nsnowstart,nsnowstop,-1
             SMASS(ji,jj)   = SMASSC(ji) + WSNOWDZ(ji,jj)
             SMASSC(ji)     = SMASSC(ji) + WSNOWDZ(ji,jj)
          ENDDO

          ! Finding the layer
          locjj = nsnowstart  ! search the layer starting from nsnowstart
          DO jj=nsnowstart,nsnowstop,-1
             IF ((SMASS(ji,jj) .LE. snowmelt_from_maxmass(ji)) .AND. (SMASS(ji,jj-1) .GE. snowmelt_from_maxmass(ji)) ) THEN
                locjj=jj-1
             ENDIF
          ENDDO

          ! Calculating the removal of snow depth
          IF (locjj .EQ. nsnowstart) THEN
             ! ZSNOWMELT_XS : Epaisseur de la couche qui a disparu partiellement 
             ZSNOWMELT_XS(ji)  = snowmelt_from_maxmass(ji)/snowrho(ji,nsnowstart)
             ZSNOWDZ(ji)     = snowdz(ji,nsnowstart) - ZSNOWMELT_XS(ji)
             snowdz(ji,nsnowstart)     = MAX(0.0, ZSNOWDZ(ji))
          ELSE IF (locjj .LT. nsnowstart) THEN
             snow_remove(ji)=0   ! masse de neige des couches qui disparaissent totalement
             DO jj=nsnowstart,locjj+1,-1
                snow_remove(ji)=snow_remove(ji)+snowdz(ji,jj)*snowrho(ji,jj)
                snowdz(ji,jj)=0
             ENDDO
             ZSNOWMELT_XS(ji)  = (snowmelt_from_maxmass(ji) - snow_remove(ji))/snowrho(ji,locjj)
             ZSNOWDZ(ji)     = snowdz(ji,locjj) - ZSNOWMELT_XS(ji)
             snowdz(ji,locjj)     = MAX(0.0, ZSNOWDZ(ji))             
          ENDIF
       ENDIF
    ENDDO

  END SUBROUTINE explicitsnow_maxmass
  
!! ================================================================================================================================
!! SUBROUTINE   : explicitsnow_subli
!!
!>\BRIEF          
!!
!! DESCRIPTION	: sublimation on snow
!!                
!!
!! RECENT CHANGE(S) :
!! 
!! MAIN OUTPUT VARIABLE(S): snow, snowdz, subsnownobio, subsinksoil, snowliq, snowtemp, vevapsno
!!                          
!! REFERENCE(S) :
!!
!! FLOWCHART    : None 
!! \n 
!_ ================================================================================================================================  
  SUBROUTINE explicitsnow_subli(kjpindex,snowrho,frac_nobio,snowdz,vevapsno, &
       snowliq, snowtemp, snow, subsnownobio, subsinksoil)

    IMPLICIT NONE

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                               :: kjpindex
    REAL(r_std),DIMENSION(kjpindex,nsnow),INTENT(in)         :: snowrho                !! Snow density (Kg/m^3)
    REAL(r_std), DIMENSION (kjpindex,nnobio), INTENT(in)     :: frac_nobio

    !! 0.2 Output variables
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(inout)   :: snowdz                 !! Snow depth
    REAL(r_std), DIMENSION (kjpindex), INTENT(inout)         :: vevapsno               !! Snow evaporation  @tex ($kg m^{-2}$) @endtex
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(inout)   :: snowliq                !! Snow liquid content (m)
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(inout)   :: snowtemp               !! Snow temperature
    REAL(r_std), DIMENSION (kjpindex), INTENT(out)           :: snow                   !! Snow mass [Kg/m^2]
    REAL(r_std), DIMENSION (kjpindex,nnobio), INTENT(out)    :: subsnownobio           !! Sublimation of snow on other surface types (ice, lakes, ...)
    REAL(r_std), DIMENSION (kjpindex), INTENT(out)           :: subsinksoil            !! Excess of sublimation as a sink for the soil

    !! 0.3 Local variables
    INTEGER(i_std)                                           :: ji, jj
    INTEGER(i_std)                                           :: nsnowstart, nsnowstop  !! start and stop index layer 
    INTEGER(i_std)                                           :: locjj
    REAL(r_std),DIMENSION(kjpindex,nsnow)                    :: WSNOWDZ                !! snow mass (kg/m2)
    REAL(r_std),DIMENSION(kjpindex)                          :: SMASSC                 !! snow mass cumulated 
    REAL(r_std),DIMENSION(kjpindex,nsnow)                    :: SMASS                  !! snow mass cumulated starting from bottom layer 
    REAL(r_std),DIMENSION(kjpindex)                          :: ZSNOWEVAPS             !! thickness of the layer that has been removed by sublimation
    REAL(r_std),DIMENSION(kjpindex)                          :: ZSNOWDZ                !! new thickness of the snow layer
    REAL(r_std), DIMENSION (kjpindex)                        :: snowacc
    REAL(r_std), DIMENSION (kjpindex)                        :: subsnowveg

    snow(:) = 0.0
    DO ji=1,kjpindex !domain size
       snow(ji) = SUM(snowrho(ji,:) * snowdz(ji,:))
    ENDDO

    subsnownobio(:,:) = zero
    subsinksoil(:) = zero

    DO ji=1, kjpindex ! domain size
       
       ! Check that sublimation is possible.
       IF (vevapsno(ji) .GT. snow(ji)) THEN
          subsinksoil (ji) = vevapsno(ji) - snow(ji)
          ! Sublimation is thus limited to what is available
          vevapsno(ji) = snow(ji)
          snow(ji) = zero
          snowdz(ji,:)  =  0
          snowliq(ji,:)   =  0
          snowtemp(ji,:) = tp_00
          IF (frac_nobio(ji,iice) .GT. min_sechiba) THEN
             subsnownobio(ji,iice) = frac_nobio(ji,iice)*vevapsno(ji)
             subsnowveg(ji) = vevapsno(ji) - subsnownobio(ji,iice)
          ELSE
             subsnownobio(ji,iice) = zero
             subsnowveg(ji) = vevapsno(ji)
          ENDIF

       ELSE
          ! Calculating the snow accumulation  
          WSNOWDZ(ji,:)= snowdz(ji,:)*snowrho(ji,:)
          SMASSC (ji)= 0.0
          DO jj=1,nsnow
             SMASS(ji,jj)   = SMASSC(ji) + WSNOWDZ(ji,jj)
             SMASSC(ji)      = SMASSC(ji) + WSNOWDZ(ji,jj)
          ENDDO
          ! Finding the layer
          locjj=1
          DO jj=1,nsnow-1
             IF ((SMASS(ji,jj) .LE. vevapsno(ji)) .AND. (SMASS(ji,jj+1) .GE. vevapsno(ji)) ) THEN
                locjj=jj+1
             ENDIF
          ENDDO

          ! Calculating the removal of snow depth
          IF (locjj .EQ. 1) THEN
             ZSNOWEVAPS(ji)  = vevapsno(ji)/snowrho(ji,1)
             ZSNOWDZ(ji)     = snowdz(ji,1) - ZSNOWEVAPS(ji)
             snowdz(ji,1)     = MAX(0.0, ZSNOWDZ(ji))
          ELSE IF (locjj .GT. 1) THEN
             snowacc(ji)=0
             DO jj=1,locjj-1
                snowacc(ji)=snowacc(ji)+snowdz(ji,jj)*snowrho(ji,jj)
                snowdz(ji,jj)=0
             ENDDO
             ZSNOWEVAPS(ji)  = (vevapsno(ji)-snowacc(ji))/snowrho(ji,locjj)
             ZSNOWDZ(ji)     = snowdz(ji,locjj) - ZSNOWEVAPS(ji)
             snowdz(ji,locjj)     = MAX(0.0, ZSNOWDZ(ji))
          ENDIF

       ENDIF
    ENDDO

  END SUBROUTINE explicitsnow_subli
  
  
!! ================================================================================================================================
!! SUBROUTINE   : explicitsnow_age
!!
!>\BRIEF          
!!
!! DESCRIPTION	: compute snow age for albedo
!!                
!!
!! RECENT CHANGE(S) :
!! 
!! MAIN OUTPUT VARIABLE(S): snowage, snownobioage
!!                          
!! REFERENCE(S) :
!!
!! FLOWCHART    : None 
!! \n 
!_ ================================================================================================================================  
  SUBROUTINE explicitsnow_age(kjpindex,snow,precip_snow,precip_rain,frac_snow_nobio, &
       temp_sol_new,snow_age,snow_nobio_age)

    IMPLICIT NONE

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                               :: kjpindex
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: snow                   !! Snow mass [Kg/m^2]
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: precip_snow            !! Snowfall
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: precip_rain            !! Rainfall
    REAL(r_std), DIMENSION (kjpindex,nnobio), INTENT(in)     :: frac_snow_nobio        !! Snow cover fraction on non-vegeted area
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: temp_sol_new           !! Surface temperature

    !! 0.2 Output variables
    REAL(r_std), DIMENSION (kjpindex), INTENT(inout)         :: snow_age               !! Snow age
    REAL(r_std), DIMENSION (kjpindex,nnobio), INTENT(inout)  :: snow_nobio_age         !! Snow age on ice, lakes, ...


    !! 0.3 Local variables
    INTEGER(i_std)                                           :: ji
    REAL(r_std), DIMENSION (kjpindex)                        :: d_age                  !! Snow age change
    REAL(r_std), DIMENSION (kjpindex)                        :: xx                     !! Temporary


    DO ji = 1, kjpindex

       !! 5.1. Snow age on land 

       IF (snow(ji) .LE. zero) THEN
          snow_age(ji) = zero
       ELSE
          snow_age(ji) =(snow_age(ji) + (un - snow_age(ji)/max_snow_age) * dt_sechiba/one_day) &
               & * EXP(-precip_snow(ji) / snow_trans)
       ENDIF

       !! 5.2. Snow age on land ice (nobio)

       !! Age of snow on ice: a little bit different because in cold regions, we really
       !! cannot negect the effect of cold temperatures on snow metamorphism any more.

       IF ( (frac_snow_nobio(ji,iice) .LE. zero) .OR. (snow(ji) .LE. zero) ) THEN
          snow_nobio_age(ji,iice) = zero
       ELSE

          d_age(ji) = ( snow_nobio_age(ji,iice) + &
               &  (un - snow_nobio_age(ji,iice)/max_snow_age) * dt_sechiba/one_day ) * &
               &  EXP(-precip_snow(ji) / snow_trans_nobio) - snow_nobio_age(ji,iice)

          IF (d_age(ji) .GT. 0. ) THEN
             xx(ji) = MAX( tp_00 - temp_sol_new(ji), zero )
             xx(ji) = ( xx(ji) / omg1 ) ** omg2
             d_age(ji) = d_age(ji) / (un+xx(ji))

             ! age increase more rapidly if it rains: 
             IF (precip_rain(ji) .GT.0.) THEN 
                d_age(ji) = d_age(ji)*2.
             ENDIF
          ENDIF

          snow_nobio_age(ji,iice) = MAX( snow_nobio_age(ji,iice) + d_age(ji), zero )
       ENDIF
    ENDDO

  END SUBROUTINE explicitsnow_age
  
END MODULE explicitsnow
