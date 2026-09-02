! =================================================================================================================================
! MODULE       : thermosoil
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        Calculates the soil temperatures by solving the heat
!! diffusion equation within the soil. This module is only used with CWRR hydrology.
!!
!!\n DESCRIPTION : General important informations about the numerical scheme and
!!                 the soil vertical discretization:\n
!!               - the soil is zmaxt deep (by default 10m) and divided into "ngrnd" layers. 
!!                 From 0-zmaxh(default 2m), the discretization is the same as for hydrology. 
!!                 From zmaxh(2m) and below, the depth increase linearly (by default) or geometrically. \n
!!               - "jg" is usually used as the index going from 1 to ngrnd to describe the
!!                  layers, from top (jg=1) to bottom (jg=ngrnd)\n
!!               - the thermal numerical scheme is implicit finite differences.\n
!!                 -- When it is resolved in thermosoil_profile at the present timestep t, the
!!                 dependancy from the previous timestep (t-1) is hidden in the
!!                 integration coefficients cgrnd and dgrnd, which are therefore
!!                 calculated at the very end of thermosoil_main (call to
!!                 thermosoil_coef) for use in the next timestep.\n
!!                 -- At timestep t, the system becomes :\n 
!!
!!                              T(k+1)=cgrnd(k)+dgrnd(k)*T(k) \n
!!                                      -- EQ1 -- \n
!!
!!                 (the bottom boundary condition has been used to obtained this equation).\n
!!                 To solve it, the uppermost soil temperature T(1) is required.
!!                 It is obtained from the surface temperature Ts, which is
!!                 considered a linear extrapolation of T(1) and T(2)\n
!!
!!                           Ts=(1+lambda)*T(1) -lambda*T(2) \n 
!!                                      -- EQ2--\n
!!
!!                 -- caveat 1 : Ts is called 'temp_soil_new' in this routine,
!!                 don' t act.\n
!!                 -- caveat 2 : actually, the surface temperature at time t Ts
!!                 depends on the soil temperature at time t through the
!!                 ground heat flux. This is again implicitly solved, with Ts(t)
!!                 expressed as :\n
!!
!!                 soilcap*(Ts(t)-Ts(t-1))/dt=soilflx+otherfluxes(Ts(t))\n 
!!                                      -- EQ3 --\n
!!
!!                 and the dependency from the previous timestep is hidden in
!!                 soilcap and soilflx (apparent surface heat capacity and heat
!!                 flux respectively). Soilcap and soilflx are therefore
!!                 calculated at the previous timestep, at the very end of thermosoil
!!                 (final call to thermosoil_coef) and stored to be used at the next time step.
!!                 At timestep t, EQ3 is solved for Ts in enerbil, and Ts
!!                 is used in thermosoil to get T(1) and solve EQ1.\n
!!
!! - lambda is the @tex $\mu$ @endtex of F. Hourdin' s PhD thesis, equation (A28); ie the
!! coefficient of the linear extrapolation of Ts (surface temperature) from T1 and T2 (ptn(jg=1) and ptn(jg=2)), so that:\n
!! Ts= (1+lambda)*T(1)-lambda*T(2) --EQ2-- \n
!! lambda = (zlt(1))/((zlt(2)-zlt(1))) \n
!!
!! RECENT CHANGE(S) : - Change soil thermal properties to consider also soil texture, rev 2922.
!!                    - Change vertical discretization, rev 2917. Note: In the revised thermosoil, 
!!                    cstgrnd and lskin are not needed any more. The depth znt, zlt and dlt
!!                    are computed in vertical_soil and are in meter
!!
!! REFERENCE(S) : None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_sechiba/thermosoil.f90 $
!! $Date: 2026-04-03 18:22:57 +0200 (ven. 03 avril 2026) $
!! $Revision: 9463 $
!! \n
!_ ================================================================================================================================

MODULE thermosoil

  ! modules used :
  USE ioipsl
  USE ioipsl_para
  USE xios_orchidee
  USE constantes
  USE time, ONLY : one_day, dt_sechiba
  USE constantes_soil
  USE sechiba_io_p
  USE grid
  USE pft_parameters
  USE interpol_help
  USE vertical_soil

  IMPLICIT NONE

  !private and public routines :
  PRIVATE
  PUBLIC :: thermosoil_main, thermosoil_clear,  thermosoil_initialize, thermosoil_finalize, thermosoil_xios_initialize

  REAL(r_std), SAVE                               :: lambda                   !! See Module description
!$OMP THREADPRIVATE(lambda)
  REAL(r_std), SAVE                               :: fz1, zalph               !! usefull constants for diverse use
!$OMP THREADPRIVATE(fz1, zalph)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:) :: ptn                    !! Vertically discretized soil temperature, per soil layer and per pft (K)
!$OMP THREADPRIVATE(ptn)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)   :: ptn_pftmean            !! Vertically discretized soil temperature, mean across all pfts
!$OMP THREADPRIVATE(ptn_pftmean)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: dz1                      !! numerical constant used in the thermal numerical
                                                                              !! scheme  @tex ($m^{-1}$) @endtex. ; it corresponds
                                                                              !! to the coefficient  @tex $d_k$ @endtex of equation
                                                                              !! (A.12) in F. Hourdin PhD thesis.
!$OMP THREADPRIVATE(dz1)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: cgrnd                    !! integration coefficient for the numerical scheme,
                                                                              !! see eq.1
!$OMP THREADPRIVATE(cgrnd)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: dgrnd                    !! integration coefficient for the numerical scheme,
                                                                              !! see eq.1
!$OMP THREADPRIVATE(dgrnd)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: pcapa                    !! volumetric vertically discretized soil heat 
                                                                              !! capacity  @tex ($J K^{-1} m^{-3}$) @endtex. 
                                                                              !! It depends on the soil
                                                                              !! moisture content (shum_ngrnd_perma) and is calculated at 
                                                                              !! each time step in thermosoil_coef.
!$OMP THREADPRIVATE(pcapa)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:) :: pcapa_per_pft          !! volumetric vertically discretized soil heat per pft 
                                                                              !! capacity  @tex ($J K^{-1} m^{-3}$) @endtex. 
                                                                              !! It depends on the soil
                                                                              !! moisture content (shum_ngrnd_perma) and is calculated at 
                                                                              !! each time step in thermosoil_coef.
!$OMP THREADPRIVATE(pcapa_per_pft)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: pkappa                   !! vertically discretized soil thermal conductivity 
                                                                              !!  @tex ($W K^{-1} m^{-1}$) @endtex. Same as pcapa.
!$OMP THREADPRIVATE(pkappa)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:) :: pkappa_per_pft         !! vertically discretized soil thermal conductivity per pft 
                                                                              !!  @tex ($W K^{-1} m^{-1}$) @endtex. Same as pcapa.
!$OMP THREADPRIVATE(pkappa_per_pft)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: pcapa_snow               !! volumetric vertically discretized snow heat 
                                                                              !! capacity @tex ($J K^{-1} m^{-3}$) @endtex. 
!$OMP THREADPRIVATE(pcapa_snow)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: pkappa_snow              !! vertically discretized snow thermal conductivity 
                                                                              !! @tex ($W K^{-1} m^{-1}$) @endtex.
!$OMP THREADPRIVATE(pkappa_snow)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: pcapa_ice                !! volumetric vertically discretized ice heat 
                                                                              !! capacity @tex ($J K^{-1} m^{-3}$) @endtex. 
!$OMP THREADPRIVATE(pcapa_ice)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: pkappa_ice              !! vertically discretized ice thermal conductivity 
                                                                              !! @tex ($W K^{-1} m^{-1}$) @endtex.
!$OMP THREADPRIVATE(pkappa_ice)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: pcapa_en                 !! heat capacity used for surfheat_incr and 
                                                                              !! coldcont_incr 
!$OMP THREADPRIVATE(pcapa_en)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:) :: pcapa_en_per_pft       !! heat capacity used for surfheat_incr and per pft 
                                                                              !! coldcont_incr 
!$OMP THREADPRIVATE(pcapa_en_per_pft)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:) :: ptn_beg                !! ptn as it is after thermosoil_profile but before thermosoil_coef, 
                                                                              !! used in thermosoil_readjust
!$OMP THREADPRIVATE(ptn_beg)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: temp_sol_beg             !! Surface temperature at previous timestep (K) 
!$OMP THREADPRIVATE(temp_sol_beg)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: surfheat_incr            !! Change in soil heat content during the timestep 
                                                                              !!  @tex ($J$) @endtex.
!$OMP THREADPRIVATE(surfheat_incr)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: coldcont_incr            !! Change in snow heat content  @tex ($J$) @endtex.
!$OMP THREADPRIVATE(coldcont_incr)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: shum_ngrnd_perma         !! Water saturation degree on the soil layers define by the thermic (0-1, dimensionless)
!$OMP THREADPRIVATE(shum_ngrnd_perma)

  !  Variables related to soil freezing
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: profil_froz              !! Frozen fraction of the soil on hydrological levels (-)
!$OMP THREADPRIVATE(profil_froz)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:):: profil_froz_pft          !! Frozen fraction of the soil on hydrological levels per pft (-)
!$OMP THREADPRIVATE(profil_froz_pft) 
  REAL(r_std), ALLOCATABLE, SAVE,DIMENSION(:,:)   :: refsoc                   !! Initialize soil organic carbon only used to calculate thermal insulating effect [gC/m3]
!$OMP THREADPRIVATE(refsoc)
  REAL(r_std), ALLOCATABLE, SAVE,DIMENSION(:,:)   :: fSOC       !! Initialize soil organic carbon only used to calculate thermal insulating effect [kg/kg]
!$OMP THREADPRIVATE(fSOC)

  REAL(r_std), ALLOCATABLE, SAVE,DIMENSION(:,:)   :: BD_FM      !! Initialize soil bulk density in fine matter only used to calculate thermal insulating effect [kg/m3]
!$OMP THREADPRIVATE(BD_FM)

  REAL(r_std), ALLOCATABLE, SAVE,DIMENSION(:,:)   :: Phi_CM     !! Initialize volume fraction of coarse matter only used to calculate thermal insulating effect [m3/m3]
!$OMP THREADPRIVATE(Phi_CM)
  REAL(r_std),ALLOCATABLE, SAVE, DIMENSION (:)    :: e_soil_lat               !! Latent heat released or consumed in the freezing/thawing processes summed vertically 
                                                                              !! for the whole soil (J/m2) and on the whole simulation to check/correct energy conservation
!$OMP THREADPRIVATE(e_soil_lat)
  REAL(r_std),ALLOCATABLE, SAVE, DIMENSION (:,:)  :: pcappa_supp              !! Additional surfacic heat capacity due to soil freezing for each soil layer (J/K/m2)
!$OMP THREADPRIVATE(pcappa_supp)    
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: dz5                      !! Used for numerical calculation [-]
!$OMP THREADPRIVATE(dz5)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: mcs                      !! Saturation humidity [m3/m3]
!$OMP THREADPRIVATE(mcs)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: QZ                       !! quartz content [-]
!$OMP THREADPRIVATE(QZ)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: so_capa_dry              !! Dry soil Heat capacity of soils,J.m^{-3}.K^{-1} 
!$OMP THREADPRIVATE(so_capa_dry)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: so_capa_ice              !! Heat capacity of saturated frozen soil (J/K/m3)
!$OMP THREADPRIVATE(so_capa_ice)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: mc_layt                  !! Volumetric soil moisture (liquid+ice) (m3/m3) on the thermodynamical levels at interface
!$OMP THREADPRIVATE(mc_layt)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: mcl_layt                 !! Volumetric soil moisture (liquid) (m3/m3) on the thermodynamical levels at interface
!$OMP THREADPRIVATE(mcl_layt)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: tmc_layt                 !! Total soil moisture content for each layer (liquid+ice) (mm) on the thermodynamical levels
!$OMP THREADPRIVATE(tmc_layt)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:) :: mc_layt_pft            !! Volumetric soil moisture (liquid+ice) (m3/m3) on the thermodynamical levels at interface per pft
!$OMP THREADPRIVATE(mc_layt_pft)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:) :: mcl_layt_pft           !! Volumetric soil moisture (liquid) (m3/m3) on the thermodynamical levels at interface per pft
!$OMP THREADPRIVATE(mcl_layt_pft)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:) :: tmc_layt_pft           !! Total soil moisture content for each layer (liquid+ice) (mm) on the thermodynamical levels per pft
!$OMP THREADPRIVATE(tmc_layt_pft)
  INTEGER(i_std), SAVE                              :: brk_flag = 0           !! Flag to consider bedrock: 0.no; 1.yes
!$OMP THREADPRIVATE(brk_flag)

!Vertical Permafrost Carbon
  LOGICAL, SAVE                                      :: use_refsoc                   !! which SOC to use in thermix: a map that we read or the SOC calculated by the model.
!$OMP THREADPRIVATE(use_refsoc)
  LOGICAL, SAVE                                      :: ok_moss                      !! Add moss in thermal conductivity and capacity calculation
!$OMP THREADPRIVATE(ok_moss)
 INTEGER(i_std), SAVE                                :: use_soilc_method             !! Method to average thermal conductivity of mineral soil and organic soil:
                                                                                     !! Possible values are SOILC_METHOD_ARITHMETIC or SOILC_METHOD_GEOMETRIC
!$OMP THREADPRIVATE(use_soilc_method)
  INTEGER(i_std), PARAMETER                          :: SOILC_METHOD_ARITHMETIC = 1  !! Index to use arithmetic mean
  INTEGER(i_std), PARAMETER                          :: SOILC_METHOD_GEOMETRIC  = 2  !! Index to use geometric mean 
                                                                                
  INTEGER(i_std), SAVE                               :: snow_cond_method                !! Method to calculate snow thermal conductivity  
                                                                                        !! Possible values are SNOW_COND_METHOD_DEFAULT and SNOW_COND_METHOD_DECHARME16 and SNOW_COND_METHOD_SNOWICE
!$OMP THREADPRIVATE(snow_cond_method)
  INTEGER(i_std), PARAMETER                          :: SNOW_COND_METHOD_DEFAULT = 1    !! Index for original method
  INTEGER(i_std), PARAMETER                          :: SNOW_COND_METHOD_DECHARME16 = 2 !! Index to follow the method by Decharme et al 2016
  INTEGER(i_std), PARAMETER                          :: SNOW_COND_METHOD_SNOWICE = 3    !! Index to follow the method by Charbit, Dumas et al 2022

  INTEGER(i_std), SAVE                               :: use_moss_cond_method         !! Method to average thermal conductivity of moss and soil:
                                                                                     !! Possible values are moss_combine or moss_replace
!$OMP THREADPRIVATE(use_moss_cond_method)
  INTEGER(i_std), PARAMETER                          :: moss_combine = 1             !! Index to combine soil and moss
  INTEGER(i_std), PARAMETER                          :: moss_replace = 2             !! Index to replace soil by moss

CONTAINS


!! =============================================================================================================================
!! SUBROUTINE:    thermosoil_xios_initialize
!!
!>\BRIEF	  Initialize xios dependant defintion before closing context defintion
!!
!! DESCRIPTION:	  Initialize xios dependant defintion before closing context defintion
!!                Reading is deactivated if the sechiba restart file exists because the
!!                variable should be in the restart file already. 
!!
!! \n
!_ ==============================================================================================================================

  SUBROUTINE thermosoil_xios_initialize

    CHARACTER(LEN=255) :: filename, name


    ! Initialization related to the reading of reftemp file
    filename = 'reftemp.nc'
    CALL getin_p('REFTEMP_FILE',filename)
    
    name = filename(1:LEN_TRIM(FILENAME)-3)
    CALL xios_orchidee_set_file_attr("reftemp_file",name=name)

    ! Check if the reftemp file will be read by XIOS, by IOIPSL or not at all    
    IF (xios_interpolation .AND. read_reftemp .AND. restname_in=='NONE') THEN
       ! The reftemp file will be read using XIOS
       IF (printlev>=2) WRITE(numout,*) 'Reading of reftemp file will be done later using XIOS. The filename is ', filename
    ELSE
       IF (.NOT. read_reftemp) THEN
          IF (printlev>=2) WRITE (numout,*) 'No reading of reftemp will be done because read_reftemp=FALSE'
       ELSE IF (restname_in=='NONE') THEN
          IF (printlev>=2) WRITE (numout,*) 'The reftemp file will be read later by IOIPSL'
       ELSE
          IF (printlev>=2) WRITE (numout,*) 'The reftemp file will not be read because the restart file exists.'
       END IF

       ! The reftemp file will not be read by XIOS. Now deactivate albedo for XIOS.
       CALL xios_orchidee_set_file_attr("reftemp_file",enabled=.FALSE.)
       CALL xios_orchidee_set_field_attr("reftemp_interp",enabled=.FALSE.)
    ENDIF


    !! Initialization related to the reading of refsoc file

    ! Get to know if the refsoc file should be read. 
    use_soilc_insulation = .FALSE.
    CALL getin_p('USE_SOILC_INSULATION', use_soilc_insulation)
    IF (ok_freeze_thermix .AND. use_soilc_insulation) THEN
       use_refsoc = .TRUE.
       CALL getin_p('USE_REFSOC',use_refsoc)
    ELSE
       use_refsoc=.FALSE.
    END IF

    ! Check if the refsoc file will be read by XIOS, by IOIPSL or not at all    
    IF (xios_interpolation .AND. use_refsoc .AND. restname_in=='NONE') THEN
       ! The refsoc file will be read using XIOS
       ! Read the file name for run.def and set it to XIOS
       filename = 'refSOC.nc'
       CALL getin_p('SOIL_REFSOC_FILE',filename)
       name = filename(1:LEN_TRIM(FILENAME)-3)
       CALL xios_orchidee_set_file_attr("refsoc_file",name=name)
       
       IF (printlev>=2) WRITE(numout,*) 'Reading of refsoc file will be done later using XIOS. The filename is ', filename
    ELSE
       IF (.NOT. use_refsoc) THEN
          IF (printlev>=2) WRITE (numout,*) 'No reading of refsoc will be done because use_refsoc=FALSE'
       ELSE IF (restname_in=='NONE') THEN
          IF (printlev>=2) WRITE (numout,*) 'The refsoc file will be read later by IOIPSL'
       ELSE
          IF (printlev>=2) WRITE (numout,*) 'The refsoc file will not be read because the restart file exists.'
       END IF

       ! The reftemp file will not be read by XIOS. Now deactivate albedo for XIOS.
       CALL xios_orchidee_set_file_attr("refsoc_file",enabled=.FALSE.)
       CALL xios_orchidee_set_field_attr("mask3D_refsoc_file",enabled=.FALSE.)
       CALL xios_orchidee_set_field_attr("refsoc_interptmp",enabled=.FALSE.)
       CALL xios_orchidee_set_field_attr("refsoc_interp",enabled=.FALSE.)
    ENDIF

    

  END SUBROUTINE thermosoil_xios_initialize

  !!  =============================================================================================================================
  !! SUBROUTINE		 		    : thermosoil_initialize
  !!
  !>\BRIEF			            Allocate module variables, read from restart file or initialize with default values
  !!
  !! DESCRIPTION			    : Allocate module variables, read from restart file or initialize with default values.
  !!                                          Call thermosoil_var_init to calculate physical constants. 
  !!                                          Call thermosoil_coef to calculate thermal soil properties.
  !!
  !! RECENT CHANGE(S)			    : None
  !!
  !! REFERENCE(S)			    : None
  !! 
  !! FLOWCHART                              : None
  !! \n
  !_ ==============================================================================================================================
  SUBROUTINE thermosoil_initialize (kjit,          kjpindex,   rest_id,                   &
                                    temp_sol_new,  snow,       shumdiag_perma,            &
                                    soilcap,       soilflx,    depth_organic_soil,        &
                                    stempdiag,     ftempdiag,  gtemp,       mc_layh,      &
                                    mcl_layh,      tmc_layh,   njsc,  soiltile,           &
                                    frac_snow_veg,frac_snow_nobio,totfrac_nobio,          &
                                    snowdz, snowrho, snowtemp, lambda_snow, cgrnd_snow,   &
                                    dgrnd_snow,    lambda_ice, cgrnd_ice,   dgrnd_ice,    &
                                    icetemp,       icedz,      ice_sheet_mask,   pb,      &
                                    som_total, veget_max, mc_layh_pft, mcl_layh_pft, tmc_layh_pft, &
                                    ptnlev1)

    !! 0. Variable and parameter declaration
    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                            :: kjit             !! Time step number (unitless) 
    INTEGER(i_std), INTENT(in)                            :: kjpindex         !! Domain size (unitless)
    INTEGER(i_std),INTENT (in)                            :: rest_id          !! Restart file identifier (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)         :: temp_sol_new     !! Surface temperature at the present time-step,
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)         :: snow             !! Snow mass (kg)
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (in)    :: shumdiag_perma   !! Soil saturation degree (0-1, unitless)
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (in)    :: mc_layh          !! Volumetric soil moisture content (liquid+ice) for hydrological layers, at node (m3/m3)
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (in)    :: mcl_layh         !! Volumetric soil moisture content (liquid) for hydrological layers, at node (m3/m3)
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (in)    :: tmc_layh         !! Total soil moisture content(liquid+ice) for hydrological layers (mm)
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)      :: njsc             !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless)
    REAL(r_std),DIMENSION (kjpindex,nstm), INTENT (in)    :: soiltile         !! Fraction of each soil tile within vegtot (0-1, unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)          :: frac_snow_veg    !! Snow cover fraction on vegeted area
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT(in)   :: frac_snow_nobio  !! Snow cover fraction on non-vegeted area
    REAL(r_std),DIMENSION (kjpindex),INTENT(in)           :: totfrac_nobio    !! Total fraction of continental ice+lakes+cities+...
                                                                              !! (unitless,0-1)
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT (in)    :: snowdz           !! Snow depth [m]
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(in)   :: snowrho          !! Snow density (Kg/m^3)
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(in)   :: snowtemp         !! Snow temperature (K)
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)        :: pb               !! Surface presure (hPa)
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm,nelements), INTENT (in) :: som_total  !! total soil organic matter for use in thermal calcs (g/m**3)
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(in)     :: veget_max         !! Fraction of vegetation type
    REAL(r_std),DIMENSION (kjpindex,nslm,nvm), INTENT (in):: mc_layh_pft       !! Volumetric soil moisture content for each layer in hydrol at nodes(liquid+ice) [m/s]
    REAL(r_std),DIMENSION (kjpindex,nslm,nvm), INTENT (in):: mcl_layh_pft      !! Volumetric soil moisture content for each layer in hydrol at nodes(liquid) [m/s]
    REAL(r_std),DIMENSION (kjpindex,nslm,nvm), INTENT (in):: tmc_layh_pft      !! Total soil moisture content for each layer in hydrol(liquid+ice) [mm]

    !! 0.2 Output variables
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)        :: soilcap          !! apparent surface heat capacity considering snow and soil surface (J m-2 K-1)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)        :: soilflx          !! apparent soil heat flux considering snow and soil surface (W m-2)
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (out)   :: stempdiag        !! temperature profile on the levels in hydrol(K)
    REAL(r_std),DIMENSION (kjpindex,ngrnd), INTENT (out)  :: ftempdiag        !! temperature profile on full depth for stream temperature (K)
    REAL(r_std),DIMENSION (kjpindex),INTENT(out)          :: gtemp            !! First soil layer temperature
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)        :: ptnlev1          !! 1st level soil temperature   

    !! 0.3 Modified variables
    REAL(r_std), DIMENSION( kjpindex), INTENT (inout)      :: depth_organic_soil!! Depth at which there is still organic matter (m) 
    REAL(r_std), DIMENSION (kjpindex), INTENT(inout)       :: lambda_snow     !! Coefficient of the linear extrapolation of surface temperature 
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT (inout):: cgrnd_snow      !! Integration coefficient for snow numerical scheme
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT (inout):: dgrnd_snow      !! Integration coefficient for snow numerical scheme
    REAL(r_std), DIMENSION (kjpindex), INTENT(inout)       :: lambda_ice      !! Coefficient of the linear extrapolation of ice surface temperature
    REAL(r_std), DIMENSION (kjpindex,nice), INTENT (out)   :: cgrnd_ice       !! Integration coefficient for ice numerical scheme
    REAL(r_std), DIMENSION (kjpindex,nice), INTENT (out)   :: dgrnd_ice       !! Integration coefficient for ice numerical scheme
    REAL(r_std), DIMENSION (kjpindex,nice), INTENT(in)     :: icetemp         !! Ice temperature (K)
    REAL(r_std), DIMENSION (kjpindex,nice), INTENT (in)    :: icedz           !! Ice depth
    INTEGER(i_std), DIMENSION (kjpindex), INTENT (in)      :: ice_sheet_mask  !! Ice sheet mask
    
        

    !! 0.4 Local variables
    INTEGER(i_std)                                        :: ier, i, jg, iv, m, jv, jsc
    LOGICAL                                               :: calculate_coef   !! Local flag to initialize variables by call to thermosoil_coef   
    CHARACTER(LEN=10)                                     :: part_str         !! string suffix indicating an index
    CHARACTER(LEN=80)                                     :: var_name
    
     
    !! Variables needed only if we impose SOC and bulk density for site simulations with flag IMPOSE_SOC_BD)
    REAL(r_std), DIMENSION(ngrnd) :: fSOC_init !! Initial mass fraction of SOC in the fine matter (discretized on 18 levels of thermosoil) [kg/kg]
    REAL(r_std), DIMENSION(ngrnd) :: BDOD_init !! Initial bulk density of fine matter [kg/m3]
    REAL(r_std), DIMENSION(ngrnd) :: CFVO_init !! Initial coarse fragments volume fraction [m3/m3]

!_ ================================================================================================================================
    !
    !  !! Flag to consider bedrock at deeper layers
    !  !! It affects heat capacity and thermal conductivity (energy balance). 
    !
    !Config Key  = BEDROCK_FLAG
    !Config Desc = Flag to consider bedrock at deeper layers.
    !Config If   = 
    !Config Def  = 0
    !Config Help = 0, no, 1, yes. 
    !Config Units = [FLAG]
    brk_flag = 0
    CALL getin_p('BEDROCK_FLAG', brk_flag)

    use_refsoc = .FALSE.
    ok_moss    = .FALSE.
    IF (ok_freeze_thermix) THEN
        IF (use_soilc_insulation) THEN
           !Config Key  = USE_REFSOC
           !Config Desc = Read a SOC map to perform the insolation effect
           !Config If   = USE_SOILC_INSULATION 
           !Config Def  = TRUE
           !Config Help = 
           !Config Units = [FLAG]
           use_refsoc = .TRUE.
           CALL getin_p('USE_REFSOC',use_refsoc)
           !
           !Config Key  = OK_MOSS
           !Config Desc = Activate moss thermal effect
           !Config If   = USE_SOILC_INSULATION
           !Config Def  = FALSE
           !Config Help = Activate moss thermal effect (thermal capacity and conductivity).
           !Config Units= [FLAG]
           ok_moss = .TRUE.
           CALL getin_p('OK_MOSS',ok_moss)
         ENDIF

         IF (ok_moss .AND. .NOT. use_soilc_insulation) THEN
            CALL ipslerr_p(3, "thermosoil_initialize", &
                 "OK_MOSS can not be activated without USE_SOILC_INSULATION=y",&
                 "Please, change parameter value in run.def","")
         END IF

         impose_soc_bd=.FALSE.
         IF (use_soilc_insulation) THEN
            IF (use_thermo_newdisc) THEN
               !Config Key = IMPOSE_SOC_BD
               !Config Desc = Flag to manually impose SOC and bulk density for site simulations (values defined in run.def)
               !Config If = USE_THERMO_NEWDISC, USE_SOILC_INSULATION
               !Config Def  = FALSE
               !Config Help =
               !Config Units = [FLAG]
               CALL getin_p('IMPOSE_SOC_BD', impose_soc_bd)
               
               IF (impose_soc_bd .AND. .NOT. use_refsoc) THEN
                  !Config Key = FSOC_INIT
                  !Config Desc = Initial values of fSOC discretized on the 18 thermosoil layers
                  !Config If = USE_THERMO_NEWDISC, IMPOSE_SOC_BD, USE_REFSOC == False
                  !Config Def  = 0.0
                  !Config Help = 
                  !Config Units = [kg SOC/kg fine matter]
                  fSOC_init(:) = 0.0
                  CALL getin_p('FSOC_INIT', fSOC_init)
                  
                  !Config Key = BDOD_INIT
                  !Config Desc = Initial values of BDOD discretized on the 18 thermosoil layers
                  !Config If = USE_THERMO_NEWDISC, IMPOSE_SOC_BD, USE_REFSOC == False
                  !Config Def  = 1620.0
                  !Config Help = 
                  !Config Units = [kg/m3]
                  BDOD_init(:) = 1620.
                  CALL getin_p('BDOD_INIT', BDOD_init)
                  
                  !Config Key = CFVO_INIT
                  !Config Desc = Initial values of CFVO discretized on the 18 thermosoil layers
                  !Config If = USE_THERMO_NEWDISC, IMPOSE_SOC_BD, USE_REFSOC == False
                  !Config Def  = 0.0
                  !Config Help = 
                  !Config Units = [m3/m3]
                  CFVO_init(:) = 0.0
                  CALL getin_p('CFVO_INIT', CFVO_init)
               ELSE IF (impose_soc_bd .AND. use_refsoc) THEN
                  CALL ipslerr_p(3,'thermosoil_initialize', 'IMPOSE_SOC_BD=TRUE and USE_REFSOC=TRUE', &
                       'This combination is not possible. Only one of IMPOSE_SOC_BD or USE_REFSOC can be true','Change in run.def')                 
               ENDIF
            ENDIF
         ENDIF !(use_soilc_insulation)
      ENDIF !(ok_freeze_thermix .AND. ok_soil_carbon_discretization)
   
    !Config Key  = USE_SOILC_METHOD
    !Config Desc = Flag to control the way to average thermal conductivity of mineral soil and organic soil
    !Config If   = OK_SOIL_CARBON_DISCRETIZATION
    !Config Def  = 1
    !Config Help = 1=arithmetic mean ; 2=geometric mean 
    !Config Units = [FLAG]
    use_soilc_method = SOILC_METHOD_ARITHMETIC
    CALL getin_p('USE_SOILC_METHOD', use_soilc_method)
    IF ( (use_soilc_method /= SOILC_METHOD_ARITHMETIC) .AND. (use_soilc_method /= SOILC_METHOD_GEOMETRIC) ) THEN
       CALL ipslerr_p(3,'thermosoil_initialize', 'Error in variable use_soilc_method','USE_SOILC_METHOD must be equal 1 or 2','')
    END IF


    !Config Key  = SNOW_COND_METHOD
    !Config Desc = Flag to choose the way to calculate snow thermal conductivity 
    !Config If   = OK_SOIL_CARBON_DISCRETIZATION
    !Config Def  = 1
    !Config Help =  1: original 2: follows Decharme et al 2016
    !Config Units = [1=original method, 2=method by Decharme et al 2016]
    snow_cond_method = SNOW_COND_METHOD_DEFAULT
    CALL getin_p('SNOW_COND_METHOD', snow_cond_method)


    IF (ok_moss) THEN
       !Config Key  = USE_MOSS_COND_METHOD
       !Config Desc = Flag to control the way to average thermal conductivity of moss and soil: 1=moss_combine, 2=moss_replace
       !Config If   = OK_SOIL_CARBON_DISCRETIZATION and OK_MOSS
       !Config Def  = 1
       !Config Help = 1=moss_combine ; 2=moss_replace
       !Config Units = [FLAG]
       use_moss_cond_method = moss_combine
       CALL getin_p('USE_MOSS_COND_METHOD', use_moss_cond_method)
       IF ( (use_moss_cond_method /= moss_combine) .AND. (use_moss_cond_method /= moss_replace) ) THEN
          CALL ipslerr_p(3,'thermosoil_initialize', 'Error in variable use_moss_cond_method','USE_MOSS_COND_METHOD must be equal 1 or 2','')
       END IF
    END IF

    IF (printlev >= 3) WRITE (numout,*) 'Start thermosoil_initialize '

    !! 1. Allocate soil temperatures variables
    IF (use_thermo_newdisc) THEN
       ALLOCATE (fSOC(kjpindex,ngrnd),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of fSOC','','')
       
       ALLOCATE (BD_FM(kjpindex,ngrnd),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of BD_FM','','')
       
       ALLOCATE (Phi_CM(kjpindex,ngrnd),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of Phi_CM','','')	
    END IF
    
    ALLOCATE (refsoc(kjpindex,ngrnd),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of refsoc','','')

    ALLOCATE (ptn(kjpindex,ngrnd,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of ptn','','')

    ALLOCATE (ptn_pftmean(kjpindex,ngrnd),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of ptn_pftmean','','')

    ALLOCATE (dz1(ngrnd),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of dz1','','')

    ALLOCATE (cgrnd(kjpindex,ngrnd-1),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of cgrnd','','')

    ALLOCATE (dgrnd(kjpindex,ngrnd-1),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of dgrnd','','')

    ALLOCATE (pcapa(kjpindex,ngrnd),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of pcapa','','')

    ALLOCATE (pkappa(kjpindex,ngrnd),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of pkappa','','')

    ALLOCATE (pkappa_per_pft(kjpindex,ngrnd,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of pkappa_per_pft','','')

    ALLOCATE (pcapa_per_pft(kjpindex,ngrnd,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of pcapa_per_pft','','')

    ALLOCATE (pcapa_en_per_pft(kjpindex,ngrnd,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of pcapa_en_per_pft','','')

    ALLOCATE (pcapa_snow(kjpindex,nsnow),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of pcapa_snow','','')

    ALLOCATE (pkappa_snow(kjpindex,nsnow),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of pkappa_snow','','')
    
    ALLOCATE (pcapa_ice(kjpindex,nice),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of pcapa_ice','','')

    ALLOCATE (pkappa_ice(kjpindex,nice),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of pkappa_ice','','')

    ! Temporary fix: Initialize following variable because they are output to xios before the first calculation
    pcapa  = 0
    pkappa = 0
    pkappa_per_pft = 0
    pcapa_per_pft = 0
    pcapa_en_per_pft = 0
    pcapa_snow  = 0
    pkappa_snow = 0
    pcapa_ice  = 0
    pkappa_ice = 0

    ALLOCATE (surfheat_incr(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of surfheat_incr','','')

    ALLOCATE (coldcont_incr(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of coldcont_incr','','')

    ALLOCATE (pcapa_en(kjpindex,ngrnd),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of pcapa_en','','')
    ! Initialization to zero used at first time step in thermosoil_energy_diag, only for diagnostic variables coldcont_incr and surfheat_incr
    pcapa_en(:,:) = 0.

    ALLOCATE (ptn_beg(kjpindex,ngrnd,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of ptn_beg','','')

    ALLOCATE (temp_sol_beg(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of temp_sol_beg','','')

    ALLOCATE (shum_ngrnd_perma(kjpindex,ngrnd),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of shum_ngrnd_perma','','')
    shum_ngrnd_perma(:,:) = 0

    ALLOCATE (profil_froz(kjpindex,ngrnd),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of profil_froz','','')

    ALLOCATE (profil_froz_pft(kjpindex,ngrnd,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of profil_froz_pft','','')

    IF (ok_freeze_thermix) THEN
       ALLOCATE (pcappa_supp(kjpindex,ngrnd),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of ok_freeze_termix','','')
       ! Initialization to zero used at first time step only for diagnostic output. 
       ! This variable is only used in thermosoil_readajust and always calculated before in thermosoil_getdiff.
       pcappa_supp(:,:) = 0.
    END IF
    IF (ok_Ecorr) THEN
       ALLOCATE (e_soil_lat(kjpindex),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of e_soil_lat','','')
    END IF

    ALLOCATE (dz5(ngrnd),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of dz5','','')

    ALLOCATE (mc_layt(kjpindex,ngrnd),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of mc_layt','','')

    ALLOCATE (mcl_layt(kjpindex,ngrnd),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of mcl_layt','','')

    ALLOCATE (tmc_layt(kjpindex,ngrnd),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of tmc_layt','','')

    ALLOCATE (mc_layt_pft(kjpindex,ngrnd,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of mc_layt_pft','','')

    ALLOCATE (mcl_layt_pft(kjpindex,ngrnd,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of mcl_layt_pft','','')

    ALLOCATE (tmc_layt_pft(kjpindex,ngrnd,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of tmc_layt_pft','','')

    ALLOCATE (mcs(nscm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of mcs','','')

    ALLOCATE (QZ(nscm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of QZ','','')

    ALLOCATE (so_capa_dry(nscm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of so_capa_dry','','')

    ALLOCATE (so_capa_ice(nscm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of so_capa_ice','','')


    !! Soil texture choose : Now useless since njsc defines the dominant texture within 13 classes whichever the soil map
    QZ(:) = QZ_usda(:)
    mcs(:) = mcs_usda(:)

    !Config Key   = DRY_SOIL_HEAT_CAPACITY
    !Config Desc  = Dry soil Heat capacity of soils
    !Config If    = 
    !Config Def   = (1.47, 1.41, 1.34, 1.27, 1.21, 1.21, 1.18, 1.32, 1.23, 1.18, 1.15, 1.09, 1.09)*e+6
    !Config Help  = Values taken from : Pielke [2002, 2013]
    !Config Units = [J.m^{-3}.K^{-1}]
    so_capa_dry(:) = so_capa_dry_usda(:)
    CALL getin_p("DRY_SOIL_HEAT_CAPACITY",so_capa_dry)

    !! Check parameter value (correct range)
    IF ( MINVAL(so_capa_dry(:)) <= zero ) THEN
       CALL ipslerr_p(3, "thermosoil_initialize", &
            &     "Wrong parameter value for DRY_SOIL_HEAT_CAPACITY", &
            &     "These parameters should be positive. ", &
            &     "Please, check parameter value in run.def or orchidee.def. ")
    END IF
    
    
    
    !! 2. Initialize variable from restart file or with default values 
    
    !! Reads restart files for soil temperatures only. If no restart file is
    !! found,  the initial soil temperature is by default set to 280K at all depths. The user
    !! can decide to initialize soil temperatures at an other value, in which case he should set the flag THERMOSOIL_TPRO
    !! to this specific value in the run.def.
    IF (printlev>=3) WRITE (numout,*) 'Read restart file for THERMOSOIL variables'

    CALL ioconf_setatt_p('UNITS', 'K')
    CALL ioconf_setatt_p('LONG_NAME','Soil Temperature profile')
    CALL restget_p (rest_id, 'ptn', nbp_glo, ngrnd, nvm, kjit, .TRUE., ptn, "gather", nbp_glo, index_g)

    ! Initialize ptn if it was not found in restart file
    IF (ALL(ptn(:,:,:)==val_exp)) THEN 
       ! ptn was not found in restart file

       IF (read_reftemp) THEN
          ! Read variable ptn from file
          CALL thermosoil_read_reftempfile(kjpindex,lalo,ptn(:,:,1))
          DO jv = 2,nvm
             ptn(:,:,jv)=ptn(:,:,1)
          ENDDO ! jv = 1,nvm
       ELSE
          ! Initialize ptn with a constant value which can be set in run.def

          !Config Key   = THERMOSOIL_TPRO
          !Config Desc  = Initial soil temperature profile if not found in restart
          !Config Def   = 280.
          !Config If    = OK_SECHIBA
          !Config Help  = The initial value of the temperature profile in the soil if 
          !Config         its value is not found in the restart file. Here
          !Config         we only require one value as we will assume a constant 
          !Config         throughout the column.
          !Config Units = Kelvin [K]
          CALL setvar_p (ptn, val_exp,'THERMOSOIL_TPRO',280._r_std)
       END IF
    END IF

    ! Read refsoc from restart file. If the variable is not found in the restart file,
    ! then read from file for the case use_refsoc and use_thermo_newdisc is false.
    CALL restget_p (rest_id, 'refsoc', nbp_glo, ngrnd, 1, kjit, .TRUE., refsoc, "gather", nbp_glo, index_g)
    IF (ALL(refsoc(:,:) == val_exp)) THEN
       IF (use_refsoc .AND. .NOT. use_thermo_newdisc) THEN
         CALL thermosoil_read_refsoc_file(kjpindex,lalo,neighbours, resolution, contfrac)
       ENDIF
    ENDIF	
	
    ! Read from restart file the variables fSOC, BD_FM and Phi_CM. Initialize all variables if fSOC not found.
    IF (use_thermo_newdisc) THEN
       CALL restget_p (rest_id, 'fSOC', nbp_glo, ngrnd, 1, kjit, .TRUE., fSOC, "gather", nbp_glo, index_g)
       CALL restget_p (rest_id, 'BD_FM', nbp_glo, ngrnd, 1, kjit, .TRUE., BD_FM, "gather", nbp_glo, index_g)
       CALL restget_p (rest_id, 'Phi_CM', nbp_glo, ngrnd, 1, kjit, .TRUE., Phi_CM, "gather", nbp_glo, index_g)
       
       IF (ALL(fSOC(:,:) == val_exp)) THEN
          IF (use_refsoc) THEN
             CALL thermosoil_read_soc_bd_file(kjpindex,lalo,neighbours,resolution, contfrac)
          ELSE IF (impose_soc_bd) THEN
             DO i = 1,kjpindex
                fSOC(i,:) = fSOC_init(:)
                BD_FM(i,:) = BDOD_init(:)
                Phi_CM(i,:) = CFVO_init(:)
             ENDDO
          ENDIF
       END IF
    ENDIF
    
    ! Initialize ptn_beg (variable needed in thermosoil_readadjust called from thermosoil_coef)
    ptn_beg(:,:,:) = ptn(:,:,:)
    
    ! Initialize temp_sol_beg with values from previous time-step
    temp_sol_beg(:) = temp_sol_new(:) 

    ! Read e_soil_lat from restart file or initialize
    IF (ok_Ecorr) THEN
       CALL restget_p (rest_id, 'e_soil_lat', nbp_glo, 1, 1, kjit, .TRUE., &
            e_soil_lat, "gather", nbp_glo, index_g)
       CALL setvar_p (e_soil_lat, val_exp,'NO_KEYWORD',zero)
    END IF

    ! Read gtemp from restart file
    CALL restget_p (rest_id, 'gtemp', nbp_glo, 1, 1, kjit, .TRUE., &
         gtemp, "gather", nbp_glo, index_g)
    CALL setvar_p (gtemp, val_exp,'NO_KEYWORD',zero)
    

    ! Read variables calculated in thermosoil_coef from restart file
    ! If the variables were not found in the restart file, the logical 
    ! calculate_coef will be true and thermosoil_coef will be called further below.
    ! These variables need to be in the restart file to avoid a time shift that
    ! would be done using thermosoil_coef at this stage.
    calculate_coef=.FALSE.
    CALL ioconf_setatt_p('UNITS', 'J m-2 K-1')
    CALL ioconf_setatt_p('LONG_NAME','Apparent surface heat capacity')
    CALL restget_p (rest_id, 'soilcap', nbp_glo, 1, 1, kjit, .TRUE., &
         soilcap, "gather", nbp_glo, index_g)
    IF (ALL(soilcap(:)==val_exp)) calculate_coef=.TRUE.

    CALL ioconf_setatt_p('UNITS', 'W m-2')
    CALL ioconf_setatt_p('LONG_NAME','Apparent soil heat flux')
    CALL restget_p (rest_id, 'soilflx', nbp_glo, 1, 1, kjit, .TRUE., &
         soilflx, "gather", nbp_glo, index_g)
    IF (ALL(soilflx(:)==val_exp)) calculate_coef=.TRUE.

    CALL ioconf_setatt_p('UNITS', '')
    CALL ioconf_setatt_p('LONG_NAME','Integration coefficient for the numerical scheme')
    CALL restget_p (rest_id, 'cgrnd', nbp_glo, ngrnd-1, 1, kjit, .TRUE., &
         cgrnd, "gather", nbp_glo, index_g)
    IF (ALL(cgrnd(:,:)==val_exp)) calculate_coef=.TRUE.

    CALL ioconf_setatt_p('UNITS', '')
    CALL ioconf_setatt_p('LONG_NAME','Integration coefficient for the numerical scheme')
    CALL restget_p (rest_id, 'dgrnd', nbp_glo, ngrnd-1, 1, kjit, .TRUE., &
         dgrnd, "gather", nbp_glo, index_g)
    IF (ALL(dgrnd(:,:)==val_exp)) calculate_coef=.TRUE.

    CALL ioconf_setatt_p('UNITS', '')
    CALL ioconf_setatt_p('LONG_NAME','Integration coefficient for the numerical scheme')
    CALL restget_p (rest_id, 'cgrnd_snow', nbp_glo, nsnow, 1, kjit, .TRUE., &
         cgrnd_snow, "gather", nbp_glo, index_g)
    IF (ALL(cgrnd_snow(:,:)==val_exp)) calculate_coef=.TRUE.

    CALL ioconf_setatt_p('UNITS', '')
    CALL ioconf_setatt_p('LONG_NAME','Integration coefficient for the numerical scheme')
    CALL restget_p (rest_id, 'dgrnd_snow', nbp_glo, nsnow, 1, kjit, .TRUE., &
         dgrnd_snow, "gather", nbp_glo, index_g)
    IF (ALL(dgrnd_snow(:,:)==val_exp)) calculate_coef=.TRUE.

    CALL ioconf_setatt_p('UNITS', '')
    CALL ioconf_setatt_p('LONG_NAME','Coefficient of the linear extrapolation of surface temperature')
    CALL restget_p (rest_id, 'lambda_snow', nbp_glo, 1, 1, kjit, .TRUE., &
         lambda_snow, "gather", nbp_glo, index_g)
    IF (ALL(lambda_snow(:)==val_exp)) calculate_coef=.TRUE.
    
    IF (ok_ice_sheet) THEN
       CALL ioconf_setatt_p('UNITS', '')
       CALL ioconf_setatt_p('LONG_NAME','Integration coefficient for the numerical scheme')
       CALL restget_p (rest_id, 'cgrnd_ice', nbp_glo, nice, 1, kjit, .TRUE., &
            cgrnd_ice, "gather", nbp_glo, index_g)
       IF (ALL(cgrnd_ice(:,:)==val_exp)) calculate_coef=.TRUE.
       
       CALL ioconf_setatt_p('UNITS', '')
       CALL ioconf_setatt_p('LONG_NAME','Integration coefficient for the numerical scheme')
       CALL restget_p (rest_id, 'dgrnd_ice', nbp_glo, nice, 1, kjit, .TRUE., &
            dgrnd_ice, "gather", nbp_glo, index_g)
       IF (ALL(dgrnd_ice(:,:)==val_exp)) calculate_coef=.TRUE.
       
       CALL ioconf_setatt_p('UNITS', '')
       CALL ioconf_setatt_p('LONG_NAME','Coefficient of the linear extrapolation of ice temperature')
       CALL restget_p (rest_id, 'lambda_ice', nbp_glo, 1, 1, kjit, .TRUE., &
            lambda_ice, "gather", nbp_glo, index_g)
       IF (ALL(lambda_ice(:)==val_exp)) calculate_coef=.TRUE.

    END IF ! ok_ice_sheet

    !! 2.2 Computes some physical constants and arrays depending on the soil vertical discretization 

    ! Calculate so_capa_ice
    DO jsc = 1, nscm
       so_capa_ice(jsc) = capa_ice*rho_ice
    END DO
    IF (printlev>=2) WRITE(numout,*) 'Calculation of so_capa_ice(:)=', so_capa_ice(:),' and capa_ice=',capa_ice
    
    ! Computing some usefull constants for the numerical scheme
    ! Use znt(depth of nodes) and zlt(depth of deeper layer interface) from vertical_soil module.  
    DO jg=1,ngrnd-1
      dz1(jg)  = un / (znt(jg+1) - znt(jg))
      dz5(jg) = (zlt(jg) - znt(jg)) * dz1(jg)
    ENDDO
    dz1(ngrnd) = 0.0
    dz5(ngrnd) = 0.0
    lambda = znt(1) * dz1(1)


    CALL ioconf_setatt_p('UNITS', 'K')
    CALL ioconf_setatt_p('LONG_NAME','Soil Temperature profile mean over pft')
    CALL restget_p (rest_id, 'ptn_pftmean', nbp_glo, ngrnd, 1, kjit, .TRUE., ptn_pftmean, "gather", nbp_glo, index_g)
    IF (ALL(ptn_pftmean(:,:)==val_exp)) THEN
       ! Initialize ptn_pftmean if not found in restart file
       IF (ok_soil_carbon_discretization) THEN
          ptn_pftmean(:,:)=0.0
          DO m=1,nvm
             DO jg = 1, ngrnd
                ptn_pftmean(:,jg) = ptn_pftmean(:,jg) + ptn(:,jg,m) * veget_max(:,m)
             ENDDO
          ENDDO
       ELSE
          ! For this case, ptn is constant over all pfts. Use here PFT=1 to initialize ptn_pftmean.
          ptn_pftmean(:,:) = ptn(:,:,1)
       END IF
    END IF
    
    ! Send out the temperature profile on the first nslm levels(the levels treated in hydrol)
    stempdiag(:,:) = ptn_pftmean(:,1:nslm)
    ftempdiag(:,:) = ptn_pftmean(:,1:ngrnd)
    
    !! Initialize output arguments to be used in sechiba
    ptnlev1(:) = ptn_pftmean(:,1)

    !! 2.3. Computes cgrnd, dgrnd, soilflx and soilcap coefficients only if they were not found in restart file.
    IF (calculate_coef) THEN
       ! Interpolate variables needed by thermosoil_coef to the thermal levels
       CALL thermosoil_humlev(kjpindex, shumdiag_perma, mc_layh, mcl_layh, tmc_layh, &
                              mc_layh_pft, mcl_layh_pft, tmc_layh_pft)

       IF (printlev>=3) WRITE (numout,*) 'thermosoil_coef will be called in the intialization phase'
       CALL thermosoil_coef (&
            kjpindex,      temp_sol_new,    snow, njsc, soiltile, &
            frac_snow_veg, frac_snow_nobio, totfrac_nobio,        &
            snowdz,        snowrho,         snowtemp,       pb,   &
            veget_max,     som_total,       depth_organic_soil,  &
            ptn, ptn_pftmean,                                     &
            soilcap,       soilflx,         cgrnd,          dgrnd,&
            lambda_snow,   cgrnd_snow,      dgrnd_snow,           &
            lambda_ice,    cgrnd_ice,       dgrnd_ice,            &
            icetemp,       icedz,           ice_sheet_mask)

    END IF

  END SUBROUTINE thermosoil_initialize


!! ================================================================================================================================
!! SUBROUTINE   : thermosoil_main
!!
!>\BRIEF        Thermosoil_main computes the soil thermal properties and dynamics, ie solves
!! the heat diffusion equation within the soil. 
!!
!! DESCRIPTION : The resolution of the soil heat diffusion equation 
!! relies on a numerical finite-difference implicit scheme
!! fully described in the reference and in the header of the thermosoil module.
!! - The dependency of the previous timestep hidden in the 
!! integration coefficients cgrnd and dgrnd (EQ1), calculated in thermosoil_coef, and 
!! called at the end of the routine to prepare for the next timestep.
!! - The effective computation of the new soil temperatures is performed in thermosoil_profile. 
!!
!! - thermosoil_coef calculates the coefficients for the numerical scheme for the very first iteration of thermosoil;
!! after that, thermosoil_coef is called only at the end of the module to calculate the coefficients for the next timestep.
!! - thermosoil_profile solves the numerical scheme.\n
!!
!! - Flags : one unique flag : THERMOSOIL_TPRO (to be set to the desired initial soil in-depth temperature in K; by default 280K)
!!
!! RECENT CHANGE(S) : Change vertical discretization (consistent with hydrology layers) and soil thermal properties (taking into account soil texture effects).
!!
!! MAIN OUTPUT VARIABLE(S): vertically discretized soil temperatures ptn, soil
!! thermal properties (pcapa, pkappa), apparent surface heat capacity (soilcap)
!! and heat flux (soilflx) to be used in enerbil at the next timestep to solve
!! the surface energy balance.
!!
!! REFERENCE(S) : 
!! - Hourdin, F. (1992). Study and numerical simulation of the general circulation of planetary atmospheres,
!!  Ph.D. thesis, Paris VII University. Remark: the part of F. Hourdin' s PhD thesis relative to the thermal
!!  integration scheme has been scanned and is provided along with the documentation, with name : 
!!  Hourdin_1992_PhD_thermal_scheme.pdf
!!
!! FLOWCHART    : 
!! \latexonly
!! \includegraphics[scale = 1]{thermosoil_flowchart.png}
!! \endlatexonly
!! 
!! \ A flag to activate the heat production by soil microbial activityn
!_ ================================================================================================================================

  SUBROUTINE thermosoil_main (kjit, kjpindex, &
       index, indexgrnd, &
       temp_sol_new, snow, soilcap, soilflx, &
       shumdiag_perma, stempdiag, ftempdiag, ptnlev1, rest_id, hist_id, hist2_id, &
       snowdz,snowrho,snowtemp,gtemp,pb,&
       mc_layh, mcl_layh, tmc_layh,  mc_layh_pft, mcl_layh_pft, tmc_layh_pft, njsc, soiltile, &
       depth_organic_soil, heat_Zimov, deeptemp_prof, deephum_prof,&
       som_total, veget_max, frac_snow_veg,frac_snow_nobio,totfrac_nobio,temp_sol_add, &
       lambda_snow, cgrnd_snow, dgrnd_snow, &
       lambda_ice,  cgrnd_ice,  dgrnd_ice, icetemp, icedz, ice_sheet_mask)

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                            :: kjit             !! Time step number (unitless) 
    INTEGER(i_std), INTENT(in)                            :: kjpindex         !! Domain size (unitless)
    INTEGER(i_std),INTENT (in)                            :: rest_id,hist_id  !! Restart_ file and history file identifier 
                                                                              !! (unitless)
    INTEGER(i_std),INTENT (in)                            :: hist2_id         !! history file 2 identifier (unitless)
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)      :: index            !! Indeces of the points on the map (unitless)
    INTEGER(i_std),DIMENSION (kjpindex*ngrnd), INTENT (in):: indexgrnd        !! Indeces of the points on the 3D map (vertical 
                                                                              !! dimension towards the ground) (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)      :: temp_sol_new     !! Surface temperature at the present time-step,
                                                                              !! Ts @tex ($K$) @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)         :: snow             !! Snow mass @tex ($kg$) @endtex.
                                                                              !! Caveat: when there is snow on the
                                                                              !! ground, the snow is integrated into the soil for
                                                                              !! the calculation of the thermal dynamics. It means
                                                                              !! that the uppermost soil layers can completely or 
                                                                              !! partially consist in snow. In the second case, zx1
                                                                              !! and zx2 are the fraction of the soil layer 
                                                                              !! consisting in snow and 'normal' soil, respectively
                                                                              !! This is calculated in thermosoil_coef.
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (in)    :: shumdiag_perma   !! Soil saturation degree (0-1, unitless)
    REAL(r_std), DIMENSION(kjpindex),INTENT (in)          :: depth_organic_soil !! Depth at which there is still organic matter (m) 
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm), INTENT (in)   :: heat_Zimov   !! Heating associated with decomposition [W/m**3 soil]
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm,nelements), INTENT (in) :: som_total  !! Total soil organic matter for use in thermal calcs (g/m**3)
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT (in)    :: veget_max        !! Fraction of vegetation type
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT (in)    :: snowdz           !! Snow depth [m]
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT (in)    :: snowrho          !! Snow density (Kg/m^3)
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT (inout) :: snowtemp         !! Snow temperature (K)
    REAL(r_std), DIMENSION (kjpindex),INTENT (in)         :: pb               !! Surface presure (hPa)
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (in)    :: mc_layh          !! Volumetric soil moisture content for each layer in hydrol at nodes(liquid + ice) (m3/m3)
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (in)    :: mcl_layh         !! Volumetric soil moisture content for each layer in hydrol at nodes(liquid) (m3/m3)
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (in)    :: tmc_layh         !! Total soil moisture content for each layer in hydrol(liquid + ice) (mm)
    REAL(r_std),DIMENSION (kjpindex,nslm,nvm), INTENT (in):: mc_layh_pft      !! Volumetric soil moisture content for each layer in hydrol at nodes(liquid + ice) (m3/m3)
    REAL(r_std),DIMENSION (kjpindex,nslm,nvm), INTENT (in):: mcl_layh_pft     !! Volumetric soil moisture content for each layer in hydrol at nodes(liquid) (m3/m3)
    REAL(r_std),DIMENSION (kjpindex,nslm,nvm), INTENT (in):: tmc_layh_pft     !! Total soil moisture content for each layer in hydrol(liquid + ice) (mm)
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)      :: njsc             !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless)
    REAL(r_std),DIMENSION (kjpindex,nstm), INTENT (in)    :: soiltile         !! Fraction of each soil tile within vegtot (0-1, unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)          :: frac_snow_veg    !! Snow cover fraction on vegeted area
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT(in)   :: frac_snow_nobio  !! Snow cover fraction on non-vegeted area
    REAL(r_std),DIMENSION (kjpindex),INTENT(in)           :: totfrac_nobio    !! Total fraction of continental ice+lakes+cities+...
                                                                              !!(unitless,0-1)
                                                                              !! @tex $(kg m^{-2})$ @endtex 
    REAL(r_std),DIMENSION (kjpindex,nice),INTENT (in)     :: icedz            !! Ice depth
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)      :: ice_sheet_mask   !! Ice sheet mask
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)      :: temp_sol_add     !! additional surface temperature due to the melt of first layer
                                                                              !! at the present time-step @tex ($K$) @endtex

    !! 0.2 Output variables

    REAL(r_std),DIMENSION (kjpindex), INTENT (out)        :: ptnlev1          !! 1st level soil temperature   
    REAL(r_std),DIMENSION (kjpindex),INTENT(out)          :: gtemp            !! First soil layer temperature
    REAL(r_std), DIMENSION (kjpindex,ngrnd,nvm), INTENT (out) :: deephum_prof !! moisture on a deep thermodynamic profile for permafrost calcs
    REAL(r_std), DIMENSION (kjpindex,ngrnd,nvm), INTENT (out) :: deeptemp_prof!! temp on a deep thermodynamic profile for permafrost calcs

    !! 0.3 Modified variables

    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)      :: soilcap          !! apparent surface heat capacity considering snow and soil surface
                                                                              !! @tex ($J m^{-2} K^{-1}$) @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)      :: soilflx          !! apparent soil heat flux  considering snow and soil surface
                                                                              !! @tex ($W m^{-2}$) @endtex
                                                                              !! , positive 
                                                                              !! towards the soil, writen as Qg (ground heat flux) 
                                                                              !! in the history files, and computed at the end of 
                                                                              !! thermosoil for the calculation of Ts in enerbil, 
                                                                              !! see EQ3.
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (out)   :: stempdiag        !! temperature profile @tex ($K$) @endtex
    REAL(r_std),DIMENSION (kjpindex,ngrnd), INTENT (out)  :: ftempdiag        !! temperature profile @tex ($K$) @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT(inout)       :: lambda_snow      !! Coefficient of the linear extrapolation of surface temperature 
    REAL(r_std),DIMENSION (kjpindex,nsnow), INTENT (inout):: cgrnd_snow       !! Integration coefficient for snow numerical scheme
    REAL(r_std),DIMENSION (kjpindex,nsnow), INTENT (inout):: dgrnd_snow       !! Integration coefficient for snow numerical scheme
    
    REAL(r_std),DIMENSION (kjpindex), INTENT(inout)       :: lambda_ice       !! Coefficient of the linear extrapolation of ice temperature 
    REAL(r_std),DIMENSION (kjpindex,nice), INTENT (inout) :: cgrnd_ice        !! Integration coefficient for snow numerical scheme
    REAL(r_std),DIMENSION (kjpindex,nice), INTENT (inout) :: dgrnd_ice        !! Integration coefficient for snow numerical scheme
    REAL(r_std),DIMENSION (kjpindex,nice),INTENT (in)     :: icetemp          !! Ice temperature (K)

    !! 0.4 Local variables

    INTEGER(i_std)                                        :: jv,ji,ii, jg, m
    REAL(r_std),DIMENSION (kjpindex)                      :: snowtemp_weighted!! Snow temperature weighted by snow density, only for diag (K)
    REAL(r_std),DIMENSION (kjpindex, nsnow)               :: pkappa_snow_diag !! Only for diag, containing xios_default_val
    REAL(r_std),DIMENSION (kjpindex, nsnow)               :: pcapa_snow_diag  !! Only for diag, containing xios_default_val
    REAL(r_std),DIMENSION (kjpindex, nsnow)               :: snowtemp_diag    !! Only for diag, containing xios_default_val
    LOGICAL                                               :: ok_zimov         !! A flag to activate the heat production by soil microbial activity
!_ ================================================================================================================================
    
  !! 3. Put the soil wetness diagnostic on the levels of the soil temperature

    !!?? this could logically be put just before the last call to
    !!thermosoil_coef, as the results are used there...
    CALL thermosoil_humlev(kjpindex, shumdiag_perma, mc_layh, mcl_layh, tmc_layh, &
                            mc_layh_pft, mcl_layh_pft, tmc_layh_pft)
    
  !! 4. Effective computation of the soil temperatures profile.
  !!    cgrnd and dgrnd have been calculated in thermosoil_coef at the previous time step 
  !!    but they are correct for the actual time-step.
    CALL thermosoil_profile (kjpindex,      temp_sol_new,                   &
                             frac_snow_veg, frac_snow_nobio, totfrac_nobio, &
                             ptn,           ptn_pftmean,     stempdiag,   snowtemp,      &
                             cgrnd_snow, dgrnd_snow, icetemp, cgrnd_ice,    &
                             dgrnd_ice, ice_sheet_mask)


  !! 5. Call to thermosoil_energy_diag for calculation of diagnostic variables
    CALL thermosoil_energy_diag(kjpindex, temp_sol_new, soilcap)

  !! Save ptn at current stage, to be used in thermosoil_readjust
    ptn_beg(:,:,:) = ptn(:,:,:)

  !! 6. Writing the history files according to the ALMA standards (or not..)

    ! Add XIOS default value where no snow
    DO ji=1,kjpindex 
       IF (snow(ji) .GT. zero) THEN
          pkappa_snow_diag(ji,:) = pkappa_snow(ji,:)
          pcapa_snow_diag(ji,:) = pcapa_snow(ji,:)
          snowtemp_diag(ji,:) = snowtemp(ji,:)
       ELSE
          pkappa_snow_diag(ji,:) = xios_default_val
          pcapa_snow_diag(ji,:) = xios_default_val
          snowtemp_diag(ji,:) = xios_default_val
       END IF
    END DO

    DO ji=1,kjpindex 
       ! Use min_sechiba instead of zero to avoid problem with division by zero
       IF (snow(ji) .GT. min_sechiba) THEN
          snowtemp_weighted(ji) = SUM(snowtemp(ji,:)*snowrho(ji,:))/SUM(snowrho(ji,:))
       ELSE
          snowtemp_weighted(ji) = xios_default_val
       END IF
    END DO
    CALL xios_orchidee_send_field("snowtemp_weighted",snowtemp_weighted)
    CALL xios_orchidee_send_field("ptn_pftmean",ptn_pftmean)
    CALL xios_orchidee_send_field("soilflx",soilflx)
    CALL xios_orchidee_send_field("surfheat_incr",surfheat_incr)
    CALL xios_orchidee_send_field("coldcont_incr",coldcont_incr)
    CALL xios_orchidee_send_field("pkappa",pkappa)
    CALL xios_orchidee_send_field("pkappa_snow",pkappa_snow_diag)
    CALL xios_orchidee_send_field("pcapa",pcapa)
    CALL xios_orchidee_send_field("pcapa_snow",pcapa_snow_diag)
    CALL xios_orchidee_send_field("snowtemp",snowtemp_diag)
    IF (ok_freeze_thermix) CALL xios_orchidee_send_field("shum_ngrnd_perma", shum_ngrnd_perma)
    IF (ok_ice_sheet) THEN
      CALL xios_orchidee_send_field("pkappa_ice",pkappa_ice)
      CALL xios_orchidee_send_field("pcapa_ice",pcapa_ice)
      CALL xios_orchidee_send_field("icetemp",icetemp)
    ENDIF

    IF ( .NOT. almaoutput ) THEN
      CALL histwrite_p(hist_id, 'ptn', kjit, ptn_pftmean, kjpindex*ngrnd, indexgrnd)
      CALL histwrite_p(hist_id, 'Qg', kjit, soilflx, kjpindex, index)
      CALL histwrite_p(hist_id, 'pkappa', kjit, pkappa, kjpindex*ngrnd, indexgrnd)
      CALL histwrite_p(hist_id, 'pcapa', kjit, pcapa, kjpindex*ngrnd, indexgrnd)

      IF (ok_freeze_thermix) THEN
         CALL histwrite_p(hist_id, 'profil_froz', kjit, profil_froz, kjpindex*ngrnd, indexgrnd)
         CALL histwrite_p(hist_id, 'pcappa_supp', kjit, pcappa_supp, kjpindex*ngrnd, indexgrnd)
      END IF
      CALL histwrite_p(hist_id, 'shum_ngrnd_perma', kjit, shum_ngrnd_perma(:,:), kjpindex*ngrnd, indexgrnd)
      
    ELSE
      CALL histwrite_p(hist_id, 'SoilTemp', kjit, ptn_pftmean, kjpindex*ngrnd, indexgrnd)
      CALL histwrite_p(hist_id, 'Qg', kjit, soilflx, kjpindex, index)
      CALL histwrite_p(hist_id, 'DelSurfHeat', kjit, surfheat_incr, kjpindex, index)
      CALL histwrite_p(hist_id, 'DelColdCont', kjit, coldcont_incr, kjpindex, index)
    ENDIF
    IF ( hist2_id > 0 ) THEN
       IF ( .NOT. almaoutput ) THEN
          CALL histwrite_p(hist2_id, 'ptn', kjit, ptn_pftmean, kjpindex*ngrnd, indexgrnd)
       ELSE
          CALL histwrite_p(hist2_id, 'SoilTemp', kjit, ptn_pftmean, kjpindex*ngrnd, indexgrnd)
          CALL histwrite_p(hist2_id, 'Qg', kjit, soilflx, kjpindex, index)
          CALL histwrite_p(hist2_id, 'DelSurfHeat', kjit, surfheat_incr, kjpindex, index)
          CALL histwrite_p(hist2_id, 'DelColdCont', kjit, coldcont_incr, kjpindex, index)
       ENDIF
    ENDIF

  !! 7. Considering the heat released by microbial respiration
    ok_zimov=.FALSE.
    IF (ok_zimov) THEN
       CALL thermosoil_add_heat_zimov(kjpindex, veget_max, ptn, heat_zimov)
    END IF
    
  !! 8. A last final call to thermosoil_coef
 
    !! A last final call to thermosoil_coef, which calculates the different
    !!coefficients (cgrnd, dgrnd, soilcap, soilflx) from this time step to be
    !!used at the next time step, either in the surface temperature calculation
    !!(soilcap, soilflx) or in the soil thermal numerical scheme.
    CALL thermosoil_coef (&
         kjpindex,      temp_sol_new,  snow, njsc, soiltile, &
         frac_snow_veg, frac_snow_nobio, totfrac_nobio,      &
         snowdz,        snowrho,         snowtemp,     pb,   &
         veget_max,     som_total,       depth_organic_soil,&
         ptn,           ptn_pftmean,                         &
         soilcap,       soilflx,         cgrnd,        dgrnd,&
         lambda_snow,   cgrnd_snow,      dgrnd_snow,         &
         lambda_ice, cgrnd_ice, dgrnd_ice, icetemp,          &
         icedz, ice_sheet_mask)
         
    ! Save variables for explicit snow model
    gtemp(:) = ptn_pftmean(:,1)

    !! Initialize output arguments to be used in sechiba
    ptnlev1(:) = ptn_pftmean(:,1)

    !! Initialize output arguments to be returned to sechiba and further used in stomate
    DO jv= 1, nvm
       deephum_prof(:,:,jv)  = shum_ngrnd_perma(:,:)
    END DO
    deeptemp_prof = ptn

    !! Surface temperature is forced to zero celcius if its value is larger than melting point, only for explicit snow scheme

    !! Write the new temperature in the full in a the diagnostic variable.
    ftempdiag(:,:) = ptn_pftmean(:,1:ngrnd)

    !! Surface temperature is forced to zero celcius if its value is larger than melting point
    DO ji=1,kjpindex
       IF  (SUM(snowdz(ji,:)) .GT. 0.0) THEN
          IF (temp_sol_new(ji) .GE. tp_00) THEN
             IF ( ok_ice_sheet .and. (ice_sheet_mask(ji).GT.0)) THEN ! Ice-sheet  point
                temp_sol_new(ji) = tp_00*frac_snow_veg(ji)*(1-totfrac_nobio(ji))+ &      ! weights related to first layer of snow cover fraction on vegetation  
                     tp_00*SUM(frac_snow_nobio(ji,:))*totfrac_nobio(ji)+ &               ! weights related to SCF on nobio
                     tp_00*(1 - SUM(frac_snow_nobio(ji,:)))*totfrac_nobio(ji)+ &         ! weights related to non snow area on nobio (ice area)  
                     temp_sol_new(ji)*(1-frac_snow_veg(ji))*(1-totfrac_nobio(ji))        ! weights related to non snow on vegetation (non ice area)

             ELSE ! no ice (standard point)
                temp_sol_new(ji) = tp_00*frac_snow_veg(ji)*(1-totfrac_nobio(ji))+ &      ! weights related to first layer of snow cover fraction on vegetation  
                     tp_00*SUM(frac_snow_nobio(ji,:))*totfrac_nobio(ji)+ &               ! weights related to SCF on nobio
                     temp_sol_new(ji)*(1-(frac_snow_veg(ji)*(1-totfrac_nobio(ji))+SUM(frac_snow_nobio(ji,:))*totfrac_nobio(ji))) ! weights related to non snow 
             ENDIF
          ENDIF
       END IF

    END DO


    IF (printlev>=3) WRITE (numout,*) ' thermosoil_main done '

  END SUBROUTINE thermosoil_main

  !!  =============================================================================================================================
  !! SUBROUTINE		 		    : thermosoil_finalize
  !!
  !>\BRIEF                                    Write to restart file
  !!
  !! DESCRIPTION			    : This subroutine writes the module variables and variables calculated in thermosoil
  !!                                          to restart file
  !! \n
  !_ ==============================================================================================================================
  SUBROUTINE thermosoil_finalize (kjit,    kjpindex, rest_id,   gtemp, &
                                  soilcap, soilflx, lambda_snow, cgrnd_snow, dgrnd_snow,&
                                  lambda_ice, cgrnd_ice, dgrnd_ice)

    !! 0. Variable and parameter declaration
    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                            :: kjit             !! Time step number (unitless) 
    INTEGER(i_std), INTENT(in)                            :: kjpindex         !! Domain size (unitless)
    INTEGER(i_std),INTENT (in)                            :: rest_id          !! Restart file identifier(unitless)
    REAL(r_std),DIMENSION (kjpindex),INTENT(in)           :: gtemp            !! First soil layer temperature
    REAL(r_std),DIMENSION (kjpindex),INTENT(in)           :: soilcap
    REAL(r_std),DIMENSION (kjpindex),INTENT(in)           :: soilflx
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)         :: lambda_snow      !! Coefficient of the linear extrapolation of surface temperature 
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT (in)  :: cgrnd_snow       !! Integration coefficient for snow numerical scheme
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT (in)  :: dgrnd_snow       !! Integration coefficient for snow numerical scheme
    
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)         :: lambda_ice      !! Coefficient of the linear extrapolation of ice temperature 
    REAL(r_std), DIMENSION (kjpindex,nice), INTENT (in)  :: cgrnd_ice       !! Integration coefficient for snow numerical scheme
    REAL(r_std), DIMENSION (kjpindex,nice), INTENT (in)  :: dgrnd_ice       !! Integration coefficient for snow numerical scheme

    !! 0.2 Local variables
    ! To store variables names for I/O
    CHARACTER(LEN=80) :: var_name
    CHARACTER(LEN=10) :: part_str
    INTEGER(i_std)    :: iv
!_ ================================================================================================================================
    
    !! 1. Write variables to restart file to be used for the next simulation
    IF (printlev>=3) WRITE (numout,*) 'Write restart file with THERMOSOIL variables'

    CALL restput_p(rest_id, 'ptn', nbp_glo, ngrnd, nvm, kjit, ptn, 'scatter', nbp_glo, index_g)
    !-
    CALL restput_p(rest_id, 'ptn_pftmean', nbp_glo, ngrnd, 1, kjit, ptn_pftmean, 'scatter', nbp_glo, index_g)
    !-
    CALL restput_p (rest_id, 'refsoc', nbp_glo, ngrnd, 1, kjit, refsoc, 'scatter', nbp_glo, index_g)

    IF (use_thermo_newdisc) THEN
       CALL restput_p(rest_id, 'fSOC', nbp_glo, ngrnd, 1, kjit, fSOC, 'scatter', nbp_glo, index_g)
       CALL restput_p(rest_id, 'BD_FM', nbp_glo, ngrnd, 1, kjit, BD_FM, 'scatter', nbp_glo, index_g)
       CALL restput_p(rest_id, 'Phi_CM', nbp_glo, ngrnd, 1, kjit, Phi_CM, 'scatter', nbp_glo, index_g)    
    END IF
    
    IF (ok_Ecorr) THEN
       CALL restput_p(rest_id, 'e_soil_lat', nbp_glo, 1 , 1, kjit, e_soil_lat, 'scatter', nbp_glo, index_g)
       !-
    END IF
    
    CALL restput_p(rest_id, 'gtemp', nbp_glo, 1, 1, kjit, gtemp, 'scatter', nbp_glo, index_g)
    !-
    CALL restput_p(rest_id, 'soilcap', nbp_glo, 1, 1, kjit, soilcap, 'scatter', nbp_glo, index_g)
    !-
    CALL restput_p(rest_id, 'soilflx', nbp_glo, 1, 1, kjit, soilflx, 'scatter', nbp_glo, index_g)
    !-
    CALL restput_p(rest_id, 'cgrnd', nbp_glo, ngrnd-1, 1, kjit, cgrnd, 'scatter', nbp_glo, index_g)
    !-
    CALL restput_p(rest_id, 'dgrnd', nbp_glo, ngrnd-1, 1, kjit, dgrnd, 'scatter', nbp_glo, index_g)
    !-
    CALL restput_p(rest_id, 'cgrnd_snow', nbp_glo, nsnow, 1, kjit, cgrnd_snow, 'scatter', nbp_glo, index_g)
    !-
    CALL restput_p(rest_id, 'dgrnd_snow', nbp_glo, nsnow, 1, kjit, dgrnd_snow, 'scatter', nbp_glo, index_g)
    !-
    CALL restput_p(rest_id, 'lambda_snow', nbp_glo, 1, 1, kjit, lambda_snow, 'scatter', nbp_glo, index_g)

    IF ( ok_ice_sheet ) THEN
      CALL restput_p(rest_id, 'cgrnd_ice', nbp_glo, nice, 1, kjit, cgrnd_ice, 'scatter', nbp_glo, index_g)
      CALL restput_p(rest_id, 'dgrnd_ice', nbp_glo, nice, 1, kjit, dgrnd_ice, 'scatter', nbp_glo, index_g)
      CALL restput_p(rest_id, 'lambda_ice', nbp_glo, 1, 1, kjit, lambda_ice, 'scatter', nbp_glo, index_g)
    ENDIF

  END SUBROUTINE thermosoil_finalize


!! ================================================================================================================================
!! SUBROUTINE   : thermosoil_clear
!!
!>\BRIEF        Deallocates the allocated arrays.
!! The call of thermosoil_clear originates from sechiba_clear but the calling sequence and 
!! its purpose require further investigation.
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S) : None 
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : None 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

 SUBROUTINE thermosoil_clear()

        IF ( ALLOCATED (ptn)) DEALLOCATE (ptn)
        IF ( ALLOCATED (ptn_pftmean)) DEALLOCATE (ptn_pftmean)
        IF ( ALLOCATED (cgrnd)) DEALLOCATE (cgrnd) 
        IF ( ALLOCATED (dgrnd)) DEALLOCATE (dgrnd) 
        IF ( ALLOCATED (pcapa)) DEALLOCATE (pcapa)
        IF ( ALLOCATED (pkappa))  DEALLOCATE (pkappa)
        IF ( ALLOCATED (pkappa_per_pft))  DEALLOCATE (pkappa_per_pft)
        IF ( ALLOCATED (pcapa_per_pft))  DEALLOCATE (pcapa_per_pft)
        IF ( ALLOCATED (pcapa_en_per_pft))  DEALLOCATE (pcapa_en_per_pft)
        IF ( ALLOCATED (pcapa_snow)) DEALLOCATE (pcapa_snow)
        IF ( ALLOCATED (pkappa_snow))  DEALLOCATE (pkappa_snow)
        IF ( ALLOCATED (pcapa_ice)) DEALLOCATE (pcapa_ice)
        IF ( ALLOCATED (pkappa_ice))  DEALLOCATE (pkappa_ice)
        IF ( ALLOCATED (pcapa_en)) DEALLOCATE (pcapa_en)
        IF ( ALLOCATED (ptn_beg)) DEALLOCATE (ptn_beg)
        IF ( ALLOCATED (temp_sol_beg)) DEALLOCATE (temp_sol_beg)
        IF ( ALLOCATED (surfheat_incr)) DEALLOCATE (surfheat_incr)
        IF ( ALLOCATED (coldcont_incr)) DEALLOCATE (coldcont_incr)
        IF ( ALLOCATED (shum_ngrnd_perma)) DEALLOCATE (shum_ngrnd_perma)
        IF ( ALLOCATED (profil_froz)) DEALLOCATE (profil_froz)
        IF ( ALLOCATED (mc_layt)) DEALLOCATE (mc_layt)
        IF ( ALLOCATED (mcl_layt)) DEALLOCATE (mcl_layt)
        IF ( ALLOCATED (tmc_layt)) DEALLOCATE (tmc_layt)
        IF ( ALLOCATED (mc_layt_pft)) DEALLOCATE (mc_layt_pft)
        IF ( ALLOCATED (mcl_layt_pft)) DEALLOCATE (mcl_layt_pft)
        IF ( ALLOCATED (tmc_layt_pft)) DEALLOCATE (tmc_layt_pft)
        IF ( ALLOCATED (profil_froz_pft)) DEALLOCATE (profil_froz_pft)
        IF ( ALLOCATED (refsoc)) DEALLOCATE (refsoc)
        IF ( ALLOCATED (fSOC)) DEALLOCATE (fSOC)
        IF ( ALLOCATED (BD_FM)) DEALLOCATE (BD_FM)
        IF ( ALLOCATED (Phi_CM)) DEALLOCATE (Phi_CM)

  END SUBROUTINE thermosoil_clear


!! ================================================================================================================================
!! SUBROUTINE   : thermosoil_coef
!!
!>\BRIEF        Calculate soil thermal properties, integration coefficients, apparent heat flux,
!! surface heat capacity,  
!!
!! DESCRIPTION	: This routine computes : \n
!!		1. the soil thermal properties. \n 
!!		2. the integration coefficients of the thermal numerical scheme, cgrnd and dgrnd,
!!              which depend on the vertical grid and on soil properties, and are used at the next 
!!              timestep.\n
!!              3. the soil apparent heat flux and surface heat capacity (soilflx
!!              and soilcap), used by enerbil to compute the surface temperature at the next
!!              timestep.\n
!!             -  The soil thermal properties depend on water content (shum_ngrnd_perma, shumdiag_perma, 
!!              mc_layt, mcl_layt, tmc_layt), dominant soil texture(njsc), and on the presence 
!!              of snow : snow is integrated into the soil for the thermal calculations, ie if there 
!!              is snow on the ground, the first thermal layer(s) consist in snow, depending on the 
!!              snow-depth. If a layer consists out of snow and soil, wheighed soil properties are 
!!              calculated\n
!!             - The coefficients cgrnd and dgrnd are the integration
!!              coefficients for the thermal scheme \n
!!                              T(k+1)=cgrnd(k)+dgrnd(k)*T(k) \n
!!                                      -- EQ1 -- \n
!!              They correspond respectively to $\beta$ and $\alpha$ from F. Hourdin\'s thesis and 
!!              their expression can be found in this document (eq A19 and A20)
!!             - soilcap and soilflx are the apparent surface heat capacity and flux
!!               used in enerbil at the next timestep to solve the surface
!!               balance for Ts (EQ3); they correspond to $C_s$ and $F_s$ in F.
!!               Hourdin\'s PhD thesis and are expressed in eq. A30 and A31. \n
!!                 soilcap*(Ts(t)-Ts(t-1))/dt=soilflx+otherfluxes(Ts(t)) \n
!!                                      -- EQ3 --\n
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): cgrnd, dgrnd, pcapa, pkappa, soilcap, soilflx
!!
!! REFERENCE(S) :
!! - Hourdin, F. (1992). Study and numerical simulation of the general circulation of planetary atmospheres,
!! Ph.D. thesis, Paris VII University. Remark: the part of F. Hourdin's PhD thesis relative to the thermal
!! integration scheme has been scanned and is provided along with the documentation, with name : 
!! Hourdin_1992_PhD_thermal_scheme.pdf
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE thermosoil_coef (kjpindex,temp_sol_new, snow,njsc,soiltile, &
       frac_snow_veg, frac_snow_nobio, totfrac_nobio,        &
       snowdz,        snowrho,         snowtemp,       pb,   &
       veget_max,     som_total,       depth_organic_soil,  &
       ptn,           ptn_pftmean,                           &
       soilcap,       soilflx,         cgrnd,          dgrnd,&
       lambda_snow,   cgrnd_snow,      dgrnd_snow,           &
       lambda_ice,    cgrnd_ice,       dgrnd_ice,   icetemp, &
       icedz, ice_sheet_mask)

    !! 0. Variables and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                             :: kjpindex     !! Domain size (unitless)
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)         :: temp_sol_new !! soil surface temperature @tex ($K$) @endtex
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)         :: snow         !! snow mass @tex ($Kg$) @endtex
    INTEGER(i_std), DIMENSION (kjpindex), INTENT (in)      :: njsc         !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless)
    REAL(r_std),DIMENSION (kjpindex,nstm), INTENT (in)     :: soiltile         !! Fraction of each soil tile within vegtot (0-1, unitless)
    !
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)           :: frac_snow_veg   !! Snow cover fraction on vegeted area
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT(in)    :: frac_snow_nobio !! Snow cover fraction on non-vegeted area
    REAL(r_std),DIMENSION (kjpindex),INTENT(in)            :: totfrac_nobio   !! Total fraction of continental ice+lakes+cities+...
    !!(unitless,0-1)
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(in)    :: snowdz          !! Snow depth (m)
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(in)    :: snowrho         !! Snow density (Kg/m^3)
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(in)    :: snowtemp        !! Snow temperature (K)
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)         :: pb              !! Surface presure (hPa)
    REAL(r_std), DIMENSION (kjpindex,nice), INTENT(in)     :: icetemp         !! Ice temperature (K)
    REAL(r_std), DIMENSION (kjpindex,nice), INTENT(in)     :: icedz           !! Ice depth (m)
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)       :: ice_sheet_mask !! Ice sheet mask
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(in)      :: veget_max       !! Fraction of vegetation type
    REAL(r_std), DIMENSION (kjpindex,ngrnd,nvm,nelements), INTENT (in) :: som_total       !! total soil carbon for use in thermal calcs (g/m**3)
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)         :: depth_organic_soil !! Depth at which there is still organic matter (m) 

    !! 0.2 Output variables

    REAL(r_std), DIMENSION (kjpindex), INTENT (out)        :: soilcap      !! surface heat capacity considering snow and soil surface
    !! @tex ($J m^{-2} K^{-1}$) @endtex
    REAL(r_std), DIMENSION (kjpindex), INTENT (out)        :: soilflx      !! surface heat flux considering snow and soil surface @tex ($W m^{-2}$) @endtex,
    !! positive towards the 
    !! soil, writen as Qg (ground heat flux) in the history 
    !! files.
    REAL(r_std), DIMENSION (kjpindex,ngrnd-1), INTENT(out) :: cgrnd        !! matrix coefficient for the computation of soil 
    !! temperatures (beta in F. Hourdin thesis)
    REAL(r_std), DIMENSION (kjpindex,ngrnd-1), INTENT(out) :: dgrnd        !! matrix coefficient for the computation of soil 
    !! temperatures (alpha in F. Hourdin thesis)


    !! 0.3 Modified variable

    REAL(r_std), DIMENSION (kjpindex,ngrnd,nvm), INTENT (inout):: ptn      !! vertically discretized soil temperatures per pft. ptn is only modified if ok_Ecorr.
    REAL(r_std), DIMENSION (kjpindex,ngrnd), INTENT (inout):: ptn_pftmean  !! vertically discretized soil temperatures. ptn is only modified if ok_Ecorr.
    REAL(r_std), DIMENSION (kjpindex), INTENT(inout)       :: lambda_snow  !! Coefficient of the linear extrapolation of surface temperature 
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT (inout):: cgrnd_snow   !! Integration coefficient for snow numerical scheme
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT (inout):: dgrnd_snow   !! Integration coefficient for snow numerical scheme

    REAL(r_std), DIMENSION (kjpindex), INTENT(inout)       :: lambda_ice   !! Coefficient of the linear extrapolation of ice temperature 
    REAL(r_std), DIMENSION (kjpindex,nice), INTENT (inout) :: cgrnd_ice    !! Integration coefficient for snow numerical scheme
    REAL(r_std), DIMENSION (kjpindex,nice), INTENT (inout) :: dgrnd_ice    !! Integration coefficient for snow numerical scheme

    !! 0.4 Local variables

    INTEGER(i_std)                                         :: ji, jg, m, jst
    REAL(r_std), DIMENSION (kjpindex,ngrnd-1)              :: zdz1         !! numerical (buffer) constant 
    !! @tex ($W m^{-1} K^{-1}$) @endtex
    REAL(r_std), DIMENSION (kjpindex,ngrnd)                :: zdz2         !! numerical (buffer) constant  
    !! @tex ($W m^{-1} K^{-1}$) @endtex
    REAL(r_std), DIMENSION (kjpindex)                      :: z1           !! numerical constant @tex ($W m^{-1} K^{-1}$) @endtex
    REAL(r_std), DIMENSION (kjpindex)                      :: soilcap_nosnow      !! surface heat capacity
    !! @tex ($J m^{-2} K^{-1}$)
    !! @endtex
    REAL(r_std), DIMENSION (kjpindex)                      :: soilflx_nosnow      !! surface heat flux @tex ($W m^{-2}$) @endtex,
    !! positive towards the soil, written as Qg
    !!(ground heat flux in the history files).
    REAL(r_std), DIMENSION (kjpindex)                      :: snowcap             !! apparent snow heat capacity @tex ($J m^{-2} K^{-1}$)
    REAL(r_std), DIMENSION (kjpindex)                      :: snowflx             !! apparent snow-atmosphere heat flux @tex ($W m^{-2}$) @endtex
    REAL(r_std), DIMENSION (kjpindex)                      :: icecap              !! apparent ice heat capacity @tex ($J m^{-2} K^{-1}$)
    REAL(r_std), DIMENSION (kjpindex)                      :: iceflx              !! apparent ice-atmosphere heat flux @tex ($W m^{-2}$) @endtex
    REAL(r_std), DIMENSION (kjpindex,nsnow)                :: dz1_snow
    REAL(r_std), DIMENSION (kjpindex,nsnow)                :: ZSNOWDZM
    REAL(r_std), DIMENSION (kjpindex,nsnow)                :: dz2_snow
    REAL(r_std), DIMENSION (kjpindex,nsnow)                :: zdz1_snow
    REAL(r_std), DIMENSION (kjpindex,nsnow)                :: zdz2_snow
    REAL(r_std), DIMENSION (kjpindex)                      :: z1_snow
    REAL(r_std), DIMENSION (kjpindex,nice)                 :: dz1_ice
    REAL(r_std), DIMENSION (kjpindex,nice)                 :: dz2_ice
    REAL(r_std), DIMENSION (kjpindex,nice)                 :: zdz1_ice
    REAL(r_std), DIMENSION (kjpindex,nice)                 :: zdz2_ice
    REAL(r_std), DIMENSION (kjpindex)                      :: z1_ice

    REAL(r_std), DIMENSION (kjpindex)                      :: snowflxtot          !! Total snow flux (including snow on vegetated and bare soil and nobio areas)
    !! @tex ($W m^{-2}$) @endtex
    !! positive towards the soil

    !_ ================================================================================================================================

    !! 1. Computation of the soil thermal properties
    !     Modify mcs if peatland exists
    IF (ok_peat_hydro .AND. (nstm == 4)) THEN
       DO ji=1,kjpindex
          ! note that jst here correspond to the index of the principal soil texture used for this grid
          ! this is different from jst in hydrol.f90, which correspond to soil tiles. xw
          jst = njsc(ji)
          
          ! the case for agripeat and tides is not yet considered here. xw 202602
          IF ((soiltile (ji,4) .GT. min_sechiba)) THEN
             mcs(jst) = mcs(jst)*(un-soiltile(ji,4))+mcs_peat*soiltile(ji,4)
          ENDIF
       ENDDO
    END IF

    ! Computation of the soil thermal properties; snow properties are also accounted for
    IF (ok_freeze_thermix) THEN
       IF (use_soilc_insulation .OR. ok_soil_carbon_discretization) THEN
          pcapa(:,:) = zero
          pcapa_en(:,:) = zero
          pkappa(:,:) = zero
          ! possible incoherence when activating peatland: organic porosity poros_org (0.92) is used, 
          ! which is different from mcs_peat (0.9) used in hydrol.f90
          ! gridcell quantity of organic matter may be different as well from peatland.
		  
          IF (use_thermo_newdisc) THEN
             CALL thermosoil_getdiff_pft_newdisc( kjpindex,  ptn,     njsc,     shum_ngrnd_perma, &
                  depth_organic_soil, som_total, snowrho, snowtemp, pb, veget_max)
          ELSE
             CALL thermosoil_getdiff_pft( kjpindex,  ptn,     njsc,     shum_ngrnd_perma, &
                  depth_organic_soil, som_total, snowrho, snowtemp, pb, veget_max, &
                  icetemp, ice_sheet_mask)
          END IF
         
          DO m=1,nvm
             DO jg = 1, ngrnd
                pcapa(:,jg) = pcapa(:,jg) + pcapa_per_pft(:,jg,m) * veget_max(:,m)
                pcapa_en(:,jg) = pcapa_en(:,jg) + pcapa_en_per_pft(:,jg,m) * veget_max(:,m)
                pkappa(:,jg) = pkappa(:,jg) + pkappa_per_pft(:,jg,m) * veget_max(:,m)
             END DO
          END DO
       ELSE
          CALL thermosoil_getdiff( kjpindex, snow, ptn_pftmean, njsc, snowrho, snowtemp, pb, totfrac_nobio, &
               icetemp, ice_sheet_mask )
       ENDIF
    ELSE
       ! Special case without soil freezing
       CALL thermosoil_getdiff_old_thermix_without_snow( kjpindex, njsc, snowrho, snowtemp, pb, icetemp, &
            ice_sheet_mask)
    ENDIF

    ! Energy conservation : Correction to make sure that the same latent heat is released and 
    ! consumed during freezing and thawing
    IF (ok_Ecorr) THEN
       CALL thermosoil_readjust(kjpindex, ptn, ptn_pftmean, veget_max)
    ENDIF


    !! 2. Computation of the coefficients of the numerical integration scheme for the soil layers

    !! 2.1 Calculate numerical coefficients zdz1 and zdz2
    DO jg=1,ngrnd
       DO ji=1,kjpindex
          zdz2(ji,jg)=pcapa(ji,jg) * dlt(jg)/dt_sechiba
       ENDDO
    ENDDO

    DO jg=1,ngrnd-1
       DO ji=1,kjpindex
          zdz1(ji,jg) = dz1(jg) * pkappa(ji,jg)
       ENDDO
    ENDDO

    !! 2.2 Calculate coefficients cgrnd and dgrnd for soil
    DO ji = 1,kjpindex
       z1(ji) = zdz2(ji,ngrnd) + zdz1(ji,ngrnd-1)
       cgrnd(ji,ngrnd-1) = zdz2(ji,ngrnd) * ptn_pftmean(ji,ngrnd) / z1(ji)
       dgrnd(ji,ngrnd-1) = zdz1(ji,ngrnd-1) / z1(ji)
    ENDDO

    DO jg = ngrnd-1,2,-1
       DO ji = 1,kjpindex
          z1(ji) = un / (zdz2(ji,jg) + zdz1(ji,jg-1) + zdz1(ji,jg) * (un - dgrnd(ji,jg)))
          cgrnd(ji,jg-1) = (ptn_pftmean(ji,jg) * zdz2(ji,jg) + zdz1(ji,jg) * cgrnd(ji,jg)) * z1(ji)
          dgrnd(ji,jg-1) = zdz1(ji,jg-1) * z1(ji)
       ENDDO
    ENDDO

    !! 3. Computation of the coefficients of the numerical integration scheme for the snow layers

    !! 3.1 Calculate numerical coefficients zdz1_snow, zdz2_snow and lambda_snow

    DO ji = 1, kjpindex
       IF ( ok_ice_sheet .and. (ice_sheet_mask(ji).GT.0)) THEN   ! Ice-sheet  point => temp in ice 
          ! For Ice
          DO jg = 1, nice-1
             dz1_ice(ji,jg)  = 2.0 / (icedz(ji,jg+1)+icedz(ji,jg))
          ENDDO
          lambda_ice(ji) = icedz(ji,1)/2.0 * dz1_ice(ji,1)
          DO jg=1,nice
             zdz2_ice(ji,jg)=pcapa_ice(ji,jg) * icedz(ji,jg)/dt_sechiba
          ENDDO
          DO jg=1,nice-1
             zdz1_ice(ji,jg) = dz1_ice(ji,jg) * pkappa_ice(ji,jg)
          ENDDO

          ! the bottom ice
          zdz1_ice(ji,nice) = pkappa_ice(ji,nice) / ( zlt(1) + icedz(ji,nice)/2 )

          ! For Snow
          ! Calculate internal values
          DO jg = 1, nsnow
             ZSNOWDZM(ji,jg) = MAX(snowdz(ji,jg),psnowdzmin)
          ENDDO
          dz2_snow(ji,:)=ZSNOWDZM(ji,:)


          DO jg = 1, nsnow-1
             dz1_snow(ji,jg)  = 2.0 / (dz2_snow(ji,jg+1)+dz2_snow(ji,jg))
          ENDDO

          lambda_snow(ji) = dz2_snow(ji,1)/2.0 * dz1_snow(ji,1)

          DO jg=1,nsnow
             zdz2_snow(ji,jg)=pcapa_snow(ji,jg) * dz2_snow(ji,jg)/dt_sechiba
          ENDDO

          DO jg=1,nsnow-1
             zdz1_snow(ji,jg) = dz1_snow(ji,jg) * pkappa_snow(ji,jg)
          ENDDO

          ! the bottom snow on ice
          zdz1_snow(ji,nsnow) = pkappa_snow(ji,nsnow) / ( icedz(ji,1) + dz2_snow(ji,nsnow)/2 )

       ELSE ! standard point (no ice)
          ! Calculate internal values
          DO jg = 1, nsnow
             ZSNOWDZM(ji,jg) = MAX(snowdz(ji,jg),psnowdzmin)
          ENDDO
          dz2_snow(ji,:)=ZSNOWDZM(ji,:)

          DO jg = 1, nsnow-1
             dz1_snow(ji,jg)  = 2.0 / (dz2_snow(ji,jg+1)+dz2_snow(ji,jg))
          ENDDO

          lambda_snow(ji) = dz2_snow(ji,1)/2.0 * dz1_snow(ji,1)

          DO jg=1,nsnow
             zdz2_snow(ji,jg)=pcapa_snow(ji,jg) * dz2_snow(ji,jg)/dt_sechiba
          ENDDO

          DO jg=1,nsnow-1
             zdz1_snow(ji,jg) = dz1_snow(ji,jg) * pkappa_snow(ji,jg)
          ENDDO

          ! the bottom snow
          zdz1_snow(ji,nsnow) = pkappa_snow(ji,nsnow) / ( zlt(1) + dz2_snow(ji,nsnow)/2 )
          lambda_ice(ji) = lambda
       ENDIF
    ENDDO

      

    !! 3.2 Calculate coefficients cgrnd_snow and dgrnd_snow for snow

    DO ji = 1,kjpindex
       IF ( ok_ice_sheet .and. (ice_sheet_mask(ji).GT.0)) THEN ! Ice-sheet  point => temp in ice
          ! ice bottom level
          z1_ice(ji) = zdz2(ji,1)+(un-dgrnd(ji,1))*zdz1(ji,1)+zdz1_ice(ji,nice)
          cgrnd_ice(ji,nice) = (zdz2(ji,1) * ptn_pftmean(ji,1) + zdz1(ji,1) * cgrnd(ji,1) ) / z1_ice(ji)
          dgrnd_ice(ji,nice) = zdz1_ice(ji,nice) / z1_ice(ji)

          ! ice next-to-bottom level
          z1_ice(ji) = zdz2_ice(ji,nice)+(un-dgrnd_ice(ji,nice))*zdz1_ice(ji,nice)+zdz1_ice(ji,nice-1)
          cgrnd_ice(ji,nice-1) = (zdz2_ice(ji,nice)*icetemp(ji,nice) +&
               zdz1_ice(ji,nice)*cgrnd_ice(ji,nice))/z1_ice(ji)
          dgrnd_ice(ji,nice-1) = zdz1_ice(ji,nice-1) / z1_ice(ji)

          DO jg = nice-1,2,-1
             z1_ice(ji) = un / (zdz2_ice(ji,jg) + zdz1_ice(ji,jg-1) + zdz1_ice(ji,jg) * (un - dgrnd_ice(ji,jg)))
             cgrnd_ice(ji,jg-1) = (icetemp(ji,jg) * zdz2_ice(ji,jg) + zdz1_ice(ji,jg) * cgrnd_ice(ji,jg)) * z1_ice(ji)
             dgrnd_ice(ji,jg-1) = zdz1_ice(ji,jg-1) * z1_ice(ji)
          ENDDO

          ! snow bottom level
          z1_snow(ji) = zdz2_ice(ji,1)+(un-dgrnd_ice(ji,1))*zdz1_ice(ji,1)+zdz1_snow(ji,nsnow)
          cgrnd_snow(ji,nsnow) = (zdz2_ice(ji,1) * icetemp(ji,1) + zdz1_ice(ji,1) * cgrnd_ice(ji,1) ) / z1_snow(ji)
          dgrnd_snow(ji,nsnow) = zdz1_snow(ji,nsnow) / z1_snow(ji)

          ! snow next-to-bottom level
          z1_snow(ji) = zdz2_snow(ji,nsnow)+(un-dgrnd_snow(ji,nsnow))*zdz1_snow(ji,nsnow)+zdz1_snow(ji,nsnow-1)
          cgrnd_snow(ji,nsnow-1) = (zdz2_snow(ji,nsnow)*snowtemp(ji,nsnow) +&
               zdz1_snow(ji,nsnow)*cgrnd_snow(ji,nsnow))/z1_snow(ji)
          dgrnd_snow(ji,nsnow-1) = zdz1_snow(ji,nsnow-1) / z1_snow(ji)

          DO jg = nsnow-1,2,-1
             z1_snow(ji) = un / (zdz2_snow(ji,jg) + zdz1_snow(ji,jg-1) + zdz1_snow(ji,jg) * (un - dgrnd_snow(ji,jg)))
             cgrnd_snow(ji,jg-1) = (snowtemp(ji,jg) * zdz2_snow(ji,jg) + zdz1_snow(ji,jg) * cgrnd_snow(ji,jg)) * z1_snow(ji)
             dgrnd_snow(ji,jg-1) = zdz1_snow(ji,jg-1) * z1_snow(ji)
          ENDDO

       ELSE ! Standard point (no ice)
          ! bottom level
          z1_snow(ji) = zdz2(ji,1)+(un-dgrnd(ji,1))*zdz1(ji,1)+zdz1_snow(ji,nsnow)
          cgrnd_snow(ji,nsnow) = (zdz2(ji,1) * ptn_pftmean(ji,1) + zdz1(ji,1) * cgrnd(ji,1) ) / z1_snow(ji)
          dgrnd_snow(ji,nsnow) = zdz1_snow(ji,nsnow) / z1_snow(ji)

          ! next-to-bottom level
          z1_snow(ji) = zdz2_snow(ji,nsnow)+(un-dgrnd_snow(ji,nsnow))*zdz1_snow(ji,nsnow)+zdz1_snow(ji,nsnow-1)
          cgrnd_snow(ji,nsnow-1) = (zdz2_snow(ji,nsnow)*snowtemp(ji,nsnow)+&
               zdz1_snow(ji,nsnow)*cgrnd_snow(ji,nsnow))/z1_snow(ji)
          dgrnd_snow(ji,nsnow-1) = zdz1_snow(ji,nsnow-1) / z1_snow(ji)

          DO jg = nsnow-1,2,-1
             z1_snow(ji) = un / (zdz2_snow(ji,jg) + zdz1_snow(ji,jg-1) + zdz1_snow(ji,jg) * (un - dgrnd_snow(ji,jg)))
             cgrnd_snow(ji,jg-1) = (snowtemp(ji,jg) * zdz2_snow(ji,jg) + zdz1_snow(ji,jg) * cgrnd_snow(ji,jg)) * z1_snow(ji)
             dgrnd_snow(ji,jg-1) = zdz1_snow(ji,jg-1) * z1_snow(ji)
          ENDDO

          cgrnd_ice(ji,:) = cgrnd(ji,1)
          dgrnd_ice(ji,:) = dgrnd(ji,1)
       ENDIF
    ENDDO




    !! 4. Computation of the apparent ground heat flux 

    !! Computation of apparent snow-atmosphere flux
    DO ji = 1,kjpindex
       snowflx(ji) = zdz1_snow(ji,1) * (cgrnd_snow(ji,1) + (dgrnd_snow(ji,1)-1.) * snowtemp(ji,1))
       snowcap(ji) = (zdz2_snow(ji,1) * dt_sechiba + dt_sechiba * (un - dgrnd_snow(ji,1)) * zdz1_snow(ji,1))
       z1_snow(ji) = lambda_snow(ji) * (un - dgrnd_snow(ji,1)) + un 
       snowcap(ji) = snowcap(ji) / z1_snow(ji)
       snowflx(ji) = snowflx(ji) + snowcap(ji) * (snowtemp(ji,1) * z1_snow(ji) - &
            & lambda_snow(ji) * cgrnd_snow(ji,1) - temp_sol_new(ji)) / dt_sechiba

       IF ( ok_ice_sheet .and. ice_sheet_mask(ji).GT.0 ) THEN ! Ice-sheet  point => temp in ice
          iceflx(ji) = zdz1_ice(ji,1) * (cgrnd_ice(ji,1) + (dgrnd_ice(ji,1)-1.) * icetemp(ji,1))
          icecap(ji) = (zdz2_ice(ji,1) * dt_sechiba + dt_sechiba * (un - dgrnd_ice(ji,1)) * zdz1_ice(ji,1))
          z1_ice(ji) = lambda_ice(ji) * (un - dgrnd_ice(ji,1)) + un 
          icecap(ji) = icecap(ji) / z1_ice(ji)
          iceflx(ji) = iceflx(ji) + icecap(ji) * (icetemp(ji,1) * z1_ice(ji) - &
               & lambda_ice(ji) * cgrnd_ice(ji,1) - temp_sol_new(ji)) / dt_sechiba
       ENDIF
    ENDDO



    !! Computation of the apparent ground heat flux (> towards the soil) and
    !! apparent surface heat capacity, used at the next timestep by enerbil to
    !! compute the surface temperature.
    DO ji = 1,kjpindex
       IF ( ok_ice_sheet .and. (ice_sheet_mask(ji).GT.0)) THEN !
          !Ice-sheet  point
          soilflx_nosnow(ji) = iceflx(ji)
          soilcap_nosnow(ji) = icecap(ji)
       ELSE ! Standard point (no ice)
          soilflx_nosnow(ji) = zdz1(ji,1) * (cgrnd(ji,1) + (dgrnd(ji&
               &,1)-1.) * ptn_pftmean(ji,1))
          soilcap_nosnow(ji) = (zdz2(ji,1) * dt_sechiba + dt_sechiba &
               &* (un - dgrnd(ji,1)) * zdz1(ji,1))
          z1(ji) = lambda * (un - dgrnd(ji,1)) + un
          soilcap_nosnow(ji) = soilcap_nosnow(ji) / z1(ji)
          soilflx_nosnow(ji) = soilflx_nosnow(ji) +&
               & soilcap_nosnow(ji) * (ptn_pftmean(ji,1) * z1(ji) -&
               & lambda * cgrnd(ji,1) - temp_sol_new(ji)) / dt_sechiba
       ENDIF
    ENDDO

    !! Add snow fraction
    ! Using an effective heat capacity and heat flux by a simple pondering of snow and soil fraction
    DO ji = 1, kjpindex
       soilcap(ji) = snowcap(ji)*frac_snow_veg(ji)*(1-totfrac_nobio(ji))+    & ! weights related to snow cover fraction on vegetation  
            snowcap(ji)*SUM(frac_snow_nobio(ji,:))*totfrac_nobio(ji)+ &        ! weights related to SCF on nobio
            soilcap_nosnow(ji)*(1-(frac_snow_veg(ji)*(1-totfrac_nobio(ji))+SUM(frac_snow_nobio(ji,:))*totfrac_nobio(ji))) ! weights related to non snow fraction
       soilflx(ji) = snowflx(ji)*frac_snow_veg(ji)*(1-totfrac_nobio(ji))+    & ! weights related to snow cover fraction on vegetation 
            snowflx(ji)*SUM(frac_snow_nobio(ji,:))*totfrac_nobio(ji)+ & ! weights related to SCF on nobio
            soilflx_nosnow(ji)*(1-(frac_snow_veg(ji)*(1-totfrac_nobio(ji))+SUM(frac_snow_nobio(ji,:))*totfrac_nobio(ji))) ! weights related to non snow fraction
    ENDDO


    ! Total snow flux (including snow on vegetated and bare soil and nobio areas)
    DO ji = 1, kjpindex
       snowflxtot(ji) = snowflx(ji)*frac_snow_veg(ji)*(1-totfrac_nobio(ji)) + &
            soilflx_nosnow(ji)*SUM(frac_snow_nobio(ji,:))*totfrac_nobio(ji)
    ENDDO

    CALL xios_orchidee_send_field("snowflxtot",snowflxtot(:))

    IF (printlev>=3) WRITE (numout,*) ' thermosoil_coef done '

  END SUBROUTINE thermosoil_coef
 
 
!! ================================================================================================================================
!! SUBROUTINE   : thermosoil_profile
!!
!>\BRIEF        In this routine solves the numerical soil thermal scheme, ie calculates the new soil temperature profile; 
!! 
!!
!! DESCRIPTION	: The calculation of the new soil temperature profile is based on
!! the cgrnd and dgrnd values from the previous timestep and the surface temperature Ts aka temp_sol_new. (see detailed
!! explanation in the header of the thermosoil module or in the reference).\n
!!                              T(k+1)=cgrnd(k)+dgrnd(k)*T(k)\n
!!                                      -- EQ1 --\n
!!                           Ts=(1+lambda)*T(1) -lambda*T(2)\n 
!!                                      -- EQ2--\n
!!
!! RECENT CHANGE(S) : None
!! 
!! MAIN OUTPUT VARIABLE(S): ptn (soil temperature profile on the thermal axis), 
!!                          stempdiag (soil temperature profile on the diagnostic axis)
!!
!! REFERENCE(S) :
!! - Hourdin, F. (1992). Study and numerical simulation of the general circulation of planetary atmospheres,
!! Ph.D. thesis, Paris VII University. Remark: the part of F. Hourdin's PhD thesis relative to the thermal
!! integration scheme has been scanned and is provided along with the documentation, with name : 
!! Hourdin_1992_PhD_thermal_scheme.pdf
!!
!! FLOWCHART    : None 
!! \n 
!_ ================================================================================================================================

  SUBROUTINE thermosoil_profile (kjpindex,      temp_sol_new,                   &
       frac_snow_veg, frac_snow_nobio, totfrac_nobio, &
       ptn,           ptn_pftmean,     stempdiag,       snowtemp,      &
       cgrnd_snow, dgrnd_snow, icetemp, cgrnd_ice,    &
       dgrnd_ice, ice_sheet_mask)

    !! 0. Variables and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                               :: kjpindex       !! Domain size (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: temp_sol_new   !! Surface temperature at the present time-step 
    !! @tex ($K$) @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)             :: frac_snow_veg  !! Snow cover fraction on vegeted area
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT(in)      :: frac_snow_nobio!! Snow cover fraction on non-vegeted area
    REAL(r_std),DIMENSION (kjpindex),INTENT(in)              :: totfrac_nobio  !! Total fraction of continental ice+lakes+cities+...
    !! (unitless,0-1)
    REAL(r_std),DIMENSION (kjpindex,nsnow), INTENT(in)       :: snowtemp       !! Snow temperature (K)
    REAL(r_std),DIMENSION (kjpindex,nsnow), INTENT(in)       :: cgrnd_snow     !! Integration coefficient for snow numerical scheme
    REAL(r_std),DIMENSION (kjpindex,nsnow), INTENT(in)       :: dgrnd_snow     !! Integration coefficient for snow numerical scheme
    REAL(r_std), DIMENSION (kjpindex,nice), INTENT(in)       :: icetemp        !! Ice temperature (K)
    REAL(r_std),DIMENSION (kjpindex,nice), INTENT(in)        :: cgrnd_ice      !! Integration coefficient for ice numerical scheme
    REAL(r_std),DIMENSION (kjpindex,nice), INTENT(in)        :: dgrnd_ice      !! Integration coefficient for ice numerical scheme
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)         :: ice_sheet_mask !! Ice sheet mask

    !! 0.3 Modified variables


    !! 0.2 Output variables
    REAL(r_std),DIMENSION (kjpindex,ngrnd,nvm), INTENT (out) :: ptn            !! vertically discretized soil temperatures per pft (K)
    !! @tex ($K$) @endtex
    REAL(r_std),DIMENSION (kjpindex,ngrnd), INTENT (out)     :: ptn_pftmean    !! vertically discretized soil temperatures (K)
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (out)      :: stempdiag      !! diagnostic temperature profile 
    !! @tex ($K$) @endtex

    !! 0.4 Local variables
    INTEGER(i_std)                                           :: ji, jg, jv
    REAL(r_std)                                              :: temp_sol_eff   !! effective surface temperature including snow and soil

    !_ ================================================================================================================================

    !! 1. Computes the soil temperatures ptn.

    !! 1.1. ptn(jg=1) using EQ1 and EQ2

    DO ji = 1,kjpindex
       IF ( ok_ice_sheet .and. ice_sheet_mask(ji).GT.0 )  THEN ! Ice-sheet  point
          temp_sol_eff = snowtemp(ji,nsnow)*frac_snow_veg(ji)*(1-totfrac_nobio(ji)) + &         ! weights related to snow cover fraction on vegetation  
               icetemp(ji,nice)*totfrac_nobio(ji) + &                                           ! weights related to SCF on nobio
               temp_sol_new(ji)*(1-(frac_snow_veg(ji)*(1-totfrac_nobio(ji))+totfrac_nobio(ji))) ! weights related to non snow fraction

          ! Soil temperature calculation with explicit snow if there is snow on the ground
          ptn_pftmean(ji,1) = cgrnd_ice(ji,nice) + dgrnd_ice(ji,nice) * temp_sol_eff

       ELSE ! Standard point (no ice)
          temp_sol_eff=snowtemp(ji,nsnow)*frac_snow_veg(ji)*(1-totfrac_nobio(ji)) + &           ! weights related to snow cover fraction on vegetation  
               snowtemp(ji,nsnow)*SUM(frac_snow_nobio(ji,:))*totfrac_nobio(ji) + &              ! weights related to SCF on nobio
               temp_sol_new(ji)*(1-(frac_snow_veg(ji)*(1-totfrac_nobio(ji)) + &                 ! weights related to non snow fraction on vegetation
               SUM(frac_snow_nobio(ji,:))*totfrac_nobio(ji)))                                   ! weights related to non snow fraction on nobio

          ! Soil temperature calculation with explicit snow if there is snow on the ground
          ptn_pftmean(ji,1) = cgrnd_snow(ji,nsnow) + dgrnd_snow(ji,nsnow) * temp_sol_eff
       ENDIF
    ENDDO

    !! 1.2. ptn_pftmean(jg=2:ngrnd) using EQ1.
    DO jg = 1,ngrnd-1
       DO ji = 1,kjpindex
          ptn_pftmean(ji,jg+1) = cgrnd(ji,jg) + dgrnd(ji,jg) * ptn_pftmean(ji,jg)
       ENDDO
    ENDDO

    !! Calculate ptn per pft
    !  Here ptn is the same at all pfts when ok_soil_carbon_discretization is activated or not. 
    !  ptn is modified in thermosoil_add_heat_zimov if ok_zimov=TRUE.
    DO jv=1,nvm
       ptn(:,:,jv) = ptn_pftmean(:,:)
    END DO


    !! 2. Assigne the soil temperature to the output variable. It is already on the right axis. 
    stempdiag(:,:) = ptn_pftmean(:,1:nslm)

    IF (printlev>=3) WRITE (numout,*) ' thermosoil_profile done '

  END SUBROUTINE thermosoil_profile



!================================================================================================================================
!! SUBROUTINE   : thermosoil_cond
!!
!>\BRIEF          Calculate soil thermal conductivity.  
!!
!! DESCRIPTION  : This routine computes soil thermal conductivity
!!                Code introduced from NOAH LSM. 
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): cnd
!!
!! REFERENCE(S) :
!!    Farouki, O.T.,1986: Thermal Properties of Soils. Series on Rock
!!            and Soil Mechanics, Vol. 11, Trans Tech, 136 PP.
!!    Johansen, O., 1975: Thermal Conductivity of Soils. Ph.D. Thesis,
!!            University of Trondheim,
!!    Peters-Lidard, C. D., Blackburn, E., Liang, X., & Wood, E. F.,
!!            1998: The effect of soil thermal conductivity 
!!            Parameterization on Surface Energy fluxes
!!            and Temperatures. J. of The Atmospheric Sciences,
!!            Vol. 55, pp. 1209-1224.
!! Modify histroy:
!!
!! FLOWCHART    : None
!! \n
!_
!================================================================================================================================

  SUBROUTINE thermosoil_cond (kjpindex, njsc, smc, qz, sh2o, cnd)

    !! 0. Variables and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                                 :: kjpindex      !! Domain size (unitless)
    INTEGER(i_std), DIMENSION (kjpindex), INTENT (in)          :: njsc          !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless)
    REAL(r_std), DIMENSION (kjpindex,ngrnd), INTENT(IN)        :: smc           !! Volumetric Soil Moisture Content (m3/m3)
    REAL(r_std), DIMENSION (nscm), INTENT(IN)                  :: qz            !! Quartz Content (Soil Type Dependent) (0-1)               
    REAL(r_std), DIMENSION (kjpindex,ngrnd), INTENT(IN)        :: sh2o          !! Unfrozen Soil Moisture Content; Frozen Soil Moisture = smc - sh2o

    !! 0.2 Output variables
    REAL(r_std), DIMENSION (kjpindex,ngrnd), INTENT(OUT)       :: cnd           !! Soil Thermal Conductivity (W/m/k)

    !! 0.3 Modified variables

    !! 0.4 Local variables
    REAL(r_std), DIMENSION (kjpindex,ngrnd)                    :: ake           !! Kersten Number (unitless)
    REAL(r_std), DIMENSION (kjpindex,ngrnd)                    :: thksat        !! Saturated Thermal Conductivity (W/m/k)
    REAL(r_std), DIMENSION (kjpindex,ngrnd)                    :: satratio      !! Degree of Saturation (0-1)
    REAL(r_std), DIMENSION (kjpindex,ngrnd)                    :: xu            !! Unfrozen Volume For Saturation (0-1)
    REAL(r_std), DIMENSION (kjpindex,ngrnd)                    :: xunfroz       !! Unfrozon Volume Fraction (0-1)
    REAL(r_std)                                                :: thko          !! Thermal Conductivity for Other Ssoil Components (W/m/k)
    REAL(r_std)                                                :: gammd         !! Dry Dendity (kg/m3)
    REAL(r_std)                                                :: thkdry        !! Dry Thermal Conductivity (W/m/k)
    REAL(r_std)                                                :: thks          !! Thermal Conductivity for the Solids Combined (Quartz + Other) (W/m/k)
    REAL(r_std), PARAMETER                                     :: THKICE = 2.2  !! Ice Thermal Conductivity (W/m/k)
    REAL(r_std), PARAMETER                                     :: THKQTZ = 7.7  !! Thermal Conductivity for Quartz (W/m/k)
    REAL(r_std), PARAMETER                                     :: THKW = 0.57   !! Water Thermal Conductivity (W/m/k)
    INTEGER(i_std)                                             :: ji, jg, jst

    !_================================================================================================================================

    !! 1. Dry and Saturated Thermal Conductivity.

    DO ji = 1,kjpindex
       jst = njsc(ji)

       !! 1.1. Dry density (Kg/m3) and Dry thermal conductivity (W.M-1.K-1)
       gammd = (1. - mcs(jst))*2700.
       thkdry = (0.135* gammd+ 64.7)/ (2700. - 0.947* gammd)

       !! 1.2. thermal conductivity of "other" soil components
       IF (qz(jst) > 0.2) THEN
          thko = 2.0
       ELSEIF (qz(jst) <= 0.2) THEN
          thko = 3.0
       ENDIF

       !! 1.3. Thermal conductivity of solids
       thks = (THKQTZ ** qz(jst))* (thko ** (1. - qz(jst)))

       DO jg = 1,ngrnd      
          !! 1.4. saturation ratio
          satratio(ji,jg) = smc(ji,jg) / mcs(jst)

          !! 1.5. Saturated Thermal Conductivity (thksat)
          IF ( smc(ji,jg) > min_sechiba ) THEN
             xunfroz(ji,jg) = sh2o(ji,jg) / smc(ji,jg)   ! Unfrozen Fraction (From i.e., 100%Liquid, to 0. (100% Frozen))
             xu(ji,jg) = xunfroz(ji,jg) * mcs(jst)  ! Unfrozen volume for saturation (porosity*xunfroz)
             thksat(ji,jg) = thks ** (1. - mcs(jst))* THKICE ** (mcs(jst) - xu(ji,jg))* THKW ** (xu(ji,jg))
          ELSE
             ! this value will not be used since ake=0 for this case
             thksat(ji,jg)=0 
          END IF
       END DO ! DO jg = 1,ngrnd

       !! 2. Kersten Number (ake)
       DO jg = 1,ngrnd
          IF ( (sh2o(ji,jg) + 0.0005) <  smc(ji,jg) ) THEN
             ! Frozen
             ake(ji,jg) = satratio(ji,jg)
          ELSE
             ! Unfrozen
             ! Eq 11 in Peters-Lidard et al., 1998
             IF ( satratio(ji,jg) >  0.1 ) THEN
                IF (jst < 4 )  THEN
                   ! Coarse 
                   ake(ji,jg) = 0.7 * LOG10 (SATRATIO(ji,jg)) + 1.0
                ELSE
                   ! Fine 
                   ake(ji,jg) = LOG10 (satratio(ji,jg)) + 1.0
                ENDIF
             ELSEIF ( satratio(ji,jg) >  0.05 .AND. satratio(ji,jg) <=  0.1 ) THEN
                IF (jst < 4 )  THEN
                   ! Coarse 
                   ake(ji,jg) = 0.7 * LOG10 (satratio(ji,jg)) + 1.0
                ELSE
                   ! Fine 
                   ake(ji,jg) = 0.0
                ENDIF
             ELSE
                ake(ji,jg) = 0.0  ! use k = kdry
             END IF
          END IF
       END DO ! DO jg = 1,ngrnd

       !! 3. Thermal conductivity (cnd)
       DO jg = 1,ngrnd
          cnd(ji,jg) = ake(ji,jg) * (thksat(ji,jg) - thkdry) + thkdry
       END DO ! DO jg = 1,ngrnd

    END DO !DO ji = 1,kjpindex

  END SUBROUTINE thermosoil_cond

!================================================================================================================================
!! SUBROUTINE   : thermosoil_cond_pft
!!
!>\BRIEF          Calculate soil thermal conductivity.  
!!
!! DESCRIPTION  : This routine computes soil thermal conductivity
!!                but considers the fact that soil organic carbon can decrease
!!                conductivity 
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): cnd
!!
!! REFERENCE(S) :
!!    Farouki, O.T.,1986: Thermal Properties of Soils. Series on Rock
!!            and Soil Mechanics, Vol. 11, Trans Tech, 136 PP.
!!    Johansen, O., 1975: Thermal Conductivity of Soils. Ph.D. Thesis,
!!            University of Trondheim,
!!    Peters-Lidard, C. D., Blackburn, E., Liang, X., & Wood, E. F.,
!!            1998: The effect of soil thermal conductivity 
!!            Parameterization on Surface Energy fluxes
!!            and Temperatures. J. of The Atmospheric Sciences,
!!            Vol. 55, pp. 1209-1224.
!!    Lawrence and Slater,2008: Incorporating organic soil into a global climate
!!            model
!! Modify histroy:
!!
!! FLOWCHART    : None
!! \n
!_
!================================================================================================================================


  SUBROUTINE thermosoil_cond_pft (kjpindex, njsc, smc, qz, sh2o,zx1,zx2,porosnet,cnd)

    !! 0. Variables and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                                 :: kjpindex      !! Domain size (unitless)
    INTEGER(i_std), DIMENSION (kjpindex), INTENT (in)          :: njsc          !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless)
    REAL(r_std), DIMENSION (kjpindex,ngrnd), INTENT(IN)        :: smc           !! Volumetric Soil Moisture Content (m3/m3)
    REAL(r_std), DIMENSION (nscm), INTENT(IN)                  :: qz            !! Quartz Content (Soil Type Dependent) (0-1)               
    REAL(r_std), DIMENSION (kjpindex,ngrnd,nvm), INTENT(IN)    :: porosnet      !! Soil Porosity (0-1)
    REAL(r_std), DIMENSION (kjpindex,ngrnd), INTENT(IN)        :: sh2o          !! Unfrozen Soil Moisture Content; Frozen Soil Moisture = smc - sh2o
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm), INTENT(IN)     :: zx1, zx2      !! proportion of organic and mineral soil

    !! 0.2 Output variables
    REAL(r_std), DIMENSION (kjpindex,ngrnd,nvm), INTENT(OUT)   :: cnd           !! Soil Thermal Conductivity (W/m/k)

    !! 0.3 Modified variables

    !! 0.4 Local variables
    REAL(r_std), DIMENSION (kjpindex,ngrnd,nvm)                :: ake           !! Kerston Number (unitless)
    REAL(r_std), DIMENSION (kjpindex,ngrnd,nvm)                :: thksat        !! Saturated Thermal Conductivity (W/m/k)
    REAL(r_std), DIMENSION (kjpindex,ngrnd)                    :: thksat_moss   !! Moss saturated thermal conductivity (W/m/K)
    REAL(r_std), DIMENSION (kjpindex,ngrnd,nvm)                :: satratio      !! Degree of Saturation (0-1)
    REAL(r_std), DIMENSION (kjpindex,ngrnd,nvm)                :: xu            !! Unfrozen Volume For Saturation (0-1)
    REAL(r_std), DIMENSION (kjpindex,ngrnd,nvm)                :: xunfroz       !! Unfrozon Volume Fraction (0-1)
    REAL(r_std)                                                :: thko          !! Thermal Conductivity for Other Ssoil Components (W/m/k)
    REAL(r_std), DIMENSION (kjpindex)                          :: gammd         !! Dry Density (kg/m3)
    REAL(r_std), DIMENSION (kjpindex)                          :: thkdry_min    !! Dry Thermal Conductivity for mineral soil (W/m/k)
    REAL(r_std), DIMENSION (kjpindex,ngrnd,nvm)                :: thkdry        !! Dry Thermal Conductivity considering organic carbon (W/m/k)
    REAL(r_std), DIMENSION (kjpindex)                          :: thks_min      !! Thermal Conductivity for the Solids Combined (Quartz + Other) (W/m/k)
    REAL(r_std), DIMENSION (kjpindex,ngrnd,nvm)                :: thks          !! Thermal Conductivity considering organic carbon (W/m/k)
    REAL(r_std), DIMENSION (kjpindex,ngrnd)                    :: thkmoss       !! Thermal conductivity of moss (W/m/K)
    REAL(r_std), DIMENSION (kjpindex,ngrnd)                    :: VWC
    INTEGER(i_std)                                             :: ji, jg, jst, jv
    REAL(r_std), PARAMETER                                     :: THKICE = 2.2 !! Ice Thermal Conductivity (W/m/k)
    REAL(r_std), PARAMETER                                     :: THKQTZ = 7.7 !! Thermal Conductivity for Quartz (W/m/k)
    REAL(r_std), PARAMETER                                     :: THKW = 0.57 !! Water Thermal Conductivity (W/m/k)
    !_================================================================================================================================

    thksat(:,:,:)=0
    thkmoss(:,:)=0
    thksat_moss(:,:)=0
    !! 1. Dry and Saturated Thermal Conductivity.

    DO ji = 1,kjpindex
       jst = njsc(ji)

       !! 1.1. Dry density (Kg/m3) and Dry thermal conductivity (W.M-1.K-1)
       gammd(ji) = (1. - mcs(jst))*2700.
       thkdry_min(ji) = (0.135* gammd(ji) + 64.7)/ (2700. - 0.947* gammd(ji))


       !! 1.2. thermal conductivity of "other" soil components
       IF (qz(jst) > 0.2) THEN
          thko = 2.0
       ELSEIF (qz(jst) <= 0.2) THEN
          thko = 3.0
       ENDIF

       !! 1.3. Thermal conductivity of solids
       thks_min(ji) = (THKQTZ ** qz(jst))* (thko ** (1. - qz(jst)))
    ENDDO

    DO jv = 1,nvm

       SELECTCASE (use_soilc_method)
       CASE (SOILC_METHOD_ARITHMETIC)
          DO jg = 1,ngrnd
             DO ji = 1,kjpindex
                thks(ji,jg,jv) = zx1(ji,jg,jv) * cond_solid_org + zx2(ji,jg,jv) * thks_min(ji)
             ENDDO
          ENDDO
       CASE (SOILC_METHOD_GEOMETRIC)
          DO jg = 1,ngrnd
             DO ji = 1,kjpindex
                ! use geometric mean rather than arithmetic mean (Decharme et al 2016)
                thks(ji,jg,jv) =(cond_solid_org**zx1(ji,jg,jv)) * (thks_min(ji)**zx2(ji,jg,jv))
             ENDDO
          ENDDO
       ENDSELECT

       DO jg = 1,ngrnd
          !! 1.4. saturation ratio
          DO ji = 1,kjpindex
             satratio(ji,jg,jv) = smc(ji,jg) / porosnet(ji, jg, jv)
          ENDDO

          !! 1.5. Saturated Thermal Conductivity (thksat)
          DO ji = 1,kjpindex
             IF ( smc(ji,jg) > min_sechiba ) THEN
                xunfroz(ji,jg,jv) = sh2o(ji,jg) / smc(ji,jg)   ! Unfrozen Fraction (From i.e., 100%Liquid, to 0. (100% Frozen))
                xu(ji,jg,jv) = xunfroz(ji,jg,jv) * porosnet(ji, jg, jv)  ! Unfrozen volume for saturation (porosity*xunfroz)
                thksat(ji,jg,jv) = thks(ji,jg,jv) ** (1. - porosnet(ji, jg, jv)) &
                     * THKICE ** (porosnet(ji, jg, jv) - xu(ji,jg,jv)) * THKW ** xu(ji,jg,jv)
                IF (ok_moss) THEN
                   IF (zlt(jg) .LE. soil_moss_depth) THEN
                      ! Saturated conductivity for moss
                      thksat_moss(ji,jg) = thksat_moss_liq ** xu(ji,jg,jv) * thksat_moss_ice ** (porosnet(ji, jg, jv) - xu(ji,jg,jv))
                   ENDIF
                ENDIF
             ELSE
                ! this value will not be used since ake=0 for this case
                thksat(ji,jg,jv)=0 
             END IF
          ENDDO
       END DO ! DO jg = 1,ngrnd

       !! 2. Kerston Number (ake)
       DO jg = 1,ngrnd
          DO ji = 1,kjpindex
             jst = njsc(ji)
             IF ( (sh2o(ji,jg) + 0.0005) <  smc(ji,jg) ) THEN
                ! Frozen
                ake(ji,jg,jv) = satratio(ji,jg,jv)
             ELSE
                ! Unfrozen
                IF ( satratio(ji,jg,jv) >  0.1 ) THEN
                   IF ( jst < 4 )  THEN
                      ! Coarse
                      ake(ji,jg,jv) = 0.7 * LOG10 (satratio(ji,jg,jv)) + 1.0
                   ELSE
                      ! Fine
                      ake(ji,jg,jv) = LOG10 (satratio(ji,jg,jv)) + 1.0
                   ENDIF
                ELSEIF ( satratio(ji,jg,jv) >  0.05 .AND. satratio(ji,jg,jv) <=  0.1 ) THEN
                   IF ( jst < 4 )  THEN
                      ! Coarse
                      ake(ji,jg,jv) = 0.7 * LOG10 (satratio(ji,jg,jv)) + 1.0
                   ELSE
                      ! Fine
                      ake(ji,jg,jv) = 0.0
                   ENDIF
                ELSE
                   ake(ji,jg,jv) = 0.0  ! use k = kdry
                END IF
             END IF
          ENDDO
       END DO ! DO jg = 1,ngrnd

       SELECTCASE (use_soilc_method)
       CASE (SOILC_METHOD_ARITHMETIC)
          DO jg = 1,ngrnd
             DO ji = 1,kjpindex
                thkdry(ji,jg,jv) = zx1(ji,jg,jv) * cond_dry_org + zx2(ji,jg,jv) * thkdry_min(ji)
             ENDDO
          ENDDO
       CASE (SOILC_METHOD_GEOMETRIC)
          DO jg = 1,ngrnd
             DO ji = 1,kjpindex
                ! use geometric mean rather than arithmetic mean (Decharme et al 2016)
                thkdry(ji,jg,jv) =(cond_dry_org**zx1(ji,jg,jv)) * (thkdry_min(ji)**zx2(ji,jg,jv))
             ENDDO
          ENDDO
       CASE DEFAULT
          CALL ipslerr_p(3,'thermosoil_cond_pft','Unsupported USE_SOILC_METHOD','','')
       ENDSELECT

       !! 3. Thermal conductivity (cnd)
       DO jg = 1,ngrnd
          DO ji = 1,kjpindex
             cnd(ji,jg,jv) = ake(ji,jg,jv) * (thksat(ji,jg,jv) - thkdry(ji, jg, jv)) + thkdry(ji, jg, jv)
             ! Add moss thermal conductivity
             IF (ok_moss) THEN
                IF (zlt(jg) .LE. soil_moss_depth) THEN
                   thkmoss(ji,jg) = thkmoss_dry + &
                        tmc_layt(ji,jg) / mille / dlt(jg) *(thksat_moss(ji,jg) - thkmoss_dry) 
                   
                   SELECTCASE(use_moss_cond_method)
                   CASE (moss_combine)
                      cnd(ji,jg,jv) = moss_frac(jv) * thkmoss(ji,jg) * cnd(ji,jg,jv) *&
                           soil_moss_depth / (moss_layer_thickness * cnd(ji,jg,jv) +&
                           soil_moss_depth * thkmoss(ji,jg)) + (1 - moss_frac(jv)) * cnd(ji,jg,jv) 
                   CASE (moss_replace)
                      cnd(ji,jg,jv) = moss_frac(jv) * thkmoss(ji,jg) + (1 - moss_frac(jv)) * cnd(ji,jg,jv)
                   ENDSELECT
                ENDIF
             ENDIF
          ENDDO
       END DO
       
    END DO !DO jv = 1,nvm 

    IF (ok_moss) CALL xios_orchidee_send_field("thkmoss",thkmoss)


  END SUBROUTINE thermosoil_cond_pft

!! ================================================================================================================================
!! SUBROUTINE   : thermosoil_humlev
!!
!>\BRIEF        Interpolates the diagnostic soil humidity profile shumdiag_perma(nslm, diagnostic axis) onto 
!!              the thermal axis, which gives shum_ngrnd_perma(ngrnd, thermal axis).
!!
!! DESCRIPTION  :  Interpolate the volumetric soil moisture content from the node to the interface of the layer. 
!!                 The values for the deep layers in thermosoil where hydrology is not existing are constant. 
!!                 No interpolation is needed for the total soil moisture content and for the soil saturation degree.
!! The depths of the diagnostic levels are diaglev(1:nslm), computed in slowproc.f90.
!! Recall that when the 11-layer hydrology is used,
!! shum_ngrnd_perma and shumdiag_perma are with reference to the moisture content (mc)
!! at the wilting point mcw : shum_ngrnd_perma=(mc-mcw)/(mcs-mcw).
!! with mcs the saturated soil moisture content.
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): mc_layt, mcl_layt, tmc_layt, shum_ngrnd_perma (soil humidity profile on the thermal axis)
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None 
!! \n 
!_ ================================================================================================================================
  SUBROUTINE thermosoil_humlev(kjpindex, shumdiag_perma, mc_layh, mcl_layh, tmc_layh, &
       mc_layh_pft, mcl_layh_pft, tmc_layh_pft)

    !! 0. Variables and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                            :: kjpindex    !! Domain size (unitless)
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (in)    :: shumdiag_perma    !! Relative soil humidity on the diagnostic axis. 
    !! (0-1, unitless). Caveats : when "hydrol" (the 11-layers
    !! hydrology) is used, this humidity is calculated with 
    !! respect to the wilting point : 
    !! shumdiag_perma= (mc-mcw)/(mcs-mcw), with mc : moisture 
    !! content; mcs : saturated soil moisture content; mcw: 
    !! soil moisture content at the wilting point. when the 2-layers
    !! hydrology "hydrolc" is used, shumdiag_perma is just
    !! a diagnostic humidity index, with no real physical 
    !! meaning.
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (in)    :: mc_layh     !! Volumetric soil moisture content for each layer in hydrol at nodes(liquid+ice) [m/s]
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (in)    :: mcl_layh    !! Volumetric soil moisture content for each layer in hydrol at nodes(liquid) [m/s]
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (in)    :: tmc_layh    !! Total soil moisture content for each layer in hydrol(liquid+ice) [mm]
    REAL(r_std),DIMENSION (kjpindex,nslm,nvm), INTENT (in) :: mc_layh_pft     !! Volumetric soil moisture content for each layer in hydrol at nodes(liquid+ice) [m/s]
    REAL(r_std),DIMENSION (kjpindex,nslm,nvm), INTENT (in) :: mcl_layh_pft    !! Volumetric soil moisture content for each layer in hydrol at nodes(liquid) [m/s]
    REAL(r_std),DIMENSION (kjpindex,nslm,nvm), INTENT (in) :: tmc_layh_pft    !! Total soil moisture content for each layer in hydrol(liquid+ice) [mm]

    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                                        :: ji, jd, jv

    !_ ================================================================================================================================

    IF (printlev >= 4) WRITE(numout,*) 'Start thermosoil_humlev'

    ! The values for the deep layers in thermosoil where hydrology is not existing are constant. 
    ! For exemple if thermosoil uses 8m, and hydrol uses 2m vertical discretization,
    ! the values between 2m and 8m are constant.
    ! The moisture computed in hydrol is at the nodes (except for the
    ! top and bottom layer which are at interfaces)
    ! A linear interpolation is applied to obtain the moisture values at
    ! the interfaces (mc_layt), from the mc_layh at the nodes

    DO ji=1,kjpindex
       DO jd = 1, nslm
          IF(jd == 1) THEN ! the moisture at the 1st interface mc_layh(1) is at the surface, no interpolation
             mc_layt(ji,jd) = mc_layh(ji,jd)
             mcl_layt(ji,jd) = mcl_layh(ji,jd)
          ELSEIF(jd == 2) THEN  !! the mc_layt at the 2nd interface is interpolated using mc_layh(1) at surface and mc_layh(2) at the node
             mc_layt(ji, jd) = mc_layh(ji,jd-1)*(znt(jd)-zlt(jd-1))/(znt(jd)-0.0) + &
                  mc_layh(ji, jd)*(zlt(jd-1)-0.0)/(znt(jd)-0.0)
             mcl_layt(ji, jd) = mcl_layh(ji,jd-1)*(znt(jd)-zlt(jd-1))/(znt(jd)-0.0) + &
                  mcl_layh(ji, jd)*(zlt(jd-1)-0.0)/(znt(jd)-0.0)
          ELSEIF(jd == nslm) THEN ! the mc_layt at the nslm interface is interpolated using mc_layh(nslm) and mc_layh(nslm-1)
             mc_layt(ji, jd) = mc_layh(ji,jd-1)*(zlt(jd)-zlt(jd-1))/(zlt(jd)-znt(jd-1))  + &
                  mc_layh(ji,jd)*(zlt(jd-1)-znt(jd-1))/(zlt(jd)-znt(jd-1))
             mcl_layt(ji, jd) = mcl_layh(ji,jd-1)*(zlt(jd)-zlt(jd-1))/(zlt(jd)-znt(jd-1))  + &
                  mcl_layh(ji,jd)*(zlt(jd-1)-znt(jd-1))/(zlt(jd)-znt(jd-1))
          ELSE ! the mc_layt at the other interfaces are interpolated using mc_layh at adjacent nodes.
             mc_layt(ji, jd) = mc_layh(ji, jd-1)*(1-dz5(jd-1)) + mc_layh(ji,jd)*dz5(jd-1)
             mcl_layt(ji, jd) = mcl_layh(ji, jd-1)*(1-dz5(jd-1)) + mcl_layh(ji,jd)*dz5(jd-1)
          ENDIF
          DO jv = 1,nvm
             IF(jd == 1) THEN
                mc_layt_pft(ji,jd,jv) = MAX(mc_layh_pft(ji,jd,jv), min_sechiba)
                mcl_layt_pft(ji,jd,jv) = MAX(mcl_layh_pft(ji,jd,jv), min_sechiba)
             ELSEIF(jd == 2) THEN
                mc_layt_pft(ji,jd,jv) = MAX(mc_layh_pft(ji,jd-1,jv)*(znt(jd)-zlt(jd-1))/(znt(jd)-0.0) + &
                     mc_layh_pft(ji,jd,jv)*(zlt(jd-1)-0.0)/(znt(jd)-0.0), min_sechiba)
                mcl_layt_pft(ji,jd,jv) = MAX(mcl_layh_pft(ji,jd-1,jv)*(znt(jd)-zlt(jd-1))/(znt(jd)-0.0) + &
                     mcl_layh_pft(ji,jd,jv)*(zlt(jd-1)-0.0)/(znt(jd)-0.0), min_sechiba)
             ELSEIF(jd == nslm) THEN
                mc_layt_pft(ji,jd,jv) = MAX(mc_layh_pft(ji,jd-1,jv)*(zlt(jd)-zlt(jd-1))/(zlt(jd)-znt(jd-1))  + &
                     mc_layh_pft(ji,jd,jv)*(zlt(jd-1)-znt(jd-1))/(zlt(jd)-znt(jd-1)),min_sechiba)
                mcl_layt_pft(ji,jd,jv) = MAX(mcl_layh_pft(ji,jd-1,jv)*(zlt(jd)-zlt(jd-1))/(zlt(jd)-znt(jd-1))  + &
                     mcl_layh_pft(ji,jd,jv)*(zlt(jd-1)-znt(jd-1))/(zlt(jd)-znt(jd-1)),min_sechiba)
             ELSE
                mc_layt_pft(ji,jd,jv) = MAX(mc_layh_pft(ji,jd-1,jv)*(1-dz5(jd-1)) + &
                     mc_layh_pft(ji,jd,jv)*dz5(jd-1), min_sechiba)
                mcl_layt_pft(ji,jd,jv) = MAX(mcl_layh_pft(ji,jd-1,jv)*(1-dz5(jd-1)) + &
                     mcl_layh_pft(ji,jd,jv)*dz5(jd-1), min_sechiba)
             ENDIF
          ENDDO ! jv
          tmc_layt(ji,jd) = tmc_layh(ji,jd)
          tmc_layt_pft(ji,jd,:) = tmc_layh_pft(ji,jd,:)
       ENDDO !jd

       ! The deep layers in thermosoil where hydro is not existing
       DO jd = nslm+1, ngrnd
          mc_layt(ji,jd) = mc_layh(ji,nslm)
          mcl_layt(ji,jd) = mcl_layh(ji,nslm)
          tmc_layt(ji,jd) = tmc_layh(ji,nslm)/dlt(nslm) *dlt(jd)
          mc_layt_pft(ji,jd,:) = mc_layh_pft(ji,nslm,:)
          mcl_layt_pft(ji,jd,:) = mcl_layh_pft(ji,nslm,:)
          tmc_layt_pft(ji,jd,:) = tmc_layh_pft(ji,nslm,:)/dlt(nslm) *dlt(jd)
       ENDDO
    ENDDO


    ! The values for the deep layers in thermosoil where hydro is not existing are constant.
    ! For exemple if thermosoil uses 8m, and hydrol uses 2m vertical discretization,
    ! the values between 2m and 8m are constant.

    DO jd = 1, nslm
       shum_ngrnd_perma(:,jd) = shumdiag_perma(:,jd)
    END DO
    DO jd = nslm+1,ngrnd
       shum_ngrnd_perma(:,jd) = shumdiag_perma(:,nslm)
    END DO

    IF (printlev >= 4) WRITE(numout,*) 'thermosoil_humlev done'

  END SUBROUTINE thermosoil_humlev



!! ================================================================================================================================
!! SUBROUTINE   : thermosoil_energy_diag
!!
!>\BRIEF         Calculate diagnostics 
!!
!! DESCRIPTION  : Calculate diagnostic variables coldcont_incr and coldcont_incr
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None 
!! \n 
!_ ================================================================================================================================

  SUBROUTINE thermosoil_energy_diag(kjpindex, temp_sol_new, soilcap)

    !! 0. Variables and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                     :: kjpindex     !! Domain size (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)  :: temp_sol_new !! Surface temperature at the present time-step, Ts 
    !! @tex ($K$) @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)  :: soilcap      !! Apparent surface heat capacity 
    !! @tex ($J m^{-2} K^{-1}$) @endtex, 
    !! see eq. A29 of F. Hourdin\'s PhD thesis.

    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                                 :: ji, jg
    !_ ================================================================================================================================

    !  Sum up the energy content of all layers in the soil.
    DO ji = 1, kjpindex

       IF (pcapa_en(ji,1) .LE. sn_capa) THEN

          ! Verify the energy conservation in the surface layer
          coldcont_incr(ji) = soilcap(ji) * (temp_sol_new(ji) - temp_sol_beg(ji))
          surfheat_incr(ji) = zero
       ELSE

          ! Verify the energy conservation in the surface layer
          surfheat_incr(ji) = soilcap(ji) * (temp_sol_new(ji) - temp_sol_beg(ji))
          coldcont_incr(ji) = zero
       ENDIF
    ENDDO

    ! Save temp_sol_new to be used at next timestep
    temp_sol_beg(:)   = temp_sol_new(:)

  END SUBROUTINE thermosoil_energy_diag



!! ================================================================================================================================
!! SUBROUTINE   : thermosoil_readjust
!!
!>\BRIEF        
!!
!! DESCRIPTION	: Energy conservation : Correction to make sure that the same latent heat is released and 
!!                consumed during freezing and thawing  
!!
!! RECENT CHANGE(S) : None
!! 
!! MAIN OUTPUT VARIABLE(S): ptn (soil temperature profile on the thermal axis per pft) and ptn_pftmean  
!!                          
!! REFERENCE(S) :
!!
!! FLOWCHART    : None 
!! \n 
!_ ================================================================================================================================

  SUBROUTINE thermosoil_readjust(kjpindex, ptn, ptn_pftmean, veget_max)

    !! 0. Variables and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                             :: kjpindex
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(in)      :: veget_max !! Fraction of vegetation type

    !! 0.2 Modified variables
    REAL(r_std),DIMENSION(kjpindex,ngrnd,nvm),INTENT(inout)    :: ptn
    REAL(r_std),DIMENSION(kjpindex,ngrnd),INTENT(inout)        :: ptn_pftmean

    !! 0.3 Local variables
    INTEGER(i_std)  :: ji, jg, m, jv
    INTEGER(i_std)  :: lev3m  !! Closest interface level to 3m
    REAL(r_std) :: ptn_tmp
    REAL(r_std), DIMENSION (kjpindex,ngrnd) :: ptn_beg_pftmean          

    ! The energy is spread over the layers down to approximatly 3m
    ! Find the closest level to 3m. It can be below or above 3m.
    lev3m=MINLOC(ABS(zlt(:)-3.0),dim=1)
    IF (printlev >= 3) WRITE(numout,*) 'In thermosoil_adjust: lev3m=',lev3m, ' zlt(lev3m)=', zlt(lev3m)

    ! Calculate ptn_beg for mean PFT
    IF (ok_soil_carbon_discretization) THEN
       ptn_beg_pftmean(:,:) = zero
       DO m=1,nvm
          DO jg = 1, ngrnd
             ptn_beg_pftmean(:,jg) = ptn_beg_pftmean(:,jg) + ptn_beg(:,jg,m) * veget_max(:,m)
          END DO
       END DO
    ELSE
       ! ptn is not depending on pft
       ptn_beg_pftmean(:,:) = ptn_beg(:,:,1)
    END IF

    DO jg=1, ngrnd
       DO ji=1, kjpindex
          ! All soil latent energy is put into e_soil_lat(ji)
          ! because the variable soil layers make it difficult to keep track of all
          ! layers in this version
          ! NOTE: pcapa has unit J/K/m3 and pcappa_supp has J/K/m2
          ! NOTE: Right now this is not effective as ptn_beg_pftmean equal ptn_pftmean
          e_soil_lat(ji)=e_soil_lat(ji)+pcappa_supp(ji,jg)*(ptn_pftmean(ji,jg)-ptn_beg_pftmean(ji,jg))
       END DO
    END DO

    DO ji=1, kjpindex
       IF (e_soil_lat(ji).GT.min_sechiba.AND.MINVAL(ptn_pftmean(ji,:)).GT.fr_center+fr_dT/2.) THEN
          ! The soil is thawed: we spread the excess of energy over the uppermost lev3m levels
          ! Here we increase the temperatures
          DO jg=1, lev3m
             ptn_tmp=ptn_pftmean(ji,jg)
             ptn_pftmean(ji,jg)=ptn_pftmean(ji,jg)+MIN(e_soil_lat(ji)/pcapa(ji,jg)/zlt(lev3m), 0.5)
             e_soil_lat(ji)=e_soil_lat(ji)-(ptn_pftmean(ji,jg)-ptn_tmp)*pcapa(ji,jg)*dlt(jg)
          ENDDO
       ELSE IF (e_soil_lat(ji).LT.-min_sechiba.AND.MINVAL(ptn_pftmean(ji,:)).GT.fr_center+fr_dT/2.) THEN
          ! The soil is thawed
          ! Here we decrease the temperatures
          DO jg=1, lev3m
             ptn_tmp=ptn_pftmean(ji,jg)
             ptn_pftmean(ji,jg)=MAX(fr_center+fr_dT/2., ptn_tmp+e_soil_lat(ji)/pcapa(ji,jg)/zlt(lev3m))
             e_soil_lat(ji)=e_soil_lat(ji)+(ptn_tmp-ptn_pftmean(ji,jg))*pcapa(ji,jg)*dlt(jg)
          END DO
       END IF
    END DO


    !! Calculate ptn per pft
    !  Here ptn is the same at all pfts when ok_soil_carbon_discretization is activated or not. 
    !  ptn is modified in thermosoil_add_heat_zimov if ok_zimov=TRUE.
    DO jv=1,nvm
       ptn(:,:,jv) = ptn_pftmean(:,:)
    END DO


  END SUBROUTINE thermosoil_readjust
   
!-------------------------------------------------------------------



!! ================================================================================================================================
!! SUBROUTINE   : thermosoil_getdiff
!!
!>\BRIEF          Computes soil and snow heat capacity and conductivity    
!!
!! DESCRIPTION	: Computation of the soil thermal properties; snow properties are also accounted for
!!
!! RECENT CHANGE(S) : None
!! 
!! MAIN OUTPUT VARIABLE(S): 
!!                          
!! REFERENCE(S) :
!!
!! FLOWCHART    : None 
!! \n 
!_ ================================================================================================================================

  SUBROUTINE thermosoil_getdiff( kjpindex, snow, ptn, njsc, snowrho, snowtemp, pb, totfrac_nobio, &
       icetemp, ice_sheet_mask )

    !! 0. Variables and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std),INTENT(in)				:: kjpindex
    REAL(r_std),DIMENSION(kjpindex),INTENT (in)	        :: snow       !! Snow mass
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)    :: njsc       !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless)
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(in) :: snowrho    !! Snow density (Kg/m^3)
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(in) :: snowtemp   !! Snow temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: pb         !! Surface presure (hPa)
    REAL(r_std),DIMENSION(kjpindex,ngrnd),INTENT(in)	:: ptn        !! Soil temperature profile
    REAL(r_std),DIMENSION (kjpindex),INTENT(in)         :: totfrac_nobio !! Total fraction of continental ice+lakes+cities+...
    REAL(r_std), DIMENSION (kjpindex,nice), INTENT(in)  :: icetemp    !! Ice temperature (K)
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)    :: ice_sheet_mask !! Ice sheet mask

    !! 0.3 Local variables
    REAL						:: xx         !! Unfrozen fraction of the soil
    REAL(r_std), DIMENSION(kjpindex)             	:: snow_h
    REAL(r_std), DIMENSION(kjpindex,ngrnd) 		:: zx1, zx2    
    INTEGER						:: ji,jg
    INTEGER                                             :: jst
    REAL(r_std), DIMENSION (kjpindex,ngrnd)             :: pcapa_tmp  !! soil heat capacity (J/m3/K)
    REAL(r_std), DIMENSION (kjpindex,ngrnd)             :: pcapa_spec !! SPECIFIC soil heat capacity (J/kg/K)
    REAL(r_std)                                         :: rho_tot    !! Soil density (kg/m3)
    REAL(r_std), DIMENSION (kjpindex,ngrnd)             :: profil_froz_old

    pcapa_tmp(:,:) = 0.0

    !! Computes soil heat capacity and conductivity
    DO ji = 1,kjpindex

       ! Since explicitsnow module is implemented set zx1=0 and zx2=1
       zx1(ji,:) = 0.
       zx2(ji,:) = 1.

       DO jg = 1, ngrnd
          jst = njsc(ji)
          pcapa_tmp(ji, jg) = so_capa_dry(jst) * (1-mcs(jst)) + water_capa * tmc_layt(ji,jg)/mille/dlt(jg)
          !
          ! 2. Calculate volumetric heat capacity with allowance for permafrost
          ! 2.1. soil heat capacity depending on temperature and humidity
          ! For SP6MIP we also diagnose a specific heat capacity (pcapa_spec), 
          ! which requires to calculate the total density of the soil (rho_tot), always >> 0

          IF (ptn(ji,jg) .LT. fr_center-fr_dT/2.) THEN
             ! frozen soil
             profil_froz(ji,jg) = 1.
             pcappa_supp(ji,jg)= 0.
             pcapa(ji, jg) = so_capa_dry(jst) * (1-mcs(jst)) + so_capa_ice(jst) * tmc_layt(ji,jg) / mille / dlt(jg)
             rho_tot = rho_soil * (1-mcs(jst)) + rho_ice * tmc_layt(ji,jg) / mille / dlt(jg) 
             pcapa_spec(ji, jg) = pcapa(ji, jg) / rho_tot

          ELSEIF (ptn(ji,jg) .GT. fr_center+fr_dT/2.) THEN
             ! unfrozen soil	 
             pcapa(ji, jg) = pcapa_tmp(ji, jg)
             profil_froz(ji,jg) = 0.
             pcappa_supp(ji,jg)= 0.
             rho_tot = rho_soil * (1-mcs(jst)) + rho_water * tmc_layt(ji,jg)/mille/dlt(jg)
             pcapa_spec(ji, jg) = pcapa(ji, jg) / rho_tot
          ELSE
             ! xx is the unfrozen fraction of soil water	   	   
             xx = (ptn(ji,jg)-(fr_center-fr_dT/2.)) / fr_dT
             profil_froz(ji,jg) = (1. - xx)

             IF (ok_freeze_thaw_latent_heat) THEN
                pcapa(ji, jg) = so_capa_dry(jst) * (1-mcs(jst)) + &
                     water_capa * tmc_layt(ji,jg)/mille / dlt(jg) * xx + &
                     so_capa_ice(jst) * tmc_layt(ji,jg) / mille/dlt(jg) * (1.-xx) + &
                     shum_ngrnd_perma(ji,jg)*mcs(jst)*lhf*rho_water/fr_dT
             ELSE
                pcapa(ji, jg) = so_capa_dry(jst) * (1-mcs(jst)) + &
                     water_capa * tmc_layt(ji,jg)/mille / dlt(jg) * xx + &
                     so_capa_ice(jst) * tmc_layt(ji,jg) / mille/dlt(jg) * (1.-xx)
             ENDIF

             rho_tot =  rho_soil* (1-mcs(jst)) + &
                  rho_water * tmc_layt(ji,jg)/mille / dlt(jg) * xx + &
                  rho_ice * tmc_layt(ji,jg) / mille/dlt(jg) * (1.-xx)
             pcapa_spec(ji, jg) = pcapa(ji, jg) / rho_tot

             ! NOTE: pcappa_supp is transformed from volumetric heat capacity into surfacic heat capacity (J/K/m2)
             pcappa_supp(ji,jg)= shum_ngrnd_perma(ji,jg)*mcs(jst)*lhf*rho_water/fr_dT*zx2(ji,jg)*dlt(jg)

          ENDIF
          !BG:is this 2.2 and 2.3 necessary since zx1=0?
          ! 
          ! 2.2. Take into account the snow and soil fractions in the layer
          !
          pcapa(ji,jg) = zx1(ji,jg) * sn_capa + zx2(ji,jg) * pcapa(ji,jg)

          ! 2.3. Calculate the heat capacity for energy conservation check 
          IF ( zx1(ji,jg).GT.0. ) THEN
             pcapa_en(ji,jg) = sn_capa
          ELSE
             pcapa_en(ji,jg) = pcapa(ji,jg)
          ENDIF
       END DO
    ENDDO

    ! Output the specific heat capcaity for SP-MIP
    CALL xios_orchidee_send_field("pcapa_spec",pcapa_spec)

    !
    ! 3. Calculate the heat conductivity with allowance for permafrost
    !
    IF (ok_freeze_thaw_latent_heat) THEN
       CALL thermosoil_cond (kjpindex, njsc, mc_layt, QZ, mcl_layt*(1-profil_froz), pkappa)
    ELSE
       CALL thermosoil_cond (kjpindex, njsc, mc_layt, QZ, mcl_layt, pkappa)
    ENDIF

    !! Computes snow heat capacity and conductivity   
    DO ji = 1,kjpindex
       pcapa_snow(ji,:) = snowrho(ji,:) * xci

       SELECTCASE (snow_cond_method)
       CASE (SNOW_COND_METHOD_DEFAULT)
          pkappa_snow(ji,:) = (ZSNOWTHRMCOND1 + ZSNOWTHRMCOND2*snowrho(ji,:)*snowrho(ji,:)) +      &
               MAX(0.0,(ZSNOWTHRMCOND_AVAP+(ZSNOWTHRMCOND_BVAP/(snowtemp(ji,:)+ &
               ZSNOWTHRMCOND_CVAP)))*(XP00/(pb(ji)*100.)))
       CASE (SNOW_COND_METHOD_DECHARME16)
          pkappa_snow(ji,:) = 2.2*((snowrho(ji,:)/1000.)**2.0) +      &
               MAX(0.0,(ZSNOWTHRMCOND_AVAP+(ZSNOWTHRMCOND_BVAP/(snowtemp(ji,:)+ &
               ZSNOWTHRMCOND_CVAP)))*(XP00/(pb(ji)*100.)))
       CASE (SNOW_COND_METHOD_SNOWICE)
          ! 1. thermal conductivity from Yen 1981 :
          pkappa_snow(ji,:) =  cond_ice * exp(1.88*log(snowrho(ji,:)/rho_water))
          ! 2. Implicit vapor diffn effects Sun et al 1999 :     
          pkappa_snow(ji,:) = pkappa_snow(ji,:) +      &
               MAX(0.0,(ZSNOWTHRMCOND_AVAP+(ZSNOWTHRMCOND_BVAP/(snowtemp(ji,:)+ &
               ZSNOWTHRMCOND_CVAP)))*(XP00/(pb(ji)*100.)))
       CASE DEFAULT
          CALL ipslerr_p(3,'thermosoil_getdiff_pft','Unsupported SNOW_COND_METHOD', &
               'Currently supported methods are ','default(1), ducharme16(2) or snowice(3)')
       ENDSELECT

       IF ( ok_ice_sheet .and. (ice_sheet_mask(ji).GT.0) ) THEN
          pkappa_ice(ji,:) = ZICETHRMCOND1 * exp(ZICETHRMCOND2*icetemp(ji,:))  !  from Yin CHao Yen (Temp > 195K)
          pcapa_ice(ji,:) = rho_ice * (ZICETHRMHEAT1 + ZICETHRMHEAT2*(icetemp(ji,:)-tp_00)) ! formulation as in GRISLI prop_th_icetemp.f90
       ELSE
          pkappa_ice(ji,:) =  ZICETHRMCOND1 * exp(ZICETHRMCOND2*(tp_00))
          pcapa_ice(ji,:) = rho_ice * ZICETHRMHEAT1 ! formulation as in GRISLI prop_th_icetemp.f90
       ENDIF
    END DO

  END SUBROUTINE thermosoil_getdiff



!! ================================================================================================================================
!! SUBROUTINE   : thermosoil_getdiff_pft
!!
!>\BRIEF          Computes soil and snow heat capacity and conductivity    
!!
!! DESCRIPTION  : Computation of the soil thermal properties; snow properties are also accounted for
!!
!! RECENT CHANGE(S) : None
!! 
!! MAIN OUTPUT VARIABLE(S): Following module variables are calculted in this subroutine: 
!!                          profil_froz_pft, pcapa_per_pft, pkappa_per_pft, pcapa_en_per_pft, tmc_layt_pft
!!                          
!! REFERENCE(S) :
!!
!! FLOWCHART    : None 
!! \n 
!_ ================================================================================================================================

  SUBROUTINE thermosoil_getdiff_pft( kjpindex,  ptn,     njsc,     shum_ngrnd_perma, &
       depth_organic_soil, som_total, snowrho, snowtemp, pb, veget_max, &
       icetemp, ice_sheet_mask)

    !! 0. Variables and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std),INTENT(in)                                          :: kjpindex
    REAL(r_std),DIMENSION(kjpindex,ngrnd),INTENT(in)                   :: shum_ngrnd_perma     !! Water saturation degree on the soil depth axes define by the thermic (0-1, dimensionless)
    REAL(r_std), DIMENSION(kjpindex), INTENT (in)                      :: depth_organic_soil   !! Depth at which there is still organic matter (m)
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm,nelements), INTENT (in)  :: som_total            !! total soil carbon for use in thermal calcs (g/m**3)
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)                   :: njsc                 !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless)
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(in)                :: snowrho              !! Snow density (Kg/m^3)
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(in)                :: snowtemp             !! Snow temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)                      :: pb                   !! Surface presure (hPa)
    REAL(r_std),DIMENSION(kjpindex,ngrnd,nvm),INTENT(in)               :: ptn                  !! Soil temperature profile 
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(in)                   :: veget_max            !! Fraction of vegetation type
    REAL(r_std), DIMENSION (kjpindex,nice), INTENT(in)                 :: icetemp              !! Ice temperature (K)
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)                   :: ice_sheet_mask       !! Ice sheet mask

    !! 0.2 Modified variables

    !! 0.3 Output variables

    !! 0.3 Local variables
    REAL(r_std)                                              :: xx                     !! Unfrozen fraction of the soil
    REAL(r_std)                                              :: p
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm)               :: so_capa_dry_net
    INTEGER(i_std)                                           :: ji,jg,jv
    INTEGER(i_std)                                           :: jst
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm)               :: poros_net
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm)               :: zx1, zx2 !! BE CAREFUL here zx1 and zx2 split between the organic and mineral soils not
    !! not snow vs. liquid water
    REAL(r_std), DIMENSION(kjpindex,ngrnd)                   :: profil_froz_mean, tmp
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm)               :: pcappa_supp_pft
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm)               :: pcapa_pft_tmp  !! soil heat capacity (J/m3/K)
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm)               :: pcapa_spec_pft !! SPECIFIC soil heat capacity (J/kg/K)
    REAL(r_std)                                              :: rho_tot    !! Soil density (kg/m3)
    REAL(r_std), DIMENSION(kjpindex,ngrnd)                   :: pcapa_spec !! SPECIFIC soil heat capacity (J/kg/K)
    REAL(r_std), DIMENSION(kjpindex,ngrnd)                   :: poros_net_glob !! Porosity based on SOC content over all PFTs
    REAL(r_std), DIMENSION(kjpindex,ngrnd)                   :: pcapa_moss !! moss heat capacity (J/m3/K)


    ! Organic and anorgaic layer fraction
    REAL(r_std), DIMENSION(kjpindex,ngrnd)		             :: som_total_pftmean !! Deep soil organic matter profiles calculated for thermosoil, mean over all PFTs in concentration g/m3
    !
    ! Default: organic layer not taken into account
    zx1(:,:,:) = 0.0
    zx2(:,:,:) = 0.0
    poros_net(:,:,:) = 0.0
    poros_net_glob(:,:) = 0.0
    pcapa_spec_pft(:,:,:) = 0.0 
    pcapa_pft_tmp(:,:,:) = 0.0
    pcapa_spec(:,:) = 0.0
    pcapa_moss(:,:) = 0.0
    som_total_pftmean(:,:) = 0.0
    
    IF ( use_soilc_insulation ) THEN
       IF (use_refsoc) THEN
          DO jv = 1,nvm
             DO jg = 1, ngrnd
                DO ji = 1,kjpindex
                   zx1(ji,jg,jv) = refsoc(ji,jg)/soilc_max   !after lawrence and slater
                ENDDO
                som_total_pftmean(:,jg) = som_total_pftmean(:,jg) + refsoc(:,jg) * veget_max(:,jv)
             ENDDO
          ENDDO
       ELSEIF (ok_soil_carbon_discretization) THEN
          DO jv = 1,nvm
             DO jg = 1, ngrnd
                DO ji = 1,kjpindex
                   zx1(ji,jg,jv) = som_total(ji,jg,jv,icarbon)/soilc_max   !after lawrence and slater
                ENDDO
                som_total_pftmean(:,jg) = som_total_pftmean(:,jg) + som_total(:,jg,jv,icarbon) * veget_max(:,jv)
             ENDDO
          ENDDO
       ELSE
          WRITE (numout,*) ' use_soilc_insulation can t be used without use_refsoc or without ok_soil_carbon_discretization. We stop.'
          STOP 'stomate_soil_carbon_discretization_deep_somcycle'
       ENDIF
       WHERE (zx1 > 1)  zx1 = 1
       ! 
    ELSE
       !
       ! level 1
       !
       DO jv = 1,nvm
          DO ji = 1,kjpindex
             IF ( depth_organic_soil(ji) .GT. zlt(1) ) THEN
                !! the 1st level is in the organic => the 1st layer is entirely organic
                zx1(ji,1,jv) = 1. !!zx1 being the fraction of each level that is organic, zx2 is the remainder
             ELSE IF ( depth_organic_soil(ji) .GT. zero ) THEN
                !! the 1st level is beyond the organic and the organic is present
                zx1(ji,1,jv) = depth_organic_soil(ji) / zlt(1)
             ELSE
                ! there is no organic at all
                zx1(ji,1,jv) = 0.
             ENDIF
          ENDDO
       ENDDO
       !
       ! other levels
       !
       DO jg = 2, ngrnd !- 2
          DO ji = 1,kjpindex
             IF ( depth_organic_soil(ji) .GT. zlt(jg) ) THEN
                ! the current level is in the organic => the current layer is
                ! entirely organic
                zx1(ji,jg,1) = 1.
             ELSE IF ( depth_organic_soil(ji) .GT. zlt(jg-1) ) THEN
                ! the current layer is partially organic
                zx1(ji,jg,1) = (depth_organic_soil(ji) - zlt(jg-1)) / (zlt(jg) - zlt(jg-1))
             ELSE
                ! both levels are out of organic => the current layer is entirely
                ! mineral soil       
                zx1(ji,jg,1) = 0.
             ENDIF
          ENDDO
       ENDDO
       DO jv = 2, nvm
          zx1(:,:,jv) = zx1(:,:,1)
       ENDDO
    
       DO jg = 1, ngrnd
          som_total_pftmean(:,jg) = zx1(:,jg,1) * soilc_max
       ENDDO
    ENDIF ! ( use_soilc_insulation ) 
    
    zx2(:,:,:) = 1.-zx1(:,:,:)
     
    CALL xios_orchidee_send_field ("deepC_thermo_pftmean", som_total_pftmean)

#ifdef STRICT_CHECK
    IF (ANY(tmc_layt_pft < min_sechiba)) CALL ipslerr_p(3, "thermosoil_getdiff_pft", "tmc_layt_pft has negative values", "", "") ! prec issues
#endif

    DO jv = 1,nvm
       DO jg = 1, ngrnd
          DO ji = 1,kjpindex
             jst = njsc(ji)
             !
             ! 1. Calculate dry heat capacity and conductivity, taking
             ! into account the organic and mineral fractions in the layer
             !
             ! Former MICT version
             !Here we take into account the new dependance of the soil heat capacity from the soil type.
             so_capa_dry_net(ji,jg,jv) = zx1(ji,jg,jv) * SO_CAPA_DRY_ORG + zx2(ji,jg,jv) * so_capa_dry(jst)

             !Here we take into account the new dependance of the porosity from the soil type.
             poros_net(ji,jg,jv) = zx1(ji,jg,jv) * poros_org + zx2(ji,jg,jv) * mcs(jst)
             !
             ! 2. Calculate heat capacity with allowance for permafrost
          ENDDO
       ENDDO
    ENDDO

    ! Output the global porosity (sum over every PFT)
    DO jv = 1,nvm
       DO jg = 1, ngrnd
          poros_net_glob(:,jg)=poros_net_glob(:,jg) + (poros_net(:,jg,jv)*veget_max(:,jv)/(SUM(veget_max(:,:),dim=2)))
       ENDDO
    END DO
    
    CALL xios_orchidee_send_field("poros_net_thermosoil",poros_net_glob)

    DO jv = 1,nvm
       DO jg = 1, ngrnd
          DO ji = 1,kjpindex
             jst = njsc(ji)
             pcapa_pft_tmp(ji, jg,jv) = so_capa_dry_net(ji,jg,jv) * (1-mcs(jst)) + water_capa * tmc_layt(ji,jg)/mille/dlt(jg)
             ! 2.1. soil heat capacity depending on temperature and humidity
             IF (ptn(ji,jg,jv) .LT. fr_center-fr_dT/2.) THEN
                ! frozen soil
                profil_froz_pft(ji,jg,jv) = 1.
                pcappa_supp_pft(ji,jg, jv)= 0.
                pcapa_per_pft(ji,jg,jv) = so_capa_dry_net(ji,jg,jv)* (1-mcs(jst)) + &
                     so_capa_ice(jst)  * tmc_layt(ji,jg) / mille / dlt(jg) 
                ! Add moss thermal capacity
                IF (ok_moss) THEN
                   IF (zlt(jg) .LE. soil_moss_depth) THEN
                      pcapa_moss(ji,jg) = pcapa_moss_dry + tmc_layt(ji,jg) / mille / dlt(jg) * (pcapa_moss_ice - pcapa_moss_dry)
                      pcapa_per_pft(ji,jg,jv) = pcapa_per_pft(ji,jg,jv) + &
                           moss_frac(jv) * pcapa_moss(ji,jg) * moss_layer_thickness / soil_moss_depth
                   ENDIF
                ENDIF
                rho_tot = rho_soil * (1-poros_net(ji,jg,jv)) + rho_ice * tmc_layt(ji,jg) / mille / dlt(jg)
                pcapa_spec_pft(ji, jg, jv) = pcapa_per_pft(ji, jg, jv) / rho_tot
             ELSEIF (ptn(ji,jg,jv) .GT. fr_center+fr_dT/2.) THEN
                ! unfrozen soil         
                profil_froz_pft(ji,jg,jv) = 0.
                pcappa_supp_pft(ji,jg,jv)= 0.
                pcapa_per_pft(ji,jg,jv) = pcapa_pft_tmp(ji, jg,jv)
                ! Add moss thermal capacity
                IF (ok_moss) THEN
                   IF (zlt(jg) .LE. soil_moss_depth) THEN
                      pcapa_moss(ji,jg) = pcapa_moss_dry + tmc_layt(ji,jg) / mille / dlt(jg) * (pcapa_moss_liq - pcapa_moss_dry)
                      pcapa_per_pft(ji,jg,jv) = pcapa_per_pft(ji,jg,jv) + &
                           moss_frac(jv) * pcapa_moss(ji,jg) * moss_layer_thickness / soil_moss_depth
                   ENDIF
                ENDIF
                rho_tot = rho_soil * (1-poros_net(ji,jg,jv)) + rho_water * tmc_layt(ji,jg)/mille/dlt(jg)
                pcapa_spec_pft(ji, jg, jv) = pcapa_per_pft(ji, jg, jv) / rho_tot
             ELSE

                pcappa_supp_pft(ji,jg,jv)= shum_ngrnd_perma(ji,jg)*mcs(jst)*lhf*rho_water/fr_dT
                IF (jg .GT. nslm) pcappa_supp_pft(ji,jg,jv)= 0.

                ! x is the unfrozen fraction of soil water              
                xx = (ptn(ji,jg,jv)-(fr_center-fr_dT/2.)) / fr_dT
                profil_froz_pft(ji,jg,jv) = (1. - xx)

                IF (ok_freeze_thaw_latent_heat) THEN
                   pcapa_per_pft(ji,jg,jv) = so_capa_dry_net(ji,jg,jv)*(1-mcs(jst)) + &
                        water_capa * tmc_layt(ji,jg)/mille / dlt(jg) * xx + &
                        so_capa_ice(jst) * tmc_layt(ji,jg) / mille/dlt(jg) * (1.-xx) + &
                        pcappa_supp_pft(ji,jg,jv) 
                ELSE
                   pcapa_per_pft(ji,jg,jv) = so_capa_dry_net(ji,jg,jv)*(1-mcs(jst)) + &
                        water_capa * tmc_layt(ji,jg)/mille / dlt(jg) * xx + &
                        so_capa_ice(jst) * tmc_layt(ji,jg) / mille/dlt(jg) * (1.-xx)
                ENDIF 

                ! Add moss thermal capacity
                IF (ok_moss) THEN
                   IF (zlt(jg) .LE. soil_moss_depth) THEN
                      pcapa_moss(ji,jg) = pcapa_moss_dry + &
                           tmc_layt(ji,jg) / mille / dlt(jg) * (xx * pcapa_moss_liq + (1.-xx) * pcapa_moss_ice - pcapa_moss_dry)
                      pcapa_per_pft(ji,jg,jv) = pcapa_per_pft(ji,jg,jv) + &
                           moss_frac(jv) * pcapa_moss(ji,jg) * moss_layer_thickness / soil_moss_depth
                   ENDIF
                ENDIF
                rho_tot =  rho_soil* (1-poros_net(ji,jg,jv)) + &
                     rho_water * tmc_layt(ji,jg)/mille / dlt(jg) * xx + &
                     rho_ice * tmc_layt(ji,jg) / mille/dlt(jg) * (1.-xx)
                pcapa_spec_pft(ji, jg,jv) = pcapa_per_pft(ji, jg,jv) / rho_tot

                ! NOTE: pcappa_supp_pft is transformed from volumetric heat capacity into surfacic heat capacity (J/K/m2)
                pcappa_supp_pft(ji,jg,jv)= pcappa_supp_pft(ji,jg,jv)*zx2(ji,jg,jv)*dlt(jg)
              ENDIF
          ENDDO
       ENDDO
    ENDDO

    ! Output the specific heat capcaity for SP-MIP
    DO jv = 1,nvm
       DO jg = 1, ngrnd
          DO ji = 1,kjpindex
             pcapa_spec(ji,jg)=pcapa_spec(ji,jg) + pcapa_spec_pft(ji,jg,jv)*veget_max(ji,jv)
          ENDDO
       ENDDO
    ENDDO

    IF (ok_moss) CALL xios_orchidee_send_field("pcapa_moss",pcapa_moss)
    CALL xios_orchidee_send_field("pcapa_spec",pcapa_spec)
    !
    ! 3. Calculate the heat conductivity with allowance for permafrost
    ! Note: mc_layt has no PFT dimention,so we calculate here profil_froz_mean. Actually, profil_froz_pft along the PFT dimention currently has no difference for each PFT.
    IF ( ANY(MAXVAL(profil_froz_pft,DIM=3)>MINVAL(profil_froz_pft,DIM=3)) ) THEN
       CALL ipslerr_p(3,'thermosoil_getdiff_pft','profil_froz_mean wrong','','')
    ENDIF
    profil_froz_mean=MINVAL(profil_froz_pft,DIM=3)
    tmp = mc_layt*(1.- profil_froz_mean)
    CALL thermosoil_cond_pft (kjpindex, njsc, mc_layt, QZ, tmp, zx1,zx2, poros_net,pkappa_per_pft)

    !! Computes snow heat capacity and conductivity    
    DO ji = 1,kjpindex
       pcapa_snow(ji,:) = snowrho(ji,:) * xci


       SELECTCASE (snow_cond_method)
       CASE (SNOW_COND_METHOD_DEFAULT)
          pkappa_snow(ji,:) = (ZSNOWTHRMCOND1 + ZSNOWTHRMCOND2*snowrho(ji,:)*snowrho(ji,:)) +      &
               MAX(0.0,(ZSNOWTHRMCOND_AVAP+(ZSNOWTHRMCOND_BVAP/(snowtemp(ji,:)+ &
               ZSNOWTHRMCOND_CVAP)))*(XP00/(pb(ji)*100.)))
       CASE (SNOW_COND_METHOD_DECHARME16)
          pkappa_snow(ji,:) = 2.2*((snowrho(ji,:)/1000.)**2.0) +      &
               MAX(0.0,(ZSNOWTHRMCOND_AVAP+(ZSNOWTHRMCOND_BVAP/(snowtemp(ji,:)+ &
               ZSNOWTHRMCOND_CVAP)))*(XP00/(pb(ji)*100.)))
       CASE (SNOW_COND_METHOD_SNOWICE)
          ! 1. thermal conductivity from Yen 1981 :
          pkappa_snow(ji,:) =  cond_ice * exp(1.88*log(snowrho(ji,:)/rho_water))
          ! 2. Implicit vapor diffn effects Sun et al 1999 :     
          pkappa_snow(ji,:) = pkappa_snow(ji,:) +      &
               MAX(0.0,(ZSNOWTHRMCOND_AVAP+(ZSNOWTHRMCOND_BVAP/(snowtemp(ji,:)+ &
               ZSNOWTHRMCOND_CVAP)))*(XP00/(pb(ji)*100.)))
       CASE DEFAULT
          CALL ipslerr_p(3,'thermosoil_getdiff_pft','Unsupported SNOW_COND_METHOD', &
               'Currently supported methods are ','default(1), ducharme16(2) or snowice(3)')
       ENDSELECT


       IF ( ok_ice_sheet .and. (ice_sheet_mask(ji).GT.0) ) THEN
          pkappa_ice(ji,:) = ZICETHRMCOND1 * exp(ZICETHRMCOND2*icetemp(ji,:))  !  from Yin CHao Yen (Temp > 195K)
          pcapa_ice(ji,:) = rho_ice * (ZICETHRMHEAT1 + ZICETHRMHEAT2*(icetemp(ji,:)-tp_00)) ! formulation as in GRISLI prop_th_icetemp.f90
       ELSE
          pkappa_ice(ji,:) =  ZICETHRMCOND1 * exp(ZICETHRMCOND2*(tp_00))
          pcapa_ice(ji,:) = rho_ice * ZICETHRMHEAT1 ! formulation as in GRISLI prop_th_icetemp.f90
       ENDIF

    END DO
  END SUBROUTINE thermosoil_getdiff_pft
                                                                                                                                                                                                                             

!! ================================================================================================================================
!! SUBROUTINE   : thermosoil_getdiff_old_thermix_without_snow
!!
!>\BRIEF          Computes soil and snow heat capacity and conductivity    
!!
!! DESCRIPTION	: Calculations of soil and snow thermal properties without effect of freezing.
!!                
!!
!! RECENT CHANGE(S) : None
!! 
!! MAIN OUTPUT VARIABLE(S):
!!                          
!! REFERENCE(S) :
!!
!! FLOWCHART    : None 
!! \n 
!_ ================================================================================================================================

  SUBROUTINE thermosoil_getdiff_old_thermix_without_snow( kjpindex, njsc, snowrho, snowtemp, pb, icetemp,&
       ice_sheet_mask)

    !! 0. Variables and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                          :: kjpindex
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)    :: njsc     !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless)
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(in) :: snowrho  !! Snow density (Kg/m^3)
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(in) :: snowtemp !! Snow temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: pb       !! Surface presure (hPa)
    REAL(r_std), DIMENSION (kjpindex,nice), INTENT(in)  :: icetemp  !! Ice temperature (K)
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)    :: ice_sheet_mask !! Ice sheet mask

    !! 0.1 Local variables
    INTEGER(i_std)    				        :: ji, jg, jst    !! Index
    REAL(r_std), DIMENSION (kjpindex,ngrnd)             :: pcapa_tmp      !! Soil heat capacity (J/m3/K) 

    !! Computes soil heat capacity and conductivity
    DO jg = 1,ngrnd
       DO ji = 1,kjpindex
          jst = njsc(ji)
          pcapa_tmp(ji, jg) = so_capa_dry(jst) * (1-mcs(jst)) + water_capa * tmc_layt(ji,jg)/mille/dlt(jg)
          pcapa(ji,jg) = pcapa_tmp(ji, jg)
          pcapa_en(ji,jg) = pcapa_tmp(ji, jg)
       ENDDO
    ENDDO


    CALL thermosoil_cond (kjpindex, njsc, mc_layt, QZ, mcl_layt, pkappa)


    IF (brk_flag == 1) THEN
       ! Bedrock flag is activated
       DO jg = ngrnd-1,ngrnd
          DO ji = 1,kjpindex
             pcapa(ji,jg) = brk_capa
             pcapa_en(ji,jg) = brk_capa 
             pkappa(ji,jg) = brk_cond
          ENDDO
       ENDDO
    ENDIF

    !! Computes snow heat capacity and conductivity
    DO ji = 1,kjpindex
       pcapa_snow(ji,:) = snowrho(ji,:) * xci

       SELECTCASE (snow_cond_method)
       CASE (SNOW_COND_METHOD_DEFAULT)
          pkappa_snow(ji,:) = (ZSNOWTHRMCOND1 + ZSNOWTHRMCOND2*snowrho(ji,:)*snowrho(ji,:)) +      &
               MAX(0.0,(ZSNOWTHRMCOND_AVAP+(ZSNOWTHRMCOND_BVAP/(snowtemp(ji,:)+ &
               ZSNOWTHRMCOND_CVAP)))*(XP00/(pb(ji)*100.)))
       CASE (SNOW_COND_METHOD_DECHARME16)
          pkappa_snow(ji,:) = 2.2*((snowrho(ji,:)/1000.)**2.0) +      &
               MAX(0.0,(ZSNOWTHRMCOND_AVAP+(ZSNOWTHRMCOND_BVAP/(snowtemp(ji,:)+ &
               ZSNOWTHRMCOND_CVAP)))*(XP00/(pb(ji)*100.)))
       CASE (SNOW_COND_METHOD_SNOWICE)
          ! 1. thermal conductivity from Yen 1981 :
          pkappa_snow(ji,:) =  cond_ice * exp(1.88*log(snowrho(ji,:)/rho_water))
          ! 2. Implicit vapor diffn effects Sun et al 1999 :     
          pkappa_snow(ji,:) = pkappa_snow(ji,:) +      &
               MAX(0.0,(ZSNOWTHRMCOND_AVAP+(ZSNOWTHRMCOND_BVAP/(snowtemp(ji,:)+ &
               ZSNOWTHRMCOND_CVAP)))*(XP00/(pb(ji)*100.)))
       CASE DEFAULT
          CALL ipslerr_p(3,'thermosoil_getdiff_pft','Unsupported SNOW_COND_METHOD', &
               'Currently supported methods are ','default(1), ducharme16(2) or snowice(3)')
       ENDSELECT

       IF ( ok_ice_sheet .and. (ice_sheet_mask(ji).GT.0) ) THEN
          pkappa_ice(ji,:) = ZICETHRMCOND1 * exp(ZICETHRMCOND2*icetemp(ji,:))  !  from Yin CHao Yen (Temp > 195K)
          pcapa_ice(ji,:) = rho_ice * (ZICETHRMHEAT1 + ZICETHRMHEAT2*(icetemp(ji,:)-tp_00)) ! formulation as in GRISLI prop_th_icetemp.f90
       ELSE
          pkappa_ice(ji,:) =  ZICETHRMCOND1 * exp(ZICETHRMCOND2*(tp_00))
          pcapa_ice(ji,:) = rho_ice * ZICETHRMHEAT1 ! formulation as in GRISLI prop_th_icetemp.f90
       ENDIF
    END DO


  END SUBROUTINE thermosoil_getdiff_old_thermix_without_snow


!! ================================================================================================================================
!! SUBROUTINE   : thermosoil_read_reftempfile
!!
!>\BRIEF          
!!
!! DESCRIPTION	: Read file with longterm soil temperature
!!                
!!
!! RECENT CHANGE(S) : None
!! 
!! MAIN OUTPUT VARIABLE(S): reftemp : Reference temerature
!!                          
!! REFERENCE(S) :
!!
!! FLOWCHART    : None 
!! \n 
!_ ================================================================================================================================
  SUBROUTINE thermosoil_read_reftempfile(kjpindex,lalo,reftemp)

    USE interpweight

    IMPLICIT NONE

    !! 0. Variables and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                           :: kjpindex
    REAL(r_std), DIMENSION(kjpindex,2), INTENT(in)       :: lalo

    !! 0.2 Output variables
    REAL(r_std), DIMENSION(kjpindex, ngrnd), INTENT(out) :: reftemp

    !! 0.3 Local variables
    INTEGER(i_std) :: ib
    CHARACTER(LEN=80) :: filename
    REAL(r_std),DIMENSION(kjpindex) :: reftemp_file                          !! Horizontal temperature field interpolated from file [C]
    INTEGER(i_std),DIMENSION(kjpindex,8) :: neighbours
    REAL(r_std)                                          :: vmin, vmax       !! min/max values to use for the
    !!   renormalization
    REAL(r_std), DIMENSION(kjpindex)                     :: areftemp         !! Availability of data for  the interpolation
    CHARACTER(LEN=80)                                    :: variablename     !! Variable to interpolate
    !!   the file
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
    !Config If    = READ_REFTEMP
    !Config Def   = reftemp.nc
    !Config Help  = 
    !Config Units = [FILE]
    filename = 'reftemp.nc'
    CALL getin_p('REFTEMP_FILE',filename)

    variablename = 'temperature'

    IF (printlev >= 1) WRITE(numout,*) "thermosoil_read_reftempfile: Read and interpolate file " &
         // TRIM(filename) //" for variable " //TRIM(variablename)

    IF (xios_interpolation) THEN

       CALL xios_orchidee_recv_field('reftemp_interp',reftemp_file)

       DO ib=1, kjpindex
          reftemp(ib,:) = reftemp_file(ib) + ZeroCelsius 
       END DO
       areftemp = 1.0
    ELSE


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
            reftemp_file, areftemp)
       IF (printlev >= 5) WRITE(numout,*)'  thermosoil_read_reftempfile after interpweight_2Dcont'

       ! Copy reftemp_file temperature to all ground levels and transform into Kelvin
       DO ib=1, kjpindex
          reftemp(ib, :) = reftemp_file(ib)+ZeroCelsius
       END DO

    END IF

    ! Write diagnostics
    CALL xios_orchidee_send_field("interp_avail_areftemp",areftemp)
    CALL xios_orchidee_send_field("interp_diag_reftemp",reftemp_file)

  END SUBROUTINE thermosoil_read_reftempfile

!! ================================================================================================================================
!! SUBROUTINE   : thermosoil_read_refsoc_file
!!
!>\BRIEF          
!!
!! DESCRIPTION : Read file of soil organic carbon to be used in thermix
!! (insulating effect)
!!                
!!
!! RECENT CHANGE(S) : None
!! 
!! MAIN OUTPUT VARIABLE(S): refsoc : soil organic carbon from data
!!                          
!! REFERENCE(S) :
!!
!! FLOWCHART    : None 
!! \n 
!_ ================================================================================================================================
    
  SUBROUTINE thermosoil_read_refsoc_file(nbpt, lalo, neighbours, resolution, contfrac)

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                    :: nbpt                  !! Number of points for which the data needs to be interpolated (unitless)             
    REAL(r_std), INTENT(in)                       :: lalo(nbpt,2)          !! Vector of latitude and longitudes (degree)        
    INTEGER(i_std), INTENT(in)                    :: neighbours(nbpt,NbNeighb)!! Vector of neighbours for each grid point (1=N,2=E,3=S,4=W)  
    REAL(r_std), INTENT(in)                       :: resolution(nbpt,2)    !! The size of each grid cell in X and Y (km)
    REAL(r_std), INTENT(in)                       :: contfrac(nbpt)        !! Fraction of land in each grid cell (unitless)   

    !! 0.4 Local variables

    INTEGER(i_std)                                :: nbvmax                !! nbvmax for interpolation (unitless) 
    CHARACTER(LEN=80)                             :: filename              !! Name of file to read
    CHARACTER(LEN=80)                             :: variablename          !! Variable to interpolate
    
    INTEGER(i_std)                                :: iml, jml, lml, tml    !! Indices
    INTEGER(i_std)                                :: fid, ib, ip, jp, fopt !! Indices
    INTEGER(i_std)                                :: ilf, ks               !! Indices
    REAL(r_std)                                   :: totarea               !! Help variable to compute average SOC
    REAL(r_std), ALLOCATABLE, DIMENSION(:)        :: lat_lu, lon_lu        !! Latitudes and longitudes read from input file
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)      :: lat_rel, lon_rel      !! Help variable to read file data and allocate memory
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)      :: mask_lu               !! Help variable to read file data and allocate memory
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:)   :: mask
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:)    :: refsoc_file           !! Help variable to read file data and allocate memory
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)      :: sub_area              !! Help variable to read file data and allocate memory
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:,:) :: sub_index             !! Help variable to read file data and allocate memory
    CHARACTER(LEN=30)                             :: callsign              !! Help variable to read file data and allocate memory
    CHARACTER(LEN=100)                            :: str                   !! Temporary string var 
    LOGICAL                                       :: ok_interpol           !! Optional return of aggregate_2d
    INTEGER                                       :: ALLOC_ERR             !! Help varialbe to count allocation error
    !_
    !================================================================================================================================

    !! 1. Open file and allocate memory

    ! Open file with SOC map

    !Config Key   = SOIL_REFSOC_FILE
    !Config Desc  = File with soil carbon stocks
    !Config If    = OK_SOIL_CARBON_DISCRETIZATION, USE_REFSOC, SOIL_CTEMPDIFF
    !Config Def   = refSOC.nc
    !Config Help  = 
    !Config Units = [FILE]
    filename = 'refSOC.nc'
    CALL getin_p('SOIL_REFSOC_FILE',filename)
    variablename = 'soil_organic_carbon'

    IF (printlev >= 1) WRITE(numout,*) "thermosoil_read_refsoc_file: Read and interpolate " &
         // TRIM(filename) // " for variable " // TRIM(variablename)

    IF (xios_interpolation) THEN
       ! Read refsoc using XIOS. Information about the input file is set in the xml files.
       CALL xios_orchidee_recv_field('refsoc_interp',refsoc)
       
    ELSE
       ! Read data from file
       IF (is_root_prc) CALL flininfo(filename, iml, jml, lml, tml, fid)
       CALL bcast(iml)
       CALL bcast(jml)
       CALL bcast(lml)
       CALL bcast(tml)
       
       IF (lml .NE. ngrnd) THEN
          WRITE(str, *) 'ngrnd=', ngrnd, ', depth found in file=', lml
          CALL ipslerr_p(3, 'thermosoil_read_refsoc_file',    &
               'depth from the file must be the same as ngrnd', &
               str, &
            filename )
       ENDIF
       
       ALLOCATE(lon_lu(iml), STAT=ALLOC_ERR)
       IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'thermosoil_read_refsoc_file','Problem in allocation of variable lon_lu','','')
       
       ALLOCATE(lat_lu(jml), STAT=ALLOC_ERR)
       IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'thermosoil_read_refsoc_file','Problem in allocation of variable lat_lu','','')
       
       ALLOCATE(mask_lu(iml,jml), STAT=ALLOC_ERR)
       IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'thermosoil_read_refsoc_file','Pb in allocation for mask_lu','','')
       
       ALLOCATE(refsoc_file(iml,jml,lml), STAT=ALLOC_ERR)
       IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'thermosoil_read_refsoc_file','Pb in allocation for refsoc_file','','')
       
       IF (is_root_prc) THEN
          CALL flinget(fid, 'longitude', iml, 0, 0, 0, 1, 1, lon_lu)
          CALL flinget(fid, 'latitude', jml, 0, 0, 0, 1, 1, lat_lu)
          CALL flinget(fid, 'mask', iml, jml, 0, 0, 1, 1, mask_lu)
          CALL flinget(fid, 'soil_organic_carbon', iml, jml, lml, tml, 1, 1, refsoc_file)
          
          CALL flinclo(fid)
       ENDIF
       
       CALL bcast(lon_lu)
       CALL bcast(lat_lu)
       CALL bcast(mask_lu)
       CALL bcast(refsoc_file)
       
       ALLOCATE(lon_rel(iml,jml), STAT=ALLOC_ERR)
       IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'thermosoil_read_refsoc_file','Pb in allocation for lon_rel','','')
       
       ALLOCATE(lat_rel(iml,jml), STAT=ALLOC_ERR)
       IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'thermosoil_read_refsoc_file','Pb in allocation for lat_rel','','')
       
       ALLOCATE(mask(iml,jml), STAT=ALLOC_ERR)
       IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'thermosoil_read_refsoc_file','Problem in allocation of variable mask','','')
       
       DO jp=1,jml
          lon_rel(:,jp) = lon_lu(:)
       ENDDO
       DO ip=1,iml
          lat_rel(ip,:) = lat_lu(:)
       ENDDO
       
       mask(:,:) = zero
       WHERE (mask_lu(:,:) > zero )
          mask(:,:) = un
       ENDWHERE
       
       ! Set nbvmax to 200 for interpolation
       ! This number is the dimension of the variables in which we store 
       ! the list of points of the source grid which fit into one grid box of the
       ! target. 
       nbvmax = 200 
       callsign = 'soil organic carbon'
       
       ! Start interpolation
       ok_interpol=.FALSE.
       DO WHILE ( .NOT. ok_interpol )
          WRITE(numout,*) "Projection arrays for ",callsign," : "
          WRITE(numout,*) "nbvmax = ",nbvmax
          
          ALLOCATE(sub_area(nbpt,nbvmax), STAT=ALLOC_ERR)
          IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'thermosoil_read_refsoc_file','Pb in allocation for sub_area','','')
          sub_area(:,:)=zero
          
          ALLOCATE(sub_index(nbpt,nbvmax,2), STAT=ALLOC_ERR)
          IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'thermosoil_read_refsoc_file','Pb in allocation for sub_index','','')
          sub_index(:,:,:)=0
          
          CALL aggregate_p(nbpt, lalo, neighbours, resolution, contfrac, &
               iml, jml, lon_rel, lat_rel, mask, callsign, &
               nbvmax, sub_index, sub_area, ok_interpol)
          
          IF ( .NOT. ok_interpol ) THEN
             DEALLOCATE(sub_area)
             DEALLOCATE(sub_index)
             nbvmax = nbvmax * 2
          ENDIF
       ENDDO
       
       ! Compute the average
       refsoc(:,:) = zero
       DO ib = 1, nbpt
          fopt = COUNT(sub_area(ib,:) > zero)
          IF ( fopt > 0 ) THEN
             totarea = zero
             DO ilf = 1, fopt
                ip = sub_index(ib,ilf,1)
                jp = sub_index(ib,ilf,2)
                refsoc(ib,:) = refsoc(ib,:) + refsoc_file(ip,jp,:) * sub_area(ib,ilf)
                totarea = totarea + sub_area(ib,ilf)
             ENDDO
             ! Normalize
             refsoc(ib,:) = refsoc(ib,:)/totarea
          ELSE
             ! Set defalut value for points where the interpolation fail
             WRITE(numout,*) 'On point ', ib, ' no points were found for interpolation data. Mean value is used.'
             WRITE(numout,*) 'Location : ', lalo(ib,2), lalo(ib,1)
             refsoc(ib,:) = 0.
          ENDIF
       ENDDO
       
       DEALLOCATE (lat_lu)
       DEALLOCATE (lat_rel)
       DEALLOCATE (lon_lu)
       DEALLOCATE (lon_rel)
       DEALLOCATE (mask_lu)
       DEALLOCATE (mask)
       DEALLOCATE (refsoc_file)
       DEALLOCATE (sub_area)
       DEALLOCATE (sub_index)
    END IF
    
    ! Write diagnostics
    CALL xios_orchidee_send_field("interp_diag_refsoc",refsoc)
   
  END SUBROUTINE thermosoil_read_refsoc_file


!! 
!================================================================================================================================ 
!! SUBROUTINE   : thermosoil_add_heat_zimov 
!! 
!>\BRIEF          heat 
!! 
!! DESCRIPTION  :  
!! 
!! RECENT CHANGE(S) : None 
!! 
!! MAIN OUTPUT VARIABLE(S): 
!! 
!! REFERENCE(S) : 
!! 
!! FLOWCHART    : None 
!! \n 
!_ 
!================================================================================================================================ 
  SUBROUTINE thermosoil_add_heat_zimov(kjpindex, veget_max, ptn, heat_Zimov)
    !! 0. Variables and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std),INTENT(in)                                 :: kjpindex
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(in)         :: veget_max    !! Fraction of vegetation type 

    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm), INTENT (in)   :: heat_Zimov   !! heating associated with decomposition [W/m**3 soil]

    !! 0.2 Modified variables
    REAL(r_std),DIMENSION(kjpindex,ngrnd,nvm),INTENT(inout)   :: ptn          !! Soil temperature profile (K)

    !! 0.3 Local variables
    INTEGER(r_std) :: ji, jg, jv 

    IF (printlev>=3) WRITE (numout,*) 'entering thermosoil_add_heat_zimov'

    DO ji = 1, kjpindex
       DO jv = 1,nvm
          DO jg = 1, ngrnd
             ! BG: should we use pcapa without pft dimension or pkappa_per_pft since they seem similar?
             ptn(ji,jg,jv) = ptn(ji,jg,jv) + heat_zimov(ji,jg,jv) * dt_sechiba / ( pcapa(ji,jg) * dlt(jg) )
          END DO
       END DO
    END DO

    ! ptn_pftmean needs to be updated to ensure consistency
    ptn_pftmean(:,:) = zero
    DO jv=1,nvm
       DO jg = 1, ngrnd
          ptn_pftmean(:,jg) = ptn_pftmean(:,jg) + ptn(:,jg,jv) * veget_max(:,jv)
       ENDDO ! jg = 1, ngrnd
    ENDDO ! m=1,nvm

    IF (printlev>=3) WRITE (numout,*) ' thermosoil_add_heat_zimov done'

  END SUBROUTINE thermosoil_add_heat_zimov

!! ================================================================================================================================
!! SUBROUTINE   : thermosoil_getdiff_pft_newdisc
!!
!>\BRIEF          Computes soil and snow heat capacity and conductivity when the soil carbon is discretized, with the new thermal param   
!!
!! DESCRIPTION  : Computation of the soil thermal properties; snow properties are also accounted for
!!
!! RECENT CHANGE(S) : None
!! 
!! MAIN OUTPUT VARIABLE(S): Following module variables are calculted in this subroutine: 
!!                          profil_froz_pft, pcapa_per_pft, pkappa_per_pft, pcapa_en_per_pft, tmc_layt_pft
!!                          
!! REFERENCE(S) :
!!
!! FLOWCHART    : None 
!! \n 
!_ ================================================================================================================================

  SUBROUTINE thermosoil_getdiff_pft_newdisc( kjpindex,  ptn,     njsc,     shum_ngrnd_perma, &
                                     depth_organic_soil, som_total, snowrho, snowtemp, pb, veget_max)

    !! 0. Variables and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std),INTENT(in)                                          :: kjpindex
    REAL(r_std),DIMENSION(kjpindex,ngrnd),INTENT(in)                   :: shum_ngrnd_perma     !! Water saturation degree on the soil depth axes define by the thermic (0-1, dimensionless)
    REAL(r_std), DIMENSION(kjpindex), INTENT (in)                      :: depth_organic_soil   !! Depth at which there is still organic matter (m)
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm,nelements), INTENT (in)  :: som_total            !! total soil carbon for use in thermal calcs (g/m**3)
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)                   :: njsc                 !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless)
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(in)                :: snowrho              !! Snow density (Kg/m^3)
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(in)                :: snowtemp             !! Snow temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)                      :: pb                   !! Surface presure (hPa)
    REAL(r_std),DIMENSION(kjpindex,ngrnd,nvm),INTENT(in)               :: ptn                  !! Soil temperature profile 
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(in)                   :: veget_max            !! Fraction of vegetation type
    !! 0.2 Modified variables

    !! 0.3 Output variables

    !! 0.3 Local variables
    INTEGER(i_std)                                           :: ji,jg,jv
    INTEGER(i_std)                                           :: jst
    REAL(r_std)                                              :: xx             !! Unfrozen fraction of the soil
    REAL(r_std)                                              :: p
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm)               :: so_capa_dry_net
    REAL(r_std), DIMENSION(nscm)                             :: so_capa_min_solid !! Solid heat capacity of mineral soil, derived from dry heat capacity and porosity, J.m^{-3}.K^{-1}
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm)               :: poros_net
    REAL(r_std), DIMENSION(kjpindex,ngrnd)                   :: profil_froz_mean, tmp
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm)               :: pcappa_supp_pft
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm)               :: pcapa_spec_pft !! SPECIFIC soil heat capacity (J/kg/K)
    REAL(r_std)                                              :: rho_tot        !! Soil density (kg/m3)
    REAL(r_std), DIMENSION(kjpindex,ngrnd)                   :: pcapa_spec     !! SPECIFIC soil heat capacity (J/kg/K)
    REAL(r_std), DIMENSION(kjpindex,ngrnd)                   :: poros_net_glob !! Porosity based on SOC content over all PFTs (A. Cuynet)
    REAL(r_std), DIMENSION(kjpindex,ngrnd)                   :: pcapa_moss     !! moss heat capacity (J/m3/K)
    REAL(r_std), DIMENSION(kjpindex,ngrnd)                   :: som_total_pftmean !! Deep soil organic matter profiles calculated for thermosoil, mean over all PFTs in concentration g/m3
    REAL(r_std), DIMENSION(kjpindex,ngrnd)                   :: fOM            !! Mass fraction of soil organic matter [kg/kg]
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm)               :: F_OM           !! Solid volume fraction of soil organic matter (does not include the air/porosity) [m3/m3]
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm)               :: F_MM           !! Solid volume fraction of mineral matter [m3/m3]
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm)               :: F_CM, F_air    !! Solid volume fractions of coarse matter and air [m3/m3]

    
    ! fSOC, BD_FM and Phi_CM come directly from the maps
    !
    ! Default: organic layer not taken into account
    fOM(:,:) = 0.0
    F_OM(:,:,:) = 0.0
    F_MM(:,:,:) = 0.0
    F_CM(:,:,:) = 0.0
    
    ! Calculation of solid mineral heat capacities from dry capacities
    so_capa_min_solid(:) = so_capa_dry(:)/(1.-mcs(:))
    
    poros_net(:,:,:) = 0.0
    poros_net_glob(:,:) = 0.0
    pcapa_spec_pft(:,:,:) = 0.0 
    pcapa_spec(:,:) = 0.0
    pcapa_moss(:,:) = 0.0
    som_total_pftmean(:,:) = 0.0
    
    ! Calculation of the volume fractions of each component
    IF ( use_soilc_insulation ) THEN
       IF (use_refsoc .OR. impose_soc_bd) THEN
          DO jv = 1,nvm
             DO jg = 1, ngrnd
                DO ji = 1,kjpindex
                   jst = njsc(ji)
                   fOM(ji,jg) = MIN(1.0, fSOC(ji,jg)/alpha_VB) ! Calculate fraction of organic matter based on SOC content in fine matter
                   F_OM(ji,jg,jv) = MIN(1.0 - poros_org_mini, fOM(ji,jg)*BD_FM(ji,jg)*(1.-Phi_CM(ji,jg))/(rho_orga)) ! Volume fraction of organic matter in the whole pixel
                   F_MM(ji,jg,jv) = MIN(1.0 - mcs(jst), (1.-fOM(ji,jg))*BD_FM(ji,jg)*(1.-Phi_CM(ji,jg))/(rho_soil)) ! Volume fraction of mineral matter in the whole pixel
                   F_CM(ji,jg,jv) = MIN(1.0 - mcs(jst), Phi_CM(ji,jg)*(1.0-mcs(jst))) ! Volume fraction of coarse matter in the whole pixel
                ENDDO
                som_total_pftmean(:,jg) = som_total_pftmean(:,jg) + F_OM(:,jg,jv) * rho_orga * veget_max(:,jv) * 1000 ! g/m3
             ENDDO
          ENDDO
       ELSE !uses som from stomate
          DO jv = 1,nvm
             DO jg = 1, ngrnd
                DO ji = 1,kjpindex
                   jst = njsc(ji)
                   F_OM(ji,jg,jv) = MIN(1.0 - poros_org, som_total(ji,jg,jv,icarbon)/(rho_orga*1000))
                   F_MM(ji,jg,jv) = MIN(1.0 - mcs(jst), (1. - mcs(jst))*(1. - (F_OM(ji,jg,jv) / (1. - poros_org))))
                ENDDO
                som_total_pftmean(:,jg) = som_total_pftmean(:,jg) + som_total(:,jg,jv,icarbon) * veget_max(:,jv) * 1000
             ENDDO
          ENDDO
       ENDIF
    ELSE
       !
       ! level 1
       !
       DO jv = 1,nvm
          DO ji = 1,kjpindex
             jst = njsc(ji)
             IF ( depth_organic_soil(ji) .GT. zlt(1) ) THEN
                !! the 1st level is in the organic => the 1st layer is entirely organic
                F_OM(ji,1,jv) = 1. - poros_org !! F_OM represents the volume occupied by organic matter and the rest corresponds to the pores
                F_MM(ji,1,jv) = 0.
             ELSE IF ( depth_organic_soil(ji) .GT. zero ) THEN
                !! the 1st level is beyond the organic and the organic is present
                F_OM(ji,1,jv) = (depth_organic_soil(ji) / zlt(1))*(1. - poros_org)
                F_MM(ji,1,jv) = (1. - mcs(jst))*(1. - (F_OM(ji,1,jv) / (1. - poros_org)))
             ELSE
                ! there is no organic at all
                F_OM(ji,1,jv) = 0.
                F_MM(ji,1,jv) = 1. - mcs(jst)
             ENDIF
          ENDDO
       ENDDO
       !
       ! other levels
       !
       DO jg = 2, ngrnd !- 2
          DO ji = 1,kjpindex
             IF ( depth_organic_soil(ji) .GT. zlt(jg) ) THEN
                ! the current level is in the organic => the current layer is
                ! entirely organic
                F_OM(ji,jg,1) = 1. - poros_org
             ELSE IF ( depth_organic_soil(ji) .GT. zlt(jg-1) ) THEN
                ! the current layer is partially organic
                F_OM(ji,jg,1) = ((depth_organic_soil(ji) - zlt(jg-1)) / (zlt(jg) - zlt(jg-1))) * (1.-poros_org)
             ELSE
                ! both levels are out of organic => the current layer is entirely
                ! mineral soil       
                F_OM(ji,jg,1) = 0.
             ENDIF
             F_MM(ji,jg,1) = (1. - mcs(jst))*(1. - (F_OM(ji,jg,1) / (1. - poros_org)))
          ENDDO
       ENDDO
       DO jv = 2, nvm
          F_OM(:,:,jv) = F_OM(:,:,1)
          F_MM(:,:,jv) = F_MM(:,:,1)
       ENDDO
       
       DO jg = 1, ngrnd
          som_total_pftmean(:,jg) = som_total_pftmean(:,jg) + F_OM(:,jg,1) * rho_orga * 1000 ! g/m3
       ENDDO
    ENDIF ! ( use_soilc_insulation )
 
     
    CALL xios_orchidee_send_field ("deepC_thermo_pftmean", som_total_pftmean)
    CALL xios_orchidee_send_field ("F_OM", F_OM(:,:,1)) ! Just get F_OM for PFT1 for the test (later on, we will have to calculate the mean value)
    
#ifdef STRICT_CHECK
    IF (ANY(tmc_layt_pft < min_sechiba)) CALL ipslerr_p(3, "thermosoil_getdiff_pft_newdisc", "tmc_layt_pft has negative values", "", "") ! prec issues
#endif

    DO jv = 1,nvm
       DO jg = 1, ngrnd
          DO ji = 1,kjpindex
             jst = njsc(ji)
             !
             ! 1. Calculate dry heat capacity and conductivity, taking
             ! into account the organic and mineral fractions in the layer
             !
             ! Former MICT version
             ! Here we take into account the new dependance of the soil heat capacity from the soil type.
             ! so_capa_dry_net(ji,jg,jv) = zx1(ji,jg,jv) * SO_CAPA_DRY_ORG + zx2(ji,jg,jv) * so_capa_dry(jst)
             ! MODIFICATION OF THE CALCULATION ! Here, instead of using zx1 that is the volume fraction including
             ! the air/porosity, we use F_OM, which represents only the solid organic matter (no pores)
             so_capa_dry_net(ji,jg,jv) = F_OM(ji,jg,jv) * SO_CAPA_DRY_ORG + (F_CM(ji,jg,jv)+F_MM(ji,jg,jv)) * so_capa_min_solid(jst)
             
             !Here we take into account the new dependance of the porosity from the soil type.
             ! poros_net(ji,jg,jv) = zx1(ji,jg,jv) * poros_org + zx2(ji,jg,jv) * mcs(jst)
             ! MODIFICATION ALSO HERE
             poros_net(ji,jg,jv) = 1. - F_OM(ji,jg,jv) - F_MM(ji,jg,jv) - F_CM(ji,jg,jv)
             
             ! 2. Calculate heat capacity with allowance for permafrost
          ENDDO
       ENDDO
    ENDDO
    
    ! Output the global porosity (sum over every PFT)
    DO jv = 1,nvm
       DO jg = 1, ngrnd
          poros_net_glob(:,jg)=poros_net_glob(:,jg) + (poros_net(:,jg,jv)*veget_max(:,jv)/(SUM(veget_max(:,:),dim=2)))
       ENDDO
    ENDDO
    
    CALL xios_orchidee_send_field("poros_net_thermosoil",poros_net_glob)
    
    DO jv = 1,nvm
       DO jg = 1, ngrnd
          DO ji = 1,kjpindex
             jst = njsc(ji)
             F_air(ji,jg,jv) = MAX(0., 1. - F_OM(ji,jg,jv) - F_MM(ji,jg,jv) - F_CM(ji,jg,jv) - tmc_layt(ji,jg) / mille / dlt(jg))
             IF (ptn(ji,jg,jv) .LT. fr_center-fr_dT/2.) THEN
                ! frozen soil
                ! xx is the unfrozen fraction of soil water (in volume)
                xx = 0.
                pcappa_supp_pft(ji,jg, jv)= 0.
             ELSEIF (ptn(ji,jg,jv) .GT. fr_center+fr_dT/2.) THEN
                xx = 1.
                pcappa_supp_pft(ji,jg,jv)= 0.
             ELSE
                ! Latent heat contribution
                ! Here, I keep mcs so that the results stay consistent with the hydrology (may need to be changed
                ! if the porosity is changed later on in hydrol)
                pcappa_supp_pft(ji,jg,jv)= shum_ngrnd_perma(ji,jg)*mcs(jst)*lhf*rho_water/fr_dT       
                xx = (ptn(ji,jg,jv)-(fr_center-fr_dT/2.)) / fr_dT
             ENDIF
             profil_froz_pft(ji,jg,jv) = (1. - xx)
             pcapa_per_pft(ji,jg,jv) = so_capa_dry_net(ji,jg,jv) + &
                  water_capa * tmc_layt(ji,jg)/mille / dlt(jg) * xx + &
                  so_capa_ice(jst) * tmc_layt(ji,jg) / mille/dlt(jg) * (1.-xx)+ &
                  air_capa * F_air(ji,jg,jv)
             
             IF (ok_freeze_thaw_latent_heat) THEN
                pcapa_per_pft(ji,jg,jv) = pcapa_per_pft(ji,jg,jv) + pcappa_supp_pft(ji,jg,jv) 
             ENDIF
             
             IF (ok_moss) THEN
                IF (zlt(jg) .LE. soil_moss_depth) THEN
                   pcapa_moss(ji,jg) = pcapa_moss_dry + &
                        tmc_layt(ji,jg) / mille / dlt(jg) * (xx * pcapa_moss_liq + (1.-xx) * pcapa_moss_ice - pcapa_moss_dry)
                   pcapa_per_pft(ji,jg,jv) = pcapa_per_pft(ji,jg,jv) + &
                        moss_frac(jv) * pcapa_moss(ji,jg) * moss_layer_thickness / soil_moss_depth
                ENDIF
             ENDIF
             
             ! rho_tot : bulk density of the soil including ice and water
             rho_tot = rho_orga * F_OM(ji,jg,jv) + rho_soil* (F_MM(ji,jg,jv) + F_CM(ji,jg,jv)) + &
                  rho_water * tmc_layt(ji,jg)/mille / dlt(jg) * xx + &
                  rho_ice * tmc_layt(ji,jg) / mille/dlt(jg) * (1.-xx)
             pcapa_spec_pft(ji, jg, jv) = pcapa_per_pft(ji, jg, jv) / rho_tot
             
             
             ! NOTE: pcappa_supp_pft is transformed from volumetric heat capacity into surfacic heat capacity (J/K/m2)
             ! zx2 here is related to the snow and this variable seems deprecated since explicit snow. zx2 is supposed to be equal to 1.
             ! It is removed.
             ! pcappa_supp_pft(ji,jg,jv)= pcappa_supp_pft(ji,jg,jv)*zx2(ji,jg,jv)*dlt(jg)
             pcappa_supp_pft(ji,jg,jv)= pcappa_supp_pft(ji,jg,jv)*dlt(jg)
          ENDDO
       ENDDO
    ENDDO
    
    ! Output the specific heat capacity for SP-MIP
    DO jv = 1,nvm
       DO jg = 1, ngrnd
          DO ji = 1,kjpindex
             pcapa_spec(ji,jg)=pcapa_spec(ji,jg) + pcapa_spec_pft(ji,jg,jv)*veget_max(ji,jv)
          ENDDO
       ENDDO
    ENDDO
    
    IF (ok_moss) CALL xios_orchidee_send_field("pcapa_moss",pcapa_moss)
    CALL xios_orchidee_send_field("pcapa_spec",pcapa_spec)
    !
    ! 3. Calculate the heat conductivity with allowance for permafrost
    ! Note: mc_layt has no PFT dimension,so we calculate here profil_froz_mean. Actually, profil_froz_pft along the PFT dimention currently has no difference for each PFT.
    IF ( ANY(MAXVAL(profil_froz_pft,DIM=3)>MINVAL(profil_froz_pft,DIM=3)) ) THEN
       CALL ipslerr_p(3,'thermosoil_getdiff_pft_newdisc','profil_froz_mean wrong','','')
    ENDIF
    profil_froz_mean=MINVAL(profil_froz_pft,DIM=3)
    ! tmp : volumetric liquid soil moisture
    ! Is the residual water also frozen ??? If not, this need to be corrected.
    tmp = mc_layt*(1.- profil_froz_mean)
    
    CALL thermosoil_cond_pft_newdisc (kjpindex, njsc, mc_layt, QZ, tmp, F_OM, F_MM, F_CM, poros_net,pkappa_per_pft)
    
    !! Computes snow heat capacity and conductivity    
    DO ji = 1,kjpindex
       pcapa_snow(ji,:) = snowrho(ji,:) * xci
       
       SELECT CASE (snow_cond_method)
       CASE (SNOW_COND_METHOD_DEFAULT)
          pkappa_snow(ji,:) = (ZSNOWTHRMCOND1 + ZSNOWTHRMCOND2*snowrho(ji,:)*snowrho(ji,:)) +      &
               MAX(0.0,(ZSNOWTHRMCOND_AVAP+(ZSNOWTHRMCOND_BVAP/(snowtemp(ji,:)+ &
               ZSNOWTHRMCOND_CVAP)))*(XP00/(pb(ji)*100.)))
       CASE (SNOW_COND_METHOD_DECHARME16)
          pkappa_snow(ji,:) = 2.2*((snowrho(ji,:)/1000.)**2.0) +      &
               MAX(0.0,(ZSNOWTHRMCOND_AVAP+(ZSNOWTHRMCOND_BVAP/(snowtemp(ji,:)+ &
               ZSNOWTHRMCOND_CVAP)))*(XP00/(pb(ji)*100.)))
       CASE DEFAULT
          CALL ipslerr_p(3,'thermosoil_getdiff_pft_newdisc','Unsupported SNOW_COND_METHOD', &
               'Currently supported methods are ','default(1) or ducharme16(2)')
       END SELECT
       
    END DO
  END SUBROUTINE thermosoil_getdiff_pft_newdisc
  
  
   
!================================================================================================================================
!! SUBROUTINE   : thermosoil_cond_pft_newdisc
!!
!>\BRIEF          Calculate soil thermal conductivity.  
!!
!! DESCRIPTION  : This routine computes soil thermal conductivity
!!                but considers the fact that soil organic carbon can decrease
!!                conductivity. New version (NewDisC 2023) 
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): cnd
!!
!! REFERENCE(S) :
!!    Farouki, O.T.,1986: Thermal Properties of Soils. Series on Rock
!!            and Soil Mechanics, Vol. 11, Trans Tech, 136 PP.
!!    Johansen, O., 1975: Thermal Conductivity of Soils. Ph.D. Thesis,
!!            University of Trondheim,
!!    Peters-Lidard, C. D., Blackburn, E., Liang, X., & Wood, E. F.,
!!            1998: The effect of soil thermal conductivity 
!!            Parameterization on Surface Energy fluxes
!!            and Temperatures. J. of The Atmospheric Sciences,
!!            Vol. 55, pp. 1209-1224.
!!    Lawrence and Slater,2008: Incorporating organic soil into a global climate
!!            model
!! Modify history:
!!
!! FLOWCHART    : None
!! \n
!_
!================================================================================================================================

  SUBROUTINE thermosoil_cond_pft_newdisc (kjpindex, njsc, smc, qz, sh2o, F_OM, F_MM, F_CM, porosnet, cnd)
    
    !! 0. Variables and parameter declaration
    
    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                                 :: kjpindex      !! Domain size (unitless)
    INTEGER(i_std), DIMENSION (kjpindex), INTENT (in)          :: njsc          !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless)
    REAL(r_std), DIMENSION (kjpindex,ngrnd), INTENT(IN)        :: smc           !! Volumetric Soil Moisture Content (m3/m3)
    REAL(r_std), DIMENSION (nscm), INTENT(IN)                  :: qz            !! Quartz Content (Soil Type Dependent) (0-1)               
    REAL(r_std), DIMENSION (kjpindex,ngrnd,nvm), INTENT(IN)    :: porosnet      !! Soil Porosity (0-1)
    REAL(r_std), DIMENSION (kjpindex,ngrnd), INTENT(IN)        :: sh2o          !! Unfrozen Soil Moisture Content; Frozen Soil Moisture = smc - sh2o
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm), INTENT(IN)     :: F_OM          !! Solid volume fraction of soil organic matter (does not include the air/porosity) [m3/m3]
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm), INTENT(IN)     :: F_MM          !! Solid volume fraction of mineral matter [m3/m3]
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm), INTENT(IN)     :: F_CM          !! Solid volume fraction of coarse matter [m3/m3]
    
    !! 0.2 Output variables
    REAL(r_std), DIMENSION (kjpindex,ngrnd,nvm), INTENT(OUT)   :: cnd           !! Soil Thermal Conductivity (W/m/k)
    
    !! 0.3 Modified variables
    
    !! 0.4 Local variables
    REAL(r_std), DIMENSION (kjpindex,ngrnd,nvm)                :: ake           !! Kerston Number (unitless)
    REAL(r_std), DIMENSION (kjpindex,ngrnd,nvm)                :: thksat        !! Saturated Thermal Conductivity (W/m/k)
    REAL(r_std), DIMENSION (kjpindex,ngrnd)                    :: thksat_moss   !! Moss saturated thermal conductivity (W/m/K)
    REAL(r_std), DIMENSION (kjpindex,ngrnd,nvm)                :: satratio      !! Degree of Saturation (0-1)
    REAL(r_std), DIMENSION (kjpindex,ngrnd,nvm)                :: xu            !! Unfrozen Volume For Saturation (0-1)
    REAL(r_std), DIMENSION (kjpindex,ngrnd,nvm)                :: xunfroz       !! Unfrozen Volume Fraction (0-1)
    REAL(r_std)                                                :: thko          !! Thermal Conductivity for Other Soil Components (W/m/k)
    REAL(r_std), DIMENSION (kjpindex)                          :: gammd         !! Dry Density (kg/m3)
    REAL(r_std), DIMENSION (kjpindex)                          :: thkdry_min    !! Dry Thermal Conductivity for mineral soil (W/m/k)
    REAL(r_std), DIMENSION (kjpindex,ngrnd,nvm)                :: thkdry        !! Dry Thermal Conductivity considering organic carbon (W/m/k)
    REAL(r_std), DIMENSION (kjpindex)                          :: thks_min      !! Thermal Conductivity for the Solids Combined (Quartz + Other) (W/m/k)
    REAL(r_std), DIMENSION (kjpindex,ngrnd,nvm)                :: thks          !! Thermal Conductivity considering organic carbon (W/m/k)
    REAL(r_std), DIMENSION (kjpindex,ngrnd)                    :: thkmoss       !! Thermal conductivity of moss (W/m/K)
    REAL(r_std), DIMENSION (kjpindex,ngrnd,nvm)                :: Phi_CM_MM, Phi_OM !! Volume fractions of coarse and mineral matter (including pores), and of organic matter (0-1) [m3/m3]
    REAL(r_std), DIMENSION (kjpindex,ngrnd)                    :: VWC
    INTEGER(i_std)                                             :: ji, jg, jst, jv
    REAL(r_std), PARAMETER                                     :: THKICE = 2.2 !! Ice Thermal Conductivity (W/m/k)
    REAL(r_std), PARAMETER                                     :: THKQTZ = 7.7 !! Thermal Conductivity for Quartz (W/m/k)
    REAL(r_std), PARAMETER                                     :: THKW = 0.57 !! Water Thermal Conductivity (W/m/k)
!_================================================================================================================================

    thksat(:,:,:)=0
    thkmoss(:,:)=0
    thksat_moss(:,:)=0
    !! 1. Dry and Saturated Thermal Conductivity.
    
    DO ji = 1,kjpindex
       jst = njsc(ji)
       
       !! Calculation of the repartition of organic matter and mineral+coarse matter including the pores (later used for thkdry)
       Phi_CM_MM(ji,:,:)=(F_CM(ji,:,:)+F_MM(ji,:,:))/(1. - mcs(jst)) ! Cumulated fractions of coarse and mineral matter
       Phi_OM(ji,:,:)=1.-Phi_CM_MM(ji,:,:)
       
       !! 1.1. Dry density (Kg/m3) and Dry thermal conductivity (W.M-1.K-1) - NO MODIFICATIONS
       gammd(ji) = (1. - mcs(jst))*2700.
       thkdry_min(ji) = (0.135* gammd(ji) + 64.7)/ (2700. - 0.947* gammd(ji))
       
       
       !! 1.2. thermal conductivity of "other" soil components
       IF (qz(jst) > 0.2) THEN
          thko = 2.0
       ELSEIF (qz(jst) <= 0.2) THEN
          thko = 3.0
       ENDIF
       
       !! 1.3. Thermal conductivity of solids
       thks_min(ji) = (THKQTZ ** qz(jst))* (thko ** (1. - qz(jst)))
    ENDDO
    
    DO jv = 1,nvm
       
       SELECT CASE (use_soilc_method)
       CASE (SOILC_METHOD_ARITHMETIC)
          DO jg = 1,ngrnd
             DO ji = 1,kjpindex
                !Add the different contributions of each element within the soil solids
                thks(ji,jg,jv) = (1. / (F_OM(ji,jg,jv) + F_MM(ji,jg,jv) + F_CM(ji,jg,jv))) * &
                     (F_OM(ji,jg,jv) * cond_solid_org + (F_MM(ji,jg,jv) + F_CM(ji,jg,jv)) * thks_min(ji))
                
                thkdry(ji,jg,jv) = Phi_OM(ji,jg,jv) * cond_dry_org + Phi_CM_MM(ji,jg,jv) * thkdry_min(ji)
             ENDDO
          ENDDO
       CASE (SOILC_METHOD_GEOMETRIC)
          DO jg = 1,ngrnd
             DO ji = 1,kjpindex
                ! use geometric mean rather than arithmetic mean (Decharme et al 2016), readapted with the correct fractions
                thks(ji,jg,jv) =((cond_solid_org**F_OM(ji,jg,jv)) * (thks_min(ji)**(F_MM(ji,jg,jv)+F_CM(ji,jg,jv)))) &
                     ** (1. / (F_OM(ji,jg,jv) + F_MM(ji,jg,jv) + F_CM(ji,jg,jv)))
                
                thkdry(ji,jg,jv) =(cond_dry_org**Phi_OM(ji,jg,jv)) * (thkdry_min(ji)**(Phi_CM_MM(ji,jg,jv)))
             ENDDO
          ENDDO
       CASE DEFAULT
          CALL ipslerr_p(3,'thermosoil_cond_pft_newdisc','Unsupported USE_SOILC_METHOD','','')
       END SELECT

       DO jg = 1,ngrnd
          !! 1.4. saturation ratio
          DO ji = 1,kjpindex
             satratio(ji,jg,jv) = smc(ji,jg) / porosnet(ji, jg, jv)
          ENDDO
          
          !! 1.5. Saturated Thermal Conductivity (thksat) - NO CHANGES
          DO ji = 1,kjpindex
             IF ( smc(ji,jg) > min_sechiba ) THEN
                xunfroz(ji,jg,jv) = sh2o(ji,jg) / smc(ji,jg)   ! Unfrozen Fraction (From i.e., 100%Liquid, to 0. (100% Frozen))
                xu(ji,jg,jv) = xunfroz(ji,jg,jv) * porosnet(ji, jg, jv)  ! Unfrozen volume for saturation (porosity*xunfroz)
                thksat(ji,jg,jv) = thks(ji,jg,jv) ** (1. - porosnet(ji, jg, jv)) &
                     * THKICE ** (porosnet(ji, jg, jv) - xu(ji,jg,jv)) * THKW ** xu(ji,jg,jv)
                IF (ok_moss) THEN
                   IF (zlt(jg) .LE. soil_moss_depth) THEN
                      ! Saturated conductivity for moss
                      thksat_moss(ji,jg) = thksat_moss_liq ** xu(ji,jg,jv) * thksat_moss_ice ** (porosnet(ji, jg, jv) - xu(ji,jg,jv))
                   ENDIF
                ENDIF
             ELSE
                ! this value will not be used since ake=0 for this case
                thksat(ji,jg,jv)=0 
             END IF
          ENDDO
       END DO ! DO jg = 1,ngrnd
       
       !! 2. Kersten Number (ake) - NO CHANGES
       DO jg = 1,ngrnd
          DO ji = 1,kjpindex
             jst = njsc(ji)
             IF ( (sh2o(ji,jg) + 0.0005) <  smc(ji,jg) ) THEN
                ! Frozen
                ake(ji,jg,jv) = satratio(ji,jg,jv)
             ELSE
                ! Unfrozen
                IF ( satratio(ji,jg,jv) >  0.1 ) THEN
                   IF ((jst < 4 .AND. soil_classif == 'usda') .OR. (jst == 1 .AND. soil_classif == 'zobler') )  THEN
                      ! Coarse
                      ake(ji,jg,jv) = 0.7 * LOG10 (satratio(ji,jg,jv)) + 1.0
                   ELSE
                      ! Fine
                      ake(ji,jg,jv) = LOG10 (satratio(ji,jg,jv)) + 1.0
                   ENDIF
                ELSEIF ( satratio(ji,jg,jv) >  0.05 .AND. satratio(ji,jg,jv) <=  0.1 ) THEN
                   IF ((jst < 4 .AND. soil_classif == 'usda') .OR. (jst == 1 .AND. soil_classif == 'zobler') )  THEN
                      ! Coarse
                      ake(ji,jg,jv) = 0.7 * LOG10 (satratio(ji,jg,jv)) + 1.0
                   ELSE
                      ! Fine
                      ake(ji,jg,jv) = 0.0
                   ENDIF
                ELSE
                   ake(ji,jg,jv) = 0.0  ! use k = kdry
                END IF
             END IF
          ENDDO
       END DO ! DO jg = 1,ngrnd
       
       !! 3. Thermal conductivity (cnd)
       
       DO jg = 1,ngrnd
          DO ji = 1,kjpindex
             cnd(ji,jg,jv) = ake(ji,jg,jv) * (thksat(ji,jg,jv) - thkdry(ji, jg, jv)) + thkdry(ji, jg, jv)
             
             ! Add moss thermal conductivity
             IF (ok_moss) THEN
                IF (zlt(jg) .LE. soil_moss_depth) THEN
                   thkmoss(ji,jg) = thkmoss_dry + &
                        tmc_layt(ji,jg) / mille / dlt(jg) *(thksat_moss(ji,jg) - thkmoss_dry) 
                   
                   SELECT CASE(use_moss_cond_method)
                   CASE (moss_combine)
                      cnd(ji,jg,jv) = moss_frac(jv) * thkmoss(ji,jg) * cnd(ji,jg,jg) *&
                           soil_moss_depth / (moss_layer_thickness * cnd(ji,jg,jv) +&
                           soil_moss_depth * thkmoss(ji,jg)) + (1 - moss_frac(jv)) * cnd(ji,jg,jv) 
                   CASE (moss_replace)
                      cnd(ji,jg,jv) = moss_frac(jv) * thkmoss(ji,jg) + (1 - moss_frac(jv)) * cnd(ji,jg,jv)
                   END SELECT
                ENDIF
             ENDIF
          ENDDO
       END DO
    END DO !DO jv = 1,nvm 
    
    IF (ok_moss) CALL xios_orchidee_send_field("thkmoss",thkmoss)


  END SUBROUTINE thermosoil_cond_pft_newdisc

!! ================================================================================================================================
!! SUBROUTINE   : thermosoil_read_soc_bd_file
!!
!>\BRIEF          
!!
!! DESCRIPTION : Reads three files of soil organic carbon to be used in thermix (insulating effect)
!! (files: SOC concentration, BD of fine matter and volume fraction of coarse matter)
!!                
!!
!! RECENT CHANGE(S) : None
!! 
!! MAIN OUTPUT VARIABLE(S):  calculates Phi_CM, fSOC and BD_FM from external data
!!                          
!! REFERENCE(S) :
!!
!! FLOWCHART    : None 
!! \n 
!_ ================================================================================================================================
    
  SUBROUTINE thermosoil_read_soc_bd_file(nbpt, lalo, neighbours, resolution, contfrac)

    !! 0. Variable and parameter declaration
    
    !! 0.1 Input variables
    
    INTEGER(i_std), INTENT(in)                    :: nbpt                  !! Number of points for which the data needs to be interpolated (unitless)
    REAL(r_std), INTENT(in)                       :: lalo(nbpt,2)          !! Vector of latitude and longitudes (degree)        
    INTEGER(i_std), INTENT(in)                    :: neighbours(nbpt,NbNeighb)!! Vector of neighbours for each grid point (1=N,2=E,3=S,4=W)  
    REAL(r_std), INTENT(in)                       :: resolution(nbpt,2)    !! The size of each grid cell in X and Y (km)
    REAL(r_std), INTENT(in)                       :: contfrac(nbpt)        !! Fraction of land in each grid cell (unitless)   
    
    !! 0.4 Local variables
    INTEGER(i_std)                                :: nbvmax                !! nbvmax for interpolation (unitless) 
    INTEGER(i_std)                                :: iml, jml, lml, tml    !! Indices
    INTEGER(i_std)                                :: fid, ib, ip, jp, fopt !! Indices
    INTEGER(i_std)                                :: ilf, ks               !! Indices
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:)   :: mask
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:,:) :: sub_index             !! Help variable to read file data and allocate memory
    INTEGER                                       :: ALLOC_ERR             !! Help varialbe to count allocation error
    INTEGER(i_std)                                :: ji, jg
    CHARACTER(LEN=80)                             :: filename              
    CHARACTER(LEN=30)                             :: callsign              !! Help variable to read file data and allocate memory
    CHARACTER(LEN=100)                            :: str                   !! Temporary string var 
    REAL(r_std)                                   :: totarea               !! Help variable to compute average SOC
    REAL(r_std), ALLOCATABLE, DIMENSION(:)        :: lat_lu, lon_lu        !! Latitudes and longitudes read from input file
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)      :: lat_rel, lon_rel      !! Help variable to read file data and allocate memory
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)      :: mask_lu               !! Help variable to read file data and allocate memory
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:)    :: refsoc_file           !! Help variable to read file data and allocate memory
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:)    :: bulk_file             !! Help variable to read file data and allocate memory
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:)    :: CM_file               !! Help variable to read file data and allocate memory
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)      :: sub_area              !! Help variable to read file data and allocate memory   
    REAL(r_std),DIMENSION(nbpt,ngrnd)             :: BD_SOC                !! Bulk density of SOC in the full volume, from the map (kg/m3)
    REAL(r_std),DIMENSION(nbpt,ngrnd)             :: BD_FM_tot             !! Mass of fine matter per volume unit (in the full volume) (kg/m3)
    LOGICAL                                       :: ok_interpol           !! Optional return of aggregate_2d
!_
!================================================================================================================================

    !! 1. Open file and allocate memory
    
    ! Open file with SOC map
    
    !Config Key   = SOC_FILE
    !Config Desc  = File with SOC concentrations (bulk density of SOC in the total volume - kg/m3)
    !Config If    = OK_SOIL_CARBON_DISCRETIZATION, USE_REFSOC, SOIL_CTEMPDIFF, USE_THERMO_NEWDISC
    !Config Def   = SOC.nc
    !Config Help  = 
    !Config Units = [FILE]
    filename = 'SOC.nc'
    CALL getin_p('SOC_FILE',filename)
    
    ! Read data from file
    IF (is_root_prc) CALL flininfo(filename, iml, jml, lml, tml, fid)
    CALL bcast(iml)
    CALL bcast(jml)
    CALL bcast(lml)
    CALL bcast(tml)
    
    IF (lml .NE. ngrnd) THEN
       WRITE(str, *) 'ngrnd=', ngrnd, ', depth found in file=', lml
       CALL ipslerr_p(3, 'thermosoil_read_soc_bd_file',    &
            'depth from the file must be the same as ngrnd', &
            str, &
            filename )
    ENDIF
    
    ALLOCATE(lon_lu(iml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'thermosoil_read_soc_bd_file','Problem in allocation of variable lon_lu','','')
    
    ALLOCATE(lat_lu(jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'thermosoil_read_soc_bd_file','Problem in allocation of variable lat_lu','','')
    
    ALLOCATE(mask_lu(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'thermosoil_read_soc_bd_file','Pb in allocation for mask_lu','','')
    
    ALLOCATE(refsoc_file(iml,jml,lml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'thermosoil_read_soc_bd_file','Pb in allocation for refsoc_file','','')
    
    IF (is_root_prc) THEN
       CALL flinget(fid, 'longitude', iml, 0, 0, 0, 1, 1, lon_lu)
       CALL flinget(fid, 'latitude', jml, 0, 0, 0, 1, 1, lat_lu)
       CALL flinget(fid, 'mask', iml, jml, 0, 0, 1, 1, mask_lu)
       CALL flinget(fid, 'Soil_Organic_Carbon', iml, jml, lml, tml, 1, 1, refsoc_file)
       
       CALL flinclo(fid)
    ENDIF
    
    CALL bcast(lon_lu)
    CALL bcast(lat_lu)
    CALL bcast(mask_lu)
    CALL bcast(refsoc_file)
    
    ALLOCATE(lon_rel(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'thermosoil_read_soc_bd_file','Pb in allocation for lon_rel','','')
    
    ALLOCATE(lat_rel(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'thermosoil_read_soc_bd_file','Pb in allocation for lat_rel','','')
    
    ALLOCATE(mask(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'thermosoil_read_soc_bd_file','Problem in allocation of variable mask','','')
    
    DO jp=1,jml
       lon_rel(:,jp) = lon_lu(:)
    ENDDO
    DO ip=1,iml
       lat_rel(ip,:) = lat_lu(:)
    ENDDO
    
    mask(:,:) = zero
    WHERE (mask_lu(:,:) > zero )
       mask(:,:) = un
    ENDWHERE
    
    ! Set nbvmax to 200 for interpolation
    ! This number is the dimension of the variables in which we store 
    ! the list of points of the source grid which fit into one grid box of the
    ! target. 
    nbvmax = 16
    callsign = 'Soil Organic Carbon'
    
    ! Start interpolation
    ok_interpol=.FALSE.
    DO WHILE ( .NOT. ok_interpol )
       WRITE(numout,*) "Projection arrays for ",callsign," : "
       WRITE(numout,*) "nbvmax = ",nbvmax
       
       ALLOCATE(sub_area(nbpt,nbvmax), STAT=ALLOC_ERR)
       IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'thermosoil_read_soc_bd_file','Pb in allocation for sub_area','','')
       sub_area(:,:)=zero
       
       ALLOCATE(sub_index(nbpt,nbvmax,2), STAT=ALLOC_ERR)
       IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'thermosoil_read_soc_bd_file','Pb in allocation for sub_index','','')
       sub_index(:,:,:)=0
       
       CALL aggregate_p(nbpt, lalo, neighbours, resolution, contfrac, &
            iml, jml, lon_rel, lat_rel, mask, callsign, &
            nbvmax, sub_index, sub_area, ok_interpol)
       
       IF ( .NOT. ok_interpol ) THEN
          DEALLOCATE(sub_area)
          DEALLOCATE(sub_index)
          nbvmax = nbvmax * 2
       ENDIF
    ENDDO
    
    ! Compute the average
    BD_SOC(:,:) = zero
    DO ib = 1, nbpt
       fopt = COUNT(sub_area(ib,:) > zero)
       IF ( fopt > 0 ) THEN
          totarea = zero
          DO ilf = 1, fopt
             ip = sub_index(ib,ilf,1)
             jp = sub_index(ib,ilf,2)
             BD_SOC(ib,:) = BD_SOC(ib,:) + refsoc_file(ip,jp,:) * sub_area(ib,ilf)
             totarea = totarea + sub_area(ib,ilf)
          ENDDO
          ! Normalize
          BD_SOC(ib,:) = BD_SOC(ib,:)/totarea
       ELSE
          ! Set defalut value for points where the interpolation fail
          WRITE(numout,*) 'On point ', ib, ' no points were found for interpolation data. Mean value is used.'
          WRITE(numout,*) 'Location : ', lalo(ib,2), lalo(ib,1)
          BD_SOC(ib,:) = 0.
       ENDIF
    ENDDO
    
    DEALLOCATE (lat_lu)
    DEALLOCATE (lat_rel)
    DEALLOCATE (lon_lu)
    DEALLOCATE (lon_rel)
    DEALLOCATE (mask_lu)
    DEALLOCATE (mask)
    DEALLOCATE (refsoc_file)
    DEALLOCATE (sub_area)
    DEALLOCATE (sub_index)
    
    ! Open file with BD map
    
    !Config Key   = BD_FILE
    !Config Desc  = File with BD of fine matter in the total volume (kg/m3)
    !Config If    = OK_SOIL_CARBON_DISCRETIZATION, USE_REFSOC, SOIL_CTEMPDIFF, USE_THERMO_NEWDISC
    !Config Def   = BD.nc
    !Config Help  = 
    !Config Units = [FILE]
    filename = 'BD.nc'
    CALL getin_p('BD_FILE',filename)
    
    ! Read data from file
    IF (is_root_prc) CALL flininfo(filename, iml, jml, lml, tml, fid)
    CALL bcast(iml)
    CALL bcast(jml)
    CALL bcast(lml)
    CALL bcast(tml)
    
    IF (lml .NE. ngrnd) THEN
       WRITE(str, *) 'ngrnd=', ngrnd, ', depth found in file=', lml
       CALL ipslerr_p(3, 'thermosoil_read_soc_bd_file',    &
            'depth from the file must be the same as ngrnd', &
            str, &
            filename )
    ENDIF
 
    ALLOCATE(lon_lu(iml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'thermosoil_read_soc_bd_file','Problem in allocation of variable lon_lu','','')
    
    ALLOCATE(lat_lu(jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'thermosoil_read_soc_bd_file','Problem in allocation of variable lat_lu','','')
    
    ALLOCATE(mask_lu(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'thermosoil_read_soc_bd_file','Pb in allocation for mask_lu','','')
    
    ALLOCATE(bulk_file(iml,jml,lml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'thermosoil_read_soc_bd_file','Pb in allocation for bulk_file','','')
    
    IF (is_root_prc) THEN
       CALL flinget(fid, 'longitude', iml, 0, 0, 0, 1, 1, lon_lu)
       CALL flinget(fid, 'latitude', jml, 0, 0, 0, 1, 1, lat_lu)
       CALL flinget(fid, 'mask', iml, jml, 0, 0, 1, 1, mask_lu)
       CALL flinget(fid, 'fine_matter_density', iml, jml, lml, tml, 1, 1, bulk_file)
       
       CALL flinclo(fid)
    ENDIF
    
    CALL bcast(lon_lu)
    CALL bcast(lat_lu)
    CALL bcast(mask_lu)
    CALL bcast(bulk_file)
    
    ALLOCATE(lon_rel(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'thermosoil_read_soc_bd_file','Pb in allocation for lon_rel','','')
    
    ALLOCATE(lat_rel(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'thermosoil_read_soc_bd_file','Pb in allocation for lat_rel','','')
    
    ALLOCATE(mask(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'thermosoil_read_soc_bd_file','Problem in allocation of variable mask','','')
    
    DO jp=1,jml
       lon_rel(:,jp) = lon_lu(:)
    ENDDO
    DO ip=1,iml
       lat_rel(ip,:) = lat_lu(:)
    ENDDO
    
    mask(:,:) = zero
    WHERE (mask_lu(:,:) > zero )
       mask(:,:) = un
    ENDWHERE
    
    ! Set nbvmax to 200 for interpolation
    ! This number is the dimension of the variables in which we store 
    ! the list of points of the source grid which fit into one grid box of the
    ! target. 
    nbvmax = 16
    callsign = 'Bulk Density Tot'
    
    ! Start interpolation
    ok_interpol=.FALSE.
    DO WHILE ( .NOT. ok_interpol )
       WRITE(numout,*) "Projection arrays for ",callsign," : "
       WRITE(numout,*) "nbvmax = ",nbvmax
       
       ALLOCATE(sub_area(nbpt,nbvmax), STAT=ALLOC_ERR)
       IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'thermosoil_read_soc_bd_file','Pb in allocation for sub_area','','')
       sub_area(:,:)=zero
       
       ALLOCATE(sub_index(nbpt,nbvmax,2), STAT=ALLOC_ERR)
       IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'thermosoil_read_soc_bd_file','Pb in allocation for sub_index','','')
       sub_index(:,:,:)=0
       
       CALL aggregate_p(nbpt, lalo, neighbours, resolution, contfrac, &
            iml, jml, lon_rel, lat_rel, mask, callsign, &
            nbvmax, sub_index, sub_area, ok_interpol)
       
       IF ( .NOT. ok_interpol ) THEN
          DEALLOCATE(sub_area)
          DEALLOCATE(sub_index)
          nbvmax = nbvmax * 2
       ENDIF
    ENDDO
    
    ! Compute the average
    BD_FM_tot(:,:) = zero
    DO ib = 1, nbpt
       fopt = COUNT(sub_area(ib,:) > zero)
       IF ( fopt > 0 ) THEN
          totarea = zero
          DO ilf = 1, fopt
             ip = sub_index(ib,ilf,1)
             jp = sub_index(ib,ilf,2)
             BD_FM_tot(ib,:) = BD_FM_tot(ib,:) + bulk_file(ip,jp,:) * sub_area(ib,ilf)
             totarea = totarea + sub_area(ib,ilf)
          ENDDO
          ! Normalize
          BD_FM_tot(ib,:) = BD_FM_tot(ib,:)/totarea
       ELSE
          ! Set defalut value for points where the interpolation fail
          WRITE(numout,*) 'On point ', ib, ' no points were found for interpolation data. Mean value is used.'
          WRITE(numout,*) 'Location : ', lalo(ib,2), lalo(ib,1)
          BD_FM_tot(ib,:) = 1440.
       ENDIF
    ENDDO
    
    DEALLOCATE (lat_lu)
    DEALLOCATE (lat_rel)
    DEALLOCATE (lon_lu)
    DEALLOCATE (lon_rel)
    DEALLOCATE (mask_lu)
    DEALLOCATE (mask)
    DEALLOCATE (bulk_file)
    DEALLOCATE (sub_area)
    DEALLOCATE (sub_index)
    
    ! Open file with CM map
    
    !Config Key   = CM_FILE
    !Config Desc  = File with volume fraction of coarse matter (m3/m3)
    !Config If    = OK_SOIL_CARBON_DISCRETIZATION, USE_REFSOC, SOIL_CTEMPDIFF, USE_THERMO_NEWDISC
    !Config Def   = CM.nc
    !Config Help  = 
    !Config Units = [FILE]
    filename = 'CM.nc'
    CALL getin_p('CM_FILE',filename)
    
    ! Read data from file
    IF (is_root_prc) CALL flininfo(filename, iml, jml, lml, tml, fid)
    CALL bcast(iml)
    CALL bcast(jml)
    CALL bcast(lml)
    CALL bcast(tml)
    
    IF (lml .NE. ngrnd) THEN
       WRITE(str, *) 'ngrnd=', ngrnd, ', depth found in file=', lml
       CALL ipslerr_p(3, 'thermosoil_read_soc_bd_file',    &
            'depth from the file must be the same as ngrnd', &
            str, &
            filename )
    ENDIF
    
    ALLOCATE(lon_lu(iml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'thermosoil_read_soc_bd_file','Problem in allocation of variable lon_lu','','')
    
    ALLOCATE(lat_lu(jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'thermosoil_read_soc_bd_file','Problem in allocation of variable lat_lu','','')
    
    ALLOCATE(mask_lu(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'thermosoil_read_soc_bd_file','Pb in allocation for mask_lu','','')
    
    ALLOCATE(CM_file(iml,jml,lml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'thermosoil_read_soc_bd_file','Pb in allocation for CM_file','','')
    
    IF (is_root_prc) THEN
       CALL flinget(fid, 'longitude', iml, 0, 0, 0, 1, 1, lon_lu)
       CALL flinget(fid, 'latitude', jml, 0, 0, 0, 1, 1, lat_lu)
       CALL flinget(fid, 'mask', iml, jml, 0, 0, 1, 1, mask_lu)
       CALL flinget(fid, 'coarse_matter', iml, jml, lml, tml, 1, 1, CM_file)
       CALL flinclo(fid)
    ENDIF
    
    CALL bcast(lon_lu)
    CALL bcast(lat_lu)
    CALL bcast(mask_lu)
    CALL bcast(CM_file)
    
    ALLOCATE(lon_rel(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'thermosoil_read_soc_bd_file','Pb in allocation for lon_rel','','')
    
    ALLOCATE(lat_rel(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'thermosoil_read_soc_bd_file','Pb in allocation for lat_rel','','')
    
    ALLOCATE(mask(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'thermosoil_read_soc_bd_file','Problem in allocation of variable mask','','')
    
    DO jp=1,jml
       lon_rel(:,jp) = lon_lu(:)
    ENDDO
    DO ip=1,iml
       lat_rel(ip,:) = lat_lu(:)
    ENDDO
    
    mask(:,:) = zero
    WHERE (mask_lu(:,:) > zero )
       mask(:,:) = un
    ENDWHERE
    
    ! Set nbvmax to 200 for interpolation
    ! This number is the dimension of the variables in which we store 
    ! the list of points of the source grid which fit into one grid box of the
    ! target. 
    nbvmax = 16
    callsign = 'Coarse Matter'
    
    ! Start interpolation
    ok_interpol=.FALSE.
    DO WHILE ( .NOT. ok_interpol )
       WRITE(numout,*) "Projection arrays for ",callsign," : "
       WRITE(numout,*) "nbvmax = ",nbvmax
       
       ALLOCATE(sub_area(nbpt,nbvmax), STAT=ALLOC_ERR)
       IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'thermosoil_read_soc_bd_file','Pb in allocation for sub_area','','')
       sub_area(:,:)=zero
       
       ALLOCATE(sub_index(nbpt,nbvmax,2), STAT=ALLOC_ERR)
       IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'thermosoil_read_soc_bd_file','Pb in allocation for sub_index','','')
       sub_index(:,:,:)=0
       
       CALL aggregate_p(nbpt, lalo, neighbours, resolution, contfrac, &
            iml, jml, lon_rel, lat_rel, mask, callsign, &
            nbvmax, sub_index, sub_area, ok_interpol)
       
       IF ( .NOT. ok_interpol ) THEN
          DEALLOCATE(sub_area)
          DEALLOCATE(sub_index)
          nbvmax = nbvmax * 2
       ENDIF
    ENDDO
    
    ! Compute the average
    Phi_CM(:,:) = zero
    DO ib = 1, nbpt
       fopt = COUNT(sub_area(ib,:) > zero)
       IF ( fopt > 0 ) THEN
          totarea = zero
          DO ilf = 1, fopt
             ip = sub_index(ib,ilf,1)
             jp = sub_index(ib,ilf,2)
             Phi_CM(ib,:) = Phi_CM(ib,:) + CM_file(ip,jp,:) * sub_area(ib,ilf)
             totarea = totarea + sub_area(ib,ilf)
          ENDDO
          ! Normalize
          Phi_CM(ib,:) = Phi_CM(ib,:)/totarea
       ELSE
          ! Set defalut value for points where the interpolation fail
          WRITE(numout,*) 'On point ', ib, ' no points were found for interpolation data. Mean value is used.'
          WRITE(numout,*) 'Location : ', lalo(ib,2), lalo(ib,1)
          Phi_CM(ib,:) = 0.
       ENDIF
    ENDDO
    
    DEALLOCATE (lat_lu)
    DEALLOCATE (lat_rel)
    DEALLOCATE (lon_lu)
    DEALLOCATE (lon_rel)
    DEALLOCATE (mask_lu)
    DEALLOCATE (mask)
    DEALLOCATE (CM_file)
    DEALLOCATE (sub_area)
    DEALLOCATE (sub_index)
    
    ! Calculate fSOC and BD_FM
    DO jg = 1,ngrnd
       DO ji = 1,nbpt
          IF (Phi_CM(ji,jg) .LT. 1.) THEN ! check that we don't divide by 0
             BD_FM(ji,jg) = BD_FM_tot(ji,jg) / (1.-Phi_CM(ji,jg))
          ELSE
             BD_FM(ji,jg) = 1620.
          ENDIF
          IF (BD_FM_tot(ji,jg) .GT. 0.) THEN ! check that we don't divide by 0
             fSOC(ji,jg) = MIN(BD_SOC(ji,jg) / BD_FM_tot(ji,jg), 1.)
          ELSE
             fSOC(ji,jg) = 0.
          ENDIF
       ENDDO
    ENDDO
    
  END SUBROUTINE thermosoil_read_soc_bd_file
  

END MODULE thermosoil
