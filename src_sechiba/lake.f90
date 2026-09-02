! ===================================================================================================\n
! MODULE        : lake
!
! CONTACT       : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE       : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
! BRIEF        This module computes the lakes energy budget
!!
!! DESCRIPTION : This module calculates the lake surface energy budget, surface temperature and fluxes 
!! as well as the thermal processes within lakes, i.e., vertical temperature structure and mixing conditions. 
!! Snow, ice and sediment profile temperatures can also be computed in specific layers. One sediment layer 
!! may be activated with keyword ok_botsed_lake set to True. A prescribed number of effective lakes may be 
!! accounted within the grid cell (3 is the default number), on which a specific energy budget is done. 
!! Output variables are lake surface temperature, radiative, momentum and turbulent fluxes, bottom temperature 
!! and shape factor of the thermocline. Ice and snow temperatures at the top and bottom of each layer are 
!! also written in the output file. 
!!
!! This module is based on the FLake lake model developed by Mironov (2008)and implemented in ORCHIDEE 
!! by Anthony Bernus during his PhD in 2019-2022 (Bernus and Ottlé, 2022). The code was revised in 2023 
!! by Zacharie Titus (Master student) and Catherine Ottlé to include an ice cover fraction based on 
!! Garnaud et al., 2022 which allowed to improve the lake surface albedo model (Titus et al., 2024).
!!
!!
!! RECENT CHANGE(S) : None
!!
!! REFERENCE(S) :
!!	Mironov, D., Ritter, B., Schulz, J. P., Buchhold, M., Lange, M., & MacHulskaya, E.: 
!!      Parameterisation of sea and lake ice in numerical weather prediction models of the German Weather Service. Tellus A: 
!!      Dynamic Meteorology and Oceanography, 64(1), 17330, doi: 10.3402/tellusa.v64i0.17330, 2012.
!!
!!	Garnaud, C., MacKay, M., & Fortin, V.: A One‐Dimensional Lake Model in ECCC's Land Surface Prediction System. 
!!      Journal of Advances in Modeling Earth Systems, 14(2), e2021MS002861, doi:10.1029/2021MS002861, 2022.
!!
!!	Bernus, A. and Ottlé, C.:  Modeling subgrid lake energy balance in ORCHIDEE terrestrial scheme using the FLake lake model, 
!!      Geosci. Model Dev., 15, 4275–4295, doi:10.5194/gmd-15-4275-2022, 2022. 
!!
!!	Titus, Z., Cuynet, A., Salmon, E., and Ottlé, C.: Brief communication: Improving lake ice modeling in 
!!      ORCHIDEE-FLake model using MODIS albedo data, EGUsphere [preprint], https://doi.org/10.5194/egusphere-2024-2907, 2024. 
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_sechiba/lake.f90 $
!! $Date: 2025-12-12 08:47:06 +0100 (ven. 12 déc. 2025) $
!! $Revision: 9253 $
!_ ===============================================================================================\n
MODULE lake

  USE ioipsl
  USE xios_orchidee
  USE constantes
  USE time, ONLY : one_day, dt_sechiba, julian_diff
  USE constantes_soil
  USE pft_parameters
  USE sechiba_io_p
  USE grid
  USE qsat_moisture
  USE interpweight
  USE function_library,  ONLY:  get_printlev

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: lake_main, lake_initialize, lake_finalize, lake_clear, lake_xios_initialize

  !
  ! Variables used inside lake module : declaration and initialisation
  !
  INTEGER(i_std), SAVE                             :: printlev_loc   !! local printlev for this module
!$OMP THREADPRIVATE(printlev_loc)
  INTEGER(i_std), SAVE                             :: lake_year      !! year for land cover and lake fraction update
!$OMP THREADPRIVATE(lake_year)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: T_mnw_flk      !! Mean temperature of the water column (in Kelvin)
!$OMP THREADPRIVATE(T_mnw_flk)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: T_snow_flk     !! Snow temperature at the atmosphere interface (in Kelvin)
!$OMP THREADPRIVATE(T_snow_flk)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: T_ice_flk      !! Ice temperature (in Kelvin)
!$OMP THREADPRIVATE(T_ice_flk)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: T_wML_flk      !! Mixed layer temperature (in Kelvin)
!$OMP THREADPRIVATE(T_wML_flk)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: T_bot_flk      !! Temperature at the lake bottom (in Kelvin)
!$OMP THREADPRIVATE(T_bot_flk)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: T_B1_flk       !! Temperature at the bottom of the sediment layer (in Kelvin)
!$OMP THREADPRIVATE(T_B1_flk)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: C_T_flk        !! Shape coefficient thermocline (unitless)
!$OMP THREADPRIVATE(C_T_flk)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: h_snow_flk     !! Snow layer thickness (in meter)            
!$OMP THREADPRIVATE(h_snow_flk)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: h_ice_flk      !! Ice layer thickness (in meter)      
!$OMP THREADPRIVATE(h_ice_flk)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: h_ML_flk       !! Mixed layer thickness (in meter)
!$OMP THREADPRIVATE(h_ML_flk)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: H_B1_flk       !! Sediment layer thickness (in meter)
!$OMP THREADPRIVATE(H_B1_flk)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: regime_flux_lake !! Flux regime (2 regimes of convection are modeled: convection driven by surface fluxes or convection driven by surface cooling
!$OMP THREADPRIVATE(regime_flux_lake)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: qs_out         !! Saturation specific humidity at surface temperature [-] (kg kg^{-1})
!$OMP THREADPRIVATE(qs_out)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: depth_w        !! Lake depth (in meter)
!$OMP THREADPRIVATE(depth_w)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: light_ext_water_lake !! water optical thickness (in meter-1)
!$OMP THREADPRIVATE(light_ext_water_lake)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: fetch          !!  wind fetch (in meter) assessed from lake area (diameter of an equivalent disc of same area)
!$OMP THREADPRIVATE(fetch)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: fraclake       !! Vector of lake surface fraction for each grid point and for each lake type
!$OMP THREADPRIVATE(fraclake)
  LOGICAL, SAVE                                    :: xios_interpolation_lakefrac  !! Flag for interpolating using XIOS 
!$OMP THREADPRIVATE(xios_interpolation_lakefrac)
  LOGICAL, SAVE                                    :: xios_interpolation_lakeparam !! Flag for interpolating using XIOS 
!$OMP THREADPRIVATE(xios_interpolation_lakeparam)

CONTAINS


  !! =============================================================================================================================
  !! SUBROUTINE:    lake_xios_initialize
  !!
  !! DESCRIPTION:	  Initialize xios dependant defintions before closing context defintion
  !! 
  !! ==============================================================================================================================
  SUBROUTINE lake_xios_initialize

    CHARACTER(LEN=255) :: filename
    CHARACTER(LEN=255) :: name

    !! Initialize local printlev
    printlev_loc=get_printlev('lake')

    xios_interpolation_lakefrac = xios_interpolation
    CALL getin_p('XIOS_INTERPOLATION_LAKEFRAC',xios_interpolation_lakefrac)

    xios_interpolation_lakeparam = xios_interpolation
    CALL getin_p('XIOS_INTERPOLATION_LAKEPARAM',xios_interpolation_lakeparam)

    !! 1. Treat the file for lake fraction

    ! Read the input file name from run.def
    filename = 'lake_frac.nc'
    CALL getin_p('LAKE_FRAC_FILE',filename) 

    ! Check if the file will be read by XIOS, by IOIPSL or not at all
    IF (xios_interpolation_lakefrac .AND. ok_lake_energy .AND. restname_in=='NONE') THEN
       IF (printlev_loc>=2) WRITE (numout,*) 'The lake fraction file will be read later by XIOS'

       ! Define the name of the input file to be red later. Remove suffix .nc from filename.
       name = filename(1:LEN_TRIM(filename)-3)
       CALL xios_orchidee_set_file_attr("lake_frac_file",name=name)
    ELSE
       IF (ok_lake_energy .AND. restname_in=='NONE') THEN
          IF (printlev_loc>=2) WRITE (numout,*) 'The lake fraction file will be read later by IOIPSL'
       ELSE
          IF (printlev_loc>=2) WRITE (numout,*) 'The lake fraction file will not be read'
       END IF

       ! Deactivate file and related fields
       CALL xios_orchidee_set_file_attr("lake_frac_file",enabled=.FALSE.)
       CALL xios_orchidee_set_field_attr("fraclake",enabled=.FALSE.)
    END IF


    !! 2. Treat the file for lake parameters

    ! Read the input file name from run.def
    filename = 'lake_param.nc'
    CALL getin_p('LAKE_PARAM_FILE',filename) 

    ! Check if the file will be read by XIOS, by IOIPSL or not at all
    IF (xios_interpolation_lakeparam .AND. ok_lake_energy .AND. restname_in=='NONE') THEN
       IF (printlev_loc>=2) WRITE (numout,*) 'The lake parameter file will be read later by XIOS'

       ! Define filename to be red later. Remove suffix .nc from filename.
       name = filename(1:LEN_TRIM(FILENAME)-3)
       CALL xios_orchidee_set_file_attr("lake_param_file",name=name)
    ELSE
       IF (ok_lake_energy .AND. restname_in=='NONE') THEN
          IF (printlev_loc>=2) WRITE (numout,*) 'The lake parameter file will be read later by IOIPSL'
       ELSE
          IF (printlev_loc>=2) WRITE (numout,*) 'The lake parameter file will not be read'
       END IF

       ! Deactivate file and related fields
       CALL xios_orchidee_set_file_attr("lake_param_file",enabled=.FALSE.)
       CALL xios_orchidee_set_field_attr("depth",enabled=.FALSE.)
       CALL xios_orchidee_set_field_attr("Cext",enabled=.FALSE.)
       CALL xios_orchidee_set_field_attr("fetch",enabled=.FALSE.)
    END IF


  END SUBROUTINE lake_xios_initialize



  !! ================================================================================================================================
  !! SUBROUTINE 	: lake_initialize
  !!
  !!
  !! DESCRIPTION : Allocate module variables and read them from restart file or initialize them by reading file or with default values
  !!
  !! MAIN OUTPUT VARIABLE(S) : 
  !!
  !! REFERENCE(S) : 
  !!
  !! FLOWCHART    : None
  !! 
  !! ================================================================================================================================

  SUBROUTINE lake_initialize ( kjit,          kjpindex,      index,         rest_id,                 &
                               lalo,          neighbours,    resolution,    contfrac,                &
                               fraclake_out,  temp_surf_lake)

    !! 0. Variable and parameter declaration
    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                                :: kjit             !! Time step number 
    INTEGER(i_std), INTENT(in)                                :: kjpindex         !! Domain size
    INTEGER(i_std), DIMENSION (kjpindex), INTENT (in)         :: index            !! Indices of the points on the map
    INTEGER(i_std), INTENT (in)                               :: rest_id          !! Restart file identifier
    REAL(r_std),DIMENSION (kjpindex,2), INTENT (in)           :: lalo             !! Geogr. coordinates (latitude,longitude) (degrees)
    INTEGER(i_std), DIMENSION (kjpindex,NbNeighb), INTENT(in) :: neighbours       !! Vector of neighbours for each grid point
    REAL(r_std), DIMENSION (kjpindex,2), INTENT(in)           :: resolution       !! size in x and y of the grid (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)             :: contfrac         !! Fraction of continent in the grid (unitless)

    !! 0.2 Output variables
    REAL(r_std), DIMENSION (kjpindex,nlake), INTENT(out)      :: fraclake_out    !! Vector of lake surface fraction for each grid point and for each lake type
    REAL(r_std),DIMENSION (kjpindex, nlake), INTENT(out)      :: temp_surf_lake  !! Lake surface temperature (K)

    !! 0.4 Local variables
    INTEGER(i_std) :: ij, ilake
    INTEGER(i_std) :: ier
    ! ================================================================================================================================

    !! 1. Initialize local printlev.
    !  This has already been done in lake_xios_initialize but is done again for the 
    !  case that subroutine was not called (case without XIOS)
    printlev_loc=get_printlev('lake')


    !! 2. Allocate module variables
    ALLOCATE (T_B1_flk(kjpindex, nlake),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'lake_initialize','Problem in allocate of variable T_B1','','')

    ALLOCATE (T_mnw_flk(kjpindex, nlake),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'lake_initialize','Problem in allocate of variable T_mnw','','')

    ALLOCATE (T_snow_flk(kjpindex, nlake),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'lake_initialize','Problem in allocate of variable T_snow','','')

    ALLOCATE (T_ice_flk(kjpindex, nlake),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'lake_initialize','Problem in allocate of variable T_ice','','')

    ALLOCATE (T_wML_flk(kjpindex, nlake),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'lake_initialize','Problem in allocate of variable T_wML','','')

    ALLOCATE (T_bot_flk(kjpindex, nlake),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'lake_initialize','Problem in allocate of variable T_bot','','')

    ALLOCATE (H_B1_flk(kjpindex, nlake),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'lake_initialize','Problem in allocate of variable H_B1','','')

    ALLOCATE (h_ML_flk(kjpindex, nlake),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'lake_initialize','Problem in allocate of variable h_ML','','')

    ALLOCATE (h_ice_flk(kjpindex, nlake),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'lake_initialize','Problem in allocate of variable h_ice','','')

    ALLOCATE (h_snow_flk(kjpindex, nlake),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'lake_initialize','Problem in allocate of variable h_snow','','')

    ALLOCATE (C_T_flk(kjpindex, nlake),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'lake_initialize','Problem in allocate of variable C_T','','')

    ALLOCATE (depth_w(kjpindex, nlake),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'lake_initialize','Problem in allocate of variable depth_w','','')

    ALLOCATE (fraclake(kjpindex, nlake),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'lake_initialize','Problem in allocate of variable fraclake','','')

    ALLOCATE (light_ext_water_lake(kjpindex, nlake),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'lake_initialize','Problem in allocate of variable light_ext_water_lake','','')

    ALLOCATE (fetch(kjpindex, nlake),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'lake_initialize','Problem in allocate of variable fetch','','')

    ALLOCATE (regime_flux_lake(kjpindex, nlake),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'lake_initialize','Problem in allocate of variable C_T','','')

    ALLOCATE (qs_out(kjpindex, nlake),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'lake_initialize','Problem in allocate of variable C_T','','')


    !! 3. Read module variables from restart file. If they are not available, initialize by reading files or using default values.

    ! Read fraclake, depth_w, light_ext_water_lake and fetch from restart file
    CALL restget_p (rest_id, 'fraclake', nbp_glo, nlake, 1, kjit, .TRUE., fraclake, "gather", nbp_glo, index_g)
    CALL restget_p (rest_id, 'depth_w', nbp_glo, nlake, 1, kjit, .TRUE., depth_w, "gather", nbp_glo, index_g)
    CALL restget_p (rest_id, 'light_ext_water_lake', nbp_glo, nlake, 1, kjit, .TRUE., light_ext_water_lake, "gather", nbp_glo, index_g)
    CALL restget_p (rest_id, 'fetch', nbp_glo, nlake, 1, kjit, .TRUE., fetch, "gather", nbp_glo, index_g)

    ! Check if all these 4 variables were found in the restart file, if not read them from file
    IF ( ALL(fraclake(:,:) .EQ. val_exp) .OR. &
         ALL(depth_w(:,:) .EQ. val_exp) .OR. &
         ALL(light_ext_water_lake(:,:) .EQ. val_exp) .OR. &
         ALL(fetch(:,:) .EQ. val_exp) ) THEN

       ! One or several variables were missing. Initialize them by reading from file
       ! Read lake fractions, fraclake from file. This map contains the actual lake fractions if the n lakes. 
       CALL lake_read_frac( kjpindex,  lalo, neighbours,  resolution, contfrac, lake_year, fraclake) 
       
       ! Read parameters depth_w, light_ext_water_lake and fetch from file
       CALL lake_read_param(kjpindex, lalo, neighbours, resolution, contfrac, lake_year, depth_w, light_ext_water_lake, fetch)
       
       ! Ajust fraclake and depth_w and check if there is a lake in the grid cell
       CALL lake_adjust_frac_and_depth(kjpindex, depth_w, fraclake)
    END IF

    ! Copy fraclake in module saved variable
    DO ij = 1, kjpindex
       DO ilake = 1, nlake
          fraclake_out(ij, ilake) = fraclake(ij, ilake)
       END DO
    END DO


    ! Read or initialize temp_surf_lake
    CALL restget_p (rest_id, 'temp_surf_lake', nbp_glo, nlake, 1, kjit, .TRUE., temp_surf_lake, "gather", nbp_glo, index_g)
    !Config Key   = TEMP_SURF_LAKE
    !Config Desc  = Initiale value of temp_surf_lake used if not found in restart file
    !Config If    = OK_LAKE_ENERGY
    !Config Def   = 273.15
    !Config Help  = The initial value of lake surface temperature if its value is not found
    !Config         in the restart file. This will only be used if the model is 
    !Config         started without a restart file.
    !Config Units = Kelvin [K]
    CALL setvar_p (temp_surf_lake, val_exp,'TEMP_SURF_LAKE', tp_00)


    ! Read or initialize T_snow_flk
    CALL restget_p (rest_id, 'T_snow_flk', nbp_glo, nlake, 1, kjit, .TRUE., T_snow_flk, "gather", nbp_glo, index_g)
    !Config Key   = T_SNOW_FLK
    !Config Desc  = Initiale value of T_snow_flk used if not found in restart file
    !Config If    = OK_LAKE_ENERGY
    !Config Def   = 273.15
    !Config Help  = The initial value if its value is not found
    !Config         in the restart file. This will only be used if the model is 
    !Config         started without a restart file.
    !Config Units = Kelvin [K]
    CALL setvar_p (T_snow_flk, val_exp,'T_SNOW_FLK', tp_00)


    ! Read or initialize T_ice_flk
    CALL restget_p (rest_id, 'T_ice_flk', nbp_glo, nlake, 1, kjit, .TRUE., T_ice_flk, "gather", nbp_glo, index_g)
    !Config Key   = T_ICE_FLK
    !Config Desc  = Initiale value of T_ice_flk used if not found in restart file
    !Config If    = OK_LAKE_ENERGY
    !Config Def   = 273.15
    !Config Help  = 
    !Config Units = Kelvin [K]
    CALL setvar_p(T_ice_flk, val_exp,'T_ICE_FLK', tp_00)


    ! Read or initialize T_wML_flk
    CALL restget_p (rest_id, 'T_wML_flk', nbp_glo, nlake, 1, kjit, .TRUE., T_wML_flk, "gather", nbp_glo, index_g)
    !Config Key   = T_WML_FLK
    !Config Desc  = Initiale value of T_wML_flk used if not found in restart file
    !Config If    = OK_LAKE_ENERGY
    !Config Def   = T_CLIM_LAKE
    !Config Help  = The initial value of mixed layer lake temperature if its value is not found
    !Config         in the restart file. This will only be used if the model is 
    !Config         started without a restart file.
    !Config Units = Kelvin [K]
    CALL setvar_p (T_wML_flk, val_exp,'T_WML_FLK', T_clim_lake)


    ! Read or initialize T_snow_flk
    CALL restget_p (rest_id, 'T_mnw_flk', nbp_glo, nlake, 1, kjit, .TRUE., T_mnw_flk, "gather", nbp_glo, index_g)
    !Config Key   = T_MNW_FLK
    !Config Desc  = Initiale value of T_mnw_flk used if not found in restart file
    !Config If    = OK_LAKE_ENERGY
    !Config Def   = T_CLIM_LAKE
    !Config Help  = The initial value of mean lake temperature if its value is not found
    !Config         in the restart file. This will only be used if the model is 
    !Config         started without a restart file.
    !Config Units = Kelvin [K]
    CALL setvar_p (T_mnw_flk, val_exp,'T_MNW_FLK', T_clim_lake)


    ! Read or initialize T_bot_flk
    CALL restget_p (rest_id, 'T_bot_flk', nbp_glo, nlake, 1, kjit, .TRUE.,T_bot_flk , "gather", nbp_glo, index_g)
    !Config Key   = T_BOT_FLK
    !Config Desc  = Initiale value of T_bot_flk used if not found in restart file
    !Config If    = OK_LAKE_ENERGY
    !Config Def   = T_CLIM_LAKE
    !Config Help  = The initial value if its value is not found
    !Config         in the restart file. This will only be used if the model is 
    !Config         started without a restart file.
    !Config Units = Kelvin [K]
    CALL setvar_p (T_bot_flk, val_exp,'T_BOT_FLK', T_clim_lake)


    ! Read or initialize T_B1_flk
    CALL restget_p (rest_id, 'T_B1_flk', nbp_glo, nlake, 1, kjit, .TRUE., T_B1_flk, "gather", nbp_glo, index_g)
    !Config Key   = T_B1_FLK
    !Config Desc  = Initiale value of T_B1_flk used if not found in restart file
    !Config If    = OK_LAKE_ENERGY
    !Config Def   = T_CLIM_LAKE
    !Config Help  = The initial value if its value is not found
    !Config         in the restart file. This will only be used if the model is 
    !Config         started without a restart file.
    !Config Units = Kelvin [K]
    CALL setvar_p (T_B1_flk, val_exp,'T_B1_FLK', T_clim_lake)


    ! Read or initialize C_T_flk
    CALL restget_p (rest_id, 'C_T_flk ', nbp_glo, nlake, 1, kjit, .TRUE., C_T_flk , "gather", nbp_glo, index_g)
    !Config Key   = C_T_FLK
    !Config Desc  = Initiale value of C_T_flk used if not found in restart file
    !Config If    = OK_LAKE_ENERGY
    !Config Def   = 0.6
    !Config Help  = The initial value if its value is not found
    !Config         in the restart file. This will only be used if the model is 
    !Config         started without a restart file.
    !Config Units = [-]
    CALL setvar_p (C_T_flk , val_exp,'C_T_FLK', 0.6_r_std)


    ! Read or initialize h_ML_flk
    CALL restget_p (rest_id, 'h_ML_flk', nbp_glo, nlake, 1, kjit, .TRUE., h_ML_flk, "gather", nbp_glo, index_g)
    !Config Key   = h_ML_FLK
    !Config Desc  = Initiale value of used if not found in restart file
    !Config If    = OK_LAKE_ENERGY
    !Config Def   = 0.1
    !Config Help  = The initial value if its value is not found
    !Config         in the restart file. This will only be used if the model is 
    !Config         started without a restart file.
    !Config Units = Meter [m]
    CALL setvar_p (h_ML_flk, val_exp,'h_ML_FLK', 0.1_r_std)


    ! Read or initialize h_ice_flk
    CALL restget_p (rest_id, 'h_ice_flk', nbp_glo, nlake, 1, kjit, .TRUE., h_ice_flk, "gather", nbp_glo, index_g)
    !Config Key   = H_ICE_FLK
    !Config Desc  = Initiale value of h_ice_flk used if not found in restart file
    !Config If    = OK_LAKE_ENERGY
    !Config Def   = 0.0
    !Config Help  = The initial value if its value is not found
    !Config         in the restart file. This will only be used if the model is 
    !Config         started without a restart file.
    !Config Units = Meter [m] 
    CALL setvar_p (h_ice_flk, val_exp,'H_ICE_FLK', zero)


    ! Read or initialize h_snow_flk
    CALL restget_p (rest_id, 'h_snow_flk', nbp_glo, nlake, 1, kjit, .TRUE., h_snow_flk, "gather", nbp_glo, index_g)
    !Config Key   = H_SNOW_FLK
    !Config Desc  = Initiale value of h_snow_flk used if not found in restart file
    !Config If    = OK_LAKE_ENERGY
    !Config Def   = 0.0
    !Config Help  = The initial value if its value is not found
    !Config         in the restart file. This will only be used if the model is 
    !Config         started without a restart file.
    !Config Units = Meter [m]
    CALL setvar_p (h_snow_flk, val_exp,'H_SNOW_FLK', zero)


    ! Read or initialize H_B1_flk
    CALL restget_p (rest_id, 'H_B1_flk', nbp_glo, nlake, 1, kjit, .TRUE., H_B1_flk, "gather", nbp_glo, index_g)
    !Config Key   = H_B1_FLK
    !Config Desc  = Initiale value of H_B1_flk used if not found in restart file
    !Config If    = OK_LAKE_ENERGY
    !Config Def   = 0.1
    !Config Help  = The initial value of if its value is not found
    !Config         in the restart file. This will only be used if the model is 
    !Config         started without a restart file.
    !Config Units = Meter [m]
    CALL setvar_p (H_B1_flk, val_exp,'H_B1_FLK', 0.1_r_std)   


  END SUBROUTINE lake_initialize


  !! ================================================================================================================================
  !! SUBROUTINE 	: lake_main
  !!
  !>\BRIEF         
  !!
  !! DESCRIPTION :
  !! - called every time step
  !! - initialization and finalization part are not done in here
  !!
  !! - 1 computes snow  ==> explicitsnow
  !! - 2 computes vegetations reservoirs  ==> lake_vegupd
  !! - 3 computes canopy  ==> lake_canop
  !! - 4 computes surface reservoir  ==> lake_flood
  !! - 5 computes soil lakeogy ==> lake_soil
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

  SUBROUTINE lake_main (kjit, lalo,         kjpindex,      index,         &
       u  ,v  , qair   , precip_rain , &
       precip_snow ,lwdown, swnet, swdown,pb ,temp_air,&
       vevapp_lake,   temp_surf_lake, fluxsens_lake, fluxlat_lake)


    !! 0. Variable and parameter declaration
    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                           :: kjit           !! Time step number 
    INTEGER(i_std), INTENT(in)                           :: kjpindex       !! Domain size
    REAL(r_std),DIMENSION (kjpindex,2), INTENT(in)       :: lalo           !! lalo 
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)     :: index          !! Indices of the points on the map
    REAL(r_std), DIMENSION (kjpindex),INTENT (in)        :: u              !! Lowest level wind speed in direction u (m/s)
    REAL(r_std), DIMENSION (kjpindex),INTENT (in)        :: v              !! Lowest level wind speed in direction v (m/s) 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: qair           !! Lowest level air specific humidity (kg/kg)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: precip_rain    !! Rain precipitation (kg/m2)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: precip_snow    !! Snow precipitation (kg/m2)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: lwdown         !! Downwelling long-wave flux (W/m2)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: swnet          !! Net surface short-wave flux (W/m2)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: swdown         !! Downwelling surface short-wave flux (W/m2) 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: temp_air       !! Air temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: pb             !! Surface pressure [hPa]

    !! 0.2 Modified variables
    REAL(r_std),DIMENSION (kjpindex, nlake), INTENT (inout)  :: temp_surf_lake     !! Lake surface temperature (K)

    !! 0.3 Output variables
    REAL(r_std),DIMENSION (kjpindex, nlake), INTENT (out)    :: fluxsens_lake     !! Sensible heat flux (W/m2)
    REAL(r_std),DIMENSION (kjpindex, nlake), INTENT (out)    :: fluxlat_lake      !! Latent heat flux (W/m2)
    REAL(r_std),DIMENSION (kjpindex, nlake), INTENT (out)    :: vevapp_lake       !! Total of evaporation (mm/day)


    !!  0.4 Local variables
    REAL(r_std) :: U_a_in
    REAL(r_std) :: P_a_in
    REAL(r_std) :: dM_snow 
    INTEGER(i_std) :: ij                !! Indice for the spatial loop
    INTEGER(i_std) :: ilake             !! Indice for the lake loop
    REAL (r_std) :: par_Coriolis

    !  Temperatures at the next time step  
    REAL (r_std) :: T_snow_flk_new
    REAL (r_std) :: T_ice_flk_new 
    REAL (r_std) :: T_mnw_flk_new 
    REAL (r_std) :: T_wML_flk_new 
    REAL (r_std) :: T_bot_flk_new 
    REAL (r_std) :: T_B1_flk_new  

    !  Thickness of various layers at the next time step 
    REAL (r_std) :: h_snow_flk_new
    REAL (r_std) :: h_ice_flk_new 
    REAL (r_std) :: h_ML_flk_new  
    REAL (r_std) :: H_B1_flk_new  
    REAL (r_std) :: temp_surf_lake_new
    REAL (r_std) :: frac_ice          !! Lake ice fraction
    REAL (r_std) :: h_crit
    REAL (r_std) :: albedo_lake

    !  The shape factor(s) at the next time step 

    REAL (r_std) :: C_T_flk_new       !! Shape factor (thermocline)
    REAL (r_std) :: C_TT_flk          !! Dimensionless parameter (thermocline)
    REAL (r_std) :: C_Q_flk           !! Shape factor with respect to the heat flux (thermocline)
    REAL (r_std) :: C_I_flk           !! Shape factor (ice)
    REAL (r_std) :: C_S_flk           !! Shape factor (snow)

    !  Heat and radiation fluxes

    REAL (r_std) :: Q_snow_flk        !! Heat flux through the air-snow interface [W/m2]
    REAL (r_std) :: Q_ice_flk         !! Heat flux through the snow-ice or air-ice interface [W/m2]
    REAL (r_std) :: Q_w_flk           !! Heat flux through the ice-water or air-water interface [W/m2]
    REAL (r_std) :: Q_bot_flk         !! Heat flux through the water-bottom sediment interface [W/m2]
    REAL (r_std) :: I_snow_flk        !! Radiation flux through the air-snow interface [W/m2]
    REAL (r_std) :: I_ice_flk         !! Radiation flux through the snow-ice or air-ice interface [W/m2]
    REAL (r_std) :: I_w_flk           !! Radiation flux through the ice-water or air-water interface [W/m2]
    REAL (r_std) :: I_h_flk           !! Radiation flux through the mixed-layer-thermocline interface [W/m2]
    REAL (r_std) :: I_bot_flk         !! Radiation flux through the water-bottom sediment interface [W/m2]
    REAL (r_std) :: I_intm_0_h_flk    !! Mean radiation flux over the mixed layer [W/m]
    REAL (r_std) :: I_intm_h_D_flk    !! Mean radiation flux over the thermocline [W/m]  


    !  Fluxes at the lake surface

    REAL (r_std) :: Q_momentum        !! Momentum flux [N/m2]
    REAL (r_std) :: Q_sensible        !! Sensible heat flux [W/m2]
    REAL (r_std) :: Q_latent          !! Latent heat flux [W/m2]
    REAL (r_std) :: Q_watvap          !! Flux of water vapour [kg/m2/s]

    REAL(r_std)   ::  Q_type

    REAL(r_std)   ::  qs



    !  Velocity scales

    REAL (r_std) :: u_star_w_flk      !! Friction velocity in the surface layer of lake water [m/s]


    ! Variables test


    REAL(r_std),DIMENSION (kjpindex, nlake) ::   albedo_water_lake_out
    REAL(r_std),DIMENSION (kjpindex, nlake) ::   frac_ice_out
    REAL(r_std),DIMENSION (kjpindex, nlake) ::   h_crit_out
    REAL(r_std),DIMENSION (kjpindex, nlake) ::   albedo_lake_out
    REAL(r_std),DIMENSION (kjpindex, nlake) ::   light_ext_water_lake_out 
    REAL(r_std),DIMENSION (kjpindex, nlake) ::   light_ext_ice_lake_out 
    REAL(r_std),DIMENSION (kjpindex, nlake) ::   light_ext_snow_lake_out


    REAL(r_std),DIMENSION (kjpindex, nlake) ::   I_h_flk_out
    REAL(r_std),DIMENSION (kjpindex, nlake) ::   I_intm_0_h_flk_out
    REAL(r_std),DIMENSION (kjpindex, nlake) ::   I_intm_h_D_flk_out    
    REAL(r_std),DIMENSION (kjpindex, nlake) ::   I_bot_flk_out  
    REAL(r_std),DIMENSION (kjpindex, nlake) ::   I_snow_flk_out    
    REAL(r_std),DIMENSION (kjpindex, nlake) ::   I_ice_flk_out 
    REAL(r_std),DIMENSION (kjpindex, nlake) ::   dM_snow_accepted

    REAL(r_std),DIMENSION (kjpindex, nlake) ::  depth_sed_lake_out
    REAL(r_std),DIMENSION (kjpindex, nlake) ::  t_sed_lake_out
    REAL(r_std),DIMENSION (kjpindex, nlake) ::  par_Coriolis_out
    REAL(r_std),DIMENSION (kjpindex, nlake) ::  Q_w_flk_out
    REAL(r_std),DIMENSION (kjpindex, nlake) ::  I_w_flk_out
    REAL(r_std),DIMENSION (kjpindex, nlake) ::  Q_momentum_out
    REAL(r_std),DIMENSION (kjpindex, nlake) ::  Q_snow_flk_out
    REAL(r_std),DIMENSION (kjpindex, nlake) ::  Q_ice_flk_out
    REAL(r_std),DIMENSION (kjpindex, nlake) ::  Q_bot_flk_out
    REAL(r_std),DIMENSION (kjpindex, nlake) ::  u_star_w_flk_out


    ! Loop on each grid point 
    DO ij = 1, kjpindex

       ! Calcul of the Coriolis parameter
       par_Coriolis = 2.* omega_earth * sin(2*3.1415*lalo(ij,1)/360)

       !1.  Transform the forcing file for FLake

       !1.1 Speed norm calculation
       U_a_in = (u(ij)**2 +v(ij)**2)**(0.5)
       !1.2 convert Pressure hPa in Pa 
       P_a_in = pb(ij)*100

       ! Loop on each lake
       DO ilake = 1, nlake
          IF (fraclake(ij,ilake)>0) THEN    
             !!          IF  (precip_snow(ij)/1000<0.0001) THEN
             dM_snow =  precip_snow(ij)/1000

             !!          ELSEIF ( (h_snow_flk(ij, ilake)+precip_snow(ij)/1000)/h_ice_flk(ij, ilake)>2) THEN
             !!              dM_snow = 0 
             !!          ELSEIF (( h_snow_flk(ij, ilake)<0.001) .AND. ( h_ice_flk(ij, ilake)<1  ) ) THEN
             !!              dM_snow = 0 
             !!          ELSE
             !!              dM_snow = MIN(precip_snow(ij)/1000, 0.2_r_std) ! Correct the forcing 
             !!          ENDIF


             !2.  Set initial values
             T_snow_flk_new = T_snow_flk(ij, ilake)
             T_ice_flk_new = T_ice_flk(ij, ilake)
             T_mnw_flk_new = T_mnw_flk(ij, ilake)
             T_wML_flk_new = T_wML_flk(ij, ilake)
             T_bot_flk_new = T_bot_flk(ij, ilake)
             T_B1_flk_new = T_B1_flk(ij, ilake)
             C_T_flk_new = C_T_flk(ij, ilake)
             h_snow_flk_new = h_snow_flk(ij, ilake)
             h_ice_flk_new = h_ice_flk(ij, ilake)
             h_ML_flk_new = h_ML_flk(ij, ilake)
             H_B1_flk_new = H_B1_flk(ij, ilake)
             temp_surf_lake_new = temp_surf_lake(ij, ilake)


             !------------------------------------------------------------------------------
             !  Set albedos of the lake water, lake ice and snow
             !------------------------------------------------------------------------------

             IF (ok_snow_ice_lake_par) THEN
                ! Use empirical formulation proposed by Mironov and Ritter (2004) for GME 
                !_nu albedo_ice_lake   = albedo_whiteice_ref 
                albedo_ice_lake   = EXP(-c_albice_MR*(tpl_T_f-temp_surf_lake(ij, ilake))/tpl_T_f)
                albedo_ice_lake   = albedo_whiteice_ref*(1.-albedo_ice_lake) + albedo_blueice_ref*albedo_ice_lake
                ! Snow is not considered
                !albedo_snow_lake  = albedo_ice_lake  

                ! Change from Mironov. 2012 Parametrisation of sea and lake ice in NWP models of GWS 
                albedo_snow_lake  =  EXP(-c_albice_MR*(tpl_T_f-temp_surf_lake(ij, ilake))/tpl_T_f)
                albedo_snow_lake = albedo_snow_max*(1.-albedo_snow_lake) + albedo_snow_min*albedo_snow_lake
             END IF

             h_crit = 0.15*fetch(ij, ilake)/27500

             !IF (h_crit>1.) THEN
             !	h_crit=1.
             !ENDIF 

             !IF (h_crit<0.10) THEN
             !	h_crit=0.10
             !ENDIF 

             ! Garnaud et al. (2022) parametrization for ice cover, https://doi.org/10.1029/2021MS002861
             ! Ice cover depends of the fetch of the lake and the properties of ice (more details in the article). 
             frac_ice = h_ice_flk(ij, ilake)/h_crit  ! x100 to have a percentage

             IF (frac_ice>1.) THEN
                frac_ice=1.
             ENDIF


             IF (h_ice_flk(ij, ilake).GE.h_Ice_min_flk) THEN             ! Ice exists
                IF (h_snow_flk(ij, ilake).GE.h_Snow_min_flk) THEN        ! There is snow above the ice
                   albedo_lake = frac_ice*albedo_snow_lake + (1-frac_ice)*albedo_water_lake
                ELSE                                                    ! No snow above the ice
                   albedo_lake = frac_ice*albedo_ice_lake + (1-frac_ice)*albedo_water_lake
                END IF
             ELSE                                                       ! No ice-snow cover
                albedo_lake = albedo_water_lake
             END IF

             !------------------------------------------------------------------------------
             !  Compute solar radiation fluxes (positive downward)
             !------------------------------------------------------------------------------


             CALL lake_radflux ( depth_w(ij, ilake), albedo_water_lake, albedo_ice_lake, albedo_snow_lake, & 
                  light_ext_water_lake(ij, ilake), light_ext_ice_lake, light_ext_snow_lake, &
                  h_snow_flk(ij, ilake), h_ice_flk(ij, ilake), h_ML_flk(ij, ilake),&
                  swdown(ij), I_w_flk, I_h_flk, I_intm_0_h_flk,   &
                  I_intm_h_D_flk, I_bot_flk, I_snow_flk, I_ice_flk,	&
                  frac_ice, h_crit, albedo_lake) 


             !------------------------------------------------------------------------------
             !  Compute long-wave radiation fluxes (positive downward)
             !------------------------------------------------------------------------------
             Q_w_flk = lwdown(ij) - c_lwrad_emis*tpsf_C_StefBoltz*temp_surf_lake(ij, ilake)**4  ! Radiation of the surface (notice the sign)


             !------------------------------------------------------------------------------
             !  Compute the surface friction velocity and fluxes of sensible and latent heat 
             !------------------------------------------------------------------------------
             CALL lake_SfcFlx_momsenlat ( height_u_lake, height_tq_lake, fetch(ij, ilake),                      &
                  U_a_in, temp_air(ij), qair(ij), temp_surf_lake(ij, ilake), P_a_in, h_ice_flk(ij, ilake),  &
                  Q_momentum, Q_sensible, Q_latent, Q_watvap, Q_type, qs )

             u_star_w_flk = SQRT(-Q_momentum/tpl_rho_w_r)

             !------------------------------------------------------------------------------
             !  Compute heat fluxes Q_snow_flk, Q_ice_flk, Q_w_flk
             !------------------------------------------------------------------------------
             Q_w_flk = Q_w_flk - Q_sensible - Q_latent  ! Add sensible and latent heat fluxes (notice the signs)
             IF( h_ice_flk(ij, ilake).GE.h_Ice_min_flk) THEN            ! Ice exists
                IF( h_snow_flk(ij, ilake).GE.h_Snow_min_flk) THEN        ! There is snow above the ice
                   Q_snow_flk = Q_w_flk
                   Q_ice_flk  = 0.
                   Q_w_flk    = 0.
                ELSE                                           ! No snow above the ice
                   Q_snow_flk = 0.
                   Q_ice_flk  = Q_w_flk
                   Q_w_flk    = 0.
                END IF
             ELSE                                             ! No ice-snow cover
                Q_snow_flk = 0.
                Q_ice_flk  = 0.
             END IF

             !------------------------------------------------------------------------------
             !  Advance FLake variables
             !------------------------------------------------------------------------------
             IF (h_ice_flk(ij, ilake)>0) THEN
                IF ((h_snow_flk(ij, ilake)+dM_snow)/h_ice_flk(ij, ilake)>1) THEN
                   dM_snow = 0
                ENDIF
             ENDIF


             CALL lake_driver ( depth_w(ij, ilake), depth_sed_lake, t_sed_lake, par_Coriolis,              &
                  light_ext_water_lake(ij, ilake),                  &
                  dt_sechiba, &                 
                  dM_snow, Q_w_flk, I_snow_flk, I_ice_flk, I_w_flk, I_h_flk, I_intm_0_h_flk, I_intm_h_D_flk, I_bot_flk, &
                  Q_snow_flk,  Q_ice_flk, Q_bot_flk, u_star_w_flk, &
                  temp_surf_lake(ij, ilake), T_snow_flk(ij, ilake), T_ice_flk(ij, ilake), T_wML_flk(ij, ilake),    &
                  T_mnw_flk(ij, ilake), T_bot_flk(ij, ilake), T_B1_flk(ij, ilake), h_snow_flk(ij, ilake), &
                  h_ice_flk(ij, ilake), h_ML_flk(ij, ilake), H_B1_flk(ij, ilake),   C_T_flk(ij, ilake),            &
                  temp_surf_lake_new,  T_snow_flk_new, T_ice_flk_new,  T_wML_flk_new,  &
                  T_mnw_flk_new, T_bot_flk_new, T_B1_flk_new, h_snow_flk_new, &
                  h_ice_flk_new, h_ML_flk_new, H_B1_flk_new, C_T_flk_new)


             !------------------------------------------------------------------------------
             !  Set output values
             !------------------------------------------------------------------------------
             T_snow_flk(ij, ilake) = T_snow_flk_new
             T_ice_flk(ij, ilake) = T_ice_flk_new
             T_mnw_flk(ij, ilake) = T_mnw_flk_new
             T_wML_flk(ij, ilake) = T_wML_flk_new
             T_bot_flk(ij, ilake) = T_bot_flk_new
             T_B1_flk(ij, ilake) = T_B1_flk_new
             C_T_flk(ij, ilake) = C_T_flk_new
             h_snow_flk(ij, ilake) = h_snow_flk_new
             h_ice_flk(ij, ilake) = h_ice_flk_new
             h_ML_flk(ij, ilake) = h_ML_flk_new
             H_B1_flk(ij, ilake) = H_B1_flk_new
             temp_surf_lake(ij, ilake) = temp_surf_lake_new
             fluxlat_lake(ij, ilake) = Q_latent 
             fluxsens_lake(ij, ilake) = Q_sensible
             vevapp_lake(ij, ilake) = Q_watvap
             regime_flux_lake(ij, ilake) = Q_type
             albedo_water_lake_out(ij, ilake) = albedo_water_lake
             light_ext_water_lake_out(ij, ilake)  = light_ext_water_lake(ij, ilake)
             light_ext_ice_lake_out(ij, ilake)  = light_ext_ice_lake
             light_ext_snow_lake_out(ij, ilake) = light_ext_snow_lake
             frac_ice_out(ij, ilake) = frac_ice
             h_crit_out(ij, ilake) = h_crit
             albedo_lake_out(ij, ilake) = albedo_lake
             qs_out(ij, ilake) = qs
             dM_snow_accepted(ij, ilake) = dM_snow*1000
             depth_sed_lake_out(ij, ilake) = depth_sed_lake
             t_sed_lake_out(ij, ilake) = t_sed_lake
             par_Coriolis_out(ij, ilake) = par_Coriolis
             Q_w_flk_out(ij, ilake) = Q_w_flk
             I_snow_flk_out(ij, ilake) = I_snow_flk
             I_ice_flk_out(ij, ilake) = I_ice_flk
             I_w_flk_out(ij, ilake) = I_w_flk
             I_h_flk_out(ij, ilake) = I_h_flk
             I_intm_0_h_flk_out(ij, ilake) = I_intm_0_h_flk
             I_intm_h_D_flk_out(ij, ilake) = I_intm_h_D_flk
             I_bot_flk_out(ij, ilake) = I_bot_flk
             Q_momentum_out(ij, ilake) = Q_momentum
             Q_snow_flk_out(ij, ilake) = Q_snow_flk
             Q_ice_flk_out(ij, ilake) = Q_ice_flk
             Q_bot_flk_out(ij, ilake) = Q_bot_flk
             u_star_w_flk_out(ij, ilake) =  u_star_w_flk 

          ELSE   

             regime_flux_lake(ij, ilake)= 0 
             dM_snow_accepted(ij, ilake)= undef_sechiba 
             T_snow_flk(ij, ilake) = undef_sechiba 
             T_ice_flk(ij, ilake) = undef_sechiba
             T_mnw_flk(ij, ilake) = undef_sechiba
             T_wML_flk(ij, ilake) = undef_sechiba
             T_bot_flk(ij, ilake) = undef_sechiba
             T_B1_flk(ij, ilake) = undef_sechiba
             C_T_flk(ij, ilake) = undef_sechiba
             h_snow_flk(ij, ilake) = undef_sechiba
             h_ice_flk(ij, ilake) = undef_sechiba
             h_ML_flk(ij, ilake) = undef_sechiba
             H_B1_flk(ij, ilake) = undef_sechiba
             temp_surf_lake(ij, ilake) = undef_sechiba
             fluxlat_lake(ij, ilake) =  undef_sechiba
             fluxsens_lake(ij, ilake) = undef_sechiba
             vevapp_lake(ij, ilake) = undef_sechiba

             albedo_water_lake_out(ij, ilake) = undef_sechiba
             light_ext_water_lake_out(ij, ilake)  = undef_sechiba
             light_ext_ice_lake_out(ij, ilake)  = undef_sechiba
             light_ext_snow_lake_out(ij, ilake) = undef_sechiba
             frac_ice_out(ij, ilake) = undef_sechiba
             h_crit_out(ij, ilake) = undef_sechiba
             albedo_lake_out(ij, ilake) = undef_sechiba
             qs_out(ij, ilake) = 0 

             depth_sed_lake_out(ij, ilake) = undef_sechiba
             t_sed_lake_out(ij, ilake) = undef_sechiba
             par_Coriolis_out(ij, ilake) = undef_sechiba
             Q_w_flk_out(ij, ilake) = undef_sechiba
             I_snow_flk_out(ij, ilake) = undef_sechiba
             I_ice_flk_out(ij, ilake) = undef_sechiba
             I_w_flk_out(ij, ilake) = undef_sechiba
             I_h_flk_out(ij, ilake) = undef_sechiba
             I_intm_0_h_flk_out(ij, ilake) = undef_sechiba
             I_intm_h_D_flk_out(ij, ilake) = undef_sechiba
             I_bot_flk_out(ij, ilake) = undef_sechiba
             Q_momentum_out(ij, ilake) = undef_sechiba
             Q_snow_flk_out(ij, ilake) = undef_sechiba
             Q_ice_flk_out(ij, ilake) = undef_sechiba
             Q_bot_flk_out(ij, ilake) = undef_sechiba
             u_star_w_flk_out(ij, ilake) =   undef_sechiba

          ENDIF

       ENDDO !ilake = 1, nlake
    ENDDO !ij = 1, kjpindex


    ! Write lake variables
    CALL xios_orchidee_send_field( "Tsfc_lake", temp_surf_lake)
    CALL xios_orchidee_send_field( "Qlatent_lake",   fluxlat_lake)
    CALL xios_orchidee_send_field( "Qsensible_lake", fluxsens_lake)  
    CALL xios_orchidee_send_field( "Qmomentum_lake", vevapp_lake)
    CALL xios_orchidee_send_field( "H_ML", h_ML_flk )
    CALL xios_orchidee_send_field( "H_ice", h_ice_flk )
    CALL xios_orchidee_send_field( "H_snow", h_snow_flk )
    CALL xios_orchidee_send_field( "C_T_thermocline", C_T_flk )
    CALL xios_orchidee_send_field( "T_snow_lake", T_snow_flk )
    CALL xios_orchidee_send_field( "T_ice_lake", T_ice_flk )
    CALL xios_orchidee_send_field( "T_ML_lake", T_wML_flk )
    CALL xios_orchidee_send_field( "T_bottom_lake", T_bot_flk )
    CALL xios_orchidee_send_field( "T_avg_lake", T_mnw_flk )
    CALL xios_orchidee_send_field( "depth_out", depth_w)
    CALL xios_orchidee_send_field( "frac_out", fraclake)
    CALL xios_orchidee_send_field( "fetch_out", fetch)
    CALL xios_orchidee_send_field( "Snowf_lake_accepted", dM_snow_accepted)
    CALL xios_orchidee_send_field( "albedo_water_lake_out", albedo_water_lake_out)
    CALL xios_orchidee_send_field( "frac_ice", frac_ice_out)
    CALL xios_orchidee_send_field( "h_crit", h_crit_out)
    CALL xios_orchidee_send_field( "albedo_lake", albedo_lake_out)
    CALL xios_orchidee_send_field( "light_ext_water_lake_out", light_ext_water_lake_out) 
    CALL xios_orchidee_send_field( "light_ext_ice_lake_out", light_ext_ice_lake_out)
    CALL xios_orchidee_send_field( "light_ext_snow_lake_out", light_ext_snow_lake_out)
    CALL xios_orchidee_send_field( "depth_sed_lake_out", depth_sed_lake_out)
    CALL xios_orchidee_send_field( "t_sed_lake_out", t_sed_lake_out)
    CALL xios_orchidee_send_field( "par_Coriolis_out", par_Coriolis_out)
    CALL xios_orchidee_send_field( "Q_w_flk_out", Q_w_flk_out)
    CALL xios_orchidee_send_field( "I_snow_flk_out", I_snow_flk_out)
    CALL xios_orchidee_send_field( "I_ice_flk_out", I_ice_flk_out)
    CALL xios_orchidee_send_field( "I_w_flk_out", I_w_flk_out)
    CALL xios_orchidee_send_field( "I_h_flk_out", I_h_flk_out)
    CALL xios_orchidee_send_field( "I_intm_0_h_flk_out", I_intm_0_h_flk_out)
    CALL xios_orchidee_send_field( "I_intm_h_D_flk_out", I_intm_h_D_flk_out)
    CALL xios_orchidee_send_field( "I_bot_flk_out", I_bot_flk_out)
    CALL xios_orchidee_send_field( "Q_momentum_out", Q_momentum_out)
    CALL xios_orchidee_send_field( "Q_snow_flk_out", Q_snow_flk_out)
    CALL xios_orchidee_send_field( "Q_ice_flk_out", Q_ice_flk_out)
    CALL xios_orchidee_send_field( "Q_bot_flk_out", Q_bot_flk_out)
    CALL xios_orchidee_send_field( "u_star_w_flk_out", u_star_w_flk_out)

  END SUBROUTINE lake_main


  !! ================================================================================================================================
  !! SUBROUTINE 	: lake_finalize
  !!
  !!
  !! DESCRIPTION : This subroutine writes the module variables and variables calculated in lake to restart file
  !!
  !! MAIN OUTPUT VARIABLE(S) : 
  !!
  !! REFERENCE(S) : 
  !!
  !! FLOWCHART    : None
  !! \n
  !_ ================================================================================================================================

  SUBROUTINE lake_finalize(kjit, kjpindex, rest_id, temp_surf_lake)

    !! 0. Variable and parameter declaration
    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                           :: kjit           !! Time step number 
    INTEGER(i_std), INTENT(in)                           :: kjpindex       !! Domain size
    INTEGER(i_std), INTENT (in)                          :: rest_id        !! Restart file identifier
    REAL(r_std),DIMENSION (kjpindex, nlake), INTENT (in) :: temp_surf_lake  !! Lake surface temperature (K)   


    !! 1. Write variables to restart file to be used for the next simulation
    IF (printlev_loc>=2) WRITE (numout,*) 'Write LAKE variables to the restart file'

    CALL restput_p(rest_id, 'fraclake', nbp_glo, nlake, 1, kjit,  fraclake, 'scatter',  nbp_glo, index_g)
    CALL restput_p(rest_id, 'depth_w', nbp_glo, nlake, 1, kjit,  depth_w, 'scatter',  nbp_glo, index_g)
    CALL restput_p(rest_id, 'light_ext_water_lake', nbp_glo, nlake, 1, kjit, light_ext_water_lake, 'scatter',  nbp_glo, index_g)
    CALL restput_p(rest_id, 'fetch', nbp_glo, nlake, 1, kjit, fetch, 'scatter',  nbp_glo, index_g)
    CALL restput_p(rest_id, 'temp_surf_lake', nbp_glo, nlake, 1, kjit,  temp_surf_lake, 'scatter',  nbp_glo, index_g)
    CALL restput_p(rest_id, 'T_snow_flk', nbp_glo, nlake, 1, kjit,  T_snow_flk, 'scatter',  nbp_glo, index_g)
    CALL restput_p(rest_id, 'T_ice_flk', nbp_glo, nlake, 1, kjit,  T_ice_flk, 'scatter',  nbp_glo, index_g)
    CALL restput_p(rest_id, 'T_wML_flk', nbp_glo, nlake, 1, kjit,  T_wML_flk, 'scatter',  nbp_glo, index_g)
    CALL restput_p(rest_id, 'T_mnw_flk', nbp_glo, nlake, 1, kjit,  T_mnw_flk, 'scatter',  nbp_glo, index_g)
    CALL restput_p(rest_id, 'T_bot_flk', nbp_glo, nlake, 1, kjit,  T_bot_flk, 'scatter',  nbp_glo, index_g)
    CALL restput_p(rest_id, 'T_B1_flk', nbp_glo, nlake, 1, kjit,  T_B1_flk, 'scatter',  nbp_glo, index_g)
    CALL restput_p(rest_id, 'C_T_flk', nbp_glo, nlake, 1, kjit,  C_T_flk, 'scatter',  nbp_glo, index_g)
    CALL restput_p(rest_id, 'h_ML_flk', nbp_glo, nlake, 1, kjit, h_ML_flk, 'scatter',  nbp_glo, index_g)
    CALL restput_p(rest_id, 'h_ice_flk', nbp_glo, nlake, 1, kjit, h_ice_flk, 'scatter',  nbp_glo, index_g)
    CALL restput_p(rest_id, 'h_snow_flk', nbp_glo,nlake, 1, kjit, h_snow_flk, 'scatter',  nbp_glo, index_g)
    CALL restput_p(rest_id, 'H_B1_flk', nbp_glo, nlake, 1, kjit, H_B1_flk, 'scatter',  nbp_glo, index_g)

    IF (printlev_loc>=3) WRITE (numout,*) 'Finish to write restart file with LAKE variables'
    ! ================================================================================================================================


  END SUBROUTINE lake_finalize



  !! ================================================================================================================================
  !! SUBROUTINE : lake_clear
  !!
  !! DESCRIPTION       Deallocate arrays 
  !!
  !! ================================================================================================================================
  SUBROUTINE lake_clear()

    ! Deallocation for lake module variables
    IF (ALLOCATED(T_mnw_flk)) DEALLOCATE (T_mnw_flk)
    IF (ALLOCATED(T_snow_flk)) DEALLOCATE (T_snow_flk)
    IF (ALLOCATED(T_ice_flk)) DEALLOCATE (T_ice_flk)
    IF (ALLOCATED(T_wML_flk)) DEALLOCATE (T_wML_flk)
    IF (ALLOCATED(T_bot_flk)) DEALLOCATE (T_bot_flk)
    IF (ALLOCATED(T_B1_flk)) DEALLOCATE (T_B1_flk)
    IF (ALLOCATED(C_T_flk)) DEALLOCATE (C_T_flk)
    IF (ALLOCATED(h_snow_flk)) DEALLOCATE (h_snow_flk)
    IF (ALLOCATED(h_ice_flk)) DEALLOCATE (h_ice_flk)
    IF (ALLOCATED(h_ML_flk)) DEALLOCATE (h_ML_flk)
    IF (ALLOCATED(H_B1_flk)) DEALLOCATE (H_B1_flk)
    IF (ALLOCATED(depth_w)) DEALLOCATE (depth_w)
    IF (ALLOCATED(fraclake)) DEALLOCATE (fraclake)
    IF (ALLOCATED(light_ext_water_lake)) DEALLOCATE (light_ext_water_lake)
    IF (ALLOCATED(fetch)) DEALLOCATE (fetch)
    IF (ALLOCATED(regime_flux_lake)) DEALLOCATE (regime_flux_lake)
    IF (ALLOCATED(qs_out)) DEALLOCATE (qs_out)


  END SUBROUTINE lake_clear



  !! ================================================================================================================================
  !! SUBROUTINE 	: lake_read_param
  !!
  !>\BRIEF         Read and process the nc file where is written the lake parameter
  !!
  !! DESCRIPTION : Read and process the nc file where is written the lake parameter 
  !!
  !! MAIN OUTPUT VARIABLE(S) : 
  !!
  !! REFERENCE(S) : 
  !!
  !! FLOWCHART    : None
  !! \n
  !! ================================================================================================================================

  SUBROUTINE lake_read_param (nbpt,  lalo, neighbours,  resolution, contfrac, lake_year, depth, light_ext_water_lake, fetch)

    !! 0. Variable and parameter declaration
    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                             :: nbpt            !! Number of points for which the data needs 
                                                                              !! to be interpolated
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in)             :: lalo            !! Vector of latitude and longitudes (beware of the order !)
    INTEGER(i_std), DIMENSION(nbpt,NbNeighb), INTENT(in)   :: neighbours      !! Vector of neighbours for each grid point
                                                                              !! (1=North and then clockwise)
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in)             :: resolution      !! The size in km of each grid-box in X and Y
    REAL(r_std), DIMENSION(nbpt), INTENT(in)               :: contfrac        !! Fraction of continent in the grid
    INTEGER(i_std), INTENT(in)                             :: lake_year       !! first year for landuse (0 == NO TIME AXIS)

    !! 0.2 Output variables
    REAL(r_std), DIMENSION(nbpt,nlake), INTENT(out)        :: depth 
    REAL(r_std), DIMENSION(nbpt,nlake), INTENT(out)        :: light_ext_water_lake
    REAL(r_std), DIMENSION(nbpt,nlake), INTENT(out)        :: fetch 

    !! 0.4 Local variables
    CHARACTER(LEN=80)                                      :: filename, name        !! Name of the netcdf file (unitless)
    INTEGER(i_std)                                         :: ier,ib,il             !! Check errors in memory allocation (unitless)
    REAL(r_std), DIMENSION(nbpt)                           :: alake                 !! Availability of the soilcol interpolation
    REAL(r_std), DIMENSION(nlake)                          :: vmin, vmax            !! min/max values to use for the renormalization
    CHARACTER(LEN=80)                                      :: variablename_depth    !! Variable to interpolate
    CHARACTER(LEN=80)                                      :: variablename_opticpar !! Variable to interpolate
    CHARACTER(LEN=80)                                      :: variablename_fetch    !! Variable to interpolate
    CHARACTER(LEN=80)                                      :: lonname, latname      !! lon, lat names in input file
    REAL(r_std), DIMENSION(nlake)                          :: variabletypevals      !! Values for all the types of the variable
                                                                                    !! (variabletypevals(1) = -un, not used)
    CHARACTER(LEN=50)                                      :: fractype              !! method of calculation of fraction
                                                                                    !!   'XYKindTime': Input values are kinds 
                                                                                    !!   of something with a temporal 
                                                                                    !!   evolution on the dx*dy matrix'
    LOGICAL                                                :: nonegative            !! whether negative values should be removed
    CHARACTER(LEN=50)                                      :: maskingtype           !! Type of masking
                                                                                    !!   'nomask': no-mask is applied
                                                                                    !!   'mbelow': take values below maskvals(1)
                                                                                    !!   'mabove': take values above maskvals(1)
                                                                                    !!   'msumrange': take values within 2 ranges;
                                                                                    !!        maskvals(2) <= SUM(vals(k)) <= maskvals(1)
                                                                                    !!        maskvals(1) < SUM(vals(k)) <= maskvals(3)
                                                                                    !!        (normalized by maskvals(3))
                                                                                    !!   'var': mask values are taken from a 
                                                                                    !!        variable inside the file (>0)
    REAL(r_std), DIMENSION(3)                              :: maskvals              !! values to use to mask (according to `maskingtype') 
    CHARACTER(LEN=250)                                     :: namemaskvar           !! name of the variable to use to mask 
    LOGICAL                                                :: new_status
    LOGICAL                                                :: old_status
    !_ ================================================================================================================================


    !! 1. Read the netcdf file to get the lake fraction

    !! 1.1 Get the name file and the dimension information
    !Config Key   = LAKE_PARAM_FILE
    !Config Desc  = Name of file which contains the spatial lake parameters
    !Config If    = OK_LAKE_ENERGY
    !Config Def   = lake_param.nc
    !Config Help  = 
    !Config Units = [FILE]
    filename = 'lake_param.nc'
    CALL getin_p('LAKE_PARAM_FILE',filename)
    name = filename(1:LEN_TRIM(FILENAME)-3)

    IF (xios_interpolation_lakeparam) THEN
       CALL xios_orchidee_recv_field('depth',depth)
       CALL xios_orchidee_recv_field('Cext',light_ext_water_lake)
       CALL xios_orchidee_recv_field('fetch', fetch)
    ELSE
       variablename_depth = 'depth'
       variablename_opticpar = 'Cext'
       variablename_fetch = 'fetch'

       IF (printlev_loc >= 1) WRITE(numout,*) "lake_read_param: Start interpolate " &
            // TRIM(filename) // " for variable " // TRIM(variablename_fetch)

       ! Assigning values to vmin, vmax
       vmin = 1
       vmax = nlake*1._r_std

       variabletypevals = -un

       !! Variables for interpweight
       ! Type of calculation of cell fractions
       fractype = 'default'
       ! Name of the longitude and latitude in the input file
       lonname = 'lon'
       latname = 'lat'
       ! Should negative values be set to zero from input file?
       nonegative = .FALSE.
       ! Type of mask to apply to the input data (see header for more details)
       maskingtype = 'nomask'
       ! Values to use for the masking
       maskvals = (/ 50.+1.e-7, 0., 51. /)
       ! Name of the variable with the values for the mask in the input file (only if maskkingtype='var') (here not used)
       namemaskvar = ''

       CALL interpweight_3D_lake(nbpt, nlake, variabletypevals, lalo, resolution, neighbours,        &
            contfrac, filename, variablename_fetch, lonname, latname, vmin, vmax, nonegative, maskingtype,        &
            maskvals, namemaskvar, nlake, 0, lake_year, fractype,                                 &
            -1., -1., fetch, alake)

       IF (printlev_loc >= 1) WRITE(numout,*) "lake_read_param: Start interpolate " &
            // TRIM(filename) // " for variable " // TRIM(variablename_depth)

       ! Assigning values to vmin, vmax
       vmin = 1
       vmax = nlake*1._r_std

       variabletypevals = -un

       !! Variables for interpweight
       ! Type of calculation of cell fractions
       fractype = 'default'
       ! Name of the longitude and latitude in the input file
       lonname = 'lon'
       latname = 'lat'
       ! Should negative values be set to zero from input file?
       nonegative = .FALSE.
       ! Type of mask to apply to the input data (see header for more details)
       maskingtype = 'nomask'
       ! Values to use for the masking
       maskvals = (/ 50.+1.e-7, 0., 51. /)
       ! Name of the variable with the values for the mask in the input file (only if maskkingtype='var') (here not used)
       namemaskvar = ''

       CALL interpweight_3D_lake(nbpt, nlake, variabletypevals, lalo, resolution, neighbours,        &
            contfrac, filename, variablename_depth, lonname, latname, vmin, vmax, nonegative, maskingtype,        &
            maskvals, namemaskvar, nlake, 0, lake_year, fractype,                                 &
            -1., -1., depth, alake)


       IF (printlev_loc >= 1) WRITE(numout,*) "lake_read_param: Start interpolate " &
            // TRIM(filename) // " for variable " // TRIM(variablename_opticpar)

       !! Variables for interpweight
       ! Type of calculation of cell fractions
       fractype = 'default'
       ! Name of the longitude and latitude in the input file
       lonname = 'lon'
       latname = 'lat'
       ! Should negative values be set to zero from input file?
       nonegative = .FALSE.
       ! Type of mask to apply to the input data (see header for more details)
       maskingtype = 'nomask'
       ! Values to use for the masking
       maskvals = (/ 10., 0., 20. /)
       ! Name of the variable with the values for the mask in the input file (only if maskkingtype='var') (here not used)
       namemaskvar = ''
       CALL interpweight_3D_lake(nbpt, nlake, variabletypevals, lalo, resolution, neighbours,        &
            contfrac, filename, variablename_opticpar, lonname, latname, vmin, vmax, nonegative, maskingtype,        &
            maskvals, namemaskvar, nlake, 0, lake_year, fractype,                                 &
            -1., -1., light_ext_water_lake, alake)
    ENDIF

    CALL xios_orchidee_send_field('interp_diag_depth', depth)  
    CALL xios_orchidee_send_field('interp_diag_Cext', light_ext_water_lake)
    CALL xios_orchidee_send_field('interp_diag_fetch', fetch)


  END SUBROUTINE lake_read_param




  !! ================================================================================================================================
  !! SUBROUTINE 	: lake_read_frac
  !!
  !>\BRIEF         Read and interpolate the netcdf file where is written the lake fraction
  !!
  !! DESCRIPTION : Read and interpolate the netcdf file where is written the lake fraction. 
  !!
  !! MAIN OUTPUT VARIABLE(S) : 
  !!
  !! REFERENCE(S) : 
  !!
  !! FLOWCHART    : None
  !! \n
  !_ ================================================================================================================================

  SUBROUTINE lake_read_frac(nbpt,  lalo, neighbours,  resolution, contfrac, lake_year, fraclake)

    !! 0. Variable and parameter declaration
    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                             :: nbpt           !! Number of points for which the data needs 
                                                                             !! to be interpolated
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in)             :: lalo           !! Vector of latitude and longitudes (beware of the order !)
    INTEGER(i_std), DIMENSION(nbpt,NbNeighb), INTENT(in)   :: neighbours     !! Vector of neighbours for each grid point
                                                                             !! (1=North and then clockwise)
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in)             :: resolution     !! The size in km of each grid-box in X and Y
    REAL(r_std), DIMENSION(nbpt), INTENT(in)               :: contfrac       !! Fraction of continent in the grid

    INTEGER(i_std), INTENT(in)                             :: lake_year      !! first year for landuse (0 == NO TIME AXIS)

    !! 0.2 Output variables
    REAL(r_std), DIMENSION(nbpt,nlake), INTENT(out)        :: fraclake

    !! 0.4 Local variables
    CHARACTER(LEN=80)                                    :: filename, name   !! Name of the netcdf file (unitless)
    REAL(r_std), DIMENSION(nbpt)                         :: alake            !! Availability of the soilcol interpolation
    REAL(r_std), DIMENSION(nlake)                        :: vmin, vmax       !! min/max values to use for the renormalization
    CHARACTER(LEN=80)                                    :: variablename     !! Variable to interpolate
    CHARACTER(LEN=80)                                    :: lonname, latname !! lon, lat names in input file
    REAL(r_std), DIMENSION(nlake)                        :: variabletypevals !! Values for all the types of the variable
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
                                                                             !!        (normalized by maskvals(3))
                                                                             !!   'var': mask values are taken from a 
                                                                             !!     variable inside the file (>0)
    REAL(r_std), DIMENSION(3)                            :: maskvals         !! values to use to mask (according to `maskingtype') 
    CHARACTER(LEN=250)                                   :: namemaskvar      !! name of the variable to use to mask 
    LOGICAL                                              :: new_status
    LOGICAL                                              :: old_status
    ! ================================================================================================================================


    !! 1. Read the netcdf file to get the lake fraction

    !! 1.1 Get the name file and the dimension information

    !Config Key   = LAKE_FRAC_FILE
    !Config Desc  = Name of file which contains the lake fractions
    !Config If    = OK_LAKE_ENERGY
    !Config Def   = lake_frac.nc
    !Config Help  = 
    !Config Units = [FILE]
    filename = 'lake_frac.nc'
    CALL getin_p('LAKE_FRAC_FILE',filename)
    name = filename(1:LEN_TRIM(FILENAME)-3)

    IF (xios_interpolation_lakefrac) THEN
       !      CALL xios_orchidee_set_file_attr("lake_frac_file",name=name)
       CALL xios_orchidee_recv_field('fraclake',fraclake)
    ELSE
       variablename = 'fraclake'
       IF (printlev_loc >= 1) WRITE(numout,*) "lake_read_frac: Start interpolate " &
            // TRIM(filename) // " for variable " // TRIM(variablename)
       ! Assigning values to vmin, vmax
       vmin = 1
       vmax = nlake*1._r_std

       variabletypevals = -un

       !! Variables for interpweight
       ! Type of calculation of cell fractions
       fractype = 'default'
       ! Name of the longitude and latitude in the input file
       lonname = 'lon'
       latname = 'lat'
       ! Should negative values be set to zero from input file?
       nonegative = .FALSE.
       ! Type of mask to apply to the input data (see header for more details)
       maskingtype = 'nomask'
       ! Values to use for the masking
       maskvals = (/ 1.-1.e-7, 0., 2. /)
       ! Name of the variable with the values for the mask in the input file (only if maskkingtype='var') (here not used)
       namemaskvar = 'nomask'


       CALL interpweight_3D_lake(nbpt, nlake, variabletypevals, lalo, resolution, neighbours,        &
            contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype,        &
            maskvals, namemaskvar, nlake, 0, lake_year, fractype,                                 &
            -1., -1., fraclake, alake)

    ENDIF

    CALL xios_orchidee_send_field("interp_diag_fraclake",fraclake)

  END SUBROUTINE lake_read_frac


  !! ================================================================================================================================
  !! SUBROUTINE 	: lake_adjust_frac_and_depth
  !!
  !! DESCRIPTION : Remove too small fractions of lake surface.
  !!               Adjust lake depth depending on the options ok_depth_lake_max and lake_depth_cons>0.
  !!
  !! MAIN OUTPUT VARIABLE(S) : 
  !!
  !! REFERENCE(S) : 
  !!
  !! FLOWCHART    : None
  !! \n
  !_ ================================================================================================================================
  SUBROUTINE lake_adjust_frac_and_depth(kjpindex, depth_w, fraclake)

    !! 0. Variable and parameter declaration
    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                                 :: kjpindex       !! Domain size (-)

    !! 0.2 Modified variables
    REAL(r_std), DIMENSION(kjpindex,nlake), INTENT(inout)      :: fraclake       !! 
    REAL(r_std), DIMENSION(kjpindex,nlake), INTENT(inout)      :: depth_w        !! 

    !! 0.4 Local variables
    REAL(r_std)                                                :: depth_min = 0.01
    REAL(r_std)                                                :: frac_err = 0.05
    INTEGER(i_std)                                             :: ilake, ij         !!

    !_ ================================================================================================================================

    DO ij = 1, kjpindex 
       DO ilake = 1, nlake

          IF (fraclake(ij, ilake)>0 .AND. depth_w(ij, ilake)>depth_min ) THEN

             IF (fraclake(ij, ilake)>1) THEN
                fraclake(ij, ilake)=1
             END IF

             ! Adjust the lake depth depending on specific options set in run.def
             IF (ok_depth_lake_max) THEN
                ! Use a maximum allowed depth for the lake
                IF (depth_w(ij,ilake) > depth_lake_max) THEN
                   depth_w(ij,ilake) = depth_lake_max
                END IF
             ELSE IF(lake_depth_cons>0) THEN
                ! Overwrite lake depth with constant value red from run.def
                depth_w(ij, ilake)=lake_depth_cons
             END IF
             
          ELSE
             ! No fraclake was found or depth_w was too low
             ! Now set fraclake to zero
             fraclake(ij,ilake) = 0
          END IF
       END DO
    END DO

  END SUBROUTINE lake_adjust_frac_and_depth



  !! ================================================================================================================================
  !! SUBROUTINE 	: lake_radflux
  !!
  !>\BRIEF         Computes the radiation fluxes
  !!
  !! DESCRIPTION : Computes the radiation fluxes 
  !!               at the snow-ice, ice-water, air-water, 
  !!               mixed layer-thermocline and water column-bottom sediment interfaces,
  !!               the mean radiation flux over the mixed layer,
  !!               and the mean radiation flux over the thermocline.
  !!
  !! MAIN OUTPUT VARIABLE(S) :  I_w_flk,       I_h_flk,   I_intm_0_h_flk,   
  !!                            I_intm_h_D_flk,    I_bot_flk,   I_snow_flk ,     I_ice_flk
  !!
  !! REFERENCE(S) : 
  !!
  !! FLOWCHART    : None
  !! \n
  !_ ================================================================================================================================
  SUBROUTINE lake_radflux ( depth_w, albedo_water_lake, albedo_ice_lake, albedo_snow_lake, & 
       light_ext_water_lake, light_ext_ice_lake, light_ext_snow_lake , h_snow_p_flk, h_ice_p_flk, h_ML_p_flk, &
       I_atm_flk ,    I_w_flk,       I_h_flk,   I_intm_0_h_flk,   &
       I_intm_h_D_flk,    I_bot_flk,   I_snow_flk ,     I_ice_flk,  &
       frac_ice,	h_crit, 	albedo_lake)     


    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    !! 0.1.1 Atmospheric variables 
    REAL (r_std), INTENT(IN) :: I_atm_flk                              !! Downwelling surface short-wave flux 

    !! 0.1.2 FLake parameteris
    REAL (r_std), INTENT(IN) :: depth_w                                !! The lake depth [m]
    REAL (r_std), INTENT(IN) :: albedo_water_lake                      !! Albedo of the water surface 
    REAL (r_std), INTENT(IN) :: albedo_ice_lake                        !! Albedo of the ice surface
    REAL (r_std), INTENT(IN) :: albedo_snow_lake                       !! Albedo of the snow surface
    REAL (r_std), INTENT(IN) :: frac_ice                               !! Ice cover of the lake according to Garnaud's article
    REAL (r_std), INTENT(IN) :: h_crit                                 !! Critical ice height
    REAL (r_std), INTENT(IN) :: albedo_lake                            !! Lake albedo considering Garnaud's ice fraction
    REAL (r_std), INTENT(IN) :: light_ext_water_lake                   !! Optical characteristics of water
    REAL (r_std), INTENT(IN) :: light_ext_ice_lake                     !! Optical characteristics of ice
    REAL (r_std), INTENT(IN) :: light_ext_snow_lake                    !! Optical characteristics of snow 


    !! 0.1.3 FLake state variables
    REAL (r_std), INTENT(IN) :: h_ice_p_flk
    REAL (r_std), INTENT(IN) :: h_snow_p_flk
    REAL (r_std), INTENT(IN) :: h_ML_p_flk


    !! 0.2 Input variables
    REAL(r_std) , INTENT(OUT) :: I_snow_flk                          !! Radiation flux through the air-snow interface [W/m2]
    REAL(r_std) , INTENT(OUT) :: I_ice_flk                           !! Radiation flux through the snow-ice or air-ice interface [W/m2]
    REAL(r_std) , INTENT(OUT) :: I_w_flk                             !! Radiation flux through the ice-water or air-water interface [W/m2]
    REAL(r_std) , INTENT(OUT) :: I_h_flk                             !! Radiation flux through the mixed-layer-thermocline interface [W/m2]
    REAL(r_std) , INTENT(OUT) :: I_bot_flk                           !! Radiation flux through the water-bottom sediment interface [W/m2]
    REAL(r_std) , INTENT(OUT) :: I_intm_0_h_flk                      !! Mean radiation flux over the mixed layer [W/m]
    REAL(r_std) , INTENT(OUT) :: I_intm_h_D_flk                      !! Mean radiation flux over the thermocline [W/m]


    !_ ================================================================================================================================

    !1. Calculation of the radiation flux at each interface: snow, ice, water 
    IF (h_ice_p_flk.GE.h_Ice_min_flk) THEN            !A) Ice exists

       IF(h_snow_p_flk.GE.h_Snow_min_flk) THEN        !B) There is snow above the ice
          ! flux penetrating the lake is moduled by ice fraction and snow albedo
          I_snow_flk = I_atm_flk*((1.-albedo_snow_lake)*frac_ice+(1.-albedo_water_lake)*(1-frac_ice))

          ! Calculation for each wave length band contribution
          I_bot_flk = EXP(-light_ext_snow_lake*h_snow_p_flk)
          I_ice_flk  = I_snow_flk*I_bot_flk

       ELSE                                           !B) No snow above the ice 
          I_snow_flk = I_atm_flk  
          ! simply ice fraction and ice albedo if no snow-fall
          I_ice_flk  = I_atm_flk*((1.-albedo_ice_lake)*frac_ice+(1.-albedo_water_lake)*(1-frac_ice)) 
       END IF

       ! Calculation for each wave length band contribution
       I_bot_flk = EXP(-light_ext_ice_lake*h_ice_p_flk)
       I_w_flk      = I_ice_flk*I_bot_flk

    ELSE                                              !A) No ice-snow cover
       I_snow_flk   = I_atm_flk
       I_ice_flk    = I_atm_flk
       I_w_flk      = I_atm_flk*(1.-albedo_water_lake)
    END IF



    !2. Calculation of the Radiation flux at the bottom of the mixed layer
    IF (h_ML_p_flk.GE.h_ML_min_flk) THEN               !A) Security: Mixed layer large enough 

       ! Calculation for each wave length band contribution
       I_bot_flk = EXP(-light_ext_water_lake*h_ML_p_flk)
       I_h_flk = I_w_flk*I_bot_flk

    ELSE                                              !A) Mixed-layer depth is less then a minimum value
       I_h_flk = I_w_flk
    END IF !A) Security: Mixed layer large enough                                       


    !3. Calculation of Radiation flux at the lake bottom
    I_bot_flk =  EXP(-light_ext_water_lake*depth_w)                        
    I_bot_flk = I_w_flk*I_bot_flk


    !4. Calculation of Integral-mean radiation flux over the mixed layer
    IF(h_ML_p_flk.GE.h_ML_min_flk) THEN               !A) Security: Mixed layer large enough           
       I_intm_0_h_flk = 1./light_ext_water_lake*  &
            (1. - EXP(-light_ext_water_lake*h_ML_p_flk))
       I_intm_0_h_flk = I_w_flk*I_intm_0_h_flk/h_ML_p_flk
    ELSE                                              !A) Mixed-layer depth is less then a minimum value
       I_intm_0_h_flk = I_h_flk !also egal to I_w_flk
    END IF !A) Security: Mixed layer large enough           


    !5. Calculation of  Integral-mean radiation flux over the thermocline
    IF(h_ML_p_flk.LE.depth_w-h_ML_min_flk) THEN       !A) Security: thermocline bigger than h_ML_min  
       I_intm_h_D_flk = 1./light_ext_water_lake*  &
            ( EXP(-light_ext_water_lake*h_ML_p_flk)             &
            - EXP(-light_ext_water_lake*depth_w))
       I_intm_h_D_flk = I_w_flk*I_intm_h_D_flk/(depth_w-h_ML_p_flk)
    ELSE                                              !A) Thermocline too small
       I_intm_h_D_flk = I_h_flk
    END IF !A) Security: thermocline bigger than h_ML_min  

  END SUBROUTINE lake_radflux


  !! ================================================================================================================================
  !! SUBROUTINE 	: lake_SfcFlx_momsenlat
  !!
  !>\BRIEF          The SfcFlx routine 
  !!                where fluxes of momentum and of sensible and latent heat 
  !!                at the air-water or air-ice (air-snow) interface are computed. 
  !!
  !! DESCRIPTION :  The SfcFlx routine 
  !!                where fluxes of momentum and of sensible and latent heat 
  !!                at the air-water or air-ice (air-snow) interface are computed. 
  !!
  !! MAIN OUTPUT VARIABLE(S) : Q_momentum, Q_sensible, Q_latent, Q_watvap
  !!
  !! REFERENCE(S) : 
  !!
  !! FLOWCHART    : None
  !! \n
  !_ ================================================================================================================================



  SUBROUTINE lake_SfcFlx_momsenlat ( height_u, height_tq, fetch,                &
       U_a, T_a, q_a, T_s, P_a, h_ice,            &
       Q_momentum, Q_sensible, Q_latent, Q_watvap, regime_flux_lake_r, q_s ) 


    !! 0. Variable and parameter declaration
    !! 0.1 Input variables
    !! 0.1.1 Atmospheric variables 
    REAL (r_std), INTENT(IN) ::  U_a               !! Wind speed [m s^{-1}]
    REAL (r_std), INTENT(IN) ::  T_a               !! Air temperature [K]
    REAL (r_std), INTENT(IN) ::  q_a               !! Air specific humidity [-]
    REAL (r_std), INTENT(IN) ::  T_s               !! Surface temperature (water, ice or snow) [K]
    REAL (r_std), INTENT(IN) ::  P_a               !! Surface air pressure [N m^{-2} = kg m^{-1} s^{-2}]

    !! 0.1.2 FLake parameters
    REAL (r_std), INTENT(IN) ::  height_u          !! Height where wind is measured [m]
    REAL (r_std), INTENT(IN) ::  height_tq         !! Height where temperature and humidity are measured [m]
    REAL (r_std), INTENT(IN) ::  fetch             !! Typical wind fetch [m]

    !! 0.1.3 FLake state variables
    REAL (r_std), INTENT(IN) ::  h_ice             !! Ice thickness [m]

    !! 0.2 Output variables
    REAL (r_std), INTENT(OUT) ::  Q_momentum       !! Momentum flux [N m^{-2}]  
    REAL (r_std), INTENT(OUT) ::  Q_sensible       !! Sensible heat flux [W m^{-2}]  
    REAL (r_std), INTENT(OUT) ::  Q_latent         !! Laten heat flux [W m^{-2}]
    REAL (r_std), INTENT(OUT) ::  Q_watvap         !! Flux of water vapout [kg m^{-2} s^{-1}]
    REAL (r_std), INTENT(OUT) ::  regime_flux_lake_r !!

    !! 0.4 Local variables
    INTEGER (i_std), PARAMETER :: n_iter_max = 24 !! Maximum number of iterations 
    LOGICAL :: l_conv_visc                        !! Switch, TRUE = viscous free convection, the Nu=C Ra^(1/3) law is used
    LOGICAL :: l_conv_cbl                         !! Switch, TRUE = CBL scale convective structures define surface fluxes 

    INTEGER (i_std) ::   i                        !! Loop index
    INTEGER (i_std) ::   n_iter                   !! Number of iterations performed 

    REAL (r_std) :: rho_a                         !! Air density [kg m^{-3}]  
    REAL (r_std) :: wvpres_s                      !! Saturation water vapour pressure at T=T_s [N m^{-2}]
    REAL (r_std), INTENT(OUT)  :: q_s                           !! Saturation specific humidity at T=T_s [-]

    REAL (r_std) :: Q_mom_tur                     !! Turbulent momentum flux [N m^{-2}]
    REAL (r_std) :: Q_sen_tur                     !! Turbulent sensible heat flux [W m^{-2}]  
    REAL (r_std) :: Q_lat_tur                     !! Turbulent laten heat flux [W m^{-2}]
    REAL (r_std) :: Q_mom_mol                     !! Molecular momentum flux [N m^{-2}]
    REAL (r_std) :: Q_sen_mol                     !! Molecular sensible heat flux [W m^{-2}]  
    REAL (r_std) :: Q_lat_mol                     !! Molecular laten heat flux [W m^{-2}]
    REAL (r_std) :: Q_mom_con                     !! Momentum flux in free convection [N m^{-2}]
    REAL (r_std) :: Q_sen_con                     !! Sensible heat flux in free convection [W m^{-2}]  
    REAL (r_std) :: Q_lat_con                     !! Latent heat flux in free convection [W m^{-2}]

    REAL (r_std) :: par_conv_visc                 !! Viscous convection stability parameter
    REAL (r_std) :: par_conv_cbl                  !! CBL convection stability parameter
    REAL (r_std) :: c_z0u_fetch                   !! Fetch-dependent Charnock parameter
    REAL (r_std) :: U_a_thresh                    !! Threshld value of the wind speed [m s^{-1}] 
    REAL (r_std) :: u_star_thresh                 !! Threshld value of friction velocity [m s^{-1}]
    REAL (r_std) :: u_star_previter               !! Friction velocity from previous iteration [m s^{-1}]
    REAL (r_std) :: u_star_n                      !! Friction velocity at neutral stratification [m s^{-1}]
    REAL (r_std) :: u_star_st                     !! Friction velocity with due regard for stratification [m s^{-1}]
    REAL (r_std) :: ZoL                           !! The z/L ratio, z=height_u
    REAL (r_std) :: Ri                            !! Gradient Richardson number 
    REAL (r_std) :: Ri_cr                         !! Critical value of Ri 
    REAL (r_std) :: R_z                           !! Ratio of "height_tq" to "height_u"
    REAL (r_std) :: Fun                           !! A function of generic variable "x"
    REAL (r_std) :: Fun_prime                     !! Derivative of "Fun" with respect to "x"
    REAL (r_std) :: Delta                         !! Relative error 
    REAL (r_std) :: psi_u                         !! The MO stability function for wind profile
    REAL (r_std) :: psi_t                         !! The MO stability function for temperature profile
    REAL (r_std) :: psi_q                         !! The MO stability function for specific humidity profile

    REAL (r_std) :: z0u_sf                        !! Roughness length with respect to wind velocity [m]
    REAL (r_std) :: z0t_sf                        !! Roughness length with respect to potential temperature [m]
    REAL (r_std) :: z0q_sf                        !! Roughness length with respect to specific humidity [m]

    INTEGER ::  regime_flux_lake_temp !!
    !_ ================================================================================================================================

    !_dm All fluxes are positive when directed upwards.

    !1.  Compute saturation specific humidity and the air density at T=T_s

    ! Saturation water vapour pressure [N m^{-2} = kg m^{-1} s^{-2}]
    ! Saturation water vapour pressure at T=T_s
    IF(h_ice.LT.h_Ice_min_flk) THEN  ! Water surface
       wvpres_s = b1_vap*EXP(b2w_vap*(T_s-b3_vap)/(T_s-b4w_vap))
    ELSE                             ! Ice surface
       wvpres_s = b1_vap*EXP(b2i_vap*(T_s-b3_vap)/(T_s-b4i_vap))
    END IF

    q_s = tpsf_Rd_o_Rv*wvpres_s/(P_a-(1.-tpsf_Rd_o_Rv)*wvpres_s)     ! Saturation specific humidity at T=T_s

    rho_a = P_a/tpsf_R_dryair/T_s/(1.+(1./tpsf_Rd_o_Rv-1.)*q_s)    ! Air density at T_s and q_s (surface values)


    !2. Compute Fluxes of momentum and of sensible and latent heat
    !2.1  Compute molecular fluxes of momentum and of sensible and latent heat

    !_dm The fluxes are in kinematic units
    Q_mom_mol = -tpsf_nu_u_a*U_a/height_u 
    Q_sen_mol = -tpsf_kappa_t_a*(T_a-T_s)/height_tq    
    Q_lat_mol = -tpsf_kappa_q_a*(q_a-q_s)/height_tq  

    !2.2  Compute fluxes in free convection
    par_conv_visc = (T_s-T_a)/T_s*SQRT(tpsf_kappa_t_a) + (q_s-q_a)*tpsf_alpha_q*SQRT(tpsf_kappa_q_a)
    IF(par_conv_visc.GT.0.) THEN   ! Viscous convection takes place
       l_conv_visc = .TRUE.
       par_conv_visc = (par_conv_visc*tpl_grav/tpsf_nu_u_a)**num_1o3_sf
       Q_sen_con = c_free_conv*SQRT(tpsf_kappa_t_a)*par_conv_visc  
       Q_sen_con = Q_sen_con*(T_s-T_a)
       Q_lat_con = c_free_conv*SQRT(tpsf_kappa_q_a)*par_conv_visc
       Q_lat_con = Q_lat_con*(q_s-q_a)
    ELSE                                  ! No viscous convection, set fluxes to zero
       l_conv_visc = .FALSE.
       Q_sen_con = 0. 
       Q_lat_con = 0.
    END IF
    Q_mom_con = 0.                 ! Momentum flux in free (viscous or CBL-scale) convection is zero  

    !2.3  Compute turbulent fluxes
    R_z   = height_tq/height_u                        ! Ratio of "height_tq" to "height_u"
    Ri_cr = c_MO_t_stab/c_MO_u_stab**2*R_z  ! Critical Ri
    Ri    = tpl_grav*((T_a-T_s)/T_s+tpsf_alpha_q*(q_a-q_s))/MAX(U_a,u_wind_min_sf)**2
    Ri    = Ri*height_u/Pr_neutral                    ! Gradient Richardson number

    Turb_Fluxes: IF(U_a.LT.u_wind_min_sf.OR.Ri.GT.Ri_cr-c_small_sf) THEN  ! Low wind or Ri>Ri_cr 
       u_star_st = 0.                       ! Set turbulent fluxes to zero 
       Q_mom_tur = 0.                       
       Q_sen_tur = 0.   
       Q_lat_tur = 0.  

    ELSE Turb_Fluxes                            ! Compute turbulent fluxes using MO similarity

       ! Compute z/L, where z=height_u
       IF(Ri.GE.0.) THEN   ! Stable stratification
          ZoL = SQRT(1.-4.*(c_MO_u_stab-R_z*c_MO_t_stab)*Ri)
          ZoL = ZoL - 1. + 2.*c_MO_u_stab*Ri
          ZoL = ZoL/2./c_MO_u_stab/c_MO_u_stab/(Ri_cr-Ri)
       ELSE                       ! Convection
          n_iter = 0
          Delta = 1.                ! Set initial error to a large value (as compared to the accuracy)
          u_star_previter = Ri*MAX(1., SQRT(R_z*c_MO_t_conv/c_MO_u_conv)) ! Initial guess for ZoL
          DO WHILE (Delta.GT.c_accur_sf.AND.n_iter.LT.n_iter_max) 
             Fun = u_star_previter**2*(c_MO_u_conv*u_star_previter-1.)  &
                  + Ri**2*(1.-R_z*c_MO_t_conv*u_star_previter)
             Fun_prime = 3.*c_MO_u_conv*u_star_previter**2              &
                  - 2.*u_star_previter - R_z*c_MO_t_conv*Ri**2
             ZoL = u_star_previter - Fun/Fun_prime
             Delta = ABS(ZoL-u_star_previter)/MAX(c_accur_sf, ABS(ZoL+u_star_previter))
             u_star_previter = ZoL
             n_iter = n_iter + 1
          END DO
       END IF

       !  Compute fetch-dependent Charnock parameter, use "u_star_min_sf"
       CALL lake_SfcFlx_roughness (fetch, U_a, u_star_min_sf, h_ice, c_z0u_fetch, u_star_thresh, z0u_sf, z0t_sf, z0q_sf)

       !  Threshold value of wind speed 
       u_star_st = u_star_thresh
       CALL lake_SfcFlx_roughness (fetch, U_a, u_star_st, h_ice, c_z0u_fetch, u_star_thresh, z0u_sf, z0t_sf, z0q_sf)
       IF(ZoL.GT.0.) THEN   ! MO function in stable stratification 
          psi_u = c_MO_u_stab*ZoL*(1.-MIN(z0u_sf/height_u, 1.))
       ELSE                        ! MO function in convection
          psi_t = (1.-c_MO_u_conv*ZoL)**c_MO_u_exp
          psi_q = (1.-c_MO_u_conv*ZoL*MIN(z0u_sf/height_u, 1.))**c_MO_u_exp
          psi_u = 2.*(ATAN(psi_t)-ATAN(psi_q))                  &
               + 2.*LOG((1.+psi_q)/(1.+psi_t))   &
               + LOG((1.+psi_q*psi_q)/(1.+psi_t*psi_t))   
       END IF
       U_a_thresh = u_star_thresh/c_Karman*(LOG(height_u/z0u_sf)+psi_u)

       !  Compute friction velocity 
       n_iter = 0
       Delta = 1.                ! Set initial error to a large value (as compared to the accuracy)
       u_star_previter = u_star_thresh  ! Initial guess for friction velocity  
       IF(U_a.LE.U_a_thresh) THEN  ! Smooth surface
          DO WHILE (Delta.GT.c_accur_sf.AND.n_iter.LT.n_iter_max) 
             CALL lake_SfcFlx_roughness (fetch, U_a, MIN(u_star_thresh, u_star_previter), h_ice,   &
                  c_z0u_fetch, u_star_thresh, z0u_sf, z0t_sf, z0q_sf)
             IF(ZoL.GE.0.) THEN  ! Stable stratification
                psi_u = c_MO_u_stab*ZoL*(1.-MIN(z0u_sf/height_u, 1.))
                Fun = LOG(height_u/z0u_sf) + psi_u
                Fun_prime = (Fun + 1. + c_MO_u_stab*ZoL*MIN(z0u_sf/height_u, 1.))/c_Karman
                Fun = Fun*u_star_previter/c_Karman - U_a
             ELSE                       ! Convection 
                psi_t = (1.-c_MO_u_conv*ZoL)**c_MO_u_exp
                psi_q = (1.-c_MO_u_conv*ZoL*MIN(z0u_sf/height_u, 1.))**c_MO_u_exp
                psi_u = 2.*(ATAN(psi_t)-ATAN(psi_q))                  &
                     + 2.*LOG((1.+psi_q)/(1.+psi_t))   &
                     + LOG((1.+psi_q*psi_q)/(1.+psi_t*psi_t))   
                Fun = LOG(height_u/z0u_sf) + psi_u
                Fun_prime = (Fun + 1./psi_q)/c_Karman
                Fun = Fun*u_star_previter/c_Karman - U_a
             END IF
             u_star_st = u_star_previter - Fun/Fun_prime
             Delta = ABS((u_star_st-u_star_previter)/(u_star_st+u_star_previter))
             u_star_previter = u_star_st
             n_iter = n_iter + 1
          END DO
       ELSE                        ! Rough surface
          DO WHILE (Delta.GT.c_accur_sf.AND.n_iter.LT.n_iter_max) 
             CALL lake_SfcFlx_roughness (fetch, U_a, MAX(u_star_thresh, u_star_previter), h_ice,   &
                  c_z0u_fetch, u_star_thresh, z0u_sf, z0t_sf, z0q_sf)
             IF(ZoL.GE.0.) THEN  ! Stable stratification
                psi_u = c_MO_u_stab*ZoL*(1.-MIN(z0u_sf/height_u, 1.))
                Fun = LOG(height_u/z0u_sf) + psi_u
                Fun_prime = (Fun - 2. - 2.*c_MO_u_stab*ZoL*MIN(z0u_sf/height_u, 1.))/c_Karman
                Fun = Fun*u_star_previter/c_Karman - U_a
             ELSE                       ! Convection 
                psi_t = (1.-c_MO_u_conv*ZoL)**c_MO_u_exp
                psi_q = (1.-c_MO_u_conv*ZoL*MIN(z0u_sf/height_u, 1.))**c_MO_u_exp
                psi_u = 2.*(ATAN(psi_t)-ATAN(psi_q))                  &
                     + 2.*LOG((1.+psi_q)/(1.+psi_t))   &
                     + LOG((1.+psi_q*psi_q)/(1.+psi_t*psi_t))   
                Fun = LOG(height_u/z0u_sf) + psi_u
                Fun_prime = (Fun - 2./psi_q)/c_Karman
                Fun = Fun*u_star_previter/c_Karman - U_a
             END IF
             IF(h_ice.GE.h_Ice_min_flk) THEN   ! No iteration is required for rough flow over ice
                u_star_st = c_Karman*U_a/MAX(c_small_sf, LOG(height_u/z0u_sf)+psi_u)
                u_star_previter = u_star_st
             ELSE                              ! Iterate in case of open water
                u_star_st = u_star_previter - Fun/Fun_prime
             END IF
             Delta = ABS((u_star_st-u_star_previter)/(u_star_st+u_star_previter))
             u_star_previter = u_star_st
             n_iter = n_iter + 1
          END DO
       END IF


       !  Momentum flux
       Q_mom_tur = -u_star_st*u_star_st

       !  Temperature and specific humidity fluxes
       CALL lake_SfcFlx_roughness (fetch, U_a, u_star_st, h_ice, c_z0u_fetch, u_star_thresh, z0u_sf, z0t_sf, z0q_sf)

       IF(ZoL.GE.0.) THEN   ! Stable stratification 
          psi_t = c_MO_t_stab*R_z*ZoL*(1.-MIN(z0t_sf/height_tq, 1.))
          psi_q = c_MO_q_stab*R_z*ZoL*(1.-MIN(z0q_sf/height_tq, 1.))
       ELSE                        ! Convection 
          psi_u = (1.-c_MO_t_conv*R_z*ZoL)**c_MO_t_exp
          psi_t = (1.-c_MO_t_conv*R_z*ZoL*MIN(z0t_sf/height_tq, 1.))**c_MO_t_exp
          psi_t = 2.*LOG((1.+psi_t)/(1.+psi_u))
          psi_u = (1.-c_MO_q_conv*R_z*ZoL)**c_MO_q_exp
          psi_q = (1.-c_MO_q_conv*R_z*ZoL*MIN(z0q_sf/height_tq, 1.))**c_MO_q_exp
          psi_q = 2.*LOG((1.+psi_q)/(1.+psi_u))

       END IF
       Q_sen_tur = -(T_a-T_s)*u_star_st*c_Karman/Pr_neutral  &
            / MAX(c_small_sf, LOG(height_tq/z0t_sf)+psi_t)
       Q_lat_tur = -(q_a-q_s)*u_star_st*c_Karman/Sc_neutral  &
            / MAX(c_small_sf, LOG(height_tq/z0q_sf)+psi_q)

    END IF Turb_Fluxes

    !4.  Decide between turbulent, molecular, and convective fluxes
    Q_momentum = MIN(Q_mom_tur, Q_mom_mol, Q_mom_con)  ! Momentum flux is negative  
    IF(l_conv_visc) THEN    ! Convection, take fluxes that are maximal in magnitude 
       IF(ABS(Q_sen_tur).GE.ABS(Q_sen_con)) THEN
          regime_flux_lake_r=1
          Q_sensible = Q_sen_tur
       ELSE
          Q_sensible = Q_sen_con
          regime_flux_lake_r=2
       END IF
       IF(ABS(Q_sensible).LT.ABS(Q_sen_mol)) THEN
          Q_sensible = Q_sen_mol
          regime_flux_lake_r=3
       END IF
       IF(ABS(Q_lat_tur).GE.ABS(Q_lat_con)) THEN
          Q_latent = Q_lat_tur
          regime_flux_lake_temp=10
       ELSE
          Q_latent = Q_lat_con
          regime_flux_lake_temp=20
       END IF
       IF(ABS(Q_latent).LT.ABS(Q_lat_mol)) THEN
          Q_latent = Q_lat_mol
          regime_flux_lake_temp=30
       END IF
    ELSE                    ! Stable or neutral stratification, chose fluxes that are maximal in magnitude 
       IF(ABS(Q_sen_tur).GE.ABS(Q_sen_mol)) THEN 
          Q_sensible = Q_sen_tur
          regime_flux_lake_r=1
       ELSE 
          Q_sensible = Q_sen_mol    
          regime_flux_lake_r=3
       END IF
       IF(ABS(Q_lat_tur).GE.ABS(Q_lat_mol)) THEN 
          Q_latent = Q_lat_tur
          regime_flux_lake_temp=10
       ELSE 
          Q_latent = Q_lat_mol  
          regime_flux_lake_temp=30
       END IF
    END IF

    regime_flux_lake_r=regime_flux_lake_r+ regime_flux_lake_temp

    !5.  Set output (notice that fluxes are no longer in kinematic units)

    Q_momentum = Q_momentum*rho_a 
    Q_sensible = Q_sensible*rho_a*tpsf_c_a_p
    Q_watvap   = Q_latent*rho_a
    Q_latent = tpsf_L_evap
    IF(h_ice.GE.h_Ice_min_flk) Q_latent = Q_latent + tpl_L_f   ! Add latent heat of fusion over ice
    Q_latent = Q_watvap*Q_latent


  END SUBROUTINE lake_SfcFlx_momsenlat


  !! ================================================================================================================================
  !! SUBROUTINE 	: lake_driver
  !!
  !>\BRIEF          The main driving routine of the lake model FLake 
  !!                where computations are performed.
  !!                Advances the surface temperature
  !!                and other FLake variables one time step.
  !!                At the moment, the Euler explicit scheme is used.
  !!
  !! DESCRIPTION :  The main driving routine of the lake model FLake 
  !!                where computations are performed.
  !!                Advances the surface temperature
  !!                and other FLake variables one time step.
  !!                At the moment, the Euler explicit scheme is used.
  !!
  !! MAIN OUTPUT VARIABLE(S) : T_sfc_n,  T_snow_n_flk, T_ice_n_flk,  T_wML_n_flk, 
  !!                           T_mnw_n_flk, T_bot_n_flk, T_B1_n_flk, h_snow_n_flk, 
  !!                           h_ice_n_flk, h_ML_n_flk, H_B1_n_flk, C_T_n_flk
  !!
  !! REFERENCE(S) : 
  !!
  !! FLOWCHART    : None
  !! \n
  !_ ================================================================================================================================


  SUBROUTINE lake_driver  ( depth_w, depth_sed_lake, t_sed_lake, par_Coriolis,           &
       extincoef_water_typ,                             &
       del_time,                                        &
       dMsnowdt_flk, Q_w_flk, I_snow_flk, I_ice_flk, &
       I_w_flk, I_h_flk, I_intm_0_h_flk, I_intm_h_D_flk, &
       I_bot_flk, Q_snow_flk,  Q_ice_flk, Q_bot_flk, u_star_w_flk,&
       T_sfc_p, T_snow_p_flk, T_ice_p_flk, T_wML_p_flk, &
       T_mnw_p_flk, T_bot_p_flk, T_B1_p_flk, h_snow_p_flk, &
       h_ice_p_flk, h_ML_p_flk, H_B1_p_flk, C_T_p_flk, &
       T_sfc_n,  T_snow_n_flk, T_ice_n_flk,  T_wML_n_flk, &
       T_mnw_n_flk, T_bot_n_flk, T_B1_n_flk, h_snow_n_flk, &
       h_ice_n_flk, h_ML_n_flk, H_B1_n_flk, C_T_n_flk)


    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    REAL (r_std), INTENT(IN) :: del_time                          !! The model time step [s]

    !! 0.1.1 Atmospheric variables 
    REAL (r_std), INTENT(IN) ::  dMsnowdt_flk

    !! 0.1.2 FLake parameters
    REAL (r_std), INTENT(IN) :: depth_w                           !! The lake depth [m]
    REAL (r_std), INTENT(IN) :: depth_sed_lake                    !! Depth of the thermally active layer of bottom sediments [m]
    REAL (r_std), INTENT(IN) :: t_sed_lake                        !! Temperature at the outer edge of 
                                                                  !! the thermally active layer of bottom sediments [K]
    REAL (r_std), INTENT(IN) :: extincoef_water_typ               !! "Typical" extinction coefficient of the lake water [m^{-1}],
                                                                  !! used to compute the equilibrium CBL depth
    REAL (r_std), INTENT(IN) :: par_Coriolis                      !! The Coriolis parameter [s^{-1}]

    !! 0.1.3 FLake state variables
    REAL (r_std), INTENT(IN) :: T_sfc_p
    REAL (r_std), INTENT(IN) :: T_snow_p_flk 
    REAL (r_std), INTENT(IN) :: T_ice_p_flk
    REAL (r_std), INTENT(IN) :: T_wML_p_flk
    REAL (r_std), INTENT(IN) :: T_mnw_p_flk 
    REAL (r_std), INTENT(IN) :: T_bot_p_flk
    REAL (r_std), INTENT(IN) :: T_B1_p_flk 
    REAL (r_std), INTENT(IN) :: h_snow_p_flk
    REAL (r_std), INTENT(IN) :: h_ice_p_flk
    REAL (r_std), INTENT(IN) :: h_ML_p_flk
    REAL (r_std), INTENT(IN) :: H_B1_p_flk 
    REAL (r_std), INTENT(IN) :: C_T_p_flk

    !! 0.2 Output variables
    REAL (r_std), INTENT(INOUT) :: T_sfc_n
    REAL (r_std), INTENT(INOUT) :: T_snow_n_flk
    REAL (r_std), INTENT(INOUT) :: T_ice_n_flk 
    REAL (r_std), INTENT(INOUT) :: T_wML_n_flk
    REAL (r_std), INTENT(INOUT) :: T_mnw_n_flk
    REAL (r_std), INTENT(INOUT) :: T_bot_n_flk 
    REAL (r_std), INTENT(INOUT) :: T_B1_n_flk
    REAL (r_std), INTENT(INOUT) :: h_snow_n_flk
    REAL (r_std), INTENT(INOUT) :: h_ice_n_flk 
    REAL (r_std), INTENT(INOUT) :: h_ML_n_flk
    REAL (r_std), INTENT(INOUT) :: H_B1_n_flk
    REAL (r_std), INTENT(INOUT) :: C_T_n_flk


    !! 0.4 Local variables
    LOGICAL ::          &
         l_ice_create    , & ! Switch, .TRUE. = ice does not exist but should be created
         l_snow_exists   , & ! Switch, .TRUE. = there is snow above the ice
         l_ice_meltabove     ! Switch, .TRUE. = snow/ice melting from above takes place

    INTEGER (i_std) :: &
         i                             ! Loop index

    REAL (r_std) ::    &
         d_T_mnw_dt             , & ! Time derivative of T_mnw [K s^{-1}] 
         d_T_ice_dt             , & ! Time derivative of T_ice [K s^{-1}] 
         d_T_bot_dt             , & ! Time derivative of T_bot [K s^{-1}] 
         d_T_B1_dt              , & ! Time derivative of T_B1 [K s^{-1}] 
         d_h_snow_dt            , & ! Time derivative of h_snow [m s^{-1}]
         d_h_ice_dt             , & ! Time derivative of h_ice [m s^{-1}]
         d_h_ML_dt              , & ! Time derivative of h_ML [m s^{-1}]
         d_H_B1_dt              , & ! Time derivative of H_B1 [m s^{-1}]
         d_C_T_dt                   ! Time derivative of C_T [s^{-1}]

    REAL (r_std) ::    &
         N_T_mean               , & ! The mean buoyancy frequency in the thermocline [s^{-1}] 
         ZM_h_scale             , & ! The ZM96 equilibrium SBL depth scale [m] 
         conv_equil_h_scale         ! The equilibrium CBL depth scale [m]

    REAL (r_std) :: &
         h_ice_threshold     , & ! If h_ice<h_ice_threshold, use quasi-equilibrium ice model 
         flk_str_1           , & ! Help storage variable
         flk_str_2           , & ! Help storage variable
         R_H_icesnow         , & ! Dimensionless ratio, used to store intermediate results
         R_rho_c_icesnow     , & ! Dimensionless ratio, used to store intermediate results
         R_TI_icesnow        , & ! Dimensionless ratio, used to store intermediate results
         R_Tstar_icesnow         ! Dimensionless ratio, used to store intermediate results

    !  The shape factor(s) at the previous time step ("p") and the updated value(s) ("n") 
    REAL (r_std) ::           &
         C_TT_flk                      , & ! Dimensionless parameter (thermocline)
         C_Q_flk                       , & ! Shape factor with respect to the heat flux (thermocline)
         C_I_flk                       , & ! Shape factor (ice)
         C_S_flk                           ! Shape factor (snow)

    !  Derivatives of the shape functions
    REAL (r_std) ::           &
         Phi_T_pr0_flk                 , & ! d\Phi_T(0)/d\zeta   (thermocline)
         Phi_I_pr0_flk                 , & ! d\Phi_I(0)/d\zeta_I (ice)
         Phi_I_pr1_flk                 , & ! d\Phi_I(1)/d\zeta_I (ice)
         Phi_S_pr0_flk                     ! d\Phi_S(0)/d\zeta_S (snow)


    !  Velocity scales
    REAL (r_std), INTENT(IN) :: u_star_w_flk ! Friction velocity in the surface layer of lake water [m s^{-1}]
    REAL (r_std) ::  w_star_sfc_flk                    ! Convective velocity scale, 
    ! using a generalized heat flux scale [m s^{-1}]

    REAL (r_std), INTENT(IN) :: I_snow_flk, I_ice_flk, I_w_flk, I_h_flk, I_intm_0_h_flk, I_intm_h_D_flk, I_bot_flk
    REAL (r_std), INTENT(INOUT) :: Q_w_flk, Q_snow_flk,  Q_ice_flk,   Q_bot_flk                     
    REAL (r_std)  Q_star_flk                   ! A generalized heat flux scale [W m^{-2}]



    REAL (r_std) :: dMlim
    INTEGER (i_std) :: Nb_stab_snow_flk 
    REAL (r_std) :: del_time_stab_snow 
    INTEGER (i_std) :: n_stab_snow_flk 



    !_dm 
    ! Security. Set time-rate-of-change of prognostic variables to zero.
    ! Set prognostic variables to their values at the previous time step.
    ! (This is to avoid spurious changes of prognostic variables 
    ! when FLake is used within a 3D model, e.g. to avoid spurious generation of ice 
    ! at the neighbouring lake points as noticed by Burkhardt Rockel.)
    !_dm 

    d_T_mnw_dt   = 0. 
    d_T_ice_dt   = 0. 
    d_T_bot_dt   = 0. 
    d_T_B1_dt    = 0. 
    d_h_snow_dt  = 0. 
    d_h_ice_dt   = 0. 
    d_h_ML_dt    = 0. 
    d_H_B1_dt    = 0. 
    d_C_T_dt     = 0.    

    !------------------------------------------------------------------------------
    !  Compute fluxes, using variables from the previous time step.
    !------------------------------------------------------------------------------

    !_dm
    ! At this point, the heat and radiation fluxes, namely,
    ! Q_snow_flk, Q_ice_flk, Q_w_flk, 
    ! I_atm_flk, I_snow_flk, I_ice_flk, I_w_flk, I_h_flk, I_bot_flk,     
    ! the mean radiation flux over the mixed layer, I_intm_0_h_flk, 
    ! and the mean radiation flux over the thermocline, I_intm_h_D_flk, 
    ! should be known.
    ! They are computed within "flake_interface" (or within the driving model)
    ! and are available to "flake_driver"
    ! through the above variables declared in the MODULE "flake".
    ! In case a lake is ice-covered, Q_w_flk is re-computed below.
    !_dm

    ! Heat flux through the ice-water interface
    IF(h_ice_p_flk.GE.h_Ice_min_flk) THEN    ! Ice exists 
       IF(h_ML_p_flk.LE.h_ML_min_flk) THEN    ! Mixed-layer depth is zero, compute flux 
          Q_w_flk = -tpl_kappa_w*(T_bot_p_flk-T_wML_p_flk)/depth_w  ! Flux with linear T(z) 
          Phi_T_pr0_flk = Phi_T_pr0_1*C_T_p_flk-Phi_T_pr0_2         ! d\Phi(0)/d\zeta (thermocline)
          Q_w_flk = Q_w_flk*MAX(Phi_T_pr0_flk, 1.)           ! Account for an increased d\Phi(0)/d\zeta 
       ELSE                    
          Q_w_flk = 0.                  ! Mixed-layer depth is greater than zero, set flux to zero
       END IF
    END IF

    ! A generalized heat flux scale 
    Q_star_flk = Q_w_flk + I_w_flk + I_h_flk - 2.*I_intm_0_h_flk

    ! Heat flux through the water-bottom sediment interface
    IF (ok_botsed_lake) THEN
       Q_bot_flk = -tpl_kappa_w*(T_B1_p_flk-T_bot_p_flk)/MAX(H_B1_p_flk, H_B1_min_flk)*Phi_B1_pr0
    ELSE  
       Q_bot_flk = 0.   ! The bottom-sediment scheme is not used
    END IF


    !------------------------------------------------------------------------------
    !  Check if ice exists or should be created.
    !  If so, compute the thickness and the temperature of ice and snow.
    !------------------------------------------------------------------------------

    !_dm
    ! Notice that a quasi-equilibrium ice-snow model is used 
    ! to avoid numerical instability when the ice is thin.
    ! This is always the case when new ice is created.
    !_dm

    !_dev
    ! The dependence of snow density and of snow heat conductivity 
    ! on the snow thickness is accounted for parametrically.
    ! That is, the time derivatives of \rho_S and \kappa_S are neglected.
    ! The exception is the equation for the snow thickness 
    ! in case of snow accumulation and no melting, 
    ! where d\rho_S/dt is incorporated.
    ! Furthermore, some (presumably small) correction terms incorporating 
    ! the snow density and the snow heat conductivity are dropped out.
    ! Those terms may be included as better formulations 
    ! for \rho_S and \kappa_S are available.
    !_dev

    ! Default values
    l_ice_create    = .FALSE.  
    l_ice_meltabove = .FALSE.  

    Phi_I_pr1_flk = Phi_I_pr1_lin + Phi_I_ast_MR*MIN(1., h_ice_n_flk/H_Ice_max)
    R_H_icesnow = Phi_I_pr1_flk/Phi_S_pr0_lin*tpl_kappa_I/lake_snowheatconduct(h_snow_n_flk) &
         *h_snow_n_flk/MAX(h_ice_n_flk, h_Ice_min_flk)

    !dMlim = 0.005

    Nb_stab_snow_flk = 1 !INT(R_H_icesnow)+1+3*INT( 0.03/(h_ice_n_flk+0.0015))   !INT(dMsnowdt_flk/dMlim)+1
    del_time_stab_snow = del_time/Nb_stab_snow_flk  

    Ice_exist: IF(h_ice_p_flk.LT.h_Ice_min_flk) THEN   ! Ice does not exist 

       l_ice_create = T_wML_p_flk.LE.(tpl_T_f+c_small_flk).AND.Q_w_flk.LT.0.

       IF(l_ice_create) THEN                            ! Ice does not exist but should be created
          DO n_stab_snow_flk = 1, Nb_stab_snow_flk
             d_h_ice_dt = -Q_w_flk/tpl_rho_I/tpl_L_f                                  
             h_ice_n_flk = h_ice_n_flk + d_h_ice_dt*del_time_stab_snow                          ! Advance h_ice 
             T_ice_n_flk = tpl_T_f + h_ice_n_flk*Q_w_flk/tpl_kappa_I/Phi_I_pr0_lin    ! Ice temperature
             d_h_snow_dt = dMsnowdt_flk/tpl_rho_S_min 
             h_snow_n_flk = h_snow_n_flk + d_h_snow_dt*del_time_stab_snow                       ! Advance h_snow
             Phi_I_pr1_flk = Phi_I_pr1_lin                                    & 
                  + Phi_I_ast_MR*MIN(1., h_ice_n_flk/H_Ice_max)       ! d\Phi_I(1)/d\zeta_I (ice)
             R_H_icesnow = Phi_I_pr1_flk/Phi_S_pr0_lin*tpl_kappa_I/lake_snowheatconduct(h_snow_n_flk) &
                  * h_snow_n_flk/MAX(h_ice_n_flk, h_Ice_min_flk)
             T_snow_n_flk = T_ice_n_flk + R_H_icesnow*(T_ice_n_flk-tpl_T_f)           ! Snow temperature
          END DO
       END IF

    ELSE Ice_exist                                     ! Ice exists

       DO n_stab_snow_flk = 1, Nb_stab_snow_flk

          l_snow_exists = h_snow_n_flk.GE.h_Snow_min_flk   ! Check if there is snow above the ice

          Melting: IF(T_snow_n_flk.GE.(tpl_T_f-c_small_flk)) THEN  ! T_sfc = T_f, check for melting from above
             ! T_snow = T_ice if snow is absent 
             IF(l_snow_exists) THEN   ! There is snow above the ice
                flk_str_1 = Q_snow_flk + I_snow_flk - I_ice_flk        ! Atmospheric forcing
                IF(flk_str_1.GE.0.) THEN  ! Melting of snow and ice from above
                   l_ice_meltabove = .TRUE.
                   d_h_snow_dt = (-flk_str_1/tpl_L_f+dMsnowdt_flk)/lake_snowdensity(h_snow_n_flk)
                   d_h_ice_dt  = -(I_ice_flk - I_w_flk - Q_w_flk)/tpl_L_f/tpl_rho_I 
                END IF
             ELSE                     ! No snow above the ice
                flk_str_1 = Q_ice_flk + I_ice_flk - I_w_flk - Q_w_flk  ! Atmospheric forcing + heating from the water
                IF(flk_str_1.GE.0.) THEN  ! Melting of ice from above, snow accumulation may occur
                   l_ice_meltabove = .TRUE.
                   d_h_ice_dt  = -flk_str_1/tpl_L_f/tpl_rho_I 
                   d_h_snow_dt = dMsnowdt_flk/tpl_rho_S_min
                END IF
             END IF
             IF(l_ice_meltabove) THEN  ! Melting from above takes place
                h_ice_n_flk  = h_ice_n_flk  + d_h_ice_dt *del_time_stab_snow  ! Advance h_ice
                h_snow_n_flk = h_snow_n_flk + d_h_snow_dt*del_time_stab_snow  ! Advance h_snow
                T_ice_n_flk  = tpl_T_f                              ! Set T_ice to the freezing point
                T_snow_n_flk = tpl_T_f                              ! Set T_snow to the freezing point
             END IF


          END IF Melting

          No_Melting: IF(.NOT.l_ice_meltabove) THEN                 ! No melting from above

             d_h_snow_dt = lake_snowdensity(h_snow_n_flk)  
             IF(d_h_snow_dt.LT.tpl_rho_S_max) THEN    ! Account for d\rho_S/dt
                flk_str_1 = h_snow_n_flk*tpl_Gamma_rho_S/tpl_rho_w_r
                flk_str_1 = flk_str_1/(1.-flk_str_1)
             ELSE                                     ! Snow density is equal to its maximum value, d\rho_S/dt=0
                flk_str_1 = 0.
             END IF
             d_h_snow_dt = dMsnowdt_flk/d_h_snow_dt/(1.+flk_str_1)       ! Snow accumulation
             h_snow_n_flk = h_snow_n_flk + d_h_snow_dt*del_time_stab_snow                         ! Advance h_snow

             Phi_I_pr0_flk = h_ice_n_flk/H_Ice_max                              ! h_ice relative to its maximum value
             C_I_flk = C_I_lin  - C_I_MR*(1.+Phi_I_ast_MR)*Phi_I_pr0_flk  ! Shape factor (ice)
             Phi_I_pr1_flk = Phi_I_pr1_lin + Phi_I_ast_MR*Phi_I_pr0_flk         ! d\Phi_I(1)/d\zeta_I (ice)
             Phi_I_pr0_flk = Phi_I_pr0_lin - Phi_I_pr0_flk                      ! d\Phi_I(0)/d\zeta_I (ice)

             h_ice_threshold = MAX(1., 2.*C_I_flk*tpl_c_I*(tpl_T_f-T_ice_n_flk)/tpl_L_f)
             h_ice_threshold = Phi_I_pr0_flk/C_I_flk*tpl_kappa_I/tpl_rho_I/tpl_c_I*h_ice_threshold
             h_ice_threshold = SQRT(h_ice_threshold*del_time_stab_snow)                   ! Threshold value of h_ice
             h_ice_threshold = MIN(0.9*H_Ice_max, MAX(h_ice_threshold, h_Ice_min_flk))
             ! h_ice(threshold) < 0.9*H_Ice_max

             IF(h_ice_n_flk.LT.h_ice_threshold) THEN  ! Use a quasi-equilibrium ice model

                IF(l_snow_exists) THEN   ! Use fluxes at the air-snow interface
                   flk_str_1 = Q_snow_flk + I_snow_flk - I_w_flk
                ELSE                     ! Use fluxes at the air-ice interface
                   flk_str_1 = Q_ice_flk + I_ice_flk - I_w_flk
                END IF
                d_h_ice_dt = -(flk_str_1-Q_w_flk)/tpl_L_f/tpl_rho_I
                h_ice_n_flk = h_ice_n_flk + d_h_ice_dt *del_time_stab_snow                         ! Advance h_ice
                T_ice_n_flk = tpl_T_f + h_ice_n_flk*flk_str_1/tpl_kappa_I/Phi_I_pr0_flk  ! Ice temperature

             ELSE                                     ! Use a complete ice model

                d_h_ice_dt = tpl_kappa_I*(tpl_T_f-T_ice_n_flk)/h_ice_n_flk*Phi_I_pr0_flk
                d_h_ice_dt = (Q_w_flk+d_h_ice_dt)/tpl_L_f/tpl_rho_I
                h_ice_n_flk = h_ice_n_flk  + d_h_ice_dt*del_time_stab_snow                         ! Advance h_ice

                R_TI_icesnow = tpl_c_I*(tpl_T_f-T_ice_n_flk)/tpl_L_f         ! Dimensionless parameter
                R_Tstar_icesnow = 1. - C_I_flk                        ! Dimensionless parameter
                IF(l_snow_exists) THEN  ! There is snow above the ice
                   R_H_icesnow = Phi_I_pr1_flk/Phi_S_pr0_lin*tpl_kappa_I/lake_snowheatconduct(h_snow_n_flk) &
                        * h_snow_n_flk/h_ice_n_flk
                   R_rho_c_icesnow = lake_snowdensity(h_snow_n_flk)*tpl_c_S/tpl_rho_I/tpl_c_I 

                   R_Tstar_icesnow = R_Tstar_icesnow*R_TI_icesnow             ! Dimensionless parameter


                   flk_str_2 = Q_snow_flk+I_snow_flk-I_w_flk                  ! Atmospheric fluxes
                   flk_str_1  = C_I_flk*h_ice_n_flk + (1.+C_S_lin*R_H_icesnow)*R_rho_c_icesnow*h_snow_p_flk
                   d_T_ice_dt = -(1.-2.*C_S_lin)*R_H_icesnow*(tpl_T_f-T_ice_n_flk)             & 
                        * tpl_c_S*dMsnowdt_flk                          ! Effect of snow accumulation
                ELSE                    ! No snow above the ice
                   R_Tstar_icesnow = R_Tstar_icesnow*R_TI_icesnow             ! Dimensionless parameter
                   flk_str_2 = Q_ice_flk+I_ice_flk-I_w_flk                    ! Atmospheric fluxes
                   flk_str_1  = C_I_flk*h_ice_n_flk
                   d_T_ice_dt = 0.
                END IF
                d_T_ice_dt = d_T_ice_dt + tpl_kappa_I*(tpl_T_f-T_ice_n_flk)/h_ice_n_flk*Phi_I_pr0_flk       &
                     * (1.-R_Tstar_icesnow)                     ! Add flux due to heat conduction
                d_T_ice_dt = d_T_ice_dt - R_Tstar_icesnow*Q_w_flk            ! Add flux from water to ice
                d_T_ice_dt = d_T_ice_dt + flk_str_2                          ! Add atmospheric fluxes
                d_T_ice_dt = d_T_ice_dt/tpl_rho_I/tpl_c_I                    ! Total forcing
                d_T_ice_dt = d_T_ice_dt/flk_str_1                            ! dT_ice/dt 
                T_ice_n_flk = T_ice_n_flk + d_T_ice_dt*del_time_stab_snow                          ! Advance T_ice
             END IF

             Phi_I_pr1_flk = MIN(1., h_ice_n_flk/H_Ice_max)          ! h_ice relative to its maximum value
             Phi_I_pr1_flk = Phi_I_pr1_lin + Phi_I_ast_MR*Phi_I_pr1_flk     ! d\Phi_I(1)/d\zeta_I (ice)
             R_H_icesnow = Phi_I_pr1_flk/Phi_S_pr0_lin*tpl_kappa_I/lake_snowheatconduct(h_snow_n_flk) &
                  *h_snow_n_flk/MAX(h_ice_n_flk, h_Ice_min_flk)
             T_snow_n_flk = T_ice_n_flk + R_H_icesnow*(T_ice_n_flk-tpl_T_f)             ! Snow temperature

          END IF No_Melting
       END DO

       !  l_snow_exists = h_snow_p_flk.GE.h_Snow_min_flk   ! Check if there is snow above the ice
       !
       !  Melting: IF(T_snow_p_flk.GE.(tpl_T_f-c_small_flk)) THEN  ! T_sfc = T_f, check for melting from above
       !                                                           ! T_snow = T_ice if snow is absent 
       !    IF(l_snow_exists) THEN   ! There is snow above the ice
       !      flk_str_1 = Q_snow_flk + I_snow_flk - I_ice_flk        ! Atmospheric forcing
       !      IF(flk_str_1.GE.0.) THEN  ! Melting of snow and ice from above
       !        l_ice_meltabove = .TRUE.
       !        d_h_snow_dt = (-flk_str_1/tpl_L_f+dMsnowdt_flk)/lake_snowdensity(h_snow_p_flk)
       !        d_h_ice_dt  = -(I_ice_flk - I_w_flk - Q_w_flk)/tpl_L_f/tpl_rho_I 
       !      END IF 
       !    ELSE                     ! No snow above the ice
       !      flk_str_1 = Q_ice_flk + I_ice_flk - I_w_flk - Q_w_flk  ! Atmospheric forcing + heating from the water
       !      IF(flk_str_1.GE.0.) THEN  ! Melting of ice from above, snow accumulation may occur
       !        l_ice_meltabove = .TRUE.
       !        d_h_ice_dt  = -flk_str_1/tpl_L_f/tpl_rho_I 
       !        d_h_snow_dt = dMsnowdt_flk/tpl_rho_S_min
       !      END IF 
       !    END IF 
       !    IF(l_ice_meltabove) THEN  ! Melting from above takes place
       !      h_ice_n_flk  = h_ice_p_flk  + d_h_ice_dt *del_time  ! Advance h_ice
       !      h_snow_n_flk = h_snow_p_flk + d_h_snow_dt*del_time  ! Advance h_snow
       !      T_ice_n_flk  = tpl_T_f                              ! Set T_ice to the freezing point
       !      T_snow_n_flk = tpl_T_f                              ! Set T_snow to the freezing point
       !    END IF
       !
       !  END IF Melting
       !
       !  No_Melting: IF(.NOT.l_ice_meltabove) THEN                 ! No melting from above
       !
       !    d_h_snow_dt = lake_snowdensity(h_snow_p_flk)  
       !    IF(d_h_snow_dt.LT.tpl_rho_S_max) THEN    ! Account for d\rho_S/dt
       !     flk_str_1 = h_snow_p_flk*tpl_Gamma_rho_S/tpl_rho_w_r
       !     flk_str_1 = flk_str_1/(1.-flk_str_1)
       !    ELSE                                     ! Snow density is equal to its maximum value, d\rho_S/dt=0
       !     flk_str_1 = 0.
       !    END IF
       !    d_h_snow_dt = dMsnowdt_flk/d_h_snow_dt/(1.+flk_str_1)       ! Snow accumulation
       !    h_snow_n_flk = h_snow_p_flk + d_h_snow_dt*del_time                         ! Advance h_snow
       !    
       !    Phi_I_pr0_flk = h_ice_p_flk/H_Ice_max                              ! h_ice relative to its maximum value
       !    C_I_flk = C_I_lin - C_I_MR*(1.+Phi_I_ast_MR)*Phi_I_pr0_flk  ! Shape factor (ice)
       !    Phi_I_pr1_flk = Phi_I_pr1_lin + Phi_I_ast_MR*Phi_I_pr0_flk         ! d\Phi_I(1)/d\zeta_I (ice)
       !    Phi_I_pr0_flk = Phi_I_pr0_lin - Phi_I_pr0_flk                      ! d\Phi_I(0)/d\zeta_I (ice)
       !
       !    h_ice_threshold = MAX(1., 2.*C_I_flk*tpl_c_I*(tpl_T_f-T_ice_p_flk)/tpl_L_f)
       !    h_ice_threshold = Phi_I_pr0_flk/C_I_flk*tpl_kappa_I/tpl_rho_I/tpl_c_I*h_ice_threshold
       !    h_ice_threshold = SQRT(h_ice_threshold*del_time)                   ! Threshold value of h_ice
       !    h_ice_threshold = MIN(0.9*H_Ice_max, MAX(h_ice_threshold, h_Ice_min_flk))
       !                                                                       ! h_ice(threshold) < 0.9*H_Ice_max
       !
       !    IF(h_ice_p_flk.LT.h_ice_threshold) THEN  ! Use a quasi-equilibrium ice model
       !
       !      IF(l_snow_exists) THEN   ! Use fluxes at the air-snow interface
       !        flk_str_1 = Q_snow_flk + I_snow_flk - I_w_flk
       !      ELSE                     ! Use fluxes at the air-ice interface
       !        flk_str_1 = Q_ice_flk + I_ice_flk - I_w_flk
       !      END IF
       !      d_h_ice_dt = -(flk_str_1-Q_w_flk)/tpl_L_f/tpl_rho_I
       !      h_ice_n_flk = h_ice_p_flk + d_h_ice_dt *del_time                         ! Advance h_ice
       !      T_ice_n_flk = tpl_T_f + h_ice_n_flk*flk_str_1/tpl_kappa_I/Phi_I_pr0_flk  ! Ice temperature
       !
       !    ELSE                                     ! Use a complete ice model
       !
       !      d_h_ice_dt = tpl_kappa_I*(tpl_T_f-T_ice_p_flk)/h_ice_p_flk*Phi_I_pr0_flk
       !      d_h_ice_dt = (Q_w_flk+d_h_ice_dt)/tpl_L_f/tpl_rho_I
       !      h_ice_n_flk = h_ice_p_flk  + d_h_ice_dt*del_time                         ! Advance h_ice
       !
       !      R_TI_icesnow = tpl_c_I*(tpl_T_f-T_ice_p_flk)/tpl_L_f         ! Dimensionless parameter
       !      R_Tstar_icesnow = 1. - C_I_flk                        ! Dimensionless parameter
       !      IF(l_snow_exists) THEN  ! There is snow above the ice
       !        R_H_icesnow = Phi_I_pr1_flk/Phi_S_pr0_lin*tpl_kappa_I/lake_snowheatconduct(h_snow_p_flk) &
       !                    * h_snow_p_flk/h_ice_p_flk
       !        R_rho_c_icesnow = lake_snowdensity(h_snow_p_flk)*tpl_c_S/tpl_rho_I/tpl_c_I 
       !!_dev 
       !!_dm 
       !! These terms should be included as an improved understanding of the snow scheme is gained, 
       !! of the effect of snow density in particular. 
       !!_dm 
       !!_nu        R_Tstar_icesnow = R_Tstar_icesnow                                                           &
       !!_nu                        + (1._ireals+C_S_lin*h_snow_p_flk/h_ice_p_flk)*R_H_icesnow*R_rho_c_icesnow
       !!_dev
       !
       !        R_Tstar_icesnow = R_Tstar_icesnow*R_TI_icesnow             ! Dimensionless parameter
       !
       !!_dev
       !!_nu        R_Tstar_icesnow = R_Tstar_icesnow                                                         &
       !!_nu                        + (1._ireals-R_rho_c_icesnow)*tpl_c_I*T_ice_p_flk/tpl_L_f
       !!_dev
       !        flk_str_2 = Q_snow_flk+I_snow_flk-I_w_flk                  ! Atmospheric fluxes
       !        flk_str_1  = C_I_flk*h_ice_p_flk + (1.+C_S_lin*R_H_icesnow)*R_rho_c_icesnow*h_snow_p_flk
       !        d_T_ice_dt = -(1.-2.*C_S_lin)*R_H_icesnow*(tpl_T_f-T_ice_p_flk)             & 
       !                   * tpl_c_S*dMsnowdt_flk                          ! Effect of snow accumulation
       !      ELSE                    ! No snow above the ice
       !        R_Tstar_icesnow = R_Tstar_icesnow*R_TI_icesnow             ! Dimensionless parameter
       !        flk_str_2 = Q_ice_flk+I_ice_flk-I_w_flk                    ! Atmospheric fluxes
       !        flk_str_1  = C_I_flk*h_ice_p_flk
       !        d_T_ice_dt = 0.
       !      END IF 
       !      d_T_ice_dt = d_T_ice_dt + tpl_kappa_I*(tpl_T_f-T_ice_p_flk)/h_ice_p_flk*Phi_I_pr0_flk       &
       !                 * (1.-R_Tstar_icesnow)                     ! Add flux due to heat conduction
       !      d_T_ice_dt = d_T_ice_dt - R_Tstar_icesnow*Q_w_flk            ! Add flux from water to ice
       !      d_T_ice_dt = d_T_ice_dt + flk_str_2                          ! Add atmospheric fluxes
       !      d_T_ice_dt = d_T_ice_dt/tpl_rho_I/tpl_c_I                    ! Total forcing
       !      d_T_ice_dt = d_T_ice_dt/flk_str_1                            ! dT_ice/dt 
       !      T_ice_n_flk = T_ice_p_flk + d_T_ice_dt*del_time                          ! Advance T_ice
       !    END IF
       !
       !    Phi_I_pr1_flk = MIN(1., h_ice_n_flk/H_Ice_max)          ! h_ice relative to its maximum value
       !    Phi_I_pr1_flk = Phi_I_pr1_lin + Phi_I_ast_MR*Phi_I_pr1_flk     ! d\Phi_I(1)/d\zeta_I (ice)
       !    R_H_icesnow = Phi_I_pr1_flk/Phi_S_pr0_lin*tpl_kappa_I/lake_snowheatconduct(h_snow_n_flk) &
       !                 *h_snow_n_flk/MAX(h_ice_n_flk, h_Ice_min_flk)
       !    T_snow_n_flk = T_ice_n_flk + R_H_icesnow*(T_ice_n_flk-tpl_T_f)             ! Snow temperature
       !
       !  END IF No_Melting


    END IF Ice_exist


    ! Security, limit h_ice by its maximum value
    h_ice_n_flk = MIN(h_ice_n_flk, H_Ice_max)      

    ! Security, limit the ice and snow temperatures by the freezing point 
    T_snow_n_flk = MIN(T_snow_n_flk, tpl_T_f)  
    T_ice_n_flk =  MIN(T_ice_n_flk,  tpl_T_f)    

    !_tmp
    ! Security, avoid too low values (these constraints are used for debugging purposes)
    !  T_snow_n_flk = MAX(T_snow_n_flk, 73.15)  
    !  T_ice_n_flk =  MAX(T_ice_n_flk,  73.15)    
    !_tmp

    ! Remove too thin ice and/or snow
    IF(h_ice_n_flk.LT.h_Ice_min_flk)  THEN        ! Check ice
       h_ice_n_flk = 0.       ! Ice is too thin, remove it, and
       T_ice_n_flk = tpl_T_f         ! set T_ice to the freezing point.
       h_snow_n_flk = 0.      ! Remove snow when there is no ice, and
       T_snow_n_flk = tpl_T_f        ! set T_snow to the freezing point.
       l_ice_create = .FALSE.        ! "Exotic" case, ice has been created but proved to be too thin
    ELSE IF(h_snow_n_flk.LT.h_Snow_min_flk) THEN  ! Ice exists, check snow
       h_snow_n_flk = 0.      ! Snow is too thin, remove it, 
       T_snow_n_flk = T_ice_n_flk    ! and set the snow temperature equal to the ice temperature.
    END IF


    !------------------------------------------------------------------------------
    !  Compute the mean temperature of the water column.
    !------------------------------------------------------------------------------

    IF(l_ice_create) Q_w_flk = 0.     ! Ice has just been created, set Q_w to zero
    d_T_mnw_dt = (Q_w_flk - Q_bot_flk + I_w_flk - I_bot_flk)/tpl_rho_w_r/tpl_c_w/depth_w
    T_mnw_n_flk = T_mnw_p_flk + d_T_mnw_dt*del_time   ! Advance T_mnw
    T_mnw_n_flk = MAX(T_mnw_n_flk, tpl_T_f)           ! Limit T_mnw by the freezing point 


    !------------------------------------------------------------------------------
    !  Compute the mixed-layer depth, the mixed-layer temperature, 
    !  the bottom temperature and the shape factor
    !  with respect to the temperature profile in the thermocline. 
    !  Different formulations are used, depending on the regime of mixing. 
    !------------------------------------------------------------------------------

    HTC_Water: IF(h_ice_n_flk.GE.h_Ice_min_flk) THEN    ! Ice exists

       T_mnw_n_flk = MIN(T_mnw_n_flk, tpl_T_r) ! Limit the mean temperature under the ice by T_r 
       T_wML_n_flk = tpl_T_f                   ! The mixed-layer temperature is equal to the freezing point 

       IF(l_ice_create) THEN                  ! Ice has just been created 
          IF(h_ML_p_flk.GE.depth_w-h_ML_min_flk) THEN    ! h_ML=D when ice is created 
             h_ML_n_flk = 0.                 ! Set h_ML to zero 
             C_T_n_flk = C_T_min                    ! Set C_T to its minimum value 
          ELSE                                          ! h_ML<D when ice is created 
             h_ML_n_flk = h_ML_p_flk                ! h_ML remains unchanged 
             C_T_n_flk = C_T_p_flk                  ! C_T (thermocline) remains unchanged 
          END IF
          T_bot_n_flk = T_wML_n_flk - (T_wML_n_flk-T_mnw_n_flk)/C_T_n_flk/(1.-h_ML_n_flk/depth_w)
          ! Update the bottom temperature 

       ELSE IF(T_bot_p_flk.LT.tpl_T_r) THEN   ! Ice exists and T_bot < T_r, molecular heat transfer 
          h_ML_n_flk = h_ML_p_flk                  ! h_ML remains unchanged 
          C_T_n_flk = C_T_p_flk                    ! C_T (thermocline) remains unchanged 
          T_bot_n_flk = T_wML_n_flk - (T_wML_n_flk-T_mnw_n_flk)/C_T_n_flk/(1.-h_ML_n_flk/depth_w)
          ! Update the bottom temperature 

       ELSE                                   ! Ice exists and T_bot = T_r, convection due to bottom heating 
          T_bot_n_flk = tpl_T_r                      ! T_bot is equal to the temperature of maximum density 
          IF(h_ML_p_flk.GE.c_small_flk) THEN   ! h_ML > 0 
             C_T_n_flk = C_T_p_flk                     ! C_T (thermocline) remains unchanged 
             h_ML_n_flk = depth_w*(1.-(T_wML_n_flk-T_mnw_n_flk)/(T_wML_n_flk-T_bot_n_flk)/C_T_n_flk)
             h_ML_n_flk = MAX(h_ML_n_flk, 0.)   ! Update the mixed-layer depth  
          ELSE                                 ! h_ML = 0 
             h_ML_n_flk = h_ML_p_flk                   ! h_ML remains unchanged 
             C_T_n_flk = (T_wML_n_flk-T_mnw_n_flk)/(T_wML_n_flk-T_bot_n_flk) 
             C_T_n_flk = MIN(C_T_max, MAX(C_T_n_flk, C_T_min)) ! Update the shape factor (thermocline)  
          END IF
       END IF

       T_bot_n_flk = MIN(T_bot_n_flk, tpl_T_r)    ! Security, limit the bottom temperature by T_r 

    ELSE HTC_Water                                      ! Open water

       ! Generalised buoyancy flux scale and convective velocity scale
       flk_str_1 =  tpl_grav*tpl_a_T*(T_wML_p_flk-tpl_T_r)*Q_star_flk/tpl_rho_w_r/tpl_c_w                    
       IF(flk_str_1.LT.0.) THEN       
          w_star_sfc_flk = (-flk_str_1*h_ML_p_flk)**(1./3.)  ! Convection     
       ELSE 
          w_star_sfc_flk = 0.                                       ! Neutral or stable stratification
       END IF

       !_dm
       ! The equilibrium depth of the CBL due to surface cooling with the volumetric heating
       ! is not computed as a solution to the transcendental equation.
       ! Instead, an algebraic formula is used
       ! that interpolates between the two asymptotic limits.
       !_dm
       conv_equil_h_scale = -Q_w_flk/MAX(I_w_flk, c_small_flk)
       IF(conv_equil_h_scale.GT.0. .AND. conv_equil_h_scale.LT.1.  &
            .AND. T_wML_p_flk.GT.tpl_T_r) THEN   ! The equilibrium CBL depth scale is only used above T_r
          conv_equil_h_scale = SQRT(6.*conv_equil_h_scale)                 &
               + 2.*conv_equil_h_scale/(1.-conv_equil_h_scale)
          conv_equil_h_scale = MIN(depth_w, conv_equil_h_scale/extincoef_water_typ)
       ELSE
          conv_equil_h_scale = 0.       ! Set the equilibrium CBL depth to zero
       END IF

       ! Mean buoyancy frequency in the thermocline
       N_T_mean = tpl_grav*tpl_a_T*((0.5*(T_wML_p_flk+T_bot_p_flk))-tpl_T_r)*(T_wML_p_flk-T_bot_p_flk)
       IF(h_ML_p_flk.LE.depth_w-h_ML_min_flk) THEN
          N_T_mean = SQRT(N_T_mean/(depth_w-h_ML_p_flk))  ! Compute N                   
       ELSE 
          N_T_mean = 0.                            ! h_ML=D, set N to zero
       END IF

       ! The rate of change of C_T
       d_C_T_dt = MAX(w_star_sfc_flk, u_star_w_flk, u_star_min_flk)**2
       d_C_T_dt = N_T_mean*(depth_w-h_ML_p_flk)**2       &
            / c_relax_C/d_C_T_dt                               ! Relaxation time scale for C_T
       d_C_T_dt = (C_T_max-C_T_min)/MAX(d_C_T_dt, c_small_flk)     ! Rate-of-change of C_T 


       !-------------------------------------------------
       IF (d_C_T_dt .GT. 0) THEN
          d_C_T_dt = MIN(d_C_T_dt, 0.00000278 )  ! ref Lemoigne2009 Crelax max = 0.01h-1
       ELSE
          d_C_T_dt = MAX(d_C_T_dt, -0.00000278 ) 
       ENDIF
       !-------------------------------------------------


       ! Compute the shape factor and the mixed-layer depth, 
       ! using different formulations for convection and wind mixing

       C_TT_flk = C_TT_1*C_T_p_flk-C_TT_2         ! C_TT, using C_T at the previous time step
       C_Q_flk = 2.*C_TT_flk/C_T_p_flk     ! C_Q using C_T at the previous time step

       Mixing_regime: IF(flk_str_1.LT.0.) THEN  ! Convective mixing 

          C_T_n_flk = C_T_p_flk + d_C_T_dt*del_time                        ! Update C_T, assuming dh_ML/dt>0
          C_T_n_flk = MIN(C_T_max, MAX(C_T_n_flk, C_T_min))                ! Limit C_T 
          d_C_T_dt = (C_T_n_flk-C_T_p_flk)/del_time                        ! Re-compute dC_T/dt

          IF(h_ML_p_flk.LE.depth_w-h_ML_min_flk) THEN       ! Compute dh_ML/dt
             IF(h_ML_p_flk.LE.h_ML_min_flk) THEN    ! Use a reduced entrainment equation (spin-up)
                d_h_ML_dt = c_cbl_1/c_cbl_2*MAX(w_star_sfc_flk, c_small_flk)


             ELSE                                   ! Use a complete entrainment equation 
                R_H_icesnow     = depth_w/h_ML_p_flk
                R_rho_c_icesnow = R_H_icesnow-1.
                R_TI_icesnow    = C_T_p_flk/C_TT_flk
                R_Tstar_icesnow = (R_TI_icesnow/2.-1.)*R_rho_c_icesnow + 1.
                d_h_ML_dt = -Q_star_flk*(R_Tstar_icesnow*(1.+c_cbl_1)-1.) - Q_bot_flk
                d_h_ML_dt = d_h_ML_dt/tpl_rho_w_r/tpl_c_w                        ! Q_* and Q_b flux terms
                flk_str_2 = (depth_w-h_ML_p_flk)*(T_wML_p_flk-T_bot_p_flk)*C_TT_2/C_TT_flk*d_C_T_dt 
                d_h_ML_dt = d_h_ML_dt + flk_str_2                                 ! Add dC_T/dt term
                flk_str_2 = I_bot_flk + (R_TI_icesnow-1.)*I_h_flk - R_TI_icesnow*I_intm_h_D_flk
                flk_str_2 = flk_str_2 + (R_TI_icesnow-2.)*R_rho_c_icesnow*(I_h_flk-I_intm_0_h_flk)
                flk_str_2 = flk_str_2/tpl_rho_w_r/tpl_c_w
                d_h_ML_dt = d_h_ML_dt + flk_str_2                                 ! Add radiation terms
                flk_str_2 = -c_cbl_2*R_Tstar_icesnow*Q_star_flk/tpl_rho_w_r/tpl_c_w/MAX(w_star_sfc_flk, c_small_flk)
                flk_str_2 = flk_str_2 + C_T_p_flk*(T_wML_p_flk-T_bot_p_flk)
                d_h_ML_dt = d_h_ML_dt/flk_str_2                                   ! dh_ML/dt = r.h.s.
             END IF
             !_dm
             ! Notice that dh_ML/dt may appear to be negative  
             ! (e.g. due to buoyancy loss to bottom sediments and/or
             ! the effect of volumetric radiation heating),
             ! although a negative generalized buoyancy flux scale indicates 
             ! that the equilibrium CBL depth has not yet been reached
             ! and convective deepening of the mixed layer should take place.
             ! Physically, this situation reflects an approximate character of the lake model.
             ! Using the self-similar temperature profile in the thermocline, 
             ! there is always communication between the mixed layer, the thermocline 
             ! and the lake bottom. As a result, the rate of change of the CBL depth
             ! is always dependent on the bottom heat flux and the radiation heating of the thermocline.
             ! In reality, convective mixed-layer deepening may be completely decoupled
             ! from the processes underneath. In order to account for this fact,
             ! the rate of CBL deepening is set to a small value
             ! if dh_ML/dt proves to be negative.
             ! This is "double insurance" however, 
             ! as a negative dh_ML/dt is encountered very rarely.
             !_dm


             d_h_ML_dt = MAX(d_h_ML_dt, c_small_flk)    
             h_ML_n_flk = h_ML_p_flk + d_h_ML_dt*del_time                       ! Update h_ML 
             h_ML_n_flk = MAX(h_ML_min_flk, MIN(h_ML_n_flk, depth_w))           ! Security, limit h_ML
          ELSE                                              ! Mixing down to the lake bottom
             h_ML_n_flk = depth_w
          END IF

       ELSE Mixing_regime                              ! Wind mixing

          d_h_ML_dt = MAX(u_star_w_flk, u_star_min_flk)                        ! The surface friction velocity
          ZM_h_scale = (ABS(par_Coriolis)/c_sbl_ZM_n + N_T_mean/c_sbl_ZM_i)*d_h_ML_dt**2
          ZM_h_scale = ZM_h_scale + flk_str_1/c_sbl_ZM_s
          ZM_h_scale = MAX(ZM_h_scale, c_small_flk)
          ZM_h_scale = d_h_ML_dt**3/ZM_h_scale 
          ZM_h_scale = MAX(h_ML_min_flk, MIN(ZM_h_scale, h_ML_max_flk))        ! The ZM96 SBL depth scale 
          ZM_h_scale = MAX(ZM_h_scale, conv_equil_h_scale)                     ! Equilibrium mixed-layer depth 

          !_dm 
          ! In order to avoid numerical discretization problems,
          ! an analytical solution to the evolution equation 
          ! for the wind-mixed layer depth is used.
          ! That is, an exponential relaxation formula is applied
          ! over the time interval equal to the model time step.
          !_dm 

          d_h_ML_dt = c_relax_h*d_h_ML_dt/ZM_h_scale*del_time
          h_ML_n_flk = ZM_h_scale - (ZM_h_scale-h_ML_p_flk)*EXP(-d_h_ML_dt)    ! Update h_ML 
          h_ML_n_flk = MAX(h_ML_min_flk, MIN(h_ML_n_flk, depth_w))             ! Limit h_ML 
          d_h_ML_dt = (h_ML_n_flk-h_ML_p_flk)/del_time                         ! Re-compute dh_ML/dt

          IF(h_ML_n_flk.LE.h_ML_p_flk)           &
               d_C_T_dt = -d_C_T_dt                 ! Mixed-layer retreat or stationary state, dC_T/dt<0
          C_T_n_flk = C_T_p_flk + d_C_T_dt*del_time                            ! Update C_T
          C_T_n_flk = MIN(C_T_max, MAX(C_T_n_flk, C_T_min))                    ! Limit C_T 
          d_C_T_dt = (C_T_n_flk-C_T_p_flk)/del_time                            ! Re-compute dC_T/dt

       END IF Mixing_regime

       ! Compute the time-rate-of-change of the the bottom temperature, 
       ! depending on the sign of dh_ML/dt 
       ! Update the bottom temperature and the mixed-layer temperature

       IF(h_ML_n_flk.LE.depth_w-h_ML_min_flk) THEN       ! Mixing did not reach the bottom 

          IF(h_ML_n_flk.GT.h_ML_p_flk) THEN   ! Mixed-layer deepening 
             R_H_icesnow     = h_ML_p_flk/depth_w
             R_rho_c_icesnow = 1.-R_H_icesnow 
             R_TI_icesnow    = 0.5*C_T_p_flk*R_rho_c_icesnow+C_TT_flk*(2.*R_H_icesnow-1.)
             R_Tstar_icesnow = (0.5+C_TT_flk-C_Q_flk)/R_TI_icesnow
             R_TI_icesnow    = (1.-C_T_p_flk*R_rho_c_icesnow)/R_TI_icesnow

             d_T_bot_dt = (Q_w_flk-Q_bot_flk+I_w_flk-I_bot_flk)/tpl_rho_w_r/tpl_c_w
             d_T_bot_dt = d_T_bot_dt - C_T_p_flk*(T_wML_p_flk-T_bot_p_flk)*d_h_ML_dt
             d_T_bot_dt = d_T_bot_dt*R_Tstar_icesnow/depth_w                   ! Q+I fluxes and dh_ML/dt term

             flk_str_2 = I_intm_h_D_flk - (1.-C_Q_flk)*I_h_flk - C_Q_flk*I_bot_flk
             flk_str_2 = flk_str_2*R_TI_icesnow/(depth_w-h_ML_p_flk)/tpl_rho_w_r/tpl_c_w
             d_T_bot_dt = d_T_bot_dt + flk_str_2                               ! Add radiation-flux term

             flk_str_2 = (1.-C_TT_2*R_TI_icesnow)/C_T_p_flk
             flk_str_2 = flk_str_2*(T_wML_p_flk-T_bot_p_flk)*d_C_T_dt
             d_T_bot_dt = d_T_bot_dt + flk_str_2                               ! Add dC_T/dt term
          ELSE                                ! Mixed-layer retreat or stationary state
             d_T_bot_dt = 0.                                            ! dT_bot/dt=0
          END IF

          T_bot_n_flk = T_bot_p_flk + d_T_bot_dt*del_time                      ! Update T_bot  
          T_bot_n_flk = MAX(T_bot_n_flk, tpl_T_f)           ! Security, limit T_bot by the freezing point
          flk_str_2 = (T_bot_n_flk-tpl_T_r)*tpl_grav*tpl_a_T*(T_mnw_n_flk-tpl_T_r)
          IF(flk_str_2.LT.0.) T_bot_n_flk = tpl_T_r  ! Security, avoid T_r crossover 
          T_wML_n_flk = C_T_n_flk*(1.-h_ML_n_flk/depth_w)
          T_wML_n_flk = (T_mnw_n_flk-T_bot_n_flk*T_wML_n_flk)/(1.-T_wML_n_flk)
          T_wML_n_flk = MAX(T_wML_n_flk, tpl_T_f)           ! Security, limit T_wML by the freezing point

       ELSE                                              ! Mixing down to the lake bottom 

          h_ML_n_flk = depth_w
          T_wML_n_flk = T_mnw_n_flk
          T_bot_n_flk = T_mnw_n_flk
          C_T_n_flk = C_T_min

       END IF

    END IF HTC_Water


    !------------------------------------------------------------------------------
    !  Compute the depth of the upper layer of bottom sediments
    !  and the temperature at that depth.
    !------------------------------------------------------------------------------

    IF (ok_botsed_lake) THEN   ! The bottom-sediment scheme is used

       IF(H_B1_p_flk.GE.depth_sed_lake-H_B1_min_flk) THEN   ! No T(z) maximum (no thermal wave) 
          H_B1_n_flk = 0.                                   ! Set H_B1_p to zero
          T_B1_n_flk = T_bot_p_flk                          ! Set T_B1_p to the bottom temperature
       END IF

       flk_str_1 = 2.*Phi_B1_pr0/(1.-C_B1)*tpl_kappa_w/tpl_rho_w_r/tpl_c_w*del_time
       h_ice_threshold = SQRT(flk_str_1)                              ! Threshold value of H_B1
       h_ice_threshold = MIN(0.9*depth_sed_lake, h_ice_threshold)     ! Limit H_B1
       flk_str_2 = C_B2/(1.-C_B2)*(t_sed_lake-T_B1_n_flk)/(depth_sed_lake-H_B1_p_flk)

       IF(H_B1_p_flk.LT.h_ice_threshold) THEN            ! Use a truncated equation for H_B1(t)
          H_B1_n_flk = SQRT(H_B1_p_flk**2+flk_str_1)     ! Advance H_B1
          d_H_B1_dt = (H_B1_n_flk-H_B1_p_flk)/del_time   ! Re-compute dH_B1/dt 
       ELSE                                              ! Use a full equation for H_B1(t)
          flk_str_1 = (Q_bot_flk+I_bot_flk)/H_B1_p_flk/tpl_rho_w_r/tpl_c_w
          flk_str_1 = flk_str_1 - (1.-C_B1)*(T_bot_n_flk-T_bot_p_flk)/del_time
          d_H_B1_dt = (1.-C_B1)*(T_bot_p_flk-T_B1_n_flk)/H_B1_n_flk + C_B1*flk_str_2
          d_H_B1_dt = flk_str_1/d_H_B1_dt
          H_B1_n_flk = H_B1_n_flk + d_H_B1_dt*del_time         ! Advance H_B1
       END IF
       d_T_B1_dt = flk_str_2*d_H_B1_dt
       T_B1_n_flk = T_B1_p_flk + d_T_B1_dt*del_time            ! Advance T_B1



       !_nu  
       ! Use a very simplistic procedure, where only the upper layer profile is used, 
       ! H_B1 is always set to depth_sed_lake, and T_B1 is always set to t_sed_lake.
       ! Then, the time derivatives are zero, and the sign of the bottom heat flux depends on 
       ! whether T_bot is smaller or greater than t_sed_lake.
       ! This is, of course, an oversimplified scheme.
       !_nu  d_H_B1_dt = 0.
       !_nu  d_T_B1_dt = 0.
       !_nu  H_B1_n_flk = H_B1_p_flk + d_H_B1_dt*del_time   ! Advance H_B1
       !_nu  T_B1_n_flk = T_B1_p_flk + d_T_B1_dt*del_time   ! Advance T_B1
       !_nu  

       l_snow_exists = H_B1_n_flk.GE.depth_sed_lake-H_B1_min_flk        & ! H_B1 reached depth_sed_lake, or
            .OR. H_B1_n_flk.LT.H_B1_min_flk                             & ! H_B1 decreased to zero, or
            .OR.(T_bot_n_flk-T_B1_n_flk)*(t_sed_lake-T_B1_n_flk).LE.0.    ! there is no T(z) maximum
       IF(l_snow_exists) THEN      
          H_B1_n_flk = depth_sed_lake                     ! Set H_B1 to the depth of the thermally active layer
          T_B1_n_flk = t_sed_lake                         ! Set T_B1 to the climatological temperature 
       END IF

    ELSE 
       ! The bottom-sediment scheme is not used

       H_B1_n_flk = depth_sed_lake                 ! H_B1 is set to a reference value 
       T_B1_n_flk = tpl_T_r                        ! T_B1 is set to the temperature of maximum density

    END IF ! ok_botsed_lake


    !------------------------------------------------------------------------------
    !  Impose additional constraints.
    !------------------------------------------------------------------------------

    ! In case of unstable stratification, force mixing down to the bottom
    flk_str_2 = (T_wML_n_flk-T_bot_n_flk)*tpl_grav*tpl_a_T*(T_mnw_n_flk-tpl_T_r)
    IF(flk_str_2.LT.0.) THEN 
       h_ML_n_flk = depth_w
       T_wML_n_flk = T_mnw_n_flk
       T_bot_n_flk = T_mnw_n_flk
       C_T_n_flk = C_T_min

    END IF


    !------------------------------------------------------------------------------
    !  Update the surface temperature.
    !------------------------------------------------------------------------------

    IF(h_snow_n_flk.GE.h_Snow_min_flk) THEN   
       T_sfc_n = T_snow_n_flk                   ! Snow exists, use the snow temperature
    ELSE IF(h_ice_n_flk.GE.h_Ice_min_flk) THEN
       T_sfc_n = T_ice_n_flk                    ! Ice exists but there is no snow, use the ice temperature
    ELSE 
       T_sfc_n = T_wML_n_flk                    ! No ice-snow cover, use the mixed-layer temperature
    END IF

  END SUBROUTINE lake_driver




  !! ================================================================================================================================
  !! SUBROUTINE 	: lake_SfcFlx_roughness
  !!
  !>\BRIEF          Computes the water-surface or the ice-surface roughness lengths
  !!                with respect to wind velocity, potential temperature and specific humidity.
  !!
  !! DESCRIPTION :  Computes the water-surface or the ice-surface roughness lengths
  !!                with respect to wind velocity, potential temperature and specific humidity.
  !!
  !!                The water-surface roughness lengths with respect to wind velocity is computed
  !!                from the Charnock formula when the surface is aerodynamically rough.
  !!                A simple empirical formulation is used to account for the dependence 
  !!                of the Charnock parameter on the wind fetch. 
  !!                When the flow is aerodynamically smooth, the roughness length with respect to 
  !!                wind velocity is proportional to the depth of the viscous sub-layer.
  !!                The water-surface roughness lengths for scalars are computed using the power-law 
  !!                formulations in terms of the roughness Reynolds number (Zilitinkevich et al. 2001).
  !!                The ice-surface aerodynamic roughness is taken to be constant.
  !!                The ice-surface roughness lengths for scalars 
  !!                are computed through the power-law formulations 
  !!                in terms of the roughness Reynolds number (Andreas 2002).
  !!
  !!
  !! MAIN OUTPUT VARIABLE(S) :  c_z0u_fetch,  u_star_thresh, z0u, z0t, z0q 
  !!
  !! REFERENCE(S) : 
  !!
  !! FLOWCHART    : None
  !! \n
  !_ ================================================================================================================================

  SUBROUTINE lake_SfcFlx_roughness (fetch, U_a, u_star, h_ice,   & 
       c_z0u_fetch, u_star_thresh, z0u, z0t, z0q)

    !  Input (procedure arguments)
    REAL (KIND = r_std), INTENT(IN) ::   &
         fetch                             , & ! Typical wind fetch [m]
         U_a                               , & ! Wind speed [m s^{-1}]
         u_star                            , & ! Friction velocity in the surface air layer [m s^{-1}]
         h_ice                                 ! Ice thickness [m]

    !  Output (procedure arguments)
    REAL (KIND = r_std), INTENT(OUT) ::   &
         c_z0u_fetch                        , & ! Fetch-dependent Charnock parameter
         u_star_thresh                      , & ! Threshold value of friction velocity [m s^{-1}]
         z0u                                , & ! Roughness length with respect to wind velocity [m]
         z0t                                , & ! Roughness length with respect to potential temperature [m]
         z0q                                    ! Roughness length with respect to specific humidity [m]

    !  Local variables of type REAL
    REAL (KIND = r_std) ::    &
         Re_s                   , & ! Surface Reynolds number 
         Re_s_thresh                ! Threshold value of Re_s

    Water_or_Ice: IF(h_ice.LT.h_Ice_min_flk) THEN  ! Water surface  

       ! The Charnock parameter as dependent on dimensionless fetch
       c_z0u_fetch = MAX(U_a, u_wind_min_sf)**2/tpl_grav/fetch  ! Inverse dimensionless fetch
       c_z0u_fetch = c_z0u_rough + c_z0u_ftch_f*c_z0u_fetch**c_z0u_ftch_ex
       c_z0u_fetch = MIN(c_z0u_fetch, c_z0u_rough_L)                      ! Limit Charnock parameter

       ! Threshold value of friction velocity
       u_star_thresh = (c_z0u_smooth/c_z0u_fetch*tpl_grav*tpsf_nu_u_a)**num_1o3_sf

       ! Surface Reynolds number and its threshold value
       Re_s = u_star**3/tpsf_nu_u_a/tpl_grav
       Re_s_thresh = c_z0u_smooth/c_z0u_fetch

       ! Aerodynamic roughness
       IF(Re_s.LE.Re_s_thresh) THEN                 
          z0u = c_z0u_smooth*tpsf_nu_u_a/u_star     ! Smooth flow
       ELSE
          z0u = c_z0u_fetch*u_star*u_star/tpl_grav  ! Rough flow
       END IF
       ! Roughness for scalars  
       z0q = c_z0u_fetch*MAX(Re_s, Re_s_thresh)
       z0t = c_z0t_rough_1*z0q**c_z0t_rough_3 - c_z0t_rough_2
       z0q = c_z0q_rough_1*z0q**c_z0q_rough_3 - c_z0q_rough_2
       z0t = z0u*EXP(-c_Karman/Pr_neutral*z0t)
       z0q = z0u*EXP(-c_Karman/Sc_neutral*z0q) 

    ELSE Water_or_Ice                              ! Ice surface

       ! The Charnock parameter is not used over ice, formally set "c_z0u_fetch" to its minimum value
       c_z0u_fetch = c_z0u_rough

       ! Threshold value of friction velocity
       u_star_thresh = c_z0u_smooth*tpsf_nu_u_a/z0u_ice_rough

       ! Aerodynamic roughness
       z0u = MAX(z0u_ice_rough, c_z0u_smooth*tpsf_nu_u_a/u_star)

       ! Roughness Reynolds number 
       Re_s = MAX(u_star*z0u/tpsf_nu_u_a, c_accur_sf)

       ! Roughness for scalars  
       IF(Re_s.LE.Re_z0s_ice_t) THEN 
          z0t = c_z0t_ice_b0t + c_z0t_ice_b1t*LOG(Re_s)
          z0t = MIN(z0t, c_z0t_ice_b0s)
          z0q = c_z0q_ice_b0t + c_z0q_ice_b1t*LOG(Re_s)
          z0q = MIN(z0q, c_z0q_ice_b0s)
       ELSE 
          z0t = c_z0t_ice_b0r + c_z0t_ice_b1r*LOG(Re_s) + c_z0t_ice_b2r*LOG(Re_s)**2
          z0q = c_z0q_ice_b0r + c_z0q_ice_b1r*LOG(Re_s) + c_z0q_ice_b2r*LOG(Re_s)**2
       END IF
       z0t = z0u*EXP(z0t)
       z0q = z0u*EXP(z0q)

    END IF Water_or_Ice

  END SUBROUTINE lake_SfcFlx_roughness



  !! ================================================================================================================================
  !! FUNCTION 	: lake_snowdensity
  !!
  !>\BRIEF          Computes the snow density,
  !!                using an empirical approximation from Heise et al. (2003).
  !!
  !! DESCRIPTION :  Computes the snow density,
  !!                using an empirical approximation from Heise et al. (2003).
  !!
  !! MAIN OUTPUT VARIABLE(S) : lake_snowdensity
  !!
  !! REFERENCE(S) : Heise et al. (2003)
  !!
  !! FLOWCHART    : None
  !! \n
  !_ ================================================================================================================================


  REAL (r_std) FUNCTION lake_snowdensity(h_snow)

    IMPLICIT NONE

    !! 0. Variable and parameter declaration
    !  Input (function argument) 
    REAL (r_std), INTENT(IN) :: h_snow          !! Snow thickness [m]

    ! Snow density [kg/m3]

    !  Security. Ensure that the expression in () does not become negative at a very large h_snow.
    !  lake_snowdensity = MAX( c_small_flk, (1. - h_snow*tpl_Gamma_rho_S/tpl_rho_w_r) )
    !  lake_snowdensity = MIN( tpl_rho_S_max, tpl_rho_S_min/lake_snowdensity )

    lake_snowdensity = 320.

  END FUNCTION lake_snowdensity



  !! ================================================================================================================================
  !! FUNCTION 	: lake_snowheatconduct
  !!
  !>\BRIEF         Computes the snow heat conductivity,
  !!               using an empirical approximation from Heise et al. (2003).
  !!
  !! DESCRIPTION :  Computes the snow heat conductivity,
  !!                using an empirical approximation from Heise et al. (2003).
  !!
  !! MAIN OUTPUT VARIABLE(S) : lake_snowheatconduct
  !!
  !! REFERENCE(S) : Heise et al. (2003)
  !!
  !! FLOWCHART    : None
  !! \n
  !_ ================================================================================================================================

  REAL (KIND = r_std) FUNCTION lake_snowheatconduct(h_snow)

    IMPLICIT NONE

    !! 0. Variable and parameter declaration

    !  Input (function argument) 
    REAL(r_std), INTENT(IN) ::  h_snow      !! Snow thickness [m]

    ! Snow heat conductivity [J m^{-1} s^{-1} K^{-1} = kg m s^{-3} K^{-1}]

    !  lake_snowheatconduct = lake_snowdensity( h_snow )   ! Compute snow density
    !  lake_snowheatconduct = MIN( tpl_kappa_S_max, tpl_kappa_S_min                      &
    !                        + h_snow*tpl_Gamma_kappa_S*lake_snowheatconduct/tpl_rho_w_r )

    lake_snowheatconduct = 0.14+(tpl_kappa_I-0.14)*exp(-5.*h_snow) 

  END FUNCTION lake_snowheatconduct

END MODULE lake
