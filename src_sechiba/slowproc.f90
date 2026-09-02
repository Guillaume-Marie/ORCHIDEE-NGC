! =================================================================================================================================
! MODULE       : slowproc
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF         Groups the subroutines that: (1) initialize all variables used in 
!! slowproc_main, (2) prepare the restart file for the next simulation, (3) Update the 
!! vegetation cover if needed, and (4) handle all slow processes if the carbon
!! cycle is activated (call STOMATE) or update the vegetation properties (LAI and 
!! fractional cover) in the case of a run with only SECHIBA.
!!
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): Allowed reading of USDA map, Nov 2014, ADucharne
!!                   November 2020: It is possible to define soil hydraulic parameters from maps,
!!                   as needed for the SP-MIP project (Tafasca Salma and Ducharne Agnes).
!!                   Changes in slowproc_xios_initialize, slowproc_soilt and slowproc_finalize
!!                   Octobre 2023: New irrigation scheme. Here interpolation of new maps for
!!                   the irrigation scheme
!!
!! REFERENCE(S)	:
!!- Tafasca S. (2020). Evaluation de l’impact des propriétés du sol sur l’hydrologie simulee dans le
!! modèle ORCHIDEE, PhD thesis, Sorbonne Universite. \n
!! SVN          :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_sechiba/slowproc.f90 $
!! $Date: 2026-05-13 12:09:00 +0200 (mer. 13 mai 2026) $
!! $Revision: 9537 $
!! \n
!_ ================================================================================================================================

MODULE slowproc

  USE defprec
  USE constantes 
  USE constantes_soil
  USE pft_parameters
  USE dynamic_parameters
  USE structures
  USE ioipsl
  USE xios_orchidee
  USE ioipsl_para
  USE sechiba_io_p
  USE interpol_help
  USE stomate
  USE stomate_data
  USE sapiens_lcchange, ONLY: check_loss_gain, check_veget, &
                              adjust_delta_veget_max, check_read_vegetmax
  USE grid
  USE solar, ONLY: solarang_noon
  USE time, ONLY : dt_sechiba, dt_stomate, one_day, FirstTsYear, LastTsDay, FirstTsMonth, FirstTsDay
  USE time, ONLY : year_start, month_start, day_start, sec_start
  USE time, ONLY : month_end, day_end, julian_diff
  USE mod_orchidee_para
  USE stomate_laieff,    ONLY: effective_lai,find_lai_per_level, &
                               calculate_z_level_photo, fitting_laieff, &
                               stomate_laieff_initialize, partial_spheroid_vol
  USE function_library,  ONLY: wood_to_qmheight, cc_to_lai, &
                               check_pixel_area, check_area_change, &
                               wood_to_height, wood_to_dia, &
                               wood_to_cn_dia, wood_to_cv, calculate_c0_alloc, &
                               get_printlev, cc_to_biomass, biomass_to_lai, &
                               lai_to_biomass

 IMPLICIT NONE

  ! Private & public routines

  PRIVATE
  PUBLIC slowproc_main, slowproc_clear, slowproc_initialize, slowproc_finalize, &
       slowproc_xios_initialize, slowproc_veget, slowproc_canopy

  !
  ! variables used inside slowproc module : declaration and initialisation
  !
  REAL(r_std), SAVE                                  :: slope_default = 0.1
!$OMP THREADPRIVATE(slope_default)
  INTEGER(i_std) , SAVE                              :: Ninput_update       !! update frequency in years for N inputs (nb of years)
!$OMP THREADPRIVATE(Ninput_update)
  REAL(r_std), PARAMETER                             :: aei_sw_default = 62.2 !! Default percentage for Area Equipped for Irrigation with Surface Water (unitless). Value of 62.2% is the mean of the GMIA v5 dataset of the FAO AQUASTAT (Siebert et al., 2010) used for the construction of the AEI_SW_pct.nc
  INTEGER, SAVE                                      :: printlev_loc        !! Local printlev in slowproc module
!$OMP THREADPRIVATE(printlev_loc)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: clayfraction        !! Clayfraction (0-1, unitless)
!$OMP THREADPRIVATE(clayfraction)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: sandfraction        !! Sandfraction (0-1, unitless)
!$OMP THREADPRIVATE(sandfraction)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: siltfraction        !! Siltfraction (0-1, unitless)
!$OMP THREADPRIVATE(siltfraction)  
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: bulk                !! Bulk density (kg/m**3)
!$OMP THREADPRIVATE(bulk)  
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: soil_ph             !! Soil pH (-)
!$OMP THREADPRIVATE(soil_ph)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:,:):: n_input_raw         !! nitrogen inputs (gN/m2/day) per points, per PFT and per type of N (Nox,NHx,Fert,Manure,BNF) - Monthly values (array of 12 elements)
!$OMP THREADPRIVATE(n_input_raw)  
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:,:,:,:) :: cc_biomass_m   !! Monthly cc_biomass to prescribe the canopy structure if not calculated by STOMATE
!$OMP THREADPRIVATE(cc_biomass_m)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:,:):: cc_n_m              !! Monthly cc_n to prescribe the canopy structure if not calculated by STOMATE
!$OMP THREADPRIVATE(cc_n_m)
  INTEGER(i_std) , SAVE                              :: ninput_year         !! year for N inputs data
!$OMP THREADPRIVATE(ninput_year)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: cn_leaf_min_2D      !! Minimal leaf CN ratio
!$OMP THREADPRIVATE(cn_leaf_min_2D)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: cn_leaf_max_2D      !! Maximal leaf CN ratio
!$OMP THREADPRIVATE(cn_leaf_max_2D)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: cn_leaf_init_2D     !! Initial leaf CN ratio
!$OMP THREADPRIVATE(cn_leaf_init_2D)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: irrigated_new       !! New year area equipped for irrigation  (m^{-2})
!$OMP THREADPRIVATE(irrigated_new)

CONTAINS


!!  =============================================================================================================================
!! SUBROUTINE:    slowproc_xios_initialize
!!
!>\BRIEF	  Initialize xios dependant defintion before closing context defintion
!!
!! DESCRIPTION:	  Initialize xios dependant defintion before closing context defintion
!!
!! RECENT CHANGE(S): Initialization of XIOS to read soil hydraulic parameters from maps,
!!                   as needed for the SP-MIP project (Tafasca Salma and Ducharne Agnes).
!! \n
!_ ==============================================================================================================================

  SUBROUTINE slowproc_xios_initialize

    CHARACTER(LEN=255) :: filename, name
    LOGICAL :: lerr
    REAL(r_std) :: slope_noreinf
    LOGICAL :: get_slope
    LOGICAL :: flag

    IF (printlev>=3) WRITE(numout,*) 'In slowproc_xios_initialize'


    !!
    !! 1. Prepare for reading of soils_param file
    !!

    ! Get the file name from run.def file and set file attributes accordingly
    filename = 'soils_param.nc'
    CALL getin_p('SOILCLASS_FILE',filename)
    name = filename(1:LEN_TRIM(FILENAME)-3)
    CALL xios_orchidee_set_file_attr("soiltext_file",name=name)

    ! Determine if soiltext_file will be read. If not, deactivate the file.    
    IF (xios_interpolation .AND. restname_in=='NONE' .AND. .NOT. impsoilt) THEN
       ! Reading will be done with XIOS later
       IF (printlev>=2) WRITE(numout,*) 'Reading of soiltext variable will be done later using XIOS. The filename is ', filename
    ELSE
       ! No reading of variable soiltext, deactivate correspoding fieldgroup and file set in context_input_orchidee.xml file
       IF (printlev>=2) WRITE(numout,*) 'Reading of soiltext variable will not be done with XIOS.'
       CALL xios_orchidee_set_fieldgroup_attr("soil_text",enabled=.FALSE.)
       CALL xios_orchidee_set_file_attr("soiltext_file",enabled=.FALSE.)
    END IF
  

    !!
    !! 2. Prepare for reading of bulk variable
    !!

    ! Get the file name from run.def file and set file attributes accordingly
    filename = 'soil_bulk_and_ph.nc' 
    CALL getin_p('SOIL_BULK_FILE',filename)
    
    name = filename(1:LEN_TRIM(FILENAME)-3)
    CALL xios_orchidee_set_file_attr("soilbulk_file",name=name)

    ! Set variables that can be used in the xml files
    lerr=xios_orchidee_setvar('bulk_default',bulk_default)
    
    ! Determine if the file will be read by XIOS. If not, deactivate reading of the file.    
    IF (xios_interpolation .AND. restname_in=='NONE' .AND. .NOT. impsoilt) THEN
       ! Reading will be done with XIOS later
       IF (printlev>=2) WRITE(numout,*) 'Reading of soilbulk file will be done later using XIOS. The filename is ', filename
    ELSE
       ! No reading by XIOS, deactivate soilbulk file and related variables declared in context_input_orchidee.xml.
       ! If this is not done, the model will crash if the file is not available in the run directory.
       IF (printlev>=2) WRITE(numout,*) 'Reading of soil_bulk file will not be done with XIOS.'
       CALL xios_orchidee_set_file_attr("soilbulk_file",enabled=.FALSE.)
       CALL xios_orchidee_set_field_attr("soilbulk",enabled=.FALSE.)
       CALL xios_orchidee_set_field_attr("soilbulk_mask",enabled=.FALSE.)
    END IF
    
    !!
    !! 3. Prepare for reading of soil ph variable
    !!

    ! Get the file name from run.def file and set file attributes accordingly
    ! soilbulk and soilph are by default in the same file but they can also be read from different files. 
    filename = 'soil_bulk_and_ph.nc' 
    CALL getin_p('SOIL_PH_FILE',filename)
    
    name = filename(1:LEN_TRIM(FILENAME)-3)
    CALL xios_orchidee_set_file_attr("soilph_file",name=name)

    ! Set variables that can be used in the xml files
    lerr=xios_orchidee_setvar('ph_default',ph_default)

    ! Determine if the file will be read by XIOS. If not, deactivate the file.    
    IF (xios_interpolation .AND. restname_in=='NONE' .AND. .NOT. impsoilt) THEN
       ! Reading will be done with XIOS later
       IF (printlev>=2) WRITE(numout,*) 'Reading of soilph file will be done later using XIOS. The filename is ', filename
    ELSE
       ! No reading by XIOS, deactivate soilph file and related variables declared in context_input_orchidee.xml.
       ! If this is not done, the model will crash if the file is not available in the run directory.
       IF (printlev>=2) WRITE(numout,*) 'Reading of soilph file will not be done with XIOS.'
       CALL xios_orchidee_set_file_attr("soilph_file",enabled=.FALSE.)
       CALL xios_orchidee_set_field_attr("soilph",enabled=.FALSE.)
       CALL xios_orchidee_set_field_attr("soilph_mask",enabled=.FALSE.)
    END IF
 
  
    !!
    !! 4. Prepare for reading of PFTmap file
    !!

    filename = 'PFTmap.nc'
    CALL getin_p('VEGETATION_FILE',filename)
    name = filename(1:LEN_TRIM(FILENAME)-3)
    CALL xios_orchidee_set_file_attr("PFTmap_file",name=name)

    ! Check if PFTmap file will be read by XIOS in this execution
    IF ( xios_interpolation .AND. .NOT. impveg .AND. &
         ((veget_update>0 .OR. restname_in=='NONE') .OR. vegetmap_reset)) THEN
       ! PFTmap will not be read if impveg=TRUE
       ! PFTmap file will be read each year if veget_update>0 
       ! PFTmap is read if the restart file do not exist and if impveg=F

       ! Reading will be done
       IF (printlev>=1) WRITE(numout,*) 'Reading of PFTmap file will be done later using XIOS. The filename is ', filename
    ELSE
       ! No reading, deactivate PFTmap file
       IF (printlev>=1) WRITE(numout,*) 'Reading of PFTmap file will not be done with XIOS.'
       
       CALL xios_orchidee_set_file_attr("PFTmap_file",enabled=.FALSE.)
       CALL xios_orchidee_set_field_attr("frac_veget",enabled=.FALSE.)
       CALL xios_orchidee_set_field_attr("frac_veget_frac",enabled=.FALSE.)
    ENDIF
    

    !!
    !! 5. Prepare for reading of topography file
    !!

    filename = 'cartepente2d_15min.nc'
    CALL getin_p('TOPOGRAPHY_SLOPE_FILE',filename)
    name = filename(1:LEN_TRIM(FILENAME)-3)
    CALL xios_orchidee_set_file_attr("topography_slope_file",name=name)
    
    ! Set default values used by XIOS for the interpolation
    slope_noreinf = 0.5 ! slope in percent
    CALL getin_p('SLOPE_NOREINF',slope_noreinf)
    lerr=xios_orchidee_setvar('slope_noreinf',slope_noreinf)
    lerr=xios_orchidee_setvar('slope_default',slope_default)
    
    get_slope = .FALSE.
    CALL getin_p('GET_SLOPE',get_slope)
    IF (xios_interpolation .AND. (restname_in=='NONE' .OR. get_slope)) THEN
       ! The slope file will be read using XIOS
       IF (printlev>=2) WRITE(numout,*) 'Reading of albedo file will be done later using XIOS. The filename is ', filename
    ELSE
       ! Deactivate slope reading
       IF (printlev>=2) WRITE(numout,*) 'The slope file will not be read by XIOS'
       CALL xios_orchidee_set_file_attr("topography_slope_file",enabled=.FALSE.)
       CALL xios_orchidee_set_field_attr("frac_slope_interp",enabled=.FALSE.)
       CALL xios_orchidee_set_field_attr("reinf_slope_interp",enabled=.FALSE.)
    END IF
     

    !!
    !! 6. Prepare for reading of lai file
    !!

    filename = 'lai2D.nc'
    CALL getin_p('LAI_FILE',filename)
    name = filename(1:LEN_TRIM(FILENAME)-3)
    CALL xios_orchidee_set_file_attr("lai_file",name=name)
    ! Determine if lai file will be read by XIOS. If not, deactivate the file.    
    IF (xios_interpolation .AND. restname_in=='NONE' .AND. (read_lai .EQ. 2)) THEN
       CALL ipslerr_p(3,'slowproc_xios_initialize', 'Option read_lai=2 is not yet implemented.','','')
       ! Reading will be done
       IF (printlev>=2) WRITE(numout,*) 'Reading of lai file will be done later using XIOS. The filename is ', filename
    ELSE
       ! No reading, deactivate lai file
       IF (printlev>=2) WRITE(numout,*) 'Reading of lai file will not be done with XIOS.'
       CALL xios_orchidee_set_file_attr("lai_file",enabled=.FALSE.)
       CALL xios_orchidee_set_field_attr("frac_lai_interp",enabled=.FALSE.)
       CALL xios_orchidee_set_field_attr("lai_interp",enabled=.FALSE.)
    END IF
    

    !!
    !! 7. Prepare for reading of woodharvest file
    !!

    filename = 'woodharvest.nc'
    CALL getin_p('WOODHARVEST_FILE',filename)
    name = filename(1:LEN_TRIM(FILENAME)-3)
    CALL xios_orchidee_set_file_attr("woodharvest_file",name=name)
    
    IF (xios_interpolation .AND. do_wood_harvest .AND. &
         (veget_update>0 .OR. restname_in=='NONE' )) THEN
       ! Woodharvest file will be read each year if veget_update>0 or if no restart file exists

       ! Reading will be done
       IF (printlev>=2) WRITE(numout,*) 'Reading of woodharvest file will be done later using XIOS. The filename is ', filename
    ELSE
       ! No reading, deactivate woodharvest file
       IF (printlev>=2) WRITE(numout,*) 'Reading of woodharvest file will not be done with XIOS.'
       CALL xios_orchidee_set_file_attr("woodharvest_file",enabled=.FALSE.)
       CALL xios_orchidee_set_field_attr("woodharvest_interp",enabled=.FALSE.)
    ENDIF

    !!
    !! 8. Prepare for reading of nitrogen maps
    !!
    flag=(ok_ncycle .AND. (.NOT. impose_CN .AND. .NOT. impose_ninput_dep))
    CALL slowproc_xios_initialize_ninput('Nammonium',flag)
    CALL slowproc_xios_initialize_ninput('Nnitrate',flag)
    CALL slowproc_xios_initialize_ninput('WETNHX',flag)
    CALL slowproc_xios_initialize_ninput('DRYNHX',flag)
    CALL slowproc_xios_initialize_ninput('WETNOY',flag)
    CALL slowproc_xios_initialize_ninput('DRYNOY',flag)

    flag=(ok_ncycle .AND. (.NOT. impose_CN .AND. .NOT. impose_ninput_fert))
    CALL slowproc_xios_initialize_ninput('Nfert',flag)
    CALL slowproc_xios_initialize_ninput('Nfert_cropland',flag)
    CALL slowproc_xios_initialize_ninput('Nfert_pasture',flag)
    CALL slowproc_xios_initialize_ninput('Nfert_ammo_cropland',flag)
    CALL slowproc_xios_initialize_ninput('Nfert_nitr_cropland',flag)
    CALL slowproc_xios_initialize_ninput('Nfert_cropC3',flag)
    CALL slowproc_xios_initialize_ninput('Nfert_cropC4',flag)
    CALL slowproc_xios_initialize_ninput('Nfert_ammo_pasture',flag)
    CALL slowproc_xios_initialize_ninput('Nfert_nitr_pasture',flag)


    flag=(ok_ncycle .AND. (.NOT. impose_CN .AND. .NOT. impose_ninput_manure))
    CALL slowproc_xios_initialize_ninput('Nmanure_cropland',flag)
    CALL slowproc_xios_initialize_ninput('Nmanure_pasture',flag)

    flag=(ok_ncycle .AND. (.NOT. impose_CN .AND. .NOT. impose_ninput_bnf))
    CALL slowproc_xios_initialize_ninput('Nbnf',flag)
    

    !! 9. Prepare for reading of irrigmap file
    filename = 'IRRIGmap.nc'
    CALL getin_p('IRRIGmap_FILE',filename)
    name = filename(1:LEN_TRIM(FILENAME)-3)
    CALL xios_orchidee_set_file_attr("irrigmap_file",name=name)
    
    IF (xios_interpolation .AND. do_irrigation .AND. &
        irrig_map_dynamic_flag) THEN
       ! IRRIGmap file will be read each year
       ! Reading will be done
       IF (printlev>=2) WRITE(numout,*) 'Reading of irrigmap file will be done later using XIOS. The filename is ', filename
    ELSEIF (xios_interpolation .AND. do_irrigation .AND. &
            restname_in=='NONE') THEN
       ! IRRIGmap file will be read only once if no restart file exists

       ! Reading will be done
       IF (printlev>=2) WRITE(numout,*) 'Reading of irrigmap file will be done later using XIOS. The filename is ', filename
    ELSE
       ! No reading, deactivate irrigmap file
       IF (printlev>=2) WRITE(numout,*) 'Reading of irrigmap file will not be done with XIOS.'
       CALL xios_orchidee_set_file_attr("irrigmap_file",enabled=.FALSE.)
       CALL xios_orchidee_set_field_attr("irrigmap_interp",enabled=.FALSE.)
    ENDIF
    
    !! 10. Prepare for reading of aei_sw file
    filename = 'AEI_SW_pct.nc'
    CALL getin_p('AEI_SW_FILE',filename)
    name = filename(1:LEN_TRIM(FILENAME)-3)
    CALL xios_orchidee_set_file_attr("aei_sw_file",name=name)
    
    ! Set default values used by XIOS for the interpolation
    lerr=xios_orchidee_setvar('aei_sw_default',aei_sw_default)
    
    IF (xios_interpolation .AND. do_irrigation .AND. &
        select_source_irrig .AND. restname_in=='NONE') THEN
       ! AEI_SW file will be read if no restart file exists

       ! Reading will be done
       IF (printlev>=2) WRITE(numout,*) 'Reading of aei_sw file will be done later using XIOS. The filename is ', filename
    ELSE
       ! No reading, deactivate aei_sw file
       IF (printlev>=2) WRITE(numout,*) 'Reading of aei_sw file will not be done with XIOS.'
       CALL xios_orchidee_set_file_attr("aei_sw_file",enabled=.FALSE.)
       CALL xios_orchidee_set_field_attr("aei_sw_interp",enabled=.FALSE.)
    ENDIF

    !! 11. Prepare for reading of soil parameter files

    !! See commented part below for the reading of params_sp_mip.nc if spmipexp='maps' 
    !! (with a bug, but helpful)
    !! This part was introduced to prepare the reading of params_sp_mip.nc if spmipexp='maps'
    !! but there are mistakes in the IF ELSE ENDIF and we go through ELSE
    !! each time xios_interpolation = T, even if we don't need to read this file
    !! and it is not provided by sechiba.card
    !! The corresponding part in context_input_orchidee.xml is also commented

    ! Get the file name from run.def file and set file attributes accordingly

!!$    filename = 'params_sp_mip.nc'
!!$    CALL getin_p('PARAM_FILE',filename)
!!$    name = filename(1:LEN_TRIM(FILENAME)-3)
!!$    IF (TRIM(filename)=='NONE') THEN
!!$      IF (printlev>=2) WRITE(numout,*) 'We won t readsoil hydraulic parameters from a file.'
!!$    ELSE
!!$      CALL xios_orchidee_set_file_attr("soilparam_file",name=name)
!!$      ! Determine if the file will be read by XIOS. If not, deactivate reading of the file.
!!$      IF (xios_interpolation .AND. restname_in=='NONE' .AND. .NOT. impsoilt) THEN
!!$         ! Reading will be done with XIOS later
!!$         IF (printlev>=2) WRITE(numout,*) 'Reading of soil hydraulic parameters file will be done later using XIOS. The filename is ', filename
!!$      ELSE
!!$         ! No reading by XIOS, deactivate soilparam_file and related variables declared in context_input_orchidee.xml.
!!$         ! If this is not done, the model will crash if the file is not available in the run directory.
!!$         IF (printlev>=2) WRITE(numout,*) 'Reading of soil parameter file will not be done with XIOS.'
!!$         CALL xios_orchidee_set_file_attr("soilparam_file",enabled=.FALSE.)
!!$         CALL xios_orchidee_set_field_attr("soilks",enabled=.FALSE.)
!!$         CALL xios_orchidee_set_field_attr("soilnvan",enabled=.FALSE.)
!!$         CALL xios_orchidee_set_field_attr("soilavan",enabled=.FALSE.)
!!$         CALL xios_orchidee_set_field_attr("soilmcr",enabled=.FALSE.)
!!$         CALL xios_orchidee_set_field_attr("soilmcs",enabled=.FALSE.)
!!$         CALL xios_orchidee_set_field_attr("soilmcfc",enabled=.FALSE.)
!!$         CALL xios_orchidee_set_field_attr("soilmcw",enabled=.FALSE.)
!!$      ENDIF
!!$   ENDIF

    IF (printlev_loc>=3) WRITE(numout,*) 'End slowproc_xios_intialize'
   
  END SUBROUTINE slowproc_xios_initialize


!! ================================================================================================================================
!! SUBROUTINE 	: slowproc_initialize
!!
!>\BRIEF         Initialize slowproc module and call initialization of stomate module
!!
!! DESCRIPTION : Allocate module variables, read from restart file or initialize with default values
!!               Call initialization of stomate module.
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_initialize (kjit,          kjpij,        kjpindex,                          &
                                  rest_id,       rest_id_stom, hist_id_stom,   hist_id_stom_IPCC, &
                                  IndexLand,     indexveg,     lalo,           neighbours,        &
                                  resolution,    contfrac,     temp_air,       coszang_noon,      &
                                  soiltile,      reinf_slope,  reinf_slope_soil,deadleaf_cover, assim_param,       &
                                  ks,             nvan,        avan,           mcr,               &
                                  mcs,            mcfc,         mcw,                              &
                                  frac_age,      height,       lai,            veget,             &
                                  frac_nobio,    njsc,         veget_max,      fraclut,           &
                                  nwdfraclut,    tot_bare_soil,totfrac_nobio,  qsintmax,          &
                                  temp_growth,   circ_class_biomass,                              &
                                  circ_class_n,  lai_per_level,laieff_fit,                        &
                                  z_array_out,   max_height_store,                                &
                                  som_total,     heat_Zimov,   altmax,         depth_organic_soil,&
                                  loss_gain,     veget_max_new, frac_nobio_new, height_dom,       &
                                  height_inv,    dia_dom,       dia_inv,       ba_inv,            &
                                  ind_inv,       Pgap_cumul,    irrigated_next,irrig_frac_next,   &
                                  fraction_aeirrig_sw,          precip_longterm)


!! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                             :: kjit           !! Time step number
    INTEGER(i_std), INTENT(in)                             :: kjpij          !! Total size of the un-compressed grid
    INTEGER(i_std),INTENT(in)                              :: kjpindex       !! Domain size - terrestrial pixels only
    INTEGER(i_std),INTENT (in)                             :: rest_id        !! Restart file identifier
    INTEGER(i_std),INTENT (in)                             :: rest_id_stom   !! STOMATE's _Restart_ file identifier
    INTEGER(i_std),INTENT (in)                             :: hist_id_stom   !! STOMATE's _history_ file identifier
    INTEGER(i_std),INTENT(in)                              :: hist_id_stom_IPCC !! STOMATE's IPCC _history_ file identifier
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)       :: IndexLand      !! Indices of the points on the land map
    INTEGER(i_std),DIMENSION (kjpindex*nvm), INTENT (in)   :: indexveg       !! Indices of the points on the vegetation (3D map ???) 
    REAL(r_std),DIMENSION (kjpindex,2), INTENT (in)        :: lalo           !! Geogr. coordinates (latitude,longitude) (degrees)
    INTEGER(i_std), DIMENSION (kjpindex,NbNeighb), INTENT(in):: neighbours   !! neighbouring grid points if land.
    REAL(r_std), DIMENSION (kjpindex,2), INTENT(in)        :: resolution     !! size in x an y of the grid (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)          :: contfrac       !! Fraction of continent in the grid (0-1, unitless)
    REAL(r_std), DIMENSION(kjpindex), INTENT(in)           :: temp_air       !! Air temperature at first atmospheric model layer (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)           :: coszang_noon   !! Solar zenith angle at noon
    
!! 0.2 Output variables 
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: temp_growth    !! Growth temperature (Â°C) - Is equal to t2m_month 
    INTEGER(i_std), DIMENSION(kjpindex), INTENT(out)       :: njsc           !! Index of the dominant soil textural class in the grid 
                                                                             !! cell (1-nscm, unitless)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(out)      :: height         !! height of vegetation (m)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(out)      :: height_dom     !! dominant height of vegetation (m)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(out)      :: dia_dom        !! dominant diameter of vegetation (m)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(out)      :: height_inv     !! height of vegetation as measured in a forest inventory (m)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(out)      :: dia_inv        !! diameter of vegetation as measured in a forest inventory (m)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(out)      :: ba_inv         !! ba of vegetation as measured in a forest inventory (m2)    
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(out)      :: ind_inv        !! stand density as measured in a forest inventory (trees m-1)
    REAL(r_std),DIMENSION(kjpindex,nvm), INTENT(out)       :: lai            !! PFT leaf area index (m^{2} m^{-2})
    REAL(r_std),DIMENSION (kjpindex,nvm,nleafages), INTENT(out):: frac_age   !! Age efficacity from STOMATE for isoprene
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)     :: veget          !! Fraction of vegetation type in the mesh (unitless)
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (out)  :: frac_nobio     !! Fraction of ice, lakes, cities etc. in the mesh (unitless)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)     :: veget_max      !! Maximum fraction of vegetation type in the mesh (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: tot_bare_soil  !! Total evaporating bare soil fraction in the mesh  (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: totfrac_nobio  !! Total fraction of ice+lakes+cities etc. in the mesh  (unitless)
    REAL(r_std),DIMENSION (kjpindex,nstm), INTENT(out)     :: soiltile       !! Fraction of each soil tile within vegtot (0-1, unitless)
    REAL(r_std),DIMENSION (kjpindex,nlut), INTENT(out)     :: fraclut        !! Fraction of each landuse tile (0-1, unitless)
    REAL(r_std),DIMENSION (kjpindex,nlut), INTENT(out)     :: nwdFraclut     !! Fraction of non-woody vegetation in each landuse tile (0-1, unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)          :: reinf_slope    !! slope coef for reinfiltration
    REAL(r_std),DIMENSION (kjpindex, nstm), INTENT(out)    :: reinf_slope_soil  !! slope coef for reinfiltration per soil tile
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: ks             !! Hydraulic conductivity at saturation (mm {-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: nvan           !! Van Genuchten coeficients n (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: avan           !! Van Genuchten coeficients a (mm-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: mcr            !! Residual volumetric water content (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: mcs            !! Saturated volumetric water content (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: mcfc           !! Volumetric water content at field capacity (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: mcw            !! Volumetric water content at wilting point (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex,nvm,npco2),INTENT (out):: assim_param    !! min+max+opt temperatures & vmax for photosynthesis 
                                                                             !! (K, \mumol m^{-2} s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: deadleaf_cover !! Fraction of soil covered by dead leaves (unitless)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)     :: qsintmax       !! Maximum water storage on vegetation from interception (mm) 
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm),INTENT (out):: heat_Zimov     !! heating associated with decomposition [W/m**3 soil]
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)     :: altmax         !! Maximul active layer thickness (m). Be careful, here active means non frozen.
                                                                             !! Not related with the active soil carbon pool.
    REAL(r_std), DIMENSION(kjpindex), INTENT(out)          :: depth_organic_soil !! Depth at which there is still organic matter (m)
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(out)         :: circ_class_biomass !! Biomass components of the model tree  
                                                                             !! within a circumference class
                                                                             !! class @tex $(g C ind^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(out)             :: circ_class_n   !! Number of trees within each circumference
                                                                             !! class @tex $(m^{-2})$ @endtex
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(out)      :: loss_gain      !! Changes in veget_max distributed over all
                                                                             !! age classes and thus taking the age-classes into 
                                                                             !! account (unitless, 0-1)
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT(out)   :: frac_nobio_new !! Fraction of ice,lakes,cities, ... (unitless)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(out)      :: veget_max_new  !! Maximum fraction of vegetation type including none
    REAL(r_std), DIMENSION (kjpindex), INTENT (out)        :: irrigated_next !! Dynamic irrig. area, calculated in slowproc and passed to routing
    REAL(r_std), DIMENSION (kjpindex), INTENT (out)        :: irrig_frac_next!! Dynamic irrig. fraction, calculated in slowproc and passed to routing
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)          :: fraction_aeirrig_sw  !! Fraction of area equipped for irrigation from surface water, of irrig_frac
                                                                                   !! 1.0 here corresponds to fraction of irrigated area, not grid cell
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)          :: precip_longterm !! Longterm annual precipitation sum mm year^{-1}


!! 0.3 Modified variables
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm,nelements), INTENT (inout) :: som_total !! total soil carbon for use in thermal (g/m**3)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)           :: lai_per_level  !! This is the LAI per vertical level
                                                                             !! @tex $(m^{2} m^{-2})$
    TYPE(laieff_type),DIMENSION (:,:,:), INTENT(inout)     :: laieff_fit     !! Fitted parameters for the effective LAI
    REAL(r_std),DIMENSION(:,:,:,:), INTENT(inout)          :: z_array_out    !! An output of h_array, to use in sechiba
    REAL(r_std),DIMENSION(:,:), INTENT(inout)              :: max_height_store!! ???
    REAL(r_std),DIMENSION(:,:,:), INTENT(inout)            :: Pgap_cumul     !! The probability of finding a gap in the in canopy from the top 
                                                                             !! of the canopy to a given level (unitless, between 0-1)

!! 0.4 Local variables
    INTEGER(i_std)                                         :: jsl, jv, ji
    REAL(r_std),DIMENSION (kjpindex,nslm)                  :: land_frac      !! To ouput the clay/sand/silt fractions with a vertical dim
    REAL(r_std), DIMENSION(ncirc)                          :: circ_height    !! temporary variable for vegetation height
    REAL(r_std), DIMENSION(ncirc)                          :: circ_dia       !! temporary variable for vegetation diameter
    REAL(r_std), DIMENSION(ncirc)                          :: circ_n         !!Temporary variable for stand density (number of trees)

!_ ================================================================================================================================

    !! Initialize local printlev
    printlev_loc=get_printlev('slowproc')
    

    IF(printlev_loc>=5) WRITE(numout,*) 'Entering slowproc_initialize'

    !! 1. Perform the allocation of all variables, define some files and some flags. 
    !     Restart file read for Sechiba.
    CALL slowproc_init (kjit, kjpindex, IndexLand, lalo, neighbours, resolution, &
         contfrac, rest_id, frac_age, veget, frac_nobio, totfrac_nobio, soiltile, &
         fraclut, nwdfraclut, reinf_slope, &
         ks,  nvan, avan, mcr, mcs, mcfc, mcw, &
         veget_max, tot_bare_soil, njsc, &
         ninput_update, ninput_year, &
         circ_class_biomass, circ_class_n, assim_param, loss_gain, veget_max_new, &
         frac_nobio_new, fraction_aeirrig_sw, precip_longterm)

    !! 2. Define time step in days for stomate
    dt_days = dt_stomate / one_day
    
    !! 3. check time step coherence between slow processes and fast processes
    IF ( dt_stomate .LT. dt_sechiba ) THEN
       WRITE(numout,*) 'slow_processes: time step smaller than forcing time step, dt_sechiba=',dt_sechiba,' dt_stomate=',dt_stomate
       CALL ipslerr_p(3,'slowproc_initialize',&
            'Coherence problem between dt_stomate and dt_sechiba',&
            'Time step smaller than forcing time step','')
    ENDIF
    
    
    !! 4. Call stomate to initialize all variables managed in stomate.
    IF ( ok_stomate ) THEN


       IF(printlev>=5) WRITE(numout,*) 'Entering stomate_initialize in slowproc'

       !  Note that some of the variables are stored in both sechiba and
       !  stomate. The values of these variables should be identical in 
       !  the stomate and sechiba restart. Anyway, the values found in
       !  sechiba will be overwritten by those found in stomate.
       CALL stomate_initialize (kjit,         kjpij,             kjpindex,     &
            rest_id_stom,   hist_id_stom,     hist_id_stom_IPCC,               &
            indexLand,      lalo,             neighbours,        resolution,   &
            contfrac,       clayfraction,     siltfraction,                    &
            bulk,           temp_air,         veget,             veget_max,    & 
            deadleaf_cover, assim_param,                                       &
            circ_class_biomass, circ_class_n, lai_per_level,     laieff_fit,   &
            temp_growth,                                                       &
            som_total,      heat_Zimov,       altmax,            depth_organic_soil,&
            cn_leaf_init_2D)

    ELSE
    
       ! The model is not using stomate but it needs to calculate a 
       ! vegetation structure anyway to calculate photosynthesis,
       ! albedo, and the energy budget
       CALL stomate_laieff_initialize()

       IF (read_lai .EQ. 1) THEN 
          lai(:,:)= 0.
          lai_per_level(:,:,:)= 0.
          laieff_fit(:,:,:)%a=zero
          laieff_fit(:,:,:)%b=zero
          laieff_fit(:,:,:)%c=zero
          laieff_fit(:,:,:)%d=zero
          laieff_fit(:,:,:)%e=zero
          circ_class_biomass(:,:,:,:,:) = 0.
          circ_class_n(:,:,:) = 0.
          z_array_out(:,:,:,:) = 0.
          height(:,:) = zero
          height_dom(:,:) = zero
          veget(:,:) = zero
       ENDIF

       ! Set a fixed value for altmax corresponding to the total soil depth for the hydrology
       altmax(:,:) = zdr(nslm)

       ! Initialize depth_organic_soil
       depth_organic_soil(:) = 0.0

    ENDIF

    !! 5. Vegetation structure can be read from sechiba and stomate
    !  All key variable should be initialized, read from a restart or imposed
    !  by the time this part of the code is reached. These variables now have
    !  to be used to calculate some derived variables that are required by
    !  sechiba. For example, circ_class_biomass and circ_class_n are used to
    !  calculate lai_effit which is required to calculate the albedo but is not
    !  stored in the sechiba restart file. veget_max may come from a map which 
    !  implies that there might be small imprecisions. This routine also calculates 
    !  veget. Update veget, since the biomass has possibly been read in from a 
    !  restart file, taken from a map or prescribed by an imposed value.
    CALL slowproc_canopy (kjpindex, circ_class_biomass, circ_class_n, &
            veget_max, lai_per_level, z_array_out, &
            max_height_store, laieff_fit, frac_age)

    CALL slowproc_veget (kjpindex, lai_per_level, z_array_out, &
         coszang_noon, circ_class_biomass, circ_class_n, frac_nobio, totfrac_nobio, &
         veget_max, veget, soiltile, tot_bare_soil, fraclut, nwdFraclut, Pgap_cumul, &
         precip_longterm)

    !! 6. Calculate height and lai from biomass
    
    ! Dynamic parameters
    CALL dynamic_parameters_calculate(kjpindex, precip_longterm)

    IF ( read_lai .NE. 1 ) THEN
       height(:,ibare_sechiba) = zero
       height_dom(:,ibare_sechiba) = zero
       height_inv(:,ibare_sechiba) = zero
       dia_dom(:,ibare_sechiba) = zero
       dia_inv(:,ibare_sechiba) = zero
       ba_inv(:,ibare_sechiba) = zero
       ind_inv(:,ibare_sechiba) = zero
       lai(:,ibare_sechiba) = zero  
       DO jv = 2,nvm
          DO ji = 1,kjpindex
             ! Skip if veget_max = 0, else calculate height and lai
             IF (veget_max(ji,jv) .EQ. zero) THEN
                height(ji,jv) = zero
                height_dom(ji,jv) = zero
                height_inv(ji,jv) = zero
                dia_dom(ji,jv) = zero
                dia_inv(ji,jv) = zero
                ba_inv(ji,jv) = zero
                ind_inv(ji,jv) = zero
                lai(ji,jv) = zero
             ELSE
                height(ji,jv) = wood_to_qmheight(circ_class_biomass(ji,jv,:,:,icarbon), &
                     circ_class_n(ji,jv,:), jv, pipe_tune2(ji,jv))
                IF(is_tree(jv))THEN
                   ! Calculate diameter and height of the dominant diameter class
                   circ_height(:) = wood_to_height(circ_class_biomass(ji,jv,:,:,icarbon),jv,pipe_tune2(ji,jv))
                   height_dom(ji,jv) = circ_height(ncirc)
                   circ_dia(:) = wood_to_dia(circ_class_biomass(ji,jv,:,:,icarbon),jv,pipe_tune2(ji,jv))
                   dia_dom(ji,jv) = circ_dia(ncirc)
                   ! Calculate quadratic mean height and diameter of the trees
                   ! of which the diameter exceeds a threshold. This is closely
                   ! matches the height and diameter reported in forest 
                   ! inventories
                   circ_n(:) = circ_class_n(ji,jv,:)
                   WHERE (circ_dia(:) .lt. dia_thresh_inv(jv))
                      circ_dia(:) = zero
                      circ_height(:) = zero
                      circ_n(:) = zero
                   END WHERE
                   IF (SUM(circ_n(:)).EQ.zero) THEN
                      dia_inv(ji,jv) =  zero
                      height_inv(ji,jv) = zero
                      ba_inv(ji,jv) = zero
                      ind_inv(ji,jv) = zero
                   ELSE
                      dia_inv(ji,jv) = (SUM(circ_dia(:)**2*circ_n(:))/SUM(circ_n(:)))**0.5
                      height_inv(ji,jv) = (SUM(circ_height(:)**2*circ_n(:))/SUM(circ_n(:)))**0.5
                      ind_inv(ji,jv) = SUM(circ_n(:))
                      ba_inv(ji,jv) = pi/4*SUM(circ_dia(:)**2*circ_n(:))*m2_to_ha 
                   END IF
                ELSE
                   ! Grasses and crops have a height but no
                   ! diameer in ORCHIDEE
                   height_dom(ji,jv) = height(ji,jv)
                   height_inv(ji,jv) = height(ji,jv)
                   dia_dom(ji,jv) = undef
                   dia_inv(ji,jv) = undef
                   ba_inv(ji,jv) = undef
                   ind_inv(ji,jv) = circ_class_n(ji,jv,1)
                ENDIF
                lai(ji,jv) = cc_to_lai(circ_class_biomass(ji,jv,:,ileaf,icarbon), &
                     circ_class_n(ji,jv,:),jv)
             ENDIF
          ENDDO
       ENDDO
    ENDIF ! read_lai .NE. 1) THEN

    !! 3. Initialize the maximum water on vegetation for interception
    qsintmax(:,:) = qsintcst * veget(:,:) * lai(:,:) 
    qsintmax(:,1) = zero

    !! 5. Specific run without the carbon cycle (STOMATE not called): 
    !!    Need to initialize some variables that will be used in SECHIBA
    IF (.NOT. ok_stomate ) THEN
       ! Initialize some missing variables
       deadleaf_cover(:) = zero
       temp_growth(:)=25.
    ENDIF

    !! 5.1  Dynamic irrigation maps as output to sechiba_end
    irrigated_next(:) = zero
    irrig_frac_next(:) = zero
    IF (do_irrigation ) THEN
      irrigated_next(:) = irrigated_new(:)
      ! irrig_frac calculation
      DO ji=1,kjpindex
          IF( area(ji)*contfrac(ji) > min_sechiba) THEN
            !SUM(routing_area(ig,:)) is totarea(ig) = m2
            irrig_frac_next(ji) = MIN( soiltile(ji,irrig_st) * SUM(veget_max(ji,:)) , &
                irrigated_new(ji) / ( area(ji)*contfrac(ji) ) )

            ! Soiltile(fraction of vegtot) * SUM(veget_max) [SUM(veget_max) is vegtot in hydrol, calculated
            ! differently in routing] = fraction of grid cell
            ! irrigated(m2)/ SUM(routing_area(ig,:) = Fraction of grid cell
            ! irrig_frac is always fraction of grid cell
          ENDIF
      ENDDO
    END IF

    !! 5.2 Calculation of reinf_slope_soil to pass to hydrol
    reinf_slope_soil(:,:) = zero
    DO ji=1,kjpindex
        reinf_slope_soil(ji,:) = reinf_slope(ji)
        IF( Reinfiltr_IrrigField .AND.  irrig_frac_next(ji) > min_sechiba ) THEN
          reinf_slope_soil(ji, irrig_st) = MAX(reinf_slope(ji), reinf_slope_cropParam)
        ENDIF
    ENDDO

    !! 6. Output with XIOS for variables done only once per run
    DO jsl=1,nslm
       land_frac(:,jsl) = clayfraction(:)
    ENDDO
    ! mean fraction of clay in grid-cell
    CALL xios_orchidee_send_field("clayfraction",land_frac)
    DO jsl=1,nslm
       land_frac(:,jsl) = sandfraction(:)
    ENDDO
    ! mean fraction of sand in grid-cell
    CALL xios_orchidee_send_field("sandfraction",land_frac) 
    DO jsl=1,nslm
       land_frac(:,jsl) = siltfraction(:)
    ENDDO
    ! mean fraction of silt in grid-cell
    CALL xios_orchidee_send_field("siltfraction",land_frac) 
    
    IF(printlev>=5) WRITE(numout,*) 'Leaving slowproc_initialize'

  END SUBROUTINE slowproc_initialize


!! ================================================================================================================================
!! SUBROUTINE   : slowproc_main
!!
!>\BRIEF         Main routine that manage variable initialisation (slowproc_init), 
!! prepare the restart file with the slowproc variables, update the time variables 
!! for slow processes, and possibly update the vegetation cover, before calling 
!! STOMATE in the case of the carbon cycle activated or just update the canopy (and possibly
!! the vegetation cover) for simulation with only SECHIBA   
!!
!!
!! DESCRIPTION  : (definitions, functional, design, flags): The subroutine manages 
!! diverses tasks:
!! (1) Initializing all variables of slowproc (first call)
!! (2) Preparation of the restart file for the next simulation with all prognostic variables
!! (3) Compute and update time variable for slow processes
!! (4) Update the vegetation cover if there is some land use change (only every years)
!! (5) Call STOMATE for the runs with the carbone cycle activated (ok_stomate) and compute the respiration
!!     and the net primary production
!! (6) Compute the LAI and possibly update the vegetation cover for run without STOMATE 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S):  ::fco2_flux_out, ::fco2_lu_out,::fco2_wh_out, ::fco2_ha_out, ::lai, ::height, ::veget, ::frac_nobio,  
!! ::veget_max, ::woodharvest, ::totfrac_nobio, ::soiltype, ::assim_param, ::deadleaf_cover, ::qsintmax,
!! and resp_maint, resp_hetero, resp_growth, npp that are calculated and stored
!! in stomate is activated.  
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : 
! \latexonly 
! \includegraphics(scale=0.5){SlowprocMainFlow.eps} !PP to be finalize!!)
! \endlatexonly
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_main (kjit, kjpij, kjpindex, njsc, &
       IndexLand, indexveg, lalo, neighbours, resolution, contfrac, soiltile, fraclut, nwdFraclut, &
       temp_air, temp_sol, stempdiag, precip_longterm, &
       vegstress, humrel, &
       shumdiag, litterhumdiag, precip_rain, precip_snow, pb, gpp, VPD, &
       tmc_pft, drainage_pft, runoff_pft, swc_pft, deadleaf_cover, &
       assim_param, qsintveg, &
       frac_age, height, lai, veget, frac_nobio, veget_max, totfrac_nobio, qsintmax, &
       rest_id, hist_id, hist2_id, rest_id_stom, hist_id_stom, hist_id_stom_IPCC, &
       fco2_flux_out, fco2_lu_out, fco2_wh_out, fco2_ha_out, &
       temp_growth, tot_bare_soil, &
       tdeep, hsdeep_long, snow, heat_Zimov, &
       sfluxCH4_deep, sfluxCO2_deep, &
       som_total,snowdz,snowrho, altmax, depth_organic_soil, &
       circ_class_biomass, circ_class_n, &
       lai_per_level, max_height_store, laieff_fit, &
       Light_Abs_Tot, Light_Tran_Tot, laieff_isotrop, &
       z_array_out, transpir, transpir_mod, &
       coszang, coszang_noon, &
       stressed, unstressed, &
       u, v, loss_gain, veget_max_new, frac_nobio_new, failed_vegfrac, &
       mcs_hydrol, mcfc_hydrol, & 
       vessel_loss, height_dom, height_inv, dia_dom, dia_inv, &
       ba_inv, ind_inv, root_profile, root_depth, us, Pgap_cumul, z0m, &
       irrigated_next, irrig_frac_next, reinf_slope, reinf_slope_soil, &
       wtp,  mc_peat_above, & 
       liqwt_ratio,  shumdiag_peat,  shumdiag_croppeat,&
       mc_croppeat_above, &
       shumdiag_man,  mc_man_above)

!! INTERFACE DESCRIPTION

!! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                          :: kjit                !! Time step number
    INTEGER(i_std), INTENT(in)                          :: kjpij               !! Total size of the un-compressed grid
    INTEGER(i_std),INTENT(in)                           :: kjpindex            !! Domain size - terrestrial pixels only
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)    :: njsc                !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless)
    INTEGER(i_std),INTENT (in)                          :: rest_id,hist_id     !! _Restart_ file and _history_ file identifier
    INTEGER(i_std),INTENT (in)                          :: hist2_id            !! _history_ file 2 identifier
    INTEGER(i_std),INTENT (in)                          :: rest_id_stom        !! STOMATE's _Restart_ file identifier
    INTEGER(i_std),INTENT (in)                          :: hist_id_stom        !! STOMATE's _history_ file identifier
    INTEGER(i_std),INTENT(in)                           :: hist_id_stom_IPCC   !! STOMATE's IPCC _history_ file identifier
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)    :: IndexLand           !! Indices of the points on the land map
    INTEGER(i_std),DIMENSION (kjpindex*nvm), INTENT (in):: indexveg            !! Indices of the points on the vegetation (3D map ???) 
    REAL(r_std),DIMENSION (kjpindex,2), INTENT (in)     :: lalo                !! Geogr. coordinates (latitude,longitude) (degrees)
    INTEGER(i_std), DIMENSION (kjpindex,NbNeighb), INTENT(in)  :: neighbours   !! neighbouring grid points if land
    REAL(r_std), DIMENSION (kjpindex,2), INTENT(in)     :: resolution          !! size in x an y of the grid (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: contfrac            !! Fraction of continent in the grid (0-1, unitless)
    REAL(r_std), DIMENSION(kjpindex), INTENT(in)        :: temp_air            !! Temperature of first model layer (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: temp_sol            !! Surface temperature (K)
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (in)  :: stempdiag           !! Soil temperature (K)
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (in)  :: shumdiag            !! Relative soil moisture (0-1, unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: litterhumdiag       !! Litter humidity  (0-1, unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: precip_rain         !! Rain precipitation (mm dt_stomate^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: precip_snow         !! Snow precipitation (mm dt_stomate^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: pb                  !! Lowest level pressure (Pa) 
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(in)    :: gpp                 !! GPP of total ground area (gC m^{-2} time step^{-1}). 
                                                                               !! Calculated in sechiba, account for vegetation cover and 
                                                                               !! effective time step to obtain gpp_d
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)       :: VPD                 !! Vapor Pressure Deficit (kPa)   
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(in)   :: tmc_pft             !! Total soil water per PFT (mm/m2) 
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(in)   :: drainage_pft        !! Drainage per PFT (mm/m2)   
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(in)   :: runoff_pft          !! Drainage per PFT (mm/m2)  
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(in)   :: swc_pft             !! Relative Soil water content [tmcr:tmcs] per pft (-)     
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)   :: transpir            !! transpiration @tex $(kg m^{-2} days^{-1})$ @endtex
    REAL(r_std), DIMENSION(kjpindex,nvm)                :: transpir_mod        !! transpir converted from mm/day to mm dt^(-1)

    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: u                   !! Lowest level wind speed in direction u (m/s)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: v                   !! Lowest level wind speed in direction v (m/s)


    REAL(r_std),DIMENSION(kjpindex),INTENT(in)                   :: coszang        !! the cosine of the solar zenith angle (unitless)  
    REAL(r_std),DIMENSION(kjpindex),INTENT(in)                   :: coszang_noon   !! the cosine of the solar zenith angle (unitless)  
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm),   INTENT (in)    :: tdeep          !! deep temperature profile (K)
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm),   INTENT (in)    :: hsdeep_long    !! deep long term soil humidity profile
    REAL(r_std), DIMENSION(kjpindex),         INTENT (in)        :: snow           !! Snow mass [Kg/m^2]
    REAL(r_std), DIMENSION(kjpindex,nsnow),INTENT(in)            :: snowdz         !! snow depth for each layer [m]
    REAL(r_std), DIMENSION(kjpindex,nsnow),INTENT(in)            :: snowrho        !! snow density for each layer (Kg/m^3)
    REAL(r_std), DIMENSION(kjpindex), INTENT(in)                 :: mcs_hydrol     !! Saturated volumetric water content output to be used in stomate_soilcarbon
    REAL(r_std), DIMENSION(kjpindex), INTENT(in)                 :: mcfc_hydrol    !! Volumetric water content at field capacity output to be used in stomate_soilcarbon
    REAL(r_std), DIMENSION(kjpindex,nvm,ncirc), INTENT(in)       :: vessel_loss    !! Proportion  of conductivity lost due to cavitation in the xylem (no unit).
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)                  :: root_profile   !! Normalized root mass/length fraction in each soil layer 
                                                                                   !! (0-1, unitless)
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)                    :: root_depth     !! Node and interface numbers at which the deepest roots
                                                                                   !! occur (1 to nslm, unitless)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)                 :: Pgap_cumul     !! The probability of finding a gap in the in canopy from the top 
                                                                                   !! of the canopy to a given level (unitless, between 0-1)
    REAL(r_std), DIMENSION(:), INTENT(in)                        :: z0m            !! Surface roughness for momentum (m)
    REAL(r_std),DIMENSION(kjpindex), INTENT(in)                  :: reinf_slope    !! slope coef for reinfiltration
    !
    ! peatland and tides
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)               :: liqwt_ratio      !! liquid water ratio for peat, compared to the saturated one
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)               :: wtp              !! water table position when peat is activated
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)               :: mc_peat_above    !! liquid soil moisture content within the first 4 layers when peat is activated
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)               :: mc_croppeat_above!! liquid soil moisture content within the first 4 layers when agri_peat is activated
    REAL(r_std), DIMENSION(kjpindex,nslm),   INTENT (in)         :: shumdiag_peat    !! soil moisture content useful for the decomposition when peat is activated   
    REAL(r_std), DIMENSION(kjpindex,nslm),   INTENT (in)         :: shumdiag_croppeat!! soil moisture content useful for the decomposition when agri_peat is activated 
    REAL(r_std), DIMENSION(kjpindex,nslm),   INTENT (in)         :: shumdiag_man     !! soil moisture content useful for the decomposition when tides is activated 
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)               :: mc_man_above     !! liquid soil moisture content within the first 4 layers when tides is activated
    !

!! 0.2 Output variables 

    REAL(r_std), DIMENSION (kjpindex), INTENT(out)      :: fco2_flux_out       !! CO2 flux per average ground area (gC m^{-2} s{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)      :: fco2_lu_out         !! CO2 flux from land-use (without forest management) (gC m^{-2} dt_stomate^{-1})
    REAL(r_std), DIMENSION (kjpindex), INTENT (out)     :: fco2_wh_out         !! CO2 Flux to Atmosphere from Wood Harvesting (gC m^{-2} dt_stomate^{-1})
    REAL(r_std), DIMENSION (kjpindex), INTENT (out)     :: fco2_ha_out         !! CO2 Flux to Atmosphere from Crop Harvesting (gC m^{-2} dt_stomate^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)      :: temp_growth         !! Growth temperature (Â°C) - Is equal to t2m_month 
    REAL(r_std), DIMENSION (kjpindex), INTENT(out)      :: tot_bare_soil       !! Total evaporating bare soil fraction in the mesh
    REAL(r_std),DIMENSION(kjpindex,nvm),INTENT(out)     :: max_height_store    !! ???

    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm),INTENT (out)   :: heat_Zimov    !! heating associated with decomposition [W/m**3 soil] 
    REAL(r_std), DIMENSION(kjpindex),     INTENT (out)  :: sfluxCH4_deep       !! surface flux of CH4 to atmosphere from permafrost
    REAL(r_std), DIMENSION(kjpindex),     INTENT (out)  :: sfluxCO2_deep       !! surface flux of CO2 to atmosphere from permafrost
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(out)   :: lai                 !! Leaf area index (m^2 m^{-2})
    REAL(r_std), DIMENSION (kjpindex), INTENT(out)      :: irrigated_next      !! Dynamic irrig. area, calculated in slowproc and passed to routing
    

!! 0.3 Modified variables

    REAL(r_std),DIMENSION (kjpindex,nvm,nleafages), INTENT(inout):: frac_age   !! Age efficacity from STOMATE for isoprene
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout)     :: veget          !! Fraction of vegetation type including none biological fractionin the mesh (unitless)
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (inout)  :: frac_nobio     !! Fraction of ice, lakes, cities etc. in the mesh
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout)     :: veget_max      !! Maximum fraction of vegetation type in the mesh (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)         :: totfrac_nobio  !! Total fraction of ice+lakes+cities etc. in the mesh
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT(inout)    :: soiltile       !! Fraction of each soil tile within vegtot (0-1, unitless)
    REAL(r_std),DIMENSION (kjpindex,nvm,npco2),INTENT (inout):: assim_param    !! vcmax, nue and leaf N for photosynthesis 
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(inout)     :: qsintveg       !! Water on vegetation due to interception @tex $(kg m^{-2})$ @endtex
    REAL(r_std), DIMENSION (kjpindex,nlut), INTENT(inout)    :: fraclut        !! Fraction of each landuse tile (0-1, unitless)
    REAL(r_std), DIMENSION (kjpindex,nlut), INTENT(inout)    :: nwdFraclut     !! Fraction of non-woody vegetation in each landuse tile (0-1, unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)         :: deadleaf_cover !! Fraction of soil covered by dead leaves (unitless)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout)     :: qsintmax       !! Maximum water storage on vegetation from interception (mm)
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm,nelements), INTENT (inout) :: som_total !! Total soil carbon for use in thermal (g/m**3)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(inout)       :: altmax        !! Maximul active layer thickness (m). Be careful, here active means non frozen.
                                                                               !! Not related with the active soil carbon pool.
    REAL(r_std), DIMENSION(kjpindex), INTENT (inout)              :: depth_organic_soil !! how deep is the organic soil?
    REAL(r_std),DIMENSION(kjpindex,nvm,ncirc,nparts,nelements), &
                                               INTENT(inout)      :: circ_class_biomass !!  Biomass components of the model tree  
                                                                               !! within a circumference class
                                                                               !! class @tex $(g C ind^{-1})$ @endtex
    REAL(r_std),DIMENSION(kjpindex,nvm,ncirc), INTENT(inout)      :: circ_class_n !! Number of trees within each circumference
                                                                               !! class @tex $(m^{-2})$ @endtex
    REAL(r_std),DIMENSION(kjpindex,nvm), INTENT(inout)            :: loss_gain !! Changes in veget_max distributed over all
                                                                               !! age classes and thus taking the age-classes into 
                                                                               !! account (unitless, 0-1)
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT(inout)        :: frac_nobio_new !! Fraction of ice,lakes,cities, ... (unitless)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(inout)           :: veget_max_new !! Maximum fraction of vegetation type including none 
    REAL(r_std),DIMENSION(:,:,:), INTENT(inout)                   :: lai_per_level !! This is the LAI per vertical level
                                                                               !! @tex $(m^{2} m^{-2})$
    TYPE(laieff_type),DIMENSION (:,:,:), INTENT(inout)            :: laieff_fit!! Fitted parameters for the effective LAI
    REAL(r_std),DIMENSION (:,:,:), INTENT (inout)                 :: Light_Abs_Tot !!Absorbed radiation per level for photosynthesis
    REAL(r_std),DIMENSION (:,:,:), INTENT (inout)                 :: Light_Tran_Tot !!Transmitted radiation per level for photosynthesis
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)                  :: laieff_isotrop!! Effective LAI
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(inout)           :: stressed  !! adjusted ecosystem functioning. Takes the unit of the variable
                                                                               !! used as a proxy for waterstress (assigned in sechiba)
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(inout)           :: unstressed!! initial ecosystem functioning after the first calculation and 
                                                                               !! before any recalculations. Takes the unit of the variable used
                                                                               !! as a proxy for unstressed.

    REAL(r_std),DIMENSION(:,:,:,:), INTENT(inout)                 :: z_array_out !! An output of h_array, to use in sechiba
    LOGICAL, DIMENSION(kjpindex,nvm), INTENT(inout)               :: failed_vegfrac !! Failed to find a PFT were some residual fraction could be added (true/false)
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT (inout)         :: vegstress !! Relative soil moisture used in stomate to calculate plant water stress (0-1, unitless)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout)          :: humrel    !! Relative humidity - not used in stomate (needed in age_class_distr)
    REAL(r_std),DIMENSION (kjpindex,nvm,nstm,nslm), INTENT(inout) :: us        !! Water stress index for transpiration
                                                                               !! (by soil layer and PFT) (0-1, unitless)
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(inout)           :: height    !! height of vegetation (m)
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(inout)           :: height_dom!! dominant height of vegetation (m)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(inout)           :: dia_dom   !! dominant diameter of vegetation (m)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(inout)           :: height_inv!! height of vegetation as measured in a forest inventory (m)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(inout)           :: dia_inv   !! diameter of vegetation as measured in a forest inventory (m)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(inout)           :: ba_inv    !! ba of vegetation as measured in a forest inventory (m2)    
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(inout)           :: ind_inv   !! stand density as measured in a forest inventory (trees m-1)
    REAL(r_std), DIMENSION (kjpindex), INTENT (inout)             :: irrig_frac_next !! Dynamic irrig. fraction, calculated in slowproc and passed to routing
    REAL(r_std),DIMENSION (kjpindex, nstm), INTENT(inout)         :: reinf_slope_soil!!  slope coef for reinfiltration per soil tile
    REAL(r_std),DIMENSION (kjpindex), INTENT(inout)               :: precip_longterm !! Longterm annual precipitation sum mm year^{-1} 


!! 0.4 Local variables
    INTEGER(i_std)                                     :: j, jv, ji,jj         !! indices (unitless)
    INTEGER(i_std)                                     :: ipts,ivm,jvm         !! indices (unitless)
    INTEGER(i_std)                                     :: igroup               !! Indices (unitless)
    REAL(r_std), DIMENSION(kjpindex,nvm)               :: resp_maint           !! Maitanance component of autotrophic respiration in (gC m^{-2} dt_stomate^{-1})
    REAL(r_std), DIMENSION(kjpindex,nvm)               :: resp_hetero          !! heterotrophic resp. (gC/(m**2 of total ground)/time step)
    REAL(r_std), DIMENSION(kjpindex,nvm)               :: resp_growth          !! Growth component of autotrophic respiration in gC m^{-2} dt_stomate^{-1})
    REAL(r_std), DIMENSION(kjpindex,nvm)               :: npp                  !! Net Ecosystem Exchange (gC/(m**2 of total ground)/time step)
    REAL(r_std),DIMENSION (kjpindex)                   :: totfrac_nobio_new    !! Total fraction for the next year
    REAL(r_std),DIMENSION (kjpindex)                   :: histvar              !! Temporary variable for output
    REAL(r_std), DIMENSION(kjpindex,nvm,12)            :: N_input_temp              
    CHARACTER(LEN=80)                                  :: fieldname            !! name of the field read in the N input map
    REAL(r_std), DIMENSION(kjpindex,nvm)               :: veget_max_temp
    REAL(r_std), DIMENSION(kjpindex,nvm)               :: veget_ny_map
    INTEGER(i_std)                                     :: nfound               !! The number of PFTs for this pixel which have veget_max
    LOGICAL(r_std),DIMENSION (nvm)                     :: lfound_veget         !! Do we have a given threshold of veget_max on for this pixel?
    REAL(r_std)                                        :: excess_veg           !! The veget_max which we have to redistribute t0 other pixels.
    REAL(r_std)                                        :: nobio_frac_sum_old   !! Total vegetated land for the current year.
    REAL(r_std)                                        :: nobio_frac_sum_new   !! Total vegetated land for next year.
    REAL(r_std)                                        :: nobio_diff           !! How much the nonvegetative land coverage changes between the current and next year.i
    REAL(r_std)                                        :: nobio_scale          !! This scaling factor is trying to redistribute the nobio_diff to new one
    CHARACTER(30)                                      :: losses               !! Loss distribution across age classes
    REAL(r_std), DIMENSION(ncirc)                      :: circ_height          !! temporary variable for vegetation height
    REAL(r_std), DIMENSION(ncirc)                      :: circ_dia             !! temporary variable for vegetation diameter
    REAL(r_std), DIMENSION(ncirc)                      :: circ_n               !! Temporary variable for stand density (number of trees)
    
!_ ================================================================================================================================

    
    !! 1. Compute and update all variables linked to the date and time
    IF (printlev_loc>=4) WRITE(numout,*) 'Entering slowproc_main, year_start, month_start, day_start, sec_start=',&
          year_start, month_start,day_start,sec_start

    !! 2. Activate slow processes if it is the end of the day
    IF ( LastTsDay ) THEN

       ! 3.2.2 Activate slow processes in the end of the day
       do_slow = .TRUE.

       ! 3.2.3 Count the number of days 
       days_since_beg = days_since_beg + 1
       IF (printlev_loc>=4) WRITE(numout,*) "New days_since_beg : ",days_since_beg

    ELSE

       do_slow = .FALSE.

    ENDIF
    
    !! 3. Update the vegetation if it is time to do so.
    !     This is done at the first sechiba time step on a new year if
    !     veget_update > 0. veget_update can only be 0 or 1. 
    !     Nothing is done if veget_update=0.
    !     Update will never be done if impveg=true because veget_update=0.
    IF ( FirstTsYear ) THEN

       IF (veget_update > 0) THEN

          IF (printlev_loc>=1) WRITE(numout,*)  'We are updating the vegetation map'

          veget_max_new(:,:) = zero
          frac_nobio_new(:,:) = zero
          loss_gain(:,:) = zero
          failed_vegfrac(:,:) = .FALSE.

          IF (hack_veget_max_new) THEN
             ! Read veget_max_new from run.def. This is intended to be used in
             ! combination with a restart file. It only works once per simulation.
             ! It can be used to debug simplified LCC test cases. Tested use:
             ! (1) run the model for 1 year with LCC = true and making use of the
             ! pft-map. For this run hack_lcc = FALSE. (2) restart the model with
             ! hack_lcc = TRUE and specify the frac_nobio_new and veget_max_new 
             ! in the run.def. The model will apply the prescribed LCC the first 
             ! day of year 2. Note that it will read veget_max_new from the run.def 
             ! in all subsequent years as well. Hence there will be no land cover 
             ! changes except in the year 2.
             CALL getin_p('FRAC_NOBIO_NEW',frac_nobio_new)
             CALL getin_p('VEGET_MAX_NEW', veget_max_new)
             WRITE(numout,*) 'WARNING: Hacked land cover change'
             WRITE(numout,*) 'frac_nobio_new, ', frac_nobio_new
             DO jv=1,nvm
                WRITE(numout,*) 'PFT, veget_max_new, ', jv, veget_max_new(test_grid,jv)
             END DO
          ELSE
             ! Read the new the vegetation from file. Output is veget_max_new and
             ! frac_nobio_new
             ! 
             IF ( .NOT. ok_dgvm .AND. (agri_peat)) THEN
                CALL slowproc_readvegetmax(kjpindex, lalo, neighbours, resolution, contfrac, &
                     veget_max, veget_max_new, frac_nobio_new, .TRUE.)
             ELSE
                CALL slowproc_readvegetmax(kjpindex, lalo, neighbours, resolution, contfrac, &
                     veget_max, veget_max_new, frac_nobio_new, .FALSE.)
             ENDIF
          ENDIF

          ! Check the map that was just read for errors. After this check the
          ! fraction (inc. frac_nobio) should add up to one. frac_nobio_new
          ! is kept as "read" from the map. 
          CALL check_read_vegetmax(kjpindex, veget_max_new, frac_nobio_new)

          ! Compare the new vegetation map with the current map and resolve
          ! possible conflicts to make sure all will go smooth in land cover change.
          losses = 'proportional'
          CALL adjust_delta_veget_max(kjpindex, veget_max, frac_nobio, &
               veget_max_new, frac_nobio_new, loss_gain, losses, failed_vegfrac)

          ! Verification and correction of veget_max_new, calculation
          ! of veget and soiltile. totfrac_nobio, veget, soiltile and
          ! tot_bare_soil are output variables.
          CALL slowproc_veget (kjpindex, lai_per_level, z_array_out, &
               coszang_noon,circ_class_biomass, circ_class_n, frac_nobio_new, totfrac_nobio_new, &
               veget_max_new, veget, soiltile, tot_bare_soil, &
               fraclut, nwdFraclut, Pgap_cumul, precip_longterm)

          ! Set the flag do_now_stomate_lcchange to activate stomate_lcchange.
          ! This flag will be kept to true until stomate_lcchange has been done. 
          ! The variable totfrac_nobio_new will only be used in stomate when this 
          ! flag is activated
          do_now_stomate_lcchange=.TRUE.
          IF ( .NOT. ok_stomate ) THEN
             ! Special case if stomate is not activated : set the variable 
             ! done_stomate_lcchange=true so that the subroutine slowproc_change_frac 
             ! will be called in the end of sechiba_main.
             done_stomate_lcchange=.TRUE.
          END IF

       ENDIF ! Veget_update>0
       
    END IF ! FirstTsYear

    !! Read the Nitrogen inputs
    IF(ok_ncycle .AND. (.NOT. impose_CN)) THEN

       IF ( (Ninput_update > 0) .AND. FirstTsYear ) THEN
          ! Update of the vegetation cover with Land Use only if 
          ! the current year match the requested condition (a multiple of "veget_update")
          Ninput_year = Ninput_year + 1
          IF ( MOD(Ninput_year - Ninput_year_orig, Ninput_update) == 0 ) THEN
             IF (printlev_loc>=1) WRITE(numout,*)  'We are updating the Ninputs map for year =' , Ninput_year
             
             IF(.NOT. impose_ninput_dep) THEN
                ! Read the new N inputs from file. Output is Ninput and frac_nobio_nextyear.
                fieldname='Nammonium'
                CALL slowproc_Ninput(kjpindex, lalo, neighbours, resolution, contfrac, fieldname, &
                     N_input_temp, Ninput_year)
                ! Conversion from mgN/m2/yr to gN/m2/day
                N_input_raw(:,:,:,iammonium)=(N_input_temp + boost_fert/2)/1000./one_year

                fieldname='DRYNHX'
                CALL slowproc_Ninput(kjpindex, lalo, neighbours, resolution, contfrac, fieldname, &
                     N_input_temp, Ninput_year)
                ! Conversion from kgN/m2/s to gN/m2/day
                N_input_raw(:,:,:,iammonium)=N_input_raw(:,:,:,iammonium)+N_input_temp*1000.*one_day

                fieldname='WETNHX'
                CALL slowproc_Ninput(kjpindex, lalo, neighbours, resolution, contfrac, fieldname, &
                     N_input_temp, Ninput_year)
                ! Conversion from kgN/m2/s to gN/m2/day
                N_input_raw(:,:,:,iammonium)=N_input_raw(:,:,:,iammonium)+N_input_temp*1000.*one_day

                fieldname='Nnitrate'
                CALL slowproc_Ninput(kjpindex, lalo, neighbours, resolution, contfrac, fieldname, &
                     N_input_temp, Ninput_year)
                ! Conversion from mgN/m2/yr to gN/m2/day
                N_input_raw(:,:,:,initrate)=(N_input_temp + boost_fert/2)/1000./one_year

                fieldname='DRYNOY'
                CALL slowproc_Ninput(kjpindex, lalo, neighbours, resolution, contfrac, fieldname, &
                     N_input_temp, Ninput_year)
                ! Conversion from kgN/m2/s to gN/m2/day
                N_input_raw(:,:,:,initrate)=N_input_raw(:,:,:,initrate)+N_input_temp*1000.*one_day

                fieldname='WETNOY'
                CALL slowproc_Ninput(kjpindex, lalo, neighbours, resolution, contfrac, fieldname, &
                     N_input_temp, Ninput_year)
                ! Conversion from kgN/m2/s to gN/m2/day
                N_input_raw(:,:,:,initrate)=N_input_raw(:,:,:,initrate)+N_input_temp*1000.*one_day

             ENDIF

             IF(.NOT. impose_ninput_fert) THEN
                fieldname='Nfert'
                CALL slowproc_Ninput(kjpindex, lalo, neighbours, resolution, contfrac, fieldname, &
                     N_input_temp, Ninput_year)
                ! Conversion from gN/m2(cropland)/yr to gN/m2/day
                N_input_raw(:,:,:,ifert_ammo) = N_input_temp(:,:,:)/one_year * ratio_nh4_fert 
                N_input_raw(:,:,:,ifert_nitr) = N_input_temp(:,:,:)/one_year * (1.-ratio_nh4_fert)

                fieldname='Nfert_cropland'
                CALL slowproc_Ninput(kjpindex, lalo, neighbours, resolution, contfrac, fieldname, &
                     N_input_temp, Ninput_year)
                ! Conversion from gN/m2(cropland)/yr to gN/m2/day
                N_input_raw(:,:,:,ifert_ammo) = N_input_raw(:,:,:,ifert_ammo)+ N_input_temp(:,:,:)/one_year * ratio_nh4_fert
                N_input_raw(:,:,:,ifert_nitr) = N_input_raw(:,:,:,ifert_nitr)+ N_input_temp(:,:,:)/one_year * (1.-ratio_nh4_fert)

                fieldname='Nfert_ammo_cropland'
                CALL slowproc_Ninput(kjpindex, lalo, neighbours, resolution, contfrac, fieldname, &
                     N_input_temp, Ninput_year)
                ! Conversion from gN/m2(cropland)/yr to gN/m2/day
                N_input_raw(:,:,:,ifert_ammo) = N_input_raw(:,:,:,ifert_ammo)+ N_input_temp(:,:,:)/one_year 

                fieldname='Nfert_nitr_cropland'
                CALL slowproc_Ninput(kjpindex, lalo, neighbours, resolution, contfrac, fieldname, &
                     N_input_temp, Ninput_year)
                ! Conversion from gN/m2(cropland)/yr to gN/m2/day
                N_input_raw(:,:,:,ifert_nitr) = N_input_raw(:,:,:,ifert_nitr)+ N_input_temp(:,:,:)/one_year 


                fieldname='Nfert_cropC3'
                CALL slowproc_Ninput(kjpindex, lalo, neighbours, resolution, contfrac, fieldname, &
                     N_input_temp, Ninput_year)
                ! Conversion from gN/m2(cropland)/yr to gN/m2/day
                N_input_raw(:,:,:,ifert_ammo) = N_input_raw(:,:,:,ifert_ammo)+ N_input_temp(:,:,:)/one_year * ratio_nh4_fert
                N_input_raw(:,:,:,ifert_nitr) = N_input_raw(:,:,:,ifert_nitr)+ N_input_temp(:,:,:)/one_year * (1.-ratio_nh4_fert)

 
                fieldname='Nfert_cropC4'
                CALL slowproc_Ninput(kjpindex, lalo, neighbours, resolution, contfrac, fieldname, &
                     N_input_temp, Ninput_year)
                ! Conversion from gN/m2(cropland)/yr to gN/m2/day
                N_input_raw(:,:,:,ifert_ammo) = N_input_raw(:,:,:,ifert_ammo)+ N_input_temp(:,:,:)/one_year * ratio_nh4_fert
                N_input_raw(:,:,:,ifert_nitr) = N_input_raw(:,:,:,ifert_nitr)+ N_input_temp(:,:,:)/one_year * (1.-ratio_nh4_fert)

                fieldname='Nfert_pasture'
                CALL slowproc_Ninput(kjpindex, lalo, neighbours, resolution, contfrac, fieldname, &
                     N_input_temp, Ninput_year)
                ! Conversion from gN/m2(pasture)/yr to gN/m2/day
                N_input_raw(:,:,:,ifert_ammo) = N_input_raw(:,:,:,ifert_ammo)+ N_input_temp(:,:,:)/one_year * ratio_nh4_fert
                N_input_raw(:,:,:,ifert_nitr) = N_input_raw(:,:,:,ifert_nitr)+ N_input_temp(:,:,:)/one_year * (1.-ratio_nh4_fert)

                fieldname='Nfert_ammo_pasture'
                CALL slowproc_Ninput(kjpindex, lalo, neighbours, resolution, contfrac, fieldname, &
                     N_input_temp, Ninput_year)
                ! Conversion from gN/m2(pasture)/yr to gN/m2/day
                N_input_raw(:,:,:,ifert_ammo) = N_input_raw(:,:,:,ifert_ammo)+ N_input_temp(:,:,:)/one_year

                fieldname='Nfert_nitr_pasture'
                CALL slowproc_Ninput(kjpindex, lalo, neighbours, resolution, contfrac, fieldname, &
                     N_input_temp, Ninput_year)
                ! Conversion from gN/m2(pasture)/yr to gN/m2/day
                N_input_raw(:,:,:,ifert_nitr) = N_input_raw(:,:,:,ifert_nitr)+ N_input_temp(:,:,:)/one_year

             ENDIF

             IF(.NOT. impose_ninput_manure) THEN
                N_input_raw(:,:,:,imanure) = zero
               fieldname='Nmanure_cropland'
               CALL slowproc_Ninput(kjpindex, lalo, neighbours, resolution, contfrac, fieldname, &
                    N_input_temp, Ninput_year)
               ! Conversion from gN/m2(cropland)/yr to gN/m2/day
                N_input_raw(:,:,:,imanure) = N_input_raw(:,:,:,imanure)+N_input_temp(:,:,:)/one_year

               fieldname='Nmanure_pasture'
               CALL slowproc_Ninput(kjpindex, lalo, neighbours, resolution, contfrac, fieldname, &
                    N_input_temp, Ninput_year)
               ! Conversion from gN/m2(cropland)/yr to gN/m2/day
               N_input_raw(:,:,:,imanure) = N_input_raw(:,:,:,imanure)+N_input_temp(:,:,:)/one_year
             ENDIF

             IF(.NOT. impose_ninput_bnf) THEN
                fieldname='Nbnf'
                CALL slowproc_Ninput(kjpindex, lalo, neighbours, resolution, contfrac, fieldname, &
                     N_input_raw(:,:,:,ibnf), Ninput_year)
                ! Conversion from mgN/m2/year to gN/m2/day
                N_input_raw(:,:,:,ibnf) = N_input_raw(:,:,:,ibnf)/1000./one_year
             ENDIF

          ENDIF
          
       ENDIF
    ENDIF

    
    !! 4. Main call to STOMATE
    IF ( ok_stomate ) THEN

       !! 4.1 Call stomate main routine that will call all c-cycle routines       !
       CALL stomate_main (kjit, kjpij, kjpindex, njsc,&
            IndexLand, lalo, neighbours, resolution, contfrac, frac_nobio, &
            totfrac_nobio, clayfraction, &
            siltfraction, bulk, temp_air, temp_sol, stempdiag, precip_longterm, &
            vegstress, humrel, & 
            shumdiag, litterhumdiag, precip_rain, precip_snow, &
            tmc_pft, drainage_pft, runoff_pft, swc_pft, gpp, VPD, &
            deadleaf_cover, &
            assim_param, qsintveg, &
            frac_age, veget, veget_max, &
            veget_max_new, loss_gain, frac_nobio_new, fraclut, &
            rest_id_stom, hist_id_stom, hist_id_stom_IPCC, &
            fco2_flux_out, fco2_lu_out, fco2_wh_out, fco2_ha_out, &
            resp_maint,resp_hetero,resp_growth,temp_growth, soil_pH, &
            pb, n_input_raw, month_end, &
            tdeep, hsdeep_long, snow, heat_Zimov, sfluxCH4_deep, sfluxCO2_deep, &
            som_total, snowdz, snowrho, altmax, depth_organic_soil, &
            cn_leaf_min_2D,cn_leaf_max_2D,cn_leaf_init_2D, &
            circ_class_biomass, &
            circ_class_n, lai_per_level, &
            laieff_fit, Light_Abs_Tot, Light_Tran_Tot, &
            laieff_isotrop, z_array_out, max_height_store, &
            transpir, transpir_mod, &
            coszang, stressed, unstressed, &
            u, v, mcs_hydrol, mcfc_hydrol, vessel_loss, &
            root_profile, root_depth, us, Pgap_cumul, z0m, &
            wtp,shumdiag_peat, mc_peat_above,&
            liqwt_ratio,  shumdiag_croppeat, &
            mc_croppeat_above, &
            shumdiag_man,  mc_man_above, soiltile)

       !+++CHECK+++
       ! Why is this written here? As a rule of thumb output should be written 
       ! where it is calculated. The same variables have just been written in 
       ! stomate_lpj. The calculation of npp is wrong because it does not account
       ! for atm_to_bm. Correct or even better delete this block.

       !! 4.2 Output the respiration terms and the net primary
       !!     production (NPP) that are calculated in STOMATE

       ! 4.2.1 Output the three respiration terms
       ! These variables could be output from stomate.
       ! Variables per pft. Note that growth and maintanance respiration
       ! are calculated only once per day but hetero_resp is calculated 
       ! every half hour. Because we use hetero_resp which is daily and
       ! has passed through stomate_lpj the output has also become daily
       ! and is therefore written in the do_slow loop.
       CALL xios_orchidee_send_field("maint_resp",resp_maint/dt_sechiba)
       CALL xios_orchidee_send_field("hetero_resp",resp_hetero/dt_sechiba)
       CALL xios_orchidee_send_field("growth_resp",resp_growth/dt_sechiba)

       ! Variables on grid-cell
       CALL xios_orchidee_send_field("rh_ipcc2",SUM(resp_hetero,dim=2)/dt_sechiba)
       histvar(:)=zero
       DO jv = 2, nvm
          IF ( .NOT. is_tree(jv) .AND. natural(jv) ) THEN
             histvar(:) = histvar(:) + resp_hetero(:,jv)
          ENDIF
       ENDDO
       CALL xios_orchidee_send_field("rhGrass",histvar/dt_sechiba)

       histvar(:)=zero
       DO jv = 2, nvm
          IF ( (.NOT. is_tree(jv)) .AND. (.NOT. natural(jv)) ) THEN
             histvar(:) = histvar(:) + resp_hetero(:,jv)
          ENDIF
       ENDDO
       CALL xios_orchidee_send_field("rhCrop",histvar/dt_sechiba)

       histvar(:)=zero
       DO jv = 2, nvm
          IF ( is_tree(jv) ) THEN
             histvar(:) = histvar(:) + resp_hetero(:,jv)
          ENDIF
       ENDDO
       CALL xios_orchidee_send_field("rhTree",histvar/dt_sechiba)

       ! 4.2.2 Compute the net primary production as the diff from
       ! Gross primary productin and the growth and maintenance respirations
       npp(:,1)=zero
       DO j = 2,nvm
          npp(:,j) = gpp(:,j) - resp_growth(:,j) - resp_maint(:,j)
       ENDDO

       
       ! 4.2.2 Output npp & respiration terms
       CALL xios_orchidee_send_field("npp",npp/dt_sechiba)

       ! Output with IOIPSL
       CALL histwrite_p(hist_id, 'npp', kjit, npp, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'maint_resp', kjit, resp_maint, &
            kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'hetero_resp', kjit, resp_hetero, &
            kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'growth_resp', kjit, resp_growth, &
            kjpindex*nvm, indexveg)

       ! Write the same information to a different history file file 
       IF ( hist2_id > 0 ) THEN
          CALL histwrite_p(hist2_id, 'maint_resp', kjit, resp_maint, &
               kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'hetero_resp', kjit, resp_hetero, &
               kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'growth_resp', kjit, resp_growth, &
               kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'npp', kjit, npp, kjpindex*nvm, indexveg)
       ENDIF

    ELSE

       ! ok_stomate is not activated
       ! Define the CO2 flux from the grid point to zero (no carbone cycle)
       fco2_flux_out(:) = zero
       fco2_lu_out(:) = zero
       fco2_wh_out(:) = zero
       fco2_ha_out(:) = zero
  
    ENDIF ! ok_stomate
    !+++++++++++++


    !! Calculate height and lai from biomass
    IF (read_lai .EQ. 1) THEN

       CALL slowproc_impose_site_level_lai(kjit,kjpindex,veget_max, coszang, temp_growth, max_height_store, height, height_dom, &
            lai, fco2_flux_out, fco2_lu_out, fco2_wh_out, fco2_ha_out, veget, qsintmax, circ_class_biomass, circ_class_n, &
            lai_per_level, z_array_out, laieff_fit)
       
    ELSE
       !! When the LAI is read at site level, the heights are already calculated by the bypass stomate.
              
       DO jv = 2,nvm
          DO ji = 1,kjpindex
             ! Skip if veget_max = 0, else calculate height and lai
             IF (veget_max(ji,jv) .EQ. zero) THEN
                height(ji,jv) = zero
                height_dom(ji,jv) = zero
                height_inv(ji,jv) = zero
                dia_dom(ji,jv) = zero
                height_inv(ji,jv) = zero
                ba_inv(ji,jv) = zero
                ind_inv(ji,jv) = zero
                lai(ji,jv) = zero
             ELSE
                height(ji,jv) = wood_to_qmheight(circ_class_biomass(ji,jv,:,:,icarbon), &
                     circ_class_n(ji,jv,:), jv, pipe_tune2(ji,jv))
                IF(is_tree(jv))THEN
                   ! Calculate diameter and height of the dominant diameter class
                   circ_height(:) = wood_to_height(circ_class_biomass(ji,jv,:,:,icarbon),jv, pipe_tune2(ji,jv))
                   height_dom(ji,jv) = circ_height(ncirc)
                   circ_dia(:) = wood_to_dia(circ_class_biomass(ji,jv,:,:,icarbon),jv,pipe_tune2(ji,jv))
                   dia_dom(ji,jv) = circ_dia(ncirc)
                   ! Calculate quadratic mean height and diameter of the trees
                   ! of which the diameter exceeds a threshold. This is closely
                   ! matches the height and diameter reported in forest 
                   ! inventories
                   circ_n(:) = circ_class_n(ji,jv,:)
                   WHERE (circ_dia(:) .lt. dia_thresh_inv(jv))
                      circ_dia(:) = zero
                      circ_height(:) = zero
                      circ_n(:) = zero
                   END WHERE
                   IF (SUM(circ_n(:)).EQ.zero) THEN
                      dia_inv(ji,jv) =  zero
                      height_inv(ji,jv) = zero
                      ba_inv(ji,jv) = zero
                      ind_inv(ji,jv) = zero
                   ELSE
                      dia_inv(ji,jv) = (SUM(circ_dia(:)**2*circ_n(:))/SUM(circ_n(:)))**0.5
                      height_inv(ji,jv) = (SUM(circ_height(:)**2*circ_n(:))/SUM(circ_n(:)))**0.5
                      ind_inv(ji,jv) = SUM(circ_n(:))
                      ba_inv(ji,jv) = pi/4*SUM(circ_dia(:)**2*circ_n(:))*m2_to_ha
                   END IF
                ELSE
                   ! Grasses and crops have a height but no
                   ! diameter in ORCHIDEE
                   height_dom(ji,jv) = height(ji,jv)
                   height_inv(ji,jv) = height(ji,jv)
                   dia_dom(ji,jv) = undef
                   dia_inv(ji,jv) = undef
                   ba_inv(ji,jv) = undef
                   ind_inv(ji,jv) = circ_class_n(ji,jv,1)
                ENDIF
                lai(ji,jv) = cc_to_lai(circ_class_biomass(ji,jv,:,ileaf,icarbon), &
                     circ_class_n(ji,jv,:),jv)
             ENDIF
          ENDDO
       ENDDO
    END IF ! read_lai .EQ. 1
    !-

    !! 5. Do daily processes if necessary
    !     This section is not done for read_lai=1
    IF (read_lai .NE. 1) THEN
       IF ( do_slow ) THEN
          
          !!  5.1 Calculate canopy structure when STOMATE is not activated
          IF ( .NOT. ok_stomate ) THEN
             
             ! Slowproc_canopy calculates the canopy structure including
             ! laieff_fit. It should therefore only be calculated once 
             ! per day
             CALL slowproc_canopy (kjpindex, circ_class_biomass, circ_class_n, &
                  veget_max, lai_per_level, z_array_out, &
                  max_height_store, laieff_fit, frac_age)
             
          ENDIF
          
          CALL slowproc_veget (kjpindex, lai_per_level, z_array_out, &
               coszang_noon, circ_class_biomass, circ_class_n, frac_nobio, &
               totfrac_nobio, veget_max, veget, soiltile, tot_bare_soil, &
               fraclut, nwdFraclut, Pgap_cumul,precip_longterm)
          
          !! 5.3 updates qsintmax and other derived variables
          IF ( .NOT. ok_stomate ) THEN
             
             ! Initialize missing variables
             deadleaf_cover(:) = zero
             temp_growth(:) = 25.
             
          ENDIF
          qsintmax(:,:) = qsintcst * veget(:,:) * lai(:,:)
          qsintmax(:,1) = zero
         
          ! Do some basic tests on the surface fractions updated above, only if
          ! slowproc_veget has been done (do_slow). No change of the variables. 
          CALL check_veget(kjpindex, frac_nobio, veget_max, veget, &
               tot_bare_soil, soiltile, failed_vegfrac)
          
       END IF

       !! 8. Write output fields
       CALL xios_orchidee_send_field("tot_bare_soil",tot_bare_soil)
       CALL xios_orchidee_send_field("HEIGHT_DOM",height_dom)
       CALL xios_orchidee_send_field("HEIGHT_INV",height_inv)
       CALL xios_orchidee_send_field("DIA_DOM",dia_dom)
       CALL xios_orchidee_send_field("DIA_INV",dia_inv)
       CALL xios_orchidee_send_field("BA_INV",ba_inv)
       CALL xios_orchidee_send_field("IND_INV",ind_inv)
       
       IF ( .NOT. almaoutput) THEN
          CALL histwrite_p(hist_id, 'tot_bare_soil', kjit, tot_bare_soil, &
               kjpindex, IndexLand)
       END IF
    
       ! Error checking
       IF(err_act.GT.1)THEN
          
          ! All initial checks should be done in slowproc right after the map
          ! is being read. If vegetation fractions or frac_nobio is adjusted
          ! afterwards, mass balance problems are unavoidable. Check whether
          ! veget_max and frac_nobio are still consistent.
          CALL check_pixel_area("End of slowproc_main", kjpindex, veget_max, frac_nobio)
          
       END IF ! err_act.GT.1
    END IF! read_lai .NE. 1


    !! 6.1. Call to interpolation of dynamic irrigation map, if time to do so (as for vegetmax interpolation)
    !! Important difference: veget map is updated every veget_update years. Irrig map_pft
    !! is updated every year for now
    !! Put here to use updates values of soiltile and veget_max
    irrigated_next(:) = zero
    IF (do_irrigation .AND. irrig_map_dynamic_flag) THEN
      ! Attention: veget_year already updated, but veget_update must be >0, I.E. it must read veget maps
      ! Seems logic to update irrigation maps if vegetation maps are updated too
      IF ( FirstTsYear) THEN
        CALL slowproc_readirrigmap_dyn(kjpindex, lalo, neighbours,  resolution, contfrac,         &
             irrigated_new)
        DO ji=1,kjpindex
           IF( area(ji)*contfrac(ji) > min_sechiba) THEN !Multiplication is total area(ig) = m2

             irrig_frac_next(ji) = MIN( soiltile(ji,irrig_st) * SUM(veget_max(ji,:)) , &
                 irrigated_new(ji) / ( area(ji)*contfrac(ji) ) )
             ! soiltile(fraction of vegtot) * SUM(veget_max)  = fraction of grid cell
             ! irrigated(m2)/ grid_cell_area = Fraction of grid cell
             ! irrig_frac is always fraction of grid cell
           ENDIF
        ENDDO
      ENDIF
     !! Here irrigated_next from sechiba = irrigated_new from slowproc!!
     !! Dynamic irrigation maps as output to sechiba_end
    ENDIF
    irrigated_next(:) = irrigated_new(:)

    !!
    !! 6.2 Calculation of reinf_slope_soil to pass to hydrol
    IF (FirstTsYear) THEN
      reinf_slope_soil(:,:) = zero
      DO ji=1,kjpindex
          reinf_slope_soil(ji,:) = reinf_slope(ji)
          IF( Reinfiltr_IrrigField .AND.  irrig_frac_next(ji) > min_sechiba ) THEN
            reinf_slope_soil(ji, irrig_st) = MAX(reinf_slope(ji), reinf_slope_cropParam)
          ENDIF
      ENDDO
    ENDIF

    IF (printlev_loc>=3)  WRITE (numout,*) ' slowproc_main done '

  END SUBROUTINE slowproc_main


!! ================================================================================================================================
!! SUBROUTINE 	: slowproc_finalize
!!
!>\BRIEF         Write to restart file variables for slowproc module and call finalization of stomate module
!!
!! DESCRIPTION : 
!!
!! RECENT CHANGE(S): Add arrays of soil hydraulic parameters to the restart file.
!!                   Linked to SP-MIP project (Tafasca Salma and Ducharne Agnes).
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_finalize (kjit,       kjpindex,  rest_id,  IndexLand,  &
                                njsc,       veget,                           &
                                frac_nobio, veget_max, reinf_slope,          &
                                ks,  nvan, avan, mcr, mcs, mcfc, mcw,        &
                                assim_param, frac_age, heat_Zimov, altmax, depth_organic_soil, &
                                circ_class_biomass, circ_class_n,            &
                                lai_per_level, laieff_fit, loss_gain,        &
                                veget_max_new, frac_nobio_new, fraction_aeirrig_sw, &
                                precip_longterm)

!! 0.1 Input variables
    INTEGER(i_std),INTENT(in)                            :: kjit           !! Time step number
    INTEGER(i_std),INTENT(in)                            :: kjpindex       !! Domain size - terrestrial pixels only
    INTEGER(i_std),INTENT(in)                            :: rest_id        !! Restart file identifier
    INTEGER(i_std),DIMENSION(kjpindex),INTENT(in)        :: IndexLand      !! Indices of the points on the land map
    INTEGER(i_std),DIMENSION(kjpindex),INTENT(in)        :: njsc           !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless)
    REAL(r_std),DIMENSION(kjpindex,nvm),INTENT(in)       :: veget          !! Fraction of vegetation type in the mesh (unitless)
    REAL(r_std),DIMENSION(kjpindex,nnobio),INTENT(in)    :: frac_nobio     !! Fraction of ice, lakes, cities etc. in the mesh (unitless)
    REAL(r_std),DIMENSION(kjpindex,nvm),INTENT(in)       :: veget_max      !! Maximum fraction of vegetation type in the mesh (unitless)
    REAL(r_std),DIMENSION(kjpindex),INTENT(in)           :: reinf_slope    !! slope coef for reinfiltration
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)         :: ks             !! Hydraulic conductivity at saturation (mm {-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)         :: nvan           !! Van Genuchten coeficients n (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)         :: avan           !! Van Genuchten coeficients a (mm-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)         :: mcr            !! Residual volumetric water content (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)         :: mcs            !! Saturated volumetric water content (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)         :: mcfc           !! Volumetric water content at field capacity (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)         :: mcw            !! Volumetric water content at wilting point (m^{3} m^{-3})
    REAL(r_std),DIMENSION(kjpindex,nvm,npco2),INTENT(in) :: assim_param    !! assimilation parameters vcmax, nue, and leaf nitrogen
    REAL(r_std),DIMENSION(kjpindex,nvm,nleafages),INTENT(in) :: frac_age   !! Age efficacity from STOMATE for isoprene
    REAL(r_std),DIMENSION(:,:,:,:,:),INTENT(in)          :: circ_class_biomass !!  Biomass components of the model tree  
                                                                           !! within a circumference class
                                                                           !! class @tex $(g C ind^{-1})$ @endtex
    REAL(r_std),DIMENSION(:,:,:),INTENT(in)              :: circ_class_n   !! Number of trees within each circumference
                                                                           !! class @tex $(m^{-2})$ @endtex
    REAL(r_std),DIMENSION(:,:,:),INTENT(in)              :: lai_per_level  !! This is the LAI per vertical level
                                                                           !! @tex $(m^{2} m^{-2})$
    TYPE(laieff_type),DIMENSION(:,:,:),INTENT(in) &
                                                         :: laieff_fit     !! Fitted parameters for the effective LAI
    REAL(r_std), DIMENSION(:,:), INTENT(in)              :: loss_gain      !! Changes in veget_max distributed over all
                                                                           !! age classes and thus taking the age-classes into 
                                                                           !! account (unitless, 0-1)
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT(in)  :: frac_nobio_new !! Fraction of ice,lakes,cities, ... (unitless)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(in)     :: veget_max_new  !! Maximum fraction of vegetation type including none 
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)         :: precip_longterm!! Longterm annual precipitation sum mm year^{-1}


!! 0.3 Modified variables
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm), INTENT(in):: heat_Zimov    !! heating associated with decomposition [W/m**3 soil]
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)    :: altmax         !! Maximul active layer thickness (m). Be careful, here active means non frozen.
                                                                           !! Not related with the active soil carbon pool.
    REAL(r_std), DIMENSION(kjpindex), INTENT (in)        :: depth_organic_soil !! how deep is the organic soil?
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)         :: fraction_aeirrig_sw    !! Fraction of area equipped for irrigation from surface water, of irrig_frac
                                                                                   !! 1.0 here corresponds to fraction of irrigated area, not grid cell
!! 0.4 Local variables
    REAL(r_std),DIMENSION(kjpindex,nmonth)               :: Ninput_ammo    !! Daily ammonium inputs (gN m-2 day-1) 
    REAL(r_std),DIMENSION(kjpindex,nmonth)               :: Ninput_nitr    !! Daily nitrate inputs (gN m-2 day-1)
    REAL(r_std),DIMENSION(kjpindex,nmonth)               :: Ninput_fert    !! Daily N fertilization (gN m-2 day-1)
    REAL(r_std),DIMENSION(kjpindex,nmonth)               :: Ninput_bnf     !! Daily biological N fixation (gN m-2 day-1)
    REAL(r_std)                                          :: tmp_day(1)     !! temporary variable for I/O
    INTEGER                                              :: ji, jv, imonth !! Indices
    INTEGER                                              :: iele,im,jf,ivm !! Indices
    CHARACTER(LEN=4)                                     :: laistring      !! Temporary character string
    CHARACTER(LEN=80)                                    :: var_name       !! To store variables names for I/O
    CHARACTER(LEN=2), DIMENSION(nelements)               :: element_str    !! Element string used to make nice variables names 
                                                                           !! in restart files
    CHARACTER(LEN=10)                                    :: part_str       !! string suffix indicating an index
    CHARACTER(LEN=2)                                     :: month_str      !! string used in variable name in restart files
    
!_ ================================================================================================================================

    IF (printlev_loc>=3) WRITE (numout,*) 'Write restart file with SLOWPROC variables '

    !! 1. Define element string
    DO iele = 1,nelements
       IF (iele == icarbon) THEN 
          element_str(iele) = 'c' 
       ELSEIF (iele == initrogen) THEN 
          element_str(iele) = 'n' 
       ELSE 
          STOP 'Define element_str' 
       ENDIF
    ENDDO

    ! 2.1 Write a series of variables controled by slowproc to the restart file
    CALL restput_p (rest_id, 'veget', nbp_glo, nvm, 1, kjit, &
         veget, 'scatter',  nbp_glo, index_g)
    !-
    CALL restput_p (rest_id, 'veget_max', nbp_glo, nvm, 1, kjit, &
         veget_max, 'scatter',  nbp_glo, index_g)
    !-
    CALL restput_p (rest_id, 'circ_class_n', nbp_glo, nvm, ncirc, kjit, &
         circ_class_n, 'scatter', nbp_glo, index_g)
    !-
    DO iele = 1,nelements
       var_name = 'cc_biomass_'//TRIM(element_str(iele))
       CALL restput_p (rest_id, var_name, nbp_glo, nvm, ncirc, nparts, kjit, &
            circ_class_biomass(:,:,:,:,iele), 'scatter', nbp_glo, index_g)
       !-
    ENDDO
    CALL restput_p (rest_id, 'assim_param',nbp_glo, nvm, npco2, kjit, &
         assim_param, 'scatter', nbp_glo, index_g)
    !-
    CALL restput_p (rest_id, 'frac_nobio', nbp_glo, nnobio, 1, kjit, frac_nobio, 'scatter',  nbp_glo, index_g)
    !-
    CALL restput_p (rest_id, 'loss_gain', nbp_glo, nvm, 1, kjit, &
         loss_gain, 'scatter',  nbp_glo, index_g)
    !-
    CALL restput_p (rest_id, 'veget_max_new', nbp_glo, nvm, 1, kjit, &
         veget_max_new, 'scatter',  nbp_glo, index_g)
    !-
    CALL restput_p (rest_id, 'frac_nobio_new', nbp_glo, nnobio, 1, kjit, frac_nobio_new, 'scatter',  nbp_glo, index_g)
    !-
    CALL restput_p (rest_id, 'precip_longterm', nbp_glo, 1, 1, kjit, precip_longterm, 'scatter',  nbp_glo, index_g)

    !+++CHECK+++
!!$    ! Are the values of frac_age what we expect for ok_bvoc? Seems that not 
!!$    ! much is happening with frac_age for the moment
!!$    CALL restput_p (rest_id, 'frac_age', nbp_glo, nvm, nleafages, kjit, &
!!$         frac_age, 'scatter',  nbp_glo, index_g)
    !+++++++++++

    IF ( do_irrigation ) THEN
          CALL restput_p (rest_id, 'irrigmap_dyn', nbp_glo, 1, 1, kjit, irrigated_new, 'scatter',  nbp_glo, index_g)
          IF ( select_source_irrig ) THEN
                CALL restput_p (rest_id, 'fraction_aeirrig_sw', nbp_glo, 1, 1, kjit, fraction_aeirrig_sw, 'scatter',  nbp_glo, index_g)
          ENDIF
    ENDIF

    ! Add the soil_classif as suffix for the variable name of njsc when it is stored in the restart file. 
    IF (soil_classif == 'zobler') THEN
       var_name= 'njsc_zobler'
    ELSE IF (soil_classif == 'usda') THEN
       var_name= 'njsc_usda'
    END IF

    CALL restput_p (rest_id, var_name, nbp_glo, 1, 1, kjit, REAL(njsc, r_std), 'scatter',  nbp_glo, index_g)
    !-
    CALL restput_p (rest_id, 'reinf_slope', nbp_glo, 1, 1, kjit, reinf_slope, 'scatter',  nbp_glo, index_g)
    !-
    CALL restput_p (rest_id, 'clay_frac', nbp_glo, 1, 1, kjit, clayfraction, 'scatter',  nbp_glo, index_g)
    !-
    CALL restput_p (rest_id, 'sand_frac', nbp_glo, 1, 1, kjit, sandfraction, 'scatter',  nbp_glo, index_g)
    !-
    CALL restput_p (rest_id, 'silt_frac', nbp_glo, 1, 1, kjit, siltfraction, 'scatter',  nbp_glo, index_g)
    !-
    CALL restput_p (rest_id, 'bulk', nbp_glo, 1, 1, kjit, bulk, 'scatter',  nbp_glo, index_g)
    !-
    CALL restput_p (rest_id, 'soil_ph', nbp_glo, 1, 1, kjit, soil_ph, 'scatter',  nbp_glo, index_g)
    !-
    CALL restput_p (rest_id, 'mcs', nbp_glo, 1, 1, kjit, mcs, 'scatter', nbp_glo, index_g)
    !-
    CALL restput_p (rest_id, 'mcr', nbp_glo, 1, 1, kjit, mcr, 'scatter', nbp_glo, index_g)
    !-
    CALL restput_p (rest_id, 'mcfc', nbp_glo, 1, 1, kjit, mcfc, 'scatter', nbp_glo, index_g)
    !-
    CALL restput_p (rest_id, 'mcw', nbp_glo, 1, 1, kjit, mcw, 'scatter', nbp_glo, index_g)
    !-
    CALL restput_p (rest_id, 'avan', nbp_glo, 1, 1, kjit, avan, 'scatter', nbp_glo, index_g)
    !-
    CALL restput_p (rest_id, 'nvan', nbp_glo, 1, 1, kjit, nvan, 'scatter', nbp_glo, index_g)
    !-
    CALL restput_p (rest_id, 'ks', nbp_glo, 1, 1, kjit, ks, 'scatter', nbp_glo, index_g)
    !-
    
    ! Specific case where the biomass is prescribed or read and not calculated by STOMATE
    ! Note that the cc_biomass_m variable has different dimensions than the standard cc_biomass.
    ! The same applies for cc_n_m and cc_n.
    IF ((read_lai .EQ. 2) .OR. (read_lai .EQ. 0)) THEN
       DO iele = 1,nelements
          DO imonth = 1,nmonth
10           FORMAT(I2)
             WRITE (month_str,10) imonth
             var_name = 'cc_biomass_m_'//TRIM(month_str)//'_'//TRIM(element_str(iele))
             CALL restput_p (rest_id, var_name, nbp_glo, nvm, ncirc, nparts, kjit, &
                  cc_biomass_m(:,:,:,:,imonth,iele), 'scatter', nbp_glo, index_g)
          ENDDO
       ENDDO
       var_name = 'cc_n_m'
       CALL restput_p (rest_id, var_name, nbp_glo, nvm, ncirc, 12, kjit, &
            cc_n_m(:,:,:,:), 'scatter', nbp_glo, index_g)
    ENDIF


    IF (ok_ncycle .AND. (.NOT. impose_CN)) THEN
       CALL restput_p (rest_id, 'Nammonium', nbp_glo, nvm , 12, kjit, N_input_raw(:,:,:,iammonium), 'scatter',  nbp_glo, index_g)
       !-
       CALL restput_p (rest_id, 'Nnitrate', nbp_glo, nvm, 12, kjit, N_input_raw(:,:,:,initrate), 'scatter',  nbp_glo, index_g)
       !-
       CALL restput_p (rest_id, 'Nfert_ammo', nbp_glo, nvm, 12, kjit, N_input_raw(:,:,:,ifert_ammo), 'scatter',  nbp_glo, index_g)
       !-
       CALL restput_p (rest_id, 'Nfert_nitr', nbp_glo, nvm, 12, kjit, N_input_raw(:,:,:,ifert_nitr), 'scatter',  nbp_glo, index_g)
       !-
       CALL restput_p (rest_id, 'Nmanure', nbp_glo, nvm, 12, kjit, N_input_raw(:,:,:,imanure), 'scatter',  nbp_glo, index_g)
       !-
       CALL restput_p (rest_id, 'Nbnf', nbp_glo, nvm, 12, kjit, N_input_raw(:,:,:,ibnf), 'scatter',  nbp_glo, index_g)
       !-
    END IF

    !
    ! If there is some N inputs change, write the year 
    CALL restput_p (rest_id, 'Ninput_year', kjit, ninput_year)
    
    ! 2.2 Write restart variables managed by STOMATE
    ! The restart files of stomate and sechiba both contain circ_class_biomass, and circ_class_n. 
    ! This allows to run a first simulation with stomate and then use its restarts
    ! for a simulation with only sechiba.
    IF ( ok_stomate ) THEN
       CALL stomate_finalize (kjit, kjpindex, indexLand, clayfraction, siltfraction, &
            bulk, assim_param, &
            heat_Zimov, altmax, depth_organic_soil, circ_class_biomass, circ_class_n, &
            lai_per_level, laieff_fit, veget_max)
    ENDIF

  END SUBROUTINE slowproc_finalize


!! ================================================================================================================================
!! SUBROUTINE   : slowproc_init
!!
!>\BRIEF         Initialisation of all variables linked to SLOWPROC
!!
!! DESCRIPTION  : (definitions, functional, design, flags): The subroutine manages 
!! diverses tasks:
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::
!! ::veget, ::frac_nobio, ::totfrac_nobio, ::veget_max, ::height, ::soiltype
!! ::Ninput_update
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_init (kjit, kjpindex, IndexLand, lalo, neighbours, resolution, contfrac, &
       rest_id, frac_age, veget, frac_nobio, totfrac_nobio, soiltile, fraclut, nwdfraclut, reinf_slope, &
       ks,  nvan, avan, mcr, mcs, mcfc, mcw, & 
       veget_max, tot_bare_soil, njsc, &
       Ninput_update, Ninput_year, &
       circ_class_biomass, circ_class_n, assim_param, loss_gain, veget_max_new, &
       frac_nobio_new, fraction_aeirrig_sw, precip_longterm)
    
    !! INTERFACE DESCRIPTION

    !! 0.1 Input variables
    INTEGER(i_std), INTENT (in)                           :: kjit           !! Time step number
    INTEGER(i_std), INTENT (in)                           :: kjpindex       !! Domain size - Terrestrial pixels only 
    INTEGER(i_std), INTENT (in)                           :: rest_id        !! Restart file identifier
    INTEGER(i_std), DIMENSION (kjpindex), INTENT (in)     :: IndexLand      !! Indices of the land points on the map
    REAL(r_std), DIMENSION (kjpindex,2), INTENT (in)      :: lalo           !! Geogr. coordinates (latitude,longitude) (degrees)
    INTEGER(i_std), DIMENSION (kjpindex,NbNeighb), INTENT(in):: neighbours  !! Vector of neighbours for each grid point
                                                                            !! (1=North and then clockwise)
    REAL(r_std), DIMENSION (kjpindex,2), INTENT(in)       :: resolution     !! size in x and y of the grid (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)         :: contfrac       !! Fraction of continent in the grid (unitless)
    
    !! 0.2 Output variables
    INTEGER(i_std), INTENT(out)                           :: Ninput_update       !! update frequency in timesteps (years) for N inputs
    INTEGER(i_std), INTENT(out)                           :: Ninput_year         !! Year for the nitrogen inputs
    INTEGER(i_std), DIMENSION(kjpindex), INTENT(out)      :: njsc                !! Index of the dominant soil textural class in the grid 
                                                                                 !! cell (1-nscm, unitless)
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT (out)   :: veget               !! Fraction of vegetation type in the mesh (unitless)
    REAL(r_std), DIMENSION (kjpindex,nnobio), INTENT (out):: frac_nobio          !! Fraction of ice,lakes,cities, in the mesh (unitless)
    REAL(r_std), DIMENSION (kjpindex), INTENT (out)       :: totfrac_nobio       !! Total fraction of ice+lakes+cities+ in the mesh (unitless) 
    REAL(r_std), DIMENSION (kjpindex), INTENT (out)       :: tot_bare_soil       !! Total evaporating bare soil fraction in the mesh (unitless)                                         
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT (out)   :: veget_max           !! Max fraction of vegetation type in the mesh (unitless)
    REAL(r_std), DIMENSION (kjpindex,nvm,nleafages), INTENT (out):: frac_age     !! Age efficacity from STOMATE for isoprene
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT(out)   :: soiltile            !! Fraction of each soil tile within vegtot 
                                                                                 !! (0-1, unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)          :: reinf_slope    !! slope coef for reinfiltration
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: ks             !! Hydraulic conductivity at saturation (mm {-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: nvan           !! Van Genuchten coeficients n (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: avan           !! Van Genuchten coeficients a (mm-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: mcr            !! Residual volumetric water content (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: mcs            !! Saturated volumetric water content (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: mcfc           !! Volumetric water content at field capacity (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: mcw            !! Volumetric water content at wilting point (m^{3} m^{-3})
    REAL(r_std),DIMENSION (:,:,:,:,:), INTENT(out)         :: circ_class_biomass  !! Biomass components of the model tree  
                                                                                 !! within a circumference class
                                                                                 !! class @tex $(g C ind^{-1})$ @endtex
    REAL(r_std), DIMENSION (:,:,:), INTENT(out)           :: circ_class_n        !! Number of trees within each circumference
                                                                                 !! class @tex $(m^{-2})$ @endtex
    REAL(r_std),DIMENSION (:,:,:), INTENT(out)            :: assim_param         !! assimilation parameters, vcmax, nue and leaf nitrogen
    REAL(r_std), DIMENSION (kjpindex,nlut), INTENT(out)   :: fraclut             !! Fraction of each landuse tile
    REAL(r_std), DIMENSION (kjpindex,nlut), INTENT(out)   :: nwdfraclut          !! Fraction of non woody vegetation in each landuse tile
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(out)    :: loss_gain           !! Fractional losses and gains following a land cover change (unitless, 0-1) 
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT(out)  :: frac_nobio_new      !! Fraction of ice,lakes,cities, ... (unitless)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(out)     :: veget_max_new       !! Maximum fraction of vegetation type including none 
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)        :: fraction_aeirrig_sw !! Fraction of area equipped for irrigation from surface water, of irrig_frac
                                                                                 !! 1.0 here corresponds to fraction of irrig. area, not grid cell
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)        :: precip_longterm     !! Fraction of area equipped for irrigation from surface water, of irrig_frac


    !! 0.3 Modified variables

    !! 0.4 Local variables 
    REAL(r_std)                                           :: zcanop              !! ???? soil depth taken for canopy
    INTEGER(i_std)                                        :: vtmp(1)             !! temporary variable
    REAL(r_std), DIMENSION(nslm)                          :: zsoil               !! soil depths at diagnostic levels
    INTEGER(i_std)                                        :: njsc_imp            !! njsc to impose nvan, ks etc. if impsoil
    CHARACTER(LEN=4)                                      :: laistring           !! Temporary character string
    INTEGER(i_std)                                        :: l, jf, im           !! Indices
    CHARACTER(LEN=80)                                     :: var_name            !! To store variables names for I/O
    INTEGER(i_std)                                        :: ji, jv, ier,jst     !! Indices 
    LOGICAL                                               :: get_slope
    REAL(r_std)                                           :: frac_nobio1         !! temporary variable for frac_nobio(see above)
    REAL(r_std)                                           :: precip_longterm_init !! Pluie annuelle supposee au cold start (mm/an)
    REAL(r_std), DIMENSION(kjpindex)                      :: tmp_real
    REAL(r_std), DIMENSION(kjpindex,nslm)                 :: stempdiag2_bid      !! matrix to store stempdiag_bid
    REAL(r_std), DIMENSION (kjpindex,nscm)                :: soilclass           !! Fractions of each soil textural class in the grid cell (0-1, unitless)
    CHARACTER(LEN=30), SAVE                               :: ninput_str          !! update frequency for N inputs
!$OMP THREADPRIVATE(ninput_str)
    CHARACTER(LEN=10)                                     :: part_str            !! string suffix indicating an index
    REAL(r_std), DIMENSION(kjpindex)                      :: frac_crop_tot       !! Total fraction occupied by crops (0-1, unitless)
    REAL(r_std), DIMENSION(kjpindex)                      :: mvan, psi_fc, psi_w  !! To calculate default wilting point and 
                                                                                  !! field capacity if impsoilt
    REAL(r_std), DIMENSION(kjpindex)                      :: mcfc_default      !! Default field capacity if impsoilt
    REAL(r_std), DIMENSION(kjpindex)                      :: mcw_default       !! Default wilting point if impsoilt
    REAL(r_std)                                           :: nvan_default      !! Default  if impsoilt
    REAL(r_std)                                           :: avan_default      !! Default  if impsoilt
    REAL(r_std)                                           :: mcr_default       !! Default  if impsoilt
    REAL(r_std)                                           :: mcs_default       !! Default  if impsoilt
    REAL(r_std)                                           :: ks_default        !! Default  if impsoilt
    REAL(r_std)                                           :: clayfraction_default  !! Default  if impsoilt
    REAL(r_std)                                           :: sandfraction_default  !! Default  if impsoilt
    LOGICAL                                               :: found_restart       !! found_restart=true if all 3 variables veget_max, 
                                                                                 !! veget and frac_nobio are read from restart file
    LOGICAL                                               :: call_slowproc_soilt !! This variables will be true if subroutine slowproc_soilt needs to be called
    CHARACTER(LEN=80)                                     :: fieldname           !! name of the field read in the N input map
    REAL(r_std)                                           :: nammonium, nnitrate !! Precribed amounts of deposition 
    REAL(r_std)                                           :: nfert, nbnf,nmanure !! Prescribed amounts of N input from fertilizer, 
                                                                                 !! biological fixation and manure
    REAL(r_std), DIMENSION(kjpindex,nvm,12)               :: N_input_temp
    INTEGER                                               :: iele, imonth, ipts  !! Indices
    INTEGER                                               :: ivmax, ivm          !! Indices
    CHARACTER(LEN=2), DIMENSION(nelements)                :: element_str         !! Element string used to make nice variables names 
                                                                                 !! in restart files
    CHARACTER(LEN=2)                                      :: month_str           !! string used in variable name in restart files
    REAL(r_std), DIMENSION(kjpindex)                      :: fracsum             !! Sum of both fracnobio and veget_max
    REAL(r_std), DIMENSION(kjpindex)                      :: count               !! pixels with an error
    REAL(r_std)                                           :: residual            !! precision error that should be corrected for
    REAL(r_std), DIMENSION(nvm)                           :: sechiba_vegmax      !! temporary variable to use impose_veg
      
!_ ================================================================================================================================


    IF (printlev_loc>=3) WRITE (numout,*) "In slowproc_init"
    
    !! 1. Allocate memory
    ALLOCATE (clayfraction(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init',&
         'Problem in allocation of variable clayfraction','','')
    clayfraction(:)=undef_sechiba
    
    ALLOCATE (sandfraction(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable sandfraction','','')
    sandfraction(:)=undef_sechiba
    
    ALLOCATE (siltfraction(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable siltfraction','','')
    siltfraction(:)=undef_sechiba

    ALLOCATE (bulk(kjpindex),stat=ier)
    IF (ier.NE.0) CALL ipslerr_p(3,'slowproc_init',&
         'Problem in allocation of variable bulk','','')
    bulk(:)=undef_sechiba

    ALLOCATE (soil_ph(kjpindex),stat=ier)
    IF (ier.NE.0) CALL ipslerr_p(3,'slowproc_init',&
         'Problem in allocation of variable soilph','','')
    soil_ph(:)=undef_sechiba

    ALLOCATE (n_input_raw(kjpindex,nvm,12,ninput),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable n_input_raw','','')

    ! Allocation of last irrig map area in case of change
    ALLOCATE(irrigated_new(kjpindex), STAT=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable irrigated_new','','')

    
    ! Allocate variables to prescribe canopy structure when
    ! stomate is not used
    IF ((read_lai .EQ. 2) .OR. (read_lai .EQ. 0)) THEN
       ! cc_biomass and cc_n are prescribed from files with monthly values
       ALLOCATE (cc_biomass_m(kjpindex,nvm,ncirc,nparts,12,nelements),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init',&
            'Problem in allocation of variable cc_biomass_m','','')
       ALLOCATE (cc_n_m(kjpindex,nvm,ncirc,12), stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init',&
            'Problem in allocation of variable cc_n_m','','')
    ELSE
       ! allocate the variables but they will never be used
       ALLOCATE (cc_biomass_m(1,1,1,1,1,1), stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init',&
            'Problem in allocation of variable cc_biomass_m(1,1,1,1,1,1)','','')
       ALLOCATE (cc_n_m(1,1,1,1), stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init',&
            'Problem in allocation of variable cc_n_m(1,1,1,1)','','')
    ENDIF

    !! 2. Read general parameters
   
    !
    ! Time step of STOMATE and LAI update when reading an LAI map
    !
    !Config Key   = DT_STOMATE
    !Config Desc  = Time step of STOMATE and other slow processes
    !Config If    = OK_STOMATE
    !Config Def   = one_day
    !Config Help  = Time step (s) of regular update of vegetation
    !Config         cover, LAI etc. This is also the time step
    !Config         of STOMATE.
    !Config Units = [seconds]
    dt_stomate = one_day
    CALL getin_p('DT_STOMATE', dt_stomate)

    !! 3. Read the soil related variables
    !  The model could start from scratch, from a restart file or the value of 
    !  specific variables could be imposed. Imposing could make use of fixed
    !  values or from maps. This results in many possible configurations some
    !  of which are not very useful. The order of the code in this subroutine 
    !  already limits the number of possible configurations to initialize the
    !  model but the real quality control on this issue is taking place in 
    !  sechiba in the subroutine check_configuration. The order implemented 
    !  here is: (1) Read the values from a restart file if available. If
    !  no restart file is found, give the variable the value val_exp (done
    !  in the subroutine restget_p). If a restart file was found, the value 
    !  for impose_soilt will be ingnored (this is taken care of in setvar_p by 
    !  checking whether the variable has the value val_exp or not). If no 
    !  restart file was found there are still two more options to initialize 
    !  the model: (2) initialize all soil variables by fixed values which can 
    !  be set in the run.def or their default values, or (3) initialize the 
    !  soil variables with values from a map or file (in this case possible 
    !  conflicts with the restart file are explicitly dealt with througf 
    !  IF-statements in this subroutine). Note that all the information 
    !  required to initialize the soil is stored in the sechiba restart. Hence, 
    !  the value of ok_stomate does not affect this part of the initialization. 

    ! The variables are read from restart file. If at least one of the related variables are not found, the  
    ! variable call_slowproc_soilt will be true and the variables will be read and interpolated from file in  
    ! subroutine slowproc_soilt. 
    call_slowproc_soilt=.FALSE. 
    
    ! Add the soil classification as suffix for the variable name of njsc to make
    ! sure that the correct njsc is read from the restart. A restart made with one
    ! soil classification cannot be used for a simulation with another soil 
    ! classification. The model will crash by saying that the specific variable 
    ! for njsc was not found in the restart file.
    IF (soil_classif == 'zobler') THEN
       var_name= 'njsc_zobler'
    ELSE IF (soil_classif == 'usda') THEN
       var_name= 'njsc_usda'
    ELSE
       CALL ipslerr_p(3,'slowproc_init',&
            'Non supported soil typeclassification','','')
    END IF

    ! Index of the dominant soil type in the grid cell
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Index of soil type')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, &
         .TRUE., tmp_real, "gather", nbp_glo, index_g)
    IF ( ALL( tmp_real(:) .EQ. val_exp) ) THEN
       njsc (:) = undef_int
       call_slowproc_soilt=.TRUE.
    ELSE
       njsc = NINT(tmp_real)
    END IF

    var_name= 'ks'
    CALL ioconf_setatt_p('UNITS', 'mm/d')
    CALL ioconf_setatt_p('LONG_NAME','Soil saturated water content')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., ks, "gather", nbp_glo, index_g)
    IF ( ALL(ks(:) .EQ. val_exp ) ) THEN 
       ! ks is not in restart file 
       call_slowproc_soilt=.TRUE. 
    END IF

    var_name= 'mcs'
    CALL ioconf_setatt_p('UNITS', 'm3/m3')
    CALL ioconf_setatt_p('LONG_NAME','')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., mcs, "gather", nbp_glo, index_g)
    IF ( ALL(mcs(:) .EQ. val_exp ) ) THEN 
       ! mcs is not in restart file 
       call_slowproc_soilt=.TRUE. 
    END IF

    var_name= 'mcr'
    CALL ioconf_setatt_p('UNITS', 'm3/m3')
    CALL ioconf_setatt_p('LONG_NAME','')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., mcr, "gather", nbp_glo, index_g)
    IF ( ALL(mcr(:) .EQ. val_exp ) ) THEN 
       ! mcr is not in restart file 
       call_slowproc_soilt=.TRUE. 
    END IF

    var_name= 'mcfc'
    CALL ioconf_setatt_p('UNITS', 'm3/m3')
    CALL ioconf_setatt_p('LONG_NAME','')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., mcfc, "gather", nbp_glo, index_g)
    IF ( ALL(mcfc(:) .EQ. val_exp ) ) THEN 
       ! mcfc is not in restart file 
       call_slowproc_soilt=.TRUE. 
    END IF

    var_name= 'mcw'
    CALL ioconf_setatt_p('UNITS', 'm3/m3')
    CALL ioconf_setatt_p('LONG_NAME','')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., mcw, "gather", nbp_glo, index_g)
    IF ( ALL(mcw(:) .EQ. val_exp ) ) THEN 
       ! mcw is not in restart file 
       call_slowproc_soilt=.TRUE. 
    END IF

    var_name= 'nvan'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., nvan, "gather", nbp_glo, index_g)
    IF ( ALL(nvan(:) .EQ. val_exp ) ) THEN 
       ! nvan is not in restart file 
       call_slowproc_soilt=.TRUE. 
    END IF

    var_name= 'avan'
    CALL ioconf_setatt_p('UNITS', 'm-1')
    CALL ioconf_setatt_p('LONG_NAME','')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., avan, "gather", nbp_glo, index_g)
    IF ( ALL(avan(:) .EQ. val_exp ) ) THEN 
       ! avan is not in restart file 
       call_slowproc_soilt=.TRUE. 
    END IF

    var_name= 'clay_frac'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Fraction of clay in each mesh')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, &
         .TRUE., clayfraction, "gather", nbp_glo, index_g)
    IF ( ALL(clayfraction(:) .EQ. val_exp ) ) THEN 
       ! clayfraction is not in restart file 
       call_slowproc_soilt=.TRUE. 
    END IF

    var_name= 'sand_frac'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Fraction of sand in each mesh')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., sandfraction, "gather", nbp_glo, index_g)	
    IF ( ALL(sandfraction(:) .EQ. val_exp ) ) THEN 
       ! sandfraction is not in restart file 
       call_slowproc_soilt=.TRUE. 
    END IF

    ! Read siltfrac instead of just recalculating it. It is already in the restart file. 
    ! Recalculating it can lead to a bitwise error unseen by looking at double precision, 
    ! which accumulates and creates a restartability problem in the future.
    var_name= 'silt_frac'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Fraction of silt in each mesh')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., siltfraction, "gather", nbp_glo, index_g)
    IF ( ALL(siltfraction(:) .EQ. val_exp ) ) THEN 
       ! siltfraction is not in restart file 
       call_slowproc_soilt=.TRUE. 
    END IF

    var_name= 'bulk'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Bulk density in each mesh')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, &
         .TRUE., bulk, "gather", nbp_glo, index_g)
    IF ( ALL(bulk(:) .EQ. val_exp ) ) THEN 
       ! bulk is not in restart file 
       call_slowproc_soilt=.TRUE. 
    END IF

    var_name= 'soil_ph'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Soil pH in each mesh')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, &
         .TRUE., soil_ph, "gather", nbp_glo, index_g)
    IF ( ALL(soil_ph(:) .EQ. val_exp ) ) THEN 
       ! soil_ph is not in restart file 
       call_slowproc_soilt=.TRUE. 
    END IF

    IF (impsoilt) THEN 

       !Config Key   = SOIL_FRACTIONS
       !Config Desc  = Areal fraction of the 13 soil USDA textures; the dominant one is selected, Loam by default
       !Config Def   = 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
       !Config If    = IMPOSE_SOILT
       !Config Help  = Determines the fraction of the 13 USDA classes with same order as in constantes_soil_var
       !Config Units = [-]
       soilclass(:,:)=val_exp
       CALL setvar_p (soilclass, val_exp, 'SOIL_FRACTIONS', soilclass_default)
       ! Simplify a heterogeneous grid-cell into an homogeneous one 
       ! with the dominant texture
       njsc(:) = 0
       DO ji = 1, kjpindex
          ! here we reduce to the dominant texture class
          njsc(ji) = MAXLOC(soilclass(ji,:),1)
       ENDDO
       njsc_imp = njsc(1) ! to prescribe the VG parameters consistely with the imposed texture    


       !Config Key   = CLAY_FRACTION
       !Config Desc  = Fraction of the clay fraction (0-dim mode)
       !Config Def   = 0.2
       !Config If    = IMPOSE_SOIL
       !Config Help  = Determines the fraction of clay in the grid box.
       !               If clayfraction was not in restart file it will be read 
       !               from the run.def file instead or initialized based on 
       !               fractions of each textural class
       !Config Units = [-] 
       clayfraction_default =  clayfrac_usda(njsc_imp)
       CALL setvar_p (clayfraction, val_exp, 'CLAY_FRACTION', clayfraction_default)

       !Config Key   = SAND_FRACTION
       !Config Desc  = Fraction of the sand fraction (0-dim mode)
       !Config Def   = 0.4
       !Config If    = IMPOSE_SOIL
       !Config Help  = Determines the fraction of sand in the grid box.
       !Config Units = [-] 
       sandfraction_default =  sandfrac_usda(njsc_imp)
       CALL setvar_p (sandfraction, val_exp, 'SAND_FRACTION',sandfraction_default)

       ! Calculate silt fraction
       siltfraction(:) = 1. - clayfraction(:) - sandfraction(:)

       !Config Key   = BULK
       !Config Desc  = Bulk density (0-dim mode)
       !Config Def   = 1000.0
       !Config If    = IMPOSE_SOIL
       !Config Help  = Determines the bulk density in the grid box.  The bulk density
       !Config         is the weight of soil in a given volume.
       !Config Units = [-] 
       CALL setvar_p (bulk, val_exp, 'BULK', bulk_default)

       !Config Key   = SOIL_PH
       !Config Desc  = Soil pH (0-dim mode)
       !Config Def   = 5.5 
       !Config If    = IMPOSE_SOIL
       !Config Help  = Determines the pH in the grid box.
       !Config Units = [-]
       CALL setvar_p (soil_ph, val_exp, 'SOIL_PH', ph_default)

       !Config Key   = NVAN_IMP
       !Config Desc  = NVAN parameter from Van genutchen equations
       !Config Def   = 1.56 if Loam
       !Config If    = IMPOSE_SOIL
       !Config Help  = Determines the nvan in the grid box.
       !Config Units = [-]
       nvan_default = nvan_usda(njsc_imp)
       CALL setvar_p (nvan, val_exp, 'NVAN_IMP', nvan_default)

       !Config Key   = AVAN_IMP
       !Config Desc  = AVAN parameter from Van genutchen equations
       !Config Def   = 0.0036 if Loam
       !Config If    = IMPOSE_SOIL
       !Config Help  = Determines the avan in the grid box.
       !Config Units = [-]
       avan_default = avan_usda(njsc_imp)
       CALL setvar_p (avan, val_exp, 'AVAN_IMP', avan_default)

       !Config Key   = MCR_IMP
       !Config Desc  = residual soil moisture
       !Config Def   = 0.078 if Loam
       !Config If    = IMPOSE_SOIL
       !Config Help  = Determines the mcr in the grid box.
       !Config Units = [-]
       mcr_default = mcr_usda(njsc_imp)
       CALL setvar_p (mcr, val_exp, 'MCR_IMP', mcr_default)

       !Config Key   = MCS_IMP
       !Config Desc  = saturation soil moisture
       !Config Def   = 0.43 if Loam
       !Config If    = IMPOSE_SOIL
       !Config Help  = Determines the mcs in the grid box.
       !Config Units = [-]
       mcs_default = mcs_usda(njsc_imp)
       CALL setvar_p (mcs, val_exp, 'MCS_IMP', mcs_default)

       !Config Key   = KS_IMP
       !Config Desc  = saturation conductivity
       !Config Def   = 249.6 if Loam
       !Config If    = IMPOSE_SOIL
       !Config Help  = Determines the ks in the grid box.
       !Config Units = [mm/d]
       ks_default = ks_usda(njsc_imp)
       CALL setvar_p (ks, val_exp, 'KS_IMP', ks_default)

       ! By default, we calculate mcf and mcw from the above values, as in slowproc_soilt,
       ! but they can be overruled by values from run.def
       DO ji=1, kjpindex
          mvan(ji) = un - (un / nvan(ji))
       ENDDO

       !
       ! Define matrix potential in mm for wilting point and field capacity (with sand vs clay-silt variation)
       ! no need to define extra value for peatland, because this variable is used to compute mcfc and mcw,
       ! for peatland mcfc and mcw are given
       psi_w(:) = 150000.
       DO ji=1, kjpindex
          IF ( ks(ji) .GE. 560 ) THEN ! Sandy soils (560 is equivalent of 2.75 at log scale of Ks, mm/d)
             psi_fc(ji) = 1000.
          ELSE ! Finer soils
             psi_fc(ji) = 3300.
          ENDIF
       ENDDO

       mcfc_default(:) = mcr(:) + (( mcs(:) - mcr(:)) / (un + ( avan(:) * psi_fc(:))** nvan(:))** mvan(:))
       mcw_default(:)  = mcr(:) + (( mcs(:) - mcr(:)) / (un + ( avan(:) *  psi_w(:))** nvan(:))** mvan(:))
       !
       !Config Key   = MCFC_IMP
       !Config Desc  = field capacity soil moisture
       !Config Def   = 0.1654 if caclulated from default 5 parameters above
       !Config If    = IMPOSE_SOIL
       !Config Help  = Determines the mcfc in the grid box.
       !Config Units = [-]
       mcfc(:) = mcfc_default(:)
       CALL setvar_p (mcfc, val_exp, 'MCFC_IMP', mcfc_default)

       !Config Key   = MCW_IMP
       !Config Desc  = wilting point soil moisture
       !Config Def   = 0.0884 if caclulated from default 5 parameters above
       !Config If    = IMPOSE_SOIL
       !Config Help  = Determines the mcw in the grid box.
       !Config Units = [-]
       mcw(:) = mcw_default(:)
       CALL setvar_p (mcw, val_exp, 'MCW_IMP', mcw_default)

    ELSE IF (.NOT. impsoilt) THEN

       IF ( call_slowproc_soilt ) THEN 
          ! At least one of the output variables from slowproc_soilt were not found in the restart file  
          ! and the user did not want to impose values by making use of the run.def or the default model settings. 
          ! Read a map and initialize. 
          !
          CALL slowproc_soilt(njsc, ks, nvan, avan, mcr, mcs, mcfc, mcw, &
               kjpindex, lalo, neighbours, resolution, &
               contfrac, soilclass, clayfraction, sandfraction, siltfraction, &
               bulk, soil_ph)
          call_slowproc_soilt=.FALSE.
       ENDIF
       
    ELSE
        
       STOP 'The flag IMPOSE_SOILT is not defined correctly'

    ENDIF

    ! XIOS export of Ks before changing the vertical profile 
    CALL xios_orchidee_send_field("ksref",ks) ! mm/d (for CMIP6, once)


    !! 3. Read the infiltration related variables
    ! This variable helps reducing surface runuff in flat areas
    
    !! 3.a Looking first in the restart files
    
    var_name= 'reinf_slope'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Slope coef for reinfiltration')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., reinf_slope, "gather", nbp_glo, index_g)

    !! 3.b We can also read from a map or prescribe, depending on IMPOSE_SLOPE

    IF (impslope) THEN

       !! Impose a constant value from the run.def (not anymore controlled by impveg)
       
       !Config Key   = REINF_SLOPE
       !Config Desc  = Fraction of reinfiltrated surface runoff 
       !Config Def   = 0.1
       !Config If    = IMPOSE_SLOPE
       !Config Help  = Determines the reinfiltration ratio in the grid box due to flat areas
       !Config Units = [-]
       slope_default=0.1
       CALL setvar_p (reinf_slope, val_exp, 'SLOPE', slope_default)

    ELSE

       !! Initialize variables by call to slowproc_slope(reading of map) if not in restart file or if get_slope=T.
     
       !Config Key   = GET_SLOPE
       !Config Desc  = Read slopes from file and do the interpolation
       !Config Def   = n
       !Config If    =
       !Config Help  = Needed for reading the slope file and doing the interpolation. This will be
       !               used by the re-infiltration parametrization
       !Config Units = [FLAG]
       get_slope = .FALSE.
       CALL getin_p('GET_SLOPE',get_slope)

       !! If not found from restart or GET_SLOPE = T, we read from a map 
       IF ( MINVAL(reinf_slope) .EQ. MAXVAL(reinf_slope) .AND. MAXVAL(reinf_slope) .EQ. val_exp .OR. get_slope) THEN
          IF (printlev_loc>=4) WRITE (numout,*) 'reinf_slope was not in restart file or get_slope=T. Now call slowproc_slope'
          
          CALL slowproc_slope(kjpindex, lalo, neighbours, resolution, contfrac, reinf_slope)
          IF (printlev_loc>=4) WRITE (numout,*) 'After slowproc_slope'
          
       ENDIF

    ENDIF

    !! 5. Read the vegetation related variables
    !  The model could start from scratch, from a restart file or the value of 
    !  specific variables could be imposed. Imposing could make use of fixed
    !  values or from maps. This results in many possible configurations some
    !  of which are not very useful. The order of the code in this subroutine 
    !  already limits the number of possible configurations to initialize the
    !  model but the real quality control on this issue is taking place in 
    !  sechiba in the subroutine check_configuration. The order implemented 
    !  here is: (1) Read the values from a restart file if available. If
    !  no restart file is found, give the variable the value val_exp (done
    !  in the subroutine restget_p). If a restart file was found, the value 
    !  for impveg will be ignored (this is taken care of in setvar_p by 
    !  checking whether the variable has the value val_exp or not). If no 
    !  restart file was found there are still two more options to initialize 
    !  the model: (2) initialize the vegetation fractions by fixed values which can 
    !  be set in the run.def or their default values, or (3) initialize the 
    !  vegetation fractions with values from a map or file. In case (2) and (3)
    !  the default of the vegetation varaibles is set to zero. This enables the
    !  users to prescribe the vegetation fractions but start the rest of the
    !  model from scratch. Finally there is 4th option in which the vegetation
    !  variables are also imposed or prescribed from a map. In this case
    !  conflicts may be introduced, for example, the map contains biomass for a
    !  given PFT but previous steps in the initialization resulted in a
    !  veget_max of zero for that PFT. Note that all the information 
    !  required to initialize the vegetation is also stored in the sechiba restart. 
    !  Hence, the value of ok_stomate does not affect this part of the initialization. 

    ! Set default value. It may get overwritten in the subsequent code
    found_restart=.TRUE.

    ! Define element string
    DO iele = 1,nelements
       IF (iele == icarbon) THEN 
          element_str(iele) = 'c' 
       ELSEIF (iele == initrogen) THEN 
          element_str(iele) = 'n' 
       ELSE 
          STOP 'Define element_str' 
       ENDIF
    ENDDO

    !! 5.1 Try to read the values from a sechiba restart file.
    var_name= 'veget'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Vegetation fraction')
    CALL restget_p (rest_id, var_name, nbp_glo, nvm, 1, kjit, .TRUE., &
         veget, "gather", nbp_glo, index_g)
    IF ( ALL( veget(:,:) .EQ. val_exp ) ) found_restart=.FALSE.

    var_name= 'veget_max'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Maximum vegetation fraction')
    CALL restget_p (rest_id, var_name, nbp_glo, nvm, 1, kjit, .TRUE., &
         veget_max, "gather", nbp_glo, index_g)
    IF ( ALL( veget_max(:,:) .EQ. val_exp ) ) found_restart=.FALSE.

    var_name= 'frac_nobio'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Special soil type fraction')
    CALL restget_p (rest_id, var_name, nbp_glo, nnobio, 1, kjit, .TRUE., &
         frac_nobio, "gather", nbp_glo, index_g)
    IF ( ALL( frac_nobio(:,:) .EQ. val_exp ) ) found_restart=.FALSE.

    var_name= 'circ_class_n'
    CALL ioconf_setatt_p('UNITS', 'trees m-2')
    CALL ioconf_setatt_p('LONG_NAME','Stand density')
    CALL restget_p (rest_id, var_name, nbp_glo, nvm, ncirc, kjit, .TRUE., &
         circ_class_n, "gather", nbp_glo, index_g)
    IF ( ALL( circ_class_n(:,:,:) .EQ. val_exp ) ) found_restart=.FALSE.

    DO iele = 1,nelements
       var_name = 'cc_biomass_'//TRIM(element_str(iele))
       CALL ioconf_setatt_p('UNITS', 'gC(N) tree-1')
       CALL ioconf_setatt_p('LONG_NAME','Carbon (or N) mass for the different &
            & biomass components of an individual tree')
       CALL restget_p (rest_id, var_name, nbp_glo, nvm, ncirc, nparts, kjit, .TRUE., &
            circ_class_biomass(:,:,:,:,iele), "gather", nbp_glo, index_g)
       IF ( ALL( circ_class_biomass(:,:,:,:,iele) .EQ. val_exp ) ) found_restart=.FALSE.
    ENDDO

    var_name= 'loss_gain'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Fraction changes during land cover changes')
    CALL restget_p (rest_id, var_name, nbp_glo, nvm, 1, kjit, .TRUE., &
         loss_gain, "gather", nbp_glo, index_g)
    IF (veget_update > 0) THEN
       ! If veget_update is used we expect to have values for loss_gain
       ! in the restart file. If that is not the case there was no restart.
       IF ( ALL( loss_gain(:,:) .EQ. val_exp ) ) found_restart=.FALSE.
       ! If veget_update=0, loss_gain is not used
    END IF

    var_name= 'veget_max_new'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Maximum vegetation fraction after land cover change')
    CALL restget_p (rest_id, var_name, nbp_glo, nvm, 1, kjit, .TRUE., &
         veget_max_new, "gather", nbp_glo, index_g)
    IF (veget_update > 0) THEN
       ! If veget_update is used we expect to have values for veget_max_new
       ! in the restart file. If that is not the case there was no restart.
       IF ( ALL( veget_max_new(:,:) .EQ. val_exp ) ) found_restart=.FALSE.
       ! If veget_update=0, veget_max_new is not used
    END IF

    var_name= 'frac_nobio_new'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Non biological fraction after land cover change')
    CALL restget_p (rest_id, var_name, nbp_glo, nnobio, 1, kjit, .TRUE., &
         frac_nobio_new, "gather", nbp_glo, index_g)
    IF (veget_update > 0) THEN
       ! If veget_update is used we expect to have values for frac_nobio_new
       ! in the restart file. If that is not the case there was no restart.
       IF ( ALL( frac_nobio_new(:,:) .EQ. val_exp ) ) found_restart=.FALSE.
       ! If veget_update=0, frac_nobio_new is not used
    END IF
  
    var_name= 'assim_param'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Assimilation parameters, Vcmax, nue and leaf nitrogen')
    CALL restget_p (rest_id, var_name, nbp_glo, nvm, npco2, kjit, .TRUE., &
         assim_param, "gather", nbp_glo, index_g)
    IF ( ALL( assim_param(:,:,:) .EQ. val_exp ) ) found_restart=.FALSE.

    ! This variable is used to calculate a paremeter that gives the  tree
    ! height when the diameter reaches 1 m. If we set this variable to 0, it 
    ! may take several years to reach a good estimate of the annual
    ! precipitation. The calculation of pipe_tune2 makes use of a minimal
    ! height (min_pipe_tune2). We can use the rain needed to give this
    ! minimal value as an initial value to speed up estimate. We only need
    ! 10mm year-1 to give a pipe_tune2 of 7 m (current value for
    ! min_pipe_tune2). That will not save lots of time. We could reason that
    ! forest are unlikley to appear below 400 mm year-1 and this as the
    ! initial value. Inital tree height in semi-arid regions will be too high
    ! now but of the PFT maps are good there shouldn be too many trees there.
    ! Note that precip_longterm is in mm/day.
    precip_longterm(:) = val_exp                 
    var_name = 'precip_longterm'                 
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., precip_longterm, "gather", nbp_glo, index_g)

    IF (ALL(precip_longterm(:) == val_exp)) THEN
        !Config Key   = PRECIP_LONGTERM_INIT
        !Config Desc  = Pluie annuelle supposee au demarrage a froid, pour initialiser precip_longterm
        !Config If    = OK_STOMATE
        !Config Def   = 700.
        !Config Help  = A zero, pipe_tune2 est rabattu sur son plancher min_pipe_tune2 (7 m) le
        !Config         premier jour. Or c'est ce jour-la que prescribe dimensionne les peuplements :
        !Config         la biomasse posee est ensuite relue avec l'allometrie climatique reelle
        !Config         (pipe_tune2 ~30), le diametre est divise par ~1.8, Nmax triple et le RDI
        !Config         tombe de 0.70 a ~0.23 sans qu'aucun arbre ne meure. Valeur par defaut :
        !Config         les forets apparaissent rarement sous 400 mm/an (cf. commentaire ci-dessus),
        !Config         mais 700 mm/an colle mieux au domaine europeen vise ici. Un scalaire global
        !Config         reste approximatif : un site atlantique a 1200 mm/an garde un residu.
        !Config         Sans effet en reprise sur fichier : cold start uniquement.
        !Config Units = [mm/an]
        precip_longterm_init = 700.
        CALL getin_p('PRECIP_LONGTERM_INIT', precip_longterm_init)
        ! Guillaume M. -- precip_longterm is held in mm/day, the namelist value in mm/yr.
        precip_longterm(:) = precip_longterm_init / 365.
        found_restart=.FALSE.
        ! Guillaume M. -- In forced mode forcing_meanprecip already holds the PER-GRID-CELL
        ! forcing mean; it wins over the scalar, which is only a fallback (coupled mode, or
        ! forcing without a readable precipitation variable). Grid cells left at zero keep
        ! the scalar: a zero value would reproduce the very defect this code fixes.
        IF ( ALLOCATED(precip_forcing_mean) ) THEN
           IF ( SIZE(precip_forcing_mean) == kjpindex ) THEN
              WHERE ( precip_forcing_mean(:) > zero ) precip_longterm(:) = precip_forcing_mean(:)
              ! Guillaume M. -- The COUNT intrinsic is unusable here: a local REAL named
              ! `count` shadows it in this routine, so COUNT(...) would be parsed as an
              ! array indexed by a logical. Hence SUM(MERGE(...)).
              WRITE(numout,*) 'slowproc_init : cold start, precip_longterm initialise PAR MAILLE ', &
                   'depuis le forcage ; mailles retombant sur le scalaire : ', &
                   SUM(MERGE(1,0,precip_forcing_mean(:) <= zero)), ' / ', kjpindex
           ELSE
              CALL ipslerr_p(2,'slowproc_init', &
                   'precip_forcing_mean a une taille inattendue : repli sur le scalaire', &
                   'PRECIP_LONGTERM_INIT','')
           ENDIF
        ENDIF
        WRITE(numout,*) 'slowproc_init : cold start, precip_longterm (mm/an) min/moy/max = ', &
             MINVAL(precip_longterm(:))*365., &
             SUM(precip_longterm(:))/REAL(MAX(kjpindex,1),r_std)*365., &
             MAXVAL(precip_longterm(:))*365.
    ENDIF

    !+++CHECK+++
!!$    ! Are the values of frac_age what we expect for ok_bvoc? Seems that not 
!!$    ! much is happening with frac_age for the moment 
!!$    ! frac_age is used in ok_bvoc which can be called if ok_stomate = F.
!!$    ! In case stomate is not used, it needs to be read from the sechiba 
!!$    ! restart file. 
!!$    CALL ioconf_setatt_p('UNITS', '-')
!!$    CALL ioconf_setatt_p('LONG_NAME','Fraction of leaves in leaf age class ')
!!$    CALL restget_p (rest_id, 'frac_age', nbp_glo, nvm, nleafages, kjit, .TRUE., &
!!$         frac_age, "gather", nbp_glo, index_g)
!!$    IF ( ALL( frac_age(:,:,:) .EQ. val_exp ) ) found_restart=.FALSE.
    !+++++++++++
    
    !! 5.2 Impose vegetation fractions or read from a PFT map
    IF (.NOT. found_restart .AND. impveg) THEN
       
       ! impveg=TRUE: there can not be any land use change, veget_update must be =0
       ! Read VEGET_UPDATE from run.def and exit if it is different from 0Y
       IF (veget_update /= 0) THEN
          WRITE(numout,*) 'veget_update=',veget_update,' is not coherent with impveg=',impveg
          CALL ipslerr_p(3,'slowproc_init','Incoherent values between impveg and veget_update', &
               'veget_update must be equal to 0 if impveg=true','')
       ENDIF

       ! If no restart file was found replace the val_exp by zero so
       ! the calculations can start.
       IF (ALL(circ_class_biomass(:,:,:,:,:).EQ.val_exp) ) THEN
          circ_class_biomass(:,:,:,:,:) = zero
          circ_class_n(:,:,:) = zero
       ENDIF

       ! Initialize the vegetation fractions by reading run.def. Previously setvar_p was
       ! used but it swapped the dimensions when kjpindex = nvm (thus when using 15 pixels).
       ! Changed to getin_p which is typically used for PFT-dependent parameters.
       !Config Key   = SECHIBA_VEGMAX
       !Config Desc  = Maximum vegetation distribution within the mesh (0-dim mode)
       !Config If    = IMPOSE_VEG
       !Config Def   = 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.8, 0.0, 0.0, 0.0
       !Config Help  = The fraction of vegetation is read from the restart file. If
       !Config         it is not found there we will use the values provided here.
       !Config Units = [-]
       sechiba_vegmax(:) = 0.
       sechiba_vegmax(1) = 1.
       CALL getin_p('SECHIBA_VEGMAX',sechiba_vegmax)
       DO ivm = 1, nvm
          veget_max(:,ivm) = sechiba_vegmax(ivm)
       END DO
       
       !Config Key   = SECHIBA_FRAC_NOBIO
       !Config Desc  = Fraction of other surface types within the mesh (0-dim mode)
       !Config If    = IMPOSE_VEG
       !Config Def   = 0.0
       !Config Help  = The fraction of ice, lakes, etc. is read from the restart file. If
       !Config         it is not found there we will use the values provided here.
       !Config         For the moment, there is only ice.
       !Config Units = [-]
       ! The routine setvar_p will only initialize the variable if it was not found in restart file.
       frac_nobio1 = frac_nobio(1,1)
       CALL setvar_p (frac_nobio1, val_exp, 'SECHIBA_FRAC_NOBIO', frac_nobio_fixed_test_1)
       frac_nobio(:,:) = frac_nobio1

       ! Check for precision errors
       count(:) = zero
       ! note that the second dimension of veget_max is nvm but for frac_nobio it is nnobio
       fracsum(:) = SUM(veget_max(:,:),2) + SUM(frac_nobio(:,:),2)
       WHERE (ABS(fracsum(:)-un).GT.10*EPSILON(un)) 
          count(:) = un
       ENDWHERE
       IF (SUM(count(:)).GT.zero) THEN
          ! At least one error was found
          DO ipts = 1,kjpindex
             IF (ABS(fracsum(ipts)-un).GT.min_stomate) THEN
                ! Too big for a precision error
                WRITE(numout,*) 'ipts, veget_max, frac_nobio, total, ', ipts, &
                     SUM(veget_max(ipts,:)), SUM(frac_nobio(ipts,:)), fracsum(ipts)
                CALL flush(numout)
                CALL ipslerr_p(3,'slowproc.f90','The prescribed vegetation does not',&
                     'add up to 1','check the run.def')
             ELSEIF ( ABS(fracsum(ipts)-un).LE.min_stomate .AND. &
                  ABS(fracsum(ipts)-un).GT.10*EPSILON(un) ) THEN
                ! Looks like a precision error so correct it
                residual = fracsum(ipts) - un
                ! Note that dim=1 refers to the dimensions of the answer
                ivmax = MAXLOC(veget_max(ipts,:),DIM=1)
                IF (veget_max(ipts,ivmax).GT.residual) THEN
                   veget_max(ipts,ivmax) = veget_max(ipts,ivmax) - residual
                   CALL ipslerr_p(2,'slowproc.f90','correcting the prescribed veget_max',&
                        'to account for a precision issue','')
                ELSE
                   CALL ipslerr_p(3,'slowproc.f90','Could not correct a precision error',&
                        'check the code because this is unecpected','') 
                END IF
             END IF 
          END DO ! ipts
       END IF ! at least one error
       
    ELSE

       ! The DGVM cannot be combined with land use change if agriculture is not 
       ! accounted for by reading a land cover change map (agriculture cannot
       ! be predicted by the dgvm and therefore needs to be prescribed by a map)
       IF (veget_update > 0 .AND. ok_dgvm .AND. .NOT. agriculture) THEN
          CALL ipslerr_p(3,'slowproc_init',&
               'The combination DGVM=TRUE, AGRICULTURE=FALSE and VEGET_UPDATE>0 is not possible', &
               'Set VEGET_UPDATE=0Y in run.def','')
       END IF  

       IF ( .NOT. found_restart .OR. vegetmap_reset) THEN

          ! Only when there is no restart and the values have not been imposed, the model
          ! will read the vegetation fractions from a file
          IF (printlev_loc>=3) WRITE(numout,*) 'Before call read_vegetmax in &
               & initialization phase without restart files'

          ! Call the routine to read the vegetation from file (output is veget_max_new)
          CALL slowproc_readvegetmax(kjpindex, lalo, neighbours, resolution, contfrac, &
               veget_max, veget_max_new, frac_nobio_new, .TRUE.)
 
          ! Check the map that was just read for errors.
          CALL check_read_vegetmax(kjpindex, veget_max_new, frac_nobio_new)
          
          ! Replace the val_exp (needed to check whether a restart was found) 
          ! by zero so the calculations can start
          IF (ALL(circ_class_biomass(:,:,:,:,:).EQ.val_exp) ) THEN
             circ_class_biomass(:,:,:,:,:) = zero
             circ_class_n(:,:,:) = zero
             veget_max(:,:) = veget_max_new(:,:)
             frac_nobio(:,:) = frac_nobio_new(:,:)
          ENDIF

          !! Reset totaly or partialy veget_max if using DGVM
          IF ( ok_dgvm  ) THEN

             ! If we are dealing with dynamic vegetation then all natural PFTs should 
             ! be set to veget_max = 0. In case no agriculture is desired, agriculture 
             ! PFTS should be set to 0 as well
             IF (agriculture) THEN

                DO jv = 2, nvm
                   IF (natural(jv)) THEN 
                      veget_max(:,jv)=zero
                   ENDIF
                ENDDO

                ! Calculate the fraction of crop for each point.
                ! Sum only on the indexes corresponding to the non_natural pfts
                frac_crop_tot(:) = zero
                DO jv = 2, nvm
                   IF(.NOT. natural(jv)) THEN
                      DO ji = 1, kjpindex
                         frac_crop_tot(ji) = frac_crop_tot(ji) + veget_max(ji,jv)
                      ENDDO
                   ENDIF
                END DO

                ! Calculate the fraction of bare soil
                DO ji = 1, kjpindex
                   veget_max(ji,1) = un - frac_crop_tot(ji) - SUM(frac_nobio(ji,:))   
                ENDDO

             ELSE
                
                ! No agriculture land in this simulation
                veget_max(:,:) = zero
                DO ji = 1, kjpindex
                   veget_max(ji,1) = un  - SUM(frac_nobio(ji,:))
                ENDDO

             END IF ! agriculture is considered in the DGVM

          END IF  ! end ok_dgvm
     
       END IF ! No restart was found

    END IF ! impose vegetation fractions

    !! Continue initializing variables not found in restart file. Case for both impveg=true and false.
    ALLOCATE(cn_leaf_min_2D(kjpindex, nvm), STAT=ier)
    ALLOCATE(cn_leaf_max_2D(kjpindex, nvm), STAT=ier)
    ALLOCATE(cn_leaf_init_2D(kjpindex, nvm), STAT=ier)

    IF (impose_cn .AND. read_cn) THEN
       ! cn_leaf_min_2D, cn_leaf_max_2D and cn_leaf_init_2D are set in slowproc_readcnleaf by reading a map in slowproc_readcnleaf
       ! Note that they are not explicitly passed to this subroutine, but they are module variables, and slowproc_readcnleaf is
       ! still in this module, so the values are changed nonetheless.
       CALL slowproc_readcnleaf(kjpindex, lalo, neighbours, resolution, contfrac)
    ELSE
       ! cn_leaf_min_2D, cn_leaf_max_2D and cn_leaf_init_2D take scalar values with constant spatial distribution
       DO ji=1,kjpindex
          cn_leaf_min_2D(ji,:)=cn_leaf_min(:)
          cn_leaf_init_2D(ji,:)=cn_leaf_init(:)
          cn_leaf_max_2D(ji,:)=cn_leaf_max(:)
       ENDDO
    ENDIF

    !! 5.3 Initialize vegetation characteristics
    !  If a restart file has been read, the vegetation characteristics are
    !  already known. Note that reading the vegetation characteristics is 
    !  only controlled by the read_lai flag.The user could still decide to 
    !  use all information from the restart but overwrite circ_class_biomass 
    !  and circ_class_n. If no restart file has been read this far, sechiba 
    !  will need a description of the canopy to be able to run. This description 
    !  should come from a so called lai map (note that the map no longer contains
    !  lai but contains information on the biomass and number of individuals)

    IF (.NOT. ok_stomate) THEN
       IF (read_lai==0 .OR. read_lai==2) THEN 
          ! Initialize
          found_restart = .TRUE.
          
          ! Read "canopy structure map". This map should be based on
          ! a previous simulation with ORCHIDEE.
          ! The difference with a normal restart is that it should read
          ! circ_class_biomass and circ_class_n for 12 months.
          var_name= 'cc_biomass_monthly'
          CALL ioconf_setatt_p('UNITS', 'g C(N) m-2 tree-1')
          CALL ioconf_setatt_p('LONG_NAME','Monthly values for biomass &
               & components per circumference class')
          DO iele = 1,nelements
             DO imonth = 1,nmonth
10              FORMAT(I2)
                WRITE (month_str,10) imonth
                var_name = 'cc_biomass_m_'//TRIM(month_str)//'_'//TRIM(element_str(iele))
                CALL restget_p (rest_id, var_name, nbp_glo, nvm, ncirc, nparts, kjit, .TRUE., &
                     cc_biomass_m(:,:,:,:,imonth,iele),"gather", nbp_glo, index_g)
                IF ( ALL( cc_biomass_m(:,:,:,:,imonth,iele) .EQ. val_exp ) ) found_restart=.FALSE.
             ENDDO
          ENDDO
          
          var_name = 'cc_n_m'
          CALL ioconf_setatt_p('UNITS', 'tree-1 m-2')
          CALL ioconf_setatt_p('LONG_NAME','trees per m2 for each circ class')
          CALL restget_p (rest_id, var_name, nbp_glo, nvm, ncirc, 12, kjit, .TRUE. , &
               cc_n_m, "gather", nbp_glo, index_g)
          IF ( ALL( cc_n_m(:,:,:,:) .EQ. val_exp) ) found_restart=.FALSE.
          
          IF ( .NOT. found_restart ) THEN           
             ! Initialize other canopy variables with zeros. 
             ! setvar checks whether the restart was found. 
             CALL setvar_p (frac_age, val_exp, 'FRAC_AGE',zero)
             CALL setvar_p (loss_gain, val_exp, 'LOSS_GAIN', zero)
             CALL setvar_p (veget_max_new, val_exp, 'VEGET_MAX_NEW', zero)
             CALL setvar_p (frac_nobio_new, val_exp, 'FRAC_NOBIO_NEW', zero)
             
    
             !  If no restart was found, the model will need to make up its assimilation
             !  parameters to calculate photosynthesis and transpiration. Note that
             !  circ_class_biomass should be available either through a restart file,
             !  or an lai map (that no longer contains lai but contains biomass and
             !  and the number of individuals.
             
             ! Initialize the variables revelant for the assimilation parameters
             DO jv = 1, nvm
                assim_param(:,jv,ivcmax) = vcmax_fix(jv)
                assim_param(:,jv,inue) = nue_opt(jv)
                assim_param(:,jv,ileafN) = SUM( &
                     circ_class_biomass(:,jv,:,ileaf,initrogen) * &
                     circ_class_n(:,jv,:),2 )
             ENDDO
          ENDIF ! .NOT. found_restart

       END IF
       SELECT CASE (read_lai)
       CASE(0)
          !! Read_lai == 0 corresponds to global level lai prescribe. 
          IF ( .NOT. found_restart ) THEN

             ! cc_biomass_m and/or cc_n_m were not found in the restart file
             ! Now prescribe them
             IF (ncirc.NE.3) THEN
                WRITE(numout,*) 'Current number of circumference classes (ncirc), ', ncirc
                CALL ipslerr_p(3,'read_lai=0 only works with ncirc=3',&
                     'data is imposed on the entire domain', '','')
             END IF

             circ_class_biomass(:,:,:,:,:) = zero
             circ_class_n(:,:,:) = zero
             DO ji = 1,kjpindex
                DO jv = 2,nvm
                   circ_class_biomass(ji,jv,1,:,icarbon) = 3*(/64.1790961402884, 1219.15235836900, 564.11952303, &
                        712.783132227959, 417.796173293409, 122.105778809076, 20.5181600482659, &
                        337.744066746780, 46.6473441173028/)     
                   circ_class_biomass(ji,jv,1,:,initrogen) = 3*(/5.38391872316575,14.5306615478275, 6.32728500214539, &
                        6.73356880989017, 3.63140992343535, 8.80379728708368, 0.909221528487578, &     
                        4.255487191806592E-004,  6.783452513144551E-003/)
                   circ_class_biomass(ji,jv,2,:,icarbon) =  3*(/285.559552359068, 8632.15052183103, 3862.94977896271, &
                        4440.95806491165, 2462.91756260235, 543.299510808887, 20.5181600482660, &     
                        337.744066746985, 99.9253437952986/)     
                   circ_class_biomass(ji,jv,2,:,initrogen) = 3*(/10.8985571313903, 36.6657154882537, 15.9658547505304, &
                        16.9910445846784, 9.16326091799929, 17.8213477283047, 3.105548985182576E-003, &
                        4.255487191808893E-004, 1.603616996966766E-002/)
                   circ_class_biomass(ji,jv,3,:,icarbon) =  3*(/1358.39825416204, 66727.3342051582, 28610.5603591189, &
                        28941.8798532254, 15160.4663937242, 2584.45953172628, 20.5181600482660, &
                        337.744066746983, 180.177367238784/)     
                   circ_class_biomass(ji,jv,3,:,initrogen) = 3*(/18.7135366977149, 74.2989491440166, 32.3530091900281, &
                        34.4304410998095, 18.5683177833507, 30.6004217526930, 1.039620165307951E-005, &
                        4.255487191808977E-004, 5.568363409712271E-002/)
                   circ_class_n(ji,jv,:) = (/0.373684875218478, 0.124561625072826, 4.152054169094205E-002/)
                END DO
             END DO
          ENDIF ! .NOT. found_restart

       CASE(1)
          !! Read_lai == 1 corresponds to site level lai reading. 
          !! Everything is read from the run.def and prescribed by the user.

          !! When LAI is read at the site level, only one circumference class is considered
          IF (ncirc.NE.1) THEN
             WRITE(numout,*) 'Current number of circumference classes (ncirc), ', ncirc
             CALL ipslerr_p(3,'instead of reading a LAI map a handfull of fixed biomass',&
                  'data is read and set for the entire domain.', &
                  'These values excpect that ncirc = 1','This seems not to be the case')
          END IF
          
          DO ji = 1, kjpindex
             DO jv = 1, nvm
                circ_class_biomass(ji,jv,ncirc,:,initrogen) = 3*(/18.7135366977149, 74.2989491440166, 32.3530091900281, &
                     34.4304410998095, 18.5683177833507, 30.6004217526930, 1.039620165307951E-005, &
                     4.255487191808977E-004, 5.568363409712271E-002/)
             ENDDO
          ENDDO
          circ_class_n(:,:,:) = canopy_density_prescribed
          !! Then, everything is calculated in the bypass stomate

       CASE(2)          
          !! Read_lai == 2 corresponds to global level biomass reading. 
           IF ( .NOT. found_restart ) THEN
             CALL ipslerr_p(3,'Option read_lai=2 is not yet implemented.',&
                  'Creating a biomass map is first needed', &
                  'Then a subroutin to read and interpolate the map needs to be implemented.', &
                  'If short of inspiration, see ORCHIDEE_2 and slowproc_interlai')
          ENDIF ! .NOT. found_restart

       END SELECT ! read_lai     
    ENDIF ! ok_stomate == FALSE

    
    !+++CHECK+++
    ! Special case for DGVM and restart file. 
    ! JG why is specific treatement needed for DGVM ? Is not the veget_max variable 
    ! correct in the end of last run ?
    IF (found_restart) THEN

       ! WITH restarts for vegetation and DGVM and NO AGRICULTURE
       IF ( ok_dgvm  .AND. .NOT. agriculture ) THEN
          ! Calculate the total fraction of crops for each point
          frac_crop_tot(:) = zero
          DO jv = 2, nvm  
             IF ( .NOT. natural (jv))  THEN
                DO ji = 1, kjpindex
                   frac_crop_tot(ji) = frac_crop_tot(ji) + veget_max(ji,jv)
                ENDDO
             ENDIF
          ENDDO
          
          ! Add the crops fraction to the bare soil fraction
          DO ji = 1, kjpindex
             veget_max(ji,1) = veget_max(ji,1) + frac_crop_tot(ji)
          ENDDO
          
          ! Set the crops fraction to zero
          DO jv = 2, nvm                  
             IF ( .NOT. natural (jv))  THEN
                veget_max(:,jv) = zero
             ENDIF
          ENDDO
       ENDIF
    ENDIF ! end found_restart
    !+++++++++++

    !! 6. Dynamic irrigation map
    !  If do_irrigation, it will look to the dynamical irrig. map in restart
    !  If not dynamic irrig. map, it will be set to zero.
    !  If not found in restart, it will try to interpolate the map
    
    irrigated_new(:) = zero !
    IF ( do_irrigation ) THEN
       ! It will look into restart file
       var_name = 'irrigmap_dyn'
       CALL ioconf_setatt_p('UNITS', 'm2')
       CALL ioconf_setatt_p('LONG_NAME','Dynamical area equipped for irrigation')
       CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., irrigated_new, "gather", nbp_glo, index_g)
       
       ! Now, if not found in the restart, read and interpolate from file
       IF ( ALL( irrigated_new(:) .EQ. val_exp ) ) THEN
          CALL slowproc_readirrigmap_dyn(kjpindex, lalo, neighbours,  resolution, contfrac,         &
               irrigated_new)
       ENDIF
       ! irrigated_next from sechiba (int output) = irrigated_new  in slowproc_initialize.
    ENDIF
    
    ! If new priorization scheme is used, it also seek into restart, or interpolate
    
    fraction_aeirrig_sw(:) = un
    IF ( do_irrigation .AND. select_source_irrig) THEN
       ! It will look into restart file
       var_name = 'fraction_aeirrig_sw'
       CALL ioconf_setatt_p('UNITS', '%')
       CALL ioconf_setatt_p('LONG_NAME','Fraction of area equipped for irrigation with surface water')
       CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., fraction_aeirrig_sw, "gather", nbp_glo, index_g)
       
       ! Now, if not found in the restart, read and interpolate from file
       IF ( ALL( fraction_aeirrig_sw(:) .EQ. val_exp ) ) THEN
          CALL slowproc_read_aeisw_map(kjpindex, lalo, neighbours,  resolution, contfrac,         &
               fraction_aeirrig_sw)
       ENDIF
       ! irrigated_next from sechiba (int output) = irrigated_new in slowproc_initialize.
    ENDIF

    
    !! 7. Some calculations always done, with and without restart files
    !  The variables veget, veget_max and frac_nobio were all read from 
    !  restart file or initialized above. Calculate now totfrac_nobio and
    !  soiltiles using these variables.
    IF (ok_bare_soil_new) THEN
       
       !+++CHECK+++
       ! We no longer want to treat the gaps in the canopy as
       ! bare soil. It needs to be tested what will happen with
       ! the evaporation in the single-layer model. The multi-
       ! layer energy budget should be able to correctly deal
       ! with the gaps in the canopy.
       tot_bare_soil(:) = veget_max(:,1)
       
       ! Total frac nobio
       totfrac_nobio(:) = SUM(frac_nobio(:,:),2)
       !+++++++++++
       
    ELSE

       ! Initialize
       fracsum(:) = zero
       tot_bare_soil(:) = veget_max(:,1)

       ! Calculate bare soil fraction
       DO ji = 1, kjpindex

          ! Total frac nobio
          totfrac_nobio(ji) = SUM(frac_nobio(ji,:))
          
          DO jv = 2, nvm

             ! Move the canopy gaps into the bare soil
             ! fraction for this time step. Calculate
             ! the total fraction of bare soil in the
             ! grid
             tot_bare_soil(ji) = tot_bare_soil(ji) + &
                  (veget_max(ji,jv) - veget(ji,jv))
             fracsum(ji) = fracsum(ji) + veget(ji,jv)

          ENDDO

          ! Consistency check
          fracsum(ji) = fracsum(ji) + tot_bare_soil(ji) + SUM(frac_nobio(ji,:)) 
          IF (fracsum(ji) .LT. 0.99999) THEN
             WRITE(numout,*)' ATTENTION, in ji, fracsum LT 1: ', ji, fracsum(ji)
             WRITE(numout,*)'  frac_nobio = ',SUM(frac_nobio(ji,:))
             WRITE(numout,*)'  veget = ',veget(ji,:)
             WRITE(numout,*)'  tot_bare_soil = ',tot_bare_soil(ji)
          ENDIF
       ENDDO
       
    ENDIF ! ok_bare_soil_new
   
    !! 8. Calculate soiltiles
    !  Soiltiles are only used in hydrol, but we fix them in here because some time 
    !  it might depend on a changing vegetation (but then some adaptation should be 
    !  made to hydrol) and be also used in the other modules to perform separated 
    !  energy balances.
    !  The sum of all soiltiles makes one, and corresponds to the bio fraction
    !  of the grid cell (called vegtot in hydrol)   
    soiltile(:,:) = zero
    DO jv = 1, nvm
       jst = pref_soil_veg(jv)
       DO ji = 1, kjpindex
             soiltile(ji,jst) = soiltile(ji,jst) + veget_max(ji,jv)
       ENDDO
    ENDDO

    DO ji = 1, kjpindex 
       IF (totfrac_nobio(ji) .LT. (1-EPSILON(un))) THEN
          ! If the calculation of 1-totfrac_nobio is correct its value
          ! should be identical to SUM(soiltile) give or take precision
          ! issues.
          IF (ABS(SUM(soiltile(ji,:))-(1-totfrac_nobio(ji))).LT.min_sechiba) THEN
             ! If the numbers are very similar, the divsion may result in 
             ! errors of 10e-10. By taking the minimum we protect against
             ! such conditions. 
             soiltile(ji,:)=MIN(soiltile(ji,:)/(1.-totfrac_nobio(ji)),un)
          ELSE
             ! Mismatch between SUM(soiltile) and 1-totfrac_nobio.
             CALL ipslerr_p (3,'slowproc_main',&
                  'Mismatch between SUM(soiltile) and 1-totfrac_nobio','','')
          ENDIF
       ELSE
          soiltile(ji,:)=zero
       ENDIF
    ENDDO  

    !! 10. Initialize different sources of nitrogen
    IF(ok_ncycle .AND. (.NOT. impose_CN)) THEN
       IF((.NOT. impose_ninput_dep) .OR. (.NOT. impose_ninput_fert) .OR. (.NOT. impose_ninput_bnf)) THEN
          var_name= 'Ninput_year'
          CALL ioconf_setatt_p('UNITS', '-')
          CALL ioconf_setatt_p('LONG_NAME','Last year get in N input file.')
          ! Read Ninput_year from restart file. For restget_p interface for scalar value, the default value 
          ! if the variable is not in the restart file is given as argument, here use REAL(Ninput_year_orig)
          CALL restget_p (rest_id, var_name, kjit, .TRUE., REAL(Ninput_year_orig), Ninput_year)

          IF (Ninput_reinit) THEN
             ! Reset Ninput_year
             Ninput_year=Ninput_year_orig
          ENDIF
       
          
          !Config Key   = NINPUT_UPDATE
          !Config Desc  = Update N input frequency
          !Config If    = ok_ncycle .AND. (.NOT. impose_cn) .AND. .NOT. impsoilt
          !Config Def   = 0Y
          !Config Help  = The veget datas will be update each this time step.
          !Config Units = [years]
          !
          ninput_update=0
          WRITE(ninput_str,'(a)') '0Y'
          CALL getin_p('NINPUT_UPDATE', ninput_str)
          l=INDEX(TRIM(ninput_str),'Y')
          READ(ninput_str(1:(l-1)),"(I2.2)") ninput_update
          IF (printlev>=1) WRITE(numout,*) "Update frequency for N inputs in years :",ninput_update
       ENDIF

       IF(.NOT. impose_Ninput_dep) THEN
          FOUND_RESTART=.TRUE.
          CALL ioconf_setatt_p('UNITS', 'kgN m-2 yr-1')
          CALL ioconf_setatt_p('LONG_NAME','N ammonium deposition')
          CALL restget_p (rest_id, 'Nammonium', nbp_glo, nvm, 12, kjit, .TRUE., N_input_raw(:,:,:,iammonium), &
                  "gather", nbp_glo, index_g)
          IF ( ALL( N_input_raw(:,:,:,iammonium) .EQ. val_exp ) ) FOUND_RESTART=.FALSE.

          CALL ioconf_setatt_p('UNITS', 'kgN m-2 yr-1')
          CALL ioconf_setatt_p('LONG_NAME','N nitrate deposition')
          CALL restget_p (rest_id, 'Nnitrate', nbp_glo, nvm, 12, kjit, .TRUE., N_input_raw(:,:,:,initrate), &
               "gather", nbp_glo, index_g)
          IF ( ALL( N_input_raw(:,:,:,initrate) .EQ. val_exp ) ) FOUND_RESTART=.FALSE.

          IF(.NOT. FOUND_RESTART) THEN
             ! Read the new N inputs from file. Output is Ninput and frac_nobio_nextyear.
             fieldname='Nammonium'
             CALL slowproc_Ninput(kjpindex, lalo, neighbours, resolution, contfrac, fieldname, &
                  N_input_temp, Ninput_year)
             ! Conversion from mgN/m2/yr to gN/m2/day
             N_input_raw(:,:,:,iammonium)=N_input_temp/1000./one_year

             fieldname='DRYNHX'
             CALL slowproc_Ninput(kjpindex, lalo, neighbours, resolution, contfrac, fieldname, &
                  N_input_temp, Ninput_year)
             ! Conversion from kgN/m2/s to gN/m2/day
             N_input_raw(:,:,:,iammonium)=N_input_raw(:,:,:,iammonium)+N_input_temp*1000.*one_day

             fieldname='WETNHX'
             CALL slowproc_Ninput(kjpindex, lalo, neighbours, resolution, contfrac, fieldname, &
                  N_input_temp, Ninput_year)
             ! Conversion from kgN/m2/s to gN/m2/day
             N_input_raw(:,:,:,iammonium)=N_input_raw(:,:,:,iammonium)+N_input_temp*1000.*one_day

             fieldname='Nnitrate'
             CALL slowproc_Ninput(kjpindex, lalo, neighbours, resolution, contfrac, fieldname, &
                  N_input_temp, Ninput_year)
             ! Conversion from mgN/m2/yr to gN/m2/day
             N_input_raw(:,:,:,initrate)=N_input_temp/1000./one_year

             fieldname='DRYNOY'
             CALL slowproc_Ninput(kjpindex, lalo, neighbours, resolution, contfrac, fieldname, &
                  N_input_temp, Ninput_year)
             ! Conversion from kgN/m2/s to gN/m2/day
             N_input_raw(:,:,:,initrate)=N_input_raw(:,:,:,initrate)+N_input_temp*1000.*one_day

             fieldname='WETNOY'
             CALL slowproc_Ninput(kjpindex, lalo, neighbours, resolution, contfrac, fieldname, &
                  N_input_temp, Ninput_year)
             ! Conversion from kgN/m2/s to gN/m2/day
             N_input_raw(:,:,:,initrate)=N_input_raw(:,:,:,initrate)+N_input_temp*1000.*one_day
          ENDIF
       ELSE
             !Config Key   = NAMMONIUM
             !Config Desc  = Amount of N ammonium deposition 
             !Config Def   = 0
             !Config If    = ok_ncycle .AND. (.NOT. impose_cn)
             !Config Help  = 
             !Config Units = [gN m-2 d-1] 
             nammonium=zero
             CALL getin_p('NAMMONIUM',nammonium)       
             n_input_raw(:,:,:,iammonium)=nammonium
             !Config Key   = NNITRATE
             !Config Desc  = Amount of N nitrate deposition 
             !Config Def   = 0
             !Config If    = ok_ncycle .AND. (.NOT. impose_cn)
             !Config Help  = 
             !Config Units = [gN m-2 d-1] 
             nnitrate=zero
             CALL getin_p ('NNITRATE',nnitrate)       
             n_input_raw(:,:,:,initrate)=nnitrate
       ENDIF


       IF(.NOT. impose_Ninput_fert) THEN
          FOUND_RESTART=.TRUE.

          CALL ioconf_setatt_p('UNITS', 'kgN m-2 yr-1')
          CALL ioconf_setatt_p('LONG_NAME','N fertilizer')
          CALL restget_p (rest_id, 'Nfert_nitr', nbp_glo, nvm, 12, kjit, .TRUE., N_input_raw(:,:,:,ifert_nitr), "gather", nbp_glo, index_g)
          IF ( ALL( N_input_raw(:,:,:,ifert_nitr) .EQ. val_exp ) ) FOUND_RESTART=.FALSE.
          CALL restget_p (rest_id, 'Nfert_ammo', nbp_glo, nvm, 12, kjit, .TRUE., N_input_raw(:,:,:,ifert_ammo), "gather", nbp_glo, index_g)
          IF ( ALL( N_input_raw(:,:,:,ifert_ammo) .EQ. val_exp ) ) FOUND_RESTART=.FALSE.


          IF(.NOT. FOUND_RESTART) THEN
             ! Read the new N inputs from file. Output is Ninput and frac_nobio_nextyear.
             fieldname='Nfert'
             CALL slowproc_Ninput(kjpindex, lalo, neighbours, resolution, contfrac, fieldname, &
                  N_input_temp, Ninput_year)
             ! Conversion from gN/m2(cropland)/yr to gN/m2/day
             N_input_raw(:,:,:,ifert_ammo) = N_input_temp(:,:,:)/one_year * ratio_nh4_fert
             N_input_raw(:,:,:,ifert_nitr) = N_input_temp(:,:,:)/one_year * (1.-ratio_nh4_fert)

             fieldname='Nfert_cropland'
             CALL slowproc_Ninput(kjpindex, lalo, neighbours, resolution, contfrac, fieldname, &
                  N_input_temp, Ninput_year)
             ! Conversion from gN/m2(cropland)/yr to gN/m2/day
             N_input_raw(:,:,:,ifert_ammo ) = N_input_raw(:,:,:,ifert_ammo)+ N_input_temp(:,:,:)/one_year * ratio_nh4_fert
             N_input_raw(:,:,:,ifert_nitr) = N_input_raw(:,:,:,ifert_nitr)+ N_input_temp(:,:,:)/one_year * (1.-ratio_nh4_fert)
             
             fieldname='Nfert_ammo_cropland'
             CALL slowproc_Ninput(kjpindex, lalo, neighbours, resolution, contfrac, fieldname, &
                  N_input_temp, Ninput_year)
             ! Conversion from gN/m2(cropland)/yr to gN/m2/day
             N_input_raw(:,:,:,ifert_ammo ) = N_input_raw(:,:,:,ifert_ammo)+ N_input_temp(:,:,:)/one_year

             fieldname='Nfert_nitr_cropland'
             CALL slowproc_Ninput(kjpindex, lalo, neighbours, resolution, contfrac, fieldname, &
                  N_input_temp, Ninput_year)
             ! Conversion from gN/m2(cropland)/yr to gN/m2/day
             N_input_raw(:,:,:,ifert_nitr) = N_input_raw(:,:,:,ifert_nitr)+ N_input_temp(:,:,:)/one_year


             fieldname='Nfert_cropC3'
             CALL slowproc_Ninput(kjpindex, lalo, neighbours, resolution, contfrac, fieldname, &
                  N_input_temp, Ninput_year)
             ! Conversion from gN/m2(cropland)/yr to gN/m2/day
             N_input_raw(:,:,:,ifert_ammo) = N_input_raw(:,:,:,ifert_ammo)+ N_input_temp(:,:,:)/one_year * ratio_nh4_fert
             N_input_raw(:,:,:,ifert_nitr) = N_input_raw(:,:,:,ifert_nitr)+ N_input_temp(:,:,:)/one_year * (1.-ratio_nh4_fert)
 
             fieldname='Nfert_cropC4'
             CALL slowproc_Ninput(kjpindex, lalo, neighbours, resolution, contfrac, fieldname, &
                  N_input_temp, Ninput_year)
             ! Conversion from gN/m2(cropland)/yr to gN/m2/day
             N_input_raw(:,:,:,ifert_ammo) = N_input_raw(:,:,:,ifert_ammo)+ N_input_temp(:,:,:)/one_year * ratio_nh4_fert
             N_input_raw(:,:,:,ifert_nitr) = N_input_raw(:,:,:,ifert_nitr)+ N_input_temp(:,:,:)/one_year * (1.-ratio_nh4_fert)

             fieldname='Nfert_pasture'
             CALL slowproc_Ninput(kjpindex, lalo, neighbours, resolution, contfrac, fieldname, &
                  N_input_temp, Ninput_year)
             ! Conversion from gN/m2(pasture)/yr to gN/m2/day
             N_input_raw(:,:,:,ifert_ammo) = N_input_raw(:,:,:,ifert_ammo)+ N_input_temp(:,:,:)/one_year * ratio_nh4_fert
             N_input_raw(:,:,:,ifert_nitr) = N_input_raw(:,:,:,ifert_nitr)+ N_input_temp(:,:,:)/one_year * (1.-ratio_nh4_fert)

             fieldname='Nfert_ammo_pasture'
             CALL slowproc_Ninput(kjpindex, lalo, neighbours, resolution, contfrac, fieldname, &
                  N_input_temp, Ninput_year)
             ! Conversion from gN/m2(pasture)/yr to gN/m2/day
             N_input_raw(:,:,:,ifert_ammo) = N_input_raw(:,:,:,ifert_ammo)+ N_input_temp(:,:,:)/one_year

             fieldname='Nfert_nitr_pasture'
             CALL slowproc_Ninput(kjpindex, lalo, neighbours, resolution, contfrac, fieldname, &
                  N_input_temp, Ninput_year)
             ! Conversion from gN/m2(pasture)/yr to gN/m2/day
             N_input_raw(:,:,:,ifert_nitr) = N_input_raw(:,:,:,ifert_nitr)+ N_input_temp(:,:,:)/one_year

          ENDIF
       ELSE
             !Config Key   = NFERT
             !Config Desc  = Amount of N fertiliser 
             !Config Def   = 0
             !Config If    = ok_ncycle .AND. (.NOT. impose_cn)
             !Config Help  = 
             !Config Units = [gN m-2 d-1] 
             nfert=zero
             CALL getin_p ('NFERT',nfert)       
             n_input_raw(:,:,:,ifert_ammo)=nfert * ratio_nh4_fert
             n_input_raw(:,:,:,ifert_nitr)=nfert * (1.-ratio_nh4_fert)
       ENDIF


       IF(.NOT. impose_Ninput_manure) THEN
          FOUND_RESTART=.TRUE.

          CALL ioconf_setatt_p('UNITS', 'kgN m-2 yr-1')
          CALL ioconf_setatt_p('LONG_NAME','N manure')
          CALL restget_p (rest_id, 'Nmanure', nbp_glo, nvm, 12, kjit, .TRUE., N_input_raw(:,:,:,imanure), "gather", nbp_glo, index_g)
          IF ( ALL( N_input_raw(:,:,:,imanure) .EQ. val_exp ) ) FOUND_RESTART=.FALSE.


          IF(.NOT. FOUND_RESTART) THEN
             ! Read the new N inputs from file. Output is Ninput and frac_nobio_nextyear.
             N_input_raw(:,:,:,imanure) = zero
             fieldname='Nmanure_cropland'
             CALL slowproc_Ninput(kjpindex, lalo, neighbours, resolution, contfrac, fieldname, &
                  N_input_temp, Ninput_year)
             ! Conversion from gN/m2(cropland)/yr to gN/m2/day
             N_input_raw(:,:,:,imanure) = N_input_raw(:,:,:,imanure)+N_input_temp(:,:,:)/one_year
                
             fieldname='Nmanure_pasture'
             CALL slowproc_Ninput(kjpindex, lalo, neighbours, resolution, contfrac, fieldname, &
                  N_input_temp, Ninput_year)
             ! Conversion from gN/m2(cropland)/yr to gN/m2/day
             N_input_raw(:,:,:,imanure) = N_input_raw(:,:,:,imanure)+N_input_temp(:,:,:)/one_year
          ENDIF
       ELSE
          !Config Key   = NMANURE
          !Config Desc  = Amount of N manure 
          !Config Def   = 0
          !Config If    = ok_ncycle .AND. (.NOT. impose_cn)
          !Config Help  = 
          !Config Units = [gN m-2 d-1] 
          nmanure=zero
          CALL getin_p ('NMANURE',nmanure)       
          n_input_raw(:,:,:,imanure)=nmanure
       ENDIF



       IF(.NOT. impose_Ninput_bnf) THEN
          FOUND_RESTART=.TRUE.
          CALL ioconf_setatt_p('UNITS', 'kgN m-2 yr-1')
          CALL ioconf_setatt_p('LONG_NAME','N bilogical fixation')
          CALL restget_p (rest_id, 'Nbnf', nbp_glo, nvm, 12, kjit, .TRUE., N_input_raw(:,:,:,ibnf), "gather", nbp_glo, index_g)
          IF ( ALL( N_input_raw(:,:,:,ibnf) .EQ. val_exp ) ) FOUND_RESTART=.FALSE.
          
          IF(.NOT. FOUND_RESTART) THEN
             fieldname='Nbnf'
             CALL slowproc_Ninput(kjpindex, lalo, neighbours, resolution, contfrac, fieldname, &
                  N_input_raw(:,:,:,ibnf), Ninput_year)

             ! Conversion from kgN/km2/yr to gN/m2/day
             N_input_raw(:,:,:,ibnf) = N_input_raw(:,:,:,ibnf)/1000./one_year
          ENDIF
       ELSE
          !Config Key   = NBNF
          !Config Desc  = Amount of N biological fixation 
          !Config Def   = 0
          !Config If    = ok_ncycle .AND. (.NOT. impose_cn)
          !Config Help  = 
          !Config Units = [gN m-2 d-1] 
          nbnf=zero
          CALL getin_p ('NBNF',nbnf)       
          n_input_raw(:,:,:,ibnf)=nbnf
       ENDIF
    ELSE
       n_input_raw(:,:,:,:)=zero
    ENDIF

    !! Calculate fraction of landuse tiles to be used only for diagnostic variables
    fraclut(:,:)=0
    nwdFraclut(:,id_psl)=0
    nwdFraclut(:,id_crp)=1.
    nwdFraclut(:,id_urb)=xios_default_val
    nwdFraclut(:,id_pst)=xios_default_val
    DO jv=1,nvm
       IF (natural(jv)) THEN
          fraclut(:,id_psl) = fraclut(:,id_psl) + veget_max(:,jv)
          IF(.NOT. is_tree(jv)) THEN
             nwdFraclut(:,id_psl) = nwdFraclut(:,id_psl) + veget_max(:,jv) 
          ENDIF
       ELSE
          fraclut(:,id_crp) = fraclut(:,id_crp) + veget_max(:,jv)
       ENDIF
    END DO
    
    WHERE (fraclut(:,id_psl) > min_sechiba)
       nwdFraclut(:,id_psl) = nwdFraclut(:,id_psl)/fraclut(:,id_psl)
    ELSEWHERE
       nwdFraclut(:,id_psl) = xios_default_val
    END WHERE    


    IF (printlev_loc>=3) WRITE (numout,*) ' slowproc_init done '
    
  END SUBROUTINE slowproc_init

!! ================================================================================================================================
!! SUBROUTINE   : slowproc_clear
!!
!>\BRIEF          Clear all variables related to slowproc and stomate modules  
!!
!_ ================================================================================================================================

  SUBROUTINE slowproc_clear 

  ! 1 clear all the variables defined as common for the routines in slowproc 

    IF (ALLOCATED (clayfraction)) DEALLOCATE (clayfraction)
    IF (ALLOCATED (sandfraction)) DEALLOCATE (sandfraction)
    IF (ALLOCATED (siltfraction)) DEALLOCATE (siltfraction)
    IF (ALLOCATED (bulk)) DEALLOCATE (bulk)
    IF (ALLOCATED (soil_ph)) DEALLOCATE (soil_ph)
    IF (ALLOCATED (cc_biomass_m)) DEALLOCATE (cc_biomass_m)
    IF (ALLOCATED (cc_n_m)) DEALLOCATE (cc_n_m)
    IF (ALLOCATED (irrigated_new)) DEALLOCATE (irrigated_new)

 ! 2. Clear all the variables in stomate 

    CALL stomate_clear 
    !
  END SUBROUTINE slowproc_clear

!!$!! ================================================================================================================================
!!$!! SUBROUTINE   : slowproc_derivvar
!!$!!
!!$!>\BRIEF         Initializes variables related to the
!!$!! parameters to be assimilated, the maximum water on vegetation, the vegetation height, 
!!$!! and the fraction of soil covered by dead leaves and the vegetation height 
!!$!!
!!$!! DESCRIPTION  : (definitions, functional, design, flags):
!!$!! (1) Initialization of the variables relevant for the assimilation parameters  
!!$!! (2) Intialization of the fraction of soil covered by dead leaves
!!$!! (3) Initialization of the Vegetation height per PFT
!!$!! (3) Initialization the maximum water on vegetation for interception with a particular treatement of the PFT no.1
!!$!!
!!$!! RECENT CHANGE(S): None
!!$!!
!!$!! MAIN OUTPUT VARIABLE(S): ::qsintmax, ::deadleaf_cover, ::height  
!!$!!
!!$!! REFERENCE(S) : None
!!$!!
!!$!! FLOWCHART    : None
!!$!! \n
!!$!_ ================================================================================================================================
!!$
!!$  SUBROUTINE slowproc_derivvar (kjpindex, veget, circ_class_biomass,circ_class_n, &
!!$       qsintmax, deadleaf_cover, height, temp_growth)
!!$
!!$    !! INTERFACE DESCRIPTION
!!$
!!$    !! 0.1 Input scalar and fields 
!!$    INTEGER(i_std),INTENT (in)                                  :: kjpindex       !! Domain size - terrestrial pixels only
!!$    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)           :: veget          !! Fraction of pixel covered by PFT in the mesh (unitless)
!!$    !! 0.2. Output scalar and fields 
!!$    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)          :: qsintmax       !! Maximum water on vegetation for interception(mm)
!!$    REAL(r_std),DIMENSION (kjpindex), INTENT (out)              :: deadleaf_cover !! fraction of soil covered by dead leaves (unitless)
!!$    REAL(r_std),DIMENSION (:,:,:,:,:), INTENT (inout)           :: circ_class_biomass
!!$    REAL(r_std),DIMENSION (:,:,:), INTENT (inout)               :: circ_class_n
!!$    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)          :: height         !! height of the vegetation or surface in general ??? (m)
!!$    REAL(r_std),DIMENSION (kjpindex), INTENT (out)              :: temp_growth    !! growth temperature (°C)  
!!$    !
!!$    !! 0.3 Local declaration
!!$    REAL(r_std),DIMENSION (kjpindex,nvm)                        :: lai            !! PFT leaf area index (m^{2} m^{-2})
!!$    INTEGER(i_std)                                              :: ji, jv         !! Local indices
!!$!_ ================================================================================================================================
!!$
!!$    !! 1. Intialize the fraction of soil covered by dead leaves 
!!$    deadleaf_cover(:) = zero
!!$
!!$    !! 2. Initialize vegetation height and LAI per PFT
!!$    height(:,1) = zero
!!$    lai(:,ibare_sechiba) = zero
!!$    DO jv = 1, nvm
!!$       DO ji = 1, kjpindex
!!$          height(:,jv) = wood_to_qmheight(circ_class_biomass(ji,jv,:,:,icarbon, &
!!$               circ_class_n(ji,jv,:),jv,pipe_tune2(ipts,ivm))
!!$          lai(ji,jv) = cc_to_lai(circ_class_biomass(ji,jv,:,ileaf,icarbon),&
!!$               circ_class_n(ji,jv,:),jv)
!!$       ENDDO
!!$    ENDDO
!!$    
!!$    !! 3. Initialize the maximum water on vegetation for interception
!!$    qsintmax(:,:) = qsintcst * veget(:,:) * lai(:,:) 
!!$    qsintmax(:,1) = zero
!!$
!!$    !! 4. Initialize the growth temperature
!!$    temp_growth(:)=25.
!!$
!!$  END SUBROUTINE slowproc_derivvar


!! ================================================================================================================================
!! SUBROUTINE   : slowproc_mean
!!
!>\BRIEF          Accumulates field_in over a period of dt_tot.
!! Has to be called at every time step (dt). 
!! Mean value is calculated if ldmean=.TRUE.
!! field_mean must be initialized outside of this routine! 
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!! (1) AcumAcuumlm 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::field_main
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_mean (kjpindex, n_dim2, dt_tot, dt, ldmean, field_in, field_mean)

    !
    !! 0 declarations

    !! 0.1 input scalar and variables 
    INTEGER(i_std), INTENT(in)                           :: kjpindex     !! Domain size- terrestrial pixels only 
    INTEGER(i_std), INTENT(in)                           :: n_dim2   !! Number of PFTs 
    REAL(r_std), INTENT(in)                              :: dt_tot   !! Time step of stomate (in days). The period over which the accumulation or the mean is computed 
    REAL(r_std), INTENT(in)                              :: dt       !! Time step in days 
    LOGICAL, INTENT(in)                                  :: ldmean   !! Flag to calculate the mean after the accumulation ???
    REAL(r_std), DIMENSION(kjpindex,n_dim2), INTENT(in)  :: field_in !! Daily field 

    !! 0.3 Modified field; The computed sum or mean field over dt_tot time period depending on the flag ldmean 
    REAL(r_std), DIMENSION(kjpindex,n_dim2), INTENT(inout)   :: field_mean !! Accumulated field at dt_tot time period or mean field over dt_tot 
 

!_ ================================================================================================================================

    !
    ! 1. Accumulation the field over dt_tot period 
    !
    field_mean(:,:) = field_mean(:,:) + field_in(:,:) * dt

    !
    ! 2. If the flag ldmean set, the mean field is computed over dt_tot period  
    !
    IF (ldmean) THEN
       field_mean(:,:) = field_mean(:,:) / dt_tot
    ENDIF

  END SUBROUTINE slowproc_mean


  
!! ================================================================================================================================
!! SUBROUTINE   : slowproc_long
!!
!>\BRIEF        Calculates a temporally smoothed field (field_long) from
!! instantaneous input fields.Time constant tau determines the strength of the smoothing.
!! For tau -> infinity??, field_long becomes the true mean value of field_inst
!! (but  the spinup becomes infinietly long, too).
!! field_long must be initialized outside of this routine! 
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!! (1) Testing the time coherence betwen the time step dt and the time tau over which
!! the rescaled of the mean is performed   
!!  (2) Computing the rescaled mean over tau period 
!! MAIN OUTPUT VARIABLE(S): field_long  
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::field_long
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_long (kjpindex, n_dim2, dt, tau, field_inst, field_long)

    !
    ! 0 declarations
    !

    ! 0.1 input scalar and fields 

    INTEGER(i_std), INTENT(in)                                 :: kjpindex        !! Domain size- terrestrial pixels only
    INTEGER(i_std), INTENT(in)                                 :: n_dim2      !! Second dimension of the fields, which represents the number of PFTs
    REAL(r_std), INTENT(in)                                    :: dt          !! Time step in days   
    REAL(r_std), INTENT(in)                                    :: tau         !! Integration time constant (has to have same unit as dt!)  
    REAL(r_std), DIMENSION(kjpindex,n_dim2), INTENT(in)            :: field_inst  !! Instantaneous field 


    ! 0.2 modified field

    ! Long-term field
    REAL(r_std), DIMENSION(kjpindex,n_dim2), INTENT(inout)         :: field_long  !! Mean value of the instantaneous field rescaled at tau time period 

!_ ================================================================================================================================

    !
    ! 1 test coherence of the time 

    IF ( ( tau .LT. dt ) .OR. ( dt .LE. zero ) .OR. ( tau .LE. zero ) ) THEN
       WRITE(numout,*) 'slowproc_long: Problem with time steps'
       WRITE(numout,*) 'dt=',dt
       WRITE(numout,*) 'tau=',tau
    ENDIF

    !
    ! 2 integration of the field over tau 

    field_long(:,:) = ( field_inst(:,:)*dt + field_long(:,:)*(tau-dt) ) / tau

  END SUBROUTINE slowproc_long


!! ================================================================================================================================
!! SUBROUTINE   : slowproc_canopy
!!
!>\BRIEF         Convert circ_class_biomass and circ_class_n in a 3D canopy used
!!               in the albedo, transpiration and energy budget calculations   
!!
!! DESCRIPTION  : 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::lai_per_level, ::z_array_out, 
!!                          ::max_height_store, ::laieff_fit, ::frac_age
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_canopy(kjpindex, circ_class_biomass, circ_class_n, &
       veget_max, lai_per_level, z_array_out, &
       max_height_store, laieff_fit, frac_age)
    !
    ! 0. Declarations
    !
    !! 0.1 Input variables 
    INTEGER(i_std), INTENT(in)                             :: kjpindex             !! Domain size - terrestrial pixels only
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(in)          :: circ_class_biomass   !! Biomass of the different components per
                                                                                   !! PFT and circ class (g C(N) tree-1 y-1)
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)              :: circ_class_n         !! Number of trees per circ_class (trees m-2)
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(in)       :: veget_max            !! Maximum fraction of vegetation type including
                                                                                   !! none biological fraction (unitless)

    !! 0.2 Modified variables 

    !! 0.3 Output
    REAL(r_std), DIMENSION(:,:,:), INTENT(out)             :: lai_per_level        !! This is the LAI per vertical level
                                                                                   !! @tex $(m^{2} m^{-2})$   
    REAL(r_std),DIMENSION(kjpindex,nvm,ncirc,nlevels_tot), INTENT(out) &
                                                           :: z_array_out          !! Height above soil of the Pgap points.
                                                                                   !! @tex $(m)$ @endtex
    REAL(r_std), DIMENSION(kjpindex, nvm), INTENT(out)     :: max_height_store     !! ???
    TYPE(laieff_type),DIMENSION (:,:,:),INTENT(out)        :: laieff_fit           !! Fitted parameters for the effective LAI
    REAL(r_std),DIMENSION(kjpindex,nvm,nleafages),INTENT(out) &
                                                           :: frac_age             !! leaf age distribution calculated in stomate 
    
    !! 0.4 Local
    REAL(r_std), DIMENSION(kjpindex,nvm,nlevels_tot)       :: z_level_photo        !! The height of the levels that we will
                                                                                   !! use to calculate the effective LAI for
                                                                                   !! the albedo routines and photosynthesis.
                                                                                   !! @tex $(m)$ @endtex
    
    INTEGER(i_std)                                         :: ji, jv, icir         !! Indices    

!_ ================================================================================================================================

    !! 1. Calculate canopy structure
    ! Compute the height of the photosynthesis levels given the height 
    ! of the energy levels and the vegetation on the grid square. The 
    ! hightest levels will be a function of the height of the vegetation 
    ! so that we don't waste computational time on empty levels.
    CALL calculate_z_level_photo(kjpindex, circ_class_biomass, circ_class_n, &
         z_level_photo)
      
    ! Finding the true LAI per level is different from finding the
    ! effective LAI per level, so we'll do that here.  This only
    ! changes once every day so hopefully it is not too expensive.
    ! It will eventually be needed by the energy budget.  It is
    ! also needed by effective_lai for grasses and crops.
    CALL find_lai_per_level(kjpindex, z_level_photo, &
         circ_class_biomass, circ_class_n, lai_per_level, &
         max_height_store)
    
    ! Change the dimensions of z_level_photo so it can be used in 
    ! sechiba, mleb, and fitting_laieff
    DO icir = 1,ncirc
       z_array_out(:,:,icir,:) = z_level_photo(:,:,:)
    END DO

    ! Now we actually find the effective LAI and fit the function 
    ! we'll use later on. Note that this fir allows us to calculate 
    ! the effective lai once per day. This should not be repeated
    ! every half hour!
    CALL fitting_laieff(kjpindex, z_array_out, circ_class_biomass, & 
         circ_class_n, veget_max, lai_per_level, laieff_fit)
    

    !! 2. Calculate frac_age
    ! The variable frac_age is only used when BVOCs are calculated.
    ! slowproc_canopy should only be used when only sechiba is used.
    frac_age(:,:,1) = un
    frac_age(:,:,2) = zero
    frac_age(:,:,3) = zero
    frac_age(:,:,4) = zero

    IF (printlev_loc>=3) WRITE(numout,*) 'Leaving slowproc_canopy'

  END SUBROUTINE slowproc_canopy


!! ================================================================================================================================
!! SUBROUTINE   : slowproc_veget
!!
!>\BRIEF        Calucate veget and soiltile.
!!
!! DESCRIPTION  : Calucate veget and soiltile.
!! (1) Calculate veget
!! (2) Calculate totfrac_nobio
!! (3) Calculate soiltile
!! (4) Calculate fraclut
!!
!! RECENT CHANGE(S):
!!  Slowproc_veget is called four times in slowproc: (1) after reading and cleaning the 
!!  land cover map in slowproc. There is no reason to do these tests again as they 
!!  were just done in the code prior to this call. (2) during the initialization. 
!!  veget_max should is either read from a file (and thus error checked) or comes 
!!  from a restart file and was error checked earlier. (3) After the call to stomate
!!  main. In stomate veget_max may change after LCC. Sapiens_lcchange avoids too
!!  small changes. No need to check again. (4) slowproc_change_frac calls slowproc_veget.
!!  slowproc_change_frac is called once in sechiba_main after slowproc_main (and thus
!!  indirectly after sapiens_lcchange.f90). The fraction should have been taken care of
!!  sapiens_lcchange.f90. 
!!  Keeping that code is a perfect way to create undetectable mass balance problems 
!!  because: (1) mass balance closure is not checked in slowproc and sechiba yet and 
!!  (2) these PFTs contain carbon and nitrogen. Simply truncating and normalizing
!!  their cover fractions should result in mass balance problems. By removing the
!!  code below, changes in veget_max occur right after reading the maps (and thus
!!  before the PFT is getting a biomass) or in sapiens_lcchange where changes in
!!  cover fractions are correctly accounted for and mass balance closure is checked.
!!  In the trunk slowproc_veget_max_limit contains the code that has been commented out
!!  in this version. When merging do not accept this code or add a mass balance check
!!  over slowproc and sechiba.
!!
!! MAIN OUTPUT VARIABLE(S): :: frac_nobio, totfrac_nobio, veget_max, veget, soiltile, fraclut
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_veget (kjpindex, lai_per_level, z_array_out, &
       coszang_noon, circ_class_biomass, circ_class_n, frac_nobio, totfrac_nobio, &
       veget_max, veget, soiltile, tot_bare_soil, fraclut, nwdFraclut, Pgap_cumul,&
       precip_longterm)
    
    !! 0.1 Input variables 
    INTEGER(i_std), INTENT(in)                             :: kjpindex             !! Domain size - terrestrial pixels only
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)              :: lai_per_level        !! This is the LAI per vertical level
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(in)          :: circ_class_biomass   !! Biomass of the different components per
                                                                                   !! PFT and circ class (g C(N) tree-1 y-1)
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)              :: circ_class_n         !! Number of trees per circ_class (trees m-2)                                                                                                               !! @tex $(m^{2} m^{-2})$ 
    REAL(r_std), DIMENSION(:,:), INTENT(in)                :: frac_nobio           !! Fraction of the mesh which is covered by ice, lakes, ...
    REAL(r_std), DIMENSION(:,:), INTENT(in)                :: veget_max            !! Maximum fraction of vegetation type including
                                                                                   !! none biological fraction (unitless)
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)            :: z_array_out          !! The physical height of the levels used in the
                                                                                   !! effective LAI routines @tex $(m)$
    REAL(r_std), DIMENSION(kjpindex), INTENT(in)           :: coszang_noon         !! The cosine of the solar zenith angle at noon
    REAL(r_std), DIMENSION(kjpindex), INTENT(in)           :: precip_longterm      !! Longterm annual precipitation sum mm year^{-1}
 
    
    !! 0.2 Output variables 
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(out)      :: veget                !! Fraction of pixel covered by PFT in the mesh (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: totfrac_nobio
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT(out)    :: soiltile             !! Fraction of each soil tile within vegtot (0-1, unitless)
    REAL(r_std), DIMENSION (kjpindex,nlut), INTENT(out)    :: fraclut              !! Fraction of each landuse tile (0-1, unitless)
    REAL(r_std), DIMENSION (kjpindex,nlut), INTENT(out)    :: nwdFraclut           !! Fraction of non-woody vegetation in each landuse tile (0-1, unitless)
    REAL(r_std), DIMENSION(kjpindex,nvm,nlevels_tot), INTENT(out) :: Pgap_cumul    !! The probability of finding a gap in the in canopy from the top 
                                                                                   !! of the canopy to a given level (unitless, between 0-1)

    !! 0.3 Modified variables   
    REAL(r_std),DIMENSION(:), INTENT (inout)               :: tot_bare_soil        !! Total evaporating bare soil fraction 

    !! 0.4 Local scalar and varaiables 
    INTEGER(i_std)                                         :: ilev, ji, jv         !! Indices
    INTEGER(i_std)                                         :: jst, ivm             !! Indices
    REAL(r_std), DIMENSION(nlevels_tot,kjpindex,nvm)       :: laieff_temp          !! The effective LAI.  Not used here, just need it as an
                                                                                   !! argument in this routine.
    REAL(r_std), DIMENSION(kjpindex,nvm,1)                 :: lai_per_level_temp   !! The total LAI found for the PFT and grid square.
    REAL(r_std),DIMENSION(kjpindex,nvm,nlevels_tot)        :: z_array4             !! Same as z_array, but one less dimension.
                                                                                   !! @tex $(m)$ @endtex (local because the call to
                                                                                   !! effective_lai requires it) 
    REAL(r_std), DIMENSION(kjpindex,nvm,nlevels_tot)       :: z_array2_out         !! Same as z_array, but one less dimension.
                                                                                   !! @tex $(m)$ @endtex (local because the 
                                                                                   !! call to effective_lai requires it)]
    REAL(r_std), DIMENSION(kjpindex)                       :: fracsum              !! Sum of both fracnobio and veget_max
    REAL(r_std), DIMENSION(kjpindex)                       :: fedge_pgap           !! Edge effect weight on canopy porosity (0-1)
                                                                                   !! See MODULE_DESIGN_AED_LIGHT.md §4.2

!_ ================================================================================================================================

    IF(printlev_loc>=4) WRITE(numout,*) 'Entering slowproc_veget'


    !! 4. Calculate veget making use of the canopy structure
    ! Determine the coverage based on stand structure and LAI. There
    ! should always be trees in the model even if we only run sechiba.
    ! The veget should be related to the effective LAI if you are looking
    ! straight down on the canopy. So let us use the same routine
    ! that we use for the effective LAI but put the solar angle equal to
    ! the zenith, zero degrees. We use the light that reaches the forest
    ! floor in the calculation of veget. For some applications, i.e.,
    ! recruitment it would be better to use the solar angle at noon, but
    ! for turbulence-related and precipitation-related application zenith
    ! is preferred. To avoid having different Pgaps we use zenith.
    veget(:,:) = zero

    CALL effective_lai(kjpindex, nlevels_tot, z_array_out,  circ_class_biomass, &
         circ_class_n, veget_max, coszang_noon, lai_per_level, &
         laieff_temp, Pgap_cumul_out=Pgap_cumul)

    !! 4bis. Non-directional canopy edge light effect
    ! Guillaume M. -- Same modulation as in fitting_laieff (stomate_laieff), applied here so
    ! that the SAVEd Pgap_cumul feeding the veget -> condveg (albedo) -> enerbil cascade stays
    ! consistent. AED_light comes from stomate_data.calculate_aed, refreshed daily by stomate.
    ! See MODULE_DESIGN_AED_LIGHT.md §4.2.
    IF (ok_aed_light .AND. ALLOCATED(AED_light)) THEN
       fedge_pgap(:) = 1.0_r_std - ((AED_light(:) - edge_distance_light/2._r_std)**2) &
                                  / (AED_light(:)**2)
       DO ivm = 1, nvm
          DO ilev = 1, nlevels_tot
             Pgap_cumul(:,ivm,ilev) = Pgap_cumul(:,ivm,ilev) &
                                    * (1.0_r_std + pgap_edge_factor * fedge_pgap(:))
          ENDDO
       ENDDO
       Pgap_cumul(:,:,:) = MIN(1.0_r_std, MAX(0.0_r_std, Pgap_cumul(:,:,:)))
       CALL xios_orchidee_send_field("AED_LIGHT_SLOWPROC", AED_light)
       CALL xios_orchidee_send_field("FEDGE_PGAP_SLOWPROC", fedge_pgap)
    ENDIF

    ! Now convert the transmission probability (1, light makes it through
    ! without hitting anything) to veget (1*veget_max, the canopy is
    ! completely opaque).
    DO ivm=1,nvm

       ! We don't want to overwrite the bare soil value, since that
       ! means something different than the other values.
       IF(ivm == ibare_sechiba) THEN
          ! Calculate veget
          veget(:,ibare_sechiba)=zero
       ELSE
          IF (hack_pgap) THEN
             ! Note that the exctinction coefficient for calculating veget
             ! was set to 1.0. Use the old approach following Lambert-Beer
             veget(:,ivm) = veget_max(:,ivm)*(un-exp(-SUM(lai_per_level(:,ivm,:),2)*1.0))
          ELSE
             ! Use the new approach based on Pgap
             veget(:,ivm)=veget_max(:,ivm)*(un-Pgap_cumul(:,ivm,1))
          END IF
       ENDIF

    ENDDO

    !+++DBEN OUTPUT+++
    CALL xios_orchidee_send_field("CA_dben",(un-Pgap_cumul(:,:,1))*m2_to_ha)
    !+++++++++++++++++

    !! 6. Calculate totfrac_nobio and tot_bare_soil making use of veget 
    IF (ok_bare_soil_new) THEN
       
       !+++CHECK+++
       ! We no longer want to treat the gaps in the canopy as
       ! bare soil. It needs to be tested what will happen with
       ! the evaporation in the single-layer model. The multi-
       ! layer energy budget should be able to correctly deal
       ! with the gaps in the canopy.
       tot_bare_soil(:) = veget_max(:,1)
       
       ! Total frac nobio
       totfrac_nobio(:) = SUM(frac_nobio(:,:),2)
       !+++++++++++
       
    ELSE

       ! Initialize
       fracsum(:) = zero
       tot_bare_soil(:) = veget_max(:,1)

       ! Calculate bare soil fraction
       DO ji = 1, kjpindex

          ! Total frac nobio
          totfrac_nobio(ji) = SUM(frac_nobio(ji,:))
          
          DO jv = 2, nvm

             ! Move the canopy gaps into the bare soil
             ! fraction for this time step. Calculate
             ! the total fraction of bare soil in the
             ! grid
             tot_bare_soil(ji) = tot_bare_soil(ji) + &
                  (veget_max(ji,jv) - veget(ji,jv))
             fracsum(ji) = fracsum(ji) + veget(ji,jv)

          ENDDO

          ! Consistency check
          fracsum(ji) = fracsum(ji) + tot_bare_soil(ji) + SUM(frac_nobio(ji,:)) 
          IF (fracsum(ji) .LT. 0.99999) THEN
             WRITE(numout,*)' ATTENTION, in ji, fracsum LT 1: ', ji, fracsum(ji)
             WRITE(numout,*)'  frac_nobio = ',SUM(frac_nobio(ji,:))
             WRITE(numout,*)'  veget = ',SUM(veget(ji,:))
             WRITE(numout,*)'  veget_max = ',SUM(veget_max(ji,:))
             WRITE(numout,*)'  tot_bare_soil = ',tot_bare_soil(ji)
          ENDIF
       ENDDO
       
    ENDIF ! ok_bare_soil_new
    
    !! 4. Calculate soiltiles
    !  Soiltiles are only used in hydrol, but we fix them in here because some time 
    !  it might depend on a changing vegetation (but then some adaptation should be 
    ! made to hydrol) and be also used in the other modules to perform separated 
    ! energy balances. The sum of all soiltiles makes one, and corresponds to the 
    ! bio fraction of the grid cell (called vegtot in hydrol).   
    soiltile(:,:) = zero
    DO jv = 1, nvm
       jst = pref_soil_veg(jv)
       DO ji = 1, kjpindex
          soiltile(ji,jst) = soiltile(ji,jst) + veget_max(ji,jv)
       ENDDO
    ENDDO

    DO ji = 1, kjpindex 
       IF (totfrac_nobio(ji) .LT. (1-EPSILON(un))) THEN
          ! If the calculation of 1-totfrac_nobio is correct its value
          ! should be identical to SUM(soiltile) give or take precision
          ! issues.
          IF (ABS(SUM(soiltile(ji,:))-(1-totfrac_nobio(ji))).LT.min_sechiba) THEN
             ! If the numbers are very similar, the divsion may result in 
             ! errors of 10e-10. By taking the minimum we protect against
             ! such conditions. Note that the new notation is expected
             ! to be a bit more precise than the old notation. It probably
             ! doesn really matter but it may help to increase surface
             ! balance closure at 10e-13 or better.
!!$             soiltile(ji,:)=MIN(soiltile(ji,:)/(1.-totfrac_nobio(ji)),un)
             soiltile(ji,:)=MIN(soiltile(ji,:)/SUM(soiltile(ji,:)),un)
          ELSE
             ! Mismatch between SUM(soiltile) and 1-totfrac_nobio.
             WRITE(numout,*) 'kjpindex, totfrac_nobio, soiltile, ', &
                  ji, totfrac_nobio(ji), SUM(soiltile(ji,:))  
             CALL ipslerr_p (3,'slowproc_main',&
                  'Mismatch between SUM(soiltile) and 1-totfrac_nobio','','')
          ENDIF
       ELSE
          soiltile(ji,:)=zero
       ENDIF
    ENDDO  

    !! 5. Calculate fraction of landuse tiles to be used only for diagnostic variables
    fraclut(:,:)=0
    nwdFraclut(:,id_psl)=0
    nwdFraclut(:,id_crp)=1.
    nwdFraclut(:,id_urb)=xios_default_val
    nwdFraclut(:,id_pst)=xios_default_val
    DO jv=1,nvm
       IF (natural(jv)) THEN
          fraclut(:,id_psl) = fraclut(:,id_psl) + veget_max(:,jv)
          IF(.NOT. is_tree(jv)) THEN
             nwdFraclut(:,id_psl) = nwdFraclut(:,id_psl) + veget_max(:,jv) 
          ENDIF
       ELSE
          fraclut(:,id_crp) = fraclut(:,id_crp) + veget_max(:,jv)
       ENDIF
    END DO
    

    WHERE (fraclut(:,id_psl) > min_sechiba)
       nwdFraclut(:,id_psl) = nwdFraclut(:,id_psl)/fraclut(:,id_psl)
    ELSEWHERE
       nwdFraclut(:,id_psl) = xios_default_val
    END WHERE    

  END SUBROUTINE slowproc_veget


!! ================================================================================================================================
!! SUBROUTINE   : slowproc_interlai

    ! Read an LAI map. In ORCHIDEE-CN-CAN lai is no longer passed from one 
    ! routine to another but more importantly it changed from a 1-D to a 3-D
    ! variable. Reading, e.g. MODIS LAI, and using that information to force 
    ! ORCHIDEE-CN-CAN would first require to downscale the observed pixel-level
    ! LAI into PFT-level LAI. Subsequently, assumptions should be made to 
    ! convert the 1-D LAI into circ_class_biomass and circ_class_n which in 
    ! turn are used to calculate a 3-D canopy structure that is used in the
    ! laieff calculation to calculate a 1-D effective LAI. The 1-D MODIS LAI
    ! cannot be assumed to be the same as the 1-D effective LAI. What was
    ! previously known as read_lai and laimap has therefore been replaces by
    ! code that reads a file containing circ_class_biomass and circ_class_n.

  SUBROUTINE slowproc_interlai(nbpt, lalo, resolution, neighbours, contfrac, laimap)

    USE interpweight

    IMPLICIT NONE

    !
    !
    !
    !  0.1 INPUT
    !
    INTEGER(i_std), INTENT(in)          :: nbpt                  !! Number of points for which the data needs to be interpolated
    REAL(r_std), INTENT(in)             :: lalo(nbpt,2)          !! Vector of latitude and longitudes 
                                                                 !! (beware of the order = 1 : latitude, 2 : longitude)
    REAL(r_std), INTENT(in)             :: resolution(nbpt,2)    !! The size in km of each grid-box in X and Y
    INTEGER(i_std), INTENT(in)          :: neighbours(nbpt,NbNeighb)!! Vector of neighbours for each grid point
                                                                 !! (1=North and then clockwise)
    REAL(r_std), INTENT(in)             :: contfrac(nbpt)        !! Fraction of land in each grid box.
    !
    !  0.2 OUTPUT
    !
    REAL(r_std), INTENT(out)    ::  laimap(nbpt,nvm,12)          !! lai read variable and re-dimensioned
    !
    !  0.3 LOCAL
    !
    CHARACTER(LEN=80) :: filename                               !! name of the LAI map read
    INTEGER(i_std) :: ib, ip, jp, it, jv
    REAL(r_std) :: lmax, lmin, ldelta
    LOGICAL ::           renormelize_lai  ! flag to force LAI renormelization
    INTEGER                  :: ier

    REAL(r_std), DIMENSION(nbpt)                         :: alaimap          !! availability of the lai interpolation 
    INTEGER, DIMENSION(4)                                :: invardims
    REAL(r_std), DIMENSION(nbpt,nvm,12)                  :: lairefrac        !! lai fractions re-dimensioned
    REAL(r_std), DIMENSION(nbpt,nvm,12)                  :: fraclaiinterp    !! lai fractions re-dimensioned
    REAL(r_std), DIMENSION(:), ALLOCATABLE               :: vmin, vmax       !! min/max values to use for the 
                                                                             !!   renormalization
    CHARACTER(LEN=80)                                    :: variablename     !! Variable to interpolate
    CHARACTER(LEN=80)                                    :: lonname, latname !! lon, lat names in input file
    REAL(r_std), DIMENSION(nvm)                          :: variabletypevals !! Values for all the types of the variable
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
    REAL(r_std), DIMENSION(3)                            :: maskvals         !! values to use to mask (according to 
                                                                             !!   `maskingtype') 
    CHARACTER(LEN=250)                                   :: namemaskvar      !! name of the variable to use to mask 
!_ ================================================================================================================================

    !
    !Config Key   = LAI_FILE
    !Config Desc  = Name of file from which the vegetation map is to be read
    !Config If    = LAI_MAP
    !Config Def   = lai2D.nc
    !Config Help  = The name of the file to be opened to read the LAI
    !Config         map is to be given here. Usualy SECHIBA runs with a 5kmx5km
    !Config         map which is derived from a Nicolas VIOVY one. 
    !Config Units = [FILE]
    !
    filename = 'lai2D.nc'
    CALL getin_p('LAI_FILE',filename)
    variablename = 'LAI'

    IF (xios_interpolation) THEN
       IF (printlev_loc >= 1) WRITE(numout,*) "slowproc_interlai: Use XIOS to read and interpolate " &
            // TRIM(filename) //" for variable " //TRIM(variablename)
    
       CALL xios_orchidee_recv_field('lai_interp',lairefrac)
       CALL xios_orchidee_recv_field('frac_lai_interp',fraclaiinterp)      
       alaimap(:) = fraclaiinterp(:,1,1)
    ELSE

      IF (printlev_loc >= 2) WRITE(numout,*) "slowproc_interlai: Start interpolate " &
           // TRIM(filename) //" for variable " //TRIM(variablename)

      ! invardims: shape of variable in input file to interpolate
      invardims = interpweight_get_var4dims_file(filename, variablename)
      ! Check coherence of dimensions read from the file
      IF (invardims(4) /= 12)  CALL ipslerr_p(3,'slowproc_interlai','Wrong dimension of time dimension in input file for lai','','')
      IF (invardims(3) /= nvm) CALL ipslerr_p(3,'slowproc_interlai','Wrong dimension of PFT dimension in input file for lai','','')

      ALLOCATE(vmin(nvm),stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3,'slowproc_interlai','Problem in allocation of variable vmin','','')

      ALLOCATE(vmax(nvm), STAT=ier)
      IF (ier /= 0) CALL ipslerr_p(3,'slowproc_interlai','Problem in allocation of variable vmax','','')


! Assigning values to vmin, vmax
      vmin = un
      vmax = nvm*un

      variabletypevals = -un

      !! Variables for interpweight
      ! Type of calculation of cell fractions
      fractype = 'default'
      ! Name of the longitude and latitude in the input file
      lonname = 'longitude'
      latname = 'latitude'
      ! Should negative values be set to zero from input file?
      nonegative = .TRUE.
      ! Type of mask to apply to the input data (see header for more details)
      maskingtype = 'mbelow'
      ! Values to use for the masking
      maskvals = (/ 20., undef_sechiba, undef_sechiba /)
      ! Name of the variable with the values for the mask in the input file (only if maskkingtype='var') (here not used)
      namemaskvar = ''

      CALL interpweight_4D(nbpt, nvm, variabletypevals, lalo, resolution, neighbours,        &
        contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype,        &
        maskvals, namemaskvar, nvm, invardims(4), -1, fractype,                            &
        -1., -1., lairefrac, alaimap)

      IF (printlev_loc >= 5) WRITE(numout,*)'  slowproc_interlai after interpweight_4D'

    ENDIF



    !
    !
    !Config Key   = RENORM_LAI
    !Config Desc  = flag to force LAI renormelization
    !Config If    = LAI_MAP
    !Config Def   = n
    !Config Help  = If true, the laimap will be renormalize between llaimin and llaimax parameters.
    !Config Units = [FLAG]
    !
    renormelize_lai = .FALSE.
    CALL getin_p('RENORM_LAI',renormelize_lai)

    !
    laimap(:,:,:) = zero
    !
    IF (printlev_loc >= 5) THEN
      WRITE(numout,*)'  slowproc_interlai before starting loop nbpt:', nbpt
    END IF 

    ! Assigning the right values and giving a value where information was not found
    DO ib=1,nbpt
      IF (alaimap(ib) < min_sechiba) THEN
        DO jv=1,nvm
          laimap(ib,jv,:) = (llaimax(jv)+llaimin(jv))/deux
        ENDDO
      ELSE
        DO jv=1, nvm
          DO it=1, 12
            laimap(ib,jv,it) = lairefrac(ib,jv,it)
          ENDDO
        ENDDO
      END IF
    ENDDO
    !
    ! Normelize the read LAI by the values SECHIBA is used to
    !
    IF ( renormelize_lai ) THEN
       DO ib=1,nbpt
          DO jv=1, nvm
             lmax = MAXVAL(laimap(ib,jv,:))
             lmin = MINVAL(laimap(ib,jv,:))
             ldelta = lmax-lmin
             IF ( ldelta < min_sechiba) THEN
                ! LAI constante ... keep it constant
                laimap(ib,jv,:) = (laimap(ib,jv,:)-lmin)+(llaimax(jv)+llaimin(jv))/deux
             ELSE
                laimap(ib,jv,:) = (laimap(ib,jv,:)-lmin)/(lmax-lmin)*(llaimax(jv)-llaimin(jv))+llaimin(jv)
             ENDIF
          ENDDO
       ENDDO
    ENDIF

    ! Write diagnostics
    CALL xios_orchidee_send_field("interp_avail_alaimap",alaimap)
    CALL xios_orchidee_send_field("interp_diag_lai",laimap)
   
    IF (printlev_loc >= 3) WRITE(numout,*) '  slowproc_interlai ended'

  END SUBROUTINE slowproc_interlai

!! ================================================================================================================================
!! SUBROUTINE   : slowproc_readvegetmax
!!
!>\BRIEF          Read and interpolate a vegetation map (by pft)
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): The subroutine was previously called slowproc_update.
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_readvegetmax(nbpt, lalo, neighbours,  resolution, contfrac, veget_last_in,         & 
       veget_next_out, frac_nobio_next, init)

    USE interpweight
    IMPLICIT NONE

    !
    !
    !
    !  0.1 INPUT
    !
    INTEGER(i_std), INTENT(in)                             :: nbpt            !! Number of points for which the data needs 
                                                                              !! to be interpolated
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in)             :: lalo            !! Vector of latitude and longitudes (beware of the order !)
    INTEGER(i_std), DIMENSION(nbpt,NbNeighb), INTENT(in)   :: neighbours      !! Vector of neighbours for each grid point
                                                                              !! (1=North and then clockwise)
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in)             :: resolution      !! The size in km of each grid-box in X and Y
    REAL(r_std), DIMENSION(nbpt), INTENT(in)               :: contfrac        !! Fraction of continent in the grid
    !
    REAL(r_std), DIMENSION(nbpt,nvm), INTENT(in)           :: veget_last_in   !! old max vegetfrac
    LOGICAL, INTENT(in)                                    :: init            !! initialisation : in case of dgvm, it forces update of all PFTs
    !
    !  0.2 OUTPUT
    !
    REAL(r_std), DIMENSION(nbpt,nvm), INTENT(out)          :: veget_next_out   !! new max vegetfrac
    REAL(r_std), DIMENSION(nbpt,nnobio), INTENT(out)       :: frac_nobio_next  !! new fraction of the mesh which is 
                                                                               !! covered by ice, lakes, ...
    
    !
    !  0.3 LOCAL
    !
    !
    CHARACTER(LEN=80) :: filename
    INTEGER(i_std) :: ib, inobio, jv, ivma, iagec
    REAL(r_std) :: sumf, err, norm
    !
    ! for DGVM case :
    REAL(r_std)                 :: sum_veg                     ! sum of vegets
    REAL(r_std)                 :: sum_nobio                   ! sum of nobios
    REAL(r_std)                 :: sumvAnthro_old, sumvAnthro  ! last an new sum of antrhopic vegets
    REAL(r_std)                 :: rapport                     ! (S-B) / (S-A)
    LOGICAL                     :: partial_update              ! if TRUE, partialy update PFT (only anthropic ones) 
                                                               ! e.g. in case of DGVM and not init (optional parameter)
    REAL(r_std), DIMENSION(nbpt,nvmap)                   :: veget_last       !! Temporary variable for veget_last_in on the same number of pfts as in the file 
    REAL(r_std), DIMENSION(nbpt,nvmap)                   :: veget_next       !! Temporary variable for veget_next_out on the same number of pfts as in the file
    REAL(r_std), DIMENSION(nbpt,nvmap)                   :: vegetrefrac      !! veget fractions re-dimensioned
    REAL(r_std), DIMENSION(nbpt,nvmap)                   :: aveget_nvmap     !! Availability of the interpolation
    REAL(r_std), DIMENSION(nbpt)                         :: aveget           !! Availability of the interpolation
    REAL(r_std), DIMENSION(nvmap)                        :: vmin, vmax       !! min/max values to use for the renormalization
    CHARACTER(LEN=80)                                    :: variablename     !! Variable to interpolate
    CHARACTER(LEN=80)                                    :: lonname, latname !! lon, lat names in input file
    REAL(r_std), DIMENSION(nvmap)                        :: variabletypevals !! Values for all the types of the variable
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
    REAL(r_std), DIMENSION(nbpt,nagec)                   :: acfrac           !! Age-class area fractions per pixel (structured cold start)
    INTEGER(i_std)                                       :: ic_acf           !! Indice de classe d'age pour AGE_CLASS_INIT_FRAC
    REAL(r_std), DIMENSION(nbpt)                         :: acnorm           !! Per-group renormalisation of acfrac
    REAL(r_std), DIMENSION(3)                            :: maskvals         !! values to use to mask (according to 
                                                                             !!   `maskingtype') 
    CHARACTER(LEN=250)                                   :: namemaskvar      !! name of the variable to use to mask 
    CHARACTER(LEN=250)                                   :: msg

!_ ================================================================================================================================

    IF (printlev_loc >= 5) PRINT *,'  In slowproc_readvegetmax'

    !Config Key   = VEGETATION_FILE
    !Config Desc  = Name of file from which the vegetation map is to be read
    !Config If    = 
    !Config Def   = PFTmap.nc
    !Config Help  = The name of the file to be opened to read a vegetation
    !Config         map (in pft) is to be given here. 
    !Config Units = [FILE]
    !
    filename = 'PFTmap.nc'
    CALL getin_p('VEGETATION_FILE',filename)
    variablename = 'maxvegetfrac'


    IF (xios_interpolation) THEN
       IF (printlev_loc >= 1) WRITE(numout,*) "slowproc_readvegetmax: Use XIOS to read and interpolate " &
            // TRIM(filename) // " for variable " // TRIM(variablename)

       CALL xios_orchidee_recv_field('frac_veget',vegetrefrac)
       CALL xios_orchidee_recv_field('frac_veget_frac',aveget_nvmap)
       aveget(:)=aveget_nvmap(:,1)

       DO ib = 1, nbpt
          IF (aveget(ib) > min_sechiba) THEN
             vegetrefrac(ib,:) = vegetrefrac(ib,:)/aveget(ib) ! intersected area normalization
          ENDIF
       ENDDO

    ELSE

      IF (printlev_loc >= 2) WRITE(numout,*) "slowproc_readvegetmax: Start interpolate " &
           // TRIM(filename) // " for variable " // TRIM(variablename)

      ! Assigning values to vmin, vmax
      vmin = 1
      vmax = nvmap*1._r_std

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
      maskingtype = 'msumrange'
      ! Values to use for the masking
      maskvals = (/ 1.-1.e-7, 0., 2. /)
      ! Name of the variable with the values for the mask in the input file (only if maskkingtype='var') (here not used)
      namemaskvar = ''
      CALL interpweight_3D(nbpt, nvmap, variabletypevals, lalo, resolution, neighbours,        &
        contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype,        &
        maskvals, namemaskvar, nvmap, 0, 1, fractype,                                 &
        -1., -1., vegetrefrac, aveget)
      IF (printlev_loc >= 5) WRITE(numout,*)'  slowproc_readvegetmax after interpeeight_3D'
   ENDIF

    !! Consider special case with age classes. 
    IF (nagec > 1)THEN
       ! If we are using age classes, we have to do this a little differently.
       ! We pass fake arrays which get the new vegetation from the maps for
       ! the PFTs ignoring the age classes, and then we set up the vegetation
       ! for the age classes using that information.
       DO ivma=1,nvmap
          veget_last(:,ivma)=SUM(veget_last_in(:,start_index(ivma):start_index(ivma)+nagec_pft(ivma)-1))
       ENDDO
    ELSE
       ! If you forget to include NAGEC in your input file, you may have an
       ! error here.
       IF(nvmap .NE. nvm)THEN
          CALL ipslerr_p (3,'slowproc_readvegetmax', 'The number of PFTs including all &
             age classes is not the same','as the number excluding age classes, despite that &
             NAGEC=1.', 'Please check the value of NAGEC in your run.def.')
       ENDIF

       ! Standard case with only 1 age class. All the PFT's are read from the
       ! file. 
       veget_last(:,:) = veget_last_in(:,:)
    END IF
    !
    ! Compute the logical for partial (only anthropic) PTFs update
    IF (ok_dgvm .AND. .NOT. init) THEN
       partial_update= .TRUE.
    ELSE
       partial_update=.FALSE.
    END IF

    IF (printlev_loc >= 5) THEN
      WRITE(numout,*)'  slowproc_readvegetmax before updating loop nbpt:', nbpt
    END IF

    IF ( .NOT. partial_update ) THEN
       ! Case for not DGVM or (DGVM and init)
       veget_next(:,:)=zero
       
       IF (printlev_loc >=3 .AND. ANY(aveget < min_sechiba)) THEN
          WRITE(numout,*) 'Some grid cells on the model grid did not have any points on the source grid in PFTmap.nc.'
          IF (init) THEN
             WRITE(numout,*) 'Initialization with full fraction of bare soil are done for the below grid cells.'
          ELSE
             WRITE(numout,*) 'Old values are kept for the below grid cells.'
          ENDIF
          WRITE(numout,*) 'List of grid cells (ib, lat, lon):'
       END IF
 
      DO ib = 1, nbpt
          IF (aveget(ib) < min_sechiba) THEN
             IF (printlev_loc >=3) WRITE(numout,*) ib,lalo(ib,1),lalo(ib,2)
             ! This grid cell didn't have any valid values on the source grid (the value was >1E10)
             ! Set this grid cell to bare soil if it is the first initialization, otherwise keep previous values
             IF (init) THEN
                veget_next(ib,1) = un
                veget_next(ib,2:nvmap) = zero
             ELSE
                veget_next(ib,:) = veget_last(ib,:)
             ENDIF
          ELSE
             veget_next(ib,:) = vegetrefrac(ib,:)           
          ENDIF
       ENDDO
    ELSE
       ! Partial update
       DO ib = 1, nbpt
          IF (aveget(ib) > min_sechiba) THEN
             ! For the case with properly interpolated grid cells (aveget>0)

             ! last veget for this point
             sum_veg=SUM(veget_last(ib,:))
             !
             ! If the DGVM is activated, only anthropic PFTs are utpdated, the others are copied from previous time-step 
             veget_next(ib,:) = veget_last(ib,:)
             
             DO jv = 2, nvmap
                IF ( .NOT. natural(jv) ) THEN       
                   veget_next(ib,jv) = vegetrefrac(ib,jv)
                ENDIF
             ENDDO

             sumvAnthro_old = zero
             sumvAnthro     = zero
             DO jv = 2, nvmap
                IF ( .NOT. natural(jv) ) THEN
                   sumvAnthro = sumvAnthro + veget_next(ib,jv)
                   sumvAnthro_old = sumvAnthro_old + veget_last(ib,jv)
                ENDIF
             ENDDO

             IF ( sumvAnthro_old < sumvAnthro ) THEN
                ! Increase of non natural vegetations (increase of agriculture)
                ! The proportion of natural PFT's must be preserved
                ! ie the sum of vegets is preserved
                !    and natural PFT / (sum of veget - sum of antropic veget)
                !    is preserved. 
                rapport = ( sum_veg - sumvAnthro ) / ( sum_veg - sumvAnthro_old )
                DO jv = 1, nvmap
                   IF ( natural(jv) ) THEN
                      veget_next(ib,jv) = veget_last(ib,jv) * rapport
                   ENDIF
                ENDDO
             ELSE
                ! Increase of natural vegetations (decrease of agriculture)
                ! The decrease of agriculture is replaced by bare soil. The DGVM will
                ! re-introduce natural PFT's.
                DO jv = 1, nvmap
                   IF ( natural(jv) ) THEN
                      veget_next(ib,jv) = veget_last(ib,jv)
                   ENDIF
                ENDDO
                veget_next(ib,1) = veget_next(ib,1) + sumvAnthro_old - sumvAnthro
             ENDIF

             ! test
             IF ( ABS( SUM(veget_next(ib,:)) - sum_veg ) > 10*EPSILON(un) ) THEN
                WRITE(numout,*) 'slowproc_readvegetmax _______'
                msg = "  No conservation of sum of veget for point "
                WRITE(numout,*) TRIM(msg), ib, ",(", lalo(ib,1),",", lalo(ib,2), ")" 
                WRITE(numout,*) "  last sum of veget ", sum_veg, " new sum of veget ",                &
                  SUM(veget_next(ib,:)), " error : ", SUM(veget_next(ib,:))-sum_veg
                WRITE(numout,*) "  Anthropic modifications : last ",sumvAnthro_old," new ",sumvAnthro     
                CALL ipslerr_p (3,'slowproc_readvegetmax',                                            &
                     &          'No conservation of sum of veget_next',                               &
                     &          "The sum of veget_next is different after reading Land Use map.",     &
                     &          '(verify the dgvm case model.)')
             ENDIF
          ELSE
             ! For the case when there was a propblem with the interpolation, aveget < min_sechiba
             WRITE(numout,*) 'slowproc_readvegetmax _______'
             WRITE(numout,*) "  No land point in the map for point ", ib, ",(", lalo(ib,1), ",",      &
               lalo(ib,2),")" 
             CALL ipslerr_p (2,'slowproc_readvegetmax',                                               &
                  &          'Problem with vegetation file for Land Use.',                            &
                  &          "No land point in the map for point",                                    & 
                  &          '(verify your land use file.)')
             veget_next(ib,:) = veget_last(ib,:)
          ENDIF
          
       ENDDO
    ENDIF
    IF (printlev_loc >= 5) WRITE(numout,*)'  slowproc_readvegetmax after updating'
    !
    frac_nobio_next (:,:) = un
    
    ! Works only for one nnobio !! (ie ice)
    DO inobio=1,nnobio
       DO jv=1,nvmap
          DO ib = 1, nbpt
             frac_nobio_next(ib,inobio) = frac_nobio_next(ib,inobio) - veget_next(ib,jv)
          ENDDO
       ENDDO
    ENDDO

    DO ib = 1, nbpt
       sum_veg = SUM(veget_next(ib,:))
       sum_nobio = SUM(frac_nobio_next(ib,:))
       IF (sum_nobio < 0.) THEN
          frac_nobio_next(ib,:) = zero
          veget_next(ib,1) = veget_next(ib,1) + sum_nobio
          sum_veg = SUM(veget_next(ib,:))
       ENDIF
       sumf = sum_veg + sum_nobio
       IF (sumf > min_sechiba) THEN
          veget_next(ib,:) = veget_next(ib,:) / sumf
          frac_nobio_next(ib,:) = frac_nobio_next(ib,:) / sumf
          norm=SUM(veget_next(ib,:))+SUM(frac_nobio_next(ib,:))
          err=norm-un
          IF (printlev_loc >=5) WRITE(numout,*) "  slowproc_readvegetmax: ib ",ib,                    &
            " SUM(veget_next(ib,:)+frac_nobio_next(ib,:))-un, sumf",err,sumf
          IF (abs(err) > -EPSILON(un)) THEN
             IF ( SUM(frac_nobio_next(ib,:)) > min_sechiba ) THEN
                frac_nobio_next(ib,1) = frac_nobio_next(ib,1) - err
             ELSE
                veget_next(ib,1) = veget_next(ib,1) - err
             ENDIF
             norm=SUM(veget_next(ib,:))+SUM(frac_nobio_next(ib,:))
             err=norm-un
             IF (printlev_loc >=5) WRITE(numout,*) "  slowproc_readvegetmax: ib ", ib,                &
               " SUM(veget_next(ib,:)+frac_nobio_next(ib,:))-un",err
             IF (abs(err) > EPSILON(un)) THEN
                WRITE(numout,*) '  slowproc_readvegetmax _______'
                WRITE(numout,*) "update : Problem with point ",ib,",(",lalo(ib,1),",",lalo(ib,2),")" 
                WRITE(numout,*) "         err(sum-1.) = ",abs(err)
                CALL ipslerr_p (2,'slowproc_readvegetmax', &
                     &          'Problem with sum vegetation + sum fracnobio for Land Use.',          &
                     &          "sum not equal to 1.", &
                     &          '(verify your land use file.)')
                aveget(ib) = -0.6
             ENDIF
          ENDIF
       ELSE
          ! sumf < min_sechiba
          WRITE(numout,*) '  slowproc_readvegetmax _______'
          WRITE(numout,*)"    No vegetation nor frac_nobio for point ", ib, ",(", lalo(ib,1), ",",    &
            lalo(ib,2),")" 
          WRITE(numout,*)"    Replaced by bare_soil !! "
          veget_next(ib,1) = un
          veget_next(ib,2:nvmap) = zero
          frac_nobio_next(ib,:) = zero
          CALL ipslerr_p (2,'slowproc_readvegetmax', &
               &          'Problem with vegetation file for Land Use.', &
               &          "No vegetation nor frac_nobio for point ", &
               &          '(verify your land use file.)')
       ENDIF
    ENDDO


    IF ( nagec > 1 ) THEN
       ! Now we need to change this into the map with age classes.        
       ! We have the total area for each PFT for next year from the map.
       ! However, we do not have the area for each age class; this is 
       ! only present in the simulation, not the  maps.  I am going to 
       ! put all of the vegetmax for a given PFT in only the youngest
       ! age class for now. This will be properly taken into account in
       ! land_cover_change_main. The reason we don't do it here is because 
       ! the age classes can change between this part of the code and between 
       ! land_cover_change main, due to growth or death.
       veget_next_out(:,:) = 0.0
       ! Guillaume M. -- The age-class split is gated on `init` ONLY. A partial restart
       ! (soil kept, veget_max dropped, PFT map re-read) still counts as a cold start for
       ! the vegetation. The ELSE branch is not a "flag off" case: it serves every annual
       ! VEGET_UPDATE call, where all-to-youngest is the convention. Taking the first
       ! branch unconditionally would reset the age structure every year.
       IF (init) THEN
          ! Guillaume M. -- Cold-start area is spread over the species age classes rather
          ! than all put in the youngest one, which would make every stand reach rotation
          ! together (synchronised clearfell wave). Fractions come from AGE_CLASS_INIT_FRAC,
          ! renormalised to 1, even by default; each class then establishes at its own
          ! diameter bin, which desynchronises the start. See MODULE_DESIGN_DIA_INV_DESYNC.md.
          DO ic_acf = 1, nagec
             acfrac(:,ic_acf) = age_class_init_frac(ic_acf)
          ENDDO
          DO ivma=1,nvmap
             ! Guillaume M. -- Renormalise PER GROUP: non-forest groups have nagec_pft = 1,
             ! and spreading 4 fractions over 1 slot would not conserve area. Normalise on
             ! the age classes the group ACTUALLY has.
             acnorm(:) = SUM(acfrac(:,1:nagec_pft(ivma)), DIM=2)
             DO iagec=0,nagec_pft(ivma)-1
                WHERE (acnorm(:) > min_sechiba)
                   veget_next_out(:,start_index(ivma)+iagec) = &
                        veget_next(:,ivma) * acfrac(:,iagec+1) / acnorm(:)
                ELSEWHERE
                   veget_next_out(:,start_index(ivma)+iagec) = &
                        veget_next(:,ivma) / REAL(nagec_pft(ivma),r_std)
                ENDWHERE
             ENDDO
          ENDDO
          WRITE(numout,*) '[DIA_INV] cold-start: veget_max reparti sur les classes', &
               ' d age pour la desync. nvmap=',nvmap,' nagec=',nagec, &
               ' fractions=',age_class_init_frac(1:nagec)
       ELSE
          DO ivma=1,nvmap
             veget_next_out(:,start_index(ivma))=veget_next(:,ivma)
          ENDDO
       ENDIF

    ELSE
       ! Standard case with 1 age class
       veget_next_out(:,:) = veget_next(:,:)
    END IF
 
    ! Write diagnostics
    CALL xios_orchidee_send_field("interp_avail_aveget",aveget)
    CALL xios_orchidee_send_field("interp_diag_vegetrefrac",vegetrefrac)
    CALL xios_orchidee_send_field("interp_diag_veget_next",veget_next)

    IF (printlev_loc >= 3) WRITE(numout,*) '  slowproc_readvegetmax ended'
    
  END SUBROUTINE slowproc_readvegetmax


!! ================================================================================================================================
!! SUBROUTINE   : slowproc_readcnleaf
!!
!>\BRIEF          Read and interpolate a map (by pft) with cn leaf ratio
!!
!! DESCRIPTION  : Note that the variables modified are not explicit INTENT(OUT), but they
!!                are module variables, so they don't need to be explicitly passed.
!!
!! RECENT CHANGE(S): 
!!
!! MAIN OUTPUT VARIABLE(S): ::cn_leaf_min_2D, ::cn_leaf_init_2D, ::cn_leaf_max_2D
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_readcnleaf(nbpt, lalo, neighbours,  resolution, contfrac)

    USE interpweight

    IMPLICIT NONE

    !
    !
    !
    !  0.1 INPUT
    !
    INTEGER(i_std), INTENT(in)                             :: nbpt            !! Number of points for which the data needs 
                                                                              !! to be interpolated
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in)             :: lalo            !! Vector of latitude and longitudes (beware of the order !)
    INTEGER(i_std), DIMENSION(nbpt,NbNeighb), INTENT(in)   :: neighbours      !! Vector of neighbours for each grid point
                                                                              !! (1=North and then clockwise)
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in)             :: resolution      !! The size in km of each grid-box in X and Y
    REAL(r_std), DIMENSION(nbpt), INTENT(in)               :: contfrac        !! Fraction of continent in the grid
    !
    !
    !  0.2 OUTPUT
    !
    
    !
    !  0.3 LOCAL
    !
    !
    CHARACTER(LEN=80) :: filename
    INTEGER(i_std) :: ib, inobio, jv
    REAL(r_std) :: sumf, err, norm
    !
    ! for DGVM case :
    REAL(r_std)                 :: sum_veg                     ! sum of vegets
    REAL(r_std)                 :: sum_nobio                   ! sum of nobios
    REAL(r_std), DIMENSION(nbpt)                         :: acnleaf          !! Availability of the soilcol interpolation
    REAL(r_std)                                          :: vmin, vmax       !! min/max values to use for the renormalization
    REAL(r_std), DIMENSION(nbpt,1)                         :: defaultvalue
    CHARACTER(LEN=80)                                    :: variablename     !! Variable to interpolate
    CHARACTER(LEN=80)                                    :: lonname, latname !! lon, lat names in input file
    REAL(r_std), DIMENSION(nvm)                          :: variabletypevals !! Values for all the types of the variable
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
    REAL(r_std), DIMENSION(3)                            :: maskvals         !! values to use to mask (according to 
                                                                             !!   `maskingtype') 
    CHARACTER(LEN=250)                                   :: namemaskvar      !! name of the variable to use to mask 
    CHARACTER(LEN=250)                                   :: msg
    REAL(r_std), DIMENSION(nbpt,nvm,1)                     :: cnleaf       !! cn leaf read

!_ ================================================================================================================================

    IF (printlev_loc >= 5) PRINT *,'  In slowproc_readcnleaf'

    !
    !Config Key   = CNLEAF_FILE
    !Config Desc  = Name of file from which the cn leaf ratio is to be read
    !Config If    = 
    !Config Def   = cnleaf_map.nc
    !Config Help  = The name of the file to be opened to read a 2D cn leaf ratio
    !Config Units = [FILE]
    !
    filename = 'cnleaf_map.nc'
    CALL getin_p('CNLEAF_FILE',filename)

    !
    !Config Key   = CNLEAF_VAR
    !Config Desc  = Name of the variable in the file from which the cn leaf ratio is to be read
    !Config If    = 
    !Config Def   = leaf_cn.nc
    !Config Help  = The name of the variable to be opened to read a 2D cn leaf ratio
    !Config Units = [VAR]
    !
    variablename = 'leaf_cn'
    CALL getin_p('CNLEAF_VAR', variablename)

    IF (printlev_loc >= 2) WRITE(numout,*) "slowproc_readcnleaf: Start interpolate " &
         // TRIM(filename) // " for variable " // TRIM(variablename)

    ! Assigning values to vmin, vmax
    vmin = 0.
    vmax = 0.

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
    maskvals = (/ 1.-1.e-7, min_sechiba, 2. /)
    ! Name of the variable with the values for the mask in the input file (only if maskkingtype='var') (here not used)
    namemaskvar = ''

!  SUBROUTINE interpweight_3D(nbpt, Nvariabletypes, variabletypes, lalo, resolution, neighbours,       &
!    contfrac, filename, varname, inlonname, inlatname, varmin, varmax, noneg, masktype,               &
!    maskvalues, maskvarname, dim1, dim2, initime, typefrac,                                           &
!    maxresollon, maxresollat, outvar3D, aoutvar)
!    CALL interpweight_3D(nbpt, nvm, variabletypevals, lalo, resolution, neighbours,        &
!      contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype,        &
!      maskvals, namemaskvar, nvm, 0, 0, fractype,                                 &
!      -1., -1., cnleaf, acnleaf)

!  SUBROUTINE interpweight_4Dcont(nbpt, dim1, dim2, lalo, resolution, neighbours,                      &
!    contfrac, filename, varname, inlonname, inlatname, varmin, varmax, noneg, masktype,               &
!    maskvalues, maskvarname, initime, typefrac, defaultvalue, defaultNOvalue,                         &
!    outvar4D, aoutvar)

    defaultvalue=0.
    CALL interpweight_4Dcont(nbpt, nvm, 1, lalo, resolution, neighbours, &
         contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype, &
         maskvals, namemaskvar, -1, fractype, defaultvalue, 0., cnleaf, acnleaf)

    IF (printlev_loc >= 5) WRITE(numout,*)'  slowproc_readcnleaf after interpeeight_3D'
    

    cn_leaf_min_2D(:,:)=cnleaf(:,:,1)
    cn_leaf_init_2D(:,:)=cnleaf(:,:,1)
    cn_leaf_max_2D(:,:)=1000.


    IF (printlev_loc >= 3) WRITE(numout,*) '  slowproc_readcnleaf ended'
    
  END SUBROUTINE slowproc_readcnleaf


!! ================================================================================================================================
!! SUBROUTINE   : slowproc_nearest
!!
!>\BRIEF         looks for nearest grid point on the fine map
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::inear
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_nearest(iml, lon5, lat5, lonmod, latmod, inear)

    !! INTERFACE DESCRIPTION
    
    !! 0.1 input variables

    INTEGER(i_std), INTENT(in)                   :: iml             !! size of the vector
    REAL(r_std), DIMENSION(iml), INTENT(in)      :: lon5, lat5      !! longitude and latitude vector, for the 5km vegmap
    REAL(r_std), INTENT(in)                      :: lonmod, latmod  !! longitude  and latitude modelled

    !! 0.2 output variables
    
    INTEGER(i_std), INTENT(out)                  :: inear           !! location of the grid point from the 5km vegmap grid
                                                                    !! closest from the modelled grid point

    !! 0.4 Local variables

    REAL(r_std)                                  :: pa, p
    REAL(r_std)                                  :: coscolat, sincolat
    REAL(r_std)                                  :: cospa, sinpa
    REAL(r_std), ALLOCATABLE, DIMENSION(:)       :: cosang
    INTEGER(i_std)                               :: i
    INTEGER(i_std), DIMENSION(1)                 :: ineartab
    INTEGER                                      :: ALLOC_ERR

!_ ================================================================================================================================

    ALLOCATE(cosang(iml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) CALL ipslerr_p(3,'slowproc_nearest','Error in allocation for cosang','','')

    pa = pi/2.0 - latmod*pi/180.0 ! dist. between north pole and the point a 
                                                      !! COLATITUDE, in radian
    cospa = COS(pa)
    sinpa = SIN(pa)

    DO i = 1, iml

       sincolat = SIN( pi/2.0 - lat5(i)*pi/180.0 ) !! sinus of the colatitude
       coscolat = COS( pi/2.0 - lat5(i)*pi/180.0 ) !! cosinus of the colatitude

       p = (lonmod-lon5(i))*pi/180.0 !! angle between a & b (between their meridian)in radians

       !! dist(i) = ACOS( cospa*coscolat + sinpa*sincolat*COS(p))
       cosang(i) = cospa*coscolat + sinpa*sincolat*COS(p) !! TL : cosang is maximum when angle is at minimal value  
!! orthodromic distance between 2 points : cosang = cosinus (arc(AB)/R), with
!R = Earth radius, then max(cosang) = max(cos(arc(AB)/R)), reached when arc(AB)/R is minimal, when
! arc(AB) is minimal, thus when point B (corresponding grid point from LAI MAP) is the nearest from
! modelled A point
    ENDDO

    ineartab = MAXLOC( cosang(:) )
    inear = ineartab(1)

    DEALLOCATE(cosang)
  END SUBROUTINE slowproc_nearest


!! ================================================================================================================================
!! SUBROUTINE   : slowproc_soilt
!!
!>\BRIEF         Interpolate the Zobler or Reynolds/USDA soil type map
!!
!! DESCRIPTION  : Read and interpolate Zobler or Reynolds/USDA soil type map. 
!!                Read and interpolate soil bulk and soil ph from file.
!!
!! RECENT CHANGE(S): Nov 2014, ADucharne
!!                   Nov 2020, Salma Tafasca and Agnes Ducharne: adding a choice for spmipexp/SPMIPEXP,
!!                             and everything needed to read all maps and assign parameter values.
!!
!! MAIN OUTPUT VARIABLE(S): ::soiltype, ::clayfraction, sandfraction, siltfraction, ::bulk, ::soilph
!!
!! REFERENCE(S) : Reynold, Jackson, and Rawls (2000). Estimating soil water-holding capacities 
!! by linking the Food and Agriculture Organization soil map of the world with global pedon
!! databases and continuous pedotransfer functions, WRR, 36, 3653-3662
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_soilt(njsc,  ks,  nvan, avan, mcr, mcs, mcfc, mcw, &
       nbpt, lalo, neighbours, resolution, contfrac, &
       soilclass, clayfraction, sandfraction, siltfraction, bulk, soil_ph)

    USE interpweight

    IMPLICIT NONE
    !
    !
    !   This subroutine should read the Zobler/Reynolds map and interpolate to the model grid. 
    !   The method is to get fraction of the three/12 main soiltypes for each grid box.
    !   For the Zobler case, also called FAO in the code, the soil fraction are going to be put 
    !   into the array soiltype in the following order : coarse, medium and fine.
    !   For the Reynolds/USDA case, the soiltype array follows the order defined in constantes_soil_var.f90
    !
    !
    !!  0.1 INPUT
    !
    INTEGER(i_std), INTENT(in)    :: nbpt                   !! Number of points for which the data needs to be interpolated
    REAL(r_std), INTENT(in)       :: lalo(nbpt,2)           !! Vector of latitude and longitudes (beware of the order !)
    INTEGER(i_std), INTENT(in)    :: neighbours(nbpt,NbNeighb)!! Vector of neighbours for each grid point
                                                              !! (1=North and then clockwise)
    REAL(r_std), INTENT(in)       :: resolution(nbpt,2)     !! The size in km of each grid-box in X and Y
    REAL(r_std), INTENT(in)       :: contfrac(nbpt)         !! Fraction of land in each grid box.
      
    !
    !  0.2 OUTPUT
    !
    INTEGER(i_std),DIMENSION (nbpt), INTENT (out)      :: njsc           !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless)
    REAL(r_std),DIMENSION (nbpt), INTENT (out)         :: ks             !! Hydraulic conductivity at saturation (mm {-1})
    REAL(r_std),DIMENSION (nbpt), INTENT (out)         :: nvan           !! Van Genuchten coeficients n (unitless)
    REAL(r_std),DIMENSION (nbpt), INTENT (out)         :: avan           !! Van Genuchten coeficients a (mm-1})
    REAL(r_std),DIMENSION (nbpt), INTENT (out)         :: mcr            !! Residual volumetric water content (m^{3} m^{-3})
    REAL(r_std),DIMENSION (nbpt), INTENT (out)         :: mcs            !! Saturated volumetric water content (m^{3} m^{-3})
    REAL(r_std),DIMENSION (nbpt), INTENT (out)         :: mcfc           !! Volumetric water content at field capacity (m^{3} m^{-3})
    REAL(r_std),DIMENSION (nbpt), INTENT (out)         :: mcw            !! Volumetric water content at wilting point (m^{3} m^{-3})

    REAL(r_std), INTENT(out)      :: soilclass(nbpt, nscm)  !! Soil type map to be created from the Zobler map
                                                            !! or a map defining the 12 USDA classes (e.g. Reynolds)
                                                            !! Holds the area of each texture class in the ORCHIDEE grid cells
                                                            !! Final unit = fraction of ORCHIDEE grid-cell (unitless)
    REAL(r_std), INTENT(out)      :: clayfraction(nbpt)     !! The fraction of clay as used by STOMATE
    REAL(r_std), INTENT(out)      :: sandfraction(nbpt)     !! The fraction of sand (for SP-MIP)
    REAL(r_std), INTENT(out)      :: siltfraction(nbpt)     !! The fraction of silt (for SP-MIP)
    REAL(r_std), INTENT(out)      :: bulk(nbpt)             !! Bulk density  as used by STOMATE
    REAL(r_std), INTENT(out)      :: soil_ph(nbpt)          !! Soil pH  as used by STOMATE
    !
    !
    !  0.3 LOCAL
    !
    REAL(r_std), DIMENSION(nbpt)        :: param            !! to be introduced in function: interpweight
    CHARACTER(LEN=80) :: filename
    INTEGER(i_std) :: ib, ilf, nbexp, i
    INTEGER(i_std) :: fopt                                  !! Nb of pts from the texture map within one ORCHIDEE grid-cell
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:) :: solt       !! Texture the different points from the input texture map 
                                                            !! in one ORCHIDEE grid cell (unitless)
    !
    ! Number of texture classes in Zobler
    !
    INTEGER(i_std), PARAMETER :: nzobler = 7                !! Nb of texture classes according in the Zobler map
    REAL(r_std),ALLOCATABLE   :: textfrac_table(:,:)        !! conversion table between the texture index
                                                            !! and the granulometric composition
    !   
    INTEGER                  :: ALLOC_ERR
    INTEGER                                              :: ntextinfile      !! number of soil textures in the in the file
    REAL(r_std), DIMENSION(:,:), ALLOCATABLE             :: textrefrac       !! text fractions re-dimensioned
    REAL(r_std), DIMENSION(nbpt)                         :: atext            !! Availability of the texture interpolation
    REAL(r_std), DIMENSION(nbpt)                         :: abulkph          !! Availability of the bulk and ph interpolation
    REAL(r_std), DIMENSION(nbpt)                         :: aparam            !! Availability of the parameter interpolation
    CHARACTER(LEN=80)                                    :: spmipexp          !! designing the number of sp-mip experiment
    CHARACTER(LEN=80)                                    :: unif_case               !! designing the model of experiment 4 (sp_mip)
    REAL(r_std)                                          :: vmin, vmax       !! min/max values to use for the 

    CHARACTER(LEN=80)                                    :: variablename     !! Variable to interpolate
    CHARACTER(LEN=80)                                    :: lonname, latname !! lon, lat name in input file
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
                                                                             !!        (normalized by maskvals(3))
                                                                             !!   'var': mask values are taken from a 
                                                                             !!     variable inside the file (>0)
    REAL(r_std), DIMENSION(3)                            :: maskvals         !! values to use to mask (according to 
                                                                             !!   `maskingtype') 
    CHARACTER(LEN=250)                                   :: namemaskvar      !! name of the variable to use to mask 
    INTEGER(i_std), DIMENSION(:), ALLOCATABLE            :: vecpos
    CHARACTER(LEN=80)                                    :: fieldname        !! name of the field read in the N input map
    REAL(r_std)                                          :: sgn              !! sum of fractions excluding glaciers and ocean

    ! For the calculation of field capacity and wilting point
    REAL(r_std),DIMENSION (nbpt)                         :: mvan             !! Van Genuchten parameter m
    REAL(r_std),DIMENSION (nbpt)                         :: psi_w            !! Matrix potential characterizing the wilting point (mm)
    REAL(r_std),DIMENSION (nbpt)                         :: psi_fc           !! Matrix potential characterizing the field capacity (mm)
    
!_ ================================================================================================================================

    IF (printlev_loc>=3) WRITE (numout,*) 'slowproc_soilt'

    ! The soil parameters are defined by several keywords in run.def:
    ! (a) soil_classif tells which kind of soil texture map you will read (mandatory):
    !    - usda for 12 USDA texture classes (Reynolds, SoilGrids, SPMIP, etc) updated to 13 classes
    !      for clay oxisols by Salma Tafasca
    !    - zobler to read teh Zobler map and reduce it to 3 classes (fine, medium, coarse)
    ! (b) spmipexp was introduced by Salma Tafasca for the SPMIP project  
    !   maps: Reading the soil parameter maps of SPMIP
    !   unif: Imposing uniform soil texture over the globe (4 texture options, with parameter values imposed by SP-MIP)
    ! Even with maps, some parameters (thermics) are defined based on texture.
    ! So we read a soil texture map in all experiments but unif, where soil texture is imposed by njsc(:).
    ! (c) unif_case to choose the soil texture assigned if spmipexp=maps (4 hard_coded possibilities)

    ! IMPORTANT: if no spmipexp is defined in run.def, the model works as before, by deriving the soil parameters 
    ! from a soil texture map, itself defined by the SOILTYPE_CLASSIF keyword, and soil_classif variable
    ! But to get a uniform texture (exp 4), you need to select a soil texture map using soil_classif, even if it's not read

    !Config Key   = SPMIPEXP
    !Config Desc  = Types of alternative hydraulic parameters
    !Config Def   = 'texture'
    !Config If    =
    !Config Help  = possible values: maps, unif
    !Config Units = [-]
    spmipexp='texture' ! default is to define parameters from soil texture, with soil_classif = 'zobler' or 'usda'
    CALL getin_p("SPMIPEXP",spmipexp)

    IF (spmipexp == 'unif') THEN
       ! case where unif=exp4 is selected: uniform soil parameters
       ! the values of the hydraulic parameters below come from SP-MIP,
       ! and correspond to the Rosetta PTF (Schaap et al., 2001)

       ! sp_mip_experiment_4: select another level of experiment: a, b, c or d in run.def

       !Config Key   = UNIF_CASE
       !Config Desc  = Types of uniform soil textures in SPMIP, possible values: a, b, c and d
       !Config Def   = 'b'
       !Config If    =
       !Config Help  = possible values: a, b, c and d
       !Config Units = [-]
       unif_case='b' ! default = loamy soil
       CALL getin_p("UNIF_CASE",unif_case)

       SELECTCASE (unif_case)

       CASE ('a') ! loamy sand
          clayfraction=0.06
          sandfraction=0.81
          siltfraction=0.13
          DO ib=1 , nbpt
             njsc(ib) = 2
             mcr(ib) = 0.049
             mcs(ib) = 0.39
             ks(ib) = (1.41e-5)*1000*24*3600
             avan(ib) = 3.475*(1e-3)
             nvan(ib) = 1.746
             mcfc(ib) = 0.1039
             mcw(ib) = 0.05221
          ENDDO

       CASE ('b') !loam
          clayfraction=0.2
          sandfraction=0.4
          siltfraction=0.4
          DO ib=1, nbpt
             njsc(ib) = 6
             mcr(ib) = 0.061
             mcs(ib) = 0.399
             ks(ib) = (3.38e-6)*1000*24*3600
             avan(ib) = 1.112*(1e-3)
             nvan(ib) = 1.472
             mcfc(ib) = 0.236
             mcw(ib) = 0.09115
          ENDDO
 
       CASE ('c') !silt
          clayfraction=0.1
          sandfraction=0.06
          siltfraction=0.84
          DO ib=1, nbpt
             njsc(ib)=5
             mcr(ib) = 0.05
             mcs(ib) = 0.489
             ks(ib) = (2.81e-6)*1000*24*3600
             avan(ib) = 0.6577*(1e-3)
             nvan(ib) = 1.679
             mcfc(ib) = 0.2854
             mcw(ib) = 0.06944
          ENDDO

       CASE ('d')!clay
          clayfraction=0.55
          sandfraction=0.15
          siltfraction=0.3
          DO ib=1, nbpt
             njsc(ib)=12
             mcr(ib) = 0.098
             mcs(ib) = 0.459
             ks(ib) = (9.74e-7)*1000*24*3600
             avan(ib) = 1.496*(1e-3)
             nvan(ib) = 1.253
             mcfc(ib) = 0.3329
             mcw(ib) = 0.1897
          ENDDO

       CASE DEFAULT

          WRITE (numout,*) 'Unsupported experiment number. Choose between a, b, c or d according to sp_mip_experiment_4 number'
          CALL ipslerr_p(3,'hydrol_init','Unsupported experiment number. ',&
               'Choose between a,b,c or d','')
       ENDSELECT

    ELSE ! spmipexp is either exp1=maps, or texture for exp2 or exp3 (or typing error!)
                                                              
    !
    !  Needs to be a configurable variable
    !
    !
    !Config Key   = SOILCLASS_FILE
    !Config Desc  = Name of file from which soil types are read
    !Config Def   = soils_param.nc
    !Config If    = NOT(IMPOSE_VEG)
    !Config Help  = The name of the file to be opened to read the soil types. 
    !Config         The data from this file is then interpolated to the grid of
    !Config         of the model. The aim is to get fractions for sand loam and
    !Config         clay in each grid box. This information is used for soil hydrology
    !Config         and respiration.
    !Config Units = [FILE]
    !
    ! soils_param.nc file is 1deg soil texture file (Zobler)
    ! The USDA map from Reynolds is soils_param_usda.nc (1/12deg resolution)


    filename = 'soils_param.nc'
    CALL getin_p('SOILCLASS_FILE',filename)

    variablename = 'soiltext'

    !! Variables for interpweight
    ! Type of calculation of cell fractions
    fractype = 'default'
    nbexp=0
    IF (xios_interpolation) THEN
       IF (printlev_loc >= 1) WRITE(numout,*) "slowproc_soilt: Use XIOS to read and interpolate " &
            // TRIM(filename) // " for variable " // TRIM(variablename)
       
       SELECT CASE(soil_classif)

       CASE('none')
          ALLOCATE(textfrac_table(nscm,ntext), STAT=ALLOC_ERR)
          IF (ALLOC_ERR/=0) CALL ipslerr_p(3,'slowproc_soilt','Error in allocation for textfrac_table','','')
          DO ib=1, nbpt
             njsc(ib) = usda_default ! 6 = Loam
             clayfraction(ib) = clayfrac_usda(usda_default) 
             sandfraction(ib) = sandfrac_usda(usda_default)
             siltfraction(ib) = 1.-clayfrac_usda(usda_default)-sandfrac_usda(usda_default)
          ENDDO

       CASE('zobler')
          IF (printlev_loc>=2) WRITE(numout,*) "Using a soilclass map with Zobler classification, to be read using XIOS"
          !
          ALLOCATE(textrefrac(nbpt,nzobler))
          ALLOCATE(textfrac_table(nzobler,ntext), STAT=ALLOC_ERR)
          IF (ALLOC_ERR/=0) CALL ipslerr_p(3,'slowproc_soilt','Error in allocation for textfrac_table','','')
          CALL get_soilcorr_zobler (nzobler, textfrac_table)
       
          CALL xios_orchidee_recv_field('soiltext1',textrefrac(:,1))
          CALL xios_orchidee_recv_field('soiltext2',textrefrac(:,2))
          CALL xios_orchidee_recv_field('soiltext3',textrefrac(:,3))
          CALL xios_orchidee_recv_field('soiltext4',textrefrac(:,4))
          CALL xios_orchidee_recv_field('soiltext5',textrefrac(:,5))
          CALL xios_orchidee_recv_field('soiltext6',textrefrac(:,6))
          CALL xios_orchidee_recv_field('soiltext7',textrefrac(:,7))
       
          CALL get_soilcorr_zobler (nzobler, textfrac_table)
        !
        !
          DO ib =1, nbpt   
            soilclass(ib,:)=0. 
            soilclass(ib,fao2usda(1))=textrefrac(ib,1) 
            soilclass(ib,fao2usda(2))=textrefrac(ib,2)+textrefrac(ib,3)+textrefrac(ib,4)+textrefrac(ib,7) 
            soilclass(ib,fao2usda(3))=textrefrac(ib,5) 
                     
            ! clayfraction is the sum of the % of clay (as a mineral of small granulometry, and not as a texture)
            ! over the zobler pixels composing the ORCHIDEE grid-cell
            clayfraction(ib) = textfrac_table(1,3) * textrefrac(ib,1)+textfrac_table(2,3) * textrefrac(ib,2) + &
                               textfrac_table(3,3) * textrefrac(ib,3)+textfrac_table(4,3) * textrefrac(ib,4) + &
                               textfrac_table(5,3) * textrefrac(ib,5)+textfrac_table(7,3) * textrefrac(ib,7)

            sandfraction(ib) = textfrac_table(1,2) * textrefrac(ib,1)+textfrac_table(2,2) * textrefrac(ib,2) + &
                               textfrac_table(3,2) * textrefrac(ib,3)+textfrac_table(4,2) * textrefrac(ib,4) + &
                               textfrac_table(5,2) * textrefrac(ib,5)+textfrac_table(7,2) * textrefrac(ib,7)

            siltfraction(ib) = textfrac_table(1,1) * textrefrac(ib,1)+textfrac_table(2,1) * textrefrac(ib,2) + &
                               textfrac_table(3,1) * textrefrac(ib,3)+textfrac_table(4,1) * textrefrac(ib,4) + &
                               textfrac_table(5,1) * textrefrac(ib,5)+textfrac_table(7,1) * textrefrac(ib,7)

            sgn=SUM(soilclass(ib,:)) ! grid-cell fraction with texture info
            
            IF (sgn < min_sechiba) THEN ! if no texture info in this grid-point, we assume that texture = Loam
              njsc(ib) = usda_default ! 6 = Loam
              clayfraction(ib) = clayfrac_usda(usda_default) 
              sandfraction(ib) = sandfrac_usda(usda_default)
              siltfraction(ib) = 1.-clayfrac_usda(usda_default)-sandfrac_usda(usda_default)  
              atext(ib)=0.
            ELSE
              atext(ib)=sgn
              clayfraction(ib) = clayfraction(ib) / sgn
              sandfraction(ib) = sandfraction(ib) / sgn
              siltfraction(ib) = siltfraction(ib) / sgn
              soilclass(ib,:)  = soilclass(ib,:) / sgn
              njsc(ib) = MAXLOC(soilclass(ib,:),1) ! Dominant texture class
            ENDIF
            
          ENDDO
          
          
       
       CASE('usda')

          IF (printlev_loc>=4) WRITE (numout,*) 'slowproc_soilt: start case usda'

          IF (printlev_loc>=1) WRITE(numout,*) "Using a soilclass map with usda classification, to be read using XIOS"
          !
          ALLOCATE(textrefrac(nbpt,nscm))
          ALLOCATE(textfrac_table(nscm,ntext), STAT=ALLOC_ERR)
          IF (ALLOC_ERR/=0) CALL ipslerr_p(3,'slowproc_soilt','Error in allocation for textfrac_table','','')
          
          CALL get_soilcorr_usda (nscm, textfrac_table)
          IF (printlev_loc>=4) WRITE (numout,*) 'slowproc_soilt: After get_soilcorr_usda'

          CALL xios_orchidee_recv_field('soiltext1',textrefrac(:,1))
          CALL xios_orchidee_recv_field('soiltext2',textrefrac(:,2))
          CALL xios_orchidee_recv_field('soiltext3',textrefrac(:,3))
          CALL xios_orchidee_recv_field('soiltext4',textrefrac(:,4))
          CALL xios_orchidee_recv_field('soiltext5',textrefrac(:,5))
          CALL xios_orchidee_recv_field('soiltext6',textrefrac(:,6))
          CALL xios_orchidee_recv_field('soiltext7',textrefrac(:,7))
          CALL xios_orchidee_recv_field('soiltext8',textrefrac(:,8))
          CALL xios_orchidee_recv_field('soiltext9',textrefrac(:,9))
          CALL xios_orchidee_recv_field('soiltext10',textrefrac(:,10))
          CALL xios_orchidee_recv_field('soiltext11',textrefrac(:,11))
          CALL xios_orchidee_recv_field('soiltext12',textrefrac(:,12))
          CALL xios_orchidee_recv_field('soiltext13',textrefrac(:,13))
       

          clayfraction(:) = 0.0
          sandfraction(:) = 0.0
          siltfraction(:) = 0.0
          nbexp=0
          DO ib =1, nbpt
            DO ilf = 1,nscm
              soilclass(ib,ilf)=textrefrac(ib,ilf)      
              clayfraction(ib) = clayfraction(ib) + textfrac_table(ilf,3)*textrefrac(ib,ilf)
              sandfraction(ib) = sandfraction(ib) + textfrac_table(ilf,2)*textrefrac(ib,ilf)
              siltfraction(ib) = siltfraction(ib) + textfrac_table(ilf,1)*textrefrac(ib,ilf)
              ! textfrac_table holds the %silt,%sand,%clay
            ENDDO


            sgn=SUM(soilclass(ib,:)) ! grid-cell fraction with texture info
            
            IF (sgn < min_sechiba) THEN ! if no texture info in this grid-point, we assume that texture = Loam
              njsc(ib) = usda_default ! 6 = Loam
              clayfraction(ib) = clayfrac_usda(usda_default) 
              sandfraction(ib) = sandfrac_usda(usda_default) 
              siltfraction(ib) = 1.-clayfrac_usda(usda_default)-sandfrac_usda(usda_default)  
              atext(ib)=0
            ELSE
              soilclass(ib,:) = soilclass(ib,:) / sgn
              clayfraction(ib) = clayfraction(ib) / sgn
              sandfraction(ib) = sandfraction(ib) / sgn
              siltfraction(ib) = siltfraction(ib) / sgn
              atext(ib)=sgn
              njsc(ib) = MAXLOC(soilclass(ib,:),1) ! Dominant texture class
            END IF          
         ENDDO

        CASE DEFAULT
             WRITE(numout,*) 'slowproc_soilt:'
             WRITE(numout,*) '  A non supported soil type classification has been chosen'
             CALL ipslerr_p(3,'slowproc_soilt','non supported soil type classification','','')
        END SELECT

    ELSE              !    xios_interpolation 
       ! Read and interpolate using stardard method with IOIPSL and aggregate
    
       IF (printlev_loc >= 1) WRITE(numout,*) "slowproc_soilt: Read and interpolate " &
            // TRIM(filename) // " for variable " // TRIM(variablename)

       ! Name of the longitude and latitude in the input file
       lonname = 'nav_lon'
       latname = 'nav_lat'

       IF (printlev_loc >= 2) WRITE(numout,*) "slowproc_soilt: Start interpolate " &
            // TRIM(filename) // " for variable " // TRIM(variablename)

       IF ( TRIM(soil_classif) /= 'none' ) THEN

          ! Define a variable for the number of soil textures in the input file
          SELECTCASE(soil_classif)
          CASE('zobler')
             ntextinfile=nzobler
          CASE('usda')
             ntextinfile=nscm
          CASE DEFAULT
             WRITE(numout,*) 'slowproc_soilt:'
             WRITE(numout,*) '  A non supported soil type classification has been chosen'
             CALL ipslerr_p(3,'slowproc_soilt','non supported soil type classification','','')
          ENDSELECT

          ALLOCATE(textrefrac(nbpt,ntextinfile), STAT=ALLOC_ERR)
          IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'slowproc_soilt','Problem in allocation of variable textrefrac',&
            '','')

          ! Assigning values to vmin, vmax
          vmin = un
          vmax = ntextinfile*un
          
          ALLOCATE(variabletypevals(ntextinfile), STAT=ALLOC_ERR)
          IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'slowproc_soilt','Problem in allocation of variabletypevals','','')
          variabletypevals = -un
          
          !! Variables for interpweight
          ! Should negative values be set to zero from input file?
          nonegative = .FALSE.
          ! Type of mask to apply to the input data (see header for more details)
          maskingtype = 'mabove'
          ! Values to use for the masking
          maskvals = (/ min_sechiba, undef_sechiba, undef_sechiba /)
          ! Name of the variable with the values for the mask in the input file (only if maskkingtype='var') ( not used)
          namemaskvar = ''

          CALL interpweight_2D(nbpt, ntextinfile, variabletypevals, lalo, resolution, neighbours,        &
             contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype,    & 
             maskvals, namemaskvar, 0, 0, -1, fractype, -1., -1., textrefrac, atext)

          ALLOCATE(vecpos(ntextinfile), STAT=ALLOC_ERR)
          IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'slowproc_soilt','Problem in allocation of variable vecpos','','')
          ALLOCATE(solt(ntextinfile), STAT=ALLOC_ERR)
          IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'slowproc_soilt','Problem in allocation of variable solt','','')
          
          IF (printlev_loc >= 5) THEN
             WRITE(numout,*)'  slowproc_soilt after interpweight_2D'
             WRITE(numout,*)'  slowproc_soilt before starting loop nbpt:', nbpt
             WRITE(numout,*)"  slowproc_soilt starting classification '" // TRIM(soil_classif) // "'..."
          END IF
       ELSE
         IF (printlev_loc >= 5) WRITE(numout,*)'  slowproc_soilt using default values all points are propertly ' // &
           'interpolated atext = 1. everywhere!'
         atext = 1.
       END IF

    nbexp = 0
    SELECTCASE(soil_classif)
    CASE('none')
       ALLOCATE(textfrac_table(nscm,ntext), STAT=ALLOC_ERR)
       IF (ALLOC_ERR/=0) CALL ipslerr_p(3,'slowproc_soilt','Error in allocation for textfrac_table','','')
       DO ib=1, nbpt
          njsc(ib) = usda_default ! 6 = Loam
          clayfraction(ib) = clayfrac_usda(usda_default) 
          sandfraction(ib) = sandfrac_usda(usda_default) 
          siltfraction(ib) = 1.-clayfrac_usda(usda_default)-sandfrac_usda(usda_default)  
       ENDDO
    CASE('zobler')
       IF (printlev_loc>=2) WRITE(numout,*) "Using a soilclass map with Zobler classification"
       !
       ALLOCATE(textfrac_table(nzobler,ntext), STAT=ALLOC_ERR)
       IF (ALLOC_ERR/=0) CALL ipslerr_p(3,'slowproc_soilt','Error in allocation for textfrac_table','','')
       CALL get_soilcorr_zobler (nzobler, textfrac_table)
       !
       IF (printlev_loc >= 5) WRITE(numout,*)'  slowproc_soilt after getting table of textures'
       DO ib =1, nbpt
          soilclass(ib,:) = zero
          clayfraction(ib) = zero
          sandfraction(ib) = zero
          siltfraction(ib) = zero
          !
          ! vecpos: List of positions where textures were not zero
          ! vecpos(1): number of not null textures found
          vecpos = interpweight_ValVecR(textrefrac(ib,:),nzobler,zero,'neq')
          fopt = vecpos(1)

          IF ( fopt .EQ. 0 ) THEN
             ! No points were found for current grid box, use default values
             nbexp = nbexp + 1
             njsc(ib) = usda_default ! 6=Loam
             clayfraction(ib) = clayfrac_usda(usda_default) 
             sandfraction(ib) = sandfrac_usda(usda_default) 
             siltfraction(ib) = 1.-clayfrac_usda(usda_default)-sandfrac_usda(usda_default)  
          ELSE
             IF (fopt == nzobler) THEN
                ! All textures are not zero
                solt=(/(i,i=1,nzobler)/)
             ELSE
               DO ilf = 1,fopt
                 solt(ilf) = vecpos(ilf+1)
               END DO
             END IF
             !
             !   Compute the fraction of each textural class
             !
             sgn = 0.
             DO ilf = 1,fopt
                   !
                   ! Here we make the correspondance between the 7 zobler textures and the 3 textures in ORCHIDEE
                   ! and soilclass correspond to surfaces covered by the 3 textures of ORCHIDEE (coase,medium,fine)
                   ! For type 6 = glacier, default values are set and it is also taken into account during the normalization 
                   ! of the fractions (done in interpweight_2D)
                   ! Note that type 0 corresponds to ocean but it is already removed using the mask above.
                   !
                IF ( (solt(ilf) .LE. nzobler) .AND. (solt(ilf) .GT. 0) .AND. & 
                     (solt(ilf) .NE. 6) ) THEN
                   SELECT CASE(solt(ilf))
                     CASE(1)
                        soilclass(ib,fao2usda(1)) = soilclass(ib,fao2usda(1)) + textrefrac(ib,solt(ilf))
                     CASE(2)
                        soilclass(ib,fao2usda(2)) = soilclass(ib,fao2usda(2)) + textrefrac(ib,solt(ilf))
                     CASE(3)
                        soilclass(ib,fao2usda(2)) = soilclass(ib,fao2usda(2)) + textrefrac(ib,solt(ilf))
                     CASE(4)
                        soilclass(ib,fao2usda(2)) = soilclass(ib,fao2usda(2)) + textrefrac(ib,solt(ilf))
                     CASE(5)
                        soilclass(ib,fao2usda(3)) = soilclass(ib,fao2usda(3)) + textrefrac(ib,solt(ilf))
                     CASE(7)
                        soilclass(ib,fao2usda(2)) = soilclass(ib,fao2usda(2)) + textrefrac(ib,solt(ilf))
                     CASE DEFAULT
                        WRITE(numout,*) 'We should not be here, an impossible case appeared'
                        CALL ipslerr_p(3,'slowproc_soilt','Bad value for solt','','')
                   END SELECT
                   ! clayfraction is the sum of the % of clay (as a mineral of small granulometry, and not as a texture)
                   ! over the zobler pixels composing the ORCHIDEE grid-cell
                   clayfraction(ib) = clayfraction(ib) + &
                        & textfrac_table(solt(ilf),3) * textrefrac(ib,solt(ilf))
                   sandfraction(ib) = sandfraction(ib) + &
                        & textfrac_table(solt(ilf),2) * textrefrac(ib,solt(ilf))
                   siltfraction(ib) = siltfraction(ib) + &
                        & textfrac_table(solt(ilf),1) * textrefrac(ib,solt(ilf))
                   ! Sum the fractions which are not glaciers nor ocean
                   sgn = sgn + textrefrac(ib,solt(ilf))
                ELSE
                   IF (solt(ilf) .GT. nzobler) THEN
                      WRITE(numout,*) 'The file contains a soil color class which is incompatible with this program'
                      CALL ipslerr_p(3,'slowproc_soilt','Problem soil color class incompatible','','')

                   ENDIF
                END IF
             ENDDO

             IF ( (sgn .LT. min_sechiba) .OR. (atext(ib) .LT. min_sechiba)) THEN
                ! Set default values if grid cells were only covered by glaciers or ocean 
                ! or if now information on the source grid was found.
                nbexp = nbexp + 1
                njsc(ib) = usda_default ! 6 = Loam
                clayfraction(ib) = clayfrac_usda(usda_default) 
                sandfraction(ib) = sandfrac_usda(usda_default) 
                siltfraction(ib) = 1.-clayfrac_usda(usda_default)-sandfrac_usda(usda_default)  
             ELSE
                ! Normalize using the fraction of surface not including glaciers and ocean
                soilclass(ib,:) = soilclass(ib,:)/sgn
                clayfraction(ib) = clayfraction(ib)/sgn
                sandfraction(ib) = sandfraction(ib)/sgn
                siltfraction(ib) = siltfraction(ib)/sgn
                njsc(ib) = MAXLOC(soilclass(ib,:),1) ! Dominant texture class
             ENDIF
          ENDIF
       ENDDO
      
       ! The "USDA" case reads a map of the 12 USDA texture classes, 
       ! such as to assign the corresponding soil properties
       CASE("usda")
          IF (printlev_loc>=2) WRITE(numout,*) "Using a soilclass map with usda classification"
 
          ALLOCATE(textfrac_table(nscm,ntext), STAT=ALLOC_ERR)
          IF (ALLOC_ERR/=0) CALL ipslerr_p(3,'slowproc_soilt','Error in allocation for textfrac_table','','')

          CALL get_soilcorr_usda (nscm, textfrac_table)

          IF (printlev_loc>=4) WRITE (numout,*) 'slowproc_soilt: After get_soilcorr_usda'
          !
          DO ib =1, nbpt
          ! GO through the point we have found
          !
          !
          ! Provide which textures were found
          ! vecpos: List of positions where textures were not zero
          !   vecpos(1): number of not null textures found
          vecpos = interpweight_ValVecR(textrefrac(ib,:),ntextinfile,zero,'neq')
          fopt = vecpos(1)
         
          !
          !    Check that we found some points
          !
          soilclass(ib,:) = 0.0
          clayfraction(ib) = 0.0
          sandfraction(ib) = 0.0
          siltfraction(ib) = 0.0
          
          IF ( fopt .EQ. 0) THEN
             ! No points were found for current grid box, use default values
             IF (printlev_loc>=3) WRITE(numout,*)'slowproc_soilt: no soil class in input file found for point=', ib
             nbexp = nbexp + 1
             njsc(ib) = usda_default ! 6 = Loam
             clayfraction(ib) = clayfrac_usda(usda_default) 
             sandfraction(ib) = sandfrac_usda(usda_default) 
             siltfraction(ib) = 1.-clayfrac_usda(usda_default)-sandfrac_usda(usda_default)  
          ELSE
             IF (fopt == nscm) THEN
                ! All textures are not zero
                solt(:) = (/(i,i=1,nscm)/)
             ELSE
               DO ilf = 1,fopt
                 solt(ilf) = vecpos(ilf+1) 
               END DO
             END IF
                !
                !   Compute the fraction of each textural class  
                DO ilf = 1,fopt
                   IF ( (solt(ilf) .LE. nscm) .AND. (solt(ilf) .GT. 0) ) THEN
                      soilclass(ib,solt(ilf)) = textrefrac(ib,solt(ilf))
                      clayfraction(ib) = clayfraction(ib) + textfrac_table(solt(ilf),3) *                &
                           textrefrac(ib,solt(ilf))
                      sandfraction(ib) = sandfraction(ib) + textfrac_table(solt(ilf),2) * &
                           textrefrac(ib,solt(ilf))
                      siltfraction(ib) = siltfraction(ib) + textfrac_table(solt(ilf),1) * &
                        textrefrac(ib,solt(ilf))
                   ELSE
                      IF (solt(ilf) .GT. nscm) THEN
                         WRITE(numout,*) 'The file contains a soil color class which is incompatible with this program'
                         CALL ipslerr_p(3,'slowproc_soilt','Problem soil color class incompatible 2','','')
                      ENDIF
                   ENDIF
                   !
                ENDDO

                
                njsc(ib) = MAXLOC(soilclass(ib,:),1) ! Dominant texture class
                
                ! Set default values if the surface in source file is too small
                ! Warning - This test is donne differently for Zobler (based on sgn, related to class 6=ice)
                IF ( atext(ib) .LT. min_sechiba) THEN
                   nbexp = nbexp + 1
                   njsc(ib) = usda_default ! 6 = Loam
                   clayfraction(ib) = clayfrac_usda(usda_default) 
                   sandfraction(ib) = sandfrac_usda(usda_default) 
                   siltfraction(ib) = 1.-clayfrac_usda(usda_default)-sandfrac_usda(usda_default)                
                ENDIF
            ENDIF

          ENDDO
          
          IF (printlev_loc>=4) WRITE (numout,*) '  slowproc_soilt: End case usda'
          
       CASE DEFAULT
          WRITE(numout,*) 'slowproc_soilt _______'
          WRITE(numout,*) '  A non supported soil type classification has been chosen'
          CALL ipslerr_p(3,'slowproc_soilt','non supported soil type classification','','')
       ENDSELECT
       IF (printlev_loc >= 5 ) WRITE(numout,*)'  slowproc_soilt end of type classification'


       IF (ALLOCATED(variabletypevals)) DEALLOCATE (variabletypevals)
       IF (ALLOCATED(textrefrac)) DEALLOCATE (textrefrac)
       IF (ALLOCATED(solt)) DEALLOCATE (solt)
       IF (ALLOCATED(textfrac_table)) DEALLOCATE (textfrac_table)
    
    ENDIF        !      xios_interpolation 


    IF ( nbexp .GT. 0 .AND. printlev>=1) THEN
       WRITE(numout,*) 'slowproc_soilt: '
       WRITE(numout,*) '  The interpolation of variable soiltext had ', nbexp
       WRITE(numout,*) '  points without data. This are either coastal points or ice covered land.'
       WRITE(numout,*) '  The problem was solved by using the default soil types.'
    ENDIF
       
       ! End of soil texture reading, for 'maps' and classical behavior
       
       IF (spmipexp == 'maps') THEN
              IF (printlev_loc>=3) WRITE (numout,*) 'slowproc_soilt: Read soil hydraulic parameters with IOIPSL'

              ! Read using IOIPSL and interpolate using aggregate tool in ORCHIDEE

              !Config Key   = PARAM_FILE
              !Config Desc  = Name of file from which soil parameter  values are read
              !Config Def   = params_sp_mip.nc
              !Config Help  = The name of the file to be opened to read values of parameters.
              !Config         The data from this file is then interpolated to the grid of
              !Config         of the model.
              !Config Units = [FILE]
              !
              ! params_sp_mip.nc file is 0.5 deg soil hydraulic parameters file provided by sp_mip

              filename = 'params_sp_mip.nc'
              CALL getin_p('PARAM_FILE',filename)

              !! Variables for interpweight
              ! Type of calculation of cell fractions
              fractype = 'default'
              ! Name of the longitude and latitude in the input file
              lonname = 'nav_lon'
              latname = 'nav_lat'
              ! Assigning values to vmin, vmax (there are not types/categories
              vmin =0.
              vmax = 99999.
              !! Variables for interpweight
              ! Should negative values be set to zero from input file?
              nonegative = .FALSE.
              ! Type of mask to apply to the input data (see header for more details)
              maskingtype = 'mabove'
              ! Values to use for the masking
              maskvals = (/ min_sechiba, undef_sechiba, undef_sechiba /)
              ! Name of the variable with the values for the mask in the input file (only if maskkingtype='var') ( not used)
              namemaskvar = ''

              variablename = 'ks'
              IF (printlev_loc >= 1) WRITE(numout,*) "slowproc_soilt: Read and interpolate " &
                   // TRIM(filename) // " for variable " // TRIM(variablename)
              CALL interpweight_2Dcont(nbpt, 0, 0, lalo, resolution, neighbours,                                &
                   contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype,     &
                   maskvals, namemaskvar, -1, fractype, 0., 0.,                              &
                   ks, aparam)
              IF (printlev_loc>=1) WRITE(numout,*) 'ks map is read _______'

              variablename = 'alpha'
              CALL interpweight_2Dcont(nbpt, 0, 0, lalo, resolution, neighbours,                                &
                   contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype,     &
                   maskvals, namemaskvar, -1, fractype, 0., 0.,                              &
                   avan, aparam)
              IF (printlev_loc>=1) WRITE(numout,*) 'avan map read _______'

              variablename = 'thetar'
              CALL interpweight_2Dcont(nbpt, 0, 0, lalo, resolution, neighbours,                                &
                   contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype,     &
                   maskvals, namemaskvar, -1, fractype, 0., 0.,                              &
                   mcr, aparam)
              IF (printlev_loc>=1) WRITE(numout,*) 'thetar map read _______'

              variablename = 'thetas'
              CALL interpweight_2Dcont(nbpt, 0, 0, lalo, resolution, neighbours,                                &
                   contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype,     &
                   maskvals, namemaskvar, -1, fractype, 0., 0.,                              &
                   mcs, aparam)
              IF (printlev_loc>=1) WRITE(numout,*) 'thetas map read _______'

              variablename = 'thetapwpvg' ! mcw
              CALL interpweight_2Dcont(nbpt, 0, 0, lalo, resolution, neighbours,                                &
                   contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype,     &
                   maskvals, namemaskvar, -1, fractype, 0., 0.,                              &
                   mcw, aparam)
              IF (printlev_loc>=1) WRITE(numout,*) 'thetapwpvg map read _______'

              variablename = 'thetafcvg' !mcfc
              CALL interpweight_2Dcont(nbpt, 0, 0, lalo, resolution, neighbours,                                &
                   contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype,     &
                   maskvals, namemaskvar, -1, fractype, 0., 0.,                              &
                   mcfc, aparam)
              IF (printlev_loc>=1) WRITE(numout,*) 'thetafcvg map read _______'

              variablename = 'nvg'
              CALL interpweight_2Dcont(nbpt, 0, 0, lalo, resolution, neighbours,                                &
                   contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype,     &
                   maskvals, namemaskvar, -1, fractype, 0., 0.,                              &
                   nvan, aparam)
              IF (printlev_loc>=1) WRITE(numout,*) 'nvan map read _______'

       ELSE ! spmipexp is not maps nor unif, then it must be texture
          IF (spmipexp == 'texture') THEN
             ! Whichever the soil texture map, we can use the USDA parameter vectors with 13 values 
             nvan(:) = nvan_usda(njsc(:))
             avan(:) = avan_usda(njsc(:))
             mcr(:) = mcr_usda(njsc(:))
             mcs(:) = mcs_usda(njsc(:))
             ks(:) = ks_usda(njsc(:))
!!$             mcfc(:) = mcf_usda(njsc(:))
!!$             mcw(:) = mcw_usda(njsc(:))
             
             !! Calculation of FC and WP based on above 5 parameters
             mvan(:) = un - (un / nvan(:))
             ! Define matrix potential in mm for wilting point and field capacity (with sand vs clay-silt variation)
             psi_w(:) = 150000.
             DO ib=1, nbpt
                IF ( ks(ib) .GE. 560 ) THEN ! Sandy soils (560 is equivalent of 2.75 at log scale of Ks, mm/d)
                   psi_fc(ib) = 1000.
                ELSE ! Finer soils
                   psi_fc(ib) = 3300. 
                ENDIF
             ENDDO
             mcfc(:) = mcr(:) + (( mcs(:) - mcr(:)) / (un + ( avan(:) * psi_fc(:))** nvan(:))** mvan(:))
             mcw(:)  = mcr(:) + (( mcs(:) - mcr(:)) / (un + ( avan(:) *  psi_w(:))** nvan(:))** mvan(:))
             
         ELSE ! if spmipexp is not among texture or maps or unif
            WRITE(numout,*) "Unsupported spmipexp=",spmipexp
            WRITE(numout,*) "Choose between texture, maps, and unif"
            CALL ipslerr_p(3,'soilproc_soilt','Bad choice of spmipexp','Choose between texture, maps, and unif','')
         ENDIF
       ENDIF
   ENDIF ! SPMIPEXP 


!!
!! Read and interpolate soil bulk and soil ph using IOIPSL or XIOS
!!
    IF (xios_interpolation) THEN
       ! Read and interpolate using XIOS

       ! Check if the restart file for sechiba is read.
       ! Reading of soilbulk and soilph with XIOS is only activated if restname==NONE.
       IF (restname_in /= 'NONE') THEN
          CALL ipslerr_p(3,'slowproc_soilt','soilbulk and soilph can not be read with XIOS if sechiba restart file exist',&
               'Remove sechiba restart file and start again','')
       END IF

       IF (printlev_loc>=3) WRITE (numout,*) 'slowproc_soilt: Read soilbulk and soilph with XIOS'
       CALL xios_orchidee_recv_field('soilbulk', bulk)
       CALL xios_orchidee_recv_field('soilph', soil_ph)

    ELSE
       ! Read using IOIPSL and interpolate using aggregate tool in ORCHIDEE
       IF (printlev_loc>=3) WRITE (numout,*) 'slowproc_soilt: Read soilbulk and soilph with IOIPSL'

       !! Read soilbulk

       !Config Key   = SOIL_BULK_FILE
       !Config Desc  = Name of file from which soil bulk should be read
       !Config Def   = soil_bulk_and_ph.nc
       !Config If    = 
       !Config Help  = 
       !Config Units = [FILE]
       
       ! By default, bulk and ph is stored in the same file but they could be separated if needed.
       filename = 'soil_bulk_and_ph.nc'
       CALL getin_p('SOIL_BULK_FILE',filename)
       
       fieldname= 'soilbulk'
       ! Name of the longitude and latitude in the input file
       lonname = 'nav_lon'
       latname = 'nav_lat'
       vmin=0  ! not used in interpweight_2Dcont
       vmax=0  ! not used in interpweight_2Dcont
       
       ! Should negative values be set to zero from input file?
       nonegative = .FALSE.
       ! Type of mask to apply to the input data (see header for more details)
       maskingtype = 'mabove'
       ! Values to use for the masking
       maskvals = (/ min_sechiba, undef_sechiba, undef_sechiba /)
       ! Name of the variable with the values for the mask in the input file (only if maskkingtype='var') ( not used)
       namemaskvar = ''
       ! Type of calculation of cell fractions
       fractype = 'default'
       CALL interpweight_2Dcont(nbpt, 0, 0, lalo, resolution, neighbours,                        &
            contfrac, filename, fieldname, lonname, latname, vmin, vmax, nonegative, maskingtype,    &
            maskvals, namemaskvar, -1, fractype, bulk_default, undef_sechiba,                        &
            bulk, abulkph)
       
       !! Read soilph
       
       !Config Key   = SOIL_PH_FILE
       !Config Desc  = Name of file from which soil ph should be read
       !Config Def   = soil_bulk_and_ph.nc
       !Config If    = 
       !Config Help  = 
       !Config Units = [FILE]
       
       filename = 'soil_bulk_and_ph.nc'
       CALL getin_p('SOIL_PH_FILE',filename)
       
       fieldname= 'soilph'
       ! Name of the longitude and latitude in the input file
       lonname = 'nav_lon'
       latname = 'nav_lat'
       CALL interpweight_2Dcont(nbpt, 0, 0, lalo, resolution, neighbours,                        &
            contfrac, filename, fieldname, lonname, latname, vmin, vmax, nonegative, maskingtype,    &
            maskvals, namemaskvar, -1, fractype, ph_default, undef_sechiba,                          &
            soil_ph, abulkph)
              
    END IF  ! xios_interpolation


    ! Write diagnostics
    CALL xios_orchidee_send_field("interp_avail_atext",atext)
    CALL xios_orchidee_send_field("interp_diag_soilclass",soilclass)
    CALL xios_orchidee_send_field("interp_diag_njsc",REAL(njsc, r_std))
    CALL xios_orchidee_send_field("interp_diag_clayfraction",clayfraction)
    CALL xios_orchidee_send_field("interp_diag_sandfraction",sandfraction)
    CALL xios_orchidee_send_field("interp_diag_siltfraction",siltfraction)
    CALL xios_orchidee_send_field("interp_diag_bulk",bulk)
    CALL xios_orchidee_send_field("interp_diag_soil_ph",soil_ph)

    IF (printlev_loc >= 3) WRITE(numout,*) '  slowproc_soilt ended'

  END SUBROUTINE slowproc_soilt

!! ================================================================================================================================
!! SUBROUTINE   : slowproc_slope
!!
!>\BRIEF         Calculate mean slope coef in each  model grid box from the slope map
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::reinf_slope
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_slope(nbpt, lalo, neighbours, resolution, contfrac, reinf_slope)

    USE interpweight

    IMPLICIT NONE

    !
    !
    !
    !  0.1 INPUT
    !
    INTEGER(i_std), INTENT(in)           :: nbpt                      !! Number of points for which the data needs to be interpolated
    REAL(r_std), INTENT(in)              :: lalo(nbpt,2)              !! Vector of latitude and longitudes (beware of the order !)
    INTEGER(i_std), INTENT(in)           :: neighbours(nbpt,NbNeighb) !! Vector of neighbours for each grid point
                                                                      !! (1=North and then clockwise)
    REAL(r_std), INTENT(in)              :: resolution(nbpt,2)        !! The size in km of each grid-box in X and Y
    REAL(r_std), INTENT (in)             :: contfrac(nbpt)            !! Fraction of continent in the grid
    !
    !  0.2 OUTPUT
    !
    REAL(r_std), INTENT(out)             ::  reinf_slope(nbpt)        !! slope coef 
    !
    !  0.3 LOCAL
    !
    !
    REAL(r_std)                                          :: slope_noreinf    !! Slope above which runoff is maximum
    CHARACTER(LEN=80)                                    :: filename
    REAL(r_std)                                          :: vmin, vmax       !! min/max values to use for the 
                                                                             !!   renormalization
    REAL(r_std), DIMENSION(nbpt)                         :: aslope           !! slope availability 

    CHARACTER(LEN=80)                                    :: variablename     !! Variable to interpolate
    CHARACTER(LEN=80)                                    :: lonname, latname !! lon, lat name in the input file
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
                                                                             !!     variable inside the file  (>0)
    REAL(r_std), DIMENSION(3)                            :: maskvals         !! values to use to mask (according to 
                                                                             !!   `maskingtype') 
    CHARACTER(LEN=250)                                   :: namemaskvar      !! name of the variable to use to mask 

!_ ================================================================================================================================
    
    
    !Config Key   = SLOPE_NOREINF
    !Config Desc  = Slope over which surface runoff does not reinfiltrate
    !Config If    = 
    !Config Def   = 0.5
    !Config Help  = The slope above which there is no reinfiltration
    !Config Units = [%]
    !
    slope_noreinf = 0.5
    CALL getin_p('SLOPE_NOREINF',slope_noreinf)
    !
    !Config Key   = TOPOGRAPHY_SLOPE_FILE
    !Config Desc  = Name of file from which the slope of the topography map is to be read
    !Config Def   = cartepente2d_15min.nc
    !Config If    = 
    !Config Help  = The name of the file to be opened to read the orography
    !Config         map is to be given here. Usualy SECHIBA runs with a 2'
    !Config         map which is derived from the NGDC one. 
    !Config Units = [FILE]
    !
    filename = 'cartepente2d_15min.nc'
    CALL getin_p('TOPOGRAPHY_SLOPE_FILE',filename)

    IF (xios_interpolation) THEN
    
      CALL xios_orchidee_recv_field('reinf_slope_interp',reinf_slope)
      CALL xios_orchidee_recv_field('frac_slope_interp',aslope)


    ELSE
    
      variablename = 'pente'
      IF (printlev_loc >= 1) WRITE(numout,*) "slowproc_slope: Read and interpolate " &
           // TRIM(filename) // " for variable " // TRIM(variablename)

      ! For this case there are not types/categories. We have 'only' a continuos field
      ! Assigning values to vmin, vmax
      vmin = 0.
      vmax = 9999.

      !! Variables for interpweight
      ! Type of calculation of cell fractions
      fractype = 'slopecalc'
      ! Name of the longitude and latitude in the input file
      lonname = 'longitude'
      latname = 'latitude'
      ! Should negative values be set to zero from input file?
      nonegative = .FALSE.
      ! Type of mask to apply to the input data (see header for more details)
      maskingtype = 'mabove'
      ! Values to use for the masking
      maskvals = (/ min_sechiba, undef_sechiba, undef_sechiba /)
      ! Name of the variable with the values for the mask in the input file (only if maskkingtype='var') (here not used)
      namemaskvar = ''

      CALL interpweight_2Dcont(nbpt, 0, 0, lalo, resolution, neighbours,                                &
        contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype,        &
        maskvals, namemaskvar, -1, fractype, slope_default, slope_noreinf,                              &
        reinf_slope, aslope)
      IF (printlev_loc >= 5) WRITE(numout,*)'  slowproc_slope after interpweight_2Dcont'

    ENDIF
    
      ! Write diagnostics
    CALL xios_orchidee_send_field("interp_avail_aslope",aslope)
    CALL xios_orchidee_send_field("interp_diag_reinf_slope",reinf_slope)

    IF (printlev_loc >= 3) WRITE(numout,*) '  slowproc_slope ended'

  END SUBROUTINE slowproc_slope


!! ================================================================================================================================
!! SUBROUTINE	: slowproc_xios_initialize_ninput
!!
!>\BRIEF        Activates or not reading of ninput file
!!
!!
!! DESCRIPTION  : This subroutine activates or not the variables in the xml files related to the reading of input files. 
!!                The subroutine is called from slowproc_xios_initialization only if xios_orchidee_ok is activated. Therefor
!!                the reading from run.def done in here are duplications and will be done again in slowproc_Ninput.
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCE(S)	: None.
!! 
!! FLOWCHART : None.
!! \n
!_ ================================================================================================================================
 
  SUBROUTINE slowproc_xios_initialize_ninput(Ninput_field, flag)

    CHARACTER(LEN=*), INTENT(in)         :: Ninput_field   !! Name of the default field reading in the map
    LOGICAL, INTENT(in)                  :: flag           !! Condition where the file will be read
    CHARACTER(LEN=80)                    :: filename
    CHARACTER(LEN=80)                    :: varname
    INTEGER                              :: ninput_update_loc
    INTEGER                              :: l
    CHARACTER(LEN=30)                    :: ninput_str     
    
    ! Read from run.def file and variable name for the current field: Ninput_field_FILE, Ninput_field_VAR
    filename = 'NONE'
    CALL getin_p(TRIM(Ninput_field)//'_FILE',filename)

    varname=Ninput_field
    CALL getin_p(TRIM(Ninput_field)//'_VAR',varname)
          
    ! Read from run.def the number of years set in the variable NINPUT_UPDATE
    ninput_update_loc=0
    WRITE(ninput_str,'(a)') '0Y'
    CALL getin_p('NINPUT_UPDATE', ninput_str)
    l=INDEX(TRIM(ninput_str),'Y')
    READ(ninput_str(1:(l-1)),"(I2.2)") ninput_update_loc

    ! Determine if reading with XIOS will be done in this executaion.
    ! Activate files and fields in the xml files if reading will be done with
    ! XIOS in this run. Otherwise deactivate the files. 
    IF (flag .AND. (restname_in=='NONE' .OR. (ninput_update_loc>0)) .AND. &
         (TRIM(filename) .NE. 'NONE') .AND. (TRIM(filename) .NE. 'none')) THEN
       IF (xios_interpolation) THEN
          ! Reading will be done with XIOS later
          IF (printlev>=1) WRITE(numout,*) 'Reading of ',TRIM(Ninput_field), &
                       ' will be done later with XIOS. File and variable name are ',filename, varname
          CALL xios_orchidee_set_file_attr(TRIM(Ninput_field)//'_file',enabled=.TRUE., name=filename(1:LEN_TRIM(filename)-3))
          CALL xios_orchidee_set_field_attr(TRIM(Ninput_field)//'_read',enabled=.TRUE., name=TRIM(varname))
          CALL xios_orchidee_set_field_attr('mask_'//TRIM(Ninput_field)//'_read',enabled=.TRUE., name=TRIM(varname))
       ELSE
          ! Reading will be done with IOIPSL later
          ! Deactivate file specification in xml files
          IF (printlev>=1) WRITE(numout,*) 'Reading of ',TRIM(Ninput_field), &
                       ' will be done with IOIPSL. File and variable name are ',filename, varname
          CALL xios_orchidee_set_file_attr(TRIM(Ninput_field)//'_file',enabled=.FALSE.)
          CALL xios_orchidee_set_field_attr(TRIM(Ninput_field)//'_interp',enabled=.FALSE.)
       END IF
    ELSE
       ! No reading will be done, deactivate corresponding file declared in context_input_orchidee.xml
       IF (printlev>=1) WRITE(numout,*) 'No reading of ',TRIM(Ninput_field),' will be done'
       CALL xios_orchidee_set_file_attr(TRIM(Ninput_field)//'_file',enabled=.FALSE.)
       CALL xios_orchidee_set_field_attr(TRIM(Ninput_field)//'_interp',enabled=.FALSE.)

       ! Deactivate controle output diagnostic field not needed since no interpolation
       CALL xios_orchidee_set_field_attr("interp_diag_"//TRIM(Ninput_field),enabled=.FALSE.)
    END IF
  
  END SUBROUTINE slowproc_xios_initialize_ninput


!! ================================================================================================================================
!! SUBROUTINE	: slowproc_Ninput
!!
!>\BRIEF        Reads in the maps containing nitrogen inputs
!!
!!
!! DESCRIPTION  : This subroutine reads in various maps containing information on the amount of nitrogen inputed
!!              to the system via manure, fertilizer, atmospheric deposition, and biological nitrogen fixation.
!!              The information is read in for a single year for all pixels present in the simulation, and
!!              interpolated to the resolution being used for the current run.
!!
!! RECENT CHANGE(S): 
!!
!! MAIN OUTPUT VARIABLE(S): Ninput_vec
!!
!! REFERENCE(S)	: None.
!! 
!! FLOWCHART : None.
!! \n
!_ ================================================================================================================================
 


  SUBROUTINE slowproc_Ninput(nbpt,         lalo,       neighbours,  resolution, contfrac, &
                             Ninput_field, Ninput_vec, Ninput_year)

    !
    !! 0. Variable and parameter declaration
    !

    !
    !! 0.1 Input variables
    !
    INTEGER(i_std), INTENT(in)                           :: nbpt           !! Number of points for which the data needs to be interpolated
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in)           :: lalo           !! Vector of latitude and longitudes (beware of the order !)
    INTEGER(i_std), DIMENSION(nbpt,8), INTENT(in)        :: neighbours     !! Vector of neighbours for each grid point 
    ! (1=N, 2=NE, 3=E, 4=SE, 5=S, 6=SW, 7=W, 8=NW)
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in)           :: resolution     !! The size in km of each grid-box in X and Y
    REAL(r_std), DIMENSION(nbpt), INTENT(in)             :: contfrac       !! Fraction of continent in the grid
    CHARACTER(LEN=80), INTENT(in)                        :: Ninput_field   !! Name of the default field reading in the map
    INTEGER(i_std), INTENT(in)                           :: Ninput_year    !! year for N inputs update

    !
    !! 0.2 Modified variables
    !

    !
    !! 0.3 Output variables
    !
    REAL(r_std), DIMENSION(nbpt, nvm,12), INTENT(out)    ::  Ninput_vec    !! Nitrogen input (kgN m-2 yr-1) 

    !! 0.4 Local variables
    !
    CHARACTER(LEN=80)                                    :: filename
    CHARACTER(LEN=30)                                    :: callsign
    INTEGER(i_std)                                       :: iml, jml, lml, tml, fid, ib, ip, jp, vid, l, im
    INTEGER(i_std)                                       :: idi, idi_last, nbvmax
    REAL(r_std)                                          :: coslat
    REAL(r_std), DIMENSION(12)                           :: Ninput_val 
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:)          :: mask
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:,:)        :: sub_index
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)             :: lat_rel, lon_rel 
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:)           :: Ninput_map,Ninput_map_temp
    REAL(r_std), ALLOCATABLE, DIMENSION(:)               :: lat_lu, lon_lu, lon_lu_temp
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)             :: sub_area
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:)           :: resol_lu
    REAL(r_std)                                          :: Ninput_read(nbpt,12)       !! Nitrogen input temporary variable  
    INTEGER(i_std)                                       :: nix, njx, iv, i
    !
    LOGICAL                                              :: ok_interpol = .FALSE.      !! optionnal return of aggregate_2d
    !
    INTEGER                                              :: ALLOC_ERR
    CHARACTER(LEN=80)                                    :: Ninput_field_read          !! Name of the field reading in the map
    CHARACTER(LEN=80)                                    :: Ninput_year_str            !! Ninput year as a string variable
    LOGICAL                                              :: latitude_exists, longitude_exists !! Test existence of variables in the input files
!_ ================================================================================================================================



    !Config Key   = NINPUT File
    !Config Desc  = Name of file from which the N-input map is to be read
    !Config Def   = NONE
    !Config If    = 
    !Config Help  = The name of the file to be opened to read the N-input map
    !Config Units = [FILE]
    !
    filename = 'NONE'
    CALL getin_p(TRIM(Ninput_field)//'_FILE',filename)

    !Config Key   = NINPUT var
    !Config Desc  = Name of the variable in the file from which the N-input map is to be read
    !Config Def   = 'Ninput_field'
    !Config If    = 
    !Config Help  = The name of the variable  to be read for the N-input map
    !Config Units = [FILE]
    !
    Ninput_field_read=Ninput_field
    CALL getin_p(TRIM(Ninput_field)//'_VAR',Ninput_field_read)
    !
    IF((TRIM(filename) .NE. 'NONE') .AND. (TRIM(filename) .NE. 'none')) THEN

       IF(Ninput_suffix_year) THEN
          l=INDEX(TRIM(filename),'.nc')
          WRITE(Ninput_year_str,'(i4)') Ninput_year
          filename=TRIM(filename(1:(l-1)))//'_'//Ninput_year_str//'.nc'
       ENDIF

       IF (xios_interpolation) THEN

          ! Read and interpolate with XIOS
          IF (TRIM(Ninput_field)=='Nammonium' .OR. TRIM(Ninput_field)=='Nnitrate' .OR. &
              TRIM(Ninput_field)=='DRYNHX' .OR. TRIM(Ninput_field)=='WETNHX' .OR. &
              TRIM(Ninput_field)=='DRYNOY' .OR. TRIM(Ninput_field)=='WETNOY' ) THEN
             ! For these 2 fields, 12 time step exist in the file
             CALL xios_orchidee_recv_field(TRIM(Ninput_field)//'_interp',Ninput_read)
          ELSE
             ! For the other fields, only 1 time step exist in the file
             CALL xios_orchidee_recv_field(TRIM(Ninput_field)//'_interp',Ninput_read(:,1))
             DO i=2,12
                Ninput_read(:,i) = Ninput_read(:,1)
             END DO
          END IF

       ELSE

          IF (printlev >= 1) WRITE(numout,*) 'Reading the variable ',TRIM(Ninput_field_read), ' from file ', TRIM(filename)
          ! Read with IOIPSL and interpolate with aggregate
          IF (is_root_prc) CALL flininfo(filename, iml, jml, lml, tml, fid)
          CALL bcast(iml)
          CALL bcast(jml)
          CALL bcast(lml)
          CALL bcast(tml)
          ALLOCATE(lat_lu(jml), STAT=ALLOC_ERR)
          IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'slowproc_ninput','Problem in allocation of variable lat_lu','','')
          
          ALLOCATE(lon_lu(iml), STAT=ALLOC_ERR)
          IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'slowproc_ninput','Problem in allocation of variable lon_lu','','')

          ALLOCATE(lon_lu_temp(iml), STAT=ALLOC_ERR)
          IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'slowproc_ninput','Problem in allocation of variable lon_lu_temp','','')
          
          ALLOCATE(Ninput_map(iml,jml,tml), STAT=ALLOC_ERR)
          IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'slowproc_ninput','Problem in allocation of variable Ninput_map','','')
          
          ALLOCATE(Ninput_map_temp(iml,jml,tml), STAT=ALLOC_ERR)
          IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'slowproc_ninput','Problem in allocation of variable Ninput_map_temp','','')
          
          ALLOCATE(resol_lu(iml,jml,2), STAT=ALLOC_ERR)
          IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'slowproc_ninput','Problem in allocation of variable resol_lu','','')
          
          
          IF (is_root_prc) THEN
             CALL flinquery_var(fid, 'longitude', longitude_exists)
             IF(longitude_exists)THEN
                CALL flinget(fid, 'longitude', iml, 0, 0, 0, 1, 1, lon_lu)
             ELSE
                CALL flinget(fid, 'lon', iml, 0, 0, 0, 1, 1, lon_lu)
             ENDIF
             CALL flinquery_var(fid, 'latitude', latitude_exists)
             IF(latitude_exists)THEN
                CALL flinget(fid, 'latitude', jml, 0, 0, 0, 1, 1, lat_lu)
             ELSE
                CALL flinget(fid, 'lat', jml, 0, 0, 0, 1, 1, lat_lu)
             ENDIF
             CALL flinget(fid, Ninput_field_read, iml, jml, 0, tml, 1, tml, Ninput_map)
             !
             CALL flinclo(fid)

             IF((ALL(lon_lu(:).GE.0.)).AND.(lon_lu(iml).GT.lon_lu(1))) THEN
                IF((lon_lu(iml/2).LT.180.).AND.(lon_lu(iml/2+1).GE.180.)) THEN
                   lon_lu_temp(iml/2+1:iml)=lon_lu(1:iml/2)
                   lon_lu_temp(1:iml/2)=lon_lu(iml/2+1:iml)
                   Ninput_map_temp(iml/2+1:iml,:,:)=Ninput_map(1:iml/2,:,:)
                   Ninput_map_temp(1:iml/2,:,:)=Ninput_map(iml/2+1:iml,:,:)
                   lon_lu=lon_lu_temp
                   Ninput_map=Ninput_map_temp
                   WHERE(lon_lu(:).GE.180.)
                      lon_lu(:)=lon_lu(:)-360.
                   ENDWHERE
                ELSE
                   CALL ipslerr_p(3,'slowproc_ninput','Problem in reading specific map 1:360','','')
                ENDIF
             ENDIF
          ENDIF
          CALL bcast(lon_lu)
          CALL bcast(lat_lu)
          CALL bcast(Ninput_map)
          
          
          ALLOCATE(lon_rel(iml,jml), STAT=ALLOC_ERR)
          IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'slowproc_slope','Problem in allocation of variable lon_rel','','')
          
          ALLOCATE(lat_rel(iml,jml), STAT=ALLOC_ERR)
          IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'slowproc_slope','Problem in allocation of variable lat_rel','','')
          
          DO ip=1,iml
             lat_rel(ip,:) = lat_lu(:)
          ENDDO
          DO jp=1,jml
             lon_rel(:,jp) = lon_lu(:)
          ENDDO
          !
          !
          ! Mask of permitted variables.
          !
          ALLOCATE(mask(iml,jml), STAT=ALLOC_ERR)
          IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'slowproc_slope','Problem in allocation of variable mask','','')
          
          mask(:,:) = zero
          DO ip=1,iml
             DO jp=1,jml
                IF (ANY(Ninput_map(ip,jp,:) .GE. 0.)) THEN
                   mask(ip,jp) = un
                ENDIF
                !
                ! Resolution in longitude
                !
                coslat = MAX( COS( lat_rel(ip,jp) * pi/180. ), mincos )     
                IF ( ip .EQ. 1 ) THEN
                   resol_lu(ip,jp,1) = ABS( lon_rel(ip+1,jp) - lon_rel(ip,jp) ) * pi/180. * R_Earth * coslat
                ELSEIF ( ip .EQ. iml ) THEN
                   resol_lu(ip,jp,1) = ABS( lon_rel(ip,jp) - lon_rel(ip-1,jp) ) * pi/180. * R_Earth * coslat
                ELSE
                   resol_lu(ip,jp,1) = ABS( lon_rel(ip+1,jp) - lon_rel(ip-1,jp) )/2. * pi/180. * R_Earth * coslat
                ENDIF
                !
                ! Resolution in latitude
                !
                IF ( jp .EQ. 1 ) THEN
                   resol_lu(ip,jp,2) = ABS( lat_rel(ip,jp) - lat_rel(ip,jp+1) ) * pi/180. * R_Earth
                ELSEIF ( jp .EQ. jml ) THEN
                   resol_lu(ip,jp,2) = ABS( lat_rel(ip,jp-1) - lat_rel(ip,jp) ) * pi/180. * R_Earth
                ELSE
                   resol_lu(ip,jp,2) =  ABS( lat_rel(ip,jp-1) - lat_rel(ip,jp+1) )/2. * pi/180. * R_Earth
                ENDIF
                !
             ENDDO
          ENDDO

          !
          !
          ! The number of maximum vegetation map points in the GCM grid is estimated.
          ! Some lmargin is taken.
          !
          IF (is_root_prc) THEN
             nix=INT(MAXVAL(resolution_g(:,1))/MAXVAL(resol_lu(:,:,1)))+2
             njx=INT(MAXVAL(resolution_g(:,2))/MAXVAL(resol_lu(:,:,2)))+2
             nbvmax = nix*njx
          ENDIF
          CALL bcast(nbvmax)
          !
          callsign="Ninput map"
          ok_interpol = .FALSE.
          DO WHILE ( .NOT. ok_interpol )
             !
             IF (printlev >= 2) WRITE(numout,*) "Projection arrays for ",callsign," : "
             IF (printlev >= 2) WRITE(numout,*) "nbvmax = ",nbvmax
             
             ALLOCATE(sub_index(nbpt,nbvmax,2), STAT=ALLOC_ERR)
             IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'slowproc_Ninput','Problem in allocation of variable sub_index','','')
             sub_index(:,:,:)=0
             
             ALLOCATE(sub_area(nbpt,nbvmax), STAT=ALLOC_ERR)
             IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'slowproc_Ninput','Problem in allocation of variable sub_area','','')
             sub_area(:,:)=zero
             
             CALL aggregate_p(nbpt, lalo, neighbours, resolution, contfrac, &
                  &                iml, jml, lon_rel, lat_rel, mask, callsign, &
                  &                nbvmax, sub_index, sub_area, ok_interpol)
             
             IF (.NOT. ok_interpol ) THEN
                IF (printlev_loc>=3) WRITE(numout,*) 'nbvmax will be increased from ',nbvmax,' to ', nbvmax*2
                DEALLOCATE(sub_area)
                DEALLOCATE(sub_index)
                nbvmax = nbvmax * 2
             END IF
          END DO
          !
          !
          DO ib = 1, nbpt
             Ninput_val(:) = zero
             
             ! Initialize last index to the highest possible 
             idi_last=nbvmax
             DO idi=1, nbvmax
                ! Leave the do loop if all sub areas are treated, sub_area <= 0
                IF ( sub_area(ib,idi) <= zero ) THEN
                   ! Set last index to the last one used
                   idi_last=idi-1
                   ! Exit do loop
                   EXIT
                END IF
                
                ip = sub_index(ib,idi,1)
                jp = sub_index(ib,idi,2)
                
                IF(tml == 12) THEN
                   Ninput_val(:) = Ninput_val(:) + Ninput_map(ip,jp,:) * sub_area(ib,idi)
                ELSE
                   Ninput_val(:) = Ninput_val(:) + Ninput_map(ip,jp,1) * sub_area(ib,idi)
                ENDIF
             ENDDO
             
             IF ( idi_last >= 1 ) THEN
                Ninput_read(ib,:) = Ninput_val(:) / SUM(sub_area(ib,1:idi_last)) 
             ELSE
                CALL ipslerr_p(2,'slowproc_ninput', '', '',&
                     &                 'No information for a point') ! Warning error
                Ninput_read(ib,:) = 0.
             ENDIF
          ENDDO
          !
          DEALLOCATE(Ninput_map)
          DEALLOCATE(Ninput_map_temp)
          DEALLOCATE(sub_index)
          DEALLOCATE(sub_area)
          DEALLOCATE(mask)
          DEALLOCATE(lon_lu)
          DEALLOCATE(lon_lu_temp)
          DEALLOCATE(lat_lu)
          DEALLOCATE(lon_rel)
          DEALLOCATE(lat_rel)
          
       END IF ! xios_interpolation 
       ! Output the variables read for control only
       IF (TRIM(Ninput_field)=='Nammonium' .OR. TRIM(Ninput_field)=='Nnitrate' .OR. &
           TRIM(Ninput_field)=='WETNHX' .OR. TRIM(Ninput_field)=='DRYNHX' .OR. &
           TRIM(Ninput_field)=='WETNOY' .OR. TRIM(Ninput_field)=='DRYNOY' ) THEN
          ! For these 2 fields, 12 time step exist in the file
          CALL xios_orchidee_send_field("interp_diag_"//TRIM(Ninput_field),Ninput_read)
       ELSE
          CALL xios_orchidee_send_field("interp_diag_"//TRIM(Ninput_field),Ninput_read(:,1))
       END IF   
       !       
       ! Initialize Ninput_vec
       Ninput_vec(:,:,:) = 0.
       SELECT CASE (Ninput_field)
          CASE ("Nammonium")
             DO iv = 1,nvm 
                Ninput_vec(:,iv,:) = Ninput_read(:,:)
             ENDDO
          CASE ("Nnitrate")
             DO iv = 1,nvm 
                Ninput_vec(:,iv,:) = Ninput_read(:,:)
             ENDDO
          CASE ("WETNHX")
             DO iv = 1,nvm 
                Ninput_vec(:,iv,:) = Ninput_read(:,:)
             ENDDO
          CASE ("DRYNHX")
             DO iv = 1,nvm 
                Ninput_vec(:,iv,:) = Ninput_read(:,:)
             ENDDO
          CASE ("WETNOY")
             DO iv = 1,nvm 
                Ninput_vec(:,iv,:) = Ninput_read(:,:)
             ENDDO
          CASE ("DRYNOY")
             DO iv = 1,nvm 
                Ninput_vec(:,iv,:) = Ninput_read(:,:)
             ENDDO
          CASE ("Nfert")
             DO iv = 2,nvm
                ! Exclude bare soil. It is not fertilized
                IF ( .NOT. natural(iv) ) THEN
                   Ninput_vec(:,iv,:) = Ninput_read(:,:)
                ENDIF
             ENDDO
          CASE ("Nfert_cropland")
             DO iv = 2,nvm
                ! Exclude bare soil. It is not fertilized
                IF ( .NOT. natural(iv) ) THEN
                   Ninput_vec(:,iv,:) = Ninput_read(:,:)
                ENDIF
             ENDDO
          CASE ("Nfert_ammo_cropland")
             DO iv = 2,nvm 
                IF ( .NOT. natural(iv) ) THEN
                   Ninput_vec(:,iv,:) = Ninput_read(:,:)
                ENDIF
             ENDDO
          CASE ("Nfert_nitr_cropland")
             DO iv = 2,nvm 
                IF ( .NOT. natural(iv) ) THEN
                   Ninput_vec(:,iv,:) = Ninput_read(:,:)
                ENDIF
             ENDDO
          CASE ("Nfert_cropC3")
             DO iv = 2,nvm 
                IF (( .NOT. natural(iv) ).AND.( .NOT. is_C4(iv))) THEN
                   Ninput_vec(:,iv,:) = Ninput_read(:,:)
                ENDIF
             ENDDO
          CASE ("Nfert_cropC4")
             DO iv = 2,nvm 
                IF (( .NOT. natural(iv) ).AND.( is_C4(iv))) THEN
                   Ninput_vec(:,iv,:) = Ninput_read(:,:)
                ENDIF
             ENDDO
          CASE ("Nmanure_cropland")
             DO iv = 2,nvm
                ! Exclude bare soil. It is not fertilized
                IF ( .NOT. natural(iv) ) THEN
                   Ninput_vec(:,iv,:) = Ninput_read(:,:)
                ENDIF
             ENDDO
          CASE ("Nfert_pasture")
             DO iv = 2,nvm
                ! Exclude bare soil. It is not fertilized
                IF ( natural(iv) .AND. (.NOT.(is_tree(iv))) ) THEN
                   Ninput_vec(:,iv,:) = Ninput_read(:,:)
                ENDIF
             ENDDO
          CASE ("Nfert_ammo_pasture")
             DO iv = 2,nvm 
                IF ( natural(iv) .AND. (.NOT.(is_tree(iv))) ) THEN
                   Ninput_vec(:,iv,:) = Ninput_read(:,:)
                ENDIF
             ENDDO
          CASE ("Nfert_nitr_pasture")
             DO iv = 2,nvm 
                IF ( natural(iv) .AND. (.NOT.(is_tree(iv))) ) THEN
                   Ninput_vec(:,iv,:) = Ninput_read(:,:)
                ENDIF
             ENDDO
          CASE ("Nmanure_pasture")
             DO iv = 2,nvm
                ! Exclude bare soil. It is not fertilized
                IF ( natural(iv) .AND. (.NOT.(is_tree(iv))) ) THEN
                   Ninput_vec(:,iv,:) = Ninput_read(:,:)
                ENDIF
             ENDDO
          CASE ("Nbnf")
             DO iv = 2,nvm
                ! Exclude bare soil. No plants = no biological N fixation
                Ninput_vec(:,iv,:) = Ninput_read(:,:)
             ENDDO
          CASE default
             WRITE (numout,*) 'This kind of Ninput_field choice is not possible. '
             CALL ipslerr_p(3,'slowproc_ninput', '', '',&
               &              'This kind of Ninput_field choice is not possible.') ! Fatal error
       END SELECT
       !
       IF (printlev >= 1) WRITE(numout,*) 'Interpolation Done in slowproc_Ninput for ',TRIM(Ninput_field)
       !
       !
    ELSE
       Ninput_vec(:,:,:)=zero
    ENDIF 
  END SUBROUTINE slowproc_Ninput

!! ================================================================================================================================
!! SUBROUTINE   : slowproc_woodharvest
!!
!>\BRIEF         
!!
!! DESCRIPTION  : 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_woodharvest(nbpt, lalo, neighbours, resolution, contfrac, woodharvest)

    USE interpweight

    IMPLICIT NONE

    !
    !
    !
    !  0.1 INPUT
    !
    INTEGER(i_std), INTENT(in)                           :: nbpt         !! Number of points for which the data needs to be interpolated
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in)           :: lalo         !! Vector of latitude and longitudes (beware of the order !)
    INTEGER(i_std), DIMENSION(nbpt,NbNeighb), INTENT(in) :: neighbours   !! Vector of neighbours for each grid point
                                                                         !! (1=North and then clockwise)
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in)           :: resolution   !! The size in km of each grid-box in X and Y
    REAL(r_std), DIMENSION(nbpt), INTENT(in)             :: contfrac     !! Fraction of continent in the grid
    !
    !  0.2 OUTPUT
    !
    REAL(r_std), DIMENSION(nbpt), INTENT(out)            ::  woodharvest !! Wood harvest
    !
    !  0.3 LOCAL
    !
    CHARACTER(LEN=80)                                    :: filename
    REAL(r_std)                                          :: vmin, vmax  
    REAL(r_std), DIMENSION(nbpt)                         :: aoutvar          !! availability of input data to
                                                                             !!   interpolate output variable 
                                                                             !!   (on the nbpt space)
    CHARACTER(LEN=80)                                    :: variablename     !! Variable to interpolate
    CHARACTER(LEN=80)                                    :: lonname, latname !! lon, lat name in the input file
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
                                                                             !!     variable inside the file  (>0)
    REAL(r_std), DIMENSION(3)                            :: maskvals         !! values to use to mask (according to 
                                                                             !!   `maskingtype') 
    CHARACTER(LEN=250)                                   :: namemaskvar      !! name of the variable to use to mask 
    REAL(r_std), DIMENSION(1)                            :: variabletypevals !! 
!    REAL(r_std), DIMENSION(nbp_mpi)                      :: woodharvest_mpi  !! Wood harvest where all thredds OMP are gatherd
!_ ================================================================================================================================
    
    
    !Config Key   = WOODHARVEST_FILE
    !Config Desc  = Name of file from which the wood harvest will be read
    !Config Def   = woodharvest.nc
    !Config If    = DO_WOOD_HARVEST
    !Config Help  = 
    !Config Units = [FILE]
    filename = 'woodharvest.nc'
    CALL getin_p('WOODHARVEST_FILE',filename)
    variablename = 'woodharvest'


    IF (xios_interpolation) THEN
       IF (printlev_loc >= 1) WRITE(numout,*) "slowproc_readwoodharvest: Use XIOS to read and interpolate " &
            // TRIM(filename) // " for variable " // TRIM(variablename)

       CALL xios_orchidee_recv_field('woodharvest_interp',woodharvest)

       aoutvar = 1.0
    ELSE

       IF (printlev_loc >= 1) WRITE(numout,*) "slowproc_readwoodharvest: Read and interpolate " &
            // TRIM(filename) // " for variable " // TRIM(variablename)

       ! For this case there are not types/categories. We have 'only' a continuos field
       ! Assigning values to vmin, vmax
       vmin = 0.
       vmax = 9999.
       
       !! Variables for interpweight
       ! Type of calculation of cell fractions
       fractype = 'default'
       ! Name of the longitude and latitude in the input file
       lonname = 'longitude'
       latname = 'latitude'
       ! Should negative values be set to zero from input file?
       nonegative = .TRUE.
       ! Type of mask to apply to the input data (see header for more details)
       maskingtype = 'nomask'
       ! Values to use for the masking
       maskvals = (/ min_sechiba, undef_sechiba, undef_sechiba /)
       ! Name of the variable with the values for the mask in the input file (only if maskkingtype='var') (here not used)
       namemaskvar = ''
       
       variabletypevals=-un
       CALL interpweight_2Dcont(nbpt, 0, 0, lalo, resolution, neighbours,                                &
            contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype,        &
            maskvals, namemaskvar, -1, fractype, 0., 0., woodharvest, aoutvar)
       IF (printlev_loc >= 5) WRITE(numout,*)'  slowproc_wodharvest after interpweight_2Dcont'
       
    END IF

    ! Write diagnostics
    CALL xios_orchidee_send_field("interp_diag_woodharvest",woodharvest)

    IF (printlev_loc >= 3) WRITE(numout,*) '  slowproc_woodharvest ended'
  END SUBROUTINE slowproc_woodharvest

!! ================================================================================================================================
!! FUNCTION 	: slowproc_impose_site_level_lai
!!
!>\BRIEF        ! This subroutine aims at bypassing stomate when the user wants it.
!!
!! DESCRIPTION   : This subroutine calculates the main output values from stomate that sechiba will use in its different modules.
!!                 The user prescribes the global LAI and height and all the different variables are calculates from this
!!
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : circ_class_biomass, heights
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_impose_site_level_lai(kjit,kjpindex,veget_max,coszang, temp_growth, max_height_store, height, height_dom, lai, fco2_flux_out, &
       fco2_lu_out, fco2_wh_out, fco2_ha_out, veget, qsintmax, circ_class_biomass, circ_class_n, lai_per_level, z_array_out, laieff_fit)

    !! 0. Variables and parameters declaration

    !! 0.1 Input variables

    INTEGER(i_std),                                                      INTENT(in)      :: kjit                     !! Time step number
    INTEGER(i_std),                                                      INTENT(in)      :: kjpindex                 !! Domain size - terrestrial pixels only
    REAL(r_std),        DIMENSION(kjpindex,nvm),                         INTENT(in)      :: veget_max                !! Maximum vegetation fraction cover (-)
    REAL(r_std),        DIMENSION(kjpindex),                             INTENT(in)      :: coszang                  !! The cosine of the solar zenith angle (unitless) 


    !! 0.2 Output variables 

    REAL(r_std),        DIMENSION(kjpindex),                             INTENT(out)     :: temp_growth              !! Growth temperature (°C) - Is equal to t2m_month 
    REAL(r_std),        DIMENSION(kjpindex,nvm),                         INTENT(out)     :: max_height_store         !! Maximum height of canopy
    REAL(r_std),        DIMENSION(kjpindex,nvm),                         INTENT(out)     :: height                   !! Height of vegetation (m)
    REAL(r_std),        DIMENSION(kjpindex,nvm),                         INTENT(out)     :: height_dom               !! Height of dominant class vegetation (m)
    REAL(r_std),        DIMENSION(kjpindex,nvm),                         INTENT(out)     :: lai                      !! Leaf area index (m^2 m^{-2})
    REAL(r_std),        DIMENSION(kjpindex),                             INTENT(out)     :: fco2_flux_out            !! CO2 flux per average ground area (gC m^{-2} dt_stomate^{-1})
    REAL(r_std),        DIMENSION(kjpindex),                             INTENT(out)     :: fco2_lu_out              !! CO2 flux from land-use (without forest management) (gC m^{-2} dt_stomate^{-1})
    REAL(r_std),        DIMENSION(kjpindex),                             INTENT(out)     :: fco2_wh_out              !! CO2 Flux to Atmosphere from Wood Harvesting (gC m^{-2} dt_stomate^{-1})
    REAL(r_std),        DIMENSION(kjpindex),                             INTENT(out)     :: fco2_ha_out              !! CO2 Flux to Atmosphere from Crop Harvesting (gC m^{-2} dt_stomate^{-1})
    

    !! 0.3 Modified variables

    REAL(r_std),        DIMENSION(kjpindex,nvm),                         INTENT(inout)   :: veget                    !! Fraction of vegetation type including none biological fractioning the mesh (unitless)
    REAL(r_std),        DIMENSION(kjpindex,nvm),                         INTENT(inout)   :: qsintmax                 !! Maximum water storage on vegetation from interception (mm)
    REAL(r_std),        DIMENSION(kjpindex,nvm,ncirc,nparts,nelements),  INTENT(inout)   :: circ_class_biomass       !! Biomass components of the model tree  
                                                                                                                     !! within a circumference class @tex $(g C ind^{-1})$ @endtex
    REAL(r_std),        DIMENSION(kjpindex,nvm,ncirc),                   INTENT(inout)   :: circ_class_n             !! Number of trees within each circumference class @tex $(m^{-2})$ @endtex
    REAL(r_std),        DIMENSION(:,:,:),                                INTENT(inout)   :: lai_per_level            !! LAI per vertical level @tex $(m^{2} m^{-2})$
    REAL(r_std),        DIMENSION(:,:,:,:),                              INTENT(inout)   :: z_array_out              !! An output of h_array, to use in sechiba
    TYPE(laieff_type),  DIMENSION(:,:,:),                                INTENT(inout)   :: laieff_fit               !! Fitted parameters for the effective LAI



    !! 0.4 Local variables


    REAL(r_std),        DIMENSION(kjpindex,nvm)                                          :: CrownVolumeL             !! Total Crown volume (m^3)
    REAL(r_std),        DIMENSION(kjpindex,nvm,nlevels_tot)                              :: PartialCrownVolumeL      !! Crown volume in each canopy layer (m^3)
    REAL(r_std),        DIMENSION(kjpindex,nvm,nlevels_tot)                              :: fitted_lai               !! fitted parameters for the LAI (-)
    REAL(r_std),        DIMENSION(ncirc)                                                 :: class_heights            !! Heights of the different circumference classes (m)
    REAL(r_std)                                                                          :: max_height               !! Maximum height of the canopy (m)
    REAL(r_std)                                                                          :: min_height               !! Starting height of the canopy (m)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                                          :: lstress_fac              !! Fraction of gaps n canopy (-)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                                          :: c0_alloc                 !! Root to sapwood tradeoff parameter (-)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                                          :: k_latosa_tmp             !! Leaf to sapwood area ratio (unitless)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                                          :: KF                       !! Scaling factor to convert sapwood mass into leaf mass (-)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                                          :: LF                       !! Scaling factor to convert sapwood mass into root mass (-)
    REAL(r_std)                                                                          :: B_sapwood                !! Sapwood biomass (gC/ind)
    REAL(r_std),        DIMENSION(kjpindex,nvm,ncirc,ndist_types)                        :: values                   !! An array which holds various canopy parameters (-)
    REAL(r_std),        DIMENSION(ncirc)                                                 :: cn_bottom                !! Height of the bottom of the canopy (m)
    REAL(r_std),        DIMENSION(kjpindex,nvm,nlevels_tot)                              :: z_levels                 !! Heights of each canopy levels (m)
    REAL(r_std)                                                                          :: upper_boundary           !! Upper boundary of the canopy level (m)
    REAL(r_std)                                                                          :: lower_boundary           !! Lower boundary of the canopy level (m)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                                          :: longevity_eff_root       !! Effective root longevity (d)
    REAL(r_std),        DIMENSION(kjpindex,nvm)                                          :: longevity_eff_sap        !! Effective sapwood longevity (d)
    REAL(r_std)                                                                          :: min_height_prescribed    !! Minimum canopy height prescribed (m)
    REAL(r_std),        DIMENSION(365)                                                   :: daily_lai_prescribed     !! Daily prescribed LAI (m^2/m^2)
    REAL(r_std)                                                                          :: lai_temp                 !! Temporary variable to calculate LAI
    INTEGER(i_std)                                                                       :: iday                     !! Day of Year
    INTEGER(i_std)                                                                       :: ipts,ivm, ilvl,icirc     !! Indives (resp: grid-cells, PFTs, canopy levels, circumference classes)
    REAL(r_std)                                                                          :: sla_est
    REAL(r_std), DIMENSION(kjpindex,nvm,nparts,nelements)                                :: tmp_bm

    !_ ================================================================================================================================


    !! 0- Initialisations

    CrownVolumeL(:,:)=0.
    PartialCrownVolumeL(:,:,:)=0.
    temp_growth(:) = 25.
    height(:,ibare_sechiba) = zero
    height_dom(:,ibare_sechiba) = zero
    lai(:,ibare_sechiba) = zero  
    fco2_flux_out(:) = zero
    fco2_lu_out(:) = zero
    fco2_wh_out(:) = zero
    fco2_ha_out(:) = zero
    lai(:,:)= 0.
    lai_per_level(:,:,:)= 0.
    laieff_fit(:,:,:)%a=zero
    laieff_fit(:,:,:)%b=zero
    laieff_fit(:,:,:)%c=zero
    laieff_fit(:,:,:)%d=zero
    laieff_fit(:,:,:)%e=zero
    circ_class_biomass(:,:,:,:,:) = 0.
    circ_class_n(:,:,:) = 0.
    z_array_out(:,:,:,:) = 0.
    height(:,:) = zero
    height_dom(:,:) = zero

    tmp_bm = cc_to_biomass(kjpindex,nvm,&
        circ_class_biomass(:,:,:,:,:),&
        circ_class_n(:,:,:))

    DO ipts=1,kjpindex
       DO ivm=1,nvm
          circ_class_n(ipts,ivm,:)=canopy_density_prescribed
       ENDDO
    ENDDO

    DO ivm=1,nvm
       longevity_eff_root(:,ivm) = longevity_root(ivm)*un
       longevity_eff_sap(:,ivm) = longevity_sap(ivm)*un
    ENDDO


    !! 1- Calculation of the daily LAI:
    !! According to the hemisphere of the site considered the calculation is a bit different but the principle stays the same. The LAI is set to its minimum value when the 
    !! DOY is in the "Low LAI period" (ie: winter in Northern Hemisphere) and set to it's maximum value in the "High LAI period" (ie: Summer in Northern Hemisphere). In the
    !! transition periods, (senescence_delta and phenology_delta), the LAI follows a constant increase (or decrease) until it reaches the end of the period. Both periods
    !! are arbitrary set to 50 days.

    daily_lai_prescribed(:) = minimum_lai
    IF (phenology_doy_start .LE. senescence_doy_end) THEN
       lai_temp = minimum_lai
       DO iday = 1, 365
          IF ((iday .GE. phenology_doy_start) .AND. (iday .LE. phenology_doy_start+phenology_delta)) THEN
             lai_temp= (maximum_lai-minimum_lai)/phenology_delta + lai_temp
             daily_lai_prescribed(iday) = lai_temp
          ELSEIF ((iday .GE. phenology_doy_start+phenology_delta) .AND. (iday .LE. senescence_doy_end-senescence_delta)) THEN
             daily_lai_prescribed(iday) = maximum_lai
          ELSEIF ((iday .GE. senescence_doy_end-senescence_delta) .AND. (iday .LE. senescence_doy_end)) THEN
             lai_temp = lai_temp - (maximum_lai-minimum_lai)/senescence_delta
             daily_lai_prescribed(iday) = lai_temp
          ELSE
             daily_lai_prescribed(iday) = minimum_lai
          ENDIF
       ENDDO
    ELSE
       lai_temp = maximum_lai
       DO iday = 1, 365
          IF ((iday .GE. senescence_doy_end-senescence_delta) .AND. (iday .LE. senescence_doy_end)) THEN
             lai_temp = lai_temp - (maximum_lai-minimum_lai)/senescence_delta
             daily_lai_prescribed(iday) = lai_temp
          ELSEIF ((iday .GE. senescence_doy_end) .AND. (iday .LE. phenology_doy_start)) THEN
             daily_lai_prescribed(iday) = minimum_lai
          ELSEIF ((iday .GE. phenology_doy_start) .AND. (iday .LE. phenology_doy_start+phenology_delta)) THEN
             lai_temp= (maximum_lai-minimum_lai)/phenology_delta + lai_temp
             daily_lai_prescribed(iday) = lai_temp
          ELSE
             daily_lai_prescribed(iday) = maximum_lai
          ENDIF
       ENDDO
    ENDIF


    DO ipts=1,kjpindex
       DO ivm=2,nvm
          IF(ALL(lvllai_prescribed(:) .LT. min_sechiba)) THEN
             IF (FLOOR((julian_diff+1)) .GE. 366) THEN
                lai(ipts,ivm)=daily_lai_prescribed(365)
             ELSE
                lai(ipts,ivm)=daily_lai_prescribed(FLOOR(julian_diff+1)) !! FLOOR returns the integer before (3.7 --> 3)
             ENDIF
          ELSE
             lai(ipts,ivm)=SUM(lvllai_prescribed(:))
          ENDIF
       ENDDO
    ENDDO


    veget(:,:) = veget_max(:,:)*(1-EXP(-lai(:,:)))
    
    !! 2- Calculations of the heights, biomass and LAI per level:
    !! Maximum height is prescribed. Minimum height is either prescribed or calculated thanks to the tree shape configured in stomate. Both cases lead to
    !! the caluclation of the canopy layers heights which have been decided to be of constant height (the interval max_height - min_height being divided 
    !! into nlevels_tot layers of constant height). 
    !! For the biomasses, the leaf biomass is calculated from the LAI and all the other biomasses are calcualated thanks to the allocation equations as in 
    !! stomate.f90. The only exception is the ratio above/below ground sapwod. It is here considered that 80% of the sapwood biomass is above ground.
    !! The final calculation is the LAI per level which is deduced from the tree shape (CrownVolume and PartialCrownVolume) as in stomate.f90.

    class_heights(:) = max_canopy_height_prescribed

    DO ipts=1 , kjpindex

       DO ivm = 2, nvm
          
          IF (veget_max(ipts,ivm) .GT. min_sechiba) THEN

             max_height = MAXVAL(class_heights(:))

             c0_alloc(ipts,ivm) = calculate_c0_alloc(ivm,longevity_eff_root(ipts,ivm), &
                  longevity_eff_sap(ipts,ivm))

             k_latosa_tmp(ipts,ivm) = k_latosa_max(ivm)

             IF (sla_dyn) THEN
                sla_est = biomass_to_lai(tmp_bm(ipts,ivm,ileaf,icarbon),ivm)/&
                        tmp_bm(ipts,ivm,ileaf,icarbon)
             ELSE
                sla_est = sla(ivm)
             ENDIF

             KF(ipts,ivm) = k_latosa_tmp(ipts,ivm) / &
                  (sla_est * pipe_density(ivm) * tree_ff(ivm))

             LF(ipts,ivm) = c0_alloc(ipts,ivm) * KF(ipts,ivm) 


             DO icirc = 1, ncirc

                !! Pour plusieurs classes de circonférences : trouver un moyen de splitter le LAI

                ! check for dynamic sla
                circ_class_biomass(ipts,ivm,icirc,ileaf,icarbon) = &
                        2.7*lai_to_biomass(lai(ipts,ivm),ivm)/circ_class_n(ipts,ivm,icirc)

                B_sapwood = circ_class_biomass(ipts,ivm,icirc,ileaf,icarbon)  * class_heights(icirc) / KF(ipts,ivm)

                circ_class_biomass(ipts,ivm,icirc,iroot,icarbon) = circ_class_biomass(ipts,ivm,icirc,ileaf,icarbon)  / LF(ipts,ivm)

                !! Calcul biomass heartabove en inversant wood_to_height
                IF (is_tree(ivm)) THEN
                   circ_class_biomass(ipts,ivm,icirc,isapabove,icarbon) = &
                        0.8*B_sapwood !! Mature forest
                   circ_class_biomass(ipts,ivm,icirc,isapbelow,icarbon) = &
                        0.2*B_sapwood
                   circ_class_biomass(ipts,ivm,icirc,iheartabove,icarbon)  = &
                        ((class_heights(icirc)**(pipe_tune3(ivm)/2+1))/(pipe_tune2(ipts,ivm)))** &
                        (2/pipe_tune3(ivm))*pi/4.*tree_ff(ivm)*pipe_density(ivm) &
                        - circ_class_biomass(ipts,ivm,icirc,isapabove,icarbon) 
                ELSE
                   ! grasses and crops have sapwood aboveground only. No
                   ! heartwood
                   circ_class_biomass(ipts,ivm,icirc,isapabove,icarbon) = &
                        B_sapwood
                   circ_class_biomass(ipts,ivm,icirc,isapbelow,icarbon) = zero
                   circ_class_biomass(ipts,ivm,icirc,iheartabove,icarbon) = zero

                ENDIF

             ENDDO


             IF (is_tree(ivm)) THEN

                ! for the stem height (aboveground only)
                values(ipts,ivm,:,iheight) = &
                     wood_to_height(circ_class_biomass(ipts,ivm,:,:,icarbon),ivm,pipe_tune2(ipts,ivm))

                ! for the stem diameter
                values(ipts,ivm,:,idiameter) = &
                     wood_to_dia(circ_class_biomass(ipts,ivm,:,:,icarbon),ivm,pipe_tune2(ipts,ivm))

                ! vertical and horizontal crown diameter
                CALL wood_to_cn_dia(circ_class_biomass(ipts,ivm,:,:,icarbon), &
                     circ_class_n(ipts,ivm,:),ivm,values(ipts,ivm,:,icndiahor), &
                     values(ipts,ivm,:,icndiaver),pipe_tune2(ipts,ivm))

                ! for the crown area
                values(ipts,ivm,:,icnarea) = &
                     pi/4*values(ipts,ivm,:,icndiahor)**2

                ! for the crown volume
                values(ipts,ivm,:,icnvol) = &
                     wood_to_cv(circ_class_biomass(ipts,ivm,:,:,icarbon),&
                     circ_class_n(ipts,ivm,:),ivm,pipe_tune2(ipts,ivm))

                cn_bottom(:)=class_heights(:)-values(ipts,ivm,:,icndiaver)

                IF(min_canopy_height_prescribed .NE. 0.) THEN
                   min_height=min_canopy_height_prescribed
                ELSE
                   min_height=MINVAL(cn_bottom(:))
                ENDIF


                z_levels(ipts,ivm,1) = min_height

                DO ilvl=1,nlevels_tot-1

                   z_levels(ipts,ivm,ilvl+1) = min_height + ilvl*(max_height - min_height)/(nlevels_tot-1)

                ENDDO


                !! The only place where z_array_out is udes is mleb.f90. In this module, only the first circumference class is 
                !! considered. However, in a matter of consistency, the lowest level of the canopy is maybe not in this class.
                !! Consequetnly, the canopy is divided into nlevels_tot leveles from the lowest level to the highest one and 
                !! this value is given to z_array_out.

                DO icirc=1,ncirc
                   z_array_out(ipts,ivm,icirc,:) = z_levels(ipts,ivm,:)
                ENDDO

                IF (ALL(lvllai_prescribed(:) .LT. min_sechiba)) THEN

                   DO icirc=1,ncirc

                      CrownVolumeL(ipts,ivm)=CrownVolumeL(ipts,ivm)+&
                           values(ipts,ivm,icirc,icnvol)*&
                           circ_class_n(ipts,ivm,icirc)/SUM(circ_class_n(ipts,ivm,:))

                      DO ilvl=1,nlevels_tot

                         IF (ilvl .NE. nlevels_tot) THEN
                            upper_boundary=z_array_out(ipts,ivm,icirc,ilvl+1)
                         ELSE
                            upper_boundary=z_array_out(ipts,ivm,icirc,ilvl)+max_height
                         ENDIF

                         lower_boundary=z_array_out(ipts,ivm,icirc,ilvl)

                         PartialCrownVolumeL(ipts,ivm,ilvl)=PartialCrownVolumeL(ipts,ivm,ilvl) + &
                              partial_spheroid_vol(upper_boundary, lower_boundary, &
                              values(ipts,ivm,icirc,icndiaver)/deux, &
                              values(ipts,ivm,icirc,icndiahor)/deux, &
                              class_heights(icirc))*circ_class_n(ipts,ivm,icirc)/SUM(circ_class_n(ipts,ivm,:))

                      ENDDO

                   ENDDO

                   DO ilvl=1,nlevels_tot

                      lai_per_level(ipts,ivm,ilvl) = MAX(0.01,lai(ipts,ivm)*PartialCrownVolumeL(ipts,ivm,ilvl)/CrownVolumeL(ipts,ivm))

                   ENDDO


                ELSE

                   lai_per_level(ipts,ivm,:) = lvllai_prescribed(:)

                ENDIF

             ELSE

                min_height=min_canopy_height_prescribed

                z_levels(ipts,ivm,1) = min_height

                DO ilvl=1,nlevels_tot-1

                   z_levels(ipts,ivm,ilvl+1) = min_height + ilvl*(max_height - min_height)/(nlevels_tot-1)

                ENDDO

                DO icirc=1,ncirc
                   z_array_out(ipts,ivm,icirc,:) = z_levels(ipts,ivm,:)
                ENDDO


                IF (ALL(lvllai_prescribed(:) .LT. min_sechiba)) THEN

                   lai_per_level(ipts,ivm,:) = MAX(0.01, lai(ipts,ivm)/nlevels_tot)

                ELSE

                   lai_per_level(ipts,ivm,:) = lvllai_prescribed(:)

                ENDIF


             ENDIF



             !! This is a highly critical point to discuss about. The fitting_laieff function is not functioning here.
             !! Thus, the effective LAI (accroding to the time of the day) is arbitrary given a value between the lai_per_level 
             !! and 2*lai_per_level when sun is low in the sky. This behaviour as no physical background but tries to mimic the 
             !! behaviour of the fitting_lai function. This mainly has an impact in "albedo_surface" and, thus, "mleb".


!!$             CALL fitting_laieff(kjpindex, z_array_out, circ_class_biomass, & 
!!$                  circ_class_n, veget_max, lai_per_level, laieff_fit)

             fitted_lai(ipts,ivm,:) = lai_per_level(ipts,ivm,:)*(1.0-coszang(ipts))+ lai_per_level(ipts,ivm,:)

             laieff_fit(ipts,ivm,:)%a=fitted_lai(ipts,ivm,:)
             laieff_fit(ipts,ivm,:)%b=zero
             laieff_fit(ipts,ivm,:)%c=zero
             laieff_fit(ipts,ivm,:)%d=zero
             laieff_fit(ipts,ivm,:)%e=zero

          ENDIF
       ENDDO

    ENDDO



    !! Last calculations needed in sechiba.

    qsintmax(:,:) = qsintcst * veget(:,:) * lai(:,:)
    qsintmax(:,1) = zero 
    DO ipts=1,kjpindex
       DO ivm=1,nvm
          IF (veget_max(ipts,ivm) .EQ. zero) THEN
             height_dom(ipts,ivm) = zero
          ELSE
             IF(is_tree(ivm))THEN
                height_dom(ipts,ivm)=MAXVAL(class_heights(:))
             ELSE
                height_dom(ipts,ivm)=height(ipts,ivm)
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    max_height_store(:,:) = height_dom(:,:)


  END SUBROUTINE slowproc_impose_site_level_lai

!! ================================================================================================================================
!! SUBROUTINE 	: get_soilcorr_zobler
!!
!>\BRIEF         The "get_soilcorr" routine defines the table of correspondence
!!               between the Zobler types and the three texture types known by SECHIBA and STOMATE :
!!               silt, sand and clay. 
!!
!! DESCRIPTION : get_soilcorr is needed if you use soils_param.nc .\n
!!               The data from this file is then interpolated to the grid of the model. \n
!!               The aim is to get fractions for sand loam and clay in each grid box.\n
!!               This information is used for soil hydrology and respiration.
!!
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S) : ::texfrac_table
!!
!! REFERENCE(S)	: 
!! - Zobler L., 1986, A World Soil File for global climate modelling. NASA Technical memorandum 87802. NASA 
!!   Goddard Institute for Space Studies, New York, U.S.A.
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================

  SUBROUTINE get_soilcorr_zobler (nzobler,textfrac_table)

    IMPLICIT NONE

    !! 0. Variables and parameters declaration
    
    INTEGER(i_std),PARAMETER :: nbtypes_zobler = 7                    !! Number of Zobler types (unitless)

    !! 0.1  Input variables
    
    INTEGER(i_std),INTENT(in) :: nzobler                              !! Size of the array (unitless)
    
    !! 0.2 Output variables 
    
    REAL(r_std),DIMENSION(nzobler,ntext),INTENT(out) :: textfrac_table !! Table of correspondence between soil texture class
                                                                       !! and granulometric composition (0-1, unitless)
    
    !! 0.4 Local variables
    
    INTEGER(i_std) :: ib                                              !! Indice (unitless)
    
!_ ================================================================================================================================

    !-
    ! 0. Check consistency
    !-  
    IF (nzobler /= nbtypes_zobler) THEN 
       CALL ipslerr_p(3,'get_soilcorr', 'nzobler /= nbtypes_zobler',&
          &   'We do not have the correct number of classes', &
          &                 ' in the code for the file.')  ! Fatal error
    ENDIF

    !-
    ! 1. Textural fraction for : silt        sand         clay
    !-
    textfrac_table(1,:) = (/ 0.12, 0.82, 0.06 /)
    textfrac_table(2,:) = (/ 0.32, 0.58, 0.10 /)
    textfrac_table(3,:) = (/ 0.39, 0.43, 0.18 /)
    textfrac_table(4,:) = (/ 0.15, 0.58, 0.27 /)
    textfrac_table(5,:) = (/ 0.34, 0.32, 0.34 /)
    textfrac_table(6,:) = (/ 0.00, 1.00, 0.00 /)
    textfrac_table(7,:) = (/ 0.39, 0.43, 0.18 /)


    !-
    ! 2. Check the mapping for the Zobler types which are going into the ORCHIDEE textures classes 
    !-
    DO ib=1,nzobler ! Loop over # classes soil
       
       IF (ABS(SUM(textfrac_table(ib,:))-1.0) > EPSILON(un)) THEN ! The sum of the textural fractions should not exceed 1 !
          WRITE(numout,*) &
               &     'Error in the correspondence table', &
               &     ' sum is not equal to 1 in', ib
          WRITE(numout,*) textfrac_table(ib,:)
          CALL ipslerr_p(3,'get_soilcorr', 'SUM(textfrac_table(ib,:)) /= 1.0',&
               &                 '', 'Error in the correspondence table') ! Fatal error
       ENDIF
       
    ENDDO ! Loop over # classes soil

    
  END SUBROUTINE get_soilcorr_zobler

!! ================================================================================================================================
!! SUBROUTINE 	: get_soilcorr_usda
!!
!>\BRIEF         The "get_soilcorr_usda" routine defines the table of correspondence
!!               between the 12 USDA textural classes and their granulometric composition, 
!!               as % of silt, sand and clay. This is used to further defien clayfraction.
!!
!! DESCRIPTION : get_soilcorr is needed if you use soils_param.nc .\n
!!               The data from this file is then interpolated to the grid of the model. \n
!!               The aim is to get fractions for sand loam and clay in each grid box.\n
!!               This information is used for soil hydrology and respiration.
!!               The default map in this case is derived from Reynolds et al 2000, \n
!!               at the 1/12deg resolution, with indices that are consistent with the \n
!!               textures tabulated below
!!
!! RECENT CHANGE(S): Created by A. Ducharne on July 02, 2014
!!
!! MAIN OUTPUT VARIABLE(S) : ::texfrac_table
!!
!! REFERENCE(S)	: 
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================

  SUBROUTINE get_soilcorr_usda (nusda,textfrac_table)

    IMPLICIT NONE

    !! 0. Variables and parameters declaration
    
    !! 0.1  Input variables
    
    INTEGER(i_std),INTENT(in) :: nusda                               !! Size of the array (unitless)
    
    !! 0.2 Output variables 
    
    REAL(r_std),DIMENSION(nusda,ntext),INTENT(out) :: textfrac_table !! Table of correspondence between soil texture class
                                                                     !! and granulometric composition (0-1, unitless)
    
    !! 0.4 Local variables

    INTEGER(i_std),PARAMETER :: nbtypes_usda = 13                    !! Number of USDA texture classes (unitless)
    INTEGER(i_std) :: n                                              !! Index (unitless)
    
!_ ================================================================================================================================

    !-
    ! 0. Check consistency
    !-  
    IF (nusda /= nbtypes_usda) THEN 
       CALL ipslerr_p(3,'get_soilcorr', 'nusda /= nbtypes_usda',&
          &   'We do not have the correct number of classes', &
          &                 ' in the code for the file.')  ! Fatal error
    ENDIF

    !! Parameters for soil type distribution :
    !! Sand, Loamy Sand, Sandy Loam, Silt Loam, Silt, Loam, Sandy Clay Loam, Silty Clay Loam, Clay Loam, Sandy Clay, Silty Clay, Clay
    ! The order comes from constantes_soil.f90
    ! The corresponding granulometric composition comes from Carsel & Parrish, 1988

    !-
    ! 1. Textural fractions for : sand, clay
    !-
    textfrac_table(1,2:3)  = (/ 0.93, 0.03 /) ! Sand
    textfrac_table(2,2:3)  = (/ 0.81, 0.06 /) ! Loamy Sand
    textfrac_table(3,2:3)  = (/ 0.63, 0.11 /) ! Sandy Loam
    textfrac_table(4,2:3)  = (/ 0.17, 0.19 /) ! Silt Loam
    textfrac_table(5,2:3)  = (/ 0.06, 0.10 /) ! Silt
    textfrac_table(6,2:3)  = (/ 0.40, 0.20 /) ! Loam
    textfrac_table(7,2:3)  = (/ 0.54, 0.27 /) ! Sandy Clay Loam
    textfrac_table(8,2:3)  = (/ 0.08, 0.33 /) ! Silty Clay Loam
    textfrac_table(9,2:3)  = (/ 0.30, 0.33 /) ! Clay Loam
    textfrac_table(10,2:3) = (/ 0.48, 0.41 /) ! Sandy Clay
    textfrac_table(11,2:3) = (/ 0.06, 0.46 /) ! Silty Clay
    textfrac_table(12,2:3) = (/ 0.15, 0.55 /) ! Clay
    textfrac_table(13,2:3) = (/ 0.15, 0.55 /) ! Clay

    ! Fraction of silt

    DO n=1,nusda
       textfrac_table(n,1) = 1. - textfrac_table(n,2) - textfrac_table(n,3)
    END DO
       
  END SUBROUTINE get_soilcorr_usda

!! ================================================================================================================================
!! FUNCTION 	: tempfunc
!!
!>\BRIEF        ! This function interpolates value between ztempmin and ztempmax
!! used for lai detection. 
!!
!! DESCRIPTION   : This subroutine calculates a scalar between 0 and 1 with the following equation :\n
!!                 \latexonly
!!                 \input{constantes_veg_tempfunc.tex}
!!                 \endlatexonly
!!
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : tempfunc_result
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  FUNCTION tempfunc (temp_in) RESULT (tempfunc_result)


    !! 0. Variables and parameters declaration

    REAL(r_std),PARAMETER    :: ztempmin=273._r_std   !! Temperature for laimin (K)
    REAL(r_std),PARAMETER    :: ztempmax=293._r_std   !! Temperature for laimax (K)
    REAL(r_std)              :: zfacteur              !! Interpolation factor   (K^{-2})

    !! 0.1 Input variables

    REAL(r_std),INTENT(in)   :: temp_in               !! Temperature (K)

    !! 0.2 Result

    REAL(r_std)              :: tempfunc_result       !! (unitless)
    
!_ ================================================================================================================================

    !! 1. Define a coefficient
    zfacteur = un/(ztempmax-ztempmin)**2
    
    !! 2. Computes tempfunc
    IF     (temp_in > ztempmax) THEN
       tempfunc_result = un
    ELSEIF (temp_in < ztempmin) THEN
       tempfunc_result = zero
    ELSE
       tempfunc_result = un-zfacteur*(ztempmax-temp_in)**2
    ENDIF !(temp_in > ztempmax)


  END FUNCTION tempfunc
  
  !! ================================================================================================================================
  !! SUBROUTINE   : slowproc_readirrigmap_dyn
  !!
  !>\BRIEF        Function to interpolate irrigation maps
  !!
  !! DESCRIPTION  : This function interpolates the irrigation maps from original resolution to simul. resolution
  !!
  !! RECENT CHANGE(S): None
  !!
  !! MAIN OUTPUT VARIABLE(S): :: irrigmap_new
  !!
  !! REFERENCE(S) : None
  !!
  !! FLOWCHART    : None
  !! \n
  !_ ================================================================================================================================

  SUBROUTINE slowproc_readirrigmap_dyn(nbpt, lalo, neighbours,  resolution, contfrac,         &
       irrigmap_new)

    USE interpweight

    IMPLICIT NONE

    !  0.1 INPUT
    !
    INTEGER(i_std), INTENT(in)                             :: nbpt            !! Number of points for which the data needs
                                                                              !! to be interpolated
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in)             :: lalo            !! Vector of latitude and longitudes (beware of the order !)
    INTEGER(i_std), DIMENSION(nbpt,NbNeighb), INTENT(in)   :: neighbours      !! Vector of neighbours for each grid point
                                                                              !! (1=North and then clockwise)
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in)             :: resolution      !! The size in km of each grid-box in X and Y
    REAL(r_std), DIMENSION(nbpt), INTENT(in)               :: contfrac        !! Fraction of continent in the grid
    !
    !  0.2 OUTPUT
    !
    REAL(r_std), DIMENSION(nbpt), INTENT(out)          :: irrigmap_new       !! new irrigation map in m2 per grid cell
    !
    !  0.3 LOCAL
    !
    !
    CHARACTER(LEN=80) :: filename
    INTEGER(i_std) :: ib
    !
    ! for irrigated_new case :

    REAL(r_std), DIMENSION(nbpt)                         :: irrigref_frac    !! irrigation fractions re-dimensioned
    REAL(r_std), DIMENSION(nbpt)                         :: airrig           !! Availability of the soilcol interpolation
    REAL(r_std)                          :: vmin, vmax       !! min/max values to use for the renormalization
    CHARACTER(LEN=80)                                    :: variablename     !! Variable to interpolate
    CHARACTER(LEN=80)                                    :: lonname, latname !! lon, lat names in input file
    REAL(r_std), DIMENSION(nvm)                          :: variabletypevals !! Values for all the types of the variable
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
    REAL(r_std), DIMENSION(3)                            :: maskvals         !! values to use to mask (according to
                                                                             !!   `maskingtype')
    CHARACTER(LEN=250)                                   :: namemaskvar      !! name of the variable to use to mask
    CHARACTER(LEN=250)                                   :: msg

  !_ ================================================================================================================================

    IF (printlev_loc >= 5) PRINT *,'  In slowproc_readirrigmap_dyn'

    !
    !Config Key   = IRRIGATION_DYN_FILE
    !Config Desc  = Name of file from which the DYNAMIC irrigation fraction map is to be read
    !Config Def   = IRRIGmap.nc
    !Config If    = IRRIG_DYN
    !Config Help  = The name of the file to be opened to read an irrigation
    !Config         map is to be given here.
    !Config Units = [FILE]
    !
    filename = 'IRRIGmap.nc'
    CALL getin_p('IRRIGATION_DYN_FILE',filename)
    variablename = 'irrig'

    IF (xios_interpolation) THEN
       IF (printlev_loc >= 1) WRITE(numout,*) "slowproc_readirrigmap_dyn: Use XIOS to read and interpolate " &
            // TRIM(filename) // " for variable " // TRIM(variablename)

       CALL xios_orchidee_recv_field('irrigmap_interp',irrigref_frac)

       airrig = 1.0
    ELSE

       IF (printlev_loc >= 1) WRITE(numout,*) "slowproc_readirrigmap_dyn: Read and interpolate " &
            // TRIM(filename) // " for variable " // TRIM(variablename)

       ! Assigning values to vmin, vmax
       vmin = 1.
       vmax = 1.

       variabletypevals = -un

       !! Variables for interpweight
       ! Type of calculation of cell fractions
       fractype = 'default'
       ! Name of the longitude and latitude in the input file
       lonname = 'lon'
       latname = 'lat'
       ! Should negative values be set to zero from input file?
       nonegative = .TRUE.
       ! Type of mask to apply to the input data (see header for more details)
       maskingtype = 'nomask'
       ! Values to use for the masking
       maskvals = (/ 0.05, 0.05 , un /)
       ! Name of the variable with the values for the mask in the input file (only if maskkingtype='var') (here not used)
       namemaskvar = ''

       CALL interpweight_2Dcont(nbpt, 0, 0, lalo, resolution, neighbours,        &
	 contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype,        &
	 maskvals, namemaskvar, -1, fractype, 0., 0.,                                 &
	 irrigref_frac, airrig)

       IF (printlev_loc >= 5) WRITE(numout,*)'  slowproc_readirrigmap_dyn after interpweight_2Dcont'

    ENDIF
    
    ! Write diagnostics
    CALL xios_orchidee_send_field("interp_diag_irrigmap_frac",irrigref_frac)

    IF (printlev_loc >= 5) THEN
      WRITE(numout,*)'  slowproc_readirrigmap_dyn before updating loop nbpt:', nbpt
    END IF
    
    ! Transform form percentage (irrigref_frac) to area (irrigmap_new)
    irrigmap_new(:) = zero
    irrigref_frac(:) = irrigref_frac(:)/100.

    DO ib=1,nbpt
      IF (irrigref_frac(ib) < 1. .AND. irrigref_frac(ib) > 0.005 ) THEN
        irrigmap_new(ib) = irrigref_frac(ib) * area(ib) * contfrac(ib)
      ELSEIF (irrigref_frac(ib) > 1.) THEN
        irrigmap_new(ib) = area(ib)*contfrac(ib)
      !ELSE  THEN irrigated_new= zero. Already set, put to lisibility
      ENDIF
    ENDDO

    ! Write diagnostics
    !CALL xios_orchidee_send_field("airrig",airrig)

    IF (printlev_loc >= 3) WRITE(numout,*) '  slowproc_readirrigmap_dyn ended'

  END SUBROUTINE slowproc_readirrigmap_dyn

  !! ================================================================================================================================
  !! SUBROUTINE   : slowproc_read_aeisw_map
  !!
  !>\BRIEF        Function to interpolate irrigation maps
  !!
  !! DESCRIPTION  : This function interpolates the maps of areas equipped for irrigation with surface water
  !!                from original resolution to simul. resolution. This is used in the new irrigation scheme, that
  !!                restrain water availability according to environmental needs and type of equippment to irrigate.
  !!
  !! RECENT CHANGE(S): None
  !!
  !! MAIN OUTPUT VARIABLE(S): :: fraction_aeirrig_sw
  !!
  !! REFERENCE(S) : None
  !!
  !! FLOWCHART    : None
  !! \n
  !_ ================================================================================================================================


  SUBROUTINE slowproc_read_aeisw_map(nbpt, lalo, neighbours,  resolution, contfrac,         &
     fraction_aeirrig_sw)

    USE interpweight

    IMPLICIT NONE

    !  0.1 INPUT
    !
    INTEGER(i_std), INTENT(in)                             :: nbpt            !! Number of points for which the data needs
                                                                              !! to be interpolated
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in)             :: lalo            !! Vector of latitude and longitudes (beware of the order !)
    INTEGER(i_std), DIMENSION(nbpt,NbNeighb), INTENT(in)   :: neighbours      !! Vector of neighbours for each grid point
                                                                              !! (1=North and then clockwise)
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in)             :: resolution      !! The size in km of each grid-box in X and Y
    REAL(r_std), DIMENSION(nbpt), INTENT(in)               :: contfrac        !! Fraction of continent in the grid
    !
    !  0.2 OUTPUT
    !
    REAL(r_std), DIMENSION(nbpt), INTENT(out)          :: fraction_aeirrig_sw       !! Fraction of area equipped for irrigation from surface water, of irrig_frac
                                                                                ! 1.0 here corresponds to fraction of irrig. area, not grid cell
    !
    !  0.3 LOCAL
    CHARACTER(LEN=80) :: filename
    INTEGER(i_std) :: ib
    !
    ! for irrigated_new case :

    REAL(r_std), DIMENSION(nbpt)                         :: irrigref_frac    !! irrigation fractions re-dimensioned (0-100, percentage)
    REAL(r_std), DIMENSION(nbpt)                         :: airrig           !! Availability of the soilcol interpolation
    REAL(r_std)                          :: vmin, vmax       !! min/max values to use for the renormalization
    CHARACTER(LEN=80)                                    :: variablename     !! Variable to interpolate
    CHARACTER(LEN=80)                                    :: lonname, latname !! lon, lat names in input file
    REAL(r_std), DIMENSION(nvm)                          :: variabletypevals !! Values for all the types of the variable
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
    REAL(r_std), DIMENSION(3)                            :: maskvals         !! values to use to mask (according to
                                                                             !!   `maskingtype')
    CHARACTER(LEN=250)                                   :: namemaskvar      !! name of the variable to use to mask
    CHARACTER(LEN=250)                                   :: msg

  !_ ================================================================================================================================

    IF (printlev_loc >= 5) PRINT *,'  In slowproc_read_aeisw_map'

    !
    !Config Key   = FRACTION_AEI_SW_FILE
    !Config Desc  = Name of file with AEI with SW
    !Config Def   = AEI_SW_pct.nc
    !Config If    = SELECT_SOURCE_IRRIG
    !Config Help  = The name of the file to be opened to read an AIE from SW
    !Config         map is to be given here.
    !Config Units = [FILE]
    !
    filename = 'AEI_SW_pct.nc'
    CALL getin_p('FRACTION_AEI_SW_FILE',filename)
    variablename = 'aeisw_pct'

    IF (xios_interpolation) THEN
       IF (printlev_loc >= 1) WRITE(numout,*) "slowproc_read_aeisw_map: Use XIOS to read and interpolate " &
            // TRIM(filename) // " for variable " // TRIM(variablename)

       CALL xios_orchidee_recv_field('aei_sw_interp',irrigref_frac)

       airrig = 1.0
    ELSE

       IF (printlev_loc >= 1) WRITE(numout,*) "slowproc_read_aeisw_map: Read and interpolate " &
            // TRIM(filename) // " for variable " // TRIM(variablename)

       ! Assigning values to vmin, vmax
       vmin = 1.
       vmax = 1.

       variabletypevals = -un

       !! Variables for interpweight
       ! Type of calculation of cell fractions
       fractype = 'default'
       ! Name of the longitude and latitude in the input file
       lonname = 'lon'
       latname = 'lat'
       ! Should negative values be set to zero from input file?
       nonegative = .TRUE.
       ! Type of mask to apply to the input data (see header for more details)
       maskingtype = 'mabove'
       ! Values to use for the masking
       maskvals = (/ 0.05, 0.05 , un /)
       ! Name of the variable with the values for the mask in the input file (only if maskkingtype='var') (here not used)
       namemaskvar = ''

       CALL interpweight_2Dcont(nbpt, 0, 0, lalo, resolution, neighbours,        &
	 contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype,        &
	 maskvals, namemaskvar, -1, fractype, aei_sw_default, aei_sw_default,                                 &
	 irrigref_frac, airrig)

       IF (printlev_loc >= 5) WRITE(numout,*)'  slowproc_read_aeisw_map after interpweight_2Dcont'
       
    END IF
    
    ! Write diagnostics
    CALL xios_orchidee_send_field("interp_diag_aeiswmap_frac",irrigref_frac) ! irrigref_frac here has a different meaning than in subroutine slowproc_readirrigmap_dyn

    IF (printlev_loc >= 5) THEN
      WRITE(numout,*)'  slowproc_read_aeisw_map before updating loop nbpt:', nbpt
    END IF
    fraction_aeirrig_sw(:) = irrigref_frac(:)/100.
    ! Write diagnostics
    !CALL xios_orchidee_send_field("airrig",airrig)
    DO ib=1,nbpt

      fraction_aeirrig_sw(ib) = MIN(fraction_aeirrig_sw(ib), 0.99 )
      fraction_aeirrig_sw(ib) = MAX(fraction_aeirrig_sw(ib), 0.01 )

    ENDDO


    IF (printlev_loc >= 3) WRITE(numout,*) 'slowproc_read_aeisw_map ended'

  END SUBROUTINE slowproc_read_aeisw_map

END MODULE slowproc
