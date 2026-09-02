! ===============================================================================================================================
! MODULE       : albedo_surface
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        Calculate the surface albedo using a variety of different schemes
!!
!! \n DESCRIPTION : This module includes the mechanisms for
!! 1. :: albedo_surface_main computes the two_stream approach of Pinty et al 2006.
!!           Requires an effective leaf area index and depends on the solar angle.
!!           Separates light into NIR and VIS, and diffuse and direct illumination
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCES(S)    : 
!!                    Chalita, S. and H Le Treut (1994), The albedo of temperate
!!                    and boreal forest and the Northern Hemisphere climate: a 
!!                    sensitivity experiment using the LMD GCM, 
!!                    Climate Dynamics, 10 231-240.
!!
!!                    B. Pinty, T. Lavergne, R.E. Dickinson, J-L. Widlowski, 
!!                    N. Gobron and M. M. Verstraete (2006). Simplifying the 
!!                    Interaction of Land Surfaces with Radiation for Relating 
!!                    Remote Sensing Products to Climate Models. Journal of 
!!                    Geophysical Research. Vol 111, D02116.
!!
!! SVN              :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_sechiba/albedo_surface.f90 $
!! $Date: 2026-04-02 12:01:32 +0200 (jeu. 02 avril 2026) $
!! $Revision: 9457 $
!! \n
!_ ================================================================================================================================

MODULE albedo_surface
  !
  USE ioipsl
  !
  ! modules used :
  USE constantes
  USE constantes_soil
  USE pft_parameters
  USE structures
  USE interpol_help
  USE ioipsl_para, ONLY : ipslerr_p
  USE stomate_laieff
  USE xios_orchidee
  USE grid
  USE time, ONLY : sec_start, sec_end
  USE function_library,    ONLY: get_printlev
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: albedo_surface_main, albedo_surface_initialize, albedo_surface_finalize, &
       albedo_surface_clear, albedo_surface_xios_initialize
  
  INTEGER, SAVE                                   :: printlev_loc   !! Local level of text output for current module
!$OMP THREADPRIVATE(printlev_loc)
  LOGICAL, SAVE                                   :: firstcall_albedo_surface=.TRUE. !! Initialization phase of this module
!$OMP THREADPRIVATE(firstcall_albedo_surface)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)  :: soilalb_dry    !! Albedo values for the dry bare soil (unitless)
!$OMP THREADPRIVATE(soilalb_dry)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)  :: soilalb_wet    !! Albedo values for the wet bare soil (unitless)
!$OMP THREADPRIVATE(soilalb_wet)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)  :: soilalb_moy    !! Albedo values for the mean bare soil (unitless)
!$OMP THREADPRIVATE(soilalb_moy)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)  :: bckgrnd_alb    !! Albedo values for the background bare soil (unitless)
!$OMP THREADPRIVATE(bckgrnd_alb)


  LOGICAL, SAVE                                   :: test_alb_rt_vegmax
!$OMP THREADPRIVATE(test_alb_rt_vegmax)

  LOGICAL, SAVE                                   :: test_alb_rt_bgsnowbare ! This is a new test that can be applied on alb1 or standard(alb2) cases
!$OMP THREADPRIVATE(test_alb_rt_bgsnowbare)


CONTAINS

  !!  =============================================================================================================================
  !! SUBROUTINE:    albedo_surface_xios_initialize
  !!
  !>\BRIEF	  Initialize xios dependant defintion before closing context defintion
  !!
  !! DESCRIPTION:	  Initialize xios dependant defintion needed for the interpolations done in albedo.
  !!                      Reading is deactivated if the sechiba restart file exists because the variable
  !!                      should be in the restart file already.
  !!                      This subruting is called before closing context with xios_orchidee_close_definition in intersurf 
  !!                      via the subroutine sechiba_xios_initialize. 
  !!
  !! \n
  !_ ==============================================================================================================================

  SUBROUTINE albedo_surface_xios_initialize
    
    CHARACTER(LEN=255)   :: filename        !! Filename read from run.def
    CHARACTER(LEN=255)   :: name            !! Filename without suffix .nc
    LOGICAL              :: lerr            !! Flag to dectect error
     
 
    ! Read the file name for the albedo input file from run.def
    filename = 'alb_bg.nc'
    CALL getin_p('ALB_BG_FILE',filename)
    
    ! Remove suffix .nc from filename
    name = filename(1:LEN_TRIM(FILENAME)-3)
    
    ! Define the file name in XIOS
    CALL xios_orchidee_set_file_attr("albedo_file",name=name)

    ! Define default values for albedo
    lerr=xios_orchidee_setvar('albbg_vis_default',0.129)
    lerr=xios_orchidee_setvar('albbg_nir_default',0.247)

    ! Check if the albedo file will be read by XIOS, by IOIPSL or not at all
    IF (xios_interpolation .AND. alb_bg_modis .AND. (restname_in=='NONE' .OR. alb_bg_modis_reinit)) THEN
       ! The albedo file will be read using XIOS
       IF (printlev>=2) WRITE(numout,*) 'Reading of albedo file will be done later using XIOS. The filename is ', filename

    ELSE
       ! The albedo file will not be read or it will be read by IOIPSL
       IF (.NOT. alb_bg_modis) THEN
          IF (printlev>=2) WRITE (numout,*) 'No reading of albedo will be done because alb_bg_modis=FALSE'
       ELSE IF (restname_in=='NONE' .OR. alb_bg_modis_reinit) THEN
          IF (printlev>=2) WRITE (numout,*) 'The albedo file will be read later by IOIPSL'
       ELSE
          IF (printlev>=2) WRITE (numout,*) 'The albedo file will not be read because the restart file exists.'
       END IF

       ! The albedo file will not be read by XIOS. Now deactivate albedo for XIOS.
       CALL xios_orchidee_set_file_attr("albedo_file",enabled=.FALSE.)
       CALL xios_orchidee_set_field_attr("bg_alb_vis_interp",enabled=.FALSE.)
       CALL xios_orchidee_set_field_attr("bg_alb_nir_interp",enabled=.FALSE.)
    END IF

    ! Get the file name from run.def file and set file attributes accordingly
    filename = 'soils_param.nc'
    CALL getin_p('SOILALB_FILE',filename)
    name = filename(1:LEN_TRIM(FILENAME)-3)
    CALL xios_orchidee_set_file_attr("soilcolor_file",name=name)

    ! Determine if reading of variable soilcolor will be done with XIOS.
    ! If not, deactivate the file and fieldgroup.  
    IF (xios_interpolation .AND. restname_in=='NONE' .AND. .NOT. alb_bg_modis) THEN
       ! Reading will be done with XIOS later
       IF (printlev>=2) WRITE(numout,*) 'Reading of variable soilcolor will be done later using XIOS. The filename is ', filename
    ELSE
       ! No reading, deactivate soilcolor_file and corresponding field group
       IF (printlev>=2) WRITE(numout,*) 'Reading of variable soilcolor will not be done with XIOS.'
       CALL xios_orchidee_set_file_attr("soilcolor_file",enabled=.FALSE.)
       CALL xios_orchidee_set_fieldgroup_attr("soil_color",enabled=.FALSE.)
    END IF
    
  END SUBROUTINE albedo_surface_xios_initialize

!!  =============================================================================================================================
!! SUBROUTINE  :  albedo_surface_initialize
!!
!>\BRIEF          Initialize albedo module
!!
!! DESCRIPTION :  Allocate module variables, read from restart file or initialize with default values.
!!
!! MAIN OUTPUT VARIABLE(S)
!!
!! REFERENCE(S)			            : None
!! 
!! \n
!_ ==============================================================================================================================
  !+++CHECK+++
  ! Variable dimensions xxx_Tot_mean are for a single pixel. Causing 1+1 problems
  ! when running a larger domain. Needs to be corrected when implementing
  ! a global use of the multi-layer energy budget.
  SUBROUTINE albedo_surface_initialize (kjit,  kjpindex, rest_id, &
                                lalo, neighbours, resolution, contfrac, &
                                drysoil_frac, veget_max, coszang, frac_nobio, & 
                                snow, snow_age, &
                                snow_nobio_age, frac_snow_veg, frac_snow_nobio,&
                                z0m, laieff_fit, &
                                albedo, &
                                Isotrop_Abs_Tot_p, Isotrop_Tran_Tot_p, &
                                Light_Abs_Tot_mean, Light_Alb_Tot_mean, &
                                laieff_isotrop, veget)
    !+++++++++++
    
    !! 0. Variable and parameter declaration
    !! 0.1 Input variables  
    INTEGER(i_std), INTENT(in)                               :: kjit                   !! Time step number 
    INTEGER(i_std), INTENT(in)                               :: kjpindex               !! Domain size
    INTEGER(i_std), INTENT(in)                               :: rest_id                !! _Restart_ file identifier
    REAL(r_std), DIMENSION(kjpindex,2), INTENT (in)          :: lalo                   !! Geographical coordinates
    INTEGER(i_std), DIMENSION(kjpindex,NbNeighb), INTENT(in) :: neighbours             !! neighoring grid points if land
    REAL(r_std), DIMENSION(kjpindex,2), INTENT(in)           :: resolution             !! size in x an y of the grid (m)
    REAL(r_std), DIMENSION(kjpindex), INTENT(in)             :: contfrac               !! Fraction of land in each grid box.
 
    REAL(r_std), DIMENSION(:), INTENT(in)                    :: drysoil_frac           !! Fraction of  dry soil (unitless)
    REAL(r_std), DIMENSION(:), INTENT(in)                    :: coszang                !! The cosine of the solar zenith angle (unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(in)                  :: veget_max              !! PFT coverage fraction of a PFT (= ind*cn_ind) (m^2 m^{-2})
    REAL(r_std),DIMENSION (:,:), INTENT(in)                  :: veget                  !! Fraction of vegetation types
    REAL(r_std), DIMENSION(:), INTENT(in)                    :: snow                   !! Snow mass in vegetation (kg m^{-2})           
    REAL(r_std), DIMENSION(:), INTENT(in)                    :: snow_age               !! Snow age (days)        
    REAL(r_std), DIMENSION(:,:), INTENT(in)                  :: snow_nobio_age         !! Snow age on continental ice, lakes, etc. (days)    
    REAL(r_std), DIMENSION(:,:), INTENT(in)                  :: frac_nobio             !! Fraction of non-vegetative surfaces
    REAL(r_std), DIMENSION(:), INTENT(in)                    :: z0m                    !! Roughness height for momentum of vegetated part (m)
    TYPE(laieff_type),DIMENSION(:,:,:), INTENT(in)           :: laieff_fit             !! Fitted parameters for the effective LAI
    REAL(r_std), DIMENSION(:,:), INTENT(in)                  :: frac_snow_nobio        !! Fraction of snow on continental ice, lakes, etc. (unitless ratio)
    REAL(r_std), DIMENSION(:), INTENT(in)                    :: frac_snow_veg          !! The fraction of ground covered with snow, between 0 and 1
    
   
   !! 0.2 Output variables
   
   REAL(r_std), DIMENSION(kjpindex,n_spectralbands), INTENT(out) :: albedo              !! Albedo (two stream radiation transfer model)
                                                                                        !! for visible and near-infrared range 
                                                                                        !! (unitless)
   REAL(r_std), DIMENSION(kjpindex,nvm,nlevels_tot), INTENT(out)  :: Isotrop_Abs_Tot_p  !! Absorbed radiation per layer for photosynthesis
   REAL(r_std), DIMENSION(kjpindex,nvm,nlevels_tot), INTENT(out)  :: Isotrop_Tran_Tot_p !! Transmitted radiation per layer for photosynthesis
   !+++CHECK+++
   ! Variable dimensions xxx_Tot_mean are for a single pixel. Causing 1+1 problems
   ! when running a larger domain. Needs to be corrected when implementing
   ! a global use of the multi-layer energy budget.
   REAL(r_std),DIMENSION(nlevels_tot), INTENT(out)                :: Light_Abs_Tot_mean !! total absorption for a given level * changed for enerbil integration
   REAL(r_std),DIMENSION(nlevels_tot), INTENT(out)                :: Light_Alb_Tot_mean !! total albedo for a given level * changed for enerbil integration
   !+++++++++++
   REAL(r_std), DIMENSION(kjpindex,nlevels_tot,nvm), INTENT(out)  :: laieff_isotrop     !! Leaf Area Index Effective converts
                                                                                        !! 3D lai into 1D lai for two stream
                                                                                        !! radiation transfer model...this is for 
                                                                                        !! isotropic light and only calculated once per day
                                                                                        !! @tex $(m^{2} m^{-2})$ @endtex      



    !! 0.4 Local variables   
    INTEGER                                                  :: ier, ji, ilevel
    LOGICAL                                                  :: init_needed             !! Variable indicating true if one or several restart 
                                                                                        !! variables were missing
    INTEGER(i_std), DIMENSION(kjpindex)                      :: index_dummy             !! Dummy array

!_ ================================================================================================================================
  
    ! Get local printlev variable
    printlev_loc=get_printlev('albedo')

    ! Set TEST_ALB_RT_VEGMAX=y in run.def to test previous version using veget (called alb1)
    test_alb_rt_vegmax=.TRUE.
    CALL getin_p("TEST_ALB_RT_VEGMAX", test_alb_rt_vegmax)
    
    ! Set TEST_ALB_RT_BGSNOWBARE=y in run.def to test a new method using snowa_veg(ibare_sechiba) at a precice moment
    test_alb_rt_bgsnowbare=.FALSE.
    CALL getin_p("TEST_ALB_RT_BGSNOWBARE", test_alb_rt_bgsnowbare)


    IF (alb_bg_modis) THEN
       ! Allocate background soil albedo
       ALLOCATE (bckgrnd_alb(kjpindex,2),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'albedo_surface_initialize','Pb in allocation for bckgrnd_alb','','')
 
       ! Read background albedo from restart file
       CALL ioconf_setatt_p('UNITS', '-')
       CALL ioconf_setatt_p('LONG_NAME','Background soil albedo for visible and near-infrared range')
       CALL restget_p (rest_id, 'bckgrnd_alb', nbp_glo, 2, 1, kjit, .TRUE., bckgrnd_alb, "gather", nbp_glo, index_g)
       
       ! Initialize by interpolating from file if the variable was not in restart file or if option alb_bg_modis_reinit=true
       IF ( alb_bg_modis_reinit .OR. ALL(bckgrnd_alb(:,:) == val_exp) ) THEN
          CALL read_background_albedo(kjpindex, lalo, neighbours, resolution, contfrac)
       END IF
       CALL xios_orchidee_send_field("bckgrnd_alb",bckgrnd_alb)

    ELSE
       ! Dry soil albedo
       ALLOCATE (soilalb_dry(kjpindex,2),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'albedo_surface_initialize','Pb in allocation for soilalb_dry','','')
       
       ! Wet soil albedo
       ALLOCATE (soilalb_wet(kjpindex,2),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'albedo_surface_initialize','Pb in allocation for soilalb_wet','','')
       
       ! Mean soil albedo
       ALLOCATE (soilalb_moy(kjpindex,2),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'albedo_surface_initialize','Pb in allocation for soilalb_moy','','')


       !! Read variables from restart file or initialize them
       ! dry soil albedo
       CALL ioconf_setatt_p('UNITS', '-')
       CALL ioconf_setatt_p('LONG_NAME','Dry bare soil albedo')
       CALL restget_p (rest_id,'soilalbedo_dry' , nbp_glo, 2, 1, kjit, .TRUE., soilalb_dry, "gather", nbp_glo, index_g)
       
       ! wet soil albedo
       CALL ioconf_setatt_p('UNITS', '-')
       CALL ioconf_setatt_p('LONG_NAME','Wet bare soil albedo')
       CALL restget_p (rest_id, 'soilalbedo_wet', nbp_glo, 2, 1, kjit, .TRUE., soilalb_wet, "gather", nbp_glo, index_g)
       
       ! mean soil aledo
       CALL ioconf_setatt_p('UNITS', '-')
       CALL ioconf_setatt_p('LONG_NAME','Mean bare soil albedo')
       CALL restget_p (rest_id, 'soilalbedo_moy', nbp_glo, 2, 1, kjit, .TRUE., soilalb_moy, "gather", nbp_glo, index_g)
       
       ! Initialize the variables if not found in restart file
       IF ( ALL(soilalb_wet(:,:) == val_exp) .OR. &
            ALL(soilalb_dry(:,:) == val_exp) .OR. &
            ALL(soilalb_moy(:,:) == val_exp)) THEN
          ! One or more of the variables were not in the restart file. 
          ! Call routine albedo_surface_soilalb to calculate them.
          CALL albedo_surface_soilalb(kjpindex, lalo, neighbours, resolution, contfrac)
          IF (printlev>=1) WRITE(numout,*) '---> val_exp ', val_exp
          IF (printlev>=1) WRITE(numout,*) '---> ALBEDO_wet VIS:', MINVAL(soilalb_wet(:,ivis)), MAXVAL(soilalb_wet(:,ivis))
          IF (printlev>=1) WRITE(numout,*) '---> ALBEDO_wet NIR:', MINVAL(soilalb_wet(:,inir)), MAXVAL(soilalb_wet(:,inir))
          IF (printlev>=1) WRITE(numout,*) '---> ALBEDO_dry VIS:', MINVAL(soilalb_dry(:,ivis)), MAXVAL(soilalb_dry(:,ivis))
          IF (printlev>=1) WRITE(numout,*) '---> ALBEDO_dry NIR:', MINVAL(soilalb_dry(:,inir)), MAXVAL(soilalb_dry(:,inir))
          IF (printlev>=1) WRITE(numout,*) '---> ALBEDO_moy VIS:', MINVAL(soilalb_moy(:,ivis)), MAXVAL(soilalb_moy(:,ivis))
          IF (printlev>=1) WRITE(numout,*) '---> ALBEDO_moy NIR:', MINVAL(soilalb_moy(:,inir)), MAXVAL(soilalb_moy(:,inir))
       ENDIF

    END IF

    !! Read variables from restart file
    !  If a variable is not found in the restart file, the flag init_needed is set to true.
    init_needed=.FALSE.
    ! albedo : DIM(kjpindex, n_spectralbands)
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','albedo')
    CALL restget_p (rest_id,'albedo' , nbp_glo, n_spectralbands, 1, kjit, .TRUE., &
         albedo, "gather", nbp_glo, index_g)
    IF (ALL(albedo(:,:) == val_exp)) init_needed=.TRUE.

    ! Isotrop_Abs_Tot_p : DIM(kjpindex,nvm,nlevels_tot)
    CALL ioconf_setatt_p('UNITS', '')
    CALL ioconf_setatt_p('LONG_NAME','Absorbed radiation per layer for photosynthesis')
    CALL restget_p (rest_id, 'Isotrop_Abs_Tot_p', nbp_glo, nvm, nlevels_tot, kjit, .TRUE., &
         Isotrop_Abs_Tot_p, "gather", nbp_glo, index_g)
    IF ( ALL(Isotrop_Abs_Tot_p(:,:,:) ==  val_exp ) ) init_needed=.TRUE.

    ! Isotrop_Tran_Tot_p : DIM(kjpindex,nvm,nlevels_tot)
    CALL ioconf_setatt_p('UNITS', '')
    CALL ioconf_setatt_p('LONG_NAME','Transmitted radiation per layer for photosynthesis')
    CALL restget_p (rest_id, 'Isotrop_Tran_Tot_p', nbp_glo, nvm, nlevels_tot, kjit, .TRUE., &
         Isotrop_Tran_Tot_p, "gather", nbp_glo, index_g)
    IF ( ALL(Isotrop_Tran_Tot_p(:,:,:) ==  val_exp ) ) init_needed=.TRUE.

    !+++CHECK+++
    ! Variable dimensions xxx_Tot_mean are for a single pixel. Causing 1+1 problems
    ! when running a larger domain. Needs to be corrected when implementing
    ! a global use of the multi-layer energy budget.
    ! Light_Abs_Tot_mean, Light_Alb_Tot_mean : DIM(nlevels_tot) : Needed in mleb_main

    IF (ok_mleb) THEN
       IF (is_root_prc) THEN

          WRITE(numout,*) "LATm", Light_Abs_Tot_mean(:)
          CALL ioconf_setatt_p('UNITS', '')
          CALL ioconf_setatt_p('LONG_NAME','')
          CALL restget (rest_id, 'Light_Abs_Tot_mean', nlevels_tot, 1, 1, kjit, .TRUE., Light_Abs_Tot_mean)

          CALL ioconf_setatt_p('UNITS', '')
          CALL ioconf_setatt_p('LONG_NAME','')
          CALL restget (rest_id, 'Light_Alb_Tot_mean', nlevels_tot, 1, 1, kjit, .TRUE., Light_Alb_Tot_mean)
       END IF
       CALL bcast(Light_Abs_Tot_mean)
       CALL bcast(Light_Alb_Tot_mean)
       IF ( ALL(Light_Abs_Tot_mean(:) ==  val_exp ) ) init_needed=.TRUE.
       IF ( ALL(Light_Alb_Tot_mean(:) ==  val_exp ) ) init_needed=.TRUE.
    ENDIF
    ! laieff_isotrop : DIM(kjpindex, nlevels_tot,nvm) : Needed in mleb_main called before condveg_main
    !+++++++++++
    CALL ioconf_setatt_p('UNITS', '')
    CALL ioconf_setatt_p('LONG_NAME','')
    CALL restget_p (rest_id, 'laieff_isotrop', nbp_glo, nlevels_tot, nvm, kjit, .TRUE., &
         laieff_isotrop, "gather", nbp_glo, index_g)
    IF ( ALL(laieff_isotrop(:,:,:) ==  val_exp ) ) init_needed=.TRUE.


    IF (init_needed) THEN
       ! One or more of the variables were not in the restart file. 
       ! Call routine albedo_surface_main to calculate them.
       ! In the initialization phase, the arguments for hist_id, hist_id2 and index will never be used. 
       ! Therefor send dummy variables.
       index_dummy(:)=0
       !+++CHECK+++
       ! Variable dimensions xxx_Tot_mean are for a single pixel. Causing 1+1 problems
       ! when running a larger domain. Needs to be corrected when implementing
       ! a global use of the multi-layer energy budget.
       CALL albedo_surface_main( &
            kjit,                kjpindex,            lalo,              0,  0,          &
            index_dummy,         drysoil_frac,        veget_max,         coszang,        &
            frac_nobio,          frac_snow_veg,       frac_snow_nobio,                   &
            snow,                snow_age,            snow_nobio_age,                    &
            z0m,                 laieff_fit,                                             &
            albedo,              Isotrop_Abs_Tot_p,   Isotrop_Tran_Tot_p,                &
            Light_Abs_Tot_mean, Light_Alb_Tot_mean,   &
            laieff_isotrop, veget)
       !+++++++++++
    END IF

    ! Initalization finished
    firstcall_albedo_surface=.FALSE.

  END SUBROUTINE albedo_surface_initialize

!! ==============================================================================================================================
!! SUBROUTINE   : albedo_surface_main
!!
!>\BRIEF        This subroutine calculates the albedo for two stream radiation transfer model.  This routine includes
!!              the effect of snow on the background albedo, through one of two methods.
!!
!! DESCRIPTION  : The albedo for two stream radiation transfer model  is calculated for both the visible and near-infrared 
!! domain. First the mean albedo of the bare soil is calculated. Two options exist: 
!! either the soil albedo depends on soil wetness (drysoil_frac variable), or the soil albedo 
!! is set to a mean soil albedo value.
!!
!!    NOTE: the main output variable, albedo, is an unweighted average of the direct and diffuse albedos
!!          for each grid point...this is done in this way right now to be consistent with the new scheme,
!!          but it doesn't have to be combined like that once we have an energy budget which can
!!          use both types of light directly
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::albedo
!!
!! REFERENCE(S) : B. Pinty, T. Lavergne, R.E. Dickinson, J-L. Widlowski, N. Gobron and M. M. Verstraete (2006).
!! Simplifying the Interaction of Land Surfaces with Radiation for Relating Remote Sensing Products to Climate Models.
!! Journal of Geophysical Research. Vol 111, D02116.
!!
!! FLOWCHART    : None
!! \n
!!
!!
!!
!_ ================================================================================================================================
  !+++CHECK+++
  ! Variable dimensions xxx_Tot_mean are for a single pixel. Causing 1+1 problems
  ! when running a larger domain. Needs to be corrected when implementing
  ! a global use of the multi-layer energy budget.
  SUBROUTINE albedo_surface_main(                                                   &
       kjit,                kjpindex,            lalo,              hist_id, hist2_id, &
       index,               drysoil_frac,        veget_max,         coszang,        &
       frac_nobio,          frac_snow_veg,       frac_snow_nobio,                   &
       snow,                snow_age,            snow_nobio_age,                    &
       z0m,                 laieff_fit,                                             &
       albedo,              Light_Abs_Tot,       Light_Tran_Tot,                    &
       Light_Abs_Tot_mean,  Light_Alb_Tot_mean,  &
       laieff_isotrop, veget)
    !+++++++++++

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables
   INTEGER(i_std), INTENT(in)                                :: kjit                   !! Time step number 
   INTEGER(i_std), INTENT(in)                                :: kjpindex               !! Domain size - Number of land pixels  (unitless)       
   REAL(r_std), DIMENSION(kjpindex,2), INTENT (in)           :: lalo                   !! Geographical coordinates
   INTEGER(i_std), INTENT (in)                               :: hist_id                !! _History_ file identifier
   INTEGER(i_std), INTENT (in)                               :: hist2_id               !! _History_ file 2 identifier
   INTEGER(i_std), DIMENSION(kjpindex), INTENT (in)          :: index                  !! Indeces of the points on the map
   REAL(r_std), DIMENSION(:), INTENT(in)                     :: drysoil_frac           !! Fraction of  dry soil (unitless)
   REAL(r_std), DIMENSION(:,:), INTENT(in)                   :: veget_max              !! PFT coverage fraction of a PFT (= ind*cn_ind) 
                                                                                       !! (m^2 m^{-2})
   REAL(r_std),DIMENSION (:,:), INTENT(in)                   :: veget                  !! Fraction of vegetation types
   REAL(r_std), DIMENSION(:), INTENT(in)                     :: coszang                !! The cosine of the solar zenith angle (unitless)
   REAL(r_std),DIMENSION (:,:), INTENT(in)                   :: frac_nobio             !! Fraction of non-vegetative surfaces, i.e. 
                                                                                       !! continental ice, lakes, etc. (unitless)     
   REAL(r_std),DIMENSION(:), INTENT(in)                      :: frac_snow_veg          !! The fraction of ground covered with snow, between 0 and 1
   REAL(r_std), DIMENSION(:,:), INTENT(in)                   :: frac_snow_nobio        !! Fraction of snow on continental ice, lakes, etc. 
                                                                                       !! (unitless ratio)
   REAL(r_std),DIMENSION (:), INTENT(in)                     :: snow                   !! Snow mass in vegetation (kg m^{-2})           
   REAL(r_std),DIMENSION (:), INTENT(in)                     :: snow_age               !! Snow age (days)        
   REAL(r_std),DIMENSION (:,:), INTENT(in)                   :: snow_nobio_age         !! Snow age on continental ice, lakes, etc. (days)    
   REAL(r_std), DIMENSION(:), INTENT(in)                     :: z0m                    !! Roughness height for momentum of vegetated part (m)
   TYPE(laieff_type),DIMENSION (:,:,:),INTENT(in)            :: laieff_fit             !! Fitted parameters for the effective LAI
   
   !! 0.2 Output variables
   
   REAL(r_std), DIMENSION (kjpindex,n_spectralbands), &
                                 INTENT (out)                :: albedo                 !! Albedo (two stream radiation transfer model)
                                                                                       !! for visible and near-infrared range 
                                                                                       !! (unitless)
   REAL(r_std),DIMENSION (kjpindex,nvm,nlevels_tot), INTENT(out) :: Light_Abs_Tot      !! Absorbed radiation per layer, averaged between light
                                                                                       !! from a collimated (direct) source and that from an
                                                                                       !! isotropic (diffuse) source.  Expressed as the fraction
                                                                                       !! of overall light hitting the canopy.
   REAL(r_std),DIMENSION (kjpindex,nvm,nlevels_tot), INTENT(out) :: Light_Tran_Tot     !! Transmitted radiation per layer, averaged between light
                                                                                       !! from a collimated (direct) source and that from an
                                                                                       !! isotropic (diffuse) source.  Expressed as the fraction
                                                                                       !! of overall light hitting the canopy.
   !+++CHECK+++
   ! Variable dimensions xxx_Tot_mean are for a single pixel. Causing 1+1 problems
   ! when running a larger domain. Needs to be corrected when implementing
   ! a global use of the multi-layer energy budget.
   REAL(r_std),DIMENSION(nlevels_tot), INTENT (out)          :: Light_Abs_Tot_mean     !! total absorption for a given level * changed for enerbil integration
   REAL(r_std),DIMENSION(nlevels_tot), INTENT (out)          :: Light_Alb_Tot_mean     !! total albedo for a given level * changed for enerbil integration
   !+++++++++++
   REAL(r_std), DIMENSION(kjpindex,nlevels_tot,nvm), INTENT(out) :: laieff_isotrop     !! Leaf Area Index Effective converts
                                                                                       !! 3D lai into 1D lai for two stream
                                                                                       !! radiation transfer model...this is for 
                                                                                       !! isotropic light and only calculated once per day
                                                                                       !! @tex $(m^{2} m^{-2})$ @endtex      



   !! 0.3 Modified variables
 
 
   !! 0.4 Local variables

   REAL(r_std),DIMENSION (nlevels_tot)   :: laieff_isotrop_pft, laieff_collim_pft      !! 
   REAL(r_std)                           :: cosine_sun_angle                           !! the cosine of the solar zenith angle
   INTEGER(i_std)                        :: ks                                         !! Index for visible and near-infraread range
   INTEGER(i_std)                        :: ivm                                        !! Index for vegetative PFTs
   INTEGER(i_std)                        :: ipts                                       !! Index for spatial domain
   INTEGER(i_std)                        :: ilevel                                     !! Index for canopy levels
   REAL(r_std)                           :: leaf_reflectance
   REAL(r_std)                           :: leaf_transmittance
   REAL(r_std)                           :: br_base_temp,br_base_temp_collim,br_base_temp_isotrop
   REAL(r_std)                           :: leaf_psd_temp
   REAL(r_std)                           :: leaf_single_scattering_albedo

   REAL(r_std), DIMENSION(kjpindex,n_spectralbands)      :: alb_bare                   !! Mean bare soil albedo for visible and near-infrared
   REAL(r_std), DIMENSION(kjpindex,n_spectralbands)      :: albedo_snow                !! Snow albedo (unitless ratio)     
   REAL(r_std), DIMENSION(kjpindex)                      :: albedo_glob                !! Mean albedo (unitless ratio)
   REAL(r_std), DIMENSION(kjpindex,nvm, n_spectralbands) :: albedo_pft                 !! Albedo (two stream radiation transfer model)
                                                                                       !! for visible and near-infrared range for each PFT (unitless)
   REAL(r_std), DIMENSION(kjpindex,nlevels_tot,nvm) :: laieff_collim                   !! Leaf Area Index Effective for direct light
   REAL(r_std),DIMENSION(nlevels_tot)    :: Collim_Tran_Uncoll                         !! collimated transmission of uncollided light for a given level 
   REAL(r_std),DIMENSION(nlevels_tot)    :: Collim_Tran_Coll                           !! collimated transmission of collided for a given level
   REAL(r_std),DIMENSION(nlevels_tot)    :: Collim_Tran_Tot                            !! collimated total transmission
   REAL(r_std),DIMENSION(nlevels_tot)    :: Collim_Abs_Tot                             !! collimated total absorption for a given level        * changed for enerbil integration
   REAL(r_std),DIMENSION(nlevels_tot)    :: Collim_Alb_Tot                             !! Collimated (direct) total albedo for a given level   * changed for enerbil integration

   REAL(r_std),DIMENSION(nlevels_tot)    :: Isotrop_Alb_Tot                            !! isotropic (diffuse) total albedo for a given level
   REAL(r_std),DIMENSION(nlevels_tot)    :: Isotrop_Tran_Uncoll                        !! isotropic transmission of uncollided light for a given level
   REAL(r_std),DIMENSION(nlevels_tot)    :: Isotrop_Tran_Coll                          !! isotropic transmission of collided for a given level
   REAL(r_std),DIMENSION(nlevels_tot)    :: Isotrop_Tran_Tot                           !! isotropic total transmission
   REAL(r_std),DIMENSION(nlevels_tot)    :: Isotrop_Abs_Tot                            !! isotropic total absporbtion for a given level
   INTEGER :: jv, inno
   REAL(r_std) :: alb_nobio
   REAL(r_std):: laieff_collim_1, laieff_isotrop_1, Collim_Alb_Tot_1,Collim_Tran_Tot_1,Collim_Abs_Tot_1,&
                 Isotrop_Alb_Tot_1,Isotrop_Tran_Tot_1,Isotrop_Abs_Tot_1,&
                 Collim_Tran_Uncollided_1, Isotrop_Tran_Uncollided_1                   !! Values for the one level solution

   REAL(r_std):: converged_albedo
   REAL(r_std):: isotrop_angle
   LOGICAL :: lconverged
   LOGICAL :: lprint                                                                   !! a flag for printing some debug statements
   REAL(r_std), DIMENSION(kjpindex,nvm,n_spectralbands) :: snowa_veg                   !! snow albedo due to vegetated surfaces, between 0 and 1
   REAL(r_std), DIMENSION(kjpindex,nnobio,n_spectralbands) :: snowa_nobio              !! Albedo of snow covered area on continental ice, 
                                                                                       !! lakes, etc. (unitless ratio)  

   REAL(r_std), DIMENSION(kjpindex,n_spectralbands)      :: albedo_diag                !! Array used to filter the values that will be send to xios
   REAL(r_std), DIMENSION(kjpindex,n_spectralbands)      :: albedo_noon                !! Array keeping albedo values only at local noon
   REAL(r_std), DIMENSION(kjpindex,n_spectralbands)      :: alb_bare_diag              !! Array used to filter the values that will be send to xios
   REAL(r_std), DIMENSION(kjpindex,n_spectralbands)      :: albedo_snow_diag           !! Array used to filter the values that will be send to xios
   REAL(r_std)                                           :: leaf_w, leaf_d             !! parameters related to the single scattering albedo
   REAL(r_std),PARAMETER                                 :: cosine_iso_angle = 0.5_r_std / 0.705_r_std
!_ ================================================================================================================================

   IF (printlev_loc>=3) WRITE(numout,*) 'Entering albedo_surface_main'

   IF (printlev_loc>=3) WRITE(numout,*) 'nlevels_tot',nlevels_tot
   IF (printlev_loc>=4) THEN
      WRITE(numout,*) 'laieff_isotrop start albedo for test_pft',test_pft,'is', &
           laieff_isotrop(test_grid,:,test_pft)

      WRITE(numout,*)'albedo.f90, coszang',coszang(:)
   ENDIF

   
   !! 1. We will now calculate the background reflectance to be used in the model.
   !  The parameters read in from the input file do not include the effect 
   !  of snow, so we need to figure out the snow's contribution for each 
   !  grid space
   !  Initialize some output variables
   albedo_pft(:,:,:)=zero
   Light_Abs_Tot(:,:,:)=zero
   Light_Tran_Tot(:,:,:)=un
   
   !+++CHECK+++
   ! Variable dimensions xxx_Tot_mean are for a single pixel. Causing 1+1 problems
   ! when running a larger domain. Needs to be corrected when implementing
   ! a global use of the multi-layer energy budget.
   Light_Abs_Tot_mean(:) = zero
   Light_Alb_Tot_mean(:) = zero

   !+++++++++++


   ! We need to calculate the LAI effective for this solar angle. The other LAI
   ! effective that we calculated was computed at an angle of 60.0 degrees for
   ! the isotropic contribution. Notice that, in many situations, the results 
   ! will be exactly the same.  One can show that Pgap is proportional to 
   ! 1/cos(theta) if the level bottom is below all canopy elements, which
   ! means that the cos(theta) in the calculation of the LAIeff will be
   ! cancelled out
   isotrop_angle=COS(60.0/180.0*pi)
   DO ipts=1,kjpindex
      DO ivm=1,nvm
         DO ilevel=1,nlevels_tot
            laieff_isotrop(ipts,ilevel,ivm) = &
                 calculate_laieff_fit(isotrop_angle,laieff_fit(ipts,ivm,ilevel))
            laieff_collim(ipts,ilevel,ivm) = &
                 calculate_laieff_fit(coszang(ipts),laieff_fit(ipts,ivm,ilevel))
            ! Because we are fitting these values, it's possible that some 
            ! will drop below zero.  This causes a serious problem here,
            ! so let's make sure this doesn't happen.  Unfortunately, by
            ! doing this, we lose a test of things going wrong (an unfitted
            ! negative LAIeff can be the sign of another problem), but we
            ! don't really have a choice.
            IF(laieff_isotrop(ipts,ilevel,ivm) .LT. min_stomate) &
                 laieff_isotrop(ipts,ilevel,ivm)=zero
            IF(laieff_collim(ipts,ilevel,ivm) .LT. min_stomate) &
                 laieff_collim(ipts,ilevel,ivm)=zero
         ENDDO
      ENDDO
   ENDDO

   ! Debug
   IF(printlev_loc>=4)THEN
      WRITE(numout,*) 'ACOS(coszang(ipts))/pi*180: ',&
           ACOS(coszang(test_grid))/pi*180.0_r_std
      WRITE(numout,*) 'coszang(ipst)', coszang
      WRITE(numout,*) 'ivel, laieff_collim, laieff_isotrop in albedo', &
           test_grid, test_pft
      DO ilevel=nlevels_tot,1,-1
         WRITE(numout,*) ilevel,laieff_collim(test_grid,ilevel,test_pft),&
              laieff_isotrop(test_grid,ilevel,test_pft)
      END DO
      DO ilevel=nlevels_tot,1,-1
         WRITE(numout,*) 'Fitting parameters: ',&
              laieff_fit(test_grid,test_pft,ilevel)%a,&
              laieff_fit(test_grid,test_pft,ilevel)%b,&
              laieff_fit(test_grid,test_pft,ilevel)%c,&
              laieff_fit(test_grid,test_pft,ilevel)%d,&
              laieff_fit(test_grid,test_pft,ilevel)%e
      ENDDO
   ENDIF
   !-

   ! Calculate the snow albedo
   CALL calculate_snow_albedo(kjpindex,  coszang,  snow, &
        snow_age, snow_nobio_age,  frac_nobio,  albedo_snow, &
        snowa_veg,  frac_snow_veg,  snowa_nobio,  frac_snow_nobio, &
        veget_max, z0m, veget)

   ! Check for possible problems in the calcualtion of the fraction
   ! of the grid cell covered by snow. Fractions should be positive
   ! numbers
   DO ipts=1,kjpindex
      IF(frac_snow_veg(ipts) .LT. 0.0)THEN
         WRITE(numout,*) 'SNOWFRAC ipts ',ipts,frac_snow_veg(ipts)
         CALL ipslerr_p (3,'albedo', 'snow frac error','','')
      END IF
     
      DO inno=1,nnobio
         IF(frac_snow_nobio(ipts,inno) .LT. 0.0)THEN
            WRITE(numout,*) 'SNOWFRAC ipts inno',ipts,inno,frac_snow_nobio(ipts,inno)
            CALL ipslerr_p (3,'albedo', 'nobio snow frac error','','')
         ENDIF
      ENDDO
   ENDDO

   ! Initialize the albedo in every square for both spectra by calculating the 
   ! bare soil albedo
   DO ipts = 1, kjpindex
      DO ks = 1, n_spectralbands

         ! If 'alb_bare_model' is set to TRUE,  the soil albedo calculation 
         ! depends on soil moisture
         ! If 'alb_bare_model' is set to FALSE, the mean soil albedo is used 
         ! without the dependance on soil moisture
         ! see subroutine 'conveg_soilalb'
         ! Note that these two methods are considered old.  The albedo 
         ! background map has been done more recently, as is now
         ! used as the default in simulations with TAG 2.2.  The 
         ! lines below ensure that the albedos over the Sahara are identical
         ! for Tags 2.2, 3.0, and 4.0 when alb_bg_modis=.TRUE.
         IF ( alb_bg_modis ) THEN
            alb_bare(ipts,ks) = bckgrnd_alb(ipts,ks)
         ELSEIF ( alb_bare_model ) THEN
            alb_bare(ipts,ks) = soilalb_wet(ipts,ks) + drysoil_frac(ipts) * &
                 (soilalb_dry(ipts,ks) -  soilalb_wet(ipts,ks))
         ELSE
            alb_bare(ipts,ks) = soilalb_moy(ipts,ks)
         ENDIF

         ! We initialise the albedo with the snow / non-snow part of Bare soil.
         albedo(ipts,ks) = veget_max(ipts,ibare_sechiba) * ( alb_bare(ipts,ks) * &
              (un-frac_snow_veg(ipts)) + frac_snow_veg(ipts) * snowa_veg(ipts,ibare_sechiba,ks) )         
      
      ENDDO ! ks = 1, n_spectralbands
      
   ENDDO ! ipts = 1, kjpindex

   !!    PFTs (except bare soil) and spectral bands
   DO ipts = 1, kjpindex ! loop over the grid squares
      !+++CHECK+++
      ! observations from the MERGE, sinang is calculated as cos of the zenith
      ! angle in solar.f90 (i.e. it is the cosine of the zenith angle). 
      ! Then maybe also the cosine_sun_angle is a redundant variable.
      cosine_sun_angle = coszang(ipts)
      !+++++++++++

      ! Debug
      IF(printlev_loc >= 4) THEN
         WRITE(numout,*)'cosine_sun_angle',cosine_sun_angle
      ENDIF
      !-

      ! For the coupled model, the angles can be negative.
      ! Offline, they are always between zero and 90 degrees.
      IF(cosine_sun_angle .LE. min_sechiba)THEN

         ! It's night, so no sunlight...don't need to calculate albedo.  
         ! Set it equal to zero, because the snow albedo has already been 
         ! calculated and it looks funny on the coupled run maps otherwise.
         albedo(ipts,:)=zero
         CYCLE

      ENDIF
      ! Are there any non-bio PFTs here that we need to take into account?
      DO ks = 1, n_spectralbands ! Loop over # of spectra
         DO jv = 1, nnobio 
            ! now update the albedo
            IF ( jv .EQ. iice ) THEN
               alb_nobio = alb_ice(ks)
            ELSE
               WRITE(numout,*) 'jv=',jv
               CALL ipslerr_p (3,'Albedo.f90', &
                    'DO NOT KNOW ALBEDO OF THIS SURFACE TYPE','','')
            ENDIF

            ! this takes into account both snow covered and non-snow 
            ! covered non-bio regios in all grid squares
            albedo(ipts,ks) = albedo(ipts,ks) + &
                 ( frac_nobio(ipts,jv) ) * &
                 ( (un-frac_snow_nobio(ipts,jv)) * alb_nobio + &
                 ( frac_snow_nobio(ipts,jv)  ) * snowa_nobio(ipts,jv,ks)   )
            
         ENDDO ! jv = 1, nnobio 
      ENDDO ! ks = 1, n_spectralbands 
      
      ! Calculate the albedo of the vegetated fractions of the pixel
      DO ivm = 2, nvm  ! Loop over # of PFTs 
         ! for this grid square and this vegetation type, we have a set LAI 
         ! effective
         laieff_collim_pft(1:nlevels_tot) = &
              laieff_collim(ipts,1:nlevels_tot,ivm)
         laieff_isotrop_pft(1:nlevels_tot) = &
              laieff_isotrop(ipts,1:nlevels_tot,ivm)

         IF(veget_max(ipts,ivm) == zero)THEN
            ! this vegetation type is not present, so no reason to do the 
            ! calculation
            CYCLE
         ENDIF
         DO ks = 1, n_spectralbands ! Loop over # of spectra
            ! set single scattering albedo and preferred scattering direction
            leaf_single_scattering_albedo = leaf_ssa(ivm,ks)
            leaf_psd_temp = leaf_psd(ivm,ks)                 
            ! calculate the rleaf and tleaf from wl=rl+tl and dl=rl/tl
            leaf_transmittance = leaf_single_scattering_albedo/ & 
                 (leaf_psd_temp+1)
            leaf_reflectance = leaf_psd_temp * &
                 leaf_transmittance
            ! We need to take into account the effect of snow on the 
            ! background reflectance here. From the snow above I can 
            ! calculate the true background reflectance for this PFT.
            ! The flag alb_bg_modis is a bit confusing because in 
            ! both cases the data source is the MODIS albedo product
!PP
!PP The IF below is not fully coherent with the above calculation of the Bare soil albedo: if alb_bg_modis = T then it is coherent 
!PP but else we have above the use of another flag "alb_bare_model" not used here... To be discussed : the difference btw bg albedo and Bare soil albedo !
!PP 
            IF (alb_bg_modis) THEN

               ! If the flag is set to TRUE the model will read a 
               ! spatially explicit map of the background albedo for
               ! each PIXEL. The background albedo of the pixel is 
               ! then used to calculate the background albedo of the
               ! PFT accounting for snow fraction and snow albedo. Note
               ! that both snow fraction and snow albedo vary over time.

               ! Note that for Tree we assume that all snow goes on the ground (Tree not covered by snow)
               ! While for short vegetation the snow stays on the top of the grass (no impact on the ground)
               IF (is_tree(ivm) .AND. test_alb_rt_bgsnowbare) THEN
                  ! This is a new test using snowa_veg(ibare_sechiba)
                  br_base_temp = (un-frac_snow_veg(ipts)) * bckgrnd_alb(ipts,ks) + &
                       frac_snow_veg(ipts) * snowa_veg(ipts,ibare_sechiba,ks)
               ELSE IF (is_tree(ivm)) THEN
                  ! This is the default case using snowa_veg(ivm)
                  br_base_temp = (un-frac_snow_veg(ipts)) * bckgrnd_alb(ipts,ks) + &
                       frac_snow_veg(ipts) * snowa_veg(ipts,ivm,ks)
               ELSE
                  br_base_temp = bckgrnd_alb(ipts,ks) 
               ENDIF

               ! Note: this case (alb_bg_modis = F) is not fully coherent with the above calculation 
               ! of background albedo for bare soil fraction
            ELSE
               
               ! If the flag is set to FALSE the model uses the original
               ! parameterization at the PFT level. Each PFT has a fixed
               ! background albedo. The background albedo is recalculated
               ! taking into account the snowfraction and the snow albedo.
               ! Note that snow fraction and snow albedo both vary over
               ! time.

               ! Note that for Tree we assume that all snow goes on the ground (Tree not covered by snow)
               ! While for short vegetation the snow stays on the top of the grass (no impact on the ground)
               IF (is_tree(ivm) .AND. test_alb_rt_bgsnowbare) THEN
                  ! This is a new test using snowa_veg(ibare_sechiba)
                  br_base_temp= (un-frac_snow_veg(ipts)) * bgd_reflectance(ivm,ks) + &
                       frac_snow_veg(ipts) * snowa_veg(ipts,ibare_sechiba,ks)
               ELSE IF (is_tree(ivm)) THEN
                  ! This is the default case using snowa_veg(ivm)
                  br_base_temp= (un-frac_snow_veg(ipts)) * bgd_reflectance(ivm,ks) + &
                       frac_snow_veg(ipts) * snowa_veg(ipts,ivm,ks)
               ELSE
                  br_base_temp = bgd_reflectance(ivm,ks) 
               ENDIF

            ENDIF

            ! At this step, we assume that there is no difference between the 
            ! direct and the diffuse background reflectance, which is true for 
            ! the background but not the canopy albedo
            br_base_temp_collim=br_base_temp
            br_base_temp_isotrop=br_base_temp
            
            ! Debug - Set flag for extra output
            IF(printlev_loc>=4 .AND. test_pft == ivm .AND. &
                 test_grid == ipts)THEN
               lprint=.TRUE.
            ELSE
               lprint=.FALSE.
            ENDIF
            !-
    
            ! Now solve the multilevel scheme
            ! The choice is between Iterative and Matrixial. 
            ! The default is the Iterative scheme, unless the user defines which_rt in the orchidee.def file.
            IF (which_rt .EQ. 'Iterative') THEN 
               CALL multilevel_albedo(cosine_sun_angle, &
                       leaf_single_scattering_albedo, &
                       leaf_psd_temp, br_base_temp_collim, br_base_temp_isotrop, &
                       laieff_collim_pft, laieff_isotrop_pft, lconverged, &
                       Collim_Alb_Tot, Collim_Tran_Coll, Collim_Abs_Tot, &
                       Isotrop_Alb_Tot, &
                       Isotrop_Tran_Coll, Isotrop_Abs_Tot, Collim_Tran_Uncoll, &
                       Isotrop_Tran_Uncoll, lprint)
               Collim_Tran_Tot(:) = Collim_Tran_Uncoll(:)+Collim_Tran_Coll(:)
               Isotrop_Tran_Tot(:) = Isotrop_Tran_Uncoll(:)+Isotrop_Tran_Coll(:)
            ELSE
               ! WHICH_RT = Matrixial 
               leaf_w = leaf_transmittance + leaf_reflectance
               leaf_d = leaf_reflectance   - leaf_transmittance
               CALL multilevel_matrix(nlevels_tot-1,cosine_sun_angle,br_base_temp_collim ,laieff_collim_pft ,leaf_w,leaf_d,&
                    Collim_Alb_Tot,Collim_Tran_Tot,Collim_Abs_Tot)
               CALL multilevel_matrix(nlevels_tot-1,cosine_iso_angle,br_base_temp_isotrop,laieff_isotrop_pft,leaf_w,leaf_d,&
                    Isotrop_Alb_Tot,Isotrop_Tran_Tot,Isotrop_Abs_Tot)
            ENDIF

            !+++CHECK+++
            ! Variable dimensions xxx_Tot_mean are for a single pixel. Causing 1+1 problems
            ! when running a larger domain. Needs to be corrected when implementing
            ! a global use of the multi-layer energy budget.
            ! Here we weight the light in a consistent way between diffuse and
            ! direct sources
            
            !! Debug ID12.
            !! Mleb purpose. Next step is to adapt them for global run.
            Light_Abs_Tot_mean(:) = Light_Abs_Tot_mean(:) + &
                 direct_light_weight*Collim_Abs_Tot(:)+ (un-direct_light_weight)*Isotrop_Abs_Tot(:)
            Light_Alb_Tot_mean(:) = Light_Alb_Tot_mean(:) + &
                 direct_light_weight*Collim_Alb_Tot(:)+ (un-direct_light_weight)*Isotrop_Alb_Tot(:)
            !! End ID12

            !+++++++++++

            ! This is a lot of information printed out for the paper
            ! on the multilevel albedo scheme (McGrath et al 2016, GMDD).  
            ! Probably not useful for anyone else, but I'm keeping it in 
            ! until the paper is published. This is why I'm protecting it 
            ! with an IF statement.  I do not use printlev_loc>=4 since it 
            ! really is too much information for normal usage.


            ! For our total albedo, we take a weighted average of the 
            ! diffuse and the direct light, which means we lose some 
            ! information. This might be changed in the future. 
            converged_albedo = direct_light_weight*Collim_Alb_Tot(nlevels_tot) + &
                 (un-direct_light_weight)*Isotrop_Alb_Tot(nlevels_tot)


            ! Addition of the different PFT contribution to the albedo
            ! i) Short vegetation (grass, crop): we assume that when there
            ! are no leaves (LAI=0), the snow aging is not impacted by
            ! the vegetation. The difference between veget and vegetmax
            ! (bare soil fraction within the PFT) should be considered in
            ! the same way as the bare soil fraction (PFT1).
            ! ii) High vegetation (trees): In this case we consider that
            ! the vegetation impacts the snow aging even when LAI=0
            ! because of the trunks, branches. So these PFTs should be
            ! treated differently without any consideration of the bare
            ! soil fraction within the PFT.
            
            IF (.NOT. test_alb_rt_vegmax) THEN
               ! test_alb_rt_vegmax=false: Corresponds to old test alb1
               IF (is_tree(ivm)) THEN
                  albedo(ipts,ks) = albedo(ipts,ks) + &
                       veget(ipts,ivm) * converged_albedo + &
                       ((veget_max(ipts,ivm) - veget(ipts,ivm)) * ((un-frac_snow_veg(ipts))*alb_bare(ipts,ks) + frac_snow_veg(ipts)*snowa_veg(ipts,ibare_sechiba,ks)))
                  albedo_pft(ipts,ivm,ks) = converged_albedo

               ELSE
                  albedo(ipts,ks) = albedo(ipts,ks) + &
                       veget(ipts,ivm) * ((un-frac_snow_veg(ipts))*converged_albedo + &
                       frac_snow_veg(ipts)*snowa_veg(ipts,ivm,ks)) + &
                       (veget_max(ipts,ivm) - veget(ipts,ivm)) * ((un-frac_snow_veg(ipts))*alb_bare(ipts,ks) + frac_snow_veg(ipts)*snowa_veg(ipts,ibare_sechiba,ks))
                  albedo_pft(ipts,ivm,ks) = (un-frac_snow_veg(ipts))*converged_albedo + &
                       frac_snow_veg(ipts)*snowa_veg(ipts,ivm,ks)
               ENDIF

            ELSE
               ! test_alb_rt_vegmax=true: This is as alb2: the version which was choosen
               IF (is_tree(ivm)) THEN
                  albedo(ipts,ks) = albedo(ipts,ks) + &
                       veget_max(ipts,ivm) * converged_albedo 
                  albedo_pft(ipts,ivm,ks) = converged_albedo

               ELSE
                  albedo(ipts,ks) = albedo(ipts,ks) + &
                       veget_max(ipts,ivm) * ((un-frac_snow_veg(ipts))*converged_albedo + &
                       frac_snow_veg(ipts)*snowa_veg(ipts,ivm,ks)) 
                  albedo_pft(ipts,ivm,ks) = (un-frac_snow_veg(ipts))*converged_albedo + &
                       frac_snow_veg(ipts)*snowa_veg(ipts,ivm,ks)
               ENDIF

            ENDIF

            ! Save the absorbed radiation for photosynthesis. We need only the
            ! visible range.  We weight the isotropic and collimated parts
            ! according to an external parameter, but 0.5 is a sensible value
            ! if we don't actually know the fraction of light coming from direct
            ! and diffuse sources.  This makes it consistent with the albedo
            ! calculated above.  This can be changed later to distinguish between 
            ! sunlit and shaded leaves.
            ! Notice that sunlight leaves will only get light from a Collimated
            ! source, whereas shaded leaves will get light from both direct
            ! sunlight (that has bounced off other leaves...Collim in this
            ! subroutine) and sunilght that has already reflected off of aerosols
            ! and clouds (Isotrop in this subroutine).  Therefore, it is not a
            ! simple matter of taking Isotrop_Abs_Tot or Collim_Abs_Tot, but 
            ! the Uncollided fraction of Collim_Abs_Tot (for the sunlit leaves)
            ! and the Collided fraction of Collum_Abs_Tot and the collided
            ! and uncollided fraction of Isotrop_Abs_Tot (for shaded leaves).
            !
            !PP Question: should the following re-ajustment for negative values be done 
            !   before the calculation of the Albedo variable hhat will be used for the Energie budget ?
            IF (ks == 1) THEN
               ! Test if Collim_Abs_Tot has negative values
               ! If so, set it to min_sechiba  
               DO  ilevel=1,nlevels_tot
                  IF (Isotrop_Abs_Tot(ilevel) .LT. zero) THEN
                     Isotrop_Abs_Tot(ilevel) =  min_sechiba
                  ENDIF
                  IF (Collim_Abs_Tot(ilevel) .LT. zero) THEN
                     Collim_Abs_Tot(ilevel) =  min_sechiba
                  ENDIF
               ENDDO
               !
               Light_Abs_Tot(ipts,ivm,:) = direct_light_weight*Collim_Abs_Tot(:) + &
                    (un-direct_light_weight)*Isotrop_Abs_Tot(:)
               Light_Tran_Tot(ipts,ivm,:)= direct_light_weight*Collim_Tran_Tot(:) + &
                    (un-direct_light_weight)*Isotrop_Tran_Tot(:)

               ! Notice that this light is actually cumulative, not per
               ! level!  This was needed for debugging purposes and running
               ! tests.  However, the photosynthesis routines need the light
               ! transmitted per level, i.e. if there is zero LAI
               ! in a level, the transmitted light will be one.
               DO ilevel=1,nlevels_tot-1
                  IF(Light_Tran_Tot(ipts,ivm,ilevel+1) .GT. alb_threshold)THEN
                     Light_Tran_Tot(ipts,ivm,ilevel)=&
                          Light_Tran_Tot(ipts,ivm,ilevel) / &
                          Light_Tran_Tot(ipts,ivm,ilevel+1)
                  ELSE

                     ! Here, we really don't know anything about how much 
                     ! light is transmitted in this layer, but there is no 
                     ! light reaching it from above so we can safely assume 
                     ! no photosynthesis takes place.  This is equivalent
                     ! to assuming it has no LAI, which means the 
                     ! transmission will be unity.
                     Light_Tran_Tot(ipts,ivm,ilevel)=un

                  ENDIF
               ENDDO

               ! Debug
               ! Notice that the sum of the transmissions may not equal one.  This
               ! is due to multiple scattering (light can be reflected up from a
               ! lower layer, and then back down).
               !-

            ENDIF ! IF ks==1

         ENDDO ! ks = 1, n_spectralbands 

      ENDDO ! ivm = 2, nvm  

   ENDDO ! ipts = 1, kjpindex 

   ! now we need to average the albedo over all our spectra, so we can 
   ! pass it to modules which do not distinguish between different 
   ! spectral bands
   albedo_glob(:) = zero
   DO ks = 1, n_spectralbands ! Loop over # of spectra
      albedo_glob(:) = albedo_glob(:) + albedo(:,ks)
   ENDDO
   albedo_glob(:)=albedo_glob(:)/REAL(n_spectralbands,r_std)

   !+++CHECK+++
   ! Variable dimensions xxx_Tot_mean are for a single pixel. Causing 1+1 problems
   ! when running a larger domain. Needs to be corrected when implementing
   ! a global use of the multi-layer energy budget.

   !! Debug ID12.
   !! For mleb purpose. Next step is to adapt those variables

   Light_Abs_Tot_mean(:) = Light_Abs_Tot_mean(:) / n_spectralbands
   Light_Alb_Tot_mean(:) = Light_Alb_Tot_mean(:) / n_spectralbands

   !! End ID12
   !+++++++++++

   ! Error checking
   IF (abs(albedo(1,1)) .gt. 1.0d0) THEN
       WRITE(numout,*) 'Albedo has to high values. Value higher than 1 has been detected.'
       WRITE(numout,*) 'albedo(1,1) = ', albedo(1,1)
       CALL print_debugging_albedo_info(1, cosine_sun_angle, &
             leaf_reflectance, leaf_transmittance, laieff_collim_1, &
             laieff_isotrop_1, br_base_temp_collim, br_base_temp_isotrop, &
             Collim_Alb_Tot_1, Collim_Tran_Tot_1, Collim_Abs_Tot_1, &
             Isotrop_Alb_Tot_1, Isotrop_Tran_Tot_1, Isotrop_Abs_Tot_1)

       CALL ipslerr_p(3,'albedo_surface_main','albedo(1,1) > 1','','')
   END IF
   !-

   ! Debug
   IF (printlev_loc>=4) THEN
         WRITE(numout,*) 'laieff_isotrop end albedo for test_pft',test_pft,'is', &
              laieff_isotrop(test_grid,:,test_pft)
   ENDIF
   !-

   !! 5. Output diagnostics
   !!    Return if current call is done from the initialization phase
   IF (firstcall_albedo_surface) RETURN

   ! Add XIOS default value where no sun
   DO ipts=1,kjpindex 
      IF(coszang(ipts) .LE. min_sechiba)THEN
         albedo_diag(ipts,:) = xios_default_val
         alb_bare_diag(ipts,:) = xios_default_val
         albedo_snow_diag(ipts,:) = xios_default_val
         albedo_noon(ipts,:) = xios_default_val
      ELSE
         albedo_diag(ipts,:) = albedo(ipts,:)
         alb_bare_diag(ipts,:) = alb_bare(ipts,:)
         albedo_snow_diag(ipts,:) = albedo_snow(ipts,:)

         ! On the Greenwich (lon=0) the time step corresponding to the solar noon is when sec_start-1 <= 43200 and sec_end+1 >= 43200
         ! (+1/-1 are added for cases when time step starts/ends precisely at noon)
         ! for other longitudes (-180 <= lon <= 180) solar noon is shifted by lon/360*86400
         ! 'modulo' function (residual after division) is used to make sure that time is never larger than day length 24*60*60=86400
         ! lalo(ipts,2) holds longitude of pixel
         IF ((MODULO(sec_start + lalo(ipts,2)/360.*86400., 86400.)-1 .LE. 43200.) .AND. (MODULO(sec_end + lalo(ipts,2)/360.*86400., 86400.)+1 .GE. 43200.)) THEN
            albedo_noon(ipts,:) = albedo(ipts,:)
         ELSE
            albedo_noon(ipts,:) = xios_default_val
         END IF
      END IF
   END DO

   IF (.NOT. impaze) THEN
      CALL xios_orchidee_send_field("soilalb_vis",alb_bare_diag(:,1))
      CALL xios_orchidee_send_field("soilalb_nir",alb_bare_diag(:,2))
   END IF
   CALL xios_orchidee_send_field("albedo_vis",albedo_diag(:,1))
   CALL xios_orchidee_send_field("albedo_nir",albedo_diag(:,2))
   CALL xios_orchidee_send_field("albedo_noon_vis",albedo_noon(:,1))
   CALL xios_orchidee_send_field("albedo_noon_nir",albedo_noon(:,2))
   CALL xios_orchidee_send_field("albedo_snow_vis",(albedo_snow_diag(:,1)))
   CALL xios_orchidee_send_field("albedo_snow_nir",(albedo_snow_diag(:,2)))
   
   IF ( almaoutput ) THEN
      CALL histwrite_p(hist_id, 'Albedo', kjit, (albedo(:,1) + albedo(:,2))/2, kjpindex, index)
      CALL histwrite_p(hist_id, 'SAlbedo', kjit, (albedo_snow(:,1)+albedo_snow(:,2))/2, kjpindex, index)
      IF ( hist2_id > 0 ) THEN
         CALL histwrite_p(hist2_id, 'Albedo', kjit, (albedo(:,1) + albedo(:,2))/2, kjpindex, index)
         CALL histwrite_p(hist2_id, 'SAlbedo', kjit, (albedo_snow(:,1)+albedo_snow(:,2))/2, kjpindex, index)
      ENDIF
   ELSE
      IF (.NOT. impaze) THEN
         CALL histwrite_p(hist_id, 'soilalb_vis', kjit, alb_bare(:,1), kjpindex, index)
         CALL histwrite_p(hist_id, 'soilalb_nir', kjit, alb_bare(:,2), kjpindex, index)
         IF ( hist2_id > 0 ) THEN
            CALL histwrite_p(hist2_id, 'soilalb_vis', kjit, alb_bare(:,1), kjpindex, index)
            CALL histwrite_p(hist2_id, 'soilalb_nir', kjit, alb_bare(:,2), kjpindex, index)
         ENDIF
      END IF
   ENDIF
   
   IF (printlev_loc>=3) WRITE(numout,*) 'Leaving albedo_surface_main'
   
 END SUBROUTINE albedo_surface_main


  !!  =============================================================================================================================
  !! SUBROUTINE		 		    : albedo_surface_finalize
  !!
  !>\BRIEF                                    Write to restart file
  !!
  !! DESCRIPTION			    : This subroutine writes the module variables and variables calculated in albedo
  !!                                          to restart file
  !!
  !! RECENT CHANGE(S)			    : None
  !!
  !! REFERENCE(S)			    : None
  !! 
  !! FLOWCHART                              : None
  !! \n
  !_ ==============================================================================================================================
 !+++CHECK+++
 ! Variable dimensions xxx_Tot_mean are for a single pixel. Causing 1+1 problems
 ! when running a larger domain. Needs to be corrected when implementing
 ! a global use of the multi-layer energy budget.
 SUBROUTINE albedo_surface_finalize (kjit,        kjpindex,            rest_id,                          &
                              albedo,              Isotrop_Abs_Tot_p,   Isotrop_Tran_Tot_p,               &
                              Light_Abs_Tot_mean, Light_Alb_Tot_mean, &
                              laieff_isotrop)
   !+++++++++++

    !! 0. Variable and parameter declaration
    !! 0.1 Input variables  
    INTEGER(i_std), INTENT(in)                                    :: kjit               !! Time step number 
    INTEGER(i_std), INTENT(in)                                    :: kjpindex           !! Domain size
    INTEGER(i_std),INTENT (in)                                    :: rest_id            !! Restart file identifier
    REAL(r_std), DIMENSION(kjpindex,n_spectralbands), INTENT(in)  :: albedo             !! Albedo (two stream radiation transfer model)
    REAL(r_std), DIMENSION(kjpindex,nvm,nlevels_tot), INTENT(in)  :: Isotrop_Abs_Tot_p  !! Absorbed radiation per layer for photosynthesis
    REAL(r_std), DIMENSION(kjpindex,nvm,nlevels_tot), INTENT(in)  :: Isotrop_Tran_Tot_p !! Transmitted radiation per layer for photosynthesis
    !+++CHECK+++
    ! Variable dimensions xxx_Tot_mean are for a single pixel. Causing 1+1 problems
    ! when running a larger domain. Needs to be corrected when implementing
    ! a global use of the multi-layer energy budget.
    REAL(r_std),DIMENSION(nlevels_tot), INTENT(in)                :: Light_Abs_Tot_mean !! total absorption for a given level * changed for enerbil integration
    REAL(r_std),DIMENSION(nlevels_tot), INTENT(in)                :: Light_Alb_Tot_mean !! total albedo for a given level * changed for enerbil integration
    !+++++++++++
    REAL(r_std), DIMENSION(kjpindex,nlevels_tot,nvm), INTENT(in)  :: laieff_isotrop     !! Leaf Area Index Effective converts
                                                                                        !! 3D lai into 1D lai for two stream

    !_ ================================================================================================================================
    
    IF ( alb_bg_modis ) THEN
       CALL restput_p (rest_id, 'bckgrnd_alb', nbp_glo, 2, 1, kjit, bckgrnd_alb, 'scatter',  nbp_glo, index_g)
    ELSE
       CALL restput_p (rest_id, 'soilalbedo_dry', nbp_glo, 2, 1, kjit, soilalb_dry, 'scatter',  nbp_glo, index_g)
       !-
       CALL restput_p (rest_id, 'soilalbedo_wet', nbp_glo, 2, 1, kjit, soilalb_wet, 'scatter',  nbp_glo, index_g)
       !-
       CALL restput_p (rest_id, 'soilalbedo_moy', nbp_glo, 2, 1, kjit, soilalb_moy, 'scatter',  nbp_glo, index_g)
    ENDIF
    
    CALL restput_p (rest_id, 'albedo' , nbp_glo, n_spectralbands, 1, kjit, albedo, "scatter", nbp_glo, index_g)
    !-
    CALL restput_p (rest_id, 'Isotrop_Abs_Tot_p', nbp_glo, nvm, nlevels_tot, kjit, Isotrop_Abs_Tot_p, "scatter", nbp_glo, index_g)
    !-
    CALL restput_p (rest_id, 'Isotrop_Tran_Tot_p', nbp_glo, nvm, nlevels_tot, kjit, Isotrop_Tran_Tot_p, "scatter", nbp_glo, index_g)

    !+++CHECK+++
    ! Variable dimensions xxx_Tot_mean are for a single pixel. Causing 1+1 problems
    ! when running a larger domain. Needs to be corrected when implementing
    ! a global use of the multi-layer energy budget.
    IF (ok_mleb .AND. is_root_prc) THEN
       CALL restput (rest_id, 'Light_Abs_Tot_mean', nlevels_tot, 1, 1, kjit, Light_Abs_Tot_mean)
       CALL restput (rest_id, 'Light_Alb_Tot_mean', nlevels_tot, 1, 1, kjit, Light_Alb_Tot_mean)
    END IF
    !+++++++++++

    CALL restput_p (rest_id, 'laieff_isotrop', nbp_glo, nlevels_tot, nvm, kjit, laieff_isotrop, "scatter", nbp_glo, index_g)


  END SUBROUTINE albedo_surface_finalize


!! ==============================================================================================================================
!! SUBROUTINE 	: albedo_surface_clear
!!
!>\BRIEF        Deallocate albedo variables
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): None 
!!
!! REFERENCES	: None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE albedo_surface_clear  ()

    IF (ALLOCATED (soilalb_dry)) DEALLOCATE (soilalb_dry)
    IF (ALLOCATED(soilalb_wet))  DEALLOCATE (soilalb_wet)
    IF (ALLOCATED(soilalb_moy))  DEALLOCATE (soilalb_moy)
    IF (ALLOCATED(bckgrnd_alb))  DEALLOCATE (bckgrnd_alb)
    
  END SUBROUTINE albedo_surface_clear


!! ==============================================================================================================================
!! SUBROUTINE   : twostream_solver
!!
!>\BRIEF        : Computes the two-stream albedo solution for a level with given single scatterer properties
!!                and a defined background albedo
!!
!! DESCRIPTION  : This solution is the two-stream solution of Pinty et al (2006) for a vegetation level above
!!                an isotropically reflecting background.  It breaks down the problem into three parts, solving
!!                each part for the case of diffuse (isotropic) and direct (collimated) light.  The three parts are,
!!                1) The term due to light that does not interact at all with the canopy (Black Canopy)
!!                2) Light that does not interact at all with the background (Black Background)
!!                3) Light that bounces between the background and the canopy
!!                This routine (and the routines it uses) were received directly from Bernard Pinty, and only some
!!                minor modifcations were made, in addition to more documentation
!! 
!!                NOTE: the single layer solution is no longer used. The code is kept here in case it is decided that
!!                      the multi layer solution is two expensive if only a single layer is used for the energy budget.
!!                      In that case a solution needs to be implemeted to calculate the light distribution within
!!                      the canopy. 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): Collim_Alb_Tot,Collim_Tran_Tot,
!!      Collim_Abs_Tot,Isotrop_Alb_Tot,Isotrop_Tran_Tot,Isotrop_Abs_To
!!
!! REFERENCE(S) :  B. Pinty, T. Lavergne, R.E. Dickinson, J-L. Widlowski, N. Gobron and M. M. Verstraete (2006).
!! Simplifying the Interaction of Land Surfaces with Radiation for Relating Remote Sensing Products to Climate Models.
!! Journal of Geophysical Research. Vol 111, D02116.
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

 SUBROUTINE twostream_solver(leaf_reflectance, leaf_transmittance, background_reflectance_collim, background_reflectance_isotrop, &
          cosine_sun_angle, Collim_Alb_Tot, Collim_Tran_Tot, Collim_Abs_Tot, &
          Isotrop_Alb_Tot, Isotrop_Tran_Tot, Isotrop_Abs_Tot, laieff_collim, &
           laieff_isotrop,  Collim_Tran_Uncollided,  Isotrop_Tran_Uncollided)


   !! 0. Variables and parameter declaration
   !! 0.1 Input variables
   REAL(r_std), INTENT(IN)    :: leaf_reflectance            !! effective leaf reflectance, between 0 and 1
                                                             !! @tex $()$ @endtex 
   REAL(r_std), INTENT(IN)    :: leaf_transmittance          !! effective leaf transmittance, between 0 and 1
                                                             !! @tex $()$ @endtex 
   REAL(r_std), INTENT(IN)    :: background_reflectance_collim      !! background reflectance for direct radiation,
                                                             !! between 0 and 1
   REAL(r_std), INTENT(IN)    :: background_reflectance_isotrop      !! background reflectance for diffuse radiation,
                                                             !! between 0 and 1
   REAL(r_std), INTENT(IN)    :: cosine_sun_angle            !! cosine of the solar zenith angle, between 0 and 1
                                                             !! @tex $()$ @endtex 
   REAL(r_std), INTENT(IN)    :: laieff_collim               !! Effective Leaf Area Index, computed at current sun angle
   REAL(r_std), INTENT(IN)    :: laieff_isotrop              !! Effective Leaf Area Index, computed at 60 degrees
                                                             !! @tex $(m^{2} m^{-2})$ @endtex 

   !! 0.2 Output variables
   !! notice that these variables (absorption + transmission + reflection) do not necessarily add up to 1 due to multiple scattering
   REAL(r_std), INTENT(OUT)   :: Collim_Alb_Tot              !! collimated total albedo from this level, between 0 and 1
                                                             !! @tex $()$ @endtex 
   REAL(r_std), INTENT(OUT)   :: Collim_Tran_Tot             !! collimated total transmission through this level, between 0 and 1
                                                             !! @tex $()$ @endtex 
   REAL(r_std), INTENT(OUT)   :: Collim_Abs_Tot              !! collimated total absorption by this level, between 0 and 1
                                                             !! @tex $()$ @endtex 
   REAL(r_std), INTENT(OUT)   :: Isotrop_Alb_Tot             !! isotropic total albedo from this level, between 0 and 1
                                                             !! @tex $()$ @endtex 
   REAL(r_std), INTENT(OUT)   :: Isotrop_Tran_Tot            !! isotropic total transmission through this level, between 0 and 1
                                                             !! @tex $()$ @endtex 
   REAL(r_std), INTENT(OUT)   :: Isotrop_Abs_Tot             !! isotropic total absorption by this level, between 0 and 1
                                                             !! @tex $()$ @endtex 
   REAL(r_std), INTENT(OUT)   :: Collim_Tran_Uncollided      !! collimated uncollied transmission through this level, between 0 and 1
                                                             !! @tex $()$ @endtex 
   REAL(r_std), INTENT(OUT)   :: Isotrop_Tran_Uncollided     !! isotropic uncollied transmission through this level, between 0 and 1
                                                             !! @tex $()$ @endtex 

   !! 0.3 Modified variables

   !! 0.4 Local variables
   LOGICAL :: ok
   REAL(r_std), DIMENSION(4)                           :: gammaCoeffs
   REAL(r_std), DIMENSION(4)                           :: gammaCoeffs_star
   REAL(r_std)                                         :: tauprimetilde
   REAL(r_std)                                         :: tauprimestar
   REAL(r_std)                                         :: sun_zenith_angle_radians
   REAL(r_std), PARAMETER                              :: isotropic_cosine_constant=0.5/0.705

   !calculated fluxes

   REAL(r_std)                :: Collim_Alb_BB 
   REAL(r_std)                :: Collim_Tran_BB
   REAL(r_std)                :: Collim_Abs_BB
   REAL(r_std)                :: Isotrop_Alb_BB
   REAL(r_std)                :: Isotrop_Tran_BB
   REAL(r_std)                :: Isotrop_Abs_BB
   REAL(r_std)                :: Collim_Tran_BC
   REAL(r_std)                :: Isotrop_Tran_BC
   REAL(r_std)                :: Collim_Tran_TotalOneWay
   REAL(r_std)                :: Isotrop_Tran_TotalOneWay
   REAL(r_std)                :: Below_reinject_rad_collim,Below_reinject_rad_isotrop

!_ ================================================================================================================================

   IF (printlev_loc>=3) WRITE(numout,*) 'Entering twostream_solver'

   ! convert angular values
   
   sun_zenith_angle_radians = acos(cosine_sun_angle)

   ! calculate the 4 gamma coefficients both in isotropic and collimated illumination 

   call gammas(leaf_reflectance, leaf_transmittance, cosine_sun_angle, gammaCoeffs)
   call gammas(leaf_reflectance,leaf_transmittance,isotropic_cosine_constant,&
        gammaCoeffs_star)

   ! estimate the effective value of the optical thickness 

   tauprimetilde = 0.5_r_std * laieff_collim
   !tauprimetilde = 1.

   ! Should become zenit angle dependent (=60°) calculated by the Monte Carlo 
   ! photon shooting routine

   tauprimestar  = 0.5_r_std * laieff_isotrop
   !tauprimestar  = 1.

   ! +++++++++++++ BLACK BACKGROUND ++++++++++++++++++++++ 
   ! * Apply the black-background 2stream solution 
   ! * These equations are written for the part of the incoming radiation 
   ! * that never hits the background but does interact with the vegetation 
   ! * canopy.
   ! * 
   ! * Note : the same routine dhrT1() is used for both the isotropic and
   ! * collimated illumination conditions but the calling arguments differ.
   ! * (especially the solar angle used).
   ! */

   !    /* 1) collimated source */
   ok=dhrT1(leaf_reflectance,leaf_transmittance,&
        gammaCoeffs(1),gammaCoeffs(2),gammaCoeffs(3),gammaCoeffs(4),&
        sun_zenith_angle_radians, tauprimetilde,&
        Collim_Alb_BB,Collim_Tran_BB,Collim_Abs_BB)                      
   !    /* 2) isotropic source */
   ok=dhrT1(leaf_reflectance,leaf_transmittance,&
        gammaCoeffs_star(1),gammaCoeffs_star(2),gammaCoeffs_star(3),&
        gammaCoeffs_star(4),&
        acos(isotropic_cosine_constant),tauprimestar,&
        Isotrop_Alb_BB,Isotrop_Tran_BB,Isotrop_Abs_BB)


   ! +++++++++++++ BLACK CANOPY ++++++++++++++++++++++ 
   ! * Apply the black canopy solution.
   ! * These equations hold for the part of the incoming radiation 
   ! * that do not interact with the vegetation, travelling through 
   ! * its gaps.
   ! */

   !    /* 1) collimated source */
   IF (cosine_sun_angle .NE. 0) THEN
      Collim_Tran_BC  = exp( - tauprimetilde/cosine_sun_angle)
      Collim_Tran_Uncollided = Collim_Tran_BC
   ENDIF


   !    /* 2) isotropic source */
   Isotrop_Tran_BC = TBarreUncoll_exact(tauprimestar)
   Isotrop_Tran_Uncollided = Isotrop_Tran_BC


   ! /* Total one-way transmissions:
   ! * The vegetation canopy is crossed (one way) by the uncollided radiation
   ! * (black canopy) and the collided one (black background). */

   !    /* 1) collimated source */
   Collim_Tran_TotalOneWay  = Collim_Tran_BC  + Collim_Tran_BB 

   !    /* 2) isotropic source */
   Isotrop_Tran_TotalOneWay = Isotrop_Tran_BC + Isotrop_Tran_BB 

   ! * The Below_reinject_rad describes the process of reflecting toward the 
   ! * background the upward travelling radiation (re-emitted from below the 
   ! * canopy). It appears in the coupling equations as the limit of the series: 
   ! *    1 + rg*rbv + (rg*rbv)^2 + (rg*rbv)^3 + ...
   ! *      where rg is the background_reflectance and rbv is Isotrop_Alb_BB
   ! *      (with Isotrop describing the Lambertian reflectance of the background).
   ! */
   ! this might be involved in Eq. 27, 
   Below_reinject_rad_collim = un / (un - background_reflectance_collim*Isotrop_Alb_BB)
   Below_reinject_rad_isotrop = un / (un - background_reflectance_isotrop*Isotrop_Alb_BB)

   !/* TOTAL ALBEDO */
   !    /* 1) collimated source */
   Collim_Alb_Tot = Collim_Alb_BB + &
        background_reflectance_collim * Collim_Tran_TotalOneWay * &
        Isotrop_Tran_TotalOneWay * Below_reinject_rad_collim
   !    /* 2) isotropic source */
   Isotrop_Alb_Tot = Isotrop_Alb_BB + & 
        background_reflectance_isotrop * Isotrop_Tran_TotalOneWay * &
        Isotrop_Tran_TotalOneWay * Below_reinject_rad_isotrop ! this seems to be
   ! exactly Eq. 33

   !/* TOTAL TRANSMITION TO THE BACKGROUND LEVEL */
   !    /* 1) collimated source */
   Collim_Tran_Tot = Collim_Tran_TotalOneWay * Below_reinject_rad_collim ;

   !    /* 2) isotropic source */
   Isotrop_Tran_Tot = Isotrop_Tran_TotalOneWay * Below_reinject_rad_isotrop ;

   !/* TOTAL ABSORPTION BY THE VEGETATION LEVEL */
   !    /* 1) collimated source */
   Collim_Abs_Tot = un - (Collim_Tran_Tot + Collim_Alb_Tot) + &
        background_reflectance_collim * Collim_Tran_Tot;
   !    /* 2) isotropic source */
   Isotrop_Abs_Tot = un - (Isotrop_Tran_Tot + Isotrop_Alb_Tot) + &
        background_reflectance_isotrop * Isotrop_Tran_Tot;

!+++ TEMP +++ 
! Combine the collimated & isotrop together as Collimate for energybudget 
   Collim_Alb_Tot = Collim_Alb_Tot + Isotrop_Alb_Tot
   Collim_Abs_Tot = Collim_Abs_Tot + Isotrop_Abs_Tot   
!+++ TEMP +++       
        

    IF (printlev_loc>=3) WRITE(numout,*) 'Exiting twostream_solver'

 END SUBROUTINE twostream_solver


!! ============================================================================\n
!! FUNCTION     : dhrT1
!!
!>\BREF         
!!
!! DESCRIPTION : Taken directly from the Fortan code provided by Pinty et al for their \n
!!               two-stream model.  The only changes are the REAL types.  This is used
!!               to compute the black background absorption, transmission, and reflectance
!!               in Pintys scheme and it's based on the older model of Meador and Weaver...
!!               Pinty 2006 also includes a short discussion of it in Appendix A 
!!
!! RECENT CHANGE(S): None\n
!!
!! RETURN VALUE : dhrT1
!!
!! REFERENCE(S) : Meador and Weaver, 'Two-stream approximations to radiative transfer in
!!    planetary atmosphere: a unified description of existing methods and a new improvement'
!!    J. Atmospheric Sciences, VOL 37, p. 630--643, 1980.
!!
!! FLOWCHART    : None
!_ =============================================================================

FUNCTION dhrT1(rl,tl,gamma1,gamma2,gamma3,gamma4,tta,tau,AlbBS,Tdif,AbsVgt)


  !! 0. Variables and parameter declaration
  !! 0.1 Input variables
  REAL(r_std), INTENT(IN) :: rl  ! the effective reflectance of a single scatterer, between 0 and 1
                                 !! @tex $()$ @endtex      
  REAL(r_std), INTENT(IN) :: tl  ! the effective transmittance of a single scatterer, between 0 and 1
                                 !! @tex $()$ @endtex      
  REAL(r_std), INTENT(IN) :: gamma1 ! the gamma coefficients from Meador and Weaver
  REAL(r_std), INTENT(IN) :: gamma2
  REAL(r_std), INTENT(IN) :: gamma3
  REAL(r_std), INTENT(IN) :: gamma4
  REAL(r_std), INTENT(IN) :: tta ! the solar zenith angle, between 0 and pi/2
                                 !! @tex $(radians)$ @endtex      

  REAL(r_std), INTENT(IN) :: tau ! Effective LAI * G(theta)
                                 !! @tex $()$ @endtex      

  !! 0.2 Output variables
  REAL(r_std), INTENT(OUT) :: AlbBS ! the albedo of level, between 0 and 1
                                 !! @tex $()$ @endtex      
  REAL(r_std), INTENT(OUT) :: Tdif ! the transmitted light (all diffuse) from this light source, between 0 and 1
                                 !! @tex $()$ @endtex      
  REAL(r_std), INTENT(OUT) :: AbsVgt ! the light absorbed by the level, between 0 and 1
                                 !! @tex $()$ @endtex      
  LOGICAL :: dhrT1 ! the output flag, which always appears to be true

  !! 0.3 Modified variables
  !! 0.4 Local variables

  REAL(r_std) :: alpha1,alpha2,ksquare,k
  REAL(r_std) :: first_term,secnd_term1,secnd_term2,secnd_term3
  REAL(r_std) :: expktau,Tdir
  REAL(r_std) :: mu0,w0 

  !_ ================================================================================================================================

  mu0=cos(tta)
  w0=rl+tl


  Tdir = exp(-tau/mu0) ! direct transmission

  ! There is a difference between conservative and non-conservative scattering 
  ! conditions */
  IF (w0 .ne. 1.0 .AND. w0 .ne. 0.0) THEN
     !NON_CONSERVATIVE SCATTERING

     ! some additional parameters, taken from Eq. 16--18 of Meador and Weaver, J. Atmospheric Sciences,
     ! Vol 37, p. 630--643 (1980)

     alpha1  = gamma1*gamma4 + gamma2*gamma3
     alpha2  = gamma1*gamma3 + gamma2*gamma4
     ksquare = gamma1*gamma1 - gamma2*gamma2
     k       = sqrt(ksquare)

     expktau = exp(k*tau)

     !Black Soil Albedo...Eq. 14 in Meador and Weaver
     first_term  = ((un-ksquare*mu0*mu0)*((k+gamma1)*expktau + (k-gamma1)/expktau))
     IF (first_term .eq. 0.0) THEN
        !we will be dividing by zero : cannot continue.
        dhrT1 = .false.
     ELSE 
        first_term = un/first_term
        secnd_term1 = (un - k*mu0)*(alpha2 + k*gamma3)*expktau
        secnd_term2 = (un + k*mu0)*(alpha2 - k*gamma3)/expktau
        secnd_term3 = 2.0_r_std * k * (gamma3 - alpha2*mu0)*Tdir

        AlbBS = (w0 * first_term * (secnd_term1 - secnd_term2 - secnd_term3))

        !Transmission...Eq. 15 in Meador and Weaver, for diffuse light?
        IF (ksquare .eq. 0.0) THEN 
           first_term = un
        ENDIF
        secnd_term1 = (un+k*mu0)*(alpha1+k*gamma4)*expktau
        secnd_term2 = (un-k*mu0)*(alpha1-k*gamma4)/expktau
        secnd_term3 = 2.0_r_std * k * (gamma4 + alpha1*mu0)
        Tdif = - w0*first_term*(Tdir*(secnd_term1 - secnd_term2) - secnd_term3)

        ! Absorption by vegetation...whatever is not transmitted or reflected 
        ! must be absorbed
        AbsVgt = (un- (Tdif+Tdir) - AlbBS)
     ENDIF ! first_term .eq. 0.0
  ELSE IF (w0 .eq. 0.) THEN
     !BLACK CANOPY
     AlbBS = zero
     Tdif  = zero
     AbsVgt = un - Tdir
  ELSE
     !CONSERATIVE SCATTERING...Eq. 24 in Meador and Weaver
     AlbBS =  (un/(un + gamma1*tau))*(gamma1*tau + (gamma3-gamma1*mu0)*&
          (un-exp(-tau/mu0)));
     Tdif   = un - AlbBS - Tdir;
     AbsVgt = zero;
  ENDIF ! w0 .ne. 1.0 .AND. w0 .ne. 0.0

! not sure what the purpose of this flag is, as it will always be set to true here
  dhrT1 = .true.

END FUNCTION dhrT1

!! ============================================================================\n
!! FUNCTION     : TBarreUncoll_exact
!!
!>\BRIEF          Computes the transmission of the diffuse black canopy light
!!
!! DESCRIPTION : Taken directly from the Fortan code provided by Pinty et al for their \n
!!               two-stream model.  The only changes are the REAL types.  This appears
!!               to be solving Eq. 16 in Pinty 2006, which is the transmission which
!!               does not collide with the canopy
!!
!! RECENT CHANGE(S): None\n
!!
!! RETURN VALUE : TBarreUncoll_exact
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!_ =============================================================================

FUNCTION TBarreUncoll_exact(tau)


  !! 0. Variables and parameter declaration

  !! 0.1 Input variables
  REAL(r_std), INTENT(IN) :: tau ! Effective LAI * G(theta)
  !! @tex $()$ @endtex      

  !! 0.2 Output variables
  REAL(r_std) :: TBarreUncoll_exact ! the isotropic light transmission uncollided
                                    ! with the canopy, between 0 and 1
  !! @tex $()$ @endtex      


  !! 0.3 Modified variables

  !! 0.4 Local variables
!  INTEGER :: j,ind
!  REAL(r_std) :: iGammaloc
!  INTEGER :: order

!_ ================================================================================================================================

  !+++++ CHECK +++++

!!$  iGammaloc = zero
!!$  order=20
!!$
!!$  ! in the case where tau is equal to zero, this crashes in the first loop where
!!$  ! iGammaloc is also zero...quick fix, and might even be accurate
!!$
!!$  IF(tau .GT. 1e-10_r_std)THEN
!!$
!!$     DO j=0,order-1
!!$        ind = order - j
!!$        iGammaloc = ind / (un + ind/(tau+iGammaloc))
!!$     END DO
!!$     iGammaloc=un/(tau + iGammaloc)
!!$  ENDIF
!!$
!!$  TBarreUncoll_exact = exp(-tau)*(un - tau + tau*tau*iGammaloc)

  ! this is a change suggested by Bernard to improve the matching between the one
  ! level and multilevel case

  TBarreUncoll_exact = exp(-tau*2.0_r_std*0.705_r_std)
  !++++++++++

END FUNCTION TBarreUncoll_exact

!! ==============================================================================================================================
!! SUBROUTINE   : gammas
!!
!>\BRIEF         Computes a set of gamma coefficients for use in the black background equations
!!
!! DESCRIPTION : Taken directly from the Fortan code provided by Pinty et al for their \n
!!               two-stream model.  The only changes are the REAL types.  This calculates
!!               the gamma coefficients based on Appendix A in the paper.  This seems to
!!               make the assumption of the Ross's G function and a spherical leaf
!!               angle distribution.
!!
!! RECENT CHANGE(S): None\n
!!
!! MAIN OUTPUT VARIABLE(S): gammaCoeff
!!
!! REFERENCE(S) :  B. Pinty, T. Lavergne, R.E. Dickinson, J-L. Widlowski, N. Gobron and M. M. Verstraete (2006).
!! Simplifying the Interaction of Land Surfaces with Radiation for Relating Remote Sensing Products to Climate Models.
!! Journal of Geophysical Research. Vol 111, D02116.
!!
!! FLOWCHART    : None
!_ ================================================================================================================================

SUBROUTINE gammas(rl, tl, mu0, gammaCoeff)


  !! 0. Variables and parameter declaration
  !! 0.1 Input variables
  REAL(r_std), INTENT(IN) :: rl  ! the effective reflectance of a single scatterer, between 0 and 1
  !! @tex $()$ @endtex      
  REAL(r_std), INTENT(IN) :: tl  ! the effective transmittance of a single scatterer, between 0 and 1
  !! @tex $()$ @endtex      
  REAL(r_std), INTENT(IN) :: mu0 ! the cosine of the solar zenith angle, between 0 and 1
  !! @tex $()$ @endtex      

  !! 0.2 Output variables
  REAL(r_std), DIMENSION(1:4), INTENT(OUT) :: gammaCoeff ! the set of gamma coefficients in the above reference
  !! @tex $()$ @endtex      

  !! 0.3 Modified variables

  !! 0.4 Local variables
  REAL(r_std) :: wd,w0,w0half,wdsixth

!_ ================================================================================================================================


  w0      = rl+tl
  wd      = rl-tl
  w0half  = w0*0.5_r_std
  wdsixth = wd/6.0_r_std

  gammaCoeff(1)=2._r_std*(un - w0half + wdsixth)
  gammaCoeff(2)=2._r_std*(w0half + wdsixth)
  gammaCoeff(3)=2._r_std*(w0half*0.5_r_std + mu0*wdsixth)/w0
  gammaCoeff(4)=un-gammaCoeff(3)

END SUBROUTINE gammas


!! ==============================================================================================================================
!! SUBROUTINE   : calculate_snow_albedo
!!
!>\BRIEF         Computes some of the information needed to calculate the effect of the snow albedo
!!               on the background reflectance in the two stream model. This is done
!!               with the snow albedo scheme from Krinner et al 2005. The calculation of
!!               snow cover fraction is taken from Yang et al. 1997
!!
!! DESCRIPTION : In order to compute the background albedo in Pinty's two stream model, we have to take
!!               into account any snow that has fallen through the canopy and landed on the ground.  In particular,
!!               we need the amount of ground covered by the snow and the albedo of this snow.  Both of these
!!               quantities are calculated here, using a function that changes the albedo of the snow based on
!!               its age. 
!!
!! RECENT CHANGE(S): None\n
!!
!! MAIN OUTPUT VARIABLE(S): ::snowa_veg, ::frac_snow_veg, ::albedo_snow
!!
!! REFERENCE(S) :  
!!                 Yang Z, Dickinson R, Robock A, Vinnikov K (1997) Validation of the snow submodel of the 
!!                 biosphere-atmosphere transfer scheme with Russian 
!!                 snow cover and meteorological observational data. Journal of Climate, 10, 353–373.
!!
!! FLOWCHART    : None
!_ ================================================================================================================================

SUBROUTINE calculate_snow_albedo(kjpindex,  coszang,  snow, &
           snow_age,  snow_nobio_age,  frac_nobio,  albedo_snow, &
           snowa_veg,  frac_snow_veg,  snowa_nobio,  frac_snow_nobio, &
           veget_max, z0m, veget)

  !! 0. Variables and parameter declaration
  !! 0.1 Input variables
   INTEGER,INTENT(in)                                  :: kjpindex        !! Domain size - Number of land pixels  (unitless) 
   REAL(r_std), DIMENSION(kjpindex), INTENT(in)        :: coszang         !!  cosine of the solar zenith angle, between 0 and 1
                                                                          !! @tex $()$ @endtex 
   REAL(r_std),DIMENSION (kjpindex), INTENT(in)        :: snow            !! Snow mass in vegetation (kg m^{-2})           
   REAL(r_std),DIMENSION (kjpindex), INTENT(in)        :: snow_age        !! Snow age (days)        
   REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT(in) :: snow_nobio_age  !! Snow age on continental ice, lakes, etc. (days)    
   REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT(in) :: frac_nobio      !! Fraction of non-vegetative surfaces, i.e. 
                                                                          !!    continental ice, lakes, etc. (unitless)     
   REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(in)    :: veget_max       !! PFT coverage fraction of a PFT (= ind*cn_ind) 
                                                                          !!    (m^2 m^{-2})
   REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(in)    :: veget           !! Fraction of vegetation types
   REAL(r_std), DIMENSION(kjpindex), INTENT(in)        :: z0m             !! Roughness height for momentum of vegetated part (m)
   REAL(r_std), DIMENSION(:),INTENT(in)                :: frac_snow_veg   !! The fraction of the surface for each PFT covered by snow
   REAL(r_std), DIMENSION(:,:),INTENT(in)              :: frac_snow_nobio !! The fraction of nonbiological land types covered by snow

   !! 0.2 Output variables
   REAL(r_std), DIMENSION (kjpindex,n_spectralbands), &
                               INTENT (out)        :: albedo_snow         !! Snow albedo (unitless ratio)   
   REAL(r_std), DIMENSION(kjpindex,nvm,n_spectralbands),&
                                INTENT(out)        :: snowa_veg           !! the albedo of snow...this is seperated into
                                                                          !! PFT types, even though the calculation is identical
                                                                          !! for all PFTs right now
   REAL(r_std), DIMENSION(kjpindex,nnobio,n_spectralbands),&
                                INTENT(out)        :: snowa_nobio         !! the albedo of snow on nonbiological land types


  !! 0.3 Modified variables

  !! 0.4 Local variables
   REAL(r_std), DIMENSION(kjpindex)                      :: agefunc_veg
   REAL(r_std), DIMENSION(kjpindex,nnobio)               :: agefunc_nobio
   REAL(r_std)                                           :: snow_alb_direct, snow_alb_diffuse,f_mu
   REAL(r_std), DIMENSION(kjpindex)                      :: fraction_veg          !! fraction of the pixel with vegetation (bare soil is considered vegetation)
   REAL(r_std), DIMENSION(kjpindex)                      :: fraction_temp         !! Temporary fraction for snow albedo averaging
   
   INTEGER :: ipts,ks,ivm,jv ! indices

   REAL(r_std), DIMENSION (nvm,2)                        :: snowa_aged_tmp        !! spectral domains (unitless) 
   REAL(r_std), DIMENSION (nvm,2)                        :: snowa_dec_tmp
   REAL(r_std), DIMENSION (2)                            :: snowa_aged_nobio_tmp  !! spectral domains (unitless) 
   REAL(r_std), DIMENSION (2)                            :: snowa_dec_nobio_tmp
  
  !_ ================================================================================================================================

  ! initialize the output
  albedo_snow(:,:) = zero
  snowa_veg(:,:,:) = zero
  snowa_nobio(:,:,:) = val_exp
  fraction_veg(:) = un - SUM(frac_nobio(:,:),2)
  fraction_temp(:) = zero
  snowa_aged_tmp(:,ivis) = snowa_aged_vis(:)
  snowa_aged_tmp(:,inir) = snowa_aged_nir(:)
  snowa_dec_tmp(:,ivis) = snowa_dec_vis(:)
  snowa_dec_tmp(:,inir) = snowa_dec_nir(:)
  snowa_aged_nobio_tmp(ivis) = snowa_aged_nobio_vis
  snowa_aged_nobio_tmp(inir) = snowa_aged_nobio_nir
  snowa_dec_nobio_tmp(ivis) = snowa_dec_nobio_vis
  snowa_dec_nobio_tmp(inir) = snowa_dec_nobio_nir

  DO ipts = 1, kjpindex ! loop over all the grid squares

     IF (SUM(veget_max(ipts,:)).LT.min_stomate) THEN
        ! The pixel is inconsistent. It is present in climate forcing
        ! but it is absent in the PFT map. If we do not want the model
        ! to crash or generate NaNs. Some variables need to be initialized
        agefunc_nobio(ipts,:) = un
     END IF
     
     DO ivm=1,nvm

        IF (ivm.NE.ibare_sechiba.AND.veget_max(ipts,ivm) == zero) THEN
           ! No vegetation so no reason to do the calculation; but for bare soil always some contrib (either from BS or from veget)
           CYCLE
        ENDIF

        DO ks=1,n_spectralbands ! loop over spectra

           ! The snow albedo could be either prescribed or 
           ! calculated following Chalita and Treut (1994)
           
           ! Check if the precribed value fixed_snow_albedo exists. 
           IF (ABS(fixed_snow_albedo - undef_sechiba) .GT. EPSILON(undef_sechiba)) THEN

              snowa_veg(ipts,ivm,ks) = fixed_snow_albedo
              snowa_nobio(ipts,ivm,ks) = fixed_snow_albedo
              
           ELSE

              ! Calculate age dependence
              ! On vegetated surfaces
              agefunc_veg(ipts) = EXP(-snow_age(ipts)/tcst_snowa)
                       

              ! Albedo of the snow
              ! under the vegetation. These values are used to calculate
              ! the background albedo which in turn is used in the
              ! two-stream solver to calculate the radiative transfer
              ! through the canopy
              snowa_veg(ipts,ivm,ks) = zero
              

              ! Same calculation of the snow albedo for each PFT but with different snow aging parameters
              ! Note that the treatment of the snow differs between vegetation type, with different paramters
              ! i) Short vegetation (grass, crop): we assume that when there
              ! are no leaves (LAI=0), the snow aging is not impacted by
              ! the vegetation. The difference between veget and vegetmax
              ! (bare soil fraction within the PFT) should be considered in
              ! the same way as the bare soil fraction (PFT1).
              ! ii) High vegetation (trees): In this case we consider that
              ! the vegetation impacts the snow aging even when LAI=0
              ! because of the trunks, branches. So these PFTs should be
              ! treated differently without any consideration of the bare
              ! soil fraction within the PFT.
              IF( fraction_veg(ipts) .GT. min_sechiba ) THEN
                    snowa_veg(ipts,ivm,ks) =  snowa_aged_tmp(ivm,ks) + &
                         snowa_dec_tmp(ivm,ks) * agefunc_veg(ipts) 
              ENDIF
              
           ENDIF ! prescribe or calculate albedo
           

           ! albedo snow (diagnostic) : snowa_veg weighted by the fraction of all PFTs (simple approach that works for now 
           ! given that the fraction of snow cover is the same for all PFTs (including bare soil);
           albedo_snow(ipts,ks) = albedo_snow(ipts,ks) + snowa_veg(ipts,ivm,ks) * &
                veget_max(ipts,ivm) * frac_snow_veg(ipts)

           ! Sum the contributing fraction for the later calculation of the mean snow albedo
           IF (ks .EQ. 1) THEN 
              fraction_temp(ipts)  = fraction_temp(ipts) + veget_max(ipts,ivm)*frac_snow_veg(ipts)
           ENDIF
           
        ENDDO ! ks=1,n_spectralbands 

     ENDDO ! ivm=1,nvm 

     ! NIR and VIS albedo for the snow on Nobio surfaces
     IF (ABS(fixed_snow_albedo - undef_sechiba) .GT. EPSILON(undef_sechiba)) THEN
        
        snowa_nobio(ipts,:,:) = fixed_snow_albedo
        
     ELSE
     
        ! On non-vegtative surfaces
        DO jv = 1, nnobio ! Loop over # nobio types
           agefunc_nobio(ipts,jv) = EXP(-snow_nobio_age(ipts,jv)/tcst_snowa_nobio)
        ENDDO

        DO ks = 1, n_spectralbands
           
           DO jv = 1, nnobio
             
              ! The albedo due to snow of this age on this nonbio surface
              ! This value is used in albedo_main to calculate
              ! the albedo of the entire pixel (= vegetation + nobio + bare soil)
              snowa_nobio(ipts,jv,ks) = ( snowa_aged_nobio_tmp(ks) + &
                   snowa_dec_nobio_tmp(ks) * agefunc_nobio(ipts,jv) )
              
           ENDDO
           
        ENDDO
        
     ENDIF ! prescribe or calculate

     ! NIR and VIS albedo for the snow on non-biological surfaces
     DO ks = 1, n_spectralbands
        
        DO jv = 1, nnobio 

           ! This value is only used as an output variable
           albedo_snow(ipts,ks) = albedo_snow(ipts,ks) + &
                ( frac_nobio(ipts,jv) ) * ( frac_snow_nobio(ipts,jv) ) * &
                snowa_nobio(ipts,jv,ks) 

           ! Sum the contributing fraction for the later calculation of the mean snow albedo
           IF (ks .EQ. 1) THEN 
              fraction_temp(ipts)  = fraction_temp(ipts) + &
                   frac_nobio(ipts,jv) * frac_snow_nobio(ipts,jv) 
           ENDIF
           
        ENDDO ! jv = 1, nnobio

        ! Final normalisation for the calculation of the mean snow albedo
        IF (fraction_temp(ipts) .GT. min_sechiba) THEN
           albedo_snow(ipts,ks) = albedo_snow(ipts,ks) / fraction_temp(ipts)
        ENDIF
        
     ENDDO ! spectralbands 
    
  ENDDO ! ipts = 1, kjpindex 

END SUBROUTINE calculate_snow_albedo

!! ==============================================================================================================================
!! SUBROUTINE   : optimize_albedo_values
!!
!>\BRIEF         Follow the radiation scattered through the canopy in the 
!!               case of multiple levels.
!!
!! DESCRIPTION : Right now, we know that using Pintys method in the case of a 
!!               single canopy level will result in different top of the canopy
!!               albedos than iwhat is found if we use multiple canopy levels.
!!               We trust that the single level case gives the 'true' values.
!!
!!               This algorithm follows light as it scatters through multiple
!!               levels in the canopy.  At each level, the scattering 
!!               solution is found by solving Pinty's two stream model.  This
!!               means that light which enters a level can either be transmitted
!!               without colliding with the vegetation, transmitted after
!!               collision with the vegetation, reflecting off the vegetation,
!!               or being absorbed.  We follow the fluxes until they get
!!               really small.
!!
!! RECENT CHANGE(S): None\n
!!
!! MAIN OUTPUT VARIABLE(S): ::lconverged, ::Collim_Alb_Tot, ::Collim_Tran_Tot, 
!!          ::Collim_Abs_Tot, ::Isotrop_Alb_Tot, ::Isotrop_Tran_Tot, ::Isotrop_Abs_Tot
!!
!! REFERENCE(S) :  B. Pinty, T. Lavergne, R.E. Dickinson, J-L. Widlowski, N. Gobron
!!          and M. M. Verstraete (2006). Simplifying the Interaction of Land 
!!          Surfaces with Radiation for Relating Remote Sensing Products to 
!!          Climate Models. Journal of Geophysical Research. Vol 111, D02116.
!!
!! FLOWCHART    : None
!_ ================================================================================================================================
SUBROUTINE multilevel_albedo(cosine_sun_angle, leaf_single_scattering_albedo_start,&
     leaf_psd_start, br_base_temp_collim,  br_base_temp_isotrop, &
     laieff_collim_pft, laieff_isotrop_pft, lconverged, &
     Collim_Alb_Coll, Collim_Tran_Coll, Collim_Abs_Tot, Isotrop_Alb_Coll, &
     Isotrop_Tran_Coll, Isotrop_Abs_Tot, Collim_Tran_Uncoll, Isotrop_Tran_Uncoll, lprint)

  !! 0. Variables and parameter declaration
  !! 0.1 Input variables
  REAL(r_std), INTENT(IN)    :: cosine_sun_angle                          !! cosine of the solar zenith angle, between 0 and 1
                                                                          !! @tex $()$ @endtex 
  REAL(r_std), INTENT(IN)    :: leaf_single_scattering_albedo_start       !! cosine of the solar zenith angle, between 0 and 1
  REAL(r_std), INTENT(IN)    :: leaf_psd_start                            !! cosine of the solar zenith angle, between 0 and 1
  REAL(r_std), INTENT(IN)    :: br_base_temp_collim                       !! cosine of the solar zenith angle, between 0 and 1
  REAL(r_std), INTENT(IN)    :: br_base_temp_isotrop                      !! cosine of the solar zenith angle, between 0 and 1
  REAL(r_std),DIMENSION(:),INTENT(IN)        :: laieff_collim_pft         !! Effective lai for a single pixel and pft to be 
  REAL(r_std),DIMENSION(:),INTENT(IN)        :: laieff_isotrop_pft        !! Effective lai for a single pixel and pft to be 
  LOGICAL,INTENT(IN)        :: lprint                                     !! A flag to print some debug statements
  !!  used in the two stream approach...every level

  !! 0.2 Output variables
  LOGICAL,INTENT(OUT) :: lconverged                                       !! did the optimization converge?
  REAL(r_std),DIMENSION(:),INTENT(OUT)          :: Collim_Alb_Coll        !! collimated total albedo for the converged case
  !! unitless, between 0 and 1
  REAL(r_std),DIMENSION(:),INTENT(OUT)          :: Collim_Tran_Coll       !! collimated total transmission
  !! unitless, between 0 and 1
  REAL(r_std),DIMENSION(:),INTENT(OUT)          :: Collim_Abs_Tot         !! collimated total absorption
  !! unitless, between 0 and 1
  REAL(r_std),DIMENSION(:),INTENT(OUT)          :: Isotrop_Alb_Coll       !! isotropic total albedo
  !! unitless, between 0 and 1
  REAL(r_std),DIMENSION(:),INTENT(OUT)          :: Isotrop_Tran_Coll      !! isotropic total transmission
  !! unitless, between 0 and 1
  REAL(r_std),DIMENSION(:),INTENT(OUT)          :: Isotrop_Abs_Tot        !! isotropic total absorption
  !! unitless, between 0 and 1
  REAL(r_std),DIMENSION(:),INTENT(OUT)          :: Collim_Tran_Uncoll     !! collimated uncollided transmission
  !! unitless, between 0 and 1
  REAL(r_std),DIMENSION(:),INTENT(OUT)          :: Isotrop_Tran_Uncoll    !! isotropic uncollided transmission
  !! unitless, between 0 and 1

  !! 0.3 Modified variables

  !! 0.4 Local variables
  ! use an extra level here...this is basically to store the transmission from the sun, which is normalized to 1.0
  REAL(r_std),DIMENSION(nlevels_tot)            :: Collim_Abs_Coll_Unscaled
  REAL(r_std),DIMENSION(nlevels_tot)            :: Collim_Alb_Coll_Unscaled
  REAL(r_std),DIMENSION(nlevels_tot)            :: Collim_Tran_Coll_Unscaled
  REAL(r_std),DIMENSION(nlevels_tot)            :: Collim_Tran_UnColl_Unscaled
  REAL(r_std),DIMENSION(nlevels_tot)            :: Isotrop_Abs_Coll_Unscaled
  REAL(r_std),DIMENSION(nlevels_tot)            :: Isotrop_Alb_Coll_Unscaled
  REAL(r_std),DIMENSION(nlevels_tot)            :: Isotrop_Tran_Coll_Unscaled
  REAL(r_std),DIMENSION(nlevels_tot)            :: Isotrop_Tran_UnColl_Unscaled
  REAL(r_std),DIMENSION(nlevels_tot)            :: Collim_Tran_Tot
  REAL(r_std),DIMENSION(nlevels_tot)            :: Isotrop_Tran_Tot

  REAL(r_std),DIMENSION(0:nlevels_tot+1)            :: &
       collim_down_cs,isotrop_down_cs,isotrop_up_cs

  INTEGER                                       :: istep,ilevel

  REAL(r_std)                                   :: leaf_reflectance
  REAL(r_std)                                   :: leaf_transmittance

  LOGICAL :: lexit
  LOGICAL :: ok
  REAL(r_std), DIMENSION(4)                     :: gammaCoeffs
  REAL(r_std), DIMENSION(4)                     :: gammaCoeffs_star
  REAL(r_std)                                   :: sun_zenith_angle_radians
  REAL(r_std), PARAMETER                        :: isotropic_cosine_constant=0.5/0.705

  !_ ================================================================================================================================
  
  ! convet the sun angle
  sun_zenith_angle_radians = acos(cosine_sun_angle)

  ! initialize some values that are used in the optimization loop
  istep=0
  lexit=.FALSE.
  lconverged=.FALSE.

  ! Initialize    
  Collim_Alb_Coll(:) = zero  
  Collim_Tran_Coll(:) = zero 
  Collim_Abs_Tot(:) = zero  
  Isotrop_Alb_Coll(:) = zero 
  Isotrop_Tran_Coll(:) = zero
  Isotrop_Abs_Tot(:) = zero  
  Collim_Tran_Uncoll(:) = zero 
  Isotrop_Tran_Uncoll(:) = zero
  
  ! calculate the albedo at our starting point
  leaf_transmittance = leaf_single_scattering_albedo_start/ & 
       ( leaf_psd_start+un)
  leaf_reflectance =  leaf_psd_start * &
       leaf_transmittance

  ! some debugging stuff
  !  DO ilevel=1,nlevels_tot
  !
  !     CALL print_debugging_albedo_info(ilevel,cosine_sun_angle,leaf_reflectance,&
  !          leaf_transmittance,&
  !          laieff_collim_pft(ilevel),laieff_isotrop_pft(ilevel), &
  !          reflectance_collim(ilevel-1),reflectance_isotrop(ilevel-1),&
  !          Collim_Alb_Tot_temp(ilevel),Collim_Tran_Tot_temp(ilevel),&
  !          Collim_Abs_Tot_temp(ilevel),&
  !          Isotrop_Alb_Tot_temp(ilevel),Isotrop_Tran_Tot_temp(ilevel),&
  !          Isotrop_Abs_Tot_temp(ilevel))
  !  ENDDO

  ! calculate the gamma coefficients used in the case of the black background 
  call gammas(leaf_reflectance, leaf_transmittance, cosine_sun_angle, gammaCoeffs)
  call gammas(leaf_reflectance, leaf_transmittance, isotropic_cosine_constant,&
       gammaCoeffs_star)
     
  !************************** step one ****** !
  ! compute all the unscaled quantities
  DO ilevel=nlevels_tot,1,-1

     Collim_Tran_UnColl_Unscaled(ilevel)=exp( - 0.5_r_std*&
          laieff_collim_pft(ilevel)/cosine_sun_angle)

     Isotrop_Tran_UnColl_Unscaled(ilevel)= &
          TBarreUncoll_exact(0.5_r_std*laieff_isotrop_pft(ilevel))


     !  /* 1) collimated source */
     ok=dhrT1(leaf_reflectance,leaf_transmittance,&
          gammaCoeffs(1),gammaCoeffs(2),gammaCoeffs(3),gammaCoeffs(4),&
          sun_zenith_angle_radians, 0.5_r_std*laieff_collim_pft(ilevel),&
          Collim_Alb_Coll_Unscaled(ilevel),Collim_Tran_Coll_Unscaled(ilevel),&
          Collim_Abs_Coll_Unscaled(ilevel))
     !  /* 2) isotropic source */
     ok=dhrT1(leaf_reflectance,leaf_transmittance,&
          gammaCoeffs_star(1),gammaCoeffs_star(2),gammaCoeffs_star(3),gammaCoeffs_star(4),&
          acos(isotropic_cosine_constant),0.5_r_std*laieff_isotrop_pft(ilevel),&
          Isotrop_Alb_Coll_Unscaled(ilevel),Isotrop_Tran_Coll_Unscaled(ilevel),&
          Isotrop_Abs_Coll_Unscaled(ilevel))
  ENDDO ! ilevel=nlevels_tot,1,-1
  
  ! Following the fate of the light at every step

  ! The downwelling array indicates the quantity of light flowing into this level
  ! from above therefore, to start our system for collimated light, we give a 
  ! unit of light coming into the top level from a collimated source.
  ! The upwelling array is light entering this level from below.
  ! cs is the current step, while ns is the next step for the iteration.
  ! Level 0 is the background. Light can enter this level from above but not below
  ! Level nlevels_tot+1 is the atomosphere. Light can enter this level from below 
  ! but not above.
  collim_down_cs(:)=zero
  collim_down_cs(nlevels_tot)=un
  isotrop_down_cs(:)=zero
  isotrop_up_cs(:)=zero

  IF(lprint)THEN
     WRITE(numout,*) "Solving fluxes for a collimated light source (i.e., the sun) in the canopy."
  ENDIF
  CALL propagate_fluxes(collim_down_cs, isotrop_down_cs, isotrop_up_cs, &
       Collim_Tran_UnColl_Unscaled, Collim_Tran_Coll_Unscaled, &
       Collim_Alb_Coll_Unscaled, Isotrop_Tran_UnColl_Unscaled, &
       Isotrop_Tran_Coll_Unscaled, Isotrop_Alb_Coll_Unscaled, &
       br_base_temp_collim, br_base_temp_isotrop, .FALSE., &
       Collim_Tran_Uncoll, Collim_Tran_Coll, Collim_Alb_Coll, lconverged, lprint)

  ! Now for the isotropic light
  collim_down_cs(:)=zero
  isotrop_down_cs(:)=zero
  isotrop_down_cs(nlevels_tot)=un
  isotrop_up_cs(:)=zero

  IF(lprint)THEN
     WRITE(numout,*) "Solving fluxes for diffuse light sources (e.g., clouds, aerosols) in the canopy."
  ENDIF
  ! The colliminated light is never used here in propagate_fluxes for the
  ! isotropic light source, but it's passed to keep things 
  ! clean between the two cases. 
  CALL propagate_fluxes(collim_down_cs, isotrop_down_cs, isotrop_up_cs, &
       Collim_Tran_UnColl_Unscaled, Collim_Tran_Coll_Unscaled, &
       Collim_Alb_Coll_Unscaled, Isotrop_Tran_UnColl_Unscaled, &
       Isotrop_Tran_Coll_Unscaled, Isotrop_Alb_Coll_Unscaled, &
       br_base_temp_collim, br_base_temp_isotrop, .TRUE., &
       Isotrop_Tran_Uncoll, Isotrop_Tran_Coll, Isotrop_Alb_Coll, lconverged, lprint)

  ! Calculate the absorption profile
  Collim_Tran_Tot(:)=Collim_Tran_Coll(:)+Collim_Tran_Uncoll(:)
  Isotrop_Tran_Tot(:)=Isotrop_Tran_Coll(:)+Isotrop_Tran_Uncoll(:)

  ! bottom level
  Collim_Abs_Tot(1)=Collim_Tran_Tot(1+1)+Collim_Tran_Tot(1)*br_base_temp_collim&
       -Collim_Tran_Tot(1)-Collim_Alb_Coll(1)
  Isotrop_Abs_Tot(1)=Isotrop_Tran_Tot(1+1)+Isotrop_Tran_Tot(1)*br_base_temp_isotrop&
       -Isotrop_Tran_Tot(1)-Isotrop_Alb_Coll(1)

  ! all middle levels
  DO ilevel=2,nlevels_tot-1
     Collim_Abs_Tot(ilevel)=Collim_Tran_Tot(ilevel+1)+Collim_Alb_Coll(ilevel-1)&
          -Collim_Tran_Tot(ilevel)-Collim_Alb_Coll(ilevel)
     Isotrop_Abs_Tot(ilevel)=Isotrop_Tran_Tot(ilevel+1)+Isotrop_Alb_Coll(ilevel-1)&
          -Isotrop_Tran_Tot(ilevel)-Isotrop_Alb_Coll(ilevel)
  ENDDO

  ! top level
  Collim_Abs_Tot(nlevels_tot)=un+Collim_Alb_Coll(nlevels_tot-1)&
       -Collim_Tran_Tot(nlevels_tot)-Collim_Alb_Coll(nlevels_tot)
  Isotrop_Abs_Tot(nlevels_tot)=un+Isotrop_Alb_Coll(nlevels_tot-1)&
       -Isotrop_Tran_Tot(nlevels_tot)-Isotrop_Alb_Coll(nlevels_tot)

END SUBROUTINE multilevel_albedo
!! ==============================================================================================================================
!! SUBROUTINE   : multilevel_matrix
!!
!>\BRIEF        : Multilevel matrix radiation transfer
!!
!! DESCRIPTION : Matrix radiation transfer where the mathematical representation
!!               of radiative transfer involves using matrices to describe the transfer
!!               of energy and properties of radiation within the medium.
!!
!! RECENT CHANGE(S): None

!!
!! MAIN OUTPUT VARIABLE(S): Upward flux, downward flux, and absorption
!!
!! REFERENCE(S) : K. Ardaneh et al. / Procedia Computer Science 255 (2025) 3–12
!!
!! FLOWCHART    :
!! 1. Initialize Inputs
!!
!! Read vegetation canopy properties (LAI, reflectance, transmittance, scattering properties).
!! Define incident solar flux, optical depth, and angular parameters.
!! Set the number of layers nl.
!!
!! 2. Compute Two-Stream Coefficients
!!
!! Compute G(mu), K(mu), mu_bar, omega, omega_beta, and omega_beta0.
!! Calculate two-stream approximation coefficients gamma1, gamma2, gamma3, gamma4.
!!
!! 3. Formulate the Matrix System
!!
!! Define the linear equations governing radiative transfer in multiple layers.
!! Construct a tridiagonal system with coefficients e_l, d_l, f_l, g_l.
!!
!! 4. Apply Boundary Conditions
!!
!! At the top: Incoming downward flux.
!! At the bottom: Upward flux determined by surface reflectance.
!!
!! 5. Solve Tridiagonal System
!!
!! Use the Thomas Algorithm for efficient computation.
!! Compute flux amplitudes A_l, B_l for all layers.
!!
!! 6. Compute Final Radiative Fluxes
!!
!! Calculate downward and upward fluxes Fdn_l and Fup_l using solved amplitudes.
!! Determine surface albedo, absorption (Fab), and transmission.
!!
!_ ================================================================================================================================
SUBROUTINE multilevel_matrix(nl, Mu, Rs, t, w, d, Fup, Fdn, Fab)
  !! 0. Variables and parameter declaration
  !! 0.1 Input variables
  INTEGER, INTENT(IN) :: nl
  REAL(KIND = r_std), INTENT(IN) :: Mu, Rs, w, d
  REAL(KIND = r_std), DIMENSION(nlevels_tot), INTENT(IN) :: t
  !! 0.2 Output variables
  REAL(KIND = r_std), DIMENSION(nlevels_tot), INTENT(OUT) :: Fup, Fdn, Fab
  !! 0.3 Modified variables
  !! 0.4 Local variables
  REAL(KIND = r_std), DIMENSION(nl * 2) :: xn, an, bn, cn, gn
  REAL(KIND = r_std), DIMENSION(2, nlevels_tot) :: data_pu
  REAL(KIND = r_std), DIMENSION(11, nl) :: data_nl
  REAL(KIND = r_std), DIMENSION(4) :: Gam
  REAL(KIND = r_std), PARAMETER :: G = 0.5_r_std
  REAL(KIND = r_std), PARAMETER :: Mubar = 1.0_r_std
  REAL(KIND = r_std) :: Bs, Lambda, Cap_Gam, Exp_trm, Cpn, Cmn
  REAL(KIND = r_std) :: Den, Exp_top, Exp_bot, Fdu, Fdb
  INTEGER :: l, le, lo
  REAL(KIND = r_std) :: Fspi, mu0, K, Beta, Beta0
  !_ ================================================================================================================================
  Fab(:) = 0.0_r_std
  K = G / Mu ! K in Eq. (12)
  Beta0 = (0.5_r_std / w) * (w + 2.0_r_std * d * Mu / 3.0_r_std) ! β0 Eq. (12)
  Beta = (0.5_r_std / w) * (w + d / 3.0_r_std) ! β in Eq. (12)

  ! Gamma coefficients (γ₁, γ₂, γ₃, and γ₄ = 1 − γ₃) according to the δ-methods
  Gam(1) = (1.0_r_std - (1.0_r_std - Beta) * w)
  Gam(2) = w * Beta
  Gam(3) = Beta0
  Gam(4) = 1.0_r_std - Gam(3)

  ! Transformation of the two-stream equations for vegetation
  mu0 = 1.0_r_std / K
  Fspi = 1.0_r_std / mu0

  data_nl(1, :) = t(nl : 1 : - 1) ! Δτₗ in Fig. 1

  Lambda = SQRT(Gam(1) ** 2 - Gam(2) ** 2) ! λ in Eq. (3)
  Cap_Gam = Gam(2) / (Gam(1) + Lambda) ! Γ in Eq. (3)

  ! Cumulative optical depths (τc)
  data_pu(1, 1) = 0.0_r_std
  DO l = 2, nlevels_tot
    data_pu(1, l) = data_pu(1, l - 1) + data_nl(1, l - 1)
  END DO

  Fdb = 0.0_r_std  ! diffuse flux at top
  data_pu(2, 1) = mu0 * Fspi ! direct flux at top

  DO l = 1, nl

    ! Eqs. (8)
    Exp_trm = EXP(- MIN(Lambda * data_nl(1, l), 30.0_r_std))

    data_nl(2, l) = 1.0_r_std + Cap_Gam * Exp_trm
    data_nl(3, l) = 1.0_r_std - Cap_Gam * Exp_trm
    data_nl(4, l) = Cap_Gam + Exp_trm
    data_nl(5, l) = Cap_Gam - Exp_trm

    ! Exponential terms in Eq. (3)
    Cpn = w * Fspi * ((Gam(1) - 1.0_r_std / mu0) * Gam(3) + Gam(4) * Gam(2))
    Cmn = w * Fspi * ((Gam(1) + 1.0_r_std / mu0) * Gam(4) + Gam(2) * Gam(3))
    Den = Lambda ** 2 - 1.0_r_std / mu0 ** 2

    Exp_top = EXP(- MIN(data_pu(1, l) / mu0, 30.0_r_std))
    Exp_bot = Exp_top * EXP(- MIN(data_nl(1, l) / mu0, 30.0_r_std))

    data_pu(2, l + 1) = mu0 * Fspi * Exp_bot

    data_nl(6, l) = Exp_top * Cpn / Den
    data_nl(7, l) = Exp_bot * Cpn / Den
    data_nl(8, l) = Exp_top * Cmn / Den
    data_nl(9, l) = Exp_bot * Cmn / Den

  END DO

  ! Boundary conditions at the top  of the structure.
  Bs = Rs * data_pu(2, nlevels_tot)
  an(1) = 0.0_r_std
  bn(1) = data_nl(2, 1)
  cn(1) = - data_nl(3, 1)
  gn(1) = Fdb - data_nl(8, 1)

  ! Matrix of coefficients in Eqs. (9)
  DO l = 1, nl - 1

    le = 2 * l
    lo = 2 * l + 1

    an(le) = &
            data_nl(3, l + 1) * data_nl(2, l) - &
            data_nl(5, l + 1) * data_nl(4, l)
    an(lo) = &
            data_nl(3, l) * data_nl(4, l) - &
            data_nl(5, l) * data_nl(2, l)

    bn(le) = &
            data_nl(3, l) * data_nl(3, l + 1) - &
            data_nl(5, l) * data_nl(5, l + 1)
    bn(lo) = &
            data_nl(2, l) * data_nl(2, l + 1) - &
            data_nl(4, l) * data_nl(4, l + 1)

    cn(le) = &
            data_nl(2, l + 1) * data_nl(5, l + 1) - &
            data_nl(3, l + 1) * data_nl(4, l + 1)
    cn(lo) = &
            data_nl(4, l) * data_nl(5, l + 1) - &
            data_nl(2, l) * data_nl(3, l + 1)

    gn(le) = & 
            data_nl(3, l + 1) * (data_nl(6, l + 1) - data_nl(7, l)) - &
            data_nl(5, l + 1) * (data_nl(8, l + 1) - data_nl(9, l))
    gn(lo) = &
            data_nl(4, l) * (data_nl(6, l + 1) - data_nl(7, l    )) + &
            data_nl(2, l) * (data_nl(9, l    ) - data_nl(8, l + 1))

  END DO

  ! Boundary conditions at the bottom of the structure.
  l = 2 * nl
  an(l) = data_nl(2, nl) - Rs * data_nl(4, nl)
  bn(l) = data_nl(3, nl) - Rs * data_nl(5, nl)
  cn(l) = 0.0_r_std
  gn(l) = Bs - data_nl(7, nl) + Rs * data_nl(9, nl)

  CALL trdiag(l, an, bn, cn, gn, xn)

  ! Matrix solutions unpacking 
  data_nl(10, :) = xn(1 : l : 2)
  data_nl(11, :) = xn(2 : l : 2)

  ! Calculation of flux 
  Fdu = &
          data_nl(10, 1) * data_nl(4, 1) - &
          data_nl(11, 1) * data_nl(5, 1) + &
          data_nl(6, 1)
  Fup(nlevels_tot) = Fdu
  Fdn(nlevels_tot) = Fdb + data_pu(2, 1)

  DO l = 1, nl
    le = nlevels_tot - l
    Fdu = &
            data_nl(10, l) * data_nl(2, l) + &
            data_nl(11, l) * data_nl(3, l) + &
            data_nl(7, l)
    Fup(le) = Fdu
    Fdb = &
            data_nl(10, l) * data_nl(4, l) + &
            data_nl(11, l) * data_nl(5, l) + &
            data_nl(9, l)
    Fdn(le) = Fdb + data_pu(2, l + 1)
  END DO

  ! Calculation of absorption
  Fab(nl : 1 : - 1) = (Fdn(nlevels_tot : 2 : - 1) - Fup(nlevels_tot : 2 : - 1)) - (Fdn(nl : 1 : - 1) - fup(nl : 1 : - 1))

END SUBROUTINE multilevel_matrix
!! ==============================================================================================================================
!! SUBROUTINE   : trdiag
!!
!>\BRIEF        : Tridiagonal matrix algorithm solver using Thomas algorithm
!!                solves Ax = g where A is a tridiagonal matrix consisting of vectors a, b, c
!! DESCRIPTION : The Thomas algorithm, also known as the tridiagonal matrix algorithm,
!!               is an efficient method for solving systems of linear equations with
!!               a tridiagonal coefficient matrix.
!!               n - number of equations
!!               a - sub-diagonal (means it is the diagonal below the main diagonal)
!!               b - the main diagonal
!!               c - sup-diagonal (means it is the diagonal above the main diagonal)
!!               g - right part
!!               x - the answer
!!
!! RECENT CHANGE(S): None

!!
!! MAIN OUTPUT VARIABLE(S): Matrix of solutions x
!!
!! REFERENCE(S) : Joe D. Hoffman, 2001, Numerical Methods for Engineers and Scientists, Second Edition
!!
!! FLOWCHART    : None
!_ ================================================================================================================================
SUBROUTINE trdiag(n, a, b, c, g, x)
  !! 0. Variables and parameter declaration
  !! 0.1 Input variables
  INTEGER, INTENT(IN) :: n
  REAL(KIND = r_std), DIMENSION(n), INTENT(IN) :: a, c
  !! 0.2 Output variables
  REAL(KIND = r_std), DIMENSION(n), INTENT(OUT) :: x
  !! 0.3 Modified variables
  REAL(KIND = r_std), DIMENSION(n), INTENT(INOUT) :: b, g
  !! 0.4 Local variables
  INTEGER :: i
  REAL(KIND = r_std) :: t
  !_ ================================================================================================================================
  ! Forward elimination phase
  DO i = 2, n
    t = a(i) / b(i - 1)
    b(i) = b(i) - c(i - 1) * t
    g(i) = g(i) - g(i - 1) * t
  END DO
  ! Backward substitution phase
  x(n) = g(n) / b(n)
  DO i = n - 1, 1, - 1
    x(i) = (g(i) - c(i) * x(i + 1)) / b(i)
  END DO
END SUBROUTINE trdiag
!! ==============================================================================================================================
!! SUBROUTINE   : print_debugging_albedo_info
!!
!>\BRIEF         Prints out some albedo information in a nice format.  
!!               Should only be used for debugging, never for production runs.
!!
!! DESCRIPTION : 
!!
!! RECENT CHANGE(S): None\n
!!
!! MAIN OUTPUT VARIABLE(S): None.
!!
!! REFERENCE(S) :  
!!
!! FLOWCHART    : None
!_ ================================================================================================================================
SUBROUTINE print_debugging_albedo_info(ilevel, cosine_sun_angle, &
     leaf_reflectance, leaf_transmittance, laieff_collim_temp, laieff_isotrop_temp, &
     br_base_temp_collim, br_base_temp_isotrop, Collim_Alb_Tot, &
     Collim_Tran_Tot, Collim_Abs_Tot, Isotrop_Alb_Tot, Isotrop_Tran_Tot, Isotrop_Abs_Tot)

  !! 0. Variables and parameter declaration
  !! 0.1 Input variables
  INTEGER,INTENT(IN)          :: ilevel
  REAL(r_std), INTENT(IN)     :: cosine_sun_angle            !! cosine of the solar zenith angle, between 0 and 1
                                                             !! @tex $()$ @endtex 
   REAL(r_std), INTENT(IN)    :: laieff_collim_temp          !! cosine of the solar zenith angle, between 0 and 1
   REAL(r_std), INTENT(IN)    :: laieff_isotrop_temp         !! cosine of the solar zenith angle, between 0 and 1
   REAL(r_std), INTENT(IN)    :: leaf_reflectance            !! cosine of the solar zenith angle, between 0 and 1
   REAL(r_std), INTENT(IN)    :: leaf_transmittance          !! cosine of the solar zenith angle, between 0 and 1
   REAL(r_std), INTENT(IN)    :: br_base_temp_collim         !! cosine of the solar zenith angle, between 0 and 1
   REAL(r_std), INTENT(IN)    :: br_base_temp_isotrop        !! cosine of the solar zenith angle, between 0 and 1
   REAL(r_std),INTENT(IN)     :: Collim_Alb_Tot
   REAL(r_std),INTENT(IN)     :: Collim_Tran_Tot
   REAL(r_std),INTENT(IN)     :: Collim_Abs_Tot
   REAL(r_std),INTENT(IN)     :: Isotrop_Alb_Tot
   REAL(r_std),INTENT(IN)     :: Isotrop_Tran_Tot
   REAL(r_std),INTENT(IN)     :: Isotrop_Abs_Tot

  !! 0.2 Output variables

  !! 0.3 Modified variables

  !! 0.4 Local variables

  !_ ================================================================================================================================
   

     WRITE (numout,'(8(11A))') '   Level   ','   Angle   ','  laieff_c  ','  laieff_i  ',&
          '   rleaf   ','   tleaf   ','   rbgd_c   ','   rbgd_i   '

     WRITE(numout,FMT='(I6,5X,7(F11.6))') &
          ilevel,180/pi*ACOS(cosine_sun_angle),&
          laieff_collim_temp,laieff_isotrop_temp,&
          leaf_reflectance,leaf_transmittance,br_base_temp_collim,&
          br_base_temp_isotrop


     WRITE (numout,'(10X,6(11A))') ' Rtot(sun) ',' Ttot(sun) ',&
          ' Atot(sun) ',' Rtot(iso) ',' Ttot(iso) ',' Atot(iso) '
     WRITE(numout,FMT='(10X,6(F11.6))') &
          Collim_Alb_Tot,Collim_Tran_Tot,Collim_Abs_Tot,&
          Isotrop_Alb_Tot,Isotrop_Tran_Tot,Isotrop_Abs_Tot


END SUBROUTINE print_debugging_albedo_info


!! ==============================================================================================================================
!! SUBROUTINE   : propagate_fluxes
!!
!>\BRIEF         Propogates the radiation fluxes through each level of the
!!               canopy.
!!
!! DESCRIPTION : This is an algorithm to follow radition fluxes through the
!!     canopy.  At each level, the path probabilities are determined by the
!!     raditional transfer scheme of Pinty et al (2006).  Notice that
!!     the fluxes given by this routine are all cumulative fluxes, not per
!!     level.
!!
!! RECENT CHANGE(S): None\n
!!
!! MAIN OUTPUT VARIABLE(S): ::Tran_Uncoll_Tot, ::Tran_Coll_Tot, ::Alb_Coll_Tot
!!
!! REFERENCE(S) :  
!!
!! FLOWCHART    : None
!_ ================================================================================================================================
SUBROUTINE propagate_fluxes(collim_down_cs, isotrop_down_cs, isotrop_up_cs, &
     Collim_Tran_UnColl_Unscaled, Collim_Tran_Coll_Unscaled, &
     Collim_Alb_Coll_Unscaled, Isotrop_Tran_UnColl_Unscaled, &
     Isotrop_Tran_Coll_Unscaled, Isotrop_Alb_Coll_Unscaled, br_base_temp_collim, &
     br_base_temp_isotrop, lisotrop, Tran_Uncoll_Tot, Tran_Coll_Tot, &
     Alb_Coll_Tot, lconverged, lprint)

  !! 0. Variables and parameter declaration
  !! 0.1 Input variables
  REAL(r_std), DIMENSION(nlevels_tot),INTENT(IN)           :: Collim_Tran_UnColl_Unscaled
  REAL(r_std), DIMENSION(nlevels_tot),INTENT(IN)           :: Collim_Tran_Coll_Unscaled
  REAL(r_std), DIMENSION(nlevels_tot),INTENT(IN)           :: Collim_Alb_Coll_Unscaled
  REAL(r_std), DIMENSION(nlevels_tot),INTENT(IN)           :: Isotrop_Tran_UnColl_Unscaled
  REAL(r_std), DIMENSION(nlevels_tot),INTENT(IN)           :: Isotrop_Tran_Coll_Unscaled
  REAL(r_std), DIMENSION(nlevels_tot),INTENT(IN)           :: Isotrop_Alb_Coll_Unscaled
  REAL(r_std), INTENT(IN)                                  :: br_base_temp_collim             !! cosine of the solar zenith angle, between 0 and 1
  REAL(r_std), INTENT(IN)                                  :: br_base_temp_isotrop            !! cosine of the solar zenith angle, between 0 and 1
  LOGICAL, INTENT(IN)                                      :: lisotrop                        !! are we dealing with an isotropic source?  only needed for correct partitioning
                                                                                              !! of the collided and uncollided transmitted light...the total is not affected
  LOGICAL, INTENT(IN)                                      :: lprint                          !! a flag to print

  !! 0.2 Output variables
  REAL(r_std), DIMENSION(nlevels_tot),INTENT(OUT)          :: Tran_Uncoll_Tot
  REAL(r_std), DIMENSION(nlevels_tot),INTENT(OUT)          :: Tran_Coll_Tot
  REAL(r_std), DIMENSION(nlevels_tot),INTENT(OUT)          :: Alb_Coll_Tot
  LOGICAL,INTENT(OUT)                                      :: lconverged

  !! 0.3 Modified variables
  REAL(r_std), DIMENSION(0:nlevels_tot+1),INTENT(INOUT)    :: collim_down_cs
  REAL(r_std), DIMENSION(0:nlevels_tot+1),INTENT(INOUT)    :: isotrop_down_cs
  REAL(r_std), DIMENSION(0:nlevels_tot+1),INTENT(INOUT)    :: isotrop_up_cs

  !! 0.4 Local variables
  INTEGER :: istep,ilevel
  REAL(r_std),DIMENSION(0:nlevels_tot+1)            :: collim_down_ns,isotrop_down_ns,&
       isotrop_up_ns,Tran_Uncoll_Tot_temp
  REAL(r_std), DIMENSION(nlevels_tot)          ::  Tran_Coll_Tot_temp, Tran_Tot
  
  !_ ================================================================================================================================

  Tran_Uncoll_Tot(:)=zero
  Tran_Tot(:)=zero
  Alb_Coll_Tot(:)=zero

  Tran_Uncoll_Tot_temp(:)=un

  istep=0

  DO

     istep=istep+1

     lconverged=.TRUE.

     ! Zero out the counters for the next step
     collim_down_ns(:)=zero
     isotrop_down_ns(:)=zero
     isotrop_up_ns(:)=zero


     ! Now we need to loop over all levels and see what light is entering the 
     ! level, and how it will propagate in the next step
     DO ilevel=1,nlevels_tot

        ! For collimated downwelling light into the level, it can be scattered 
        ! up, down, or pass through uncollided
        IF(collim_down_cs(ilevel) .GT. zero)THEN
           collim_down_ns(ilevel-1)=collim_down_ns(ilevel-1)+&
                collim_down_cs(ilevel)*Collim_Tran_UnColl_Unscaled(ilevel)


           ! This statement checks to see if this light has been previously
           ! scattered or not.  This term is only present in level nlevels_tot
           ! for the first step, level nlevels_tot-1 for the second step, 
           ! level nlevels_tot-2 for the third step, etc., and it only happens
           ! in the case of a collimated light source.
           IF(ilevel == nlevels_tot-istep+1 .AND. .NOT. lisotrop) THEN
              Tran_Uncoll_Tot_temp(ilevel)=Tran_Uncoll_Tot_temp(ilevel+1)*&
                   Collim_Tran_UnColl_Unscaled(ilevel)
           ENDIF

           isotrop_down_ns(ilevel-1)=isotrop_down_ns(ilevel-1)+&
                collim_down_cs(ilevel)*Collim_Tran_Coll_Unscaled(ilevel)
           isotrop_up_ns(ilevel+1)=isotrop_up_ns(ilevel+1)+&
                collim_down_cs(ilevel)*Collim_Alb_Coll_Unscaled(ilevel)

        ENDIF ! collim_down_cs(ilevel) .GT. zero

        ! For isotropic downwelling light, it can also be scattered up, down,
        ! or pass through uncollided
        IF(isotrop_down_cs(ilevel) .GT. zero)THEN
           isotrop_down_ns(ilevel-1)=isotrop_down_ns(ilevel-1)+&
                isotrop_down_cs(ilevel)*Isotrop_Tran_UnColl_Unscaled(ilevel)

           ! This is the same check as above, but this time the light has
           ! an isotropic source and not a collimated source
           IF(ilevel == nlevels_tot-istep+1 .AND. lisotrop) THEN
              Tran_Uncoll_Tot_temp(ilevel)=Tran_Uncoll_Tot_temp(ilevel+1)&
                   *Isotrop_Tran_UnColl_Unscaled(ilevel)
           ENDIF

           isotrop_down_ns(ilevel-1)=isotrop_down_ns(ilevel-1)+&
                isotrop_down_cs(ilevel)*Isotrop_Tran_Coll_Unscaled(ilevel)
           isotrop_up_ns(ilevel+1)=isotrop_up_ns(ilevel+1)+&
                isotrop_down_cs(ilevel)*Isotrop_Alb_Coll_Unscaled(ilevel)
        ENDIF ! isotrop_down_cs(ilevel) .GT. zero

        ! Isotropic upwelling light can pass through upwards collided or 
        ! uncollided with vegetation in this level, or it can be reflected downwards
        IF(isotrop_up_cs(ilevel) .GT. zero)THEN

           isotrop_up_ns(ilevel+1)=isotrop_up_ns(ilevel+1)+isotrop_up_cs(ilevel)&
                *Isotrop_Tran_UnColl_Unscaled(ilevel)
           isotrop_up_ns(ilevel+1)=isotrop_up_ns(ilevel+1)+isotrop_up_cs(ilevel)&
                *Isotrop_Tran_Coll_Unscaled(ilevel)
           isotrop_down_ns(ilevel-1)=isotrop_down_ns(ilevel-1)+&
                isotrop_up_cs(ilevel)*Isotrop_Alb_Coll_Unscaled(ilevel)
        ENDIF

        ! there can be no collimated upwards light, since upwards light must 
        ! always have reflected off something


     ENDDO ! ilevel=1,nlevels_tot

     ! The background is a bit special. there is no transmission, but there is 
     ! reflection, which leads to a source term to the bottom level from below.
     isotrop_up_ns(1)=isotrop_up_ns(1)+collim_down_cs(0)*br_base_temp_collim
     isotrop_up_ns(1)=isotrop_up_ns(1)+isotrop_down_cs(0)*br_base_temp_isotrop


     ! Now we add the light generated here to our cumulative counters to track
     ! the total amount that leaves the canopy, either through being absorbed
     ! by the background or being reflected back into the atmosphere.
     ! We keep track of the uncoll above, so here we track the total light that
     ! is transmitted to the soil, and then at the end we taken the difference 
     ! to get the collided radation.
     Tran_Tot(1:nlevels_tot)=Tran_Tot(1:nlevels_tot)+&
          isotrop_down_ns(0:nlevels_tot-1)+collim_down_ns(0:nlevels_tot-1)
     Alb_Coll_Tot(1:nlevels_tot)=Alb_Coll_Tot(1:nlevels_tot)+isotrop_up_ns(2:nlevels_tot+1)

     ! now we update the light we are currently tracking
     collim_down_cs(:)=collim_down_ns(:)
     isotrop_down_cs(:)=isotrop_down_ns(:)
     isotrop_up_cs(:)=isotrop_up_ns(:)


     ! check for convegence...if all the values of light are currently below a 
     ! threshold, we're probably good assuming the threshold is low enough.
     IF(MAXVAL(collim_down_cs,1) .GT. alb_threshold) lconverged=.FALSE.
     IF(MAXVAL(isotrop_down_cs,1) .GT. alb_threshold) lconverged=.FALSE.
     IF(MAXVAL(isotrop_up_cs,1) .GT. alb_threshold) lconverged=.FALSE.

     IF(lconverged)THEN
        IF(lprint) WRITE(numout,*) 'Converged after this many steps: ',istep
        EXIT
     ENDIF

     ! This number here could also be externalized.
     IF(istep .GE. 1000)THEN
        WRITE(numout,*) '*********************************************************'
        WRITE(numout,'(A,I6,F14.10)') 'Albedo not converging!',istep,alb_threshold
        WRITE(numout,'(A)') '                  collim_down_cs ' // &
             'isotrop_down_cs  isotrop_up_cs'
        DO ilevel=0,nlevels_tot+1
           WRITE(numout,'(I4,3F14.10)') ilevel, &
                collim_down_cs(ilevel),isotrop_down_cs(ilevel),isotrop_up_cs(ilevel)

        ENDDO
        WRITE(numout,*) 'You should increase either the number' // &
             'of steps or the alb_threshold.'
        WRITE(numout,*) '*********************************************************'
        EXIT
     ENDIF ! istep .GE. 1000
  ENDDO ! convergence loop


  ! now separate the collided from the uncollided light.  Notice that
  ! this is not really needed for any purposes other than debugging, as the
  ! important quantity is the total amount of light striking the ground.
  Tran_UnColl_Tot(1:nlevels_tot)=Tran_UnColl_Tot_temp(1:nlevels_tot)
  ! added lines since the ticket 850
  Tran_Coll_Tot_temp(:)=Tran_Tot(:)-Tran_UnColl_Tot(:)
  WHERE((Tran_Coll_Tot_temp<0) .AND. (Tran_UnColl_Tot .eq. 1))
        Tran_UnColl_Tot=0
  ENDWHERE 
  Tran_Coll_Tot(1:nlevels_tot)=Tran_Tot(1:nlevels_tot)-Tran_UnColl_Tot(1:nlevels_tot)
  
  ! Some debugging information
  IF(lprint)THEN
     WRITE(numout,'(7X,3(A15,3X))') 'Tran_Uncoll_Tot','  Tran_Coll_Tot','   Alb_Coll_Tot'
     DO ilevel=nlevels_tot,1,-1
        WRITE(numout,'(I4,3X,3(F15.6,3X))') ilevel, &
             Tran_Uncoll_Tot(ilevel),Tran_Coll_Tot(ilevel),Alb_Coll_Tot(ilevel)
     ENDDO
  ENDIF

END SUBROUTINE propagate_fluxes


!! ==============================================================================================================================
!! SUBROUTINE   : read_background_albedo
!!
!>\BRIEF        This subroutine reads the background albedo
!!
!! DESCRIPTION  This subroutine reads the background albedo map in 0.5 x 0.5 deg resolution 
!! derived from JRCTIP product. These values are then interpolated to the resolution of the
!! simulation. For deserts and fallow croplands, the background albedo will
!! be similar to the bare soil albedo. For all other vegetated PFTs, the background albedo 
!! will be determined by the understory and/or the litter layer.
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): bckgrnd_alb for visible and near-infrared range 
!!
!! REFERENCES   : None
!!
!! FLOWCHART    : None
!! \n
!_ ==============================================================================================================================

  SUBROUTINE read_background_albedo(nbpt, lalo, neighbours, resolution, contfrac)

    USE interpweight

    IMPLICIT NONE

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                    :: nbpt                  !! Number of points for which the data needs to be 
                                                                           !! interpolated (unitless)             
    REAL(r_std), INTENT(in)                       :: lalo(nbpt,2)          !! Vector of latitude and longitudes (degree)        
    INTEGER(i_std), INTENT(in)                    :: neighbours(nbpt,NbNeighb)!! Vector of neighbours for each grid point 
                                                                           !! (1=N, 2=E, 3=S, 4=W)  
    REAL(r_std), INTENT(in)                       :: resolution(nbpt,2)    !! The size of each grid cell in X and Y (km)
    REAL(r_std), INTENT(in)                       :: contfrac(nbpt)        !! Fraction of land in each grid cell (unitless)   

    !! 0.4 Local variables

    CHARACTER(LEN=80)                             :: filename              !! Filename of background albedo
    REAL(r_std), DIMENSION(nbpt)                  :: aalb_bg               !! Availability of the interpolation
    INTEGER                                       :: ALLOC_ERR             !! Help varialbe to count allocation error
    REAL(r_std)                                   :: vmin, vmax            !! min/max values to use for the 
                                                                           !!   renormalization
    CHARACTER(LEN=80)                             :: variablename          !! Variable to interpolate
    CHARACTER(LEN=80)                             :: lonname, latname      !! lon, lat names in input file
    CHARACTER(LEN=50)                             :: fractype              !! method of calculation of fraction
                                                                           !!   'XYKindTime': Input values are kinds 
                                                                           !!     of something with a temporal 
                                                                           !!     evolution on the dx*dy matrix'
    LOGICAL                                       :: nonegative            !! whether negative values should be removed
    CHARACTER(LEN=50)                             :: maskingtype           !! Type of masking
                                                                           !!   'nomask': no-mask is applied
                                                                           !!   'mbelow': take values below maskvals(1)
                                                                           !!   'mabove': take values above maskvals(1)
                                                                           !!   'msumrange': take values within 2 ranges;
                                                                           !!      maskvals(2) <= SUM(vals(k)) <= maskvals(1)
                                                                           !!      maskvals(1) < SUM(vals(k)) <= maskvals(3)
                                                                           !!        (normalized by maskedvals(3))
                                                                           !!   'var': mask values are taken from a 
                                                                           !!     variable inside the file (>0)
    REAL(r_std), DIMENSION(3)                     :: maskvals              !! values to use to mask (according to 
                                                                           !!   `maskingtype') 
    CHARACTER(LEN=250)                            :: namemaskvar           !! name of the variable to use to mask 
    REAL(r_std)                                   :: albbg_norefinf        !! No value
    REAL(r_std), ALLOCATABLE, DIMENSION(:)        :: albbg_default         !! Default value

!_ ================================================================================================================================

  !! 1. Open file and allocate memory

  ! Open file with background albedo

  !Config Key   = ALB_BG_FILE
  !Config Desc  = Name of file from which the background albedo is read 
  !Config Def   = alb_bg.nc
  !Config If    = ALB_BG_MODIS
  !Config Help  = The name of the file to be opened to read background albedo 
  !Config Units = [FILE]
  !
  filename = 'alb_bg.nc'
  CALL getin_p('ALB_BG_FILE',filename)

 
  IF (xios_interpolation) THEN
     IF (printlev_loc >= 2) WRITE(numout,*)'Now start reading background albedo with XIOS'
     ! Read and interpolation background albedo using XIOS
     CALL xios_orchidee_recv_field('bg_alb_vis_interp',bckgrnd_alb(:,ivis))
     CALL xios_orchidee_recv_field('bg_alb_nir_interp',bckgrnd_alb(:,inir))
     
     aalb_bg(:)=1
     
  ELSE
     IF (printlev_loc >= 2) WRITE(numout,*)'Now start reading background albedo with IOIPSL'
     ALLOCATE(albbg_default(2), STAT=ALLOC_ERR)
     IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'read_background_albedo','Pb in allocation for albbg_default','','')
     
     ! For this case there are not types/categories. We have 'only' a continuous field
     ! Assigning values to vmin, vmax

     vmin = 0.
     vmax = 9999.
  
     !! Variables for interpweight
     ! Type of calculation of cell fractions (not used here)
     fractype = 'default'
     ! Name of the longitude and latitude in the input file
     lonname = 'longitude'
     latname = 'latitude'
     ! Default value when no value is get from input file
     albbg_default(ivis) = 0.129
     albbg_default(inir) = 0.247
     ! Reference value when no value is get from input file (not used here)
     albbg_norefinf = undef_sechiba
     ! Should negative values be set to zero from input file?
     nonegative = .FALSE.
     ! Type of mask to apply to the input data (see header for more details)
     maskingtype = 'var'
     ! Values to use for the masking (here not used)
     maskvals = (/ undef_sechiba, undef_sechiba, undef_sechiba /)
     ! Name of the variable with the values for the mask in the input file (only if maskkingtype='var')
     namemaskvar = 'mask'
     
     ! There is a variable for each chanel 'infrared' and 'visible'
     ! Interpolate variable bg_alb_vis
     variablename = 'bg_alb_vis'
     IF (printlev_loc >= 2) WRITE(numout,*) "read_background_albedo: Start interpolate " &
          // TRIM(filename) // " for variable " // TRIM(variablename)
     CALL interpweight_2Dcont(nbpt, 0, 0, lalo, resolution, neighbours,                                  &
          contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype,          &
          maskvals, namemaskvar, -1, fractype, albbg_default(ivis), albbg_norefinf,                         &
          bckgrnd_alb(:,ivis), aalb_bg)
     IF (printlev_loc >= 5) WRITE(numout,*)"  read_background_albedo after InterpWeight2Dcont for '" //   &
          TRIM(variablename) // "'"
     
     ! Interpolate variable bg_alb_nir in the same file
     variablename = 'bg_alb_nir'
     IF (printlev_loc >= 2) WRITE(numout,*) "read_background_albedo: Start interpolate " &
          // TRIM(filename) // " for variable " // TRIM(variablename)
     CALL interpweight_2Dcont(nbpt, 0, 0, lalo, resolution, neighbours,                                  &
          contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype,          &
          maskvals, namemaskvar, -1, fractype, albbg_default(inir), albbg_norefinf,                         &
          bckgrnd_alb(:,inir), aalb_bg)
     IF (printlev_loc >= 5) WRITE(numout,*)"  read_background_albedo after InterpWeight2Dcont for '" //   &
          TRIM(variablename) // "'"
     
     IF (ALLOCATED(albbg_default)) DEALLOCATE(albbg_default)
     
  ENDIF
  
  ! Write diagnostics
  CALL xios_orchidee_send_field("interp_diag_alb_vis",bckgrnd_alb(:,ivis))
  CALL xios_orchidee_send_field("interp_diag_alb_nir",bckgrnd_alb(:,inir)) 
  CALL xios_orchidee_send_field("interp_avail_aalb_bg",aalb_bg)

  IF (printlev_loc >= 3) WRITE(numout,*)'  read_background_albedo ended'


  END SUBROUTINE read_background_albedo


!! ==============================================================================================================================
!! SUBROUTINE   : albedo_surface_soilalb
!!
!>\BRIEF        This subroutine calculates the albedo of soil (without snow).
!!
!! DESCRIPTION  This subroutine reads the soil colour maps in 1 x 1 deg resolution 
!! from the Henderson-Sellers & Wilson database. These values are interpolated to 
!! the model's resolution and transformed into 
!! dry and wet albedos.\n
!!
!! If the soil albedo is calculated without the dependence of soil moisture, the
!! soil colour values are transformed into mean soil albedo values.\n 
!!
!! The calculations follow the assumption that the grid of the data is regular and 
!! it covers the globe. The calculation for the model grid are based on the borders 
!! of the grid of the resolution.
!!
!! RECENT CHANGE(S): None
!!
!! CALCULATED MODULE VARIABLE(S): soilalb_dry for visible and near-infrared range,
!!                                soilalb_wet for visible and near-infrared range, 
!!                                soilalb_moy for visible and near-infrared range 
!!
!! REFERENCE(S) : 
!! -Wilson, M.F., and A. Henderson-Sellers, 1985: A global archive of land cover and
!!  soils data for use in general circulation climate models. J. Clim., 5, 119-143.
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  
  SUBROUTINE albedo_surface_soilalb(nbpt, lalo, neighbours, resolution, contfrac)
    
    USE interpweight
    
    IMPLICIT NONE
    
    
    !! 0. Variable and parameter declaration
    
    !! 0.1 Input variables
    
    INTEGER(i_std), INTENT(in)                    :: nbpt                  !! Number of points for which the data needs to be 
                                                                           !! interpolated (unitless)             
    REAL(r_std), INTENT(in)                       :: lalo(nbpt,2)          !! Vector of latitude and longitudes (degree)        
    INTEGER(i_std), INTENT(in)                    :: neighbours(nbpt,NbNeighb)!! Vector of neighbours for each grid point 
                                                                           !! (1=N, 2=E, 3=S, 4=W)  
    REAL(r_std), INTENT(in)                       :: resolution(nbpt,2)    !! The size of each grid cell in X and Y (km)
    REAL(r_std), INTENT(in)                       :: contfrac(nbpt)        !! Fraction of land in each grid cell (unitless)   

    !! 0.4 Local variables

    CHARACTER(LEN=80)                             :: filename              !! Filename of soil colour map
    INTEGER(i_std)                                :: i, ib, ip, nbexp      !! Indices
    INTEGER                                       :: ALLOC_ERR             !! Help varialbe to count allocation error
    REAL(r_std), DIMENSION(nbpt)                  :: asoilcol              !! Availability of the soilcol interpolation
    REAL(r_std), DIMENSION(:), ALLOCATABLE        :: variabletypevals      !! Values for all the types of the variable
                                                                           !!   (variabletypevals(1) = -un, not used)
    REAL(r_std), DIMENSION(:,:), ALLOCATABLE      :: soilcolrefrac         !! soilcol fractions re-dimensioned
    REAL(r_std)                                   :: vmin, vmax            !! min/max values to use for the 
                                                                           !!   renormalization
    CHARACTER(LEN=80)                             :: variablename          !! Variable to interpolate
    CHARACTER(LEN=80)                             :: lonname, latname      !! lon, lat names in input file
    CHARACTER(LEN=50)                             :: fractype              !! method of calculation of fraction
                                                                           !!   'XYKindTime': Input values are kinds 
                                                                           !!     of something with a temporal 
                                                                           !!     evolution on the dx*dy matrix'
    LOGICAL                                       :: nonegative            !! whether negative values should be removed
    CHARACTER(LEN=50)                             :: maskingtype           !! Type of masking
                                                                           !!   'nomask': no-mask is applied
                                                                           !!   'mbelow': take values below maskvals(1)
                                                                           !!   'mabove': take values above maskvals(1)
                                                                           !!   'msumrange': take values within 2 ranges;
                                                                           !!      maskvals(2) <= SUM(vals(k)) <= maskvals(1)
                                                                           !!      maskvals(1) < SUM(vals(k)) <= maskvals(3)
                                                                           !!        (normalized by maskvals(3))
                                                                           !!   'var': mask values are taken from a 
                                                                           !!     variable inside the file (>0)
    REAL(r_std), DIMENSION(3)                     :: maskvals              !! values to use to mask (according to 
                                                                           !!   `maskingtype') 
    CHARACTER(LEN=250)                            :: namemaskvar           !! name of the variable to use to mask 
    CHARACTER(LEN=250)                            :: msg
    INTEGER                                       :: fopt
    INTEGER(i_std), DIMENSION(:), ALLOCATABLE     :: vecpos
    INTEGER(i_std), DIMENSION(:), ALLOCATABLE     :: solt

!_ ================================================================================================================================


    ! Open file with soil colours 

    !Config Key   = SOILALB_FILE
    !Config Desc  = Name of file from which the bare soil albedo
    !Config Def   = soils_param.nc
    !Config If    = NOT(IMPOSE_AZE)
    !Config Help  = The name of the file to be opened to read the soil types from 
    !Config         which we derive then the bare soil albedos. This file is 1x1 
    !Config         deg and based on the soil colors defined by Wilson and Henderson-Seller.
    !Config Units = [FILE]
    filename = 'soils_param.nc'
    CALL getin_p('SOILALB_FILE',filename)
    
    ! Read and interpolate variable soilcolor
    variablename = 'soilcolor'
    
    IF (xios_interpolation) THEN
       IF (printlev_loc >= 1) WRITE(numout,*) "albedo_surface_soilalb: Use XIOS to read and interpolate " &
            // TRIM(filename) // " for variable " // TRIM(variablename)
       
       ALLOCATE(soilcolrefrac(nbpt, classnb), STAT=ALLOC_ERR) 
       IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'albedo_surface_soilalb','Problem in allocation of variable soilcolrefrac','','')
       
       CALL xios_orchidee_recv_field('soilcolor1',soilcolrefrac(:,1))
       CALL xios_orchidee_recv_field('soilcolor2',soilcolrefrac(:,2))
       CALL xios_orchidee_recv_field('soilcolor3',soilcolrefrac(:,3))
       CALL xios_orchidee_recv_field('soilcolor4',soilcolrefrac(:,4))
       CALL xios_orchidee_recv_field('soilcolor5',soilcolrefrac(:,5))
       CALL xios_orchidee_recv_field('soilcolor6',soilcolrefrac(:,6))
       CALL xios_orchidee_recv_field('soilcolor7',soilcolrefrac(:,7))
       CALL xios_orchidee_recv_field('soilcolor8',soilcolrefrac(:,8))
       CALL xios_orchidee_recv_field('soilcolor9',soilcolrefrac(:,9))
       
       nbexp = 0
       soilalb_dry(:,:) = zero
       soilalb_wet(:,:) = zero
       soilalb_moy(:,:) = zero
       DO ib=1,nbpt ! Loop over domain size
          asoilcol(ib) = SUM(soilcolrefrac(ib,:))
          
          IF (asoilcol(ib) .LT. min_sechiba) THEN
             IF (printlev>=1) WRITE(numout,*)'  albedo_surface_soilalb: for point=', ib, ' no soil class read from file! still frac read=',asoilcol(ib)
             
             ! Initialize with mean value if no points were interpolated or if no data was found
             nbexp = nbexp + 1
             soilalb_dry(ib,ivis) = (SUM(vis_dry)/classnb + SUM(vis_wet)/classnb)/deux
             soilalb_dry(ib,inir) = (SUM(nir_dry)/classnb + SUM(nir_wet)/classnb)/deux
             soilalb_wet(ib,ivis) = (SUM(vis_dry)/classnb + SUM(vis_wet)/classnb)/deux
             soilalb_wet(ib,inir) = (SUM(nir_dry)/classnb + SUM(nir_wet)/classnb)/deux
             soilalb_moy(ib,ivis) = SUM(albsoil_vis)/classnb
             soilalb_moy(ib,inir) = SUM(albsoil_nir)/classnb
          ELSE
             ! Take into account all different fractons of different soilcolors read from file. Divide by asoilcol corresponding the total fracion in the grid cell which has information. 
             DO ip=1,classnb
                soilalb_dry(ib,ivis) = soilalb_dry(ib,ivis) + vis_dry(ip)*soilcolrefrac(ib,ip)/asoilcol(ib)
                soilalb_dry(ib,inir) = soilalb_dry(ib,inir) + nir_dry(ip)*soilcolrefrac(ib,ip)/asoilcol(ib)
                soilalb_wet(ib,ivis) = soilalb_wet(ib,ivis) + vis_wet(ip)*soilcolrefrac(ib,ip)/asoilcol(ib)
                soilalb_wet(ib,inir) = soilalb_wet(ib,inir) + nir_wet(ip)*soilcolrefrac(ib,ip)/asoilcol(ib)
                soilalb_moy(ib,ivis) = soilalb_moy(ib,ivis) + albsoil_vis(ip)*soilcolrefrac(ib,ip)/asoilcol(ib)
                soilalb_moy(ib,inir) = soilalb_moy(ib,inir) + albsoil_nir(ip)*soilcolrefrac(ib,ip)/asoilcol(ib)
             END DO
          END IF
          
       END DO
       
       IF ( (nbexp .GT. 0) .AND. printlev>=1) THEN
          WRITE(numout,*) 'albedo_surface_soilalb _______'
          WRITE(numout,*) 'albedo_surface_soilalb: The interpolation of the bare soil albedo had ', nbexp
          WRITE(numout,*) 'albedo_surface_soilalb: points without data. This are either coastal points or'
          WRITE(numout,*) 'albedo_surface_soilalb: ice covered land.'
          WRITE(numout,*) 'albedo_surface_soilalb: The problem was solved by using the average of all soils'
          WRITE(numout,*) 'albedo_surface_soilalb: in dry and wet conditions'
          WRITE(numout,*) 'albedo_surface_soilalb: Use the diagnostic output field interp_avail_asoilcol to see location of these points'
       ENDIF
       
       DEALLOCATE(soilcolrefrac)
       
    ELSE
       
       ! Use standard method with IOIPSL and aggregate to read and interpolate (without XIOS)
       ALLOCATE(soilcolrefrac(nbpt, classnb), STAT=ALLOC_ERR) 
       IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'albedo_surface_soilalb','Problem in allocation of variable soilcolrefrac','','')
       ALLOCATE(vecpos(classnb), STAT=ALLOC_ERR) 
       IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'albedo_surface_soilalb','Problem in allocation of variable vecpos','','')
       ALLOCATE(solt(classnb), STAT=ALLOC_ERR)
       IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'albedo_surface_soilalb','Problem in allocation of variable solt','','')
       
       ! Assigning values to vmin, vmax
       vmin = 1.0
       vmax = classnb
       
       ALLOCATE(variabletypevals(classnb),STAT=ALLOC_ERR)
       IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variabletypevals','','')
       variabletypevals = -un
       
       !! Variables for interpweight
       ! Type of calculation of cell fractions
       fractype = 'default'
       ! Name of the longitude and latitude in the input file
       lonname = 'nav_lon'
       latname = 'nav_lat'
       ! Should negative values be set to zero from input file?
       nonegative = .FALSE.
       ! Type of mask to apply to the input data (see header for more details)
       maskingtype = 'mabove'
       ! Values to use for the masking
       maskvals = (/ min_sechiba, undef_sechiba, undef_sechiba /)
       ! Name of the variable with the values for the mask in the input file (only if maskkingtype='var') (here not used)
       namemaskvar = ''
       
       IF (printlev_loc >= 2) WRITE(numout,*) "albedo_surface_soilalb: Start interpolate " &
            // TRIM(filename) // " for variable " // TRIM(variablename)
       CALL interpweight_2D(nbpt, classnb, variabletypevals, lalo, resolution, neighbours,          &
            contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype,          &
            maskvals, namemaskvar, 0, 0, -1, fractype,                                                        &
            -1., -1., soilcolrefrac, asoilcol)
       IF (printlev_loc >= 5) WRITE(numout,*)'  albedo_surface_soilalb after interpweight_2D'
       
       ! Check how many points with soil information are found
       nbexp = 0
       
       soilalb_dry(:,:) = zero
       soilalb_wet(:,:) = zero
       soilalb_moy(:,:) = zero
       IF (printlev_loc >= 5) THEN
          WRITE(numout,*)'  albedo_surface_soilalb before starting loop nbpt:', nbpt
          WRITE(numout,*)'  albedo_surface_soilalb initial values classnb: ',classnb
          WRITE(numout,*)'  albedo_surface_soilalb vis_dry. SUM:',SUM(vis_dry),' vis_dry= ',vis_dry
          WRITE(numout,*)'  albedo_surface_soilalb nir_dry. SUM:',SUM(nir_dry),' nir_dry= ',nir_dry
          WRITE(numout,*)'  albedo_surface_soilalb vis_wet. SUM:',SUM(vis_wet),' vis_wet= ',vis_wet
          WRITE(numout,*)'  albedo_surface_soilalb nir_wet. SUM:',SUM(nir_wet),' nir_wet= ',nir_wet
       END IF
       
       DO ib=1,nbpt ! Loop over domain size
          
          ! vecpos: List of positions where textures were not zero
          ! vecpos(1): number of not null textures found
          vecpos = interpweight_ValVecR(soilcolrefrac(ib,:),classnb,zero,'neq')
          fopt = vecpos(1)
          IF (fopt == classnb) THEN
             ! All textures are not zero
             solt(:) = (/(i,i=1,classnb)/)
          ELSE IF (fopt == 0) THEN
             IF (printlev>=1) WRITE(numout,*)'  albedo_surface_soilalb: for point=', ib, ' no soil class!'
          ELSE
             DO ip = 1,fopt
                solt(ip) = vecpos(ip+1)
             END DO
          END IF
          
          !! 3. Compute the average bare soil albedo parameters
          
          IF ( (fopt .EQ. 0) .OR. (asoilcol(ib) .LT. min_sechiba)) THEN
             ! Initialize with mean value if no points were interpolated or if no data was found
             nbexp = nbexp + 1
             soilalb_dry(ib,ivis) = (SUM(vis_dry)/classnb + SUM(vis_wet)/classnb)/deux
             soilalb_dry(ib,inir) = (SUM(nir_dry)/classnb + SUM(nir_wet)/classnb)/deux
             soilalb_wet(ib,ivis) = (SUM(vis_dry)/classnb + SUM(vis_wet)/classnb)/deux
             soilalb_wet(ib,inir) = (SUM(nir_dry)/classnb + SUM(nir_wet)/classnb)/deux
             soilalb_moy(ib,ivis) = SUM(albsoil_vis)/classnb
             soilalb_moy(ib,inir) = SUM(albsoil_nir)/classnb
          ELSE          
             ! If points were interpolated
             DO ip=1, fopt
                IF ( solt(ip) .LE. classnb) THEN
                   ! Set to zero if the value is below min_sechiba
                   IF (soilcolrefrac(ib,solt(ip)) < min_sechiba) soilcolrefrac(ib,solt(ip)) = zero
                   
                   soilalb_dry(ib,ivis) = soilalb_dry(ib,ivis) + vis_dry(solt(ip))*soilcolrefrac(ib,solt(ip))
                   soilalb_dry(ib,inir) = soilalb_dry(ib,inir) + nir_dry(solt(ip))*soilcolrefrac(ib,solt(ip))
                   soilalb_wet(ib,ivis) = soilalb_wet(ib,ivis) + vis_wet(solt(ip))*soilcolrefrac(ib,solt(ip))
                   soilalb_wet(ib,inir) = soilalb_wet(ib,inir) + nir_wet(solt(ip))*soilcolrefrac(ib,solt(ip))
                   soilalb_moy(ib,ivis) = soilalb_moy(ib,ivis) + albsoil_vis(solt(ip))*                    &
                        soilcolrefrac(ib,solt(ip))
                   soilalb_moy(ib,inir) = soilalb_moy(ib,inir) + albsoil_nir(solt(ip))*                    &
                        soilcolrefrac(ib,solt(ip))
                ELSE
                   msg = 'The file contains a soil color class which is incompatible with this program'
                   CALL ipslerr_p(3,'albedo_surface_soilalb',TRIM(msg),'','')
                ENDIF
             ENDDO
          ENDIF
          
       ENDDO
       
       IF ( nbexp .GT. 0 .AND. printlev>=1) THEN
          WRITE(numout,*) 'albedo_surface_soilalb _______'
          WRITE(numout,*) 'albedo_surface_soilalb: The interpolation of the bare soil albedo had ', nbexp
          WRITE(numout,*) 'albedo_surface_soilalb: points without data. This are either coastal points or'
          WRITE(numout,*) 'albedo_surface_soilalb: ice covered land.'
          WRITE(numout,*) 'albedo_surface_soilalb: The problem was solved by using the average of all soils'
          WRITE(numout,*) 'albedo_surface_soilalb: in dry and wet conditions'
          WRITE(numout,*) 'albedo_surface_soilalb: Use the diagnostic output field interp_avail_asoilcol to see location of these points'
       ENDIF
       
       DEALLOCATE (soilcolrefrac)
       DEALLOCATE (variabletypevals)
       
    END IF
    
    ! Write diagnostics
    CALL xios_orchidee_send_field("interp_avail_asoilcol",asoilcol)
    CALL xios_orchidee_send_field("soilalb_dry",soilalb_dry)
    CALL xios_orchidee_send_field("soilalb_wet",soilalb_wet)
    CALL xios_orchidee_send_field("soilalb_moy",soilalb_moy)
    
    IF (printlev_loc >= 3) WRITE(numout,*)'  albedo_surface_soilalb ended'
    
  END SUBROUTINE albedo_surface_soilalb
  
! ----------------------------------------------------------------------------------------

END MODULE albedo_surface
