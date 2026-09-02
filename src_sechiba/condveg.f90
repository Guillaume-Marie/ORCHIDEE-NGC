! ===============================================================================================================================
! MODULE       : condveg
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        Initialise, compute and update the surface parameters emissivity and roughness. 
!!
!! \n DESCRIPTION : This module calculates the emissivity and roughness. It calls albedo_surface_main for the albedo calculations.
!! The module uses 2 settings to control its flow:\n
!! 1. :: rough_dyn to choose between two methods to calculate 
!!    the roughness height. If set to false: the roughness height is 
!!    calculated by the old formulation which does not distinguish between 
!!    z0m and z0h and which does not vary with LAI. If set to true: the 
!!    grid average is calculated by the formulation proposed by Su et al. (2001)
!! 2. :: impaze for choosing surface parameters. If set to false, the 
!!    values for the soil albedo, emissivity and roughness height are 
!!    set to default values which are read from the run.def. If set to 
!!    true, the user imposes its own values, fixed for the grid point. 
!!    This is useful if one performs site simulations, however, 
!!    it is not recommended to do so for spatialized simulations.
!!    roughheight_scal imposes the roughness height in (m) , 
!!    same for emis_scal (in %), albedo_scal (in %), zo_scal (in m)                       
!!    Note that these values are only used if 'impaze' is true.\n
!!
!!    The surface fluxes are calculated between two levels: the 
!!    atmospheric level reference and the effective roughness height 
!!    defined as the difference between the mean height of the vegetation 
!!    and the displacement height (zero wind level). Over bare soils, the 
!!    zero wind level is equal to the soil roughness. Over vegetation, 
!!    the zero wind level is increased by the displacement height
!!    which depends on the height of the vegetation. For a grid point
!!    composed of different types of vegetation, an effective surface 
!!    roughness has to be calculated
!!
!! RECENT CHANGE(S): Added option rough_dyn and subroutine condveg_z0cdrag_dyn. Removed 
!!                   subroutine condveg_z0logz. June 2016.
!!                   The condveg_albedo scheme has been replaced by albedo_surface_main, 
!!                   developed for the DOFOCO branch. The albedo_surface_main is contained
!!                   in albedo.f90 module.
!! REFERENCES(S)    : None
!!

!! SVN              :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_sechiba/condveg.f90 $
!! $Date: 2025-05-26 17:13:29 +0200 (lun. 26 mai 2025) $
!! $Revision: 8972 $
!! \n
!_ ================================================================================================================================

MODULE condveg

  USE ioipsl
  USE xios_orchidee
  USE constantes
  USE constantes_soil
  USE pft_parameters
  USE qsat_moisture
  USE interpol_help
  USE mod_orchidee_para
  USE ioipsl_para 
  USE sechiba_io_p
  USE albedo_surface
  USE structures
  USE grid
  USE function_library,    ONLY: get_printlev

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: condveg_xios_initialize, condveg_main, condveg_initialize, condveg_finalize, condveg_clear 

  !
  ! Variables used inside condveg module
  !
  LOGICAL, SAVE                :: l_first_condveg=.TRUE.              !! To keep first call's trace
!$OMP THREADPRIVATE(l_first_condveg)

  INTEGER, SAVE                                     :: printlev_loc   !! Output debug level
!$OMP THREADPRIVATE(printlev_loc)

CONTAINS

  !!  =============================================================================================================================
  !! SUBROUTINE:    condveg_xios_initialize
  !!
  !>\BRIEF	  Initialize xios dependant defintion before closing context defintion
  !!
  !! DESCRIPTION:	  Initialize xios dependant defintion needed for the interpolations done in condveg.
  !!                      Reading is deactivated if the sechiba restart file exists because the variable
  !!                      should be in the restart file already.
  !!                      This subruting is called before closing context with xios_orchidee_close_definition in intersurf 
  !!                      via the subroutine sechiba_xios_initialize. 
  !!
  !! \n
  !_ ==============================================================================================================================

  SUBROUTINE condveg_xios_initialize
    
    CALL albedo_surface_xios_initialize()

  END SUBROUTINE condveg_xios_initialize

!!  =============================================================================================================================
!! SUBROUTINE		 		    : condveg_initialize
!!
!>\BRIEF			            Allocate module variables, read from restart file or initialize with default values
!!
!! DESCRIPTION			    : Allocate module variables, read from restart file or initialize with default values.
!!                                          condveg_snow is called to initialize corresponding variables.
!!
!! RECENT CHANGE(S)			    : None
!!
!! MAIN OUTPUT VARIABLE(S)
!!
!! REFERENCE(S)			            : None
!! 
!! FLOWCHART                                : None
!! \n
!_ ==============================================================================================================================
  !+++CHECK+++
 ! Variable dimensions xxx_Tot_mean are for a single pixel. Causing 1+1 problems
 ! when running a larger domain. Needs to be corrected when implementing
 ! a global use of the multi-layer energy budget.
  SUBROUTINE condveg_initialize (kjit, kjpindex, index, rest_id, &
       lalo, neighbours, resolution, contfrac, veget, veget_max, frac_nobio, totfrac_nobio, &
       zlev, snow, snow_age, snow_nobio_age, &
       drysoil_frac, height, height_dom, snowdz, snowrho, tot_bare_soil, &
       temp_air, pb, u, v, &
       lai, & 
       emis, albedo, z0m, z0h, roughheight, &
       frac_snow_veg,frac_snow_nobio, &
       coszang, &
       Light_Abs_Tot, Light_Tran_Tot, laieff_fit, &
       Light_Abs_Tot_mean, Light_Alb_Tot_mean, &
       laieff_isotrop)
    !+++++++++++
    
    !! 0. Variable and parameter declaration
    !! 0.1 Input variables  
    INTEGER(i_std), INTENT(in)                       :: kjit             !! Time step number 
    INTEGER(i_std), INTENT(in)                       :: kjpindex         !! Domain size
    INTEGER(i_std),INTENT (in)                       :: rest_id          !! _Restart_ file identifier
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in) :: index            !! Indeces of the points on the map
    REAL(r_std),DIMENSION (kjpindex,2), INTENT (in)  :: lalo             !! Geographical coordinates
    INTEGER(i_std),DIMENSION (kjpindex,NbNeighb), INTENT(in):: neighbours!! neighoring grid points if land
    REAL(r_std), DIMENSION (kjpindex,2), INTENT(in)  :: resolution       !! size in x an y of the grid (m)
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)    :: contfrac         ! Fraction of land in each grid box.
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(in) :: veget            !! Fraction of vegetation types
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(in) :: veget_max        !! Fraction of vegetation type
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT(in) :: frac_nobio    !! Fraction of continental ice, lakes, ...
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)     :: totfrac_nobio    !! total fraction of continental ice+lakes+...
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)    :: zlev             !! Height of first layer
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)     :: snow             !! Snow mass [Kg/m^2]
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)     :: snow_age         !! Snow age
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT(in) :: snow_nobio_age   !! Snow age on ice, lakes, ...
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)     :: drysoil_frac     !! Fraction of visibly Dry soil(between 0 and 1)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(in) :: height           !! Vegetation Height (m)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(in) :: height_dom       !! Dominant vegetation height (m)
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT(in):: snowdz           !! Snow depth at each snow layer (m)
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT(in):: snowrho          !! Snow density at each snow layer (Kg/m^3) 
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)    :: tot_bare_soil    !! Total evaporating bare soil fraction 
    REAL(r_std),DIMENSION(kjpindex),INTENT(in)       :: temp_air         !! Air temperature
    REAL(r_std), DIMENSION(kjpindex), INTENT(in)     :: pb               !! Surface pressure (hPa)
    REAL(r_std),DIMENSION(kjpindex),INTENT(in)       :: u                !! Horizontal wind speed, u direction 
    REAL(r_std),DIMENSION(kjpindex),INTENT(in)       :: v                !! Horizontal wind speed, v direction
                                                                         !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(in) :: lai              !! Leaf area index (m2[leaf]/m2[ground])
    REAL(r_std), DIMENSION(kjpindex), INTENT(in)     :: coszang          !! the cosine of the solar zenith angle (unitless)
    TYPE(laieff_type),DIMENSION (:,:,:),INTENT(in)   :: laieff_fit       !! Fitted parameters for the effective LAI

    !! 0.2 Output variables
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)   :: emis             !! Emissivity
    REAL(r_std),DIMENSION (kjpindex,2), INTENT (out) :: albedo           !! Albedo, vis(1) and nir(2)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)   :: z0m              !! Roughness for momentum (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)   :: z0h              !! Roughness for heat (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)   :: roughheight      !! Effective height for roughness
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)    :: frac_snow_veg    !! Snow cover fraction on vegeted area
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT(out):: frac_snow_nobio  !! Snow cover fraction on non-vegeted area
   REAL(r_std),DIMENSION (kjpindex,nvm,nlevels_tot), &
                                         INTENT(out) :: Light_Abs_Tot    !! Absorbed radiation per layer, averaged between light
                                                                         !! from a collimated (direct) source and that from an
                                                                         !! isotropic (diffuse) source.  Expressed as the fraction
                                                                         !! of overall light hitting the canopy.
   REAL(r_std),DIMENSION (kjpindex,nvm,nlevels_tot), &
                                         INTENT(out) :: Light_Tran_Tot   !! Transmitted radiation per layer, averaged between light
                                                                         !! from a collimated (direct) source and that from an
                                                                         !! isotropic (diffuse) source.  Expressed as the fraction
                                                                         !! of overall light hitting the canopy.
   !+++CHECK+++
   ! Variable dimensions xxx_Tot_mean are for a single pixel. Causing 1+1 problems
   ! when running a larger domain. Needs to be corrected when implementing
   ! a global use of the multi-layer energy budget.
   REAL(r_std),DIMENSION(nlevels_tot), INTENT (out)  :: Light_Abs_Tot_mean  !! total absorption for a given level 
   REAL(r_std),DIMENSION(nlevels_tot), INTENT (out)  :: Light_Alb_Tot_mean  !! total albedo for a given level 
   !+++++++++++
    REAL(r_std), DIMENSION(:,:,:), INTENT(out)       :: laieff_isotrop   !! Leaf Area Index Effective converts
                                                                         !! 3D lai inot 1D lai for two steam
                                                                         !! radiation transfer model...this is for 
                                                                         !! isotropic light and only calculated once per day
                                                                         !! @tex $(m^{2} m^{-2})$ @endtex

    !! 0.4 Local variables   
    INTEGER                                          :: ier
    INTEGER                                          :: ji,jv

!_ ================================================================================================================================
  
    IF (.NOT. l_first_condveg) CALL ipslerr_p(3,'condveg_initialize','Error: initialization already done','','')
    l_first_condveg=.FALSE.

    !! Initialize local printlev
    printlev_loc=get_printlev('condveg')
    

    IF (printlev>=3) WRITE (numout,*) 'Start condveg_initialize'
    
    !! 1. Allocate module variables and read from restart or initialize

    ! z0m
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Roughness for momentum')
    CALL restget_p (rest_id, 'z0m', nbp_glo, 1, 1, kjit, .TRUE., z0m, "gather", nbp_glo, index_g)

    ! z0h
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Roughness for heat')
    CALL restget_p (rest_id, 'z0h', nbp_glo, 1, 1, kjit, .TRUE., z0h, "gather", nbp_glo, index_g)

    ! roughness height
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Roughness height')
    CALL restget_p (rest_id, 'roughheight', nbp_glo, 1, 1, kjit, .TRUE., roughheight, "gather", nbp_glo, index_g)

    !! Initialize emissivity
    IF ( impaze ) THEN
       ! Use parameter CONDVEG_EMIS from run.def
       emis(:) = emis_scal
    ELSE
       ! Set emissivity to 1.
       emis_scal = un
       emis(:) = emis_scal
    ENDIF


    !! 3. Calculate the fraction of snow on vegetation and nobio
    CALL condveg_frac_snow(kjpindex, snowrho, snowdz, &
       frac_snow_veg, frac_snow_nobio, totfrac_nobio)

    !! 4. Calculate roughness height if it was not found in the restart file
    IF ( ALL(z0m(:) == val_exp) .OR. ALL(z0h(:) == val_exp) .OR. ALL(roughheight(:) == val_exp)) THEN
       !! Calculate roughness height
       ! Chooses between two methods to calculate the grid average of the 
       ! roughness. If impaze set to true:  The grid average is calculated by 
       ! averaging the drag coefficients over PFT. If impaze set to false: 
       ! The grid average is calculated by averaging the logarithm of the 
       ! roughness length per PFT.
       IF ( impaze ) THEN
          ! Use parameter CONDVEG_Z0 and ROUGHHEIGHT from run.def
          z0m(:) = z0_scal
          z0h(:) = z0_scal
          roughheight(:) = roughheight_scal
       ELSE
          ! Caluculate roughness height
          IF( rough_dyn ) THEN
             CALL condveg_z0cdrag_dyn(kjpindex, veget, veget_max, frac_nobio, totfrac_nobio, zlev, &
                  &               height, height_dom, temp_air, pb, u, v, lai, frac_snow_veg, z0m, z0h, roughheight)
          ELSE
             CALL condveg_z0cdrag(kjpindex, veget, veget_max, frac_nobio, totfrac_nobio, zlev, &
                  height, height_dom, tot_bare_soil, frac_snow_veg, z0m, z0h, roughheight)
          ENDIF
       END IF
    END IF


    !! 5. Initialze albedo module
    !+++CHECK+++
    ! Variable dimensions xxx_Tot_mean are for a single pixel. Causing 1+1 problems
    ! when running a larger domain. Needs to be corrected when implementing
    ! a global use of the multi-layer energy budget.
    CALL albedo_surface_initialize(kjit,  kjpindex,  rest_id, &
                           lalo, neighbours, resolution, contfrac, &
                           drysoil_frac, veget_max, coszang, frac_nobio, & 
                           snow, snow_age, &
                           snow_nobio_age, frac_snow_veg, frac_snow_nobio,&
                           z0m, laieff_fit, &
                           albedo, &
                           Light_Abs_Tot, Light_Tran_Tot, &
                           Light_Abs_Tot_mean, Light_Alb_Tot_mean, &
                           laieff_isotrop, veget)
    !+++++++++++
         
    IF (printlev>=3) WRITE (numout,*) 'condveg_initialize done '
    
  END SUBROUTINE condveg_initialize



!! ==============================================================================================================================
!! SUBROUTINE   : condveg_main
!!
!>\BRIEF        Calls the subroutines update the variables for current time step
!!
!!
!! MAIN OUTPUT VARIABLE(S):  emis (emissivity), albedo (albedo of 
!! vegetative PFTs in visible and near-infrared range), z0 (surface roughness height),
!! roughheight (grid effective roughness height), soil type (fraction of soil types) 
!! 
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!!
!! REVISION(S)  : None
!!
!_ ================================================================================================================================

  !+++CHECK+++
  ! Variable dimensions xxx_Tot_mean are for a single pixel. Causing 1+1 problems
  ! when running a larger domain. Needs to be corrected when implementing
  ! a global use of the multi-layer energy budget.
  SUBROUTINE condveg_main (kjit, kjpindex, index, rest_id, hist_id, hist2_id, &
       lalo, neighbours, resolution, contfrac, veget, veget_max, frac_nobio, totfrac_nobio, &
       zlev, snow, snow_age, snow_nobio_age, &
       drysoil_frac, height, height_dom, snowdz, snowrho, tot_bare_soil, &
       Temp_air, pb, u, v, lai,& 
       emis, albedo, z0m, z0h, roughheight, &
       frac_snow_veg, frac_snow_nobio, coszang, &
       Light_Abs_Tot, Light_Tran_Tot, laieff_fit, &
       Light_Abs_Tot_mean, Light_Alb_Tot_mean, &
       laieff_isotrop )
    !+++++++++++

     !! 0. Variable and parameter declaration

    !! 0.1 Input variables  

    INTEGER(i_std), INTENT(in)                       :: kjit             !! Time step number 
    INTEGER(i_std), INTENT(in)                       :: kjpindex         !! Domain size
    INTEGER(i_std),INTENT (in)                       :: rest_id          !! _Restart_ file identifier
    INTEGER(i_std),INTENT (in)                       :: hist_id          !! _History_ file identifier
    INTEGER(i_std), OPTIONAL, INTENT (in)            :: hist2_id          !! _History_ file 2 identifier
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in) :: index            !! Indeces of the points on the map
    REAL(r_std),DIMENSION (kjpindex,2), INTENT (in)  :: lalo             !! Geographical coordinates
    INTEGER(i_std),DIMENSION (kjpindex,NbNeighb), INTENT(in):: neighbours!! neighoring grid points if land
    REAL(r_std), DIMENSION (kjpindex,2), INTENT(in)  :: resolution       !! size in x an y of the grid (m)
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)    :: contfrac         ! Fraction of land in each grid box.
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(in) :: veget            !! Fraction of vegetation types
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(in) :: veget_max        !! Fraction of vegetation type
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT(in) :: frac_nobio    !! Fraction of continental ice, lakes, ...
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)     :: totfrac_nobio    !! total fraction of continental ice+lakes+...
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)    :: zlev             !! Height of first layer
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)     :: snow             !! Snow mass [Kg/m^2]
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)     :: snow_age         !! Snow age
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT(in) :: snow_nobio_age   !! Snow age on ice, lakes, ...
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)     :: drysoil_frac     !! Fraction of visibly Dry soil(between 0 and 1)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(in) :: height           !! Vegetation Height (m)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(in) :: height_dom       !! Dominant vegetation height (m)
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT(in):: snowdz           !! Snow depth at each snow layer (m)
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT(in):: snowrho          !! Snow density at each snow layer (Kg/m^3) 
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)    :: tot_bare_soil    !! Total evaporating bare soil fraction 
    REAL(r_std),DIMENSION(kjpindex),INTENT(in)       :: temp_air         !! Air temperature
    REAL(r_std), DIMENSION(kjpindex), INTENT(in)     :: pb               !! Surface pressure (hPa)
    REAL(r_std),DIMENSION(kjpindex),INTENT(in)       :: u                !! Horizontal wind speed, u direction 
    REAL(r_std),DIMENSION(kjpindex),INTENT(in)       :: v                !! Horizontal wind speed, v direction
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)     :: coszang          !! the cosine of the solar enith angle (unitless)
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(in) :: lai              !! Leaf area index (m2[leaf]/m2[ground])
    TYPE(laieff_type),DIMENSION (:,:,:),INTENT(in)   :: laieff_fit       !! Fitted parameters for the effective LAI


    !! 0.2 Output variables

    REAL(r_std),DIMENSION (kjpindex), INTENT (out)   :: emis             !! Emissivity
    REAL(r_std),DIMENSION (kjpindex,2), INTENT (out) :: albedo           !! Albedo, vis(1) and nir(2)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)   :: z0m              !! Roughness for momentum (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)   :: z0h              !! Roughness for heat (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)   :: roughheight      !! Effective height for roughness
    REAL(r_std),DIMENSION (:,:,:), &
                                      INTENT (out)   :: Light_Abs_Tot !!Absorbed radiation per level for photosynthesis
    REAL(r_std),DIMENSION (:,:,:), &
                                      INTENT (out)   :: Light_Tran_Tot !!Transmitted radiation per level for photosynthesis
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)    :: frac_snow_veg    !! Snow cover fraction on vegeted area
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT(out):: frac_snow_nobio  !! Snow cover fraction on non-vegeted area

    !! 0.3 Modified variables
    !+++CHECK+++
    ! Variable dimensions xxx_Tot_mean are for a single pixel. Causing 1+1 problems
    ! when running a larger domain. Needs to be corrected when implementing
    ! a global use of the multi-layer energy budget.
    REAL(r_std),DIMENSION(:), INTENT(inout)          :: Light_Abs_Tot_mean    !! total light absorption for a given level
    REAL(r_std),DIMENSION(:), INTENT(inout)          :: Light_Alb_Tot_mean   !! total albedo for a given level
    !+++++++++++
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)     :: laieff_isotrop   !! Leaf Area Index Effective converts
                                                                         !! 3D lai into 1D lai for two stream
                                                                         !! radiation transfer model...this is for
                                                                         !! isotropic light and only calculated once per day
                                                                         !! @tex $(m^{2} m^{-2})$ @endtex

    !! 0.4 Local variables 
    CHARACTER(LEN=80)                                :: var_name         !! To store variables names for I/O  
    INTEGER(i_std)                                   :: ji,ilevel,ivm

!_ ================================================================================================================================
  
    !! 1. Calculate the fraction of snow on vegetation and nobio
    CALL condveg_frac_snow(kjpindex, snowrho, snowdz, &
            frac_snow_veg, frac_snow_nobio, totfrac_nobio)

    !! 2. Calculate emissivity
    emis(:) = emis_scal
    
    !! 3. Calculate roughness height
    ! If TRUE read in prescribed values for roughness height
    IF ( impaze ) THEN

       DO ji = 1, kjpindex
         z0m(ji) = z0_scal
         z0h(ji) = z0_scal
         roughheight(ji) = roughheight_scal
      ENDDO

    ELSE

       ! Calculate roughness height
       IF ( rough_dyn ) THEN
          CALL condveg_z0cdrag_dyn (kjpindex, veget, veget_max, &
               frac_nobio, totfrac_nobio, zlev, height, height_dom, &
               temp_air, pb, u, v, lai, frac_snow_veg, &
               z0m, z0h, roughheight)
       ELSE
          CALL condveg_z0cdrag (kjpindex, veget, veget_max, &
               frac_nobio, totfrac_nobio, zlev, &
               height, height_dom, tot_bare_soil, frac_snow_veg, z0m, z0h, roughheight)
       ENDIF
     
    ENDIF ! impaze


    !! 4. Calculate albedo
    !+++CHECK+++
    ! Variable dimensions xxx_Tot_mean are for a single pixel. Causing 1+1 problems
    ! when running a larger domain. Needs to be corrected when implementing
    ! a global use of the multi-layer energy budget.
    CALL albedo_surface_main(                                                          &
         kjit,                kjpindex,            lalo,            hist_id, hist2_id, &
         index,               drysoil_frac,        veget_max,       coszang,           &
         frac_nobio,          frac_snow_veg,       frac_snow_nobio,                    &
         snow,                snow_age,            snow_nobio_age,                     &
         z0m,                 laieff_fit,                                              &
         albedo,              Light_Abs_Tot,       Light_Tran_Tot,                     &
         Light_Abs_Tot_mean,  Light_Alb_Tot_mean,                                      &
         laieff_isotrop, veget)
    !+++++++++++

    ! Debug
    IF(printlev >= 4)THEN
       DO ivm = 1,nvm
          WRITE(numout,*) 'laieff_isotrop condveg_main for ivm', &
              ivm,'is', laieff_isotrop(:,:,ivm)
       ENDDO
    ENDIF
    !-


    IF (printlev>=3) WRITE (numout,*)' condveg_main done '

  END SUBROUTINE condveg_main

  !!  =============================================================================================================================
  !! SUBROUTINE		 		    : condveg_finalize
  !!
  !>\BRIEF                                    Write to restart file
  !!
  !! DESCRIPTION			    : This subroutine writes the module variables and variables calculated in condveg
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
  SUBROUTINE condveg_finalize (kjit,                kjpindex,            rest_id,                          &
                               z0m,                 z0h,                 roughheight,                      &
                               albedo,              Light_Abs_Tot,       Light_Tran_Tot,               &
                               Light_Abs_Tot_mean,  Light_Alb_Tot_mean,  &
                               laieff_isotrop)
    !+++++++++++

    !! 0. Variable and parameter declaration
    !! 0.1 Input variables  
    INTEGER(i_std), INTENT(in)                  :: kjit             !! Time step number 
    INTEGER(i_std), INTENT(in)                  :: kjpindex         !! Domain size
    INTEGER(i_std),INTENT (in)                  :: rest_id          !! Restart file identifier
    REAL(r_std),DIMENSION(kjpindex), INTENT(in) :: z0m              !! Roughness for momentum
    REAL(r_std),DIMENSION(kjpindex), INTENT(in) :: z0h              !! Roughness for heat
    REAL(r_std),DIMENSION(kjpindex), INTENT(in) :: roughheight      !! Grid effective roughness height (m)     
    REAL(r_std), DIMENSION(kjpindex,n_spectralbands), &
                                    INTENT(in)  :: albedo           !! Albedo (two stream radiation transfer model)
    REAL(r_std),DIMENSION (kjpindex,nvm,nlevels_tot), &
                                    INTENT(in) :: Light_Abs_Tot     !! Absorbed radiation per layer, averaged between light
                                                                    !! from a collimated (direct) source and that from an
                                                                    !! isotropic (diffuse) source.  Expressed as the fraction
                                                                    !! of overall light hitting the canopy.
    REAL(r_std),DIMENSION (kjpindex,nvm,nlevels_tot), &
                                    INTENT(in) :: Light_Tran_Tot    !! Transmitted radiation per layer, averaged between light
                                                                    !! from a collimated (direct) source and that from an
                                                                    !! isotropic (diffuse) source.  Expressed as the fraction
                                                                    !! of overall light hitting the canopy.
    !+++CHECK+++
    ! Variable dimensions xxx_Tot_mean are for a single pixel. Causing 1+1 problems
    ! when running a larger domain. Needs to be corrected when implementing
    ! a global use of the multi-layer energy budget.
    REAL(r_std),DIMENSION(nlevels_tot), INTENT(in)                :: Light_Abs_Tot_mean    !! total light absorption for a given level
    REAL(r_std),DIMENSION(nlevels_tot), INTENT(in)                :: Light_Alb_Tot_mean    !! total albedo for a given level
    !+++++++++++
    REAL(r_std), DIMENSION(kjpindex,nlevels_tot,nvm), INTENT(in)  :: laieff_isotrop     !! Leaf Area Index Effective converts
                                                                                        !! 3D lai into 1D lai for two stream
    
    !_ ================================================================================================================================
    CALL restput_p (rest_id, 'z0m', nbp_glo, 1, 1, kjit, z0m, 'scatter',  nbp_glo, index_g)
    !-
    CALL restput_p (rest_id, 'z0h', nbp_glo, 1, 1, kjit, z0h, 'scatter',  nbp_glo, index_g)
    !-
    CALL restput_p (rest_id, 'roughheight', nbp_glo, 1, 1, kjit, roughheight, 'scatter',  nbp_glo, index_g)

    ! Finalize albedo module
    !+++CHECK+++
    ! Variable dimensions xxx_Tot_mean are for a single pixel. Causing 1+1 problems
    ! when running a larger domain. Needs to be corrected when implementing
    ! a global use of the multi-layer energy budget.
    CALL albedo_surface_finalize (kjit,                kjpindex,            rest_id,                          &
                          albedo,              Light_Abs_Tot,   Light_Tran_Tot,               &
                          Light_Abs_Tot_mean, Light_Alb_Tot_mean, &
                          laieff_isotrop)
    !+++++++++++
 
  END SUBROUTINE condveg_finalize

!! ==============================================================================================================================
!! SUBROUTINE 	: condveg_clear
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

  SUBROUTINE condveg_clear  ()

      l_first_condveg=.TRUE.
       
      CALL albedo_surface_clear()

  END SUBROUTINE condveg_clear

!! ==============================================================================================================================
!! SUBROUTINE   : condveg_frac_snow
!!
!>\BRIEF        This subroutine calculates the fraction of snow on vegetation and nobio
!!
!! DESCRIPTION  
!!
!! RECENT CHANGE(S): These calculations were previously done in condveg_snow.
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  
  SUBROUTINE condveg_frac_snow(kjpindex, snowrho, snowdz, &
                 frac_snow_veg, frac_snow_nobio, totfrac_nobio)

    !! 0. Variable and parameter declaration
    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                          :: kjpindex        !! Domain size
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT(in)   :: snowrho         !! Snow density at each snow layer (Kg/m^3) 
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT(in)   :: snowdz          !! Snow depth at each snow layer (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)        :: totfrac_nobio    !! total fraction of continental ice+lakes+...

    !! 0.2 Output variables
    REAL(r_std), DIMENSION(kjpindex), INTENT(out)       :: frac_snow_veg   !! Fraction of snow on vegetation (unitless ratio)
    REAL(r_std), DIMENSION(kjpindex,nnobio), INTENT(out):: frac_snow_nobio !! Fraction of snow on continental ice, lakes, etc. 

    !! 0.3 Local variables
    REAL(r_std), DIMENSION(kjpindex)                    :: snowrho_ave     !! Average snow density
    REAL(r_std), DIMENSION(kjpindex)                    :: snowdepth       !! Snow depth
    REAL(r_std), DIMENSION(kjpindex)                    :: snowrho_snowdz       !! Snow rho time snowdz
    INTEGER(i_std)                                      :: jv
   

!_ ================================================================================================================================

    ! Calculate snow cover fraction for both total vegetated and 
    ! total non-vegetative surfaces.
    snowdepth=sum(snowdz,2)
    snowrho_snowdz=sum(snowrho*snowdz,2)
    WHERE(snowdepth(:) .LT. min_sechiba)
       frac_snow_veg(:) = 0.
    ELSEWHERE
       snowrho_ave(:)=snowrho_snowdz(:)/snowdepth(:)
       frac_snow_veg(:) = tanh(snowdepth(:)/(0.025*(snowrho_ave(:)/50.)))
    END WHERE

    WHERE (totfrac_nobio(:) .EQ. 1) 
	frac_snow_veg(:)=0.
    END WHERE

    frac_snow_nobio(:,:)=0.
    snowdepth=SUM(snowdz,2)
    snowrho_snowdz=sum(snowrho*snowdz,2)
    WHERE(snowdepth(:) .LT. min_sechiba)
      frac_snow_nobio(:,iice) = 0.
    ELSEWHERE
      snowrho_ave(:)=snowrho_snowdz(:)/snowdepth(:)
      frac_snow_nobio(:,iice) = tanh(snowdepth(:)/(0.025*(snowrho_ave(:)/50.)))
    ENDWHERE

    IF (printlev>=3) WRITE (numout,*) ' condveg_frac_snow done '
    
  END SUBROUTINE condveg_frac_snow

  

!! ==============================================================================================================================
!! SUBROUTINE   : condveg_z0cdrag
!!
!>\BRIEF        Computation of grid average of roughness length by calculating 
!! the drag coefficient.
!!
!! DESCRIPTION  : This routine calculates the mean roughness height and mean 
!! effective roughness height over the grid cell. The mean roughness height (z0) 
!! is computed by averaging the drag coefficients  \n
!!
!! \latexonly 
!! \input{z0cdrag1.tex}
!! \endlatexonly
!! \n 
!!
!! where C is the drag coefficient at the height of the vegetation, kappa is the 
!! von Karman constant, z (Ztmp) is the height at which the fluxes are estimated and z0 the roughness height. 
!! The reference level for z needs to be high enough above the canopy to avoid 
!! singularities of the LOG. This height is set to  minimum 10m above ground. 
!! The drag coefficient increases with roughness height to represent the greater 
!! turbulence generated by rougher surfaces. 
!! The roughenss height is obtained by the inversion of the drag coefficient equation.\n
!!
!! The roughness height for the non-vegetative surfaces is calculated in a second step. 
!! In order to calculate the transfer coefficients the 
!! effective roughness height is calculated. This effective value is the difference
!! between the height of the vegetation and the zero plane displacement height.\nn
!!
!! RECENT CHANGE(S): None
!! 
!! MAIN OUTPUT VARIABLE(S):  :: roughness height(z0) and grid effective roughness height(roughheight)
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE condveg_z0cdrag (kjpindex, veget, veget_max, frac_nobio,&
    totfrac_nobio, zlev, height, height_dom, tot_bare_soil, frac_snow_veg,&
    z0m, z0h, roughheight)

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables
   
    INTEGER(i_std), INTENT(in)                          :: kjpindex      !! Domain size - Number of land pixels  (unitless)
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(in)    :: veget         !! PFT coverage fraction of a PFT (= ind*cn_ind) 
                                                                         !! (m^2 m^{-2})
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(in)    :: veget_max     !! PFT "Maximal" coverage fraction of a PFT 
                                                                         !! (= ind*cn_ind) (m^2 m^{-2})
    REAL(r_std), DIMENSION(kjpindex,nnobio), INTENT(in) :: frac_nobio    !! Fraction of non-vegetative surfaces, 
                                                                         !! i.e. continental ice, lakes, etc. (unitless)
    REAL(r_std), DIMENSION(kjpindex), INTENT(in)        :: totfrac_nobio !! Total fraction of non-vegetative surfaces, 
                                                                         !! i.e. continental ice, lakes, etc. (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: zlev          !! Height of first layer (m)           
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(in)    :: height        !! Vegetation height (m)
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(in)    :: height_dom    !! Dominant vegetation height (m)
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)       :: tot_bare_soil !! Total evaporating bare soil fraction 
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)        :: frac_snow_veg !! Snow cover fraction on vegeted area

    !! 0.2 Output variables

    REAL(r_std), DIMENSION(kjpindex), INTENT(out)       :: z0m           !! Roughness height for momentum (m)
    REAL(r_std), DIMENSION(kjpindex), INTENT(out)       :: z0h           !! Roughness height for heat (m) 
    REAL(r_std), DIMENSION(kjpindex), INTENT(out)       :: roughheight   !! Grid effective roughness height (m) 
    
    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                                      :: jv            !! Loop index over PFTs (unitless)
    REAL(r_std), DIMENSION(kjpindex)                    :: sumveg        !! Fraction of bare soil (unitless)
    REAL(r_std), DIMENSION(kjpindex)                    :: ztmp          !! Max height of the atmospheric level (m)
    REAL(r_std), DIMENSION(kjpindex)                    :: ave_height    !! Average vegetation height (m)
    REAL(r_std), DIMENSION(kjpindex)                    :: d_veg         !! PFT coverage of vegetative PFTs 
                                                                         !! (= ind*cn_ind) (m^2 m^{-2})
    REAL(r_std), DIMENSION(kjpindex)                    :: zhdispl       !! Zero plane displacement height (m)
    REAL(r_std)                                         :: z0_nobio      !! Roughness height of non-vegetative fraction (m),  
                                                                         !! i.e. continental ice, lakes, etc. 
    REAL(r_std), DIMENSION(kjpindex)                    :: dragm         !! Dragcoefficient for momentum
    REAL(r_std), DIMENSION(kjpindex)                    :: dragh         !! Dragcoefficient for heat
    REAL(r_std), DIMENSION(kjpindex)                    :: z0_ground     !! z0m value used for ground surface

!_ ================================================================================================================================
    
    !! 1. Preliminary calculation

    ! Set maximal height of first layer
    ztmp(:) = MAX(10., zlev(:))

    z0_ground(:) = (1.-frac_snow_veg(:))*z0_bare + frac_snow_veg(:)*z0_bare/10.

    ! Calculate roughness for non-vegetative surfaces
    ! with the von Karman constant 
    dragm(:) = tot_bare_soil(:) * (ct_karman/LOG(ztmp(:)/z0_ground))**2
    dragh(:) = tot_bare_soil(:) * (ct_karman/LOG(ztmp(:)/(z0_ground/ratio_z0m_z0h(1))))*(ct_karman/LOG(ztmp(:)/z0_ground))

    ! Fraction of bare soil
    sumveg(:) = tot_bare_soil(:)

    ! Set average vegetation height to zero
    ave_height(:) = zero
    
    !! 2. Calculate the mean roughness height 
    
    ! Calculate the mean roughness height of
    ! vegetative PFTs over the grid cell
    DO jv = 2, nvm

       ! In the case of forest, use parameter veget_max because 
       ! tree trunks influence the roughness even when there are no leaves
       IF ( is_tree(jv) ) THEN
          ! In the case of grass, use parameter veget because grasses 
          ! only influence the roughness during the growing season
          d_veg(:) = veget_max(:,jv)
       ELSE
          ! grasses only have an influence if they are really there!
          d_veg(:) = veget(:,jv)
       ENDIF
       
       ! Calculate the average roughness over the grid cell:
       ! The unitless drag coefficient is per vegetative PFT
       ! calculated by use of the von Karman constant, the height 
       ! of the first layer and the roughness. The roughness
       ! is calculated as the vegetation height  per PFT 
       ! multiplied by the roughness  parameter 'z0_over_height= 1/16'. 
       ! If this scaled value is lower than 0.01 then the value for 
       ! the roughness of bare soil (0.01) is used. 
       ! The sum over all PFTs gives the average roughness 
       ! per grid cell for the vegetative PFTs.
       IF(use_height_dom)THEN
          dragm(:) = dragm(:) + d_veg(:) *(ct_karman/LOG(ztmp(:)/MAX(height_dom(:,jv)*z0_over_height(jv),z0_ground)))**2
          dragh(:) = dragh(:) + d_veg(:) *(ct_karman/LOG(ztmp(:)/(MAX(height_dom(:,jv)*z0_over_height(jv),z0_ground) / &
               ratio_z0m_z0h(jv)))) * (ct_karman/LOG(ztmp(:)/MAX(height_dom(:,jv)*z0_over_height(jv),z0_ground)))
       ELSE
          dragm(:) = dragm(:) + d_veg(:) *(ct_karman/LOG(ztmp(:)/MAX(height(:,jv)*z0_over_height(jv),z0_ground)))**2
          dragh(:) = dragh(:) + d_veg(:) *(ct_karman/LOG(ztmp(:)/(MAX(height(:,jv)*z0_over_height(jv),z0_ground) / &
               ratio_z0m_z0h(jv)))) * (ct_karman/LOG(ztmp(:)/MAX(height(:,jv)*z0_over_height(jv),z0_ground)))
       ENDIF

       ! Sum of bare soil and fraction vegetated fraction
       sumveg(:) = sumveg(:) + d_veg(:)
       
       ! Weigh height of vegetation with maximal cover fraction
       IF(use_height_dom)THEN
          ave_height(:) = ave_height(:) + veget_max(:,jv)*height_dom(:,jv)
       ELSE
          ave_height(:) = ave_height(:) + veget_max(:,jv)*height(:,jv)
       ENDIF
    ENDDO
    
    !! 3. Calculate the mean roughness height of vegetative PFTs over the grid cell
    
    !  Search for pixels with vegetated part to normalise 
    !  roughness height
    WHERE ( sumveg(:) .GT. min_sechiba ) 
       dragm(:) = dragm(:) / sumveg(:)
       dragh(:) = dragh(:) / sumveg(:)
    ENDWHERE
    ! Calculate fraction of roughness for vegetated part 
    dragm(:) = (un - totfrac_nobio(:)) * dragm(:)
    dragh(:) = (un - totfrac_nobio(:)) * dragh(:)

    DO jv = 1, nnobio ! Loop over # of non-vegative surfaces

       ! Set rougness for ice
       IF ( jv .EQ. iice ) THEN
          z0_nobio = z0_ice
       ELSE
          WRITE(numout,*) 'jv=',jv
          WRITE(numout,*) 'DO NOT KNOW ROUGHNESS OF THIS SURFACE TYPE'
          CALL ipslerr_p(3,'condveg_z0cdrag','DO NOT KNOW ROUGHNESS OF THIS SURFACE TYPE','','')
       ENDIF
       
       ! Sum of vegetative roughness length and non-vegetative roughness length
       dragm(:) = dragm(:) + frac_nobio(:,jv) *(ct_karman/LOG(ztmp(:)/z0_nobio))**2
       dragh(:) = dragh(:) + frac_nobio(:,jv) *(ct_karman/LOG(ztmp(:)/(z0_nobio/ratio_z0m_z0h(1))))*(ct_karman/LOG(ztmp(:)/z0_nobio))

    ENDDO ! Loop over # of non-vegative surfaces
    
    !! 4. Calculate the zero plane displacement height and effective roughness length

    !  Take the exponential of the roughness
    z0m(:) = ztmp(:) / EXP(ct_karman/SQRT(dragm(:)))
    z0h(:) = ztmp(:) / EXP((ct_karman**2.)/(dragh(:)*LOG(ztmp(:)/z0m(:))))


    ! Compute the zero plane displacement height which
    ! is an equivalent height for the absorption of momentum
    zhdispl(:) = ave_height(:) * height_displacement

    ! In order to calculate the fluxes we compute what we call the grid effective roughness height.
    ! This is the height over which the roughness acts. It combines the
    ! zero plane displacement height and the vegetation height.
    roughheight(:) = ave_height(:) - zhdispl(:)

  END SUBROUTINE condveg_z0cdrag


!! ==============================================================================================================================
!! SUBROUTINE   : condveg_z0cdrag_dyn
!!
!>\BRIEF        Computation of grid average of roughness length by calculating 
!! the drag coefficient based on formulation proposed by Su et al. (2001). 
!!
!! DESCRIPTION  : This routine calculates the mean roughness height and mean 
!! effective roughness height over the grid cell. The mean roughness height (z0) 
!! is computed by averaging the drag coefficients  \n
!!
!! \latexonly 
!! \input{z0cdrag1.tex}
!! \endlatexonly
!! \n 
!!
!! where C is the drag coefficient at the height of the vegetation, kappa is the 
!! von Karman constant, z (Ztmp) is the height at which the fluxes are estimated and z0 the roughness height. 
!! The reference level for z needs to be high enough above the canopy to avoid 
!! singularities of the LOG. This height is set to  minimum 10m above ground. 
!! The drag coefficient increases with roughness height to represent the greater 
!! turbulence generated by rougher surfaces. 
!! The roughenss height is obtained by the inversion of the drag coefficient equation.\n
!! In the formulation of Su et al. (2001), one distinguishes the roughness height for
!! momentum (z0m) and the one for heat (z0h). 
!! z0m is computed as a function of LAI (z0m increases with LAI) and z0h is computed  
!! with a so-called kB-1 term (z0m/z0h=exp(kB-1))
!!
!! RECENT CHANGE(S): Written by N. Vuichard (2016)
!! 
!! MAIN OUTPUT VARIABLE(S):  :: roughness height(z0) and grid effective roughness height(roughheight)
!!
!! REFERENCE(S) : 
!! - Su, Z., Schmugge, T., Kustas, W.P., Massman, W.J., 2001. An Evaluation of Two Models for 
!! Estimation of the Roughness Height for Heat Transfer between the Land Surface and the Atmosphere. J. Appl. 
!! Meteorol. 40, 1933–1951. doi:10.1175/1520-0450(2001)
!! - Ershadi, A., McCabe, M.F., Evans, J.P., Wood, E.F., 2015. Impact of model structure and parameterization 
!! on Penman-Monteith type evaporation models. J. Hydrol. 525, 521–535. doi:10.1016/j.jhydrol.2015.04.008
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE condveg_z0cdrag_dyn (kjpindex,veget,veget_max,frac_nobio,totfrac_nobio,zlev, height, height_dom, &
       &                      temp_air, pb, u, v, lai, frac_snow_veg, z0m, z0h, roughheight)

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables
   
    INTEGER(i_std), INTENT(in)                          :: kjpindex      !! Domain size - Number of land pixels  (unitless)
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(in)    :: veget         !! PFT coverage fraction of a PFT (= ind*cn_ind) 
                                                                         !! (m^2 m^{-2})
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(in)    :: veget_max     !! PFT "Maximal" coverage fraction of a PFT 
                                                                         !! (= ind*cn_ind) (m^2 m^{-2})
    REAL(r_std), DIMENSION(kjpindex,nnobio), INTENT(in) :: frac_nobio    !! Fraction of non-vegetative surfaces, 
                                                                         !! i.e. continental ice, lakes, etc. (unitless)
    REAL(r_std), DIMENSION(kjpindex), INTENT(in)        :: totfrac_nobio !! Total fraction of non-vegetative surfaces, 
                                                                         !! i.e. continental ice, lakes, etc. (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: zlev          !! Height of first layer (m)           
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(in)    :: height        !! Vegetation height (m)    
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(in)    :: height_dom    !! Dominant vegetation height (m)    
    REAL(r_std), DIMENSION(kjpindex), INTENT(in)        :: temp_air      !! 2m air temperature (K)
    REAL(r_std), DIMENSION(kjpindex), INTENT(in)        :: pb            !! Surface pressure (hPa)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: u             !! Lowest level wind speed in direction u 
                                                                         !! @tex $(m.s^{-1})$ @endtex 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: v             !! Lowest level wind speed in direction v 
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(in)    :: lai           !! Leaf area index (m2[leaf]/m2[ground])
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)        :: frac_snow_veg    !! Snow cover fraction on vegeted area
    !! 0.2 Output variables

    REAL(r_std), DIMENSION(kjpindex), INTENT(out)       :: z0m           !! Roughness height for momentum (m)
    REAL(r_std), DIMENSION(kjpindex), INTENT(out)       :: z0h           !! Roughness height for heat (m)
    REAL(r_std), DIMENSION(kjpindex), INTENT(out)       :: roughheight   !! Grid effective roughness height (m) 
    
    !! 0.3 Modified variables

    !! 0.4 Local variables
    INTEGER(i_std)                                      :: ji            !! Loop index over grid points
    INTEGER(i_std)                                      :: jv            !! Loop index over PFTs (unitless)
    REAL(r_std)                                         :: threshold     !! Treshold for calculation of kB_m1
    REAL(r_std), DIMENSION(kjpindex)                    :: sumveg        !! Fraction of bare soil (unitless)
    REAL(r_std), DIMENSION(kjpindex)                    :: ztmp          !! Max height of the atmospheric level (m)
    REAL(r_std), DIMENSION(kjpindex)                    :: ave_height    !! Average vegetation height (m)
    REAL(r_std), DIMENSION(kjpindex)                    :: zhdispl       !! Zero plane displacement height (m)
    REAL(r_std)                                         :: z0_nobio      !! Roughness height of non-vegetative fraction (m),  
    REAL(r_std), DIMENSION(kjpindex)                    :: z0m_pft       !! Roughness height for momentum for a specific PFT
    REAL(r_std), DIMENSION(kjpindex)                    :: z0h_pft       !! Roughness height for heat for a specific PFT
    REAL(r_std), DIMENSION(kjpindex)                    :: dragm         !! Dragcoefficient for momentum
    REAL(r_std), DIMENSION(kjpindex)                    :: dragh         !! Dragcoefficient for heat
    REAL(r_std), DIMENSION(kjpindex)                    :: eta           !! Ratio of friction velocity to the wind speed at the canopy top - See Ershadi et al. (2015)
    REAL(r_std), DIMENSION(kjpindex)                    :: eta_ec        !! Within-canopy wind speed profile estimation coefficient - See Ershadi et al. (2015)
    REAL(r_std), DIMENSION(kjpindex)                    :: Ct_star       !! Heat transfer coefficient of the soil - see Su et al. (2001)
    REAL(r_std), DIMENSION(kjpindex)                    :: kBs_m1        !! Canopy model of Brutsaert (1982) for a bare soil surface - used in the calculation of kB_m1 (see Ershadi et al. (2015))
    REAL(r_std), DIMENSION(kjpindex)                    :: kB_m1         !! kB**-1: Term used in the calculation of z0h where B-1 is the inverse Stanton number (see Ershadi et al. (2015))
    REAL(r_std), DIMENSION(kjpindex)                    :: fc            !! fractional canopy coverage
    REAL(r_std), DIMENSION(kjpindex)                    :: fs            !! fractional soil coverage
    REAL(r_std), DIMENSION(kjpindex)                    :: Reynolds      !! Reynolds number
    REAL(r_std), DIMENSION(kjpindex)                    :: wind          !! wind Speed (m)
    REAL(r_std), DIMENSION(kjpindex)                    :: u_star        !! friction velocity
    REAL(r_std), DIMENSION(kjpindex)                    :: z0_ground     !! z0m value used for ground surface
    REAL(r_std), DIMENSION(kjpindex,nvm)                :: loc_height    !! Vegetation height (m)  

!_ ================================================================================================================================
    
    !! 1. Preliminary calculation

    ! Set maximal height of first layer
    ztmp(:) = MAX(10., zlev(:))
    
    z0_ground(:) = (1.-frac_snow_veg(:))*z0_bare + frac_snow_veg(:)*z0_bare/10.

    ! Calculate roughness for non-vegetative surfaces
    ! with the von Karman constant 
    dragm(:) = veget_max(:,1) * (ct_karman/LOG(ztmp(:)/z0_ground(:)))**2

    wind(:) = SQRT(u(:)*u(:)+v(:)*v(:))
    u_star(:)= ct_karman * MAX(min_wind,wind(:)) / LOG(zlev(:)/z0_ground(:))
    Reynolds(:) = z0_ground(:) * u_star(:) &
         / (1.327*1e-5 * (pb_std/pb(:)) * (temp_air(:)/ZeroCelsius)**(1.81))
    
    kBs_m1(:) = 2.46 * reynolds**(1./4.) - LOG(7.4)

    IF (use_ratio_z0m_z0h) THEN
       ! Do not use exp(kBs_m1(:) but and use ratio_z0m_z0h instead
       ! Same as rough_dyn=F except veget_max(:,1) is used instead of tot_bare_soil
       dragh(:) = veget_max(:,1) * (ct_karman/LOG(ztmp(:)/z0_ground(:)))*(ct_karman/LOG(ztmp(:)/(z0_ground(:)/ratio_z0m_z0h(1)) ))
    ELSE
       dragh(:) = veget_max(:,1) * (ct_karman/LOG(ztmp(:)/z0_ground(:)))*(ct_karman/LOG(ztmp(:)/(z0_ground(:)/exp(kBs_m1(:))) ))
    END IF
    
    ! Fraction of bare soil
    sumveg(:) = veget_max(:,1)

    ! Set average vegetation height to zero
    ave_height(:) = zero

    !+++HACK+++
    ! Use higher threshold for small lai in the calculation of kB_m1 in condveg
    IF (hack_kb_m1) THEN
       threshold=1.E-2
    ELSE
       threshold=min_sechiba
    END IF
    !+++HACK+++
    
    !! 2. Calculate the mean roughness height 
    ! height is an input variable write it to a local variable
    ! such that the value can be changed where needed
    IF(use_height_dom)THEN
       loc_height(:,:) = height_dom(:,:)
    ELSE
       loc_height(:,:) = height(:,:)
    ENDIF

    ! Calculate the mean roughness height of
    ! vegetative PFTs over the grid cell
    DO jv = 2, nvm
     DO ji = 1, kjpindex
      IF (veget_max(ji,jv) .GT. zero) THEN

          IF (lai(ji,jv).GT.min_sechiba .AND. height(ji,jv).LT.min_sechiba) THEN

             ! This is possible when the trees are under coppicing. Coppicing removes 
             ! the aboveground biomass but retains the belowground components. LAI
             ! is calculated making use of the sum of above and belowground biomass
             ! height only makes use of the aboverground biomass. Following 
             ! coppincing we thus have an LAI but no height. I propose a quick fix
             ! by prescribing a height of 10 cm (= the stumps of the trees)
             loc_height(ji,jv) = 0.1

          END IF

          ! Calculate the average roughness over the grid cell:
          ! The unitless drag coefficient is per vegetative PFT
          ! calculated by use of the von Karman constant, the height 
          ! of the first layer and the roughness. The roughness
          ! is calculated as the vegetation height  per PFT 
          ! multiplied by the roughness  parameter 'z0_over_height= 1/16'. 
          ! If this scaled value is lower than 0.01 then the value for 
          ! the roughness of bare soil (0.01) is used. 
          ! The sum over all PFTs gives the average roughness 
          ! per grid cell for the vegetative PFTs.
          eta(ji) = c1 - c2 * exp(-c3 * Cdrag_foliage * lai(ji,jv))

          z0m_pft(ji) = (loc_height(ji,jv)*(1-height_displacement)*(exp(-ct_karman/eta(ji))-exp(-ct_karman/(c1-c2)))) &
               + z0_ground(ji)

          dragm(ji) = dragm(ji) + veget_max(ji,jv) *(ct_karman/LOG(ztmp(ji)/z0m_pft(ji)))**2
   
          fc(ji) = veget(ji,jv)/veget_max(ji,jv)
          fs(ji) = 1. - fc(ji)

          eta_ec(ji) = ( Cdrag_foliage * lai(ji,jv)) / (2 * eta(ji)*eta(ji))
          wind(ji) = SQRT(u(ji)*u(ji)+v(ji)*v(ji))
          u_star(ji)= ct_karman * MAX(min_wind,wind(ji)) / LOG((zlev(ji)+(loc_height(ji,jv)*(1-height_displacement)))/z0m_pft(ji))
          Reynolds(ji) = z0_ground(ji) * u_star(ji) &
               / (1.327*1e-5 * (pb_std/pb(ji)) * (temp_air(ji)/ZeroCelsius)**(1.81))
                 
          kBs_m1(ji) = 2.46 * reynolds(ji)**(1./4.) - LOG(7.4)
          Ct_star(ji) = Prandtl**(-2./3.) * SQRT(1./Reynolds(ji))

          
          IF (lai(ji,jv) .GT. threshold .AND. loc_height(ji,jv).GT.min_sechiba) THEN
             kB_m1(ji) = (ct_karman * Cdrag_foliage) / (4 * Ct * eta(ji) * (1 - exp(-eta_ec(ji)/2.))) * fc(ji)**2. &
                  + 2*fc(ji)*fs(ji) * (ct_karman * eta(ji) * z0m_pft(ji) / loc_height(ji,jv)) / Ct_star(ji) &
                  + kBs_m1(ji) * fs(ji)**2. 
          ELSE
             kB_m1(ji) = kBs_m1(ji) * fs(ji)**2. 
          END IF
          
          IF (use_ratio_z0m_z0h) THEN
             ! Do not use exp(kBs_m1(ji) but and use ratio_z0m_z0h instead
             z0h_pft(ji) = z0m_pft(ji) / ratio_z0m_z0h(jv)
          ELSE
             z0h_pft(ji) = z0m_pft(ji) / exp(kB_m1(ji))
          END IF
          dragh(ji) = dragh(ji) + veget_max(ji,jv) * (ct_karman/LOG(ztmp(ji)/z0m_pft(ji)))*(ct_karman/LOG(ztmp(ji)/z0h_pft(ji)))
          
          ! Sum of bare soil and fraction vegetated fraction
          sumveg(ji) = sumveg(ji) + veget_max(ji,jv)

          ! Weigh height of vegetation with maximal cover fraction
          ave_height(ji) =ave_height(ji) + veget_max(ji,jv)*loc_height(ji,jv)

       END IF
      END DO
    END DO
    
    !! 3. Calculate the mean roughness height of vegetative PFTs over the grid cell
    
    !  Search for pixels with vegetated part to normalise 
    !  roughness height
    WHERE ( sumveg(:) .GT. min_sechiba ) 
       dragh(:) = dragh(:) / sumveg(:)
       dragm(:) = dragm(:) / sumveg(:)
    ENDWHERE

    ! Calculate fraction of roughness for vegetated part 
    dragh(:) = (un - totfrac_nobio(:)) * dragh(:)
    dragm(:) = (un - totfrac_nobio(:)) * dragm(:)

    DO jv = 1, nnobio ! Loop over # of non-vegative surfaces

       ! Set rougness for ice
       IF ( jv .EQ. iice ) THEN          
         z0_nobio = z0_ice
       ELSE
          WRITE(numout,*) 'jv=',jv
          WRITE(numout,*) 'DO NOT KNOW ROUGHNESS OF THIS SURFACE TYPE'
          CALL ipslerr_p(3,'condveg_z0cdrag_dyn','DO NOT KNOW ROUGHNESS OF THIS SURFACE TYPE','','')
       ENDIF
       
       ! Sum of vegetative roughness length and non-vegetative roughness length
       dragm(:) = dragm(:) + frac_nobio(:,jv) * (ct_karman/LOG(ztmp(:)/z0_nobio))**2
       
       u_star(:)= ct_karman * MAX(min_wind,wind(:)) / LOG(zlev(:)/z0_nobio)
       Reynolds(:) = z0_nobio * u_star(:) &
            / (1.327*1e-5 * (pb_std/pb(:)) * (temp_air(:)/ZeroCelsius)**(1.81))
       
       kBs_m1(:) = 2.46 * reynolds**(1./4.) - LOG(7.4)

       IF (use_ratio_z0m_z0h) THEN
          ! Do not use exp(kBs_m1(:) but and use ratio_z0m_z0h instead
          dragh(:) = dragh(:) + frac_nobio(:,jv) * (ct_karman/LOG(ztmp(:)/z0_nobio)) * &
               (ct_karman/LOG(ztmp(:)/(z0_nobio/ ratio_z0m_z0h(1)) ))
       ELSE
          dragh(:) = dragh(:) + frac_nobio(:,jv) * (ct_karman/LOG(ztmp(:)/z0_nobio)) * &
               (ct_karman/LOG(ztmp(:)/(z0_nobio/ exp(kBs_m1(:))) ))
       END IF
       
    ENDDO ! Loop over # of non-vegative surfaces
    
    !! 4. Calculate the zero plane displacement height and effective roughness length
    !  Take the exponential of the roughness
    z0m(:) = ztmp(:) / EXP(ct_karman/SQRT(dragm(:)))
    z0h(:) = ztmp(:) / EXP((ct_karman**2.)/(dragh(:)*LOG(ztmp(:)/z0m(:))))


    ! Compute the zero plane displacement height which
    ! is an equivalent height for the absorption of momentum
    zhdispl(:) = ave_height(:) * height_displacement

    ! In order to calculate the fluxes we compute what we call the grid effective roughness height.
    ! This is the height over which the roughness acts. It combines the
    ! zero plane displacement height and the vegetation height.
    roughheight(:) = ave_height(:) - zhdispl(:)

  END SUBROUTINE condveg_z0cdrag_dyn


END MODULE condveg
