! =================================================================================================================================
! PROGRAM       : forcesoil
!
! CONTACT	: orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      	: IPSL (2006). This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        This subroutine runs the soilcarbon submodel using specific initial conditions 
!! and driving variables in order to obtain soil carbon stocks closed to the steady-state values 
!! quicker than when using the ''full'' ORCHIDEE.  
!!	
!!\n DESCRIPTION: None
!! This subroutine computes the soil carbon stocks by calling the soilcarbon routine at each time step. \n
!! The aim is to obtain soil carbon stocks closed to the steady-state values and ultimately to create  
!! an updated stomate restart file for the stomate component. The state variables of the subsystem are the clay content 
!! (fixed value) and the soil carbon stocks. Initial conditions for the state variables are read in an  
!! input stomate restart file. Driving variables are Soil carbon input, Water and Temperature stresses on 
!! Organic Matter decomposition. Driving variables are read from a specific forcing file produced by a former run of ORCHIDEE
!! (SECHIBA+STOMATE). \n 
!! The FORCESOIL program first consists in reading a set of input files, allocating variables and 
!! preparing output stomate restart file. \n                                                             
!! Then, a loop over time is performed in which the soilcarbon routine is called at each time step. \n
!! Last, final values of the soil carbon stocks are written into the output stomate restart file. \n
!!
!! Flags related:
!! OK_SOIL_CARBON_DISCRETIZATION yes/no(default)
!! STOMATE_CFORCING netcdf filename(default: stomate_cforcing.nc) 
!!              Only if OK_SOIL_CARBON_DISCRETIZATION = y.
!!              If set to NONE no output file will be created.
!! FORCESOIL_STEP_PER_YEAR 1 to 366
!! FORCESOIL_NB_YEAR  1 to *
!!
!! In order to create a carbon forcing file:
!!   Run an orchidee_ol simulation.
!!     Enable OK_SOIL_CARBON_DISCRETIZATION
!!     Set a filename for the flag STOMATE_CFORCING
!!     If required redefine FORCESOIL_STEP_PER_YEAR and FORCESOIL_NB_YEAR
!!
!!   The STOMATE_CFORCING file will keep the data for the latest period given by FORCESOIL_NB_YEAR.
!!   By default it's 1 year.
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S) : None    
!!
!! FLOWCHART    : None
!!
!! SVN		:
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_driver/forcesoil.f90 $ 
!! $Date: 2023-07-04 15:51:03 +0200 (mar. 04 juil. 2023) $
!! $Revision: 8070 $
!! \n
!_ =================================================================================================================================

PROGRAM forcesoil
 
  USE netcdf
  !-
  USE utils
  USE defprec
  USE constantes
  USE constantes_var
  USE constantes_mtc
  USE constantes_soil
  USE pft_parameters 
  USE stomate_data
  USE ioipsl_para
  USE mod_orchidee_para
  USE stomate_soil_carbon_discretization
  USE stomate_io_soil_carbon_discretization
  USE constantes_soil_var
  USE stomate
  USE vertical_soil , ONLY: vertical_soil_init
  USE grid , ONLY : grid_init, regular_lonlat
  USE function_library, ONLY: get_printlev

#ifdef CPP_PARA
  USE mpi
#endif
  !-
  IMPLICIT NONE
  !-
  !-
  CHARACTER(LEN=80)                          :: sto_restname_in,sto_restname_out
  INTEGER(i_std)                             :: iim,jjm                !! Indices (unitless)

  INTEGER(i_std),PARAMETER                   :: llm = 1                !! Vertical Layers (requested by restini routine) (unitless)
  INTEGER(i_std)                             :: kjpindex               !! Domain size (unitless)

  INTEGER(i_std)                             :: itau_dep,itau_len      !! Time step read in the restart file (?) 
                                                                       !! and number of time steps of the simulation (unitless) 
  CHARACTER(LEN=30)                          :: time_str               !! Length of the simulation (year)
  REAL(r_std)                                :: dt_files               !! time step between two successive itaus (?) 
                                                                       !! (requested by restini routine) (seconds)
  REAL(r_std)                                :: date0                  !! Time at which itau = 0 (requested by restini routine) (?)
  INTEGER(i_std)                             :: rest_id_sto            !! ID of the input restart file (unitless)
  CHARACTER(LEN=20), SAVE                    :: thecalendar = 'noleap' !! Type of calendar defined in the input restart file 
                                                                       !! (unitless)
  !-
  CHARACTER(LEN=100)                         :: Cforcing_name          !! Name of the forcing file (unitless)
  INTEGER                                    :: Cforcing_id            !! ID of the forcing file (unitless)
  INTEGER                                    :: v_id                   !! ID of the variable 'Index' stored in the forcing file 
                                                                       !! (unitless)
  REAL(r_std)                                :: dt_forcesoil           !! Time step at which soilcarbon routine is called (days) 
  INTEGER                                    :: nparan                 !! Number of values stored per year in the forcing file 
                                                                       !! (unitless)
  INTEGER                                    :: nbyear
  INTEGER(i_std),DIMENSION(:),ALLOCATABLE    :: indices                !! Grid Point Index used per processor (unitless)
  INTEGER(i_std),DIMENSION(:),ALLOCATABLE    :: indices_g              !! Grid Point Index for all processor (unitless)
  REAL(r_std),DIMENSION(:),ALLOCATABLE       :: x_indices_g            !! Grid Point Index for all processor (unitless)
  REAL(r_std),DIMENSION(:,:),ALLOCATABLE     :: lon, lat               !! Longitude and Latitude of each grid point defined 
                                                                       !! in lat/lon (2D) (degrees)
  REAL(r_std),DIMENSION(llm)                 :: lev                    !! Number of level (requested by restini routine) (unitless)


  INTEGER                                    :: i,m,iatt,iv,iyear      !! counters (unitless)
                                                                       
  CHARACTER(LEN=100)                          :: var_name              
  CHARACTER(LEN=8000)                        :: taboo_vars             !! string used for storing the name of the variables 
                                                                       !! of the stomate restart file that are not automatically 
                                                                       !! duplicated from input to output restart file (unitless)
  REAL(r_std),DIMENSION(1)                   :: xtmp                   !! scalar read/written in restget/restput routines (unitless)
  INTEGER(i_std),PARAMETER                   :: nbvarmax=1000           !! maximum # of variables assumed in the stomate restart file 
                                                                       !! (unitless)
  INTEGER(i_std)                             :: nbvar                  !! # of variables effectively present 
                                                                       !! in the stomate restart file (unitless)
  CHARACTER(LEN=1000),DIMENSION(nbvarmax)      :: varnames              !! list of the names of the variables stored 
                                                                       !! in the stomate restart file (unitless)
  INTEGER(i_std)                             :: varnbdim               !! # of dimensions of a given variable 
                                                                       !! of the stomate restart file
  INTEGER(i_std),PARAMETER                   :: varnbdim_max=20        !! maximal # of dimensions assumed for any variable 
                                                                       !! of the stomate restart file 
  INTEGER,DIMENSION(varnbdim_max)            :: vardims                !! length of each dimension of a given variable 
                                                                       !! of the stomate restart file
  LOGICAL                                    :: l1d                    !! boolean : TRUE if all dimensions of a given variable 
                                                                       !! of the stomate restart file are of length 1 (ie scalar) 
                                                                       !! (unitless)
  REAL(r_std)                                :: x_tmp                  !! temporary variable used to store return value 
  INTEGER(i_std)                               :: orch_vardims         !! Orchidee dimensions (different to IOIPSL -exclude time dim-)
                                                                       !! from nf90_get_att (unitless)
  CHARACTER(LEN=10)  :: part_str                                       !! string suffix indicating the index of a PFT 
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)  :: clay_g                 !! clay fraction (nbpglo) (unitless)
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)  :: depth_organic_soil_g   !! Depth to organic soil (\f$m\f$)
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE   :: control_temp_g         !! Temperature control (nbp_glo,above/below,time) on OM decomposition 
                                                                       !! (unitless)
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE   :: control_moist_g        !! Moisture control (nbp_glo,abo/below,time) on OM decomposition 
                                                                       !! ?? Should be defined per PFT as well (unitless)
  REAL(r_std),DIMENSION(:,:,:,:),ALLOCATABLE :: som_g                  !! Soil organic matter stocks (nbp_glo,ncarb,nvm,nelements) (\f$gC m^{-2}\f$)
                                                                       
  REAL(r_std),ALLOCATABLE :: clay(:)                                   !! clay fraction (nbp_loc) (unitless)
  REAL(r_std),ALLOCATABLE :: depth_organic_soil(:)                     !! Depth to organic soil (\f$m\f$)
  REAL(r_std),ALLOCATABLE :: som_input(:,:,:,:,:)                      !! soil organic matter input (nbp_loc,ncarb,nvm,nelements,time) 
                                                                       !! (\f$gC m^{-2} dt_forcesoil^{-1}\f$) 
  REAL(r_std),ALLOCATABLE :: control_temp(:,:,:)                       !! Temperature control (nbp_loc,above/below,time) on OM decomposition 
                                                                       !! (unitless)
  REAL(r_std),ALLOCATABLE :: control_moist(:,:,:)                      !! Moisture control (nbp_loc,abo/below,time) on OM decomposition 
                                                                       !! ?? Should be defined per PFT as well (unitless)
  REAL(r_std),ALLOCATABLE :: som(:,:,:,:)                              !! Soil organic matter stocks (nbp_loc,ncarb,nvm,nelements) (\f$gC m^{-2}\f$)
  REAL(r_std),DIMENSION(:,:),ALLOCATABLE     :: resp_hetero_soil       !! Heterotrophic respiration (\f$gC m^{-2} dt_forcesoil^{-1}\f$) 
                                                                       !! (requested by soilcarbon routine but not used here) 
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE   :: n_mineralisation       !! net nitrogen mineralisation of decomposing SOM 
  INTEGER(i_std)                             :: printlev_loc           !! Local write level                                                                     
  INTEGER(i_std)                             :: ier,iret               !! Used for hangling errors 
                                                                       
  CHARACTER(LEN=50) :: temp_name                                       
  CHARACTER(LEN=100) :: msg3                                       
  LOGICAL :: debug                                                     !! boolean used for printing messages
  LOGICAL :: l_error                                                   !! boolean for memory allocation
  ! allocateable arrays needed for permafrost carbon 
  REAL(r_std),DIMENSION(:,:,:,:),ALLOCATABLE  :: deepSOM_a_g           !! active organic matter concentration
  REAL(r_std),DIMENSION(:,:,:,:),ALLOCATABLE  :: deepSOM_s_g           !! slow organic matter concentration  
  REAL(r_std),DIMENSION(:,:,:,:),ALLOCATABLE  :: deepSOM_p_g           !! passive matter concentration
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE  :: O2_soil_g               !! oxygen in the soil
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE  :: CH4_soil_g              !! methane in the soil
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE  :: O2_snow_g               !! oxygen in the snow
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE  :: CH4_snow_g              !! methane in the snow 
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE :: snowdz_g                 !! snow depth at each layer
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE :: snowrho_g                !! snow density at each layer
  REAL(r_std),DIMENSION(:,:,:,:),ALLOCATABLE  :: deepSOM_a             !! active organic matter concentration
  REAL(r_std),DIMENSION(:,:,:,:),ALLOCATABLE  :: deepSOM_s             !! slow organic matter concentration
  REAL(r_std),DIMENSION(:,:,:,:),ALLOCATABLE  :: deepSOM_p             !! passive organic matter concentration
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE  :: O2_soil                 !! oxygen in the soil
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE  :: CH4_soil                !! methane in the soil
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE  :: O2_snow                 !! oxygen in the snow
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE  :: CH4_snow                !! methane in the snow 
  REAL(r_std),DIMENSION(:,:),ALLOCATABLE  :: pb                        !! surface pressure
  REAL(r_std),DIMENSION(:,:),ALLOCATABLE  :: snow                      !! snow mass
  REAL(r_std),DIMENSION(:,:,:,:),ALLOCATABLE  :: tprof                 !! deep soil temperature profile
  REAL(r_std),DIMENSION(:,:,:,:),ALLOCATABLE  :: fbact                 !! factor for soil organic matter decomposition
  REAL(r_std),DIMENSION(:,:,:,:),ALLOCATABLE  :: hslong                !! deep soil humidity
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE  :: veget_max               !! maximum vegetation fraction
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE  :: rprof                   !! PFT rooting depth
  REAL(r_std),DIMENSION(:,:),ALLOCATABLE  :: tsurf                     !! surface temperature
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE :: snowrho                  !! snow density
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE :: snowdz                   !! snow depth
  REAL(r_std),DIMENSION(:,:),ALLOCATABLE      :: lalo                  !! Geogr. coordinates (latitude,longitude) (degrees)    
  REAL(r_std), DIMENSION(:,:,:), ALLOCATABLE    :: heat_Zimov                 !! heating associated with decomposition  [W/m**3 soil]
  REAL(R_STD), DIMENSION(:), ALLOCATABLE      :: sfluxCH4_deep, sfluxCO2_deep !! [g / m**2]
  REAL(R_STD), DIMENSION(:,:), ALLOCATABLE      :: altmax                     !! active layer thickness (m)
  REAL(R_STD), DIMENSION(:,:), ALLOCATABLE      :: altmax_g                   !! global active layer thickness (m)
  REAL(r_std), DIMENSION(:,:,:,:),  ALLOCATABLE  :: som_surf_g
  REAL(r_std), DIMENSION(:,:,:,:),  ALLOCATABLE  :: som_surf                 !! vertically-integrated (diagnostic) soil organic matter pools: active, slow, or passive, (gC/(m**2 of ground))
  REAL(R_STD), ALLOCATABLE, DIMENSION(:,:)      :: fixed_cryoturbation_depth  !! depth to hold cryoturbation to for fixed runs
  REAL(R_STD), ALLOCATABLE, DIMENSION(:,:,:,:)  :: CN_target  !! C to N ratio of SOM flux from one pool to another (gN m-2 dt-1) 
  LOGICAL, SAVE                             :: satsoil = .FALSE.
  LOGICAL                                   :: reset_soilc = .false.

  INTEGER(i_std)                            :: start_2d(2), count_2d(2) 
  INTEGER(i_std)                            :: start_5d(5), count_5d(5), start_4d(4), count_4d(4), start_3d(3), count_3d(3)
!_ =================================================================================================================================
 
  CALL Init_orchidee_para
  CALL init_timer

! Set specific write level to forcesoil using PRINTLEV_forcesoil=[0-4] in run.def. 
! The global printlev is used as default value. 
  printlev_loc=get_printlev('forcesoil')

!-
! Configure the number of PFTS 
!-
  ok_soil_carbon_discretization=.FALSE. 
  CALL getin_p('OK_SOIL_CARBON_DISCRETIZATION',ok_soil_carbon_discretization)

  IF (.NOT. ok_soil_carbon_discretization) THEN
     CALL ipslerr_p(3,'forcesoil','OK_SOIL_CARBON_DISCRETIZATION must be enabled','But found:','True')
  ENDIF

  ! 1. Read the number of PFTs
  !
  !Config Key   = NVM
  !Config Desc  = number of PFTs  
  !Config If    = OK_SECHIBA or OK_STOMATE
  !Config Def   = 13
  !Config Help  = The number of vegetation types define by the user
  !Config Units = [-]
  CALL getin_p('NVM',nvm)

  ! 2. Allocation
  ALLOCATE(pft_to_mtc(nvm),stat=ier)
  IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'pft_to_mtc : error in memory allocation', '', '')

  ! 3. Initialisation of the correspondance table
  pft_to_mtc(:) = undef_int
  
  ! 4.Reading of the conrrespondance table in the .def file
  !
  !Config Key   = PFT_TO_MTC
  !Config Desc  = correspondance array linking a PFT to MTC
  !Config if    = OK_SECHIBA or OK_STOMATE
  !Config Def   = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13
  !Config Help  =
  !Config Units = [-]
  CALL getin_p('PFT_TO_MTC',pft_to_mtc)

  ! 4.1 if nothing is found, we use the standard configuration
  IF(nvm <= nvmc ) THEN
     IF(pft_to_mtc(1) == undef_int) THEN
        WRITE(numout,*) 'Note to the user : we will use ORCHIDEE to its standard configuration'
        pft_to_mtc(:) = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 /)
     ENDIF
  ELSE   
     IF(pft_to_mtc(1) == undef_int) THEN
        WRITE(numout,*)' The array PFT_TO_MTC is empty : we stop'
     ENDIF
  ENDIF
  
  ! 4.2 What happened if pft_to_mtc(j) > nvmc (if the mtc doesn't exist)?
  DO i = 1, nvm
     IF(pft_to_mtc(i) > nvmc) THEN
        CALL ipslerr_p(3, 'forcesoil', 'the MTC you chose doesnt exist', 'we stop reading pft_to_mtc', '')
     ENDIF
  ENDDO
  
  ! 4.3 Check if pft_to_mtc(1) = 1 
  IF(pft_to_mtc(1) /= 1) THEN
     CALL ipslerr_p(3, 'forcesoil', 'the first pft has to be the bare soil', 'we stop reading next values of pft_to_mtc', '')
  ENDIF

  DO i = 2,nvm
     IF(pft_to_mtc(i) == 1) THEN
        CALL ipslerr_p(3, 'forcesoil', 'only pft_to_mtc(1) has to be the bare soil', 'we stop reading next values of pft_to_mtc', '')
     ENDIF
  ENDDO
  
  ! 5. Allocate and initialize natural and is_c4
  
  ! 5.1 Memory allocation
  ALLOCATE(natural(nvm),stat=ier)
  IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'natural : error in memory allocation', '', '')

  ALLOCATE(is_c4(nvm),stat=ier)
  IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'is_c4 : error in memory allocation', '', '')

!  ALLOCATE(permafrost_veg_exists(nvm),stat=ier)
!  IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'permafrost_veg_exists : error in memory allocation', '', '')

  ! 5.2 Initialisation
  DO i = 1, nvm
     natural(i) = natural_mtc(pft_to_mtc(i))
     is_c4(i) = is_c4_mtc(pft_to_mtc(i))
  ENDDO

!  DO i = 1, nvm
!     permafrost_veg_exists(i) = permafrost_veg_exists_mtc(pft_to_mtc(i))
!  ENDDO
  !!- 
  !! 1. Initialisation stage
  !! Reading a set of input files, allocating variables and preparing output restart file.     
  !!-
  ! Define restart file name
  ! for reading initial conditions (sto_restname_in var) and for writting final conditions (sto_restname_out var). 
  ! User values are used if present in the .def file.
  ! If not present, default values (stomate_start.nc and stomate_rest_out.c) are used.
  !-
  IF (is_root_prc) THEN
     sto_restname_in = 'stomate_start.nc'
     CALL getin ('STOMATE_RESTART_FILEIN',sto_restname_in)
     WRITE(numout,*) 'STOMATE INPUT RESTART_FILE: ',TRIM(sto_restname_in)
     sto_restname_out = 'stomate_rest_out.nc'
     CALL getin ('STOMATE_RESTART_FILEOUT',sto_restname_out)
     WRITE(numout,*) 'STOMATE OUTPUT RESTART_FILE: ',TRIM(sto_restname_out)
     CALL getin ('satsoil', satsoil)
     !-
     ! Open the input file and Get some Dimension and Attributes ID's 
     !-
     CALL nccheck( NF90_OPEN (sto_restname_in, NF90_NOWRITE, rest_id_sto))
     CALL nccheck( NF90_INQUIRE_DIMENSION (rest_id_sto,1,len=iim_g))
     CALL nccheck( NF90_INQUIRE_DIMENSION (rest_id_sto,2,len=jjm_g))
     CALL nccheck( NF90_INQ_VARID (rest_id_sto, "time", iv))
     CALL nccheck( NF90_GET_ATT (rest_id_sto, iv, 'calendar',thecalendar))
     CALL nccheck( NF90_CLOSE (rest_id_sto))
     i=INDEX(thecalendar,ACHAR(0))
     IF ( i > 0 ) THEN
        thecalendar(i:20)=' '
     ENDIF
     !-
     ! Allocate longitudes and latitudes
     !-
     ALLOCATE (lon(iim_g,jjm_g), stat=ier)
     IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'lon : error in memory allocation', '', '')
     ALLOCATE (lat(iim_g,jjm_g), stat=ier)
     IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'lat : error in memory allocation', '', '')
     lon(:,:) = zero
     lat(:,:) = zero
     lev(1)   = zero
     !-
     CALL restini &
          & (sto_restname_in, iim_g, jjm_g, lon, lat, llm, lev, &
          &  sto_restname_out, itau_dep, date0, dt_files, rest_id_sto, &
          &  use_compression=nc_restart_compression )
  ENDIF

  CALL bcast(date0)
  CALL bcast(thecalendar)
  WRITE(numout,*) "calendar = ",thecalendar
  !-
  ! calendar
  !-
  CALL ioconf_calendar (thecalendar)
  CALL ioget_calendar  (one_year,one_day)
  CALL ioconf_startdate(date0)
  !
  !-
  ! define forcing file's name (Cforcing_name var)
  ! User value is used if present in the .def file
  ! If not, default (NONE) is used
  !-
  Cforcing_name = 'stomate_cforcing.nc'
  CALL getin ('STOMATE_CFORCING_NAME',Cforcing_name)
  !
  IF (TRIM(Cforcing_name) .EQ. 'NONE') THEN
    CALL ipslerr_p(3,'forcesoil','STOMATE_CFORCING_NAME key must be defined with a filename',&
            'But found:',Cforcing_name)
  ENDIF
  !-
  ! Initailize ngrnd
  !-
  CALL vertical_soil_init()
  !-
  ! Disable XIOS outputs for spinup
  !-
  soilc_isspinup = .TRUE.
  !
  !! For master process only
  !
  IF (is_root_prc) THEN
     !-
     ! Open FORCESOIL's forcing file to read some basic info (dimensions, variable ID's)
     ! and allocate variables.
     !-
#ifdef CPP_PARA
     CALL nccheck( NF90_OPEN (TRIM(Cforcing_name),IOR(NF90_NOWRITE, NF90_MPIIO),Cforcing_id, &
         & comm = MPI_COMM_ORCH, info = MPI_INFO_NULL ))
#else
     CALL nccheck( NF90_OPEN (TRIM(Cforcing_name),NF90_NOWRITE,Cforcing_id))
#endif
     !-
     ! Total Domain size is stored in nbp_glo variable
     !-
     CALL nccheck( NF90_GET_ATT (Cforcing_id,NF90_GLOBAL,'kjpindex',x_tmp))
     nbp_glo = NINT(x_tmp)
     !-
     ! Number of values stored per year in the forcing file is stored in nparan var.
     !-
     CALL nccheck( NF90_GET_ATT (Cforcing_id,NF90_GLOBAL,'nparan',x_tmp))
     nparan = NINT(x_tmp)
     CALL nccheck( NF90_GET_ATT (Cforcing_id,NF90_GLOBAL,'nbyear',x_tmp))
     nbyear = NINT(x_tmp)
     !-
     ALLOCATE (indices_g(nbp_glo), stat=ier)
     IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'indices_g : error in memory allocation', '', '')
     ALLOCATE (clay_g(nbp_glo), stat=ier)
     IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'clay_g : error in memory allocation', '', '')
     ALLOCATE (depth_organic_soil_g(nbp_glo), stat=ier)
     IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'depth_organic_soil_g : error in memory allocation', '', '')
     !-
     ALLOCATE (x_indices_g(nbp_glo),stat=ier)
     IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'x_indices_g : error in memory allocation', '', '')
     CALL nccheck( NF90_INQ_VARID (Cforcing_id,'index',v_id))
     CALL nccheck( NF90_GET_VAR   (Cforcing_id,v_id,x_indices_g))
     indices_g(:) = NINT(x_indices_g(:))
     WRITE(numout,*) mpi_rank,"indices globaux : ",indices_g
     DEALLOCATE (x_indices_g)
     !-
     CALL nccheck( NF90_INQ_VARID (Cforcing_id,'clay',v_id))
     CALL nccheck( NF90_GET_VAR   (Cforcing_id,v_id,clay_g))
     CALL nccheck( NF90_INQ_VARID (Cforcing_id,'depth_organic_soil',v_id))
     CALL nccheck( NF90_GET_VAR   (Cforcing_id,v_id,depth_organic_soil_g))
     !-
     ! time step of forcesoil program (in days)
     !-
     dt_forcesoil = one_year / FLOAT(nparan)
     WRITE(numout,*) 'time step (d): ',dt_forcesoil
     WRITE(numout,*) 'nparan: ',nparan
     WRITE(numout,*) 'nbyear: ',nbyear    
     !-
     ! read and write the variables in the output restart file we do not modify within the Forcesoil program
     ! ie all variables stored in the input restart file except those stored in taboo_vars
     !-
     !-
     taboo_vars = '$nav_lon$ $nav_lat$ $nav_lev$ $time$ $time_steps$ '// &
     &            '$day_counter$ $dt_days$ $date$ $deepSOM_a$ $deepSOM_s$ '// &
     &            '$deepSOM_p$ $O2_soil$ $CH4_soil$ $O2_snow$ $CH4_snow$ '// &
     &            '$altmax$ '  
     !-
     !-
     CALL ioget_vname(rest_id_sto, nbvar, varnames)
     !-
     ! read and write some special variables (1D or variables that we need)
     !-
     CALL restransfer (rest_id_sto, 'day_counter', itau_dep)
     !-
     CALL restransfer (rest_id_sto, 'dt_days', itau_dep)
     !-
     CALL restransfer (rest_id_sto, 'date', itau_dep)
     !-
     DO iv=1,nbvar
        !-- check if the variable is to be written here
        IF (INDEX(taboo_vars,'$'//TRIM(varnames(iv))//'$') == 0 ) THEN

           CALL restransfer(rest_id_sto, TRIM(varnames(iv)), itau_dep, nbp_glo, indices_g)
         
        ENDIF ! INDEX(taboo_vars,'$'//TRIM(varnames(iv))//'$') == 0 )
     ENDDO
     ! Length of the run (in Years)
     ! User value is used if present in the .def file
     ! If not, default value (10000 Years) is used
     !-
     WRITE(time_str,'(a)') '10000Y'
     CALL getin('TIME_LENGTH', time_str)
     write(numout,*) 'Number of years for som spinup : ',time_str
     ! transform into itau
     CALL tlen2itau(time_str, dt_forcesoil*one_day, date0, itau_len)
     write(numout,*) 'Number of time steps to do: ',itau_len

     ! read soil organic matter stocks values stored in the input restart file
     !-
     !-
     ! Permafrost soil organic matter
     !-
      ALLOCATE(som_g(nbp_glo,ncarb,nvm,nelements), stat=ier)
      IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'som_g : error in memory allocation', '', '')
      som_g(:,:,:,:) = 0.
      ALLOCATE(som_surf_g(nbp_glo,ncarb,nvm,nelements), stat=ier)
      IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'som_surf_g : error in memory allocation', '', '')
      som_surf_g(:,:,:,:) = 0.
      ALLOCATE(deepSOM_a_g(nbp_glo,ngrnd,nvm,nelements), stat=ier)
      IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'deepSOM_a_g : error in memory allocation', '', '')
      ALLOCATE(deepSOM_s_g(nbp_glo,ngrnd,nvm,nelements), stat=ier)
      IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'deepSOM_s_g : error in memory allocation', '', '')
      ALLOCATE(deepSOM_p_g(nbp_glo,ngrnd,nvm,nelements), stat=ier)
      IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'deepSOM_p_g : error in memory allocation', '', '')
      ALLOCATE(O2_soil_g(nbp_glo,ngrnd,nvm), stat=ier)
      IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'O2_soil_g : error in memory allocation', '', '')
      ALLOCATE(CH4_soil_g(nbp_glo,ngrnd,nvm), stat=ier)
      IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'CH4_soil_g : error in memory allocation', '', '')
      ALLOCATE(O2_snow_g(nbp_glo,nsnow,nvm), stat=ier)
      IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'O2_snow_g : error in memory allocation', '', '')
      ALLOCATE(CH4_snow_g(nbp_glo,nsnow,nvm), stat=ier)
      IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'CH4_snow_g : error in memory allocation', '', '')
      ALLOCATE(altmax_g(nbp_glo,nvm), stat=ier)
      IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'altmax_g : error in memory allocation', '', '')

      deepSOM_a_g(:,:,:,:) = val_exp
      CALL restget (rest_id_sto, 'deepSOM_a', nbp_glo, ngrnd, nvm, nelements, itau_dep, &
            &               .TRUE., deepSOM_a_g, 'gather', nbp_glo, indices_g)
      IF (ALL(deepSOM_a_g == val_exp)) deepSOM_a_g = zero

      deepSOM_s_g(:,:,:,:) = val_exp
      CALL restget (rest_id_sto, 'deepSOM_s', nbp_glo, ngrnd, nvm, nelements, itau_dep, &
            &               .TRUE., deepSOM_s_g, 'gather', nbp_glo, indices_g)
      IF (ALL(deepSOM_s_g == val_exp)) deepSOM_s_g = zero
     
      deepSOM_p_g(:,:,:,:) = val_exp
      CALL restget (rest_id_sto, 'deepSOM_p', nbp_glo, ngrnd, nvm, nelements, itau_dep, &
           &               .TRUE., deepSOM_p_g, 'gather', nbp_glo, indices_g)
      IF (ALL(deepSOM_p_g == val_exp)) deepSOM_p_g = zero

      var_name= 'altmax'
      altmax_g(:,:) = val_exp
      CALL restget (rest_id_sto, var_name, nbp_glo, nvm, 1, itau_dep, .TRUE., altmax_g, "gather", nbp_glo, indices_g)
      IF ( ALL( altmax_g(:,:) .EQ. val_exp ) ) THEN
          CALL ipslerr(3, 'forcesoil', 'altmax is not found in stomate restart file', '', '')
      END IF

      CALL getin('reset_soilc', reset_soilc)
      IF (reset_soilc) THEN
          CALL ipslerr(1, 'forcesoil', 'deepSOM_a, deepSOM_s and deeSOM_p',  & 
                       'are ignored and set to zero value due to', 'reset_soilc option')
         deepSOM_a_g(:,:,:,:) = zero
         deepSOM_s_g(:,:,:,:) = zero
         deepSOM_p_g(:,:,:,:) = zero
      ENDIF
     
      O2_soil_g(:,:,:) = val_exp
      CALL restget (rest_id_sto, 'O2_soil', nbp_glo, ngrnd, nvm, itau_dep, &
           &               .TRUE., O2_soil_g, 'gather', nbp_glo, indices_g)
      IF (ALL(O2_soil_g == val_exp)) O2_soil_g = O2_init_conc

      CH4_soil_g(:,:,:) = val_exp
      CALL restget (rest_id_sto, 'CH4_soil', nbp_glo, ngrnd, nvm, itau_dep, &
          &               .TRUE., CH4_soil_g, 'gather', nbp_glo, indices_g)
      IF (ALL(CH4_soil_g == val_exp)) CH4_soil_g =  CH4_init_conc

      O2_snow_g(:,:,:) = val_exp
      CALL restget (rest_id_sto, 'O2_snow', nbp_glo, nsnow, nvm, itau_dep, &
           &               .TRUE., O2_snow_g, 'gather', nbp_glo, indices_g)
      IF (ALL(O2_snow_g == val_exp)) O2_snow_g =  O2_init_conc

      CH4_snow_g(:,:,:) = val_exp
      CALL restget (rest_id_sto, 'CH4_snow', nbp_glo, nsnow, nvm, itau_dep, &
          &               .TRUE., CH4_snow_g, 'gather', nbp_glo, indices_g)
      IF (ALL(CH4_snow_g == val_exp)) CH4_snow_g = CH4_init_conc
  ENDIF ! is_root_prc
  !
  CALL bcast(nbp_glo)
  CALL bcast(iim_g)
  CALL bcast(jjm_g)
  IF (.NOT. ALLOCATED(indices_g)) ALLOCATE (indices_g(nbp_glo))
  CALL bcast(indices_g)
  CALL bcast(nparan)
  CALL bcast(nbyear)
  CALL bcast(dt_forcesoil)
  CALL bcast(itau_dep)
  CALL bcast(itau_len)
  !
  ! we must initialize data_para :
  CALL init_orchidee_data_para_driver(nbp_glo,indices_g)

  kjpindex=nbp_loc
  jjm=jj_nb
  iim=iim_g
  IF (printlev_loc>=3) WRITE(numout,*) "Local grid : ",kjpindex,iim,jjm
  !-
  ! Initialize grid type
  !-
  CALL grid_init ( kjpindex, 4, regular_lonlat, "ForcingGrid" )
  !-
  ! Analytical spinup is set to false
  !
  spinup_analytic = .FALSE.
  !-
  ! read soil organic matter inputs, water and temperature stresses on OM
  ! decomposition 
  ! into the forcing file - We read an average year.
  !-
  !-
  ! Read permafrost-related soil organic matter from Cforcing_id
  !-
  ALLOCATE(pb(kjpindex,nparan*nbyear), stat=ier)
  IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'pb : error in memory allocation', '', '')
  ALLOCATE(snow(kjpindex,nparan*nbyear), stat=ier)
  IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'snow : error in memory allocation', '', '')
  ALLOCATE(tprof(kjpindex,ngrnd,nvm,nparan*nbyear), stat=ier)
  IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'tprof : error in memory allocation', '', '')
  ALLOCATE(fbact(kjpindex,ngrnd,nvm,nparan*nbyear), stat=ier)
  IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'fbact : error in memory allocation', '', '')
  ALLOCATE(hslong(kjpindex,ngrnd,nvm,nparan*nbyear), stat=ier)
  IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'hslong : error in memory allocation', '', '')
  ALLOCATE(veget_max(kjpindex,nvm,nparan*nbyear), stat=ier)
  IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'veget_max : error in memory allocation', '', '')
  ALLOCATE(rprof(kjpindex,nvm,nparan*nbyear), stat=ier)
  IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'rprof : error in memory allocation', '', '')
  ALLOCATE(tsurf(kjpindex,nparan*nbyear), stat=ier)
  IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'tsurf : error in memory allocation', '', '')
  ALLOCATE(lalo(kjpindex,2), stat=ier)
  IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'lalo : error in memory allocation', '', '')
  ALLOCATE(snowdz(kjpindex,nsnow,nparan*nbyear), stat=ier)
  IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'snowdz_ : error in memory allocation', '', '')
  ALLOCATE(snowrho(kjpindex,nsnow,nparan*nbyear), stat=ier)
  IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'snowrho : error in memory allocation', '', '')
  ALLOCATE(som_input(kjpindex,ncarb,nvm,nelements,nparan*nbyear), stat=ier)
  IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'som_input : error in memory allocation', '', '')
  ALLOCATE(CN_target(kjpindex,nvm,ncarb,nparan*nbyear), stat=ier)
  IF (ier /= 0) CALL ipslerr(3, 'forcesoil', 'CN_target : error in memory allocation', 'Error code:', '')
  ALLOCATE(n_mineralisation(kjpindex, nvm, nparan*nbyear), stat=ier)
  IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'n_mineralisation : error in memory allocation', '', '')    
  !-
  CALL stomate_io_soil_carbon_discretization_read(Cforcing_name,  nparan,      nbyear,&
                nbp_mpi_para_begin(mpi_rank),   nbp_mpi_para(mpi_rank),         &
                som_input,               pb,         snow,       tsurf,  &
                tprof,                          fbact,      hslong,     rprof,  &
                lalo,                           snowdz,     snowrho,    veget_max, &
                CN_target, n_mineralisation )
  !---
  !--- Create the index table
  !---
  !--- This job returns a LOCAL kindex.
  !---
  ALLOCATE (indices(kjpindex),stat=ier)
  IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'indices : error in memory allocation', '', '')
  !
  !! scattering to all processes in parallel mode
  !
  CALL scatter(indices_g,indices)
  indices(1:kjpindex)=indices(1:kjpindex)-(jj_begin-1)*iim_g
  IF (printlev_loc>=3) WRITE(numout,*) mpi_rank,"indices locaux = ",indices(1:kjpindex)
  !-
  ! Allocation of the variables for a processor
  !-
  ALLOCATE(clay(kjpindex), stat=ier)
  IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'clay : error in memory allocation', '', '')
  ALLOCATE(depth_organic_soil(kjpindex), stat=ier)
  IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'depth_organic_soil : error in memory allocation', '', '')
  ALLOCATE(som(kjpindex,ncarb,nvm,nelements), stat=ier)
  IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'indices : error in memory allocation', '', '')
  !-
  ALLOCATE(som_surf(kjpindex,ncarb,nvm,nelements), stat=ier)
  IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'som_surf : error in memory allocation', '', '')
  ALLOCATE(deepSOM_a(kjpindex,ngrnd,nvm,nelements), stat=ier)
  IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'deepSOM_a : error in memory allocation', '', '')
  ALLOCATE(deepSOM_s(kjpindex,ngrnd,nvm,nelements), stat=ier)
  IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'deepSOM_s : error in memory allocation', '', '')
  ALLOCATE(deepSOM_p(kjpindex,ngrnd,nvm,nelements), stat=ier)
  IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'deepSOM_p : error in memory allocation', '', '')
  ALLOCATE(O2_soil(kjpindex,ngrnd,nvm), stat=ier)
  IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'O2_soil : error in memory allocation', '', '')
  ALLOCATE(CH4_soil(kjpindex,ngrnd,nvm), stat=ier)
  IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'CH4_soil : error in memory allocation', '', '')
  ALLOCATE(O2_snow(kjpindex,nsnow,nvm), stat=ier)
  IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'O2_snow : error in memory allocation', '', '')
  ALLOCATE(CH4_snow(kjpindex,nsnow,nvm), stat=ier)
  IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'CH4_snow : error in memory allocation', '', '')
  ALLOCATE(altmax(kjpindex,nvm), stat=ier)
  IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'altmax : error in memory allocation', '', '')
  !-
  ! Initialization of the variables for a processor
  !-
  CALL Scatter(clay_g,clay)
  CALL Scatter(depth_organic_soil_g,depth_organic_soil)
  CALL Scatter(som_g,som)
  CALL Scatter(som_surf_g,som_surf)
  CALL Scatter(deepSOM_a_g,deepSOM_a)
  CALL Scatter(deepSOM_s_g,deepSOM_s)
  CALL Scatter(deepSOM_p_g,deepSOM_p)
  CALL Scatter(O2_soil_g,O2_soil)
  CALL Scatter(CH4_soil_g,CH4_soil)
  CALL Scatter(O2_snow_g,O2_snow)
  CALL Scatter(CH4_snow_g,CH4_snow)
  CALL Scatter(altmax_g,altmax)
!-
! Configuration of the parameters
!-
   !-
    ! soilcarbon parameters 
    !-
    !Config Key   = ACTIVE_TO_PASS_CLAY_FRAC
    !Config Desc  = 
    !Config if    = OK_STOMATE 
    !Config Def   = 0.68  
    !Config Help  =
    !Config Units = [-] 
    CALL getin_p('ACTIVE_TO_PASS_CLAY_FRAC',active_to_pass_clay_frac)


    !Config Key   = ACTIVE_TO_PASS_REF_FRAC
    !Config Desc  = Fixed fraction from Active to Passive pool
    !Config if    = OK_STOMATE 
    !Config Def   = 0.003
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('ACTIVE_TO_PASS_REF_FRAC',active_to_pass_ref_frac)
    !
    !Config Key   = SURF_TO_SLOW_REF_FRAC
    !Config Desc  = Fixed fraction from Surface to Slow pool
    !Config if    = OK_STOMATE 
    !Config Def   = 0.4
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('SURF_TO_SLOW_REF_FRAC',surf_to_slow_ref_frac)
    !
    !Config Key   = ACTIVE_TO_CO2_REF_FRAC
    !Config Desc  = Fixed fraction from Active pool to CO2 emission
    !Config if    = OK_STOMATE 
    !Config Def   = 0.85
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('ACTIVE_TO_CO2_REF_FRAC',active_to_co2_ref_frac)
    !
    !Config Key   = SLOW_TO_PASS_REF_FRAC
    !Config Desc  = Fixed fraction from Slow to Passive pool
    !Config if    = OK_STOMATE 
    !Config Def   = 0.003
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('SLOW_TO_PASS_REF_FRAC',slow_to_pass_ref_frac)
    !
    !Config Key   = SLOW_TO_CO2_REF_FRAC
    !Config Desc  = Fixed fraction from Slow pool to CO2 emission
    !Config if    = OK_STOMATE 
    !Config Def   = 0.55
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('SLOW_TO_CO2_REF_FRAC',slow_to_co2_ref_frac)
    !
    !Config Key   = PASS_TO_ACTIVE_REF_FRAC
    !Config Desc  = Fixed fraction from Passive to Active pool
    !Config if    = OK_STOMATE 
    !Config Def   = 0.45
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('PASS_TO_ACTIVE_REF_FRAC',pass_to_active_ref_frac)
    !
    !Config Key   = PASS_TO_SLOW_REF_FRAC
    !Config Desc  = Fixed fraction from Passive to Slow pool
    !Config if    = OK_STOMATE 
    !Config Def   = 0.
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('PASS_TO_SLOW_REF_FRAC',pass_to_slow_ref_frac)
    !
    !Config Key   = ACTIVE_TO_CO2_CLAY_SILT_FRAC
    !Config Desc  = Clay-Silt-dependant fraction from Active pool to CO2 emission
    !Config if    = OK_STOMATE 
    !Config Def   = 0.68
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('ACTIVE_TO_CO2_CLAY_SILT_FRAC',active_to_co2_clay_silt_frac)
    !
    !Config Key   = SLOW_TO_PASS_CLAY_FRAC
    !Config Desc  = Clay-dependant fraction from Slow to Passive pool
    !Config if    = OK_STOMATE 
    !Config Def   = -0.009
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('SLOW_TO_PASS_CLAY_FRAC',slow_to_pass_clay_frac)
    !
    !Config Key   = SOM_TURN_IACTIVE
    !Config Desc  = turnover in active pool
    !Config if    = OK_STOMATE 
    !Config Def   = 7.3
    !Config Help  =
    !Config Units =  [year-1] 
    CALL getin_p('SOM_TURN_IACTIVE',som_turn_iactive)
    !
    !Config Key   = SOM_TURN_ISLOW
    !Config Desc  = turnover in slow pool
    !Config if    = OK_STOMATE 
    !Config Def   = 0.2
    !Config Help  =
    !Config Units = [year-1]
    CALL getin_p('SOM_TURN_ISLOW',som_turn_islow)
    !
    !Config Key   = SOM_TURN_IPASSIVE
    !Config Desc  = turnover in passive pool
    !Config if    = OK_STOMATE 
    !Config Def   = 0.0045
    !Config Help  = 
    !Config Units = [year-1] 
    CALL getin_p('SOM_TURN_IPASSIVE',som_turn_ipassive)
    !
    !Config Key   = SOM_TURN_IACTIVE_CLAY_FRAC
    !Config Desc  = clay-dependant parameter impacting on turnover rate of active pool - Tm parameter of Parton et al. 1993 (-)
    !Config if    = OK_STOMATE 
    !Config Def   = 0.75
    !Config Help  = 
    !Config Units = [-] 
    CALL getin_p('SOM_TURN_IACTIVE_CLAY_FRAC',som_turn_iactive_clay_frac)
    !
    !Config Key   = CN_TARGET_IACTIVE_REF
    !Config Desc  = CN target ratio of active pool for soil min N = 0
    !Config if    = OK_STOMATE 
    !Config Def   = 15.
    !Config Help  = 
    !Config Units = [-] 
    CALL getin_p('CN_TARGET_IACTIVE_REF',CN_target_iactive_ref)
    !
    !Config Key   = CN_TARGET_ISLOW_REF
    !Config Desc  = CN target ratio of slow pool for soil min N = 0
    !Config if    = OK_STOMATE 
    !Config Def   = 20.
    !Config Help  = 
    !Config Units = [-] 
    CALL getin_p('CN_TARGET_ISLOW_REF',CN_target_islow_ref)
    !
    !Config Key   = CN_TARGET_IPASSIVE_REF
    !Config Desc  = CN target ratio of passive pool for soil min N = 0
    !Config if    = OK_STOMATE 
    !Config Def   = 10.
    !Config Help  = 
    !Config Units = [-] 
    CALL getin_p('CN_TARGET_IPASSIVE_REF',CN_target_ipassive_ref)
    !
    !Config Key   = CN_TARGET_IACTIVE_NMIN
    !Config Desc  = CN target ratio change per mineral N unit (g m-2) for active pool 
    !Config if    = OK_STOMATE 
    !Config Def   = -6.
    !Config Help  = 
    !Config Units = [(g m-2)-1] 
    CALL getin_p('CN_TARGET_IACTIVE_NMIN',CN_target_iactive_Nmin)
    !
    !Config Key   = CN_TARGET_ISLOW_NMIN
    !Config Desc  = CN target ratio change per mineral N unit (g m-2) for slow pool 
    !Config if    = OK_STOMATE 
    !Config Def   = -4.
    !Config Help  = 
    !Config Units = [(g m-2)-1] 
    CALL getin_p('CN_TARGET_ISLOW_NMIN',CN_target_islow_Nmin)
    !
    !Config Key   = CN_TARGET_IPASSIVE_NMIN
    !Config Desc  = CN target ratio change per mineral N unit (g m-2) for passive pool 
    !Config if    = OK_STOMATE 
    !Config Def   = -1.5
    !Config Help  = 
    !Config Units = [(g m-2)-1] 
    CALL getin_p('CN_TARGET_IPASSIVE_NMIN',CN_target_ipassive_Nmin)

  !
  !! 2. Computational step
  !! Loop over time - Call of soilcarbon routine at each time step 
  !! Updated soil carbon stocks are stored into carbon variable
  !! We only keep the last value of carbon variable (no time dimension).
  !!-
   IF ( satsoil )  hslong(:,:,:,:) = 1.
   !these variables are only ouputs from deep_carbcycle (thus not necessary for
   !Gather and Scatter)
   ALLOCATE(heat_Zimov(kjpindex,ngrnd,nvm), stat=ier)
   IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'heat_Zimov : error in memory allocation', '', '')
   ALLOCATE(sfluxCH4_deep(kjpindex), stat=ier)
   IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'sfluxCH4_deep : error in memory allocation', '', '')
   ALLOCATE(sfluxCO2_deep(kjpindex), stat=ier)
   IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'sfluxCO2_deep : error in memory allocation', '', '')
   ALLOCATE(fixed_cryoturbation_depth(kjpindex,nvm), stat=ier)
   IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'fixed_cryoturbation_depth : error in memory allocation', '', '')
   ALLOCATE(resp_hetero_soil(kjpindex, nvm), stat=ier)
   IF (ier /= 0) CALL ipslerr_p(3, 'forcesoil', 'resp_hetero_soil : error in memory allocation', '', '')
   iatt = 0
   iyear=1
   DO i=1,itau_len
      iatt = iatt+1
      IF (iatt > nparan*nbyear) THEN
            IF (printlev>=3) WRITE(numout,*) iyear
            iatt = 1
            iyear=iyear+1
      ENDIF
      WRITE(numout, *) "Forcesoil:: deep_carbcycle, iyear=", iyear

      CALL stomate_soil_carbon_discretization_deep_somcycle(kjpindex, indices, iatt, dt_forcesoil*one_day, lalo, clay, &
         tsurf(:,iatt), tprof(:,:,:,iatt), hslong(:,:,:,iatt), snow(:,iatt), heat_Zimov, pb(:,iatt), &
         sfluxCH4_deep, sfluxCO2_deep,  &
         deepSOM_a, deepSOM_s, deepSOM_p, O2_soil, CH4_soil, O2_snow, CH4_snow, &
         depth_organic_soil, som_input(:,:,:,:,iatt), &
         veget_max(:,:,iatt), rprof(:,:,iatt), altmax,  som, som_surf, resp_hetero_soil, &
         fbact(:,:,:,iatt), CN_target(:,:,:,iatt), fixed_cryoturbation_depth, snowdz(:,:,iatt), snowrho(:,:,iatt),n_mineralisation(:,:,iatt))
   ENDDO
  !!-
  !! 3. write new soil organic matter stocks into the ouput restart file
  !!-
  CALL restput_p (rest_id_sto, 'som', nbp_glo, ncarb , nvm, nelements, itau_dep, &
         &     som, 'scatter', nbp_glo, indices_g)
  CALL restput_p (rest_id_sto, 'deepSOM_a', nbp_glo, ngrnd, nvm, nelements, itau_dep, &
        &               deepSOM_a, 'scatter', nbp_glo, indices_g)
  CALL restput_p (rest_id_sto, 'deepSOM_s', nbp_glo, ngrnd, nvm, nelements, itau_dep, &
        &               deepSOM_s, 'scatter', nbp_glo, indices_g)
  CALL restput_p (rest_id_sto, 'deepSOM_p', nbp_glo, ngrnd, nvm, nelements, itau_dep, &
        &               deepSOM_p, 'scatter', nbp_glo, indices_g)
  CALL restput_p (rest_id_sto, 'O2_soil', nbp_glo, ngrnd, nvm, itau_dep, &
        &               O2_soil, 'scatter', nbp_glo, indices_g)
  CALL restput_p (rest_id_sto, 'CH4_soil', nbp_glo, ngrnd, nvm, itau_dep, &
        &               CH4_soil, 'scatter', nbp_glo, indices_g)
  CALL restput_p (rest_id_sto, 'O2_snow', nbp_glo, nsnow, nvm, itau_dep, &
        &               O2_snow, 'scatter', nbp_glo, indices_g)
  CALL restput_p (rest_id_sto, 'CH4_snow', nbp_glo, nsnow, nvm, itau_dep, &
        &               CH4_snow, 'scatter', nbp_glo, indices_g)
  CALL restput_p (rest_id_sto, 'altmax', nbp_glo, nvm, 1, itau_dep,     &
        &               altmax, 'scatter',  nbp_glo, indices_g)
  !-
  IF (is_root_prc) THEN
        !- Close restart files
        CALL getin_dump
        CALL restclo
  ENDIF
  !-
#ifdef CPP_PARA
  CALL MPI_FINALIZE(ier)
#endif
  WRITE(numout,*) "End of forcesoil."
  !--------------------
END PROGRAM forcesoil
