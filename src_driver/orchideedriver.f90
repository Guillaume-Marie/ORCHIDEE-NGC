! =================================================================================================================================
! PROGRAM       : orchideedriver
!
! CONTACT       : jan.polcher@lmd.jussieu.fr
!
! LICENCE      : IPSL (2016)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF      This is the main program for the new driver. This only organises the data and calls sechiba_main.
!!            The main work is done in glogrid.f90 and forcing_tools.f90.
!!
!!\n DESCRIPTION: Call the various modules to get the forcing data and provide it to SECHIBA. The only complexity
!!                is setting-up the domain decomposition and distributing the grid information.
!!                The code is parallel from tip to toe using the domain decomposition inherited from LMDZ.
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S)	:
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_driver/orchideedriver.f90 $
!! $Date: 2025-02-10 14:38:24 +0100 (lun. 10 févr. 2025) $
!! $Revision: 8824 $
!! \n
!_
!================================================================================================================================
! 
PROGRAM orchideedriver
  !---------------------------------------------------------------------
  !-
  !- 
  !---------------------------------------------------------------------
  USE defprec
  USE netcdf
  !
  !
  USE ioipsl_para
  USE mod_orchidee_para
  !
  USE grid
  USE time
  USE timer
  USE constantes
  USE constantes_soil
  USE forcing_tools
  USE globgrd
  !
  USE sechiba
  USE control
  USE ioipslctrl
  USE xios_orchidee
  USE function_library, ONLY: get_printlev

  !added in 202409
  USE grid_var, ONLY : grid_type, unstructured, regular_lonlat, regular_xy
  !
  !-
  IMPLICIT NONE
  !-
  CHARACTER(LEN=80) :: gridfilename
  CHARACTER(LEN=80), DIMENSION(100) :: forfilename
  INTEGER(i_std) :: nb_forcefile
  CHARACTER(LEN=8)  :: model_guess
  INTEGER(i_std)    ::  file_id
  INTEGER(i_std)    :: iim_glo, jjm_glo
  INTEGER(i_std)    :: nbseg
  INTEGER(i_std) :: nbindex_g
  REAL(r_std), DIMENSION(2) :: zoom_lon, zoom_lat
  LOGICAL           :: is_spinupacc
  !
  ! Variables for the global grid available on all procs and used
  ! to fill the ORCHIDEE variable on the root_proc
  !
  REAL(r_std), ALLOCATABLE, DIMENSION(:,:)    :: lalo_glo
  REAL(r_std), ALLOCATABLE, DIMENSION(:)      :: contfrac_glo
  REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: lon_glo, lat_glo, area_glo
  REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: mask_glo
  REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:,:) :: corners_glo
  INTEGER(i_std), ALLOCATABLE, DIMENSION(:) :: kindex_g
  CHARACTER(LEN=20)                           :: calendar
  !-
  !- Variables local to each processors.
  !-
  INTEGER(i_std) :: i, j, ik, jk, in, nbdt, first_point
  INTEGER(i_std) :: itau, itau_offset, itau_sechiba
  REAL(r_std)    :: date0, date0_shifted, dt, julian
  REAL(r_std)    :: date0_tmp, dt_tmp
  INTEGER(i_std) :: nbdt_tmp
  REAL(r_std)    :: timestep_interval(2), timestep_int_next(2)
  !The index of smallest lat lon of the zoomed region, in referene to the full region, for wrf safran grids
  ! safran latitudes are organized from south to north, longitudes from west to east
  ! the smallest lat and lon corresponds to the southern west corner of the region
  INTEGER, DIMENSION(2)          :: zoom_index_init  ! array organized in the order of (lon, lat) 
  !
  INTEGER(i_std) :: rest_id, rest_id_stom
  INTEGER(i_std) ::  hist_id, hist2_id, hist_id_stom, hist_id_stom_IPCC
  REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: lalo_loc
  INTEGER(i_std) :: iim, jjm, ier
  REAL(r_std), ALLOCATABLE, DIMENSION(:,:)   :: lon, lat
  REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:) :: corners_lon, corners_lat
  REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: bounds_lon_loc, bounds_lat_loc
  !REAL(r_std), ALLOCATABLE, DIMENSION(:) :: areas
  INTEGER(i_std) :: kjpindex
  INTEGER(i_std), ALLOCATABLE, DIMENSION(:) :: kindex
  !-
  !- input fields
  !-
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: u             !! Lowest level wind speed
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: v             !! Lowest level wind speed 
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: zlev_uv       !! Height of first layer
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: zlev_tq       !! Height of first layer
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: qair          !! Lowest level specific humidity
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: precip_rain   !! Rain precipitation
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: precip_snow   !! Snow precipitation
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: lwdown        !! Down-welling long-wave flux 
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: swdown        !! Downwelling surface short-wave flux
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: sinang        !! cosine of solar zenith angle
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: temp_air      !! Air temperature in Kelvin
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: epot_air      !! Air potential energy
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: ccanopy       !! CO2 concentration in the canopy
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: petAcoef      !! Coeficients A from the PBL resolution
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: peqAcoef      !! One for T and another for q
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: petBcoef      !! Coeficients B from the PBL resolution
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: peqBcoef      !! One for T and another for q
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: cdrag         !! Cdrag
  REAL(r_std), ALLOCATABLE, DIMENSION (:)             :: pb            !! Lowest level pressure
  !-
  !- output fields
  !-
  REAL(r_std), ALLOCATABLE, DIMENSION (:)            :: z0m           !! Surface roughness for momentum (m)
  REAL(r_std), ALLOCATABLE, DIMENSION (:)            :: z0h           !! Surface roughness for heat (m)
  REAL(r_std), ALLOCATABLE, DIMENSION (:)            :: coastalflow   !! Diffuse flow of water into the ocean (m^3/dt)
  REAL(r_std), ALLOCATABLE, DIMENSION (:)            :: riverflow     !! Largest rivers flowing into the ocean (m^3/dt)
  REAL(r_std), ALLOCATABLE, DIMENSION (:)            :: tsol_rad      !! Radiative surface temperature
  REAL(r_std), ALLOCATABLE, DIMENSION (:)            :: vevapp        !! Total of evaporation
  REAL(r_std), ALLOCATABLE, DIMENSION (:)            :: temp_sol_new  !! New soil temperature
  REAL(r_std), ALLOCATABLE, DIMENSION (:)            :: qsurf         !! Surface specific humidity
  REAL(r_std), ALLOCATABLE, DIMENSION (:,:)          :: albedo        !! Albedo
  REAL(r_std), ALLOCATABLE, DIMENSION (:)            :: fluxsens      !! Sensible chaleur flux
  REAL(r_std), ALLOCATABLE, DIMENSION (:)            :: fluxlat       !! Latent chaleur flux
  REAL(r_std), ALLOCATABLE, DIMENSION (:)            :: emis          !! Emissivity
  REAL(r_std), ALLOCATABLE, DIMENSION (:)            :: netco2        !! netco2flux: Sum CO2 flux over PFTs
  REAL(r_std), ALLOCATABLE, DIMENSION (:)            :: carblu        !! fco2_lu: Land Cover Change CO2 flux
  REAL(r_std), ALLOCATABLE, DIMENSION (:)            :: carbwh        !! fco2_wh: Wood harvest CO2 flux
  REAL(r_std), ALLOCATABLE, DIMENSION (:)            :: carbha        !! fco2_ha: Crop harvest CO2 flux
  REAL(r_std), ALLOCATABLE, DIMENSION (:,:)          :: veget_diag    !! Fraction of vegetation type (unitless, 0-1)
  REAL(r_std), ALLOCATABLE, DIMENSION (:,:)          :: lai_diag      !! Leaf area index (m^2 m^{-2}
  REAL(r_std), ALLOCATABLE, DIMENSION (:,:)          :: height_diag   !! Vegetation Height (m)
  REAL(r_std), ALLOCATABLE, DIMENSION (:)            :: run_off_lic      !! Calving, melting and liquid precipitation, needed when coupled to atmosphere
  REAL(r_std), ALLOCATABLE, DIMENSION (:)            :: run_off_lic_frac !! Cell fraction corresponding to run_off_lic

  !-
  !-
  !-
  REAL(r_std) :: atmco2
  REAL(r_std), ALLOCATABLE, DIMENSION (:)  :: u_tq, v_tq, swnet
  LOGICAL :: lrestart_read = .TRUE. !! Logical for _restart_ file to read
  LOGICAL :: lrestart_write = .FALSE. !! Logical for _restart_ file to write'
  !
  ! Timer variables
  !
  LOGICAL, PARAMETER :: timemeasure=.TRUE.
  REAL(r_std) :: waitput_cputime=0.0, waitget_cputime=0.0, orchidee_cputime=0.0
  REAL(r_std) :: waitput_walltime=0.0, waitget_walltime=0.0, orchidee_walltime=0.0
  !
  !
  ! Print point
  !
!!  REAL(r_std), DIMENSION(2) :: testpt=(/44.8,-25.3/)
!!  REAL(r_std), DIMENSION(2) :: testpt=(/44.8,-18.3/)
!!  REAL(r_std), DIMENSION(2) :: testpt=(/-60.25,-5.25/)
!!  REAL(r_std), DIMENSION(2) :: testpt=(/46.7,10.3/)
!!  REAL(r_std), DIMENSION(2) :: testpt=(/0.25,49.25/)
  ! Case when no ouput is desired.
  REAL(r_std), DIMENSION(2) :: testpt=(/9999.99,9999.99/)
  INTEGER(i_std) :: ktest, alloc_stat
  INTEGER :: printlev_loc  !! local write level
  INTEGER(i_std) :: ierr

  OFF_LINE_MODE = .TRUE. 


  !-
  !---------------------------------------------------------------------------------------
  !- 
  !- Define MPI communicator
  !- 
  !---------------------------------------------------------------------------------------
  !-
  !
  ! Set parallel processing in ORCHIDEE
  !
  CALL Init_orchidee_para()
  !
  !====================================================================================
  !
  ! Start timer now that the paralelisation is initialized.
  !
  IF ( timemeasure ) THEN
     CALL init_timer
     CALL start_timer(timer_global)
     CALL start_timer(timer_mpi)
  ENDIF
  !
  !
  !---------------------------------------------------------------------------------------
  !-
  !- Start the getconf processes 
  !-
  !---------------------------------------------------------------------------------------
  !-

  ! Set specific write level to orchideedriver using PRINTLEV_orchideedriver=[0-4] 
  ! in run.def. The global printlev is used as default value. 
  printlev_loc=get_printlev('orchideedriver')

  !-
  !Config Key   = GRID_FILE
  !Config Desc  = Name of file containing the forcing data
  !Config If    = [-]
  !Config Def   = grid_file.nc
  !Config Help  = This is the name of the file from which we will read
  !Config         or write into it the description of the grid from 
  !Config         the forcing file.
  !Config         compliant.
  !Config Units = [FILE] 
  !- 
  gridfilename='NONE'
  CALL getin_p('GRID_FILE', gridfilename)
  !-
  forfilename(:)=" "
  forfilename(1)='forcing_file.nc'
  CALL getin_p('FORCING_FILE', forfilename)


  !Config Key   = OK_SORT_LON_LAT
  !Config Desc  = Flag to order longitudes and latitudes in increasing order
  !Config If    = [orchideedriver]
  !Config Def   = n
  !Config Help  = If this flag is true, the latitudes and longitudes in the forcing file wil be
  !Config         changed to always be in increasing order. 
  !Config Units = [y/n]
  ok_sort_lon_lat=.FALSE.
  CALL getin_p('OK_SORT_LON_LAT', ok_sort_lon_lat)


  !-
  !- Define the zoom
  !-
  zoom_lon=(/-180,180/)
  zoom_lat=(/-90,90/)
  !
  !Config Key   = LIMIT_WEST
  !Config Desc  = Western limit of region
  !Config If    = [-]
  !Config Def   = -180.
  !Config Help  = Western limit of the region we are 
  !Config         interested in. Between -180 and +180 degrees
  !Config         The model will use the smalest regions from
  !Config         region specified here and the one of the forcing file.
  !Config Units = [Degrees] 
  !- 
  CALL getin_p('LIMIT_WEST',zoom_lon(1))
  !-
  !Config Key   = LIMIT_EAST
  !Config Desc  = Eastern limit of region
  !Config If    = [-]
  !Config Def   = 180.
  !Config Help  = Eastern limit of the region we are
  !Config         interested in. Between -180 and +180 degrees
  !Config         The model will use the smalest regions from
  !Config         region specified here and the one of the forcing file.
  !Config Units = [Degrees] 
  !-
  CALL getin_p('LIMIT_EAST',zoom_lon(2))
  !-
  !Config Key   = LIMIT_NORTH
  !Config Desc  = Northern limit of region
  !Config If    = [-]
  !Config Def   = 90.
  !Config Help  = Northern limit of the region we are
  !Config         interested in. Between +90 and -90 degrees
  !Config         The model will use the smalest regions from
  !Config         region specified here and the one of the forcing file.
  !Config Units = [Degrees]
  !-
  CALL getin_p('LIMIT_NORTH',zoom_lat(2))
  !-
  !Config Key   = LIMIT_SOUTH
  !Config Desc  = Southern limit of region
  !Config If    = [-]
  !Config Def   = -90.
  !Config Help  = Southern limit of the region we are
  !Config         interested in. Between 90 and -90 degrees
  !Config         The model will use the smalest regions from
  !Config         region specified here and the one of the forcing file.
  !Config Units = [Degrees]
  !-
  CALL getin_p('LIMIT_SOUTH',zoom_lat(1))
  IF ( (zoom_lon(1)+180 < EPSILON(zoom_lon(1))) .AND. (zoom_lon(2)-180 < EPSILON(zoom_lon(2))) .AND.&
       &(zoom_lat(1)+90 < EPSILON(zoom_lat(1))) .AND. (zoom_lat(2)-90 < EPSILON(zoom_lat(2))) ) THEN

     ! We are here only if zoom_lon and zoom_lat have there original values which
     ! means that they have not been modified by the getin LIMIT_ above. 
     ! Read WEST_EAST and SOUTH_NORTH from run.def.
     
     !Config Key   = WEST_EAST
     !Config Desc  = Longitude interval to use from the forcing data
     !Config If    = [-]
     !Config Def   = -180, 180
     !Config Help  = This function allows to zoom into the forcing data
     !Config Units = [degrees east] 
     !- 
     CALL getin_p('WEST_EAST', zoom_lon)
     !
     !Config Key   = SOUTH_NORTH
     !Config Desc  = Latitude interval to use from the forcing data
     !Config If    = [-]
     !Config Def   = -90, 90
     !Config Help  = This function allows to zoom into the forcing data
     !Config Units = [degrees north] 
     !- 
     CALL getin_p('SOUTH_NORTH', zoom_lat)
  ENDIF
  !-
  !-
  !- Get some basic variables from the run.def
  !-
  atmco2=350.
  CALL getin_p('ATM_CO2',atmco2)
  !
  !====================================================================================
  !-
  !-
  !- Get the grid on all processors.
  !-
  !---------------------------------------------------------------------------------------
  !-
  !- Read the grid, only on the root proc. from the forcing file using tools in the globgrd module. 
  !- The grid is then broadcast to all other broadcast.
  !
  nb_forcefile = 0
  DO ik=1,100
     IF ( INDEX(forfilename(ik), '.nc') > 0 ) nb_forcefile = nb_forcefile+1
  ENDDO
  !
  ! the output variables are defined only on root proc
  ! globgrd_getdomsz calls forcing_getglogrid in the case of regular grids only
  ! forcing_getglogrid defines _glo variables with SAVE attribute.
  ! in the case of unstructured/wrf grids, forcing_getglogrid is not called globgrd_getdomsz.
  IF ( is_root_prc) THEN
     CALL globgrd_getdomsz(gridfilename, iim_glo, jjm_glo, nbindex_g, nbseg, model_guess, zoom_index_init, &
          file_id, forfilename, zoom_lon, zoom_lat)
  ENDIF
  ! broadcast some to all processors: nbindex_g represents the number of valid land pixels
  CALL bcast(iim_glo)
  CALL bcast(jjm_glo)
  CALL bcast(nbindex_g)
  CALL bcast(nbseg)
  !-
  !- Allocation of memory
  !- variables over the entire grid (thus in x,y). Glo represents the processor 0, that to be broadcasted to all the processors
  ! here lon_glo etc have no SAVE attribute. They contains all pixels including oceans.
  ! It seems better to define such variables to _full_proc0, not _glo, to avoid confusion with global Saved variables on root by forcing_xxx
  ALLOCATE(lon_glo(iim_glo, jjm_glo),stat=alloc_stat)
  ALLOCATE(lat_glo(iim_glo, jjm_glo),stat=alloc_stat)
  ALLOCATE(mask_glo(iim_glo, jjm_glo),stat=alloc_stat)
  ALLOCATE(area_glo(iim_glo, jjm_glo),stat=alloc_stat)
  ALLOCATE(corners_glo(iim_glo, jjm_glo, nbseg, 2),stat=alloc_stat)
  !
  ! Gathered variables: nbindex_g should be only the land pixels size
  ! contfra_glo, kindex_g contains only land points
  ALLOCATE(kindex_g(nbindex_g),stat=alloc_stat)
  ALLOCATE(contfrac_glo(nbindex_g),stat=alloc_stat)
  !
  !-kindex_g: 1D index of land grids for zoomed region (if simulate zoomed region)
  ! or the land point index in the full map in case of full land domain (which is the case of unstructured grid).
  ! kindex_g: the values of indices are reference to the total pixels in the domain, including ocean
  ! The shape of output variables of globgrd_getgrid must be consistent to the allocation definition above
  ! otherwise the simulation will be freezed.
  ! call globgrd_getgrid: to get spatial domain variables. The variables _glo here are not global public variables used for the whole module, instead, local variables in orchideedriver.f90.
  ! contfrac_glo: contains only land pixels
  ! lon_glo and lat_glo: 2D array, in case of unstructured code, the first dimension(lon) is 1, because orchidee parallize in latitude. They contain the full pixels including oceans.
  ! corners_glo, mask_glo, area_glo also contains the ocean points and land points
  ! globgrd_getgrid gives corners_glo which is not a global variable with SAVE
  ! lat_glo and lon_glo contains missing values here
  ! corners_glo(:,:,:,1): longitude
  ! corners_glo(:,:,:,2): latitude
  ! calendar
  IF ( is_root_prc) THEN
     CALL globgrd_getgrid(file_id, iim_glo, jjm_glo, nbindex_g, nbseg, model_guess, zoom_index_init, &
          &               lon_glo, lat_glo, mask_glo, area_glo, corners_glo,&
          &               kindex_g, contfrac_glo, calendar, zoom_lon, zoom_lat)
  ENDIF
  !
  CALL bcast(lon_glo)
  CALL bcast(lat_glo)
  CALL bcast(mask_glo)
  CALL bcast(area_glo)
  CALL bcast(corners_glo)
  CALL bcast(kindex_g)
  CALL bcast(contfrac_glo)
  CALL bcast(calendar)
  CALL bcast(model_guess)  
  ! calendar will be broadcast to all processors. 
  !
  ALLOCATE(lalo_glo(nbindex_g,2), stat=alloc_stat)
  IF ( model_guess .EQ. "unstruct") THEN
     DO ik=1,nbindex_g
        i = iim_glo
        j = kindex_g(ik)
        lalo_glo(ik,1) = lat_glo(i,j)
        lalo_glo(ik,2) = lon_glo(i,j)
     ENDDO
  ELSE
     DO ik=1,nbindex_g
        ! j: lat, i: lon
        j = ((kindex_g(ik)-1)/iim_glo)+1
        i = (kindex_g(ik)-(j-1)*iim_glo)
        !
        IF ( i > iim_glo .OR. j > jjm_glo ) THEN
           WRITE(100+mpi_rank,*) "Error in the indexing (ik, kindex, i, j) : ", ik, kindex(ik), i, j
           STOP "ERROR in orchideedriver"
        ENDIF
        !
        lalo_glo(ik,1) = lat_glo(i,j)
        lalo_glo(ik,2) = lon_glo(i,j)
        !
     ENDDO
  ENDIF
  
  WRITE(*,*) "Rank", mpi_rank, " Before parallel region All land points : ",  nbindex_g
  WRITE(*,*) "Rank", mpi_rank, " from ", iim_glo, " point in Lon. and ", jjm_glo, "in Lat."
  !-
  !---------------------------------------------------------------------------------------
  !-
  !- Now that the grid is distributed on all procs we can start
  !- initialise the ORCHIDEE domain on each proc (longitude, latitude, indices)
  !-
  !---------------------------------------------------------------------------------------
  !-
  !- init_data_para also transfers kindex_g to index_g (the variable used in ORCHIDEE) 
  !- nbindex_g is used to define nbp_glo used by grid_set_glo
  !     In grid_allocate_glo:
  !     ALLOCATE(lalo_g(nbp_glo,2), contfrac_g(nbp_glo),index_g(nbp_glo))
  !     ALLOCATE(lon_g(iim_g, jjm_g), lat_g(iim_g, jjm_g))
  ! so index_g has the dimension of the number of land points
  CALL grid_set_glo(iim_glo, jjm_glo, nbindex_g)
  CALL grid_allocate_glo(nbseg)
  !
  ! Define is_spinupacc flag when the model_guess is unstruct
  ! ! If the number of land point is less than total points:
  ! then it is unstructured grids with ocean and land. is_spinupacc is false.
  ! If the number of land point is equal to total points:
  ! then it is the case of spinup accelation and optimization, with only land pixels
  is_spinupacc = .FALSE.
  IF ( model_guess == "unstruct" .AND. nbindex_g < jjm_glo ) THEN
     is_spinupacc = .FALSE.
  ELSE IF ( model_guess == "unstruct" .AND. nbindex_g .EQ. jjm_glo ) THEN
     is_spinupacc = .TRUE.
  ENDIF
  ! Copy the list of indexes of land points into index_g used by ORCHIDEE and then broacast to all
  ! processors
  ! index_g: defined in src_parallel/mod_orchidee_para_var.F90:  INTEGER(i_std),ALLOCATABLE,DIMENSION(:),SAVE   :: index_g
  ! allocated in grid.f90
  ! index_g: index of land points, referenced to the total pixels on the domain (ocean included)
  CALL bcast(nbindex_g)
  IF ( is_root_prc) index_g = kindex_g
  CALL bcast(index_g)
  
  ! nbindex_g: number of valid land grids in the domain
  WRITE(numout,*) "Rank", mpi_rank, "Into Init_orchidee_data_para_driver with ", nbindex_g
  ! to print out the index of the first and the last land point in the global scale
  WRITE(numout,*) "Rank", mpi_rank, "Into ", index_g(1), index_g(nbindex_g)
  CALL flush(numout)
  !
  ! the driver uses the land points for computation for the parallization
  CALL Init_orchidee_data_para_driver(nbindex_g,index_g)
  CALL init_ioipsl_para 
  WRITE(numout,*) "Rank", mpi_rank, "After init_data_para global size : ",  nbp_glo, SIZE(index_g), iim_g, iim_glo, jjm_g, jjm_glo
  IF ( model_guess .NE. "unstruct") THEN
     WRITE(numout,'("After init_data_para local : ij_nb, jj_nb",2I4)') iim_glo, jj_nb
  ELSE IF ( model_guess .EQ. "unstruct") THEN 
     WRITE(numout,'("After init_data_para local : ij_nb, jj_nb",2I4)') jj_nb, jj_nb
  ENDIF
  CALL flush(numout)
  !
  ! Allocate grid varibles on the local processor by using nbp_loc:
  !
  IF ( model_guess == "regular") THEN
     CALL grid_init (nbp_loc, nbseg, regular_lonlat, "ForcingGrid")
  ELSE IF ( model_guess == "WRF") THEN
     CALL grid_init (nbp_loc, nbseg, regular_xy, "WRFGrid")
  ELSE IF ( model_guess == "unstruct") THEN
     CALL grid_init (nbp_loc, nbseg, unstructured, "Unstructured")
  ELSE
     CALL ipslerr(3, "orchideedriver", "The grid found in the GRID_FILE is not supported by ORCHIDEE", "", "")
  ENDIF
  !
  !
  ! Transfer the global grid variables to the ORCHIDEE version on the root proc
  ! *_glo -> *_g
  ! Variables *_g were allocated with the CALL init_grid
  !
  !lalo_g: not really useful. can be removed 202311
  lalo_g(:,:) = lalo_glo(:,:)
  contfrac_g(:) = contfrac_glo(:)
  lon_g(:,:) = lon_glo(:,:)
  lat_g(:,:) = lat_glo(:,:)
  !
  !
  ! Set the local dimensions of the fields on each processor
  !
  iim = iim_glo
  jjm = jj_nb
  ! kjpindex = number of land pixels on each processor, used by the driver
  kjpindex = nbp_loc
  !
  WRITE(numout,*) mpi_rank, "DIMENSIONS of grid on processor : iim, jjm, kjpindex = ", iim, jjm, kjpindex
  CALL flush(numout)
  !
  ALLOCATE(lon(iim,jjm), lat(iim,jjm))
  ALLOCATE(corners_lon(nbseg,iim,jj_para_nb(mpi_rank)), corners_lat(nbseg,iim,jj_para_nb(mpi_rank)))
  ALLOCATE(kindex(kjpindex))
  IF ( .NOT. Allocated(contfrac)) ALLOCATE(contfrac(kjpindex))
  !
  ! lat lon used for xios_orchidee_init and grid_init_unstructuured
  lon=lon_glo(:,jj_para_begin(mpi_rank):jj_para_end(mpi_rank))
  lat=lat_glo(:,jj_para_begin(mpi_rank):jj_para_end(mpi_rank))
  ! corners_lon is used in xios_orchidee_init
  DO in=1,nbseg
     DO j=1,jj_para_nb(mpi_rank)
        jk=jj_para_begin(mpi_rank)+(j-1)
        DO i=1,iim
           ! Copy the corners/bounds of the meshes which are local to this processor.
           ! Serves for OASIS and XIOS.
           corners_lon(in,i,j)=corners_glo(i,jk,in,1)
           corners_lat(in,i,j)=corners_glo(i,jk,in,2)
        ENDDO
     ENDDO
  ENDDO
  !
  ! Redistribute the indeces on all procs (appel distribution of land points) 
  ! lon_g and lag_g: allocated in grid.f90 ALLOCATE(lon_g(iim_g, jjm_g), lat_g(iim_g, jjm_g))
  ! iim_g and jjm_g equal to the input variables in grid_set_glo, which is iim_glo, jjm_glo
  CALL bcast(lon_g)
  CALL bcast(lat_g)
  ! distribute index_g to multiple processors, with each processor having the kindex
  CALL scatter(index_g, kindex)
  !
  ! Apply the offset needed so that kindex refers to the index of the land point 
  ! on the current region, i.e. the local lon lat domain.
  !
  kindex(1:kjpindex)=kindex(1:kjpindex)-(jj_begin-1)*iim_glo
  !
  !
  ! This routine transforms the global grid into a series of polygons for all land
  ! points identified by index_g.
  !
  IF ( model_guess /= "unstruct") THEN
     ! it scatters useful grid informations to the processors, such as corner, nieghbour, area, segment etc
     CALL grid_stuff(nbindex_g, iim_g, jjm_g, lon_g, lat_g, index_g, contfrac_glo)
  ELSE
     CALL scatter(contfrac_glo, contfrac)
  ENDIF
  !
  ! Distribute the global lalo to the local processor level lalo
  ! lalo is same as lat and lon
  ALLOCATE(lalo_loc(kjpindex,2))
  CALL scatter(lalo_glo, lalo_loc)
  lalo(:,:) = lalo_loc(:,:)
  !
  !====================================================================================
  !-
  !- Prepare the time for the simulation
  !-
  !- Set the calendar and get some information
  !-
  CALL ioconf_calendar(calendar)
  CALL ioget_calendar(one_year, one_day)
  !-
  !- get the time period for the run
  !-
  CALL forcing_integration_time(date0, dt, nbdt)
  !
  !-
  !- Set the start date in IOIPSL for the calendar and initialize the module time
  !-
  CALL ioconf_startdate(date0)
  CALL time_initialize(0, date0, dt, "END")
  !
  !
  !====================================================================================
  !-
  !- Initialize the forcing files and prepare the time stepping through the data.
  !-
  !  This subroutine calls forcing_getglogrid and forcing_zoomgrid, forcing_vertical, forcing_integrationtime etc
  !  subroutine forcing_getglogrid defines several variables _glo, which have SAVE attribute, and can be shared within the whole module
  !   would these variables be in conflict with get_globgrid ? or getdomsz ?
  ! forcing_open call forcing_getglogrid and forcing_zoomgrid
  ! forcing_getglogrid will also compute nbland_glo, iim_glo, jjm_glo, lat_glo, lon_glo,
  ! contfrac_glo, lindex_glo for the regular grids.
  !
  ! forcing_zoomgrid will compute iim_loc, jjm_loc, nbpoint_loc for all types of grid
  ! It computes also corners_loc, areas_loc for wrf and regular grids based on other subroutines
  ! It finds computes corners_loc, areas_loc for for unstructured grids based on corners_glo, areas_glo variables
  !
  CALL forcing_open(forfilename, iim_glo,  jjm_glo, lon_glo, lat_glo, nbindex_g, zoom_lon, zoom_lat, &
       &            index_g, kjpindex, numout, model_guess, zoom_index_init)
  !
  ! Moyenne des precipitations du forcage, par maille. DOIT rester ici : entre forcing_open
  ! (qui ouvre les fichiers et fixe la grille) et sechiba_initialize (dont slowproc_init
  ! consomme le resultat au demarrage a froid pour initialiser precip_longterm). Deplacer
  ! cet appel apres sechiba_initialize le rendrait silencieusement inoperant.
  IF ( .NOT. ALLOCATED(precip_forcing_mean) ) ALLOCATE(precip_forcing_mean(kjpindex))
  CALL forcing_meanprecip(model_guess, precip_forcing_mean)
  !
  !-
  !====================================================================================
  !-
  !- Initialise the ORCHIDEE system in 4 steps :
  !- 1 The control flags,
  !- 2 Allocate memory (needs to be done after initializing the control flags because of nvm). 
  !- 2 the restart system of IOIPSL
  !- 3 The history mechanism
  !- 4 Finally the first call to SECHIBA will initialise all the internal variables
  !
  ! 1 Setting flags and some parameters
  ! 
  ! Fixed PFT dependent parameters (nvm)
  CALL control_initialize(kjpindex)
  !
  !
  ! 2 Allocation
  !
  ALLOCATE(zlev_tq(kjpindex), zlev_uv(kjpindex))
  ALLOCATE(u(kjpindex), v(kjpindex), pb(kjpindex))
  ALLOCATE(temp_air(kjpindex))
  ALLOCATE(qair(kjpindex))
  ALLOCATE(petAcoef(kjpindex), peqAcoef(kjpindex), petBcoef(kjpindex), peqBcoef(kjpindex))
  ALLOCATE(ccanopy(kjpindex))
  ALLOCATE(cdrag(kjpindex))
  ALLOCATE(precip_rain(kjpindex))
  ALLOCATE(precip_snow(kjpindex))
  ALLOCATE(swdown(kjpindex))
  ALLOCATE(swnet(kjpindex))
  ALLOCATE(lwdown(kjpindex))
  ALLOCATE(sinang(kjpindex))
  ALLOCATE(vevapp(kjpindex))
  ALLOCATE(fluxsens(kjpindex))
  ALLOCATE(fluxlat(kjpindex))
  ALLOCATE(coastalflow(kjpindex))
  ALLOCATE(riverflow(kjpindex))
  ALLOCATE(netco2(kjpindex))
  ALLOCATE(carblu(kjpindex))
  ALLOCATE(carbwh(kjpindex))
  ALLOCATE(carbha(kjpindex))
  ALLOCATE(tsol_rad(kjpindex))
  ALLOCATE(temp_sol_new(kjpindex))
  ALLOCATE(qsurf(kjpindex))
  ALLOCATE(albedo(kjpindex,2))
  ALLOCATE(emis(kjpindex))
  ALLOCATE(epot_air(kjpindex))
  ALLOCATE(u_tq(kjpindex), v_tq(kjpindex))
  ALLOCATE(z0m(kjpindex))
  ALLOCATE(z0h(kjpindex))
  ALLOCATE(veget_diag(kjpindex,nvm))
  ALLOCATE(lai_diag(kjpindex,nvm))
  ALLOCATE(height_diag(kjpindex,nvm))
  ALLOCATE(run_off_lic(kjpindex))
  ALLOCATE(run_off_lic_frac(kjpindex))
  !-
  !---------------------------------------------------------------------------------------
  !-
  !- Get a first set of forcing data
  !-
  !---------------------------------------------------------------------------------------
  !-
  !- Some default values so that the operations before the ORCHIDEE initialisation do not fail.
  !-
  z0m(:) = 0.1
  albedo(:,:) = 0.13
  !
  itau = 0
  !
  CALL ioipslctrl_restini(itau, date0, dt, rest_id, rest_id_stom, itau_offset, date0_shifted)
  !
  ! To ensure that itau starts with 0 at date0 for the restart, we have to set an off-set to achieve this. 
  ! itau_offset will get used to prduce itau_sechiba.
  !
  itau_offset=-itau_offset-1
  !
  ! Get the date of the first time step
  !
  WRITE(*,*) "itau_offset : date0 : ", year_end, month_end, day_end, sec_end
  !
  !!- Initialize module for output with XIOS
  CALL xios_orchidee_init( MPI_COMM_ORCH,                         &  
       date0,    year_end, month_end,   day_end,    julian_diff,  &
       lon,      lat,      corners_lon, corners_lat, is_spinupacc)
  CALL sechiba_xios_initialize
  CALL xios_orchidee_close_definition
  IF (printlev_loc >= 2) WRITE(numout,*) 'After xios_orchidee_close_definition'

  !- Initialize IOIPSL sechiba output files
  itau_sechiba = itau+itau_offset
  CALL ioipslctrl_history(iim, jjm, lon, lat,  kindex, kjpindex, itau_sechiba, &
       date0, dt, hist_id, hist2_id, hist_id_stom, hist_id_stom_IPCC)
  WRITE(*,*) "HISTORY : Defined for ", itau_sechiba, date0, dt
  !-
  !---------------------------------------------------------------------------------------
  !-
  !- Go into the time loop
  !-
  !---------------------------------------------------------------------------------------
  !-
  DO itau = 1,nbdt
     !
     CALL time_nextstep( itau )
     !     
     timestep_interval(1) = julian_start
     timestep_interval(2) = julian_end
     julian = (julian_start + julian_end) /2.0  !julian_end
     !
     ! Get the forcing data
     !
     CALL forcing_getvalues(timestep_interval, dt, zlev_tq, zlev_uv, temp_air, qair, &
          &                 precip_rain, precip_snow, swdown, lwdown, sinang, u, v, pb, model_guess, is_spinupacc)
     !-
     IF ( itau == nbdt ) lrestart_write = .TRUE.
     !
     ! Adaptation of the forcing data to SECHIBA's needs
     !
     ! Contrary to what the documentation says, ORCHIDEE expects surface pressure in hPa.
     pb(:) = pb(:)/100.
     epot_air(:) = cp_air*temp_air(:)+cte_grav*zlev_tq(:)
     ccanopy(:) = atmco2
     cdrag(:) = 0.0
     !
     petBcoef(:) = epot_air(:)
     peqBcoef(:) = qair(:)
     petAcoef(:) = zero
     peqAcoef(:) = zero
     !
     ! Interpolate the wind (which is at hight zlev_uv) to the same height 
     ! as the temperature and humidity (at zlev_tq).
     !
     u_tq(:) = u(:)*LOG(zlev_tq(:)/z0m(:))/LOG(zlev_uv(:)/z0m(:))
     v_tq(:) = v(:)*LOG(zlev_tq(:)/z0m(:))/LOG(zlev_uv(:)/z0m(:))
     !
     swnet(:) =(1.-(albedo(:,1)+albedo(:,2))/2.)*swdown(:)
     !
     CALL forcing_printpoint(julian, testpt(1), testpt(2), kjpindex, lalo_loc, temp_air, "RECEIVED Air temperature")
     CALL forcing_printpoint(julian, testpt(1), testpt(2), kjpindex, lalo_loc, qair, "RECEIVED Air humidity")
     CALL forcing_printpoint(julian, testpt(1), testpt(2), kjpindex, lalo_loc, precip_rain*one_day, "RECEIVED Rainfall")
     CALL forcing_printpoint(julian, testpt(1), testpt(2), kjpindex, lalo_loc, precip_snow*one_day, "RECEIVED Snowfall")
     CALL forcing_printpoint(julian, testpt(1), testpt(2), kjpindex, lalo_loc, swnet, "RECEIVED net solar")
     CALL forcing_printpoint(julian, testpt(1), testpt(2), kjpindex, lalo_loc, lwdown, "RECEIVED lwdown")
     CALL forcing_printpoint(julian, testpt(1), testpt(2), kjpindex, lalo_loc, u, "RECEIVED East-ward wind")
     CALL forcing_printpoint(julian, testpt(1), testpt(2), kjpindex, lalo_loc, v, "RECEIVED North-ward wind")
     CALL forcing_printpoint(julian, testpt(1), testpt(2), kjpindex, lalo_loc, pb*100, "RECEIVED surface pressure")
     CALL forcing_printpoint(julian, testpt(1), testpt(2), kjpindex, lalo_loc, zlev_uv, "RECEIVED UV height")
     CALL forcing_printpoint(julian, testpt(1), testpt(2), kjpindex, lalo_loc, zlev_tq, "RECEIVED TQ height")
     CALL forcing_printpoint(julian, testpt(1), testpt(2), kjpindex, lalo_loc, sinang, "RECEIVED sinang")
     !
     IF ( itau .NE. 1 ) THEN
        IF ( timemeasure ) THEN
           waitget_cputime = waitget_cputime + Get_cpu_Time(timer_global)
           waitget_walltime = waitget_walltime + Get_real_Time(timer_global)
           CALL stop_timer(timer_global)
           CALL start_timer(timer_global)
        ENDIF
     ENDIF
     !
     !---------------------------------------------------------------------------------------
     !-
     !- IF first time step : Call to SECHIBA_initialize to set-up ORCHIDEE before doing an actual call
     !- which will provide the first fluxes.
     !-
     !---------------------------------------------------------------------------------------
     !
     itau_sechiba = itau+itau_offset
     !
     ! Update the calendar in xios by sending the new time step
     !CALL xios_orchidee_update_calendar(itau_sechiba)
     CALL xios_orchidee_update_calendar(itau_sechiba)
     !
     IF ( itau == 1 ) THEN
        !
        IF ( timemeasure ) THEN
           WRITE(numout,*) '------> CPU Time for start-up of main : ',Get_cpu_Time(timer_global)
           WRITE(numout,*) '------> Real Time for start-up of main : ',Get_real_Time(timer_global)
           CALL stop_timer(timer_global)
           CALL start_timer(timer_global)
        ENDIF
        ! contfrac: contfrac values distributed on one procesor.
        ! defined by grid_scatter
        ! lalo_loc: only land points
        CALL sechiba_initialize( &
             itau_sechiba,  iim*jjm,      kjpindex,      kindex,                   &
             lalo_loc,     contfrac,     neighbours,    resolution,  zlev_tq,      &
             u_tq,         v_tq,         qair,          temp_air,                  &
             petAcoef,     peqAcoef,     petBcoef,      peqBcoef,                  &
             precip_rain,  precip_snow,  lwdown,        swnet,       swdown,       &
             pb,           rest_id,      hist_id,       hist2_id,                  &
             rest_id_stom, hist_id_stom, hist_id_stom_IPCC,                        &
             coastalflow,  riverflow,    tsol_rad,      vevapp,      qsurf,        &
             z0m,          z0h,          albedo,        fluxsens,    fluxlat,  emis, &
             temp_sol_new, cdrag,       sinang)
        !
        CALL forcing_printpoint(julian, testpt(1), testpt(2), kjpindex, lalo_loc, temp_sol_new, "Init temp_sol_new")
        !
        ! Net solar and the wind at the right hight are recomputed with the correct values.
        !
        swnet(:) =(1.-(albedo(:,1)+albedo(:,2))/2.)*swdown(:)
        u_tq(:) = u(:)*LOG(zlev_tq(:)/z0m(:))/LOG(zlev_uv(:)/z0m(:))
        v_tq(:) = v(:)*LOG(zlev_tq(:)/z0m(:))/LOG(zlev_uv(:)/z0m(:))
        !
        lrestart_read = .FALSE.
        !
        CALL histwrite_p(hist_id, 'LandPoints',  itau+1, (/ REAL(kindex) /), kjpindex, kindex)
        !CALL histwrite_p(hist_id, 'Areas',  itau+1, area, kjpindex, kindex)
        CALL histwrite_p(hist_id, 'Contfrac',  itau+1, contfrac, kjpindex, kindex)
        !
        IF ( timemeasure ) THEN
           WRITE(numout,*) '------> CPU Time for set-up of ORCHIDEE : ',Get_cpu_Time(timer_global)
           WRITE(numout,*) '------> Real Time for set-up of ORCHIDEE : ',Get_real_Time(timer_global)
           CALL stop_timer(timer_global)
           CALL start_timer(timer_global)
        ENDIF
        !
     ENDIF
     !
     !---------------------------------------------------------------------------------------
     !-
     !- Main call to SECHIBA  
     !-
     !---------------------------------------------------------------------------------------
     !
     !
     !
     CALL sechiba_main (itau_sechiba, iim*jjm, kjpindex, kindex, &
          & lrestart_read, lrestart_write, &
          & lalo_loc, contfrac, neighbours, resolution, &
          ! First level conditions
          & zlev_tq, u_tq, v_tq, qair, temp_air, epot_air, ccanopy, &
          ! Variables for the implicit coupling
          & cdrag, petAcoef, peqAcoef, petBcoef, peqBcoef, &
          ! Rain, snow, radiation and surface pressure
          & precip_rain ,precip_snow, lwdown, swnet, swdown, sinang, pb, &
          ! Output : Fluxes
          & vevapp, fluxsens, fluxlat, coastalflow, riverflow, &
          ! CO2 fluxes
          & netco2, carblu, carbwh, carbha, &
          ! Surface temperatures and surface properties
          & tsol_rad, temp_sol_new, qsurf, albedo, emis, z0m, z0h, &
          ! Vegetation, lai and vegetation height
          & veget_diag, lai_diag, height_diag, &
          ! File ids
          & rest_id, hist_id, hist2_id, rest_id_stom, hist_id_stom, hist_id_stom_IPCC, &
          ! Variables related to continental ice
          & run_off_lic, run_off_lic_frac)
     !
     CALL forcing_printpoint(julian, testpt(1), testpt(2), kjpindex, lalo_loc, temp_sol_new, "Produced temp_sol_new")
     CALL forcing_printpoint(julian, testpt(1), testpt(2), kjpindex, lalo_loc, fluxsens, "Produced fluxsens")
     CALL forcing_printpoint(julian, testpt(1), testpt(2), kjpindex, lalo_loc, fluxlat, "Produced fluxlat")
     !
     IF ( timemeasure ) THEN
        orchidee_cputime = orchidee_cputime + Get_cpu_Time(timer_global)
        orchidee_walltime = orchidee_walltime + Get_real_Time(timer_global)
        CALL stop_timer(timer_global)
        CALL start_timer(timer_global)
     ENDIF
     !
     !---------------------------------------------------------------------------------------
     !-
     !- Write diagnostics
     !-
     !---------------------------------------------------------------------------------------
     !
     CALL xios_orchidee_send_field("LandPoints" ,(/ ( REAL(ik), ik=1,kjpindex ) /))
     CALL xios_orchidee_send_field("areas", area)
     CALL xios_orchidee_send_field("contfrac",contfrac)
     CALL xios_orchidee_send_field("temp_air",temp_air)
     CALL xios_orchidee_send_field("qair",qair)
     CALL xios_orchidee_send_field("swnet",swnet)
     CALL xios_orchidee_send_field("swdown",swdown)
     ! zpb in hPa, output in Pa
     CALL xios_orchidee_send_field("pb",pb)
     !
     IF ( .NOT. almaoutput ) THEN
        ! 
        !  ORCHIDEE INPUT variables
        !
        CALL histwrite_p (hist_id, 'swdown',   itau_sechiba, swdown,   kjpindex, kindex)
        CALL histwrite_p (hist_id, 'tair',     itau_sechiba, temp_air, kjpindex, kindex)
        CALL histwrite_p (hist_id, 'qair',     itau_sechiba, qair, kjpindex, kindex)
        CALL histwrite_p (hist_id, 'evap',     itau_sechiba, vevapp, kjpindex, kindex)
        CALL histwrite_p (hist_id, 'coastalflow',itau_sechiba, coastalflow, kjpindex, kindex)
        CALL histwrite_p (hist_id, 'riverflow',itau_sechiba, riverflow, kjpindex, kindex)
        ! 
        CALL histwrite_p (hist_id, 'temp_sol', itau_sechiba, temp_sol_new, kjpindex, kindex)
        CALL histwrite_p (hist_id, 'tsol_max', itau_sechiba, temp_sol_new, kjpindex, kindex)
        CALL histwrite_p (hist_id, 'tsol_min', itau_sechiba, temp_sol_new, kjpindex, kindex)
        CALL histwrite_p (hist_id, 'fluxsens', itau_sechiba, fluxsens, kjpindex, kindex)
        CALL histwrite_p (hist_id, 'fluxlat',  itau_sechiba, fluxlat,  kjpindex, kindex)
        CALL histwrite_p (hist_id, 'swnet',    itau_sechiba, swnet,    kjpindex, kindex)
        CALL histwrite_p (hist_id, 'alb_vis',  itau_sechiba, albedo(:,1), kjpindex, kindex)
        CALL histwrite_p (hist_id, 'alb_nir',  itau_sechiba, albedo(:,2), kjpindex, kindex)
        !
        IF ( hist2_id > 0 ) THEN
           CALL histwrite_p (hist2_id, 'swdown',   itau_sechiba, swdown, kjpindex, kindex)
           CALL histwrite_p (hist2_id, 'tair',     itau_sechiba, temp_air, kjpindex, kindex)
           CALL histwrite_p (hist2_id, 'qair',     itau_sechiba, qair, kjpindex, kindex)
           !
           CALL histwrite_p (hist2_id, 'evap',     itau_sechiba, vevapp, kjpindex, kindex)
           CALL histwrite_p (hist2_id, 'coastalflow',itau_sechiba, coastalflow, kjpindex, kindex)
           CALL histwrite_p (hist2_id, 'riverflow',itau_sechiba, riverflow, kjpindex, kindex)
           ! 
           CALL histwrite_p (hist2_id, 'temp_sol', itau_sechiba, temp_sol_new, kjpindex, kindex)
           CALL histwrite_p (hist2_id, 'tsol_max', itau_sechiba, temp_sol_new, kjpindex, kindex)
           CALL histwrite_p (hist2_id, 'tsol_min', itau_sechiba, temp_sol_new, kjpindex, kindex)
           CALL histwrite_p (hist2_id, 'fluxsens', itau_sechiba, fluxsens, kjpindex, kindex)
           CALL histwrite_p (hist2_id, 'fluxlat',  itau_sechiba, fluxlat,  kjpindex, kindex)
           CALL histwrite_p (hist2_id, 'swnet',    itau_sechiba, swnet,    kjpindex, kindex)
           !
           CALL histwrite_p (hist2_id, 'alb_vis',  itau_sechiba, albedo(:,1), kjpindex, kindex)
           CALL histwrite_p (hist2_id, 'alb_nir',  itau_sechiba, albedo(:,2), kjpindex, kindex)
        ENDIF
     ELSE
        !
        ! Input variables
        !
        CALL histwrite_p (hist_id, 'SinAng', itau_sechiba, sinang, kjpindex, kindex)
        CALL histwrite_p (hist_id, 'LWdown', itau_sechiba, lwdown, kjpindex, kindex)
        CALL histwrite_p (hist_id, 'SWdown', itau_sechiba, swdown, kjpindex, kindex)
        CALL histwrite_p (hist_id, 'Tair', itau_sechiba, temp_air, kjpindex, kindex)
        CALL histwrite_p (hist_id, 'Qair', itau_sechiba, qair, kjpindex, kindex)
        CALL histwrite_p (hist_id, 'SurfP', itau_sechiba, pb, kjpindex, kindex)
        CALL histwrite_p (hist_id, 'Windu', itau_sechiba, u_tq, kjpindex, kindex)
        CALL histwrite_p (hist_id, 'Windv', itau_sechiba, v_tq, kjpindex, kindex)
        !
        CALL histwrite_p (hist_id, 'Evap', itau_sechiba, vevapp, kjpindex, kindex)
        CALL histwrite_p (hist_id, 'SWnet',    itau_sechiba, swnet, kjpindex, kindex)
        CALL histwrite_p (hist_id, 'Qh', itau_sechiba, fluxsens, kjpindex, kindex)
        CALL histwrite_p (hist_id, 'Qle',  itau_sechiba, fluxlat, kjpindex, kindex)
        CALL histwrite_p (hist_id, 'AvgSurfT', itau_sechiba, temp_sol_new, kjpindex, kindex)
        CALL histwrite_p (hist_id, 'RadT', itau_sechiba, temp_sol_new, kjpindex, kindex)
        !
        ! There is a mess with the units passed to the coupler. To be checked with Marc
        !
        IF ( river_routing ) THEN
           CALL histwrite_p (hist_id, 'CoastalFlow',itau_sechiba, coastalflow, kjpindex, kindex)
           CALL histwrite_p (hist_id, 'RiverFlow',itau_sechiba, riverflow, kjpindex, kindex)
        ENDIF
        !
        IF ( hist2_id > 0 ) THEN
           CALL histwrite_p (hist2_id, 'Evap', itau_sechiba, vevapp, kjpindex, kindex)
           CALL histwrite_p (hist2_id, 'SWnet',    itau_sechiba, swnet, kjpindex, kindex)
           CALL histwrite_p (hist2_id, 'Qh', itau_sechiba, fluxsens, kjpindex, kindex)
           CALL histwrite_p (hist2_id, 'Qle',  itau_sechiba, fluxlat, kjpindex, kindex)
           CALL histwrite_p (hist2_id, 'AvgSurfT', itau_sechiba, temp_sol_new, kjpindex, kindex)
           CALL histwrite_p (hist2_id, 'RadT', itau_sechiba, temp_sol_new, kjpindex, kindex)
        ENDIF
     ENDIF
     !
     !
  ENDDO
  !-
  !-
  !---------------------------------------------------------------------------------------
  !-
  !- Close eveything
  !-
  !--
  ! Close IOIPSL history files
  CALL histclo
  IF(is_root_prc) THEN
     ! Close restart files
     CALL restclo
     CALL getin_dump
  ENDIF
  !-
  !- Deallocate all variables and reset initialization variables
  !-
  CALL orchideedriver_clear()
  ! close forcing files which was open by forcing_open
  IF ( itau == nbdt ) THEN
     DO ik=1,nb_forcefile
       ierr = NF90_CLOSE(force_id(ik))
    ENDDO
  ENDIF
  !
  WRITE(numout,*) "End at proc ", mpi_rank
  !
  !
  !---------------------------------------------------------------------------------------
  !-
  !- Get time and close IOIPSL, OASIS and MPI
  !-
  !---------------------------------------------------------------------------------------
  !-
  IF ( timemeasure ) THEN
     WRITE(numout,*) '------> Total CPU Time waiting to get forcing : ',waitget_cputime
     WRITE(numout,*) '------> Total Real Time waiting to get forcing : ',waitget_walltime
     WRITE(numout,*) '------> Total CPU Time for ORCHIDEE : ', orchidee_cputime
     WRITE(numout,*) '------> Total Real Time for ORCHIDEE : ', orchidee_walltime
     WRITE(numout,*) '------> Total CPU Time waiting to put fluxes : ',waitput_cputime
     WRITE(numout,*) '------> Total Real Time waiting to put fluxes : ',waitput_walltime
     WRITE(numout,*) '------> Total CPU Time for closing : ',  Get_cpu_Time(timer_global)
     WRITE(numout,*) '------> Total Real Time for closing : ', Get_real_Time(timer_global)
     WRITE(numout,*) '------> Total without MPI : CPU Time  : ', Get_cpu_Time(timer_mpi)
     WRITE(numout,*) '------> Total without MPI : Real Time : ', Get_real_Time(timer_mpi)
     CALL stop_timer(timer_global)
     CALL stop_timer(timer_mpi)
  ENDIF
  !
  ! Finalize MPI and XIOS
  CALL Finalize_mpi
 
  IF (printlev>=1) WRITE(numout,*) 'END of orchideedriver'
 
CONTAINS
!
!! ================================================================================================================================
!! SUBROUTINE   : orchideedriver_clear
!!
!>\BRIEF         Clear orchideedriver
!!
!! DESCRIPTION  :  Deallocate memory and reset initialization variables to there original values
!!                 This subroutine calls forcing_tools_clear and sechiba_clear.
!!
!_ ================================================================================================================================
!
  SUBROUTINE orchideedriver_clear
    !- Deallocate all variables existing on all procs
    !-
    !- Deallocate all variables existing on all procs (list still incomplete)
    !-
    IF ( ALLOCATED(lon_glo) ) DEALLOCATE(lon_glo)
    IF ( ALLOCATED(lat_glo) ) DEALLOCATE(lat_glo)
    IF ( ALLOCATED(mask_glo) ) DEALLOCATE(mask_glo)
    IF ( ALLOCATED(area_glo) ) DEALLOCATE(area_glo)
    IF ( ALLOCATED(corners_glo) ) DEALLOCATE(corners_glo)
    IF ( ALLOCATED(corners_lon) ) DEALLOCATE(corners_lon)
    IF ( ALLOCATED(corners_lat) ) DEALLOCATE(corners_lat)
    IF ( ALLOCATED(kindex_g) ) DEALLOCATE(kindex_g)
    IF ( ALLOCATED(contfrac_glo) ) DEALLOCATE(contfrac_glo)
    IF ( ALLOCATED(lalo_glo) ) DEALLOCATE(lalo_glo)
    IF ( ALLOCATED(lon) ) DEALLOCATE(lon)
    IF ( ALLOCATED(lat) ) DEALLOCATE(lat)
    IF ( ALLOCATED(kindex) ) DEALLOCATE(kindex)
    IF ( ALLOCATED(lalo_loc) ) DEALLOCATE(lalo_loc)
    IF ( ALLOCATED(zlev_tq) ) DEALLOCATE(zlev_tq)
    IF ( ALLOCATED(zlev_uv) ) DEALLOCATE(zlev_uv)
    IF ( ALLOCATED(u) ) DEALLOCATE(u)
    IF ( ALLOCATED(v) ) DEALLOCATE(v)
    IF ( ALLOCATED(pb) ) DEALLOCATE(pb)
    IF ( ALLOCATED(temp_air) ) DEALLOCATE(temp_air)
    IF ( ALLOCATED(qair) ) DEALLOCATE(qair)
    IF ( ALLOCATED(precip_rain) ) DEALLOCATE(precip_rain)
    IF ( ALLOCATED(precip_snow) ) DEALLOCATE(precip_snow)
    IF ( ALLOCATED(swdown) ) DEALLOCATE(swdown)
    IF ( ALLOCATED(swnet) ) DEALLOCATE(swnet)
    IF ( ALLOCATED(lwdown) ) DEALLOCATE(lwdown)
    IF ( ALLOCATED(sinang) ) DEALLOCATE(sinang)
    IF ( ALLOCATED(epot_air) ) DEALLOCATE(epot_air)
    IF ( ALLOCATED(ccanopy) ) DEALLOCATE(ccanopy)
    IF ( ALLOCATED(cdrag) ) DEALLOCATE(cdrag)
    IF ( ALLOCATED(swnet) ) DEALLOCATE(swnet)
    IF ( ALLOCATED(petAcoef) ) DEALLOCATE(petAcoef)
    IF ( ALLOCATED(peqAcoef) ) DEALLOCATE(peqAcoef)
    IF ( ALLOCATED(petBcoef) ) DEALLOCATE(petBcoef)
    IF ( ALLOCATED(peqBcoef) ) DEALLOCATE(peqBcoef)
    IF ( ALLOCATED(u_tq) ) DEALLOCATE(u_tq)
    IF ( ALLOCATED(v_tq) ) DEALLOCATE(v_tq)
    IF ( ALLOCATED(vevapp) ) DEALLOCATE(vevapp)
    IF ( ALLOCATED(fluxsens) ) DEALLOCATE(fluxsens)
    IF ( ALLOCATED(fluxlat) ) DEALLOCATE(fluxlat)
    IF ( ALLOCATED(coastalflow) ) DEALLOCATE(coastalflow)
    IF ( ALLOCATED(riverflow) ) DEALLOCATE(riverflow)
    IF ( ALLOCATED(netco2) ) DEALLOCATE(netco2)
    IF ( ALLOCATED(carblu) ) DEALLOCATE(carblu)
    IF ( ALLOCATED(carbwh) ) DEALLOCATE(carbwh)
    IF ( ALLOCATED(carbha) ) DEALLOCATE(carbha)
    IF ( ALLOCATED(tsol_rad) ) DEALLOCATE(tsol_rad)
    IF ( ALLOCATED(temp_sol_new) ) DEALLOCATE(temp_sol_new)
    IF ( ALLOCATED(qsurf) ) DEALLOCATE(qsurf)
    IF ( ALLOCATED(albedo) ) DEALLOCATE(albedo)
    IF ( ALLOCATED(emis) ) DEALLOCATE(emis)
    IF ( ALLOCATED(z0m) ) DEALLOCATE(z0m)
    IF ( ALLOCATED(z0h) ) DEALLOCATE(z0h)
    !
    WRITE(numout,*) "Memory deallocated"
    !
  END SUBROUTINE orchideedriver_clear
  !
END PROGRAM orchideedriver
