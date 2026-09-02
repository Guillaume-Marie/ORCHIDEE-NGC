MODULE routing_interp_topo_flow_mod
#ifdef XIOS
  USE ioipsl   
  USE xios_orchidee
  USE ioipsl_para 
  USE constantes
  USE constantes_soil
  USE pft_parameters
  USE sechiba_io_p
  USE interpol_help
  USE grid, not_used => contfrac
  USE function_library
  USE mod_orchidee_para


  IMPLICIT NONE
  PRIVATE
 
  PUBLIC :: routing_flow_xios_initialize, routing_flow_set, routing_flow_get, routing_flow_main
  PUBLIC :: routing_flow_initialize, routing_flow_finalize, routing_flow_clear, routing_flow_make_mean
  PUBLIC :: routing_flow_diags, routing_flow_get_water_balance, routing_flow_retreive_outflow
  PUBLIC :: compute_coastline
    
  INTEGER,SAVE :: nbpt
  !$OMP THREADPRIVATE(nbpt)
   
  !! PARAMETERS
  REAL(r_std), SAVE                         :: fast_tcst = 9e-3             !! Property of the fast reservoir (day/km)
  !$OMP THREADPRIVATE(fast_tcst)
  REAL(r_std), SAVE                         :: slow_tcst = 1.2e-2           !! Property of the slow reservoir (day/km)
  !$OMP THREADPRIVATE(slow_tcst)
  REAL(r_std), SAVE                         :: stream_tcst = 3e-5           !! Property of the stream reservoir (day/km)
  !$OMP THREADPRIVATE(stream_tcst)
  

  REAL(r_std),SAVE,ALLOCATABLE :: runoff_mean(:)               
  !$OMP THREADPRIVATE(runoff_mean)

  REAL(r_std),SAVE,ALLOCATABLE :: drainage_mean(:)               
  !$OMP THREADPRIVATE(drainage_mean)

  REAL(r_std),SAVE,ALLOCATABLE,PUBLIC,PROTECTED,TARGET :: fast_reservoir_r(:)         !! Water amount in the fast reservoir (kg) - (local routing grid)
  !$OMP THREADPRIVATE(fast_reservoir_r)

  REAL(r_std),SAVE,ALLOCATABLE,PUBLIC,PROTECTED,TARGET:: slow_reservoir_r(:)          !! Water amount in the slow reservoir (kg) - (local routing grid)
  !$OMP THREADPRIVATE(slow_reservoir_r)

  REAL(r_std),SAVE,ALLOCATABLE,PUBLIC,PROTECTED,TARGET :: stream_reservoir_r(:)       !! Water amount in the stream reservoir (kg) - (local routing grid)
  !$OMP THREADPRIVATE(stream_reservoir_r)

  
  REAL(r_std),SAVE,ALLOCATABLE,PUBLIC,PROTECTED,TARGET :: fast_diag(:)         !! Diag on water amount in the fast reservoir (kg/m^2) - (orchidee grid)
  !$OMP THREADPRIVATE(fast_diag)

  REAL(r_std),SAVE,ALLOCATABLE,PUBLIC,PROTECTED,TARGET:: slow_diag(:)          !! Diag on water amount in the slow reservoir (kg/m^2) - (orchidee grid)
  !$OMP THREADPRIVATE(slow_diag)

  REAL(r_std),SAVE,ALLOCATABLE,PUBLIC,PROTECTED,TARGET :: stream_diag(:)       !! Diag on water amount in the stream reservoir (kg/m^2) - (orchidee grid)
  !$OMP THREADPRIVATE(stream_diag) 

  REAL(r_std),SAVE,ALLOCATABLE,PUBLIC,PROTECTED,TARGET :: hydrographs_diag(:)       !! Diag on water amount in the stream reservoir (kg/m^2) - (orchidee grid)
  !$OMP THREADPRIVATE(hydrographs_diag) 

  
  REAL(r_std),SAVE,ALLOCATABLE,PUBLIC,PROTECTED :: riverflow_mean(:)         !! Water amount in the stream reservoir (kg) - (local routing grid)
  !$OMP THREADPRIVATE(riverflow_mean) 

  REAL(r_std),SAVE,ALLOCATABLE,PUBLIC,PROTECTED :: coastalflow_mean(:)       !! Water amount in the stream reservoir (kg) - (local routing grid)
  !$OMP THREADPRIVATE(coastalflow_mean) 

  REAL(r_std),SAVE,ALLOCATABLE,PUBLIC,PROTECTED :: lakeinflow_mean(:)        !! Water amount in the stream reservoir (kg) - (local routing grid)
  !$OMP THREADPRIVATE(lakeinflow_mean) 
  
  REAL(r_std),SAVE :: integration_time
  !$OMP THREADPRIVATE(integration_time)

  !
  ! Specific variables for interp_topo routing
  !
  REAL(r_std),SAVE,ALLOCATABLE :: topoind_r(:)                !! Topographic index of the retention time (m) index - (local routing grid)
  !$OMP THREADPRIVATE(topoind_r)
  REAL(r_std),SAVE             :: topoind_factor              !! conversion factor for topographic index (merit : m, standard : km), topo index must be in m
  !$OMP THREADPRIVATE(topoind_factor)
  INTEGER,SAVE,ALLOCATABLE     :: route_flow_rp1(:)           !! flow index from cell to neighboring cell following the trip direction - (local routing grid + halo)   
  !$OMP THREADPRIVATE(route_flow_rp1)
  LOGICAL,SAVE,ALLOCATABLE     :: is_lakeinflow_r(:)          !! is lake inflow point  - (local routing grid)
  !$OMP THREADPRIVATE(is_lakeinflow_r)
  LOGICAL,SAVE,ALLOCATABLE     :: is_coastalflow_r(:)         !! is coastal flow point - (local routing grid)
  !$OMP THREADPRIVATE(is_coastalflow_r)
  LOGICAL,SAVE,ALLOCATABLE     :: is_riverflow_r(:)           !! is river flow point - (local routing grid) 
  !$OMP THREADPRIVATE(is_riverflow_r)
  LOGICAL,SAVE,ALLOCATABLE     :: is_coastline(:)             !! is coastline point - (local native grid) 
  !$OMP THREADPRIVATE(is_coastline)
  LOGICAL,SAVE,ALLOCATABLE     :: is_streamflow_r(:)          !! is stream flow point - (local routing grid) 
  !$OMP THREADPRIVATE(is_streamflow_r)
  LOGICAL,SAVE,ALLOCATABLE,PUBLIC,PROTECTED    :: routing_mask_r(:)    !! valid routing point - (local routing grid) 
  !$OMP THREADPRIVATE(routing_mask_r)
  LOGICAL,SAVE,ALLOCATABLE     :: coast_mask(:)               !! is a coast point - (local native grid)
  !$OMP THREADPRIVATE(coast_mask)
  INTEGER,SAVE                 :: total_coast_points        !! global number of coast point - (local native grid)                   
  !$OMP THREADPRIVATE(total_coast_points)
  INTEGER,SAVE                 :: nbpt_r                      !! number of point in local routing grid
  !$OMP THREADPRIVATE(nbpt_r)
  INTEGER,SAVE                 :: nbpt_rp1                    !! number of point in local routing grid with halo of 1
  !$OMP THREADPRIVATE(nbpt_rp1)
  REAL(r_std),SAVE,ALLOCATABLE :: routing_weight(:)           !! Weight to transform runoff and drainage flux to water quantity in a conservative way (local native grid -> local routing grid)
  !$OMP THREADPRIVATE(routing_weight)
  REAL(r_std),SAVE,ALLOCATABLE :: routing_weight_in(:)         !! Weight to transform runoff and drainage flux to water quantity in a conservative way (local native grid -> local routing grid)
  !$OMP THREADPRIVATE(routing_weight_in)
  REAL(r_std),SAVE,ALLOCATABLE :: unrouted_weight(:)           !! Weight to transform runoff and drainage flux to water quantity in a conservative way (local native grid -> local routing grid)
  !$OMP THREADPRIVATE(unrouted_weight)

  REAL(r_std),SAVE,ALLOCATABLE :: weight_coast_to_coast_r(:)   
  !$OMP THREADPRIVATE(weight_coast_to_coast_r)
  
  REAL(r_std),SAVE,ALLOCATABLE :: weight_coast_to_lake_r(:)   
  !$OMP THREADPRIVATE(weight_coast_to_lake_r)
  
  REAL(r_std),SAVE,ALLOCATABLE :: weight_lake_to_coast_r(:)   
  !$OMP THREADPRIVATE(weight_lake_to_coast_r)
  
  REAL(r_std),SAVE,ALLOCATABLE :: weight_lake_to_lake_r(:)   
  !$OMP THREADPRIVATE(weight_lake_to_lake_r)

  REAL(r_std),SAVE,ALLOCATABLE :: basins_area_r(:)            !! upstream basin area of the routing cells (m2) - (local routing grid) 
  !$OMP THREADPRIVATE(basins_area_r)

  !! when doing interpolation from local routing grid to local native grid (river flow+coastal flow)
  REAL(r_std),SAVE,ALLOCATABLE :: basins_extended_r(:)        !! basins riverflow id (local routing grid) 
  !$OMP THREADPRIVATE(basins_extended_r)
  INTEGER(i_std)               :: basins_count                !! number of basins (local routing grid) 
  !$OMP THREADPRIVATE(basins_count)
  INTEGER(i_std)               :: basins_out                  !! number of basins to output for diag 
  !$OMP THREADPRIVATE(basins_out)

  INTEGER(i_std)               :: split_routing               !! time spliting for routing
  !$OMP THREADPRIVATE(split_routing)


  ! Variables for station hydrographs (diagnostics)
  INTEGER(i_std),SAVE                :: nb_station            !! number of stations
  CHARACTER(LEN=60),SAVE,ALLOCATABLE :: station(:)            !! station names
  INTEGER(i_std),SAVE,ALLOCATABLE    :: station_index(:)      !! station index on the routing grid
  INTEGER(i_std),SAVE                :: station_ts = 0        !! timestep index for station output file
  
  INTEGER,SAVE :: printlev_loc
  !$OMP THREADPRIVATE(printlev_loc)

!  INTEGER(i_std), PARAMETER :: nb_stations=14
!  REAL(r_std),PARAMETER  :: station_lon(nb_stations) = &
!       (/ -162.8830, -90.9058, -55.5110, -49.3242, -133.7447, -63.6000,  28.7167, &
!       15.3000,  66.5300,  89.6700,  86.5000,  127.6500,   3.3833, 117.6200 /)
!  REAL(r_std),PARAMETER  :: station_lat(nb_stations) = &
!       (/  61.9340,  32.3150, -1.9470, -5.1281, 67.4583,  8.1500, 45.2167, &
!       -4.3000,  66.5700, 25.1800, 67.4800, 70.7000, 11.8667, 30.7700 /)
!  CHARACTER(LEN=17),PARAMETER  :: station_name(nb_stations) = &
!       (/ "Pilot station    ", "Vicksburg        ", "Obidos           ", &
!       "Itupiranga       ", "Arctic red river ", "Puente Angostura ", &
!       "Ceatal Izmail    ", "Kinshasa         ", "Salekhard        ", &
!       "Bahadurabad      ", "Igarka           ", "Kusur            ", &
!       "Malanville       ", "Datong           " /)

CONTAINS


  SUBROUTINE routing_flow_get(slow_reservoir_r, fast_reservoir_r, stream_reservoir_r, riverflow_mean, coastalflow_mean, &
                              lakeinflow_mean, runoff_mean, drainage_mean)
  IMPLICIT NONE
    REAL(r_std),OPTIONAL, INTENT(OUT) :: slow_reservoir_r(:)
    REAL(r_std),OPTIONAL, INTENT(OUT) :: fast_reservoir_r(:)
    REAL(r_std),OPTIONAL, INTENT(OUT) :: stream_reservoir_r(:)
    REAL(r_std),OPTIONAL, INTENT(OUT) :: riverflow_mean(:)
    REAL(r_std),OPTIONAL, INTENT(OUT) :: coastalflow_mean(:)
    REAL(r_std),OPTIONAL, INTENT(OUT) :: lakeinflow_mean(:)
    REAL(r_std),OPTIONAL, INTENT(OUT) :: runoff_mean(:)
    REAL(r_std),OPTIONAL, INTENT(OUT) :: drainage_mean(:)

    CALL routing_flow_get_(slow_reservoir_r, fast_reservoir_r, stream_reservoir_r, riverflow_mean, coastalflow_mean, lakeinflow_mean, &
                           runoff_mean, drainage_mean)
  END SUBROUTINE routing_flow_get

  SUBROUTINE routing_flow_get_(slow_reservoir_r_, fast_reservoir_r_, stream_reservoir_r_, riverflow_mean_, &
                               coastalflow_mean_, lakeinflow_mean_, runoff_mean_, drainage_mean_)
  IMPLICIT NONE
    REAL(r_std),OPTIONAL, INTENT(OUT) :: slow_reservoir_r_(:)
    REAL(r_std),OPTIONAL, INTENT(OUT) :: fast_reservoir_r_(:)
    REAL(r_std),OPTIONAL, INTENT(OUT) :: stream_reservoir_r_(:)
    REAL(r_std),OPTIONAL, INTENT(OUT) :: riverflow_mean_(:)
    REAL(r_std),OPTIONAL, INTENT(OUT) :: coastalflow_mean_(:)
    REAL(r_std),OPTIONAL, INTENT(OUT) :: lakeinflow_mean_(:)
    REAL(r_std),OPTIONAL, INTENT(OUT) :: runoff_mean_(:)
    REAL(r_std),OPTIONAL, INTENT(OUT) :: drainage_mean_(:)

    IF (PRESENT(slow_reservoir_r_))   slow_reservoir_r_ = slow_reservoir_r
    IF (PRESENT(fast_reservoir_r_))   fast_reservoir_r_ = fast_reservoir_r
    IF (PRESENT(stream_reservoir_r_)) stream_reservoir_r_ = stream_reservoir_r
    IF (PRESENT(riverflow_mean_))     riverflow_mean_ = riverflow_mean
    IF (PRESENT(coastalflow_mean_))   coastalflow_mean_ = coastalflow_mean
    IF (PRESENT(lakeinflow_mean_))    lakeinflow_mean_ = lakeinflow_mean
    IF (PRESENT(runoff_mean_))        runoff_mean_ = runoff_mean
    IF (PRESENT(drainage_mean_))      drainage_mean_ = drainage_mean

  END SUBROUTINE routing_flow_get_

  SUBROUTINE routing_flow_set(slow_reservoir_r, fast_reservoir_r, stream_reservoir_r, riverflow_mean, coastalflow_mean, &
                             lakeinflow_mean, runoff_mean, drainage_mean)
  IMPLICIT NONE
    REAL(r_std),OPTIONAL, INTENT(IN) :: slow_reservoir_r(:)
    REAL(r_std),OPTIONAL, INTENT(IN) :: fast_reservoir_r(:)
    REAL(r_std),OPTIONAL, INTENT(IN) :: stream_reservoir_r(:)
    REAL(r_std),OPTIONAL, INTENT(IN) :: riverflow_mean(:)
    REAL(r_std),OPTIONAL, INTENT(IN) :: coastalflow_mean(:)
    REAL(r_std),OPTIONAL, INTENT(IN) :: lakeinflow_mean(:)
    REAL(r_std),OPTIONAL, INTENT(IN) :: runoff_mean(:)
    REAL(r_std),OPTIONAL, INTENT(IN) :: drainage_mean(:)
    
    CALL routing_flow_set_(slow_reservoir_r, fast_reservoir_r, stream_reservoir_r, riverflow_mean, coastalflow_mean, lakeinflow_mean)

  END SUBROUTINE routing_flow_set

  SUBROUTINE routing_flow_set_(slow_reservoir_r_, fast_reservoir_r_, stream_reservoir_r_, riverflow_mean_, &
                               coastalflow_mean_, lakeinflow_mean_, runoff_mean_, drainage_mean_)
  IMPLICIT NONE
    REAL(r_std),OPTIONAL, INTENT(IN) :: slow_reservoir_r_(:)
    REAL(r_std),OPTIONAL, INTENT(IN) :: fast_reservoir_r_(:)
    REAL(r_std),OPTIONAL, INTENT(IN) :: stream_reservoir_r_(:)
    REAL(r_std),OPTIONAL, INTENT(IN) :: riverflow_mean_(:)
    REAL(r_std),OPTIONAL, INTENT(IN) :: coastalflow_mean_(:)
    REAL(r_std),OPTIONAL, INTENT(IN) :: lakeinflow_mean_(:)
    REAL(r_std),OPTIONAL, INTENT(IN) :: runoff_mean_(:)
    REAL(r_std),OPTIONAL, INTENT(IN) :: drainage_mean_(:)

    IF (PRESENT(slow_reservoir_r_))   slow_reservoir_r = slow_reservoir_r_
    IF (PRESENT(fast_reservoir_r_))   fast_reservoir_r = fast_reservoir_r_
    IF (PRESENT(stream_reservoir_r_)) stream_reservoir_r = stream_reservoir_r_
    IF (PRESENT(riverflow_mean_))     riverflow_mean = riverflow_mean_
    IF (PRESENT(coastalflow_mean_))   coastalflow_mean = coastalflow_mean_
    IF (PRESENT(lakeinflow_mean_))     lakeinflow_mean = lakeinflow_mean_
    IF (PRESENT(runoff_mean_))        runoff_mean = runoff_mean_
    IF (PRESENT(drainage_mean_))      drainage_mean = drainage_mean_

  END SUBROUTINE routing_flow_set_  
!!  =============================================================================================================================
!! SUBROUTINE:    routing_flow_xios_initialize
!!
!>\BRIEF	  Initialize xios dependant definition before closing context definition
!!
!! DESCRIPTION:	  Initialize xios dependant definition before closing context definition.
!!                This subroutine is called before the xios context is closed. 
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S): None
!! 
!! FLOWCHART: None
!! \n
!_ ==============================================================================================================================

  SUBROUTINE routing_flow_xios_initialize
    USE xios
    USE routing_subgrid_halfdeg, ONLY : routing_names
    IMPLICIT NONE     

    INTEGER(i_std) ::ib
    INTEGER :: nbasmax=1
    !! 0 Variable and parameter description
    CHARACTER(LEN=60),ALLOCATABLE :: label(:)
    LOGICAL :: file_exists

    IF (is_omp_root) THEN   
       CALL xios_get_axis_attr("basins", n_glo=basins_out)    ! get nb basins to output
       ALLOCATE(label(basins_out))
       CALL routing_names(basins_out,label)
       CALL xios_set_axis_attr("basins", label=label)         ! set riverflow basins name
       INQUIRE(FILE="routing_start.nc", EXIST=file_exists)
       IF (file_exists) CALL xios_set_file_attr("routing_start", enabled=.TRUE.)  
    ENDIF

    !! Define XIOS axis size needed for the model output
    ! Add axis for homogeneity between all routing schemes, these dimensions are currently not used in this scheme
    CALL xios_orchidee_addaxis("nbhtu", nbasmax, (/(REAL(ib,r_std),ib=1,nbasmax)/))
    CALL xios_orchidee_addaxis("nbasmon", 1, (/(REAL(ib,r_std),ib=1,1)/))

  END SUBROUTINE routing_flow_xios_initialize


  SUBROUTINE routing_flow_initialize(kjit, rest_id, nbpt_, dt_routing, contfrac, nbpt_r_, riverflow, coastalflow)
  USE grid, ONLY : area
  USE xios
  IMPLICIT NONE
    INTEGER ,INTENT(IN)                     :: kjit           
    INTEGER ,INTENT(IN)                     :: rest_id           
    INTEGER, INTENT(IN)                     :: nbpt_                 !! nb points orchidee grid
    REAL(r_std), INTENT(IN)                 :: dt_routing                       
    REAL(r_std), INTENT(IN)                 :: contfrac(nbpt_)       !! fraction of land
    INTEGER, INTENT(OUT)                    :: nbpt_r_               !! nb points routing grid
    REAL(r_std), INTENT(OUT)                :: riverflow(nbpt_)      !! Outflow of the major rivers. The flux will be located on the continental grid but this should be a coastal point (kg/dt)
    REAL(r_std), INTENT(OUT)                :: coastalflow(nbpt_)    !! Outflow on coastal points by small basins. This is the water which flows in a disperse way into the ocean (kg/dt)
    INTEGER :: ier    
    CHARACTER(LEN=80)   :: var_name       !! To store variables names for I/O (unitless)
     
    nbpt = nbpt_ 
    CALL routing_flow_init_local(contfrac, nbpt_r_, dt_routing)
    nbpt_r = nbpt_r_
    CALL routing_flow_init_mean(kjit, rest_id)
    CALL initialize_stations(dt_routing)
    !
    ! Put into the restart file the fluxes so that they can be regenerated at restart.
    !
    ALLOCATE (lakeinflow_mean(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_flow_initialize','Pb in allocate for lakeinflow_mean','','')
    var_name = 'lakeinflow'
    CALL ioconf_setatt_p('UNITS', 'Kg/dt')
    CALL ioconf_setatt_p('LONG_NAME','Lake inflow')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., lakeinflow_mean, "gather", nbp_glo, index_g)
    CALL setvar_p (lakeinflow_mean, val_exp, 'NO_KEYWORD', zero)

    ALLOCATE (riverflow_mean(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_flow_initialize','Pb in allocate for riverflow_mean','','')
    var_name = 'riverflow'
    CALL ioconf_setatt_p('UNITS', 'Kg/dt')
    CALL ioconf_setatt_p('LONG_NAME','River flux into the sea')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., riverflow_mean, "gather", nbp_glo, index_g)
    CALL setvar_p (riverflow_mean, val_exp, 'NO_KEYWORD', zero)

    ALLOCATE (coastalflow_mean(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_flow_initialize','Pb in allocate for coastalflow_mean','','')
    var_name = 'coastalflow'
    CALL ioconf_setatt_p('UNITS', 'Kg/dt')
    CALL ioconf_setatt_p('LONG_NAME','Diffuse flux into the sea')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., coastalflow_mean, "gather", nbp_glo, index_g)
    CALL setvar_p (coastalflow_mean, val_exp, 'NO_KEYWORD', zero)
    
    riverflow(:)   =  riverflow_mean(:)
    coastalflow(:) =  coastalflow_mean(:)

    CALL routing_flow_initialize_diag(kjit, rest_id, contfrac)

  END SUBROUTINE routing_flow_initialize

  SUBROUTINE routing_flow_initialize_diag(kjit, rest_id, contfrac)
    USE xios
    USE grid, ONLY : area 
    IMPLICIT NONE
    INTEGER ,INTENT(IN)                     :: kjit           
    INTEGER ,INTENT(IN)                     :: rest_id           
    REAL(r_std), INTENT(IN)                 :: contfrac(nbpt)       !! fraction of land
    INTEGER :: ier    
    CHARACTER(LEN=80)   :: var_name       !! To store variables names for I/O (unitless)
     
    
    REAL(r_std) :: contfrac_mpi(nbp_mpi)
    REAL(r_std) :: area_mpi(nbp_mpi)
    REAL(r_std) :: fast_diag_mpi(nbp_mpi)
    REAL(r_std) :: slow_diag_mpi(nbp_mpi)
    REAL(r_std) :: stream_diag_mpi(nbp_mpi)
    REAL(r_std) :: hydrographs_diag_mpi(nbp_mpi)

    CALL gather_omp(contfrac, contfrac_mpi)
    CALL gather_omp(area, area_mpi)
  
    ALLOCATE(fast_diag(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_flow_initialize_diag','Pb in allocate for fast_diag','','')
    ALLOCATE(slow_diag(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_flow_initialize_diag','Pb in allocate for slow_diag','','')
    ALLOCATE(stream_diag(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_flow_initialize_diag','Pb in allocate for stream_diag','','')

    IF (is_omp_root) THEN
      CALL xios_send_field("routing_fast_diag0_r", fast_reservoir_r)
      CALL xios_recv_field("routing_fast_diag0", fast_diag_mpi)
      CALL xios_send_field("routing_slow_diag0_r", slow_reservoir_r)
      CALL xios_recv_field("routing_slow_diag0", slow_diag_mpi)
      CALL xios_send_field("routing_stream_diag0_r", stream_reservoir_r)
      CALL xios_recv_field("routing_stream_diag0", stream_diag_mpi)

      fast_diag_mpi=fast_diag_mpi/(area_mpi*contfrac_mpi)      !! kg => kg/m^2
      slow_diag_mpi=slow_diag_mpi/(area_mpi*contfrac_mpi)      !! kg => kg/m^2
      stream_diag_mpi=stream_diag_mpi/(area_mpi*contfrac_mpi)  !! kg => kg/m^2
    ENDIF
    
    CALL scatter_omp(fast_diag_mpi,fast_diag)
    CALL scatter_omp(slow_diag_mpi,slow_diag)
    CALL scatter_omp(stream_diag_mpi,stream_diag)
    
    ALLOCATE(hydrographs_diag(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_flow_initialize_diag','Pb in allocate for hydrographs_diag','','')
 
    var_name = 'hydrographs'
    CALL ioconf_setatt_p('UNITS', 'kg/dt_sechiba')
    CALL ioconf_setatt_p('LONG_NAME','Hydrograph at outlow of grid')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., hydrographs_diag, "gather", nbp_glo, index_g)
    CALL setvar_p (hydrographs_diag, val_exp, 'NO_KEYWORD', zero)

  END SUBROUTINE routing_flow_initialize_diag


  SUBROUTINE compute_coastline(contfrac, coastline)
  USE mod_orchidee_para
  USE grid, not_used => contfrac
  IMPLICIT NONE
    REAL(r_std),INTENT(IN)     :: contfrac(nbpt)    ! INPUT   : fraction of continent (unitless)
    LOGICAL,INTENT(OUT)        :: coastline(nbpt)   ! OUTPUT  : coastline mask (true : on coastaline, else false
    
    REAL(r_std) :: contfrac_glo(nbp_glo)
    REAL(r_std) :: contfrac2D_glo(iim_g*jjm_g)
    LOGICAL     :: coastline_glo(nbp_glo)
    LOGICAL     :: coastline2D_glo(iim_g*jjm_g)
    LOGICAL     :: mask1d(nbpt)
    LOGICAL     :: mask1d_glo(nbp_glo)
    LOGICAL     :: mask2d_glo(iim_g*jjm_g)
    INTEGER(i_std) :: i,j,ij,next_i,next_j,next_ij,m,k
    REAL(r_std), PARAMETER :: epsilon=1e-5
    LOGICAL :: is_periodic
    
    is_periodic = global
    mask1d=.TRUE. ;
    CALL gather(mask1d,mask1d_glo)
    CALL gather(contfrac,contfrac_glo)
    
    
    IF (is_mpi_root .AND. is_omp_root) THEN
      contfrac2D_glo(:)=0
      mask2d_glo=.FALSE.
      DO i=1,nbp_glo
        contfrac2D_glo(index_g(i)) = contfrac_glo(i)
        mask2d_glo(index_g(i))=mask1d_glo(i)
      ENDDO

      coastline2D_glo(:)=.FALSE.
      DO j=1,jjm_g
        DO i=1,iim_g
          ij=(j-1)*iim_g+i
          IF (grid_type==unstructured) THEN
            IF (contfrac2D_glo(ij)>epsilon .AND. contfrac2D_glo(ij)<1-epsilon) THEN
              coastline2D_glo(ij)=.TRUE.
            ENDIF
          ELSE
            DO k=-1,1
              DO m=-1,1
                IF (k==0 .AND. m==0) CYCLE
                next_i=i+k
                next_j=j+m
                
                IF (next_i==0) THEN ! manage periodicity
                  IF (is_periodic) THEN
                    next_i=iim_g
                  ELSE 
                    coastline2D_glo(ij)=.TRUE.
                    CYCLE
                  ENDIF
                ENDIF 

                IF (next_i==iim_g+1) THEN! manage periodicity
                  IF (is_periodic) THEN
                    next_i=1 !! periodic but not true for limited area 
                  ELSE
                    coastline2D_glo(ij)=.TRUE.
                    CYCLE
                  ENDIF
                ENDIF
                IF (next_j==0 .OR. next_j==jjm_g+1) THEN
                  IF (.NOT. is_periodic) coastline2D_glo(ij)=.TRUE.
                  CYCLE
                ENDIF
                next_ij =  (next_j-1)*iim_g+next_i  
                IF (.NOT. mask2d_glo(next_ij)) coastline2D_glo(ij)=.TRUE.
    !            IF (contfrac2D_glo(next_ij)<=epsilon) coastline2D_glo(ij)=.TRUE.
              ENDDO
            ENDDO
          ENDIF    
        ENDDO
      ENDDO
    
      DO i=1,nbp_glo
        coastline_glo(i) = coastline2D_glo(index_g(i))
      ENDDO

    ENDIF
    
    CALL scatter(coastline_glo, coastline)
    
  END SUBROUTINE compute_coastline

 
 SUBROUTINE routing_flow_correct_riverflow(ni, nj, contfrac, coastline, trip_r, trip_extended_r, topoind_r)
  USE xios
  USE grid, ONLY : global
  IMPLICIT NONE
    INCLUDE "mpif.h"
    INTEGER, INTENT(IN)        :: ni                      ! INPUT : size i (longitude) of the local routing domain
    INTEGER, INTENT(IN)        :: nj                      ! INPUT : size i (latitude)  of the local routing domain
    REAL(r_std), INTENT(IN)    :: contfrac(nbp_mpi)       ! INPUT : continental fraction (unitless)
    LOGICAL, INTENT(IN)        :: coastline(nbp_mpi)      ! INPUT/OUTPUT : coastline mask
    REAL(r_std), INTENT(INOUT) :: trip_r(ni,nj)           ! INPUT/OUTPUT : diection of flow, which will be modified by the routine
    REAL(r_std), INTENT(INOUT)    :: trip_extended_r(ni,nj)  ! INPUT  : direction of flow extended to ocean points
    REAL(r_std), INTENT(INOUT) :: topoind_r(ni,nj)        ! INPUT/OUTPUT : topographic index which will be modified by the routine if extended to ocean
    
    REAL(r_std)                :: mask(nbp_mpi)
    REAL(r_std)                :: contfrac_r(ni,nj)
    REAL(r_std)                :: mask_coast(nbp_mpi)
    REAL(r_std)                :: frac_coast_r(ni,nj)
    REAL(r_std)                :: trip_extended_rp1(0:ni+1,0:nj+1)
    REAL(r_std)                :: trip_rp1(0:ni+1,0:nj+1)
    REAL(r_std)                :: state_r(ni,nj)
    REAL(r_std)                :: state_rp1(0:ni+1,0:nj+1)
    REAL(r_std)                :: done_rp1(0:ni+1,0:nj+1)
    REAL(r_std)                :: trip_out_r(0:ni+1,0:nj+1)

    INTEGER, PARAMETER ::       is_ter=0
    INTEGER, PARAMETER ::       is_coast=1
    INTEGER, PARAMETER ::       is_oce=2
    INTEGER, PARAMETER ::       riverflow=99
    INTEGER, PARAMETER ::       coastalflow=98
    INTEGER, PARAMETER ::       lakeinflow=97
    INTEGER, PARAMETER ::       norouting=100

    INTEGER            :: updated
    INTEGER            :: it, i,j,k,m, nexti,nextj,previ,prevj,ierr
    INTEGER            :: ibegin, jbegin, ni_glo, nj_glo
    REAL(r_std)        :: lonvalue_1d(ni), latvalue_1d(nj)
    REAL(r_std),PARAMETER :: epsilon=1e-5
 
    TYPE(xios_duration) :: ts 
    TYPE(xios_domain) :: domain_hdl
    TYPE(xios_domaingroup) :: domain_def_hdl
    TYPE(xios_field) :: field_hdl
    TYPE(xios_fieldgroup) :: field_def_hdl
    TYPE(xios_file) :: file_hdl
    TYPE(xios_filegroup) :: file_def_hdl
    TYPE(xios_expand_domain) :: domain_expand_hdl
    TYPE(xios_date) :: start_date, time_origin
    CHARACTER(LEN=20)  :: calendar_type
    LOGICAL :: ok

      CALL xios_get_domain_attr("routing_domain", ibegin=ibegin, jbegin=jbegin, ni_glo=ni_glo,nj_glo=nj_glo)    ! get routing domain dimension
      CALL xios_get_domain_attr("routing_domain", lonvalue_1d=lonvalue_1d, latvalue_1d=latvalue_1d)    ! get routing domain dimension
      
      ! manage non periodic boundary in case of LAM
      IF (.NOT. global) THEN
        IF (ibegin==0) THEN
          i=1 ;
          DO j=1,nj
            IF (trip_r(i,j)>=1 .AND. trip_r(i,j)<=8) trip_r(i,j) = coastalflow
            IF (trip_extended_r(i,j)>=1 .AND. trip_extended_r(i,j)<=8) trip_extended_r(i,j) = coastalflow
          ENDDO
        ENDIF

        IF (ibegin+ni==ni_glo) THEN
          i=ni ;
          DO j=1,nj
            IF (trip_r(i,j)>=1 .AND. trip_r(i,j)<=8) trip_r(i,j) = coastalflow
            IF (trip_extended_r(i,j)>=1 .AND. trip_extended_r(i,j)<=8) trip_extended_r(i,j) = coastalflow
          ENDDO
        ENDIF
        
        IF (jbegin==0) THEN
          j=1 ;
          DO i=1,ni
            IF (trip_r(i,j)>=1 .AND. trip_r(i,j)<=8) trip_r(i,j) = coastalflow
            IF (trip_extended_r(i,j)>=1 .AND. trip_extended_r(i,j)<=8) trip_extended_r(i,j) = coastalflow
          ENDDO
        ENDIF

        IF (jbegin+nj==nj_glo) THEN
          j=nj ;
          DO i=1,ni
            IF (trip_r(i,j)>=1 .AND. trip_r(i,j)<=8) trip_r(i,j) = coastalflow
            IF (trip_extended_r(i,j)>=1 .AND. trip_extended_r(i,j)<=8) trip_extended_r(i,j) = coastalflow
          ENDDO
        ENDIF
      ENDIF

      mask = 1
      CALL xios_send_field("routing_contfrac", mask) 
      CALL xios_recv_field("routing_contfrac_r", contfrac_r)

      CALL xios_send_field("trip_ext_r",trip_extended_r)
      CALL xios_recv_field("trip_ext_rp1",trip_extended_rp1)

      mask_coast=0
      WHERE(coastline) mask_coast=1
      CALL xios_send_field("mask_coastline",mask_coast)
      CALL xios_recv_field("frac_coastline_r",frac_coast_r)
      
      DO j=1,nj
        DO i=1,ni
          IF (frac_coast_r(i,j)>epsilon) THEN
            state_r(i,j)=is_coast
          ELSE 
            IF (contfrac_r(i,j)> epsilon) THEN 
             state_r(i,j)=is_ter
            ELSE
              state_r(i,j)=is_oce
            ENDIF
          ENDIF
        ENDDO
      ENDDO


      DO j=1,nj
        DO i=1,ni
          IF (trip_r(i,j)>=100) THEN
            trip_r(i,j)=norouting
          ENDIF
        ENDDO
      ENDDO          

      CALL xios_context_initialize("orchidee_init_routing",MPI_COMM_ORCH)
      CALL xios_orchidee_change_context("orchidee_init_routing")
      CALL xios_set_domain_attr("routing_domain", ibegin=ibegin, jbegin=jbegin, ni=ni, nj=nj, ni_glo=ni_glo, nj_glo=nj_glo)    ! get routing domain dimension
      CALL xios_set_domain_attr("routing_domain", lonvalue_1d=lonvalue_1d, latvalue_1d=latvalue_1d)                           ! get routing domain dimension
      CALL xios_close_context_definition()
      CALL xios_update_calendar(1)

      CALL xios_send_field("trip_r_init",trip_r)
      CALL xios_send_field("topoind_r_init",topoind_r)
      
      CALL xios_send_field("trip_r",trip_r)
      CALL xios_recv_field("trip_rp1",trip_rp1)
     
      CALL xios_send_field("state_r",state_r)
      CALL xios_recv_field("state_rp1",state_rp1)
      
!      DO j=1,nj
!        DO i=1,ni
!          IF (trip_rp1(i,j)<=8) THEN
!            CALL next_trip(trip_rp1, ni, nj, i, j, nexti, nextj)
!            IF (trip_rp1(nexti,nextj)>=100) trip_rp1(i,j)=98 ! unrouted point becomes coastal flow
!          ENDIF
!        ENDDO
!      ENDDO          

      state_r(1:ni,1:nj) = state_rp1(1:ni,1:nj)
      trip_r(1:ni,1:nj)  = trip_rp1(1:ni,1:nj)

      updated=1
      it=2
      DO WHILE(updated>0)
        updated=0

        CALL xios_update_calendar(it)
        
        CALL xios_send_field("trip_r", trip_r)
        CALL xios_recv_field("trip_rp1",trip_rp1)
     
        CALL xios_send_field("state_r",state_r)
        CALL xios_recv_field("state_rp1",state_rp1)

        state_r(1:ni,1:nj) = state_rp1(1:ni,1:nj)
        trip_r(1:ni,1:nj)  = trip_rp1(1:ni,1:nj)
        
        ! first => riverflow and coastalflow routing to ocean must be routed to coast => to long 
        DO j=1,nj
          DO i=1,ni
            CALL next_trip(trip_rp1, ni, nj, i, j, nexti, nextj, ok)
            IF (ok) THEN
              IF (state_rp1(nexti,nextj)==is_oce .AND. (trip_rp1(nexti,nextj)==riverflow .OR. trip_rp1(nexti,nextj)==coastalflow .OR. trip_rp1(nexti,nextj)==lakeinflow)) THEN
                trip_r(i,j) = trip_rp1(nexti,nextj)
                updated=updated + 1
              ENDIF
            ENDIF
          ENDDO
        ENDDO

        DO j=1,nj
          DO i=1,ni
            IF (state_rp1(i,j)==is_oce) THEN
              IF (trip_rp1(i,j)==riverflow .OR. trip_rp1(i,j)==coastalflow .OR. trip_rp1(i,j)==lakeinflow) THEN
                trip_r(i,j) = norouting
                updated=updated+1
              ENDIF 
            ENDIF
          ENDDO
        ENDDO     

        CALL MPI_ALLREDUCE(MPI_IN_PLACE , updated, 1, MPI_INT_ORCH, MPI_SUM, MPI_COMM_ORCH, ierr)
        it = it +1
      ENDDO

      ! second prolongate riverflow and coastalflow to coast if they are in land 
      updated=1
      DO WHILE(updated>0)
        updated=0

        CALL xios_update_calendar(it)
        
        CALL xios_send_field("trip_r", trip_r)
        CALL xios_recv_field("trip_rp1",trip_rp1)
     
        CALL xios_send_field("state_r",state_r)
        CALL xios_recv_field("state_rp1",state_rp1)

        state_r(1:ni,1:nj) = state_rp1(1:ni,1:nj)
        trip_r(1:ni,1:nj)  = trip_rp1(1:ni,1:nj)
        
        DO j=1,nj
          DO i=1,ni
            IF (state_rp1(i,j)==is_ter .AND. (trip_rp1(i,j)==riverflow .OR. trip_rp1(i,j)==coastalflow .OR. trip_rp1(i,j)==lakeinflow)) THEN
              CALL next_trip(trip_extended_rp1, ni, nj, i, j, nexti, nextj, ok)
              IF (ok) THEN
                IF (state_rp1(nexti,nextj)==is_ter .OR. state_rp1(nexti,nextj)==is_coast) THEN
                  trip_r(i,j) = trip_extended_rp1(i,j)
                  topoind_r(i,j) = 1e-10  ! flow instantaneously for now but better to take the topoind of the incoming flux
                  updated=updated+1
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDDO

        DO j=1,nj
          DO i=1,ni
            IF (state_rp1(i,j)==is_ter .OR. state_rp1(i,j)==is_coast) THEN
              DO k=-1,1
                DO m=-1,1
                  IF (k==0 .AND. m==0) CYCLE
                  prevj = j+k ; previ = i+m
                  CALL next_trip(trip_extended_rp1, ni, nj, previ, prevj, nexti, nextj, ok)
                  IF (ok) THEN
                    IF (nexti==i .AND. nextj==j) THEN
                      IF (state_rp1(previ,prevj)==is_ter .AND. (trip_rp1(previ,prevj)==riverflow .OR. trip_rp1(previ,prevj)==coastalflow .OR. trip_rp1(previ,prevj)==lakeinflow)) THEN
                        IF (trip_r(i,j)/=riverflow) THEN
                          trip_r(i,j)=trip_rp1(previ,prevj)
                        ELSE 
                          trip_r(i,j)=riverflow
                        ENDIF
                        updated=updated+1
                      ENDIF
                    ENDIF
                  ENDIF
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDDO
        
        CALL MPI_ALLREDUCE(MPI_IN_PLACE , updated, 1, MPI_INT_ORCH, MPI_SUM, MPI_COMM_ORCH, ierr)
        it = it +1
      ENDDO

      CALL xios_update_calendar(it)
        
      CALL xios_send_field("trip_r", trip_r)
      CALL xios_recv_field("trip_rp1",trip_rp1)
     
      CALL xios_send_field("state_r",state_r)
      CALL xios_recv_field("state_rp1",state_rp1)      
      DO j=1,nj
        DO i=1,ni
          IF (state_rp1(i,j)==is_coast) THEN
            IF (trip_r(i,j)==norouting) trip_r(i,j) = coastalflow
          ENDIF
        ENDDO
      ENDDO
      
!      it = it +1
!      CALL xios_update_calendar(it)
!        
!      CALL xios_send_field("trip_r", trip_r)
!      CALL xios_recv_field("trip_rp1",trip_rp1)
!     
!      CALL xios_send_field("state_r",state_r)
!      CALL xios_recv_field("state_rp1",state_rp1)
!
!      DO j=1,nj
!        DO i=1,ni
!          IF (state_rp1(i,j)==is_oce) THEN
!            trip_r(i,j) = norouting
!          ENDIF
!        ENDDO
!      ENDDO
      
      it = it +1
      CALL xios_update_calendar(it)
        
      CALL xios_send_field("trip_r", trip_r)
      CALL xios_recv_field("trip_rp1",trip_rp1)
     
      CALL xios_send_field("state_r",state_r)
      CALL xios_recv_field("state_rp1",state_rp1)
! check
      updated=0
      DO j=1,nj
        DO i=1,ni
          IF (state_rp1(i,j)==is_ter) THEN
            IF (trip_r(i,j) == coastalflow .OR. trip_r(i,j) == riverflow .OR. trip_r(i,j) == norouting) THEN
              updated=updated+1
              trip_r(i,j) = lakeinflow
            ENDIF
          ELSE IF (state_rp1(i,j)==is_coast) THEN
            IF (trip_r(i,j) == norouting) THEN
              updated=updated+1
              trip_r(i,j) = coastalflow
            ENDIF
          ENDIF
        ENDDO
      ENDDO


      CALL MPI_ALLREDUCE(MPI_IN_PLACE , updated, 1, MPI_INT_ORCH, MPI_SUM, MPI_COMM_ORCH, ierr)
      IF (updated/=0 .AND. is_mpi_root)  WRITE (numout,*) "WARNING : ",updated," riverflow or coastal flow are not on the coast, ie in middle of land ==> converted into lakeinflow"
      CALL MPI_BARRIER(MPI_COMM_ORCH, ierr)

      updated=0
      DO j=1,nj
        DO i=1,ni
          CALL next_trip(trip_rp1,ni,nj,i,j,nexti,nextj,ok)
          IF (ok) THEN
            IF (trip_rp1(nexti,nextj)==norouting) THEN
              updated=updated+1
              WRITE (numout,*) "point",ibegin+i,jbegin+j," (trip(i,j)=",trip_rp1(i,j),") is routed to nothing => ",ibegin+nexti,jbegin+nextj, &
                     " (trip(i,j)=",trip_rp1(nexti,nextj),")"
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      
      CALL MPI_BARRIER(MPI_COMM_ORCH, ierr)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE , updated, 1, MPI_INT_ORCH, MPI_SUM, MPI_COMM_ORCH, ierr)
      IF (updated/=0) THEN
        IF (is_mpi_root) WRITE (numout,*) "ERROR : ",updated," points are not correctly routed (routed to nothing)"
        CALL xios_context_finalize
        STOP
      ENDIF


      updated=0
      DO j=1,nj
        DO i=1,ni
          CALL next_trip(trip_rp1,ni,nj,i,j,nexti,nextj,ok)
          IF (ok) THEN
            IF (state_rp1(nexti,nextj)==is_oce .AND. (trip_rp1(nexti,nextj)==riverflow .OR. trip_rp1(nexti,nextj)==coastalflow .OR. trip_rp1(nexti,nextj)==lakeinflow )) THEN
              updated=updated+1
              WRITE (numout,*) "point",ibegin+i,jbegin+j," (trip(i,j)=",trip_rp1(i,j),") is routed to coastalflow/riverflow/lakeinflow  point that is not on land grid => ",ibegin+nexti,jbegin+nextj, &
                     " (trip(i,j)=",trip_rp1(nexti,nextj),")"
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      
      CALL MPI_BARRIER(MPI_COMM_ORCH, ierr)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE , updated, 1, MPI_INT_ORCH, MPI_SUM, MPI_COMM_ORCH, ierr)
      IF (updated/=0) THEN
        IF (is_mpi_root) WRITE (numout,*) "ERROR : ",updated," points are routed to coastalflow/riverflow/lakeinflow points that are not on land grid"
        CALL xios_context_finalize
        STOP
      ENDIF



      it = it +1
      CALL xios_update_calendar(it)
        
      CALL xios_send_field("trip_r", trip_r)
      CALL xios_recv_field("trip_rp1",trip_rp1)
     
      CALL xios_send_field("state_r",state_r)
      CALL xios_recv_field("state_rp1",state_rp1)
      state_r(1:ni,1:nj) = state_rp1(1:ni,1:nj)
      trip_r(1:ni,1:nj)  = trip_rp1(1:ni,1:nj)

      CALL xios_send_field("trip_r_final",trip_r)
      
      ! Check if routing cells have bad topo_index and fix it
      updated=0
      DO j=1,nj
        DO i=1,ni
          IF (topoind_r(i,j) <= 0.) THEN
            updated=updated+1
            topoind_r(i,j) = 1e-10  ! flow instantaneously
          ENDIF
        ENDDO
      ENDDO
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, updated, 1, MPI_INT_ORCH, MPI_SUM, MPI_COMM_ORCH, ierr)
      IF (updated/=0 .AND. is_mpi_root)  WRITE(numout,*) "WARNING : ",updated," routing cells have a zero value for topographic index, which have been replaced by 1e-10"
      
      CALL xios_send_field("topoind_r_final",topoind_r)
      
      CALL xios_context_finalize
      CALL xios_orchidee_change_context("orchidee")
    
  CONTAINS

    SUBROUTINE next_trip(trip, ni, nj, i, j, nexti, nextj, ok)
    IMPLICIT NONE
      REAL, INTENT(IN)     :: trip(0:ni+1,0:nj+1)
      INTEGER, INTENT(IN)  :: ni
      INTEGER, INTENT(IN)  :: nj
      INTEGER, INTENT(IN)  :: i
      INTEGER, INTENT(IN)  :: j
      INTEGER, INTENT(OUT) :: nexti
      INTEGER, INTENT(OUT) :: nextj
      LOGICAL, INTENT(OUT) :: ok

      ok=.TRUE. 
      
      SELECT CASE(NINT(trip(i,j)))
        CASE (1)
           nextj=j-1 ; nexti=i    ! north
        CASE (2)
           nextj=j-1 ; nexti=i+1  ! north-east
        CASE (3)
           nextj=j ;   nexti=i+1  ! east
        CASE (4)
           nextj=j+1 ;  nexti=i+1 ! south-east
        CASE (5)
           nextj=j+1 ;  nexti=i   ! south
        CASE(6)
           nextj=j+1 ;  nexti=i-1 ! south-west
        CASE (7)
           nextj=j ;   nexti=i-1  ! west
        CASE (8)
           nextj=j-1 ;  nexti=i-1 ! north-west
        CASE DEFAULT
          ok=.FALSE.
       END SELECT
       
     END SUBROUTINE next_trip
     
  END SUBROUTINE routing_flow_correct_riverflow


  
  SUBROUTINE routing_flow_init_mean(kjit, rest_id)
  USE ioipsl   
  IMPLICIT NONE
  INTEGER(i_std), INTENT(in)    :: kjit
  INTEGER(i_std), INTENT(in)    :: rest_id        !! Restart file identifier (unitless)
 
  INTEGER :: ier    
  CHARACTER(LEN=80)   :: var_name       !! To store variables names for I/O (unitless)
    

    ! Get from the restart the fluxes we accumulated.
    !
  
    ALLOCATE (runoff_mean(nbpt), stat=ier)
    runoff_mean(:) = 0
    IF (ier /= 0) CALL ipslerr_p(3,'routing_flow_init_mean','Pb in allocate for runoff_mean','','')
    var_name = 'runoff_route'
    CALL ioconf_setatt_p('UNITS', 'Kg')
    CALL ioconf_setatt_p('LONG_NAME','Accumulated runoff for routing')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., runoff_mean, "gather", nbp_glo, index_g)
    CALL setvar_p (runoff_mean, val_exp, 'NO_KEYWORD', zero)

    ALLOCATE(drainage_mean(nbpt), stat=ier)
    drainage_mean(:) = 0
    IF (ier /= 0) CALL ipslerr_p(3,'routing_flow_init_mean','Pb in allocate for drainage_mean','','')
    var_name = 'drainage_route'
    CALL ioconf_setatt_p('UNITS', 'Kg')
    CALL ioconf_setatt_p('LONG_NAME','Accumulated drainage for routing')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., drainage_mean, "gather", nbp_glo, index_g)
    CALL setvar_p (drainage_mean, val_exp, 'NO_KEYWORD', zero)
    
  END SUBROUTINE routing_flow_init_mean

  SUBROUTINE routing_flow_make_mean(runoff, drainage)
  IMPLICIT NONE
    REAL(r_std),INTENT(IN) :: runoff(:)
    REAL(r_std),INTENT(IN) :: drainage(:)
    
    runoff_mean(:) = runoff_mean(:) + runoff
    drainage_mean(:) =  drainage_mean(:) + drainage
    
  END SUBROUTINE routing_flow_make_mean
  

  SUBROUTINE routing_flow_reset_mean
  IMPLICIT NONE
    runoff_mean(:) = 0
    drainage_mean(:) =  0
    
  END SUBROUTINE routing_flow_reset_mean


  SUBROUTINE routing_flow_finalize_mean(kjit, rest_id)
  USE ioipsl
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: kjit
  INTEGER,INTENT(IN) :: rest_id
 
    CALL restput_p (rest_id, 'runoff_route', nbp_glo, 1, 1, kjit, runoff_mean, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'drainage_route', nbp_glo, 1, 1, kjit, drainage_mean, 'scatter',  nbp_glo, index_g)
     
    DEALLOCATE(runoff_mean)
    DEALLOCATE(drainage_mean)
    
  END SUBROUTINE routing_flow_finalize_mean


  SUBROUTINE routing_flow_init_local(contfrac, nbpt_r, dt_routing)
    USE xios
    USE grid, ONLY : area
    IMPLICIT NONE
    INCLUDE "mpif.h"

    !! 0 Variable and parameter description
    !! 0.1 Input variables
    REAL(r_std),INTENT(IN)                 :: contfrac(nbpt)   !! fraction of land
    INTEGER,INTENT(OUT)                    :: nbpt_r           !! nb points routing grid
    REAL(r_std),INTENT(IN)                     :: dt_routing       !! routing timestep (s)

    !! 0.2 Local variables
    INTEGER :: ni                                         !! longitude dimension of local routing grid
    INTEGER :: nj                                         !! latitude dimension of local routing grid
    REAL(r_std)             :: contfrac_mpi(nbp_mpi)
    REAL(r_std),ALLOCATABLE :: trip_rp1(:,:)             !! direction of flow (1-8) or river flow (99) or coastal flow (98) or lake inflow (97) - local routing grid + halo (0:ni+1,0:nj+1)
    REAL(r_std),ALLOCATABLE :: trip_extended_r(:,:)       !! direction of flow (1-8) or river flow (99) or coastal flow (98) or lake inflow (97) - local routing grid (ni,nj)
    !! routing is artificially computed on sea and endorheic basins
    REAL(r_std),ALLOCATABLE :: diag_r(:)    !! fraction of routing cell intersected by coastal cells of native grid
    REAL(r_std),ALLOCATABLE :: frac_lake_r(:)          !! fraction of routing cell intersected by cells of native grid
    REAL(r_std),ALLOCATABLE :: frac_coast_r(:)    
  
    REAL(r_std)             :: diag(nbp_mpi)    !! fraction of routing cell intersected by coastal cells of native grid
    
    LOGICAL :: coastline(nbpt)
    LOGICAL :: coastline_mpi(nbp_mpi)
    REAL(r_std) :: area_mpi(nbp_mpi)

    INTEGER :: ij, ij_r, ij_rp1, i,j,jp1,jm1,ip1,im1,jr,ir
    INTEGER :: basins_count_mpi
    INTEGER :: nb_coast_points
    INTEGER :: ierr
    LOGICAL :: file_exists
    REAL(r_std) :: epsilon = 1e-5
    CHARACTER(LEN=255) :: routing_file_type


    printlev_loc = get_printlev("routing_interp_topo_flow")
!_ ================================================================================================================================
    !
    !> A value for property of each reservoir (in day/m) is given to compute a time constant (in day)
    !> for each reservoir (product of tcst and topo_resid).
    !> The value of tcst has been calibrated for the three reservoirs over the Senegal river basin only,
    !> during the 1 degree NCEP Corrected by Cru (NCC) resolution simulations (Ngo-Duc et al., 2005, Ngo-Duc et al., 2006) and
    !> generalized for all the basins of the world. The "slow reservoir" and the "fast reservoir"
    !> have the highest value in order to simulate the groundwater. 
    !> The "stream reservoir", which represents all the water of the stream, has the lowest value.
    !> Those figures are the same for all the basins of the world.

    CALL getin_p('SLOW_TCST', slow_tcst)
    !
    !Config Key   = FAST_TCST
    !Config Desc  = Time constant for the fast reservoir 
    !Config If    = RIVER_ROUTING 
    !Config Def   = 3.0
    !Config Help  = This parameters allows the user to fix the 
    !Config         time constant (in days) of the fast reservoir
    !Config         in order to get better river flows for 
    !Config         particular regions.
    !Config Units = [days]

    CALL getin_p('FAST_TCST', fast_tcst)
    !
    !Config Key   = STREAM_TCST
    !Config Desc  = Time constant for the stream reservoir 
    !Config If    = RIVER_ROUTING
    !Config Def   = 0.24
    !Config Help  = This parameters allows the user to fix the 
    !Config         time constant (in days) of the stream reservoir
    !Config         in order to get better river flows for 
    !Config         particular regions.
    !Config Units = [days]

    CALL getin_p('STREAM_TCST', stream_tcst)
    !
    !Config Key   = FLOOD_TCST
    !Config Desc  = Time constant for the flood reservoir 
    !Config If    = RIVER_ROUTING
    !Config Def   = 4.0
    !Config Help  = This parameters allows the user to fix the 
    !Config         time constant (in days) of the flood reservoir
    !Config         in order to get better river flows for 
    !Config         particular regions.
    !Config Units = [days]

    CALL getin_p('SLOW_TCST', slow_tcst)

    routing_file_type = "standard"
    CALL getin_p("routing_file_type", routing_file_type) 

    IF (TRIM(routing_file_type)=="standard") THEN
      topoind_factor=1000
    ELSE IF (TRIM(routing_file_type)=="merit") THEN
      topoind_factor=1
    ELSE
      CALL ipslerr_p(3,'routing_flow_init_local', &
        'Error in call routing_flow_init_local','getin routing_file_type : bad value, must be <standard> or <merit>','')
    ENDIF

    split_routing=1
    CALL getin_p("SPLIT_ROUTING",split_routing)
    
    CALL compute_coastline(contfrac, coastline)
    
    CALL gather_omp(contfrac, contfrac_mpi)
    CALL gather_omp(area, area_mpi)
    CALL gather_omp(coastline, coastline_mpi)
    
    ALLOCATE(is_coastline(nbp_mpi))
    is_coastline=coastline_mpi

    IF (is_omp_root) THEN
       WHERE(contfrac_mpi <= epsilon) contfrac_mpi=0
       WHERE(contfrac_mpi >= 1-epsilon) contfrac_mpi=1

       CALL xios_get_domain_attr("routing_domain", ni=ni, nj=nj)    ! get routing domain dimension

       nbpt_r= ni*nj                                                
       nbpt_rp1= (ni+2)*(nj+2)

       ! Allocate module variable     
       ALLOCATE(fast_reservoir_r(ni*nj))       
       ALLOCATE(slow_reservoir_r(ni*nj))       
       ALLOCATE(stream_reservoir_r(ni*nj))       

       ALLOCATE(is_lakeinflow_r(ni*nj))       
       ALLOCATE(is_coastalflow_r(ni*nj))       
       ALLOCATE(is_riverflow_r(ni*nj)) 
       ALLOCATE(is_streamflow_r(ni*nj)) 

       ALLOCATE(topoind_r(ni*nj))       
       ALLOCATE(routing_mask_r(ni*nj)) 
       ALLOCATE(basins_extended_r(nbpt_r))
       ALLOCATE(route_flow_rp1((ni+2)*(nj+2)))

       ALLOCATE(routing_weight(nbp_mpi))
       ALLOCATE(routing_weight_in(nbp_mpi))
       ALLOCATE(unrouted_weight(nbp_mpi))
      
       ALLOCATE(diag_r(nbpt_r)) ! for diags
       ALLOCATE(frac_coast_r(nbpt_r))
       ALLOCATE(frac_lake_r(nbpt_r))
       ALLOCATE(weight_coast_to_coast_r(nbpt_r))
       ALLOCATE(weight_coast_to_lake_r(nbpt_r))
       ALLOCATE(weight_lake_to_coast_r(nbpt_r))
       ALLOCATE(weight_lake_to_lake_r(nbpt_r))

       ALLOCATE(trip_extended_r(ni,nj))

      ! correction on coastline in case of LAM, the cells near the frontier are considered as coastline 
      ! => interpolate model grid to routing grid and reinterpolate to model grid. Fractionnal cells will be considered as coastline too
    !  diag=1
    !  CALL xios_send_field("one", diag)
    !  CALL xios_recv_field("tmp1_coastline", diag_r)
    !  CALL xios_send_field("tmp2_coastline", diag_r)
    !  CALL xios_recv_field("lam_coastline", diag)
    !  WHERE(diag<1-epsilon) is_coastline=.TRUE.

       diag=0
       WHERE(is_coastline) diag=1
       CALL xios_send_field("is_coastline", diag)

       is_lakeinflow_r(:) = .FALSE.
       is_coastalflow_r(:) = .FALSE.
       is_riverflow_r(:) = .FALSE.
   
       
       INQUIRE(FILE="routing_start.nc", EXIST=file_exists)

       IF (file_exists) THEN  
          CALL xios_recv_field("fast_reservoir_start",fast_reservoir_r)
          CALL xios_recv_field("slow_reservoir_start",slow_reservoir_r)
          CALL xios_recv_field("stream_reservoir_start",stream_reservoir_r)
          CALL xios_recv_field("integration_time_start",integration_time)
       ELSE
          fast_reservoir_r(:)=0
          slow_reservoir_r(:)=0
          stream_reservoir_r(:)=0
          integration_time = dt_routing   ! to be saved in restart
       ENDIF
       
       ALLOCATE(trip_rp1(0:ni+1,0:nj+1))
       
       trip_rp1(:,:)=1e10
       CALL xios_recv_field("trip_r",trip_rp1(1:ni,1:nj))          ! recv trip array with halo of 1
       CALL xios_recv_field("trip_extended_r",trip_extended_r)     ! recv extended trip array from file
       CALL xios_recv_field("topoind_r",topoind_r)                 ! recv topo index array from file
       CALL xios_recv_field("basins_extended_r",basins_extended_r) ! recv basins index from file
 
       CALL routing_flow_correct_riverflow(ni, nj, contfrac_mpi, is_coastline, trip_rp1(1:ni,1:nj), trip_extended_r, topoind_r) 
   
       CALL xios_send_field("trip_update_r",trip_rp1(1:ni,1:nj))        ! send to xios trip array to update for halo
       CALL xios_recv_field("trip_rp1",trip_rp1)                        ! recv trip array with halo of 1
       


       !! Compute the routing
       !! Loop on all point point of the local routing grid + halo
  
       route_flow_rp1(:)=-1

       DO j=0,nj+1                                                 
          jp1=j+1
          jm1=j-1
          DO i=0,ni+1
             ij_rp1=i+(ni+2)*j+1
             ij_r=i+ni*(j-1)
             ip1=i+1
             im1=i-1

             IF (trip_rp1(i,j) < 100) THEN                    
                
                IF (i>=1 .AND. i<=ni .AND. j>=1 .AND. j<=nj)   routing_mask_r(ij_r)=.TRUE.
                
                ir=-1 ; jr=-1  !  -> -1 for 97,98,99
                SELECT CASE (NINT(trip_rp1(i,j)))                     ! get the trip value for each points
                CASE (1)
                   jr=jm1 ; ir=i    ! north
                CASE (2)
                   jr=jm1 ; ir=ip1  ! north-east
                CASE (3)
                   jr=j ;   ir=ip1  ! east
                CASE (4)
                   jr=jp1 ;  ir=ip1 ! south-east
                CASE (5)
                   jr=jp1 ;  ir=i   ! south
                CASE(6)
                   jr=jp1 ;  ir=im1 ! south-west
                CASE (7)
                   jr=j ;   ir=im1  ! west
                CASE (8)
                   jr=jm1 ;  ir=im1 ! north-west
                CASE (97)
                   IF ( i>0 .AND. i<ni+1 .AND. j>0 .AND. j<nj+1)  THEN         ! if inside my local domain
                      is_lakeinflow_r(ij_r)=.TRUE.                             ! I am a lakeinflow point and route to myself
                      jr=j ;   ir=i                                
                   ENDIF
                CASE (98)
                   IF ( i>0 .AND. i<ni+1 .AND. j>0 .AND. j<nj+1) THEN         ! if inside my local domain
                      is_coastalflow_r(ij_r)=.TRUE.                           ! I am a coastal flow point and route to myself           
                      jr=j ;   ir=i
                   ENDIF
                CASE (99)
                   IF ( i>0 .AND. i<ni+1 .AND. j>0 .AND. j<nj+1) THEN         ! if inside my local domain
                      is_riverflow_r(ij_r)=.TRUE.                             ! I am a riverflow point and route to myself
                      jr=j ;   ir=i
                   ENDIF
                END SELECT

                IF (ir<0 .OR. ir>ni+1 .OR. jr<0 .OR. jr>nj+1) THEN
                   route_flow_rp1(ij_rp1)=-1                                     ! if route outside my local domain+halo, no routing (will be done by other process)
                ELSE 
                  IF (ir<1 .OR. ir>ni .OR. jr<1 .OR. jr>nj) THEN  
                      route_flow_rp1(ij_rp1)=-1                                  ! if route outside my local domain, no routing (will be done by other process)
                   ELSE
                      route_flow_rp1(ij_rp1)=ir+ni*(jr-1)                        ! define the cell where to flow
                      IF (trip_rp1(ir,jr)>99) STOP 'Pb point not routed to outflow'
                   ENDIF
                ENDIF
             ELSE
               IF (i>=1 .AND. i<=ni .AND. j>=1 .AND. j<=nj)  routing_mask_r(ij_r)=.FALSE.             
             ENDIF
          ENDDO
       ENDDO
       is_streamflow_r(:) = .NOT. (is_lakeinflow_r(:) .OR. is_coastalflow_r(:) .OR. is_riverflow_r(:)) .AND. routing_mask_r(:)
       
       diag_r(:)=0
       WHERE(routing_mask_r) diag_r=1
       CALL xios_send_field("routing_mask_r", diag_r)
       CALL xios_recv_field("frac_routing", diag)
       WHERE (diag < epsilon) diag=0
       routing_weight_in (:) = 0
       unrouted_weight (:) = 0
       WHERE (diag > 0) routing_weight_in = contfrac_mpi*area_mpi/diag
       WHERE (diag == 0 ) unrouted_weight = contfrac_mpi*area_mpi
       

      ! compute the number of coast cells to distribute lakeinflow onto
       diag=0
       WHERE (is_coastline) diag=1
       nb_coast_points=SUM(diag)
       CALL reduce_sum_mpi(nb_coast_points, total_coast_points)
       CALL bcast_mpi(total_coast_points)


       DO ij=1,nbp_mpi
          routing_weight(ij)=contfrac_mpi(ij)*area_mpi(ij)  
       ENDDO

       ! looking for basins
       basins_count_mpi=0
       DO ij=1,nbpt_r
          IF (basins_extended_r(ij) > basins_count_mpi) basins_count_mpi=basins_extended_r(ij)
       ENDDO
       CALL MPI_ALLREDUCE(basins_count_mpi,basins_count,1,MPI_INT_ORCH, MPI_MAX,MPI_COMM_ORCH,ierr)

       diag(:) = 0
       WHERE (is_coastline) diag=1
       CALL xios_send_field("mask_coastal",diag)
       CALL xios_recv_field("frac_coastal_r",diag_r)
       WHERE (diag_r > 100) ! missing_value
         diag_r=0
       ELSE WHERE (diag_r < epsilon) 
         diag_r=0
       ELSEWHERE (diag_r > 1-epsilon) 
         diag_r=1
       END WHERE
       frac_coast_r=diag_r ;

       diag(:) = 0
       WHERE (.NOT. is_coastline) diag=1
       CALL xios_send_field("mask_lake",diag)
       CALL xios_recv_field("frac_lake_r",diag_r)
       WHERE (diag_r > 100) ! missing_value
         diag_r=0.
       ELSEWHERE (diag_r < epsilon) 
         diag_r=0.
       ELSEWHERE (diag_r > 1-epsilon) 
         diag_r=1.
       END WHERE
       frac_lake_r = diag_r

       weight_coast_to_coast_r(:)=0
       weight_coast_to_lake_r(:)=0
       weight_lake_to_coast_r(:)=0
       weight_lake_to_lake_r(:)=0

       DO ij=1,nbpt_r
         IF (is_riverflow_r(ij) .OR. is_coastalflow_r(ij) ) THEN
           IF (frac_coast_r(ij)==0 .AND. frac_lake_r(ij)==0) THEN
             WRITE (numout,*) "riverflow or costalflow point on routing grid can not be interpolated on orchidee grid, this might not happen"
             STOP
           ELSE IF (frac_coast_r(ij)==0 .AND. frac_lake_r(ij)>0) THEN
             weight_coast_to_lake_r(ij) = 1./frac_lake_r(ij)
             weight_coast_to_coast_r(ij) = 0
           ELSE IF (frac_coast_r(ij)>0) THEN
             weight_coast_to_coast_r(ij) = 1./(frac_coast_r(ij))
             weight_coast_to_lake_r(ij) = 0
           ELSE
           ENDIF
         ELSE IF (is_lakeinflow_r(ij)) THEN
           IF ( frac_coast_r(ij)==0 .AND. frac_lake_r(ij)==0) THEN
             WRITE (numout,*) "lakeinflow on routing grid can not be interpolated on orchidee grid, this might not happen"
             STOP
           ELSE IF (frac_lake_r(ij)==0 .AND. frac_coast_r(ij)>0) THEN
             weight_lake_to_coast_r(ij) = 1./frac_coast_r(ij)
             weight_lake_to_lake_r(ij) = 0
           ELSE IF (frac_lake_r(ij)>0 ) THEN
             weight_lake_to_lake_r(ij) = 1./(frac_lake_r(ij))
             weight_lake_to_coast_r(ij) = 0
           ENDIF
         ENDIF
       ENDDO
       
       CALL compute_basins_area(nbpt_r,nbpt_rp1)

    ELSE

       nbpt_r=0

    ENDIF !is_omp_root

    CALL bcast_omp(integration_time)
    
  END SUBROUTINE routing_flow_init_local     
  
  
  !!  =============================================================================================================================
  !! SUBROUTINE:         compute_basins_area
  !!
  !>\BRIEF	         Compute the basin upstream area of each cell of the routing grid
  !!
  !! DESCRIPTION:        Works as a cascade of areas :
  !!                     Each routing cell iteratively sends area information to the adjacent downstream cells,
  !!                     and receives area information from the adjacent upstream cells to update its basins_area_r.
  !!                     When no more area information is sent (sum_area=0), the variable basins_area_r has the right values.
  !!                     Must be called from an omp_root thread.
  !!
  !! RECENT CHANGE(S)
  !!
  !! REFERENCE(S)
  !! 
  !! FLOWCHART   
  !! \n
  !_ ==============================================================================================================================
  SUBROUTINE compute_basins_area(nbpt_r,nbpt_rp1 )
  USE xios
  USE routing_interp_topo_para
  IMPLICIT NONE
    INCLUDE 'mpif.h'
    INTEGER, INTENT(IN) :: nbpt_r
    INTEGER, INTENT(IN) :: nbpt_rp1
    REAL(r_std),PARAMETER :: Pi=atan(1.)*4
    REAL(r_std),PARAMETER :: earth_radius=6.371009E6 ! in meter

    INTEGER ::  ni, nj
    REAL(r_std),ALLOCATABLE :: lon(:), lat(:)
    REAL(r_std)  :: area_rp1(nbpt_rp1)
    REAL(r_std)  :: send_area_rp1(nbpt_rp1)
    REAL(r_std)  :: recv_area_r(nbpt_rp1)
    REAL :: delta_lon, delta_lat
    INTEGER :: i,j,ij,it,ig,ierr
    REAL(r_std) :: sum_area
   
    ALLOCATE(basins_area_r(nbpt_r))

    CALL xios_get_domain_attr("routing_domain", ni=ni, nj=nj)
    ALLOCATE(lon(ni),lat(nj))
    CALL xios_get_domain_attr("routing_domain",lonvalue_1d=lon, latvalue_1d=lat)
    
   !! Compute the routing cell areas (regular lon-lat cells)
   !! and initialize area_rp1 and basins_area_r in the local domain only
    delta_lon=ABS(lon(2)-lon(1))*Pi/180.
    delta_lat=ABS(lat(2)-lat(1))*Pi/180.

    DO j=1,nj 
      DO i=1,ni
        ig=(j+1-1)*(ni+2)+i+1 
        area_rp1(ig) = delta_lon*delta_lat*cos(lat(j)*Pi/180.)*earth_radius*earth_radius  ! pretty good approximation
        ij = (j-1)*ni+i
        basins_area_r(ij) = delta_lon*delta_lat*cos(lat(j)*Pi/180.)*earth_radius*earth_radius  ! pretty good approximation
      ENDDO
    ENDDO
    
    !! Update area_rp1 in the halo (from adjacent local domains)
    CALL update_halo(area_rp1)

    !! The first area "sent" by the routing cells is their own area
    send_area_rp1(:) = area_rp1(:)
    !! ... except the riverflow/coastalflow/lakeinflow cells, which are routed to themselves
    DO j=1,nj
      DO i=1,ni
        ij=(j-1)*ni+i
        ig=(j+1-1)*(ni+2)+i+1
        IF (is_riverflow_r(ij) .OR. is_coastalflow_r(ij) .OR. is_lakeinflow_r(ij)) send_area_rp1(ig)=0
      ENDDO
    ENDDO 
    
    sum_area=1.
    DO WHILE (sum_area/=0) ! while there is something to send downstream

      !! (Re)set recv_area_r at the beginning of each iteration
      recv_area_r(:) = 0
      
      !! Update basins_area_r and recv_area_r
      DO ig=1,nbpt_rp1  ! loop on local domain+halo
        IF ( route_flow_rp1(ig) > 0 ) THEN ! if the cell (ig) flows into a cell of this local domain
          !! Add the area sent by the flowing cell to the basins_area_r of the downstream cell
          basins_area_r(route_flow_rp1(ig)) = basins_area_r(route_flow_rp1(ig)) + send_area_rp1(ig) 
          recv_area_r(route_flow_rp1(ig)) = recv_area_r(route_flow_rp1(ig)) + send_area_rp1(ig)
        ENDIF
      ENDDO  ! end loop on local domain+halo
      
      !! Update the sent area (send_area_rp1) in each cell of the local domain
      DO j=1,nj
        DO i=1,ni
          ij=(j-1)*ni+i
          ig=(j+1-1)*(ni+2)+i+1 
          send_area_rp1(ig) = recv_area_r(ij)
          IF (is_riverflow_r(ij) .OR. is_coastalflow_r(ij) .OR. is_lakeinflow_r(ij)) send_area_rp1(ig)=0 ! riverflow/coastalflow/lakeinflow cells dont send area
        ENDDO
      ENDDO
      !! Update send_area_rp1 in the halo (from adjacent local domains)
      CALL update_halo(send_area_rp1)
      
      !! Compute the total sent area in the local domain+halo and store it in sum_area
      sum_area = sum(send_area_rp1)
      !! Sum these total sent areas from all local domains+halos and overwrite sum_area with it
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, sum_area, 1, MPI_REAL_ORCH, MPI_SUM, MPI_COMM_ORCH, ierr)

    ENDDO

  END SUBROUTINE compute_basins_area

  
  !!  =============================================================================================================================
  !! SUBROUTINE:         initialize_stations
  !!
  !>\BRIEF	         Snap the user's stations on the routing grid and initialize the output file that will contain the hydrographs (stations.nc)
  !!
  !! DESCRIPTION:        Read the list of stations in run.def with their information (coordinates and upstream area)
  !!                     Try to place them at best on the routing grid using a research perimeter around the coordinates, with 2 possible methods :
  !!                     1. if the station's area is provided, find the routing cell that has the closest upstream area (computed by compute_basins_area)
  !!                        within the research perimeter
  !!                     2. otherwise (or if the routing cell upstream areas are too different from the provided station area),
  !!                        find the biggest river (biggest upstream area) within the research perimeter
  !!
  !! RECENT CHANGE(S)
  !!
  !! REFERENCE(S)
  !! 
  !! FLOWCHART   
  !! \n
  !_ ==============================================================================================================================
  SUBROUTINE initialize_stations(dt_routing)
  USE xios
  USE mod_orchidee_para, ONLY : is_omp_root
  IMPLICIT NONE
    INCLUDE 'mpif.h'
    REAL(r_std), INTENT (in)                     :: dt_routing                !! Routing time step (s)

    REAL(r_std),PARAMETER :: Pi=atan(1.)*4
    REAL(r_std),PARAMETER :: earth_radius=6.371009E6 ! in meter
    INTEGER ::  ni, nj
    REAL(r_std),ALLOCATABLE :: lon(:), lat(:) ! routing grid axes (degrees)
  
    REAL,ALLOCATABLE :: station_lonlat(:,:)  ! station coordinates (lon,lat) (degrees)
    REAL,ALLOCATABLE :: station_prec(:)  ! research perimeter around the station coordinates (m)
    REAL,ALLOCATABLE :: station_area(:)  ! basin area upstream of the stations (km2)
    INTEGER :: station_index_prec, station_index_prec_i, station_index_prec_j
    INTEGER :: station_index_area, station_index_area_i, station_index_area_j
    INTEGER :: max_station_index_area  ! to find the global maximum station_index_area
    LOGICAL :: has_lonlat
    LOGICAL :: has_area
    INTEGER ::i,j,r,k
  
    TYPE(xios_duration) :: ts 
    TYPE(xios_scalar) :: scalar_hdl
    TYPE(xios_scalargroup) :: scalar_def_hdl
    TYPE(xios_field) :: field_hdl
    TYPE(xios_fieldgroup) :: field_def_hdl, fieldgroup_hdl
    TYPE(xios_file) :: file_hdl
    TYPE(xios_filegroup) :: file_def_hdl
    TYPE(xios_date) :: start_date, time_origin
    CHARACTER(LEN=20)  :: calendar_type

    CHARACTER(LEN=256) :: str_station_ind
    REAL :: lon_a,lat_a,lon_b,lat_b,dist,max_area, max_max_area, min_dist, min_area, min_min_area
    INTEGER :: min_dist_index, min_dist_index_i, min_dist_index_j
    INTEGER :: ierr

    IF (is_omp_root) THEN

      !! 1. Read and snap the stations given by the user
      
       CALL xios_get_domain_attr("routing_domain", ni=ni, nj=nj)
       ALLOCATE(lon(ni),lat(nj))
       CALL xios_get_domain_attr("routing_domain",lonvalue_1d=lon, latvalue_1d=lat)

       nb_station=0
       CALL getin("nb_station",nb_station)
       
       IF (nb_station/=0) THEN
         ALLOCATE(station(nb_station)) 
         ALLOCATE(station_lonlat(nb_station,2)) 
         ALLOCATE(station_index(nb_station)) 
         ALLOCATE(station_prec(nb_station)) 
         ALLOCATE(station_area(nb_station))

         DO k=1,nb_station
           WRITE(str_station_ind,*) k
           str_station_ind=ADJUSTL(str_station_ind)
           CALL getin("station"//TRIM(str_station_ind)//"_id",station(k))
           station_lonlat(k,:) = 1.1e10
           has_lonlat=.TRUE.
           CALL getin("station"//TRIM(str_station_ind)//"_coor",station_lonlat(k,:))
           IF (station_lonlat(k,1)>1e10 .OR. station_lonlat(k,2)>1e10) has_lonlat=.FALSE.
           station_prec(k)=50000 ! meters ! NB: if routing_file_type=standard (0.5deg resolution), 50km is ~1 routing cell ! one may want to look wider
           CALL getin("station"//TRIM(str_station_ind)//"_prec",station_prec(k))
           station_area(k) = -1
           has_area=.TRUE.
           CALL getin("station"//TRIM(str_station_ind)//"_area",station_area(k))
           IF (station_area(k)==-1) has_area=.FALSE.

           lon_a = station_lonlat(k,1)*Pi/180.
           lat_a = station_lonlat(k,2)*Pi/180.

           max_area=0
           min_dist=HUGE(min_dist)
           min_area=HUGE(min_area)

           DO j=1,nj 
             DO i=1,ni
               r=((j-1)*ni)+i
               IF (routing_mask_r(r)) THEN
                 lon_b=lon(i)*Pi/180.
                 lat_b=lat(j)*Pi/180.
                 dist=earth_radius*acos(sin(lat_a)*sin(lat_b)+cos(lat_a)*cos(lat_b)*cos(lon_b-lon_a)) ! Spherical Law of Cosines
                 IF (dist<station_prec(k)) THEN ! if the routing cell is within the research perimeter
                   ! Find the closest routing cell to the station's coordinates (not used)
                   IF (dist<min_dist) THEN
                     min_dist=dist
                     min_dist_index = r
                     min_dist_index_i = i
                     min_dist_index_j = j
                   ENDIF

                   ! Find the routing cell with the biggest upstream area in the research perimeter (biggest river, most downstream routing cell)
                   IF (basins_area_r(r) > max_area) THEN
                     max_area = basins_area_r(r)
                     station_index_prec = r
                     station_index_prec_i = i
                     station_index_prec_j = j
                   ENDIF

                   ! Find the routing cell with the closest upstream area value to the station's area in the research perimeter
                   IF (has_area) THEN
                     IF (ABS(basins_area_r(r)-station_area(k)*1e6) < min_area) THEN
                       min_area = ABS(basins_area_r(r)-station_area(k)*1e6)
                       station_index_area = r
                       station_index_area_i = i
                       station_index_area_j = j
                     ENDIF
                   ENDIF

                 ENDIF
               ENDIF
             ENDDO
           ENDDO      

           ! Compare with other domains (MPI)
           CALL MPI_ALLREDUCE(max_area, max_max_area, 1, MPI_REAL_ORCH, MPI_MAX, MPI_COMM_ORCH, ierr)
           IF (max_area /= max_max_area) station_index_prec=-1
           IF (max_area == 0) station_index_prec=-1

           IF (has_area) THEN
             CALL MPI_ALLREDUCE(min_area, min_min_area, 1, MPI_REAL_ORCH, MPI_MIN, MPI_COMM_ORCH, ierr)
             IF (min_area /= min_min_area) station_index_area=-1
             IF (min_area == HUGE(min_area)) station_index_area=-1
           ELSE
             station_index_area=-1
           ENDIF

           ! Do the snapping
           station_index(k)=-1 ! initialization
           IF (station_index_prec/=-1 .OR. station_index_area/=-1) THEN
             IF (printlev>=2) THEN
               WRITE(numout,*)"Station ",TRIM(station(k))
               WRITE(numout,*)"  - coordinates lon-lat: ",station_lonlat(k,:)
               IF (has_area) WRITE(numout,*)"  - upstream basin area (km2): ", station_area(k)
               WRITE(numout,*)"=> Based on a research perimeter of ", station_prec(k)/1000, "km:"
             ENDIF
           ENDIF

           IF (station_index_area/=-1) THEN
             ! Do the "closest upstream area" snapping
             IF (printlev>=2) THEN
               WRITE(numout,*)"- the corresponding routing cell with the closest upstream area (", &
                      basins_area_r(station_index_area)/1000/1000, "km2) is at coordinates lon-lat", &
                      lon(station_index_area_i),lat(station_index_area_j)
               WRITE(numout,*)"  The station is snapped on this routing cell."
             ENDIF
             station_index(k) = station_index_area
           ENDIF
           
           ! Check on all domains (MPI) if we did a "closest upstream area" snapping (max_station_index_area/=-1)
           CALL MPI_ALLREDUCE(station_index_area, max_station_index_area, 1, MPI_INT_ORCH, MPI_MAX, MPI_COMM_ORCH, ierr)
             
           IF (station_index_prec/=-1) THEN
             IF (printlev>=2) THEN
               WRITE(numout,*)"- the corresponding routing cell with the biggest upstream area (", &
                      basins_area_r(station_index_prec)/1000/1000, "km2) is at coordinates lon-lat", &
                      lon(station_index_prec_i),lat(station_index_prec_j)
             ENDIF
             IF (max_station_index_area==-1) THEN ! if there was no snapping by "closest upstream area" done anywhere
               ! Do the "biggest river not too far" snapping
               IF (printlev>=2) WRITE(numout,*)"  The station is snapped on this routing cell."
               station_index(k) = station_index_prec
             ENDIF
           ENDIF

        ENDDO ! k=1,nb_station

        !! 2. Initialize the output file that will contain the hydrographs of the stations

         CALL xios_get_calendar_type(calendar_type)
         CALL xios_get_start_date(start_date)
         CALL xios_get_time_origin(time_origin)

         CALL xios_context_initialize("orchidee_routing_out",MPI_COMM_ORCH)
         CALL xios_orchidee_change_context("orchidee_routing_out")
         ts%second=dt_routing
         calendar_type="gregorian" 
         CALL xios_define_calendar(type=calendar_type,start_date=start_date,time_origin=time_origin, timestep=ts )

         CALL xios_get_handle("scalar_definition",scalar_def_hdl)

         CALL xios_get_handle("field_definition",field_def_hdl)
         CALL xios_add_child(field_def_hdl, fieldgroup_hdl, "stations")
         DO k=1,nb_station
           CALL xios_add_child(scalar_def_hdl,scalar_hdl,TRIM(station(k))//"_scalar")
           CALL xios_set_attr(scalar_hdl,label=TRIM(station(k)))
           CALL xios_add_child(fieldgroup_hdl,field_hdl, TRIM(station(k)))
           CALL xios_set_attr(field_hdl,scalar_ref=TRIM(station(k))//"_scalar")
         ENDDO

         CALL xios_get_handle("file_definition",file_def_hdl)
         CALL xios_add_child(file_def_hdl, file_hdl, "stations")
         CALL xios_set_attr(file_hdl, type="one_file", output_freq=ts, sync_freq=ts, enabled=.TRUE.)
         CALL xios_add_child(file_hdl, fieldgroup_hdl)
         CALL xios_set_attr(fieldgroup_hdl, group_ref="stations", operation="average")
         CALL xios_close_context_definition()
         CALL xios_orchidee_change_context("orchidee")
      END IF ! nb_station/=0
    END IF ! is_omp_root
    
  END SUBROUTINE initialize_stations
  

  SUBROUTINE routing_flow_main(dt_routing, contfrac) !
    
    USE xios
    USE grid, ONLY : area 
    USE routing_interp_topo_para
    IMPLICIT NONE
    INCLUDE "mpif.h"
    
    !! 0 Variable and parameter description
    !! 0.1 Input variables
    REAL(r_std), INTENT (in)                     :: dt_routing                !! Routing time step (s)
    REAL(r_std), INTENT (in)                     :: contfrac(nbpt)                !! Routing time step (s)

    !! 0.4 Local variables
    REAL(r_std)                                  :: runoff(nbp_mpi)              !! Grid-point runoff (kg/dt)
    REAL(r_std)                                  :: runoff_in(nbp_mpi)              !! Grid-point runoff (kg/dt)
    REAL(r_std)                                  :: drainage(nbp_mpi)            !! Grid-point drainage (kg/dt)
    REAL(r_std)                                  :: drainage_in(nbp_mpi)            !! Grid-point drainage (kg/dt)
    REAL(r_std)                                  :: riverflow(nbp_mpi)             
    REAL(r_std)                                  :: coastalflow(nbp_mpi)           
    REAL(r_std)                                  :: lakeinflow(nbp_mpi)            
    REAL(r_std)                                  :: area_mpi(nbp_mpi) ! cell area         
    REAL(r_std)                                  :: contfrac_mpi(nbp_mpi) ! cell area         
    REAL(r_std)                                  :: flow_coast(nbp_mpi)          
    REAL(r_std)                                  :: flow_lake(nbp_mpi)          
    
    ! diags
    REAL(r_std)                                  :: fast_diag_mpi(nbp_mpi)  
    REAL(r_std)                                  :: slow_diag_mpi(nbp_mpi)  
    REAL(r_std)                                  :: stream_diag_mpi(nbp_mpi)  
    REAL(r_std)                                  :: hydrographs_diag_mpi(nbp_mpi)  
   
    ! from input model -> routing_grid
    REAL(r_std)                      :: runoff_r(nbpt_r)              !! Grid-point runoff (kg/m^2/dt)
    REAL(r_std)                      :: drainage_r(nbpt_r)            !! Grid-point drainage (kg/m^2/dt)

    REAL(r_std), DIMENSION(nbpt_r)        :: fast_flow_r                 !! Outflow from the fast reservoir (kg/dt)
    REAL(r_std), DIMENSION(nbpt_r)        :: slow_flow_r                 !! Outflow from the slow reservoir (kg/dt)
    REAL(r_std), DIMENSION(nbpt_r)        :: stream_flow_r               !! Outflow from the stream reservoir (kg/dt)
    REAL(r_std), DIMENSION(nbpt_r)        :: hydrographs_r                !! hydrograph (kg/dt)
    REAL(r_std), DIMENSION(nbpt_r)        :: transport_r                 !! Water transport between basins (kg/dt)

    INTEGER(i_std)                               :: ig        !! Indices (unitless)
    INTEGER(i_std)                               :: isplit
    
    LOGICAL, PARAMETER                           :: check_reservoir = .TRUE. !! Logical to choose if we write informations when a negative amount of water is occurring in a reservoir (true/false)

    REAL(r_std), DIMENSION(nbpt_r)        :: lakeinflow_r   
    REAL(r_std), DIMENSION(nbpt_r)        :: coastalflow_r   
    REAL(r_std), DIMENSION(nbpt_r)        :: riverflow_r   

    REAL(r_std), DIMENSION(nbpt_r)        :: flow_r   
    REAL(r_std), DIMENSION(nbpt_rp1)      :: flow_rp1   
    REAL(r_std)                           :: flow                      !! Outflow computation for the reservoirs (kg/dt)
    REAL(r_std)                           :: riverflow_omp(nbpt)
    REAL(r_std)                           :: coastalflow_omp(nbpt)
    REAL(r_std)                           :: lakeinflow_omp(nbpt)
    REAL(r_std) :: basins_riverflow_mpi(0:basins_count)
    REAL(r_std) :: basins_riverflow(0:basins_count)
    REAL(r_std) :: water_balance_before, water_balance_after
    REAL(r_std) :: sum_water_before, sum_water_after
    REAL(r_std)  :: value
    INTEGER(i_std) :: k
    INTEGER :: ierr

    !_ ================================================================================================================================
    !> The outflow fluxes from the three reservoirs are computed. 
    !> The outflow of volume of water Vi into the reservoir i is assumed to be linearly related to its volume.
    !> The water travel simulated by the routing scheme is dependent on the water retention index topo_resid
    !> given by a 0.5 degree resolution map for each pixel performed from a simplification of Manning's formula
    !> (Dingman, 1994; Ducharne et al., 2003).
    !> The resulting product of tcst (in day/m) and topo_resid (in m) represents the time constant (day)
    !> which is an e-folding time, the time necessary for the water amount
    !> in the stream reservoir to decrease by a factor e. Hence, it gives an order of
    !> magnitude of the travel time through this reservoir between
    !> the sub-basin considered and its downstream neighbor.



    CALL gather_omp(runoff_mean,runoff)
    CALL gather_omp(drainage_mean, drainage)
    CALL gather_omp(area, area_mpi)
    CALL gather_omp(contfrac, contfrac_mpi)

    IF (is_omp_root) THEN
      CALL xios_send_field("routing_basins_area",basins_area_r)

      hydrographs_r(:)=0
      ! water balance before

      sum_water_before = sum((runoff(:)+drainage(:))*routing_weight(:))+sum(fast_reservoir_r(:)+slow_reservoir_r(:)+stream_reservoir_r(:))
      CALL MPI_ALLREDUCE(sum_water_before, water_balance_before,1,MPI_REAL_ORCH,MPI_SUM,MPI_COMM_ORCH,ierr)


      runoff_in=runoff*routing_weight_in
      CALL xios_send_field("routing_runoff",runoff_in)   ! interp conservative model -> routing
      CALL xios_recv_field("routing_runoff_r",runoff_r)
      WHERE(.NOT. routing_mask_r) runoff_r = 0

      drainage_in=drainage*routing_weight_in
      CALL xios_send_field("routing_drainage",drainage_in)   ! interp conservative model -> routing
      CALL xios_recv_field("routing_drainage_r",drainage_r)
      WHERE(.NOT. routing_mask_r) drainage_r = 0

      CALL MPI_ALLREDUCE(sum(runoff*(routing_weight-unrouted_weight)),sum_water_before,1,MPI_REAL_ORCH,MPI_SUM,MPI_COMM_ORCH,ierr)
      CALL MPI_ALLREDUCE(sum(runoff_r),sum_water_after,1,MPI_REAL_ORCH,MPI_SUM,MPI_COMM_ORCH,ierr)
      IF (sum_water_after/=0) THEN
        IF (is_mpi_root .AND. printlev_loc >= 3) WRITE (numout,*) "runoff fixer : ", sum_water_before/sum_water_after
        DO ig=1,nbpt_r
          runoff_r(ig) =  runoff_r(ig) * sum_water_before/sum_water_after
        ENDDO
      ENDIF

      CALL MPI_ALLREDUCE(sum(drainage*(routing_weight-unrouted_weight)),sum_water_before,1,MPI_REAL_ORCH,MPI_SUM,MPI_COMM_ORCH,ierr)
      CALL MPI_ALLREDUCE(sum(drainage_r),sum_water_after,1,MPI_REAL_ORCH,MPI_SUM,MPI_COMM_ORCH,ierr)
      IF (sum_water_after/=0) THEN
        IF (is_mpi_root .AND. printlev_loc >= 3) WRITE (numout,*) "drainage fixer : ", sum_water_before/sum_water_after
        DO ig=1,nbpt_r
          drainage_r(ig) =  drainage_r(ig) * sum_water_before/sum_water_after
        ENDDO
      ENDIF
          

      CALL MPI_ALLREDUCE(sum((runoff+drainage)*(routing_weight-unrouted_weight)),sum_water_before,1,MPI_REAL_ORCH,MPI_SUM,MPI_COMM_ORCH,ierr)
      CALL MPI_ALLREDUCE(sum(runoff_r+drainage_r),sum_water_after,1,MPI_REAL_ORCH,MPI_SUM,MPI_COMM_ORCH,ierr)

      IF (is_mpi_root .AND. printlev_loc >= 3) WRITE (numout,*) "Loose water by interpolation ;  before : ", sum_water_before," ; after : ",sum_water_after,  &
                                 " ; delta : ", 100.*(sum_water_after-sum_water_before)/(0.5*(sum_water_after+sum_water_before)),"%"


      runoff_in(:) = runoff(:)*unrouted_weight(:)
      drainage_in(:) = drainage(:)*unrouted_weight(:)
      runoff_r(:) = runoff_r(:) / split_routing
      drainage_r(:) = drainage_r(:) / split_routing
      hydrographs_r(:) = 0
      transport_r(:)=0

      DO isplit=1,split_routing

        DO ig=1,nbpt_r
          IF ( routing_mask_r(ig) ) THEN

            flow = MIN(fast_reservoir_r(ig)/((topoind_r(ig)*topoind_factor/1000.)*fast_tcst*one_day/(dt_routing/split_routing)),&
                   & fast_reservoir_r(ig)-min_sechiba)
            fast_flow_r(ig) = MAX(flow, zero)
          !
            flow = MIN(slow_reservoir_r(ig)/((topoind_r(ig)*topoind_factor/1000.)*slow_tcst*one_day/(dt_routing/split_routing)),&
                   & slow_reservoir_r(ig)-min_sechiba)
            slow_flow_r(ig) = MAX(flow, zero)
          !
            flow = MIN(stream_reservoir_r(ig)/((topoind_r(ig)*topoind_factor/1000.)*stream_tcst*one_day/(dt_routing/split_routing)),&
                   & stream_reservoir_r(ig)-min_sechiba)
            stream_flow_r(ig) = MAX(flow, zero)
          ! 
          !
          ELSE
            fast_flow_r(ig) = zero
            slow_flow_r(ig) = zero
            stream_flow_r(ig) = zero
          ENDIF
        ENDDO


        !-
        !- Compute the transport
        !-

        DO ig=1,nbpt_r
          flow_rp1(r_to_rp1(ig))=fast_flow_r(ig) + slow_flow_r(ig) + stream_flow_r(ig)
        ENDDO

        CALL update_halo(flow_rp1) ! transfert halo

        DO ig=1,nbpt_rp1
          IF ( route_flow_rp1(ig) > 0 ) THEN
            transport_r(route_flow_rp1(ig))=transport_r(route_flow_rp1(ig))+ flow_rp1(ig)
          ENDIF
        ENDDO


        DO ig=1,nbpt_r
          IF ( routing_mask_r(ig) ) THEN
            fast_reservoir_r(ig) =  fast_reservoir_r(ig) + runoff_r(ig) - fast_flow_r(ig)
            slow_reservoir_r(ig) = slow_reservoir_r(ig) + drainage_r(ig) -  slow_flow_r(ig)
            stream_reservoir_r(ig) = stream_reservoir_r(ig) - stream_flow_r(ig)
            IF (is_streamflow_r(ig)) THEN
              stream_reservoir_r(ig)= stream_reservoir_r(ig) +  transport_r(ig)
              transport_r(ig) = 0
            ENDIF  

            IF ( stream_reservoir_r(ig) .LT. zero ) THEN
              IF ( check_reservoir ) THEN
                WRITE(numout,*) "WARNING : negative stream reservoir at :", ig, ". Problem is being corrected."
                WRITE(numout,*) "stream_reservoir, transport, stream_flow,: ", &
                                   stream_reservoir_r(ig), transport_r(ig), stream_flow_r(ig)
              ENDIF
              fast_reservoir_r(ig) =  fast_reservoir_r(ig) + stream_reservoir_r(ig)
              stream_reservoir_r(ig) = zero
            ENDIF
            !
            IF ( fast_reservoir_r(ig) .LT. zero ) THEN
              IF ( check_reservoir ) THEN
                WRITE(numout,*) "WARNING : negative fast reservoir at :", ig, ". Problem is being corrected."
                WRITE(numout,*) "fast_reservoir, runoff, fast_flow : ", fast_reservoir_r(ig), &
                                 &runoff_r(ig), fast_flow_r(ig)
              ENDIF
              slow_reservoir_r(ig) =  slow_reservoir_r(ig) + fast_reservoir_r(ig)
              fast_reservoir_r(ig) = zero
            ENDIF

            IF ( slow_reservoir_r(ig) .LT. - min_sechiba ) THEN
              WRITE(numout,*) 'WARNING : There is a negative reservoir at :', ig
              WRITE(numout,*) 'WARNING : slowr, slow_flow, drainage', &
                               & slow_reservoir_r(ig), slow_flow_r(ig), drainage_r(ig)
              CALL ipslerr_p(2, 'routing_interp_topo_flow', 'WARNING negative slow_reservoir.','','')
            ENDIF
          ENDIF
        ENDDO

        DO ig=1,nbpt_r
          IF ( routing_mask_r(ig) ) THEN
            hydrographs_r(ig)=hydrographs_r(ig)+fast_flow_r(ig)+slow_flow_r(ig)+stream_flow_r(ig)
          ENDIF
        ENDDO

      ENDDO ! isplit

      lakeinflow_r(:)=0
      coastalflow_r(:)=0
      riverflow_r(:)=0
      basins_riverflow_mpi(:)=0

      DO ig=1,nbpt_r
        IF ( routing_mask_r(ig) ) THEN
          IF (is_lakeinflow_r(ig))  THEN
            lakeinflow_r(ig) = transport_r(ig)
            basins_riverflow_mpi(basins_extended_r(ig)) = basins_riverflow_mpi(basins_extended_r(ig))+lakeinflow_r(ig)
          ENDIF

          IF (is_coastalflow_r(ig)) THEN
            coastalflow_r(ig) = transport_r(ig)
            basins_riverflow_mpi(basins_extended_r(ig)) = &
            basins_riverflow_mpi(basins_extended_r(ig))+coastalflow_r(ig)
          ENDIF

          IF (is_riverflow_r(ig)) THEN
            riverflow_r(ig) = transport_r(ig)
            basins_riverflow_mpi(basins_extended_r(ig)) = &
            basins_riverflow_mpi(basins_extended_r(ig))+riverflow_r(ig)
          ENDIF
        ENDIF
      ENDDO

      ! now send riverflow, coastalflow and lakeinflow on orchidee grid

      coastalflow(:)=0.
      riverflow(:)=0.
      lakeinflow(:)=0.

      CALL xios_send_field("routing_coastalflow_to_coast_r" ,coastalflow_r*weight_coast_to_coast_r)
      CALL xios_recv_field("routing_coastalflow_to_coast" ,flow_coast)
      WHERE(.NOT.is_coastline) flow_coast=0

      CALL xios_send_field("routing_coastalflow_to_lake_r" ,coastalflow_r*weight_coast_to_lake_r)
      CALL xios_recv_field("routing_coastalflow_to_lake" ,flow_lake)
      WHERE(is_coastline) flow_lake=0.

      CALL MPI_ALLREDUCE(sum(coastalflow_r),sum_water_before,1,MPI_REAL_ORCH,MPI_SUM,MPI_COMM_ORCH,ierr)
      CALL MPI_ALLREDUCE(sum(flow_coast+flow_lake),sum_water_after,1,MPI_REAL_ORCH,MPI_SUM,MPI_COMM_ORCH,ierr)
      IF (sum_water_after/=0) THEN
        IF (is_mpi_root .AND. printlev_loc >= 3) WRITE (numout,*) "coastalflow fixer : ", sum_water_before/sum_water_after
        flow_coast(:)=flow_coast(:)*(sum_water_before/sum_water_after)
        flow_lake(:)=flow_lake(:)*(sum_water_before/sum_water_after)
      ENDIF
       
      coastalflow=coastalflow+flow_coast
      lakeinflow=lakeinflow+flow_lake

      CALL xios_send_field("routing_riverflow_to_coast_r" ,riverflow_r*weight_coast_to_coast_r)
      CALL xios_recv_field("routing_riverflow_to_coast" ,flow_coast)
      WHERE(.NOT.is_coastline) flow_coast=0

      CALL xios_send_field("routing_riverflow_to_lake_r" ,riverflow_r*weight_coast_to_lake_r)
      CALL xios_recv_field("routing_riverflow_to_lake" ,flow_lake)
      WHERE(is_coastline) flow_lake=0.
      
      CALL MPI_ALLREDUCE(sum(riverflow_r),sum_water_before,1,MPI_REAL_ORCH,MPI_SUM,MPI_COMM_ORCH,ierr)
      CALL MPI_ALLREDUCE(sum(flow_coast+flow_lake),sum_water_after,1,MPI_REAL_ORCH,MPI_SUM,MPI_COMM_ORCH,ierr)
      IF (sum_water_after/=0) THEN
        IF (is_mpi_root .AND. printlev_loc >= 3) WRITE (numout,*) "riverflow fixer : ", sum_water_before/sum_water_after
        flow_coast(:)=flow_coast(:)*(sum_water_before/sum_water_after)
        flow_lake(:)=flow_lake(:)*(sum_water_before/sum_water_after)
      ENDIF   
       
      riverflow=riverflow+flow_coast
      lakeinflow=lakeinflow+flow_lake

      CALL xios_send_field("routing_lakeinflow_to_coast_r" ,lakeinflow_r*weight_lake_to_coast_r)
      CALL xios_recv_field("routing_lakeinflow_to_coast" ,flow_coast)
      WHERE(.NOT.is_coastline) flow_coast=0

      CALL xios_send_field("routing_lakeinflow_to_lake_r" ,lakeinflow_r*weight_lake_to_lake_r)
      CALL xios_recv_field("routing_lakeinflow_to_lake" ,flow_lake)
      WHERE(is_coastline) flow_lake=0.

      CALL MPI_ALLREDUCE(sum(lakeinflow_r),sum_water_before,1,MPI_REAL_ORCH,MPI_SUM,MPI_COMM_ORCH,ierr)
      CALL MPI_ALLREDUCE(sum(flow_coast+flow_lake),sum_water_after,1,MPI_REAL_ORCH,MPI_SUM,MPI_COMM_ORCH,ierr)
      IF (sum_water_after/=0) THEN
        IF (is_mpi_root.AND. printlev_loc >= 3) WRITE (numout,*) "lakeinflow fixer : ", sum_water_before/sum_water_after
        flow_coast(:)=flow_coast(:)*(sum_water_before/sum_water_after)
        flow_lake(:)=flow_lake(:)*(sum_water_before/sum_water_after)
      ENDIF   
       
      coastalflow=coastalflow+flow_coast+(runoff+drainage)*unrouted_weight
      lakeinflow=lakeinflow+flow_lake


!      WHERE(is_coastline) coastalflow = coastalflow + runoff_in + drainage_in
!      WHERE(.NOT.is_coastline) lakeinflow = lakeinflow + runoff_in + drainage_in

      sum_water_after = sum(coastalflow(:) + riverflow(:) + lakeinflow(:))+sum(fast_reservoir_r(:)+slow_reservoir_r(:)+stream_reservoir_r(:))
      CALL MPI_ALLREDUCE(sum_water_after, water_balance_after,1,MPI_REAL_ORCH,MPI_SUM,MPI_COMM_ORCH,ierr)
      
      IF (is_mpi_root .AND. printlev_loc >= 3) WRITE (numout,*) "routing flow water Balance ;  before : ", water_balance_before," ; after : ",water_balance_after,  &
                               " ; delta : ", 100*(water_balance_after-water_balance_before)/(0.5*(water_balance_after+water_balance_before)),"%"
  ! diag
      CALL xios_send_field("routing_fast_reservoir_r"  , fast_reservoir_r)
      CALL xios_send_field("routing_slow_reservoir_r"  , slow_reservoir_r)
      CALL xios_send_field("routing_stream_reservoir_r"  , stream_reservoir_r)
      CALL xios_send_field("routing_riverflow_r"  , riverflow_r)
      CALL xios_send_field("routing_coastalflow_r"  , coastalflow_r)
      CALL xios_send_field("routing_lakeinflow_r"  , lakeinflow_r)
      CALL xios_send_field("out_flow",lakeinflow+coastalflow+riverflow)
      CALL xios_send_field("routing_hydrographs_r", (hydrographs_r+lakeinflow_r+coastalflow_r+riverflow_r)/1000./dt_routing)
      CALL xios_send_field("routing_riverflow"  , riverflow)
      CALL xios_send_field("routing_coastalflow"  , coastalflow)
      CALL xios_send_field("routing_lakeinflow"  , lakeinflow)
!
!     stations
      IF (nb_station/=0) THEN
        CALL xios_orchidee_change_context("orchidee_routing_out")
        station_ts = station_ts +1
        CALL xios_update_calendar(station_ts)
        DO k=1,nb_station
          IF (station_index(k) /=-1) THEN
            value = hydrographs_r(station_index(k))/1000./dt_routing
          ELSE
            value = 0
          ENDIF
          CALL MPI_ALLREDUCE(MPI_IN_PLACE,value,1,MPI_REAL_ORCH,MPI_SUM,MPI_COMM_ORCH,ierr)
          CALL xios_send_field(TRIM(station(k)),value)
        ENDDO

        CALL xios_orchidee_change_context("orchidee")  !! return on orchidee context
      END IF !nb_station/=0
      
      !!! for orchidee history diag

      CALL xios_send_field("routing_fast_diag_r", fast_reservoir_r)
      CALL xios_recv_field("routing_fast_diag", fast_diag_mpi)
      CALL xios_send_field("routing_slow_diag_r", slow_reservoir_r)
      CALL xios_recv_field("routing_slow_diag", slow_diag_mpi)
      CALL xios_send_field("routing_stream_diag_r", stream_reservoir_r)
      CALL xios_recv_field("routing_stream_diag", stream_diag_mpi)
      CALL xios_send_field("routing_hydrographs_diag_r", hydrographs_r+lakeinflow_r+coastalflow_r+riverflow_r)
      CALL xios_recv_field("routing_hydrographs_diag", hydrographs_diag_mpi)

      fast_diag_mpi=fast_diag_mpi/(area_mpi*contfrac_mpi)      !! kg => kg/m^2
      slow_diag_mpi=slow_diag_mpi/(area_mpi*contfrac_mpi)      !! kg => kg/m^2
      stream_diag_mpi=stream_diag_mpi/(area_mpi*contfrac_mpi)  !! kg => kg/m^2
      hydrographs_diag_mpi = hydrographs_diag_mpi/1000.

    ENDIF ! is_omp_root

    CALL scatter_omp(riverflow, riverflow_omp)
    riverflow_mean = riverflow_mean + riverflow_omp
    CALL scatter_omp(coastalflow, coastalflow_omp)
    coastalflow_mean = coastalflow_mean + coastalflow_omp
    CALL scatter_omp(lakeinflow, lakeinflow_omp)
    lakeinflow_mean = lakeinflow_mean + lakeinflow_omp
   
    CALL scatter_omp(fast_diag_mpi, fast_diag)
    CALL scatter_omp(slow_diag_mpi, slow_diag)
    CALL scatter_omp(stream_diag_mpi, stream_diag)
    CALL scatter_omp(hydrographs_diag_mpi, hydrographs_diag)
    
    CALL routing_flow_reset_mean
    integration_time = integration_time + dt_routing
  END SUBROUTINE routing_flow_main

  SUBROUTINE routing_flow_retreive_outflow(dt, riverflow, coastalflow, lakeinflow)
  IMPLICIT NONE
    REAL(r_std),INTENT(IN) :: dt                    ! part of time over retreiving outflow mean
    REAL(r_std),INTENT(OUT) :: riverflow(nbpt)
    REAL(r_std),INTENT(OUT) :: coastalflow(nbpt)
    REAL(r_std),INTENT(OUT) :: lakeinflow(nbpt)
    REAL(r_std)            :: fract

    fract = dt/integration_time
    riverflow   = riverflow_mean * fract 
    coastalflow = coastalflow_mean * fract
    lakeinflow  = lakeinflow_mean * fract

    riverflow_mean = riverflow_mean - riverflow
    coastalflow_mean = coastalflow_mean - coastalflow
    lakeinflow_mean = lakeinflow_mean - lakeinflow
    
    integration_time = integration_time-dt

  END SUBROUTINE routing_flow_retreive_outflow


  FUNCTION routing_flow_get_water_balance(contfrac)
  USE grid, ONLY : area
  USE mod_orchidee_para, ONLY : is_omp_root, is_mpi_root
  IMPLICIT NONE
    REAL(r_std),INTENT(IN) :: contfrac(nbpt)
    REAL(r_std) :: routing_flow_get_water_balance   !  Amount of water (local) in kg
    REAL(r_std) :: water
    
    water = 0
    water = water + SUM((runoff_mean + drainage_mean)*contfrac*area)
    IF (is_omp_root) THEN
      water = water + SUM(fast_reservoir_r + slow_reservoir_r + stream_reservoir_r)
    ENDIF
    water = water + SUM(riverflow_mean + coastalflow_mean + lakeinflow_mean) 
    routing_flow_get_water_balance = water
  
  END FUNCTION routing_flow_get_water_balance

  SUBROUTINE routing_flow_diags(dt_sechiba)
  USE xios
  IMPLICIT NONE
    REAL(r_std) :: dt_sechiba

    CALL xios_orchidee_send_field("fastr",fast_diag)  
    CALL xios_orchidee_send_field("slowr",slow_diag)
    CALL xios_orchidee_send_field("streamr",stream_diag)
    CALL xios_orchidee_send_field("hydrographs",hydrographs_diag*dt_sechiba)

  END SUBROUTINE routing_flow_diags


  !!  =============================================================================================================================
  !! SUBROUTINE:         routing_flow_finalize
  !!
  !>\BRIEF	         Write to restart file
  !!
  !! DESCRIPTION:        Write module variables to restart file
  !!
  !! RECENT CHANGE(S)
  !!
  !! REFERENCE(S)
  !! 
  !! FLOWCHART   
  !! \n
  !_ ==============================================================================================================================

  SUBROUTINE routing_flow_finalize(kjit, rest_id)
    USE xios
    USE ioipsl
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: kjit
    INTEGER, INTENT(IN) :: rest_id
    !_ ================================================================================================================================

    IF (is_omp_root) THEN
       CALL xios_send_field("fast_reservoir_restart",fast_reservoir_r)
       CALL xios_send_field("slow_reservoir_restart",slow_reservoir_r)
       CALL xios_send_field("stream_reservoir_restart",stream_reservoir_r)
       CALL xios_send_field("integration_time_restart",integration_time)
    ENDIF
    
    CALL restput_p (rest_id, 'riverflow', nbp_glo, 1, 1, kjit, riverflow_mean, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'coastalflow', nbp_glo, 1, 1, kjit, coastalflow_mean, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'lakeinflow', nbp_glo, 1, 1, kjit, lakeinflow_mean, 'scatter',  nbp_glo, index_g)
    CALL routing_flow_finalize_mean(kjit, rest_id)
    
    CALL routing_flow_finalize_diag(kjit, rest_id)

    ! Station context
    IF (nb_station/=0) THEN
      CALL xios_orchidee_change_context("orchidee_routing_out")
      CALL xios_orchidee_context_finalize()
      CALL xios_orchidee_change_context("orchidee")
    ENDIF
 
  END SUBROUTINE routing_flow_finalize

  SUBROUTINE routing_flow_finalize_diag(kjit, rest_id)
    USE ioipsl
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: kjit
    INTEGER, INTENT(IN) :: rest_id
    
     CALL restput_p (rest_id, 'hydrographs', nbp_glo, 1, 1, kjit, hydrographs_diag, 'scatter',  nbp_glo, index_g)

  END SUBROUTINE routing_flow_finalize_diag

  !! ================================================================================================================================
  !! SUBROUTINE 	: routing_flow_clear
  !!
  !>\BRIEF         This subroutine deallocates the block memory previously allocated.
  !!
  !! DESCRIPTION:  This subroutine deallocates the block memory previously allocated.
  !!
  !! RECENT CHANGE(S): None
  !!
  !! MAIN OUTPUT VARIABLE(S):
  !!
  !! REFERENCES   : None
  !!
  !! FLOWCHART    :None
  !! \n
  !_ ================================================================================================================================

  SUBROUTINE routing_flow_clear
  IMPLICIT NONE

    IF (is_omp_root) THEN
       IF (ALLOCATED(topoind_r)) DEALLOCATE(topoind_r)
       IF (ALLOCATED(route_flow_rp1)) DEALLOCATE(route_flow_rp1)
       IF (ALLOCATED(routing_mask_r)) DEALLOCATE(routing_mask_r)
       IF (ALLOCATED(fast_reservoir_r)) DEALLOCATE(fast_reservoir_r)
       IF (ALLOCATED(slow_reservoir_r)) DEALLOCATE(slow_reservoir_r)
       IF (ALLOCATED(is_lakeinflow_r)) DEALLOCATE(is_lakeinflow_r) 
       IF (ALLOCATED(is_coastalflow_r)) DEALLOCATE(is_coastalflow_r)    
       IF (ALLOCATED(is_riverflow_r)) DEALLOCATE(is_riverflow_r)   
    ENDIF

  END SUBROUTINE routing_flow_clear

#endif
END MODULE routing_interp_topo_flow_mod
