! =================================================================================================================================
! MODULE       : routing_interp_topo_v0
!
!  CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       This module routes the water over the continents into the oceans and computes the water
!!             stored in floodplains or taken for irrigation.
!!
!!\n DESCRIPTION: The subroutines in this subroutine is only called when ROUTING_METHOD=interp_topo_v0 is set in run.def.
!!                The method can be used for regular latitude-longitude grid or for unstructured grid.
!!                It computes the routing on a separate grid and exchanges data with the ORCHIDEE grid using interpolations
!!                performed via XIOS.
!!
!! RECENT CHANGE(S): July 2024: Changed routing method's name "simple" to "interp_topo_v0".
!!
!! REFERENCE(S)	:
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_sechiba/routing_interp_topo_v0.f90 $
!! $Date: 2025-09-16 12:13:09 +0200 (mar. 16 sept. 2025) $
!! $Revision: 9095 $
!! \n
!_ ================================================================================================================================

 

! Histoire Salee
!---------------
! La douce riviere
! Sortant de son lit
! S'est jetee ma chere
! dans les bras mais oui
! du beau fleuve
!
! L'eau coule sous les ponts
! Et puis les flots s'emeuvent
! - N'etes vous pas au courant ?
! Il parait que la riviere 
! Va devenir mer
!                       Roland Bacri
!


MODULE routing_interp_topo_v0
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
  USE mod_orchidee_para


  IMPLICIT NONE
  PRIVATE
  PUBLIC :: routing_interp_topo_v0_main, routing_interp_topo_v0_xios_initialize
  PUBLIC :: routing_interp_topo_v0_initialize, routing_interp_topo_v0_finalize, routing_interp_topo_v0_clear

  !! PARAMETERS
  INTEGER(i_std), PARAMETER                 :: nbasmax=5                   !! The maximum number of basins we wish to have per grid box (truncation of the model) (unitless)
  INTEGER(i_std), SAVE                      :: nbvmax                      !! The maximum number of basins we can handle at any time during the generation of the maps (unitless)
  !$OMP THREADPRIVATE(nbvmax)
  REAL(r_std), SAVE                         :: fast_tcst = 3.0             !! Property of the fast reservoir (day/m)
  !$OMP THREADPRIVATE(fast_tcst)
  REAL(r_std), SAVE                         :: slow_tcst = 25.0            !! Property of the slow reservoir (day/m)
  !$OMP THREADPRIVATE(slow_tcst)
  REAL(r_std), SAVE                         :: stream_tcst = 0.24          !! Property of the stream reservoir (day/m)
  !$OMP THREADPRIVATE(stream_tcst)
  REAL(r_std), SAVE                         :: flood_tcst = 4.0            !! Property of the floodplains reservoir (day/m)
  !$OMP THREADPRIVATE(flood_tcst)
  REAL(r_std), SAVE                         :: swamp_cst = 0.2             !! Fraction of the river transport that flows to the swamps (unitless;0-1)
  !$OMP THREADPRIVATE(swamp_cst)
  !
  !  Relation between volume and fraction of floodplains
  !
  REAL(r_std), SAVE                         :: beta = 2.0                  !! Parameter to fix the shape of the floodplain (>1 for convex edges, <1 for concave edges) (unitless)
  !$OMP THREADPRIVATE(beta)
  REAL(r_std), SAVE                         :: betap = 0.5                 !! Ratio of the basin surface intercepted by ponds and the maximum surface of ponds (unitless;0-1)
  !$OMP THREADPRIVATE(betap)
  REAL(r_std), SAVE                         :: floodcri = 2000.0           !! Potential height for which all the basin is flooded (mm)
  !$OMP THREADPRIVATE(floodcri)
  !
  !  Relation between maximum surface of ponds and basin surface, and drainage (mm/j) to the slow_res
  !
  REAL(r_std), PARAMETER                    :: pond_bas = 50.0             !! [DISPENSABLE] - not used
  REAL(r_std), SAVE                         :: pondcri = 2000.0            !! Potential height for which all the basin is a pond (mm)
  !$OMP THREADPRIVATE(pondcri)

  REAL(r_std), PARAMETER                    :: maxevap_lake = 7.5/86400.   !! Maximum evaporation rate from lakes (kg/m^2/s)
  REAL(r_std),SAVE                          :: dt_routing                  !! Routing time step (s)
  !$OMP THREADPRIVATE(dt_routing)
  INTEGER(i_std), SAVE                      :: diagunit = 87               !! Diagnostic file unit (unitless)
  !$OMP THREADPRIVATE(diagunit)

  ! Logicals to control model configuration
  !
  LOGICAL, SAVE                             :: dofloodinfilt = .FALSE.     !! Logical to choose if floodplains infiltration is activated or not (true/false)
  !$OMP THREADPRIVATE(dofloodinfilt)
  LOGICAL, SAVE                             :: doswamps = .FALSE.          !! Logical to choose if swamps are activated or not (true/false)
  !$OMP THREADPRIVATE(doswamps)
  LOGICAL, SAVE                             :: doponds = .FALSE.           !! Logical to choose if ponds are activated or not (true/false)
  !$OMP THREADPRIVATE(doponds)

  INTEGER(i_std), SAVE                                       :: num_largest                 !! Number of largest river basins which should be treated as independently as rivers
                                                                                            !! (not flow into ocean as diffusion coastal flow) (unitless)
  !$OMP THREADPRIVATE(num_largest)
  REAL(r_std), SAVE                                          :: time_counter                !! Time counter (s)
  !$OMP THREADPRIVATE(time_counter)
  REAL(r_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:,:)     :: routing_area_loc            !! Surface of basin (m^2)
  !$OMP THREADPRIVATE(routing_area_loc)
  REAL(r_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:,:)     :: topo_resid_loc              !! Topographic index of the retention time (m)
  !$OMP THREADPRIVATE(topo_resid_loc)
  INTEGER(i_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:,:)  :: route_togrid_loc            !! Grid into which the basin flows (unitless)
  !$OMP THREADPRIVATE(route_togrid_loc)
  INTEGER(i_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:,:)  :: route_tobasin_loc           !! Basin in to which the water goes (unitless)
  !$OMP THREADPRIVATE(route_tobasin_loc)
  INTEGER(i_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:,:)  :: route_nbintobas_loc         !! Number of basin into current one (unitless)
  !$OMP THREADPRIVATE(route_nbintobas_loc)
  INTEGER(i_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:,:)  :: global_basinid_loc          !! ID of basin (unitless)
  !$OMP THREADPRIVATE(global_basinid_loc)
  INTEGER(i_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:,:)  :: hydrodiag_loc               !! Variable to diagnose the hydrographs
  !$OMP THREADPRIVATE(hydrodiag_loc)
  REAL(r_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:)       :: hydroupbasin_loc            !! The area upstream of the gauging station (m^2)
  !$OMP THREADPRIVATE(hydroupbasin_loc)

  ! parallelism
  REAL(r_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:,:)     :: routing_area_glo            !! Surface of basin (m^2)
  !$OMP THREADPRIVATE(routing_area_glo)
  REAL(r_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:,:)     :: topo_resid_glo              !! Topographic index of the retention time (m)
  !$OMP THREADPRIVATE(topo_resid_glo)
  INTEGER(i_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:,:)  :: route_togrid_glo            !! Grid into which the basin flows (unitless)
  !$OMP THREADPRIVATE(route_togrid_glo)
  INTEGER(i_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:,:)  :: route_tobasin_glo           !! Basin in to which the water goes (unitless)
  !$OMP THREADPRIVATE(route_tobasin_glo)
  INTEGER(i_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:,:)  :: route_nbintobas_glo         !! Number of basin into current one (unitless)
  !$OMP THREADPRIVATE(route_nbintobas_glo)
  INTEGER(i_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:,:)  :: global_basinid_glo          !! ID of basin (unitless)
  !$OMP THREADPRIVATE(global_basinid_glo)
  INTEGER(i_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:,:)  :: hydrodiag_glo               !! Variable to diagnose the hydrographs
  !$OMP THREADPRIVATE(hydrodiag_glo)
  REAL(r_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:)       :: hydroupbasin_glo            !! The area upstream of the gauging station (m^2)
  !$OMP THREADPRIVATE(hydroupbasin_glo)

  REAL(r_std), SAVE, POINTER, DIMENSION(:,:)                 :: routing_area                !! Surface of basin (m^2)
  !$OMP THREADPRIVATE(routing_area)
  REAL(r_std), SAVE, POINTER, DIMENSION(:,:)                 :: topo_resid                  !! Topographic index of the retention time (m)
  !$OMP THREADPRIVATE(topo_resid)
  INTEGER(i_std), SAVE, POINTER, DIMENSION(:,:)              :: route_togrid                !! Grid into which the basin flows (unitless)
  !$OMP THREADPRIVATE(route_togrid)
  INTEGER(i_std), SAVE, POINTER, DIMENSION(:,:)              :: route_tobasin               !! Basin in to which the water goes (unitless)
  !$OMP THREADPRIVATE(route_tobasin)
  INTEGER(i_std), SAVE, POINTER, DIMENSION(:,:)              :: route_nbintobas             !! Number of basin into current one (unitless)
  !$OMP THREADPRIVATE(route_nbintobas)
  INTEGER(i_std), SAVE, POINTER, DIMENSION(:,:)              :: global_basinid              !! ID of basin (unitless)
  !$OMP THREADPRIVATE(global_basinid)
  INTEGER(i_std), SAVE, POINTER, DIMENSION(:,:)              :: hydrodiag                   !! Variable to diagnose the hydrographs
  !$OMP THREADPRIVATE(hydrodiag)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: slowflow_diag               !! Diagnostic slow flow hydrographs (kg/dt)
  !$OMP THREADPRIVATE(slowflow_diag)  
  REAL(r_std), SAVE, POINTER, DIMENSION(:)                   :: hydroupbasin                !! The area upstream of the gauging station (m^2)
  !$OMP THREADPRIVATE(hydroupbasin)

  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: irrigated                   !! Area equipped for irrigation in each grid box (m^2)
  !$OMP THREADPRIVATE(irrigated)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: floodplains                 !! Maximal surface which can be inundated in each grid box (m^2)
  !$OMP THREADPRIVATE(floodplains)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: swamp                       !! Maximal surface of swamps in each grid box (m^2)
  !$OMP THREADPRIVATE(swamp)
  !
  ! The reservoirs, variables kept in the restart file.
  !
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)             :: fast_reservoir              !! Water amount in the fast reservoir (kg)
  !$OMP THREADPRIVATE(fast_reservoir)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)             :: slow_reservoir              !! Water amount in the slow reservoir (kg)
  !$OMP THREADPRIVATE(slow_reservoir)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)             :: stream_reservoir            !! Water amount in the stream reservoir (kg)
  !$OMP THREADPRIVATE(stream_reservoir)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)             :: flood_reservoir             !! Water amount in the floodplains reservoir (kg)
  !$OMP THREADPRIVATE(flood_reservoir)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: lake_reservoir              !! Water amount in the lake reservoir (kg)
  !$OMP THREADPRIVATE(lake_reservoir)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: pond_reservoir              !! Water amount in the pond reservoir (kg)
  !$OMP THREADPRIVATE(pond_reservoir)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)             :: flood_frac_bas              !! Flooded fraction per basin (unitless;0-1)
  !$OMP THREADPRIVATE(flood_frac_bas)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: pond_frac                   !! Pond fraction per grid box (unitless;0-1)
  !$OMP THREADPRIVATE(pond_frac)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: flood_height                !! Floodplain height (mm)
  !$OMP THREADPRIVATE(flood_height)
  !
  ! The accumulated fluxes.
  !
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: floodout_mean               !! Accumulated flow out of floodplains (kg/m^2/dt)
  !$OMP THREADPRIVATE(floodout_mean)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: runoff_mean                 !! Accumulated runoff (kg/m^2/dt)
  !$OMP THREADPRIVATE(runoff_mean)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: drainage_mean               !! Accumulated drainage (kg/m^2/dt)
  !$OMP THREADPRIVATE(drainage_mean)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: transpot_mean               !! Mean potential transpiration from the plants (kg/m^2/dt)
  !$OMP THREADPRIVATE(transpot_mean)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: precip_mean                 !! Accumulated precipitation (kg/m^2/dt)
  !$OMP THREADPRIVATE(precip_mean)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: humrel_mean                 !! Mean soil moisture stress, mean root extraction potential (unitless)
  !$OMP THREADPRIVATE(humrel_mean)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: totnobio_mean               !! Mean last total fraction of no bio (unitless;0-1)
  !$OMP THREADPRIVATE(totnobio_mean)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: vegtot_mean                 !! Mean potentially vegetated fraction (unitless;0-1)
  !$OMP THREADPRIVATE(vegtot_mean)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: k_litt_mean                 !! Mean averaged conductivity for saturated infiltration in the 'litter' layer (kg/m^2/dt)
  !$OMP THREADPRIVATE(k_litt_mean)
  !
  ! The averaged outflow fluxes.
  !
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: lakeinflow_mean              !! Mean lake inflow (kg/m^2/dt)
  !$OMP THREADPRIVATE(lakeinflow_mean)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: returnflow_mean              !! Mean water flow from lakes and swamps which returns to the grid box.
  !! This water will go back into the hydrol or hydrolc module to allow re-evaporation (kg/m^2/dt)
  !$OMP THREADPRIVATE(returnflow_mean)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: reinfiltration_mean          !! Mean water flow which returns to the grid box (kg/m^2/dt)
  !$OMP THREADPRIVATE(reinfiltration_mean)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: irrigation_mean              !! Mean irrigation flux.
  !! This is the water taken from the reservoirs and beeing put into the upper layers of the soil (kg/m^2/dt)
  !$OMP THREADPRIVATE(irrigation_mean)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: riverflow_mean               !! Mean Outflow of the major rivers.
  !! The flux will be located on the continental grid but this should be a coastal point (kg/dt)
  !$OMP THREADPRIVATE(riverflow_mean)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: coastalflow_mean             !! Mean outflow on coastal points by small basins.
  !! This is the water which flows in a disperse way into the ocean (kg/dt)
  !$OMP THREADPRIVATE(coastalflow_mean)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: floodtemp                    !! Temperature to decide if floodplains work (K)
  !$OMP THREADPRIVATE(floodtemp)
  INTEGER(i_std), SAVE                                       :: floodtemp_lev                !! Temperature level to decide if floodplains work (K)
  !$OMP THREADPRIVATE(floodtemp_lev)
  !
  ! Diagnostic variables ... well sort of !
  !
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: irrig_netereq                !! Irrigation requirement (water requirements by the crop for its optimal growth (kg/m^2/dt)
  !$OMP THREADPRIVATE(irrig_netereq)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: hydrographs                  !! Hydrographs at the outflow of the grid box for major basins (kg/dt)
  !$OMP THREADPRIVATE(hydrographs)
  !
  ! Diagnostics for the various reservoirs we use (Kg/m^2)
  !
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: fast_diag                    !! Diagnostic for the fast reservoir (kg/m^2)
  !$OMP THREADPRIVATE(fast_diag)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: slow_diag                    !! Diagnostic for the slow reservoir (kg/m^2)
  !$OMP THREADPRIVATE(slow_diag)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: stream_diag                  !! Diagnostic for the stream reservoir (kg/m^2)
  !$OMP THREADPRIVATE(stream_diag)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: flood_diag                   !! Diagnostic for the floodplain reservoir (kg/m^2)
  !$OMP THREADPRIVATE(flood_diag)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: pond_diag                    !! Diagnostic for the pond reservoir (kg/m^2)
  !$OMP THREADPRIVATE(pond_diag)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: lake_diag                    !! Diagnostic for the lake reservoir (kg/m^2)
  !$OMP THREADPRIVATE(lake_diag)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: delsurfstor                  !! Diagnostic of the change in surface water storage (flood, pond and lake reservoirs) (kg/m^2)
  !$OMP THREADPRIVATE(delsurfstor)
  
  !
  ! Specific variables for interp_topo_v0 routing
  !
  REAL(r_std),SAVE,ALLOCATABLE :: topoind_r(:)                !! Topographic index of the retention time (m) index - (local routing grid)
  !$OMP THREADPRIVATE(topoind_r)
  INTEGER,SAVE,ALLOCATABLE     :: route_flow_rp1(:)           !! flow index from cell to neighboring cell following the trip direction - (local routing grid + halo)   
  !$OMP THREADPRIVATE(route_flow_rp1)
  REAL(r_std),SAVE,ALLOCATABLE :: fast_reservoir_r(:)         !! Water amount in the fast reservoir (kg) - (local routing grid)
  !$OMP THREADPRIVATE(fast_reservoir_r)
  REAL(r_std),SAVE,ALLOCATABLE :: slow_reservoir_r(:)         !! Water amount in the slow reservoir (kg) - (local routing grid)
  !$OMP THREADPRIVATE(slow_reservoir_r)
  REAL(r_std),SAVE,ALLOCATABLE :: stream_reservoir_r(:)       !! Water amount in the stream reservoir (kg) - (local routing grid)
  !$OMP THREADPRIVATE(stream_reservoir_r)
  LOGICAL,SAVE,ALLOCATABLE     :: is_lakeinflow_r(:)          !! is lake inflow point  - (local routing grid)
  !$OMP THREADPRIVATE(is_lakeinflow_r)
  LOGICAL,SAVE,ALLOCATABLE     :: is_coastalflow_r(:)         !! is coastal flow point - (local routing grid)
  !$OMP THREADPRIVATE(is_coastalflow_r)
  LOGICAL,SAVE,ALLOCATABLE     :: is_riverflow_r(:)           !! is river flow point - (local routing grid) 
  !$OMP THREADPRIVATE(is_riverflow_r)
  LOGICAL,SAVE,ALLOCATABLE     :: is_streamflow_r(:)          !! is stream flow point - (local routing grid) 
  !$OMP THREADPRIVATE(is_streamflow_r)
  LOGICAL,SAVE,ALLOCATABLE     :: routing_mask_r(:)           !! valid routing point - (local routing grid) 
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
  REAL(r_std),SAVE,ALLOCATABLE :: routing_weight_r(:)         !! Weight to transfer lost water into intersecting cell when doing interpolation from local routing grid to local native grid (lakeinflow)
  !$OMP THREADPRIVATE(routing_weight_r)
  REAL(r_std),SAVE,ALLOCATABLE :: routing_weight_coast_r(:)   !! Weight to transfer lost water into intersecting cell on coast 
  !$OMP THREADPRIVATE(routing_weight_coast_r)
  !! when doing interpolation from local routing grid to local native grid (river flow+coastal flow)
  REAL(r_std),SAVE,ALLOCATABLE :: basins_extended_r(:)        !! basins riverflow id (local routing grid) 
  !$OMP THREADPRIVATE(basins_extended_r)
  INTEGER(i_std)               :: basins_count                !! number of basins (local routing grid) 
  !$OMP THREADPRIVATE(basins_count)
  INTEGER(i_std)               :: basins_out                  !! number of basins to output for diag 
  !$OMP THREADPRIVATE(basins_out)

  INTEGER(i_std)               :: split_routing               !! time spliting for routing
  !$OMP THREADPRIVATE(split_routing)

  REAL(r_std), SAVE                                          :: max_lake_reservoir           !! Maximum limit of water in lake_reservoir [kg/m2]
  !$OMP THREADPRIVATE(max_lake_reservoir)


  INTEGER(i_std), PARAMETER :: nb_stations=14
  REAL(r_std),PARAMETER  :: station_lon(nb_stations) = &
       (/ -162.8830, -90.9058, -55.5110, -49.3242, -133.7447, -63.6000,  28.7167, &
       15.3000,  66.5300,  89.6700,  86.5000,  127.6500,   3.3833, 117.6200 /)
  REAL(r_std),PARAMETER  :: station_lat(nb_stations) = &
       (/  61.9340,  32.3150, -1.9470, -5.1281, 67.4583,  8.1500, 45.2167, &
       -4.3000,  66.5700, 25.1800, 67.4800, 70.7000, 11.8667, 30.7700 /)
  CHARACTER(LEN=17),PARAMETER  :: station_name(nb_stations) = &
       (/ "Pilot station    ", "Vicksburg        ", "Obidos           ", &
       "Itupiranga       ", "Arctic red river ", "Puente Angostura ", &
       "Ceatal Izmail    ", "Kinshasa         ", "Salekhard        ", &
       "Bahadurabad      ", "Igarka           ", "Kusur            ", &
       "Malanville       ", "Datong           " /)

CONTAINS

!!  =============================================================================================================================
!! SUBROUTINE:    routing_interp_topo_v0_xios_initialize
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

  SUBROUTINE routing_interp_topo_v0_xios_initialize
    USE xios
    USE routing_subgrid_halfdeg, ONLY : routing_names
    IMPLICIT NONE     

    INTEGER(i_std) ::ib

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

  END SUBROUTINE routing_interp_topo_v0_xios_initialize



!!  =============================================================================================================================
!! SUBROUTINE:    routing_interp_topo_v0_init_2
!!
!>\BRIEF	  Private subroutine
!!
!! DESCRIPTION:	  Private subroutine to the module routing_interp_topo_v0. This subroutine is called in the end of routing_interp_topo_v0_initialize
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S): None
!! 
!! FLOWCHART: None
!! \n
!_ ==============================================================================================================================

  SUBROUTINE routing_interp_topo_v0_init_2(nbpt, contfrac)
    USE xios
    USE grid, ONLY : area
    IMPLICIT NONE
    INCLUDE "mpif.h"

    !! 0 Variable and parameter description
    !! 0.1 Input variables
    INTEGER,INTENT(in) :: nbpt             !! nb points native grid
    REAL,INTENT(in)    :: contfrac(nbpt)   !! fraction of land

    !! 0.2 Local variables
    INTEGER :: ni                                         !! longitude dimension of local routing grid
    INTEGER :: nj                                         !! latitude dimension of local routing grid

    REAL(r_std),ALLOCATABLE :: trip_rp1(:,: )             !! direction of flow (1-8) or river flow (99) or coastal flow (98) or lake inflow (97) - local routing grid + halo (0:ni+1,0:nj+1)
    REAL(r_std),ALLOCATABLE :: trip_extended_r(:,:)       !! direction of flow (1-8) or river flow (99) or coastal flow (98) or lake inflow (97) - local routing grid (ni,nj)
    !! routing is artificially computed on sea and endoric basins
    LOGICAL,ALLOCATABLE     :: routing_mask_rp1(:,:)      !! valid routing point - local routing grid+halo (0:ni+1,0:nj+1) 
    REAL(r_std),ALLOCATABLE :: mask_native(:)             !! some mask on native grid (nbpt)
    REAL(r_std),ALLOCATABLE :: frac_routing_r(:)          !! fraction of routing cell intesected by cells of native grid
    REAL(r_std),ALLOCATABLE :: frac_routing_coast_r(:)    !! fraction of routing cell intesected by coastal cells of native grid
    INTEGER :: ij, ij_r, ij_rp1, i,j,jp1,jm1,ip1,im1,jr,ir
    REAL(r_std) :: sum1,sum2,sum_frac_routing_r
    REAL(r_std) :: contfrac_mpi(nbp_mpi)
    REAL(r_std) :: area_mpi(nbp_mpi)
    INTEGER :: basins_count_mpi
    INTEGER :: nb_coast_points
    INTEGER :: ierr
    LOGICAL :: file_exists

!_ ================================================================================================================================

    split_routing=1
    CALL getin_p("SPLIT_ROUTING",split_routing)

    CALL gather_omp(contfrac,contfrac_mpi)
    CALL gather_omp(area, area_mpi)
 
    IF (is_omp_root) THEN

       CALL xios_get_domain_attr("routing_domain", ni=ni, nj=nj)    ! get routing domain dimension

       nbpt_r= ni*nj                                                
       nbpt_rp1= (ni+2)*(nj+2)

       ALLOCATE(trip_extended_r(ni,nj))
       ALLOCATE(topoind_r(ni*nj))       
       ALLOCATE(is_lakeinflow_r(ni*nj))       
       ALLOCATE(is_coastalflow_r(ni*nj))       
       ALLOCATE(is_riverflow_r(ni*nj)) 
       ALLOCATE(is_streamflow_r(ni*nj)) 
       ALLOCATE(routing_mask_r(ni*nj)) 
       ALLOCATE(fast_reservoir_r(ni*nj))       
       ALLOCATE(slow_reservoir_r(ni*nj))       
       ALLOCATE(stream_reservoir_r(ni*nj))       

       ALLOCATE(mask_native(nbp_mpi))
       ALLOCATE(coast_mask(nbp_mpi))
       ALLOCATE(frac_routing_r(nbpt_r))
       ALLOCATE(frac_routing_coast_r(nbpt_r))
       ALLOCATE(routing_weight(nbp_mpi))
       ALLOCATE(routing_weight_r(nbpt_r))
       ALLOCATE(routing_weight_coast_r(nbpt_r))
       ALLOCATE(basins_extended_r(nbpt_r))

       coast_mask=.FALSE.
       WHERE( contfrac_mpi(:)< 1.-1.e-5) coast_mask(:)=.TRUE.            ! create mask for coastal cells on native grid
       
       INQUIRE(FILE="routing_start.nc", EXIST=file_exists)

       IF (file_exists) THEN  
          CALL xios_recv_field("fast_reservoir_start",fast_reservoir_r)
          CALL xios_recv_field("slow_reservoir_start",slow_reservoir_r)
          CALL xios_recv_field("stream_reservoir_start",stream_reservoir_r)
       ELSE
          fast_reservoir_r(:)=0
          slow_reservoir_r(:)=0
          stream_reservoir_r(:)=0
       ENDIF

       ALLOCATE(trip_rp1(0:ni+1,0:nj+1))
       trip_rp1(:,:)=1e10
       ALLOCATE(routing_mask_rp1(0:ni+1,0:nj+1))
       ALLOCATE(route_flow_rp1((ni+2)*(nj+2)))

       CALL xios_send_field("routing_contfrac",contfrac_mpi)       ! diag
       CALL xios_recv_field("trip_r",trip_rp1(1:ni,1:nj))          ! recv trip array with halo of 1
       CALL xios_recv_field("trip_extended_r",trip_extended_r)     ! recv extended trip array from file
       CALL xios_recv_field("topoind_r",topoind_r)                 ! recv topo index array from file
       CALL xios_recv_field("basins_extended_r",basins_extended_r) ! recv basins index from file

       is_lakeinflow_r(:) = .FALSE.
       is_coastalflow_r(:) = .FALSE.
       is_riverflow_r(:) = .FALSE.

       route_flow_rp1(:)=-1

       routing_mask_rp1=.TRUE.
       DO j=1,nj
          DO i=1,ni
             IF (trip_rp1(i,j)>99) routing_mask_rp1(i,j)=.FALSE.         ! if ocean or endoric basins => masked point
          ENDDO
       ENDDO

       mask_native(:)=0
       frac_routing_r(:)=0
       WHERE( .NOT. coast_mask(:)) mask_native(:)=1

       nb_coast_points=SUM(1-mask_native)
       CALL reduce_sum_mpi(nb_coast_points, total_coast_points)
       CALL bcast_mpi(total_coast_points)
       
       CALL xios_send_field("mask_native_lake",mask_native)              ! send full land point to XIOS (native grid)
       CALL xios_recv_field("frac_routing_lake_r",frac_routing_r)        ! receive fraction of intersected cell by full land, on routing grid 


       DO j=1,nj
          DO i=1,ni
             ij=i+ni*(j-1)
             IF (frac_routing_r(ij)>1-1e-5) THEN                                                       ! if full land point (no part of coast)
                IF ((.NOT. routing_mask_rp1(i,j)) .OR. trip_rp1(i,j)==98 .OR. trip_rp1(i,j)==99) THEN   ! if masked point or river flow or coastal flow
                   routing_mask_rp1(i,j)=.TRUE.                                                          ! transform to non masked
                   trip_rp1(i,j)=trip_extended_r(i,j)                                                  ! replace by extended trip value         
                ENDIF
             ENDIF
          ENDDO
       ENDDO

       mask_native(:)=0
       frac_routing_coast_r(:)=0
       WHERE( coast_mask(:)) mask_native(:)=1
       CALL xios_send_field("mask_native_coast",mask_native)               ! send only coastal point to XIOS (native grid)
       CALL xios_recv_field("frac_routing_coast_r",frac_routing_coast_r)   ! receive fraction of intersected cell by coastal point, on routing grid 

       DO j=1,nj
          DO i=1,ni
             ij=i+ni*(j-1)
             IF (frac_routing_coast_r(ij) > 0) THEN                        ! fraction of coastal point ?
                IF (.NOT. routing_mask_rp1(i,j)) THEN                      ! if masked, transform to coastal flow to conserve the water balance
                   routing_mask_rp1(i,j)=.TRUE.
                   trip_rp1(i,j)=98     ! transform to coastal flow
                   topoind_r(ij)=1e-10  ! no residence time
                ENDIF
             ENDIF
          ENDDO
       ENDDO


       mask_native(:)=1
       frac_routing_r(:)=0
       CALL xios_send_field("mask_native",mask_native)              ! send only all point to XIOS (native grid)
       CALL xios_recv_field("frac_routing_r",frac_routing_r)        ! receive fraction of intersected cell, (routing grid)
       DO j=1,nj
          DO i=1,ni
             ij=i+ni*(j-1)
             IF (frac_routing_r(ij)>0) THEN
                IF (.NOT. routing_mask_rp1(i,j)) THEN
                   routing_mask_rp1(i,j)=.TRUE.                         ! may never happen
                ENDIF
             ELSE
                routing_mask_rp1(i,j)=.FALSE.                           ! if no intersection from native cells, mask the point (no routage) 
                trip_rp1(i,j)=1e10
             ENDIF
          ENDDO
       ENDDO

       CALL xios_send_field("trip_update_r",trip_rp1(1:ni,1:nj))        ! send to xios trip array to update for halo
       CALL xios_recv_field("trip_rp1",trip_rp1)                        ! recv trip array with halo of 1


       routing_mask_rp1=.TRUE.
       DO j=0,nj+1
          DO i=0,ni+1
             IF (trip_rp1(i,j)>99) routing_mask_rp1(i,j)=.FALSE.        !  update mask with halo of 1
          ENDDO
       ENDDO

       !! Compute the routing
       !! Loop on all point point of the local routing grid + halo
       DO j=0,nj+1                                                 
          jp1=j+1
          jm1=j-1
          DO i=0,ni+1
             ij_rp1=i+(ni+2)*j+1
             ij_r=i+ni*(j-1)
             ip1=i+1
             im1=i-1

             IF (trip_rp1(i,j) < 100) THEN                    

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

                   IF ( .NOT. routing_mask_rp1(ir,jr)) THEN                      ! if routing to a masked point, cut the flow
                      jr=j ; ir=i ;                                              ! and route to myself
                      IF (i>0 .AND. i<ni+1 .AND. j>0 .AND. j<nj+1) THEN          ! if some flow come from neighbour cell 
                         IF (     trip_rp1(i-1,j-1)==4 .OR. trip_rp1(i,j-1)==5 .OR. trip_rp1(i+1,j)==6     &
                              .OR. trip_rp1(i-1,j)==3   .OR. trip_rp1(i+1,j)==7                             & 
                              .OR. trip_rp1(i-1,j+1)==2 .OR. trip_rp1(i-1,j+1)==1 .OR. trip_rp1(i+1,j+1)==8) THEN
                            is_riverflow_r(ij_r)=.TRUE.                          !  => transform to riverflow
                         ELSE
                            is_coastalflow_r(ij_r)=.TRUE.                        !  else transform to coastalflow
                         ENDIF
                      ENDIF
                   ENDIF

                   IF (ir<1 .OR. ir>ni .OR. jr<1 .OR. jr>nj) THEN  
                      route_flow_rp1(ij_rp1)=-1                                  ! if route outside my local domain, no routing (will be done by other process)
                   ELSE
                      route_flow_rp1(ij_rp1)=ir+ni*(jr-1)                        ! define the cell where to flow
                   ENDIF
                ENDIF

             ENDIF
          ENDDO
       ENDDO

       routing_mask_r(:)=reshape(routing_mask_rp1(1:ni,1:nj),(/ni*nj/))

       is_streamflow_r(:)=.NOT. (is_lakeinflow_r(:) .OR. is_coastalflow_r(:) .OR. is_riverflow_r(:)) .AND. routing_mask_r(:)

       DO ij=1,nbp_mpi
          routing_weight(ij)=contfrac_mpi(ij)*area_mpi(ij)  
       ENDDO

       routing_weight_r(:)=0
       DO ij=1,nbpt_r
          IF (frac_routing_r(ij)>0) routing_weight_r(ij)=1./frac_routing_r(ij)
       ENDDO

       routing_weight_coast_r(:)=0.
       DO ij=1,nbpt_r
          IF (frac_routing_coast_r(ij)>0) routing_weight_coast_r(ij)=1./frac_routing_coast_r(ij)
       ENDDO

       CALL xios_send_field("routing_weight_coast_r",routing_weight_coast_r)
       ! looking for basins
       basins_count_mpi=0
       DO ij=1,nbpt_r
          IF (basins_extended_r(ij) > basins_count_mpi) basins_count_mpi=basins_extended_r(ij)
       ENDDO
       CALL MPI_ALLREDUCE(basins_count_mpi,basins_count,1,MPI_INT_ORCH, MPI_MAX,MPI_COMM_ORCH,ierr)

    ELSE

       nbpt_r=0

    ENDIF

  END SUBROUTINE routing_interp_topo_v0_init_2


  !! ================================================================================================================================
  !! SUBROUTINE 	: routing_interp_topo_v0_flow
  !!
  !>\BRIEF         This subroutine computes the transport of water in the various reservoirs
  !!                (including ponds and floodplains) and the water withdrawals from the reservoirs for irrigation.
  !!
  !! DESCRIPTION (definitions, functional, design, flags) :
  !! This will first compute the amount of water which flows out of each of the 3 reservoirs using the assumption of an 
  !! exponential decrease of water in the reservoir (see Hagemann S and Dumenil L. (1998)). Then we compute the fluxes 
  !! for floodplains and ponds. All this will then be used in order to update each of the basins : taking water out of 
  !! the up-stream basin and adding it to the down-stream one.
  !! As this step happens globaly we have to stop the parallel processing in order to exchange the information. Once 
  !! all reservoirs are updated we deal with irrigation. The final step is to compute diagnostic fluxes. Among them
  !! the hydrographs of the largest rivers we have chosen to monitor.
  !!
  !! RECENT CHANGE(S): None
  !!
  !! MAIN OUTPUT VARIABLE(S): lakeinflow, returnflow, reinfiltration, irrigation, riverflow, coastalflow, hydrographs, flood_frac, flood_res
  !!
  !! REFERENCES   :
  !! - Ngo-Duc, T., K. Laval, G. Ramillien, J. Polcher, and A. Cazenave (2007)
  !!   Validation of the land water storage simulated by Organising Carbon and Hydrology in Dynamic Ecosystems (ORCHIDEE) with Gravity Recovery and Climate Experiment (GRACE) data.
  !!   Water Resour. Res., 43, W04427, doi:10.1029/2006WR004941.
  !! * Irrigation:
  !! - de Rosnay, P., J. Polcher, K. Laval, and M. Sabre (2003)
  !!   Integrated parameterization of irrigation in the land surface model ORCHIDEE. Validation over Indian Peninsula.
  !!   Geophys. Res. Lett., 30(19), 1986, doi:10.1029/2003GL018024.
  !! - A.C. Vivant (2003)
  !!   Les plaines d'inondations et l'irrigation dans ORCHIDEE, impacts de leur prise en compte.
  !!   , , 51pp.
  !! - N. Culson (2004)
  !!   Impact de l'irrigation sur le cycle de l'eau
  !!   Master thesis, Paris VI University, 55pp.
  !! - X.-T. Nguyen-Vinh (2005)
  !!   Analyse de l'impact de l'irrigation en Amerique du Nord - plaine du Mississippi - sur la climatologie regionale
  !!   Master thesis, Paris VI University, 33pp.
  !! - M. Guimberteau (2006)
  !!   Analyse et modifications proposees de la modelisation de l'irrigation dans un modele de surface.
  !!   Master thesis, Paris VI University, 46pp.
  !! - Guimberteau M. (2010)
  !!   Modelisation de l'hydrologie continentale et influences de l'irrigation sur le cycle de l'eau.
  !!   Ph.D. thesis, Paris VI University, 195pp.
  !! - Guimberteau M., Laval K., Perrier A. and Polcher J. (2011).
  !!   Global effect of irrigation and its impact on the onset of the Indian summer monsoon.
  !!   In press, Climate Dynamics, doi: 10.1007/s00382-011-1252-5.
  !! * Floodplains:
  !! - A.C. Vivant (2002)
  !!   L'ecoulement lateral de l'eau sur les surfaces continentales. Prise en compte des plaines d'inondations dans ORCHIDEE.
  !!   Master thesis, Paris VI University, 46pp.
  !! - A.C. Vivant (2003)
  !!   Les plaines d'inondations et l'irrigation dans ORCHIDEE, impacts de leur prise en compte.
  !!   , , 51pp.
  !! - T. d'Orgeval (2006)
  !!   Impact du changement climatique sur le cycle de l'eau en Afrique de l'Ouest: modelisation et incertitudes.
  !!   Ph.D. thesis, Paris VI University, 188pp.
  !! - T. d'Orgeval, J. Polcher, and P. de Rosnay (2008)
  !!   Sensitivity of the West African hydrological cycle in ORCHIDEE to infiltration processes.
  !!   Hydrol. Earth Syst. Sci., 12, 1387-1401
  !! - M. Guimberteau, G. Drapeau, J. Ronchail, B. Sultan, J. Polcher, J.-M. Martinez, C. Prigent, J.-L. Guyot, G. Cochonneau,
  !!   J. C. Espinoza, N. Filizola, P. Fraizy, W. Lavado, E. De Oliveira, R. Pombosa, L. Noriega, and P. Vauchel (2011)
  !!   Discharge simulation in the sub-basins of the Amazon using ORCHIDEE forced by new datasets.
  !!   Hydrol. Earth Syst. Sci. Discuss., 8, 11171-11232, doi:10.5194/hessd-8-11171-2011
  !!
  !! FLOWCHART    :None
  !! \n
  !_ ================================================================================================================================

  SUBROUTINE routing_interp_topo_v0_flow(nbpt, dt_routing, lalo, floodout, runoff_omp, drainage_omp, &
                                         vegtot, totnobio, transpot_mean, precip, humrel, k_litt, floodtemp, reinf_slope, &
                                         lakeinflow_omp, returnflow, reinfiltration, irrigation, riverflow_omp, &
                                         coastalflow_omp, hydrographs, slowflow_diag, flood_frac, flood_res)
    
    USE xios
    USE grid, ONLY : area 
    IMPLICIT NONE
    INCLUDE "mpif.h"
    
    !! 0 Variable and parameter description
    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                   :: nbpt                      !! Domain size (unitless)
    REAL(r_std), INTENT (in)                     :: dt_routing                !! Routing time step (s)
    REAL(r_std), INTENT(in)                      :: lalo(nbpt,2)              !! Vector of latitude and longitudes
    REAL(r_std), INTENT(in)                      :: runoff_omp(nbpt)              !! Grid-point runoff (kg/m^2/dt)
    REAL(r_std), INTENT(in)                      :: floodout(nbpt)            !! Grid-point flow out of floodplains (kg/m^2/dt)
    REAL(r_std), INTENT(in)                      :: drainage_omp(nbpt)            !! Grid-point drainage (kg/m^2/dt)
    REAL(r_std), INTENT(in)                      :: vegtot(nbpt)              !! Potentially vegetated fraction (unitless;0-1)
    REAL(r_std), INTENT(in)                      :: totnobio(nbpt)            !! Other areas which can not have vegetation
    REAL(r_std), INTENT(in)                      :: transpot_mean(nbpt)       !! Mean potential transpiration of the vegetation (kg/m^2/dt)
    REAL(r_std), INTENT(in)                      :: precip(nbpt)              !! Rainfall (kg/m^2/dt)
    REAL(r_std), INTENT(in)                      :: humrel(nbpt)              !! Soil moisture stress, root extraction potential (unitless)
    REAL(r_std), INTENT(in)                      :: k_litt(nbpt)              !! Averaged conductivity for saturated infiltration in the 'litter' layer (kg/m^2/dt)
    REAL(r_std), INTENT(in)                      :: floodtemp(nbpt)           !! Temperature to decide if floodplains work (K)
    REAL(r_std), INTENT(in)                      :: reinf_slope(nbpt)         !! Coefficient which determines the reinfiltration ratio in the grid box due to flat areas (unitless;0-1)

    !! 0.2 Output variables

    REAL(r_std), INTENT(out)                     :: lakeinflow_omp(nbpt)          !! Water inflow to the lakes (kg/dt)
    REAL(r_std), INTENT(out)                     :: returnflow(nbpt)          !! The water flow from lakes and swamps which returns into the grid box.
    !! This water will go back into the hydrol or hydrolc module to allow re-evaporation (kg/m^2/dt)
    REAL(r_std), INTENT(out)                     :: reinfiltration(nbpt)      !! Water flow from ponds and floodplains which returns to the grid box (kg/m^2/dt)
    REAL(r_std), INTENT(out)                     :: irrigation(nbpt)          !! Irrigation flux. This is the water taken from the reservoirs and beeing put into the upper layers of the soil (kg/m^2/dt)
    REAL(r_std), INTENT(out)                     :: riverflow_omp(nbpt)           !! Outflow of the major rivers. The flux will be located on the continental grid but this should be a coastal point (kg/dt)
    REAL(r_std), INTENT(out)                     :: coastalflow_omp(nbpt)         !! Outflow on coastal points by small basins. This is the water which flows in a disperse way into the ocean (kg/dt)
    REAL(r_std), INTENT(out)                     :: hydrographs(nbpt)         !! Hydrographs at the outflow of the grid box for major basins (kg/dt)
    REAL(r_std), INTENT(out)                     :: slowflow_diag(nbpt)       !! Hydrographs of slow_flow = routed slow_flow for major basins (kg/dt)
    REAL(r_std), INTENT(out)                     :: flood_frac(nbpt)          !! Flooded fraction of the grid box (unitless;0-1)
    REAL(r_std), INTENT(out)                     :: flood_res(nbpt)           !! Diagnostic of water amount in the floodplains reservoir (kg)

    !! 0.4 Local variables
    REAL(r_std)                                  :: runoff(nbp_mpi)              !! Grid-point runoff (kg/dt)
    REAL(r_std)                                  :: drainage(nbp_mpi)            !! Grid-point drainage (kg/dt)
    REAL(r_std)                                  :: riverflow(nbp_mpi)             
    REAL(r_std)                                  :: coastalflow(nbp_mpi)           
    REAL(r_std)                                  :: lakeinflow(nbp_mpi)            
    REAL(r_std)                                  :: fast_diag_mpi(nbp_mpi)          
    REAL(r_std)                                  :: slow_diag_mpi(nbp_mpi)        
    REAL(r_std)                                  :: stream_diag_mpi(nbp_mpi)         
    REAL(r_std)                                  :: area_mpi(nbp_mpi) ! cell area         
    REAL(r_std)                                  :: lake_reservoir_mpi(nbp_mpi) ! cell area         

    ! from input model -> routing_grid
    REAL(r_std)                      :: runoff_r(nbpt_r)              !! Grid-point runoff (kg/m^2/dt)
    REAL(r_std)                      :: floodout_r(nbpt_r)            !! Grid-point flow out of floodplains (kg/m^2/dt)
    REAL(r_std)                      :: drainage_r(nbpt_r)            !! Grid-point drainage (kg/m^2/dt)
    REAL(r_std)                      :: vegtot_r(nbpt_r)              !! Potentially vegetated fraction (unitless;0-1)
    REAL(r_std)                      :: totnobio_r(nbpt_r)            !! Other areas which can not have vegetation
    REAL(r_std)                      :: transpot_mean_r(nbpt_r)       !! Mean potential transpiration of the vegetation (kg/m^2/dt)
    REAL(r_std)                      :: precip_r(nbpt_r)              !! Rainfall (kg/m^2/dt)
    REAL(r_std)                      :: humrel_r(nbpt_r)              !! Soil moisture stress, root extraction potential (unitless)
    REAL(r_std)                      :: k_litt_r(nbpt_r)              !! Averaged conductivity for saturated infiltration in the 'litter' layer (kg/m^2/dt)
    REAL(r_std)                      :: floodtemp_r(nbpt_r)           !! Temperature to decide if floodplains work (K)
    REAL(r_std)                      :: reinf_slope_r(nbpt_r)         !! Coefficient which determines the reinfiltration ratio in the grid box due to flat areas (unitless;0-1)



    REAL(r_std), DIMENSION(nbpt_r)        :: fast_flow_r                 !! Outflow from the fast reservoir (kg/dt)
    REAL(r_std), DIMENSION(nbpt_r)        :: slow_flow_r                 !! Outflow from the slow reservoir (kg/dt)
    REAL(r_std), DIMENSION(nbpt_r)        :: stream_flow_r               !! Outflow from the stream reservoir (kg/dt)
    REAL(r_std), DIMENSION(nbpt_r)        :: hydrographs_r                !! hydrograph (kg/dt)
    REAL(r_std), DIMENSION(nbpt_r)        :: flood_flow_r                !! Outflow from the floodplain reservoir (kg/dt)
    REAL(r_std), DIMENSION(nbpt_r)        :: pond_inflow_r               !! Inflow to the pond reservoir (kg/dt)
    REAL(r_std), DIMENSION(nbpt_r)        :: pond_drainage_r             !! Drainage from pond (kg/m^2/dt)
    REAL(r_std), DIMENSION(nbpt_r)        :: flood_drainage_r            !! Drainage from floodplains (kg/m^2/dt)
    REAL(r_std), DIMENSION(nbpt_r)        :: return_swamp_r              !! Inflow to the swamp (kg/dt)
    !
    ! Irrigation per basin
    !
    REAL(r_std), DIMENSION(nbpt_r)        :: irrig_needs_r               !! Total irrigation requirement (water requirements by the crop for its optimal growth) (kg)
    REAL(r_std), DIMENSION(nbpt_r)        :: irrig_actual_r              !! Possible irrigation according to the water availability in the reservoirs (kg)
    REAL(r_std), DIMENSION(nbpt_r)        :: irrig_deficit_r             !! Amount of water missing for irrigation (kg)
    REAL(r_std), DIMENSION(nbpt_r)        :: irrig_adduct_r              !! Amount of water carried over from other basins for irrigation (kg)
    !
    REAL(r_std), DIMENSION(nbpt_r)        :: transport_r                 !! Water transport between basins (kg/dt)
    REAL(r_std), DIMENSION(nbpt_r)        :: floods_r                    !! Water flow in to the floodplains (kg/dt)
    REAL(r_std), DIMENSION(nbpt_r)        :: potflood_r                  !! Potential inflow to the swamps (kg/dt)
    REAL(r_std), DIMENSION(nbpt)                 :: tobeflooded               !! Maximal surface which can be inundated in each grid box (m^2)
    REAL(r_std), DIMENSION(nbpt)                 :: totarea                   !! Total area of basin (m^2)
    REAL(r_std), DIMENSION(nbpt)                 :: totflood                  !! Total amount of water in the floodplains reservoir (kg) (not used)
    REAL(r_std), DIMENSION(nbasmax)              :: pond_excessflow           !! 
    REAL(r_std)                                  :: flow                      !! Outflow computation for the reservoirs (kg/dt)
    REAL(r_std)                                  :: floodindex                !! Fraction of grid box area inundated (unitless;0-1)
    REAL(r_std)                                  :: pondex                    !! 
    REAL(r_std)                                  :: flood_frac_pot            !! Total fraction of the grid box which is flooded at optimum repartition (unitless;0-1)
    REAL(r_std)                                  :: stream_tot                !! Total water amount in the stream reservoirs (kg)
    REAL(r_std)                                  :: adduction                 !! Importation of water from a stream reservoir of a neighboring grid box (kg)
    REAL(r_std), DIMENSION(8,nbasmax)            :: streams_around            !! Stream reservoirs of the neighboring grid boxes (kg)
    INTEGER(i_std), DIMENSION(8)                 :: igrd                      !! 
    INTEGER(i_std), DIMENSION(2)                 :: ff                        !! 
    INTEGER(i_std), DIMENSION(1)                 :: fi                        !! 
    INTEGER(i_std)                               :: ig, ib, ib2, ig2          !! Indices (unitless)
    INTEGER(i_std)                               :: rtg, rtb, in              !! Indices (unitless)
    INTEGER(i_std)                               :: ier                       !! Error handling 
    INTEGER(i_std)                               :: isplit
    
    LOGICAL, PARAMETER                           :: check_reservoir = .TRUE. !! Logical to choose if we write informations when a negative amount of water is occurring in a reservoir (true/false)
    REAL(r_std), DIMENSION(nbpt_r)        :: flood_frac_bas_r   
    REAL(r_std), DIMENSION(nbpt_r)        :: flood_reservoir_r   
    REAL(r_std), DIMENSION(nbpt_r)        :: pond_reservoir_r   

    REAL(r_std), DIMENSION(nbpt_r)        :: lakeinflow_r   
    REAL(r_std), DIMENSION(nbpt_r)        :: coastalflow_r   
    REAL(r_std), DIMENSION(nbpt_r)        :: riverflow_r   

    REAL(r_std), DIMENSION(nbpt_r)        :: flow_r   
    REAL(r_std), DIMENSION(nbpt_rp1)      :: flow_rp1   
    REAL(r_std) :: water_balance_before, water_balance_after
    REAL(r_std) :: basins_riverflow_mpi(0:basins_count)
    REAL(r_std) :: basins_riverflow(0:basins_count)
    REAL(r_std) :: lake_overflow,sum_lake_overflow, total_lake_overflow
    INTEGER :: ierr

    !_ ================================================================================================================================
    !


    irrig_netereq(:) = zero
    totarea(:) = zero
    totflood(:) = zero

    flood_frac_bas_r(:)=zero
    !
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
    CALL gather_omp(runoff_omp,runoff)
    CALL gather_omp(drainage_omp, drainage)
    CALL gather_omp(area, area_mpi)
    CALL gather_omp(lake_reservoir, lake_reservoir_mpi)

    IF (is_omp_root) THEN

       hydrographs_r(:)=0

       DO isplit=1,split_routing


          irrig_needs_r(:) = zero
          irrig_actual_r(:) = zero
          irrig_deficit_r(:) = zero
          irrig_adduct_r(:) = zero

          DO ig=1,nbpt_r
             IF ( routing_mask_r(ig) ) THEN
                !
                ! Each of the fluxes is limited by the water in the reservoir and a small margin 
                ! (min_reservoir) to avoid rounding errors.
                !
                !               IF (fast_reservoir_r(ig)>1e-30) THEN
                !                 PRINT*,"fast_reservoir > 0  ",ig, fast_reservoir_r(ig), topoind_r(ig)
                !               ENDIF 
                flow = MIN(fast_reservoir_r(ig)/((topoind_r(ig)/1000.)*fast_tcst*one_day/(dt_routing/split_routing)),&
                     & fast_reservoir_r(ig)-min_sechiba)
                fast_flow_r(ig) = MAX(flow, zero)
                !
                flow = MIN(slow_reservoir_r(ig)/((topoind_r(ig)/1000.)*slow_tcst*one_day/(dt_routing/split_routing)),&
                     & slow_reservoir_r(ig)-min_sechiba)
                slow_flow_r(ig) = MAX(flow, zero)
                !
                flow = MIN(stream_reservoir_r(ig)/((topoind_r(ig)/1000.)*stream_tcst* & 
                     & MAX(un-SQRT(flood_frac_bas_r(ig)),min_sechiba)*one_day/(dt_routing/split_routing)),&
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

          DO ig=1,nbpt_r
             flood_drainage_r(ig) = zero
             flood_flow_r(ig) = zero
             flood_reservoir_r(ig) = zero
          ENDDO

          DO ig=1,nbpt_r
             pond_inflow_r(ig) = zero
             pond_drainage_r(ig) = zero
             pond_reservoir_r(ig) = zero
          ENDDO


          !-
          !- Compute the transport
          !-

          flow_r(:)=fast_flow_r(:) + slow_flow_r(:) + stream_flow_r(:)

          CALL xios_send_field("flow_r",flow_r)       ! transfer halo
          CALL xios_recv_field("flow_rp1",flow_rp1)

          transport_r(:)=0

          DO ig=1,nbpt_rp1
             IF ( route_flow_rp1(ig) > 0 ) THEN
                transport_r(route_flow_rp1(ig))=transport_r(route_flow_rp1(ig))+ flow_rp1(ig)
             ENDIF
          ENDDO


          !-
          !- Do the floodings - First initialize
          !-
          return_swamp_r(:)=zero
          floods_r(:)=zero

          !
          ! Update all reservoirs
          !> The slow and deep reservoir (slow_reservoir) collect the deep drainage whereas the
          !> fast_reservoir collects the computed surface runoff. Both discharge into a third reservoir
          !> (stream_reservoir) of the next sub-basin downstream.
          !> Water from the floodplains reservoir (flood_reservoir) flows also into the stream_reservoir of the next sub-basin downstream.
          !> Water that flows into the pond_reservoir is withdrawn from the fast_reservoir.
          !
          runoff=runoff*routing_weight
          CALL xios_send_field("routing_runoff",runoff)   ! interp conservative model -> routing
          CALL xios_recv_field("routing_runoff_r",runoff_r)
          drainage=drainage*routing_weight
          CALL xios_send_field("routing_drainage",drainage)   ! interp conservative model -> routing
          CALL xios_recv_field("routing_drainage_r",drainage_r)

          CALL xios_send_field("fast_reservoir_r",fast_reservoir_r)
          CALL xios_send_field("slow_reservoir_r",slow_reservoir_r)
          CALL xios_send_field("stream_reservoir_r",stream_reservoir_r)

          WHERE (.NOT. routing_mask_r) runoff_r(:)=0
          WHERE (.NOT. routing_mask_r) drainage_r(:)=0

          CALL MPI_ALLREDUCE(sum(runoff+drainage),water_balance_before,1,MPI_REAL_ORCH,MPI_SUM,MPI_COMM_ORCH,ierr)
          CALL MPI_ALLREDUCE(sum(runoff_r+drainage_r),water_balance_after,1,MPI_REAL_ORCH,MPI_SUM,MPI_COMM_ORCH,ierr)

          PRINT *,"routing water Balance ;  before : ", water_balance_before," ; after : ",water_balance_after,  &
               " ; delta : ", 100.*(water_balance_after-water_balance_before)/(0.5*(water_balance_after+water_balance_before)),"%"

          CALL xios_recv_field("water_balance_before",water_balance_before)


          DO ig=1,nbpt_r
             IF ( routing_mask_r(ig) ) THEN

                !
                fast_reservoir_r(ig) =  fast_reservoir_r(ig) + runoff_r(ig) - &
                     & fast_flow_r(ig) - pond_inflow_r(ig)
                !
                slow_reservoir_r(ig) = slow_reservoir_r(ig) + drainage_r(ig) - &
                     & slow_flow_r(ig)
                !

                stream_reservoir_r(ig) = stream_reservoir_r(ig) + flood_flow_r(ig)  - &
                     & stream_flow_r(ig) - return_swamp_r(ig) - floods_r(ig)

                IF (is_streamflow_r(ig)) stream_reservoir_r(ig)= stream_reservoir_r(ig) +  transport_r(ig)  
                !
                flood_reservoir_r(ig) = flood_reservoir_r(ig) + floods_r(ig) - &
                     & flood_flow_r(ig) 
                !
                pond_reservoir_r(ig) = pond_reservoir_r(ig) + pond_inflow_r(ig) - pond_drainage_r(ig)
                !
                IF ( flood_reservoir_r(ig) .LT. zero ) THEN
                   IF ( check_reservoir ) THEN
                      WRITE(numout,*) "WARNING : negative flood reservoir at :", ig, ". Problem is being corrected."
                      WRITE(numout,*) "flood_reservoir, floods, flood_flow : ", flood_reservoir_r(ig), floods_r(ig), &
                           & flood_flow_r(ig) 
                   ENDIF
                   stream_reservoir_r(ig) = stream_reservoir_r(ib) + flood_reservoir_r(ib)
                   flood_reservoir_r(ib) = zero
                ENDIF
                !
                IF ( stream_reservoir_r(ig) .LT. zero ) THEN
                   IF ( check_reservoir ) THEN
                      WRITE(numout,*) "WARNING : negative stream reservoir at :", ig, ". Problem is being corrected."
                      WRITE(numout,*) "stream_reservoir, flood_flow, transport : ", stream_reservoir_r(ig), flood_flow_r(ig), &
                           &  transport_r(ig)
                      WRITE(numout,*) "stream_flow, return_swamp, floods :", stream_flow_r(ig), return_swamp_r(ig), floods_r(ig)
                   ENDIF
                   fast_reservoir_r(ig) =  fast_reservoir_r(ig) + stream_reservoir_r(ig)
                   stream_reservoir_r(ig) = zero
                ENDIF
                !
                IF ( fast_reservoir_r(ig) .LT. zero ) THEN
                   IF ( check_reservoir ) THEN
                      WRITE(numout,*) "WARNING : negative fast reservoir at :", ig, ". Problem is being corrected."
                      WRITE(numout,*) "fast_reservoir, runoff, fast_flow, ponf_inflow  : ", fast_reservoir_r(ig), &
                           &runoff_r(ig), fast_flow_r(ig), pond_inflow_r(ig)
                   ENDIF
                   slow_reservoir_r(ig) =  slow_reservoir_r(ig) + fast_reservoir_r(ig)
                   fast_reservoir_r(ig) = zero
                ENDIF

                IF ( slow_reservoir_r(ig) .LT. - min_sechiba ) THEN
                   WRITE(numout,*) 'WARNING : There is a negative reservoir at :', ig
                   WRITE(numout,*) 'WARNING : slowr, slow_flow, drainage', &
                        & slow_reservoir_r(ig), slow_flow_r(ig), drainage_r(ig)
                   WRITE(numout,*) 'WARNING : pondr, pond_inflow, pond_drainage', &
                        & pond_reservoir_r(ig), pond_inflow_r(ig), pond_drainage_r(ig)
                   CALL ipslerr_p(2, 'routing_interp_topo_v0_flow', 'WARNING negative slow_reservoir.','','')
                ENDIF

             ENDIF
          ENDDO

          DO ig=1,nbpt_r
             IF ( routing_mask_r(ig) ) THEN
                hydrographs_r(ig)=hydrographs_r(ig)+fast_flow_r(ig)+slow_flow_r(ig)+stream_flow_r(ig)
             ENDIF
          ENDDO

       ENDDO ! isplit

       !
       !
       !
       ! The three output types for the routing : endorheic basins, rivers and 
       ! diffuse coastal flow.
       !
       lakeinflow_r(:)=0
       coastalflow_r(:)=0
       riverflow_r(:)=0
       basins_riverflow_mpi(:)=0

       DO ig=1,nbpt_r
          IF ( routing_mask_r(ig) ) THEN

             IF (is_lakeinflow_r(ig))  THEN
                lakeinflow_r(ig) = transport_r(ig)
                basins_riverflow_mpi(basins_extended_r(ig)) = &
                     basins_riverflow_mpi(basins_extended_r(ig))+lakeinflow_r(ig)
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


       CALL MPI_ALLREDUCE(basins_riverflow_mpi,basins_riverflow,basins_count+1,MPI_REAL_ORCH,MPI_SUM,MPI_COMM_ORCH,ierr)
       CALL xios_send_field("basins_riverflow",basins_riverflow(1:basins_out)/1000./dt_routing)

       !
       flood_res = flood_diag + pond_diag
       !
       CALL xios_send_field("routing_lakeinflow_r"  ,lakeinflow_r*routing_weight_r)
       CALL xios_recv_field("routing_lakeinflow"    ,lakeinflow)
       CALL xios_send_field("routing_coastalflow_r" ,coastalflow_r*routing_weight_coast_r)
       CALL xios_recv_field("routing_coastalflow_temp"   ,coastalflow)
       WHERE(.NOT. coast_mask(:)) coastalflow(:)=0
       CALL xios_send_field("routing_coastalflow"   ,coastalflow)

       CALL xios_send_field("routing_riverflow_r"   ,riverflow_r*routing_weight_coast_r)
       CALL xios_recv_field("routing_riverflow_temp"     ,riverflow)
       WHERE(.NOT. coast_mask(:)) riverflow(:)=0
       CALL xios_send_field("routing_riverflow"     ,riverflow)

       CALL xios_send_field("out_flow",lakeinflow+coastalflow+riverflow)
       CALL xios_send_field("routing_hydrographs_r", hydrographs_r/1000./dt_routing)

       ! diag
       CALL xios_send_field("routing_fast_reservoir_r"  , fast_reservoir_r)
       CALL xios_recv_field("routing_fast_reservoir"  , fast_diag_mpi)

       CALL xios_send_field("routing_slow_reservoir_r"  , slow_reservoir_r)
       CALL xios_recv_field("routing_slow_reservoir"  ,   slow_diag_mpi)

       CALL xios_send_field("routing_stream_reservoir_r"  , stream_reservoir_r)
       CALL xios_recv_field("routing_stream_reservoir"    , stream_diag_mpi)

       CALL xios_recv_field("water_balance_after"  ,water_balance_after)

       PRINT *,"routing water Balance ;  before : ", water_balance_before," ; after : ",water_balance_after,  &
            " ; delta : ", 100*(water_balance_after-water_balance_before)/(0.5*(water_balance_after+water_balance_before)),"%"

   
   !! Remove water from lake reservoir if it exceeds the maximum limit and distribute it 
    !! uniformly over all possible the coastflow gridcells
    
    ! Calculate lake_overflow and remove it from lake_reservoir
      sum_lake_overflow=0
      DO ig=1,nbp_mpi
         lake_overflow = MAX(0., lake_reservoir_mpi(ig) - max_lake_reservoir*area_mpi(ig))
         lake_reservoir_mpi(ig) = lake_reservoir_mpi(ig) - lake_overflow
         sum_lake_overflow = sum_lake_overflow+lake_overflow
      END DO

    ! Calculate the sum of the lake_overflow and distribute it uniformly over all gridboxes
      CALL reduce_sum_mpi(sum_lake_overflow,total_lake_overflow)
      CALL bcast_mpi(total_lake_overflow)

      WHERE(coast_mask) coastalflow = coastalflow + total_lake_overflow/total_coast_points

    ENDIF ! is_omp_root

    CALL scatter_omp(riverflow,riverflow_omp)
    CALL scatter_omp(coastalflow,coastalflow_omp)
    CALL scatter_omp(lakeinflow,lakeinflow_omp)
    CALL scatter_omp(lake_reservoir_mpi,lake_reservoir)
    CALL scatter_omp(fast_diag_mpi,fast_diag)
    CALL scatter_omp(slow_diag_mpi,slow_diag)
    CALL scatter_omp(stream_diag_mpi,stream_diag)

    DO ig=1,nbpt
       fast_diag(ig)=fast_diag(ig)/area(ig)  
       slow_diag(ig)=slow_diag(ig)/area(ig)    
       stream_diag(ig)=stream_diag(ig)/area(ig)  
    ENDDO


    flood_frac(:) = zero
    flood_height(:) = zero
    flood_frac_bas_r(:) = zero

    returnflow(:) = zero
    reinfiltration(:) = zero

    !
    !
    ! Compute the fluxes which leave the routing scheme
    !
    ! Lakeinflow is in Kg/dt
    ! returnflow is in Kg/m^2/dt
    !
    delsurfstor(:) = -flood_diag(:)-pond_diag(:)-lake_diag(:)
    hydrographs(:) = zero
    slowflow_diag(:) = zero
    flood_diag(:) =  zero
    pond_diag(:) =  zero
    irrigation(:) = zero



  END SUBROUTINE routing_interp_topo_v0_flow


  !!  =============================================================================================================================
  !! SUBROUTINE:         routing_interp_topo_v0_initialize
  !!
  !>\BRIEF	         Initialize the routing_interp_topo_v0 module
  !!
  !! DESCRIPTION:        Initialize the routing_interp_topo_v0 module. Read from restart file or read the routing.nc file
  !!                     to initialize the routing scheme. 
  !!
  !! RECENT CHANGE(S)
  !!
  !! REFERENCE(S)
  !! 
  !! FLOWCHART   
  !! \n
  !_ ==============================================================================================================================

  SUBROUTINE routing_interp_topo_v0_initialize( kjit,       nbpt,           index,                 &
                                                rest_id,     hist_id,        hist2_id,   lalo,      &
                                                neighbours,  resolution,     contfrac,   stempdiag, &
                                                returnflow,  reinfiltration, irrigation, riverflow, &
                                                coastalflow, flood_frac,     flood_res )

    IMPLICIT NONE

    !! 0 Variable and parameter description
    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)     :: kjit                 !! Time step number (unitless)
    INTEGER(i_std), INTENT(in)     :: nbpt                 !! Domain size (unitless)
    INTEGER(i_std), INTENT(in)     :: index(nbpt)          !! Indices of the points on the map (unitless)
    INTEGER(i_std),INTENT(in)      :: rest_id              !! Restart file identifier (unitless)
    INTEGER(i_std),INTENT(in)      :: hist_id              !! Access to history file (unitless)
    INTEGER(i_std),INTENT(in)      :: hist2_id             !! Access to history file 2 (unitless)
    REAL(r_std), INTENT(in)        :: lalo(nbpt,2)         !! Vector of latitude and longitudes (beware of the order !)

    INTEGER(i_std), INTENT(in)     :: neighbours(nbpt,8)   !! Vector of neighbours for each grid point 
                                                           !! (1=N, 2=NE, 3=E, 4=SE, 5=S, 6=SW, 7=W, 8=NW) (unitless)
    REAL(r_std), INTENT(in)        :: resolution(nbpt,2)   !! The size of each grid box in X and Y (m)
    REAL(r_std), INTENT(in)        :: contfrac(nbpt)       !! Fraction of land in each grid box (unitless;0-1)
    REAL(r_std), INTENT(in)        :: stempdiag(nbpt,nslm) !! Diagnostic soil temperature profile

    !! 0.2 Output variables
    REAL(r_std), INTENT(out)       :: returnflow(nbpt)     !! The water flow from lakes and swamps which returns to the grid box.
                                                           !! This water will go back into the hydrol or hydrolc module to allow re-evaporation (kg/m^2/dt)
    REAL(r_std), INTENT(out)       :: reinfiltration(nbpt) !! Water flow from ponds and floodplains which returns to the grid box (kg/m^2/dt)
    REAL(r_std), INTENT(out)       :: irrigation(nbpt)     !! Irrigation flux. This is the water taken from the reservoirs and beeing put into the upper layers of the soil (kg/m^2/dt)
    REAL(r_std), INTENT(out)       :: riverflow(nbpt)      !! Outflow of the major rivers. The flux will be located on the continental grid but this should be a coastal point (kg/dt)

    REAL(r_std), INTENT(out)       :: coastalflow(nbpt)    !! Outflow on coastal points by small basins. This is the water which flows in a disperse way into the ocean (kg/dt)
    REAL(r_std), INTENT(out)       :: flood_frac(nbpt)     !! Flooded fraction of the grid box (unitless;0-1)
    REAL(r_std), INTENT(out)       :: flood_res(nbpt)      !! Diagnostic of water amount in the floodplains reservoir (kg)

    !! 0.3 Local variables
    LOGICAL                        :: init_irrig           !! Logical to initialize the irrigation (true/false)
    LOGICAL                        :: init_flood           !! Logical to initialize the floodplains (true/false)
    LOGICAL                        :: init_swamp           !! Logical to initialize the swamps (true/false)

    !_ ================================================================================================================================

    !
    ! do initialisation
    !
    nbvmax = 440
    ! Here we will allocate the memory and get the fixed fields from the restart file.
    ! If the info is not found then we will compute the routing map.
    !
    CALL routing_interp_topo_v0_init_1 (kjit, nbpt, index, returnflow, reinfiltration, irrigation, &
                                        riverflow, coastalflow, flood_frac, flood_res, stempdiag, rest_id)

    routing_area => routing_area_loc  
    topo_resid => topo_resid_loc
    route_togrid => route_togrid_loc
    route_tobasin => route_tobasin_loc
    global_basinid => global_basinid_loc
    hydrodiag => hydrodiag_loc

    ! This routine computes the routing map if the route_togrid_glo is undefined. This means that the
    ! map has not been initialized during the restart process..
    !
    !! Reads in the map of the basins and flow directions to construct the catchments of each grid box
    !
    IF ( COUNT(route_togrid_glo .GE. undef_int) .GT. 0 ) THEN
       !       CALL routing_basins_p(nbpt, lalo, neighbours, resolution, contfrac)
    ENDIF
    !
    ! Do we have what we need if we want to do irrigation
    !! Initialisation of flags for irrigated land, flood plains and swamps
    !
    init_irrig = .FALSE.
    IF ( do_irrigation ) THEN 
       IF (COUNT(irrigated .GE. undef_sechiba-1) > 0) init_irrig = .TRUE.
    END IF

    init_flood = .FALSE.
    IF ( do_floodplains ) THEN
       IF (COUNT(floodplains .GE. undef_sechiba-1) > 0) init_flood = .TRUE.
    END IF

    init_swamp = .FALSE.
    IF ( doswamps ) THEN
       IF (COUNT(swamp .GE. undef_sechiba-1) > 0 ) init_swamp = .TRUE.
    END IF

    !! If we have irrigated land, flood plains or swamps then we need to interpolate the 0.5 degree
    !! base data set to the resolution of the model.

    IF ( init_irrig .OR. init_flood .OR. init_swamp ) THEN
       !       CALL routing_irrigmap(nbpt, index, lalo, neighbours, resolution, &
       !            contfrac, init_irrig, irrigated, init_flood, floodplains, init_swamp, swamp, hist_id, hist2_id)
    ENDIF

    IF ( do_irrigation ) THEN 
       CALL xios_orchidee_send_field("irrigmap",irrigated)

       WRITE(numout,*) 'Verification : range of irrigated : ', MINVAL(irrigated), MAXVAL(irrigated) 
       IF ( .NOT. almaoutput ) THEN
          CALL histwrite_p(hist_id, 'irrigmap', 1, irrigated, nbpt, index)
       ELSE
          CALL histwrite_p(hist_id, 'IrrigationMap', 1, irrigated, nbpt, index)
       ENDIF
       IF ( hist2_id > 0 ) THEN
          IF ( .NOT. almaoutput ) THEN
             CALL histwrite_p(hist2_id, 'irrigmap', 1, irrigated, nbpt, index)
          ELSE
             CALL histwrite_p(hist2_id, 'IrrigationMap', 1, irrigated, nbpt, index)
          ENDIF
       ENDIF
    ENDIF

    IF ( do_floodplains ) THEN
       CALL xios_orchidee_send_field("floodmap",floodplains)

       WRITE(numout,*) 'Verification : range of floodplains : ', MINVAL(floodplains), MAXVAL(floodplains) 
       IF ( .NOT. almaoutput ) THEN
          CALL histwrite_p(hist_id, 'floodmap', 1, floodplains, nbpt, index)
       ELSE
          CALL histwrite_p(hist_id, 'FloodplainsMap', 1, floodplains, nbpt, index)
       ENDIF
       IF ( hist2_id > 0 ) THEN
          IF ( .NOT. almaoutput ) THEN
             CALL histwrite_p(hist2_id, 'floodmap', 1, floodplains, nbpt, index)
          ELSE
             CALL histwrite_p(hist2_id, 'FloodplainsMap', 1, floodplains, nbpt, index)
          ENDIF
       ENDIF
    ENDIF

    IF ( doswamps ) THEN
       CALL xios_orchidee_send_field("swampmap",swamp)

       WRITE(numout,*) 'Verification : range of swamp : ', MINVAL(swamp), MAXVAL(swamp) 
       IF ( .NOT. almaoutput ) THEN
          CALL histwrite_p(hist_id, 'swampmap', 1, swamp, nbpt, index)
       ELSE
          CALL histwrite_p(hist_id, 'SwampMap', 1, swamp, nbpt, index)
       ENDIF
       IF ( hist2_id > 0 ) THEN
          IF ( .NOT. almaoutput ) THEN
             CALL histwrite_p(hist2_id, 'swampmap', 1, swamp, nbpt, index)
          ELSE
             CALL histwrite_p(hist2_id, 'SwampMap', 1, swamp, nbpt, index)
          ENDIF
       ENDIF
    ENDIF

    !! This routine gives a diagnostic of the basins used.
    !    CALL routing_diagnostic_p(nbpt, index, lalo, resolution, contfrac, hist_id, hist2_id)

    CALL routing_interp_topo_v0_init_2(nbpt, contfrac)

  END SUBROUTINE routing_interp_topo_v0_initialize


  !! ================================================================================================================================
  !! SUBROUTINE   : routing_interp_topo_v0_main 
  !!
  !>\BRIEF          This module routes the water over the continents (runoff and
  !!                drainage produced by the hydrolc or hydrol module) into the oceans. 
  !!
  !! DESCRIPTION (definitions, functional, design, flags):
  !! The routing scheme (Polcher, 2003) carries the water from the runoff and drainage simulated by SECHIBA
  !! to the ocean through reservoirs, with some delay. The routing scheme is based on
  !! a parametrization of the water flow on a global scale (Miller et al., 1994; Hagemann
  !! and Dumenil, 1998). Given the global map of the main watersheds (Oki et al., 1999;
  !! Fekete et al., 1999; Vorosmarty et al., 2000) which delineates the boundaries of subbasins
  !! and gives the eight possible directions of water flow within the pixel, the surface
  !! runoff and the deep drainage are routed to the ocean. The time-step of the routing is one day.
  !! The scheme also diagnoses how much water is retained in the foodplains and thus return to soil
  !! moisture or is taken out of the rivers for irrigation. \n
  !!
  !! RECENT CHANGE(S): None
  !!
  !! MAIN OUTPUT VARIABLE(S): 
  !! The result of the routing are 3 fluxes :
  !! - riverflow   : The water which flows out from the major rivers. The flux will be located
  !!                 on the continental grid but this should be a coastal point.
  !! - coastalflow : This is the water which flows in a disperse way into the ocean. Essentially these
  !!                 are the outflows from all of the small rivers.
  !! - returnflow  : This is the water which flows into a land-point - typically rivers which end in
  !!                 the desert. This water will go back into the hydrol module to allow re-evaporation.
  !! - irrigation  : This is water taken from the reservoir and is being put into the upper 
  !!                 layers of the soil.
  !! The two first fluxes are in kg/dt and the last two fluxes are in kg/(m^2dt).\n
  !!
  !! REFERENCE(S) :
  !! - Miller JR, Russell GL, Caliri G (1994)
  !!   Continental-scale river flow in climate models.
  !!   J. Clim., 7:914-928
  !! - Hagemann S and Dumenil L. (1998)
  !!   A parametrization of the lateral waterflow for the global scale.
  !!   Clim. Dyn., 14:17-31
  !! - Oki, T., T. Nishimura, and P. Dirmeyer (1999)
  !!   Assessment of annual runoff from land surface models using total runoff integrating pathways (TRIP)
  !!   J. Meteorol. Soc. Jpn., 77, 235-255
  !! - Fekete BM, Charles V, Grabs W (2000)
  !!   Global, composite runoff fields based on observed river discharge and simulated water balances.
  !!   Technical report, UNH/GRDC, Global Runoff Data Centre, Koblenz
  !! - Vorosmarty, C., B. Fekete, B. Meybeck, and R. Lammers (2000)
  !!   Global system of rivers: Its role in organizing continental land mass and defining land-to-ocean linkages
  !!   Global Biogeochem. Cycles, 14, 599-621
  !! - Vivant, A-C. (?? 2002)
  !!   Développement du schéma de routage et des plaines d'inondation, MSc Thesis, Paris VI University
  !! - J. Polcher (2003)
  !!   Les processus de surface a l'echelle globale et leurs interactions avec l'atmosphere
  !!   Habilitation a diriger les recherches, Paris VI University, 67pp.
  !!
  !! FLOWCHART    :
  !! \latexonly 
  !! \includegraphics[scale=0.75]{routing_main_flowchart.png}
  !! \endlatexonly
  !! \n
  !_ ================================================================================================================================


  SUBROUTINE routing_interp_topo_v0_main(kjit, nbpt, index, &
                                         lalo, neighbours, resolution, contfrac, totfrac_nobio, veget_max, floodout, runoff, &
                                         drainage, transpot, precip_rain, humrel, k_litt, flood_frac, flood_res, &
                                         stempdiag, reinf_slope, returnflow, reinfiltration, irrigation, riverflow, coastalflow, &
                                         rest_id, hist_id, hist2_id)

    IMPLICIT NONE

    !! 0 Variable and parameter description
    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)     :: kjit                 !! Time step number (unitless)
    INTEGER(i_std), INTENT(in)     :: nbpt                 !! Domain size (unitless)
    INTEGER(i_std),INTENT(in)      :: rest_id              !! Restart file identifier (unitless)
    INTEGER(i_std),INTENT(in)      :: hist_id              !! Access to history file (unitless)
    INTEGER(i_std),INTENT(in)      :: hist2_id             !! Access to history file 2 (unitless)
    INTEGER(i_std), INTENT(in)     :: index(nbpt)          !! Indices of the points on the map (unitless)
    REAL(r_std), INTENT(in)        :: lalo(nbpt,2)         !! Vector of latitude and longitudes (beware of the order !)
    INTEGER(i_std), INTENT(in)     :: neighbours(nbpt,8)   !! Vector of neighbours for each grid point (1=N, 2=NE, 3=E, 4=SE, 5=S, 6=SW, 7=W, 8=NW) (unitless)
    REAL(r_std), INTENT(in)        :: resolution(nbpt,2)   !! The size of each grid box in X and Y (m)
    REAL(r_std), INTENT(in)        :: contfrac(nbpt)       !! Fraction of land in each grid box (unitless;0-1)
    REAL(r_std), INTENT(in)        :: totfrac_nobio(nbpt)  !! Total fraction of no-vegetation (continental ice, lakes ...) (unitless;0-1)
    REAL(r_std), INTENT(in)        :: veget_max(nbpt,nvm)  !! Maximal fraction of vegetation (unitless;0-1)
    REAL(r_std), INTENT(in)        :: floodout(nbpt)       !! Grid-point flow out of floodplains (kg/m^2/dt)
    REAL(r_std), INTENT(in)        :: runoff(nbpt)         !! Grid-point runoff (kg/m^2/dt)
    REAL(r_std), INTENT(in)        :: drainage(nbpt)       !! Grid-point drainage (kg/m^2/dt)
    REAL(r_std), INTENT(in)        :: transpot(nbpt,nvm)   !! Potential transpiration of the vegetation (kg/m^2/dt)
    REAL(r_std), INTENT(in)        :: precip_rain(nbpt)    !! Rainfall (kg/m^2/dt)
    REAL(r_std), INTENT(in)        :: k_litt(nbpt)         !! Averaged conductivity for saturated infiltration in the 'litter' layer (kg/m^2/dt)
    REAL(r_std), INTENT(in)        :: humrel(nbpt,nvm)     !! Soil moisture stress, root extraction potential (unitless)
    REAL(r_std), INTENT(in)        :: stempdiag(nbpt,nslm) !! Diagnostic soil temperature profile
    REAL(r_std), INTENT(in)        :: reinf_slope(nbpt)    !! Coefficient which determines the reinfiltration ratio in the grid box due to flat areas (unitless;0-1)

    !! 0.2 Output variables
    REAL(r_std), INTENT(out)       :: returnflow(nbpt)     !! The water flow from lakes and swamps which returns to the grid box.
    !! This water will go back into the hydrol or hydrolc module to allow re-evaporation (kg/m^2/dt)
    REAL(r_std), INTENT(out)       :: reinfiltration(nbpt) !! Water flow from ponds and floodplains which returns to the grid box (kg/m^2/dt)
    REAL(r_std), INTENT(out)       :: irrigation(nbpt)     !! Irrigation flux. This is the water taken from the reservoirs and beeing put into the upper layers of the soil (kg/m^2/dt)
    REAL(r_std), INTENT(out)       :: riverflow(nbpt)      !! Outflow of the major rivers. The flux will be located on the continental grid but this should be a coastal point (kg/dt)
    REAL(r_std), INTENT(out)       :: coastalflow(nbpt)    !! Outflow on coastal points by small basins. This is the water which flows in a disperse way into the ocean (kg/dt)
    REAL(r_std), INTENT(out)       :: flood_frac(nbpt)     !! Flooded fraction of the grid box (unitless;0-1)
    REAL(r_std), INTENT(out)       :: flood_res(nbpt)      !! Diagnostic of water amount in the floodplains reservoir (kg)

    !! 0.3 Local variables
    CHARACTER(LEN=30)              :: var_name             !! To store variables names for I/O (unitless)
    REAL(r_std), DIMENSION(1)      :: tmp_day              !! 
    REAL(r_std), DIMENSION(nbpt)   :: return_lakes         !! Water from lakes flowing back into soil moisture (kg/m^2/dt)
    INTEGER(i_std)                 :: ig, jv               !! Indices (unitless)
    REAL(r_std), DIMENSION(nbpt)   :: tot_vegfrac_nowoody  !! Total fraction occupied by grass (0-1,unitless)
    REAL(r_std),PARAMETER :: delta_lon=(360./144.)/2.
    REAL(r_std),PARAMETER :: delta_lat=(180./144.)/2.

    !_ ================================================================================================================================

    !
    !! Computes the variables averaged between routing time steps and which will be used in subsequent calculations
    !
    floodout_mean(:) = floodout_mean(:) + floodout(:)
    runoff_mean(:) = runoff_mean(:) + runoff(:)
    drainage_mean(:) = drainage_mean(:) + drainage(:)
    floodtemp(:) = stempdiag(:,floodtemp_lev)
    precip_mean(:) =  precip_mean(:) + precip_rain(:)
    !
    !! Computes the total fraction occupied by the grasses and the crops for each grid cell
    tot_vegfrac_nowoody(:) = zero
    DO jv  = 1, nvm
       IF ( (jv /= ibare_sechiba) .AND. .NOT.(is_tree(jv)) ) THEN
          tot_vegfrac_nowoody(:) = tot_vegfrac_nowoody(:) + veget_max(:,jv) 
       END IF
    END DO

    DO ig = 1, nbpt
       IF ( tot_vegfrac_nowoody(ig) .GT. min_sechiba ) THEN
          DO jv = 1,nvm
             IF ( (jv /= ibare_sechiba) .AND. .NOT.(is_tree(jv)) ) THEN
                transpot_mean(ig) = transpot_mean(ig) + transpot(ig,jv) * veget_max(ig,jv)/tot_vegfrac_nowoody(ig)  
             END IF
          END DO
       ELSE
          IF (MAXVAL(veget_max(ig,2:nvm)) .GT. min_sechiba) THEN
             DO jv = 2, nvm
                transpot_mean(ig) = transpot_mean(ig) + transpot(ig,jv) * veget_max(ig,jv)/ SUM(veget_max(ig,2:nvm))
             ENDDO
          ENDIF
       ENDIF
    ENDDO

    !
    ! Averaged variables (i.e. *dt_sechiba/dt_routing). This accounts for the difference between the shorter
    ! timestep dt_sechiba of other parts of the model and the long dt_routing timestep (set to one day at present)
    !
    totnobio_mean(:) = totnobio_mean(:) + totfrac_nobio(:)*dt_sechiba/dt_routing
    k_litt_mean(:) = k_litt_mean(:) + k_litt(:)*dt_sechiba/dt_routing
    !
    ! Only potentially vegetated surfaces are taken into account. At the start of
    ! the growing seasons we will give more weight to these areas.
    !
    DO jv=2,nvm
       DO ig=1,nbpt
          humrel_mean(ig) = humrel_mean(ig) + humrel(ig,jv)*veget_max(ig,jv)*dt_sechiba/dt_routing
          vegtot_mean(ig) = vegtot_mean(ig) + veget_max(ig,jv)*dt_sechiba/dt_routing
       ENDDO
    ENDDO
    !
    time_counter = time_counter + dt_sechiba 
    !
    ! If the time has come we do the routing.
    !
    IF ( NINT(time_counter) .GE. NINT(dt_routing) ) THEN 
       !
       ! Check the water balance if needed
       !
       !       IF ( check_waterbal ) THEN
       !          CALL routing_waterbal(nbpt, .TRUE., floodout_mean, runoff_mean, drainage_mean, returnflow_mean, &
       !               & reinfiltration_mean, irrigation_mean, riverflow_mean, coastalflow_mean)
       !       ENDIF
       !
       ! Make sure we do not flood north of 49N as there freezing processes start to play a role and they
       ! are not yet well treated in ORCHIDEE.
       !
       DO ig=1,nbpt
          IF ( lalo(ig,1) > 49.0 ) THEN
             floodtemp(ig) = tp_00 - un
          ENDIF
       ENDDO
       !
       !! Computes the transport of water in the various reservoirs
       !
       !ym       runoff_mean=0
       !ym       drainage_mean=0
       !ym      DO ig=1,nbpt
       !ym         IF ( lalo(ig,1)-delta_lat < 0 .AND. lalo(ig,1)+delta_lat > 0 .AND. lalo(ig,2)-delta_lon < 32. .AND. lalo(ig,2)+delta_lon > 32.) runoff_mean(ig)=one_day/dt_routing        
       !ym       ENDDO

       CALL routing_interp_topo_v0_flow(nbpt, dt_routing, lalo, floodout_mean, runoff_mean, drainage_mean, &
            & vegtot_mean, totnobio_mean, transpot_mean, precip_mean, humrel_mean, k_litt_mean, floodtemp, reinf_slope, &
            & lakeinflow_mean, returnflow_mean, reinfiltration_mean, irrigation_mean, riverflow_mean, &
            & coastalflow_mean, hydrographs, slowflow_diag, flood_frac, flood_res)
       !
       !! Responsible for storing the water in lakes
       !
       CALL routing_interp_topo_v0_lake(nbpt, dt_routing, lakeinflow_mean, humrel_mean, contfrac, return_lakes)
       !
       returnflow_mean(:) = returnflow_mean(:) + return_lakes(:)
       !
       !! Check the water balance in the routing scheme
       !
       !       IF ( check_waterbal ) THEN
       !          CALL routing_waterbal(nbpt, .FALSE., floodout_mean, runoff_mean, drainage_mean, returnflow_mean, &
       !               & reinfiltration_mean, irrigation_mean, riverflow_mean, coastalflow_mean)
       !       ENDIF
       !
       time_counter = zero
       !
       floodout_mean(:) = zero
       runoff_mean(:) = zero
       drainage_mean(:) = zero
       transpot_mean(:) = zero
       precip_mean(:) = zero
       !
       humrel_mean(:) = zero
       totnobio_mean(:) = zero
       k_litt_mean(:) = zero
       vegtot_mean(:) = zero
       !
       ! Change the units of the routing fluxes from kg/dt_routing into kg/dt_sechiba
       !                                    and from m^3/dt_routing into m^3/dt_sechiba
       !
       returnflow_mean(:) = returnflow_mean(:)/dt_routing*dt_sechiba
       reinfiltration_mean(:) = reinfiltration_mean(:)/dt_routing*dt_sechiba
       irrigation_mean(:) = irrigation_mean(:)/dt_routing*dt_sechiba
       irrig_netereq(:) = irrig_netereq(:)/dt_routing*dt_sechiba
       !
       ! Change units as above but at the same time transform the kg/dt_sechiba to m^3/dt_sechiba
       !
       riverflow_mean(:) = riverflow_mean(:)/dt_routing*dt_sechiba/mille
       coastalflow_mean(:) = coastalflow_mean(:)/dt_routing*dt_sechiba/mille
       hydrographs(:) = hydrographs(:)/dt_routing*dt_sechiba/mille
       slowflow_diag(:) = slowflow_diag(:)/dt_routing*dt_sechiba/mille
    ENDIF

    !
    ! Return the fraction of routed water for this time step.
    !
    returnflow(:) = returnflow_mean(:)
    reinfiltration(:) = reinfiltration_mean(:)
    irrigation(:) = irrigation_mean(:)
    riverflow(:) = riverflow_mean(:)
    coastalflow(:) = coastalflow_mean(:)

    !
    ! Write diagnostics
    !

    CALL xios_orchidee_send_field("riversret",returnflow*one_day/dt_sechiba)
    CALL xios_orchidee_send_field("hydrographs",hydrographs/dt_sechiba)
    CALL xios_orchidee_send_field("slowflow",slowflow_diag/dt_sechiba) ! Qb in m3/s
    CALL xios_orchidee_send_field("reinfiltration",reinfiltration)
    CALL xios_orchidee_send_field("fastr",fast_diag)
    CALL xios_orchidee_send_field("slowr",slow_diag)
    CALL xios_orchidee_send_field("streamr",stream_diag)
    CALL xios_orchidee_send_field("laker",lake_diag)
    CALL xios_orchidee_send_field("pondr",pond_diag)
    CALL xios_orchidee_send_field("irrigation",irrigation*one_day/dt_sechiba)
    CALL xios_orchidee_send_field("netirrig",irrig_netereq*one_day/dt_sechiba)
    CALL xios_orchidee_send_field("floodh",flood_height)
    CALL xios_orchidee_send_field("floodr",flood_diag)
    delsurfstor = delsurfstor + flood_diag + pond_diag + lake_diag
    CALL xios_orchidee_send_field("SurfStor",flood_diag+pond_diag+lake_diag)

    ! Transform from kg/dt_sechiba into m^3/s
    CALL xios_orchidee_send_field("hydrographs",hydrographs/mille/dt_sechiba)
    CALL xios_orchidee_send_field("slowflow",slowflow_diag/mille/dt_sechiba) ! previous id name: Qb
    CALL xios_orchidee_send_field("coastalflow",coastalflow/dt_sechiba)
    CALL xios_orchidee_send_field("riverflow",riverflow/dt_sechiba)


    IF ( .NOT. almaoutput ) THEN
       !
       CALL histwrite_p(hist_id, 'riversret', kjit, returnflow, nbpt, index)
       IF (do_floodplains .OR. doponds) THEN
          CALL histwrite_p(hist_id, 'reinfiltration', kjit, reinfiltration, nbpt, index)
       ENDIF
       CALL histwrite_p(hist_id, 'hydrographs', kjit, hydrographs, nbpt, index)
       !
       CALL histwrite_p(hist_id, 'fastr', kjit, fast_diag, nbpt, index)
       CALL histwrite_p(hist_id, 'slowr', kjit, slow_diag, nbpt, index)
       CALL histwrite_p(hist_id, 'streamr', kjit, stream_diag, nbpt, index)
       IF ( do_floodplains ) THEN
          CALL histwrite_p(hist_id, 'floodr', kjit, flood_diag, nbpt, index)
          CALL histwrite_p(hist_id, 'floodh', kjit, flood_height, nbpt, index)
       ENDIF
       CALL histwrite_p(hist_id, 'pondr', kjit, pond_diag, nbpt, index)
       CALL histwrite_p(hist_id, 'lakevol', kjit, lake_diag, nbpt, index)
       !
       IF ( do_irrigation ) THEN
          CALL histwrite_p(hist_id, 'irrigation', kjit, irrigation, nbpt, index)
          CALL histwrite_p(hist_id, 'returnflow', kjit, returnflow, nbpt, index)
          CALL histwrite_p(hist_id, 'netirrig', kjit, irrig_netereq, nbpt, index)
       ENDIF
       !
    ELSE
       !
       delsurfstor = delsurfstor + flood_diag + pond_diag + lake_diag
       CALL histwrite_p(hist_id, 'DelSurfStor', kjit, delsurfstor, nbpt, index)
       CALL histwrite_p(hist_id, 'SurfStor', kjit, flood_diag+pond_diag+lake_diag, nbpt, index)
       CALL histwrite_p(hist_id, 'Dis', kjit, hydrographs, nbpt, index)
       IF ( do_irrigation ) THEN
          CALL histwrite_p(hist_id, 'Qirrig', kjit, irrigation, nbpt, index)
          CALL histwrite_p(hist_id, 'Qirrig_req', kjit, irrig_netereq, nbpt, index)
       ENDIF
       !
    ENDIF
    IF ( hist2_id > 0 ) THEN
       IF ( .NOT. almaoutput ) THEN
          !
          CALL histwrite_p(hist2_id, 'riversret', kjit, returnflow, nbpt, index)
          IF (do_floodplains .OR. doponds) THEN
             CALL histwrite_p(hist2_id, 'reinfiltration', kjit, reinfiltration, nbpt, index)
          ENDIF
          CALL histwrite_p(hist2_id, 'hydrographs', kjit, hydrographs, nbpt, index)
          !
          CALL histwrite_p(hist2_id, 'fastr', kjit, fast_diag, nbpt, index)
          CALL histwrite_p(hist2_id, 'slowr', kjit, slow_diag, nbpt, index)
          IF ( do_floodplains ) THEN
             CALL histwrite_p(hist2_id, 'floodr', kjit, flood_diag, nbpt, index)
             CALL histwrite_p(hist2_id, 'floodh', kjit, flood_height, nbpt, index)
          ENDIF
          CALL histwrite_p(hist2_id, 'pondr', kjit, pond_diag, nbpt, index)
          CALL histwrite_p(hist2_id, 'streamr', kjit, stream_diag, nbpt, index)
          CALL histwrite_p(hist2_id, 'lakevol', kjit, lake_diag, nbpt, index)
          !
          IF ( do_irrigation ) THEN
             CALL histwrite_p(hist2_id, 'irrigation', kjit, irrigation, nbpt, index)
             CALL histwrite_p(hist2_id, 'returnflow', kjit, returnflow, nbpt, index)
             CALL histwrite_p(hist2_id, 'netirrig', kjit, irrig_netereq, nbpt, index)
          ENDIF
          !
       ELSE
          !
          delsurfstor=delsurfstor + flood_diag + pond_diag + lake_diag
          CALL histwrite_p(hist2_id, 'DelSurfStor', kjit, delsurfstor, nbpt, index)
          CALL histwrite_p(hist2_id, 'SurfStor', kjit, flood_diag+pond_diag+lake_diag, nbpt, index)
          CALL histwrite_p(hist2_id, 'Dis', kjit, hydrographs, nbpt, index)
          !
       ENDIF
    ENDIF
    !
    !
  END SUBROUTINE routing_interp_topo_v0_main

  !!  =============================================================================================================================
  !! SUBROUTINE:         routing_interp_topo_v0_finalize
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

  SUBROUTINE routing_interp_topo_v0_finalize( kjit, nbpt, rest_id, flood_frac, flood_res )
    USE xios
    IMPLICIT NONE

    !! 0 Variable and parameter description
    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)     :: kjit                 !! Time step number (unitless)
    INTEGER(i_std), INTENT(in)     :: nbpt                 !! Domain size (unitless)
    INTEGER(i_std),INTENT(in)      :: rest_id              !! Restart file identifier (unitless)
    REAL(r_std), INTENT(in)        :: flood_frac(nbpt)     !! Flooded fraction of the grid box (unitless;0-1)
    REAL(r_std), INTENT(in)        :: flood_res(nbpt)      !! Diagnostic of water amount in the floodplains reservoir (kg)

    !! 0.2 Local variables
    REAL(r_std), DIMENSION(1)      :: tmp_day              

    !_ ================================================================================================================================

    IF (is_omp_root) THEN
       CALL xios_send_field("fast_reservoir_restart",fast_reservoir_r)
       CALL xios_send_field("slow_reservoir_restart",slow_reservoir_r)
       CALL xios_send_field("stream_reservoir_restart",stream_reservoir_r)
    ENDIF

    !
    ! Write restart variables
    !
    tmp_day(1) = time_counter
    IF (is_root_prc) CALL restput (rest_id, 'routingcounter', 1, 1, 1, kjit, tmp_day)

    CALL restput_p (rest_id, 'routingarea', nbp_glo, nbasmax, 1, kjit, routing_area, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'routetogrid', nbp_glo, nbasmax, 1, kjit, REAL(route_togrid,r_std), 'scatter', &
         nbp_glo, index_g)
    CALL restput_p (rest_id, 'routetobasin', nbp_glo, nbasmax, 1, kjit, REAL(route_tobasin,r_std), 'scatter', &
         nbp_glo, index_g)
    CALL restput_p (rest_id, 'basinid', nbp_glo, nbasmax, 1, kjit, REAL(global_basinid,r_std), 'scatter', &
         nbp_glo, index_g)
    CALL restput_p (rest_id, 'topoindex', nbp_glo, nbasmax, 1, kjit, topo_resid, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'fastres', nbp_glo, nbasmax, 1, kjit, fast_reservoir, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'slowres', nbp_glo, nbasmax, 1, kjit, slow_reservoir, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'streamres', nbp_glo, nbasmax, 1, kjit, stream_reservoir, 'scatter',nbp_glo,index_g)
    CALL restput_p (rest_id, 'floodres', nbp_glo, nbasmax, 1, kjit, flood_reservoir, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'floodh', nbp_glo, 1, 1, kjit, flood_height, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'flood_frac_bas', nbp_glo, nbasmax, 1, kjit, flood_frac_bas, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'pond_frac', nbp_glo, 1, 1, kjit, pond_frac, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'flood_frac', nbp_glo, 1, 1, kjit, flood_frac, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'flood_res', nbp_glo, 1, 1, kjit, flood_res, 'scatter', nbp_glo, index_g)

    CALL restput_p (rest_id, 'lakeres', nbp_glo, 1, 1, kjit, lake_reservoir, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'pondres', nbp_glo, 1, 1, kjit, pond_reservoir, 'scatter',  nbp_glo, index_g)

    CALL restput_p (rest_id, 'lakeinflow', nbp_glo, 1, 1, kjit, lakeinflow_mean, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'returnflow', nbp_glo, 1, 1, kjit, returnflow_mean, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'reinfiltration', nbp_glo, 1, 1, kjit, reinfiltration_mean, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'riverflow', nbp_glo, 1, 1, kjit, riverflow_mean, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'coastalflow', nbp_glo, 1, 1, kjit, coastalflow_mean, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'hydrographs', nbp_glo, 1, 1, kjit, hydrographs, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'slowflow_diag', nbp_glo, 1, 1, kjit, slowflow_diag, 'scatter',  nbp_glo, index_g)
    !
    ! Keep track of the accumulated variables
    !
    CALL restput_p (rest_id, 'floodout_route', nbp_glo, 1, 1, kjit, floodout_mean, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'runoff_route', nbp_glo, 1, 1, kjit, runoff_mean, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'drainage_route', nbp_glo, 1, 1, kjit, drainage_mean, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'transpot_route', nbp_glo, 1, 1, kjit, transpot_mean, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'precip_route', nbp_glo, 1, 1, kjit, precip_mean, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'humrel_route', nbp_glo, 1, 1, kjit, humrel_mean, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'totnobio_route', nbp_glo, 1, 1, kjit, totnobio_mean, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'k_litt_route', nbp_glo, 1, 1, kjit, k_litt_mean, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'vegtot_route', nbp_glo, 1, 1, kjit, vegtot_mean, 'scatter',  nbp_glo, index_g)

    IF ( do_irrigation ) THEN
       CALL restput_p (rest_id, 'irrigated', nbp_glo, 1, 1, kjit, irrigated, 'scatter',  nbp_glo, index_g)
       CALL restput_p (rest_id, 'irrigation', nbp_glo, 1, 1, kjit, irrigation_mean, 'scatter',  nbp_glo, index_g)
    ENDIF

    IF ( do_floodplains ) THEN
       CALL restput_p (rest_id, 'floodplains', nbp_glo, 1, 1, kjit, floodplains, 'scatter',  nbp_glo, index_g)
    ENDIF
    IF ( doswamps ) THEN
       CALL restput_p (rest_id, 'swamp', nbp_glo, 1, 1, kjit, swamp, 'scatter',  nbp_glo, index_g)
    ENDIF



  END SUBROUTINE routing_interp_topo_v0_finalize

  !! ================================================================================================================================
  !! SUBROUTINE 	: routing_interp_topo_v0_init_1
  !!
  !>\BRIEF         This subroutine allocates the memory and get the fixed fields from the restart file.
  !!
  !! DESCRIPTION:	  Privat subroutine to the module routing_interp_topo_v0. This subroutine is called in the begining 
  !!                      of routing_interp_topo_v0_initialize
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

  SUBROUTINE routing_interp_topo_v0_init_1(kjit, nbpt, index, returnflow, reinfiltration, irrigation, &
                                           riverflow, coastalflow, flood_frac, flood_res, stempdiag, rest_id)
    
    IMPLICIT NONE
        
    !! 0 Variable and parameter description
    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                   :: kjit           !! Time step number (unitless)
    INTEGER(i_std), INTENT(in)                   :: nbpt           !! Domain size (unitless)
    INTEGER(i_std), DIMENSION (nbpt), INTENT(in) :: index          !! Indices of the points on the map (unitless)
    REAL(r_std), DIMENSION(nbpt,nslm),INTENT(in) :: stempdiag      !! Temperature profile in soil
    INTEGER(i_std), INTENT(in)                   :: rest_id        !! Restart file identifier (unitless)
    
    !! 0.2 Output variables
    REAL(r_std), DIMENSION (nbpt),INTENT(out)    :: returnflow     !! The water flow from lakes and swamps which returns into the grid box.
                                                                   !! This water will go back into the hydrol or hydrolc module to allow re-evaporation (kg/m^2/dt)
    REAL(r_std), DIMENSION (nbpt),INTENT(out)    :: reinfiltration !! Water flow from ponds and floodplains which returns to the grid box (kg/m^2/dt)
    REAL(r_std), DIMENSION (nbpt),INTENT(out)    :: irrigation     !! Irrigation flux. This is the water taken from the reservoirs and beeing put into the upper layers of the soil.(kg/m^2/dt)
    REAL(r_std), DIMENSION (nbpt),INTENT(out)    :: riverflow      !! Outflow of the major rivers. The flux will be located on the continental grid but this should be a coastal point (kg/dt)
    REAL(r_std), DIMENSION (nbpt),INTENT(out)    :: coastalflow    !! Outflow on coastal points by small basins. This is the water which flows in a disperse way into the ocean (kg/dt)
    REAL(r_std), DIMENSION (nbpt),INTENT(out)    :: flood_frac     !! Flooded fraction of the grid box (unitless;0-1)
    REAL(r_std), DIMENSION (nbpt),INTENT(out)    :: flood_res      !! Diagnostic of water amount in the floodplains reservoir (kg)
    
    !! 0.3 Local variables
    CHARACTER(LEN=80)                            :: var_name       !! To store variables names for I/O (unitless)
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)     :: tmp_real_g     !! A temporary real array for the integers
    REAL(r_std), DIMENSION(1)                    :: tmp_day        !!
    REAL(r_std)                                  :: ratio          !! Diagnostic ratio to check that dt_routing is a multiple of dt_sechiba (unitless)
    REAL(r_std)                                  :: totarea        !! Total area of basin (m^2)
    INTEGER(i_std)                               :: ier, ig, ib, ipn(1) !! Indices (unitless)

    !_ ================================================================================================================================
    !
    !
    ! These variables will require the configuration infrastructure
    !
    !Config Key   = DT_ROUTING 
    !Config If    = RIVER_ROUTING
    !Config Desc  = Time step of the routing scheme
    !Config Def   = one_day
    !Config Help  = This values gives the time step in seconds of the routing scheme. 
    !Config         It should be multiple of the main time step of ORCHIDEE. One day
    !Config         is a good value.
    !Config Units = [seconds]
    !
    dt_routing = one_day
    CALL getin_p('DT_ROUTING', dt_routing)
    !
    !Config Key   = ROUTING_RIVERS
    !Config If    = RIVER_ROUTING
    !Config Desc  = Number of rivers 
    !Config Def   = 50
    !Config Help  = This parameter chooses the number of largest river basins
    !Config         which should be treated as independently as rivers and not
    !Config         flow into the oceans as diffusion coastal flow.
    !Config Units = [-]
    num_largest = 50
    CALL getin_p('ROUTING_RIVERS', num_largest)
    !
    !Config Key   = DO_FLOODINFILT
    !Config Desc  = Should floodplains reinfiltrate into the soil 
    !Config If    = RIVER_ROUTING
    !Config Def   = n
    !Config Help  = This parameters allows the user to ask the model
    !Config         to take into account the flood plains reinfiltration 
    !Config         into the soil moisture. It then can go 
    !Config         back to the slow and fast reservoirs
    !Config Units = [FLAG]
    !
    dofloodinfilt = .FALSE.
    CALL getin_p('DO_FLOODINFILT', dofloodinfilt)
    !
    !Config Key   = DO_SWAMPS
    !Config Desc  = Should we include swamp parameterization 
    !Config If    = RIVER_ROUTING
    !Config Def   = n
    !Config Help  = This parameters allows the user to ask the model
    !Config         to take into account the swamps and return 
    !Config         the water into the bottom of the soil. It then can go 
    !Config         back to the atmopshere. This tried to simulate 
    !Config         internal deltas of rivers.
    !Config Units = [FLAG]
    !
    doswamps = .FALSE.
    CALL getin_p('DO_SWAMPS', doswamps)
    !
    !Config Key   = DO_PONDS
    !Config Desc  = Should we include ponds 
    !Config If    = RIVER_ROUTING
    !Config Def   = n
    !Config Help  = This parameters allows the user to ask the model
    !Config         to take into account the ponds and return 
    !Config         the water into the soil moisture. It then can go 
    !Config         back to the atmopshere. This tried to simulate 
    !Config         little ponds especially in West Africa.
    !Config Units = [FLAG]
    !
    doponds = .FALSE.
    CALL getin_p('DO_PONDS', doponds)

    !Config Key   = SLOW_TCST
    !Config Desc  = Time constant for the slow reservoir 
    !Config If    = RIVER_ROUTING 
    !Config Def   = 25.0
    !Config Help  = This parameters allows the user to fix the 
    !Config         time constant (in days) of the slow reservoir
    !Config         in order to get better river flows for 
    !Config         particular regions.
    !Config Units = [days]
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
    CALL getin_p('FLOOD_TCST', flood_tcst)

    !Config Key   = SWAMP_CST
    !Config Desc  = Fraction of the river that flows back to swamps 
    !Config If    = RIVER_ROUTING
    !Config Def   = 0.2
    !Config Help  = This parameters allows the user to fix the 
    !Config         fraction of the river transport
    !Config         that flows to swamps
    !Config Units = [-]
    CALL getin_p('SWAMP_CST', swamp_cst)

    !Config Key   = FLOOD_BETA
    !Config Desc  = Parameter to fix the shape of the floodplain  
    !Config If    = RIVER_ROUTING
    !Config Def   = 2.0
    !Config Help  = Parameter to fix the shape of the floodplain
    !Config         (>1 for convex edges, <1 for concave edges)
    !Config Units = [-] 
    CALL getin_p("FLOOD_BETA", beta)

    !Config Key   = POND_BETAP
    !Config Desc  = Ratio of the basin surface intercepted by ponds and the maximum surface of ponds
    !Config If    = RIVER_ROUTING
    !Config Def   = 0.5
    !Config Help  = 
    !Config Units = [-] 
    CALL getin_p("POND_BETAP", betap)    

    !Config Key   = FLOOD_CRI
    !Config Desc  = Potential height for which all the basin is flooded
    !Config If    = DO_FLOODPLAINS or DO_PONDS
    !Config Def   = 2000.
    !Config Help  = 
    !Config Units = [mm] 
    CALL getin_p("FLOOD_CRI", floodcri)

    !Config Key   = POND_CRI
    !Config Desc  = Potential height for which all the basin is a pond
    !Config If    = DO_FLOODPLAINS or DO_PONDS
    !Config Def   = 2000.
    !Config Help  = 
    !Config Units = [mm] 
    CALL getin_p("POND_CRI", pondcri)

    !Config Key   = MAX_LAKE_RESERVOIR
    !Config Desc  = Maximum limit of water in lake_reservoir
    !Config If    = RIVER_ROUTING
    !Config Def   = 7000
    !Config Help  = 
    !Config Units = [kg/m2(routing area)] 
    max_lake_reservoir = 7000
    CALL getin_p("MAX_LAKE_RESERVOIR", max_lake_reservoir)

    ! In order to simplify the time cascade check that dt_routing
    ! is a multiple of dt_sechiba
    !
    ratio = dt_routing/dt_sechiba
    IF ( ABS(NINT(ratio) - ratio) .GT. 10*EPSILON(ratio)) THEN
       WRITE(numout,*) 'WARNING -- WARNING -- WARNING -- WARNING'
       WRITE(numout,*) "The chosen time step for the routing is not a multiple of the"
       WRITE(numout,*) "main time step of the model. We will change dt_routing so that"
       WRITE(numout,*) "this condition os fulfilled"
       dt_routing = NINT(ratio) * dt_sechiba
       WRITE(numout,*) 'THE NEW DT_ROUTING IS : ', dt_routing
    ENDIF
    !
    IF ( dt_routing .LT. dt_sechiba) THEN
       WRITE(numout,*) 'WARNING -- WARNING -- WARNING -- WARNING'
       WRITE(numout,*) 'The routing timestep can not be smaller than the one'
       WRITE(numout,*) 'of the model. We reset its value to the model''s timestep.'
       WRITE(numout,*) 'The old DT_ROUTING is : ', dt_routing
       dt_routing = dt_sechiba
       WRITE(numout,*) 'THE NEW DT_ROUTING IS : ', dt_routing
    ENDIF
    !
    var_name ="routingcounter"
    IF (is_root_prc) THEN
       CALL ioconf_setatt('UNITS', 's')
       CALL ioconf_setatt('LONG_NAME','Time counter for the routing scheme')
       CALL restget (rest_id, var_name, 1, 1, 1, kjit, .TRUE., tmp_day)
       IF (tmp_day(1) == val_exp) THEN
          time_counter = zero
       ELSE
          time_counter = tmp_day(1) 
       ENDIF
    ENDIF
    CALL bcast(time_counter)


    ALLOCATE (routing_area_loc(nbpt,nbasmax), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for routing_area_loc','','')

    ALLOCATE (routing_area_glo(nbp_glo,nbasmax))
    IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for routing_area_glo','','')
    var_name = 'routingarea'
    IF (is_root_prc) THEN
       CALL ioconf_setatt('UNITS', 'm^2')
       CALL ioconf_setatt('LONG_NAME','Area of basin')
       CALL restget (rest_id, var_name, nbp_glo, nbasmax, 1, kjit, .TRUE., routing_area_glo, "gather", nbp_glo, index_g)
    ENDIF
    CALL scatter(routing_area_glo,routing_area_loc)
    routing_area=>routing_area_loc

    ALLOCATE (tmp_real_g(nbp_glo,nbasmax), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for tmp_real_g','','')

    ALLOCATE (route_togrid_loc(nbpt,nbasmax), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for route_togrid_loc','','')
    ALLOCATE (route_togrid_glo(nbp_glo,nbasmax), stat=ier)      ! used in global in routing_interp_topo_v0_flow
    IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for route_togrid_glo','','')

    IF (is_root_prc) THEN
       var_name = 'routetogrid'
       CALL ioconf_setatt('UNITS', '-')
       CALL ioconf_setatt('LONG_NAME','Grid into which the basin flows')
       CALL restget (rest_id, var_name, nbp_glo, nbasmax, 1, kjit, .TRUE., tmp_real_g, "gather", nbp_glo, index_g)
       route_togrid_glo(:,:) = undef_int
       WHERE ( tmp_real_g .LT. val_exp )
          route_togrid_glo = NINT(tmp_real_g)
       ENDWHERE
    ENDIF
    CALL bcast(route_togrid_glo)                      ! used in global in routing_interp_topo_v0_flow
    CALL scatter(route_togrid_glo,route_togrid_loc)
    route_togrid=>route_togrid_loc
    !
    ALLOCATE (route_tobasin_loc(nbpt,nbasmax), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for route_tobasin_loc','','')

    ALLOCATE (route_tobasin_glo(nbp_glo,nbasmax), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for route_tobasin_glo','','')

    IF (is_root_prc) THEN
       var_name = 'routetobasin'
       CALL ioconf_setatt('UNITS', '-')
       CALL ioconf_setatt('LONG_NAME','Basin in to which the water goes')
       CALL restget (rest_id, var_name, nbp_glo, nbasmax, 1, kjit, .TRUE., tmp_real_g, "gather", nbp_glo, index_g)
       route_tobasin_glo = undef_int
       WHERE ( tmp_real_g .LT. val_exp )
          route_tobasin_glo = NINT(tmp_real_g)
       ENDWHERE
    ENDIF
    CALL scatter(route_tobasin_glo,route_tobasin_loc)
    route_tobasin=>route_tobasin_loc
    !
    ! nbintobasin
    !
    ALLOCATE (route_nbintobas_loc(nbpt,nbasmax), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for route_nbintobas_loc','','')
    ALLOCATE (route_nbintobas_glo(nbp_glo,nbasmax), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for route_nbintobas_glo','','')

    IF (is_root_prc) THEN
       var_name = 'routenbintobas'
       CALL ioconf_setatt('UNITS', '-')
       CALL ioconf_setatt('LONG_NAME','Number of basin into current one')
       CALL restget (rest_id, var_name, nbp_glo, nbasmax, 1, kjit, .TRUE., tmp_real_g, "gather", nbp_glo, index_g)
       route_nbintobas_glo = undef_int
       WHERE ( tmp_real_g .LT. val_exp )
          route_nbintobas_glo = NINT(tmp_real_g)
       ENDWHERE
    ENDIF
    CALL scatter(route_nbintobas_glo,route_nbintobas_loc)
    route_nbintobas=>route_nbintobas_loc
    !
    ALLOCATE (global_basinid_loc(nbpt,nbasmax), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for global_basinid_loc','','')
    ALLOCATE (global_basinid_glo(nbp_glo,nbasmax), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for global_basinid_glo','','')

    IF (is_root_prc) THEN
       var_name = 'basinid'
       CALL ioconf_setatt('UNITS', '-')
       CALL ioconf_setatt('LONG_NAME','ID of basin')
       CALL restget (rest_id, var_name, nbp_glo, nbasmax, 1, kjit, .TRUE., tmp_real_g, "gather", nbp_glo, index_g)
       global_basinid_glo = undef_int
       WHERE ( tmp_real_g .LT. val_exp )
          global_basinid_glo = NINT(tmp_real_g)
       ENDWHERE
    ENDIF
    CALL scatter(global_basinid_glo,global_basinid_loc)
    global_basinid=>global_basinid_loc
    !
    ALLOCATE (topo_resid_loc(nbpt,nbasmax), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for topo_resid_loc','','')
    ALLOCATE (topo_resid_glo(nbp_glo,nbasmax), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for topo_resid_glo','','')

    IF (is_root_prc) THEN
       var_name = 'topoindex'
       CALL ioconf_setatt('UNITS', 'm')
       CALL ioconf_setatt('LONG_NAME','Topographic index of the residence time')
       CALL restget (rest_id, var_name, nbp_glo, nbasmax, 1, kjit, .TRUE., topo_resid_glo, "gather", nbp_glo, index_g)
    ENDIF
    CALL scatter(topo_resid_glo,topo_resid_loc)
    topo_resid=>topo_resid_loc

    ALLOCATE (fast_reservoir(nbpt,nbasmax), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for fast_reservoir','','')
    var_name = 'fastres'
    CALL ioconf_setatt_p('UNITS', 'Kg')
    CALL ioconf_setatt_p('LONG_NAME','Water in the fast reservoir')
    CALL restget_p (rest_id, var_name, nbp_glo, nbasmax, 1, kjit, .TRUE., fast_reservoir, "gather", nbp_glo, index_g)
    CALL setvar_p (fast_reservoir, val_exp, 'NO_KEYWORD', zero)

    ALLOCATE (slow_reservoir(nbpt,nbasmax), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for slow_reservoir','','')
    var_name = 'slowres'
    CALL ioconf_setatt_p('UNITS', 'Kg')
    CALL ioconf_setatt_p('LONG_NAME','Water in the slow reservoir')
    CALL restget_p (rest_id, var_name, nbp_glo, nbasmax, 1, kjit, .TRUE., slow_reservoir, "gather", nbp_glo, index_g)
    CALL setvar_p (slow_reservoir, val_exp, 'NO_KEYWORD', zero)

    ALLOCATE (stream_reservoir(nbpt,nbasmax), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for stream_reservoir','','')
    var_name = 'streamres'
    CALL ioconf_setatt_p('UNITS', 'Kg')
    CALL ioconf_setatt_p('LONG_NAME','Water in the stream reservoir')
    CALL restget_p (rest_id, var_name, nbp_glo, nbasmax, 1, kjit, .TRUE., stream_reservoir, "gather", nbp_glo, index_g)
    CALL setvar_p (stream_reservoir, val_exp, 'NO_KEYWORD', zero)

    ALLOCATE (flood_reservoir(nbpt,nbasmax), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for flood_reservoir','','')
    var_name = 'floodres'
    CALL ioconf_setatt_p('UNITS', 'Kg')
    CALL ioconf_setatt_p('LONG_NAME','Water in the flood reservoir')
    CALL restget_p (rest_id, var_name, nbp_glo, nbasmax, 1, kjit, .TRUE., flood_reservoir, "gather", nbp_glo, index_g)
    CALL setvar_p (flood_reservoir, val_exp, 'NO_KEYWORD', zero)

    ALLOCATE (flood_frac_bas(nbpt,nbasmax), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for flood_frac_bas','','')
    var_name = 'flood_frac_bas'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Flooded fraction per basin')
    CALL restget_p (rest_id, var_name, nbp_glo, nbasmax, 1, kjit, .TRUE., flood_frac_bas, "gather", nbp_glo, index_g)
    CALL setvar_p (flood_frac_bas, val_exp, 'NO_KEYWORD', zero)

    ALLOCATE (flood_height(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for flood_height','','')
    var_name = 'floodh'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., flood_height, "gather", nbp_glo, index_g)
    CALL setvar_p (flood_height, val_exp, 'NO_KEYWORD', zero)

    ALLOCATE (pond_frac(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for pond_frac','','')
    var_name = 'pond_frac'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Pond fraction per grid box')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., pond_frac, "gather", nbp_glo, index_g)
    CALL setvar_p (pond_frac, val_exp, 'NO_KEYWORD', zero)

    var_name = 'flood_frac'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Flooded fraction per grid box')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., flood_frac, "gather", nbp_glo, index_g)
    CALL setvar_p (flood_frac, val_exp, 'NO_KEYWORD', zero)

    var_name = 'flood_res'
    CALL ioconf_setatt_p('UNITS','mm')
    CALL ioconf_setatt_p('LONG_NAME','Flooded quantity (estimation)')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., flood_res, "gather", nbp_glo, index_g)
    CALL setvar_p (flood_res, val_exp, 'NO_KEYWORD', zero)
    !    flood_res = zero

    ALLOCATE (lake_reservoir(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for lake_reservoir','','')
    var_name = 'lakeres'
    CALL ioconf_setatt_p('UNITS', 'Kg')
    CALL ioconf_setatt_p('LONG_NAME','Water in the lake reservoir')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., lake_reservoir, "gather", nbp_glo, index_g)
    CALL setvar_p (lake_reservoir, val_exp, 'NO_KEYWORD', zero)

    ALLOCATE (pond_reservoir(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for pond_reservoir','','')
    var_name = 'pondres'
    CALL ioconf_setatt_p('UNITS', 'Kg')
    CALL ioconf_setatt_p('LONG_NAME','Water in the pond reservoir')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., pond_reservoir, "gather", nbp_glo, index_g)
    CALL setvar_p (pond_reservoir, val_exp, 'NO_KEYWORD', zero)
    !
    ! Map of irrigated areas
    !
    IF ( do_irrigation ) THEN
       ALLOCATE (irrigated(nbpt), stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for irrigated','','')
       var_name = 'irrigated'
       CALL ioconf_setatt_p('UNITS', 'm^2')
       CALL ioconf_setatt_p('LONG_NAME','Surface of irrigated area')
       CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., irrigated, "gather", nbp_glo, index_g)
       CALL setvar_p (irrigated, val_exp, 'NO_KEYWORD', undef_sechiba)
    ENDIF

    IF ( do_floodplains ) THEN
       ALLOCATE (floodplains(nbpt), stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for floodplains','','')
       var_name = 'floodplains'
       CALL ioconf_setatt_p('UNITS', 'm^2')
       CALL ioconf_setatt_p('LONG_NAME','Surface which can be flooded')
       CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., floodplains, "gather", nbp_glo, index_g)
       CALL setvar_p (floodplains, val_exp, 'NO_KEYWORD', undef_sechiba)
    ENDIF
    IF ( doswamps ) THEN
       ALLOCATE (swamp(nbpt), stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for swamp','','')
       var_name = 'swamp'
       CALL ioconf_setatt_p('UNITS', 'm^2')
       CALL ioconf_setatt_p('LONG_NAME','Surface which can become swamp')
       CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., swamp, "gather", nbp_glo, index_g)
       CALL setvar_p (swamp, val_exp, 'NO_KEYWORD', undef_sechiba)
    ENDIF
    !
    ! Put into the restart file the fluxes so that they can be regenerated at restart.
    !
    ALLOCATE (lakeinflow_mean(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for lakeinflow_mean','','')
    var_name = 'lakeinflow'
    CALL ioconf_setatt_p('UNITS', 'Kg/dt')
    CALL ioconf_setatt_p('LONG_NAME','Lake inflow')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., lakeinflow_mean, "gather", nbp_glo, index_g)
    CALL setvar_p (lakeinflow_mean, val_exp, 'NO_KEYWORD', zero)

    ALLOCATE (returnflow_mean(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for returnflow_mean','','')
    var_name = 'returnflow'
    CALL ioconf_setatt_p('UNITS', 'Kg/m^2/dt')
    CALL ioconf_setatt_p('LONG_NAME','Deep return flux')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., returnflow_mean, "gather", nbp_glo, index_g)
    CALL setvar_p (returnflow_mean, val_exp, 'NO_KEYWORD', zero)
    returnflow(:) = returnflow_mean(:)

    ALLOCATE (reinfiltration_mean(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for reinfiltration_mean','','')
    var_name = 'reinfiltration'
    CALL ioconf_setatt_p('UNITS', 'Kg/m^2/dt')
    CALL ioconf_setatt_p('LONG_NAME','Top return flux')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., reinfiltration_mean, "gather", nbp_glo, index_g)
    CALL setvar_p (reinfiltration_mean, val_exp, 'NO_KEYWORD', zero)
    reinfiltration(:) = reinfiltration_mean(:)

    ALLOCATE (irrigation_mean(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for irrigation_mean','','')
    ALLOCATE (irrig_netereq(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for irrig_netereq','','')
    irrig_netereq(:) = zero

    IF ( do_irrigation ) THEN
       var_name = 'irrigation'
       CALL ioconf_setatt_p('UNITS', 'Kg/dt')
       CALL ioconf_setatt_p('LONG_NAME','Artificial irrigation flux')
       CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., irrigation_mean, "gather", nbp_glo, index_g)
       CALL setvar_p (irrigation_mean, val_exp, 'NO_KEYWORD', zero)
    ELSE
       irrigation_mean(:) = zero
    ENDIF
    irrigation(:) = irrigation_mean(:) 

    ALLOCATE (riverflow_mean(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for riverflow_mean','','')
    var_name = 'riverflow'
    CALL ioconf_setatt_p('UNITS', 'Kg/dt')
    CALL ioconf_setatt_p('LONG_NAME','River flux into the sea')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., riverflow_mean, "gather", nbp_glo, index_g)
    CALL setvar_p (riverflow_mean, val_exp, 'NO_KEYWORD', zero)
    riverflow(:) = riverflow_mean(:)

    ALLOCATE (coastalflow_mean(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for coastalflow_mean','','')
    var_name = 'coastalflow'
    CALL ioconf_setatt_p('UNITS', 'Kg/dt')
    CALL ioconf_setatt_p('LONG_NAME','Diffuse flux into the sea')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., coastalflow_mean, "gather", nbp_glo, index_g)
    CALL setvar_p (coastalflow_mean, val_exp, 'NO_KEYWORD', zero)
    coastalflow(:) = coastalflow_mean(:)

    ! Locate it at the 2m level
    ipn = MINLOC(ABS(zlt-2))
    floodtemp_lev = ipn(1)
    ALLOCATE (floodtemp(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for floodtemp','','')
    floodtemp(:) = stempdiag(:,floodtemp_lev)

    ALLOCATE(hydrographs(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for hydrographs','','')
    var_name = 'hydrographs'
    CALL ioconf_setatt_p('UNITS', 'm^3/dt')
    CALL ioconf_setatt_p('LONG_NAME','Hydrograph at outlow of grid')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., hydrographs, "gather", nbp_glo, index_g)
    CALL setvar_p (hydrographs, val_exp, 'NO_KEYWORD', zero)

    ALLOCATE(slowflow_diag(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for slowflow_diag','','')
    var_name = 'slowflow_diag'
    CALL ioconf_setatt_p('UNITS', 'm^3/dt')
    CALL ioconf_setatt_p('LONG_NAME','Slowflow hydrograph at outlow of grid')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE.,slowflow_diag, "gather", nbp_glo, index_g)
    CALL setvar_p (slowflow_diag, val_exp, 'NO_KEYWORD', zero)

    !
    ! The diagnostic variables, they are initialized from the above restart variables.
    !
    ALLOCATE(fast_diag(nbpt), slow_diag(nbpt), stream_diag(nbpt), flood_diag(nbpt), &
         & pond_diag(nbpt), lake_diag(nbpt), delsurfstor(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for fast_diag,..','','')

    fast_diag(:) = zero
    slow_diag(:) = zero
    stream_diag(:) = zero
    flood_diag(:) = zero
    pond_diag(:) = zero
    lake_diag(:) = zero
    delsurfstor(:) = zero

    DO ig=1,nbpt
       totarea = zero
       DO ib=1,nbasmax
          totarea = totarea + routing_area(ig,ib)
          fast_diag(ig) = fast_diag(ig) + fast_reservoir(ig,ib)
          slow_diag(ig) = slow_diag(ig) + slow_reservoir(ig,ib)
          stream_diag(ig) = stream_diag(ig) + stream_reservoir(ig,ib)
          flood_diag(ig) = flood_diag(ig) + flood_reservoir(ig,ib)
       ENDDO
       !
       fast_diag(ig) = fast_diag(ig)/totarea
       slow_diag(ig) = slow_diag(ig)/totarea
       stream_diag(ig) = stream_diag(ig)/totarea
       flood_diag(ig) = flood_diag(ig)/totarea
       !
       ! This is the volume of the lake scaled to the entire grid.
       ! It would be batter to scale it to the size of the lake
       ! but this information is not yet available.
       !
       lake_diag(ig) = lake_reservoir(ig)/totarea
       !
    ENDDO
    !
    ! Get from the restart the fluxes we accumulated.
    !
    ALLOCATE (floodout_mean(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for floodout_mean','','')
    var_name = 'floodout_route'
    CALL ioconf_setatt_p('UNITS', 'Kg')
    CALL ioconf_setatt_p('LONG_NAME','Accumulated flow out of floodplains for routing')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., floodout_mean, "gather", nbp_glo, index_g)
    CALL setvar_p (floodout_mean, val_exp, 'NO_KEYWORD', zero)

    ALLOCATE (runoff_mean(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for runoff_mean','','')
    var_name = 'runoff_route'
    CALL ioconf_setatt_p('UNITS', 'Kg')
    CALL ioconf_setatt_p('LONG_NAME','Accumulated runoff for routing')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., runoff_mean, "gather", nbp_glo, index_g)
    CALL setvar_p (runoff_mean, val_exp, 'NO_KEYWORD', zero)

    ALLOCATE(drainage_mean(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for drainage_mean','','')
    var_name = 'drainage_route'
    CALL ioconf_setatt_p('UNITS', 'Kg')
    CALL ioconf_setatt_p('LONG_NAME','Accumulated drainage for routing')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., drainage_mean, "gather", nbp_glo, index_g)
    CALL setvar_p (drainage_mean, val_exp, 'NO_KEYWORD', zero)

    ALLOCATE(transpot_mean(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for transpot_mean','','')
    var_name = 'transpot_route'
    CALL ioconf_setatt_p('UNITS', 'Kg/m^2')
    CALL ioconf_setatt_p('LONG_NAME','Accumulated potential transpiration for routing/irrigation')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., transpot_mean, "gather", nbp_glo, index_g)
    CALL setvar_p (transpot_mean, val_exp, 'NO_KEYWORD', zero)

    ALLOCATE(precip_mean(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for precip_mean','','')
    var_name = 'precip_route'
    CALL ioconf_setatt_p('UNITS', 'Kg/m^2')
    CALL ioconf_setatt_p('LONG_NAME','Accumulated rain precipitation for irrigation')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., precip_mean, "gather", nbp_glo, index_g)
    CALL setvar_p (precip_mean, val_exp, 'NO_KEYWORD', zero)

    ALLOCATE(humrel_mean(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for humrel_mean','','')
    var_name = 'humrel_route'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Mean humrel for irrigation')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., humrel_mean, "gather", nbp_glo, index_g)
    CALL setvar_p (humrel_mean, val_exp, 'NO_KEYWORD', un)

    ALLOCATE(k_litt_mean(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for k_litt_mean','','')
    var_name = 'k_litt_route'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Mean cond. for litter')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., k_litt_mean, "gather", nbp_glo, index_g)
    CALL setvar_p (k_litt_mean, val_exp, 'NO_KEYWORD', zero)

    ALLOCATE(totnobio_mean(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for totnobio_mean','','')
    var_name = 'totnobio_route'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Last Total fraction of no bio for irrigation')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., totnobio_mean, "gather", nbp_glo, index_g)
    CALL setvar_p (totnobio_mean, val_exp, 'NO_KEYWORD', zero)

    ALLOCATE(vegtot_mean(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for vegtot_mean','','')
    var_name = 'vegtot_route'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Last Total fraction of vegetation')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., vegtot_mean, "gather", nbp_glo, index_g)
    CALL setvar_p (vegtot_mean, val_exp, 'NO_KEYWORD', un)
    !
    !
    DEALLOCATE(tmp_real_g)
    !
    ! Allocate diagnostic variables
    !
    ALLOCATE(hydrodiag_loc(nbpt,nbasmax),hydrodiag_glo(nbp_glo,nbasmax),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for hydrodiag_glo','','')
    hydrodiag=>hydrodiag_loc

    ALLOCATE(hydroupbasin_loc(nbpt),hydroupbasin_glo(nbp_glo), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_v0_init_1','Pb in allocate for hydroupbasin_glo','','')
    hydroupbasin=>hydroupbasin_loc

  END SUBROUTINE routing_interp_topo_v0_init_1


  !! ================================================================================================================================
  !! SUBROUTINE 	: routing_interp_topo_v0_clear
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

  SUBROUTINE routing_interp_topo_v0_clear()

    IF (is_omp_root) THEN
       IF (ALLOCATED(topoind_r)) DEALLOCATE(topoind_r)
       IF (ALLOCATED(route_flow_rp1)) DEALLOCATE(route_flow_rp1)
       IF (ALLOCATED(fast_reservoir_r)) DEALLOCATE(fast_reservoir_r)
       IF (ALLOCATED(slow_reservoir_r)) DEALLOCATE(slow_reservoir_r)
       IF (ALLOCATED(is_lakeinflow_r)) DEALLOCATE(is_lakeinflow_r) 
       IF (ALLOCATED(is_coastalflow_r)) DEALLOCATE(is_coastalflow_r)    
       IF (ALLOCATED(is_riverflow_r)) DEALLOCATE(is_riverflow_r)   
    ENDIF


    IF (ALLOCATED(routing_area_loc)) DEALLOCATE(routing_area_loc)
    IF (ALLOCATED(route_togrid_loc)) DEALLOCATE(route_togrid_loc)
    IF (ALLOCATED(route_tobasin_loc)) DEALLOCATE(route_tobasin_loc)
    IF (ALLOCATED(route_nbintobas_loc)) DEALLOCATE(route_nbintobas_loc)
    IF (ALLOCATED(global_basinid_loc)) DEALLOCATE(global_basinid_loc)
    IF (ALLOCATED(topo_resid_loc)) DEALLOCATE(topo_resid_loc)
    IF (ALLOCATED(routing_area_glo)) DEALLOCATE(routing_area_glo)
    IF (ALLOCATED(route_togrid_glo)) DEALLOCATE(route_togrid_glo)
    IF (ALLOCATED(route_tobasin_glo)) DEALLOCATE(route_tobasin_glo)
    IF (ALLOCATED(route_nbintobas_glo)) DEALLOCATE(route_nbintobas_glo)
    IF (ALLOCATED(global_basinid_glo)) DEALLOCATE(global_basinid_glo)
    IF (ALLOCATED(topo_resid_glo)) DEALLOCATE(topo_resid_glo)
    IF (ALLOCATED(fast_reservoir)) DEALLOCATE(fast_reservoir)
    IF (ALLOCATED(slow_reservoir)) DEALLOCATE(slow_reservoir)
    IF (ALLOCATED(stream_reservoir)) DEALLOCATE(stream_reservoir)
    IF (ALLOCATED(flood_reservoir)) DEALLOCATE(flood_reservoir)
    IF (ALLOCATED(flood_frac_bas)) DEALLOCATE(flood_frac_bas)
    IF (ALLOCATED(flood_height)) DEALLOCATE(flood_height)
    IF (ALLOCATED(pond_frac)) DEALLOCATE(pond_frac)
    IF (ALLOCATED(lake_reservoir)) DEALLOCATE(lake_reservoir)
    IF (ALLOCATED(pond_reservoir)) DEALLOCATE(pond_reservoir)
    IF (ALLOCATED(returnflow_mean)) DEALLOCATE(returnflow_mean)
    IF (ALLOCATED(reinfiltration_mean)) DEALLOCATE(reinfiltration_mean)
    IF (ALLOCATED(riverflow_mean)) DEALLOCATE(riverflow_mean)
    IF (ALLOCATED(coastalflow_mean)) DEALLOCATE(coastalflow_mean)
    IF (ALLOCATED(lakeinflow_mean)) DEALLOCATE(lakeinflow_mean)
    IF (ALLOCATED(runoff_mean)) DEALLOCATE(runoff_mean)
    IF (ALLOCATED(floodout_mean)) DEALLOCATE(floodout_mean)
    IF (ALLOCATED(drainage_mean)) DEALLOCATE(drainage_mean)
    IF (ALLOCATED(transpot_mean)) DEALLOCATE(transpot_mean)
    IF (ALLOCATED(precip_mean)) DEALLOCATE(precip_mean)
    IF (ALLOCATED(humrel_mean)) DEALLOCATE(humrel_mean)
    IF (ALLOCATED(k_litt_mean)) DEALLOCATE(k_litt_mean)
    IF (ALLOCATED(totnobio_mean)) DEALLOCATE(totnobio_mean)
    IF (ALLOCATED(vegtot_mean)) DEALLOCATE(vegtot_mean)
    IF (ALLOCATED(floodtemp)) DEALLOCATE(floodtemp)
    IF (ALLOCATED(hydrodiag_loc)) DEALLOCATE(hydrodiag_loc)
    IF (ALLOCATED(hydrodiag_glo)) DEALLOCATE(hydrodiag_glo)
    IF (ALLOCATED(hydroupbasin_loc)) DEALLOCATE(hydroupbasin_loc)    
    IF (ALLOCATED(hydroupbasin_glo)) DEALLOCATE(hydroupbasin_glo)
    IF (ALLOCATED(hydrographs)) DEALLOCATE(hydrographs)
    IF (ALLOCATED(slowflow_diag)) DEALLOCATE(slowflow_diag)
    IF (ALLOCATED(irrigation_mean)) DEALLOCATE(irrigation_mean)
    IF (ALLOCATED(irrigated)) DEALLOCATE(irrigated)
    IF (ALLOCATED(floodplains)) DEALLOCATE(floodplains)
    IF (ALLOCATED(swamp)) DEALLOCATE(swamp)
    IF (ALLOCATED(fast_diag)) DEALLOCATE(fast_diag)
    IF (ALLOCATED(slow_diag)) DEALLOCATE(slow_diag)
    IF (ALLOCATED(stream_diag)) DEALLOCATE(stream_diag)
    IF (ALLOCATED(flood_diag)) DEALLOCATE(flood_diag)
    IF (ALLOCATED(pond_diag)) DEALLOCATE(pond_diag)
    IF (ALLOCATED(lake_diag)) DEALLOCATE(lake_diag)
    IF (ALLOCATED(delsurfstor)) DEALLOCATE(delsurfstor)
    !
  END SUBROUTINE routing_interp_topo_v0_clear
  !

  !! ================================================================================================================================
  !! SUBROUTINE 	: routing_interp_topo_v0_lake
  !!
  !>\BRIEF        : This subroutine stores water in lakes so that it does not cycle through the runoff.
  !!                For the moment it only works for endorheic lakes but it can be extended in the future.
  !!
  !! DESCRIPTION (definitions, functional, design, flags): The return flow to the soil moisture reservoir
  !! is based on a maximum lake evaporation rate (maxevap_lake). \n
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

  SUBROUTINE routing_interp_topo_v0_lake(nbpt, dt_routing, lakeinflow, humrel, contfrac, return_lakes)

    USE grid, ONLY : area
    
    IMPLICIT NONE
    !! 0 Variable and parameter description
    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in) :: nbpt               !! Domain size (unitless)
    REAL(r_std), INTENT (in)   :: dt_routing         !! Routing time step (s)
    REAL(r_std), INTENT(in)    :: lakeinflow(nbpt)   !! Water inflow to the lakes (kg/dt)
    REAL(r_std), INTENT(in)    :: humrel(nbpt)       !! Soil moisture stress, root extraction potential (unitless)
    REAL(r_std), INTENT(in)    :: contfrac(nbpt)     !! Fraction of land in each grid box (unitless;0-1)
    
    !! 0.2 Output variables
    REAL(r_std), INTENT(out)   :: return_lakes(nbpt) !! Water from lakes flowing back into soil moisture (kg/m^2/dt)
    
    !! 0.3 Local variables 
    INTEGER(i_std)             :: ig                 !! Indices (unitless)
    REAL(r_std)                :: refill             !!
    REAL(r_std)                :: total_area         !! Sum of all the surfaces of the basins (m^2)

    !_ ================================================================================================================================

    
    DO ig=1,nbpt
       !
       total_area = area(ig)*contfrac(ig)
       !
       lake_reservoir(ig) = lake_reservoir(ig) + lakeinflow(ig)

       IF ( doswamps ) THEN
          ! uptake in Kg/dt
          refill = MAX(zero, maxevap_lake * (un - humrel(ig)) * dt_routing * total_area)
          return_lakes(ig) = MIN(refill, lake_reservoir(ig))
          lake_reservoir(ig) = lake_reservoir(ig) - return_lakes(ig)
          ! Return in Kg/m^2/dt
          return_lakes(ig) = return_lakes(ig)/total_area
       ELSE
          return_lakes(ig) = zero
       ENDIF
       !
       ! This is the volume of the lake scaled to the entire grid.
       ! It would be batter to scale it to the size of the lake
       ! but this information is not yet available.
       lake_diag(ig) = lake_reservoir(ig)/total_area
       !
    ENDDO

  END SUBROUTINE routing_interp_topo_v0_lake
#endif
END MODULE routing_interp_topo_v0
