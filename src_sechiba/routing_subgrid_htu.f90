! ================================================================================================================================
! MODULE       : routing_subgrid_htu
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       This module routes the water over the continents into the oceans and computes the water
!!             stored in floodplains or taken for irrigation.
!!
!!\n DESCRIPTION: This method computes the routing on the ORCHIDEE grid using high resolution subgrid-basins
!!               (Hydrological Transfer Units HTUs).
!!
!! RECENT CHANGE(S): July 2024: Changed routing method's name "highres" to "subgrid_htu".
!!                   Now works together with the routing pre-processor : https://gitlab.in2p3.fr/ipsl/lmd/intro/routingpp
!!
!! REFERENCE(S)	:
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_sechiba/routing_subgrid_htu.f90 $
!! $Date: 2026-04-19 16:35:39 +0200 (dim. 19 avril 2026) $
!! $Revision: 9494 $
!! \n
!_ ================================================================================================================================
! 
!
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


MODULE routing_subgrid_htu
  
  USE ioipsl   
  USE xios_orchidee
  USE ioipsl_para 
  USE constantes
  USE constantes_var
  USE constantes_soil
  USE pft_parameters
  USE sechiba_io_p
  USE interpol_help
  USE grid, not_used => contfrac
  USE mod_orchidee_para

  USE haversine

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: routing_subgrid_htu_main, routing_subgrid_htu_initialize, routing_subgrid_htu_finalize, routing_subgrid_htu_clear, routing_subgrid_htu_xios_initialize

  INTERFACE routing_hr_landgather
     MODULE PROCEDURE routing_hr_landgather_i1, routing_hr_landgather_i2, routing_hr_landgather_r
  END INTERFACE routing_hr_landgather

  INTEGER(i_std),PARAMETER                                   :: WaterCp=1000.*4.1813        !! water heat capacity in J/Kg/K
  
!! PARAMETERS
  INTEGER(i_std), SAVE                                       :: nbasmax=-1                  !! The maximum number of basins we wish to have per grid box (truncation of the model) (unitless)
  INTEGER(i_std), SAVE                                       :: nbasmon = 4                 !! Number of basins to be monitored
  INTEGER(i_std), SAVE                                       :: inflows=-1                  !! The maximum number of inflows (unitless)
  INTEGER(i_std), SAVE                                       :: nbvmax                      !! The maximum number of basins we can handle at any time during the generation of the maps (unitless)
!$OMP THREADPRIVATE(nbvmax)
  REAL(r_std), SAVE                                          :: fast_tcst = -1.             !! Property of the fast reservoir, (s/km)
!$OMP THREADPRIVATE(fast_tcst)
  REAL(r_std), SAVE                                          :: slow_tcst = -1.             !! Property of the slow reservoir, (s/km)
!$OMP THREADPRIVATE(slow_tcst)
  REAL(r_std), SAVE                                          :: stream_tcst = -1.           !! Property of the stream reservoir, (s/km)
!$OMP THREADPRIVATE(stream_tcst)
  REAL(r_std), SAVE                                          :: flood_tcst = -1.            !! Property of the floodplains reservoir, (s/km)
!$OMP THREADPRIVATE(flood_tcst)
  REAL(r_std), SAVE                                          :: swamp_cst = -1.             !! Fraction of the river transport that flows to the swamps (unitless;0-1)
!$OMP THREADPRIVATE(swamp_cst)
  REAL(r_std), SAVE                                          :: lim_floodcri = -1.          !! Minimal orog diff between two consecutive floodplains htu (m)
!$OMP THREADPRIVATE(lim_floodcri)
  !
  !  Relation between volume and fraction of floodplains
  !
  REAL(r_std), SAVE                                          :: betap = 0.5                 !! Ratio of the basin surface intercepted by ponds and the maximum surface of ponds (unitless;0-1)
!$OMP THREADPRIVATE(betap)
  REAL(r_std), SAVE                                          :: rfloodmax = 0.5             !! Maximal discharge reducer when there are floodplains
!$OMP THREADPRIVATE(rfloodmax)
  REAL(r_std), SAVE                                          :: overflow_tcst = 5           !! Maximal discharge reducer when there are floodplains
  !$OMP THREADPRIVATE(overflow_tcst)
  INTEGER(i_std), SAVE                                       :: overflow_repetition = 1     !! Number of repetition of overflow for each routing step
  !$OMP THREADPRIVATE(overflow_repetition)
  ! Soil temperature depth to be used to estimate runoff and drainage temperatures
  !
  REAL(r_std), PARAMETER, DIMENSION(2) :: runofftempdepth = (/ 0.0, 0.3 /)                  !! Layer which will determine the temperature of runoff
  REAL(r_std), PARAMETER, DIMENSION(2) :: drainagetempdepth = (/ 3.0, 90.0 /)               !! Layer which will determine the temperature of runoff
  !
  !
  !  Relation between maximum surface of ponds and basin surface, and drainage (mm/j) to the slow_res
  !
  REAL(r_std), PARAMETER                                     :: pond_bas = 50.0             !! [DISPENSABLE] - not used
  REAL(r_std), SAVE                                          :: pondcri = 2000.0            !! Potential height for which all the basin is a pond (mm)
!$OMP THREADPRIVATE(pondcri)
  !
  REAL(r_std), PARAMETER                                     :: maxevap_lake = 7.5/86400.   !! Maximum evaporation rate from lakes (kg/m^2/s)
  !
  REAL(r_std),SAVE                                           :: dt_routing                  !! Routing time step (s)
!$OMP THREADPRIVATE(dt_routing)
  !
  INTEGER(i_std), SAVE                                       :: ntemp_layer = 4             !! Number of layers to be taken to determine the ground water temperature.
!$OMP THREADPRIVATE(ntemp_layer)
  INTEGER(i_std), SAVE                                       :: diagunit = 87               !! Diagnostic file unit (unitless)
!$OMP THREADPRIVATE(diagunit)
  !
  ! Logicals to control model configuration
  !
  LOGICAL, SAVE                                              :: dofloodinfilt = .FALSE.     !! Logical to choose if floodplains infiltration is activated or not (true/false)
!$OMP THREADPRIVATE(dofloodinfilt)
  LOGICAL, SAVE                                              :: dofloodoverflow = .FALSE.   !! Logical to choose if floodplains overflow is activated or not (true/false)
!$OMP THREADPRIVATE(dofloodoverflow)
  LOGICAL, SAVE                                              :: doswamps = .FALSE.          !! Logical to choose if swamps are activated or not (true/false)
!$OMP THREADPRIVATE(doswamps)
  LOGICAL, SAVE                                              :: doponds = .FALSE.           !! Logical to choose if ponds are activated or not (true/false)
!$OMP THREADPRIVATE(doponds)
  REAL(r_std), SAVE                                          :: conduct_factor = 1.         !! Adjustment factor for floodplains reinfiltration
!$OMP THREADPRIVATE(conduct_factor)
  !
  ! The variables describing the basins and their routing, need to be in the restart file.
  !
  INTEGER(i_std), SAVE                                       :: num_largest = 200           !! Number of largest river basins which should be treated as independently as rivers
  !! (not flow into ocean as diffusion coastal flow) (unitless)
!$OMP THREADPRIVATE(num_largest)
  !
  CHARACTER(LEN=80),SAVE                                     :: graphfilename="routing_graph.nc"
!$OMP THREADPRIVATE(graphfilename)
  REAL(r_std), SAVE                                          :: undef_graphfile
!$OMP THREADPRIVATE(undef_graphfile)
  REAL(r_std), SAVE                                          :: graphfile_version = 0.0
!$OMP THREADPRIVATE(graphfile_version)
  REAL(r_std), SAVE                                          :: maxtimestep = 1800.0        !! A reasonalble maximum time step. Actual value to be read from graphfile.
!$OMP THREADPRIVATE(maxtimestep)
  REAL(r_std), SAVE                                          :: time_counter                !! Time counter (s)
!$OMP THREADPRIVATE(time_counter)
  REAL(r_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:,:)     :: routing_area_loc            !! Surface of basin (m^2)
!$OMP THREADPRIVATE(routing_area_loc)
  REAL(r_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:,:)     :: topo_resid_loc              !! Topographic index of the retention time (m)
!$OMP THREADPRIVATE(topo_resid_loc)
  REAL(r_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:,:)     :: stream_resid_loc              !! Topographic index of the retention time (m)
!$OMP THREADPRIVATE(stream_resid_loc)
  INTEGER(i_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:,:)  :: route_togrid_loc            !! Grid into which the basin flows (unitless)
!$OMP THREADPRIVATE(route_togrid_loc)
  INTEGER(i_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:,:)  :: route_tobasin_loc           !! Basin in to which the water goes (unitless)
!$OMP THREADPRIVATE(route_tobasin_loc)
  INTEGER(i_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:,:)  :: route_nbintobas_loc         !! Number of basin into current one (unitless)
!$OMP THREADPRIVATE(route_nbintobas_loc)
  INTEGER(i_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:,:)  :: global_basinid_loc          !! ID of basin (unitless)
!$OMP THREADPRIVATE(global_basinid_loc)
  INTEGER(i_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:)    :: hydrodiag_loc               !! Variable to diagnose the hydrographs
!$OMP THREADPRIVATE(hydrodiag_loc)
  INTEGER(i_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:,:)  :: HTUdiag_loc                 !! Variable to diagnose the hydrographs
!$OMP THREADPRIVATE(HTUdiag_loc)
  INTEGER(i_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:,:)  :: HTUdiag_glo                 !! Variable to diagnose the hydrographs
!$OMP THREADPRIVATE(HTUdiag_glo)
  LOGICAL, SAVE                                              :: MonitoringinGraph=.FALSE.
  LOGICAL, SAVE                                              :: ReadGraph=.FALSE.
  LOGICAL, SAVE                                              :: ReadMonitoring=.FALSE.
  REAL(r_std), SAVE                                          :: stream_maxresid
!$OMP THREADPRIVATE(stream_maxresid)
  !
  ! parallelism
  REAL(r_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:,:)     :: routing_area_glo            !! Surface of basin (m^2)
!$OMP THREADPRIVATE(routing_area_glo)
  REAL(r_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:,:)     :: topo_resid_glo              !! Topographic index of the retention time (m)
!$OMP THREADPRIVATE(topo_resid_glo)
  REAL(r_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:,:)     :: stream_resid_glo            !! Topographic index of the retention time (m)
!$OMP THREADPRIVATE(stream_resid_glo)
  INTEGER(i_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:,:)  :: route_togrid_glo            !! Grid into which the basin flows (unitless)
!$OMP THREADPRIVATE(route_togrid_glo)
  INTEGER(i_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:,:)  :: route_tobasin_glo           !! Basin in to which the water goes (unitless)
!$OMP THREADPRIVATE(route_tobasin_glo)
  INTEGER(i_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:,:)  :: route_nbintobas_glo         !! Number of basin into current one (unitless)
!$OMP THREADPRIVATE(route_nbintobas_glo)
  INTEGER(i_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:,:)  :: global_basinid_glo          !! ID of basin (unitless)
!$OMP THREADPRIVATE(global_basinid_glo)
  INTEGER(i_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:)    :: hydrodiag_glo               !! Variable to diagnose the hydrographs
!$OMP THREADPRIVATE(hydrodiag_glo)
  !
  REAL(r_std), SAVE, POINTER, DIMENSION(:,:)                 :: routing_area                !! Surface of basin (m^2)
!$OMP THREADPRIVATE(routing_area)
  REAL(r_std), SAVE, POINTER, DIMENSION(:,:)                 :: topo_resid                  !! Topographic index of the retention time (m)
!$OMP THREADPRIVATE(topo_resid)
  REAL(r_std), SAVE, POINTER, DIMENSION(:,:)                 :: stream_resid                  !! Topographic index of the retention time (m)
!$OMP THREADPRIVATE(stream_resid)
  INTEGER(i_std), SAVE, POINTER, DIMENSION(:,:)              :: route_togrid                !! Grid into which the basin flows (unitless)
!$OMP THREADPRIVATE(route_togrid)
  INTEGER(i_std), SAVE, POINTER, DIMENSION(:,:)              :: route_tobasin               !! Basin in to which the water goes (unitless)
!$OMP THREADPRIVATE(route_tobasin)
  INTEGER(i_std), SAVE, POINTER, DIMENSION(:,:)              :: route_nbintobas             !! Number of basin into current one (unitless)
!$OMP THREADPRIVATE(route_nbintobas)
  INTEGER(i_std), SAVE, POINTER, DIMENSION(:,:)              :: global_basinid              !! ID of basin (unitless)
!$OMP THREADPRIVATE(global_basinid)
  INTEGER(i_std), SAVE, POINTER, DIMENSION(:)                :: hydrodiag                   !! Variable to diagnose the hydrographs
!$OMP THREADPRIVATE(hydrodiag)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: slowflow_diag               !! Diagnostic slow flow hydrographs (kg/dt)
!$OMP THREADPRIVATE(slowflow_diag)  
  !
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: irrigated                   !! Area equipped for irrigation in each grid box (m^2)
!$OMP THREADPRIVATE(irrigated)
  REAL(r_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:,:)     :: floodplains_glo             !! Maximal surface which can be inundated in each grid box (m^2)
!$OMP THREADPRIVATE(floodplains_glo)
  REAL(r_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:,:)     :: floodplains_loc             !! Maximal surface which can be inundated in each grid box (m^2)
!$OMP THREADPRIVATE(floodplains_loc)
  REAL(r_std), SAVE, POINTER, DIMENSION(:,:)                 :: floodplains                 !! Maximal surface which can be inundated in each grid box (m^2)
!$OMP THREADPRIVATE(floodplains)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: floodmap
!! Floodplains Fraction for each grid point.
!$OMP THREADPRIVATE(floodmap)

!!!

  REAL(r_std),  SAVE, ALLOCATABLE, DIMENSION(:,:)            :: tempdiag_mean              !! Averaged soil temperatures
!$OMP THREADPRIVATE(tempdiag_mean)
!
! FLOOD OVERFLOW
  REAL(r_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:,:)  :: orog_min_glo !!            
!$OMP THREADPRIVATE(orog_min_glo)!
  REAL(r_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:,:)  :: orog_min_loc             !! 
!$OMP THREADPRIVATE(orog_min_loc)
  REAL(r_std), SAVE, POINTER, DIMENSION(:,:)              :: orog_min                 !!
!$OMP THREADPRIVATE(orog_min)
  !
  INTEGER(i_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:,:)  :: route_innum_glo             !! 
!$OMP THREADPRIVATE(route_innum_glo)
  INTEGER(i_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:,:)  :: route_innum_loc             !! 
!$OMP THREADPRIVATE(route_innum_loc)
  INTEGER(i_std), SAVE, POINTER, DIMENSION(:,:)              :: route_innum                 !! 
!$OMP THREADPRIVATE(route_innum)
  !
  INTEGER(i_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: route_ingrid_glo             !! 
!$OMP THREADPRIVATE(route_ingrid_glo)
  INTEGER(i_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: route_ingrid_loc             !! 
!$OMP THREADPRIVATE(route_ingrid_loc)
  INTEGER(i_std), SAVE, POINTER, DIMENSION(:,:,:)             :: route_ingrid                 !! 
!$OMP THREADPRIVATE(route_ingrid)
  !
  INTEGER(i_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: route_inbasin_glo             !! 
!$OMP THREADPRIVATE(route_inbasin_glo)
  INTEGER(i_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: route_inbasin_loc             !! 
!$OMP THREADPRIVATE(route_inbasin_loc)
  INTEGER(i_std), SAVE, POINTER, DIMENSION(:,:,:)             :: route_inbasin                 !! 
!$OMP THREADPRIVATE(route_inbasin)
  !  
!!!
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: swamp                       !! Maximal surface of swamps in each grid box (m^2)
!$OMP THREADPRIVATE(swamp)
  REAL(r_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:,:)     :: fp_beta_glo                 !! Parameter to fix the shape of the floodplain (>1 for convex edges, <1 for concave edges) (unitless)
!$OMP THREADPRIVATE(fp_beta_glo)
  REAL(r_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:,:)     :: fp_beta_loc                 !! Parameter to fix the shape of the floodplain (>1 for convex edges, <1 for concave edges) (unitless)
!$OMP THREADPRIVATE(fp_beta_loc)
  REAL(r_std), SAVE, POINTER, DIMENSION(:,:)                 :: fp_beta                     !! Parameter to fix the shape of the floodplain (>1 for convex edges, <1 for concave edges) (unitless)
!$OMP THREADPRIVATE(fp_beta)
  REAL(r_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:,:)     :: floodcri_glo                !! Potential height for which all the basin is a pond (mm)
!$OMP THREADPRIVATE(floodcri_glo)
  REAL(r_std), SAVE, ALLOCATABLE, TARGET, DIMENSION(:,:)     :: floodcri_loc                !! Potential height for which all the basin is a pond (mm)
!$OMP THREADPRIVATE(floodcri_loc)
  REAL(r_std), SAVE, POINTER, DIMENSION(:,:)                 :: floodcri                    !! Potential height for which all the basin is a pond (mm)
!$OMP THREADPRIVATE(floodcri)
  !
  ! The reservoirs, also to be put into the restart file.
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
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)             :: flood_height                !! Floodplain height (mm)
!$OMP THREADPRIVATE(flood_height)
  !
  ! Reservoir temperatures
  !
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)             :: fast_temp                   !! Water temperature in the fast reservoir (K)
!$OMP THREADPRIVATE(fast_temp)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)             :: slow_temp                   !! Water temperature in the slow reservoir (K)
!$OMP THREADPRIVATE(slow_temp)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)             :: stream_temp                 !! Water temperature in the stream reservoir (K)
!$OMP THREADPRIVATE(stream_temp)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: streamlimit                 !!
!$OMP THREADPRIVATE(streamlimit)
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
                                                                                             !! This water will go back into the hydrol module to allow re-evaporation (kg/m^2/dt)
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
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: hydrotemp                    !! Temperature of the largest river (in the HTUdiag sense) in the grid (K)
!$OMP THREADPRIVATE(hydrotemp)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)             :: HTUhgmon                     !! Hydrographs to be monitored on specific HTUs (kg/dt)
!$OMP THREADPRIVATE(HTUhgmon)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)             :: HTUhgmon_glo                 !! Hydrographs to be monitored on specific HTUs (kg/dt)
!$OMP THREADPRIVATE(HTUhgmon_glo)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)             :: HTUtempmon                   !! Temperature to be monitored on specific HTUs (K)
!$OMP THREADPRIVATE(HTUtempmon)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)             :: HTUtempmon_glo                !! Temperature to be monitored on specific HTUs (K)
!$OMP THREADPRIVATE(HTUtempmon_glo)
  !
  ! Diagnostics for the various reservoirs we use (Kg/m^2)
  !
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: fast_diag                    !! Diagnostic for the fast reservoir (kg/m^2)
!$OMP THREADPRIVATE(fast_diag)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: slow_diag                    !! Diagnostic for the slow reservoir (kg/m^2)
!$OMP THREADPRIVATE(slow_diag)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: stream_diag                  !! Diagnostic for the stream reservoir (kg/m^2)
!$OMP THREADPRIVATE(stream_diag)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: flood_diag                   !! Diagnostic for the floodplain reservoir over grid (kg/m^2)
!$OMP THREADPRIVATE(flood_diag)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: flood_height_diag            !! Diagnostic for the floodplain height over grid (kg/m^2)
!$OMP THREADPRIVATE(flood_height_diag)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: pond_diag                    !! Diagnostic for the pond reservoir (kg/m^2)
!$OMP THREADPRIVATE(pond_diag)
  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: lake_diag                    !! Diagnostic for the lake reservoir (kg/m^2)
!$OMP THREADPRIVATE(lake_diag)

  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: mask_coast                   !! Mask with coastal gridcells on local grid(1/0)
!$OMP THREADPRIVATE(mask_coast)
  REAL(r_std), SAVE                                          :: max_lake_reservoir           !! Maximum limit of water in lake_reservoir [kg/m2]
  !$OMP THREADPRIVATE(max_lake_reservoir)
  INTEGER(i_std), SAVE                                       :: nb_coast_gridcells           !! Number of gridcells which can receive coastalflow
!$OMP THREADPRIVATE(nb_coast_gridcells)


CONTAINS
  !!  =============================================================================================================================
  !! SUBROUTINE:         routing_subgrid_htu_initialize
  !!
  !>\BRIEF	         Initialize the routing module
  !!
  !! DESCRIPTION:        Initialize the routing module. Read from restart file or read the routing.nc file to initialize the
  !!                     routing scheme. 
  !!
  !! RECENT CHANGE(S)
  !!
  !! REFERENCE(S)
  !! 
  !! FLOWCHART   
  !! \n
  !_ ==============================================================================================================================

  SUBROUTINE routing_subgrid_htu_initialize( kjit,        nbpt,           index,                &
                                             rest_id,     hist_id,        hist2_id,   lalo,     &
                                             neighbours,  resolution,     contfrac,   tempdiag, &
                                             returnflow,  reinfiltration, irrigation, riverflow,&
                                             coastalflow, flood_frac,     flood_res )
       
    IMPLICIT NONE
    
    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)     :: kjit                 !! Time step number (unitless)
    INTEGER(i_std), INTENT(in)     :: nbpt                 !! Domain size (unitless)
    INTEGER(i_std), INTENT(in)     :: index(nbpt)          !! Indices of the points on the map (unitless)
    INTEGER(i_std),INTENT(in)      :: rest_id              !! Restart file identifier (unitless)
    INTEGER(i_std),INTENT(in)      :: hist_id              !! Access to history file (unitless)
    INTEGER(i_std),INTENT(in)      :: hist2_id             !! Access to history file 2 (unitless)
    REAL(r_std), INTENT(in)        :: lalo(nbpt,2)         !! Vector of latitude and longitudes (beware of the order !)

    INTEGER(i_std), INTENT(in)     :: neighbours(nbpt,NbNeighb) !! Vector of neighbours for each grid point 
                                                           !! (1=N, 2=NE, 3=E, 4=SE, 5=S, 6=SW, 7=W, 8=NW) (unitless)
    REAL(r_std), INTENT(in)        :: resolution(nbpt,2)   !! The size of each grid box in X and Y (m)
    REAL(r_std), INTENT(in)        :: contfrac(nbpt)       !! Fraction of land in each grid box (unitless;0-1)
    REAL(r_std), INTENT(in)        :: tempdiag(nbpt,ngrnd) !! Diagnostic soil temperature profile

    !! 0.2 Output variables
    REAL(r_std), INTENT(out)       :: returnflow(nbpt)     !! The water flow from lakes and swamps which returns to the grid box.
                                                           !! This water will go back into the hydrol module to allow re-evaporation (kg/m^2/dt)
    REAL(r_std), INTENT(out)       :: reinfiltration(nbpt) !! Water flow from ponds and floodplains which returns to the grid box (kg/m^2/dt)
    REAL(r_std), INTENT(out)       :: irrigation(nbpt)     !! Irrigation flux. This is the water taken from the reservoirs and beeing put into the upper layers of the soil (kg/m^2/dt)
    REAL(r_std), INTENT(out)       :: riverflow(nbpt)      !! Outflow of the major rivers. The flux will be located on the continental grid but this should be a coastal point (kg/dt)

    REAL(r_std), INTENT(out)       :: coastalflow(nbpt)    !! Outflow on coastal points by small basins. This is the water which flows in a disperse way into the ocean (kg/dt)
    REAL(r_std), INTENT(out)       :: flood_frac(nbpt)     !! Flooded fraction of the grid box (unitless;0-1)
    REAL(r_std), INTENT(out)       :: flood_res(nbpt)      !! Diagnostic of water amount in the floodplains reservoir (kg)
    
    !! 0.3 Local variables
    REAL(r_std), DIMENSION(nbp_glo)        :: mask_coast_glo       !! Mask with coastal gridcells on global grid (1/0)
    LOGICAL                                :: init_irrig           !! Logical to initialize the irrigation (true/false)
    LOGICAL                                :: init_flood           !! Logical to initialize the floodplains (true/false)
    LOGICAL                                :: init_swamp           !! Logical to initialize the swamps (true/false)
    INTEGER                                :: ig, ib, rtg, rtb     !! Index
    REAL(r_std)                            :: stream_tcst_orig, totarea
    INTEGER                                :: ier                  !! Error handeling
!_ ================================================================================================================================

    IF ( do_floodplains ) THEN
       !
       ! Force ReadGraph when do_floodplains
       !
       ReadGraph = .TRUE.
    ENDIF
    !
    ! do initialisation
    !
    nbvmax = 440
    ! Here we will allocate the memory and get the fixed fields from the restart file.
    ! If the info is not found then we will compute the routing map.
    !

    CALL routing_hr_init (kjit, nbpt, index, returnflow, reinfiltration, irrigation, &
         riverflow, coastalflow, flood_frac, flood_res, tempdiag, rest_id)
    
    ! This routine computes the routing map if the route_togrid_glo is undefined. This means that the
    ! map has not been initialized during the restart process..
    !
    !! Reads in the map of the basins and flow directions to construct the catchments of each grid box
    !
    IF ( ReadGraph .OR. ReadMonitoring) THEN
       CALL routing_hr_basins_p(nbpt, lalo, neighbours, resolution, contfrac)
    ENDIF
    
    ! Keep the information so we can check the time step.
    stream_tcst_orig = stream_tcst
    !
    IF (stream_tcst .LE. 0 .OR. fast_tcst .LE. 0 .OR. slow_tcst .LE. 0 .OR. flood_tcst .LE. 0 ) THEN
       CALL ipslerr(3,'routing_subgrid_htu_initialize',' The time constants of the routing reservoirs were not initialized. ', &
            'Please check if they are present in the HTU graph file', ' ')
    ELSE
       !
!> The time constants for the various reservoirs should be read from the graph file
!> produced by routingpp (https://gitlab.in2p3.fr/ipsl/lmd/intro/routingpp). They are
!> also saved in the restart file so that we do not need to read the graph file at each restart.
!> But once they are set in the model the user can changed them through the run.def.
!> This is a useful option to test values but should not be an operational solution. The
!> correct value should be given to the model through the graph file.
!> The getin_p operation cannot be done earlier as in routing_hr_init above these constant
!> might not yet be known.
       !
       !Config Key   = SLOW_TCST
       !Config Desc  = Time constant for the slow reservoir 
       !Config If    = RIVER_ROUTING 
       !Config Def   = 25.0 
       !Config Help  = This parameters allows the user to fix the 
       !Config         time constant (s/km) of the slow reservoir
       !Config         in order to get better river flows for 
       !Config         particular regions.
       !Config Units = [days]
       !
       CALL getin_p('SLOW_TCST', slow_tcst)
       !
       !Config Key   = FAST_TCST
       !Config Desc  = Time constant for the fast reservoir 
       !Config If    = RIVER_ROUTING 
       !Config Def   = 3.0 
       !Config Help  = This parameters allows the user to fix the 
       !Config         time constant (s/km) of the fast reservoir
       !Config         in order to get better river flows for 
       !Config         particular regions.
       !Config Units = [days]
       CALL getin_p('FAST_TCST', fast_tcst)

       !Config Key   = STREAM_TCST
       !Config Desc  = Time constant for the stream reservoir 
       !Config If    = RIVER_ROUTING
       !Config Def   = 0.24
       !Config Help  = This parameters allows the user to fix the 
       !Config         time constant (s/km) of the stream reservoir
       !Config         in order to get better river flows for 
       !Config         particular regions.
       !Config Units = [days]
       CALL getin_p('STREAM_TCST', stream_tcst)

       !Config Key   = FLOOD_TCST
       !Config Desc  = Time constant for the flood reservoir 
       !Config If    = RIVER_ROUTING
       !Config Def   = 4.0
       !Config Help  = This parameters allows the user to fix the 
       !Config         time constant (s/km) of the flood reservoir
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
       !
       !Config Key   = LIM_FLOODCRI
       !Config Desc  = Difference of orography between floodplains HTUs.
       !Config If    = RIVER_ROUTING
       !Config Def   = 0.3
       !Config Help  = This parameters allows the user to fix the 
       !Config         minimal difference of orography between two consecutive
       !Config         floodplains HTU.
       !Config Units = [meter]
       CALL getin_p('LIM_FLOODCRI', lim_floodcri)
       !
    ENDIF
    !
    ! Final set of pointer associations. These variables can
    ! either come from the restart or the graph files.
    !
    routing_area => routing_area_loc  
    hydrodiag => hydrodiag_loc
    topo_resid=>topo_resid_loc
    stream_resid=>stream_resid_loc
    !
    route_togrid=>route_togrid_loc
    route_tobasin=>route_tobasin_loc
    global_basinid=>global_basinid_loc
    route_nbintobas=>route_nbintobas_loc
    !
    floodplains=>floodplains_loc
    fp_beta=>fp_beta_loc
    floodcri=>floodcri_loc
    orog_min=>orog_min_loc
    !
    route_innum=>route_innum_loc
    route_ingrid=>route_ingrid_loc
    route_inbasin=>route_inbasin_loc
    !
    ! Initialize the diagnostic variables before they are used.
    !
    fast_diag(:) = zero
    slow_diag(:) = zero
    stream_diag(:) = zero
    flood_diag(:) = zero
    flood_height_diag(:) = zero
    pond_diag(:) = zero
    lake_diag(:) = zero
    
    DO ig=1,nbpt
       totarea = zero
       DO ib=1,nbasmax
          totarea = totarea + routing_area(ig,ib)
          fast_diag(ig) = fast_diag(ig) + fast_reservoir(ig,ib)
          slow_diag(ig) = slow_diag(ig) + slow_reservoir(ig,ib)
          stream_diag(ig) = stream_diag(ig) + stream_reservoir(ig,ib)
          IF (floodplains(ig,ib) .GT. min_sechiba) THEN
             flood_diag(ig) = flood_diag(ig) + flood_reservoir(ig,ib)
             flood_height_diag(ig) = flood_height_diag(ig) + flood_height(ig,ib)*routing_area(ig,ib)
          ENDIF
       ENDDO
       !
       fast_diag(ig) = fast_diag(ig)/totarea
       slow_diag(ig) = slow_diag(ig)/totarea
       stream_diag(ig) = stream_diag(ig)/totarea
       flood_diag(ig) = flood_diag(ig)/totarea
       flood_height_diag(ig) = flood_height_diag(ig)/totarea
       !
       ! This is the volume of the lake scaled to the entire grid.
       ! It would be better to scale it to the size of the lake
       ! but this information is not yet available.
       !
       lake_diag(ig) = lake_reservoir(ig)/totarea
       !
    ENDDO
    !
    ! Verify that the time step is compatible with the graph file.
    ! If the user has changed the time constant of the stream reservoir then
    ! the maximum time step needs to be adjusted.
    !
    IF ( stream_tcst_orig == 0 ) THEN
       WRITE(*,*) "routing_subgrid_htu_initialize : Update stream_tcst ", stream_tcst_orig, stream_tcst
       stream_tcst_orig = stream_tcst
    ENDIF
    IF ( dt_routing > maxtimestep/stream_tcst_orig*stream_tcst ) THEN
       WRITE(*,*) "routing_subgrid_htu_initialize : Chosen time step : ", dt_routing
       WRITE(*,*) "routing_subgrid_htu_initialize : Recommended time step : ", maxtimestep/stream_tcst_orig*stream_tcst
       CALL ipslerr_p(2,'routing_subgrid_htu_initialize','The chosen time step is larger than the value recommended','in the graph file.','')
    ENDIF
    !
    !
    !
    IF (dofloodoverflow) THEN
      CALL routing_hr_inflows(nbp_glo, nbasmax, inflows, floodplains_glo,route_innum_glo,route_ingrid_glo,route_inbasin_glo)
    END IF 

    !! Create a mask containing all possible coastal gridcells and count total number of coastal gridcells
    IF (is_root_prc) THEN
       mask_coast_glo(:)=0
       DO ib=1,nbasmax
          DO ig=1,nbp_glo
             rtg = route_togrid_glo(ig,ib)
             rtb = route_tobasin_glo(ig,ib)
             ! Coastal gridcells are stored in nbasmax+2
             IF ( rtb == nbasmax+2) THEN
                mask_coast_glo(rtg) = 1
             END IF
          END DO
       END DO
       nb_coast_gridcells=SUM(mask_coast_glo)
       IF (printlev>=3) WRITE(numout,*) 'Number of coastal gridcells = ', nb_coast_gridcells

       IF (nb_coast_gridcells == 0)THEN
          CALL ipslerr(3,'routing_subgrid_htu_initialize',&
               'Number of coastal gridcells is zero for routing. ', &
               'If this is a global run, this is an error.',&
               'If this is a regional run, please check to make sure your region includes a full basin or turn routing off.')
       ENDIF

    ENDIF
    CALL bcast(nb_coast_gridcells)
   
    ALLOCATE(mask_coast(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_inititalize','Pb in allocate for mask_coast','','')
    CALL scatter(mask_coast_glo, mask_coast)
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
    
    !IF ( init_irrig .OR. init_flood .OR. init_swamp ) THEN
    !   CALL routing_hr_irrigmap(nbpt, index, lalo, neighbours, resolution, &
    !        contfrac, init_irrig, irrigated, init_flood, floodplains, init_swamp, swamp, hist_id, hist2_id)
    !ENDIF
    
    IF (printlev >= 5) WRITE(numout,*) 'End of routing_subgrid_htu_initialize'

  END SUBROUTINE routing_subgrid_htu_initialize

  !!  =============================================================================================================================
  !! SUBROUTINE:    routing_subgrid_htu_xios_initialize
  !!
  !>\BRIEF          Initialize xios dependant definition before closing context definition
  !!
  !! DESCRIPTION:   Initialize xios dependant definition before closing context definition.
  !!                This subroutine is called before the xios context is closed. 
  !!
  !! RECENT CHANGE(S): None
  !!
  !! REFERENCE(S): None
  !! 
  !! FLOWCHART: None
  !! \n
  !_ ==============================================================================================================================

  SUBROUTINE routing_subgrid_htu_xios_initialize
    USE xios_orchidee
    IMPLICIT NONE

    INTEGER(i_std) ::ib

    !
    ! If the routing_graph file is available we will extract the information in the dimensions
    ! and parameters.
    !
    !Config Key   = ROUTING_FILE
    !Config Desc  = Name of file which contains the routing information graph on the model grid
    !Config If    = RIVER_ROUTING
    !Config Def   = routing.nc
    !Config Help  = The file provided here should allows to route the water from one HTU
    !Config         to another. The RoutingPP code needs to be used in order to generate
    !Config         the routing graph for the model grid.
    !Config         More details on : https://gitlab.in2p3.fr/ipsl/lmd/intro/routingpp
    !Config Units = [FILE]
    !
    graphfilename = 'routing_graph.nc'
    CALL getin('ROUTING_FILE',graphfilename)
    CALL routing_hr_graphinfo(graphfilename, nbasmax, inflows, nbasmon, undef_graphfile, stream_tcst, fast_tcst, slow_tcst, &
         &                 flood_tcst, swamp_cst, lim_floodcri)

    CALL xios_orchidee_addaxis("nbhtu", nbasmax, (/(REAL(ib,r_std),ib=1,nbasmax)/))
    CALL xios_orchidee_addaxis("nbasmon", nbasmon, (/(REAL(ib,r_std),ib=1,nbasmon)/))

  END SUBROUTINE routing_subgrid_htu_xios_initialize

!! ================================================================================================================================
!! SUBROUTINE   : routing_subgrid_htu_main 
!!
!>\BRIEF          This module routes the water over the continents (runoff and
!!                drainage produced by the hydrol module) into the oceans. 
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

SUBROUTINE routing_subgrid_htu_main(kjit, nbpt, index, &
       & lalo, neighbours, resolution, contfrac, totfrac_nobio, veget_max, floodout, runoff, &
       & drainage, transpot, precip_rain, humrel, k_litt, flood_frac, flood_res, &
       & tempdiag, reinf_slope, returnflow, reinfiltration, irrigation, riverflow, coastalflow, rest_id, hist_id, hist2_id)

    IMPLICIT NONE

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)     :: kjit                 !! Time step number (unitless)
    INTEGER(i_std), INTENT(in)     :: nbpt                 !! Domain size (unitless)
    INTEGER(i_std),INTENT(in)      :: rest_id              !! Restart file identifier (unitless)
    INTEGER(i_std),INTENT(in)      :: hist_id              !! Access to history file (unitless)
    INTEGER(i_std),INTENT(in)      :: hist2_id             !! Access to history file 2 (unitless)
    INTEGER(i_std), INTENT(in)     :: index(nbpt)          !! Indices of the points on the map (unitless)
    REAL(r_std), INTENT(in)        :: lalo(nbpt,2)         !! Vector of latitude and longitudes (beware of the order !)
    INTEGER(i_std), INTENT(in)     :: neighbours(nbpt,NbNeighb)   !! Vector of neighbours for each grid point (1=N, 2=NE, 3=E, 4=SE, 5=S, 6=SW, 7=W, 8=NW) (unitless)
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
    REAL(r_std), INTENT(in)        :: tempdiag(nbpt,ngrnd) !! Diagnostic soil temperature profile
    REAL(r_std), INTENT(in)        :: reinf_slope(nbpt)    !! Coefficient which determines the reinfiltration ratio in the grid box due to flat areas (unitless;0-1)

    !! 0.2 Output variables
    REAL(r_std), INTENT(out)       :: returnflow(nbpt)     !! The water flow from lakes and swamps which returns to the grid box.
                                                           !! This water will go back into the hydrol module to allow re-evaporation (kg/m^2/dt)
    REAL(r_std), INTENT(out)       :: reinfiltration(nbpt) !! Water flow from ponds and floodplains which returns to the grid box (kg/m^2/dt)
    REAL(r_std), INTENT(out)       :: irrigation(nbpt)     !! Irrigation flux. This is the water taken from the reservoirs and beeing put into the upper layers of the soil (kg/m^2/dt)
    REAL(r_std), INTENT(out)       :: riverflow(nbpt)      !! Outflow of the major rivers. The flux will be located on the continental grid but this should be a coastal point (kg/dt)
    REAL(r_std), INTENT(out)       :: coastalflow(nbpt)    !! Outflow on coastal points by small basins. This is the water which flows in a disperse way into the ocean (kg/dt)
    REAL(r_std), INTENT(out)       :: flood_frac(nbpt)     !! Flooded fraction of the grid box (unitless;0-1)
    REAL(r_std), INTENT(out)       :: flood_res(nbpt)      !! Diagnostic of water amount in the floodplains reservoir (kg)

    !! 0.3 Local variables
    REAL(r_std), DIMENSION(nbpt)   :: return_lakes         !! Water from lakes flowing back into soil moisture (kg/m^2/dt)

    INTEGER(i_std)                 :: ig, jv, ib           !! Indices (unitless)
    REAL(r_std), DIMENSION(nbpt)   :: tot_vegfrac_nowoody  !! Total fraction occupied by grass (0-1,unitless)

    REAL(r_std), DIMENSION(nbpt)   :: fast_diag_old        !! Reservoir in the beginning of the time step
    REAL(r_std), DIMENSION(nbpt)   :: slow_diag_old        !! Reservoir in the beginning of the time step
    REAL(r_std), DIMENSION(nbpt)   :: stream_diag_old      !! Reservoir in the beginning of the time step
    REAL(r_std), DIMENSION(nbpt)   :: lake_diag_old        !! Reservoir in the beginning of the time step
    REAL(r_std), DIMENSION(nbpt)   :: pond_diag_old        !! Reservoir in the beginning of the time step
    REAL(r_std), DIMENSION(nbpt)   :: flood_diag_old       !! Reservoir in the beginning of the time step

    !! For water budget check in the three routing reservoirs (positive if input > output)
    !! Net fluxes averaged over each grid cell in kg/m^2/dt
    REAL(r_std), DIMENSION(nbpt)   :: netflow_stream_diag  !! Input - Output flow to stream reservoir
    REAL(r_std), DIMENSION(nbpt)   :: netflow_fast_diag    !! Input - Output flow to fast reservoir
    REAL(r_std), DIMENSION(nbpt)   :: netflow_slow_diag    !! Input - Output flow to slow reservoir
    !
    REAL(r_std), DIMENSION(nbpt,nbasmax)   :: stemp_total_tend, stemp_advec_tend, stemp_relax_tend
    !
!_ ================================================================================================================================

    ! Save reservoirs in beginning of time step to calculate the water budget
    fast_diag_old   = fast_diag
    slow_diag_old   = slow_diag
    stream_diag_old = stream_diag
    lake_diag_old   = lake_diag
    pond_diag_old   = pond_diag
    flood_diag_old  = flood_diag

    !
    !! Computes the variables averaged between routing time steps and which will be used in subsequent calculations
    !
    floodout_mean(:) = floodout_mean(:) + floodout(:)
    runoff_mean(:) = runoff_mean(:) + runoff(:)
    drainage_mean(:) = drainage_mean(:) + drainage(:)
    floodtemp(:) = tempdiag(:,floodtemp_lev)
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
    tempdiag_mean(:,:) = tempdiag_mean(:,:) + tempdiag(:,:)*dt_sechiba/dt_routing
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
       !! Computes the transport of water in the various reservoirs
       !
       CALL routing_hr_flow(nbpt, dt_routing, lalo, floodout_mean, runoff_mean, drainage_mean, &
            & vegtot_mean, totnobio_mean, transpot_mean, precip_mean, humrel_mean, k_litt_mean, floodtemp, &
            & tempdiag_mean, reinf_slope, lakeinflow_mean, returnflow_mean, reinfiltration_mean, &
            & irrigation_mean, riverflow_mean, coastalflow_mean, hydrographs, slowflow_diag, flood_frac, &
            & flood_res, netflow_stream_diag, netflow_fast_diag, netflow_slow_diag, &
            & stemp_total_tend, stemp_advec_tend, stemp_relax_tend)
       !
       !! Responsible for storing the water in lakes
       !
       CALL routing_hr_lake(nbpt, dt_routing, lakeinflow_mean, humrel_mean, return_lakes)
       !
       returnflow_mean(:) = returnflow_mean(:) + return_lakes(:)

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
       tempdiag_mean(:,:) = zero
       vegtot_mean(:) = zero

       ! Change the units of the routing fluxes from kg/dt_routing into kg/dt_sechiba
       hydrographs(:) = hydrographs(:)/dt_routing*dt_sechiba
       HTUhgmon(:,:) = HTUhgmon(:,:)/dt_routing*dt_sechiba
       slowflow_diag(:) = slowflow_diag(:)/dt_routing*dt_sechiba

       ! Change the units of the routing fluxes from kg/m^2/dt_routing into kg/m^2/dt_sechiba
       returnflow_mean(:) = returnflow_mean(:)/dt_routing*dt_sechiba
       reinfiltration_mean(:) = reinfiltration_mean(:)/dt_routing*dt_sechiba
       irrigation_mean(:) = irrigation_mean(:)/dt_routing*dt_sechiba
       irrig_netereq(:) = irrig_netereq(:)/dt_routing*dt_sechiba
       
       ! Change units as above but at the same time transform the kg/dt_routing to m^3/dt_sechiba
       riverflow_mean(:) = riverflow_mean(:)/dt_routing*dt_sechiba/mille
       coastalflow_mean(:) = coastalflow_mean(:)/dt_routing*dt_sechiba/mille

       ! Water budget residu of the three routing reservoirs (in kg/m^2/s)
       ! Note that these diagnostics are done using local variables only calculated 
       ! during the time steps when the routing is calculated
       CALL xios_orchidee_send_field("wbr_stream",(stream_diag - stream_diag_old - netflow_stream_diag)/dt_routing)
       CALL xios_orchidee_send_field("wbr_fast",  (fast_diag   - fast_diag_old - netflow_fast_diag)/dt_routing)
       CALL xios_orchidee_send_field("wbr_slow",  (slow_diag   - slow_diag_old - netflow_slow_diag)/dt_routing)
       CALL xios_orchidee_send_field("wbr_lake",  (lake_diag   - lake_diag_old - &
            lakeinflow_mean + return_lakes)/dt_routing)
       CALL xios_orchidee_send_field("StreamT_TotTend", stemp_total_tend)
       CALL xios_orchidee_send_field("StreamT_AdvTend", stemp_advec_tend)
       CALL xios_orchidee_send_field("StreamT_RelTend", stemp_relax_tend)
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
    !
    CALL xios_orchidee_send_field("mask_coast",mask_coast)

    IF ( do_irrigation ) THEN 
       CALL xios_orchidee_send_field("irrigmap",irrigated)
    ENDIF
       
    IF ( do_floodplains ) THEN
       !! May be improved by performing the operation with XIOS
       floodmap(:) = 0.0
       DO ig=1,nbpt
          floodmap(ig) = SUM(floodplains(ig,:)) / (area(ig)*contfrac(ig))
       END DO
       CALL xios_orchidee_send_field("floodmap",floodmap)
    ENDIF
       
    IF ( doswamps ) THEN
       CALL xios_orchidee_send_field("swampmap",swamp)
    ENDIF
       
    !
    ! Water storage in reservoirs [kg/m^2]
    CALL xios_orchidee_send_field("fastr",fast_diag)
    CALL xios_orchidee_send_field("slowr",slow_diag)
    CALL xios_orchidee_send_field("streamr",stream_diag)
    CALL xios_orchidee_send_field("laker",lake_diag)
    CALL xios_orchidee_send_field("pondr",pond_diag)
    CALL xios_orchidee_send_field("floodr",flood_diag)
    CALL xios_orchidee_send_field("floodh",flood_height_diag)

    ! FLOODPLAINS
    CALL xios_orchidee_send_field("flood_frac",flood_frac)

    ! Difference between the end and the beginning of the routing time step [kg/m^2]
    CALL xios_orchidee_send_field("delfastr",   fast_diag   - fast_diag_old)
    CALL xios_orchidee_send_field("delslowr",   slow_diag   - slow_diag_old)
    CALL xios_orchidee_send_field("delstreamr", stream_diag - stream_diag_old)
    CALL xios_orchidee_send_field("dellaker",   lake_diag   - lake_diag_old)
    CALL xios_orchidee_send_field("delpondr",   pond_diag   - pond_diag_old)
    CALL xios_orchidee_send_field("delfloodr",  flood_diag  - flood_diag_old)

    ! Water fluxes converted from kg/m^2/dt_sechiba into kg/m^2/s 
    CALL xios_orchidee_send_field("irrigation",irrigation/dt_sechiba)
    CALL xios_orchidee_send_field("netirrig",irrig_netereq/dt_sechiba)
    CALL xios_orchidee_send_field("riversret",returnflow/dt_sechiba)
    CALL xios_orchidee_send_field("reinfiltration",reinfiltration/dt_sechiba)

    ! Transform from kg/dt_sechiba into m^3/s
    CALL xios_orchidee_send_field("hydrographs",hydrographs/mille/dt_sechiba)
    CALL xios_orchidee_send_field("htuhgmon",HTUhgmon/mille/dt_sechiba)
    CALL xios_orchidee_send_field("htutempmon",HTUtempmon)
    CALL xios_orchidee_send_field("hydrotemp", hydrotemp)
    CALL xios_orchidee_send_field("streamlimit", streamlimit)

    CALL xios_orchidee_send_field("slowflow",slowflow_diag/mille/dt_sechiba) ! previous id name: Qb
    CALL xios_orchidee_send_field("coastalflow",coastalflow/dt_sechiba)
    CALL xios_orchidee_send_field("riverflow",riverflow/dt_sechiba)

    IF ( .NOT. xios_orchidee_ok) THEN
       IF ( .NOT. almaoutput ) THEN
          !
          CALL histwrite_p(hist_id, 'riversret', kjit, returnflow, nbpt, index)
          IF (do_floodplains .OR. doponds) THEN
             CALL histwrite_p(hist_id, 'reinfiltration', kjit, reinfiltration, nbpt, index)
          ENDIF
          CALL histwrite_p(hist_id, 'hydrographs', kjit, hydrographs/mille, nbpt, index)
          !
          CALL histwrite_p(hist_id, 'fastr', kjit, fast_diag, nbpt, index)
          CALL histwrite_p(hist_id, 'slowr', kjit, slow_diag, nbpt, index)
          CALL histwrite_p(hist_id, 'streamr', kjit, stream_diag, nbpt, index)
          IF ( do_floodplains ) THEN
             CALL histwrite_p(hist_id, 'floodr', kjit, flood_diag, nbpt, index)
             CALL histwrite_p(hist_id, 'floodh', kjit, flood_height_diag, nbpt, index)
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
          CALL histwrite_p(hist_id, 'SurfStor', kjit, flood_diag+pond_diag+lake_diag, nbpt, index)
          CALL histwrite_p(hist_id, 'Dis', kjit, hydrographs/mille, nbpt, index)
          !
          CALL histwrite_p(hist_id, 'slowr', kjit, slow_diag, nbpt, index)
          CALL histwrite_p(hist_id, 'fastr', kjit, fast_diag, nbpt, index)
          CALL histwrite_p(hist_id, 'streamr', kjit, stream_diag, nbpt, index)
          CALL histwrite_p(hist_id, 'lakevol', kjit, lake_diag, nbpt, index)
          CALL histwrite_p(hist_id, 'pondr', kjit, pond_diag, nbpt, index)
          !
          IF ( do_irrigation ) THEN
             CALL histwrite_p(hist_id, 'Qirrig', kjit, irrigation, nbpt, index)
             CALL histwrite_p(hist_id, 'Qirrig_req', kjit, irrig_netereq, nbpt, index)
          ENDIF
          !
       ENDIF
       IF ( hist2_id > 0 ) THEN
          IF ( .NOT. almaoutput) THEN
             !
             CALL histwrite_p(hist2_id, 'riversret', kjit, returnflow, nbpt, index)
             IF (do_floodplains .OR. doponds) THEN
                CALL histwrite_p(hist2_id, 'reinfiltration', kjit, reinfiltration, nbpt, index)
             ENDIF
             CALL histwrite_p(hist2_id, 'hydrographs', kjit, hydrographs/mille, nbpt, index)
             !
             CALL histwrite_p(hist2_id, 'fastr', kjit, fast_diag, nbpt, index)
             CALL histwrite_p(hist2_id, 'slowr', kjit, slow_diag, nbpt, index)
             IF ( do_floodplains ) THEN
                CALL histwrite_p(hist2_id, 'floodr', kjit, flood_diag, nbpt, index)
                CALL histwrite_p(hist2_id, 'floodh', kjit, flood_height_diag, nbpt, index)
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
             CALL histwrite_p(hist2_id, 'SurfStor', kjit, flood_diag+pond_diag+lake_diag, nbpt, index)
             CALL histwrite_p(hist2_id, 'Dis', kjit, hydrographs/mille, nbpt, index)
             !
          ENDIF
       ENDIF
    ENDIF
    !
    !
  END SUBROUTINE routing_subgrid_htu_main
  
  !!  =============================================================================================================================
  !! SUBROUTINE:         routing_subgrid_htu_finalize
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

  SUBROUTINE routing_subgrid_htu_finalize( kjit, nbpt, rest_id, flood_frac, flood_res )
    
    IMPLICIT NONE
    
    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)     :: kjit                 !! Time step number (unitless)
    INTEGER(i_std), INTENT(in)     :: nbpt                 !! Domain size (unitless)
    INTEGER(i_std),INTENT(in)      :: rest_id              !! Restart file identifier (unitless)
    REAL(r_std), INTENT(in)        :: flood_frac(nbpt)     !! Flooded fraction of the grid box (unitless;0-1)
    REAL(r_std), INTENT(in)        :: flood_res(nbpt)      !! Diagnostic of water amount in the floodplains reservoir (kg)
    
    !! 0.2 Local variables      

!_ ================================================================================================================================
    
    !
    ! Write restart variables
    !
    CALL restput_p (rest_id, 'routingcounter', kjit, time_counter)

    CALL restput_p (rest_id, 'streamtcst', kjit, stream_tcst)
    CALL restput_p (rest_id, 'slowtcst', kjit, slow_tcst)
    CALL restput_p (rest_id, 'fasttcst', kjit, fast_tcst)
    CALL restput_p (rest_id, 'floodtcst', kjit, flood_tcst)
    CALL restput_p (rest_id, 'swampcst', kjit, swamp_cst)

    CALL restput_p (rest_id, 'lim_floodcri', kjit, lim_floodcri)
    
    CALL restput_p (rest_id, 'nbasmax', kjit, nbasmax)
    CALL restput_p (rest_id, 'nbasmon', kjit, nbasmon)
    CALL restput_p (rest_id, 'inflows', kjit, inflows)
    
    CALL restput_p (rest_id, 'routingarea', nbp_glo, nbasmax, 1, kjit, routing_area, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'routetogrid', nbp_glo, nbasmax, 1, kjit, REAL(route_togrid,r_std), 'scatter', &
         nbp_glo, index_g)
    CALL restput_p (rest_id, 'routetobasin', nbp_glo, nbasmax, 1, kjit, REAL(route_tobasin,r_std), 'scatter', &
         nbp_glo, index_g)
    CALL restput_p (rest_id, 'routenbintobas', nbp_glo, nbasmax, 1, kjit, REAL(route_nbintobas,r_std), 'scatter', &
         nbp_glo, index_g)
    CALL restput_p (rest_id, 'basinid', nbp_glo, nbasmax, 1, kjit, REAL(global_basinid,r_std), 'scatter', &
         nbp_glo, index_g)
    CALL restput_p (rest_id, 'topoindex', nbp_glo, nbasmax, 1, kjit, topo_resid, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'topoindex_stream', nbp_glo, nbasmax, 1, kjit, stream_resid, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'fastres', nbp_glo, nbasmax, 1, kjit, fast_reservoir, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'slowres', nbp_glo, nbasmax, 1, kjit, slow_reservoir, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'streamres', nbp_glo, nbasmax, 1, kjit, stream_reservoir, 'scatter',nbp_glo,index_g)
    CALL restput_p (rest_id, 'floodres', nbp_glo, nbasmax, 1, kjit, flood_reservoir, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'floodh', nbp_glo, nbasmax, 1, kjit, flood_height, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'flood_frac_bas', nbp_glo, nbasmax, 1, kjit, flood_frac_bas, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'pond_frac', nbp_glo, 1, 1, kjit, pond_frac, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'flood_frac', nbp_glo, 1, 1, kjit, flood_frac, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'flood_res', nbp_glo, 1, 1, kjit, flood_res, 'scatter', nbp_glo, index_g)

    CALL restput_p (rest_id, 'fasttemp', nbp_glo, nbasmax, 1, kjit, fast_temp, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'slowtemp', nbp_glo, nbasmax, 1, kjit, slow_temp, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'streamtemp', nbp_glo, nbasmax, 1, kjit, stream_temp, 'scatter',nbp_glo,index_g)
 
    
    CALL restput_p (rest_id, 'lakeres', nbp_glo, 1, 1, kjit, lake_reservoir, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'pondres', nbp_glo, 1, 1, kjit, pond_reservoir, 'scatter',  nbp_glo, index_g)

    CALL restput_p (rest_id, 'lakeinflow', nbp_glo, 1, 1, kjit, lakeinflow_mean, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'returnflow', nbp_glo, 1, 1, kjit, returnflow_mean, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'reinfiltration', nbp_glo, 1, 1, kjit, reinfiltration_mean, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'riverflow', nbp_glo, 1, 1, kjit, riverflow_mean, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'coastalflow', nbp_glo, 1, 1, kjit, coastalflow_mean, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'hydrographs', nbp_glo, 1, 1, kjit, hydrographs, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'htuhgmon', nbp_glo, nbasmon, 1, kjit, HTUhgmon, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'slowflow_diag', nbp_glo, 1, 1, kjit, slowflow_diag, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'hydrotemp', nbp_glo, 1, 1, kjit, hydrotemp, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'htutempmon', nbp_glo, nbasmon, 1, kjit, HTUtempmon, 'scatter',  nbp_glo, index_g)
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
    CALL restput_p (rest_id, 'tempdiag_route', nbp_glo, ngrnd, 1, kjit, tempdiag_mean, 'scatter',  nbp_glo, index_g)

    CALL restput_p (rest_id, 'gridrephtu', nbp_glo, 1, 1, kjit, REAL(hydrodiag,r_std), 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'htudiag', nbp_glo, nbasmon, 1, kjit, REAL(HTUdiag_loc,r_std), 'scatter',  nbp_glo, index_g)
    
    IF ( do_irrigation ) THEN
       CALL restput_p (rest_id, 'irrigated', nbp_glo, 1, 1, kjit, irrigated, 'scatter',  nbp_glo, index_g)
       CALL restput_p (rest_id, 'irrigation', nbp_glo, 1, 1, kjit, irrigation_mean, 'scatter',  nbp_glo, index_g)
    ENDIF

    IF ( do_floodplains ) THEN
       CALL restput_p (rest_id, 'floodplains', nbp_glo, nbasmax, 1, kjit, floodplains, 'scatter',  nbp_glo, index_g)
       CALL restput_p (rest_id, 'floodcri', nbp_glo, nbasmax, 1, kjit, floodcri, 'scatter',  nbp_glo, index_g)
       CALL restput_p (rest_id, 'floodp_beta', nbp_glo, nbasmax, 1, kjit, fp_beta, 'scatter',  nbp_glo, index_g)
    ENDIF
    IF ( dofloodoverflow ) THEN
       CALL restput_p (rest_id, 'orog_min', nbp_glo, nbasmax, 1,kjit,orog_min, 'scatter',  nbp_glo, index_g)
    END IF
    IF ( doswamps ) THEN
       CALL restput_p (rest_id, 'swamp', nbp_glo, 1, 1, kjit, swamp, 'scatter',  nbp_glo, index_g)
    ENDIF
  
  END SUBROUTINE routing_subgrid_htu_finalize

!! ================================================================================================================================
!! SUBROUTINE 	: routing_hr_init
!!
!>\BRIEF         This subroutine allocates the memory and get the fixed fields from the restart file.
!!
!! DESCRIPTION (definitions, functional, design, flags) : None
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

  SUBROUTINE routing_hr_init(kjit, nbpt, index, returnflow, reinfiltration, irrigation, &
       &                  riverflow, coastalflow, flood_frac, flood_res, tempdiag, rest_id)
    !
    IMPLICIT NONE
    !
    ! interface description
    !
!! INPUT VARIABLES
    INTEGER(i_std), INTENT(in)                   :: kjit           !! Time step number (unitless)
    INTEGER(i_std), INTENT(in)                   :: nbpt           !! Domain size (unitless)
    INTEGER(i_std), DIMENSION (nbpt), INTENT(in) :: index          !! Indices of the points on the map (unitless)
    REAL(r_std), DIMENSION(nbpt,ngrnd),INTENT(in) :: tempdiag      !! Temperature profile in soil
    INTEGER(i_std), INTENT(in)                   :: rest_id        !! Restart file identifier (unitless)
    !
!! OUTPUT VARIABLES
    REAL(r_std), DIMENSION (nbpt),INTENT(out)    :: returnflow     !! The water flow from lakes and swamps which returns into the grid box.
                                                                   !! This water will go back into the hydrol module to allow re-evaporation (kg/m^2/dt)
    REAL(r_std), DIMENSION (nbpt),INTENT(out)    :: reinfiltration !! Water flow from ponds and floodplains which returns to the grid box (kg/m^2/dt)
    REAL(r_std), DIMENSION (nbpt),INTENT(out)    :: irrigation     !! Irrigation flux. This is the water taken from the reservoirs and beeing put into the upper layers of the soil.(kg/m^2/dt)
    REAL(r_std), DIMENSION (nbpt),INTENT(out)    :: riverflow      !! Outflow of the major rivers. The flux will be located on the continental grid but this should be a coastal point (kg/dt)
    REAL(r_std), DIMENSION (nbpt),INTENT(out)    :: coastalflow    !! Outflow on coastal points by small basins. This is the water which flows in a disperse way into the ocean (kg/dt)
    REAL(r_std), DIMENSION (nbpt),INTENT(out)    :: flood_frac     !! Flooded fraction of the grid box (unitless;0-1)
    REAL(r_std), DIMENSION (nbpt),INTENT(out)    :: flood_res      !! Diagnostic of water amount in the floodplains reservoir (kg)
    !
!! LOCAL VARIABLES
    CHARACTER(LEN=80)                            :: var_name       !! To store variables names for I/O (unitless)
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)     :: tmp_real_g     !! A temporary real array for the integers
    REAL(r_std), ALLOCATABLE, DIMENSION(:)       :: tmp_real       !
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)     :: tmp_real_g2
    REAL(r_std)                                  :: ratio          !! Diagnostic ratio to check that dt_routing is a multiple of dt_sechiba (unitless)
    REAL(r_std)                                  :: totarea        !! Total area of basin (m^2)
    INTEGER(i_std)                               :: ier, ig, im, ib, ipn(1), nbhtumon !! Indices (unitless)
    REAL(r_std)                                  :: nbasmon_tmp, nbasmax_tmp, inflows_tmp

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
    dt_routing = dt_sechiba
    CALL getin_p('DT_ROUTING', dt_routing)
    !
    !
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
    IF ( do_floodplains ) CALL getin_p('DO_FLOODINFILT', dofloodinfilt)
    !
    !Config Key   = CONDUCT_FACTOR_FP
    !Config Desc  = Adjustment factor for floodplains reinfiltration
    !Config If    = RIVER_ROUTING
    !Config Def   = n
    !Config Help  = Factor used to reduce the infiltration from the
    !Config         floodplains. For a value of 1, the infiltration is
    !Config         unchanged, for a value of 0 there is no infiltration.
    !Config Units = -
    !
    conduct_factor = 1.0
    IF ( do_floodplains ) CALL getin_p('CONDUCT_FACTOR_FP', conduct_factor)
    !
    !
    !Config Key   = DO_FLOODOVERFLOW
    !Config Desc  = Should floodplains overflow to upstream HTUs floodplains 
    !Config If    = RIVER_ROUTING
    !Config Def   = n
    !Config Help  = This parameters allows the user to ask the model
    !Config         to take into account the overflow of the  
    !Config         floodplains. The water can flow to the upstream
    !Config         floodplains reservoir if the current flood height
    !Config         is higher than the upstream one.
    !Config Units = [FLAG]
    !
    dofloodoverflow = .FALSE.
    IF ( do_floodplains ) CALL getin_p('DO_FLOODOVERFLOW', dofloodoverflow)
    !
    !Config Key   = OVERFLOW_REPETITION
    !Config Desc  = Repetition of overflow at each routing time step
    !Config If    = RIVER_ROUTING
    !Config Def   = n
    !Config Help  = This parameters allows the user to ask the model
    !Config         repeat the overflow a certain amount of time 
    !Config         in order to have more stability with lower
    !Config         overflow time step.
    !Config Units = [FLAG]
    !
    overflow_repetition = 1
    IF ( do_floodplains ) CALL getin_p('OVERFLOW_REPETITION', overflow_repetition)
    !
    !Config Key   = R_FLOODMAX
    !Config Desc  = Maximal values for R factor 
    !Config If    = DO_FLOODPLAINS
    !Config Def   = 0.5
    !Config Help  = R is the factor of reduction of the stream discharge
    !Config         if there is floodplains. This is the maximal value 
    !Config         when the HTU is fully filled.
    !Config         R = 1 -> discharge = 0
    !Config         R = 0 -> Maximal discharge
    !
    rfloodmax = 0.5
    IF ( do_floodplains ) CALL getin_p('R_FLOODMAX', rfloodmax)
    !
    !Config Key   = OVERFLOW_TCST
    !Config Desc  = Time Constant for overflow in day 
    !Config If    = DO_FLOODPLAINS
    !Config Def   = 1
    !Config Help  = OVERFLOW_TCST is the time constant  
    !Config         For the floodplains overflow
    !
    overflow_tcst = 1
    IF ( do_floodplains ) CALL getin_p('OVERFLOW_TCST', overflow_tcst)
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
    !
    !Config Key   = FLOOD_BETA
    !Config Desc  = Parameter to fix the shape of the floodplain  
    !Config If    = RIVER_ROUTING
    !Config Def   = 2.0
    !Config Help  = Parameter to fix the shape of the floodplain
    !Config         (>1 for convex edges, <1 for concave edges)
    !Config Units = [-] 

    ! ANTHONY OLD FLOODPLAINS
    !CALL getin_p("FLOOD_BETA", beta)
    !
    !Config Key   = POND_BETAP
    !Config Desc  = Ratio of the basin surface intercepted by ponds and the maximum surface of ponds
    !Config If    = RIVER_ROUTING
    !Config Def   = 0.5
    !Config Help  = 
    !Config Units = [-] 
    CALL getin_p("POND_BETAP", betap)    
    !
    !Config Key   = FLOOD_CRI
    !Config Desc  = Potential height for which all the basin is flooded
    !Config If    = DO_FLOODPLAINS or DO_PONDS
    !Config Def   = 2000.
    !Config Help  = 
    !Config Units = [mm] 

    ! ANTHONY OLD FLOODPLAINS
    !CALL getin_p("FLOOD_CRI", floodcri)
    !
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

    !
    !
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
    ! If the routing_graph file is available we will extract the information in the dimensions
    ! and parameters.
    !
    !Config Key   = ROUTING_FILE
    !Config Desc  = Name of file which contains the routing information graph on the model grid
    !Config If    = RIVER_ROUTING
    !Config Def   = routing.nc
    !Config Help  = The file provided here should allows to route the water from one HTU
    !Config         to another. The RoutingPP code needs to be used in order to generate 
    !Config         the routing graph for the model grid.
    !Config         More details on : https://gitlab.in2p3.fr/ipsl/lmd/intro/routingpp
    !Config Units = [FILE]
    !
    !graphfilename = 'routing_graph.nc'
    !CALL getin('ROUTING_FILE',graphfilename)
    !CALL routing_hr_graphinfo(graphfilename, nbasmax, inflows, nbasmon, undef_graphfile, stream_tcst, fast_tcst, slow_tcst, &
    !                          &                 flood_tcst, swamp_cst, lim_floodcri)
    ! At this stage we could have an option to force reading of graph
    !
    ! Constants which can be in the restart file
    !
    var_name ="routingcounter"
    CALL ioconf_setatt_p('UNITS', 's')
    CALL ioconf_setatt_p('LONG_NAME','Time counter for the routing scheme')
    CALL restget_p (rest_id, var_name, kjit, .TRUE., zero, time_counter)
    CALL setvar_p (time_counter, val_exp, 'NO_KEYWORD', zero)

    ! Parameters which are in the restart file
    IF (stream_tcst .LE. 0 ) THEN
       var_name ="streamtcst"
       CALL ioconf_setatt_p('UNITS', 's/km')
       CALL ioconf_setatt_p('LONG_NAME','Time constant for the stream reservoir')
       CALL restget_p (rest_id, var_name, kjit, .TRUE., zero, stream_tcst)
    ENDIF
    IF (slow_tcst .LE. 0 ) THEN
       var_name ="slowtcst"
       CALL ioconf_setatt_p('UNITS', 's/km')
       CALL ioconf_setatt_p('LONG_NAME','Time constant for the slow reservoir')
       CALL restget_p (rest_id, var_name, kjit, .TRUE., zero, slow_tcst)
    ENDIF
    IF (fast_tcst .LE. 0 ) THEN
       var_name ="fasttcst"
       CALL ioconf_setatt_p('UNITS', 's/km')
       CALL ioconf_setatt_p('LONG_NAME','Time constant for the fast reservoir')
       CALL restget_p (rest_id, var_name, kjit, .TRUE., zero, fast_tcst)
    ENDIF
    IF (flood_tcst .LE. 0 ) THEN
       var_name ="floodtcst"
       CALL ioconf_setatt_p('UNITS', 's/km')
       CALL ioconf_setatt_p('LONG_NAME','Time constant for the flood reservoir')
       CALL restget_p (rest_id, var_name, kjit, .TRUE., zero, flood_tcst)
    ENDIF
    IF (swamp_cst .LE. 0 ) THEN
       var_name ="swampcst"
       CALL ioconf_setatt_p('UNITS', '-')
       CALL ioconf_setatt_p('LONG_NAME','Fraction of the river transport that flows to the swamps')
       CALL restget_p (rest_id, var_name, kjit, .TRUE., zero, swamp_cst)
    ENDIF
    IF (lim_floodcri .LE. 0 ) THEN
       var_name ="lim_floodcri"
       CALL ioconf_setatt_p('UNITS', 'm')
       CALL ioconf_setatt_p('LONG_NAME','Minimal difference of orography consecutive floodplains HTUs')
       CALL restget_p (rest_id, var_name, kjit, .TRUE., zero, lim_floodcri)
    ENDIF
    !
    ! Number of HTUs
    !
    var_name ="nbasmax"
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Number of HTU per grid box')
    CALL restget_p (rest_id, var_name, kjit, .TRUE., zero, nbasmax_tmp)
    CALL routing_hr_restartconsistency(var_name, nbasmax, nbasmax_tmp)
    !
    ! Number of inflows
    !
    var_name ="inflows"
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Maximum number of inflows per HTU')
    CALL restget_p (rest_id, var_name, kjit, .TRUE., zero, inflows_tmp)
    CALL routing_hr_restartconsistency(var_name, inflows, inflows_tmp)
    !
    ! Dimension of HTU monitoring variable
    !
    var_name ="nbasmon"
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Number of HTU to be monitored')
    CALL restget_p (rest_id, var_name, kjit, .TRUE., zero, nbasmon_tmp)
    CALL routing_hr_restartconsistency(var_name, nbasmon, nbasmon_tmp)
    !
    ! Continuation of extraction from restart file.
    !
    ALLOCATE (routing_area_loc(nbpt,nbasmax), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for routing_area_loc','','')

    ALLOCATE (routing_area_glo(nbp_glo,nbasmax))
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for routing_area_glo','','')
    IF ( .NOT. ReadGraph ) THEN
       var_name = 'routingarea'
       IF (is_root_prc) THEN
          CALL ioconf_setatt('UNITS', 'm^2')
          CALL ioconf_setatt('LONG_NAME','Area of basin')
          CALL restget (rest_id, var_name, nbp_glo, nbasmax, 1, kjit, .TRUE., routing_area_glo, "gather", nbp_glo, index_g)
       ENDIF
       CALL scatter(routing_area_glo,routing_area_loc)
       routing_area=>routing_area_loc
    ENDIF
    CALL scatter(routing_area_glo,routing_area_loc)
    routing_area=>routing_area_loc
    
    
    IF ( do_floodplains ) THEN
      !!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! ANTHONY - BETA 
      ALLOCATE (fp_beta_loc(nbpt,nbasmax), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for fp_beta_loc','','')

      ALLOCATE (fp_beta_glo(nbp_glo,nbasmax))
      IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for fp_beta_glo','','')

      IF ( .NOT. ReadGraph ) THEN
         IF (is_root_prc) THEN
            var_name = 'floodp_beta'
            CALL ioconf_setatt('UNITS', '-')
            CALL ioconf_setatt('LONG_NAME','Beta parameter for floodplains')
            CALL restget (rest_id, var_name, nbp_glo, nbasmax, 1, kjit, .TRUE., fp_beta_glo, "gather", nbp_glo, index_g)
         ENDIF
         CALL scatter(fp_beta_glo,fp_beta_loc)
      END IF
      !!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! ANTHONY - h0 - floodcri 
      ALLOCATE (floodcri_loc(nbpt,nbasmax), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for floodcri_loc','','')

      ALLOCATE (floodcri_glo(nbp_glo,nbasmax))
      IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for floodcri_glo','','')

      IF ( .NOT. ReadGraph ) THEN
         IF (is_root_prc) THEN
            var_name = 'floodcri'
            CALL ioconf_setatt('UNITS', 'mm')
            CALL ioconf_setatt('LONG_NAME','Height of complete flood')
            CALL restget (rest_id, var_name, nbp_glo, nbasmax, 1, kjit, .TRUE., floodcri_glo, "gather", nbp_glo, index_g)
         END IF
         CALL scatter(floodcri_glo,floodcri_loc)
      ENDIF
    END IF

    ALLOCATE (tmp_real_g(nbp_glo,nbasmax), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for tmp_real_g','','')

    ALLOCATE (route_togrid_loc(nbpt,nbasmax), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for route_togrid_loc','','')
    ALLOCATE (route_togrid_glo(nbp_glo,nbasmax), stat=ier)      ! used in global in routing_hr_flow
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for route_togrid_glo','','')

    IF ( .NOT. ReadGraph ) THEN
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
       CALL bcast(route_togrid_glo)                      ! used in global in routing_hr_flow
       CALL scatter(route_togrid_glo,route_togrid_loc)
    ENDIF
    !
    ALLOCATE (route_tobasin_loc(nbpt,nbasmax), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for route_tobasin_loc','','')

    ALLOCATE (route_tobasin_glo(nbp_glo,nbasmax), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for route_tobasin_glo','','')

    IF ( .NOT. ReadGraph ) THEN
       IF (is_root_prc) THEN
          var_name = 'routetobasin'
          CALL ioconf_setatt('UNITS', '-')
          CALL ioconf_setatt('LONG_NAME','Basin in to which the water goes')
          CALL restget (rest_id, var_name, nbp_glo, nbasmax, 1, kjit, .TRUE., tmp_real_g, "gather", nbp_glo, index_g)
          route_tobasin_glo = undef_int
          WHERE ( tmp_real_g .LT. val_exp )
             route_tobasin_glo = NINT(tmp_real_g)
          ENDWHERE
          num_largest = COUNT(route_tobasin_glo .EQ. nbasmax+3)
       ENDIF
       CALL scatter(route_tobasin_glo,route_tobasin_loc)
       CALL bcast(num_largest)
    ENDIF
    !
    ! nbintobasin
    !
    ALLOCATE (route_nbintobas_loc(nbpt,nbasmax), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for route_nbintobas_loc','','')
    ALLOCATE (route_nbintobas_glo(nbp_glo,nbasmax), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for route_nbintobas_glo','','')

    IF ( .NOT. ReadGraph ) THEN
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
    ENDIF
    !
    ALLOCATE (global_basinid_loc(nbpt,nbasmax), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for global_basinid_loc','','')
    ALLOCATE (global_basinid_glo(nbp_glo,nbasmax), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for global_basinid_glo','','')

    IF ( .NOT. ReadGraph ) THEN
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
    ENDIF
    !
    ALLOCATE (topo_resid_loc(nbpt,nbasmax), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for topo_resid_loc','','')
    ALLOCATE (topo_resid_glo(nbp_glo,nbasmax), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for topo_resid_glo','','')

    IF ( .NOT. ReadGraph ) THEN
       IF (is_root_prc) THEN
          var_name = 'topoindex'
          CALL ioconf_setatt('UNITS', 'km')
          CALL ioconf_setatt('LONG_NAME','Topographic index of the residence time')
          CALL restget (rest_id, var_name, nbp_glo, nbasmax, 1, kjit, .TRUE., topo_resid_glo, "gather", nbp_glo, index_g)
       ENDIF
       CALL scatter(topo_resid_glo,topo_resid_loc)
    ENDIF
    !
    ALLOCATE (stream_resid_loc(nbpt,nbasmax), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for stream_resid_loc','','')
    ALLOCATE (stream_resid_glo(nbp_glo,nbasmax), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for stream_resid_glo','','')

    IF ( .NOT. ReadGraph ) THEN
       IF (is_root_prc) THEN
          var_name = 'topoindex_stream'
          CALL ioconf_setatt('UNITS', 'km')
          CALL ioconf_setatt('LONG_NAME','Topographic index of the residence time')
          CALL restget (rest_id, var_name, nbp_glo, nbasmax, 1, kjit, .TRUE., stream_resid_glo, "gather", nbp_glo, index_g)
          stream_maxresid=MAXVAL(stream_resid_glo, MASK=stream_resid_glo .LT. undef_graphfile)
       ENDIF
       CALL bcast(stream_maxresid)
       CALL scatter(stream_resid_glo,stream_resid_loc)
    ENDIF
    !
    ALLOCATE (fast_reservoir(nbpt,nbasmax), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for fast_reservoir','','')
    var_name = 'fastres'
    CALL ioconf_setatt_p('UNITS', 'Kg')
    CALL ioconf_setatt_p('LONG_NAME','Water in the fast reservoir')
    CALL restget_p (rest_id, var_name, nbp_glo, nbasmax, 1, kjit, .TRUE., fast_reservoir, "gather", nbp_glo, index_g)
    CALL setvar_p (fast_reservoir, val_exp, 'NO_KEYWORD', zero)

    ALLOCATE (slow_reservoir(nbpt,nbasmax), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for slow_reservoir','','')
    var_name = 'slowres'
    CALL ioconf_setatt_p('UNITS', 'Kg')
    CALL ioconf_setatt_p('LONG_NAME','Water in the slow reservoir')
    CALL restget_p (rest_id, var_name, nbp_glo, nbasmax, 1, kjit, .TRUE., slow_reservoir, "gather", nbp_glo, index_g)
    CALL setvar_p (slow_reservoir, val_exp, 'NO_KEYWORD', zero)

    ALLOCATE (stream_reservoir(nbpt,nbasmax), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for stream_reservoir','','')
    var_name = 'streamres'
    CALL ioconf_setatt_p('UNITS', 'Kg')
    CALL ioconf_setatt_p('LONG_NAME','Water in the stream reservoir')
    CALL restget_p (rest_id, var_name, nbp_glo, nbasmax, 1, kjit, .TRUE., stream_reservoir, "gather", nbp_glo, index_g)
    CALL setvar_p (stream_reservoir, val_exp, 'NO_KEYWORD', zero)

    ALLOCATE (flood_reservoir(nbpt,nbasmax), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for flood_reservoir','','')
    var_name = 'floodres'
    CALL ioconf_setatt_p('UNITS', 'Kg')
    CALL ioconf_setatt_p('LONG_NAME','Water in the flood reservoir')
    CALL restget_p (rest_id, var_name, nbp_glo, nbasmax, 1, kjit, .TRUE., flood_reservoir, "gather", nbp_glo, index_g)
    CALL setvar_p (flood_reservoir, val_exp, 'NO_KEYWORD', zero)

    ALLOCATE (flood_frac_bas(nbpt,nbasmax), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for flood_frac_bas','','')
    var_name = 'flood_frac_bas'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Flooded fraction per basin')
    CALL restget_p (rest_id, var_name, nbp_glo, nbasmax, 1, kjit, .TRUE., flood_frac_bas, "gather", nbp_glo, index_g)
    CALL setvar_p (flood_frac_bas, val_exp, 'NO_KEYWORD', zero)

    ALLOCATE (flood_height(nbpt, nbasmax), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for flood_height','','')
    var_name = 'floodh'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','')
    CALL restget_p (rest_id, var_name, nbp_glo, nbasmax, 1, kjit, .TRUE., flood_height, "gather", nbp_glo, index_g)
    CALL setvar_p (flood_height, val_exp, 'NO_KEYWORD', zero)
    
    ALLOCATE (pond_frac(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for pond_frac','','')
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
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for lake_reservoir','','')
    var_name = 'lakeres'
    CALL ioconf_setatt_p('UNITS', 'Kg')
    CALL ioconf_setatt_p('LONG_NAME','Water in the lake reservoir')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., lake_reservoir, "gather", nbp_glo, index_g)
    CALL setvar_p (lake_reservoir, val_exp, 'NO_KEYWORD', zero)
    
    ALLOCATE (pond_reservoir(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for pond_reservoir','','')
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
       IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for irrigated','','')
       var_name = 'irrigated'
       CALL ioconf_setatt_p('UNITS', 'm^2')
       CALL ioconf_setatt_p('LONG_NAME','Surface of irrigated area')
       CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., irrigated, "gather", nbp_glo, index_g)
       CALL setvar_p (irrigated, val_exp, 'NO_KEYWORD', undef_sechiba)
    ENDIF

    ALLOCATE (floodmap(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for floodmap','','')
    ALLOCATE (floodplains_loc(nbpt,nbasmax), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for floodplains_loc','','')
    ALLOCATE (floodplains_glo(nbp_glo,nbasmax), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for floodplains_glo','','')
    
    IF ( do_floodplains ) THEN    
       IF ( .NOT. ReadGraph ) THEN
          IF (is_root_prc) THEN
             var_name = 'floodplains'
             CALL ioconf_setatt_p('UNITS', 'm^2')
             CALL ioconf_setatt_p('LONG_NAME','Surface which can be flooded')
             CALL restget (rest_id, var_name, nbp_glo, nbasmax, 1, kjit, .TRUE., floodplains_glo, "gather", nbp_glo, index_g)
          END IF
          CALL scatter(floodplains_glo,floodplains_loc)
       END IF
    ELSE
       floodmap(:) = zero
       floodplains_glo(:,:) = zero
       floodplains_loc(:,:) = zero
    ENDIF
    !!! 
    !!! ANTHONY : OVERFLOW
    !!!
    IF ( dofloodoverflow) THEN
      ALLOCATE (orog_min_loc(nbpt,nbasmax), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for orog_min_loc','','')
      ALLOCATE (orog_min_glo(nbp_glo,nbasmax), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for orog_min_glo','','') 
      !!
      IF ( .NOT. ReadGraph ) THEN
         IF (is_root_prc) THEN
            var_name = 'orog_min'
            CALL ioconf_setatt('UNITS', 'm')
            CALL ioconf_setatt('LONG_NAME','HTU minimum orography')
            CALL restget (rest_id, var_name, nbp_glo, nbasmax, 1, kjit, .TRUE., orog_min_glo, "gather", nbp_glo, index_g)
         END IF 
         CALL scatter(orog_min_glo,orog_min_loc)
      END IF
      !
      ALLOCATE (route_innum_loc(nbpt,nbasmax), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for route_innum_loc','','')
      ALLOCATE (route_innum_glo(nbp_glo,nbasmax), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for route_innum_glo','','')
      CALL scatter(route_innum_glo,route_innum_loc)
      !
      ALLOCATE (route_ingrid_loc(nbpt,nbasmax, inflows), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for route_ingrid_loc','','')
      ALLOCATE (route_ingrid_glo(nbp_glo,nbasmax,inflows), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for route_ingrid_glo','','') 
      CALL scatter(route_ingrid_glo,route_ingrid_loc)
      !
      ALLOCATE (route_inbasin_loc(nbpt,nbasmax, inflows), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for route_inbasin_loc','','')
      ALLOCATE (route_inbasin_glo(nbp_glo,nbasmax, inflows), stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for route_inbasin_glo','','')
      CALL scatter(route_inbasin_glo,route_inbasin_loc)
   END IF
    !!!  
    !!! 
    !!!      
    IF ( doswamps ) THEN
       ALLOCATE (swamp(nbpt), stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for swamp','','')
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
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for lakeinflow_mean','','')
    var_name = 'lakeinflow'
    CALL ioconf_setatt_p('UNITS', 'Kg/dt')
    CALL ioconf_setatt_p('LONG_NAME','Lake inflow')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., lakeinflow_mean, "gather", nbp_glo, index_g)
    CALL setvar_p (lakeinflow_mean, val_exp, 'NO_KEYWORD', zero)
    
    ALLOCATE (returnflow_mean(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for returnflow_mean','','')
    var_name = 'returnflow'
    CALL ioconf_setatt_p('UNITS', 'Kg/m^2/dt')
    CALL ioconf_setatt_p('LONG_NAME','Deep return flux')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., returnflow_mean, "gather", nbp_glo, index_g)
    CALL setvar_p (returnflow_mean, val_exp, 'NO_KEYWORD', zero)
    returnflow(:) = returnflow_mean(:)
    
    ALLOCATE (reinfiltration_mean(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for reinfiltration_mean','','')
    var_name = 'reinfiltration'
    CALL ioconf_setatt_p('UNITS', 'Kg/m^2/dt')
    CALL ioconf_setatt_p('LONG_NAME','Top return flux')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., reinfiltration_mean, "gather", nbp_glo, index_g)
    CALL setvar_p (reinfiltration_mean, val_exp, 'NO_KEYWORD', zero)
    reinfiltration(:) = reinfiltration_mean(:)
    
    ALLOCATE (irrigation_mean(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for irrigation_mean','','')
    ALLOCATE (irrig_netereq(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for irrig_netereq','','')
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
    
    ALLOCATE (fast_temp(nbpt,nbasmax), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for fast_temp','','')
    var_name = 'fasttemp'
    CALL ioconf_setatt_p('UNITS', 'K')
    CALL ioconf_setatt_p('LONG_NAME','Water temperature in the fast reservoir')
    CALL restget_p (rest_id, var_name, nbp_glo, nbasmax, 1, kjit, .TRUE., fast_temp, "gather", nbp_glo, index_g)

    ALLOCATE (slow_temp(nbpt,nbasmax), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for slow_temp','','')
    var_name = 'slowtemp'
    CALL ioconf_setatt_p('UNITS', 'K')
    CALL ioconf_setatt_p('LONG_NAME','Water temperature in the slow reservoir')
    CALL restget_p (rest_id, var_name, nbp_glo, nbasmax, 1, kjit, .TRUE., slow_temp, "gather", nbp_glo, index_g)
       
    IF ( COUNT(fast_temp == val_exp) == nbpt*nbasmax ) THEN
       CALL groundwatertemp(nbpt, nbasmax, ngrnd, tempdiag, znt, dlt, fast_temp, slow_temp)
    ENDIF
       
    ALLOCATE (stream_temp(nbpt,nbasmax), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for stream_temp','','')
    var_name = 'streamtemp'
    CALL ioconf_setatt_p('UNITS', 'K')
    CALL ioconf_setatt_p('LONG_NAME','Water temperature in the stream reservoir')
    CALL restget_p (rest_id, var_name, nbp_glo, nbasmax, 1, kjit, .TRUE., stream_temp, "gather", nbp_glo, index_g)

    
    IF ( COUNT(stream_temp == val_exp) == nbpt*nbasmax ) THEN
       DO ig=1,nbpt 
         stream_temp(ig,:) = tempdiag(ig,1)
      ENDDO
    ENDIF
    
    ALLOCATE (riverflow_mean(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for riverflow_mean','','')
    var_name = 'riverflow'
    CALL ioconf_setatt_p('UNITS', 'Kg/dt')
    CALL ioconf_setatt_p('LONG_NAME','River flux into the sea')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., riverflow_mean, "gather", nbp_glo, index_g)
    CALL setvar_p (riverflow_mean, val_exp, 'NO_KEYWORD', zero)
    riverflow(:) = riverflow_mean(:)
    
    ALLOCATE (coastalflow_mean(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for coastalflow_mean','','')
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
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for floodtemp','','')
    floodtemp(:) = tempdiag(:,floodtemp_lev)
    
    ALLOCATE(hydrographs(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for hydrographs','','')
    var_name = 'hydrographs'
    CALL ioconf_setatt_p('UNITS', 'kg/dt_sechiba')
    CALL ioconf_setatt_p('LONG_NAME','Hydrograph at outlow of grid')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., hydrographs, "gather", nbp_glo, index_g)
    CALL setvar_p (hydrographs, val_exp, 'NO_KEYWORD', zero)
    
    ALLOCATE(hydrotemp(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for hydrotemp','','')
    var_name = 'hydrotemp'
    CALL ioconf_setatt_p('UNITS', 'K')
    CALL ioconf_setatt_p('LONG_NAME','Temperature of most significant river of grid')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., hydrotemp, "gather", nbp_glo, index_g)
    CALL setvar_p (hydrotemp, val_exp, 'NO_KEYWORD', ZeroCelsius)
    
    ALLOCATE(slowflow_diag(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for slowflow_diag','','')
    var_name = 'slowflow_diag'
    CALL ioconf_setatt_p('UNITS', 'kg/dt_sechiba')
    CALL ioconf_setatt_p('LONG_NAME','Slowflow hydrograph at outlow of grid')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE.,slowflow_diag, "gather", nbp_glo, index_g)
    CALL setvar_p (slowflow_diag, val_exp, 'NO_KEYWORD', zero)
    !
    ! Grid diagnostic at representative HTU
    !
    ALLOCATE(hydrodiag_loc(nbpt),hydrodiag_glo(nbp_glo),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for hydrodiag_glo','','')
    IF ( .NOT. ReadGraph ) THEN
       IF (is_root_prc) THEN
          ALLOCATE(tmp_real(nbp_glo))
          var_name = 'gridrephtu'
          CALL ioconf_setatt('UNITS', '-')
          CALL ioconf_setatt('LONG_NAME','Representative HTU for the grid')
          CALL restget(rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE.,tmp_real, "gather", nbp_glo, index_g)
          hydrodiag_glo(:) = 1
          WHERE ( tmp_real .LT. val_exp )
             hydrodiag_glo = NINT(tmp_real)
          ENDWHERE
          DEALLOCATE(tmp_real)
       ENDIF
       CALL scatter(hydrodiag_glo, hydrodiag_loc)
    ENDIF
    !
    ! Station diagnostics
    !
    ALLOCATE(HTUdiag_loc(nbpt,nbasmon), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for HTUdiag_loc','','')
    ALLOCATE(HTUdiag_glo(nbp_glo,nbasmon), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for HTUdiag_glo','','')
    ALLOCATE(tmp_real_g2(nbp_glo,nbasmon), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for tmp_real_g2','','')
    !
    IF (is_root_prc) THEN
       var_name = 'htudiag'
       CALL ioconf_setatt('UNITS', 'index')
       CALL ioconf_setatt('LONG_NAME','Index of HTU to be monitored')
       CALL restget(rest_id, var_name, nbp_glo, nbasmon, 1, kjit, .TRUE., tmp_real_g2, "gather", nbp_glo, index_g)
       HTUdiag_glo(:,:) = -1
       WHERE ( tmp_real_g2 .LT. val_exp )
          HTUdiag_glo  = NINT(tmp_real_g2)
       ENDWHERE
       nbhtumon = 0
       DO ig=1,nbp_glo
          DO im=1,nbasmon
             IF ( HTUdiag_glo(ig,im) > 0 ) THEN
                nbhtumon = nbhtumon + 1
             ENDIF
          ENDDO
       ENDDO
       WRITE(numout,*) "After restget : Found a total of ", nbhtumon, " HTUs to be monitored and written into HTUhgmon"
    ENDIF
    CALL scatter(HTUdiag_glo, HTUdiag_loc)
    CALL bcast(nbhtumon)
    DEALLOCATE(tmp_real_g2)
    !
    ALLOCATE(HTUhgmon(nbpt,nbasmon), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for HTUhgmon','','')
    ALLOCATE(HTUhgmon_glo(nbp_glo,nbasmon), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for HTUhgmon_glo','','')
    !
    IF (is_root_prc) THEN
       var_name = 'htuhgmon'
       CALL ioconf_setatt('UNITS', 'kg/dt_sechiba')
       CALL ioconf_setatt('LONG_NAME','Hydrograph at selected HTU of grid')
       CALL restget(rest_id, var_name, nbp_glo, nbasmon, 1, kjit, .TRUE., HTUhgmon_glo, "gather", nbp_glo, index_g)
       WHERE ( HTUhgmon_glo .GE. val_exp )
          HTUhgmon_glo = zero
       ENDWHERE
    ENDIF
    CALL scatter(HTUhgmon_glo, HTUhgmon)
    DEALLOCATE(HTUhgmon_glo)
    !
    ! Restart of the temperature monitoring
    !
    ALLOCATE(HTUtempmon(nbpt,nbasmon), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for HTUtempmon','','')
    ALLOCATE(HTUtempmon_glo(nbp_glo,nbasmon), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for HTUtempmon_glo','','')
    !
    IF (is_root_prc) THEN
       var_name = 'htutempmon'
       CALL ioconf_setatt('UNITS', 'K')
       CALL ioconf_setatt('LONG_NAME','Temperature at selected HTU of grid')
       CALL restget(rest_id, var_name, nbp_glo, nbasmon, 1, kjit, .TRUE., HTUtempmon_glo, "gather", nbp_glo, index_g)
       WHERE ( HTUtempmon_glo .GE. val_exp )
          HTUtempmon_glo = ZeroCelsius
       ENDWHERE
       HTUtempmon_glo(:,:) = ZeroCelsius
    ENDIF
    CALL scatter(HTUtempmon_glo, HTUtempmon)
    DEALLOCATE(HTUtempmon_glo)
    !
    ! The diagnostic variables, they are initialized from the above restart variables.
    !
    ALLOCATE(fast_diag(nbpt), slow_diag(nbpt), stream_diag(nbpt), flood_diag(nbpt), &
         & pond_diag(nbpt), lake_diag(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for fast_diag,..','','')
    ALLOCATE(flood_height_diag(nbpt))
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for fast_height_diag','','')
    !
    !
    ! Get from the restart the fluxes we accumulated.
    !
    ALLOCATE (floodout_mean(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for floodout_mean','','')
    var_name = 'floodout_route'
    CALL ioconf_setatt_p('UNITS', 'Kg')
    CALL ioconf_setatt_p('LONG_NAME','Accumulated flow out of floodplains for routing')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., floodout_mean, "gather", nbp_glo, index_g)
    CALL setvar_p (floodout_mean, val_exp, 'NO_KEYWORD', zero)
    
    ALLOCATE (runoff_mean(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for runoff_mean','','')
    var_name = 'runoff_route'
    CALL ioconf_setatt_p('UNITS', 'Kg')
    CALL ioconf_setatt_p('LONG_NAME','Accumulated runoff for routing')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., runoff_mean, "gather", nbp_glo, index_g)
    CALL setvar_p (runoff_mean, val_exp, 'NO_KEYWORD', zero)
    
    ALLOCATE(drainage_mean(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for drainage_mean','','')
    var_name = 'drainage_route'
    CALL ioconf_setatt_p('UNITS', 'Kg')
    CALL ioconf_setatt_p('LONG_NAME','Accumulated drainage for routing')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., drainage_mean, "gather", nbp_glo, index_g)
    CALL setvar_p (drainage_mean, val_exp, 'NO_KEYWORD', zero)
    
    ALLOCATE(transpot_mean(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for transpot_mean','','')
    var_name = 'transpot_route'
    CALL ioconf_setatt_p('UNITS', 'Kg/m^2')
    CALL ioconf_setatt_p('LONG_NAME','Accumulated potential transpiration for routing/irrigation')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., transpot_mean, "gather", nbp_glo, index_g)
    CALL setvar_p (transpot_mean, val_exp, 'NO_KEYWORD', zero)

    ALLOCATE(precip_mean(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for precip_mean','','')
    var_name = 'precip_route'
    CALL ioconf_setatt_p('UNITS', 'Kg/m^2')
    CALL ioconf_setatt_p('LONG_NAME','Accumulated rain precipitation for irrigation')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., precip_mean, "gather", nbp_glo, index_g)
    CALL setvar_p (precip_mean, val_exp, 'NO_KEYWORD', zero)
    
    ALLOCATE(humrel_mean(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for humrel_mean','','')
    var_name = 'humrel_route'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Mean humrel for irrigation')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., humrel_mean, "gather", nbp_glo, index_g)
    CALL setvar_p (humrel_mean, val_exp, 'NO_KEYWORD', un)
    
    ALLOCATE(k_litt_mean(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for k_litt_mean','','')
    var_name = 'k_litt_route'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Mean cond. for litter')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., k_litt_mean, "gather", nbp_glo, index_g)
    CALL setvar_p (k_litt_mean, val_exp, 'NO_KEYWORD', zero)

    ALLOCATE(totnobio_mean(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for totnobio_mean','','')
    var_name = 'totnobio_route'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Last Total fraction of no bio for irrigation')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., totnobio_mean, "gather", nbp_glo, index_g)
    CALL setvar_p (totnobio_mean, val_exp, 'NO_KEYWORD', zero)
    
    ALLOCATE(vegtot_mean(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for vegtot_mean','','')
    var_name = 'vegtot_route'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Last Total fraction of vegetation')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., vegtot_mean, "gather", nbp_glo, index_g)
    CALL setvar_p (vegtot_mean, val_exp, 'NO_KEYWORD', un)
    !
    ALLOCATE(tempdiag_mean(nbpt,ngrnd), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_init','Pb in allocate for tempdiag_mean','','')
    var_name = 'tempdiag_route'
    CALL ioconf_setatt_p('UNITS', 'K')
    CALL ioconf_setatt_p('LONG_NAME','Mean temperature profile')
    CALL restget_p (rest_id, var_name, nbp_glo, ngrnd, 1, kjit, .TRUE., tempdiag_mean, "gather", nbp_glo, index_g)
    CALL setvar_p (tempdiag_mean, val_exp, 'NO_KEYWORD', Zero)
    !
    DEALLOCATE(tmp_real_g)
    !
    ! Other variables
    !
    ALLOCATE(streamlimit(nbpt), stat=ier)
    !
  END SUBROUTINE routing_hr_init
  !
!! ================================================================================================================================
!! SUBROUTINE 	: routing_hr_restartconsistency
!!
!>\BRIEF        : This subroutine will verify that the important dimensions for the routing exist in the restart file
!!                or can be read from the routing_graph.nc file. This ensures that if the restart and graph file are
!!                not consistent the latter is read and overwrite whatever was in the restart. Then the user will be
!!                using the new routing graph from the routing_graph.nc and not whatever is in the restart.
!! \n
!_ ================================================================================================================================

  SUBROUTINE routing_hr_restartconsistency(varname, dimgraph, dimrestart)
    !
    IMPLICIT NONE
    !
    !! INPUT VARIABLES
    CHARACTER(LEN=80), INTENT(in)       :: varname      !! Name of dimension
    INTEGER(i_std), INTENT(inout)       :: dimgraph     !! Dimension read in the graph file and also final result.
    REAL(r_std), INTENT(in)             :: dimrestart   !! Dimension read in the restart file.
    !
    ! Case when the routing_hr_graphinfo could not get any information from the routing_graph.nc
    IF ( dimgraph < un ) THEN
       IF ( dimrestart > zero ) THEN
          ! Only information from the restart.
          dimgraph = NINT(dimrestart)
       ELSE
          WRITE(*,*) "Problem : No information in the routing_graph file and no routing information in restart ", TRIM(varname)
          CALL ipslerr(3,'routing_hr_restartconsistency',&
               'No routing_graph file availble and no information in restart.', &
               'Cannot perform routing in ORCHIDEE.', ' ')
       ENDIF
       
    ! Information from the routing_graph.nc file exists !
    ELSE
       IF ( dimgraph .NE. NINT(dimrestart) ) THEN
          WRITE(*,*) "Problem for ", TRIM(varname)," in restart is not the same as in routing_graph.nc "
          WRITE(*,*) "Value of ", TRIM(varname), " in restart file : ", dimrestart
          WRITE(*,*) "Value of ", TRIM(varname), " in routing_graph.nc file : ", dimgraph
          CALL ipslerr(2,'routing_hr_restartconsistency',&
               'The value of dimension provided is not consistant with the one in routing_graph file.', &
               'We will read a new graph from the given file.', ' ')
          ReadGraph = .TRUE.
       ELSE
          !! Nothing to do
       ENDIF
    ENDIF
  END SUBROUTINE routing_hr_restartconsistency
!! ================================================================================================================================
!! SUBROUTINE 	: routing_subgrid_htu_clear
!!
!>\BRIEF        : This subroutine deallocates the block memory previously allocated.
!! \n
!_ ================================================================================================================================

  SUBROUTINE routing_subgrid_htu_clear()

    IF (ALLOCATED(routing_area_loc)) DEALLOCATE(routing_area_loc)
    IF (ALLOCATED(route_togrid_loc)) DEALLOCATE(route_togrid_loc)
    IF (ALLOCATED(route_tobasin_loc)) DEALLOCATE(route_tobasin_loc)
    IF (ALLOCATED(route_nbintobas_loc)) DEALLOCATE(route_nbintobas_loc)
    IF (ALLOCATED(global_basinid_loc)) DEALLOCATE(global_basinid_loc)
    IF (ALLOCATED(topo_resid_loc)) DEALLOCATE(topo_resid_loc)
    IF (ALLOCATED(stream_resid_loc)) DEALLOCATE(stream_resid_loc)
    IF (ALLOCATED(routing_area_glo)) DEALLOCATE(routing_area_glo)
    IF (ALLOCATED(route_togrid_glo)) DEALLOCATE(route_togrid_glo)
    IF (ALLOCATED(route_tobasin_glo)) DEALLOCATE(route_tobasin_glo)
    IF (ALLOCATED(route_nbintobas_glo)) DEALLOCATE(route_nbintobas_glo)
    IF (ALLOCATED(global_basinid_glo)) DEALLOCATE(global_basinid_glo)
    IF (ALLOCATED(topo_resid_glo)) DEALLOCATE(topo_resid_glo)
    IF (ALLOCATED(stream_resid_glo)) DEALLOCATE(stream_resid_glo)
    IF (ALLOCATED(fast_reservoir)) DEALLOCATE(fast_reservoir)
    IF (ALLOCATED(slow_reservoir)) DEALLOCATE(slow_reservoir)
    IF (ALLOCATED(stream_reservoir)) DEALLOCATE(stream_reservoir)

    IF (ALLOCATED(fast_temp)) DEALLOCATE(fast_temp)
    IF (ALLOCATED(slow_temp)) DEALLOCATE(slow_temp)
    IF (ALLOCATED(stream_temp)) DEALLOCATE(stream_temp)
    
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
    IF (ALLOCATED(tempdiag_mean)) DEALLOCATE(tempdiag_mean)
    IF (ALLOCATED(totnobio_mean)) DEALLOCATE(totnobio_mean)
    IF (ALLOCATED(vegtot_mean)) DEALLOCATE(vegtot_mean)
    IF (ALLOCATED(floodtemp)) DEALLOCATE(floodtemp)
    IF (ALLOCATED(hydrodiag_loc)) DEALLOCATE(hydrodiag_loc)
    IF (ALLOCATED(hydrodiag_glo)) DEALLOCATE(hydrodiag_glo)
    IF (ALLOCATED(hydrographs)) DEALLOCATE(hydrographs)
    IF (ALLOCATED(hydrotemp)) DEALLOCATE(hydrotemp)
    IF (ALLOCATED(HTUhgmon)) DEALLOCATE(HTUhgmon)
    IF (ALLOCATED(HTUhgmon_glo)) DEALLOCATE(HTUhgmon_glo)
    IF (ALLOCATED(HTUtempmon)) DEALLOCATE(HTUtempmon)
    IF (ALLOCATED(HTUtempmon_glo)) DEALLOCATE(HTUtempmon_glo)
    IF (ALLOCATED(slowflow_diag)) DEALLOCATE(slowflow_diag)
    IF (ALLOCATED(irrigation_mean)) DEALLOCATE(irrigation_mean)
    IF (ALLOCATED(irrigated)) DEALLOCATE(irrigated)
    IF (ALLOCATED(floodplains_glo)) DEALLOCATE(floodplains_glo)
    IF (ALLOCATED(floodplains_loc)) DEALLOCATE(floodplains_loc)
    IF (ALLOCATED(swamp)) DEALLOCATE(swamp)
    IF (ALLOCATED(fast_diag)) DEALLOCATE(fast_diag)
    IF (ALLOCATED(slow_diag)) DEALLOCATE(slow_diag)
    IF (ALLOCATED(stream_diag)) DEALLOCATE(stream_diag)
    IF (ALLOCATED(flood_diag)) DEALLOCATE(flood_diag)
    IF (ALLOCATED(flood_height_diag)) DEALLOCATE(flood_height_diag)
    IF (ALLOCATED(pond_diag)) DEALLOCATE(pond_diag)
    IF (ALLOCATED(lake_diag)) DEALLOCATE(lake_diag)
    !
    IF (ALLOCATED(route_innum_loc)) DEALLOCATE(route_innum_loc)
    IF (ALLOCATED(route_ingrid_loc)) DEALLOCATE(route_ingrid_loc)
    IF (ALLOCATED(route_inbasin_loc)) DEALLOCATE(route_inbasin_loc)
    IF (ALLOCATED(route_innum_glo)) DEALLOCATE(route_innum_glo)
    IF (ALLOCATED(route_ingrid_glo)) DEALLOCATE(route_ingrid_glo)
    IF (ALLOCATED(route_inbasin_glo)) DEALLOCATE(route_inbasin_glo)
    !
    IF (ALLOCATED(orog_min_loc)) DEALLOCATE(orog_min_loc)
    IF (ALLOCATED(orog_min_glo)) DEALLOCATE(orog_min_glo)
    !
    IF (ALLOCATED(floodcri_loc)) DEALLOCATE(floodcri_loc)
    IF (ALLOCATED(floodcri_glo)) DEALLOCATE(floodcri_glo)
    IF (ALLOCATED(fp_beta_loc)) DEALLOCATE(fp_beta_loc)
    IF (ALLOCATED(fp_beta_glo)) DEALLOCATE(fp_beta_glo)

  END SUBROUTINE routing_subgrid_htu_clear
  !

!! ================================================================================================================================
!! SUBROUTINE 	: routing_hr_flow
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

  SUBROUTINE routing_hr_flow(nbpt, dt_routing, lalo, floodout, runoff, drainage, &
       &                  vegtot, totnobio, transpot_mean, precip, humrel, k_litt, floodtemp, tempdiag, &
       &                  reinf_slope, lakeinflow, returnflow, reinfiltration, irrigation, riverflow, &
       &                  coastalflow, hydrographs, slowflow_diag, flood_frac, flood_res, &
       &                  netflow_stream_diag, netflow_fast_diag, netflow_slow_diag, &
       &                  stemp_total_tend, stemp_advec_tend, stemp_relax_tend)
    !
    IMPLICIT NONE
    !
!! INPUT VARIABLES
    INTEGER(i_std), INTENT(in)                   :: nbpt                      !! Domain size (unitless)
    REAL(r_std), INTENT (in)                     :: dt_routing                !! Routing time step (s)
    REAL(r_std), INTENT(in)                      :: lalo(nbpt,2)              !! Vector of latitude and longitudes
    REAL(r_std), INTENT(in)                      :: runoff(nbpt)              !! Grid-point runoff (kg/m^2/dt)
    REAL(r_std), INTENT(in)                      :: floodout(nbpt)            !! Grid-point flow out of floodplains (kg/m^2/dt)
    REAL(r_std), INTENT(in)                      :: drainage(nbpt)            !! Grid-point drainage (kg/m^2/dt)
    REAL(r_std), INTENT(in)                      :: vegtot(nbpt)              !! Potentially vegetated fraction (unitless;0-1)
    REAL(r_std), INTENT(in)                      :: totnobio(nbpt)            !! Other areas which can not have vegetation
    REAL(r_std), INTENT(in)                      :: transpot_mean(nbpt)       !! Mean potential transpiration of the vegetation (kg/m^2/dt)
    REAL(r_std), INTENT(in)                      :: precip(nbpt)              !! Rainfall (kg/m^2/dt)
    REAL(r_std), INTENT(in)                      :: humrel(nbpt)              !! Soil moisture stress, root extraction potential (unitless)
    REAL(r_std), INTENT(in)                      :: k_litt(nbpt)              !! Averaged conductivity for saturated infiltration in the 'litter' layer (kg/m^2/dt)
    REAL(r_std), INTENT(in)                      :: floodtemp(nbpt)           !! Temperature to decide if floodplains work (K)
    REAL(r_std), INTENT(in)                      :: tempdiag(nbpt,ngrnd)      !! Soil temperature profiles (K)
    REAL(r_std), INTENT(in)                      :: reinf_slope(nbpt)         !! Coefficient which determines the reinfiltration ratio in the grid box due to flat areas (unitless;0-1)
    REAL(r_std), INTENT(out)                     :: lakeinflow(nbpt)          !! Water inflow to the lakes (kg/dt)
    !
!! OUTPUT VARIABLES
    REAL(r_std), INTENT(out)                     :: returnflow(nbpt)          !! The water flow from lakes and swamps which returns into the grid box.
                                                                              !! This water will go back into the hydrol module to allow re-evaporation (kg/m^2/dt_routing)
    REAL(r_std), INTENT(out)                     :: reinfiltration(nbpt)      !! Water flow from ponds and floodplains which returns to the grid box (kg/m^2/dt)
    REAL(r_std), INTENT(out)                     :: irrigation(nbpt)          !! Irrigation flux. This is the water taken from the reservoirs and beeing put into the upper layers of the soil (kg/m^2/dt_routing)
    REAL(r_std), INTENT(out)                     :: riverflow(nbpt)           !! Outflow of the major rivers. The flux will be located on the continental grid but this should be a coastal point (kg/dt_routing)
    REAL(r_std), INTENT(out)                     :: coastalflow(nbpt)         !! Outflow on coastal points by small basins. This is the water which flows in a disperse way into the ocean (kg/dt_routing)
    REAL(r_std), INTENT(out)                     :: hydrographs(nbpt)         !! Hydrographs at the outflow of the grid box for major basins (kg/dt)
    REAL(r_std), INTENT(out)                     :: slowflow_diag(nbpt)       !! Hydrographs of slow_flow = routed slow_flow for major basins (kg/dt)
    REAL(r_std), INTENT(out)                     :: flood_frac(nbpt)          !! Flooded fraction of the grid box (unitless;0-1)
    REAL(r_std), INTENT(out)                     :: flood_res(nbpt)           !! Diagnostic of water amount in the floodplains reservoir (kg)

    REAL(r_std), INTENT(out)                     :: netflow_stream_diag(nbpt) !! Input - Output flow to stream reservoir
    REAL(r_std), INTENT(out)                     :: netflow_fast_diag(nbpt)   !! Input - Output flow to fast reservoir
    REAL(r_std), INTENT(out)                     :: netflow_slow_diag(nbpt)   !! Input - Output flow to slow reservoir
    REAL(r_std), INTENT(out)                     :: stemp_total_tend(nbpt, nbasmax)  !! Total tendency in GJ/s computed for the stream reservoir.
    REAL(r_std), INTENT(out)                     :: stemp_advec_tend(nbpt, nbasmax)  !! Tendency (GJ/s) produced by advection
    REAL(r_std), INTENT(out)                     :: stemp_relax_tend(nbpt, nbasmax)  !! Tendency (GJ/s) produced by relaxation
    !
!! LOCAL VARIABLES
    REAL(r_std), DIMENSION(nbpt, nbasmax)        :: fast_flow                 !! Outflow from the fast reservoir (kg/dt)
    REAL(r_std), DIMENSION(nbpt, nbasmax)        :: slow_flow                 !! Outflow from the slow reservoir (kg/dt)
    REAL(r_std), DIMENSION(nbpt, nbasmax)        :: stream_flow               !! Outflow from the stream reservoir (kg/dt)
    REAL(r_std), DIMENSION(nbpt, nbasmax)        :: flood_flow                !! Outflow from the floodplain reservoir (kg/dt)
    REAL(r_std), DIMENSION(nbpt, nbasmax)        :: pond_inflow               !! Inflow to the pond reservoir (kg/dt)
    REAL(r_std), DIMENSION(nbpt, nbasmax)        :: pond_drainage             !! Drainage from pond (kg/m^2/dt)
    REAL(r_std), DIMENSION(nbpt, nbasmax)        :: flood_drainage            !! Drainage from floodplains (kg/m^2/dt)
    REAL(r_std), DIMENSION(nbpt, nbasmax)        :: return_swamp              !! Inflow to the swamp (kg/dt)
    REAL(r_std), DIMENSION(nbpt, nbasmax)        :: source
    REAL(r_std), DIMENSION(nbpt, nbasmax)        :: ewh
    !
    ! Irrigation per basin
    !
    REAL(r_std), DIMENSION(nbpt, nbasmax)        :: irrig_needs               !! Total irrigation requirement (water requirements by the crop for its optimal growth) (kg)
    REAL(r_std), DIMENSION(nbpt, nbasmax)        :: irrig_actual              !! Possible irrigation according to the water availability in the reservoirs (kg)
    REAL(r_std), DIMENSION(nbpt, nbasmax)        :: irrig_deficit             !! Amount of water missing for irrigation (kg)
    REAL(r_std), DIMENSION(nbpt, nbasmax)        :: irrig_adduct              !! Amount of water carried over from other basins for irrigation (kg)
    !
    ! The transport terms are over a larger indexing space so that outlfows to ocean and lakes do not generate out of bounds issues.
    ! Non existing HTU have their index set to zero and their memory will end-up in index 0 of transport.
    !
    REAL(r_std), DIMENSION(nbpt, 0:nbasmax+3)    :: transport                 !! Water transport between basins (kg/dt)
    REAL(r_std), DIMENSION(nbp_glo, 0:nbasmax+3) :: transport_glo             !! Water transport between basins (kg/dt)
    REAL(r_std), DIMENSION(nbpt, 0:nbasmax+3)    :: transport_temp            !! Temperature transport between grids
    REAL(r_std), DIMENSION(nbp_glo, 0:nbasmax+3) :: transport_temp_glo        !! Temperature transport global for transfers
    !
    REAL(r_std)                                  :: oldtemp
    REAL(r_std)                                  :: oldstream
    INTEGER(i_std), SAVE                         :: nbunpy=0
    !
    REAL(r_std), DIMENSION(nbpt, nbasmax)        :: floods                    !! Water flow in to the floodplains (kg/dt)
    REAL(r_std), DIMENSION(nbpt, nbasmax)        :: potflood                  !! Potential inflow to the swamps (kg/dt)
    REAL(r_std), DIMENSION(nbpt)                 :: tobeflooded               !! Maximal surface which can be inundated in each grid box (m^2)
    REAL(r_std), DIMENSION(nbpt)                 :: totarea                   !! Total area of basin (m^2)
    REAL(r_std), DIMENSION(nbpt)                 :: totflood                  !! Total amount of water in the floodplains reservoir (kg)
    REAL(r_std), DIMENSION(nbasmax)              :: pond_excessflow           !! 
    REAL(r_std)                                  :: flow                      !! Outflow computation for the reservoirs (kg/dt)
    REAL(r_std)                                  :: floodindex                !! Fraction of grid box area inundated (unitless;0-1)
    REAL(r_std)                                  :: pondex                    !! 
    REAL(r_std)                                  :: stream_tot                !! Total water amount in the stream reservoirs (kg)
    REAL(r_std)                                  :: adduction                 !! Importation of water from a stream reservoir of a neighboring grid box (kg)
    REAL(r_std), DIMENSION(nbp_glo)              :: lake_overflow_g           !! Removed water from lake reservoir on global grid (kg/gridcell/dt_routing)
    REAL(r_std), DIMENSION(nbpt)                 :: lake_overflow             !! Removed water from lake reservoir on local grid (kg/gridcell/dt_routing)
    REAL(r_std), DIMENSION(nbpt)                 :: lake_overflow_coast       !! lake_overflow distributed on coast gridcells, only diag(kg/gridcell/dt_routing)
    REAL(r_std)                                  :: total_lake_overflow       !! Sum of lake_overflow over full grid (kg)
    REAL(r_std), DIMENSION(8,nbasmax)            :: streams_around            !! Stream reservoirs of the neighboring grid boxes (kg)
    INTEGER(i_std), DIMENSION(8)                 :: igrd                      !! 
    INTEGER(i_std), DIMENSION(2)                 :: ff                        !! 
    INTEGER(i_std), DIMENSION(1)                 :: fi                        !! 
    INTEGER(i_std)                               :: ig, ib, ib2, ig2, im      !! Indices (unitless)
    INTEGER(i_std)                               :: rtg, rtb, in, ing, inb,inf!! Indices (unitless)
    !INTEGER(i_std)                               :: numflood                  !!
    INTEGER(i_std)                               :: ier, negslow              !! Error handling
    INTEGER(i_std), DIMENSION(20)                :: negig, negib
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)     :: fast_flow_g               !! Outflow from the fast reservoir (kg/dt)
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)     :: slow_flow_g               !! Outflow from the slow reservoir (kg/dt)
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)     :: stream_flow_g             !! Outflow from the stream reservoir (kg/dt)
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)     :: fast_temp_g               !! Temperature of the fast reservoir (K)
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)     :: slow_temp_g               !! Temperature of the slow reservoir (K)
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)     :: stream_temp_g             !! Temperature of the stream reservoir (K)
    !REAL(r_std), ALLOCATABLE, DIMENSION(:,:)     :: flood_height_g            !! Floodplains height (m)
    !REAL(r_std), ALLOCATABLE, DIMENSION(:,:)     :: flood_frac_bas_g          !! Fraction of the HTU flooded 
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)     :: irrig_deficit_glo         !! Amount of water missing for irrigation (kg)
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)     :: stream_reservoir_glo      !! Water amount in the stream reservoir (kg)
    !REAL(r_std), ALLOCATABLE, DIMENSION(:,:)     :: flood_reservoir_glo      !! Water amount in the stream reservoir (kg)
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)     :: irrig_adduct_glo          !! Amount of water carried over from other basins for irrigation (kg)

    REAL(r_std)                                  :: reduced                   !! Discharge reduction due to floodplains
    REAL(r_std)                                  :: htmp, hscale              !! Water height scalingfor temperature relaxation
    REAL(r_std)                                  :: krelax, den
    !! PARAMETERS
    LOGICAL, PARAMETER                           :: check_reservoir = .FALSE. !! Logical to choose if we write informations when a negative amount of water is occurring in a reservoir (true/false)
!_ ================================================================================================================================
    !
    !
    hscale = 1.
    CALL getin_p('ROUTING_HSCALEKH',hscale)
    !
    transport(:,:) = zero
    transport_glo(:,:) = zero
    transport_temp(:,:) = zero !tp_00  its a transport, not a temperature !!
    transport_temp_glo(:,:) = zero !tp_00
    
    irrig_netereq(:) = zero
    irrig_needs(:,:) = zero
    irrig_actual(:,:) = zero
    irrig_deficit(:,:) = zero
    irrig_adduct(:,:) = zero
    totarea(:) = zero
    totflood(:) = zero
    !
    ! Compute all the fluxes
    !    
    DO ib=1,nbasmax
       DO ig=1,nbpt
          !
          totarea(ig) = totarea(ig) + routing_area(ig,ib)
          totflood(ig) = totflood(ig) + flood_reservoir(ig,ib)
       ENDDO
    ENDDO
          !
!> The outflow fluxes from the three reservoirs are computed. 
!> The outflow of volume of water Vi into the reservoir i is assumed to be linearly related to its volume.
!> The water travel simulated by the routing scheme is dependent on the water retention index topo_resid
!> given by a 0.5 degree resolution map for each pixel performed from a simplification of Manning's formula
!> (Dingman, 1994; Ducharne et al., 2003).
!> The resulting product of tcst (in s/km) and topo_resid (in km) represents the time constant (s)
!> which is an e-folding time, the time necessary for the water amount
!> in the stream reservoir to decrease by a factor e. Hence, it gives an order of
!> magnitude of the travel time through this reservoir between
!> the sub-basin considered and its downstream neighbor.
    !
    CALL groundwatertemp(nbpt, nbasmax, ngrnd, tempdiag, znt, dlt, fast_temp, slow_temp)
    !
    streamlimit(:) = zero
    !
    DO ib=1,nbasmax
       DO ig=1,nbpt
          IF ( route_tobasin(ig,ib) .GT. 0 ) THEN
             !
             ! Each of the fluxes is limited by the water in the reservoir and a small margin 
             ! (min_reservoir) to avoid rounding errors.
             !
             flow = MIN(fast_reservoir(ig,ib)/(topo_resid(ig,ib)*fast_tcst/dt_routing),&
                  & fast_reservoir(ig,ib)-min_sechiba)
             fast_flow(ig,ib) = MAX(flow, zero)

             flow = MIN(slow_reservoir(ig,ib)/(topo_resid(ig,ib)*slow_tcst/dt_routing),&
                  & slow_reservoir(ig,ib)-min_sechiba)
             slow_flow(ig,ib) = MAX(flow, zero)

             ! Need to adjust the reduction of the flow
             reduced = MAX(1-SQRT(MIN(flood_frac_bas(ig,ib),rfloodmax)), min_sechiba) ! Add the reduction flow parameter
             flow = stream_reservoir(ig,ib)/(stream_resid(ig,ib)*stream_tcst/dt_routing)*reduced
             flow = MIN(flow, stream_reservoir(ig,ib)-min_sechiba)
             stream_flow(ig,ib) = MAX(flow, zero)
             IF ( stream_flow(ig,ib) .GE. stream_reservoir(ig,ib)-min_sechiba .AND. stream_flow(ig,ib) > zero .AND. &
                  & routing_area(ig,ib) > zero ) THEN
                streamlimit(ig) = streamlimit(ig)+1.0
             ENDIF
             !
          ELSE
             fast_flow(ig,ib) = zero
             slow_flow(ig,ib) = zero
             stream_flow(ig,ib) = zero
          ENDIF
       ENDDO
    ENDDO
    !-
    !- Compute the fluxes out of the floodplains and ponds if they exist.
    !-
    IF (do_floodplains .OR. doponds) THEN
       DO ig=1,nbpt
          IF (flood_frac(ig) .GT. min_sechiba) THEN
             !!!! 
             ! PONDS : not actualized
             !
             !flow = MIN(floodout(ig)*totarea(ig)*pond_frac(ig)/flood_frac(ig), pond_reservoir(ig)+totflood(ig))
             !pondex = MAX(flow - pond_reservoir(ig), zero)
             !pond_reservoir(ig) = pond_reservoir(ig) - (flow - pondex) 
             !
             ! If demand was over reservoir size, we will take it out from floodplains
             !
             !pond_excessflow(:) = zero
             !DO ib=1,nbasmax
             !   pond_excessflow(ib) = MIN(pondex*flood_frac_bas(ig,ib)/(flood_frac(ig)-pond_frac(ig)),&
             !        &                    flood_reservoir(ig,ib))
             !   pondex = pondex - pond_excessflow(ib)
             !ENDDO
             !
             !IF ( pondex .GT. min_sechiba) THEN
             !   WRITE(numout,*) "Unable to redistribute the excess pond outflow over the water available in the floodplain."
             !   WRITE(numout,*) "Pondex = ", pondex
             !   WRITE(numout,*) "pond_excessflow(:) = ", pond_excessflow(:)
             !ENDIF
             !
             DO ib=1,nbasmax
                !
                ! when ponds actualized : add pond_excessflow to flow
                ! This is the flow out of the reservoir due to ET (+ pond excessflow(ig), suppressed here) 
                !flow = floodout(ig)*routing_area(ig,ib)*flood_frac_bas(ig,ib)/flood_frac(ig)
                flow = floodout(ig)*routing_area(ig,ib)*flood_frac_bas(ig,ib)
                !
                flood_reservoir(ig,ib) = flood_reservoir(ig,ib) - flow
                !
                !
                IF (flood_reservoir(ig,ib) .LT. min_sechiba) THEN
                   flood_reservoir(ig,ib) = zero
                ENDIF
                IF (pond_reservoir(ig) .LT. min_sechiba) THEN
                   pond_reservoir(ig) = zero
                ENDIF
             ENDDO
          ENDIF
       ENDDO
    ENDIF

    !-
    !- Computing the drainage and outflow from floodplains
!> Drainage from floodplains is depending on a averaged conductivity (k_litt) 
!> for saturated infiltration in the 'litter' layer. Flood_drainage will be
!> a component of the total reinfiltration that leaves the routing scheme.
    !-
    IF (do_floodplains) THEN
       IF (dofloodinfilt) THEN
          DO ib=1,nbasmax
             DO ig=1,nbpt
                flood_drainage(ig,ib) = MAX(zero, MIN(flood_reservoir(ig,ib), &
                & flood_frac_bas(ig,ib)* routing_area(ig,ib) * k_litt(ig) * &
                & conduct_factor * dt_routing/one_day))
                flood_reservoir(ig,ib) = flood_reservoir(ig,ib) - flood_drainage(ig,ib)
             ENDDO
          ENDDO
       ELSE
          DO ib=1,nbasmax
             DO ig=1,nbpt
                flood_drainage(ig,ib) = zero 
             ENDDO
          ENDDO
       ENDIF
!> Outflow from floodplains is computed depending a delay. This delay is characterized by a time constant
!> function of the surface of the floodplains and the product of topo_resid and flood_tcst. flood_tcst
!> has been calibrated through observations in the Niger Inner Delta (D'Orgeval, 2006).
!
       DO ib=1,nbasmax
          DO ig=1,nbpt
             IF ( route_tobasin(ig,ib) .GT. 0 ) THEN
                IF (flood_reservoir(ig,ib) .GT. min_sechiba) THEN
                   flow = MIN(flood_reservoir(ig,ib)/(stream_resid(ig,ib)*flood_tcst/dt_routing),&
                  & flood_reservoir(ig,ib)-min_sechiba)
                   flow = MAX(flow, zero)
                ELSE
                   flow = zero
                ENDIF
                flood_flow(ig,ib) = flow
             ELSE
                flood_flow(ig,ib) = zero
             ENDIF
          ENDDO
       ENDDO
    ELSE
       DO ib=1,nbasmax
          DO ig=1,nbpt
             flood_drainage(ig,ib) = zero
             flood_flow(ig,ib) = zero
             flood_reservoir(ig,ib) = zero
          ENDDO
       ENDDO
    ENDIF

    !-
    !- Computing drainage and inflow for ponds
!> Drainage from ponds is computed in the same way than for floodplains.
!> Reinfiltrated fraction from the runoff (i.e. the outflow from the fast reservoir)
!> is the inflow of the pond reservoir.
    !-
    IF (doponds) THEN
       ! If used, the slope coef is not used in hydrol for water2infilt
       DO ib=1,nbasmax
          DO ig=1,nbpt
             pond_inflow(ig,ib) = fast_flow(ig,ib) * reinf_slope(ig)
             pond_drainage(ig,ib) = MIN(pond_reservoir(ig)*routing_area(ig,ib)/totarea(ig), &
                  & pond_frac(ig)*routing_area(ig,ib)*k_litt(ig)*dt_routing/one_day)
             fast_flow(ig,ib) = fast_flow(ig,ib) - pond_inflow(ig,ib) 
          ENDDO
       ENDDO
    ELSE
       DO ib=1,nbasmax
          DO ig=1,nbpt
             pond_inflow(ig,ib) = zero
             pond_drainage(ig,ib) = zero
             pond_reservoir(ig) = zero
          ENDDO
       ENDDO
    ENDIF

    source(:,:) = fast_flow(:,:) + slow_flow(:,:) + stream_flow(:,:)
    CALL downstreamsum(nbpt, nbasmax, source, transport)
    source(:,:) = fast_flow(:,:)*fast_temp(:,:) + slow_flow(:,:)*slow_temp(:,:) +  &
            &                  stream_flow(:,:)*stream_temp(:,:)
    CALL downstreamsum(nbpt, nbasmax, source, transport_temp)
    !-
    !- Do the floodings - First initialize
    !-
    return_swamp(:,:)=zero
    floods(:,:)=zero
    !-
!> Over swamp areas, a fraction of water (return_swamp) is withdrawn from the river depending on the
!> parameter swamp_cst.
!> It will be transferred into soil moisture and thus does not return directly to the river.
    !
    !- 1. Swamps: Take out water from the river to put it to the swamps
    !-
    !
    IF ( doswamps ) THEN
       tobeflooded(:) = swamp(:)
       DO ib=1,nbasmax
          DO ig=1,nbpt
             potflood(ig,ib) = transport(ig,ib) 
             !
             IF ( tobeflooded(ig) > 0. .AND. potflood(ig,ib) > 0. .AND. floodtemp(ig) > tp_00 ) THEN
                !
                IF (routing_area(ig,ib) > tobeflooded(ig)) THEN
                   floodindex = tobeflooded(ig) / routing_area(ig,ib)
                ELSE
                   floodindex = 1.0
                ENDIF
                return_swamp(ig,ib) = swamp_cst * potflood(ig,ib) * floodindex
                !
                tobeflooded(ig) = tobeflooded(ig) - routing_area(ig,ib) 
                !
             ENDIF
          ENDDO
       ENDDO
    ENDIF
    !-
    !- 2. Floodplains: Update the reservoir with the flux computed above.
    !-
    IF ( do_floodplains ) THEN
       DO ig=1,nbpt
          DO ib=1,nbasmax
            IF (floodplains(ig, ib) .GT. min_sechiba .AND. floodtemp(ig) .GT. tp_00) THEN
                floods(ig,ib) = transport(ig,ib) - return_swamp(ig,ib) 
            ENDIF
          ENDDO
       ENDDO
    ENDIF
    !
    ! Update all reservoirs
!> The slow and deep reservoir (slow_reservoir) collect the deep drainage whereas the
!> fast_reservoir collects the computed surface runoff. Both discharge into a third reservoir
!> (stream_reservoir) of the next sub-basin downstream.
!> Water from the floodplains reservoir (flood_reservoir) flows also into the stream_reservoir of the next sub-basin downstream.
!> Water that flows into the pond_reservoir is withdrawn from the fast_reservoir.
    !
    negslow = 0
    DO ig=1,nbpt
       DO ib=1,nbasmax
          !
          fast_reservoir(ig,ib) =  fast_reservoir(ig,ib) + runoff(ig)*routing_area(ig,ib) - &
               & fast_flow(ig,ib) - pond_inflow(ig,ib)
          !
          slow_reservoir(ig,ib) = slow_reservoir(ig,ib) + drainage(ig)*routing_area(ig,ib) - &
               & slow_flow(ig,ib)
          !
          oldstream = stream_reservoir(ig, ib) * stream_temp(ig,ib)
          !
          stream_reservoir(ig,ib) = stream_reservoir(ig,ib) + flood_flow(ig,ib) + transport(ig,ib) - &
               & stream_flow(ig,ib) - return_swamp(ig,ib) - floods(ig,ib)
          !
          ! Diagnostics of the stream reservoir 
          !
          IF ( routing_area(ig,ib) > zero ) THEN
             ! 1000 to transform kg into m^3
             htmp = stream_reservoir(ig,ib)*1000/routing_area(ig,ib)
             ewh(ig,ib) = 1.0/(1.0+htmp*hscale)
          ELSE
             ewh(ig,ib) = 1.0
          ENDIF
          !
          !reste du calcul
          !
          krelax = ewh(ig,ib)
          !
          den = 1.0/(1.0+dt_routing*krelax)
          IF ( stream_reservoir(ig,ib) > 1.e-6 ) THEN
             oldtemp = stream_temp(ig,ib)
             stream_temp(ig,ib) = den * dt_routing * krelax * fast_temp(ig,ib) + &
                  & den * oldstream/stream_reservoir(ig,ib) + &
                  & den * transport_temp(ig, ib)/stream_reservoir(ig,ib) - &
                  & den * oldtemp*stream_flow(ig,ib)/stream_reservoir(ig,ib)
             !
             !Stream_temp [K], stream_reservoir [kg], WaterCp [J/g/K] yields tendencies in GJ/s
             ! 
             stemp_total_tend(ig,ib) = WaterCp*1.e-6*(stream_temp(ig,ib)*stream_reservoir(ig,ib) - oldstream)/dt_routing
             stemp_advec_tend(ig,ib) = WaterCp*1.e-6*(transport_temp(ig, ib) - oldtemp*stream_flow(ig,ib))/dt_routing
             stemp_relax_tend(ig,ib) = WaterCp*1.e-6*stream_reservoir(ig,ib)*krelax*(fast_temp(ig,ib)-stream_temp(ig,ib))
          ELSE
             stream_temp(ig,ib) = MAX(fast_temp(ig,ib), ZeroCelsius)
             stemp_total_tend(ig,ib) = zero
             stemp_advec_tend(ig,ib) = zero
             stemp_relax_tend(ig,ib) = zero
          ENDIF
          !
          flood_reservoir(ig,ib) = flood_reservoir(ig,ib) + floods(ig,ib) - &
               & flood_flow(ig,ib) 
          !
          pond_reservoir(ig) = pond_reservoir(ig) + pond_inflow(ig,ib) - pond_drainage(ig,ib)
          !
          IF ( flood_reservoir(ig,ib) .LT. zero ) THEN
             IF ( check_reservoir ) THEN
                WRITE(numout,*) "WARNING : negative flood reservoir at :", ig, ib, ". Problem is being corrected."
                WRITE(numout,*) "flood_reservoir, floods, flood_flow : ", flood_reservoir(ig,ib), floods(ig,ib), &
                     & flood_flow(ig,ib) 
             ENDIF
             stream_reservoir(ig,ib) = stream_reservoir(ig,ib) + flood_reservoir(ig,ib)
             flood_reservoir(ig,ib) = zero
          ENDIF
          !
          IF ( stream_reservoir(ig,ib) .LT. zero ) THEN
             IF ( check_reservoir ) THEN
                WRITE(numout,*) "WARNING : negative stream reservoir at :", ig, ib, ". Problem is being corrected."
                WRITE(numout,*) "stream_reservoir, flood_flow, transport : ", stream_reservoir(ig,ib), flood_flow(ig,ib), &
                     &  transport(ig,ib)
                WRITE(numout,*) "stream_flow, return_swamp, floods :", stream_flow(ig,ib), return_swamp(ig,ib), floods(ig,ib)
             ENDIF
             fast_reservoir(ig,ib) =  fast_reservoir(ig,ib) + stream_reservoir(ig,ib)
             stream_reservoir(ig,ib) = zero
          ENDIF
          !
          IF ( fast_reservoir(ig,ib) .LT. zero ) THEN
             IF ( check_reservoir ) THEN
                WRITE(numout,*) "WARNING : negative fast reservoir at :", ig, ib, ". Problem is being corrected."
                WRITE(numout,*) "fast_reservoir, runoff, fast_flow, ponf_inflow  : ", fast_reservoir(ig,ib), &
                     &runoff(ig), fast_flow(ig,ib), pond_inflow(ig,ib)
             ENDIF
             slow_reservoir(ig,ib) =  slow_reservoir(ig,ib) + fast_reservoir(ig,ib)
             fast_reservoir(ig,ib) = zero
          ENDIF

          IF ( slow_reservoir(ig,ib) .LT. - min_sechiba ) THEN
             IF ( negslow < 20 ) THEN
                negslow = negslow + 1
                negig(negslow) = ig
                negib(negslow) = ib
             ENDIF
          ENDIF

       ENDDO
    ENDDO

    IF ( negslow > 0 ) THEN
       DO ier = 1,negslow 
          ig = negig(ier)
          ib = negib(ier)
          WRITE(numout,*) 'WARNING : There is a negative reservoir at :', ig, ib,lalo(ig,:)
          WRITE(numout,*) 'WARNING : slowr, slow_flow, drainage', &
               & slow_reservoir(ig,ib), slow_flow(ig,ib), drainage(ig)
          WRITE(numout,*) 'WARNING : pondr, pond_inflow, pond_drainage', &
               & pond_reservoir(ig), pond_inflow(ig,ib), pond_drainage(ig,ib)
          CALL ipslerr_p(2, 'routing_hr_flow', 'WARNING negative slow_reservoir.','','')
       ENDDO
    ENDIF

    totflood(:) = zero
    DO ig=1,nbpt
       DO ib=1,nbasmax
          totflood(ig) = totflood(ig) + flood_reservoir(ig,ib)
       ENDDO
    ENDDO
    !
    ! ESTIMATE the flooded fraction
    !
    IF (do_floodplains .OR. doponds) THEN
      CALL routing_hr_flood(nbpt, flood_frac, totarea, totflood)
    ELSE
       flood_frac(:) = zero
       flood_height(:,:) = zero
       flood_frac_bas(:,:) = zero
    ENDIF
    

   !! ANTHONY : OVERFLOW
   !! CALCULATE TRANSFER BETWEEN FLOODPLAINS RESERVOIR
    IF (do_floodplains .AND. dofloodoverflow) Then
      ! The overflow is repeated "overflow_repetition" times
      ! This is in order to have more stability and
      ! be able to use lower "overflow_tcst".
         DO ier = 1,overflow_repetition
           CALL routing_hr_overflow(nbpt, nbasmax)
         END DO
       ! Once done we update the floodplains fraction and the floodplains height
       CALL routing_hr_flood(nbpt, flood_frac, totarea, totflood)
    END IF 


!-
!- Compute the total reinfiltration and returnflow to the grid box
!> A term of returnflow is computed including the water from the swamps that does not return directly to the river
!> but will be put into soil moisture (see hydrol module).
!> A term of reinfiltration is computed including the water that reinfiltrated from the ponds and floodplains areas.
!> It will be put into soil moisture (see hydrol module).
    !-
    IF (do_floodplains .OR. doswamps .OR. doponds) THEN
       returnflow(:) = zero
       reinfiltration(:) = zero
       !
       DO ib=1,nbasmax
          DO ig=1,nbpt
             returnflow(ig) =  returnflow(ig) + return_swamp(ig,ib)
             reinfiltration(ig) =  reinfiltration(ig) + pond_drainage(ig,ib) + flood_drainage(ig,ib) 
          ENDDO
       ENDDO
       !
       DO ig=1,nbpt
          returnflow(ig) = returnflow(ig)/totarea(ig)
          reinfiltration(ig) = reinfiltration(ig)/totarea(ig)
       ENDDO
    ELSE
       returnflow(:) = zero
       reinfiltration(:) = zero
    ENDIF

    !
    ! Compute the net irrigation requirement from Univ of Kassel
    !
    ! This is a very low priority process and thus only applies if
    ! there is some water left in the reservoirs after all other things.
    !
!> The computation of the irrigation is performed here.
!> * First step
!> In a first time, the water requirements (irrig_netereq) by the crops for their optimal growth are calculated
!> over each irrigated fraction (irrigated(ig)/totarea(ig)). It is the difference
!> between the maximal water loss by the crops (transpot_mean) and the net water amount kept by the soil
!> (precipitation and reinfiltration). Transpot_mean is computed in the routines enerbil and diffuco. It
!> is derived from the effective transpiration parametrization under stress-free conditions, called potential transpiration.
!> Crop_coef was used by a previous parametrization of irrigation in the code. Here, its value is equal to one.
!> The crop coefficient was constant in space and time to represent a mean resistance of the vegetation to the potential evaporation.
!> Now, the term crop_coef*Epot is substituted by transpot_mean (see Guimberteau et al., 2011).
!> * Second step
!> We compute irrigation needs in order to supply Irrig_netereq. Water for irrigation (irrig_actual) is withdrawn
!> from the reservoirs. The amount of water is withdrawn in priority from the stream reservoir.
!> If the irrigation requirement is higher than the water availability of the reservoir, water is withdrawn
!> from the fast reservoir or, in the extreme case, from the slow reservoir.
!> * Third step
!> We compute a deficit in water for irrigation. If it is positive, irrigation (depending on water availibility in the reservoirs)
!> has not supplied the crops requirements.
!
    IF ( do_irrigation ) THEN
       DO ig=1,nbpt
          !
          IF ((vegtot(ig) .GT. min_sechiba) .AND. (humrel(ig) .LT. un-min_sechiba) .AND. &
               & (runoff(ig) .LT. min_sechiba) ) THEN
             
             irrig_netereq(ig) = (irrigated(ig) / totarea(ig) ) * MAX(zero, transpot_mean(ig) - &
                  & (precip(ig)+reinfiltration(ig)) )
             
          ENDIF
          !
          DO ib=1,nbasmax
             IF ( routing_area(ig,ib) .GT. 0 ) THEN
             
                irrig_needs(ig,ib) = irrig_netereq(ig) * routing_area(ig,ib)

                irrig_actual(ig,ib) = MIN(irrig_needs(ig,ib),&
                     &   stream_reservoir(ig,ib) + fast_reservoir(ig,ib) + slow_reservoir(ig,ib) )
                
                slow_reservoir(ig,ib) = MAX(zero, slow_reservoir(ig,ib) + &
                     & MIN(zero, fast_reservoir(ig,ib) + MIN(zero, stream_reservoir(ig,ib)-irrig_actual(ig,ib))))

                fast_reservoir(ig,ib) = MAX( zero, &
                     &  fast_reservoir(ig,ib) + MIN(zero, stream_reservoir(ig,ib)-irrig_actual(ig,ib)))

                stream_reservoir(ig,ib) = MAX(zero, stream_reservoir(ig,ib)-irrig_actual(ig,ib) )

                irrig_deficit(ig,ib) = irrig_needs(ig,ib)-irrig_actual(ig,ib)

             ENDIF
          ENDDO
          !
          ! Check if we cannot find the missing water in another basin of the same grid (stream reservoir only).
          ! If we find that then we create some adduction from that subbasin to the one where we need it for
          ! irrigation.
          !
!> If crops water requirements have not been supplied (irrig_deficit>0), we check if we cannot find the missing water
!> in another basin of the same grid. If there is water in the stream reservoir of this subbasin, we create some adduction
!> from that subbasin to the one where we need it for irrigation.
!> 
          DO ib=1,nbasmax

             stream_tot = SUM(stream_reservoir(ig,:))

             DO WHILE ( irrig_deficit(ig,ib) > min_sechiba .AND. stream_tot > min_sechiba)
                
                fi = MAXLOC(stream_reservoir(ig,:))
                ib2 = fi(1)

                irrig_adduct(ig,ib) = MIN(irrig_deficit(ig,ib), stream_reservoir(ig,ib2))
                stream_reservoir(ig,ib2) = stream_reservoir(ig,ib2)-irrig_adduct(ig,ib)
                irrig_deficit(ig,ib) = irrig_deficit(ig,ib)-irrig_adduct(ig,ib)
             
                stream_tot = SUM(stream_reservoir(ig,:))
                
             ENDDO
             
          ENDDO
          !
       ENDDO
       !
       ! If we are at higher resolution we might need to look at neighboring grid boxes to find the streams
       ! which can feed irrigation
!
!> At higher resolution (grid box smaller than 100x100km), we can import water from neighboring grid boxes
!> to the one where we need it for irrigation.
       !
       IF (is_root_prc) THEN
          ALLOCATE(irrig_deficit_glo(nbp_glo, nbasmax), stream_reservoir_glo(nbp_glo, nbasmax), &
               &        irrig_adduct_glo(nbp_glo, nbasmax), stat=ier)
       ELSE
          ALLOCATE(irrig_deficit_glo(0, 0), stream_reservoir_glo(0, 0), &
               &        irrig_adduct_glo(0, 0), stat=ier)
       ENDIF
       IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_flow','Pb in allocate for irrig_deficit_glo, stream_reservoir_glo,...','','')

       CALL gather(irrig_deficit, irrig_deficit_glo)
       CALL gather(stream_reservoir,  stream_reservoir_glo)
       CALL gather(irrig_adduct, irrig_adduct_glo)

       IF (is_root_prc) THEN
          !
          DO ig=1,nbp_glo
             ! Only work if the grid box is smaller than 100x100km. Else the piplines we build
             ! here would be too long to be reasonable.
             IF ( resolution_g(ig,1) < 100000. .AND. resolution_g(ig,2) < 100000. ) THEN
                DO ib=1,nbasmax
                   !
                   IF ( irrig_deficit_glo(ig,ib)  > min_sechiba ) THEN
                      !
                      streams_around(:,:) = zero
                      !
                      DO in=1,NbNeighb
                         ig2 = neighbours_g(ig,in)
                         IF (ig2 .GT. 0 ) THEN
                            streams_around(in,:) = stream_reservoir_glo(ig2,:)
                            igrd(in) = ig2
                         ENDIF
                      ENDDO
                      !
                      IF ( MAXVAL(streams_around) .GT. zero ) THEN
                         !
                         ff=MAXLOC(streams_around)
                         ig2=igrd(ff(1))
                         ib2=ff(2)
                         !
                         IF ( routing_area_glo(ig2,ib2) .GT. 0 .AND. stream_reservoir_glo(ig2,ib2) > zero ) THEN
                            adduction = MIN(irrig_deficit_glo(ig,ib), stream_reservoir_glo(ig2,ib2))
                            stream_reservoir_glo(ig2,ib2) = stream_reservoir_glo(ig2,ib2) - adduction
                            irrig_deficit_glo(ig,ib) = irrig_deficit_glo(ig,ib) - adduction
                            irrig_adduct_glo(ig,ib) = irrig_adduct_glo(ig,ib) + adduction
                         ENDIF
                         !
                      ENDIF
                      !
                   ENDIF
                   !
                ENDDO
             ENDIF
          ENDDO
          !
       ENDIF
       !

       CALL scatter(irrig_deficit_glo, irrig_deficit)
       CALL scatter(stream_reservoir_glo,  stream_reservoir)
       CALL scatter(irrig_adduct_glo, irrig_adduct)

       DEALLOCATE(irrig_deficit_glo, stream_reservoir_glo, irrig_adduct_glo)

    ENDIF

    !! Calculate the net water flow to each routing reservoir (in kg/dt)
    !! to further diagnose the corresponding water budget residu
    !! in routing_subgrid_htu_main

    netflow_fast_diag(:) = zero
    netflow_slow_diag(:) = zero
    netflow_stream_diag(:) = zero

    DO ib=1,nbasmax
       DO ig=1,nbpt
          netflow_fast_diag(ig) = netflow_fast_diag(ig) + runoff(ig)*routing_area(ig,ib) &
               - fast_flow(ig,ib) - pond_inflow(ig,ib)
          netflow_slow_diag(ig) = netflow_slow_diag(ig) + drainage(ig)*routing_area(ig,ib) &
               - slow_flow(ig,ib)
          netflow_stream_diag(ig) = netflow_stream_diag(ig) + flood_flow(ig,ib) + transport(ig,ib) &
               - stream_flow(ig,ib) - return_swamp(ig,ib) - floods(ig,ib)
       ENDDO
    ENDDO

    !! Grid cell averaging
    DO ig=1,nbpt
       netflow_fast_diag(ig) = netflow_fast_diag(ig)/totarea(ig)
       netflow_slow_diag(ig) = netflow_slow_diag(ig)/totarea(ig)
       netflow_stream_diag(ig) = netflow_stream_diag(ig)/totarea(ig)
    ENDDO

    !
    !
    ! Compute the fluxes which leave the routing scheme
    !
    ! Lakeinflow is in Kg/dt
    ! returnflow is in Kg/m^2/dt
    !
    hydrographs(:) = zero
    hydrotemp(:) = zero
    HTUhgmon(:,:) = zero
    HTUtempmon(:,:) = zero
    slowflow_diag(:) = zero
    fast_diag(:) = zero
    slow_diag(:) = zero
    stream_diag(:) = zero
    flood_diag(:) =  zero
    flood_height_diag(:) =  zero
    pond_diag(:) =  zero
    irrigation(:) = zero
    !
    !
    DO ib=1,nbasmax
       !
       DO ig=1,nbpt
          !
          DO im=1,nbasmon
             IF (HTUdiag_loc(ig,im) > 0 .AND. HTUdiag_loc(ig,im) .EQ. ib ) THEN
                HTUhgmon(ig,im) = fast_flow(ig,ib) + slow_flow(ig,ib) + stream_flow(ig,ib)
                HTUtempmon(ig,im) = stream_temp(ig,ib)
             ENDIF
          ENDDO
          !
          IF (hydrodiag(ig) == ib) THEN
             hydrographs(ig) = fast_flow(ig,ib) + slow_flow(ig,ib) + stream_flow(ig,ib)
             hydrotemp(ig) = stream_temp(ig,ib)
             slowflow_diag(ig) = slowflow_diag(ig) + slow_flow(ig,ib)
          ENDIF
          fast_diag(ig) = fast_diag(ig) + fast_reservoir(ig,ib)
          slow_diag(ig) = slow_diag(ig) + slow_reservoir(ig,ib)
          stream_diag(ig) = stream_diag(ig) + stream_reservoir(ig,ib)
          IF (floodplains(ig,ib) .GT. min_sechiba) THEN
             flood_diag(ig) = flood_diag(ig) + flood_reservoir(ig,ib)
             flood_height_diag(ig) = flood_height_diag(ig) + flood_height(ig,ib)*routing_area(ig,ib)
          ENDIF
          irrigation (ig) = irrigation (ig) + irrig_actual(ig,ib) + irrig_adduct(ig,ib)
       ENDDO
    ENDDO
    !
    DO ig=1,nbpt
       fast_diag(ig) = fast_diag(ig)/totarea(ig)
       slow_diag(ig) = slow_diag(ig)/totarea(ig)
       stream_diag(ig) = stream_diag(ig)/totarea(ig)
       flood_diag(ig) = flood_diag(ig)/totarea(ig)
       flood_height_diag(ig) = flood_height_diag(ig)/totarea(ig)
       pond_diag(ig) = pond_reservoir(ig)/totarea(ig)
       !
       irrigation(ig) = irrigation(ig)/totarea(ig)
       !
       ! The three output types for the routing : endorheic basins, rivers and 
       ! diffuse coastal flow.
       !
       lakeinflow(ig) = transport(ig,nbasmax+1)
       coastalflow(ig) = transport(ig,nbasmax+2)
       riverflow(ig) = transport(ig,nbasmax+3)
       !
    ENDDO
    !
    flood_res = flood_diag + pond_diag
    

    !! Remove water from lake reservoir if it exceeds the maximum limit and distribute it 
    !! uniformly over all possible the coastflow gridcells
    
    ! Calculate lake_overflow and remove it from lake_reservoir
    DO ig=1,nbpt
       lake_overflow(ig) = MAX(0., lake_reservoir(ig) - max_lake_reservoir*totarea(ig))
       lake_reservoir(ig) = lake_reservoir(ig) - lake_overflow(ig)
    END DO
    ! Transform lake_overflow from kg/grid-cell/dt_routing into kg/m^2/s
    CALL xios_orchidee_send_field("lake_overflow",lake_overflow(:)/totarea(:)/dt_routing)

    ! Calculate the sum of the lake_overflow and distribute it uniformly over all gridboxes
    CALL gather(lake_overflow,lake_overflow_g)
    IF (is_root_prc) THEN
       total_lake_overflow=SUM(lake_overflow_g)
    END IF
    CALL bcast(total_lake_overflow)

    ! Distribute the lake_overflow uniformly over all coastal gridcells
    ! lake_overflow_coast is only calculated to be used as diagnostics if needed
    DO ig=1,nbpt
       coastalflow(ig) = coastalflow(ig) + total_lake_overflow/nb_coast_gridcells * mask_coast(ig)
       lake_overflow_coast(ig) = total_lake_overflow/nb_coast_gridcells * mask_coast(ig)
    END DO
    ! Transform from kg/grid-cell/dt_routing into m^3/grid-cell/s to match output unit of coastalflow
    CALL xios_orchidee_send_field("lake_overflow_coast",lake_overflow_coast/mille/dt_routing)
   

  END SUBROUTINE routing_hr_flow
  !
!! ================================================================================================================================
!! SUBROUTINE 	: groundwatertemp
!!
!>\BRIEF        : This subroutine computes the temperature of the groundwater leaving the HTU
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
    !-
  SUBROUTINE groundwatertemp(nbpt, nbasmax, nl, tempdiag, lev, dlz, fast_temp, slow_temp)
    ! INPUT
    INTEGER(i_std), INTENT(in)                   :: nbpt, nbasmax, nl
    REAL(r_std), INTENT(in)                      :: tempdiag(nbpt,nl)
    REAL(r_std), INTENT(in)                      :: lev(nl), dlz(nl)
    REAL(r_std), INTENT(inout)                   :: slow_temp(nbpt,nbasmax), fast_temp(nbpt,nbasmax)
    ! OUTPUT
    ! LOCAL
    INTEGER(i_std)                               :: ig, ib, im
    REAL(r_std)                                  :: sw
    REAL(r_std)                                  :: rw(nl), dw(nl)
    LOGICAL, SAVE                                :: alltop=.FALSE.
    LOGICAL, SAVE                                :: FirstCall=.TRUE.
    !
    IF ( FirstCall ) THEN
       !Config Key   = ROUTING_ALLTOPT
       !Config Desc  = Should drainage have the temperature of the top soil (0.3m) ?
       !Config Def   = False
       !Config Help  = The default behaviour of the scheme is that runoff has the temperature
       !Config Help    of the top 30 cm of soil. Drainage will have the temperature of the lowest
       !Config Help    soil layer (3-17m). If set to True this flag will give drainage the same
       !Config Help    temperature as runoff.
       !Config Units = Logical
       alltop=.FALSE.
       CALL getin_p('ROUTING_ALLTOPT', alltop)
       !
       WRITE(numout,*) "Runoff will have the average soil temperature of layers from ", runofftempdepth(1),&
            &          " to ", runofftempdepth(2), "[m]"
       !
       IF ( alltop ) THEN
          WRITE(numout,*) "Drainage will have the average soil temperature of layers from ", runofftempdepth(1),&
            &          " to ", runofftempdepth(2), "[m]"
       ELSE
          WRITE(numout,*) "Drainage will have the average soil temperature of layers from ", drainagetempdepth(1),&
               &          " to ", MIN(drainagetempdepth(2), SUM(dlz)), "[m]"
       ENDIF
       FirstCall=.FALSE.
    ENDIF
    !
    CALL tempdepthweight(nl, dlz, runofftempdepth(1), runofftempdepth(2), rw)
    CALL tempdepthweight(nl, dlz, drainagetempdepth(1), MIN(drainagetempdepth(2), SUM(dlz)), dw)
    !
    slow_temp(:,:) = zero
    fast_temp(:,:) = zero
    ! Compute for each HTU the temperature of runoff and drainage water.
    DO im = 1,nl
       DO ib=1,nbasmax
          DO ig=1,nbpt
             fast_temp(ig,ib) = fast_temp(ig,ib) + tempdiag(ig,im)*rw(im)
             ! The option to have drainage water at the same temperature as runoff
             IF ( alltop ) THEN
                slow_temp(ig,ib) = slow_temp(ig,ib) + tempdiag(ig,im)*rw(im)
             ELSE
                slow_temp(ig,ib) = slow_temp(ig,ib) + tempdiag(ig,im)*dw(im)
             ENDIF
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE groundwatertemp

  SUBROUTINE tempdepthweight(n, dz, top, bot, w)
    ! Input
    INTEGER(i_std), INTENT(in)                   :: n
    REAL(r_std), INTENT(in)                      :: dz(n)
    REAL(r_std), INTENT(in)                      :: top, bot
    ! Output
    REAL(r_std), INTENT(out)                     :: w(n)
    ! Local
    INTEGER(i_std)                               :: i
    REAL(r_std)                                  :: sw
    w(:) = zero
    sw = zero
    DO i=1,n
       w(i) = MAX(zero, MIN(sw+dz(i), bot) - MAX(top, sw))
       sw = sw + dz(i)
    ENDDO
    w(:) = w(:)/(bot-top)
  END SUBROUTINE tempdepthweight

!! ================================================================================================================================
!! SUBROUTINE 	: downstreamsum
!!
!>\BRIEF        : This subroutine sums the input variables onto the downstream HTU in the river graph.
!!
!! DESCRIPTION  : We assume that the downstream HTU is defined by route_togrid and route_tobas. As these
!!                donwstream HTU can be on another processor we do this job on the root processor. So before we need to
!!                transfer all the data onto that processor and then redistribute the result.
!!                Keep in mind that if an HTU does not exit then route_tobas = 0. So the result array needs
!!                to have this index. The end of the rivers are between nbmax+1 and nbmax+3 so this indexing space is also
!!                needed in the result array.
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
    !-
  SUBROUTINE downstreamsum(nbpt, nbmax, v, t)
    ! Input
    INTEGER(i_std), INTENT(in)                           :: nbpt, nbmax
    REAL(r_std), INTENT(in), DIMENSION(nbpt, nbmax)      :: v
    ! Output
    REAL(r_std), INTENT(out), DIMENSION(nbpt, 0:nbmax+3) :: t
    !
    ! Local
    !
    INTEGER(i_std)                                       :: ig, ib, rtg, rtb
    INTEGER(i_std)                                       :: ier
    REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)       :: v_g, t_g
    !
    ! Allocate memory if needed. Should only happen only once in order to reduce computing time.
    !
    IF ( .NOT. ALLOCATED(v_g) ) THEN
       IF (is_root_prc)  THEN
          ALLOCATE(v_g(nbp_glo,nbmax), stat=ier)
          IF (ier /= 0) CALL ipslerr_p(3,'downstreamsum','Pb in allocate for v_g','','')
       ELSE
          ALLOCATE(v_g(1,1))
       ENDIF
    ENDIF
    IF ( .NOT. ALLOCATED(t_g) ) THEN
       IF (is_root_prc)  THEN
          ALLOCATE(t_g(nbp_glo,0:nbmax+3), stat=ier)
          IF (ier /= 0) CALL ipslerr_p(3,'downstreamsum','Pb in allocate for t_g','','')
       ELSE
          ALLOCATE(t_g(1,1))
       ENDIF
    ENDIF
    !
    ! Gather the source variable on the root processor.
    !
    CALL gather(v, v_g)
    !
    ! The downstream sum is performed only on the root processor.
    !
    IF (is_root_prc) THEN
       t_g(:,:) = zero
       DO ib=1,nbmax
          DO ig=1,nbp_glo
             rtg = route_togrid_glo(ig,ib)
             rtb = route_tobasin_glo(ig,ib)
             t_g(rtg,rtb) = t_g(rtg,rtb) + v_g(ig,ib)
          ENDDO
       ENDDO
    ENDIF
    !
    ! Redistribute the downstream field to the all processors.
    !
    CALL scatter(t_g, t)
    !
  END SUBROUTINE downstreamsum
!! ================================================================================================================================
!! SUBROUTINE 	: routing_hr_flood
!!
!>\BRIEF        : This subroutine estimate the flood fraction and the flood height for each HTU 
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
    !-
  SUBROUTINE routing_hr_flood(nbpt, flood_frac, totarea, totflood)
    !
   IMPLICIT NONE
   !
   !! INPUT VARIABLES
   INTEGER(i_std), INTENT(in)                   :: nbpt                      !! Domain size (unitless)
   REAL(r_std), INTENT(in), DIMENSION(nbpt)                 :: totflood                  !! Total amount of water in the floodplains reservoir (kg)
   REAL(r_std), INTENT(in), DIMENSION(nbpt)                 :: totarea                   !! Total area of basin (m^2)
   !! Flooded fraction of the grid box (unitless;0-1)
   !
   !! OUTPUT VARIABLES
   REAL(r_std), INTENT(inout)                   :: flood_frac(nbpt)
   
   !
   !! LOCAL VARIABLES
   INTEGER(i_std)                               :: ig, ib                    !! Indices (unitless)
   REAL(r_std)                                  :: diff, voltemp             !! Discharge reduction due to floodplains   
   !_ ================================================================================================================================
   !
   ! 
   ! Initialize the variables
   flood_frac(:) = zero
   flood_height(:,:) = zero
   flood_frac_bas(:,:) = zero
   DO ig=1, nbpt
      IF (totflood(ig) .GT. min_sechiba) THEN
         DO ib=1,nbasmax
            IF (floodplains(ig,ib) .GT. min_sechiba) THEN
              ! We have to convert h0 to m and the flood_reservoir in m^3  
               flood_frac_bas(ig,ib) = ((fp_beta(ig,ib)+un) * flood_reservoir(ig,ib) / 1000) / ( floodcri(ig,ib) / 1000 * floodplains(ig,ib))
               flood_frac_bas(ig,ib) = (flood_frac_bas(ig,ib)) ** (fp_beta(ig,ib)/(fp_beta(ig,ib)+1))
               flood_frac_bas(ig,ib) = MIN(flood_frac_bas(ig,ib), floodplains(ig,ib)/ routing_area(ig,ib) )

               ! flood_height is in mm
               ! there is two cases: flood_height < h0, flood_height >= h0 (this corresponds to flood_frac_bas = 1 )
               IF ( flood_frac_bas(ig,ib) .EQ. floodplains(ig,ib) / routing_area(ig,ib) ) THEN
                 ! voltemp is on m^3
                 ! Calculation of volume corresponding to h0
                 voltemp = floodplains(ig,ib)/(fp_beta(ig,ib)+un) * ( floodcri(ig,ib) / 1000 )
                 voltemp = flood_reservoir(ig,ib) / 1000 - voltemp
                 ! flood height is in mm
                 flood_height(ig, ib) = voltemp / floodplains(ig,ib) * 1000 + floodcri(ig,ib)
               ELSE
                 ! flood height is in mm
                 flood_height(ig, ib) = (flood_frac_bas(ig,ib)) ** (1/fp_beta(ig,ib)) * floodcri(ig,ib)
               END IF 
            ENDIF
         ENDDO
       ENDIF
         
       DO ib=1,nbasmax
            flood_frac(ig) = flood_frac(ig) + flood_frac_bas(ig,ib) * routing_area(ig,ib) / totarea(ig)
       END DO
       flood_frac(ig) = flood_frac(ig) + pond_frac(ig)
       !
   ENDDO
   
   END SUBROUTINE routing_hr_flood
   !
!! ================================================================================================================================
!! SUBROUTINE 	: routing_hr_overflow
!!
!>\BRIEF        : This subroutine performs the overflow fluxes
!!
!! DESCRIPTION (definitions, functional, design, flags): \n
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
   !-
   SUBROUTINE routing_hr_overflow(nbpt, nbasmax)
      !
     IMPLICIT NONE
     !
     !! INPUT VARIABLES
     INTEGER(i_std), INTENT(in)                   :: nbpt,nbasmax              !! Domain size (unitless)
     !
     !! LOCAL VARIABLES
     REAL(r_std), DIMENSION(nbpt,nbasmax)         :: transport_overflow        !! Water transport between floodplains - flood overflow (kg/dt)
     REAL(r_std), DIMENSION(nbp_glo,nbasmax)      :: transport_overflow_glo    !! Water transport between floodplains - flood overflow (kg/dt)
     REAL(r_std), DIMENSION(nbpt,nbasmax)         :: overflow_loss             !! Water loss from flood overflow (kg/dt)
     REAL(r_std), DIMENSION(nbp_glo,nbasmax)      :: overflow_loss_glo         !! Water loss from flood overflow (kg/dt)
     !
     REAL(r_std), ALLOCATABLE, DIMENSION(:,:)     :: flood_height_g            !! Floodplains height (m)
     REAL(r_std), ALLOCATABLE, DIMENSION(:,:)     :: flood_frac_bas_g          !! Fraction of the HTU flooded 
     REAL(r_std), ALLOCATABLE, DIMENSION(:,:)     :: flood_reservoir_glo       !! Water amount in the stream reservoir (kg)  
     !
     REAL(r_std), ALLOCATABLE, DIMENSION(:)       :: DH,DH_temp                !! Difference of height - flood overflow (kg/dt)
     !
     INTEGER(i_std)                               :: numflood                  !!
     !
     INTEGER(i_std)                               :: ig, ib, inf,inb,ing       !! Indices (unitless)
     REAL(r_std)                                  :: diff                      !! Discharge reduction due to floodplains   
     REAL(r_std)                                  :: flow                      !! Outflow computation for the reservoirs (kg/dt)
     REAL(r_std)                                  :: dorog                     !! Discharge reduction due to floodplains
     INTEGER(i_std)                               :: ier                       !! Error handling

  
     !_ ================================================================================================================================
     !
     !! ANTHONY : OVERFLOW
     !! CALCULATE TRANSFER BETWEEN FLOODPLAINS RESERVOIR
     IF (is_root_prc)  THEN
        ALLOCATE( flood_height_g(nbp_glo, nbasmax), flood_frac_bas_g(nbp_glo, nbasmax), stat=ier) 
        ALLOCATE( flood_reservoir_glo(nbp_glo, nbasmax), stat=ier) 
     ELSE
        ALLOCATE( flood_height_g(1,1), flood_frac_bas_g(1,1), stat=ier)
        ALLOCATE( flood_reservoir_glo(1, 1), stat=ier) 
     ENDIF
     !
     IF (ier /= 0) CALL ipslerr_p(3,'routing_hr_flow','Pb in allocate for flood_height_glo/floog_frac_glo','','')
     !
     CALL gather(flood_height,flood_height_g)
     CALL gather(flood_frac_bas,flood_frac_bas_g)
     CALL gather(flood_reservoir,flood_reservoir_glo)
     !
     IF (is_root_prc) THEN
        transport_overflow_glo(:,:) = 0
        overflow_loss_glo(:,:) = 0
        DO ib=1,nbasmax
           DO ig=1,nbp_glo
              IF ( floodplains_glo(ig,ib)/routing_area_glo(ig,ib) .GT. 0.5) THEN
                 numflood = 0 ! Number of inflows for overflow
                 ALLOCATE(DH(route_innum_glo(ig,ib)))
                 DH(:) = 0
                 DH_temp(:) = -1
                 DO inf=1,route_innum_glo(ig,ib)
                    ing  = route_ingrid_glo(ig,ib,inf)
                    inb  = route_inbasin_glo(ig,ib,inf)
                    IF ( floodplains_glo(ing,inb)/routing_area_glo(ing,inb) .GT. 0 ) THEN
                       ! Minimum of deltaorog is defined at lim_floodcri (0.3 m
                       ! can be used).
                       dorog = MAX(orog_min_glo(ing,inb)- orog_min_glo(ig,ib), lim_floodcri)
                       ! flood_height is in mm and orog min in m
                       diff = (flood_height_g(ig,ib)- flood_height_g(ing,inb))/1000 - dorog
                       DH(inf) = max(diff, 0.)
                       ! 
                       ! Flux is estimated via floodplains_glo
                       ! Then factor 1000 is to convert m^3 to kg
                       ! OVERFLOW_TCST is in seconds
                       flow = DH(inf) * (floodplains_glo(ig,ib)* floodplains_glo(ing,inb))/(floodplains_glo(ig,ib)+floodplains_glo(ing,inb))*1000 / overflow_tcst * dt_routing / one_day
                       transport_overflow_glo(ing,inb) = transport_overflow_glo(ing,inb) + flow
                       overflow_loss_glo(ig,ib) = overflow_loss_glo(ig,ib) + flow
                    END IF
                 END DO
                 DEALLOCATE(DH)
              END IF
           ENDDO
        ENDDO
     END IF
     ! Send to local variables
     CALL scatter(transport_overflow_glo, transport_overflow)
     CALL scatter(overflow_loss_glo, overflow_loss)
     ! Apply the volume changes
     DO ig=1,nbpt
        DO ib=1,nbasmax
           IF ( floodplains(ig,ib) .GT. 0 ) THEN
                 flood_reservoir(ig,ib) = flood_reservoir(ig,ib) + transport_overflow(ig,ib) - overflow_loss(ig,ib)
                 ! NEED to check if flood reservoir is less than 0, this may be a critical issue
                 ! Solved by an adequate use of an higher overflow time constant 
                 ! To obtain the same result as with a lower overflow parameter 
                 ! -> repeat a few time the operation with and higher overflow parameter
                 IF ( flood_reservoir(ig,ib) .LT. 0 ) THEN

                     WRITE(*,*) "Issue of flood reservoir < 0 due to overflow at ", ig, ib
                     stream_reservoir(ig,ib) =  stream_reservoir(ig,ib) + stream_reservoir(ig,ib) ! + because negative !
                     flood_reservoir(ig,ib) = 0
                 END IF
           END IF
        END DO
     END DO
     DEALLOCATE( flood_height_g, flood_frac_bas_g) 
     
     END SUBROUTINE routing_hr_overflow
   !
!! ================================================================================================================================
!! SUBROUTINE 	: routing_hr_lake
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

  SUBROUTINE routing_hr_lake(nbpt, dt_routing, lakeinflow, humrel, return_lakes)
    !
    IMPLICIT NONE
    !
!! INPUT VARIABLES
    INTEGER(i_std), INTENT(in) :: nbpt               !! Domain size (unitless)
    REAL(r_std), INTENT (in)   :: dt_routing         !! Routing time step (s)
    REAL(r_std), INTENT(out)    :: lakeinflow(nbpt)   !! Water inflow to the lakes (kg/dt)
    REAL(r_std), INTENT(in)    :: humrel(nbpt)       !! Soil moisture stress, root extraction potential (unitless)
    !
!! OUTPUT VARIABLES
    REAL(r_std), INTENT(out)   :: return_lakes(nbpt) !! Water from lakes flowing back into soil moisture (kg/m^2/dt)
    !
!! LOCAL VARIABLES
    INTEGER(i_std)             :: ig                 !! Indices (unitless)
    REAL(r_std)                :: refill             !!
    REAL(r_std)                :: total_area         !! Sum of all the surfaces of the basins (m^2)

!_ ================================================================================================================================
    !
    !
    DO ig=1,nbpt
       !
       total_area = SUM(routing_area(ig,:))
       !
       lake_reservoir(ig) = lake_reservoir(ig) + lakeinflow(ig)
       
       IF ( doswamps ) THEN
          ! Calculate a return flow that will be extracted from the lake reservoir and reinserted in the soil in hydrol
          ! Uptake in Kg/dt
          refill = MAX(zero, maxevap_lake * (un - humrel(ig)) * dt_routing * total_area)
          return_lakes(ig) = MIN(refill, lake_reservoir(ig))
          lake_reservoir(ig) = lake_reservoir(ig) - return_lakes(ig)
          ! Return in Kg/m^2/dt
          return_lakes(ig) = return_lakes(ig)/total_area
       ELSE
          return_lakes(ig) = zero
       ENDIF

       ! This is the volume of the lake scaled to the entire grid.
       ! It would be better to scale it to the size of the lake
       ! but this information is not yet available.
       lake_diag(ig) = lake_reservoir(ig)/total_area

       lakeinflow(ig) = lakeinflow(ig)/total_area

    ENDDO
    !
  END SUBROUTINE routing_hr_lake
  !
!! ================================================================================================================================
!! SUBROUTINE 	: routing_hr_basins_p
!!
!>\BRIEF        This routing read the file created by RoutingPreProc : https://gitlab.in2p3.fr/ipsl/lmd/intro/routingpp
!!
!! DESCRIPTION (definitions, functional, design, flags) : None
!!  Once the atmospheric grid is defined and the land/sea mask set, RoutingPreProc has to used to generate the
!!  HTU graphs for the domain. This can be done either on the basis of the HydroSHEDS, MERIT or the old Vörösmarty map
!!  of catchments. During this step all the information will be created to allow ORCHIDEE to route the water and
!!  and monitor the flows at given stations.
!!  For the moment the ROUTING_FILE (Perhaps to renamed RoutingGraph) is read using IOIPSL but that should evolve toward XIOS.
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S):
!!
!! REFERENCES   : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE routing_hr_basins_p(nbpt, lalo, neighbours, resolution, contfrac)
    !
    IMPLICIT NONE
    !
!! INPUT VARIABLES
    INTEGER(i_std), INTENT(in) :: nbpt               !! Domain size (unitless)
    REAL(r_std), INTENT(in)    :: lalo(nbpt,2)       !! Vector of latitude and longitudes (beware of the order !)
    INTEGER(i_std), INTENT(in) :: neighbours(nbpt,NbNeighb) !! Vector of neighbours for each grid point (1=North and then clockwise) (unitless)
    REAL(r_std), INTENT(in)    :: resolution(nbpt,2) !! The size of each grid box in X and Y (m)
    REAL(r_std), INTENT(in)    :: contfrac(nbpt)     !! Fraction of land in each grid box (unitless;0-1)
    !
    ! LOCAL
    !
    INTEGER(i_std)    :: iml, jml, lml, tml
    INTEGER(i_std)    :: i, j, ni, fid, ib, ig, ic, ign, ibn, og, ob, ier, im
    REAL(r_std)       :: corr
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:)    :: tmpvar_glo
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)      :: tmpvar
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)      :: lon, lat, landindex
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:)   :: indextab
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:)   :: landfileindex
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:)     :: land2land
    INTEGER(i_std)                                :: nbhtumon
    REAL(r_std), ALLOCATABLE, DIMENSION(:)        :: floodmap_glo
    !
!_ ================================================================================================================================
    !
    !
    !
    IF (is_root_prc) THEN
       !
       CALL flininfo(graphfilename, iml, jml, lml, tml, fid)
       !
       IF (iml .NE. iim_g .AND. jml .NE. jjm_g ) THEN
          CALL ipslerr(3,'routing_hr_basins_p',&
               'The routing graph file does not have the right dimensions for the model.', &
               'Are you sure you are using the right routing graph file ?', '  ')
       ENDIF
       !
       !
       ALLOCATE(tmpvar_glo(iml,jml,nbasmax))
       ALLOCATE(tmpvar(iml,jml))
       ALLOCATE(lon(iml,jml))
       ALLOCATE(lat(iml,jml))
       ALLOCATE(landindex(iml,jml))
       ALLOCATE(indextab(iml,jml))
       ALLOCATE(landfileindex(iml,jml))
       ALLOCATE(land2land(iml*jml))
       !      
       CALL flinget(fid, 'lon', iml, jml, 1, tml, 1, 0, lon)
       CALL flinget(fid, 'lat', iml, jml, 1, tml, 1, 0, lat)
       CALL flinget(fid, 'nbpt_glo', iml, jml, 1, tml, 1, 0, landindex)
       !
       ! Replace NaN and other undef values
       !
       DO i=1,iml
          DO j=1,jml
             IF ( landindex(i,j) /= landindex(i,j) .OR. landindex(i,j) < 0 .OR. landindex(i,j) >= undef_graphfile) THEN
                landindex(i,j) = -1
             ENDIF
          ENDDO
       ENDDO
       !
       ! Compute land index for file data. Information could be in file !
       !
       ni=NINT(MAXVAL(landindex))
       IF ( ni .NE. nbp_glo) THEN
          WRITE(numout,*) "Error routing_hr_basins_p : ni, nbp_glo : ", ni, nbp_glo, undef_graphfile
          CALL ipslerr(3,'routing_hr_basins_p',&
               'The routing graph file does not have the same number', &
               'of land points as the model.',&
               '  ')
       ENDIF
       !
       CALL routing_hr_indexfilegrid(iml, jml, nbp_glo, lon, lat, landindex, indextab, land2land)
       !
       CALL flinget(fid, 'basin_area', iml, jml, nbasmax, tml, 1, 0, tmpvar_glo)
       CALL routing_hr_landgather(iml, jml, nbasmax, nbp_glo, indextab, tmpvar_glo, routing_area_glo, &
            &                    zero)
       
       IF ( do_floodplains ) THEN
          CALL flinget(fid, 'basin_floodp', iml, jml, nbasmax, tml, 1, 0, tmpvar_glo)
          CALL routing_hr_landgather(iml, jml, nbasmax, nbp_glo, indextab, tmpvar_glo, floodplains_glo, &
            &                    zero)
          !
          CALL flinget(fid, 'floodcri', iml, jml, nbasmax, tml, 1, 0, tmpvar_glo)
          CALL routing_hr_landgather(iml, jml, nbasmax, nbp_glo, indextab, tmpvar_glo, floodcri_glo, &
            &                    un)
          !
          CALL flinget(fid, 'basin_beta_fp', iml, jml, nbasmax, tml, 1, 0, tmpvar_glo)
          CALL routing_hr_landgather(iml, jml, nbasmax, nbp_glo, indextab, tmpvar_glo, fp_beta_glo, &
            &                    un)
       END IF
       
       CALL flinget(fid, 'topoindex', iml, jml, nbasmax, tml, 1, 0, tmpvar_glo)
       CALL routing_hr_landgather(iml, jml, nbasmax, nbp_glo, indextab, tmpvar_glo, topo_resid_glo, &
            &                    undef_graphfile)

       IF ( graphfile_version >= 2.0) THEN
          CALL flinget(fid, 'topoindex_stream', iml, jml, nbasmax, tml, 1, 0, tmpvar_glo)
          CALL routing_hr_landgather(iml, jml, nbasmax, nbp_glo, indextab, tmpvar_glo, stream_resid_glo, &
               &                    undef_graphfile)
          CALL ipslerr(1,'routing_hr_basins_p',&
               'The topoindex_stream variable was found in routing_graph.nc', &
               'It will be used the topographic index of the stream store.',&
               '  ')
       ELSE
          stream_resid_glo(:,:) = topo_resid_glo(:,:)
       ENDIF
       stream_maxresid=MAXVAL(stream_resid_glo, MASK=stream_resid_glo .LT. undef_graphfile)
       
       CALL flinget(fid, 'basinid', iml, jml, nbasmax, tml, 1, 0, tmpvar_glo)
       CALL routing_hr_landgather(iml, jml, nbasmax, nbp_glo, indextab, tmpvar_glo, global_basinid_glo, &
            &                    undef_int)

       CALL flinget(fid, 'routetogrid', iml, jml, nbasmax, tml, 1, 0, tmpvar_glo)
       CALL routing_hr_landgather(iml, jml, nbasmax, nbp_glo, indextab, tmpvar_glo, route_togrid_glo, &
            &                    undef_int)
       CALL routing_hr_convertlandpts(nbp_glo, nbasmax, land2land, route_togrid_glo)
       
       CALL flinget(fid, 'routetobasin', iml, jml, nbasmax, tml, 1, 0, tmpvar_glo)
       CALL routing_hr_landgather(iml, jml, nbasmax, nbp_glo, indextab, tmpvar_glo, route_tobasin_glo, 0)
       
       CALL flinget(fid, 'routenbintobas', iml, jml, nbasmax, tml, 1, 0, tmpvar_glo)
       CALL routing_hr_landgather(iml, jml, nbasmax, nbp_glo, indextab, tmpvar_glo, route_nbintobas_glo, 0)
       
       !!
       IF ( dofloodoverflow ) THEN
          CALL flinget(fid, 'basin_orog_min', iml, jml, nbasmax, tml, 1, 0, tmpvar_glo)
          CALL routing_hr_landgather(iml, jml, nbasmax, nbp_glo, indextab, tmpvar_glo, orog_min_glo, un)
       END IF 
       !!
       IF ( graphfile_version >= 2.6) THEN
          CALL flinget(fid, 'gridrephtu', iml, jml, 1, tml, 1, 0, tmpvar)
          CALL routing_hr_landgather(iml, jml, nbp_glo, indextab, tmpvar, hydrodiag_glo, -1)
       ELSE
          hydrodiag_glo(:) = 1
       ENDIF
       !!
       IF ( MonitoringinGraph ) THEN
          CALL flinget(fid, 'HTUmonitor', iml, jml, nbasmon, tml, 1, 0, tmpvar_glo)
          CALL routing_hr_landgather(iml, jml, nbasmon, nbp_glo, indextab, tmpvar_glo, HTUdiag_glo, -1)
       ELSE
          HTUdiag_glo(:,:) = -1
       ENDIF
       !
       CALL flinclo(fid)
       DEALLOCATE(indextab)
       DEALLOCATE(lon)
       DEALLOCATE(lat)
       DEALLOCATE(tmpvar_glo)
       DEALLOCATE(tmpvar)
       !
       ! Convert floodplains fraction into floodplains surface
       IF ( do_floodplains ) THEN
          !floodplains_glo(:, :) = 0
          DO ig = 1,nbp_glo
              DO ib = 1,nbasmax
                 floodplains_glo(ig, ib) = routing_area_glo(ig,ib) * floodplains_glo(ig, ib)
                 IF ( floodplains_glo(ig,ib) .GT. min_sechiba ) THEN
                    IF (floodcri_glo(ig,ib) .LT. min_sechiba ) THEN
                       WRITE(*,*) "routing_hr_basins_p : Have problem with floodcri : ", floodcri_glo(ig,ib), " on ", ig, ib
                       WRITE(*,*) "routing_hr_basins_p setting : floodcri_glo(ig,ib) = 100. [mm]"
                       ! The default value is 100mm
                       floodcri_glo(ig,ib) = 100.
                    ENDIF
                 ENDIF
              END DO
           END DO
           ! Compute global floodmap
           ALLOCATE(floodmap_glo(nbp_glo))
           floodmap_glo(:) = 0.0
           DO ig = 1,nbp_glo
              floodmap_glo(ig) = SUM(floodplains_glo(ig,:))
           ENDDO
       END IF
       !
       ! Verifications of the routing graph.
       !
       nbhtumon = 0
       DO ig = 1,nbp_glo
          ! Noramlize the areas so that differences in precision of area compution by RoutingPP do not affect the model
          !
          corr = contfrac_g(ig)*area_g(ig)/SUM(routing_area_glo(ig,:))
          IF (ABS(1 - corr) > 0.0002 ) THEN
             WRITE(numout ,*) "Correcting the HTU area to take into account contfrac", corr
             IF ( ABS(1 - corr) > 0.1) THEN
                WRITE(numout,*) "Coordinates : ", lalo_g(ig,1), lalo_g(ig,2)
                WRITE(numout,*) "Contfrac and area in model : ", contfrac_g(ig), area_g(ig)
                WRITE(numout,*) "Total grid area in graph file : ", SUM(routing_area_glo(ig,:))
                WRITE(numout,*) "The new areas are : ", SUM(routing_area_glo(ig,:)), contfrac_g(ig)*area_g(ig)
                WRITE(numout,*) "Correction factor : ", corr
                CALL ipslerr(3,'routing_hr_basins_p',&
                     'There is a mismatch in the  area of the grid', &
                     'Either there are issues with the projection of the grid ',' or contfrac mismatches.')
             ELSE
                CALL ipslerr(2,'routing_hr_basins_p',&
                     'The area of the grid had to be adjusted by less than 10% :', &
                     ' ','  ')
             ENDIF
          ENDIF
          DO ib = 1,nbasmax
             routing_area_glo(ig,ib) = corr*routing_area_glo(ig,ib)
          ENDDO
          !     
          !
          DO ib = 1,nbasmax
             !
             IF (topo_resid_glo(ig,ib) <= zero .AND. route_tobasin_glo(ig, ib) .LE. nbasmax+3) THEN
                ! If the basin has no surface we change silently as it does not matter.
                IF ( routing_area_glo(ig,ib) > zero ) THEN
                   CALL ipslerr(2,'routing_hr_basins_p',&
                        'Some zero topo_resid (topoindex) values were encoutered and replaced here :', &
                        ' ','  ')
                   WRITE(*,*) "routing_hr_basins_p : topo_resid_glo : ", topo_resid_glo(ig,ib), routing_area_glo(ig,ib)
                   WRITE(*,*) "routing_hr_basins_p : Coordinates : ", lalo_g(ig,1), lalo_g(ig,2)
                   topo_resid_glo(ig,ib) = 10
                   stream_resid_glo(ig,ib) = 10
                   WRITE(*,*) "routing_hr_basins_p : New topo_resid_glo : ", topo_resid_glo(ig,ib)
                ELSE
                   topo_resid_glo(ig,ib) = 10
                   stream_resid_glo(ig,ib) = 10
                ENDIF
             ENDIF
             !
             !
             IF ( route_togrid_glo(ig, ib) > nbp_glo ) THEN
                IF ( route_tobasin_glo(ig,ib) <= nbasmax+3 ) THEN
                   WRITE(*,*) "Issues with the global grid : ", ig, ib, route_togrid_glo(ig, ib), route_tobasin_glo(ig,ib)
                   CALL ipslerr(3,'routing_hr_basins_p','route_togrid is not compatible with the model configuration', &
                        ' ','  ')
                ENDIF
             ELSE
                ic = 0
                ign = ig
                ibn = ib
                ! Locate outflow point
                DO WHILE (ibn .GT. 0 .AND. ibn .LE. nbasmax .AND. ic .LT. nbasmax*nbp_glo)
                   ic = ic + 1
                   og = ign
                   ob = ibn
                   ign = route_togrid_glo(og, ob)
                   ibn = route_tobasin_glo(og, ob)
                   !
                   IF (ibn .GT. nbasmax+3 .OR. ign .GT. nbp_glo) THEN
                      WRITE(*,*) "Reached point ", ign, ibn, " on condition ", nbasmax+3, nbp_glo
                      WRITE(*,*) "Why do we flow into basin :",  route_tobasin_glo(og, ob), " at ", og,ob
                      WRITE(*,*) "Coordinates : ", lalo_g(ob,1), lalo_g(ob,2)
                      WRITE(*,*) "neighbours_g : ", MINVAL(neighbours_g(ob,:))
                      CALL ipslerr(3,'routing_hr_basins_p','The river flows into a place outside of the grid.', &
                           ' ','  ')
                   ENDIF
                ENDDO
                IF ( ic .GE. nbasmax*nbp_glo) THEN
                   WRITE(*,*) "Some river did not converge on point ", ig, ib, ic
                   WRITE(*,*) "The start point in the graph was : ", lalo_g(ig,2), lalo_g(ig,1), ib
                   WRITE(*,*) "The last point we passed through was : ", lalo_g(og,2), lalo_g(og,1), ob
                   WRITE(*,*) "The next one would be : ", lalo_g(ign,2), lalo_g(ign,1), ibn
                   og = route_togrid_glo(ign, ibn)
                   ob = route_tobasin_glo(ign, ibn)
                   WRITE(*,*) "The after next HTU would be : ", lalo_g(og,2), lalo_g(og,1), ob
                   WRITE(*,*) "Last information : ", ign, ibn
                   CALL ipslerr(3,'routing_hr_basins_p','The river never flows into an outflow point.', &
                        ' ','  ')
                ENDIF
             ENDIF
          ENDDO
          !
          ! Count stations to be monitored
          !
          DO im=1,nbasmon
             IF ( HTUdiag_glo(ig,im) > 0 ) THEN
                nbhtumon = nbhtumon + 1
             ENDIF
          ENDDO
       ENDDO
       WRITE(numout,*) "Found a total of ", nbhtumon, " HTUs to be monitored and written into HTUhgmon"
       !
       ! Compute num_largest
       !
       num_largest = COUNT(route_tobasin_glo .EQ. nbasmax+3)
       WRITE(numout,*) "After _basins_p : Number of largest rivers : ", COUNT(route_tobasin_glo .EQ. nbasmax+3)
    ENDIF
    !
    CALL bcast(num_largest)
    CALL bcast(nbasmax)
    CALL bcast(nbasmon)
    CALL bcast(inflows)
    !
    CALL scatter(routing_area_glo,routing_area_loc)
    IF ( do_floodplains ) THEN
       CALL scatter(floodplains_glo,floodplains_loc)
       CALL scatter(floodcri_glo, floodcri_loc)
       CALL scatter(fp_beta_glo, fp_beta_loc)
    END IF
    CALL scatter(global_basinid_glo, global_basinid_loc)
    CALL scatter(topo_resid_glo, topo_resid_loc)
    CALL scatter(stream_resid_glo, stream_resid_loc)
    CALL scatter(route_togrid_glo, route_togrid_loc)
    CALL scatter(route_tobasin_glo, route_tobasin_loc)
    CALL scatter(route_nbintobas_glo, route_nbintobas_loc)
    CALL scatter(hydrodiag_glo, hydrodiag_loc)
    CALL scatter(HTUdiag_glo, HTUdiag_loc)
    IF ( do_floodplains .AND. dofloodoverflow ) THEN
       CALL scatter(orog_min_glo, orog_min_loc) 
    END IF
    !
    CALL bcast(stream_tcst)
    CALL bcast(fast_tcst)
    CALL bcast(slow_tcst)
    CALL bcast(flood_tcst)
    CALL bcast(swamp_cst)
    CALL bcast(lim_floodcri)
    CALL bcast(stream_maxresid)
    !    
  END SUBROUTINE routing_hr_basins_p
  !
!! ================================================================================================================================
!! SUBROUTINE 	: routing_hr_graphinfo
!!
!>\BRIEF Extract some basic information from the routing graph file which cannot be obtained through IOIPSL.
!!
!! ================================================================================================================================
  SUBROUTINE routing_hr_graphinfo(filename, basmax, infmax, basmon, undef, tstream, tfast, tslow, tflood, cswamp, lfpcri)
    !
    USE netcdf
    !
    IMPLICIT NONE
    !
    !! 0. Variables and parameter declaration
    !! 0.1 Input variables
    CHARACTER(LEN=*), INTENT(in)                         :: filename      !! filename: name of the file to open
    INTEGER(i_std), INTENT(inout)                        :: basmax        !! maximum number of HTUs
    INTEGER(i_std), INTENT(inout)                        :: basmon        !! Number of HTUs to be monitored by grid box.
    INTEGER(i_std), INTENT(inout)                        :: infmax        !! Maximum number of inflows.
    REAL(r_std), INTENT(out)                             :: undef
    REAL(r_std), INTENT(out)                             :: tstream, tfast, tslow, tflood, cswamp !! Time constants to be extracted
    REAL(r_std), INTENT(out)                             :: lfpcri !! Constant lim_floodcri to be taken from graph file.
    !
    INTEGER(i_std)                                       :: rcode, nid, dimid, ndims, nvars
    INTEGER(i_std), DIMENSION(4)                         :: dimids
    INTEGER(i_std)                                       :: iv, ndimsvar
    CHARACTER(LEN=20)                                    :: dname, varname
    !
    !
    IF (is_root_prc) THEN
       !
       rcode = nf90_open(TRIM(filename), NF90_NOWRITE, nid)
       IF (rcode == NF90_NOERR) THEN
          !
          ! Get graph file version
          !
          rcode = nf90_get_att(nid, NF90_GLOBAL, "RoutingPPVersion", graphfile_version)
          IF (rcode /= NF90_NOERR) THEN
             graphfile_version = 0.0
          ENDIF
          !
          ! Assumes that the number of HTUs is in the dimension z
          !
          rcode = nf90_inq_dimid(nid, "z", dimid)
          IF (rcode /= NF90_NOERR) CALL ipslerr_p(3, 'routing_hr_graphinfo', 'Error in nf90_inq_dimid for z', &
               TRIM(nf90_strerror(rcode)),'')
          rcode = nf90_inquire_dimension(nid, dimid, dname, basmax)
          IF (rcode /= NF90_NOERR) CALL ipslerr_p(3, 'routing_hr_graphinfo', 'Error in nf90_inquire_dimension basmax', &
               TRIM(nf90_strerror(rcode)),'')
          !
          rcode = nf90_inq_dimid(nid, "inflow", dimid)
          IF (rcode /= NF90_NOERR) CALL ipslerr_p(3, 'routing_hr_graphinfo', 'Error in nf90_inq_dimid for inflow', &
               TRIM(nf90_strerror(rcode)),'')
          rcode = nf90_inquire_dimension(nid, dimid, dname, infmax)
          IF (rcode /= NF90_NOERR) CALL ipslerr_p(3, 'routing_hr_graphinfo', 'Error in nf90_inquire_dimension inflows', &
               TRIM(nf90_strerror(rcode)),'')
          !
          rcode = nf90_inq_dimid(nid, "htumon", dimid)
          IF (rcode /= NF90_NOERR) THEN
             MonitoringinGraph = .FALSE.
             basmon = 1
          ELSE
             rcode = nf90_inquire_dimension(nid, dimid, dname, basmon)
             IF (rcode /= NF90_NOERR) CALL ipslerr_p(3, 'routing_hr_graphinfo', 'Error in nf90_inquire_dimension for basmon', &
                  TRIM(nf90_strerror(rcode)),'')
             MonitoringinGraph = .TRUE.
          ENDIF
          !    
          !
          rcode = NF90_INQUIRE (nid, nDimensions=ndims, nVariables=nvars)
          IF (rcode /= NF90_NOERR) CALL ipslerr_p(3, 'routing_hr_graphinfo', 'Error in nf90_inquire', &
               TRIM(nf90_strerror(rcode)),'')
          !
          DO iv=1,nvars
             !
             rcode = NF90_INQUIRE_VARIABLE(nid, iv, name=varname, ndims=ndimsvar, dimids=dimids)
             !
             SELECT CASE (varname) 
             CASE ("basin_area")
                rcode = NF90_GET_ATT(nid, iv, "missing_value", undef)
                IF (rcode /= NF90_NOERR) THEN
                   IF ( rcode == NF90_ENOTATT ) THEN 
                      rcode = NF90_GET_ATT(nid, iv, "_FillValue", undef)
                      IF ( undef /= undef ) THEN
                         CALL ipslerr_p(2, 'routing_hr_graphinfo', "_FillValue Detected as NaN and set to undef_sechiba", "", "")
                         undef = undef_sechiba
                      ENDIF
                      IF (rcode /= NF90_NOERR) CALL ipslerr_p(3, 'routing_hr_graphinfo', 'Did not get FillValue with nf90_get_att', &
                           TRIM(nf90_strerror(rcode)),'')
                   ELSE 
                      IF (rcode /= NF90_NOERR) CALL ipslerr_p(3, 'routing_hr_graphinfo', 'Error in nf90_get_att', &
                           TRIM(nf90_strerror(rcode)),'')
                   ENDIF
                ENDIF
             CASE("StreamTimeCst")
                rcode = NF90_GET_VAR(nid,iv,tstream)
                IF (rcode /= NF90_NOERR) CALL ipslerr_p(3, 'routing_hr_graphinfo', 'Error in nf90_get_var for variable StreamTimeCst', '','')
                ! If in an old version convert 10^3d/km into s/km
                IF (graphfile_version < 1.0) THEN
                   tstream = tstream/1000*one_day
                ENDIF
             CASE("FastTimeCst")
                rcode = NF90_GET_VAR(nid,iv,tfast)
                IF (rcode /= NF90_NOERR) CALL ipslerr_p(3, 'routing_hr_graphinfo', 'Error in nf90_get_var for variable FastTimeCst', '','')
                ! If in an old version convert 10^3d/km into s/km
                IF (graphfile_version < 1.0) THEN
                   tfast = tfast/1000*one_day
                ENDIF
             CASE("SlowTimeCst")
                rcode = NF90_GET_VAR(nid,iv,tslow)
                IF (rcode /= NF90_NOERR) CALL ipslerr_p(3, 'routing_hr_graphinfo', 'Error in nf90_get_var for variable SlowTimeCst', '','')
                ! If in an old version convert 10^3d/km into s/km
                IF (graphfile_version < 1.0) THEN
                   tslow = tslow/1000*one_day
                ENDIF
             CASE("FloodTimeCst")
                rcode = NF90_GET_VAR(nid,iv,tflood)
                IF (rcode /= NF90_NOERR) CALL ipslerr_p(3, 'routing_hr_graphinfo', 'Error in nf90_get_var for variable FloodTimeCst', '','')
                ! If in an old version convert 10^3d/km into s/km
                IF (graphfile_version < 1.0) THEN
                   tflood = tflood/1000*one_day
                ENDIF
             CASE("SwampCst")
                rcode = NF90_GET_VAR (nid,iv,cswamp)
                IF (rcode /= NF90_NOERR) CALL ipslerr_p(3, 'routing_hr_graphinfo', 'Error in nf90_get_var for variable SwampCst', '','')
             CASE("MaxTimeStep")
                rcode = NF90_GET_VAR (nid,iv,maxtimestep)
                IF (rcode /= NF90_NOERR) CALL ipslerr_p(3, 'routing_hr_graphinfo', 'Error in nf90_get_var for variable MaxTimeStep', '','')
             CASE("LimFloodcri")
                rcode = NF90_GET_VAR(nid,iv,lfpcri)
                IF (rcode /= NF90_NOERR) CALL ipslerr_p(3, 'routing_hr_graphinfo', 'Error in nf90_get_var for variable LimFloodcri', '','')
             END SELECT
          ENDDO
          rcode = NF90_CLOSE(nid)
          IF (rcode /= NF90_NOERR) CALL ipslerr_p(3, 'routing_hr_graphinfo', 'Error in nf90_close', &
               TRIM(nf90_strerror(rcode)),'')
          !
          ! Before RoutingGraph version 2.5 the lfpcri parameter was hardcoded in routing.f90 and set to 2m.
          ! This information is preserved here.
          IF ( graphfile_version < 2.5 ) THEN
             lfpcri = 2.0
          ENDIF
          !
       ELSE
          ! Case without Graphfile. So we consider that the information will be found in the restart.
          CALL ipslerr_p(2, 'routing_hr_graphinfo', 'Could not open the routing_graph.nc file', &
               "Expect to find all the information needed in the restart file.",'')
          !
          MonitoringinGraph = .FALSE.
          basmax = -1
          basmon = -1
          infmax = -1
          undef = undef_sechiba
          tstream = -1
          tfast = -1
          tslow = -1
          tflood = -1
          cswamp = -1
          lfpcri = -1
       ENDIF
    ENDIF
    !!
    CALL bcast(MonitoringinGraph)
    CALL bcast(basmax)
    CALL bcast(basmon)
    CALL bcast(infmax)
    CALL bcast(undef)
    CALL bcast(tstream)
    CALL bcast(tfast)
    CALL bcast(tslow)
    CALL bcast(tflood)
    CALL bcast(cswamp)
    CALL bcast(lfpcri)
    !!
    !!
  END SUBROUTINE routing_hr_graphinfo
  !
!! ================================================================================================================================
!! SUBROUTINE 	: routing_hr_indexfilegrid
!!
!>\BRIEF Locates all the points of the routing graph file on the model grid. This ensure that no assumption is made on the
!!       orientation of the grid in the file.
!!
!! ================================================================================================================================
  SUBROUTINE routing_hr_indexfilegrid(im, jm, nbl, lon, lat, landindex, indextab, il2il)

    USE grid_var
    
    INTEGER(i_std), INTENT(IN) :: im, jm, nbl
    REAL(r_std), INTENT(IN) :: lon(im,jm), lat(im,jm)
    REAL(r_std), INTENT(IN) :: landindex(im,jm)
    INTEGER(i_std), INTENT(INOUT) :: indextab(im,jm)
    INTEGER(i_std), INTENT(OUT) :: il2il(nbl)

    REAL(r_std) :: dirlon, dirlat
    
    IF ( grid_type == regular_xy ) THEN
       !
       ! The grid read in the routing_graph file could have generated with a different Earth radius to the
       ! one used in ORCHIDEE. If the directions of the projected grid are the same we will use
       ! it instead (We already know they are the same dimension). If that is not the case then there is a problem somewhere.
       !
       dirlon=(lon_proj(im,1) - lon_proj(1,1))/ABS(lon_proj(im,1) - lon_proj(1,1))
       dirlat=(lat_proj(1,jm) - lat_proj(1,1))/ABS(lat_proj(1,jm) - lat_proj(1,1))
       IF ( (lon(im,1) - lon(1,1))/ABS(lon(im,1) - lon(1,1)) == dirlon .AND. &
            & (lat(1,jm) - lat(1,1))/ABS(lat(1,jm) - lat(1,1)) == dirlat ) THEN
          CALL routing_hr_indexfilegrid_internal(im, jm, nbl, lon_proj, lat_proj, landindex, indextab, il2il)
       ELSE
          CALL ipslerr_p(3, 'routing_hr_indexfilegrid', 'Grid in routing_graph.nc file does', &
               "not match the one of ORCHIDEE's projection of a regular XY grid.",'')
       ENDIF
    ELSE
       ! Normal procedure
       CALL routing_hr_indexfilegrid_internal(im, jm, nbl, lon, lat, landindex, indextab, il2il)

    ENDIF

  END SUBROUTINE routing_hr_indexfilegrid
  !!
  !! Internal version of Subroutine
  !!
  SUBROUTINE routing_hr_indexfilegrid_internal(im, jm, nbl, lon, lat, landindex, indextab, il2il)
    INTEGER(i_std), INTENT(IN) :: im, jm, nbl
    REAL(r_std), INTENT(IN) :: lon(im,jm), lat(im,jm)
    REAL(r_std), INTENT(IN) :: landindex(im,jm)
    INTEGER(i_std), INTENT(INOUT) :: indextab(im,jm)
    INTEGER(i_std), INTENT(OUT) :: il2il(nbl)
    !
    INTEGER(i_std) :: il,i,j
    INTEGER(i_std) :: f(2)
    INTEGER(i_std) :: ih, jh, ir, jr
    REAL(r_std) :: nd
    REAL(r_std), DIMENSION(im,jm) :: dist
    REAL(r_std) :: mindist = 1000. !! Minimum distance in m between two points to be matched.
    !
    indextab(:,:) = -1
    ih = INT(im/2.)
    ir = NINT(im/2.)+1
    jh = INT(jm/2.)
    jr = NINT(jm/2.)+1
    !
    DO il=1,nbl
       dist(:,:) = undef_sechiba
       DO i=MAX(1,ih-ir),MIN(im,ih+ir)
          DO j=MAX(1,jh-jr),MIN(jm,jh+jr)
             dist(i,j) = haversine_distance(lon(i,j), lat(i,j), lalo_g(il,2), lalo_g(il,1))
          ENDDO
       ENDDO
       f=MINLOC(dist)
       IF ( dist(f(1),f(2)) < mindist ) THEN
          indextab(f(1),f(2)) = il
          il2il(NINT(landindex(f(1),f(2)))) = il
          dist(f(1),f(2)) = undef_sechiba
       ELSE
          CALL ipslerr(3,'routing_hr_indexfilegrid',&
               'Distance of the closest point in the two grids is too large. ', &
               'Are you sure the routing graph file is on the correct grid ?',&
               '  ')
       ENDIF
       !
       ! See if the next point is close by
       !
       nd = haversine_distance(lalo_g(il,2), lalo_g(il,1), lalo_g(MIN(il+1,nbl),2), lalo_g(MIN(il+1,nbl),1))
       IF ( nd < MINVAL(dist)*3 ) THEN
          ! The next point is close so zoom in
          ih = f(1)
          ir = 4
          jh = f(2)
          jr = 4
       ELSE
          ! Back to starting conditions as the next point is far away
          ih = INT(im/2.)
          ir = NINT(im/2.)+1
          jh = INT(jm/2.)
          jr = NINT(jm/2.)+1
       ENDIF
    ENDDO
  END SUBROUTINE routing_hr_indexfilegrid_internal
  !
!! ================================================================================================================================
!! SUBROUTINE 	: routing_hr_convertlandpts
!!
!>\BRIEF In case the order of land points was different in RoutingPreProc and the model. The route_togrid is corrected.
!!
!! ================================================================================================================================
  !
  SUBROUTINE routing_hr_convertlandpts(nbl, nbas, land2land, route_togrid)
    INTEGER(i_std), INTENT(IN)                          :: nbl, nbas
    INTEGER(i_std), INTENT(IN), DIMENSION(nbl)          :: land2land
    INTEGER(i_std), INTENT(INOUT), DIMENSION(nbl,nbas)  :: route_togrid
    !
    INTEGER(i_std) :: ip, ib
    !
    DO ip=1,nbl
       DO ib=1,nbas
          IF ( route_togrid(ip,ib) < undef_int .AND. route_togrid(ip,ib) > 0 ) THEN
             route_togrid(ip,ib) = land2land(route_togrid(ip,ib))
          ELSE
             route_togrid(ip,ib) = ip
          ENDIF
       ENDDO
    ENDDO
    !
  END SUBROUTINE routing_hr_convertlandpts
  !
  !! ================================================================================================================================
  !! SUBROUTINE 	: routing_hr_inflows
  !!
  !>\BRIEF Calculate the inflows from the outflows information.
  !!
  !! ================================================================================================================================
    !
    SUBROUTINE routing_hr_inflows(nbl, nbas, inf, floodplains_glo,route_innum_glo,route_ingrid_glo,route_inbasin_glo)
  
      IMPLICIT None

     INTEGER(i_std), INTENT(IN)                          :: nbl, nbas, inf
     REAL(r_std), INTENT(IN), DIMENSION(nbl,nbas)  :: floodplains_glo
     INTEGER(i_std), INTENT(INOUT), DIMENSION(nbl,nbas)  :: route_innum_glo
     INTEGER(i_std), INTENT(INOUT), DIMENSION(nbl,nbas, inf)  :: route_ingrid_glo, route_inbasin_glo
     !
     INTEGER(i_std) :: ig, ib, og, ob
     !
     route_innum_glo(:,:) = 0
     route_ingrid_glo(:,:,:) = 0
     route_inbasin_glo(:,:,:) = 0
     DO ig=1,nbl
        DO ib=1,nbas
           IF (floodplains_glo(ig,ib) .GT. 0) THEN
              og = route_togrid_glo(ig,ib)
              ob = route_tobasin_glo(ig,ib)
              IF (ob .LE. nbasmax) THEN
                 IF  (floodplains_glo(og,ob) .GT. 0) THEN
                   route_innum_glo(og, ob) = route_innum_glo(og, ob) + 1
                   route_ingrid_glo(og,ob,route_innum_glo(og, ob)) = ig
                   route_inbasin_glo(og,ob,route_innum_glo(og, ob)) = ib
                 END IF 
              END IF
           END IF
        ENDDO
     ENDDO
   END SUBROUTINE routing_hr_inflows
   !
  !!
!! ================================================================================================================================
!! SUBROUTINE 	: routing_hr_landgather
!!
!>\BRIEF Gathers the routing information onto landpoints, i.e. goes from an X/Y grid to a list of land points.
!!
!! ================================================================================================================================  
!
  SUBROUTINE routing_hr_landgather_r(im,jm,nbas,nbl,indextab,ijfield,landfield,def)
    !
    INTEGER(i_std), INTENT(IN) :: im,jm,nbas,nbl
    INTEGER(i_std), INTENT(IN) :: indextab(im,jm)
    REAL(r_std), INTENT(IN), DIMENSION(im,jm,nbas) :: ijfield
    REAL(r_std), INTENT(OUT), DIMENSION(nbl,nbas) :: landfield
    REAL(r_std), INTENT(IN) :: def
    !
    INTEGER(i_std) :: i,j,k
    !
    DO i=1,im
       DO j=1,jm
          IF ( indextab(i,j) > 0 ) THEN
             DO k=1,nbas
                ! Catch undef or NaN values
                IF (ijfield(i,j,k) >= undef_graphfile .OR. ijfield(i,j,k) /= ijfield(i,j,k)) THEN
                   landfield(indextab(i,j),k) = def
                ELSE
                   landfield(indextab(i,j),k) = ijfield(i,j,k)
                ENDIF
             ENDDO
          ENDIF
       ENDDO
    ENDDO
  END SUBROUTINE routing_hr_landgather_r
  SUBROUTINE routing_hr_landgather_i2(im,jm,nbas,nbl,indextab,ijfield,landfield,def)
    !
    INTEGER(i_std), INTENT(IN) :: im,jm,nbas,nbl
    INTEGER(i_std), INTENT(IN) :: indextab(im,jm)
    REAL(r_std), INTENT(IN), DIMENSION(im,jm,nbas) :: ijfield
    INTEGER(i_std), INTENT(OUT), DIMENSION(nbl,nbas) :: landfield
    INTEGER(i_std), INTENT(IN) :: def
    !
    INTEGER(i_std) :: i,j,in
    !
    DO i=1,im
       DO j=1,jm
          IF ( indextab(i,j) > 0 ) THEN
             DO in=1,nbas
                IF (ijfield(i,j,in) .GE. undef_int .OR. ijfield(i,j,in) /= ijfield(i,j,in)) THEN
                   landfield(indextab(i,j),in) = def
                ELSE
                   landfield(indextab(i,j),in) = NINT(ijfield(i,j,in))
                ENDIF
             ENDDO
          ENDIF
       ENDDO
    ENDDO
  END SUBROUTINE routing_hr_landgather_i2
  SUBROUTINE routing_hr_landgather_i1(im,jm,nbl,indextab,ijfield,landfield,def)
    !
    INTEGER(i_std), INTENT(IN) :: im,jm,nbl
    INTEGER(i_std), INTENT(IN) :: indextab(im,jm)
    REAL(r_std), INTENT(IN), DIMENSION(im,jm) :: ijfield
    INTEGER(i_std), INTENT(OUT), DIMENSION(nbl) :: landfield
    INTEGER(i_std), INTENT(IN) :: def
    !
    INTEGER(i_std) :: i,j
    !
    DO i=1,im
       DO j=1,jm
          IF ( indextab(i,j) > 0 ) THEN
             IF (ijfield(i,j) .GE. undef_int .OR. ijfield(i,j) /= ijfield(i,j)) THEN
                landfield(indextab(i,j)) = def
             ELSE
                landfield(indextab(i,j)) = NINT(ijfield(i,j))
             ENDIF
          ENDIF
       ENDDO
    ENDDO
  END SUBROUTINE routing_hr_landgather_i1
  !
  !
!! ================================================================================================================================
!! SUBROUTINE 	: routing_hr_irrigmap
!!
!>\BRIEF         This  subroutine interpolates the 0.5x0.5 degree based map of irrigated areas to the resolution of the model.
!!
!! DESCRIPTION (definitions, functional, design, flags) : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S):
!!
!! REFERENCES   : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

SUBROUTINE routing_hr_irrigmap (nbpt, index, lalo, neighbours, resolution, contfrac, &
       &                       init_irrig, irrigated, init_flood, init_swamp, swamp, hist_id, hist2_id)
    !
    IMPLICIT NONE
    !
!! PARAMETERS
    INTEGER(i_std), PARAMETER                      :: ilake = 1             !! Number of type of lakes area (unitless)
    INTEGER(i_std), PARAMETER                      :: idam = 2              !! Number of type of dams area (unitless)
    INTEGER(i_std), PARAMETER                      :: iflood = 3            !! Number of type of floodplains area (unitless)
    INTEGER(i_std), PARAMETER                      :: iswamp = 4            !! Number of type of swamps area (unitless)
    INTEGER(i_std), PARAMETER                      :: isal = 5              !! Number of type of salines area (unitless)
    INTEGER(i_std), PARAMETER                      :: ipond = 6             !! Number of type of ponds area (unitless)
    INTEGER(i_std), PARAMETER                      :: ntype = 6             !! Number of types of flooded surfaces (unitless)

!! INPUT VARIABLES
    INTEGER(i_std), INTENT(in)                     :: nbpt                  !! Domain size  (unitless)
    INTEGER(i_std), INTENT(in)                     :: index(nbpt)           !! Index on the global map.
    REAL(r_std), INTENT(in)                        :: lalo(nbpt,2)          !! Vector of latitude and longitudes (beware of the order !)
    INTEGER(i_std), INTENT(in)                     :: neighbours(nbpt,NbNeighb)!! Vector of neighbours for each grid point
    REAL(r_std), INTENT(in)                        :: resolution(nbpt,2)    !! The size of each grid box in X and Y (m)
    REAL(r_std), INTENT(in)                        :: contfrac(nbpt)        !! Fraction of land in each grid box (unitless;0-1)
    INTEGER(i_std), INTENT(in)                     :: hist_id               !! Access to history file (unitless)
    INTEGER(i_std), INTENT(in)                     :: hist2_id              !! Access to history file 2 (unitless)
    LOGICAL, INTENT(in)                            :: init_irrig            !! Logical to initialize the irrigation (true/false)
    LOGICAL, INTENT(in)                            :: init_flood            !! Logical to initialize the floodplains (true/false)
    LOGICAL, INTENT(in)                            :: init_swamp            !! Logical to initialize the swamps (true/false)
    !
!! OUTPUT VARIABLES
    REAL(r_std), INTENT(out)                       :: irrigated(:)          !! Irrigated surface in each grid box (m^2)
!!    REAL(r_std), INTENT(out)                       :: floodplains(:)        !! Surface which can be inundated in each grid box (m^2)
    REAL(r_std), INTENT(out)                       :: swamp(:)              !! Surface which can be swamp in each grid box (m^2)
    !
!! LOCAL VARIABLES
    ! Interpolation variables
    ! 
    INTEGER(i_std)                                 :: nbpmax, nix, njx, fopt !!
    CHARACTER(LEN=30)                              :: callsign              !!
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:)     :: resol_lu              !! Resolution read on the map
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:)    :: mask                  !! Mask to exclude some points (unitless)
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)       :: irrsub_area           !! Area on the fine grid (m^2)
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:,:)  :: irrsub_index          !! Indices of the points we need on the fine grid (unitless)
    INTEGER                                        :: ALLOC_ERR             !!
    LOGICAL                                        :: ok_interpol = .FALSE. !! Flag for interpolation (true/false)
    !
    CHARACTER(LEN=80)                              :: filename              !! Name of the netcdf file (unitless)
    INTEGER(i_std)                                 :: iml, jml, lml, tml, fid, ib, ip, jp, itype !! Indices (unitless)
    REAL(r_std)                                    :: lev(1), date, dt, coslat !!
    INTEGER(i_std)                                 :: itau(1)               !!
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)       :: latrel                !! Latitude
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)       :: lonrel                !! Longitude
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)       :: irrigated_frac        !! Irrigated fraction of the grid box (unitless;0-1)
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:)     :: flood_fracmax         !! Maximal flooded fraction of the grid box (unitless;0-1)
    REAL(r_std)                                    :: area_irrig            !! Irrigated surface in the grid box (m^2)
    REAL(r_std)                                    :: area_flood(ntype)     !! Flooded surface in the grid box (m^2)
!!$    REAL(r_std)                                :: irrigmap(nbpt)
!!$    REAL(r_std)                                :: swampmap(nbpt)

!_ ================================================================================================================================

    !
    !Config Key   = IRRIGATION_FILE
    !Config Desc  = Name of file which contains the map of irrigated areas
    !Config Def   = floodplains.nc
    !Config If    = DO_IRRIGATION OR DO_FLOODPLAINS
    !Config Help  = The name of the file to be opened to read the field
    !Config         with the area in m^2 of the area irrigated within each
    !Config         0.5 0.5 deg grid box. The map currently used is the one
    !Config         developed by the Center for Environmental Systems Research 
    !Config         in Kassel (1995).
    !Config Units = [FILE]
    !
    filename = 'floodplains.nc'
    CALL getin_p('IRRIGATION_FILE',filename)
    !
    IF (is_root_prc) THEN
       CALL flininfo(filename,iml, jml, lml, tml, fid)
       CALL flinclo(fid)
    ELSE
       iml = 0
       jml = 0
       lml = 0
       tml = 0
    ENDIF
    !
    CALL bcast(iml)
    CALL bcast(jml)
    CALL bcast(lml)
    CALL bcast(tml)
    !
    !
    !
    ALLOCATE (latrel(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'routing_hr_irrigmap','Pb in allocate for latrel','','')

    ALLOCATE (lonrel(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'routing_hr_irrigmap','Pb in allocate for lonrel','','')

    ALLOCATE (irrigated_frac(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'routing_hr_irrigmap','Pb in allocate for irrigated_frac','','')

    ALLOCATE (flood_fracmax(iml,jml,ntype), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'routing_hr_irrigmap','Pb in allocate for flood_fracmax','','')

    IF (is_root_prc) CALL flinopen(filename, .FALSE., iml, jml, lml, lonrel, latrel, lev, tml, itau, date, dt, fid)

    CALL bcast(lonrel)
    CALL bcast(latrel)
    !
    IF (is_root_prc) CALL flinget(fid, 'irrig', iml, jml, lml, tml, 0, 0, irrigated_frac)
    CALL bcast(irrigated_frac)
    IF (is_root_prc) CALL flinget(fid, 'lake', iml, jml, lml, tml, 0, 0, flood_fracmax(:,:,ilake))
    IF (is_root_prc) CALL flinget(fid, 'dam', iml, jml, lml, tml, 0, 0, flood_fracmax(:,:,idam))
    IF (is_root_prc) CALL flinget(fid, 'flood', iml, jml, lml, tml, 0, 0, flood_fracmax(:,:,iflood))
    IF (is_root_prc) CALL flinget(fid, 'swamp', iml, jml, lml, tml, 0, 0, flood_fracmax(:,:,iswamp))
    IF (is_root_prc) CALL flinget(fid, 'saline', iml, jml, lml, tml, 0, 0, flood_fracmax(:,:,isal))
    IF (is_root_prc) CALL flinget(fid, 'pond', iml, jml, lml, tml, 0, 0, flood_fracmax(:,:,ipond))
    CALL bcast(flood_fracmax)
    !
    IF (is_root_prc) CALL flinclo(fid)
    !
    ! Set to zero all fraction which are less than 0.5%
    !
    DO ip=1,iml
       DO jp=1,jml
          !
          IF ( irrigated_frac(ip,jp) .LT. undef_sechiba-un) THEN
             irrigated_frac(ip,jp) = irrigated_frac(ip,jp)/100.
             IF ( irrigated_frac(ip,jp) < 0.005 ) irrigated_frac(ip,jp) = zero
          ENDIF
          !
          DO itype=1,ntype
             IF ( flood_fracmax(ip,jp,itype) .LT. undef_sechiba-1.) THEN
                flood_fracmax(ip,jp,itype) = flood_fracmax(ip,jp,itype)/100
                IF ( flood_fracmax(ip,jp,itype) < 0.005 )  flood_fracmax(ip,jp,itype) = zero
             ENDIF
          ENDDO
          !
       ENDDO
    ENDDO
    
    IF (printlev>=2) THEN
       WRITE(numout,*) 'lonrel : ', MAXVAL(lonrel), MINVAL(lonrel)
       WRITE(numout,*) 'latrel : ', MAXVAL(latrel), MINVAL(latrel)
       WRITE(numout,*) 'irrigated_frac : ', MINVAL(irrigated_frac, MASK=irrigated_frac .GT. 0), &
            MAXVAL(irrigated_frac, MASK=irrigated_frac .LT. undef_sechiba)
       WRITE(numout,*) 'flood_fracmax : ', MINVAL(flood_fracmax, MASK=flood_fracmax .GT. 0), &
            MAXVAL(flood_fracmax, MASK=flood_fracmax .LT. undef_sechiba)
    END IF

    ! Consider all points a priori
    !
    ALLOCATE(resol_lu(iml,jml,2), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'routing_hr_irrigmap','Pb in allocate for resol_lu','','')

    ALLOCATE(mask(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'routing_hr_irrigmap','Pb in allocate for mask','','')
    mask(:,:) = 0

    DO ip=1,iml
       DO jp=1,jml
          !
          ! Exclude the points where we are close to the missing value.
          !
!MG This condition cannot be applied in floodplains/swamps configuration because
!   the same mask would be used for the interpolation of irrigation, floodplains and swamps maps.
!          IF ( irrigated_frac(ip,jp) < undef_sechiba ) THEN
             mask(ip,jp) = 1
!          ENDIF
          !
          ! Resolution in longitude
          !
          coslat = MAX( COS( latrel(ip,jp) * pi/180. ), mincos )     
          IF ( ip .EQ. 1 ) THEN
             resol_lu(ip,jp,1) = ABS( lonrel(ip+1,jp) - lonrel(ip,jp) ) * pi/180. * R_Earth * coslat
          ELSEIF ( ip .EQ. iml ) THEN
             resol_lu(ip,jp,1) = ABS( lonrel(ip,jp) - lonrel(ip-1,jp) ) * pi/180. * R_Earth * coslat
          ELSE
             resol_lu(ip,jp,1) = ABS( lonrel(ip+1,jp) - lonrel(ip-1,jp) )/2. * pi/180. * R_Earth * coslat
          ENDIF
          !
          ! Resolution in latitude
          !
          IF ( jp .EQ. 1 ) THEN
             resol_lu(ip,jp,2) = ABS( latrel(ip,jp) - latrel(ip,jp+1) ) * pi/180. * R_Earth
          ELSEIF ( jp .EQ. jml ) THEN
             resol_lu(ip,jp,2) = ABS( latrel(ip,jp-1) - latrel(ip,jp) ) * pi/180. * R_Earth
          ELSE
             resol_lu(ip,jp,2) =  ABS( latrel(ip,jp-1) - latrel(ip,jp+1) )/2. * pi/180. * R_Earth
          ENDIF
          !
       ENDDO
    ENDDO
    !
    ! The number of maximum vegetation map points in the GCM grid is estimated.
    ! Some lmargin is taken.
    !
    callsign = 'Irrigation map'
    ok_interpol = .FALSE.
    IF (is_root_prc) THEN
       nix=INT(MAXVAL(resolution_g(:,1))/MAXVAL(resol_lu(:,:,1)))+2
       njx=INT(MAXVAL(resolution_g(:,2))/MAXVAL(resol_lu(:,:,2)))+2
       nbpmax = nix*njx*2
       IF (printlev>=1) THEN
          WRITE(numout,*) "Projection arrays for ",callsign," : "
          WRITE(numout,*) "nbpmax = ",nbpmax, nix, njx
       END IF
    ENDIF
    CALL bcast(nbpmax)

    ALLOCATE(irrsub_index(nbpt, nbpmax, 2), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'routing_hr_irrigmap','Pb in allocate for irrsub_index','','')
    irrsub_index(:,:,:)=0

    ALLOCATE(irrsub_area(nbpt, nbpmax), STAT=ALLOC_ERR)
    IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'routing_hr_irrigmap','Pb in allocate for irrsub_area','','')
    irrsub_area(:,:)=zero

    CALL aggregate_p(nbpt, lalo, neighbours, resolution, contfrac, &
         &                iml, jml, lonrel, latrel, mask, callsign, &
         &                nbpmax, irrsub_index, irrsub_area, ok_interpol)
    !
    !
    WHERE (irrsub_area < 0) irrsub_area=zero
    !  
    ! Test here if not all sub_area are larger than 0 if so, then we need to increase nbpmax
    !
    DO ib=1,nbpt
       !
       area_irrig = 0.0
       area_flood = 0.0
       !
       DO fopt=1,COUNT(irrsub_area(ib,:) > zero)
          !
          ip = irrsub_index(ib, fopt, 1)
          jp = irrsub_index(ib, fopt, 2)
          !
          IF (irrigated_frac(ip,jp) .LT. undef_sechiba-1.) THEN
             area_irrig = area_irrig + irrsub_area(ib,fopt)*irrigated_frac(ip,jp)
          ENDIF
          !
          DO itype=1,ntype
             IF (flood_fracmax(ip,jp,itype) .LT. undef_sechiba-1.) THEN
                area_flood(itype) = area_flood(itype) + irrsub_area(ib,fopt)*flood_fracmax(ip,jp,itype)
             ENDIF
          ENDDO
       ENDDO
       !
       ! Put the total irrigated and flooded areas in the output variables
       !
       IF ( init_irrig ) THEN
          irrigated(ib) = MIN(area_irrig, resolution(ib,1)*resolution(ib,2)*contfrac(ib))
          IF ( irrigated(ib) < 0 ) THEN
             WRITE(numout,*) 'We have a problem here : ', irrigated(ib) 
             WRITE(numout,*) 'resolution :', resolution(ib,1), resolution(ib,2)
             WRITE(numout,*) area_irrig
             CALL ipslerr_p(3,'routing_hr_irrigmap','Problem with irrigated...','','')
          ENDIF
!!$          ! Compute a diagnostic of the map.
!!$          IF(contfrac(ib).GT.zero) THEN
!!$             irrigmap (ib) = irrigated(ib) / ( resolution(ib,1)*resolution(ib,2)*contfrac(ib) )
!!$          ELSE
!!$             irrigmap (ib) = zero
!!$          ENDIF
          !
       ENDIF
       !
       !
       !
       IF ( init_swamp ) THEN
          swamp(ib) = MIN(area_flood(iswamp), resolution(ib,1)*resolution(ib,2)*contfrac(ib))
          IF ( swamp(ib) < 0 ) THEN
             WRITE(numout,*) 'We have a problem here : ', swamp(ib) 
             WRITE(numout,*) 'resolution :', resolution(ib,1), resolution(ib,2)
             WRITE(numout,*) area_flood
             CALL ipslerr_p(3,'routing_hr_irrigmap','Problem with swamp...','','')
          ENDIF
!!$          ! Compute a diagnostic of the map.
!!$          IF(contfrac(ib).GT.zero) THEN
!!$             swampmap(ib) = swamp(ib) / ( resolution(ib,1)*resolution(ib,2)*contfrac(ib) )
!!$          ELSE
!!$             swampmap(ib) = zero
!!$          ENDIF
       ENDIF
       !
       !
    ENDDO
    !
    !
    
    IF (printlev>=1) THEN
       IF ( init_irrig ) WRITE(numout,*) "Diagnostics irrigated :", MINVAL(irrigated), MAXVAL(irrigated)
       IF ( init_flood ) WRITE(numout,*) "Diagnostics floodplains :", MINVAL(floodplains), MAXVAL(floodplains)
       IF ( init_swamp ) WRITE(numout,*) "Diagnostics swamp :", MINVAL(swamp), MAXVAL(swamp)
    END IF

! No compensation is done for overlapping floodplains, swamp and irrig. At least overlapping will not
! happen between floodplains and swamp alone
!    IF ( init_irrig .AND. init_flood ) THEN
!       DO ib = 1, nbpt
!          surp = (floodplains(ib)+swamp(ib)+irrigated(ib)) / (resolution(ib,1)*resolution(ib,2)*contfrac(ib))
!          IF ( surp .GT. un ) THEN
!             floodplains(ib) = floodplains(ib) / surp
!             swamp(ib) = swamp(ib) / surp
!             irrigated(ib) = irrigated(ib) / surp
!          ENDIF
!       ENDDO
!    ENDIF
    !
    DEALLOCATE (irrsub_area)
    DEALLOCATE (irrsub_index)
    !
    DEALLOCATE (mask)
    DEALLOCATE (resol_lu)
    !
    DEALLOCATE (lonrel)
    DEALLOCATE (latrel)
    !
  END SUBROUTINE routing_hr_irrigmap
  !
!! ================================================================================================================================
!! SUBROUTINE 	: routing_hr_waterbal
!!
!>\BRIEF         This subroutine checks the water balance in the routing module.
!!
!! DESCRIPTION (definitions, functional, design, flags) : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S):
!!
!! REFERENCES   : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

SUBROUTINE routing_hr_waterbal(nbpt, reinit, floodout, runoff, drainage, returnflow, &
               & reinfiltration, irrigation, riverflow, coastalflow)
    !
    IMPLICIT NONE
    !
!! INPUT VARIABLES
    INTEGER(i_std), INTENT(in) :: nbpt                 !! Domain size  (unitless)
    LOGICAL, INTENT(in)        :: reinit               !! Controls behaviour (true/false)
    REAL(r_std), INTENT(in)    :: floodout(nbpt)       !! Grid-point flow out of floodplains (kg/m^2/dt)
    REAL(r_std), INTENT(in)    :: runoff(nbpt)         !! Grid-point runoff (kg/m^2/dt)
    REAL(r_std), INTENT(in)    :: drainage(nbpt)       !! Grid-point drainage (kg/m^2/dt)
    REAL(r_std), INTENT(in)    :: returnflow(nbpt)     !! The water flow from lakes and swamps which returns to the grid box.
                                                       !! This water will go back into the hydrol module to allow re-evaporation (kg/m^2/dt)
    REAL(r_std), INTENT(in)    :: reinfiltration(nbpt) !! Water flow from ponds and floodplains which returns to the grid box (kg/m^2/dt)
    REAL(r_std), INTENT(in)    :: irrigation(nbpt)     !! Irrigation flux. This is the water taken from the reservoirs and beeing put into the upper layers of the soil (kg/m^2/dt)
    REAL(r_std), INTENT(in)    :: riverflow(nbpt)      !! Outflow of the major rivers. The flux will be located on the continental grid but this should be a coastal point (kg/dt)
    REAL(r_std), INTENT(in)    :: coastalflow(nbpt)    !! Outflow on coastal points by small basins. This is the water which flows in a disperse way into the ocean (kg/dt)
    !
    ! We sum-up all the water we have in the warious reservoirs
    !
    REAL(r_std), SAVE          :: totw_flood           !! Sum of all the water amount in the floodplains reservoirs (kg)
!$OMP THREADPRIVATE(totw_flood)
    REAL(r_std), SAVE          :: totw_stream          !! Sum of all the water amount in the stream reservoirs (kg)
!$OMP THREADPRIVATE(totw_stream)
    REAL(r_std), SAVE          :: totw_fast            !! Sum of all the water amount in the fast reservoirs (kg)
!$OMP THREADPRIVATE(totw_fast)
    REAL(r_std), SAVE          :: totw_slow            !! Sum of all the water amount in the slow reservoirs (kg)
!$OMP THREADPRIVATE(totw_slow)
    REAL(r_std), SAVE          :: totw_lake            !! Sum of all the water amount in the lake reservoirs (kg)
!$OMP THREADPRIVATE(totw_lake)
    REAL(r_std), SAVE          :: totw_pond            !! Sum of all the water amount in the pond reservoirs (kg)
!$OMP THREADPRIVATE(totw_pond)
    REAL(r_std), SAVE          :: totw_in              !! Sum of the water flow in to the routing scheme
!$OMP THREADPRIVATE(totw_in)
    REAL(r_std), SAVE          :: totw_out             !! Sum of the water flow out to the routing scheme
!$OMP THREADPRIVATE(totw_out)
    REAL(r_std), SAVE          :: totw_return          !! 
!$OMP THREADPRIVATE(totw_return)
    REAL(r_std), SAVE          :: totw_irrig           !! 
!$OMP THREADPRIVATE(totw_irrig)
    REAL(r_std), SAVE          :: totw_river           !! 
!$OMP THREADPRIVATE(totw_river)
    REAL(r_std), SAVE          :: totw_coastal         !! 
!$OMP THREADPRIVATE(totw_coastal)
    REAL(r_std)                :: totarea              !! Total area of basin (m^2)
    REAL(r_std)                :: area                 !! Total area of routing (m^2)
    INTEGER(i_std)             :: ig                   !! 
    !
    ! Just to make sure we do not get too large numbers !
    !
!! PARAMETERS
    REAL(r_std), PARAMETER     :: scaling = 1.0E+6     !!
    REAL(r_std), PARAMETER     :: allowed_err = 50.    !!

!_ ================================================================================================================================
    !
    IF ( reinit ) THEN
       !
       totw_flood = zero
       totw_stream = zero
       totw_fast = zero
       totw_slow = zero
       totw_lake = zero
       totw_pond = zero 
       totw_in = zero
       !
       DO ig=1,nbpt
          !
          totarea = SUM(routing_area(ig,:))
          !
          totw_flood = totw_flood + SUM(flood_reservoir(ig,:)/scaling)
          totw_stream = totw_stream + SUM(stream_reservoir(ig,:)/scaling)
          totw_fast = totw_fast + SUM(fast_reservoir(ig,:)/scaling)
          totw_slow = totw_slow + SUM(slow_reservoir(ig,:)/scaling)
          totw_lake = totw_lake + lake_reservoir(ig)/scaling
          totw_pond = totw_pond + pond_reservoir(ig)/scaling
          !
          totw_in = totw_in + (runoff(ig)*totarea + drainage(ig)*totarea - floodout(ig)*totarea)/scaling
          !
       ENDDO
       !
    ELSE
       !
       totw_out = zero
       totw_return = zero
       totw_irrig = zero
       totw_river = zero
       totw_coastal = zero
       area = zero
       !
       DO ig=1,nbpt
          !
          totarea = SUM(routing_area(ig,:))
          !
          totw_flood = totw_flood - SUM(flood_reservoir(ig,:)/scaling)
          totw_stream = totw_stream - SUM(stream_reservoir(ig,:)/scaling)
          totw_fast = totw_fast - SUM(fast_reservoir(ig,:)/scaling)
          totw_slow = totw_slow - SUM(slow_reservoir(ig,:)/scaling)
          totw_lake = totw_lake - lake_reservoir(ig)/scaling
          totw_pond = totw_pond - pond_reservoir(ig)/scaling
          !
          totw_return = totw_return + (reinfiltration(ig)+returnflow(ig))*totarea/scaling
          totw_irrig = totw_irrig + irrigation(ig)*totarea/scaling
          totw_river = totw_river + riverflow(ig)/scaling
          totw_coastal = totw_coastal + coastalflow(ig)/scaling
          !
          area = area + totarea
          !
       ENDDO
       totw_out = totw_return + totw_irrig + totw_river + totw_coastal
       !
       ! Now we have all the information to balance our water
       !
       IF ( ABS((totw_flood + totw_stream + totw_fast + totw_slow + totw_lake + totw_pond) - &
            & (totw_out - totw_in)) > allowed_err ) THEN
          WRITE(numout,*) 'WARNING : Water not conserved in routing. Limit at ', allowed_err, ' 10^6 kg'
          WRITE(numout,*) '--Water-- change : flood stream fast ', totw_flood, totw_stream, totw_fast
          WRITE(numout,*) '--Water-- change : slow, lake ', totw_slow, totw_lake
          WRITE(numout,*) '--Water>>> change in the routing res. : ', totw_flood + totw_stream + totw_fast + totw_slow + totw_lake
          WRITE(numout,*) '--Water input : ', totw_in
          WRITE(numout,*) '--Water output : ', totw_out
          WRITE(numout,*) '--Water output : return, irrig ', totw_return, totw_irrig
          WRITE(numout,*) '--Water output : river, coastal ',totw_river, totw_coastal
          WRITE(numout,*) '--Water>>> change by fluxes : ', totw_out - totw_in, ' Diff [mm/dt]: ',   &
               & ((totw_flood + totw_stream + totw_fast + totw_slow + totw_lake) - (totw_out - totw_in))/area

          ! Stop the model
          CALL ipslerr_p(3, 'routing_hr_waterbal', 'Water is not conserved in routing.','','')
       ENDIF
       !
    ENDIF
    !
  END SUBROUTINE routing_hr_waterbal
  !
  !
END MODULE routing_subgrid_htu
