! =================================================================================================================================
! MODULE       : routing_interp_topo_mod
!
!  CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       This module routes the water over the continents into the oceans and computes the water
!!             stored in floodplains or taken for irrigation.
!!
!!\n DESCRIPTION: The subroutines in this subroutine is only called when ROUTING_METHOD=interp_topo is set in run.def.
!!                The method can be used for regular latitude-longitude grid or for unstructured grid.
!!                It computes the routing on a separate grid and exchanges data with the ORCHIDEE grid using interpolations
!!                performed via XIOS.
!!                It is a complete and more robust rewriting of routing_interp_topo_v0.f90
!!
!! RECENT CHANGE(S): July 2024: Changed routing method's name "native" to "interp_topo".
!!
!! REFERENCE(S)	:
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_sechiba/routing_interp_topo.f90 $
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


MODULE routing_interp_topo_mod
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
  PUBLIC :: routing_interp_topo_main, routing_interp_topo_xios_initialize
  PUBLIC :: routing_interp_topo_initialize, routing_interp_topo_finalize, routing_interp_topo_clear

  REAL,SAVE :: dt_routing
  !$OMP THREADPRIVATE(dt_routing) 

  INTEGER,SAVE :: nbpt
  !$OMP THREADPRIVATE(nbpt) 
  
  INTEGER,SAVE :: nbpt_r
  !$OMP THREADPRIVATE(nbpt_r) 
  
  LOGICAL,SAVE :: dofloodinfilt
  !$OMP THREADPRIVATE(dofloodinfilt)   

  LOGICAL,SAVE :: doswamps
  !$OMP THREADPRIVATE(doswamps)  

  LOGICAL,SAVE :: doponds
  !$OMP THREADPRIVATE(doponds)  


  REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:)               :: returnflow_mean              !! Mean water flow from lakes and swamps which returns to the grid box.
 !$OMP THREADPRIVATE(returnflow_mean)

 REAL, SAVE :: time_counter = 0
 !$OMP THREADPRIVATE(time_counter)  
   
 LOGICAL,SAVE :: check_water_balance =.FALSE.
 !$OMP THREADPRIVATE(check_water_balance) 

CONTAINS

!!  =============================================================================================================================
!! SUBROUTINE:    routing_interp_topo_xios_initialize
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

  SUBROUTINE routing_interp_topo_xios_initialize
    USE xios
    USE routing_interp_topo_flow_mod
    IMPLICIT NONE     
     
      CALL xios_orchidee_set_field_attr("basinmap",enabled=.FALSE.)
      CALL xios_orchidee_set_field_attr("nbrivers",enabled=.FALSE.)
      CALL xios_orchidee_set_field_attr("riversret",enabled=.FALSE.)
      CALL xios_orchidee_set_field_attr("hydrographs",enabled=.TRUE.)
      CALL xios_orchidee_set_field_attr("htuhgmon",enabled=.FALSE.)
      CALL xios_orchidee_set_field_attr("fastr",enabled=.TRUE.)
      CALL xios_orchidee_set_field_attr("slowr",enabled=.TRUE.)
      CALL xios_orchidee_set_field_attr("streamr",enabled=.TRUE.)
      CALL xios_orchidee_set_field_attr("laker",enabled=.FALSE.)
      CALL xios_orchidee_set_field_attr("lake_overflow",enabled=.FALSE.)
      CALL xios_orchidee_set_field_attr("mask_coast",enabled=.FALSE.)
      CALL xios_orchidee_set_field_attr("pondr",enabled=.FALSE.)
      CALL xios_orchidee_set_field_attr("floodr",enabled=.FALSE.)
      CALL xios_orchidee_set_field_attr("slowflow",enabled=.FALSE.)
      CALL xios_orchidee_set_field_attr("delfastr",enabled=.FALSE.)
      CALL xios_orchidee_set_field_attr("delslowr",enabled=.FALSE.)
      CALL xios_orchidee_set_field_attr("delstreamr",enabled=.FALSE.)
      CALL xios_orchidee_set_field_attr("dellaker",enabled=.FALSE.)
      CALL xios_orchidee_set_field_attr("delpondr",enabled=.FALSE.)
      CALL xios_orchidee_set_field_attr("delfloodr",enabled=.FALSE.)
      CALL xios_orchidee_set_field_attr("irrigmap",enabled=.FALSE.)
      CALL xios_orchidee_set_field_attr("swampmap",enabled=.FALSE.)
      CALL xios_orchidee_set_field_attr("wbr_stream",enabled=.FALSE.)
      CALL xios_orchidee_set_field_attr("wbr_fast",enabled=.FALSE.)
      CALL xios_orchidee_set_field_attr("wbr_slow",enabled=.FALSE.)
      CALL xios_orchidee_set_field_attr("wbr_lake",enabled=.FALSE.)
      CALL xios_orchidee_set_field_attr("reinfiltration",enabled=.FALSE.)
      CALL xios_orchidee_set_field_attr("SurfStor",enabled=.FALSE.)
      CALL xios_orchidee_set_field_attr("htutempmon",enabled=.FALSE.)
      CALL xios_orchidee_set_field_attr("streamlimit",enabled=.FALSE.)
      CALL xios_orchidee_set_field_attr("StreamT_TotTend",enabled=.FALSE.)
      CALL xios_orchidee_set_field_attr("StreamT_AdvTend",enabled=.FALSE.)
      CALL xios_orchidee_set_field_attr("StreamT_RelTend",enabled=.FALSE.)

     CALL routing_flow_xios_initialize

     
   END SUBROUTINE routing_interp_topo_xios_initialize







  !!  =============================================================================================================================
  !! SUBROUTINE:         routing_interp_topo_initialize
  !!
  !>\BRIEF	         Initialize the routing_interp_topo module
  !!
  !! DESCRIPTION:        Initialize the routing_interp_topo module. Read from restart file or read the routing.nc file
  !!                     to initialize the routing scheme. 
  !!
  !! RECENT CHANGE(S)
  !!
  !! REFERENCE(S)
  !! 
  !! FLOWCHART   
  !! \n
  !_ ==============================================================================================================================

  SUBROUTINE routing_interp_topo_initialize( kjit,       nbpt,           index,                 &
                                             rest_id,     hist_id,        hist2_id,   lalo,      &
                                             neighbours,  resolution,     contfrac,   stempdiag, &
                                             returnflow,  reinfiltration, irrigation, riverflow, &
                                             coastalflow, flood_frac,     flood_res , irrigated_next)
    
    USE routing_interp_topo_flow_mod, ONLY : routing_flow_initialize
    USE routing_interp_topo_irrig_mod, ONLY : irrigation_initialize
    USE routing_interp_topo_lake_mod, ONLY : routing_lake_initialize
    USE routing_interp_topo_para, ONLY: routing_para_initialize
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
    REAL(r_std), INTENT(in)        :: irrigated_next (nbpt)  !! Dynamic irrig. area, calculated in slowproc and passed to routing!

    !_ ================================================================================================================================
    
    
      CALL routing_para_initialize 
      CALL routing_interp_topo_init_local(kjit, rest_id, nbpt, contfrac)
      CALL routing_flow_initialize(kjit, rest_id, nbpt, dt_routing, contfrac, nbpt_r, riverflow, coastalflow)
      CALL routing_lake_initialize(kjit, rest_id, nbpt, contfrac)
      CALL irrigation_initialize(kjit, rest_id, nbpt, nbpt_r,irrigated_next )
      
      reinfiltration(:)=0
      irrigation(:)=0
      flood_frac(:)=0
      flood_res(:)=0
      returnflow(:)=returnflow_mean(:)

    END SUBROUTINE routing_interp_topo_initialize
    
    
  !! ================================================================================================================================
  !! SUBROUTINE 	: routing_interp_topo_initialize
  !!
  !>\BRIEF         This subroutine allocates the memory and get the fixed fields from the restart file.
  !!
  !! DESCRIPTION:	  Private subroutine to the module routing_interp_topo. This subroutine is called in the begining 
  !!                      of routing_interp_topo_initialize
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

  SUBROUTINE routing_interp_topo_init_local(kjit, rest_id, nbpt_, contfrac)
  USE time, ONLY : dt_sechiba 
  USE routing_interp_topo_flow_mod   
  IMPLICIT NONE
    INTEGER(i_std),INTENT(in)      :: kjit
    INTEGER(i_std),INTENT(in)      :: rest_id              !! Restart file identifier (unitless)
    INTEGER,INTENT(in) :: nbpt_             !! nb points native grid
    REAL,INTENT(in)    :: contfrac(nbpt)   !! fraction of land

    INTEGER :: ier    
    CHARACTER(LEN=80)   :: var_name       !! To store variables names for I/O (unitless)
  
    REAL(r_std) ::ratio
    
    nbpt = nbpt_

    ALLOCATE (returnflow_mean(nbpt), stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_interp_topo_init_local','Pb in allocate for returnflow_mean','','')
    var_name = 'returnflow'
    CALL ioconf_setatt_p('UNITS', 'Kg/m^2/dt')
    CALL ioconf_setatt_p('LONG_NAME','Deep return flux')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., returnflow_mean, "gather", nbp_glo, index_g)
    CALL setvar_p (returnflow_mean, val_exp, 'NO_KEYWORD', zero)
    
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

    check_water_balance=.FALSE.
    CALL getin_p('ROUTING_CHECK_WATER_BALANCE', check_water_balance)

    ratio = dt_routing/dt_sechiba
    IF ( ABS(NINT(ratio) - ratio) .GT. 10*EPSILON(ratio)) THEN
       WRITE(numout,*) 'WARNING -- WARNING -- WARNING -- WARNING'
       WRITE(numout,*) "The chosen time step for the routing is not a multiple of the"
       WRITE(numout,*) "main time step of the model. We will change dt_routing so that"
       WRITE(numout,*) "this condition is fulfilled"
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
    
  END SUBROUTINE routing_interp_topo_init_local
  


  !! ================================================================================================================================
  !! SUBROUTINE   : routing_interp_topo_main 
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


  SUBROUTINE routing_interp_topo_main(kjit, nbpt, index, &
                                      lalo, neighbours, resolution, contfrac, totfrac_nobio, veget_max, floodout, runoff, &
                                      drainage, transpot, precip_rain, humrel, k_litt, flood_frac, flood_res, &
                                      stempdiag, reinf_slope, returnflow, reinfiltration, irrigation, riverflow, coastalflow, &
                                      rest_id, hist_id, hist2_id ,&
                                      soiltile, root_deficit, irrigated_next, irrig_frac, fraction_aeirrig_sw) 

    USE routing_interp_topo_flow_mod,ONLY : routing_flow_make_mean, routing_flow_main, riverflow_mean, coastalflow_mean, lakeinflow_mean, &
                                            routing_flow_get, routing_flow_set, routing_flow_diags, routing_flow_retreive_outflow, routing_flow_get_water_balance
    USE routing_interp_topo_irrig_mod, ONLY: irrigation_mean_make, irrigation_main, irrigation_get
    USE routing_interp_topo_lake_mod, ONLY: routing_lake_mean_make, routing_lake_route_coast, routing_lake_main, routing_lake_get_water_balance
    USE xios
    USE grid, ONLY : area 
    USE mod_orchidee_para, ONLY : reduce_sum, is_mpi_root, is_omp_root
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
    REAL(r_std), INTENT(in)        :: root_deficit(nbpt)   !! soil water deficit
    REAL(r_std), INTENT(in)        :: soiltile(nbpt,nstm)  !! Fraction of each soil tile within vegtot (0-1, unitless)
    REAL(r_std), INTENT(in)        :: irrig_frac(nbpt)     !! Irrig. fraction interpolated in routing, and saved to pass to slowproc if irrigated_soiltile = .TRUE.
    REAL(r_std), INTENT(in)        :: irrigated_next (nbpt)!! Dynamic irrig. area, calculated in slowproc and passed to routing
    REAL(r_std), INTENT(in)        :: fraction_aeirrig_sw(nbpt) !! Fraction of area equipped for irrigation from surface water, of irrig_frac

    !! 0.2 Output variables
    REAL(r_std), INTENT(out)       :: returnflow(nbpt)     !! The water flow from lakes and swamps which returns to the grid box.
    !! This water will go back into the hydrol or hydrolc module to allow re-evaporation (kg/m^2/dt)
    REAL(r_std), INTENT(out)       :: reinfiltration(nbpt) !! Water flow from ponds and floodplains which returns to the grid box (kg/m^2/dt)
    REAL(r_std), INTENT(out)       :: irrigation(nbpt)     !! Irrigation flux. This is the water taken from the reservoirs to be put into the upper layers of the soil (kg/m^2/dt)
    REAL(r_std), INTENT(out)       :: riverflow(nbpt)      !! Outflow of the major rivers. The flux will be located on the continental grid but this should be a coastal point (kg/dt)
    REAL(r_std), INTENT(out)       :: coastalflow(nbpt)    !! Outflow on coastal points by small basins. This is the water which flows in a disperse way into the ocean (kg/dt)
    REAL(r_std), INTENT(out)       :: flood_frac(nbpt)     !! Flooded fraction of the grid box (unitless;0-1)
    REAL(r_std), INTENT(out)       :: flood_res(nbpt)      !! Diagnostic of water amount in the floodplains reservoir (kg)

    !! 0.3 Local variables
    REAL(r_std), DIMENSION(nbpt)   :: return_lakes         !! Water from lakes flowing back into soil moisture (kg/m^2/dt)
    REAL(r_std), DIMENSION(nbpt)   :: lakeinflow
    REAL(r_std), DIMENSION(nbpt)   :: irrig_netereq        !! diag, Irrigation requirement flux (kg/m^2/dt)
    REAL(r_std)                    :: water_balance_before, before
    REAL(r_std)                    :: water_balance_after, after

    IF (check_water_balance) THEN
      water_balance_before = 0
      water_balance_before = water_balance_before + SUM( (runoff + drainage) *contfrac*area )
      water_balance_before = water_balance_before + routing_flow_get_water_balance(contfrac)
      water_balance_before = water_balance_before + routing_lake_get_water_balance(contfrac)
    ENDIF

    CALL routing_flow_make_mean(runoff, drainage)
    CALL irrigation_mean_make(dt_routing, veget_max, humrel, transpot, runoff , precip_rain)
    CALL routing_lake_mean_make(dt_routing, humrel, veget_max)
       
    time_counter = time_counter + dt_sechiba
    
    reinfiltration(:)= 0 ! for now
    
    ! If the time has come we do the routing.
    IF ( NINT(time_counter) .GE. NINT(dt_routing) ) THEN
       returnflow_mean(:)=0
 
       CALL routing_flow_main(dt_routing, contfrac)

       CALL routing_flow_get(coastalflow_mean=coastalflow)
       CALL routing_lake_route_coast(contfrac,coastalflow)
       CALL routing_flow_set(coastalflow_mean=coastalflow)
      
       CALL irrigation_main( dt_routing, contfrac, reinfiltration, irrigated_next, irrig_frac, root_deficit, soiltile, &
                             fraction_aeirrig_sw )
       
       CALL routing_flow_get(lakeinflow_mean=lakeinflow)
       CALL routing_lake_main(dt_routing, contfrac, lakeinflow, return_lakes)
       lakeinflow(:) = 0. ! lake reservoir has been updated in routing_lake_main, lakeinflow can thus be put to 0 to ensure consistency of the water budget
       CALL routing_flow_set(lakeinflow_mean=lakeinflow)
       
       returnflow_mean(:)=returnflow_mean(:)+return_lakes(:)
       
       time_counter = 0
    ENDIF

    CALL routing_flow_retreive_outflow(dt_sechiba, riverflow, coastalflow, lakeinflow)
    CALL irrigation_get(routing_irrigation_mean=irrigation, routing_irrig_netereq_mean=irrig_netereq)

    
    CALL routing_flow_diags(dt_sechiba) ! output routing flow diagnostics on orchidee grid

    ! Convert fluxes from kg/m^2/dt_routing to kg/m^2/dt(_sechiba)
    returnflow  = returnflow/dt_routing*dt_sechiba
    irrigation = irrigation/dt_routing*dt_sechiba
    irrig_netereq = irrig_netereq/dt_routing*dt_sechiba
    
    CALL xios_orchidee_send_field("irrigation", irrigation/dt_sechiba)
    CALL xios_orchidee_send_field("netirrig", irrig_netereq/dt_sechiba)
    
    returnflow(:) = 0
    reinfiltration(:) =0
    flood_frac(:) = 0
    flood_res(:) = 0
    
    IF (check_water_balance) THEN
      water_balance_after = 0
      water_balance_after = water_balance_after + SUM( riverflow + coastalflow + lakeinflow)
      water_balance_after = water_balance_after + routing_flow_get_water_balance(contfrac)
      water_balance_after = water_balance_after + routing_lake_get_water_balance(contfrac)

      CALL reduce_sum(water_balance_before, before)
      CALL reduce_sum(water_balance_after, after)
      IF (is_mpi_root .AND. is_omp_root) THEN
        WRITE (numout,*) "Routing interp topo :  Water balance_before : ", before
        WRITE (numout,*) "Routing interp topo :  Water balance after : ", after
        WRITE (numout,*) "Routing interp topo :  Water balance relative : ", (after-before)/(0.5*(after+before+1.e-100))
      ENDIF
    ENDIF

    ! convert coastalflow and riverflow from kg/dt_sechiba to m3/dt_sechiba
    coastalflow(:)=coastalflow(:)/1000.
    riverflow(:)=riverflow(:)/1000.
    
  END SUBROUTINE routing_interp_topo_main

  !!  =============================================================================================================================
  !! SUBROUTINE:         routing_interp_topo_finalize
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

  SUBROUTINE routing_interp_topo_finalize(kjit, nbpt_, rest_id, flood_frac, flood_res )
    USE routing_interp_topo_flow_mod
    USE routing_interp_topo_lake_mod
    USE routing_interp_topo_irrig_mod
    IMPLICIT NONE
    INTEGER, INTENT(IN)          :: kjit
    INTEGER,INTENT(IN)           :: nbpt_
    INTEGER, INTENT(IN)          :: rest_id
    REAL(r_std), INTENT(in)      :: flood_frac(nbpt)     !! Flooded fraction of the grid box (unitless;0-1)
    REAL(r_std), INTENT(in)      :: flood_res(nbpt)      !! Diagnostic of water amount in the floodplains reservoir (kg)

    CALL routing_interp_topo_local_finalize(kjit, rest_id)
    CALL routing_flow_finalize(kjit, rest_id)
    CALL routing_lake_finalize(kjit, rest_id)
    CALL irrigation_finalize(kjit, rest_id)

  END SUBROUTINE routing_interp_topo_finalize

  SUBROUTINE routing_interp_topo_local_finalize(kjit,rest_id)
    USE ioipsl_para
    USE grid
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: kjit
    INTEGER, INTENT(IN) :: rest_id

    CALL restput_p (rest_id, 'returnflow', nbp_glo, 1, 1, kjit, returnflow_mean, 'scatter',  nbp_glo, index_g)
    DEALLOCATE(returnflow_mean)

  END SUBROUTINE routing_interp_topo_local_finalize


  SUBROUTINE routing_interp_topo_clear
    USE xios
    IMPLICIT NONE



  END SUBROUTINE routing_interp_topo_clear
#endif  
END MODULE routing_interp_topo_mod
