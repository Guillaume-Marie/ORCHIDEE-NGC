! ================================================================================================================================
!  MODULE       : routing_wrapper
!
!  CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
!  LICENCE      : IPSL (2006)
!  This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF          Interfaces to all routing schemes.
!!
!!\n DESCRIPTION: This module contains uniformed subroutines called from sechiba. These subroutines make the swich the between 
!!                the different existing routing modules.
!!                
!!                Depending on the key world ROUTING_METHOD set in run.def, this module calls one of the 
!!                available routing modules: 
!!                - ROUTING_METHOD=subgrid_halfdeg for the routing scheme available in module routing_subgrid_halfdeg.
!!                - ROUTING_METHOD=subgrid_htu for the high resolution routing scheme in module routing_subgrid_htu. 
!!                - ROUTING_METHOD=interp_topo_v0 for the routing scheme in module routing_interp_topo_v0.
!!                - ROUTING_METHOD=interp_topo for the routing scheme in module routing_interp_topo.
!!
!! RECENT CHANGE(S): July 2024: Changed names of the different routing methods.
!!
!! REFERENCE(S) : None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_sechiba/routing_wrapper.f90 $
!! $Date: 2024-09-04 16:47:32 +0200 (mer. 04 sept. 2024) $
!! $Revision: 8648 $
!! \n
!_ ================================================================================================================================

MODULE routing_wrapper

  USE defprec
  USE pft_parameters
  USE grid, not_used => contfrac
  USE routing_subgrid_halfdeg
  USE routing_subgrid_htu
  USE routing_interp_topo_v0
  USE routing_interp_topo_mod
  USE constantes_var, ONLY : routing_method

  IMPLICIT NONE

  PUBLIC :: routing_wrapper_xios_initialize, routing_wrapper_initialize, &
            routing_wrapper_main, routing_wrapper_finalize, routing_wrapper_clear 
  PRIVATE

CONTAINS

!!  =============================================================================================================================
!! SUBROUTINE:    routing_wrapper_xios_initialize
!!
!>\BRIEF	  First initialization phase of the choosen routing module
!!
!! DESCRIPTION:	  Call the xios initialization subroutine from corresponding routing module.
!!                ROUTING_METHOD is read previously from run.def.
!!                This subroutine is called before the xios context is closed. 
!!                It is called from sechiba_initialize only if 1 is activated.
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S): None
!! 
!! FLOWCHART: None
!! \n
!_ ==============================================================================================================================
  SUBROUTINE routing_wrapper_xios_initialize()

    IF(routing_method=='subgrid_halfdeg') THEN
       CALL routing_subgrid_halfdeg_xios_initialize
       
    ELSEIF(routing_method=='subgrid_htu') THEN
       CALL routing_subgrid_htu_xios_initialize
       
    ELSEIF(routing_method=='interp_topo_v0') THEN  
#ifdef XIOS
       CALL routing_interp_topo_v0_xios_initialize
#else
       CALL ipslerr_p(3,'routing_wrapper_xios_initialize','ROUTING_METHOD interp_topo_v0 needs XIOS',&
                        'You must compile XIOS and then ORCHIDEE with cpp key XIOS','')
#endif

    ELSEIF(routing_method=='interp_topo') THEN
#ifdef XIOS
       CALL routing_interp_topo_xios_initialize
#else
       CALL ipslerr_p(3,'routing_wrapper_xios_initialize','ROUTING_METHOD interp_topo needs XIOS',&
                        'You must compile XIOS and then ORCHIDEE with cpp key XIOS','')
#endif

    ENDIF

  END SUBROUTINE routing_wrapper_xios_initialize




!!  =============================================================================================================================
!! SUBROUTINE:    routing_wrapper_initialize
!!
!>\BRIEF	  Initialize the choosen routing module
!!
!! DESCRIPTION:	  Call the initialization subroutine from corresponding routing module
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S): None
!! 
!! FLOWCHART: None
!! \n
!_ ==============================================================================================================================
  SUBROUTINE routing_wrapper_initialize( &
       kjit,        nbpt,           index,                 &
       rest_id,     hist_id,        hist2_id,   lalo,      &
       neighbours,  resolution,     contfrac,   stempdiag, ftempdiag, &
       soiltile,    irrig_frac,     veget_max,  irrigated_next, &    
       returnflow,  reinfiltration, irrigation, riverflow, &
       coastalflow, flood_frac,     flood_res )
       


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
    REAL(r_std), INTENT(in)        :: ftempdiag(nbpt,ngrnd)!! Diagnostic soil temperature profile over full column
    REAL(r_std), INTENT(in)        :: soiltile(nbpt,nstm)  !! Fraction of each soil tile within vegtot (0-1, unitless)
    REAL(r_std), INTENT(in)        :: veget_max(nbpt,nvm)  !! Maximal fraction of vegetation (unitless;0-1) !
    REAL(r_std), INTENT(in)        :: irrigated_next (nbpt)!! Dynamic irrig. area, calculated in slowproc and passed to routing!
    REAL(r_std), INTENT(in)        :: irrig_frac(nbpt)     !! Irrig. fraction interpolated in routing, and saved to pass to slowproc if irrigated_soiltile = .TRUE.

    !! 0.2 Output variables
    REAL(r_std), INTENT(out)       :: returnflow(nbpt)     !! The water flow from lakes and swamps which returns to the grid box.
                                                           !! This water will go back into the hydrol or hydrolc module to allow re-evaporation (kg/m^2/dt)
    REAL(r_std), INTENT(out)       :: reinfiltration(nbpt) !! Water flow from ponds and floodplains which returns to the grid box (kg/m^2/dt)
    REAL(r_std), INTENT(out)       :: irrigation(nbpt)     !! Irrigation flux. This is the water taken from the reservoirs and beeing put into the upper layers of the soil (kg/m^2/dt)
    REAL(r_std), INTENT(out)       :: riverflow(nbpt)      !! Outflow of the major rivers. The flux will be located on the continental grid but this should be a coastal point (kg/dt)
    REAL(r_std), INTENT(out)       :: coastalflow(nbpt)    !! Outflow on coastal points by small basins. This is the water which flows in a disperse way into the ocean (kg/dt)
    REAL(r_std), INTENT(out)       :: flood_frac(nbpt)     !! Flooded fraction of the grid box (unitless;0-1)
    REAL(r_std), INTENT(out)       :: flood_res(nbpt)      !! Diagnostic of water amount in the floodplains reservoir (kg)


    !_ ================================================================================================================================

    !! 1. Get routing_method from run.def
    !!    ROUTING_METHOD is now read in control_initialize
        
    !! 2. Initialize the choosen routing module
    IF (routing_method == 'subgrid_halfdeg') THEN

       CALL routing_subgrid_halfdeg_initialize( kjit,        nbpt,           index,                 &
                                                rest_id,     hist_id,        hist2_id,   lalo,      &
                                                neighbours,  resolution,     contfrac,   stempdiag, &
                                                returnflow,  reinfiltration, irrigation, riverflow, &
                                                coastalflow, flood_frac,     flood_res,  soiltile,  &
                                                irrig_frac,  veget_max,      irrigated_next)

    ELSE IF (routing_method == 'subgrid_htu') THEN

       CALL routing_subgrid_htu_initialize(     kjit,        nbpt,           index,                 &
                                                rest_id,     hist_id,        hist2_id,   lalo,      &
                                                neighbours,  resolution,     contfrac,   ftempdiag, &
                                                returnflow,  reinfiltration, irrigation, riverflow, &
                                                coastalflow, flood_frac,     flood_res )

    ELSE IF(routing_method== 'interp_topo_v0') THEN 
#ifdef XIOS
       CALL routing_interp_topo_v0_initialize(  kjit,        nbpt,           index,                 &
                                                rest_id,     hist_id,        hist2_id,   lalo,      &
                                                neighbours,  resolution,     contfrac,   stempdiag, &
                                                returnflow,  reinfiltration, irrigation, riverflow, &
                                                coastalflow, flood_frac,     flood_res )
#endif
       riverflow(:) = zero
       coastalflow(:) = zero
       returnflow(:) = zero
       reinfiltration(:) = zero
       irrigation(:) = zero
       flood_frac(:) = zero
       flood_res(:) = zero
   
    ELSE IF(routing_method== 'interp_topo') THEN
#ifdef XIOS
       CALL routing_interp_topo_initialize(     kjit,        nbpt,           index,                 &
                                                rest_id,     hist_id,        hist2_id,   lalo,      &
                                                neighbours,  resolution,     contfrac,   stempdiag, &
                                                returnflow,  reinfiltration, irrigation, riverflow, &
                                                coastalflow, flood_frac,     flood_res , irrigated_next)
#endif
       riverflow(:) = zero
       coastalflow(:) = zero
       returnflow(:) = zero
       reinfiltration(:) = zero
       irrigation(:) = zero
       flood_frac(:) = zero
       flood_res(:) = zero

    ELSE
       ! Bad choice of routing_method. Exit the model now.
       ! Redundant with a check done in control_initialize (when ROUTING_METHOD is read)
       WRITE(numout,*) 'Following routing method is not implemented, ROUTING_METHOD=',routing_method
       CALL ipslerr_p(3,'routing_wrapper_initialize',&
                        'ROUTING_METHOD can only be subgrid_halfdeg, subgrid_htu, interp_topo_v0 or interp_topo',&
                        'Error in run.def','')
    ENDIF

  END SUBROUTINE routing_wrapper_initialize



!!  =============================================================================================================================
!! SUBROUTINE:    routing_wrapper_main
!!
!>\BRIEF	  Call the main subroutine for the choosen routing module
!!
!! DESCRIPTION:	  Call the main subroutine for the choosen routing module according to ROUTING_METHOD
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S): None
!! 
!! FLOWCHART: None
!! \n
!_ ==============================================================================================================================
  SUBROUTINE routing_wrapper_main(kjit, nbpt, index, &
       lalo, neighbours, resolution, contfrac, totfrac_nobio, veget_max, floodout, runoff, &
       drainage, transpot, precip_rain, humrel, k_litt, flood_frac, flood_res, stempdiag, &
       ftempdiag, reinf_slope, returnflow, reinfiltration, irrigation, riverflow, coastalflow, rest_id, hist_id, hist2_id, &
       soiltile, root_deficit, irrigated_next, irrig_frac, fraction_aeirrig_sw) 

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
    REAL(r_std), INTENT(in)        :: stempdiag(nbpt,nslm) !! Diagnostic soil temperature profile
    REAL(r_std), INTENT(in)        :: ftempdiag(nbpt,ngrnd)!! Diagnostic soil temperature profile over full column
    REAL(r_std), INTENT(in)        :: reinf_slope(nbpt)    !! Coefficient which determines the reinfiltration ratio in the grid box due to flat areas (unitless;0-1)
    REAL(r_std), INTENT(in)        :: root_deficit(nbpt)   !! soil water deficit
    REAL(r_std), INTENT(in)        :: soiltile(nbpt,nstm)  !! Fraction of each soil tile within vegtot (0-1, unitless)
    REAL(r_std), INTENT(in)        :: irrig_frac(nbpt)     !! Irrig. fraction interpolated in routing, and saved to pass to slowproc if irrigated_soiltile = .TRUE.
    REAL(r_std), INTENT(in)        :: irrigated_next (nbpt)!! Dynamic irrig. area, calculated in slowproc and passed to routing
    REAL(r_std), INTENT(in)        :: fraction_aeirrig_sw(nbpt) !! Fraction of area equipped for irrigation from surface water, of irrig_frac

    !! 0.2 Output variables
    REAL(r_std), INTENT(out)       :: returnflow(nbpt)     !! The water flow from lakes and swamps which returns to the grid box.
    !! This water will go back into the hydrol module to allow re-evaporation (kg/m^2/dt)
    REAL(r_std), INTENT(out)       :: reinfiltration(nbpt) !! Water flow from ponds and floodplains which returns to the grid box (kg/m^2/dt)
    REAL(r_std), INTENT(out)       :: irrigation(nbpt)     !! Irrigation flux. This is the water taken from the reservoirs and beeing put into the upper layers of the soil (kg/m^2/dt)
    REAL(r_std), INTENT(out)       :: riverflow(nbpt)      !! Outflow of the major rivers. The flux will be located on the continental grid but this should be a coastal point (kg/dt)
    REAL(r_std), INTENT(out)       :: coastalflow(nbpt)    !! Outflow on coastal points by small basins. This is the water which flows in a disperse way into the ocean (kg/dt)
    REAL(r_std), INTENT(out)       :: flood_frac(nbpt)     !! Flooded fraction of the grid box (unitless;0-1)
    REAL(r_std), INTENT(out)       :: flood_res(nbpt)      !! Diagnostic of water amount in the floodplains reservoir (kg)

    !_ ================================================================================================================================

    !! 1. Call the main subroutine from the routing module corresponding to the choice of ROUTING_METHOD

    IF (routing_method=='subgrid_halfdeg') THEN

       CALL routing_subgrid_halfdeg_main (kjit, nbpt, index, &
            lalo, neighbours, resolution, contfrac, totfrac_nobio, veget_max, floodout, runoff, &
            drainage, transpot, precip_rain, humrel, k_litt, flood_frac, flood_res, &
            stempdiag, reinf_slope, returnflow, reinfiltration, irrigation, riverflow, coastalflow, rest_id, hist_id, hist2_id, &
            soiltile, root_deficit, irrigated_next, irrig_frac, fraction_aeirrig_sw)

    ELSE IF (routing_method=='subgrid_htu') THEN

       CALL routing_subgrid_htu_main (kjit, nbpt, index, &
            lalo, neighbours, resolution, contfrac, totfrac_nobio, veget_max, floodout, runoff, &
            drainage, transpot, precip_rain, humrel, k_litt, flood_frac, flood_res, &
            ftempdiag, reinf_slope, returnflow, reinfiltration, irrigation, riverflow, coastalflow, rest_id, hist_id, hist2_id)

    ELSE IF(routing_method=='interp_topo_v0') THEN 
#ifdef XIOS
       CALL routing_interp_topo_v0_main (kjit, nbpt, index, &
            lalo, neighbours, resolution, contfrac, totfrac_nobio, veget_max, floodout, runoff, &
            drainage, transpot, precip_rain, humrel, k_litt, flood_frac, flood_res, &
            stempdiag, reinf_slope, returnflow, reinfiltration, irrigation, riverflow, coastalflow, &
            rest_id, hist_id, hist2_id) 
#endif
    ELSE IF(routing_method=='interp_topo') THEN
#ifdef XIOS
       CALL routing_interp_topo_main (kjit, nbpt, index, &
            lalo, neighbours, resolution, contfrac, totfrac_nobio, veget_max, floodout, runoff, &
            drainage, transpot, precip_rain, humrel, k_litt, flood_frac, flood_res, &
            stempdiag, reinf_slope, returnflow, reinfiltration, irrigation, riverflow, coastalflow, &
            rest_id, hist_id, hist2_id, &
            soiltile, root_deficit, irrigated_next, irrig_frac, fraction_aeirrig_sw) 
#endif
    ENDIF

  END SUBROUTINE routing_wrapper_main


!!  =============================================================================================================================
!! SUBROUTINE:    routing_wrapper_finalize
!!
!>\BRIEF	  Call the finalization subroutine for the choosen routing module
!!
!! DESCRIPTION:	  Call the subroutine for finalization for the choosen routing module according to ROUTING_METHOD
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S): None
!! 
!! FLOWCHART: None
!! \n
!_ ==============================================================================================================================
  SUBROUTINE routing_wrapper_finalize( kjit, nbpt, rest_id, flood_frac, flood_res )

    IMPLICIT NONE
    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)     :: kjit                 !! Time step number (unitless)
    INTEGER(i_std), INTENT(in)     :: nbpt                 !! Domain size (unitless)
    INTEGER(i_std),INTENT(in)      :: rest_id              !! Restart file identifier (unitless)
    REAL(r_std), INTENT(in)        :: flood_frac(nbpt)     !! Flooded fraction of the grid box (unitless;0-1)
    REAL(r_std), INTENT(in)        :: flood_res(nbpt)      !! Diagnostic of water amount in the floodplains reservoir (kg)

    !_ ================================================================================================================================

    !! 1. Call the finalization subroutine from the routing module corresponding to the choice of ROUTING_METHOD

    IF (routing_method=='subgrid_halfdeg') THEN

       CALL routing_subgrid_halfdeg_finalize( kjit, nbpt, rest_id, flood_frac, flood_res )

    ELSE IF (routing_method=='subgrid_htu') THEN

       CALL routing_subgrid_htu_finalize( kjit, nbpt, rest_id, flood_frac, flood_res )

    ELSE IF(routing_method=='interp_topo_v0') THEN 
#ifdef XIOS
       CALL routing_interp_topo_v0_finalize( kjit, nbpt, rest_id, flood_frac, flood_res )
#endif
    ELSE IF(routing_method=='interp_topo') THEN
#ifdef XIOS
       CALL routing_interp_topo_finalize( kjit, nbpt, rest_id, flood_frac, flood_res )
#endif
    ENDIF

  END SUBROUTINE routing_wrapper_finalize


!!  =============================================================================================================================
!! SUBROUTINE:    routing_wrapper_clear
!!
!>\BRIEF	  Call the clear subroutine for the choosen routing module
!!
!! DESCRIPTION:	  Call the clear subroutine for the choosen routing module according to ROUTING_METHOD
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S): None
!! 
!! FLOWCHART: None
!! \n
!_ ==============================================================================================================================
  SUBROUTINE routing_wrapper_clear

    IF (routing_method=='subgrid_halfdeg') THEN

       CALL routing_subgrid_halfdeg_clear

    ELSE IF (routing_method=='subgrid_htu') THEN

       CALL routing_subgrid_htu_clear

    ELSE IF(routing_method=='interp_topo_v0') THEN 
#ifdef XIOS
       CALL routing_interp_topo_v0_clear
#endif    
    ELSE IF(routing_method=='interp_topo') THEN
#ifdef XIOS
       CALL routing_interp_topo_clear
#endif
    ENDIF

  END SUBROUTINE routing_wrapper_clear

END MODULE routing_wrapper
