MODULE routing_interp_topo_lake_mod
#ifdef XIOS
  USE constantes
   
  PRIVATE
  
  PUBLIC :: routing_lake_initialize, routing_lake_mean_make, routing_lake_main, routing_lake_route_coast, routing_lake_finalize, &
            routing_lake_get_water_balance
  
  REAL(r_std), SAVE, ALLOCATABLE :: lake_reservoir(:)
!$OMP THREADPRIVATE(lake_reservoir)

  REAL(r_std), SAVE, ALLOCATABLE :: humrel_mean(:)
!$OMP THREADPRIVATE(humrel_mean)
  
  INTEGER,SAVE  :: nbpt
!$OMP THREADPRIVATE( nbpt)
 
  LOGICAL,SAVE :: do_swamps
!$OMP THREADPRIVATE(do_swamps)

  INTEGER,SAVE :: nb_coast_cells
!$OMP THREADPRIVATE(nb_coast_cells)

REAL(r_std), PARAMETER   :: maxevap_lake = 7.5/86400.   !! Maximum evaporation rate from lakes (kg/m^2/s)

REAL(r_std), SAVE    :: max_lake_reservoir           !! Maximum limit of water in lake_reservoir [kg/m2]
!$OMP THREADPRIVATE(max_lake_reservoir)
  
LOGICAL, ALLOCATABLE :: is_coastline(:)
!$OMP THREADPRIVATE(is_coastline)  

CONTAINS
 
 SUBROUTINE routing_lake_initialize(kjit, rest_id, nbpt_, contfrac)
 IMPLICIT NONE
   INTEGER, INTENT(IN)     :: kjit
   INTEGER, INTENT(IN)     :: rest_id
   INTEGER,INTENT(IN)      :: nbpt_
   REAL(r_std), INTENT(IN)    :: contfrac(nbpt)     !! Fraction of land in each grid box (unitless;0-1)


   CALL routing_lake_init_local(kjit, rest_id, nbpt_, contfrac)
   CALL routing_lake_mean_init(kjit, rest_id)

 END SUBROUTINE routing_lake_initialize


 SUBROUTINE routing_lake_finalize(kjit, rest_id)
 IMPLICIT NONE
   INTEGER, INTENT(IN)     :: kjit
   INTEGER, INTENT(IN)     :: rest_id

   CALL routing_lake_finalize_local(kjit, rest_id)
   CALL routing_lake_mean_finalize(kjit, rest_id)

 END SUBROUTINE routing_lake_finalize

 
 SUBROUTINE routing_lake_init_local(kjit, rest_id, nbpt_, contfrac)
 USE mod_orchidee_para, ONLY : reduce_sum, bcast
 USE ioipsl_para
 USE grid, ONLY : nbp_glo, index_g
 USE sechiba_io_p
 USE routing_interp_topo_flow_mod, ONLY : compute_coastline
 IMPLICIT NONE
   INTEGER, INTENT(IN) :: kjit
   INTEGER, INTENT(IN) :: rest_id
   INTEGER,INTENT(IN) :: nbpt_
   REAL(r_std), INTENT(IN)    :: contfrac(nbpt)     !! Fraction of land in each grid box (unitless;0-1)
   
   INTEGER :: ier    
   CHARACTER(LEN=80)   :: var_name       !! To store variables names for I/O (unitless)
 
   INTEGER  :: ig
   INTEGER  :: nb_cells
   
   nbpt=nbpt_
   
   
    ALLOCATE(lake_reservoir(nbpt),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'routing_lake_init_local','Pb in allocate for lake_reservoir','','')
    var_name = 'lakeres'
    CALL ioconf_setatt_p('UNITS', 'Kg')
    CALL ioconf_setatt_p('LONG_NAME','Water in the lake reservoir')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., lake_reservoir, "gather", nbp_glo, index_g)
    CALL setvar_p (lake_reservoir, val_exp, 'NO_KEYWORD', zero)
       
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
    do_swamps = .FALSE.
    CALL getin_p('DO_SWAMPS', do_swamps)

    !Config Key   = MAX_LAKE_RESERVOIR
    !Config Desc  = Maximum limit of water in lake_reservoir
    !Config If    = RIVER_ROUTING
    !Config Def   = 7000
    !Config Help  =
    !Config Units = [kg/m2(routing area)]
    max_lake_reservoir = 7000
    CALL getin_p("MAX_LAKE_RESERVOIR", max_lake_reservoir)
        
    
    ! compute  number of coast cells
    ! WARNING : contfrac is fraction of ter, but 1-contfrac is not fraction of ocean because of landice
    ! => to be corrected...
    ALLOCATE(is_coastline(nbpt))
    CALL compute_coastline(contfrac, is_coastline)
    nb_cells=0
    
    DO ig=1,nbpt
      IF (is_coastline(ig)) nb_cells=nb_cells+1
    ENDDO
    
    CALL reduce_sum(nb_cells, nb_coast_cells)
    CALL bcast(nb_coast_cells)
  
  END SUBROUTINE routing_lake_init_local


  SUBROUTINE routing_lake_finalize_local(kjit, rest_id)
  USE ioipsl_para
  USE grid, ONLY : nbp_glo, index_g
  IMPLICIT NONE
    INTEGER, INTENT(IN) :: kjit
    INTEGER, INTENT(IN) :: rest_id
 
    CALL restput_p (rest_id, 'lakeres', nbp_glo, 1, 1, kjit, lake_reservoir, 'scatter',  nbp_glo, index_g)
    DEALLOCATE(lake_reservoir)

  END SUBROUTINE  routing_lake_finalize_local
  
  

 SUBROUTINE routing_lake_mean_init(kjit, rest_id)
 USE ioipsl_para
 USE grid
 USE sechiba_io_p
 IMPLICIT NONE
    INTEGER, INTENT(IN) :: kjit
    INTEGER, INTENT(IN) :: rest_id

    INTEGER :: ier    
    CHARACTER(LEN=80)   :: var_name       !! To store variables names for I/O (unitless)
  
    ALLOCATE(humrel_mean(nbpt), stat=ier)

    IF (ier /= 0) CALL ipslerr_p(3,'routing_lake_mean_init','Pb in allocate for humrel_mean','','')
    var_name = 'humrel_lake'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Mean humrel for irrigation')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., humrel_mean, "gather", nbp_glo, index_g)
    CALL setvar_p (humrel_mean, val_exp, 'NO_KEYWORD', un)
   
 END SUBROUTINE  routing_lake_mean_init
 
 SUBROUTINE routing_lake_mean_make(dt_routing, humrel, veget_max)
 USE constantes
 USE pft_parameters
 IMPLICIT NONE
   REAL(r_std), INTENT(IN)    :: dt_routing 
   REAL(r_std), INTENT(IN)    :: humrel(:,:)       !! Soil moisture stress, root extraction potential (unitless)
   REAL(r_std), INTENT(IN)    :: veget_max(:,:)       !! Soil moisture stress, root extraction potential (unitless)
   INTEGER :: jv
   
   IF ( .NOT. old_irrig_scheme ) THEN
     DO jv=1,nvm
       humrel_mean(:) = humrel_mean(:) + humrel(:,jv)*veget_max(:,jv)*dt_sechiba/dt_routing
     ENDDO
   ELSE
     DO jv=2,nvm
       humrel_mean(:) = humrel_mean(:) + humrel(:,jv)*veget_max(:,jv)*dt_sechiba/dt_routing
     ENDDO
   ENDIF
    
 END SUBROUTINE routing_lake_mean_make  

  
 SUBROUTINE routing_lake_mean_reset
 IMPLICIT NONE
   humrel_mean(:) = 0
 END SUBROUTINE routing_lake_mean_reset 

 
 SUBROUTINE routing_lake_mean_finalize(kjit, rest_id)
 USE ioipsl_para
 USE grid
 IMPLICIT NONE
   INTEGER, INTENT(IN) :: kjit
   INTEGER, INTENT(IN) :: rest_id
   
   CALL restput_p (rest_id, 'humrel_lake', nbp_glo, 1, 1, kjit, humrel_mean, 'scatter',  nbp_glo, index_g)
   DEALLOCATE(humrel_mean)
 
 END SUBROUTINE routing_lake_mean_finalize
 
 
 SUBROUTINE routing_lake_route_coast(contfrac, coastalflow)
 USE grid, ONLY : area
 USE mod_orchidee_para, ONLY : reduce_sum, bcast

 IMPLICIT NONE
   REAL(r_std), INTENT(IN)    :: contfrac(nbpt)     !! Fraction of land in each grid box (unitless;0-1)
   REAL(r_std), INTENT(INOUT) :: coastalflow(nbpt)   !! Water inflow to the lakes (kg/dt)
   
   REAL(r_std) :: sum_lake_overflow, lake_overflow, total_lake_overflow
   INTEGER :: ig
    
   !! Remove water from lake reservoir if it exceeds the maximum limit and distribute it 
   !! uniformly over all possible the coastflow gridcells
    
   ! Calculate lake_overflow and remove it from lake_reservoir
    sum_lake_overflow=0
    DO ig=1,nbpt
      lake_overflow = MAX(0., lake_reservoir(ig) - max_lake_reservoir*area(ig)*contfrac(ig))
      lake_reservoir(ig) = lake_reservoir(ig) - lake_overflow
      sum_lake_overflow = sum_lake_overflow + lake_overflow
    END DO

   ! Calculate the sum of the lake_overflow and distribute it uniformly over all gridboxes
    CALL reduce_sum(sum_lake_overflow, total_lake_overflow)
    CALL bcast(total_lake_overflow)

    WHERE (is_coastline) coastalflow = coastalflow + total_lake_overflow/nb_coast_cells 
 
 END SUBROUTINE routing_lake_route_coast



 
 SUBROUTINE routing_lake_main(dt_routing, contfrac, lakeinflow, return_lakes)

    USE grid, ONLY : area
    
    IMPLICIT NONE
    !! 0 Variable and parameter description
    !! 0.1 Input variables

    REAL(r_std), INTENT(IN)    :: dt_routing         !! Routing time step (s)
    REAL(r_std), INTENT(IN)    :: lakeinflow(nbpt)   !! Water inflow to the lakes (kg/dt)
    REAL(r_std), INTENT(IN)    :: contfrac(nbpt)     !! Fraction of land in each grid box (unitless;0-1)
    
    !! 0.2 Output variables
    REAL(r_std), INTENT(OUT)   :: return_lakes(nbpt) !! Water from lakes flowing back into soil moisture (kg/m^2/dt)
    
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
       !uptake in Kg/dt
       IF ( do_swamps ) THEN
         ! Calculate a return flow that will be extracted from the lake reservoir and reinserted in the soil in hydrol
         ! Uptake in Kg/dt
         refill = MAX(zero, maxevap_lake * (un - humrel_mean(ig)) * dt_routing * total_area)
         return_lakes(ig) = MIN(refill, lake_reservoir(ig))
         lake_reservoir(ig) = lake_reservoir(ig) - return_lakes(ig)
         !Return in Kg/m^2/dt
          return_lakes(ig) = return_lakes(ig)/total_area
       ELSE
          return_lakes(ig) = zero
       ENDIF
       !
       ! This is the volume of the lake scaled to the entire grid.
       ! It would be batter to scale it to the size of the lake
       ! but this information is not yet available.
!ym for now       lake_diag(ig) = lake_reservoir(ig)/total_area
       !
    ENDDO

    CALL routing_lake_mean_reset
 
  END SUBROUTINE routing_lake_main

  FUNCTION routing_lake_get_water_balance(contfrac)
  USE grid, ONLY : area
  USE mod_orchidee_para, ONLY : is_omp_root
  IMPLICIT NONE
    REAL(r_std) :: contfrac(nbpt)
    REAL(r_std) :: routing_lake_get_water_balance   !  Amount of water (local) in kg
    
    routing_lake_get_water_balance = SUM(lake_reservoir)
  
  END FUNCTION routing_lake_get_water_balance

#endif
END MODULE routing_interp_topo_lake_mod
