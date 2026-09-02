MODULE routing_interp_topo_irrig_mod
#ifdef XIOS
  USE constantes

  
  PRIVATE
  
  INTEGER, SAVE :: nbpt
  !$OMP THREADPRIVATE(nbpt)
  
  INTEGER, SAVE :: nbpt_r
  !$OMP THREADPRIVATE(nbpt_r)

  REAL(r_std), ALLOCATABLE, SAVE :: irrig_gw_source_r(:)    ! diag, irrig withdrawal from gw reservoir (on routing grid)
  !$OMP THREADPRIVATE(irrig_gw_source_r)
  REAL(r_std), ALLOCATABLE, SAVE :: irrig_fast_source_r(:)  ! diag, irrig withdrawal from fast reservoir (on routing grid)
  !$OMP THREADPRIVATE(irrig_fast_source_r)
  REAL(r_std), ALLOCATABLE, SAVE :: irrig_str_source_r(:)   ! diag, irrig withdrawal from stream reservoir (on routing grid)
  !$OMP THREADPRIVATE(irrig_str_source_r)
 
  REAL(r_std), ALLOCATABLE, SAVE :: vegtot_mean(:)
  !$OMP THREADPRIVATE(vegtot_mean)
  
  REAL(r_std), ALLOCATABLE, SAVE :: humrel_mean(:)
  !$OMP THREADPRIVATE(humrel_mean)
  
  REAL(r_std), ALLOCATABLE, SAVE :: transpot_mean(:)
  !$OMP THREADPRIVATE(transpot_mean)
  
  REAL(r_std), ALLOCATABLE, SAVE :: runoff_mean(:)
  !$OMP THREADPRIVATE(runoff_mean)
  
  REAL(r_std), ALLOCATABLE, SAVE :: precip_mean(:)
  !$OMP THREADPRIVATE(precip_mean)

  REAL(r_std), ALLOCATABLE, SAVE :: irrigation_mean(:)  ! irrigation flux taken from the reservoirs to be put at the surface (kg/m^2/dt_routing) (on the ORC grid)
  !$OMP THREADPRIVATE(irrigation_mean)

  REAL(r_std), ALLOCATABLE, SAVE :: irrig_netereq_mean(:)  ! diag, irrigation requirement flux (kg/m^2/dt_routing) (on the ORC grid)
  !$OMP THREADPRIVATE(irrig_netereq_mean)
 
  REAL(r_std), ALLOCATABLE, SAVE :: irrigated(:)
  !$OMP THREADPRIVATE(irrigated)
 
  PUBLIC irrigation_initialize, irrigation_main, irrigation_mean_make, irrigation_get, irrigation_finalize

CONTAINS


  !! Public subroutine to ask for the module variables
  SUBROUTINE irrigation_get(routing_irrigation_mean, routing_irrig_netereq_mean)
  IMPLICIT NONE
    REAL(r_std),OPTIONAL, INTENT(OUT) :: routing_irrigation_mean(:)
    REAL(r_std),OPTIONAL, INTENT(OUT) :: routing_irrig_netereq_mean(:)

    CALL irrigation_get_(irrigation_mean_=routing_irrigation_mean, &
                         irrig_netereq_mean_=routing_irrig_netereq_mean)
  END SUBROUTINE irrigation_get

  !! Private subroutine (called by the public one) to get the module variables
  SUBROUTINE irrigation_get_(irrigation_mean_, irrig_netereq_mean_)
  IMPLICIT NONE
    REAL(r_std),OPTIONAL, INTENT(OUT) :: irrigation_mean_(:)
    REAL(r_std),OPTIONAL, INTENT(OUT) :: irrig_netereq_mean_(:)
    IF (PRESENT(irrigation_mean_))      irrigation_mean_ = irrigation_mean        !! assign local variable with the module variable irrigation_mean
    IF (PRESENT(irrig_netereq_mean_))   irrig_netereq_mean_ = irrig_netereq_mean  !! assign local variable with the module variable irrig_netereq_mean

  END SUBROUTINE irrigation_get_


  SUBROUTINE irrigation_initialize(kjit, rest_id, nbpt_, nbpt_r_, irrigated_next )
  USE constantes
  IMPLICIT NONE
    INTEGER, INTENT(IN)         :: kjit
    INTEGER, INTENT(IN)         :: rest_id
    INTEGER, INTENT(IN)         :: nbpt_
    INTEGER, INTENT(IN)         :: nbpt_r_
    REAL(r_std), INTENT(IN)     :: irrigated_next(nbpt_)
  
    CALL irrigation_local_init(kjit, rest_id, nbpt_, nbpt_r_,irrigated_next)
    CALL irrigation_mean_init(kjit, rest_id)
    
  END SUBROUTINE irrigation_initialize


  SUBROUTINE irrigation_finalize(kjit, rest_id)
  IMPLICIT NONE
    INTEGER, INTENT(IN)     :: kjit
    INTEGER, INTENT(IN)     :: rest_id
    
    CALL irrigation_mean_finalize(kjit, rest_id)
    
  END SUBROUTINE irrigation_finalize
  
  
  SUBROUTINE irrigation_local_init(kjit, rest_id, nbpt_, nbpt_r_, irrigated_next)
  IMPLICIT NONE
    INTEGER, INTENT(IN)         :: kjit
    INTEGER, INTENT(IN)         :: rest_id
    INTEGER, INTENT(IN)         :: nbpt_
    INTEGER, INTENT(IN)         :: nbpt_r_
    REAL(r_std), INTENT(IN)     :: irrigated_next(nbpt_)

     nbpt = nbpt_
     nbpt_r = nbpt_r_
     IF (do_irrigation) THEN
       ALLOCATE(irrigated(nbpt))
       irrigated(:)=irrigated_next(:)
       
       ! Diagnostics
       ALLOCATE(irrig_gw_source_r(nbpt_r))
       ALLOCATE(irrig_fast_source_r(nbpt_r))
       ALLOCATE(irrig_str_source_r(nbpt_r))
       
     ENDIF
     
  END SUBROUTINE irrigation_local_init


  
  SUBROUTINE irrigation_mean_init(kjit, rest_id)
   USE ioipsl_para
   USE grid
   USE sechiba_io_p
   IMPLICIT NONE
    INTEGER, INTENT(IN)     :: kjit
    INTEGER, INTENT(IN)     :: rest_id     !! Restart file identifier (unitless)
    CHARACTER(LEN=80)       :: var_name    !! To store variables names for I/O (unitless)
    INTEGER(i_std)          :: ier
    
    ! Get from the restart the flux we accumulated.
    ALLOCATE(irrigation_mean(nbpt), stat=ier)
    irrigation_mean(:) = 0
    IF (do_irrigation) THEN
      IF (ier /= 0) CALL ipslerr_p(3,'irrigation_mean_init','Pb in allocate for irrigation_mean','','')
      var_name = 'irrigation'
      CALL ioconf_setatt_p('UNITS', 'Kg/dt')
      CALL ioconf_setatt_p('LONG_NAME','Artificial irrigation flux')
      CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., irrigation_mean, "gather", nbp_glo, index_g)
      CALL setvar_p (irrigation_mean, val_exp, 'NO_KEYWORD', zero)
    
      ALLOCATE(vegtot_mean(nbpt), stat=ier)
      ALLOCATE(humrel_mean(nbpt), stat=ier)
      ALLOCATE(transpot_mean(nbpt), stat=ier)
      ALLOCATE(runoff_mean(nbpt), stat=ier)
      ALLOCATE(precip_mean(nbpt), stat=ier)
      CALL irrigation_mean_reset
    ENDIF
    
    ! Allocate diag (necessary even without do_irrigation, since it is used by routine irrigation_get)
    ALLOCATE(irrig_netereq_mean(nbpt), stat=ier)
    irrig_netereq_mean(:) = 0 ! maybe should be put and read in restart instead, like irrigation
                                ! but consistent with the reset of the other variables,
                                ! + it's only diagnostic

  END SUBROUTINE irrigation_mean_init

  
  SUBROUTINE irrigation_mean_reset
  IMPLICIT NONE

    vegtot_mean(:) = 0
    humrel_mean(:) = 0
    transpot_mean(:) = 0
    runoff_mean(:) = 0
    precip_mean(:) = 0

  END SUBROUTINE irrigation_mean_reset
  
  SUBROUTINE irrigation_mean_make(dt_routing, veget_max, humrel, transpot, runoff, precip )
  USE pft_parameters
  IMPLICIT NONE
    REAL(r_std),INTENT(IN) :: dt_routing
    REAL(r_std),INTENT(IN) :: veget_max(:,:)
    REAL(r_std),INTENT(IN) :: humrel(:,:)
    REAL(r_std),INTENT(IN) :: transpot(:,:)
    REAL(r_std),INTENT(IN) :: runoff(:)
    REAL(r_std),INTENT(IN) :: precip(:)

    REAL(r_std), DIMENSION(nbpt)   :: tot_vegfrac_nowoody  !! Total fraction occupied by grass (0-1,unitless)
    INTEGER :: jv, ig

    IF (do_irrigation) THEN
      runoff_mean(:) = runoff_mean(:) + runoff(:)
      precip_mean(:) = precip_mean(:) + precip(:)

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

      ! New irrigation scheme uses mean of vegtot with jv 1 to nvm
      ! Old scheme keeps jv 2 to nvm, even if possibly an error
      IF ( .NOT. old_irrig_scheme ) THEN
        DO jv=1,nvm
           DO ig=1,nbpt
              humrel_mean(ig) = humrel_mean(ig) + humrel(ig,jv)*veget_max(ig,jv)*dt_sechiba/dt_routing
              vegtot_mean(ig) = vegtot_mean(ig) + veget_max(ig,jv)*dt_sechiba/dt_routing
           ENDDO
        ENDDO
      ELSE
        DO jv=2,nvm
           DO ig=1,nbpt
              humrel_mean(ig) = humrel_mean(ig) + humrel(ig,jv)*veget_max(ig,jv)*dt_sechiba/dt_routing
              vegtot_mean(ig) = vegtot_mean(ig) + veget_max(ig,jv)*dt_sechiba/dt_routing
           ENDDO
        ENDDO
      ENDIF
    ENDIF !do_irrigation
  END SUBROUTINE irrigation_mean_make
 
  SUBROUTINE irrigation_mean_finalize(kjit, rest_id)
  USE grid
  USE ioipsl_para
  IMPLICIT NONE
    INTEGER, INTENT(IN)     :: kjit
    INTEGER, INTENT(IN)     :: rest_id

    IF (do_irrigation) THEN
      CALL restput_p (rest_id, 'irrigation', nbp_glo, 1, 1, kjit, irrigation_mean, 'scatter',  nbp_glo, index_g)  
    ENDIF

  END SUBROUTINE irrigation_mean_finalize
  

  SUBROUTINE irrigation_main(dt_routing, contfrac, reinfiltration, irrigated_next, irrig_frac, root_deficit, soiltile, &
                             fraction_aeirrig_sw )
  USE constantes
  USE constantes_soil
  USE grid, ONLY : area
  USE mod_orchidee_para, ONLY : nbp_mpi
  USE xios
  IMPLICIT NONE 
    REAL(r_std), INTENT(IN)       :: dt_routing               
    REAL(r_std), INTENT(IN)       :: contfrac(nbpt)               
    REAL(r_std), INTENT(IN)       :: reinfiltration(nbpt)               
    REAL(r_std), INTENT(IN)       :: irrigated_next(nbpt)           
    REAL(r_std), INTENT(IN)       :: irrig_frac(nbpt)     !! Irrig. fraction interpolated in routing, and saved to pass to slowproc if irrigated_soiltile = .TRUE.
    REAL(r_std), INTENT(IN)       :: root_deficit(nbpt)   !! soil water deficit
    REAL(r_std), INTENT(IN)       :: soiltile(nbpt,nstm)  !! Fraction of each soil tile within vegtot (0-1, unitless)
    REAL(r_std), INTENT(IN)       :: fraction_aeirrig_sw(nbpt) !! Fraction of area equipped for irrigation from surface water, of irrig_frac

    REAL(r_std)                   :: irrig_needs_r(nbpt_r)   !! Total irrigation requirement (water requirements by the crop for its optimal growth) (kg) (on routing grid)
    REAL(r_std)                   :: irrig_actual_r(nbpt_r)  !! Possible irrigation according to the water availability in the reservoirs (kg) (on routing grid)
    REAL(r_std)                   :: irrig_actual_mpi(nbp_mpi)  !! Possible irrigation according to the water availability - gathered on each MPI (kg) (on ORC grid)

    IF (do_irrigation) THEN
      IF (irrig_map_dynamic_flag) irrigated(:)=irrigated_next(:)
   
      IF (old_irrig_scheme) THEN
        CALL irrigation_compute_requested_old(contfrac, reinfiltration, irrigated, irrig_needs_r)
        CALL irrigation_old(irrig_needs_r, irrig_actual_r)
      ELSE
        CALL irrigation_compute_requested_new(dt_routing, contfrac, root_deficit, soiltile, irrig_frac, irrig_needs_r)
        CALL irrigation_new(contfrac,fraction_aeirrig_sw, irrig_needs_r, irrig_actual_r)
      ENDIF
      
      IF (is_omp_root) THEN
        CALL xios_send_field("routing_irrigation_r",irrig_actual_r)
        CALL xios_recv_field("routing_irrigation",irrig_actual_mpi)
      ENDIF
      CALL scatter_omp(irrig_actual_mpi,irrigation_mean)
      irrigation_mean(:)=irrigation_mean(:)/(area(:)*contfrac(:))
      CALL irrigation_mean_reset()
    ELSE
      irrigation_mean(:) = 0
    ENDIF
    
  END SUBROUTINE irrigation_main
  
  
 


  SUBROUTINE irrigation_new(contfrac, fraction_aeirrig_sw, irrig_needs_r, irrig_actual_r)
  USE constantes
  USE routing_interp_topo_flow_mod, ONLY : routing_mask_r, routing_flow_get, routing_flow_set 
  USE mod_orchidee_para
  USE xios
  IMPLICIT NONE
    REAL(r_std),INTENT(IN)        :: contfrac(nbpt)           
    REAL(r_std),INTENT(IN)        :: fraction_aeirrig_sw(nbpt)           
    REAL(r_std),INTENT(IN)        :: irrig_needs_r(nbpt_r)               !! Total irrigation requirement (water requirements by the crop for its optimal growth) (kg)
    REAL(r_std),INTENT(OUT)       :: irrig_actual_r(nbpt_r)              !! Possible irrigation according to the water availability in the reservoirs (kg)

    REAL(r_std)    :: irrig_deficit_r(nbpt_r)             !! Amount of water missing for irrigation (kg)

    LOGICAL        :: IsFail_slow        !! Logical to ask if slow reserv. does not fit irrigation demand
    LOGICAL        :: IsFail_fast        !! Logical to ask if fast reserv. does not fit irrigation demand
    LOGICAL        :: IsFail_stre        !! Logical to ask if stream reserv. does not fit irrigation demand

    REAL(r_std)    :: Count_failure_slow(nbpt_r)       !! Counter times slow reserv. does not fit irrigation demand
    REAL(r_std)    :: Count_failure_fast(nbpt_r)        !! Counter times fast reserv. does not fit irrigation demand
    REAL(r_std)    :: Count_failure_stre(nbpt_r)        !! Counter times stream reserv. does not fit irrigation demand
    
    REAL(r_std)    :: pot_slow_wdr_dummy, pot_fast_wdr_dummy, pot_stre_wdr_dummy  !! Maximum available water for irrigation for each reservoir (kg)
    REAL(r_std)    :: slow_wdr_dummy, fast_wdr_dummy, stre_wdr_dummy              !! Temporary variables for computation of diagnostics

    REAL(r_std)    :: slow_reservoir_r(nbpt_r)            !! Local slow_reservoir, that is retrieved from the routing and given back after update from irrigation (kg)
    REAL(r_std)    :: fast_reservoir_r(nbpt_r)            !! Local fast_reservoir, that is retrieved from the routing and given back after update from irrigation (kg)
    REAL(r_std)    :: stream_reservoir_r(nbpt_r)          !! Local stream_reservoir, that is retrieved from the routing and given back after update from irrigation (kg)
    REAL(r_std)    :: fraction_aeirrig_sw_mpi(nbp_mpi)
    REAL(r_std)    :: fraction_aeirrig_sw_r(nbpt_r)
    REAL(r_std)    :: pcent_vol_irrig                     !! Fraction of total avail. surface water / total avail water (surface+GW)
    INTEGER :: ig

    !! Initialize diagnostics
    irrig_gw_source_r(:) = 0
    irrig_fast_source_r(:) = 0
    irrig_str_source_r(:) = 0


    CALL gather_omp(fraction_aeirrig_sw, fraction_aeirrig_sw_mpi)
    !! fraction moyenne => need contfrac ??
      
    IF (is_omp_root) THEN
      
      CALL xios_send_field("fraction_aeirrig_sw", fraction_aeirrig_sw_mpi)
      CALL xios_recv_field("fraction_aeirrig_sw_r", fraction_aeirrig_sw_r) ! ==> need conservative quantity interp => no flux => to check !!!


      Count_failure_slow(:) = zero !
      Count_failure_fast(:) = zero !
      Count_failure_stre(:) = zero !

      CALL routing_flow_get(slow_reservoir_r=slow_reservoir_r, fast_reservoir_r=fast_reservoir_r, stream_reservoir_r=stream_reservoir_r)      

      DO ig=1,nbpt_r

        IF (routing_mask_r(ig)) THEN
        
          IF (select_source_irrig) THEN
            ! Priority order is :
            !   - for surface water, withdrawals in stream then fast
            !   - for ground water,  withdrawal in slow. GW withdrawal is disconnected from surface water withdrawal

            ! For irrig. scheme, avail_reserve gives the amount of water available for irrigation in every reservoir
            ! --> avail_reserve is a vector of dimension=3, BY DEFINITION i=1 for streamflow, i=2 fast, and i=3 slow reservoir

            ! The new priorization scheme (select_source_irrig) takes into account irrig. infrastructure according to GMIA map
            ! It also withdraws water according to availability (avail_reserve, between 0-1),
            ! meaning that it wont seek for all the water in the stream reservoir, even if this one could respond to the demand by itself

            pot_slow_wdr_dummy = ( 1 - fraction_aeirrig_sw_r(ig)) * avail_reserve(3) * slow_reservoir_r(ig)
            pot_fast_wdr_dummy = fraction_aeirrig_sw_r(ig) * avail_reserve(2) * fast_reservoir_r(ig)
            pot_stre_wdr_dummy = fraction_aeirrig_sw_r(ig) * avail_reserve(1) * stream_reservoir_r(ig)
            pcent_vol_irrig = zero
            
            IsFail_slow = .FALSE. !
            IsFail_fast = .FALSE. !
            IsFail_stre = .FALSE. !

            irrig_actual_r(ig) = MIN(irrig_needs_r(ig), pot_stre_wdr_dummy + pot_fast_wdr_dummy + pot_slow_wdr_dummy)
            ! irrig_actual set to zero if there is no available water.
            ! Put to avoid negative values
            irrig_actual_r(ig) = MAX(irrig_actual_r(ig), zero)

            !! pcent_vol_irrig corresponds to the fraction of available water in surface,
            !! considering environmental needs and irrigation equipment by source from map. It controls
            !! how the source of water withdrawal is distributed, especially when requirements < available water
            !! To calculate pcent_vol_irrig, we need to check that the total avail. water is not zero
            IF (  (pot_stre_wdr_dummy + pot_fast_wdr_dummy + pot_slow_wdr_dummy)  .GT. min_sechiba ) THEN

              pcent_vol_irrig = ( pot_stre_wdr_dummy + pot_fast_wdr_dummy ) / &
                                ( pot_stre_wdr_dummy + pot_fast_wdr_dummy + pot_slow_wdr_dummy )            
            !ELSE
              !Already zero because pcent_vol_irrig initialized to zero
              !Put here to readability but not necessary
              !    pcent_vol_irrig = zero
            ENDIF

            ! Update the reservoirs.
            ! Computation of irrig_gw_source(ig): first we add the slow_reservoir volume. Then we substract the updated slow_reservoir.
            ! What results is the volume used for irrigation that comes from GW
            ! Idem for irrig_fast_source and irrig_str_source

            slow_wdr_dummy = slow_reservoir_r(ig) ! slow reservoir before withdrawal for irrigation
            !! Slow reservoir update:
             !! we take all the relevant ground water ((1-pcent_vol_irrig)*irrig_actual_r) from the slow reservoir.
             !! if the resulting withdrawal exceeds the slow reservoir's maximum water availability,
             !! we only withdraw up to this availability
            slow_reservoir_r(ig) = MAX( slow_reservoir_r(ig) - pot_slow_wdr_dummy, &
                                        slow_reservoir_r(ig) - (un - pcent_vol_irrig) * irrig_actual_r(ig) )
            slow_wdr_dummy = slow_wdr_dummy - slow_reservoir_r(ig)
            irrig_gw_source_r(ig) = irrig_gw_source_r(ig) + slow_wdr_dummy

            fast_wdr_dummy = fast_reservoir_r(ig)
            !! Fast reservoir update:
             !! if the stream reservoir can not supply for all the surface water needs (pcent_vol_irrig*irrig_actual_r),
             !! we complete with the fast reservoir, up to the reservoir's maximum water availability
            fast_reservoir_r(ig) = MAX( fast_reservoir_r(ig) - pot_fast_wdr_dummy, &
                                        fast_reservoir_r(ig) + MIN(zero, pot_stre_wdr_dummy - pcent_vol_irrig * irrig_actual_r(ig) ) )
            fast_wdr_dummy = fast_wdr_dummy - fast_reservoir_r(ig)
            irrig_fast_source_r(ig) = irrig_fast_source_r(ig) + fast_wdr_dummy

            stre_wdr_dummy = stream_reservoir_r(ig)
            !! Stream reservoir update:
             !! we take all the relevant surface water (pcent_vol_irrig*irrig_actual_r) from the stream reservoir,
             !! except if it exceeds the reservoir's maximum water availability (then, see the fast reservoir update above)
            stream_reservoir_r(ig) = MAX( stream_reservoir_r(ig) - pot_stre_wdr_dummy, &
                                          stream_reservoir_r(ig) - pcent_vol_irrig * irrig_actual_r(ig) )
            stre_wdr_dummy = stre_wdr_dummy - stream_reservoir_r(ig)
            irrig_str_source_r(ig) = irrig_str_source_r(ig) + stre_wdr_dummy

            irrig_deficit_r(ig) = irrig_needs_r(ig)-irrig_actual_r(ig)

            ! A reservoir is failing to give water for irrigation if pot. req > pot. withdrawal
            ! We assume that the pot. requirement = Needs * fraction of area equipped for SW/GW
            ! In the case of surface, we also sustract the withdrawal from Fast/Stream, because both are
            ! considered as surface water

            IsFail_slow = ( ( irrig_needs_r(ig)*(un - fraction_aeirrig_sw_r(ig)) ) > pot_slow_wdr_dummy )
            IsFail_fast = ( ( irrig_needs_r(ig)*fraction_aeirrig_sw_r(ig) - stre_wdr_dummy ) > pot_fast_wdr_dummy )
            IsFail_stre = ( ( irrig_needs_r(ig)*fraction_aeirrig_sw_r(ig) - fast_wdr_dummy ) > pot_stre_wdr_dummy )

            IF( IsFail_stre ) Count_failure_stre(ig) = un
            IF( IsFail_fast ) Count_failure_fast(ig) = un
            IF( IsFail_slow ) Count_failure_slow(ig) = un


          ELSE IF (.NOT. select_source_irrig) THEN
            ! Priority order is :
            !   - withdrawal in stream
            !   - then withdrawal in fast
            !   - then withdrawal in slow

            ! For irrig. scheme, avail_reserve gives the amount of water available for irrigation in every reservoir
            ! --> avail_reserve is a vector of dimension=3, BY DEFINITION i=1 for streamflow, i=2 fast, and i=3 slow reservoir

            pot_slow_wdr_dummy = avail_reserve(3)*slow_reservoir_r(ig)
            pot_fast_wdr_dummy = avail_reserve(2)*fast_reservoir_r(ig)
            pot_stre_wdr_dummy = avail_reserve(1)*stream_reservoir_r(ig)
            IsFail_slow = .FALSE. !
            IsFail_fast = .FALSE. !
            IsFail_stre = .FALSE. !

            irrig_actual_r(ig) = MIN(irrig_needs_r(ig), pot_stre_wdr_dummy + pot_fast_wdr_dummy + pot_slow_wdr_dummy )
            ! irrig_actual set to zero if there is no available water.
            ! Put to avoid negative values
            irrig_actual_r(ig) = MAX(irrig_actual_r(ig), zero)

            ! Computation of irrig_gw_source(ig): first we add the slow_reservoir volume. Then we substract the updated slow_reservoir.
            ! It should result the volume used for irrigation that comes from GW
            ! Idem for irrig_fast_source and irrig_str_source

            slow_wdr_dummy = slow_reservoir_r(ig)
            !! Fast reservoir update:
             !! if the stream and fast reservoirs can not supply for all the water withdrawal (irrig_actual_r),
             !! we complete with the slow reservoir, up to the reservoir's maximum water availability
            slow_reservoir_r(ig) = MAX( slow_reservoir_r(ig) - pot_slow_wdr_dummy,       &
                                        slow_reservoir_r(ig)                             &
                                    + MIN(zero, pot_fast_wdr_dummy + MIN(zero, pot_stre_wdr_dummy-irrig_actual_r(ig))))
            slow_wdr_dummy = slow_wdr_dummy - slow_reservoir_r(ig)
            irrig_gw_source_r(ig) = irrig_gw_source_r(ig) + slow_wdr_dummy

            fast_wdr_dummy = fast_reservoir_r(ig)
            !! Fast reservoir update:
             !! if the stream reservoir can not supply for all the water withdrawal (irrig_actual_r),
             !! we complete with the fast reservoir, up to the reservoir's maximum water availability (then, see the slow reservoir update above)
            fast_reservoir_r(ig) = MAX( fast_reservoir_r(ig) - pot_fast_wdr_dummy ,  &
                                        fast_reservoir_r(ig)                         &
                                    + MIN(zero, pot_stre_wdr_dummy-irrig_actual_r(ig)))
            fast_wdr_dummy = fast_wdr_dummy - fast_reservoir_r(ig)
            irrig_fast_source_r(ig) = irrig_fast_source_r(ig) + fast_wdr_dummy

            stre_wdr_dummy = stream_reservoir_r(ig)
            !! Stream reservoir update:
             !! we take all the relevant water (irrig_actual_r) from the stream reservoir,
             !! except if it exceeds the reservoir's maximum water availability (then, see the fast reservoir update above)
            stream_reservoir_r(ig) = MAX( stream_reservoir_r(ig) - pot_stre_wdr_dummy, &
                                          stream_reservoir_r(ig) - irrig_actual_r(ig) )
            stre_wdr_dummy = stre_wdr_dummy - stream_reservoir_r(ig)
            irrig_str_source_r(ig) = irrig_str_source_r(ig) + stre_wdr_dummy

            irrig_deficit_r(ig) = irrig_needs_r(ig)-irrig_actual_r(ig)

            ! A reservoir is failing to give water for irrigation if pot. req > pot. withdrawal
            ! Because it follows the old scheme, we do not separate between surface/gw, but consider that
            ! priority is given in this order: River (Stream), then Fast, then Slow reservoirs.

            IsFail_slow = ( ( irrig_needs_r(ig) - stre_wdr_dummy - fast_wdr_dummy  ) > pot_slow_wdr_dummy )
            IsFail_fast = ( ( irrig_needs_r(ig) - stre_wdr_dummy ) > pot_fast_wdr_dummy )
            IsFail_stre = ( irrig_needs_r(ig) > pot_stre_wdr_dummy )

            IF( IsFail_stre ) Count_failure_stre(ig) = un
            IF( IsFail_fast ) Count_failure_fast(ig) = un
            IF( IsFail_slow ) Count_failure_slow(ig) = un

          ENDIF !IF (select_source_irrig) 
        ELSE
          irrig_actual_r(ig) = 0
        ENDIF !IF (routing_mask_r(ig))
      ENDDO
  
      CALL routing_flow_set(slow_reservoir_r=slow_reservoir_r, fast_reservoir_r=fast_reservoir_r, stream_reservoir_r=stream_reservoir_r)

    ENDIF !IF (is_omp_root)
   
  END SUBROUTINE irrigation_new


  SUBROUTINE irrigation_old(irrig_needs_r, irrig_actual_r)
  USE constantes
  USE mod_orchidee_para
  USE routing_interp_topo_flow_mod, ONLY : routing_mask_r, routing_flow_get, routing_flow_set 
  IMPLICIT NONE
    REAL(r_std),INTENT(IN)        :: irrig_needs_r(nbpt_r)               !! Total irrigation requirement (water requirements by the crop for its optimal growth) (kg)
    REAL(r_std),INTENT(OUT)       :: irrig_actual_r(nbpt_r)             !! Possible irrigation according to the water availability in the reservoirs (kg)
 
    REAL(r_std)    :: irrig_deficit_r(nbpt_r)             !! Amount of water missing for irrigation (kg)

    LOGICAL        :: IsFail_slow        !! Logical to ask if slow reserv. does not fit irrigation demand
    LOGICAL        :: IsFail_fast        !! Logical to ask if fast reserv. does not fit irrigation demand
    LOGICAL        :: IsFail_stre        !! Logical to ask if stream reserv. does not fit irrigation demand

    REAL(r_std)    :: Count_failure_slow(nbpt_r)       !! Counter times slow reserv. does not fit irrigation demand
    REAL(r_std)    :: Count_failure_fast(nbpt_r)        !! Counter times fast reserv. does not fit irrigation demand
    REAL(r_std)    :: Count_failure_stre(nbpt_r)        !! Counter times stream reserv. does not fit irrigation demand
    
    REAL(r_std)    :: pot_slow_wdr_dummy, pot_fast_wdr_dummy, pot_stre_wdr_dummy
    REAL(r_std)    :: slow_wdr_dummy, fast_wdr_dummy, stre_wdr_dummy

    REAL(r_std)    :: slow_reservoir_r(nbpt_r)
    REAL(r_std)    :: fast_reservoir_r(nbpt_r)
    REAL(r_std)    :: stream_reservoir_r(nbpt_r)

    INTEGER :: ig

    !! Initialize diagnostics
    irrig_gw_source_r(:) = 0
    irrig_fast_source_r(:) = 0
    irrig_str_source_r(:) = 0

    IF (is_omp_root) THEN

      Count_failure_slow(:) = zero !
      Count_failure_fast(:) = zero !
      Count_failure_stre(:) = zero !

      CALL routing_flow_get(slow_reservoir_r=slow_reservoir_r, fast_reservoir_r=fast_reservoir_r, stream_reservoir_r=stream_reservoir_r)      

      DO ig=1,nbpt_r

        IF (routing_mask_r(ig)) THEN

          ! Old irrigation scheme as in tag 2.0
          ! Note for irrig_gw_source(ig): first we add the slow_reservoir volume. Then we substract the updated slow_reservoir. It should be the
          ! Volume used for irrigation that comes from GW
          ! Idem for irrig_fast_source and irrig_str_source

          pot_slow_wdr_dummy = slow_reservoir_r(ig)
          pot_fast_wdr_dummy = fast_reservoir_r(ig)
          pot_stre_wdr_dummy = stream_reservoir_r(ig)
          IsFail_slow = .FALSE. !
          IsFail_fast = .FALSE. !
          IsFail_stre = .FALSE. !
           
          irrig_actual_r(ig) = MIN(irrig_needs_r(ig), stream_reservoir_r(ig) + fast_reservoir_r(ig) + slow_reservoir_r(ig) )

          slow_wdr_dummy = slow_reservoir_r(ig)
          slow_reservoir_r(ig) = MAX(zero, slow_reservoir_r(ig) + MIN(zero, fast_reservoir_r(ig) &
                                 + MIN(zero, stream_reservoir_r(ig)-irrig_actual_r(ig))))
          slow_wdr_dummy = slow_wdr_dummy - slow_reservoir_r(ig)
          irrig_gw_source_r(ig) = irrig_gw_source_r(ig) + slow_wdr_dummy

          fast_wdr_dummy = fast_reservoir_r(ig)
          fast_reservoir_r(ig) = MAX( zero, fast_reservoir_r(ig) + MIN(zero, stream_reservoir_r(ig)-irrig_actual_r(ig)))
          fast_wdr_dummy = fast_wdr_dummy - fast_reservoir_r(ig)
          irrig_fast_source_r(ig) = irrig_fast_source_r(ig) + fast_wdr_dummy

          stre_wdr_dummy = stream_reservoir_r(ig)
          stream_reservoir_r(ig) = MAX(zero, stream_reservoir_r(ig)-irrig_actual_r(ig) )
          stre_wdr_dummy = stre_wdr_dummy - stream_reservoir_r(ig)
          irrig_str_source_r(ig) = irrig_str_source_r(ig) + stre_wdr_dummy

          irrig_deficit_r(ig) = irrig_needs_r(ig)-irrig_actual_r(ig)

          ! A reservoir is failing to give water for infiltration if pot. req > pot. withdrawal
          ! Because it follows the old scheme, we do not separate between surface/gw, but consider that
          ! priority is given in this order: River, Fast and Slow reservoir.

          IsFail_slow = ( ( irrig_needs_r(ig) - stre_wdr_dummy - fast_wdr_dummy  ) > pot_slow_wdr_dummy )
          IsFail_fast = ( ( irrig_needs_r(ig) - stre_wdr_dummy ) > pot_fast_wdr_dummy )
          IsFail_stre = ( irrig_needs_r(ig) > pot_stre_wdr_dummy )

          IF( IsFail_stre ) Count_failure_stre(ig) = un
          IF( IsFail_fast ) Count_failure_fast(ig) = un
          IF( IsFail_slow ) Count_failure_slow(ig) = un
        ELSE
          irrig_actual_r(ig) = 0
        ENDIF
      ENDDO
  
      CALL routing_flow_set(slow_reservoir_r=slow_reservoir_r, fast_reservoir_r=fast_reservoir_r, stream_reservoir_r=stream_reservoir_r)

    ENDIF
  
  END SUBROUTINE irrigation_old  
    
 
  SUBROUTINE irrigation_compute_requested_new(dt_routing, contfrac, root_deficit, soiltile, irrig_frac, irrig_netereq_r)
  USE constantes
  USE constantes_soil
  USE xios
  USE mod_orchidee_para
  USE grid, ONLY : area
  IMPLICIT NONE
    REAL(r_std),INTENT(IN)   :: dt_routing          
    REAL(r_std),INTENT(IN)   :: contfrac(nbpt)          
    REAL(r_std),INTENT(IN)   :: root_deficit(nbpt)               
    REAL(r_std),INTENT(IN)   :: soiltile(nbpt,nstm)  !! Fraction of each soil tile within vegtot (0-1, unitless)
    REAL(r_std),INTENT(IN)   :: irrig_frac(nbpt)      !! Irrig. fraction interpolated in routing, and saved to pass to slowproc if irrigated_soiltile = .TRUE.
    REAL(r_std),INTENT(OUT)  :: irrig_netereq_r(nbpt_r)  !! Irrigation requirement on the routing grid (kg/dt_routing)

    REAL(r_std)              :: irrig_netereq(nbpt)   !! Irrigation requirement flux on the ORC grid (water requirements by the plants for their optimal growth) (kg/m^2 of grid cell/dt_routing)
    REAL(r_std)              :: irrig_netereq_mpi(nbp_mpi)  !! Irrigation requirement on the ORC grid - gathered on each MPI (kg/dt_routing)
   
    INTEGER  :: ig
   
    irrig_netereq(:)=0
   
    DO ig=1,nbpt
      !It enters to the new irrigation module only if there is an irrigated fraction, if not irrig_netereq = zero for that cell
      IF ( irrig_frac(ig) .GT. min_sechiba ) THEN

        irrig_netereq(ig) = irrig_netereq(ig) + MIN( irrig_dosmax/3600*dt_routing, &
                            root_deficit(ig) ) * soiltile(ig, irrig_st) * vegtot_mean(ig)
       ! By definition, irrig_dosmax is in kg/m2 of soil tile/hour,dividing by 3600(seconds/hour) * DT_ROUTING  !
       ! = kg/m2 of soil tile/(routing timestep)
       ! irrig_netereq(kg/m2 of grid cell / routing timestep ) is equal to
       ! root_deficit (kg/m2 of soil tile) * soiltile*vegtot (fraction of soil tile at cell level) = kg/m2 of grid cell

       IF (.NOT. irrigated_soiltile .AND. ( soiltile(ig,irrig_st) .GT. min_sechiba ) .AND. (vegtot_mean(ig) .GT. min_sechiba) ) THEN
         ! Irrigated_soiltile asks if there is an independent soil tile for irrigated crops. If not,
         ! actual volume calculated for irrig_netereq assumed that the whole SOILTILE was irrigated, but in this case
         ! just a fraction of the irrig_st (irrigated soil tile, by default = 3) is actually irrigated,
         ! and irrig_netereq needs to be reduced by irrig_frac/( soiltile * vegtot ) (note that it is max = 1 thanks to irrig_frac calculation in l. 424)
         ! Demand(ST3)*irrig_frac/(soiltile(3)*vegtot) = irrig_netereq_In_ST3, then
         !irrig_netereq_In_ST3 * (soiltile(3)*vegtot) = irrig_netereq at grid scale = Demand(ST3)*irrig_frac.
         irrig_netereq(ig) = irrig_netereq(ig) * irrig_frac(ig) / ( soiltile(ig,irrig_st) * vegtot_mean(ig) )
          !irrig_netereq = kg/m2 of grid cell
        ENDIF
      ENDIF
    ENDDO
    
    ! Update diagnostics for irrigation requirement
    irrig_netereq_mean(:)=irrig_netereq(:)
    
    ! Interpolation from ORCHIDEE grid to the routing grid
    CALL gather_omp(irrig_netereq*area*contfrac, irrig_netereq_mpi)
    
    IF (is_omp_root) THEN
      CALL xios_send_field("irrig_netereq",irrig_netereq_mpi) !! contfrac => ok ?
      CALL xios_recv_field("irrig_netereq_r",irrig_netereq_r) ! ==> need conservative quantity interp
    ENDIF
      
  END SUBROUTINE irrigation_compute_requested_new
  
 
  SUBROUTINE irrigation_compute_requested_old(contfrac, reinfiltration, irrigated, irrig_need_r)
  USE constantes
  USE xios
  USE mod_orchidee_para
  USE grid, ONLY : area
  IMPLICIT NONE
    REAL(r_std),INTENT(IN)   :: contfrac(nbpt)          
    REAL(r_std),INTENT(IN)   :: reinfiltration(nbpt)               
    REAL(r_std),INTENT(IN)   :: irrigated(nbpt)               
    REAL(r_std),INTENT(OUT)  :: irrig_need_r(nbpt_r)

    REAL(r_std)              :: irrig_netereq(nbpt)
    REAL(r_std)              :: irrig_netereq_mpi(nbp_mpi)
   
    INTEGER  :: ig, ir
  
      irrig_netereq(:)=0
  
      DO ig=1,nbpt
        IF ((vegtot_mean(ig) .GT. min_sechiba) .AND. (humrel_mean(ig) .LT. un-min_sechiba) .AND. (runoff_mean(ig) .LT. min_sechiba) ) THEN
          irrig_netereq(ig) = (irrigated(ig)/area(ig)) * MAX(zero, transpot_mean(ig) - (precip_mean(ig)+reinfiltration(ig))) ! ==> ok kg
        ELSE
          irrig_netereq(ig) = 0
        ENDIF
      ENDDO
          
      ! Update diagnostics for irrigation requirement
      irrig_netereq_mean(:)=irrig_netereq(:)
    
      ! Interpolation from ORCHIDEE grid to the routing grid
      CALL gather_omp(irrig_netereq*area*contfrac, irrig_netereq_mpi)
     
      IF (is_omp_root) THEN
        CALL xios_send_field("irrig_netereq",irrig_netereq_mpi) ! contfrac ??
        CALL xios_recv_field("irrig_netereq_r",irrig_need_r) ! ==> need conservative quantity interp
      ENDIF
      
  END SUBROUTINE irrigation_compute_requested_old



!for now  SUBROUTINE abduction
!for now    REAL(r_std), DIMENSION(nbpt)        :: irrig_adduct              !! Amount of water carried over from other basins for irrigation (kg)
!for now  
!for now    DO ig=1,nbpt

!for now      !
!for now      ! Check if we cannot find the missing water in another basin of the same grid (stream reservoir only).
!for now      ! If we find that then we create some adduction from that subbasin to the one where we need it for
!for now      ! irrigation.
!for now      !

!for now      !> If crops water requirements have not been supplied (irrig_deficit>0), we check if we cannot find the missing water
!for now      !> in another basin of the same grid. If there is water in the stream reservoir of this subbasin, we create some adduction
!for now      !> from that subbasin to the one where we need it for irrigation.
!for now      !>
!for now    
!for now      DO ib=1,nbasmax
!for now        stream_tot = a_stream_adduction * SUM(stream_reservoir(ig,:))

!for now        DO WHILE ( irrig_deficit(ig,ib) > min_sechiba .AND. stream_tot > min_sechiba)
!for now          fi = MAXLOC(stream_reservoir(ig,:))
!for now          ib2 = fi(1)

!for now          irrig_adduct(ig,ib) = MIN(irrig_deficit(ig,ib), a_stream_adduction * stream_reservoir(ig,ib2))
!for now          stream_reservoir(ig,ib2) = stream_reservoir(ig,ib2)-irrig_adduct(ig,ib)
!for now          irrig_deficit(ig,ib) = irrig_deficit(ig,ib)-irrig_adduct(ig,ib)
!for now          stream_tot = a_stream_adduction * SUM(stream_reservoir(ig,:))

!for now        ENDDO
!for now      ENDDO
!for now
!for now    ENDDO

!for now    !
!for now    ! If we are at higher resolution we might need to look at neighboring grid boxes to find the streams
!for now    ! which can feed irrigation
!for now    !
!for now    !> At higher resolution (grid box smaller than 100x100km), we can import water from neighboring grid boxes
!for now    !> to the one where we need it for irrigation.
!for now    !
!for now    IF (is_root_prc) THEN
!for now      ALLOCATE(irrig_deficit_glo(nbp_glo, nbasmax), stream_reservoir_glo(nbp_glo, nbasmax), &
!for now               irrig_adduct_glo(nbp_glo, nbasmax), stat=ier)
!for now    ELSE
!for now      ALLOCATE(irrig_deficit_glo(0, 0), stream_reservoir_glo(0, 0), irrig_adduct_glo(0, 0), stat=ier)
!for now    ENDIF
!for now    IF (ier /= 0) CALL ipslerr_p(3,'routing_flow','Pb in allocate for irrig_deficit_glo, stream_reservoir_glo,...','','')

!for now    CALL gather(irrig_deficit, irrig_deficit_glo)
!for now    CALL gather(stream_reservoir,  stream_reservoir_glo)
!for now    CALL gather(irrig_adduct, irrig_adduct_glo)

!for now    IF (is_root_prc) THEN
!for now      !
!for now      DO ig=1,nbp_glo
!for now        ! Only work if the grid box is smaller than 100x100km. Else the piplines we build
!for now        ! here would be too long to be reasonable.
!for now        IF ( resolution_g(ig,1) < 100000. .AND. resolution_g(ig,2) < 100000. ) THEN
!for now          
!for now          DO ib=1,nbasmax
!for now          !
!for now            IF ( irrig_deficit_glo(ig,ib)  > min_sechiba ) THEN
!for now            !
!for now              streams_around(:,:) = zero
!for now              !
!for now              DO in=1,NbNeighb
!for now                ig2 = neighbours_g(ig,in)
!for now                IF (ig2 .GT. 0 ) THEN
!for now                  streams_around(in,:) = a_stream_adduction * stream_reservoir_glo(ig2,:)
!for now                  igrd(in) = ig2
!for now                ENDIF
!for now              ENDDO
!for now              !
!for now              IF ( MAXVAL(streams_around) .GT. zero ) THEN
!for now              !
!for now                ff=MAXLOC(streams_around)
!for now                ig2=igrd(ff(1))
!for now                ib2=ff(2)
!for now                !
!for now                IF ( routing_area_glo(ig2,ib2) .GT. 0 .AND. a_stream_adduction * stream_reservoir_glo(ig2,ib2) > zero ) THEN
!for now                  adduction = MIN(irrig_deficit_glo(ig,ib), a_stream_adduction * stream_reservoir_glo(ig2,ib2))
!for now                  stream_reservoir_glo(ig2,ib2) = stream_reservoir_glo(ig2,ib2) - adduction
!for now                  irrig_deficit_glo(ig,ib) = irrig_deficit_glo(ig,ib) - adduction
!for now                  irrig_adduct_glo(ig,ib) = irrig_adduct_glo(ig,ib) + adduction
!for now                ENDIF
!for now                !
!for now              ENDIF
!for now              !
!for now            ENDIF
!for now            !
!for now          ENDDO
!for now      
!for now        ENDIF

!for now      ENDDO
!for now          !
!for now    ENDIF

!for now    CALL scatter(irrig_deficit_glo, irrig_deficit)
!for now    CALL scatter(stream_reservoir_glo,  stream_reservoir)
!for now    CALL scatter(irrig_adduct_glo, irrig_adduct)

!for now    DEALLOCATE(irrig_deficit_glo, stream_reservoir_glo, irrig_adduct_glo)  
!for now  
!for now  
!for now  END SUBROUTINE abduction
  

#endif
END MODULE routing_interp_topo_irrig_mod
