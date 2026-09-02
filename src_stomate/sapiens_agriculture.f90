! ================================================================================================================================
! MODULE       : sapiens_agriculture
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       Simple approach to harvesting croplands
!!
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S) : None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_stomate/sapiens_agriculture.f90 $
!! $Date: 2026-01-14 12:52:55 +0100 (mer. 14 janv. 2026) $
!! $Revision: 9305 $
!! \n
!_ ================================================================================================================================

MODULE sapiens_agriculture

  ! modules used:
  USE stomate_data
  USE constantes
  USE grid
  USE pft_parameters
  USE function_library,    ONLY: cc_to_biomass, lai_to_biomass, &
                           biomass_to_lai, get_printlev

  IMPLICIT NONE

  ! private & public routines

  PRIVATE
  PUBLIC crop_planting, crop_harvest, sapiens_agriculture_initialize

  ! Variable declaration 

  LOGICAL, SAVE             :: firstcall_sapiens_agriculture = .TRUE.   !! first call flag
!$OMP THREADPRIVATE(firstcall_sapiens_agriculture)
  INTEGER(i_std), SAVE       :: printlev_loc                            !! Local level of text output for current module
!$OMP THREADPRIVATE(printlev_loc)

CONTAINS

  !! ================================================================================================================================
  !! SUBROUTINE   : sapiens_agriculture_initialize
  !!
  !>\BRIEF        Initialize module sapiens_agriculture
  !!
  !! DESCRIPTION  : This subroutine is called from stomate_initialize. It will initialize printlev_loc variable by 
  !!               reading PRINTLEV_sapiens_agriculture from run.def.
  !!
  !! RECENT CHANGE(S) : 
  !!
  !! MAIN OUTPUT VARIABLE(S): 
  !!
  !! REFERENCE(S) :
  !!
  !! FLOWCHART    : None
  !! \n
  !_ =============================================================================================================================

  SUBROUTINE sapiens_agriculture_initialize()

    !! Initialize local printlev
      ! Guillaume M. -- get_printlev calls getin_p, a COLLECTIVE bcast. Called here,
      ! inside the compute loop under a local firstcall, it hangs the run as soon as the
      ! ranks stop reaching it together. Use printlev, read once at startup by every rank.
    printlev_loc=printlev
    firstcall_sapiens_agriculture=.FALSE.

  END SUBROUTINE sapiens_agriculture_initialize


  !! ================================================================================================================================
  !! SUBROUTINE   : crop_planting
  !!
  !>\BRIEF        Plant croplands
  !!
  !! DESCRIPTION  : Plant the crops when the weather is suitable. Uses begin_leaves to decide
  !! when to plant a crop. To be used within the DOFOCO branch.
  !!
  !! RECENT CHANGE(S) : 
  !!
  !! MAIN OUTPUT VARIABLE(S): 
  !!
  !! REFERENCE(S) :
  !!
  !! FLOWCHART    : None
  !! \n
  !_ =============================================================================================================================

  SUBROUTINE crop_planting(npts, dt, ipts, ivm, &
       veget_max, PFTpresent, c0_alloc, when_growthinit, &
       time_hum_min, everywhere, plant_status, &
       circ_class_n, KF, leaf_frac, age, &
       npp_longterm, lm_lastyearmax, circ_class_biomass, &
       atm_to_bm, k_latosa_adapt, soil_n_min, som, &
       deepSOM_a, deepSOM_s, deepSOM_p, zf_soil, & 
       n_estab_atm, n_estab_soil )


  !! 0. Variable and parameter declaration
  
  !! 0.1 Input variables
  
    INTEGER(i_std), INTENT(in)                         :: npts                 !! Domain size - number of grid cells (unitless)
    INTEGER(i_std), INTENT(in)                         :: ipts                 !! specific pixel to plant (unitless)
    INTEGER(i_std), INTENT(in)                         :: ivm                  !! specific pft to plant in the specific pixel
                                                                               !! (unitless)
    REAL(r_std), INTENT(in)                            :: dt                   !! time step (dt_days)  
    REAL(r_std), DIMENSION(:,:), INTENT(in)            :: veget_max            !! "maximal" coverage fraction of a 
                                                                               !! PFT (LAI -> infinity) on ground 
                                                                               !! (0-1, unitless)
    REAL(r_std), INTENT(in)                            :: c0_alloc             !! Root to sapwood tradeoff parameter
    REAL(r_std), DIMENSION(0:ngrnd), INTENT(in)        :: zf_soil

  !! 0.2 Ouput variables


  !! 0.3 Modified variables
    
    LOGICAL, DIMENSION(:,:), INTENT(inout)             :: PFTpresent           !! PFT exists (true/false)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: plant_status         !! Growth and phenological status of the plant
                                                                               !! Different stati are listed in constantes_var 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: when_growthinit      !! how many days since the 
                                                                               !! beginning of the growing season (days) 
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: time_hum_min         !! time elapsed since strongest 
                                                                               !! moisture availability (days)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: everywhere           !! is the PFT everywhere in the grid box or 
                                                                               !! very localized (after its introduction) (?)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: circ_class_n         !! Number of individuals in each circ class
                                                                               !! @tex $(ind m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: KF                   !! Scaling factor to convert sapwood mass
                                                                               !! into leaf mass (m)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: k_latosa_adapt       !! Leaf to sapwood area adapted for long 
                                                                               !! term water stress (m)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: leaf_frac            !! fraction of leaves in leaf age 
                                                                               !! class (0-1, unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: age                  !! mean age (years)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: npp_longterm         !! "long term" net primary productivity
                                                                               !! @tex ($gC m^{-2} year^{-1}$) @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)         :: lm_lastyearmax       !! last year's maximum leaf mass for each PFT 
                                                                               !! @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: atm_to_bm            !! co2 taken up by carbohydrate 
                                                                               !! reserve at the beginning of the 
                                                                               !! growing season @tex ($gC m^{-2} 
                                                                               !! of total ground/day$) @endtex
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)   :: circ_class_biomass   !! Biomass components of the model tree  
                                                                               !! within a circumference class
                                                                               !! class @tex $(g C ind^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)       :: soil_n_min           !!mineral nitrogen in the soil (gN/m**2)
 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: som                  !! Carbon pool: active, slow, or passive, 
                                                                               !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)    :: n_estab_atm          !! N taken from atm at establishment [g N m-2 day-1] 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)    :: n_estab_soil         !! N taken from soil at establishment [g N m-2 day-1] 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: deepSOM_a            !! Soil carbon discretized with depth active (g/m**3) 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: deepSOM_s            !! Soil carbon discretized with depth slow (g/m**3) 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)     :: deepSOM_p            !! Soil carbon discretized with depth passive (g/m**3) 

!! 0.4 Local variables
    
    REAL(r_std), DIMENSION(npts,nvm)                   :: LF                   !! Scaling factor to convert sapwood mass
                                                                               !! into root mass (unitless)
    REAL(r_std), DIMENSION(npts,nvm)                   :: histvar              !! controls the history output 
                                                                               !! level - 0: nothing is written; 
                                                                               !! 10: everything is written 
                                                                               !! (0-10, unitless)
    REAL(r_std), DIMENSION(npts,nvm)                   :: lstress_fac          !! Light stress factor, based on total
                                                                               !! transmitted light (unitless, 0-1)
    REAL(r_std)                                        :: Cs_grass             !! Individual plant, sapwood compartment
                                                                               !! @tex $(gC. ind^{-1})$ @endtex
    REAL(r_std)                                        :: Cl_init              !! Initial leaf carbon required to start
                                                                               !! the growing season
                                                                               !! @tex $(gC tree^{-1})$ @endtex
    REAL(r_std)                                        :: Cr_init              !! Initial root carbon required to start
                                                                               !! the growing season
                                                                               !! @tex $(gC tree^{-1})$ @endtex
                                                                               !! relationships (m) 
    INTEGER(i_std)                                     :: igrn                 !! index (unitless)
    REAL(r_std)                                        :: k_latosa_tmp         !! Temporaray variable in the calculation of 
                                                                               !! k_latosa_adapt                                                        
    REAL(r_std)                                        :: cn_leaf,cn_root,cn_wood !! C/N ratio of leaves, root and wood pool
                                                                               !! (gC/gN)
    REAL(r_std), DIMENSION(nelements)                  :: cn_tmp               !! Temporary variable used in the planting of
                                                                               !! saplings. C added to atm_to_bm. N taken from
                                                                               !! soil (g m**2)       
    REAL(r_std)                                        :: sla_est              !! The value of SLA using dyn_sla = y or n
    REAL(r_std), DIMENSION(nelements,ngrnd)            :: cn_tmp_discrete      !! Discretized temporary variable used in the planting of
                                                                               !! saplings. C added to atm_to_bm. N taken from
                                                                               !! soil or added atm_to_bm, dependeing on 
                                                                               !! ok_crop_planting_n_soil flag. (g m**2) 
    REAL(r_std), DIMENSION(npts,nvm,nelements)         :: tot_act              !! Total soil carbon and nitrogen for active
                                                                               !!  pool when soil carbon is discretized   
    REAL(r_std), DIMENSION(npts,nvm,nelements)         :: tot_slo              !! Total soil carbon and nitrogen for slow
                                                                               !!  pool when soil carbon is discretized 
    REAL(r_std), DIMENSION(npts,nvm,nelements)         :: tot_pas              !! Total soil carbon and nitrogen for passive
                                                                               !!  pool when soil carbon is discretized       
    REAL(r_std), DIMENSION(npts,ngrnd,nvm,nelements)   :: temp_deepSOM_a       !! Temporary variable for soil carbon discretized 
                                                                               !! with depth active (g/m**2) 
    REAL(r_std), DIMENSION(npts,ngrnd,nvm,nelements)   :: temp_deepSOM_s       !! Temporary variable for soil carbon discretized 
                                                                               !! with depth slow (g/m**2) 
    REAL(r_std), DIMENSION(npts,ngrnd,nvm,nelements)   :: temp_deepSOM_p       !! Temporary variable for soil carbon discretized 
                                                                               !! with depth passive (g/m**2)
!_ ================================================================================================================================


    IF (printlev_loc>=3) WRITE(numout,*) 'Entering sapiens_crop_planting'
        
    IF ( veget_max(ipts,ivm) .GT. min_stomate .AND. &
         SUM(circ_class_n(ipts,ivm,:)) .LT. min_stomate) THEN

       !! 1. Prescribe the crops
 
       ! Lightstress varies from 0 to 1 and is calculated from the canopy structure (veget)
       ! Given that there is no vegetation at this point, the lightstress cannot be 
       ! calculated and is set to 1. However as soon as we grow a canopy it will 
       ! experience light stress and so KF should be adjusted. If this is not done, we 
       ! can't prescribe tall vegetation which is a useful feature to speed up 
       ! optimisation and testing. We will have iterate over light stress and KF to 
       ! prescribe vegetation that is in balance with its light environment. 
       lstress_fac(ipts,ivm) = un

       ! initialize everything to make sure there are not random values floating
       ! around
       ! for ncirc = 1
       circ_class_biomass(ipts,ivm,:,:,:) = zero
       cn_tmp(:) = zero

       ! Similar as for trees, the initial height of the vegetation can be
       ! calculated from prescribed parameters. Use the function lai_to_biomass
       ! to account for dynamic sla. 
       ! Calculating leaf biomass here already to use in calculcation of KF/LF.
       ! Root and sapwood mass calculated below.
       circ_class_biomass(ipts,ivm,1,ileaf,icarbon) = lai_to_biomass(&
                height_init(ivm)/lai_to_height(ivm),ivm)
       
       ! For crops we assume that an individual crop is 1 m2 of crop. This is set
       ! in nmaxplants in pft_parameters.f90. It could be set to another value but 
       ! with the current code this should not have any meaning. The cropland
       ! does not necessarily covers the whole 1 m2 so adjust for the canopy cover
       circ_class_n(ipts,ivm,1) = nmaxplants(ivm) * ha_to_m2 * canopy_cover(ivm)

       ! Now we generate the size of a single crop sapling. We do not deal with 
       ! circumference classes for crops, but we want to keep the arrays the 
       ! same as for the trees so we put the sapling information into the first 
       ! circumference class

       ! Calculate the sapwood to leaf mass in a similar way as has been done for trees.
       ! For trees this approach had been justified by observations. For crops such
       ! justification is not supported by observations but we didn't try to find it.
       ! Needs more work by someone interested in crops. There might be a more elegant
       ! solution making use of a well observed parameter. 
       ! The mass of the structural carbon relates to the mass of the leaves through
       ! a prescribed parameter ::k_latosa
       k_latosa_tmp = k_latosa_adapt(ipts,ivm) + (lstress_fac(ipts,ivm) * &
            (k_latosa_max(ivm)-k_latosa_min(ivm)))

       IF (sla_dyn) THEN 
          sla_est = &
                biomass_to_lai(circ_class_biomass(ipts,ivm,1,ileaf,icarbon),ivm)/&
                circ_class_biomass(ipts,ivm,1,ileaf,icarbon)
       ELSE
          sla_est = sla(ivm)
       END IF
          
       KF(ipts,ivm) = k_latosa_tmp / & 
                (sla_est * pipe_density(ivm) * tree_ff(ivm))

       ! Calculate leaf to root area 
       LF(ipts,ivm) = c0_alloc * KF(ipts,ivm)

       ! Debug
       IF(printlev_loc>=4 .AND. ivm==test_pft)THEN
          WRITE(numout,*) 'sapiens_agriculture'
          WRITE(numout,*) 'KF, ', KF(ipts,ivm)
          WRITE(numout,*) 'LF, c0_alloc, ', c0_alloc * KF(ipts,ivm), &
               c0_alloc
       ENDIF
       !-

       ! initialize everything to make sure there are not random values floating around
       ! for ncirc = 1
       circ_class_biomass(ipts,ivm,:,:,:) = zero
       cn_tmp(:) = zero
       cn_tmp_discrete(:,:) =zero
 
       ! Similar as for trees, the initial height of the vegetation can be calculated 
       ! from prescribed parameters. Use the function lai_to_biomass to account
       ! for dynamic sla.
       IF (test_biomass_init) THEN
          circ_class_biomass(ipts,ivm,1,ileaf,icarbon) = biomass_init(ivm)
       ELSE
          circ_class_biomass(ipts,ivm,1,ileaf,icarbon) = lai_to_biomass(&
               height_init(ivm)/lai_to_height(ivm),ivm)
       ENDIF
       
       ! Use allometric relationships to define the root mass based on leaf mass. Some 
       ! sapwood mass is needed to store the reserves.
       circ_class_biomass(ipts,ivm,1,iroot,icarbon) = &
            circ_class_biomass(ipts,ivm,1,ileaf,icarbon) / LF(ipts,ivm)
       circ_class_biomass(ipts,ivm,1,isapabove,icarbon) = &
            circ_class_biomass(ipts,ivm,1,ileaf,icarbon) / KF(ipts,ivm)

       ! Pools that are defined in the same way for trees and grasses      
       circ_class_biomass(ipts,ivm,1,ilabile,icarbon) = labile_to_total * &
            ( circ_class_biomass(ipts,ivm,1,ileaf,icarbon) + &
            fcn_root(ivm) *  circ_class_biomass(ipts,ivm,1,iroot,icarbon) + &
            fcn_wood(ivm) * ( circ_class_biomass(ipts,ivm,1,isapabove,icarbon) +  &
            circ_class_biomass(ipts,ivm,1,isapbelow,icarbon) + &
            circ_class_biomass(ipts,ivm,1,icarbres,icarbon)))
       
       ! Allocate the required N for the phenological growth
       ! Similar as to the carbon this is taken from the atmosphere
       ! That is dealt with later in the code.
       ! Calculate the C/N ratio
       cn_leaf=cn_leaf_init(ivm)
       cn_wood=cn_leaf/fcn_wood(ivm)
       cn_root=cn_leaf/fcn_root(ivm) 
       
       ! Use the C/N ratio to calculate the N content for
       ! all the biomass components.
       circ_class_biomass(ipts,ivm,1,:,initrogen) = &
            circ_class_biomass(ipts,ivm,1,:,icarbon) / cn_root

       ! Overwrite the N content for the leaves and wood
       ! components.
       circ_class_biomass(ipts,ivm,1,ileaf,initrogen) = &
            circ_class_biomass(ipts,ivm,1,ileaf,icarbon) / cn_leaf
       circ_class_biomass(ipts,ivm,1,isapabove,initrogen) = &
            circ_class_biomass(ipts,ivm,1,isapabove,icarbon) / cn_wood
      
       ! Some of the biomass components that exist for trees are undefined for grasses
       circ_class_biomass(ipts,ivm,1,isapbelow,:) = zero
       circ_class_biomass(ipts,ivm,1,iheartabove,:) = zero
       circ_class_biomass(ipts,ivm,1,iheartbelow,:) = zero
       circ_class_biomass(ipts,ivm,1,ifruit,:) = zero
       circ_class_biomass(ipts,ivm,1,icarbres,:) = zero     
            
       ! Write initial values
       IF (printlev_loc>=4) THEN
          WRITE(numout,*) 'Crop prescribed biomass, carbon ', &
               circ_class_biomass(ipts,ivm,1,:,icarbon)
          WRITE(numout,*) 'Crop prescribed biomass, nitrogen ', &
               circ_class_biomass(ipts,ivm,1,:,initrogen)
       ENDIF
       
       ! Set leaf age classes -> all leaves will be current year leaves
       leaf_frac(ipts,ivm,:) = zero
       leaf_frac(ipts,ivm,1) = un
       
       ! Set time since last beginning of growing season but only
       ! for the first day of the whole simulation. When the model
       ! is initialized when_growthinit is set to undef. In subsequent
       ! time steps it should have a value.
       IF (when_growthinit(ipts,ivm) .EQ. undef) THEN
          when_growthinit(ipts,ivm) = 200
          time_hum_min(ipts,ivm) = 200
       ENDIF
      
       ! The biomass to build the saplings is taken from the atmosphere, 
       ! keep track of
       ! the amount to calculate the C-balance closure
       
        cn_tmp(:) = SUM(cc_to_biomass(npts, ivm, &
            circ_class_biomass(ipts,ivm,:,:,:), &
            circ_class_n(ipts,ivm,:)),1)

       atm_to_bm(ipts,ivm,icarbon) = &
            atm_to_bm(ipts,ivm,icarbon) + ( cn_tmp(icarbon) / dt )
     
       ! The nitrogen required is sourced from a soil pool (soil_n_min or som)
       ! If enough soil nitrogen available, 
       ! take nitrogen from either ammonium or
       ! nitrate pool in soil, otherwise take from som pool.
       IF (ok_soil_carbon_discretization) THEN
          ! NOTE: "deepSOM_x" does not exit in tag 2.1
          ! When soil is discreitzed we consider that the N taken in the som
          ! is only for the 5 top layers
          tot_act(ipts,ivm,initrogen) = zero
          tot_slo(ipts,ivm,initrogen) = zero
          tot_pas(ipts,ivm,initrogen) = zero
          temp_deepSOM_a(ipts,:,ivm,initrogen) = zero
          temp_deepSOM_s(ipts,:,ivm,initrogen) = zero
          temp_deepSOM_p(ipts,:,ivm,initrogen) = zero
          ! We convert deepSOM variables from g/m3 to g/m2 to be coherent with the 
          ! unit of cn_tmp
          DO igrn = 1,ngrnd
             temp_deepSOM_a(ipts,igrn,ivm,initrogen) = deepSOM_a(ipts,igrn,ivm,initrogen) * &
                  (zf_soil(igrn)-zf_soil(igrn-1))
             temp_deepSOM_s(ipts,igrn,ivm,initrogen) = deepSOM_s(ipts,igrn,ivm,initrogen) * &
                  (zf_soil(igrn)-zf_soil(igrn-1))
             temp_deepSOM_p(ipts,igrn,ivm,initrogen) = deepSOM_p(ipts,igrn,ivm,initrogen) * &
                  (zf_soil(igrn)-zf_soil(igrn-1))
          ENDDO
          
          DO igrn = 1,5
             tot_act(ipts,ivm,initrogen) = tot_act(ipts,ivm,initrogen) + &
                  deepSOM_a(ipts,igrn,ivm,initrogen)* (zf_soil(igrn)-zf_soil(igrn-1))
             tot_slo(ipts,ivm,initrogen) = tot_slo(ipts,ivm,initrogen) + &
                  deepSOM_s(ipts,igrn,ivm,initrogen)* (zf_soil(igrn)-zf_soil(igrn-1))
             tot_pas(ipts,ivm,initrogen) = tot_pas(ipts,ivm,initrogen) + &
                  deepSOM_p(ipts,igrn,ivm,initrogen)* (zf_soil(igrn)-zf_soil(igrn-1))
                cn_tmp_discrete(:,igrn) = cn_tmp(:) * (zf_soil(igrn)-zf_soil(igrn-1))/zf_soil(5)
          END DO

          IF ( cn_tmp(initrogen) .LE. &
               MAX(soil_n_min(ipts,ivm,initrate), &
               soil_n_min(ipts,ivm,iammonium)) ) THEN
             IF (soil_n_min(ipts,ivm,initrate) .GT. &
                  soil_n_min(ipts,ivm,iammonium)) THEN
                soil_n_min(ipts,ivm,initrate) = &
                     soil_n_min(ipts,ivm,initrate) - cn_tmp(initrogen)
             ELSE
                soil_n_min(ipts,ivm,iammonium) = &
                     soil_n_min(ipts,ivm,iammonium) - cn_tmp(initrogen)
             ENDIF

             n_estab_soil(ipts,ivm) = &
                n_estab_soil(ipts,ivm) + (cn_tmp(initrogen) / dt)
             
          ELSEIF ( cn_tmp(initrogen) .LT. tot_act(ipts,ivm,initrogen) ) THEN

              DO igrn = 1,5
                 IF ((temp_deepSOM_a(ipts,igrn,ivm,initrogen)- cn_tmp_discrete(initrogen,igrn)) .GE. zero) THEN
                    temp_deepSOM_a(ipts,igrn,ivm,initrogen)= &
                         temp_deepSOM_a(ipts,igrn,ivm,initrogen)- cn_tmp_discrete(initrogen,igrn)                         
                 ELSE
                      cn_tmp_discrete(initrogen,igrn+1) = cn_tmp_discrete(initrogen,igrn+1) + &
                                cn_tmp_discrete(initrogen,igrn) - temp_deepSOM_a(ipts,igrn,ivm,initrogen)
                      cn_tmp_discrete(initrogen,igrn) = temp_deepSOM_a(ipts,igrn,ivm,initrogen)
                      temp_deepSOM_a(ipts,igrn,ivm,initrogen)= zero
                   ENDIF
             ENDDO 

             n_estab_soil(ipts,ivm) = &
                 n_estab_soil(ipts,ivm) + (cn_tmp(initrogen) / dt)             

                IF (cn_tmp_discrete(initrogen,6) .GT. zero) THEN
                     ! Sometimes there is not enough nitrogen in either the
                     ! soil_n_min pools or n the first layers of the som pools. 
                     ! cn_temp_discrete at layer 5 is still higher than the som pools. 
                     ! Then we put some N in layer 6. In this case, we take the 
                     ! nitrogen from the atmosphere using atm_to_bm. We keep a record
                     ! of these instances using n_estab_atm.

                     ! Debug
                     IF(printlev_loc>=4 .AND. ivm==test_pft)THEN
                             WRITE(numout,*) 'ipts, ivm', ipts, ivm
                             WRITE(numout,*) 'cn_tmp(initrogen)', cn_tmp(initrogen)
                             WRITE(numout,*) ' tot_act(ipts,ivm,initrogen)', &
                                   tot_act(ipts,ivm,initrogen)
                             WRITE(numout,*) ' tot_pas(ipts,ivm,initrogen)', &
                                   tot_pas(ipts,ivm,initrogen)
                             WRITE(numout,*) ' tot_slo(ipts,ivm,initrogen)', &
                                   tot_slo(ipts,ivm,initrogen)
                     ENDIF

                     atm_to_bm(ipts,ivm,initrogen) = &
                          atm_to_bm(ipts,ivm,initrogen) + ( cn_tmp_discrete(initrogen,6) / dt )

                     n_estab_atm(ipts,ivm) = &
                          n_estab_atm(ipts,ivm) + (cn_tmp_discrete(initrogen,6) / dt)
               ENDIF

          ELSEIF ( cn_tmp(initrogen) .LT. tot_pas(ipts,ivm,initrogen)  ) THEN

              DO igrn = 1,5
                 IF ((temp_deepSOM_p(ipts,igrn,ivm,initrogen)- cn_tmp_discrete(initrogen,igrn)) .GE. zero) THEN
                    temp_deepSOM_p(ipts,igrn,ivm,initrogen)= &
                         temp_deepSOM_p(ipts,igrn,ivm,initrogen)- cn_tmp_discrete(initrogen,igrn)
                 ELSE
                      cn_tmp_discrete(initrogen,igrn+1) = cn_tmp_discrete(initrogen,igrn+1) + &
                                cn_tmp_discrete(initrogen,igrn) - temp_deepSOM_p(ipts,igrn,ivm,initrogen)
                      cn_tmp_discrete(initrogen,igrn) = temp_deepSOM_p(ipts,igrn,ivm,initrogen)
                      temp_deepSOM_p(ipts,igrn,ivm,initrogen)= zero
                   ENDIF
             ENDDO

             n_estab_soil(ipts,ivm) = &
                n_estab_soil(ipts,ivm) + (cn_tmp(initrogen) / dt)   

                IF (cn_tmp_discrete(initrogen,6) .GT. zero) THEN
                     ! Sometimes there is not enough nitrogen in either the
                     ! soil_n_min pools or n the first layers of the som pools. 
                     ! cn_temp_discrete at layer 5 is still higher than the som pools. 
                     ! Then we put some N in layer 6. In this case, we take the 
                     ! nitrogen from the atmosphere using atm_to_bm. We keep a record
                     ! of these instances using n_estab_atm.

                     ! Debug
                     IF(printlev_loc>=4 .AND. ivm==test_pft)THEN
                             WRITE(numout,*) 'ipts, ivm', ipts, ivm
                             WRITE(numout,*) 'cn_tmp(initrogen)', cn_tmp(initrogen)
                             WRITE(numout,*) ' tot_act(ipts,ivm,initrogen)', &
                                   tot_act(ipts,ivm,initrogen)
                             WRITE(numout,*) ' tot_pas(ipts,ivm,initrogen)', &
                                   tot_pas(ipts,ivm,initrogen)
                             WRITE(numout,*) ' tot_slo(ipts,ivm,initrogen)', &
                                   tot_slo(ipts,ivm,initrogen)
                     ENDIF

                     atm_to_bm(ipts,ivm,initrogen) = &
                          atm_to_bm(ipts,ivm,initrogen) + ( cn_tmp_discrete(initrogen,6) / dt )

                     n_estab_atm(ipts,ivm) = &
                          n_estab_atm(ipts,ivm) + (cn_tmp_discrete(initrogen,6) / dt)
               ENDIF

          ELSEIF ( cn_tmp(initrogen) .LE. tot_slo(ipts,ivm,initrogen)   ) THEN

              DO igrn = 1,5
                 IF ((temp_deepSOM_s(ipts,igrn,ivm,initrogen)- cn_tmp_discrete(initrogen,igrn)) .GE. zero) THEN
                    temp_deepSOM_s(ipts,igrn,ivm,initrogen)= &
                         temp_deepSOM_s(ipts,igrn,ivm,initrogen)- cn_tmp_discrete(initrogen,igrn)
                 ELSE
                      cn_tmp_discrete(initrogen,igrn+1) = cn_tmp_discrete(initrogen,igrn+1) + &
                                cn_tmp_discrete(initrogen,igrn) - temp_deepSOM_s(ipts,igrn,ivm,initrogen)
                      cn_tmp_discrete(initrogen,igrn) = temp_deepSOM_s(ipts,igrn,ivm,initrogen)
                      temp_deepSOM_s(ipts,igrn,ivm,initrogen)= zero
                   ENDIF
             ENDDO

             n_estab_soil(ipts,ivm) = &
                n_estab_soil(ipts,ivm) + (cn_tmp(initrogen) / dt)

                IF (cn_tmp_discrete(initrogen,6) .GT. zero) THEN
                     ! Sometimes there is not enough nitrogen in either the
                     ! soil_n_min pools or n the first layers of the som pools. 
                     ! cn_temp_discrete at layer 5 is still higher than the som pools. 
                     ! Then we put some N in layer 6. In this case, we take the 
                     ! nitrogen from the atmosphere using atm_to_bm. We keep a record
                     ! of these instances using n_estab_atm.

                     ! Debug
                     IF(printlev_loc>=4 .AND. ivm==test_pft)THEN
                             WRITE(numout,*) 'ipts, ivm', ipts, ivm
                             WRITE(numout,*) 'cn_tmp(initrogen)', cn_tmp(initrogen)
                             WRITE(numout,*) ' tot_act(ipts,ivm,initrogen)', &
                                   tot_act(ipts,ivm,initrogen)
                             WRITE(numout,*) ' tot_pas(ipts,ivm,initrogen)', &
                                   tot_pas(ipts,ivm,initrogen)
                             WRITE(numout,*) ' tot_slo(ipts,ivm,initrogen)', &
                                   tot_slo(ipts,ivm,initrogen)
                     ENDIF
     
                     atm_to_bm(ipts,ivm,initrogen) = &
                          atm_to_bm(ipts,ivm,initrogen) + ( cn_tmp_discrete(initrogen,6) / dt )

                     n_estab_atm(ipts,ivm) = &
                          n_estab_atm(ipts,ivm) + (cn_tmp_discrete(initrogen,6) / dt)
               ENDIF

            ELSE
                ! Sometimes there is not enough nitrogen in either the
                ! soil_n_min pools or the som pools. In this case, we take the 
                ! nitrogen from the atmosphere using atm_to_bm. We keep a record
                ! of these instances using n_estab_atm.

                ! Debug
                IF(printlev_loc>=4 .AND. ivm==test_pft)THEN
                        WRITE(numout,*) 'ipts, ivm', ipts, ivm
                        WRITE(numout,*) 'cn_tmp(initrogen)', cn_tmp(initrogen)
                        WRITE(numout,*) ' tot_act(ipts,ivm,initrogen)', &
                              tot_act(ipts,ivm,initrogen)
                        WRITE(numout,*) ' tot_pas(ipts,ivm,initrogen)', &
                              tot_pas(ipts,ivm,initrogen)
                        WRITE(numout,*) ' tot_slo(ipts,ivm,initrogen)', &
                              tot_slo(ipts,ivm,initrogen)
                ENDIF

                atm_to_bm(ipts,ivm,initrogen) = &
                     atm_to_bm(ipts,ivm,initrogen) + ( cn_tmp(initrogen) / dt )

                n_estab_atm(ipts,ivm) = &
                     n_estab_atm(ipts,ivm) + (cn_tmp(initrogen) / dt)

            ENDIF

             DO igrn = 1,ngrnd
                deepSOM_a(ipts,igrn,ivm,initrogen) = temp_deepSOM_a(ipts,igrn,ivm,initrogen) / &
                                                      (zf_soil(igrn)-zf_soil(igrn-1))
                deepSOM_s(ipts,igrn,ivm,initrogen) = temp_deepSOM_s(ipts,igrn,ivm,initrogen) / &
                                                      (zf_soil(igrn)-zf_soil(igrn-1))
                deepSOM_p(ipts,igrn,ivm,initrogen) = temp_deepSOM_p(ipts,igrn,ivm,initrogen) / &
                                                      (zf_soil(igrn)-zf_soil(igrn-1))
             ENDDO

             n_estab_soil(ipts,ivm) = &
                n_estab_soil(ipts,ivm) + (cn_tmp(initrogen) / dt) 

       ELSE ! not ok_soil_carbon_discretization
          IF ( cn_tmp(initrogen) .LE. &
               MAX(soil_n_min(ipts,ivm,initrate), &
               soil_n_min(ipts,ivm,iammonium)) ) THEN
             
             IF (soil_n_min(ipts,ivm,initrate) .GT. &
                  soil_n_min(ipts,ivm,iammonium)) THEN
                soil_n_min(ipts,ivm,initrate) = & 
                     soil_n_min(ipts,ivm,initrate) - cn_tmp(initrogen)
             ELSE
                soil_n_min(ipts,ivm,iammonium) = &
                     soil_n_min(ipts,ivm,iammonium) - cn_tmp(initrogen) 
             ENDIF

             n_estab_soil(ipts,ivm) = &
                n_estab_soil(ipts,ivm) + (cn_tmp(initrogen) / dt)
             
          ELSEIF ( cn_tmp(initrogen) .LT. som(ipts,iactive,ivm,initrogen) ) THEN
             som(ipts,iactive,ivm,initrogen) = &
                  som(ipts,iactive,ivm,initrogen) - cn_tmp(initrogen) 

             n_estab_soil(ipts,ivm) = &
                n_estab_soil(ipts,ivm) + (cn_tmp(initrogen) / dt)
             
          ELSEIF ( cn_tmp(initrogen) .LT. som(ipts,ipassive,ivm,initrogen) ) THEN
             som(ipts,ipassive,ivm,initrogen) = &
                  som(ipts,ipassive,ivm,initrogen) - cn_tmp(initrogen) 
             n_estab_soil(ipts,ivm) = &
                n_estab_soil(ipts,ivm) + (cn_tmp(initrogen) / dt)
             
          ELSEIF ( cn_tmp(initrogen) .LE. som(ipts,islow,ivm,initrogen)  ) THEN
             som(ipts,islow,ivm,initrogen) = &
                  som(ipts,islow,ivm,initrogen) - cn_tmp(initrogen)
             n_estab_soil(ipts,ivm) = &
                 n_estab_soil(ipts,ivm) + (cn_tmp(initrogen) / dt)
             
          ELSE  
             ! Sometimes there is not enough nitrogen in either the
             ! soil_n_min pools or the som pools. In this case, we take the 
             ! nitrogen from the atmosphere using atm_to_bm. We keep a record
             ! of these instances using n_estab_atm.
             
             ! Debug
             IF(printlev_loc>=4 .AND. ivm==test_pft)THEN
                WRITE(numout,*) 'ipts, ivm', ipts, ivm
                WRITE(numout,*) 'cn_tmp(initrogen)', cn_tmp(initrogen)
                WRITE(numout,*) 'som(ipts,iactive,ivm,initrogen)', &
                     som(ipts,iactive,ivm,initrogen)
                WRITE(numout,*) 'som(ipts,ipassive,ivm,initrogen)', &
                     som(ipts,ipassive,ivm,initrogen)
                WRITE(numout,*) 'som(ipts,islow,ivm,initrogen)', &
                     som(ipts,islow,ivm,initrogen)
             ENDIF
             
             atm_to_bm(ipts,ivm,initrogen) = &
                  atm_to_bm(ipts,ivm,initrogen) + ( cn_tmp(initrogen) / dt )
             
             n_estab_atm(ipts,ivm) = & 
                  n_estab_atm(ipts,ivm) + (cn_tmp(initrogen) / dt)
             
          ENDIF
       ENDIF
       

       cn_tmp(:) = zero

       !! 2.3 Declare PFT present 
       !! Now that the PFT has biomass it should be declared 'present'
       !! everywhere in that grid box. Assign some additional properties
       PFTpresent(ipts,ivm) = .TRUE.
       everywhere(ipts,ivm) = un
       age(ipts,ivm) = zero
       npp_longterm(ipts,ivm) = npp_longterm_init
       lm_lastyearmax(ipts,ivm) = zero

       !! 4.3 reset when_growthinit counter: start of the growing season
       when_growthinit(ipts,ivm) = zero
       plant_status(ipts,ivm) = icanopy
    
    ELSE

       ! PFT is not present, initialize the variables
       ! At the plant level
       circ_class_n(ipts,ivm,:) = zero
       circ_class_biomass(ipts,ivm,:,:,:) = zero
       
    ENDIF ! veget_max gt min_stomate

    IF (printlev_loc>=4) WRITE(numout,*) 'Leaving sapiens_crop_planting'

  END SUBROUTINE crop_planting


  !! ================================================================================================================================
  !! SUBROUTINE   : crop_harvest
  !!
  !>\BRIEF        Harvest of croplands
  !!
  !! DESCRIPTION  : To take into account biomass harvest from crop (mainly 
  !! to take into account for the reduced litter input and then decreased 
  !! soil carbon. it is a constant (40\%) fraction of turnover. To be used 
  !! within the DOFOCO branch.
  !!
  !! RECENT CHANGE(S) : Harvest is stored in harvest_pool for further 
  !! modelling of product use 
  !!
  !! MAIN OUTPUT VARIABLE(S): ::harvest_pool
  !!
  !! REFERENCE(S) :
  !! - Piao, S., P. Ciais, P. Friedlingstein, N. de Noblet-Ducoudre, P. Cadule, 
  !!   N. Viovy, and T. Wang. 2009. Spatiotemporal patterns of terrestrial 
  !!   carbon cycle during the 20th century. Global Biogeochemical Cycles 
  !!   23:doi:10.1029/2008GB003339.
  !!
  !! FLOWCHART    : None
  !! \n
  !_ ================================================================================================================================

  SUBROUTINE crop_harvest(npts, ipts, ivm, dt_days, &
       veget_max, turnover, harvest_pool, harvest_type, &
       harvest_cut, harvest_area, circ_class_biomass, &
       circ_class_n, leaf_meanage) 


    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER, INTENT(in)                                :: npts              !! Domain size (unitless)
    INTEGER(i_std), INTENT(in)                         :: ipts              !! specific pixel to plant (unitless)
    INTEGER(i_std), INTENT(in)                         :: ivm               !! specific pft to plant in the specific pixel
                                                                            !! (unitless)
    REAL(r_std), INTENT(in)                            :: dt_days           !! Time step (days)
    REAL(r_std), DIMENSION(:,:), INTENT(in)            :: veget_max         !! Fraction of the pixel covered by a PFT 
                                                                            !! @tex $(m^2 m^{-2})$ @endtex

    !! 0.2 Output variables
   
    !! 0.3 Modified variables

    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)  :: circ_class_biomass !! biomass @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(:,:,:,:),INTENT(inout)     :: turnover           !! Turnover rates 
                                                                            !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)  :: harvest_pool       !! The wood and biomass that have been
                                                                            !! havested by humans @tex $(gC)$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)        :: harvest_type       !! Type of management that resulted
                                                                            !! in the harvest (unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)        :: harvest_cut        !! Type of cutting that was used for the harvest
                                                                            !! (unitless)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)      :: harvest_area       !! Harvested area (m^{2})
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)      :: circ_class_n       !! Density of individuals 
                                                                            !! @tex $(m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)        :: leaf_meanage       !! mean age of the leaves (days)


    !! 0.4 Local variables
               

 !_ ================================================================================================================================

    IF(printlev_loc>=3) WRITE(numout,*) 'Entering sapiens_crop_harvest'

    !! 1. Initialize
    IF(printlev>=4 .AND. ipts == test_grid .AND. ivm == test_pft)THEN
       WRITE(numout,*) "Entering the crop_harvest routine: ", &
            SUM(SUM(SUM(harvest_pool(ipts,ivm,:,:,:),3),2))
    ENDIF

    ! We should not reset!  Sometimes there is harvest in the harvest 
    ! pool, and our mass balance will not be closed if we eliminate it 
    ! here. It seems possible that we hit senesence twice in a single 
    ! year. This happened on a grid point near the equator.
    ! harvest_pool(ipts,ivm,:,:,:) = zero


    !! 2. Harvest croplands

    ! NOTE that this is a quick fix rather than a final solution. The turnover
    ! of crop biomass has for long be calculated and was assumed to entirely 
    ! return to the litter pools. This approach ignored crop harvest. Rather 
    ! than harvesting the aboveground biomass (or below ground for some
    ! crops i.e. potatoes, beet, etc which is currently not accounted for) this
    ! routine makes use of the turnover and diverts part of the turnover into 
    ! the harvest pool. By doing so it does not end up in the litter.

    ! Assume all aboveground biomass is harvested and put it all in the first 
    ! diameter class.
    harvest_pool(ipts,ivm,1,:,iharvest) = harvest_pool(ipts,ivm,1,:,iharvest) + &
         harvest_ratio(ivm) * circ_class_n(ipts,ivm,1) * (circ_class_biomass(ipts,ivm,1,ileaf,:) + &
         circ_class_biomass(ipts,ivm,1,icarbres,:) + circ_class_biomass(ipts,ivm,1,ilabile,:) + &
         circ_class_biomass(ipts,ivm,1,ifruit,:) + circ_class_biomass(ipts,ivm,1,isapabove,:)) * &
         veget_max(ipts,ivm) * area(ipts) * contfrac(ipts)

    ! Associated harvest variables
    ! Multiplying with 10 allows to record all harvest types because the ifm_x value
    ! ranges between 1 and 7. When the model is run without disturbances harvest_cut
    ! also works correctly. If disturbances are uses icut_x is 10 or more and the simple
    ! trick to store all changes no longer works.
    harvest_type(ipts,ivm) = harvest_type(ipts,ivm) * 10 + ifm_crop
    harvest_cut(ipts,ivm) = harvest_cut(ipts,ivm) * 10 + icut_crop
    harvest_area(ipts,ivm,iharvest) = harvest_area(ipts,ivm,iharvest) + &
         veget_max(ipts,ivm) * area(ipts) * contfrac(ipts)

    ! The belowground becomes litter. All roots remain on site. 
    ! Half (::harvest_ratio) of the other biomas pools was harvested.
    turnover(ipts,ivm,iroot,:) = turnover(ipts,ivm,iroot,:) + &
         circ_class_biomass(ipts,ivm,1,iroot,:)*circ_class_n(ipts,ivm,1)
    turnover(ipts,ivm,ileaf,:) = turnover(ipts,ivm,ileaf,:) + &
         (un - harvest_ratio(ivm)) * circ_class_biomass(ipts,ivm,1,ileaf,:)*&
         circ_class_n(ipts,ivm,1)
    turnover(ipts,ivm,icarbres,:) = turnover(ipts,ivm,icarbres,:) + &
         (un - harvest_ratio(ivm)) * circ_class_biomass(ipts,ivm,1,icarbres,:)*&
         circ_class_n(ipts,ivm,1)
    turnover(ipts,ivm,ilabile,:) = turnover(ipts,ivm,ilabile,:) + &
         (un - harvest_ratio(ivm)) * circ_class_biomass(ipts,ivm,1,ilabile,:)*&
         circ_class_n(ipts,ivm,1)
    turnover(ipts,ivm,ifruit,:) = turnover(ipts,ivm,ifruit,:) + &
         (un - harvest_ratio(ivm)) * circ_class_biomass(ipts,ivm,1,ifruit,:)*&
         circ_class_n(ipts,ivm,1)
    turnover(ipts,ivm,isapabove,:) = turnover(ipts,ivm,isapabove,:) + &
         (un - harvest_ratio(ivm)) * circ_class_biomass(ipts,ivm,1,isapabove,:)*&
         circ_class_n(ipts,ivm,1)


    ! Account for the pools that should not contain biomass for the moment
    ! Do so to make sure the mass balance will be closed 
    turnover(ipts,ivm,isapbelow,:) = turnover(ipts,ivm,isapbelow,:) +  &
         circ_class_biomass(ipts,ivm,1,isapbelow,:)*&
         circ_class_n(ipts,ivm,1)
    turnover(ipts,ivm,iheartabove,:) = turnover(ipts,ivm,iheartabove,:) +  &
         circ_class_biomass(ipts,ivm,1,iheartabove,:)*&
         circ_class_n(ipts,ivm,1)
    turnover(ipts,ivm,iheartbelow,:) = turnover(ipts,ivm,iheartbelow,:) +  &
         circ_class_biomass(ipts,ivm,1,iheartbelow,:)*&
         circ_class_n(ipts,ivm,1)

    ! Reset the biomass pools
    circ_class_biomass(ipts,ivm,:,isapabove,:) = zero
    circ_class_biomass(ipts,ivm,:,ileaf,:) = zero
    circ_class_biomass(ipts,ivm,:,iroot,:) = zero
    circ_class_biomass(ipts,ivm,:,ifruit,:) = zero
    circ_class_biomass(ipts,ivm,:,icarbres,:) = zero
    circ_class_biomass(ipts,ivm,:,ilabile,:) = zero
                
    ! These pools should always be empty for crops - just to be sure
    circ_class_biomass(ipts,ivm,:,iheartbelow,:) = zero
    circ_class_biomass(ipts,ivm,:,iheartabove,:) = zero
    circ_class_biomass(ipts,ivm,:,isapbelow,:) = zero

    ! No individuals left
     circ_class_n(ipts,ivm,:) = zero

    ! reset leaf age
    leaf_meanage(ipts,ivm) = zero

    IF(printlev_loc>=4) WRITE(numout,*) 'Leaving sapiens_crop_harvest'

  END SUBROUTINE crop_harvest

END MODULE sapiens_agriculture
