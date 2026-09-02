! ================================================================================================================================
! MODULE       : stomate_litter
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       Update litter and lignine content after litter fall and 
!! calculating litter decomposition.      
!!
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S)	: None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_stomate/stomate_litter.f90 $
!! $Date: 2026-04-03 18:22:57 +0200 (ven. 03 avril 2026) $
!! $Revision: 9463 $
!! \n
!_ ================================================================================================================================

MODULE stomate_litter

  ! modules used:

  USE ioipsl_para
  USE stomate_data
  USE constantes
  USE constantes_soil
  USE pft_parameters
  USE function_library, ONLY : check_mass_balance, biomass_to_lai, cc_to_biomass, get_printlev
  USE grid, ONLY : area, contfrac

  IMPLICIT NONE

  ! private & public routines

  PRIVATE
  PUBLIC littercalc,littercalc_clear, deadleaf

  LOGICAL, SAVE             :: firstcall_litter = .TRUE.       !! first call
!$OMP THREADPRIVATE(firstcall_litter)
  INTEGER(i_std), SAVE       :: printlev_loc                   !! Local level of text output for current module
!$OMP THREADPRIVATE(printlev_loc)

CONTAINS

!! ================================================================================================================================
!! SUBROUTINE   : littercalc_clear
!!
!!\BRIEF        Set the flag ::firstcall_litter to .TRUE. and as such activate section
!! 1.1 of the subroutine littercalc (see below).
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : None
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE littercalc_clear
    firstcall_litter =.TRUE.
  END SUBROUTINE littercalc_clear


!! ================================================================================================================================
!! SUBROUTINE   : littercalc
!!
!!\BRIEF        Calculation of the litter decomposition and therefore of the 
!! heterotrophic respiration from litter following Parton et al. (1987).
!!
!! DESCRIPTION  : The littercal routine splits the litter in 4 pools: 
!! aboveground metaboblic, aboveground structural, belowground metabolic and 
!! belowground structural. the fraction (F) of plant material going to metabolic 
!! and structural is defined following Parton et al. (1987)
!! \latexonly
!! \input{littercalc1.tex}
!! \endlatexonly
!! \n
!! where L is the lignin content of the plant carbon pools considered and CN 
!! its CN ratio. L and CN are fixed parameters for each plant carbon pools,
!! therefore it is the ratio between each plant carbon pool within a PFT, which
!! controlled the part of the total litter, that will be considered as
!! recalcitrant (i.e. structural litter) or labile (i.e. metabolic litter),
!! woody or snag pool .\n  
!! 
!! The routine calculates the fraction of aboveground litter which is metabolic
!! structural woody or snag (the litterpart variable) which is then used in lpj_fire.f90.\n 
!! 
!! In the section 2, the routine calculate the new plant material entering the
!! litter pools by phenological death of plants organs (corresponding to the
!! variable turnover) and by fire, herbivory and others non phenological causes
!! (variable bm_to_litter). This calculation is first done for each PFT and then
!! the values calculated for each PFT are added up. Following the same approach
!! the lignin content of the total structural litter is calculated and will be
!! then used as a factor control of the decomposition of the structural, woody
!! and snag litter in the section 5.1.2. A test is performed to avoid that we add
!! more lignin than structural litter. Finally, the variable litterpart is
!! updated.\n
!! 
!! In the section 3 and 4 the temperature and the moisture controlling the
!! decomposition are calculated for above and belowground. For aboveground
!! litter, air temperature and litter moisture are calculated in sechiba and used 
!! directly. For belowground, soil temperature and moisture are also calculated 
!! in sechiba but are modulated as a function of the soil depth. The modulation 
!! is a multiplying factor exponentially distributed between 0 (in depth) and 1
!! in surface.\n  
!! 
!! Then, in the section 5, the routine calculates the structural litter decomposition 
!! (C) following first order kinetics following Parton et al. (1987).
!! \latexonly
!! \input{littercalc2.tex}
!! \endlatexonly
!! \n
!! with k the decomposition rate of the structural litter. 
!! k corresponds to
!! \latexonly
!! \input{littercalc3.tex}
!! \endlatexonly
!! \n
!! with littertau the turnover rate, T a function of the temperature and M a function of
!! the moisture described below.\n
!!  
!! Then, the fraction of dead leaves (DL) composed by aboveground structural litter is
!! calculated as following
!! \latexonly
!! \input{littercalc4.tex}
!! \endlatexonly
!! \n
!! with k the decomposition rate of the structural litter previously
!! described.\n
!! ?? below redundant with 5. ? 
!! In the section 5.1, the fraction of decomposed structural litter
!! incorporated to the soil (Input) and its associated heterotrophic respiration are
!! calculated. For structural litter, the C decomposed could go in the active
!! soil carbon pool or in the slow carbon, as described in 
!! stomate_soilcarbon.f90.\n
!! \latexonly
!! \input{littercalc5.tex}
!! \endlatexonly
!! \n
!! with f a parameter describing the fraction of structural litter incorporated
!! into the considered soil carbon pool, C the amount of litter decomposed and L 
!! the amount of lignin in the litter. The litter decomposed which is not
!! incorporated into the soil is respired.\n
!!
!! In the section 5.2, the fraction of decomposed metabolic litter
!! incorporated to the soil and its associated heterotrophic respiration are
!! calculated with the same approaches presented for 5.1 but no control factor
!! depending on the lignin content are used.\n
!! 
!! In the section 6 the dead leaf cover is calculated through a call to the 
!! deadleaf subroutine presented below.\n
!!
!! In the section 7, if the flag SPINUP_ANALYTIC is set to true, we fill MatrixA
!! and VectorB following Lardy(2011).
!!
!! MAIN OUTPUT VARIABLES: ::deadleaf_cover, ::resp_hetero_litter
!! ::control_temp, ::control_moist
!!
!! REFERENCES:
!! - Parton, WJ, Schimel, DS, Cole, CV, and Ojima, DS. 1987. Analysis
!! of factors controlling soil organic matter levels in Great Plains
!! grasslands. Soil Science Society of America journal (USA)
!! (51):1173-1179.
!! - Lardy, R, et al., A new method to determine soil organic carbon equilibrium,
!! Environmental Modelling & Software (2011), doi:10.1016|j.envsoft.2011.05.016
!!
!! FLOWCHART    :
!! \latexonly
!! \includegraphics(scale=0.5){littercalcflow.jpg}
!! \endlatexonly
!! \n
!_ ================================================================================================================================

  SUBROUTINE littercalc (npts, turnover, bm_to_litter, tree_bm_to_litter, &
       veget_max, tsurf, stempdiag, shumdiag, litterhum, som, clay, silt, & 
       soil_n_min, n_input, month, harvest_pool_acc, litter, dead_leaves, lignin_struc, &
       lignin_wood, lignin_snag, n_mineralisation, deadleaf_cover, resp_hetero_litter, &
       litterfuel, som_input, control_temp, control_moist, &
       MatrixA, VectorB, CN_target, CN_som_litter_longterm, tau_CN_longterm, &
       circ_class_biomass, circ_class_n, tsoil_decomp, snag_to_wood_flux, &
       shumdiag_peat, mc_peat_above, shumdiag_croppeat, mc_croppeat_above, &
       shumdiag_man, mc_man_above, soiltile)

    !! 0. Variable and parameter declaration
    
    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                           :: npts               !! Domain size - number of grid pixels
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)          :: turnover           !! Turnover rates of plant biomass 
                                                                               !! @tex $(gC m^{-2} dt\_slow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)          :: bm_to_litter       !! Conversion of biomass to litter 
                                                                               !! @tex $(gC m^{-2} dt\_slow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)          :: tree_bm_to_litter  !! Conversion of biomass to litter 
                                                                               !! @tex $(gC m^{-2} dt\_slow^{-1})$ @endtex 
    REAL(r_std),DIMENSION(:,:),INTENT(in)                :: veget_max          !! PFT "Maximal" coverage fraction of a PFT 
                                                                               !! defined in the input vegetation map 
                                                                               !! @tex $(m^2 m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(:), INTENT(in)                :: tsurf              !! Temperature (K) at the surface
    REAL(r_std), DIMENSION(:,:), INTENT(in)              :: stempdiag          !! Soil temperature (K)
    REAL(r_std), DIMENSION(:,:), INTENT(in)              :: shumdiag           !! Daily soil humidity of each soil layer 
                                                                               !! (unitless)
    REAL(r_std), DIMENSION(:), INTENT(in)                :: litterhum          !! Daily litter humidity (unitless)
    REAL(r_std), DIMENSION(:), INTENT(in)                :: clay               !! Clay fraction (unitless, 0-1)
    REAL(r_std), DIMENSION(:), INTENT(in)                :: silt               !! Silt fraction (unitless, 0-1)
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(in)        :: circ_class_biomass !! Biomass components of the model tree  
                                                                               !! within a circumference class 
                                                                               !! @tex $(g C ind^{=1})$ @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)            :: circ_class_n       !! Number of individuals in each circ class
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(inout)       :: n_input            !! nitrogen inputs into the soil  (gN/m**2/day)
    INTEGER(i_std), INTENT(in)                           :: month              !! month number required for n_input (1-12)

    !!!peatland
    REAL(r_std),DIMENSION (:,:), INTENT(in)              :: soiltile           !! soil tile fraction
    REAL(r_std),DIMENSION (npts,nslm),INTENT(in)         :: shumdiag_peat      !! soil moisture content for the decomposition when peat is activated
    REAL(r_std), DIMENSION(npts,nlevs)                   :: control_moist_peat !! Moisture control of heterotrophic 
                                                                               !! respiration (0.25-1, unitless), when peat is activated
    REAL(r_std), DIMENSION(npts),INTENT(in)              :: mc_peat_above      !! Liquid soil moisture content within the first 4 layers when peat is activated
    REAL(r_std), DIMENSION(npts),INTENT(in)              :: mc_croppeat_above  !! Liquid soil moisture content within the first 4 layers when agri_peat is activated
    REAL(r_std), DIMENSION(npts)                         :: mc_peat_below      !! soil humidity weighted by the assumed presence of the decomposers for peatland activated
    REAL(r_std),DIMENSION (npts,nslm),INTENT(in)         :: shumdiag_croppeat  !! soil moisture content for the decomposition when agri_peat is activated
    REAL(r_std), DIMENSION(npts,nlevs)                   :: control_moist_croppeat !! Moisture control of heterotrophic 
                                                                                   !! respiration (0.25-1, unitless), when agri_peat is activated
    REAL(r_std), DIMENSION(npts)                         :: mc_croppeat_below  !! soil humidity weighted by the assumed presence of the decomposers for agri_peat activated

    !!! tides
    REAL(r_std),DIMENSION (npts,nslm),INTENT(in)         :: shumdiag_man       !! soil moisture content for the decomposition when tides is activated
    REAL(r_std), DIMENSION(npts,nlevs)                   :: control_moist_man  !! Moisture control of heterotrophic 
                                                                               !! respiration (0.25-1, unitless), when agri_peat is activated
    REAL(r_std), DIMENSION(npts),INTENT(in)              :: mc_man_above       !! Liquid soil moisture content within the first 4 layers when tides is activated
    REAL(r_std), DIMENSION(npts)                         :: mc_man_below       !! soil temperature weighted by the assumed presence of the decomposers for tides being activated

    !! 0.2 Output variables
    
    REAL(r_std), DIMENSION(:), INTENT(out)               :: deadleaf_cover     !! Fraction of soil covered by dead leaves 
                                                                               !! over all PFTs (0-1, unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(out)             :: resp_hetero_litter !! Litter heterotrophic respiration. The unit
                                                                               !! is given by m^2 of ground.  
                                                                               !! @tex $(gC dt_sechiba one\_day^{-1}) m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(out)         :: som_input          !! Quantity of Carbon (or Nitrogen) going into SOM pools
                                                                               !! from litter decomposition. The unit is  
                                                                               !! given by m^2 of ground 
                                                                               !! @tex $(gC(orN) m^{-2} dt\_slow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(:,:), INTENT(out)             :: control_temp       !! Temperature control of heterotrophic 
                                                                               !! respiration, above and below (0-1, 
                                                                               !! unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(out)             :: control_moist      !! Moisture control of heterotrophic 
                                                                               !! respiration (0.25-1, unitless)
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(out)         :: MatrixA            !! Matrix containing the fluxes between the
                                                                               !! carbon pools per sechiba time step 
                                                                               !! @tex $(gC.m^2.day^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(out)           :: VectorB            !! Vector containing the litter increase per
                                                                               !! sechiba time step
                                                                               !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(out)           :: CN_target          !! C to N ratio of SOM flux from one pool to another (gN m-2 dt-1) 
    REAL(r_std),  DIMENSION(npts,nvm)                    :: ld_redistribute    !! logical set to redistribute som and litter  
    REAL(r_std), DIMENSION(:), INTENT(out)               :: tsoil_decomp       !! Temperature used for decompostition in soil (K)
    REAL(r_std), DIMENSION(npts,nvm,nlevs,nelements), INTENT(out) :: snag_to_wood_flux  !! Snag fall, litter goes from snag to woody pool (gC.m-2)


    !! 0.3 Modified variables
   
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(in)          :: som                !! SOM pools: active, slow, or passive, \f$(gC m^{2})$\f 
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)     :: litter             !! Metabolic and structural litter,above and
                                                                               !! below ground. The unit is given by m^2 of 
                                                                               !! ground @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)     :: litterfuel         !! Metabolic and structural litter,above and
                                                                               !! below ground. The unit is given by m2 of ground
                                                                               !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)         :: dead_leaves        !! Dead leaves per ground unit area, per PFT,
                                                                               !! metabolic and structural in 
                                                                               !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)         :: lignin_struc       !! Ratio Lignin content in structural litter,
                                                                               !! above and below ground, (0-1, unitless)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)         :: lignin_wood        !! Ratio Lignin/Carbon in woody litter,	
                                                                               !! above and below ground, (0-1, unitless)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)         :: lignin_snag        !! Ratio Lignin/Carbon in snag litter,	
                                                                               !! above and below ground, (0-1, unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)           :: n_mineralisation   !! Mineral N pool (gN m-2)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)         :: CN_som_litter_longterm !! Longterm CN ratio of litter and som pools (gC/gN)
    REAL(r_std), INTENT(inout)                           :: tau_CN_longterm    !! Counter used for calculating the longterm CN_ratio of som and litter pools [seconds]
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)     :: harvest_pool_acc   !! The wood and biomass havested by humans
    REAL(r_std), DIMENSION(:,:,:),INTENT(inout)          :: soil_n_min         !! mineral nitrogen in the soil (gN/m**2)

    !! 0.4 Local variables
 
    REAL(r_std)                                                 :: dt                 !! Number of day per sechiba(fast processes) time-step  
    REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:)          :: litterfrac         !! The fraction of leaves, wood, etc. that 
                                                                                      !! goes into metabolic and structural 
                                                                                      !! litterpools (0-1, unitless)
!$OMP THREADPRIVATE(litterfrac)
    REAL(r_std)                                                 :: rpc                !! Integration constant for vertical root 
                                                                                      !! profiles (unitless)
    REAL(r_std), SAVE, DIMENSION(nlitt)                         :: litter_dec_fac     !! Turnover time in litter pools (days)
!$OMP THREADPRIVATE(litter_dec_fac)
    REAL(r_std), SAVE                                           :: snag_fall_fac      !! Turnover time in snags from fall (days)
!$OMP THREADPRIVATE(snag_fall_fac)
    REAL(r_std), SAVE, DIMENSION(nlitt,ncarb,nlevs)             :: frac_soil          !! Fraction of litter that goes into soil 
                                                                                      !! (litter -> carbon, above and below). The
                                                                                      !! remaining part goes to the atmosphere
!$OMP THREADPRIVATE(frac_soil)
    REAL(r_std)                                                 :: soilhum_decomp     !! Humidity used for decompostition in soil
                                                                                      !! (unitless)
    REAL(r_std)                                                 :: fd                 !! Fraction of structural or metabolic litter
                                                                                      !! decomposed (unitless)
    REAL(r_std)                                                 :: qd_icarbon, qd_initrogen !! Quantity of structural or metabolic litter
                                                                                      !! decomposed @tex $(gC m^{-2})$ @endtex
                                                                                      !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std)                                                 :: old_struc          !! Old structural litter, above and below 
                                                                                      !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std)                                                 :: old_woody          !! Old woody litter, above and below 
                                                                                      !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std)                                                 :: old_snag           !! Old snag litter, above and below 
                                                                                      !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nlitt,nvm,nlevs,nelements)      :: litter_inc         !! Increase of metabolic and structural 
                                                                                      !! litter, above and below ground. The unit 
                                                                                      !! is given by m^2 of ground. 
                                                                                      !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nlitt,nvm,nhour,nelements)      :: litterfuel_inc     !! Increase of metabolic and structural
                                                                                      !! litter fuel
    REAL(r_std), DIMENSION(npts,nlitt,nvm,nelements)            :: litterfuel_all     !! The sum of litter fuel of all classes.

    REAL(r_std), DIMENSION(npts,nvm,nlevs)                      :: lignin_struc_inc   !! Lignin increase in structural litter, 
                                                                                      !! above and below ground. The unit is given 
                                                                                      !! by m^2 of ground. 
                                                                                      !! @tex $(gC m^{-2})$ @endtex                                             
    REAL(r_std), DIMENSION(npts,nvm,nlevs)                      :: lignin_wood_inc    !! Lignin increase in woody litter, 
                                                                                      !! above and below ground. The unit is given 
                                                                                      !! by m^2 of ground. 
                                                                                      !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm,nlevs)                      :: lignin_snag_inc    !! Lignin increase in woody litter, 
                                                                                      !! above and below ground. The unit is given 
                                                                                      !! by m^2 of ground. 
                                                                                      !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts)                                :: zdiff_min          !! Intermediate field for looking for minimum
                                                                                      !! of what? this is not used in the code. 
                                                                                      !! [??CHECK] could we delete it?
    CHARACTER(LEN=10), DIMENSION(nlitt)                         :: litter_str         !! Messages to write output information about
                                                                                      !! the litter
    CHARACTER(LEN=22), DIMENSION(nparts)                        :: part_str           !! Messages to write output information about
                                                                                      !! the plant
    CHARACTER(LEN=7), DIMENSION(ncarb)                          :: carbon_str         !! Messages to write output information about
                                                                                      !! the soil carbon
    CHARACTER(LEN=5), DIMENSION(nlevs)                          :: level_str          !! Messages to write output information about
                                                                                      !! the level (aboveground or belowground litter)
    INTEGER(i_std)                                              :: ilitt,ilev,ibdl    !! Indices (unitless)
    INTEGER(i_std)                                              :: ipts,ivm,j,ihour   !! Indices (unitless)
    INTEGER(i_std)                                              :: icarb,inspec,iele  !! Indices (unitless)
    INTEGER(i_std)                                              :: ipar,imbc,iicarb   !! Indices (unitless)
    REAL(r_std)                                                 :: f_soil_n_min       !! soil_n_min function response used for defing C to N target ratios (-) 
    REAL(r_std)                                                 :: f_pnc              !! plant nitrogen concentration function response
                                                                                      !! used for defining C to N target ratios (-)  
    INTEGER(i_std)                                              :: itarget            !! target som pool 
    REAL(r_std)                                                 :: CN                 !! CN ratio of the litter pools 
    REAL(r_std), DIMENSION(npts,nvm,nmbcomp,nelements)          :: check_intern       !! Contains the components of the internal
                                                                                      !! mass balance check for this routine
                                                                                      !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)                  :: closure_intern     !! Check closure of internal mass balance
                                                                                      !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)                  :: pool_start         !! Start and end pool of this routine 
                                                                                      !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)                  :: pool_end           !! Start and end pool of this routine 
                                                                                      !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)                            :: temp
    REAL(r_std), DIMENSION(npts)                                :: err                     !! Array for error checking 
    REAL(r_std), DIMENSION(ncarb,nelements)                     :: som_input_old           !! Used to store som_input values
                                                                                           !! in case of inadequate soil N
    REAL(r_std), DIMENSION(ncarb,nelements)                     :: delta_som_input         !! Used to store som_input values
                                                                                           !! in case of inadequate soil N
    REAL(r_std), DIMENSION(npts,nvm,nelements)                  :: som_input_new           !! Used to recalculate som_input in
                                                                                           !! case of inadequate soil N
    REAL(r_std)                                                 :: CN_input_total          !! overall C/N ratio of som_input_total
    REAL(r_std), DIMENSION(npts,nvm,nelements)                  :: som_input_total         !! Used to recalculate som_input in
                                                                                           !! case of inadequate soil N
    REAL(r_std), DIMENSION(npts,nvm,nelements)                  :: litter_total_new        !! Updated litter pool in case of nitrogen
                                                                                           !! deficit
    REAL(r_std), DIMENSION(npts,nvm)                            :: total_init_nitrogen     !! Litter pool at the very start of the 
                                                                                           !! calculations. Required to correct for the
                                                                                           !! nitrogen deficit
    REAL(r_std)                                                 :: n_mineralisation_init   !! temp variable in case of overspending
    REAL(r_std), DIMENSION(npts,nvm)                            :: delta_hetero_litter     !! reduction factor of resp_hetero
                                                                                           !! in case of nitrogen shortage
    REAL(r_std)                                                 :: share_litter_struc_a    !! Fractions of litter in various pools
    REAL(r_std)                                                 :: share_litter_struc_b    !! Fractions of litter in various pools
    REAL(r_std)                                                 :: share_litter_met_a      !! Fractions of litter in various pools
    REAL(r_std)                                                 :: share_litter_met_b      !! Fractions of litter in various pools
    REAL(r_std)                                                 :: share_litter_woody_a    !! Fractions of litter in various pools
    REAL(r_std)                                                 :: share_litter_woody_b    !! Fractions of litter in various pools
    REAL(r_std)                                                 :: share_litter_snag_a     !! Fractions of litter in various pools
    REAL(r_std)                                                 :: share_litter_snag_b     !! Fractions of litter in various pools
    REAL(r_std)                                                 :: share_som_input_active  !! Active pool fraction of 
    REAL(r_std)                                                 :: share_som_input_slow    !! Slow pool fraction of
    REAL(r_std)                                                 :: share_som_input_surface !! Surface pool fraction of
    REAL(r_std), DIMENSION(npts,nvm)                            :: n_min_old               !! Temp variable
    INTEGER(i_std), DIMENSION(2)                                :: ix                      !! Row/column number of the maximal value
    INTEGER(i_std)                                              :: ier                     !! Error handling
    REAL(r_std), DIMENSION(npts)                                :: error_count             !! Count the number of errors
    REAL(r_std), DIMENSION(nlitt)                               :: tree_litterfrac
    REAL(r_std)                                                 :: decomp_residual         !! temporary variable to ensure mass balance closure for som_input and hetero_resp
    REAL(r_std)                                                 :: temp_flux               !! temporary variable for litter to som fluxes
    REAL(r_std)                                                 :: total                   !! Total of share_litter used as temporary variable (unitless,0-1)
    REAL(r_std)                                                 :: err_min_1,err_min_2,err_min_3,err_min_4,err_min_5,err_min_6,err_min_7,err_min_8,err_max_1
    REAL(r_std)                                                 :: err_sum_1,err_sum_2,err_sum_3,err_sum_4,err_sum_5,err_sum_6,err_sum_7,err_sum_8
    INTEGER(i_std)                                              :: istm
!_ ================================================================================================================================
    
    IF ( firstcall_litter ) THEN
      ! 1.1 Initialize local printlev
      ! Guillaume M. -- get_printlev calls getin_p, a COLLECTIVE bcast. Called here,
      ! inside the compute loop under a local firstcall, it hangs the run as soon as the
      ! ranks stop reaching it together. Use printlev, read once at startup by every rank.
      printlev_loc=printlev
    END IF

    IF (printlev_loc>=3) WRITE(numout,*) 'Entering littercalc'

  !! 1. Initialisations of the different fields during the first call of the routine 
    dt = dt_sechiba/one_day

    IF ( firstcall_litter ) THEN

       IF ( .NOT. ALLOCATED(litterfrac) ) THEN
          ALLOCATE(litterfrac(npts,nvm,nparts,nlitt))
       ENDIF

       !! 1.1.2 exponential decay rates conversion (day-1)
       !  The corresponding 95% residence times in litter pools 
       !  are respectively of .2 years, 0.75 years, 2.25 and 225 years
       litter_dec_fac(imetabolic) = turn_metabolic / one_year      
       litter_dec_fac(istructural) = turn_struct / one_year        
       litter_dec_fac(iwoody) = turn_woody / one_year  
       litter_dec_fac(isnag) = turn_snag / one_year  
      
       !! 1.1.3 exponential decay rate for snag pool from fall (day-1)
       !  95% residence time of snags regarding fall is 22.5 years 
       snag_fall_fac = turn_snag_fall / one_year      
       
       !! 1.1.5 decomposition flux fraction that goes into soil
       !  (litter -> carbon, above and below)
       !  1-frac_soil goes into atmosphere
       frac_soil(:,:,:) = zero

       ! structural litter: lignin fraction goes into slow pool + respiration,
       ! rest into active and surface pool + respiration
       frac_soil(istructural,isurface,iabove) = frac_soil_struct_sua
       frac_soil(istructural,iactive,iabove) = frac_soil_struct_aa
       frac_soil(istructural,iactive,ibelow) = frac_soil_struct_ab
       frac_soil(istructural,islow,iabove) = frac_soil_struct_sa
       frac_soil(istructural,islow,ibelow) = frac_soil_struct_sb

       ! metabolic litter: all goes into active pool + respiration.
       ! Nothing into slow or passive pool.
       frac_soil(imetabolic,isurface,iabove) = frac_soil_metab_sua
       frac_soil(imetabolic,iactive,iabove) = frac_soil_metab_aa
       frac_soil(imetabolic,iactive,ibelow) = frac_soil_metab_ab

       ! woody litter: lignin fraction goes into slow pool + respiration,
       ! rest into active and surface pool + respiration
       ! Remark : woody litter lignin fraction should be lower than that of structural
       ! litter, else woody litter will not decompose quickly enough.
       ! In effect, this decreases the carbon use efficiency of the system to  
       ! allow for faster litter decomposition. The parameter frac_woody is not
       ! well-constrained and was tested for a boreal PFT in northern Sweden.       
       frac_soil(iwoody,isurface,iabove) = frac_woody * frac_soil_struct_sua
       frac_soil(iwoody,iactive,iabove) = frac_woody * frac_soil_struct_aa
       frac_soil(iwoody,iactive,ibelow) = frac_woody * frac_soil_struct_ab
       frac_soil(iwoody,islow,iabove) = frac_woody * frac_soil_struct_sa
       frac_soil(iwoody,islow,ibelow) = frac_woody * frac_soil_struct_sb
       
       ! snag litter: lignin fraction goes into slow pool + respiration,
       ! rest into active and surface pool + respiration
       frac_soil(isnag,isurface,iabove) = frac_snag * frac_soil_struct_sua
       frac_soil(isnag,iactive,iabove) = frac_snag * frac_soil_struct_aa
       frac_soil(isnag,iactive,ibelow) = frac_snag * frac_soil_struct_ab
       frac_soil(isnag,islow,iabove) = frac_snag * frac_soil_struct_sa
       frac_soil(isnag,islow,ibelow) = frac_snag * frac_soil_struct_sb

       !! 1.3 messages
       litter_str(imetabolic) = 'metabolic'
       litter_str(istructural) = 'structural'
       litter_str(iwoody) = 'woody'
       litter_str(isnag) = 'snag'

       carbon_str(iactive) = 'active'
       carbon_str(isurface) = 'surface'
       carbon_str(islow) = 'slow'
       carbon_str(ipassive) = 'passive'

       level_str(iabove) = 'above'
       level_str(ibelow) = 'below'

       part_str(ileaf) = 'leaves'
       part_str(isapabove) = 'sap above ground'
       part_str(isapbelow) = 'sap below ground'
       part_str(iheartabove) = 'heartwood above ground'
       part_str(iheartbelow) = 'heartwood below ground'
       part_str(iroot) = 'roots'
       part_str(ifruit) = 'fruits'
       part_str(icarbres) = 'carbohydrate reserve'
       part_str(ilabile) = 'labile reserve'

       ! Debug
       IF (printlev_loc >= 4) THEN
          WRITE(numout,*) 'litter:'

          WRITE(numout,*) '   > C/N ratios: '
          DO ipar = 1, nparts
             WRITE(numout,*) '       ', part_str(ipar), ': ',CN_fix(ipar)
          ENDDO
          
          WRITE(numout,*) '   > scaling depth for decomposition (m): ',z_decomp

          WRITE(numout,*) '   > decay rates of litter pools (d-1):'
          DO ilitt = 1, nlitt
             WRITE(numout,*) '       ',litter_str(ilitt),':',litter_dec_fac(ilitt)
          ENDDO
          
          WRITE(numout,*) '   > litter decomposition flux fraction that really goes '
          WRITE(numout,*) '     into carbon pools (rest into the atmosphere):'
          DO ilitt = 1, nlitt
             DO ilev = 1, nlevs
                DO icarb = 1, ncarb
                   WRITE(numout,*) '       ',litter_str(ilitt),' ',level_str(ilev),' -> ',&
                        carbon_str(icarb),':', frac_soil(ilitt,icarb,ilev)
                ENDDO
             ENDDO
          ENDDO

       ENDIF
       !-
          
       firstcall_litter = .FALSE.

    ENDIF ! firstcall_litter

    !! 1.1.3 litter fractions:
    !  What fraction of leaves, wood, etc. goes into metabolic and 
    !  structural litterpools
    ! Initialised for when everything is printed below
    !+++CHECK+++
    ! Add documentation
    ! Calculate the CN ratio. If the CN ratio cannot be calculated
    ! then use the prescribed value
    DO ipar = 1, nparts
       ! Note that all the carbon and nitrogen contained in tree_bm_to_litter
       ! is also contained in bm_to_litter. tree_bm_to_litter is only needed
       ! to ensure that some of the woody litter that moved to bare soil, crop/grass
       ! PFTs during LCC is correctly dealt with. tree_bm_to_litter should
       ! NOT be accounted for in most of the calculations.
       DO ivm = 1,nvm
          DO ipts = 1,npts
             litterfrac(ipts,ivm,ipar,:) = zero 
             CN = zero
             IF((bm_to_litter(ipts,ivm,ipar,initrogen)+turnover(ipts,ivm,ipar,initrogen)).GT.min_stomate) THEN
                     CN = (bm_to_litter(ipts,ivm,ipar,icarbon)+turnover(ipts,ivm,ipar,icarbon)) / &
                             (bm_to_litter(ipts,ivm,ipar,initrogen)+turnover(ipts,ivm,ipar,initrogen))
             ELSE
                     CN = CN_fix(ipar)  
             ENDIF
             IF ( ((ipar == isapabove) .OR. (ipar == isapbelow) .OR. &
                     (ipar == iheartabove) .OR. (ipar == iheartbelow)) .AND. &
                     ivm == ibare_sechiba) THEN
                     IF (bm_to_litter(ipts,ivm,ipar,icarbon).GT.10*EPSILON(un))THEN
                        ! PFT 1 with woody litter will be treated as
                        ! a forest. The threshold needs to be much lower than min_stomate.
                        ! If not, mass is being lost. Because this calculation is repeated
                        ! 48 times, small errors (1e-9) soon exceed 1e-8 and make the model
                        ! crash in the consistency cross-check for nbp.
                        ! First we hard-code a fifty-fifty division of woody
                        ! pool in woody and snag pools
                        IF ( ok_snags ) THEN
                           litterfrac(ipts,ivm,ipar,iwoody) = 1 - snag_frac 
                           litterfrac(ipts,ivm,ipar,isnag) = snag_frac
                        ELSE 
                           litterfrac(ipts,ivm,ipar,iwoody) = 1.0
                           litterfrac(ipts,ivm,ipar,isnag) = zero
                        ENDIF
                        litterfrac(ipts,ivm,ipar,imetabolic) = zero
                        litterfrac(ipts,ivm,ipar,istructural) = zero
                     ENDIF
             ELSEIF ( ((ipar == isapabove) .OR. (ipar == isapbelow) .OR. &
                             (ipar == iheartabove) .OR. (ipar == iheartbelow)) &
                             .AND. is_tree(ivm) ) THEN
                ! Woody components of forest PFTs
                IF ( ok_snags ) THEN
                   litterfrac(ipts,ivm,ipar,iwoody) = 1 - snag_frac
                   litterfrac(ipts,ivm,ipar,isnag) = snag_frac
                ELSE 
                   litterfrac(ipts,ivm,ipar,iwoody) = 1.0
                   litterfrac(ipts,ivm,ipar,isnag) = zero
                ENDIF
                litterfrac(ipts,ivm,ipar,imetabolic) = zero
                litterfrac(ipts,ivm,ipar,istructural) = zero
             ELSE 
                ! PFT 1 without woody litter and all non-forest PFTs
                IF( (ipar == icarbres) .OR. (ipar == ilabile)) THEN
                   litterfrac(ipts,ivm,ipar,imetabolic) = 1.0
                ELSE
                ! This repartition is derived from eq 1 of Parton et al. 1993
                   litterfrac(ipts,ivm,ipar,imetabolic) = &
                           MAX(metabolic_ref_frac - metabolic_LN_ratio  * &
                           LC(ivm,ipar) * CN,zero)
                ENDIF
                litterfrac(ipts,ivm,ipar,istructural) = &
                        1. - litterfrac(ipts,ivm,ipar,imetabolic)
                litterfrac(ipts,ivm,ipar,iwoody) = zero
                litterfrac(ipts,ivm,ipar,isnag) = zero
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    ! Error checking
    err_max_1 = zero
    IF(err_act.GT.1)THEN
       DO ipar = 1, nparts 
          DO ivm = 1,nvm
             DO ipts = 1,npts
                IF((bm_to_litter(ipts,ivm,ipar,initrogen) + &
                     turnover(ipts,ivm,ipar,initrogen)).GT.min_stomate)THEN
                   err_max_1 = MAX(err_max_1,ABS(SUM(litterfrac(ipts,ivm,ipar,:))-1.0))
                ENDIF
             ENDDO
          ENDDO
       ENDDO

       IF (err_max_1.GT.10*EPSILON(un)) THEN
          CALL ipslerr_p(3,'stomate_litter',&
                  'litterfrac do not add up to 1','','')
       ENDIF
    ENDIF
    
    ! Error checking
    ! soil_n_min should never be zero. It is a nitrogen pool
    ! expressed in gN m-2
    ! Error checking
    ! This code is OK only when the land cover change and the age
    ! class distribution code are working properly. Better to test
    ! it in those subroutines than to enforce it here. If there are
    ! problems with land cover change or age class distribution the
    ! mass balance for soil_n_min may be violated here.
    err_sum_1 = zero
    err_sum_2 = zero
    err_sum_3 = zero

    DO inspec = 1, nnspec
       DO ivm = 1, nvm
          DO ipts = 1, npts
             IF((veget_max(ipts,ivm).LT.min_stomate) .AND. &
                  (soil_n_min(ipts,ivm,inspec).GT. min_stomate)) THEN
                err_sum_3 = err_sum_3 + un  
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    DO ivm = 1, nvm      
       DO ipts = 1, npts 
          IF(soil_n_min(ipts,ivm,initrate).LT.zero) THEN
                  err_sum_1 = err_sum_1 + un
          ENDIF
          IF(soil_n_min(ipts,ivm,iammonium).LT.zero)THEN
                  err_sum_2 = err_sum_2 + un
          ENDIF
       ENDDO
    ENDDO

    IF (err_sum_1.GT.zero) THEN
       CALL ipslerr_p(3,'stomate_litter',&
            'At least one soil_n_min value initrate',&
            'is negative in stomate_litter','')
    ENDIF
    IF (err_sum_2.GT.zero) THEN
       CALL ipslerr_p(3,'stomate_litter iammonium',&
            'At least one soil_n_min value',&
            'is negative in stomate_litter','')
    ENDIF
    IF(err_sum_3.GT.zero) THEN
            CALL ipslerr_p(3,'stomate_litter',&
                    'soil_n_min present in pixels without veget_max',&
                    '','')
    ENDIF
    ! f_soil_n_min is derived from Fig.4 of Parton et al. 1993, which shows CN
    ! ratio as an affine function of the sum of soil ammonium and nitrate which
    ! becomes constant above 2.0
    ! It should not be < 0 but the check is already done above for ammonium and
    ! nitrate from which it is calculated

    DO ivm = 1,nvm
       DO ipts = 1,npts
          f_soil_n_min = &
                  MIN((soil_n_min(ipts,ivm,iammonium) + soil_n_min(ipts,ivm,initrate)), 2.)
          ! C to N target ratios of different pools
          CN_target(ipts,ivm,:)=zero
          ! CN ratio ranges between 3 and 15 for active pool
          ! Figure 4 of Parton et al. (1993) show lower CN values for active 
          ! pool of 2 rather than 3
          CN_target(ipts,ivm,iactive)= CN_target_iactive_ref + &
                  CN_target_iactive_Nmin * f_soil_n_min

          ! CN ratio ranges between 12 and 20 for slow pool
          CN_target(ipts,ivm,islow)= CN_target_islow_ref + &
                  CN_target_islow_Nmin * f_soil_n_min
          ! CN ratio ranges between 7 and 10 for passive pool
          ! OCN uses a fixed value of 9. (don't know why)
          ! Figure 4a of Parton et al. (1993) show lower CN values 
          ! for passive pool of 3 rather than 7
          CN_target(ipts,ivm,ipassive)= CN_target_ipassive_ref + &
                  CN_target_ipassive_Nmin * f_soil_n_min
          ! Litter nitrogen content (%)
          IF(litter(ipts,imetabolic,ivm,iabove,icarbon)+&
                  litter(ipts,istructural,ivm,iabove,icarbon)+&
                  litter(ipts,iwoody,ivm,iabove,icarbon) .GT. zero) THEN


             !+++CHECK+++
             ! Checking the fig. 4b of Parton 1993 it seems that f_pnc is 
             ! an estimation of the N content (%) of the plant. However it would
             ! need to be documented why this estimation relies on the N to 2*C
             ! ratio of litter. Even in Parton originally this 2.0 seems somehow
             ! arbitrary as it is found in both figures (4a and b) though the units are not
             ! the same (g.m-2 and %). Snags are not added here as they are not
             ! mixed with the rest of the litter pools and thus their CN should
             ! nor influence the CN of the ground litter 
             f_pnc=MIN(( litter(ipts,imetabolic,ivm,iabove,initrogen)+&
                     litter(ipts,istructural,ivm,iabove,initrogen)+&
                     litter(ipts,iwoody,ivm,iabove,initrogen) ) / &
                     (2.*(litter(ipts,imetabolic,ivm,iabove,icarbon)+&
                     litter(ipts,istructural,ivm,iabove,icarbon)+&
                     litter(ipts,iwoody,ivm,iabove,icarbon) )) * 100., 2.) 
             !+++++++++++
          ELSE
             f_pnc= zero
          ENDIF    
          ! 
          CN_target(ipts,ivm,isurface) = CN_target_isurface_ref + &
                  CN_target_isurface_pnc * f_pnc
          !+++++++++++
          !! 1.4 set output to zero
          resp_hetero_litter(ipts,ivm) = zero
          som_input(ipts,:,ivm,:) = zero
       ENDDO
    ENDDO

    
    !! 1.5 Initialize check for mass balance closure
    !  The mass balance is calculated at the end of this routine
    !  in section 6. Initial biomass and harvest pool all other
    !  relevant pools were just set to zero.
    IF (err_act.GT.1) THEN
       DO ivm = 1,nvm
          DO ipts = 1,npts
             pool_start(ipts,ivm,:) = zero
             pool_start(ipts,ivm,initrogen) = pool_start(ipts,ivm,initrogen) + &
                        n_mineralisation(ipts,ivm) * veget_max(ipts,ivm)
             pool_start(ipts,ivm,initrogen) = pool_start(ipts,ivm,initrogen) + &
                        (soil_n_min(ipts,ivm,iammonium) + soil_n_min(ipts,ivm,initrate))*veget_max(ipts,ivm)
                   ! The litter pool
             DO iele = 1,nelements
                DO ilitt = 1,nlitt
                   DO ilev = 1,nlevs
                      pool_start(ipts,ivm,iele) = pool_start(ipts,ivm,iele) + &
                              (litter(ipts,ilitt,ivm,ilev,iele) * veget_max(ipts,ivm))
                   ENDDO
                ENDDO
                ! Pools added to the litter
                ! Note that all the carbon and nitrogen conatined in tree_bm_to_litter
                ! is also contained in bm_to_litter. tree_bm_to_litter is only needed
                ! to ensure that some of the woody litter that moved to a crop/grass
                ! PFT during LCC is correctly dealt with. tree_bm_to_litter should
                ! NOT be accounted for in most of the calculations.
                DO ipar = 1,nparts
                   pool_start(ipts,ivm,iele) = pool_start(ipts,ivm,iele) + &
                           (turnover(ipts,ivm,ipar,iele) + &
                           bm_to_litter(ipts,ivm,ipar,iele)) * veget_max(ipts,ivm)
                ENDDO
                pool_start(ipts,ivm,iele) = pool_start(ipts,ivm,iele) + &
                        SUM(SUM(harvest_pool_acc(ipts,ivm,:,iele,:),2))/(area(ipts)*contfrac(ipts))
             ENDDO
          ENDDO
       ENDDO
    ENDIF ! err_act.GT.1
    ! Store the initial available N. In case there is overspending of nitrogen
    ! this value can be used to adjust the spending

    DO ivm = 1,nvm
       DO ipts = 1,npts
          IF (veget_max(ipts,ivm).GT.min_stomate)THEN
             total_init_nitrogen(ipts,ivm) = n_mineralisation(ipts,ivm) +  &
                     SUM(litter(ipts,:,ivm,:,initrogen)) + SUM(turnover(ipts,ivm,:,initrogen)) + &
                     SUM(bm_to_litter(ipts,ivm,:,initrogen))
          ELSE
             total_init_nitrogen(ipts,ivm) = zero
          ENDIF
       ENDDO
    ENDDO

    DO ilev = 1, nlevs !Loop over litter levels (above and below ground)
       DO ivm = 1,nvm !Loop over PFTs
          DO ipts = 1,npts
             lignin_struc_inc(ipts,ivm,ilev) = zero
             lignin_wood_inc(ipts,ivm,ilev) = zero
             lignin_snag_inc(ipts,ivm,ilev) = zero
             litter_inc(ipts,:,ivm,ilev,:) = zero
          ENDDO
       ENDDO
    ENDDO
    
    litterfuel_inc(:,:,:,:,:) = zero
    IF ( ok_snags ) THEN
       tree_litterfrac(iwoody) = 1 - snag_frac
       tree_litterfrac(isnag) = snag_frac
    ELSE
       tree_litterfrac(iwoody) = 1.0
       tree_litterfrac(isnag) = zero
    ENDIF
    tree_litterfrac(imetabolic) = zero
    tree_litterfrac(istructural) = zero

    DO ivm = 1,nvm !Loop over PFTs
       DO ipts = 1,npts
          !! 2.2.1 litter
          DO ilitt = 1, nlitt        !Loop over litter pools (metabolic and structural)
             DO iele = 1, nelements  ! Loop over element pools (carbon and nitrogen)
                litter_inc(ipts,ilitt,ivm,iabove,iele) = & 
                        litterfrac(ipts,ivm,ileaf,ilitt) * bm_to_litter(ipts,ivm,ileaf,iele) + & 
                        litterfrac(ipts,ivm,isapabove,ilitt) * &
                        (bm_to_litter(ipts,ivm,isapabove,iele)-tree_bm_to_litter(ipts,ivm,isapabove,iele)) + &  
                        litterfrac(ipts,ivm,iheartabove,ilitt) * &
                        (bm_to_litter(ipts,ivm,iheartabove,iele)-tree_bm_to_litter(ipts,ivm,iheartabove,iele)) + &  
                        tree_litterfrac(ilitt) * tree_bm_to_litter(ipts,ivm,isapabove,iele) + &  
                        tree_litterfrac(ilitt) * tree_bm_to_litter(ipts,ivm,iheartabove,iele) + & 
                        litterfrac(ipts,ivm,ifruit,ilitt) * bm_to_litter(ipts,ivm,ifruit,iele) + & 
                        litterfrac(ipts,ivm,icarbres,ilitt) * bm_to_litter(ipts,ivm,icarbres,iele) + & 
                        litterfrac(ipts,ivm,ilabile,ilitt) * bm_to_litter(ipts,ivm,ilabile,iele) + & 
                        litterfrac(ipts,ivm,ileaf,ilitt) * turnover(ipts,ivm,ileaf,iele) + & 
                        litterfrac(ipts,ivm,isapabove,ilitt) * turnover(ipts,ivm,isapabove,iele) + & 
                        litterfrac(ipts,ivm,iheartabove,ilitt) * turnover(ipts,ivm,iheartabove,iele) + & 
                        litterfrac(ipts,ivm,ifruit,ilitt) * turnover(ipts,ivm,ifruit,iele) + & 
                        litterfrac(ipts,ivm,icarbres,ilitt) * turnover(ipts,ivm,icarbres,iele)+ &  
                        litterfrac(ipts,ivm,ilabile,ilitt) * turnover(ipts,ivm,ilabile,iele)
                IF (litter_inc(ipts,ilitt,ivm,iabove,iele).LT.10*EPSILON(un)) litter_inc(ipts,ilitt,ivm,iabove,iele) = zero
                litter_inc(ipts,ilitt,ivm,ibelow,iele) = & 
                        litterfrac(ipts,ivm,isapbelow,ilitt) * &
                        (bm_to_litter(ipts,ivm,isapbelow,iele)-tree_bm_to_litter(ipts,ivm,isapbelow,iele)) + &
                        litterfrac(ipts,ivm,iheartbelow,ilitt) * &
                        (bm_to_litter(ipts,ivm,iheartbelow,iele)-tree_bm_to_litter(ipts,ivm,iheartbelow,iele)) + &  
                        tree_litterfrac(ilitt) * tree_bm_to_litter(ipts,ivm,isapbelow,iele) + &  
                        tree_litterfrac(ilitt) * tree_bm_to_litter(ipts,ivm,iheartbelow,iele) + &        
                        litterfrac(ipts,ivm,iroot,ilitt) * bm_to_litter(ipts,ivm,iroot,iele) + & 
                        litterfrac(ipts,ivm,isapbelow,ilitt) * turnover(ipts,ivm,isapbelow,iele) + & 
                        litterfrac(ipts,ivm,iheartbelow,ilitt) * turnover(ipts,ivm,iheartbelow,iele) + & 
                        litterfrac(ipts,ivm,iroot,ilitt) * turnover(ipts,ivm,iroot,iele)    
                IF (litter_inc(ipts,ilitt,ivm,ibelow,iele).LT.10*EPSILON(un)) litter_inc(ipts,ilitt,ivm,ibelow,iele) = zero

                ! From the calculation of litter_inc above, the labile,reserve and fruit
                ! mass were treated completely as aboveground, we will do the same for fuel.
                ! Note bm_to_litter is used here only for diagnositc purpose to calculate
                ! surface litter fuel of different fuel categories. Nothing below concern with 
                ! mass balance check. For the same reason, litterfuel is not included in mass balance 
                ! check.

                IF (is_tree(ivm)) THEN
                  ! we treat tree and herbaceous differently. For tree, the components of
                  ! leaf and fruit go 100% to 1h fuel category. Other components were
                  ! splited into the four classes of 1h/10h/100h/1000h.
                  ! alloc_firefuel of different fuel classes, defined as constant, add exactly to one
                  DO ihour = 1,nhour
                      IF ( ok_snags ) THEN
                         ! Snag should be excluded from the litter fuel pool,
                         ! as they do not contribute to fire spread or ignition.
                         ! Hence, we set the tree_litterfrac(isnag) as zero
                         ! before calculating the litter fuel. when we get the litter
                         ! fuel, we set the tree_litterfrac(isnag) as previous
                         ! value. 
                         tree_litterfrac(isnag) = zero
                         litterfuel_inc(ipts,ilitt,ivm,ihour,iele) =alloc_firefuel(ihour) *&
                              ( litterfrac(ipts,ivm,isapabove,ilitt) * &
                              (bm_to_litter(ipts,ivm,isapabove,iele)-tree_bm_to_litter(ipts,ivm,isapabove,iele))+ &
                              litterfrac(ipts,ivm,iheartabove,ilitt) * &
                              (bm_to_litter(ipts,ivm,iheartabove,iele)-tree_bm_to_litter(ipts,ivm,iheartabove,iele)) + &
                              tree_litterfrac(ilitt) * tree_bm_to_litter(ipts,ivm,isapabove,iele) + &
                              tree_litterfrac(ilitt) * tree_bm_to_litter(ipts,ivm,iheartabove,iele) + &
                              litterfrac(ipts,ivm,icarbres,ilitt) * bm_to_litter(ipts,ivm,icarbres,iele) + &
                              litterfrac(ipts,ivm,ilabile,ilitt) * bm_to_litter(ipts,ivm,ilabile,iele) + &
                              litterfrac(ipts,ivm,isapabove,ilitt) * turnover(ipts,ivm,isapabove,iele) + &
                              litterfrac(ipts,ivm,iheartabove,ilitt) * turnover(ipts,ivm,iheartabove,iele) + &
                              litterfrac(ipts,ivm,icarbres,ilitt) * turnover(ipts,ivm,icarbres,iele)+ &
                              litterfrac(ipts,ivm,ilabile,ilitt) * turnover(ipts,ivm,ilabile,iele) )
                         tree_litterfrac(isnag) = snag_frac
                      ELSE
                         litterfuel_inc(ipts,ilitt,ivm,ihour,iele) = alloc_firefuel(ihour) *&
                              ( litterfrac(ipts,ivm,isapabove,ilitt) * &
                              (bm_to_litter(ipts,ivm,isapabove,iele)-tree_bm_to_litter(ipts,ivm,isapabove,iele)) + &
                              litterfrac(ipts,ivm,iheartabove,ilitt) * &
                              (bm_to_litter(ipts,ivm,iheartabove,iele)-tree_bm_to_litter(ipts,ivm,iheartabove,iele)) + &
                              tree_litterfrac(ilitt) * tree_bm_to_litter(ipts,ivm,isapabove,iele) + &
                              tree_litterfrac(ilitt) * tree_bm_to_litter(ipts,ivm,iheartabove,iele) + &
                              litterfrac(ipts,ivm,icarbres,ilitt) * bm_to_litter(ipts,ivm,icarbres,iele) + &
                              litterfrac(ipts,ivm,ilabile,ilitt) * bm_to_litter(ipts,ivm,ilabile,iele) + &
                              litterfrac(ipts,ivm,isapabove,ilitt) * turnover(ipts,ivm,isapabove,iele) + &
                              litterfrac(ipts,ivm,iheartabove,ilitt) * turnover(ipts,ivm,iheartabove,iele) + &
                              litterfrac(ipts,ivm,icarbres,ilitt) * turnover(ipts,ivm,icarbres,iele)+ &
                              litterfrac(ipts,ivm,ilabile,ilitt) * turnover(ipts,ivm,ilabile,iele) )
                      ENDIF
                      IF (litterfuel_inc(ipts,ilitt,ivm,ihour,iele).LT.10*EPSILON(un)) THEN
                          litterfuel_inc(ipts,ilitt,ivm,ihour,iele) = zero
                      ENDIF
                  ENDDO
                  
                  ! leaf and fruit go 100% to 1h fuel class.
                  litterfuel_inc(ipts,ilitt,ivm,ihour1,iele) = litterfuel_inc(ipts,ilitt,ivm,ihour1,iele) + &
                      litterfrac(ipts,ivm,ileaf,ilitt) * bm_to_litter(ipts,ivm,ileaf,iele) + &
                      litterfrac(ipts,ivm,ifruit,ilitt) * bm_to_litter(ipts,ivm,ifruit,iele) + &
                      litterfrac(ipts,ivm,ileaf,ilitt) * turnover(ipts,ivm,ileaf,iele) + &
                      litterfrac(ipts,ivm,ifruit,ilitt) * turnover(ipts,ivm,ifruit,iele)
                  IF (litterfuel_inc(ipts,ilitt,ivm,ihour1,iele).LT.10*EPSILON(un)) THEN
                      litterfuel_inc(ipts,ilitt,ivm,ihour1,iele) = zero
                  ENDIF
                ELSE
                   ! for herbaceous PFTs, all aboveground litter goes to 1h fuel.
                   litterfuel_inc(ipts,ilitt,ivm,ihour1,iele) = litter_inc(ipts,ilitt,ivm,iabove,iele)
                   IF (litterfuel_inc(ipts,ilitt,ivm,ihour1,iele).LT.10*EPSILON(un)) THEN
                       litterfuel_inc(ipts,ilitt,ivm,ihour1,iele) = zero
                   ENDIF
                ENDIF
             ENDDO
             
             IF ( ilitt /= iwoody .AND. ilitt /= isnag ) THEN
                !! 2.2.3 dead leaves, for soil cover.
                dead_leaves(ipts,ivm,ilitt) = &
                        dead_leaves(ipts,ivm,ilitt) + litterfrac(ipts,ivm,ileaf,ilitt) * &
                        (bm_to_litter(ipts,ivm,ileaf,icarbon) + &
                        turnover(ipts,ivm,ileaf,icarbon) )
             ENDIF
             IF(ivm == ibare_sechiba.OR.is_tree(ivm))THEN
                ! Note that a bare soil pixel may have some left 
                ! overs of the litter of a previous forest. So
                ! bare soil is treated as forest here.
                IF ( ilitt == istructural ) THEN
                   !! 2.2.4 lignin increase in structural litter
                   lignin_struc_inc(ipts,ivm,iabove) = &
                           LC(ivm,ileaf) * bm_to_litter(ipts,ivm,ileaf,icarbon) + &
                           LC(ivm,ifruit) * bm_to_litter(ipts,ivm,ifruit,icarbon) + &
                           LC(ivm,icarbres) * bm_to_litter(ipts,ivm,icarbres,icarbon) + &
                           LC(ivm,ilabile) * bm_to_litter(ipts,ivm,ilabile,icarbon) + &
                           LC(ivm,ileaf) * turnover(ipts,ivm,ileaf,icarbon) + &
                           LC(ivm,ifruit) * turnover(ipts,ivm,ifruit,icarbon) + &
                           LC(ivm,icarbres) * turnover(ipts,ivm,icarbres,icarbon) + &
                           LC(ivm,ilabile) * turnover(ipts,ivm,ilabile,icarbon)
                   lignin_struc_inc(ipts,ivm,ibelow) = &
                           LC(ivm,iroot) * bm_to_litter(ipts,ivm,iroot,icarbon) + &
                           LC(ivm,iroot) * turnover(ipts,ivm,iroot,icarbon)
                ELSEIF ( ilitt == iwoody ) THEN       
                   !! 2.2.4 lignin increase in woody litter
                   IF ( ok_snags ) THEN 
                      lignin_wood_inc(ipts,ivm,iabove) = (1 - snag_frac) * &
                              ( LC(ivm,isapabove) * turnover(ipts,ivm,isapabove,icarbon) + &
                              LC(ivm,iheartabove) * turnover(ipts,ivm,iheartabove,icarbon) )
                      lignin_wood_inc(ipts,ivm,ibelow) = (1 - snag_frac) * &
                              ( LC(ivm,isapbelow) * turnover(ipts,ivm,isapbelow,icarbon) + &
                              LC(ivm,iheartbelow) * turnover(ipts,ivm,iheartbelow,icarbon) )
                      ELSE 
                      lignin_wood_inc(ipts,ivm,iabove) = &
                              LC(ivm,isapabove) * turnover(ipts,ivm,isapabove,icarbon) + &
                              LC(ivm,iheartabove) * turnover(ipts,ivm,iheartabove,icarbon) 
                      lignin_wood_inc(ipts,ivm,ibelow) = &
                              LC(ivm,isapbelow) * turnover(ipts,ivm,isapbelow,icarbon) + &
                              LC(ivm,iheartbelow) * turnover(ipts,ivm,iheartbelow,icarbon) 
                   ENDIF
                ELSEIF ( ilitt == isnag ) THEN       
                   !! 2.2.4 lignin increase in snag litter
                   IF ( ok_snags ) THEN 
                      lignin_snag_inc(ipts,ivm,iabove) = snag_frac * &
                              ( LC(ivm,isapabove) * turnover(ipts,ivm,isapabove,icarbon) + &
                              LC(ivm,iheartabove) * turnover(ipts,ivm,iheartabove,icarbon) )
                      lignin_snag_inc(ipts,ivm,ibelow) = snag_frac * &
                              ( LC(ivm,isapbelow) * turnover(ipts,ivm,isapbelow,icarbon) + &
                              LC(ivm,iheartbelow) * turnover(ipts,ivm,iheartbelow,icarbon) )
                   ELSE 
                      lignin_snag_inc(ipts,ivm,iabove) = zero
                      lignin_snag_inc(ipts,ivm,ibelow) = zero
                   ENDIF
                ENDIF ! istructural or iwoody or isnag
             ELSE 
                ! None woody vegetation 
                IF (ilitt == istructural) THEN
                   !! 2.2.4 lignin increase in structural litter
                   lignin_struc_inc(ipts,ivm,iabove) = lignin_struc_inc(ipts,ivm,iabove) + &
                           LC(ivm,isapabove) * turnover(ipts,ivm,isapabove,icarbon) + &
                           LC(ivm,iheartabove) * turnover(ipts,ivm,iheartabove,icarbon) 
                   lignin_struc_inc(ipts,ivm,ibelow) = lignin_struc_inc(ipts,ivm,ibelow) + &
                           LC(ivm,isapbelow)*turnover(ipts,ivm,isapbelow,icarbon) + &
                           LC(ivm,iheartbelow)*turnover(ipts,ivm,iheartbelow,icarbon)
                ENDIF ! istructural
             ENDIF ! is_tree
             IF (ilitt == istructural) THEN
                IF ( ok_snags ) THEN
                   ! above woody / snag
                   temp_flux = & 
                        LC(ivm,isapabove) * tree_bm_to_litter(ipts,ivm,isapabove,icarbon) + &
                        LC(ivm,iheartabove) * tree_bm_to_litter(ipts,ivm,iheartabove,icarbon) 
                   lignin_wood_inc(ipts,ivm,iabove) =  (1 - snag_frac) * temp_flux
                   lignin_snag_inc(ipts,ivm,iabove) = snag_frac * temp_flux
                   ! below woody / snag
                   temp_flux = &
                        LC(ivm,isapbelow) * tree_bm_to_litter(ipts,ivm,isapbelow,icarbon) + &
                        LC(ivm,iheartbelow) * tree_bm_to_litter(ipts,ivm,iheartbelow,icarbon) 
                   lignin_wood_inc(ipts,ivm,ibelow) = (1 - snag_frac)  * temp_flux
                   lignin_snag_inc(ipts,ivm,ibelow) = snag_frac * temp_flux
                ELSE 
                   lignin_wood_inc(ipts,ivm,iabove) = &
                        LC(ivm,isapabove) * tree_bm_to_litter(ipts,ivm,isapabove,icarbon) + &
                        LC(ivm,iheartabove) * tree_bm_to_litter(ipts,ivm,iheartabove,icarbon)
                   lignin_wood_inc(ipts,ivm,ibelow) = &
                        LC(ivm,isapbelow) * tree_bm_to_litter(ipts,ivm,isapbelow,icarbon) + &
                        LC(ivm,iheartbelow) * tree_bm_to_litter(ipts,ivm,iheartbelow,icarbon) 
                   lignin_snag_inc(ipts,ivm,iabove) = zero
                   lignin_snag_inc(ipts,ivm,ibelow) = zero
                ENDIF
                lignin_struc_inc(ipts,ivm,iabove) = lignin_struc_inc(ipts,ivm,iabove) + &
                        LC(ivm,isapabove) * (bm_to_litter(ipts,ivm,isapabove,icarbon)-tree_bm_to_litter(ipts,ivm,isapabove,icarbon)) + &
                        LC(ivm,iheartabove) * (bm_to_litter(ipts,ivm,iheartabove,icarbon)-tree_bm_to_litter(ipts,ivm,iheartabove,icarbon))
                lignin_struc_inc(ipts,ivm,ibelow) = lignin_struc_inc(ipts,ivm,ibelow) + &
                        LC(ivm,isapbelow) * (bm_to_litter(ipts,ivm,isapbelow,icarbon)-tree_bm_to_litter(ipts,ivm,isapbelow,icarbon)) + &
                        LC(ivm,iheartbelow) * (bm_to_litter(ipts,ivm,iheartbelow,icarbon)-tree_bm_to_litter(ipts,ivm,iheartbelow,icarbon))
             END IF
          ENDDO
       ENDDO 
    ENDDO !#nvm

    err_max_1 = zero
    DO iele = 1,nelements
       DO ivm = 1,nvm
          DO ipts = 1,npts
             err_max_1 = MAX(err_max_1,ABS(SUM(litter_inc(ipts,:,ivm,:,iele)) - &
                     SUM(bm_to_litter(ipts,ivm,:,iele)) - &
                     SUM(turnover(ipts,ivm,:,iele))))
          END DO 
       END DO 
    ENDDO

    IF(err_max_1.GT.10000*EPSILON(un)) THEN
       CALL ipslerr_p(plev,'littercalc','more than 1e-10 carbon or nitrogen', &
            'has gone missing','this might hint at a mass balance problem')
    ENDIF
    
    litterfuel(:,:,:,:,:) = litterfuel(:,:,:,:,:) + litterfuel_inc(:,:,:,:,:)
    !! 2.2.5 add new litter (struct/met, above/below)
    ! litter_inc should be a positive number but due to all the multiplications
    ! above it is possible that very small negative numbers are generated (10e-21)
    ! When added to a very small number for litter, the total litter could become
    ! a negative number representing zero. That could fail some of the quality 
    ! tests further down. This is a safe place to enhance numerical consistency of 
    ! the code.
    !! 2.2.6 for security: can't add more lignin than structural 
    !! litter (above/below)

    DO ilev = 1, nlevs !Loop over litter levels (above and below ground)
       DO ivm = 1,nvm !Loop over PFTs
          DO ipts = 1,npts
             old_struc = litter(ipts,istructural,ivm,ilev,icarbon)
             old_woody = litter(ipts,iwoody,ivm,ilev,icarbon)
             old_snag = litter(ipts,isnag,ivm,ilev,icarbon)
             litter(ipts,:,ivm,ilev,:) = litter(ipts,:,ivm,ilev,:) + litter_inc(ipts,:,ivm,ilev,:)
             lignin_struc_inc(ipts,ivm,ilev) = MIN( lignin_struc_inc(ipts,ivm,ilev), &
                     litter_inc(ipts,istructural,ivm,ilev,icarbon) )
             lignin_wood_inc(ipts,ivm,ilev) = MIN( lignin_wood_inc(ipts,ivm,ilev), &
                     litter_inc(ipts,iwoody,ivm,ilev,icarbon) )
             lignin_snag_inc(ipts,ivm,ilev) = MIN( lignin_snag_inc(ipts,ivm,ilev), &
                     litter_inc(ipts,isnag,ivm,ilev,icarbon) )
             !! 2.2.7 new lignin content: add old lignin and lignin increase, divide by 
             !!       total structural litter (above/below)
             IF( litter(ipts,istructural,ivm,ilev,icarbon) .GT. min_stomate ) THEN
                lignin_struc(ipts,ivm,ilev) = ( lignin_struc(ipts,ivm,ilev) * &
                        old_struc + lignin_struc_inc(ipts,ivm,ilev) ) / &
                        litter(ipts,istructural,ivm,ilev,icarbon)
             ELSE
                lignin_struc(ipts,ivm,ilev) = zero
             ENDIF
             IF ( litter(ipts,iwoody,ivm,ilev,icarbon) .GT. min_stomate ) THEN
                lignin_wood(ipts,ivm,ilev) = ( lignin_wood(ipts,ivm,ilev) * &
                        old_woody + lignin_wood_inc(ipts,ivm,ilev) ) / &
                        litter(ipts,iwoody,ivm,ilev,icarbon)
             ELSE
                lignin_wood(ipts,ivm,ilev) = zero
             ENDIF
             IF ( litter(ipts,isnag,ivm,ilev,icarbon) .GT. min_stomate ) THEN
                lignin_snag(ipts,ivm,ilev) = ( lignin_snag(ipts,ivm,ilev) * &
                        old_snag + lignin_snag_inc(ipts,ivm,ilev) ) / &
                        litter(ipts,isnag,ivm,ilev,icarbon)
             ELSE
                lignin_snag(ipts,ivm,ilev) = zero
             ENDIF
          ENDDO
       ENDDO ! ivm
    ENDDO ! ilev
    
    !! 2.2.8 Add the carbon of the manure
    ! The N-maps only contain the nitrogen contained in the manure. In reality
    ! manure contains a lot of carbon. Here this carbon, which is absent from
    ! the maps, is added to the litter pools. Add manure into the metabolic
    ! above ground pool if enough C is harvested over the pixel. If not enough
    ! C is harvested a part goes as litter and the other part goes at mineral
    ! N. Here we assume that N is not limiting, only the C harvested is limiting
    ! Note that here we only treat input for imanure because it is added to the 
    ! litter pool. The other nitrogen inputs are added to the soil and therefore 
    ! taken care of in stomate_soilcarbon (nitrogen_dynamics).
   
    ! SpitfireNote: we tentatively add manure into 1h fuel class as manure
    ! often has small size. But in fact currently spitfire handles only fires over 
    ! natural vegetation (dictated by the parameter `burnable`). Fires over
    ! intensively managed land is possible but these fires are also well under
    ! human control and not suitable for fire propagation simulation developed 
    ! based on naturally uncontroled burning.

    DO ivm=1,nvm
       DO ipts = 1,npts
          IF((veget_max(ipts,ivm).GT.min_stomate) .AND. &
               (SUM(harvest_pool_acc(ipts,ivm,1,icarbon,:)) .GE. &
               n_input(ipts,ivm,month,imanure)*cn_ratio_manure*dt*veget_max(ipts,ivm)*area(ipts)*contfrac(ipts)))THEN
             ! The harvest pools can provide all the carbon
             litter(ipts,imetabolic,ivm,iabove,icarbon) = litter(ipts,imetabolic,ivm,iabove,icarbon) + &
                  n_input(ipts,ivm,month,imanure)*cn_ratio_manure*dt
             litter(ipts,imetabolic,ivm,iabove,initrogen) = litter(ipts,imetabolic,ivm,iabove,initrogen) + &
                  n_input(ipts,ivm,month,imanure)*dt
             litterfuel(ipts,imetabolic,ivm,ihour1,icarbon) = litterfuel(ipts,imetabolic,ivm,ihour1,icarbon) + &
                  n_input(ipts,ivm,month,imanure)*cn_ratio_manure*dt
             litterfuel(ipts,imetabolic,ivm,ihour1,initrogen) =litterfuel(ipts,imetabolic,ivm,ihour1,initrogen) + &
                  n_input(ipts,ivm,month,imanure)*dt

             
             ! Take the carbon from the harvest pools
             IF(harvest_pool_acc(ipts,ivm,1,icarbon,iharvest) .GT. &
                  n_input(ipts,ivm,month,imanure)*cn_ratio_manure*dt*veget_max(ipts,ivm)*area(ipts)*contfrac(ipts))THEN
                ! There is enough carbon in the harvest pool. Take all the carbon from iharvest because this is likely
                ! the largest pool.
                harvest_pool_acc(ipts,ivm,1,icarbon,iharvest) = harvest_pool_acc(ipts,ivm,1,icarbon,iharvest) - &
                     n_input(ipts,ivm,month,imanure)*cn_ratio_manure*dt*veget_max(ipts,ivm)*area(ipts)*contfrac(ipts)
             ELSE
                ! There is not enough carbon in iharvest so take from iharvest and ilcc
                harvest_pool_acc(ipts,ivm,1,icarbon,ilcc) = harvest_pool_acc(ipts,ivm,1,icarbon,ilcc) - &
                     (n_input(ipts,ivm,month,imanure)*cn_ratio_manure*dt*veget_max(ipts,ivm)*area(ipts)*contfrac(ipts) - &
                     harvest_pool_acc(ipts,ivm,1,icarbon,iharvest))
                harvest_pool_acc(ipts,ivm,1,icarbon,iharvest) = zero
             END IF
             ! Account for the manure when litter decomposition will have to be recalculated
             ! further down this subroutine
             total_init_nitrogen(ipts,ivm) = total_init_nitrogen(ipts,ivm) + n_input(ipts,ivm,month,imanure)*dt
          ELSEIF((veget_max(ipts,ivm).GT.min_stomate) .AND. &
               (SUM(harvest_pool_acc(ipts,ivm,1,icarbon,:)) .LT. &
               n_input(ipts,ivm,month,imanure)*cn_ratio_manure*dt*veget_max(ipts,ivm)*area(ipts)*contfrac(ipts)))THEN
             ! There is not enough carbon in the harvest pools. Empty the harvest pools and take the re
             

             litter(ipts,imetabolic,ivm,iabove,icarbon) = litter(ipts,imetabolic,ivm,iabove,icarbon) + &
                  SUM(harvest_pool_acc(ipts,ivm,1,icarbon,:))/veget_max(ipts,ivm)/area(ipts)/contfrac(ipts)
             litter(ipts,imetabolic,ivm,iabove,initrogen) = litter(ipts,imetabolic,ivm,iabove,initrogen) + &
                  SUM(harvest_pool_acc(ipts,ivm,1,icarbon,:))/veget_max(ipts,ivm)/area(ipts)/contfrac(ipts)/cn_ratio_manure
             litterfuel(ipts,imetabolic,ivm,ihour1,icarbon) = litterfuel(ipts,imetabolic,ivm,ihour1,icarbon) + &
                  SUM(harvest_pool_acc(ipts,ivm,1,icarbon,:))/veget_max(ipts,ivm)/area(ipts)/contfrac(ipts)
             litterfuel(ipts,imetabolic,ivm,ihour1,initrogen) = litterfuel(ipts,imetabolic,ivm,ihour1,initrogen) + &
                  SUM(harvest_pool_acc(ipts,ivm,1,icarbon,:))/veget_max(ipts,ivm)/area(ipts)/contfrac(ipts)/cn_ratio_manure

             soil_n_min(ipts,ivm,iammonium) = soil_n_min(ipts,ivm,iammonium) + &
                  (n_input(ipts,ivm,month,imanure)*dt - &
                  SUM(harvest_pool_acc(ipts,ivm,1,icarbon,:))/veget_max(ipts,ivm)/area(ipts)/contfrac(ipts)/cn_ratio_manure)*&
                  ratio_nh4_fert
             soil_n_min(ipts,ivm,initrate) = soil_n_min(ipts,ivm,initrate) + &
                  (n_input(ipts,ivm,month,imanure)*dt - &
                  SUM(harvest_pool_acc(ipts,ivm,1,icarbon,:))/veget_max(ipts,ivm)/area(ipts)/contfrac(ipts)/cn_ratio_manure)*&
                  (1.-ratio_nh4_fert)
             ! Account for the manure when litter decomposition will have to be recalculated
             ! further down this subroutine
             total_init_nitrogen(ipts,ivm) = total_init_nitrogen(ipts,ivm) + &
                  SUM(harvest_pool_acc(ipts,ivm,1,icarbon,:))/veget_max(ipts,ivm)/area(ipts)/contfrac(ipts)/cn_ratio_manure 
             harvest_pool_acc(ipts,ivm,1,icarbon,:) = zero
          ELSE
             ! Too small fraction, the n_input is not used
             N_input(ipts,ivm,month,imanure)=zero
          ENDIF
       ENDDO
    ENDDO
     
    !! 2.2.9 Add a given fraction of the snag pool to the woody pool according
    ! to an exponential model decay (Hararuk et al. 2020)
    snag_to_wood_flux(:,:,:,:) = zero
    IF ( ok_snags ) THEN
       DO ivm=1,nvm
          DO ipts = 1,npts
             IF (( veget_max(ipts,ivm).GT.min_stomate ) .AND. is_tree(ivm)) THEN
                snag_to_wood_flux(ipts,ivm,:,:) = dt * snag_fall_fac * litter(ipts,isnag,ivm,:,:)
                WHERE ( snag_to_wood_flux(ipts,ivm,:,:).LE.litter(ipts,isnag,ivm,:,:) )
                      litter(ipts,iwoody,ivm,:,:) = litter(ipts,iwoody,ivm,:,:) + snag_to_wood_flux(ipts,ivm,:,:)
                      litter(ipts,isnag,ivm,:,:) = litter(ipts,isnag,ivm,:,:) - snag_to_wood_flux(ipts,ivm,:,:)
                ENDWHERE
                WHERE ( snag_to_wood_flux(ipts,ivm,:,:).GT.litter(ipts,isnag,ivm,:,:) )
                      snag_to_wood_flux(ipts,ivm,:,:) = litter(ipts,isnag,ivm,:,:)
                      litter(ipts,iwoody,ivm,:,:) = litter(ipts,iwoody,ivm,:,:) + snag_to_wood_flux(ipts,ivm,:,:)
                      litter(ipts,isnag,ivm,:,:) = zero
                ENDWHERE
             ENDIF
          ENDDO
       ENDDO
    ENDIF
  
    !! 3. Temperature control on decay: Factor between 0 and 1
    !! 3.1 above: surface temperature
     
    DO ipts = 1,npts
        control_temp(ipts,iabove) = exp( soil_Q10 * ( tsurf(ipts) - (ZeroCelsius+tsoil_ref)) / Q10 )
        control_temp(ipts,iabove) = MIN( un, control_temp(ipts,iabove) )
        !! 3.2 below: convolution of temperature and decomposer profiles
        !! exponential decomposer profile supposed)
        !! 3.2.1 rpc is an integration constant such that the integral of 
        !  the decomposer profile becomes one. This looks like the calculation
        !  of a root profile but it uses a different constant; z_decomp instead 
        !  of humncste for the root profile. The integral goes from the top
        !  of the soil all the way to the bottom.
        !  All pixels have the same zdr profile and all PFTs and pixels
        !  are assumed to have the same z_decomp so rpc is globally fixed.
        rpc = un / (EXP(-zdr(0)/z_decomp) - EXP(-zdr(nslm)/z_decomp))

        !! 3.2.2 integrate over the nslm levels
        tsoil_decomp(ipts) = zero
        DO ibdl = 1, nslm
        
           ! Calculate the soil temperature weighted by the assumed presence of
           ! the decomposers. Presence of the decomposers is prescribes by the 
           ! presence = e^(-soil_depth/z_decomp) where z_decomp described the 
           ! shape of an exponentially decaying function with soil depth. stempdiag
           ! is descritized along the nodes of the soil layers following the scheme
           ! of De Rosnay (PhD, figure C.2 page 156), this has been calculated in 
           ! vertical_soil.f90 as zdr.
           tsoil_decomp(ipts) = &
                   tsoil_decomp(ipts) + stempdiag(ipts,ibdl) * rpc * &
                   ( EXP(-zdr(ibdl-1)/z_decomp) - EXP(-zdr(ibdl)/z_decomp) )

        ENDDO
        control_temp(ipts,ibelow) = exp( soil_Q10 * ( tsoil_decomp(ipts) - (ZeroCelsius+tsoil_ref)) / Q10 )
        control_temp(ipts,ibelow) = MIN( un, control_temp(ipts,ibelow) ) 
        !! 4. Moisture control. Factor between 0 and 1
        !! 4.1 above the ground: litter humidity
        control_moist(ipts,iabove) = -moist_coeff(1) * litterhum(ipts) * litterhum(ipts) &
                +moist_coeff(2) * litterhum(ipts) - moist_coeff(3)
        control_moist(ipts,iabove) = MAX( moistcont_min, MIN( un, control_moist(ipts,iabove) ) )
     
        !! 4.2 below: convolution of humidity and decomposer profiles
        !            (exponential decomposer profile supposed)

        !! 4.2.2 integrate over the nslm levels
        soilhum_decomp = zero
        DO ibdl = 1, nslm !Loop over soil levels
        
           ! Calculate the soil humidity weighted by the assumed presence of
           ! the decomposers. Presence of the decomposers is prescribes by the 
           ! presence = e^(-soil_depth/z_decomp) where z_decomp described the 
           ! shape of an exponentially decaying function with soil depth. stempdiag
           ! is descritized along the nodes of the soil layers following the scheme
           ! of De Rosnay (PhD, figure C.2 page 156), this has been calculated in 
           ! vertical_soil.f90 as zdr
           soilhum_decomp = &
                   soilhum_decomp + shumdiag(ipts,ibdl) * rpc * &
                   ( EXP(-zdr(ibdl-1)/z_decomp) - EXP(-zdr(ibdl)/z_decomp) )
        
        ENDDO
     
        control_moist(ipts,ibelow) = -moist_coeff(1) * soilhum_decomp * soilhum_decomp &
                +moist_coeff(2) * soilhum_decomp - moist_coeff(3)
        control_moist(ipts,ibelow) = MAX( moistcont_min, MIN( un, control_moist(ipts,ibelow) ) )

     ENDDO

     ! peatland
     ! The below codes are not put into npts loop as above
     ! because we call several different functions below, and those functions are more complex than 
     ! the two lines of code above for computing control_moist
     IF (perma_peat) THEN
        !peatland : use water content instead of relative hum, xw
        control_moist_peat(:,iabove)=control_moist_func_peat(npts, mc_peat_above)
     ENDIF

     IF (agri_peat) THEN
       control_moist_croppeat(:,iabove)=control_moist_func_peat(npts, mc_croppeat_above)
     ENDIF

     IF (perma_peat .AND. tides) THEN
       control_moist_man(:,iabove)=control_moist_func_man(npts, mc_man_above)
     ENDIF

     ! peatland: use water content instead of relative hum
     IF (perma_peat) THEN
        ! mc_peat_below has the dimension of npts
        mc_peat_below(:) =zero
        ! added rpc definition again
        rpc = un / (EXP(-zdr(0)/z_decomp) - EXP(-zdr(nslm)/z_decomp))
        ! 
        DO ibdl = 1, nslm !Loop over soil levels
           mc_peat_below(:) = &
                mc_peat_below(:)+shumdiag_peat(:,ibdl)  * rpc * & 
                ( EXP( -zdr(ibdl-1)/z_decomp ) - EXP( -zdr(ibdl)/z_decomp ) )
        ENDDO
        control_moist_peat(:,ibelow)=control_moist_func_peat(npts, mc_peat_below)
        
     ENDIF
     
     IF (perma_peat .AND. tides ) THEN
        mc_man_below(:) =zero
        rpc = un / (EXP(-zdr(0)/z_decomp) - EXP(-zdr(nslm)/z_decomp))
        DO ibdl = 1, nslm !Loop over soil levels
           mc_man_below(:) = &
                mc_man_below(:)+shumdiag_man(:,ibdl)  * rpc * &  
                ( EXP( -zdr(ibdl-1)/z_decomp ) - EXP( -zdr(ibdl)/z_decomp ) )
        ENDDO
        control_moist_man(:,ibelow)=control_moist_func_man(npts, mc_man_below)
     ENDIF

     ! crops
     IF (agri_peat) THEN
        mc_croppeat_below(:)=zero
        rpc = un / (EXP(-zdr(0)/z_decomp) - EXP(-zdr(nslm)/z_decomp))
        DO ibdl = 1, nslm !Loop over soil levels
           mc_croppeat_below(:)=&
                mc_croppeat_below(:)+shumdiag_croppeat(:,ibdl)  * rpc * & 
                ( EXP( -zdr(ibdl-1)/z_decomp ) - EXP( -zdr(ibdl)/z_decomp ) )
        ENDDO
        control_moist_croppeat(:,ibelow)=control_moist_func_peat(npts,mc_croppeat_below)
     ENDIF
   
     ! 5. fluxes from litter to carbon pools and respiration
     err_min_1 = zero
     err_min_2 = zero
     err_min_3 = zero
     err_min_4 = zero
     err_min_5 = zero
     err_min_6 = zero
     err_min_7 = zero
     err_min_8 = zero
     DO ivm = 1,nvm 
        ! define the soil tile type
        istm = pref_soil_veg(ivm)
        DO ipts = 1,npts
           DO ilev = 1, nlevs
              IF ( ok_soil_carbon_discretization ) THEN
                 itarget=iactive
              ELSE
                 IF (ilev.EQ.iabove)THEN 
                    itarget=isurface 
                 ELSE 
                    itarget=iactive 
                 ENDIF
              END IF

              !! 5.1 structural litteripts goes into active and slow carbon 
              !  pools + respiration

              !! 5.1.1 total quantity of structural litter which is decomposed
              !
              IF (perma_peat .AND. is_peat(ivm) .AND. soiltile(ipts,istm) .GT. zero ) THEN
                 fd = dt*litter_dec_fac(istructural) * &
                      control_temp(ipts,ilev) * control_moist_peat(ipts,ilev) * &
                      exp( -litter_struct_coef * lignin_struc(ipts,ivm,ilev) )
              ELSE
                 fd = dt*litter_dec_fac(istructural) * &
                      control_temp(ipts,ilev) * control_moist(ipts,ilev) * &
                      exp( -litter_struct_coef * lignin_struc(ipts,ivm,ilev) )
              ENDIF
              ! Ensure that fd is between 0 and 1
              fd = MAX(zero,MIN(fd,un))

              !! decompose same fraction of structural part of dead leaves. 
              ! Not exact as lignin content is not the same as that of the total 
              ! structural litter. to avoid a multiple (for ibelow and iabove) 
              ! modification of dead_leaves, we do this test to do this calculate 
              ! only ones in 1,nlev loop
              IF (ilev == iabove)  THEN
              
                 dead_leaves(ipts,ivm,istructural) = &
                         dead_leaves(ipts,ivm,istructural) * ( un - fd )
              
              ENDIF
           
              !! 5.1.2 Calculate the amount of structural litter
              !  that was decomposed
              
              ! Calculate the quantity of structural litter that should
              ! be decomposed
              qd_icarbon = litter(ipts,istructural,ivm,ilev,icarbon) * fd
              qd_initrogen = litter(ipts,istructural,ivm,ilev,initrogen) * fd

              ! Keep track of how many errors there are
              err_min_1 = MIN(err_min_1, MIN(litter(ipts,istructural,ivm,ilev,icarbon), &
                   litter(ipts,istructural,ivm,ilev,initrogen)))
            
              ! Subtract the mineralized fraction from the litter but make sure
              ! the mineralisation cannot become negative
              
              ! qd is calculated as litter(:,istructural,ivm,ilev,iele)*fd(:)
              ! as long as fd is positive (which was enforced above) the
              ! subtraction below should never result in a number that is
              ! less than zero (except for precision issues). fd is the same
              ! for carbon and nitrogen and therefore the ratio is preserved
              ! during mineralization. The test conditions should be
              ! satisfied for C and N at the same time.
              err_min_2 = MIN(err_min_2,MIN((litter(ipts,istructural,ivm,ilev,icarbon) - qd_icarbon), &
                   (litter(ipts,istructural,ivm,ilev,initrogen)-qd_initrogen)))
              litter(ipts,istructural,ivm,ilev,icarbon) = MAX(zero, &
                   litter(ipts,istructural,ivm,ilev,icarbon) - qd_icarbon)
              litter(ipts,istructural,ivm,ilev,initrogen) = MAX(zero, &
                   litter(ipts,istructural,ivm,ilev,initrogen) - qd_initrogen)
              n_mineralisation(ipts,ivm) = n_mineralisation(ipts,ivm) + qd_initrogen

              ! Debug
              IF ( ( printlev_loc >= 4 ).AND.( ivm == test_pft ).AND.( ipts == test_grid ) ) THEN
                 WRITE(numout,*) 'het_resp before struc decomp', resp_hetero_litter(test_grid,test_pft)
              ENDIF

              !! 5.1.3 non-lignin fraction of structural litter goes into
              !!       active (or surface) carbon pool + respiration
              decomp_residual = qd_icarbon
              som_input(ipts,itarget,ivm,icarbon) = som_input(ipts,itarget,ivm,icarbon) + frac_soil(istructural,itarget,ilev) * &
                   qd_icarbon * ( 1. - lignin_struc(ipts,ivm,ilev) ) / dt 
              decomp_residual = decomp_residual - frac_soil(istructural,itarget,ilev) * &
                   qd_icarbon * ( 1. - lignin_struc(ipts,ivm,ilev) )  

              ! Here resp_hetero_litter is divided by dt to 
              ! have a value which corresponds to the sechiba time step but 
              ! then in stomate.f90 resp_hetero_litter is multiplied by dt. 
              ! Perhaps it could be simplified.
              resp_hetero_litter(ipts,ivm) = resp_hetero_litter(ipts,ivm) + ( 1. - frac_soil(istructural,itarget,ilev) ) * &
                   qd_icarbon * ( 1. - lignin_struc(ipts,ivm,ilev) ) / dt
              decomp_residual = decomp_residual - ( 1. - frac_soil(istructural,itarget,ilev) ) * &
                   qd_icarbon * ( 1. - lignin_struc(ipts,ivm,ilev) ) 
      
              !! 5.1.4 lignin fraction of structural litter goes into
              !!       slow carbon pool + respiration
              som_input(ipts,islow,ivm,icarbon) = som_input(ipts,islow,ivm,icarbon) + frac_soil(istructural,islow,ilev) * &
                   qd_icarbon * lignin_struc(ipts,ivm,ilev) / dt 
              decomp_residual = decomp_residual - frac_soil(istructural,islow,ilev) * &
                   qd_icarbon * lignin_struc(ipts,ivm,ilev)  
              
              ! BE CAREFUL: Here resp_hetero_litter is divided by dt to have a
              ! value which corresponds to the sechiba time step but then in 
              ! stomate.f90 resp_hetero_litter is multiplied by dt. Perhaps it 
              ! could be simplified. Moreover, we must totally adapt the 
              ! routines to the dt_sechiba/one_day time step and avoid some 
              ! constructions that could create bug during future developments.
              ! The last component of the heterothropic respiration could be
              ! calculated as below but the residual will be used to avoid
              ! mass balance problems.
              !resp_hetero_litter(ipts,ivm) = resp_hetero_litter(ipts,ivm) + &
              !     ( 1. - frac_soil(istructural,islow,ilev) ) * qd_icarbon * &
              !     lignin_struc(ipts,ivm,ilev) / dt
              resp_hetero_litter(ipts,ivm) = resp_hetero_litter(ipts,ivm) + (decomp_residual / dt)

              ! LG
              IF ( ( printlev_loc >= 4 ).AND.( ivm == test_pft ).AND.( ipts == test_grid ) ) THEN
                 WRITE(numout,*) 'het_resp before metab decomp', resp_hetero_litter(test_grid,test_pft)
              ENDIF

              !! 5.2 metabolic litter goes into active carbon pool + respiration
              !! 5.2.1 total quantity of metabolic litter that is decomposed
              IF (perma_peat .AND. is_peat(ivm) .AND. soiltile(ipts,istm) .GT. zero) THEN
                 fd = dt*litter_dec_fac(imetabolic) * control_temp(ipts,ilev) * &
                      control_moist_peat(ipts,ilev)
              ELSEIF (agri_peat .AND. (ivm==15 .OR. ivm==16)) THEN
                 fd =  dt*litter_dec_fac(imetabolic) * control_temp(ipts,ilev) * control_moist_croppeat(ipts,ilev)*flux_tot_coeff(1)
              ELSE
                 fd = dt*litter_dec_fac(imetabolic) * control_temp(ipts,ilev) * &
                      control_moist(ipts,ilev)
              ENDIF

              ! Ensure that fd is between 0 and 1
              fd = MAX(zero,MIN(fd,un))

              ! Calculate qd, make sure it has a positive value
              ! Calculate the quantity of metabolic litter that should
              ! be decomposed
              qd_icarbon = litter(ipts,imetabolic,ivm,ilev,icarbon) * fd

              ! Accumulate qd(initrogen) for later use
              qd_initrogen = litter(ipts,imetabolic,ivm,ilev,initrogen) * fd

              ! Keep track of how many errors there are
              err_min_3 = MIN(err_min_3, MIN(litter(ipts,imetabolic,ivm,ilev,icarbon), &
                   litter(ipts,imetabolic,ivm,ilev,initrogen))) 

              ! Subtract the mineralization from the litter pool but make
              ! sure that no negative values for the litter pool are created

              ! qd is calculated as litter(:,imetabolic,ivm,ilev,iele)*fd(:)
              ! as long as fd is positive (which was enforced above) the
              ! subtraction below should never result in a number that is
              ! less than zero (except for precision issues). fd is the same
              ! for carbon and nitrogen and therefore the ratio is preserved
              ! during mineralization. The test conditions should be
              ! satisfied for C and N at the same time.
              err_min_4 = MIN(err_min_4,MIN((litter(ipts,imetabolic,ivm,ilev,icarbon)  -qd_icarbon),&
                   (litter(ipts,imetabolic,ivm,ilev,initrogen)-qd_initrogen)))
              litter(ipts,imetabolic,ivm,ilev,icarbon) = MAX(zero, &
                   litter(ipts,imetabolic,ivm,ilev,icarbon) - qd_icarbon)
              litter(ipts,imetabolic,ivm,ilev,initrogen) = MAX(zero, &
                   litter(ipts,imetabolic,ivm,ilev,initrogen) - qd_initrogen)

              
              !! 5.2.2 decompose same fraction of metabolic part of dead leaves.
              !  to avoid a multiple (for ibelow and iabove) modification of 
              !  dead_leaves, we do this test to do this calcul only once in 
              !  1,nlev loop
              IF (ilev == iabove) THEN
                      dead_leaves(ipts,ivm,imetabolic) = &
                              dead_leaves(ipts,ivm,imetabolic) * ( 1. - fd )

              ENDIF

              ! The commentary below is not true anymore (fungivores were
              ! removed) :
              ! "After large-scale dieback events, so much soil mineral N becomes immobilized 
              ! to decompose litter that too little N is left for plant regrowth. To address 
              ! this, we implicitly represent the action of fungivores, which release N for 
              ! the plants and increase N turnover rates.  We set aside a fraction of qd, 
              ! n_fungivores, which becomes available for plant uptake in nitrogen_dynamics."      
              n_mineralisation(ipts,ivm) = n_mineralisation(ipts,ivm) + qd_initrogen

              !! 5.2.3 put decomposed litter into active 
              !  (or surface) pool + respiration
              decomp_residual = qd_icarbon
              som_input(ipts,itarget,ivm,icarbon) = som_input(ipts,itarget,ivm,icarbon) + &
                   frac_soil(imetabolic,itarget,ilev) * qd_icarbon / dt 
              decomp_residual = decomp_residual - frac_soil(imetabolic,itarget,ilev) * qd_icarbon 
              
              ! BE CAREFUL: Here resp_hetero_litter is divided by dt to have a 
              ! value which corresponds to the sechiba time step but then in 
              ! stomate.f90 resp_hetero_litter is multiplied by dt. 
              ! Perhaps it could be simplified. Moreover, we must totally adapt 
              ! the routines to the dtradia/one_day time step and avoid some 
              ! constructions that could create bugs during future developments.
              ! The last component of the heterothropic respiration could be
              ! calculated as below but the residual will be used to avoid
              ! mass balance problems.
              !resp_hetero_litter(:,ivm) = resp_hetero_litter(:,ivm) + &
              !     ( 1. - frac_soil(imetabolic,itarget,ilev) ) * qd(:,icarbon) / dt
              resp_hetero_litter(ipts,ivm) = resp_hetero_litter(ipts,ivm) + &
                   (decomp_residual / dt)
              
              ! LG
              IF ( ( printlev_loc >= 4 ).AND.( ivm == test_pft ) .AND. &
                   ( ipts == test_grid ) ) THEN
                 WRITE(numout,*) 'het_resp before wood decomp', &
                      resp_hetero_litter(test_grid,test_pft)
              ENDIF

              !! 5.3 woody litter: goes into active and slow carbon 
              !  pools + respiration
           
              !! 5.3.1 total quantity of woody litter which is decomposed
              fd = dt*litter_dec_fac(iwoody) * &
                   control_temp(ipts,ilev) * control_moist(ipts,ilev) * &
                   EXP( -3. * lignin_wood(ipts,ivm,ilev) )
           
              ! Ensure that fd is between 0 and 1
              fd = MAX(zero,MIN(fd,un))

              ! Calculate qd for woody litter

              ! Calculate the quantity of woody litter that should
              ! be decomposed
              qd_icarbon = litter(ipts,iwoody,ivm,ilev,icarbon) * fd
              qd_initrogen = litter(ipts,iwoody,ivm,ilev,initrogen) * fd

              ! Keep track of how many errors there are
              err_min_5 = MIN(err_min_5, MIN(litter(ipts,iwoody,ivm,ilev,icarbon), &
                   litter(ipts,iwoody,ivm,ilev,initrogen)))

              ! Subtract the mineralization from the litter pool but make
              ! sure that no negative values for the litter pool are created

              ! qd is calculated as litter(:,iwoody,ivm,ilev,iele)*fd(:)
              ! as long as fd is positive (which was enforced above) the
              ! subtraction below should never result in a number that is
              ! less than zero (except for precision issues). fd is the same
              ! for carbon and nitrogen and therefore the ratio is preserved
              ! during mineralization. The test conditions should be
              ! satisfied for C and N at the same time.
              err_min_6 = MIN(err_min_6,MIN((litter(ipts,iwoody,ivm,ilev,icarbon) - &
                   qd_icarbon),&
                   (litter(ipts,iwoody,ivm,ilev,initrogen)-qd_initrogen)))
              litter(ipts,iwoody,ivm,ilev,icarbon) = MAX(zero, &
                   litter(ipts,iwoody,ivm,ilev,icarbon) - qd_icarbon)
              litter(ipts,iwoody,ivm,ilev,initrogen) = MAX(zero, &
                   litter(ipts,iwoody,ivm,ilev,initrogen) - qd_initrogen)
              n_mineralisation(ipts,ivm) = n_mineralisation(ipts,ivm) + qd_initrogen
              
              !! 5.3.2 non-lignin fraction of woody litter goes into
              !!       active/structural carbon pool + respiration (per time unit)
              decomp_residual = qd_icarbon
              som_input(ipts,itarget,ivm,icarbon) = som_input(ipts,itarget,ivm,icarbon) + &
                   frac_soil(iwoody,itarget,ilev) * &
                   qd_icarbon * ( 1. - lignin_wood(ipts,ivm,ilev) ) / dt 
              decomp_residual = decomp_residual - frac_soil(iwoody,itarget,ilev) * &
                   qd_icarbon * ( 1. - lignin_wood(ipts,ivm,ilev) )  
              
              resp_hetero_litter(ipts,ivm) = resp_hetero_litter(ipts,ivm) + &
                   ( 1. - frac_soil(iwoody,itarget,ilev) ) * &
                   qd_icarbon * ( 1. - lignin_wood(ipts,ivm,ilev) ) / dt
              decomp_residual = decomp_residual - &
                   ( 1. - frac_soil(iwoody,itarget,ilev) ) * qd_icarbon * &
                   ( 1. - lignin_wood(ipts,ivm,ilev) )
    
              !! 5.3.3 lignin fraction of woody litter goes into
              !!       slow carbon pool + respiration (per time unit)
              som_input(ipts,islow,ivm,icarbon) = som_input(ipts,islow,ivm,icarbon) + frac_soil(iwoody,islow,ilev) * &
                   qd_icarbon * lignin_wood(ipts,ivm,ilev) / dt 
              decomp_residual = decomp_residual - frac_soil(iwoody,islow,ilev) * qd_icarbon * lignin_wood(ipts,ivm,ilev)  

              ! The last component of the heterothropic respiration could be
              ! calculated as below but the residual will be used to avoid
              ! mass balance problems.
              ! resp_hetero_litter(ipts,ivm) = resp_hetero_litter(ipts,ivm) + &
              !         ( 1. - frac_soil(iwoody,islow,ilev) ) * qd_icarbon * &
              !         lignin_wood(ipts,ivm,ilev) / dt
              resp_hetero_litter(ipts,ivm) = resp_hetero_litter(ipts,ivm) + (decomp_residual / dt)
              
              !! 5.4 snag litter: goes into active and slow carbon 
              !  pools + respiration
           
              !! 5.4.1 total quantity of snag litter which is decomposed
              fd = dt*litter_dec_fac(isnag) * &
                   control_temp(ipts,ilev) * control_moist(ipts,ilev) * &
                   EXP( -3. * lignin_snag(ipts,ivm,ilev) )
           
              ! Ensure that fd is between 0 and 1
              fd = MAX(zero,MIN(fd,un))

              ! Calculate qd for snag litter

              ! Calculate the quantity of snag litter that should
              ! be decomposed
              qd_icarbon = litter(ipts,isnag,ivm,ilev,icarbon) * fd
              qd_initrogen = litter(ipts,isnag,ivm,ilev,initrogen) * fd

              ! Keep track of how many errors there are
              err_min_7 = MIN(err_min_7, MIN(litter(ipts,isnag,ivm,ilev,icarbon), &
                   litter(ipts,isnag,ivm,ilev,initrogen)))
              err_min_8 = MIN(err_min_8,MIN((litter(ipts,isnag,ivm,ilev,icarbon)  -qd_icarbon),&
                   (litter(ipts,isnag,ivm,ilev,initrogen)-qd_initrogen)))
              litter(ipts,isnag,ivm,ilev,icarbon) = MAX(zero, &
                   litter(ipts,isnag,ivm,ilev,icarbon) - qd_icarbon)
              litter(ipts,isnag,ivm,ilev,initrogen) = MAX(zero, &
                   litter(ipts,isnag,ivm,ilev,initrogen) - qd_initrogen)
              n_mineralisation(ipts,ivm) = n_mineralisation(ipts,ivm) + qd_initrogen
              
              !! 5.4.2 non-lignin fraction of snag litter goes into
              !!       active/structural carbon pool + respiration (per time unit)
              decomp_residual = qd_icarbon
              som_input(ipts,itarget,ivm,icarbon) = som_input(ipts,itarget,ivm,icarbon) + frac_soil(isnag,itarget,ilev) * &
                                                qd_icarbon * ( 1. - lignin_snag(ipts,ivm,ilev) ) / dt 
              decomp_residual = decomp_residual - frac_soil(isnag,itarget,ilev) * &
                            qd_icarbon * ( 1. - lignin_snag(ipts,ivm,ilev) ) 
              
              resp_hetero_litter(ipts,ivm) = resp_hetero_litter(ipts,ivm) + ( 1. - frac_soil(isnag,itarget,ilev) ) * &
                                         qd_icarbon * ( 1. - lignin_snag(ipts,ivm,ilev) ) / dt
              decomp_residual = decomp_residual - ( 1. - frac_soil(isnag,itarget,ilev) ) * &
                            qd_icarbon * ( 1. - lignin_snag(ipts,ivm,ilev) )
    
              !! 5.4.3 lignin fraction of snag litter goes into
              !!       slow carbon pool + respiration (per time unit)
              som_input(ipts,islow,ivm,icarbon) = som_input(ipts,islow,ivm,icarbon) + frac_soil(isnag,islow,ilev) * &
                                              qd_icarbon * lignin_snag(ipts,ivm,ilev) / dt 
              decomp_residual = decomp_residual - frac_soil(isnag,islow,ilev) * &
                            qd_icarbon * lignin_snag(ipts,ivm,ilev) 

              ! The last component of the heterothropic respiration could be
              ! calculated as below but the residual will be used to avoid
              ! mass balance problems.
              !resp_hetero_litter(ipts,ivm) = resp_hetero_litter(ipts,ivm) + &         
              !                           ( 1. - frac_soil(isnag,islow,ilev) ) * qd_icarbon * lignin_snag(ipts,ivm,ilev) / dt
              resp_hetero_litter(ipts,ivm) = resp_hetero_litter(ipts,ivm) + (decomp_residual / dt)         

              ! LG
              IF ( ( printlev_loc >= 4 ).AND.( ivm == test_pft ).AND.( ipts == test_grid ) ) THEN
                 WRITE(numout,*) 'het_resp after wood decomp', resp_hetero_litter(test_grid,test_pft)
              ENDIF
           ENDDO ! nvm
        ENDDO
     ENDDO ! nlevs

     ! Check for errors
     ! Guillaume M. -- threshold is -min_stomate, not exact zero: emptying an age-class
     ! slot leaves a rounding residue (~1e-12 gC/m2) that a .LT. zero test reads as a
     ! negative pool. err_min_* hold the MINIMUM pool value, not a count, despite the
     ! message below. A real negative pool is orders of magnitude above the threshold.
     IF (err_min_1 .LT. -min_stomate) THEN
        WRITE(numout,*) 'Number of errors found: ,', err_min_1
        WRITE(numout,*) 'Negative structural litter pools in stomate_litter'
        WRITE(numout,*) 'That does not make any sense'
        CALL ipslerr_p(plev,'stomate_litter','Negative litter pools',&
                'This should not be the case','')
     END IF
     IF (err_min_2 .LT. -min_stomate) THEN
        WRITE(numout,*) 'Number of errors found: ,', err_min_2
        WRITE(numout,*) 'Mineralization exceeds litter pools in stomate_litter'
        WRITE(numout,*) 'This will cause problems later on'
        CALL ipslerr_p(plev,'stomate_litter','Structural mineralization exceeds litter pool',&
                'This will cause problems later on','Fix the problem now')
     ENDIF     
     IF (err_min_3 .LT. -min_stomate) THEN
        WRITE(numout,*) 'Number of errors found in icarbon: ,', err_min_3
        WRITE(numout,*) 'Negative litter pools in stomate_litter'
        WRITE(numout,*) 'That does not make any sense'
        CALL ipslerr_p(plev,'stomate_litter','Negative metabolic litter pools',&
                'This should not be the case','')
     END IF
     IF (err_min_4 .LT. -min_stomate) THEN
        WRITE(numout,*) 'Number of errors found: ,', err_min_4
        WRITE(numout,*) 'Mineralization exceeds litter pools in stomate_litter'
        WRITE(numout,*) 'This will cause problems latter on'
        CALL ipslerr_p(plev,'stomate_litter','Metabolic mineralization exceeds litter pool',&
                'This will cause problems later on','Fix the problem now')
     ENDIF
     IF (err_min_5 .LT. -min_stomate) THEN
        WRITE(numout,*) 'Number of errors found: ,', err_min_5
        WRITE(numout,*) 'Negative litter pools in stomate_litter'
        WRITE(numout,*) 'That does not make any sense'
        CALL ipslerr_p(plev,'stomate_litter','Negative woody litter pools',&
                'This should not be the case','')
     END IF
     IF (err_min_6 .LT. -min_stomate) THEN
        WRITE(numout,*) 'Number of errors found: ,', err_min_6
        WRITE(numout,*) 'Mineralization exceeds litter pools in stomate_litter'
        WRITE(numout,*) 'This will cause problems latter on'
        CALL ipslerr_p(plev,'stomate_litter','Woody mineralization exceeds litter pool',&
                'This will cause problems later on','Fix the problem now')
     ENDIF
     IF (err_min_7 .LT. -min_stomate) THEN
        WRITE(numout,*) 'Number of errors found: ,', err_min_7
        WRITE(numout,*) 'Negative litter pools in stomate_litter'
        WRITE(numout,*) 'That does not make any sense'
        CALL ipslerr_p(plev,'stomate_litter','Negative snag litter pools',&
                'This should not be the case','')
     END IF
     IF (err_min_8 .LT. -min_stomate) THEN
        WRITE(numout,*) 'Number of errors found: ,', err_min_8
        WRITE(numout,*) 'Mineralization exceeds litter pools in stomate_litter'
        WRITE(numout,*) 'This will cause problems latter on'
        CALL ipslerr_p(plev,'stomate_litter','Snag mineralization exceeds litter pool',&
                'This will cause problems later on','Fix the problem now')
     ENDIF
     ! 5.4 Calculate som_input for nitrogen
     ! This far the som input for icarbon has been recorded. Now
     ! use the CN ratios to calculate the som input for nitrogen.
     ! Note the difference in units between ::som_input (gN m-2 s-1) 
     ! and ::n_mineralisation, (gN m-2).
     DO ivm = 1, nvm
        DO ipts = 1, npts
           som_input(ipts,iactive,ivm,initrogen) = &
                   som_input(ipts,iactive,ivm,icarbon)/CN_target(ipts,ivm,iactive)  
           som_input(ipts,islow,ivm,initrogen) = &
                   som_input(ipts,islow,ivm,icarbon)/CN_target(ipts,ivm,islow) 
           som_input(ipts,isurface,ivm,initrogen) = &
                   som_input(ipts,isurface,ivm,icarbon)/CN_target(ipts,ivm,isurface)
           ! Ideally we take the nitrogen contained in som_input from the
           ! n_mineralisation pool. Before we can do so we have to check whether
           ! there is enough nitrogen to do so. If not we will have to decompose
           ! less litter than we would like to. That is what is being checked
           ! in the rest of this subroutine. If we have to adjust som_input, we 
           ! may need the following factors later:
           n_min_old(ipts,ivm) = n_mineralisation(ipts,ivm)
           som_input_total(ipts,ivm,initrogen) = &
                   som_input(ipts,iactive,ivm,initrogen) + &
                   som_input(ipts,islow,ivm,initrogen) + &
                   som_input(ipts,isurface,ivm,initrogen)
        ENDDO
     ENDDO

     IF (.NOT.test_adjust_rhetero) THEN

        ! Follow the ORC3 approach
        ! Calculate the required mineralisation, possible mass issues
        ! will be dealt with in nitrogen_dynamics.
        ! Multiply by dt !!
        n_mineralisation(:,:) = n_mineralisation(:,:) - & 
             ( som_input(:,iactive,:,initrogen) + & 
             som_input(:,islow,:,initrogen)   + & 
             som_input(:,isurface,:,initrogen))*dt 

     ELSE

        ! Follow a mass balance approach. A pool should not be negative, in
        ! case n_mineralisation. So litter decomposition will be adjusted
        ! such that n_mineralisation never becomes negative.
        
        ! At this point we calculated the fluxes that will occur if we use
        ! the som_input calculated above in stomate_litter. We will now check
        ! whether this will result in conflicts. If it results in conflicts
        ! som_input will be adjusted. if it does not result in conflicts we
        ! keep som_input as calculated above.
        err_sum_1 = zero
        err_sum_2 = zero
        err_sum_3 = zero
        err_sum_4 = zero
        err_sum_5 = zero
        err_sum_6 = zero
        err_sum_7 = zero
        err_sum_8 = zero
        err_min_1 = zero
        err_max_1 = zero
        ! Calculate n_mineralisation for som_input
        DO ivm = 1, nvm
           DO ipts = 1, npts
              ! No redistrubtion needed unless specified in the code
              ld_redistribute(ipts,ivm) = zero
              delta_hetero_litter(ipts,ivm) = zero
              ! Only check for existing PFTs
              IF (veget_max(ipts,ivm).LE.zero) CYCLE
              ! To decompose the som + som_input from the litter we 
              ! need to take nitrogen from the mineralised pool. Note that 
              ! litter decomposition does not contribute directly to the 
              ! passive component of the soil organic matter [gN m-2 dt-1] 
              ! multiply with dt to obtain [gN m-2].
              som_input_total(ipts,ivm,:) = (som_input(ipts,iactive,ivm,:) + & 
                   som_input(ipts,islow,ivm,:) + & 
                   som_input(ipts,isurface,ivm,:))*dt

              ! Store the positive component of n_mineralisation calculated
              ! above. This n_mineralisation is consistent with som_input.
              ! This is the mineralisation we need to decompose the amount
              ! of litter initialily estimated in stomate_litter
              n_mineralisation_init = n_mineralisation(ipts,ivm)

              ! Calculate the n_mineralisation pool after all internal
              ! n fluxes from som_input(initrogen) are accounted for. This
              ! looks like a worst case scenario. 
              ! n_mineralisation can be both positive or negative
              n_mineralisation(ipts,ivm) = n_mineralisation_init - &
                   som_input_total(ipts,ivm,initrogen)

              ! Check whether the decomposition estimates above are realistic in 
              ! terms of the nitrogen available for immobilisation if 
              ! n_mineralisation is negative.
              IF ((soil_n_min(ipts,ivm,initrate)+soil_n_min(ipts,ivm,iammonium)) + &
                   n_mineralisation_init.LT.zero) THEN

                 err_sum_1 = err_sum_1 + un

                 !+++CHECK+++
                 ! Possible solution but it has never been checked
                 ! All the C and N contained in som_input will need to be moved 
                 ! back into the litter pool. Furthermore the n_mineralisation
                 ! needs to be decreased as is the heterotrophic respiration
                 DO iele = 1,nelements
                    litter_total_new(ipts,ivm,iele) = &
                         SUM(litter(ipts,:,ivm,:,iele)) + som_input_total(ipts,ivm,iele)
                    som_input_total(ipts,ivm,iele) = zero
                 ENDDO

                 ! Adjust resp_hetero_litter proportionaly.
                 delta_hetero_litter(ipts,ivm) = (resp_hetero_litter(ipts,ivm)- &
                      (resp_hetero_litter(ipts,ivm)/n_mineralisation_init) * &
                      (soil_n_min(ipts,ivm,initrate)+soil_n_min(ipts,ivm,iammonium)))*dt

                 ! Update the litter pools
                 litter_total_new(ipts,ivm,initrogen) = litter_total_new(ipts,ivm,initrogen) + &
                      n_mineralisation_init - (soil_n_min(ipts,ivm,initrate)+soil_n_min(ipts,ivm,iammonium))
                 litter_total_new(ipts,ivm,icarbon) = litter_total_new(ipts,ivm,icarbon) + delta_hetero_litter(ipts,ivm)

                 ! Set flag to redistribute the new pools
                 ld_redistribute(ipts,ivm) = 1.
                 !+++++++++++

              ELSEIF ((((soil_n_min(ipts,ivm,initrate)+soil_n_min(ipts,ivm,iammonium)) + &
                   n_mineralisation(ipts,ivm).LT.zero) .AND. &
                   ((soil_n_min(ipts,ivm,initrate)+soil_n_min(ipts,ivm,iammonium)) + &
                   n_mineralisation_init.GT.zero)))THEN

                 ! If the above conditions are satified there exists a solution.
                 ! The new n_mineralisation should get a value between its initial
                 ! and its current value. 

                 ! Initialize
                 litter_total_new(ipts,ivm,:) = zero

                 ! If we have persistent negative values, however, we 
                 ! can reach a situation where we completely deplete our soil_n_min 
                 ! pools to satisfy the demand for immobilization. To remedy this, 
                 ! we'll truncate som_input based on the current 
                 ! n_mineralisation pool, but we'll limit the occurrence of negative
                 ! values. Instead we'll set a minimum (min_n) for n_mineralisation 
                 ! (e.g. .00001 g/m^2) so we can keep some N in the soil. 
                 IF(n_min_old(ipts,ivm).GT.min_n) THEN  

                    ! Later in this code n_mineralisation will be calculated as 
                    ! n_mineralisation(ipts,ivm) = n_min_old(ipts,ivm) - som_input_new(ipts,ivm,initrogen)
                    ! So, if the n_min_old is GT then min_n we will truncate the 
                    ! n_mineralisation at min_n. This requires that som_input_new is 
                    ! calculated as below.
                    som_input_new(ipts,ivm,initrogen) = n_min_old(ipts,ivm) - &
                         min_n

                 ELSE

                    ! If n_mineralisation is less than min_n, take a fraction of the
                    ! remaining n_mineralisation for som_input. This is arbitrarily 
                    ! set to 50% for now. We are just trying to avoid getting values
                    ! that equal to zero.
                    som_input_new(ipts,ivm,initrogen) = 0.5*n_min_old(ipts,ivm)

                 ENDIF

                 ! Error checking
                 IF (som_input_new(ipts,ivm,initrogen).LT.zero) THEN

                    err_min_1 = err_min_1 + un

                 ENDIF

                 ! The adjusted som_input carbon [gC m-2] should be
                 ! som_input_new(ipts,ivm,icarbon) = som_input_new(ipts,ivm,initrogen) / CN_input_total
                 ! However, When too little N is available to reach the target C:N ratio of
                 ! the receiving pool, the equation above results in values that are too 
                 ! low for som_input_total(:,icarbon), meaning carbon litter can not 
                 ! decompose. som_input_total(:,icarbon) still needs to be constrainted, 
                 ! however, or else the C:N ratio in the surface soil pool will become 
                 ! unrealistic.  
                 som_input_new(ipts,ivm,icarbon) = som_input_total(ipts,ivm,icarbon)

                 ! Protect against dividing by zero.  som_input_new(ipts,ivm,initrogen) may be zero here in
                 ! cases where n_mineralisation is zero at the beginning of the run, and then there
                 ! is some new litter input which forces us to this point in the code.  This can
                 ! happen if there is no litter for a long time, but to date we have only seen
                 ! this in cases where the GPP is zero for a long time.  So I will flag it as a 
                 ! warning.
                 IF(som_input_new(ipts,ivm,initrogen) .NE. zero)THEN
                    IF((som_input_new(ipts,ivm,icarbon)/som_input_new(ipts,ivm,initrogen)) .GT. max_cn) THEN
                       som_input_new(ipts,ivm,icarbon) = som_input_new(ipts,ivm,initrogen) * max_cn
                    ENDIF
                 ELSE
                    ! leave the new soil organic matter carbon input as it was, but 
                    ! write a warning. This warning is not very useful for crops
                    ! so it is only written for trees.
                    IF (natural(ivm)) THEN

                       err_sum_2 = err_sum_2 + un

                    ENDIF !natural
                 ENDIF

                 ! Calculate the n_mineralisation pool after all internal
                 ! n fluxes from som_input(initrogen) are accounted for. 
                 ! n_mineralisation can be both positive or negative
                 n_mineralisation(ipts,ivm) = n_min_old(ipts,ivm) - &
                      som_input_new(ipts,ivm,initrogen)

                 ! The nitrogen that was not mineralized should be added back 
                 ! into the litter pool. The carbon no longer included in 
                 ! som_input should also be added back into the litter pool.
                 litter_total_new(ipts,ivm,initrogen) = total_init_nitrogen(ipts,ivm) - &
                      n_mineralisation(ipts,ivm) - som_input_new(ipts,ivm,initrogen)

                 litter_total_new(ipts,ivm,icarbon) = SUM(litter(ipts,:,ivm,:,icarbon)) + &
                      som_input_total(ipts,ivm,icarbon) - som_input_new(ipts,ivm,icarbon)

                 ! Som_input was adjusted to ensure that after accounting for
                 ! n_mineralisation there is enough soil nitrogen left to
                 ! account for nitrogen immobilisation. When litter is 
                 ! transformed into som_input, part of the carbon is lost as
                 ! heterotrophic respiration. For each unit of som_input a
                 ! certain number given by ::resp_ratio is respired
                 IF ((som_input_total(ipts,ivm,icarbon) .GT. zero) .AND. &
                      (resp_hetero_litter(ipts,ivm).GT.zero)) THEN

                    ! The change will be used to adjust the litter pool. The units should
                    ! be [gC m-2]
                    delta_hetero_litter(ipts,ivm) = &
                         (som_input_total(ipts,ivm,icarbon) - som_input_new(ipts,ivm,icarbon)) * &
                         (resp_hetero_litter(ipts,ivm)/som_input_total(ipts,ivm,icarbon))*dt

                    ! Respiration is a flux, the units should be [gC m-2 dt-1]
                    resp_hetero_litter(ipts,ivm) = som_input_new(ipts,ivm,icarbon) * &
                         resp_hetero_litter(ipts,ivm) / som_input_total(ipts,ivm,icarbon)

                 ELSE

                    err_sum_3 = err_sum_3 + un

                 ENDIF

                 ! The carbon that is not respired stays in the litter pool.
                 litter_total_new(ipts,ivm,icarbon) = litter_total_new(ipts,ivm,icarbon) + &
                      delta_hetero_litter(ipts,ivm)

                 ! Set flag to redistribute the new pools
                 ld_redistribute(ipts,ivm) = 1.


              ELSEIF ((soil_n_min(ipts,ivm,initrate)+soil_n_min(ipts,ivm,iammonium)) + &
                   n_mineralisation(ipts,ivm).GE.zero) THEN

                 ! Set flag to redistribute the new pools
                 ld_redistribute(ipts,ivm) = 0

              ELSE

                 err_sum_4 = err_sum_4 + un

              ENDIF
           ENDDO
        ENDDO

        DO iele = 1,nelements 
           DO ivm = 1, nvm  
              DO ipts = 1, npts 
                 IF (ld_redistribute(ipts,ivm) == 1.) THEN
                    ! Share of each component to the total flux
                    IF (som_input_total(ipts,ivm,iele).GT.zero) THEN

                       share_som_input_active = som_input(ipts,iactive,ivm,iele) / &
                            (som_input_total(ipts,ivm,iele) / dt)
                       share_som_input_slow = som_input(ipts,islow,ivm,iele) / &
                            (som_input_total(ipts,ivm,iele) / dt)
                       share_som_input_surface = MAX(1 - share_som_input_active - &
                            share_som_input_slow, zero)
                       ! Truncate som_input values
                       som_input(ipts,iactive,ivm,iele) = (share_som_input_active * &
                            som_input_new(ipts,ivm,iele))
                       som_input(ipts,islow,ivm,iele) = (share_som_input_slow * &
                            som_input_new(ipts,ivm,iele))
                       som_input(ipts,isurface,ivm,iele) = (share_som_input_surface * &
                            som_input_new(ipts,ivm,iele))
                       ! Error checking
                       IF(err_act.GT.1)THEN

                          IF ( abs(share_som_input_surface + share_som_input_active + &
                               share_som_input_slow - 1) .GT. 10*EPSILON(un)) THEN

                             err_max_1 = err_max_1 + un

                          ENDIF
                       ENDIF ! err_act.GT.1

                    ELSE

                       err_sum_5 = err_sum_5 + un

                    ENDIF

                    ! Convert to the initial units to avoid problems when closing the mass balance
                    som_input(ipts,:,ivm,iele) = som_input(ipts,:,ivm,iele) / dt

                    ! Distribute delta_som_input back over the litter pools
                    ! Calculate the share of each litter pool.
                    ! Calculate the share of each litter pool. In principle all values should be
                    ! larger than zero so the max is to prevent against zero stored as very small
                    ! negative numbers, e.g. -1e-17
                    share_litter_struc_a = MAX(litter(ipts,istructural,ivm,iabove,iele) / &
                         SUM(litter(ipts,:,ivm,:,iele)),zero)
                    share_litter_struc_b = MAX(litter(ipts,istructural,ivm,ibelow,iele) / &
                         SUM(litter(ipts,:,ivm,:,iele)),zero)
                    share_litter_woody_a = MAX(litter(ipts,iwoody,ivm,iabove,iele) / &
                         SUM(litter(ipts,:,ivm,:,iele)),zero)
                    share_litter_woody_b = MAX(litter(ipts,iwoody,ivm,ibelow,iele) / &
                         SUM(litter(ipts,:,ivm,:,iele)),zero)
                    share_litter_snag_a = MAX(litter(ipts,isnag,ivm,iabove,iele) / &
                         SUM(litter(ipts,:,ivm,:,iele)),zero)
                    share_litter_snag_b = MAX(litter(ipts,isnag,ivm,ibelow,iele) / &
                         SUM(litter(ipts,:,ivm,:,iele)),zero)

                    ! All PFTs should have metabolic litter. Use the metabolic litter to ensure
                    ! mass balance closure. Not all PFTs have woody litter so the woody litter 
                    ! cannot be used as the residual term. Note that the initial calculation of
                    ! share_litter_met_b as the residual of the other components could result
                    ! in a small negative number. 
                    share_litter_met_a = MAX(litter(ipts,imetabolic,ivm,iabove,iele) / &
                         SUM(litter(ipts,:,ivm,:,iele)),zero)
                    share_litter_met_b = un - (share_litter_struc_a + &
                         share_litter_struc_b + share_litter_met_a + &
                         share_litter_snag_b + share_litter_snag_a + &
                         share_litter_woody_a + share_litter_woody_b)

                    ! Correct precision errors caused by very small negative values for share_litter_met_b
                    IF (share_litter_met_b.LT.zero) THEN

                       IF (share_litter_met_b.LT.-min_stomate) THEN
                          err_sum_6 = err_sum_6 + un
                       ELSE

                          ! Value between 0 and -min_stomate. Correct it. The aboveground wood
                          ! pool is expected to be larger than the belowground pool so transfer
                          ! the precision issue to the aboverground pool.
                          err_sum_7 = err_sum_7 + un

                          !-

                          ! If share_litter_met_a _ share_litter_met_b is positive it is easy to solve the 
                          ! the precision problem
                          IF(share_litter_met_a+share_litter_met_b.GT.zero)THEN

                             ! Move the precision error into share_litter_met_a
                             share_litter_met_a = share_litter_met_a + &
                                  share_litter_met_b
                             share_litter_met_b = zero

                          ELSE

                             ! Difficult case. Nevertheless mass needs to be conserved. Set the negative value
                             ! to zero and then try to rescale all other share_litter_xxxs
                             share_litter_met_b = zero
                             total = MAX(zero, share_litter_struc_a + &
                                  share_litter_struc_b + share_litter_woody_a + &
                                  share_litter_snag_b + share_litter_snag_a + &
                                  share_litter_woody_b + share_litter_met_a)

                             IF(total.GT.zero)THEN

                                ! Rescaling is possible
                                share_litter_struc_a=share_litter_struc_a/total
                                share_litter_struc_b=share_litter_struc_b/total 
                                share_litter_woody_a=share_litter_woody_a/total 
                                share_litter_woody_b=share_litter_woody_b/total
                                share_litter_snag_a=share_litter_snag_a/total 
                                share_litter_snag_b=share_litter_snag_b/total
                                share_litter_met_a=share_litter_met_a/total
                             ELSE

                                err_sum_8 = err_sum_8 + un

                             END IF ! total.GT.zero

                          END IF ! share_litter_met_a+share_litter_met_b.EQ.zero

                       END IF ! share_litter_met_b.LT.-min_stomate

                    END IF ! share_litter_met_b.LT.zero

                    ! Calculate the new litter value. When being decomposed, these values will not result 
                    ! in negative soil_n_min values that cannot be accounted for through immobilisation.
                    litter(ipts,istructural,ivm,iabove,iele) = share_litter_struc_a * &
                         litter_total_new(ipts,ivm,iele)
                    litter(ipts,istructural,ivm,ibelow,iele) = share_litter_struc_b * &
                         litter_total_new(ipts,ivm,iele)
                    litter(ipts,imetabolic,ivm,iabove,iele) = share_litter_met_a * &
                         litter_total_new(ipts,ivm,iele) 
                    litter(ipts,imetabolic,ivm,ibelow,iele) = share_litter_met_b * &
                         litter_total_new(ipts,ivm,iele)
                    litter(ipts,iwoody,ivm,iabove,iele) = share_litter_woody_a * &
                         litter_total_new(ipts,ivm,iele)
                    litter(ipts,iwoody,ivm,ibelow,iele) = share_litter_woody_b * &
                         litter_total_new(ipts,ivm,iele)
                    litter(ipts,isnag,ivm,iabove,iele) = share_litter_snag_a * &
                         litter_total_new(ipts,ivm,iele)
                    litter(ipts,isnag,ivm,ibelow,iele) = share_litter_snag_b * &
                         litter_total_new(ipts,ivm,iele)
                    ! below snag should be empty thus share eq 0
                    IF ( (.NOT. ok_snags) .AND. &  
                         ( (share_litter_snag_a .GT. min_stomate) .OR. &
                         (share_litter_snag_b .GT. min_stomate) ) ) THEN
                       WRITE(numout,*) 'littercalc : snag_ab litter share :', &
                            share_litter_snag_a, &
                            'snag_be litter share :', &
                            share_litter_snag_b
                       CALL ipslerr_p(3,'stomate_litter',&
                            'Snags litter fraction is not zero &
                            though the pool is supposed empty','','')
                       share_litter_woody_a = share_litter_woody_a + &
                            share_litter_snag_a
                       share_litter_woody_b = share_litter_woody_b + &
                            share_litter_snag_b
                       share_litter_snag_a = zero
                       share_litter_snag_b = zero
                       litter(ipts,iwoody,ivm,iabove,iele) = share_litter_woody_a * &
                            litter_total_new(ipts,ivm,iele)
                       litter(ipts,iwoody,ivm,ibelow,iele) = share_litter_woody_b * &
                            litter_total_new(ipts,ivm,iele)
                       litter(ipts,isnag,ivm,iabove,iele) = share_litter_snag_a * &
                            litter_total_new(ipts,ivm,iele)
                       litter(ipts,isnag,ivm,ibelow,iele) = share_litter_snag_b * &
                            litter_total_new(ipts,ivm,iele)
                    ENDIF
                 ENDIF !ld_redistribute 
              ENDDO !npts
           ENDDO !nvm
        ENDDO !nelements

        IF (err_sum_1.GT.zero) THEN
           CALL ipslerr_p(3,'stomate_litter',&
                'First time that this problem is encountered',&
                'Debug the code first before continuing','Really!')
        ENDIF
        IF (err_min_1.GT.zero) THEN
           CALL ipslerr_p(3,'stomate_litter',&
                'Seems that the nitrogen can be satisfied',&
                'Should not have ended up in this part of the code','')
        ENDIF
        IF (err_sum_2.GT.zero) THEN
           CALL ipslerr_p(2,'stomate_litter','Perhaps no n_mineralization &
                & pools, but new litter','Previously only seen in cases &
                & where GPP is zero for long term','')
        ENDIF
        IF (err_sum_3.GT.zero) THEN
           WRITE(numout,*) 'ERROR: surprise when recalculating &
                & resp_hetero_litter, '
           CALL ipslerr_p(3,'stomate_litter','problems with the sign',&
                'when calculating delta_hetero_litter','')
        ENDIF
        IF (err_sum_4.GT.zero) THEN
           WRITE(numout,*) 'ERROR: unexpected case'
           CALL ipslerr_p(3,'stomate_litter','Unexpected case',&
                'Check the logic of IF-statement','')
        ENDIF
        IF (err_sum_5.GT.zero) THEN
           WRITE(numout,*) 'ERROR: new problem. som_input = zero &
                & cannot redistribute the overspended nitrogen'
           CALL ipslerr_p(3,'stomate_litter','som_input = zero',&
                'cannot redistribute the overspended N or C','')
        ENDIF
        IF ( err_max_1 .GT. zero) THEN
           WRITE(numout,*) 'ERROR: sum of shares should be 0 '
           CALL ipslerr_p(3,'stomate_litter','sum of shares',&
                'share_som_input_xxx','')
        ENDIF
        IF (err_sum_6.GT.zero) THEN
           CALL ipslerr_p(3,'stomate_litter','share_litter_met_b is negative',&
                'the value is too negative to be a precision error','')
        ENDIF
        IF (err_sum_7.GT.zero) THEN
           CALL ipslerr_p(2,'stomate_litter','share_litter_met_b is slightly negative',&
                'if this happens rarely it is OK','If it happens a lot debug!')
        ENDIF
        IF (err_sum_8.GT.zero) THEN
           CALL ipslerr_p(3,'stomate_litter','all share_litter_xxx_a/b turn out to be zero',&
                'This is very strange and should have resulted in a divide by zero earlier',&
                '')
        ENDIF

     END IF ! test_adjust_rhetero
    
     !! Reset variables
     ! We assume no distinction exisits in decomposition among different fuel
     ! classes.
     ! Therefore we will just re-adjust litter fuel according to new litter pool
     ! after decomposition.
     ! There is no need to add litterfuel in analytic spinup because of the lines
     ! below.
     litterfuel_all = SUM(litterfuel,DIM=4)
     DO ihour = 1,nhour
       WHERE(litterfuel_all .GE. min_stomate)
         litterfuel(:,:,:,ihour,:) = litterfuel(:,:,:,ihour,:)/litterfuel_all * &
                                     litter(:,:,:,iabove,:)
       ELSEWHERE
         litterfuel(:,:,:,ihour,:) = 0.25 * litter(:,:,:,iabove,:)
 
       ENDWHERE
    ENDDO
    
    !! 6. Check mass balance closure

    !! 6.1 Calculate components of the mass balance
    !  The carbon in turnover and bm_to_litter is moved to different pools but
    !  ::turnover and ::bm_to_litter are never updated. They should either be
    !  set to zero or not be included in the calculation of ::pool_end 
    IF (err_act.GT.1) THEN
       DO ivm = 1,nvm
          DO ipts = 1,npts
             pool_end(ipts,ivm,:) = zero
             closure_intern(ipts,ivm,:) = zero
             check_intern(ipts,ivm,:,:) = zero
             pool_end(ipts,ivm,initrogen) = pool_end(ipts,ivm,initrogen) + &
                     n_mineralisation(ipts,ivm) * veget_max(ipts,ivm)
             pool_end(ipts,ivm,initrogen) = pool_end(ipts,ivm,initrogen) + &
                     (soil_n_min(ipts,ivm,iammonium) + soil_n_min(ipts,ivm,initrate))*veget_max(ipts,ivm)
             check_intern(ipts,ivm,iland2atm,icarbon) = -un * resp_hetero_litter(ipts,ivm) * &
                     veget_max(ipts,ivm) * dt
             DO iele = 1,nelements
             ! The litter pool
                DO ilitt = 1,nlitt
                   DO ilev = 1,nlevs
                      pool_end(ipts,ivm,iele) = pool_end(ipts,ivm,iele) + &
                              (litter(ipts,ilitt,ivm,ilev,iele) * veget_max(ipts,ivm))
                   ENDDO
                ENDDO
                pool_end(ipts,ivm,iele) = pool_end(ipts,ivm,iele) + &
                        SUM(SUM(harvest_pool_acc(ipts,ivm,:,iele,:),2))/(area(ipts)*contfrac(ipts))
             ENDDO
          ENDDO
       ENDDO
       
       !! 6.2 Calculate mass balance
       DO iele = 1,nelements
          DO ivm = 1,nvm
             DO ipts = 1,npts
                check_intern(ipts,ivm,ilat2out,iele) = -un * &
                        SUM(som_input(ipts,:,ivm,iele)) * veget_max(ipts,ivm) * dt
                check_intern(ipts,ivm,ipoolchange,iele) = -un * (pool_end(ipts,ivm,iele) - &
                        pool_start(ipts,ivm,iele))
             ENDDO
          ENDDO
       ENDDO
       DO ivm = 1,nvm
          DO ipts = 1,npts
             check_intern(ipts,ivm,iatm2land,initrogen) = &
                     check_intern(ipts,ivm,iatm2land,initrogen) + &
                     n_input(ipts,ivm,month,imanure) * dt * veget_max(ipts,ivm)
             DO imbc = 1,nmbcomp
                DO iele=1,nelements
                   closure_intern(ipts,ivm,iele) = closure_intern(ipts,ivm,iele) + &
                        check_intern(ipts,ivm,imbc,iele)
                ENDDO
             ENDDO
          ENDDO
       ENDDO     
       
    ENDIF ! err_act.GT.1

    IF (err_act.GT.1) THEN
       CALL check_mass_balance("littercalc", closure_intern, npts, pool_end, &
            pool_start, veget_max, 'pft')
    ENDIF ! err_act.GT.1


  !! 7. calculate fraction of total soil covered by dead leaves

    CALL deadleaf (npts, veget_max, dead_leaves, deadleaf_cover)

    !! 8. (Quasi-)Analytical Spin-up : Start filling MatrixA

    IF (spinup_analytic) THEN

       MatrixA(:,:,:,:) = zero
       VectorB(:,:,:) = zero
       
       ! During the spinup PFT1 (bare soil) should still be empty
       ! for the moment the only way to get carbon/nitrogen in PFT1
       ! is after a LCC. LCC should not be used during spinup.
      
       DO ivm = 2,nvm
          
          !- MatrixA : carbon fluxes leaving the litter
          
          MatrixA(:,ivm ,istructural_above,istructural_above)= - dt*litter_dec_fac(istructural) * &
               control_temp(:,iabove) * control_moist(:,iabove) * exp( -litter_struct_coef * lignin_struc(:,ivm ,iabove) )
          
          MatrixA(:,ivm ,istructural_below,istructural_below) = - dt*litter_dec_fac(istructural) * &
               control_temp(:,ibelow) * control_moist(:,ibelow) * exp( -litter_struct_coef * lignin_struc(:,ivm ,ibelow) )
          
          MatrixA(:,ivm ,imetabolic_above,imetabolic_above) = - dt*litter_dec_fac(imetabolic) * & 
               control_temp(:,iabove) * control_moist(:,iabove)
          
          MatrixA(:,ivm ,imetabolic_below,imetabolic_below) = - dt*litter_dec_fac(imetabolic) * & 
               control_temp(:,ibelow) * control_moist(:,ibelow)
          
          ! Flux leaving the woody above litter pool :
          MatrixA(:, ivm , iwoody_above, iwoody_above) = - dt * litter_dec_fac(iwoody) * control_temp(:,iabove) * &
               control_moist(:,iabove) * exp( -3. * lignin_wood(:,ivm ,iabove) )
          
          ! Flux leaving the woody below litter pool :
          MatrixA(:, ivm , iwoody_below, iwoody_below) = - dt * litter_dec_fac(iwoody) * control_temp(:,ibelow) * &
               control_moist(:,ibelow) * exp( -3. * lignin_wood(:,ivm ,ibelow))

          IF ( ok_soil_carbon_discretization ) THEN  
             ! Flux received by the carbon active from the woody above litter pool :
             MatrixA(:, ivm , iactive_pool, iwoody_above) = frac_soil(iwoody, iactive, iabove) * &
                  dt *litter_dec_fac(iwoody) * control_temp(:,iabove) * control_moist(:,iabove) *  &
                  exp( -3. * lignin_wood(:,ivm ,iabove) ) * ( 1. -  lignin_wood(:,ivm ,iabove) )
          ELSE        
             ! Flux received by the carbon surface from the woody above litter pool :
             MatrixA(:, ivm , isurface_pool, iwoody_above) = frac_soil(iwoody, isurface, iabove) * & 
                  dt *litter_dec_fac(iwoody) * control_temp(:,iabove) * control_moist(:,iabove) *  &
                  exp( -3. * lignin_wood(:,ivm ,iabove) ) * ( 1. -  lignin_wood(:,ivm ,iabove) ) 
          ENDIF

          ! Flux received by the carbon active from the woody below litter pool :
          MatrixA(:, ivm , iactive_pool, iwoody_below) = frac_soil(iwoody, iactive, ibelow) * & 
               dt *litter_dec_fac(iwoody) * control_temp(:,ibelow) * control_moist(:,ibelow) * &
               exp( -3. * lignin_wood(:,ivm ,ibelow) ) * ( 1. -  lignin_wood(:,ivm ,ibelow) ) 
          
          ! Flux received by the carbon slow from the woody above litter pool :
          MatrixA(:, ivm , islow_pool, iwoody_above) = frac_soil(iwoody, islow, iabove) * &
               dt *litter_dec_fac(iwoody) * control_temp(:,iabove) * control_moist(:,iabove) * &
               exp( -3. * lignin_wood(:,ivm ,iabove) ) * lignin_wood(:,ivm ,iabove)
             
          ! Flux received by the carbon slow from the woody below litter pool :
          MatrixA(:, ivm , islow_pool, iwoody_below) =  frac_soil(iwoody, islow, ibelow) * & 
               dt *litter_dec_fac(iwoody) * control_temp(:,ibelow) * control_moist(:,ibelow) * &
               exp( -3. * lignin_wood(:,ivm ,ibelow) ) * lignin_wood(:,ivm ,ibelow)

          ! First order snag fall model from Hararuk et al. 2020
          IF ( ok_snags ) THEN 
             ! Flux leaving the snag above litter pool :
             MatrixA(:, ivm , isnag_above, isnag_above) = - dt * litter_dec_fac(isnag) * control_temp(:,iabove) * &
                  control_moist(:,iabove) * exp( -3. * lignin_snag(:,ivm ,iabove) ) - dt * snag_fall_fac
          
             ! Flux leaving the snag below litter pool :
             MatrixA(:, ivm , isnag_below, isnag_below) = - dt * litter_dec_fac(isnag) * control_temp(:,ibelow) * &
                  control_moist(:,ibelow) * exp( -3. * lignin_snag(:,ivm ,ibelow)) - dt * snag_fall_fac
             
             ! Flux received by the wood above litter pool from the snag above litter pool
             MatrixA(:, ivm , iwoody_above, isnag_above) = dt * snag_fall_fac
          
             ! Flux received by the wood below litter pool from the snag below litter pool
             MatrixA(:, ivm , iwoody_below, isnag_below) = dt * snag_fall_fac

          ELSE  
             ! Flux leaving the snag above litter pool :
             MatrixA(:, ivm , isnag_above, isnag_above) = - dt * litter_dec_fac(isnag) * control_temp(:,iabove) * &
                  control_moist(:,iabove) * exp( -3. * lignin_snag(:,ivm ,iabove) )
          
             ! Flux leaving the snag below litter pool :
             MatrixA(:, ivm , isnag_below, isnag_below) = - dt * litter_dec_fac(isnag) * control_temp(:,ibelow) * &
                  control_moist(:,ibelow) * exp( -3. * lignin_snag(:,ivm ,ibelow))
          
          ENDIF

          IF ( ok_soil_carbon_discretization ) THEN
             ! Flux received by the carbon active from the snag above litter pool :
             MatrixA(:, ivm , iactive_pool, isnag_above) = frac_soil(isnag, iactive, iabove) * &
                  dt *litter_dec_fac(isnag) * control_temp(:,iabove) * control_moist(:,iabove) *  &
                  exp( -3. * lignin_snag(:,ivm ,iabove) ) * ( 1. -  lignin_snag(:,ivm ,iabove) )
          ELSE
             ! Flux received by the carbon surface from the snag above litter pool :
             MatrixA(:, ivm , isurface_pool, isnag_above) = frac_soil(isnag, isurface, iabove) * & 
                  dt *litter_dec_fac(isnag) * control_temp(:,iabove) * control_moist(:,iabove) *  &
                  exp( -3. * lignin_snag(:,ivm ,iabove) ) * ( 1. -  lignin_snag(:,ivm ,iabove) ) 
          ENDIF

          ! Flux received by the carbon active from the snag below litter pool :
          MatrixA(:, ivm , iactive_pool, isnag_below) = frac_soil(isnag, iactive, ibelow) * & 
               dt *litter_dec_fac(isnag) * control_temp(:,ibelow) * control_moist(:,ibelow) * &
               exp( -3. * lignin_snag(:,ivm ,ibelow) ) * ( 1. -  lignin_snag(:,ivm ,ibelow) ) 
          
          ! Flux received by the carbon slow from the snag above litter pool :
          MatrixA(:, ivm , islow_pool, isnag_above) = frac_soil(isnag, islow, iabove) * &
               dt *litter_dec_fac(isnag) * control_temp(:,iabove) * control_moist(:,iabove) * &
               exp( -3. * lignin_snag(:,ivm ,iabove) ) * lignin_snag(:,ivm ,iabove)
             
          ! Flux received by the carbon slow from the snag below litter pool :
          MatrixA(:, ivm , islow_pool, isnag_below) =  frac_soil(isnag, islow, ibelow) * & 
               dt *litter_dec_fac(isnag) * control_temp(:,ibelow) * control_moist(:,ibelow) * &
               exp( -3. * lignin_snag(:,ivm ,ibelow) ) * lignin_snag(:,ivm ,ibelow)
                    
          IF ( ok_soil_carbon_discretization ) THEN
             !- MatrixA : carbon fluxes between the litter and the pools (the rest of the matrix is filled in stomate_soilcarbon.f90)
             MatrixA(:,ivm ,iactive_pool,istructural_above) = frac_soil(istructural,iactive,iabove) * &
                  dt*litter_dec_fac(istructural) * control_temp(:,iabove) * control_moist(:,iabove) * &
                  exp( -litter_struct_coef * lignin_struc(:,ivm ,iabove) ) * &
                  ( 1. - lignin_struc(:,ivm ,iabove) )
          ELSE
             !- MatrixA : carbon fluxes between the litter and the pools (the rest of the matrix is filled in stomate_soilcarbon.f90)
             MatrixA(:,ivm ,isurface_pool,istructural_above) = frac_soil(istructural,isurface,iabove) * &
                  dt*litter_dec_fac(istructural) * control_temp(:,iabove) * control_moist(:,iabove) * & 
                  exp( -litter_struct_coef * lignin_struc(:,ivm ,iabove) ) * &
                  ( 1. - lignin_struc(:,ivm ,iabove) ) 
          ENDIF

          MatrixA(:,ivm ,iactive_pool,istructural_below) = frac_soil(istructural,iactive,ibelow) * &
               dt*litter_dec_fac(istructural) * control_temp(:,ibelow) * control_moist(:,ibelow) * & 
               exp( -litter_struct_coef * lignin_struc(:,ivm ,ibelow) ) * &
               ( 1. - lignin_struc(:,ivm ,ibelow) ) 
             
          IF ( ok_soil_carbon_discretization ) THEN
             MatrixA(:,ivm ,iactive_pool,imetabolic_above) =  frac_soil(imetabolic,iactive,iabove) * &
                  dt*litter_dec_fac(imetabolic) * control_temp(:,iabove) * control_moist(:,iabove)
          ELSE
             MatrixA(:,ivm ,isurface_pool,imetabolic_above) =  frac_soil(imetabolic,isurface,iabove) * &
                  dt*litter_dec_fac(imetabolic) * control_temp(:,iabove) * control_moist(:,iabove) 
          ENDIF
       
          MatrixA(:,ivm ,iactive_pool,imetabolic_below) =  frac_soil(imetabolic,iactive,ibelow) * &
               dt*litter_dec_fac(imetabolic) * control_temp(:,ibelow) * control_moist(:,ibelow)          
                    
          MatrixA(:,ivm ,islow_pool,istructural_above) = frac_soil(istructural,islow,iabove) * &
               dt*litter_dec_fac(istructural) * control_temp(:,iabove) * control_moist(:,iabove) * &
               exp( -litter_struct_coef * lignin_struc(:,ivm ,iabove) )* &
               lignin_struc(:,ivm ,iabove)
          
          MatrixA(:,ivm ,islow_pool,istructural_below) = frac_soil(istructural,islow,ibelow) * &
               dt*litter_dec_fac(istructural) * control_temp(:,ibelow) * control_moist(:,ibelow) *  &
               exp( -litter_struct_coef * lignin_struc(:,ivm ,ibelow) )* &
               lignin_struc(:,ivm ,ibelow) 
          
          
          !- VectorB : carbon input -
          
          VectorB(:,ivm ,istructural_above) = litter_inc(:,istructural,ivm ,iabove,icarbon)
          VectorB(:,ivm ,istructural_below) = litter_inc(:,istructural,ivm ,ibelow,icarbon)
          VectorB(:,ivm ,imetabolic_above) = litter_inc(:,imetabolic,ivm ,iabove,icarbon)
          VectorB(:,ivm ,imetabolic_below) = litter_inc(:,imetabolic,ivm ,ibelow,icarbon)
          VectorB(:,ivm ,iwoody_above) = litter_inc(:,iwoody,ivm ,iabove,icarbon)
          VectorB(:,ivm ,iwoody_below) = litter_inc(:,iwoody,ivm ,ibelow,icarbon)
          VectorB(:,ivm ,isnag_above) = litter_inc(:,isnag,ivm ,iabove,icarbon)
          VectorB(:,ivm ,isnag_below) = litter_inc(:,isnag,ivm ,ibelow,icarbon)

          IF (printlev>=4) WRITE(numout,*) 'We filled MatrixA and VectorB' 
          
!!$          ! Error checking
!!$          ! If there is no veget_max, the litter pools should be empty
!!$          error_count(:) = zero
!!$          WHERE ( SUM(SUM(ABS(litter(:,:,ivm,:,initrogen)),3),2) .GT. min_stomate .AND. &
!!$               veget_max(:,ivm) .LT. min_stomate)
!!$             error_count(:) = un
!!$          END WHERE
!!$          IF (SUM(error_count(:)).GT.zero) THEN
!!$             DO ipts = 1,npts
!!$                IF (error_count(ipts).GT.zero) THEN
!!$                   WRITE(numout,*) 'There is litter but no veget_max in pixel and PFT, ', ipts, ivm
!!$                   WRITE(numout,*) 'veget_max, ', veget_max(ipts,ivm)
!!$                   DO ilitt = 1,nlitt
!!$                      DO ilev = 1,nlevs
!!$                         WRITE(numout,*) 'litter C, ', ipts, ilitt, ivm, ilev, litter(ipts,ilitt,ivm,ilev,icarbon)
!!$                         WRITE(numout,*) 'litter N, ', ipts, ilitt, ivm, ilev, litter(ipts,ilitt,ivm,ilev,initrogen)
!!$                      END DO
!!$                   END DO
!!$                   CALL ipslerr_p(3, 'stomate_litter','Problem when preparing for the analytical spinup',&
!!$                     'veget_max indicates that the PFT is not available','however there is litter')
!!$                END IF
!!$             END DO
!!$          END IF
!!$          !-

          WHERE(litter(:,istructural,ivm ,iabove,initrogen) .GT. min_stomate)
             CN_som_litter_longterm(:,ivm ,istructural_above) = &
                  ( CN_som_litter_longterm(:,ivm ,istructural_above) * (tau_CN_longterm-dt) &
                  + litter(:,istructural,ivm ,iabove,icarbon)/litter(:,istructural,ivm ,iabove,initrogen) * dt)/ (tau_CN_longterm)
          ENDWHERE
          
          WHERE(litter(:,istructural,ivm ,ibelow,initrogen) .GT. min_stomate)
             CN_som_litter_longterm(:,ivm ,istructural_below) = &
                  ( CN_som_litter_longterm(:,ivm ,istructural_below) * (tau_CN_longterm-dt) &
                  + litter(:,istructural,ivm ,ibelow,icarbon)/litter(:,istructural,ivm ,ibelow,initrogen) * dt)/ (tau_CN_longterm)
          ENDWHERE

          WHERE(litter(:,imetabolic,ivm ,iabove,initrogen) .GT. min_stomate)
             CN_som_litter_longterm(:,ivm ,imetabolic_above) = &
                  ( CN_som_litter_longterm(:,ivm ,imetabolic_above) * (tau_CN_longterm-dt) &
                  + litter(:,imetabolic,ivm ,iabove,icarbon)/litter(:,imetabolic,ivm ,iabove,initrogen) * dt)/ (tau_CN_longterm)
          ENDWHERE

          WHERE(litter(:,imetabolic,ivm ,ibelow,initrogen) .GT. min_stomate)
             CN_som_litter_longterm(:,ivm ,imetabolic_below) = &
                  ( CN_som_litter_longterm(:,ivm ,imetabolic_below) * (tau_CN_longterm-dt) &
                  + litter(:,imetabolic,ivm ,ibelow,icarbon)/litter(:,imetabolic,ivm ,ibelow,initrogen) * dt)/ (tau_CN_longterm)
          ENDWHERE

          WHERE(litter(:,iwoody,ivm ,iabove,initrogen) .GT. min_stomate)
             CN_som_litter_longterm(:,ivm ,iwoody_above) = &
                  ( CN_som_litter_longterm(:,ivm ,iwoody_above) * (tau_CN_longterm-dt) &
                  + litter(:,iwoody,ivm ,iabove,icarbon)/litter(:,iwoody,ivm ,iabove,initrogen) * dt)/ (tau_CN_longterm)
          ENDWHERE
          
          WHERE(litter(:,iwoody,ivm ,ibelow,initrogen) .GT. min_stomate)
             CN_som_litter_longterm(:,ivm ,iwoody_below) = &
                  ( CN_som_litter_longterm(:,ivm ,iwoody_below) * (tau_CN_longterm-dt) &
                  + litter(:,iwoody,ivm ,ibelow,icarbon)/litter(:,iwoody,ivm ,ibelow,initrogen) * dt)/ (tau_CN_longterm)
          ENDWHERE            

          WHERE(litter(:,isnag,ivm ,iabove,initrogen) .GT. min_stomate)
             CN_som_litter_longterm(:,ivm ,isnag_above) = &
                  ( CN_som_litter_longterm(:,ivm ,isnag_above) * (tau_CN_longterm-dt) &
                  + litter(:,isnag,ivm ,iabove,icarbon)/litter(:,isnag,ivm ,iabove,initrogen) * dt)/ (tau_CN_longterm)
          ENDWHERE
          
          WHERE(litter(:,isnag,ivm ,ibelow,initrogen) .GT. min_stomate)
             CN_som_litter_longterm(:,ivm ,isnag_below) = &
                  ( CN_som_litter_longterm(:,ivm ,isnag_below) * (tau_CN_longterm-dt) &
                  + litter(:,isnag,ivm ,ibelow,icarbon)/litter(:,isnag,ivm ,ibelow,initrogen) * dt)/ (tau_CN_longterm)
          ENDWHERE

       ENDDO ! Loop over # PFTs

          
    ENDIF ! spinup analytic 

    IF (printlev_loc>=4) WRITE(numout,*) 'Leaving littercalc'

  END SUBROUTINE littercalc


!! ==============================================================================================================================\n
!! SUBROUTINE   : deadleaf
!!
!>\BRIEF        This routine calculates the deadleafcover. 
!!
!! DESCRIPTION  : It first calculates the lai corresponding to the dead leaves (LAI) using 
!! the dead leaves carbon content (DL) the specific leaf area (sla) and the 
!! maximal coverage fraction of a PFT (veget_max) using the following equations:
!! \latexonly
!! \input{deadleaf1.tex}
!! \endlatexonly
!! \n
!! Then, the dead leaf cover (DLC) is calculated as following:\n
!! \latexonly
!! \input{deadleaf2.tex}
!! \endlatexonly
!! \n
!! 
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE: ::deadleaf_cover
!! 
!! REFERENCE(S) : None
!!
!! FLOWCHART : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE deadleaf (npts, veget_max, dead_leaves, deadleaf_cover)

  !! 0. Variable and parameter declaration
    
    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                          :: npts           !! Domain size - number of grid pixels (unitless)
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)           :: dead_leaves    !! Dead leaves per ground unit area, per PFT, 
                                                                          !! metabolic and structural  
                                                                          !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std),DIMENSION(:,:),INTENT(in)               :: veget_max      !! PFT "Maximal" coverage fraction of a PFT defined in 
                                                                          !! the input vegetation map 
                                                                          !! @tex $(m^2 m^{-2})$ @endtex 
    
    !! 0.2 Output variables
    
    REAL(r_std), DIMENSION(:), INTENT(out)              :: deadleaf_cover !! Fraction of soil covered by dead leaves over all PFTs
                                                                          !! (0-1, unitless)

    !! 0.3 Modified variables

    !! 0.4 Local variables

    REAL(r_std), DIMENSION(npts)                        :: dead_lai       !! LAI of dead leaves @tex $(m^2 m^{-2})$ @endtex
    INTEGER(i_std)                                      :: ivm,ipts      !! Index (unitless)
    REAL(r_std)                                         :: total_dead_leaves,biomass_to_lai_sca !! imetabolic + istructural
!_ ================================================================================================================================
  !! 1. LAI of dead leaves
  
    IF (printlev>=3) WRITE(numout,*) 'Entering deadleaf'
    
    DO ipts = 1,npts
       dead_lai(ipts) = zero
       DO ivm = 1,nvm !Loop over PFTs
          total_dead_leaves = dead_leaves(ipts,ivm,imetabolic) + dead_leaves(ipts,ivm,istructural)
          IF(sla_dyn) THEN
             biomass_to_lai_sca = log(1.+ ext_coeff_N(ivm) * total_dead_leaves * slainit(ivm))/(ext_coeff_N(ivm))
          ELSE
             biomass_to_lai_sca = total_dead_leaves * sla(ivm)
          ENDIF
          dead_lai(ipts) = dead_lai(ipts) + biomass_to_lai_sca &
                  * veget_max(ipts,ivm)
       ENDDO
       !! 2. fraction of soil covered by dead leaves
       deadleaf_cover(ipts) = un - exp( - 0.5 * dead_lai(ipts) )
    ENDDO

    IF (printlev>=4) WRITE(numout,*) 'Leaving deadleaf'

  END SUBROUTINE deadleaf   

  
!! ==============================================================================================================================\n
!! SUBROUTINE   : control_moist_func_peat
!!
!>\BRIEF        This routine calculates the moist control factor over peatland for litter and soil C decomposition 
!!
!! ==============================================================================================================================\n

  FUNCTION control_moist_func_peat (npts, mc_peat) RESULT (moistfunc_result)

  !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                       :: npts      !! Domain size - number of grid pixel (unitless)
    REAL(r_std), DIMENSION(npts), INTENT(in)         :: mc_peat   !! relative humidity (unitless)
    REAL(r_std),DIMENSION(dim_moyanno_moist)         :: mc        !! used for Moyano et al., 2012, volumetric moisture, 0.01 interval    
    REAL(r_std),DIMENSION(dim_moyanno_moist)         :: pcsr      !! hetero respiration response to moisture accorcing to moyano
    REAL(r_std),DIMENSION(dim_moyanno_moist)         :: sr        !! temporary normalized response of soil heteo respiration to soil moisture by using the max
    REAL(r_std),DIMENSION(dim_moyanno_moist)         :: corgmat   !! normalized respiration response to moisture
    INTEGER(i_std)                                   :: ind       !! the index for the max of sr
    INTEGER(i_std)                                   :: mc_ind    !! index correspind to model soil moisture

    !! 0.2 Output variables

    REAL(r_std), DIMENSION(npts)                     :: moistfunc_result    !! Moisture control factor (0.25-1, unitless)

    !! 0.3 Modified variables
    INTEGER                                          :: ii,jj
    !! 0.4 Local variables

!_
!================================================================================================================================

    DO jj=1,dim_moyanno_moist
       mc(jj)=0.01+0.02*(jj-1)
    ENDDO

    pcsr(:)=0.97509-0.48212*mc(:)+1.83997*(mc(:)**2)-1.56379*(mc(:)**3)+ &
          0.09867*1.2+1.39944*0.05+0.17938*0.3-0.30307*mc(:)*1.2-0.30885*mc(:)*0.3

    DO jj=1,dim_moyanno_moist
       IF (jj==1) THEN
          sr(jj) = pcsr(jj)
       ELSE
          sr(jj)=sr(jj-1)* pcsr(jj)
       ENDIF
    ENDDO
    sr(:)=sr(:)/ MAXVAL(sr)

    corgmat(:)=sr(:)
    ind= MAXLOC(corgmat,1)
    corgmat(1:ind)=corgmat(1:ind)-MINVAL(corgmat(1:ind))
    corgmat(1:ind)=corgmat(1:ind)/MAXVAL(corgmat(1:ind))

    DO ii=1,npts
       mc_ind = MIN(dim_moyanno_moist, MAX(1, INT(mc_peat(ii)/0.02)+1))
       moistfunc_result(ii)= corgmat(mc_ind)
       moistfunc_result(ii)= MIN(un,MAX(zero,moistfunc_result(ii)))
    ENDDO

  END FUNCTION control_moist_func_peat
!_

!! ==============================================================================================================================\n
!! SUBROUTINE   : control_moist_func_man
!!
!>\BRIEF        This routine calculates the moist control factor over mangroves for litter and soil C decomposition
!!
!! ==============================================================================================================================\n


  FUNCTION control_moist_func_man (npts, mc_man) RESULT (moistfunc_result)

  !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                                      :: npts              !! Domain size - number of grid pixel (unitless)
    REAL(r_std), DIMENSION(npts), INTENT(in)                        :: mc_man            !! humidity (unitless) from orchidee model
    REAL(r_std),DIMENSION(dim_moyanno_moist)                        :: mc                !! humidity for Moyano et al., 2012, created with 0.02 interval    
    REAL(r_std),DIMENSION(dim_moyanno_moist)                        :: pcsr              !! response of soil heteo respiration to soil moisture 
    REAL(r_std),DIMENSION(dim_moyanno_moist)                        :: sr                !! normalized response of soil heteo respiration to soil moisture by using the max
    REAL(r_std),DIMENSION(dim_moyanno_moist)                        :: corgmat           !! response of soil heteo respiration to soil moisture by using the min-max
    INTEGER(i_std)                                                  :: ind               !! the index for the max of sr
    INTEGER(i_std)                                                  :: mc_ind            !! index of real soil moisture corresponding to mc array used for moyano

    !! 0.2 Output variables

    REAL(r_std), DIMENSION(npts)                                    :: moistfunc_result    !! Moisture control factor (0.25-1, unitless)

    !! 0.3 Modified variables
    INTEGER                                                         :: ii,jj
    !! 0.4 Local variables

!_
!================================================================================================================================
    DO jj=1,dim_moyanno_moist
       mc(jj)=0.01+0.02*(jj-1)
    ENDDO

    pcsr(:)=0.97509-0.48212*mc(:)+1.83997*(mc(:)**2)-1.56379*(mc(:)**3)+ &
          0.09867*1.2+1.39944*0.05+0.17938*0.3-0.30307*mc(:)*1.2-0.30885*mc(:)*0.3

    DO jj=1,dim_moyanno_moist
       IF (jj==1) THEN
          sr(jj) = pcsr(jj)
       ELSE
          sr(jj)=sr(jj-1)* pcsr(jj)
       ENDIF
    ENDDO
    sr(:)=sr(:)/ MAXVAL(sr)

    corgmat(:)=sr(:)
    ind= MAXLOC(corgmat,1)
    corgmat(1:ind)=corgmat(1:ind)-MINVAL(corgmat(1:ind))
    corgmat(1:ind)=corgmat(1:ind)/MAXVAL(corgmat(1:ind))

    DO ii=1,npts
       mc_ind = MIN(dim_moyanno_moist, MAX(1, INT(mc_man(ii)/0.02)+1))
       moistfunc_result(ii)= corgmat(mc_ind)
       moistfunc_result(ii)= MIN(un,MAX(zero,moistfunc_result(ii)))
    ENDDO

  END FUNCTION control_moist_func_man

END MODULE stomate_litter
