!  ==============================================================================================================================\n
!  MODULE 	: sechiba
! 
!  CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
!  LICENCE      : IPSL (2006)
!  This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        Structures the calculation of atmospheric and hydrological 
!! variables by calling diffuco_main, enerbil_main, hydrol_main,
!! condveg_main and thermosoil_main. Note that sechiba_main
!! calls slowproc_main and thus indirectly calculates the biogeochemical
!! processes as well.
!!
!!\n DESCRIPTION  : :: shumdiag, :: litterhumdiag and :: stempdiag :: ftempdiag are not 
!! saved in the restart file because at the first time step because they 
!! are recalculated. However, they must be saved as they are in slowproc 
!! which is called before the modules which calculate them.
!! 
!! RECENT CHANGE(S): November 2020: It is possible to define soil hydraulic parameters from maps, 
!!                    as needed for the SP-MIP project (Tafasca Salma and Ducharne Agnes). 
!!                    Here, it leads to declare and allocate global variables.
!! 
!! REFERENCE(S) : None
!!   
!! SVN     :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_sechiba/sechiba.f90 $ 
!! $Date: 2026-05-13 15:33:20 +0200 (mer. 13 mai 2026) $
!! $Revision: 9541 $
!! \n
!_ ================================================================================================================================
 
MODULE sechiba
 
  USE ioipsl
  USE xios_orchidee
  
  ! modules used :
  USE constantes
  USE time, ONLY : one_day, dt_sechiba, FirstTsDay
  USE solar, ONLY: solarang_noon
  USE constantes_soil
  USE pft_parameters
  USE dynamic_parameters
  USE grid
  USE structures
  USE diffuco
  USE condveg
  USE enerbil
  USE mleb
  USE hydrol
  USE lake
  USE thermosoil
  USE sechiba_io_p
  USE slowproc
  USE routing_wrapper
  USE ioipsl_para
  USE chemistry
  USE stomate_laieff
  USE sapiens_lcchange,  ONLY : check_veget
  USE sapiens_forestry,  ONLY : sapiens_forestry_xios_initialize
  USE function_library,  ONLY : get_printlev, wood_to_height, wood_to_qmheight, wood_to_dia
  USE explicitsnow, ONLY : explicitsnow_xios_initialize

  IMPLICIT NONE

  PRIVATE
  PUBLIC sechiba_main, sechiba_initialize, sechiba_clear, &
       sechiba_interface_orchidee_inca, sechiba_xios_initialize

  INTEGER(i_std), SAVE                             :: printlev_loc   !! local printlev for this module
!$OMP THREADPRIVATE(printlev_loc)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION (:) :: indexveg       !! indexing array for the 3D fields of vegetation
!$OMP THREADPRIVATE(indexveg)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION (:) :: indexlai       !! indexing array for the 3D fields of vegetation
!$OMP THREADPRIVATE(indexlai)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION (:) :: indexnobio     !! indexing array for the 3D fields of other surfaces (ice,
                                                                     !! lakes, ...)
!$OMP THREADPRIVATE(indexnobio)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION (:) :: indexsoil      !! indexing array for the 3D fields of soil types (kjpindex*nstm)
!$OMP THREADPRIVATE(indexsoil)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION (:) :: indexgrnd      !! indexing array for the 3D ground heat profiles (kjpindex*ngrnd)
!$OMP THREADPRIVATE(indexgrnd)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION (:) :: indexlayer     !! indexing array for the 3D fields of soil layers in CWRR (kjpindex*nslm)
!$OMP THREADPRIVATE(indexlayer)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION (:) :: indexnslm      !! indexing array for the 3D fields of diagnostic soil layers (kjpindex*nslm)
!$OMP THREADPRIVATE(indexnslm)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION (:) :: indexalb       !! indexing array for the 2 fields of albedo
!$OMP THREADPRIVATE(indexalb)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION (:) :: indexsnow      !! indexing array for the 3D fields snow layers
!$OMP THREADPRIVATE(indexsnow)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION (:) :: indexcan       !! indexing array for the level fields of the canopy
!$OMP THREADPRIVATE(indexcan)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: veget          !! Fraction of vegetation type (unitless, 0-1)       
!$OMP THREADPRIVATE(veget)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: veget_max      !! Max. fraction of vegetation type (LAI -> infty, unitless)
!$OMP THREADPRIVATE(veget_max)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: tot_bare_soil  !! Total evaporating bare soil fraction 
!$OMP THREADPRIVATE(tot_bare_soil)
   REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: height        !! Vegetation Height (m)
!$OMP THREADPRIVATE(height)
   REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: height_dom    !! Dominant vegetation height (m)
!$OMP THREADPRIVATE(height_dom)
   REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: dia_dom       !! Dominant tree diameter (m)
!$OMP THREADPRIVATE(dia_dom)
   REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: height_inv    !! Equivalent tree height to obs (m)
!$OMP THREADPRIVATE(height_inv)
   REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: dia_inv       !! Equivalent tree diameter to obs (m)
!$OMP THREADPRIVATE(dia_inv)
   REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: ba_inv        !! Equivalent basal area to obs (m2)
!$OMP THREADPRIVATE(ba_inv)
   REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: ind_inv       !! Equivalent stand density to obs (trees m-1)
!$OMP THREADPRIVATE(ind_inv)     

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: totfrac_nobio  !! Total fraction of continental ice+lakes+cities+...
                                                                     !! (unitless, 0-1)
!$OMP THREADPRIVATE(totfrac_nobio)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: floodout       !! Flow out of floodplains from hydrol
!$OMP THREADPRIVATE(floodout)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: runoff         !! Surface runoff calculated by hydrol
                                                                     !! @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(runoff)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: drainage       !! Deep drainage calculatedd by hydrol
                                                                     !! @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(drainage)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: returnflow     !! Water flow from lakes and swamps which returns to 
                                                                     !! the grid box @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(returnflow)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: reinfiltration !! Routed water which returns into the soil
!$OMP THREADPRIVATE(reinfiltration)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: irrigation     !! Irrigation flux taken from the routing reservoirs and 
                                                                     !! being put into the upper layers of the soil 
                                                                     !! @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(irrigation)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: emis           !! Surface emissivity (unitless)
!$OMP THREADPRIVATE(emis)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: z0h           !! Surface roughness for heat (m)
!$OMP THREADPRIVATE(z0h)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: z0m           !! Surface roughness for momentum (m)
!$OMP THREADPRIVATE(z0m)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: roughheight    !! Effective height for roughness (m)
!$OMP THREADPRIVATE(roughheight)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: reinf_slope      !! slope coefficient (reinfiltration)
!$OMP THREADPRIVATE(reinf_slope)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: reinf_slope_soil !! slope coefficient (reinfiltration) per soil tile
!$OMP THREADPRIVATE(reinf_slope_soil)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)     :: ks              !! Saturated soil conductivity (mm d^{-1})
!$OMP THREADPRIVATE(ks)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)     :: nvan            !! Van Genushten n parameter (unitless)
!$OMP THREADPRIVATE(nvan)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)     :: avan            !! Van Genushten alpha parameter (mm ^{-1})
!$OMP THREADPRIVATE(avan)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)     :: mcr             !! Residual soil moisture (m^{3} m^{-3})
!$OMP THREADPRIVATE(mcr)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)     :: mcs             !! Saturated soil moisture (m^{3} m^{-3})
!$OMP THREADPRIVATE(mcs)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)     :: mcfc            !! Volumetric water content at field capacity (m^{3} m^{-3})
!$OMP THREADPRIVATE(mcfc)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)     :: mcw             !! Volumetric water content at wilting point (m^{3} m^{-3})
!$OMP THREADPRIVATE(mcw)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:,:):: us             !! Water stress index for transpiration 
                                                                       !! (by soil layer and PFT) (0-1, unitless)
!$OMP THREADPRIVATE(us)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: shumdiag       !! Mean relative soil moisture in the different levels used 
                                                                     !! by thermosoil.f90 (unitless, 0-1)
!$OMP THREADPRIVATE(shumdiag)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: shumdiag_perma !! Saturation degree of the soil 
!$OMP THREADPRIVATE(shumdiag_perma)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: k_litt         !! litter cond.
!$OMP THREADPRIVATE(k_litt)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: litterhumdiag  !! Litter dryness factor (unitless, 0-1)
!$OMP THREADPRIVATE(litterhumdiag)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: stempdiag      !! Temperature which controls canopy evolution (K)
!$OMP THREADPRIVATE(stempdiag)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: ftempdiag      !! Temperature over the full soil column for river temperature (K)
!$OMP THREADPRIVATE(ftempdiag)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: qsintveg       !! Water on vegetation due to interception 
                                                                     !! @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(qsintveg)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: resist_intercept     !! Sum of the resistances applied to intercted water evaporation flux (s m^{-1})
!$OMP THREADPRIVATE(resist_intercept)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: frac_evap_intercept  !! Fraction of the gridcell where intercted water evaporation occurs (-)
!$OMP THREADPRIVATE(frac_evap_intercept)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: frac_evap_transp     !! Fraction of the gridcell where transpiration occurs (-)
!$OMP THREADPRIVATE(frac_evap_transp)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: resist_transp        !! Sum of the resistances applied to transpiration flux (s m^{-1})
!$OMP THREADPRIVATE(resist_transp)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: vbeta3pot      !! Potential vegetation resistance
!$OMP THREADPRIVATE(vbeta3pot)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: gsmean         !! Mean stomatal conductance for CO2 (mol m-2 s-1) 
!$OMP THREADPRIVATE(gsmean) 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: cimean         !! mean intercellular CO2 concentration (ppm)
!$OMP THREADPRIVATE(cimean)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: assimtot         !! Total assimilation per pft
!$OMP THREADPRIVATE(assimtot)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: Rdtot_pft         !! Total respiration per pft
!$OMP THREADPRIVATE(Rdtot_pft)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: vevapwet       !! Interception loss over each PFT 
                                                                     !! @tex $(kg m^{-2} days^{-1})$ @endtex
!$OMP THREADPRIVATE(vevapwet)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: transpir       !! Transpiration @tex $(kg m^{-2} days^{-1})$ @endtex
!$OMP THREADPRIVATE(transpir)

! energy budget variables -----------------------------------------------
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:) :: netrad        !! Net radiation (W m^{-2})
!$OMP THREADPRIVATE(netrad)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:) :: lwabs         !! LW radiation absorbed by the surface (W m^{-2})
!$OMP THREADPRIVATE(lwabs)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:) :: lwnet         !! Net Long-wave radiation (W m^{-2})
!$OMP THREADPRIVATE(lwnet)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:) :: fluxsubli     !! Energy of sublimation (mm day^{-1})
!$OMP THREADPRIVATE(fluxsubli)


! hydraulic stress variables -------------------------------------

   REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: transpir_supply  !! Supply of water for transpiration @tex$$(mm dt^{-1})$ @endtex
!$OMP THREADPRIVATE(transpir_supply)
   REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: vir_transpir_supply !! Supply of water for transpiration @tex$$(mm dt^{-1})$ @endtex
!$OMP THREADPRIVATE(vir_transpir_supply)
   REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: stressed      !! Adjusted ecosystem functioning. Takes the unit of the variable
                                                                     !! used as a proxy for waterstress
!$OMP THREADPRIVATE(stressed)
   REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: unstressed    !! Initial ecosystem functioning after the first calculation and 
                                                                     !! before any recalculations. Takes the unit of the variable used
                                                                     !! as a proxy for unstressed.
!$OMP THREADPRIVATE(unstressed)
   REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: vessel_loss   !! Conductivity lost due to cavitation in the xylem (no unit).
!$OMP THREADPRIVATE(vessel_loss)
   REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) :: e_frac      !! Fraction of water transpired supplied by individual layers (no units)
!$OMP THREADPRIVATE(e_frac)

!$   REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: vir_transpir_mod !! potential transpiration (transpot) divided by veget_max
!!$!$OMP THREADPRIVATE(vir_transpir_mod)
   REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: transpir_mod      !! transpir divided by veget_max
!$OMP THREADPRIVATE(transpir_mod)

   REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:, :, :)    :: transpir_column            !! Supply of water for transpiration 
!$OMP THREADPRIVATE(transpir_column)
   REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:, :, :)    :: transpir_supply_column     !! Supply of water for transpiration 
!$OMP THREADPRIVATE(transpir_supply_column)
   REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:, :, :)    :: transpir_mod_column        !! Supply of water for transpiration 
!$OMP THREADPRIVATE(transpir_mod_column)


   REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: transpot      !! Potential transpiration (needed for irrigation)
!$OMP THREADPRIVATE(transpot)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: qsintmax       !! Maximum amount of water in the canopy interception 
                                                                     !! reservoir @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(qsintmax)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: rveget         !! Surface resistance for the vegetation 
                                                                     !! @tex $(s m^{-1})$ @endtex
!$OMP THREADPRIVATE(rveget)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: rstruct        !! Vegetation structural resistance
!$OMP THREADPRIVATE(rstruct)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: warnings     !! Holds a count of how many warnings we run into
                                                                     !! of different types.  It will get reset at the
                                                                     !! end of each period since we don't wish to restart it.
                                                                     !! Purely a technical diagnostic variable.
!$OMP THREADPRIVATE(warnings)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: snow_nobio_age !! Snow age on non-vegetative surfaces (days)
!$OMP THREADPRIVATE(snow_nobio_age)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: frac_nobio     !! Fraction of non-vegetative surfaces (continental ice, 
                                                                     !! lakes, ...) (unitless, 0-1)
!$OMP THREADPRIVATE(frac_nobio)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:):: assim_param    !! vcmax, nue, and leaf N for photosynthesis  for photosynthesis
!$OMP THREADPRIVATE(assim_param)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: lai            !! Surface foliaire
!$OMP THREADPRIVATE(lai)
  TYPE(laieff_type), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: laieff_fit   !! The parameters for fitting the effective
                                                                     !! LAI function
!$OMP THREADPRIVATE(laieff_fit)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: gpp            !! STOMATE: GPP. gC/m**2 of total area
!$OMP THREADPRIVATE(gpp)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)     :: temp_growth    !! Growth temperature (C) - Is equal to t2m_month 
!$OMP THREADPRIVATE(temp_growth) 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: humrel         !! Relative soil moisture used in stomate to calculate plant water stress (0-1, unitless)
!$OMP THREADPRIVATE(humrel)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: psi_leaf              !! Leaf Water Potential (MPa)
!$OMP THREADPRIVATE(psi_leaf)  
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: psi_leaf_next     !! Approximated Leaf Water Potential at time step n+1 (MPa) (=psi_leaf if no stress)
!$OMP THREADPRIVATE(psi_leaf_next)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: psi_leaf_midday     !! Leaf Water Potential at midday (MPa) 
!$OMP THREADPRIVATE(psi_leaf_midday)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: psi_sto_leaf_save     !! Leaf storage Water Potential at timestep n-1(MPa)
!$OMP THREADPRIVATE(psi_sto_leaf_save) 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: psi_sto_wood_save     !! Wood storage Water Potential at timestep n-1(MPa)
!$OMP THREADPRIVATE(psi_sto_wood_save) 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: psi_root_sup          !! Superficial roots Water Potential (MPa)
!$OMP THREADPRIVATE(psi_root_sup) 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: psi_root_inf          !! Inferior roots Water Potential (MPa)
!$OMP THREADPRIVATE(psi_root_inf)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: psi_soil_sup          !! Superficial soil Water Potential (MPa)
!$OMP THREADPRIVATE(psi_soil_sup)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: psi_soil_inf          !! Inferior soil Water Potential (MPa)
!$OMP THREADPRIVATE(psi_soil_inf)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: psi_xylem_trunk       !! Xylem (trunk level) Water Potential (MPa)
!$OMP THREADPRIVATE(psi_xylem_trunk)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: psi_xylem_leaf        !! Xylem (leaf level) Water Potential (MPa)
!$OMP THREADPRIVATE(psi_xylem_leaf)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: psi_xylem_collar      !! Xylem (collar level) Water Potential (MPa)
!$OMP THREADPRIVATE(psi_xylem_collar)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: psi_sto_wood          !! Wood storage Water Potential (MPa)
!$OMP THREADPRIVATE(psi_sto_wood)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: psi_sto_leaf          !! Leaf storage Water Potential (MPa)
!$OMP THREADPRIVATE(psi_sto_leaf)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:,:)  :: mc_i_sup            !! Water content at each node of the muff of the superficial soil layer (m^3.m^-3)
!$OMP THREADPRIVATE(mc_i_sup) 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:,:)  :: mc_i_inf            !! Water content at each node of the muff of the inferior soil layer (m^3.m^-3)
!$OMP THREADPRIVATE(mc_i_inf)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: F_absorption          !! Total root absorption flux (m^3/s)
!$OMP THREADPRIVATE(F_absorption) 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: vegstress      !! Vegetation moisture stress (only for vegetation growth)
!$OMP THREADPRIVATE(vegstress)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:):: frac_age       !! Age efficacity from STOMATE for isoprene 
!$OMP THREADPRIVATE(frac_age)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: soiltile       !! Fraction of vegtot occupied by each soil tile (0-1, unitless) - so vegtot*soiltile = fraction of each soiltile on gridcell
!$OMP THREADPRIVATE(soiltile)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: fraclut        !! Fraction of each landuse tile (0-1, unitless)
!$OMP THREADPRIVATE(fraclut)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: nwdFraclut     !! Fraction of non-woody vegetation in each landuse tile (0-1, unitless)
!$OMP THREADPRIVATE(nwdFraclut)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION (:) :: njsc           !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless)
!$OMP THREADPRIVATE(njsc)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: frac_sub_snow   !! Fraction of the gridcell where snow sublimation occurs (-)
!$OMP THREADPRIVATE(frac_sub_snow)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: frac_evap_bare  !! Fraction of the gridcell where bare soil evaporation occurs (-)
!$OMP THREADPRIVATE(frac_evap_bare)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: frac_evap_flood !! Fraction of the gridcell where floodplains evaporation occurs (-)
!$OMP THREADPRIVATE(frac_evap_flood)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: soilcap        !!
!$OMP THREADPRIVATE(soilcap)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: soilflx        !!
!$OMP THREADPRIVATE(soilflx)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: temp_sol       !! Surface temperature
!$OMP THREADPRIVATE(temp_sol)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: qsurf          !! near soil air moisture
!$OMP THREADPRIVATE(qsurf)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: flood_res      !! flood reservoir estimate
!$OMP THREADPRIVATE(flood_res)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: flood_frac     !! flooded fraction
!$OMP THREADPRIVATE(flood_frac)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: snow           !! Snow mass [Kg/m^2]
!$OMP THREADPRIVATE(snow)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: snow_age       !! Snow age @tex ($d$) @endtex
!$OMP THREADPRIVATE(snow_age)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: drysoil_frac   !! Fraction of visibly (albedo) Dry soil (Between 0 and 1)
!$OMP THREADPRIVATE(drysoil_frac)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: evap_bare_lim  !! Bare soil stress
!$OMP THREADPRIVATE(evap_bare_lim)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:) ::swc            !! Soil water content (a copy of mc) m3 m-3
!$OMP THREADPRIVATE(swc)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:) ::ksoil          !! Soil conductivity (copy of k mm/d)
!$OMP THREADPRIVATE(ksoil)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: evap_bare_lim_ns !! Bare soil stress
!$OMP THREADPRIVATE(evap_bare_lim_ns)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: co2_to_bm      !! virtual CO2 flux (gC/m**2 of average ground/s)
!$OMP THREADPRIVATE(co2_to_bm)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: evapot         !! Soil Potential Evaporation
!$OMP THREADPRIVATE(evapot)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: evapot_corr    !! Soil Potential Evaporation Correction (Milly 1992)
!$OMP THREADPRIVATE(evapot_corr)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: vevapflo       !! Floodplains evaporation
!$OMP THREADPRIVATE(vevapflo)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: vevapsno       !! Snow evaporation
!$OMP THREADPRIVATE(vevapsno)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: vevapnu        !! Bare soil evaporation
!$OMP THREADPRIVATE(vevapnu)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: tot_melt       !! Total melt
!$OMP THREADPRIVATE(tot_melt)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: resist_equiv_evap     !! Equivalent resistance to the fluxes involved in the evaporation flux (bare soil evaporation + interception loss + transpiration + floodplains evaporation) (s m^{-1})
!$OMP THREADPRIVATE(resist_equiv_evap)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: resist_equiv_evap_veg !! Equivalent resistance to the fluxes involved in the evaporation flux from vegetation (bare soil evaporation + interception loss + transpiration) (s m^{-1})
!$OMP THREADPRIVATE(resist_equiv_evap_veg)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: resist_equiv_tot      !! Equivalent resistance to the fluxes involved in the total evaporative flux from vegetation (all evaporations + snow sublimation) (s m^{-1})
!$OMP THREADPRIVATE(resist_equiv_tot)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: raero          !! Aerodynamic resistance
!$OMP THREADPRIVATE(raero)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: rau            !! Density
!$OMP THREADPRIVATE(rau)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: deadleaf_cover !! Fraction of soil covered by dead leaves
!$OMP THREADPRIVATE(deadleaf_cover)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: ptnlev1        !! 1st level Different levels soil temperature
!$OMP THREADPRIVATE(ptnlev1)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:):: laieff_isotrop !! Effective LAI
!$OMP THREADPRIVATE(laieff_isotrop)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:):: Light_Abs_Tot  !! Absorbed radiation per level for photosynthesis
!$OMP THREADPRIVATE(Light_Abs_Tot)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:):: Light_Tran_Tot !! Transmitted radiation per level for photosynthesis
!$OMP THREADPRIVATE(Light_Tran_Tot)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: Light_Abs_Tot_mean  !! Absorbed radiation per level for photosynthesis
!$OMP THREADPRIVATE(Light_Abs_Tot_mean)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: Light_Alb_Tot_mean !! Refelcted radiation per level for photosynthesis
!$OMP THREADPRIVATE(Light_Alb_Tot_mean)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: mc_layh        !! Volumetric soil moisture for each layer in hydrol(liquid + ice) (m3/m3)
!$OMP THREADPRIVATE(mc_layh)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: mcl_layh       !! Volumetric soil moisture for each layer in hydrol(liquid) (m3/m3)
!$OMP THREADPRIVATE(mcl_layh)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: soilmoist      !! Total soil moisture content for each layer in hydrol(liquid + ice) (mm)
!$OMP THREADPRIVATE(soilmoist)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:):: mc_layh_s      !! Volumetric soil moisture for each layer in hydrol per soiltile (liquid + ice) (m3/m3)
!$OMP THREADPRIVATE(mc_layh_s)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:):: mcl_layh_s     !! Volumetric soil moisture for each layer in hydrol per soiltile (liquid) (m3/m3)
!$OMP THREADPRIVATE(mcl_layh_s)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:):: soilmoist_s    !! Total soil moisture content for each layer in hydrol per soiltile (liquid + ice) (mm)
!$OMP THREADPRIVATE(soilmoist_s)

! multi-layer variables --------------------------------
 REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)      :: u_speed                    !! Canopy wind speed profile
!$OMP THREADPRIVATE(u_speed)
 REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:,:) :: circ_class_biomass        !! Stem diameter @tex $(m)$ @endtex
!$OMP THREADPRIVATE(circ_class_biomass)
 REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  :: circ_class_n                 !! Number of trees within each circumference
!$OMP THREADPRIVATE(circ_class_n)
 REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: loss_gain                   !! loss and gains in each PFT due to LCC
!$OMP THREADPRIVATE(loss_gain)
 REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: veget_max_new               !! New year fraction of vegetation type (0-1, unitless)
!$OMP THREADPRIVATE(veget_max_new)
 REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: frac_nobio_new              !! New year fraction of ice+lakes+cities+... (0-1, unitless)
!$OMP THREADPRIVATE(frac_nobio_new)
 REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)      :: nlevels_loc                  !! The number of physical levels in the canopy
!$OMP THREADPRIVATE(nlevels_loc)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: sum_veget_diff          !! The is the difference between the total
                                                                              !! vegetation fraction in the new and old
                                                                              !! land cover maps.  Useful for land cover
                                                                              !! changes. [-]
!$OMP THREADPRIVATE(sum_veget_diff)
 
REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:,:)      :: z_array_out       !! heights of tree levels from stomate           
!$OMP THREADPRIVATE(z_array_out)

 REAL(r_std),ALLOCATABLE, SAVE, DIMENSION (:, :, :)      :: profile_vbeta3
!$OMP THREADPRIVATE(profile_vbeta3)

 REAL(r_std),ALLOCATABLE, SAVE, DIMENSION (:, :, :)      :: profile_rveget
!$OMP THREADPRIVATE(profile_rveget)

 REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)         :: max_height_store     !! Same as z_array, but one less dimension.
                                                                                 !! @tex $(m)$ @endte
!$OMP THREADPRIVATE(max_height_store)

 REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)         :: delta_c13_assim      !! C13 concentration in delta notation
                                                                                 !! @tex $ permille $ @endtex (per thousand)   
!$OMP THREADPRIVATE(delta_c13_assim)
  
 REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)         :: leaf_ci_out          !! Ci Leaf internal CO2 concentration ()      
!$OMP THREADPRIVATE(leaf_ci_out)
 
 REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)         :: gpp_day              !! Number of time steps when there is gpp
!$OMP THREADPRIVATE(gpp_day)

 REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:)        :: lai_per_level        !! The LAI per vertical level
                                                                                 !! @tex $(m^2 / m^2)$ @endtex
!$OMP THREADPRIVATE(lai_per_level)
 REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)            :: frac_snow_pix        !! The fraction of the whole pixel covered
                                                                                 !! by snow. This is computed from the above
                                                                                 !! two. @tex $-$ @endtex
!$OMP THREADPRIVATE(frac_snow_pix)


  LOGICAL, SAVE                                    :: l_first_sechiba = .TRUE. !! Flag controlling the intialisation (true/false)
!$OMP THREADPRIVATE(l_first_sechiba)

  ! Variables related to snow processes calculations  
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)    :: frac_snow_veg   !! Snow cover fraction on vegetation (unitless)  
!$OMP THREADPRIVATE(frac_snow_veg)  
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)  :: frac_snow_nobio !! Snow cover fraction on continental ice, lakes, etc (unitless)  
!$OMP THREADPRIVATE(frac_snow_nobio) 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)  :: snowrho      !! snow density for each layer (Kg/m^3)
!$OMP THREADPRIVATE(snowrho)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)  :: snowheat     !! snow heat content for each layer (J/m2)
!$OMP THREADPRIVATE(snowheat)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)  :: snowgrain    !! snow grain size (m)
!$OMP THREADPRIVATE(snowgrain)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)  :: snowtemp     !! snow temperature profile (K)
!$OMP THREADPRIVATE(snowtemp)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)  :: snowdz       !! snow layer thickness (m)
!$OMP THREADPRIVATE(snowdz)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)    :: gtemp        !! soil surface temperature
!$OMP THREADPRIVATE(gtemp)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: pgflux       !! net energy into snow pack
!$OMP THREADPRIVATE(pgflux)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)  :: cgrnd_snow   !! Integration coefficient for snow numerical scheme
!$OMP THREADPRIVATE(cgrnd_snow)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)  :: dgrnd_snow   !! Integration coefficient for snow numerical scheme
!$OMP THREADPRIVATE(dgrnd_snow)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)    :: lambda_snow  !! Coefficient of the linear extrapolation of surface temperature 
                                                                  !! from the first and second snow layers
!$OMP THREADPRIVATE(lambda_snow)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)    :: temp_sol_add !! Additional energy to melt snow for snow ablation case (K)
!$OMP THREADPRIVATE(temp_sol_add)

  ! Variables related to ice processes calculations
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)  :: icetemp     !! ice temperature profile (K)
!$OMP THREADPRIVATE(icetemp)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)  :: icedz       !! ice layer thickness (m)
!$OMP THREADPRIVATE(icedz)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)  :: cgrnd_ice   !! Integration coefficient for ice numerical scheme
!$OMP THREADPRIVATE(cgrnd_ice)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)  :: dgrnd_ice   !! Integration coefficient for ice numerical scheme
!$OMP THREADPRIVATE(dgrnd_ice)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)    :: lambda_ice  !! Coefficient of the linear extrapolation of surface temperature 
                                                                  !! from the first and second ice layers
!$OMP THREADPRIVATE(lambda_ice)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION (:) :: ice_sheet_mask !! Ice sheet mask O=no ice, 1=ice unitless)                                                          
!$OMP THREADPRIVATE(ice_sheet_mask)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)    :: qsol_sat_new !! New saturated surface air moisture (kg kg^{-1})
!$OMP THREADPRIVATE(qsol_sat_new)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)    :: qair_new     !! New specific humidity at lowest level (kg kg^{-1})
!$OMP THREADPRIVATE(qair_new)

  !+++CHECK+++
  ! Variables defined in CN-CAN but no longer present in CN
  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:)        :: osfcmelt         !! Indicate snow melting in each gridcell
!$OMP THREADPRIVATE(osfcmelt)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)    :: snowflx         !! Snow flux (W m^{-2})
!$OMP THREADPRIVATE(snowflx)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)    :: snowcap         !! Snow calorific capacity (J K^{-1])
!$OMP THREADPRIVATE(snowcap)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)    :: cgrnd_soil      !! matrix coefficient for the computation of soil, from thermosoil
!$OMP THREADPRIVATE(cgrnd_soil)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)    :: dgrnd_soil      !! matrix coefficient for the computation of soil, from thermosoil
!$OMP THREADPRIVATE(dgrnd_soil)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)    :: zdz1_soil       !! numerical constant from thermosoil
!$OMP THREADPRIVATE(zdz1_soil)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)    :: zdz2_soil       !! numerical constant from thermosoil
!$OMP THREADPRIVATE(zdz2_soil)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: albedo_undersnow !! albedo under the snowpack
!$OMP THREADPRIVATE(albedo_undersnow)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: pkappa_snow     !! snow thermal conductivity
!$OMP THREADPRIVATE(pkappa_snow)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)    :: gthick          !! soil surface layer thickness
!$OMP THREADPRIVATE(gthick)

 !+++++++++++

! Variables related to carbon soil discretization
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:)     :: tdeep            !! Deep temperature profile (K)
!$OMP THREADPRIVATE(tdeep)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:)     :: hsdeep           !! Deep soil humidity profile (unitless)
!$OMP THREADPRIVATE(hsdeep)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:)     :: heat_Zimov       !! Heating associated with decomposition [W/m**3 soil]
!$OMP THREADPRIVATE(heat_Zimov)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)         :: sfluxCH4_deep    !! Surface flux of CH4 to atmosphere from permafrost
!$OMP THREADPRIVATE(sfluxCH4_deep)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)         :: sfluxCO2_deep    !! Surface flux of CO2 to atmosphere from permafrost
!$OMP THREADPRIVATE(sfluxCO2_deep)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:)   :: som_total        !! total soil organic matter for use in thermal calcs (g/m**3)
!$OMP THREADPRIVATE(som_total)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)      :: altmax           !! Maximul active layer thickness (m). Be careful, here active means non frozen.
                                                                           !! Not related with the active soil carbon pool.
!$OMP THREADPRIVATE(altmax) 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)         :: depth_organic_soil !! Depth at which there is still organic matter (m)
!$OMP THREADPRIVATE(depth_organic_soil)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:,:) :: root_profile      !! Normalized root length fraction in each soil layer 
                                                                           !! (0-1, unitless)
!$OMP THREADPRIVATE(root_profile)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)   :: root_depth        !! Node and interface numbers at which the deepest roots
                                                                           !! occur (1 to nslm, unitless) 
!$OMP THREADPRIVATE(root_depth)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)       :: coszang_noon      !! Solar zenith angle at noon (dimensionless) 
!$OMP THREADPRIVATE(coszang_noon)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:)    :: Pgap_cumul        !! The probability of finding a gap in the in canopy from the top 
                                                                           !! of the canopy to a given level (unitless, between 0-1)
!$OMP THREADPRIVATE(Pgap_cumul)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: root_deficit !! water deficit to reach IRRIGATION SOIL MOISTURE TARGET
!$OMP THREADPRIVATE(root_deficit)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: irrig_frac   !! Irrig. fraction interpolated in routing, and saved to pass to slowproc if irrigated_soiltile = .TRUE.
!$OMP THREADPRIVATE(irrig_frac)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: irrigated_next         !! Dynamic irrig. area, calculated in slowproc and passed to routing
                                                                            !!  @tex $(m^{-2})$ @endtex
!$OMP THREADPRIVATE(irrigated_next)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: fraction_aeirrig_sw    !! Fraction of area equipped for irrigation from surface water, of irrig_frac
                                                                            !! 1.0 here corresponds to irrig_frac, not grid cell
!$OMP THREADPRIVATE(fraction_aeirrig_sw)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: precip_longterm        !! Longterm annual precipitation sum mm year^{-1}
!$OMP THREADPRIVATE(precip_longterm)

! peatland and tides
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: wtp                   !! water table position when the peatland is activated
!$OMP THREADPRIVATE(wtp)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: liqwt_ratio           !! Liquid water ratio in the peat, compared to the satured one, when peat is activated
!$OMP THREADPRIVATE(liqwt_ratio)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: mc_peat_above         !! Liquid soil moisture content within the first 4 layers when peat is activated
!$OMP THREADPRIVATE(mc_peat_above)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: mc_croppeat_above     !! Liquid soil moisture content within the first 4 layers when agri_peat is activated
!$OMP THREADPRIVATE(mc_croppeat_above)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: mc_man_above          !! Liquid soil moisture content within the first 4 layers when tides is activated
!$OMP THREADPRIVATE(mc_man_above)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: shumdiag_peat         !! soil moisture content useful for the decomposition, when peat is activated
!$OMP THREADPRIVATE(shumdiag_peat)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: shumdiag_croppeat     !! soil moisture content useful for the decomposition, when agri_peat is activated
!$OMP THREADPRIVATE(shumdiag_croppeat)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: shumdiag_man          !! soil moisture content useful for the decomposition, when tides is activated
!$OMP THREADPRIVATE(shumdiag_man)
  
  ! Variables related to lake energy budget

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:) :: frac_lake       !! Lake cover fraction on continental (unitless)
!$OMP THREADPRIVATE(frac_lake)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:) :: vevapp_lake     !! Lake water evaporation for each tile and grid cell 
                                                                     !! @tex $(kg m^{-2} days^{-1})$ @endtex
!$OMP THREADPRIVATE(vevapp_lake)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:) :: temp_sol_lake   !! Lake temperature (K)
!$OMP THREADPRIVATE(temp_sol_lake)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:) :: fluxsens_lake   !! Sensible heat flux 
                                                                    !! @tex $(W m^{-2})$ @endtex
!$OMP THREADPRIVATE(fluxsens_lake)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:) :: fluxlat_lake    !! Latent heat flux 
                                                                    !! @tex $(W m^{-2})$ @endtex
!$OMP THREADPRIVATE(fluxlat_lake)

CONTAINS


!!  =============================================================================================================================
!! SUBROUTINE:    sechiba_xios_initialize
!!
!>\BRIEF	  Initialize xios dependant defintion before closing context defintion
!!
!! DESCRIPTION:	  Initialize xios dependant defintion before closing context defintion
!!
!! \n
!_ ==============================================================================================================================

  SUBROUTINE sechiba_xios_initialize
    LOGICAL :: lerr
    
    IF (xios_orchidee_ok) THEN
       lerr=xios_orchidee_setvar('min_sechiba',min_sechiba)
       CALL lake_xios_initialize
       CALL slowproc_xios_initialize
       CALL explicitsnow_xios_initialize
       CALL condveg_xios_initialize
       CALL chemistry_xios_initialize
       CALL thermosoil_xios_initialize
       CALL routing_wrapper_xios_initialize
       CALL sapiens_forestry_xios_initialize
    END IF
    IF (printlev_loc>=3) WRITE(numout,*) 'End sechiba_xios_initialize'

  END SUBROUTINE sechiba_xios_initialize




!!  =============================================================================================================================
!! SUBROUTINE:    sechiba_initialize
!!
!>\BRIEF	  Initialize all prinicipal modules by calling their "_initialize" subroutines
!!
!! DESCRIPTION:	  Initialize all prinicipal modules by calling their "_initialize" subroutines
!!
!! \n
!_ ==============================================================================================================================

  SUBROUTINE sechiba_initialize( &
       kjit,         kjpij,        kjpindex,     index,                   &
       lalo,         contfrac,     neighbours,   resolution, zlev,        &
       u,            v,            qair,         temp_air,    &
       petAcoef,     peqAcoef,     petBcoef,     peqBcoef,                &
       precip_rain,  precip_snow,  lwdown,       swnet,      swdown,      &
       pb,           rest_id,      hist_id,      hist2_id,                &
       rest_id_stom, hist_id_stom, hist_id_stom_IPCC,                     &
       coastalflow,  riverflow,    tsol_rad,     vevapp,       qsurf_out, &
       z0m_out,      z0h_out,      albedo,       fluxsens,     fluxlat,      emis_out,  &
       temp_sol_new, tq_cdrag,     coszang)

!! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                               :: kjit              !! Time step number (unitless)
    INTEGER(i_std), INTENT(in)                               :: kjpij             !! Total size of the un-compressed grid 
                                                                                  !! (unitless)
    INTEGER(i_std), INTENT(in)                               :: kjpindex          !! Domain size - terrestrial pixels only 
                                                                                  !! (unitless)
    INTEGER(i_std),INTENT (in)                               :: rest_id           !! _Restart_ file identifier (unitless)
    INTEGER(i_std),INTENT (in)                               :: hist_id           !! _History_ file identifier (unitless)
    INTEGER(i_std),INTENT (in)                               :: hist2_id          !! _History_ file 2 identifier (unitless)
    INTEGER(i_std),INTENT (in)                               :: rest_id_stom      !! STOMATE's _Restart_ file identifier 
                                                                                  !! (unitless)
    INTEGER(i_std),INTENT (in)                               :: hist_id_stom      !! STOMATE's _History_ file identifier 
                                                                                  !! (unitless)
    INTEGER(i_std),INTENT(in)                                :: hist_id_stom_IPCC !! STOMATE's IPCC _history_ file file 
                                                                                  !! identifier (unitless)
    REAL(r_std),DIMENSION (kjpindex,2), INTENT (in)          :: lalo              !! Geographic coordinates (latitude,longitude)
                                                                                  !! for grid cells (degrees)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: contfrac          !! Fraction of continent in the grid 
                                                                                  !! (unitless, 0-1)
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)         :: index             !! Indices of the pixels on the map. 
                                                                                  !! Sechiba uses a reduced grid excluding oceans
                                                                                  !! ::index contains the indices of the 
                                                                                  !! terrestrial pixels only! (unitless)
    INTEGER(i_std), DIMENSION (kjpindex,NbNeighb), INTENT(in):: neighbours        !! Neighboring grid points if land!(unitless)
    REAL(r_std), DIMENSION (kjpindex,2), INTENT(in)          :: resolution        !! Size in x and y of the grid (m)
    
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: u                 !! Lowest level wind speed in direction u 
                                                                                  !! @tex $(m.s^{-1})$ @endtex 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: v                 !! Lowest level wind speed in direction v 
                                                                                  !! @tex $(m.s^{-1})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: zlev              !! Height of first layer (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: qair              !! Lowest level specific humidity 
                                                                                  !! @tex $(kg kg^{-1})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: precip_rain       !! Rain precipitation 
                                                                                  !! @tex $(kg m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: precip_snow       !! Snow precipitation 
                                                                                  !! @tex $(kg m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: lwdown            !! Down-welling long-wave flux 
                                                                                  !! @tex $(W m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: swnet             !! Net surface short-wave flux 
                                                                                  !! @tex $(W m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: swdown            !! Down-welling surface short-wave flux 
                                                                                  !! @tex $(W m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: temp_air          !! Air temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: petAcoef          !! Coefficients A for T from the Planetary 
                                                                                  !! Boundary Layer
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: peqAcoef          !! Coefficients A for q from the Planetary 
                                                                                  !! Boundary Layer
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: petBcoef          !! Coefficients B for T from the Planetary 
                                                                                  !! Boundary Layer
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: peqBcoef          !! Coefficients B for q from the Planetary 
                                                                                  !! Boundary Layer
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: pb                !! Surface pressure (hPa)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: coszang           !! cosine of solar zenith angle


!! 0.2 Output variables
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: coastalflow       !! Outflow on coastal points by small basins.
                                                                                  !! This is the water which flows in a disperse 
                                                                                  !! way into the ocean
                                                                                  !! @tex $(kg dt_routing^{-1})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: riverflow         !! Outflow of the major rivers.
                                                                                  !! The flux will be located on the continental 
                                                                                  !! grid but this should be a coastal point  
                                                                                  !! @tex $(kg dt_routing^{-1})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: tsol_rad          !! Radiative surface temperature 
                                                                                  !! @tex $(W m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: vevapp            !! Total of evaporation 
                                                                                  !! @tex $(kg m^{-2} days^{-1})$ @endtex
    
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: qsurf_out         !! Surface specific humidity 
                                                                                  !! @tex $(kg kg^{-1})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: z0m_out           !! Surface roughness momentum (output diagnostic, m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: z0h_out           !! Surface roughness heat (output diagnostic, m)
    REAL(r_std),DIMENSION (kjpindex,2), INTENT (out)         :: albedo            !! Surface albedo for visible and near-infrared
                                                                                  !! (unitless, 0-1)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: fluxsens          !! Sensible heat flux 
                                                                                  !! @tex $(W m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: fluxlat           !! Latent heat flux 
                                                                                  !! @tex $(W m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: emis_out          !! Emissivity (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: temp_sol_new      !! New surface temperature (K)

!! 0.3 Modified
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)         :: tq_cdrag          !! Surface drag coefficient (-)

!! 0.4 Local variables
    INTEGER(i_std)                                           :: ji, jv, ilev      !! Index (unitless)
    REAL(r_std), DIMENSION(kjpindex)                         :: zmaxh_glo         !! 2D field of constant soil depth (zmaxh) (m)
    CHARACTER(LEN=80)                                        :: var_name          !! To store variables names for I/O (unitless)
    REAL(r_std),DIMENSION (kjpindex)                         :: epot_air          !! Air potential energy (??J)
    !+++CHECK+++
    ! Variable dimensions xxx_Tot_mean are for a single pixel. Causing 1+1 problems
    ! when running a larger domain. Needs to be corrected when implementing
    ! a global use of the multi-layer energy budget.
!!$    REAL(r_std),DIMENSION(nlevels_tot)                       :: Light_Abs_Tot_mean!! total light absorption for a given canopy level
!!$    REAL(r_std),DIMENSION(nlevels_tot)                       :: Light_Alb_Tot_mean!! total albedo for a given level
    !+++++++++++
    INTEGER(i_std)                                           :: init_config       !! Identifer of the configuration used to 
                                                                                  !! initialize stomate/or sechiba
    REAL(r_std), DIMENSION(kjpindex,nslm,nvm)                :: mc_layh_pft       !! As mc_layh per pft
    REAL(r_std), DIMENSION(kjpindex,nslm,nvm)                :: mcl_layh_pft      !! As mcl_layh per pft
    REAL(r_std), DIMENSION(kjpindex,nslm,nvm)                :: soilmoist_pft     !! As soilmoist per pft
!_ ================================================================================================================================

    !! Initialize local printlev
    printlev_loc=get_printlev('sechiba')

    IF (printlev_loc>=3) WRITE(numout,*) 'Start sechiba_initialize'

    !! 1. Initialize variables on first call
    
    !! 1.2 Initialize most of sechiba's variables
    CALL sechiba_init (kjit, kjpij, kjpindex, index, rest_id, lalo)
    
    !! 1.3 Initialize stomate's variables
    CALL slowproc_initialize (kjit,          kjpij,        kjpindex,                          &
                              rest_id,       rest_id_stom, hist_id_stom,   hist_id_stom_IPCC, &
                              index,         indexveg,     lalo,           neighbours,        &
                              resolution,    contfrac,     temp_air,       coszang_noon,      &
                              soiltile,      reinf_slope,  reinf_slope_soil, deadleaf_cover, assim_param,       &
                              ks,             nvan,        avan,           mcr,               &
                              mcs,            mcfc,        mcw,                               &
                              frac_age,      height,       lai,            veget,             &
                              frac_nobio,    njsc,         veget_max,      fraclut,           &
                              nwdfraclut,    tot_bare_soil,totfrac_nobio,  qsintmax,          &
                              temp_growth,   circ_class_biomass,                              &
                              circ_class_n,  lai_per_level,laieff_fit,                        &
                              z_array_out,   max_height_store,                                &
                              som_total,     heat_Zimov,   altmax,         depth_organic_soil,&
                              loss_gain,     veget_max_new,frac_nobio_new, height_dom,        &
                              height_inv,    dia_dom,      dia_inv,        ba_inv,            &
                              ind_inv,       Pgap_cumul,   irrigated_next, irrig_frac,   fraction_aeirrig_sw,   &
                              precip_longterm)
    

    !! 1.4 Initialize diffusion coefficients
    CALL diffuco_initialize (kjit,    kjpindex, index,                  &
                             rest_id, lalo,     neighbours, resolution, &
                             rstruct, tq_cdrag)

    
    !! 1.5 Initialize variables for the energy budget
    IF (ok_mleb) THEN
       !! Multi-layer energy budget
       CALL mleb_initialize (kjit, kjpindex, rest_id,                    &
                             temp_air, qair,                             & 
                             temp_sol, temp_sol_new, tsol_rad,           &
                             evapot,   evapot_corr,  qsurf,    fluxsens, &
                             fluxlat,  vevapp,                           &
                             u_speed, z_array_out,  & 
                             max_height_store )   
    ELSE
       !! Energy budget using enerbil module
       CALL enerbil_initialize (kjit,     kjpindex,     index,    rest_id,  &
                                qair,                                       &
                                temp_sol, temp_sol_new, tsol_rad,           &
                                evapot,   evapot_corr,  qsurf,    fluxsens, &
                                fluxlat,  vevapp )
    END IF


    !! 1.6 Initialize variables for lake energy budget   
    IF (ok_lake_energy) THEN
       CALL lake_initialize (kjit,         kjpindex,      index,         rest_id,       &
                             lalo,         neighbours,    resolution,    contfrac,      &
                             frac_lake,    temp_sol_lake)
    ENDIF

    
    !! 1.7 Initialize remaining hydrological variables

    CALL hydrol_initialize (ks,  nvan, avan, mcr, mcs, mcfc, mcw,    &
         kjit,           kjpindex,  index,         rest_id,          &
         njsc,           soiltile,  veget,         veget_max,        &
         frac_nobio,     altmax,                                     &
         humrel,         vegstress, drysoil_frac,                    &
         shumdiag_perma,    qsintveg,                                &
         evap_bare_lim, evap_bare_lim_ns,  snow,   snow_age,         &
         snow_nobio_age, snowrho,   snowtemp,      snowgrain,        &
         snowdz,                   snowheat,                         &
         mc_layh,        mcl_layh,  soilmoist, mc_layh_s,  mcl_layh_s,   &
         soilmoist_s,    swc,       ksoil,         root_profile,     &
         us,             icetemp,   icedz, ice_sheet_mask,           &
         psi_leaf, psi_leaf_next, psi_leaf_midday, psi_sto_leaf_save, &
         psi_sto_wood_save,         psi_root_sup,  psi_root_inf, &
         psi_soil_sup, psi_soil_inf, psi_xylem_trunk, psi_xylem_leaf, &
         psi_xylem_collar, psi_sto_wood, psi_sto_leaf,&
         mc_i_sup, mc_i_inf, F_absorption, &
         wtp,  liqwt_ratio,   shumdiag_peat,  shumdiag_croppeat, shumdiag_man)
       


    !! 1.9 Initialize surface parameters (emissivity, albedo and roughness)
    !+++CHECK+++
    ! Variable dimensions xxx_Tot_mean are for a single pixel. Causing 1+1 problems
    ! when running a larger domain. Needs to be corrected when implementing
    ! a global use of the multi-layer energy budget.
    CALL condveg_initialize (kjit, kjpindex, index, rest_id, &
         lalo, neighbours, resolution, contfrac, &
         veget, veget_max, frac_nobio, totfrac_nobio, &
         zlev, snow, snow_age, snow_nobio_age, &
         drysoil_frac, height, height_dom, snowdz,snowrho, tot_bare_soil, &
         temp_air, pb, u, v, &
         lai, &
         emis, albedo, z0m, z0h, roughheight, &
         frac_snow_veg,frac_snow_nobio, &
         coszang, &
         Light_Abs_Tot, Light_Tran_Tot, laieff_fit, &
         Light_Abs_Tot_mean, Light_Alb_Tot_mean, &
         laieff_isotrop)
    !+++++++++++
    
    !! 1.10 Initialization of soil thermodynamics
    DO jv = 1,nvm
       mc_layh_pft(:,:,jv) = mc_layh_s(:,:,pref_soil_veg(jv))
       mcl_layh_pft(:,:,jv) = mcl_layh_s(:,:,pref_soil_veg(jv))
       soilmoist_pft(:,:,jv) = soilmoist_s(:,:,pref_soil_veg(jv))
    END DO
    ! 
    CALL thermosoil_initialize (kjit, kjpindex, rest_id,  &
         temp_sol_new, snow,       shumdiag_perma,        &
         soilcap,      soilflx,    depth_organic_soil,    & 
         stempdiag,    ftempdiag,  gtemp,                 &
         mc_layh,  mcl_layh,   soilmoist, njsc,soiltile,   &
         frac_snow_veg,frac_snow_nobio,totfrac_nobio,     &
         snowdz, snowrho, snowtemp, lambda_snow, cgrnd_snow, &
         dgrnd_snow,   lambda_ice,  cgrnd_ice,   dgrnd_ice,  &
         icetemp,      icedz,       ice_sheet_mask,     pb,  &
         som_total, &
         veget_max, mc_layh_pft, mcl_layh_pft, soilmoist_pft, &
         ptnlev1)

    
    !! 1.12 Initialize river routing
    IF ( river_routing .AND. nbp_glo .GT. 1) THEN

       !! 1.12.1 Initialize river routing
       CALL routing_wrapper_initialize( &
            kjit,        kjpindex,       index,                 &
            rest_id,     hist_id,        hist2_id,   lalo,      &
            neighbours,  resolution,     contfrac,   stempdiag, ftempdiag, &
            soiltile,    irrig_frac,     veget_max,  irrigated_next, &
            returnflow,  reinfiltration, irrigation, riverflow, &
            coastalflow, flood_frac,     flood_res )
    ELSE
       !! 1.12.2 No routing, set variables to zero
       riverflow(:) = zero
       coastalflow(:) = zero
       returnflow(:) = zero
       reinfiltration(:) = zero
       irrigation(:) = zero
       flood_frac(:) = zero
       flood_res(:) = zero

       CALL xios_orchidee_send_field("coastalflow",coastalflow/dt_sechiba)  
       CALL xios_orchidee_send_field("riverflow",riverflow/dt_sechiba)  

    ENDIF
    
    !! 1.13 Write internal variables to output fields
    z0m_out(:) = z0m(:)
    z0h_out(:) = z0h(:)
    emis_out(:) = emis(:) 
    qsurf_out(:) = qsurf(:)

    !! 2. Output variables only once
    zmaxh_glo(:) = zmaxh
    CALL xios_orchidee_send_field("zmaxh",zmaxh_glo)

    IF (printlev_loc>=3) WRITE(numout,*) 'sechiba_initialize done'

  END SUBROUTINE sechiba_initialize

!! ==============================================================================================================================\n
!! SUBROUTINE 	: sechiba_main
!!
!>\BRIEF        Main routine for the sechiba module performing three functions:
!! calculating temporal evolution of all variables and preparation of output and 
!! restart files (during the last call only)
!!
!!\n DESCRIPTION : Main routine for the sechiba module. 
!! One time step evolution consists of:
!! - call sechiba_var_init to do some initialization,
!! - call slowproc_main to do some daily calculations
!! - call diffuco_main for diffusion coefficient calculation,
!! - call enerbil_main for energy budget calculation,
!! - call hydrol_main for hydrologic processes calculation,
!! - call condveg_main for surface conditions such as roughness, albedo, and emmisivity,
!! - call thermosoil_main for soil thermodynamic calculation,
!! - call sechiba_end to swap previous to new fields.
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): Hydrological variables (:: coastalflow and :: riverflow),
!! components of the energy budget (:: tsol_rad, :: vevapp, :: fluxsens, 
!! :: temp_sol_new and :: fluxlat), surface characteristics (:: z0_out, :: emis_out, 
!! :: tq_cdrag and :: albedo) and land use related CO2 fluxes :: fco2_flux_out 
!! :: fco2_lu_out, :: fco2_wh_out, ::fco2_ha_out)            
!!
!! REFERENCE(S)	: 
!!
!! FLOWCHART    : 
!! \latexonly 
!! \includegraphics[scale = 0.5]{sechibamainflow.png}
!! \endlatexonly
!! \n
!_ ================================================================================================================================

  SUBROUTINE sechiba_main (kjit, kjpij, kjpindex, index, &
       & ldrestart_read, ldrestart_write, &
       & lalo, contfrac, neighbours, resolution,&
       & zlev, u, v, qair, temp_air, epot_air, ccanopy, &
       & tq_cdrag, petAcoef, peqAcoef, petBcoef, peqBcoef, &
       & precip_rain, precip_snow, lwdown, swnet, swdown, coszang, pb, &
       & vevapp, fluxsens, fluxlat, coastalflow, riverflow, &
       & fco2_flux_out, fco2_lu_out, fco2_wh_out, fco2_ha_out, &
       & tsol_rad, temp_sol_new, qsurf_out, albedo, emis_out, z0m_out, z0h_out,&
       & veget_out, lai_out, height_out, &
       & rest_id, hist_id, hist2_id, rest_id_stom, hist_id_stom, hist_id_stom_IPCC, &
       & run_off_lic, run_off_lic_frac)

!! 0.1 Input variables
    
    INTEGER(i_std), INTENT(in)                               :: kjit              !! Time step number (unitless)
    INTEGER(i_std), INTENT(in)                               :: kjpij             !! Total size of the un-compressed grid 
                                                                                  !! (unitless)
    INTEGER(i_std), INTENT(in)                               :: kjpindex          !! Domain size - terrestrial pixels only 
                                                                                  !! (unitless)
    INTEGER(i_std),INTENT (in)                               :: rest_id           !! _Restart_ file identifier (unitless)
    INTEGER(i_std),INTENT (in)                               :: hist_id           !! _History_ file identifier (unitless)
    INTEGER(i_std),INTENT (in)                               :: hist2_id          !! _History_ file 2 identifier (unitless)
    INTEGER(i_std),INTENT (in)                               :: rest_id_stom      !! STOMATE's _Restart_ file identifier 
                                                                                  !! (unitless)
    INTEGER(i_std),INTENT (in)                               :: hist_id_stom      !! STOMATE's _History_ file identifier 
                                                                                  !! (unitless)
    INTEGER(i_std),INTENT(in)                                :: hist_id_stom_IPCC !! STOMATE's IPCC _history_ file file 
                                                                                  !! identifier (unitless)
    LOGICAL, INTENT(in)                                      :: ldrestart_read    !! Logical for _restart_ file to read 
                                                                                  !! (true/false)
    LOGICAL, INTENT(in)                                      :: ldrestart_write   !! Logical for _restart_ file to write 
                                                                                  !! (true/false)
    REAL(r_std),DIMENSION (kjpindex,2), INTENT (in)          :: lalo              !! Geographic coordinates (latitude,longitude)
                                                                                  !! for grid cells (degrees)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: contfrac          !! Fraction of continent in the grid 
                                                                                  !! (unitless, 0-1)
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)         :: index             !! Indices of the pixels on the map. 
                                                                                  !! Sechiba uses a reduced grid excluding oceans
                                                                                  !! ::index contains the indices of the 
                                                                                  !! terrestrial pixels only! (unitless)
    INTEGER(i_std), DIMENSION(kjpindex,NbNeighb), INTENT(in) :: neighbours        !! Neighboring grid points if land!(unitless)
    REAL(r_std), DIMENSION (kjpindex,2), INTENT(in)          :: resolution        !! Size in x and y of the grid (m)
    
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: u                 !! Lowest level wind speed in direction u 
                                                                                  !! @tex $(m.s^{-1})$ @endtex 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: v                 !! Lowest level wind speed in direction v 
                                                                                  !! @tex $(m.s^{-1})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: zlev              !! Height of first layer (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: qair              !! Lowest level specific humidity 
                                                                                  !! @tex $(kg kg^{-1})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: precip_snow       !! Snow precipitation 
                                                                                  !! @tex $(kg m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: lwdown            !! Down-welling long-wave flux 
                                                                                  !! @tex $(W m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: coszang           !! Cosine of the solar zenith angle (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: swnet             !! Net surface short-wave flux 
                                                                                  !! @tex $(W m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: swdown            !! Down-welling surface short-wave flux 
                                                                                  !! @tex $(W m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: temp_air          !! Air temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: epot_air          !! Air potential energy (??J)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: ccanopy           !! CO2 concentration in the canopy (ppm)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: petAcoef          !! Coefficients A for T from the Planetary 
                                                                                  !! Boundary Layer
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: peqAcoef          !! Coefficients A for q from the Planetary 
                                                                                  !! Boundary Layer
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: petBcoef          !! Coefficients B for T from the Planetary 
                                                                                  !! Boundary Layer
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: peqBcoef          !! Coefficients B for q from the Planetary 
                                                                                  !! Boundary Layer
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: pb                !! Surface pressure (hPa)


!! 0.2 Output variables

    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: coastalflow       !! Outflow on coastal points by small basins.
                                                                                  !! This is the water which flows in a disperse 
                                                                                  !! way into the ocean
                                                                                  !! @tex $(kg dt_routing^{-1})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: riverflow         !! Outflow of the major rivers.
                                                                                  !! The flux will be located on the continental 
                                                                                  !! grid but this should be a coastal point  
                                                                                  !! @tex $(kg dt_routing^{-1})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: tsol_rad          !! Radiative surface temperature 
                                                                                  !! @tex $(W m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: vevapp            !! Total of evaporation 
                                                                                  !! @tex $(kg m^{-2} days^{-1})$ @endtex
    
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: qsurf_out         !! Surface specific humidity 
                                                                                  !! @tex $(kg kg^{-1})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: z0m_out           !! Surface roughness momentum (output diagnostic, m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: z0h_out           !! Surface roughness heat (output diagnostic, m)
    REAL(r_std),DIMENSION (kjpindex,2), INTENT (out)         :: albedo            !! Surface albedo for visible and near-infrared
                                                                                  !! (unitless, 0-1)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: fluxsens          !! Sensible heat flux 
                                                                                  !! @tex $(W m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: fluxlat           !! Latent heat flux 
                                                                                  !! @tex $(W m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: emis_out          !! Emissivity (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: fco2_flux_out     !! Sum CO2 flux over PFTs 
                                                                                  !! (gC/m2/s)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: fco2_lu_out       !! Land Cover Change CO2 flux (gC/m2/one_day)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: fco2_wh_out       !! Wood harvest CO2 flux (gC/m2/one_day)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: fco2_ha_out       !! Crop harvest CO2 flux (gC/m2/one_day)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)       :: veget_out         !! Fraction of vegetation type (unitless, 0-1)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)       :: lai_out           !! Leaf area index (m^2 m^{-2}) 
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)       :: height_out        !! Vegetation Height (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: temp_sol_new      !! New surface temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: run_off_lic       !! Contains calving, melting and liquid precipitation
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: run_off_lic_frac  !! Contains cell fraction corresponding to run_off_lic

!! 0.3 Modified

    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)         :: tq_cdrag          !! Surface drag coefficient  (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)         :: precip_rain       !! Rain precipitation 
                                                                                  !! @tex $(kg m^{-2})$ @endtex

!! 0.4 local variables

    INTEGER(i_std)                                           :: ji, jv, ilevel    !! Index (unitless)
    REAL(r_std), DIMENSION(kjpindex)                         :: histvar           !! Computations for history files (unitless)
    REAL(r_std), DIMENSION(kjpindex,nlut)                    :: histvar2          !! Computations for history files (unitless)
    CHARACTER(LEN=80)                                        :: var_name          !! To store variables names for I/O (unitless)
    REAL(r_std), DIMENSION(kjpindex)                         :: sum_treefrac      !! Total fraction occupied by trees (0-1, uniless) 
    REAL(r_std), DIMENSION(kjpindex)                         :: sum_height_treefrac    !! Sum of product height*treefrac over tree tiles (meter*landfraction) 
    REAL(r_std), DIMENSION(kjpindex)                         :: sum_grassfracC3   !! Total fraction occupied by C3 grasses (0-1, unitless)
    REAL(r_std), DIMENSION(kjpindex)                         :: sum_grassfracC4   !! Total fraction occupied by C4 grasses (0-1, unitless)
    REAL(r_std), DIMENSION(kjpindex)                         :: sum_height_grassfrac    !! Sum of product height*(grassfracC3+grassfracC4) over tree tiles (meter*landfraction) 
    REAL(r_std), DIMENSION(kjpindex)                         :: vegheighttree     !! Mean height for tree land 
    REAL(r_std), DIMENSION(kjpindex)                         :: vegheightgrass    !! Mean height for grass land 

    REAL(r_std), DIMENSION(kjpindex)                         :: sum_cropfracC3    !! Total fraction occupied by C3 crops (0-1, unitess)
    REAL(r_std), DIMENSION(kjpindex)                         :: sum_cropfracC4    !! Total fraction occupied by C4 crops (0-1, unitess)
    REAL(r_std), DIMENSION(kjpindex)                         :: sum_treeFracNdlEvg!! Total fraction occupied by treeFracNdlEvg (0-1, unitess)
    REAL(r_std), DIMENSION(kjpindex)                         :: sum_treeFracBdlEvg!! Total fraction occupied by treeFracBdlEvg (0-1, unitess)
    REAL(r_std), DIMENSION(kjpindex)                         :: sum_treeFracNdlDcd!! Total fraction occupied by treeFracNdlDcd (0-1, unitess)
    REAL(r_std), DIMENSION(kjpindex)                         :: sum_treeFracBdlDcd!! Total fraction occupied by treeFracBdlDcd (0-1, unitess)
    REAL(r_std),DIMENSION (kjpindex,nvm)                     :: tmc_pft           !! Total soil water per PFT (mm/m2) 
    REAL(r_std),DIMENSION (kjpindex,nvm)                     :: drainage_pft      !! Drainage per PFT (mm/m2)   
    REAL(r_std),DIMENSION (kjpindex,nvm)                     :: runoff_pft        !! Runoff per PFT (mm/m2)   
    REAL(r_std),DIMENSION (kjpindex,nvm)                     :: swc_pft           !! Relative Soil water content [tmcr:tmcs] per pft (-)    
    REAL(r_std), DIMENSION(kjpindex)                         :: grndflux          !! Net energy into soil (W/m2)
    REAL(r_std), DIMENSION(kjpindex,nsnow)                   :: snowliq           !! Liquid water content (m)
    REAL(r_std), DIMENSION (kjpindex)                        :: snowmelt          !! Snow melt [mm/dt_sechiba]
!cdc ajout zrainfall pour calcul smb    
    REAL(r_std), DIMENSION (kjpindex)                        :: zrainfall         !! Rain precipitation on snow 
                                                                                  !! @tex $(kg m^{-2})$ @endtex
    !+++CHECK+++
    ! Variable dimensions xxx_Tot_mean are for a single pixel. Causing 1+1 problems
    ! when running a larger domain. Needs to be corrected when implementing
    ! a global use of the multi-layer energy budget.
!!$    REAL(r_std),DIMENSION(nlevels_tot)                       :: Light_Abs_Tot_mean!! total light absorption for a given canopy level
!!$    REAL(r_std),DIMENSION(nlevels_tot)                       :: Light_Alb_Tot_mean!! total albedo for a given level
    !+++++++++++
    REAL(r_std), DIMENSION(kjpindex)                         :: snow_age_diag     !! Only for diag, contains xios_default_val
    REAL(r_std), DIMENSION(kjpindex,nnobio)                  :: snow_nobio_age_diag !! Only for diag, contains xios_default_val
    REAL(r_std), DIMENSION(kjpindex)                         :: snowage_glob      !! Snow age on total area including snow on vegetated and bare soil and nobio area @tex ($d$) @endtex 
    REAL(r_std), DIMENSION(kjpindex,nlut)                    :: gpplut            !! GPP on landuse tile, only for diagnostics
    REAL(r_std), DIMENSION(kjpindex,nslm,nvm)                :: mc_layh_pft       !! As mc_layh per pft
    REAL(r_std), DIMENSION(kjpindex,nslm,nvm)                :: mcl_layh_pft      !! As mcl_layh per pft
    REAL(r_std), DIMENSION(kjpindex,nslm,nvm)                :: soilmoist_pft     !! As soilmoist per pft
    REAL(r_std), DIMENSION(kjpindex)                         :: vbeta2sum         !! sum of vbeta2 coefficients across all PFTs (-)
    REAL(r_std), DIMENSION(kjpindex)                         :: vbeta3sum         !! sum of vbeta3 coefficients across all PFTs (-)
    REAL(r_std), DIMENSION(kjpindex)                         :: netrad            !! Net radiation (W m^{-2})
    REAL(r_std), DIMENSION(kjpindex)                         :: lwabs             !! LW radiation absorbed by the surface (W m^{-2})
    REAL(r_std), DIMENSION(kjpindex)                         :: lwnet             !! Net Long-wave radiation (W m^{-2})
    REAL(r_std), DIMENSION(kjpindex)                         :: fluxsubli         !! Energy of sublimation (mm day^{-1})
    REAL(r_std),DIMENSION (kjpindex,nvm)                     :: vbeta23           !! Beta for fraction of wetted foliage that will
                                                                                  !! transpire once intercepted water has evaporated (-)
    REAL(r_std), DIMENSION (kjpindex,nvm,nlai)               :: leaf_ci           !! intercellular CO2 concentration (ppm)
    REAL(r_std),DIMENSION (kjpindex)                         :: vbeta1            !! Temporary variable : Beta for snow sublimation
    REAL(r_std),DIMENSION (kjpindex,nvm)                     :: vbeta2            !! Temporary variable : Beta for interception loss
    REAL(r_std),DIMENSION (kjpindex,nvm)                     :: vbeta3            !! Temporary variable : Beta for transpiration
    REAL(r_std),DIMENSION (kjpindex)                         :: vbeta4            !! Temporary variable : Beta for bare soil evaporation
    REAL(r_std),DIMENSION (kjpindex)                         :: vbeta5            !! Temporary variable : Beta for floodplains evaporation
    REAL(r_std),DIMENSION (kjpindex)                         :: vbeta             !! Temporary variable : Total beta 
    REAL(r_std), DIMENSION (kjpindex,nvm,nlevels_tot)        :: gs_distribution
    REAL(r_std), DIMENSION (kjpindex,nvm,nlevels_tot)        :: gs_diffuco_output

    REAL(r_std), DIMENSION (kjpindex)                        :: VPD               !! Vapor Pressure Deficit (kPa)
    REAL(r_std), DIMENSION (kjpindex)                        :: air_relhum        !! Air relative humidity (%)
    REAL(r_std), DIMENSION (kjpindex,nvm)                    :: f_gco2_out

    REAL(r_std), DIMENSION (kjpindex,nvm,nlevels_tot)        :: gstot_component
    REAL(r_std), DIMENSION (kjpindex,nvm,nlevels_tot)        :: gstot_frac
    REAL(r_std),DIMENSION (kjpindex,nvm,nlevels_tot)         :: profile_vbeta3
    REAL(r_std),DIMENSION (kjpindex,nvm,nlevels_tot)         :: profile_rveget
    INTEGER                                                  :: ilev,ivm,ipts  
    LOGICAL, DIMENSION(kjpindex,nvm)                         :: failed_vegfrac     !! Pixels and pfts for which the model failed to 
                                                                                   !! find a PFT were some residual fraction could be added (true/false)
    REAL(r_std),DIMENSION (kjpindex,nvm)                     :: cresist            !!coefficient for resistances (??)
    REAL(r_std), DIMENSION(kjpindex)                         :: mcs_hydrol         !! Saturated volumetric water content output to be used in stomate_soilcarbon
    REAL(r_std), DIMENSION(kjpindex)                         :: mcfc_hydrol        !! Volumetric water content at field capacity output to be used in stomate_soilcarbon

    REAL(r_std), DIMENSION (kjpindex,nvm,nlevels_tot+1)      :: info_limitphoto    !!Save information of limitation for photosynthesis
    REAL(r_std), DIMENSION (kjpindex,nvm,nlevels_tot)        :: JJ_out
    REAL(r_std), DIMENSION (kjpindex,nvm,nlevels_tot)        :: assimi_lev
    LOGICAL                                                  :: hydrol_flag        !! flag that 'trips' the energy budget for each grid square,
                                                                                   !! but should be false when mleb is called from sechiba   
    LOGICAL, DIMENSION (kjpindex)                            :: hydrol_flag2       !! flag that 'trips' the energy budget for each grid square,
                                                                                   !! hydrol_flag2 on the full global domain for all processors 
                                                                                   !! but should be false when mleb is called from sechiba 
    LOGICAL, DIMENSION (kjpindex,nvm)                        :: hydrol_flag3       !! flag that 'trips' the energy budget for each grid square and PFT
                                                                                   !! but should be false when mleb is called from sechiba

    REAL(r_std), DIMENSION (kjpindex,nvm)                    :: temporary_array    !! A temporary array
    REAL(r_std), DIMENSION (kjpindex)                        :: zenith_angle
    REAL(r_std), DIMENSION (ncirc)                           :: circ_height
    REAL(r_std), DIMENSION (ncirc)                           :: circ_dia
    REAL(r_std), DIMENSION (ncirc)                           :: circ_n
 
!_ ================================================================================================================================

    IF (printlev_loc>=3) WRITE(numout,*) 'Start sechiba_main kjpindex =',kjpindex
    IF (FirstTsDay) THEN
       CALL solarang_noon(kjpindex, lalo, coszang_noon)
    END IF
    zenith_angle(:) = ACOS(coszang_noon)*(180/pi)
    WHERE (zenith_angle(:).LE.zero)
       zenith_angle(:) = xios_default_val
    ENDWHERE
    CALL xios_orchidee_send_field("zenith_ang_noon",zenith_angle(:))
    CALL xios_orchidee_send_field("nav_lat",lalo(:,1))
    CALL xios_orchidee_send_field("nav_lon",lalo(:,2))
 
    ! Daily initialization
    IF (FirstTsDay) THEN

       ! Dynamic parameters that are calculated once per day
       CALL dynamic_parameters_calculate(kjpindex, precip_longterm)

       ! As dynamic parameter changes, all the variables uses dynamic parameters
       ! needs to be recalculated before other calculations start. (Ticket XXX)
       CALL slowproc_canopy (kjpindex, circ_class_biomass, circ_class_n, &
            veget_max, lai_per_level, z_array_out, &
            max_height_store, laieff_fit, frac_age)
       !! Calculation of veget and soiltile.
       CALL slowproc_veget (kjpindex, lai_per_level, z_array_out, &
            coszang_noon, circ_class_biomass, circ_class_n, frac_nobio, totfrac_nobio, &
            veget_max, veget, soiltile, tot_bare_soil, &
            fraclut, nwdFraclut, Pgap_cumul, precip_longterm)

       IF (read_lai .NE. 1) THEN

          DO jv = 2,nvm
             DO ji = 1,kjpindex
             ! Skip if veget_max = 0, else calculate height and lai
                IF (veget_max(ji,jv) .EQ. zero) THEN
                   height(ji,jv) = zero
                   height_dom(ji,jv) = zero
                   height_inv(ji,jv) = zero
                   dia_dom(ji,jv) = zero
                   ba_inv(ji,jv) = zero
                   ind_inv(ji,jv) = zero
                ELSE
                   height(ji,jv) = wood_to_qmheight(circ_class_biomass(ji,jv,:,:,icarbon), &
                     circ_class_n(ji,jv,:), jv, pipe_tune2(ji,jv))
                   IF(is_tree(jv))THEN
                      ! Calculate diameter and height of the dominant diameter
                      ! class
                      circ_height(:) = wood_to_height(circ_class_biomass(ji,jv,:,:,icarbon), &
                         jv, pipe_tune2(ji,jv))
                      height_dom(ji,jv) = circ_height(ncirc)
                      circ_dia(:) = wood_to_dia(circ_class_biomass(ji,jv,:,:,icarbon),jv, &
                         pipe_tune2(ji,jv))
                      dia_dom(ji,jv) = circ_dia(ncirc)
                      ! Calculate quadratic mean height and diameter of the trees
                      ! of which the diameter exceeds a threshold. This is closely
                      ! matches the height and diameter reported in forest 
                      ! inventories
                      circ_n(:) = circ_class_n(ji,jv,:)
                      WHERE (circ_dia(:) .lt. dia_thresh_inv(jv))
                         circ_dia(:) = zero
                         circ_height(:) = zero
                         circ_n(:) = zero
                      END WHERE
                      IF (SUM(circ_n(:)).EQ.zero) THEN
                         dia_inv(ji,jv) =  zero
                         height_inv(ji,jv) = zero
                         ba_inv(ji,jv) = zero
                         ind_inv(ji,jv) = zero
                      ELSE
                         dia_inv(ji,jv) = (SUM(circ_dia(:)**2*circ_n(:))/SUM(circ_n(:)))**0.5
                         height_inv(ji,jv) = (SUM(circ_height(:)**2*circ_n(:))/SUM(circ_n(:)))**0.5
                         ind_inv(ji,jv) = SUM(circ_n(:))
                         ba_inv(ji,jv) = pi/4*SUM(circ_dia(:)**2*circ_n(:))*m2_to_ha
                      END IF
                   ELSE
                      ! Grasses and crops have a height but no
                      ! diameter in ORCHIDEE
                      height_dom(ji,jv) = height(ji,jv)
                      height_inv(ji,jv) = height(ji,jv)
                      dia_dom(ji,jv) = undef
                      dia_inv(ji,jv) = undef
                      ba_inv(ji,jv) = undef
                      ind_inv(ji,jv) = circ_class_n(ji,jv,1)
                   ENDIF ! IF is tree
                ENDIF
             ENDDO
          ENDDO ! Do nvm
       END IF ! read_lai .EQ. 1

       qsintmax(:,:) = qsintcst * veget(:,:) * lai(:,:)
       qsintmax(:,1) = zero

    ENDIF ! IF FirstTsDay

    !! 1. Initialize variables at each time step
    CALL sechiba_var_init (kjpindex, rau, pb, temp_air) 
   
    !! 1.1 For rain exclusion experiments
    CALL rain_exclusion(precip_rain)

    !! 2. Compute diffusion coefficients
    CALL diffuco_main (kjit, kjpindex, index, indexveg, indexlai, u, v, &
         zlev, z0m, z0h, roughheight, temp_sol, temp_air, temp_growth, rau, tq_cdrag, qsurf, qair, pb, &
         evap_bare_lim, evap_bare_lim_ns, evapot, evapot_corr, snow, flood_frac, flood_res, &
         frac_nobio, totfrac_nobio, &
         swnet, swdown, coszang, ccanopy, humrel, veget, veget_max, lai, &
         qsintveg, qsintmax, assim_param, &
         resist_equiv_evap, resist_equiv_evap_veg, resist_equiv_tot, frac_sub_snow, frac_evap_intercept, resist_intercept, &
         frac_evap_transp, resist_transp, vbeta3pot, frac_evap_bare, frac_evap_flood, raero, &
         gsmean, rveget, rstruct, cimean, gpp, cresist, &
         lalo, neighbours, resolution, ptnlev1, precip_rain, frac_age, tot_bare_soil, frac_snow_veg, frac_snow_nobio, &
         hist_id, hist2_id, vbeta2sum, vbeta3sum, Light_Abs_Tot, &
         Light_Tran_Tot, lai_per_level, vbeta23, leaf_ci, &
         gs_distribution, gs_diffuco_output, gstot_component, gstot_frac,&
         VPD, air_relhum, f_gco2_out, warnings, u_speed, profile_vbeta3, profile_rveget, &
         delta_c13_assim, leaf_ci_out, info_limitphoto, JJ_out, assimi_lev, &
         psi_leaf_next, assimtot, Rdtot_pft)

    IF ( ok_c13 ) THEN 
                 
      ! When there is no photosynthesis, carbon isotopic values is not calulated.  
      ! Because daily mean of carbon isotopic value in ORCHIDEE contained night values,  
      ! to avoide biased daily mean value by zeros, daily mean of carbon isotopic  
      ! discriminations should be calculated only when there is photosynthesis 
      
       gpp_day(:,:) = zero 
      
       DO ipts = 1,kjpindex 
          DO jv = 1,nvm 
             IF (gpp(ipts,jv) .LT. min_sechiba) THEN 
                gpp_day(ipts,jv) = zero 
             ELSE 
                gpp_day(ipts,jv) = gpp_day(ipts,jv) + 1 
             END IF
          END DO
       END DO
       
    END IF ! ok_c13


    !! 3. Compute energy balance

    IF (.NOT. ok_mleb) THEN

       CALL enerbil_main (kjit, kjpindex, netrad, lwabs, lwnet, fluxsubli, index, &
            indexveg, lwdown, swnet, epot_air, temp_air, u, v, petAcoef, petBcoef, qair, &
            peqAcoef, peqBcoef, pb, rau, resist_equiv_evap, resist_equiv_evap_veg, resist_equiv_tot, &
            frac_sub_snow, frac_evap_intercept, resist_intercept, &
            frac_evap_transp, resist_transp, vbeta3pot, frac_evap_bare, frac_evap_flood, raero, &
            emis, soilflx, soilcap, tq_cdrag, humrel, fluxsens, fluxlat, &
            vevapp, transpir, transpot, vevapnu, vevapwet, vevapsno, vevapflo, temp_sol, tsol_rad, &
            temp_sol_new, qsurf, evapot, evapot_corr, precip_rain, pgflux, &
            snowdz, temp_sol_add, qair_new, qsol_sat_new)
       
       ! Write enerbil output. Because enerbil_main is called twice in the case
       ! with water stress when ok_hydrol_arch, calling this module in the 
       ! enrebil_main can cause a problem. For now, the module is placed separately.
       CALL enerbil_write (kjit, kjpindex, hist_id, hist2_id, index, &
            evapot, evapot_corr, fluxsubli, lwdown, lwnet, netrad, transpir, &
            vevapflo, vevapnu, vevapp, vevapsno, vevapwet)
    ELSE
       
       ! These hydrol flag will always be false, when the model is not run with
       ! the hydraulic architecture.
       hydrol_flag = .FALSE.
       hydrol_flag2 = .FALSE.
       hydrol_flag3 = .FALSE.
       
       !! Temporary variables used only in mleb. As soon as mleb is clean, those vbeta will be removed.
       vbeta1(:) = frac_sub_snow(:)
       vbeta4(:) = frac_evap_bare(:)
       DO ji=1,kjpindex
          vbeta2(ji,:) = frac_evap_intercept(ji,:) * raero(ji) / resist_intercept(ji,:)
          vbeta3(ji,:) = frac_evap_transp(ji,:) * raero(ji) / resist_transp(ji,:)
          vbeta(ji) = vbeta4(ji) + SUM(vbeta2(ji,:) + vbeta3(ji,:))
       ENDDO
       vbeta5(:) = frac_evap_flood(:)
       
       !+++CHECK+++
       ! Variable dimensions xxx_Tot_mean are for a single pixel. Causing 1+1 problems
       ! when running a larger domain. Needs to be corrected when implementing
       ! a global use of the multi-layer energy budget.
       CALL mleb_main (kjit, kjpindex, ldrestart_read, ldrestart_write, &
            index, indexveg, lwdown, swnet, swdown, epot_air, temp_air, &
            u, v, petAcoef, petBcoef, qair, peqAcoef, peqBcoef, pb, rau, &
            vbeta, vbeta1, vbeta2, vbeta3, vbeta3pot, vbeta4, vbeta5, &
            emis, soilflx, soilcap, tq_cdrag, humrel, &
            fluxsens, fluxlat, vevapp, transpir, transpot, vevapnu,vevapwet, &
            vevapsno, vevapflo, temp_sol, tsol_rad, temp_sol_new, qsurf, &
            evapot, evapot_corr, rest_id, hist_id, hist2_id, &
            ok_mleb_history_file, &
            transpir_supply, hydrol_flag, hydrol_flag2, hydrol_flag3,vbeta2sum, &
            vbeta3sum, veget_max, qsol_sat_new, qair_new, &
            veget, &
            Light_Abs_Tot_mean, Light_Alb_Tot_mean, &
            laieff_isotrop, lai_per_level,&
            z_array_out, transpir_supply_column, u_speed, &
            profile_vbeta3, profile_rveget, max_height_store, &
            precip_rain,  pgflux, snowdz, temp_sol_add)
       !+++++++++++
       
    ENDIF !! .NOT. ok_mleb
    



    IF (ok_lake_energy) THEN
       !! 3.2 lake only
       CALL lake_main (kjit,    lalo,      kjpindex,      index,         &
            u  ,v  , qair   , precip_rain , precip_snow ,lwdown,swnet,swdown, pb,temp_air,&                      
            vevapp_lake,   temp_sol_lake, fluxsens_lake, fluxlat_lake)
    
       !! 3.3 Calculation on all the tile
       CALL sechiba_energy_merge_diag (kjpindex, temp_sol_lake, fluxsens_lake, fluxlat_lake,&
            temp_sol_new, fluxsens, fluxlat, frac_lake)   
    ENDIF
       
    
    
     !! 4. Compute hydrology
     !! 4.1 Water balance from CWRR module (11 soil layers)
     CALL hydrol_main (ks,  nvan, avan, mcr, mcs, mcfc, mcw, kjit, kjpindex, &
          & index, indexveg, indexsoil, indexlayer, indexnslm, &
          & temp_sol_new, floodout, runoff, drainage, frac_nobio, totfrac_nobio, frac_snow_nobio, vevapwet, veget, veget_max, njsc, &
          & qsintmax, qsintveg, vevapnu, vevapsno, vevapflo, snow, snow_age, snow_nobio_age,  &
          & tot_melt, transpir, precip_rain, precip_snow, returnflow, reinfiltration, irrigation, &
          & humrel, vegstress, drysoil_frac, evapot, evapot_corr, evap_bare_lim, evap_bare_lim_ns, flood_frac, flood_res, z0m, &
          & shumdiag,shumdiag_perma, k_litt, litterhumdiag, soilcap, soiltile, fraclut, reinf_slope_soil,&
          & rest_id, hist_id, hist2_id,&
          & contfrac, stempdiag, &
          & temp_air, pb, u, v, tq_cdrag, swnet, pgflux, &
          & snowrho, snowtemp, snowgrain, snowdz, snowheat, snowliq, &
          & grndflux,gtemp,tot_bare_soil, &
          & lambda_snow,cgrnd_snow,dgrnd_snow,frac_snow_veg,temp_sol_add, &
          & lambda_ice, cgrnd_ice, dgrnd_ice, ice_sheet_mask, icetemp, icedz, snowmelt, zrainfall,&
          & mc_layh, mcl_layh, tmc_pft, drainage_pft, runoff_pft, swc_pft, soilmoist, &
          & mc_layh_s, mcl_layh_s, soilmoist_s, swc, e_frac, ksoil, &
          & mcs_hydrol, mcfc_hydrol, altmax, root_profile, root_depth, root_deficit, &
          & circ_class_biomass, us, run_off_lic, run_off_lic_frac, circ_class_n, gsmean,&
          & psi_leaf, psi_leaf_next, psi_leaf_midday, psi_sto_leaf_save, psi_sto_wood_save, &
          & psi_root_sup, psi_root_inf, psi_soil_sup, psi_soil_inf, psi_xylem_trunk, psi_xylem_leaf, &
          & psi_xylem_collar, psi_sto_wood, psi_sto_leaf,&
          & mc_i_sup, mc_i_inf, F_absorption, lalo, z_array_out, &
          & qair, qsol_sat_new, rau, frac_evap_flood, frac_sub_snow, frac_evap_transp, cimean, assimtot, Rdtot_pft, vessel_loss,&
          & wtp,mc_peat_above,liqwt_ratio,shumdiag_peat,shumdiag_croppeat,mc_croppeat_above, &
          & shumdiag_man,mc_man_above)

     !! 6. Compute surface variables (emissivity, albedo and roughness)
     !+++CHECK+++
     ! Variable dimensions xxx_Tot_mean are for a single pixel. Causing 1+1 problems
     ! when running a larger domain. Needs to be corrected when implementing
     ! a global use of the multi-layer energy budget.
     CALL condveg_main (kjit, kjpindex, index, rest_id, hist_id, hist2_id, &
         lalo, neighbours, resolution, contfrac, &
         veget, veget_max, frac_nobio, totfrac_nobio, &
         zlev, snow, snow_age, snow_nobio_age, &
         drysoil_frac, height, height_dom, snowdz, snowrho, tot_bare_soil, &
         temp_air, pb, u, v, &
         lai, &
         emis, albedo, z0m, z0h, roughheight, &
         frac_snow_veg, frac_snow_nobio, coszang, &
         Light_Abs_Tot, Light_Tran_Tot, laieff_fit, &
         Light_Abs_Tot_mean, Light_Alb_Tot_mean, &
         laieff_isotrop)
    !+++++++++++

    !! 7. Compute soil thermodynamics

    DO jv = 1,nvm
      mc_layh_pft(:,:,jv) = mc_layh_s(:,:,pref_soil_veg(jv))
      mcl_layh_pft(:,:,jv) = mcl_layh_s(:,:,pref_soil_veg(jv))
      soilmoist_pft(:,:,jv) = soilmoist_s(:,:,pref_soil_veg(jv))
    END DO

    CALL thermosoil_main (kjit, kjpindex, &
         index, indexgrnd, &
         temp_sol_new, snow, soilcap, soilflx, &
         shumdiag_perma, stempdiag, ftempdiag, ptnlev1, rest_id, hist_id, hist2_id, &
         snowdz,snowrho,snowtemp,gtemp,pb,&
         mc_layh, mcl_layh, soilmoist, &
         mc_layh_pft, mcl_layh_pft,soilmoist_pft, njsc, &
         soiltile, &
         depth_organic_soil, heat_Zimov, tdeep, hsdeep,&
         som_total, veget_max, &
         frac_snow_veg,frac_snow_nobio,totfrac_nobio,temp_sol_add, &
         lambda_snow, cgrnd_snow, dgrnd_snow, &
         lambda_ice,  cgrnd_ice,  dgrnd_ice, icetemp, icedz, ice_sheet_mask)
    
    !! 8. Compute river routing 
    IF ( river_routing .AND. nbp_glo .GT. 1) THEN
       !! 8.1 River routing
       CALL routing_wrapper_main (kjit, kjpindex, index, &
            & lalo, neighbours, resolution, contfrac, totfrac_nobio, veget_max, floodout, runoff, &
            & drainage, transpot, precip_rain, humrel, k_litt, flood_frac, flood_res, &
            & stempdiag, ftempdiag, reinf_slope, returnflow, reinfiltration, irrigation, riverflow, coastalflow, &
            & rest_id, hist_id, hist2_id, soiltile, root_deficit, irrigated_next, irrig_frac, fraction_aeirrig_sw)
    ELSE
       !! 8.2 No routing, set variables to zero
       riverflow(:) = zero
       coastalflow(:) = zero
       returnflow(:) = zero
       reinfiltration(:) = zero
       irrigation(:) = zero
       flood_frac(:) = zero
       flood_res(:) = zero

       CALL xios_orchidee_send_field("coastalflow",coastalflow/dt_sechiba)
       CALL xios_orchidee_send_field("riverflow",riverflow/dt_sechiba)
    ENDIF


    !! 9. Compute slow processes (i.e. 'daily' and annual time step)
    CALL slowproc_main (kjit, kjpij, kjpindex, njsc, &
         index, indexveg, lalo, neighbours, resolution, contfrac, soiltile, fraclut, nwdFraclut, &
         temp_air, temp_sol, stempdiag, precip_longterm, &
         vegstress, humrel, &
         shumdiag, litterhumdiag, precip_rain, precip_snow, pb, gpp, VPD, &
         tmc_pft, drainage_pft, runoff_pft, swc_pft, deadleaf_cover, &
         assim_param, qsintveg, &
         frac_age, height, lai, veget, frac_nobio, veget_max, totfrac_nobio, qsintmax, &
         rest_id, hist_id, hist2_id, rest_id_stom, hist_id_stom, hist_id_stom_IPCC, &
         fco2_flux_out, fco2_lu_out, fco2_wh_out, fco2_ha_out, &
         temp_growth, tot_bare_soil, &
         tdeep, hsdeep, snow, heat_Zimov, &
         sfluxCH4_deep, sfluxCO2_deep, &
         som_total, snowdz, snowrho, altmax, depth_organic_soil, &
         circ_class_biomass, circ_class_n, &
         lai_per_level, max_height_store, laieff_fit, &
         Light_Abs_Tot, Light_Tran_Tot, laieff_isotrop, &
         z_array_out, transpir, transpir_mod, &
         coszang, coszang_noon, &
         stressed, unstressed, &
         u, v, loss_gain, veget_max_new, frac_nobio_new, failed_vegfrac, &
         mcs_hydrol, mcfc_hydrol,vessel_loss, height_dom, height_inv, dia_dom, dia_inv, &
         ba_inv, ind_inv, root_profile, root_depth, us, Pgap_cumul, z0m, &
         irrigated_next, irrig_frac, reinf_slope, reinf_slope_soil, &
         wtp,  mc_peat_above,& 
         liqwt_ratio,  shumdiag_peat,  shumdiag_croppeat,&
         mc_croppeat_above, &
         shumdiag_man, mc_man_above)

    !+++CHECK+++
    ! slowproc_main basically ends with a CALL to slowproc_veget and 
    ! check_veget if do_slow. LCC takes place at the end of the first day 
    ! so one could expect do_slow = TRUE. Are the following CALLs really
    ! needed. Seems that they should not do too much. Consider deleting.

    ! Update variables such as veget, frac_nobio, soiltile, lai_per_level
    ! after LCC has been done in stomatelpj. Some of these variables are
    ! used in calculating variables for the history files. 
    IF (done_stomate_lcchange) THEN
       
       !! Calculation of veget and soiltile.
       CALL slowproc_veget (kjpindex, lai_per_level, z_array_out, &
            coszang_noon, circ_class_biomass, circ_class_n, frac_nobio, totfrac_nobio, &
            veget_max, veget, soiltile, tot_bare_soil, &
            fraclut, nwdFraclut, Pgap_cumul, precip_longterm)

       !! Do some basic tests on the surface fractions updated above
       CALL check_veget(kjpindex, frac_nobio, veget_max, &
            veget, tot_bare_soil, soiltile, failed_vegfrac)

       done_stomate_lcchange = .FALSE.
       
    END IF

    !! 10. Update the temperature (temp_sol) with newly computed values
    CALL sechiba_end (kjpindex, temp_sol_new, temp_sol)

    !! 14. If it is the last time step, write restart files
    IF (ldrestart_write) THEN
       !+++CHECK+++
       ! Variable dimensions xxx_Tot_mean are for a single pixel. Causing 1+1 problems
       ! when running a larger domain. Needs to be corrected when implementing
       ! a global use of the multi-layer energy budget.
       CALL sechiba_finalize( &
            kjit,     kjpij,  kjpindex, index,   rest_id, &
            tq_cdrag, vevapp, fluxsens, fluxlat, tsol_rad, &
            albedo, &
            Light_Abs_Tot_mean,   Light_Alb_Tot_mean, &
            wtp,liqwt_ratio)
       !+++++++++++
    END IF
    
    !! 11. Write internal variables to output fields
    z0m_out(:) = z0m(:)
    z0h_out(:) = z0h(:)
    emis_out(:) = emis(:)
    qsurf_out(:) = qsurf(:)
    veget_out(:,:)  = veget(:,:)
    lai_out(:,:) = lai(:,:)
    height_out(:,:) = height(:,:)
 
    !! 12. Write global variables to history files
    sum_height_treefrac(:) = zero
    sum_treefrac(:) = zero
    sum_height_grassfrac(:) = zero
    sum_grassfracC3(:) = zero
    sum_grassfracC4(:) = zero
    sum_cropfracC3(:) = zero
    sum_cropfracC4(:) = zero
    sum_treeFracNdlEvg(:) = zero
    sum_treeFracBdlEvg(:) = zero
    sum_treeFracNdlDcd(:) = zero
    sum_treeFracBdlDcd(:) = zero
    DO jv = 2, nvm 
       IF (is_tree(jv) .AND. natural(jv)) THEN
          sum_treefrac(:) = sum_treefrac(:) + veget_max(:,jv)
          sum_height_treefrac(:) = sum_height_treefrac(:) + veget_max(:,jv)*height(:,jv)
       ELSE IF ((.NOT. is_tree(jv))  .AND. natural(jv)) THEN
          ! Grass
          IF (is_c4(jv)) THEN
             sum_grassfracC4(:) = sum_grassfracC4(:) + veget_max(:,jv)
          ELSE
             sum_grassfracC3(:) = sum_grassfracC3(:) + veget_max(:,jv)
          END IF
          sum_height_grassfrac(:) = sum_height_grassfrac(:) + veget_max(:,jv)*height(:,jv)
       ELSE 
          ! Crop and trees not natural
          IF (is_c4(jv)) THEN
             sum_cropfracC4(:) = sum_cropfracC4(:) + veget_max(:,jv)
          ELSE
             sum_cropfracC3(:) = sum_cropfracC3(:) + veget_max(:,jv)
          END IF
       ENDIF

       IF (is_tree(jv)) THEN
          IF (is_evergreen(jv)) THEN
             IF (is_needleleaf(jv)) THEN
                ! Fraction for needleleaf evergreen trees (treeFracNdlEvg)
                sum_treeFracNdlEvg(:) = sum_treeFracNdlEvg(:) + veget_max(:,jv)
             ELSE
                ! Fraction for broadleaf evergreen trees (treeFracBdlEvg)
                sum_treeFracBdlEvg(:) = sum_treeFracBdlEvg(:) + veget_max(:,jv)
             END IF
          ELSE IF (is_deciduous(jv)) THEN
             IF (is_needleleaf(jv)) THEN
                ! Fraction for needleleaf deciduous trees (treeFracNdlDcd)
                sum_treeFracNdlDcd(:) = sum_treeFracNdlDcd(:) + veget_max(:,jv)
             ELSE 
                ! Fraction for broadleafs deciduous trees (treeFracBdlDcd)
                sum_treeFracBdlDcd(:) = sum_treeFracBdlDcd(:) + veget_max(:,jv)
             END IF
          END IF
       END IF
    ENDDO          

    ! === DBG-CHAIN : dump chaine causale SECHIBA au pas dt_sechiba (comparaison Julia) ===
    ! ji=1, PFT=7. Unites alignees Julia : precip/gpp en /s, transpir/evapnu/vevapwet en /dt,
    ! pb en hPa, tmc=SUM(soilmoist) kg/m2. Grep 'CHAINCSV' dans le log puis strip tag -> CSV.
    ! Colonnes = step,temp_air,qair,precip,swdown,lwdown,pb,wind,tmc,humrel,vegstress,temp_sol,
    !            netrad,fluxsens,fluxlat,soilflx,transpir,evapnu,vevapwet,gpp,raero,rveget
    ! Harnais 1 MOIS (amazon PFT=2) : garde kjit<=1488 (31j x 48 pas), PFT=2 (tropical evergreen).
    IF (kjit <= 1488) THEN
    WRITE(numout,'(A,22(",",ES16.8))') 'CHAINCSV', &
         REAL(kjit,r_std), &
         temp_air(1), qair(1), precip_rain(1)/dt_sechiba, swdown(1), lwdown(1), pb(1), &
         SQRT(u(1)**2+v(1)**2), &
         SUM(soilmoist(1,:)), humrel(1,2), vegstress(1,2), &
         temp_sol_new(1), netrad(1), fluxsens(1), fluxlat(1), soilflx(1), &
         transpir(1,2), vevapnu(1), vevapwet(1,2), &
         gpp(1,2)/dt_sechiba, raero(1), rveget(1,2)
    ENDIF

    histvar(:)=zero
    DO jv = 2, nvm
       IF (is_deciduous(jv)) THEN
          histvar(:) = histvar(:) + veget_max(:,jv)*100*contfrac
       ENDIF
    ENDDO
    CALL xios_orchidee_send_field("treeFracPrimDec",histvar)

    histvar(:)=zero
    DO jv = 2, nvm
       IF (is_evergreen(jv)) THEN
          histvar(:) = histvar(:) + veget_max(:,jv)*100*contfrac
       ENDIF
    ENDDO
    CALL xios_orchidee_send_field("treeFracPrimEver",histvar)

    histvar(:)=zero
    DO jv = 2, nvm
       IF ( .NOT.(is_c4(jv)) ) THEN
          histvar(:) = histvar(:) + veget_max(:,jv)*100*contfrac
       ENDIF
    ENDDO
    CALL xios_orchidee_send_field("c3PftFrac",histvar)

    histvar(:)=zero
    DO jv = 2, nvm
       IF ( is_c4(jv) ) THEN
          histvar(:) = histvar(:) + veget_max(:,jv)*100*contfrac
       ENDIF
    ENDDO
    CALL xios_orchidee_send_field("c4PftFrac",histvar)

    CALL xios_orchidee_send_field("temp_sol_new",temp_sol_new)
    CALL xios_orchidee_send_field("fluxsens",fluxsens)
    CALL xios_orchidee_send_field("fluxlat",fluxlat)
    ! Add XIOS default value where no snow
    DO ji=1,kjpindex 
       IF (snow(ji) .GT. zero) THEN
          snow_age_diag(ji) = snow_age(ji)
          snow_nobio_age_diag(ji,:) = snow_nobio_age(ji,:)
       
          snowage_glob(ji) = snow_age(ji)*frac_snow_veg(ji)*(1-totfrac_nobio(ji)) + &
               SUM(snow_nobio_age(ji,:)*frac_snow_nobio(ji,:)*frac_nobio(ji,:))
          IF (snowage_glob(ji) .NE. 0) snowage_glob(ji) = snowage_glob(ji) / &
               (frac_snow_veg(ji)*(1-totfrac_nobio(ji)) + SUM(frac_snow_nobio(ji,:)*frac_nobio(ji,:)))
       ELSE
          snow_age_diag(ji) = xios_default_val
          snow_nobio_age_diag(ji,:) = xios_default_val
          snowage_glob(ji) = xios_default_val
       END IF
    END DO
    
    CALL xios_orchidee_send_field("snow",snow)
    CALL xios_orchidee_send_field("snowage",snow_age_diag)
    CALL xios_orchidee_send_field("snownobioage",snow_nobio_age_diag)
    CALL xios_orchidee_send_field("snowage_glob",snowage_glob)

    CALL xios_orchidee_send_field("frac_snow", SUM(frac_snow_nobio,2)*totfrac_nobio+frac_snow_veg*(1-totfrac_nobio))
    CALL xios_orchidee_send_field("frac_snow_veg", frac_snow_veg)
    CALL xios_orchidee_send_field("frac_snow_nobio", frac_snow_nobio)
    CALL xios_orchidee_send_field("frac_snow", SUM(frac_snow_nobio,2)*totfrac_nobio+frac_snow_veg*(1-totfrac_nobio))
    CALL xios_orchidee_send_field("frac_snow_veg", frac_snow_veg)
    CALL xios_orchidee_send_field("frac_snow_nobio", frac_snow_nobio)
    CALL xios_orchidee_send_field("pgflux",pgflux)
    CALL xios_orchidee_send_field("reinf_slope",reinf_slope)
    CALL xios_orchidee_send_field("njsc",REAL(njsc, r_std))
    CALL xios_orchidee_send_field("vegetfrac",veget)
    CALL xios_orchidee_send_field("maxvegetfrac",veget_max)
    CALL xios_orchidee_send_field("irrigmap_dyn",irrigated_next)
    CALL xios_orchidee_send_field("aei_sw",fraction_aeirrig_sw)
    CALL xios_orchidee_send_field("nobiofrac",frac_nobio)
    CALL xios_orchidee_send_field("soiltile",soiltile)
    CALL xios_orchidee_send_field("rstruct",rstruct)

    !-
    ! Write effective LAI to the output file.  We choose the isotropic LAI
    ! here since that is integrated over all solar angles.  Comparisons
    ! come from satellies and are dependent on the viewing angle, which is not
    ! available from observations, so this becomes the LAI as viewed from
    ! 60 degrees instead of the typical LAI, which is viewed from 90 degrees.
    ! laieff_isotrop is calculated for every level of LAI, but we just
    ! want to print out a single value for the canopy.
    ! This is a variable that is part stomate, part sechiba.  Print it
    ! out in sechiba because it's used here, even though it depends
    ! completely on stomate.
    temporary_array(:,:)=0.0d0
    DO ilevel = 1,nlevels_tot
       temporary_array(:,:)=temporary_array(:,:)+laieff_isotrop(:,ilevel,:)
    ENDDO
    CALL xios_orchidee_send_field("EFFECTIVE_LAI",temporary_array(:,:))

    ! NOTE! These values do not account for atm_to_bm. They are supposed to
    ! be different from gpp send to the history files from stomate. 
    CALL xios_orchidee_send_field("gpp",gpp/dt_sechiba)
    CALL xios_orchidee_send_field("gpp_ipcc2",SUM(gpp,dim=2)/dt_sechiba)

    histvar(:)=zero
    DO jv = 2, nvm
       IF ( .NOT. is_tree(jv) .AND. natural(jv) ) THEN
          histvar(:) = histvar(:) + gpp(:,jv)
       ENDIF
    ENDDO
    CALL xios_orchidee_send_field("gppgrass",histvar/dt_sechiba)

    histvar(:)=zero
    DO jv = 2, nvm
       IF ( (.NOT. is_tree(jv)) .AND. (.NOT. natural(jv)) ) THEN
          histvar(:) = histvar(:) + gpp(:,jv)
       ENDIF
    ENDDO
    CALL xios_orchidee_send_field("gppcrop",histvar/dt_sechiba)

    histvar(:)=zero
    DO jv = 2, nvm
       IF ( is_tree(jv) ) THEN
          histvar(:) = histvar(:) + gpp(:,jv)
       ENDIF
    ENDDO
    CALL xios_orchidee_send_field("gpptree",histvar/dt_sechiba)
    !++++++++++++
    IF (printlev_loc .GT. 3) THEN
       WRITE(numout,*)'resist_transp(test_grid/pft):',resist_transp(test_grid,test_pft)
    ENDIF

    CALL xios_orchidee_send_field("drysoil_frac",drysoil_frac)
    CALL xios_orchidee_send_field("vevapflo",vevapflo/dt_sechiba)
    CALL xios_orchidee_send_field("k_litt",k_litt)
    CALL xios_orchidee_send_field("resist_equiv_evap",resist_equiv_evap)
    CALL xios_orchidee_send_field("resist_equiv_evap_veg",resist_equiv_evap_veg)
    CALL xios_orchidee_send_field("resist_equiv_tot",resist_equiv_tot)
    CALL xios_orchidee_send_field("frac_sub_snow",frac_sub_snow)
    CALL xios_orchidee_send_field("frac_evap_intercept",frac_evap_intercept)
    CALL xios_orchidee_send_field("resist_intercept",resist_intercept)
    CALL xios_orchidee_send_field("frac_evap_transp",frac_evap_transp)
    CALL xios_orchidee_send_field("resist_transp",resist_transp)
    CALL xios_orchidee_send_field("frac_evap_bare",frac_evap_bare)
    CALL xios_orchidee_send_field("frac_evap_flood",frac_evap_flood)
    CALL xios_orchidee_send_field("VPD",VPD)
    CALL xios_orchidee_send_field("air_relhum",air_relhum)
    CALL xios_orchidee_send_field("f_gco2",f_gco2_out/nlevels_tot)
    CALL xios_orchidee_send_field("gsmean",gsmean)
    CALL xios_orchidee_send_field("cimean",cimean)
    CALL xios_orchidee_send_field("rveget",rveget)
 
    IF ( ok_c13 ) THEN
       WHERE (ABS(delta_c13_assim) < min_sechiba)
           delta_c13_assim = xios_default_val
       ENDWHERE
       WHERE (leaf_ci_out < min_sechiba)
           leaf_ci_out = xios_default_val
       ENDWHERE
       WHERE (gpp_day < min_sechiba)
           gpp_day = xios_default_val
       ENDWHERE
    ELSE
       delta_c13_assim(:,:)=-9999
       leaf_ci_out(:,:)=-9999
       gpp_day(:,:)=-9999
    ENDIF

    CALL xios_orchidee_send_field('delta_c13_assim', delta_c13_assim)
    CALL xios_orchidee_send_field('leaf_ci_out', leaf_ci_out)
    CALL xios_orchidee_send_field('gpp_day', gpp_day)
 
    histvar(:)=SUM(vevapwet(:,:),dim=2)
    CALL xios_orchidee_send_field("evspsblveg",histvar/dt_sechiba)
    histvar(:)= vevapnu(:)+vevapsno(:)
    CALL xios_orchidee_send_field("evspsblsoi",histvar/dt_sechiba)
    histvar(:)=SUM(transpir(:,:),dim=2)
    CALL xios_orchidee_send_field("tran",histvar/dt_sechiba)

    ! For CMIP6 data request: the following fractions are fractions of the total grid-cell,
    ! which explains the multiplication by contfrac
    CALL xios_orchidee_send_field("treeFrac",sum_treefrac(:)*100*contfrac(:))

    WHERE (sum_treefrac .GT. zero)
       vegheighttree=sum_height_treefrac(:)/sum_treefrac(:)
    ELSEWHERE
       vegheighttree=xios_default_val
    END WHERE
    CALL xios_orchidee_send_field("vegHeightTree",vegheighttree)

    WHERE ((sum_grassfracC3 .GT. zero) .OR. (sum_grassfracC4 .GT. zero))
       vegheightgrass=sum_height_grassfrac/(sum_grassfracC3+sum_grassfracC4)
    ELSEWHERE
       vegheightgrass=xios_default_val
    END WHERE   
    CALL xios_orchidee_send_field("vegHeightGrass",vegheightgrass)

    CALL xios_orchidee_send_field("vegHeight",SUM(veget_max(:,:)*height(:,:), dim=2))
    
    CALL xios_orchidee_send_field("grassFracC3",sum_grassfracC3(:)*100*contfrac(:))
    CALL xios_orchidee_send_field("grassFracC4",sum_grassfracC4(:)*100*contfrac(:))
    CALL xios_orchidee_send_field("cropFracC3",sum_cropfracC3(:)*100*contfrac(:))
    CALL xios_orchidee_send_field("cropFracC4",sum_cropfracC4(:)*100*contfrac(:))
    CALL xios_orchidee_send_field("treeFracNdlEvg",sum_treeFracNdlEvg(:)*100*contfrac(:))
    CALL xios_orchidee_send_field("treeFracBdlEvg",sum_treeFracBdlEvg(:)*100*contfrac(:))
    CALL xios_orchidee_send_field("treeFracNdlDcd",sum_treeFracNdlDcd(:)*100*contfrac(:))
    CALL xios_orchidee_send_field("treeFracBdlDcd",sum_treeFracBdlDcd(:)*100*contfrac(:))

    histvar(:)=veget_max(:,1)*100*contfrac(:)
    CALL xios_orchidee_send_field("baresoilFrac",histvar)
    histvar(:)=SUM(frac_nobio(:,1:nnobio),dim=2)*100*contfrac(:)
    CALL xios_orchidee_send_field("residualFrac",histvar)

    ! For CMIP6 data request: cnc = canopy cover fraction over land area
    histvar(:)=zero
    DO jv=2,nvm
       histvar(:) = histvar(:) + veget_max(:,jv)*100
    END DO
    CALL xios_orchidee_send_field("cnc",histvar)
    
    CALL xios_orchidee_send_field("tsol_rad",tsol_rad-273.15)
    CALL xios_orchidee_send_field("qsurf",qsurf)
    CALL xios_orchidee_send_field("emis",emis)
    CALL xios_orchidee_send_field("z0m",z0m)
    CALL xios_orchidee_send_field("z0h",z0h)
    CALL xios_orchidee_send_field("roughheight",roughheight)
    CALL xios_orchidee_send_field("u",u)
    CALL xios_orchidee_send_field("v",v)
    CALL xios_orchidee_send_field("zlev",zlev)

    ! calculate laimean (pixel including no_bio)
    histvar(:)=zero   
    DO ji = 1, kjpindex
       IF (SUM(veget_max(ji,:)) > zero) THEN
         DO jv=2,nvm
            histvar(ji) = histvar(ji) + veget_max(ji,jv)*lai(ji,jv)/ &
                 SUM(veget_max(ji,:))
         END DO
       END IF
    END DO
    CALL xios_orchidee_send_field("lai",lai)
    CALL xios_orchidee_send_field("LAImean",histvar)

    CALL xios_orchidee_send_field("vevapsno",vevapsno/dt_sechiba)
    CALL xios_orchidee_send_field("vevapp",vevapp/dt_sechiba)
    CALL xios_orchidee_send_field("vevapnu",vevapnu/dt_sechiba)
    CALL xios_orchidee_send_field("transpir",transpir*one_day/dt_sechiba)
    CALL xios_orchidee_send_field("transpot",transpot*one_day/dt_sechiba)
    CALL xios_orchidee_send_field("inter",vevapwet*one_day/dt_sechiba)

    histvar(:)=zero
    DO jv=1,nvm
      histvar(:) = histvar(:) + vevapwet(:,jv)
    ENDDO
    CALL xios_orchidee_send_field("ECanop",histvar/dt_sechiba)
    histvar(:)=zero
    DO jv=1,nvm
      histvar(:) = histvar(:) + transpir(:,jv)
    ENDDO
    CALL xios_orchidee_send_field("TVeg",histvar/dt_sechiba)

    !! Calculate diagnostic variables on Landuse tiles for LUMIP/CMIP6

    ! Calculate fraction of landuse tiles related to the whole grid cell
    DO jv=1,nlut
       histvar2(:,jv) = fraclut(:,jv) * contfrac(:)
    END DO
    CALL xios_orchidee_send_field("fraclut",histvar2)

    CALL xios_orchidee_send_field("nwdFraclut",nwdFraclut(:,:))
   
    ! Calculate GPP on landuse tiles
    ! val_exp is used as missing value where the values are not known i.e. where the tile is not represented 
    ! or for pasture (id_pst) or urban land (id_urb). 
    gpplut(:,:)=0
    DO jv=1,nvm
       IF (natural(jv)) THEN
          gpplut(:,id_psl) = gpplut(:,id_psl) + gpp(:,jv)
       ELSE
          gpplut(:,id_crp) = gpplut(:,id_crp) + gpp(:,jv)
       ENDIF
    END DO

    ! Transform from gC/m2/s into kgC/m2/s
    WHERE (fraclut(:,id_psl)>min_sechiba)
       gpplut(:,id_psl) = gpplut(:,id_psl)/fraclut(:,id_psl)/1000
    ELSEWHERE
       gpplut(:,id_psl) = xios_default_val
    END WHERE
    WHERE (fraclut(:,id_crp)>min_sechiba)
       gpplut(:,id_crp) = gpplut(:,id_crp)/fraclut(:,id_crp)/1000
    ELSEWHERE
       gpplut(:,id_crp) = xios_default_val
    END WHERE
    gpplut(:,id_pst) = xios_default_val
    gpplut(:,id_urb) = xios_default_val

    CALL xios_orchidee_send_field("gpplut",gpplut)

    ! +++ DBEN Output +++
    CALL xios_orchidee_send_field("z0m_dben",z0m)
    CALL xios_orchidee_send_field("z0h_dben",z0h)
    CALL xios_orchidee_send_field("albedo_vis_dben",albedo(:,1))
    CALL xios_orchidee_send_field("albedo_nir_dben",albedo(:,2))
    ! +++++++++++++++++++

    ! IOIPSL output
    IF ( .NOT. almaoutput ) THEN
       ! Write history file in IPSL-format
       CALL histwrite_p(hist_id, 'z0m', kjit, z0m, kjpindex, index)
       CALL histwrite_p(hist_id, 'z0h', kjit, z0h, kjpindex, index)
       CALL histwrite_p(hist_id, 'roughheight', kjit, roughheight, kjpindex, index)
       CALL histwrite_p(hist_id, 'vegetfrac', kjit, veget, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'maxvegetfrac', kjit, veget_max, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'nobiofrac', kjit, frac_nobio, kjpindex*nnobio, indexnobio)

       ! calculate lai and laimean (pixel including no_bio fraction)
       histvar(:)=zero   
       DO ji = 1, kjpindex
          IF (SUM(veget_max(ji,:)) > zero) THEN
             DO jv=2,nvm
                histvar(ji) = histvar(ji) + veget_max(ji,jv)*lai(ji,jv)/ &
                     SUM(veget_max(ji,:))
             END DO
          END IF
       END DO
       CALL histwrite_p(hist_id, 'lai', kjit, lai, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'LAImean', kjit, histvar, kjpindex, index)
       CALL histwrite_p(hist_id, 'subli', kjit, vevapsno, kjpindex, index)
       CALL histwrite_p(hist_id, 'evapnu', kjit, vevapnu, kjpindex, index)
       CALL histwrite_p(hist_id, 'transpir', kjit, transpir, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'inter', kjit, vevapwet, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'drysoil_frac', kjit, drysoil_frac, kjpindex, index)
       CALL histwrite_p(hist_id, 'rveget', kjit, rveget, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'rstruct', kjit, rstruct, kjpindex*nvm, indexveg)

       CALL histwrite_p(hist_id, 'snow', kjit, snow, kjpindex, index)
       CALL histwrite_p(hist_id, 'snowage', kjit, snow_age, kjpindex, index)
       CALL histwrite_p(hist_id, 'snownobioage', kjit, snow_nobio_age, kjpindex*nnobio, indexnobio)

       CALL histwrite_p(hist_id, 'grndflux', kjit, grndflux, kjpindex,index)
       CALL histwrite_p(hist_id, 'snowtemp',kjit,snowtemp,kjpindex*nsnow,indexsnow)
       CALL histwrite_p(hist_id, 'snowliq', kjit,snowliq,kjpindex*nsnow,indexsnow)
       CALL histwrite_p(hist_id, 'snowdz', kjit,snowdz,kjpindex*nsnow,indexsnow)
       CALL histwrite_p(hist_id, 'snowrho', kjit,snowrho,kjpindex*nsnow,indexsnow)
       CALL histwrite_p(hist_id, 'snowgrain',kjit,snowgrain,kjpindex*nsnow,indexsnow)
       CALL histwrite_p(hist_id, 'snowheat',kjit,snowheat,kjpindex*nsnow,indexsnow)
       
       CALL histwrite_p(hist_id,'frac_snow_veg',kjit,frac_snow_veg,kjpindex,index)
       CALL histwrite_p(hist_id, 'frac_snow_nobio', kjit,frac_snow_nobio,kjpindex*nnobio, indexnobio)
       CALL histwrite_p(hist_id, 'pgflux',kjit,pgflux,kjpindex,index)
       CALL histwrite_p(hist_id, 'soiltile',  kjit, soiltile, kjpindex*nstm, indexsoil)

       CALL histwrite_p(hist_id, 'soilindex',  kjit, REAL(njsc, r_std), kjpindex, index)
       CALL histwrite_p(hist_id, 'reinf_slope',  kjit, reinf_slope, kjpindex, index)
       CALL histwrite_p(hist_id, 'k_litt', kjit, k_litt, kjpindex, index)
       
       IF ( do_floodplains ) THEN
          CALL histwrite_p(hist_id, 'evapflo', kjit, vevapflo, kjpindex, index)
          CALL histwrite_p(hist_id, 'flood_frac', kjit, flood_frac, kjpindex, index)
       ENDIF
       
       CALL histwrite_p(hist_id, 'gsmean', kjit, gsmean, kjpindex*nvm, indexveg)    
       CALL histwrite_p(hist_id, 'gpp', kjit, gpp, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'cimean', kjit, cimean, kjpindex*nvm, indexveg)    
       
       IF ( ok_c13 ) THEN 
          CALL histwrite_p(hist_id, 'delta_c13_assim', kjit, delta_c13_assim, & 
          kjpindex*nvm, indexveg) 
          CALL histwrite_p(hist_id, 'leaf_ci_out', kjit, leaf_ci_out, kjpindex*nvm, & 
          indexveg) 
          CALL histwrite_p(hist_id, 'gpp_day', kjit, gpp_day, kjpindex*nvm, & 
          indexveg)             
       ENDIF

       histvar(:)=SUM(vevapwet(:,:),dim=2)
       CALL histwrite_p(hist_id, 'evspsblveg', kjit, histvar, kjpindex, index)

       histvar(:)= vevapnu(:)+vevapsno(:)
       CALL histwrite_p(hist_id, 'evspsblsoi', kjit, histvar, kjpindex, index)

       histvar(:)=SUM(transpir(:,:),dim=2)
       CALL histwrite_p(hist_id, 'tran', kjit, histvar, kjpindex, index)

       histvar(:)= sum_treefrac(:)*100*contfrac(:)
       CALL histwrite_p(hist_id, 'treeFrac', kjit, histvar, kjpindex, index) 

       histvar(:)= (sum_grassfracC3(:)+sum_grassfracC4(:))*100*contfrac(:)
       CALL histwrite_p(hist_id, 'grassFrac', kjit, histvar, kjpindex, index) 

       histvar(:)= (sum_cropfracC3(:)+sum_cropfracC4(:))*100*contfrac(:)
       CALL histwrite_p(hist_id, 'cropFrac', kjit, histvar, kjpindex, index)

       histvar(:)=veget_max(:,1)*100*contfrac(:)
       CALL histwrite_p(hist_id, 'baresoilFrac', kjit, histvar, kjpindex, index)

       histvar(:)=SUM(frac_nobio(:,1:nnobio),dim=2)*100*contfrac(:)
       CALL histwrite_p(hist_id, 'residualFrac', kjit, histvar, kjpindex, index)
    ELSE
       ! Write history file in ALMA format 
       CALL histwrite_p(hist_id, 'vegetfrac', kjit, veget, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'maxvegetfrac', kjit, veget_max, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'nobiofrac', kjit, frac_nobio, kjpindex*nnobio, indexnobio) 
       CALL histwrite_p(hist_id, 'lai', kjit, lai, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'ESoil', kjit, vevapnu, kjpindex, index)
       CALL histwrite_p(hist_id, 'EWater', kjit, vevapflo, kjpindex, index)
       CALL histwrite_p(hist_id, 'SWE', kjit, snow, kjpindex, index)
       histvar(:)=zero
       DO jv=1,nvm
          histvar(:) = histvar(:) + transpir(:,jv)
       ENDDO
       CALL histwrite_p(hist_id, 'TVeg', kjit, histvar, kjpindex, index)
       histvar(:)=zero
       DO jv=1,nvm
          histvar(:) = histvar(:) + vevapwet(:,jv)
       ENDDO
       CALL histwrite_p(hist_id, 'ECanop', kjit, histvar, kjpindex, index)
       !
       CALL histwrite_p(hist_id, 'Z0m', kjit, z0m, kjpindex, index)
       CALL histwrite_p(hist_id, 'Z0h', kjit, z0h, kjpindex, index)
       CALL histwrite_p(hist_id, 'EffectHeight', kjit, roughheight, kjpindex, index)
       !
       IF ( do_floodplains ) THEN
          CALL histwrite_p(hist_id, 'Qflood', kjit, vevapflo, kjpindex, index)
          CALL histwrite_p(hist_id, 'FloodFrac', kjit, flood_frac, kjpindex, index)
       ENDIF

       CALL histwrite_p(hist_id, 'gsmean', kjit, gsmean, kjpindex*nvm, indexveg)    
       CALL histwrite_p(hist_id, 'cimean', kjit, cimean, kjpindex*nvm, indexveg)    
       CALL histwrite_p(hist_id, 'GPP', kjit, gpp, kjpindex*nvm, indexveg)
       
    ENDIF ! almaoutput

    !! 13. Write additional output file with higher frequency
    IF ( hist2_id > 0 ) THEN
       IF ( .NOT. almaoutput ) THEN
!!$          WRITE(numout,*)'vbeta, what is happening',vbeta
          ! Write history file in IPSL-format
          CALL histwrite_p(hist2_id, 'tsol_rad', kjit, tsol_rad, kjpindex, index)
          CALL histwrite_p(hist2_id, 'qsurf', kjit, qsurf, kjpindex, index)
          CALL histwrite_p(hist2_id, 'albedo', kjit, albedo, kjpindex*2, indexalb)
          CALL histwrite_p(hist2_id, 'emis', kjit, emis, kjpindex, index)
          CALL histwrite_p(hist2_id, 'z0m', kjit, z0m, kjpindex, index)
          CALL histwrite_p(hist2_id, 'z0h', kjit, z0h, kjpindex, index)
          CALL histwrite_p(hist2_id, 'roughheight', kjit, roughheight, kjpindex, index)
          CALL histwrite_p(hist2_id, 'vegetfrac', kjit, veget, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'maxvegetfrac', kjit, veget_max, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'nobiofrac', kjit, frac_nobio, kjpindex*nnobio, indexnobio)
          CALL histwrite_p(hist2_id, 'lai', kjit, lai, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'subli', kjit, vevapsno, kjpindex, index)
          IF ( do_floodplains ) THEN
             CALL histwrite_p(hist2_id, 'vevapflo', kjit, vevapflo, kjpindex, index)
             CALL histwrite_p(hist2_id, 'flood_frac', kjit, flood_frac, kjpindex, index)
          ENDIF
          CALL histwrite_p(hist2_id, 'vevapnu', kjit, vevapnu, kjpindex, index)
          CALL histwrite_p(hist2_id, 'transpir', kjit, transpir, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'inter', kjit, vevapwet, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'drysoil_frac', kjit, drysoil_frac, kjpindex, index)
          CALL histwrite_p(hist2_id, 'rveget', kjit, rveget, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'rstruct', kjit, rstruct, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'snow', kjit, snow, kjpindex, index)
          CALL histwrite_p(hist2_id, 'snowage', kjit, snow_age, kjpindex, index)
          CALL histwrite_p(hist2_id, 'snownobioage', kjit, snow_nobio_age, kjpindex*nnobio, indexnobio)

          CALL histwrite_p(hist2_id, 'soilindex',  kjit, REAL(njsc, r_std), kjpindex, index)
          CALL histwrite_p(hist2_id, 'reinf_slope',  kjit, reinf_slope, kjpindex, index)
          
          CALL histwrite_p(hist2_id, 'gsmean', kjit, gsmean, kjpindex*nvm, indexveg)    
          CALL histwrite_p(hist2_id, 'gpp', kjit, gpp, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'cimean', kjit, cimean, kjpindex*nvm, indexveg)    
          
       ELSE
          ! Write history file in ALMA format
          CALL histwrite_p(hist2_id, 'vegetfrac', kjit, veget, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'maxvegetfrac', kjit, veget_max, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'nobiofrac', kjit, frac_nobio, kjpindex*nnobio, indexnobio)
          CALL histwrite_p(hist2_id, 'ESoil', kjit, vevapnu, kjpindex, index)
          IF ( do_floodplains ) THEN
             CALL histwrite_p(hist2_id, 'EWater', kjit, vevapflo, kjpindex, index)
             CALL histwrite_p(hist2_id, 'FloodFrac', kjit, flood_frac, kjpindex, index)
          ENDIF
          CALL histwrite_p(hist2_id, 'SWE', kjit, snow, kjpindex, index)
          histvar(:)=zero
          DO jv=1,nvm
             histvar(:) = histvar(:) + transpir(:,jv)
          ENDDO
          CALL histwrite_p(hist2_id, 'TVeg', kjit, histvar, kjpindex, index)
          histvar(:)=zero
          DO jv=1,nvm
             histvar(:) = histvar(:) + vevapwet(:,jv)
          ENDDO
          CALL histwrite_p(hist2_id, 'ECanop', kjit, histvar, kjpindex, index)
          CALL histwrite_p(hist2_id, 'GPP', kjit, gpp, kjpindex*nvm, indexveg)
          
       ENDIF ! almaoutput
    ENDIF ! hist2_id

  END SUBROUTINE sechiba_main

!!  =============================================================================================================================
!! SUBROUTINE:    sechiba_finalize
!!
!>\BRIEF	  Finalize all modules by calling their "_finalize" subroutines.
!!
!! DESCRIPTION:	  Finalize all modules by calling their "_finalize" subroutines. These subroutines will write variables to 
!!                restart file. 
!!
!! \n
!_ ==============================================================================================================================
  !+++CHECK+++
  ! Variable dimensions xxx_Tot_mean are for a single pixel. Causing 1+1 problems
  ! when running a larger domain. Needs to be corrected when implementing
  ! a global use of the multi-layer energy budget.
  SUBROUTINE sechiba_finalize( &
       kjit,     kjpij,  kjpindex, index,   rest_id, &
       tq_cdrag, vevapp, fluxsens, fluxlat, tsol_rad, &
       albedo, &
       Light_Abs_Tot_mean,   Light_Alb_Tot_mean, &
       wtp, liqwt_ratio)

!! 0.1 Input variables    
    INTEGER(i_std), INTENT(in)                               :: kjit              !! Time step number (unitless)
    INTEGER(i_std), INTENT(in)                               :: kjpij             !! Total size of the un-compressed grid 
                                                                                  !! (unitless)
    INTEGER(i_std), INTENT(in)                               :: kjpindex          !! Domain size - terrestrial pixels only 
                                                                                  !! (unitless)
    INTEGER(i_std),INTENT (in)                               :: rest_id           !! _Restart_ file identifier (unitless)
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)         :: index             !! Indices of the pixels on the map. 
                                                                                  !! Sechiba uses a reduced grid excluding oceans
                                                                                  !! ::index contains the indices of the 
                                                                                  !! terrestrial pixels only! (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)           :: tsol_rad           !! Radiative surface temperature 
                                                                                  !! @tex $(W m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)           :: vevapp             !! Total of evaporation 
                                                                                  !! @tex $(kg m^{-2} days^{-1})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)           :: fluxsens           !! Sensible heat flux 
                                                                                  !! @tex $(W m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)           :: fluxlat            !! Latent heat flux 
                                                                                  !! @tex $(W m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)           :: tq_cdrag           !! Surface drag coefficient (-)
    REAL(r_std),DIMENSION (kjpindex,2), INTENT (in)         :: albedo             !! Surface albedo for visible and near-infrared (unitless, 0-1)
    !+++CHECK+++
    ! Variable dimensions xxx_Tot_mean are for a single pixel. Causing 1+1 problems
    ! when running a larger domain. Needs to be corrected when implementing
    ! a global use of the multi-layer energy budget.
    REAL(r_std),DIMENSION(nlevels_tot)                      :: Light_Abs_Tot_mean !! total light absorption for a given canopy level
    REAL(r_std),DIMENSION(nlevels_tot)                      :: Light_Alb_Tot_mean !! total albedo for a given level
    REAL(r_std),DIMENSION (kjpindex)                        :: wtp
    REAL(r_std),DIMENSION (kjpindex)                        :: liqwt_ratio
    !+++++++++++

!! 0.2 Local variables
    INTEGER(i_std)                                          :: ji, jv		  !! Index (unitless)
    REAL(r_std), DIMENSION(kjpindex)                        :: histvar            !! Computations for history files (unitless)
    CHARACTER(LEN=80)                                       :: var_name           !! To store variables names for I/O (unitless)

! ==============================================================================================================================

    !! Write restart file for the next simulation from SECHIBA and other modules
    IF (printlev_loc>=3) WRITE (numout,*) 'Start sechiba_finalize for writing restart files'

    !! 1. Call diffuco_finalize to write restart files
    CALL diffuco_finalize (kjit, kjpindex, rest_id, rstruct )


    !! 2. Call energy budget to write restart files
    IF (ok_mleb) THEN
       ! Multi-layer energy budget
       CALL  mleb_finalize (kjit, kjpindex, rest_id,                    &
                            evapot, evapot_corr, temp_sol, tsol_rad,    &
                            qsurf,  fluxsens,    fluxlat,  vevapp,      &
                            u_speed, z_array_out,  &
                            max_height_store )
    ELSE
       ! Energy budget using enerbil module
       CALL enerbil_finalize(kjit,   kjpindex,    rest_id,            &
                             evapot, evapot_corr, temp_sol, tsol_rad, &
                             qsurf,  fluxsens,    fluxlat,  vevapp )
    ENDIF


    !! 2.2 Energy Lake only

    IF ( ok_lake_energy) THEN
       CALL lake_finalize (kjit,   kjpindex,    rest_id, temp_sol_lake)
    ENDIF


    !! 3. Call hydrology to write restart files
    CALL hydrol_finalize(kjit,  kjpindex, rest_id,  vegstress,  qsintveg,  humrel,   snow,     snow_age,          &
                         snow_nobio_age,  snowrho,  snowtemp,   snowdz,    snowheat, snowgrain,                   &
                         drysoil_frac,    evap_bare_lim, evap_bare_lim_ns, swc,      ksoil,    root_profile,      &
                         us,    icetemp,  psi_leaf, psi_leaf_next, psi_leaf_midday,  psi_sto_leaf_save,  psi_sto_wood_save, &
                         psi_root_sup,   psi_root_inf, psi_soil_sup,       psi_soil_inf,       psi_xylem_trunk,   &
                         psi_xylem_leaf, psi_xylem_collar,   psi_sto_wood, psi_sto_leaf,       mc_i_sup,          &
                         mc_i_inf,       F_absorption, &
                         wtp,            liqwt_ratio)

     !! 4. Call condveg to write surface variables to restart files
    !+++CHECK+++
    ! Variable dimensions xxx_Tot_mean are for a single pixel. Causing 1+1 problems
    ! when running a larger domain. Needs to be corrected when implementing
    ! a global use of the multi-layer energy budget.
    CALL condveg_finalize (kjit,           kjpindex,            rest_id,                          &
                           z0m,            z0h,                 roughheight,                      &
                           albedo,         Light_Abs_Tot,   Light_Tran_Tot,               &
                           Light_Abs_Tot_mean, Light_Alb_Tot_mean,      &
                           laieff_isotrop)
    !+++++++++++


    !! 5. Call soil thermodynamic to write restart files
    CALL thermosoil_finalize (kjit,    kjpindex, rest_id,   gtemp, &
         soilcap, soilflx, lambda_snow, cgrnd_snow, dgrnd_snow,    &
         lambda_ice, cgrnd_ice, dgrnd_ice)

    !! 6. Add river routing to restart files  
    IF ( river_routing .AND. nbp_glo .GT. 1) THEN
       !! 6.1 Call river routing to write restart files 
       CALL routing_wrapper_finalize( kjit, kjpindex, rest_id, flood_frac, flood_res )
    ELSE
       !! 6.2 No routing, set variables to zero
       reinfiltration(:) = zero
       returnflow(:) = zero
       irrigation(:) = zero
       flood_frac(:) = zero
       flood_res(:) = zero
    ENDIF

    !! 7. Call slowproc_main to add 'daily' and annual variables to restart file
    CALL slowproc_finalize (kjit,      kjpindex,    rest_id,         index,      &
                            njsc,      veget,       frac_nobio,                  &
                            veget_max, reinf_slope, ks,              nvan,       &
                            avan,      mcr,         mcs,             mcfc,       &
                            mcw,       assim_param, frac_age,                    &
                            heat_Zimov,altmax,      depth_organic_soil,          &
                            circ_class_biomass,     circ_class_n,                &
                            lai_per_level,          laieff_fit,      loss_gain,  &
                            veget_max_new,          frac_nobio_new,  fraction_aeirrig_sw, &
                            precip_longterm)
    

    !! 8. Write variables from sechiba module to restart file
    !!    Here are only the variables set which are not already restarted from other modules.
    CALL restput_p (rest_id, 'coszang_noon', nbp_glo, 1, 1, kjit, coszang_noon, 'scatter',  nbp_glo, index_g)

    IF (printlev_loc>=3) WRITE (numout,*) 'sechiba_finalize done'
    
  END SUBROUTINE sechiba_finalize


!! ==============================================================================================================================\n
!! SUBROUTINE 	: sechiba_init
!!
!>\BRIEF        Dynamic allocation of the variables, the dimensions of the 
!! variables are determined by user-specified settings. 
!! 
!! DESCRIPTION  : The domain size (:: kjpindex) is used to allocate the correct
!! dimensions to all variables in sechiba. Depending on the variable, its 
!! dimensions are also determined by the number of PFT's (::nvm), number of 
!! soil types (::nstm), number of non-vegetative surface types (::nnobio),
!! number of soil levels (::ngrnd), number of soil layers in the hydrological 
!! model (i.e. cwrr) (::nslm). Values for these variables are set in
!! constantes_soil.f90 and constantes_veg.f90.\n
!!
!! Memory is allocated for all Sechiba variables and new indexing tables
!! are build making use of both (::kjpij) and (::kjpindex). New indexing tables 
!! are needed because a single pixel can contain several PFTs, soil types, etc.
!! The new indexing tables have separate indices for the different
!! PFTs, soil types, etc.\n
!!
!! RECENT CHANGE(S): None
!! 
!! MAIN OUTPUT VARIABLE(S): Strictly speaking the subroutine has no output 
!! variables. However, the routine allocates memory and builds new indexing 
!! variables for later use.
!!
!! REFERENCE(S)	: None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================ 

  SUBROUTINE sechiba_init (kjit, kjpij, kjpindex, index, rest_id, lalo)

!! 0.1 Input variables
 
    INTEGER(i_std), INTENT (in)                         :: kjit               !! Time step number (unitless)
    INTEGER(i_std), INTENT (in)                         :: kjpij              !! Total size of the un-compressed grid (unitless)
    INTEGER(i_std), INTENT (in)                         :: kjpindex           !! Domain size - terrestrial pixels only (unitless)
    INTEGER(i_std), INTENT (in)                         :: rest_id            !! _Restart_ file identifier (unitless)
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)    :: index              !! Indeces of the points on the map (unitless)
    REAL(r_std),DIMENSION (kjpindex,2), INTENT (in)     :: lalo               !! Geographical coordinates (latitude,longitude) 
                                                                              !! for pixels (degrees)
!! 0.2 Output variables

!! 0.3 Modified variables

!! 0.4 Local variables

    INTEGER(i_std)                                      :: ier                !! Check errors in memory allocation (unitless)
    INTEGER(i_std)                                      :: ji,jv,ilevel,ipts  !! Indeces (unitless)
    CHARACTER(LEN=80)                                   :: var_name           !! To store variables names in restart file
!_ ==============================================================================================================================

!! 1. Initialize variables 
  
    ! Debug
    ! It is good to leave this in here.  It is only written out once,
    ! it takes almost no time, and it's necessary for identifying a 
    ! problem pixel.
    DO ipts=1,kjpindex
       IF (printlev>=1) WRITE(numout,'(A,I6,10F20.10)') 'pixel number to lat/lon: ',ipts,lalo(ipts,1:2)
    ENDDO
    !-
  
    ! Dynamic allocation with user-specified dimensions on first call
    IF (l_first_sechiba) THEN 
       l_first_sechiba=.FALSE.
    ELSE 
       CALL ipslerr_p(3,'sechiba_init',' l_first_sechiba false . we stop ','','')
    ENDIF

    !! 1.1 Initialize 3D vegetation indexation table
    ALLOCATE (indexveg(kjpindex*nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for indexveg','','')

    ALLOCATE (indexlai(kjpindex*(nlai+1)),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for indexlai','','')

    ALLOCATE (indexsoil(kjpindex*nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for indexsoil','','')

    ALLOCATE (indexnobio(kjpindex*nnobio),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for indexnobio','','')

    ALLOCATE (indexgrnd(kjpindex*ngrnd),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for indexgrnd','','')

    ALLOCATE (indexsnow(kjpindex*nsnow),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for indexsnow','','')

    ALLOCATE (indexlayer(kjpindex*nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for indexlayer','','')

    ALLOCATE (indexnslm(kjpindex*nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for indexnslm','','')

    ALLOCATE (indexalb(kjpindex*2),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for indexalb','','')

    !! 1.2  Initialize 1D array allocation with restartable value
    ALLOCATE (flood_res(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for flood_res','','')
    flood_res(:) = undef_sechiba

    ALLOCATE (flood_frac(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for kjpindex','','')
    flood_frac(:) = undef_sechiba

    ALLOCATE (snow(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for snow','','')
    snow(:) = undef_sechiba

    ALLOCATE (snow_age(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for snow_age','','')
    snow_age(:) = undef_sechiba

    ALLOCATE (drysoil_frac(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for drysoil_frac','','')

    ALLOCATE (evap_bare_lim(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for evap_bare_lim','','')

    ALLOCATE (evap_bare_lim_ns(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for evap_bare_lim_ns','','')

    ALLOCATE (evapot(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for evapot','','')
    evapot(:) = undef_sechiba

    ALLOCATE (evapot_corr(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for evapot_corr','','')

    ALLOCATE (humrel(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for humrel','','')
    humrel(:,:) = undef_sechiba

    ALLOCATE (psi_leaf(kjpindex,nvm,ncirc),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for psi_leaf','','')
    psi_leaf(:,:,:) = undef_sechiba

    ALLOCATE (psi_leaf_next(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for psi_leaf_next','','')
    psi_leaf_next(:,:) = undef_sechiba

    ALLOCATE (psi_leaf_midday(kjpindex,nvm,ncirc),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for psi_leaf_midday','','')
    psi_leaf_midday(:,:,:) = undef_sechiba

    ALLOCATE (psi_sto_leaf_save(kjpindex,nvm,ncirc),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for psi_sto_leaf_save','','')
    psi_sto_leaf_save(:,:,:) = undef_sechiba

    ALLOCATE (psi_sto_wood_save(kjpindex,nvm,ncirc),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for psi_sto_wood_save','','')
    psi_sto_wood_save(:,:,:) = undef_sechiba

    ALLOCATE (psi_root_sup(kjpindex,nvm,ncirc),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for psi_root_sup','','')
    psi_root_sup(:,:,:) = undef_sechiba

    ALLOCATE (psi_root_inf(kjpindex,nvm,ncirc),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for psi_root_inf','','')
    psi_root_inf(:,:,:) = undef_sechiba

    ALLOCATE (psi_soil_sup(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for psi_soil_sup','','')
    psi_soil_sup(:,:) = undef_sechiba

    ALLOCATE (psi_soil_inf(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for psi_soil_inf','','')
    psi_soil_inf(:,:) = undef_sechiba

    ALLOCATE (psi_xylem_trunk(kjpindex,nvm,ncirc),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for psi_xylem_trunk','','')
    psi_xylem_trunk(:,:,:) = undef_sechiba

    ALLOCATE (psi_xylem_leaf(kjpindex,nvm,ncirc),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for psi_xylem_leaf','','')
    psi_xylem_leaf(:,:,:) = undef_sechiba

    ALLOCATE (psi_xylem_collar(kjpindex,nvm,ncirc),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for psi_xylem_collar','','')
    psi_xylem_collar(:,:,:) = undef_sechiba

    ALLOCATE (psi_sto_wood(kjpindex,nvm,ncirc),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for psi_sto_wood','','')
    psi_sto_wood(:,:,:) = undef_sechiba

    ALLOCATE (psi_sto_leaf(kjpindex,nvm,ncirc),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for psi_sto_leaf','','')
    psi_sto_leaf(:,:,:) = undef_sechiba

    ALLOCATE (mc_i_sup(kjpindex,nvm,ncirc,nrp),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for mc_i_sup','','')
    mc_i_sup(:,:,:,:) = undef_sechiba

    ALLOCATE (mc_i_inf(kjpindex,nvm,ncirc,nrp),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for mc_i_inf','','')
    mc_i_inf(:,:,:,:) = undef_sechiba

    ALLOCATE (F_absorption(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for F_absorption','','')
    F_absorption(:,:) = undef_sechiba

    ALLOCATE (vegstress(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for vegstress','','')
    vegstress(:,:) = undef_sechiba

    ALLOCATE (njsc(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for njsc','','')
    njsc(:)= undef_int
    
    ALLOCATE (soiltile(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for soiltile','','')

    ALLOCATE (fraclut(kjpindex,nlut),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for fraclut','','')

    ALLOCATE (nwdFraclut(kjpindex,nlut),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for nwdFraclut','','')

    ALLOCATE (reinf_slope(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for reinf_slope','','')

    ALLOCATE (reinf_slope_soil(kjpindex, nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for reinf_slope_soil','','') !

    ALLOCATE (ks(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for ks','','')

    ALLOCATE (nvan(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for nvan ','','')

    ALLOCATE (avan(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for avan','','')

    ALLOCATE (mcr(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for mcr','','')

    ALLOCATE (mcs(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for mcs','','')

    ALLOCATE (mcfc(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for mcfc','','')

    ALLOCATE (mcw(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for mcw','','')

    ALLOCATE (frac_sub_snow(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for frac_sub_snow','','')

    ALLOCATE (us(kjpindex,nvm,nstm,nslm),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable us','','')

    ALLOCATE (frac_evap_bare(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for frac_evap_bare','','')

    ALLOCATE (frac_evap_flood(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for frac_evap_flood','','')

    ALLOCATE (soilcap(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for soilcap','','')

    ALLOCATE (soilflx(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for soilflx','','')

    ALLOCATE (temp_sol(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for temp_sol','','')
    temp_sol(:) = undef_sechiba

    ALLOCATE (qsurf(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for qsurf','','')
    qsurf(:) = undef_sechiba

    !! 1.3 Initialize 2D array allocation with restartable value
    ALLOCATE (qsintveg(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for qsintveg','','')
    qsintveg(:,:) = undef_sechiba

    ALLOCATE (frac_evap_intercept(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for frac_evap_intercept','','')

    ALLOCATE (resist_intercept(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for resist_intercept','','')

    ALLOCATE (frac_evap_transp(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for frac_evap_transp','','')

    ALLOCATE (resist_transp(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for resist_transp','','')

    ALLOCATE (vbeta3pot(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for vbeta3pot','','')

    ALLOCATE (gsmean(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for gsmean','','')

    ALLOCATE (assimtot(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for assimtot','','')

    ALLOCATE (Rdtot_pft(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for Rdtot_pft','','')

    ALLOCATE (cimean(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for cimean','','')

    ALLOCATE (gpp(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for gpp','','')
    gpp(:,:) = undef_sechiba
 
    ALLOCATE (temp_growth(kjpindex),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for temp_growth','','')
    temp_growth(:) = undef_sechiba 

    ALLOCATE (veget(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for veget','','')
    veget(:,:)=undef_sechiba

    ALLOCATE (veget_max(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for veget_max','','')

    ALLOCATE (tot_bare_soil(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for tot_bare_soil','','')

    ALLOCATE (lai(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for lai','','')
    lai(:,:)=undef_sechiba

    ALLOCATE (laieff_fit(kjpindex,nvm,nlevels_tot),stat=ier)
    IF (ier.NE.0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for laieff_fit','','')
    CALL laieff_type_init(kjpindex, nlevels_tot, laieff_fit)

    ALLOCATE (frac_age(kjpindex,nvm,nleafages),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for frac_age','','')
    frac_age(:,:,:)=undef_sechiba

    ALLOCATE (height(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for height','','')
    height(:,:)=undef_sechiba

    ALLOCATE (height_dom(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for height_dom','','')
    height_dom(:,:)=undef_sechiba

    ALLOCATE (dia_dom(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for dia_dom','','')
    dia_dom(:,:)=undef_sechiba

    ALLOCATE (height_inv(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for height_inv','','')
    height_inv(:,:)=undef_sechiba

    ALLOCATE (dia_inv(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for dia_inv','','')
    dia_inv(:,:)=undef_sechiba

    ALLOCATE (ba_inv(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for ba_inv','','')
    ba_inv(:,:)=undef_sechiba

    ALLOCATE (ind_inv(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for ind_inv','','')
    ind_inv(:,:)=undef_sechiba

    ALLOCATE (frac_nobio(kjpindex,nnobio),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for frac_nobio','','')
    frac_nobio(:,:) = undef_sechiba

    ALLOCATE (snow_nobio_age(kjpindex,nnobio),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for snow_nobio_age','','')
    snow_nobio_age(:,:) = undef_sechiba

    ALLOCATE (assim_param(kjpindex,nvm,npco2),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for assim_param','','')

    !! 1.4 Initialize 1D array allocation 
    ALLOCATE (vevapflo(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for vevapflo','','')
    vevapflo(:)=zero

    ALLOCATE (vevapsno(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for vevapsno','','')

    ALLOCATE (vevapnu(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for vevapnu','','')

    ALLOCATE (totfrac_nobio(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for totfrac_nobio','','')

    ALLOCATE (floodout(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for floodout','','')

    ALLOCATE (runoff(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for runoff','','')

    ALLOCATE (drainage(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for drainage','','')

    ALLOCATE (returnflow(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for returnflow','','')
    returnflow(:) = zero

    ALLOCATE (reinfiltration(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for reinfiltration','','')
    reinfiltration(:) = zero

    ALLOCATE (irrigation(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for irrigation','','')
    irrigation(:) = zero

    ALLOCATE (z0h(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for z0h','','')

    ALLOCATE (z0m(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for z0m','','')

    ALLOCATE (roughheight(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for roughheight','','')

    ALLOCATE (emis(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for emis','','')

    ALLOCATE (tot_melt(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for tot_melt','','')

    ALLOCATE (resist_equiv_evap(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for resist_equiv_evap','','')

    ALLOCATE (resist_equiv_evap_veg(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for resist_equiv_evap_veg','','')

    ALLOCATE (resist_equiv_tot(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for resist_equiv_tot','','')

    ALLOCATE (raero(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for raero','','')

    ALLOCATE (rau(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for rau','','')

    ALLOCATE (deadleaf_cover(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for deadleaf_cover','','')

    ALLOCATE (stempdiag(kjpindex, nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for stempdiag','','')

    ALLOCATE (co2_to_bm(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for co2_to_bm','','')
    
    ALLOCATE (ftempdiag(kjpindex, ngrnd),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for ftempdiag','','')

    ALLOCATE (shumdiag(kjpindex,nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for shumdiag','','')
    
    ALLOCATE (shumdiag_perma(kjpindex,nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for shumdiag_perma','','')

    ALLOCATE (litterhumdiag(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for litterhumdiag','','')

    ALLOCATE (ptnlev1(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for ptnlev1','','')

    ALLOCATE (k_litt(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for k_litt','','')

    ALLOCATE (netrad(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for netrad','','')

    ALLOCATE (lwabs(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for lwabs','','')

    ALLOCATE (lwnet(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for lwnet','','')

    ALLOCATE (fluxsubli(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for fluxsubli','','')

    ALLOCATE (precip_longterm(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for precip_longterm','','')

    !! 1.5 Initialize 2D array allocation
    ALLOCATE (vevapwet(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for vevapwet','','')
    vevapwet(:,:)=undef_sechiba

    ALLOCATE (transpir(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for transpir','','')

    ALLOCATE (transpot(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for transpot','','')

    ALLOCATE (transpir_mod(kjpindex,nvm),stat=ier)
    IF (ier.NE.0) CALL ipslerr_p (3,'sechiba_init','Pb in alloc for transpir_mod','','')
    transpir_mod(:,:) = zero 


    ALLOCATE (stressed(kjpindex,nvm),stat=ier)
    IF (ier.NE.0) CALL ipslerr_p (3,'sechiba_init','Pb in alloc for stressed','','')
    stressed(:,:) = zero

    ALLOCATE (unstressed(kjpindex,nvm),stat=ier)
    IF (ier.NE.0) CALL ipslerr_p (3,'sechiba_init','Pb in alloc for unstressed','','')
    unstressed(:,:) = zero

    ALLOCATE (transpir_supply(kjpindex,nvm),stat=ier)
    IF (ier.NE.0) CALL ipslerr_p (3,'sechiba_init','Pb in alloc for transpir_supply','','')
    transpir_supply(:,:) = zero

    ALLOCATE (vir_transpir_supply(kjpindex,nvm),stat=ier)
    IF (ier.NE.0) CALL ipslerr_p (3,'sechiba_init','Pb in alloc for vir_transpir_supply','','')
    vir_transpir_supply(:,:) = zero

    ALLOCATE (transpir_supply_column(nlevels_tot,kjpindex,nvm),stat=ier)
    IF (ier.NE.0) CALL ipslerr_p (3,'sechiba_init','Pb in alloc for transpir_supply_column','','')
    transpir_supply_column(:,:,:) = zero

    ALLOCATE (vessel_loss(kjpindex,nvm,ncirc),stat=ier)
    IF (ier.NE.0) CALL ipslerr_p (3,'sechiba_init','Pb in alloc for vessel_loss','','')

    ALLOCATE (e_frac(kjpindex,nvm,nslm,nstm),stat=ier)
    IF (ier.NE.0) CALL ipslerr_p (3,'sechiba_init','Pb in alloc for e_frac','','')
    e_frac(:,:,:,:) = zero

    ALLOCATE (qsintmax(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for qsintmax','','')

    ALLOCATE (rveget(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for rveget','','')

    ALLOCATE (rstruct(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for rstruct','','')

    ALLOCATE (pgflux(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for pgflux','','')
    pgflux(:)= 0.0

    ALLOCATE (cgrnd_snow(kjpindex,nsnow),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for cgrnd_snow','','')
    cgrnd_snow(:,:) = 0

    ALLOCATE (dgrnd_snow(kjpindex,nsnow),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for dgrnd_snow','','')
    dgrnd_snow(:,:) = 0

    ALLOCATE (lambda_snow(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for lambda_snow','','')
    lambda_snow(:) = 0
    
    ALLOCATE (cgrnd_ice(kjpindex,nice),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for cgrnd_ice','','')

    ALLOCATE (dgrnd_ice(kjpindex,nice),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for dgrnd_ice','','')

    ALLOCATE (lambda_ice(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for lambda_ice','','')
    
    ALLOCATE (ice_sheet_mask(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for ice_sheet_mask','','')
    
    ALLOCATE (icetemp(kjpindex,nice),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for icetemp','','')
    
    ALLOCATE (icedz(kjpindex,nice),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for icedz','','')

    ALLOCATE (temp_sol_add(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for temp_sol_add','','')

    ALLOCATE (qsol_sat_new(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for qsol_sat_new','','')

    ALLOCATE (qair_new(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for qair_new','','')

    ALLOCATE (gtemp(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for gtemp','','')

    ALLOCATE (frac_snow_veg(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for frac_snow_veg','','')

    ALLOCATE (frac_snow_nobio(kjpindex,nnobio),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for frac_snow_nobio','','')

    ALLOCATE (snowrho(kjpindex,nsnow),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for snowrho','','')

    ALLOCATE (snowheat(kjpindex,nsnow),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for snowheat','','')

    ALLOCATE (snowgrain(kjpindex,nsnow),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for snowgrain','','')

    ALLOCATE (snowtemp(kjpindex,nsnow),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for snowtemp','','')

    ALLOCATE (snowdz(kjpindex,nsnow),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for snowdz','','')

    ALLOCATE(max_height_store(kjpindex,nvm),stat=ier)
    IF (ier.NE.0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for max_height_store','','')

    ALLOCATE (warnings(kjpindex,nvm,nwarns),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in warnings allocation. We stop. We need kjpindex x nvm x nwarns words = ',&
            & kjpindex,' x ' ,nvm, ' x ',nwarns,' = ',kjpindex*nvm*nwarns
      CALL ipslerr_p (3,'sechiba','sechiba_init','','')
    END IF
    ! We aren't restarting the warnings at the moment.
    warnings(:,:,:)=zero

    !! 1.5b Initialize 3D array allocation for albedo
    ALLOCATE (Light_Abs_Tot(kjpindex,nvm,nlevels_tot),stat=ier)
    IF (ier.NE.0) CALL ipslerr_p (3,'sechiba_init', &
         'Pb in alloc for Light_Abs_Tot','','')

    ALLOCATE (Light_Tran_Tot(kjpindex,nvm,nlevels_tot),stat=ier)
    IF (ier.NE.0) CALL ipslerr_p (3,'sechiba_init', &
         'Pb in alloc for Light_Tran_Tot','','')

    ALLOCATE (Light_Abs_Tot_mean(nlevels_tot),stat=ier)
    IF (ier.NE.0) CALL ipslerr_p (3,'sechiba_init', &
         'Pb in alloc for Light_Abs_Tot_mean','','')

    ALLOCATE (Light_Alb_Tot_mean(nlevels_tot),stat=ier)
    IF (ier.NE.0) CALL ipslerr_p (3,'sechiba_init', &
         'Pb in alloc for Light_Alb_Tot_mean','','')

    ALLOCATE (laieff_isotrop(kjpindex,nlevels_tot,nvm),stat=ier)
    IF (ier.NE.0) CALL ipslerr_p (3,'sechiba_init', &
         'Pb in alloc for laieff_isotrop','','')

    ALLOCATE (mc_layh(kjpindex, nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for mc_layh','','')

    ALLOCATE (mcl_layh(kjpindex, nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for mcl_layh','','')

    ALLOCATE (soilmoist(kjpindex, nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for soilmoist','','')

    ALLOCATE (mcl_layh_s(kjpindex, nslm, nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for mcl_layh_s','','')

    ALLOCATE (mc_layh_s(kjpindex, nslm, nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for mc_layh_s','','')

    ALLOCATE (soilmoist_s(kjpindex, nslm, nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for soilmoist_s','','')

    ALLOCATE(tdeep(kjpindex,ngrnd,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for tdeep','','')

    ALLOCATE(hsdeep(kjpindex,ngrnd,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for hsdeep','','')

    ALLOCATE(heat_Zimov(kjpindex,ngrnd,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for heat_Zimov','','')

    ALLOCATE(sfluxCH4_deep(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for sfluxCH4_deep','','')

    ALLOCATE(sfluxCO2_deep(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for sfluxCO2_deep','','')

    ALLOCATE(som_total(kjpindex,ngrnd,nvm,nelements),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for som_total','','')

    ALLOCATE (altmax(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for altmax','','')
    ALLOCATE(depth_organic_soil(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for depth_organic_soil','','')


    !! 1.5 Irrigation related variables
    ALLOCATE (root_deficit(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for root_deficit','','') !

    ALLOCATE (irrig_frac(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for irrig_frac','','') !
    irrigation(:) = zero

    ALLOCATE (irrigated_next(kjpindex),stat=ier) !
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable irrigated_next','','') !

    ALLOCATE (fraction_aeirrig_sw(kjpindex),stat=ier) !
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for fraction_aeirrig_sw','','')


    !! 1.6 Initialize lake variable allocation

    IF (ok_lake_energy) THEN
       ALLOCATE (frac_lake(kjpindex,nlake),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for frac_lake','','')
       
       ALLOCATE (vevapp_lake(kjpindex,nlake),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for vevapp_lake','','')
       
       ALLOCATE (temp_sol_lake(kjpindex,nlake),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for temp_sol_lake','','')
       
       ALLOCATE (fluxsens_lake(kjpindex,nlake),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for fluxsens_lake','','')
       
       ALLOCATE (fluxlat_lake(kjpindex,nlake),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for fluxlat_lake','','')
    END IF

    !! 1.6 Initialize indexing table for the vegetation fields. 
    ! In SECHIBA we work on reduced grids but to store in the full 3D filed vegetation variable 
    ! we need another index table : indexveg, indexsoil, indexnobio and indexgrnd
    DO ji = 1, kjpindex
       !
       DO jv = 1, nlai+1
          indexlai((jv-1)*kjpindex + ji) = INDEX(ji) + (jv-1)*kjpij + offset_omp - offset_mpi
       ENDDO
       !
       DO jv = 1, nvm
          indexveg((jv-1)*kjpindex + ji) = INDEX(ji) + (jv-1)*kjpij + offset_omp - offset_mpi
       ENDDO
       !      
       DO jv = 1, nstm
          indexsoil((jv-1)*kjpindex + ji) = INDEX(ji) + (jv-1)*kjpij + offset_omp - offset_mpi
       ENDDO
       !      
       DO jv = 1, nnobio
          indexnobio((jv-1)*kjpindex + ji) = INDEX(ji) + (jv-1)*kjpij + offset_omp - offset_mpi
       ENDDO
       !
       DO jv = 1, ngrnd
          indexgrnd((jv-1)*kjpindex + ji) = INDEX(ji) + (jv-1)*kjpij + offset_omp - offset_mpi
       ENDDO
       !
       DO jv = 1, nsnow
          indexsnow((jv-1)*kjpindex + ji) = INDEX(ji) + (jv-1)*kjpij
       ENDDO

       DO jv = 1, nslm
          indexnslm((jv-1)*kjpindex + ji) = INDEX(ji) + (jv-1)*kjpij
       ENDDO

       DO jv = 1, nslm
          indexlayer((jv-1)*kjpindex + ji) = INDEX(ji) + (jv-1)*kjpij + offset_omp - offset_mpi
       ENDDO
       !
       DO jv = 1, 2
          indexalb((jv-1)*kjpindex + ji) = INDEX(ji) + (jv-1)*kjpij + offset_omp - offset_mpi
       ENDDO
       !
    ENDDO

    ! 1.7 now we create and initialize some plant variables that need to be passed
    ! around sechiba and stomate (from the merge, there might be more to
    ! put here)
    ALLOCATE(nlevels_loc(kjpindex),stat=ier)
    IF (ier .NE. 0) CALL ipslerr_p (3,'sechiba_init','Pb in alloc for nlevels_loc','','')
    nlevels_loc(:) = val_exp

    ALLOCATE (z_array_out(kjpindex,nvm,ncirc,nlevels_tot),STAT=ier)
    IF (ier.NE.0) CALL ipslerr_p (3,'sechiba_init','Pb in alloc for z_array_out','','')

    ALLOCATE (profile_vbeta3(kjpindex,nvm,nlevels_tot),STAT=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in profile_vbeta3 allocation.'
       WRITE (numout,*) 'We stop. We need kjpindex*nvm words =',kjpindex*nvm*nlevels_tot
       STOP
    END IF

    ALLOCATE (profile_rveget(kjpindex,nvm,nlevels_tot),STAT=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in profile_rveget allocation.'
       WRITE (numout,*) 'We stop. We need kjpindex*nvm words =',kjpindex*nvm*nlevels_tot
       STOP
    END IF

        ALLOCATE (delta_c13_assim(kjpindex,nvm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in delta_c13_assim allocation. We stop. We need kjpindex x nvm words = ',&
            & kjpindex,' x ' ,nvm, ' = ',kjpindex*nvm
       STOP 'sechiba_init'
    END IF
    delta_c13_assim(:,:)=undef_sechiba

    ALLOCATE (leaf_ci_out(kjpindex,nvm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in leaf_ci_out allocation. We stop. We need kjpindex x nvm words = ',&
            & kjpindex,' x ' ,nvm, ' = ',kjpindex*nvm
       STOP 'sechiba_init'
    END IF
    leaf_ci_out(:,:)=undef_sechiba

    ALLOCATE (gpp_day(kjpindex,nvm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in gpp_day allocation. We stop. We need kjpindex x nvm words = ',&
            & kjpindex,' x ' ,nvm, ' = ',kjpindex*nvm
       STOP 'sechiba_init'
    END IF
    gpp_day(:,:)=undef_sechiba

    ALLOCATE(circ_class_n(kjpindex,nvm,ncirc),stat=ier)
    IF (ier .NE. 0) THEN
       WRITE(numout,*) 'Memory allocation error for circ_class_n. We stop. We need kjpindex*nvm*ncirc words', &
       &      kjpindex,nvm,ncirc
      CALL ipslerr_p (3,'sechiba','sechiba_init','','')
    ENDIF
    circ_class_n(:,:,:) = val_exp

    ALLOCATE(circ_class_biomass(kjpindex,nvm,ncirc,nparts,nelements),stat=ier)
    IF (ier .NE. 0) THEN
       WRITE(numout,*) 'Memory allocation error for circ_class_biomass.'
       WRITE(numout,*) 'We stop. We need kjpindex*nvm*nparts*ncirc*nelmements words', &
       &      kjpindex,nvm,ncirc,nparts,nelements
      CALL ipslerr_p (3,'sechiba','sechiba_init','','')
    ENDIF
    circ_class_biomass(:,:,:,:,:) = val_exp

    ALLOCATE(lai_per_level(kjpindex,nvm,nlevels_tot),stat=ier)
    IF (ier .NE. 0) THEN
       WRITE(numout,*) 'Memory allocation error for lai_per_level. We stop. '//&
            'We need kjpindex*nvm*nlevels_tot words', &
       &      kjpindex,nvm,nlevels_tot
      CALL ipslerr_p (3,'sechiba','sechiba_init','','')
    ENDIF
    lai_per_level(:,:,:)=val_exp

    ALLOCATE (loss_gain(kjpindex,nvm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in loss_gain allocation. We stop. We need kjpindex x nvm words = ',&
            & kjpindex,' x ' ,nvm, ' = ',kjpindex*nvm
       STOP 'sechiba_init'
    END IF
    loss_gain(:,:)=undef_sechiba

    ALLOCATE (veget_max_new(kjpindex,nvm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in veget_max_new allocation. We stop. We need kjpindex x nvm words = ',&
            & kjpindex,' x ' ,nvm, ' = ',kjpindex*nvm
       STOP 'sechiba_init'
    END IF
    veget_max_new(:,:)=undef_sechiba

    ALLOCATE (frac_nobio_new(kjpindex,nnobio),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in frac_nobio_new allocation. We stop. We need kjpindex x nnobio words = ',&
            & kjpindex,' x ' ,nnobio, ' = ',kjpindex*nnobio
       STOP 'sechiba_init'
    END IF
    frac_nobio_new(:,:)=undef_sechiba
    
    ! multilayer allocations
    ALLOCATE (u_speed(kjpindex,jnlvls),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in u_speed allocation. We stop. We need kjpindex words = ',kjpindex,jnlvls
       CALL ipslerr_p (3,'sechiba','sechiba_init','','')
    END IF

    ALLOCATE (swc(kjpindex,nslm,nstm))
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in swc allocation. We stop. We need kjpindex words = ',kjpindex*nslm*nstm
       CALL ipslerr_p (3,'sechiba','sechiba_init','','')
    END IF

    ALLOCATE (ksoil(kjpindex,nslm,nstm))
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in ksoil allocation. We stop. We need kjpindexwords = ',kjpindex*nslm*nstm
       CALL ipslerr_p (3,'sechiba','sechiba_init','','')
    END IF

    ALLOCATE (root_profile(kjpindex,nvm,nslm,nroot_prof))
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in root_profile allocation. We stop. We need kjpindexwords = ',kjpindex*nvm*nslm*nroot_prof 
       CALL ipslerr_p (3,'sechiba','sechiba_init','','')
    END IF

    ALLOCATE (root_depth(kjpindex,nvm,ndepths))
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in root_depth allocation. We stop. We need kjpindexwords = ',kjpindex*nvm*ndepths 
       CALL ipslerr_p (3,'sechiba','sechiba_init','','')
    END IF

    ALLOCATE (coszang_noon(kjpindex))
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in coszang_noon allocation. We stop. We need kjpindexwords = ',kjpindex 
       CALL ipslerr_p (3,'sechiba','sechiba_init','','')
    END IF

    ! Read coszang_noon from restart file
    var_name= 'coszang_noon'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Solar angel at noon')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, &
         .TRUE., coszang_noon, "gather", nbp_glo, index_g)
    IF ( ALL( coszang_noon(:) .EQ. val_exp ) ) THEN
       ! coszang_noon was not found in restart file.
       IF (printlev>=1) WRITE (numout,*) 'coszang_noon was not found in restart file. Call solarang_noon to calculate it.'
       CALL solarang_noon(kjpindex, lalo, coszang_noon)
    ELSEIF ( ANY( ABS(coszang_noon(:)) .GT. 2 ) ) THEN
       ! This is a test if coszang_noon is out of range. This can be the case if
       ! there is a mismatch in the land-sea mask. coszang_noon is the first
       ! restart variable read in ORCHIDEE. 
       CALL ipslerr_p(3, 'sechiba_init', 'coszang_noon from restart file is out of range.',&
                         'There is probably a mismatch in the land/sea mask in the model and in the restart file.',&
                         'Start the model again without restart files for ORCHIDEE.')
    END IF

    ALLOCATE (Pgap_cumul(kjpindex,nvm,nlevels_tot))
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in Pgap_cumul allocation. We stop. We need kjpindexwords = ', kjpindex,nvm,nlevels_tot
       CALL ipslerr_p (3,'sechiba','sechiba_init','','')
    END IF

    ALLOCATE (wtp(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for wtp','','')
    wtp(:) = undef_sechiba

    ALLOCATE (liqwt_ratio(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for liqwt_ratio','','')

    ALLOCATE (mc_peat_above(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for mc_peat_above','','')

    ALLOCATE (mc_croppeat_above(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for mc_croppeat_above','','')

    ALLOCATE (shumdiag_peat(kjpindex, nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for shumdiag_peat','','')

    ALLOCATE (shumdiag_croppeat(kjpindex, nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for shumdiag_croppeat','','')

    ALLOCATE (shumdiag_man(kjpindex, nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for shumdiag_man','','')

    ALLOCATE (mc_man_above(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for mc_man_above','','')


!! 2. Read the default value that will be put into variable which are not in the restart file
    CALL ioget_expval(val_exp)
    
    IF (printlev_loc>=3) WRITE (numout,*) ' sechiba_init done '

  END SUBROUTINE sechiba_init
  

!! ==============================================================================================================================\n
!! SUBROUTINE 	: sechiba_clear
!!
!>\BRIEF        Deallocate memory of sechiba's variables
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): None 
!!
!! REFERENCE(S)	: None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================ 

  SUBROUTINE sechiba_clear()

!! 1. Initialize first run

    l_first_sechiba=.TRUE.

!! 2. Deallocate dynamic variables of sechiba

    IF ( ALLOCATED (indexveg)) DEALLOCATE (indexveg)
    IF ( ALLOCATED (indexlai)) DEALLOCATE (indexlai)
    IF ( ALLOCATED (indexsoil)) DEALLOCATE (indexsoil)
    IF ( ALLOCATED (indexnobio)) DEALLOCATE (indexnobio)
    IF ( ALLOCATED (indexsnow)) DEALLOCATE (indexsnow)
    IF ( ALLOCATED (indexgrnd)) DEALLOCATE (indexgrnd)
    IF ( ALLOCATED (indexlayer)) DEALLOCATE (indexlayer)
    IF ( ALLOCATED (indexnslm)) DEALLOCATE (indexnslm)
    IF ( ALLOCATED (indexalb)) DEALLOCATE (indexalb)
    IF ( ALLOCATED (flood_res)) DEALLOCATE (flood_res)
    IF ( ALLOCATED (flood_frac)) DEALLOCATE (flood_frac)
    IF ( ALLOCATED (snow)) DEALLOCATE (snow)
    IF ( ALLOCATED (snow_age)) DEALLOCATE (snow_age)
    IF ( ALLOCATED (drysoil_frac)) DEALLOCATE (drysoil_frac)
    IF ( ALLOCATED (evap_bare_lim)) DEALLOCATE (evap_bare_lim)
    IF ( ALLOCATED (evap_bare_lim_ns)) DEALLOCATE (evap_bare_lim_ns)
    IF ( ALLOCATED (evapot)) DEALLOCATE (evapot)
    IF ( ALLOCATED (netrad)) DEALLOCATE (netrad)
    IF ( ALLOCATED (lwabs)) DEALLOCATE (lwabs)
    IF ( ALLOCATED (lwnet)) DEALLOCATE (lwnet)
    IF ( ALLOCATED (fluxsubli)) DEALLOCATE (fluxsubli)
    IF ( ALLOCATED (evapot_corr)) DEALLOCATE (evapot_corr)
    IF ( ALLOCATED (humrel)) DEALLOCATE (humrel)
    IF ( ALLOCATED (vessel_loss)) DEALLOCATE (vessel_loss)
    IF ( ALLOCATED (psi_leaf)) DEALLOCATE (psi_leaf)
    IF ( ALLOCATED (psi_leaf_midday)) DEALLOCATE (psi_leaf_midday)
    IF ( ALLOCATED (psi_leaf_next)) DEALLOCATE (psi_leaf_next)
    IF ( ALLOCATED (psi_sto_leaf_save)) DEALLOCATE (psi_sto_leaf_save)
    IF ( ALLOCATED (psi_sto_wood_save)) DEALLOCATE (psi_sto_wood_save)
    IF ( ALLOCATED (psi_root_sup)) DEALLOCATE (psi_root_sup)
    IF ( ALLOCATED (psi_root_inf)) DEALLOCATE (psi_root_inf)
    IF ( ALLOCATED (psi_soil_sup)) DEALLOCATE (psi_soil_sup)
    IF ( ALLOCATED (psi_soil_inf)) DEALLOCATE (psi_soil_inf)
    IF ( ALLOCATED (psi_xylem_trunk)) DEALLOCATE (psi_xylem_trunk)
    IF ( ALLOCATED (psi_xylem_leaf)) DEALLOCATE (psi_xylem_leaf)
    IF ( ALLOCATED (psi_xylem_collar)) DEALLOCATE (psi_xylem_collar)
    IF ( ALLOCATED (psi_sto_wood)) DEALLOCATE (psi_sto_wood)
    IF ( ALLOCATED (psi_sto_leaf)) DEALLOCATE (psi_sto_leaf)
    IF ( ALLOCATED (mc_i_sup)) DEALLOCATE (mc_i_sup)
    IF ( ALLOCATED (mc_i_inf)) DEALLOCATE (mc_i_inf)
    IF ( ALLOCATED (F_absorption)) DEALLOCATE (F_absorption)
    IF ( ALLOCATED (vegstress)) DEALLOCATE (vegstress)
    IF ( ALLOCATED (soiltile)) DEALLOCATE (soiltile)
    IF ( ALLOCATED (fraclut)) DEALLOCATE (fraclut)
    IF ( ALLOCATED (nwdFraclut)) DEALLOCATE (nwdFraclut)
    IF ( ALLOCATED (njsc)) DEALLOCATE (njsc)
    IF ( ALLOCATED (reinf_slope)) DEALLOCATE (reinf_slope)
    IF ( ALLOCATED (reinf_slope_soil)) DEALLOCATE (reinf_slope_soil)
    IF ( ALLOCATED (ks)) DEALLOCATE (ks)
    IF ( ALLOCATED (nvan)) DEALLOCATE (nvan)
    IF ( ALLOCATED (avan)) DEALLOCATE (avan)
    IF ( ALLOCATED (mcr)) DEALLOCATE (mcr)
    IF ( ALLOCATED (mcs)) DEALLOCATE (mcs)
    IF ( ALLOCATED (mcfc)) DEALLOCATE (mcfc)
    IF ( ALLOCATED (mcw)) DEALLOCATE (mcw)
    IF ( ALLOCATED (us)) DEALLOCATE (us)
    IF ( ALLOCATED (frac_sub_snow)) DEALLOCATE (frac_sub_snow)
    IF ( ALLOCATED (frac_evap_bare)) DEALLOCATE (frac_evap_bare)
    IF ( ALLOCATED (frac_evap_flood)) DEALLOCATE (frac_evap_flood)
    IF ( ALLOCATED (soilcap)) DEALLOCATE (soilcap)
    IF ( ALLOCATED (soilflx)) DEALLOCATE (soilflx)
    IF ( ALLOCATED (temp_sol)) DEALLOCATE (temp_sol)
    IF ( ALLOCATED (qsurf)) DEALLOCATE (qsurf)
    IF ( ALLOCATED (qsintveg)) DEALLOCATE (qsintveg)
    IF ( ALLOCATED (frac_evap_intercept))  DEALLOCATE (frac_evap_intercept)
    IF ( ALLOCATED (resist_intercept))  DEALLOCATE (resist_intercept)
    IF ( ALLOCATED (frac_evap_transp))  DEALLOCATE (frac_evap_transp)
    IF ( ALLOCATED (resist_transp))  DEALLOCATE (resist_transp)
    IF ( ALLOCATED (vbeta3pot)) DEALLOCATE (vbeta3pot)
    IF ( ALLOCATED (gsmean)) DEALLOCATE (gsmean)
    IF ( ALLOCATED (cimean)) DEALLOCATE (cimean)
    IF ( ALLOCATED (assimtot)) DEALLOCATE (assimtot)
    IF ( ALLOCATED (Rdtot_pft)) DEALLOCATE (Rdtot_pft)
    IF ( ALLOCATED (gpp)) DEALLOCATE (gpp)
    IF ( ALLOCATED (temp_growth)) DEALLOCATE (temp_growth) 
    IF ( ALLOCATED (veget)) DEALLOCATE (veget)
    IF ( ALLOCATED (veget_max)) DEALLOCATE (veget_max)
    IF ( ALLOCATED (tot_bare_soil)) DEALLOCATE (tot_bare_soil)
    IF ( ALLOCATED (lai)) DEALLOCATE (lai)
    IF ( ALLOCATED (frac_age)) DEALLOCATE (frac_age)
    IF ( ALLOCATED (height)) DEALLOCATE (height)
    IF ( ALLOCATED (height_dom)) DEALLOCATE (height_dom)
    IF ( ALLOCATED (dia_dom)) DEALLOCATE (dia_dom)
    IF ( ALLOCATED (height_inv)) DEALLOCATE (height_inv)
    IF ( ALLOCATED (dia_inv)) DEALLOCATE (dia_inv)
    IF ( ALLOCATED (ba_inv)) DEALLOCATE (ba_inv)
    IF ( ALLOCATED (ind_inv)) DEALLOCATE (ind_inv)
    IF ( ALLOCATED (roughheight)) DEALLOCATE (roughheight)
    IF ( ALLOCATED (frac_nobio)) DEALLOCATE (frac_nobio)
    IF ( ALLOCATED (snow_nobio_age)) DEALLOCATE (snow_nobio_age)
    IF ( ALLOCATED (assim_param)) DEALLOCATE (assim_param)
    IF ( ALLOCATED (vevapflo)) DEALLOCATE (vevapflo)
    IF ( ALLOCATED (vevapsno)) DEALLOCATE (vevapsno)
    IF ( ALLOCATED (vevapnu)) DEALLOCATE (vevapnu)
    IF ( ALLOCATED (totfrac_nobio)) DEALLOCATE (totfrac_nobio)
    IF ( ALLOCATED (floodout)) DEALLOCATE (floodout)
    IF ( ALLOCATED (runoff)) DEALLOCATE (runoff)
    IF ( ALLOCATED (drainage)) DEALLOCATE (drainage)
    IF ( ALLOCATED (reinfiltration)) DEALLOCATE (reinfiltration)
    IF ( ALLOCATED (irrigation)) DEALLOCATE (irrigation)
    IF ( ALLOCATED (tot_melt)) DEALLOCATE (tot_melt)
    IF ( ALLOCATED (resist_equiv_evap)) DEALLOCATE (resist_equiv_evap)
    IF ( ALLOCATED (resist_equiv_evap_veg)) DEALLOCATE (resist_equiv_evap_veg)
    IF ( ALLOCATED (resist_equiv_tot)) DEALLOCATE (resist_equiv_tot)
    IF ( ALLOCATED (raero)) DEALLOCATE (raero)
    IF ( ALLOCATED (rau)) DEALLOCATE (rau)
    IF ( ALLOCATED (deadleaf_cover)) DEALLOCATE (deadleaf_cover)
    IF ( ALLOCATED (stempdiag)) DEALLOCATE (stempdiag)
    IF ( ALLOCATED (co2_to_bm)) DEALLOCATE (co2_to_bm)
    IF ( ALLOCATED (ftempdiag)) DEALLOCATE (ftempdiag)
    IF ( ALLOCATED (shumdiag)) DEALLOCATE (shumdiag)
    IF ( ALLOCATED (shumdiag_perma)) DEALLOCATE (shumdiag_perma)
    IF ( ALLOCATED (litterhumdiag)) DEALLOCATE (litterhumdiag)
    IF ( ALLOCATED (ptnlev1)) DEALLOCATE (ptnlev1)
    IF ( ALLOCATED (k_litt)) DEALLOCATE (k_litt)
    IF ( ALLOCATED (vevapwet)) DEALLOCATE (vevapwet)
    IF ( ALLOCATED (transpir)) DEALLOCATE (transpir)
    IF ( ALLOCATED (stressed)) DEALLOCATE (stressed)
    IF ( ALLOCATED (unstressed)) DEALLOCATE (unstressed)
    IF ( ALLOCATED (transpir_mod)) DEALLOCATE (transpir_mod)
    IF ( ALLOCATED (transpir_supply)) DEALLOCATE (transpir_supply)
    IF ( ALLOCATED (vir_transpir_supply)) DEALLOCATE (vir_transpir_supply)
    IF ( ALLOCATED (e_frac)) DEALLOCATE (e_frac)
    IF ( ALLOCATED (transpot)) DEALLOCATE (transpot)
    IF ( ALLOCATED (qsintmax)) DEALLOCATE (qsintmax)
    IF ( ALLOCATED (rveget)) DEALLOCATE (rveget)
    IF ( ALLOCATED (rstruct)) DEALLOCATE (rstruct)
    IF ( ALLOCATED (frac_snow_veg)) DEALLOCATE (frac_snow_veg)
    IF ( ALLOCATED (frac_snow_nobio)) DEALLOCATE (frac_snow_nobio)
    IF ( ALLOCATED (snowrho)) DEALLOCATE (snowrho)
    IF ( ALLOCATED (ice_sheet_mask)) DEALLOCATE (ice_sheet_mask)
    IF ( ALLOCATED (snowgrain)) DEALLOCATE (snowgrain)
    IF ( ALLOCATED (snowtemp)) DEALLOCATE (snowtemp)
    IF ( ALLOCATED (snowdz)) DEALLOCATE (snowdz)
    IF ( ALLOCATED (snowheat)) DEALLOCATE (snowheat)
    IF ( ALLOCATED (cgrnd_snow)) DEALLOCATE (cgrnd_snow)
    IF ( ALLOCATED (dgrnd_snow)) DEALLOCATE (dgrnd_snow)
    IF ( ALLOCATED (lambda_snow)) DEALLOCATE(lambda_snow)
    IF ( ALLOCATED (temp_sol_add)) DEALLOCATE(temp_sol_add)
    IF ( ALLOCATED (qsol_sat_new)) DEALLOCATE(qsol_sat_new)
    IF ( ALLOCATED (qair_new)) DEALLOCATE(qair_new)
    IF ( ALLOCATED (gtemp)) DEALLOCATE (gtemp)
    IF ( ALLOCATED (pgflux)) DEALLOCATE (pgflux)
    IF ( ALLOCATED (u_speed)) DEALLOCATE(u_speed)
    IF ( ALLOCATED (warnings)) DEALLOCATE (warnings)
    IF ( ALLOCATED (lai_per_level)) DEALLOCATE(lai_per_level)
    IF ( ALLOCATED (laieff_fit)) DEALLOCATE(laieff_fit)
    IF ( ALLOCATED (max_height_store)) DEALLOCATE(max_height_store)
    IF ( ALLOCATED (z_array_out)) DEALLOCATE(z_array_out)
    IF ( ALLOCATED (Light_Abs_Tot)) DEALLOCATE(Light_Abs_Tot)
    IF ( ALLOCATED (Light_Tran_Tot)) DEALLOCATE(Light_Tran_Tot)
    IF ( ALLOCATED (Light_Abs_Tot_mean)) DEALLOCATE(Light_Abs_Tot_mean)
    IF ( ALLOCATED (Light_Alb_Tot_mean)) DEALLOCATE(Light_Alb_Tot_mean)
    IF ( ALLOCATED (laieff_isotrop)) DEALLOCATE(laieff_isotrop)
    IF ( ALLOCATED (circ_class_biomass)) DEALLOCATE(circ_class_biomass)
    IF ( ALLOCATED (circ_class_n)) DEALLOCATE(circ_class_n)
    IF ( ALLOCATED (loss_gain)) DEALLOCATE(loss_gain)
    IF ( ALLOCATED (veget_max_new)) DEALLOCATE(veget_max_new)
    IF ( ALLOCATED (frac_nobio_new)) DEALLOCATE(frac_nobio_new)
    IF ( ALLOCATED (mc_layh)) DEALLOCATE (mc_layh)
    IF ( ALLOCATED (mcl_layh)) DEALLOCATE (mcl_layh)
    IF ( ALLOCATED (soilmoist)) DEALLOCATE (soilmoist)
    IF ( ALLOCATED (profile_vbeta3)) DEALLOCATE(profile_vbeta3)
    IF ( ALLOCATED (profile_rveget)) DEALLOCATE(profile_rveget)
    IF ( ALLOCATED (delta_c13_assim)) DEALLOCATE (delta_c13_assim)
    IF ( ALLOCATED (leaf_ci_out)) DEALLOCATE (leaf_ci_out)
    IF ( ALLOCATED (mcl_layh_s)) DEALLOCATE (mc_layh_s)
    IF ( ALLOCATED (mc_layh_s)) DEALLOCATE (mc_layh_s)
    IF ( ALLOCATED (soilmoist_s)) DEALLOCATE (soilmoist_s)
    IF ( ALLOCATED (tdeep)) DEALLOCATE (tdeep)
    IF ( ALLOCATED (hsdeep)) DEALLOCATE (hsdeep)
    IF ( ALLOCATED (heat_Zimov)) DEALLOCATE (heat_Zimov)
    IF ( ALLOCATED (sfluxCH4_deep)) DEALLOCATE (sfluxCH4_deep)
    IF ( ALLOCATED (sfluxCO2_deep)) DEALLOCATE (sfluxCO2_deep)
    IF ( ALLOCATED (som_total)) DEALLOCATE (som_total)
    IF ( ALLOCATED (altmax)) DEALLOCATE (altmax)
    IF ( ALLOCATED (swc)) DEALLOCATE(swc)
    IF ( ALLOCATED (ksoil)) DEALLOCATE(ksoil)
    IF ( ALLOCATED (root_profile)) DEALLOCATE (root_profile)
    IF ( ALLOCATED (coszang_noon)) DEALLOCATE (coszang_noon)
    IF ( ALLOCATED (icetemp)) DEALLOCATE (icetemp)
    IF ( ALLOCATED (icedz)) DEALLOCATE (icedz)
    IF ( ALLOCATED (cgrnd_ice)) DEALLOCATE (cgrnd_ice)
    IF ( ALLOCATED (dgrnd_ice)) DEALLOCATE (dgrnd_ice)
    IF ( ALLOCATED (lambda_ice)) DEALLOCATE(lambda_ice)
    IF ( ALLOCATED (root_deficit)) DEALLOCATE (root_deficit)
    IF ( ALLOCATED (irrig_frac)) DEALLOCATE (irrig_frac)
    IF ( ALLOCATED (irrigated_next)) DEALLOCATE (irrigated_next)
    IF ( ALLOCATED (precip_longterm)) DEALLOCATE (precip_longterm)
    
    ! peat
    IF ( ALLOCATED (wtp)) DEALLOCATE (wtp)
    IF ( ALLOCATED (liqwt_ratio)) DEALLOCATE (liqwt_ratio)
    IF ( ALLOCATED (mc_peat_above)) DEALLOCATE (mc_peat_above)
    IF ( ALLOCATED (mc_croppeat_above)) DEALLOCATE (mc_croppeat_above)
    IF ( ALLOCATED (shumdiag_peat)) DEALLOCATE (shumdiag_peat)
    IF ( ALLOCATED (shumdiag_croppeat)) DEALLOCATE (shumdiag_croppeat)
    IF ( ALLOCATED (shumdiag_man)) DEALLOCATE (shumdiag_man)
    IF ( ALLOCATED (mc_man_above)) DEALLOCATE (mc_man_above)

    ! lake
    IF ( ALLOCATED (frac_lake)) DEALLOCATE (frac_lake)
    IF ( ALLOCATED (vevapp_lake)) DEALLOCATE (vevapp_lake)
    IF ( ALLOCATED (temp_sol_lake)) DEALLOCATE (temp_sol_lake)
    IF ( ALLOCATED (fluxsens_lake)) DEALLOCATE (fluxsens_lake)
    IF ( ALLOCATED (fluxlat_lake)) DEALLOCATE (fluxlat_lake)

!! 3. Clear all allocated memory

    CALL pft_parameters_clear
    CALL dynamic_parameters_clear
    CALL slowproc_clear 
    CALL diffuco_clear 
    IF (ok_mleb) THEN
       CALL mleb_clear
    ELSE
       CALL enerbil_clear  
    ENDIF
    CALL lake_clear
    CALL hydrol_clear 
    CALL thermosoil_clear
    CALL condveg_clear 
    CALL routing_wrapper_clear

  END SUBROUTINE sechiba_clear


!! ==============================================================================================================================\n
!! SUBROUTINE 	: sechiba_var_init
!!
!>\BRIEF        Calculate air density as a function of air temperature and 
!! pressure for each terrestrial pixel.
!! 
!! RECENT CHANGE(S): None
!! 
!! MAIN OUTPUT VARIABLE(S): air density (::rau, kg m^{-3}).
!! 
!! REFERENCE(S)	: None
!! 
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE sechiba_var_init (kjpindex, rau, pb, temp_air) 

!! 0.1 Input variables

    INTEGER(i_std), INTENT (in)                    :: kjpindex        !! Domain size - terrestrial pixels only (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)  :: pb              !! Surface pressure (hPa)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)  :: temp_air        !! Air temperature (K)
    
!! 0.2 Output variables

    REAL(r_std),DIMENSION (kjpindex), INTENT (out) :: rau             !! Air density @tex $(kg m^{-3})$ @endtex

!! 0.3 Modified variables

!! 0.4 Local variables

    INTEGER(i_std)                                 :: ji              !! Indices (unitless)
!_ ================================================================================================================================
    
!! 1. Calculate intial air density (::rau)
   
    DO ji = 1,kjpindex
       rau(ji) = pa_par_hpa * pb(ji) / (cte_molr*temp_air(ji))
    END DO

    IF (printlev_loc>=3) WRITE (numout,*) ' sechiba_var_init done '

  END SUBROUTINE sechiba_var_init


!! ==============================================================================================================================\n
!! SUBROUTINE 	: sechiba_end
!!
!>\BRIEF        Swap old for newly calculated soil temperature.
!! 
!! RECENT CHANGE(S): None
!! 
!! MAIN OUTPUT VARIABLE(S): soil temperature (::temp_sol; K)
!! 
!! REFERENCE(S)	: None
!! 
!! FLOWCHART    : None
!! \n
!! ================================================================================================================================ 

  SUBROUTINE sechiba_end (kjpindex, temp_sol_new, temp_sol)
                         

!! 0.1 Input variables

    INTEGER(i_std), INTENT (in)                       :: kjpindex           !! Domain size - terrestrial pixels only (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)     :: temp_sol_new       !! New surface temperature (K)
    
    !! 0.2 Output variables
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)    :: temp_sol           !! Surface temperature (K)

!_ ================================================================================================================================
    
!! 1. Swap temperature

    temp_sol(:) = temp_sol_new(:)
    
    IF (printlev_loc>=3) WRITE (numout,*) ' sechiba_end done '

  END SUBROUTINE sechiba_end


!! ==============================================================================================================================\n
!! SUBROUTINE 	: sechiba_interface_orchidee_inca
!!
!>\BRIEF        make the interface between surface and atmospheric chemistry
!! 
!! DESCRIPTION  : This subroutine is called from INCA, the atmospheric chemistry model. It is used to transfer variables from ORCHIDEE to INCA. 
!!
!! RECENT CHANGE(S): move from chemistry module to be more generic (feb - 2017)
!! 
!! MAIN OUTPUT VARIABLE(S): emission COV to be transport by orchidee to inca in fields_out array 
!! 
!! REFERENCE(S)	: None
!! 
!! FLOWCHART    : None
!! \n
!! ================================================================================================================================ 
  SUBROUTINE sechiba_interface_orchidee_inca( &
       nvm_out, veget_max_out, veget_frac_out, lai_out, snow_out, &
       field_out_COV_names, fields_out_COV, field_in_COV_names, fields_in_COV, &
       field_out_Nsoil_names, fields_out_Nsoil, nnspec_out)


    INTEGER, INTENT(out)                      :: nvm_out            !! Number of vegetation types
    REAL(r_std), DIMENSION (:,:), INTENT(out) :: veget_max_out      !! Max. fraction of vegetation type (LAI -> infty)
    REAL(r_std), DIMENSION (:,:), INTENT(out) :: veget_frac_out     !! Fraction of vegetation type (unitless, 0-1)  
    REAL(r_std), DIMENSION (:,:), INTENT(out) :: lai_out            !! Surface foliere
    REAL(r_std), DIMENSION (:)  , INTENT(out) :: snow_out           !! Snow mass [Kg/m^2]

    !
    ! Optional arguments
    !
    ! Names and fields for emission variables : to be transport by Orchidee to Inca
    CHARACTER(LEN=*),DIMENSION(:), OPTIONAL, INTENT(IN) :: field_out_COV_names
    REAL(r_std),DIMENSION(:,:,:), OPTIONAL, INTENT(OUT) :: fields_out_COV
    !
    ! Names and fields for deposit variables : to be transport from chemistry model by INCA to ORCHIDEE.
    CHARACTER(LEN=*),DIMENSION(:), OPTIONAL, INTENT(IN) :: field_in_COV_names
    REAL(r_std),DIMENSION(:,:), OPTIONAL, INTENT(IN)    :: fields_in_COV

    CHARACTER(LEN=*),DIMENSION(:), OPTIONAL, INTENT(IN) :: field_out_Nsoil_names  !! Not used before ORCHIDEE_3
    REAL(r_std),DIMENSION(:,:,:), OPTIONAL, INTENT(out) :: fields_out_Nsoil       !! Not used before ORCHIDEE_3
    INTEGER, OPTIONAL, INTENT(out)                      :: nnspec_out             !! Not used before ORCHIDEE_3


    IF (PRESENT(field_out_Nsoil_names) .OR. PRESENT(fields_out_Nsoil) .OR. PRESENT(nnspec_out)) THEN
       CALL ipslerr_p(3,'sechiba_interface_orchidee_inca','This version of Orchidee is not usable for coupling nitrogen with atmosphere',&
            'Please use Orchidee_3 or modify coupling between atm and surf','')
    ENDIF

    ! Variables always transmitted from sechiba to inca
    nvm_out = nvm 
    veget_max_out(:,:)  = veget_max(:,:) 
    veget_frac_out(:,:) = veget(:,:) 
    lai_out(:,:)  = lai(:,:) 
    snow_out(:)  = snow(:) 

    ! Call chemistry_flux_interface if at least one of variables field_out_names or
    ! field_in_names is present in the argument list of sechiba_interface_orchidee_inca when called from inca.
    IF (PRESENT(field_out_COV_names) .AND. .NOT. PRESENT(field_in_COV_names)) THEN 
       CALL chemistry_flux_interface(field_out_names=field_out_COV_names, fields_out=fields_out_COV)
    ELSE IF (.NOT. PRESENT(field_out_COV_names) .AND. PRESENT(field_in_COV_names)) THEN 
       CALL chemistry_flux_interface(field_in_names=field_in_COV_names, fields_in=fields_in_COV)
    ELSE IF (PRESENT(field_out_COV_names) .AND. PRESENT(field_in_COV_names)) THEN 
       CALL chemistry_flux_interface(field_out_names=field_out_COV_names, fields_out=fields_out_COV, &
            field_in_names=field_in_COV_names, fields_in=fields_in_COV)
    ENDIF

  END SUBROUTINE sechiba_interface_orchidee_inca



 !! ==============================================================================================================================\n
 !! SUBROUTINE 	: sechiba_energy_merge_diag
 !!
 !! 
 !! DESCRIPTION  : This subroutine will merge in the grid cell the fluxes on the lake fractions with the vegetation and
 !!                bare soil fractions. The results will be written to output diagnostic file.
 !!
 !! RECENT CHANGE(S): 
 !! 
 !! MAIN OUTPUT VARIABLE(S):  
 !! 
 !! REFERENCE(S)	: None
 !! 
 !! FLOWCHART    : None
 !! \n
 !! ================================================================================================================================ 
  SUBROUTINE sechiba_energy_merge_diag(kjpindex, temp_sol_lake, fluxsens_lake, fluxlat_lake, &
       temp_sol, fluxsens, fluxlat, frac_lake)
    
    !! 0. Variable and parameter declaration
    !! 0.1 Input variables
    
    INTEGER(i_std), INTENT (in)                          ::  kjpindex
    REAL(r_std),DIMENSION (kjpindex,nlake), INTENT (in)  :: temp_sol_lake
    REAL(r_std),DIMENSION (kjpindex,nlake), INTENT (in)  :: fluxsens_lake
    REAL(r_std),DIMENSION (kjpindex,nlake), INTENT (in)  :: fluxlat_lake
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: temp_sol
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: fluxsens
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: fluxlat
    REAL(r_std),DIMENSION (kjpindex,nlake), INTENT (in)  :: frac_lake
    
    
    !! 0.2 Local variables
    REAL(r_std),DIMENSION (kjpindex) :: temp_sol_grid
    REAL(r_std),DIMENSION (kjpindex) :: fluxsens_grid
    REAL(r_std),DIMENSION (kjpindex) :: fluxlat_grid  
    REAL(r_std),DIMENSION (kjpindex) :: frac_nolake_tot
    REAL(r_std)                      :: frac_lake_tot
    INTEGER(i_std)                   :: ij, ilake



    !! 1. Calcul on the tile
   
    DO ij = 1, kjpindex
       
       frac_lake_tot = SUM(frac_lake(ij,1:nlake))
       
       IF (frac_lake_tot >0) THEN
          temp_sol_grid(ij) = (1-frac_lake_tot) * temp_sol(ij)**4.
          fluxsens_grid(ij) = (1-frac_lake_tot) * fluxsens(ij)
          fluxlat_grid(ij)  = (1-frac_lake_tot) * fluxlat(ij)
          frac_nolake_tot(ij) = (1-frac_lake_tot)
          
          DO ilake = 1, nlake
             IF (frac_lake(ij, ilake)>0) THEN
                temp_sol_grid(ij) = temp_sol_grid(ij) + frac_lake(ij, ilake) * temp_sol_lake(ij, ilake)**4. 
                fluxsens_grid(ij) = fluxsens_grid(ij) + frac_lake(ij, ilake) * fluxsens_lake(ij, ilake) 
                fluxlat_grid(ij)  = fluxlat_grid(ij) + frac_lake(ij, ilake) * fluxlat_lake(ij, ilake)
                frac_nolake_tot = frac_nolake_tot(ij) + frac_lake(ij, ilake) 
             ENDIF
          ENDDO
          
          temp_sol_grid(ij) = temp_sol_grid(ij)**0.25 
       ELSE
          frac_nolake_tot = 1
          temp_sol_grid(ij) = temp_sol(ij)
          fluxsens_grid(ij) = fluxsens(ij)
          fluxlat_grid(ij) = fluxlat(ij)
       ENDIF
       
    ENDDO
    
    CALL xios_orchidee_send_field("frac_nolake", frac_nolake_tot)
    CALL xios_orchidee_send_field("temp_sol_grid",temp_sol_grid)
    CALL xios_orchidee_send_field("fluxsens_grid",fluxsens_grid)
    CALL xios_orchidee_send_field("fluxlat_grid",fluxlat_grid)
    
  END SUBROUTINE sechiba_energy_merge_diag
  
 !! ==============================================================================================================================\n
 !! SUBROUTINE 	: rain_exclusion
 !!
 !! 
 !! DESCRIPTION  : This subroutine allows to quickly play with rain (cut rain,
 !!                double the rain) 
 !!
 !! RECENT CHANGE(S): 
 !! 
 !! MAIN OUTPUT VARIABLE(S):  
 !! 
 !! REFERENCE(S)	: None
 !! 
 !! FLOWCHART    : None
 !! \n
 !! ================================================================================================================================ 
  SUBROUTINE rain_exclusion(precip_rain)
    
    !! 0. Variable and parameter declaration
    !! 0.1 Modified variables
    
    REAL(r_std),DIMENSION(:),INTENT(inout)      :: precip_rain

    !! 1. Modifies precipitation based on a parameter (default 1)
  
    precip_rain(:) = precip_rain(:) * rain_fac
 
  END SUBROUTINE rain_exclusion
  
  
END MODULE sechiba

