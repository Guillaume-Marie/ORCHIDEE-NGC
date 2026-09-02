! =================================================================================================================================
! MODULE 	: constantes_soil_var
!
! CONTACT       : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE       : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF         "constantes_soil_var" module contains the parameters related to soil and hydrology.
!!
!!\n DESCRIPTION : The non saturated hydraulic properties are defined from the  
!!                 formulations of van Genuchten (1980) and Mualem (1976), combined as  
!!                 explained in d'Orgeval (2006). \n
!!                 The related parameters for main soil textures (coarse, medium and fine if "fao", 
!!                 12 USDA testures if "usda") come from Carsel and Parrish (1988).
!!
!! RECENT CHANGE(S): AD: mcw and mcf depend now on soil texture, based on Van Genuchten equations 
!!                   and classical matric potential values, and pcent is adapted
!!                   November 2020 by Salma Tafasca and Agnes Ducharne : we introduce a new texture class 
!!                   for clay oxisols (cf. Tafasca, 2020, PhD thesis; Tafasca et al., in prep for GRL). 
!!                   It makes no change if we read a soil texture map with only 12 USDA classes.
!!                   Lookup tables for Zobler replaces by pointer to read the corresponding values in the
!!                   13-value USDA tables
!!
!! REFERENCE(S)	:
!!- Roger A.Pielke, (2002), Mesoscale meteorological modeling, Academic Press Inc. 
!!- Polcher, J., Laval, K., Dümenil, L., Lean, J., et Rowntree, P. R. (1996).
!! Comparing three land surface schemes used in general circulation models. Journal of Hydrology, 180(1-4), 373--394.
!!- Ducharne, A., Laval, K., et Polcher, J. (1998). Sensitivity of the hydrological cycle
!! to the parametrization of soil hydrology in a GCM. Climate Dynamics, 14, 307--327. 
!!- Rosnay, P. de et Polcher, J. (1999). Modelling root water uptake in a complex land surface
!! scheme coupled to a GCM. Hydrol. Earth Syst. Sci., 2(2/3), 239--255.
!!- d'Orgeval, T. et Polcher, J. (2008). Impacts of precipitation events and land-use changes
!! on West African river discharges during the years 1951--2000. Climate Dynamics, 31(2), 249--262. 
!!- Carsel, R. and Parrish, R.: Developing joint probability distributions of soil water
!! retention characteristics, Water Resour. Res.,24, 755–769, 1988.
!!- Mualem Y (1976). A new model for predicting the hydraulic conductivity  
!! of unsaturated porous media. Water Resources Research 12(3):513-522
!!- Van Genuchten M (1980). A closed-form equation for predicting the  
!! hydraulic conductivity of unsaturated soils. Soil Sci Soc Am J, 44(5):892-898
!!- Tafasca S. (2020). Evaluation de l impact des proprietes du sol sur l hydrologie simulee dans le 
!! modele ORCHIDEE, PhD thesis, Sorbonne Universite.  
!!- Tafasca S., Ducharne A. and Valentin C. Accounting for soil structure in pedo-transfer functions: 
!!  swelling vs non swelling clays. In prep for GRL.
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_parameters/constantes_soil_var.f90 $
!! $Date: 2026-04-06 15:26:48 +0200 (lun. 06 avril 2026) $
!! $Revision: 9464 $
!! \n
!_ ================================================================================================================================

MODULE constantes_soil_var

  USE defprec
  USE vertical_soil_var

  IMPLICIT NONE

  LOGICAL, SAVE             :: check_cwrr               !! Calculate diagnostics to check the water balance in hydrol (true/false)
!$OMP THREADPRIVATE(check_cwrr)

  !! Number of soil classes

  INTEGER(i_std), PARAMETER :: ntext=3                  !! Number of soil textures (Silt, Sand, Clay)
  INTEGER(i_std), SAVE      :: nstm=3                   !! Number of soil tiles (unitless). nstm will be changed only for case OK_PEAT_HYDRO(nstm=4), 
                                                        !! TIDES(nstm=6) else nstm=3
!$OMP THREADPRIVATE(nstm)
  
  INTEGER(i_std), PARAMETER :: lim_layer=9              !! Index of the layer that limitate the superficial of the inferior layers (ok_hydrol_arch_muff)
  INTEGER(i_std), PARAMETER :: nrp=15                   !! Number of nodes when solving Richards equation within the muffs (ok_hydrol_arch_muff) (unitless)
  CHARACTER(LEN=30)         :: soil_classif             !! Type of classification used for the map of soil types.
                                                        !! It must be consistent with soil file given by 
                                                        !! SOILCLASS_FILE parameter.
!$OMP THREADPRIVATE(soil_classif)
  INTEGER(i_std), PARAMETER :: nscm_fao=3               !! For FAO Classification (unitless)
  INTEGER(i_std), PARAMETER :: nscm_usda=13             !! For USDA Classification (unitless)
  INTEGER(i_std), SAVE :: nscm                          !! Default value for nscm
!$OMP THREADPRIVATE(nscm)

  !! Soil types for stomate_windthrow
  INTEGER(i_std), PARAMETER :: n_soil_types = 4         !! Total number of soil types used in
                                                        !! windthrow according to
                                                        !! referencer data
  INTEGER(i_std), PARAMETER :: ifree_draining = 1       !! Free-draining mineral soil
  INTEGER(i_std), PARAMETER :: igleyed = 2              !! Gleyed mineral soil
  INTEGER(i_std), PARAMETER :: ipeaty = 3               !! Peaty mineral soil
  INTEGER(i_std), PARAMETER :: ipeat = 4                !! Deep peat
  INTEGER(i_std), PARAMETER :: n_soil_depths = 3        !! Total number of soil depths used in
                                                        !! windthrow according to
                                                        !! reference data
  INTEGER(i_std), PARAMETER :: ideep = 1                !! Deep-rooted
  INTEGER(i_std), PARAMETER :: ishallow = 2             !! Shallow-rooted
  INTEGER(i_std), PARAMETER :: iaverage = 3             !! Average-rooted

  !! Parameters for soil thermodynamics
  REAL(r_std), SAVE :: sn_cond                          !! Thermal Conductivity of snow 
                                                        !! @tex $(W.m^{-2}.K^{-1})$ @endtex  
!$OMP THREADPRIVATE(sn_cond)
  REAL(r_std), SAVE :: sn_dens                          !! Snow density for the soil thermodynamics
                                                        !! (kg/m3)
!$OMP THREADPRIVATE(sn_dens)
  REAL(r_std), SAVE :: sn_capa                          !! Volumetric heat capacity for snow 
                                                        !! @tex $(J.m^{-3}.K^{-1})$ @endtex
!$OMP THREADPRIVATE(sn_capa)
  REAL(r_std), PARAMETER :: capa_ice = 2.228*1.E3       !! Specific heat capacity of ice (J/kg/K)

  REAL(r_std), PARAMETER :: poros_org = 0.92            !! Organic soil porosity [m3/m3] it is just a number from Dmitry's code
                                                        !! but it is consistent with range given by Rezanezhad et al., 2016
                                                        !! Chem. Geol. [0.71 - 0.951]
  REAL(r_std), PARAMETER :: poros_org_mini = 0.71       !! Minimum possible organic soil porosity [m3/m3], from range given by Rezanezhad et al., 2016  
                                                        !! Chem. Geol. 
                                                        !! Also consistent with Robinson et al., 2022, https://doi.org/10.1038/s41598-022-11099-7 
  REAL(r_std), PARAMETER :: cond_solid_org = 0.25       !! W/m/K from Farouki via Lawrence and Slater
  REAL(r_std), PARAMETER :: cond_dry_org = 0.05         !! W/m/K from Lawrence and Slater
  
  REAL(r_std), PARAMETER :: thkmoss_dry = 0.05          !! Dry moss thermal conductivity (W/m/K) estimated from O'Donnel et al., 2009  
  REAL(r_std), PARAMETER :: thksat_moss_liq = 0.56      !! Saturated liquid moss thermal conductivity (W/M/K) estimated from O'Donnel et al., 2009 
  REAL(r_std), PARAMETER :: thksat_moss_ice = 1.4       !! Saturated ice moss thermal conductivity (W/M/K) estimated from Porada et al., 2016 (TO BE CHECKED)
  
  REAL(r_std), PARAMETER :: so_capa_dry_org = 2.5e6     !! J/K/m^3 from Farouki via Lawrence and Slater

  REAL(r_std), PARAMETER :: pcapa_moss_dry = 0.29e6     !! Dry moss thermal capacity (J/K/m^3) from Soudzilovskaia et al., 2013
  REAL(r_std), PARAMETER :: pcapa_moss_liq = 4.29e6     !! Saturated liquid moss thermal capacity (J/K/m^3) from Soudzilovskaia et al., 2013
  REAL(r_std), PARAMETER :: pcapa_moss_ice = 3.26e6     !! Saturated ice moss thermal conductivity (J/K/m^3) from Druel et al., 2017
  REAL(r_std), PARAMETER :: water_capa = 4.18e+6        !! Volumetric water heat capacity 
                                                        !! @tex $(J.m^{-3}.K^{-1})$ @endtex
  REAL(r_std), PARAMETER :: air_capa = 1.25e+3          !! Volumetric air heat capacity, Hillel 1982 
                                                        !! @tex $(J.m^{-3}.K^{-1})$ @endtex 
  REAL(r_std), PARAMETER :: brk_capa = 2.0e+6           !! Volumetric heat capacity of generic rock
                                                        !! @tex $(J.m^{-3}.K^{-1})$ @endtex
  REAL(r_std), PARAMETER :: brk_cond = 3.0              !! Thermal conductivity of saturated granitic rock
                                                        !! @tex $(W.m^{-1}.K^{-1})$ @endtex
  REAL(r_std), SAVE   :: soilc_max                      !! g/m^3 from lawrence and slater
!$OMP THREADPRIVATE(soilc_max)

  REAL(r_std), SAVE :: qsintcst                         !! Transforms leaf area index into size of interception reservoir
                                                        !! (unitless)
!$OMP THREADPRIVATE(qsintcst)

  !! Parameters for muff resolution

  REAL(r_std), SAVE :: mcr_sup_param                    !! Residual moisture content in the superficial soil layer
                                                        !! (m^3/m^3)
!$OMP THREADPRIVATE(mcr_sup_param)
  REAL(r_std), SAVE :: mcr_inf_param                    !! Residual moisture content in the inferior soil layer
                                                        !! (m^3/m^3)
!$OMP THREADPRIVATE(mcr_inf_param)
  REAL(r_std), SAVE :: mcs_sup_param                    !! Saturated moisture content in the superficial soil layer
                                                        !! (m^3/m^3)
!$OMP THREADPRIVATE(mcs_sup_param)
  REAL(r_std), SAVE :: mcs_inf_param                    !! Saturated moisture content in the inferior soil layer
                                                        !! (m^3/m^3)
!$OMP THREADPRIVATE(mcs_inf_param)

  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: b_muff_param = &               !! Campbell coefficient b (unitless)
 & (/ 2.708_r_std, 2.708_r_std, 2.708_r_std, 2.708_r_std, 2.708_r_std, 2.708_r_std, &
 &    2.708_r_std, 2.708_r_std, 2.708_r_std, 2.708_r_std, 2.708_r_std, 2.708_r_std, &
 &    2.708_r_std/) 

  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: psi_air_entry_param = &        !! Air entry water potential (unitless)
 & (/ -0.015_r_std, -0.015_r_std, -0.015_r_std,-0.015_r_std, -0.015_r_std, -0.015_r_std, &
 &    -0.015_r_std, -0.015_r_std, -0.015_r_std,-0.015_r_std, -0.015_r_std, -0.015_r_std, &
 &    -0.015_r_std /) 

  !
  !! Parameters for vertical soil discretization
  !
  REAL(r_std), PARAMETER :: minaltmax = 0.1             !! Minimum active layer thickness (m)
  REAL(r_std), PARAMETER :: maxaltmax = 2.              !! Maximum active layer thickness (m)


  !
  !! Parameters specific for the CWRR hydrology.
  !
  !!  1. Parameters for FAO-Zobler Map

  INTEGER(i_std), PARAMETER,DIMENSION(nscm_fao) :: fao2usda = (/ 3,6,9 /) !! To find the values of Coarse, Medium, Fine in Zobler map
                                                                          !! from the USDA lookup tables
  
  !!  2. Parameters for USDA Classification

  !! Parameters for soil type distribution :
  !! Sand, Loamy Sand, Sandy Loam, Silt Loam, Silt, Loam, Sandy Clay Loam, Silty Clay Loam, Clay Loam, Sandy Clay, Silty Clay, Clay
  INTEGER(i_std), SAVE      :: usda_default = 6                          !! Default USDA texture class if no value found from map
!$OMP THREADPRIVATE(usda_default) 

  REAL(r_std), PARAMETER, DIMENSION(nscm_usda) :: soilclass_default = (/0.0, &
       0.0, 0.0, 0.0, 0.0, 1.0, 0.0, &  !! Areal fraction of the 13 soil USDA textures; 
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)   !! the dominant one will selected


  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: nvan_usda = &            !! Van Genuchten coefficient n (unitless)
 & (/ 2.68_r_std, 2.28_r_std, 1.89_r_std, 1.41_r_std, &                   !  RK: 1/n=1-m
 &    1.37_r_std, 1.56_r_std, 1.48_r_std, 1.23_r_std, &
 &    1.31_r_std, 1.23_r_std, 1.09_r_std, 1.09_r_std, & 
 &    1.552_r_std    /) ! oxisols  

  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: avan_usda = &            !! Van Genuchten coefficient a 
 & (/ 0.0145_r_std, 0.0124_r_std, 0.0075_r_std, 0.0020_r_std, &          !!  @tex $(mm^{-1})$ @endtex
 &    0.0016_r_std, 0.0036_r_std, 0.0059_r_std, 0.0010_r_std, &
 &    0.0019_r_std, 0.0027_r_std, 0.0005_r_std, 0.0008_r_std, & 
 &    0.0132_r_std /) ! oxisols 

  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: mcr_usda = &             !! Residual volumetric water content 
 & (/ 0.045_r_std, 0.057_r_std, 0.065_r_std, 0.067_r_std, &              !!  @tex $(m^{3} m^{-3})$ @endtex
 &    0.034_r_std, 0.078_r_std, 0.100_r_std, 0.089_r_std, &
 &    0.095_r_std, 0.100_r_std, 0.070_r_std, 0.068_r_std, & 
 &    0.068_r_std /) ! oxisols 

  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: mcs_usda = &             !! Saturated volumetric water content 
 & (/ 0.43_r_std, 0.41_r_std, 0.41_r_std, 0.45_r_std, &                  !!  @tex $(m^{3} m^{-3})$ @endtex
 &    0.46_r_std, 0.43_r_std, 0.39_r_std, 0.43_r_std, &
 &    0.41_r_std, 0.38_r_std, 0.36_r_std, 0.38_r_std, & 
 &    0.503_r_std  /) ! oxisols 

  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: ks_usda = &              !! Hydraulic conductivity at saturation
 & (/ 7128.0_r_std, 3501.6_r_std, 1060.8_r_std, 108.0_r_std, &           !!  @tex $(mm d^{-1})$ @endtex
 &    60.0_r_std, 249.6_r_std, 314.4_r_std, 16.8_r_std, &
 &    62.4_r_std, 28.8_r_std, 4.8_r_std, 48.0_r_std, & 
 &    6131.4_r_std  /) ! oxisols  

! The max available water content is smaller when mcw and mcf depend on texture,
! so we increase pcent to a classical value of 80%  
  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: pcent_usda = &           !! Fraction of saturated volumetric soil moisture
 & (/ 0.8_r_std, 0.8_r_std, 0.8_r_std, 0.8_r_std, &                      !! above which transpir is max (0-1, unitless)
 &    0.8_r_std, 0.8_r_std, 0.8_r_std, 0.8_r_std, &
 &    0.8_r_std, 0.8_r_std, 0.8_r_std, 0.8_r_std, & 
 &    0.8_r_std /) ! oxisols 

  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: free_drain_max_usda = &  !! Max=default value of the permeability coeff 
 & (/ 1.0_r_std, 1.0_r_std, 1.0_r_std, 1.0_r_std, &                      !! at the bottom of the soil (0-1, unitless)
 &    1.0_r_std, 1.0_r_std, 1.0_r_std, 1.0_r_std, &
 &    1.0_r_std, 1.0_r_std, 1.0_r_std, 1.0_r_std,  & 
 &    1.0_r_std /) 
  
!! We use the VG relationships to derive mcw and mcf depending on soil texture
!! assuming that the matric potential for wilting point and field capacity is
!! -150m (permanent WP) and -3.3m respectively
!! (-1m for FC for the three sandy soils following Richards, L.A. and Weaver, L.R. (1944)
!! Note that mcw GE mcr
  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: mcf_usda = &             !! Volumetric water content at field capacity
 & (/ 0.0493_r_std, 0.0710_r_std, 0.1218_r_std, 0.2402_r_std, &          !!  @tex $(m^{3} m^{-3})$ @endtex
      0.2582_r_std, 0.1654_r_std, 0.1695_r_std, 0.3383_r_std, &
      0.2697_r_std, 0.2672_r_std, 0.3370_r_std, 0.3469_r_std, & 
      0.172_r_std  /) ! oxisols
  
  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: mcw_usda = &             !! Volumetric water content at wilting point
 & (/ 0.0450_r_std, 0.0570_r_std, 0.0657_r_std, 0.1039_r_std, &          !!  @tex $(m^{3} m^{-3})$ @endtex
      0.0901_r_std, 0.0884_r_std, 0.1112_r_std, 0.1967_r_std, &
      0.1496_r_std, 0.1704_r_std, 0.2665_r_std, 0.2707_r_std, & 
      0.075_r_std  /) ! oxisols 

  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: mc_awet_usda = &         !! Vol. wat. cont. above which albedo is cst
 & (/ 0.25_r_std, 0.25_r_std, 0.25_r_std, 0.25_r_std, &                  !!  @tex $(m^{3} m^{-3})$ @endtex
 &    0.25_r_std, 0.25_r_std, 0.25_r_std, 0.25_r_std, &
 &    0.25_r_std, 0.25_r_std, 0.25_r_std, 0.25_r_std, & 
 &    0.25_r_std /) 

  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: mc_adry_usda = &         !! Vol. wat. cont. below which albedo is cst
 & (/ 0.1_r_std, 0.1_r_std, 0.1_r_std, 0.1_r_std, &                      !!  @tex $(m^{3} m^{-3})$ @endtex
 &    0.1_r_std, 0.1_r_std, 0.1_r_std, 0.1_r_std, &
 &    0.1_r_std, 0.1_r_std, 0.1_r_std, 0.1_r_std, & 
 &    0.1_r_std /) ! oxisols

  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: QZ_usda = &              !! QUARTZ CONTENT (SOIL TYPE DEPENDENT)
 & (/ 0.92_r_std, 0.82_r_std, 0.60_r_std, 0.25_r_std, &
 &    0.10_r_std, 0.40_r_std, 0.60_r_std, 0.10_r_std, &
 &    0.35_r_std, 0.52_r_std, 0.10_r_std, 0.25_r_std, &                  !! Peters et al [1998]
 &     0.25_r_std /)  ! oxisols                  

  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: so_capa_dry_usda = &     !! Dry soil Volumetric heat capacity of soils,J.m^{-3}.K^{-1}
 & (/ 1.47e+6_r_std, 1.41e+6_r_std, 1.34e+6_r_std, 1.27e+6_r_std, &
 &    1.21e+6_r_std, 1.21e+6_r_std, 1.18e+6_r_std, 1.32e+6_r_std, &
 &    1.23e+6_r_std, 1.18e+6_r_std, 1.15e+6_r_std, 1.09e+6_r_std, &      !! Pielke [2002, 2013]
 &    1.09e+6_r_std /) ! oxisols  

  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: clayfrac_usda = &        !! % clay particles in the 13 USDA texture classes 
      (/ 0.03_r_std, 0.06_r_std, 0.11_r_std, 0.19_r_std , &              !! values taken from get_soilcorr_usda in slowproc 
         0.10_r_std, 0.20_r_std, 0.27_r_std, 0.33_r_std, & 
         0.33_r_std, 0.41_r_std, 0.46_r_std, 0.55_r_std, & 
         0.55_r_std /) ! oxisols                                                                   
     
  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: sandfrac_usda = &        !! % sand particles in the 13 USDA texture classes  
       (/ 0.93_r_std, 0.81_r_std, 0.63_r_std, 0.17_r_std, &              !! values taken from get_soilcorr_usda in slowproc    
          0.06_r_std, 0.40_r_std, 0.54_r_std, 0.08_r_std, & 
          0.30_r_std, 0.48_r_std, 0.06_r_std, 0.15_r_std, & 
          0.15_r_std /) ! oxisols 


  !! Parameters for the numerical scheme used by CWRR

  INTEGER(i_std), PARAMETER :: imin = 1                                 !! Start for CWRR linearisation (unitless)
  INTEGER(i_std), PARAMETER :: nbint = 50                               !! Number of interval for CWRR linearisation (unitless)
  INTEGER(i_std), PARAMETER :: imax = nbint+1                           !! Number of points for CWRR linearisation (unitless)
  REAL(r_std), PARAMETER    :: w_time = 1.0_r_std                       !! Time weighting for CWRR numerical integration (unitless)


  !! Variables related to soil freezing, in thermosoil : 
  LOGICAL, SAVE        :: ok_Ecorr                    !! Flag for energy conservation correction
!$OMP THREADPRIVATE(ok_Ecorr)
  LOGICAL, SAVE        :: ok_freeze_thermix           !! Flag to activate thermal part of the soil freezing scheme
!$OMP THREADPRIVATE(ok_freeze_thermix)
  LOGICAL, SAVE        :: ok_freeze_thaw_latent_heat  !! Flag to activate latent heat part of the soil freezing scheme
!$OMP THREADPRIVATE(ok_freeze_thaw_latent_heat)
  LOGICAL, SAVE        :: read_reftemp                !! Flag to initialize soil temperature using climatological temperature
!$OMP THREADPRIVATE(read_reftemp)
  LOGICAL, SAVE        :: read_reftempice             !! Flag to initialize ice temperature using equilibrated temperature
!$OMP THREADPRIVATE(read_reftempice)
  REAL(r_std), SAVE    :: fr_dT                       !! Freezing window width (K)
!$OMP THREADPRIVATE(fr_dT)
  REAL(r_std), SAVE    :: fr_center                   !! Freezing window center (K)
!$OMP THREADPRIVATE(fr_center)
  LOGICAL, SAVE        :: use_soilc_insulation        !! Do we want to activate the soil C effect on thermix
!$OMP THREADPRIVATE(use_soilc_insulation)

  !! Variables related to soil freezing, in hydrol : 
  LOGICAL, SAVE        :: ok_freeze_cwrr              !! CWRR freezing scheme by I. Gouttevin
!$OMP THREADPRIVATE(ok_freeze_cwrr)
  LOGICAL, SAVE        :: ok_thermodynamical_freezing !! Calculate frozen fraction thermodynamically
!$OMP THREADPRIVATE(ok_thermodynamical_freezing)
  REAL(r_std), SAVE    :: moss_layer_thickness        !! Thickness of moss layer used in thermosoil (m) 
!$OMP THREADPRIVATE(moss_layer_thickness)
  REAL(r_std), SAVE    :: soil_moss_depth             !! Maximum depth for computation of moss thermal effects (m). 
!$OMP THREADPRIVATE(soil_moss_depth)

  !! Variables related to SOM initialization
  LOGICAL, SAVE        :: use_initsom       !! Initialize SOM using SOM map if no restart file
!$OMP THREADPRIVATE(use_initsom)
  REAL(r_std), SAVE    :: initsom_frac_a    !! Fraction of initSOM to put in active C pool
!$OMP THREADPRIVATE(initsom_frac_a)  
  REAL(r_std), SAVE    :: initsom_frac_s    !! Fraction of initSOM to put in slow C pool
!$OMP THREADPRIVATE(initsom_frac_s)
  REAL(r_std), SAVE    :: initsom_frac_p    !! Fraction of initSOM to put in passive C pool
!$OMP THREADPRIVATE(initsom_frac_p)
  
  ! PEATLAND
  !
  INTEGER(i_std), PARAMETER ::  dim_moyanno_moist = 45 !! The dimension used in Moyano et al., 2012, volumetric moisture, 0.01 interval 
  REAL(r_std), PARAMETER    ::  nvan_peat=1.38    !! Van Genuchten coeficients for peat
  REAL(r_std), PARAMETER    ::  avan_peat=0.00507 !! Van Genuchten coeficients for peat
  REAL(r_std), PARAMETER    ::  mcr_peat=0.15     !! Residual volumetric water content (m^{3} m^{-3}) 
  REAL(r_std), PARAMETER    ::  mcs_peat=0.90     !! Saturated volumetric water content (m^{3} m^{-3})
  REAL(r_std), PARAMETER    ::  ks_peat=2120_r_std ! Hydraulic conductivity at saturation (mm {-1})
  REAL(r_std), PARAMETER    ::  mcw_peat=0.210    !! wilting-point volumetric water content (m^{3} m^{-3})
  REAL(r_std), PARAMETER    ::  mcf_peat=0.406    !! volumetric water content related to field capacity (m^{3} m^{-3})
  REAL(r_std), PARAMETER    ::  mc_awet_peat =0.25 ! Vol. wat. cont. above which albedo is cst
  REAL(r_std), PARAMETER    ::  mc_adry_peat =0.1  ! Vol. wat. cont. below which albedo is cst
  REAL(r_std), PARAMETER    ::  pcent_peat=0.8 ! Fraction of saturated volumetric soil moisture above which transpir is max
  REAL(r_std), PARAMETER, DIMENSION(18) :: peat_bulk_density =  &    !!in g/cm3, median of core measurement
  & (/ 0.063, 0.063, 0.063, 0.047, 0.064, 0.066, 0.079, 0.091, 0.101, 0.109, 0.109, &
  &    0.103, 0.112, 0.110, 0.104, 0.073, 0.092, 0.092 /)
  
  REAL(r_std), SAVE ::  tau_peat     !!k0 = 0.1yr**-1, tau_peat in s, turnover rate for peat
!$OMP THREADPRIVATE(tau_peat)
  REAL(r_std), SAVE ::  z_tau        !!the e-folding depth of turnover rates
!$OMP THREADPRIVATE(z_tau)
  

END MODULE constantes_soil_var
