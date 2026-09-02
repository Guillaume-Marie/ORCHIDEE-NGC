!=================================================================================================================================
! MODULE       : constantes_mtc
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2011)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF         This module contains the standard values of the parameters for the 13 metaclasses of vegetation used by ORCHIDEE.
!!
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): 
!!
!! REFERENCE(S)	: 
!! - Kuppel, S. (2012): Doctoral Thesis, Assimilation de mesures de flux turbulents d'eau et de carbone dans un modèle de la biosphère 
!! continentale 
!! - Kuppel, S., Peylin, P., Chevallier, F., Bacour, C., Maignan, F., and Richardson, A. D. (2012). Constraining a global ecosystem
!! model with multi-site eddy-covariance data, Biogeosciences, 9, 3757-3776, DOI 10.5194/bg-9-3757-2012.
!! - Wohlfahrt, G., M. Bahn, E. Haubner, I. Horak, W. Michaeler, K.Rottmar, U. Tappeiner, and A. Cemusca, 1999: Inter-specific 
!! variation of the biochemical limitation to photosynthesis and related leaf traits of 30 species from mountain grassland 
!! ecosystems under different land use. Plant Cell Environ., 22, 12811296.
!! - Malhi, Y., Doughty, C., and Galbraith, D. (2011). The allocation of ecosystem net primary productivity in tropical forests, 
!! Philosophical Transactions of the Royal Society B-Biological Sciences, 366, 3225-3245, DOI 10.1098/rstb.2011.0062.
!! - Earles, J. M., Yeh, S., and Skog, K. E. (2012). Timing of carbon emissions from global forest clearance, Nature Climate Change, 2, 
!! 682-685, Doi 10.1038/Nclimate1535.
!! - Piao, S. L., Luyssaert, S., Ciais, P., Janssens, I. A., Chen, A. P., Cao, C., Fang, J. Y., Friedlingstein, P., Luo, Y. Q., and 
!! Wang, S. P. (2010). Forest annual carbon cost: A global-scale analysis of autotrophic respiration, Ecology, 91, 652-661, 
!! Doi 10.1890/08-2176.1.
!! - Verbeeck, H., Peylin, P., Bacour, C., Bonal, D., Steppe, K., and Ciais, P. (2011). Seasonal patterns of co2 fluxes in amazon 
!! forests: Fusion of eddy covariance data and the orchidee model, Journal of Geophysical Research-Biogeosciences, 116, 
!! Artn G02018, Doi 10.1029/2010jg001544.
!! - MacBean, N., Maignan, F., Peylin, P., Bacour, C., Breon, F. M., & Ciais, P. (2015). Using satellite data to improve the leaf 
!! phenology of a global terrestrial biosphere model. Biogeosciences, 12(23), 7185-7208.
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_parameters/constantes_mtc.f90 $
!! $Date: 2026-05-13 15:11:27 +0200 (mer. 13 mai 2026) $
!! $Revision: 9540 $
!! \n
!_ ================================================================================================================================

MODULE constantes_mtc

  USE defprec
  USE constantes_var

  IMPLICIT NONE

  !
  ! METACLASSES CHARACTERISTICS
  !

  INTEGER(i_std), PARAMETER :: nvmc = 13                         !! Number of MTCS fixed in the code (unitless)

  CHARACTER(len=35), PARAMETER, DIMENSION(nvmc) :: MTC_name = &  !! description of the MTC (unitless)
  & (/ 'Bare ground                        ', &          !  1
  &    'Tropical Evergreen Broadleaf tree  ', &          !  2
  &    'Tropical Deciduous Broadleaf tree  ', &          !  3
  &    'Temperate Evergreen Needleleaf tree', &          !  4
  &    'Temperate Evergreen Broadleaf tree ', &          !  5
  &    'Temperate Deciduous Broadleaf tree ', &          !  6
  &    'Boreal Evergreen Needleleaf tree   ', &          !  7
  &    'Boreal Deciduous Broadleaf tree    ', &          !  8
  &    'Boreal Deciduous Needleleaf tree   ', &          !  9
  &    'C3 grass                           ', &          ! 10
  &    'C4 grass                           ', &          ! 11
  &    'C3 crop                            ', &          ! 12
  &    'C4 crop                            '  /)         ! 13
  

  !
  ! VEGETATION STRUCTURE
  !
  INTEGER(i_std),PARAMETER, DIMENSION(nvmc) :: leaf_tab_mtc  =  &                 !! leaf type (1-4, unitless)
  & (/  4,   1,   1,   2,   1,   1,   2,   &                                      !! 1=broad leaved tree, 2=needle leaved tree
  &     1,   2,   3,   3,   3,   3   /)                                           !! 3=grass 4=bare ground

  CHARACTER(len=6), PARAMETER, DIMENSION(nvmc) :: pheno_model_mtc  =  &  !! which phenology model is used? (tabulated) 
  & (/  'none  ',   'none  ',   'moi   ',   'none  ',   'none  ',  &
  &     'ncdgdd',   'none  ',   'ncdgdd',   'ngd   ',   'moigdd',  &
  &     'moigdd',   'moigdd',   'moigdd'  /) 

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: gdd_threshold_mtc  =  &    !! threshold temperature above which gdd is counted 
  & (/ 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  &                !! (degrees C)
  &    0.0,   0.0,   0.0,   5.0,   5.0,   5.0  /)  

  LOGICAL, PARAMETER, DIMENSION(nvmc) :: is_tropical_mtc  =  &                       !! Is PFT tropical ? (true/false)
  & (/ .FALSE.,   .TRUE.,    .TRUE.,    .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,  &
  &    .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE. /)    

  LOGICAL, PARAMETER, DIMENSION(nvmc) :: is_temperate_mtc  =  &         !! Is PFT temperate ? (true/false)
  & (/ .FALSE.,   .FALSE.,   .FALSE.,   .TRUE.,    .TRUE.,    .TRUE.,   .FALSE.,  &
  &    .FALSE.,   .FALSE.,   .TRUE.,    .TRUE.,    .TRUE.,    .TRUE. /)

  LOGICAL, PARAMETER, DIMENSION(nvmc) :: is_boreal_mtc = &              !! Is PFT boreal ? (true/false)
  & (/ .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .TRUE.,  &
  &    .TRUE.,    .TRUE.,    .FALSE.,   .FALSE.,   .FALSE.,   .FALSE. /)

  CHARACTER(LEN=5), PARAMETER, DIMENSION(nvmc) :: type_of_lai_mtc  =  &  !! Type of behaviour of the LAI evolution algorithm
  & (/ 'inter', 'inter', 'inter', 'inter', 'inter',  &                   !! for each vegetation type. (unitless)
  &    'inter', 'inter', 'inter', 'inter', 'inter',  &                   !! Value of type_of_lai : mean or interp
  &    'inter', 'inter', 'inter' /)

  LOGICAL, PARAMETER, DIMENSION(nvmc) :: natural_mtc =  &                         !! natural?  (true/false)
  & (/ .TRUE.,   .TRUE.,   .TRUE.,   .TRUE.,   .TRUE.,    .TRUE.,   .TRUE.,  &
  &    .TRUE.,   .TRUE.,   .TRUE.,   .TRUE.,   .FALSE.,   .FALSE.  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: llaimax_mtc  =  &          !! laimax for maximum
  & (/ 0.0,   8.0,   8.0,   4.0,   4.5,   4.5,   4.0,  &                !! See also type of lai interpolation
  &    4.5,   4.0,   2.0,   2.0,   2.0,   2.0  /)                       !! @tex $(m^2.m^{-2})$ @endtex

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: llaimin_mtc  = &           !! laimin for minimum lai
  & (/ 0.0,   8.0,   0.0,   4.0,   4.5,   0.0,   4.0,  &                !! See also type of lai interpolation (m^2.m^{-2})
  &    0.0,   0.0,   0.0,   0.0,   0.0,   0.0  /)                       !! @tex $(m^2.m^{-2})$ @endtex

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: height_presc_mtc  =  &     !! prescribed height of vegetation (m)
  & (/  0.0,   30.0,   30.0,   20.0,   20.0,   20.0,   15.0,  &         !! Value for height_presc : one for each vegetation type
  &    15.0,   15.0,    0.5,    0.6,    1.0,    1.0  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: z0_over_height_mtc = &         !! Factor to calculate roughness height from 
  & (/  0.0, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625,  &         !! vegetation height (unitless)   
  &  0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: ratio_z0m_z0h_mtc = &      !! Ratio between z0m and z0h values (roughness height for momentum and for heat)
  & (/ 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0,  &         
  &    10.0, 10.0, 10.0, 10.0, 10.0, 10.0  /)


  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: rveg_mtc  =  &             !! Potentiometer to set vegetation resistance (unitless)
  & (/ 1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &                !! Nathalie on March 28th, 2006 - from Fred Hourdin,
  &    1.0,   1.0,   1.0,   1.0,   1.0,   1.0   /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: sla_mtc  =  &              !! specif leaf area @tex $(m^2.gC^{-1})$ @endtex
  & (/ 1.5E-2,   1.53E-2,   2.6E-2,   9.26E-3,     2E-2,   2.6E-2,   9.26E-3,  &
  &    2.6E-2,    1.9E-2,   2.6E-2,    2.6E-2,   2.6E-2,   2.6E-2  /) 

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: slainit_mtc  =  &          !! specif leaf area @tex $(m^2.gC^{-1})$ @endtex
  & (/ 0.026,     0.03963,   0.044,   0.00575,  0.03604,  0.06131, 0.008153, &
  &    0.06822,   0.03684,   0.031,   0.031,    0.02,     0.02  /) 

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: availability_fact_mtc  = & !! calculate mortality in lpj_gap
  & (/ undef,   0.14,  0.14,   0.10,   0.10,   0.10,   0.05,  &
  &     0.05,   0.05,  undef,  undef,  undef,  undef  /)

  !
  ! EVAPOTRANSPIRATION (sechiba)
  !
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: rstruct_const_mtc  =  &  !! Structural resistance.
  & (/ 0.0,   25.0,   25.0,   25.0,   25.0,   25.0,   25.0,  &        !! @tex $(s.m^{-1})$ @endtex
  &   25.0,   25.0,    2.5,    2.0,    2.0,    2.0   /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: kzero_mtc  =  &                  !! A vegetation dependent constant used in the 
  & (/    0.0,   12.E-5,   12.E-5,   12.E-5,   12.E-5,   25.E-5,   12.E-5,  & !! calculation  of the surface resistance. 
  &    25.E-5,   25.E-5,   30.E-5,   30.E-5,   30.E-5,   30.E-5  /)           !! @tex $(kg.m^2.s^{-1})$ @endtex

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: alpha_watstress_pft_mtc = &          !!parameter of water stress to control the exp function in each PFT
  & (/    1.0,  1.0,      1.0,      1.0,      1.0,       1.0,      1.0, &
          1.0,  1.0,      1.0,      1.0,      1.0,       1.0     /) 

  !
  ! WATER (sechiba)
  !
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: wmax_veg_mtc  =  &        !! Volumetric available soil water capacity in each PFT
  & (/ 150.0,   150.0,   150.0,   150.0,   150.0,   150.0,   150.0,  & !! @tex $(kg.m^{-3} of soil)$ @endtex
  &    150.0,   150.0,   150.0,   150.0,   150.0,   150.0  /)         
                                                                      

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: humcste_ref4m  =  &       !! Root profile description for the different 
  & (/ 5.0,   0.4,   0.4,   1.0,   0.8,   0.8,   1.0,  &               !! vegetations types. @tex $(m^{-1})$ @endtex
  &    1.0,   0.8,   4.0,   1.0,   4.0,   1.0  /)                      !! These are the factor in the exponential which gets       
                                                                       !! the root density as a function of depth
                                                                       !! Values for zmaxh = 4.0  
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: max_root_depth_ref4m_mtc=&!! Maximum rooting depth for the PFT irrespective of other 
  & (/ 4.0,   4.0,   4.0,   4.0,   4.0,   4.0,   4.0,  &               !! constraints from the active layer thickness @tex $(m)$ @endtex
  &    4.0,   4.0,   2.0,   2.0,   0.8,   0.8  /)                      !! Values for zmaxh = 4.0 
  
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: humcste_ref2m  =  &       !! Root profile description for the different
  & (/ 1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &               !! vegetations types.  @tex $(m^{-1})$ @endtex
  &    1.0,   1.0,   0.6,   0.6,   0.6,   0.6  /)                      !! These are the factor in the exponential which gets       
                                                                       !! the root density as a function of depth
                                                                       !! Values for zmaxh = 2.0

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: max_root_depth_ref2m_mtc=&!! Maximum rooting depth for the PFT irrespective of other 
  & (/ 2.0,   2.0,   2.0,   2.0,   2.0,   2.0,   2.0,  &               !! constraints from the active layer thickness @tex $(m)$ @endtex
  &    2.0,   2.0,   1.0,   1.0,   0.8,   0.8  /)      

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: throughfall_by_mtc  =  &  !! Percent by PFT of precip that is not intercepted by the canopy
  & (/ 30.0,   35.0,   35.0,   30.0,   30.0,   30.0,   30.0,  &        !! (0-100, unitless)
  &    30.0,   30.0,   30.0,   30.0,   30.0,   30.0  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: rsoil_scale_mtc  =  &     !! Scaling factor to the soil resistance to evaporate per mtc
  & (/ 1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &               
  &    1.0,   1.0,   1.0,   1.0,   1.0,   1.0  /)
 

  !
  ! ALBEDO (sechiba)
  !
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: snowa_aged_vis_mtc  =  &  !! Minimum snow albedo value for each vegetation type
  & (/ 0.68,    0.24,    0.07,  0.5,    0.5,    0.125,  0.377,  &      !! after aging (dirty old snow) (unitless), visible albedo
  &    0.5,     0.303,   0.68,  0.68,   0.68,   0.68  /)               !! Source : Values optimized for ORCHIDEE2.0

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: snowa_aged_nir_mtc  =  &  !! Minimum snow albedo value for each vegetation type
  & (/ 0.50,    0.37,   0.08,   0.10,   0.37,   0.08,   0.16,  &       !! after aging (dirty old snow) (unitless), near infrared albedo
  &    0.17,    0.27,   0.44,   0.44,   0.44,   0.44  /)               !! Source : Values optimized for ORCHIDEE2.0

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: snowa_dec_vis_mtc  =  &   !! Decay rate of snow albedo value for each vegetation type
  & (/ 0.3,    0.14,    0.08,   0.01,   0.01,   0.01,   0.01,  &       !! as it will be used in condveg_snow (unitless), visible albedo
  &    0.3,    0.01,    0.3,    0.3,    0.3,    0.3  /)                !! Source : Values optimized for ORCHIDEE2.0

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: snowa_dec_nir_mtc  =  &   !! Decay rate of snow albedo value for each vegetation type
  & (/ 0.13,    0.10,   0.16,   0.10,   0.10,   0.16,   0.04,  &       !! as it will be used in condveg_snow (unitless), near infrared albedo
  &    0.07,    0.08,   0.12,   0.12,   0.12,   0.12 /)                !! Source : Values optimized for ORCHIDEE2.0

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: alb_leaf_vis_mtc  =  &    !! leaf albedo of vegetation type, visible albedo, optimized on 04/07/2016
  & (/ 0.00,   0.04, 0.04, 0.04, 0.04, 0.03, 0.03,  &                  !! (unitless)
  &    0.03,   0.03, 0.06, 0.06, 0.06, 0.06  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: alb_leaf_nir_mtc  =  &    !! leaf albedo of vegetation type, near infrared albedo, optimized on 04/07/2016
  & (/ 0.00,   0.23,  0.18,  0.18,  0.20,  0.24,  0.15,  &             !! (unitless)
  &    0.26,   0.20,  0.24,  0.27,  0.28,  0.26  /)

  ! albedo values for albedo type 'pinty'
  ! these next values were determined by fitting to global MODIS data and using the inversion scheme of
  ! Pinty et al (see Pinty B,Andredakis I, Clerici M, et al. (2011) ! 'Exploiting the MODIS albedos 
  ! with the Two-stream Inversion Package (JRC-TIP): 1. Effective leaf area index,
  ! vegetation, and soil properties'. Journal of Geophysical Research.
  !
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: leaf_ssa_vis_mtc  =  &           !! leaf single scattering albedo, visible light (unitless)
  (/ 0.173, 0.0914, 0.0543, 0.0658, 0.0831,  0.05, 0.0997, &
     0.0956, 0.101, 0.175,  0.242,  0.222,  0.221 /)
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: leaf_ssa_nir_mtc  =  &           !! leaf single scattering albedo, near infrared (unitless)
  (/ 0.703, 0.606, 0.680,  0.617,  0.611,  0.702, 0.600, &
     0.741, 0.607,  0.697, 0.662,  0.699, 0.790/)
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: leaf_psd_vis_mtc  =  &           !! leaf preferred scattering direction, visible light (unitless)
  (/ 1.00, 0.855,  0.5, 0.5, 0.5,  0.5, 0.517, &
     0.51, 0.625, 1.71, 1.67, 2.0,  1.68 /)
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: leaf_psd_nir_mtc  =  &           !! leaf preferred scattering direction, NIR light (unitless)
  (/ 2.01,  1.82,  1.7,  1.47,  1.09,  1.67,  0.5, &
     1.93,  0.744, 1.7,  1.54,  2.22,  2.06 /)
    
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: bgd_reflectance_vis_mtc  =  &    !! background reflectance, visible light (unitless)
  (/ 0.13128,  0.08002,  0.10057,  0.05125,  0.0559,  0.08873,  0.03511, &
     0.0597,   0.05174,  0.13085,  0.1177,   0.11008, 0.124 /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: bgd_reflectance_nir_mtc  =  &    !! background reflectance, NIR light (unitless)
  (/ 0.24556,  0.13978,  0.17651,  0.08851,  0.0957,  0.14427,  0.04983, &
     0.09275,  0.08372,  0.24711,  0.22773,  0.20472, 0.23701 /)


  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: leaf_to_shoot_clumping_mtc  =  & !! The clumping factor for leaves to shoots in the 
  & (/ un,   un,   un,   un,   un,   un,   un,  &                    !! effective LAI calculation...notice this should be
  &   un,   un,   un,   un,   un,   un /)                            !! equal to unity for grasslands/croplands

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: lai_correction_factor_mtc  =  &  !! see the note about this variable in pft_parameters
  & (/ un,   un,   un,   un,   un,   un,   un,  &  
  &   un,   un,   un,   un,   un,   un /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: min_level_sep_mtc  =  & !! The minimum level thickness for photosynthesis [m]
  & (/ un,   0.1,   0.1,   0.1,   0.1,   0.1,   0.1,  &              !! This number is arbitrary at the moment.  The idea
  &   0.1,   0.1,   0.1,   0.1,   0.1,   0.1 /)                      !! is to have a small number to make as many levels
                                                                     !! as possible in the canopies, but not too small which
                                                                     !! results in too little LAI in all the levels.  If all your
                                                                     !! levels have less than 0.1 LAI in them, that's probably 
                                                                     !! too small.

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: lai_top_mtc = &                  !! Diffuco.f90 calculates the stomatal conductance of the
  (/ 1.0,   0.1,   0.1,   0.2,   0.1,   0.1,   0.2,  &                         !! top layer of the canopy. Because the top layer can contain
     0.1,   0.1,   5.0,   5.0,   5.0,   5.0 /)                                !! different amounts of LAI depending on the crown diameter
                                                                              !! we had to define top layer in terms of the LAI it contains.
                                                                              !! stomatal conductance in the top layer contributes to the 
                                                                              !! transpiration (m2 m-2). Arbitrary values.

  !
  ! SOIL - VEGETATION
  !
  INTEGER(i_std), PARAMETER, DIMENSION(nvmc) :: pref_soil_veg_mtc  =  &       !! The soil tile number for each vegetation
  & (/ 1,   2,   2,   2,   2,   2,   2,  &                                    
  &    2,   2,   3,   3,   3,   3  /)                                         

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: moss_frac_mtc  =  &            !! Moss fraction
  & (/ 0.,   0.,   0.,   0.,   0.,   0.,   1.,  &
  &    1.,   1.,   0.,   0.,   0.,   0.  /)

  !
  ! VEGETATION - AGE CLASSES
  !
  INTEGER(i_std), PARAMETER, DIMENSION(nvmc) :: agec_group_mtc  =  &       !! The age class group that each PFT belongs to.
       (/ 1,   2,   3,   4,   5,   6,   7,  &                                    
       8,   9,   10,   11,   12,   13  /)

  !
  ! PHOTOSYNTHESIS
  !
  !-
  ! 1 .CO2
  !-
  LOGICAL, PARAMETER, DIMENSION(nvmc) :: is_c4_mtc  =  &                            !! flag for C4 vegetation types (true/false)
  & (/ .FALSE.,  .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,  &
  &    .FALSE.,  .FALSE.,   .FALSE.,   .TRUE.,    .FALSE.,   .TRUE.  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: vcmax_fix_mtc  =  &     !! values used for vcmax when STOMATE is not
  & (/  0.0,   40.0,   50.0,   30.0,   35.0,   40.0,   30.0,  &      !! activated @tex $(\mu mol.m^{-2}.s^{-1})$ @endtex
  &    40.0,   35.0,   60.0,   60.0,   70.0,   70.0  /)

! For C4 plant we define a very small downregulation effect as C4 plant are
! currently saturate with respect to CO2 impact on vcmax
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: downregulation_co2_coeff_mtc  =  &  !! coefficient for CO2 downregulation 
  & (/  0.0,   0.38,   0.38,   0.28,   0.28,   0.28,   0.22,  &
  &     0.22,  0.22,   0.26,   0.03,   0.26,   0.03 /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: E_KmC_mtc  = &            !! Energy of activation for KmC (J mol-1)
  & (/undef,  79430.,  79430.,  79430.,  79430.,  79430.,  79430.,  &  !! See Medlyn et al. (2002) 
  &  79430.,  79430.,  79430.,  79430.,  79430.,  79430.  /)           !! from Bernacchi al. (2001)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: E_KmO_mtc  = &            !! Energy of activation for KmO (J mol-1)
  & (/undef,  36380.,  36380.,  36380.,  36380.,  36380.,  36380.,  &  !! See Medlyn et al. (2002) 
  &  36380.,  36380.,  36380.,  36380.,  36380.,  36380.  /)           !! from Bernacchi al. (2001)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: E_Sco_mtc  = &            !! Energy of activation for Sco (J mol-1)
  & (/undef, -24460., -24460., -24460., -24460., -24460., -24460.,  &  !! See Table 2 of Yin et al. (2009)
  & -24460., -24460., -24460., -24460., -24460., -24460.  /)           !! Value for C4 plants is not mentioned - We use C3 for all plants

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: E_gamma_star_mtc  = &     !! Energy of activation for gamma_star (J mol-1)
  & (/undef,  37830.,  37830.,  37830.,  37830.,  37830.,  37830.,  &  !! See Medlyn et al. (2002) from Bernacchi al. (2001) 
  &  37830.,  37830.,  37830.,  37830.,  37830.,  37830.  /)           !! for C3 plants - We use the same values for C4 plants

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: E_Vcmax_mtc  = &          !! Energy of activation for Vcmax (J mol-1)
  & (/undef,  71513.,  71513.,  71513.,  71513.,  71513.,  71513.,  &  !! See Table 2 of Yin et al. (2009) for C4 plants
  &  71513.,  71513.,  71513.,  67300.,  71513.,  67300.  /)           !! and Kattge & Knorr (2007) for C3 plants (table 3)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: E_Jmax_mtc  = &            !! Energy of activation for Jmax (J mol-1)
  & (/undef,  49884.,  49884.,  49884.,  49884.,  49884.,  49884.,  &   !! See Table 2 of Yin et al. (2009) for C4 plants
  &  49884.,  49884.,  49884.,  77900.,  49884.,  77900.  /)            !! and Kattge & Knorr (2007) for C3 plants (table 3)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: aSV_mtc     = &            !! a coefficient of the linear regression (a+bT) defining the Entropy term for Vcmax (J K-1 mol-1)
  & (/undef,  668.39,  668.39,  668.39,  668.39,  668.39,  668.39,  &   !! See Table 3 of Kattge & Knorr (2007)
  &  668.39,  668.39,  668.39,  641.64,  668.39,  641.64  /)            !! For C4 plants, we assume that there is no
                                                                        !! acclimation and that at for a temperature of 25°C, aSV is the same for both C4 and C3 plants (no strong jusitification - need further parametrization)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: bSV_mtc     = &            !! b coefficient of the linear regression (a+bT) defining the Entropy term for Vcmax (J K-1 mol-1 °C-1)
  & (/undef,   -1.07,   -1.07,   -1.07,   -1.07,   -1.07,   -1.07,  &   !! See Table 3 of Kattge & Knorr (2007)
  &   -1.07,   -1.07,   -1.07,      0.,   -1.07,      0.  /)            !! We assume No acclimation term for C4 plants

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: tphoto_min_mtc  =  &       !! minimum photosynthesis temperature (deg C) 
  & (/  undef,   -4.0,    -4.0,   -4.0,   -4.0,   -4.0,   -4.0,  & 
  &      -4.0,   -20.0,   -4.0,   -4.0,   -4.0,   -4.0  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: tphoto_max_mtc  =  &       !! maximum photosynthesis temperature (deg C) 
  & (/  undef,   55.0,    55.0,   55.0,   55.0,   55.0,   55.0,  & 
  &      55.0,   55.0,    55.0,   55.0,   55.0,   55.0  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: aSJ_mtc     = &            !! a coefficient of the linear regression (a+bT) defining the Entropy term for Jmax (J K-1 mol-1)
  & (/undef,  659.70,  659.70,  659.70,  659.70,  659.70,  659.70,  &   !! See Table 3 of Kattge & Knorr (2007)
  &  659.70,  659.70,  659.70,    630.,  659.70,    630.  /)            !! and Table 2 of Yin et al. (2009) for C4 plants

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: bSJ_mtc     = &            !! b coefficient of the linear regression (a+bT) defining the Entropy term for Jmax (J K-1 mol-1 °C-1)
  & (/undef,   -0.75,   -0.75,   -0.75,   -0.75,   -0.75,   -0.75,  &   !! See Table 3 of Kattge & Knorr (2007)
  &   -0.75,   -0.75,   -0.75,      0.,   -0.75,      0.  /)            !! We assume no acclimation term for C4 plants

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: D_Vcmax_mtc  = &           !! Energy of deactivation for Vcmax (J mol-1)
  & (/undef, 200000., 200000., 200000., 200000., 200000., 200000.,  &   !! Medlyn et al. (2002) also uses 200000. for C3 plants (same value than D_Jmax)
  & 200000., 200000., 200000., 192000., 200000., 192000.  /)            !! 'Consequently', we use the value of D_Jmax for C4 plants

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: D_Jmax_mtc  = &            !! Energy of deactivation for Jmax (J mol-1)
  & (/undef, 200000., 200000., 200000., 200000., 200000., 200000.,  &   !! See Table 2 of Yin et al. (2009)
  & 200000., 200000., 200000., 192000., 200000., 192000.  /)            !! Medlyn et al. (2002) also uses 200000. for C3 plants

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: E_gm_mtc  = &              !! Energy of activation for gm (J mol-1) 
  & (/undef,  49600.,  49600.,  49600.,  49600.,  49600.,  49600.,  &   !! See Table 2 of Yin et al. (2009) 
  &  49600.,  49600.,  49600.,   undef,  49600.,   undef  /)            
	 	 
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: S_gm_mtc  = &              !! Entropy term for gm (J K-1 mol-1) 
  & (/undef,   1400.,   1400.,   1400.,   1400.,   1400.,   1400.,  &   !! See Table 2 of Yin et al. (2009) 
  &   1400.,   1400.,   1400.,   undef,   1400.,   undef  /) 
	 	 
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: D_gm_mtc  = &              !! Energy of deactivation for gm (J mol-1) 
  & (/undef, 437400., 437400., 437400., 437400., 437400., 437400.,  &   !! See Table 2 of Yin et al. (2009) 
  & 437400., 437400., 437400.,   undef, 437400.,   undef  /)            

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: E_Rd_mtc  = &              !! Energy of activation for Rd (J mol-1)
  & (/undef,  46390.,  46390.,  46390.,  46390.,  46390.,  46390.,  &   !! See Table 2 of Yin et al. (2009)
  &  46390.,  46390.,  46390.,  46390.,  46390.,  46390.  /)           

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: Vcmax25_mtc  =  &          !! Maximum rate of Rubisco activity-limited carboxylation at 25°C
  & (/ undef,   45.0,    45.0,    35.0,   40.0,   50.0,   45.0,  &      !! @tex $(\mu mol.m^{-2}.s^{-1})$ @endtex
  &     35.0,   35.0,    50.0,    50.0,   60.0,   60.0  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: arJV_mtc    = &            !! a coefficient of the linear regression (a+bT) defining the Jmax25/Vcmax25 ratio (mu mol e- (mu mol CO2)-1)
  & (/undef,    2.1,    2.5,    2.04,    2.22,    2.18,    2.28,  &   !! See Table 3 of Kattge & Knorr (2007)
  &    2.18,    2.0,    2.0,    2.0,     2.0,     2.0  /)            !! For C4 plants, we assume that there is no
                                                                        !! acclimation and that for a temperature of 25°C, aSV is the same for both C4 and C3 plants (no strong jusitification - need further parametrization)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: brJV_mtc    = &            !! b coefficient of the linear regression (a+bT) defining the Jmax25/Vcmax25 ratio ((mu mol e- (mu mol CO2)-1) (°C)-1)
  & (/undef,  -0.035,  -0.035,  -0.035,  -0.035,  -0.035,  -0.035,  &   !! See Table 3 of Kattge & Knorr (2007)
  &  -0.035,  -0.035,  -0.035,      0.,  -0.035,      0.  /)            !! We assume No acclimation term for C4 plants

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: KmC25_mtc  = &             !! Michaelis–Menten constant of Rubisco for CO2 at 25°C (ubar)
  & (/undef,   404.9,   404.9,   404.9,   404.9,  404.9,   404.9,  &    !! See Table 2 of Yin et al. (2009) for C4
  &   404.9,   404.9,   404.9,    650.,   404.9,   650.  /)             !! and Medlyn et al (2002) for C3

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: KmO25_mtc  = &             !! Michaelis–Menten constant of Rubisco for O2 at 25°C (ubar)
  & (/undef, 278400., 278400., 278400., 278400., 278400., 278400.,  &   !! See Table 2 of Yin et al. (2009) for C4 plants and Medlyn et al. (2002) for C3
  & 278400., 278400., 278400., 450000., 278400., 450000.  /)           

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: Sco25_mtc  = &             !! Relative CO2 /O2 specificity factor for Rubisco at 25Â°C (bar bar-1)
  & (/undef,   2800.,   2800.,   2800.,   2800.,   2800.,   2800.,  &   !! See Table 2 of Yin et al. (2009)
  &   2800.,   2800.,   2800.,   2590.,   2800.,   2590.  /)           

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: gm25_mtc  = &              !! Mesophyll diffusion conductance at 25Â°C (mol m-2 s-1 bar-1) 
  & (/undef,     0.4,     0.4,     0.4,     0.4,    0.4,      0.4,  &   !! See legend of Figure 6 of Yin et al. (2009) 
  &     0.4,     0.4,     0.4,   undef,     0.4,  undef  /)             !! and review by Flexas et al. (2008) - gm is not used for C4 plants 

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: gamma_star25_mtc  = &      !! Ci-based CO2 compensation point in the absence of Rd at 25°C (ubar)
  & (/undef,   42.75,   42.75,   42.75,   42.75,   42.75,   42.75,  &   !! See Medlyn et al. (2002) for C3 plants - For C4 plants, we use the same value (probably uncorrect)
  &   42.75,   42.75,   42.75,   42.75,   42.75,   42.75  /)    

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: a1_mtc  = &                !! Empirical factor involved in the calculation of fvpd (-)
  & (/undef,    0.94,    0.92,    0.90,    0.93,    0.95,  0.90,  &     !! Adjusted from Table 2 of Yin et al. (2009)
  &    0.95,    0.90,    0.95,    0.85,    0.96,    0.95  /)           

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: b1_mtc  = &                !! Empirical factor involved in the calculation of fvpd (-)
  & (/undef,    0.048,   0.053,   0.056,   0.051,   0.044, 0.056,  &    !! Adjusted from Table 2 of Yin et al. (2009)
  &    0.044,   0.056,   0.045,   0.058,   0.039,   0.046  /)           

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: a_fpsi_mtc = &             !!Empirical factor involved in the calculation of fpsi (-)
  & (/undef,    9.1,   3.67,    3.14,   6.81,    2.27,     3.02,  &     !! See Table 2 of Tuzet et al. (2003) Optimized by JAlleon
  &    6.85,    2.,    6.96,    3.09,   6.74,    7.04  /)

   REAL(r_std), PARAMETER, DIMENSION(nvmc) :: sf_mtc  = &               !! Sensitivity parameter involved in the calculation of f_psi
  & (/undef,    2.5,    2.82,    2.89,  2.74,    3.76,     2.72,   &    !! See Table 1 of Tuzet et al. (2017) Optimized by JAlleon
  &    3.18,    4.3,    3.08,    2.64,  3.63,    3.32  /)
  
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: Rroot_sup_mtc  = &         !! Superficial Roots Water resistance involved in the calculation of psi_leaf
  & (/undef, 500000., 500000., 500000., 500000., 500000., 500000.,  &   !! See Table 1 of Tuzet et al. (2017)
  &  500000., 500000., 500000., 500000., 500000., 500000. /)           

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: Rroot_inf_mtc  = &         !! Inferior Roots Water resistance involved in the calculation of psi_leaf
  & (/undef, 2000000., 2000000., 2000000., 2000000., 2000000., 2000000., &    !! See Table 1 of Tuzet et al. (2017)
  &  2000000., 2000000., 2000000., 2000000., 2000000., 2000000. /)          

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: Rleaf_mtc  = &             !! Leaf Water resistance involved in the calculation of psi_leaf
  & (/undef, 500000., 500000., 500000., 500000., 500000., 500000.,  &   !! See Table 1 of Tuzet et al. (2017)
  &  500000., 500000., 500000., 500000., 500000., 500000. /)          

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: Rxylem_mtc  = &            !! Xylem Water resistance involved in the calculation of psi_leaf
  & (/undef, 2500000., 2500000., 2500000., 2500000., 2500000., 2500000., &    !! See Table 1 of Tuzet et al. (2017)
  &  2500000., 2500000., 2500000., 2500000., 2500000., 2500000. /)           

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: Rroot_soil_mtc  = &        !! Soil to Roots Water resistance involved in the calculation of psi_leaf
  & (/undef,    150.,    150.,    150.,    150.,   150.,   150.,  &     !! See Table 1 of Liu et al. (2015)
  &    150.,    150.,    150.,    150.,    150.,    150.  /)           

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: Rsto_leaf_mtc  = &         !! Leaf water storage Water resistance involved in the calculation of psi_leaf
  & (/undef, 1000000., 1000000., 1000000., 1000000., 1000000., 1000000.,  &           !! See Table 1 of Tuzet et al. (2017)
  &  1000000., 1000000., 1000000., 1000000., 1000000., 1000000.  /)           

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: Rsto_wood_mtc  = &         !! Wood water storage Water resistance involved in the calculation of psi_leaf
  & (/undef, 60000000., 60000000., 60000000., 60000000., 60000000., 60000000.,  &           !! See Table 1 of Tuzet et al. (2017)
  &  60000000., 60000000., 60000000., 60000000., 60000000., 60000000.  /)           

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: wood_wat_store_fac_mtc  = &   !! Empirical factor to increase/decrease wood water storage
  & (/undef,    2.18,   7.35,   6.79,   2.25,  8.5,   8.34,   &     !! 
  &    2.38,    5.0,    9.07,   9.35,   0.8,   0.9  /)
 
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: leaf_wat_store_fac_mtc  = &   !! Empirical factor to increase/decrease leaf water storage
  & (/undef,    5.84,   6.89,   3.9,    7.6,   1.82,  6.87,  &     !! 
  &    8.24,    6.0,    4.07,   6.0,    2.94,  6.15 /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: psi_ref_g0_mtc  = &        !! Reference water potential involved in the calculation of f_psi in the calculation of stomatal conductance
  & (/undef,    -1.97,   -1.05,   -1.17,   -0.81,   -1.67,   -1.15,  &     !! See Table 1 of Tuzet et al. (2017)
  &     -1.46,   -1.6,  -1.31,   -1.6,   -1.45,   -1.53  /) 
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: psi_ref_sto_leaf_mtc  = &  !! Reference water potential involved in the calculation of leaf storage water potential
  & (/undef,    -3.0,    -3.0,    -3.0,    -3.0,   -3.0,   -3.0,  &     !! See Table 1 of Tuzet et al. (2017)
  &    -3.0,    -3.0,    -3.0,    -3.0,    -3.0,    -3.0  /)           

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: psi_ref_sto_wood_mtc  = &  !! Reference water potential involved in the calculation of wood storage water potential
  & (/undef,    -1.4,    -1.4,    -1.4,    -1.4,   -1.4,   -1.4,  &     !! See Table 1 of Tuzet et al. (2017)
  &    -1.4,    -1.4,    -1.4,    -1.4,    -1.4,    -1.4  /)           

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: psi_50_res_leaf_mtc  = &   !! Reference water potential involved in the calculation of leaf water resistance
  & (/undef,    -1.1,    -1.1,    -1.1,    -1.1,   -1.1,   -1.1,  &     !! See Parameterisation of Yao et al. (2022)
  &    -1.1,    -1.1,    -1.1,    -1.1,    -1.1,    -1.1  /)           

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: psi_50_res_wood_mtc  = &   !! Reference water potential involved in the calculation of wood water resistance
  & (/undef,    -1.3,    -1.3,    -1.3,    -1.3,   -1.3,   -1.3,  &     !! See Parameterisation of Yao et al. (2021) 
  &    -1.3,    -1.3,    -1.3,    -1.3,    -1.3,    -1.3  /)           
  !& (/ undef,    -3.0,    -3.0,    -3.52,   -3.83,  -2.4,   -2.8, &              !! conductivity through cavitation. @tex $(MPa)$ @endtex
   !     -3.15,   -3.66,   undef,   undef,   undef,   undef /)
  
 !REAL(r_std), PARAMETER, DIMENSION(nvmc) :: psi_50_res_wood_mtc  = &   !! IF USED SHOULD BE WITH PSI_50_RES OF LEAF AND ROOT ADAPTED
  !& (/ undef,    -2.38,    -2.34,    -4.89,   -2.59,  -2.18,   -3.78, &  
   !     -1.75,   -3.92,   undef,   undef,   undef,   undef /)          !! derived from XFT database analysis 

 !REAL(r_std), PARAMETER, DIMENSION(nvmc) :: psi_ref_g0_mtc  = &       
  !& (/ undef,    -1.46,    -1.39,    -1.72,   -1.57,  -1.46,   -1.41, &!! psi50gs derived from the combination of sf and psi90, both being functions of the XFT psi50 of embolism.
   !     -1.31,   -1.18,   undef,   undef,   undef,   undef /)          !! The functions are derived from MartinStPaul2016. Biggest uncertainty (increased by the exponential) from sf
 
 !REAL(r_std), PARAMETER, DIMENSION(nvmc) :: sf_mtc  = &               
  !& (/ undef,     8.54,     6.38,     3.37,    7.31,   9.88,    2.88, &  !! sf calculated from XFT data (psi50_wood and sf) and an exponential best fit inspired from MartinStPaul2016
   !      9.83,    2.17,   undef,   undef,   undef,   undef /)          !!    however the exponential then averaged seems to give a lot of weight to outliers
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: psi_50_res_root_mtc  = &   !! Reference water potential involved in the calculation of root water resistance
  & (/undef,    -1.1,    -1.1,    -1.1,    -1.1,   -1.1,   -1.1,  &     !! See Parameterisation of Yao et al. (2022)
  &    -1.1,    -1.1,    -1.1,    -1.1,    -1.1,    -1.1  /)           

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: a_leaf_mtc  = &            !! Parameter involved in the calculation of leaf water resistance
  & (/undef,    -2.5,    -2.5,    -2.5,    -2.5,   -2.5,   -2.5,  &     !! See Parameterisation of Yao et al. (2022)
  &    -2.5,    -2.5,    -2.5,    -2.5,    -2.5,    -2.5  /)           

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: a_wood_mtc  = &            !! Parameter involved in the calculation of wood water resistance
  & (/undef,    -3.0,    -3.0,    -3.0,    -3.0,   -3.0,   -3.0,  &     !! See Parameterisation of Yao et al. (2022)
  &    -3.0,    -3.0,    -3.0,    -3.0,    -3.0,    -3.0  /)           

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: a_root_mtc  = &            !! Parameter involved in the calculation of root water resistance
  & (/undef,    -2.3,    -2.3,    -2.3,    -2.3,   -2.3,   -2.3,  &     !! See Parameterisation of Yao et al. (2022)
  &    -2.3,    -2.3,    -2.3,    -2.3,    -2.3,    -2.3  /)           

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: k_leaf_max_mtc  = &           !! Maximum leaf water conductance involved in the calculation of leaf water resistance
  & (/undef,0.000000525,0.000000525,0.000000525, 0.000000525, 0.000000525, 0.000000525,  &     !! See Parameterisation of Yao et al. (2022)
  &  0.000000525,0.000000525, 0.000000525,0.000000525,0.000000525,0.000000525  /)           

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: k_wood_max_mtc  = &           !! Maximum wood water conductance involved in the calculation of wood water resistance
  & (/undef,0.000000175,0.000000175,0.000000175, 0.000000175, 0.000000175, 0.000000175,  &     !! See Parameterisation of Yao et al. (2022)
  &  0.000000175,0.000000175, 0.000000175,0.000000175,0.000000175,0.000000175  /)      

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: sapwood_density_mtc = &                  !! Sapood density @tex $(gC.m^{-3})$ @endtex
  &(/undef,   5.3e5,  5.3e5,  5.3e5,  6.3e5,  5.9e5,  4.9e5,  4.1e5,  5.2e5,  &       !! From sapwood specific gravity (dry mass/green volume) of XFT database classified by PFT. 
  & undef,    undef,  undef,  undef/)                                                 !! By lack of PFT3 we took PFT2 value for it
                                                                                     
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: ks_max_mtc  = &               !! Maximum sapwood conductivity involved in the calculation of wood water resistance (m3/m/MPa/s)
  & (/undef,  0.00357,  0.00613,  0.00132,  0.00327,  0.00327,  0.00130,  0.00268,       &
  &  0.00030, undef,    undef,    undef,    undef/)                        !! XFT average per PFT, except PFT3 from Janssen et al. (2020) values

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: k_root_max_mtc  = &           !! Maximum root water conductance involved in the calculation of root water resistance
  & (/undef,0.000000175,0.000000175,0.000000175, 0.000000175, 0.000000175, 0.000000175,  &     !! See Parameterisation of Yao et al. (2022)
  &  0.000000175,0.000000175, 0.000000175,0.000000175,0.000000175,0.000000175  /)           

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: lambda_leaf_mtc  = &           !! Parameter involved in the calculation of leaf storage water potential
  & (/undef,    1.3,    1.3,    1.3,    1.3,   1.3,   1.3,  &               !!  See Table 1 of Tuzet et al. (2017)
  &    1.3,    1.3,    1.3,    1.3,    1.3,    1.3  /)           

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: lambda_wood_mtc  = &           !! Parameter involved in the calculation of wood storage water potential
  & (/undef,   0.9,     0.9,     0.9,     0.9,    0.9,    0.9,  &           !!  See Table 1 of Tuzet et al. (2017)
  &   0.9,     0.9,     0.9,     0.9,     0.9,     0.9  /)           

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: d_l_mtc  = &                !! PArameter for gs sensitivity to VPD (kPa) from Anderegg et al. (2017)
  & (/undef,   3.0,     3.0,     3.0,     3.0,     3.0,     3.0, &   
  &   3.0,     3.0,     3.0,     3.0,     3.0,     3.0  /)           
  
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: g0_mtc  = &                 !! Residual stomatal conductance when irradiance approaches zero (mol CO2 m−2 s−1 bar−1)
  & (/undef, 0.00625, 0.00625, 0.00625, 0.00625, 0.00625, 0.00625,  &   !! Value from ORCHIDEE - No other reference.
  & 0.00625, 0.00625, 0.00625, 0.01875, 0.00625, 0.01875  /)            !! modofy to account for the conversion for conductance to H2O to CO2 
  !& (/undef,  0.0049,  0.0049,  0.0049,  0.0049,  0.0049, 0.0049,  &     !! For tree PFTs, from Duursma et al. (2019)(n=221) Other values from ORCHIDEE - No other reference. 
  !& 0.0049,   0.0049,  0.0049,  0.01875, 0.00625, 0.01875  /)            !! modify to account for the conversion for conductance to H2O to CO2 

  !& (/undef, 0.00625, 0.00625, 0.00625, 0.00625, 0.00625, 0.00625,  &   !! Value from ORCHIDEE - No other reference.
  !& 0.00625, 0.00625, 0.00625, 0.01875, 0.00625, 0.01875  /)            !! modofy to account for the conversion for conductance to H2O to CO2 

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: h_protons_mtc  = &         !! Number of protons required to produce one ATP (mol mol-1)
  & (/undef,      4.,      4.,      4.,      4.,      4.,    4.,  &     !! See Table 2 of Yin et al. (2009) - h parameter
  &      4.,      4.,      4.,      4.,      4.,      4.  /)           

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: fpsir_mtc = &              !! Fraction of PSII e− transport rate 
  & (/undef,   undef,   undef,   undef,   undef,  undef,  undef,  &     !! partitioned to the C4 cycle (-)
  &   undef,   undef,   undef,     0.4,   undef,    0.4  /)             !! See Table 2 of Yin et al. (2009) - x parameter        
 
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: fQ_mtc = &                 !! Fraction of electrons at reduced plastoquinone 
  & (/undef,   undef,   undef,   undef,   undef,  undef,  undef,  &     !! that follow the Q-cycle (-) - Values for C3 platns are not used
  &   undef,   undef,   undef,      1.,   undef,     1.  /)             !! See Table 2 of Yin et al. (2009)          

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: fpseudo_mtc = &            !! Fraction of electrons at PSI that follow 
  & (/undef,   undef,   undef,   undef,   undef,  undef,  undef,  &     !! pseudocyclic transport (-) - Values for C3 platns are not used
  &   undef,   undef,   undef,     0.1,   undef,    0.1  /)             !! See Table 2 of Yin et al. (2009)    

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: kp_mtc = &                 !! Initial carboxylation efficiency of the PEP carboxylase (mol m−2 s−1 bar−1) 
  & (/undef,   undef,   undef,   undef,   undef,  undef,  undef,  &     !! See Table 2 of Yin et al. (2009)
  &   undef,   undef,   undef,     0.7,   undef,    0.7  /)                 

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: alpha_mtc = &              !! Fraction of PSII activity in the bundle sheath (-)
  & (/undef,   undef,   undef,   undef,   undef,  undef,  undef,  &     !! See legend of Figure 6 of Yin et al. (2009)
  &   undef,   undef,   undef,     0.1,   undef,    0.1  /)                 

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: gbs_mtc = &                !! Bundle-sheath conductance (mol m−2 s−1 bar−1)
  & (/undef,   undef,   undef,   undef,   undef,  undef,  undef,  &     !! See legend of Figure 6 of Yin et al. (2009)
  &   undef,   undef,   undef,   0.003,   undef,  0.003  /)    

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: theta_mtc = &              !! Convexity factor for response of J to irradiance (-)
  & (/undef,     0.7,     0.7,     0.7,     0.7,    0.7,    0.7,  &     !! See Table 2 of Yin et al. (2009) 
  &     0.7,     0.7,     0.7,     0.7,     0.7,    0.7  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: alpha_LL_mtc = &           !! Conversion efficiency of absorbed light into J at strictly limiting light (mol e− (mol photon)−1)
  & (/undef,     0.53,    0.53,    0.3,     0.3,    0.3,    0.3,  &     !! See comment from Yin et al. (2009) after eq. 4
  &     0.3,     0.3,     0.3,     0.3,     0.3,    0.3  /)             !! alpha value from Medlyn et al. (2002)   
                                                                        !! 0.093 mol CO2 fixed per mol absorbed photons
                                                                        !! times 4 mol e- per mol CO2 produced

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: stress_vcmax_mtc = &       !! Water stress on vcmax
  & (/    1.,     1.,     1.,       1.,      1.,     1.,      1., &
  &      1.,     1.,     1.,       1.,      1.,     1.  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: stress_gs_mtc = &          !! Water stress on gs
  & (/    1.,     1.,     1.,       1.,      1.,     1.,      1., &
  &      1.,     1.,     1.,       1.,      1.,     1.  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: stress_gm_mtc = &          !! Water stress on gm
  & (/    1.,     1.,     1.,       1.,      1.,     1.,      1., &
  &      1.,     1.,     1.,       1.,      1.,     1.  /)
    
  !-
  ! 2 .Stomate
  !-
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: ext_coeff_mtc = &             !! extinction coefficient of the Monsi&Saeki
  & (/ 0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,  &                   !! relationship (1953) ((m2[ground]) (m-2[leaf])) 
  &    0.5,   0.5,   0.5,   0.5,   0.5,   0.5  /)
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: ext_coeff_vegetfrac_mtc = &   !! extinction coefficient used for defining the fraction
  & (/ 1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &                   !!  of bare soil (unitless)
  &    1.0,   1.0,   1.0,   1.0,   1.0,   1.0  /)
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: ext_coeff_N_mtc = &           !! extinction coefficient of the leaf N content profile within the canopy
  & (/ 0.15,  0.15,  0.15,  0.15,  0.15,  0.15,  0.15,  &                  !! ((m2[ground]) (m-2[leaf]))
  &    0.15,  0.15,  0.15,  0.15,  0.15,  0.15  /)                         !! based on Dewar et al. (2012, value of 0.18), on Carswell et al. (2000, value of 0.11 used in OCN)
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: nue_opt_mtc = &               !! Nitrogen use efficiency of Vcmax 
  & (/ undef,   14.08,    30.0,   10.5,  18.88,  56.61,  10.5, &           !! ((mumol[CO2] s-1) (gN[leaf])-1)
  &    45.00,   20.79,     45.,    45.,    60.,    60.     /)              !! based on the work of Kattge et al. (2009, GCB)
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: vmax_uptake_nh4_mtc = &       !! Vmax of nitrogen uptake by plants for ammonium (umol (g DryWeight_root)-1 h-1)). 
  & (/ undef,   13.6,     12.0,  13.66,  23.96,   10.18,   10.65, &        !! from  Kronzucker et al. (1995, 1996) but externalized and tuned for ORC4
  &    14.06,   16.51,     9.,    9.,     9.,    9. /)
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: vmax_uptake_no3_mtc = &       !! Vmax of nitrogen uptake by plants for nitrate (umol (g DryWeight_root)-1 h-1))
  & (/ undef,   13.13,  12.0,   4.854,   5.03,   23.51,   6.608, &         !! from Zaehle & Friend (2010) but externalized and tuned for ORC4
  &     10.6,   17.97,    9.,    9.,    9.,      9. /)

    
  !
  ! RESPIRATION (stomate)
  !
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: frac_growthresp_mtc  =  &              !! fraction of GPP which is lost as growth respiration
  & (/   0.28,   0.28,   0.28,   0.28,   0.28,   0.28,   0.28,  &
  &      0.28,   0.28,   0.28,   0.28,   0.28,   0.28  /)

   REAL(r_std), PARAMETER, DIMENSION(nvmc) :: maint_resp_slope_c_mtc  =  &          !! slope of maintenance respiration coefficient (1/K),
  & (/  undef,   0.44, 0.6, 0.2, 0.6, 0.3, 0.3, &                                   !! constant c of aT^2+bT+c, tabulated
  &     0.4,     0.4,  0.5, 0.5, 0.5, 0.7  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: maint_resp_slope_b_mtc  =  &           !! slope of maintenance respiration coefficient (1/K),
  & (/  undef,   -0.01, -0.01, -0.04, -0.02, -0.02, -0.01, &                        !! constant b of aT^2+bT+c, tabulated
  &     -0.01,   -0.01, -0.01, -0.0,  -0.0, -0.01  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: maint_resp_slope_a_mtc  =  &           !! slope of maintenance respiration coefficient (1/K),
  & (/  undef,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  &                         !! constant a of aT^2+bT+c, tabulated
  &       0.0,   0.0,   0.0,   0.0,   0.0,   0.0  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: maint_resp_slope_c_cn_mtc  =  &        !! slope of maintenance respiration coefficient (1/K),
  & (/  undef,   0.12,   0.12,   0.16,   0.16,   0.16,   0.25,  &                   !! constant c of aT^2+bT+c, tabulated
  &      0.25,   0.25,   0.16,   0.12,   0.16,   0.12  /)
  
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: maint_resp_slope_b_cn_mtc  =  &        !! slope of maintenance respiration coefficient (1/K),
  & (/  undef,   0.0,        0.0,   0.0,        0.0,   0.0,   0.0,  &               !! constant b of aT^2+bT+c, tabulated
  &       0.0,   0.0,   -0.00133,   0.0,   -0.00133,   0.0  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: coeff_maint_init_mtc  =   &            !! maintenance respiration coefficient
  & (/  undef,  0.04531,  0.04,  0.04615,  0.0227,  0.226,  0.022, &              !! at 10 deg C - from Sitch et al. 2003 and Zaehle (OCN)
  &    0.09333,  0.2049,   0.06,  0.05,     0.05,     0.08  /)                       !! @tex $(gC.gN^{-1}.day^{-1})$ @endtex

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: coeff_maint_init_cn_mtc  =   &         !! maintenance respiration coefficient
       & (/   undef,   7.0E-2,   7.0E-2,   7.0E-2,   7.0E-2,   7.0E-2,  7.0E-2,  &  !! at 10 deg C - from Sitch et al. 2003 and Zaehle (OCN)
       &    7.0E-2,   7.0E-2,   7.0E-2,   7.0E-2,   7.0E-2,   7.0E-2  /)            !! @tex $(gC.gN^{-1}.day^{-1})$ @endtex

  
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: tref_maint_resp_mtc  =   &             !! maintenance respiration Temperature coefficient (deg C)
  & (/   undef,     10.02,     10.02,     10.02,     10.02,     10.02,     10.02,  &  
  &      10.02,     10.02,     10.02,     10.02,     10.02,     10.02  /)             
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: tmin_maint_resp_mtc  =   &             !! maintenance respiration Temperature coefficient (deg C)
  & (/   undef,    -46.02,    -46.02,    -46.02,    -46.02,    -46.02,    -46.02,  &  
  &     -46.02,    -46.02,    -46.02,    -46.02,    -46.02,    -46.02  /)             
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: e0_maint_resp_mtc  =   &               !! maintenance respiration Temperature coefficient (unitless)
  & (/   undef,   308.56,     308.56,   308.56,     308.56,   308.56,     308.56,  &  
  &     308.56,   308.56,     308.56,   308.56,     308.56,   308.56  /) 


  !
  ! Allocation (stomate)
  !
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: tref_labile_mtc  =   &                  !! Growth from labile pool - temperature at which all labile C will be allocated to growth (deg C)
  & (/   undef,   5.0,      5.0,      5.0,      5.0,      5.0,      5.0,  &  
  &       5.0,    5.0,      5.0,      5.0,      5.0,      5.0  /)             
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: tmin_labile_mtc  =   &                  !! Growth from labile pool  - temperature above which labile will be allocated to growth (deg C)
  & (/   undef,    -2.0,    -2.0,    -2.0,    -2.0,    -2.0,    -2.0,  &  
  &      -2.0,     -2.0,    -2.0,    -2.0,    -2.0,    -2.0  /) 
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: e0_labile_mtc  =   &                  !! Growth temperature coefficient - tuned see stomate_growth_fun_all.f90 (unitless)
  & (/   undef,   15.0,     15.0,   15.0,     15.0,   15.0,     15.0,  &  
  &      15.0,    15.0,     15.0,   15.0,     15.0,   15.0  /)
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: always_labile_mtc = &                 !! share of the labile pool that will remain in the labile pool (unitless)
  & (/   undef,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01, &
  &      0.01,   0.01,  0.01,  0.01,  0.01/)

  !
  ! SOM decomposition (stomate)
  !
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: LC_leaf_mtc = &        !! Lignin/C ratio of leaf pool (unitless)
  & (/   0.15,   0.18,   0.18,   0.24,   0.18,   0.18,   0.24,  &   !! based on CN from White et al. (2000)       
  &      0.18,   0.24,   0.09,   0.09,   0.09,   0.09  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: LC_sapabove_mtc = &    !! Lignin/C ratio of sapabove pool (unitless)
  & (/   0.15,   0.23,   0.23,   0.29,   0.23,   0.23,   0.29,  &   !! based on CN from White et al. (2000)       
  &      0.23,   0.29,   0.09,   0.09,   0.09,   0.09  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: LC_sapbelow_mtc = &    !! Lignin/C ratio of sapbelow pool (unitless)
  & (/   0.15,   0.23,   0.23,   0.29,   0.23,   0.23,   0.29,  &   !! based on CN from White et al. (2000)       
  &      0.23,   0.29,   0.09,   0.09,   0.09,   0.09  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: LC_heartabove_mtc = &  !! Lignin/C ratio of heartabove pool (unitless)
  & (/   0.15,   0.23,   0.23,   0.29,   0.23,   0.23,   0.29,  &   !! based on CN from White et al. (2000)       
  &      0.23,   0.29,   0.09,   0.09,   0.09,   0.09  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: LC_heartbelow_mtc = &  !! Lignin/C ratio of heartbelow pool (unitless)
  & (/   0.15,   0.23,   0.23,   0.29,   0.23,   0.23,   0.29,  &   !! based on CN from White et al. (2000)       
  &      0.23,   0.29,   0.09,   0.09,   0.09,   0.09  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: LC_fruit_mtc = &       !! Lignin/C ratio of fruit pool (unitless)
  & (/   0.09,   0.09,   0.09,   0.09,   0.09,   0.09,   0.09,  &   !! based on CN from White et al. (2000)       
  &      0.09,   0.09,   0.09,   0.09,   0.09,   0.09  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: LC_root_mtc = &        !! Lignin/C ratio of root pool (unitless)
  & (/   0.22,   0.22,   0.22,   0.22,   0.22,   0.22,   0.22,  &   !! based on CN from White et al. (2000)       
  &      0.22,   0.22,   0.22,   0.22,   0.22,   0.22  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: LC_carbres_mtc = &     !! Lignin/C ratio of carbres pool (unitless)
  & (/   0.15,   0.18,   0.18,   0.24,   0.18,   0.18,   0.24,  &   !! based on CN from White et al. (2000)       
  &      0.18,   0.24,   0.09,   0.09,   0.09,   0.09  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: LC_labile_mtc = &      !! Lignin/C ratio of labile pool (unitless)
  & (/   0.15,   0.18,   0.18,   0.24,   0.18,   0.18,   0.24,  &   !! based on CN from White et al. (2000)       
  &      0.18,   0.24,   0.09,   0.09,   0.09,   0.09  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: decomp_factor_mtc  =  &  !! Multpliactive factor modifying the standard decomposition factor for each SOM pool
  & (/     1.,     1.,     1.,     1.,     1.,     1.,     1.,  &          
  &        1.,     1.,     1.,     1.,    1.2,    1.4  /)
  !
  ! STAND STRUCTURE
  !
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: crown_to_height_mtc = &  !! Crown depth as a function of tree height (unitless)
  &(/ undef,  0.5,    0.5,   0.7,   0.5,   0.66,   0.75, &
  &     0.66,  0.66,  1.0, 1.0, 1.0, 1.0 /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: crown_vertohor_dia_mtc = & !! Crown diameter a function of crown depth (and thus tree height) (unitless).
  &(/ undef,  2.0,    2.0,   0.5,   1.5,   1.0,   0.4, &
  &     1.0,  1.0,  1.0, 1.0, 1.0, 1.0 /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: tree_ff_mtc = &                          !! Tree form factor to reduce
  &(/ undef, 0.4681,  0.4234,  0.7717, 0.793,  0.8552, 0.6036, &                           !! the volume of a cylinder
  &    0.9442,  0.6209,  1.0,     1.0,  1.0,  1.0 /)                                     !! to the volume of the real tree shape

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: pipe_density_mtc = &                     !! Wood density @tex $(gC.m^{-3})$ @endtex
  &(/   undef,   2.87458e5,   2.68753e5,   2.08333e5,   3.0e5,   2.38e5,   1.95e5,  & !! Current values are taken from the trunk. 
  &    2.38e5,  2.4875e5,     2.0e5,       2.0e5,       2.0e5,   2.0e5  /)            !! forestry-branch has more realistic values 
                                                                                      !! in it. Source: AFOCEL 2006

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: pipe_tune1_mtc = &                       !! cn_area = pipe_tune1*...
  &(/ undef,  undef,  undef,  undef,  undef,  undef,  undef, &                        !!    stem diameter**pipe_tune_exp_coeff
  &    undef,  undef, undef, undef, undef, undef /)                                   !! for consistency reason pipe_tune1 is calculated as 0.66*pi/4*pipe_tune2**2
                                                                                      !! see pft_parameters.f90

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: pipe_tune2_mtc = &                       !! height=pipe_tune2 * diameter**pipe_tune3
  &(/ undef,   55.,   55.,   45.,   14.,   50.,   30., &
  &     30.,   30., undef, undef, undef, undef /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: min_pipe_tune2_mtc = &                    !! initial value for pipe_tune2. pipe_tune2 is calculated as 
  &(/ undef,   7.,    7.,    7.,    7.,    7.,   7., &                                 !! function of precipitation 
  &     7.,    7., undef, undef, undef, undef /)
 
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: pipe_tune3a_mtc = &                       !! height=pipe_tune2 * diameter**pipe_tune3
  &(/ undef,   0.45,   0.45,   0.45,   0.45,   0.45,  0.50, &                          !! when 'test_dynamic_max_height'
  &    0.52,   0.52,  undef,  undef,  undef,  undef /)
  
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: pipe_tune3b_mtc = &                       !! height=pipe_tune2 * diameter**pipe_tune3
  &(/ undef,   0.65,   0.65,   0.57,   0.33,   0.66,  0.58, &
  &    0.52,   0.52,  undef,  undef,  undef,  undef /)
      
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: pipe_tune4_mtc = &                       !! CHECK - needed for stem diameter
  &(/ undef,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3, &
        0.3,   0.3, undef, undef, undef, undef /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: pipe_k1_mtc = &                          !! CHECK 
  &(/ undef,  8.e3,  8.e3,  8.e3,  8.e3,  8.e3,  8.e3, &
  &    8.e3,  8.e3, undef, undef, undef, undef /) 

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: pipe_tune_exp_coeff_mtc = &              !! cn_area = pipe_tune1*...  
  &(/ undef,   undef,   undef,   undef,   undef,   undef,   undef, &                  !!    stem diameter**pipe_tune_exp_coeff
  &     undef,   undef, undef, undef, undef, undef /)                                 !! for consistency reasons pipe_tune_exp_coeff is calculated as 2*pipe_tune3
                                                                                      !! see pft_parameters.f90
    
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: mass_ratio_heart_sap_mtc = &             !! mass ratio (heartwood+sapwood)/heartwood 
  &(/ undef,    3.,    3.,    3.,    3.,    3.,    3., &
  &      3.,    3., undef, undef, undef, undef /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: canopy_cover_mtc = &                     !! Prescribed canopy cover (1-gap fraction) 
  & (/ undef,    0.9,   0.9,   0.7,   0.7,   0.7,   0.6, &                            !! of a canopy (unitless)
  &      0.5,    0.5,   0.9,   0.9,   0.9,   0.9 /)  

  INTEGER(i_std), PARAMETER, DIMENSION(nvmc) :: nmaxplants_mtc = &                    !! Initial number of trees per ha. This parameter is 
  & (/ undef,  65000.,  65000.,  65000.,  65000.,  65000., 65000.,  &                 !! used at .firstcall. and after clearcuts
  &    65000.,  65000.,  10000.,  10000.,  10000.,  10000. /)                         !! the value is used by the allometric allocation
                                                                                      !! and forestry subroutines.

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: height_init_mtc = &                      !! The height (m) of a grass or crop when the vegetation is established. 
  &(/ undef,  undef,  undef,  undef,  undef,  undef,  undef, &                        !! In combination with the parameter lai_to_height and the allometric relationships
  &   undef,  undef,   0.3,   0.3,   0.2,   0.2 /)                                    !! this setting determines all biomass components of a newly establised vegetation

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: biomass_init_mtc = &                     !! The biomass (g) when the plant is established.          
  &(/ undef,  undef,  undef,  undef,  undef,  undef,  undef, &                        !! Must be used with OK_BIOMASS_INIT flag.        
  &   undef,  undef,   20.,   20.,   70.,   80. /)                                        !! An alternative to using height_init and lai_to_height. 
  
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: qmd_init_mtc = &                         !! The minimum (aboveground) diameter (m) of saplings when a forest 
  &(/ undef,    0.025,    0.025,   0.025,    0.025,   0.025,   0.025, &                     !! stand is established. 
  &    0.025,    0.025,    undef,  undef,   undef,  undef /)
                                            
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: dia_thresh_inv_mtc = &                   !! Diameter (m) below wich trees are not measured in forest  
  &(/ undef,    0.07,    0.07,    0.07,    0.07,   0.07,    0.07, &                   !! inventories. This value only affects the history variables
  &    0.07,    0.07,    undef,   undef,   undef,  undef /)                           !! dia_inv and height_inv

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: alpha_self_thinning_mtc = &              !! Coefficient of the self-thinning relationship D=alpha*N^beta 
  &(/ undef,  2827.,  3700.,  880.,  2430.,  2000.,  2507., &                         !! estimated from German, French, Spanish and Swedish 
  &     770.,  800.,  undef,  undef,  undef,  undef/)                                 !! forest inventories
 
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: beta_self_thinning_mtc = &               !! Exponent of the self-thinning relationship D=alpha*N^beta 
  &(/ undef,  -0.73,  -0.67,  -0.57,  -0.69,  -0.67,  -0.73, &                        !! estimated from German, French, Spanish and Swedish 
  &    -0.59, -0.59,  undef,  undef,  undef,  undef/)                                 !! forest inventories

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: ref_alpha_self_thin_mtc = &              !! Basic change rate of alpha self thinning coefficient.
  &(/ undef,  2827.,  3700.,  880.,  1300.,  1220.,  4050., &                        
  &    1000.,  800.,  undef,  undef,  undef,  undef/)                                  

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: init_pre_indust_ref_gpp_mtc = &          !! The pre-industrial reference gpp that is used in the
  &(/ undef,  8.21,    8.21,   4.11,   3.40,   4.11,   2.74, &                        !! calculations of the temporal and spatial dynamics of 
  &    2.47,  2.47,   undef,  undef,  undef,  undef/)                                 !! alpha_self_thinning. 
  
  !
  ! RECRUITMENT (stomate)  
  !                                                                       
  LOGICAL, PARAMETER, DIMENSION(nvmc) :: recruitment_pft_mtc = &                      !! Do recruitment? (true/false)  
  & (/ .FALSE.,   .TRUE.,   .TRUE.,    .TRUE.,    .TRUE.,    .TRUE.,   .TRUE.,  &  
  &    .TRUE.,    .TRUE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE. /)    

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: recruitment_height_mtc = &               !! Prescribed height (above) for the recruited stems @tex $(m)$ @endtex  
  & (/ undef, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, &
  &    undef, undef, undef, undef/)        

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: recruitment_alpha_mtc = &                !! Alpha parameter for the power function used to model recruitment 
  & (/ undef, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, &                       !! from light and density @tex $(unitless)$ @endtex. 
  &    undef, undef, undef, undef/)                                                   

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: recruitment_beta_mtc = &                 !! Beta parameter for the power function used to model recruitment 
  &(/ undef, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, &        !! from light and density @tex $(unitless)$ @endtex
  &   undef, undef, undef, undef/)

  
  ! BARK BEETLES (ok_pest = y)
  LOGICAL, PARAMETER, DIMENSION(nvmc) :: beetle_pft_mtc = &                            !! Do bark beetle attack? (true/false)
  & (/ .FALSE.,   .FALSE.,   .FALSE.,   .TRUE.,   .FALSE.,   .FALSE.,.TRUE.,  &
  &    .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE. /)

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: S_activity_mtc  =  &                       !! 
  & (/  undef,   undef,   undef,   -500.0,   undef,   undef,   -500.0, & 
  &      undef,   undef,   undef,  undef,  undef,  undef  /)

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: act_limit_mtc  =  &                        !! 
  & (/  undef,   undef,   undef,   0.06,   undef,   undef,   0.06, & 
  &      undef,   undef,   undef,  undef,  undef,  undef  /)

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: S_competition_mtc  =  &                    !! 
  & (/  undef,   undef,   undef,   -5.0,   undef,   undef,   -5.0, &
  &      undef,   undef,   undef,  undef,  undef,  undef  /)

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: S_susceptibility_mtc  =  &                       !! 
  & (/  undef,   undef,   undef,   -5.0,   undef,   undef,   -5.0, &
  &      undef,   undef,   undef,  undef,  undef,  undef  /)

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: i_rd_susceptibility_mtc  =  &                     !! 
  & (/  undef,   undef,   undef,   0.6,   undef,   undef,   0.6, & 
  &      undef,   undef,   undef,  undef,  undef,  undef  /)

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: S_share_mtc  =  &                          !! 
  & (/  undef,   undef,   undef,  0.55,   undef,   undef,   0.55, & 
  &      undef,   undef,   undef,  undef,  undef,  undef  /)

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: SH_limit_mtc  =  &                         !! 
  & (/  undef,   undef,   undef,   15.5,   undef,   undef,   15.5, & 
  &      undef,   undef,   undef,  undef,  undef,  undef  /)

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: S_defence_mtc  =  &                        !! 
  & (/  undef,   undef,   undef,   -9.5,   undef,   undef,   -9.5, & 
  &      undef,   undef,   undef,  undef,  undef,  undef  /)

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: PWS_limit_mtc  =  &                        !! 
  & (/  undef,   undef,   undef,   0.4,   undef,   undef,   0.4, &  
  &      undef,   undef,   undef,  undef,  undef,  undef  /)

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: beetle_generation_a_mtc  =&                !! a parameter for the calculation of the number of beetle generation per year 
  & (/  undef,   undef,   undef,   1.0,   undef,   undef, 1.0, & 
  &      undef,   undef,   undef,  undef,  undef,  undef  /)

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: beetle_generation_b_mtc  =&               !! b parameter for the calculation of the number of beetle generation per year 
  & (/  undef,   undef,   undef, 557.0,   undef,   undef, 557.0, & 
  &      undef,   undef,   undef,  undef,  undef,  undef  /)

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: beetle_generation_c_mtc  =&               !! c parameter for the calculation of the number of beetle generation per year
  & (/  undef,   undef,   undef, 1.0,   undef,   undef, 1.0, & 
  &      undef,   undef,   undef,  undef,  undef,  undef  /)

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: min_temp_beetle_mtc  =&                   !! temperature threshold below which Teff is not calculated (*C) 
  & (/  undef,   undef,   undef, 8.3,   undef,   undef, 8.3, & 
  &      undef,   undef,   undef,  undef,  undef,  undef  /)

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: max_temp_beetle_mtc  =&                   !! temperature threshold above which Teff is not calculated (*C) 
  & (/  undef,   undef,   undef, 38.4,   undef,   undef, 38.4, & 
  &      undef,   undef,   undef,  undef,  undef,  undef  /)

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: opt_temp_beetle_mtc  =&                   !! optimal temperature to breed bark beetle (*C) 
  & (/  undef,   undef,   undef, 30.3,   undef,   undef, 30.3, &
  &      undef,   undef,   undef,  undef,  undef,  undef  /)

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: eff_temp_beetle_a_mtc  =&                 !! a parameter for the calculation of the effective temperature used in beetle phenology
  & (/  undef,   undef,   undef,  0.02876507,   undef,   undef, 0.02876507, & 
  &      undef,   undef,   undef,  undef,  undef,  undef  /)

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: eff_temp_beetle_b_mtc  =&                 !! b parameter for the calculation of the effective temperature used in beetle phenology
  & (/  undef,   undef,   undef, 40.9958913,   undef,   undef, 40.9958913, & 
  &      undef,   undef,   undef,  undef,  undef,  undef  /)

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: eff_temp_beetle_c_mtc  =&                 !! c parameter for the calculation of the effective temperature used in beetle phenology
  & (/  undef,   undef,   undef, 3.5922336,   undef,   undef, 3.5922336, &
  &      undef,   undef,   undef,  undef,  undef,  undef  /)

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: eff_temp_beetle_d_mtc  =&                 !! d parameter for the calculation of the effective temperature used in beetle phenology 
  & (/  undef,   undef,   undef, 1.24657367,   undef,   undef, 1.24657367, & 
  &      undef,   undef,   undef,  undef,  undef,  undef  /)

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: diapause_thres_daylength_mtc  =&          !! daylength in hour above which bark beetle start diapause
  & (/  undef,   undef,   undef,   14.5,   undef,   undef,   14.5, & 
  &      undef,   undef,   undef,  undef,  undef,  undef  /)

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: S_massattack_mtc  =&                      !!
  & (/  undef,   undef,   undef,   -30.0,   undef,   undef,   -30.0, & 
  &      undef,   undef,   undef,  undef,  undef,  undef  /)

  !! SCALE DEPENDENT BARK BEETLE PARAMETERS.
  ! The choice of the scale is linked to the resolution of the forcing files
  ! large scale: CRU_JRA, ISIMIP, LMDzOR, zoomed LMDzOR
  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: max_Nwood_degree_mtc  =  &                !! Value of Nwood at which ihosts dead = 1.0
  & (/  undef,   undef,   undef,   0.2,   undef,   undef,   0.2, & 
  &      undef,   undef,   undef,  undef,  undef,  undef  /)

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: BP_limit_degree_mtc  =&                   !! Beetle pressure at which beetles can mass attack
  & (/  undef,   undef,   undef,   0.12,   undef,   undef,   0.12, & 
  &      undef,   undef,   undef,  undef,  undef,  undef  /)

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: RDi_limit_degree_mtc  =  &                !! Relative density index at which ihosts competition = 0.5
  & (/  undef,   undef,   undef,    0.4,   undef,   undef,   0.4, & 
  &      undef,   undef,   undef,  undef,  undef,  undef  /)

  ! medium scale: SAFRAN
  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: max_Nwood_km2_mtc  =  &                   !! Value of Nwood at which ihosts dead = 1.0
  & (/  undef,   undef,   undef,   0.2,   undef,   undef,   0.2, & 
  &      undef,   undef,   undef,  undef,  undef,  undef  /)

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: BP_limit_km2_mtc  =&                      !! Beetle pressure at which beetles can mass attack
  & (/  undef,   undef,   undef,   0.12,   undef,   undef,   0.12, & 
  &      undef,   undef,   undef,  undef,  undef,  undef  /)

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: RDi_limit_km2_mtc  =  &                   !! Relative density index at which ihosts competition = 0.5
  & (/  undef,   undef,   undef,    0.4,   undef,   undef,   0.4, & 
  &      undef,   undef,   undef,  undef,  undef,  undef  /)

  ! small scale: FLUXNET, site level
  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: max_Nwood_ha_mtc  =  &                    !! Value of Nwood at which ihosts dead = 1.0
  & (/  undef,   undef,   undef,   0.2,   undef,   undef,   0.2, & 
  &      undef,   undef,   undef,  undef,  undef,  undef  /)

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: BP_limit_ha_mtc  =&                       !! Beetle pressure at which beetles can mass attack
  & (/  undef,   undef,   undef,   0.12,   undef,   undef,   0.12, & 
  &      undef,   undef,   undef,  undef,  undef,  undef  /)

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: RDi_limit_ha_mtc  =  &                    !! Relative density index at which ihosts competition = 0.5
  & (/  undef,   undef,   undef,    0.4,   undef,   undef,   0.4, & 
  &      undef,   undef,   undef,  undef,  undef,  undef  /)

  !
  ! WINDFALL (stomate)
  !

  ! At the moment, species-related parameters are scarce (more tree pulling tests are needed), therefore defining parameters for metaclasses
  ! is far from being straightforward.

  ! NOTE:  All the parameter values below were originated from six needle-leaf species:
  ! SS: Sitka Spruce 
  ! NS: Norway Spruce
  ! SP: Scots Pine
  ! LP: Lodgepople Pine
  ! CP: Corsican Pine  
  ! EL: European Larch
  !     and three broad leaved species:
  ! BI: Beech (Fagus)
  ! BE: Birch (Betula)
  ! OK: Oak   (Quercus)
  ! see the parameter table for the detail of combination of species to PFT

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: streamlining_c_leaf_mtc  =  &   !! Streamlining parameter. @tex $(unitless)$ @endtex. Streamlining is the
  & (/  undef,   2.34,   2.34,   2.70,   2.66,   2.34,   2.71, &            !! change of shape of the crowns due to wind. Used in the calculation of the
  &      2.15,   3.07,   undef,  undef,  undef,  undef  /)                   !! critical wind speed according to the GALES (Hale et al. 2015) model.
                                                                            !! In this case, the tree is in leaf.

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: streamlining_c_leafless_mtc = & !! Streamlining parameter. @tex $(unitless)$ @endtex. Streamlining is the
  & (/  undef,   2.34,   2.34,   2.70,   2.66,   2.34,   2.71, &            !! change of shape of the crowns due to wind. Used in the calculation of the
  &      2.15,   3.07,   undef,  undef,  undef,  undef  /)                   !! critical wind speed according to the GALES (Hale et al. 2015) model.
                                                                            !! In this case, the tree is leafless.

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: streamlining_n_leaf_mtc  =  &   !! Streamlining parameter. @tex $(unitless)$ @endtex. Streamlining is the
  & (/  undef,   0.88,   0.88,   0.64,   0.85,   0.88,   0.63, &            !! change of shape of the crowns due to wind. Used in the calculation of the
  &      0.88,   0.75,   undef,  undef,  undef,  undef  /)                   !! critical wind speed according to the GALES (Hale et al. 2015) model.
                                                                            !! In this case, the tree is in leaf.

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: streamlining_n_leafless_mtc = & !! Streamlining parameter. @tex $(unitless)$ @endtex. Streamlining is the
  & (/  undef,   0.88,   0.88,   0.64,   0.85,   0.88,   0.63, &            !! change of shape of the crowns due to wind. Used in the calculation of the
  &      0.88,   0.75,   undef,  undef,  undef,  undef  /)                   !! critical wind speed according to the GALES (Hale et al. 2015) model.
                                                                            !! In this case, the tree is leafless.

!  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: streamlining_rb_leaf_mtc  =  &  !! Streamlining parameter. @tex $(unitless)$ @endtex. Streamlining is the
!  & (/  undef,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, &                  !! change of shape of the crowns due to wind. Used in the calculation of the
!  &      0.0,   0.0,   undef,  undef, undef, undef  /)                         !! critical wind speed according to the GALES (Hale et al. 2015) model.
                                                                            !! In this case, the tree is in leaf.

!  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: streamlining_rb_leafless_mtc =& !! Streamlining parameter. @tex $(unitless)$ @endtex. Streamlining is the
!  & (/  undef,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, &                  !! change of shape of the crowns due to wind. Used in the calculation of the
!  &      0.0,   0.0,   undef,  undef, undef, undef  /)                         !! critical wind speed according to the GALES (Hale et al. 2015) model.
                                                                            !! In this case, the tree is leafless.

!  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: canopy_density_leaf_mtc  =  &   !! Density of the tree canopy @tex $(kg/m^{3})$ @endtex, i.e. of the 
!  & (/  undef,   2.5,   2.5,   2.5,   2.5,   2.5,   2.5, &                  !! branches and leaves. Used in the calculation of the critical wind speed
!  &      2.5,   2.5,   undef, undef,  undef, undef  /)                         !! according to the GALES (Hale et al. 2015) model.

!  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: canopy_density_leafless_mtc = & !! Density of the tree canopy @tex $(kg/m^{3})$ @endtex, i.e. of the 
!  & (/  undef,   2.5,   2.5,   2.5,   2.5,   2.5,   2.5, &                  !! branches (no leaves). Used in the calculation of the critical wind speed
!  &      2.5,   2.5,   undef, undef,  undef, undef  /)                         !! according to the GALES (Hale et al. 2015) model.

!  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: intercept_breadth_mtc  =  &      !! Intercept in the equation for calculating the canopy breadth.
!  & (/  undef,   0.5824,   0.5824,   0.5824,   0.5824,   0.5824,   0.5824, & !! @tex $(unitless)$ @endtex
!  &      0.5824,   0.5824,  undef,  undef,   undef,  undef /)

!  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: slope_breadth_mtc  =  &          !! Slope in the equation for calculating the canopy breadth.
!  & (/  undef,   0.115,   0.115,   0.115,   0.115,   0.115,   0.115, &       !! @tex $(unitless)$ @endtex
!  &      0.115,   0.115,  undef,  undef,  undef, undef  /)

!  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: intercept_depth_mtc  =  &        !! Intercept in the equation for calculating the canopy depth.
!  & (/  undef,   0.4206,   0.4206,   0.4206,   0.4206,   0.4206,   0.4206, & !! @tex $(unitless)$ @endtex
!  &      0.4206,   0.4206,  undef,   undef,   undef,   undef  /)

!  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: slope_depth_mtc  =  &            !! Slope in the equation for calculating the canopy depth.
!  & (/  undef,   0.4368,   0.4368,   0.4368,   0.4368,   0.4368,   0.4368, & !! @tex $(unitless)$ @endtex
!  &      0.4368,   0.4368,   undef,  undef,  undef,  undef /)

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: modulus_rupture_mtc  =  &  !! Modulus of rupture @tex $(Pa)$ @endtex. The measure of a species’ strength
  & (/  undef,   6.23E7, 6.23E7, 4.13E7, 5.90E7, 6.23E7,  4.10E7, &    !! before rupture when being bent. Used in the calculation of the critical
  &      6.27E7, 5.30E7,  undef,  undef,  undef,  undef  /)            !! wind speed according to the GALES (Hale et al. 2015) model. IMPORTANT:
                                                                       !! greenwood values are used and not the more frequently available drywood
                                                                       !! modulus of rupture.

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: f_knot_mtc  =  &           !! Knot Factor @tex $(unitless)$ @endtex. This modifier represents the knots
  & (/  undef,   1.0,   1.0,   0.87, 1.0,  1.0,   0.88, &              !! in the wood, and hence the decrease in structural strength. Used in the
  &      1.0,   0.85,   undef, undef, undef, undef /)                  !! calculation of the critical wind speed according to the GALES
                                                                       !! (Hale et al. 2015) model.

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: overturning_free_draining_shallow_mtc  =  &  !! Overturning moment multiplier @tex $(Nm/kg)$ @endtex.
  & (/  undef,   175.3,   175.3,   134.7,  198.5,   175.3,  132.6, &   !! This is derived from the generic soil type (free_draining mineral soils;
  &      152.0,  145.2,   undef,   undef,  undef,   undef /)           !! Gleyed mineral soils; Peaty mineral soils; Deep peats) and the soil depth
                                                                       !! (shallow, deep, average). Used in the calculation of the critical wind speed
                                                                       !! according to the GALES (Hale et al. 2015) model.
                                                                       !! In this case, the tree is in leaf.

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: overturning_free_draining_shallow_leafless_mtc  =  &  !! Overturning moment multiplier @tex $(Nm/kg)$ @endtex.
  & (/  undef,   175.3,  175.3,   134.7,   198.5,   175.3,  132.6,   & !! This is derived from the generic soil type (free_draining mineral soils;
  &      152.0,  145.2,  undef,   undef,   undef,   undef /)           !! Gleyed mineral soils; Peaty mineral soils; Deep peats) and the soil depth
                                                                       !! (shallow, deep, average). Used in the calculation of the critical wind speed
                                                                       !! according to the GALES (Hale et al. 2015) model.
                                                                       !! In this case, the tree is leafless.

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: overturning_free_draining_deep_mtc  =  &  !! Overturning moment multiplier @tex $(Nm/kg)$ @endtex.
  & (/  undef,   203.8,   203.8,   157.2,   230.8,   230.8,   154.8, & !! This is derived from the generic soil type (free_draining mineral soils;
  &     176.7,   169.4,  undef,   undef,   undef,   undef  /)          !! Gleyed mineral soils; Peaty mineral soils; Deep peats) and the soil depth
                                                                       !! (shallow, deep, average). Used in the calculation of the critical wind speed
                                                                       !! according to the GALES (Hale et al. 2015) model.
                                                                       !! In this case, the tree is in leaf.

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: overturning_free_draining_deep_leafless_mtc  =  &  !! Overturning moment multiplier @tex $(Nm/kg)$ @endtex.
  & (/  undef,   203.8,   203.8,   157.2,   230.8,   230.8,   154.8, & !! This is derived from the generic soil type (free_draining mineral soils;
  &      176.7,  169.4,   undef,  undef,   undef,  undef  /)           !! Gleyed mineral soils; Peaty mineral soils; Deep peats) and the soil depth
                                                                       !! (shallow, deep, average). Used in the calculation of the critical wind speed
                                                                       !! according to the GALES (Hale et al. 2015) model.
                                                                       !! In this case, the tree is leafless.

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: overturning_free_draining_average_mtc  =  &  !! Overturning moment multiplier @tex $(Nm/kg)$ @endtex.
  & (/  undef,   178.7,   178.7,   137.8,   202.4,   178.7,   135.7, & !! This is derived from the generic soil type (free_draining mineral soils;
  &      155.0,  148.6,   undef,   undef,  undef,   undef  /)          !! Gleyed mineral soils; Peaty mineral soils; Deep peats) and the soil depth
                                                                       !! (shallow, deep, average). Used in the calculation of the critical wind speed
                                                                       !! according to the GALES (Hale et al. 2015) model.
                                                                       !! In this case, the tree is in leaf.

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: overturning_free_draining_average_leafless_mtc  =  &  !! Overturning moment multiplier @tex $(Nm/kg)$ @endtex.
  & (/  undef,   178.7,   178.7,   137.8,   202.4,   178.7,   135.7, & !! This is derived from the generic soil type (free_draining mineral soils;
  &      155.0,  148.6,   undef,   undef,   undef,   undef  /)         !! Gleyed mineral soils; Peaty mineral soils; Deep peats) and the soil depth
                                                                       !! (shallow, deep, average). Used in the calculation of the critical wind speed
                                                                       !! according to the GALES (Hale et al. 2015) model.
                                                                       !! In this case, the tree is leafless.

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: overturning_gleyed_shallow_mtc  =  &  !! Overturning moment multiplier @tex $(Nm/kg)$ @endtex.
  & (/  undef,   155.4,   155.4,   119.4,   176.0,   155.4,   117.6, & !! This is derived from the generic soil type (free_draining mineral soils;
  &      134.8,  128.7,   undef,   undef,   undef,   undef  /)         !! Gleyed mineral soils; Peaty mineral soils; Deep peats) and the soil depth
                                                                       !! (shallow, deep, average). Used in the calculation of the critical wind speed
                                                                       !! according to the GALES (Hale et al. 2015) model.
                                                                       !! In this case, the tree is in leaf.

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: overturning_gleyed_shallow_leafless_mtc  =  &  !! Overturning moment multiplier @tex $(Nm/kg)$ @endtex.
  & (/  undef,   155.4,   155.4,   119.4,   176.0,   155.4,   117.6, & !! This is derived from the generic soil type (free_draining mineral soils;
  &      134.8,  128.7,   undef,   undef,   undef,   undef  /)         !! Gleyed mineral soils; Peaty mineral soils; Deep peats) and the soil depth
                                                                       !! (shallow, deep, average). Used in the calculation of the critical wind speed
                                                                       !! according to the GALES (Hale et al. 2015) model.
                                                                       !! In this case, the tree is leafless.

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: overturning_gleyed_deep_mtc  =  &  !! Overturning moment multiplier @tex $(Nm/kg)$ @endtex.
  & (/  undef,   180.6,   180.6,   139.3,   204.6,   180.6,   137.2, & !! This is derived from the generic soil type (free_draining mineral soils;
  &     156.7,   150.2,   undef,   undef,   undef,   undef  /)         !! Gleyed mineral soils; Peaty mineral soils; Deep peats) and the soil depth
                                                                       !! (shallow, deep, average). Used in the calculation of the critical wind speed
                                                                       !! according to the GALES (Hale et al. 2015) model.
                                                                       !! In this case, the tree is in leaf.

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: overturning_gleyed_deep_leafless_mtc  =  &  !! Overturning moment multiplier @tex $(Nm/kg)$ @endtex.
  & (/  undef,   180.6,   180.6,   139.3,   204.6,   180.6,   137.2, & !! This is derived from the generic soil type (free_draining mineral soils;
  &     156.7,   150.2,   undef,   undef,   undef,   undef  /)         !! Gleyed mineral soils; Peaty mineral soils; Deep peats) and the soil depth
                                                                       !! (shallow, deep, average). Used in the calculation of the critical wind speed
                                                                       !! according to the GALES (Hale et al. 2015) model.
                                                                       !! In this case, the tree is leafless.

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: overturning_gleyed_average_mtc  =  &  !! Overturning moment multiplier @tex $(Nm/kg)$ @endtex.
  & (/  undef,   158.5,   158.5,   122.2,   179.5,   158.5,  120.3, &  !! This is derived from the generic soil type (free_draining mineral soils;
  &     137.4,   131.7,   undef,   undef,   undef,   undef  /)         !! Gleyed mineral soils; Peaty mineral soils; Deep peats) and the soil depth
                                                                       !! (shallow, deep, average). Used in the calculation of the critical wind speed
                                                                       !! according to the GALES (Hale et al. 2015) model.
                                                                       !! In this case, the tree is in leaf.

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: overturning_gleyed_average_leafless_mtc  =  &  !! Overturning moment multiplier @tex $(Nm/kg)$ @endtex.
  & (/  undef,   158.8,   158.5,   122.2,   179.5,   158.5,   120.3, & !! This is derived from the generic soil type (free_draining mineral soils;
  &     137.4,   131.7,   undef,   undef,   undef,   undef /)          !! Gleyed mineral soils; Peaty mineral soils; Deep peats) and the soil depth
                                                                       !! (shallow, deep, average). Used in the calculation of the critical wind speed
                                                                       !! according to the GALES (Hale et al. 2015) model.
                                                                       !! In this case, the tree is leafless.

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: overturning_peaty_shallow_mtc  =  &  !! Overturning moment multiplier @tex $(Nm/kg)$ @endtex.
  & (/  undef,   169.7,   169.7,   130.4,   192.2,   169.7,   128.4, & !! This is derived from the generic soil type (free_draining mineral soils;
  &      147.2,  140.6,   undef,   undef,   undef,   undef  /)         !! Gleyed mineral soils; Peaty mineral soils; Deep peats) and the soil depth
                                                                       !! (shallow, deep, average). Used in the calculation of the critical wind speed
                                                                       !! according to the GALES (Hale et al. 2015) model.
                                                                       !! In this case, the tree is in leaf.

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: overturning_peaty_shallow_leafless_mtc  =  &  !! Overturning moment multiplier @tex $(Nm/kg)$ @endtex.
  & (/  undef,   169.7,   169.7,   130.4,   192.2,   169.7,   128.4, & !! This is derived from the generic soil type (free_draining mineral soils;
  &      147.2,  140.6,   undef,   undef,   undef,   undef  /)         !! Gleyed mineral soils; Peaty mineral soils; Deep peats) and the soil depth
                                                                       !! (shallow, deep, average). Used in the calculation of the critical wind speed
                                                                       !! according to the GALES (Hale et al. 2015) model.
                                                                       !! In this case, the tree is leafless.

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: overturning_peaty_deep_mtc  =  &  !! Overturning moment multiplier @tex $(Nm/kg)$ @endtex.
  & (/  undef,   191.4,   191.4,   152.1,   223.5,   191.4,   141.9, & !! This is derived from the generic soil type (free_draining mineral soils;
  &      159.2,  164.0,   undef,   undef,   undef,   undef  /)         !! Gleyed mineral soils; Peaty mineral soils; Deep peats) and the soil depth
                                                                       !! (shallow, deep, average). Used in the calculation of the critical wind speed
                                                                       !! according to the GALES (Hale et al. 2015) model.
                                                                       !! In this case, the tree is in leaf.

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: overturning_peaty_deep_leafless_mtc  =  &  !! Overturning moment multiplier @tex $(Nm/kg)$ @endtex.
  & (/  undef,   191.4,   191.4,   152.1,   223.5,   191.4,   141.9, & !! This is derived from the generic soil type (free_draining mineral soils;
  &      159.2,  164.0,   undef,   undef,   undef,   undef  /)         !! Gleyed mineral soils; Peaty mineral soils; Deep peats) and the soil depth
                                                                       !! (shallow, deep, average). Used in the calculation of the critical wind speed
                                                                       !! according to the GALES (Hale et al. 2015) model.
                                                                       !! In this case, the tree is leafless.

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: overturning_peaty_average_mtc  =  &  !! Overturning moment multiplier @tex $(Nm/kg)$ @endtex.
  & (/  undef,   178.9,   178.9,   133.4,   195.9,   178.9,  131.4, &  !! This is derived from the generic soil type (free_draining mineral soils;
  &     162.0,   143.8,   undef,   undef,   undef,   undef  /)         !! Gleyed mineral soils; Peaty mineral soils; Deep peats) and the soil depth
                                                                       !! (shallow, deep, average). Used in the calculation of the critical wind speed
                                                                       !! according to the GALES (Hale et al. 2015) model.
                                                                       !! In this case, the tree is in leaf.

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: overturning_peaty_average_leafless_mtc  =  &  !! Overturning moment multiplier @tex $(Nm/kg)$ @endtex.
  & (/  undef,   178.9,   178.9,   133.4,   195.9,   178.9,  131.4, &  !! This is derived from the generic soil type (free_draining mineral soils;
  &     162.0,   143.8,   undef,   undef,   undef,   undef  /)         !! Gleyed mineral soils; Peaty mineral soils; Deep peats) and the soil depth
                                                                       !! (shallow, deep, average). Used in the calculation of the critical wind speed
                                                                       !! according to the GALES (Hale et al. 2015) model.
                                                                       !! In this case, the tree is leafless.

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: overturning_peat_shallow_mtc  =  &  !! Overturning moment multiplier @tex $(Nm/kg)$ @endtex.
  & (/  undef,   193.0,   193.0,   148.3,  218.6,  193.0,   146.0, &   !! This is derived from the generic soil type (free_draining mineral soils;
  &     167.4,   159.9,   undef,   undef,  undef,  undef  /)           !! Gleyed mineral soils; Peaty mineral soils; Deep peats) and the soil depth
                                                                       !! (shallow, deep, average). Used in the calculation of the critical wind speed
                                                                       !! according to the GALES (Hale et al. 2015) model.
                                                                       !! In this case, the tree is in leaf.

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: overturning_peat_shallow_leafless_mtc  =  &  !! Overturning moment multiplier @tex $(Nm/kg)$ @endtex.
  & (/  undef,   193.0,   193.0,   148.3,   218.6, 193.0,   146.0, &   !! This is derived from the generic soil type (free_draining mineral soils;
  &     167.4,   159.9,   undef,   undef,   undef, undef  /)           !! Gleyed mineral soils; Peaty mineral soils; Deep peats) and the soil depth
                                                                       !! (shallow, deep, average). Used in the calculation of the critical wind speed
                                                                       !! according to the GALES (Hale et al. 2015) model.
                                                                       !! In this case, the tree is leafless.

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: overturning_peat_deep_mtc  =  &  !! Overturning moment multiplier @tex $(Nm/kg)$ @endtex.
  & (/  undef,   224.4,  224.4,  173.1,   254.2,  224.4,  170.4, &     !! This is derived from the generic soil type (free_draining mineral soils;
  &     194.7,   186.6,  undef,  undef,   undef,  undef  /)            !! Gleyed mineral soils; Peaty mineral soils; Deep peats) and the soil depth
                                                                       !! (shallow, deep, average). Used in the calculation of the critical wind speed
                                                                       !! according to the GALES (Hale et al. 2015) model.
                                                                       !! In this case, the tree is in leaf.

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: overturning_peat_deep_leafless_mtc  =  &  !! Overturning moment multiplier @tex $(Nm/kg)$ @endtex.
  & (/  undef,  224.4,   224.4,  173.1,  254.2,  224.4,   170.4, &     !! This is derived from the generic soil type (free_draining mineral soils;
  &     194.7,  186.6,   undef,  undef,  undef,  undef  /)             !! Gleyed mineral soils; Peaty mineral soils; Deep peats) and the soil depth
                                                                       !! (shallow, deep, average). Used in the calculation of the critical wind speed
                                                                       !! according to the GALES (Hale et al. 2015) model.
                                                                       !! In this case, the tree is leafless.

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: overturning_peat_average_mtc  =  &  !! Overturning moment multiplier @tex $(Nm/kg)$ @endtex.
  & (/  undef,   196.9,   196.9,   151.8,   223.0,   196.9,   149.4, & !! This is derived from the generic soil type (free_draining mineral soils;
  &     170.8,   163.6,   undef,   undef,   undef,   undef  /)         !! Gleyed mineral soils; Peaty mineral soils; Deep peats) and the soil depth
                                                                       !! (shallow, deep, average). Used in the calculation of the critical wind speed
                                                                       !! according to the GALES (Hale et al. 2015) model.
                                                                       !! In this case, the tree is in leaf.

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: overturning_peat_average_leafless_mtc  =  &  !! Overturning moment multiplier @tex $(Nm/kg)$ @endtex.
  & (/  undef,  196.9,   196.9,   151.8,   223.0,   196.9,    149.4, & !! This is derived from the generic soil type (free_draining mineral soils;
  &     170.8,  163.6,   undef,   undef,   undef,   undef  /)          !! Gleyed mineral soils; Peaty mineral soils; Deep peats) and the soil depth
                                                                       !! (shallow, deep, average). Used in the calculation of the critical wind speed
                                                                       !! according to the GALES (Hale et al. 2015) model.
                                                                       !! In this case, the tree is leafless.

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: max_damage_further_mtc = & !! A tunning parameter for determining wind damage rate/level @text $(unitless)$ @endtex.  
  & (/    undef,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8, &           !! The value of this tunning parameter is suggested by filed observation data for
  &       0.8,   0.8,  undef,   undef,   undef, undef /)               !! different PFTs.
  
  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: max_damage_closer_mtc = &  !! A tunning parameter for determining wind damage rate/level @text $(unitless)$ @endtex. 
  & (/    undef,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8, &           !! The value of this tunning parameter is suggested by filed observation data for
  &       0.8,   0.8,  undef,   undef,  undef, undef /)                !! different PFTs.

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: sfactor_further_mtc  =  &  !! A tunning parameter for determining wind damage rate/level @text $(unitless)$ @endtex. 
  & (/    undef,   5.0,   5.0,   5.0,   5.0,   5.0,   5.0, &           !! The value of this tunning parameter is suggested by filed observation data for
  &       5.0,   5.0,  undef,   undef,  undef,  undef /)               !! different PFTs.

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: sfactor_closer_mtc  =  &   !! A tunning parameter for determining wind damage rate/level @text $(unitless)$ @endtex. 
  & (/    undef,   5.0,   5.0,   5.0,   5.0,   5.0,   5.0, &           !! The value of this tunning parameter is suggested by filed observation data for
  &       5.0,   5.0,  undef,  undef,  undef,  undef /)                !! different PFTs.

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: green_density_mtc  =  &    !! The green wood density of the meta classes sepcies @text $(kg/m^3)$ @endtex.
  & (/    undef,  1007.,   1007.,   985.,  1060., 1007., 990.,    &    !!  
  &       968.,    900.,   undef,  undef,  undef, undef /)             !! 

  
  ! Spitfire
  !
  ! Note: the sources of these parameters are Thonicke et al. (2010) BG except
  ! for
  ! those otherwise indicated.
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: bulkdensity_deadfuel_mtc = & !!Fuel bulk density for litter fuel (kg dry mass m^{-3})
  & (/ 0.0,   25.0,   25.0,   25.0,   10.0,   22.0,   25.0,  &            !!Source: Thonicke et al. (2010) BG Table 1
  &   22.0,   22.0,    2.0,    2.0,    2.0,    2.0 /)
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: bulkdensity_livefuel_mtc = & !!Bulk density for live fuel (kg dry mass m^{-3}) 
  & (/ undef,  undef,  undef,  undef,  undef,  undef,  undef,  &          !! 
  &    undef,  undef,    4.0,    4.0,  undef,  undef /)
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: f_sh_mtc = &                 !!Scorch height parameter 
  & (/  0.0,  0.1487,  0.061,  0.1,  0.371,  0.094,  0.11,  &
  &   0.094,   0.094,    0.0,  0.0,    0.0,    0.0 /)
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: BTpar1_mtc = &               !!Bark thickness parameter
  & (/   0.0,  0.0301,  0.1085,  0.0670,  0.0451,  0.0347,  0.0292,  &
  &   0.0347,  0.0347,     0.0,     0.0,     0.0,     0.0 /)
  REAL(r_std) , PARAMETER, DIMENSION(nvmc) :: BTpar2_mtc = &              !!Bark thickness parameter
  & (/   0.0,  0.0281,  0.212,  0.559,  0.1412,  0.1086,  0.2632,  &
  &   0.1086,  0.1086,    0.0,    0.0,     0.0,     0.0 /)
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: r_ck_mtc = &                 !!Crown damage parameter   
  & (/ 0.0,    0.95,    0.05,    0.95,    0.95,    0.95,    0.95,  &
  &   0.95,    0.95,     0.0,     0.0,     0.0,     0.0 /)
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: p_ck_mtc = &                 !!Crown damage parameter
  & (/ 0.0,    3.0,    3.0,    3.75,    3.0,    3.0,    3.0,  &
  &    3.0,    3.0,    0.0,     0.0,    0.0,    0.0 /)
  REAL(r_std) , PARAMETER, DIMENSION(nvmc) :: me_mtc = &                  !!Moisture of extinction. The source is unclear.
  & (/ 1.0,    0.2,    0.3,    0.3,    0.3,    0.3,    0.35,  &
  &   0.35,   0.35,    0.2,    0.2,    0.2,    0.2 /)
  REAL(r_std) , PARAMETER, DIMENSION(nvmc) :: fdi_sf_pft_mtc = &          !!PFT-dependent scaling factor for fire danger index 
  & (/ 0.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0,  &
  &   1.0,   1.0,    1.0,    1.0,    0.0,    0.0 /)

  ! Emissions factors that convert dry mass emissions into different species (g (kg dry mass)^{-1})
  REAL(r_std) , PARAMETER, DIMENSION(nvmc) :: ef_CO2_mtc = &
  & (/   0.0,  1580.0,  1664.0,  1568.0,  1568.0,  1568.0,  1568.0,  &
  &   1568.0,  1568.0,  1568.0,  1664.0,  1568.0,  1664.0 /)
  REAL(r_std) , PARAMETER, DIMENSION(nvmc) :: ef_CO_mtc = &
  & (/  0.0,    103.0,    63.0,    106.0,    106.0,    106.0,    106.0,  &
  &   106.0,    106.0,   106.0,     63.0,    106.0,     63.0 /)
  REAL(r_std) , PARAMETER, DIMENSION(nvmc) :: ef_CH4_mtc = &
  & (/ 0.0,    6.8,    2.2,    4.8,    4.8,    4.8,    4.8,  &
  &    4.8,    4.8,    4.8,    2.2,    4.8,    2.2 /)
  REAL(r_std) , PARAMETER, DIMENSION(nvmc) :: ef_VOC_mtc = &
  & (/ 0.0,    8.1,    3.4,    5.7,    5.7,    5.7,    5.7,  &
  &    5.7,    5.7,    5.7,    3.4,    5.7,    3.4 /)
  REAL(r_std) , PARAMETER, DIMENSION(nvmc) :: ef_TPM_mtc = &
  & (/ 0.0,    8.5,    8.5,    17.6,    17.6,    17.6,    17.6,  &
  &   17.6,   17.6,   17.6,     8.5,    17.6,     8.5 /)
  REAL(r_std) , PARAMETER, DIMENSION(nvmc) :: ef_NOx_mtc = &
  & (/ 0.0,    1.999,    2.54,    3.24,    3.24,    3.24,    3.24,  &
  &   3.24,     3.24,    3.24,    2.54,    3.24,    2.54 /)


  !! The 100hr and 1000hr maximum combustion fractions (CF) are from Table 2 of Leeuwen et al. (2014),
  !! with 100hr corresponding to woody debris of 2.55-7.6 cm and 1000hr corresponding
  !! to woody debris of 7.2-20.5 cm. The CF for tropical tress are taken from
  !! Table 2 (b); The CF for temperate trees are taken from Table 2 (c). As the
  !! The equivalent CF for woody debris were not reported for boreal forests in 
  !! Table 2 (d), we tentatively use those for temperate forests. The CF for
  !! grassland are taken from Table 3 by averaging the CFs for savanna and grassland savanna 
  !! as we don't make distinction between savanna and grassland savanna in our
  !! model yet. The CFs for crops are set tentatively as zero to indicate the fact
  !! that cropland fires are not simulated in the model.
  REAL(r_std) , PARAMETER, DIMENSION(nvmc) :: fire_max_cf_100hr_mtc = &
  & (/ 0.0,    0.65,    0.65,    0.73,    0.73,    0.73,    0.73,  &
  &   0.73,    0.73,    0.76,    0.76,     0.0,     0.0 /)
  REAL(r_std) , PARAMETER, DIMENSION(nvmc) :: fire_max_cf_1000hr_mtc = &
  & (/ 0.0,    0.41,    0.41,    0.38,    0.38,    0.38,    0.38,  &
  &   0.38,    0.38,    0.76,    0.76,     0.0,     0.0 /)
  LOGICAL, PARAMETER, DIMENSION(nvmc) :: burnable_mtc =  &                       !! Is this PFT allowed to be burned? 
  & (/ .FALSE.,   .TRUE.,   .TRUE.,   .TRUE.,   .TRUE.,    .TRUE.,  .TRUE.,  &   !! (Boolean)
  &     .TRUE.,   .TRUE.,   .TRUE.,   .TRUE.,  .FALSE.,   .FALSE. /)




  !
  ! FLUX - LUC
  !
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: coeff_lcchange_1_mtc  =  &   !! Coeff of biomass export for the year
  & (/  undef,   0.897,   0.897,   0.597,   0.597,   0.597,   0.597,  &   !! (0-1, unitless)
  &     0.597,   0.597,   1.000,   1.000,   1.000,   1.000  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: coeff_lcchange_10_mtc  =  &  !! Coeff of biomass export for the decade 
  & (/  undef,   0.103,   0.103,   0.299,   0.299,   0.299,   0.299,  &   !! (0-1, unitless)
  &     0.299,   0.299,   0.000,   0.000,   0.000,   0.000  /) 

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: coeff_lcchange_100_mtc  =  & !! Coeff of biomass export for the century
  & (/  undef,     0.0,     0.0,   0.104,   0.104,   0.104,   0.104,  &   !! (0-1, unitless)
  &     0.104,   0.104,     0.0,   0.000,   0.000,   0.000  /)

  !
  ! FLUX - LUC
  !
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: coeff_lcchange_s_mtc  =  &   !! Coeff of biomass export for the year
  & (/  undef,   0.897,   0.897,   0.597,   0.597,   0.597,   0.597,  &   !! (0-1, unitless)
  &     0.597,   0.597,   1.000,   1.000,   1.000,   1.000  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: coeff_lcchange_m_mtc  =  &   !! Coeff of biomass export for the decade
  & (/  undef,   0.103,   0.103,   0.299,   0.299,   0.299,   0.299,  &   !! (0-1, unitless)
  &     0.299,   0.299,   0.000,   0.000,   0.000,   0.000  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: coeff_lcchange_l_mtc  =  &   !! Coeff of biomass export for the century
  & (/  undef,     0.0,     0.0,   0.104,   0.104,   0.104,   0.104,  &   !! (0-1, unitless)
  &     0.104,   0.104,     0.0,   0.000,   0.000,   0.000  /)



  !
  ! PHENOLOGY
  !
  ! The latest modifications regarding senescence_temp_c, leaffall, hum_min_time and nosenescence_hum are inspired by
  ! MacBean et al. (2015), following the optimization of phenology parameters using MODIS NDVI (FM/PP).
  !-
  ! 1. Stomate
  !-
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: lai_max_to_happy_mtc  =  &  !! threshold of LAI below which plant uses carbohydrate reserves
  & (/  undef,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,  &
  &       0.5,   0.5,   0.5,   0.5,   0.5,   0.5  /)

  REAL(r_std), PARAMETER, DIMENSION (nvmc) :: lai_max_mtc  =  &          !! maximum LAI, PFT-specific 
  & (/ undef,   7.0,   5.0,   5.0,   4.0,   5.0,   3.5,  &               !! @tex $(m^2.m^{-2})$ @endtex
  &      4.0,   3.0,   2.5,   2.0,   5.0,   5.0  /)

  INTEGER(i_std), PARAMETER, DIMENSION(nvmc) :: pheno_type_mtc  =  &     !! type of phenology (0-4, unitless)
  & (/  0,   1,   3,   1,   1,   2,   1,  &                              !! 0=bare ground 1=evergreen,  2=summergreen, 
  &     2,   2,   4,   4,   2,   3  /)                                   !! 3=raingreen,  4=perennial
  !-
  ! 2. Leaf Onset
  !-
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: force_pheno_mtc = &           !! Number of days after the mean 
  & (/ undef,  undef,  42.,  undef,  undef,  42.,  undef, &                  !! doy at which budbreak occurs 
  &    28.,  28.,  35.,  35.,  28.,  28. /)                                       !! at which phenology will be forced
  
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: pheno_gdd_crit_c_mtc  =  &    !! critical gdd, tabulated (C),
  & (/  undef,   undef,   undef,   undef,   undef,   undef,   undef,  &    !! constant c of aT^2+bT+c
  &     undef,   undef,   640.0,   400.0,   320.0,   300.0  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: pheno_gdd_crit_b_mtc  =  &    !! critical gdd, tabulated (C),
  & (/  undef,   undef,   undef,   undef,   undef,   undef,   undef,  &    !! constant b of aT^2+bT+c
  &     undef,   undef,    32.,     0.0,    6.25,     6.0  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: pheno_gdd_crit_a_mtc  =  &    !! critical gdd, tabulated (C),
  & (/  undef,   undef,     undef,   undef,   undef,   undef,   undef,  &  !! constant a of aT^2+bT+c
  &     undef,   undef,    0.1,   0.0,   0.0315,   0.01  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: pheno_moigdd_t_crit_mtc  = &  !! temperature threshold for C4 grass(C)
  & (/  undef,   undef,     undef,   undef,   undef,   undef,   undef,  &  
  &     undef,   undef,     undef,    22.0,   undef,   undef  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: ngd_crit_mtc  =  &            !! critical ngd, tabulated. 
  & (/  undef,   undef,   undef,   undef,   undef,   undef,   undef,  &    !! Threshold -5 degrees (days)
  &     undef,    23.9,   undef,   undef,   undef,   undef  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: ncdgdd_temp_mtc  =  &         !! critical temperature for the ncd vs. gdd 
  & (/  undef,   undef,   undef,   undef,   undef,     5.0,   undef,  &    !! function in phenology (C)
  &       0.0,   -10.0,   undef,   undef,   undef,   undef  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: hum_frac_mtc  =  &            !! critical humidity (relative to min/max) 
  & (/  undef,   undef,   0.5,   undef,   undef,   undef,   undef, &       !! for phenology (unitless)
  &     undef,   undef,   0.5,     0.5,     0.5,     0.5  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: hum_min_time_mtc  =  &        !! minimum time elapsed since
  & (/  undef,   undef,   50.0,   undef,   undef,   undef,   undef,  &     !! moisture minimum (days)
  &     undef,   undef,   36.0,    35.0,    75.0,    75.0  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: relsoilmoist_always_mtc  =  & !! moisture monthly availability above which moisture tendency
  & (/  undef,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,  &    !! doesn't matter (0-1, unitless)    
  &       1.0,     1.0,    0.98,     0.6,     0.6,     0.6  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: t_always_mtc  = &             !! temperature monthly above which temperature tendency (K)
  & (/  undef,  283.15,  283.15,  283.15,  283.15,  283.15,  283.15,  &    !! doesn't matter
  &    283.15,  283.15,    280.,  283.15,  283.15,  283.15  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: harvest_time_a_mtc  =  &      !! Parameter to calculate harvest time of crops (days)
  & (/  undef,   undef,   undef,   undef,   undef,   undef,   undef,  &     
  &     undef,   undef,   undef,   undef,   160.0,   160.0  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: harvest_time_b_mtc  =  &      !! Parameter to calculate harvest time of crops (K)
  & (/  undef,   undef,   undef,   undef,   undef,   undef,   undef,  &     
  &     undef,   undef,   undef,   undef,   282.0,   282.0  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: longevity_root_mtc  =  &      !! roots longevity (days). This parameter describes the root turnover
  & (/  undef,  200.,   200.,    300.,    300.,   300.,   300.,  &         !! within the growing season. A longevity of 1000 days implies that every
  &      300.,  300.,   100.,    100.,    100.,   100.  /)                 !! day 1/1000 of the root mass will be replaced. For a temperate PFT this
                                                                           !! implies that that 18% of the root mass dies (and thus needs to be
                                                                           !! replaced) within the growing season. This definition is straightforward 
                                                                           !! for evergreen species but needs the above nuance for deciduous PFTs 

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: longevity_sap_mtc  =  &       !! time (days). Because the sapwood is always there, the definition of this  
  & (/  undef, 7300., 3766., 3934., 3534., 1563., 4016., &                 !! parameter is straightforward for evergeens and deciduous species. A longevity
  &     1937., 1570., 180.,  180.,  480.,  480. /)                         !! of 7300 days means that 5% of the sapwood will turnover every year.

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: longevity_leaf_mtc  =  &      !! leaf longevity (days). This parameter describes the leaf
  & (/  undef,   730.,   180.,    910.,   730.,   180.,   2000.,  &       !! longevity at a location with the average temperature for the range of the PFT
  &      180.,   180.,   180.,    180.,   200.,   200.  /)               !! It is used to calculate leaf_age_crit (as a function of temperature.
                                                                           !! The variable leaf_age_crit describes the maximum leaf age. If the 
                                                                           !! actual leaf age exceeds leaf_age_crit, senescence will start.
                                                                           !! This definition is straightforward for evergreen species but needs some 
                                                                           !! nuance for deciduous PFTs because it is not clear when within season turnover 
                                                                           !! (describing stochastic processes of leaf mortality due to insects, wind, and 
                                                                           !! self-shading) are overruled by climatological senescence.

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: leaf_age_crit_tref_mtc = &    !! Reference temperature of the PFT (degrees Celsius)
  & (/ undef,   25.,   25.,   15.,   20.,   15.,   5.,    &                !! Used to calculate the leaf_age_crit as a function of
  &       5.,    5.,   15.,   20.,   15.,   20.    /)                      !! longevity_leaf

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: leaf_age_crit_coeff1_mtc = &  !! Coeff1 (unitless) to link leaf_age_crit to leaf_age_crit_tref
  & (/ undef,   1.,   1.,   1.,   1.,   1.,   2.0,    &                
  &       1.,   1.,   1.5,   1.5,   1.5,   1.5  /)      
  
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: leaf_age_crit_coeff2_mtc = &  !! Coeff2 (unitless) to link leaf_age_crit to leaf_age_crit_tref
  & (/ undef,   1.,   1.,   1.,   1.,   1.,   0.365,    &                
  &       1.,   1.,   0.75,   0.75,   0.75,   0.75 /) 

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: leaf_age_crit_coeff3_mtc = &  !! Coeff3 (unitless) to link leaf_age_crit to leaf_age_crit_tref
  & (/ undef,   0.,   0.,   0.,   0.,   0.,   120.,    &                
  &       0.,   0.,   10.,   10.,   10.,   10.  /) 

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: longevity_fruit_mtc  =  &     !! fruit lifetime (days)
  & (/  undef,  90.0,    90.0,    90.0,    90.0,   90.0,   90.0,  &
  &      90.0,  90.0,    undef,   undef,   undef,  undef  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: ecureuil_mtc  =  &            !! fraction of primary leaf and root allocation
  & (/  undef,   0.0,   1.0,   0.0,   0.0,   1.0,   0.0,  &                !! put into reserve (0-1, unitless)
  &       1.0,   1.0,   1.0,   1.0,   1.0,   1.0  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: alloc_min_mtc  =  &           !! allocation above/below = f(age) 
  & (/    0.2,   0.2,   0.2,   0.2,   0.2,  0.2,   0.35,  &
  &       0.7,   0.2,   0.2,   0.2,   0.2,  0.2  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: alloc_max_mtc  =  &           !! allocation above/below = f(age) 
  & (/    0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,  &
  &       0.8,   0.8,   0.8,   0.8,   0.8,   0.8  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: demi_alloc_mtc  =  &          !! allocation above/below = f(age) 
  & (/  undef,   5.0,     5.0,     5.0,     5.0,    5.0,   5.0,  &         
  &       5.0,   5.0,   undef,   undef,   undef,   undef  /)

  !-
  ! 3. Senescence
  !-
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: leaffall_mtc  =  &             !! length of death of leaves (days)
  & (/  undef,   undef,   5.0,    undef,    undef,   15.0,   undef,  &
  &      25.0,    25.0,    5.0,     5.0,    undef,   undef  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: presenescence_ratio_mtc  =  &  !! The ratio of maintenance respiration to
  & (/  undef,   0.6,     0.6,     0.6,     0.6,     0.6,     0.6,  &       !! gpp beyond which presenescence stage of
  &       0.6,   0.6,     0.6,     0.6,     0.6,     0.6  /)                !! plant phenology is declared to begin (0-1, unitless)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: senescence_ratio_mtc  =  &     !! The ratio of maintenance respiration to
  & (/  undef,   1.0,     1.0,     1.0,     1.0,     1.0,     1.0,  &       !! gpp beyond which senescence stage of
  &      0.87,  0.89,    0.75,     1.0,     1.0,     1.0  /)                !! plant phenology is declared to begin (0-1, unitless)

  CHARACTER(LEN=6), PARAMETER, DIMENSION(nvmc) :: senescence_type_mtc  =  & !! type of senescence, tabulated (unitless)
  & (/  'none  ',  'none  ',   'dry   ',  'none  ',  'none  ',  &
  &     'cold  ',  'none  ',   'cold  ',  'cold  ',  'mixed ',  &
  &     'mixed ',  'crop  ',   'crop  '            /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: senescence_hum_mtc  =  &       !! critical relative moisture availability
  & (/  undef,   undef,   0.3,   undef,   undef,   undef,   undef,  &       !! for senescence (0-1, unitless)
  &     undef,   undef,   0.2,     0.2,     0.3,     0.2  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: nosenescence_hum_mtc  =  &     !! relative moisture availability above which 
  & (/  undef,   undef,   0.8,   undef,   undef,   undef,   undef,  &       !! there is no humidity-related senescence
  &     undef,   undef,   0.3,     0.3,     0.3,     0.3  /)                !! (0-1, unitless)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: max_turnover_time_mtc  =  &    !! maximum turnover time for grasses (days)
  & (/  undef,   undef,   undef,   undef,   undef,   undef,   undef,  &
  &     undef,   undef,    80.0,    80.0,    80.0,    80.0  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: min_turnover_time_mtc  =  &    !! minimum turnover time for grasses (days)
  & (/  undef,   undef,   undef,   undef,   undef,   undef,   undef,  &
  &     undef,   undef,    10.0,    10.0,    10.0,    10.0  /)
 
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: recycle_leaf_mtc = &           !! Fraction of N leaf that is recycled when leaves are senescent
  & (/  undef,     0.5,     0.5,     0.5,     0.5,     0.5,     0.5,  &
  &       0.5,     0.5,     0.5,     0.5,     0.5,     0.5  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: recycle_root_mtc = &           !! Fraction of N leaf that is recycled when leaves are senescent
  & (/  undef,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2,  &
  &       0.2,     0.2,     0.2,     0.2,     0.2,     0.2  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: min_leaf_age_for_senescence_mtc  =  &  !! minimum leaf age to allow 
  & (/  undef,   undef,   90.0,   undef,   undef,   90.0,   undef,  &               !! senescence g (days)
  &      60.0,    60.0,   30.0,    30.0,    30.0,   30.0  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: senescence_temp_c_mtc  =  &    !! critical temperature for senescence (C)
  & (/  undef,   undef,    undef,   undef,   undef,   16.0,   undef,  &     !! constant c of aT^2+bT+c, tabulated
  &      14.0,    13.0,      5.0,     5.0,     12.0,    13.0  /)             !! (unitless)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: senescence_temp_b_mtc  =  &    !! critical temperature for senescence (C), 
  & (/  undef,   undef,   undef,   undef,   undef,   0.0,   undef,  &       !! constant b of aT^2+bT+c, tabulated
  &       0.0,     0.0,     0.1,     0.0,     0.0,   0.0  /)                !! (unitless)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: senescence_temp_a_mtc  =  &    !! critical temperature for senescence (C), 
  & (/  undef,   undef,     undef,   undef,   undef,   0.0,   undef,  &     !! constant a of aT^2+bT+c, tabulated
  &       0.0,     0.0,   0.00375,     0.0,     0.0,   0.0  /)              !! (unitless)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: gdd_senescence_mtc  =  &       !! minimum gdd to allow senescence of crops (days)
  & (/  undef,   undef,    undef,   undef,     undef,    undef,    undef,  &
  &     undef,   undef,    undef,   undef,     2500.,    2500.  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: pre_gdd_mtc  =  &              !! Subtracted from gdd_senescence to allow 
  & (/  undef,   undef,    undef,   undef,     undef,    undef,    undef, & !! presenescence in crops    
  &     undef,   undef,    undef,   undef,     1000.,    1000.  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: pre_harv_mtc  =  &             !! Subtracted from harvest_time to allow presenescence 
  & (/  undef,   undef,    undef,   undef,     undef,    undef,    undef, & !! presenescence in crops    
  &     undef,   undef,    undef,   undef,     500.,    500.  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: gr_gt_mtc  =  &                !! The minimum number of days after growth is 
  & (/  undef,   undef,    undef,   undef,     undef,    undef,    undef, & !! initialised that plants can enter presenescence 
  &     undef,   undef,    undef,   undef,     60.,    60.  /)              !! (currently only for crops) 

  LOGICAL, PARAMETER, DIMENSION(nvmc) :: always_init_mtc  =  &              !! take carbon from atmosphere if carbohydrate reserve too small (true/false)
  & (/ .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE., &!! default is true for all pfts except pft=11 C4 grass
  &    .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE. /)    

  !-
  ! 4. N cycle
  !-
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: cn_leaf_min_mtc  = &            !! minimum CN ratio of leaves 
  & (/  undef,     16.,      16.,     28.2,       16.,      16.,      28.2, &  !! (gC/gN)
  &       16.,     16.,      16.,     16.,       16.,      16.   /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: cn_leaf_max_mtc  = &            !! maximum CN ratio of leaves 
 & (/   undef,     45.5,      45.5,     74.8,       45.5,      45.5,      74.8,  & !! (gC/gN) 
 &        45.5,     45.5,      45.5,     45.5,       45.5,      45.5   /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: max_soil_n_bnf_mtc   = &        !! Value of total N (NH4+NO3) 
 & (/     0.0,     1.5,      1.5,     1.5,       1.5,      1.5,      1.5,  & 
 &        1.5,     1.5,       2.,      2.,        2.,       2.   /)          !! above which we stop adding N via BNF
                                                                             !! (gN/m**2)
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: manure_pftweight_mtc = &        !! Weight of the distribution of manure over the PFT surface
 & (/   0.,     0.,      0.,     0.,       0.,      0.,      0.,  &          !!(to a same number correspond the same concentration) 
 &        0.,     0.,      1.,     1.,       1.,      1.   /)

  !
  ! DGVM
  !
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: residence_time_mtc  =  &        !! residence time of trees (years). Functions as a back-up. Trees should
  & (/  undef,   1000.0,   1000.0,   1000.0,   1000.0,   1000.0,   1000.0, & !! die from other causes well before this value is reached.
  &     1000.0,  1000.0,   undef,    undef,    undef,    undef  /) 

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: tmin_crit_mtc  =  &
  & (/  undef,     0.0,     0.0,   -30.0,   -14.0,   -30.0,   -45.0,  &  !! critical tmin, tabulated (C)
  &     -45.0,   -60.0,   undef,   undef,   undef,   undef  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: tcm_crit_mtc  =  &
  & (/  undef,   undef,   undef,     5.0,    15.5,    15.5,   -8.0,  &   !! critical tcm, tabulated (C)
  &      -8.0,    -8.0,   undef,   undef,   undef,   undef  /)

REAL(r_std), PARAMETER, DIMENSION(nvmc) :: mortality_min_mtc = &       !! Asymptotic mortality if plant growth exceeds long term
  & (/  undef,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  &               !! NPP (thus a strongly growing PFT)
  &      0.01,  0.01,  0.01,  0.01,  0.01,  0.01  /)                      !! @tex $(year^{-1})$ @endtex
  
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: mortality_max_mtc = &       !! Maximum mortality if plants hardly grows thus 
  & (/  undef,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1,  &                      !! NPP << NPPlongterm @tex $(year^{-1})$ @endtex
  &       0.1,  0.1,  0.1,  0.1,  0.1,  0.1  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: ref_mortality_mtc = &       !! Reference mortality rate used to calculate mortality
  & (/  undef,  0.035,  0.035,  0.035,  0.035,  0.035,  0.035,  &        !! as a function of the plant vigor 
  &     0.035,  0.035,  0.035,  0.035,  0.035,  0.035  /)                !! @tex $(year^{-1})$ @endtex

  !
  ! Biogenic Volatile Organic Compounds
  !
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: em_factor_isoprene_mtc = &     !! Isoprene emission factor 
  & (/  0.,    24.,   24.,    8.,   16.,   45.,   8.,  &                    !!
  &    18.,    0.5,   12.,   18.,    5.,    5.  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: em_factor_monoterpene_mtc = &  !! Monoterpene emission factor
  & (/   0.,   2.0,    2.0,   1.8,    1.4,    1.6,    1.8,  &               !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
  &    1.4,    1.8,    0.8,   0.8,    0.22,     0.22  /)

  REAL(r_std), PARAMETER :: LDF_mono_mtc = 0.6                                  !! monoterpenes fraction dependancy to light
  REAL(r_std), PARAMETER :: LDF_sesq_mtc = 0.5                                  !! sesquiterpenes fraction dependancy to light
  REAL(r_std), PARAMETER :: LDF_meth_mtc = 0.8                                  !! methanol fraction dependancy to light
  REAL(r_std), PARAMETER :: LDF_acet_mtc = 0.2                                  !! acetone fraction dependancy to light

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: em_factor_apinene_mtc = &      !! Alfa pinene emission factor percentage
  & (/   0.,   0.395,   0.395,   0.354,   0.463,   0.326,   0.354, &        !! ATTENTION: for each PFT they are PERCENTAGE of monoterpene EF
  &   0.316,   0.662,   0.231,   0.200,   0.277,   0.277 /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: em_factor_bpinene_mtc = &      !! Beta pinene emission factor  percentage      
  & (/   0.,   0.110,   0.110,   0.146,   0.122,   0.087,   0.146, &        !! ATTENTION: for each PFT they are PERCENTAGE of monoterpene EF
  &   0.063,   0.150,   0.123,   0.080,   0.154,   0.154  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: em_factor_limonene_mtc = &     !! Limonene emission factor percentage
  & (/   0.,   0.092,   0.092,   0.083,   0.122,   0.061,   0.083, &        !! ATTENTION: for each PFT they are PERCENTAGE of monoterpene EF
  &   0.071,   0.037,   0.146,   0.280,   0.092,   0.092  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: em_factor_myrcene_mtc = &      !! Myrcene emission factor percentage
  & (/   0.,   0.073,   0.073,   0.050,   0.054,   0.028,   0.050, &        !! ATTENTION: for each PFT they are PERCENTAGE of monoterpene EF
  &   0.019,   0.025,   0.062,   0.057,   0.046,   0.046  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: em_factor_sabinene_mtc = &     !! Sabinene emission factor percentage
  & (/   0.,   0.073,   0.073,   0.050,   0.083,   0.304,   0.050, &        !! ATTENTION: for each PFT they are PERCENTAGE of monoterpene EF
  &   0.263,   0.030,   0.065,   0.050,   0.062,   0.062  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: em_factor_camphene_mtc = &     !! Camphene emission factor percentage
  & (/   0.,   0.055,   0.055,   0.042,   0.049,   0.004,   0.042, &        !! ATTENTION: for each PFT they are PERCENTAGE of monoterpene EF
  &   0.005,   0.023,   0.054,   0.053,   0.031,   0.031  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: em_factor_3carene_mtc = &      !! 3-carene emission factor percentage
  & (/   0.,   0.048,   0.048,   0.175,   0.010,   0.024,   0.175, &        !! ATTENTION: for each PFT they are PERCENTAGE of monoterpene EF
  &   0.013,   0.042,   0.065,   0.057,   0.200,   0.200  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: em_factor_tbocimene_mtc = &    !! T-beta-ocimene emission factor percentage
  & (/   0.,   0.092,   0.092,   0.054,   0.044,   0.113,   0.054, &        !! ATTENTION: for each PFT they are PERCENTAGE of monoterpene EF
  &   0.105,   0.028,   0.138,   0.120,   0.031,   0.031  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: em_factor_othermonot_mtc = &   !! Other monoterpenes emission factor percentage
  & (/   0.,   0.062,   0.062,   0.046,   0.054,   0.052,   0.046, &        !! ATTENTION: for each PFT they are PERCENTAGE of monoterpene EF
  &   0.144,   0.003,   0.115,   0.103,   0.108,   0.108  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: em_factor_sesquiterp_mtc = &   !! Sesquiterpene emission factor
  & (/   0.,  0.45,   0.45,   0.13,   0.30,   0.36,   0.15, &               !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
  &    0.30,  0.25,   0.60,   0.60,   0.08,   0.08  /)

  REAL(r_std), PARAMETER :: beta_mono_mtc = 0.10                            !! Monoterpenes temperature dependency coefficient 
  REAL(r_std), PARAMETER :: beta_sesq_mtc = 0.17                            !! Sesquiterpenes temperature dependency coefficient 
  REAL(r_std), PARAMETER :: beta_meth_mtc = 0.08                            !! Methanol temperature dependency coefficient 
  REAL(r_std), PARAMETER :: beta_acet_mtc = 0.10                            !! Acetone temperature dependency coefficient 
  REAL(r_std), PARAMETER :: beta_oxyVOC_mtc = 0.13                          !! Other oxygenated BVOC temperature dependency coefficient 


  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: em_factor_ORVOC_mtc = &        !! ORVOC emissions factor 
  &  (/  0.,    1.5,    1.5,    1.5,    1.5,   1.5,    1.5,  &              !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
  &     1.5,    1.5,    1.5,    1.5,    1.5,   1.5  /) 

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: em_factor_OVOC_mtc = &         !! OVOC emissions factor 
  &  (/  0.,    1.5,    1.5,    1.5,    1.5,   1.5,    1.5,  &              !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
  &     1.5,    1.5,    1.5,    1.5,    1.5,   1.5  /)
  
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: em_factor_MBO_mtc = &          !! MBO emissions factor
  & (/     0., 2.e-5, 2.e-5,   1.4, 2.e-5, 2.e-5, 0.14,  &                  !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
  &     2.e-5, 2.e-5, 2.e-5, 2.e-5, 2.e-5, 2.e-5  /)  
  
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: em_factor_methanol_mtc = &     !! Methanol emissions factor 
  & (/  0.,    0.8,   0.8,   1.8,   0.9,   1.9,   1.8,  &                   !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
  &    1.8,    1.8,   0.7,   0.9,    2.,     2.  /)  
  
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: em_factor_acetone_mtc = &      !! Acetone emissions factor
  & (/  0.,   0.25,   0.25,   0.30,   0.20,   0.33,   0.30,  &              !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
  &   0.25,   0.25,   0.20,   0.20,   0.08,   0.08  /)
  
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: em_factor_acetal_mtc = &       !! Acetaldehyde emissions factor 
  & (/  0.,   0.2,    0.2,     0.2,   0.2,   0.25,   0.25,   0.16,   &      !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
  &   0.16,   0.12,   0.12,   0.035,   0.020  /)  
  
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: em_factor_formal_mtc = &       !! Formaldehyde emissions factor
  & (/  0.,   0.04,   0.04,  0.08,    0.04,    0.04,  0.04,  &              !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
  &   0.04,   0.04,  0.025, 0.025,   0.013,   0.013  /)  

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: em_factor_acetic_mtc = &       !! Acetic Acid emissions factor
  & (/   0.,   0.025,   0.025,   0.025,   0.022,   0.08,   0.025,   &      !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
  &   0.022,   0.013,   0.012,   0.012,   0.008,   0.008  /)  

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: em_factor_formic_mtc = &       !! Formic Acid emissions factor
  & (/  0.,  0.015,  0.015,   0.02,    0.02,   0.025,  0.025,  &            !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
  &  0.015,  0.015,  0.010,  0.010,   0.008,   0.008  /)  

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: em_factor_DMS_mtc = &          !! DMS emissions factor
  & (/  0.,   0.037,   0.037,   0.027,   0.027,   0.027,   0.027,  &        !! @tex $(\mu g.g^{-1}.h^{-1})$ @endtex
  &  0.027,   0.027,   0.027,   0.027,   0.027,   0.027  /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: em_factor_H2S_mtc = &          !! H2S emissions factor
  & (/  0.,   0.0095,   0.0095,   0.0064,   0.0064,   0.0064,   0.0064,  &  !! @tex $(\mu g.g^{-1}.h^{-1})$ @endtex
  &  0.0064,   0.0064,   0.0064,   0.0064,   0.0064,   0.0064  /)

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: em_factor_no_wet_mtc = &        !! NOx emissions factor soil emissions and exponential
  & (/  0.,   2.6,   0.06,   0.03,   0.03,   0.03,   0.03,  &               !! dependancy factor for wet soils
  &  0.03,   0.03,   0.36,   0.36,   0.36,   0.36  /)                       !! @tex $(ngN.m^{-2}.s^{-1})$ @endtex

  REAL(r_std),PARAMETER, DIMENSION(nvmc) :: em_factor_no_dry_mtc = &        !! NOx emissions factor soil emissions and exponential
  & (/  0.,   8.60,   0.40,   0.22,   0.22,   0.22,   0.22,  &              !! dependancy factor for dry soils
  &   0.22,   0.22,   2.65,   2.65,   2.65,   2.65  /)                      !! @tex $(ngN.m^{-2}.s^{-1})$ @endtex

  !
  ! MORTALITY (stomate)
  !

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: death_distribution_factor_mtc = &    !!  The scale factor between the smallest and largest
       (/ undef,  1.,  1.,    1.,    1.,    1.,    1.,  &                         !! circ class for tree mortality in stomate_mark_kill.
          1.,     1.,  undef, undef, undef, undef /)                              !! (unitless)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: npp_reset_value_mtc = &  !! The value of the NPP that the long-term value is
       (/ undef,  undef,  undef,  undef,  undef,  undef,  undef,  &   !! reset to after a PFT dies in stomate_kill.  This
          undef,  undef,   500.,   500.,   500.,   500. /)            !! only seems to be used for non-trees.
                                                                      !! @tex $(gC m^{-2})$ @endtex 

  !
  ! ALLOCATION and related
  !

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: Larch_mtc = &                  !! Larcher 1991 SAI/LAI ratio (unitless)
  & (/   0.,   0.015,   0.015,   0.003,   0.005,   0.005,   0.003,  &
  &   0.005,   0.003,   0.005,   0.005,   0.008,   0.008  /)  

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: k_latosa_max_mtc = &  !! Maximum leaf-to-sapwood area ratio as defined in McDowell et al
  & (/ undef,  13600., 10000., 4300.,  5800.,  20000., 6700., &    !! 2002, Oecologia and compiled in Hickler et al 2006, Appendix S2 
  &    19800., 30049., 7000.,   5500.,  9500.,  11000. /)          !! The values for grasses and crops are tuned. More work is needed
                                                                   !! to fully justify this approach for the herbacuous PFTs (unitless)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: k_latosa_min_mtc = &  !! Minimum leaf-to-sapwood area ratio as defined in McDowell et al
  & (/ undef,   7300.,  5500.,  3300.,  2800.,  7300.,  3700., &  !! 2002, Oecologia and compiled in Hickler et al 2006, Appendix S2 
  &    13000., 14730.,  7000.,  5500.,  9500.,  11000. /)          !! The values for grasses and crops are tuned. More work is needed
                                                                   !! to fully justify this approach for the herbacuous PFTs (unitless)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: lai_to_height_mtc = &                    !! Convertion from lai to height for grasses
  &(/ undef, undef, undef, undef, undef, undef, undef, &                              !! and cropland. Convert lai because that way a dynamic
  &   undef, undef,   0.1,   0.1,   0.1,   0.1 /)                                     !! sla is accounted for  

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: deleuze_a_mtc = &                          !! intercept of the intra-tree competition within a stand
  & (/ undef,  0.23,  0.23,  0.23,  0.23,  0.23,  0.23, &                               !! based on the competion rule of Deleuze and Dhote 2004
  &     0.23,  0.23, undef, undef, undef, undef /)                                      !! Used when n_circ > 6

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: deleuze_b_mtc = &                          !! slope of the intra-tree competition within a stand
  & (/ undef,  0.58,  0.58,  0.58,  0.58,  0.58,  0.58, &                               !! based on the competion rule of Deleuze and Dhote 2004
  &     0.58,  0.58, undef, undef, undef, undef /)                                      !! Used when n_circ > 6

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: deleuze_p_all_mtc = &                      !! Percentile of the circumferences that receives photosynthates
  & (/ undef,  0.50,  0.50,  0.70,  0.70,  0.70,  0.70, &                               !! based on the competion rule of Deleuze and Dhote 2004
  &     0.70,  0.70, undef, undef, undef, undef /)                                      !! Used when n_circ > 6 for FM1, FM2 and FM4

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: deleuze_p_coppice_mtc = &                  !! Percentile of the circumferences that receives photosynthates
  & (/ undef,  0.50,  0.50,  0.50,  0.50,  0.50,  0.50, &                               !! based on the competion rule of Deleuze and Dhote 2004
  &     0.50,  0.50, undef, undef, undef, undef /)                                      !! Used when n_circ > 6 for FM3

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: deleuze_power_a_mtc = &                    !! slope to calculate divisor of the power for the slope
  & (/ undef,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, &                                     !! of intra-tree competition whithin a stand
  &     0.0,  0.0, undef, undef, undef, undef /)                                        !! based on the competition rule of Deleuze and Dhote 2004

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: m_dv_mtc = &                               !! Parameter in the Deleuze & Dhote allocation 
  & (/ undef,  1.05,  1.05,  1.05,  1.05,  1.05,  1.05, &                               !! rule that relaxes the cut-off imposed by
  &     1.05,  1.05,    0.,    0.,    0.,    0. /)                                      !! ::sigma. Owing to m_relax trees still grow
                                                                                        !! a little when their ::circ is below ::sigma
  
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: fruit_alloc_mtc = &   !! Fraction of biomass allocated to fruit production (0-1)
  & (/ undef,   0.1,    0.1,    0.1,     0.1,     0.1,    0.1, &   !! currently only parameterized for forest PFTs 
  &      0.1,   0.1,    0.1,    0.1,     0.2,     0.2 /)  

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: frac_growthresp_res_lim_mtc = &          !! Fraction of growth respiration expressed as 
  &(/  0.28,   0.28,   0.28,   0.28,   0.28,   0.28,   0.28, &                        !! share of the total C that is to be allocated
  &    0.28,   0.28,   0.28,   0.28,   0.28,   0.28 /) 

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: frac_growthresp_fun_all_mtc = &          !! Fraction of growth respiration expressed as
  &(/  0.28,   0.28,   0.28,   0.28,   0.28,   0.28,   0.28, &                        !! share of the total C that is to be allocated
  &    0.28,   0.28,   0.28,   0.28,   0.28,   0.28 /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: labile_reserve_mtc = &                   !! Parameter to calculate labile reserve pool 
  &(/  0.,   2.0,  2.0,  2.0, 2.0,  2.0,  2.0, &                                        
  &    2.0,  2.0,  1.0,   1.0,  4.0,   4.0  /)                                         

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: evergreen_reserve_mtc = &                !! Fraction of sapwood mass stored in the reserve pool of evergreen 
  &(/ undef, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, &                                    !! trees (unitless, 0-1)
  &    0.05,  0.05, 0.05, 0.05, 0.05, 0.05 /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: senescense_reserve_mtc = &               !! Fraction of sapwood mass stored in the reserve pool of deciduous 
  &(/ undef, 0.15, 0.15, 0.15, 0.25, 0.25, 0.15, &                                    !! trees during senescense(unitless, 0-1)
  &    0.25,  0.15, 0.15, 0.15, 0.15, 0.15 /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: deciduous_reserve_mtc = &                !! Fraction of sapwood mass stored in the reserve pool of deciduous 
  &(/ undef, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, &                                    !! trees during the growing season (unitless, 0-1) 
  &    0.24,  0.24, 0.3, 0.3, 0.3, 0.3 /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: root_reserve_mtc = &                     !! Fraction of max root biomass which are 
                                                                                
  &(/ undef, 0.3, 1., 0.3, 0.3, 1., 0.3, &                                            !! covered by carbon reserve 
  &    1.,  1., 1., 1., 1., 1. /)                                                     !! for deciduous species we keep the whole root mass. 
                                                                                      !! For evergreens we are happy with 30%
 
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: fcn_root_mtc = &      !! N/C of "root" for allocation relative to leaf N/C  according 
  & (/ undef,   .86,    .86,    .86,    .86,    .86,   .86,   &    !! to Sitch et al 2003 (https://doi.org/10.1046/j.1365-2486.2003.00569.x)
  &      .86,   .86,    .86,    .86,    .86,    .86   /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: fcn_wood_mtc = &      !! N/C of "wood" for allocation relative to leaf N/C  according 
  & (/ undef,    .087,   .087,   .087,   .087,   .087,  .087,  &   !! to Sitch et al 2003 (https://doi.org/10.1046/j.1365-2486.2003.00569.x)
  &     .087,    .087,   .087,   .087,   .087,   .087   /)

  !
  ! CROP MANAGEMENT
  !
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: harvest_ratio_mtc = &                      !! Share of biomass that is removed from the site during harvest
  & (/ undef,  undef,  undef,  undef,  undef,  undef,  undef, &                         !! A high value indicates a high harvest efficiency and thus a 
  &    undef,  undef,  0.5,  0.5,    0.5,    0.5 /)                                     !! low input of residuals. (unitless, 0-1). 

  !
  ! FOREST MANAGEMENT
  !
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: forest_managed_forced_mtc =  &             !! The user decided against using a PFT-map and will prescribe
  & (/ undef, 1., 1., 1., 1., 1., 1., 1., 1., undef, undef, undef, undef /)             !! forest management 

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: thinstrat_mtc =  &                         !! Thinning strategy. The FM code distinguished
  & (/   0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &                              !! thinning from above (<0) or from below (>0).  
  &    1.0,   1.0,     0.0,     0.0,     0.0,     0.0   /)                              

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: taumin_mtc =  &                            !! Minimum probability that a tree get thinned (unitless)
  & (/   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  &
  &    0.0,   0.0,     0.0,     0.0,     0.0,     0.0   /) 

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: taumax_mtc =  &                            !! Maximum probability that a tree get thinned (unitless)
  & (/   0.0,   0.30,   0.30,   0.30,   0.30,   0.30,   0.30,  &
  &    0.30,   0.30,     0.0,     0.0,     0.0,     0.0   /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: circ_dist_shape_none_mtc = &               !! Circumference class distribution shape
   &(/ undef,    1.2,    1.2,    1.2,  1.2,    1.2,    1.2, &                           !! parameter for unmanaged forests (Weibull distribution) 
   &     1.2,    1.2,  undef,  undef,  undef,  undef/) 

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: circ_dist_shape_thin_mtc = &               !! Circumference class distribution shape
   &(/ undef,    -0.09,    -0.09,    -0.09, -0.09,    -0.09,    -0.09, &                !! parameter for rotational forests (Weibull distribution)  
   &     -0.09,    -0.09,  undef,  undef,  undef,  undef/)
 
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: circ_dist_shape_uneven_mtc = &             !! Circumference class distribution shape                 
   &(/ undef,    1.2,    1.2,    1.2,  1.2,    1.2,    1.2, &                           !! parameter for uneven-aged forests (Weibull distribution) 
   &     1.2,    1.2,  undef,  undef,  undef,  undef/) 

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: circ_dist_shape_cop_mtc = &                !! Circumference class distribution shape           
   &(/ undef,    1.2,    1.2,    1.2,  1.2,    1.2,    1.2, &                           !! parameter for coppice forests (Weibull distribution) 
   &     1.2,    1.2,  undef,  undef,  undef,  undef/) 

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: circ_dist_shape_src_mtc = &                !! Circumference class distribution shape  
   &(/ undef,    1.2,    1.2,    1.2, 1.2,    1.2,    1.2, &                            !! for short rotation cycle forests (Weibull distribution) 
   &     1.2,    1.2,  undef,  undef,  undef,  undef/)  

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: final_dia_none_mtc =  &                    !! Maximal tree diameter (m) controlling RDI
  & (/ undef,  0.65,     0.65,    0.55,    0.50,    0.60,   0.42,  &                    !! for unmanaged forests 
  &     0.50,  0.50,   undef,  undef,  undef,  undef/)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: final_dia_thin_mtc =  &                    !! Maximal tree diameter (m) controlling RDI 
  & (/ undef,  0.45,   0.45,    0.35,    0.30,    0.38,   0.30,  &                      !! for rotational forests  
  &     0.20,  0.30,   undef,  undef,  undef,  undef/) 

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: final_dia_uneven_mtc =  &                  !! Maximal tree diameter (m) controlling RDI
  & (/ undef,  0.45,   0.45,    0.35,    0.30,    0.40,   0.30,  &                      !! for uneven-aged forests 
  &     0.30,  0.30,   undef,  undef,  undef,  undef/) 

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: final_dia_cop_mtc =  &                     !! Maximal tree diameter (m) controlling RDI
  & (/ undef,  0.35,     0.3,    0.25,    0.25,    0.3,   0.25,  &                      !! for coppice forests
  &      0.3,  0.25,   undef,  undef,  undef,  undef/)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: final_dia_src_mtc =  &                     !! Maximal tree diameter (m) controlling RDI
  & (/ undef,  0.20,     0.15,    0.20,    0.20,    0.20,   0.20,  &                    !! for short rotation cycle forests                     
  &     0.20,  0.20,   undef,  undef,  undef,  undef/)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: fuelwood_diameter_mtc = &                  !! Diameter below which the wood harvest is used as fuelwood (m)
  &(/ undef,   0.3,   0.3,   0.2,   0.2,   0.2,   0.2, &                                !! Affects the way the wood is used in the dim_product_use   
  &     0.1,   0.1, undef, undef, undef, undef/)                                        !! subroutine 

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: coppice_kill_be_wood_mtc = &               !! The fraction of the belowground wood killed during coppicing.
  &(/ undef,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, &                                !! (unitless)
  &     0.0,   0.0, undef, undef, undef, undef/) 

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: branch_ratio_mtc =  &                      !! Ratio of branches to total woody biomass (unitless)
  & (/  0.0,   0.38,   0.38,   0.25,   0.38,   0.38,   0.25,  &
  &    0.38,   0.25,    0.0,    0.0,    0.0,    0.0 /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: coppice_diameter_mtc = &                   !! The trunk diameter above which one coppices
  & (/ undef,  0.15,  0.15,  0.15,  0.15,  0.15,  0.15, &                               !! trees. (m)
  &     0.15,  0.15,    0.,    0.,    0.,    0. /)                                      

  INTEGER(i_std), PARAMETER, DIMENSION(nvmc) :: shoots_per_stool_mtc = &                !! The number of shoots which regrow on a stool after
  (/ 9999,  6,     6,     6,      6,     6,     6, &                                    !! coppicing. 
  &     6,  6,  9999,  9999,   9999,  9999 /)   

  INTEGER(i_std), PARAMETER, DIMENSION(nvmc) :: src_rot_length_mtc = &                  !! The number of years between SRC cuttings.
  (/ 9999,   3,     3,     3,     3,     3,     3, &                                    !! (-)
  &     3,   3,  9999,  9999,  9999,  9999 /)                            

  INTEGER(i_std), PARAMETER, DIMENSION(nvmc) :: src_nrots_mtc = &                       !! The number of SRC rotations before the whole stand
  (/ 9999,   10,    10,    10,    10,    10,    10, &                                   !! is harvested (-)
  &    10,   10,  9999,  9999,  9999,  9999 /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: ndying_year_mtc = &                        !! Number of years during which an unmanaged forest will die after    
  (/ undef, 15.,   15.,   25.,   15.,   15.,   15., &                                   !! reaching the tree stem density threshold
  &    15., 15.,   15.,   15.,   15.,   15.  /)

  !! Forest management parameters for hectare-scale runs
  !  - Distinguishes between biomass and timber production for PFT4 and PFT6
  !  - Realistic thinning: infrequent but strong when it happens
  !  - The results should be representative for the stand but become unrealistic 
  !     when used at large pixels.
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: largest_tree_dia_ha_mtc =  &               !! Maximal tree diameter (m). If this diameter is exceeded a
  & (/ undef,  0.45,   0.45,    0.35,    0.30,    0.40,   0.30,  &                      !! a clearcut will happen
  &     0.30,  0.30,   undef,  undef,  undef,  undef/) 

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: thin_intensity_ha_mtc = &                  !! Intensity of a thinnig event when a forest reach the upper
  &(/ undef,   0.25,   0.25,  0.25,   0.25,    0.25,   0.25, &                          !! rdi limit
  &    0.25,   0.25,  undef,  undef,  undef,  undef/)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: rdi_dia_ref_ha_mtc = &                     !! RDI at a reference tree diameter used to shape the end of 
  &(/ undef,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8, &                                !! the RDI/diameter relationship
  &    0.8,   0.8,  undef,  undef,  undef,  undef/)   

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: rdi_start_ha_mtc = &                       !! The rdi value at the start of a forest rotation, use
  &(/ undef,    0.4, 0.4,  0.4,  0.4,  0.4,  0.4, &                                     !! 0.4 to obtain a reallistic growth rate during
  &   0.4,  0.4,  undef,     undef,    undef,    undef/)                                !! the early stage of the forest

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: rdi_max_none_ha_mtc = &                    !! The maximum RDI that a forest will reach during 
  &(/ undef,   0.9,   0.9,   0.9,   0.9,   0.9,   0.9, &                                !! its rotation : unmanaged forest
  &    0.9,   0.9,  undef,  undef,  undef,  undef/)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: rdi_max_thin_ha_mtc = &                    !! The maximum RDI that a forest will reach during 
  &(/ undef,    0.7,    0.7,    0.7,    0.7,    0.7,    0.7, &                          !! its rotation : rotational managed forest
  &     0.7,    0.7,  undef,  undef,  undef,  undef/)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: rdi_max_uneven_ha_mtc = &                  !! The maximum RDI that a forest will reach during 
  &(/ undef,  0.6,  0.6,  0.6,  0.6,  0.6,  0.6, &                                      !! its rotation : uneven-aged managed forest
  &   0.6,  0.6,  undef,  undef,  undef,  undef/)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: rdi_max_cop_ha_mtc = &                     !! The maximum RDI that a forest will reach during 
  &(/ undef,  1.2,  1.2,  1.2,  1.2,  1.2,  1.2, &                                      !! its rotation : coppice managed forest
  &  1.2, 1.2,   undef,     undef,    undef,   undef/)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: rdi_max_src_ha_mtc = &                     !! The maximum RDI that a forest will reach during 
  &(/ undef,  1.2,  1.2,  1.2,  1.2,  1.2,  1.2, &                                      !! its rotation : src managed forest
  &  1.2, 1.2,   undef,     undef,    undef,   undef/)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: dia_boost_end_ha_mtc = &                   !! Diameter threshold at which the rdi is
  &(/ undef,   10.0,   10.0,   10.0,   10.0,   10.0,   10.0, &                          !! not artificially maintained close to rdi_start
  &    10.0,   10.0,  undef,  undef,  undef,  undef/)

 REAL(r_std), PARAMETER, DIMENSION(nvmc) :: dens_target_none_ha_mtc = &                 !! Minimal density. Below this density the forest will
  &(/ undef,   200.0,   200.0,   200.0,   200.0,   200.0,   200.0, &                    !! be clearcut (trees.ha-1) : unmanaged forest
  &    200.0,   200.0,  undef,  undef,  undef,  undef/)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: dens_target_thin_ha_mtc = &                !! Minimal density. Below this density the forest will
  &(/ undef,   200.0,   200.0,   200.0,   200.0,   200.0,   200.0, &                    !! be clearcut (trees.ha-1) : rotational managed forest
  &    200.0,    200.0,  undef,  undef,  undef,  undef/)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: dens_target_uneven_ha_mtc = &              !! Minimal density. Below this density the forest will
  &(/ undef,   200.0,   200.0,   200.0,   200.0,   200.0,   200.0, &                    !! be clearcut (trees.ha-1) : uneven-aged managed forest
  &    200.0,   200.0,  undef,  undef,  undef,  undef/)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: dens_target_cop_ha_mtc = &                 !! Minimal density. Below this density the forest will
  &(/ undef,   200.0,   200.0,   200.0,   200.0,   200.0,   200.0, &                    !! be clearcut (trees.ha-1) : coppice managed forest
  &  200.0, 200.0,   undef,     undef,    undef,   undef/)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: dens_target_src_ha_mtc = &                 !! Minimal density. Below this density the forest will
  &(/ undef,  200.0,  200.0,  200.0,  200.0,  200.0,  200.0, &                          !! be clearcut (trees.ha-1) : src managed forest
  &  200.0, 200.0,   undef,     undef,    undef,   undef/)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: branch_harvest_none_ha_mtc =  &            !! Ratio of branches harvested in unmanaged management.
  & (/  0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &                               !! CHECK whether this parameter is already used for this FM in the code
  &     1.0,   1.0,   0.0,   0.0,   0.0,   0.0 /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: branch_harvest_thin_ha_mtc =  &            !! Ratio of branches harvested in FM2 management.
  & (/  0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &                               !! This parameter is used for FM2
  &     1.0,   1.0,   0.0,   0.0,   0.0,   0.0 /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: branch_harvest_uneven_ha_mtc =  &          !! Ratio of branches harvested in FM5 management.
  & (/  0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &                               !! This parameter is used for FM5
  &     1.0,   1.0,   0.0,   0.0,   0.0,   0.0 /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: branch_harvest_cop_ha_mtc =  &             !! Ratio of branches harvested in FM3 management.
  & (/  0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &                               !! CHECK whether this parameter is already used for this FM in the code
  &     1.0,   1.0,   0.0,   0.0,   0.0,   0.0 /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: branch_harvest_src_ha_mtc =  &             !! Ratio of branches harvested in FM4 management.
  & (/  0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &                               !! CHECK whether this parameter is already used for this FM in the code
  &     1.0,   1.0,   0.0,   0.0,   0.0,   0.0 /)
  
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: clearcut_residue_none_ha_mtc =  &          !! Ratio of wood harvest during a thinning/clearcut. the wood left is considered unselable or too expensive to remove.
  & (/  0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &                               !! CHECK whether this parameter is already used for this FM in the code
  &     1.0,   1.0,   0.0,   0.0,   0.0,   0.0 /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: clearcut_residue_thin_ha_mtc =  &          !! Ratio of wood harvest during a thinning/clearcut. the wood left is considered unselable or too expensive to remove.
  & (/  0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &                               !! This parameter is used for FM2
  &     1.0,   1.0,   0.0,   0.0,   0.0,   0.0 /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: clearcut_residue_uneven_ha_mtc =  &        !! Ratio of wood harvest during a thinning/clearcut. the wood left is considered unselable or too expensive to remove.
  & (/  0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &                               !! This parameter is used for FM5
  &     1.0,   1.0,   0.0,   0.0,   0.0,   0.0 /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: clearcut_residue_cop_ha_mtc =  &           !! Ratio of wood harvest during a thinning/clearcut. the wood left is considered unselable or too expensive to remove.
  & (/  0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &                               !! CHECK whether this parameter is already used for this FM in the code
  &     1.0,   1.0,   0.0,   0.0,   0.0,   0.0 /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: clearcut_residue_src_ha_mtc =  &           !! Ratio of wood harvest during a thinning/clearcut. the wood left is considered unselable or too expensive to remove.
  & (/  0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &                               !! CHECK whether this parameter is already used for this FM in the code
  &     1.0,   1.0,   0.0,   0.0,   0.0,   0.0 /)
  
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: salvage_dist_none_ha_mtc =  &              !! Ratio of wood remove after a mortality event the wood left is considered unselable or too expensive to remove.
  & (/  0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &                               !! CHECK whether this parameter is already used for this FM in the code. If used, it would represent
  &     1.0,   1.0,   0.0,   0.0,   0.0,   0.0 /)                                       !! sanitation cuttings rather than salvage logging. Different term but similar effect on the C and N cycles.

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: salvage_dist_thin_ha_mtc =  &              !! Ratio of wood remove after a mortality event the wood left is considered unselable or too expensive to remove.
  & (/  0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &
  &     1.0,   1.0,   0.0,   0.0,   0.0,   0.0 /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: salvage_dist_uneven_ha_mtc =  &            !! Ratio of wood remove after a mortality event the wood left is considered unselable or too expensive to remove.
  & (/  0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &
  &     1.0,   1.0,   0.0,   0.0,   0.0,   0.0 /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: salvage_dist_cop_ha_mtc =  &               !! Ratio of wood remove after a mortality event the wood left is considered unselable or too expensive to remove.
  & (/  0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &
  &     1.0,   1.0,   0.0,   0.0,   0.0,   0.0 /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: salvage_dist_src_ha_mtc =  &               !! Ratio of wood remove after a mortality event the wood left is considered unselable or too expensive to remove.
  & (/  0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &
  &     1.0,   1.0,   0.0,   0.0,   0.0,   0.0 /)

  !! Forest management parameters for global runs
  !  - All PFTs are managed for biomass production
  !  - Pixel-scale thinning: little wood is thinned but very frequently
  !  - When ORCHIDEE read the FM map and finds FM2 but if the scale is global this is set to FM5
  !  - No clear cuts so when integrating over space, FM2 should be simulated as FM5
  !  - The results should be representative for the entire pixel, not for the stand
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: largest_tree_dia_global_mtc =  &           !! Maximal tree diameter (m). If this diameter is exceeded a
  & (/ undef,  0.32,   0.23,    0.35,    0.30,    0.38,   0.24,  &                      !! a clearcut will happen
  &     0.20,  0.30,   undef,  undef,  undef,  undef/)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: thin_intensity_global_mtc = &              !! Intensity of a thinnig event when a forest reach the upper
  &(/ undef,   0.01,   0.01,  0.01,   0.01,    0.01,   0.01, &                          !! rdi limit
  &    0.01,   0.01,  undef,  undef,  undef,  undef/)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: rdi_dia_ref_global_mtc = &                 !! RDI at a reference tree diameter used to shape the end of  
  &(/ undef,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8, &                                !! the RDI/diameter relationship
  &    0.8,   0.8,  undef,  undef,  undef,  undef/) 

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: rdi_start_global_mtc = &                   !! The rdi value at the start of a forest rotation, use
  &(/ undef,    0.4, 0.4,  0.4,  0.4,  0.4,  0.4, &                                     !! 0.4 to obtain a reallistic growth rate during
  &   0.4,  0.4,  undef,     undef,    undef,    undef/)                                !! the early stage of the forest

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: rdi_max_none_global_mtc = &                !! The maximum RDI that a forest will reach during 
  &(/ undef,   0.9,   0.9,   0.9,   0.9,   0.9,   0.9, &                                !! its rotation : unmanaged forest
  &    0.9,    0.9, undef,  undef,  undef,  undef/)                                     !! its rotation : coppice managed forest

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: rdi_max_thin_global_mtc = &                !! The maximum RDI that a forest will reach during
  &(/ undef,    0.7,    0.7,    0.7,    0.7,    0.7,    0.7, &                          !! its rotation : evenaged managed forest
  &     0.7,    0.7,  undef,  undef,  undef,  undef/)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: rdi_max_uneven_global_mtc = &              !! The maximum RDI that a forest will reach during
  &(/ undef,  0.7,  0.7,  0.7,  0.7,  0.7,  0.7, &                                      !! its rotation : unevenaged managed forest
  &   0.7,  0.7,  undef,  undef,  undef,  undef/)   

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: rdi_max_cop_global_mtc = &                 !! The maximum RDI that a forest will reach during
  &(/ undef,  1.2,  1.2,  1.2,  1.2,  1.2,  1.2, &                                      !! its rotation : coppice managed forest
  &  1.2, 1.2,   undef,     undef,    undef,   undef/)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: rdi_max_src_global_mtc = &                 !! The maximum RDI that a forest will reach during
  &(/ undef,  1.2,  1.2,  1.2,  1.2,  1.2,  1.2, &                                      !! its rotation : src managed forest
  &  1.2, 1.2,   undef,     undef,    undef,   undef/)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: dia_boost_end_global_mtc = &               !! Diameter threshold at which the rdi is
  &(/ undef,   10.0,   10.0,   10.0,   10.0,   10.0,   10.0, &                          !! not artificially maintained close to rdi_start
  &    10.0,   10.0,  undef,  undef,  undef,  undef/)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: dens_target_none_global_mtc = &            !! Minimal density. Below this density the forest will
  &(/ undef,   200.0,   200.0,   200.0,   200.0,   200.0,   200.0, &                    !! be clearcut (trees.ha-1) : unmanaged forest
  &    200.0,   200.0,  undef,  undef,  undef,  undef/)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: dens_target_thin_global_mtc = &            !! Minimal density. Below this density the forest will
  &(/ undef,   200.0,   200.0,   200.0,   200.0,   200.0,   200.0, &                    !! be clearcut (trees.ha-1) : rotational managed forest
  &    200.0,   200.0,  undef,  undef,  undef,  undef/)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: dens_target_uneven_global_mtc = &          !! Minimal density. Below this density the forest will
  &(/ undef,   200.0,   200.0,   200.0,   200.0,   200.0,   200.0, &                    !! be clearcut (trees.ha-1) : uneven-aged managed forest
  &    200.0,   200.0,  undef,  undef,  undef,  undef/)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: dens_target_cop_global_mtc = &             !! Minimal density. Below this density the forest will
  &(/ undef,   200.0,   200.0,   200.0,   200.0,   200.0,   200.0, &                    !! be clearcut (trees.ha-1) : coppice managed forest
  &    200.0,  200.0,   undef,  undef,  undef,  undef/)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: dens_target_src_global_mtc = &             !! Minimal density. Below this density the forest will
  &(/ undef,   200.0,   200.0,   200.0,   200.0,   200.0,   200.0, &                    !! be clearcut (trees.ha-1) : src managed forest
  &    200.0,  200.0,   undef,  undef,  undef,  undef/)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: branch_harvest_none_global_mtc =  &        !! Ratio of branches harvested in unmanaged management.
  & (/  0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &                               !! CHECK whether this parameter is already used for this FM in the code
  &     1.0,   1.0,   0.0,   0.0,   0.0,   0.0 /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: branch_harvest_thin_global_mtc =  &        !! Ratio of branches harvested in FM2 management.
  & (/  0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &                               !! This parameter is used for FM2
  &     1.0,   1.0,   0.0,   0.0,   0.0,   0.0 /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: branch_harvest_uneven_global_mtc =  &      !! Ratio of branches harvested in FM5 management.
  & (/  0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &                               !! This parameter is used for FM5
  &     1.0,   1.0,   0.0,   0.0,   0.0,   0.0 /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: branch_harvest_cop_global_mtc =  &         !! Ratio of branches harvested in FM3 management.
  & (/  0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &                               !! CHECK whether this parameter is already used for this FM in the code
  &     1.0,   1.0,   0.0,   0.0,   0.0,   0.0 /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: branch_harvest_src_global_mtc =  &         !! Ratio of branches harvested in FM4 management.
  & (/  0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &                               !! CHECK whether this parameter is already used for this FM in the code
  &     1.0,   1.0,   0.0,   0.0,   0.0,   0.0 /)
  
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: clearcut_residue_none_global_mtc =  &      !! Ratio of wood harvest during a thinning/clearcut. the wood left is considered unselable or too expensive to remove.
  & (/  0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &                               !! CHECK whether this parameter is already used for this FM in the code
  &     1.0,   1.0,   0.0,   0.0,   0.0,   0.0 /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: clearcut_residue_thin_global_mtc =  &      !! Ratio of wood harvest during a thinning/clearcut. the wood left is considered unselable or too expensive to remove.
  & (/  0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &                               !! This parameter is used for FM2
  &     1.0,   1.0,   0.0,   0.0,   0.0,   0.0 /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: clearcut_residue_uneven_global_mtc =  &    !! Ratio of wood harvest during a thinning/clearcut. the wood left is considered unselable or too expensive to remove.
  & (/  0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &                               !! This parameter is used for FM5. In the global set-up it is used to implement selective logging in the tropics.
  &     1.0,   1.0,   1.0,   0.0,   0.0,   0.0 /)                                       !! In the temperate and boreal forest it is used to account for a harvest ratio < 1. In reality the sink will be in
                                                                                        !! be in the biomass (due to land abandonment) in ORCHIDEE the sink will be in the soil = first order approximation.

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: clearcut_residue_cop_global_mtc =  &       !! Ratio of wood harvest during a thinning/clearcut. the wood left is considered unselable or too expensive to remove.
  & (/  0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &                               !! CHECK whether this parameter is already used for this FM in the code
  &     1.0,   1.0,   0.0,   0.0,   0.0,   0.0 /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: clearcut_residue_src_global_mtc =  &       !! Ratio of wood harvest during a thinning/clearcut. the wood left is considered unselable or too expensive to remove.
  & (/  0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &                               !! CHECK whether this parameter is already used for this FM in the code
  &     1.0,   1.0,   0.0,   0.0,   0.0,   0.0 /)
  
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: salvage_dist_none_global_mtc =  &          !! Ratio of wood remove after a mortality event the wood left is considered unselable or too expensive to remove.
  & (/  0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &                               !! CHECK whether this parameter is already used for this FM in the code. If used, it would represent
  &     1.0,   1.0,   0.0,   0.0,   0.0,   0.0 /)                                       !! sanitation cuttings rather than salvage logging. Different term but similar effect on the C and N cycles.

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: salvage_dist_thin_global_mtc =  &          !! Ratio of wood remove after a mortality event the wood left is considered unselable or too expensive to remove.
  & (/  0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &
  &     1.0,   1.0,   0.0,   0.0,   0.0,   0.0 /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: salvage_dist_uneven_global_mtc =  &        !! Ratio of wood remove after a mortality event the wood left is considered unselable or too expensive to remove.
  & (/  0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &
  &     1.0,   1.0,   0.0,   0.0,   0.0,   0.0 /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: salvage_dist_cop_global_mtc =  &           !! Ratio of wood remove after a mortality event the wood left is considered unselable or too expensive to remove.
  & (/  0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &
  &     1.0,   1.0,   0.0,   0.0,   0.0,   0.0 /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: salvage_dist_src_global_mtc =  &           !! Ratio of wood remove after a mortality event the wood left is considered unselable or too expensive to remove.
  & (/  0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &
  &     1.0,   1.0,   0.0,   0.0,   0.0,   0.0 /)

  !! Forest management parameters for European runs
  !  - Distinguishes between biomass and timber production for PFT4 and PFT6
  !  - Use these settings with age classes
  !  - Probably not the best set-up to plot NBP maps because the results may 
  !     show too much heterogeneity due to forest management.
  !  - Realistic thinning: infrequent but strong when it happens
  ! MTC 7 REVISE 0.24 -> 0.32 m. Guillaume M. -- The boreal conifer was felled at 24 cm,
  ! well under what Nordic forestry actually harvests: the Swedish NFI (SC Sweden, pine +
  ! spruce) gives a dominant diameter of 31.2 cm for stands aged 60-80 yr, the real felling
  ! window, and 34.9 cm at 80-100 yr. At 24 cm a stand reached the mature class in ~38 yr
  ! against the ~50 yr observed. /!\ This raises the threshold; it does NOT fix the ~1.3x
  ! excess in simulated diameter growth (0.53 vs 0.41 cm/yr observed), a separate defect.
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: largest_tree_dia_eu_mtc =  &               !! Maximal tree diameter (m). If this diameter is exceeded a
  & (/ undef,  0.32,   0.23,    0.35,    0.30,    0.38,   0.32,  &                      !! a clearcut will happen
  &     0.32,  0.30,   undef,  undef,  undef,  undef/)

  ! ROTATION DE REFERENCE PAR MTC (2026-07-31). Duree de revolution d'un peuplement
  ! conduit en gestion MULTIFONCTIONNELLE : c'est la reference que MI_ROTATION_FACTOR
  ! module ensuite par classe d'intensite. Le taux de recolte en decoule,
  ! taux = 1/(rotation_ref * facteur), et REMPLACE l'estimateur EFDA -- une serie de
  ! 39 ans ne peut pas mesurer une rotation de 100-200 ans.
  ! Sources : itineraires techniques CNPF (France), PROF Lisboa e Vale do Tejo
  ! (ICNF 2018) pour l'eucalyptus iberique.
  !   MTC 4 conif. temp.  : douglas 50-70, pin maritime 35-60, pin sylvestre -> 60
  !     /!\ Cette valeur est TEMPEREE (sources francaises) et le MTC 4 court de 34,5 a
  !     68,5 N. La sylviculture nordique est autre : ~80 ans en Suede du Sud, minimum legal
  !     de 45 (station fertile, sud) a 100 ans (station pauvre, nord). Ne PAS l'absorber en
  !     allongeant cette valeur, ce qui casserait la France : c'est un probleme de GROUPE.
  !   MTC 5 sclerophylle  : chene vert / chene liege, >100 ans               -> 120
  !   MTC 6 feuillu decidu temp. : chene 100-150, hetre 70-120               -> 110
  !   MTC 7 conif. boreal : plus lent que le tempere                         -> 110
  !   MTC 8 feuillu decidu boreal (bouleau)                                  -> 80
  !   MTC 9 meleze : 40-60                                                   -> 50
  !   MTC 2-3 tropicaux : AUCUNE source, valeur temperee la plus proche par defaut.
  ! /!\ L'EUCALYPTUS n'a PAS de MTC propre (nvmc=13) : il emprunte un MTC existant et
  ! ne s'en distingue que par des surcharges de namelist. Sa rotation (10-14 ans en
  ! taillis de trituration, 25-35 ans en futaie de sciage) doit donc etre posee par
  ! ROTATION_REF sur ses creneaux PFT ; elle ne peut pas venir de cette table.
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: rotation_ref_mtc =  &                      !! Reference rotation length (yr) under
  & (/ undef,  110.,   80.,    60.,    120.,    110.,   110.,  &                        !! MULTIFUNCTIONAL management
  &     80.,   50.,   undef,  undef,  undef,  undef/)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: thin_intensity_eu_mtc = &                  !! Intensity of a thinnig event when a forest reach the upper
  &(/ undef,   0.25,   0.25,  0.25,   0.25,    0.25,   0.25, &                          !! rdi limit
  &    0.25,   0.25,  undef,  undef,  undef,  undef/)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: rdi_dia_ref_eu_mtc = &                     !! RDI at a reference tree diameter used to shape the end of
  &(/ undef,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8, &                                !! the RDI/diameter relationship
  &    0.8,   0.8,  undef,  undef,  undef,  undef/)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: rdi_start_eu_mtc = &                       !! The rdi value at the start of a forest rotation, use
  &(/ undef,    0.4, 0.4,  0.4,  0.4,  0.4,  0.4, &                                     !! 0.4 to obtain a reallistic growth rate during
  &   0.4,  0.4,  undef,     undef,    undef,    undef/)                                !! the early stage of the forest

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: rdi_max_none_eu_mtc = &                    !! The maximum RDI that a forest will reach during 
  &(/ undef,   0.9,   0.9,   0.9,   0.9,   0.9,   0.9, &                                !! its rotation : unmanaged forest
  &    0.9,   0.9,  undef,  undef,  undef,  undef/)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: rdi_max_thin_eu_mtc = &                    !! The maximum RDI that a forest will reach during 
  &(/ undef,    0.7,    0.7,    0.7,    0.7,    0.7,    0.7, &                          !! its rotation : rotational managed forest
  &     0.7,    0.7,  undef,  undef,  undef,  undef/)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: rdi_max_uneven_eu_mtc = &                  !! The maximum RDI that a forest will reach during 
  &(/ undef,  0.6,  0.6,  0.6,  0.6,  0.6,  0.6, &                                      !! its rotation : uneven-aged managed forest
  &   0.6,  0.6,  undef,  undef,  undef,  undef/)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: rdi_max_cop_eu_mtc = &                     !! The maximum RDI that a forest will reach during 
  &(/ undef,  1.2,  1.2,  1.2,  1.2,  1.2,  1.2, &                                      !! its rotation : coppice managed forest
  &  1.2, 1.2,   undef,     undef,    undef,   undef/)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: rdi_max_src_eu_mtc = &                     !! The maximum RDI that a forest will reach during
  &(/ undef,  1.2,  1.2,  1.2,  1.2,  1.2,  1.2, &                                      !! its rotation : src managed forest
  &  1.2, 1.2,   undef,     undef,    undef,   undef/)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: dia_boost_end_eu_mtc = &                   !! Diameter threshold at which the rdi is
  &(/ undef,   10.0,   10.0,   10.0,   10.0,   10.0,   10.0, &                          !! not artificially maintained close to rdi_start
  &    10.0,   10.0,  undef,  undef,  undef,  undef/)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: dens_target_none_eu_mtc = &                !! Minimal density. Below this density the forest will
  &(/ undef,   200.0,   200.0,   200.0,   200.0,   200.0,   200.0, &                    !! be clearcut (trees.ha-1) : unmanaged forest
  &    200.0,   200.0,  undef,  undef,  undef,  undef/)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: dens_target_thin_eu_mtc = &                !! Minimal density. Below this density the forest will
  &(/ undef,   200.0,   200.0,   200.0,   200.0,   200.0,   200.0, &                    !! be clearcut (trees.ha-1) : rotational managed forest
  &    200.0,   200.0,  undef,  undef,  undef,  undef/)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: dens_target_uneven_eu_mtc = &              !! Minimal density. Below this density the forest will
  &(/ undef,   200.0,   200.0,   200.0,   200.0,   200.0,   200.0, &                    !! be clearcut (trees.ha-1) : uneven-aged managed forest
  &    200.0,   200.0,  undef,  undef,  undef,  undef/)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: dens_target_cop_eu_mtc = &                 !! Minimal density. Below this density the forest will
  &(/ undef,   200.0,   200.0,   200.0,   200.0,   200.0,   200.0, &                    !! be clearcut (trees.ha-1) : coppice managed forest
  &    200.0,   200.0,  undef,  undef,  undef,  undef/)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: dens_target_src_eu_mtc = &                 !! Minimal density. Below this density the forest will
  &(/ undef,   200.0,   200.0,   200.0,   200.0,   200.0,   200.0, &                    !! be clearcut (trees.ha-1) : src managed forest
  &    200.0,   200.0,  undef,  undef,  undef,  undef/)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: branch_harvest_none_eu_mtc =  &            !! Ratio of branches harvested in unmanaged management.
  & (/  0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &                               !! CHECK whether this parameter is already used for this FM in the code
  &     1.0,   1.0,   0.0,   0.0,   0.0,   0.0 /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: branch_harvest_thin_eu_mtc =  &            !! Ratio of branches harvested in FM2 management.
  & (/  0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &                               !! This parameter is used for FM2
  &     1.0,   1.0,   0.0,   0.0,   0.0,   0.0 /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: branch_harvest_uneven_eu_mtc =  &          !! Ratio of branches harvested in FM5 management.
  & (/  0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &                               !! This parameter is used for FM5
  &     1.0,   1.0,   0.0,   0.0,   0.0,   0.0 /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: branch_harvest_cop_eu_mtc =  &             !! Ratio of branches harvested in FM3 management.
  & (/  0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &                               !! CHECK whether this parameter is already used for this FM in the code
  &     1.0,   1.0,   0.0,   0.0,   0.0,   0.0 /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: branch_harvest_src_eu_mtc =  &             !! Ratio of branches harvested in FM4 management.
  & (/  0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &                               !! CHECK whether this parameter is already used for this FM in the code
  &     1.0,   1.0,   0.0,   0.0,   0.0,   0.0 /)
  
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: clearcut_residue_none_eu_mtc =  &          !! Ratio of wood harvest during a thinning/clearcut. the wood left is considered unselable or too expensive to remove.
  & (/  0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &                               !! CHECK whether this parameter is already used for this FM in the code
  &     1.0,   1.0,   0.0,   0.0,   0.0,   0.0 /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: clearcut_residue_thin_eu_mtc =  &          !! Ratio of wood harvest during a thinning/clearcut. the wood left is considered unselable or too expensive to remove.
  & (/  0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &                               !! This parameter is used for FM2
  &     1.0,   1.0,   0.0,   0.0,   0.0,   0.0 /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: clearcut_residue_uneven_eu_mtc =  &        !! Ratio of wood harvest during a thinning/clearcut. the wood left is considered unselable or too expensive to remove.
  & (/  0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &                               !! This parameter is used for FM5
  &     1.0,   1.0,   0.0,   0.0,   0.0,   0.0 /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: clearcut_residue_cop_eu_mtc =  &           !! Ratio of wood harvest during a thinning/clearcut. the wood left is considered unselable or too expensive to remove.
  & (/  0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &                               !! CHECK whether this parameter is already used for this FM in the code
  &     1.0,   1.0,   0.0,   0.0,   0.0,   0.0 /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: clearcut_residue_src_eu_mtc =  &           !! Ratio of wood harvest during a thinning/clearcut. the wood left is considered unselable or too expensive to remove.
  & (/  0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &                               !! CHECK whether this parameter is already used for this FM in the code
  &     1.0,   1.0,   0.0,   0.0,   0.0,   0.0 /)
  
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: salvage_dist_none_eu_mtc =  &              !! Ratio of wood remove after a mortality event the wood left is considered unselable or too expensive to remove.
  & (/  0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &                               !! CHECK whether this parameter is already used for this FM in the code. If used, it would represent
  &     1.0,   1.0,   0.0,   0.0,   0.0,   0.0 /)                                       !! sanitation cuttings rather than salvage logging. Different term but similar effect on the C and N cycles.

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: salvage_dist_thin_eu_mtc =  &              !! Ratio of wood remove after a mortality event the wood left is considered unselable or too expensive to remove.
  & (/  0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &
  &     1.0,   1.0,   0.0,   0.0,   0.0,   0.0 /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: salvage_dist_uneven_eu_mtc =  &            !! Ratio of wood remove after a mortality event the wood left is considered unselable or too expensive to remove.
  & (/  0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &
  &     1.0,   1.0,   0.0,   0.0,   0.0,   0.0 /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: salvage_dist_cop_eu_mtc =  &               !! Ratio of wood remove after a mortality event the wood left is considered unselable or too expensive to remove.
  & (/  0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &
  &     1.0,   1.0,   0.0,   0.0,   0.0,   0.0 /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: salvage_dist_src_eu_mtc =  &               !! Ratio of wood remove after a mortality event the wood left is considered unselable or too expensive to remove.
  & (/  0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  &
  &     1.0,   1.0,   0.0,   0.0,   0.0,   0.0 /)

  !! +++ End of scale dependent forestry parameters


  !
  ! ALLOCATION
  !
  !+++CHECK+++
  ! Ideally k_root (see below) and a value for k_soil_to_root are used to calculate
  ! k_belowground which is used in the allocation. Problem is that we need the root biomass
  ! to calculate k_soil_to_root and that we need k_soil_to_root to calculate the root
  ! biomass. Now allocation starts from Cs, it could be written to start from Cr but that
  ! is not an easy task. The benefits would be full consistency between allocation
  ! and plant water stress and a dynamic root allocation (as has been observed). The
  ! latter has been tried in 2014 but it turned out to be much more difficult than 
  ! expected. This development needs to be thought throught more carefully. In addition
  ! ::k in hydrology is the effective soil conductivity. To use that value we would
  ! have reconsider the hydraulic architecture. The solution for the moment is to
  ! decouple allocation and hydraulic architecture by defining two related parameters
  ! independently.
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: k_belowground_mtc = &                       !! Belowground (roots + soil) specific conductivity.     
  & (/-9999.,     8.62E-08,  6.6E-08,  1.929E-08,  2.82E-08,  3.0E-07,  5.11E-09, &  !! @tex $(m^{3} kg^{-1} s^{-1} MPa^{-1})$ @endtex    
  &   3.287E-07,  4.85E-08,  4.E-07,   4.E-07,     4.E-07,    4.E-07 /)  

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: cn_leaf_init_mtc = &                       !! C/N of leaves according to Sitch et al 2003
  (/ undef,  25.,  25.,  41.7,  25.,  25.,  43.,  &                                     !! (https://doi.org/10.1046/j.1365-2486.2003.00569.x)
       25.,  25.,  25.,  25.,  25.,  25.   /)
  !
  ! HYDRAULIC ARCHITECTURE
  !

  !++++++++++

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: k_root_mtc = &                    !! Fine root conductivity. Values based on Bonan et al 2014      
  & (/ undef,  7.02E-4,  7.02E-4,  7.02E-4,  7.02E-4,  7.02E-4,  7.02E-4, &    !! et al. 2006. @tex $(m^{3} kg^{-1} s^{-1} MPa^{-1})$ @endtex    
  &       7.02E-4,  7.02E-4,   7.02E-4,  7.02E-4,  7.02E-4,  7.02E-4   /)         
                                                              
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: k_sap_mtc = &                     !! Maximal sapwood specific conductivity. Values compiled in T. Hickler 
  & (/-9999.,  0.009,  0.002,  0.0011,  0.0002,  0.0012,  0.0013, &            !! et al. 2006. @tex $(m^{2} s^{-1} MPa^{-1})$ @endtex 
  &   0.0012,  0.0012, 0.0006, 0.0006,  0.0006,  0.0006   /)                   !! Values from DOFOCO run.def

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: k_leaf_mtc = &        !! Leaf conductivity. Values compiled in T. Hickler et al 2006
  & (/ undef,  2.5,  2.5,  1.5,  1.5,  2.5,  1.5,         &        !! @tex $(m s^{-1} MPa^{-1})$ @endtex
         2.5,  2.5,  3.0,  3.0,  3.0,  3.0         /)*1.E-7

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: psi_leaf_min_mtc = &  !! Minimal leaf water potential. Values in T. Hickler et al 2006
  & (/ undef, -2.2, -2.2, -1.95, -4.48, -2.2, -1.78,   &           !! @tex $(MPa)$ @endtex
        -2.2, -3.0, -3.0, -2.2,  -3.0,  -2.2        /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: psi_soil_tune_mtc = & !! Additive tuning parameter to account for soil-root interactions
  & (/ undef,  0.,  0.,  0.,  0.,  0.,  0.,  &                     !! @tex $(MPa)$ @endtex
          0.,  0.,  0.,  0.,  0.,  0. /)                                  

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: srl_mtc = &           !! Specific root length (m g⁻1). Values are obtained from Metcalfe 
   & (/ undef,  18.3,   18.3,   18.3,  18.3,  18.3,  18.3,  &      !! et al. 2008 and Ostonen et al. 2007  
         18.3,  18.3,   18.3,   18.3,  18.3,  18.3 /) 

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: r_froot_mtc = &       !! Fine root radius (m). Values are obtained from Bonan  
   & (/ undef,  0.29E-3,   0.29E-3,  0.29E-3,  0.29E-3, &          !! et al.2014 and Ostonen et al. 2007 (Tree physiology)   
   &  0.21E-3,  0.24E-3,   0.21E-3,  0.21E-3,  0.075E-3, &
   &  0.075E-3, 0.168E-3,  0.168E-3 /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: psi_root_min_mtc = &      !! The minimum root water potential. Tested by Emilie joetzjer
   & (/ undef,  -5.,   -5.,   -5.,  -5.,   -5.,  -5.,  &           !! 
   &      -5.,  -5.,   -5.,   -5.,  -5.,   -5. /)                  

  !
  ! SOIL DECOMPOSITION
  !
  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: som_turn_ipassive_mtc = &  !! turnover of the passive pool (year-1)
   & (/ 0.0025, 0.0025, 0.0025, 0.0025, 0.0025, 0.0025, 0.0025,  &            
   &    0.0025, 0.0025, 0.0025, 0.0025, 0.0025, 0.0025 /)

  REAL(r_std), PARAMETER, DIMENSION(nvmc) :: hack_som_turn_ipassive_mtc = &  !! turnover of the passive pool (year-1)
   & (/ 0.0025, 0.0025, 0.0025, 0.0160, 0.0027, 0.0058, 0.0060,  &            
   &    0.0025, 0.0025, 0.0060, 0.0012, 0.0045, 0.0100 /)

  !
  ! Peatland
  ! 
  LOGICAL, PARAMETER, DIMENSION(nvmc) :: is_peat_mtc  =  & !!natural peatland vegetations
  & (/ .FALSE.,  .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE., &
  &    .FALSE.,  .FALSE.,   .FALSE.,   .FALSE.,    .FALSE.,   .FALSE. /)  

  LOGICAL, PARAMETER, DIMENSION(nvmc) :: is_croppeat_mtc  =  & !!crops on peatland
  & (/ .FALSE.,  .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE.,   .FALSE., &
  &    .FALSE.,  .FALSE.,   .FALSE.,   .FALSE.,    .FALSE.,  .FALSE. /)  



END MODULE constantes_mtc
