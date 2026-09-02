! =================================================================================================================================
! MODULE       : pft_parameters_var
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2011)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        This module contains the variables in function of plant funtional type (pft).
!!
!!\n DESCRIPTION: This module contains the declarations for the externalized variables in function of the 
!!                plant foncional type(pft). \n
!!                The module is already USE in module pft_parameters. Therefor no need to USE it seperatly except
!!                if the subroutines in module pft_parameters are not needed.\n
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S)	: None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_parameters/pft_parameters_var.f90 $
!! $Date: 2026-05-13 15:11:27 +0200 (mer. 13 mai 2026) $
!! $Revision: 9540 $
!! \n
!_ ================================================================================================================================

MODULE pft_parameters_var

  USE defprec
  
  IMPLICIT NONE


  !
  ! PFT GLOBAL
  !
  INTEGER(i_std), SAVE :: nvm = 13                                       !! Number of vegetation types (2-N, unitless)
!$OMP THREADPRIVATE(nvm)

  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:) :: pft_to_mtc          !! Table of conversion : we associate one pft to one metaclass 
                                                                         !! (1-13, unitless)
!$OMP THREADPRIVATE(pft_to_mtc)

  CHARACTER(LEN=35), ALLOCATABLE, SAVE, DIMENSION(:) :: PFT_name         !! Description of the PFT (unitless)
!$OMP THREADPRIVATE(PFT_name)

  LOGICAL, SAVE   :: l_first_pft_parameters = .TRUE.                     !! To keep first call trace of the module (true/false)
!$OMP THREADPRIVATE(l_first_pft_parameters)

  !
  ! VEGETATION STRUCTURE
  !
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:) :: leaf_tab            !! leaf type (1-4, unitless)
                                                                         !! 1=broad leaved tree, 2=needle leaved tree, 
                                                                         !! 3=grass 4=bare ground 
!$OMP THREADPRIVATE(leaf_tab)

  CHARACTER(len=6), ALLOCATABLE, SAVE, DIMENSION(:) :: pheno_model       !! which phenology model is used? (tabulated) (unitless)
!$OMP THREADPRIVATE(pheno_model)
  
  REAL, ALLOCATABLE, SAVE, DIMENSION(:) :: gdd_threshold                 !! Temperature above which growing degrees are accounted for (Celcius)
!$OMP THREADPRIVATE(gdd_threshold)

  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:) :: is_tree                    !! Is the vegetation type a tree ? (true/false)
!$OMP THREADPRIVATE(is_tree)
     
  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:) :: is_deciduous               !! Is PFT deciduous ? (true/false)
!$OMP THREADPRIVATE(is_deciduous)

  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:) :: is_evergreen               !! Is PFT evegreen ? (true/false)
!$OMP THREADPRIVATE(is_evergreen)

  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:) :: is_needleleaf              !! Is PFT needleleaf ? (true/false)
!$OMP THREADPRIVATE(is_needleleaf)
 
  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:) :: is_tropical                !! Is PFT tropical ? (true/false)
!$OMP THREADPRIVATE(is_tropical)

  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:) :: is_temperate               !! Is PFT temperate ? (true/false)
!$OMP THREADPRIVATE(is_temperate)
   
  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:) :: is_boreal                  !! Is PFT boreal ? (true/false)
!$OMP THREADPRIVATE(is_boreal)

  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:) :: natural                    !! natural? (true/false)
!$OMP THREADPRIVATE(natural)

  CHARACTER(len=5), ALLOCATABLE, SAVE, DIMENSION(:) :: type_of_lai       !! Type of behaviour of the LAI evolution algorithm 
                                                                         !! for each vegetation type.
                                                                         !! Value of type_of_lai, one for each vegetation type :
                                                                         !! mean or interp
!$OMP THREADPRIVATE(type_of_lai)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: llaimax                !! laimax for maximum lai see also type of lai 
                                                                         !! interpolation
                                                                         !! @tex $(m^2.m^{-2})$ @endtex 
!$OMP THREADPRIVATE(llaimax)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: llaimin                !! laimin for minimum lai see also type of lai 
                                                                         !! interpolation 
                                                                         !! @tex $(m^2.m^{-2})$ @endtex
!$OMP THREADPRIVATE(llaimin)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: height_presc           !! prescribed height of vegetation.(m)
                                                                         !! Value for height_presc : one for each vegetation type
!$OMP THREADPRIVATE(height_presc)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: z0_over_height        !! Factor to calculate roughness height from 
                                                                        !! vegetation height (unitless)   
!$OMP THREADPRIVATE(z0_over_height)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: ratio_z0m_z0h         !! Ratio between z0m and z0h
!$OMP THREADPRIVATE(ratio_z0m_z0h)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) ::  rveg_pft              !! Potentiometer to set vegetation resistance (unitless)
                                                                         !! Nathalie on March 28th, 2006 - from Fred Hourdin,
!$OMP THREADPRIVATE(rveg_pft)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: sla                    !! specif leaf area @tex $(m^2.gC^{-1})$ @endtex 
!$OMP THREADPRIVATE(sla)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: slainit                !! specif leaf area @tex $(m^2.gC^{-1})$ @endtex 
!$OMP THREADPRIVATE(slainit)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: availability_fact      !! calculate dynamic mortality in lpj_gap
!$OMP THREADPRIVATE(availability_fact)

  !
  ! EVAPOTRANSPIRATION (sechiba)
  !
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: rstruct_const          !! Structural resistance. 
                                                                         !! Value for rstruct_const : one for each vegetation type
                                                                         !! @tex $(s.m^{-1})$ @endtex
!$OMP THREADPRIVATE(rstruct_const)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: kzero                  !! A vegetation dependent constant used in the calculation
                                                                         !! of the surface resistance. 
                                                                         !! Value for kzero one for each vegetation type
                                                                         !! @tex $(kg.m^2.s^{-1})$ @endtex
!$OMP THREADPRIVATE(kzero)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: alpha_watstress_pft     !! PFT-dependent water stress to decide exponential function
!$OMP THREADPRIVATE(alpha_watstress_pft)
  !
  ! WATER (sechiba)
  !
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: wmax_veg               !! Volumetric available soil water capacity in each PFT
                                                                         !! @tex $(kg.m^{-3} of soil)$ @endtex
!$OMP THREADPRIVATE(wmax_veg)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: humcste                !! Root profile description for the different vegetation types.
                                                                         !! These are the factor in the exponential which gets
                                                                         !! the root density as a function of depth 
                                                                         !! @tex $(m^{-1})$ @endtex
!$OMP THREADPRIVATE(humcste)
  
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: max_root_depth         !! Maximum rooting depth for the PFT irrespective of other 
                                                                         !! constraints from the active layer thickness @tex $(m)$ @endtex
!$OMP THREADPRIVATE(max_root_depth)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: throughfall_by_pft     !! Percent by PFT of precip that is not intercepted by the canopy
!$OMP THREADPRIVATE(throughfall_by_pft)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: rsoil_scale            !! Scaling factor to the soil resistance to evaporate per pft
!$OMP THREADPRIVATE(rsoil_scale)

  !
  ! ALBEDO (sechiba)
  !
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: snowa_aged_vis         !! Minimum snow albedo value for each vegetation type
                                                                         !! after aging (dirty old snow), visible albedo (unitless)
                                                                         !! Source : Values are from the Thesis of S. Chalita (1992)
!$OMP THREADPRIVATE(snowa_aged_vis)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: snowa_aged_nir         !! Minimum snow albedo value for each vegetation type
                                                                         !! after aging (dirty old snow), near infrared albedo (unitless)
                                                                         !! Source : Values are from the Thesis of S. Chalita (1992)
!$OMP THREADPRIVATE(snowa_aged_nir)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: snowa_dec_vis          !! Decay rate of snow albedo value for each vegetation type
                                                                         !! as it will be used in condveg_snow, visible albedo (unitless)
                                                                         !! Source : Values are from the Thesis of S. Chalita (1992)
!$OMP THREADPRIVATE(snowa_dec_vis)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: snowa_dec_nir          !! Decay rate of snow albedo value for each vegetation type
                                                                         !! as it will be used in condveg_snow, near infrared albedo (unitless)
                                                                         !! Source : Values are from the Thesis of S. Chalita (1992)
!$OMP THREADPRIVATE(snowa_dec_nir)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: alb_leaf_vis           !! leaf albedo of vegetation type, visible albedo (unitless)
!$OMP THREADPRIVATE(alb_leaf_vis)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: alb_leaf_nir           !! leaf albedo of vegetation type, near infrared albedo (unitless)
!$OMP THREADPRIVATE(alb_leaf_nir)

   REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:) :: leaf_ssa            !! leaf single scattering albedo of all 
!$OMP THREADPRIVATE(leaf_ssa)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:) :: leaf_psd             !! leaf prefered scattering direction of all 
                                                                         !! vegetation types and spectra (unitless)
!$OMP THREADPRIVATE(leaf_psd)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:) :: bgd_reflectance      !! background reflectance of all vegetation types and spectra (unitless)
!$OMP THREADPRIVATE(bgd_reflectance)



  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: leaf_to_shoot_clumping !! The clumping factor for leaves to shoots in the 
                                                                         !! effective LAI calculation...notice this should be
                                                                         !! equal to unity for grasslands/croplands
!$OMP THREADPRIVATE(leaf_to_shoot_clumping)


! NOTE: this next variable originally was a plant-to-stand clumping factor to be used
! in describing how grasses and crops clump together at the plant level (but not trees, 
! as their plant-to-stand clumping is calculated directly) to calculate abledo. However, we
! get the effective spectral parameters for the albedo calculation from inverting satelite data,
! which includes all clumping.  Therefore we do not wish to account for this effect twice.
! What will be incorrect about our grassland and crop albedo is the lack of management options
! for these PFTs in ORCHIDEE, which will lead to LAI values which are wrong.  Thus we 
! will include an LAI correction factor in the calculation of the effective LAI which
! allows us to compensate for this via tuning. 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: lai_correction_factor  !! see note above
!$OMP THREADPRIVATE(lai_correction_factor)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: min_level_sep          !! This is used in determining the levels
                                                                         !! for photosynthesis.  This is the thinnest
                                                                         !! that the levels are allowed to be, in 
                                                                         !! vertical thickness.
                                                                         !! @tex $(m)$ @endtex
!$OMP THREADPRIVATE(min_level_sep)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: lai_top                !!Diffuco.f90 calculates the stomatal conductance of the
                                                                         !! top layer of the canopy. Because the top layer can contain
                                                                         !! different amounts of LAI depending on the crown diameter
                                                                         !! we had to define top layer in terms of the LAI it contains.
                                                                         !! stomatal conductance in the top layer contributes to the
                                                                         !! transpiration (m2 m-2). Arbitrary values.
!$OMP THREADPRIVATE(lai_top) 

  !
  ! SOIL - VEGETATION
  !
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:) :: pref_soil_veg      !! Table which contains the links between the soil
                                                                        !! tiles and the PFTs.
!$OMP THREADPRIVATE(pref_soil_veg)

  !
  ! PHOTOSYNTHESIS
  !
  !-
  ! 1. CO2
  !-
  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:) :: is_c4             !! flag for C4 vegetation types (true/false)
!$OMP THREADPRIVATE(is_c4)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: vcmax_fix     !! values used for vcmax when STOMATE is not activated
                                                                !! @tex $(\mu mol.m^{-2}.s^{-1})$ @endtex
!$OMP THREADPRIVATE(vcmax_fix)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: downregulation_co2_coeff !! Coefficient for CO2 downregulation (unitless)
!$OMP THREADPRIVATE(downregulation_co2_coeff)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: E_KmC         !! Energy of activation for KmC (J mol-1)
!$OMP THREADPRIVATE(E_KmC)                                                               
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: E_KmO         !! Energy of activation for KmO (J mol-1)
!$OMP THREADPRIVATE(E_KmO)         
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: E_Sco         !! Energy of activation for Sco (J mol-1) 
!$OMP THREADPRIVATE(E_Sco) 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: E_gamma_star  !! Energy of activation for gamma_star (J mol-1)
!$OMP THREADPRIVATE(E_gamma_star)    
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: E_Vcmax       !! Energy of activation for Vcmax (J mol-1) 
!$OMP THREADPRIVATE(E_Vcmax)                                                              
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: E_Jmax        !! Energy of activation for Jmax (J mol-1)
!$OMP THREADPRIVATE(E_Jmax) 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: aSV           !! a coefficient of the linear regression (a+bT) defining the Entropy term for Vcmax (J K-1 mol-1)
!$OMP THREADPRIVATE(aSV)    
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: bSV           !! b coefficient of the linear regression (a+bT) defining the Entropy term for Vcmax (J K-1 mol-1 °C-1)
!$OMP THREADPRIVATE(bSV)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: tphoto_min   !! minimum photosynthesis temperature (deg C)
!$OMP THREADPRIVATE(tphoto_min) 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: tphoto_max   !! maximum photosynthesis temperature (deg C)
!$OMP THREADPRIVATE(tphoto_max) 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: aSJ           !! a coefficient of the linear regression (a+bT) defining the Entropy term for Jmax (J K-1 mol-1)
!$OMP THREADPRIVATE(aSJ)    
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: bSJ           !! b coefficient of the linear regression (a+bT) defining the Entropy term for Jmax (J K-1 mol-1 °C-1)
!$OMP THREADPRIVATE(bSJ)    
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: D_Vcmax       !! Energy of deactivation for Vcmax (J mol-1)
!$OMP THREADPRIVATE(D_Vcmax)                     
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: D_Jmax        !! Energy of deactivation for Jmax (J mol-1)
!$OMP THREADPRIVATE(D_Jmax)                           

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: E_gm          !! Energy of activation for gm (J mol-1) 
!$OMP THREADPRIVATE(E_gm)                                       
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: S_gm          !! Entropy term for gm (J K-1 mol-1) 
!$OMP THREADPRIVATE(S_gm)                                       
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: D_gm          !! Energy of deactivation for gm (J mol-1) 
!$OMP THREADPRIVATE(D_gm)                                       
         
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: E_Rd          !! Energy of activation for Rd (J mol-1)
!$OMP THREADPRIVATE(E_Rd)                                      
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: Vcmax25       !! Maximum rate of Rubisco activity-limited carboxylation at 25°C
                                                                !! @tex $(\mu mol.m^{-2}.s^{-1})$ @endtex
!$OMP THREADPRIVATE(Vcmax25)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: arJV          !! a coefficient of the linear regression (a+bT) defining the Jmax25/Vcmax25 ratio (mu mol e- (mu mol CO2)-1)
!$OMP THREADPRIVATE(arJV)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: brJV          !! b coefficient of the linear regression (a+bT) defining the Jmax25/Vcmax25 ratio (mu mol e- (mu mol CO2)-1)
!$OMP THREADPRIVATE(brJV)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: KmC25         !! Michaelis–Menten constant of Rubisco for CO2 at 25°C (ubar)
!$OMP THREADPRIVATE(KmC25)                                     
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: KmO25         !! Michaelis–Menten constant of Rubisco for O2 at 25°C (ubar)
!$OMP THREADPRIVATE(KmO25)     
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: Sco25         !! Relative CO2 /O2 specificity factor for Rubisco at 25°C (bar bar-1) 
!$OMP THREADPRIVATE(Sco25)           
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: gamma_star25  !! Ci-based CO2 compensation point in the absence of Rd at 25°C (ubar)
!$OMP THREADPRIVATE(gamma_star25)   
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: gm25          !! Mesophyll diffusion conductance at 25°C (mol mâ2 sâ1 barâ1) 
!$OMP THREADPRIVATE(gm25)    
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: a1            !! Empirical factor involved in the calculation of fvpd (-)
!$OMP THREADPRIVATE(a1)                                        
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: b1            !! Empirical factor involved in the calculation of fvpd (-)
!$OMP THREADPRIVATE(b1)    
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: a_fpsi        !! Empirical factor involved in the calculation of fpsi (-)
!$OMP THREADPRIVATE(a_fpsi)    
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: sf            !! Sensitivity parameter involved in the calculation of fpsi (-)
!$OMP THREADPRIVATE(sf)    
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: Rroot_sup     !! Superficial roots water resistance involved in the calculation of psi_leaf
!$OMP THREADPRIVATE(Rroot_sup)     
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: Rroot_inf     !! Inferior roots water resistance involved in the calculation of psi_leaf
!$OMP THREADPRIVATE(Rroot_inf)     
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: Rleaf         !! Leaf water resistance involved in the calculation of psi_leaf
!$OMP THREADPRIVATE(Rleaf)     
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: Rxylem        !! Xylem water resistance involved in the calculation of psi_leaf
!$OMP THREADPRIVATE(Rxylem)     
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: Rroot_soil    !! Soil to roots water resistance involved in the calculation of psi_leaf
!$OMP THREADPRIVATE(Rroot_soil)    
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: Rsto_leaf     !! Leaf storage water resistance involved in the calculation of psi_leaf
!$OMP THREADPRIVATE(Rsto_leaf)     
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: Rsto_wood     !! Wood storage water resistance involved in the calculation of psi_leaf
!$OMP THREADPRIVATE(Rsto_wood)  
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: leaf_wat_store_fac !! Empirical factor involved in calculating leaf water storage (-)
!$OMP THREADPRIVATE(leaf_wat_store_fac)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: wood_wat_store_fac !! Empirical factor involved in calculating wood water storage (-)
!$OMP THREADPRIVATE(wood_wat_store_fac) 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: psi_ref_g0    !! Reference water potential involved in the calculation of fpsi in the calculation of stomatal conductance (-)
!$OMP THREADPRIVATE(psi_ref_g0)    
  REAL(r_std), SAVE                            :: psi_min = -20.00    !! Minimum water potential for all organs
!$OMP THREADPRIVATE(psi_min)    
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: psi_ref_sto_leaf    !! Reference water potential involved in the calculation of leaf storage water potential (-)
!$OMP THREADPRIVATE(psi_ref_sto_leaf)     
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: psi_ref_sto_wood    !! Reference water potential involved in the calculation of wood storage water potential (-)
!$OMP THREADPRIVATE(psi_ref_sto_wood)    
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: psi_50_res_leaf    !! Reference water potential involved in the calculation of leaf water resistance (-)
!$OMP THREADPRIVATE(psi_50_res_leaf)       
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: psi_50_res_wood    !! Reference water potential involved in the calculation of wood water resistance (-)
!$OMP THREADPRIVATE(psi_50_res_wood)     
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: psi_50_res_root    !! Reference water potential involved in the calculation of root water resistance (-)
!$OMP THREADPRIVATE(psi_50_res_root)     
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: a_leaf     !! Parameter involved in the calculation of leaf water resistance (-)
!$OMP THREADPRIVATE(a_leaf)              
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: a_wood     !! Parameter involved in the calculation of wood water resistance (-)
!$OMP THREADPRIVATE(a_wood)      
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: a_root     !! Parameter involved in the calculation of root water resistance (-)
!$OMP THREADPRIVATE(a_root)           
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: k_leaf_max     !! Leaf maximum water conductance involved in the calculation of leaf water resistance (-)
!$OMP THREADPRIVATE(k_leaf_max)               
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: k_wood_max     !! Wood maximum water conductance involved in the calculation of wood water resistance (-)
!$OMP THREADPRIVATE(k_wood_max)                  
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: sapwood_density !! Sapwood specific gravity (dry mass/green volume, gC/m3)
!$OMP THREADPRIVATE(sapwood_density)                  
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: ks_max         !! Wood maximum water conductivity involved in the calculation of wood water resistance (m^3/s/MPa/m)
!$OMP THREADPRIVATE(ks_max)                  
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: k_root_max     !! Root maximum water conductance involved in the calculation of root water resistance (m^3/s/MPa)
!$OMP THREADPRIVATE(k_root_max)                     
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: lambda_leaf     !! Parameter involved in the calculation of leaf storage water potential (-)
!$OMP THREADPRIVATE(lambda_leaf)                  
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: lambda_wood     !! Parameter involved in the calculation of wood storage water potential (-)
!$OMP THREADPRIVATE(lambda_wood)                                                                          
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: g0            !! Residual stomatal conductance when irradiance approaches zero (mol m−2 s−1 bar−1)
!$OMP THREADPRIVATE(g0)                                        
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: d_l           !! Parameter for VPD sensitivity of stomatal conductance (kPa)
!$OMP THREADPRIVATE(d_l)                                        
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: h_protons     !! Number of protons required to produce one ATP (mol mol-1)
!$OMP THREADPRIVATE(h_protons)                                 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: fpsir         !! Fraction of PSII e− transport rate partitioned to the C4 cycle (-)
!$OMP THREADPRIVATE(fpsir)                                         
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: fQ            !! Fraction of electrons at reduced plastoquinone that follow the Q-cycle (-) - Values for C3 platns are not used
!$OMP THREADPRIVATE(fQ)                                        
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: fpseudo       !! Fraction of electrons at PSI that follow  pseudocyclic transport (-) - Values for C3 platns are not used
!$OMP THREADPRIVATE(fpseudo)                                   
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: kp            !! Initial carboxylation efficiency of the PEP carboxylase (mol m−2 s−1 bar−1) 
!$OMP THREADPRIVATE(kp)                                        
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: alpha         !! Fraction of PSII activity in the bundle sheath (-)
!$OMP THREADPRIVATE(alpha)                                     
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: gbs           !! Bundle-sheath conductance (mol m−2 s−1 bar−1)
!$OMP THREADPRIVATE(gbs)                                       
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: theta         !! Convexity factor for response of J to irradiance (-)
!$OMP THREADPRIVATE(theta)                                     
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: alpha_LL      !! Conversion efficiency of absorbed light into J at strictly limiting light (mol e− (mol photon)−1)
!$OMP THREADPRIVATE(alpha_LL)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: stress_vcmax  !! Stress on vcmax
!$OMP THREADPRIVATE(stress_vcmax)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: stress_gs     !! Stress on vcmax
!$OMP THREADPRIVATE(stress_gs)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: stress_gm     !! Stress on vcmax
!$OMP THREADPRIVATE(stress_gm)

  !-
  ! 2. Stomate
  !-
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: ext_coeff     !! extinction coefficient of the Monsi&Saeki relationship (1953)
                                                                !! (unitless)
!$OMP THREADPRIVATE(ext_coeff)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: ext_coeff_vegetfrac     !! extinction coefficient used for the calculation of the 
                                                                !! bare soil fraction (unitless)
!$OMP THREADPRIVATE(ext_coeff_vegetfrac)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: ext_coeff_N   !! extinction coefficient of the leaf N content profile within the canopy
                                                                !! ((m2[ground]) (m-2[leaf]))
                                                                !! based on Dewar et al. (2012, value of 0.18), on Carswell et al. (2000, value of 0.11 used in OCN) 
!$OMP THREADPRIVATE(ext_coeff_N)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: nue_opt       !! Nitrogen use efficiency of Vcmax 
                                                                !! ((mumol[CO2] s-1) (gN[leaf])-1)
                                                                !! based on the work of Kattge et al. (2009, GCB)
!$OMP THREADPRIVATE(nue_opt)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:) :: vmax_uptake !! Vmax of nitrogen uptake by plants 
                                                                !! for Ammonium (ind.1) and Nitrate (ind.2)
                                                                !! (in umol (g DryWeight_root)-1 h-1)
                                                                !! from  Kronzucker et al. (1995, 1996)
!$OMP THREADPRIVATE(vmax_uptake)


  !
  ! RESPIRATION (stomate)
  !
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: frac_growthresp  !! fraction of GPP which is lost as growth respiration

!$OMP THREADPRIVATE(frac_growthresp)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: coeff_maint_init    !! maintenance respiration coefficient at 10 deg C,
                                                                      !! @tex $(gC.gN^{-1}.day^{-1})$ @endtex
!$OMP THREADPRIVATE(coeff_maint_init)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: tref_maint_resp     !! maintenance respiration Temperature coefficient,
                                                                      !! @tex $(degC)$ @endtex
!$OMP THREADPRIVATE(tref_maint_resp)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: tmin_maint_resp     !! maintenance respiration Temperature coefficient,
                                                                      !! @tex $(degC)$ @endtex
!$OMP THREADPRIVATE(tmin_maint_resp)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: e0_maint_resp       !! maintenance respiration Temperature coefficient,
                                                                      !! @tex $(unitless)$ @endtex
!$OMP THREADPRIVATE(e0_maint_resp)

 
  !
  ! ALLOCATION
  ! 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: tref_labile         !! Growth from labile pool - temperature at which all labile C will be allocated to growth 
                                                                      !! @tex $(degC)$ @endtex
!$OMP THREADPRIVATE(tref_labile)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: tmin_labile         !! Growth from labile pool  - temperature above which labile will be allocated to growth
                                                                      !! @tex $(degC)$ @endtex
!$OMP THREADPRIVATE(tmin_labile)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: e0_labile           !! Growth temperature coefficient - tuned see stomate_growth_fun_all.f90
                                                                      !! @tex $(unitless)$ @endtex
!$OMP THREADPRIVATE(e0_labile)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: always_labile       !! share of the labile pool that will remain in the labile pool (unitless) 
  
!$OMP THREADPRIVATE(always_labile)

  ! SPITFIRE
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: bulkdensity_deadfuel !! fuel bulk density
!$OMP THREADPRIVATE(bulkdensity_deadfuel)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: bulkdensity_livefuel !! Livefuel bulk density
!$OMP THREADPRIVATE(bulkdensity_livefuel)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: f_sh              !! scorch height parameter for crown fire
!$OMP THREADPRIVATE(f_sh)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: BTpar1            !! Bark thickness parameter
!$OMP THREADPRIVATE(BTpar1)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: BTpar2            !! Bark thickness parameter
!$OMP THREADPRIVATE(BTpar2)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: r_ck              !! parameter for postfire mortality as a result of crown damage
!$OMP THREADPRIVATE(r_ck)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: p_ck              !! parameter for postfire mortality as a result of crown damage
!$OMP THREADPRIVATE(p_ck)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: fdi_sf_pft
!$OMP THREADPRIVATE(fdi_sf_pft)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: ef_CO2            !! emissions factors
!$OMP THREADPRIVATE(ef_CO2)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: ef_CO
!$OMP THREADPRIVATE(ef_CO)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: ef_CH4
!$OMP THREADPRIVATE(ef_CH4)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: ef_VOC
!$OMP THREADPRIVATE(ef_VOC)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: ef_TPM
!$OMP THREADPRIVATE(ef_TPM)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: ef_NOx
!$OMP THREADPRIVATE(ef_NOx)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: me                 !! flammability threshold
!$OMP THREADPRIVATE(me)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: fire_max_cf_100hr       !! Maximum combustion fraction for 100hr and 1000hr fuel for different PFTs
!$OMP THREADPRIVATE(fire_max_cf_100hr)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: fire_max_cf_1000hr       !! Maximum combustion fraction for 100hr and 1000hr fuel for different PFTs
!$OMP THREADPRIVATE(fire_max_cf_1000hr)
  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:) :: burnable
!$OMP THREADPRIVATE(burnable)


  !
  ! FLUX - LUC (Land Use Change)
  !
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: coeff_lcchange_1   !! Coeff of biomass export for the year (unitless)
!$OMP THREADPRIVATE(coeff_lcchange_1)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: coeff_lcchange_10  !! Coeff of biomass export for the decade (unitless)
!$OMP THREADPRIVATE(coeff_lcchange_10)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: coeff_lcchange_100 !! Coeff of biomass export for the century (unitless)
!$OMP THREADPRIVATE(coeff_lcchange_100)
 
  !
  ! FLUX - LUC (Land Use Change)
  !
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: coeff_lcchange_s   !! Coeff of biomass export for the year (unitless)
!$OMP THREADPRIVATE(coeff_lcchange_s)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: coeff_lcchange_m  !! Coeff of biomass export for the decade (unitless)
!$OMP THREADPRIVATE(coeff_lcchange_m)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: coeff_lcchange_l  !! Coeff of biomass export for the century (unitless)
!$OMP THREADPRIVATE(coeff_lcchange_l)
 
  !
  ! PHENOLOGY
  !
  !-
  ! 1. Stomate
  !-
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: lai_max_to_happy  !! threshold of LAI below which plant uses carbohydrate reserves
!$OMP THREADPRIVATE(lai_max_to_happy)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: lai_max           !! maximum LAI, PFT-specific @tex $(m^2.m^{-2})$ @endtex
!$OMP THREADPRIVATE(lai_max)

  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:) :: pheno_type     !! type of phenology (0-4, unitless)
                                                                    !! 0=bare ground 1=evergreen,  2=summergreen, 
                                                                    !! 3=raingreen,  4=perennial
                                                                    !! For the moment, the bare ground phenotype is not managed, 
                                                                    !! so it is considered as "evergreen"
!$OMP THREADPRIVATE(pheno_type)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: force_pheno       !! Number of days after the mean
                                                                    !! doy at which budbreak occurs
                                                                    !! at which phenology will be forced
!$OMP THREADPRIVATE(force_pheno)
  !-
  ! 2. Leaf Onset
  !-
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:) :: pheno_gdd_crit   !! critical gdd,tabulated (C), used in the code
!$OMP THREADPRIVATE(pheno_gdd_crit)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: pheno_gdd_crit_c   !! critical gdd,tabulated (C), 
                                                                     !! constant c of aT^2+bT+c (unitless)
!$OMP THREADPRIVATE(pheno_gdd_crit_c)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: pheno_gdd_crit_b   !! critical gdd,tabulated (C), 
                                                                     !! constant b of aT^2+bT+c (unitless)
!$OMP THREADPRIVATE(pheno_gdd_crit_b)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: pheno_gdd_crit_a   !! critical gdd,tabulated (C), 
                                                                     !! constant a of aT^2+bT+c (unitless)
!$OMP THREADPRIVATE(pheno_gdd_crit_a)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: pheno_moigdd_t_crit!! Monthly avearage temperature treashold used for C4 grass (C)
!$OMP THREADPRIVATE(pheno_moigdd_t_crit)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: ngd_crit           !! critical ngd,tabulated. Threshold -5 degrees (days)
!$OMP THREADPRIVATE(ngd_crit)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: ncdgdd_temp        !! critical temperature for the ncd vs. gdd function
                                                                     !! in phenology (C)
!$OMP THREADPRIVATE(ncdgdd_temp)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: hum_frac           !! critical humidity (relative to min/max) for phenology
                                                                     !! (0-1, unitless)
!$OMP THREADPRIVATE(hum_frac)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: hum_min_time       !! minimum time elapsed since moisture minimum (days)
!$OMP THREADPRIVATE(hum_min_time)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: relsoilmoist_always!! moisture monthly availability above which moisture tendency
                                                                     !! doesn't matter (0-1, unitless)
!$OMP THREADPRIVATE(relsoilmoist_always)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: t_always           !! temperature monthly above which temperature tendency (K)
                                                                     !! doesn't matter
!$OMP THREADPRIVATE(t_always)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: harvest_time_a     !! Parameter to calculate harvest time of crops (day)
!$OMP THREADPRIVATE(harvest_time_a)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) ::  harvest_time_b    !! Parameter to calculate harvest time of crops (K)
!$OMP THREADPRIVATE(harvest_time_b)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: longevity_sap      !! sapwood -> heartwood conversion time (days)
!$OMP THREADPRIVATE(longevity_sap)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: longevity_leaf     !! leaf turnover (1/years)
!$OMP THREADPRIVATE(longevity_leaf)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: leaf_age_crit_tref !! Reference temperature of the PFT (degrees Celsius)
                                                                     !! Used to calculate the leaf_age_crit as a function of longevity
!$OMP THREADPRIVATE(leaf_age_crit_tref)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: leaf_age_crit_coeff1 !! Coeff1 (unitless) to link leaf_age_crit to leaf_age_crit_tref
!$OMP THREADPRIVATE(leaf_age_crit_coeff1)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: leaf_age_crit_coeff2 !! Coeff2 (unitless) to link leaf_age_crit to leaf_age_crit_tref
!$OMP THREADPRIVATE(leaf_age_crit_coeff2)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: leaf_age_crit_coeff3 !! Coeff3 (unitless) to link leaf_age_crit to leaf_age_crit_tref
!$OMP THREADPRIVATE(leaf_age_crit_coeff3)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: longevity_fruit    !! fruit lifetime (days)
!$OMP THREADPRIVATE(longevity_fruit)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: longevity_root     !! root turnover (1/days)
!$OMP THREADPRIVATE(longevity_root)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: ecureuil           !! fraction of primary leaf and root allocation put
                                                                     !! into reserve (0-1, unitless)
!$OMP THREADPRIVATE(ecureuil)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: alloc_min          !! NEW - allocation above/below = f(age) - 30/01/04 NV/JO/PF
!$OMP THREADPRIVATE(alloc_min)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: alloc_max          !! NEW - allocation above/below = f(age) - 30/01/04 NV/JO/PF
!$OMP THREADPRIVATE(alloc_max)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: demi_alloc         !! NEW - allocation above/below = f(age) - 30/01/04 NV/JO/PF
!$OMP THREADPRIVATE(demi_alloc)

  !-
  ! 3. Senescence
  !-
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: leaffall              !! length of death of leaves,tabulated (days)
!$OMP THREADPRIVATE(leaffall)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: presenescence_ratio   !! The ratio of maintenance respiration to
                                                                        !! gpp beyond which presenescence stage of
                                                                        !! plant phenology is declared to begin (0-1, unitless).
!$OMP THREADPRIVATE(presenescence_ratio)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: senescence_ratio      !! The ratio of maintenance respiration to
                                                                        !! gpp beyond which senescence stage of
                                                                        !! plant phenology is declared to begin (0-1, unitless).
!$OMP THREADPRIVATE(senescence_ratio)

  CHARACTER(len=6), ALLOCATABLE, SAVE, DIMENSION(:) :: senescence_type  !! type of senescence,tabulated (unitless)
                                                                        !! List of avaible types of senescence :
                                                                        !! 'cold  ', 'dry   ', 'mixed ', 'none  '
!$OMP THREADPRIVATE(senescence_type)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: senescence_hum        !! critical relative moisture availability for senescence
                                                                        !! (0-1, unitless)
!$OMP THREADPRIVATE(senescence_hum)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: nosenescence_hum      !! relative moisture availability above which there is
                                                                        !! no humidity-related senescence (0-1, unitless)
!$OMP THREADPRIVATE(nosenescence_hum)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: max_turnover_time     !! maximum turnover time for grasses (days)
!$OMP THREADPRIVATE(max_turnover_time)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: min_turnover_time     !! minimum turnover time for grasses (days)
!$OMP THREADPRIVATE(min_turnover_time)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: recycle_leaf          !! Fraction of N leaf that is recycled when leaves are senescent
!$OMP THREADPRIVATE(recycle_leaf)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: recycle_root          !! Fraction of N root that is recycled when leaves are senescent
!$OMP THREADPRIVATE(recycle_root)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: min_leaf_age_for_senescence  !! minimum leaf age to allow senescence g (days)
!$OMP THREADPRIVATE(min_leaf_age_for_senescence)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:) :: senescence_temp     !! critical temperature for senescence (C),
                                                                        !! used in the code
!$OMP THREADPRIVATE(senescence_temp)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: senescence_temp_c     !! critical temperature for senescence (C), 
                                                                        !! constant c of aT^2+bT+c , tabulated (unitless)
!$OMP THREADPRIVATE(senescence_temp_c)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: senescence_temp_b     !! critical temperature for senescence (C), 
                                                                        !! constant b of aT^2+bT+c , tabulated (unitless)
!$OMP THREADPRIVATE(senescence_temp_b)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: senescence_temp_a     !! critical temperature for senescence (C),
                                                                        !! constant a of aT^2+bT+c , tabulated (unitless)
!$OMP THREADPRIVATE(senescence_temp_a)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: gdd_senescence        !! minimum gdd to allow senescence of crops (days)
!$OMP THREADPRIVATE(gdd_senescence)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: pre_gdd        !! minimum gdd to allow senescence of crops (days)
!$OMP THREADPRIVATE(pre_gdd)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: pre_harv        !! minimum gdd to allow senescence of crops (days)
!$OMP THREADPRIVATE(pre_harv)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: gr_gt        !! minimum gdd to allow senescence of crops (days)
!$OMP THREADPRIVATE(gr_gt)

  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:) :: always_init               !! take carbon from atmosphere if carbohydrate reserve too small? (true/false)
!$OMP THREADPRIVATE(always_init)

  !-
  ! 4. N cycle
  !-
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: cn_leaf_min           !! minimum CN ratio of leaves (gC/gN)
!$OMP THREADPRIVATE(cn_leaf_min)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: cn_leaf_max           !! maximum CN ratio of leaves (gC/gN) 
!$OMP THREADPRIVATE(cn_leaf_max)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: max_soil_n_bnf        !! Value of total N (NH4+NO3) 
                                                                        !! above which we stop adding N via BNF
                                                                        !! (gN/m**2)
!$OMP THREADPRIVATE(max_soil_n_bnf)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: manure_pftweight      ! Weight of the distribution of manure over the PFT surface
!$OMP THREADPRIVATE(manure_pftweight)

  !
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)    :: harvest_ratio  !! Share of biomass that is removed from the site during harvest
                                                                    !! A high value indicates a high harvest efficiency and thus a 
                                                                    !! input of residuals. (unitless, 0-1).

!$OMP THREADPRIVATE(harvest_ratio)

  !
  ! STOMATE - Age classes
  !

  INTEGER(i_std), SAVE                            :: nvmap          !! The number of PFTs we have if we ignore age classes.
                                                                    !! @tex $-$ @endtex
!$OMP THREADPRIVATE(nvmap) 
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:) :: agec_group     !! The age class group that this PFT belongs to.
                                                                    !! If you're not using age classes, this will just be
                                                                    !! set to the number of the PFT and should be ignored
                                                                    !! in the code.
                                                                    !! @tex $-$ @endtex
!$OMP THREADPRIVATE(agec_group) 
! I do not like the location of these next two variables.  They are computed
! after agec_group is read in.  Ideally, they would be passed around
! as arguments or in a structure, since they are not really
! parameters read in from the input file.
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:) :: start_index    !! Gives the index that this real PFT starts
                                                                    !! on, ignoring age classes
                                                                    !! @tex $-$ @endtex
!$OMP THREADPRIVATE(start_index) 
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:) :: nagec_pft      !! The number of age classes for each PFT.
                                                                    !! Only 1 or nagec are supported right now.
                                                                    !! @tex $-$ @endtex
!$OMP THREADPRIVATE(nagec_pft) 

  !
  ! DGVM
  !

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: residence_time        !! residence time of trees (y) 
!$OMP THREADPRIVATE(residence_time)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: tmin_crit             !! critical tmin, tabulated (C)
!$OMP THREADPRIVATE(tmin_crit)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: tcm_crit              !! critical tcm, tabulated (C)
!$OMP THREADPRIVATE(tcm_crit)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)   :: mortality_min     !! Asymptotic mortality if plant growth exceeds long term
                                                                    !! NPP @tex $(year^{-1})$ @endtex
!$OMP THREADPRIVATE(mortality_min) 

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)   :: mortality_max     !! Maximum mortality if plants hardly grows
                                                                    !! @tex $(year^{-1})$ @endtex
!$OMP THREADPRIVATE(mortality_max) 

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)   :: ref_mortality     !! Reference mortality rate used to calculate mortality
                                                                    !! as a function of the plant vigor @tex $(year^{-1})$ @endtex
!$OMP THREADPRIVATE(ref_mortality)

  !
  ! Biogenic Volatile Organic Compounds
  !

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_isoprene       !! Isoprene emission factor
                                                                           !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_isoprene)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_monoterpene    !! Monoterpene emission factor
                                                                           !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_monoterpene)

  REAL(r_std), SAVE :: LDF_mono                                            !! monoterpenes fraction dependancy to light
!$OMP THREADPRIVATE(LDF_mono) 
  REAL(r_std), SAVE :: LDF_sesq                                            !! sesquiterpenes fraction dependancy to light
!$OMP THREADPRIVATE(LDF_sesq)
  REAL(r_std), SAVE :: LDF_meth                                            !! methanol fraction dependancy to light
!$OMP THREADPRIVATE(LDF_meth)
  REAL(r_std), SAVE :: LDF_acet                                            !! acetone fraction dependancy to light
!$OMP THREADPRIVATE(LDF_acet)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_apinene        !! Alfa pinene emission factor
                                                                           !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_apinene)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_bpinene        !! Beta pinene emission factor
                                                                           !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_bpinene)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_limonene       !! Limonene emission factor
                                                                           !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_limonene)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_myrcene        !! Myrcene emission factor
                                                                           !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_myrcene)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_sabinene       !! Sabinene emission factor
                                                                           !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_sabinene)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_camphene       !! Camphene emission factor
                                                                           !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_camphene)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_3carene        !! 3-carene emission factor
                                                                           !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_3carene)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_tbocimene      !! T-beta-ocimene emission factor
                                                                           !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_tbocimene)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_othermonot     !! Other monoterpenes emission factor
                                                                           !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_othermonot)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_sesquiterp     !! Sesquiterpene emission factor
                                                                           !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_sesquiterp)

  REAL(r_std), SAVE :: beta_mono                                           !! Monoterpenes temperature dependency coefficient 
!$OMP THREADPRIVATE(beta_mono)
  REAL(r_std), SAVE :: beta_sesq                                           !! Sesquiterpenes temperature dependency coefficient
!$OMP THREADPRIVATE(beta_sesq)
  REAL(r_std), SAVE :: beta_meth                                           !! Methanol temperature dependency coefficient 
!$OMP THREADPRIVATE(beta_meth)
  REAL(r_std), SAVE :: beta_acet                                           !! Acetone temperature dependency coefficient
!$OMP THREADPRIVATE(beta_acet)
  REAL(r_std), SAVE :: beta_oxyVOC                                         !! Other oxygenated BVOC temperature dependency coefficient
!$OMP THREADPRIVATE(beta_oxyVOC)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_ORVOC          !! ORVOC emissions factor
                                                                           !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_ORVOC)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_OVOC           !! OVOC emissions factor
                                                                           !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_OVOC)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_MBO            !! MBO emissions factor 
                                                                           !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_MBO)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_methanol       !! Methanol emissions factor
                                                                           !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_methanol)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_acetone        !! Acetone emissions factor 
                                                                           !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_acetone)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_acetal         !! Acetaldehyde emissions factor
                                                                           !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_acetal)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_formal         !! Formaldehyde emissions factor
                                                                           !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_formal)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_acetic         !! Acetic Acid emissions factor 
                                                                           !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_acetic)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_formic         !! Formic Acid emissions factor 
                                                                           !! @tex $(\mu gC.g^{-1}.h^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_formic)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_DMS            !! DMS emissions factor
                                                                           !! @tex $(\mu g.g^{-1}.h^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_DMS)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_H2S            !! H2S emissions factor
                                                                           !! @tex $(\mu g.g^{-1}.h^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_H2S)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_no_wet         !! NOx emissions factor soil emissions and 
                                                                           !! exponential dependancy factor for wet soils
                                                                           !! @tex $(ngN.m^{-2}.s^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_no_wet)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: em_factor_no_dry         !! NOx emissions factor soil emissions and
                                                                           !! exponential dependancy factor for dry soils
                                                                           !! @tex $(ngN.m^{-2}.s^{-1})$ @endtex
!$OMP THREADPRIVATE(em_factor_no_dry)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: Larch                    !! Larcher 1991 SAI/LAI ratio (unitless)
!$OMP THREADPRIVATE(Larch)

  !
  ! INTERNAL PARAMETERS USED IN STOMATE_DATA
  !

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: lai_initmin   !! Initial lai for trees/grass 
                                                                !! @tex $(m^2.m^{-2})$ @endtex
!$OMP THREADPRIVATE(lai_initmin)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: bm_sapl   !! sapling biomass @tex $(gC.ind^{-1})$ @endtex
!$OMP THREADPRIVATE(bm_sapl)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: migrate       !! migration speed @tex $(m.year^{-1})$ @endtex
!$OMP THREADPRIVATE(migrate)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: cn_sapl       !! crown of tree when sapling  @tex $(m^2$)$ @endtex
!$OMP THREADPRIVATE(cn_sapl)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: k_latosa_max       !! Maximum leaf-to-sapwood area ratio (unitless)
!$OMP THREADPRIVATE(k_latosa_max)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: k_latosa_min       !! Minimum leaf-to-sapwood area ratio (unitless)
!$OMP THREADPRIVATE(k_latosa_min)

 REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:) :: LC                 !! Lignin/C ratio of the different biomass pools and PFTs (unitless)
                                                                      !! based on CN from White et al. (2000) 
!$OMP THREADPRIVATE(LC)

 REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: LC_leaf              !! Lignin/C ratio of leaf pool (unitless)
                                                                      !! based on CN from White et al. (2000) 
!$OMP THREADPRIVATE(LC_leaf)

 REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: LC_sapabove          !! Lignin/C ratio of sapabove pool (unitless)
                                                                      !! based on CN from White et al. (2000) 
!$OMP THREADPRIVATE(LC_sapabove)

 REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: LC_sapbelow          !! Lignin/C ratio of sapbelow pool (unitless)
                                                                      !! based on CN from White et al. (2000) 
!$OMP THREADPRIVATE(LC_sapbelow)

 REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: LC_heartabove        !! Lignin/C ratio of heartabove pool (unitless)
                                                                      !! based on CN from White et al. (2000) 
!$OMP THREADPRIVATE(LC_heartabove)

 REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: LC_heartbelow        !! Lignin/C ratio of heartbelow pool (unitless)
                                                                      !! based on CN from White et al. (2000) 
!$OMP THREADPRIVATE(LC_heartbelow)

 REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: LC_fruit             !! Lignin/C ratio of fruit pool (unitless)
                                                                      !! based on CN from White et al. (2000) 
!$OMP THREADPRIVATE(LC_fruit)

 REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: LC_root              !! Lignin/C ratio of root pool (unitless)
                                                                      !! based on CN from White et al. (2000) 
!$OMP THREADPRIVATE(LC_root)

 REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: LC_carbres           !! Lignin/C ratio of carbres pool (unitless)
                                                                      !! based on CN from White et al. (2000) 
!$OMP THREADPRIVATE(LC_carbres)

 REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: LC_labile            !! Lignin/C ratio of labile pool (unitless)
                                                                      !! based on CN from White et al. (2000) 
!$OMP THREADPRIVATE(LC_labile)

 REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: decomp_factor        !! Multpliactive factor modifying 
                                                                      !! the standard decomposition factor for each SOM pool
!$OMP THREADPRIVATE(decomp_factor)

  !
  ! STAND STRUCTURE (stomate) 
  !
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: pipe_density        !! Wood density in @tex $(gC.m^{-3})$ @endtex
!$OMP THREADPRIVATE(pipe_density)
  
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: min_pipe_tune2      !! Initial value for pipe_tune2 when it is calculated as a function of precipitation (m)
!$OMP THREADPRIVATE(min_pipe_tune2)
  
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: pipe_tune3          !! height=pipe_tune2 * diameter**pipe_tune3
!$OMP THREADPRIVATE(pipe_tune3)      
  
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: pipe_tune4          !! ???needed for stem diameter
!$OMP THREADPRIVATE(pipe_tune4)      
  
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: pipe_k1             !! ???
!$OMP THREADPRIVATE(pipe_k1)       
  
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: tree_ff             !! Volume reduction factor from cylinder to real tree shape (inc.branches)
!$OMP THREADPRIVATE(tree_ff)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: crown_to_height    !! Ratio between tree height and the vertical crown diameter. If this value is changed check beforehand that the crown diameter will never exceed the tree height.
!$OMP THREADPRIVATE(crown_to_height)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: crown_vertohor_dia !!Ratio between the vertical and horizontal crown diameter height, so indirectly the horizonatl crown diameter also depends on crown diameter
!$OMP THREADPRIVATE(crown_vertohor_dia)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: mass_ratio_heart_sap!! mass ratio (heartwood+sapwood)/heartwood 
!$OMP THREADPRIVATE(mass_ratio_heart_sap)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: canopy_cover        !! Canopy cover - current values are guesses for testing
                                                                      !! could tune this variable to match MODIS albedo
!$OMP THREADPRIVATE(canopy_cover)
   
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: nmaxplants           !! Intial number of seedlings per hectare. Used 
                                                                      !! in prescribe to initialize the model and after 
                                                                      !! every clearcut
!$OMP THREADPRIVATE(nmaxplants)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: height_init         !! The height when a new grassland or cropland is established  
                                                                      !! @tex $(m)$ @endtex
!$OMP THREADPRIVATE(height_init)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: biomass_init        !! The biomass when a new grassland or cropland is established
!$OMP THREADPRIVATE(biomass_init)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: qmd_init        !! Min diameter use to initilize stand in prescribe
!$OMP THREADPRIVATE(qmd_init)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: dia_thresh_inv      !! Diameter below which trees are not measured in a forest inventory (m)
!$OMP THREADPRIVATE(dia_thresh_inv)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: k_root              !! Fine root specific conductivity 
                                                                      !! @tex $(m^{3} kg^{-1} s^{-1} MPa^{-1})$ @endtex 
!$OMP THREADPRIVATE(k_root)
 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: k_belowground       !! Belowground (root + soil) specific conductivity 
                                                                      !! @tex $(m^{3} kg^{-1} s^{-1} MPa^{-1})$ @endtex 
!$OMP THREADPRIVATE(k_belowground)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: k_sap               !! Sapwood specific conductivity
                                                                      !! @tex $(m^{3} kg^{-1} s^{-1} MPa^{-1})$ @endtex 
!$OMP THREADPRIVATE(k_sap)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: psi_leaf_min            !! Minimal leaf water potential @tex $(m s^{-1} MPa^{-1})$ @endtex
!$OMP THREADPRIVATE(psi_leaf_min)

 REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: srl                  !! Specific root length  @tex $(m g^{-1})$ @endtex       
!$OMP THREADPRIVATE(srl)

 REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: r_froot              !! Radius of fine roots  @tex $(m)$ @endtex       
!$OMP THREADPRIVATE(r_froot)

 REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: psi_root_min         !! Minimum root water potential  @tex $(MPa)$ @endtex       
!$OMP THREADPRIVATE(psi_root_min)

  !
  ! SOIL DECOMPOSITION
  ! 
REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: som_turn_ipassive     !! Turnover of the passive pool (year-1)
!$OMP THREADPRIVATE(som_turn_ipassive)
 
  !  
  ! RECRUITMENT (stomate) 
  !  
  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:) :: recruitment_pft                 !! Do recruitment? (true/false)  
!$OMP THREADPRIVATE(recruitment_pft)   

  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:) :: beetle_pft                 !!Do beetle? (true/false)  
!$OMP THREADPRIVATE(beetle_pft)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: recruitment_height          !! Height of recruits  
!$OMP THREADPRIVATE(recruitment_height)  

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: recruitment_alpha           !! Alpha parameter for recruitment model  
!$OMP THREADPRIVATE(recruitment_alpha)  

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: recruitment_beta            !! Beta parameter for recruitment model  
!$OMP THREADPRIVATE(recruitment_beta) 

  !
  ! PRESCRIBE (stomate)
  !
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: tune_reserves_in_sapling  !! A factor to scale the reserve pool of newly
                                                                            !!planted saplings.  This is required by some deciduous
                                                                            !! trees in order to survive the first year until budburst,
                                                                            !!but it has no physical basis.
                                                                            !!(unitless)
    
!$OMP THREADPRIVATE(tune_reserves_in_sapling)

!
  ! MORTALITY (stomate)
  !
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) &
                          :: death_distribution_factor              !! The scale factor between the smallest and largest 
                                                                    !! circ class for tree mortality in lpj_kill.
                                                                    !! (unitless)
!$OMP THREADPRIVATE(death_distribution_factor)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) &
                          :: npp_reset_value                        !! The value of the NPP that the long-term value is
                                                                    !! reset to after a PFT dies in stomate_kill.  This
                                                                    !! only seems to be used for non-trees.
                                                                    !! @tex $(gC m^{-2})$ @endtex
!$OMP THREADPRIVATE(npp_reset_value)

REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)   :: ndying_year       !! Reference number of year in which the forest will disappear after reaching the the stem density threshold.  
!$OMP THREADPRIVATE(ndying_year)

  ! BEETLE (stomate_pest)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: S_activity  !! fraction of the beetle population which remaing on the stand after beettle departure 
  !$OMP THREADPRIVATE(S_activity)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: act_limit  !! parameter wich increase/reduce the importance of the previous infestation in the calculation of the BPI (beetle pressure index)
  !$OMP THREADPRIVATE(act_limit)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: S_competition  !! a parameter for the relationship between stand rdi and beetle susceptibility
  !$OMP THREADPRIVATE(S_competition)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: RDi_limit  !! 
  !$OMP THREADPRIVATE(RDi_limit)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: S_susceptibility  !! 
  !$OMP THREADPRIVATE(S_susceptibility)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: i_rd_susceptibility  !! b parameter for the relationship between stand rdi and beetle susceptibility
  !$OMP THREADPRIVATE(i_rd_susceptibility)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: S_share  !! a parameter for the relationship between stand share and beetle susceptibility
  !$OMP THREADPRIVATE(S_share)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: SH_limit  !! b parameter for the relationship between stand share and beetle susceptibility
  !$OMP THREADPRIVATE(SH_limit)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: S_defence  !! a parameter for the relationship between drought and beetle susceptibility
  !$OMP THREADPRIVATE(S_defence)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: PWS_limit  !! b parameter for the relationship between drought and beetle susceptibility
  !$OMP THREADPRIVATE(PWS_limit)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: max_Nwood  !! parameter for the relationship between windthrow and beetle susceptibility
  !$OMP THREADPRIVATE(max_Nwood)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: beetle_generation_a  !! a parameter for the calculation of the number of beetle generation per year
  !$OMP THREADPRIVATE(beetle_generation_a)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: beetle_generation_b  !! b parameter for the calculation of the number of beetle generation per year
  !$OMP THREADPRIVATE(beetle_generation_b)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: beetle_generation_c  !! c parameter for the calculation of the number of beetle generation per year
  !$OMP THREADPRIVATE(beetle_generation_c)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: min_temp_beetle  !! temperature threshold below which Teff is not calculated (*C)
  !$OMP THREADPRIVATE(min_temp_beetle)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: max_temp_beetle  !! temperature threshold above which Teff is not calculated (*C)
  !$OMP THREADPRIVATE(max_temp_beetle)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: opt_temp_beetle  !! optimal temperature to breed bark beetle (*C)
  !$OMP THREADPRIVATE(opt_temp_beetle)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: eff_temp_beetle_a  !! a parameter for the calculation of the effective temperature used in beetle phenology 
  !$OMP THREADPRIVATE(eff_temp_beetle_a)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: eff_temp_beetle_b  !! b parameter for the calculation of the effective temperature used in beetle phenology 
  !$OMP THREADPRIVATE(eff_temp_beetle_b)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: eff_temp_beetle_c  !! c parameter for the calculation of the effective temperature used in beetle phenology
  !$OMP THREADPRIVATE(eff_temp_beetle_c)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: eff_temp_beetle_d  !! d parameter for the calculation of the effective temperature used in beetle phenology
  !$OMP THREADPRIVATE(eff_temp_beetle_d)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: diapause_thres_daylength  !! daylength in hour above which bark beetle start diapause 
  !$OMP THREADPRIVATE(diapause_thres_daylength)
 REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: S_massattack ! a parameter of the weight of rdi susceptibility in the calculation of landscape susceptibility
 !$OMP THREADPRIVATE(S_massattack)
 REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: BP_limit ! b parameter of the weight of rdi susceptibility in the calculation of landscape susceptibility
 !$OMP THREADPRIVATE(BP_limit)

  !
  ! WINDFALL (stomate_windthrow)
  !
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: streamlining_c_leaf                         !! Modulus of Rupture (unitless)
!$OMP THREADPRIVATE(streamlining_c_leaf)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: streamlining_c_leafless                     !! Modulus of Rupture (unitless)
!$OMP THREADPRIVATE(streamlining_c_leafless)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: streamlining_n_leaf                         !! Modulus of Rupture (unitless)
!$OMP THREADPRIVATE(streamlining_n_leaf)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: streamlining_n_leafless                     !! Modulus of Rupture (unitless)
!$OMP THREADPRIVATE(streamlining_n_leafless)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: modulus_rupture                             !! Modulus of Rupture (Pa)
!$OMP THREADPRIVATE(modulus_rupture)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: f_knot                                      !! Modulus of Rupture (unitless)
!$OMP THREADPRIVATE(f_knot)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: overturning_free_draining_shallow           !! Modulus of Rupture (Nm.kg-1)
!$OMP THREADPRIVATE(overturning_free_draining_shallow)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: overturning_free_draining_shallow_leafless  !! Modulus of Rupture (Nm.kg-1)
!$OMP THREADPRIVATE(overturning_free_draining_shallow_leafless)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: overturning_free_draining_deep              !! Modulus of Rupture (Nm.kg-1)
!$OMP THREADPRIVATE(overturning_free_draining_deep)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: overturning_free_draining_deep_leafless     !! Modulus of Rupture (Nm.kg-1)
!$OMP THREADPRIVATE(overturning_free_draining_deep_leafless)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: overturning_free_draining_average           !! Modulus of Rupture (Nm.kg-1)
!$OMP THREADPRIVATE(overturning_free_draining_average)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: overturning_free_draining_average_leafless  !! Modulus of Rupture (Nm.kg-1)
!$OMP THREADPRIVATE(overturning_free_draining_average_leafless)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: overturning_gleyed_shallow                  !! Modulus of Rupture (Nm.kg-1)
!$OMP THREADPRIVATE(overturning_gleyed_shallow)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: overturning_gleyed_shallow_leafless         !! Modulus of Rupture (Nm.kg-1)
!$OMP THREADPRIVATE(overturning_gleyed_shallow_leafless)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: overturning_gleyed_deep                     !! Modulus of Rupture (Nm.kg-1)
!$OMP THREADPRIVATE(overturning_gleyed_deep)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: overturning_gleyed_deep_leafless            !! Modulus of Rupture (Nm.kg-1)
!$OMP THREADPRIVATE(overturning_gleyed_deep_leafless)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: overturning_gleyed_average                  !! Modulus of Rupture (Nm.kg-1)
!$OMP THREADPRIVATE(overturning_gleyed_average)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: overturning_gleyed_average_leafless         !! Modulus of Rupture (Nm.kg-1)
!$OMP THREADPRIVATE(overturning_gleyed_average_leafless)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: overturning_peaty_shallow                   !! Modulus of Rupture (Nm.kg-1)
!$OMP THREADPRIVATE(overturning_peaty_shallow)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: overturning_peaty_shallow_leafless          !! Modulus of Rupture (Nm.kg-1)
!$OMP THREADPRIVATE(overturning_peaty_shallow_leafless)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: overturning_peaty_deep                      !! Modulus of Rupture (Nm.kg-1)
!$OMP THREADPRIVATE(overturning_peaty_deep)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: overturning_peaty_deep_leafless             !! Modulus of Rupture (Nm.kg-1)
!$OMP THREADPRIVATE(overturning_peaty_deep_leafless)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: overturning_peaty_average                   !! Modulus of Rupture (Nm.kg-1)
!$OMP THREADPRIVATE(overturning_peaty_average)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: overturning_peaty_average_leafless          !! Modulus of Rupture (Nm.kg-1)
!$OMP THREADPRIVATE(overturning_peaty_average_leafless)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: overturning_peat_shallow                    !! Modulus of Rupture (Nm.kg-1)
!$OMP THREADPRIVATE(overturning_peat_shallow)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: overturning_peat_shallow_leafless           !! Modulus of Rupture (Nm.kg-1)
!$OMP THREADPRIVATE(overturning_peat_shallow_leafless)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: overturning_peat_deep                       !! Modulus of Rupture (Nm.kg-1)
!$OMP THREADPRIVATE(overturning_peat_deep)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: overturning_peat_deep_leafless              !! Modulus of Rupture (Nm.kg-1)
!$OMP THREADPRIVATE(overturning_peat_deep_leafless)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: overturning_peat_average                    !! Modulus of Rupture (Nm.kg-1)
!$OMP THREADPRIVATE(overturning_peat_average)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: overturning_peat_average_leafless           !! Modulus of Rupture (Nm.kg-1)
!$OMP THREADPRIVATE(overturning_peat_average_leafless)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: max_damage_further                          !! Maximum damage rate for inner area (unitless)  
!$OMP THREADPRIVATE(max_damage_further)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: max_damage_closer                           !! Maximum damage rate for forest edge (unitless)
!$OMP THREADPRIVATE(max_damage_closer)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: sfactor_further                             !! Scaling coefficient s for inner forest (unitless)
!$OMP THREADPRIVATE(sfactor_further)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: sfactor_closer                              !! Scaling coefficient s for forest edge (unitless)
!$OMP THREADPRIVATE(sfactor_closer)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: green_density                               !! (kg/m^3)
!$OMP THREADPRIVATE(green_density)

 

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: lai_to_height       !! Covert lai into vegetation height for grasses and crops
!$OMP THREADPRIVATE(lai_to_height)    

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)    :: deleuze_a      !! intercept of the intra-tree competition within a stand
                                                                    !! based on the competion rule of Deleuze and Dhote 2004
                                                                    !! Used when n_circ > 6
!$OMP THREADPRIVATE(deleuze_a)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)    :: deleuze_b      !! slope of the intra-tree competition within a stand
                                                                    !! based on the competion rule of Deleuze and Dhote 2004
                                                                    !! Used when n_circ > 6
!$OMP THREADPRIVATE(deleuze_b)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)    :: deleuze_p_all  !! Percentile of the circumferences that receives photosynthates
                                                                    !! based on the competion rule of Deleuze and Dhote 2004
                                                                    !! Used when n_circ > 6 for FM1, FM2 and FM4 
!$OMP THREADPRIVATE(deleuze_p_all)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)    :: deleuze_p_coppice  !! Percentile of the circumferences that receives photosynthates
                                                                    !! based on the competion rule of Deleuze and Dhote 2004
                                                                    !! Used when n_circ > 6 for FM3
!$OMP THREADPRIVATE(deleuze_p_coppice)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)    :: deleuze_power_a!! Divisor of the power for the slope of intra-tree competition within a stand
                                                                    !! based on the competition rule of Deleuze and Dhote 2004
!$OMP THREADPRIVATE(deleuze_power_a)


  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: m_dv                !! Parameter in the Deleuze & Dhote allocation rule that
                                                                      !! relaxes the cut-off imposed by ::sigma. Owing to m_relax 
                                                                      !! trees still grow a little when their ::circ is below 
                                                                      !! ::sigma
!$OMP THREADPRIVATE(m_dv) 
 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: beta_self_thinning!! Exponent of the self-thinning relationship D=alpha*N^beta 
                                                                    !! estimated from German, French, Spanish and swedish
                                                                    !! forest inventories
!$OMP THREADPRIVATE(beta_self_thinning)
  
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: ref_alpha_self_thin !! Initial alpha self thinning parameter
!$OMP THREADPRIVATE(ref_alpha_self_thin)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:) :: pre_indust_ref_gpp !! pre-industrial reference gpp that used to scale alpha_self_thinning
!$OMP THREADPRIVATE(pre_indust_ref_gpp)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: init_pre_indust_ref_gpp !! Initial pre-industrial reference gpp that used to scale alpha_self_thinning
!$OMP THREADPRIVATE(init_pre_indust_ref_gpp)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: thinstrat         !! The thinning strategy used for forest management.
                                                                    !! Comes from Eq. 12 of Bellassen et al (2010)
                                                                    !! @tex $(unitless)$ @endtex
!$OMP THREADPRIVATE(thinstrat) 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: taumin            !! Minimum tree death probability. stomate_forest.f90
                                                                    !! Comes from Eq. 12 of Bellassen et al (2010)
                                                                    !! @tex $(unitless)$ @endtex
!$OMP THREADPRIVATE(taumin) 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: taumax            !! Maximum tree death probability. stomate_forest.f90
                                                                    !! Comes from Eq. 12 of Bellassen et al (2010)
                                                                    !! @tex $(unitless)$ @endtex
!$OMP THREADPRIVATE(taumax)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) :: poly_rdi    !! A table that gather the polynimial parameters of 
                                                                    !! the 4 management style X Upper/lower boudaries
                                                                    !! of an RDI/quadratic diameter relationship

!$OMP THREADPRIVATE(poly_rdi)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: rdi_dia_ref       !!RDI at a reference tree diameter used to shape the end of
                                                                    !!the RDI/diameter relationship


!$OMP THREADPRIVATE(rdi_dia_ref)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: rdi_start         !!The rdi value at the start of a forest rotation, use
                                                                    !!0.4 to obtain a reallistic growth rate during 
                                                                    !!the early stage of the forest )

!$OMP THREADPRIVATE(rdi_start)  
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:) :: rdi_max         !!The maximum RDI that a forest will reach during
                                                                    !!its rotation
!$OMP THREADPRIVATE(rdi_max)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:) :: final_dia       !!Maximal tree diameter (m). If this diameter is exceeded a
                                                                    !!cut will happen
!$OMP THREADPRIVATE(final_dia)
 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: dia_boost_end     !!Shape parameter of the unmanaged forest circular
                                                                    !!distance distributi
!$OMP THREADPRIVATE(dia_boost_end)  
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: rdi_ramp_out  !!RDI vise en SORTIE de la classe d'amelioration, PAR CRENEAU (-)
!$OMP THREADPRIVATE(rdi_ramp_out)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: forced_thin_intensity  !! Fraction de RDI retiree chaque annee par l'eclaircie forcee, PAR CRENEAU (-)
!$OMP THREADPRIVATE(forced_thin_intensity)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: rdi_band_low  !!Borne BASSE de la bande de RDI cible, PAR CRENEAU (-)
!$OMP THREADPRIVATE(rdi_band_low)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: rdi_band_high !!Borne HAUTE de la bande de RDI cible, PAR CRENEAU (-)
!$OMP THREADPRIVATE(rdi_band_high)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: thin_intensity    !!Intensity of a thinnig event when a forest reach the upper
                                                                    !!rdi limit 
!$OMP THREADPRIVATE(thin_intensity)  
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: circ_dist_shape_none !! Shape parameter for the distribution of the 
                                                                   !! distance between trees in unmanaged forests 
!$OMP THREADPRIVATE(circ_dist_shape_none) 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: circ_dist_shape_thin   !! Shape parameter for the distribution of the 
                                                                    !! diameter class in managed forests
!$OMP THREADPRIVATE(circ_dist_shape_thin)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: circ_dist_shape_cop !! Shape parameter for the distribution of the 
                                                                    !! diameter class in coppice forests
!$OMP THREADPRIVATE(circ_dist_shape_cop)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: circ_dist_shape_src !! Shape parameter for the distribution of the 
                                                                    !! diameter class in short rotation coppice forests
!$OMP THREADPRIVATE(circ_dist_shape_src)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: circ_dist_shape_uneven !! Shape parameter for the distribution of the 
  !! diameter class in uneven-age forests
!$OMP THREADPRIVATE(circ_dist_shape_uneven) 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:) :: dens_target     !! The minimum density of trees in a stand before
                                                                    !! they all die off and we replant.  This is to prevent
                                                                    !! the stand from becoming just one large tree.
                                                                    !! @tex $(trees ha{-1})$ @endtex
!$OMP THREADPRIVATE(dens_target)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: forest_managed_forced !! Prescribed management system. This will only be used if
                                                                    !! no management is read
!$OMP THREADPRIVATE(forest_managed_forced)  
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: largest_tree_dia  !! The diameter at which we decide to clearcut
                                                                    !! a stand because our equipment cannot handle
                                                                    !! trees larger than this.
                                                                    !! @tex $(cm)$ @endtex
!$OMP THREADPRIVATE(largest_tree_dia)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: rotation_ref      !! Rotation de reference PAR PFT (an), gestion multifonctionnelle.
                                                                    !! Defaut par MTC (rotation_ref_mtc) ; surchargeable par ROTATION_REF.
                                                                    !! <= 0 => pas de rotation definie pour ce PFT.
!$OMP THREADPRIVATE(rotation_ref)
  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:) :: fm_src_pft            !! Force forest_managed = ifm_src (taillis courte rotation) sur ces PFT.
                                                                    !! Defaut FALSE partout. Pilote PAR PFT et non en dur sur "eucalyptus" :
                                                                    !! il n'existe pas de MTC eucalyptus, l'espece n'est qu'une affectation
                                                                    !! PFT_TO_MTC assortie de surcharges de namelist.
!$OMP THREADPRIVATE(fm_src_pft)
  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:) :: fm_cop_pft            !! Force forest_managed = ifm_cop (taillis simple) sur ces PFT.
                                                                    !! Defaut FALSE partout. Meme patron que fm_src_pft, autre regime :
                                                                    !! ifm_cop est la rotation longue (chene vert, 35-60 ans), ifm_src la
                                                                    !! courte (eucalyptus, 10-14 ans). Exclusifs l'un de l'autre.
!$OMP THREADPRIVATE(fm_cop_pft)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: fuelwood_diameter !! Diameter below which the wood harvest is used as fuelwood (m)
                                                                    !! Affects the way the wood is used in the dim_product_use
                                                                    !! subroutine         
!$OMP THREADPRIVATE(fuelwood_diameter)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: coppice_kill_be_wood !! Diameter below which the wood harvest is used as fuelwood (m)
!$OMP THREADPRIVATE(coppice_kill_be_wood)
                                 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: branch_ratio      !! branches/total aboveground biomass ratio
                                                                    !! (cf carbofor for CITEPA inventory, these 
                                                                    !! Guerric, Lim 2004, Peischl 2007, 
                                                                    !! Schulp 2008: 15-30% slash after harvest,
                                                                    !! Zaehle 2007: 30% slash after harvest)
!$OMP THREADPRIVATE(branch_ratio)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:) :: branch_harvest  !! The fraction of branches which are harvested
                                                                    !! during thinning and clearcut operations on 
                                                                    !! forests.  1.0 means all branches are taken offsite,
                                                                    !! 0.0 means all branches are left onsite and go
                                                                    !! into the litter pool.  This number is not
                                                                    !! based on any data.
!$OMP THREADPRIVATE(branch_harvest)               
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:) :: clearcut_residue!! Fraction of residues left on side after a thinning or harvest
!$OMP THREADPRIVATE(clearcut_residue) 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:) :: salvage_dist    !! The fraction of dead wood that is remove after disturbance (the
                                                                    !! rest are left onsite)

!$OMP THREADPRIVATE(salvage_dist) 
    

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)    :: coppice_diameter !! The trunk diameter at which a coppice will be cut.
                                                                    !! @tex $(m)$ @endtex
!$OMP THREADPRIVATE(coppice_diameter)  
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:) :: shoots_per_stool !! The number of shoots that regrow per stool after
                                                                    !! the first coppice cut
                                                                    !! @tex $-$ @endtex
!$OMP THREADPRIVATE(shoots_per_stool) 
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:) :: src_rot_length !! The number of years between cuttings for short
                                                                    !! rotation coppices.
                                                                    !! @tex $-$ @endtex
!$OMP THREADPRIVATE(src_rot_length) 
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:) :: src_nrots      !! The number of rotations for short rotations coppices
                                                                    !! before the roots are killed and replanted.
                                                                    !! @tex $-$ @endtex
!$OMP THREADPRIVATE(src_nrots) 

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: fruit_alloc         !! Fraction of biomass allocated to fruit production (0-1)

!$OMP THREADPRIVATE(fruit_alloc) 

 REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: labile_reserve       !! The size of the labile pool as a fraction of the
                                                                      !! weekly gpp (-). For example, 3 indicates that the 
                                                                      !! is 3 times the weekly gpp. 
!$OMP THREADPRIVATE(labile_reserve)

 REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: deciduous_reserve    !! Fraction of sapwood mass stored in the reserve pool of deciduous 
                                                                      !! trees during the growing season (unitless, 0-1) 

!$OMP THREADPRIVATE(deciduous_reserve)

 REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: evergreen_reserve    !! Fraction of sapwood mass stored in the reserve pool of evergreen 
                                                                      !! trees (unitless, 0-1) 

!$OMP THREADPRIVATE(evergreen_reserve)

 REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: senescense_reserve   !! Fraction of sapwood mass stored in the reserve pool of deciduous 
                                                                      !! trees during senescense(unitless, 0-1)

!$OMP THREADPRIVATE(senescense_reserve)

 REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: root_reserve         !! Fraction of sapwood mass stored in the reserve pool of deciduous 
                                                                      !! trees during the growing season (unitless, 0-1) 
!$OMP THREADPRIVATE(root_reserve)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: fcn_wood            !! CN ratio of wood for allocation, relative to leaf CN according
                                                                      !! to Sitch et al 2003 (https://doi.org/10.1046/j.1365-2486.2003.00569.x)
!$OMP THREADPRIVATE(fcn_wood)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: fcn_root            !! CN of roots for allocation, relative to leaf CN according
                                                                      !! to Sitch et al 2003 (https://doi.org/10.1046/j.1365-2486.2003.00569.x)
!$OMP THREADPRIVATE(fcn_root)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: cn_leaf_init        !! CN of foliage for allocation, according to Sitch et al 2003 
                                                                      !! (https://doi.org/10.1046/j.1365-2486.2003.00569.x)
!$OMP THREADPRIVATE(cn_leaf_init)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:) :: maint_resp_slope  !! slope of maintenance respiration coefficient  
                                                                      !! (1/K, 1/K^2, 1/K^3), used in the code 
!$OMP THREADPRIVATE(maint_resp_slope) 
 	 		 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: maint_resp_slope_c  !! slope of maintenance respiration coefficient (1/K), 
                                                                      !! constant c of aT^2+bT+c , tabulated 
!$OMP THREADPRIVATE(maint_resp_slope_c) 
  
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: maint_resp_slope_b  !! slope of maintenance respiration coefficient (1/K),  
                                                                      !! constant b of aT^2+bT+c , tabulated 
!$OMP THREADPRIVATE(maint_resp_slope_b) 
 	 	 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: maint_resp_slope_a  !! slope of maintenance respiration coefficient (1/K),  
                                                                      !! constant a of aT^2+bT+c , tabulated 
!$OMP THREADPRIVATE(maint_resp_slope_a)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:) :: moss_frac           !! Moss fraction by PFT to compute thermal capacity and conductivity 
!$OMP THREADPRIVATE(moss_frac)

  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:) :: is_peat                 !! flag by PFT for peatland vegetation types (true/false)
!$OMP THREADPRIVATE(is_peat)

  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:) :: is_croppeat             !! flag by PFT for peatland with agriculture (true/false)
!$OMP THREADPRIVATE(is_croppeat)

END MODULE pft_parameters_var
