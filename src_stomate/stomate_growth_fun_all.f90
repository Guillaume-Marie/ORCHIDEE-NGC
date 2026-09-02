!=================================================================================================================================
! MODULE       : stomate_allocation
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
!                This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF         Plant growth and C-allocation among the biomass components (leaves, wood, roots, fruit, reserves, labile) 
!!               is calculated making use of functional allocation which combines the pipe model and allometric relationships 
!!               proposed by Sitch et al 2003 and adjusted by Zaehle et al 2010.
!!
!!\n DESCRIPTION: This module calculates three processes: (1) daily maintenance respiration based on the half-hourly
!!               respiration calculated in stomate_resp.f90, (2) the absolute allocation to the different biomass
!!               components based on functional allocation approach and (3) the allocatable biomass as the residual
!!               of GPP-Ra. Multiplication of the allocation fractions and allocatable biomass given the changes in
!!               biomass pools.
!!
!! RECENT CHANGE(S): Until 1.9.6 only one allocation scheme was available (now contained in stomate_grwoth_res_lim.f90). 
!!               This module consists of an alternative formalization of plant growth.
!!                             
!! REFERENCE(S)	: - Sitch, S., Smith, B., Prentice, I.C., Arneth, A., Bondeau, A., Cramer, W.,
!!               Kaplan, J.O., Levis, S., Lucht, W., Sykes, M.T., Thonicke, K., Venevsky, S. (2003), Evaluation of 
!!               ecosystem dynamics, plant geography and terrestrial carbon cycling in the LPJ Dynamic Global Vegetation 
!!               Model, Global Change Biology, 9, 161-185.\n 
!!               - Zaehle, S. and Friend, A.D. (2010), Carbon and nitrogen cycle dynamics in the O-CN land surface model: 1. 
!!               Model description, site-scale evaluation, and sensitivity to parameter estimates, Global Biogeochemical 
!!               Cycles, 24, GB1005.\n
!!               - Magnani F., Mencuccini M. & Grace J. 2000. Age-related decline in stand productivity: the role of 
!!               structural acclimation under hydraulic constraints Plant, Cell and Environment 23, 251–263.\n
!!               - Bloom A.J., Chapin F.S. & Mooney H.A. (1985) Resource limitation in plants. An economic analogy.  
!!               Annual Review Ecology Systematics 16, 363–392.\n
!!               - Case K.E. & Fair R.C. (1989) Principles of Economics. Prentice Hall, London.\n
!!               - McDowell, N., Barnard, H., Bond, B.J., Hinckley, T., Hubbard, R.M., Ishii, H., Köstner, B., 
!!               Magnani, F. Marshall, J.D., Meinzer, F.C., Phillips, N., Ryan, M.G., Whitehead D. 2002. The 
!!               relationship between tree height and leaf area: sapwood area ratio. Oecologia, 132:12–20.\n
!!               - Novick, K., Oren, R., Stoy, P., Juang, F.-Y., Siqueira, M., Katul, G. 2009. The relationship between
!!               reference canopy conductance and simplified hydraulic architecture. Advances in water resources 32,
!!               809-819.     
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_stomate/stomate_growth_fun_all.f90 $
!! $Date: 2026-03-03 14:02:01 +0100 (mar. 03 mars 2026) $
!! $Revision: 9417 $
!! \n
!_ ===============================================================================================================================

MODULE stomate_growth_fun_all

  ! Modules used:
  USE ioipsl_para
  USE xios_orchidee
  USE pft_parameters
  USE dynamic_parameters
  USE stomate_data
  USE constantes
  USE constantes_soil
  USE function_library,    ONLY: wood_to_qmdia, wood_to_qmheight, &
                           wood_to_ba_eff, cc_to_lai, lai_to_biomass, &
                           biomass_to_lai, cc_to_biomass, biomass_to_cc, &
                           calculate_c0_alloc, wood_to_ba, &
                           check_vegetation_area, check_mass_balance, &
                           wood_to_height, intermediate_mass_balance_check, &
                           get_printlev, wood_to_dia

  IMPLICIT NONE

  ! Private & public routines

  PRIVATE
  PUBLIC growth_fun_all_clear, growth_fun_all

 ! Variables shared by all subroutines in this module

  LOGICAL, SAVE                                             :: firstcall_growth_fun_all = .TRUE.  !! Is this the first call? (true/false)
!$OMP THREADPRIVATE(firstcall_growth_fun_all)

  !+++TEMP+++
  INTEGER, SAVE                                             :: istep = 0
!$OMP THREADPRIVATE(istep)
  INTEGER(i_std), SAVE       :: printlev_loc       !! Local level of text output for current module
!$OMP THREADPRIVATE(printlev_loc)

CONTAINS


!! ================================================================================================================================
!! SUBROUTINE   : growth_fun_all_clear
!!
!>\BRIEF          Set the flag ::firstcall to .TRUE. and as such activate section 
!! 1.1 of the subroutine alloc (see below).\n
!!
!_ ================================================================================================================================

  SUBROUTINE growth_fun_all_clear
    firstcall_growth_fun_all = .TRUE.
  END SUBROUTINE growth_fun_all_clear



!! ================================================================================================================================
!! SUBROUTINE 	: growth_fun_all
!!
!>\BRIEF          Allocate net primary production (= photosynthesis
!!                minus autothrophic respiration) to: labile carbon pool carbon reserves, aboveground sapwood,
!!                belowground sapwood, root, fruits and leaves following the pipe model and allometric constraints.
!!
!! DESCRIPTION  : Total maintenance respiration for the whole plant is calculated by summing maintenance 
!!                respiration of the different plant compartments. Maintenance respiration is subtracted 
!!                from whole-plant allocatable biomass (up to a maximum fraction of the total allocatable biomass). 
!!                Growth respiration is then calculated as a prescribed (0.75) fraction of the allocatable
!!                biomass. Subsequently NPP is calculated by substracting total autotrophic  respiration from GPP 
!!                i.e. NPP = GPP - maintenance resp - growth resp.
!!
!!                The pipe model assumes that one unit of leaf mass requires a proportional amount of sapwood to 
!!                transport water from the roots to the leaves. Also a proportional fraction of roots is needed to 
!!                take up the water from the soil. The proportional amounts between leaves, sapwood and roots are 
!!                given by allocation factors. These allocation factors are PFT specific and depend on a parameter 
!!                quantifying the leaf to sapwood area (::k_latosa_target), the specific leaf area (::sla), wood 
!!                density (::pipe_density) and a scaling parameter between leaf and root mass.
!! 
!!                Lai is optimised for mean annual radiation use efficiency and the C cost for producing the 
!!                canopy. The cost-benefit ratio is optimised when the marginal gain / marginal cost = 1 lai target 
!!                is used to calculate whether the reserves are used. This approach allows plants to get out of 
!!                senescence and to start developping a canopy in early spring. 
!!                  
!!                As soon as a canopy has emerged, C (b_inc_tot) becomes available at the stand level through 
!!                photosynthesis and, C is allocated at the tree level (b_inc) following both the pipe model and 
!!                allometric constraints. Mass conservation requires:
!!                (1) Cs_inc + Cr_inc + Cl_inc = b_inc
!!                (2) sum(b_inc) = b_inc_tot
!!
!!                Wood allocation depends on tree basal area following the rule of Deleuze & Dhote
!!                delta_ba = gammas*(circ - m*sigmas + sqrt((m*sigmas + circ).^2 - (4*sigmas*circ)))/2
!!                (3) <=> delta_ba = circ_class_dba*gammas
!!                Where circ_class_dba = (circ - m*sigmas + sqrt((m*sigmas + circ).^2 - (4*sigmas*circ)))/2
!! 
!!                Allometric relationships
!!                height = pipe_tune2*(dia.^pipe_tune3)
!!                Re-write this relationship as a function of ba 
!!                (4) height = pipe_tune2 * (4/pi*ba)^(pipe_tune3/2)
!!                (5a) Cl/Cs = KF/height for trees
!!                (5b) Cs = Cl / KF
!!                (6) Cl = Cr * LF
!!
!!                Use a linear approximation to avoid iterations. Given that allocation is calculated daily, a 
!!                local lineair assumption is fair. Eq (4) can thus be rewritten as:
!!                s = step/(pipe_tune2*(4/pi*(ba+step)).^(pipe_tune3/2)-pipe_tune2*(4/pi*ba).^(pipe_tune3/2))
!!                Where step is a small but realistic (for the time step) change in ba
!!                (7)  <=> delta_height = delta_ba/s
!!
!!                Calculate Cs_inc from allometric relationships
!!                Cs_inc = tree_ff*pipe_density*(ba+delta_ba)*(height+delta_height) - Cs - Ch         
!!                Cs_inc = tree_ff*pipe_density*(ba+delta_ba)*(height+delta_ba/s) - Cs - Ch
!!                (8)  <=> Cs_inc = tree_ff*pipe_density*(ba+a*gammas)*(height+(a/s*gammas)) - Cs - Ch
!!
!!                Rewrite (5) as 
!!                Cl_inc = KF*(Cs_inc+Cs)/(height+delta_height) - Cl
!!                Substitute (7) in (4) and solve for Cl_inc
!!                Cl_inc = KF*(tree_ff*pipe_density*(ba+circ_class_dba*gammas)*(height+(circ_class_dba/s*gammas)) - Ch)/ & 
!!                   (height+(circ_class_dba/s*gammas)) - Cl  
!!                (9)  <=> Cl_inc = KF*tree_ff*pipe_density*(ba+circ_class_dba*gammas) - &
!!                            (KF*Ch)/(height+(circ_class_dba/s*gammas)) - Cl
!!
!!                Rewrite (6) as 
!!                Cr_inc = (Cl_inc+Cl)/LF - Cr
!!                Substitute (9) in (6)
!!                (10)  <=> Cr_inc = KF/LF*tree_ff*pipe_density*(ba+circ_class_dba*gammas) - &
!!                            (KF*Ch/LF)/(height+(circ_class_dba/s*gammas)) - Cr
!!
!!                Depending on the specific case that needs to be solved equations (1) takes one of the following forms:
!!                (a) b_inc = Cl_inc + Cr_inc + Cs_inc, (b) b_inc = Cl_inc + Cr_inc, (c) b_inc = Cl_inc + Cs_inc or 
!!                (d) b_inc = Cr_inc + Cs_inc. One of these alternative forms of eq. 1 are then combined with 
!!                eqs 8, 9 and 10 and solved for gammas. The details for the solution of these four cases are given in the
!!                code. Once gammas is know, eqs 6 - 10 are used to calculate the allocation to leaves (Cl_inc), 
!!                roots (Cr_inc) and sapwood (Cs_inc).
!!
!!                Because of turnover, biomass pools are not all the time in balance according to rules prescribed
!!                by the pipe model. To test whether biomass pools are balanced, the target biomasses are calculated 
!!                and balance is restored whenever needed up to the level that the biomass pools for leaves, sapwood
!!                and roots are balanced according to the pipe model. Once the balance is restored C is allocated to
!!                fruits, leaves, sapwood and roots by making use of the pipe model (below this called ordinary 
!!                allocation).
!!
!!                Although strictly speaking allocation factors are not necessary in this scheme (Cl_inc could simply 
!!                be added to biomass(:,:,ileaf,icarbon), Cr_inc to biomass(:,:,iroot,icarbon), etc.), they are 
!!                nevertheless calculated because using allocation factors facilitates comparison to the resource 
!!                limited allocation scheme (stomate_growth_res_lim.f90) and it comes in handy for model-data comparison.
!!
!!                Effective basal area, height and circumferences are use in the allocation scheme because their 
!!                calculations make use of the total (above and belowground) biomass. In forestry the same measures
!!                exist (and they are also calculated in ORCHIDEE) but only account for the aboverground biomass.
!! 
!!
!! RECENT CHANGE(S): - The code by Sonke Zaehle made use of ::Cl_target that was derived from ::lai_target which in turn 
!!                was a function of ::rue_longterm. Cl_target was then used as a threshold value to decide whether there 
!!                was only phenological growth (just leaves and roots) or whether there was full allometric growth to the 
!!                leaves, roots and sapwood. This approach was inconsistent with the pipe model because full allometric 
!!                growth can only occur if all three biomass pools are in balance. ::lai_target is no longer used as a 
!!                criterion to switch between phenological and full allometric growth. Its use is now restricted to trigger
!!                the use of reserves in spring. 
!!
!!                Early 2019: Nitrogen limitations to growth were added, primarily based on the ratio of carbon
!!                            to nitrogen in the leaf.
!!
!! MAIN OUTPUT VARIABLE(S): ::npp and :: biomass. Seven different biomass compartments (leaves, roots, above and 
!!                belowground wood, carbohydrate reserves, labile and fruits).
!!
!! REFERENCE(S)	:- Sitch, S., Smith, B., Prentice, I.C., Arneth, A., Bondeau, A., Cramer, W.,
!!                Kaplan, J.O., Levis, S., Lucht, W., Sykes, M.T., Thonicke, K., Venevsky, S. (2003), Evaluation of 
!!                ecosystem dynamics, plant geography and terrestrial carbon cycling in the LPJ Dynamic Global Vegetation 
!!                Model, Global Change Biology, 9, 161-185.\n 
!!                - Zaehle, S. and Friend, A.D. (2010), Carbon and nitrogen cycle dynamics in the O-CN land surface model: 1. 
!!                Model description, site-scale evaluation, and sensitivity to parameter estimates, Global Biogeochemical 
!!                Cycles, 24, GB1005.\n
!!                - Magnani F., Mencuccini M. & Grace J. 2000. Age-related decline in stand productivity: the role of 
!!                structural acclimation under hydraulic constraints Plant, Cell and Environment 23, 251–263.
!!                - Bloom A.J., Chapin F.S. & Mooney H.A. (1985) Resource limitation in plants. An economic analogy. Annual 
!!                Review Ecology Systematics 16, 363–392.
!!                - Case K.E. & Fair R.C. (1989) Principles of Economics. Prentice Hall, London.
!!                - McDowell, N., Barnard, H., Bond, B.J., Hinckley, T., Hubbard, R.M., Ishii, H., Köstner, B., 
!!                Magnani, F. Marshall, J.D., Meinzer, F.C., Phillips, N., Ryan, M.G., Whitehead D. 2002. The 
!!                relationship between tree height and leaf area: sapwood area ratio. Oecologia, 132:12–20
!!                - Novick, K., Oren, R., Stoy, P., Juang, F.-Y., Siqueira, M., Katul, G. 2009. The relationship between
!!                reference canopy conductance and simplified hydraulic architecture. Advances in water resources 32,
!!                809-819.
!!                - Jefferey Amthor. 2000. The McCree-de Wit-Penning de Vries-Thornley Respiration Paradigms: 30 Years Later. 
!!                Annals of Botany 86: 1-20, doi:10.1006/anbo.2000.1175
!!                - JEFFREY Q. CHAMBERS, EDGARD S. TRIBUZY, LIGIA C. TOLEDO, BIANCA F. CRISPIM, NIRO HIGUCHI, JOAQUIM DOS SANTOS, 
!!                ALESSANDRO C. ARAU´ JO, BART KRUIJT, ANTONIO D. NOBRE, AND SUSAN E. TRUMBORE. 2004. RESPIRATION FROM A TROPICAL 
!!                FOREST ECOSYSTEM: PARTITIONING OF SOURCES AND LOW CARBON USE EFFICIENCY. Ecological Applications, 14(4) 
!!                Supplement, 2004, pp. S72S88
!!                -  Teemu Hölttä, Anna Lintunen, Tommy Chan, Annikki Mäkelä and Eero Nikinmaa. 2017. A steady-state stomatal 
!!                model of balanced leaf gas exchange, hydraulics and maximal sourcesink flux. Tree Physiology 37, 851-868. 
!!                doi:10.1093/treephys/tpx011
!!
!! FLOWCHART    : 
!!
!_ ================================================================================================================================

  SUBROUTINE growth_fun_all (npts, dt, veget_max, veget, &
       PFTpresent, plant_status ,when_growthinit, t2m, &
       nstress_season, vegstress_season, &
       gpp_daily, gpp_week, resp_maint_part, resp_maint, &
       resp_growth, npp, bm_alloc, age, &
       leaf_age, leaf_frac, &
       rue_longterm, circ_class_n, &
       circ_class_biomass, KF, sigma, &
       gammas, longevity_eff_leaf, longevity_eff_sap, &
       longevity_eff_root, k_latosa_adapt, forest_managed, &
       cn_leaf_min_season, atm_to_bm,  &
       cn_leaf_min_2D, cn_leaf_max_2D, sugar_load, &
       n_reserve_balance, n_reserve_longterm)       

 !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                        :: npts                   !! Domain size - number of grid cells 
                                                                                !! (unitless)
    REAL(r_std), INTENT(in)                           :: dt                     !! Time step of the simulations for stomate
                                                                                !! (days)
    REAL(r_std), DIMENSION(:), INTENT(in)             :: t2m                    !! Temperature at 2 meter (K)   
    REAL(r_std), DIMENSION(:,:), INTENT(in)           :: veget_max              !! PFT "Maximal" coverage fraction of a PFT 
                                                                                !! (= ind*cn_ind) 
                                                                                !! @tex $(m^2 m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(in)           :: veget                  !! Fraction of forest floor covered by vegetation (unitless, 0-1)   
    REAL(r_std), DIMENSION(:,:), INTENT(in)           :: when_growthinit        !! Days since beginning of growing season 
                                                                                !! (days)
    REAL(r_std), DIMENSION(:,:), INTENT(in)           :: rue_longterm           !! Longterm "radiation use efficicency"
                                                                                !! calculated as the ratio of GPP over 
                                                                                !! the fraction of radiation absorbed
                                                                                !! by the canopy 
                                                                                !! @tex $(gC.m^{-2}day^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(in)           :: longevity_eff_root     !! Effective root turnover time that accounts
                                                                                !! waterstress (days)
    REAL(r_std), DIMENSION(:,:), INTENT(in)           :: longevity_eff_sap      !! Effective sapwood turnover time that accounts
                                                                                !! waterstress (days)
    REAL(r_std), DIMENSION(:,:), INTENT(in)           :: longevity_eff_leaf     !! Effective leaf turnover time that accounts
                                                                                !! waterstress (days)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)      :: circ_class_n           !! Number of individuals in each circ class
                                                                                !! @tex $(m^{-2})$ @endtex
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)         :: resp_maint_part        !! Maintenance respiration of different  
                                                                                !! plant parts
                                                                                !! @tex $(gC.m^{-2}dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(in)           :: plant_status           !! Growth and phenological status of the plant

    LOGICAL, DIMENSION(:,:), INTENT(in)               :: PFTpresent             !! PFT exists (true/false)    
    REAL(r_std), DIMENSION(:,:), INTENT(in)           :: nstress_season         !! N-related seasonal stress (used for allocation)    
    REAL(r_std), DIMENSION(:,:), INTENT(in)           :: vegstress_season    !! mean growing season moisture availability
    INTEGER(i_std), DIMENSION(:,:), INTENT(in)        :: forest_managed         !! Forest management flag: 0 = orchidee 
                                                                                !! standard, 1= self-thinning only, 2= 
                                                                                !! high-stand, 3= high-stand smoothed, 4= 
                                                                                !! coppices
    REAL(r_std), DIMENSION(:,:), INTENT(in)           :: gpp_week               !! PFT gross primary productivity 
    REAL(r_std),DIMENSION(npts,nvm), INTENT(in)       :: cn_leaf_min_2D         !! minimal leaf C/N ratio 
    REAL(r_std),DIMENSION(npts,nvm), INTENT(in)       :: cn_leaf_max_2D         !! maximal leaf C/N ratio 
    REAL(r_std), DIMENSION(:,:), INTENT(in)           :: cn_leaf_min_season     !! Min leaf nitrogen concentration (C:N) of the growing season
                                                                                !! (gC/gN)  
    REAL(r_std), DIMENSION(:,:), INTENT(in)           :: n_reserve_longterm     !! "longer term" actual to potential  N reserve pool
                                                                                !! (0-1, unitless) 
   
    !! 0.2 Output variables

    REAL(r_std), DIMENSION(:,:), INTENT(out)          :: resp_maint             !! PFT maintenance respiration 
                                                                                !! @tex $(gC.m^{-2}dt^{-1})$ @endtex    
    REAL(r_std), DIMENSION(:,:), INTENT(out)          :: resp_growth            !! PFT growth respiration 
                                                                                !! @tex $(gC.m^{-2}dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(out)          :: npp                    !! PFT net primary productivity 
                                                                                !! @tex $(gC.m^{-2}dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(out)      :: bm_alloc               !! PFT biomass increase, i.e. NPP per plant  
                                                                                !! part @tex $(gC.m^{-2}dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(out)          :: sigma                  !! Threshold for indivudal tree growth in 
                                                                                !! the equation of Deleuze & Dhote (2004)(m).
                                                                                !! Trees whose circumference is smaller than 
                                                                                !! sigma don't grow much
    REAL(r_std), DIMENSION(:,:), INTENT(out)          :: gammas                 !! Slope for individual tree growth in the 
                                                                                !! equation of Deleuze & Dhote (2004) (m)

    REAL(r_std), DIMENSION(:,:), INTENT(out)          :: sugar_load             !! Relative sugar loading of the labile pool (unitless)
    REAL(r_std), DIMENSION(:,:), INTENT(out)          :: n_reserve_balance      !! Actual to potential N reserve pool (unitless)

    !! 0.3 Modified variables

    REAL(r_std), DIMENSION(:,:), INTENT(inout)        :: gpp_daily              !! PFT gross primary productivity 
                                                                                !! @tex $(gC.m^{-2}dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)        :: age                    !! PFT age (days)      
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)      :: leaf_age               !! PFT age of different leaf classes 
                                                                                !! (days)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)      :: leaf_frac              !! PFT fraction of leaves in leaf age class
                                                                                !! (0-1, unitless)
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(inout)  :: circ_class_biomass     !! Biomass components of the model tree  
                                                                                !! within a circumference class
                                                                                !! class @tex $(g C ind^{-1})$ @endtex
    REAL(r_std), DIMENSION(:,:), INTENT(inout)        :: KF                     !! Scaling factor to convert sapwood mass
                                                                                !! into leaf mass (m)
    REAL(r_std), DIMENSION(:,:), INTENT(inout)        :: k_latosa_adapt         !! Leaf to sapwood area adapted for long 
                                                                                !! term water stress (m)
    REAL(r_std), DIMENSION(:,:,:), INTENT(inout)      :: atm_to_bm              !! Nitrogen and carbon which is added to the ecosystem to 
                                                                                !! support vegetation growth (gC or gN/m2/day)
    
    !! 0.4 Local variables

    CHARACTER(30)                                     :: var_name               !! To store variable names for I/O
    REAL(r_std), DIMENSION(npts,nvm)                  :: c0_alloc               !! Root to sapwood tradeoff parameter
    LOGICAL                                           :: grow_wood              !! Flag to grow wood
    INTEGER(i_std)                                    :: ipts,j,k,l,m           !! Indices(unitless)
    INTEGER(i_std)                                    :: icir,imed,ipool        !! Indices(unitless)
    INTEGER(i_std)                                    :: ifm,icut               !! Indices
    INTEGER(i_std)                                    :: ipar,iele,imbc         !! Indices(unitless)
    INTEGER(i_std)                                    :: ilev                   !! Indices(unitless)
    REAL(r_std)                                       :: frac                   !! No idea??
    REAL(r_std)                                       :: a,b,c                  !! Temporary variables to solve a
                                                                                !! quadratic equation (unitless)
    REAL(r_std),DIMENSION(npts,nvm)                   :: gtemp                  !! Turnover coefficient of labile C pool
                                                                                !! (0-1)
    REAL(r_std),DIMENSION(npts,nvm,nelements)         :: reserve_target         !! Intentional size of the reserve pool
                                                                                !! @tex $(gC/N.m^{-2})$ @endtex
    REAL(r_std),DIMENSION(npts,nvm,nelements)         :: labile_target          !! Intentional size of the labile pool
                                                                                !! @tex $(gC/N.m^{-2})$ @endtex
    REAL(r_std)                                       :: reserve_scal           !! Protection of the reserve against 
                                                                                !! overuse (unitless)
    REAL(r_std)                                       :: use_lab                !! Availability of labile biomass 
                                                                                !! @tex $(gC.m^{-2})$ @endtex
    REAL(r_std)                                       :: use_res                !! Availability of resource biomass 
                                                                                !! @tex $(gC.m^{-2})$ @endtex
    REAL(r_std)                                       :: use_max                !! Maximum use of labile and resource pool 
                                                                                !! @tex $(gC.m^{-2})$ @endtex
    REAL(r_std)                                       :: leaf_meanage           !! Mean age of the leaves (days?)
    REAL(r_std)                                       :: reserve_time           !! Maximum number of days during which 
                                                                                !! carbohydrate reserve may be used (days)
    REAL(r_std)                                       :: b_inc_tot              !! Carbon that needs to allocated in the
                                                                                !! fixed number of trees (gC)
    REAL(r_std)                                       :: b_inc_temp             !! Temporary b_inc at the stand-level
                                                                                !! @tex $(gC.plant^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)                  :: scal                   !! Scaling factor between average 
                                                                                !! individual and individual plant
                                                                                !! @tex $(plant.m^{-2})$ @endtex
    REAL(r_std)                                       :: total_inc              !! Total biomass increase
                                                                                !! @tex $(gC.plant^{-1})$ @endtex
    REAL(r_std)                                       :: KF_old                 !! Scaling factor to convert sapwood mass
                                                                                !! into leaf mass (m) at the previous
                                                                                !! time step
    REAL(r_std)                                       :: sla_est                !! A first estimate of sla in case its calculation is 
                                                                                !! dynamic @tex $(m^2.gC^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)                  :: sla_out                !! variable to write the value of SLA to the history file 
    REAL(r_std), DIMENSION(nvm)                       :: lai_happy              !! Lai threshold below which carbohydrate 
                                                                                !! reserve may be used 
                                                                                !! @tex $(m^2 m^{-2})$ @endtex
    REAL(r_std), DIMENSION(nvm)                       :: deleuze_p              !! Percentile of trees that will receive 
                                                                                !! photosynthates. The proxy for intra stand
                                                                                !! competition. Depends on the management
                                                                                !! strategy when ncirc < 6
    REAL(r_std), DIMENSION(npts)                      :: tl                     !! Long term annual mean temperature (C)
    REAL(r_std), DIMENSION(npts)                      :: bm_add                 !! Biomass increase 
		                                                                !! @tex $(gC.m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts)                      :: bm_new                 !! New biomass @tex $(gC.m^{-2})$ @endtex
    REAL(r_std)                                       :: alloc_sap_above        !! Fraction allocated to sapwood above 
                                                                                !! ground
    REAL(r_std), DIMENSION(npts,nvm)                  :: residual               !! Copy of b_inc_tot after all C has been
                                                                                !! allocated @tex $(gC.m^{-2})$ @endtex
                                                                                !! if all went well the value should be zero
    REAL(r_std), DIMENSION(npts,nvm)                  :: residual_write         !! Copy of b_inc_tot after all C has been
                                                                                !! allocated @tex $(gC.m^{-2})$ @endtex
                                                                                !! if all went well the value should be zero. This
                                                                                !! value is written to the history file to better
                                                                                !! monitor the residuals
    REAL(r_std), DIMENSION(npts,nvm)                  :: lai_target             !! Target LAI @tex $(m^{2}m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)                  :: ltor                   !! Leaf to root ratio (unitless)   
    REAL(r_std), DIMENSION(npts,nvm)                  :: lstress_fac            !! Light stress factor, based on total
                                                                                !! transmitted light (unitless, 0-1)
    REAL(r_std), DIMENSION(npts,nvm)                  :: LF                     !! Scaling factor to convert sapwood mass
                                                                                !! into root mass (unitless)
    REAL(r_std), DIMENSION(npts,nvm)                  :: lm_old                 !! Variable to store leaf biomass from 
                                                                                !! previous time step 
                                                                                !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)                  :: bm_alloc_tot           !! Allocatable biomass for the whole plant
                                                                                !! @tex $(gC.m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)                  :: temp_bm_alloc_tot      !! Allocatable biomass for the whole plant
                                                                                !! @tex $(gC.m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)                  :: resid_bm_alloc_tot     !! Allocatable biomass for the whole plant
                                                                                !! @tex $(gC.m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)                  :: leaf_mass_young        !! Leaf biomass in youngest leaf age class 
                                                                                !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)                  :: lai                    !! PFT leaf area index 
                                                                                !! @tex $(m^2 m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm)                  :: qm_dia                 !! Quadratic mean diameter of the stand (m)
    REAL(r_std), DIMENSION(npts,nvm)                  :: qm_height              !! Height of a tree with the quadratic mean
                                                                                !! diameter (m)
    REAL(r_std), DIMENSION(npts,nvm)                  :: ba                     !! Basal area. variable for histwrite only (m2)
    REAL(r_std), DIMENSION(npts,nvm,nparts)           :: f_alloc                !! PFT fraction of NPP that is allocated to
                                                                                !! the different components (0-1, unitless)
    REAL(r_std), DIMENSION(npts,ncirc,nparts)         :: f_alloc_circ           !! Fraction of that is allocated to each circc_class
                                                                                !! the different components (0-1, unitless)
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements) :: tmp_bm                 !! temporary variable to indicate biomass for each PFT
                                                                                !! over per unit PFT area.
                                                                                !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements) :: tmp_init_bm            !! temporary variable to use site-level biomass 
                                                                                !! @tex $(gC m^{-2})$ @endtex

    ! Tree level	
    REAL(r_std), DIMENSION(ncirc)                     :: step                   !! Temporary variables to solve a
                                                                                !! quadratic equation (unitless)
    REAL(r_std), DIMENSION(ncirc)                     :: s                      !! tree-level linear relationship between
                                                                                !! basal area and height. This variable is
                                                                                !! used to linearize the allocation scheme
    REAL(r_std), DIMENSION(ncirc)                     :: circ_class_n_safe      !! circ_class_n with empty classes set to
                                                                                !! one so the element-wise division below
                                                                                !! stays finite (@tex $(trees m^{-2})$ @endtex)
    REAL(r_std), DIMENSION(ncirc)                     :: Cs_inc_est             !! Initial value estimate for Cs_inc. The
                                                                                !! value is used to linearize the ba~height
                                                                                !! relationship 
                                                                                !! @tex $(gC.plant^{-1})$ @endtex
    REAL(r_std), DIMENSION(ncirc)                     :: Cl                     !! Individual plant, leaf compartment
                                                                                !! @tex $(gC.plant^{-1})$ @endtex
    REAL(r_std), DIMENSION(ncirc)                     :: Cr                     !! Individual plant, root compartment
                                                                                !! @tex $(gC.plant^{-1})$ @endtex
    REAL(r_std), DIMENSION(ncirc)                     :: Cs                     !! Individual plant, sapwood compartment
                                                                                !! @tex $(gC.plant^{-1})$ @endtex
    REAL(r_std), DIMENSION(ncirc)                     :: Ch                     !! Individual plant, heartwood compartment
                                                                                !! @tex $(gC.plant^{-1})$ @endtex
    REAL(r_std), DIMENSION(ncirc)                     :: Cl_inc                 !! Individual plant increase in leaf 
                                                                                !! compartment
                                                                                !! @tex $(gC.plant^{-1})$ @endtex 
    REAL(r_std), DIMENSION(ncirc)                     :: Cr_inc                 !! Individual plant increase in root
                                                                                !! compartment
                                                                                !! @tex $(gC.plant^{-1})$ @endtex 
    REAL(r_std), DIMENSION(ncirc)                     :: Cs_inc                 !! Individual plant increase in sapwood
                                                                                !! compartment
                                                                                !! @tex $(gC.plant^{-1})$ @endtex 
    REAL(r_std), DIMENSION(ncirc)                     :: Cf_inc                 !! Individual plant increase in fruit
                                                                                !! compartment
                                                                                !! @tex $(gC.plant^{-1})$ @endtex
    REAL(r_std), DIMENSION(ncirc)                     :: Cl_incp                !! Phenology related individual plant 
                                                                                !! increase in leaf compartment
                                                                                !! @tex $(gC.plant^{-1})$ @endtex
    REAL(r_std), DIMENSION(ncirc)                     :: Cr_incp                !! Phenology related individual plant 
                                                                                !! increase in leaf compartment
                                                                                !! @tex $(gC.plant^{-1})$ @endtex
    REAL(r_std), DIMENSION(ncirc)                     :: Cs_incp                !! Phenology related individual plant 
                                                                                !! increase in sapwood compartment
                                                                                !! @tex $(gC.plant^{-1})$ @endtex
    REAL(r_std), DIMENSION(ncirc)                     :: Cl_target              !! Individual plant maximum leaf mass given 
                                                                                !! its current sapwood mass
                                                                                !! @tex $(gC.plant^{-1})$ @endtex
    REAL(r_std), DIMENSION(ncirc)                     :: Cr_target              !! Individual plant maximum root mass given 
                                                                                !! its current sapwood mass
                                                                                !! @tex $(gC.plant^{-1})$ @endtex
    REAL(r_std), DIMENSION(ncirc)                     :: Cs_target              !! Individual plant maximum sapwood mass 
                                                                                !! given its current leaf mass or root mass
                                                                                !! @tex $(gC.plant^{-1})$ @endtex     
    REAL(r_std), DIMENSION(ncirc)                     :: delta_ba               !! Change in basal area for a unit
                                                                                !! investment into sapwood mass (m)
    REAL(r_std), DIMENSION(ncirc)                     :: delta_height           !! Change in height for a unit
                                                                                !! investment into sapwood mass (m)
    REAL(r_std), DIMENSION(ncirc)                     :: circ_class_ba          !! Basal area (forestry definition) of the model 
                                                                                !! tree in each circ class 
                                                                                !! @tex $(m^{2} m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(ncirc)                     :: circ_class_ba_eff      !! Effective basal area of the model tree in each
                                                                                !! circ class @tex $(m^{2} m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(ncirc)                     :: circ_class_dba         !! Share of an individual tree in delta_ba
                                                                                !! thus, circ_class_dba*gammas = delta_ba
    REAL(r_std), DIMENSION(ncirc)                     :: circ_class_height_eff  !! Effective tree height calculated from allometric 
                                                                                !! relationships (m)
    REAL(r_std), DIMENSION(ncirc)                     :: circ_class_circ_eff    !! Effective circumference of individual trees (m)
    REAL(r_std)                                       :: woody_biomass          !! Woody biomass. Temporary variable to 
                                                                                !! calculate wood volume (gC m-2)
    REAL(r_std), DIMENSION(ncirc)                     :: share_ncirc            !! Temporary variable to store the share
                                                                                !! of biomass of each circumference class
                                                                                !! to the total biomass
    REAL(r_std)                                       :: temp_share             !! Temporary variable to store the share
                                                                                !! of biomass of each circumference class
                                                                                !! to the total biomass        
    REAL(r_std)                                       :: temp_class_biomass     !! Biomass across parts for a single circ
                                                                                !! class @tex $(gC m^{-2})$ @endtex
    REAL(r_std)                                       :: temp_total_biomass     !! Biomass across parts and circ classes
                                                                                !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,ncirc)            :: store_delta_ba_eff     !! Store effective delta_ba in this variable before writing
                                                                                !! to the output file (m). Adding this variable
                                                                                !! was faster than changing the dimensions
                                                                                !! of delta_ba which would have been the same
    REAL(r_std), DIMENSION(npts,nvm,ncirc)            :: store_delta_ba
    REAL(r_std), DIMENSION(npts,nvm,ncirc)            :: store_circ_class_ba    !! Store circ_class_ba in this variable before
                                                                                !! writing to the output file (m). Adding this
                                                                                !! variable was faster than changing the 
                                                                                !! dimensions of circ_class_ba_ba which would 
                                                                                !! have been the same
    REAL(r_std), DIMENSION(npts,nvm,ncirc)            :: circ_class_ba_init     !! Basal area per diameter class before growth
    REAL(r_std), DIMENSION(npts,nvm,ncirc)            :: ring_width             !! Increase of radius of trunk. Calculated using
                                                                                !! store_delta_ba which is always positive
    REAL(r_std), DIMENSION(npts,nvm,ncirc)            :: circ_height            !! Height of trees per diameter class
    REAL(r_std), DIMENSION(npts,nvm,nmbcomp,nelements):: check_intern           !! Contains the components of the internal
                                                                                !! mass balance chech for this routine
                                                                                !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nmbcomp,nelements):: check_intern_init      !! Contains the components of the internal
                                                                                !! mass balance chech for this routine
                                                                                !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)        :: closure_intern         !! Check closure of internal mass balance
                                                                                !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm,nelements)        :: pool_start, pool_end   !! Start and end pool of this routine 
                                                                                !! @tex $(gC pixel^{-1} dt^{-1})$ @endtex
    REAL(r_std)                                       :: median_circ            !! Median circumference (m)
    REAL(r_std)                                       :: deficit                !! Carbon that needs to be respired in 
                                                                                !! excess of todays gpp  
                                                                                !! @tex $(gC.m^{-2}dt^{-1})$ @endtex
    REAL(r_std)                                       :: excess                 !! Carbon that needs to be re-allocated
                                                                                !! after the needs of the reserve and 
                                                                                !! labile pool are satisfied  
                                                                                !! @tex $(gC.m^{-2}dt^{-1})$ @endtex
    REAL(r_std)                                       :: shortage               !! Shortage in the reserves that needs to 
                                                                                !! be re-allocated after to minimise the
                                                                                !! tension between required and available
                                                                                !! reserves
                                                                                !! @tex $(gC.m^{-2}dt^{-1})$ @endtex
    INTEGER                                           :: i,tempi                !   (temp variables for impose intraseasonal LAI dynamic)   

    INTEGER                                           :: month_id               !! index of month 
    REAL(r_std)                                       :: ratio_move             !! temperal variable to move the allocatable carbon
                                                                                !! from leaf to sapwood
    REAL(r_std), DIMENSION(13)                        :: lai_scale              !! monthly lai scaling facter    
    REAL(r_std)                                       :: daily_lai              !! Daily LAI value interpolated by impose lai & lai_scale 
    CHARACTER(len=256)                                :: temp_text              !! dummy text variable exchange

    ! Nitrogen cycle 
    REAL(r_std), DIMENSION(npts,nvm)                  :: n_alloc_tot            !! nitrogen growth (gN/m2/dt) 
    REAL(r_std) , DIMENSION(npts,nvm)                 :: cn_leaf                !! nitrogen concentration in leaves (gC/gN) 
    REAL(r_std) , DIMENSION(npts,nvm)                 :: transloc               !! Transloc variables 
    REAL(r_std), DIMENSION(npts)                      :: alloc_c,alloc_d,alloc_e!! allocation coefficients of nitrogen to leaves, roots and wood 
    REAL(r_std), DIMENSION(npts)                      :: sum_sap,sum_oth        !! carbon growth of wood and root+fruits (gC/m2) 
    REAL(r_std)                                       :: costf                  !! nitrogen cost of a unit carbon growth given current C partitioning and nitrogen concentration 
    REAL(r_std)                                       :: deltacn,deltacnmax     !! (maximum) change in leaf nitrogen concentration 
    REAL(r_std)                                       :: n_avail                !! nitrogen available for growth (dummy) 
    REAL(r_std)                                       :: bm_supply_n            !! carbon growth sustainable by n_avail, considering costf 
    CHARACTER(LEN=2), DIMENSION(nelements)            :: element_str            !! string suffix indicating element 
    REAL(r_std)                                       :: frac_growthresp_dyn    !! Fraction of gpp used for growth respiration (-)
                                                                                !! considers the special case at the leaf onset : frac_growthresp_dyn=0
    REAL(r_std), DIMENSION(npts,nvm)                  :: veget_max_begin        !! temporary storage of veget_max to check area conservation
    REAL(r_std), DIMENSION(npts,nvm)                  :: temp
    REAL(r_std), DIMENSION(ncirc)                     :: ba1                    !! for calculating basal area after phenological
    REAL(r_std), DIMENSION(ncirc)                     :: ba2                    !! growth
    REAL(r_std)                                       :: d_mean                 !! temporal value to calculate deleuze_power
    REAL(r_std)                                       :: deleuze_power          !! denominator of power of delueze-dhote eq.
    REAL(r_std), DIMENSION(ncirc,nparts)              :: temp_mass              !! same with ba1, ba2
    REAL(r_std)                                       :: k_latosa_tmp           !! Temporaray variable in the calculation of k_latosa_adapt
    REAL(r_std)                                       :: optimal_share          !! optimal share between the labile and carbres pools (-)
    REAL(r_std)                                       :: total_reserves         !! Temporary variable in the calculation of the labile and carbres pools.
    REAL(r_std)                                       :: update_sugar_load      !! Instantaneous Relative sugar loading of the labile pool (unitless)
    REAL(r_std)                                       :: n_deficit              !! The amount of nitrogen deficiency in labile nitrogen (dummy)
                                                                                !! only used when quantifying the nitrogen limitation before allocation
    REAL(r_std), DIMENSION(npts,nvm)                  :: height_rel             !! relative height used to make KF dynamic when the stand height increases
    REAL(r_std), DIMENSION(npts,nvm)                  :: residual10b            !! contains the residual (gC tree-1) for warning 10b 
    REAL(r_std), DIMENSION(npts,nvm,ncirc)            :: histvar2               !! temporary variable for xios send
    REAL(r_std), DIMENSION(npts,nvm)                  :: error_count            !! temporary variable to check for errors
!_ ================================================================================================================================

!! 1. Initialize

    !! 1.1 First call only
    IF (firstcall_growth_fun_all) THEN
      !! Initialize local printlev
      ! Guillaume M. -- get_printlev calls getin_p, a COLLECTIVE bcast. Called here,
      ! inside the compute loop under a local firstcall, it hangs the run as soon as the
      ! ranks stop reaching it together. Use printlev, read once at startup by every rank.
      printlev_loc=printlev
      firstcall_growth_fun_all=.FALSE.
    END IF

    IF (printlev_loc.GE.2) WRITE(numout,*) 'Entering functional allocation growth'

    !! 1.3 Initialize variables at every call
    bm_alloc(:,:,:,:) = zero 
    n_alloc_tot(:,:) = zero 
    qm_height(:,:) = zero
    delta_ba = zero
    lai_target(:,:) = zero
    resp_maint(:,:) = zero
    resp_growth(:,:) = zero
    lstress_fac(:,:) = zero
    sigma(:,:) = zero
    gammas(:,:) = zero
    bm_alloc_tot(:,:) = zero
    store_circ_class_ba(:,:,:) = zero
    store_delta_ba_eff(:,:,:) = zero
    store_delta_ba(:,:,:) = zero
    check_intern(:,:,:,:) = zero
    check_intern_init(:,:,:,:) = zero
    excess = zero
    residual(:,:) = zero
    residual_write(:,:) = zero
    residual10b(:,:) = zero
    n_reserve_balance(:,:) = un
    reserve_target(:,:,:) = zero
    labile_target(:,:,:) = zero
    gtemp(:,:) = zero
    circ_class_ba_init(:,:,:) = zero
    height_rel(:,:)= zero
    sla_out(:,:) = zero

    ! If npp is not initialized, bare soil value is never set.
    npp(:,:) = zero 

    ! bare soil never gets set here
    c0_alloc(:,1)=zero
    

    !! 1.4 Initialize check for mass balance closure
    !  The mass balance is calculated at the end of this routine
    !  in section 8
    IF (err_act.GT.1) THEN

       pool_start(:,:,:) = zero
       DO iele = 1,nelements

          ! atm_to_bm has as intent inout, the variable 
          ! accumulates carbon over the course of a day. 
          ! Use the difference between start and the end of 
          ! this routine
          check_intern_init(:,:,iatm2land,iele) = - un * &
               atm_to_bm(:,:,iele) * veget_max(:,:) * dt

          DO ipar = 1,nparts
             DO icir = 1,ncirc
                ! Initial biomass pool
                pool_start(:,:,iele) = pool_start(:,:,iele) + &
                     (circ_class_biomass(:,:,icir,ipar,iele) * &
                     circ_class_n(:,:,icir) * veget_max(:,:))
             ENDDO
          ENDDO
       
       ENDDO      
 
       !! 1.5 Initialize check for surface area conservation
       !  Veget_max is a INTENT(in) variable and can therefore
       !  not be changed during the course of this subroutine
       !  Check it anyway, in case the intent get changed.
       veget_max_begin(:,:) = veget_max(:,:)

    ENDIF ! err_act.GT.1
    
    !! 1.6 Calculate LAI threshold below which carbohydrate reserve is used.
    !  Lai_max and lai_max_to_happy are PFT-dependent parameter specified in 
    !  stomate_constants.f90
    ! +++CHECK+++
    ! Can we make this a function of Cs or rue_longterm? this double prescribed 
    ! value does not make too much sense to me. It is not really dynamic.
    lai_happy(:) = lai_max(:) * lai_max_to_happy(:)
    ! +++++++++++
   
    !! 1.7 Store the biomass pools at the beginning of allocation
    !  These values will be used to calculate the increment in each pool
    !  at the end of the code. tmp_init_bm should not be changed, updated,
    !  or overwritten in this subroutine
    tmp_init_bm(:,:,:,:) = cc_to_biomass(npts,nvm,&
         circ_class_biomass(:,:,:,:,:),&
         circ_class_n(:,:,:))
          
    ! Store basal area and wood volume before the growth to 
    ! calculate basal area increment.    
    DO j = 2,nvm
       IF ( is_tree(j) ) THEN
          DO ipts = 1,npts
             circ_class_ba_init(ipts,j,:) = wood_to_ba(&
                  circ_class_biomass(ipts,j,:,:,icarbon),&
                  j,pipe_tune2(ipts,j))
          ENDDO
       ENDIF
    ENDDO  

    !! 1.8 Calculate C/N ratio of the leaves 
    !  Nitrogen concentration in leaves as CN. First calculate biomass at 
    !  the stand level as that seems quicker than doing these calculations 
    !  on circ_class_biomass.
    tmp_bm = tmp_init_bm
    WHERE( tmp_bm(:,:,ileaf,initrogen).GT.min_stomate .AND. &
         tmp_bm(:,:,ileaf,icarbon).GT.min_stomate)

       ! Calculate the C:N ratio
       cn_leaf(:,:)=tmp_bm(:,:,ileaf,icarbon)/tmp_bm(:,:,ileaf,initrogen)

    ELSEWHERE

       ! Prescribe the C:N ratio
       cn_leaf(:,:)=cn_leaf_min_season(:,:)

    ENDWHERE
   

    !! 1.8 Save old leaf mass
    !  biomass got last updated in stomate_phenology.f90
    lm_old(:,:)=SUM(circ_class_biomass(:,:,:,ileaf,icarbon)*&
         circ_class_n(:,:,:),3)

    !! 1.9 Lai for bare soil is by definition zero
    lai(:,ibare_sechiba) = zero


    !! 2. Use carbohydrate reserve to support growth

    DO j = 2, nvm ! Loop over # PFTs 

       !! 2.1 Calculate demand for carbohydrate reserve to support leaf and root growth.
       !  Maximum time (days) since start of the growing season during which carbohydrate 
       !  may be used
       IF ( is_tree(j) ) THEN

          reserve_time = reserve_time_tree    

       ELSE

          reserve_time = reserve_time_grass

       ENDIF

       !! 2.2 Calculate lai and c0_alloc
       !  The current functions require a loop over npts
       !  Calculate lai

       DO ipts = 1,npts
 
          lai(ipts,j) = cc_to_lai(circ_class_biomass(ipts,j,:,ileaf,icarbon),&
               circ_class_n(ipts,j,:),j)

          ! We might need the c0_alloc factor, so let's calculate it.
          c0_alloc(ipts,j) = calculate_c0_alloc(j, longevity_eff_root(ipts,j), &
               longevity_eff_sap(ipts,j))
       ENDDO

    ENDDO ! loop over # PFTs


 !! 3. Initialize allocation

    DO j = 2, nvm ! Loop over # PFTs 
       
       !! 3.1 Calculate scaling factors, temperature sensitivity, target 
       !  lai to decide on reserve use, labile fraction, labile biomass 
       !  and total allocatable biomass. Convert temperature from K to C
       tl(:) = t2m(:) - ZeroCelsius

       DO ipts = 1, npts

          IF (veget_max(ipts,j) .LE. min_stomate .OR. &
               SUM(circ_class_n(ipts,j,:)) .LE. min_stomate) THEN

              ! This vegetation type is not present, so no reason to do the 
              ! calculation. CYCLE will take us out of the innermost DO loop
               CYCLE

          ENDIF

          !! 3.1 Water stress
          !  The waterstress factor varies between 0.1 and 1 and is calculated
          !  from ::vegstress_season. The latter is only used in the allometric 
          !  allocation and its time integral is determined by longevity_sap for trees
          !  (see constantes_mtc.f90 for longevity_sap and see pft_constantes.f90 for 
          !  the definition of tau_hum_growingseason). The time integral for 
          !  grasses and crops is a prescribed constant (see constantes.f90). For
          !  trees KF (and indirecrtly LF) and for grasses LF are multiplied
          !  by wstress. Because the calculated values are too low for its purpose
          !  Sonke Zhaele multiply it by two in the N-branch (see stomate_season.f90). 
          !  This approach maintains the physiological basis of KF while combining it 
          !  with a simple multiplicative factor for water stress. Clearly after 
          !  multiplication with 2, wstress is closer to 1 and will thus result in a 
          !  KF values closer to the physiologically expected KF. We did not see the 
          !  need to multiply by 2 because the way we now calculate ::vegstress_season
          !  is less volatile than before. Before it ranged between 0 and 1, now the 
          !  range is more like 0.4 to 0.9.

          ! Veget is now calculated from Pgap to be fully consistent within the model. Hence
          ! dividing by veget_max gives a value between 0 and 1 that denotes the amount of 
          ! light reaching the forest floor.
          IF (veget_max(ipts,j) .GT. min_stomate) THEN

             ! Basically recalculate Pgap and use it as the lstress. We do not use
             ! light_tran_to_floor_season here because that variables takes the annual
             ! mean and we want a quicker response here. veget is recalculated daily 
             ! starting from Pgap_cumul (in slowproc.f90)
             lstress_fac(ipts,j) = un - (veget(ipts,j) / veget_max(ipts,j))
             
             ! +++CHECK+++
             ! This is not rocket science so there are a couple
             ! of alternative functions. Some of these functions
             ! try to account for the gaps in the canopy which  has
             ! already been taken care of in Pgap and is reflected
             ! in the ORCHIDEE-CAN way of calculating ::veget. I did
             ! follow these changes.
             ! Alternative 1
             ! lstress_fac(ipts,j) = (un - (veget(ipts,j) / veget_max(ipts,j)))**(0.5)
             ! Alternative 2
             ! lai_temp=-LOG(1-veget(ipts,j))/0.5
             ! veget_temp=1-exp(-0.5*(lai_temp**1.5))
             ! lstress_fac(ipts,j) = (un - (veget_temp ))**(0.3)
             ! Alternative 3
             ! veget_temp=1-exp(-0.5*(lai_temp**3.0))
             ! lstress_fac(ipts,j) = (un - (veget_temp ))**(0.05)
             ! ++++++++++

          ELSE

             lstress_fac(ipts,j) = zero

          ENDIF

          !! 3.2 Initialize scaling factors
          ! Stand level scaling factors
          LF(ipts,j) = 1._r_std

          ! Tree level scaling factors
          ltor(ipts,j) = 1._r_std
          circ_class_height_eff(:) = 1._r_std

          !! 3.3 Calculate structural characteristics
          !  Target lai is calculated at the stand level for the tree 
          !  height of a virtual tree with the mean basal area or the 
          !  so called quadratic mean diameter
          qm_dia(ipts,j) = &
               wood_to_qmdia(circ_class_biomass(ipts,j,:,:,icarbon), &
               circ_class_n(ipts,j,:), j, pipe_tune2(ipts,j))
          qm_height(ipts,j) = &
               wood_to_qmheight(circ_class_biomass(ipts,j,:,:,icarbon), &
               circ_class_n(ipts,j,:), j, pipe_tune2(ipts,j))

          !! 3.4 Calculate allocation factors for trees and grasses 
          IF ( SUM(SUM(circ_class_biomass(ipts,j,:,:,icarbon),1)) .GT. min_stomate ) THEN

             ! Note that KF may already be calculated in stomate_prescribe.f90 (if called)
             ! it is recalculated because the biomass pools for grasses and crops
             ! may have been changed in stomate_phenology.f90. Trees were added to this
             ! calculation just to be consistent.
             
             ! Scaling factor to convert sapwood mass into leaf mass (KF)
             ! derived from
             ! LA_ind = k1 * SA_ind, k1=latosa (pipe-model)
             ! <=> Cl * vm/ind * sla = k1 * Cs * vm/ind / wooddens / tree_ff / height_new
             ! <=> Cl = Cs * k1 / wooddens / tree_ff/ height_new /sla
             ! <=> Cl = Cs * KF / height_new, where KF = k1 / (wooddens * sla * tree_ff)
             ! (1) Cl = Cs * KF / height_new
             KF_old = KF(ipts,j)
             
             ! To be fully consistent with the hydraulic limitations and pipe theory,
             ! k_latosa_zero should be calculated from equation (18) in Magnani et al.
             ! To do so, total hydraulic resistance and tree height need to be known. This
             ! poses a problem as the resistance depends on the leaf area and the leaf 
             ! area on the resistance. There is no independent equation and equations 12
             ! and 18 depend on each other and substitution would be circular. Hence 
             ! prescribed k_latosa_adapt values were obtained from observational records 
             ! and are given in mtc_parameters.f90

             ! The most simple approach to estimate k_latosa is by prescribing it. Note
             ! that for the moment lstress = 0. We decided to keep k_latosa_min and
             ! k_latosa_max just in case we want to test more complex relationships. Note
             ! as well that in the parameter files k_latosa_max = k_latosa_min. 
             ! This approach is not fully able to compensate for the increase in height
             ! Cl = KF*Cs/height. If height increases, KF should increase as well
             ! to maintain the lai. Lstress = Pgap and saturates above an
             ! lai of 4-5. If Lai drops from 7 to 6, this approach does not respond
             ! sufficiently. Part of this drop was found to be due to a quick drop in
             ! N-availability during the spinup (that is the purpose of the spinup). If the
             ! model is resarted after a clear cut (r7250), this big drop in lai largely 
             ! disappears and lai decreases 0.5 to 1.5 units over a 200 year long simulation. 
             ! That is considered acceptable. 
             k_latosa_tmp = (k_latosa_adapt(ipts,j) + (lstress_fac(ipts,j) * &
                  (k_latosa_max(j)-k_latosa_min(j))))
             
             !+++ALTERNATIVES+++
             ! At one point it looked like a good idea to take the max of two
             ! options but by doing so we cannot recalculate KF in phenology.
             ! It has not been confirmed that this is really a problem. Use the
             ! simpelest approach but leave the alternative in the code as a
             ! suggestion of a possible solution in case something goes wrong.
             !!$             k_latosa_tmp = MAX(k_latosa_min(ivm),k_latosa_adapt(ipts,ivm) + &
             !!$                 (lstress_fac(ipts,ivm) * &
             !!$                 (k_latosa_max(ivm)-k_latosa_min(ivm))))
             
             ! The relationship between height and k_latosa as reported in McDowell 
             ! et al 2002 and Novick et al 2009 is implemented to adjust k_latosa for 
             ! the height of the stand. The slope of the relationship is calculated in 
             ! stomate_data.f90 This did NOT result in a realistic model behavior.
             !!$             k_latosa(ipts,j) = wstress_fac(ipts,j) * & 
             !!$            (k_latosa_max(j) - latosa_height(j) * qm_height(ipts,j))
             ! Another relationship with height was implemented. This resulted in acceptable
             ! model behavior (r7250) but had very little impact on the temporal patterns.
             ! For that reason the most simple formulation was favored over this approach
             !!$             height_rel(ipts,j) = MAX(MIN((qm_height(ipts,j)/pipe_tune2(j))**
             !!$                  (1/exp_kf),un),zero)
             !!$             k_latosa_tmp = k_latosa_adapt(ipts,j) + height_rel(ipts,j) * &
             !!$                  (k_latosa_max(j)-k_latosa_min(j))
             
             ! Alternatively, k_latosa is also reported to be a function of diameter 
             ! (i.e. stand thinning, Simonin et al 2006, Tree Physiology, 26:493-503).
             ! Here the relationship with thinning was interpreted as a realtionship with
             ! light stress. This is the same formulation as we use now but to make it 
             ! function l_stress should be calculated and the parameters for k_latosa_min
             ! and k_latosa_max should differ from each other.
             ! k_latosa(ipts,j) = (k_latosa_adapt(ipts,j) + &
             !     (lstress_fac(ipts,j) * &
             !     (k_latosa_max(j)-k_latosa_min(j))))
    
             ! Also k_latosa has been reported to be a function of CO2 concentration 
             ! (Atwell et al. 2003, Tree Physiology, 23:13-21 and Pakati et al. 2000, 
             ! Global Change Biology, 6:889-897). This effect is not accounted for in 
             ! the current code

             ! How dow we want to account for waterstress?
             !!$             k_latosa(ipts,j) = k_latosa_min(j) + (wstress_fac(ipts,j) * &
             !!$             lstress_fac(ipts,j) * &
             !!$             (k_latosa_max(j)-k_latosa_min(j)))
             !!$             k_latosa(ipts,j) = wstress_fac(ipts,j) * (k_latosa_min(j) + &
             !!$             (lstress_fac(ipts,j) * &
             !!$             (k_latosa_max(j)-k_latosa_min(j))))
             !++++++++++++++++++

             ! Calculate the sla for the current amount of leaf biomass. Use a trick.
             ! use biomass_to_lai to calculate the lai with a dynamic or a 
             ! static sla calculation. Then divide the lai by the biomass to
             ! obtain the actual value for sla (m2 g-1).
             IF (sla_dyn) THEN
                IF (tmp_init_bm(ipts,j,ileaf,icarbon).GT.min_stomate) THEN
                   ! Calculate a dynamic sla
                   sla_est = biomass_to_lai(tmp_init_bm(ipts,j,ileaf,icarbon),j) / &
                        tmp_init_bm(ipts,j,ileaf,icarbon)
                ELSE
                   ! Nothing changes, calculate sla_est such that KF will remain
                   ! the same in the KF calculation below this IF-statement.
                   sla_est = k_latosa_tmp / &
                        (KF_old *  pipe_density(j) * tree_ff(j))
                ENDIF
             ELSE
                ! Use the prescribed fixed sla
                sla_est = sla(j)
             ENDIF

             ! save SLA as output
             sla_out(ipts,j) = sla_est
             
             ! Calculate the actual KF
             KF(ipts,j) = k_latosa_tmp / &
                  (sla_est * pipe_density(j) * tree_ff(j))
             
             ! KF of the previous time step was stored in ::KF_old to check its absolute 
             ! change. If this absolute change is too big the whole allocation will crash
             ! because it will calculate negative increments which are compensated by 
             ! positive increments that exceed the available carbon for allocation. This 
             ! would suggest that for example the plant destroys leaves and uses the 
             ! available carbon to produce more roots. This would represent an unwanted 
             ! outcome. Large changes from one time step to another makes it difficult for 
             ! the scheme to ever reach allometric balance. This balance is needed for the
             ! allocation scheme to allow 'ordinary allocation', which in turn is needed
             ! to make use of the allocation rule of Dhote and Deleuze. It needs to be 
             ! avoided that the code spends too much time in phenological growth and the 
             ! if-then statements that help to restore allometric balance. For this reason
             ! the absolute changes in KF from one time step to another are truncated.
             IF (KF_old - KF(ipts,j) .GT. max_delta_KF ) THEN
                
                IF(printlev_loc>=4)THEN
                   WRITE(numout,*) 'WARNING 2: KF was truncated'
                   WRITE(numout,*) 'WARNING 2: PFT, ipts: ',j,ipts
                   WRITE(numout,'(A,3F20.10)') 'WARNING 2: KF_old, KF(ipts,j), '//&
                        'max_delta_KF: ', KF_old, KF(ipts,j), max_delta_KF
                ENDIF
                
                ! Add maximum absolute change
                KF(ipts,j) = KF_old - max_delta_KF
                
                IF(printlev_loc>=4)THEN
                   WRITE(numout,'(A,3F20.10)') 'WARNING 2: Reset, KF_old, KF(ipts,j): ',&
                        KF_old, KF(ipts,j)
                ENDIF
                
             ELSEIF (KF_old - KF(ipts,j) .LT. -max_delta_KF) THEN
                
                IF(printlev_loc>=4)THEN
                   WRITE(numout,*) 'WARNING 3: KF was truncated'
                   WRITE(numout,*) 'WARNING 3: PFT, ipts: ',j,ipts
                   WRITE(numout,'(A,3F20.10)') 'WARNING 3: KF_old, KF(ipts,j), '//&
                        'max_delta_KF: ', KF_old, KF(ipts,j), -max_delta_KF
                ENDIF
                
                ! Remove maximum absolute change
                KF(ipts,j) = KF_old + max_delta_KF
                
                IF(printlev_loc>=4)THEN
                   WRITE(numout,'(A,3F20.10)') 'WARNING 3: Reset, KF_old, KF(ipts,j): ',&
                        KF_old, KF(ipts,j)
                ENDIF
                
             ELSE
                
                ! The change in KF is acceptable no action required
                
             ENDIF
             
             ! Scaling factor to convert sapwood mass into root mass  (LF) 
             ! derived from 
             ! Cs = c0 * height * Cr (Magnani 2000)
             ! Cr = Cs / c0 / height_new
             ! scaling parameter between leaf and root mass, derived from
             ! Cr = Cs / c0 / height_new
             ! let Cs = Cl / KF * height_new
             ! <=> Cr = ( Cl * height_new / KF ) / ( c0 * height_new )
             ! <=> Cl = Cr * KF * c0
             ! <=> Cl = Cr * LF, where LF = KF * c0
             LF(ipts,j) = c0_alloc(ipts,j) * KF(ipts,j) 
             
             ! Calculate non-nitrogen stressed leaf to root ratio to calculate the
             ! allocation to the reserves. Should be multiplied by a nitrogen stress
             ! have a look in OCN. This code should be considered as a placeholder
             ltor(ipts,j) = c0_alloc(ipts,j) * KF(ipts,j)
             
             ! Debug
             IF (j.EQ.test_pft .AND. printlev_loc.GE.4 .AND. ipts.EQ.test_grid) THEN
                WRITE(numout,*) 'Updating KF and related variables'
                WRITE(numout,*) 'KF, LF, ', KF(ipts,j), LF(ipts,j)
                WRITE(numout,*) 'c0_alloc, ', c0_alloc(ipts,j)
                WRITE(numout,*) 'longevity_root, longevity_sap, ', longevity_eff_root(ipts,j), &
                     longevity_eff_sap(ipts,j)
                WRITE(numout,*) 'k_belowground, k_sap, ', k_belowground(j), k_sap(j)
                WRITE(numout,*) 'ltor, ', ltor(ipts,j)
             ENDIF
             !-

          ENDIF ! SUM(circ_class_biomass) .gt. zero

          !+++CHECK+++
          !! 3.5 Calculate optimal LAI
          !  The calculation of the optimal LAI was copied and adjusted from O-CN. 
          !  In O-CN it was also used in the allocation but that seems to be 
          !  inconsistent with the allometric rules that are implemented. Say that 
          !  the actual LAI is below the optimal LAI. Then the O-CN approach will 
          !  keep pumping carbon to grow the optimal LAI. If we would apply
          !  the same method it means that during this phase the rule of Deleuze 
          !  and Dhote would not be used. For that reason we dropped the use of 
          !  LAI_optimal and replaced it by an allometric-based Cl_target value. 
          !  Initially, lai_target was still calculated as described below and used 
          !  in the calculation of the reserves. Further testing showed that for 
          !  some parameter sets lai_target was over 8 whereas the realized lai was 
          !  close to 4. This leaves us with a frustrated plant that will invest a
          !  lot in its reserves but can never use them because it is constrained by 
          !  the allometric rules. To grow an LAI of 8 it would need to have a crazy 
          !  sapwoodmass. At a more fundamental level it is clear why the plant's 
          !  LAI should not exceed lai_target because then it costs more to produce 
          !  and maintain the leaf than what the new leaf can produce but there is no 
          !  reason why the plant should try to reach lai_target. For these
          !  reasons it was decided to abandon this approach to lai_target and simply 
          !  replace lai_target by Cl_target * sla

          !! 3.5.1 Scaling factor
          !  Scaling factor to convert variables to the individual plant
          !  Different approach between the DGVM and statitic approach
          IF (ok_dgvm) THEN

             ! The DGVM does currently NOT work with the new allocation, consider this as
             ! placeholder. The original code had two different transformations to 
             ! calculate the scalars. Both could be used but the units will differ.
             ! When fixing the DGVM check which quantities need to be multiplied by scal 
             ! scal = ind(ipts,j) * cn_ind(ipts,j) / veget_max(ipts,j)
             scal(ipts,j) = veget_max(ipts,j) / SUM(circ_class_n(ipts,j,:)) 
 
          ELSE

             ! circ_class_biomass contain the data at the tree level
             ! no conversion required
             scal(ipts,j) = 1.

          ENDIF

          !! 3.5.2 Calculate lai_target based on the allometric rules
          IF ( SUM(SUM(circ_class_biomass(ipts,j,:,:,icarbon),1)) .GT. min_stomate ) THEN
          
             IF ( is_tree(j)) THEN

                ! Basal area at the tree level (m2 tree-1)
                circ_class_ba_eff(:) = wood_to_ba_eff(circ_class_biomass(ipts,j,:,:,icarbon),&
                     j, pipe_tune2(ipts,j))

                ! Current biomass pools per tree (gC tree^-1) 
                ! We will have different trees so this has to be calculated from the 
                ! diameter relationships            
                Cs(:) = ( circ_class_biomass(ipts,j,:,isapabove,icarbon) + &
                     circ_class_biomass(ipts,j,:,isapbelow,icarbon) ) * scal(ipts,j)
                Cr(:) = circ_class_biomass(ipts,j,:,iroot,icarbon) * scal(ipts,j)
                Cl(:) = circ_class_biomass(ipts,j,:,ileaf,icarbon) * scal(ipts,j)
                Ch(:) = ( circ_class_biomass(ipts,j,:,iheartabove,icarbon) + &
                     circ_class_biomass(ipts,j,:,iheartbelow,icarbon) ) * scal(ipts,j)

                DO l = 1,ncirc 

                   !  Calculate tree height
                   circ_class_height_eff(l) = pipe_tune2(ipts,j)*(4/pi*circ_class_ba_eff(l))**&
                        (pipe_tune3(j)/2)

                   ! Debug
                   IF(printlev>=4 .AND. ipts.EQ.test_grid .AND. j.EQ.test_pft)THEN
                      WRITE(numout,*) 'circ_class_height, ',circ_class_height_eff(l)
                      WRITE(numout,*) 'KF, ',  KF(ipts,j)
                   ENDIF
                   !-

                   !  Do the biomass pools respect the pipe model?
                   !  Do the current leaf, sapwood and root components respect the allometric 
                   !  constraints? Due to plant phenology it is possible that we have too much 
                   !  sapwood compared to the leaf and root mass (i.e. in early spring). 
                   !  Calculate the optimal root and leaf mass, given the current wood mass 
                   !  by using the basic allometric relationships. Calculate the optimal sapwood
                   !  mass as a function of the current leaf and root mass.
                   Cl_target(l) = MAX( KF(ipts,j) * Cs(l) / circ_class_height_eff(l), &
                        Cr(l) * LF(ipts,j) , Cl(l))
                   Cs_target(l) = MAX( Cl(l) / KF(ipts,j) * circ_class_height_eff(l), &
                        Cr(l) * LF(ipts,j) / KF(ipts,j) * circ_class_height_eff(l) , Cs(l))

                   ! Check dimensions of the trees
                   ! If Cs = Cs_target then ba and height are correct, else calculate the 
                   ! correct dimensions

                   IF ( Cs_target(l) - Cs(l) .GT. min_stomate ) THEN

                      ! If Cs = Cs_target then dia and height are correct. However, if 
                      ! Cl = Cl_target or Cr = Cr_target then dia and height need to be 
                      ! re-estimated. Cs_target should satify the relationship 
                      ! Cl/Cs = KF/height where height is a function of Cs_target
                      ! Search Cs needed to sustain the max of Cl or Cr.
                      ! Search max of Cl and Cr first
                      !
                      ! [UPDATE] After the code passes through turnover or mortality
                      ! we may end up in a situation where we have lost more
                      ! sapwood than leaves and roots (i.e. sapwood turnover).
                      ! The model would then suggest that at time=t+1 the tree
                      ! should be smaller than at time=t0.  From a physiological
                      ! standpoint this is not possible for the heartwood. If
                      ! we now calculate Cs_target on the basis of the actual Cl
                      ! or Cr, we find that Cs_target > Cs. The first priority
                      ! of the allocation scheme will be to allocate C to Cs.
                      ! Because we don't know yet whether the actual Cr or Cl is
                      ! what drives the need to allocate to Cs, we calculate
                      ! Cl_target first. 
                      Cl_target(l) = MAX(Cl(l), Cr(l)*LF(ipts,j))

                      ! Debug
                      IF (ipts.EQ.test_grid .AND. j.EQ.test_pft .AND. printlev_loc.GE.4) THEN
                         WRITE(numout,*) 'Does the tree need reshaping? ipts, class: ', &
                              ipts,l
                         WRITE(numout,*) 'circ_class_height_eff, ', circ_class_height_eff(l)
                         WRITE(numout,*) 'KF, LF, ', KF(ipts,j), LF(ipts,j)
                         WRITE(numout,*) 'Cl_target-Cl, ', Cl_target(l)-Cl(l), Cl_target(l), Cl(l)
                         WRITE(numout,*) 'Cs, ', Cs(l)
                         WRITE(numout,*) 'Cr, ', Cr(l)
                         WRITE(numout,*) 'Ch, ', Ch(l)
                      ENDIF
                      !-

                      ! We now have the Cl_target that we will use to calculate
                      ! Cs_target. Given the allometric relationships we can
                      ! calculate Cs_target as Cl_target*height/KF. 
                      ! height is a function of ba, which in turn is a function 
                      ! of woodmass (Woodmass = Cs+Ch (sapwood+heartwood) ). We 
                      ! therefore substitute the following equations into one 
                      ! another:
                      ! (1) Cs_target = Cl_target*height/KF
                      ! (2) height = as a function of ba
                      ! (3) ba = as a function of woodmass_ind
                      !
                      ! This gives:
                      ! (4) Cl_target = (KF*Cs_target)/(pipe_tune2*(Cs_target+Ch)/ &
                      !                 & pi/4)**(pipe_tune3/(2+pipe_tune3))
                      !
                      ! The function newX searches for the value for Cs_target
                      ! that satisfies this equation (4).
                      Cs_target(l) =  newX(KF(ipts,j), Ch(l),pipe_tune2(ipts,j), &
                           & pipe_tune3(j), Cl_target(l), Cs_target(l),&
                           & tree_ff(j)*pipe_density(j)*pi/4*pipe_tune2(ipts,j), Cs(l),&
                           & 2*Cs(l), 2, j, ipts)

                      ! Recalculate height and ba from the correct
                      ! Cs_target
                      circ_class_height_eff(l) = Cs_target(l)*KF(ipts,j)&
                           &/Cl_target(l)
                      circ_class_ba_eff(l) = pi/4*(circ_class_height_eff(l)&
                           &/pipe_tune2(ipts,j))**(2/pipe_tune3(j))
                      Cl_target(l) = KF(ipts,j) * Cs_target(l) /&
                           & circ_class_height_eff(l)
                      Cr_target(l) = Cl_target(l) / LF(ipts,j)

                      ! Debug
                      IF (printlev_loc.GE.4 .AND. ipts.EQ.test_grid .AND. j.EQ.test_pft) THEN
                         WRITE(numout,*) 'New Cl_target, ', Cl_target(l)
                         WRITE(numout,*) 'New Cs_target, ', Cs_target(l)
                         WRITE(numout,*) 'New Cr_target, ', Cr_target(l)       
                      ENDIF
                      !- 

                   ENDIF

                ENDDO

                ! Calculate lai_target
                lai_target(ipts,j) = cc_to_lai(Cl_target(:),circ_class_n(ipts,j,:),j)

             ELSEIF ( .NOT. is_tree(j)) THEN
   
                ! Grasses and croplands
                ! Current biomass pools per grass/crop (gC ind^-1)
                ! Cs has too many dimensions for grass/crops. To have a consistent 
                ! notation the same variables are used as for trees but the dimension 
                ! of Cs, Cl and Cr i.e. ::ncirc should be ignored            
                Cs(1) = circ_class_biomass(ipts,j,1,isapabove,icarbon) * scal(ipts,j)
                Cr(1) = circ_class_biomass(ipts,j,1,iroot,icarbon) * scal(ipts,j)
                Cl(1) = circ_class_biomass(ipts,j,1,ileaf,icarbon) * scal(ipts,j)
                Ch(1) = zero

                ! Do the biomass pools respect the pipe model?
                ! Do the current leaf, sapwood and root components respect the allometric 
                ! constraints? Calculate the optimal root and leaf mass, given the current 
                ! wood mass by using the basic allometric relationships. Calculate the 
                ! optimal sapwood mass as a function of the current leaf and root mass.
                Cl_target(1) = MAX( Cs(1) * KF(ipts,j) , Cr(1) * LF(ipts,j), Cl(1) )
                Cs_target(1) = MAX( Cl_target(1) / KF(ipts,j), &
                     Cr(1) * LF(ipts,j) / KF(ipts,j), Cs(1) ) 
                Cr_target(1) = MAX( Cl_target(1) / LF(ipts,j), &
                     Cs_target(1) * KF(ipts,j) / LF(ipts,j), Cr(1) )

                ! Calculate lai_target
                lai_target(ipts,j) = cc_to_lai(Cl_target(:),circ_class_n(ipts,j,:), j)

                ! Debug
                IF (printlev_loc.GE.4 .AND. ipts.EQ.test_grid .AND. j.EQ.test_pft) THEN
                   WRITE(numout,*) 'KF, LF, ', KF(ipts,j), LF(ipts,j)
                ENDIF
                !-

             ENDIF

          ELSE

             ! circ_class_biomass is empty  
             lai_target(ipts,j) = zero

          ENDIF ! SUM(circ_class_biomass) .GT. min_stomate

!!$          !! 3.5 Calculate optimum LAI 
!!$          !  Lai is optimised for mean annual radiation use efficiency and the C costs
!!$          !  for producing the canopy. The cost-benefit ratio is optimised when the 
!!$          !  marginal gain / marginal cost = 1
!!$          !  Investing 1 gC in the canopy comes at a total cost that is composed by the 
!!$          !  C required for the canopy in addition to the roots and the sapwood to support
!!$          !  the canopy. The total cost (C) is thus calculated as C: 
!!$          !  LAI/sla * ( (one_year/longevity_leaf) + (one_year/longevity_root)/LF + & 
!!$          !  (one_year/longevity_sap)*height/KF))
!!$          !  The marginal cost for one unit of LAI is then dC/dLAI : 
!!$          !  (one_year/longevity_leaf)/sla + (one_year/longevity_root)/LF/sla + é
!!$          !  (one_year/longevity_sap)*height/KF/sla)
!!$          !  Where, longevity_leaf is given by ::longevity_leaf in days, longevity_root by ::longevity_root in
!!$          !  days and longevity_sap by ::longevity_sap in days. LF is unitless, KF is expressed in meters
!!$          !  and sla in m^2.gC^{-1}. The unit of dC/dLAI is thus gC.m^{-2} but all turnover
!!$          !  times need to be expressed on an annual scale.  
!!$          !  Investing 1gC in the canopy enables the plant to assimilate more carbon
!!$          !  The gain (G) can be approximated by using the 'radiation use efficiency' as 
!!$          !  follows: RUE * one_year ( 1. - exp (-0.5 * LAI ))
!!$          !  Where, 0.5 is the extinction factor that accounts for the fact the lower parts 
!!$          !  of the canopy receive less light. Note that RUE has a peculiar definition 
!!$          !  and is calculated as the ratio of GPP over the fraction of radiation 
!!$          !  absorbed by the canopy.
!!$          !  Hence the unit of RUE is gC.m^{-2}.day^{-1}. The marginal gain of one 
!!$          !  unit of LAI is dG/dLAI: 
!!$          !  0.5 * one_year * RUE * exp (-0.5 * LAI). 
!!$          !  Subsequently, the optimal LAI is approximated by
!!$          !  LAI_opt = -2. * log(2*(dC/dt)/(RUE*one_year))           
!!$          !  Added the qm_height requirement since for a grass, it had no biomass
!!$          !  but it did have individuals.  This caused qm_height to be zero and a crash
!!$          !  in the calculation of lai_target.
!!$          IF ( (rue_longterm(ipts,j) .GT. min_stomate) .AND. (ind(ipts,j) .NE. zero &
!!$               .AND. qm_height(ipts,j) .NE. 0) ) THEN
!!$
!!$             ! Scheme in line with the documentation
!!$             lai_target(ipts,j) = -deux* log( (deux * (one_year/longevity_leaf(j))/sla(j) + &
!!$                  ((one_year/longevity_root(j))/LF(ipts,j))/sla(j) + &
!!$                  ((one_year/longevity_sap(j))*qm_height(ipts,j)/KF(ipts,j))/sla(j)) / &
!!$                  (rue_longterm(ipts,j)*one_year))
!!$             lai_target(ipts,j) = MAX(MIN(lai_target(ipts,j),12.),.5)
!!$
!!$          ELSE
!!$
!!$             lai_target(ipts,j) = 0.5
!!$
!!$          ENDIF
          !++++++++++

          !! 3.6 Calculate mean leaf age
          leaf_meanage = zero
          DO m = 1,nleafages
          
             leaf_meanage = leaf_meanage + &
                  leaf_age(ipts,j,m) * leaf_frac(ipts,j,m)
          
          ENDDO

          !! 3.8 Calculate total allocatable biomass during this time step determined from GPP.
          ! It is bit easier to deal with this issue at the stand level because
          ! gpp_daily, the reserve, and labile pool are calculated at the 
          ! stand level. After making stand level calculations the updated
          ! pools will have to be converted to the plant level.
          ! Calculate stand level biomass
          tmp_bm(ipts,j,:,:) = cc_to_biomass(npts,j,&
               circ_class_biomass(ipts,j,:,:,:),&
               circ_class_n(ipts,j,:))
          
          ! Debug
          IF(printlev_loc.GE.3 .AND. ipts == test_grid .AND. j == test_pft)THEN
             WRITE(numout,*) 'Initial labile and carbres pools'
             WRITE(numout,*) 'gpp_daily(ipts,j)',gpp_daily(ipts,j)
             WRITE(numout,*) 'biomass(ipts,j,ilabile,icarbon)',&
                  tmp_bm(ipts,j,ilabile,icarbon)
             WRITE(numout,*) 'biomass(ipts,j,icarbres,icarbon)',&
                  tmp_bm(ipts,j,icarbres,icarbon)
             WRITE(numout,*) 'biomass(ipts,j,ilabile,initrogen)',&
                  tmp_bm(ipts,j,ilabile,initrogen)
             WRITE(numout,*) 'biomass(ipts,j,icarbres,initrogen)',&
                  tmp_bm(ipts,j,icarbres,initrogen)
          ENDIF
          !-

          ! If plant goes to senescence or ipresenescence, gpp
          ! goes to reserve pool. This is to make the fresh GPP not readily
          ! available for allocation, in order to preserve the reserve pool.
          IF ( plant_status(ipts,j).EQ.isenescent .OR. &
               plant_status(ipts,j).EQ.ipresenescence ) THEN

             ! The plant is in senescence or pre-senescence: icarbres should be used
   
             ! GPP was calculated as CO2 assimilation in enerbil.f90
             ! Under some exceptional conditions :gpp could be negative when
             ! the dark respiration exceeds the photosynthesis. When this happens
             ! the dark respiration is paid for by the labile and carbres pools
             ! Account for dark respiration if needed
             IF ( (tmp_bm(ipts,j,icarbres,icarbon) + &
                  gpp_daily(ipts,j) * dt) .LT. zero ) THEN
                
                deficit = (tmp_bm(ipts,j,icarbres,icarbon) + gpp_daily(ipts,j) * dt)

                ! The deficit is less than the carbon reserve
                IF (-deficit .LE. tmp_bm(ipts,j,ilabile,icarbon)) THEN
                   
                   ! Pay the deficit from the reserve pool
                   tmp_bm(ipts,j,ilabile,icarbon) = &
                        tmp_bm(ipts,j,ilabile,icarbon) + deficit
                   tmp_bm(ipts,j,icarbres,icarbon) = &
                        tmp_bm(ipts,j,icarbres,icarbon) - deficit

                ELSE

                   ! Not enough carbon to pay the deficit, the individual 
                   ! is going to die at the end of this day
                   tmp_bm(ipts,j,icarbres,icarbon) = &
                        tmp_bm(ipts,j,ilabile,icarbon) + &
                        tmp_bm(ipts,j,icarbres,icarbon) 
                   tmp_bm(ipts,j,ilabile,icarbon) = zero
                   
                   ! Truncate the dark respiration to the available carbon.  Now we
                   ! should use up all the reserves.  If the plant has no leaves, it
                   ! will die quickly after this.
                   gpp_daily(ipts,j) = - tmp_bm(ipts,j,icarbres,icarbon)/dt 

                ENDIF

             ENDIF ! labile pool is empty
             
             ! Not senescent add GPP (irrespective of whether it is positive or 
             ! negativeto labile pool
             tmp_bm(ipts,j,icarbres,icarbon) = tmp_bm(ipts,j,icarbres,icarbon) + &
                  gpp_daily(ipts,j) * dt

          ELSE
            
             ! The plant is still growing: gpp should go into ilabile

             ! GPP was calculated as CO2 assimilation in enerbil.f90
             ! Under some exceptional conditions :gpp could be negative when
             ! the dark respiration exceeds the photosynthesis. When this happens
             ! the dark respiration is paid for by the labile and carbres pools
             ! Account for dark respiration if needed
             IF ( (tmp_bm(ipts,j,ilabile,icarbon) + &
                  gpp_daily(ipts,j) * dt) .LT. zero ) THEN
                
                deficit = (tmp_bm(ipts,j,ilabile,icarbon) + gpp_daily(ipts,j) * dt)

                ! The deficit is less than the carbon reserve
                IF (-deficit .LE. tmp_bm(ipts,j,icarbres,icarbon)) THEN
                   
                   ! Pay the deficit from the reserve pool
                   tmp_bm(ipts,j,icarbres,icarbon) = &
                        tmp_bm(ipts,j,icarbres,icarbon) + deficit
                   tmp_bm(ipts,j,ilabile,icarbon) = &
                        tmp_bm(ipts,j,ilabile,icarbon) - deficit

                ELSE

                   ! Not enough carbon to pay the deficit, the individual 
                   ! is going to die at the end of this day
                   tmp_bm(ipts,j,ilabile,icarbon) = &
                        tmp_bm(ipts,j,ilabile,icarbon) + &
                        tmp_bm(ipts,j,icarbres,icarbon) 
                   tmp_bm(ipts,j,icarbres,icarbon) = zero
                   
                   ! Truncate the dark respiration to the available carbon.  Now we
                   ! should use up all the reserves.  If the plant has no leaves, it
                   ! will die quickly after this.
                   gpp_daily(ipts,j) = - tmp_bm(ipts,j,ilabile,icarbon)/dt 

                ENDIF

             ENDIF ! labile pool is empty 
             
             ! Not senescent add GPP (irrespective of whether it is positive or 
             ! negativeto labile pool
             tmp_bm(ipts,j,ilabile,icarbon) = tmp_bm(ipts,j,ilabile,icarbon) + &
                  gpp_daily(ipts,j) * dt

          ENDIF ! check plant_status 
          
          ! This would be a good place to update the plant level
          ! labile and carbres pool but as it may be subject to
          ! further modifications in the next block of code, it 
          ! will be done later 

          ! Debug
          IF(printlev_loc.GE.3 .AND. ipts == test_grid .AND. j == test_pft)THEN
             WRITE(numout,*) 'Added gpp to labile pool'
             WRITE(numout,*) 'gpp_daily(ipts,j)',gpp_daily(ipts,j)
             WRITE(numout,*) 'biomass(ipts,j,ilabile,icarbon)',&
                  tmp_bm(ipts,j,ilabile,icarbon)
             WRITE(numout,*) 'biomass(ipts,j,ilabile,initrogen)',&
                  tmp_bm(ipts,j,ilabile,initrogen)
          ENDIF
          !-

          !! 3.9 Calculate activity of labile carbon pool  
             
          ! Similar relationship as that used for the temperature 
          ! response of maintenance respiration but the parameters were 
          ! tuned to reflect a temperature-growth relationships. 
          ! The parameters in the equation were calibrated to 
          ! give no growth below -2, at 0 degrees only 3% of the labile
          ! pool can be allocated to growth and at 5 degrees 100% of the 
          ! labile pool can be allocated to growth if there is enough 
          ! nitrogen. This equation partly decouples growth and gpp for
          ! days that it is warm enough of gpp (above 0 degrees) but 
          ! probably too cold to grow (above 5 degrees). This is our 
          ! partial answer to the sink/source discussion bu Fatichi et al 
          ! 2013 in New Phytologist. Note that this has little effect on 
          ! the evergreen species and results in a couple of sudden 
          ! small dips in for example the LAI in winter in the temperate zone.
          ! This approach has even less effect on the deciduous species 
          ! because for those species phenology typically happens after
          ! the 5 degree threshold has been passed. At the time the 
          ! deciduous trees get their leaves, gpp and growth are coupled
          ! to the extent that there is enough nitrogen to support the 
          ! growth.
          IF (tl(ipts) .GT. tmin_labile(j)) THEN 
             gtemp(ipts,j) = EXP((e0_labile(j))*(1.0/(tref_labile(j)-tmin_labile(j)) - &
                  1.0/(tl(ipts)-tmin_labile(j))))
          ELSE
             ! Too cold to grow
             gtemp(ipts,j) = zero            
          ENDIF
          
          !  If there is a plant, and we are either at the very start or in 
          !  the growing season not during senescences, calculate labile 
          !  pool use for growth.
          ! CYmark: we allow calculation of bm_alloc_tot as well for ipresenescence
          ! stage.
          IF ( SUM(circ_class_n(ipts,j,:)) .GT. min_stomate .AND. &
               ( plant_status(ipts,j) .EQ. ibudbreak .OR. &
                 plant_status(ipts,j) .EQ. icanopy .OR. &
                 plant_status(ipts,j) .EQ. ipresenescence ) .AND. &
               SUM(circ_class_biomass(ipts,j,:,ileaf,icarbon)) .GT. min_stomate ) THEN

             IF ( (tmp_bm(ipts,j,ilabile,icarbon) .GT. min_stomate) .OR. &
                  (tmp_bm(ipts,j,icarbres,icarbon) .GT. min_stomate) ) THEN
                
                ! Truncate gtemp between zero and 1. If we set the upper
                ! bound to one, we may run into numerical (precision) problems
                ! later caused by very small (10e-15) negative values. Rather
                ! than dealing with the precision issues it is easier to use
                ! 0.99 instead. 0.99 may be too high so this parameter was
                ! externalized and is pft-dependent.
                gtemp(ipts,j) = MAX(MIN(gtemp(ipts,j), un-always_labile(j)), zero)           

             ELSE

                ! There is nothing to allocate so we could as well set
                ! gtemp to zero
                gtemp(ipts,j) = zero

             ENDIF
            
             ! Prioritize the use of the carbohydrate pool. Move
             ! carbohydrates to the labile pool.
             ! CYmark: Such moving reserve  to labile pool is not allowed for 
             ! ipresenescence stage.
             IF ( ( plant_status(ipts,j) .EQ. ibudbreak .OR. &
                    plant_status(ipts,j) .EQ. icanopy ) .AND. &
                  (tmp_bm(ipts,j,icarbres,icarbon) .GT. min_stomate) ) THEN

                tmp_bm(ipts,j,ilabile,icarbon) = &
                     tmp_bm(ipts,j,ilabile,icarbon) + 0.05 * &
                     tmp_bm(ipts,j,icarbres,icarbon)
                tmp_bm(ipts,j,icarbres,icarbon) = & 
                     tmp_bm(ipts,j,icarbres,icarbon) * 0.95
             ENDIF

             ! This would be a good place to update the plant level
             ! labile and carbres pool but as it may be subject to
             ! further modifications in the next block of code, it 
             ! will be done later

          ELSE
             
             ! The plant is absent, senescent or dead so there will be 
             ! no allocation. Set gtemp to zero to keep the output files 
             ! clean. Due to this line, gtemp will be zero as soon as
             ! plant_status becomes isenescent or idormant.
             gtemp(ipts,j) = zero

          ENDIF

          ! Since the plant is in ipresenescence, we want to stop allocating gpp (i.e.,
          ! actually allocatable biomass) to tissue growth. This means: that
          ! if allocatable biomass is higher than Rm, we have to give this
          ! surplus back to reserve pool. 

          !! 3.10 Calculate allocatable part of the labile pool
          !  If there is a plant and not in senescence or dormancy phase,
          !  we calculate labile pool use for growth.
          IF (SUM(circ_class_n(ipts,j,:)) .GT. min_stomate .AND. &
               ( plant_status(ipts,j) .EQ. ibudbreak .OR. &
                 plant_status(ipts,j) .EQ. icanopy .OR. & 
                 plant_status(ipts,j) .EQ. ipresenescence) .AND. &
               SUM(circ_class_biomass(ipts,j,:,ileaf,icarbon)) .GT. min_stomate ) THEN

             ! Use carbon from the labile pool to allocate. The allometric (or 
             ! functional) allocation scheme transfers gpp to the labile pool 
             ! (see above) and then uses the labile pool (gpp + labile(t-1)) to sustain 
             ! growth. The fraction of the labile pool that can be used is a 
             ! function is given by gtemp (see above). bm_alloc_tot is in 
             ! gC m-2 dt-1
             bm_alloc_tot(ipts,j) = gtemp(ipts,j)*tmp_bm(ipts,j,ilabile,icarbon)

             ! Avoid issues with small estimates for bm_alloc_tot. Such small 
             ! number issues could result in mass balance problems.
             IF (bm_alloc_tot(ipts,j) .LE. min_stomate) THEN
                
                ! Not enough C to calculate the allocation. Keep this carbon 
                ! in the labile pool and try again later
                bm_alloc_tot(ipts,j) = zero

             END IF
             
             ! Update the labile carbon pool
             tmp_bm(ipts,j,ilabile,icarbon) = tmp_bm(ipts,j,ilabile,icarbon) - &
                  bm_alloc_tot(ipts,j)
             
          ELSE

             ! The conditions do not support growth
             bm_alloc_tot(ipts,j) = zero

          ENDIF

          ! Debug
          IF(printlev_loc.GE.3 .AND. ipts == test_grid .AND. j == test_pft)THEN 
             WRITE(numout,*) "First bm_alloc_tot ", bm_alloc_tot(ipts,j)
             WRITE(numout,*) "plant_status(ipts,j) ", plant_status(ipts,j)
             WRITE(numout,*) "gtemp ", gtemp(ipts,j)
             WRITE(numout,*) "biomass(ipts,j,ilabile,icarbon) ", &
                  tmp_bm(ipts,j,ilabile,icarbon)
             WRITE(numout,*) "biomass(ipts,j,icarbres,icarbon) ", &
                  tmp_bm(ipts,j,icarbres,icarbon)
          ENDIF
          !-

          !! 3.11 Maintenance respiration
          !  First, total maintenance respiration for the whole plant is 
          !  calculated by summing maintenance respiration of the 
          !  different plant compartments. This simply recalculates the 
          !  maintenance respiration from stomate_resp.f90. Maintenance 
          !  respiration of the different plant parts is calculated in 
          !  stomate_resp.f90 as a function of the plant's temperature, 
          !  the long term temperature and plant coefficients:
          !  The unit of ::resp_maint is gC m-2 dt-1
          resp_maint(ipts,j) = resp_maint(ipts,j) + &
               SUM(resp_maint_part(ipts,j,:))

          ! Following the calculation of hourly maintenance respiration, 
          ! verify that the PFT has not been killed after calcul of 
          ! resp_maint_part in stomate. Can this generaly calculated 
          ! ::resp_maint be use under the given conditions? Surpress 
          ! the respiration for deciduous PFTs as long as they haven't 
          ! carried leaves at least once. When starting from scratch 
          ! there is no budburst in the first year because the longterm 
          ! phenological parameters are not initialized yet. If not 
          ! surpressed respiration consumes all the reserves before the
          ! PFT can start growing. The code would establish a new PFT
          ! but it was decided to surpress this respiration because  
          ! it has no physiological bases.
          IF (SUM(circ_class_n(ipts,j,:)) .GT. min_stomate .AND. &
               rue_longterm(ipts,j) .NE. un) THEN

             !+++CHECK+++
             ! Can the calculated maintenance respiration be used ? Or 
             ! does it have to be adjusted for special cases. Maintenance 
             ! respiration should be positive. In case it is very low, use 20%
             ! (::maint_from_labile) of the active labile carbon pool 
             ! (gC m-2 dt-1)
             ! resp_maint(ipts,j) = MAX(zero, MAX(maint_from_labile * gtemp * 
             ! tmp_bm(ipts,j,ilabile,icarbon), resp_maint(ipts,j)))
             ! Calculate resp_maint for the labile pool as well, 
             ! no need to have the above threshold. Make sure resp_maint 
             ! is not zero
             resp_maint(ipts,j) = MAX(zero, resp_maint(ipts,j))
             !+++++++++++

          ELSE

             ! No plants, no respiration
             resp_maint(ipts,j) = zero

          ENDIF
          
          ! The calculation of ::resp_maint is solely based on the demand i.e.
          ! given the biomass and the condition of the plant, how much should be
          ! respired. It is not sure that this demand can be satisfied i.e. the 
          ! calculated maintenance respiration may exceed the available carbon.
          IF ( bm_alloc_tot(ipts,j) - resp_maint(ipts,j) .LT. zero ) THEN

             IF (plant_status(ipts,j) .EQ. isenescent .OR. & 
                 plant_status(ipts,j) .EQ. idormant .OR. &
                 plant_status(ipts,j) .EQ. ibudsavail) THEN

                ! Under these conditions, bm_alloc_tot will be zero.
                ! this line essentially sets resp_maint as zero during these
                ! stages. This is to not lose the accumulated reserve during
                ! active growing phase.
                resp_maint(ipts,j) = bm_alloc_tot(ipts,j)

             ELSE

                ! the deficit in Rm will be paid by reserve when plant is in stages
                ! of ibudbreak, icanopy, and ipresenescence.
                deficit = bm_alloc_tot(ipts,j) - resp_maint(ipts,j)
                ! The deficit is less than the carbon reserve
                IF (-deficit .LE. tmp_bm(ipts,j,icarbres,icarbon)) THEN
                   
                   ! Pay the deficit from the reserve pool
                   tmp_bm(ipts,j,icarbres,icarbon) = &
                        tmp_bm(ipts,j,icarbres,icarbon) + deficit
                   bm_alloc_tot(ipts,j) = bm_alloc_tot(ipts,j) - deficit

                ELSE
                   
                   ! Not enough carbon to pay the deficit, the individual 
                   ! is going to die at the end of this day
                   bm_alloc_tot(ipts,j) = bm_alloc_tot(ipts,j) + &
                        tmp_bm(ipts,j,icarbres,icarbon) 
                   tmp_bm(ipts,j,icarbres,icarbon) = zero
                   
                   ! Truncate the maintenance respiration to the available carbon
                   resp_maint(ipts,j) = bm_alloc_tot(ipts,j)
                ENDIF
             
             ENDIF 
          ENDIF
          
          ! Final ::resp_maint is known
          bm_alloc_tot(ipts,j) = bm_alloc_tot(ipts,j) - resp_maint(ipts,j)

          ! CYmark: if plant_status is ipresenescence, we don't 
          ! want any allocation to tissue growth. Therefore we put it back to 
          ! reserve pool.
          IF ( plant_status(ipts,j) .EQ. ipresenescence ) THEN
            tmp_bm(ipts,j,icarbres,icarbon) = tmp_bm(ipts,j,icarbres,icarbon) + &
                bm_alloc_tot(ipts,j)
            bm_alloc_tot(ipts,j) = zero
          ENDIF

          ! Debug
          IF(printlev_loc.GE.3 .AND. ipts == test_grid .AND. j == test_pft)THEN
             WRITE(numout,*) 'remaining bm_alloc_tot, ',bm_alloc_tot(ipts,j)
             WRITE(numout,*) "resp_maint ", resp_maint(ipts,j)
          ENDIF
          !-

          !+++CHECK+++
          ! It is more logical to deal with all the respiration terms at the
          ! same time (as being done right here) but there are good reason to calculate
          ! growth respiration in the end. Especialy if in the future we want
          ! to have different growth respiration factors for different tissues.
          
          !  Surpress the respiration for deciduous PFTs as long as they haven't 
          !  carried leaves at least once. If not surpressed respiration consumes 
          !  all the reserves before the PFT can start growing. The code would 
          !  establish a new PFT but it was decided to surpress this respiration 
          !  because it has no physiological bases (in reality the new PFT does 
          !  not start to grow on January 1st as in the model but will be established
          !  at the beginning of the growing season).
          IF (SUM(circ_class_n(ipts,j,:)) .GT. min_stomate .AND. &
               rue_longterm(ipts,j) .NE. un) THEN
             
             frac_growthresp_dyn = frac_growthresp(j)

          ELSE

             frac_growthresp_dyn = zero

          ENDIF

          !! 3.12 Growth respiration
          !  Reserve enough carbon to pay for growth respiration in case all
          !  the available carbon can be allocated. Ideally, growth respiration 
          !  should be included as an additional component in the allocation. 
          !  That way the model could even have different growth respiration 
          !  costs for the different plant organs and/or tissues. The unit of 
          !  resp_growth is gC m-2 dt-1. Calculate resp_growth such that it is 
          !  28% of bm_alloc_tot after resp_growth has been subtracted from 
          !  bm_alloc_tot.
          !  resp_growth = (bm_alloc_tot - resp_growth) * frac_growthresp
          resp_growth(ipts,j) = MAX(zero, bm_alloc_tot(ipts,j)) * &
               (frac_growthresp_dyn / (1 + frac_growthresp_dyn))

          !+++CHECK+++
          ! Set to zero to follow the trunk but it would more straightforward
          ! to delete this parameter from the equations where it is now used but 
          ! set to zero
          frac_growthresp_dyn = 0.
          !+++++++++++

          ! First estimate of ::resp_growth is known. If there is enough
          ! nitrogen to allocate all the C there is no need to recalculate
          ! bm_alloc_tot and resp_growth and so this may be the final
          ! estimate.
          bm_alloc_tot(ipts,j) = bm_alloc_tot(ipts,j) - resp_growth(ipts,j)

          ! Debug
          IF(printlev_loc.GE.3 .AND. ipts == test_grid .AND. j == test_pft)THEN
             WRITE(numout,*) 'remaining bm_alloc_tot, ',bm_alloc_tot(ipts,j)
             WRITE(numout,*) 'resp_growth ', resp_growth(ipts,j)
             IF(bm_alloc_tot(ipts,j).GT.min_stomate)THEN
                WRITE(numout,*) 'ratio resp_growth/bm_alloc_tot, ', &
                     resp_growth(ipts,j)/bm_alloc_tot(ipts,j)
             ENDIF
          ENDIF
          !-

          ! Occasionally, there is a very special situation which arises, where
          ! bm_alloc_tot is greater than min_stomate before accounting for growth 
          ! respiration, but not afterwards.  This causes a mass balance error 
          ! because growth respiration is non-zero but bm_alloc_tot is too small 
          ! to trigger loops below, so nothing is done with that carbon. In this 
          ! situation, the amout of carbon to allocate is so low that nothing 
          ! really changes.  We set the growth respiration to zero in this special 
          ! case to avoid mass imbalance, even though this will not effect the 
          ! trajectory of the plant.  It seems to happen on the same day as leaves 
          ! start growing, before any GPP is calculated.
          IF(((bm_alloc_tot(ipts,j) + resp_growth(ipts,j)) .GT. min_stomate) &
               .AND. (bm_alloc_tot(ipts,j) .LT. min_stomate)) THEN

             bm_alloc_tot(ipts,j) = bm_alloc_tot(ipts,j) + resp_growth(ipts,j)
             resp_growth(ipts,j) = zero

             ! Debug
             IF(printlev_loc.GE.3 .AND. j == test_pft .AND. ipts == test_grid)THEN
                WRITE(numout,*) 'stomate_allocation - Hit exception 25'
                WRITE(numout,*) 'bm_alloc_tot-resp_growth, ',bm_alloc_tot(ipts,j)
                WRITE(numout,*) 'resp_growth ', resp_growth(ipts,j)
                IF (bm_alloc_tot(ipts,j).GT.min_stomate)THEN
                   WRITE(numout,*) 'ratio resp_growth/bm_alloc_tot, ', &
                        resp_growth(ipts,j)/bm_alloc_tot(ipts,j)
                ENDIF
             ENDIF
             !-

          ENDIF

          !! 3.13 Distribute stand level ilabile and icarbres at the tree level
          !  The labile and carbres pools are calculated at the stand level but
          !  are then redistributed at the tree level. Tree level biomass is the
          !  prognostic variable in ORCHIDEE. Biomass is here used as a local
          !  variable to deal with the reserve and labile pools.
          circ_class_biomass(ipts,j,:,ilabile,icarbon) = &
               biomass_to_cc(tmp_bm(ipts,j,ilabile,icarbon),&
               circ_class_biomass(ipts,j,:,ilabile,icarbon),&
               circ_class_n(ipts,j,:))
          circ_class_biomass(ipts,j,:,icarbres,icarbon) = &
               biomass_to_cc(tmp_bm(ipts,j,icarbres,icarbon),&
               circ_class_biomass(ipts,j,:,icarbres,icarbon),&
               circ_class_n(ipts,j,:))

          ! Intermediate mass balance check. Note that this part of
          ! the code is within a DO-loop over nvm so the ipft check
          ! should be used.
          IF (err_act.EQ.4) THEN
             
             ! Reset pool_end
             pool_end(:,:,:) = zero
             
             ! Add bm_alloc_tot into the pool
             pool_end(:,:,icarbon) = pool_end(:,:,icarbon) + &
                  bm_alloc_tot(:,:) * veget_max(:,:)
             
             ! Check mass balance closure. The code above intializes many of the
             ! variables/parameters used in allocation. gpp_daily is allocated
             ! maintenance respiration is calculated and growth respiration is 
             ! reserved. There has been no allocation to leaves, roots and stems yet.
             CALL intermediate_mass_balance_check(pool_start, pool_end, circ_class_biomass, &
                  circ_class_n, veget_max, bm_alloc_tot, gpp_daily, atm_to_bm, dt, npts, &
                  resp_maint, resp_growth, check_intern_init, ipts, j, '1', 'ipft')
          
          ENDIF ! err_act.GT.4

       ENDDO ! npts 

       
 !! 5. Allometric allocation

       DO ipts = 1, npts

          !!  5.1 Initialize allocated biomass pools
          f_alloc(ipts,j,:) = zero
          f_alloc_circ(ipts,:,:) = zero
          Cl_inc(:) = zero
          Cs_inc(:) = zero
          Cr_inc(:) = zero
          Cf_inc(:) = zero
          Cl_incp(:) = zero
          Cs_incp(:) = zero
          Cr_incp(:) = zero 
          Cs_inc_est(:) = zero
          Cl_target(:) = zero
          Cr_target(:) = zero
          Cs_target(:) = zero
          ba1(:) = zero
          ba2(:) = zero
          b_inc_tot = zero

          IF (veget_max(ipts,j) .LE. min_stomate .OR. &
               SUM(circ_class_n(ipts,j,:)) .LE. min_stomate) THEN
             
             ! This vegetation type is not present, so no reason to do the 
             ! calculation. CYCLE will take us out of the innermost DO loop
             CYCLE
             
          ENDIF
 
          !! 5.2 Calculate allocated biomass pools for trees
          !! 5.2.1 Stand to tree allocation rule of Deleuze & Dhote
          IF ( is_tree(j) .AND. bm_alloc_tot(ipts,j) .GT. min_stomate ) THEN

             !  Basal area at the tree level (m2 tree-1)
             circ_class_ba_eff(:) = wood_to_ba_eff(circ_class_biomass(ipts,j,:,:,icarbon),&
                  j, pipe_tune2(ipts,j))
             circ_class_circ_eff(:) = 2 * pi * SQRT(circ_class_ba_eff(:)/pi)

             ! According to equation (-) in Bellasen et al 2010. 
             ! ln(sigmas) = a_sig * ln(circ_med) + b_sig
             ! sigmas = exp(a_sig*log(median(circ_med))+b_sig);
             ! However, in the code (sapiens_forestry.f90) a different expression was used
             ! sigmas = 0.023+0.58*prctile(circ_med,0.05);
             ! Any of these implementations could work but seem to be more suited for 
             ! continues or nearly continuous diameter distributions, say n_circ > 10
             ! For a small number of diameter classes sigma depends on a prescribed
             ! circumference percentile.
             IF (ncirc .GE. 6) THEN          
               
                ! Calculate the median circumference
                DO l = 1,ncirc
                   
                   IF (SUM(circ_class_n(ipts,j,1:l)) .GE. &
                        0.5 * SUM(circ_class_n(ipts,j,:))) THEN
                      
                      median_circ  = circ_class_circ_eff(l) - 5 * min_stomate
                      EXIT
                         
                   ENDIF
                   
                ENDDO
                
                sigma(ipts,j) = deleuze_a(j) + deleuze_b(j) * median_circ
                
             ELSE
                
                ! The X percentile of the trees that will receive the photosynthates
                ! depends on the FM type. In a coppice stand there is a lot of 
                ! competition between the shoots and only the top half of the shoots
                ! will receive GPP, the other half receives only little GPP. This was 
                ! implemnted to get a reasonable diameter growth of coppice stands. 
                ! If deleuze_p is independent from FM, FM strategies with high densities
                ! have very slow diameter growth because the GPP has to be distributed 
                ! over a large number of individuals. 
                IF (forest_managed(ipts,j) == ifm_cop) THEN

                   deleuze_p(j) = deleuze_p_coppice(j)

                ELSEIF (forest_managed(ipts,j) == ifm_none .OR. &
                     forest_managed(ipts,j) == ifm_thin .OR. &
                     forest_managed(ipts,j) == ifm_src .OR. &
                     forest_managed(ipts,j) == ifm_uneven) THEN
                   
                   deleuze_p(j) = deleuze_p_all(j)
                   
                ELSE
                   
                   WRITE(numout, *) 'forest management, ',ipts,j,forest_managed(ipts,j)
                   CALL ipslerr_p (3,'growth_fun_all', &
                        'Forest management strategy does not exist','','')

                ENDIF

                ! Search for the X percentile, where X is given by ::deleuze_p
                ! Substract a very small number (5*min_stomate) just to be sure that 
                ! the circ_class will be corectly accounted for in GE or LE statements
                DO l = 1,ncirc
                      IF (SUM(circ_class_n(ipts,j,1:l)) .GE. deleuze_p(j) * &
                          SUM(circ_class_n(ipts,j,:))) THEN
                         sigma(ipts,j) = circ_class_circ_eff(l) - 5 * min_stomate
                         EXIT
                      ENDIF
                ENDDO
             ENDIF
          
             !! 5.2 Calculate allocated biomass pools for trees
             !  Only possible if there is biomass to allocate
             !  Use sigma and m_dv to calculate a single coefficient that can be
             !  used in the subsequent allocation scheme.

             ! In the original deleuze-dhote equation, basal area increment
             ! linearly increases by the size of trees. But as the diameter and
             ! crown volume increment have saturation points, it can be hypothesized
             ! that basal increment has the saturation point, as well.
             ! Based on this assumption, decreasing power of deleuze-dhote equation
             ! deleuze-dhote equation is implemented here, the simple function of
             ! mean-diameter. The range of delueze_power was set empirically.
             ! Growth diversity between size classes can be highly sensitive to 
             ! deleuze_power_a, which determines a degree of decrease of power
             ! of deleuze_dhote eq. If deleuze_power_a the equation work as 
             ! its original but if it is bigger than 0 growth diversity will decrease. 
             d_mean = 0
             DO l = 1,ncirc
                d_mean = d_mean + ((circ_class_circ_eff(l)/pi)*circ_class_n(ipts,j,l))
             ENDDO
             d_mean = d_mean/SUM(circ_class_n(ipts,j,:))

             deleuze_power = 1.8 + deleuze_power_a(j)*d_mean

             IF (deleuze_power .LE. 2.0) THEN
                deleuze_power = 2.0
             ELSEIF (deleuze_power .GE. 3.5) THEN
                deleuze_power = 3.5
             ENDIF

             circ_class_dba(:) = (circ_class_circ_eff(:) - m_dv(j)*sigma(ipts,j) + &
                  ((m_dv(j)*sigma(ipts,j) + circ_class_circ_eff(:))**2 - &
                  (4*sigma(ipts,j)*circ_class_circ_eff(:)))**(1/deleuze_power))/ 2

             IF (printlev_loc >= 4 .AND. j==test_pft .AND. ipts==test_grid) THEN
                WRITE(numout,*) 'circ_class_dba(:) =', circ_class_dba(:) 
                WRITE(numout,*) 'sigma =', sigma(test_grid,test_pft) 
                WRITE(numout,*) 'median_circ =', median_circ 
             ENDIF    
      
             !! 5.2.2 Scaling factor to convert variables to the individual plant
             !  Allocation is on an individual basis. Stand-level variables need to 
             !  convert to a single individual. Different approach between the DGVM 
             !  and statitic approach
             IF (ok_dgvm) THEN

                ! The DGVM does currently NOT work with the new allocation, consider this as
                ! placeholder. The original code had two different transformations to 
                ! calculate the scalars. Both could be used but the units will differ.
                ! When fixing the DGVM check which quantities need to be multiplied by scal 
                ! scal = ind(ipts,j) * cn_ind(ipts,j) / veget_max(ipts,j)
                scal(ipts,j) = veget_max(ipts,j) / SUM(circ_class_n(ipts,j,:))  
             
             ELSE
             
                ! circ_class_biomass contain the data at the tree level
                ! no conversion required
                scal(ipts,j) = 1.
             
             ENDIF       

             !! 5.2.3 Current biomass pools per tree (gC tree^-1) 
             ! We will have different trees so this has to be calculated from the 
             ! diameter relationships            
             Cs(:) = ( circ_class_biomass(ipts,j,:,isapabove,icarbon) + &
                  circ_class_biomass(ipts,j,:,isapbelow,icarbon) ) * scal(ipts,j)
             Cr(:) = circ_class_biomass(ipts,j,:,iroot,icarbon) * scal(ipts,j)
             Cl(:) = circ_class_biomass(ipts,j,:,ileaf,icarbon) * scal(ipts,j)
             Ch(:) = ( circ_class_biomass(ipts,j,:,iheartabove,icarbon) + &
                  circ_class_biomass(ipts,j,:,iheartbelow,icarbon) ) * scal(ipts,j)

             ! Make a crude estimate of how much carbon can be allocated given 
             ! the available nitrogen. The same code as in the section 5.4 is used exept
             ! that we don't use allocation coeff to modulate n_avail. So
             ! costf=1. It is a strong assumption compared to previous versions. 
             ! It means that ordinary allocation can only happens when phenological
             ! allocation is ok. In other case no wood growth is allowed.
             ! In the case of a strong limitation by Nitrogen, the growth period
             ! for sapwood will be shorten because we reach allometry late in
             ! the growing season.
             n_avail = MAX(tmp_bm(ipts,j,ilabile,initrogen)*0.9,0.0)

             ! Calculate how much carbon could be allocated with the available nitrogen
             bm_supply_n = n_avail  / (1.-frac_growthresp_dyn) * &
                  cn_leaf(ipts,j)

             ! If there is not enough nitrogen, move nitrogen from the reserve
             ! as much as needed, keeping 10% of reserve (arbitral portion)
             IF(bm_alloc_tot(ipts,j) .GT. bm_supply_n &
                  .AND. n_avail .GT. zero) THEN

                ! Calculate the deficit
                n_deficit = bm_alloc_tot(ipts,j) * (1.-frac_growthresp_dyn) / &
                     cn_leaf(ipts,j) - n_avail

                IF(n_deficit .LE. tmp_bm(ipts,j,icarbres,initrogen) * 0.9) THEN

                   ! Enougn N in the reserve pools to fill the labile pool
                   n_avail = n_avail + n_deficit
                   bm_supply_n = n_avail / (1.-frac_growthresp_dyn) * &
                        cn_leaf(ipts,j)
                   tmp_bm(ipts,j,icarbres,initrogen) = tmp_bm(ipts,j,icarbres,initrogen) - &
                        (n_avail/0.9 - tmp_bm(ipts,j,ilabile,initrogen))
                   tmp_bm(ipts,j,ilabile,initrogen) = n_avail/0.9

                   ! tmp_bm is a temporary varaiable so the prognostic variable, i.e.,
                   ! circ_class_biomass also needs to be updated.
                   circ_class_biomass(ipts,j,:,icarbres,initrogen) = &
                        biomass_to_cc(tmp_bm(ipts,j,icarbres,initrogen),&
                        circ_class_biomass(ipts,j,:,icarbres,initrogen),&
                        circ_class_n(ipts,j,:))
                   circ_class_biomass(ipts,j,:,ilabile,initrogen) = &
                        biomass_to_cc(tmp_bm(ipts,j,ilabile,initrogen),&
                        circ_class_biomass(ipts,j,:,ilabile,initrogen),&
                        circ_class_n(ipts,j,:))

                ELSE

                   ! Deficit exceeds 90% of reserve. fill labile as much as
                   ! possible
                   tmp_bm(ipts,j,ilabile,initrogen) = tmp_bm(ipts,j,ilabile,initrogen) + &
                        tmp_bm(ipts,j,icarbres,initrogen) * 0.9
                   tmp_bm(ipts,j,icarbres,initrogen) = tmp_bm(ipts,j,icarbres,initrogen) - &
                        tmp_bm(ipts,j,icarbres,initrogen) * 0.9

                   ! tmp_bm is a temporary varaiable so the prognostic variable, i.e.,
                   ! circ_class_biomass also needs to be updated.
                   circ_class_biomass(ipts,j,:,icarbres,initrogen) = &
                        biomass_to_cc(tmp_bm(ipts,j,icarbres,initrogen),&
                        circ_class_biomass(ipts,j,:,icarbres,initrogen),&
                        circ_class_n(ipts,j,:))
                   circ_class_biomass(ipts,j,:,ilabile,initrogen) = &
                        biomass_to_cc(tmp_bm(ipts,j,ilabile,initrogen),&
                        circ_class_biomass(ipts,j,:,ilabile,initrogen),&
                        circ_class_n(ipts,j,:))

                   ! Update the available nitrogen and the carbon that could be allocated
                   ! with that amount of nitrogen
                   n_avail = MAX(tmp_bm(ipts,j,ilabile,initrogen)*0.9,0.0)
                   bm_supply_n = n_avail / (1.-frac_growthresp_dyn) * &
                        cn_leaf(ipts,j)
                ENDIF

             ENDIF
       
             deltacnmax = 1. - exp(-((1.6 * MIN((1./cn_leaf(ipts,j))-&
                          (1./cn_leaf_min_2D(ipts,j)),0.) / &
                          ( (1./(cn_leaf_max_2D(ipts,j))) - &
                          (1./cn_leaf_min_2D(ipts,j)) ) )**4.1))
             
             IF ( bm_alloc_tot(ipts,j) .GT. bm_supply_n ) THEN

                IF (impose_cn) THEN

                   ! Calculate how much nitrogen is missing to allocate all the
                   ! carbon contained in bm_alloc_tot
                   n_deficit = (bm_alloc_tot(ipts,j)-bm_supply_n) * &
                        (1.-frac_growthresp_dyn) / cn_leaf(ipts,j)/0.9
                   
                   ! The nitrogen missing to allocate the entire bm_alloc_tot will be taken
                   ! from the atmosphere and put in the labile pool. 
                   atm_to_bm(ipts,j,initrogen) = atm_to_bm(ipts,j,initrogen) + &
                        n_deficit/dt
                   tmp_bm(ipts,j,ilabile,initrogen) = &
                        tmp_bm(ipts,j,ilabile,initrogen) + n_deficit

                   ! tmp_bm is a temporary varaiable so the prognostic variable, i.e.,
                   ! circ_class_biomass also needs to be updated.
                   circ_class_biomass(ipts,j,:,ilabile,initrogen) = &
                        biomass_to_cc(tmp_bm(ipts,j,ilabile,initrogen),&
                        circ_class_biomass(ipts,j,:,ilabile,initrogen),&
                        circ_class_n(ipts,j,:))

                   ! Estimate the nitrogen pool that is required to allocate all the
                   ! carbon in bm_alloc_tot.
                   n_alloc_tot(ipts,j) = bm_alloc_tot(ipts,j) * &
                        (1.-frac_growthresp_dyn)/cn_leaf(ipts,j)

                ELSE

                   IF (printlev_loc .GE. 4 .AND. ipts.EQ.test_grid .AND. j.EQ.test_pft) THEN
                      WRITE(numout,*) 'N-limitation before allocation'
                   ENDIF 
                   deltacnmax = Dmax * (1.-deltacnmax) 
                   deltacn = n_avail /  ( bm_alloc_tot(ipts,j) * &
                        (1.-frac_growthresp_dyn) * 1./cn_leaf(ipts,j) )
                   deltacn = MIN(MAX(deltacn,1.0-deltacnmax),1.0)
                    
                   n_alloc_tot(ipts,j) =  MIN( n_avail , & 
                        bm_alloc_tot(ipts,j) * (1.-frac_growthresp_dyn) * &
                        MAX(MIN( 1./cn_leaf(ipts,j)*deltacn, 1./cn_leaf_min_2D(ipts,j)), &
                        1./cn_leaf_max_2D(ipts,j)) ) 
                   
                   tmp_bm(ipts,j,ilabile,icarbon) = &
                        tmp_bm(ipts,j,ilabile,icarbon) + &
                        bm_alloc_tot(ipts,j)
                   
                   bm_alloc_tot(ipts,j) = MIN( bm_alloc_tot(ipts,j) , &
                        n_alloc_tot(ipts,j) / (1.-frac_growthresp_dyn) / &
                        MAX(MIN(1./cn_leaf(ipts,j)*deltacn, &
                        1./cn_leaf_min_2D(ipts,j)), 1./cn_leaf_max_2D(ipts,j)) )
                   
                   tmp_bm(ipts,j,ilabile,icarbon) = tmp_bm(ipts,j,ilabile,icarbon) - &
                        bm_alloc_tot(ipts,j) 

                ENDIF ! if impose_cn

             ELSE

                !+++CHECK+++
                ! This code was copied from below but think it could be deleted.
                ! If there is enough N bm_alloc_tot(ipts,j) does not need to be 
                ! updated. Nothing is done with these calculations.
                IF (printlev_loc .GE. 4 .AND. ipts.EQ.test_grid .AND. j.EQ.test_pft) THEN
                      WRITE(numout,*) 'Sufficient nitrogen before allocation'
                ENDIF
               
                deltacnmax=Dmax * deltacnmax
                deltacn = n_avail / ( bm_alloc_tot(ipts,j) * &
                     (1.-frac_growthresp_dyn) * 1./cn_leaf(ipts,j) )
                deltacn=MIN(MAX(deltacn,1.0),1.+deltacnmax)

                n_alloc_tot(ipts,j) =  MIN( n_avail , & 
                     bm_alloc_tot(ipts,j) * (1.-frac_growthresp_dyn) * &
                     MAX(MIN(1./cn_leaf(ipts,j)*deltacn, & 
                     1./cn_leaf_min_2D(ipts,j)),1./cn_leaf_max_2D(ipts,j)) )
                !+++++++++++
              
             ENDIF

             ! Total amount of carbon that needs to ba allocated (::bm_alloc_tot). 
             ! bm_alloc_tot is in gC m-2 day-1. At 1 m2 there are ::ind number of 
             ! trees. We calculate the allocation for ::ncirc trees. Hence b_inc_tot 
             ! needs to be scaled in the allocation routines. For all cases were 
             ! allocation takes place for a single circumference class, scaling 
             ! could be done before the allocation. In the ordinary allocation 
             ! allocation takes place to all circumference classes at the same time. 
             ! Hence scaling takes place in that step for consistency we scale during 
             ! allocation. Note that b_inc (the carbon allocated to an individual 
             ! circumference class cannot be estimates at this point.
             IF (bm_alloc_tot(ipts,j).GT.min_stomate) THEN

                ! There is enough carbon to allocate
                b_inc_tot = bm_alloc_tot(ipts,j)

             ELSE

                ! There is so little carbon that it is not worth the hassle
                ! to allocate. Allocating very small amounts increases the
                ! risk to run into precision errors.
                tmp_bm(ipts,j,ilabile,icarbon) = tmp_bm(ipts,j,ilabile,icarbon) + &
                     bm_alloc_tot(ipts,j)
                b_inc_tot = zero

             ENDIF

              ! Labile carbon is updated in consequence
             circ_class_biomass(ipts,j,:,ilabile,icarbon) = &
                  biomass_to_cc(tmp_bm(ipts,j,ilabile,icarbon),&
                  circ_class_biomass(ipts,j,:,ilabile,icarbon),&
                  circ_class_n(ipts,j,:))

          END IF

          ! Intermediate mass balance check. Note that this part of
          ! the code is in DO-loops over nvm and npts so the
          ! 'ipts' label is used in the mass balance check
          IF(err_act.EQ.4 .AND.is_tree(j)) THEN

             ! Reset pool_end
             pool_end(:,:,:) = zero
             
             ! Add bm_alloc_tot into the pool
             pool_end(ipts,j,icarbon) = pool_end(ipts,j,icarbon) + &
                  b_inc_tot * veget_max(ipts,j)

             ! Check mass balance closure. Between intermediate check 1 and 2a
             ! bm_inc_tot was recalculated by accounting for the available nitrogen
             CALL intermediate_mass_balance_check(pool_start, pool_end, circ_class_biomass, &
                  circ_class_n, veget_max, bm_alloc_tot, gpp_daily, atm_to_bm, dt, npts, &
                  resp_maint, resp_growth, check_intern_init, ipts, j, '2a', 'ipft')
             
          END IF ! err_act.EQ.4

          ! The initial estimate of bm_alloc_tot was high enough to consider allocation
          ! but after accounting for the available nitrogen bm_alloc_tot may have
          ! dropped below the min_stomate threshold so it needs to be tested again
          IF ( is_tree(j) .AND. bm_alloc_tot(ipts,j) .GT. min_stomate ) THEN
             
             !! 5.2.4 C-allocation for trees
             !  The mass conservation equations are detailed in the header of this subroutine.
             !  The scheme assumes a functional relationships between leaves, sapwood and 
             !  roots. When carbon is added to the leaf biomass pool, an increase in the root
             !  biomass is to be expected to sustain water transport from the roots to the 
             !  leaves. Also sapwood is needed to sustain this water transport and to support
             !  the leaves.
             DO l = 1,ncirc 

                !! 5.2.4.1 Calculate tree height
                circ_class_height_eff(l) = pipe_tune2(ipts,j)* & 
                     (4/pi*circ_class_ba_eff(l))**(pipe_tune3(j)/2)

                !! 5.2.4.2 Do the biomass pools respect the pipe model?
                !  Do the current leaf, sapwood and root components respect the allometric 
                !  constraints? Due to plant phenology it is possible that we have too much 
                !  sapwood compared to the leaf and root mass (i.e. in early spring). 
                !  Calculate the optimal root and leaf mass, given the current wood mass 
                !  by using the basic allometric relationships. Calculate the optimal sapwood
                !  mass as a function of the current leaf and root mass.
                Cl_target(l) = MAX( KF(ipts,j) * Cs(l) / circ_class_height_eff(l), &
                     Cr(l) * LF(ipts,j) , Cl(l))
                Cr_target(l) = MAX( Cl_target(l) / LF(ipts,j), &
                     Cs(l) * KF(ipts,j) / LF(ipts,j) / circ_class_height_eff(l) , Cr(l))
                Cs_target(l) = MAX( Cl(l) / KF(ipts,j) * circ_class_height_eff(l), &
                     Cr(l) * LF(ipts,j) / KF(ipts,j) * circ_class_height_eff(l) , Cs(l))

                ! Debug
                IF (printlev_loc.GE.4 .AND. ipts.EQ.test_grid .AND. j.EQ.test_pft) THEN
                   WRITE(numout,*) 'bm_alloc_tot, ', bm_alloc_tot(ipts,j)
                   WRITE(numout,*) 'Does the tree need reshaping? Class: ',l
                   WRITE(numout,*) 'circ_class_height_eff, ', circ_class_height_eff(l)
                   WRITE(numout,*) 'KF, LF, ', KF(ipts,j), LF(ipts,j)
                   WRITE(numout,*) 'Cl_target-Cl, ', Cl_target(l)-Cl(l), Cl_target(l), Cl(l)
                   WRITE(numout,*) 'Cs_target-Cs, ', Cs_target(l)-Cs(l), Cs_target(l), Cs(l)
                   WRITE(numout,*) 'Cr_target-Cr, ', Cr_target(l)-Cr(l), Cr_target(l), Cr(l)
                ENDIF

                !-

                !! 5.2.4.2 Check dimensions of the trees
                ! If Cs = Cs_target then ba and height are correct, else calculate 
                ! the correct dimensions
                IF ( Cs_target(l) - Cs(l) .GT. min_stomate ) THEN

                   ! If Cs = Cs_target then dia and height are correct. However,
                   ! if Cl = Cl_target or Cr = Cr_target then dia and height
                   ! need to be re-estimated. Cs_target should satify the relationship 
                   ! Cl/Cs = KF/height where height is a function of Cs_target
                   ! Search Cs needed to sustain the max of Cl or Cr.
                   ! Search max of Cl and Cr first
                   !
                   ! [UPDATE] After the code passes through turnover or mortality
                   ! we may end up in a situation where we have lost more sapwood 
                   ! than leaves and roots (i.e. sapwood turnover). The model would 
                   ! then suggest that at time=t+1 the tree should be smaller than 
                   ! at time=t0.  From a physiological standpoint this is not 
                   ! possible for the heartwood. If we now calculate Cs_target on 
                   ! the basis of Cl or Cr, we find that Cs_target > Cs. The first 
                   ! priority of the allocation scheme will be to allocate C to Cs.
                   ! Because we don't know yet whether the actual Cr or Cl is what 
                   ! drives the need to allocate to Cs, we calculate Cl_target first. 
                   Cl_target(l) = MAX(Cl(l), Cr(l)*LF(ipts,j))

                   ! Debug
                   IF (printlev_loc.GE.4 .AND. ipts.EQ.test_grid .AND. j.EQ.test_pft) THEN
                      WRITE(numout,*) 'Does the tree need reshaping? ipts, class:',ipts, l
                      WRITE(numout,*) 'circ_class_height_eff, ', circ_class_height_eff(l)
                      WRITE(numout,*) 'KF, LF, ', KF(ipts,j), LF(ipts,j)
                      WRITE(numout,*) 'Cl_target-Cl, ', Cl_target(l)-Cl(l), Cl_target(l), Cl(l)
                      WRITE(numout,*) 'Cs, ', Cs(l)
                      WRITE(numout,*) 'Cs_target', Cs_target(l)
                      WRITE(numout,*) 'Cr, ', Cr(l)
                      WRITE(numout,*) 'Ch, ', Ch(l)
                   ENDIF
                   !-

                   ! We now have the Cl_target that we will use to calculate
                   ! Cs_target. Given the allometric relationships we can
                   ! calculate Cs_target as Cl_target*height/KF. 
                   ! height is a function of ba, which in turn is a function 
                   ! of woodmass (Woodmass = Cs+Ch (sapwood+heartwood) ). We 
                   ! therefore substitute the following equations into one another:
                   !
                   ! (1) Cs_target = Cl_target*height/KF
                   ! (2) height = as a function of ba
                   ! (3) ba = as a function of woodmass_ind
                   !
                   ! This gives:
                   !
                   ! (4) Cl_target = (KF*Cs_target)/(pipe_tune2*(Cs_target+Ch)/ &
                   !                 & pi/4)**(pipe_tune3/(2+pipe_tune3))
                   !
                   ! The function newX searches for the value for Cs_target
                   ! that satisfies this equation (4).
                   Cs_target(l) =  newX(KF(ipts,j), Ch(l), pipe_tune2(ipts,j), &
                        pipe_tune3(j), Cl_target(l), Cs_target(l), &
                        tree_ff(j)*pipe_density(j)*pi/4*pipe_tune2(ipts,j), &
                        Cs(l), 2*Cs(l), 2, j, ipts)

                   ! Recalculate height and ba from the correct Cs_target
                   circ_class_height_eff(l) = Cs_target(l)*KF(ipts,j)/Cl_target(l)
                   circ_class_ba_eff(l) = pi/4*(circ_class_height_eff(l)/ & 
                        pipe_tune2(ipts,j))**(2/pipe_tune3(j))
                   Cl_target(l) = KF(ipts,j) * Cs_target(l) / circ_class_height_eff(l)
                   Cr_target(l) = Cl_target(l) / LF(ipts,j)

                ENDIF

                ! Debug
                IF (printlev_loc.GE.4 .AND. ipts.EQ.test_grid .AND. j.EQ.test_pft) THEN
                   WRITE(numout,*) 'Target values were adjusted if needed, ',ipts,j,l
                   WRITE(numout,*) 'height_fin, ba_fin, ', circ_class_height_eff(l), &
                        circ_class_ba_eff(l)
                   WRITE(numout,*) 'Cl_target, Cs_target, Cr_target, ', Cl_target(l), &
                        Cs_target(l), Cr_target(l)
                   WRITE(numout,*) 'New target values'
                   WRITE(numout,*) 'Cl_target-Cl, ', Cl_target(l)-Cl(l), Cl_target(l), Cl(l)
                   WRITE(numout,*) 'Cs_target-Cs, ', Cs_target(l)-Cs(l), Cs_target(l), Cs(l)
                   WRITE(numout,*) 'Cr_target-Cr, ', Cr_target(l)-Cr(l), Cr_target(l), Cr(l)
                ENDIF
                !-

             ENDDO !ncirc

             ! The step estimate is used to linearalize the diameter vs height 
             ! relationship. Use a prior to distribute b_inc_tot over the individual 
             ! trees. The share of the total sapwood mass is used as a prior. 
             ! Subsequently, estimate the change in diameter by assuming all the 
             ! available C for allocation will be used in Cs. Hence, this represents 
             ! the maximum possible diameter increase. It was not tested whether 
             ! this is the best prior but it seems to work OK although it often 
             ! results in very small (1e-8) negative values, with even more rare 
             ! 1e-6 negative values. A C-balance closure check could reveal 
             ! whether this is a real issue and requires to change the prior or not.
             ! Calculate the linear slope (::s) of the relationship between ba and h as 
             ! (1) s = (ba2-ba)/(height2-height). 
             ! The goal is to approximate the ba2 that is predicted through the 
             ! non-linear ordinary allocation approach, as this will keep the 
             ! trees in allometric balance. In the next time step, allometric 
             ! balance is recalculated and can be corrected through the so-called 
             ! phenological growth; hence, small deviations resulting from the 
             ! linearization will not accumulate with time.
             ! Note that ba2 = ba + delta_ba and that height and ba are related as
             ! (2) height = k2*(4*ba/pi)**(k3/2)
             ! At this stage the only information we have is that there is b_inc_tot 
             ! (gC m-2) available for allocation. There are two obvious approximations 
             ! both making use of the same assumption, i.e. that for the initial 
             ! estimate of delta_ba height is constant. The first approximation is 
             ! crude and assumes that all the available C is used in Cs_inc 
             ! (thus Cs_inc = b_inc_tot / ind ). The second approximation,
             ! implemented here, makes use of the allometric rules and thus accounts 
             ! for the knowledge that allocating one unit the sapwood comes with a cost 
             ! in leaves and roots thus:
             ! b_inc_temp = Cs_inc+Cl_inc+Cr_inc
             ! (3) <=> b_inc_temp ~= (Cs_inc_est+Cs) + KF*(Cs_inc_est+Cs)/H + ...
             !    KF/LF*(Cs_inc_est+Cs)/H - Cs - Cl - Cr
             ! b_inc_temp is the amount of carbon that can be allocated to each diameter 
             ! class. However, only the total amount i.e. b_inc_tot is known. Total 
             ! allocatable carbon is distributed over the different diameter classes 
             ! proportional to their share of the total wood biomass. Divide by 
             ! circ_class_n to get the correct units (gC tree-1)
             ! (4)  b_inc_temp ~= b_inc_tot / circ_class_n * (circ_class_n * ba**(1+k3)) / ...
             !    sum(circ_class_n * ba**(1+k3))
             ! By substituting (4) in (3) an expression is obtained to approximate the 
             ! carbon that will be allocated to sapwood growth per diameter class 
             ! ::Cs_inc_est. This estimate is then used to calculate delta_ba 
             ! (called ::step) 
             ! step = (Cs+Ch+Cs_inc_set)/(tree_ff*pipe_density*height) - ba 
             ! where height is calculated from (2) after replacing ba by ba+delta_ba

             ! Keep it simple - as described in the documentation
             ! Cs_inc_est(:) = ( b_inc_tot / circ_class_n(ipts,j,:) * &
             !     (circ_class_n(ipts,j,:) * circ_class_ba_eff(:)**(un+pipe_tune3(j))) / &
             !     (SUM(circ_class_n(ipts,j,:) * circ_class_ba_eff(:)**(un+pipe_tune3(j)))))

             ! We implemented a more precise approach following the same principles
             ! Guillaume M. -- The guard at the top of this loop tests the SUM of circ_class_n,
             ! so a single empty class used to reach this element-wise division: Inf * 0 = NaN,
             ! silently poisoning gammas and every class of the PFT. An empty class carries a
             ! zero weight (circ_class_n * ba**(1+k3)) anyway, so dividing it by one instead of
             ! zero leaves every populated class bit-identical and the empty class finite.
             circ_class_n_safe(:) = circ_class_n(ipts,j,:)
             WHERE (circ_class_n_safe(:) .LE. zero) circ_class_n_safe(:) = un
             Cs_inc_est(:) = ( b_inc_tot / circ_class_n_safe(:) * &
                  (circ_class_n(ipts,j,:) * circ_class_ba_eff(:)**(un+pipe_tune3(j))) / &
                  (SUM(circ_class_n(ipts,j,:) * circ_class_ba_eff(:)**(un+pipe_tune3(j)))) + &
                  Cs(:) + Cl(:) + Cr(:)) * circ_class_height_eff(:) / &
                  (circ_class_height_eff(:) + KF(ipts,j) + KF(ipts,j)/LF(ipts,j)) - Cs(:)
             step(:) = ((Ch(:)+Cs(:)+Cs_inc_est(:)) / (tree_ff(j)*pipe_density(j)* &
                  circ_class_height_eff(:))) - circ_class_ba_eff(:)

             ! Debug
             IF (printlev_loc.GE.4 .AND. ipts.EQ.test_grid .AND. j.EQ.test_pft) THEN
                WRITE(numout,*) 'ipts, j, ', ipts, j
                WRITE(numout,*) 'initial guess for step, ', step(:)
             END IF

             ! It can happen that step is equal to zero sometimes.  I'm not sure why, but 
             ! there was a case where it was nonzero for circ classes 1 and 3, and zero 
             ! for 2.  This causes s to be zero and provokes a divide by zero error later 
             ! on.  What if we make it not zero?  This might cause a small mass balance
             ! error for this timestep, but I would rather have that than getting an 
             ! infinite biomass, which is what happened in the other case.  These limits 
             ! are arbitrary and adjusted by hand.  If the output file doesn't show this 
             ! warning very often, I think we're okay, since the amount of carbon is really
             ! small.
             DO l=1,ncirc
                IF(step(l) .LT. min_stomate*0.01 .AND. step(l) .GE. zero)THEN
                   step(l)=min_stomate*0.02
                   IF (printlev_loc.GE.4) THEN
                      WRITE(numout,*) 'WARNING: Might cause mass balance problems '//&
                           'in fun_all, position 1'
                      WRITE(numout,*) 'WARNING: ipts,j ',ipts,j
                   END IF
                ELSEIF(step(l) .GT. -min_stomate*0.01 .AND. step(l) .LT. zero)THEN
                   step(l)=-min_stomate*0.02
                   IF (printlev_loc.GE.4) THEN
                      WRITE(numout,*) 'WARNING: Might cause mass balance problems '//&
                           'in fun_all, position 2'
                      WRITE(numout,*) 'WARNING: ipts,j ',ipts,j
                   END IF
                ENDIF
             ENDDO
             s(:) = step(:)/(pipe_tune2(ipts,j)*(4.0_r_std/pi*(circ_class_ba_eff(:)+step(:)))**&
                  (pipe_tune3(j)/deux) - &
                  pipe_tune2(ipts,j)*(4.0_r_std/pi*circ_class_ba_eff(:))**(pipe_tune3(j)/deux))
             
             ! Debug
             IF (printlev_loc.GE.4 .AND. ipts.EQ.test_grid .AND. j.EQ.test_pft) THEN
                WRITE(numout,*) 'ipts, j, ', ipts, j
                WRITE(numout,*) 'final value for step, ', step(:)
                WRITE(numout,*) 's, ', s(:)
             END IF

             !! 5.2.4.3 Phenological growth
             !  Phenological growth and reshaping of the tree in line with the pipe model. 
             !  Turnover removes C from the different plant components but at a 
             !  component-specific rate, as such the allometric constraints are distorted 
             !  at every time step and should be restored before ordinary growth can 
             !  take place
             l = ncirc
             DO WHILE ((l .GT. zero) .AND. (b_inc_tot .GT. min_stomate))

                !! 5.2.4.3.1 The available wood can sustain the available leaves and roots
                !  Calculate whether the wood is in allometric balance. The target values 
                !  should always be larger than the current pools so the use of ABS is 
                !  redundant but was used to be on the safe side (here and in the rest 
                !  of the module) as it could help to find logical flaws.
                IF ( ABS(Cs_target(l) - Cs(l)) .LT. min_stomate ) THEN

                   ! Use the difference between the target and the actual to
                   ! ensure mass balance closure because l times a values 
                   ! smaller than min_stomate can still add up to a value 
                   ! exceeding min_stomate.
                   Cs_incp(l) = MAX(zero, Cs_target(l) - Cs(l))

                   ! Enough leaves and wood, only grow roots
                   IF ( ABS(Cl_target(l) - Cl(l))  .LT. min_stomate ) THEN

                      ! Allocate at the tree level to restore allometric balance
                      ! Some carbon may have been used for Cs_incp and Cl_incp
                      ! adjust the total allocatable carbon
                      Cl_incp(l) = MAX(zero, Cl_target(l) - Cl(l))
                      Cr_incp(l) = MAX( MIN(b_inc_tot / circ_class_n(ipts,j,l) - &
                           Cs_incp(l) - Cl_incp(l), Cr_target(l) - Cr(l)), zero )

                      ! Write debug comments to output file
                      IF ((j.EQ.test_pft .AND. ipts.EQ.test_grid .AND. printlev_loc.GE.4) .OR. printlev_loc>=5) THEN
                         CALL comment(npts, Cl_target, Cl, Cs_target, Cs, Cr_target, Cr, &
                              delta_ba, ipts, j, l, b_inc_tot, Cl_incp, Cs_incp, Cr_incp, &
                              KF, LF, Cl_inc, Cs_inc, Cr_inc, Cf_inc, grow_wood, &
                              circ_class_n, 1)
                      ENDIF

                   ! Sufficient wood and roots, allocate C to leaves
                   ELSEIF ( ABS(Cr_target(l) - Cr(l)) .LT. min_stomate ) THEN

                      ! Allocate at the tree level to restore allometric balance
                      ! Some carbon may have been used for Cs_incp and Cr_incp
                      ! adjust the total allocatable carbon
                      Cr_incp(l) = MAX(zero, Cr_target(l) - Cr(l))
                      Cl_incp(l) = MAX( MIN(b_inc_tot / circ_class_n(ipts,j,l) - &
                           Cs_incp(l) - Cr_incp(l), Cl_target(l) - Cl(l)), zero )

                      ! Write debug comments to output file
                      IF ((j.EQ.test_pft .AND. ipts.EQ.test_grid .AND. printlev_loc.GE.4) .OR. printlev_loc>=5) THEN
                         CALL comment(npts, Cl_target, Cl, Cs_target, Cs, Cr_target, Cr, &
                              delta_ba, ipts, j, l, b_inc_tot, Cl_incp, Cs_incp, Cr_incp, &
                              KF, LF, Cl_inc, Cs_inc, Cr_inc, Cf_inc, &
                              grow_wood, circ_class_n, 2)
                      ENDIF
                      
                   ! Both leaves and roots are needed to restore the allometric relationships
                   ELSEIF ( ABS(Cl_target(l) - Cl(l)) .GT. min_stomate .AND. &
                        ABS(Cr_target(l) - Cr(l)) .GT. min_stomate ) THEN                 

                      ! Allocate at the tree level to restore allometric balance
                      !  The equations can be rearanged and written as
                      !  (i) b_inc = Cl_inc + Cr_inc
                      !  (ii) Cr_inc = (Cl_inc+Cl)/LF - Cr
                      !  Substitue (ii) in (i) and solve for Cl_inc
                      !  <=> Cl_inc = (LF*(b_inc+Cr)-Cl)/(1+LF)
                      Cl_incp(l) = MIN( ((LF(ipts,j) * ((b_inc_tot/circ_class_n(ipts,j,l) - &
                           Cs_incp(l)) + Cr(l))) - Cl(l)) / & 
                           (1 + LF(ipts,j)), Cl_target(l) - Cl(l) )
                      Cr_incp(l) = MIN ( ((Cl_incp(l) + Cl(l)) / LF(ipts,j)) - Cr(l), &
                           Cr_target(l) - Cr(l))

                      ! The imbalance between Cr and Cl can be so big that (Cl+Cl_inc)/LF 
                      ! is still less then the available root carbon (observed!). This would 
                      ! result in a negative Cr_incp
                      IF ( Cr_incp(l) .LT. zero ) THEN

                         Cl_incp(l) = MIN( b_inc_tot/circ_class_n(ipts,j,l) - Cs_incp(l), &
                              Cl_target(l) - Cl(l) )
                         Cr_incp(l) = (b_inc_tot/circ_class_n(ipts,j,l)) - Cs_incp(l) - &
                              Cl_incp(l)

                      ELSEIF (Cl_incp(l) .LT. zero) THEN

                         Cr_incp(l) = MIN( b_inc_tot/circ_class_n(ipts,j,l) - Cs_incp(l), &
                              Cr_target(l) - Cr(l) )
                         Cl_incp(l) = (b_inc_tot/circ_class_n(ipts,j,l)) - &
                              Cs_incp(l) - Cr_incp(l)

                      ENDIF                          

                      ! Write debug comments to output file
                      IF ((j.EQ.test_pft .AND. ipts.EQ.test_grid .AND. printlev_loc.GE.4) .OR. printlev_loc>=5) THEN
                         CALL comment(npts, Cl_target, Cl, Cs_target, Cs, Cr_target, Cr, &
                              delta_ba, ipts, j, l, b_inc_tot, Cl_incp, Cs_incp, Cr_incp, &
                              KF, LF, Cl_inc, Cs_inc, Cr_inc, Cf_inc, &
                              grow_wood, circ_class_n, 3)
                      ENDIF

                   ELSE

                      WRITE(numout,*) 'Exc 1-3: unexpected exception' 
                      IF(err_act.GT.1)THEN
                         CALL ipslerr_p (3,'growth_fun_all',&
                              'Exc 1-3: unexpected exception','','')
                      ENDIF

                   ENDIF

                !! 5.2.4.3.2 Enough leaves to sustain the wood and roots
                ELSEIF ( ABS(Cl_target(l) - Cl(l)) .LT. min_stomate ) THEN

                   ! Use the difference between the target and the actual to
                   ! ensure mass balance closure because l times a values 
                   ! smaller than min_stomate can still add up to a value 
                   ! exceeding min_stomate.
                   Cl_incp(l) = MAX(zero, Cl_target(l) - Cl(l))

                   ! Enough leaves and wood, only grow roots
                   ! This duplicates Exc 1 and these lines should never be called 
                   IF ( ABS(Cs_target(l) - Cs(l)) .LT. min_stomate ) THEN

                      ! Allocate at the tree level to restore allometric balance
                      ! Some carbon may have been used for Cs_incp and Cl_incp
                      ! adjust the total allocatable carbon
                      Cs_incp(l) = MAX(zero, ABS(Cs_target(l) - Cs(l)))
                      Cr_incp(l) = MAX( MIN(b_inc_tot/circ_class_n(ipts,j,l) - &
                           Cl_incp(l) - Cs_incp(l), Cr_target(l) - Cr(l)), zero )

                      ! Write debug comments to output file
                      IF ((j.EQ.test_pft .AND. ipts.EQ.test_grid .AND. printlev_loc.GE.4) .OR. printlev_loc>=5) THEN
                         CALL comment(npts, Cl_target, Cl, Cs_target, Cs, Cr_target, Cr, &
                              delta_ba, ipts, j, l, b_inc_tot, Cl_incp, Cs_incp, Cr_incp, &
                              KF, LF, Cl_inc, Cs_inc, Cr_inc, Cf_inc, grow_wood, &
                              circ_class_n, 4)
                      ENDIF

                   ! Enough leaves and roots. Need to grow sapwood to support the available 
                   ! canopy and roots
                   ELSEIF ( ABS(Cr_target(l) - Cr(l)) .LT. min_stomate ) THEN

                      ! In truth, there might be a little root carbon to allocate here,
                      ! since min_stomate is not equal to zero.  If there is
                      ! enough of this small carbon in every circ class, and there
                      ! are enough circ classes, ordinary allocation will be skipped
                      ! below and we might try to force allocation, which is silly
                      ! if the different in the root masses is around 1e-8. This
                      ! means we will allocate a tiny amount to the roots to make
                      ! sure they are exactly in balance.                  
                      Cr_incp(l) = MAX(zero, ABS(Cr_target(l) - Cr(l)))
                      Cs_incp(l) = MAX( MIN(b_inc_tot/circ_class_n(ipts,j,l) - &
                           Cl_incp(l) - Cr_incp(l), Cs_target(l) - Cs(l)), zero )

                      ! Write debug comments to output file
                      IF ((j.EQ.test_pft .AND. ipts.EQ.test_grid .AND. printlev_loc.GE.4) .OR. printlev_loc>=5) THEN
                         CALL comment(npts, Cl_target, Cl, Cs_target, Cs, Cr_target, Cr, &
                              delta_ba, ipts, j, l, b_inc_tot, Cl_incp, Cs_incp, Cr_incp, &
                              KF, LF, Cl_inc, Cs_inc, Cr_inc, Cf_inc, &
                              grow_wood, circ_class_n, 5)
                      ENDIF                     

                   ! Need both wood and roots to restore the allometric relationships
                   ELSEIF ( ABS(Cs_target(l) - Cs(l) ) .GT. min_stomate .AND. &
                        ABS(Cr_target(l) - Cr(l)) .GT. min_stomate ) THEN

                      ! circ_class_ba_eff and circ_class_height_eff are already calculated
                      ! for a tree in balance. It would be rather complicated to follow
                      ! the allometric rules for wood allocation (implying changes in height 
                      ! and basal area) because the tree is not in balance yet. First try 
                      ! if we can simply satisfy the allocation needs
                      IF (Cs_target(l) - Cs(l) + Cr_target(l) - Cr(l) .LE. &
                           b_inc_tot/circ_class_n(ipts,j,l) - Cl_incp(l)) THEN
                         
                         Cr_incp(l) = Cr_target(l) - Cr(l)
                         Cs_incp(l) = Cs_target(l) - Cs(l)

                      ! Try to satisfy the need for roots
                      ELSEIF (Cr_target(l) - Cr(l) .LE. b_inc_tot/circ_class_n(ipts,j,l) - &
                           Cl_incp(l)) THEN

                         Cr_incp(l) = Cr_target(l) - Cr(l)
                         Cs_incp(l) = b_inc_tot/circ_class_n(ipts,j,l) - &
                              Cl_incp(l) - Cr_incp(l)
                         
                      ! There is not enough use whatever is available
                      ELSE
                         
                         Cr_incp(l) = b_inc_tot/circ_class_n(ipts,j,l) - Cl_incp(l)
                         Cs_incp(l) = zero
                         
                      ENDIF

                      ! Write debug comments to output file
                      IF ((j.EQ.test_pft .AND. ipts.EQ.test_grid .AND. printlev_loc.GE.4) .OR. printlev_loc>=5) THEN
                         CALL comment(npts, Cl_target, Cl, Cs_target, Cs, Cr_target, Cr, &
                              delta_ba, ipts, j, l, b_inc_tot, Cl_incp, Cs_incp, Cr_incp, &
                              KF, LF, Cl_inc, Cs_inc, Cr_inc, Cf_inc, grow_wood, &
                              circ_class_n, 6)
                      ENDIF   

                   ELSE

                      WRITE(numout,*) 'Exc 4-6: unexpected exception'
                      IF(err_act.GT.1)THEN
                         CALL ipslerr_p (3,'growth_fun_all',&
                              'Exc 4-6: unexpected exception','','')
                      ENDIF
                      
                   ENDIF

                !! 5.2.4.3.3 Enough roots to sustain the wood and leaves
                ELSEIF ( ABS(Cr_target(l) - Cr(l)) .LT. min_stomate ) THEN

                   ! Use the difference between the target and the actual to
                   ! ensure mass balance closure because l times a values 
                   ! smaller than min_stomate can still add up to a value 
                   ! exceeding min_stomate.
                   Cr_incp(l) = MAX(zero, Cr_target(l) - Cr(l))

                   ! Enough roots and wood, only grow leaves
                   ! This duplicates Exc 2 and these lines should thus never be called 
                   IF ( ABS(Cs_target(l) - Cs(l)) .LT. min_stomate ) THEN

                      ! Allocate at the tree level to restore allometric balance
                      Cs_incp(l) = MAX(zero, Cs_target(l) - Cs(l))
                      Cl_incp(l) = MAX( MIN(b_inc_tot/circ_class_n(ipts,j,l) - &
                           Cs_incp(l) - Cr_incp(l), &
                           Cl_target(l) - Cl(l)), zero )

                      ! Write debug comments to output file
                      IF ((j.EQ.test_pft .AND. ipts.EQ.test_grid .AND. printlev_loc.GE.4) .OR. printlev_loc>=5) THEN
                         CALL comment(npts, Cl_target, Cl, Cs_target, Cs, Cr_target, Cr, &
                              delta_ba, ipts, j, l, b_inc_tot, Cl_incp, Cs_incp, Cr_incp, &
                              KF, LF, Cl_inc, Cs_inc, Cr_inc, Cf_inc, grow_wood, &
                              circ_class_n, 7)
                      ENDIF
                   ! Enough leaves and roots. Need to grow sapwood to support the 
                   ! available canopy and roots. Duplicates Exc. 4 and these lines 
                   ! should thus never be called 
                   ELSEIF ( ABS(Cl_target(l) - Cl(l)) .LT. min_stomate ) THEN

                      ! Allocate at the tree level to restore allometric balance
                      Cl_incp(l) = MAX(zero, Cl_target(l) - Cl(l))
                      Cs_incp(l) = MAX( MIN(b_inc_tot/circ_class_n(ipts,j,l) - &
                           Cr_incp(l) - Cl_incp(l), Cs_target(l) - Cs(l) ), zero )

                      ! Write debug comments to output file
                      IF ((j.EQ.test_pft .AND. ipts.EQ.test_grid .AND. printlev_loc.GE.4) .OR. printlev_loc>=5) THEN
                         CALL comment(npts, Cl_target, Cl, Cs_target, Cs, Cr_target, Cr, &
                              delta_ba, ipts, j, l, b_inc_tot, Cl_incp, Cs_incp, Cr_incp, &
                              KF, LF, Cl_inc, Cs_inc, Cr_inc, Cf_inc, grow_wood, &
                              circ_class_n, 8)
                      ENDIF

                   ! Need both wood and leaves to restore the allometric relationships
                   ELSEIF ( ABS(Cs_target(l) - Cs(l)) .GT. min_stomate .AND. &
                        ABS(Cl_target(l) - Cl(l)) .GT. min_stomate ) THEN

                      ! circ_class_ba_eff and circ_class_height_eff are already calculated
                      ! for a tree in balance. It would be rather complicated to follow
                      ! the allometric rules for wood allocation (implying changes in height 
                      ! and basal area) because the tree is not in balance.First try if we 
                      ! can simply satisfy the allocation needs
                      IF (Cs_target(l) - Cs(l) + Cl_target(l) - Cl(l) .LE. &
                           b_inc_tot/circ_class_n(ipts,j,l) - Cr_incp(l)) THEN

                         Cl_incp(l) = Cl_target(l) - Cl(l)
                         Cs_incp(l) = Cs_target(l) - Cs(l)

                      ! Try to satisfy the need for leaves
                      ELSEIF (Cl_target(l) - Cl(l) .LE. b_inc_tot/circ_class_n(ipts,j,l) - &
                           Cr_incp(l)) THEN

                         Cl_incp(l) = Cl_target(l) - Cl(l)
                         Cs_incp(l) = b_inc_tot/circ_class_n(ipts,j,l) - &
                              Cr_incp(l) - Cl_incp(l)

                      ! There is not enough use whatever is available
                      ELSE

                         Cl_incp(l) = b_inc_tot/circ_class_n(ipts,j,l) - Cr_incp(l)
                         Cs_incp(l) = zero

                      ENDIF

                      ! Write debug comments to output file
                      IF ((j.EQ.test_pft .AND. ipts.EQ.test_grid .AND. printlev_loc.GE.4) .OR. printlev_loc>=5) THEN
                         CALL comment(npts, Cl_target, Cl, Cs_target, Cs, Cr_target, Cr, &
                              delta_ba, ipts, j, l, b_inc_tot, Cl_incp, Cs_incp, Cr_incp, &
                              KF, LF, Cl_inc, Cs_inc, Cr_inc, Cf_inc, grow_wood, &
                              circ_class_n, 9)
                      ENDIF

                   ELSE

                      WRITE(numout,*) 'Exc 7-9: unexpected exception'
                      IF(err_act.GT.1)THEN
                         CALL ipslerr_p (3,'growth_fun_all',&
                              'Exc 7-9: unexpected exception','','')
                      ENDIF

                   ENDIF

                
                ELSE

                   ! Either Cl_target, Cs_target or Cr_target should be zero
                   ! Something possibly important was overlooked                   
                   WRITE(numout,*) 'WARNING 4: logical flaw in the phenological '//&
                        'allocation, PFT, class: ', j, l
                   WRITE(numout,*) 'WARNING 4: PFT, ipts: ',j,ipts
                   WRITE(numout,*) 'Cs - Cs_target', Cs(l), Cs_target(l)
                   WRITE(numout,*) 'Cl - Cl_target', Cl(l), Cl_target(l)
                   WRITE(numout,*) 'Cr - Cr_target', Cr(l), Cr_target(l)
                   IF(err_act.GT.1)THEN
                       CALL ipslerr_p (3,'growth_fun_all',&
                            'WARNING 4: logical flaw in the phenological allocation','','')
                   ENDIF

                ENDIF

                IF ( Cl_incp(l) .GE. zero .AND. Cr_incp(l) .GE. zero .AND. &
                     Cs_incp(l) .GE. zero) THEN

                   ! Prevent overspending for leaves
                   IF (b_inc_tot - circ_class_n(ipts,j,l) * Cl_incp(l) .LT. zero) THEN
                      Cl_incp(l) = b_inc_tot/circ_class_n(ipts,j,l)
                   ENDIF
                   b_inc_tot = MAX(zero,b_inc_tot - circ_class_n(ipts,j,l) * Cl_incp(l))

                   ! Prevent overspending for roots
                   IF (b_inc_tot - circ_class_n(ipts,j,l) * Cr_incp(l) .LT. zero) THEN  
                      Cr_incp(l) = b_inc_tot/circ_class_n(ipts,j,l)
                   ENDIF
                   b_inc_tot = MAX(zero,b_inc_tot - circ_class_n(ipts,j,l) * Cr_incp(l))

                   ! Prevent overspending for sapwood
                   IF (b_inc_tot - circ_class_n(ipts,j,l) * Cs_incp(l) .LT. zero) THEN  
                      Cs_incp(l) = b_inc_tot/circ_class_n(ipts,j,l)
                   ENDIF
                   b_inc_tot = MAX(zero,b_inc_tot - circ_class_n(ipts,j,l) * Cs_incp(l))

                   ! Fake allocation for less messy equations in next case, 
                   ! incp needs to be added to inc at the end. 
                   Cl(l) = Cl(l) + Cl_incp(l)
                   Cr(l) = Cr(l) + Cr_incp(l)
                   Cs(l) = Cs(l) + Cs_incp(l)
 
                   IF (b_inc_tot .LT. zero) THEN
                      WRITE(numout,*) 'WARNING 5: numerical problem, '//&
                           'overspending in phenological allocation'
                      WRITE(numout,*) 'WARNING 5: PFT, ipts: ',j,ipts
                      WRITE(numout,*) 'b_inc_tot, ',b_inc_tot
                      WRITE(numout,*) 'Cl_incp(l), Cr_incp(l), Cs_incp(l), ', &
                           l, Cl_incp(l), Cr_incp(l), Cs_incp(l)
                      CALL ipslerr_p (3,'growth_fun_all',&
                           'WARNING 5: numerical problem, overspending in',&
                           'phenological allocation','')
                   ENDIF

                ELSE

                   ! The code was written such that the increment pools should be 
                   ! greater than or equal to zero. If this is not the case, something 
                   ! fundamental is wrong with the if-then constructs under §5.2.4.3
                   WRITE(numout,*) 'WARNING 6: PFT, ipts: ',j,ipts
                   CALL ipslerr_p (3,'growth_fun_all',&
                        'WARNING 6: numerical problem,',&
                        'one of the increment pools is less than zero','')
                ENDIF

                ! Set counter for next circumference class
                l = l-1

             ENDDO ! DO WHILE l.GE.1 .AND. b_inc_tot .GT. min_stomate
            
             ! Intermediate mass balance check. Note that this part of
             ! the code is in DO-loops over nvm and npts so the
             ! 'ipts' label is used in the mass balance check
             IF(err_act.EQ.4) THEN
                
                ! Reset pool_end
                pool_end(:,:,:) = zero
                
                ! Add allocated pools to pool_end
                pool_end(ipts,j,icarbon) = pool_end(ipts,j,icarbon) + &
                     (SUM((Cl_incp(:) + Cs_incp(:) + Cr_incp(:)) * circ_class_n(ipts,j,:)) + &
                     b_inc_tot) * veget_max(ipts,j)
                
                ! Check mass balance closure. Between intermediate check 2a and 2b
                ! phenological allocation was accounted for.
                CALL intermediate_mass_balance_check(pool_start, pool_end, circ_class_biomass, &
                     circ_class_n, veget_max, bm_alloc_tot, gpp_daily, atm_to_bm, dt, npts, &
                     resp_maint, resp_growth, check_intern_init, ipts, j, '2b', 'ipft')
                
             END IF ! err_act.EQ.4

             !! 5.2.4 Record basal area growth during phenological growth 
             ! During phenological growth, some carbon may have been allocated
             ! to the sapwood which then resulted in an increase in basal area
             ! This increase in basal area has not been recorded yet. This
             ! is done below. Later in the code, this increase is then added to
             ! the increase in basal area due to ordinary growth to obtain the 
             ! total increase which is an output variable and is used to 
             ! calculate the tree ring width.
             ba1(:) = wood_to_ba_eff(circ_class_biomass(ipts,j,:,:,icarbon),&
                  j, pipe_tune2(ipts,j))
             temp_mass(:,:) = circ_class_biomass(ipts,j,:,:,icarbon)
             temp_mass(:,isapabove) = temp_mass(:,isapabove) + Cs_incp(:)
             ba2(:) = wood_to_ba_eff(temp_mass,j,pipe_tune2(ipts,j))
             store_delta_ba_eff(ipts,j,:) = ba2(:) - ba1(:)

             !! 5.2.5 Calculate the expected size of the reserve pool
             !  Reserve and labile pools are calculated at the stand level so 
             !  first calculate the stand level biomass from the tree level 
             !  biomass and the number of trees. Note that this value might
             !  be different from the previous values calculated in biomass
             !  because phenological growth could have increased the 
             !  sapwood biomass. 
             tmp_bm(ipts,j,:,:) = cc_to_biomass(ipts,j,&
                  circ_class_biomass(ipts,j,:,:,:),&
                  circ_class_n(ipts,j,:))

             ! use the minimum of either (1) 2% of the total sapwood biomass 
             ! or (2) the amount of carbon needed to develop the optimal LAI 
             ! and the roots. This reserve pool estimate is only used to decide 
             ! whether wood should be grown or not. When really dealing with 
             ! the reserves the reserve pool is recalculated (and the fraction
             ! is 12% rather than the 2% used here). See further below §7.1. 
             reserve_target(ipts,j,icarbon) = &
                  MIN( 0.02 * ( tmp_bm(ipts,j,isapabove,icarbon) + & 
                  tmp_bm(ipts,j,isapbelow,icarbon)), &
                  lai_to_biomass(lai_target(ipts,j),j) * &
                  (1.+root_reserve(j)/ltor(ipts,j)))               
             

             ! If the carbohydrate pool is too small, don't grow wood
             IF ( (pheno_type(j) .NE. 1) .AND. &
                  (tmp_bm(ipts,j,icarbres,icarbon) .LE. reserve_target(ipts,j,icarbon)) ) THEN
                grow_wood = .FALSE.
                ! Debug
                IF(printlev_loc.GE.4 .AND. ipts.EQ.test_grid .AND. j.EQ.test_pft)THEN
                   WRITE(numout,*) 'Not enough carbres to develop the optimal LAI ',j,ipts
                   WRITE(numout,*) 'Reserve pool:',tmp_bm(ipts,j,icarbres,icarbon)
                ENDIF
                !-
             ELSE
                grow_wood = .TRUE.   
             ENDIF

             ! Write debug comments to output file
             IF ((j.EQ.test_pft .AND. ipts.EQ.test_grid .AND. printlev_loc.GE.4) .OR. printlev_loc>=5) THEN
                CALL comment(npts, Cl_target, Cl, Cs_target, Cs, Cr_target, Cr, &
                     delta_ba, ipts, j, l, b_inc_tot, Cl_incp, Cs_incp, Cr_incp, &
                     KF, LF, Cl_inc, Cs_inc, Cr_inc, Cf_inc, grow_wood, &
                     circ_class_n, 10)
             ENDIF

             !! 5.2.6 Ordinary growth
             !  Allometric relationship between components is respected, sustain 
             !  ordinary growth and allocate biomass to leaves, wood, roots and
             !  fruits.
             IF ( (SUM(ABS(Cl_target(:)-Cl(:))) .LE. min_stomate) .AND. &
                  (SUM(ABS(Cr_target(:)-Cr(:))) .LE. min_stomate) .AND. &
                  (SUM(ABS(Cs_target(:)-Cs(:))) .LE. min_stomate) .AND. &
                  grow_wood .AND. b_inc_tot .GT. min_stomate  ) THEN

                ! Allocate fraction of carbon to fruit production (at the tree level)
                Cf_inc(:) = b_inc_tot / SUM(circ_class_n(ipts,j,:)) * fruit_alloc(j)

                ! Residual carbon is allocated to the other components (b_inc_tot is
                ! at the stand level)
                b_inc_tot = b_inc_tot * (un-fruit_alloc(j))

                ! Substitute (7), (8) and (9) in (1) 
                ! b_inc = tree_ff*pipe_density*(ba+circ_class_dba*gammas)*...
                ! (height+(circ_class_dba/s*gammas)) - Cs - Ch + ...
                !    KF*tree_ff*pipe_density*(ba+circ_class_dba*gammas) - ... 
                !    (KF*Ch)/(height+(circ_class_dba/s*gammas)) - Cl + ...
                !    KF/LF*tree_ff*pipe_density*(ba+circ_class_dba*gammas) - ...
                !    (KF*Ch/LF)/(height+(circ_class_dba/s*gammas)) - Cr
                !
                ! b_inc+Cs+Ch+Cl+Cr = tree_ff*pipe_density*(ba+circ_class_dba*gammas)*...
                !    (height+(circ_class_dba/s*gammas))  + ...
                !    KF*tree_ff*pipe_density*(ba+circ_class_dba*gammas) - ...
                !    (KF*Ch)/(height+(circ_class_dba/s*gammas)) + ...
                !    KF/LF*tree_ff*pipe_density*(ba+circ_class_dba*gammas) - ...
                !    (KF*Ch/LF)/(height+(circ_class_dba/s*gammas))
                ! <=> b_inc+Cs+Ch+Cl+Cr = circ_class_dba^2/s*tree_ff*...
                !    pipe_density*gammas^2 + circ_class_dba/s*ba*tree_ff*...
                !    pipe_density*gammas + ...
                !    circ_class_dba*height*tree_ff*pipe_density*gammas + ...
                !    bcirc_class_dba*height*tree_ff*pipe_density - ...
                !    (Ch*KF*s)/(circ_class_dba*gammas+height*s) + ...
                !    circ_class_dba*KF*tree_ff*pipe_density*gammas + ...
                !    ba*KF*tree_ff*pipe_density - ...
                !    (Ch*KF*s)/(LF*(circ_class_dba*gammas+height*s)) + ...
                !    circ_class_dba*KF/LF*tree_ff*pipe_density*gammas + ...
                !    ba*KF/LF*tree_ff*pipe_density
                ! (10) b_inc+Cs+Ch+Cl+Cr = (circ_class_dba^2/s*tree_ff*...
                !    pipe_density)*gammas^2 + ...
                !    (circ_class_dba/s*ba*tree_ff*pipe_density + ...
                !    circ_class_dba*height*tree_ff*pipe_density + ...
                !    circ_class_dba*KF*tree_ff*pipe_density + ...
                !    circ_class_dba*KF/LF*tree_ff*pipe_density)*gammas - ...
                !    (Ch*KF*s)(1+1/LF)/(circ_class_dba*gammas+height*s) + ...
                !    bcirc_class_dba*height*tree_ff*pipe_density + ...
                !    ba*KF*tree_ff*pipe_density + ba*KF/LF*tree_ff*pipe_density
                !
                ! Note that b_inc is not known, only b_inc_tot (= sum(b_inc) is known. 
                ! The above equations are for individual trees, at the stand level we 
                ! have to take the sum over the individuals which is 
                ! equivalant to substituting (10) in (2) 
                ! (11) sum(b_inc) + sum(Cs+Ch+Cl+Cr) = ...
                !    sum(circ_class_dba^2/s*tree_ff*pipe_density) * gammas^2 + ...
                !    sum(circ_class_dba/s*ba*tree_ff*pipe_density + ...
                !    circ_class_dba*height*tree_ff*pipe_density + ...
                !    circ_class_dba*KF*tree_ff*pipe_density + ...
                !    circ_class_dba*KF/LF*tree_ff*pipe_density) * gammas - ...
                !    sum[(Ch*KF*s)(1+1/LF)/(circ_class_dba*gammas+height*s)] + ...
                !    sum(bcirc_class_dba*height*tree_ff*pipe_density + ...
                !    ba*KF*tree_ff*pipe_density + ba*KF/LF*tree_ff*pipe_density)
                !
                ! The term sum[(Ch*KF*s)(1+1/LF)/(circ_class_dba*gammas+height*s)] 
                ! can be approximated by a series expansion
                ! (12) sum((Ch*KF*s)(1+1/LF)/(height*s) + ...
                !    sum((Ch*KF*s)(1+1/LF)*circ_class_dba/(height*s)^2)*gammas + ...
                !    sum((Ch*KF*s)(1+1/LF)*circ_class_dba^2/(height*s)^3)*gammas^2
                !
                ! Substitute (12) in (11)
                ! sum(b_inc) + sum(Cs+Ch+Cl+Cr) = ...
                !    sum(circ_class_dba^2/s*tree_ff*pipe_density - ...
                !    (Ch*KF*s)*(1+1/LF)*circ_class_dba^2/(height*s)^3) * gammas^2 + ...
                !    sum(circ_class_dba/s*ba*tree_ff*pipe_density + ...
                !    circ_class_dba*height*tree_ff*pipe_density + ...
                !    circ_class_dba*KF*tree_ff*pipe_density + ...
                !    circ_class_dba*KF/LF*tree_ff*pipe_density + ...
                !    (Ch*KF*s)*(1+1/LF)*circ_class_dba/(height*s)^2) * gammas + ...
                !    sum(bcirc_class_dba*height*tree_ff*pipe_density + ...
                !    ba*KF*tree_ff*pipe_density + ba*KF/LF*tree_ff*pipe_density - ...
                !    (Ch*KF*s)*(1+1/LF)/(height*s))
                !
                ! Solve this quadratic equation for gammas.
                a = SUM( circ_class_n(ipts,j,:) * &
                     (circ_class_dba(:)**2/s(:)*tree_ff(j)*pipe_density(j) - &
                     (Ch(:)*KF(ipts,j)*s(:))*(1+1/LF(ipts,j))*&
                     (circ_class_dba(:)**2/(circ_class_height_eff(:)*s(:))**3)) )
                b = SUM( circ_class_n(ipts,j,:) * &
                     (circ_class_dba(:)/s(:)*circ_class_ba_eff(:)*tree_ff(j)*pipe_density(j) + &
                     circ_class_dba(:)*circ_class_height_eff(:)*tree_ff(j)*pipe_density(j) + &
                     circ_class_dba(:)*KF(ipts,j)*tree_ff(j)*pipe_density(j) + &
                     circ_class_dba(:)*KF(ipts,j)/LF(ipts,j)*tree_ff(j)*pipe_density(j) + &
                     (Ch(:)*KF(ipts,j)*s(:))*(1+1/LF(ipts,j))*circ_class_dba(:)/&
                     (circ_class_height_eff(:)*s(:))**2) )
                c = SUM( circ_class_n(ipts,j,:) * &
                     (circ_class_ba_eff(:)*circ_class_height_eff(:)*&
                     tree_ff(j)*pipe_density(j) + &
                     circ_class_ba_eff(:)*KF(ipts,j)*tree_ff(j)*pipe_density(j) + &
                     circ_class_ba_eff(:)*KF(ipts,j)/LF(ipts,j)*tree_ff(j)*pipe_density(j) - &
                     (Ch(:)*KF(ipts,j)*s(:))*(1+1/LF(ipts,j))/&
                     (circ_class_height_eff(:)*s(:)) - &
                     (Cs(:) + Ch(:) + Cl(:) + Cr(:))) ) - b_inc_tot

                ! Solve the quadratic equation a*gammas2 + b*gammas + c = 0, for gammas.
                gammas(ipts,j) = (-b + sqrt(b**2-4*a*c)) / (2*a)
                
                ! After thousands of simulation years we had a single pixel where 
                ! some of the three circ_class got a negative growth. This was because
                ! both roots of the quadratic equation were negative. If both roots are 
                ! negative, we don't allocate and simply leave the carbon in the labile 
                ! pool. We will try again with more carbon the next day.
                IF (gammas(ipts,j).LT.zero) THEN

                   ! Move the unallocatable carbon back into the labile
                   ! pool. Update related variables to pass the mass balance
                   ! check. Put the fruit allocation back first. That will give
                   ! more carbon at the next time step.
                   b_inc_tot = b_inc_tot + SUM(Cf_inc(:) * circ_class_n(ipts,j,:))
                   Cf_inc(:) = zero
                   tmp_bm(ipts,j,ilabile,icarbon) = &
                        tmp_bm(ipts,j,ilabile,icarbon) + b_inc_tot
                   bm_alloc_tot(ipts,j) = bm_alloc_tot(ipts,j) - b_inc_tot
 
                   ! Calculate C that was not allocated (b_inc_tot), the 
                   ! equation should read b_inc_tot = b_inc_tot - b_inc_tot 
                   b_inc_tot = zero

                   IF (printlev_loc.GE.4 .AND. ipts.EQ.test_grid .AND. j.EQ.test_pft) THEN
                      WRITE(numout,*) 'Both roots are negative for PFT, ', j
                      WRITE(numout,*) 'bm_alloc_tot, ', bm_alloc_tot(ipts,j)
                   ENDIF
                  
                ELSE

                   ! One gammas is positive. The solution for gammas is now used to 
                   ! calculate delta_ba (eq. 3), delta_height (eq. 6), Cs_inc (eq. 7), 
                   ! Cl_inc (eq. 8) and Cr_inc (eq. 9). See comment on the calculation 
                   ! of delta_height and its implications on numerical consistency at 
                   ! the similar statement in §5.2.4.3.1. 
                   ! Tree rings: delta_ba is a sum of phenological and ordinary growth 
                   ! but for further calculation, re-calculate delta_ba considering 
                   ! only ordinary growth. Note that calculated delta_ba is effectvie 
                   ! basal area increment.
                   delta_ba(:) = circ_class_dba(:) * gammas(ipts,j)                
                   store_delta_ba_eff(ipts,j,:) = store_delta_ba_eff(ipts,j,:) + delta_ba(:)
                   delta_height(:) = delta_ba(:)/s(:)              
                   Cs_inc(:) = tree_ff(j)*pipe_density(j)*(circ_class_ba_eff(:) + &
                        delta_ba(:))*(circ_class_height_eff(:) + &
                        delta_height(:)) - Cs(:) - Ch(:)
                   Cl_inc(:) = KF(ipts,j)*tree_ff(j)*pipe_density(j)*&
                        (circ_class_ba_eff(:)+delta_ba(:)) - &
                        (KF(ipts,j)*Ch(:))/(circ_class_height_eff(:)+delta_height(:)) - Cl(:)
                   Cr_inc(:) = KF(ipts,j)/LF(ipts,j)*tree_ff(j)*pipe_density(j)*&
                        (circ_class_ba_eff(:)+delta_ba(:)) - &
                        (KF(ipts,j)*Ch(:)/LF(ipts,j))/(circ_class_height_eff(:)+&
                        delta_height(:)) - Cr(:)

                   ! Write the initial residual to the history file to check
                   ! whether all goes well (or to see how often and where it goes
                   ! wrong).
                   residual_write(ipts,j) = b_inc_tot - SUM(circ_class_n(ipts,j,:)* &
                        (Cl_inc(:) + Cr_inc(:) + Cs_inc(:)))

                   ! Debug
                   IF (printlev_loc.GE.4 .AND. ipts.EQ.test_grid .AND. j.EQ.test_pft) THEN
                      WRITE(numout,*) 'One gamma is positive, ', j, gammas(ipts,j)
                      WRITE(numout,*) 'delta_ba, ', delta_ba(:)
                      WRITE(numout,*) 'Cl_inc, Cr_inc, Cs_inc, ', &
                           Cl_inc(:), Cr_inc(:), Cs_inc(:)
                   ENDIF
                   !-

                   ! There are two possible problems: (1) one of the Cx_inc is negative or
                   ! (2) all Cx_inc are positive but we are (slightly) overspending.
                   IF (MINVAL(Cs_inc(:)) .LT. zero .OR. MINVAL(Cr_inc(:)) .LT. zero .OR. &
                        MINVAL(Cl_inc(:)) .LT. zero) THEN

                      ! The first (rare) problem we need to catch is when one of the increment 
                      ! pools is negative. This is an undesired outcome (see comment where 
                      ! ::KF_old is calculated in this routine. In that case we write a 
                      ! warning, set all increment pools to zero and try it again at the 
                      ! next time step. A likely cause of this problem is a too large change 
                      ! in KF from one time step to another (note that at this point both
                      ! roots should be positive - so that can no longer be the cause). Try 
                      ! decreasing the acceptable value for an absolute increase in KF.
                      
                      ! Do not allocate - save the carbon for the next time step

                      ! Debug
                      IF(err_act.GT.1) THEN
                         WRITE(numout,*) 'WARNING 10a: numerical problem, '//&
                              'one of the increment pools is less than zero'
                         WRITE(numout,*) 'WARNING 10a: PFT, ipts: ',j,ipts
                         WRITE(numout,*) 'WARNING 10a: Cl_inc(:): ',Cl_inc(:)
                         WRITE(numout,*) 'WARNING 10a: Cr_inc(:): ',Cr_inc(:)
                         WRITE(numout,*) 'WARNING 10a: Cs_inc(:): ',Cs_inc(:)
                         WRITE(numout,*) 'WARNING 10a: We will revert the allocation'
                         WRITE(numout,*) ' and save the carbon for the next day'
                      END IF
                      !-

                      ! Move the unallocatable carbon back into the labile
                      ! pool. Update related variables to pass the mass balance check.
                      b_inc_tot = b_inc_tot + SUM(Cf_inc(:) * circ_class_n(ipts,j,:))
                      bm_alloc_tot(ipts,j) = bm_alloc_tot(ipts,j) - b_inc_tot
                      tmp_bm(ipts,j,ilabile,icarbon) = &
                           tmp_bm(ipts,j,ilabile,icarbon) + b_inc_tot 
                      
                      ! Revert the allocation
                      store_delta_ba_eff(ipts,j,:) = store_delta_ba_eff(ipts,j,:) - delta_ba(:)
                      delta_ba(:) = zero
                      delta_height(:) = zero
                      Cl_inc(:) = zero
                      Cs_inc(:) = zero
                      Cr_inc(:) = zero
                      Cf_inc(:) = zero

                      ! Calculate C that was not allocated (b_inc_tot), the 
                      ! equation should read b_inc_tot = b_inc_tot - b_inc_tot 
                      ! note that Cf_inc was already accounted for.
                      b_inc_tot = zero

                      IF (printlev_loc.GE.4 .AND. ipts.EQ.test_grid .AND. j.EQ.test_pft) THEN
                         WRITE(numout,*) 'Negative increment pools, ', j
                         WRITE(numout,*) 'Allocation was reverted' 
                         WRITE(numout,*) 'bm_alloc_tot, ', bm_alloc_tot(ipts,j)
                      ENDIF
                                          
                   ELSEIF (b_inc_tot - SUM(circ_class_n(ipts,j,:)* &
                        (Cl_inc(:) + Cr_inc(:) + Cs_inc(:))).LT. zero) THEN

                      ! We should only be here if there is a positive root and all
                      ! increment pools are positive. Overspending should thus be
                      ! a numerical issue and is expected to be small (less than 10-8)
                      ! If the residual is larger than expected, reduce gamma a bit and 
                      ! recalculate the allocation. The residual is added then back into
                      ! the labile pool. Do not allocate - save the carbon for the next 
                      ! time step
                      residual10b(ipts,j) = b_inc_tot - SUM(circ_class_n(ipts,j,:) * &
                           (Cl_inc(:) + Cr_inc(:) + Cs_inc(:)))
                      
!!$                      ! Debug
!!$                      IF(err_act.EQ.3) THEN
!!$                         WRITE(numout,*) 'WARNING 10b: numerical problem, '//&
!!$                              'residual is negative we are overspending'
!!$                         WRITE(numout,*) 'WARNING 10b: PFT, ipts: ',j,ipts
!!$                         WRITE(numout,*) 'WARNING 10b: residual, ', residual10b
!!$                         WRITE(numout,*) 'WARNING 10b: allocation is being adjusted'
!!$                      ENDIF
!!$                      !-

                      ! It is considered too large so we try to reduce the residual with 
                      ! some brute force. First revert the previous store_delta_ba_eff 
                      ! else we will double count tree ring width growth.
                      store_delta_ba_eff(ipts,j,:) = store_delta_ba_eff(ipts,j,:)-delta_ba(:)
                      ! Calculate new delta_ba (note the reduction factor 0.99) 
                      ! and update the other variables. The reduction factor of 0.99
                      ! could be too large, In that case we won't allocate.
                      delta_ba(:) = circ_class_dba(:) * gammas(ipts,j) * 0.99
                      store_delta_ba_eff(ipts,j,:) = store_delta_ba_eff(ipts,j,:)+delta_ba(:)
                      delta_height(:) = delta_ba(:)/s(:)              
                      Cs_inc(:) = MAX(zero,tree_ff(j)*pipe_density(j)*(circ_class_ba_eff(:) + &
                           delta_ba(:))*(circ_class_height_eff(:) + &
                           delta_height(:)) - Cs(:) - Ch(:))
                      Cl_inc(:) = MAX(zero,KF(ipts,j)*tree_ff(j)*pipe_density(j)*&
                           (circ_class_ba_eff(:)+delta_ba(:)) - &
                           (KF(ipts,j)*Ch(:)) / &
                           (circ_class_height_eff(:)+delta_height(:)) - Cl(:))
                      Cr_inc = MAX(zero,KF(ipts,j)/LF(ipts,j)*tree_ff(j)*pipe_density(j)*&
                           (circ_class_ba_eff(:)+delta_ba(:)) - &
                           (KF(ipts,j)*Ch(:)/LF(ipts,j))/(circ_class_height_eff(:)+&
                           delta_height(:)) - Cr(:))

                      ! Debug
                      IF (printlev_loc.GE.4 .AND. ipts.EQ.test_grid .AND. j.EQ.test_pft) THEN
                         WRITE(numout,*) 'adjusted after overspending, ', j
                         WRITE(numout,*) 'delta_ba, ', delta_ba(:)
                         WRITE(numout,*) 'Cl_inc, Cr_inc, Cs_inc, ', &
                              Cl_inc(:), Cr_inc(:), Cs_inc(:)
                      ENDIF
                      !-
                      
                      ! Check whether the recalculation worked
                      IF (b_inc_tot - SUM(circ_class_n(ipts,j,:) * &
                        (Cl_inc(:) + Cr_inc(:) + Cs_inc(:))) .LT. zero) THEN
                        
                         ! Debug
                         IF(err_act.GT.1) THEN
                            WRITE(numout,*) 'WARNING 10c: numerical problem, '//&
                                 'residual is still negative we are still overspending'
                            WRITE(numout,*) 'WARNING 10c: PFT, ipts: ',j,ipts
                            WRITE(numout,*) 'WARNING 10c: residual, ',b_inc_tot - &
                                 SUM(circ_class_n(ipts,j,:)* &
                                 (Cl_inc(:) + Cr_inc(:) + Cs_inc(:)))
                            WRITE(numout,*) 'WARNING 10b: allocation is being adjusted'
                         ENDIF
                         !-

                         ! Move the unallocatable carbon back into the labile
                         ! pool. Update related variables to pass the mass balance
                         ! check.
                         b_inc_tot = b_inc_tot + SUM(Cf_inc(:) * circ_class_n(ipts,j,:))
                         tmp_bm(ipts,j,ilabile,icarbon) = &
                              tmp_bm(ipts,j,ilabile,icarbon) + b_inc_tot 
                         bm_alloc_tot(ipts,j) = bm_alloc_tot(ipts,j) - b_inc_tot
                        
                         ! The residual is still negative. We will give up.
                         ! Revert the allocation
                         store_delta_ba_eff(ipts,j,:) = store_delta_ba_eff(ipts,j,:) - delta_ba(:)
                         delta_ba(:) = zero
                         delta_height(:) = zero
                         Cl_inc(:) = zero
                         Cs_inc(:) = zero
                         Cr_inc(:) = zero
                         Cf_inc(:) = zero

                         ! Calculate C that was not allocated (b_inc_tot), the 
                         ! equation should read b_inc_tot = b_inc_tot - b_inc_tot 
                         ! note that Cf_inc was already accounted for.
                         b_inc_tot = zero

                         IF (printlev_loc.GE.4 .AND. ipts.EQ.test_grid .AND. j.EQ.test_pft) THEN
                            WRITE(numout,*) 'Initial solution did not work, ', j
                            WRITE(numout,*) 'Allocation was reverted' 
                            WRITE(numout,*) 'bm_alloc_tot, ', bm_alloc_tot(ipts,j)
                         ENDIF

                      ELSE

                         ! Debug
                         IF (printlev_loc.GE.4 .AND. ipts.EQ.test_grid .AND. j.EQ.test_pft) THEN
                            WRITE(numout,*) 'Initial b_inc_tot, ', j, b_inc_tot
                            WRITE(numout,*) 'Initial tmp_bm, ', tmp_bm(ipts,j,ilabile,icarbon)
                            WRITE(numout,*) 'circ_class_n, ',circ_class_n(ipts,j,:)
                            WRITE(numout,*) 'Cl_inc, Cr_inc, Cs_inc, ', &
                                 Cl_inc(:), Cr_inc(:), Cs_inc(:)
                         ENDIF
                         !-

                         ! The adjustment was succesful. Finish the allocation.
                         ! Reduce b_inc_tot and move the difference back into the 
                         ! labile pool where it comes from. Thanks to the IF we know
                         ! for sure that the new b_inc_tot (calculated on the next 
                         ! line) is positive.
                         b_inc_tot =  b_inc_tot - SUM(circ_class_n(ipts,j,:) * &
                              (Cl_inc(:) + Cr_inc(:) + Cs_inc(:)))
                         tmp_bm(ipts,j,ilabile,icarbon) = tmp_bm(ipts,j,ilabile,icarbon) + &
                              b_inc_tot
 
                         ! Debug
                         IF (printlev_loc.GE.4 .AND. ipts.EQ.test_grid .AND. j.EQ.test_pft) THEN
                            WRITE(numout,*) 'Recalculated b_inc_tot, ', j, b_inc_tot
                            WRITE(numout,*) 'Recalculated tmp_bm, ', tmp_bm(ipts,j,ilabile,icarbon)
                            WRITE(numout,*) 'Initial bm_alloc_tot, ', bm_alloc_tot(ipts,j) 
                         ENDIF
                         !-

                         ! What is left in b_inc_tot was not allocated so adjust
                         ! the total allocation.
                         bm_alloc_tot(ipts,j) = bm_alloc_tot(ipts,j) - b_inc_tot

                         ! There is nothing left to allocate
                         b_inc_tot = zero

                         ! Debug
                         IF (printlev_loc.GE.4 .AND. ipts.EQ.test_grid .AND. j.EQ.test_pft) THEN
                            WRITE(numout,*) 'Recalculated bm_alloc_tot, ', bm_alloc_tot(ipts,j)
                         ENDIF
                         !-

                      END IF

                   ELSE

                      ! All is well, wrap up the allocation
                      ! Wrap-up ordinary growth calculate C that was not allocated, note 
                      ! that Cf_inc was already subtracted. We know from the IF
                      ! loops above that the new b_inc_tot (calculated at the next line
                      ! of code) will be positive
                      b_inc_tot = b_inc_tot - SUM(circ_class_n(ipts,j,:) * &
                           (Cl_inc(:) + Cr_inc(:) + Cs_inc(:)))

                      ! If all went well b_inc_tot should now be very close to zero. Due
                      ! to numerical approximations some C may be left. Whatever is
                      ! left is moved back into the labile pool to conserve mass.
                      tmp_bm(ipts,j,ilabile,icarbon) = tmp_bm(ipts,j,ilabile,icarbon) + &
                           b_inc_tot
                      bm_alloc_tot(ipts,j) = bm_alloc_tot(ipts,j) - b_inc_tot
                      b_inc_tot = zero
                      
                      ! Debug
                      IF ((j.EQ.test_pft .AND. ipts.EQ.test_grid .AND. printlev_loc.GE.4) .OR. printlev_loc>=5) THEN
                         WRITE(numout,*) 'wrap-up ordinary allocation, left b_in_tot, ', &
                              b_inc_tot 
                         WRITE(numout,*) 'a, b, c, gammas, ', a, b, c, gammas(ipts,j)
                         WRITE(numout,*) 'delta_height, ', delta_height(:)
                         CALL comment(npts, Cl_target, Cl, Cs_target, Cs, Cr_target, Cr, &
                              delta_ba, ipts, j, l, b_inc_tot, Cl_incp, Cs_incp, Cr_incp, &
                              KF, LF, Cl_inc, Cs_inc, Cr_inc, Cf_inc, grow_wood, &
                              circ_class_n, 11)
                      ENDIF
                      !-    

                   END IF

                END IF ! postive root for quadratic equation

                ! Intermediate mass balance check. Note that this part of
                ! the code is in DO-loops over nvm and npts so the
                ! 'ipts' label is used in the mass balance check
                IF(err_act.EQ.4) THEN

                   ! Update circ_class_biomass
                   circ_class_biomass(ipts,j,:,ilabile,icarbon) = &
                        biomass_to_cc(tmp_bm(ipts,j,ilabile,icarbon),&
                        circ_class_biomass(ipts,j,:,ilabile,icarbon),&
                        circ_class_n(ipts,j,:))
                   circ_class_biomass(ipts,j,:,icarbres,icarbon) = &
                        biomass_to_cc(tmp_bm(ipts,j,icarbres,icarbon),&
                        circ_class_biomass(ipts,j,:,icarbres,icarbon),&
                        circ_class_n(ipts,j,:))

                   ! All carbon should have been allocated and the remainder was moved
                   ! back into the labile pool. b_inc_tot should be zero. If not, the 
                   ! calculation of pool_end is wrong.
                   IF(ABS(b_inc_tot).GT.min_stomate)THEN
                      WRITE(numout,*) 'b_inc_tot differs from zero, ', ipts,j,b_inc_tot
                      CALL ipslerr_p(3,'stomate_growth_fun_all.f90','intermediate mbcheck 2c',&
                           'b_inc_tot differs from zero','')
                   END IF
                   
                   ! Reset pool_end
                   pool_end(:,:,:) = zero
                 
                   ! Update pool_end
                   pool_end(ipts,j,icarbon) = pool_end(ipts,j,icarbon) + &
                        SUM((Cl_incp(:) + Cs_incp(:) + Cr_incp(:) + &
                        Cl_inc(:) + Cs_inc(:) + Cr_inc(:) + Cf_inc(:)) * circ_class_n(ipts,j,:)) * &
                        veget_max(ipts,j)
                   
                   ! Check mass balance closure. Between intermediate check 2a and 2b
                   ! phenological allocation was accounted for.
                   CALL intermediate_mass_balance_check(pool_start, pool_end, circ_class_biomass, &
                        circ_class_n, veget_max, bm_alloc_tot, gpp_daily, atm_to_bm, dt, npts, &
                        resp_maint, resp_growth, check_intern_init, ipts, j, '2c', 'ipft')
                   
                END IF ! err_act.EQ.4

             !! 5.2.7 Don't grow wood, use C to fill labile pool
             ELSEIF ( ((.NOT. grow_wood) .AND. (b_inc_tot .GT. min_stomate)) ) THEN

                ! Calculate the C that needs to be distributed to the 
                ! labile pool. The fraction is proportional to the ratio 
                ! between the total allocatable biomass and the unallocated 
                ! biomass per tree (b_inc now contains the unallocated 
                ! biomass). At the end of the allocation scheme bm_alloc_tot 
                ! is substracted from the labile biomass pool to update the 
                ! biomass pool (tmp_bm(:,:,ilabile) = tmp_bm(:,:,ilabile) - 
                ! bm_alloc_tot(:,:)). At that point, the scheme puts the 
                ! unallocated b_inc into the labile pool. What we 
                ! want is that the unallocated fraction is removed from 
                ! ::bm_alloc_tot such that only the allocated C is removed 
                ! from the labile pool. b_inc_tot will be moved back into
                ! the labile pool in 5.2.11
                bm_alloc_tot(ipts,j) = bm_alloc_tot(ipts,j) - b_inc_tot
                tmp_bm(ipts,j,ilabile,icarbon) = &
                     tmp_bm(ipts,j,ilabile,icarbon) + b_inc_tot 

                ! Wrap-up ordinary growth  
                ! Calculate C that was not allocated (b_inc_tot), the 
                ! equation should read b_inc_tot = b_inc_tot - b_inc_tot 
                ! note that Cf_inc was already substracted
                b_inc_tot = zero 

                ! Debug
                IF (printlev_loc.GE.4 .AND. ipts.EQ.test_grid .AND. j.EQ.test_pft) THEN
                   WRITE(numout,*) 'No wood growth, move remaining C to labile pool'
                   WRITE(numout,*) 'bm_alloc_tot_new, ',bm_alloc_tot(ipts,j)
                   WRITE(numout,*) 'wrap-up ordinary allocation, left b_inc_tot, ', &
                        b_inc_tot 

                ENDIF
                !-

             !! 5.2.8 Error - the allocation scheme is overspending 
             ELSEIF (b_inc_tot .LT. min_stomate) THEN

                IF (b_inc_tot .LT. -10*EPSILON(zero)) THEN

                   ! Something is wrong with the calculations
                   WRITE(numout,*) 'WARNING 7: numerical problem overspending '//&
                        'in ordinary allocation'
                   WRITE(numout,*) 'WARNING 7: PFT, ipts: ',j,ipts
                   WRITE(numout,*) 'WARNING 7: b_inc_tot', b_inc_tot
                   IF(err_act.GT.1)THEN
                      CALL ipslerr_p (3,'growth_fun_all',&
                           'WARNING 7: numerical problem',&
                           'overspending in ordinary allocation','')
                   ENDIF

                ELSE

                   IF (printlev_loc.GE.4 .AND. ipts.EQ.test_grid .AND. j.EQ.test_pft) THEN

                      ! Succesful allocation
                      WRITE(numout,*) 'Insufficient carbon for ordinary &
                                        allocation'

                   ENDIF

                ENDIF

                ! Although the biomass components respect the allometric 
                ! relationships, there is no carbon left to allocate                      
                b_inc_tot = zero

             ENDIF !End ordinary allocation

             ! Update circ_class_biomass
             circ_class_biomass(ipts,j,:,ilabile,icarbon) = &
                  biomass_to_cc(tmp_bm(ipts,j,ilabile,icarbon),&
                  circ_class_biomass(ipts,j,:,ilabile,icarbon),&
                  circ_class_n(ipts,j,:))
             circ_class_biomass(ipts,j,:,icarbres,icarbon) = &
                  biomass_to_cc(tmp_bm(ipts,j,icarbres,icarbon),&
                  circ_class_biomass(ipts,j,:,icarbres,icarbon),&
                  circ_class_n(ipts,j,:))
          
             !! 5.2.9 Error checking
             IF ( b_inc_tot .GT. min_stomate) THEN

                ! Although this should not happen. In case the functional 
                ! allocation did not consume all the allocatable carbon, 
                ! the remaining C is left for the next day. The numerical 
                ! precision of the allocation scheme (i.e. the linearisation) 
                ! is similar to min_stomate (i.e. 10-8) resulting in 'false' 
                ! warnings.
                WRITE(numout,*) 'WARNING 8: b_inc_tot greater than min_stomate '//&
                     'force allocation'
                WRITE(numout,*) 'WARNING 8: PFT, ipts: ',j,ipts
                WRITE(numout,*) 'WARNING 8: b_inc_tot, ', b_inc_tot
                IF(err_act.GT.1)THEN
                   CALL ipslerr_p (3,'growth_fun_all',&
                        'WARNING 8: b_inc_tot greater than min_stomate',&
                        'force allocation','')
                ENDIF

             ELSEIF ( (b_inc_tot .LT. min_stomate) .AND. (b_inc_tot .GE. zero) ) THEN

                ! Successful allocation
                IF (printlev_loc.GE.4 .AND. ipts.EQ.test_grid .AND. j.EQ.test_pft) THEN
                   WRITE(numout,*) 'Successful allocation'
                ENDIF

             ELSE

                ! Something possibly important was overlooked
                IF ( (b_inc_tot .LT. zero) .AND. &
                     (b_inc_tot .GE. -100*min_stomate) ) THEN
                   IF (printlev_loc.GE.4 .AND. ipts.EQ.test_grid .AND. j.EQ.test_pft) THEN
                      WRITE(numout,*) 'Marginally successful allocation - '//&
                           'precision better than 5 10-6'
                      WRITE(numout,*) 'PFT, b_inc_tot', j, b_inc_tot
                   ENDIF
                ELSE
                   WRITE(numout,*) 'WARNING 9: Logical flaw '//&
                        'unexpected result in the ordinary allocation'
                   WRITE(numout,*) 'WARNING 9: b_inc_tot, ',b_inc_tot
                   WRITE(numout,*) 'WARNING 9: PFT, ipts: ',j,ipts
                   IF(err_act.GT.1)THEN
                      CALL ipslerr_p (3,'growth_fun_all',&
                           'WARNING 9: Logical flaw',&
                           'unexpected result in the ordinary allocation','')
                   ENDIF
                ENDIF

             ENDIF

             !Debug
             IF (printlev_loc.GE.4 .AND. ipts.EQ.test_grid .AND. j.EQ.test_pft) THEN
                WRITE(numout,*) 'Final allocation', ipts, j
                WRITE(numout,*) 'Cl, Cs, Cr', Cl(:), Cs(:), Cr(:) 
                WRITE(numout,*) 'Cl_incp, Cs_incp, Cr_incp, ', &
                     Cl_incp(:), Cs_incp(:), Cr_incp(:)
                WRITE(numout,*) 'Cl_inc, Cs_ins, Cr_inc, Cf_inc, ', &
                     Cl_inc(:), Cs_inc(:), Cr_inc(:), Cf_inc(:)
                WRITE(numout,*) 'b_inc_tot, ', b_inc_tot
                WRITE(numout,*) 'Old ba, delta_ba, new ba, ', circ_class_ba_eff(:), &
                     delta_ba(:), circ_class_ba_eff(:)+delta_ba(:)
                DO l=1,ncirc
                   WRITE(numout,*) 'Circ_class_biomass, ', &
                        circ_class_biomass(ipts,j,l,:,icarbon)
                ENDDO
             ENDIF
             !-

             !! 5.2.10 Wrap-up phenological and ordinary allocation
             Cl_inc(:) = Cl_inc(:) + Cl_incp(:)
             Cr_inc(:) = Cr_inc(:) + Cr_incp(:)
             Cs_inc(:) = Cs_inc(:) + Cs_incp(:)
             residual(ipts,j) = b_inc_tot

             !+++CHECK+++
             ! All options to leave allocation have a b_inc_tot of zero.
             ! This code may no longer be needed. As we are working towards
             ! a dealine it was left in. It should not be harmful, it may 
             ! just complexify the code and slowdown the model a bit. 
  
             !! 5.2.11 Account for the residual
             !  The residual is usually around ::min_stomate but we deal 
             !  with it anyway to make sure the mass balance is closed
             !  and as a way to detect errors. Move the unallocated carbon
             !  back into the labile pool
             IF (tmp_bm(ipts,j,ilabile,icarbon) + residual(ipts,j) .LE. min_stomate) THEN

                deficit = tmp_bm(ipts,j,ilabile,icarbon) + residual(ipts,j)

                ! The deficit is less than the carbon reserve
                IF (-deficit .LE. tmp_bm(ipts,j,icarbres,icarbon)) THEN

                   ! Pay the deficit from the reserve pool
                   tmp_bm(ipts,j,icarbres,icarbon) = &
                        tmp_bm(ipts,j,icarbres,icarbon) + deficit
                   tmp_bm(ipts,j,ilabile,icarbon)  = &
                        tmp_bm(ipts,j,ilabile,icarbon) - deficit

                ELSE

                   ! Not enough carbon to pay the deficit
                   ! There is likely a bigger problem somewhere in
                   ! this routine
                   WRITE(numout,*) 'WARNING 11: PFT, ipts: ',j,ipts
                   WRITE(numout,*) 'resiudal, labile, deficit, ', &
                        residual(ipts,j), tmp_bm(ipts,j,ilabile,icarbon), &
                        deficit, tmp_bm(ipts,j,icarbres,icarbon) 
                   CALL ipslerr_p (3,'growth_fun_all',&
                        'WARNING 11: numerical problem overspending ',&
                        'when trying to account for unallocatable C ','')

                ENDIF

             ELSE

                ! Move the unallocated carbon back into the labile pool
                tmp_bm(ipts,j,ilabile,icarbon) = &
                     tmp_bm(ipts,j,ilabile,icarbon) + residual(ipts,j)

             ENDIF
             !+++++++++++
            
             !! 5.2.12 Distribute stand level ilabile and icarbres at the tree level
             !  The labile and carbres pools are calculated at the stand level but
             !  have to be redistributed at the tree level. Tree level biomass is the
             !  prognostic variable in ORCHIDEE. Biomass is sometimes used as a local
             !  variable mainly to deal with the reserve and labile pools.
             circ_class_biomass(ipts,j,:,ilabile,icarbon) = &
                  biomass_to_cc(tmp_bm(ipts,j,ilabile,icarbon),&
                  circ_class_biomass(ipts,j,:,ilabile,icarbon),&
                  circ_class_n(ipts,j,:))
             circ_class_biomass(ipts,j,:,icarbres,icarbon) = &
                  biomass_to_cc(tmp_bm(ipts,j,icarbres,icarbon),&
                  circ_class_biomass(ipts,j,:,icarbres,icarbon),&
                  circ_class_n(ipts,j,:))
            
             !! 5.2.13 Standardise allocation factors
             !  Strictly speaking the allocation factors do not need to be 
             !  calculated because the functional allocation scheme allocates 
             !  absolute amounts of carbon. Hence, Cl_inc could simply be 
             !  added to tmp_bm(:,:,ileaf,icarbon), Cr_inc to 
             !  tmp_bm(:,:,iroot,icarbon), etc. However, using allocation 
             !  factors bears some elegance in respect to distributing the 
             !  growth respiration if this would be required. Further it 
             !  facilitates comparison to the resource limited allocation 
             !  scheme (stomate_growth_res_lim.f90) and it comes in handy 
             !  for model-data comparison. This allocation takes place at 
             !  the tree level - note that ::biomass is the only prognostic 
             !  variable from the tree-based allocation
             !  WARNING: the reserves pools are ignored when calculating
             !  the allocation factors. Make sure this is OK for you before
             !  using these factors.

             !  Allocation   
             Cl_inc(:) = MAX(zero, circ_class_n(ipts,j,:) * Cl_inc(:))
             Cr_inc(:) = MAX(zero, circ_class_n(ipts,j,:) * Cr_inc(:))
             Cs_inc(:) = MAX(zero, circ_class_n(ipts,j,:) * Cs_inc(:))
             Cf_inc(:) = MAX(zero, circ_class_n(ipts,j,:) * Cf_inc(:))

             ! Total_inc is based on the updated Cl_inc, Cr_inc, Cs_inc and 
             ! Cf_inc. Therefore, do not multiply
             ! circ_class_n(ipts,j,:) again
             total_inc = SUM(Cf_inc(:) + Cl_inc(:) + Cs_inc(:) + Cr_inc(:))

             ! Relative allocation
             IF ( total_inc .GT. min_stomate ) THEN

                Cl_inc(:) = Cl_inc(:) / total_inc
                Cs_inc(:) = Cs_inc(:) / total_inc
                Cr_inc(:) = Cr_inc(:) / total_inc
                Cf_inc(:) = Cf_inc(:) / total_inc

             ELSE

                bm_alloc_tot(ipts,j) = zero
                Cl_inc(:) = zero
                Cs_inc(:) = zero
                Cr_inc(:) = zero
                Cf_inc(:) = zero

             ENDIF

             !! 5.2.13 Convert allocation to allocation facors
             !  Convert allocation of individuals to ORCHIDEE's allocation 
             !  factors - see comment for 5.2.5. Aboveground sapwood 
             !  allocation is age dependent in trees. ::alloc_min and 
             !  ::alloc_max must range between 0 and 1. 
             alloc_sap_above = alloc_min(j) + ( alloc_max(j) - alloc_min(j) ) * &
                  ( 1. - EXP( -age(ipts,j) / demi_alloc(j) ) )

             ! Leaf, wood, root and fruit allocation. Note that the X_inc
             ! are normalized before being used here.
             f_alloc(ipts,j,ileaf) = SUM(Cl_inc(:))
             f_alloc(ipts,j,isapabove) = SUM(Cs_inc(:)*alloc_sap_above)
             f_alloc(ipts,j,isapbelow) = SUM(Cs_inc(:)*(1.-alloc_sap_above))
             f_alloc(ipts,j,iroot) = SUM(Cr_inc(:))
             f_alloc(ipts,j,ifruit) = SUM(Cf_inc(:))

             ! Store f_alloc per circ_class to calculate the allocation
             ! in circ_class_biomass after bm_alloc_tot has been checked
             ! for the N availability. Note that the X_inc
             ! are normalized before being used here.
             f_alloc_circ(ipts,:,ileaf) = Cl_inc(:)
             f_alloc_circ(ipts,:,isapabove) = Cs_inc(:)*alloc_sap_above
             f_alloc_circ(ipts,:,isapbelow) = Cs_inc(:)*(1.-alloc_sap_above)
             f_alloc_circ(ipts,:,iroot) = Cr_inc(:)
             f_alloc_circ(ipts,:,ifruit) = Cf_inc(:)

             ! Debug
             IF (printlev_loc.GE.4 .AND. ipts.EQ.test_grid .AND. j.EQ.test_pft) THEN
                tempi = zero
                DO icir = 1,ncirc
                   IF (Cl_inc(icir) .LT. zero .OR. &
                        Cs_inc(icir) * alloc_sap_above .LT. zero .OR. &
                        Cs_inc(icir) * (un - alloc_sap_above) .LT. zero .OR. &
                        Cr_inc(icir) .LT. zero .OR. & 
                        Cf_inc(icir) .LT. zero .OR. &
                        total_inc .LT. zero .OR. &
                        circ_class_n(ipts,j,icir) .LT. zero) THEN
                      WRITE(numout,*) 'Cl_inc, ', j, Cl_inc(icir)
                      WRITE(numout,*) 'Cs_inc aboveground, ', j, &
                           Cs_inc(icir) * alloc_sap_above
                      WRITE(numout,*) 'Cs_inc aboveground, ', j, &
                           Cs_inc(icir) * (un-alloc_sap_above)
                      WRITE(numout,*) 'Cr_inc, ', j, Cr_inc(icir)
                      WRITE(numout,*) 'Cf_inc, ', j, Cf_inc(icir)
                      WRITE(numout,*) 'total_inc, ', j, total_inc
                      WRITE(numout,*) 'circ_class_n, ', j, circ_class_n(ipts,j,icir)
                      CALL ipslerr_p (3,'growth_fun_all',&
                           'WARNING 11bis: the solution has negative values',&
                           'None of these variables should be negative','')
                   ENDIF
                ENDDO
             ENDIF
             !-

          ELSEIF (is_tree(j)) THEN
       
             ! bm_alloc_tot was less than min_stomate. No effort
             ! to allocate but this little bit of carbon should be
             ! correctly accounted for.
             residual(ipts,j) = bm_alloc_tot(ipts,j)
             
             ! Debug
             IF(printlev_loc.GE.4 .AND. ipts.EQ.test_grid .AND. j.EQ.test_pft) WRITE(numout,*) &
                  'there is no tree biomass to allocate, PFT, ', j
             !-

          ENDIF ! Is there biomass to allocate (§5.2 - far far up)

          !! 5.3 Calculate allocated biomass pools for grasses and crops
          !  Only possible if there is biomass to allocate
          IF ( .NOT. is_tree(j) .AND. bm_alloc_tot(ipts,j) .GT. min_stomate ) THEN

             !! 5.3.1 Scaling factor to convert variables to the individual plant
             !  Allocation is on an individual basis (gC ind-1). Stand-level variables 
             !  need to convert to a single individual. The absence of sapwood makes 
             !  this irrelevant because the allocation reduces to a linear function 
             !  (contrary to the non-linearity of tree allocation). For the
             !  beauty of consistency, the transformations will be implemented.
             !  Different approach between the DGVM and statitic approach
             IF (ok_dgvm) THEN

                ! The DGVM does NOT work with the functional allocation. Consider
                ! this code as a placeholder. The original code had two different 
                ! transformations to calculate the scalars. Both could be used but 
                ! the units will differ. For consistency only one was retained 
                ! scal = ind(ipts,j) * cn_ind(ipts,j) / veget_max(ipts,j)
                scal(ipts,j) = veget_max(ipts,j) / circ_class_n(ipts,j,1)

             ELSE

                ! By dividing the actual biomass by the number of individuals
                ! the biomass of an individual is obtained. Note that a grass/crop 
                ! individual was defined as 1m-2 of vegetation 
                scal(ipts,j) = 1.

             ENDIF

             !! 5.3.2 Current biomass pools per grass/crop (gC ind^-1)
             !  Cs has too many dimensions for grass/crops. To have a consistent 
             !  notation the same variables are used as for trees but the dimension 
             !  of Cs, Cl and Cr i.e. ::ncirc should be ignored            
             Cs(:) = circ_class_biomass(ipts,j,1,isapabove,icarbon) * scal(ipts,j)
             Cr(:) = circ_class_biomass(ipts,j,1,iroot,icarbon) * scal(ipts,j)
             Cl(:) = circ_class_biomass(ipts,j,1,ileaf,icarbon) * scal(ipts,j)
             Ch(:) = zero

             ! Quantify and account for nitrogen limitation on allometric
             ! allocation. The same code as in the section 5.4 is used exept
             ! that we don't use allocation coeff to modulate n_avail. So
             ! costf=1. It is a strong assamption compared to the previous version. 
             ! It means that ordinary allocation can only happens when allometric
             ! allocation is is ok. In other case no wood growth is allowed.
             ! In the case of a strong limitation by Nitrogen, the growth period
             ! for sapwood will be shorten because we reach allometry late in
             ! the growing season.
             n_avail = MAX(tmp_bm(ipts,j,ilabile,initrogen)*0.9,0.0)
             bm_supply_n = n_avail  / (1.-frac_growthresp_dyn) * &
                           cn_leaf(ipts,j)

             ! Calculate how much carbon could be allocated with the available nitrogen
             bm_supply_n = n_avail  / (1.-frac_growthresp_dyn) * &
                  cn_leaf(ipts,j)

             ! If there is not enough nitrogen, move nitrogen from the reserve
             ! as much as needed, keeping 10% of reserve (arbitral portion)
             IF(bm_alloc_tot(ipts,j) .GT. bm_supply_n &
                  .AND. n_avail .GT. zero) THEN

                ! Calculate the deficit
                n_deficit = bm_alloc_tot(ipts,j) * (1.-frac_growthresp_dyn) / &
                     cn_leaf(ipts,j) - n_avail

                IF(n_deficit .LE. tmp_bm(ipts,j,icarbres,initrogen) * 0.9) THEN

                   ! Enougn N in the reserve pools to fill the labile pool
                   n_avail = n_avail + n_deficit
                   bm_supply_n = n_avail / (1.-frac_growthresp_dyn) * &
                        cn_leaf(ipts,j)
                   tmp_bm(ipts,j,icarbres,initrogen) = tmp_bm(ipts,j,icarbres,initrogen) - &
                        (n_avail/0.9 - tmp_bm(ipts,j,ilabile,initrogen))
                   tmp_bm(ipts,j,ilabile,initrogen) = n_avail/0.9

                   ! tmp_bm is a temporary varaiable so the prognostic variable, i.e.,
                   ! circ_class_biomass also needs to be updated.
                   circ_class_biomass(ipts,j,:,icarbres,initrogen) = &
                        biomass_to_cc(tmp_bm(ipts,j,icarbres,initrogen),&
                        circ_class_biomass(ipts,j,:,icarbres,initrogen),&
                        circ_class_n(ipts,j,:))
                   circ_class_biomass(ipts,j,:,ilabile,initrogen) = &
                        biomass_to_cc(tmp_bm(ipts,j,ilabile,initrogen),&
                        circ_class_biomass(ipts,j,:,ilabile,initrogen),&
                        circ_class_n(ipts,j,:))

                ELSE

                   ! Deficit exceeds 90% of reserve. fill labile as much as
                   ! possible
                   tmp_bm(ipts,j,ilabile,initrogen) = tmp_bm(ipts,j,ilabile,initrogen) + &
                        tmp_bm(ipts,j,icarbres,initrogen) * 0.9
                   tmp_bm(ipts,j,icarbres,initrogen) = tmp_bm(ipts,j,icarbres,initrogen) - &
                        tmp_bm(ipts,j,icarbres,initrogen) * 0.9

                   ! tmp_bm is a temporary varaiable so the prognostic variable, i.e.,
                   ! circ_class_biomass also needs to be updated.
                   circ_class_biomass(ipts,j,:,icarbres,initrogen) = &
                        biomass_to_cc(tmp_bm(ipts,j,icarbres,initrogen),&
                        circ_class_biomass(ipts,j,:,icarbres,initrogen),&
                        circ_class_n(ipts,j,:))
                   circ_class_biomass(ipts,j,:,ilabile,initrogen) = &
                        biomass_to_cc(tmp_bm(ipts,j,ilabile,initrogen),&
                        circ_class_biomass(ipts,j,:,ilabile,initrogen),&
                        circ_class_n(ipts,j,:))

                   ! Update the available nitrogen and the carbon that could be allocated
                   ! with that amount of nitrogen
                   n_avail = MAX(tmp_bm(ipts,j,ilabile,initrogen)*0.9,0.0)
                   bm_supply_n = n_avail / (1.-frac_growthresp_dyn) * &
                        cn_leaf(ipts,j)
                ENDIF

             ENDIF

             deltacnmax = 1. - exp(-((1.6 * MIN((1./cn_leaf(ipts,j))-&
                  (1./cn_leaf_min_2D(ipts,j)),0.) / &
                  ( (1./(cn_leaf_max_2D(ipts,j))) - &
                  (1./cn_leaf_min_2D(ipts,j)) ) )**4.1))
             
             IF ( bm_alloc_tot(ipts,j) .GT. bm_supply_n ) THEN

                IF (impose_cn) THEN

                   ! Calculate how much nitrogen is missing to allocate all the
                   ! carbon contained in bm_alloc_tot
                   n_deficit = (bm_alloc_tot(ipts,j)-bm_supply_n) * &
                        (1.-frac_growthresp_dyn) / cn_leaf(ipts,j)/0.9

                   ! The nitrogen missing to allocate the entire bm_alloc_tot will be taken
                   ! from the atmosphere and put in the labile pool. 
                   atm_to_bm(ipts,j,initrogen) = atm_to_bm(ipts,j,initrogen) + &
                        n_deficit/dt
                   tmp_bm(ipts,j,ilabile,initrogen) = &
                        tmp_bm(ipts,j,ilabile,initrogen) + n_deficit

                   ! tmp_bm is a temporary varaiable so the prognostic variable, i.e.,
                   ! circ_class_biomass also needs to be updated.
                   circ_class_biomass(ipts,j,:,ilabile,initrogen) = &
                        biomass_to_cc(tmp_bm(ipts,j,ilabile,initrogen),&
                        circ_class_biomass(ipts,j,:,ilabile,initrogen),&
                        circ_class_n(ipts,j,:))

                   ! Estimate the nitrogen pool that is required to allocate all the
                   ! carbon in bm_alloc_tot.
                   n_alloc_tot(ipts,j) = bm_alloc_tot(ipts,j) * &
                        (1.-frac_growthresp_dyn)/cn_leaf(ipts,j)
                  
                ELSE

                   IF (printlev_loc .GE. 4 .AND. ipts.EQ.test_grid .AND. j.EQ.test_pft) THEN
                      WRITE(numout,*) 'N-limitation before allocation'
                   ENDIF 
                   deltacnmax = Dmax * (1.-deltacnmax) 
                   deltacn = n_avail /  ( bm_alloc_tot(ipts,j) * &
                        (1.-frac_growthresp_dyn) * 1./cn_leaf(ipts,j) )
                   deltacn = MIN(MAX(deltacn,1.0-deltacnmax),1.0)
                    
                   n_alloc_tot(ipts,j) =  MIN( n_avail , & 
                        bm_alloc_tot(ipts,j) * (1.-frac_growthresp_dyn) * &
                        MAX(MIN( 1./cn_leaf(ipts,j)*deltacn, 1./cn_leaf_min_2D(ipts,j)), &
                        1./cn_leaf_max_2D(ipts,j)) ) 
                   
                   tmp_bm(ipts,j,ilabile,icarbon) = &
                        tmp_bm(ipts,j,ilabile,icarbon) + &
                        bm_alloc_tot(ipts,j)
                   
                   bm_alloc_tot(ipts,j) = MIN( bm_alloc_tot(ipts,j) , &
                        n_alloc_tot(ipts,j) / (1.-frac_growthresp_dyn) / &
                        MAX(MIN(1./cn_leaf(ipts,j)*deltacn, &
                        1./cn_leaf_min_2D(ipts,j)), 1./cn_leaf_max_2D(ipts,j)) )
                   
                   tmp_bm(ipts,j,ilabile,icarbon) = tmp_bm(ipts,j,ilabile,icarbon) - &
                        bm_alloc_tot(ipts,j) 

                ENDIF ! if impose_cn

             ELSE

                deltacnmax=Dmax * deltacnmax
                deltacn = n_avail / ( bm_alloc_tot(ipts,j) * &
                     (1.-frac_growthresp_dyn) * 1./cn_leaf(ipts,j) ) 
                deltacn=MIN(MAX(deltacn,1.0),1.+deltacnmax)

                n_alloc_tot(ipts,j) =  MIN( n_avail , & 
                     bm_alloc_tot(ipts,j) * (1.-frac_growthresp_dyn) * &
                     MAX(MIN(1./cn_leaf(ipts,j)*deltacn, & 
                     1./cn_leaf_min_2D(ipts,j)),1./cn_leaf_max_2D(ipts,j)) )
              
             ENDIF

             ! Total amount of carbon that needs to ba allocated (::bm_alloc_tot). 
             ! bm_alloc_tot is in gC m-2 day-1. At 1 m2 there are ::ind number of 
             ! trees. We calculate the allocation for ::ncirc trees. Hence b_inc_tot 
             ! needs to be scaled in the allocation routines. For all cases were 
             ! allocation takes place for a single circumference class, scaling 
             ! could be done before the allocation. In the ordinary allocation 
             ! allocation takes place to all circumference classes at the same time. 
             ! Hence scaling takes place in that step for consistency we scale during 
             ! allocation. Note that b_inc (the carbon allocated to an individual 
             ! circumference class cannot be estimates at this point.
             IF (bm_alloc_tot(ipts,j).GT.min_stomate) THEN

                ! There is enough carbon to allocate
                b_inc_tot = bm_alloc_tot(ipts,j)

             ELSE

                ! There is so little carbon that it is not worth the hassle
                ! to allocate. Allocating very small amounts increases the
                ! risk to run into precision errors.
                tmp_bm(ipts,j,ilabile,icarbon) = tmp_bm(ipts,j,ilabile,icarbon) + &
                     bm_alloc_tot(ipts,j)
                b_inc_tot = zero

             ENDIF

             ! Labile carbon is updated in consequence
             circ_class_biomass(ipts,j,:,ilabile,icarbon) = &
                  biomass_to_cc(tmp_bm(ipts,j,ilabile,icarbon),&
                  circ_class_biomass(ipts,j,:,ilabile,icarbon),&
                  circ_class_n(ipts,j,:))

          END IF

          ! Intermediate mass balance check
          IF (err_act.EQ.4 .AND. .NOT.is_tree(j)) THEN

             ! Reset entire array to zero to calculate mass balance for each
             ! pixel x pft separatly
             pool_end(:,:,:) = zero
                
             ! Add bm_alloc_tot into the pool
             pool_end(ipts,j,icarbon) = pool_end(ipts,j,icarbon) + &
                     b_inc_tot * veget_max(ipts,j)
                            
             ! Check mass balance closure. Between intermediate check 1 and 3a
             ! bm_inc_tot was recalculated by accounting for the available nitrogen
             CALL intermediate_mass_balance_check(pool_start, pool_end, circ_class_biomass, &
                  circ_class_n, veget_max, bm_alloc_tot, gpp_daily, atm_to_bm, dt, npts, &
                  resp_maint, resp_growth, check_intern_init, ipts, j, '3a', 'ipft')
             
          ENDIF ! err_act.EQ.4

          ! The initial estimate of bm_alloc_tot was high enough to consider allocation
          ! but after accounting for the available nitrogen bm_alloc_tot may have
          ! dropped below the min_stomate threshold so it needs to be tested again
          IF ( .NOT. is_tree(j) .AND. bm_alloc_tot(ipts,j) .GT. min_stomate ) THEN

             !! 5.3.3 C-allocation for crops and grasses
             !  The mass conservation equations are detailed in the header of 
             !  this subroutine. The scheme assumes a functional relationships 
             !  between leaves and roots for grasses and crops. When carbon is 
             !  added to the leaf biomass pool, an increase in the root biomass 
             !  is to be expected to sustain water transport from the roots to 
             !  the leaves.

             !! 5.3.3.1 Do the biomass pools respect the pipe model?
             !  Do the current leaf, sapwood and root components respect the 
             !  allometric constraints? Calculate the optimal root and leaf mass, 
             !  given the current wood mass by using the basic allometric 
             !  relationships. Calculate the optimal sapwood mass as a function 
             !  of the current leaf and root mass.
             Cl_target(1) = MAX( Cs(1) * KF(ipts,j) , Cr(1) * LF(ipts,j), Cl(1) )
             Cs_target(1) = MAX( Cl_target(1) / KF(ipts,j), &
                  Cr(1) * LF(ipts,j) / KF(ipts,j), Cs(1) ) 
             Cr_target(1) = MAX( Cl_target(1) / LF(ipts,j), &
                  Cs_target(1) * KF(ipts,j) / LF(ipts,j), Cr(1) )
             
             ! Debug
             IF (printlev_loc.GE.4 .AND. ipts.EQ.test_grid .AND. j.EQ.test_pft) THEN
                WRITE(numout,*) 'bm_alloc_tot, ',bm_alloc_tot(ipts,j)
                WRITE(numout,*) 'Does the grass/crop needs reshaping?'
                WRITE(numout,*) 'KF, LF, ', KF(ipts,j), LF(ipts,j)
                WRITE(numout,*) 'Cl_target-Cl, ', Cl_target(1)-Cl(1), Cl_target(1), Cl(1)
                WRITE(numout,*) 'Cs_target-Cs, ', Cs_target(1)-Cs(1), Cs_target(1), Cs(1)
                WRITE(numout,*) 'Cr_target-Cr, ', Cr_target(1)-Cr(1), Cr_target(1), Cr(1)
             ENDIF
             !-

             !! 5.3.3.2 Phenological growth
             !  Phenological growth and reshaping of the grass/crop in line with 
             !  the pipe model. Turnover removes C from the different plant components 
             !  but at a component-specific rate, as such the allometric constraints 
             !  are distorted at every time step and should be restored before ordinary 
             !  growth can take place

             !! 5.3.3.2.1 The available C can sustain the present leaves and roots
             !  Calculate whether the structural c is in allometric balance. The target 
             !  values should always be larger than the current pools so the use of ABS 
             !  is redundant but was used to be on the safe side (here and in the rest 
             !  of the module) as it could help to find logical flaws.        
             IF ( ABS(Cs_target(1) - Cs(1)) .LT. min_stomate ) THEN

                Cs_incp(1) = MAX(zero, Cs_target(1) - Cs(1))

                ! Enough leaves and structural biomass, only grow roots
                IF ( ABS(Cl_target(1) - Cl(1))  .LT. min_stomate ) THEN

                   ! Allocate at the tree level to restore allometric balance
                   Cl_incp(1) = MAX(zero, Cl_target(1) - Cl(1))
                   Cr_incp(1) = MAX( MIN(b_inc_tot / circ_class_n(ipts,j,1) - &
                        Cs_incp(1) - Cl_incp(1), Cr_target(1) - Cr(1)), zero )

                   ! Write debug comments to output file
                   IF ((j.EQ.test_pft .AND. ipts.EQ.test_grid .AND. printlev_loc.GE.4) .OR. printlev_loc>=5) THEN
                      CALL comment(npts, Cl_target, Cl, Cs_target, Cs, Cr_target, Cr, &
                           delta_ba, ipts, j, l, b_inc_tot, Cl_incp, Cs_incp, Cr_incp, &
                           KF, LF, Cl_inc, Cs_inc, Cr_inc, Cf_inc, grow_wood, &
                           circ_class_n, 12)
                   ENDIF

                ! Sufficient structural C and roots, allocate C to leaves
                ELSEIF ( ABS(Cr_target(1) - Cr(1)) .LT. min_stomate ) THEN

                   ! Allocate at the tree level to restore allometric balance 
                   Cr_incp(1) = MAX(zero, Cr_target(1) - Cr(1))
                   Cl_incp(1) = MAX( MIN(b_inc_tot / circ_class_n(ipts,j,1) - &
                        Cs_incp(1) - Cr_incp(1), Cl_target(1) - Cl(1)), zero )

                   ! Update vegetation height
                   qm_height = biomass_to_lai(Cl(1) + Cl_incp(1),j) * lai_to_height(j)
               
                   ! Write debug comments to output file
                   IF ((j.EQ.test_pft .AND. ipts.EQ.test_grid .AND. printlev_loc.GE.4) .OR. printlev_loc>=5) THEN
                      CALL comment(npts, Cl_target, Cl, Cs_target, Cs, Cr_target, Cr, &
                           delta_ba, ipts, j, l, b_inc_tot, Cl_incp, Cs_incp, Cr_incp, &
                           KF, LF, Cl_inc, Cs_inc, Cr_inc, Cf_inc, grow_wood, &
                           circ_class_n, 13)
                   ENDIF

                ! Both leaves and roots are needed to restore the allometric relationships
                ELSEIF ( ABS(Cl_target(1) - Cl(1)) .GT. min_stomate .AND. &
                     ABS(Cr_target(1) - Cr(1)) .GT. min_stomate ) THEN                 

                   ! Allocate at the tree level to restore allometric balance
                   !  The equations can be rearanged and written as
                   !  (i) b_inc = Cl_inc + Cr_inc
                   !  (ii) Cr_inc = (Cl_inc+Cl)/LF - Cr
                   !  Substitue (ii) in (i) and solve for Cl_inc
                   !  <=> Cl_inc = (LF*(b_inc+Cr)-Cl)/(1+LF)
                   Cl_incp(1) = MIN( ((LF(ipts,j) * ((b_inc_tot/circ_class_n(ipts,j,1)) - &
                        Cs_incp(1) + Cr(1))) - Cl(1)) / (1 + LF(ipts,j)), &
                        Cl_target(1) - Cl(1) )
                   Cr_incp(1) = MIN ( ((Cl_incp(1) + Cl(1)) / LF(ipts,j)) - Cr(1), &
                        Cr_target(1) - Cr(1))

                   ! The imbalance between Cr and Cl can be so big that (Cl+Cl_inc)/LF 
                   ! is still less then the available root carbon (observed!). This 
                   ! would result in a negative Cr_incp
                   IF ( Cr_incp(1) .LT. zero ) THEN

                      Cl_incp(1) = MIN( b_inc_tot/circ_class_n(ipts,j,1) - &
                           Cs_incp(1), Cl_target(1) - Cl(1) )
                      Cr_incp(1) = b_inc_tot/circ_class_n(ipts,j,1) - &
                           Cs_incp(1) - Cl_incp(1)

                   ELSEIF (Cl_incp(1) .LT. zero) THEN

                      Cr_incp(1) = MIN( b_inc_tot/circ_class_n(ipts,j,1) - &
                           Cs_incp(1), Cr_target(1) - Cr(1) )
                      Cl_incp(1) = (b_inc_tot/circ_class_n(ipts,j,1)) - &
                           Cs_incp(1) - Cr_incp(1)

                   ENDIF

                   ! Update vegetation height
                   qm_height = biomass_to_lai(Cl(1) + Cl_incp(1),j) * lai_to_height(j)

                   ! Write debug comments to output file
                   IF ((j.EQ.test_pft .AND. ipts.EQ.test_grid .AND. printlev_loc.GE.4) .OR. printlev_loc>=5) THEN
                      CALL comment(npts, Cl_target, Cl, Cs_target, Cs, Cr_target, Cr, &
                           delta_ba, ipts, j, l, b_inc_tot, Cl_incp, Cs_incp, Cr_incp, &
                           KF, LF, Cl_inc, Cs_inc, Cr_inc, Cf_inc, grow_wood, &
                           circ_class_n, 14)
                   ENDIF    

                ELSE

                   WRITE(numout,*) 'WARNING 12: Exc 1-3 unexpected exception'
                   WRITE(numout,*) 'WARNING 12: PFT, ipts: ',j,ipts
                   IF(err_act.GT.1)THEN
                      CALL ipslerr_p (3,'growth_fun_all',&
                          'WARNING 12: Exc 1-3 unexpected exception','','') 
                   ENDIF

                ENDIF

             !! 5.3.3.3.2 Enough leaves to sustain the structural C and roots
             ELSEIF ( ABS(Cl_target(1) - Cl(1)) .LT. min_stomate ) THEN
                
                Cl_incp(1) = MAX(zero, Cl_target(1) - Cl(1))

                ! Enough leaves and structural C, only grow roots
                ! This duplicates Exc 1 and these lines should never be called 
                IF ( ABS(Cs_target(1) - Cs(1)) .LT. min_stomate ) THEN

                   ! Allocate at the tree level to restore allometric balance
                   Cs_incp(1) = MAX(zero, Cs_target(1) - Cs(1))
                   Cr_incp(1) = MAX( MIN(b_inc_tot/circ_class_n(ipts,j,1) - &
                        Cl_incp(1) - Cs_incp(1), Cr_target(1) - Cr(1)), zero )

                   ! Write debug comments to output file 
                   IF ((j.EQ.test_pft .AND. ipts.EQ.test_grid .AND. printlev_loc.GE.4) .OR. printlev_loc>=5) THEN
                      CALL comment(npts, Cl_target, Cl, Cs_target, Cs, Cr_target, Cr, &
                           delta_ba, ipts, j, l, b_inc_tot, Cl_incp, Cs_incp, Cr_incp, &
                           KF, LF, Cl_inc, Cs_inc, Cr_inc, Cf_inc, grow_wood, &
                           circ_class_n, 15)
                   ENDIF 

                ! Enough leaves and roots. Need to grow structural C to support 
                ! the available canopy and roots
                ELSEIF ( ABS(Cr_target(1) - Cr(1)) .LT. min_stomate ) THEN

                   Cr_incp(1) = MAX(zero, Cr_target(1) - Cr(1))
                   Cs_incp(1) = MAX( MIN(b_inc_tot/circ_class_n(ipts,j,1) - &
                        Cr_incp(1) - Cl_incp(1), Cs_target(1) - Cs(1)), zero )

                   ! Write debug comments to output file 
                   IF ((j.EQ.test_pft .AND. ipts.EQ.test_grid .AND. printlev_loc.GE.4) .OR. printlev_loc>=5) THEN
                      CALL comment(npts, Cl_target, Cl, Cs_target, Cs, Cr_target, Cr, &
                           delta_ba, ipts, j, l, b_inc_tot, Cl_incp, Cs_incp, Cr_incp, &
                           KF, LF, Cl_inc, Cs_inc, Cr_inc, Cf_inc, grow_wood, &
                           circ_class_n, 16)
                   ENDIF

                ! Need both structural C and roots to restore the allometric relationships
                ELSEIF ( ABS(Cs_target(1) - Cs(1) ) .GT. min_stomate .AND. &
                     ABS(Cr_target(1) - Cr(1)) .GT. min_stomate ) THEN

                   !  First try if we can simply satisfy the allocation needs
                   IF (Cs_target(1) - Cs(1) + Cr_target(1) - Cr(1) .LE. &
                        b_inc_tot/circ_class_n(ipts,j,1) - Cl_incp(1)) THEN
                         
                      Cr_incp(1) = Cr_target(1) - Cr(1)
                      Cs_incp(1) = Cs_target(1) - Cs(1)

                   ! Try to satisfy the need for the roots
                   ELSEIF (Cr_target(1) - Cr(1) .LE. &
                        b_inc_tot/circ_class_n(ipts,j,1) - Cl_incp(1)) THEN

                      Cr_incp(1) = Cr_target(1) - Cr(1)
                      Cs_incp(1) = b_inc_tot/circ_class_n(ipts,j,1) - &
                           Cl_incp(1) - Cr_incp(1)
                      

                   ! There is not enough use whatever is available
                   ELSE
                         
                      Cr_incp(1) = b_inc_tot/circ_class_n(ipts,j,1) - Cl_incp(1)
                      Cs_incp(1) = zero
                         
                   ENDIF

                   ! Write debug comments to output file 
                   IF ((j.EQ.test_pft .AND. ipts.EQ.test_grid .AND. printlev_loc.GE.4) .OR. printlev_loc>=5) THEN
                      CALL comment(npts, Cl_target, Cl, Cs_target, Cs, Cr_target, Cr, &
                           delta_ba, ipts, j, l, b_inc_tot, Cl_incp, Cs_incp, Cr_incp, &
                           KF, LF, Cl_inc, Cs_inc, Cr_inc, Cf_inc, grow_wood, &
                           circ_class_n, 17)
                   ENDIF 

                ELSE

                   WRITE(numout,*) 'WARNING 13: Exc 4-6 unexpected exception'
                   WRITE(numout,*) 'WARNING 13: PFT, ipts: ',j,ipts
                   IF(err_act.GT.1)THEN
                      CALL ipslerr_p (3,'growth_fun_all',&
                           'WARNING 13: Exc 4-6 unexpected exception','','')
                   ENDIF

                ENDIF

             !! 5.3.3.3.3 Enough roots to sustain the wood and leaves
             ELSEIF ( ABS(Cr_target(1) - Cr(1)) .LT. min_stomate ) THEN

                Cr_incp(1) = MAX(zero, Cr_target(1) - Cr(1)) 

                ! Enough roots and wood, only grow leaves
                ! This duplicates Exc 2 and these lines should thus never be called 
                IF ( ABS(Cs_target(1) - Cs(1)) .LT. min_stomate ) THEN

                   ! Allocate at the tree level to restore allometric balance
                   Cs_incp(1) = MAX(zero, Cs_target(1) - Cs(1)) 
                   Cl_incp(1) = MAX( MIN(b_inc_tot/circ_class_n(ipts,j,1) - &
                        Cr_incp(1) - Cs_incp(1), Cl_target(1) - Cl(1)), zero )

                   ! Update vegetation height
                   qm_height = biomass_to_lai(Cl(1) + Cl_incp(1),j) * lai_to_height(j)

                   ! Write debug comments to output file 
                   IF ((j.EQ.test_pft .AND. ipts.EQ.test_grid .AND. printlev_loc.GE.4) .OR. printlev_loc>=5) THEN
                      CALL comment(npts, Cl_target, Cl, Cs_target, Cs, Cr_target, Cr, &
                           delta_ba, ipts, j, l, b_inc_tot, Cl_incp, Cs_incp, Cr_incp, &
                           KF, LF, Cl_inc, Cs_inc, Cr_inc, Cf_inc, grow_wood, &
                           circ_class_n, 18)
                   ENDIF 

                ! Enough leaves and roots. Need to grow sapwood to support the 
                ! available canopy and roots. Duplicates Exc. 4 and these lines 
                ! should thus never be called 
                ELSEIF ( ABS(Cl_target(1) - Cl(1)) .LT. min_stomate ) THEN

                   ! Allocate at the tree level to restore allometric balance
                   Cl_incp(1) = MAX(zero, Cl_target(1) - Cl(1)) 
                   Cs_incp(1) = MAX( MIN(b_inc_tot/circ_class_n(ipts,j,1) - &
                        Cr_incp(1) - Cl_incp(1), Cs_target(1) - Cs(1) ), zero )

                   ! Write debug comments to output file 
                   IF ((j.EQ.test_pft .AND. ipts.EQ.test_grid .AND. printlev_loc.GE.4) .OR. printlev_loc>=5) THEN
                      CALL comment(npts, Cl_target, Cl, Cs_target, Cs, Cr_target, Cr, &
                           delta_ba, ipts, j, l, b_inc_tot, Cl_incp, Cs_incp, Cr_incp, &
                           KF, LF, Cl_inc, Cs_inc, Cr_inc, Cf_inc, grow_wood, &
                           circ_class_n, 19)
                   ENDIF                       

                ! Need both wood and leaves to restore the allometric relationships
                ELSEIF ( ABS(Cs_target(1) - Cs(1)) .GT. min_stomate .AND. &
                     ABS(Cl_target(1) - Cl(1)) .GT. min_stomate ) THEN

                   ! circ_class_ba_eff and circ_class_height_eff are already calculated
                   ! for a tree in balance. It would be rather complicated to follow
                   ! the allometric rules for wood allocation (implying changes in height 
                   ! and basal area) because the tree is not in balance.First try if we 
                   ! can simply satisfy the allocation needs
                   IF (Cs_target(1) - Cs(1) + Cl_target(1) - Cl(1) .LE. &
                        b_inc_tot/circ_class_n(ipts,j,1) - Cr_incp(1)) THEN
                      
                      Cl_incp(1) = Cl_target(1) - Cl(1)
                      Cs_incp(1) = Cs_target(1) - Cs(1)
                      
                   ! Try to satisfy the need for leaves
                   ELSEIF (Cl_target(1) - Cl(1) .LE. &
                        b_inc_tot/circ_class_n(ipts,j,1) - Cr_incp(1)) THEN

                      Cl_incp(1) = Cl_target(1) - Cl(1)
                      Cs_incp(1) = b_inc_tot/circ_class_n(ipts,j,1) - &
                           Cr_incp(1) - Cl_incp(1)

                   ! There is not enough use whatever is available
                   ELSE

                      Cl_incp(1) = b_inc_tot/circ_class_n(ipts,j,1) - Cr_incp(1)
                      Cs_incp(1) = zero
                      
                   ENDIF
                      
                   ! Calculate the height of the expanded canopy
                   qm_height(ipts,j) = biomass_to_lai(Cl(1) + Cl_inc(1),j) * lai_to_height(j)

                   ! Write debug comments to output file 
                   IF ((j.EQ.test_pft .AND. ipts.EQ.test_grid .AND. printlev_loc.GE.4) .OR. printlev_loc>=5) THEN
                      CALL comment(npts, Cl_target, Cl, Cs_target, Cs, Cr_target, Cr, &
                           delta_ba, ipts, j, l, b_inc_tot, Cl_incp, Cs_incp, Cr_incp, &
                           KF, LF, Cl_inc, Cs_inc, Cr_inc, Cf_inc, grow_wood, &
                           circ_class_n, 20)
                   ENDIF

                ELSE

                   WRITE(numout,*) 'WARNING 14: Exc 7-9 unexpected exception'
                   WRITE(numout,*) 'WARNING 14: PFT, ipts: ',j, ipts
                   IF(err_act.GT.1)THEN
                      CALL ipslerr_p (3,'growth_fun_all',&
                           'WARNING 14: Exc 7-9 unexpected exception','','')
                   ENDIF

                ENDIF

             ! Either Cl_target, Cs_target or Cr_target should be zero
             ELSE

                ! Something possibly important was overlooked
                WRITE(numout,*) 'WARNING 15: Logical flaw in phenological allocation '
                WRITE(numout,*) 'WARNING 15: PFT, ipts: ',j, ipts
                WRITE(numout,*) 'Cs - Cs_target', Cs(1), Cs_target(1)
                WRITE(numout,*) 'Cl - Cl_target', Cl(1), Cl_target(1)
                WRITE(numout,*) 'Cr - Cr_target', Cr(1), Cr_target(1)
                IF(err_act.GT.1)THEN
                   CALL ipslerr_p (3,'growth_fun_all',&
                        'WARNING 15: Logical flaw in phenological allocation','','')
                ENDIF

             ENDIF

             !! 5.3.4 Wrap-up phenological allocation
             IF ( Cl_incp(1) .GE. zero .OR. Cr_incp(1) .GE. zero .OR. &
                  Cs_incp(1) .GE. zero) THEN

                ! Prevent overspending for leaves
                IF (b_inc_tot - circ_class_n(ipts,j,1) * Cl_incp(1) .LT. zero) THEN
                   Cl_incp(1) = b_inc_tot/circ_class_n(ipts,j,1)
                ENDIF
                b_inc_tot = MAX(zero,b_inc_tot - circ_class_n(ipts,j,1) * Cl_incp(1))
                
                ! Prevent overspending for roots
                IF (b_inc_tot - circ_class_n(ipts,j,1) * Cr_incp(1) .LT. zero) THEN  
                   Cr_incp(1) = b_inc_tot/circ_class_n(ipts,j,1)
                ENDIF
                b_inc_tot = MAX(zero,b_inc_tot - circ_class_n(ipts,j,1) * Cr_incp(1))
                
                ! Prevent overspending for sapwood
                IF (b_inc_tot - circ_class_n(ipts,j,1) * Cs_incp(1) .LT. zero) THEN  
                   Cs_incp(1) = b_inc_tot/circ_class_n(ipts,j,1)
                ENDIF
                b_inc_tot = MAX(zero,b_inc_tot - circ_class_n(ipts,j,1) * Cs_incp(1))
                   
                ! Fake allocation for less messy equations in next case, 
                ! incp needs to be added to inc at the end. 
                Cl(1) = Cl(1) + Cl_incp(1)
                Cr(1) = Cr(1) + Cr_incp(1)
                Cs(1) = Cs(1) + Cs_incp(1)
                
             ELSE

                ! The code was written such that the increment pools should be greater 
                ! than or equal to zero. If this is not the case, something fundamental 
                ! is wrong with the if-then constructs under §5.3.3.2
                WRITE(numout,*) 'WARNING 16: numerical problem, '//&
                     'one of the increment pools is less than zero'
                WRITE(numout,*) 'WARNING 16: Cl_incp(1), Cr_incp(1), Cs_incp(1), j, ipts',&
                     Cl_incp(1), Cr_incp(1), Cs_incp(1), j, ipts
                IF(err_act.GT.1)THEN
                   CALL ipslerr_p (3,'growth_fun_all',&
                        'WARNING 16: numerical problem',&
                        'one of the increment pools is less than zero','')
                ENDIF

             ENDIF

             ! Something is wrong with the calculations
             IF (b_inc_tot .LT. zero) THEN

                WRITE(numout,*) 'WARNING 17: numerical problem overspending '//&
                     'in the phenological allocation'
                WRITE(numout,*) 'WARNING 17: b_inc_tot, j, ipts',b_inc_tot, j, ipts 
                WRITE(numout,*) 'WARNING 17: Cl_incp, Cr_incp, Cs_incp, ', &
                     Cl_incp(1), Cr_incp(1), Cs_incp(1)
                IF(err_act.GT.1)THEN
                    CALL ipslerr_p (3,'growth_fun_all',&
                         'WARNING 17: numerical problem',&
                         'overspending in the phenological allocation','')
                ENDIF

             ENDIF

             ! Intermediate mass balance check. Note that this part of
             ! the code is in DO-loops over nvm and npts so the
             ! 'ipts' label is used in the mass balance check
             IF(err_act.EQ.4) THEN

                ! Reset pool_end
                pool_end(:,:,:) = zero
                
                ! Add bm_alloc_tot into the pool
                pool_end(ipts,j,icarbon) = pool_end(ipts,j,icarbon) + &
                     (SUM((Cl_incp(1) + Cs_incp(1) + Cr_incp(1)) * circ_class_n(ipts,j,:)) + &
                     b_inc_tot) * veget_max(ipts,j)
                
                ! Check mass balance closure. Between intermediate check 2b and 3b
                CALL intermediate_mass_balance_check(pool_start, pool_end, circ_class_biomass, &
                     circ_class_n, veget_max, bm_alloc_tot, gpp_daily, atm_to_bm, dt, npts, &
                     resp_maint, resp_growth, check_intern_init, ipts, j, '3b', 'ipft')
                
             END IF ! err_act.EQ.4

             ! Height depends on Cl, so update height when Cl gets updated
             qm_height(ipts,j) = biomass_to_lai(Cl(1),j) * lai_to_height(j) 
             grow_wood = .TRUE.
             
             ! Write debug comments to output file 
             IF ((j.EQ.test_pft .AND. ipts.EQ.test_grid .AND. printlev_loc.GE.4) .OR. printlev_loc>=5) THEN
                CALL comment(npts, Cl_target, Cl, Cs_target, Cs, Cr_target, Cr, &
                     delta_ba, ipts, j, l, b_inc_tot, Cl_incp, Cs_incp, Cr_incp, &
                     KF, LF, Cl_inc, Cs_inc, Cr_inc, Cf_inc, grow_wood, &
                     circ_class_n, 21)
             ENDIF

             !! 5.3.6 Ordinary growth
             !  Allometric relationship between components is respected, sustain 
             !  ordinary growth and allocate biomass to leaves, wood, roots and fruits.
             IF ( (ABS(Cl_target(1) - Cl(1) ) .LE. min_stomate) .AND. &
                  (ABS(Cs_target(1) - Cs(1) ) .LE. min_stomate) .AND. &
                  (ABS(Cr_target(1) - Cr(1) ) .LE. min_stomate) .AND. &
                  (grow_wood) .AND. &
                  (b_inc_tot .GT. min_stomate) ) THEN 

                ! Allocate fraction of carbon to fruit production (at the plant level)
                Cf_inc(1) = b_inc_tot * fruit_alloc(j)/circ_class_n(ipts,j,1)

                ! Residual carbon is allocated to the other components (b_inc_tot is 
                ! at the stand level)
                b_inc_tot = b_inc_tot * (1-fruit_alloc(j))

                ! Following allometric allocation
                ! (i) b_inc = Cl_inc + Cr_inc + Cs_inc
                ! (ii) Cr_inc = (Cl + Cl_inc)/LF - Cr
                ! (iii) Cs_inc = (Cl + Cl_inc) / KF - Cs
                ! Substitue (ii) and (iii) in (i) and solve for Cl_inc 
                ! <=> b_inc = Cl_inc + ( Cl_inc + Cl ) / KF - Cs + ( Cl_inc + Cl ) / LF - Cr
	        ! <=> b_inc = Cl_inc * ( 1.+ 1/KF + 1./LF ) + Cl/LF - Cs - Cr
	        ! <=> Cl_inc = ( b_inc - Cl/LF + Cs + Cr ) / ( 1.+ 1/KF + 1./LF )
                Cl_inc(1) = MAX( (b_inc_tot/circ_class_n(ipts,j,1) - Cl(1)/LF(ipts,j) - &
                     Cl(1)/KF(ipts,j) + Cs(1) + Cr(1)) / &
                     (1. + 1./KF(ipts,j) + 1./LF(ipts,j)), zero)
               
                IF (Cl_inc(1) .LE. zero) THEN

                   Cr_inc(:) = zero
                   Cs_inc(:) = zero

                ELSE

                   ! Wrap-up ordinary growth. Calculate C that was not allocated, note 
                   ! that Cf_inc was already substracted
                   ! Prevent overspending for leaves
                   IF (b_inc_tot - circ_class_n(ipts,j,1) * Cl_inc(1) .LT. zero) THEN
                      
                      Cl_inc(1) = b_inc_tot/circ_class_n(ipts,j,1)
                      b_inc_tot = MAX(zero,b_inc_tot - circ_class_n(ipts,j,1) * Cl_inc(1))

                      ! All carbon was used for leaves. No new growth for the roots and the stems
                      Cs_inc(1) = zero
                      Cr_inc(1) = zero
                                      
                   ELSE

                      ! If we end up here we can allocate the calculated Cl_inc
                      b_inc_tot = b_inc_tot - circ_class_n(ipts,j,1) * Cl_inc(1)

                      ! Calculate the height of the expanded canopy
                      qm_height(ipts,j) = biomass_to_lai(Cl(1) + Cl_inc(1),j) * &
                           lai_to_height(j)

                      ! Use the solution for Cl_inc to calculate Cr_inc and 
                      ! Cs_inc according to (ii) and (iii)
                      Cr_inc(1) = (Cl(1) + Cl_inc(1)) / LF(ipts,j) - Cr(1)

                      IF (b_inc_tot - circ_class_n(ipts,j,1) * Cr_inc(1) .LT. zero) THEN
                         
                         Cr_inc(1) = b_inc_tot/circ_class_n(ipts,j,1)
                         b_inc_tot = MAX(zero, b_inc_tot - circ_class_n(ipts,j,1) * Cr_inc(1))

                         ! No C left to grow new stems
                         Cs_inc(1) = zero

                      ELSE
                         
                         ! If we end up here we can allocate the calculated Cr_inc
                         b_inc_tot = b_inc_tot - circ_class_n(ipts,j,1) * Cr_inc(1)

                         ! Still carbon left so allocate it to the stems
                         ! Cs_inc(1) can be calculated as follows
                         ! Cs_inc(1) = (Cl(1) + Cl_inc(1)) / KF(ipts,j) - Cs(1)
                         ! It is easier and better for the mass balance closure to move
                         ! all the remaining b_inc_tot in the Cs_inc. Should be the
                         ! same except for the rounding and precision issues
                         Cs_inc(1) = MAX(zero,b_inc_tot/circ_class_n(ipts,j,1))
                         b_inc_tot = zero

                      END IF

                   END IF
                   
                END IF ! Cl_inc(1). LE. zero

                ! Write debug comments to output file 
                IF ((j.EQ.test_pft .AND. ipts.EQ.test_grid .AND. printlev_loc.GE.4) .OR. printlev_loc>=5) THEN
                   CALL comment(npts, Cl_target, Cl, Cs_target, Cs, Cr_target, Cr, &
                        delta_ba, ipts, j, l, b_inc_tot, Cl_incp, Cs_incp, Cr_incp, &
                        KF, LF, Cl_inc, Cs_inc, Cr_inc, Cf_inc, grow_wood, &
                        circ_class_n, 22)
                ENDIF

                ! Debug
                IF (printlev_loc.GE.4 .AND. ipts.EQ.test_grid .AND. j.EQ.test_pft) THEN
                   WRITE(numout,*) 'wrap-up ordinary allocation, left b_in_tot, ', &
                        b_inc_tot 
                ENDIF
                !-

                ! Intermediate mass balance check. Note that this part of
                ! the code is in DO-loops over nvm and one over npts so the
                ! 'ipts' label is used in the mass balance check
                IF (err_act.EQ.4) THEN

                   ! All carbon should have been allocated and the remainder was moved
                   ! back into the labile pool. b_inc_tot should be zero. If not, the 
                   ! calculation of pool_end is wrong.
                   IF(ABS(b_inc_tot).GT.min_stomate)THEN
                      WRITE(numout,*) 'b_inc_tot differs from zero, ', ipts,j,b_inc_tot
                      CALL ipslerr_p(3,'stomate_growth_fun_all.f90','intermediate mbcheck 3c',&
                           'b_inc_tot differs from zero','')
                   END IF

                   ! Reset pool_end
                   pool_end(:,:,:) = zero
                   WRITE(numout,*) 'Cl, ', Cl_incp(1)*circ_class_n(ipts,j,1)*veget_max(ipts,j), Cl_inc(1)*circ_class_n(ipts,j,1)*veget_max(ipts,j)
                   WRITE(numout,*) 'Cs, ', Cs_incp(1)*circ_class_n(ipts,j,1)*veget_max(ipts,j), Cs_inc(1)*circ_class_n(ipts,j,1)*veget_max(ipts,j)
                   WRITE(numout,*) 'Cr, ', Cr_incp(1)*circ_class_n(ipts,j,1)*veget_max(ipts,j), Cr_inc(1)*circ_class_n(ipts,j,1)*veget_max(ipts,j)
                   WRITE(numout,*) 'Cf, ', Cf_inc(1)*circ_class_n(ipts,j,1)
                   ! Add bm_alloc_tot into the pool
                   pool_end(ipts,j,icarbon) = ((Cl_incp(1) + Cr_incp(1) + Cs_incp(1) + &
                        Cl_inc(1) + Cs_inc(1) + Cr_inc(1) + Cf_inc(1)) * &
                        circ_class_n(ipts,j,1)) * veget_max(ipts,j)

                   ! Check mass balance closure. Between intermediate check 3b and 3c ordinary
                   ! allocation was accounted for. However, ordinary allocation was calculated in 
                   ! temporary variables but has not been accounted for yet. This check comes at 
                   ! the end of the allocation for grasses and crops.
                   CALL intermediate_mass_balance_check(pool_start, pool_end, circ_class_biomass, &
                        circ_class_n, veget_max, bm_alloc_tot, gpp_daily, atm_to_bm, dt, npts, &
                        resp_maint, resp_growth, check_intern_init, ipts, j, '3c', 'ipft')
                
                ENDIF ! err_act.EQ.4

             !! 5.3.7 Don't grow wood, use C to fill labile pool
             ELSEIF ( (.NOT. grow_wood) .AND. (b_inc_tot .GT. min_stomate) ) THEN

                ! grow_wood is always .TRUE. see 5.3.5 around line 3652. Is this
                ! intended or did we delete an if-statement?
                ! Calculate the C that needs to be distributed to the 
                ! labile pool. The fraction is proportional to the ratio 
                ! between the total allocatable biomass and the unallocated 
                ! biomass per tree (b_inc now contains the unallocated 
                ! biomass). At the end of the allocation scheme bm_alloc_tot 
                ! is substracted from the labile biomass pool to update the 
                ! biomass pool (tmp_bm(:,:,ilabile) = tmp_bm(:,:,ilabile) - 
                ! bm_alloc_tot(:,:)). At that point, the scheme puts the 
                ! unallocated b_inc into the labile pool. What we 
                ! want is that the unallocated fraction is removed from 
                ! ::bm_alloc_tot such that only the allocated C is removed 
                ! from the labile pool. b_inc_tot will be moved back into
                ! the labile pool in 5.2.11. resp_growth will be adjusted
                ! later in the code.
                bm_alloc_tot(ipts,j) = bm_alloc_tot(ipts,j) - b_inc_tot
                circ_class_biomass(ipts,j,1,ilabile,icarbon) = &
                     circ_class_biomass(ipts,j,1,ilabile,icarbon) + &
                     b_inc_tot

                ! Wrap-up ordinary growth  
                ! Calculate C that was not allocated (b_inc_tot), the 
                ! equation should read b_inc_tot = b_inc_tot - b_inc_tot 
                ! note that Cf_inc was already substracted
                b_inc_tot = zero
                

                ! Debug
                IF (printlev_loc.GE.3 .AND. ipts.EQ.test_grid .AND. j.EQ.test_pft) THEN
                   WRITE(numout,*) 'No wood growth, move remaining C to labile pool'
                   WRITE(numout,*) 'bm_alloc_tot, ',bm_alloc_tot(ipts,j)
                   WRITE(numout,*) 'wrap-up ordinary allocation, left b_inc_tot, ', &
                        b_inc_tot 
                ENDIF
                !- 

             !! 5.3.8 Error - the allocation scheme is overspending 
             ELSEIF (b_inc_tot .LE. min_stomate) THEN  
      
                IF (b_inc_tot .LT. zero) THEN

                   ! Something is wrong with the calculations
                   WRITE(numout,*) 'WARNING 18: numerical problem'//&
                        'overspending in ordinary allocation'
                   WRITE(numout,*) 'WARNING 18: PFT, ipts, b_inc_tot: ', &
                        j, ipts,b_inc_tot
                   IF(err_act.GT.1)THEN
                      CALL ipslerr_p (3,'growth_fun_all',&
                           'WARNING 18: numerical problem',&
                           'overspending in ordinary allocation','')
                   ENDIF

                ELSE

                   IF (j .EQ. test_pft .AND. printlev_loc.GE.4) THEN

                      ! Succesful allocation
                      WRITE(numout,*) 'Successful allocation'

                   ENDIF

                ENDIF
              
                ! Althought the biomass components respect the allometric 
                ! relationships, there is less than min_stomate carbon left 
                ! to allocate. Put this little carbon in the leaves to 
                ! preserve mass balance closure.
                Cl_inc(1) = Cl_inc(1) + b_inc_tot/circ_class_n(ipts,j,1)
                b_inc_tot = zero
                Cs_inc(1) = zero
                Cr_inc(1) = zero
                Cf_inc(1) = zero

             ELSE

                WRITE(numout,*) 'WARNING 19: Logical flaw'//&
                     'unexpected result in ordinary allocation'
                WRITE(numout,*) 'WARNING 19: PFT, ipts: ', j, ipts
                WRITE(numout,*) 'WARNING 19: ',ABS(Cl_target(1) - Cl(1) ) , Cl(1)
                WRITE(numout,*) 'WARNING 19: ',ABS(Cs_target(1) - Cs(1) )  , Cs(1)
                WRITE(numout,*) 'WARNING 19: ',ABS(Cr_target(1) - Cr(1) )  , Cr(1)
                WRITE(numout,*) 'WARNING 19: ',grow_wood
                WRITE(numout,*) 'WARNING 19: ',b_inc_tot,circ_class_n(ipts,j,1),&
                     b_inc_tot/circ_class_n(ipts,j,1)
                IF(err_act.GT.1)THEN
                   CALL ipslerr_p (3,'growth_fun_all',&
                        'WARNING 19: Logical flaw',&
                        'unexpected result in ordinary allocation','')
                ENDIF
                
             ENDIF ! Ordinary allocation

             !! 5.3.9 Error checking
             IF ( b_inc_tot .GT. min_stomate ) THEN 

                ! This should not happen, in case the functional allocation 
                ! did not consume all the allocatable carbon, the remaining C 
                ! is left for the next day.
                WRITE(numout,*) 'WARNING 20: unexpected outcome force allocation'
                WRITE(numout,*) 'WARNING 20: grow_wood, b_inc_tot: ', grow_wood, b_inc_tot
                WRITE(numout,*) 'WARNING 20: PFT, ipts: ',j,ipts
                IF(err_act.GT.1)THEN
                   CALL ipslerr_p (3,'growth_fun_all',&
                        'WARNING 20: unexpected outcome force allocation','','')
                ENDIF

             ELSEIF ( (b_inc_tot .LT. min_stomate) .AND. &
                  (b_inc_tot .GE. zero) ) THEN

                ! Successful allocation
                ! Debug
                IF (printlev_loc.GE.4 .AND. ipts.EQ.test_grid .AND. j.EQ.test_pft) THEN
                   WRITE(numout,*) 'Successful allocation'
                ENDIF
                !----------

             ELSE

                ! Something possibly important was overlooked
                IF ( (b_inc_tot .LT. zero) .AND. &
                     (b_inc_tot .GE. -100*min_stomate) ) THEN
                   IF (printlev_loc.GE.4 .AND. ipts.EQ.test_grid .AND. j.EQ.test_pft) THEN
                      WRITE(numout,*) 'Marginally successful allocation - '//&
                           'precision is better than 10-6', j
                   ENDIF
                ELSE
                   WRITE(numout,*) 'WARNING 21: Logical flaw '//&
                        'unexpected result in ordinary allocation'
                   WRITE(numout,*) 'WARNING 21: b_inc_tot',b_inc_tot
                   WRITE(numout,*) 'WARNING 21: PFT, ipts: ',j,ipts
                   CALL ipslerr_p (3,'growth_fun_all',&
                        'WARNING 21: Logical flaw unexpected result',&
                        'in ordinary allocation','')
                ENDIF

             ENDIF

             ! The second problem we need to catch is when one of the increment 
             ! pools is negative. This is an undesired outcome (see comment where 
             ! ::KF_old is calculated in this routine. In that case we write a 
             ! warning, set all increment pools to zero and try it again at the 
             ! next time step. A likely cause of this problem is a too large change 
             ! in KF from one time step to another. Try decreasing the acceptable 
             ! value for an absolute increase in KF.
             IF (Cs_inc(1) .LT. zero .OR. & 
                Cr_inc(1) .LT. zero .OR. &
                Cs_inc(1) .LT. zero) THEN
             
                ! Do not allocate - save the carbon for the next time step
                Cl_inc(1) = zero
                Cr_inc(1) = zero
                Cs_inc(1) = zero
                WRITE(numout,*) 'WARNING 22: numerical problem, one of the increment '//&
                     'pools is less than zero'
                WRITE(numout,*) 'WARNING 22: PFT, ipts: ',j,ipts
               
             ENDIF

             !! 5.3.10 Wrap-up phenological and ordinary allocation
             Cl_inc(1) = Cl_inc(1) + Cl_incp(1)
             Cr_inc(1) = Cr_inc(1) + Cr_incp(1)
             Cs_inc(1) = Cs_inc(1) + Cs_incp(1)
             residual(ipts,j) = b_inc_tot

             ! Debug
             IF (printlev_loc.GE.4 .AND. ipts.EQ.test_grid .AND. j.EQ.test_pft) THEN
                WRITE(numout,*) 'Final allocation', ipts, j 
                WRITE(numout,*) 'Cl, Cs, Cr', Cl(1), Cs(1), Cr(1) 
                WRITE(numout,*) 'Cl_incp, Cs_incp, Cr_incp, ', &
                     Cl_incp(1), Cs_incp(1), Cr_incp(1)
                WRITE(numout,*) 'Cl_inc, Cs_ins, Cr_inc, Cf_inc, ', &
                     Cl_inc(1), Cs_inc(1), Cr_inc(1), Cf_inc(1)
                WRITE(numout,*) 'unallocated/residual, ', b_inc_tot
             ENDIF
             !-

             !! 5.3.11 Account for the residual
             !  The residual is usually around ::min_stomate but we deal 
             !  with it anyway to make sure the mass balance is closed
             !  and as a way to detect errors. Move the unallocated carbon
             !  back into the labile pool
             IF (circ_class_biomass(ipts,j,1,ilabile,icarbon) + &
                  residual(ipts,j) .LE. min_stomate) THEN

                deficit = circ_class_biomass(ipts,j,1,ilabile,icarbon) + residual(ipts,j)

                ! The deficit is less than the carbon reserve
                IF (-deficit .LE. circ_class_biomass(ipts,j,1,icarbres,icarbon)) THEN

                   ! Pay the deficit from the reserve pool
                   circ_class_biomass(ipts,j,1,icarbres,icarbon) = &
                        circ_class_biomass(ipts,j,1,icarbres,icarbon) + deficit
                   circ_class_biomass(ipts,j,1,ilabile,icarbon)  = &
                        circ_class_biomass(ipts,j,1,ilabile,icarbon) - deficit

                ELSE

                   ! Not enough carbon to pay the deficit
                   ! There is likely a bigger problem somewhere in
                   ! this routine
                   WRITE(numout,*) 'WARNING 23: PFT, ipts: ',j,ipts
                   CALL ipslerr_p (3,'growth_fun_all',&
                        'WARNING 23: numerical problem overspending ',&
                        'when trying to account for unallocatable C ','')

                ENDIF

             ELSE
                
                ! Move the unallocated carbon back into the labile pool
                circ_class_biomass(ipts,j,1,ilabile,icarbon) = &
                     circ_class_biomass(ipts,j,1,ilabile,icarbon) + residual(ipts,j)
                
             ENDIF

       
             !! 5.3.12 Standardise allocation factors
             !  Strictly speaking the allocation factors do not need to be 
             !  calculated because the functional allocation scheme allocates 
             !  absolute amounts of carbon. Hence, Cl_inc could simply be added to
             !  tmp_bm(:,:,ileaf,icarbon), Cr_inc to tmp_bm(:,:,iroot,icarbon),
             !  etc. However, using allocation factors bears some elegance in 
             !  respect to distributing the growth respiration if this would be 
             !  required. Further it facilitates comparison to the resource 
             !  limited allocation scheme (stomate_growth_res_lim.f90) and it 
             !  comes in handy for model-data comparison. This allocation
             !  takes place at the tree level - note that ::biomass is the only 
             !  prognostic variable from the tree-based allocation
                          
             !  Allocation   
             Cl_inc(1) = MAX(zero, circ_class_n(ipts,j,1) * Cl_inc(1))
             Cr_inc(1) = MAX(zero, circ_class_n(ipts,j,1) * Cr_inc(1))
             Cs_inc(1) = MAX(zero, circ_class_n(ipts,j,1) * Cs_inc(1))
             Cf_inc(1) = MAX(zero, circ_class_n(ipts,j,1) * Cf_inc(1))
             
             ! Total_inc is based on the updated Cl_inc, Cr_inc, 
             ! Cs_inc and Cf_inc. Therefore, do not multiply
             ! ind(ipts,j) again
             total_inc = (Cf_inc(1) + Cl_inc(1) + Cs_inc(1) + Cr_inc(1))
             
             ! Relative allocation
             IF ( total_inc .GT. min_stomate ) THEN

                Cl_inc(1) = Cl_inc(1) / total_inc
                Cs_inc(1) = Cs_inc(1) / total_inc
                Cr_inc(1) = Cr_inc(1) / total_inc
                Cf_inc(1) = Cf_inc(1) / total_inc

             ELSE

                bm_alloc_tot(ipts,j) = zero
                Cl_inc(1) = zero
                Cs_inc(1) = zero
                Cr_inc(1) = zero
                Cf_inc(1) = zero

             ENDIF

             !! 5.3.13 Convert allocation to allocation facors
             !  Convert allocation of individuals to ORCHIDEE's allocation 
             !  factors - see comment for 5.2.5
             !  Aboveground sapwood allocation is age dependent in trees, 
             !  but there is only aboveground allocation in grasses 
             alloc_sap_above = un

             ! Leaf, wood, root and fruit allocation. Note that the X_inc
             ! are normalized before being used here. Calculate f_alloc(fruit) 
             ! as the residual to enhance mass balance closure
             f_alloc(ipts,j,ileaf) = Cl_inc(1)
             f_alloc(ipts,j,isapabove) = Cs_inc(1)*alloc_sap_above
             f_alloc(ipts,j,isapbelow) = Cs_inc(1)*(1.-alloc_sap_above)
             f_alloc(ipts,j,iroot) = Cr_inc(1)
             f_alloc(ipts,j,ifruit) = Cf_inc(1)

             ! Store f_alloc per circ_class to calculate the allocation
             ! in circ_class_biomass after bm_alloc_tot has been checked
             ! for the N availability. Note that the X_inc
             ! are normalized before being used here. Calculate f_alloc_circ(fruit) 
             ! as the residual to enhance mass balance closure
             f_alloc_circ(ipts,1,ileaf) = Cl_inc(1)
             f_alloc_circ(ipts,1,isapabove) = Cs_inc(1)*alloc_sap_above
             f_alloc_circ(ipts,1,isapbelow) = Cs_inc(1)*(1.-alloc_sap_above)
             f_alloc_circ(ipts,1,iroot) = Cr_inc(1)
             f_alloc_circ(ipts,1,ifruit) = Cf_inc(1)

          ELSEIF (.NOT. is_tree(j)) THEN

             ! The first option is IF ( .NOT. is_tree(j) .AND. &
             ! bm_alloc_tot(ipts,j) .GT. min_stomate ) THEN
             ! If we end up here there is not enough biomass to allocate
             f_alloc(ipts,j,ileaf) = zero
             f_alloc(ipts,j,isapabove) = zero
             f_alloc(ipts,j,isapbelow) = zero
             f_alloc(ipts,j,iroot) = zero
             f_alloc(ipts,j,ifruit) = zero
             residual(ipts,j) = zero

             ! Debug
             IF(printlev_loc.GE.4 .AND. ipts.EQ.test_grid .AND. j.EQ.test_pft) WRITE(numout,*) &
                  'there is no non-tree biomass '//&
                  'to allocate, PFT, ', ipts, j
             !-

          ENDIF ! .NOT. is_tree(j) and there is biomass to allocate (§5.3 - far far up)

          ! Intermediate mass balance check. Note that this part of
          ! the code is in DO-loops over nvm and npts so the
          ! 'ipts' label is used in the mass balance check
          IF(err_act.EQ.4) THEN
             
             ! Reset pool_end
             pool_end(:,:,:) = zero
             
             ! The error check makes use of Cx_inc and Cx_incp so it should be
             ! done in the DO-loop for npts.
             pool_end(ipts,j,icarbon) = pool_end(ipts,j,icarbon) + &
                  bm_alloc_tot(ipts,j) * veget_max(ipts,j)             
             
             ! Check mass balance closure. Between intermediate check 3a/b and 4
             ! allocation factors were calculated but not used. Residuals
             ! was moved back into the labile pool.
             CALL intermediate_mass_balance_check(pool_start, pool_end, circ_class_biomass, &
                  circ_class_n, veget_max, bm_alloc_tot, gpp_daily, atm_to_bm, dt, npts, &
                  resp_maint, resp_growth, check_intern_init, ipts, j, '4', 'ipft')
             
          END IF ! err_act.EQ.4

       ENDDO ! npts

       ! Account for the residual (carbon that could not be allocated
       ! during phenological or ordinary growth) before using bm_alloc_tot.
       bm_alloc_tot(:,j)= bm_alloc_tot(:,j) - residual(:,j)

       ! Debug
       IF (printlev_loc.GE.4 .AND. j.EQ.test_pft) THEN
          WRITE(numout,*) 'Accounted for residual'
          WRITE(numout,*) 'bm_alloc_tot_new, ', test_grid, j, &
               bm_alloc_tot(test_grid,j)
       ENDIF
       !-   
       
       !! 5.4 Quantify and account for nitrogen limitation on growth
       DO ipts = 1 , npts 

          ! This far we calculated how we would like to allocate the available
          ! carbon (::f_alloc) and how much carbon of the available carbon we
          ! can allocate (typically except for some exceptional cases and numerical
          ! residuals). Nothing has yet been allocated. Allocation itself, meaning 
          ! the updating of biomass pools is taken care off in the subsequent parts
          ! of the code. circ_class_biomass contains the latest information for
          ! all the biomass pools. The next section is generic for trees, grasses and
          ! crops so we first have to update the information in the biomass variable
          ! so it can be used later
          !+++CHECK+++
          ! replace by biomass_to_cc_2d
!!$       tmp_bm(:,:,:,:) = cc_to_biomass(npts,nvm,&
!!$            circ_class_biomass(:,:,:,:,:),&
!!$            circ_class_n(:,:,:))
          DO iele = 1,nelements
             DO ipar = 1,nparts
                tmp_bm(ipts,j,ipar,iele) = &
                     SUM(circ_class_biomass(ipts,j,:,ipar,iele)*&
                     circ_class_n(ipts,j,:))
             ENDDO
          ENDDO

          !++++++++++++

          ! Initialize
          deltacn=1.0

          ! Nitrogen cost, given required N for a given allocatable 
          ! biomass C, and an intended leaf CN as N = C * costf / C:N
          ! Note that fcn is not a classic c/n ratio but is the c/n ratio
          ! compared to the c/n ratio for leaves. When bm_supply_n is
          ! calculated this is accounted for through multiplying
          ! with cn_leaf. The unit of costf is gN required per gN in the leaf
          costf = f_alloc(ipts,j,ileaf) + fcn_wood(j) * &
               (f_alloc(ipts,j,isapabove)+f_alloc(ipts,j,isapbelow)) + &
               fcn_root(j) * ( f_alloc(ipts,j,iroot) + f_alloc(ipts,j,ifruit))
          
          ! Debug
          IF (printlev_loc.GE.4 .AND. ipts.EQ.test_grid .AND. j.EQ.test_pft) THEN
             WRITE(numout,*) 'costf, ', ipts,j, costf
             WRITE(numout,*) 'bm_alloc_tot, ',bm_alloc_tot(ipts,j)
             WRITE(numout,*) 'f_alloc, ', f_alloc(ipts,j,ileaf),&
               f_alloc(ipts,j,isapabove), f_alloc(ipts,j,isapbelow),&
               f_alloc(ipts,j,iroot),f_alloc(ipts,j,ifruit)
             WRITE(numout,*) 'fcn, ',fcn_wood(j), fcn_root(j)
          ENDIF
          !-
          
          ! Only check if there is biomass growth
          IF ( costf.GT.min_stomate ) THEN

             ! fraction of labile N allocatable for growth
             ! no growth respiration calculated here!
             n_avail = MAX(tmp_bm(ipts,j,ilabile,initrogen)*0.9,0.0)

             ! carbon growth possible given nitrogen availability and 
             ! current nitrogen concentration
             bm_supply_n = n_avail / costf / (1.-frac_growthresp_dyn) * &
                  cn_leaf(ipts,j) 

             ! elasticity of leaf nitrogen concentration
             ! deltacnmax=exp(-(1./(cn_leaf(ipts,j)*0.5*(1./cn_leaf_max(j)+1./&
             ! cn_leaf_min(j))))**8)
             deltacnmax = 1. - exp(-((1.6 * MIN((1./cn_leaf(ipts,j))-(1./cn_leaf_min_2D(ipts,j)),0.) / &
                  ( (1./(cn_leaf_max_2D(ipts,j))) - (1./cn_leaf_min_2D(ipts,j)) ) )**4.1))

             ! Debug
             IF (printlev_loc.GE.3) THEN
                IF((test_grid == ipts).AND.(test_pft==j)) THEN 
                   WRITE(numout,*) 'cn_leaf: ', cn_leaf(ipts,j)
                   WRITE(numout,*) 'bm_supply_n: ', bm_supply_n
                   WRITE(numout,*) 'bm_alloc_tot: ', bm_alloc_tot(ipts,j)
                ENDIF
             ENDIF
             !-
            
             ! Check whether we can allocate all the carbon or whether we
             ! have N-limitation.
             IF ( bm_alloc_tot(ipts,j) .GT. bm_supply_n ) THEN

                IF (impose_cn) THEN
                
                   ! Debug
                   IF (printlev_loc.GE.4) THEN
                      IF((test_grid == ipts) .AND. (test_pft==j)) THEN
                         WRITE(numout,*) "atm_to_bm(initrogen)",atm_to_bm(ipts,j,initrogen)
                         WRITE(numout,*) 'bm_alloc_tot, ',bm_alloc_tot(ipts,j)
                         WRITE(numout,*) "bm_supply_n:",bm_supply_n
                         WRITE(numout,*) 'frac_growthresp_dyn,  ', &
                              1.-frac_growthresp_dyn
                         WRITE(numout,*) 'biomass ilabile N L4090, ', j,test_pft, &
                              tmp_bm(test_grid,test_pft,ilabile,initrogen)
                      ENDIF
                   ENDIF
                   !-
                   
                   ! Impose_cn = y so just take the required nitrogen from
                   ! the atmosphere and add it to the labile pool of the plant
                   atm_to_bm(ipts,j,initrogen) = atm_to_bm(ipts,j,initrogen) + &
                        ((bm_alloc_tot(ipts,j)-bm_supply_n)&
                        *costf*(1.-frac_growthresp_dyn) / &
                        cn_leaf(ipts,j)/0.9)/dt
                   tmp_bm(ipts,j,ilabile,initrogen) = &
                        tmp_bm(ipts,j,ilabile,initrogen) + &
                        (bm_alloc_tot(ipts,j)-bm_supply_n)&
                        *costf*(1.-frac_growthresp_dyn) / &
                        cn_leaf(ipts,j)/0.9
                                    
                   ! tmp_bm(ilabile) has changed, update circ_class_biomass 
                   circ_class_biomass(ipts,j,:,ilabile,initrogen) = &
                        biomass_to_cc(tmp_bm(ipts,j,ilabile,initrogen),&
                        circ_class_biomass(ipts,j,:,ilabile,initrogen),&
                        circ_class_n(ipts,j,:))
                   
                   n_alloc_tot(ipts,j) =  bm_alloc_tot(ipts,j) * &
                        (1.-frac_growthresp_dyn) * costf / cn_leaf(ipts,j)
             
                   ! Debug
                   IF (printlev_loc.GE.4) THEN   
                      IF((test_grid == ipts).AND.(test_pft==j)) THEN
                         WRITE(numout,*) "atm_to_bm(nitrogen)",atm_to_bm(ipts,j,initrogen)
                         WRITE(numout,*) 'biomass ilabile N L4121, ', j,test_pft, &
                              tmp_bm(test_grid,test_pft,ilabile,initrogen)
                      ENDIF
                   ENDIF
                      
                ELSE 

                   !Do not impose_cn thus use the dynamic N-cycle
                   IF (printlev_loc .GE. 4 .AND. ipts.EQ.test_grid .AND. j.EQ.test_pft) THEN
                      WRITE(numout,*) 'N-limitation is a fact'
                   ENDIF
             
                   ! case of not enough nitrogen to sustain intended growth, 
                   ! reduce carbon allocation to meet nitrogen availability
                   ! taking into account the maximal change of nitrogen 
                   ! concentrations. delta of nitrogen concentrations in 
                   ! response to nitrogen deficit
                   deltacnmax=Dmax * (1.-deltacnmax)
                   deltacn = n_avail /  ( bm_alloc_tot(ipts,j) * &
                        (1.-frac_growthresp_dyn) * costf * &
                        1./cn_leaf(ipts,j) ) 
                   deltacn=MIN(MAX(deltacn,1.0-deltacnmax),1.0)
                               
                   ! nitrogen demand given possible nitrogen concentration change
                   n_alloc_tot(ipts,j) =  MIN( n_avail , & 
                        bm_alloc_tot(ipts,j) * (1.-frac_growthresp_dyn) * costf * &
                        MAX(MIN( 1./cn_leaf(ipts,j)*deltacn, 1./cn_leaf_min_2D(ipts,j)), &
                        1./cn_leaf_max_2D(ipts,j)) ) 

                   ! if not successful, reduce growth
                   ! constrain carbon used for growth dependent on available 
                   ! nitrogen under the assumption that the f_allocs are 
                   ! piecewise linear with bm_alloc_tot, which is first-order 
                   ! correct. Remember that at the start of this module we
                   ! took bm_alloc_tot from the labile pool. Put it back, 
                   ! recalculate bm_alloc_tot and than take it back from the 
                   ! labile pool.
                   tmp_bm(ipts,j,ilabile,icarbon) = &
                        tmp_bm(ipts,j,ilabile,icarbon) + &
                        bm_alloc_tot(ipts,j) 
                   bm_alloc_tot(ipts,j) = MIN( bm_alloc_tot(ipts,j) , &
                        n_alloc_tot(ipts,j) / costf / (1.-frac_growthresp_dyn) / &
                        MAX(MIN(1./cn_leaf(ipts,j)*deltacn, &
                        1./cn_leaf_min_2D(ipts,j)), 1./cn_leaf_max_2D(ipts,j)) )

                   ! If bm_alloc_tot did not change, the labile pool will not
                   ! have changed either. In case bm_alloc_tot was adjusted in 
                   ! line with n_alloc_tot (line above), then the excess 
                   ! carbon is stored in the labile pool
                   tmp_bm(ipts,j,ilabile,icarbon) = &
                        tmp_bm(ipts,j,ilabile,icarbon) - &
                        bm_alloc_tot(ipts,j) 

                ENDIF !impose_cn

             ELSE

                IF(printlev_loc.GE.4.AND.ipts.EQ.test_grid.AND.j.EQ.test_pft) THEN
                    WRITE(numout,*) 'Sufficient nitrogen'
                ENDIF
                ! Sufficient nitrogen, increase of leaf nitrogen concentration 
                ! dependent on distance to maximal leaf nitrogen concentration
                ! cannot change leaf C:N in bud burst period nitrogen 
                ! constrained such that nitrogen concentration can only 
                ! increase 1% per day
                deltacnmax=Dmax * deltacnmax
                deltacn = n_avail /  &
                     ( bm_alloc_tot(ipts,j) * (1.-frac_growthresp_dyn) * costf * 1./cn_leaf(ipts,j) ) 
                
                deltacn=MIN(MAX(deltacn,1.0),1.+deltacnmax)

                ! Debug
                IF(printlev_loc.GE.4.AND.test_grid.EQ.ipts.AND.test_pft.EQ.j) THEN
                   WRITE(numout,*) 'biomass ilabile N L4028, ', j,test_pft, &
                        tmp_bm(test_grid,test_pft,ilabile,initrogen)
                ENDIF
                !-

                n_alloc_tot(ipts,j) =  MIN( n_avail , & 
                     bm_alloc_tot(ipts,j) * (1.-frac_growthresp_dyn) * &
                     costf * MAX(MIN(1./cn_leaf(ipts,j)*deltacn, & 
                     1./cn_leaf_min_2D(ipts,j)),1./cn_leaf_max_2D(ipts,j)) )

             ENDIF

          ENDIF ! costf.GT.min_stomate

          !! 5.X Calculate final growth respiration
          !  Since growth respiration was estimated at the start
          !  of this routine, bm_alloc_tot and thus the respiration
          !  associated to the growth may have changed. Here
          !  growth respiration is recalculated. This is the final
          !  calculation. Move the initial estimate back into the
          !  labile pool. Then recalculate resp_growth and finally
          !  take it out of the labile pool again. The labile pool
          !  may have changed so it needs to be updated.
          tmp_bm(ipts,j,ilabile,icarbon) = tmp_bm(ipts,j,ilabile,icarbon) + &
               resp_growth(ipts,j)
                  
          ! Resp_growth may have been set to zero (see exception 25). So use
          ! the lowest estimate for resp_growth
          resp_growth(ipts,j) = MIN(resp_growth(ipts,j), &
               frac_growthresp(j) * bm_alloc_tot(ipts,j))

          ! Take resp_growth from the labile pool
          tmp_bm(ipts,j,ilabile,icarbon) = tmp_bm(ipts,j,ilabile,icarbon) - &
               resp_growth(ipts,j)
                  
          ! tmp_bm(ilabile) has changed, update circ_class_biomass 
          circ_class_biomass(ipts,j,:,ilabile,icarbon) = &
               biomass_to_cc(tmp_bm(ipts,j,ilabile,icarbon),&
               circ_class_biomass(ipts,j,:,ilabile,icarbon),&
               circ_class_n(ipts,j,:))
                  
      ENDDO ! # domain npts

      ! Debug
      IF(printlev_loc.GE.4.AND.j.EQ.test_pft) THEN
         WRITE(numout,*) 'stomate_allocation - bm_alloc_tot may have been adjusted'
         WRITE(numout,*) 'bm_alloc_tot, ',bm_alloc_tot(test_grid,j)
         WRITE(numout,*) 'resp_growth ', resp_growth(test_grid,j)
         IF(bm_alloc_tot(test_grid,j).GT.min_stomate)THEN
            WRITE(numout,*) 'ratio resp_growth/bm_alloc_tot, ', &
                  resp_growth(test_grid,j)/bm_alloc_tot(test_grid,j)
         ENDIF
       ENDIF
       !-


       !! 5.5 Allocate allocatable biomass to different plant compartments   
       !  Absolute allocation at the tree level and for an individual tree (gC tree-1)
       !  The labile and reserve pools are not allocated at the tree level. However,
       !  stand level ilabile and icarbres biomass will be redistributed at the tree
       !  level later in this subroutine. This is done after the relative allocation
       !  beacuse now ::alloc_sap_above is known
       DO ipts = 1,npts
          DO icir = 1,ncirc
             DO ipar = 1,nparts
                IF (ipar.EQ.ileaf .OR. ipar.EQ.isapabove .OR. &
                     ipar.EQ.isapbelow .OR. ipar.EQ.iroot .OR. &
                     ipar.EQ.ifruit) THEN

                   ! Check whether the allocation factor is defined
                   ! and whether there are trees in this circ_class
                   IF (f_alloc_circ(ipts,icir,ipar).GE.zero .AND. &
                        f_alloc_circ(ipts,icir,ipar).LE.un .AND. &
                        circ_class_n(ipts,j,icir).GT.min_stomate) THEN

                      ! Calculate circ_class biomass with the latest and
                      ! final bm_alloc_tot
                      circ_class_biomass(ipts,j,icir,ipar,icarbon) = &
                           circ_class_biomass(ipts,j,icir,ipar,icarbon) + &
                           ( f_alloc_circ(ipts,icir,ipar) * bm_alloc_tot(ipts,j) / &
                           circ_class_n(ipts,j,icir) )

                   ENDIF ! f_alloc is defined
                ENDIF ! select plant parts
             ENDDO ! ipar
          ENDDO ! icir
       ENDDO ! ipts

       !  The amount of allocatable carbon biomass to each compartment is a 
       !  fraction ::f_alloc of the total allocatable biomass
       DO k = 1, nparts

          bm_alloc(:,j,k,icarbon) = f_alloc(:,j,k) * bm_alloc_tot(:,j) 

       ENDDO

       ! All the carbon contained in bm_alloc_tot has been allocated
       bm_alloc_tot(:,j) = zero

       ! Zero the array for PFT 1, since it has not been calculated but it
       ! is used in implict loops below
       bm_alloc_tot(:,1) = zero
       bm_alloc(:,1,:,:) = zero

    ENDDO ! # End Loop over # of PFTs

    ! Intermediate mass balance check. N has not been allocated yet. It is 
    ! tested but the test is not very informative.
    IF (err_act.EQ.4) THEN
          
       ! Reset pool_end
       pool_end(:,:,:) = zero
          
       ! Check mass balance closure. Between intermediate check 4 and 5
       ! the carbon allocation was checked against the available nitrogen
       ! and was finally allocated.
       CALL intermediate_mass_balance_check(pool_start, pool_end, circ_class_biomass, &
            circ_class_n, veget_max, bm_alloc_tot, gpp_daily, atm_to_bm, dt, npts, &
            resp_maint, resp_growth, check_intern_init, ipts, j, '5', 'pft')
       
    ENDIF ! err_act.EQ.4
 
    !! 6. Calculate nitrogen fluxes associated with biomass growth 
    
    ! this is the case of dynamic allocation of nitrogen taken up from the soil
    ! The principles of this allocation are
    ! 1) nitrogen allocated to the plant tissue is dependent on the 
    !    labile nitrogen/carbon ratio i.e. carbon_alloc*(N/C)_labile
    ! 2) the proportions of C/N ratios between different compartments 
    !    are prescribed as in Hybrid 3
    ! 3) since different to Hybrid, grasses have sapwood (...) the 
    !    proportions have be adjusted for grasses, since their tillers 
    !    are not lignified...
    DO j = 2, nvm

       DO ipts = 1, npts

          IF (veget_max(ipts,j) .LE. min_stomate) THEN

             ! This vegetation type is not present, so no reason to do the 
             ! calculation. CYCLE will take us out of the innermost DO loop
             CYCLE
          
          ENDIF

          ! only allocate nitrogen when there is construction of new biomass
          IF(SUM(bm_alloc(ipts,j,:,icarbon)).GT.min_stomate) THEN

             alloc_c(ipts)=0.0
             alloc_d(ipts)=0.0
             alloc_e(ipts)=0.0
             
             ! pool sapwood and roots+fruits into each on pool
             sum_sap(ipts) = bm_alloc(ipts,j,isapabove,icarbon) + &
                  bm_alloc(ipts,j,isapbelow,icarbon)
             sum_oth(ipts) = bm_alloc(ipts,j,iroot,icarbon) + &
                  bm_alloc(ipts,j,ifruit,icarbon)

             IF((sum_sap(ipts)+sum_oth(ipts)).GT.min_stomate) THEN
                
                IF(bm_alloc(ipts,j,ileaf,icarbon).GT.min_stomate) THEN
                 
                   ! in case there is new allocation to leaves
                   alloc_c(ipts) = 1./(1.+(fcn_wood(j)*sum_sap(ipts) + &
                        fcn_root(j)*sum_oth(ipts))/bm_alloc(ipts,j,ileaf,icarbon))
                   alloc_d(ipts) = fcn_wood(j)*sum_sap(ipts) / &
                        bm_alloc(ipts,j,ileaf,icarbon)*alloc_c(ipts)
                   alloc_e(ipts) = fcn_root(j)*sum_oth(ipts) / &
                        bm_alloc(ipts,j,ileaf,icarbon)*alloc_c(ipts)
 
                ELSE

                   ! otherwise add no nitrogen to leaves and alocate the N 
                   ! to whatever is constructed
                   alloc_c(ipts)=0.0
                   alloc_d(ipts)=sum_sap(ipts)/(sum_sap(ipts)+fcn_root(j) / &
                        fcn_wood(j)*sum_oth(ipts))
                   alloc_e(ipts)=1.-alloc_d(ipts)

                ENDIF

             ELSEIF(bm_alloc(ipts,j,ileaf,icarbon).GT.min_stomate)THEN
                
                ! case of only allocation to leaves!
                alloc_c(ipts)=1.0

             ENDIF

             ! Error checking
             ! Sum of the allocation facors c,d and e should not exceed 1
             IF (err_act.EQ.4) THEN
                temp_share = alloc_c(ipts)+alloc_d(ipts)+alloc_e(ipts)
                IF (temp_share-un.GT.zero) THEN
                   ! Total allocation is more than 100% of the available N
                   ! something may have gone wrong.
                   IF (temp_share-un.GT.10000*EPSILON(un)) THEN
                      ! Too large to be a precision error. Fix the problem upstream
                      WRITE(numout,*) 'Overspending nitrogen (%), ', (temp_share-1)*100.
                      CALL ipslerr_p(3,'Sum of alloc_c, alloc_d and alloc_e',&
                           'is more than 100%','N is not preserved','')
                   ELSE
                      ! Looks like a precision error. Fix it.
                      alloc_c(ipts) = alloc_c(ipts)/temp_share
                      alloc_d(ipts) = alloc_d(ipts)/temp_share
                      alloc_e(ipts) = un-alloc_c(ipts)-alloc_d(ipts)
                      ! Check again
                      temp_share = alloc_c(ipts)+alloc_d(ipts)+alloc_e(ipts)
                      IF (temp_share-un.GT.10000*EPSILON(un)) THEN
                         WRITE(numout,*) 'Still overspending nitrogen (%), ', (temp_share-1)*100.
                         CALL ipslerr_p(3,'The fix didn work','revist the fix','','')
                      END IF
                   END IF
                END IF
             END IF
             !-

             !calculate allocation 
             bm_alloc(ipts,j,ileaf,initrogen) = alloc_c(ipts)*n_alloc_tot(ipts,j)

             IF((test_grid == ipts).AND.(test_pft==j).AND. printlev_loc.GE.4)THEN
                WRITE(numout,*) 'alloc_c(ipts)=',alloc_c(ipts) 
                WRITE(numout,*) 'n_alloc_tot(ipts,j)=',n_alloc_tot(ipts,j)
                WRITE(numout,*) 'bm_alloc(ipts,j,ileaf,initrogen)=', &
                     bm_alloc(ipts,j,ileaf,initrogen)
                WRITE(numout,*) 'bm_alloc(ipts,j,ileaf,icarbon)=', &
                     bm_alloc(ipts,j,ileaf,icarbon)
                WRITE(numout,*) 'bm_alloc(ipts,j,:,icarbon)=',&
                     SUM(bm_alloc(ipts,j,:,icarbon))
                IF (bm_alloc(ipts,j,ileaf,initrogen).GT.zero) THEN
                   WRITE(numout,*) "((bm_alloc(ipts,j,ileaf,icarbon)/'&
                        'bm_alloc(ipts,j,ileaf,initrogen)):",&
                        (bm_alloc(ipts,j,ileaf,icarbon)/&
                        bm_alloc(ipts,j,ileaf,initrogen))
                ENDIF
             ENDIF
             !-

             IF(sum_sap(ipts).GT.min_stomate) THEN
                bm_alloc(ipts,j,isapabove,initrogen)=alloc_d(ipts)* & 
                     bm_alloc(ipts,j,isapabove,icarbon)/sum_sap(ipts)* &
                     n_alloc_tot(ipts,j)
                bm_alloc(ipts,j,isapbelow,initrogen) = alloc_d(ipts)* &
                     bm_alloc(ipts,j,isapbelow,icarbon)/sum_sap(ipts)* &
                     n_alloc_tot(ipts,j)
             ENDIF

             IF(sum_oth(ipts).GT.min_stomate) THEN
                bm_alloc(ipts,j,iroot,initrogen) = alloc_e(ipts) * & 
                     bm_alloc(ipts,j,iroot,icarbon)/sum_oth(ipts) * &
                     n_alloc_tot(ipts,j)
                bm_alloc(ipts,j,ifruit,initrogen) = alloc_e(ipts) * & 
                     bm_alloc(ipts,j,ifruit,icarbon)/sum_oth(ipts) * &
                     n_alloc_tot(ipts,j)
             ENDIF

             ! Necessary because bm_alloc_tot can be positive in deciduous, 
             ! but all C is put into the reserve. In that case, all fractions 
             ! are set to zero, thus no nitrogen is allocated alternatively 
             ! formulation is minus (c+d+e)*n_alloc_tot
             n_alloc_tot(ipts,j) = (alloc_c(ipts) + &
                  alloc_d(ipts) + alloc_e(ipts)) * &
                  n_alloc_tot(ipts,j)
                   
          ELSE !IF bm_alloc carbon > min_stomate

             n_alloc_tot(ipts,j) = zero

          ENDIF

          ! Error checking
          ! bm_alloc is the nitrogen that will be allocated. 
          ! It should equal to n_alloc_tot
          IF (err_act.EQ.4) THEN
             temp_share = SUM(bm_alloc(ipts,j,:,initrogen))
             IF (ABS(temp_share-n_alloc_tot(ipts,j)).GT.10000*EPSILON(un)) THEN
                ! Unlikley to be a precision error
                WRITE(numout,*) 'N mass is not preserved, ', &
                     temp_share-n_alloc_tot(ipts,j)
                CALL ipslerr_p(3,'n_alloc_tot does not equal bm_alloc(introgen)', &
                     'This will result in a mass balance problem',&
                     'Fix the problem upstream from the error message','')
             END IF
          END IF
          !-
          
       ENDDO ! # npts
  
    ENDDO ! #PFT

    !! 6.3 Retrieve allocated biomass from labile pool
    !  Only now the allocatable nitrogen is known 
    !  so we will take it from the labile pool.
    tmp_bm(:,:,ilabile,initrogen) = tmp_bm(:,:,ilabile,initrogen) - &
         n_alloc_tot(:,:)

    ! Error checking
    error_count(:,:) = zero
    IF (err_act.EQ.4) THEN
       WHERE (tmp_bm(:,:,ilabile,initrogen).LT.zero)
          error_count(:,:) = un
       END WHERE
       IF (SUM(error_count(:,:)).GT.zero) THEN
          CALL ipslerr_p(3,'Overspending labile nitrogen',&
               'the labile pool is negative','','')
       END IF
       WRITE(numout,*) 'PASSED check: tmp_bm(ilabile) is positive'
    END IF
    !-

    ! Nitrogen allocation is now accounted for in tmp_bm(ilabile,initrogen)
    ! n_alloc_tot could be set to zero
    n_alloc_tot(:,:) = zero

    ! Some temporary variables to simplify the calculations
    tmp_bm(:,:,ileaf,:)=zero
    tmp_bm(:,:,iroot,:)=zero
    DO icir = 1,ncirc
       DO iele = 1,nelements
          tmp_bm(:,:,ileaf,iele) = tmp_bm(:,:,ileaf,iele) + &
               circ_class_biomass(:,:,icir,ileaf,iele) * &
               circ_class_n(:,:,icir)
          tmp_bm(:,:,iroot,iele) = tmp_bm(:,:,iroot,iele) + &
               circ_class_biomass(:,:,icir,iroot,iele) * &
               circ_class_n(:,:,icir)
       ENDDO
    ENDDO

    ! Calculate the nitrogen that is being translocated from the leaves
    ! back to the labile pool. If ORCHIDEE can allocate more than the
    ! current N/C ratio it will irrespective of whether we are below
    ! or above the optimal C/N ratio. If leaf N/C is already low but
    ! bm_all_N/bm_alloc_C is even lower we will take N out of the leaves
    ! and make the N-limitation even worse. Although seems counterintuitive
    ! it will decrease NUE, Vcmax and thus the GPP at the next time step.
    ! transloc is thus a short term control of N-allocation reflecting the
    ! higher reactivity of nitrogen compared to carbon. As a test transloc
    ! was set to zero after it was calculated. This resulted in 0 to 5%
    ! changes in GPP for all PFTs except PFT5 and PFT15. Without transloc
    ! PFT5 did not grow very well. For PFT15 the effect of transloc was
    ! mixed showing pixels where the PFT grew better as well as pixel where
    ! it grew worse than with transloc. transloc seems to be a more
    ! short term response that is more or less trying to do the same as
    ! sugar_load. sugar_load is much slower but also much more intrusive.
    ! Although it appears that this block of code is intended to represent
    ! a real process (and was unlikely added as a fix afterward), it is
    ! explicitly described in Zaehle et al 2010 (incl the SI). 
    
    ! Set to zero for ileaf
    transloc(:,:) = zero
    WHERE((tmp_bm(:,:,ileaf,icarbon).GT.min_stomate).AND.&
         (bm_alloc(:,:,ileaf,icarbon).GT.min_stomate))
       
       transloc(:,:) = tmp_bm(:,:,ileaf,icarbon) * 0.05 * & 
            (bm_alloc(:,:,ileaf,initrogen)/bm_alloc(:,:,ileaf,icarbon) - & 
            tmp_bm(:,:,ileaf,initrogen) / &
            tmp_bm(:,:,ileaf,icarbon))
       transloc(:,:) = MAX(MIN(tmp_bm(:,:,ilabile,initrogen)*0.7,&
            transloc(:,:)), -bm_alloc(:,:,ileaf,initrogen))
       tmp_bm(:,:,ilabile,initrogen) = tmp_bm(:,:,ilabile,initrogen) - &
            transloc(:,:)
       bm_alloc(:,:,ileaf,initrogen) = bm_alloc(:,:,ileaf,initrogen) + &
            transloc(:,:) 

    ENDWHERE

    ! Error checking
    IF (err_act.EQ.4) THEN
       error_count(:,:) = zero
       WHERE (tmp_bm(:,:,ilabile,initrogen).LT.zero)
          error_count(:,:) = un
       ENDWHERE
       IF (SUM(error_count(:,:)).GT.zero) & 
            CALL ipslerr_p(3,'test 1 tranloc exceeds tmp_bm(ilabile)','','','')
       WHERE (bm_alloc(:,:,ileaf,initrogen).LT.zero)
          error_count(:,:) = un
       ENDWHERE
       IF (SUM(error_count(:,:)).GT.zero) & 
            CALL ipslerr_p(3,'test 2 tranloc exceeds bm_alloc(ileaf)','','','')
       WRITE(numout,*) 'PASSED check: transloc leaves OK'
    END IF
    !-
    
    ! Set to zero for iroot
    transloc(:,:) = zero
    
    WHERE((tmp_bm(:,:,iroot,icarbon).GT.min_stomate).AND. &
         (bm_alloc(:,:,iroot,icarbon).GT.min_stomate))
       
       transloc(:,:) = tmp_bm(:,:,iroot,icarbon) * 0.05 * &
            (bm_alloc(:,:,iroot,initrogen)/bm_alloc(:,:,iroot,icarbon) - &
            tmp_bm(:,:,iroot,initrogen) / &
            tmp_bm(:,:,iroot,icarbon))
       transloc(:,:) = MAX(MIN(tmp_bm(:,:,ilabile,initrogen)*0.7,&
            transloc(:,:)), -bm_alloc(:,:,iroot,initrogen))
       tmp_bm(:,:,ilabile,initrogen) = tmp_bm(:,:,ilabile,initrogen) - &
            transloc(:,:)
       bm_alloc(:,:,iroot,initrogen) = bm_alloc(:,:,iroot,initrogen) + &
            transloc(:,:)
    ENDWHERE

    ! Error checking
    IF (err_act.EQ.4) THEN
       WHERE (tmp_bm(:,:,ilabile,initrogen).LT.zero)
          error_count(:,:) = un
       ENDWHERE
       IF (SUM(error_count(:,:)).GT.zero) & 
            CALL ipslerr_p(3,'test 3 tranloc exceeds tmp_bm(ilabile)','','','')
       WHERE (bm_alloc(:,:,iroot,initrogen).LT.zero)
          error_count(:,:) = un
       ENDWHERE
       IF (SUM(error_count(:,:)).GT.zero) & 
            CALL ipslerr_p(3,'test 4 tranloc exceeds bm_alloc(iroot)','','','')
       WRITE(numout,*) 'PASSED check: transloc roots OK'
    END IF
    !-

 !! 7. Update the biomass with newly allocated nitrogen biomass 

    ! For icarbon, circ_class_biomass contains the latest information
    ! for all its pools, For nitrogen this is also the case except for 
    ! the labile pool. Account for the latest changes in tmp_bm(ilabile) 
    ! and update circ_class_biomass so it contains the latest values 
    ! for all pools for both elements.
    circ_class_biomass(:,:,:,ilabile,initrogen) = &
         biomass_to_cc(tmp_bm(:,:,ilabile,initrogen),&
         circ_class_biomass(:,:,:,ilabile,initrogen),&
         circ_class_n(:,:,:),npts,nvm)

    ! Now convert circ_class_biomass to biomass to do all subsequent
    ! calculations. This needs to recalculated because it is possible
    ! that circ_class_biomass has changed since the last time we 
    ! calculated biomass (see 5.4) 
    tmp_bm(:,:,:,:) = cc_to_biomass(npts,nvm,&
         circ_class_biomass(:,:,:,:,:),&
         circ_class_n(:,:,:)) 
    tmp_bm(:,:,:,initrogen) = tmp_bm(:,:,:,initrogen) + &
         bm_alloc(:,:,:,initrogen)

    ! All the biomass pools may have increased for nitrogen so update 
    ! circ_class_biomass. After this operation circ_class_biomass and
    ! biomass are in sync for carbon and nitrogen
    DO ipar = 1,nparts
       circ_class_biomass(:,:,:,ipar,initrogen) = &
            biomass_to_cc(tmp_bm(:,:,ipar,initrogen),&
            circ_class_biomass(:,:,:,ipar,initrogen),&
            circ_class_n(:,:,:),npts,nvm)
    ENDDO

    ! Intermediate mass balance check. Both C and N can be checked.
    IF (err_act.EQ.4) THEN

       ! calculate pool_end
       pool_end(:,:,:) = zero
       
       ! Check mass balance closure. Between intermediate check 5 and 6
       ! n has been allocated to the different pools. Both C and N are now
       ! really being checked.
       CALL intermediate_mass_balance_check(pool_start, pool_end, circ_class_biomass, &
            circ_class_n, veget_max, bm_alloc_tot, gpp_daily, atm_to_bm, dt, npts, &
            resp_maint, resp_growth, check_intern_init, ipts, j, '6', 'pft')
       
    ENDIF ! err_act.EQ.4

    !! 8. Use or fill reserve pools depending on relative size of the labile and reserve C pool

    ! +++ CHECK +++
    ! Externalize all the hard coded values i.e. 0.3
    ! Calculate the labile pool for all plants and also the reserve pool for trees      
    DO j = 2,nvm

       DO ipts = 1,npts

          ! Initialize 
          labile_target(ipts,j,:) = zero 
          reserve_target(ipts,j,:) = zero 

          IF ( (veget_max(ipts,j) .LE. min_stomate) .OR. &
               (SUM(tmp_bm(ipts,j,:,icarbon)) .LT. min_stomate) ) THEN

             ! This vegetation type is not present, so no reason to do the 
             ! calculation. CYCLE will take us out of the innermost DO loop
             CYCLE

          ENDIF

          !! 8.2 Calculate the reserves
          !  There is vegetation present and it has started growing. The second and 
          !  third condition required to make the PFT survive the first year during 
          !  which the long term climate variables are initialized for the phenology. 
          !  If these conditions are not added, the reserves are respired well 
          !  before growth ever starts
          IF ( veget_max(ipts,j) .GT. min_stomate .AND. &
               rue_longterm(ipts,j) .GE. zero .AND. &                  
               rue_longterm(ipts,j) .NE. un) THEN

             !! 8.3 Calculate the optimal size of the pools   
             ! We had an endless series of problems which were often difficult to
             ! understand but which always seemed to be related to a sudden drop
             ! in tmp_bm(ilabile). This drop was often the result of a sudden 
             ! change in labile_target. Given that there is not much science behind
             ! this approach it seems a good idea to remove this max statement to
             ! avoid sudden changes. Rather than using the actual biomass we propose
             ! to use the target biomass. This assumes that the tree would like to 
             ! fill its labile pool to be optimal when it would be in allometric
             ! balance.
             IF (is_tree(j)) THEN

                ! We will make use of the actual sapwood, heartwood and effective height
                ! and then calculate the target leaves and roots. This approach gives
                ! us a target for a labile_target of a tree in allometric balance. 
                ! Basal area at the tree level (m2 tree-1)
                circ_class_ba_eff(:) = &
                     wood_to_ba_eff(circ_class_biomass(ipts,j,:,:,icarbon),&
                     j, pipe_tune2(ipts,j))             

                ! Current biomass pools per tree (gC tree^-1) 
                ! We will have different trees so this has to be calculated from the 
                ! diameter relationships            
                Cs(:) = (circ_class_biomass(ipts,j,:,isapabove,icarbon) + &
                     circ_class_biomass(ipts,j,:,isapbelow,icarbon)) * scal(ipts,j)
                Cr(:) = (circ_class_biomass(ipts,j,:,iroot,icarbon)) * scal(ipts,j)
                Cl(:) = (circ_class_biomass(ipts,j,:,ileaf,icarbon)) * scal(ipts,j)

                DO l = 1,ncirc 

                   !  Calculate tree height
                   circ_class_height_eff(l) = pipe_tune2(ipts,j) * &
                        (4/pi*circ_class_ba_eff(l))**(pipe_tune3(j)/2)       

                ENDDO

                ! Use the pipe model to calculate the target leaf and root
                ! biomasses
                DO l = 1,ncirc
                   IF(circ_class_height_eff(l) .GT. min_stomate) THEN
                      Cl_target(l) = KF(ipts,j) * Cs(l) / circ_class_height_eff(l)

                   ELSE
                      Cl_target(l) = MAX((KF(ipts,j) * Cs(l) / circ_class_height_eff(l)), &
                           (MAX(Cs(1)*KF(ipts,j), Cr(1)*LF(ipts,j), Cl(1))))

                   ENDIF
                ENDDO

                DO l = 1,ncirc
                   Cr_target(l) = Cl_target(l) / LF(ipts,j)
                ENDDO
                
             ELSEIF ( .NOT. is_tree(j)) THEN

                ! grasses and crops
                ! Initialize
                Cs(:) = zero
                Cl_target(:) = zero 
                Cr_target(:) = zero

                ! Current biomass pools per grass/crop (gC ind^-1)
                ! Cs has too many dimensions for grass/crops. To have a consistent 
                ! notation the same variables are used as for trees but the dimension 
                ! of Cs, Cl and Cr i.e. ::ncirc should be ignored            
                Cs(1) = tmp_bm(ipts,j,isapabove,icarbon) * scal(ipts,j)

                ! Use the pipe model to calculate the target leaf and root
                ! biomasses
                Cl_target(1) = Cs(1) * KF(ipts,j)
                Cr_target(1) = Cl_target(1) / LF(ipts,j)

             ENDIF !is_tree

             ! Accounting for the N-concentration of the tissue as a proxy
             ! of tissue activity. There were some problems for the deciduous trees, when the 
             ! targeted labile pool was calculated based on the targeted values from the allometric 
             ! allocation (see previous versions). There was an inconsistency in the units, Therefore
             ! a corrected approach has been introduced. With the corrected approach the labile_target
             ! typically exceeds 10*gpp_week during the growing season. With a too low labile_target
             ! the feedback through sugar_load was way too strong. 
             labile_target(ipts,j,icarbon)=gtemp(ipts,j)*labile_reserve(j)*(tmp_bm(ipts,j,ileaf,initrogen)+ &
                                           tmp_bm(ipts,j,iroot,initrogen) + tmp_bm(ipts,j,ifruit,initrogen) + &
                                           tmp_bm(ipts,j,isapabove,initrogen)+ tmp_bm(ipts,j,isapbelow,initrogen))      
             labile_target(ipts,j,icarbon) = MAX ( labile_target(ipts,j,icarbon), 10. * gpp_week(ipts,j) )

             ! Debug 
             IF (printlev_loc.GE.4 .AND. ipts.EQ.test_grid .AND. j.EQ.test_pft) THEN
                WRITE(numout,*) 'circ_class_biomass(ipts,j,:,isapbelow,icarbon)=', &
                     circ_class_biomass(ipts,j,:,isapbelow,icarbon)
                WRITE(numout,*) 'circ_class_biomass(ipts,j,:,isapabove,icarbon)=', &
                     circ_class_biomass(ipts,j,:,isapabove,icarbon)
                WRITE(numout,*) 'scal=',scal(ipts,j)
                WRITE(numout,*) 'labile_reserve(j)=',labile_reserve(j)
                WRITE(numout,*) 'Cl_target(:)=',Cl_target(:)
                WRITE(numout,*) 'Cr_target(:)=',Cr_target(:)
                WRITE(numout,*) 'Cs(:)=',Cs(:)
                WRITE(numout,*) 'fcn_root(j)=',fcn_root(j)
                WRITE(numout,*) 'fcn_wood(j)=',fcn_wood(j)
                WRITE(numout,*) 'cn_leaf(ipts,j)=',cn_leaf(ipts,j)
             ENDIF
             !-

             ! The max size of reserve pool is proportional to the size of the 
             ! storage organ (the sapwood) and a the leaf functional trait of the 
             ! PFT (::phene_type_tab). The reserve pool is constrained by the mass 
             ! needed to replace foliage and roots. This constraint prevents the
             ! scheme from putting too much reserves in big trees (which have a lot 
             ! of sapwood compared to small trees). Exessive storage would hamper 
             ! tree growth and would make mortality less likely.
             IF(is_tree(j)) THEN

                IF (pheno_type(j).EQ.1) THEN 

                   ! Evergreen trees are not very conservative with respect to 
                   ! C-storage. Therefore, only 5% of their sapwood mass is stored 
                   ! in their reserve pool.
                   reserve_target(ipts,j,icarbon) = MIN(evergreen_reserve(j) * &
                        ( tmp_bm(ipts,j,isapabove,icarbon) + &
                        tmp_bm(ipts,j,isapbelow,icarbon)), &
                        lai_to_biomass(lai_target(ipts,j),j)*&
                        (1.+root_reserve(j)/ltor(ipts,j)))

                   ! Debug
                   IF (printlev_loc.GE.4 .AND. ipts.EQ.test_grid .AND. j.EQ.test_pft) THEN
                      WRITE(numout,*) 'What happens to the reserve and labile &
                           & pools? Evergreen'
                      WRITE(numout,*) 'carbres, reserve_target: ',&
                           tmp_bm(ipts,j,icarbres,icarbon),reserve_target(ipts,j,icarbon)
                      WRITE(numout,*) 'ilabile, labile_target: ',&
                           tmp_bm(ipts,j,ilabile,icarbon),labile_target(ipts,j,icarbon)
                      WRITE(numout,*) 'evergreen_reserve(j): ',&
                           evergreen_reserve(j)
                      WRITE(numout,*) 'isapabove,isapbelow: ',&
                           tmp_bm(ipts,j,isapabove,icarbon), &
                           tmp_bm(ipts,j,isapbelow,icarbon)
                   ENDIF
                   !-

                ELSE

                   ! Deciduous trees are more conservative and 12% of their sapwood mass 
                   ! is stored in the reserve pool. The scheme avoids that during the 
                   ! growing season too much reserve are accumulated (which would hamper 
                   ! growth), therefore, the reduced rate of 12% is used until scenecence.
                   IF (SUM(bm_alloc(ipts,j,:,icarbon)) .GT. min_stomate) THEN

                      reserve_target(ipts,j,icarbon) = deciduous_reserve(j) * &
                           ( tmp_bm(ipts,j,isapabove,icarbon) + & 
                           tmp_bm(ipts,j,isapbelow,icarbon) )

                      ! Debug
                      IF (printlev_loc.GE.4 .AND. ipts.EQ.test_grid .AND. j.EQ.test_pft) THEN
                         WRITE(numout,*) 'What happens to the reserve and '//&
                              'labile pools? Deciduous'
                         WRITE(numout,*) 'carbres, reserve_target: ',&
                              tmp_bm(ipts,j,icarbres,icarbon),reserve_target(ipts,j,icarbon)
                         WRITE(numout,*) 'deciduous_reserve: ',&
                              deciduous_reserve(j)
                         WRITE(numout,*) 'isapabove,isapbelow: ',&
                              tmp_bm(ipts,j,isapabove,icarbon), &
                              tmp_bm(ipts,j,isapbelow,icarbon)
                      ENDIF
                      !-

                   ELSE  

                      ! If the plant is senescent, allow for a higher reserve mass. Plants 
                      ! can then use the excess labile C, that is no longer used for growth 
                      ! and would be respired otherwise, to regrow leaves after the dormant 
                      ! period. This code is a more stable alternative by Nicolas
                      reserve_target(ipts,j,icarbon) = senescense_reserve(j) * &
                           ( tmp_bm(ipts,j,isapabove,icarbon) + & 
                           tmp_bm(ipts,j,isapbelow,icarbon))

                      ! Debug
                      IF (j .EQ. test_pft .AND. printlev_loc.GE.4 .AND. &
                           ipts == test_grid) THEN
                         WRITE(numout,*) 'What happens to the reserve and labile pools? '//&
                              'Senescent'
                         WRITE(numout,*) 'carbres, reserve_target: ',&
                              tmp_bm(ipts,j,icarbres,icarbon),reserve_target(ipts,j,icarbon)
                         WRITE(numout,*) 'senescense_reserve(j): ',&
                              senescense_reserve(j)
                         WRITE(numout,*) 'isapabove,isapbelow: ',&
                              tmp_bm(ipts,j,isapabove,icarbon), &
                              tmp_bm(ipts,j,isapbelow,icarbon)
                      ENDIF

                   ENDIF ! Scenecent

                ENDIF ! Phenology type

                
             ELSE

                ! Grasses
                ! The min criterion results in the reserves being zero because 
                ! isapabove goes to zero when the reserves are most needed. Use
                ! lai_to_biomass to account for a dynamic sla. Some high latitude
                ! pixels were found to have very low biomasses which result in 
                ! a zero lai_target (10e-15) which in turn results for a reserve_target
                ! of zero. This may cause problems later on.
                reserve_target(ipts,j,icarbon) = &
                     MIN(deciduous_reserve(j) * (tmp_bm(ipts,j,iroot,icarbon) + &
                     tmp_bm(ipts,j,isapabove,icarbon) + & 
                     tmp_bm(ipts,j,isapbelow,icarbon)), &
                     lai_to_biomass(lai_target(ipts,j),j) * &
                     (1.+root_reserve(j)/ltor(ipts,j)))

                ! Debug 
                IF (j .EQ. test_pft .AND. printlev_loc.GE.4 .AND. ipts == test_grid) THEN
                   WRITE(numout,*) 'reserve target, ', reserve_target(ipts,j,icarbon)
                ENDIF
                !-

             ENDIF

             !! 8.5 Move carbon between the reserve and labile pools
             !  Fill the reserve pools up to their optimal level or until the min/max 
             !  limits are reached. The original approach in OCN resulted in 
             !  instabilities and sometimes oscilations. For this reason a more 
             !  simple and straightforward transfer between the pools has been 
             !  implemented. After sugar loading was implemeted this simplified
             !  approach was found to come with some sudden regime switches because
             !  C was only transfered from one pool to another if one pool was full.
             !  The latest approach tries to fill both pools at the same time but 
             !  with a different (arbitrary speed). At present there are only two 
             !  important differences between the reserve and the labile pool in 
             !  ORCHIDEE: (1) only the labile pool is used in the allocation (hence,
             !  we should try to store as much carbon as possible in the labile pool) 
             !  and (2) the reserve pool comes without autotrophic respiration (hence,
             !  we should try to store as much carbon in the reserve pool as possible
             !  because that will enable us to grow more leaves in spring). These
             !  conflicting considerations were translated in the idea that
             !  reserve_target X**2 + labile_target X = total reserves. Where X is the share
             !  of the labile pool in the total reserves (= labile + reservers) In other 
             !  words the labile pool gets priority over the reserve pool but the model 
             !  will try to fill both pools at the same time. The optimal pools can 
             !  therefore be calculated by solving a simple quadratic eaqutions: 
             !  X = -b - sqrt(b² - 4ac)/(2a) where b = labile_target, and 
             !  c = tmp_bm(labile+reserve).
             !  Another way of looking at this approach is assuming that the labile pool
             !  fills up linearly and the reserve pool quadratically. As long as the labile
             !  pool is below its optimal, it will thus have priority. As soon as the labile
             !  pool is reached up to its target value (=considered the optimum) the reserve
             !  pool will start taking in more carbon. For simplicity it was assumed that 
             !  carbon can move freely between both pools. 
             total_reserves = tmp_bm(ipts,j,icarbres,icarbon) + tmp_bm(ipts,j,ilabile,icarbon)
             optimal_share = zero

             ! Only calculate a solution if there is carbon in the reserve pool
             IF (total_reserves.GT.min_stomate) THEN

                ! Avoid divide by zero in case the reserve_target = zero
                IF (reserve_target(ipts,j,icarbon).GT.min_stomate) THEN

                   ! First solution of a quadratic equation. Note: the quadratic solution
                   ! is based on aX2+bX+c=0. c is thus -total_reserves 
                   optimal_share =(-labile_target(ipts,j,icarbon) + SQRT(labile_target(ipts,j,icarbon)**2 + &
                           4*reserve_target(ipts,j,icarbon)*total_reserves)) / (2*reserve_target(ipts,j,icarbon))

                   IF (optimal_share.LT.zero) THEN
                      
                      ! Second solution of a quadartic equation in case the first solution
                      ! turned out to be the negative solution which we don want.
                      optimal_share =(-labile_target(ipts,j,icarbon) - SQRT(labile_target(ipts,j,icarbon)**2 + &
                           4*reserve_target(ipts,j,icarbon)*total_reserves)) / (2*reserve_target(ipts,j,icarbon))
                   ENDIF

                END IF
                
             END IF

             ! Error checking
             IF (optimal_share.LT.zero) THEN
                WRITE(numout,*) 'Tried both roots of the quadratic equation and both were negative', &
                     optimal_share
                CALL ipslerr_p(3,'growth_fun_all', 'Both solutions are negative',&
                     'something must be wrong with the calculation of the labile', &
                     'and reserve pools')
             END IF
             !-

             ! Calculate the optimal distribution between the reserve and the
             ! labile pool. Assume full mobility between both pools and ignore
             ! possible constraints during dormacy. Calculate ilabile as the 
             ! residual term to avoid mass balance problems (the alternative way
             ! to calculate it is tmp_bm(ilabile) = labile_target*optimal_share.     
             ! Optimal share was calculated with the equation above.
             ! This solution can result in precision issues and it can deplete the
             ! carbres pool. Add a safety valve.
             tmp_bm(ipts,j,icarbres,icarbon) = MAX(total_reserves * 0.1, &
                  MIN(reserve_target(ipts,j,icarbon)*(optimal_share**2), 0.9 * total_reserves))
             tmp_bm(ipts,j,ilabile,icarbon) = MAX(zero,total_reserves - &
                  tmp_bm(ipts,j,icarbres,icarbon))

             ! For evergreen PFTs with total reserves well below the optimal the code 
             ! above may results in reducing the reserves pools every day. If this
             ! happens, the labile carbon becomes e-68 within 3 years and soon after 
             ! causes an overflow error with e-302. The lines below should avoid this
             ! to happen.
             total_reserves = tmp_bm(ipts,j,icarbres,icarbon) + tmp_bm(ipts,j,ilabile,icarbon)
             IF (total_reserves .LT. 1000 * EPSILON(zero)) THEN
                tmp_bm(ipts,j,icarbres,icarbon) = zero
                tmp_bm(ipts,j,ilabile,icarbon) = zero
                circ_class_biomass(ipts,j,1,isapabove,icarbon) = &
                     circ_class_biomass(ipts,j,1,isapabove,icarbon) + &
                     total_reserves
             END IF

             ! Error checking
             IF (tmp_bm(ipts,j,ilabile,icarbon).LT.zero .OR. &
                  tmp_bm(ipts,j,icarbres,icarbon).LT.zero) THEN
                WRITE(numout,*) 'The reserve pool is negative after &
                     & re-allocation.  Not good!'
                WRITE(numout,*) 'tmp_bm(ipts,j,icarbres,icarbon): ',&
                     tmp_bm(ipts,j,icarbres,icarbon)
                WRITE(numout,*) 'ipts,j : ',ipts,j
                WRITE(numout,*) 'total_reserves : ',total_reserves
                WRITE(numout,*) 'optimal_share : ',optimal_share
                CALL ipslerr_p(3,'growth_fun_all', 'Negative reserves after re-distruting',&
                     'carbon between the labile and reserve pools.  Not sure', &
                     'how this happened.')
             ENDIF
             !-

          ELSEIF ( veget_max(ipts,j) .GT. min_stomate .AND. &
               rue_longterm(ipts,j) .EQ. un) THEN

             ! There hasn't been any photosynthesis yet. This happens when a 
             ! new vegetation is prescribed and the longterm phenology variables 
             ! are not initialized yet. These conditions happen when the model is 
             ! started from scratch (no restart files). Because the plants are 
             ! very small, they contain little reserves. We increased the amount 
             ! of reserves by a factor ::tune_r_in_sapling where r stands for 
             ! reserves. However, this amount gets simply respired before it is 
             ! needed because the reserve_target is calculated as a function of the 
             ! sapwood biomass which is very low because the plants are really 
             ! small. Here we skip recalculating the reserve_target until the day 
             ! we start using it.

          ELSE

             ! No reason to be here
             WRITE(numout,*) 'Error: unexpected condition for the reserve pools, pft, ',j
             WRITE(numout,*) 'veget_max, rue_longterm, ', veget_max(ipts,j), &
                  rue_longterm(ipts,j)

          ENDIF ! rue_longterm

          ! The model code does not control the C/N ratio of the labile pool hence,  
          ! even if there is a strong N-limitation, the model can accumulate lots  
          ! of carbon in the labile pool. The first CN-version was indeed doing this  
          ! the plant could easily store several 1000 gC m-2. As this was considered  
          ! unrealistic the excess C in the labile pool was burned-off by some excess  
          ! respiration. Although this luxury/wastage respiration has been suggested  
          ! in the literature (see Amthor et al 2000 and Chamber et al 2004) it is  
          ! not confirmed by many observations. We first tried to control the C/N ratio  
          ! of the labile pool but ran into several numerical issues with small numbers  
          ! and some of the dynamics of the pool during phenology and senescence  
          ! (N-resorption). It was then decided to simply control the size of the  
          ! labile pool. The model already had an estimate of the optimal pool size of  
          ! the labile and carbres pools. If the plant has more labile carbon than the  
          ! optimal, GPP is downregulated (too much sugars in the leaves will increase  
          ! the viscosity and hamper the sapflow in the phloem. The viscosity can be  
          ! decreased again by closing the stomata and transpiring less of the sapflow  
          ! in the xylem. By closing the stomata, GPP will be downregulated. See Holtta  
          ! et al 2017). Because ORCHIDEE has no sapflow, turgor and viscosity yet, we  
          ! used a simple ratio to downregulate NUE. The regulation is smoothened by  
          ! setting boundaries to avoid sudden decreases in GPP (which are not apparent  
          ! in the data). Smoothing is taken care of in stomate_vmax.f90. If the plant  
          ! has less carbon in its labile and carbres pools than wanted, the NUE is  
          ! upregulated. Up regulation is also capped to avoid crazy NUE values and high  
          ! frequency changes between up and downregulation. Up and downregulation are 
          ! done in stomate_vmax.f90.
          ! CYmark: we think the active regulation of sugar load only occurs in active growth
          ! stage, i.e., when plant_status is equal to icanopy.
          IF (tmp_bm(ipts,j,ilabile,icarbon)+tmp_bm(ipts,j,icarbres,icarbon).GT.zero .AND.&  
               plant_status(ipts,j).EQ.icanopy) THEN

             ! First priority: too much sugar in the labile pool-> downregulate
             update_sugar_load = (labile_target(ipts,j,icarbon)+reserve_target(ipts,j,icarbon)) / &
                  (tmp_bm(ipts,j,ilabile,icarbon) + tmp_bm(ipts,j,icarbres,icarbon))
             
             sugar_load(ipts,j) = (sugar_load(ipts,j) * (tau_sugarload_week - dt ) + & 
                  max(sugar_load_min,min(update_sugar_load,sugar_load_max)) &
                  * dt) /tau_sugarload_week
             
          ELSEIF (tmp_bm(ipts,j,ileaf,icarbon).GT.zero) THEN
             sugar_load(ipts,j)  = (sugar_load(ipts,j) * (tau_sugarload_week - dt ) + & 
                  sugar_load_max * dt) /tau_sugarload_week
          ELSE
             ! Out of growing season, too much labile but not enough reserves, etc
             ! -> do nothing
             sugar_load(ipts,j) = un             
          ENDIF
          

          !! 8.1 Calculate NPP
          !  Calculate the NPP @tex $(gC.m^{-2}dt^{-1})$ @endtex as the difference 
          !  between GPP and the two components of autotrophic respiration 
          !  (maintenance and growth respiration). GPP, R_maint and R_growth 
          !  are prognostic variables, NPP is calculated as the residuals and is 
          !  thus a diagnostic variable. Note that NPP is not used in the 
          !  allocation scheme, instead bm_alloc_tot is allocated. The 
          !  physiological difference between both is that bm_alloc_tot does no 
          !  longer contain the reserves and labile pools and is only the carbon 
          !  that needs to go into the biomass pools. NPP contains the reserves
          !  and labile carbon. Note that GPP is in gC m-2 s-1 whereas the 
          !  respiration terms were calculated in gC m-2 dt-1
          npp(ipts,j) = gpp_daily(ipts,j) - resp_growth(ipts,j)/dt - &
               resp_maint(ipts,j)/dt 

          !! 8.6 Use or fill reserve pools depending on relative size of the 
          !  labile and reserve N pool
          IF(veget_max(ipts,j) .GT. min_stomate) THEN 

             costf = f_alloc(ipts,j,ileaf) + fcn_wood(j) * &
                  (f_alloc(ipts,j,isapabove)+f_alloc(ipts,j,isapbelow)) + &
                  fcn_root(j) * ( f_alloc(ipts,j,iroot) + f_alloc(ipts,j,ifruit) )

             IF (costf.EQ.0.0) costf=1.

             ! In case the PFT is planted for the first time through LCC cn_leaf
             ! is undefined. So we use the fixed parameter cn_leaf_min_2D instead.
             IF (cn_leaf(ipts,j).LE.zero) THEN
                cn_leaf(ipts,j) = cn_leaf_min_2D(ipts,j)
             END IF
             
             !+++CHECK+++
             ! Should we calculate the target for labile nitrogen based on the
             ! target labile carbon or the actual labile carbon?
!!$             labile_target(ipts,j,initrogen) = tmp_bm(ipts,j,ilabile,icarbon) / &
!!$                  cn_leaf(ipts,j) * costf
             labile_target(ipts,j,initrogen) = labile_target(ipts,j,icarbon) / &
                     cn_leaf(ipts,j) * costf
             !+++++++++++

             ! excess or deficit of nitrogen in the labile pool
             use_lab = tmp_bm(ipts,j,ilabile,initrogen) - labile_target(ipts,j,initrogen)

             !+++CHECK+++
             ! CN-CAN proposed a simplified approach that seems more consistent.
             ! Rather than recalculating the carbres pool the calculation starts
             ! from the the carbon that is in carbres pool.
!!$          IF(is_tree(j))THEN
!!$             reserve_target = 0.12 * ( tmp_bm(ipts,j,isapabove,icarbon) + &
!!$                  tmp_bm(ipts,j,isapbelow,icarbon))/cn_leaf(ipts,j) * & 
!!$                  (1.+fcn_root(j)/ltor(ipts,j))/(1.+0.3/ltor(ipts,j))
!!$          ELSE
!!$             reserve_target = 0.3 * ( tmp_bm(ipts,j,iroot,icarbon) + &
!!$                  tmp_bm(ipts,j,isapabove,icarbon) + &
!!$                  tmp_bm(ipts,j,isapbelow,icarbon))/cn_leaf(ipts,j) * & 
!!$                  (1.+fcn_root(j)/ltor(ipts,j))/(1.+0.3/ltor(ipts,j))
!!$          ENDIF
             
             !+++CHECK+++
             ! Should we calculate the target for reserve nitrogen based on the
             ! target reserve carbon or the actual reserve carbon?
!!$             reserve_target(ipts,j,initrogen) = tmp_bm(ipts,j,icarbres,icarbon)/cn_leaf(ipts,j) * &
!!$                  (1.+fcn_root(j)/ltor(ipts,j))/(1.+root_reserve(j)/ltor(ipts,j))
             reserve_target(ipts,j,initrogen) = reserve_target(ipts,j,icarbon)/cn_leaf(ipts,j) * &
                  (1.+fcn_root(j)/ltor(ipts,j))/(1.+root_reserve(j)/ltor(ipts,j))
             !+++++++++++

             ! It needs to be avoided that the labile N can increase during the 
             ! dormancy (which includes the period that the buds are available 
             ! but still closed). If labile N increases, nitrogen uptake from 
             ! plant can be zero followed by unrealistic peak because it is the 
             ! function of labile nitrogen and carbon. Since this issue was fixed,
             ! the code does no longer allow the labile N to exceed the reserve 
             ! nitrogen during dormancy. This avoids earlier on in the code that 
             ! the labile pool can increase. This code is no longer needed and 
             ! was therefore commented out but it was left as a reminder of the issue.
!             IF(plant_status(ipts,j) .EQ. idormant .OR. &
!                plant_status(ipts,j) .EQ. ibudsavail) THEN
!                     use_max=-use_lab
!             ENDIF

             ! Sudden growth in diameter can occur if all of the following conditions
             ! are satisfied: 1) leaf C:N ratio and little allocation to root decreased
             ! the estimated N targets for labileN and reserveN; 2) high reserveN; 
             ! 3) then both labileN and reserveN were above the target; 4) no N 
             ! movement between pools; 5) extra N accumulates in labile; 6) n_alloc_tot
             ! increases; and 7) a lot of growth. To solve the problem the following 
             ! approach was implemented: 1) IF labileN is over the target, move extra 
             ! labileN to reserve.  Because N Growth is a lot dependent on the labile 
             ! pool in the current version so it is better to be cleaned. At a single 
             ! test site, remaining labile N as just much as the target pool resulted 
             ! in a slower initial increase in GPP but differences in endpoints were 
             ! negligible  2) IF labileN is under the target, move from reserve pool as
             ! much as extra N the reserve pool has. Move N from
             ! reserve pool much as needed but less than 90% of total reserve nitrogen. 
             ! Note that as there is a weak scientific basis for the movement of
             ! non-structural nitrogen, this redistribution of N is a numerical
             ! solution that aims at stabilizing the model.
             IF(use_lab .GT. 0.0) THEN 
                ! Enough N in labile: move N from labile to reserve
                use_max = -use_lab
             ELSE   
                ! depleted N in labile
                ! Fill the labile N as much as needed unless it acceeds 90% of
                ! reserveN. Before r6825 there is no N mobility when either
                ! labile and nitrogen are depleted. This change was made to use
                ! reserve N and reduce N limitation when the plant has nitrogen
                ! in its reserves
                use_max = MIN(-use_lab,0.9*tmp_bm(ipts,j,icarbres,initrogen))
             ENDIF

             tmp_bm(ipts,j,icarbres,initrogen) = tmp_bm(ipts,j,icarbres,initrogen)-use_max
             tmp_bm(ipts,j,ilabile,initrogen) = tmp_bm(ipts,j,ilabile,initrogen)+use_max
             bm_alloc(ipts,j,icarbres,initrogen) = bm_alloc(ipts,j,icarbres,initrogen)-use_max

             ! We need to keep the reserve nitrogen pool filled according to
             ! the filling status of the carbon reserve to allow decidiuous PFTs 
             ! to grow when we initialize an impose_cn=n run with the restarts 
             ! from an impose_cn=y simulation; if we don't ensure that there is 
             ! enough reserve N the PFT will get extinct.
             IF(is_tree(j)) THEN
                IF (pheno_type(j).NE.1) THEN
                   IF (impose_cn) THEN

                      ! take what is needed to keep nitrogen reserves optimal
                      atm_to_bm(ipts,j,initrogen) = atm_to_bm(ipts,j,initrogen) + &
                            MAX((tmp_bm(ipts,j,icarbres,icarbon)/cn_leaf(ipts,j) * & 
                            (1.+fcn_root(j)/ltor(ipts,j)) / &
                            (1.+root_reserve(j)/ltor(ipts,j))) - &
                            tmp_bm(ipts,j,icarbres,initrogen),zero)

                      ! fill the pool
                      tmp_bm(ipts,j,icarbres,initrogen) = &
                           MAX(tmp_bm(ipts,j,icarbres,initrogen), &
                           tmp_bm(ipts,j,icarbres,icarbon)/cn_leaf(ipts,j) * & 
                           (1.+fcn_root(j)/ltor(ipts,j))/(1.+root_reserve(j)/ltor(ipts,j)))

                      ! tmp_bm(icarbres) has changed, circ_class_biomass needs to be updated 
                      circ_class_biomass(ipts,j,:,icarbres,initrogen) = &
                           biomass_to_cc(tmp_bm(ipts,j,icarbres,initrogen),&
                           circ_class_biomass(ipts,j,:,icarbres,initrogen),&
                           circ_class_n(ipts,j,:))

                   ENDIF !impose_cn
                ENDIF
             ENDIF !is_tree
            
          ENDIF !IF_vegetmax

          ! Calculate demand or excess of nitrogen in the reserve
          ! reserve_target is the amount of nitrogen that is being targeted,
          ! icrabes is the actual amount of nitrogen. The balance is the
          ! ratio between the actual and the target and its long-term mean is
          ! used to modulate n-uptake by the roots.
          IF (reserve_target(ipts,j,initrogen) .GT. min_stomate) THEN
             n_reserve_balance(ipts,j) = tmp_bm(ipts,j,icarbres,initrogen) / &
                  reserve_target(ipts,j,initrogen) 
          ELSE
             ! If the reserver pool is zero the model might just be at an
             ! early time step. In that case n_reserve_balance should 
             ! probably be neutral (thus 1). If this happens later in the
             ! simulation, the model might be really short of nitrogen so 
             ! stimulating root uptake (n_reserve_balance << 1) is probably
             ! what is needed. Unless reserve_target suddenly dropped to zero,
             ! n_reserve_longeterm should be adjusted and is likley a good
             ! estimate for n_reserve_balance.
             n_reserve_balance(ipts,j) = n_reserve_longterm(ipts,j)
             ! As an alternative: compromise between the two options described 
             ! above but chose in favor of n-uptake
             !n_reserve_balance(ipts,j) = 0.5
          ENDIF

          ! tmp_bm can be very low. Truncate the distribution of n_reserve_balance
          n_reserve_balance(ipts,j) = MAX(zero, MIN(2.,n_reserve_balance(ipts,j)))

          ! biomass has changed so update circ_class_biomass
          DO ipar = 1,nparts
             circ_class_biomass(ipts,j,:,ipar,icarbon) = &
                  biomass_to_cc(tmp_bm(ipts,j,ipar,icarbon),&
                  circ_class_biomass(ipts,j,:,ipar,icarbon),&
                  circ_class_n(ipts,j,:))
             circ_class_biomass(ipts,j,:,ipar,initrogen) = &
                  biomass_to_cc(tmp_bm(ipts,j,ipar,initrogen),&
                  circ_class_biomass(ipts,j,:,ipar,initrogen),&
                  circ_class_n(ipts,j,:))
          ENDDO
       ENDDO ! npts

    ENDDO ! PFTs


 !! 9. Check mass balance closure

   IF (err_act.GT.1) THEN

      ! 9.2 Check surface area
      CALL check_vegetation_area("stomate_allocation", npts, veget_max_begin, &
           veget_max,'pft')

      ! 9.3 Mass balance closure 
      ! 9.3.1 Calculate final biomass
      pool_end(:,:,:) = zero
      DO ipar = 1,nparts
         DO iele = 1,nelements
            DO icir = 1,ncirc
               pool_end(:,:,iele) = pool_end(:,:,iele) + &
                    (circ_class_biomass(:,:,icir,ipar,iele) * &
                    circ_class_n(:,:,icir) * veget_max(:,:))
            ENDDO
         ENDDO
      ENDDO

      ! 9.3.2 Calculate mass balance 
      ! Specific processes for carbon
      check_intern(:,:,:,:) = zero
      check_intern(:,:,iatm2land,icarbon) = &
           check_intern_init(:,:,iatm2land,icarbon) + &
           (gpp_daily(:,:) + atm_to_bm(:,:,icarbon)) * dt * veget_max(:,:)
      check_intern(:,:,iatm2land,initrogen) = &
           check_intern_init(:,:,iatm2land,initrogen) + & 
           atm_to_bm(:,:,initrogen) * dt * veget_max(:,:)
      check_intern(:,:,iland2atm,icarbon) = -un * &
           (resp_maint(:,:) + resp_growth(:,:)) * &
           veget_max(:,:)

      ! Common processes for icarbon and initrogen
      DO iele=1,nelements
         check_intern(:,:,ipoolchange,iele) = -un * (pool_end(:,:,iele) - &
              pool_start(:,:,iele))
      ENDDO

      closure_intern = zero

      DO imbc = 1,nmbcomp
         DO iele=1,nelements
            ! Debug
            DO j=1,nvm
               IF (printlev_loc>=4 .AND. j.EQ.test_pft) THEN
                  WRITE(numout,*) 'check_intern, imbc, iele, ', imbc, &
                       iele, check_intern(:,test_pft,imbc,iele)
               END IF
            END DO
            !-
            closure_intern(:,:,iele) = closure_intern(:,:,iele) + &
                 check_intern(:,:,imbc,iele)
         END DO
      END DO

      ! 9.3.3 Check mass balance closure
      CALL check_mass_balance("stomate_allocation", closure_intern, npts, pool_end, &
           pool_start, veget_max, 'pft')

   END IF ! err_act.GT.1


   !! 10. Update leaf age
   !  Leaf age is needed to calculate the turnover and vmax in the 
   !  stomate_turnover.f90 and stomate_vmax.f90 routines. Leaf biomass 
   !  is distributed according to its age into several "age classes" 
   !  with age class=1 representing the youngest class, and consisting 
   !  of the most newly allocated leaf biomass. 
   
   !  Update biomass first
   tmp_bm(:,:,:,:) = cc_to_biomass(npts,nvm,circ_class_biomass(:,:,:,:,:),&
        circ_class_n(:,:,:))
   
   !! 9.1 Update quantity and age of the leaf biomass in the youngest class
   !  The new amount of leaf biomass in the youngest age class 
   !  (leaf_mass_young) is the sum of :
   !  - the leaf biomass that was already in the youngest age class 
   !  (leaf_frac(:,j,1) * lm_old(:,j)) with the leaf age given in 
   !  leaf_age(:,j,1) 
   !  - and the new biomass allocated to leaves 
   !  (bm_alloc(:,j,ileaf,icarbon)) with a leaf age of zero.
   leaf_mass_young(:,:) = leaf_frac(:,:,1) * lm_old(:,:) + bm_alloc(:,:,ileaf,icarbon)

   ! The age of the updated youngest age class is the average of the ages of 
   ! its 2 components: bm_alloc(leaf) of age '0', and leaf_frac * &
   ! lm_old(=leaf_mass_young-bm_alloc) of age 'leaf_age(:,:,1)' 
   DO ipts=1,npts

      DO j=1,nvm

         ! IF(veget_max(ipts,j) == zero)THEN
         !     ! this vegetation type is not present, so no reason to do the 
         !     ! calculation
         !     CYCLE
         ! ENDIF

         IF( (bm_alloc(ipts,j,ileaf,icarbon) .GT. min_stomate ) .AND. &
              ( leaf_mass_young(ipts,j) .GT. min_stomate ) )THEN

            leaf_age(ipts,j,1) = MAX ( zero, leaf_age(ipts,j,1) * &
                 ( leaf_mass_young(ipts,j) - bm_alloc(ipts,j,ileaf,icarbon) ) / &
                 & leaf_mass_young(ipts,j) )

         ENDIF

      ENDDO

   ENDDO

   !! 11 Update leaf age
   !  Update fractions of leaf biomass in each age class (fraction 
   !  in youngest class increases)

   !! 11.1 Update age of youngest leaves
   !  For age class 1 (youngest class), because we have added biomass 
   !  to the youngest class, we need to update the fraction of total 
   !  leaf biomass that belongs to the youngest age class : updated mass 
   !  in class divided by new total leaf mass
   WHERE ( tmp_bm(:,:,ileaf,icarbon) .GT. min_stomate )

      leaf_frac(:,:,1) = leaf_mass_young(:,:) / tmp_bm(:,:,ileaf,icarbon)

   ENDWHERE


   !! 11.2 Update age of other age classes
   !  Because the total leaf biomass has changed, we need to update the 
   !  fraction of leaves in each age class: mass in leaf age class (from 
   !  previous fraction of leaves in this class and previous total leaf 
   !  biomass) divided by new total mass
   DO m = 2, nleafages ! Loop over # leaf age classes

      WHERE ( tmp_bm(:,:,ileaf,icarbon) .GT. min_stomate )

         leaf_frac(:,:,m) = leaf_frac(:,:,m) * lm_old(:,:) / tmp_bm(:,:,ileaf,icarbon)

      ENDWHERE

   ENDDO       ! Loop over # leaf age classes
   !-----


   !! 12. Update whole-plant age 
   !! 12.1 PFT age
   !  At every time step, increase age of the biomass that was already 
   !  present at previous time step. Age is expressed in years, and 
   !  the time step 'dt' in days so age increase is: dt divided by number 
   !  of days in a year.
   WHERE ( PFTpresent(:,:) )

      age(:,:) = age(:,:) + dt/one_year

   ELSEWHERE

      age(:,:) = zero

   ENDWHERE


   !! 12.2 Age of grasses and crops
   !  For grasses and crops, biomass with age 0 has been added to the 
   !  whole plant with age 'age'. New biomass is the sum of the current 
   !  total biomass in all plant parts (bm_new), bm_new(:) = 
   !  SUM( tmp_bm(:,j,:), DIM=2 ). The biomass that has just been added 
   !  is the sum of the allocatable biomass of all plant parts (bm_add), 
   !  its age is zero. bm_add(:) = SUM( bm_alloc(:,j,:,icarbon), DIM=2 ). 
   !  Before allocation, the plant biomass is bm_new-bm_add, its age is 
   !  "age(:,j)". The age of the new biomass is the average of the ages 
   !  of previous and added biomass. For trees, age is treated in 
   !  "establish" if vegetation is dynamic, and in turnover routines if 
   !  it is static (in this case, only the age of the heartwood is 
   !  accounted for).
   DO j = 2,nvm

      IF ( .NOT. is_tree(j) ) THEN

         bm_new(:) = tmp_bm(:,j,ileaf,icarbon) + tmp_bm(:,j,isapabove,icarbon) + &
              tmp_bm(:,j,iroot,icarbon) + tmp_bm(:,j,ifruit,icarbon)
         bm_add(:) = bm_alloc(:,j,ileaf,icarbon) + bm_alloc(:,j,isapabove,icarbon) + &
              bm_alloc(:,j,iroot,icarbon) + bm_alloc(:,j,ifruit,icarbon)

         WHERE ( ( bm_new(:) .GT. min_stomate ) .AND. ( bm_add(:) .GT. min_stomate ) )

            age(:,j) = age(:,j) * ( bm_new(:) - bm_add(:) ) / bm_new(:)

         ENDWHERE

      ENDIF ! is .NOT. tree

   ENDDO  ! Loop over #PFTs 

   !! 13. Write history files

   ! Save in history file the variables describing the biomass allocated to the plant parts
   ! Calculate the change in biomass at the end of the routine
   tmp_bm(:,:,:,:) = cc_to_biomass(npts,nvm,circ_class_biomass(:,:,:,:,:),&
         circ_class_n(:,:,:)) - tmp_init_bm(:,:,:,:)

    DO l=1,nelements
       IF     (l == icarbon) THEN
          element_str(l) = '_c'
       ELSEIF (l == initrogen) THEN
          element_str(l) = '_n'
       ELSE
          STOP 'Define element_str'
       ENDIF

       CALL histwrite_p (hist_id_stomate, 'BM_ALLOC_LEAF'//TRIM(element_str(l)), itime, &
            tmp_bm(:,:,ileaf,l), npts*nvm, horipft_index)
       CALL histwrite_p (hist_id_stomate, 'BM_ALLOC_SAP_AB'//TRIM(element_str(l)), itime, &
            tmp_bm(:,:,isapabove,l), npts*nvm, horipft_index)
       CALL histwrite_p (hist_id_stomate, 'BM_ALLOC_SAP_BE'//TRIM(element_str(l)), itime, &
            tmp_bm(:,:,isapbelow,l), npts*nvm, horipft_index)
       CALL histwrite_p (hist_id_stomate, 'BM_ALLOC_ROOT'//TRIM(element_str(l)), itime, &
            tmp_bm(:,:,iroot,l), npts*nvm, horipft_index)
       CALL histwrite_p (hist_id_stomate, 'BM_ALLOC_FRUIT'//TRIM(element_str(l)), itime, &
            tmp_bm(:,:,ifruit,l), npts*nvm, horipft_index)
       CALL histwrite_p (hist_id_stomate, 'BM_ALLOC_RES'//TRIM(element_str(l)), itime, &
            tmp_bm(:,:,icarbres,l), npts*nvm, horipft_index)
       CALL histwrite_p (hist_id_stomate, 'BM_ALLOC_LABILE'//TRIM(element_str(l)), itime, &
            tmp_bm(:,:,ilabile,l), npts*nvm, horipft_index)
      
        CALL xios_orchidee_send_field('BM_ALLOC_LEAF'//TRIM(element_str(l)), tmp_bm(:,:,ileaf,l))
        CALL xios_orchidee_send_field('BM_ALLOC_SAP_AB'//TRIM(element_str(l)), tmp_bm(:,:,isapabove,l))
        CALL xios_orchidee_send_field('BM_ALLOC_SAP_BE'//TRIM(element_str(l)), tmp_bm(:,:,isapbelow,l))
        CALL xios_orchidee_send_field('BM_ALLOC_ROOT'//TRIM(element_str(l)), tmp_bm(:,:,iroot,l))
        CALL xios_orchidee_send_field('BM_ALLOC_FRUIT'//TRIM(element_str(l)), tmp_bm(:,:,ifruit,l))
        CALL xios_orchidee_send_field('BM_ALLOC_RES'//TRIM(element_str(l)), tmp_bm(:,:,icarbres,l))
        CALL xios_orchidee_send_field('BM_ALLOC_LABILE'//TRIM(element_str(l)), tmp_bm(:,:,ilabile,l))       

    ENDDO

    CALL histwrite_p (hist_id_stomate, 'RUE_LONGTERM', itime, &
         rue_longterm(:,:), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'KF', itime, &
         KF(:,:), npts*nvm, horipft_index)
    
    CALL xios_orchidee_send_field('RUE_LONGTERM', rue_longterm(:,:))
    CALL xios_orchidee_send_field('KF', KF(:,:))
    CALL xios_orchidee_send_field('REL_HEIGHT', height_rel(:,:)) 
    CALL xios_orchidee_send_field('C0_ALLOC', c0_alloc(:,:))
    CALL xios_orchidee_send_field('RESERVE_TARGET_c', reserve_target(:,:,icarbon))
    CALL xios_orchidee_send_field('LABILE_TARGET_c', labile_target(:,:,icarbon))
    CALL xios_orchidee_send_field('RESERVE_TARGET_n', reserve_target(:,:,initrogen))
    CALL xios_orchidee_send_field('LABILE_TARGET_n', labile_target(:,:,initrogen))
    CALL xios_orchidee_send_field('GTEMP_ALLOC', gtemp(:,:))
    CALL xios_orchidee_send_field('N_RESERVE_BALANCE',n_reserve_balance(:,:))
    CALL xios_orchidee_send_field('SLA', sla_out)

    ! Initilaize
    ring_width(:,:,:) = zero
    circ_height(:,:,:) = zero
    
    DO ipts = 1,npts
       DO j = 1,nvm
          IF(is_tree(j))THEN
  
             ! Calculate the forestry basal area (thus NOT the effective ba)
             circ_class_ba(:) = wood_to_ba(circ_class_biomass(ipts,j,:,:,icarbon),&
                  j, pipe_tune2(ipts,j))
             ba(ipts,j) = SUM(circ_class_ba(:)*circ_class_n(ipts,j,:)) * m2_to_ha
             store_circ_class_ba(ipts,j,:) = circ_class_ba(:)
             store_delta_ba(ipts,j,:) = circ_class_ba(:) - circ_class_ba_init(ipts,j,:)
             
             ! Calculate radius increment using store_delta_ba (which is alway
             ! positive) 
             ring_width(ipts,j,:) = SQRT(circ_class_ba(:)/pi) - &
                  SQRT(circ_class_ba_init(ipts,j,:)/pi)

             ! Calculate height per diameter class (above ground, not eff)
             circ_height(ipts,j,:) = wood_to_height(circ_class_biomass(ipts,j,:,:,icarbon),&
                  j, pipe_tune2(ipts,j)) 

          ELSE
             ba(ipts,j) = undef
             store_circ_class_ba(ipts,j,:) = undef
          ENDIF
       ENDDO
    ENDDO


    CALL histwrite_p (hist_id_stomate, 'BA', itime, &
         ba(:,:), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'RESIDUAL', itime, &
         residual_write(:,:), npts*nvm, horipft_index)

    ! Some of these variables should only be used when the model
    ! does not account for land cover changes and does not make
    ! use of age classes. When using land cover changes and
    ! age classes the biomass pool may get diluted which results
    ! in shrinking trees. This doesn't make sense at the tree
    ! level but is an acceptable approach at the landscape
    ! level.  
    DO icir = 1,ncirc
       WRITE(var_name,'(A,I3.3)') 'CCBA_',icir
       CALL histwrite_p (hist_id_stomate, var_name, itime, &
            store_circ_class_ba(:,:,icir), npts*nvm, horipft_index)
       WRITE(var_name,'(A,I3.3)') 'CCDELTABA_',icir
       CALL histwrite_p (hist_id_stomate, VAR_NAME, itime, &
            store_delta_ba(:,:,icir), npts*nvm, horipft_index)
       WRITE(var_name,'(A,I3.3)') 'CCN_',icir
       CALL histwrite_p (hist_id_stomate, VAR_NAME, itime, &
            circ_class_n(:,:,icir), npts*nvm, horipft_index)
       WRITE(var_name,'(A,I3.3)') 'CCTRW_',icir
       CALL histwrite_p (hist_id_stomate, VAR_NAME, itime, &
            ring_width(:,:,icir), npts*nvm, horipft_index)
       WRITE(var_name,'(A,I3.3)') 'CCH_',icir
       CALL histwrite_p (hist_id_stomate, VAR_NAME, itime, &
            circ_height(:,:,icir), npts*nvm, horipft_index)
    ENDDO

    CALL xios_orchidee_send_field('BA',ba(:,:))
    CALL xios_orchidee_send_field('RESIDUAL',residual_write(:,:))
    CALL xios_orchidee_send_field('CCBA', store_circ_class_ba(:,:,:))
    CALL xios_orchidee_send_field('CCDELTABA', store_delta_ba(:,:,:)) 
    CALL xios_orchidee_send_field('CCIND', circ_class_n(:,:,:)) 
    CALL xios_orchidee_send_field('CCTRW', ring_width(:,:,:))
    CALL xios_orchidee_send_field('CCHEIGHT', circ_height(:,:,:))
    WHERE (store_circ_class_ba(:,:,:).EQ.undef)
        histvar2(:,:,:) = xios_default_val
    ELSE WHERE
        histvar2(:,:,:) = SQRT(4/pi*store_circ_class_ba(:,:,:))
    END WHERE
    CALL xios_orchidee_send_field('CCDIAMETER', histvar2(:,:,:))
    CALL xios_orchidee_send_field('CCSAP_M_AB_c', circ_class_biomass(:,:,:,isapabove,icarbon))
    CALL xios_orchidee_send_field('CCHEART_M_AB_c',circ_class_biomass(:,:,:,iheartabove,icarbon))
    CALL xios_orchidee_send_field('CCSAP_M_BE_c', circ_class_biomass(:,:,:,isapbelow,icarbon))
    CALL xios_orchidee_send_field('CCTOTAL_M_c', SUM(circ_class_biomass(:,:,:,:,icarbon),DIM=4))

    ! Send value that caused a warning to xios
    CALL xios_orchidee_send_field('MBC_alloc10b_c', residual10b(:,:))

    ! +++ DBEN Output +++
    CALL xios_orchidee_send_field('BA_dben',ba(:,:))
    CALL xios_orchidee_send_field('BAgrowth_dben', &
         SUM(store_delta_ba(:,:,:)*circ_class_n(:,:,:),3)*m2_to_ha*365)
    ! +++++++++++++++++++

    ! Debug
    IF (printlev_loc.GE.4) THEN
       WRITE(numout,*) 'leaf_biomass C:', tmp_bm(test_grid,test_pft,ileaf,icarbon)
       WRITE(numout,*) 'leaf_biomass N:', tmp_bm(test_grid,test_pft,ileaf,initrogen)
       WRITE(numout,*) 'reserve_biomass C:', tmp_bm(test_grid,test_pft,icarbres,icarbon)
       WRITE(numout,*) 'reserve_biomass N:', tmp_bm(test_grid,test_pft,icarbres,initrogen)
       WRITE(numout,*) 'circ_height(test_grid,test_pft,:) =', circ_height(test_grid,test_pft,:) 
    ENDIF 
    !-

    IF (printlev.GE.3) WRITE(numout,*) 'Leaving functional allocation growth'

    
END SUBROUTINE growth_fun_all



!! ================================================================================================================================
!! FUNCTION	: func_derfunc
!!
!>\BRIEF        Calculate value for a function and its derivative
!!
!!
!! DESCRIPTION  : the routine describes the function and its derivative. Both function and derivative are used
!!              by the optimisation scheme. Hence, this function is part of the optimisation scheme and is only 
!!              called by the optimisation 
!!
!! RECENT CHANGE(S): 
!!
!! MAIN OUTPUT VARIABLE(S): f, df
!!
!! REFERENCE(S)	: Numerical recipies in Fortran 77
!! 
!! FLOWCHART : 
!! \n
!_ ================================================================================================================================
 
 SUBROUTINE func_derfunc(x, n, o, p, q, r, t, eq_num, f, df)

!! 0. Variable and parameter declaration

    !! 0.1 Input variables
    REAL(r_std), INTENT(in)                :: x           !! x value for which the function f(x) will be evaluated
    REAL(r_std), INTENT(in)                :: n,o,p,q,r,t !! Coefficients of the equation. Not all equations use all coefficients
    INTEGER(i_std), INTENT(in)             :: eq_num      !! Function i.e. f(x), g(x), ...

    !! 0.2 Output variables
    REAL(r_std), INTENT(out)               :: f           !! Value y for f(x)
    REAL(r_std), INTENT(out)               :: df          !! Value y for derivative[f(x)]   

    !! 0.3 Modified variables

    !! 0.4 Local variables
!_ ================================================================================================================================

!! 1. Calculate f(x) and df(x)

    IF (eq_num .EQ. 1) THEN

       !f = n*x**4 + o*x**3 + p*x**2 + q*x + r  
       !df = 4*n*x**3 + 3*o*x**2 + 2*p*x + q
    
    ELSEIF (eq_num .EQ. 2) THEN
    
       f = ( (n*x)/(p*((x+o)/t)**(q/(2+q))) ) - r
       df = ( n*(o*(q+2)+2*x)*((o+x)/t)**(-q/(q+2)) ) / ( p*(q+2)*(o+x) )
    
    ENDIF

 END SUBROUTINE func_derfunc


!! ================================================================================================================================
!! FUNCTION	: iterative_solver
!!
!>\BRIEF        find best fitting x for f(x)
!!
!!
!! DESCRIPTION  : The function makes use of an iterative approach to optimise the value for X. The solver
!!              splits the search region in two but there is an additional check to ensure that bounds are not
!!              exceeded. 
!!
!!              Use the derivative (df) of the function (f) calculated in func_derfunc to narrow down the 
!!              search range for X using the Newton-Raphson method for convergence as described in Numerical
!!              Recipes in Fortran 77 (page 355-360).            
!!
!! RECENT CHANGE(S): 
!!
!! MAIN OUTPUT VARIABLE(S): x
!!
!! REFERENCE(S)	: Numerical recipies in Fortran 77
!! 
!! FLOWCHART : 
!! \n
!_ ================================================================================================================================

  FUNCTION newX(n, o, p, q, r, s, t, x1, x2, eq_num, j, ipts)

!! 0. Variable and parameter declaration

    !! 0.1 Input variables
    REAL(r_std), INTENT(in)        :: n,o,p,q,r,s,t    !! Coefficients of the equation. Not all 
                                                       !! equations use all coefficients
    REAL(r_std), INTENT(in)        :: x1               !! Lower boundary off search range
    REAL(r_std), INTENT(in)        :: x2               !! Upper boundary off search range
    INTEGER(i_std), INTENT(in)     :: eq_num           !! Function for which an iterative solution is 
                                                       !! searched
    INTEGER(i_std), INTENT(in)     :: j                !! Number of PFT
    INTEGER(i_std), INTENT(in)     :: ipts             !! Number of grdi square...for debugging
    
    !! 0.2 Output variables
   
    !! 0.3 Modified variables

    !! 0.4 Local variables
    INTEGER(i_std), PARAMETER      :: maxit = 20       !! Maximum number of iterations
    INTEGER(i_std), PARAMETER      :: max_attempt = 5  !! Maximum number of iterations
    INTEGER(i_std)                 :: i, attempt       !! Index
    REAL(r_std)                    :: newX             !! New estimate for X
    REAL(r_std)                    :: fl, fh, f        !! Value of the function for the lower bound (x1), 
                                                       !! upper bound (x2) and the new value (newX)
    REAL(r_std)                    :: xh, xl           !! Checked lower and upper bounds
    REAL(r_std)                    :: df               !! Value of the derivative of the function for newX
    REAL(r_std)                    :: dx, dxold        !! Slope of improvement
    REAL(r_std)                    :: temp             !! Dummy variable for value swaps
    REAL(r_std)                    :: low, high        !! temporary variables for x1 and x2 to avoid
                                                       !! intent in/out conflicts with Cs
    LOGICAL                        :: found_range      !! Flag indicating whether the range in which
                                                       !! a solution exists was identified.
    
    
!_ ================================================================================================================================   

!! 1. Find solution for X
 
    ! Not sure whether our initial range is large enough. We will
    ! start with a narrow range so we are more likely to find the
    ! solution witin ::maxit iterations. If there is no solution
    ! in the initial range we will expand the range and try again

    ! Initilaze flags and counters
    attempt = 1
    found_range = .FALSE.
    low = x1
    high = x2

    ! Calculate y for the upper and lower bound
    DO WHILE (.NOT. found_range .AND. attempt .LE. max_attempt)
 
       CALL func_derfunc(low, n, o, p, q, r, t, eq_num, fl, df)
       CALL func_derfunc(high, n, o, p, q, r, t, eq_num, fh, df)
 
       IF ((fl .GT. 0.0 .AND. fh .GT. 0.0) .OR. &
            (fl .LT. 0.0 .AND. fh .LT. 0.0)) THEN
       
          ! Update counter
          attempt = attempt + 1

          IF (attempt .GT. max_attempt) THEN

             ! If the sign of y does not changes between the upper 
             ! and lower bound there no solution with the specified range
             found_range = .FALSE.

          ELSE
             
             IF (fl.GE.zero) THEN

                ! Both values are positive
                IF (fh.GT.fl) THEN

                   ! the only option to get a difference in sign
                   ! is by decreasing the lower boundary such that
                   ! fl can become negative. Both fh and fl are too
                   ! large so the current lower boundary could be
                   ! used as the upper boundary in the next attampt.
                   temp = low
                   low = low / 2
                   high = temp

                ELSE

                   ! Strange. We should decrease the higher bound
                   ! to find a solution. the uper bound will thus approach
                   ! the lower bound. The solution should alreday be in the
                   ! initial range. Unless we are searching for a local
                   ! minimum in the range. This is possible when the initial
                   ! range is too large. Try swapping the values and hope
                   ! that the next attempt will be more logic.
                   temp = low
                   low = high
                   high = temp
                   
                ENDIF
                
             ELSEIF (fl.LT.zero) THEN

                IF (fh.GT.fl) THEN

                   ! the only option to get a difference in sign
                   ! is by increasing the upper boundary such that
                   ! fh can become psitive. Both boundaries were
                   ! too low so it is safe to give the lower boundary
                   ! the value of the higher boundary and to
                   ! increase the higher boundary. Increase the
                   ! higher boundary by a factor of 10.
                   temp = high
                   high = high * 10
                   low = temp

                ELSE

                   ! Srange. Swap the values. See above for more
                   ! details on the reasoning.
                   temp = low
                   low = high
                   high = temp

                ENDIF

             ELSE

                ! Overlooked something
                CALL ipslerr_p(3,'logical flaw in an IF-statement', &
                     'case 1 in newX', '','')
             ENDIF ! fl.GE.zero
             
             ! Enlarge the search range
             IF(printlev_loc.GE.4)THEN
                WRITE(numout,*) 'Iterative procedure - enlarge the search range'
                WRITE(numout,*) 'New range: ', low, high
                WRITE(numout,*) 'PFT, grid square, range: ',j,ipts,attempt
             ENDIF

          ENDIF ! attempt.GT.max_attempt
       
       ELSE

          found_range = .TRUE.
          
       ENDIF ! fl and fh have the same sign
      
    ENDDO ! .NOT. found_range .AND. attempt .LE. max_attempt

    ! Only when we found a range we will search for the solution 
    IF (found_range) THEN

       ! If the sign of y changes between the upper and lower bound there is a solution
       IF ( ABS(fl) .LT. min_stomate ) THEN          

          ! The lower bound is the solution
          newX = low
          RETURN

       ELSEIF ( ABS(fh) .LT. min_stomate ) THEN
          
          ! The upper bound is the solution
          newX = high
          RETURN

       ELSEIF (fl .LT. 0.0) THEN
          
          ! Accept the lower and upper bounds as specified
          xl = low
          xh = high
          
       ELSE

          ! Lower and upper bounds were swapped, correct their ranking 
          xh = low
          xl = high
          
       ENDIF

       ! Estimate the initial newX value 	
       newX = 0.5 * (low+high)
       dxold = ABS(high-low)
       dx = dxold

       ! Calculate y=f(x) and df(x) for initial guess of newX
       CALL func_derfunc(newX, n, o, p, q, r, t, eq_num, f, df)

       ! Evaluate for the maximum number of iterations  
       DO  i = 1,maxit
          
          IF ( ((newX-xh)*df-f)*((newX-xl)*df-f) .GT. 0.0 .OR. ABS(deux*f) > ABS(dxold*df) ) THEN
          
             ! Bisection
             dxold = dx
             dx = 0.5 * (xh-xl)
             newX = xl+dx
             IF (xl .EQ. newX) RETURN
             
          ELSE
             
             ! Newton
             dxold = dx
             dx = f/df
             temp = newX
             newX = newX-dx
             IF (temp .EQ. newX) RETURN
          
          ENDIF
       
          ! Precision reached
          IF ( ABS(dx) .LT. min_stomate) RETURN
          
          ! Precision was not reached calculate f(x) and df(x) for newX
          CALL func_derfunc(newX, n, o, p, q, r, t, eq_num, f, df)
          
          ! Narrow down the range
          IF (f .LT. 0.0) then
             xl = newX
          ELSE
             xh = newX
          ENDIF

       ENDDO ! maximum number of iterations

    ELSE

       ! The fact that we did not find a solution in the initial range or
       ! a slightly adjusted range suggests that the new value of X (= Cs)
       ! is very different from its current value. That is a bit strange;
       ! something crazy may have happened with Cs. No optimal solution was
       ! found but we will just adjust newX to our first guess of Cl_target.
       ! Note that if we get too many suboptimal solution in a row, the error
       ! will propagate in the code resulting in unrealistic biomasses. At one
       ! point LAIs of 5000 were observed. It is important that this suboptimal 
       ! solution is as realistic as possible.
       IF (fl.GT.zero .AND. fh.GT.fl) THEN
          newX = s
       ELSEIF (fl.GT.zero .AND. fh.LE.fl) THEN
          newX = s 
       ELSEIF (fl.LT.zero .AND. fh.LE.fl) THEN
          newX = s 
       ELSEIF (fl.LT.zero .AND. fh.GT.fl) THEN
          newX = s 
       ELSE
          ! Overlooked something
          CALL ipslerr_p(3,'logical flaw in an IF-statement', &
               'case 1 in newX', '','')
       ENDIF

       WRITE(numout,*) 'PFT, grid square: ',j,ipts
       CALL ipslerr_p (2,'growth_fun_all','newX',&
            'Iterative procedure - tried really hard but failed',&
            'had to use a suboptimal solution instead')

    ENDIF

  END FUNCTION newX


!! ================================================================================================================================
!! SUBROUTINE	: comments
!!
!>\BRIEF        Contains all comments to check the code
!!
!!
!! DESCRIPTION  : contains all comments to check the code. By setting pft_test to 0, this routine is not called 
!!
!! RECENT CHANGE(S): 
!!
!! MAIN OUTPUT VARIABLE(S): none
!!
!! REFERENCE(S)	: none
!! 
!! FLOWCHART : 
!! \n
!_ ================================================================================================================================

  SUBROUTINE comment(npts, Cl_target, Cl, Cs_target, & 
       Cs, Cr_target, Cr, delta_ba, &
       ipts, j, l, b_inc_tot, & 
       Cl_incp, Cs_incp, Cr_incp, KF, LF, &
       Cl_inc, Cs_inc, Cr_inc, Cf_inc, &
       grow_wood, circ_class_n, n_comment)

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                         :: npts		                    !! Defined in stomate_growth_fun_all
    REAL(r_std), DIMENSION(:), INTENT(in)              :: Cl_target, Cs_target, Cr_target   !! Defined in stomate_growth_fun_all
    REAL(r_std), DIMENSION(:), INTENT(in)              :: Cl_incp, Cs_incp, Cr_incp         !! Defined in stomate_growth_fun_all
    REAL(r_std), DIMENSION(:), INTENT(in)              :: Cl_inc, Cs_inc, Cr_inc, Cf_inc    !! Defined in stomate_growth_fun_all
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)          :: circ_class_n                      !! Defined in stomate_growth_fun_all
    REAL(r_std), DIMENSION(:), INTENT(in)              :: Cl, Cs, Cr                        !! Defined in stomate_growth_fun_all
    REAL(r_std), DIMENSION(:), INTENT(in)              :: delta_ba                          !! Defined in stomate_growth_fun_all
    REAL(r_std), DIMENSION(:,:), INTENT(in)            :: KF, LF                            !! Defined in stomate_growth_fun_all
    REAL(r_std), INTENT(in)                            :: b_inc_tot                         !! Defined in stomate_growth_fun_all
    INTEGER(i_std), INTENT(in)                         :: ipts, j, l                        !! Defined in stomate_growth_fun_all
    LOGICAL, INTENT(in)                                :: grow_wood                         !! Defined in stomate_growth_fun_all
                                                                                            !! Note that grow_wood is only calculated for n_comment=10 and higher

    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables
    INTEGER(i_std)                                     :: n_comment                         !! Comment number  
    !_ ================================================================================================================================

    SELECT CASE (n_comment)
    CASE (1)
       ! Enough leaves and wood, grow roots
       WRITE(numout,*) 'Exc 1: Cl_incp(=0), Cs_incp (=0), Cr_incp (<>0), unallocated, class, '
       WRITE(numout,*) Cl_incp(l), Cs_incp(l), Cr_incp(l), &
            b_inc_tot - (circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ), l
       IF (b_inc_tot - (circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ) .LT. -10*EPSILON(zero)) THEN
          WRITE(numout,*) 'Exc 1.1: unallocated less then 0: overspending, ', b_inc_tot - &
               (circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)))
       ELSE
          IF ( (b_inc_tot - ( circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ) .GE. -10*EPSILON(zero)) .AND. &
               ((circ_class_n(ipts,j,l) * ABS(Cr_target(l)-Cr(l)-Cr_incp(l))) .LE. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 1.2: unallocated <>= 0 but tree is in good shape: successful allocation'
          ELSEIF ( (b_inc_tot - ( circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ) &
               .LE. min_stomate) .AND. &
               (circ_class_n(ipts,j,l) * ABS(Cr_target(l)-Cr(l)-Cr_incp(l)) .GT. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 1.3: unallocated = 0, but the tree needs more reshaping: successful allocation'
          ELSE
             WRITE(numout,*) 'WARNING 24: Exc 1.4 unexpected result'
             WRITE(numout,*) 'WARNING 24: PFT, ipts: ',j,ipts
          ENDIF
       ENDIF

    CASE (2)
       ! Enough wood and roots, grow leaves 
       WRITE(numout,*) 'Exc 2: Cl_incp(<>0), Cs_incp (=0), Cr_incp (=0), unallocated, class, '
       WRITE(numout,*) Cl_incp(l), Cs_incp(l), Cr_incp(l), &
            b_inc_tot - (circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ), l
       IF (b_inc_tot - (circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ) .LT. -10*EPSILON(zero)) THEN
          WRITE(numout,*) 'Exc 2.1: unallocated less then 0: overspending, ', b_inc_tot - &
               (circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)))
       ELSE
          IF ( (b_inc_tot - ( circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ) .GE. -10*EPSILON(zero)) .AND. &
               ((circ_class_n(ipts,j,l) * ABS(Cl_target(l)-Cl(l)-Cl_incp(l))) .LE. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 2.2: unallocated <>= 0 but tree is in good shape: successful allocation'
          ELSEIF ( (b_inc_tot - ( circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ) &
               .LE. min_stomate) .AND. &
               ((circ_class_n(ipts,j,l) * ABS(Cl_target(l)-Cl(l)-Cl_incp(l))) .GT. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 2.3: unallocated = 0, but the tree needs more reshaping: successful allocation'
          ELSE
             WRITE(numout,*) 'WARNING 25: Exc 2.4 unexpected result'
             WRITE(numout,*) 'WARNING 25: PFT, ipts: ',j,ipts
          ENDIF
       ENDIF


    CASE (3)

       ! Enough wood, grow leaves and roots
       WRITE(numout,*) 'Exc 3: Cl_incp(<>0), Cs_incp(=0), Cr_incp(<>0), unallocated, class, '
       WRITE(numout,*) Cl_incp(l), Cs_incp(l), Cr_incp(l), b_inc_tot - & 
            (circ_class_n(ipts,j,l) * (Cl_incp(l)+Cs_incp(l)+Cr_incp(l))), l
       IF (b_inc_tot - circ_class_n(ipts,j,l) * (Cl_incp(l) + Cs_incp(l) + Cr_incp(l))  &
            .LT. -10*EPSILON(zero)) THEN
          WRITE(numout,*) 'Exc 3.1: unallocated less then 0: overspending, ', b_inc_tot - &
               (circ_class_n(ipts,j,l) * (Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) )
       ELSE
          IF ( (b_inc_tot - (circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l))) &
               .GE. min_stomate) .AND. &
               ((circ_class_n(ipts,j,l) * ABS(Cl_target(l)-Cl(l)-Cl_incp(l)) ) .LE. min_stomate) .AND. &
               ((circ_class_n(ipts,j,l) * ABS(Cr_target(l)-Cr(l)-Cr_incp(l)) ) .LE. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 3.2: unallocated <>= 0 but tree is in good shape: successful allocation'
          ELSEIF ( (b_inc_tot - (circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ) &
               .LE. min_stomate) .AND. &
               ((circ_class_n(ipts,j,l) * ABS(Cl_target(l)-Cl(l)-Cl_incp(l)) ) .GT. min_stomate) .AND. &
               ((circ_class_n(ipts,j,l) * ABS(Cr_target(l)-Cr(l)-Cr_incp(l)) ) .GT. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 3.3: unallocated = 0, but the tree needs more reshaping: successful allocation'
          ELSE
             WRITE(numout,*) 'WARNING 26: Exc 3.4 unexpected result'
             WRITE(numout,*) 'WARNING 26: PFT, ipts: ',j,ipts
          ENDIF
       ENDIF

    CASE(4)
       ! Enough leaves and wood, grow roots
       WRITE(numout,*) 'Exc 4: Cl_incp(=0), Cs_incp (=0), Cr_incp (<>0), unallocated, class, '
       WRITE(numout,*) Cl_incp(l), Cs_incp(l), Cr_incp(l), &
            b_inc_tot - (circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ), l
       IF (b_inc_tot - (circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ) .LT. -10*EPSILON(zero)) THEN
          WRITE(numout,*) 'Exc 4.1: unallocated less then 0: overspending, ', b_inc_tot - &
               (circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)))
       ELSE
          IF ( (b_inc_tot - ( circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ) .GE. -10*EPSILON(zero)) .AND. &
               ((circ_class_n(ipts,j,l) * ABS(Cr_target(l)-Cr(l)-Cr_incp(l))) .LE. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 4.2: unallocated <>= 0 but tree is in good shape: successful allocation'
          ELSEIF ( (b_inc_tot - ( circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ) &
               .LE. min_stomate) .AND. &
               ((circ_class_n(ipts,j,l) * ABS(Cr_target(l)-Cr(l)-Cr_incp(l))) .GT. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 4.3: unallocated = 0, but the tree needs more reshaping: successful allocation'
          ELSE
             WRITE(numout,*) 'WARNING 27: Exc 4.4 unexpected result'
             WRITE(numout,*) 'WARNING 27: PFT, ipts: ',j,ipts
          ENDIF
       ENDIF

    CASE(5)
       ! Enough leaves and roots, grow wood
       WRITE(numout,*) 'Exc 5: Cl_incp(=0), Cs_incp (<>0), Cr_incp (=0), unallocated, class, '
       WRITE(numout,*) Cl_incp(l), Cs_incp(l), Cr_incp(l), &
            b_inc_tot - (circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ), l
       IF (b_inc_tot - (circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ) .LT. -10*EPSILON(zero)) THEN
          WRITE(numout,*) 'Exc 5.1: unallocated less then 0: overspending, ', b_inc_tot - &
               (circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)))
       ELSE
          IF ( (b_inc_tot - ( circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ) .GE. -10*EPSILON(zero)) .AND. &
               ((circ_class_n(ipts,j,l) * ABS(Cs_target(l)-Cs(l)-Cs_incp(l))) .LE. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 5.2: unallocated <>= 0 but tree is in good shape: successful allocation'
          ELSEIF ( (b_inc_tot - ( circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ) &
               .LE. min_stomate) .AND. &
               ((circ_class_n(ipts,j,l) * ABS(Cs_target(l)-Cs(l)-Cs_incp(l))) .GT. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 5.3: unallocated = 0, but the tree needs more reshaping: successful allocation'
          ELSE
             WRITE(numout,*) 'WARNING 28: Exc 5.4 unexpected result'
             WRITE(numout,*) 'WARNING 28: PFT, ipts: ',j,ipts
          ENDIF
       ENDIF

    CASE(6)
       ! Enough leaves, grow wood and roots
       WRITE(numout,*) 'Exc 6: Cl_incp(=0), Cs_incp(<>0), Cr_incp(<>0), unallocated'
       WRITE(numout,*) Cl_incp(l), Cs_incp(l), Cr_incp(l), &
            b_inc_tot - (circ_class_n(ipts,j,l) * (Cl_incp(l)+Cs_incp(l)+Cr_incp(l)))
       IF (b_inc_tot - (circ_class_n(ipts,j,l) * (Cl_incp(l)+Cs_incp(l)+Cr_incp(l))) .LT. -10*EPSILON(zero)) THEN
          WRITE(numout,*) 'Exc 6.1: unallocated less then 0: overspending, ', &
               b_inc_tot - (circ_class_n(ipts,j,l) * (Cl_incp(l)+Cs_incp(l)+Cr_incp(l)))
       ELSE
          IF ( (b_inc_tot - (circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l))) .GE. -10*EPSILON(zero)) .AND. &
               ((circ_class_n(ipts,j,l) * ABS(Cs_target(l)-Cs(l)-Cs_incp(l))) .LE. min_stomate) .AND. &
               ((circ_class_n(ipts,j,l) * ABS(Cr_target(l)-Cr(l)-Cr_incp(l))) .LE. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 6.2: unallocated <>= 0 but tree is in good shape: successful allocation'
          ELSEIF ( (b_inc_tot - (circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l))) .LE. min_stomate) .AND. &
               ((circ_class_n(ipts,j,l) * ABS(Cs_target(l)-Cs(l)-Cs_incp(l))) .GT. min_stomate) .OR. &
               ((circ_class_n(ipts,j,l) * ABS(Cr_target(l)-Cr(l)-Cr_incp(l))) .GT. min_stomate) ) THEN
             WRITE(numout,*) &
                  'Exc 6.3: unallocated = 0, but the tree needs more reshaping: successful allocation'
          ELSE
             WRITE(numout,*) 'WARNING 29: Exc 6.4 unexpected result'
             WRITE(numout,*) 'WARNING 29: PFT, ipts: ',j,ipts
          ENDIF
       ENDIF

    CASE(7)
       ! Enough leaves and wood, grow roots
       WRITE(numout,*) 'Exc 7: Cl_incp(=0), Cs_incp (=0), Cr_incp (<>0), unallocated, class, '
       WRITE(numout,*) Cl_incp(l), Cs_incp(l), Cr_incp(l), &
            b_inc_tot - (circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ), l
       IF (b_inc_tot - (circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ) .LT. -10*EPSILON(zero)) THEN
          WRITE(numout,*) 'Exc 7.1: unallocated less then 0: overspending, ', b_inc_tot - &
               (circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)))
       ELSE
          IF ( (b_inc_tot - ( circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ) .GE. -10*EPSILON(zero)) .AND. &
               ((circ_class_n(ipts,j,l) * ABS(Cl_target(l)-Cl(l)-Cl_incp(l))) .LE. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 7.2: unallocated <>= 0 but tree is in good shape: successful allocation'
          ELSEIF ( (b_inc_tot - ( circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ) &
               .LE. min_stomate) .AND. &
               ((circ_class_n(ipts,j,l) * ABS(Cl_target(l)-Cl(l)-Cl_incp(l))) .GT. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 7.3: unallocated = 0, but the tree needs more reshaping: successful allocation'
          ELSE
             WRITE(numout,*) 'WARNING 30: Exc 7.4 unexpected result'
             WRITE(numout,*) 'WARNING 30: PFT, ipts: ',j,ipts
          ENDIF
       ENDIF

    CASE(8)
       ! Enough leaves and roots, grow wood
       WRITE(numout,*) 'Exc 8: Cl_incp(=0), Cs_incp (<>0), Cr_incp (=0), unallocated, class, '
       WRITE(numout,*) Cl_incp(l), Cs_incp(l), Cr_incp(l), &
            b_inc_tot - (circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ), l
       IF (b_inc_tot - (circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ) .LT. -10*EPSILON(zero)) THEN
          WRITE(numout,*) 'Exc 8.1: unallocated less then 0: overspending, ', b_inc_tot - &
               (circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)))
       ELSE
          IF ( (b_inc_tot - ( circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ) .GE. -10*EPSILON(zero)) .AND. &
               ((circ_class_n(ipts,j,l) * ABS(Cs_target(l)-Cs(l)-Cs_incp(l))) .LE. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 8.2: unallocated <>= 0 but tree is in good shape: successful allocation'
          ELSEIF ( (b_inc_tot - ( circ_class_n(ipts,j,l)*(Cl_incp(l)+Cs_incp(l)+Cr_incp(l)) ) &
               .LE. min_stomate) .AND. &
               ((circ_class_n(ipts,j,l) * ABS(Cs_target(l)-Cs(l)-Cs_incp(l))) .GT. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 8.3: unallocated = 0, but the tree needs more reshaping: successful allocation'
          ELSE
             WRITE(numout,*) 'WARNING 31: Exc 8.4 unexpected result'
             WRITE(numout,*) 'WARNING 31: PFT, ipts: ',j,ipts
          ENDIF
       ENDIF

    CASE(9)
       ! Enough roots, grow leaves and wood
       WRITE(numout,*) 'Exc 9: delta_ba, Cl_incp(<>0), Cs_incp(<>0), Cr_incp(=0), unallocated, class, '
       WRITE(numout,*) delta_ba(:), Cl_incp(l), Cs_incp(l), Cr_incp(l), &
            b_inc_tot - (circ_class_n(ipts,j,l) * (Cl_incp(l)+Cs_incp(l)+Cr_incp(l))), l
       IF (b_inc_tot - (circ_class_n(ipts,j,l) * (Cl_incp(l)+Cs_incp(l)+Cr_incp(l))) .LT. -10*EPSILON(zero)) THEN
          WRITE(numout,*) 'Exc 9.1: unallocated less then 0: overspending, ', &
               b_inc_tot - (circ_class_n(ipts,j,l) * (Cl_incp(l)+Cs_incp(l)+Cr_incp(l)))
       ELSE
          IF ( (b_inc_tot - (circ_class_n(ipts,j,l) * (Cl_incp(l)+Cs_incp(l)+Cr_incp(l))) .GE. -10*EPSILON(zero)) .AND. &
               ((circ_class_n(ipts,j,l) * ABS(Cl_target(l)-Cl(l)-Cl_incp(l))) .LE. min_stomate) .AND. &
               ((circ_class_n(ipts,j,l) * ABS(Cs_target(l)-Cs(l)-Cs_incp(l))) .LE. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 9.2: unallocated <>= 0 but tree is in good shape: successful allocation'
          ELSEIF ( (b_inc_tot - (circ_class_n(ipts,j,l) * (Cl_incp(l)+Cs_incp(l)+Cr_incp(l))) &
               .LE. min_stomate) .AND. &
               ((circ_class_n(ipts,j,l) * ABS(Cl_target(l)-Cl(l)-Cl_incp(l))) .GT. min_stomate) .OR. &
               ((circ_class_n(ipts,j,l) * ABS(Cs_target(l)-Cs(l)-Cs_incp(l))) .GT. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 9.3: unallocated = 0, but the tree needs more reshaping: successful allocation'
          ELSE
             WRITE(numout,*) 'WARNING 32: Exc 9.4 unexpected result'
             WRITE(numout,*) 'WARNING 32: PFT, ipts: ',j,ipts
          ENDIF
       ENDIF

    CASE(10)
       ! Ready for ordinary allocation
       WRITE(numout,*) 'Ready for ordinary allocation?'
       WRITE(numout,*) 'KF, LF, ', KF(ipts,j), LF(ipts,j)
       WRITE(numout,*) 'b_inc_tot, ', b_inc_tot
       WRITE(numout,*) 'Cl, Cs, Cr', Cl(:), Cs(:), Cr(:)
       WRITE(numout,*) 'Cl_target-Cl, ', Cl_target(:)-Cl(:)
       WRITE(numout,*) 'Cs_target-Cs, ', Cs_target(:)-Cs(:)
       WRITE(numout,*) 'Cr_target-Cr, ', Cr_target(:)-Cr(:)
       IF (b_inc_tot .GT. min_stomate) THEN
          IF (SUM(ABS(Cl_target(:)-Cl(:))) .LE. min_stomate) THEN
             IF (SUM(ABS(Cs_target(:)-Cs(:))) .LE. min_stomate) THEN
                IF (SUM(ABS(Cr_target(:)-Cr(:))) .LE. min_stomate) THEN
                   IF (grow_wood) THEN
                      WRITE(numout,*) 'should result in exc 10.1 or 10.2'
                   ELSE
                      WRITE(numout,*) 'No wood growth.  Not a problem!  Just an observation.'
                   ENDIF
                ELSE
                   WRITE(numout,*) 'WARNING 34: problem with Cr_target'
                   WRITE(numout,*) 'WARNING 34: PFT, ipts: ',j,ipts
                ENDIF
             ELSE
                WRITE(numout,*) 'WARNING 35: problem with Cs_target'
                WRITE(numout,*) 'WARNING 35: PFT, ipts: ',j,ipts
             ENDIF
          ELSE
             WRITE(numout,*) 'WARNING 36: problem with Cl_target'
             WRITE(numout,*) 'WARNING 36: PFT, ipts: ',j,ipts
          ENDIF
       ELSEIF(b_inc_tot .LT. -min_stomate) THEN 
          WRITE(numout,*) 'WARNING 37: problem with b_inc_tot'
          WRITE(numout,*) 'WARNING 37: PFT, ipts: ',j,ipts
       ELSE
          WRITE(numout,*) 'no unallocated fraction'
       ENDIF

    CASE(11)
       ! Ordinary allocation
       WRITE(numout,*) 'delta_ba, ', delta_ba
       IF ( (SUM(Cl_inc(:)) .GE. zero) .AND. (SUM(Cs_inc(:)) .GE. zero) .AND. &
            (SUM(Cr_inc(:)) .GE. zero) .AND. &
            ( b_inc_tot - SUM(circ_class_n(ipts,j,:) * (Cl_inc(:)+Cs_inc(:)+Cr_inc(:))) .GT. -1*min_stomate) .AND. &
            ( b_inc_tot - SUM(circ_class_n(ipts,j,:) * (Cl_inc(:)+Cs_inc(:)+Cr_inc(:))) .LT. min_stomate ) ) THEN
          WRITE(numout,*) 'Exc 10.1: Ordinary allocation was succesful'
          WRITE(numout,*) 'Cl_inc, Cs_inc, Cr_inc, unallocated', Cl_inc(:), Cs_inc(:), Cr_inc(:), & 
               b_inc_tot - SUM(circ_class_n(ipts,j,:) * (Cl_inc(:)+Cs_inc(:)+Cr_inc(:)))
       ELSE
          WRITE(numout,*) 'WARNING 38: Exc 10.2 problem with ordinary allocation'
          WRITE(numout,*) 'WARNING 38: PFT, ipts: ',j,ipts
          WRITE(numout,*) 'Cl_inc, Cs_inc, Cr_inc, unallocated', Cl_inc(:), Cs_inc(:), Cr_inc(:), & 
               b_inc_tot - SUM(circ_class_n(ipts,j,:) * (Cl_inc(:)+Cs_inc(:)+Cr_inc(:)))
       ENDIF

    CASE(12)
       ! Enough leaves and structure, grow roots
       WRITE(numout,*) 'Exc 1: Cl_incp(=0), Cs_incp (=0), Cr_incp (<>0), unallocated, '
       WRITE(numout,*) Cl_incp(1), Cs_incp(1), Cr_incp(1), &
            b_inc_tot - (SUM(circ_class_n(ipts,j,:))*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) )
       IF (b_inc_tot - (SUM(circ_class_n(ipts,j,:))*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) ) .LT. -10*EPSILON(zero)) THEN
          WRITE(numout,*) 'Exc 1.1: unallocated less then 0: overspending, ', b_inc_tot - &
               (SUM(circ_class_n(ipts,j,:))*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)))
       ELSE
          IF ( (b_inc_tot - ( SUM(circ_class_n(ipts,j,:))*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) ) .GE. -10*EPSILON(zero)) .AND. &
               ((SUM(circ_class_n(ipts,j,:)) * ABS(Cr_target(1)-Cr(1)-Cr_incp(1))) .LE. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 1.2: unallocated <>= 0 but tree is in good shape: successful allocation'
          ELSEIF ( (b_inc_tot - ( SUM(circ_class_n(ipts,j,:))*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) ) &
               .LE. min_stomate) .AND. &
               (SUM(circ_class_n(ipts,j,:)) * ABS(Cr_target(1)-Cr(1)-Cr_incp(1)) .GT. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 1.3: unallocated = 0, but the tree needs more reshaping: successful allocation'
          ELSE
             WRITE(numout,*) 'WARNING 39: Exc 1.4 unexpected result'
             WRITE(numout,*) 'WARNING 39: PFT, ipts: ',j,ipts
          ENDIF
       ENDIF

    CASE(13)
       ! Enough structural C and roots, grow leaves
       WRITE(numout,*) 'Exc 2: Cl_incp(<>0), Cs_incp (=0), Cr_incp (=0), unallocated, '
       WRITE(numout,*) Cl_incp(1), Cs_incp(1), Cr_incp(1), &
            b_inc_tot - (SUM(circ_class_n(ipts,j,:))*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) )
       IF (b_inc_tot - (SUM(circ_class_n(ipts,j,:))*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) ) .LT. -10*EPSILON(zero)) THEN
          WRITE(numout,*) 'Exc 2.1: unallocated less then 0: overspending, ', b_inc_tot - &
               (SUM(circ_class_n(ipts,j,:))*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)))
       ELSE
          IF ( (b_inc_tot - ( SUM(circ_class_n(ipts,j,:))*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) ) .GE. -10*EPSILON(zero)) .AND. &
               ((SUM(circ_class_n(ipts,j,:)) * ABS(Cl_target(1)-Cl(1)-Cl_incp(1))) .LE. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 2.2: unallocated <>= 0 but tree is in good shape: successful allocation'
          ELSEIF ( (b_inc_tot - ( SUM(circ_class_n(ipts,j,:))*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) ) &
               .LE. min_stomate) .AND. &
               ((SUM(circ_class_n(ipts,j,:)) * ABS(Cl_target(1)-Cl(1)-Cl_incp(1))) .GT. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 2.3: unallocated = 0, but the tree needs more reshaping: successful allocation'
          ELSE
             WRITE(numout,*) 'WARNING 40: Exc l.4 unexpected result'
             WRITE(numout,*) 'WARNING 40: PFT, ipts: ',j,ipts
          ENDIF
       ENDIF

    CASE(14)
       ! Enough structural C and root, grow leaves
       WRITE(numout,*) 'Exc 3: Cl_incp(<>0), Cs_incp(=0), Cr_incp(<>0), unallocated, '
       WRITE(numout,*) Cl_incp(1), Cs_incp(1), Cr_incp(1), b_inc_tot - & 
            (SUM(circ_class_n(ipts,j,:)) * (Cl_incp(1)+Cs_incp(1)+Cr_incp(1)))
       IF (b_inc_tot - SUM(circ_class_n(ipts,j,:)) * (Cl_incp(1) + Cs_incp(1) + Cr_incp(1))  &
            .LT. -10*EPSILON(zero)) THEN
          WRITE(numout,*) 'Exc 3.1: unallocated less then 0: overspending, ', b_inc_tot - &
               (SUM(circ_class_n(ipts,j,:)) * (Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) )
       ELSE
          IF ( (b_inc_tot - (SUM(circ_class_n(ipts,j,:))*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1))) &
               .GE. min_stomate) .AND. &
               ((SUM(circ_class_n(ipts,j,:)) * ABS(Cl_target(1)-Cl(1)-Cl_incp(1)) ) .LE. min_stomate) .AND. &
               ((SUM(circ_class_n(ipts,j,:)) * ABS(Cr_target(1)-Cr(1)-Cr_incp(1)) ) .LE. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 3.2: unallocated <>= 0 but tree is in good shape: successful allocation'
          ELSEIF ( (b_inc_tot - (SUM(circ_class_n(ipts,j,:))*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) ) &
               .LE. min_stomate) .AND. &
               ((SUM(circ_class_n(ipts,j,:)) * ABS(Cl_target(1)-Cl(1)-Cl_incp(1)) ) .GT. min_stomate) .AND. &
               ((SUM(circ_class_n(ipts,j,:)) * ABS(Cr_target(1)-Cr(1)-Cr_incp(1)) ) .GT. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 3.3: unallocated = 0, but the tree needs more reshaping: successful allocation'
          ELSE
             WRITE(numout,*) 'WARNING 41: Exc 3.4 unexpected result'
             WRITE(numout,*) 'WARNING 41: PFT, ipts: ',j,ipts
          ENDIF
       ENDIF

    CASE(15)
       ! Enough leaves and structural C, grow roots
       WRITE(numout,*) 'Exc 4: Cl_incp(=0), Cs_incp (=0), Cr_incp (<>0), unallocated, '
       WRITE(numout,*) Cl_incp(1), Cs_incp(1), Cr_incp(1), &
            b_inc_tot - (SUM(circ_class_n(ipts,j,:))*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) )
       IF (b_inc_tot - (SUM(circ_class_n(ipts,j,:))*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) ) .LT. -10*EPSILON(zero)) THEN
          WRITE(numout,*) 'Exc 4.1: unallocated less then 0: overspending, ', b_inc_tot - &
               (SUM(circ_class_n(ipts,j,:))*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)))
       ELSE
          IF ( (b_inc_tot - ( SUM(circ_class_n(ipts,j,:))*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) ) .GE. -10*EPSILON(zero)) .AND. &
               ((SUM(circ_class_n(ipts,j,:)) * ABS(Cr_target(1)-Cr(1)-Cr_incp(1))) .LE. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 4.2: unallocated <>= 0 but tree is in good shape: successful allocation'
          ELSEIF ( (b_inc_tot - ( SUM(circ_class_n(ipts,j,:))*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) ) &
               .LE. min_stomate) .AND. &
               ((SUM(circ_class_n(ipts,j,:)) * ABS(Cr_target(1)-Cr(1)-Cr_incp(1))) .GT. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 4.3: unallocated = 0, but the tree needs more reshaping: successful allocation'
          ELSE
             WRITE(numout,*) 'WARNING 42: Exc 4.4 unexpected result'
             WRITE(numout,*) 'WARNING 42: PFT, ipts: ',j,ipts
          ENDIF
       ENDIF

    CASE(16)
       ! Enough leaves and roots, grow structural C            
       WRITE(numout,*) 'Exc 5: Cl_incp(=0), Cs_incp (<>0), Cr_incp (=0), unallocated, '
       WRITE(numout,*) Cl_incp(1), Cs_incp(1), Cr_incp(1), &
            b_inc_tot - (SUM(circ_class_n(ipts,j,:))*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) )
       IF (b_inc_tot - (SUM(circ_class_n(ipts,j,:))*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) ) .LT. -10*EPSILON(zero)) THEN
          WRITE(numout,*) 'Exc 5.1: unallocated less then 0: overspending, ', b_inc_tot - &
               (SUM(circ_class_n(ipts,j,:))*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)))
       ELSE
          IF ( (b_inc_tot - ( SUM(circ_class_n(ipts,j,:))*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) ) .GE. -10*EPSILON(zero)) .AND. &
               ((SUM(circ_class_n(ipts,j,:)) * ABS(Cs_target(1)-Cs(1)-Cs_incp(1))) .LE. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 5.2: unallocated <>= 0 but tree is in good shape: successful allocation'
          ELSEIF ( (b_inc_tot - ( SUM(circ_class_n(ipts,j,:))*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) ) &
               .LE. min_stomate) .AND. &
               ((SUM(circ_class_n(ipts,j,:)) * ABS(Cs_target(1)-Cs(1)-Cs_incp(1))) .GT. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 5.3: unallocated = 0, but the tree needs more reshaping: successful allocation'
          ELSE
             WRITE(numout,*) 'WARNING 43: Exc 5.4 unexpected result'
             WRITE(numout,*) 'WARNING 43: PFT, ipts: ',j,ipts
          ENDIF
       ENDIF

    CASE(17)
       ! Enough leaves, grow structural C and roots
       WRITE(numout,*) 'Exc 6: Cl_incp(=0), Cs_incp(<>0), Cr_incp(<>0), unallocated'
       WRITE(numout,*) Cl_incp(1), Cs_incp(1), Cr_incp(1), &
            b_inc_tot - (SUM(circ_class_n(ipts,j,:)) * (Cl_incp(1)+Cs_incp(1)+Cr_incp(1)))
       IF (b_inc_tot - (SUM(circ_class_n(ipts,j,:)) * (Cl_incp(1)+Cs_incp(1)+Cr_incp(1))) .LT. -10*EPSILON(zero)) THEN
          WRITE(numout,*) 'Exc 6.1: unallocated less then 0: overspending, ', &
               b_inc_tot - (SUM(circ_class_n(ipts,j,:)) * (Cl_incp(1)+Cs_incp(1)+Cr_incp(1)))
       ELSE
          IF ( (b_inc_tot - SUM(circ_class_n(ipts,j,:))*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) .GE. -10*EPSILON(zero)) .AND. &
               ((SUM(circ_class_n(ipts,j,:)) * ABS(Cs_target(1)-Cs(1)-Cs_incp(1))) .LE. min_stomate) .AND. &
               ((SUM(circ_class_n(ipts,j,:)) * ABS(Cr_target(1)-Cr(1)-Cr_incp(1))) .LE. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 6.2: unallocated <>= 0 but tree is in good shape: successful allocation'
          ELSEIF ( (b_inc_tot - SUM(circ_class_n(ipts,j,:))*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) .LE. min_stomate) .AND. &
               ((SUM(circ_class_n(ipts,j,:)) * ABS(Cs_target(1)-Cs(1)-Cs_incp(1))) .GT. min_stomate) .AND. &
               ((SUM(circ_class_n(ipts,j,:)) * ABS(Cr_target(1)-Cr(1)-Cr_incp(1))) .GT. min_stomate) ) THEN
             WRITE(numout,*) &
                  'Exc 6.3: unallocated = 0, but the tree needs more reshaping: successful allocation'
          ELSE
             WRITE(numout,*) 'WARNING 44: Exc 6.4 unexpected result'
             WRITE(numout,*) 'WARNING 44: PFT, ipts: ',j,ipts
          ENDIF
       ENDIF

    CASE(18)
       ! Enough leaves and structural C, grow roots
       WRITE(numout,*) 'Exc 7: Cl_incp(=0), Cs_incp (=0), Cr_incp (<>0), unallocated, '
       WRITE(numout,*) Cl_incp(1), Cs_incp(1), Cr_incp(1), &
            b_inc_tot - (SUM(circ_class_n(ipts,j,:))*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) )
       IF (b_inc_tot - (SUM(circ_class_n(ipts,j,:))*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) ) .LT. -10*EPSILON(zero)) THEN
          WRITE(numout,*) 'Exc 7.1: unallocated less then 0: overspending, ', b_inc_tot - &
               (SUM(circ_class_n(ipts,j,:))*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)))
       ELSE
          IF ( (b_inc_tot - ( SUM(circ_class_n(ipts,j,:))*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) ) .GE. -10*EPSILON(zero)) .AND. &
               ((SUM(circ_class_n(ipts,j,:)) * ABS(Cl_target(1)-Cl(1)-Cl_incp(1))) .LE. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 7.2: unallocated <>= 0 but tree is in good shape: successful allocation'
          ELSEIF ( (b_inc_tot - ( SUM(circ_class_n(ipts,j,:))*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) ) &
               .LE. min_stomate) .AND. &
               ((SUM(circ_class_n(ipts,j,:)) * ABS(Cl_target(1)-Cl(1)-Cl_incp(1))) .GT. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 7.3: unallocated = 0, but the tree needs more reshaping: successful allocation'
          ELSE
             WRITE(numout,*) 'WARNING 45: Exc 7.4 unexpected result'
             WRITE(numout,*) 'WARNING 45: PFT, ipts: ',j,ipts
          ENDIF
       ENDIF

    CASE(19)
       ! Enough leaves and roots, grow structural C
       WRITE(numout,*) 'Exc 8: Cl_incp(=0), Cs_incp (<>0), Cr_incp (=0), unallocated, '
       WRITE(numout,*) Cl_incp(1), Cs_incp(1), Cr_incp(1), &
            b_inc_tot - (SUM(circ_class_n(ipts,j,:))*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) )
       IF (b_inc_tot - (SUM(circ_class_n(ipts,j,:))*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) ) .LT. -10*EPSILON(zero)) THEN
          WRITE(numout,*) 'Exc 8.1: unallocated less then 0: overspending, ', b_inc_tot - &
               (SUM(circ_class_n(ipts,j,:))*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)))
       ELSE
          IF ( (b_inc_tot - ( SUM(circ_class_n(ipts,j,:))*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) ) .GE. -10*EPSILON(zero)) .AND. &
               ((SUM(circ_class_n(ipts,j,:)) * ABS(Cs_target(1)-Cs(1)-Cs_incp(1))) .LE. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 8.2: unallocated <>= 0 but tree is in good shape: successful allocation'
          ELSEIF ( (b_inc_tot - ( SUM(circ_class_n(ipts,j,:))*(Cl_incp(1)+Cs_incp(1)+Cr_incp(1)) ) &
               .LE. min_stomate) .AND. &
               ((SUM(circ_class_n(ipts,j,:)) * ABS(Cs_target(1)-Cs(1)-Cs_incp(1))) .GT. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 8.3: unallocated = 0, but the tree needs more reshaping: successful allocation'
          ELSE
             WRITE(numout,*) 'WARNING 46: Exc 8.4 unexpected result'
             WRITE(numout,*) 'WARNING 46: PFT, ipts: ',j,ipts
          ENDIF
       ENDIF

    CASE(20)
       ! Enough roots, grow structural C and leaves
       WRITE(numout,*) 'Exc 9: Cl_incp(<>0), Cs_incp(<>0), Cr_incp(=0), unallocated, '
       WRITE(numout,*) Cl_incp(1), Cs_incp(1), Cr_incp(1), &
            b_inc_tot - (SUM(circ_class_n(ipts,j,:)) * (Cl_incp(1)+Cs_incp(1)+Cr_incp(1)))
       WRITE(numout,*) 'term 1', b_inc_tot - (SUM(circ_class_n(ipts,j,:)) * (Cl_incp(1)+Cs_incp(1)+Cr_incp(1)))
       WRITE(numout,*) 'term 2', (SUM(circ_class_n(ipts,j,:)) * ABS(Cl_target(1)-Cl(1)-Cl_incp(1)))
       WRITE(numout,*) 'term 3', (SUM(circ_class_n(ipts,j,:)) * ABS(Cs_target(1)-Cs(1)-Cs_incp(1)))
       IF (b_inc_tot - (SUM(circ_class_n(ipts,j,:)) * (Cl_incp(1)+Cs_incp(1)+Cr_incp(1))) .LT. -10*EPSILON(zero)) THEN
          WRITE(numout,*) 'Exc 9.1: unallocated less then 0: overspending, ', &
               b_inc_tot - (SUM(circ_class_n(ipts,j,:)) * (Cl_incp(1)+Cs_incp(1)+Cr_incp(1)))
       ELSE
          IF ( (b_inc_tot - (SUM(circ_class_n(ipts,j,:)) * (Cl_incp(1)+Cs_incp(1)+Cr_incp(1))) .GE. -10*EPSILON(zero)) .AND. &
               ((SUM(circ_class_n(ipts,j,:)) * ABS(Cl_target(1)-Cl(1)-Cl_incp(1))) .LE. min_stomate) .AND. &
               ((SUM(circ_class_n(ipts,j,:)) * ABS(Cs_target(1)-Cs(1)-Cs_incp(1))) .LE. min_stomate) ) THEN
             WRITE(numout,*) 'Exc 9.2: unallocated <>= 0 but tree is in good shape: successful allocation'
          ELSEIF ( (b_inc_tot - (SUM(circ_class_n(ipts,j,:)) * (Cl_incp(1)+Cs_incp(1)+Cr_incp(1))) .LE. min_stomate) .AND. &
               (((SUM(circ_class_n(ipts,j,:)) * ABS(Cl_target(1)-Cl(1)-Cl_incp(1))) .GT. min_stomate) .OR. &
               ((SUM(circ_class_n(ipts,j,:)) * ABS(Cs_target(1)-Cs(1)-Cs_incp(1))) .GT. min_stomate) ) ) THEN
             WRITE(numout,*) 'Exc 9.3: unallocated = 0, but the tree needs more reshaping: successful allocation'
          ELSE
             WRITE(numout,*) 'WARNING 47: Exc 9.4 unexpected result'
             WRITE(numout,*) 'WARNING 47: PFT, ipts: ',j,ipts
          ENDIF
       ENDIF

    CASE(21)
       ! Ready for ordinary allocation
       WRITE(numout,*) 'Ready for ordinary allocation?'
       WRITE(numout,*) 'KF, LF, ', KF(ipts,j), LF(ipts,j)
       WRITE(numout,*) 'b_inc_tot, ', b_inc_tot
       WRITE(numout,*) 'Cl, Cs, Cr', Cl(1), Cs(1), Cr(1)
       WRITE(numout,*) 'Cl_target-Cl, ', Cl_target(1)-Cl(1)
       WRITE(numout,*) 'Cs_target-Cs, ', Cs_target(1)-Cs(1)
       WRITE(numout,*) 'Cr_target-Cr, ', Cr_target(1)-Cr(1)
       IF (b_inc_tot .GT. min_stomate) THEN
          IF (ABS(Cl_target(1)-Cl(1)) .LE. min_stomate) THEN
             IF (ABS(Cs_target(1)-Cs(1)) .LE. min_stomate) THEN
                IF (ABS(Cr_target(1)-Cr(1)) .LE. min_stomate) THEN
                   IF (b_inc_tot .GT. min_stomate) THEN
                      IF (grow_wood) THEN
                         WRITE(numout,*) 'should result in exc 10.1 or 10.2'
                      ELSE
                         WRITE(numout,*) 'WARNING 48: no wood growth'
                         WRITE(numout,*) 'WARNING 48: PFT, ipts: ',j,ipts
                      ENDIF
                   ENDIF
                ELSE
                   WRITE(numout,*) 'WARNING 49: problem with Cr_target'
                   WRITE(numout,*) 'WARNING 49: PFT, ipts: ',j,ipts
                ENDIF
             ELSE
                WRITE(numout,*) 'WARNING 50: problem with Cs_target'
                WRITE(numout,*) 'WARNING 50: PFT, ipts: ',j,ipts
             ENDIF
          ELSE
             WRITE(numout,*) 'WARNING 51: problem with Cl_target'
             WRITE(numout,*) 'WARNING 51: PFT, ipts: ',j,ipts
          ENDIF
       ELSEIF(b_inc_tot .LT. -min_stomate) THEN 
          WRITE(numout,*) 'WARNING 52: problem with b_inc_tot'
          WRITE(numout,*) 'WARNING 52: PFT, ipts: ',j,ipts
       ELSE
          WRITE(numout,*) 'no unallocated fraction'
       ENDIF

    CASE(22)
       ! Ordinary allocation
       IF ( ((Cl_inc(1)) .GE. zero) .AND. ((Cs_inc(1)) .GE. zero) .AND. &
            ((Cr_inc(1)) .GE. zero) .AND. &
            ( b_inc_tot - (SUM(circ_class_n(ipts,j,:)) * (Cl_inc(1)+Cs_inc(1)+Cr_inc(1))) .LT. min_stomate ) ) THEN
          WRITE(numout,*) 'Exc 10.1: Ordinary allocation was succesful'
          WRITE(numout,*) 'Cl_inc, Cs_inc, Cr_inc, unallocated', Cl_inc(1), Cs_inc(1), Cr_inc(1), & 
               b_inc_tot - (SUM(circ_class_n(ipts,j,:)) * (Cl_inc(1)+Cs_inc(1)+Cr_inc(1)))
       ELSE
          WRITE(numout,*) 'WARNING 53: Exc 10.2 problem with ordinary allocation'
          WRITE(numout,*) 'WARNING 53: PFT, ipts: ',j,ipts
          WRITE(numout,*) 'Cl_inc, Cs_inc, Cr_inc, unallocated', Cl_inc(1), Cs_inc(1), Cr_inc(1), & 
               b_inc_tot - (SUM(circ_class_n(ipts,j,:)) * (Cl_inc(1)+Cs_inc(1)+Cr_inc(1)))
       ENDIF

    END SELECT

  END SUBROUTINE comment

!! ================================================================================================================================
!! SUBROUTINE	: add_fluxes_to_pools
!!
!>\BRIEF        Add carbon and nitrogen fluxes to pools before veget_max changes
!!
!!
!! DESCRIPTION  : GPP, plant_n_uptake and n_deposition represent incoming
!!  carbon and nitrogen fluxes that were calculated sechiba or
!!  stomate. They were calculated for veget_max. Excessive
!!  complexity can be avoided by adding these fluxes to their
!!  target pools before accounting for the processes in which
!!  veget_max may change, i.e., land cover change, stand replacing
!!  disturbances and forestry (when using age classes).
!!
!! RECENT CHANGE(S): 
!!
!! MAIN OUTPUT VARIABLE(S): circ_class_biomass, soil_n_min
!!
!! REFERENCE(S)	: none
!! 
!! FLOWCHART : 
!! \n
!_ ================================================================================================================================


END MODULE stomate_growth_fun_all
