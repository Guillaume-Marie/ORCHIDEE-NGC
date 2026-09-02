! =================================================================================================================================
! MODULE       : pft_parameters
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2011)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        This module initializes all the pft parameters in function of the
!!              number of vegetation types and of the values chosen by the user.
!!
!!\n DESCRIPTION:  This module allocates and initializes the pft parameters in function of the number of pfts
!!                 and the values of the parameters. \n
!!                 The number of PFTs is read in control.f90 (subroutine control_initialize). \n
!!                 Then we can initialize the parameters. \n
!!                 This module is the result of the merge of constantes_co2, constantes_veg, stomate_constants.\n
!!
!! RECENT CHANGE(S): Josefine Ghattas 2013 : The declaration part has been extracted and moved to module pft_parameters_var
!!
!! REFERENCE(S)	: None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE/src_parameters/pft_parameters.f90 $
!! $Date: 2026-05-12 21:40:35 +0200 (mar. 12 mai 2026) $
!! $Revision: 9536 $
!! \n
!_ ================================================================================================================================

MODULE pft_parameters

  USE pft_parameters_var
  USE vertical_soil_var
  USE constantes_soil_var
  USE constantes_mtc
  USE constantes_var
  USE ioipsl
  USE ioipsl_para 
  USE defprec

  IMPLICIT NONE

CONTAINS
  

!! ================================================================================================================================
!! SUBROUTINE   : pft_parameters_main
!!
!>\BRIEF          This subroutine initializes all the pft parameters in function of the
!! number of vegetation types chosen by the user.
!!
!! DESCRIPTION  : This subroutine is called after the reading of the number of PFTS and the options 
!!                activated by the user in the configuration files. \n
!!                The allocation is done just before reading the correspondence table  between PFTs and MTCs
!!                defined by the user in the configuration file.\n
!!                With the correspondence table, the subroutine can initialize the pft parameters in function
!!                of the flags activated (ok_sechiba, ok_stomate, routing,...) in order to
!!                optimize the memory allocation. \n
!!                If the number of PFTs and pft_to_mtc are not found, the standard configuration will be used
!!                (13 PFTs, PFT = MTC). \n 
!!                Some restrictions : the pft 1 can only be the bare soil and it is unique. \n
!!                Algorithm : Build new PFT from 13 generic-PFT or meta-classes.
!!                1. Read the number of PFTs in "run.def". If nothing is found, it is assumed that the user intend to use 
!!                   the standard of PFTs (13).
!!                2. Read the index vector in "run.def". The index vector associates one PFT to one meta-classe (or generic PFT).
!!                   When the association is done, the PFT defined by the user inherited the default values from the meta classe.
!!                   If nothing is found, it is assumed to use the standard index vector (PFT = MTC).
!!                3. Check consistency
!!                4. Memory allocation and initialization.
!!                5. The parameters are read in the configuration file in config_initialize (control module).
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE pft_parameters_main()

    IMPLICIT NONE

    !! 0. Variables and parameters declaration

    !! 0.4 Local variables  

    INTEGER(i_std) :: j                             !! Index (unitless)
    

    !_ ================================================================================================================================ 

    !
    ! PFT global
    !

    IF(l_first_pft_parameters) THEN

       !! 1. First time step
       IF(printlev>=3) THEN
          WRITE(numout,*) 'l_first_pft_parameters :we read the parameters from the def files'
       ENDIF

       !! 2. Memory allocation for the pfts-parameters
       CALL pft_parameters_alloc()

       !! 3. Correspondance table 

       !! 3.1 Initialisation of the correspondance table
       !! Initialisation of the correspondance table
       IF (nvm == nvmc) THEN
          pft_to_mtc = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 /)
       ELSE
          pft_to_mtc(:) = undef_int
       ENDIF !(nvm  == nvmc)

       !! 3.2 Reading of the conrrespondance table in the .def file
       !
       !Config Key   = PFT_TO_MTC
       !Config Desc  = correspondance array linking a PFT to MTC
       !Config if    = OK_SECHIBA or OK_STOMATE
       !Config Def   = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13
       !Config Help  =
       !Config Units = [-]
       CALL getin_p('PFT_TO_MTC',pft_to_mtc)

       !! 3.3 If the user want to use the standard configuration, he needn't to fill the correspondance array
       !!     If the configuration is wrong, send a error message to the user.
       IF(nvm /= nvmc ) THEN
          !
          IF(pft_to_mtc(1) == undef_int) THEN
             CALL ipslerr_p (3,'The array PFT_TO_MTC is empty : we stop','','','')
          ENDIF !(pft_to_mtc(1) == undef_int)
          !
       ENDIF !(nvm /= nvmc )

       !! 3.4 Some error messages

       !! 3.4.1 What happened if pft_to_mtc(j) > nvmc or pft_to_mtc(j) <=0 (if the mtc doesn't exist)?
       DO j = 1, nvm ! Loop over # PFTs  
          !
          IF( (pft_to_mtc(j) > nvmc) .OR. (pft_to_mtc(j) <= 0) ) THEN
             CALL ipslerr_p(3,'the metaclass chosen does not exist', &
                  'we stop reading pft_to_mtc','','') 
          ENDIF !( (pft_to_mtc(j) > nvmc) .OR. (pft_to_mtc(j) <= 0) )
          !
       ENDDO  ! Loop over # PFTs  


       !! 3.4.2 Check if pft_to_mtc(1) = 1 
       IF(pft_to_mtc(1) /= 1) THEN
          !
          CALL ipslerr_p(3,'the first pft has to be the bare soil', & 
               'we stop reading next values of pft_to_mtc','','') 
          !
       ELSE
          !
          DO j = 2,nvm ! Loop over # PFTs different from bare soil
             !
             IF(pft_to_mtc(j) == 1) THEN
                CALL ipslerr_p(3,'only pft_to_mtc(1) has to be the bare soil',&
                     'we stop reading pft_to_mtc','','')
             ENDIF ! (pft_to_mtc(j) == 1)
             !
          ENDDO ! Loop over # PFTs different from bare soil
          !
       ENDIF !(pft_to_mtc(1) /= 1)


       !! 4.Initialisation of the pfts-parameters
       CALL pft_parameters_init()

       !! 5. Useful data

       !! 5.1 Read the name of the PFTs given by the user
       !
       !Config Key   = PFT_NAME
       !Config Desc  = Name of a PFT
       !Config if    = OK_SECHIBA or OK_STOMATE
       !Config Def   = bare ground, tropical broad-leaved evergreen, tropical broad-leaved raingreen, 
       !Config         temperate needleleaf evergreen, temperate broad-leaved evergreen temperate broad-leaved summergreen,
       !Config         boreal needleleaf evergreen, boreal broad-leaved summergreen, boreal needleleaf summergreen,
       !Config         C3 grass, C4 grass, C3 agriculture, C4 agriculture    
       !Config Help  = the user can name the new PFTs he/she introducing for new species
       !Config Units = [-]
       CALL getin_p('PFT_NAME',pft_name)

       !! 5.2 A useful message to the user: correspondance between the number of the pft
       !! and the name of the associated mtc 
       IF (printlev >=3 ) THEN
          WRITE(numout,*) ''
          DO j = 2,nvm ! Loop over # PFTs
             WRITE(numout,*) 'The PFT',j, 'called ', trim(PFT_name(j)),' corresponds to the MTC : ',trim(MTC_name(pft_to_mtc(j)))
             CALL flush(numout)
          END DO
          WRITE(numout,*) ''
       END IF


       !! 6. End message
       IF (printlev>=3) THEN
          WRITE(numout,*) 'pft_parameters_done'
          CALL flush(numout)
       END IF
          
       !! 8. Reset flag
       l_first_pft_parameters = .FALSE.

    ELSE 

       RETURN

    ENDIF !(l_first_pft_parameters)

  END SUBROUTINE pft_parameters_main


!! ================================================================================================================================
!! SUBROUTINE   : pft_parameters_init 
!!
!>\BRIEF          This subroutine initializes all the pft parameters by the default values
!! of the corresponding metaclasse. 
!!
!! DESCRIPTION  : This subroutine is called after the reading of the number of PFTS and the correspondence
!!                table defined by the user in the configuration files. \n
!!                With the correspondence table, the subroutine can search the default values for the parameter
!!                even if the PFTs are classified in a random order (except bare soil). \n
!!                With the correspondence table, the subroutine can initialize the pft parameters in function
!!                of the flags activated (ok_sechiba, ok_stomate, routing,...).\n
!!
!! RECENT CHANGE(S): Didier Solyga : Simplified PFT loops : use vector notation. 
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE pft_parameters_init()

    IMPLICIT NONE

    !! 0. Variables and parameters declaration

    !! 0.1 Input variables

    !! 0.4 Local variables

    INTEGER(i_std)                :: jv                !! Index (unitless)
    
    !_ ================================================================================================================================ 

    !
    ! 1. Correspondance between the PFTs values and thes MTCs values 
    !


    ! 1.1 For parameters used anytime

    PFT_name(:) = MTC_name(pft_to_mtc(:))
    !
    ! Vegetation structure 
    !
    llaimax(:) = llaimax_mtc(pft_to_mtc(:))
    llaimin(:) = llaimin_mtc(pft_to_mtc(:))
    height_presc(:) = height_presc_mtc(pft_to_mtc(:))
    z0_over_height(:) = z0_over_height_mtc(pft_to_mtc(:))
    ratio_z0m_z0h(:) = ratio_z0m_z0h_mtc(pft_to_mtc(:))
    type_of_lai(:) = type_of_lai_mtc(pft_to_mtc(:))
    natural(:) = natural_mtc(pft_to_mtc(:))
    !
    ! Water - sechiba
    !
    IF (zmaxh == 2.0) THEN
       IF (printlev>=2) WRITE(numout,*)'Initialize humcst using reference values for 2m soil depth'
       humcste(:) = humcste_ref2m(pft_to_mtc(:))  ! values for 2m soil depth
       max_root_depth(:) = max_root_depth_ref2m_mtc(pft_to_mtc(:))
    ELSE IF (zmaxh == 4.0) THEN
       IF (printlev>=2) WRITE(numout,*)'Initialize humcst using reference values for 4m soil depth'
       humcste(:) = humcste_ref4m(pft_to_mtc(:))  ! values for 4m soil depth
       max_root_depth(:) = max_root_depth_ref4m_mtc(pft_to_mtc(:))
    ELSE
       IF (printlev>=2) WRITE(numout,*)'Note that humcste is initialized with values for 2m soil depth bur zmaxh=', zmaxh
       humcste(:) = humcste_ref2m(pft_to_mtc(:))  ! values for 2m soil depth
    END IF

    !
    ! Soil - vegetation
    !
    pref_soil_veg(:) = pref_soil_veg_mtc(pft_to_mtc(:))

    !
    ! Vegetation - age classes
    !
    agec_group(:) = agec_group_mtc(pft_to_mtc(:))
    !
    ! Peatland
    !
    is_peat(:)=is_peat_mtc(pft_to_mtc(:))
    is_croppeat(:)=is_croppeat_mtc(pft_to_mtc(:))
    !
    ! Photosynthesis
    !
    is_c4(:) = is_c4_mtc(pft_to_mtc(:))
    vcmax_fix(:) = vcmax_fix_mtc(pft_to_mtc(:))
    downregulation_co2_coeff(:) = downregulation_co2_coeff_mtc(pft_to_mtc(:))
    E_KmC(:)      = E_KmC_mtc(pft_to_mtc(:))
    E_KmO(:)      = E_KmO_mtc(pft_to_mtc(:))
    E_Sco(:)      = E_Sco_mtc(pft_to_mtc(:))
    E_gamma_star(:) = E_gamma_star_mtc(pft_to_mtc(:))
    E_Vcmax(:)    = E_Vcmax_mtc(pft_to_mtc(:))
    E_Jmax(:)     = E_Jmax_mtc(pft_to_mtc(:))
    aSV(:)        = aSV_mtc(pft_to_mtc(:))
    bSV(:)        = bSV_mtc(pft_to_mtc(:))
    tphoto_min(:) = tphoto_min_mtc(pft_to_mtc(:))
    tphoto_max(:) = tphoto_max_mtc(pft_to_mtc(:))
    aSJ(:)        = aSJ_mtc(pft_to_mtc(:))
    bSJ(:)        = bSJ_mtc(pft_to_mtc(:))
    D_Vcmax(:)     = D_Vcmax_mtc(pft_to_mtc(:))
    D_Jmax(:)     = D_Jmax_mtc(pft_to_mtc(:))
    E_gm(:)       = E_gm_mtc(pft_to_mtc(:)) 
    S_gm(:)       = S_gm_mtc(pft_to_mtc(:)) 
    D_gm(:)       = D_gm_mtc(pft_to_mtc(:)) 
    E_Rd(:)       = E_Rd_mtc(pft_to_mtc(:))
    Vcmax25(:)    = Vcmax25_mtc(pft_to_mtc(:))
    arJV(:)       = arJV_mtc(pft_to_mtc(:))
    brJV(:)       = brJV_mtc(pft_to_mtc(:))
    KmC25(:)      = KmC25_mtc(pft_to_mtc(:))
    KmO25(:)      = KmO25_mtc(pft_to_mtc(:))
    Sco25(:)      = Sco25_mtc(pft_to_mtc(:))
    gm25(:)       = gm25_mtc(pft_to_mtc(:)) 
    gamma_star25(:)  = gamma_star25_mtc(pft_to_mtc(:))
    wood_wat_store_fac(:) = wood_wat_store_fac_mtc(pft_to_mtc(:))
    leaf_wat_store_fac(:) = leaf_wat_store_fac_mtc(pft_to_mtc(:))
    a1(:)         = a1_mtc(pft_to_mtc(:))
    b1(:)         = b1_mtc(pft_to_mtc(:))
    a_fpsi(:)     = a_fpsi_mtc(pft_to_mtc(:))
    sf(:)         = sf_mtc(pft_to_mtc(:))
    Rroot_sup(:)  = Rroot_sup_mtc(pft_to_mtc(:))
    Rroot_inf(:)  = Rroot_inf_mtc(pft_to_mtc(:))
    Rleaf(:)      = Rleaf_mtc(pft_to_mtc(:))
    Rxylem(:)     = Rxylem_mtc(pft_to_mtc(:))
    Rroot_soil(:) = Rroot_soil_mtc(pft_to_mtc(:))
    Rsto_leaf(:) = Rsto_leaf_mtc(pft_to_mtc(:))
    Rsto_wood(:) = Rsto_wood_mtc(pft_to_mtc(:))
    psi_ref_g0(:) = psi_ref_g0_mtc(pft_to_mtc(:))
    psi_ref_sto_leaf(:) = psi_ref_sto_leaf_mtc(pft_to_mtc(:))
    psi_ref_sto_wood(:) = psi_ref_sto_wood_mtc(pft_to_mtc(:))
    psi_50_res_leaf(:) = psi_50_res_leaf_mtc(pft_to_mtc(:))
    psi_50_res_wood(:) = psi_50_res_wood_mtc(pft_to_mtc(:))
    psi_50_res_root(:) = psi_50_res_root_mtc(pft_to_mtc(:))
    a_leaf(:) = a_leaf_mtc(pft_to_mtc(:))
    a_wood(:) = a_wood_mtc(pft_to_mtc(:))
    a_root(:) = a_root_mtc(pft_to_mtc(:))
    k_leaf_max(:) = k_leaf_max_mtc(pft_to_mtc(:))
    k_wood_max(:) = k_wood_max_mtc(pft_to_mtc(:))
    sapwood_density(:) = sapwood_density_mtc(pft_to_mtc(:))
    ks_max(:) = ks_max_mtc(pft_to_mtc(:))
    k_root_max(:) = k_root_max_mtc(pft_to_mtc(:))
    lambda_leaf(:) = lambda_leaf_mtc(pft_to_mtc(:))
    lambda_wood(:) = lambda_wood_mtc(pft_to_mtc(:))
    g0(:)         = g0_mtc(pft_to_mtc(:))
    d_l(:)        = d_l_mtc(pft_to_mtc(:))
    h_protons(:)  = h_protons_mtc(pft_to_mtc(:))
    fpsir(:)      = fpsir_mtc(pft_to_mtc(:))
    fQ(:)         = fQ_mtc(pft_to_mtc(:))     
    fpseudo(:)    = fpseudo_mtc(pft_to_mtc(:))    
    kp(:)         = kp_mtc(pft_to_mtc(:))
    alpha(:)      = alpha_mtc(pft_to_mtc(:))
    gbs(:)        = gbs_mtc(pft_to_mtc(:))
    theta(:)      = theta_mtc(pft_to_mtc(:))        
    alpha_LL(:)   = alpha_LL_mtc(pft_to_mtc(:))
    stress_vcmax(:) = stress_vcmax_mtc(pft_to_mtc(:))
    stress_gs(:)    = stress_gs_mtc(pft_to_mtc(:))
    stress_gm(:)    = stress_gm_mtc(pft_to_mtc(:))
    ext_coeff(:) = ext_coeff_mtc(pft_to_mtc(:))
    ext_coeff_vegetfrac(:) = ext_coeff_vegetfrac_mtc(pft_to_mtc(:))
    !
    !! Define labels from physiologic characteristics 
    !
    leaf_tab(:) = leaf_tab_mtc(pft_to_mtc(:)) 
    pheno_model(:) = pheno_model_mtc(pft_to_mtc(:))
    gdd_threshold(:) = gdd_threshold_mtc(pft_to_mtc(:))
    !
    is_tree(:) = .FALSE.
    DO jv = 1,nvm
       IF ( leaf_tab(jv) <= 2 ) is_tree(jv) = .TRUE.
    END DO
    !
    is_deciduous(:) = .FALSE.
    DO jv = 1,nvm
       IF ( is_tree(jv) .AND. (pheno_model(jv) /= "none") ) is_deciduous(jv) = .TRUE.
    END DO
    !
    is_evergreen(:) = .FALSE.
    DO jv = 1,nvm
       IF ( is_tree(jv) .AND. (pheno_model(jv) == "none") ) is_evergreen(jv) = .TRUE.
    END DO
    !
    is_needleleaf(:) = .FALSE.
    DO jv = 1,nvm
       IF ( leaf_tab(jv) == 2 ) is_needleleaf(jv) = .TRUE.
    END DO

    is_tropical(:) = is_tropical_mtc(pft_to_mtc(:)) 
    is_temperate(:) = is_temperate_mtc(pft_to_mtc(:)) 
    is_boreal(:) = is_boreal_mtc(pft_to_mtc(:)) 

    ! 1.2 For sechiba parameters
    IF (ok_sechiba) THEN
       !
       ! Vegetation structure - sechiba
       !
       rveg_pft(:) = rveg_mtc(pft_to_mtc(:))
       !
       ! Evapotranspiration -  sechiba
       !
       rstruct_const(:) = rstruct_const_mtc(pft_to_mtc(:))
       kzero(:) = kzero_mtc(pft_to_mtc(:))
       alpha_watstress_pft(:) = alpha_watstress_pft_mtc(pft_to_mtc(:))
       !
       ! Water - sechiba
       !
       wmax_veg(:) = wmax_veg_mtc(pft_to_mtc(:))
       IF ( OFF_LINE_MODE ) THEN
          throughfall_by_pft(:) = 0.
       ELSE
          throughfall_by_pft(:) = throughfall_by_mtc(pft_to_mtc(:))
       ENDIF

       rsoil_scale(:) = rsoil_scale_mtc(pft_to_mtc(:))
       
       !
       ! Albedo - sechiba
       !
       snowa_aged_vis(:) = snowa_aged_vis_mtc(pft_to_mtc(:))
       snowa_aged_nir(:) = snowa_aged_nir_mtc(pft_to_mtc(:))
       snowa_dec_vis(:) = snowa_dec_vis_mtc(pft_to_mtc(:)) 
       snowa_dec_nir(:) = snowa_dec_nir_mtc(pft_to_mtc(:)) 
       alb_leaf_vis(:) = alb_leaf_vis_mtc(pft_to_mtc(:))  
       alb_leaf_nir(:) = alb_leaf_nir_mtc(pft_to_mtc(:))
       
       leaf_ssa(:,ivis) = leaf_ssa_vis_mtc(pft_to_mtc(:))
       leaf_ssa(:,inir) = leaf_ssa_nir_mtc(pft_to_mtc(:))
       leaf_psd(:,ivis) = leaf_psd_vis_mtc(pft_to_mtc(:))
       leaf_psd(:,inir) = leaf_psd_nir_mtc(pft_to_mtc(:))
       bgd_reflectance(:,ivis) = bgd_reflectance_vis_mtc(pft_to_mtc(:))
       bgd_reflectance(:,inir) = bgd_reflectance_nir_mtc(pft_to_mtc(:))
       leaf_to_shoot_clumping(:) = leaf_to_shoot_clumping_mtc(pft_to_mtc(:))
       lai_correction_factor(:) = lai_correction_factor_mtc(pft_to_mtc(:))
       min_level_sep(:) = min_level_sep_mtc(pft_to_mtc(:))
       !
       ! Diffuco and hydrol_arch
       !
       lai_top(:) = lai_top_mtc(pft_to_mtc(:))
       k_root(:) = k_root_mtc(pft_to_mtc(:))
       k_belowground(:) = k_belowground_mtc(pft_to_mtc(:))
       k_sap(:) = k_sap_mtc(pft_to_mtc(:))
       psi_leaf_min(:) = psi_leaf_min_mtc(pft_to_mtc(:))
       srl(:) = srl_mtc(pft_to_mtc(:))
       r_froot(:) = r_froot_mtc(pft_to_mtc(:))
       psi_root_min(:) = psi_root_min_mtc(pft_to_mtc(:))
       !
       ! Laieff
       !
       crown_to_height(:) = crown_to_height_mtc(pft_to_mtc(:))
       crown_vertohor_dia(:) = crown_vertohor_dia_mtc(pft_to_mtc(:))
       pipe_density(:) = pipe_density_mtc(pft_to_mtc(:))
       tree_ff(:) = tree_ff_mtc(pft_to_mtc(:))
       min_pipe_tune2(:) = min_pipe_tune2_mtc(pft_to_mtc(:))
       ! The dynamic pipe_tune2 tends to be lower than its
       ! static value height growth slowed down too much.
       ! Trying to compensate by changing the diameter-height
       ! relationship. 
       pipe_tune3(:) = pipe_tune3a_mtc(pft_to_mtc(:))
       pipe_tune4(:) = pipe_tune4_mtc(pft_to_mtc(:))
       pipe_k1(:) = pipe_k1_mtc(pft_to_mtc(:)) 
       sla(:) = sla_mtc(pft_to_mtc(:))
       slainit(:) = slainit_mtc(pft_to_mtc(:))
       lai_to_height(:) = lai_to_height_mtc(pft_to_mtc(:)) 

    ENDIF !(ok_sechiba)


    ! For mosses in sechiba
    moss_frac(:) = moss_frac_mtc(pft_to_mtc(:))

    ! 1.3 For BVOC parameters

    IF (ok_bvoc) THEN
       !
       ! Biogenic Volatile Organic Compounds
       !
       em_factor_isoprene(:) = em_factor_isoprene_mtc(pft_to_mtc(:))
       em_factor_monoterpene(:) = em_factor_monoterpene_mtc(pft_to_mtc(:))
       LDF_mono = LDF_mono_mtc 
       LDF_sesq = LDF_sesq_mtc 
       LDF_meth = LDF_meth_mtc 
       LDF_acet = LDF_acet_mtc 

       em_factor_apinene(:) = em_factor_apinene_mtc(pft_to_mtc(:))
       em_factor_bpinene(:) = em_factor_bpinene_mtc(pft_to_mtc(:))
       em_factor_limonene(:) = em_factor_limonene_mtc(pft_to_mtc(:))
       em_factor_myrcene(:) = em_factor_myrcene_mtc(pft_to_mtc(:))
       em_factor_sabinene(:) = em_factor_sabinene_mtc(pft_to_mtc(:))
       em_factor_camphene(:) = em_factor_camphene_mtc(pft_to_mtc(:))
       em_factor_3carene(:) = em_factor_3carene_mtc(pft_to_mtc(:))
       em_factor_tbocimene(:) = em_factor_tbocimene_mtc(pft_to_mtc(:))
       em_factor_othermonot(:) = em_factor_othermonot_mtc(pft_to_mtc(:))
       em_factor_sesquiterp(:) = em_factor_sesquiterp_mtc(pft_to_mtc(:))

       beta_mono = beta_mono_mtc
       beta_sesq = beta_sesq_mtc
       beta_meth = beta_meth_mtc
       beta_acet = beta_acet_mtc
       beta_oxyVOC = beta_oxyVOC_mtc

       em_factor_ORVOC(:) = em_factor_ORVOC_mtc(pft_to_mtc(:)) 
       em_factor_OVOC(:) = em_factor_OVOC_mtc(pft_to_mtc(:))
       em_factor_MBO(:) = em_factor_MBO_mtc(pft_to_mtc(:))
       em_factor_methanol(:) = em_factor_methanol_mtc(pft_to_mtc(:))
       em_factor_acetone(:) = em_factor_acetone_mtc(pft_to_mtc(:)) 
       em_factor_acetal(:) = em_factor_acetal_mtc(pft_to_mtc(:))
       em_factor_formal(:) = em_factor_formal_mtc(pft_to_mtc(:))
       em_factor_acetic(:) = em_factor_acetic_mtc(pft_to_mtc(:))
       em_factor_formic(:) = em_factor_formic_mtc(pft_to_mtc(:))
       em_factor_DMS(:) = em_factor_DMS_mtc(pft_to_mtc(:))
       em_factor_H2S(:) = em_factor_H2S_mtc(pft_to_mtc(:))
       em_factor_no_wet(:) = em_factor_no_wet_mtc(pft_to_mtc(:))
       em_factor_no_dry(:) = em_factor_no_dry_mtc(pft_to_mtc(:))
       Larch(:) = Larch_mtc(pft_to_mtc(:)) 
       !-
    ENDIF !(ok_bvoc)

    ! 1.4 For stomate parameters
    
    IF (ok_stomate) THEN

       !
       ! Vegetation structure - stomate
       !
       availability_fact(:) = availability_fact_mtc(pft_to_mtc(:))
       !
       ! Respiration - stomate
       !
       frac_growthresp(:) = frac_growthresp_mtc(pft_to_mtc(:))
       IF (maint_resp_control == 'cn') THEN
          coeff_maint_init(:) = coeff_maint_init_cn_mtc(pft_to_mtc(:))
       ELSE
          coeff_maint_init(:) = coeff_maint_init_mtc(pft_to_mtc(:))
       END IF
       tref_maint_resp(:) = tref_maint_resp_mtc(pft_to_mtc(:))
       tmin_maint_resp(:) = tmin_maint_resp_mtc(pft_to_mtc(:))
       e0_maint_resp(:) = e0_maint_resp_mtc(pft_to_mtc(:))
       !
       ! Allocation - stomate
       ! 
       tref_labile(:) = tref_labile_mtc(pft_to_mtc(:))
       tmin_labile(:) = tmin_labile_mtc(pft_to_mtc(:))
       e0_labile(:) = e0_labile_mtc(pft_to_mtc(:))
       always_labile(:) = always_labile_mtc(pft_to_mtc(:))

       !
       ! Spitfire - stomate
       !
       bulkdensity_deadfuel(:) = bulkdensity_deadfuel_mtc(pft_to_mtc(:))
       bulkdensity_livefuel(:) = bulkdensity_livefuel_mtc(pft_to_mtc(:))
       f_sh(:) = f_sh_mtc(pft_to_mtc(:))
       BTpar1(:) = BTpar1_mtc(pft_to_mtc(:))
       BTpar2(:) = BTpar2_mtc(pft_to_mtc(:))
       r_ck(:) = r_ck_mtc(pft_to_mtc(:))
       p_ck(:) = p_ck_mtc(pft_to_mtc(:))
       ef_CO2(:) = ef_CO2_mtc(pft_to_mtc(:))
       ef_CO(:) = ef_CO_mtc(pft_to_mtc(:))
       ef_CH4(:) = ef_CH4_mtc(pft_to_mtc(:))
       ef_VOC(:) = ef_VOC_mtc(pft_to_mtc(:))
       ef_TPM(:) = ef_TPM_mtc(pft_to_mtc(:))
       ef_NOx(:) = ef_NOx_mtc(pft_to_mtc(:))
       me(:) = me_mtc(pft_to_mtc(:))
       fire_max_cf_100hr(:) = fire_max_cf_100hr_mtc(pft_to_mtc(:))
       fire_max_cf_1000hr(:) = fire_max_cf_1000hr_mtc(pft_to_mtc(:))
       burnable(:) = burnable_mtc(pft_to_mtc(:))
       fdi_sf_pft(:) = fdi_sf_pft_mtc(pft_to_mtc(:))

       !
       ! Flux - LUC
       !
       coeff_lcchange_s(:) = coeff_lcchange_s_mtc(pft_to_mtc(:))
       coeff_lcchange_m(:) = coeff_lcchange_m_mtc(pft_to_mtc(:))
       coeff_lcchange_l(:) = coeff_lcchange_l_mtc(pft_to_mtc(:))

       !
       ! Phenology
       !
       !
       ! 1. Stomate
       !
       lai_max_to_happy(:) = lai_max_to_happy_mtc(pft_to_mtc(:))  
       lai_max(:) = lai_max_mtc(pft_to_mtc(:))
       pheno_type(:) = pheno_type_mtc(pft_to_mtc(:))
       !
       ! 2. Leaf Onset
       !
       force_pheno(:) = force_pheno_mtc(pft_to_mtc(:))
       pheno_gdd_crit_c(:) = pheno_gdd_crit_c_mtc(pft_to_mtc(:))
       pheno_gdd_crit_b(:) = pheno_gdd_crit_b_mtc(pft_to_mtc(:))         
       pheno_gdd_crit_a(:) = pheno_gdd_crit_a_mtc(pft_to_mtc(:))
       pheno_moigdd_t_crit(:) = pheno_moigdd_t_crit_mtc(pft_to_mtc(:))
       ngd_crit(:) =  ngd_crit_mtc(pft_to_mtc(:))
       ncdgdd_temp(:) = ncdgdd_temp_mtc(pft_to_mtc(:)) 
       hum_frac(:) = hum_frac_mtc(pft_to_mtc(:))
       hum_min_time(:) = hum_min_time_mtc(pft_to_mtc(:))
       relsoilmoist_always(:) = relsoilmoist_always_mtc(pft_to_mtc(:))
       t_always(:) = t_always_mtc(pft_to_mtc(:))
       harvest_time_a(:) = harvest_time_a_mtc(pft_to_mtc(:))
       harvest_time_b(:) = harvest_time_b_mtc(pft_to_mtc(:))
       leaf_age_crit_tref(:) = leaf_age_crit_tref_mtc(pft_to_mtc(:))
       leaf_age_crit_coeff1(:) = leaf_age_crit_coeff1_mtc(pft_to_mtc(:))
       leaf_age_crit_coeff2(:) = leaf_age_crit_coeff2_mtc(pft_to_mtc(:))
       leaf_age_crit_coeff3(:) = leaf_age_crit_coeff3_mtc(pft_to_mtc(:))
       longevity_fruit(:) = longevity_fruit_mtc(pft_to_mtc(:))
       longevity_leaf(:) = longevity_leaf_mtc(pft_to_mtc(:))
       ecureuil(:) = ecureuil_mtc(pft_to_mtc(:))
       alloc_min(:) = alloc_min_mtc(pft_to_mtc(:))
       alloc_max(:) = alloc_max_mtc(pft_to_mtc(:))
       demi_alloc(:) = demi_alloc_mtc(pft_to_mtc(:))
      
       !
       ! 3. Senescence
       !
       leaffall(:) = leaffall_mtc(pft_to_mtc(:))
       presenescence_ratio(:) = presenescence_ratio_mtc(pft_to_mtc(:)) 
       senescence_ratio(:) = senescence_ratio_mtc(pft_to_mtc(:))
       senescence_hum(:) = senescence_hum_mtc(pft_to_mtc(:)) 
       nosenescence_hum(:) = nosenescence_hum_mtc(pft_to_mtc(:)) 
       max_turnover_time(:) = max_turnover_time_mtc(pft_to_mtc(:))
       min_turnover_time(:) = min_turnover_time_mtc(pft_to_mtc(:))
       recycle_leaf(:) = recycle_leaf_mtc(pft_to_mtc(:))
       recycle_root(:) = recycle_root_mtc(pft_to_mtc(:))
       min_leaf_age_for_senescence(:) = min_leaf_age_for_senescence_mtc(pft_to_mtc(:))
       senescence_temp_c(:) = senescence_temp_c_mtc(pft_to_mtc(:))
       senescence_temp_b(:) = senescence_temp_b_mtc(pft_to_mtc(:))
       senescence_temp_a(:) = senescence_temp_a_mtc(pft_to_mtc(:))
       gdd_senescence(:) = gdd_senescence_mtc(pft_to_mtc(:))
       pre_gdd(:) = pre_gdd_mtc(pft_to_mtc(:))
       pre_harv(:) = pre_harv_mtc(pft_to_mtc(:))
       gr_gt(:) = gr_gt_mtc(pft_to_mtc(:))
       always_init(:) = always_init_mtc(pft_to_mtc(:))

       !-
       ! 4. N cycle
       !-
       max_soil_n_bnf(:) = max_soil_n_bnf_mtc(pft_to_mtc(:))
       manure_pftweight(:) =  manure_pftweight_mtc(pft_to_mtc(:))       

       !
       ! DGVM
       !
       residence_time(:) = residence_time_mtc(pft_to_mtc(:))
       tmin_crit(:) = tmin_crit_mtc(pft_to_mtc(:))
       tcm_crit(:) = tcm_crit_mtc(pft_to_mtc(:))
       !-
       
       !
       ! Recruitment (stomate)  
       !  
       recruitment_pft(:) = recruitment_pft_mtc(pft_to_mtc(:))  
       recruitment_height(:) = recruitment_height_mtc(pft_to_mtc(:))  
       recruitment_alpha(:) = recruitment_alpha_mtc(pft_to_mtc(:))  
       recruitment_beta(:) = recruitment_beta_mtc(pft_to_mtc(:))  
       
       !
       ! Mortality - stomate_kill
       !
       beetle_pft(:) = beetle_pft_mtc(pft_to_mtc(:))
       death_distribution_factor(:) = death_distribution_factor_mtc(pft_to_mtc(:))
       npp_reset_value(:) = npp_reset_value_mtc(pft_to_mtc(:))
       ndying_year(:) = ndying_year_mtc(pft_to_mtc(:))
       !
       ! Bark beetle module (stomate)
       !
       S_activity(:) = S_activity_mtc(pft_to_mtc(:))
       act_limit(:) = act_limit_mtc(pft_to_mtc(:))
       S_competition(:) = S_competition_mtc(pft_to_mtc(:))
       S_susceptibility(:) = S_susceptibility_mtc(pft_to_mtc(:))
       i_rd_susceptibility(:) = i_rd_susceptibility_mtc(pft_to_mtc(:))
       S_share(:) = S_share_mtc(pft_to_mtc(:))
       SH_limit(:) = SH_limit_mtc(pft_to_mtc(:))
       S_defence(:) = S_defence_mtc(pft_to_mtc(:))
       PWS_limit(:) = PWS_limit_mtc(pft_to_mtc(:))
       beetle_generation_a(:) = beetle_generation_a_mtc(pft_to_mtc(:))
       beetle_generation_b(:) = beetle_generation_b_mtc(pft_to_mtc(:))
       beetle_generation_c(:) = beetle_generation_c_mtc(pft_to_mtc(:))
       min_temp_beetle(:) = min_temp_beetle_mtc(pft_to_mtc(:))
       max_temp_beetle(:) = max_temp_beetle_mtc(pft_to_mtc(:))
       opt_temp_beetle(:) = opt_temp_beetle_mtc(pft_to_mtc(:))
       eff_temp_beetle_a(:) = eff_temp_beetle_a_mtc(pft_to_mtc(:))
       eff_temp_beetle_b(:) = eff_temp_beetle_b_mtc(pft_to_mtc(:))
       eff_temp_beetle_c(:) = eff_temp_beetle_c_mtc(pft_to_mtc(:))
       eff_temp_beetle_d(:) = eff_temp_beetle_d_mtc(pft_to_mtc(:))
       diapause_thres_daylength(:) = diapause_thres_daylength_mtc(pft_to_mtc(:))
       S_massattack(:) = S_massattack_mtc(pft_to_mtc(:))
       ! Scale dependent bark beetle parameters

       SELECT CASE (forcing_resolution)
       CASE('2degrees_6hourly','1degree_daily','0.5degree_6hourly','LMDzOR_zoomed','LMDzOR')
          ! Parameters for forcings with a coarse spatial resolution
          BP_limit(:) = BP_limit_degree_mtc(pft_to_mtc(:))
          max_Nwood(:) = max_Nwood_degree_mtc(pft_to_mtc(:))
          RDi_limit(:) = RDi_limit_degree_mtc(pft_to_mtc(:))
       CASE('8km_hourly')
          ! Parameters for SAFRAN
          BP_limit(:) = BP_limit_km2_mtc(pft_to_mtc(:))
          max_Nwood(:) = max_Nwood_km2_mtc(pft_to_mtc(:))
          RDi_limit(:) = RDi_limit_km2_mtc(pft_to_mtc(:))
          !+++DELETE+++
          ! Delete once the parameters have been set
          CALL ipslerr_p(3,'Remember to specify the resolution of the forcing', &
               'Cannot select the parameters to calculate pest damage','','')
          !++++++++++++
       CASE('site')
          ! Parameters for FLUXNET 
          BP_limit(:) = BP_limit_ha_mtc(pft_to_mtc(:))
          max_Nwood(:) = max_Nwood_ha_mtc(pft_to_mtc(:))
          RDi_limit(:) = RDi_limit_ha_mtc(pft_to_mtc(:))
          !+++DELETE+++
          ! Delete once the parameters have been set
          CALL ipslerr_p(3,'Remember to specify the resolution of the forcing', &
               'Cannot select the parameters to calculate pest damage','','')
          !++++++++++++
       CASE DEFAULT
          CALL ipslerr_p(3,'Remember to specify the resolution of the forcing', &
               'Cannot select the parameters to calculate pest damage','','')
       ENDSELECT
       
       !
       ! Windfall - stomate_windthrow
       !
       IF (ok_windthrow) THEN
          streamlining_c_leaf(:) = streamlining_c_leaf_mtc(pft_to_mtc(:))
          streamlining_c_leafless(:) = streamlining_c_leafless_mtc(pft_to_mtc(:))
          streamlining_n_leaf(:) = streamlining_n_leaf_mtc(pft_to_mtc(:))
          streamlining_n_leafless(:) = streamlining_n_leafless_mtc(pft_to_mtc(:))
          modulus_rupture(:) = modulus_rupture_mtc(pft_to_mtc(:))
          f_knot(:) = f_knot_mtc(pft_to_mtc(:))
          overturning_free_draining_shallow(:) = overturning_free_draining_shallow_mtc(pft_to_mtc(:))
          overturning_free_draining_shallow_leafless(:) = overturning_free_draining_shallow_leafless_mtc(pft_to_mtc(:))
          overturning_free_draining_deep(:) = overturning_free_draining_deep_mtc(pft_to_mtc(:))
          overturning_free_draining_deep_leafless(:) = overturning_free_draining_deep_leafless_mtc(pft_to_mtc(:))
          overturning_free_draining_average(:) = overturning_free_draining_average_mtc(pft_to_mtc(:))
          overturning_free_draining_average_leafless(:) = overturning_free_draining_average_leafless_mtc(pft_to_mtc(:))
          overturning_gleyed_shallow(:) = overturning_gleyed_shallow_mtc(pft_to_mtc(:))
          overturning_gleyed_shallow_leafless(:) = overturning_gleyed_shallow_leafless_mtc(pft_to_mtc(:))
          overturning_gleyed_deep(:) = overturning_gleyed_deep_mtc(pft_to_mtc(:))
          overturning_gleyed_deep_leafless(:) = overturning_gleyed_deep_leafless_mtc(pft_to_mtc(:))
          overturning_gleyed_average(:) = overturning_gleyed_average_mtc(pft_to_mtc(:))
          overturning_gleyed_average_leafless(:) = overturning_gleyed_average_leafless_mtc(pft_to_mtc(:))
          overturning_peaty_shallow(:) = overturning_peaty_shallow_mtc(pft_to_mtc(:))
          overturning_peaty_shallow_leafless(:) = overturning_peaty_shallow_leafless_mtc(pft_to_mtc(:))
          overturning_peaty_deep(:) = overturning_peaty_deep_mtc(pft_to_mtc(:))
          overturning_peaty_deep_leafless(:) = overturning_peaty_deep_leafless_mtc(pft_to_mtc(:))
          overturning_peaty_average(:) = overturning_peaty_average_mtc(pft_to_mtc(:))
          overturning_peaty_average_leafless(:) = overturning_peaty_average_leafless_mtc(pft_to_mtc(:))
          overturning_peat_shallow(:) = overturning_peat_shallow_mtc(pft_to_mtc(:))
          overturning_peat_shallow_leafless(:) = overturning_peat_shallow_leafless_mtc(pft_to_mtc(:))
          overturning_peat_deep(:) = overturning_peat_deep_mtc(pft_to_mtc(:))
          overturning_peat_deep_leafless(:) = overturning_peat_deep_leafless_mtc(pft_to_mtc(:))
          overturning_peat_average(:) = overturning_peat_average_mtc(pft_to_mtc(:))
          overturning_peat_average_leafless(:) = overturning_peat_average_leafless_mtc(pft_to_mtc(:))
          max_damage_further(:) = max_damage_further_mtc(pft_to_mtc(:))
          max_damage_closer(:) = max_damage_closer_mtc(pft_to_mtc(:))
          sfactor_further(:) = sfactor_further_mtc(pft_to_mtc(:))
          sfactor_closer(:) = sfactor_closer_mtc(pft_to_mtc(:))
          green_density(:) = green_density_mtc(pft_to_mtc(:))
       END IF


       !
       ! SOM decomposition (stomate)
       !
       
       LC_leaf(:)       = LC_leaf_mtc(pft_to_mtc(:))
       LC_sapabove(:)   = LC_sapabove_mtc(pft_to_mtc(:))
       LC_sapbelow(:)   = LC_sapbelow_mtc(pft_to_mtc(:))
       LC_heartabove(:) = LC_heartabove_mtc(pft_to_mtc(:))
       LC_heartbelow(:) = LC_heartbelow_mtc(pft_to_mtc(:))
       LC_fruit(:)      = LC_fruit_mtc(pft_to_mtc(:))
       LC_root(:)       = LC_root_mtc(pft_to_mtc(:))
       LC_carbres(:)    = LC_carbres_mtc(pft_to_mtc(:))
       LC_labile(:)     = LC_labile_mtc(pft_to_mtc(:))
 
       decomp_factor(:) = decomp_factor_mtc(pft_to_mtc(:))
 
       !
       ! Stand structure
       !
       mass_ratio_heart_sap(:) = mass_ratio_heart_sap_mtc(pft_to_mtc(:))
       canopy_cover = canopy_cover_mtc(pft_to_mtc(:))
       nmaxplants(:) = nmaxplants_mtc(pft_to_mtc(:))
       height_init(:) = height_init_mtc(pft_to_mtc(:))
       biomass_init(:) = biomass_init_mtc(pft_to_mtc(:))
       qmd_init(:) = qmd_init_mtc(pft_to_mtc(:))
       deleuze_a(:) = deleuze_a_mtc(pft_to_mtc(:))
       deleuze_b(:) = deleuze_b_mtc(pft_to_mtc(:))
       deleuze_p_all(:) = deleuze_p_all_mtc(pft_to_mtc(:))
       deleuze_p_coppice(:) = deleuze_p_coppice_mtc(pft_to_mtc(:))
       deleuze_power_a(:) = deleuze_power_a_mtc(pft_to_mtc(:))
       beta_self_thinning(:) = beta_self_thinning_mtc(pft_to_mtc(:))
       ref_alpha_self_thin(:) = ref_alpha_self_thin_mtc(pft_to_mtc(:))
       init_pre_indust_ref_gpp(:) = init_pre_indust_ref_gpp_mtc(pft_to_mtc(:))
       
       fuelwood_diameter(:) = fuelwood_diameter_mtc(pft_to_mtc(:))
       coppice_kill_be_wood(:) = coppice_kill_be_wood_mtc(pft_to_mtc(:))
       thinstrat(:) = thinstrat_mtc(pft_to_mtc(:))
       taumin(:) = taumin_mtc(pft_to_mtc(:))
       taumax(:) = taumax_mtc(pft_to_mtc(:))            
       final_dia(:,ifm_none) = final_dia_none_mtc(pft_to_mtc(:))
       final_dia(:,ifm_thin) = final_dia_thin_mtc(pft_to_mtc(:))
       final_dia(:,ifm_uneven)  = final_dia_uneven_mtc(pft_to_mtc(:))
       final_dia(:,ifm_cop)  = final_dia_cop_mtc(pft_to_mtc(:))
       final_dia(:,ifm_src)  = final_dia_src_mtc(pft_to_mtc(:))
       
       ! Select the right scale-dependent parameters
       SELECT CASE (concept_scale)
          CASE ('hectare')
             largest_tree_dia(:) = largest_tree_dia_ha_mtc(pft_to_mtc(:))
             rotation_ref(:)     = rotation_ref_mtc(pft_to_mtc(:))
             rdi_dia_ref(:) = rdi_dia_ref_ha_mtc(pft_to_mtc(:))
             rdi_start(:) = rdi_start_ha_mtc(pft_to_mtc(:))
             rdi_max(:,ifm_none) = rdi_max_none_ha_mtc(pft_to_mtc(:))
             rdi_max(:,ifm_thin) = rdi_max_thin_ha_mtc(pft_to_mtc(:))
             rdi_max(:,ifm_uneven)  = rdi_max_uneven_ha_mtc(pft_to_mtc(:))
             rdi_max(:,ifm_cop)  = rdi_max_cop_ha_mtc(pft_to_mtc(:))
             rdi_max(:,ifm_src)  = rdi_max_src_ha_mtc(pft_to_mtc(:))
             dens_target(:,ifm_none) = dens_target_none_ha_mtc(pft_to_mtc(:))
             dens_target(:,ifm_thin) = dens_target_thin_ha_mtc(pft_to_mtc(:))
             dens_target(:,ifm_uneven)  = dens_target_uneven_ha_mtc(pft_to_mtc(:))
             dens_target(:,ifm_cop)  = dens_target_cop_ha_mtc(pft_to_mtc(:))
             dens_target(:,ifm_src)  = dens_target_src_ha_mtc(pft_to_mtc(:))             
             dens_target(:,0) = 0.0_r_std
             branch_harvest(:,ifm_none) = branch_harvest_none_ha_mtc(pft_to_mtc(:))
             branch_harvest(:,ifm_thin) = branch_harvest_thin_ha_mtc(pft_to_mtc(:))
             branch_harvest(:,ifm_uneven)  = branch_harvest_uneven_ha_mtc(pft_to_mtc(:))
             branch_harvest(:,ifm_cop)  = branch_harvest_cop_ha_mtc(pft_to_mtc(:))
             branch_harvest(:,ifm_src)  = branch_harvest_src_ha_mtc(pft_to_mtc(:))             
             branch_harvest(:,0) = 0.0_r_std
             clearcut_residue(:,ifm_none) = clearcut_residue_none_ha_mtc(pft_to_mtc(:))
             clearcut_residue(:,ifm_thin) = clearcut_residue_thin_ha_mtc(pft_to_mtc(:))
             clearcut_residue(:,ifm_uneven)  = clearcut_residue_uneven_ha_mtc(pft_to_mtc(:))
             clearcut_residue(:,ifm_cop)  = clearcut_residue_cop_ha_mtc(pft_to_mtc(:))
             clearcut_residue(:,ifm_src)  = clearcut_residue_src_ha_mtc(pft_to_mtc(:))             
             clearcut_residue(:,0) = 0.0_r_std
             salvage_dist(:,ifm_none) = salvage_dist_none_ha_mtc(pft_to_mtc(:))
             salvage_dist(:,ifm_thin) = salvage_dist_thin_ha_mtc(pft_to_mtc(:))
             salvage_dist(:,ifm_uneven)  = salvage_dist_uneven_ha_mtc(pft_to_mtc(:))
             salvage_dist(:,ifm_cop)  = salvage_dist_cop_ha_mtc(pft_to_mtc(:))
             salvage_dist(:,ifm_src)  = salvage_dist_src_ha_mtc(pft_to_mtc(:))             
             salvage_dist(:,0) = 0.0_r_std
             dia_boost_end(:) = dia_boost_end_ha_mtc(pft_to_mtc(:))
             thin_intensity(:) = thin_intensity_ha_mtc(pft_to_mtc(:))
             ! HACK for soil carbon - should be removed after proper tuning
             ! Tuning the passive pool is convenient because it has little
             ! impact on the N-cycle but it scientifically difficult to
             ! justify why the passive pool in a CENTURY model should differ
             ! between PFTs. Litter quality and the turnover of the surface
             ! and more active litter and soil pools makes more sense but
             ! will come with a strong impact on the N-cycle.
             som_turn_ipassive(:) = hack_som_turn_ipassive_mtc(pft_to_mtc(:))
          CASE ('europe')
             largest_tree_dia(:) = largest_tree_dia_eu_mtc(pft_to_mtc(:))
             rotation_ref(:)     = rotation_ref_mtc(pft_to_mtc(:))
             rdi_dia_ref(:) = rdi_dia_ref_eu_mtc(pft_to_mtc(:))
             rdi_start(:) = rdi_start_eu_mtc(pft_to_mtc(:))
             rdi_max(:,ifm_none) = rdi_max_none_eu_mtc(pft_to_mtc(:))
             rdi_max(:,ifm_thin) = rdi_max_thin_eu_mtc(pft_to_mtc(:))
             rdi_max(:,ifm_uneven)  = rdi_max_uneven_eu_mtc(pft_to_mtc(:))
             rdi_max(:,ifm_cop)  = rdi_max_cop_eu_mtc(pft_to_mtc(:))
             rdi_max(:,ifm_src)  = rdi_max_src_eu_mtc(pft_to_mtc(:))
             dens_target(:,ifm_none) = dens_target_none_eu_mtc(pft_to_mtc(:))
             dens_target(:,ifm_thin) = dens_target_thin_eu_mtc(pft_to_mtc(:))
             dens_target(:,ifm_uneven)  = dens_target_uneven_eu_mtc(pft_to_mtc(:))
             dens_target(:,ifm_cop)  = dens_target_cop_eu_mtc(pft_to_mtc(:))
             dens_target(:,ifm_src)  = dens_target_src_eu_mtc(pft_to_mtc(:))
             dens_target(:,0) = 0.0_r_std
             dia_boost_end(:) = dia_boost_end_eu_mtc(pft_to_mtc(:))
             thin_intensity(:) = thin_intensity_eu_mtc(pft_to_mtc(:))
             branch_harvest(:,ifm_none) = branch_harvest_none_eu_mtc(pft_to_mtc(:))
             branch_harvest(:,ifm_thin) = branch_harvest_thin_eu_mtc(pft_to_mtc(:))
             branch_harvest(:,ifm_uneven)  = branch_harvest_uneven_eu_mtc(pft_to_mtc(:))
             branch_harvest(:,ifm_cop)  = branch_harvest_cop_eu_mtc(pft_to_mtc(:))
             branch_harvest(:,ifm_src)  = branch_harvest_src_eu_mtc(pft_to_mtc(:))             
             branch_harvest(:,0) = 0.0_r_std
             clearcut_residue(:,ifm_none) = clearcut_residue_none_eu_mtc(pft_to_mtc(:))
             clearcut_residue(:,ifm_thin) = clearcut_residue_thin_eu_mtc(pft_to_mtc(:))
             clearcut_residue(:,ifm_uneven)  = clearcut_residue_uneven_eu_mtc(pft_to_mtc(:))
             clearcut_residue(:,ifm_cop)  = clearcut_residue_cop_eu_mtc(pft_to_mtc(:))
             clearcut_residue(:,ifm_src)  = clearcut_residue_src_eu_mtc(pft_to_mtc(:))             
             clearcut_residue(:,0) = 0.0_r_std
             salvage_dist(:,ifm_none) = salvage_dist_none_eu_mtc(pft_to_mtc(:))
             salvage_dist(:,ifm_thin) = salvage_dist_thin_eu_mtc(pft_to_mtc(:))
             salvage_dist(:,ifm_uneven)  = salvage_dist_uneven_eu_mtc(pft_to_mtc(:))
             salvage_dist(:,ifm_cop)  = salvage_dist_cop_eu_mtc(pft_to_mtc(:))
             salvage_dist(:,ifm_src)  = salvage_dist_src_eu_mtc(pft_to_mtc(:))             
             salvage_dist(:,0) = 0.0_r_std
             ! HACK for soil carbon - should be removed after proper tuning
             ! Tuning the passive pool is convenient because it has little
             ! impact on the N-cycle but it scientifically difficult to
             ! justify why the passive pool in a CENTURY model should differ
             ! between PFTs. Litter quality and the turnover of the surface
             ! and more active litter and soil pools makes more sense but
             ! will come with a strong impact on the N-cycle.
             som_turn_ipassive(:) = hack_som_turn_ipassive_mtc(pft_to_mtc(:))
          CASE ('global')
             largest_tree_dia(:) = largest_tree_dia_global_mtc(pft_to_mtc(:))
             rotation_ref(:)     = rotation_ref_mtc(pft_to_mtc(:))
             rdi_dia_ref(:) = rdi_dia_ref_global_mtc(pft_to_mtc(:))
             rdi_start(:) = rdi_start_global_mtc(pft_to_mtc(:))
             rdi_max(:,ifm_none) = rdi_max_none_global_mtc(pft_to_mtc(:))
             rdi_max(:,ifm_thin) = rdi_max_thin_global_mtc(pft_to_mtc(:))
             rdi_max(:,ifm_uneven)  = rdi_max_uneven_global_mtc(pft_to_mtc(:))
             rdi_max(:,ifm_cop)  = rdi_max_cop_global_mtc(pft_to_mtc(:))
             rdi_max(:,ifm_src)  = rdi_max_src_global_mtc(pft_to_mtc(:))
             dens_target(:,ifm_none) = dens_target_none_global_mtc(pft_to_mtc(:))
             dens_target(:,ifm_thin) = dens_target_thin_global_mtc(pft_to_mtc(:))
             dens_target(:,ifm_uneven)  = dens_target_uneven_global_mtc(pft_to_mtc(:))
             dens_target(:,ifm_cop)  = dens_target_cop_global_mtc(pft_to_mtc(:))
             dens_target(:,ifm_src)  = dens_target_src_global_mtc(pft_to_mtc(:))             
             dens_target(:,0) = 0.0_r_std
             dia_boost_end(:) = dia_boost_end_global_mtc(pft_to_mtc(:))
             thin_intensity(:) = thin_intensity_global_mtc(pft_to_mtc(:))
             branch_harvest(:,ifm_none) = branch_harvest_none_global_mtc(pft_to_mtc(:))
             branch_harvest(:,ifm_thin) = branch_harvest_thin_global_mtc(pft_to_mtc(:))
             branch_harvest(:,ifm_uneven)  = branch_harvest_uneven_global_mtc(pft_to_mtc(:))
             branch_harvest(:,ifm_cop)  = branch_harvest_cop_global_mtc(pft_to_mtc(:))
             branch_harvest(:,ifm_src)  = branch_harvest_src_global_mtc(pft_to_mtc(:))             
             branch_harvest(:,0) = 0.0_r_std
             clearcut_residue(:,ifm_none) = clearcut_residue_none_global_mtc(pft_to_mtc(:))
             clearcut_residue(:,ifm_thin) = clearcut_residue_thin_global_mtc(pft_to_mtc(:))
             clearcut_residue(:,ifm_uneven)  = clearcut_residue_uneven_global_mtc(pft_to_mtc(:))
             clearcut_residue(:,ifm_cop)  = clearcut_residue_cop_global_mtc(pft_to_mtc(:))
             clearcut_residue(:,ifm_src)  = clearcut_residue_src_global_mtc(pft_to_mtc(:))             
             clearcut_residue(:,0) = 0.0_r_std
             salvage_dist(:,ifm_none) = salvage_dist_none_global_mtc(pft_to_mtc(:))
             salvage_dist(:,ifm_thin) = salvage_dist_thin_global_mtc(pft_to_mtc(:))
             salvage_dist(:,ifm_uneven)  = salvage_dist_uneven_global_mtc(pft_to_mtc(:))
             salvage_dist(:,ifm_cop)  = salvage_dist_cop_global_mtc(pft_to_mtc(:))
             salvage_dist(:,ifm_src)  = salvage_dist_src_global_mtc(pft_to_mtc(:))             
             salvage_dist(:,0) = 0.0_r_std
             ! HACK for soil carbon - should be removed after proper tuning
             ! Tuning the passive pool is convenient because it has little
             ! impact on the N-cycle but it scientifically difficult to
             ! justify why the passive pool in a CENTURY model should differ
             ! between PFTs. Litter quality and the turnover of the surface
             ! and more active litter and soil pools makes more sense but
             ! will come with a strong impact on the N-cycle.
             som_turn_ipassive(:) = som_turn_ipassive_mtc(pft_to_mtc(:))
          CASE DEFAULT
            WRITE(numout,*) 'The conceptual scale of the simulation was set to ', concept_scale
            CALL ipslerr_p(3, 'The value for concept_scale is not known',&
                 'define its value in pft_parameters','','')
       END SELECT

       circ_dist_shape_none(:) = circ_dist_shape_none_mtc(pft_to_mtc(:))  
       circ_dist_shape_thin(:) = circ_dist_shape_thin_mtc(pft_to_mtc(:))
       circ_dist_shape_cop(:) = circ_dist_shape_cop_mtc(pft_to_mtc(:))
       circ_dist_shape_src(:) = circ_dist_shape_src_mtc(pft_to_mtc(:))
       circ_dist_shape_uneven(:) = circ_dist_shape_uneven_mtc(pft_to_mtc(:))
       branch_ratio(:) = branch_ratio_mtc(pft_to_mtc(:))
       coppice_diameter(:) = coppice_diameter_mtc(pft_to_mtc(:))
       shoots_per_stool(:) = shoots_per_stool_mtc(pft_to_mtc(:))
       src_rot_length(:) = src_rot_length_mtc(pft_to_mtc(:))
       src_nrots(:) = src_nrots_mtc(pft_to_mtc(:))

       forest_managed_forced(:) = forest_managed_forced_mtc(pft_to_mtc(:))
       m_dv(:) = m_dv_mtc(pft_to_mtc(:))
       fruit_alloc(:) = fruit_alloc_mtc(pft_to_mtc(:))
 
       labile_reserve(:) = labile_reserve_mtc(pft_to_mtc(:))
       evergreen_reserve(:) = evergreen_reserve_mtc(pft_to_mtc(:))
       deciduous_reserve(:) = deciduous_reserve_mtc(pft_to_mtc(:))
       senescense_reserve(:) = senescense_reserve_mtc(pft_to_mtc(:))
       root_reserve(:) = root_reserve_mtc(pft_to_mtc(:))

       fcn_wood(:) = fcn_wood_mtc(pft_to_mtc(:))
       fcn_root(:) = fcn_root_mtc(pft_to_mtc(:))

       ! 
       ! Cropland management 
       ! 
       harvest_ratio(:) = harvest_ratio_mtc(pft_to_mtc(:))

    ENDIF !(ok_stomate)
    
    !! Following parameters are used with and without ok_stomate
    nue_opt(:) = nue_opt_mtc(pft_to_mtc(:))
    vmax_uptake(:,iammonium) = vmax_uptake_nh4_mtc(pft_to_mtc(:))
    vmax_uptake(:,initrate) = vmax_uptake_no3_mtc(pft_to_mtc(:))
    cn_leaf_min(:) = cn_leaf_min_mtc(pft_to_mtc(:))
    cn_leaf_max(:) = cn_leaf_max_mtc(pft_to_mtc(:))
    cn_leaf_init(:) = cn_leaf_init_mtc(pft_to_mtc(:))
    ext_coeff_N(:) = ext_coeff_N_mtc(pft_to_mtc(:))
    IF (maint_resp_control == 'cn') THEN
       maint_resp_slope_c(:) = maint_resp_slope_c_cn_mtc(pft_to_mtc(:))
       maint_resp_slope_b(:) = maint_resp_slope_b_cn_mtc(pft_to_mtc(:)) 
    ELSE
       maint_resp_slope_c(:) = maint_resp_slope_c_mtc(pft_to_mtc(:))
       maint_resp_slope_b(:) = maint_resp_slope_b_mtc(pft_to_mtc(:)) 
    END IF
    maint_resp_slope_a(:) = maint_resp_slope_a_mtc(pft_to_mtc(:))
    dia_thresh_inv(:) = dia_thresh_inv_mtc(pft_to_mtc(:))
    longevity_root(:) = longevity_root_mtc(pft_to_mtc(:))
    longevity_sap(:) = longevity_sap_mtc(pft_to_mtc(:))
    k_latosa_max(:) = k_latosa_max_mtc(pft_to_mtc(:))
    k_latosa_min(:) = k_latosa_min_mtc(pft_to_mtc(:))
    senescence_type(:) = senescence_type_mtc(pft_to_mtc(:))

  END SUBROUTINE pft_parameters_init


!! ================================================================================================================================
!! SUBROUTINE   : pft_parameters_alloc
!!
!>\BRIEF         This subroutine allocates memory needed for the PFT parameters 
!! in function  of the flags activated.  
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE pft_parameters_alloc()

    IMPLICIT NONE

    !! 0. Variables and parameters declaration

    !! 0.1 Input variables 

    !! 0.4 Local variables

    LOGICAL :: l_error                             !! Diagnostic boolean for error allocation (true/false) 
    INTEGER :: ier                                 !! Return value for memory allocation (0-N, unitless)

    !_ ================================================================================================================================


    !
    ! 1. Parameters used anytime
    !

    l_error = .FALSE.

    ALLOCATE(pft_to_mtc(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for pft_to_mtc. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(PFT_name(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for PFT_name. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(height_presc(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for height_presc. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(z0_over_height(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for z0_over_height. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(ratio_z0m_z0h(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for ratio_z0m_z0h. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(is_tree(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for is_tree. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(natural(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for natural. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(is_c4(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for is_c4. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(humcste(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for humcste. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(max_root_depth(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for max_root_depth. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(downregulation_co2_coeff(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for downregulation_co2_coeff. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(E_KmC(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for E_KmC. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(E_KmO(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for E_KmO. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(E_Sco(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for E_Sco. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(E_gamma_star(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for E_gamma_star. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(E_vcmax(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for E_Vcmax. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(E_Jmax(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for E_Jmax. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(aSV(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for aSV. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(bSV(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for bSV. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(tphoto_min(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for tphoto_min. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(tphoto_max(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for tphoto_max. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(aSJ(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for aSJ. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(bSJ(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for bSJ. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(D_Vcmax(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for D_Vcmax. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(D_Jmax(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for D_Jmax. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(E_gm(nvm),stat=ier) 
    l_error = l_error .OR. (ier /= 0) 
    IF (l_error) THEN 
       WRITE(numout,*) ' Memory allocation error for E_gm. We stop. We need nvm words = ',nvm 
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','') 
    END IF
    
    ALLOCATE(S_gm(nvm),stat=ier) 
    l_error = l_error .OR. (ier /= 0) 
    IF (l_error) THEN 
       WRITE(numout,*) ' Memory allocation error for S_gm. We stop. We need nvm words = ',nvm 
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','') 
    END IF
    
    ALLOCATE(D_gm(nvm),stat=ier) 
    l_error = l_error .OR. (ier /= 0) 
    IF (l_error) THEN 
       WRITE(numout,*) ' Memory allocation error for D_gm. We stop. We need nvm words = ',nvm 
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','') 
    END IF
    
    ALLOCATE(E_Rd(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for E_Rd. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(Vcmax25(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for Vcmax25. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(arJV(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for arJV. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(brJV(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for brJV. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(KmC25(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for KmC25. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(KmO25(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for KmO25. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(Sco25(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for Sco25. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF
    
    ALLOCATE(gm25(nvm),stat=ier) 
    l_error = l_error .OR. (ier /= 0) 
    IF (l_error) THEN 
       WRITE(numout,*) ' Memory allocation error for gm25. We stop. We need nvm words = ',nvm 
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','') 
    END IF

    ALLOCATE(gamma_star25(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for gamma_star25. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(a1(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for a1. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(b1(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for b1. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(a_fpsi(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for a_fpsi. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(sf(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for sf. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(leaf_wat_store_fac(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for leaf_wat_store_fac. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(wood_wat_store_fac(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for wood_wat_store_fac. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(Rroot_sup(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for Rroot_sup. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(Rroot_inf(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for Rroot_inf. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(Rleaf(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for Rleaf. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(Rxylem(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for Rxylem. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(Rroot_soil(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for Rroot_soil. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(Rsto_leaf(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for Rsto_leaf. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(Rsto_wood(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for Rsto_wood. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(psi_ref_g0(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for psi_ref_g0. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(psi_ref_sto_leaf(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for psi_ref_sto_leaf. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(psi_ref_sto_wood(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for psi_ref_sto_wood. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(psi_50_res_leaf(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for psi_50_res_leaf. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(psi_50_res_wood(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for psi_50_res_wood. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(psi_50_res_root(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for psi_50_res_root. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(a_leaf(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for a_leaf. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(a_wood(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for a_wood. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(a_root(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for a_root. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(k_leaf_max(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for k_leaf_max. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(sapwood_density(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for sapwood_density. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(ks_max(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for ks_max. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(k_wood_max(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for k_wood_max. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(k_root_max(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for k_root_max. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(lambda_leaf(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for lambda_leaf. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(lambda_wood(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for lambda_wood. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(d_l(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for d_l. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(g0(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for g0. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(h_protons(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for h_protons. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(fpsir(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for fpsir. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(fQ(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for fQ. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(fpseudo(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for fpseudo. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(kp(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for kp. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(alpha(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for alpha. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(gbs(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for gbs. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(theta(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for theta. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(alpha_LL(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for alpha_LL. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(stress_vcmax(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for stress_vcmax. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF
    
    ALLOCATE(stress_gs(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for stress_gs. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF
    
    ALLOCATE(stress_gm(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for stress_gm. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(ext_coeff(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for ext_coeff. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(ext_coeff_vegetfrac(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for ext_coeff_vegetfrac. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(llaimax(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for llaimax. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(llaimin(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for llaimin. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(type_of_lai(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for type_of_lai. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(vcmax_fix(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for vcmax_fix. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(pref_soil_veg(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for pref_soil_veg. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(agec_group(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for agec_group. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF
    
    ALLOCATE(start_index(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for start_index. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF
    
    ALLOCATE(nagec_pft(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for nagec_pft. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(leaf_tab(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for leaf_tab. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(pheno_model(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for pheno_model. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(gdd_threshold(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for gdd_threshold. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(is_deciduous(nvm),stat=ier) 
    l_error = l_error .OR. (ier /= 0) 
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for is_deciduous. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(is_evergreen(nvm),stat=ier) 
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for is_evergreen. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(is_needleleaf(nvm),stat=ier)  
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for is_needleleaf. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(is_tropical(nvm),stat=ier)   
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for is_tropical. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(is_temperate(nvm),stat=ier)   
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for is_temperate. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF
    
    ALLOCATE(is_boreal(nvm),stat=ier)   
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for is_boreal. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    !
    ! 2. Parameters used if ok_sechiba only
    !
    IF ( ok_sechiba ) THEN

       l_error = .FALSE.

       ALLOCATE(rstruct_const(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for rstruct_const. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(kzero(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for kzero. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(alpha_watstress_pft(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for alpha_watstress_pft. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF


       ALLOCATE(rveg_pft(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for rveg_pft. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(wmax_veg(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for wmax_veg. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(throughfall_by_pft(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for throughfall_by_pft. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(rsoil_scale(nvm),stat=ier)
       IF (ier /= 0) THEN
          WRITE(numout,*) ' Memory allocation error for rsoil_scale. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(snowa_aged_vis(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for snowa_aged_vis. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(snowa_aged_nir(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for snowa_aged_nir. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(snowa_dec_vis(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for snowa_dec_vis. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(snowa_dec_nir(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for snowa_dec_nir. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(alb_leaf_vis(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for alb_leaf_vis. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(alb_leaf_nir(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for alb_leaf_nir. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(leaf_ssa(nvm,n_spectralbands),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for leaf_ssa. We stop. We need nvm*n_spectralbands words = ',&
               nvm*n_spectralbands
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(leaf_psd(nvm,n_spectralbands),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for leaf_psd. We stop. We need nvm*n_spectralbands words = ',&
               nvm*n_spectralbands
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(bgd_reflectance(nvm,n_spectralbands),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for bgd_reflectance. We need nvm*n_spectralbands words = ',&
               nvm*n_spectralbands
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(leaf_to_shoot_clumping(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for leaf_to_shoot_clumping. We need nvm words = ',&
               nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(lai_correction_factor(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for lai_correction_factor. We need nvm words = ',&
               nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(min_level_sep(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for min_level_sep. We need nvm words = ',&
               nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(lai_top(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for lai_top. We need nvm words = ',&
               nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       !
       ! Hydraulic architecture
       !       

         ALLOCATE(k_root(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for k_root. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(k_belowground(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for k_belowground. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(k_sap(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for k_sap. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(psi_leaf_min(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for psi_leaf_min. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(srl(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for srl. We stop.We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(r_froot(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for r_froot. We stop.We need nvmwords = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(psi_root_min(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for psi_root_min. We stop.We need nvmwords = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(crown_to_height(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for crown_to_height. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(crown_vertohor_dia(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for crown_vertohor_dia. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(pipe_density(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for pipe_density. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(tree_ff(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for tree_ff. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(min_pipe_tune2(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for min_pipe_tune2. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(pipe_tune3(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for pipe_tune3. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

        ALLOCATE(pipe_tune4(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for pipe_tune4. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(pipe_k1(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for pipe_k1. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(sla(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for sla. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(slainit(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for slainit. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

        ALLOCATE(lai_to_height(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for lai_to_height. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       IF( ok_bvoc ) THEN

          l_error = .FALSE.

          ALLOCATE(em_factor_isoprene(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_isoprene. We stop. We need nvm words = ',nvm
             CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
          END IF

          ALLOCATE(em_factor_monoterpene(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_monoterpene. We stop. We need nvm words = ',nvm
             CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
          END IF

          ALLOCATE(em_factor_apinene(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_apinene. We stop. We need nvm words = ',nvm
             CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
          END IF

          ALLOCATE(em_factor_bpinene(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_bpinene. We stop. We need nvm words = ',nvm
             CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
          END IF

          ALLOCATE(em_factor_limonene(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_limonene. We stop. We need nvm words = ',nvm
             CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
          END IF

          ALLOCATE(em_factor_myrcene(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_myrcene. We stop. We need nvm words = ',nvm
             CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
          END IF

          ALLOCATE(em_factor_sabinene(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_sabinene. We stop. We need nvm words = ',nvm
             CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
          END IF

          ALLOCATE(em_factor_camphene(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_camphene. We stop. We need nvm words = ',nvm
             CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
          END IF

          ALLOCATE(em_factor_3carene(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_3carene. We stop. We need nvm words = ',nvm
             CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
          END IF

          ALLOCATE(em_factor_tbocimene(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_tbocimene. We stop. We need nvm words = ',nvm
             CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
          END IF

          ALLOCATE(em_factor_othermonot(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_othermonot. We stop. We need nvm words = ',nvm
             CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
          END IF

          ALLOCATE(em_factor_sesquiterp(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_sesquiterp. We stop. We need nvm words = ',nvm
             CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
          END IF


          ALLOCATE(em_factor_ORVOC(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_ORVOC. We stop. We need nvm words = ',nvm
             CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
          END IF

          ALLOCATE(em_factor_OVOC(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0)       
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_OVOC. We stop. We need nvm words = ',nvm
             CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
          END IF

          ALLOCATE(em_factor_MBO(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_MBO. We stop. We need nvm words = ',nvm
             CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
          END IF

          ALLOCATE(em_factor_methanol(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_methanol. We stop. We need nvm words = ',nvm
             CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
          END IF

          ALLOCATE(em_factor_acetone(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_acetone. We stop. We need nvm words = ',nvm
             CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
          END IF

          ALLOCATE(em_factor_acetal(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_acetal. We stop. We need nvm words = ',nvm
             CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
          END IF

          ALLOCATE(em_factor_formal(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_formal. We stop. We need nvm words = ',nvm
             CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
          END IF

          ALLOCATE(em_factor_acetic(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0)       
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_acetic. We stop. We need nvm words = ',nvm
             CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
          END IF

          ALLOCATE(em_factor_formic(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_formic. We stop. We need nvm words = ',nvm
             CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
          END IF

          ALLOCATE(em_factor_DMS(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0)
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_DMS. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

          ALLOCATE(em_factor_H2S(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0)
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_H2S. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

          ALLOCATE(em_factor_no_wet(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0)
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_no_wet. We stop. We need nvm words = ',nvm
             CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
          END IF

          ALLOCATE(em_factor_no_dry(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0)       
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_no_dry. We stop. We need nvm words = ',nvm
             CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
          END IF

          ALLOCATE(Larch(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for Larch. We stop. We need nvm words = ',nvm
             CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
          END IF

       ENDIF ! (ok_bvoc) 

    ENDIF !(ok_sechiba)

    !
    ! 3. Parameters used if ok_stomate only
    !



    IF ( ok_stomate ) THEN

       l_error = .FALSE.

       ALLOCATE(availability_fact(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for availability_fact. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(pheno_gdd_crit_c(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for pheno_gdd_crit_c. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(pheno_gdd_crit_b(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for pheno_gdd_crit_b. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(pheno_gdd_crit_a(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for pheno_gdd_crit_a. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(pheno_gdd_crit(nvm,3),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for pheno_gdd_crit. We stop. We need nvm words = ',nvm*3
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       pheno_gdd_crit(:,:) = zero

       ALLOCATE(pheno_moigdd_t_crit(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for pheno_moigdd_t_crit. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(ngd_crit(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for ngd_crit. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(ncdgdd_temp(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for ncdgdd_temp. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(hum_frac(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for hum_frac. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(hum_min_time(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for hum_min_time. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(relsoilmoist_always(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for reslsoilmoist_always. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(t_always(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for t_always. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(harvest_time_a(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for harvest_time_a. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(harvest_time_b(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for harvest_time_b. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(longevity_leaf(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for longevity_leaf. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(leaf_age_crit_tref(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for leaf_age_crit_tref. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(leaf_age_crit_coeff1(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for leaf_age_crit_coeff1. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(leaf_age_crit_coeff2(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for leaf_age_crit_coeff2. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(leaf_age_crit_coeff3(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for leaf_age_crit_coeff3. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(longevity_fruit(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for longevity_fruit. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(ecureuil(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for ecureuil. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(alloc_min(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for alloc_min. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(alloc_max(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for alloc_max. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(demi_alloc(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for . We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(frac_growthresp(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for frac_growthresp. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(coeff_maint_init(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for coeff_maint_init. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
 
       ALLOCATE(tref_maint_resp(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for tref_maint_resp. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(tmin_maint_resp(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for tmin_maint_resp. We stop. We need nvm  words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(e0_maint_resp(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for e0_maint_resp. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(tref_labile(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for tref_labile. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(tmin_labile(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for tmin_labile. We stop. We need nvm  words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(e0_labile(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for e0_labile. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(always_labile(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for always_labile. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(bulkdensity_deadfuel(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for bulkdensity_deadfuel. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(bulkdensity_livefuel(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for bulkdensity_livefuel. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(f_sh(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for f_sh. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(BTpar1(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for BTpar1. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(BTpar2(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for BTpar2. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(r_ck(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for r_ck. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(p_ck(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for p_ck. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(fdi_sf_pft(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for fdi_sf_pft. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
    
       ALLOCATE(ef_CO2(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for ef_CO2. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(ef_CO(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for ef_CO. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(ef_CH4(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for ef_CH4. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(ef_VOC(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for ef_VOC. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(ef_TPM(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for ef_TPM. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(ef_NOx(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for ef_NOx. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(me(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for me. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(fire_max_cf_100hr(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for fire_max_cf_100hr. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(fire_max_cf_1000hr(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for fire_max_cf_1000hr. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(burnable(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for burnable. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(coeff_lcchange_s(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for coeff_lcchange_s. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(coeff_lcchange_m(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for coeff_lcchange_m. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(coeff_lcchange_l(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for coeff_lcchange_l. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(lai_max_to_happy(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for lai_max_to_happy. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(lai_max(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for lai_max. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(pheno_type(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for pheno_type. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(force_pheno(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for force_pheno. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(leaffall(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for leaffall. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(presenescence_ratio(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for presenescence_ratio. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(senescence_ratio(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for senescence_ratio. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(senescence_hum(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for senescence_hum. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(nosenescence_hum(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for nosenescence_hum. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(max_turnover_time(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for max_turnover_time. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(min_turnover_time(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for min_turnover_time. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(recycle_leaf(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for recycle_leaf. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(recycle_root(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for recycle_root. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(min_leaf_age_for_senescence(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for min_leaf_age_for_senescence. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(senescence_temp_c(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for senescence_temp_c. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(senescence_temp_b(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for senescence_temp_b. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(senescence_temp_a(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for senescence_temp_a. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(senescence_temp(nvm,3),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for senescence_temp. We stop. We need nvm*3 words = ',nvm*3
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       senescence_temp(:,:) = zero

       ALLOCATE(gdd_senescence(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for gdd_senescence. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(pre_gdd(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for pre_gdd. We stop. We need nvm words = ',nvm
             CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(pre_harv(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for pre_harv. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(gr_gt(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for gr_gt. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(always_init(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for always_init. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(max_soil_n_bnf(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for max_soil_n_bnf. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(manure_pftweight(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for manure_pftweight. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(residence_time(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for residence_time. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(tmin_crit(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for tmin_crit. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(tcm_crit(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for tcm_crit. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(lai_initmin(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for . We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(bm_sapl(nvm,nparts,nelements),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for bm_sapl. We stop. We need nvm*nparts*nelements words = ',& 
               &  nvm*nparts*nelements
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(migrate(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for migrate. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(cn_sapl(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for cn_sapl. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       !
       ! SOM decomposition (stomate)
       !
       ALLOCATE(LC(nvm,nparts),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for LC. We stop. We need nvm*nparts words = ',nvm,nparts
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(LC_leaf(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for LC_leaf. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(LC_sapabove(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for LC_sapabove. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(LC_sapbelow(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for LC_sapbelow. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(LC_heartabove(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for LC_heartabove. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(LC_heartbelow(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for LC_heartbelow. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(LC_fruit(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for LC_fruit. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(LC_root(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for LC_root. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(LC_carbres(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for LC_carbres. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(LC_labile(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for LC_labile. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(decomp_factor(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for decomp_factor. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(mass_ratio_heart_sap(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for mass_ratio_heart_sap. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(canopy_cover(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for canopy_cover. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(nmaxplants(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for nmaxplants. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(height_init(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for height_init. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(biomass_init(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
           WRITE(numout,*) ' Memory allocation error for biomass_init. We stop. We need nvm words = ',nvm
           CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       ENDIF
       
       ALLOCATE(qmd_init(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for qmd_init. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(deleuze_a(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for deleuze_a. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(deleuze_b(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for deleuze_b. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(deleuze_p_all(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for deleuze_p_all. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(beta_self_thinning(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for beta_self_thinning. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(ref_alpha_self_thin(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for ref_alpha_self_thin. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(init_pre_indust_ref_gpp(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for init_pre_indust_ref_gpp. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(fuelwood_diameter(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for fuelwood_diameter. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(coppice_kill_be_wood(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for coppice_kill_be_wood. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(deleuze_p_coppice(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for deleuze_p_coppice. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(deleuze_power_a(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for deleuze_power_a. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(m_dv(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for m_dv. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(forest_managed_forced(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for forest_managed_forced. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(thinstrat(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for thinstrat. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(taumin(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for taumin. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(taumax(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for taumax. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(poly_rdi(nvm, nfm_types, 2, 4),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for poly_rdi. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(rdi_dia_ref(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for rdi_dia_ref. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(rdi_start(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for rdi_start. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(rdi_max(nvm, nfm_types),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for rdi_max. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(dens_target(nvm, 0:nfm_types), stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for dens_target. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(final_dia(nvm, nfm_types),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for final_dia. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(dia_boost_end(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for dia_boost_end. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(rdi_ramp_out(nvm),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'pft_parameters_alloc', &
            'Memory allocation error for rdi_ramp_out','','')
       ALLOCATE(forced_thin_intensity(nvm),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'pft_parameters_alloc', &
            'Memory allocation error for forced_thin_intensity','','')
       ALLOCATE(rdi_band_low(nvm),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'pft_parameters_alloc', &
            'Memory allocation error for rdi_band_low','','')
       ALLOCATE(rdi_band_high(nvm),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'pft_parameters_alloc', &
            'Memory allocation error for rdi_band_high','','')
       ALLOCATE(thin_intensity(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for thin_intensity. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(circ_dist_shape_none(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for circ_dist_shape_none. Westop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF


       ALLOCATE(circ_dist_shape_uneven(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for circ_dist_shape_uneven. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF        
       
       ALLOCATE(circ_dist_shape_thin(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for circ_dist_shape_thin. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF  

       ALLOCATE(circ_dist_shape_cop(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for circ_dist_shape_cop. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF 

       ALLOCATE(circ_dist_shape_src(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for circ_dist_shape_src. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF  

       ALLOCATE(largest_tree_dia(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for largest_tree_dia. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(branch_ratio(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for branch_ratio. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(branch_harvest(nvm,0:nfm_types),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for branch_harvest. We stop. We need nvm * nfm_types words = ',nvm, nfm_types
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
      
       ALLOCATE(clearcut_residue(nvm,0:nfm_types),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for clearcut_residue. We stop.We need nvm * nfm_types words = ',nvm, nfm_types
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(salvage_dist(nvm,0:nfm_types),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for salvage_dist. We stop.Weneed nvm * nfm_types words = ',nvm, nfm_types 
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(coppice_diameter(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for coppice_diameter. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(shoots_per_stool(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for shoots_per_stool. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(src_rot_length(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for src_rot_length. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(src_nrots(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for src_nrots. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
      
       ALLOCATE(fruit_alloc(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for fruit_alloc. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(labile_reserve(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for labile_reserve. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(evergreen_reserve(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for evergreen_reserve. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(deciduous_reserve(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for deciudous_reserve. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(senescense_reserve(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for senescense_reserve. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(root_reserve(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for root_reserve. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(fcn_wood(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for fcn_wood. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       ALLOCATE(fcn_root(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for fcn_root. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF
       
       !
       ! RECRUITMENT 
       ! 
       ALLOCATE(recruitment_pft(nvm),stat=ier)     
       l_error = l_error .OR. (ier /= 0)  
       IF (l_error) THEN  
          WRITE(numout,*) ' Memory allocation error for recruitment_pft. We stop. We need nvm words = ',nvm  
          CALL ipslerr_p (3, 'pft_parameters','pft_parameters_alloc','','')  
       END IF

       ALLOCATE(recruitment_height(nvm),stat=ier)     
       l_error = l_error .OR. (ier /= 0)  
       IF (l_error) THEN  
          WRITE(numout,*) ' Memory allocation error for recruitment_height. We stop. We need nvm words = ',nvm  
          CALL ipslerr_p (3, 'pft_parameters','pft_parameters_alloc','','')  
       END IF

       ALLOCATE(recruitment_alpha(nvm),stat=ier)     
       l_error = l_error .OR. (ier /= 0) 
       IF (l_error) THEN  
          WRITE(numout,*) ' Memory allocation error for recruitment_alpha. We stop. We need nvm words = ',nvm  
          CALL ipslerr_p (3, 'pft_parameters','pft_parameters_alloc','','') 
       END IF

       ALLOCATE(recruitment_beta(nvm),stat=ier)     
       l_error = l_error .OR. (ier /= 0)  
       IF (l_error) THEN  
          WRITE(numout,*) ' Memory allocation error for recruitment_beta. We stop. We need nvm words = ',nvm  
          CALL ipslerr_p (3, 'pft_parameters','pft_parameters_alloc','','')  
       END IF

       !
       ! MORTALITY
       !

       ALLOCATE(beetle_pft(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for beetle_pft. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3, 'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(death_distribution_factor(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for death_distribution_factor. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(npp_reset_value(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for npp_reset_value. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(ndying_year(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for ndying_year. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       !BEETLE

       ALLOCATE(S_activity(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for S_activity. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(act_limit(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for act_limit. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(S_competition(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for S_competition. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(RDi_limit(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for rdi_limt. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(S_susceptibility(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for S_susceptibility. We stop.We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(i_rd_susceptibility(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for i_rd_susceptibility. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(S_share(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for S_share. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(SH_limit(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for SH_limit. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(S_defence(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for S_defence. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(PWS_limit(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for PWS_limit. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(max_Nwood(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for max_Nwood. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

      ALLOCATE(beetle_generation_a(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for beetle_generation_a. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

      ALLOCATE(beetle_generation_b(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for beetle_generation_b. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

      ALLOCATE(beetle_generation_c(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for beetle_generation_c. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

      ALLOCATE(min_temp_beetle(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for min_temp_beetle. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

      ALLOCATE(max_temp_beetle(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for max_temp_beetle. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

      ALLOCATE(opt_temp_beetle(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for opt_temp_beetle. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF


      ALLOCATE(eff_temp_beetle_a(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for eff_temp_beetle_a. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

      ALLOCATE(eff_temp_beetle_b(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for eff_temp_beetle_b. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

      ALLOCATE(eff_temp_beetle_c(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for eff_temp_beetle_c. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

      ALLOCATE(eff_temp_beetle_d(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for eff_temp_beetle_d. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

      ALLOCATE(diapause_thres_daylength(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for diapause_thres_daylength. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(S_massattack(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for wgth_sirdi_a. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(BP_limit(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for BP_limit. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       !
       ! WINDTHROW
       !

       ALLOCATE(streamlining_c_leaf(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for streamlining_c_leaf. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(streamlining_c_leafless(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for streamlining_c_leafless. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(streamlining_n_leaf(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for streamlining_n_leaf. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(streamlining_n_leafless(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for streamlining_n_leafless. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(modulus_rupture(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for modulus_rupture. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(f_knot(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for f_knot. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(overturning_free_draining_shallow(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for overturning_free_draining_shallow. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(overturning_free_draining_shallow_leafless(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for overturning_free_draining_shallow_leafless.',&
               ' We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(overturning_free_draining_deep(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for overturning_free_draining_deep. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(overturning_free_draining_deep_leafless(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for overturning_free_draining_deep_leafles. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(overturning_free_draining_average(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for overturning_free_draining_average. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(overturning_free_draining_average_leafless(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for overturning_free_draining_average_leafless. ', &
               ' We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(overturning_gleyed_shallow(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for overturning_gleyed_shallow. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(overturning_gleyed_shallow_leafless(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for overturning_gleyed_shallow_leafless. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(overturning_gleyed_deep(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for overturning_gleyed_deep. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(overturning_gleyed_deep_leafless(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for overturning_gleyed_deep_leafless. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(overturning_gleyed_average(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for overturning_gleyed_average. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(overturning_gleyed_average_leafless(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for overturning_gleyed_average_leafless. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(overturning_peaty_shallow(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for overturning_peaty_shallow. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(overturning_peaty_shallow_leafless(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for overturning_peaty_shallow_leafless. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(overturning_peaty_deep(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for overturning_peaty_deep. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(overturning_peaty_deep_leafless(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for overturning_peaty_deep_leafless. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(overturning_peaty_average(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for overturning_peaty_average. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(overturning_peaty_average_leafless(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for overturning_peaty_average_leafless. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(overturning_peat_shallow(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for overturning_peat_shallow. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(overturning_peat_shallow_leafless(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for overturning_peat_shallow_leafless. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(overturning_peat_deep(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for overturning_peat_deep. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(overturning_peat_deep_leafless(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for overturning_peat_deep_leafless. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(overturning_peat_average(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for overturning_peat_average. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(overturning_peat_average_leafless(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for overturning_peat_average_leafles. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(max_damage_further(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for max_damage_further. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(max_damage_closer(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for max_damage_cloeer. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(sfactor_further(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for sfactor_further. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(sfactor_closer(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for sfactor_closer. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(green_density(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for green_density. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

       ALLOCATE(som_turn_ipassive(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for som_turn_ipassive. We stop.We need nvmwords = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF


       !
       ! CROPLAND MANAGEMENT
       !
       ALLOCATE(harvest_ratio(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for harvest_ratio. We stop. We need nvm words = ',nvm
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       END IF

    ENDIF ! (ok_stomate)

    ALLOCATE(rotation_ref(nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'pft_parameters','Pb alloc rotation_ref','','')
    rotation_ref(:) = -un                      ! sentinelle : inactif => bit-neutre
    ALLOCATE(fm_src_pft(nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'pft_parameters','Pb alloc fm_src_pft','','')
    fm_src_pft(:) = .FALSE.
    ALLOCATE(fm_cop_pft(nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'pft_parameters','Pb alloc fm_cop_pft','','')
    fm_cop_pft(:) = .FALSE.

    !! Following parameters are used with and without ok_stomate

    ALLOCATE(nue_opt(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for nue_opt. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(vmax_uptake(nvm,nionspec),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for nue_opt. We stop. We need nvm*nionspec words = ',nvm,nionspec
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(cn_leaf_min(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for cn_leaf_min. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF
    
    ALLOCATE(cn_leaf_max(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for cn_leaf_max. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(cn_leaf_init(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for cn_leaf_init. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF
  
    ALLOCATE(ext_coeff_N(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for ext_coeff_N. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(maint_resp_slope(nvm,3),stat=ier) 
    l_error = l_error .OR. (ier /= 0) 
    IF (l_error) THEN 
       WRITE(numout,*) ' Memory allocation error for maint_resp_slope. We stop. We need nvm*3 words = ',nvm*3 
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF
    maint_resp_slope(:,:) = zero 
    
    ALLOCATE(maint_resp_slope_c(nvm),stat=ier) 
    l_error = l_error .OR. (ier /= 0) 
    IF (l_error) THEN 
       WRITE(numout,*) ' Memory allocation error for maint_resp_slope_c. We stop. We need nvm words = ',nvm 
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF
    
    ALLOCATE(maint_resp_slope_b(nvm),stat=ier) 
    l_error = l_error .OR. (ier /= 0) 
    IF (l_error) THEN 
       WRITE(numout,*) ' Memory allocation error for maint_resp_slope_b. We stop. We need nvm words = ',nvm 
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF
    
    ALLOCATE(maint_resp_slope_a(nvm),stat=ier) 
    l_error = l_error .OR. (ier /= 0) 
    IF (l_error) THEN 
       WRITE(numout,*) ' Memory allocation error for maint_resp_slope_a. We stop. We need nvm words = ',nvm 
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(dia_thresh_inv(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for dia_thresh_inv. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(longevity_root(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for longevity_root. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(longevity_sap(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for longevity_sap. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(k_latosa_min(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for k_latosa_min. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(k_latosa_max(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for k_latosa_max. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF

    ALLOCATE(senescence_type(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for . We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF


    ALLOCATE(moss_frac(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for moss_frac. We stop. We need nvm words = ',nvm
       CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
    END IF
    
    ALLOCATE(is_peat(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for is_peat. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF
    ALLOCATE(is_croppeat(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for is_croppeat. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

  END SUBROUTINE pft_parameters_alloc

!! ================================================================================================================================
!! SUBROUTINE   : config_pft_parameters 
!!
!>\BRIEF          This subroutine will read the imposed values for the global pft
!! parameters (sechiba + stomate). It is not called if IMPOSE_PARAM is set to NO.
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE config_pft_parameters

    IMPLICIT NONE

    !! 0. Variables and parameters declaration

    !! 0.4 Local variable

    INTEGER(i_std) :: jv, ivm                   !! Index (untiless)

    !_ ================================================================================================================================ 


    !
    ! Vegetation structure
    !

    !Config Key   = LEAF_TAB
    !Config Desc  = leaf type : 1=broad leaved tree, 2=needle leaved tree, 3=grass 4=bare ground
    !Config if    = OK_STOMATE
    !Config Def   = 4, 1, 1, 2, 1, 1, 2, 1, 2, 3, 3, 3, 3 
    !Config Help  = 
    !Config Units = [-] 
    CALL getin_p('LEAF_TAB',leaf_tab)

    !Config Key   = PHENO_MODEL
    !Config Desc  = which phenology model is used? (tabulated) 
    !Config if    = OK_STOMATE
    !Config Def   = none, none, moi, none, none, ncdgdd, none, ncdgdd, ngd, moigdd, moigdd, moigdd, moigdd
    !Config Help  =
    !Config Units = [-] 
    CALL getin_p('PHENO_MODEL',pheno_model)

    !! Redefine the values for is_tree, is_deciduous, is_needleleaf, is_evergreen if values have been modified
    !! in run.def

    is_tree(:) = .FALSE.
    DO jv = 1,nvm
       IF ( leaf_tab(jv) <= 2 ) is_tree(jv) = .TRUE.
    END DO
    !
    is_deciduous(:) = .FALSE.
    DO jv = 1,nvm
       IF ( is_tree(jv) .AND. (pheno_model(jv) /= "none") ) is_deciduous(jv) = .TRUE.
    END DO
    !
    is_evergreen(:) = .FALSE.
    DO jv = 1,nvm
       IF ( is_tree(jv) .AND. (pheno_model(jv) == "none") ) is_evergreen(jv) = .TRUE.
    END DO
    !
    is_needleleaf(:) = .FALSE.
    DO jv = 1,nvm
       IF ( leaf_tab(jv) == 2 ) is_needleleaf(jv) = .TRUE.
    END DO

    !Config Key   = GDD_THRESHOLD 
    !Config Desc  = Temperature above which GDD is accounted for
    !Config If    = OK_STOMATE 
    !Config Def   = 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 5, 5, 5 
    !Config Help  = GDD : Growing-Degree-Day
    !Config Units = [Celcius] 
    CALL getin_p('GDD_THRESHOLD',gdd_threshold)

    !Config Key   = SECHIBA_LAI
    !Config Desc  = laimax for maximum lai(see also type of lai interpolation)
    !Config if    = OK_SECHIBA or IMPOSE_VEG
    !Config Def   = 0., 8., 8., 4., 4.5, 4.5, 4., 4.5, 4., 2., 2., 2., 2.
    !Config Help  = Maximum values of lai used for interpolation of the lai map
    !Config Units = [m^2/m^2]
    CALL getin_p('SECHIBA_LAI',llaimax)

    !Config Key   = LLAIMIN
    !Config Desc  = laimin for minimum lai(see also type of lai interpolation)
    !Config if    = OK_SECHIBA or IMPOSE_VEG
    !Config Def   = 0., 8., 0., 4., 4.5, 0., 4., 0., 0., 0., 0., 0., 0.
    !Config Help  = Minimum values of lai used for interpolation of the lai map
    !Config Units = [m^2/m^2]
    CALL getin_p('LLAIMIN',llaimin)

    !Config Key   = SLOWPROC_HEIGHT
    !Config Desc  = prescribed height of vegetation 
    !Config if    = OK_SECHIBA
    !Config Def   = 0., 30., 30., 20., 20., 20., 15., 15., 15., .5, .6, 1., 1.
    !Config Help  =
    !Config Units = [m] 
    CALL getin_p('SLOWPROC_HEIGHT',height_presc)

    !Config Key   = Z0_OVER_HEIGHT
    !Config Desc  = factor to calculate roughness height from height of canopy 
    !Config if    = OK_SECHIBA
    !Config Def   = 0., 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625
    !Config Help  =
    !Config Units = [-] 
    CALL getin_p('Z0_OVER_HEIGHT',z0_over_height)

    !
    !Config Key   = RATIO_Z0M_Z0H
    !Config Desc  = Ratio between z0m and z0h
    !Config Def   = 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0 
    !Config if    = OK_SECHIBA
    !Config Help  = 
    !Config Units = [-]
    CALL getin_p('RATIO_Z0M_Z0H',ratio_z0m_z0h)


    !Config Key   = TYPE_OF_LAI
    !Config Desc  = Type of behaviour of the LAI evolution algorithm 
    !Config if    = OK_SECHIBA
    !Config Def   = inter, inter, inter, inter, inter, inter, inter, inter, inter, inter, inter, inter, inter
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('TYPE_OF_LAI',type_of_lai)

    !Config Key   = NATURAL
    !Config Desc  = natural? 
    !Config if    = OK_SECHIBA, OK_STOMATE
    !Config Def   = y, y, y, y, y, y, y, y, y, y, y, n, n 
    !Config Help  =
    !Config Units = [BOOLEAN]
    CALL getin_p('NATURAL',natural)

    !Config Key   = IS_TROPICAL
    !Config Desc  = PFT IS TROPICAL 
    !Config if    = OK_STOMATE
    !Config Def   = FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('IS_TROPICAL',is_tropical)

    !Config Key   = IS_TEMPERATE
    !Config Desc  = PFT IS TEMPERATE 
    !Config if    = OK_STOMATE
    !Config Def   = FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE 
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('IS_TEMPERATE',is_temperate)
    
    !Config Key   = IS_BOREAL
    !Config Desc  = PFT IS BOREAL 
    !Config if    = OK_STOMATE
    !Config Def   = FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE 
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('IS_BOREAL',is_boreal)


    !
    ! Peatland
    !

    !Config Key   = OK_PEAT_HYDRO
    !Config Desc  = Activate peatland in sechiba
    !Config if    = 
    !Config Def   = FALSE
    !Config Help  =
    !Config Units = [BOOLEAN]
    ok_peat_hydro= .FALSE.
    CALL getin_p('OK_PEAT_HYDRO', ok_peat_hydro)

    !Config Key   = PEAT_TIDES
    !Config Desc  = Activate peatland in sechiba
    !Config if    = OK_PEAT_HYDRO
    !Config Def   = FALSE
    !Config Help  =
    !Config Units = [BOOLEAN]
    tides=.FALSE.
    CALL getin_p('PEAT_TIDES', tides)

    !Config Key   = IS_PEAT
    !Config Desc  = flag by PFT for peatland
    !Config if    = 
    !Config Def   = FALSE (if not OK_PEAT_HYDRO)
    !Config Help  =
    !Config Units = [BOOLEAN]
    IF (ok_peat_hydro .AND. .NOT. tides .AND. .NOT. agri_peat) THEN
       is_peat(:)=.FALSE.
       is_peat(nvm)=.TRUE.
    END IF
    CALL getin_p('IS_PEAT',is_peat)

    !Config Key   = IS_CROPPEAT
    !Config Desc  = flag by PFT for peatland agriculture
    !Config if    = 
    !Config Def   = FALSE (if not AGRI_PEAT)
    !Config Help  =
    !Config Units = [BOOLEAN]
    CALL getin_p('IS_CROPPEAT',is_croppeat)
    
    ! Default parameter setting for peatland
    IF (ok_peat_hydro .AND. .NOT. tides .AND. .NOT. agri_peat) THEN
       nstm= 4
       ! By defaut, the last PFT is dedicated to peatland
       pref_soil_veg(nvm) = 4 
       Vcmax25(nvm)= 46.7
       sla(nvm) = 0.0153
       arJV(nvm) = 2.805

       perma_peat = .TRUE.
    ENDIF

    IF (tides) THEN
       ! Warning, for this case, the parameters need to be set in run.def(or orchidee_pft.def), 
       ! the default settings here are not enough
       nstm= 6
       pref_soil_veg(nvm) = 6
       Vcmax25(nvm)= 46.7
       sla(nvm) = 0.0153
       arJV(nvm) = 2.805
       perma_peat = .TRUE.
    ENDIF


    !
    ! Photosynthesis
    !

    !Config Key   = IS_C4
    !Config Desc  = flag for C4 vegetation types
    !Config if    = OK_SECHIBA or OK_STOMATE
    !Config Def   = n, n, n, n, n, n, n, n, n, n, n, y, n, y
    !Config Help  =
    !Config Units = [BOOLEAN]
    CALL getin_p('IS_C4',is_c4)

    !Config Key   = VCMAX_FIX
    !Config Desc  = values used for vcmax when STOMATE is not activated
    !Config if    = OK_SECHIBA and NOT(OK_STOMATE)
    !Config Def   = 0., 40., 50., 30., 35., 40.,30., 40., 35., 60., 60., 70., 70.
    !Config Help  =
    !Config Units = [micromol/m^2/s] 
    CALL getin_p('VCMAX_FIX',vcmax_fix)

    !Config Key   = DOWNREG_CO2
    !Config Desc  = coefficient for CO2 downregulation (unitless)
    !Config if    = 
    !Config Def   = 0., 0.38, 0.38, 0.28, 0.28, 0.28, 0.22, 0.22, 0.22, 0.26, 0.26, 0.26, 0.26
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('DOWNREG_CO2',downregulation_co2_coeff)

    !Config Key   = E_KmC
    !Config Desc  = Energy of activation for KmC
    !Config if    = 
    !Config Def   = undef,  79430., 79430., 79430., 79430., 79430., 79430., 79430., 79430., 79430., 79430., 79430., 79430.
    !Config Help  = See Medlyn et al. (2002) 
    !Config Units = [J mol-1]
    CALL getin_p('E_KMC',E_KmC)

    !Config Key   = E_KmO
    !Config Desc  = Energy of activation for KmO
    !Config if    = 
    !Config Def   = undef, 36380.,  36380.,  36380.,  36380.,  36380., 36380., 36380., 36380., 36380., 36380., 36380., 36380.
    !Config Help  = See Medlyn et al. (2002) 
    !Config Units = [J mol-1]
    CALL getin_p('E_KMO',E_KmO)

    !Config Key   = E_Sco
    !Config Desc  = Energy of activation for Sco
    !Config if    = 
    !Config Def   = undef, -24460., -24460., -24460., -24460., -24460., -24460., -24460., -24460., -24460., -24460., -24460., -24460.
    !Config Help  = See Table 2 of Yin et al. (2009) - Value for C4 plants is not mentioned - We use C3 for all plants
    !Config Units = [J mol-1]
    CALL getin_p('E_SCO',E_Sco)
    
    !Config Key   = E_gamma_star
    !Config Desc  = Energy of activation for gamma_star
    !Config if    = 
    !Config Def   = undef, 37830.,  37830.,  37830.,  37830.,  37830., 37830., 37830., 37830., 37830., 37830., 37830., 37830.
    !Config Help  = See Medlyn et al. (2002) from Bernacchi al. (2001) 
    !Config Units = [J mol-1]
    CALL getin_p('E_GAMMA_STAR',E_gamma_star)

    !Config Key   = E_Vcmax
    !Config Desc  = Energy of activation for Vcmax
    !Config if    = 
    !Config Def   = undef, 71513., 71513., 71513., 71513., 71513., 71513., 71513., 71513., 71513., 67300., 71513., 67300.
    !Config Help  = See Table 2 of Yin et al. (2009) for C4 plants and Kattge & Knorr (2007) for C3 plants (table 3)
    !Config Units = [J mol-1]
    CALL getin_p('E_VCMAX',E_Vcmax)

    !Config Key   = E_Jmax
    !Config Desc  = Energy of activation for Jmax
    !Config if    = 
    !Config Def   = undef, 49884., 49884., 49884., 49884., 49884., 49884., 49884., 49884., 49884., 77900., 49884., 77900. 
    !Config Help  = See Table 2 of Yin et al. (2009) for C4 plants and Kattge & Knorr (2007) for C3 plants (table 3)
    !Config Units = [J mol-1]
    CALL getin_p('E_JMAX',E_Jmax)

    !Config Key   = aSV
    !Config Desc  = a coefficient of the linear regression (a+bT) defining the Entropy term for Vcmax
    !Config if    = 
    !Config Def   = undef, 668.39, 668.39, 668.39, 668.39, 668.39, 668.39, 668.39, 668.39, 668.39, 641.64, 668.39, 641.64 
    !Config Help  = See Table 3 of Kattge & Knorr (2007) - For C4 plants, we assume that there is no acclimation and that at for a temperature of 25°C, aSV is the same for both C4 and C3 plants (no strong jusitification - need further parametrization)
    !Config Units = [J K-1 mol-1]
    CALL getin_p('ASV',aSV)

    !Config Key   = bSV
    !Config Desc  = b coefficient of the linear regression (a+bT) defining the Entropy term for Vcmax
    !Config if    = 
    !Config Def   = undef, -1.07, -1.07, -1.07, -1.07, -1.07, -1.07, -1.07, -1.07, -1.07, 0., -1.07, 0. 
    !Config Help  = See Table 3 of Kattge & Knorr (2007) - For C4 plants, we assume that there is no acclimation
    !Config Units = [J K-1 mol-1 °C-1]
    CALL getin_p('BSV',bSV)

    !Config Key   = TPHOTO_MIN
    !Config Desc  = minimum photosynthesis temperature (deg C)
    !Config if    = OK_STOMATE
    !Config Def   = undef,  -4., -4., -4., -4.,-4.,-4., -4., -4., -4., -4., -4., -4.
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('TPHOTO_MIN',tphoto_min)

    !Config Key   = TPHOTO_MAX
    !Config Desc  = maximum photosynthesis temperature (deg C)
    !Config if    = OK_STOMATE
    !Config Def   = undef, 55., 55., 55., 55., 55., 55., 55., 55., 55., 55., 55., 55.
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('TPHOTO_MAX',tphoto_max)

    !Config Key   = aSJ
    !Config Desc  = a coefficient of the linear regression (a+bT) defining the Entropy term for Jmax
    !Config if    = 
    !Config Def   = undef, 659.70, 659.70, 659.70, 659.70, 659.70, 659.70, 659.70, 659.70, 659.70, 630., 659.70, 630. 
    !Config Help  = See Table 3 of Kattge & Knorr (2007) - and Table 2 of Yin et al. (2009) for C4 plants
    !Config Units = [J K-1 mol-1]
    CALL getin_p('ASJ',aSJ)

    !Config Key   = bSJ
    !Config Desc  = b coefficient of the linear regression (a+bT) defining the Entropy term for Jmax
    !Config if    = 
    !Config Def   = undef, -0.75, -0.75, -0.75, -0.75, -0.75, -0.75, -0.75, -0.75, -0.75, 0., -0.75, 0. 
    !Config Help  = See Table 3 of Kattge & Knorr (2007) - For C4 plants, we assume that there is no acclimation
    !Config Units = [J K-1 mol-1 °C-1]
    CALL getin_p('BSJ',bSJ)

    !Config Key   = D_Vcmax
    !Config Desc  = Energy of deactivation for Vcmax
    !Config if    = 
    !Config Def   = undef, 200000., 200000., 200000., 200000., 200000., 200000., 200000., 200000., 200000., 192000., 200000., 192000.
    !Config Help  = Medlyn et al. (2002) also uses 200000. for C3 plants (same value than D_Jmax). 'Consequently', we use the value of D_Jmax for C4 plants.
    !Config Units = [J mol-1]
    CALL getin_p('D_VCMAX',D_Vcmax)

    !Config Key   = D_Jmax
    !Config Desc  = Energy of deactivation for Jmax
    !Config if    = 
    !Config Def   = undef, 200000., 200000., 200000., 200000., 200000., 200000., 200000., 200000., 200000., 192000., 200000., 192000.
    !Config Help  = See Table 2 of Yin et al. (2009)
    !Config Units = [J mol-1]
    CALL getin_p('D_JMAX',D_Jmax)
    
    !Config Key   = E_gm 
    !Config Desc  = Energy of activation for gm 
    !Config if    =  
    !Config Def   = undef, 49600., 49600., 49600., 49600., 49600., 49600., 49600., 49600., 49600., undef, 49600., undef 
    !Config Help  = See Table 2 of Yin et al. (2009) 
    !Config Units = [J mol-1] 
    CALL getin_p('E_GM',E_gm) 
    
    !Config Key   = S_gm 
    !Config Desc  = Entropy term for gm 
    !Config if    =  
    !Config Def   = undef, 1400., 1400., 1400., 1400., 1400., 1400., 1400., 1400., 1400., undef, 1400., undef 
    !Config Help  = See Table 2 of Yin et al. (2009) 
    !Config Units = [J K-1 mol-1] 
    CALL getin_p('S_GM',S_gm) 
    
    !Config Key   = D_gm 
    !Config Desc  = Energy of deactivation for gm 
    !Config if    =  
    !Config Def   = undef, 437400., 437400., 437400., 437400., 437400., 437400., 437400., 437400., 437400., undef, 437400., undef 
    !Config Help  = See Table 2 of Yin et al. (2009) 
    !Config Units = [J mol-1] 
    CALL getin_p('D_GM',D_gm) 
    
    !Config Key   = E_Rd
    !Config Desc  = Energy of activation for Rd
    !Config if    = 
    !Config Def   = undef, 46390., 46390., 46390., 46390., 46390., 46390., 46390., 46390., 46390., 46390., 46390., 46390.
    !Config Help  = See Table 2 of Yin et al. (2009)
    !Config Units = [J mol-1]
    CALL getin_p('E_RD',E_Rd)

    !Config Key   = VCMAX25
    !Config Desc  = Maximum rate of Rubisco activity-limited carboxylation at 25°C
    !Config if    = OK_STOMATE
    !Config Def   = undef, 45.0, 45.0, 35.0, 40.0, 50.0, 45.0, 35.0, 35.0, 50.0, 50.0, 60.0, 60.0
    !Config Help  = Notice that, with the introduction of the nitrogen cycle, this
    !Config         parameter is no longer used to influence the simulation.  It is kept solely as
    !Config         a way to compare to old revisions (nue_opt is the new parameter that controls
    !Config         photosynthesis in this way).
    !Config Units = [micromol/m^2/s]
    CALL getin_p('VCMAX25',Vcmax25)

    !Config Key   = ARJV
    !Config Desc  = a coefficient of the linear regression (a+bT) defining the Jmax25/Vcmax25 ratio 
    !Config if    = OK_STOMATE
    !Config Def   = undef, 2.59, 2.59, 2.59, 2.59, 2.59, 2.59, 2.59, 2.59, 2.59, 1.715, 2.59, 1.715
    !Config Help  = See Table 3 of Kattge & Knorr (2007) - For C4 plants, we assume that there is no acclimation and that for a temperature of 25°C, aSV is the same for both C4 and C3 plants (no strong jusitification - need further parametrization)
    !Config Units = [mu mol e- (mu mol CO2)-1]
    CALL getin_p('ARJV',arJV)

    !Config Key   = BRJV
    !Config Desc  = b coefficient of the linear regression (a+bT) defining the Jmax25/Vcmax25 ratio 
    !Config if    = OK_STOMATE
    !Config Def   = undef, -0.035, -0.035, -0.035, -0.035, -0.035, -0.035, -0.035, -0.035, -0.035, 0., -0.035, 0.
    !Config Help  = See Table 3 of Kattge & Knorr (2007) -  We assume No acclimation term for C4 plants
    !Config Units = [(mu mol e- (mu mol CO2)-1) (°C)-1]
    CALL getin_p('BRJV',brJV)

    !Config Key   = KmC25
    !Config Desc  = Michaelis–Menten constant of Rubisco for CO2 at 25°C
    !Config if    = 
    !Config Def   = undef, 404.9, 404.9, 404.9, 404.9, 404.9, 404.9, 404.9, 404.9, 404.9, 650., 404.9, 650.
    !Config Help  = See Table 2 of Yin et al. (2009) for C4 plants and Medlyn et al. (2002) for C3 plants
    !Config Units = [ubar]
    CALL getin_p('KMC25',KmC25)

    !Config Key   = KmO25
    !Config Desc  = Michaelis–Menten constant of Rubisco for O2 at 25°C
    !Config if    = 
    !Config Def   = undef, 278400., 278400., 278400., 278400., 278400., 278400., 278400., 278400., 278400., 450000., 278400., 450000.
    !Config Help  = See Table 2 of Yin et al. (2009) for C4 plants and Medlyn et al. (2002) for C3 plants
    !Config Units = [ubar]
    CALL getin_p('KMO25',KmO25)

    !Config Key   = Sco25
    !Config Desc  = Relative CO2 /O2 specificity factor for Rubisco at 25Â°C
    !Config if    = 
    !Config Def   = undef, 2800., 2800., 2800., 2800., 2800., 2800., 2800., 2800., 2800., 2590., 2800., 2590.
    !Config Help  = See Table 2 of Yin et al. (2009)
    !Config Units = [bar bar-1]
    CALL getin_p('SCO25',Sco25)
    
    !Config Key   = gm25 
    !Config Desc  = Mesophyll diffusion conductance at 25ÃÂ°C 
    !Config if    =  
    !Config Def   = undef, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, undef, 0.4, undef 
    !Config Help  = See legend of Figure 6 of Yin et al. (2009) and review by Flexas et al. (2008) - gm is not used for C4 plants 
    !Config Units = [mol m-2 s-1 bar-1] 
    CALL getin_p('GM25',gm25) 
    
    !Config Key   = gamma_star25
    !Config Desc  = Ci-based CO2 compensation point in the absence of Rd at 25°C (ubar)
    !Config if    = 
    !Config Def   = undef, 42.75, 42.75, 42.75, 42.75, 42.75, 42.75, 42.75, 42.75, 42.75, 42.75, 42.75, 42.75
    !Config Help  = See Medlyn et al. (2002) for C3 plants - For C4 plants, we use the same value (probably uncorrect)
    !Config Units = [ubar]
    CALL getin_p('gamma_star25',gamma_star25)

    !Config Key   = a1
    !Config Desc  = Empirical factor involved in the calculation of fvpd
    !Config if    = 
    !Config Def   = undef, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.72, 0.85, 0.72
    !Config Help  = See Table 2 of Yin et al. (2009)
    !Config Units = [-]
    CALL getin_p('A1',a1)

    !Config Key   = b1
    !Config Desc  = Empirical factor involved in the calculation of fvpd
    !Config if    = 
    !Config Def   = undef, 0.14, 0.14, 0.14, 0.14, 0.14, 0.14, 0.14, 0.14, 0.14, 0.20, 0.14, 0.20
    !Config Help  = See Table 2 of Yin et al. (2009)
    !Config Units = [-]
    CALL getin_p('B1',b1)

    !Config Key   = a_fpsi
    !Config Desc  = Empirical factor involved in the calculation of fpsi
    !Config if    = 
    !Config Def   = 1. / ( 1. / a1 - 1 )
    !Config Help  = 
    !Config Units = [-]
    CALL getin_p('A_FPSI',a_fpsi)

    !Config Key   = sf
    !Config Desc  = Sensibility parameter involved in the calculation of fpsi
    !Config if    = 
    !Config Def   = undef, 4.3, 4.3, 4.3, 4.3, 4.3, 4.3, 4.3, 4.3, 4.3, 4.3, 4.3, 4.3
    !Config Help  = See Table 1 of Tuzet et al. (2017)
    !Config Units = [-]
    CALL getin_p('SF',sf)

    !Config Key   = Rroot_sup
    !Config Desc  = Superficial roots water resistance involved in the calculation of psi_leaf
    !Config if    = 
    !Config Def   = undef, 500000, 500000, 500000, 500000, 500000, 500000, 500000, 500000, 500000, 500000, 500000, 500000
    !Config Help  = See Table 1 of Tuzet et al. (2017)
    !Config Units = [-]
    CALL getin_p('RROOT_SUP',Rroot_sup)

    !Config Key   = Rroot_inf
    !Config Desc  = Inferior roots water resistance involved in the calculation of psi_leaf
    !Config if    = 
    !Config Def   = undef, 2000000, 2000000, 2000000, 2000000, 2000000, 2000000, 2000000, 2000000, 2000000, 2000000, 2000000, 2000000
    !Config Help  = See Table 1 of Tuzet et al. (2017)
    !Config Units = [-]
    CALL getin_p('RROOT_INF',Rroot_inf)

    !Config Key   = Rleaf
    !Config Desc  = Leaf water resistance involved in the calculation of psi_leaf
    !Config if    = 
    !Config Def   = undef, 500000, 500000, 500000, 500000, 500000, 500000, 500000, 500000, 500000, 500000, 500000, 500000
    !Config Help  = See Table 1 of Tuzet et al. (2017)
    !Config Units = [-]
    CALL getin_p('RLEAF',RLEAF)

    !Config Key   = Rxylem
    !Config Desc  = Xylem water resistance involved in the calculation of psi_leaf
    !Config if    = 
    !Config Def   = undef, 2500000, 2500000, 2500000, 2500000, 2500000, 2500000, 2500000, 2500000, 2500000, 2500000, 2500000, 2500000
    !Config Help  = See Table 1 of Tuzet et al. (2017)
    !Config Units = [-]
    CALL getin_p('RXYLEM',Rxylem)

    !Config Key   = Rroot_soil
    !Config Desc  = Soil to roots water resistance involved in the calculation of psi_leaf
    !Config if    = 
    !Config Def   = undef, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150
    !Config Help  = See Table 1 of Liu et al. (2015)
    !Config Units = [-]
    CALL getin_p('RROOT_SOIL',Rroot_soil)

    !Config Key   = Rsto_leaf
    !Config Desc  = Leaf water storage resistance involved in the calculation of psi_leaf
    !Config if    = 
    !Config Def   = undef, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000
    !Config Help  = See Table 1 of Liu et al. (2015)
    !Config Units = [-]
    CALL getin_p('RSTO_LEAF',Rsto_leaf)

    !Config Key   = Rsto_wood
    !Config Desc  = Wood water storage resistance involved in the calculation of psi_leaf
    !Config if    = 
    !Config Def   = undef, 60000000, 60000000, 60000000, 60000000, 60000000, 60000000, 60000000, 60000000, 60000000, 60000000, 60000000, 60000000
    !Config Help  = See Table 1 of Liu et al. (2015)
    !Config Units = [-]
    CALL getin_p('RSTO_WOOD',Rsto_wood)

    !Config Key   = psi_ref_g0
    !Config Desc  = Reference Potential involved in the calculation of fpsi in the calculation of stomatal conductance
    !Config if    = 
    !Config Def   = undef, -1.6, -1.6, -1.6, -1.6, -1.6, -1.6, -1.6, -1.6, -1.6, -1.6, -1.6, -1.6
    !Config Help  = See Table 1 of Tuzet et al. (2017)
    !Config Units = [-]
    CALL getin_p('PSI_REF_G0',psi_ref_g0)

    !Config Key   = psi_ref_sto_leaf
    !Config Desc  = Reference Potential involved in the calculation of leaf water storage potential
    !Config if    = 
    !Config Def   = undef, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0
    !Config Help  = See Table 1 of Tuzet et al. (2017)
    !Config Units = [-]
    CALL getin_p('PSI_REF_STO_LEAF',psi_ref_sto_leaf)

    !Config Key   = psi_ref_sto_wood
    !Config Desc  = Reference Potential involved in the calculation of wood water storage potential
    !Config if    = 
    !Config Def   = undef, -1.4, -1.4, -1.4, -1.4, -1.4, -1.4, -1.4, -1.4, -1.4, -1.4, -1.4, -1.4
    !Config Help  = See Table 1 of Tuzet et al. (2017)
    !Config Units = [-]
    CALL getin_p('PSI_REF_STO_WOOD',psi_ref_sto_wood)
    
    !Config Key   = psi_min
    !Config Desc  = It is used to avoid numerical problems with resistances and potentials
    !Config if    = 
    !Config Def   = -20.
    !Config Help  = 
    !Config Units = [-]
    CALL getin_p('PSI_MIN',psi_min)

    !Config Key   = psi_50_res_leaf
    !Config Desc  = Reference Potential involved in the calculation of leaf water resistance
    !Config if    = 
    !Config Def   = undef, -1.1, -1.1, -1.1, -1.1, -1.1, -1.1, -1.1, -1.1, -1.1, -1.1, -1.1, -1.1
    !Config Help  = See Parameterization of Yao et al. (2021) (Under review at the time this part is developed)
    !Config Units = [-]
    CALL getin_p('PSI_50_RES_LEAF',psi_50_res_leaf)

    !Config Key   = psi_50_res_wood
    !Config Desc  = Reference Potential involved in the calculation of wood water resistance
    !Config if    = 
    !Config Def   = undef, -1.3, -1.3, -1.3, -1.3, -1.3, -1.3, -1.3, -1.3, -1.3, -1.3, -1.3, -1.3
    !Config Help  = See Parameterization of Yao et al. (2021) (Under review at the time this part is developed)
    !Config Units = [-]
    CALL getin_p('PSI_50_RES_WOOD',psi_50_res_wood)

    !Config Key   = psi_50_res_root
    !Config Desc  = Reference Potential involved in the calculation of root water resistance
    !Config if    = 
    !Config Def   = undef, -1.1, -1.1, -1.1, -1.1, -1.1, -1.1, -1.1, -1.1, -1.1, -1.1, -1.1, -1.1
    !Config Help  = See Parameterization of Yao et al. (2021) (Under review at the time this part is developed)
    !Config Units = [-]
    CALL getin_p('PSI_50_RES_ROOT',psi_50_res_root)

    !Config Key   = a_leaf
    !Config Desc  = Parameter involved in the calculation of leaf water resistance
    !Config if    = 
    !Config Def   = undef, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5, -2.5
    !Config Help  = See Parameterization of Yao et al. (2021) (Under review at the time this part is developed)
    !Config Units = [-]
    CALL getin_p('A_LEAF',a_leaf)

    !Config Key   = a_wood
    !Config Desc  = Parameter involved in the calculation of wood water resistance
    !Config if    = 
    !Config Def   = undef, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0
    !Config Help  = See Parameterization of Yao et al. (2021) (Under review at the time this part is developed)
    !Config Units = [-]
    CALL getin_p('A_WOOD',a_wood)

    !Config Key   = a_root
    !Config Desc  = Parameter involved in the calculation of root water resistance
    !Config if    = 
    !Config Def   = undef, -2.3, -2.3, -2.3, -2.3, -2.3, -2.3, -2.3, -2.3, -2.3, -2.3, -2.3, -2.3
    !Config Help  = See Parameterization of Yao et al. (2021) (Under review at the time this part is developed)
    !Config Units = [-]
    CALL getin_p('A_ROOT',a_root)

    !Config Key   = k_leaf_max
    !Config Desc  = Maximum leaf water conductance involved in the calculation of leaf water resistance
    !Config if    = 
    !Config Def   = undef, 0.000000525, 0.000000525, 0.000000525, 0.000000525, 0.000000525, 0.000000525, 0.000000525, 0.000000525, 0.000000525, 0.000000525, 0.000000525, 0.000000525
    !Config Help  = See Parameterization of Yao et al. (2021) (Under review at the time this part is developed)
    !Config Units = [-]
    CALL getin_p('K_LEAF_MAX',k_leaf_max)

    !Config Key   = sapwood_density
    !Config Desc  = Maximum sapwood water conductivity
    !Config if    = 
    !Config Def   = /undef,   5.3e5,  5.3e5,  5.3e5,  6.3e5,  5.9e5,  4.9e5,  4.1e5,  5.2e5,
    !Config Def   =  undef,    undef,  undef,  undef/
    !Config Help  = based on XFT database classified by PFT
    !Config Units = [-]
    CALL getin_p('SAPWOOD_DENSITY',sapwood_density)

    !Config Key   = ks_max
    !Config Desc  = Maximum sapwood water conductivity
    !Config if    = 
    !Config Def   = /undef,  0.00357,  0.00613,  0.00132,  0.00327,  0.00327,  0.00130,  0.00268,
    !Config Def   = 0.00030, undef,    undef,    undef,    undef
    !Config Help  = from XFT database classified by PFT
    !Config Units = [-]
    CALL getin_p('KS_MAX',ks_max)

    !Config Key   = k_wood_max
    !Config Desc  = Maximum wood water conductance involved in the calculation of wood water resistance
    !Config if    = 
    !Config Def   = undef, 0.000000175, 0.000000175, 0.000000175, 0.000000175, 0.000000175, 0.000000175, 0.000000175, 0.000000175, 0.000000175, 0.000000175, 0.000000175, 0.000000175
    !Config Help  = See Parameterization of Yao et al. (2021) (Under review at the time this part is developed)
    !Config Units = [-]
    CALL getin_p('K_WOOD_MAX',k_wood_max)

    !Config Key   = k_root_max
    !Config Desc  = Maximum root water conductance involved in the calculation of root water resistance
    !Config if    = 
    !Config Def   = undef, 0.000000175, 0.000000175, 0.000000175, 0.000000175, 0.000000175, 0.000000175, 0.000000175, 0.000000175, 0.000000175, 0.000000175, 0.000000175, 0.000000175
    !Config Help  = See Parameterization of Yao et al. (2021) (Under review at the time this part is developed)
    !Config Units = [-]
    CALL getin_p('K_ROOT_MAX',k_root_max)

    !Config Key   = lambda_leaf
    !Config Desc  = Parameter involved in the calculation of leaf water storage potential
    !Config if    = 
    !Config Def   = undef, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3
    !Config Help  = See Table 1 of Tuzet et al. (2017)
    !Config Units = [-]
    CALL getin_p('LAMBDA_LEAF',lambda_leaf)

    !Config Key   = lambda_wood
    !Config Desc  = Parameter involved in the calculation of wood water storage potential
    !Config if    = 
    !Config Def   = undef, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9
    !Config Help  = See Table 1 of Tuzet et al. (2017)
    !Config Units = [-]
    CALL getin_p('LAMBDA_WOOD',lambda_wood)

    !Config Key   = d_l
    !Config Desc  = Parameter for stomatal conductance sensitivity to VPD (kPa) 
    !Config if    = 
    !Config Def   = undef, 5.5, 5.5, 5.5, 5.5, 5.5, 5.5, 5.5, 5.5, 5.5, 5.5, 5.5, 5.5, 5.5, 5.5
    !Config Help  = Value from Anderegg et al. (2017) 
    !Config Units = [kPa]
    CALL getin_p('D_L',d_l)
    
    !Config Key   = g0
    !Config Desc  = Residual stomatal conductance when irradiance approaches zero 
    !Config if    = 
    !Config Def   = undef, 0.0049, 0.0049, 0.0049, 0.0049, 0.0049, 0.0049, 0.0049, 0.0049, 0.0049, 0.01875, 0.00625, 0.01875
    !Config Help  = Value for tree PFTs from Duursma et al. (2019) others from ORCHIDEE - No other reference.
    !Config Units = [mol m−2 s−1 bar−1]
    CALL getin_p('G0',g0)

    !Config Key   = h_protons
    !Config Desc  = Number of protons required to produce one ATP
    !Config if    = 
    !Config Def   = undef, 4., 4., 4., 4., 4., 4., 4., 4., 4., 4., 4., 4. 
    !Config Help  = See Table 2 of Yin et al. (2009) - h parameter
    !Config Units = [mol mol-1]
    CALL getin_p('H_PROTONS',h_protons)

    !Config Key   = fpsir
    !Config Desc  = Fraction of PSII e− transport rate partitioned to the C4 cycle
    !Config if    = 
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, 0.4, undef, 0.4 
    !Config Help  = See Table 2 of Yin et al. (2009)
    !Config Units = [-]
    CALL getin_p('FPSIR',fpsir)

    !Config Key   = fQ
    !Config Desc  = Fraction of electrons at reduced plastoquinone that follow the Q-cycle
    !Config if    = 
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, 1., undef, 1.
    !Config Help  = See Table 2 of Yin et al. (2009) - Values for C3 plants are not used
    !Config Units = [-]
    CALL getin_p('FQ',fQ)

    !Config Key   = fpseudo
    !Config Desc  = Fraction of electrons at PSI that follow pseudocyclic transport 
    !Config if    = 
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, 0.1, undef, 0.1
    !Config Help  = See Table 2 of Yin et al. (2009) - Values for C3 plants are not used
    !Config Units = [-]
    CALL getin_p('FPSEUDO',fpseudo)

    !Config Key   = kp
    !Config Desc  = Initial carboxylation efficiency of the PEP carboxylase
    !Config if    = 
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, 0.7, undef, 0.7
    !Config Help  = See Table 2 of Yin et al. (2009) 
    !Config Units = [mol m−2 s−1 bar−1]
    CALL getin_p('KP',kp)

    !Config Key   = alpha
    !Config Desc  = Fraction of PSII activity in the bundle sheath
    !Config if    = 
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, 0.1, undef, 0.1
    !Config Help  = See legend of Figure 6 of Yin et al. (2009)
    !Config Units = [-]
    CALL getin_p('ALPHA',alpha)

    !Config Key   = gbs
    !Config Desc  = Bundle-sheath conductance
    !Config if    = 
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, 0.003, undef, 0.003
    !Config Help  = See legend of Figure 6 of Yin et al. (2009)
    !Config Units = [mol m−2 s−1 bar−1]
    CALL getin_p('GBS',gbs)

    !Config Key   = theta
    !Config Desc  = Convexity factor for response of J to irradiance
    !Config if    = 
    !Config Def   = undef, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7
    !Config Help  = See Table 2 of Yin et al. (2009)   
    !Config Units = [−]
    CALL getin_p('THETA',theta)

    !Config Key   = alpha_LL
    !Config Desc  = Conversion efficiency of absorbed light into J at strictly limiting light
    !Config if    = 
    !Config Def   = undef, 0.372, 0.372, 0.372, 0.372, 0.372, 0.372, 0.372, 0.372, 0.372, 0.372, 0.372, 0.372
    !Config Help  = See comment from Yin et al. (2009) after eq. 4
    !Config Units = [mol e− (mol photon)−1]
    CALL getin_p('ALPHA_LL',alpha_LL)

    !Config Key   = STRESS_VCMAX
    !Config Desc  = Stress on vcmax
    !Config if    = OK_SECHIBA or OK_STOMATE
    !Config Def   = 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('STRESS_VCMAX', stress_vcmax)
    
    !Config Key   = STRESS_GS
    !Config Desc  = Stress on gs
    !Config if    = OK_SECHIBA or OK_STOMATE
    !Config Def   = 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('STRESS_GS', stress_gs)
    
    !Config Key   = STRESS_GM
    !Config Desc  = Stress on gm
    !Config if    = OK_SECHIBA or OK_STOMATE
    !Config Def   = 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('STRESS_GM', stress_gm)

    !Config Key   = EXT_COEFF
    !Config Desc  = extinction coefficient of the Monsi&Seaki relationship (1953)
    !Config if    = OK_SECHIBA or OK_STOMATE
    !Config Def   = .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('EXT_COEFF',ext_coeff)

    !Config Key   = EXT_COEFF_VEGETFRAC
    !Config Desc  = extinction coefficient used for the calculation of the bare soil fraction 
    !Config if    = OK_SECHIBA or OK_STOMATE
    !Config Def   = 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('EXT_COEFF_VEGETFRAC',ext_coeff_vegetfrac)

    !
    ! Water-hydrology - sechiba
    !

    !Config Key   = HYDROL_HUMCSTE
    !Config Desc  = Parameter to describe the shape of the structural root profile
    !Config Def   = humcste_ref2m or humcste_ref4m depending on zmaxh
    !Config if    = OK_SECHIBA
    !Config Help  = See module constantes_mtc for different default values
    !Config Units = [-]
    CALL getin_p('HYDROL_HUMCSTE',humcste)

    !Config Key   = MAX_ROOT_DEPTH
    !Config Desc  = Maximum depth of the root profile
    !Config Def   = Maximum depth of the root profile irrespective of the active layer thickness
    !Config if    = OK_SECHIBA
    !Config Help  = See module constantes_mtc for different default values
    !Config Units = [m]
    CALL getin_p('MAX_ROOT_DEPTH',max_root_depth)

    !
    ! Soil - vegetation
    !

    !Config Key   = PREF_SOIL_VEG
    !Config Desc  = The soil tile number for each vegetation
    !Config if    = OK_SECHIBA or OK_STOMATE
    !Config Def   = 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3
    !Config Help  = Gives the number of the soil tile on which we will
    !Config         put each vegetation. This allows to divide the hydrological column
    !Config Units = [-]        
    CALL getin_p('PREF_SOIL_VEG',pref_soil_veg)

    !Config Key   = MAINT_RESP_SLOPE_C 
    !Config Desc  = slope of maintenance respiration coefficient (1/K), constant c of aT^2+bT+c , tabulated 
    !Config if    = OK_STOMATE 
    !Config Def   = undef, 0.44, 0.6, 0.2, 0.6, 0.3, 0.3, 0.4, 0.4, 0.5, 0.5, 0.5, 0.7
    !Config Help  = 
    !Config Units = [-] 
    CALL getin_p('MAINT_RESP_SLOPE_C',maint_resp_slope_c)  
    
    !Config Key   = MAINT_RESP_SLOPE_B 
    !Config Desc  = slope of maintenance respiration coefficient (1/K), constant b of aT^2+bT+c , tabulated 
    !Config if    = OK_STOMATE 
    !Config Def   = undef, -0.01, -0.01, -0.04, -0.02, -0.02, -0.01, -0.01, -0.01, -0.01, -0.0, -0.0, -0.01
    !Config Help  = 
    !Config Units = [-] 
    CALL getin_p('MAINT_RESP_SLOPE_B',maint_resp_slope_b) 
    
    !Config Key   = MAINT_RESP_SLOPE_A 
    !Config Desc  = slope of maintenance respiration coefficient (1/K), constant a of aT^2+bT+c , tabulated 
    !Config if    = OK_STOMATE 
    !Config Def   = undef, .0, .0, .0, .0, .0, .0, .0, .0, .0, .0, .0, .0     
    !Config Help  = 
    !Config Units = [-] 
    CALL getin_p('MAINT_RESP_SLOPE_A',maint_resp_slope_a)

    
    !
    ! Vegetation - Age classes
    !
    
    !Config Key   = NVMAP
    !Config Desc  = The number of PFTs if we ignore age classes.  
    !Config if    = OK_SECHIBA or OK_STOMATE
    !Config Def   = nvm
    !Config Help  = Gives the total number of PFTs ignoring age classes. 
    !               If nagec equals to 1, nvmap is just nvm.
    !Config Units = [-]  
    nvmap=nvm
    CALL getin_p('NVMAP',nvmap)
    IF (printlev>=1) WRITE(numout,*)'the number of pfts for nvmap used by the model: ',nvmap
    IF(nagec > 1 .AND. nvmap == nvm)THEN
       WRITE(numout,*) 'WARNING: The number of age classes is greater than one, but'
       WRITE(numout,*) '         the input file indicates that none of the PFTs have age classes.'
       WRITE(numout,*) '         You should change either nagec or nvmap.'
    ENDIF
    
    !Config Key   = AGEC_GROUP
    !Config Desc  = The species group that each PFT belongs to.  
    !Config if    = OK_SECHIBA or OK_STOMATE
    !Config Def   = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13
    !Config Help  = The species group that each PFT belongs to.
    !Config         A single species/MTC can be represented by a fixed numeber 
    !               of age classes. All age classes of the same species make up the
    !               species group. Within a species group it is assumed that the
    !               age classes are sorted from young to old. Note that it was chosen
    !               to represent the age classes by a diameter threshold rather than
    !               an age threshold. Diameter was picked because this relates better
    !               to height and stand structure and thus the biophysics of the stand
    !               than age (in different climate zone the canopy of stands with the
    !               same age is thought to differ much more than of stands with the
    !               same diameter).
    !Config Units = [-]   
    DO ivm=1,nvm
       agec_group(ivm)=ivm
    ENDDO
    CALL getin_p('AGEC_GROUP',agec_group)
   
    
  END SUBROUTINE config_pft_parameters


!! ================================================================================================================================
!! SUBROUTINE   : config_sechiba_pft_parameters
!!
!>\BRIEF        This subroutine will read the imposed values for the sechiba pft
!! parameters. It is not called if IMPOSE_PARAM is set to NO. 
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE config_sechiba_pft_parameters()

    IMPLICIT NONE

    !! 0. Variables and parameters declaration

    !! 0.1 Input variables

    !! 0.4 Local variable

    !_ ================================================================================================================================ 

    !
    ! Evapotranspiration -  sechiba
    !

    !Config Key   = RSTRUCT_CONST
    !Config Desc  = Structural resistance 
    !Config if    = OK_SECHIBA
    !Config Def   = 0.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0,  2.5,  2.0,  2.0,  2.0
    !Config Help  =
    !Config Units = [s/m]
    CALL getin_p('RSTRUCT_CONST',rstruct_const)

    !Config Key   = KZERO
    !Config Desc  = A vegetation dependent constant used in the calculation of the surface resistance.
    !Config if    = OK_SECHIBA
    !Config Def   = 0.0, 12.E-5, 12.E-5, 12.e-5, 12.e-5, 25.e-5, 12.e-5,25.e-5, 25.e-5, 30.e-5, 30.e-5, 30.e-5, 30.e-5 
    !Config Help  =
    !Config Units = [kg/m^2/s]
    CALL getin_p('KZERO',kzero)

    !Config Key   = ALPHA_WATSTRESS_PFT
    !Config Desc  = PFT dependent parameter to decide exponential function 
    !Config if    = OK_SECHIBA
    !Config Def   = 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    !1.0
    !Config Help  =
    !Config Units = [kg/m^2/s]
    CALL getin_p('ALPHA_WATSTRESS_PFT',alpha_watstress_pft)

    !Config Key   = RVEG_PFT
    !Config Desc  = Artificial parameter to increase or decrease canopy resistance.
    !Config if    = OK_SECHIBA
    !Config Def   = 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.
    !Config Help  = This parameter is set by PFT.
    !Config Units = [-]
    CALL getin_p('RVEG_PFT',rveg_pft)    

    !
    ! Water-hydrology - sechiba
    !

    !Config Key   = WMAX_VEG
    !Config Desc  = Maximum field capacity for each of the vegetations (Temporary): max quantity of water
    !Config if    = OK_SECHIBA
    !Config Def   = 150., 150., 150., 150., 150., 150., 150.,150., 150., 150., 150., 150., 150.
    !Config Help  =
    !Config Units = [kg/m^3]
    CALL getin_p('WMAX_VEG',wmax_veg)

    !Config Key   = PERCENT_THROUGHFALL_PFT
    !Config Desc  = Percent by PFT of precip that is not intercepted by the canopy. Default value depend on run mode.
    !Config if    = OK_SECHIBA
    !Config Def   = Case offline [0. 0. 0....] else [30. 30. 30.....]
    !Config Help  = During one rainfall event, PERCENT_THROUGHFALL_PFT% of the incident rainfall
    !Config         will get directly to the ground without being intercepted, for each PFT.
    !Config Units = [%]
    CALL getin_p('PERCENT_THROUGHFALL_PFT',throughfall_by_pft)
    throughfall_by_pft(:) = throughfall_by_pft(:) / 100. 


    !Config Key   = RSOIL_SCALE
    !Config Desc  = Applay a sclaing factor to the soil resistance to evaporate, value per PFT
    !Config if    = 
    !Config Def   =  1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1.
    !Config Help  = 
    !Config         
    !Config Units = [-]
    CALL getin_p('RSOIL_SCALE', rsoil_scale)

    !! Check that all RSOIL_SCALE are positive
    IF ( ANY(rsoil_scale .LT. 0.0) ) THEN
       WRITE(numout,*) 'RSOIL_SCALE=', rsoil_scale
       CALL ipslerr_p(3,'config_sechiba_pft_parameters','One or several values for RSOIL_SCALE are negative',&
            'Change RSOIL_SCALE in run.def or orchidee.def to zero or higher','')
    END IF
    

    !
    ! Albedo - sechiba
    !

    !Config Key   = SNOWA_AGED_VIS
    !Config Desc  = Minimum snow albedo value for each vegetation type after aging (dirty old snow), visible albedo
    !Config if    = OK_SECHIBA
    !Config Def   = 0.74, 0.0, 0.0, 0.08, 0.24, 0.07, 0.18, 0.18, 0.33, 0.57, 0.57, 0.57, 0.57
    !Config Help  = Values optimized for ORCHIDEE2.0
    !Config Units = [-]
    CALL getin_p('SNOWA_AGED_VIS',snowa_aged_vis)

    !Config Key   = SNOWA_AGED_NIR
    !Config Desc  = Minimum snow albedo value for each vegetation type after aging (dirty old snow), near infrared albedo
    !Config if    = OK_SECHIBA
    !Config Def   = 0.50, 0.0, 0.0, 0.10, 0.37, 0.08, 0.16, 0.17, 0.27, 0.44, 0.44, 0.44, 0.44  
    !Config Help  = Values optimized for ORCHIDEE2.0
    !Config Units = [-]
    CALL getin_p('SNOWA_AGED_NIR',snowa_aged_nir)

    !Config Key   = SNOWA_DEC_VIS
    !Config Desc  = Decay rate of snow albedo value for each vegetation type as it will be used in condveg_snow, visible albedo
    !Config if    = OK_SECHIBA
    !Config Def   = 0.21, 0.0, 0.0, 0.14, 0.08, 0.17, 0.05, 0.06, 0.09, 0.15, 0.15, 0.15, 0.15 
    !Config Help  = Values optimized for ORCHIDEE2.0
    !Config Units = [-]
    CALL getin_p('SNOWA_DEC_VIS',snowa_dec_vis)

    !Config Key   = SNOWA_DEC_NIR
    !Config Desc  = Decay rate of snow albedo value for each vegetation type as it will be used in condveg_snow, near infrared albedo
    !Config if    = OK_SECHIBA
    !Config Def   = 0.13, 0.0, 0.0, 0.10, 0.10, 0.16, 0.04, 0.07, 0.08, 0.12, 0.12, 0.12, 0.12
    !Config Help  = Values optimized for ORCHIDEE2.0
    !Config Units = [-]
    CALL getin_p('SNOWA_DEC_NIR',snowa_dec_nir)

    !Config Key   = ALB_LEAF_VIS
    !Config Desc  = leaf albedo of vegetation type, visible albedo
    !Config if    = OK_SECHIBA
    !Config Def   = 0.00, 0.04, 0.04, 0.04, 0.04, 0.03, 0.03, 0.03, 0.03, 0.06, 0.06, 0.06, 0.06
    !Config Help  = Values optimized for ORCHIDEE2.0
    !Config Units = [-]
    CALL getin_p('ALB_LEAF_VIS',alb_leaf_vis)

    !Config Key   = ALB_LEAF_NIR
    !Config Desc  = leaf albedo of vegetation type, near infrared albedo
    !Config if    = OK_SECHIBA
    !Config Def   = 0.00, 0.23, 0.18, 0.18, 0.20, 0.24, 0.15, 0.26, 0.20, 0.24, 0.27, 0.28, 0.26
    !Config Help  = Values optimized for ORCHIDEE2.0
    !Config Units = [-]
    CALL getin_p('ALB_LEAF_NIR',alb_leaf_nir)

    !Config Key  = LEAF_SSA_VIS
    !Config Desc = Leaf_single_scattering_albedo_vis values
    !Config If   = ALBEDO_TYPE is Pinty
    !Config Def  = 0.17192, 0.12560, 0.16230, 0.13838, 0.13202, 0.14720, 0.14680, 0.14415, 0.15485, 0.17544, 0.17384, 0.17302, 0.17116 
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('LEAF_SSA_VIS',leaf_ssa(:,ivis))

    !Config Key  = LEAF_SSA_NIR
    !Config Desc = Leaf_single_scattering_albedo_nir values
    !Config If   = ALBEDO_TYPE is Pinty
    !Config Def  = 0.70253, 0.68189, 0.69684, 0.68778, 0.68356, 0.69533, 0.69520, 0.69195, 0.69180, 0.71236, 0.71904, 0.71220, 0.71190
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('LEAF_SSA_NIR',leaf_ssa(:,inir))
    !
    !Config Key  = LEAF_PSD_VIS
    !Config Desc = Preferred scattering direction values in the visibile spectra
    !Config If   = ALBEDO_TYPE is Pinty
    !Config Def  = 1.00170, 0.96776, 0.99250, 0.97170, 0.97119, 0.98077, 0.97672, 0.97810, 0.98605, 1.00490, 1.00360, 1.00320, 1.00130
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('LEAF_PSD_VIS',leaf_psd(:,ivis))
    !
    !Config Key  = LEAF_PSD_NIR
    !Config Desc =  Preferred scattering direction values in the near infrared spectra
    !Config If   = ALBEDO_TYPE is Pinty
    !Config Def  = 2.00520, 1.95120, 1.98990, 1.97020, 1.95900, 1.98190, 1.98890, 1.97400, 1.97780, 2.02430, 2.03350, 2.02070, 2.02150
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('LEAF_PSD_NIR',leaf_psd(:,inir)) 
    !
    !Config Key  = BGRD_REF_VIS
    !Config Desc = Background reflectance values in the visibile spectra
    !Config If   = ALBEDO_TYPE is Pinty
    !Config Def   = 0.2300000,   0.0866667,   0.0800000,   0.0533333,   0.0700000,   0.0933333,   0.0533333, 0.0833333,   0.0633333,   0.1033330,   0.1566670,   0.1166670,   0.1200000
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('BGRD_REF_VIS',bgd_reflectance(:,ivis)) 
    !
    !Config Key  = BGRD_REF_NIR
    !Config Desc = Background reflectance values in the near infrared spectra
    !Config If   = ALBEDO_TYPE is Pinty
    !Config Def   = 0.4200000,   0.1500000,   0.1300000,   0.0916667,   0.1066670,   0.1650000,   0.0900000, 0.1483330,   0.1066670,   0.1900000,   0.3183330,   0.2200000,   0.2183330
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('BGRD_REF_NIR',bgd_reflectance(:,inir)) 

    !Config Key  = LEAF_TO_SHOOT_CLUMPING
    !Config Desc = The leaf-to-shoot clumping factor
    !Config If   = ALBEDO_TYPE is Pinty
    !Config Def   = un, un, un, un, un, un, un, un, un, un, un, un, un
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('LEAF_TO_SHOOT_CLUMPING',leaf_to_shoot_clumping(:)) 
    !
    !Config Key  = LAI_CORRECTION_FACTOR
    !Config Desc = The correction factor for the LAI for grasslands 
    !              and crops (see note in pft_parameters)
    !Config If   = ALBEDO_TYPE is Pinty
    !Config Def   = un, un, un, un, un, un, un, un, un, un, un, un, un
    !Config Help  =
    !Config Units = [-] 
    CALL getin_p('LAI_CORRECTION_FACTOR',lai_correction_factor(:)) 

    !Config Key  = MIN_LEVEL_SEP
    !Config Desc = The minimum level thickness we use for photosynthesis 
    !Config If   = ALBEDO_TYPE is Pinty
    !Config Def  = un, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1
    !Config Help  =
    !Config Units = [m]
    CALL getin_p('MIN_LEVEL_SEP',min_level_sep(:)) 

    !Config Key  = LAI_TOP
    !Config Desc = Definition, in terms of LAI of the top layer 
    !              (used to calculate one of the resistences of 
    !              vbeta3) to calculate transpiration 
    !Config If   = OK_SECHIBA
    !Config Def  = un, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1
    !Config Help  = Coupling of the canopy to the atmosphere. See ADVANCES IN ECOLOGICAL 
    !               RESEARCH, VOLUME 15. Stomatal Control of Transpiration: Scaling Up 
    !               from Leaf to Region. 1986. P. G. JARVIS and K. G. MCNAUGHTON 
    !Config Units = [m2 m2]
    CALL getin_p('LAI_TOP',lai_top(:)) 

    !
    ! Hydraulic architecture
    !
    !Config Key   = K_ROOT
    !Config Desc  = Fine root specific conductivity
    !Config if    = OK_STOMATE
    !Config Def   = (undef, 7.02, 7.02, 7.02, 7.02, 7.02, 7.02, 7.02, 7.02, 7.02, 7.02, 7.02, 7.02)*1.e-4 
    !Config Help  =
    !Config Units = [m^{3} kg^{-1} s^{-1} MPa^{-1}] 
    CALL getin_p('K_ROOT',k_root) 

    !Config Key   = K_BELOWGROUND
    !Config Desc  = Belowground (roots + soil) specific conductivity used in allocation 
    !Config if    = OK_STOMATE
    !Config Def   = (undef, 7., 7., 7., 7., 7., 7., 7., 7., 42., 42., 42., 42.)*1.e-7 
    !Config Help  =
    !Config Units = [m^{3} kg^{-1} s^{-1} MPa^{-1}] 
    CALL getin_p('K_BELOWGROUND',k_belowground)

    !Config Key   = K_SAP
    !Config Desc  = Sapwood specific conductivity
    !Config if    = OK_STOMATE
    !Config Def   = (undef, 50., 10., 8., 5., 30., 8., 20., 8., undef, undef, undef, undef)*1.e-4
    !Config Help  =
    !Config Units = [m^{2} s^{-1} MPa^{-1}] 
    CALL getin_p('K_SAP',k_sap)

    !Config Key   = PSI_LEAF_MIN
    !Config Desc  = Minimal leaf potential
    !Config if    = OK_STOMATE, 11-LAYERS, FUNCTIONAL ALLOCATION
    !Config Def   = undef, -2.2, -2.2, -2.2, -3.5, -2.2, -2.2, -2.2, -2.2, -2.2, -2.2, -2.2, -2.2
    !Config Help  =
    !Config Units = [MPa] 
    CALL getin_p('PSI_LEAF_MIN',psi_leaf_min)

    !Config Key   = SRL
    !Config Desc  = Specific root length
    !Config if    = OK_STOMATE, 11-LAYERS, FUNCTIONAL ALLOCATION,
    !Config if      HYDRAULIC_ARCHITECTURE
    !Config Def   = undef, 10, 10, 9.2, 9.2, 14, 18.3, 18.3, 18.3, undef, undef, undef, undef 
    !Config Help  = Specific root length for the calculations of soil to root resistance.
    !Config Units = [m g^(-1)] 
    CALL getin_p('SRL',srl)

    !Config Key   = R_FROOT
    !Config Desc  = Fine root radius 
    !Config if    = OK_STOMATE, 11-LAYERS, FUNCTIONAL ALLOCATION,
    !Config if      HYDRAULIC_ARCHITECTURE
    !Config Def   = undef,  0.29E-3, 0.29E-3,  0.29E-3,  0.29E-3, 0.29E-3,  0.24E-3, 0.21E-3,  0.21E-3, undef, undef,  undef,  undef 
    !Config Help  = Fine root radius for the calculations of soil to root resistance
    !Config Units = [m] 
    CALL getin_p('R_FROOT',r_froot)

    !Config Key   = PSI_ROOT_MIN
    !Config Desc  = Minimum root water potential 
    !Config if    = OK_STOMATE, 11-LAYERS, FUNCTIONAL ALLOCATION,
    !Config if      HYDRAULIC_ARCHITECTURE
    !Config Def   = undef, -4, -4, -4, -4, -4, -4, -4, -4, undef, undef, undef, undef 
    !Config Help  = Minimum root water potential for the calculations of 
    !               soil to root resistance.
    !Config Units = [MPa] 
    CALL getin_p('PSI_ROOT_MIN',psi_root_min)
    
    !
    ! Laieff - .NOT. ok_stomate
    !
    !Config Key   = CROWN_TO_HEIGHT
    !Config Desc  = Ratio between tree height and the vertical crown diameter.
    !Config If    = OK_STOMATE 
    !Config Def   = undef, 0.6, 0.6, 0.6, 0.6, 0.6, 0.8, 0.8, 0.8, 0., 0., 0.,
    !0.
    !Config Help  = Ratio between tree height and the vertical crown diameter.
    !If this value is changed check beforehand that the crown diameter will
    !never exceed the tree height.
    !Config Units = [-]  
    CALL getin_p('CROWN_TO_HEIGHT',crown_to_height)

    !Config Key   = CROWN_VERTOHOR_DIA
    !Config Desc  = Ratio between the vertical and horizontal crown diameter height.
    !diameter 
    !Config If    = OK_STOMATE 
    !Config Def   = undef, 1.0, 1.0, 0.66, 1.0, 1.0, 0.66, 1.0, 1.0, 1.0, 1.0, 1.0,
    !1.0
    !Config Help  = Ratio between the vertical and horizontal crown diameter
    !height, so indirectly the horizontal crown diameter also depends on crown
    !diameter 
    !Config Units = [-]  
    CALL getin_p('CROWN_VERTOHOR_DIA',crown_vertohor_dia)

    !Config Key   = PIPE_DENSITY 
    !Config Desc  = 
    !Config if    = 
    !Config Def   = undef, 3.e5, 3.e5, 2.e5, 3.e5, 3.e5, 2.e5, 3.e5, 2.e5, 2.e5, 2.e5, 2.e5, 2.e5
    !Config Help  = 
    !Config Units = 
    CALL getin_p("PIPE_DENSITY",pipe_density)
    !
    !Config Key   = TREE_FF
    !Config Desc  = Tree form factor reducing the volume of a cylinder
    !               to the real volume of the tree shape (including the 
    !               branches)
    !Config If    = OK_STOMATE 
    !Config Def   = undef, 0.6, 0.6, 0.6, 0.6, 0.6, 0.8, 0.8, 0.8, 0., 0., 0., 0.
    !Config Help  = 
    !Config Units = [-]  
    CALL getin_p('TREE_FF',tree_ff)
    !
    !Config Key   = MIN_PIPE_TUNE2
    !Config Desc  = Initial value when PIPE_TUNE2 is calculated as a function of precipitation
    !Config If    = OK_STOMATE, HACK_DYNAMIC_HEIGHT 
    !Config Def   = undef, 7, 7, 7, 7, 7, 7, 7, 7, 0., 0., 0., 0.   
    !Config Help  = 
    !Config Units = [m]
    CALL getin_p('MIN_PIPE_TUNE2',min_pipe_tune2)
    !
    !Config Key   = PIPE_TUNE3
    !Config Desc  = height=pipe_tune2 * diameter**pipe_tune3
    !Config If    = OK_STOMATE 
    !Config Def   = undef, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0., 0., 0., 0.   
    !Config Help  = 
    !Config Units = [-]    
    CALL getin_p('PIPE_TUNE3',pipe_tune3)
    !
    !Config Key   = PIPE_TUNE4
    !Config Desc  = needed for stem diameter
    !Config If    = OK_STOMATE 
    !Config Def   = undef, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0., 0., 0., 0.
    !Config Help  = 
    !Config Units = [-]  
    CALL getin_p('PIPE_TUNE4',pipe_tune4)
    !
    !Config Key   = PIPE_K1 
    !Config Desc  = 
    !Config If    = OK_STOMATE 
    !Config Def   = undef, 8.e3, 8.e3, 8.e3, 8.e3, 8.e3, 8.e3, 8.e3, 8.e3, 0., 0., 0., 0. 
    !Config Help  = 
    !Config Units = [-]   
    CALL getin_p('PIPE_K1',pipe_k1)
    !
    !Config Key   = SLA
    !Config Desc  = specif leaf area 
    !Config if    = OK_STOMATE
    !Config Def   = 1.5E-2, 1.53E-2, 2.6E-2, 9.26E-3, 2E-2, 2.6E-2, 9.26E-3, 2.6E-2, 1.9E-2, 2.6E-2, 2.6E-2, 2.6E-2, 2.6E-2
    !Config Help  =
    !Config Units = [m^2/gC]
    CALL getin_p('SLA',sla)
    !
    !Config Key   = SLAINIT
    !Config Desc  = initial specif leaf area at (ie at bottom of canopy eq. lai=0) 
    !Config if    = OK_STOMATE
    !Config Def   = 2.6E-2, 2.6E-2, 4.4E-2, 1.4E-2, 3.0E-2, 3.9E-2, 1.3E-2, 3.7E-2, 2.4E-2, 3.1E-2, 3.1E-2, 3.9E-2, 3.9E-2
    !Config Help  =
    !Config Units = [m^2/gC]
    CALL getin_p('SLAINIT',slainit)
    !
    !Config Key   = LAI_TO_HEIGHT
    !Config Desc  = Convertion factor from lai to vegetation height for grasses and crops
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, 0.2, 0.5, 0.2, 0.5
    !Config Help  = 
    !Config Units = [m m2 m-2] 
    CALL getin_p('LAI_TO_HEIGHT',lai_to_height)

    IF ( ok_bvoc ) THEN
       !
       ! BVOC
       !

       !Config Key   = ISO_ACTIVITY
       !Config Desc  = Biogenic activity for each age class : isoprene
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0.5, 1.5, 1.5, 0.5
       !Config Help  =
       !Config Units = [-]
       iso_activity = (/0.5, 1.5, 1.5, 0.5/) 
       CALL getin_p('ISO_ACTIVITY',iso_activity)

       !Config Key   = METHANOL_ACTIVITY
       !Config Desc  = Isoprene emission factor for each age class : methanol
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 1., 1., 0.5, 0.5
       !Config Help  =
       !Config Units = [-]
       methanol_activity = (/1., 1., 0.5, 0.5/)
       CALL getin_p('METHANOL_ACTIVITY',methanol_activity)

       !Config Key   = EM_FACTOR_ISOPRENE
       !Config Desc  = Isoprene emission factor
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0., 24., 24., 8., 16., 45., 8., 18., 0.5, 12., 18., 5., 5.
       !Config Help  =
       !Config Units = [ugC/g/h] 
       CALL getin_p('EM_FACTOR_ISOPRENE',em_factor_isoprene)

       !Config Key   = EM_FACTOR_MONOTERPENE
       !Config Desc  = Monoterpene emission factor 
       !Config if    = CHEMISTRY_BVOC 
       !Config Def   = 0., 2.0, 2.0, 1.8, 1.4, 1.6, 1.8, 1.4, 1.8, 0.8, 0.8,  0.22, 0.22
       !Config Help  =
       !Config Units = [ugC/g/h] 
       CALL getin_p('EM_FACTOR_MONOTERPENE',em_factor_monoterpene)

       !Config Key   = C_LDF_MONO 
       !Config Desc  = Monoterpenes fraction dependancy to light
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0.6
       !Config Help  =
       !Config Units = []
       CALL getin_p('C_LDF_MONO',LDF_mono)

       !Config Key   = C_LDF_SESQ 
       !Config Desc  = Sesquiterpenes fraction dependancy to light
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0.5
       !Config Help  =
       !Config Units = []
       CALL getin_p('C_LDF_SESQ',LDF_sesq)

       !Config Key   = C_LDF_METH 
       !Config Desc  = Methanol fraction dependancy to light
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0.8
       !Config Help  =
       !Config Units = []
       CALL getin_p('C_LDF_METH',LDF_meth)

       !Config Key   = C_LDF_ACET 
       !Config Desc  = Acetone fraction dependancy to light
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0.2
       !Config Help  =
       !Config Units = []
       CALL getin_p('C_LDF_ACET',LDF_acet)

       !Config Key   = EM_FACTOR_APINENE 
       !Config Desc  = Alfa pinene  emission factor 
       !Config if    = CHEMISTRY_BVOC 
       !Config Def   = 0., 1.35, 1.35, 0.85, 0.95, 0.75, 0.85, 0.60, 1.98, 0.30, 0.30, 0.09, 0.09
       !Config Help  =
       !Config Units = [ugC/g/h] 
       CALL getin_p('EM_FACTOR_APINENE',em_factor_apinene)

       !Config Key   = EM_FACTOR_BPINENE
       !Config Desc  = Beta pinene  emission factor
       !Config if    = CHEMISTRY_BVOC 
       !Config Def   = 0., 0.30, 0.30, 0.35, 0.25, 0.20, 0.35, 0.12, 0.45, 0.16, 0.12, 0.05, 0.05
       !Config Help  =
       !Config Units = [ugC/g/h] 
       CALL getin_p('EM_FACTOR_BPINENE',em_factor_bpinene)

       !Config Key   = EM_FACTOR_LIMONENE
       !Config Desc  = Limonene  emission factor
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0., 0.25, 0.25, 0.20, 0.25, 0.14, 0.20, 0.135, 0.11, 0.19, 0.42, 0.03, 0.03
       !Config Help  =
       !Config Units = [ugC/g/h] 
       CALL getin_p('EM_FACTOR_LIMONENE',em_factor_limonene)

       !Config Key   = EM_FACTOR_MYRCENE
       !Config Desc  = Myrcene  emission factor
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0., 0.20, 0.20, 0.12, 0.11, 0.065, 0.12, 0.036, 0.075, 0.08,  0.085, 0.015, 0.015
       !Config Help  =
       !Config Units = [ugC/g/h] 
       CALL getin_p('EM_FACTOR_MYRCENE',em_factor_myrcene)

       !Config Key   = EM_FACTOR_SABINENE
       !Config Desc  = Sabinene  emission factor
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0., 0.20, 0.20, 0.12, 0.17, 0.70, 0.12, 0.50, 0.09, 0.085, 0.075, 0.02, 0.02
       !Config Help  =
       !Config Units = [ugC/g/h] 
       CALL getin_p('EM_FACTOR_SABINENE',em_factor_sabinene)

       !Config Key   = EM_FACTOR_CAMPHENE 
       !Config Desc  = Camphene  emission factor 
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0., 0.15, 0.15, 0.10, 0.10, 0.01, 0.10, 0.01, 0.07, 0.07, 0.08, 0.01, 0.01
       !Config Help  =
       !Config Units = [ugC/g/h] 
       CALL getin_p('EM_FACTOR_CAMPHENE',em_factor_camphene)

       !Config Key   = EM_FACTOR_3CARENE 
       !Config Desc  = 3-Carene  emission factor
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0., 0.13, 0.13, 0.42, 0.02, 0.055, 0.42,0.025, 0.125, 0.085, 0.085, 0.065, 0.065
       !Config Help  =
       !Config Units = [ugC/g/h] 
       CALL getin_p('EM_FACTOR_3CARENE',em_factor_3carene)

       !Config Key   = EM_FACTOR_TBOCIMENE
       !Config Desc  = T-beta-ocimene  emission factor
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0., 0.25, 0.25, 0.13, 0.09, 0.26, 0.13, 0.20, 0.085, 0.18, 0.18, 0.01, 0.01
       !Config Help  =
       !Config Units = [ugC/g/h] 
       CALL getin_p('EM_FACTOR_TBOCIMENE', em_factor_tbocimene)

       !Config Key   = EM_FACTOR_OTHERMONOT
       !Config Desc  = Other monoterpenes  emission factor
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0., 0.17, 0.17, 0.11, 0.11, 0.125, 0.11, 0.274, 0.01, 0.15, 0.155, 0.035, 0.035
       !Config Help  =
       !Config Units = [ugC/g/h] 
       CALL getin_p('EM_FACTOR_OTHERMONOT',em_factor_othermonot)

       !Config Key   = EM_FACTOR_SESQUITERP 
       !Config Desc  = Sesquiterpenes  emission factor 
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0., 0.45, 0.45, 0.13, 0.3, 0.36, 0.15, 0.3, 0.25, 0.6, 0.6, 0.08, 0.08
       !Config Help  =
       !Config Units = [ugC/g/h] 
       CALL getin_p('EM_FACTOR_SESQUITERP',em_factor_sesquiterp)



       !Config Key   = C_BETA_MONO 
       !Config Desc  = Monoterpenes temperature dependency coefficient
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0.1
       !Config Help  =
       !Config Units = []
       CALL getin_p('C_BETA_MONO',beta_mono)

       !Config Key   = C_BETA_SESQ 
       !Config Desc  = Sesquiterpenes temperature dependency coefficient
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0.17
       !Config Help  =
       !Config Units = []
       CALL getin_p('C_BETA_SESQ',beta_sesq)

       !Config Key   = C_BETA_METH 
       !Config Desc  = Methanol temperature dependency coefficient
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0.08
       !Config Help  =
       !Config Units = []
       CALL getin_p('C_BETA_METH',beta_meth)

       !Config Key   = C_BETA_ACET 
       !Config Desc  = Acetone temperature dependency coefficient
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0.1
       !Config Help  =
       !Config Units = []
       CALL getin_p('C_BETA_ACET',beta_acet)

       !Config Key   = C_BETA_OXYVOC 
       !Config Desc  = Other oxygenated BVOC temperature dependency coefficient
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0.13
       !Config Help  =
       !Config Units = []
       CALL getin_p('C_BETA_OXYVOC',beta_oxyVOC)

       !Config Key   = EM_FACTOR_ORVOC
       !Config Desc  = ORVOC emissions factor 
       !Config if    = CHEMISTRY_BVOC 
       !Config Def   = 0., 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5
       !Config Help  =
       !Config Units = [ugC/g/h]  
       CALL getin_p('EM_FACTOR_ORVOC',em_factor_ORVOC)

       !Config Key   = EM_FACTOR_OVOC
       !Config Desc  = OVOC emissions factor
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0., 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5
       !Config Help  =
       !Config Units = [ugC/g/h]        
       CALL getin_p('EM_FACTOR_OVOC',em_factor_OVOC)

       !Config Key   = EM_FACTOR_MBO
       !Config Desc  = MBO emissions factor 
       !Config if    = CHEMISTRY_BVOC 
       !Config Def   = 0., 2.e-5, 2.e-5, 1.4, 2.e-5, 2.e-5, 0.14, 2.e-5, 2.e-5, 2.e-5, 2.e-5, 2.e-5, 2.e-5
       !Config Help  =
       !Config Units = [ugC/g/h]  
       CALL getin_p('EM_FACTOR_MBO',em_factor_MBO)

       !Config Key   = EM_FACTOR_METHANOL
       !Config Desc  = Methanol emissions factor 
       !Config if    = CHEMISTRY_BVOC 
       !Config Def   = 0., 0.8, 0.8, 1.8, 0.9, 1.9, 1.8, 1.8, 1.8, 0.7, 0.9, 2., 2.
       !Config Help  =
       !Config Units = [ugC/g/h]  
       CALL getin_p('EM_FACTOR_METHANOL',em_factor_methanol)

       !Config Key   = EM_FACTOR_ACETONE
       !Config Desc  = Acetone emissions factor
       !Config if    = CHEMISTRY_BVOC 
       !Config Def   = 0., 0.25, 0.25, 0.3, 0.2, 0.33, 0.3, 0.25, 0.25, 0.2, 0.2, 0.08, 0.08
       !Config Help  =
       !Config Units = [ugC/g/h]     
       CALL getin_p('EM_FACTOR_ACETONE',em_factor_acetone)

       !Config Key   = EM_FACTOR_ACETAL
       !Config Desc  = Acetaldehyde emissions factor 
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0., 0.2, 0.2, 0.2, 0.2, 0.25, 0.25, 0.16, 0.16, 0.12, 0.12, 0.035, 0.02
       !Config Help  =
       !Config Units = [ugC/g/h]  
       CALL getin_p('EM_FACTOR_ACETAL',em_factor_acetal)

       !Config Key   = EM_FACTOR_FORMAL
       !Config Desc  = Formaldehyde emissions factor
       !Config if    = CHEMISTRY_BVOC 
       !Config Def   = 0., 0.04, 0.04, 0.08, 0.04, 0.04, 0.04, 0.04, 0.04, 0.025, 0.025, 0.013, 0.013
       !Config Help  = 
       !Config Units = [ugC/g/h]  
       CALL getin_p('EM_FACTOR_FORMAL',em_factor_formal)

       !Config Key   = EM_FACTOR_ACETIC
       !Config Desc  = Acetic Acid emissions factor
       !Config if    = CHEMISTRY_BVOC 
       !Config Def   = 0., 0.025, 0.025,0.025,0.022,0.08,0.025,0.022,0.013,0.012,0.012,0.008,0.008
       !Config Help  =
       !Config Units = [ugC/g/h]  
       CALL getin_p('EM_FACTOR_ACETIC',em_factor_acetic)

       !Config Key   = EM_FACTOR_FORMIC
       !Config Desc  = Formic Acid emissions factor
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0., 0.015, 0.015, 0.02, 0.02, 0.025, 0.025, 0.015, 0.015,0.010,0.010,0.008,0.008
       !Config Help  =
       !Config Units = [ugC/g/h]  
       CALL getin_p('EM_FACTOR_FORMIC',em_factor_formic)

       !Config Key   = EM_FACTOR_DMS
       !Config Desc  = DMS emissions factor
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0., 0.037, 0.037, 0.027, 0.027, 0.027, 0.027, 0.027, 0.027, 0.027, 0.027, 0.027, 0.027
       !Config Help  =
       !Config Units = [ug/g/h]
       CALL getin_p('EM_FACTOR_DMS',em_factor_DMS)

       !Config Key   = EM_FACTOR_H2S
       !Config Desc  = H2S emissions factor
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0., 0.0095, 0.0095, 0.0064, 0.0064, 0.0064, 0.0064, 0.0064, 0.0064, 0.0064, 0.0064, 0.0064, 0.0064
       !Config Help  =
       !Config Units = [ug/g/h]
       CALL getin_p('EM_FACTOR_H2S',em_factor_H2S)

       !Config Key   = EM_FACTOR_NO_WET
       !Config Desc  = NOx emissions factor wet soil emissions and exponential dependancy factor 
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0., 2.6, 0.06, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.36, 0.36, 0.36, 0.36
       !Config Help  =
       !Config Units = [ngN/m^2/s]
       CALL getin_p('EM_FACTOR_NO_WET',em_factor_no_wet)

       !Config Key   = EM_FACTOR_NO_DRY
       !Config Desc  = NOx emissions factor dry soil emissions and exponential dependancy factor 
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0., 8.60, 0.40, 0.22, 0.22, 0.22, 0.22, 0.22, 0.22, 2.65, 2.65, 2.65, 2.65
       !Config Help  =
       !Config Units = [ngN/m^2/s] 
       CALL getin_p('EM_FACTOR_NO_DRY',em_factor_no_dry)

       !Config Key   = LARCH
       !Config Desc  = Larcher 1991 SAI/LAI ratio
       !Config if    = CHEMISTRY_BVOC 
       !Config Def   = 0., 0.015, 0.015, 0.003, 0.005, 0.005, 0.003, 0.005, 0.003, 0.005, 0.005, 0.008, 0.008
       !Config Help  =
       !Config Units = [-]  
       CALL getin_p('LARCH',Larch)

    ENDIF ! (ok_bvoc)


    !Config Key   = NUE_OPT
    !Config Desc  = Nitrogen use efficiency of Vcmax 
    !Config if    = OK_STOMATE
    !Config Def   = undef,  14.,  30., 20., 33.,  38., 15., 38., 22.,  45.,  45.,  60.,  60.  
    !Config Help  =
    !Config Units = [(mumol[CO2] s-1) (gN[leaf])-1]
    CALL getin_p('NUE_OPT',nue_opt)

    !Config Key   = VMAX_UPTAKE_NH4
    !Config Desc  = Vmax of ammonium uptake by plant roots
    !Config if    = OK_STOMATE
    !Config Def   = undef,  9., 9., 9., 9., 9., 9., 9., 9., 9., 9., 9., 9.  
    !Config Help  =
    !Config Units = umol (g DryWeight_root)-1 h-1
    CALL getin_p('VMAX_UPTAKE_NH4',vmax_uptake(:,iammonium))

    !Config Key   = VMAX_UPTAKE_NO3
    !Config Desc  = Vmax of nitrate uptake by plant roots
    !Config if    = OK_STOMATE
    !Config Def   = undef,  9., 9., 9., 9., 9., 9., 9., 9., 9., 9., 9., 9.  
    !Config Help  =
    !Config Units = umol (g DryWeight_root)-1 h-1
    CALL getin_p('VMAX_UPTAKE_NO3',vmax_uptake(:,initrate))
    
    !Config Key   = CN_LEAF_MIN
    !Config Desc  = minimum CN ratio of leaves  
    !Config if    = OK_STOMATE
    !Config Def   = undef, 16., 16., 28., 16., 16., 28., 16., 16., 16., 16., 16., 16. 
    !Config Help  =
    !Config Units = [gC/gN] 
    CALL getin_p("CN_LEAF_MIN", cn_leaf_min)
    
    !Config Key   = CN_LEAF_MAX
    !Config Desc  = maximum CN ratio of leaves  
    !Config if    = OK_STOMATE
    !Config Def   = undef, 45., 45., 75., 45., 45., 75., 45., 45., 45., 45., 45., 45. 
    !Config Help  =
    !Config Units = [gC/gN] 
    CALL getin_p("CN_LEAF_MAX", cn_leaf_max)
    
    !Config Key   = CN_LEAF_INIT
    !Config Desc  = 
    !Config if    = 
    !Config Def   = undef, 25.,  25.,  41.7,  25.,  25.,  43., 25.,  25.,  25.,  25.,  25.,  25.
    !Config Help  = Comes from Sitch et al 2003 (https://doi.org/10.1046/j.1365-2486.2003.00569.x),
    !Config         although the defaults have changed for an unknown reason.  In Sitch et al,
    !Config         the leaf ratio is 29.
    !Config Units = 
    CALL getin_p("CN_LEAF_INIT",cn_leaf_init)

    !Config Key   = EXT_COEFF_N
    !Config Desc  = Extinction coefficient of the leaf N content profile within the canopy
    !Config if    = OK_STOMATE
    !Config Def   =  0.15, 0.15, 0.15,0.15,0.15, 0.15,0.15,0.15,0.15, 0.15, 0.15, 0.15, 0.15
    !Config Help  =
    !Config Units = [(m2[ground]) (m-2[leaf])]
    CALL getin_p('EXT_COEFF_N',ext_coeff_N)

    !Config Key   = MOSS_FRAC
    !Config Desc  = Moss fraction for each PFT
    !Config if    = OK_MOSS
    !Config Def   = 0., 0., 0., 0.,  0.,  0.,  1., 1.,  1.,  0.,  0.,  0., 0.
    !Config Help  = Default value is one for boreal pfts
    !Config Units = [0-1]
    CALL getin_p("MOSS_FRAC",moss_frac)

    !! Below is a stomate variable but it can be a diagnotic variable if running without ok_stomate
    !! or they are used by READ_LAI=1 (option to bypass stomate)

    !Config Key   = DIA_THRESH_INV
    !Config Desc  = Diameter below which trees are not measured in a forest inventory
    !Config if    = OK_STOMATE
    !Config Def   = undef, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, undef, undef, undef,undef
    !Config Help  = Forest inventories do not sample the smallest trees. The cut-off  diameter depends 
    !               on the country but typically varies between 7 and 10 cm. If ORCHIDEE output is compared
    !               against forest inventory data. HEIGHT_INV and DIA_INV should be calculated with the correct
    !              threshold diameter
    !Config Units = [m]
    CALL getin_p("DIA_THRESH_INV",dia_thresh_inv)

    !Config Key   = SENESCENCE_TYPE
    !Config Desc  = type of senescence, tabulated
    !Config if    = OK_STOMATE
    !Config Def   = none, none, dry, none, none, cold, none, cold, cold, mixed, mixed, crop, crop 
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('SENESCENCE_TYPE',senescence_type) 

    !Config Key   = LONGEVITY_SAP
    !Config Desc  = sapwood -> heartwood conversion time
    !Config if    = OK_STOMATE
    !Config Def   = undef, 730., 730., 730., 730., 730., 730., 730., 730., undef, undef, undef, undef
    !Config Help  =
    !Config Units = [days]
    CALL getin_p('LONGEVITY_SAP',longevity_sap)

    !Config Key   = LONGEVITY_ROOT
    !Config Desc  = root longevity
    !Config if    = OK_STOMATE
    !Config Def   = undef, 256., 256., 256., 256., 256., 256., 256., 256., 256., 256., 256., 256.
    !Config Help  =
    !Config Units = [days]
    CALL getin_p('LONGEVITY_ROOT',longevity_root)

    !Config Key   = K_LATOSA_MAX
    !Config Desc  = Maximum leaf-to-sapwood area ratio
    !Config if    = OK_STOMATE
    !Config Def   = (undef, 5., 5., 5., 3., 5., 5., 5., 5., undef, undef, undef, undef)*1.e3
    !Config Help  =
    !Config Units = [-] 
    CALL getin_p('K_LATOSA_MAX',k_latosa_max)

    !Config Key   = K_LATOSA_MIN
    !Config Desc  = Minimum leaf-to-sapwood area ratio
    !Config if    = OK_STOMATE
    !Config Def   = (undef, 5., 5., 5., 3., 5., 5., 5., 5., undef, undef, undef, undef)*1.e3
    !Config Help  =
    !Config Units = [-] 
    CALL getin_p('K_LATOSA_MIN',k_latosa_min) 

  END SUBROUTINE config_sechiba_pft_parameters


!! ================================================================================================================================
!! SUBROUTINE   : config_stomate_pft_parameters 
!!
!>\BRIEF         This subroutine will read the imposed values for the stomate pft
!! parameters. It is not called if IMPOSE_PARAM is set to NO.
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE config_stomate_pft_parameters

    IMPLICIT NONE

    !! 0. Variables and parameters declaration

    !! 0.4 Local variable
    INTEGER(i_std)  ::  ivma, ivm             !! indices for number of pfts with and without age classes

    !_ ================================================================================================================================

    !
    ! Vegetation structure
    !
    !Config Key   = AVAILABILITY_FACT 
    !Config Desc  = Calculate dynamic mortality in lpj_gap, pft dependent parameter
    !Config If    = OK_STOMATE 
    !Config Def   = undef, 0.14, 0.14, 0.10, 0.10, 0.10, 0.05, 0.05, 0.05, undef, undef, undef, undef 
    !Config Help  = 
    !Config Units = [-]   
    CALL getin_p('AVAILABILITY_FACT',availability_fact)

    !
    ! Respiration - stomate
    !

    !Config Key   = FRAC_GROWTHRESP
    !Config Desc  = fraction of GPP which is lost as growth respiration
    !Config if    = OK_STOMATE
    !Config Def   = undef, 0.35, 0.35, 0.28, 0.28, 0.28, 0.35, 0.35, 0.35, 0.28, 0.28, 0.28, 0.28
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('FRAC_GROWTHRESP',frac_growthresp) 
    
    !Config Key   = COEFF_MAINT_INIT
    !Config Desc  = maintenance respiration coefficient at 10 deg C
    !Config if    = OK_STOMATE
    !Config Def   = undef, 3.06E-2, 3.06E-2, 6.46E-2, 6.46E-2, 6.46E-2, 6.46E-2, 6.46E-2, 6.46E-2, 6.46E-2, 6.46E-2, 6.46E-2, 6.46E-2
    !Config Help  =
    !Config Units = [gC/gN/day]
    CALL getin_p('COEFF_MAINT_INIT',coeff_maint_init)

    !Config Key   = TREF_MAINT_RESP
    !Config Desc  = maintenance respiration Temperature coefficient
    !Config if    = OK_STOMATE
    !Config Def   =   &  undef, 56.02, 56.02, 56.02, 56.02, 56.02, 56.02, 56.02, 56.02, 56.02, 56.02, 56.02, 56.02   
    !Config Help  =
    !Config Units = [degC]
    CALL getin_p('TREF_MAINT_RESP',tref_maint_resp)

    !Config Key   = TMIN_MAINT_RESP
    !Config Desc  = maintenance respiration Temperature coefficient
    !Config if    = OK_STOMATE
    !Config Def   =  undef, 46.02, 46.02, 46.02, 46.02, 46.02, 46.02, 46.02, 46.02, 46.02, 46.02, 46.02, 46.02   
    !Config Help  =
    !Config Units = [degC]
    CALL getin_p('TMIN_MAINT_RESP',tmin_maint_resp)

    !Config Key   = E0_MAINT_RESP
    !Config Desc  = maintenance respiration Temperature coefficient
    !Config if    = OK_STOMATE
    !Config Def   = undef, 308.56, 308.56, 308.56, 308.56, 308.56, 308.56, 308.56, 308.56, 308.56, 308.56, 308.56, 308.56   
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('E0_MAINT_RESP',e0_maint_resp)

    !
    ! Allocation
    !
    !Config Key   = TREF_LABILE
    !Config Desc  = Growth from labile pool - temperature at which all labile Cmaintenance respiration Temperature coefficient
    !Config if    = OK_STOMATE
    !Config Def   = undef, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5    
    !Config Help  =
    !Config Units = [degC]
    CALL getin_p('TREF_LABILE',tref_labile)

    !Config Key   = TMIN_LABILE
    !Config Desc  = Growth from labile pool  - temperature above which labile will be allocated to growth
    !Config if    = OK_STOMATE
    !Config Def   = undef, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2
    !Config Help  =
    !Config Units = [degC]
    CALL getin_p('TMIN_LABILE',tmin_labile)

    !Config Key   = E0_LABILE
    !Config Desc  = Growth temperature coefficient - tuned see stomate_growth_fun_all.f90
    !Config if    = OK_STOMATE
    !Config Def   = undef, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('E0_LABILE',e0_labile)

    !Config Key   = ALWAYS_LABILE
    !Config Desc  = share of the labile pool that will remain in the labile pool
    !Config if    = OK_STOMATE
    !Config Def   = undef, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('ALWAYS_LABILE',always_labile)


    !
    ! Flux - LUC
    !
    !Config Key   = COEFF_LCCHANGE_s
    !Config Desc  = Coeff of biomass export for the year
    !Config if    = OK_STOMATE
    !Config Def   = undef, 0.897, 0.897, 0.597, 0.597, 0.597, 0.597, 0.597, 0.597, 0.597, 0.597, 0.597, 0.597 
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('COEFF_LCCHANGE_s',coeff_lcchange_s)
    
    !Config Key   = COEFF_LCCHANGE_m
    !Config Desc  = Coeff of biomass export for the decade
    !Config if    = OK_STOMATE
    !Config Def   = undef, 0.103, 0.103, 0.299, 0.299, 0.299, 0.299, 0.299, 0.299, 0.299, 0.403, 0.299, 0.403
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('COEFF_LCCHANGE_m',coeff_lcchange_m)
    
    !Config Key   = COEFF_LCCHANGE_l
    !Config Desc  = Coeff of biomass export for the century
    !Config if    = OK_STOMATE
    !Config Def   = undef, 0., 0., 0.104, 0.104, 0.104, 0.104, 0.104, 0.104, 0.104, 0., 0.104, 0.
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('COEFF_LCCHANGE_l',coeff_lcchange_l)

    !
    ! Phenology
    !

    !Config Key   = LAI_MAX_TO_HAPPY
    !Config Desc  = threshold of LAI below which plant uses carbohydrate reserves
    !Config if    = OK_STOMATE
    !Config Def   = undef, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5 
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('LAI_MAX_TO_HAPPY',lai_max_to_happy) 

    !Config Key   = LAI_MAX
    !Config Desc  = maximum LAI, PFT-specific
    !Config if    = OK_STOMATE
    !Config Def   = undef, 7.0, 5.0, 5.0, 4.0, 5.0, 3.5, 4.0, 3.0, 2.5, 2.0, 5.0, 5.0
    !Config Help  =
    !Config Units = [m^2/m^2]
    CALL getin_p('LAI_MAX',lai_max)

    !Config Key   = PHENO_TYPE
    !Config Desc  = type of phenology, 0=bare ground 1=evergreen,  2=summergreen,  3=raingreen,  4=perennial
    !Config if    = OK_STOMATE
    !Config Def   = 0, 1, 3, 1, 1, 2, 1, 2, 2, 4, 4, 2, 3
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('PHENO_TYPE',pheno_type)

    !
    ! Phenology : Leaf Onset
    !
    !Config Key   = FORCE_PHENO
    !Config Desc  = Offset from mean doy at which phenology will be forced
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, 42, undef, undef, 42, undef, 28, 28, 35, 35, 28, 28
    !Config Help  =
    !Config Units = [days]
    CALL getin_p('FORCE_PHENO',force_pheno)

    !Config Key   = PHENO_GDD_CRIT_C
    !Config Desc  = critical gdd, tabulated (C), constant c of aT^2+bT+c
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef, 640., 400., 400., 450.
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('PHENO_GDD_CRIT_C',pheno_gdd_crit_c)

    !Config Key   = PHENO_GDD_CRIT_B
    !Config Desc  = critical gdd, tabulated (C), constant b of aT^2+bT+c
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef, 32., 0., 6.25, 0.
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('PHENO_GDD_CRIT_B',pheno_gdd_crit_b)

    !Config Key   = PHENO_GDD_CRIT_A
    !Config Desc  = critical gdd, tabulated (C), constant a of aT^2+bT+c
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef, 0.03125,  0., 0., 0.
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('PHENO_GDD_CRIT_A',pheno_gdd_crit_a)

    !Config Key   = PHENO_MOIGDD_T_CRIT
    !Config Desc  = Average temperature threashold for C4 grass used in pheno_moigdd
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, 22.0, undef, undef
    !Config Help  =
    !Config Units = [C]
    CALL getin_p('PHENO_MOIGDD_T_CRIT',pheno_moigdd_t_crit)

    !Config Key   = NGD_CRIT
    !Config Desc  = critical ngd, tabulated. Threshold -5 degrees
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, 23.9, undef, undef, undef, undef, undef
    !Config Help  = NGD : Number of Growing Days.
    !Config Units = [days]
    CALL getin_p('NGD_CRIT',ngd_crit)

    !Config Key   = NCDGDD_TEMP
    !Config Desc  = critical temperature for the ncd vs. gdd function in phenology
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, undef, undef, undef, 5., undef, 0., undef, undef, undef, undef, undef
    !Config Help  =
    !Config Units = [C] 
    CALL getin_p('NCDGDD_TEMP',ncdgdd_temp)

    !Config Key   = HUM_FRAC
    !Config Desc  = critical humidity (relative to min/max) for phenology
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, .5, undef, undef, undef, undef, undef,  undef, .5, .5, .5,.5     
    !Config Help  =
    !Config Units = [%]
    CALL getin_p('HUM_FRAC',hum_frac)

    !Config Key   = HUM_MIN_TIME
    !Config Desc  = minimum time elapsed since moisture minimum
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, 50., undef, undef, undef, undef, undef, undef, 35., 35., 75., 75.
    !Config Help  =
    !Config Units = [days]
    CALL getin_p('HUM_MIN_TIME',hum_min_time)

    !Config Key   = RELSOILMOIST_ALWAYS
    !Config Desc  = relative soil moisture availability above which moisture tendency doesn't matter 
    !Config If    = OK_STOMATE 
    !Config Def   = undef, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.98, 0.6, 0.6, 0.6 
    !Config Help  = 
    !Config Units = [-]   
    CALL getin_p('RELSOILMOIST_ALWAYS',relsoilmoist_always)
    !
    !Config Key   = T_ALWAYS
    !Config Desc  = monthly temp. above which temp. tendency doesn't matter 
    !Config If    = OK_STOMATE 
    !Config Def   = undef, 283.15, 283.15, 283.15, 283.15, 283.15, 283.15, 283.15, 283.15, 280., 283.15, 283.15, 283.15 
    !Config Help  = 
    !Config Units = [K]    
    CALL getin_p('T_ALWAYS',t_always)
    !
    !Config Key   = HARVEST_TIME_A
    !Config Desc  = Parameter to calculate harvest time of crops
    !Config If    = OK_STOMATE 
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, 160, 160
    !Config Help  = 
    !Config Units = [day]    
    CALL getin_p('HARVEST_TIME_A',harvest_time_a)
    !
    !Config Key   = HARVEST_TIME_B
    !Config Desc  = Parameter to calculate harvest time of crops
    !Config If    = OK_STOMATE 
    !Config Def   = undef, 
    !Config Help  = undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, 282, 282
    !Config Units = [K]    
    CALL getin_p('HARVEST_TIME_B',harvest_time_b)
    !
    !Config Key   = LEAF_AGE_CRIT_TREF
    !Config Desc  = Reference temperature
    !Config if    = OK_STOMATE
    !Config Def   = undef, 25., 25., 15., 20., 15., 5., 5., 5., 15., 20., 15., 20.
    !Config Help  = Reference temperature of the PFT (degrees Celsius)
    !               Used to calculate the leaf_age_crit as a function of longevity_leaf
    !Config Units = [degrees C]
    CALL getin_p('LEAF_AGE_CRIT_TREF',leaf_age_crit_tref)

    !Config Key   = LEAF_AGE_CRIT_COEFF1
    !Config Desc  = Coeff1 (unitless) to link leaf_age_crit to leaf_age_crit_tref
    !Config if    = OK_STOMATE
    !Config Def   = undef, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5
    !Config Help  = 
    !Config Units = [-]
    CALL getin_p('LEAF_AGE_CRIT_COEFF1',leaf_age_crit_coeff1)

    !Config Key   = LEAF_AGE_CRIT_COEFF2
    !Config Desc  = Coeff1 (unitless) to link leaf_age_crit to leaf_age_crit_tref
    !Config if    = OK_STOMATE
    !Config Def   = undef, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75
    !Config Help  = 
    !Config Units = [-]
    CALL getin_p('LEAF_AGE_CRIT_COEFF2',leaf_age_crit_coeff2)

    !Config Key   = LEAF_AGE_CRIT_COEFF3
    !Config Desc  = Coeff1 (unitless) to link leaf_age_crit to leaf_age_crit_tref
    !Config if    = OK_STOMATE
    !Config Def   = undef, 10., 10., 10., 10., 10., 10., 10., 10., 10., 10., 10., 10.
    !Config Help  = 
    !Config Units = [-]
    CALL getin_p('LEAF_AGE_CRIT_COEFF3',leaf_age_crit_coeff3)
    

    !Config Key   = LONGEVITY_FRUIT
    !Config Desc  = fruit lifetime
    !Config if    = OK_STOMATE
    !Config Def   = undef, 90., 90., 90., 90., 90., 90., 90., 90., undef, undef, undef, undef
    !Config Help  =
    !Config Units = [days]
    CALL getin_p('LONGEVITY_FRUIT',longevity_fruit)

    !Config Key   = LONGEVITY_LEAF
    !Config Desc  = leaf longevity
    !Config if    = OK_STOMATE
    !Config Def   = undef, 730., 180., 910., 730., 180., 910., 180., 180., 120., 120., 90., 90. 
    !Config Help  =
    !Config Units = [days]
    CALL getin_p('LONGEVITY_LEAF',longevity_leaf)
    
    !Config Key   = ECUREUIL
    !Config Desc  = fraction of primary leaf and root allocation put into reserve
    !Config if    = OK_STOMATE
    !Config Def   = undef, .0, 1., .0, .0, 1., .0, 1., 1., 1., 1., 1., 1.
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('ECUREUIL',ecureuil)

    !Config Key   = ALLOC_MIN
    !Config Desc  = minimum allocation above/below = f(age) - 30/01/04 NV/JO/PF
    !Config if    = OK_STOMATE
    !Config Def   = undef, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, undef, undef, undef, undef 
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('ALLOC_MIN',alloc_min)

    !Config Key   = ALLOC_MAX
    !Config Desc  = maximum allocation above/below = f(age) - 30/01/04 NV/JO/PF
    !Config if    = OK_STOMATE
    !Config Def   = undef, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, undef, undef, undef, undef
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('ALLOC_MAX',alloc_max)

    !Config Key   = DEMI_ALLOC 
    !Config Desc  = mean allocation above/below = f(age) - 30/01/04 NV/JO/PF
    !Config if    = OK_STOMATE
    !Config Def   = undef, 5., 5., 5., 5., 5., 5., 5., 5., undef, undef, undef, undef
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('DEMI_ALLOC',demi_alloc)

    !
    ! SOM decomposition (stomate)
    !

    !Config Key   = LC_leaf 
    !Config Desc  = Lignine/C ratio of leaf pool
    !Config If    = OK_STOMATE 
    !Config Def   = undef, 0.18, 0.18, 0.24, 0.18, 0.18, 0.24, 0.18, 0.24, 0.09, 0.09, 0.09, 0.09
    !Config Help  = 
    !Config Units = [-]   
    CALL getin_p('LC_leaf',LC_leaf)

    !Config Key   = LC_sapabove 
    !Config Desc  = Lignine/C ratio of sapabove pool
    !Config If    = OK_STOMATE 
    !Config Def   = undef, 0.23, 0.23, 0.29, 0.23, 0.23, 0.29, 0.23, 0.29, 0.09, 0.09, 0.09, 0.09
    !Config Help  = 
    !Config Units = [-]   
    CALL getin_p('LC_sapabove',LC_sapabove)

    !Config Key   = LC_sapbelow 
    !Config Desc  = Lignine/C ratio of sapbelow pool
    !Config If    = OK_STOMATE 
    !Config Def   = undef, 0.23, 0.23, 0.29, 0.23, 0.23, 0.29, 0.23, 0.29, 0.09, 0.09, 0.09, 0.09
    !Config Help  = 
    !Config Units = [-]   
    CALL getin_p('LC_sapbelow',LC_sapbelow)

    !Config Key   = LC_heartabove 
    !Config Desc  = Lignine/C ratio of heartabove pool
    !Config If    = OK_STOMATE 
    !Config Def   = undef, 0.23, 0.23, 0.29, 0.23, 0.23, 0.29, 0.23, 0.29, 0.09, 0.09, 0.09, 0.09
    !Config Help  = 
    !Config Units = [-]   
    CALL getin_p('LC_heartabove',LC_heartabove)

    !Config Key   = LC_heartbelow 
    !Config Desc  = Lignine/C ratio of heartbelow pool
    !Config If    = OK_STOMATE 
    !Config Def   = undef, 0.23, 0.23, 0.29, 0.23, 0.23, 0.29, 0.23, 0.29, 0.09, 0.09, 0.09, 0.09
    !Config Help  = 
    !Config Units = [-]   
    CALL getin_p('LC_heartbelow',LC_heartbelow)

    !Config Key   = LC_fruit
    !Config Desc  = Lignine/C ratio of fruit pool
    !Config If    = OK_STOMATE 
    !Config Def   = undef, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09
    !Config Help  = 
    !Config Units = [-]   
    CALL getin_p('LC_fruit',LC_fruit)

    !Config Key   = LC_root
    !Config Desc  = Lignine/C ratio of fruit pool
    !Config If    = OK_STOMATE 
    !Config Def   = undef, 0.22, 0.22, 0.22, 0.22, 0.22, 0.22, 0.22, 0.22, 0.22, 0.22, 0.22, 0.22
    !Config Help  = 
    !Config Units = [-]   
    CALL getin_p('LC_root',LC_root)

    !Config Key   = LC_carbres
    !Config Desc  = Lignine/C ratio of carbres pool
    !Config If    = OK_STOMATE 
    !Config Def   = undef, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.
    !Config Help  = 
    !Config Units = [-]   
    CALL getin_p('LC_carbres',LC_carbres)

    !Config Key   = LC_labile
    !Config Desc  = Lignine/C ratio of labile pool
    !Config If    = OK_STOMATE 
    !Config Def   = undef, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.
    !Config Help  = 
    !Config Units = [-]   
    CALL getin_p('LC_labile',LC_labile)

    !Config Key   = DECOMP_FACTOR 
    !Config Desc  = Multpliactive factor modifying the standard decomposition factor for each SOM pool 
    !Config if    = 
    !Config Def   = undef, 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.2, 1.4 
    !Config Help  = 
    !Config Units = 
    CALL getin_p("DECOMP_FACTOR",decomp_factor)

    !
    ! Stand structure - stomate
    !
    !Config Key   = MASS_RATIO_HEART_SAP 
    !Config Desc  = mass ratio (heartwood+sapwood)/heartwood 
    !Config If    = OK_STOMATE 
    !Config Def   = undef, 3., 3., 3., 3., 3., 3., 3., 3., 0., 0., 0., 0.  
    !Config Help  = 
    !Config Units = [-]   
    CALL getin_p('MASS_RATIO_HEART_SAP',mass_ratio_heart_sap)

    !Config Key   = CANOPY_COVER
    !Config Desc  = Test values for canopy cover
    !Config if    = OK_STOMATE
    !Config Def   = undef, 0.9, 0.9, 0.7, 0.7, 0.7, 0.6, 0.5, 0.5, 0.9, 0.9, 0.9, 0.9
    !Config Help  =
    !Config Units = [-] 
    CALL getin_p('CANOPY_COVER',canopy_cover)

    !Config Key   = NMAXPLANTS
    !Config Desc  = number of grasses and crops planted at the start of a rotation
    !Config if    = OK_STOMATE 
    !Config Def   = (undef, 15, 15, 15, 15, 15, 15, 15, 15, 10., 10., 10., 10.)*1.e3
    !Cofig Help   = An individual grass or crop is 1m2. Number of grasses and crops 
    !               planted at the start of a rotation
    !Config Units = [trees ha-1]
    CALL getin_p("NMAXPLANTS",nmaxplants)

    !Config Key   = HEIGHT_INIT
    !Config Desc  = height of a newly established vegetation
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef, 0.3, 0.3, 0.3, 0.3
    !Config Help  = Initially defined for trees, grasses and crops. Now
    !               only used for grasses and crops.
    !Config Units = [m]
    CALL getin_p("HEIGHT_INIT",height_init)

    !Config Key   = BIOMASS_INIT
    !Config Desc  = biomass of newly established vegetation
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef, 20, 30, 70, 80
    !Config Help  = Used for crops.
    !Config Units = [g]
    CALL getin_p("BIOMASS_INIT",biomass_init)

    !Config Key   = QMD_INIT
    !Config Desc  = minimum diameter of a newly established forest stand
    !Config if    = OK_STOMATE
    !Config Def   = undef, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, undef, undef, undef, undef
    !Config Help  = 
    !Config Units = [m]
    CALL getin_p("QMD_INIT",qmd_init)
    
    !Config Key   = BETA_SELF_THINNING
    !Config Desc  = beta coefficient of the self thinning relationship
    !Config if    = OK_STOMATE 
    !Config Def   = undef, -0.57, -0.57, -0.55, -0.61, -0.58, -0.55, -0.56, -0.56, undef, undef, undef, undef
    !Config Help  = 
    !Config Units = [-]
    CALL getin_p("BETA_SELF_THINNING",beta_self_thinning)

    !Config Key   = REF_ALPHA_SELF_THIN
    !Config Desc  = Basic change rate of the alpha self thinning parameter over time
    !Config If    = TEST_DYNAMIC_ALPHA_SELF_THIN, OK_STOMATE
    !Config Def   = undef, 2827, 3700, 1348, 1220, 2000, 2827, 800, 800, undef, undef, undef, undef
    !Config Help  = Basic change rate of the alpha self thinning parameter over time
    !Config Units = [-]
    CALL getin_p('REF_ALPHA_SELF_THIN',ref_alpha_self_thin)

    !Config Key   = INIT_PRE_INDUST_REF_GPP
    !Config Desc  = Reference pre-industrial gpp.
    !Config If    = TEST_DYNAMIC_ALPHA_SELF_THIN, OK_STOMATE
    !Config Def   = undef, 8.21, 8.21, 4.11, 3.40, 4.11, 2.74, 2.47, 2.47, undef, undef, undef, undef
    !Config Help  = The pre-industrial reference gpp that is used in the calculations of the temporal and spatial dynamics of alpha_self_thinning
    !Config Units = [gC m-2 day-1]
    CALL getin_p('INIT_PRE_INDUST_REF_GPP',init_pre_indust_ref_gpp)
    
    !Config Key   = FUELWOOD_DIAMETER
    !Config Desc  = Diameter below which harvest will be used as fuelwood
    !Config if    = OK_STOMATE, DIMENSIONAL WOOD PRODUCTS
    !Config Def   = undef, 0.3, 0.3, 0.2, 0.3, 0.3, 0.2, 0.2, 0.2, undef, undef, undef, undef
    !Config Help  = 
    !Config Units = [m]
    CALL getin_p("FUELWOOD_DIAMETER",fuelwood_diameter)
    
    !Config Key   = COPPICE_KILL_BE_WOOD
    !Config Desc  = The fraction of belowground wood killed during coppicing
    !Config if    = FOREST_MANAGED equals to 3 (Coppice)
    !Config Def   = undef, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, undef, undef, undef, undef
    !Config Help  = 
    !Config Units = [m]
    CALL getin_p("COPPICE_KBEW",coppice_kill_be_wood)

    !Config Key   = DELEUZE_A
    !Config Desc  = intercept of the intra-tree competition within a stand
    !Config if    = OK_STOMATE, NCIRC>6
    !Config Def   = undef, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23, undef, undef, undef, undef
    !Config Help  = intercept of the intra-tree competition within a stand
    !               based on the competion rule of Deleuze and Dhote 2004
    !               Used when n_circ > 6
    !Config Units = [-]
    CALL getin_p("DELEUZE_A",deleuze_a)

    !Config Key   = DELEUZE_B
    !Config Desc  = slope of the intra-tree competition within a stand
    !Config if    = OK_STOMATE, NCIRC>6
    !Config Def   = undef, 0.58, 0.58, 0.58, 0.58, 0.58, 0.58, 0.58, 0.58, 0.58, undef, undef, undef, undef
    !Config Help  = slope of the intra-tree competition within a stand
    !               based on the competion rule of Deleuze and Dhote 2004
    !               Used when n_circ > 6
    !Config Units = [-]
    CALL getin_p("DELEUZE_B",deleuze_b)

    !Config Key   = DELEUZE_P_ALL
    !Config Desc  = Percentile of the circumferences that receives photosynthates
    !Config if    = OK_STOMATE, NCIRC>1 AND NCIRC<6
    !Config Def   = undef, 0.5, 0.5, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, undef, undef, undef, undef
    !Config Help  = Percentile of the circumferences that receives photosynthates
    !               based on the competion rule of Deleuze and Dhote 2004
    !               Used when n_circ < 6 for FM 1, FM2 and FM4
    !Config Units = [0-1]
    CALL getin_p("DELEUZE_P_ALL",deleuze_p_all)
 
    !Config Key   = DELEUZE_P_COPPICE
    !Config Desc  = Percentile of the circumferences that receives photosynthates
    !Config if    = OK_STOMATE, functional allocation  
    !Config Def   = undef, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, undef, undef, undef, undef
    !Config Help  = Percentile of the circumferences that receives photosynthates
    !               based on the competion rule of Deleuze and Dhote 2004
    !               Used when n_circ < 6 for FM3
    !Config Units = [0-1]
    CALL getin_p("DELEUZE_P_COPPICE",deleuze_p_coppice)
    
    !Config Key   = DELEUZE_POWER_A
    !Config Desc  = Slope parameter for intra-specific competition
    !Config if    = OK_STOMATE
    !Config Def   = undef, 0, 0, 0, 0, 0, 0, 0, 0, 0, undef, undef, undef, undef
    !Config Help  = Divisor of the power for the slope of the intra-tree competition within a stand
    !               based on the competion rule of Deleuze and Dhote 2004.
    !Config Units = [-]
    CALL getin_p("DELEUZE_POWER_A",deleuze_power_a)
    
    !Config Key   = M_DV
    !Config Desc  = Relaxation factor of deleuze relationship
    !Config if    = OK_STOMATE, NCIRC>1 
    !Config Def   = undef, 1.05, 1.05, 1.05, 1.05, 1.05, 1.05, 1.05, 1.05, 1.05, undef, undef, undef, undef
    !Config Help  = Allows some allocation to trees below the threshold (sigma) if 
    !               the value exceeds 0.
    !Config Units = [-]
    CALL getin_p("M_DV",m_dv)
    
    !Config Key   = FOREST_MANAGED_FORCED
    !Config Desc  = Prescribe forest management
    !Config if    = OK_STOMATE
    !Config Def   = undef, 1, 1, 1, 1, 1, 1, 1, 1, undef, undef, undef, undef
    !Config Help  = If no forest management map is used, forest management has to
    !               be prescribed 
    !Config Units = [-]
    CALL getin_p("FOREST_MANAGED_FORCED",forest_managed_forced)

    !Config Key   = LARGEST_TREE_DIA
    !Config Desc  = Maximum tree diameter of a stand
    !Config if    = OK_STOMATE
    !Config Def   = 0.0, .45, .60, .45, .45, .60, .45, .60, .45, 0.0, 0.0, 0.0, 0.0
    !Config Help  = If the mean diameter of the largest diameter class
    !               exceeds this threshold, the stand will be killed and replanted.
    !Config Units = [m]
    CALL getin_p("LARGEST_TREE_DIA",largest_tree_dia)
    !Config Key   = ROTATION_REF
    !Config Desc  = Rotation de reference par PFT (an), gestion multifonctionnelle
    !Config If    = OK_MANAGEMENT_INTENSITY
    !Config Def   = rotation_ref_mtc(pft_to_mtc)
    !Config Help  = Le taux de recolte vaut 1/(ROTATION_REF * MI_ROTATION_FACTOR(classe)).
    !Config         Un taux estime sur EFDA a ete essaye puis retire : 39 ans de serie ne
    !Config         peuvent pas mesurer une rotation de 100-200 ans. Valeurs issues des
    !Config         itineraires techniques CNPF et du PROF Lisboa e Vale do Tejo (ICNF).
    !Config         /!\ L'EUCALYPTUS n'a pas de MTC propre : poser ici sa rotation sur
    !Config         ses creneaux PFT (10-14 ans en taillis, 25-35 ans en futaie).
    !Config         Valeur <= 0 => pas de rotation definie pour ce PFT.
    !Config Units = [an]
    CALL getin_p("ROTATION_REF",rotation_ref)
    !Config Key   = FM_SRC_PFT
    !Config Desc  = Force le regime taillis a courte rotation (ifm_src) sur ces PFT
    !Config If    = OK_READ_FM_MAP
    !Config Def   = n (pour tous les PFT)
    !Config Help  = ifm_src est entierement implemente (parametres *_SRC, coupes
    !Config         icut_cop1/2/3, first/second_coppice_cut) mais n'etait jamais affecte.
    !Config         Motif : l'eucalyptus iberique se conduit en taillis par rejet de
    !Config         souche, 1 plantation puis 2-3 revolutions de 10-14 ans (PROF LVT,
    !Config         ICNF 2018) -- ni plantation ni regeneration naturelle au sens du modele.
    !Config Units = [FLAG]
    CALL getin_p("FM_SRC_PFT",fm_src_pft)

    !Config Key   = FM_COP_PFT
    !Config Desc  = Force le taillis simple (ifm_cop) sur ces PFT, la ou la carte
    !Config         de gestion les dit geres.
    !Config If    = OK_READ_FM_MAP
    !Config Def   = n (pour tous les PFT)
    !Config Help  = Meme patron que FM_SRC_PFT, autre regime : ifm_cop est le taillis
    !Config         a rotation LONGUE. Motif : le chene vert mediterraneen se conduit
    !Config         en taillis simple de 35-60 ans (CNPF, itineraire I5_CV), pas en
    !Config         futaie reguliere. Exclusif de FM_SRC_PFT sur un meme PFT.
    !Config Units = [FLAG]
    CALL getin_p("FM_COP_PFT",fm_cop_pft)

    ! Guillaume M. -- The two coppice flags name mutually exclusive regimes on one PFT:
    ! long rotation (ifm_cop) versus short rotation (ifm_src). Setting both leaves the
    ! outcome to the order of two WHERE blocks, which is not a user decision. Refuse it.
    IF (ANY(fm_src_pft(:) .AND. fm_cop_pft(:))) THEN
       WRITE(numout,*) 'FM_SRC_PFT et FM_COP_PFT sur un meme PFT : ', &
            PACK((/(ivm, ivm=1,nvm)/), fm_src_pft(:) .AND. fm_cop_pft(:))
       CALL ipslerr_p(3,'pft_parameters', &
            'FM_SRC_PFT et FM_COP_PFT sont exclusifs sur un meme PFT', &
            'ifm_src est le taillis a courte rotation, ifm_cop le taillis simple', &
            'choisir l un des deux')
    ENDIF

    ! Guillaume M. -- The six SPITFIRE fire parameters below are read per PFT from the
    ! namelist; defaults stay the *_mtc values assigned above, so a run without a namelist
    ! entry is bit-neutral. Two PFTs sharing one MTC can now differ in fire behaviour.
    ! /!\ Do not set a thin-bark species by expert guess: post-fire resprouting is not
    ! represented, so a realistic thin bark over-kills it. Calibrate on observed burnt area.

    !Config Key   = F_SH
    !Config Desc  = Scorch height parameter for crown fire (SPITFIRE)
    !Config if    = OK_SPITFIRE
    !Config Def   = f_sh_mtc(pft_to_mtc)
    !Config Help  = Higher value = taller scorch height for a given fire intensity,
    !               hence more crown damage. MTC5 already carries the highest value
    !               of the 13 MTCs (0.371).
    !Config Units = [-]
    CALL getin_p("F_SH",f_sh)

    !Config Key   = BTPAR1
    !Config Desc  = Bark thickness allometry parameter 1 (SPITFIRE)
    !Config if    = OK_SPITFIRE
    !Config Def   = BTpar1_mtc(pft_to_mtc)
    !Config Help  = Bark thickness drives cambial insulation; time to kill the
    !               cambium scales with bark thickness squared. Eucalypts have a
    !               NEGATIVE bark allometry (thin bark) - see the warning above.
    !Config Units = [-]
    CALL getin_p("BTPAR1",BTpar1)

    !Config Key   = BTPAR2
    !Config Desc  = Bark thickness allometry parameter 2 (SPITFIRE)
    !Config if    = OK_SPITFIRE
    !Config Def   = BTpar2_mtc(pft_to_mtc)
    !Config Units = [-]
    CALL getin_p("BTPAR2",BTpar2)

    !Config Key   = R_CK
    !Config Desc  = Post-fire mortality parameter from crown damage (SPITFIRE)
    !Config if    = OK_SPITFIRE
    !Config Def   = r_ck_mtc(pft_to_mtc)
    !Config Units = [-]
    CALL getin_p("R_CK",r_ck)

    !Config Key   = P_CK
    !Config Desc  = Post-fire mortality parameter from crown damage (SPITFIRE)
    !Config if    = OK_SPITFIRE
    !Config Def   = p_ck_mtc(pft_to_mtc)
    !Config Units = [-]
    CALL getin_p("P_CK",p_ck)

    !Config Key   = ME
    !Config Desc  = Moisture of extinction / flammability threshold (SPITFIRE)
    !Config if    = OK_SPITFIRE
    !Config Def   = me_mtc(pft_to_mtc)
    !Config Help  = Lower value = the fuel stops burning at a lower moisture,
    !               i.e. a LESS flammable PFT. Raising it makes the PFT burn under
    !               wetter conditions.
    !Config Units = [-]
    CALL getin_p("ME",me)

    !Config Key   = THINSTRAT
    !Config Desc  = Thinnig management style (from below or from above)
    !Config if    = FOREST_MANAGEMENT 
    !Config Def   = 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0
    !Config Help  = 
    !Config Units = [-]
    CALL getin_p("THINSTRAT",thinstrat)

    !Config Key   = TAUMIN
    !Config Desc  = Minimum probability that a tree get thinned
    !Config if    = FOREST_MANAGEMENT 
    !Config Def   = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    !Config Help  = 
    !Config Units = [-]
    CALL getin_p("TAUMIN",taumin)

    !Config Key   = TAUMAX
    !Config Desc  = Maximum probability that a tree get thinned
    !Config if    = FOREST_MANAGEMENT 
    !Config Def   = 0.0, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0
    !Config Help  = 
    !Config Units = [-]
    CALL getin_p("TAUMAX",taumax)

    !Config Key   = RDI_DIA_REF
    !Config Desc  = RDI at a reference tree diameter used to shape the end of
    !the RDI/diameter relationship
    !Config if    = OK_STOMATE 
    !Config Def   = undef, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, undef, undef, undef, undef
    !Config Help  = 
    !Config Units = [-]
    CALL getin_p("RDI_DIA_REF", rdi_dia_ref)

    !Config Key   = THIN_INTENSITY
    !Config Desc  = Intensity of a thinnig event when a forest reach the upper
    !rdi limit
    !Config if    = OK_STOMATE 
    !Config Def   = undef, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, undef, undef, undef, undef
    !Config Help  = 
    !Config Units = [-]
    !Config Key   = RDI_BAND_LOW / RDI_BAND_HIGH
    !Config Desc  = Bande de RDI cible des classes intermediaires, PAR CRENEAU
    !Config If    = OK_CLEARCUT_LAST_CLASS
    !Config Def   = 0.70 / 0.80
    !Config Help  = Court-circuite les polynomes calculate_rdi_boundaries pour les classes
    !Config         2..nagec-1 : au-dela de la classe 1 la cible n'est plus une propriete
    !Config         du peuplement mais une CONSIGNE de gestion.
    !Config         /!\ INDEXES PAR CRENEAU depuis le 2026-08-12. Ils etaient SCALAIRES, donc
    !Config         incapables de distinguer la classe 2 de la classe 3, et a fortiori de
    !Config         varier a l'interieur d'une classe. Or la conduite en arbres d'avenir
    !Config         demande une consigne qui DESCEND pendant la classe 3 : la classe 2 doit
    !Config         rester dense (education) pendant que la classe 3 s'ouvre (amelioration).
    !Config         Une valeur unique au namelist remplit tous les creneaux, donc les
    !Config         configurations existantes sont inchangees.
    !Config         Cela rend aussi possible l'extension annoncee ici de longue date :
    !Config         indexer la bande sur l'intensite de gestion.
    !Config Units = [-]
    rdi_band_low(:)  = 0.70_r_std
    rdi_band_high(:) = 0.80_r_std
    CALL getin_p('RDI_BAND_LOW',  rdi_band_low)
    CALL getin_p('RDI_BAND_HIGH', rdi_band_high)
    ! Guillaume M. -- Guard applied slot by slot: an inverted band on a single slot is
    ! enough to make thinning inconsistent there, and a global guard would not see it.
    DO ivm = 1, nvm
       IF (rdi_band_low(ivm) <= zero .OR. rdi_band_high(ivm) > un .OR. &
            rdi_band_low(ivm) >= rdi_band_high(ivm)) THEN
          WRITE(numout,*) 'creneau, RDI_BAND_LOW, RDI_BAND_HIGH : ', &
               ivm, rdi_band_low(ivm), rdi_band_high(ivm)
          CALL ipslerr_p(3,'pft_parameters',&
               'il faut 0 < RDI_BAND_LOW < RDI_BAND_HIGH <= 1 sur CHAQUE creneau',&
               'la bande encadre le RDI sur pied observe (~0,686 Pucher)','')
       ENDIF
    ENDDO

    !Config Key   = RDI_RAMP_OUT
    !Config Desc  = RDI vise en SORTIE de la classe d'amelioration (rampe intra-classe)
    !Config If    = OK_RDI_RAMP
    !Config Def   = 0.04
    !Config Help  = La consigne de RDI DESCEND pendant la classe d'amelioration, de
    !Config         RDI_BAND_LOW a l'entree jusqu'a RDI_RAMP_OUT a la sortie. Motif : la
    !Config         densite y est en TRANSITION -- une consigne constante posee a l'entree
    !Config         serait un echelon (mesure : 91 % du peuplement retire en une annee,
    !Config         soit une coupe rase deguisee), alors que la meme reduction etalee sur
    !Config         les ~18 ans de la classe vaut 12,5 %/an, un regime d'eclaircie ordinaire.
    !Config         0.04 correspond a ~59 tiges/ha a 16 cm, soit les 50-70 arbres d'avenir
    !Config         des itineraires publies, et donne RDI ~0,30 en classe terminale au
    !Config         diametre d'exploitabilite. Voir design/MODULE_DESIGN_CTN_BASAL_AREA.md §7.
    !Config Units = [-]
    rdi_ramp_out(:) = 0.04_r_std
    CALL getin_p('RDI_RAMP_OUT', rdi_ramp_out)

    !Config Key   = FORCED_THIN_INTENSITY
    !Config Desc  = Fraction de RDI retiree chaque annee par l'eclaircie forcee
    !Config If    = OK_FORCED_THINNING
    !Config Def   = 0.02
    !Config Help  = La cible passee a `thinning` vaut rdi_courant * (1 - cette fraction).
    !Config         0.02 est CALCULE : la classe terminale porte 1480 tiges/ha pour une
    !Config         cible d'ordre 700, et le peuplement passe ~38 ans dans les classes 1-2,
    !Config         d'ou (1-r)^38 = 0.47 soit r = 1.97 %/an. Ignore la retroaction --
    !Config         eclaircir fait grossir les survivants, l'effet sur le QMD sera plus fort
    !Config         que proportionnel. Voir design/MODULE_DESIGN_FORCED_THINNING.md
    !Config Units = [-]
    forced_thin_intensity(:) = 0.02_r_std
    CALL getin_p('FORCED_THIN_INTENSITY', forced_thin_intensity)
    IF (ok_forced_thinning) THEN
       DO ivm = 1, nvm
          IF (forced_thin_intensity(ivm) < zero .OR. forced_thin_intensity(ivm) >= un) THEN
             WRITE(numout,*) 'creneau, FORCED_THIN_INTENSITY : ', ivm, forced_thin_intensity(ivm)
             CALL ipslerr_p(3,'pft_parameters','FORCED_THIN_INTENSITY hors [0,1[', &
                  'une intensite negative ou totale n a pas de sens','')
          ENDIF
       ENDDO
    ENDIF
    IF (ok_rdi_ramp) THEN
       DO ivm = 1, nvm
          IF (rdi_ramp_out(ivm) <= zero .OR. rdi_ramp_out(ivm) >= rdi_band_low(ivm)) THEN
             WRITE(numout,*) 'creneau, RDI_RAMP_OUT, RDI_BAND_LOW : ', &
                  ivm, rdi_ramp_out(ivm), rdi_band_low(ivm)
             CALL ipslerr_p(3,'pft_parameters',&
                  'il faut 0 < RDI_RAMP_OUT < RDI_BAND_LOW sur CHAQUE creneau',&
                  'la rampe DESCEND : une sortie au-dessus de l entree la renverserait','')
          ENDIF
       ENDDO
    ENDIF

    CALL getin_p("THIN_INTENSITY", thin_intensity)

    !Config Key   = RDI_START
    !Config Desc  = The rdi value at the start of a forest rotation, use
    !0.4 to obtain a reallistic growth rate during the early stage of the forest ) 
    !Config if    = OK_STOMATE 
    !Config Def   = undef, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, undef, undef, undef, undef
    !Config Help  = 
    !Config Units = [-]
    CALL getin_p("RDI_START", rdi_start)

    !Config Key   = FINAL_DIA_NONE
    !Config Desc  = The maximum Diameter before the decline in RDI
    !Config if    = OK_STOMATE
    !Config Def   = undef, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, undef, undef, undef, undef
    !Config Help  =
    !Config Units = [-]
    CALL getin_p("FINAL_DIA_NONE", final_dia(:,ifm_none))

    !Config Key   = FINAL_DIA_THIN
    !Config Desc  = The maximum Diameter before the decline in RD
    !Config if    = OK_STOMATE
    !Config Def   = undef, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, undef, undef, undef, undef
    !Config Help  =
    !Config Units = [-]
    CALL getin_p("FINAL_DIA_THIN", final_dia(:,ifm_thin))

    !Config Key   = FINAL_DIA_UNEVEN
    !Config Desc  = The maximum Diameter before the decline in RD
    !Config if    = OK_STOMATE
    !Config Def   = undef, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, undef, undef, undef, undef
    !Config Help  =
    !Config Units = [-]
    CALL getin_p("FINAL_DIA_UNEVEN", final_dia(:,ifm_uneven))

    !Config Key   = FINAL_DIA_COP
    !Config Desc  = The maximum Diameter before the decline in RD
    !Config if    = OK_STOMATE
    !Config Def   = undef, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, undef, undef, undef, undef
    !Config Help  =
    !Config Units = [-]
    CALL getin_p("FINAL_DIA_COP", final_dia(:,ifm_cop))

    !Config Key   = FINAL_DIA_SRC
    !Config Desc  = The maximum Diameter before the decline in RD
    !Config if    = OK_STOMATE
    !Config Def   = undef, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, undef, undef, undef, undef
    !Config Help  =
    !Config Units = [-]
    CALL getin_p("FINAL_DIA_SRC", final_dia(:,ifm_src))

    !Config Key   = RDI_MAX_NONE
    !Config Desc  = The maximum RDI that a forest will reach during the rotation
    !Config if    = OK_STOMATE 
    !Config Def   = undef, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, undef, undef, undef, undef
    !Config Help  = 
    !Config Units = [-]
    CALL getin_p("RDI_MAX_NONE", rdi_max(:,ifm_none))

    !Config Key   = RDI_MAX_THIN
    !Config Desc  = The maximum RDI that a forest will reach during the rotation
    !Config if    = OK_STOMATE 
    !Config Def   = undef, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, undef, undef,
    !undef, undef
    !Config Help  = 
    !Config Units = [-]
    CALL getin_p("RDI_MAX_THIN", rdi_max(:,ifm_thin))

    !Config Key   = RDI_MAX_UNEVEN
    !Config Desc  = The maximum RDI that a forest will reach during the rotation
    !Config if    = OK_STOMATE 
    !Config Def   = undef, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, undef, undef,
    !undef, undef
    !Config Help  = 
    !Config Units = [-]
    CALL getin_p("RDI_MAX_UNEVEN", rdi_max(:,ifm_uneven))

    !Config Key   = RDI_MAX_COP
    !Config Desc  = The maximum RDI that a forest will reach during the rotation
    !Config if    = OK_STOMATE 
    !Config Def   = undef, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, undef, undef,
    !undef, undef
    !Config Help  = 
    !Config Units = [-]
    CALL getin_p("RDI_MAX_COP", rdi_max(:,ifm_cop))

    !Config Key   = RDI_MAX_SRC
    !Config Desc  = The maximum RDI that a forest will reach during the rotation
    !Config if    = OK_STOMATE 
    !Config Def   = undef, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, undef, undef,
    !undef, undef
    !Config Help  = 
    !Config Units = [-]
    CALL getin_p("RDI_MAX_SRC", rdi_max(:,ifm_src))

    !Config Key   = DENS_TARGET_NONE
    !Config Desc  = Minimal density. Below this density the forest will be clearcut (trees.ha-1)
    !Config if    = OK_STOMATE
    !Config Def   = undef, 200, 200, 200, 200, 200, 200, 200, 200, undef, undef, undef, undef
    !Config Help  =
    !Config Units = [-]
    CALL getin_p("DENS_TARGET_NONE", dens_target(:,ifm_none))

    !Config Key   = DENS_TARGET_THIN
    !Config Desc  = Minimal density. Below this density the forest will be clearcut (trees.ha-1)
    !Config if    = OK_STOMATE
    !Config Def   = undef, 200, 200, 200, 200, 200, 200, 200, 200, undef, undef, undef, undef
    !Config Help  =
    !Config Units = [-]
    CALL getin_p("DENS_TARGET_THIN", dens_target(:,ifm_thin))

    !Config Key   = DENS_TARGET_UNEVEN
    !Config Desc  = Minimal density. Below this density the forest will be clearcut (trees.ha-1)
    !Config if    = OK_STOMATE
    !Config Def   = undef, 200, 200, 200, 200, 200, 200, 200, 200, undef, undef, undef, undef
    !Config Help  =
    !Config Units = [-]
    CALL getin_p("DENS_TARGET_UNEVEN", dens_target(:,ifm_uneven))

    !Config Key   = DENS_TARGET_COP
    !Config Desc  = Minimal density. Below this density the forest will be clearcut (trees.ha-1)
    !Config if    = OK_STOMATE
    !Config Def   = undef, 200, 200, 200, 200, 200, 200, 200, 200, undef, undef, undef, undef
    !Config Help  =
    !Config Units = [-]
    CALL getin_p("DENS_TARGET_COP", dens_target(:,ifm_cop))

    !Config Key   = DENS_TARGET_SRC
    !Config Desc  = Minimal density. Below this density the forest will be clearcut (trees.ha-1)
    !Config if    = OK_STOMATE
    !Config Def   = undef, 200, 200, 200, 200, 200, 200, 200, 200, undef, undef, undef, undef
    !Config Help  =
    !Config Units = [-]
    CALL getin_p("DENS_TARGET_SRC", dens_target(:,ifm_src))

    !Config Key   = DIA_BOOST_END
    !Config Desc  = Diameter at which the model stop maintaining a low rdi
    !(based on rdi_start)
    !Config if    = OK_STOMATE 
    !Config Def   = undef, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, undef, undef, undef, undef
    !Config Help  = 
    !Config Units = [-]
    CALL getin_p("DIA_BOOST_END", dia_boost_end)

    !Config Key   = CIRC_DIST_SHAPE_NONE 
    !Config Desc  = Shape parameter of the unmanaged forest circular distance distribution
    !Config if    = FOREST_MANAGEMENT
    !Config Def   = undef, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, undef, undef, undef, undef
    !Config Help  =
    !Config Units = [-]
    CALL getin_p("CIRC_DIST_SHAPE_NONE",circ_dist_shape_none)

    !Config Key   = CIRC_DIST_SHAPE_THIN
    !Config Desc  = Shape parameter of the managed forest circ class distribution
    !Config if    = FOREST_MANAGEMENT
    !Config Def   = undef, -0.09, -0.09, -0.09, -0.09, -0.09, -0.09, -0.09, -0.09, undef, undef, undef, undef
    !Config Help  =
    !Config Units = [-]
    CALL getin_p("CIRC_DIST_SHAPE_THIN",circ_dist_shape_thin) 

    !Config Key   = CIRC_DIST_SHAPE_COP
    !Config Desc  = Shape parameter of the coppice forest circ class distribution
    !Config if    = FOREST_MANAGEMENT
    !Config Def   = undef, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, undef, undef, undef, undef
    !Config Help  =
    !Config Units = [-]
    CALL getin_p("CIRC_DIST_SHAPE_COP",circ_dist_shape_cop)

    !Config Key   = CIRC_DIST_SHAPE_SRC
    !Config Desc  = Shape parameter of short rotation cycle for the circ class distribution
    !Config if    = FOREST_MANAGEMENT
    !Config Def   = undef, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, undef, undef, undef, undef
    !Config Help  =
    !Config Units = [-]
    CALL getin_p("CIRC_DIST_SHAPE_SRC",circ_dist_shape_src)

    !Config Key   = CIRC_DIST_SHAPE_UNEVEN
    !Config Desc  = Shape parameter of Uneven-age forest for the circ class distribution
    !Config if    = FOREST_MANAGEMENT
    !Config Def   = undef, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, undef, undef, undef, undef
    !Config Help  =
    !Config Units = [-]
    CALL getin_p("CIRC_DIST_SHAPE_UNEVEN",circ_dist_shape_uneven)    

    !Config Key   = BRANCH_HARVEST_NONE
    !Config Desc  = The fraction of branches which are harvested from unmanaged forests
    !Config if    = FOREST_MANAGEMENT 
    !Config Def   = 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0  
    !Config Help  = CHECK whether this parameter is already used in the code for this FM.
    !Config Units = [-]
    CALL getin_p("BRANCH_HARVEST_NONE",branch_harvest(:,ifm_none))

    !Config Key   = BRANCH_HARVEST_THIN
    !Config Desc  = The fraction of branches which are harvested from thin and fell forests
    !Config if    = FOREST_MANAGEMENT 
    !Config Def   = 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0  
    !Config Help  = FM2 has branch harvest implemented
    !Config Units = [-]
    CALL getin_p("BRANCH_HARVEST_THIN",branch_harvest(:,ifm_thin))

    !Config Key   = BRANCH_HARVEST_UNEVEN
    !Config Desc  = The fraction of branches which are harvested from uneven-aged managed forests
    !Config if    = FOREST_MANAGEMENT 
    !Config Def   = 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0  
    !Config Help  = CHECK whether this parameter is already used in the code for this FM. 
    !Config Units = [-]
    CALL getin_p("BRANCH_HARVEST_UNEVEN",branch_harvest(:,ifm_uneven))

    !Config Key   = BRANCH_HARVEST_COP
    !Config Desc  = The fraction of branches which are harvested from coppice forests
    !Config if    = FOREST_MANAGEMENT 
    !Config Def   = 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0  
    !Config Help  = CHECK whether this parameter is already used in the code for this FM.
    !Config Units = [-]
    CALL getin_p("BRANCH_HARVEST_COP",branch_harvest(:,ifm_cop))

    !Config Key   = BRANCH_HARVEST_SRC
    !Config Desc  = The fraction of branches which are harvested from unmanaged forests
    !Config if    = FOREST_MANAGEMENT 
    !Config Def   = 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0  
    !Config Help  = CHECK whether this parameter is alraedy used in the code for this FM.
    !Config Units = [-]
    CALL getin_p("BRANCH_HARVEST_SRC",branch_harvest(:,ifm_src))

    !Config Key   = CLEARCUT_RESIDUE_NONE
    !Config Desc  = The fraction of wood which are harvested during thinning and felling rest are left on site)
    !Config if    = FOREST_MANAGEMENT 
    !Config Def   = 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0  
    !Config Help  = CHECK whether this parameter is alraedy used in the code for this FM.
    !Config Units = [-]
    CALL getin_p("CLEARCUT_RESIDUE_NONE",clearcut_residue(:,ifm_none))

    !Config Key   = CLEARCUT_RESIDUE_THIN
    !Config Desc  = The fraction of wood which are harvested during thin and fell (the rest are left on site)
    !Config if    = FOREST_MANAGEMENT 
    !Config Def   = 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0  
    !Config Help  = FM2 has clearcut_residue implemented.
    !Config Units = [-]
    CALL getin_p("CLEARCUT_RESIDUE_THIN",clearcut_residue(:,ifm_thin))

    !Config Key   = CLEARCUT_RESIDUE_UNEVEN
    !Config Desc  = The fraction of wood which are harvested during uneven-aged forestry (the rest are left on site)
    !Config if    = FOREST_MANAGEMENT 
    !Config Def   = 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0  
    !Config Help  = FM5 has clearcut_residue implemented.
    !Config Units = [-]
    CALL getin_p("CLEARCUT_RESIDUE_UNEVEN",clearcut_residue(:,ifm_uneven))

    !Config Key   = CLEARCUT_RESIDUE_COP
    !Config Desc  = The fraction of wood which are harvested during coppicing (the rest are left on site)
    !Config if    = FOREST_MANAGEMENT 
    !Config Def   = 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0  
    !Config Help  = CHECK whether this parameter is alraedy used in the code for this FM.
    !Config Units = [-]
    CALL getin_p("CLEARCUT_RESIDUE_COP",clearcut_residue(:,ifm_cop))

    !Config Key   = CLEARCUT_RESIDUE_SRC
    !Config Desc  = The fraction of wood which are harvested during src (the rest are left on site)
    !Config if    = FOREST_MANAGEMENT 
    !Config Def   = 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0  
    !Config Help  = CHECK whether this parameter is alraedy used in the code for this FM.
    !Config Units = [-]
    CALL getin_p("CLEARCUT_RESIDUE_SRC",clearcut_residue(:,ifm_src))

    !Config Key   = SALVAGE_DIST_NONE
    !Config Desc  = The fraction of dead wood that is remove after disturbance (the rest are left on site)
    !Config if    = FOREST_MANAGEMENT 
    !Config Def   = 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0  
    !Config Help  = CHECK whether this parameter is alraedy used in the code for this FM.
    !Config Units = [-]
    CALL getin_p("SALVAGE_DIST_NONE",salvage_dist(:,ifm_none))

    !Config Key   = SALVAGE_DIST_THIN
    !Config Desc  = The fraction of dead wood that is remove after disturbance (the rest are left on site)
    !Config if    = FOREST_MANAGEMENT 
    !Config Def   = 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0  
    !Config Help  = CHECK whether this parameter is alraedy used in the code for this FM.
    !Config Units = [-]
    CALL getin_p("SALVAGE_DIST_THIN",salvage_dist(:,ifm_thin))

    !Config Key   = SALVAGE_DIST_UNEVEN
    !Config Desc  = The fraction of dead wood that is remove after disturbance (the rest are left on site)
    !Config if    = FOREST_MANAGEMENT 
    !Config Def   = 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0  
    !Config Help  = CHECK whether this parameter is alraedy used in the code for this FM.
    !Config Units = [-]
    CALL getin_p("SALVAGE_DIST_UNEVEN",salvage_dist(:,ifm_uneven))

    !Config Key   = SALVAGE_DIST_COP
    !Config Desc  = The fraction of dead wood that is remove after disturbance (the rest are left on site)
    !Config if    = FOREST_MANAGEMENT 
    !Config Def   = 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0  
    !Config Help  = CHECK whether this parameter is alraedy used in the code for this FM.
    !Config Units = [-]
    CALL getin_p("SALVAGE_DIST_COP",salvage_dist(:,ifm_cop))

    !Config Key   = SALVAGE_DIST_SRC
    !Config Desc  = The fraction of dead wood that is remove after disturbance (the rest are left on site)
    !Config if    = FOREST_MANAGEMENT 
    !Config Def   = 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0  
    !Config Help  = CHECK whether this parameter is alraedy used in the code for this FM.
    !Config Units = [-]
    CALL getin_p("SALVAGE_DIST_SRC",salvage_dist(:,ifm_src))

    !Config Key   = COPPICE_DIAMETER
    !Config Desc  = The trunk diameter at which a coppice will be cut
    !Config if    = FOREST_MANAGEMENT 
    !Config Def   = undef, 0.2, 0.2, 0.2, 0.2, 0.1, 0.2, 0.2, 0.2, undef, undef, undef, undef
    !Config Help  = 
    !Config Units = [m]
    CALL getin_p("COPPICE_DIAMETER",coppice_diameter)

    !Config Key   = SHOOTS_PER_STOOL
    !Config Desc  = The number of shoots that will regrow per stool after the first coppice cut
    !Config if    = FOREST_MANAGEMENT 
    !Config Def   = undef, 6, 6, 6, 6, 6, 6, 6, 6, undef, undef, undef, undef
    !Config Help  = 
    !Config Units = [shoots.stool-1]
    CALL getin_p("SHOOTS_PER_STOOL",shoots_per_stool)

    !Config Key   = SRC_ROT_LENGTH
    !Config Desc  = The number of years between cuttings for short rotation coppices
    !Config if    = FOREST_MANAGEMENT 
    !Config Def   = undef, 3, 3, 3, 3, 3, 3, 3, 3, undef, undef, undef, undef
    !Config Help  = 
    !Config Units = [years]
    CALL getin_p("SRC_ROT_LENGTH",src_rot_length)

    !Config Key   = SRC_NROTS
    !Config Desc  = Number of rotations before  afinal cut
    !Config if    = FOREST_MANAGEMENT 
    !Config Def   = undef, 10, 10, 10, 10, 10, 10, 10, 10, undef, undef, undef, undef
    !Config Help  = The number of rotations for short rotations coppices
    !               after which the roots/stools are supposed to be exhausted. The 
    !               stool is killed and replanted.
    !Config Units = [-]
    CALL getin_p("SRC_NROTS",src_nrots)

    !Config Key   = FRUIT_ALLOC
    !Config Desc  = Fraction of allocatable carbon that will go to fruit production
    !Config if    = OK_STOMATE
    !Config Def   = (undef, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.2)
    !Config Help  = Guestimates - should be confirmed
    !Config Units = [-] 
    CALL getin_p('FRUIT_ALLOC',fruit_alloc)

    !Config Key   = LABILE_RESERVE 
    !Config Desc  = Determines labile reserve pool 
    !Config if    = OK_STOMATE
    !Config Def   = undef, 60, 30, 60, 60, 30, 60, 10, 10, 2, 2, 2, 2 
    !Config Help  =  
    !Config Units = [-]
    CALL getin_p("LABILE_RESERVE",labile_reserve)

    !Config Key   = EVERGREEN_RESERVE
    !Config Desc  = Fraction of sapwood mass stored in the reserve pool of evergreen trees
    !Config If    = OK_STOMATE
    !Config Def   = undef, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05 
    !Config Help  = 
    !Config Units = [-]  
    CALL getin_p('EVERGREEN_RESERVE',evergreen_reserve)

    !Config Key   = DECIDUOUS_RESERVE
    !Config Desc  = Fraction of sapwood mass stored in the reserve pool
    !Config If    = OK_STOMATE
    !Config Def   = undef, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12
    !Config Help  = Fraction of sapwood mass stored in the reserve pool of 
    !               deciduous trees during the growing season
    !Config Units = [-]  
    CALL getin_p('DECIDUOUS_RESERVE',deciduous_reserve)

    !Config Key   = SENESCENSE_RESERVE
    !Config Desc  = Fraction of sapwood mass stored in the reserve pool of 
    !               deciduous trees during the senescense
    !Config If    = OK_STOMATE
    !Config Def   = undef, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15 
    !Config Help  = Fraction of sapwood mass stored in the reserve pool of 
    !               deciduous trees during the senescense
    !Config Units = [-]  
    CALL getin_p('SENESCENSE_RESERVE',senescense_reserve)

    !Config Key   = ROOT_RESERVE
    !Config Desc  = Fraction of max root biomass which are covered by the carbon reserve
    !Config If    = OK_STOMATE
    !Config Def   = undef, 0.3, 1., 0.3, 0.3, 1., 0.3, 1.,  1., 1., 1., 1., 1. 
    !Config Help  = Fraction of max root biomass which are covered by the carbon reserve.
    !               For evergreens we are happy with 30%, for deciduous we use 100%. In
    !               other words the reserves contain enough C to regrow 100% of the root
    !               biomass for deciduous species.
    !Config Units = [-]  
    CALL getin_p('ROOT_RESERVE',root_reserve)

    !Config Key   = FCN_WOOD
    !Config Desc  = CN of wood for allocation, relative to leaf CN 
    !Config if    = OK_STOMATE
    !Config Def   = undef, .087, .087, .087, .087, .087, .087, .087, .087, .087, .087, .087
    !Config Help  = Comes from Sitch et al 2003 (https://doi.org/10.1046/j.1365-2486.2003.00569.x),
    !Config         although the variables are respresented a bit differntly here.  The sapwood CN
    !Config         ratio in Sitch et al is 330, and that for the leaves and roots is 29. 29/330=0.088.
    !Config Units = [-] 
    CALL getin_p('FCN_WOOD',fcn_wood)  

    !Config Key   = FCN_ROOT
    !Config Desc  = CN roots for allocation, relative to leaf CN 
    !Config if    = OK_STOMATE
    !Config Def   = undef, 0.86, 0.86, 0.86, 0.86, 0.86, 0.86, 0.86, 0.86, 0.86, 0.86, 0.86
    !Config Help  = Comes from Sitch et al 2003 (https://doi.org/10.1046/j.1365-2486.2003.00569.x),
    !Config         although the variables are respresented a bit differntly here.  The root CN
    !Config         ratio in Sitch et al is 29, the same as leaves. 29/29=1.0.  Unclear why
    !Config         the default changed to 0.86.
    !Config Units = [-] 
    CALL getin_p('FCN_ROOT',fcn_root) 

    !Config Key   = BRANCH_RATIO
    !Config Desc  = Share of the sapwood and heartwood that is used for branches
    !Config if    = FOREST_MANAGEMENT 
    !Config Def   = 0.0, 0.38, 0.38, 0.25, 0.38, 0.38, 0.25, 0.38, 0.25, 0.0, 0.0, 0.0, 0.0  
    !Config Help  = 
    !Config Units = [-]
    CALL getin_p("BRANCH_RATIO",branch_ratio)

    !
    ! Recruitment - stomate_prescribe  
    !     
    !Config Key   = RECRUITMENT_PFT  
    !Config Desc  = Logical recruitment flag for each pft
    !Config if    = OK_STOMATE  
    !Config Def   = FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE  
    !Config Help  = Logical recruitment flag for each pft. So PFTs with recruitment or no recruitment can coexist  
    !Config Units = [FLAG]  
    CALL getin_p("RECRUITMENT_PFT",recruitment_pft)  

    !Config Key   = RECRUITMENT_HEIGHT  
    !Config Desc  = Prescribed height for tree recruits (m)   
    !Config if    = OK_STOMATE  
    !Config Def   = undef, 1, 1, 1, 1, 1, 1, 1, 1, 1, undef, undef, undef  
    !Config Help  =
    !Config Units = [m]  
    CALL getin_p("RECRUITMENT_HEIGHT",recruitment_height)  

    !Config Key   = RECRUITMENT_ALPHA  
    !Config Desc  = Intercept of power model relating light and recruitment numbers  
    !Config if    = OK_STOMATE  
    !Config Def   = undef, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0, undef, undef, undef, undef  
    !Config Help  =   
    !Config Units = [-]  
    CALL getin_p("RECRUITMENT_ALPHA",recruitment_alpha)  

    !Config Key   = RECRUITMENT_BETA  
    !Config Desc  = Slope of power model relating light and recruitment numbers  
    !Config if    = OK_STOMATE  
    !Config Def   = undef, 0.8, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, undef, undef, undef, undef  
    !Config Help  =   
    !Config Units = [-]  
    CALL getin_p("RECRUITMENT_BETA",recruitment_beta) 

    !
    ! Mortality
    !
    !Config Key   = DEATH_DISTRIBUTION_FACTOR
    !Config Desc  = Shape parameter for tree mortality
    !Config if    = OK_STOMATE, FUNCTIONAL ALLOCATION
    !Config Def   = undef, 100., 100., 100., 100., 100., 100., 100., 100., undef, undef, undef, undef  
    !Config Help  = The scale factor between the smallest and largest
    !               circ class for tree mortality in stomate_mark_kill.
    !Config Units = [-] 
    CALL getin_p('DEATH_DF',death_distribution_factor) 

    !Config Key   = NPP_RESET_VALUE
    !Config Desc  = The value longterm NPP is reset to npp_reset_value after a non-tree stand dies. 
    !Config if    = OK_STOMATE, FUNCTIONAL ALLOCATION
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef, 500., 500., 500., 500.  
    !Config Help  =
    !Config Units = [gC m-2 y-1] 
    CALL getin_p('NPP_RESET_VALUE',npp_reset_value) 

    !Config Key   = NDYING_YEAR
    !Config Desc  = Number of year for a forest to die
    !Config if    = OK_STOMATE
    !Config Def   = undef, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0,
    !15.0, 15.0, 15.0, 15.0, 15.0  
    !Config Help  = Number of year during which an unmanaged forest will
    ! decay and eventually die after reaching stem density threshold
    !Config Units = [year] 
    CALL getin_p('NDYING_YEAR',ndying_year)

    ! Recruitment - stomate_prescribe  
    !     
    !Config Key   = BEETLE_PFT  
    !Config Desc  = Logical bark beetle mortality flag for each pft
    !Config if    = OK_STOMATE  
    !Config Def   = FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
    !FALSE, FALSE, FALSE, FALSE, FALSE  
    !Config Help  = Logical bark beetle mortality flag for each pft. So PFTs with
    !               bark beetle attack or not bark beetle can coexist.    
    !Config Units = [FLAG]  
    CALL getin_p("BEETLE_PFT",beetle_pft)



       !Config Key   = S_ACTIVITY
       !Config Desc  = beetles population fraction that remain on the stand
       !after beetle departure at the end of an epidemic phase
       !Config if    = OK_STOMATE, OK_PEST
       !Config Def   = undef, undef, undef, 0.5, undef, undef, 0.5, undef,
       !undef, undef, undef, undef, undef
       !Config Help  =
       !Config Units = [-] 
       CALL getin_p('S_ACTIVITY', S_activity)

       !Config Key   = ACT_LIMIT
       !Config Desc  = parameter wich control the feedback of the previous
       !instation on the calculation of the BPI (Beetle pressure index)
       !Config if    = OK_STOMATE, OK_PEST
       !Config Def   = undef, undef, undef, 0.75, undef, undef, 0.75, undef,
       !undef, undef, undef, undef, undef
       !Config Help  =
       !Config Units = [-] 
       CALL getin_p('ACT_LIMIT', act_limit)

       !Config Key   = S_COMPETITION
       !Config Desc  = a parameter for the relationship between rdi and
       !beetle susceptibility 
       !Config if    = OK_STOMATE, OK_PEST
       !Config Def   = undef, undef, undef, 15.5, undef, undef, 15.5, undef,
       !undef, undef, undef, undef, undef
       !Config Help  =
       !Config Units = [-] 
       CALL getin_p('S_COMPETITION',S_competition)

       !Config Key   = RDI_LIMIT
       !Config Desc  = b parameter for the relationship between rdi and
       !beetle susceptibility 
       !Config if    = OK_STOMATE, OK_PEST
       !Config Def   = undef, undef, undef, 0.6, undef, undef,
       !0.6, undef, undef, undef, undef, undef, undef
       !Config Help  =
       !Config Units = [-]
       CALL getin_p('RDI_LIMIT', RDi_limit)

       !Config Key   = S_WEAKNESS
       !Config Desc  = target susceptibility when we come back to endemic
       !Config if    = OK_STOMATE, OK_PEST
       !Config Def   = undef, undef, undef, -5.0, undef, undef,
       !-5.0, undef, undef, undef, undef, undef, undef
       !Config Help  =
       !Config Units = [-]
       CALL getin_p('S_SUSCEPTIBILITY',S_susceptibility)

       !Config Key   = RDI_WEAKNESS
       !Config Desc  = target susceptibility when we come back to endemic
       !Config if    = OK_STOMATE, OK_PEST
       !Config Def   = undef, undef, undef, 0.6, undef, undef,
       !0.6, undef, undef, undef, undef, undef, undef
       !Config Help  =
       !Config Units = [-]
       CALL getin_p('I_RD_SUSCEPTIBILITY',i_rd_susceptibility)

       !Config Key   = S_SHARE
       !Config Desc  = a parameter for the relationship between share and
       !beetle susceptibility 
       !Config if    = OK_STOMATE, OK_PEST
       !Config Def   = undef, undef, undef, 1.5, undef, undef, 15.5, undef,
       !undef, undef, undef, undef, undef
       !Config Help  =
       !Config Units = [-] 
       CALL getin_p('S_SHARE',S_share)

       !Config Key   = SH_LIMIT
       !Config Desc  = b parameter for the relationship between share and
       !beetle susceptibility 
       !Config if    = OK_STOMATE, OK_PEST
       !Config Def   = undef, undef, undef, 0.6, undef, undef,
       !0.6, undef, undef, undef, undef, undef, undef
       !Config Help  =
       !Config Units = [-]
       CALL getin_p('SH_LIMIT',SH_limit)

       !Config Key   = S_DEFENCE
       !Config Desc  = a parameter for the relationship between drought and beetle susceptibility
       !Config if    = OK_STOMATE, OK_PEST
       !Config Def   = undef, undef, undef, -9.5, undef, undef, -9.5, undef, undef, undef, undef, undef, undef
       !Config Help  =
       !Config Units = [-]
       CALL getin_p('S_DEFENCE', S_defence)

       !Config Key   = PWS_LIMIT
       !Config Desc  = b parameter for the relationship between drought and beetle susceptibility
       !Config if    = OK_STOMATE, OK_PEST
       !Config Def   = undef, undef, undef, 0.4, undef, undef, 0.4, undef, undef, undef, undef, undef, undef
       !Config Help  =
       !Config Units = [-]
       CALL getin_p('PWS_LIMIT',PWS_limit)

       !Config Key   = MAX_NWOOD
       !Config Desc  = tune parameter for the relationship between woodleftover and beetle susceptibility
       !Config if    = OK_STOMATE, OK_PEST
       !Config Def   = undef, undef, undef, 1.0, undef, undef, 0.5, undef, undef, undef, undef, undef, undef
       !Config Help  =
       !Config Units = [-]
       CALL getin_p('MAX_NWOOD', max_Nwood)

       !Config Key   = BEETLE_GENERATION_A
       !Config Desc  = a parameter for the calculation of the number of beetle generation per year
       !Config if    = OK_STOMATE
       !Config Def   = undef, undef, undef, 3.307963, undef, undef, 3.307963, undef, undef, undef, undef, undef, undef
       !Config Help  =
       !Config Units = [-] 
       CALL getin_p('BEETLE_GENERATION_A',beetle_generation_a)

       !Config Key   = BEETLE_GENERATION_B
       !Config Desc  = b parameter for the calculation of the number of beetle generation per year
       !Config if    = OK_STOMATE
       !Config Def   = undef, undef, undef, 557.0, undef, undef, 557.0, undef, undef, undef, undef, undef, undef
       !Config Help  = This parameter is the one we have to change across species. It represent the number of DD needed to breed 1 generation 
       !Config Units = [degrees day] 
       CALL getin_p('BEETLE_GENERATION_B',beetle_generation_b)

       !Config Key   = BEETLE_GENERATION_C
       !Config Desc  = c parameter for the calculation of the number of beetle generation per year
       !Config if    = OK_STOMATE
       !Config Def   = undef, undef, undef, 1.980938, undef, undef, 1.980938, undef, undef, undef, undef, undef, undef
       !Config Help  =
       !Config Units = [-] 
       CALL getin_p('BEETLE_GENERATION_C',beetle_generation_c)

       !Config Key   = MIN_TEMP_BEETLE
       !Config Desc  = temperature threshold below which Teff is not calculated
       !Config if    = OK_STOMATE
       !Config Def   = undef, undef, undef, 38.4, undef, undef, 38.4, undef, undef, undef, undef, undef, undef 
       !Config Help  =
       !Config Units = [degree celcius] 
       CALL getin_p('MIN_TEMP_BEETLE',min_temp_beetle)

       !Config Key   = MAX_TEMP_BEETLE
       !Config Desc  = temperature threshold above which Teff is not calculated
       !Config if    = OK_STOMATE
       !Config Def   = undef, undef, undef, 38.4, undef, undef, 38.4, undef, undef, undef, undef, undef, undef
       !Config Help  =
       !Config Units = [ degree celcius] 
       CALL getin_p('MAX_TEMP_BEETLE',max_temp_beetle)

       !Config Key   = OPT_TEMP_BEETLE
       !Config Desc  = a parameter for the calculation of the effective temperature used in beetle phenology
       !Config if    = OK_STOMATE
       !Config Def   = undef, undef, undef, 30.3, undef, undef, 30.3, undef, undef, undef, undef, undef, undef
       !Config Help  =
       !Config Units = [-] 
       CALL getin_p('OPT_TEMP_BEETLE',opt_temp_beetle)

       !Config Key   = EFF_TEMP_BEETLE_A
       !Config Desc  = a parameter for the calculation of the effective temperature used in beetle phenology
       !Config if    = OK_STOMATE
       !Config Def   = undef, undef, undef, 0.02876507, undef, undef, 0.02876507, undef, undef, undef, undef, undef, undef
       !Config Help  =
       !Config Units = [-] 
       CALL getin_p('EFF_TEMP_BEETLE_A',eff_temp_beetle_a)

       !Config Key   = EFF_TEMP_BEETLE_B
       !Config Desc  = b parameter for the calculation of the effective temperature used in beetle phenology
       !Config if    = OK_STOMATE
       !Config Def   = undef, undef, undef, 40.9958913, undef, undef, 40.9958913, undef, undef, undef, undef, undef, undef
       !Config Help  =
       !Config Units = [-] 
       CALL getin_p('EFF_TEMP_BEETLE_B',eff_temp_beetle_b)

       !Config Key   = EFF_TEMP_BEETLE_C
       !Config Desc  = c parameter for the calculation of the effective temperature used in beetle phenology
       !Config if    = OK_STOMATE
       !Config Def   = undef, undef, undef, 3.5922336, undef, undef, 3.5922336, undef, undef, undef, undef, undef, undef
       !Config Help  =
       !Config Units = [-] 
       CALL getin_p('EFF_TEMP_BEETLE_C',eff_temp_beetle_c)

       !Config Key   = EFF_TEMP_BEETLE_D
       !Config Desc  = d parameter for the calculation of the effective temperature used in beetle phenology
       !Config if    = OK_STOMATE
       !Config Def   = undef, undef, undef, 1.24657367, undef, undef, 1.24657367, undef, undef, undef, undef, undef, undef
       !Config Help  =
       !Config Units = [-] 
       CALL getin_p('EFF_TEMP_BEETLE_D',eff_temp_beetle_d)

       !Config Key   = DIAPAUSE_THRES_DAYLENGTH
       !Config Desc  = daylength in hour above which bark beetle start diapause
       !Config if    = OK_STOMATE
       !Config Def   = undef, undef, undef, 14.5, undef, undef, 14.5, undef, undef, undef, undef, undef, undef 
       !Config Help  =
       !Config Units = [hour] 
       CALL getin_p('DIAPAUSE_THRES_DAYLENGTH',diapause_thres_daylength)

       !Config Key   = S_MASSATTACK
       !Config Desc  = ""
       !Config if    = OK_STOMATE
       !Config Def   = undef, undef, undef, 15.5, undef, undef, 15.5, undef,
       !undef, undef, undef, undef, undef 
       !Config Help  =
       !Config Units = [-] 
       CALL getin_p('S_MASSATTACK', S_massattack)

       !Config Key   = BP_LIMIT
       !Config Desc  = ""
       !Config if    = OK_STOMATE
       !Config Def   = undef, undef, undef, 0.5, undef, undef, 0.5, undef,
       !undef, undef, undef, undef, undef 
       !Config Help  =
       !Config Units = [-] 
       CALL getin_p('BP_LIMIT', BP_limit)

    !
    ! Windthrow - stomate
    !
    IF (ok_windthrow) THEN

       !Config Key   = STREAMLINING_C_LEAF
       !Config Desc  = streamlining parameter for crown with leaves
       !Config if    = OK_STOMATE, OK_WINDTHROW
       !Config Def   = undef, 2.34, 2.34, 2.70, 2.66, 2.34, 2.71, 2.15, 3.07, undef, undef, undef, undef
       !Config Help  =
       !Config Units = [-] 
       CALL getin_p('STREAMLINING_C_LEAF',streamlining_c_leaf)

       !Config Key   = STREAMLINING_C_LEAFLESS
       !Config Desc  = streamlining parameter for crown without leaves
       !Config if    = OK_STOMATE, OK_WINDTHROW
       !Config Def   = undef, 2.34, 2.34, 2.70, 2.66, 2.34, 2.71, 2.15, 3.07, undef, undef, undef, undef
       !Config Help  =
       !Config Units = [-]
       CALL getin_p('STREAMLINING_C_LEAFLESS',streamlining_c_leafless)

       !Config Key   = STREAMLINING_N_LEAF
       !Config Desc  = streamlining parameter for crown with leaves
       !Config if    = OK_STOMATE, OK_WINDTHROW
       !Config Def   = undef, 0.88, 0.88, 0.64, 0.85, 0.88, 0.63, 0.88, 0.75, undef, undef, undef, undef
       !Config Help  =
       !Config Units = [-]
       CALL getin_p('STREAMLINING_N_LEAF',streamlining_n_leaf)

       !Config Key   = STREAMLINING_N_LEAFLESS
       !Config Desc  = streamlining parameter for crown without leaves
       !Config if    = OK_STOMATE, OK_WINDTHROW
       !Config Def   = undef, 0.88, 0.88, 0.64, 0.85, 0.88, 0.63, 0.88, 0.75, undef, undef, undef, undef
       !Config Help  =
       !Config Units = [-]
       CALL getin_p('STREAMLINING_N_LEAFLESS',streamlining_n_leafless)
       
       !Config Key   = MODULUS_RUPTURE
       !Config Desc  = Modulus of rupture
       !Config if    = OK_STOMATE, OK_WINDTHROW
       !Config Def   = undef, 6.23E7, 6.23E7, 4.13E7, 5.90E7, 6.23E7, 4.10E7, 6.27E7, 5.30E7, undef, undef, undef, undef
       !Config Help  = The measure of a species’ strength before rupture 
       !               when being bent. Used in the calculation of the critical
       !               wind speed according to the GALES (Hale et al. 2015) model. 
       !               IMPORTANT: greenwood values are used and not the more frequently 
       !               available drywood modulus of rupture.
       !Config Units = [Pa]
       CALL getin_p('MODULUS_RUPTURE',modulus_rupture)

       !Config Key   = F_KNOT
       !Config Desc  = Knot factor 
       !Config if    = OK_STOMATE, OK_WINDTHROW
       !Config Def   = undef, 1.0, 1.0, 0.87, 1.0, 1.0, 0.88, 1.0, 0.85, undef, undef, undef, undef
       !Config Help  = This modifier represents the knot in the wood, and hence 
       !               the decrease in structural strength. Used in the calculation 
       !               of the critical wind speed according to the GALES 
       !               (Hale et al. 2015) model.
       !Config Units = [unitless]
       CALL getin_p('F_KNOT',f_knot)

       !Config Key   = GREEN_DENSITY
       !Config Desc  = Green density of the tree
       !Config if    = OK_STOMATE, OK_WINDTHROW
       !Config Def   = undef, 1007, 1007, 985, 1060, 1007, 990, 968, 900, undef, undef, undef, undef
       !Config Help  = 
       !Config Units = [kg.m-3]
       CALL getin_p('GREEN_DENSITY',green_density)

       !Config Key   = OV_FD_SHALLOW 
       !Config Desc  = Regression coefficient for overturning in free draining and shallow soil type
       !Config if    = OK_STOMATE, OK_WINDTHROW
       !Config Def   = undef, 175.3, 175.3, 134.7, 198.5, 175.3, 132.6, 152.0, 145.2, undef, undef, undef, undef 
       !Config Help  = Values derived from generic soil types (free_draining mineral
       !               soils; Gleyed mineral soils; Peaty mineral soils; Deep peats) 
       !               and the soil depth (shallow, deep, average)
       !Config Units = [Nm/kg]
       CALL getin_p('OV_FD_SHALLOW',overturning_free_draining_shallow)

       !Config Key   = OV_FD_SHALLOW_LESS
       !Config Desc  = Regression coefficient for overturning in free draining and shallow soil type leafless
       !Config if    = OK_STOMATE, OK_WINDTHROW
       !Config Def   = undef, 175.3, 175.3, 134.7, 198.5, 175.3, 132.6, 152.0, 145.2, undef, undef, undef, undef
       !Config Help  = Values derived from generic soil types (free_draining mineral
       !               soils; Gleyed mineral soils; Peaty mineral soils; Deep peats) 
       !               and the soil depth (shallow, deep, average)
       !Config Units = [Nm/Kg]
       CALL getin_p('OV_FD_SHALLOW_LESS',overturning_free_draining_shallow_leafless)

       !Config Key   = OV_FD_DEEP
       !Config Desc  = Regression coefficient for overturning in free draining and deep soil type
       !Config if    = OK_STOMATE, OK_WINDTHROW
       !Config Def   = undef, 203.8, 203.8, 157.2, 230.8, 230.8, 154.8, 176.7, 169.4, undef, undef, undef, undef 
       !Config Help  = Values derived from generic soil types (free_draining mineral
       !               soils; Gleyed mineral soils; Peaty mineral soils; Deep peats) 
       !               and the soil depth (shallow, deep, average)
       !Config Units = [Nm/Kg]
       CALL getin_p('OV_FD_DEEP',overturning_free_draining_deep)

       !Config Key   = OV_FD_DEEP_LESS
       !Config Desc  = Regression coefficient for overturning in free draining and deep soil type leafless
       !Config if    = OK_STOMATE, OK_WINDTHROW
       !Config Def   = undef, 203.8, 203.8, 157.2, 230.8, 230.8, 154.8, 176.7, 169.4, undef, undef, undef, undef  
       !Config Help  = Values derived from generic soil types (free_draining mineral
       !               soils; Gleyed mineral soils; Peaty mineral soils; Deep peats) 
       !               and the soil depth (shallow, deep, average)
       !Config Units = [Nm/Kg]
       CALL getin_p('OV_FD_DEEP_LESS',overturning_free_draining_deep_leafless)

       !Config Key   = OV_FD_AVERAGE
       !Config Desc  = Regression coefficient for overturning in free draining and medium soil type  
       !Config if    = OK_STOMATE, OK_WINDTHROW
       !Config Def   = undef, 178.7, 178.7, 137.8, 202.4, 178.7, 135.7, 155.0, 148.6, undef, undef, undef, undef  
       !Config Help  = Values derived from generic soil types (free_draining mineral
       !               soils; Gleyed mineral soils; Peaty mineral soils; Deep peats) 
       !               and the soil depth (shallow, deep, average)
       !Config Units = [Nm/Kg]
       CALL getin_p('OV_FD_AVERAGE',overturning_free_draining_average)

       !Config Key   = OV_FD_AVERAGE_LESS 
       !Config Desc  = Regression coefficient for overturning in free draining and medium soil type leafless 
       !Config if    = OK_STOMATE, OK_WINDTHROW
       !Config Def   = undef, 178.7, 178.7, 137.8, 202.4, 178.7, 135.7, 155.0, 148.6, undef, undef, undef, undef 
       !Config Help  = Values derived from generic soil types (free_draining mineral
       !               soils; Gleyed mineral soils; Peaty mineral soils; Deep peats) 
       !               and the soil depth (shallow, deep, average)
       !Config Units = [Nm/Kg]
       CALL getin_p('OV_FD_AVERAGE_LESS',overturning_free_draining_average_leafless)

       !Config Key   = OV_GLEYED_SHALLOW
       !Config Desc  = Regression coefficient for overturning in gleyed and shallow soil type
       !Config if    = OK_STOMATE, OK_WINDTHROW
       !Config Def   = undef, 155.4, 155.4, 119.4, 176.0, 155.4, 117.6, 134.8, 128.7, undef, undef, undef, undef 
       !Config Help  = Values derived from generic soil types (free_draining mineral
       !               soils; Gleyed mineral soils; Peaty mineral soils; Deep peats) 
       !               and the soil depth (shallow, deep, average)
       !Config Units = [Nm/Kg]
       CALL getin_p('OV_GLEYED_SHALLOW',overturning_gleyed_shallow)

       !Config Key   = OV_GLEYED_SHALLOW_LESS
       !Config Desc  = Regression coefficient for overturning in gleyed and shallow soil type leafless
       !Config if    = OK_STOMATE, OK_WINDTHROW
       !Config Def   = undef, 155.4, 155.4, 119.4, 176.0, 155.4, 117.6, 134.8, 128.7, undef, undef, undef, undef  
       !Config Help  = Values derived from generic soil types (free_draining mineral
       !               soils; Gleyed mineral soils; Peaty mineral soils; Deep peats) 
       !               and the soil depth (shallow, deep, average)
       !Config Units = [Nm/Kg]
       CALL getin_p('OV_GLEYED_SHALLOW_LESS',overturning_gleyed_shallow_leafless)

       !Config Key   = OV_GLEYED_DEEP
       !Config Desc  = Regression coefficient for overturning in gleyed and deep soil type
       !Config if    = OK_STOMATE, OK_WINDTHROW
       !Config Def   = undef, 180.6, 180.6, 139.3, 204.6, 180.6, 137.2, 156.7, 150.2, undef, undef, undef, undef 
       !Config Help  = Values derived from generic soil types (free_draining mineral
       !               soils; Gleyed mineral soils; Peaty mineral soils; Deep peats) 
       !               and the soil depth (shallow, deep, average)
       !Config Units = [Nm/Kg
       CALL getin_p('OV_GLEYED_DEEP',overturning_gleyed_deep)

       !Config Key   = OV_GLEYED_DEEP_LESS 
       !Config Desc  = Regression coefficient for overturning in gleyed and deep soil type leafless
       !Config if    = OK_STOMATE, OK_WINDTHROW
       !Config Def   = undef, 180.6, 180.6, 139.3, 204.6, 180.6, 137.2, 156.7, 150.2, undef, undef, undef, undef  
       !Config Help  = Values derived from generic soil types (free_draining mineral
       !               soils; Gleyed mineral soils; Peaty mineral soils; Deep peats) 
       !               and the soil depth (shallow, deep, average)
       !Config Units = [Nm/Kg]
       CALL getin_p('OV_GLEYED_DEEP_LESS',overturning_gleyed_deep_leafless)

       !Config Key   = OV_GLEYED_AVERAGE
       !Config Desc  = Regression coefficient for overturning in gleyed and medium soil type
       !Config if    = OK_STOMATE, OK_WINDTHROW
       !Config Def   = undef, 158.5, 158.5, 122.2, 179.5, 158.5, 120.3, 137.4, 131.7, undef, undef, undef, undef  
       !Config Help  = Values derived from generic soil types (free_draining mineral
       !               soils; Gleyed mineral soils; Peaty mineral soils; Deep peats) 
       !               and the soil depth (shallow, deep, average)
       !Config Units = [Nm/Kg]
       CALL getin_p('OV_GLEYED_AVERAGE',overturning_gleyed_average)

       !Config Key   = OV_GLEYED_AVERAGE_LESS
       !Config Desc  = Regression coefficient for overturning in gleyed and medium soil type leafless
       !Config if    = OK_STOMATE, OK_WINDTHROW
       !Config Def   = undef, 158.5, 158.5, 122.2, 179.5, 158.5, 120.3, 137.4, 131.7, undef, undef, undef, undef  
       !Config Help  = Values derived from generic soil types (free_draining mineral
       !               soils; Gleyed mineral soils; Peaty mineral soils; Deep peats) 
       !               and the soil depth (shallow, deep, average)
       !Config Units = [Nm/Kg]
       CALL getin_p('OV_GLEYED_AVERAGE_LESS',overturning_gleyed_average_leafless)

       !Config Key   = OV_PEATY_SHALLOW
       !Config Desc  = Regression coefficient for overturning in peaty and shallow soil type
       !Config if    = OK_STOMATE, OK_WINDTHROW
       !Config Def   = undef, 169.7, 169.7, 130.4, 192.2, 169.7, 128.4, 147.2, 140.6, undef, undef, undef, undef  
       !Config Help  = Values derived from generic soil types (free_draining mineral
       !               soils; Gleyed mineral soils; Peaty mineral soils; Deep peats) 
       !               and the soil depth (shallow, deep, average)
       !Config Units = [Nm/Kg]
       CALL getin_p('OV_PEATY_SHALLOW',overturning_peaty_shallow)

       !Config Key   = OV_PEATY_SHALLOW_LESS
       !Config Desc  = Regression coefficient for overturning in peaty and shallow soil type leafless
       !Config if    = OK_STOMATE, OK_WINDTHROW
       !Config Def   = undef, 169.7, 169.7, 130.4, 192.2, 169.7, 128.4, 147.2, 140.6, undef, undef, undef, undef  
       !Config Help  = Values derived from generic soil types (free_draining mineral
       !               soils; Gleyed mineral soils; Peaty mineral soils; Deep peats) 
       !               and the soil depth (shallow, deep, average)
       !Config Units = [Nm/Kg]
       CALL getin_p('OV_PEATY_SHALLOW_LESS',overturning_peaty_shallow_leafless)

       !Config Key   = OV_PEATY_DEEP
       !Config Desc  = Regression coefficient for overturning in peaty and deep soil type
       !Config if    = OK_STOMATE, OK_WINDTHROW
       !Config Def   = undef, 191.4, 191.4, 152.1, 223.5, 191.4, 141.9, 159.2, 164.0, undef, undef, undef, undef   
       !Config Help  = Values derived from generic soil types (free_draining mineral
       !               soils; Gleyed mineral soils; Peaty mineral soils; Deep peats) 
       !               and the soil depth (shallow, deep, average)
       !Config Units = [Nm/Kg]
       CALL getin_p('OV_PEATY_DEEP',overturning_peaty_deep)

       !Config Key   = OV_PEATY_DEEP_LESS
       !Config Desc  = Regression coefficient for overturning in peaty and deep soil type leafless  
       !Config if    = OK_STOMATE, OK_WINDTHROW
       !Config Def   = undef, 191.4, 191.4, 152.1, 223.5, 191.4, 141.9, 159.2, 164.0, undef, undef, undef, undef  
       !Config Help  = Values derived from generic soil types (free_draining mineral
       !               soils; Gleyed mineral soils; Peaty mineral soils; Deep peats) 
       !               and the soil depth (shallow, deep, average)
       !Config Units = [Nm/Kg]
       CALL getin_p('OV_PEATY_DEEP_LESS',overturning_peaty_deep_leafless)

       !Config Key   = OV_PEATY_AVERAGE
       !Config Desc  = Regression coefficient for overturning in peaty and medium soil type 
       !Config if    = OK_STOMATE, OK_WINDTHROW
       !Config Def   = undef, 178.9, 178.9, 133.4, 195.9, 178.9, 131.4, 162.0, 143.8, undef, undef, undef, undef 
       !Config Help  = Values derived from generic soil types (free_draining mineral
       !               soils; Gleyed mineral soils; Peaty mineral soils; Deep peats) 
       !               and the soil depth (shallow, deep, average)
       !Config Units = [Nm/Kg]
       CALL getin_p('OV_PEATY_AVERAGE',overturning_peaty_average)

       !Config Key   = OV_PEATY_AVERAGE_LESS
       !Config Desc  = Regression coefficient for overturning in peaty and medium soil type leafless 
       !Config if    = OK_STOMATE, OK_WINDTHROW
       !Config Def   = undef, 178.9, 178.9, 133.4, 195.9, 178.9, 131.4, 162.0, 143.8, undef, undef, undef, undef   
       !Config Help  = Values derived from generic soil types (free_draining mineral
       !               soils; Gleyed mineral soils; Peaty mineral soils; Deep peats) 
       !               and the soil depth (shallow, deep, average)
       !Config Units = [Nm/Kg]
       CALL getin_p('OV_PEATY_AVERAGE_LESS',overturning_peaty_average_leafless)

       !Config Key   = OV_PEAT_SHALLOW
       !Config Desc  = Regression coefficient for overturning in shallow peat soil type
       !Config if    = OK_STOMATE, OK_WINDTHROW
       !Config Def   = undef, 193.0, 193.0, 148.3, 218.6, 193.0, 146.0, 167.4, 159.9, undef, undef, undef, undef  
       !Config Help  = Values derived from generic soil types (free_draining mineral
       !               soils; Gleyed mineral soils; Peaty mineral soils; Deep peats) 
       !               and the soil depth (shallow, deep, average)
       !Config Units = [Nm/Kg]
       CALL getin_p('OV_PEAT_SHALLOW',overturning_peat_shallow)

       !Config Key   = OV_PEAT_SHALLOW_LESS
       !Config Desc  = Regression coefficient for overturning in shallow peat soil leafless
       !Config if    = OK_STOMATE, OK_WINDTHROW
       !Config Def   = undef, 193.0, 193.0, 148.3, 218.6, 193.0, 146.0, 167.4, 159.9, undef, undef, undef, undef 
       !Config Help  = Values derived from generic soil types (free_draining mineral
       !               soils; Gleyed mineral soils; Peaty mineral soils; Deep peats) 
       !               and the soil depth (shallow, deep, average)
       !Config Units = [Nm/Kg]
       CALL getin_p('OV_PEAT_SHALLOW_LESS',overturning_peat_shallow_leafless)

       !Config Key   = OV_PEAT_DEEP
       !Config Desc  = Regression coefficient for overturning in deep peat soil
       !Config if    = OK_STOMATE, OK_WINDTHROW
       !Config Def   = undef, 224.4, 224.4, 173.1, 254.2, 224.4, 170.4, 194.7, 186.6, undef, undef, undef, undef
       !Config Help  = Values derived from generic soil types (free_draining mineral
       !               soils; Gleyed mineral soils; Peaty mineral soils; Deep peats) 
       !               and the soil depth (shallow, deep, average)
       !Config Units = [Nm/Kg]
       CALL getin_p('OV_PEAT_DEEP',overturning_peat_deep)

       !Config Key   = OV_PEAT_DEEP_LESS
       !Config Desc  = Regression coefficient for overturning in deep peat soil leafless
       !Config if    = OK_STOMATE, OK_WINDTHROW
       !Config Def   = undef, 224.4, 224.4, 173.1, 254.2, 224.4, 170.4, 194.7, 186.6, undef, undef, undef, undef
       !Config Help  = Values derived from generic soil types (free_draining mineral
       !               soils; Gleyed mineral soils; Peaty mineral soils; Deep peats) 
       !               and the soil depth (shallow, deep, average)
       !Config Units = [Nm/Kg]
       CALL getin_p('OV_PEAT_DEEP_LESS',overturning_peat_deep_leafless)

       !Config Key   = OV_PEAT_AVERAGE
       !Config Desc  = Regression coefficient for overturning in medium peat soil
       !Config if    = OK_STOMATE, OK_WINDTHROW
       !Config Def   = undef, 196.9, 196.9, 151.8, 223.0, 196.9, 149.4, 170.8, 163.6, undef, undef, undef, undef  
       !Config Help  = Values derived from generic soil types (free_draining mineral
       !               soils; Gleyed mineral soils; Peaty mineral soils; Deep peats) 
       !               and the soil depth (shallow, deep, average)
       !Config Units = [Nm/Kg]
       CALL getin_p('OV_PEAT_AVERAGE',overturning_peat_average)

       !Config Key   = OV_PEAT_AVERAGE_LESS
       !Config Desc  = Regression coefficient for overturning in medium peat soil leafless
       !Config if    = OK_STOMATE, OK_WINDTHROW
       !Config Def   = undef, 196.9, 196.9, 151.8, 223.0, 196.9, 149.4, 170.8, 163.6, undef, undef, undef, undef 
       !Config Help  = Values derived from generic soil types (free_draining mineral
       !               soils; Gleyed mineral soils; Peaty mineral soils; Deep peats) 
       !               and the soil depth (shallow, deep, average)
       !Config Units = [Nm/Kg]
       CALL getin_p('OV_PEAT_AVERAGE_LESS',overturning_peat_average_leafless)

       !Config Key   = MDF
       !Config Desc  = Maximum damage rate away from the forest edge
       !Config if    = OK_STOMATE, OK_WINDTHROW
       !Config Def   = undef, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8 
       !Config Help  = A tuning parameter for determining wind damage rate/level
       !Config Units = [unitless]
       CALL getin_p('MAX_DAMAGE_FURTHER',max_damage_further)

       !Config Key   = MDC
       !Config Desc  = Maximum damage rate nearby the forest edge
       !Config if    = OK_STOMATE, OK_WINDTHROW
       !Config Def   = undef, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8 
       !Config Help  = A tuning parameter for determining wind damage rate/level
       !Config Units = [unitless]
       CALL getin_p('MAX_DAMAGE_CLOSER',max_damage_closer)

       !Config Key   = SFF
       !Config Desc  = Scaling factor for maximum damage rate away from the forest edge 
       !Config if    = OK_STOMATE, OK_WINDTHROW
       !Config Def   = undef, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8 
       !Config Help  = A tuning parameter for determining wind damage rate/level
       !Config Units = [unitless]
       CALL getin_p('SFACTOR_FURTHER',sfactor_further)

       !Config Key   = SFC
       !Config Desc  = Scaling factor for maximum damage rate nearby the forest edge
       !Config if    = OK_STOMATE, OK_WINDTHROW
       !Config Def   = undef, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8 
       !Config Help  = A tuning parameter for determining wind damage rate/level
       !Config Units = [unitless]
       CALL getin_p('SFACTOR_CLOSER',sfactor_closer)

    ENDIF
    
    !
    ! Phenology : Senescence
    !
    !
    !Config Key   = LEAFFALL
    !Config Desc  = length of death of leaves, tabulated 
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, 10., undef, undef, 10., undef, 10., 10., 10., 10., 10., 10. 
    !Config Help  =
    !Config Units = [days]
    CALL getin_p('LEAFFALL',leaffall)

    !Config Key   = PRESENESCENCE_RATIO
    !Config Desc  = The ratio of maintenance respiration to gpp beyond which presenescence 
    !               stage of plant phenology is declared to begin.
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, .6, undef, undef, .6, undef, .6, .6, .6, .6, .6, .6
    !Config Help  =
    !Config Units = [0-1, unitless]
    CALL getin_p('PRESENESCENCE_RATIO',presenescence_ratio)

    !Config Key   = SENESCENCE_RATIO
    !Config Desc  = The ratio of maintenance respiration to gpp beyond which senescence 
    !               stage of plant phenology is declared to begin.
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, 1, undef, undef, 1, undef, 0.87, 0.89, 0.75, 1, 1, 1 
    !Config Help  =
    !Config Units = [0-1, unitless]
    CALL getin_p('SENESCENCE_RATIO',senescence_ratio)

    !Config Key   = SENESCENCE_HUM
    !Config Desc  = critical relative moisture availability for senescence
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, .3, undef, undef, undef, undef, undef, undef, .2, .2, .3, .2 
    !Config Help  =
    !Config Units = [-] 
    CALL getin_p('SENESCENCE_HUM',senescence_hum)

    !Config Key   = NOSENESCENCE_HUM
    !Config Desc  = relative moisture availability above which there is no humidity-related senescence
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, .8, undef, undef, undef, undef, undef, undef, .3, .3, .3, .3 
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('NOSENESCENCE_HUM',nosenescence_hum) 

    !Config Key   = MAX_TURNOVER_TIME
    !Config Desc  = maximum turnover time for grasse
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef,  80.,  80., 80., 80. 
    !Config Help  =
    !Config Units = [days]
    CALL getin_p('MAX_TURNOVER_TIME',max_turnover_time)

    !Config Key   = MIN_TURNOVER_TIME
    !Config Desc  = minimum turnover time for grasse 
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef, 10., 10., 10., 10. 
    !Config Help  =
    !Config Units = [days]
    CALL getin_p('MIN_TURNOVER_TIME',min_turnover_time)

    !Config Key   = RECYCLE_LEAF
    !Config Desc  = Fraction of N leaf that is recycled when leaves are senescent 
    !Config if    = OK_STOMATE
    !Config Def   = undef, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('RECYCLE_LEAF',recycle_leaf)

    !Config Key   = RECYCLE_ROOT
    !Config Desc  = Fraction of N root that is recycled when roots are senescent 
    !Config if    = OK_STOMATE
    !Config Def   = undef, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2 
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('RECYCLE_ROOT',recycle_root)

    !Config Key   = MIN_LEAF_AGE_FOR_SENESCENCE
    !Config Desc  = minimum leaf age to allow senescence g
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, 90., undef, undef, 90., undef, 60., 60., 30., 30., 30., 30.
    !Config Help  =
    !Config Units = [days] 
    CALL getin_p('MIN_LEAF_AGE_FOR_SENESCENCE',min_leaf_age_for_senescence)

    !Config Key   = SENESCENCE_TEMP_C
    !Config Desc  = critical temperature for senescence (C), constant c of aT^2+bT+c, tabulated
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, undef, undef, undef, 16., undef, 14., 13., 5, 5., 12., 13.
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('SENESCENCE_TEMP_C',senescence_temp_c)

    !Config Key   = SENESCENCE_TEMP_B
    !Config Desc  = critical temperature for senescence (C), constant b of aT^2+bT+c ,tabulated
    !Config if    = OK_STOMATE 
    !Config Def   = undef, undef, undef, undef, undef, 0., undef, 0., 0., .1, 0., 0., 0.
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('SENESCENCE_TEMP_B',senescence_temp_b)

    !Config Key   = SENESCENCE_TEMP_A
    !Config Desc  = critical temperature for senescence (C), constant a of aT^2+bT+c , tabulated
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, undef, undef, undef, 0., undef, 0., 0.,.00375, 0., 0., 0. 
    !Config Help  =
    !Config Units = [-] 
    CALL getin_p('SENESCENCE_TEMP_A',senescence_temp_a)

    !Config Key   = GDD_SENESCENCE
    !Config Desc  = minimum gdd to allow senescence of crops  
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, 950., 4000.
    !Config Help  =
    !Config Units = [days] 
    CALL getin_p("GDD_SENESCENCE", gdd_senescence)

    !Config Key   = PRE_GDD
    !Config Desc  = minimum gdd to allow senescence of crops  
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef,
    !undef, undef, undef, 60., 60.
    !Config Help  =
    !Config Units = [days] 
    CALL getin_p("PRE_GDD", pre_gdd)

    !Config Key   = PRE_HARV
    !Config Desc  = minimum gdd to allow senescence of crops  
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef,
    !undef, undef, undef, 60., 60.
    !Config Help  =
    !Config Units = [days] 
    CALL getin_p("PRE_HARV", pre_harv)

    !Config Key   = GR_GT
    !Config Desc  = minimum gdd to allow senescence of crops  
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef,
    !undef, undef, undef, 60., 60.
    !Config Help  =
    !Config Units = [days] 
    CALL getin_p("GR_GT", gr_gt)

    !Config Key   = ALWAYS_INIT
    !Config Desc  = Take carbon from atmosphere if carbohydrate reserve too small
    !Config if    = OK_STOMATE
    !Config Def   = y, y, y, y, y, y, y, y, y, y, n, y, y
    !Config Help  =
    !Config Units = [BOOLEAN]
    CALL getin_p('ALWAYS_INIT',always_init)

    !
    ! N cycle
    !
    !Config Key   = MAX_SOIL_N_BNF
    !Config Desc  = Value of total N (NH4+NO3) above which we stop adding N via BNF (gN/m**2)  
    !Config if    = OK_STOMATE
    !Config Def   = 0.0, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 2., 2., 2., 2.
    !Config Help  =
    !Config Units = [gN/m**2] 
    CALL getin_p("MAX_SOIL_N_BNF", max_soil_n_bnf)

    !Config Key   = MANURE_PFTWEIGHT
    !Config Desc  = Weight of the distribution of manure over the PFT surface
    !Config if    = OK_STOMATE
    !Config Def   = 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 1., 1., 1. 
    !Config Help  =
    !Config Units = [gC/gN] 
    CALL getin_p("MANURE_PFTWEIGHT", manure_pftweight)

    !
    ! CROPLAND MANAGEMENT
    !
    !Config Key   = HARVEST_RATIO
    !Config Desc  = Share of biomass that is harvested
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, 0.5, 0.5
    !Config Help  = Share of biomass that is harvested. This residual equal to 1 minus harvest_ratio
    !Config Units = [unitless] 
    CALL getin_p("HARVEST_RATIO", harvest_ratio)
    
    !
    ! DGVM
    !

    !Config Key   = RESIDENCE_TIME
    !Config Desc  = residence time of trees
    !Config if    = OK_DGVM and NOT(LPJ_GAP_CONST_MORT)
    !Config Def   = undef, 30.0, 30.0, 40.0, 40.0, 40.0, 80.0, 80.0, 80.0, 0.0, 0.0, 0.0, 0.0 
    !Config Help  =
    !Config Units = [years]
    CALL getin_p('RESIDENCE_TIME',residence_time)

    !Config Key   = TMIN_CRIT
    !Config Desc  = critical tmin, tabulated
    !Config if    = OK_STOMATE
    !Config Def   = undef,  0.0, 0.0, -30.0, -14.0, -30.0, -45.0, -45.0, undef, undef, undef, undef, undef
    !Config Help  = 
    !Config Units = [C]
    CALL getin_p('TMIN_CRIT',tmin_crit)

    !Config Key   = TCM_CRIT
    !Config Desc  = critical tcm, tabulated 
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, undef, 5.0, 15.5, 15.5, -8.0, -8.0, -8.0, undef, undef, undef, undef
    !Config Help  =
    !Config Units = [C]
    CALL getin_p('TCM_CRIT',tcm_crit)

    !
    ! SOIL DECOMPOSITION
    !
    !Config Key   = SOM_TURN_IPASSIVE
    !Config Desc  = Turnover of the passive pool (year-1) 
    !Config if    = OK_STOMATE
    !Config Def   = .0025, .0025, .0025, .0025, .0025, .0025, .0025, .0025, .0025, .0025, .0025, .0025, .0025
    !Config Help  = Turnover of the passive soil pool. Conceptually it is
    !               more correct to have this as a constant for all PFTs,
    !               but changing it at the PFT level could help to tune the
    !               total soil carbon
    !Config Units = [year-1] 
    CALL getin_p('SOM_TURN_IPASSIVE',som_turn_ipassive)

    ! Age classes
    ! I want to create a temporary array that indicates which "real" PFT starts
    ! on which index.  This could probably be put somewhere else, but this
    ! routine is only called once a year and this loop is not expensive.
    start_index(:)=-1
    nagec_pft(:)=-1
    DO ivma=1,nvmap
       ! The start index is just the first place we find this real PFT.
       DO ivm=1,nvm
          IF(agec_group(ivm) .EQ. ivma)THEN
             start_index(ivma)=ivm
             ! It is possible that not all forests will have multiple age 
             ! classes.  For example, the species might have age classes 
             ! but metaclasses (running outside Europe) might not. Let's 
             ! check to see how many age classes each PFT has. Right now, 
             ! the only options are 1 or nagec, but this could be changed 
             ! without too much difficulty.
             IF((ivm+nagec-1) .LE. nvm)THEN
                ! This first if loop prevents an out of bounds error.
                ! Guillaume M. -- .LE. and not .LT.: when the LAST group of the list
                ! carries nagec classes, ivm+nagec-1 is exactly nvm and agec_group(nvm)
                ! is in bounds. With .LT. that group silently got nagec_pft=1, so the
                ! allow_mgmt_cut guard (last slot of the group) fell on its CLASS 1,
                ! which then performed rotation cuts outside the terminal class.
                IF(agec_group(ivm+nagec-1) == ivma)THEN
                   nagec_pft(ivma)=nagec
                ELSE
                   nagec_pft(ivma)=1
                ENDIF
             ELSE
                nagec_pft(ivma)=1
             ENDIF
             EXIT
          ENDIF
       ENDDO
    ENDDO
    ! Check to see if the calculation worked and we found indices for all of them.
    DO ivma=1,nvmap
       IF(start_index(ivma) .LT. 0)THEN
          WRITE(numout,*) 'Could not find a start index for one age class group!'
          WRITE(numout,*) 'Check the input file and ',&
               'make sure the following ivma appears in agec_group'
          WRITE(numout,*) 'ivma,nvmap',ivma,nvmap
          WRITE(numout,*) 'agec_group',agec_group(:)
          CALL ipslerr_p (3,'pft_parameters','pft_parameters_alloc','','')
       ENDIF
    ENDDO

  END SUBROUTINE config_stomate_pft_parameters


!! ================================================================================================================================
!! SUBROUTINE   : derive_age_class_bounds
!!
!>\BRIEF         Derives the per-PFT age class bounds from each PFT clearcut diameter.
!!
!! DESCRIPTION : design/MODULE_DESIGN_AGE_CLASS_BOUNDS_PFT.md (Marie 2026).
!!
!! age_class_bound(c,jv) = AGE_CLASS_FRAC(c) * largest_tree_dia(jv)   for c = 1..nagec-1
!! age_class_bound(nagec,jv) = unchanged open upper bound (5 m)
!!
!! WHY A SEPARATE ROUTINE, CALLED FROM control.f90 AND NOT FROM constantes.f90 :
!! LARGEST_TREE_DIA is only read from the namelist in config_stomate_pft_parameters
!! (control.f90:1465), i.e. AFTER config_stomate_parameters (control.f90:1459) where
!! age_class_bound is allocated. Deriving the bounds there would have silently used the
!! MTC defaults and ignored the 0.20 m set for eucalyptus in orchidee_pft.def_51pft.4ac,
!! that is, produced wrong bounds for the very PFT that motivates the development.
!!
!! Only tree PFTs are touched: largest_tree_dia is meaningless for grasses and crops
!! (a zero would collapse every bound to zero and push any stand up a class at once).
!!
!! No-op when OK_AGE_CLASS_BOUND_PFT is false => bit-identical to the historical code.
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): age_class_bound
!!
!! REFERENCE(S) : design/MODULE_DESIGN_AGE_CLASS_BOUNDS_PFT.md
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE derive_age_class_bounds()

    IMPLICIT NONE

    !! 0.4 Local variables
    INTEGER(i_std) :: jv        !! PFT index
    INTEGER(i_std) :: ic        !! Age class index
    REAL(r_std)    :: cut_dia   !! Effective clearcut diameter of this PFT (m)

    !_ ================================================================================================================================

    IF (.NOT. ok_age_class_bound_pft) RETURN
    IF (.NOT. ALLOCATED(age_class_bound)) RETURN
    IF (nagec < 2) RETURN

    DO jv = 1, nvm

       IF (.NOT. is_tree(jv)) CYCLE
       IF (largest_tree_dia(jv) <= zero) CYCLE

       DO ic = 1, nagec-1
          age_class_bound(ic,jv) = age_class_frac(ic) * largest_tree_dia(jv)
       ENDDO
       ! Guillaume M. -- The upper bound of the last class stays open: it is a safeguard
       ! against aberrant diameters (ipslerr in age_class_distr), not a physical limit.

       ! Guillaume M. -- Bound versus clearcut threshold consistency, checked PFT by PFT
       ! and not only on the fractions: a LARGEST_TREE_DIA set in the namelist can break
       ! what the global guard of constantes.f90 validated on AGE_CLASS_FRAC alone.
       cut_dia = largest_tree_dia(jv)
       ! Guillaume M. -- The rotation tolerance is applied unconditionally, no flag gates it.
       cut_dia = cut_dia * (un - dia_rotation_tol)
       ! Guillaume M. -- Equality is admitted: the lower bound of the last class is derived
       ! from dia_rotation_tol (constantes.f90), so bound == clearcut threshold is the normal
       ! case. Only a threshold STRICTLY below the bound, i.e. outside the class, is rejected;
       ! per-PFT rounding of a namelist LARGEST_TREE_DIA can still produce it. The tolerance
       ! is relative to the cut diameter, not absolute: LTD ranges from 0.20 to 1 m.
       IF (age_class_bound(nagec-1,jv) > cut_dia * (un + 1.e-9)) THEN
          WRITE(numout,*) 'PFT, last bound, clearcut diameter: ', &
               jv, age_class_bound(nagec-1,jv), cut_dia
          CALL ipslerr_p(3,'derive_age_class_bounds',&
               'the clearcut threshold falls outside the last age class',&
               'lower AGE_CLASS_FRAC(nagec-1) or raise LARGEST_TREE_DIA for this PFT','')
       ENDIF

    ENDDO

    IF (printlev >= 1) THEN
       WRITE(numout,*) '[AGE_CLASS_BOUND_PFT] bounds derived from largest_tree_dia, ', &
            'AGE_CLASS_FRAC = ', age_class_frac(:)
       DO jv = 1, nvm
          IF (is_tree(jv) .AND. largest_tree_dia(jv) > zero) THEN
             WRITE(numout,*) '[AGE_CLASS_BOUND_PFT] PFT ', jv, ' LTD ', &
                  largest_tree_dia(jv), ' bounds ', age_class_bound(:,jv)
          ENDIF
       ENDDO
    ENDIF

  END SUBROUTINE derive_age_class_bounds


!! ================================================================================================================================
!! SUBROUTINE   : pft_parameters_clear
!!
!>\BRIEF         This subroutine deallocates memory at the end of the simulation. 
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE pft_parameters_clear

    l_first_pft_parameters = .TRUE.

    IF (ALLOCATED(pft_to_mtc)) DEALLOCATE(pft_to_mtc)
    IF (ALLOCATED(PFT_name)) DEALLOCATE(PFT_name)
    IF (ALLOCATED(llaimax)) DEALLOCATE(llaimax)
    IF (ALLOCATED(llaimin)) DEALLOCATE(llaimin)
    IF (ALLOCATED(height_presc)) DEALLOCATE(height_presc)   
    IF (ALLOCATED(z0_over_height)) DEALLOCATE(z0_over_height)   
    IF (ALLOCATED(ratio_z0m_z0h)) DEALLOCATE(ratio_z0m_z0h)   
    IF (ALLOCATED(type_of_lai)) DEALLOCATE(type_of_lai)
    IF (ALLOCATED(is_tree)) DEALLOCATE(is_tree)
    IF (ALLOCATED(natural)) DEALLOCATE(natural)
    IF (ALLOCATED(is_deciduous)) DEALLOCATE(is_deciduous)
    IF (ALLOCATED(is_evergreen)) DEALLOCATE(is_evergreen)
    IF (ALLOCATED(is_needleleaf)) DEALLOCATE(is_needleleaf)
    IF (ALLOCATED(is_tropical)) DEALLOCATE(is_tropical)
    IF (ALLOCATED(is_temperate)) DEALLOCATE(is_temperate)
    IF (ALLOCATED(is_boreal)) DEALLOCATE(is_boreal)
    IF (ALLOCATED(agec_group)) DEALLOCATE(agec_group)
    IF (ALLOCATED(start_index)) DEALLOCATE(start_index)
    IF (ALLOCATED(nagec_pft)) DEALLOCATE(nagec_pft)
    IF (ALLOCATED(humcste)) DEALLOCATE(humcste)
    IF (ALLOCATED(max_root_depth)) DEALLOCATE(max_root_depth)
    IF (ALLOCATED(pref_soil_veg)) DEALLOCATE(pref_soil_veg)
    IF (ALLOCATED(is_c4)) DEALLOCATE(is_c4)  
    IF (ALLOCATED(vcmax_fix)) DEALLOCATE(vcmax_fix)
    IF (ALLOCATED(downregulation_co2_coeff)) DEALLOCATE(downregulation_co2_coeff) 
    IF (ALLOCATED(E_KmC)) DEALLOCATE(E_KmC)
    IF (ALLOCATED(E_KmO)) DEALLOCATE(E_KmO)
    IF (ALLOCATED(E_Sco)) DEALLOCATE(E_Sco)
    IF (ALLOCATED(E_gamma_star)) DEALLOCATE(E_gamma_star)
    IF (ALLOCATED(E_Vcmax)) DEALLOCATE(E_Vcmax)
    IF (ALLOCATED(E_Jmax)) DEALLOCATE(E_Jmax)
    IF (ALLOCATED(aSV)) DEALLOCATE(aSV)
    IF (ALLOCATED(bSV)) DEALLOCATE(bSV)
    IF (ALLOCATED(tphoto_min)) DEALLOCATE(tphoto_min)
    IF (ALLOCATED(tphoto_max)) DEALLOCATE(tphoto_max)
    IF (ALLOCATED(aSJ)) DEALLOCATE(aSJ)
    IF (ALLOCATED(bSJ)) DEALLOCATE(bSJ)
    IF (ALLOCATED(D_Vcmax)) DEALLOCATE(D_Vcmax)
    IF (ALLOCATED(D_Jmax)) DEALLOCATE(D_Jmax)
    IF (ALLOCATED(E_gm)) DEALLOCATE(E_gm) 
    IF (ALLOCATED(S_gm)) DEALLOCATE(S_gm) 
    IF (ALLOCATED(D_gm)) DEALLOCATE(D_gm) 
    IF (ALLOCATED(E_Rd)) DEALLOCATE(E_Rd)
    IF (ALLOCATED(Vcmax25)) DEALLOCATE(Vcmax25)
    IF (ALLOCATED(arJV)) DEALLOCATE(arJV)
    IF (ALLOCATED(brJV)) DEALLOCATE(brJV)
    IF (ALLOCATED(KmC25)) DEALLOCATE(KmC25)
    IF (ALLOCATED(KmO25)) DEALLOCATE(KmO25)
    IF (ALLOCATED(Sco25)) DEALLOCATE(Sco25)
    IF (ALLOCATED(gm25)) DEALLOCATE(gm25) 
    IF (ALLOCATED(gamma_star25)) DEALLOCATE(gamma_star25)
    IF (ALLOCATED(a1)) DEALLOCATE(a1)
    IF (ALLOCATED(b1)) DEALLOCATE(b1)
    IF (ALLOCATED(a_fpsi)) DEALLOCATE(a_fpsi)
    IF (ALLOCATED(sf)) DEALLOCATE(sf)
    IF (ALLOCATED(wood_wat_store_fac)) DEALLOCATE(wood_wat_store_fac)
    IF (ALLOCATED(leaf_wat_store_fac)) DEALLOCATE(leaf_wat_store_fac)
    IF (ALLOCATED(Rroot_sup)) DEALLOCATE(Rroot_sup)
    IF (ALLOCATED(Rroot_inf)) DEALLOCATE(Rroot_inf)
    IF (ALLOCATED(Rleaf)) DEALLOCATE(Rleaf)
    IF (ALLOCATED(Rxylem)) DEALLOCATE(Rxylem)
    IF (ALLOCATED(Rroot_soil)) DEALLOCATE(Rroot_soil)
    IF (ALLOCATED(Rsto_leaf)) DEALLOCATE(Rsto_leaf)
    IF (ALLOCATED(Rsto_wood)) DEALLOCATE(Rsto_wood)
    IF (ALLOCATED(psi_ref_g0)) DEALLOCATE(psi_ref_g0)
    IF (ALLOCATED(psi_ref_sto_leaf)) DEALLOCATE(psi_ref_sto_leaf)
    IF (ALLOCATED(psi_ref_sto_wood)) DEALLOCATE(psi_ref_sto_wood)
    IF (ALLOCATED(psi_50_res_leaf)) DEALLOCATE(psi_50_res_leaf)
    IF (ALLOCATED(psi_50_res_wood)) DEALLOCATE(psi_50_res_wood)
    IF (ALLOCATED(psi_50_res_root)) DEALLOCATE(psi_50_res_root)
    IF (ALLOCATED(a_leaf)) DEALLOCATE(a_leaf)
    IF (ALLOCATED(a_wood)) DEALLOCATE(a_wood)
    IF (ALLOCATED(a_root)) DEALLOCATE(a_root)
    IF (ALLOCATED(k_leaf_max)) DEALLOCATE(k_leaf_max)
    IF (ALLOCATED(sapwood_density)) DEALLOCATE(sapwood_density)
    IF (ALLOCATED(ks_max)) DEALLOCATE(ks_max)
    IF (ALLOCATED(d_l)) DEALLOCATE(d_l)
    IF (ALLOCATED(k_wood_max)) DEALLOCATE(k_wood_max)
    IF (ALLOCATED(k_root_max)) DEALLOCATE(k_root_max)
    IF (ALLOCATED(lambda_leaf)) DEALLOCATE(lambda_leaf)
    IF (ALLOCATED(lambda_wood)) DEALLOCATE(lambda_wood)
    IF (ALLOCATED(g0)) DEALLOCATE(g0)
    IF (ALLOCATED(h_protons)) DEALLOCATE(h_protons)
    IF (ALLOCATED(fpsir)) DEALLOCATE(fpsir)
    IF (ALLOCATED(fQ)) DEALLOCATE(fQ)
    IF (ALLOCATED(fpseudo)) DEALLOCATE(fpseudo)
    IF (ALLOCATED(kp)) DEALLOCATE(kp)
    IF (ALLOCATED(alpha)) DEALLOCATE(alpha)
    IF (ALLOCATED(gbs)) DEALLOCATE(gbs)
    IF (ALLOCATED(theta)) DEALLOCATE(theta)
    IF (ALLOCATED(alpha_LL)) DEALLOCATE(alpha_LL)
    IF (ALLOCATED(stress_vcmax)) DEALLOCATE(stress_vcmax)
    IF (ALLOCATED(stress_gs)) DEALLOCATE(stress_gs)
    IF (ALLOCATED(stress_gm)) DEALLOCATE(stress_gm)
    IF (ALLOCATED(ext_coeff)) DEALLOCATE(ext_coeff)
    IF (ALLOCATED(ext_coeff_vegetfrac)) DEALLOCATE(ext_coeff_vegetfrac)
    IF (ALLOCATED(rveg_pft)) DEALLOCATE(rveg_pft)
    IF (ALLOCATED(rstruct_const)) DEALLOCATE(rstruct_const)
    IF (ALLOCATED(kzero)) DEALLOCATE(kzero)
    IF (ALLOCATED(alpha_watstress_pft)) DEALLOCATE(alpha_watstress_pft)
    IF (ALLOCATED(wmax_veg)) DEALLOCATE(wmax_veg)
    IF (ALLOCATED(throughfall_by_pft)) DEALLOCATE(throughfall_by_pft)
    IF (ALLOCATED(snowa_aged_vis)) DEALLOCATE(snowa_aged_vis)
    IF (ALLOCATED(snowa_aged_nir)) DEALLOCATE(snowa_aged_nir)
    IF (ALLOCATED(snowa_dec_vis)) DEALLOCATE(snowa_dec_vis)
    IF (ALLOCATED(snowa_dec_nir)) DEALLOCATE(snowa_dec_nir)
    IF (ALLOCATED(alb_leaf_vis)) DEALLOCATE(alb_leaf_vis)
    IF (ALLOCATED(alb_leaf_nir)) DEALLOCATE(alb_leaf_nir) 
    IF (ALLOCATED(leaf_ssa)) DEALLOCATE(leaf_ssa)   
    IF (ALLOCATED(leaf_psd)) DEALLOCATE(leaf_psd)   
    IF (ALLOCATED(bgd_reflectance)) DEALLOCATE(bgd_reflectance)
    IF (ALLOCATED(leaf_to_shoot_clumping)) DEALLOCATE(leaf_to_shoot_clumping) 
    IF (ALLOCATED(lai_correction_factor)) DEALLOCATE(lai_correction_factor) 
    IF (ALLOCATED(min_level_sep)) DEALLOCATE(min_level_sep)
    IF (ALLOCATED(lai_top)) DEALLOCATE(lai_top)   
    IF (ALLOCATED(em_factor_isoprene)) DEALLOCATE(em_factor_isoprene)
    IF (ALLOCATED(em_factor_monoterpene)) DEALLOCATE(em_factor_monoterpene)
    IF (ALLOCATED(em_factor_apinene)) DEALLOCATE(em_factor_apinene)
    IF (ALLOCATED(em_factor_bpinene)) DEALLOCATE(em_factor_bpinene)
    IF (ALLOCATED(em_factor_limonene)) DEALLOCATE(em_factor_limonene)
    IF (ALLOCATED(em_factor_myrcene)) DEALLOCATE(em_factor_myrcene)
    IF (ALLOCATED(em_factor_sabinene)) DEALLOCATE(em_factor_sabinene)
    IF (ALLOCATED(em_factor_camphene)) DEALLOCATE(em_factor_camphene)
    IF (ALLOCATED(em_factor_3carene)) DEALLOCATE(em_factor_3carene)
    IF (ALLOCATED(em_factor_tbocimene)) DEALLOCATE(em_factor_tbocimene)
    IF (ALLOCATED(em_factor_othermonot)) DEALLOCATE(em_factor_othermonot)
    IF (ALLOCATED(em_factor_sesquiterp)) DEALLOCATE(em_factor_sesquiterp)
    IF (ALLOCATED(em_factor_ORVOC)) DEALLOCATE(em_factor_ORVOC)
    IF (ALLOCATED(em_factor_OVOC)) DEALLOCATE(em_factor_OVOC)
    IF (ALLOCATED(em_factor_MBO)) DEALLOCATE(em_factor_MBO)
    IF (ALLOCATED(em_factor_methanol)) DEALLOCATE(em_factor_methanol)
    IF (ALLOCATED(em_factor_acetone)) DEALLOCATE(em_factor_acetone)
    IF (ALLOCATED(em_factor_acetal)) DEALLOCATE(em_factor_acetal)
    IF (ALLOCATED(em_factor_formal)) DEALLOCATE(em_factor_formal)
    IF (ALLOCATED(em_factor_acetic)) DEALLOCATE(em_factor_acetic)
    IF (ALLOCATED(em_factor_formic)) DEALLOCATE(em_factor_formic)
    IF (ALLOCATED(em_factor_DMS)) DEALLOCATE(em_factor_DMS)
    IF (ALLOCATED(em_factor_H2S)) DEALLOCATE(em_factor_H2S)
    IF (ALLOCATED(em_factor_no_wet)) DEALLOCATE(em_factor_no_wet)
    IF (ALLOCATED(em_factor_no_dry)) DEALLOCATE(em_factor_no_dry)
    IF (ALLOCATED(Larch)) DEALLOCATE(Larch)
    IF (ALLOCATED(leaf_tab)) DEALLOCATE(leaf_tab)
    IF (ALLOCATED(sla)) DEALLOCATE(sla)
    IF (ALLOCATED(slainit)) DEALLOCATE(slainit)
    IF (ALLOCATED(availability_fact)) DEALLOCATE(availability_fact)
    IF (ALLOCATED(nue_opt)) DEALLOCATE(nue_opt)
    IF (ALLOCATED(vmax_uptake)) DEALLOCATE(vmax_uptake)
    IF (ALLOCATED(ext_coeff_N)) DEALLOCATE(ext_coeff_N)
    IF (ALLOCATED(frac_growthresp)) DEALLOCATE(frac_growthresp)
    IF (ALLOCATED(coeff_maint_init)) DEALLOCATE(coeff_maint_init)
    IF (ALLOCATED(tref_maint_resp)) DEALLOCATE(tref_maint_resp)
    IF (ALLOCATED(tmin_maint_resp)) DEALLOCATE(tmin_maint_resp)
    IF (ALLOCATED(e0_maint_resp)) DEALLOCATE(e0_maint_resp)
    IF (ALLOCATED(tref_labile)) DEALLOCATE(tref_labile)
    IF (ALLOCATED(tmin_labile)) DEALLOCATE(tmin_labile)
    IF (ALLOCATED(e0_labile)) DEALLOCATE(e0_labile)
    IF (ALLOCATED(always_labile)) DEALLOCATE(always_labile)
    IF (ALLOCATED(coeff_lcchange_s)) DEALLOCATE(coeff_lcchange_s)
    IF (ALLOCATED(coeff_lcchange_m)) DEALLOCATE(coeff_lcchange_m)
    IF (ALLOCATED(coeff_lcchange_l)) DEALLOCATE(coeff_lcchange_l)
    IF (ALLOCATED(lai_max_to_happy)) DEALLOCATE(lai_max_to_happy)
    IF (ALLOCATED(lai_max)) DEALLOCATE(lai_max)
    IF (ALLOCATED(pheno_model)) DEALLOCATE(pheno_model)
    IF (ALLOCATED(gdd_threshold)) DEALLOCATE(gdd_threshold)
    IF (ALLOCATED(pheno_type)) DEALLOCATE(pheno_type)
    IF (ALLOCATED(force_pheno)) DEALLOCATE(force_pheno)
    IF (ALLOCATED(pheno_gdd_crit_c)) DEALLOCATE(pheno_gdd_crit_c)
    IF (ALLOCATED(pheno_gdd_crit_b)) DEALLOCATE(pheno_gdd_crit_b)
    IF (ALLOCATED(pheno_gdd_crit_a)) DEALLOCATE(pheno_gdd_crit_a)
    IF (ALLOCATED(pheno_gdd_crit)) DEALLOCATE(pheno_gdd_crit)
    IF (ALLOCATED(pheno_moigdd_t_crit)) DEALLOCATE(pheno_moigdd_t_crit)
    IF (ALLOCATED(ngd_crit)) DEALLOCATE(ngd_crit)
    IF (ALLOCATED(ncdgdd_temp)) DEALLOCATE(ncdgdd_temp)
    IF (ALLOCATED(hum_frac)) DEALLOCATE(hum_frac)
    IF (ALLOCATED(hum_min_time)) DEALLOCATE(hum_min_time)
    IF (ALLOCATED(relsoilmoist_always)) DEALLOCATE(relsoilmoist_always)
    IF (ALLOCATED(t_always)) DEALLOCATE(t_always)
    IF (ALLOCATED(harvest_time_a)) DEALLOCATE(harvest_time_a)
    IF (ALLOCATED(harvest_time_b)) DEALLOCATE(harvest_time_b)
    IF (ALLOCATED(longevity_sap)) DEALLOCATE(longevity_sap)
    IF (ALLOCATED(longevity_leaf)) DEALLOCATE(longevity_leaf)
    IF (ALLOCATED(leaf_age_crit_tref)) DEALLOCATE(leaf_age_crit_tref)
    IF (ALLOCATED(leaf_age_crit_coeff1)) DEALLOCATE(leaf_age_crit_coeff1)
    IF (ALLOCATED(leaf_age_crit_coeff2)) DEALLOCATE(leaf_age_crit_coeff2)
    IF (ALLOCATED(leaf_age_crit_coeff3)) DEALLOCATE(leaf_age_crit_coeff3)
    IF (ALLOCATED(longevity_fruit)) DEALLOCATE(longevity_fruit)
    IF (ALLOCATED(longevity_root)) DEALLOCATE(longevity_root)
    IF (ALLOCATED(ecureuil)) DEALLOCATE(ecureuil)
    IF (ALLOCATED(alloc_min)) DEALLOCATE(alloc_min)
    IF (ALLOCATED(alloc_max)) DEALLOCATE(alloc_max)
    IF (ALLOCATED(demi_alloc)) DEALLOCATE(demi_alloc)
    IF (ALLOCATED(leaffall)) DEALLOCATE(leaffall)
    IF (ALLOCATED(presenescence_ratio)) DEALLOCATE(presenescence_ratio)
    IF (ALLOCATED(senescence_ratio)) DEALLOCATE(senescence_ratio)
    IF (ALLOCATED(senescence_type)) DEALLOCATE(senescence_type)
    IF (ALLOCATED(senescence_hum)) DEALLOCATE(senescence_hum)
    IF (ALLOCATED(nosenescence_hum)) DEALLOCATE(nosenescence_hum)
    IF (ALLOCATED(max_turnover_time)) DEALLOCATE(max_turnover_time)
    IF (ALLOCATED(min_turnover_time)) DEALLOCATE(min_turnover_time)
    IF (ALLOCATED(recycle_leaf)) DEALLOCATE(recycle_leaf)
    IF (ALLOCATED(recycle_root)) DEALLOCATE(recycle_root)
    IF (ALLOCATED(min_leaf_age_for_senescence)) DEALLOCATE(min_leaf_age_for_senescence)
    IF (ALLOCATED(senescence_temp_c)) DEALLOCATE(senescence_temp_c)
    IF (ALLOCATED(senescence_temp_b)) DEALLOCATE(senescence_temp_b)
    IF (ALLOCATED(senescence_temp_a)) DEALLOCATE(senescence_temp_a)
    IF (ALLOCATED(senescence_temp)) DEALLOCATE(senescence_temp)
    IF (ALLOCATED(gdd_senescence)) DEALLOCATE(gdd_senescence)
    IF (ALLOCATED(pre_gdd)) DEALLOCATE(pre_gdd)
    IF (ALLOCATED(pre_harv)) DEALLOCATE(pre_harv)
    IF (ALLOCATED(gr_gt)) DEALLOCATE(gr_gt)
    IF (ALLOCATED(always_init)) DEALLOCATE(always_init)
    IF (ALLOCATED(cn_leaf_min)) DEALLOCATE(cn_leaf_min)
    IF (ALLOCATED(cn_leaf_max)) DEALLOCATE(cn_leaf_max)
    IF (ALLOCATED(max_soil_n_bnf)) DEALLOCATE(max_soil_n_bnf)
    IF (ALLOCATED(manure_pftweight)) DEALLOCATE(manure_pftweight)
    IF (ALLOCATED(residence_time)) DEALLOCATE(residence_time)
    IF (ALLOCATED(tmin_crit)) DEALLOCATE(tmin_crit)
    IF (ALLOCATED(tcm_crit)) DEALLOCATE(tcm_crit)
    IF (ALLOCATED(lai_initmin)) DEALLOCATE(lai_initmin)
    IF (ALLOCATED(bm_sapl)) DEALLOCATE(bm_sapl)
    IF (ALLOCATED(migrate)) DEALLOCATE(migrate)
    IF (ALLOCATED(cn_sapl)) DEALLOCATE(cn_sapl)
    IF (ALLOCATED(k_latosa_max)) DEALLOCATE(k_latosa_max)
    IF (ALLOCATED(k_latosa_min)) DEALLOCATE(k_latosa_min)
    IF (ALLOCATED(LC)) DEALLOCATE(LC)
    IF (ALLOCATED(LC_leaf)) DEALLOCATE(LC_leaf)
    IF (ALLOCATED(LC_sapabove)) DEALLOCATE(LC_sapabove)
    IF (ALLOCATED(LC_sapbelow)) DEALLOCATE(LC_sapbelow)
    IF (ALLOCATED(LC_heartabove)) DEALLOCATE(LC_heartabove)
    IF (ALLOCATED(LC_heartbelow)) DEALLOCATE(LC_heartbelow)
    IF (ALLOCATED(LC_fruit)) DEALLOCATE(LC_fruit)
    IF (ALLOCATED(LC_root)) DEALLOCATE(LC_root)
    IF (ALLOCATED(LC_carbres)) DEALLOCATE(LC_carbres)
    IF (ALLOCATED(LC_labile)) DEALLOCATE(LC_labile)
    IF (ALLOCATED(decomp_factor)) DEALLOCATE(decomp_factor)
    IF (ALLOCATED(crown_vertohor_dia)) DEALLOCATE(crown_vertohor_dia)
    IF (ALLOCATED(crown_to_height)) DEALLOCATE(crown_to_height)   
    IF (ALLOCATED(pipe_density)) DEALLOCATE(pipe_density)
    IF (ALLOCATED(tree_ff)) DEALLOCATE(tree_ff)
    IF (ALLOCATED(min_pipe_tune2)) DEALLOCATE(min_pipe_tune2)
    IF (ALLOCATED(pipe_tune3)) DEALLOCATE(pipe_tune3)
    IF (ALLOCATED(pipe_tune4)) DEALLOCATE(pipe_tune4)
    IF (ALLOCATED(pipe_k1)) DEALLOCATE(pipe_k1)
    IF (ALLOCATED(mass_ratio_heart_sap)) DEALLOCATE(mass_ratio_heart_sap)
    IF (ALLOCATED(canopy_cover)) DEALLOCATE(canopy_cover)
    IF (ALLOCATED(nmaxplants)) DEALLOCATE(nmaxplants)
    IF (ALLOCATED(height_init)) DEALLOCATE(height_init)
    IF (ALLOCATED(biomass_init)) DEALLOCATE(biomass_init)
    IF (ALLOCATED(qmd_init)) DEALLOCATE(qmd_init)
    IF (ALLOCATED(dia_thresh_inv)) DEALLOCATE(dia_thresh_inv)
    IF (ALLOCATED(beta_self_thinning)) DEALLOCATE(beta_self_thinning)
    IF (ALLOCATED(ref_alpha_self_thin)) DEALLOCATE(ref_alpha_self_thin)
    IF (ALLOCATED(init_pre_indust_ref_gpp)) DEALLOCATE(init_pre_indust_ref_gpp)
    IF (ALLOCATED(fuelwood_diameter)) DEALLOCATE(fuelwood_diameter)
    IF (ALLOCATED(coppice_kill_be_wood)) DEALLOCATE(coppice_kill_be_wood)
    IF (ALLOCATED(lai_to_height)) DEALLOCATE(lai_to_height)
    IF (ALLOCATED(deleuze_a)) DEALLOCATE(deleuze_a)
    IF (ALLOCATED(deleuze_b)) DEALLOCATE(deleuze_b)
    IF (ALLOCATED(deleuze_p_all)) DEALLOCATE(deleuze_p_all)
    IF (ALLOCATED(deleuze_power_a)) DEALLOCATE(deleuze_power_a)
    IF (ALLOCATED(m_dv)) DEALLOCATE(m_dv)
    IF (ALLOCATED(forest_managed_forced)) DEALLOCATE(forest_managed_forced)
    IF (ALLOCATED(dens_target)) DEALLOCATE(dens_target)
    IF (ALLOCATED(thinstrat)) DEALLOCATE(thinstrat)
    IF (ALLOCATED(taumin)) DEALLOCATE(taumin)
    IF (ALLOCATED(taumax)) DEALLOCATE(taumax)
    IF (ALLOCATED(poly_rdi)) DEALLOCATE(poly_rdi)
    IF (ALLOCATED(rdi_dia_ref)) DEALLOCATE(rdi_dia_ref)
    IF (ALLOCATED(rdi_start)) DEALLOCATE(rdi_start)
    IF (ALLOCATED(rdi_max)) DEALLOCATE(rdi_max)
    IF (ALLOCATED(final_dia)) DEALLOCATE(final_dia)
    IF (ALLOCATED(rdi_ramp_out)) DEALLOCATE(rdi_ramp_out)
    IF (ALLOCATED(forced_thin_intensity)) DEALLOCATE(forced_thin_intensity)
    IF (ALLOCATED(rdi_band_low)) DEALLOCATE(rdi_band_low)
    IF (ALLOCATED(rdi_band_high)) DEALLOCATE(rdi_band_high)
    IF (ALLOCATED(thin_intensity)) DEALLOCATE(thin_intensity)
    IF (ALLOCATED(dia_boost_end)) DEALLOCATE(dia_boost_end)
    IF (ALLOCATED(circ_dist_shape_thin)) DEALLOCATE(circ_dist_shape_thin) 
    IF (ALLOCATED(circ_dist_shape_none)) DEALLOCATE(circ_dist_shape_none) 
    IF (ALLOCATED(circ_dist_shape_cop)) DEALLOCATE(circ_dist_shape_cop) 
    IF (ALLOCATED(circ_dist_shape_src)) DEALLOCATE(circ_dist_shape_src)
    IF (ALLOCATED(circ_dist_shape_uneven)) DEALLOCATE(circ_dist_shape_uneven)      
    IF (ALLOCATED(largest_tree_dia)) DEALLOCATE(largest_tree_dia)
    IF (ALLOCATED(coppice_diameter)) DEALLOCATE(coppice_diameter)
    IF (ALLOCATED(shoots_per_stool)) DEALLOCATE(shoots_per_stool)
    IF (ALLOCATED(src_rot_length)) DEALLOCATE(src_rot_length)
    IF (ALLOCATED(src_nrots)) DEALLOCATE(src_nrots)
    IF (ALLOCATED(fruit_alloc)) DEALLOCATE(fruit_alloc)
    IF (ALLOCATED(labile_reserve)) DEALLOCATE(labile_reserve)
    IF (ALLOCATED(evergreen_reserve)) DEALLOCATE(evergreen_reserve)
    IF (ALLOCATED(deciduous_reserve)) DEALLOCATE(deciduous_reserve)
    IF (ALLOCATED(senescense_reserve)) DEALLOCATE(senescense_reserve)
    IF (ALLOCATED(root_reserve)) DEALLOCATE(root_reserve)
    IF (ALLOCATED(fcn_wood)) DEALLOCATE(fcn_wood)
    IF (ALLOCATED(fcn_root)) DEALLOCATE(fcn_root)
    IF (ALLOCATED(branch_ratio)) DEALLOCATE(branch_ratio)
    IF (ALLOCATED(cn_leaf_init)) DEALLOCATE(cn_leaf_init)
    IF (ALLOCATED(moss_frac)) DEALLOCATE(moss_frac)
    IF (ALLOCATED(k_root)) DEALLOCATE(k_root)
    IF (ALLOCATED(k_belowground)) DEALLOCATE(k_belowground)
    IF (ALLOCATED(k_sap)) DEALLOCATE(k_sap)
    IF (ALLOCATED(psi_leaf_min)) DEALLOCATE(psi_leaf_min)
    IF (ALLOCATED(srl)) DEALLOCATE(srl)
    IF (ALLOCATED(r_froot)) DEALLOCATE(r_froot)
    IF (ALLOCATED(psi_root_min)) DEALLOCATE(psi_root_min)
    IF (ALLOCATED(som_turn_ipassive)) DEALLOCATE(som_turn_ipassive)
    IF (ALLOCATED(recruitment_pft)) DEALLOCATE(recruitment_pft) 
    IF (ALLOCATED(beetle_pft)) DEALLOCATE(beetle_pft) 
    IF (ALLOCATED(recruitment_height)) DEALLOCATE(recruitment_height)  
    IF (ALLOCATED(recruitment_alpha)) DEALLOCATE(recruitment_alpha)  
    IF (ALLOCATED(recruitment_beta)) DEALLOCATE(recruitment_beta) 
    IF (ALLOCATED(harvest_ratio)) DEALLOCATE(harvest_ratio)
    IF (ALLOCATED(death_distribution_factor)) DEALLOCATE(death_distribution_factor)
    IF (ALLOCATED(npp_reset_value)) DEALLOCATE(npp_reset_value)
    IF (ALLOCATED(ndying_year)) DEALLOCATE(ndying_year)
    IF (ALLOCATED(S_activity)) DEALLOCATE(S_activity)
    IF (ALLOCATED(act_limit)) DEALLOCATE(act_limit)
    IF (ALLOCATED(S_competition)) DEALLOCATE(S_competition)
    IF (ALLOCATED(RDi_limit)) DEALLOCATE(RDi_limit)
    IF (ALLOCATED(S_susceptibility)) DEALLOCATE(S_susceptibility)
    IF (ALLOCATED(i_rd_susceptibility)) DEALLOCATE(i_rd_susceptibility)
    IF (ALLOCATED(S_share)) DEALLOCATE(S_share)
    IF (ALLOCATED(SH_limit)) DEALLOCATE(SH_limit)
    IF (ALLOCATED(S_defence)) DEALLOCATE(S_defence)
    IF (ALLOCATED(PWS_limit)) DEALLOCATE(PWS_limit)
    IF (ALLOCATED(max_Nwood)) DEALLOCATE(max_Nwood)
    IF (ALLOCATED(beetle_generation_a))DEALLOCATE(beetle_generation_a)
    IF (ALLOCATED(beetle_generation_b))DEALLOCATE(beetle_generation_b)
    IF (ALLOCATED(beetle_generation_c))DEALLOCATE(beetle_generation_c)
    IF (ALLOCATED(min_temp_beetle))DEALLOCATE(min_temp_beetle)
    IF (ALLOCATED(max_temp_beetle))DEALLOCATE(max_temp_beetle)
    IF (ALLOCATED(opt_temp_beetle))DEALLOCATE(opt_temp_beetle)
    IF (ALLOCATED(eff_temp_beetle_a))DEALLOCATE(eff_temp_beetle_a)
    IF (ALLOCATED(eff_temp_beetle_b))DEALLOCATE(eff_temp_beetle_b)
    IF (ALLOCATED(eff_temp_beetle_c))DEALLOCATE(eff_temp_beetle_c)
    IF (ALLOCATED(eff_temp_beetle_d))DEALLOCATE(eff_temp_beetle_d)
    IF (ALLOCATED(diapause_thres_daylength))DEALLOCATE(diapause_thres_daylength)
    IF (ALLOCATED(S_massattack)) DEALLOCATE(S_massattack)
    IF (ALLOCATED(BP_limit)) DEALLOCATE(BP_limit)
    IF (ALLOCATED(streamlining_c_leaf)) DEALLOCATE(streamlining_c_leaf)
    IF (ALLOCATED(streamlining_c_leafless)) DEALLOCATE(streamlining_c_leafless)
    IF (ALLOCATED(streamlining_n_leaf)) DEALLOCATE(streamlining_n_leaf)
    IF (ALLOCATED(streamlining_n_leafless)) DEALLOCATE(streamlining_n_leafless)
    IF (ALLOCATED(modulus_rupture)) DEALLOCATE(modulus_rupture)
    IF (ALLOCATED(f_knot)) DEALLOCATE(f_knot)
    IF (ALLOCATED(overturning_free_draining_shallow)) DEALLOCATE(overturning_free_draining_shallow)
    IF (ALLOCATED(overturning_free_draining_shallow_leafless)) DEALLOCATE(overturning_free_draining_shallow_leafless)
    IF (ALLOCATED(overturning_free_draining_deep)) DEALLOCATE(overturning_free_draining_deep)
    IF (ALLOCATED(overturning_free_draining_deep_leafless)) DEALLOCATE(overturning_free_draining_deep_leafless)
    IF (ALLOCATED(overturning_free_draining_average)) DEALLOCATE(overturning_free_draining_average)
    IF (ALLOCATED(overturning_free_draining_average_leafless)) DEALLOCATE(overturning_free_draining_average_leafless)
    IF (ALLOCATED(overturning_gleyed_shallow)) DEALLOCATE(overturning_gleyed_shallow)
    IF (ALLOCATED(overturning_gleyed_shallow_leafless)) DEALLOCATE(overturning_gleyed_shallow_leafless)
    IF (ALLOCATED(overturning_gleyed_deep)) DEALLOCATE(overturning_gleyed_deep)
    IF (ALLOCATED(overturning_gleyed_deep_leafless)) DEALLOCATE(overturning_gleyed_deep_leafless)
    IF (ALLOCATED(overturning_gleyed_average)) DEALLOCATE(overturning_gleyed_average)
    IF (ALLOCATED(overturning_gleyed_average_leafless)) DEALLOCATE(overturning_gleyed_average_leafless)
    IF (ALLOCATED(overturning_peaty_shallow)) DEALLOCATE(overturning_peaty_shallow)
    IF (ALLOCATED(overturning_peaty_shallow_leafless)) DEALLOCATE(overturning_peaty_shallow_leafless)
    IF (ALLOCATED(overturning_peaty_deep)) DEALLOCATE(overturning_peaty_deep)
    IF (ALLOCATED(overturning_peaty_deep_leafless)) DEALLOCATE(overturning_peaty_deep_leafless)
    IF (ALLOCATED(overturning_peaty_average)) DEALLOCATE(overturning_peaty_average)
    IF (ALLOCATED(overturning_peaty_average_leafless)) DEALLOCATE(overturning_peaty_average_leafless)
    IF (ALLOCATED(overturning_peat_shallow)) DEALLOCATE(overturning_peat_shallow)
    IF (ALLOCATED(overturning_peat_shallow_leafless)) DEALLOCATE(overturning_peat_shallow_leafless)
    IF (ALLOCATED(overturning_peat_deep)) DEALLOCATE(overturning_peat_deep)
    IF (ALLOCATED(overturning_peat_deep_leafless)) DEALLOCATE(overturning_peat_deep_leafless)
    IF (ALLOCATED(overturning_peat_average)) DEALLOCATE(overturning_peat_average)
    IF (ALLOCATED(overturning_peat_average_leafless)) DEALLOCATE(overturning_peat_average_leafless)
    IF (ALLOCATED(max_damage_further)) DEALLOCATE(max_damage_further)
    IF (ALLOCATED(max_damage_closer)) DEALLOCATE(max_damage_closer)
    IF (ALLOCATED(sfactor_further)) DEALLOCATE(sfactor_further)
    IF (ALLOCATED(sfactor_closer)) DEALLOCATE(sfactor_closer)
    IF (ALLOCATED(green_density)) DEALLOCATE(green_density)
    IF (ALLOCATED(maint_resp_slope)) DEALLOCATE(maint_resp_slope) 
    IF (ALLOCATED(maint_resp_slope_c)) DEALLOCATE(maint_resp_slope_c) 
    IF (ALLOCATED(maint_resp_slope_b)) DEALLOCATE(maint_resp_slope_b) 
    IF (ALLOCATED(maint_resp_slope_a)) DEALLOCATE(maint_resp_slope_a)
    IF (ALLOCATED(is_peat)) DEALLOCATE(is_peat)
    IF (ALLOCATED(is_croppeat)) DEALLOCATE(is_croppeat)

  END SUBROUTINE pft_parameters_clear

END MODULE pft_parameters
